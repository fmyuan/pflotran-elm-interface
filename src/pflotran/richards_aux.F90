module Richards_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

#ifdef BUFFER_MATRIX
  use Matrix_Buffer_module
#endif

  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none
  
  private 

  PetscReal, public :: richards_itol_scaled_res = 1.d-5
  PetscReal, public :: richards_itol_rel_update = UNINITIALIZED_DOUBLE
  PetscInt, public :: richards_ni_count
  PetscInt, public :: richards_ts_cut_count
  PetscInt, public :: richards_ts_count

  type, public :: richards_auxvar_type
  
    PetscReal :: pc
    PetscReal :: kvr
    PetscReal :: kr
    PetscReal :: dkvr_dp
    PetscReal :: dsat_dp
    PetscReal :: dden_dp
#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
    PetscReal :: bc_alpha  ! Brooks Corey parameterization: alpha
    PetscReal :: bc_lambda ! Brooks Corey parameterization: lambda    
#endif

    PetscReal :: d2sat_dp2
    PetscReal :: d2den_dp2
    PetscReal :: mass
    PetscReal :: dpres_dtime
    PetscReal :: dmass_dtime

  end type richards_auxvar_type
  
  type, public :: richards_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: auxvars_cell_pressures_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(richards_auxvar_type), pointer :: auxvars(:)
    type(richards_auxvar_type), pointer :: auxvars_bc(:)
    type(richards_auxvar_type), pointer :: auxvars_ss(:)
#ifdef BUFFER_MATRIX
    type(matrix_buffer_type), pointer :: matrix_buffer
#endif
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type richards_type

  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public :: RichardsAuxCreate, RichardsAuxDestroy, &
            RichardsAuxVarCompute, RichardsAuxVarInit, &
            RichardsAuxVarCopy, &
            RichardsAuxVarCompute2ndOrderDeriv

contains

! ************************************************************************** !

function RichardsAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(richards_type), pointer :: RichardsAuxCreate
  
  type(richards_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%auxvars_cell_pressures_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
#ifdef BUFFER_MATRIX
  nullify(aux%matrix_buffer)
#endif
  nullify(aux%matrix_zeroing)

  RichardsAuxCreate => aux
  
end function RichardsAuxCreate

! ************************************************************************** !

subroutine RichardsAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(richards_auxvar_type) :: auxvar
  type(option_type) :: option
  
  auxvar%pc = 0.d0

  auxvar%kvr = 0.d0
  auxvar%kr = 0.d0
  auxvar%dkvr_dp = 0.d0

  auxvar%dsat_dp = 0.d0
  auxvar%dden_dp = 0.d0

  auxvar%d2sat_dp2 = 0.d0
  auxvar%d2den_dp2 = 0.d0
  auxvar%mass = 0.d0
  auxvar%dpres_dtime = 0.d0
  auxvar%dmass_dtime = 0.0d0

#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
  auxvar%bc_alpha  = 0.0d0
  auxvar%bc_lambda  = 0.0d0
#endif 
  
end subroutine RichardsAuxVarInit

! ************************************************************************** !

subroutine RichardsAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Option_module

  implicit none
  
  type(richards_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pc = auxvar%pc

  auxvar2%kvr = auxvar%kvr
  auxvar2%kr = auxvar%kr
  auxvar2%dkvr_dp = auxvar%dkvr_dp

  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dden_dp = auxvar%dden_dp
 
  auxvar2%d2sat_dp2 = auxvar%d2sat_dp2
  auxvar2%d2den_dp2 = auxvar%d2den_dp2
  auxvar2%mass = auxvar%mass
  auxvar2%dpres_dtime = auxvar%dpres_dtime
  auxvar2%dmass_dtime = auxvar%dmass_dtime

#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
  auxvar2%bc_alpha  = auxvar%bc_alpha
  auxvar2%bc_lambda = auxvar%bc_lambda
#endif

end subroutine RichardsAuxVarCopy

! ************************************************************************** !

subroutine RichardsAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, &
                                 characteristic_curves,natural_id, &
                                 update_porosity,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(richards_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  PetscBool :: update_porosity
  
  PetscInt :: i
  PetscBool :: saturated
  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: pert, pw_pert, dw_kg_pert
  PetscReal :: fs, ani_A, ani_B, ani_C, ani_n, ani_coef
  PetscReal :: dkr_sat
  PetscReal :: aux(1)
  PetscReal, parameter :: tol = 1.d-3
  PetscReal :: compressed_porosity, dcompressed_porosity_dp
  
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0
  
  auxvar%kvr = 0.d0
  auxvar%kr = 0.d0

  kr = 0.d0
 
  global_auxvar%pres = x(1)
  global_auxvar%temp = option%flow%reference_temperature

  if (update_porosity) then
    call MaterialAuxVarCompute(material_auxvar,global_auxvar%pres(1))
  endif
 
  ! For a very large negative liquid pressure (e.g. -1.d18), the capillary 
  ! pressure can go near infinite, resulting in ds_dp being < 1.d-40 below 
  ! and flipping the cell to saturated, when it is really far from saturated.
  ! The large negative liquid pressure is then passed to the EOS causing it 
  ! to blow up.  Therefore, we truncate to the max capillary pressure here.
  auxvar%pc = min(option%flow%reference_pressure - global_auxvar%pres(1), &
                  characteristic_curves%saturation_function%pcmax)
  
!***************  Liquid phase properties **************************
  pw = option%flow%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0

  if (auxvar%pc > 0.d0) then
#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
    if (auxvar%bc_alpha > 0.d0) then
      select type(sf => characteristic_curves%saturation_function)
        class is(sat_func_VG_type)
          sf%m     = auxvar%bc_lambda
          sf%alpha = auxvar%bc_alpha
        class is(sat_func_BC_type)
            sf%lambda = auxvar%bc_lambda
            sf%alpha  = auxvar%bc_alpha
        class default
          option%io_buffer = 'CLM-PFLOTRAN only supports ' // &
            'sat_func_VG_type and sat_func_BC_type'
          call PrintErrMsg(option)
      end select

      select type(rpf => characteristic_curves%liq_rel_perm_function)
        class is(rpf_Mualem_VG_liq_type)
          rpf%m = auxvar%bc_lambda
        class is(rpf_Burdine_BC_liq_type)
          rpf%lambda = auxvar%bc_lambda
        class is(rpf_Mualem_BC_liq_type)
          rpf%lambda = auxvar%bc_lambda
        class is(rpf_Burdine_VG_liq_type)
          rpf%m = auxvar%bc_lambda
        class default
          option%io_buffer = 'Unsupported LIQUID-REL-PERM-FUNCTION'
          call PrintErrMsg(option)
      end select
    endif
#endif
    saturated = PETSC_FALSE
    call characteristic_curves%saturation_function% &
                               Saturation(auxvar%pc,global_auxvar%sat(1), &
                                          ds_dp,option)
    ! if ds_dp is 0, we consider the cell saturated.
    if (ds_dp < 1.d-40) then
      saturated = PETSC_TRUE
    else
      call characteristic_curves%liq_rel_perm_function% &
                       RelativePermeability(global_auxvar%sat(1), &
                                            kr,dkr_sat,option)
      dkr_dp = ds_dp * dkr_sat
    endif
  else
    saturated = PETSC_TRUE
  endif  
  
  ! the purpose for splitting this condition from the 'else' statement
  ! above is due to SaturationFunctionCompute switching a cell to
  ! saturated to prevent unstable (potentially infinite) derivatives when 
  ! capillary pressure is very small
  if (saturated) then
    auxvar%pc = 0.d0
    global_auxvar%sat = 1.d0  
    kr = 1.d0    
    pw = global_auxvar%pres(1)
  endif

  if (.not.option%flow%density_depends_on_salinity) then
    call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol, &
                         dw_dp,dw_dt,ierr)
    if (ierr /= 0) then
      call PrintMsgByCell(option,natural_id, &
                          'Error in RichardsAuxVarCompute->EOSWaterDensity')
    endif
    ! may need to compute dpsat_dt to pass to VISW
    call EOSWaterSaturationPressure(global_auxvar%temp,sat_pressure,ierr)
    !geh: 0.d0 passed in for derivative of pressure w/respect to temp
    call EOSWaterViscosity(global_auxvar%temp,pw,sat_pressure,0.d0, &
                           visl,dvis_dt,dvis_dp,ierr) 
  else
    aux(1) = global_auxvar%m_nacl(1)
    call EOSWaterDensityExt(global_auxvar%temp,pw,aux, &
                            dw_kg,dw_mol,dw_dp,dw_dt,ierr)
    if (ierr /= 0) then
      call PrintMsgByCell(option,natural_id, &
                          'Error in RichardsAuxVarCompute->EOSWaterDensityExt')
    endif
    call EOSWaterViscosityExt(global_auxvar%temp,pw,sat_pressure,0.d0,aux, &
                              visl,dvis_dt,dvis_dp,ierr) 
  endif
  if (.not.saturated) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif
 
  global_auxvar%den = dw_mol
  global_auxvar%den_kg = dw_kg
  auxvar%dsat_dp = ds_dp
  auxvar%dden_dp = dw_dp
  auxvar%kr = kr  ! stored solely for output purposes
  auxvar%kvr = kr/visl
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp

  if (size(global_auxvar%sat) > 1) then
    global_auxvar%sat(2) = 1.d0 - global_auxvar%sat(1)
  endif
  
end subroutine RichardsAuxVarCompute

! ************************************************************************** !
subroutine RichardsAuxVarCompute2ndOrderDeriv(rich_auxvar,global_auxvar, &
                                              material_auxvar, &
                                              characteristic_curves, &
                                              option)

  ! Computes 2nd order derivatives auxiliary variables for each grid cell
  ! 
  ! Author: Gautam Bisht, Satish Karra
  ! Date: 07/02/18, 06/26/19
  ! 

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  type(richards_auxvar_type) :: rich_auxvar, rich_auxvar_pert
  type(global_auxvar_type) :: global_auxvar, global_auxvar_pert
  class(material_auxvar_type) :: material_auxvar
  type(material_auxvar_type) :: material_auxvar_pert
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), pert
  PetscInt :: ideriv
  PetscErrorCode :: ierr

  rich_auxvar%d2sat_dp2 = 0.d0
  rich_auxvar%d2den_dp2 = 0.d0

  call GlobalAuxVarInit(global_auxvar_pert,option)  
  call MaterialAuxVarInit(material_auxvar_pert,option)  
  call RichardsAuxVarCopy(rich_auxvar,rich_auxvar_pert,option)
  call GlobalAuxVarCopy(global_auxvar,global_auxvar_pert,option)
  call MaterialAuxVarCopy(material_auxvar,material_auxvar_pert,option)
  x(1) = global_auxvar%pres(1)
  
  ideriv = 1
  pert = max(dabs(x(ideriv)*perturbation_tolerance),0.1d0)
  x_pert = x
  if (x_pert(ideriv) < option%flow%reference_pressure) pert = -1.d0*pert
  x_pert(ideriv) = x_pert(ideriv) + pert

  call RichardsAuxVarCompute(x_pert(1),rich_auxvar_pert,global_auxvar_pert, &
                       material_auxvar_pert, &
                       characteristic_curves, &
                       -999, PETSC_TRUE, &
                       option)   

  rich_auxvar%d2den_dp2 = (rich_auxvar_pert%dden_dp - rich_auxvar%dden_dp)/pert
  if (rich_auxvar%pc > 0.d0) then
    call characteristic_curves%saturation_function% &
                               D2SatDP2(rich_auxvar%pc, &
                                          rich_auxvar%d2sat_dp2,option)
  endif

    
  call GlobalAuxVarStrip(global_auxvar_pert)  
  call MaterialAuxVarStrip(material_auxvar_pert)  

end subroutine RichardsAuxVarCompute2ndOrderDeriv

! ************************************************************************** !

subroutine AuxVarDestroy(auxvar)
  ! 
  ! Deallocates a richards auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  implicit none

  type(richards_auxvar_type) :: auxvar
  
end subroutine AuxVarDestroy

! ************************************************************************** !

subroutine RichardsAuxDestroy(aux)
  ! 
  ! Deallocates a richards auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none

  type(richards_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call AuxVarDestroy(aux%auxvars(iaux))
    enddo  
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call AuxVarDestroy(aux%auxvars_bc(iaux))
    enddo  
    deallocate(aux%auxvars_bc)
  endif
  nullify(aux%auxvars_bc)
  if (associated(aux%auxvars_ss)) then
    do iaux = 1, aux%num_aux_ss
      call AuxVarDestroy(aux%auxvars_ss(iaux))
    enddo  
    deallocate(aux%auxvars_ss)
  endif
  nullify(aux%auxvars_ss)
  
  call MatrixZeroingDestroy(aux%matrix_zeroing)

#ifdef BUFFER_MATRIX
  if (associated(aux%matrix_buffer)) then
    call MatrixBufferDestroy(aux%matrix_buffer)
  endif
  nullify(aux%matrix_buffer)
#endif
  
  deallocate(aux)
  nullify(aux)
    
end subroutine RichardsAuxDestroy

end module Richards_Aux_module

