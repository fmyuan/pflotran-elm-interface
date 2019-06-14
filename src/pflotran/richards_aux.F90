module Richards_Aux_module

#ifdef BUFFER_MATRIX
  use Matrix_Buffer_module
#endif

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

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

    ! OLD-VAR-NAMES            = NEW-VAR
    ! ------------------------------------------------
    ! P_min                    = vars_for_sflow(1)
    ! P_max                    = vars_for_sflow(2)
    ! coeff_for_cubic_approx   = vars_for_sflow(3:6)
    ! range_for_linear_approx  = vars_for_sflow(7:10)
    ! bcflux_default_scheme    = vars_for_sflow(11)
    PetscReal, pointer :: vars_for_sflow(:)

  end type richards_auxvar_type
  
  type, public :: richards_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

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
  end type richards_type

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
#include "petsc/finclude/petscsys.h"
  use petscsys

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
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)
#ifdef BUFFER_MATRIX
  nullify(aux%matrix_buffer)
#endif

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

  if (option%surf_flow_on) then
    allocate(auxvar%vars_for_sflow(11))
    auxvar%vars_for_sflow(:) = 0.d0
  else
    nullify(auxvar%vars_for_sflow)
  endif

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

  if (option%surf_flow_on) &
    auxvar2%vars_for_sflow(:) = auxvar%vars_for_sflow(:)

#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
  auxvar2%bc_alpha  = auxvar%bc_alpha
  auxvar2%bc_lambda = auxvar%bc_lambda
#endif

end subroutine RichardsAuxVarCopy

! ************************************************************************** !

subroutine RichardsAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, &
                                 characteristic_curves,natural_id,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
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
  
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0
  
  auxvar%kvr = 0.d0
  auxvar%kr = 0.d0

  kr = 0.d0
 
  global_auxvar%pres = x(1)
  global_auxvar%temp = option%reference_temperature
 
  ! For a very large negative liquid pressure (e.g. -1.d18), the capillary 
  ! pressure can go near infinite, resulting in ds_dp being < 1.d-40 below 
  ! and flipping the cell to saturated, when it is really far from saturated.
  ! The large negative liquid pressure is then passed to the EOS causing it 
  ! to blow up.  Therefore, we truncate to the max capillary pressure here.
  auxvar%pc = min(option%reference_pressure - global_auxvar%pres(1), &
                  characteristic_curves%saturation_function%pcmax)
  
!***************  Liquid phase properties **************************
  pw = option%reference_pressure
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
          call printErrMsg(option)
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
          call printErrMsg(option)
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
      call printMsgByCell(option,natural_id, &
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
      call printMsgByCell(option,natural_id, &
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
  
end subroutine RichardsAuxVarCompute

! ************************************************************************** !
subroutine RichardsAuxVarCompute2ndOrderDeriv(auxvar,characteristic_curves,option)

  ! Computes 2nd order derivatives auxiliary variables for each grid cell
  ! 
  ! Author: Gautam Bisht
  ! Date: 07/02/18
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  type(richards_auxvar_type) :: auxvar

  PetscReal, parameter :: dp = 1.d-4
  PetscReal :: d2s_dp2
  PetscReal :: d2den_dp2
  PetscReal :: pw1,pw2
  PetscReal :: dw_dp_1,dw_dp_2
  PetscReal :: dw_kg,dw_mol,dw_dt
  PetscBool :: saturated
  PetscErrorCode :: ierr

  auxvar%d2sat_dp2 = 0.d0
  d2s_dp2 = 0.d0
  d2den_dp2 = 0.d0

  pw1 = option%reference_pressure
  saturated = PETSC_FALSE

  if (auxvar%pc > 0.d0) then
    call characteristic_curves%saturation_function% &
                               D2SatDP2(auxvar%pc, &
                                          d2s_dp2,option)
  endif

  pw2 = pw1 + dp
  call EOSWaterDensity(option%reference_temperature,pw1,dw_kg,dw_mol, &
                       dw_dp_1,dw_dt,ierr)
  call EOSWaterDensity(option%reference_temperature,pw2,dw_kg,dw_mol, &
                       dw_dp_2,dw_dt,ierr)
  d2den_dp2 = (dw_dp_2 - dw_dp_2)/dp

  auxvar%d2sat_dp2 = d2s_dp2
  auxvar%d2den_dp2 = d2den_dp2

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
  
  call DeallocateArray(aux%zero_rows_local)
  call DeallocateArray(aux%zero_rows_local_ghosted)

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
