module ZFlow_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none

  private

  PetscReal, parameter, public :: zflow_density_kg = 998.32d0
  PetscReal, parameter, public :: zflow_density_kmol = zflow_density_kg / FMWH2O
  PetscReal, parameter, public :: zflow_viscosity = 8.9d-4

  PetscReal, public :: zflow_pres_rel_pert = 1.d-8
  PetscReal, public :: zflow_pres_min_pert = 1.d-2

  PetscBool, public :: zflow_calc_accum = PETSC_TRUE
  PetscBool, public :: zflow_calc_flux = PETSC_TRUE
  PetscBool, public :: zflow_calc_bcflux = PETSC_TRUE
  PetscBool, public :: zflow_calc_chem = PETSC_TRUE

  PetscBool, public :: zflow_numerical_derivatives = PETSC_FALSE
  PetscBool, public :: zflow_simult_function_evals = PETSC_TRUE

  ! debugging
  PetscInt, public :: zflow_ni_count
  PetscInt, public :: zflow_ts_cut_count
  PetscInt, public :: zflow_ts_count

  PetscInt, parameter, public :: ZFLOW_LIQUID_PRESSURE_DOF = 1

  PetscInt, parameter, public :: ZFLOW_LIQUID_EQUATION_INDEX = 1

  PetscInt, parameter, public :: ZFLOW_LIQUID_PRESSURE_INDEX = 1
  PetscInt, parameter, public :: ZFLOW_LIQUID_FLUX_INDEX = 1
  PetscInt, parameter, public :: ZFLOW_LIQUID_CONDUCTANCE_INDEX = 2
  PetscInt, parameter, public :: ZFLOW_MAX_INDEX = 2

  PetscInt, parameter, public :: ZFLOW_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: ZFLOW_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: ZFLOW_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: ZFLOW_UPDATE_FOR_BOUNDARY = 2

  PetscInt, parameter, public :: ZFLOW_ADJOINT_PERMEABILITY = 1
  PetscInt, parameter, public :: ZFLOW_ADJOINT_POROSITY = 2
  PetscBool, public :: zflow_calc_adjoint = PETSC_FALSE
  PetscInt, public :: zflow_adjoint_parameter = ZFLOW_ADJOINT_PERMEABILITY

  type, public :: zflow_auxvar_type
    PetscReal :: pres ! liquid pressure
    PetscReal :: sat  ! liquid saturation
    PetscReal :: pc   ! capillary pressure
    PetscReal :: kr   ! relative permeability
    PetscReal :: effective_porosity
    PetscReal :: dpor_dp
    PetscReal :: dsat_dp ! derivative of saturation wrt pressure
    PetscReal :: dkr_dp  ! derivative of rel. perm. wrt pressure
    PetscReal :: pert
  end type zflow_auxvar_type

  type, public :: zflow_parameter_type
    PetscBool :: check_post_converged
  end type zflow_parameter_type

  type, public :: zflow_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(zflow_parameter_type), pointer :: zflow_parameter
    type(zflow_auxvar_type), pointer :: auxvars(:,:)
    type(zflow_auxvar_type), pointer :: auxvars_bc(:)
    type(zflow_auxvar_type), pointer :: auxvars_ss(:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type zflow_type

  interface ZFlowAuxVarDestroy
    module procedure ZFlowAuxVarSingleDestroy
    module procedure ZFlowAuxVarArray1Destroy
    module procedure ZFlowAuxVarArray2Destroy
  end interface ZFlowAuxVarDestroy

  interface ZFlowOutputAuxVars
    module procedure ZFlowOutputAuxVars1
  end interface ZFlowOutputAuxVars

  public :: ZFlowAuxCreate, &
            ZFlowAuxDestroy, &
            ZFlowAuxVarCompute, &
            ZFlowAuxVarInit, &
            ZFlowAuxVarCopy, &
            ZFlowAuxVarDestroy, &
            ZFlowAuxVarStrip, &
            ZFlowAuxVarPerturb, &
            ZFlowPrintAuxVars, &
            ZFlowOutputAuxVars

contains

! ************************************************************************** !

function ZFlowAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Option_module

  implicit none

  type(option_type) :: option

  type(zflow_type), pointer :: ZFlowAuxCreate

  type(zflow_type), pointer :: aux

  allocate(aux)
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%matrix_zeroing)

  allocate(aux%zflow_parameter)
  aux%zflow_parameter%check_post_converged = PETSC_FALSE

  ZFlowAuxCreate => aux

end function ZFlowAuxCreate

! ************************************************************************** !

subroutine ZFlowAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Option_module

  implicit none

  type(zflow_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%pc = 0.d0
  auxvar%kr = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%dpor_dp = 0.d0
  auxvar%dsat_dp = 0.d0
  auxvar%dkr_dp = 0.d0
  auxvar%pert = 0.d0

end subroutine ZFlowAuxVarInit

! ************************************************************************** !

subroutine ZFlowAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Option_module

  implicit none

  type(zflow_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%sat = auxvar%sat
  auxvar2%pc = auxvar%pc
  auxvar2%kr = auxvar%kr
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%dpor_dp = auxvar%dpor_dp
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dkr_dp = auxvar%dkr_dp
  auxvar2%pert = auxvar%pert

end subroutine ZFlowAuxVarCopy

! ************************************************************************** !

subroutine ZFlowAuxVarCompute(x,zflow_auxvar,global_auxvar, &
                              material_auxvar,characteristic_curves, &
                              natural_id,update_porosity,option)
  !
  ! Computes auxiliary variables for each grid cell
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Option_module
  use Global_Aux_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use Variables_module, only : SOIL_REFERENCE_PRESSURE

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(1)
  type(zflow_auxvar_type) :: zflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscBool :: update_porosity
  PetscInt :: natural_id

  PetscBool :: saturated
  PetscReal :: dkr_dsat

  zflow_auxvar%pres = x(ZFLOW_LIQUID_PRESSURE_DOF)
  global_auxvar%temp = option%flow%reference_temperature

  if (update_porosity .and. soil_compressibility_index > 0) then
    call MaterialCompressSoil(material_auxvar,zflow_auxvar%pres, &
                              zflow_auxvar%effective_porosity, &
                              zflow_auxvar%dpor_dp)
  else
    zflow_auxvar%effective_porosity = material_auxvar%porosity
    zflow_auxvar%dpor_dp = 0.d0
  endif
!  if (option%iflag /= ZFLOW_UPDATE_FOR_DERIVATIVE) then
!    material_auxvar%porosity = zflow_auxvar%effective_porosity
!  endif

  ! For a very large negative liquid pressure (e.g. -1.d18), the capillary
  ! pressure can go near infinite, resulting in ds_dp being < 1.d-40 below
  ! and flipping the cell to saturated, when it is really far from saturated.
  ! The large negative liquid pressure is then passed to the EOS causing it
  ! to blow up.  Therefore, we truncate to the max capillary pressure here.
  zflow_auxvar%pc = min(option%flow%reference_pressure - zflow_auxvar%pres, &
                        characteristic_curves%saturation_function%pcmax)

  if (zflow_auxvar%pc > 0.d0) then
    saturated = PETSC_FALSE
    call characteristic_curves%saturation_function% &
                               Saturation(zflow_auxvar%pc, &
                                          zflow_auxvar%sat, &
                                          zflow_auxvar%dsat_dp,option)
    ! if ds_dp is 0, we consider the cell saturated.
    if (zflow_auxvar%dsat_dp < 1.d-40) then
      saturated = PETSC_TRUE
    else
      dkr_dsat = 0.d0
      call characteristic_curves%liq_rel_perm_function% &
                       RelativePermeability(zflow_auxvar%sat, &
                                            zflow_auxvar%kr, &
                                            dkr_dsat,option)
      zflow_auxvar%dkr_dp = zflow_auxvar%dsat_dp * dkr_dsat
    endif
  else
    saturated = PETSC_TRUE
  endif

  ! the purpose for splitting this condition from the 'else' statement
  ! above is due to SaturationFunctionCompute switching a cell to
  ! saturated to prevent unstable (potentially infinite) derivatives when
  ! capillary pressure is very small
  if (saturated) then
    zflow_auxvar%pc = 0.d0
    zflow_auxvar%sat = 1.d0
    zflow_auxvar%kr = 1.d0
    zflow_auxvar%dsat_dp = 0.d0
    zflow_auxvar%dkr_dp = 0.d0
  endif

  if (option%iflag /= ZFLOW_UPDATE_FOR_DERIVATIVE) then
    if (size(global_auxvar%sat) > 1) then
      global_auxvar%sat(1) = zflow_auxvar%sat
      global_auxvar%sat(2) = 1.d0 - global_auxvar%sat(1)
    endif
  endif

end subroutine ZFlowAuxVarCompute

! ************************************************************************** !

subroutine ZFlowAuxVarPerturb(zflow_auxvar,global_auxvar, &
                              material_auxvar, &
                              characteristic_curves,natural_id, &
                              option)
  ! Calculates auxiliary variables for perturbed system
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(zflow_auxvar_type) :: zflow_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal :: x, x_pert(1), pert

  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-10

  ! ZFLOW_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = ZFLOW_UPDATE_FOR_DERIVATIVE
  x = zflow_auxvar(ZERO_INTEGER)%pres
  pert = x*zflow_pres_rel_pert+zflow_pres_min_pert
  zflow_auxvar(1)%pert = pert
  x_pert(1) = x + pert
  call ZFlowAuxVarCompute(x_pert,zflow_auxvar(ONE_INTEGER),global_auxvar, &
                          material_auxvar, &
                          characteristic_curves,natural_id, &
                          PETSC_TRUE,option)

end subroutine ZFlowAuxVarPerturb

! ************************************************************************** !

subroutine ZFlowPrintAuxVars(zflow_auxvar,global_auxvar,material_auxvar, &
                             natural_id,string,option)
  !
  ! Prints out the contents of an auxvar
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(zflow_auxvar_type) :: zflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  print *, '         liquid pressure: ', zflow_auxvar%pres
  print *, '      capillary pressure: ', zflow_auxvar%pc
  print *, '       liquid saturation: ', zflow_auxvar%sat
  print *, '      liquid sat (deriv): ', zflow_auxvar%dsat_dp
  print *, '         liquid rel perm: ', zflow_auxvar%kr
  print *, ' liquid rel perm (deriv): ', zflow_auxvar%dkr_dp
  print *, '      effective porosity: ', zflow_auxvar%effective_porosity
  print *, '   eff. porosity (deriv): ', zflow_auxvar%dpor_dp
  print *, '--------------------------------------------------------'

end subroutine ZFlowPrintAuxVars

! ************************************************************************** !

subroutine ZFlowOutputAuxVars1(zflow_auxvar,global_auxvar,material_auxvar, &
                               natural_id,string,append,option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(zflow_auxvar_type) :: zflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  PetscBool :: append
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string2

  write(string2,*) natural_id
  string2 = trim(adjustl(string)) // '_' // trim(adjustl(string2)) // '.txt'
  if (append) then
    open(unit=IUNIT_TEMP,file=string2,position='append')
  else
    open(unit=IUNIT_TEMP,file=string2)
  endif

  write(IUNIT_TEMP,*) '--------------------------------------------------------'
  write(IUNIT_TEMP,*) trim(string)
  write(IUNIT_TEMP,*) '             cell id: ', natural_id
  write(IUNIT_TEMP,*) '                 cell id: ', natural_id
  write(IUNIT_TEMP,*) '         liquid pressure: ', zflow_auxvar%pres
  write(IUNIT_TEMP,*) '      capillary pressure: ', zflow_auxvar%pc
  write(IUNIT_TEMP,*) '       liquid saturation: ', zflow_auxvar%sat
  write(IUNIT_TEMP,*) '      liquid sat (deriv): ', zflow_auxvar%dsat_dp
  write(IUNIT_TEMP,*) '         liquid rel perm: ', zflow_auxvar%kr
  write(IUNIT_TEMP,*) ' liquid rel perm (deriv): ', zflow_auxvar%dkr_dp
  write(IUNIT_TEMP,*) '      effective porosity: ', zflow_auxvar%effective_porosity
  write(IUNIT_TEMP,*) '   eff. porosity (deriv): ', zflow_auxvar%dpor_dp
  write(IUNIT_TEMP,*) '--------------------------------------------------------'

  close(IUNIT_TEMP)

end subroutine ZFlowOutputAuxVars1

! ************************************************************************** !

subroutine ZFlowAuxVarSingleDestroy(auxvar)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  implicit none

  type(zflow_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call ZFlowAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine ZFlowAuxVarSingleDestroy

! ************************************************************************** !

subroutine ZFlowAuxVarArray1Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  implicit none

  type(zflow_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call ZFlowAuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine ZFlowAuxVarArray1Destroy

! ************************************************************************** !

subroutine ZFlowAuxVarArray2Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  implicit none

  type(zflow_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call ZFlowAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine ZFlowAuxVarArray2Destroy

! ************************************************************************** !

subroutine ZFlowAuxVarStrip(auxvar)
  !
  ! ZFlowAuxVarDestroy: Deallocates a general auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(zflow_auxvar_type) :: auxvar

end subroutine ZFlowAuxVarStrip

! ************************************************************************** !

subroutine ZFlowAuxDestroy(aux)
  !
  ! Deallocates a general auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(zflow_type), pointer :: aux

  if (.not.associated(aux)) return

  call ZFlowAuxVarDestroy(aux%auxvars)
  call ZFlowAuxVarDestroy(aux%auxvars_bc)
  call ZFlowAuxVarDestroy(aux%auxvars_ss)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%zflow_parameter)) then
  endif
  nullify(aux%zflow_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine ZFlowAuxDestroy

end module ZFlow_Aux_module
