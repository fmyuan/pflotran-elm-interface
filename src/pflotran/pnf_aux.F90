module PNF_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none

  private

  PetscReal, parameter, public :: pnf_density_kg = 998.32d0
!  PetscReal, parameter, public :: pnf_viscosity = 8.9d-4
  PetscReal, parameter, public :: pnf_viscosity = 1.d-3

  ! debugging
  PetscInt, public :: pnf_ts_cut_count
  PetscInt, public :: pnf_ts_count

  PetscInt, parameter, public :: PNF_LIQUID_PRESSURE_DOF = 1

  PetscInt, parameter, public :: PNF_LIQUID_EQUATION_INDEX = 1

  PetscInt, parameter, public :: PNF_LIQUID_PRESSURE_INDEX = 1
  PetscInt, parameter, public :: PNF_LIQUID_FLUX_INDEX = 1
  PetscInt, parameter, public :: PNF_LIQUID_CONDUCTANCE_INDEX = 2
  PetscInt, parameter, public :: PNF_MAX_INDEX = 2

  type, public :: pnf_auxvar_type
    PetscReal :: head ! liquid head
  end type pnf_auxvar_type

  type, public :: pnf_parameter_type
!    PetscBool :: check_post_converged
  end type pnf_parameter_type

  type, public :: pnf_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(pnf_parameter_type), pointer :: pnf_parameter
    type(pnf_auxvar_type), pointer :: auxvars(:)
    type(pnf_auxvar_type), pointer :: auxvars_bc(:)
    type(pnf_auxvar_type), pointer :: auxvars_ss(:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type pnf_type

  interface PNFAuxVarDestroy
    module procedure PNFAuxVarSingleDestroy
    module procedure PNFAuxVarArray1Destroy
    module procedure PNFAuxVarArray2Destroy
  end interface PNFAuxVarDestroy

  public :: PNFAuxCreate, &
            PNFAuxDestroy, &
            PNFAuxVarCompute, &
            PNFAuxVarInit, &
            PNFAuxVarCopy, &
            PNFAuxVarDestroy, &
            PNFAuxVarStrip

contains

! ************************************************************************** !

function PNFAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use Option_module

  implicit none

  type(option_type) :: option

  type(pnf_type), pointer :: PNFAuxCreate

  type(pnf_type), pointer :: aux

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

  allocate(aux%pnf_parameter)
!  aux%pnf_parameter%check_post_converged = PETSC_FALSE

  PNFAuxCreate => aux

end function PNFAuxCreate

! ************************************************************************** !

subroutine PNFAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use Option_module

  implicit none

  type(pnf_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%head = 0.d0

end subroutine PNFAuxVarInit

! ************************************************************************** !

subroutine PNFAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use Option_module

  implicit none

  type(pnf_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%head = auxvar%head

end subroutine PNFAuxVarCopy

! ************************************************************************** !

subroutine PNFAuxVarCompute(x,pnf_auxvar,global_auxvar, &
                               natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use Option_module
  use Global_Aux_module
  use Characteristic_Curves_module
  use Material_Aux_module
  use Variables_module, only : SOIL_REFERENCE_PRESSURE

  implicit none

  type(option_type) :: option
  PetscReal :: x(1)
  type(pnf_auxvar_type) :: pnf_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: natural_id

  PetscBool :: saturated
  PetscReal :: dkr_dsat

  pnf_auxvar%head = x(PNF_LIQUID_PRESSURE_DOF)
  global_auxvar%temp = option%flow%reference_temperature


end subroutine PNFAuxVarCompute

! ************************************************************************** !

subroutine PNFAuxVarSingleDestroy(auxvar)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  implicit none

  type(pnf_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call PNFAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine PNFAuxVarSingleDestroy

! ************************************************************************** !

subroutine PNFAuxVarArray1Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  implicit none

  type(pnf_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call PNFAuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine PNFAuxVarArray1Destroy

! ************************************************************************** !

subroutine PNFAuxVarArray2Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  implicit none

  type(pnf_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call PNFAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine PNFAuxVarArray2Destroy

! ************************************************************************** !

subroutine PNFAuxVarStrip(auxvar)
  !
  ! PNFAuxVarDestroy: Deallocates a general auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(pnf_auxvar_type) :: auxvar

end subroutine PNFAuxVarStrip

! ************************************************************************** !

subroutine PNFAuxDestroy(aux)
  !
  ! Deallocates a general auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(pnf_type), pointer :: aux

  if (.not.associated(aux)) return

  call PNFAuxVarDestroy(aux%auxvars)
  call PNFAuxVarDestroy(aux%auxvars_bc)
  call PNFAuxVarDestroy(aux%auxvars_ss)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%pnf_parameter)) then
  endif
  nullify(aux%pnf_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine PNFAuxDestroy

end module PNF_Aux_module
