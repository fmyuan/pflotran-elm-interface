module ERT_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none

  private

  type, public :: ert_auxvar_type
    PetscReal, pointer :: potential(:)    ! ERT potential for all electrodes
    PetscReal, pointer :: jacobian(:)     ! ERT jacobian for all measurements
    PetscReal, pointer :: delM(:)         ! system matrix derivative dM/dcond
  end type ert_auxvar_type

  type, public :: ert_type
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux
    ! ert auxvars for local and ghosted cells
    type(ert_auxvar_type), pointer :: auxvars(:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type ert_type

  public :: ERTAuxCreate, ERTAuxDestroy, &
            ERTAuxVarCompute, ERTAuxVarInit, &
            ERTAuxCheckElectrodeBounds, &
            ERTAuxVarCopy

  ! conductivity mapping types
  PetscInt, parameter, public :: ARCHIE = 1
  PetscInt, parameter, public :: WAXMAN_SMITS = 2

contains

! ************************************************************************** !

function ERTAuxCreate()
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  !

  use Option_module

  implicit none

  type(ert_type), pointer :: ERTAuxCreate

  type(ert_type), pointer :: aux

  allocate(aux)
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0

  nullify(aux%auxvars)
  nullify(aux%matrix_zeroing)

  ERTAuxCreate => aux

end function ERTAuxCreate

! ************************************************************************** !

subroutine ERTAuxVarInit(auxvar,survey,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  !

  use Option_module
  use Survey_module

  implicit none

  type(ert_auxvar_type) :: auxvar
  type(survey_type) :: survey
  type(option_type) :: option

  allocate(auxvar%potential(survey%num_electrode))
  auxvar%potential = 0.d0

  nullify(auxvar%jacobian)
  if (option%geophysics%compute_jacobian) then
    allocate(auxvar%jacobian(survey%num_measurement))
    auxvar%jacobian = 0.d0
  endif
  nullify(auxvar%delM)

end subroutine ERTAuxVarInit

! ************************************************************************** !

subroutine ERTAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  !

  use Option_module

  implicit none

  type(ert_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%potential = auxvar%potential

end subroutine ERTAuxVarCopy

! ************************************************************************** !

subroutine ERTAuxVarCompute(x,ert_auxvar,global_auxvar,rt_auxvar, &
                            material_auxvar,natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell
  ! Could be Archie's law
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  !
  use Option_module
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  type(ert_auxvar_type) :: ert_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reactive_transport_auxvar_type) :: rt_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id

  ! calculate bulk_conductivity = f(global,rt,material-auxars)

end subroutine ERTAuxVarCompute

! ************************************************************************** !

subroutine ERTAuxCheckElectrodeBounds(num_electrodes,id,id2,option)
  !
  ! Checks to ensure that electrode requests (ivar, ivar2) are valid
  ! electrode IDs.
  !
  ! Author: Glenn Hammond
  ! Date: 12/09/21
  !
  use Option_module
  use String_module

  type(option_type) :: option
  PetscInt :: num_electrodes
  PetscInt :: id, id2

  PetscInt :: imax

  imax = max(id,id2)
  if (imax > num_electrodes) then
    option%io_buffer = 'Electrode ID (' // trim(StringWrite(imax)) // &
      ') specified for output is greater than the maximum electrode ID (' // &
      trim(StringWrite(num_electrodes)) // ').'
    call PrintErrMsg(option)
  endif

end subroutine ERTAuxCheckElectrodeBounds

! ************************************************************************** !

subroutine ERTAuxVarDestroy(auxvar)
  !
  ! Deallocates an ert auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  !
  use Utility_module, only: DeallocateArray

  implicit none

  type(ert_auxvar_type) :: auxvar

  call DeallocateArray(auxvar%potential)
  call DeallocateArray(auxvar%jacobian)
  call DeallocateArray(auxvar%delM)

end subroutine ERTAuxVarDestroy

! ************************************************************************** !

subroutine ERTAuxDestroy(aux)
  !
  ! Deallocates an ert auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  !

  implicit none

  type(ert_type), pointer :: aux
  PetscInt :: iaux

  if (.not.associated(aux)) return

  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call ERTAuxVarDestroy(aux%auxvars(iaux))
    enddo
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  deallocate(aux)
  nullify(aux)

end subroutine ERTAuxDestroy

end module ERT_Aux_module
