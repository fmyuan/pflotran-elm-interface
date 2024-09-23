module Matrix_Zeroing_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: matrix_zeroing_type
    PetscBool :: zero_rows_exist
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:)         ! 1-based indexing
    ! zero_rows_local /= zero_rows_local_ghosted + 1
    PetscInt, pointer :: zero_rows_local_ghosted(:) ! 0-based indexing
  end type matrix_zeroing_type

  public :: MatrixZeroingCreate, &
            MatrixZeroingAllocateArray, &
            MatrixZeroingZeroVecEntries, &
            MatrixZeroingZeroArrayEntries, &
            MatrixZeroingZeroMatEntries, &
            MatrixZeroingDestroy

contains

! ************************************************************************** !

function MatrixZeroingCreate()
  !
  ! MatrixZeroingCreate: Allocate and initialize zeroing object
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/19
  !
  implicit none

  type(matrix_zeroing_type), pointer :: matrix_zeroing

  type(matrix_zeroing_type), pointer :: MatrixZeroingCreate

  allocate(matrix_zeroing)
  matrix_zeroing%zero_rows_exist = PETSC_FALSE
  matrix_zeroing%n_zero_rows = 0
  nullify(matrix_zeroing%zero_rows_local)
  nullify(matrix_zeroing%zero_rows_local_ghosted)

  MatrixZeroingCreate => matrix_zeroing

end function MatrixZeroingCreate

! ************************************************************************** !

subroutine MatrixZeroingAllocateArray(matrix_zeroing,n_zero_rows)
  !
  ! Initialize zero arrays in matrix zeroing object
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/19

  use Utility_module, only : DeallocateArray

  implicit none

  type(matrix_zeroing_type), pointer :: matrix_zeroing
  PetscInt :: n_zero_rows

  if (.not.associated(matrix_zeroing)) matrix_zeroing => MatrixZeroingCreate()
  call DeallocateArray(matrix_zeroing%zero_rows_local)
  call DeallocateArray(matrix_zeroing%zero_rows_local_ghosted)

  matrix_zeroing%n_zero_rows = n_zero_rows
  allocate(matrix_zeroing%zero_rows_local(n_zero_rows))
  matrix_zeroing%zero_rows_local = 0
  allocate(matrix_zeroing%zero_rows_local_ghosted(n_zero_rows))
  matrix_zeroing%zero_rows_local_ghosted = 0

end subroutine MatrixZeroingAllocateArray

! ************************************************************************** !

subroutine MatrixZeroingZeroVecEntries(matrix_zeroing,vec)
  !
  ! Zeros select entries in PETSc Vec
  !
  ! Author: Glenn Hammond
  ! Date: 09/19/24

  implicit none

  type(matrix_zeroing_type) :: matrix_zeroing
  Vec :: vec

  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  if (.not.matrix_zeroing%zero_rows_exist) return

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  call MatrixZeroingZeroArrayEntries(matrix_zeroing,vec_ptr)
  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine MatrixZeroingZeroVecEntries

! ************************************************************************** !

subroutine MatrixZeroingZeroArrayEntries(matrix_zeroing,array)
  !
  ! Zeros select entries in real array
  !
  ! Author: Glenn Hammond
  ! Date: 09/19/24

  implicit none

  type(matrix_zeroing_type) :: matrix_zeroing
  PetscReal :: array(:)
  PetscInt :: i

  do i=1,matrix_zeroing%n_zero_rows
    array(matrix_zeroing%zero_rows_local(i)) = 0.d0
  enddo

end subroutine MatrixZeroingZeroArrayEntries

! ************************************************************************** !

subroutine MatrixZeroingZeroMatEntries(matrix_zeroing,mat)
  !
  ! Zeros select entries in real array
  !
  ! Author: Glenn Hammond
  ! Date: 09/19/24

  implicit none

  type(matrix_zeroing_type) :: matrix_zeroing
  Mat :: mat

  PetscReal, parameter :: diagonal_value = 1.d0
  PetscErrorCode :: ierr

  if (.not.matrix_zeroing%zero_rows_exist) return

  call MatZeroRowsLocal(mat,matrix_zeroing%n_zero_rows, &
                        matrix_zeroing%zero_rows_local_ghosted, &
                        diagonal_value,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                        ierr);CHKERRQ(ierr)

end subroutine MatrixZeroingZeroMatEntries

! ************************************************************************** !

subroutine MatrixZeroingDestroy(matrix_zeroing)
  !
  ! Deallocates a matrix zeroing object
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/19
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(matrix_zeroing_type), pointer :: matrix_zeroing

  if (.not.associated(matrix_zeroing)) return

  call DeallocateArray(matrix_zeroing%zero_rows_local)
  call DeallocateArray(matrix_zeroing%zero_rows_local_ghosted)

  deallocate(matrix_zeroing)
  nullify(matrix_zeroing)

end subroutine MatrixZeroingDestroy

end module Matrix_Zeroing_module
