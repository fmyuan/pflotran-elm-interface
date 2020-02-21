module Matrix_Zeroing_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  type, public :: matrix_zeroing_type
    PetscInt :: n_inactive_rows
    PetscInt, pointer :: inactive_rows_local(:)
    PetscInt, pointer :: inactive_rows_local_ghosted(:)
    PetscInt, pointer :: row_zeroing_array(:)
  end type matrix_zeroing_type
  
  public :: MatrixZeroingCreate, &
            MatrixZeroingInitInactive, &
            MatrixZeroingInitRowZeroing, &
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
  matrix_zeroing%n_inactive_rows = 0
  nullify(matrix_zeroing%inactive_rows_local)         ! inactive
  nullify(matrix_zeroing%inactive_rows_local_ghosted) ! inactive
  nullify(matrix_zeroing%row_zeroing_array)           ! disabled (e.g. isotherm)
  
  MatrixZeroingCreate => matrix_zeroing
  
end function MatrixZeroingCreate

! ************************************************************************** !

subroutine MatrixZeroingInitInactive(matrix_zeroing,n_inactive_rows)
  ! 
  ! Initialize inactive arrays in matrix zeroing object
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/19
  
  use Utility_module, only : DeallocateArray

  implicit none

  type(matrix_zeroing_type), pointer :: matrix_zeroing
  PetscInt :: n_inactive_rows

  if (.not.associated(matrix_zeroing)) matrix_zeroing => MatrixZeroingCreate()
  call DeallocateArray(matrix_zeroing%inactive_rows_local)
  call DeallocateArray(matrix_zeroing%inactive_rows_local_ghosted)

  matrix_zeroing%n_inactive_rows = n_inactive_rows
  allocate(matrix_zeroing%inactive_rows_local(n_inactive_rows))
  matrix_zeroing%inactive_rows_local = 0
  allocate(matrix_zeroing%inactive_rows_local_ghosted(n_inactive_rows))
  matrix_zeroing%inactive_rows_local_ghosted = 0
  
end subroutine MatrixZeroingInitInactive

! ************************************************************************** !

subroutine MatrixZeroingInitRowZeroing(matrix_zeroing,nrows)
  ! 
  ! Initialize inactive arrays in matrix zeroing object
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/19

  use Utility_module, only : DeallocateArray

  implicit none

  type(matrix_zeroing_type), pointer :: matrix_zeroing
  PetscInt :: nrows

  if (.not.associated(matrix_zeroing)) matrix_zeroing => MatrixZeroingCreate()
  call DeallocateArray(matrix_zeroing%row_zeroing_array)

  allocate(matrix_zeroing%row_zeroing_array(nrows))
  matrix_zeroing%row_zeroing_array = 0
  
end subroutine MatrixZeroingInitRowZeroing

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
  
  call DeallocateArray(matrix_zeroing%inactive_rows_local)
  call DeallocateArray(matrix_zeroing%inactive_rows_local_ghosted)
  call DeallocateArray(matrix_zeroing%row_zeroing_array)
  
  deallocate(matrix_zeroing)
  nullify(matrix_zeroing)
  
end subroutine MatrixZeroingDestroy

end module Matrix_Zeroing_module
