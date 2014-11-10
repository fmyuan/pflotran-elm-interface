module EPICData

  use EPICConstants

  implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  private

  type, public :: epic_data_type

     character(len=EPIC_MAXWORDLENGTH)    :: name              ! Data name
     character(len=EPIC_MAXWORDLENGTH)    :: units             ! units

     Vec                     :: data              ! PETSc vectors storing the data
     PetscInt                :: nsize             ! local size of PETSc vector
     PetscInt                :: var_id            ! hksat_x/hksat_y/….
     PetscBool               :: is_seq            ! Seq or MPI PETSc vector

    type(epic_data_type), pointer :: next

  end type epic_data_type

  type, public :: epic_data_list_type
    PetscInt                         :: ndata
    type(epic_data_type), pointer    :: first
    type(epic_data_type), pointer    :: last
  end type epic_data_list_type

  interface EPICDataCreate
    module procedure EPICDataCreate1
    module procedure EPICDataCreate2
  end interface EPICDataCreate

  interface EPICDataAddToList
    module procedure EPICDataAddToList1
    module procedure EPICDataAddToList2
  end interface EPICDataAddToList

  public :: EPICDataCreate
  public :: EPICDataDestroy
  public :: EPICDataListCreate
  public :: EPICDataAddToList
  public :: EPICDataGetPointer
  public :: EPICDataListDestroy

contains

! ************************************************************************** !

function EPICDataCreate1()

  implicit none
  
  type(epic_data_type),pointer :: EPICDataCreate1

  type(epic_data_type),pointer :: extphysidata

  allocate(extphysidata)
  extphysidata%name = ''
  extphysidata%units = ''

  extphysidata%data = 0
  extphysidata%nsize = 0
  extphysidata%var_id = EPIC_VAR_NULL
  extphysidata%is_seq = PETSC_FALSE
  nullify(extphysidata%next)

  EPICDataCreate1 => extphysidata

end function EPICDataCreate1

! ************************************************************************** !

function EPICDataCreate2(comm, name, units, nsize, var_id, is_seq)

  implicit none
  
  PetscMPIInt :: comm
  character(len=EPIC_MAXWORDLENGTH)    :: name
  character(len=EPIC_MAXWORDLENGTH)    :: units
  PetscInt :: nsize
  PetscInt :: var_id
  PetscBool :: is_seq

  type(epic_data_type),pointer :: EPICDataCreate2

  type(epic_data_type),pointer :: extphysidata
  PetscErrorCode :: ierr

  allocate(extphysidata)
  extphysidata%name = trim(adjustl(name))
  extphysidata%units = trim(adjustl(units))

  extphysidata%nsize = nsize
  extphysidata%var_id = var_id
  extphysidata%is_seq = is_seq

  if (is_seq) then
    call VecCreateSeq(PETSC_COMM_SELF, nsize, extphysidata%data, ierr)
  else
    call VecCreateMPI(comm, nsize, PETSC_DECIDE, extphysidata%data, ierr)
  endif

  nullify(extphysidata%next)

  EPICDataCreate2 => extphysidata

end function EPICDataCreate2

! ************************************************************************** !

function EPICDataListCreate()

  implicit none
  
  type(epic_data_list_type), pointer :: EPICDataListCreate
  
  type(epic_data_list_type), pointer :: data_list
  
  allocate(data_list)
  nullify(data_list%first)
  nullify(data_list%last)
  data_list%ndata = 0
  
  EPICDataListCreate => data_list
  
end function EPICDataListCreate

! ************************************************************************** !

subroutine EPICDataAddToList1(list, data)

  implicit none
  
  type(epic_data_list_type) :: list
  type(epic_data_type), pointer :: data
  
  if (.not. associated(list%first)) then
    list%first => data
  else
    list%last%next => data
  endif
  list%last => data
  
  list%ndata = list%ndata + 1
  
end subroutine EPICDataAddToList1

! ************************************************************************** !

subroutine EPICDataAddToList2(list, comm, name, units, nsize, var_id, &
                              is_seq)

  implicit none
  
  type(epic_data_list_type) :: list
  PetscMPIInt :: comm
  character(len=EPIC_MAXWORDLENGTH)    :: name
  character(len=EPIC_MAXWORDLENGTH)    :: units
  PetscInt :: nsize
  PetscInt :: var_id
  PetscBool :: is_seq

  type(epic_data_type), pointer :: data

  data => EPICDataCreate(comm, name, units, nsize, var_id, is_seq)
  call EPICDataAddToList1(list, data)

end subroutine EPICDataAddToList2

! ************************************************************************** !

function EPICDataGetPointer(list, data_type)

  implicit none

  type(epic_data_list_type), pointer :: list
  PetscInt :: data_type

  type(epic_data_type),pointer :: curr_epic_data
  Vec, pointer :: EPICDataGetPointer
  PetscBool :: found

  found = PETSC_FALSE
  curr_epic_data => list%first
  do 
    if (.not.associated(curr_epic_data)) exit
    if (curr_epic_data%var_id == data_type) then
      found = PETSC_TRUE
      exit
    endif
    curr_epic_data => curr_epic_data%next
  enddo

  EPICDataGetPointer => curr_epic_data%data

end function EPICDataGetPointer

! ************************************************************************** !

subroutine EPICDataListDestroy(list)

  implicit none
  
  type(epic_data_list_type), pointer :: list
  
  nullify(list%last)
  call EPICDataDestroy(list%first)
  
  deallocate(list)
  nullify(list)
  
end subroutine EPICDataListDestroy

! ************************************************************************** !

recursive subroutine EPICDataDestroy(epic_data)

  implicit none
  
  type(epic_data_type), pointer :: epic_data
  
  PetscErrorCode :: ierr

  if (.not.associated(epic_data)) return

  call VecDestroy(epic_data%data,ierr)
  epic_data%data = 0
  epic_data%nsize = 0
  epic_data%var_id = EPIC_VAR_NULL
  epic_data%is_seq = PETSC_FALSE

  call EPICDataDestroy(epic_data%next)  

  deallocate(epic_data)
  nullify(epic_data)
  
end subroutine EPICDataDestroy

end module EPICData
