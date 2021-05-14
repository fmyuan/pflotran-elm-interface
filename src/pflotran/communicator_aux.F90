module Communicator_Aux_module

  use petscsys

  implicit none

#include "petsc/finclude/petscsys.h"

  private

  type, public :: comm_type
    PetscMPIInt :: global_comm      ! MPI_COMM_WORLD
    PetscMPIInt :: global_rank      ! rank in MPI_COMM_WORLD
    PetscMPIInt :: global_commsize  ! size of MPI_COMM_WORLD
    PetscMPIInt :: global_group     ! id of group in MPI_COMM_WORLD

    PetscMPIInt :: mycomm           ! PETSC_COMM_WORLD
    PetscMPIInt :: myrank           ! rank in PETSC_COMM_WORLD
    PetscMPIInt :: mycommsize       ! size of PETSC_COMM_WORLD
    PetscMPIInt :: mygroup          ! id of group in PETSC_COMM_WORLD
    PetscMPIInt :: mygroup_id
  end type comm_type

  public :: CommCreate, &
            CommInit, &
            CommCreateProcessorGroups, &
            CommDestroy
  
contains

! ************************************************************************** !

function CommCreate()
  !
  ! Creates a comm object that holds global and local communicators, sizes 
  ! and ranks
  !
  ! Author: Glenn Hammond
  ! Date: 05/12/21
  
  implicit none

  type(comm_type), pointer :: comm
  type(comm_type), pointer :: CommCreate

  allocate(comm)
  comm%global_comm = 0
  comm%global_rank = 0
  comm%global_commsize = 0
  comm%global_group = 0

  comm%mycomm = 0
  comm%myrank = 0
  comm%mycommsize = 0
  comm%mygroup = 0

  comm%mygroup_id = 0

  CommCreate => comm

end function CommCreate

! ************************************************************************** !

subroutine CommInit(comm,communicator)

  implicit none

  type(comm_type) :: comm

  PetscMPIInt, optional :: communicator

  PetscErrorCode :: ierr

  if (present(communicator)) PETSC_COMM_WORLD = communicator
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)

  comm%global_comm = PETSC_COMM_WORLD
  call MPI_Comm_rank(comm%global_comm,comm%global_rank, ierr)
  call MPI_Comm_size(comm%global_comm,comm%global_commsize,ierr)
  call MPI_Comm_group(comm%global_comm,comm%global_group,ierr)
  comm%mycomm = comm%global_comm
  comm%myrank = comm%global_rank
  comm%mycommsize = comm%global_commsize
  comm%mygroup = comm%global_group

end subroutine CommInit

! ************************************************************************** !

subroutine CommCreateProcessorGroups(comm,num_groups)
  ! 
  ! Splits MPI_COMM_WORLD into N separate
  ! processor groups
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/09
  ! 
  implicit none

  type(comm_type) :: comm
  PetscInt :: num_groups

  PetscInt :: local_commsize
  PetscInt :: offset, delta, remainder
  PetscInt :: igroup
  PetscMPIInt :: mycolor_mpi, mykey_mpi
  PetscErrorCode :: ierr

  local_commsize = comm%global_commsize / num_groups
  remainder = comm%global_commsize - num_groups * local_commsize
  offset = 0
  do igroup = 1, num_groups
    delta = local_commsize
    if (igroup < remainder) delta = delta + 1
    if (comm%global_rank >= offset .and. &
        comm%global_rank < offset + delta) exit
    offset = offset + delta
  enddo
  mycolor_mpi = igroup
  comm%mygroup_id = igroup
  mykey_mpi = comm%global_rank - offset
  call MPI_Comm_split(comm%global_comm,mycolor_mpi,mykey_mpi,comm%mycomm,ierr)
  call MPI_Comm_group(comm%mycomm,comm%mygroup,ierr)

  call MPI_Comm_rank(comm%mycomm,comm%myrank, ierr)
  call MPI_Comm_size(comm%mycomm,comm%mycommsize,ierr)

end subroutine CommCreateProcessorGroups

! ************************************************************************** !

subroutine CommDestroy(comm)

  implicit none

  type(comm_type), pointer :: comm

  PetscErrorCode :: ierr

  deallocate(comm)
  nullify(comm)
  call PetscFinalize(ierr);CHKERRQ(ierr)

end subroutine CommDestroy

end module Communicator_Aux_module
