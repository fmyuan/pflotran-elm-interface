module Communicator_Aux_module

  use petscsys

  implicit none

#include "petsc/finclude/petscsys.h"

  private

  type, public :: comm_type
    PetscMPIInt :: global_comm
    PetscMPIInt :: global_rank
    PetscMPIInt :: global_commsize
    PetscMPIInt :: global_group

    PetscMPIInt :: mycomm
    PetscMPIInt :: myrank
    PetscMPIInt :: mycommsize
    PetscMPIInt :: mygroup

    PetscMPIInt :: mygroup_id
  end type comm_type

  interface CommCreate
    module procedure CommCreate1
    module procedure CommCreate2
  end interface

  public :: CommCreate, &
            CommInit, &
            CommCreateProcessorGroups, &
            CommDestroy
  
contains

! ************************************************************************** !

function CommCreate1()
  !
  ! Creates a comm object that holds global and local communicators, sizes 
  ! and ranks
  !
  ! Author: Glenn Hammond
  ! Date: 05/12/21
  
  implicit none

  type(comm_type), pointer :: comm
  type(comm_type), pointer :: CommCreate1

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

  CommCreate1 => comm

end function CommCreate1

! ************************************************************************** !

function CommCreate2(orig_comm)
  !
  ! Creates a comm object that holds global and local communicators, sizes 
  ! and ranks
  !
  ! Author: Glenn Hammond
  ! Date: 05/12/21
  
  implicit none

  type(comm_type) :: orig_comm

  type(comm_type), pointer :: newcomm
  type(comm_type), pointer :: CommCreate2

  allocate(newcomm)
  newcomm%global_comm = orig_comm%global_comm
  newcomm%global_rank = orig_comm%global_rank
  newcomm%global_commsize = orig_comm%global_commsize
  newcomm%global_group = orig_comm%global_group

  newcomm%mycomm = orig_comm%mycomm
  newcomm%myrank = orig_comm%myrank
  newcomm%mycommsize = orig_comm%mycommsize
  newcomm%mygroup = orig_comm%mygroup

  newcomm%mygroup_id = orig_comm%mygroup_id

  CommCreate2 => newcomm

end function CommCreate2

! ************************************************************************** !

subroutine CommInit(comm)

  implicit none

  type(comm_type) :: comm

  PetscMPIInt :: communicator
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)
  communicator = PETSC_COMM_WORLD

  comm%global_comm = communicator
  call MPI_Comm_rank(communicator,comm%global_rank, ierr)
  call MPI_Comm_size(communicator,comm%global_commsize,ierr)
  call MPI_Comm_group(communicator,comm%global_group,ierr)
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

  deallocate(comm)
  nullify(comm)

end subroutine CommDestroy

end module Communicator_Aux_module
