module Communicator_Aux_module

  use petscsys

  implicit none

#include "petsc/finclude/petscsys.h"

  private

  type, public :: comm_type
    PetscMPIInt :: communicator ! communicator
    PetscMPIInt :: rank         ! rank in communicator
    PetscMPIInt :: size         ! size of communicator
    PetscMPIInt :: group        ! id of group in communicator

    PetscLogDouble :: start_time
    PetscInt :: group_id     ! PFLOTRAN centric group id
    PetscMPIInt :: io_rank

    type(comm_type), pointer :: parent
  end type comm_type

  interface CommIsIORank
    module procedure :: CommIsIORank1
    module procedure :: CommIsIORank2
  end interface

  public :: CommCreate, &
            CommInitPetsc, &
            CommResetStartTime, &
            CommPopulate, &
            CommIsIORank, &
            CommCreateProcessGroups, &
            CommStrip, &
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
  comm%communicator = MPI_COMM_NULL
  comm%rank = 0
  comm%size = 0
  comm%group = MPI_GROUP_NULL

  comm%start_time = 0.d0
  comm%group_id = 0
  comm%io_rank = 0

  nullify(comm%parent)

  CommCreate => comm

end function CommCreate

! ************************************************************************** !

subroutine CommInitPetsc(comm,communicator)

  implicit none

  type(comm_type), pointer :: comm
  PetscMPIInt, optional :: communicator

  PetscErrorCode :: ierr

  if (.not.associated(comm)) comm => CommCreate()
  if (present(communicator)) PETSC_COMM_WORLD = communicator
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)

  comm%communicator = PETSC_COMM_WORLD
  call CommResetStartTime(comm)
  call CommPopulate(comm)

end subroutine CommInitPetsc

! ************************************************************************** !

subroutine CommResetStartTime(comm)

  implicit none

  type(comm_type) :: comm

  PetscErrorCode :: ierr

  call PetscTime(comm%start_time,ierr);CHKERRQ(ierr)

end subroutine CommResetStartTime

! ************************************************************************** !

subroutine CommPopulate(comm)

  implicit none

  type(comm_type) comm

  PetscErrorCode :: ierr

  call MPI_Comm_rank(comm%communicator,comm%rank,ierr);CHKERRQ(ierr)
  call MPI_Comm_size(comm%communicator,comm%size,ierr);CHKERRQ(ierr)
  call MPI_Comm_group(comm%communicator,comm%group,ierr);CHKERRQ(ierr)

end subroutine CommPopulate

! ************************************************************************** !

function CommIsIORank1(comm)

  implicit none

  type(comm_type) comm

  PetscBool :: CommIsIORank1

  CommIsIORank1 = (comm%rank == comm%io_rank)

end function CommIsIORank1

! ************************************************************************** !

function CommIsIORank2(comm,rank)

  implicit none

  type(comm_type) comm
  PetscMPIInt :: rank

  PetscBool :: CommIsIORank2

  CommIsIORank2 = (rank == comm%io_rank)

end function CommIsIORank2

! ************************************************************************** !

subroutine CommCreateProcessGroups(parent_comm,num_groups, &
                                   force_equal_sized_groups,child_comm,ierr)
  !
  ! Splits an MPI communicator into N separate processor groups
  !
  ! Author: Glenn Hammond
  ! Date: 08/11/09
  !
  implicit none

  type(comm_type), target :: parent_comm
  PetscInt :: num_groups
  PetscBool :: force_equal_sized_groups
  type(comm_type), pointer :: child_comm
  PetscErrorCode :: ierr

  PetscInt :: local_commsize
  PetscInt :: offset, delta, remainder
  PetscInt :: igroup
  PetscMPIInt :: mycolor_mpi, mykey_mpi

  ierr = 0
  if (associated(child_comm)) then
    ierr = 1
    return
  endif

  local_commsize = parent_comm%size / num_groups
  remainder = parent_comm%size - num_groups * local_commsize
  if (force_equal_sized_groups .and. remainder > 0) then
    ierr = 1
    return
  endif
  offset = 0
  do igroup = 1, num_groups
    delta = local_commsize
    if (igroup < remainder) delta = delta + 1
    if (parent_comm%rank >= offset .and. &
        parent_comm%rank < offset + delta) exit
    offset = offset + delta
  enddo
  mycolor_mpi = igroup
  mykey_mpi = parent_comm%rank - offset
  child_comm => CommCreate()
  child_comm%group_id = igroup
  call MPI_Comm_split(parent_comm%communicator,mycolor_mpi,mykey_mpi, &
                      child_comm%communicator,ierr);CHKERRQ(ierr)
  call MPI_Comm_rank(child_comm%communicator,child_comm%rank, &
                     ierr);CHKERRQ(ierr)
  call MPI_Comm_size(child_comm%communicator,child_comm%size, &
                     ierr);CHKERRQ(ierr)
  call MPI_Comm_group(child_comm%communicator,child_comm%group, &
                      ierr);CHKERRQ(ierr)

  child_comm%parent => parent_comm

end subroutine CommCreateProcessGroups

! ************************************************************************** !

subroutine CommStrip(comm)

  ! deallocate MPI objects prior to PetscFinalize

  implicit none

  type(comm_type) :: comm

  PetscErrorCode :: ierr

  if (comm%group /= MPI_GROUP_NULL) then
    call MPI_Group_free(comm%group,ierr);CHKERRQ(ierr)
    comm%group = MPI_GROUP_NULL
  endif

  if (comm%communicator /= PETSC_COMM_WORLD .and. &
      comm%communicator /= MPI_COMM_NULL) then
    call MPI_Comm_free(comm%communicator,ierr);CHKERRQ(ierr)
    comm%communicator = MPI_COMM_NULL
  endif

end subroutine CommStrip

! ************************************************************************** !

subroutine CommDestroy(comm)

  implicit none

  type(comm_type), pointer :: comm

  if (.not.associated(comm)) return

  call CommStrip(comm)

  deallocate(comm)
  nullify(comm)

end subroutine CommDestroy

end module Communicator_Aux_module
