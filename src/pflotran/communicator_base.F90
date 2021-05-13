module Communicator_Base_class

#include "petsc/finclude/petscvec.h"
   use petscvec

   use PFLOTRAN_Constants_module

  implicit none

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

  type, abstract, public :: communicator_type
  contains
    procedure(SetDM), public, deferred :: SetDM 
    procedure(VecToVec), public, deferred :: GlobalToLocal 
    procedure(VecToVec), public, deferred :: LocalToGlobal
    procedure(VecToVec), public, deferred :: LocalToLocal 
    procedure(VecToVec), public, deferred :: GlobalToNatural 
    procedure(VecToVec), public, deferred :: NaturalToGlobal 
    procedure(MapArray), public, deferred :: AONaturalToPetsc 
    procedure(BaseDestroy), public, deferred :: Destroy 
  end type communicator_type
  
  abstract interface
  
#ifdef SIMPLIFY    
    subroutine SetDM(this)
      import communicator_type
      implicit none
      class(communicator_type) :: this
#else
    subroutine SetDM(this,dm_ptr)
#include "petsc/finclude/petscdm.h"
      use petscdm
      use DM_Kludge_module
      import communicator_type
      implicit none
      class(communicator_type) :: this
      type(dm_ptr_type) :: dm_ptr
#endif    
    end subroutine
  
    subroutine VecToVec(this,source,destination)
#include "petsc/finclude/petscvec.h"
      use petscvec
      import communicator_type
      implicit none
      class(communicator_type) :: this
      Vec :: source
      Vec :: destination
    end subroutine VecToVec

    subroutine MapArray(this,array)
      import communicator_type
      implicit none
      class(communicator_type) :: this
      PetscInt :: array(:)
    end subroutine MapArray

    subroutine BaseDestroy(this)
      import communicator_type
      implicit none
      class(communicator_type) :: this
    end subroutine BaseDestroy

  end interface

  interface CommCreate
    module procedure CommCreate1
    module procedure CommCreate2
  end interface
  
  public :: CommCreate, &
            CommCreateProcessorGroups
  
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

subroutine CommCreateProcessorGroups(global_comm,global_commsize,global_rank, &
                                     mycomm,mycommsize,myrank, &
                                     mygroup,mygroup_id,num_groups)
  ! 
  ! Splits MPI_COMM_WORLD into N separate
  ! processor groups
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/09
  ! 
  implicit none
  
  PetscMPIInt :: global_comm
  PetscMPIInt :: global_commsize
  PetscMPIInt :: global_rank
  PetscMPIInt :: mycomm
  PetscMPIInt :: mycommsize
  PetscMPIInt :: myrank
  PetscMPIInt :: mygroup
  PetscMPIInt :: mygroup_id
  PetscInt :: num_groups

  PetscInt :: local_commsize
  PetscInt :: offset, delta, remainder
  PetscInt :: igroup
  PetscMPIInt :: mycolor_mpi, mykey_mpi
  PetscErrorCode :: ierr

  local_commsize = global_commsize / num_groups
  remainder = global_commsize - num_groups * local_commsize
  offset = 0
  do igroup = 1, num_groups
    delta = local_commsize
    if (igroup < remainder) delta = delta + 1
    if (global_rank >= offset .and. &
        global_rank < offset + delta) exit
    offset = offset + delta
  enddo
  mycolor_mpi = igroup
  mygroup_id = igroup
  mykey_mpi = global_rank - offset
  call MPI_Comm_split(MPI_COMM_WORLD,mycolor_mpi,mykey_mpi,mycomm,ierr)
  call MPI_Comm_group(mycomm,mygroup,ierr)

  call MPI_Comm_rank(mycomm,myrank, ierr)
  call MPI_Comm_size(mycomm,mycommsize,ierr)

end subroutine CommCreateProcessorGroups
  
end module Communicator_Base_class
