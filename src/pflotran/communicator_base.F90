module Communicator_Base_class

#include "petsc/finclude/petscvec.h"
   use petscvec

   use PFLOTRAN_Constants_module

  implicit none

  private

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

end module Communicator_Base_class
