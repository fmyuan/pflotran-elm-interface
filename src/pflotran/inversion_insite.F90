module Inversion_INSITE_class

#include "petsc/finclude/petscvec.h"
  use petscvec

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Realization_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_base_type) :: inversion_insite_type
    class(realization_subsurface_type), pointer :: realization
    Vec :: quantity_of_interest
    PetscInt :: iqoi
  contains
    procedure, public :: Init => InversionINSITEInit
    procedure, public :: Initialize => InversionINSITEInitialize
    procedure, public :: UpdateParameters => InversionINSITEUpdateParameters
    procedure, public :: CalculateInverse => InversionINSITECalculateInverse
    procedure, public :: Finalize => InversionINSITEFinalize
    procedure, public :: Strip => InversionINSITEStrip
  end type inversion_insite_type

  public :: InversionINSITECreate, &
            InversionINSITEStrip, &
            InversionINSITEDestroy

contains

! ************************************************************************** !

function InversionINSITECreate(driver)
  !
  ! Allocates and initializes a new inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Driver_module

  class(driver_type), pointer :: driver

  class(inversion_insite_type), pointer :: InversionINSITECreate

  allocate(InversionINSITECreate)
  call InversionINSITECreate%Init(driver)

end function InversionINSITECreate

! ************************************************************************** !

subroutine InversionINSITEInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, PERMEABILITY
  use Driver_module

  class(inversion_insite_type) :: this
  class(driver_type), pointer :: driver

  this%quantity_of_interest = PETSC_NULL_VEC
  this%iqoi = PERMEABILITY
  nullify(this%realization)
  call InversionBaseInit(this,driver)

end subroutine InversionINSITEInit

! ************************************************************************** !

subroutine InversionINSITEInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Discretization_module
  use Material_module

  class(inversion_insite_type) :: this

  PetscErrorCode :: ierr

  ! non-ghosted Vec
  call VecDuplicate(this%realization%field%work, &
                    this%quantity_of_interest,ierr);CHKERRQ(ierr)
  ! ghosted Vec
  call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                               this%realization%field%work_loc, &
                               this%iqoi,ZERO_INTEGER)
  call DiscretizationLocalToGlobal(this%realization%discretization, &
                                   this%realization%field%work_loc, &
                                   this%quantity_of_interest,ONEDOF)

end subroutine InversionINSITEInitialize

! ************************************************************************** !

subroutine InversionINSITEUpdateParameters(this)
  !
  ! Updates input parameters
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Discretization_module
  use Material_module

  class(inversion_insite_type) :: this

  if (this%quantity_of_interest == PETSC_NULL_VEC) then
    call this%Initialize()
  else
    call DiscretizationGlobalToLocal(this%realization%discretization, &
                                     this%quantity_of_interest, &
                                     this%realization%field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi, &
                                 ZERO_INTEGER)
  endif

end subroutine InversionINSITEUpdateParameters

! ************************************************************************** !

subroutine InversionINSITECalculateInverse(this)
  !
  ! Calculates update to input parameters
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_insite_type) :: this

  PetscErrorCode :: ierr

  if (this%quantity_of_interest /= PETSC_NULL_VEC) then
    call VecScale(this%quantity_of_interest,2.d0,ierr);CHKERRQ(ierr)
  endif

end subroutine InversionINSITECalculateInverse

! ************************************************************************** !

subroutine InversionINSITEFinalize(this)
  !
  ! Finalizes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_insite_type) :: this

  call InversionBaseFinalize(this)

end subroutine InversionINSITEFinalize

! ************************************************************************** !

subroutine InversionINSITEStrip(this)
  !
  ! Deallocates members of inversion insite
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_insite_type) :: this

  PetscErrorCode :: ierr

  call InversionBaseStrip(this)

  nullify(this%realization)
  if (this%quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%quantity_of_interest,ierr);CHKERRQ(ierr)
  endif

end subroutine InversionINSITEStrip

! ************************************************************************** !

subroutine InversionINSITEDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_insite_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionINSITEDestroy

end module Inversion_INSITE_class
