module Inversion_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none

  private
    
  type, public :: inversion_type
    PetscReal, pointer :: p(:)       ! vector of dim -> num of inv cells
    PetscReal, pointer :: r(:)       ! vector of dim -> num of measurements
    PetscReal, pointer :: Jp(:)      ! Product of Jacobian with p
    PetscReal, pointer :: Jtr(:)     ! Product of Jacobian transpose with r
  end type inversion_type

contains 
  
  ! ************************************************************************** !

    function InversionCreate()
        !
        ! Creates inversion type
        !
        ! Author: Piyoosh Jaysaval
        ! Date: 05/03/21
        !
        implicit none
      
        type(inversion_type), pointer :: InversionCreate
        type(inversion_type), pointer :: inversion

        nullify(inversion%p)
        nullify(inversion%r)
        nullify(inversion%Jp)
        nullify(inversion%Jtr)

        InversionCreate => inversion

    end function InversionCreate

! ************************************************************************** !

    subroutine InversionInit(inversion,survey,grid)
        !
        ! Initialize inversion object
        !
        ! Author: Piyoosh Jaysaval
        ! Date: 05/03/21
        !

        use Survey_module
        use Grid_module
      
        implicit none
      
        type(inversion_type) :: inversion
        type(survey_type) :: survey
        type(grid_type), pointer :: grid
        
        allocate(inversion%p(grid%nlmax))
        allocate(inversion%r(survey%num_measurement))
        allocate(inversion%Jp(survey%num_measurement))
        allocate(inversion%Jtr(grid%nlmax))

        inversion%p = 0.d0
        inversion%r = 0.d0
        inversion%Jp = 0.d0
        inversion%Jtr = 0.d0

    end subroutine InversionInit

! ************************************************************************** !

    subroutine InversionDestroy(inversion)
        !
        ! Deallocates inversion
        !
        ! Author: Piyoosh Jaysaval
        ! Date: 05/03/21
        !
        use Utility_module, only : DeallocateArray
      
        implicit none
      
        type(inversion_type), pointer :: inversion
      
        if (.not.associated(inversion)) return
        
        call DeallocateArray(inversion%p)
        call DeallocateArray(inversion%r)
        call DeallocateArray(inversion%Jp)
        call DeallocateArray(inversion%Jtr)

        deallocate(inversion)
        nullify(inversion)

    end subroutine InversionDestroy

end module Inversion_module
