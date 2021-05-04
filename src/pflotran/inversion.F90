module Inversion_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none

  private
    
  type, public :: inversion_type
    PetscInt :: miniter,maxiter      ! min/max CGLS iterations
    
    PetscReal :: beta                ! regularization parameter
    PetscReal :: min_beta_red        ! minimum beta reduction 
    PetscReal :: mincond,maxcond     ! min/max conductivity
    PetscReal :: target_chi2         ! target CHI^2 norm 

    ! arrays for CGLS algorithm
    PetscReal, pointer :: p(:)       ! vector of dim -> num of inv cells
    PetscReal, pointer :: q(:)       ! product of Jacobian with p = Jp
    PetscReal, pointer :: r(:)       ! vector of dim -> num of measurements
    PetscReal, pointer :: s(:)       ! product of Jacobian transpose with r
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

  allocate(inversion)

  ! Default inversion parameters
  inversion%miniter = 10
  inversion%maxiter = 50

  inversion%beta = 100.d0
  inversion%min_beta_red = 0.5d0
  inversion%mincond = 0.00001d0
  inversion%maxcond = 10.d0
  inversion%target_chi2 = 1.d0

  nullify(inversion%p)
  nullify(inversion%q)
  nullify(inversion%r) 
  nullify(inversion%s)

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
  allocate(inversion%q(survey%num_measurement))
  allocate(inversion%r(survey%num_measurement))        
  allocate(inversion%s(grid%nlmax))

  inversion%p = 0.d0
  inversion%q = 0.d0
  inversion%r = 0.d0
  inversion%s = 0.d0

end subroutine InversionInit

! ************************************************************************** !

subroutine InversionOptionRead(inversion)
  !
  ! Read inversion options
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/04/21
      
  implicit none
      
  type(inversion_type) :: inversion

  ! Here we read inversion options from input file
        
end subroutine InversionOptionRead

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
  call DeallocateArray(inversion%q)
  call DeallocateArray(inversion%r)
  call DeallocateArray(inversion%s)

  deallocate(inversion)
  nullify(inversion)

end subroutine InversionDestroy

end module Inversion_module
