module Inversion_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Survey_module

  implicit none

  private
    
  type, public :: inversion_type
    PetscInt :: iteration            ! iteration number
    PetscInt :: miniter,maxiter      ! min/max CGLS iterations
    
    PetscReal :: beta                ! regularization parameter
    PetscReal :: beta_red_factor     ! beta reduction factor
    PetscReal :: mincond,maxcond     ! min/max conductivity
    PetscReal :: target_chi2         ! target CHI^2 norm
    PetscReal :: current_chi2
    
    ! Cost/objective functions
    PetscReal :: min_phi_red         ! min change in cost function
    PetscReal :: phi_total_0,phi_total
    PetscReal :: phi_data_0,phi_data
    PetscReal :: phi_model_0,phi_model

    PetscBool :: cull_flag           ! flag to ignore data outliers
    PetscReal :: cull_dev            ! data culling cutoff (std. deviation)

    PetscBool :: converg_flag        ! convergence flag

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
  inversion%beta_red_factor = 0.5d0
  inversion%mincond = 0.00001d0
  inversion%maxcond = 10.d0
  inversion%target_chi2 = 1.d0
  inversion%min_phi_red = 0.2d0

  inversion%iteration = UNINITIALIZED_INTEGER
  inversion%current_chi2 = UNINITIALIZED_DOUBLE
  inversion%phi_total_0 = UNINITIALIZED_DOUBLE
  inversion%phi_data_0 = UNINITIALIZED_DOUBLE
  inversion%phi_model_0 = UNINITIALIZED_DOUBLE
  inversion%phi_total = UNINITIALIZED_DOUBLE
  inversion%phi_data = UNINITIALIZED_DOUBLE
  inversion%phi_model = UNINITIALIZED_DOUBLE

  inversion%cull_flag = PETSC_FALSE
  inversion%cull_dev = UNINITIALIZED_DOUBLE

  inversion%converg_flag = PETSC_FALSE

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

subroutine InversionCheckConvergence(inversion,survey)
  !
  ! Check Inversion convergence
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/05/21
      
  implicit none
  
  type(inversion_type) :: inversion
  type(survey_type) :: survey

  inversion%converg_flag = PETSC_FALSE
  call InversionEvaluateCostFunctions(inversion,survey)
  if (inversion%current_chi2 <= inversion%target_chi2) &
                           inversion%converg_flag = PETSC_TRUE
      
end subroutine InversionCheckConvergence

! ************************************************************************** !

subroutine InversionEvaluateCostFunctions(inversion,survey)
  !
  ! Evaluates cost functions for inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/05/21
        
  implicit none
        
  type(inversion_type) :: inversion
  type(survey_type) :: survey
  
  PetscInt :: idata,ndata,ncull
  PetscReal :: err_mean,err_sdev
  PetscReal, allocatable :: data_vector(:)

  ndata = survey%num_measurement
  allocate(data_vector(ndata))
  if (.not.associated(survey%Wd_cull)) allocate(survey%Wd_cull(ndata))
  data_vector = 0.d0
  survey%Wd_cull = 1

  data_vector = survey%Wd * (survey%dobs - survey%dsim)
    
  ncull = 0
  if (inversion%cull_flag) then
    err_mean = sum(data_vector) / ndata
    err_sdev = sqrt( dot_product( (data_vector - err_mean), &
                                  (data_vector - err_mean) ) / ndata )
    do idata=1,ndata
      if (data_vector(idata) < (err_mean - err_sdev*inversion%cull_dev) .or. &
          data_vector(idata) > (err_mean + err_sdev*inversion%cull_dev)) then
        survey%Wd_cull(idata) = 0
        ncull = ncull + 1
      endif
    enddo    
  endif

  data_vector = survey%Wd_cull * data_vector
    
  inversion%phi_data = dot_product(data_vector,data_vector)
  inversion%current_chi2 = inversion%phi_data / (ndata - ncull)
    
  deallocate(data_vector)

  ! TODO: compute phi_model
  inversion%phi_total = inversion%phi_data + inversion%phi_model

  if (inversion%iteration == 1) then
    inversion%phi_data_0 = inversion%phi_data
    inversion%phi_model_0 = inversion%phi_model
    inversion%phi_total_0 = inversion%phi_total
  endif  

end subroutine InversionEvaluateCostFunctions
  
! ************************************************************************** !

subroutine InversionCheckBeta(inversion)
  !
  ! Check Beta if it needs cooling/reduction
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/05/21
        
  implicit none
        
  type(inversion_type) :: inversion

  if ( abs(inversion%phi_total_0 - inversion%phi_total) <= &
       inversion%min_phi_red ) then
    inversion%beta = inversion%beta * inversion%beta_red_factor
    inversion%phi_model = inversion%beta_red_factor * inversion%phi_model
  endif

  ! update the cost functions
  inversion%phi_data_0 = inversion%phi_data
  inversion%phi_model_0 = inversion%phi_model
  inversion%phi_total_0 = inversion%phi_total

end subroutine InversionCheckBeta

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
