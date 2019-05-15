module PM_TH_TS_class

#include "petsc/finclude/petscts.h"
#include "petsc/finclude/petscvec.h"
  use petscts
  use petscvec
  use TH_module
  use TH_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_TH_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_th_type) :: pm_th_ts_type
  contains
    procedure, public :: UpdateAuxVars => PMTHTSUpdateAuxVars
    procedure, public :: IFunction => PMTHTSIFunction
    procedure, public :: IJacobian => PMTHTSIJacobian
    procedure, public :: InitializeTimestep => PMTHTSInitializeTimestep
    procedure, public :: Destroy => PMTHTSDestroy
    procedure, public :: CheckConvergence => PMTHTSCheckConvergence
  end type pm_th_ts_type

  public :: PMTHTSCreate, &
            PMTHTSUpdateAuxVarsPatch
  	
  
contains

! ************************************************************************** !

function PMTHTSCreate()
  ! 
  ! Creates TH TS process models shell
  ! 
  ! Author: Satish Karra
  ! Date: 05/08/19
  ! 

  implicit none
  
  class(pm_th_ts_type), pointer :: PMTHTSCreate

  class(pm_th_ts_type), pointer :: this

  PetscReal, parameter :: pres_abs_inf_tol = 1.d0
  PetscReal, parameter :: temp_abs_inf_tol = 1.d-5
  PetscReal, parameter :: abs_update_inf_tol(2) = &
                            [pres_abs_inf_tol,temp_abs_inf_tol]
  PetscReal, parameter :: pres_rel_inf_tol = 1.d-5
  PetscReal, parameter :: temp_rel_inf_tol = 1.d-5
  PetscReal, parameter :: rel_update_inf_tol(2) = &
                            [pres_rel_inf_tol,temp_rel_inf_tol]
  PetscReal, parameter :: residual_abs_inf_tol(2) = 1.d-5
  PetscReal, parameter :: residual_scaled_inf_tol(2) = 1.d-5
  
#ifdef PM_TH_TS_DEBUG
  print *, 'PMTHTSCreate()'
#endif  

  allocate(this)

  nullify(this%commN)

  call PMSubsurfaceFlowCreate(this)
  this%name = 'TH_TS Flow'
  this%header = 'TH_TS FLOW'

  this%residual_abs_inf_tol = residual_abs_inf_tol
  this%residual_scaled_inf_tol = residual_scaled_inf_tol
  this%abs_update_inf_tol = abs_update_inf_tol
  this%rel_update_inf_tol = rel_update_inf_tol

  PMTHTSCreate => this
  
end function PMTHTSCreate

! ************************************************************************** !

subroutine PMTHTSUpdateAuxVars(this)
  ! 
  ! Author: Satish Karra
  ! Date: 05/08/19

  implicit none
  
  class(pm_th_ts_type) :: this

  call PMTHTSUpdateAuxVarsPatch(this%realization)

end subroutine PMTHTSUpdateAuxVars


! ************************************************************************** !

subroutine PMTHTSUpdateAuxVarsPatch(realization)
  ! 
  ! Author: Satish Karra
  ! Date: 05/08/19


  use Realization_Subsurface_class
  use Field_module
  use Grid_module
  use Patch_module
  use Option_module
  
  implicit none

  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(TH_auxvar_type), pointer :: TH_auxvars(:) 
  PetscInt :: ghosted_id,local_id,istart,iend
  PetscReal, pointer :: xx_loc_p(:),xxdot_loc_p(:)
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  TH_auxvars => patch%aux%TH%auxvars

  ! 1. Update auxvars based on new values of pressure, temperature
  call THUpdateAuxVars(realization)

  ! 2. Update auxvars based on new value of dpressure_dtime, mass, and 
  !    dmass_dtime
  call VecGetArrayReadF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_xxdot_loc,xxdot_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle 
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    th_auxvars(ghosted_id)%dpres_dtime = xxdot_loc_p((ghosted_id-1)*option%nflowdof+1)
    th_auxvars(ghosted_id)%dtemp_dtime = xxdot_loc_p((ghosted_id-1)*option%nflowdof+2)
!    call RichardsAuxVarCompute2ndOrderDeriv(rich_auxvars(ghosted_id), &
!                                       patch%characteristic_curves_array( &
!                                         patch%sat_func_id(ghosted_id))%ptr, &
!                                       option)
  enddo

  call VecRestoreArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_xxdot_loc,xxdot_loc_p, &
                              ierr);CHKERRQ(ierr)
  

end subroutine PMTHTSUpdateAuxVarsPatch

! ************************************************************************** !

subroutine PMTHTSIFunction(this,ts,time,U,Udot,F,ierr)
  !
  !
  ! Author: Satish Karra
  ! Date: 05/08/19
  
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module

  implicit none
  
  class(pm_th_ts_type) :: this

  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  Vec :: F
  PetscErrorCode :: ierr

end subroutine PMTHTSIFunction


! ************************************************************************** !
subroutine IFunctionAccumulation(F,realization,ierr)
  !
  !
  !
  ! Author: Satish Karra
  ! Date: 06/07/18
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Patch_module
  use Option_module
  use Material_Aux_class

  implicit none
  
  class(realization_subsurface_type), pointer :: realization
  Vec :: F
  PetscErrorCode :: ierr
  
end subroutine IFunctionAccumulation

! ************************************************************************** !

subroutine PMTHTSIJacobian(this,ts,time,U,Udot,shift,A,B,ierr)
  !
  !
  !
  ! Author: Satish Karra
  ! Date: 06/07/18
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module

  implicit none
  
  class(pm_th_ts_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  PetscReal :: shift
  Mat :: A, B
  PetscErrorCode :: ierr  
  
end subroutine PMTHTSIJacobian


! ************************************************************************** !
subroutine IJacobianAccumulation(J,shift,realization,ierr)
  !
  !
  !
  ! Author: Satish Karra
  ! Date: 06/07/18
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Patch_module
  use Option_module
  use Material_Aux_class

  implicit none
  
  class(realization_subsurface_type), pointer :: realization
  PetscReal :: shift
  Mat :: J
  PetscErrorCode :: ierr
  
    
end subroutine IJacobianAccumulation

! ************************************************************************** !

subroutine PMTHTSInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/08/19
  !

  use TH_module, only : THInitializeTimestep
  
  implicit none
  
  class(pm_th_ts_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMTHTSInitializeTimestep

! ************************************************************************** !

subroutine PMTHTSCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)
  ! 
  ! Adds a convergence check for the nonlinear problem
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/08/19
  ! 

#include "petsc/finclude/petscsnes.h"
  use petscsnes

  implicit none
  
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm ! 2-norm of updated solution
  PetscReal :: unorm ! 2-norm of update. PETSc refers to this as snorm
  PetscReal :: fnorm ! 2-norm of updated residual
  SNESConvergedReason :: reason
  class(pm_th_ts_type) :: this
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string

  call SNESConvergedDefault(snes,it,xnorm,unorm,fnorm,reason, &
                            0,ierr);CHKERRQ(ierr)

  if (this%option%print_screen_flag) then
    select case(int(reason))
      case(2)
        string = 'atol'
      case(3)
        string = 'rtol'
      case(4)
        string = 'stol'
      case default
        write(string,'(i3)') reason
    end select

    write(*,'(i3," 2r:",es9.2, &
            & " 2x:",es9.2, &
            & " 2u:",es9.2, &
            & " rsn: ",a)') &
            it, fnorm, xnorm, unorm, trim(string)
  endif

end subroutine PMTHTSCheckConvergence

! ************************************************************************** !

subroutine PMTHTSDestroy(this)
  ! 
  ! Destroys TH process model
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/08/19
  ! 
  use TH_module, only : THDestroy

  implicit none
  
  class(pm_th_ts_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call PMTHDestroy(this)
  
end subroutine PMTHTSDestroy

end module PM_TH_TS_class

