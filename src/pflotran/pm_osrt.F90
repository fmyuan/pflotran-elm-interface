module PM_OSRT_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes

  use PM_RT_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_rt_type) :: pm_osrt_type
    Vec :: fixed_accum
    Vec :: rhs
    PetscLogDouble :: cumulative_transport_time
    PetscLogDouble :: cumulative_reaction_time
  contains
    procedure, public :: InitializeRun => PMOSRTInitializeRun
    procedure, public :: InitializeTimestep => PMOSRTInitializeTimestep
    procedure, public :: FinalizeTimestep => PMOSRTFinalizeTimestep
    procedure, public :: UpdateTimestep => PMOSRTUpdateTimestep
    procedure, public :: PreSolve => PMOSRTPreSolve
    procedure, public :: PostSolve => PMOSRTPostSolve
    procedure, public :: AcceptSolution => PMOSRTAcceptSolution
    procedure, public :: FinalizeRun => PMOSRTFinalizeRun
    procedure, public :: Destroy => PMOSRTDestroy
  end type pm_osrt_type
  
  public :: PMOSRTCreate

contains

! ************************************************************************** !

function PMOSRTCreate()
  ! 
  ! Creates operator split reactive transport process model shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 

  implicit none
  
  class(pm_osrt_type), pointer :: PMOSRTCreate

  class(pm_osrt_type), pointer :: pm_osrt
  
  allocate(pm_osrt)
  call PMOSRTInit(pm_osrt)
  pm_osrt%name = 'Oper.-Split Reactive Transport'
  pm_osrt%header = 'OPER.-SPLIT REACTIVE TRANSPORT'
  pm_osrt%cumulative_transport_time = 0.d0
  pm_osrt%cumulative_reaction_time = 0.d0

  PMOSRTCreate => pm_osrt
  
end function PMOSRTCreate

! ************************************************************************** !

subroutine PMOSRTInit(pm_osrt)
  ! 
  ! Initializes reactive transport process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 
  implicit none

  class(pm_osrt_type) :: pm_osrt

  call PMRTInit(pm_osrt)

  pm_osrt%fixed_accum = PETSC_NULL_VEC
  pm_osrt%rhs = PETSC_NULL_VEC
  
  pm_osrt%operator_split = PETSC_TRUE

end subroutine PMOSRTInit

! ************************************************************************** !

recursive subroutine PMOSRTInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/19
  ! 
  use Discretization_module

  implicit none

  class(pm_osrt_type) :: this

  call PMRTInitializeRun(this)
  call DiscretizationDuplicateVector(this%realization%discretization, &
                                     this%realization%field%tran_xx, &
                                     this%fixed_accum)
  call DiscretizationDuplicateVector(this%realization%discretization, &
                                     this%realization%field%tran_xx, &
                                     this%rhs)

  ! check to ensure zero tran_dt performed at end of PMRTInitializeRun

end subroutine PMOSRTInitializeRun

! ************************************************************************** !

subroutine PMOSRTInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 

  use Realization_Subsurface_class
  use Reactive_Transport_module
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_TIMESTEP
  use Global_module
  use PM_Base_class
  use Material_module
  use Option_module

  implicit none
  
  class(pm_osrt_type) :: this

  class(realization_subsurface_type), pointer :: realization
  PetscErrorCode :: ierr

  realization => this%realization
 
  this%option%tran_dt = this%option%dt

  call PMRTWeightFlowParameters(this,TIME_T)
  
  call VecCopy(realization%field%tran_xx,realization%field%tran_yy, &
               ierr);CHKERRQ(ierr)
  ! should I be updating BCs at this point?
                                 ! cells      bcs        act coefs
  call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_FALSE)

  if (this%realization%reaction%act_coef_update_frequency == &
      ACT_COEF_FREQUENCY_TIMESTEP) then
    call RTUpdateActivityCoefficients(this%realization,PETSC_TRUE,PETSC_TRUE)
  endif

end subroutine PMOSRTInitializeTimestep

! ************************************************************************** !

subroutine PMOSRTPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 

  use Reactive_Transport_module, only : RTUpdateTransportCoefs
  use Global_module  
  use Material_module
  use Data_Mediator_module

  implicit none
  
  class(pm_osrt_type) :: this
  
  PetscErrorCode :: ierr
  
  ! set densities and saturations to t+dt
  if (this%option%nflowdof > 0 .and. .not. this%steady_flow) then
    if (this%option%flow%transient_porosity) then
      ! weight material properties (e.g. porosity)
      call MaterialWeightAuxVars(this%realization%patch%aux%Material, &
                                 this%tran_weight_t1, &
                                 this%realization%field,this%comm1)
    endif
    call GlobalWeightAuxVars(this%realization,this%tran_weight_t1)
  else if (this%transient_porosity) then
    this%tran_weight_t1 = 1.d0
    call MaterialWeightAuxVars(this%realization%patch%aux%Material, &
                               this%tran_weight_t1, &
                               this%realization%field,this%comm1)
  endif

  call RTUpdateTransportCoefs(this%realization)
  
#if 0
  ! the problem here is that activity coefficients will be updated every time
  ! presolve is called, regardless of TS vs NI.  We need to split this out.
  if (this%realization%reaction%act_coef_update_frequency /= &
      ACT_COEF_FREQUENCY_OFF) then
    call RTUpdateAuxVars(this%realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)
!       The below is set within RTUpdateAuxVarsPatch() when 
!         PETSC_TRUE,PETSC_TRUE,* are passed
!       patch%aux%RT%auxvars_up_to_date = PETSC_TRUE 
  endif
#endif

  if (this%realization%reaction%use_log_formulation) then
    call VecCopy(this%realization%field%tran_xx, &
                 this%realization%field%tran_log_xx,ierr);CHKERRQ(ierr)
    call VecLog(this%realization%field%tran_log_xx,ierr);CHKERRQ(ierr)
  endif
  
  call DataMediatorUpdate(this%realization%tran_data_mediator_list, &
                          this%realization%field%tran_mass_transfer, &
                          this%realization%option)
  
end subroutine PMOSRTPreSolve

! ************************************************************************** !

subroutine PMOSRTPostSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 

  implicit none
  
  class(pm_osrt_type) :: this
  
end subroutine PMOSRTPostSolve

! ************************************************************************** !

subroutine PMOSRTFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 

  use Reactive_Transport_module, only : RTMaxChange
  use Variables_module, only : POROSITY
  use Material_module, only : MaterialGetAuxVarVecLoc
  use Material_Aux_class, only : POROSITY_BASE 
  use Global_module

  implicit none
  
  class(pm_osrt_type) :: this
  PetscReal :: time  
  PetscErrorCode :: ierr

#if 0
  if (this%transient_porosity) then
    call VecCopy(this%realization%field%porosity_tpdt, &
                 this%realization%field%porosity_t,ierr);CHKERRQ(ierr)
    call RealizationUpdatePropertiesTS(this%realization)
    call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 POROSITY,POROSITY_BASE)
    call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                  this%realization%field%porosity_tpdt)
  endif
  
  call RTMaxChange(this%realization,this%max_concentration_change, &
                   this%max_volfrac_change)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dcmx= ",1pe12.4,"  dc/dt= ",1pe12.4, &
            &" [mol/s]")') &
      maxval(this%max_concentration_change), &
      maxval(this%max_concentration_change)/this%option%tran_dt
    if (this%realization%reaction%mineral%nkinmnrl > 0) then
      write(*,'("               dvfmx= ",1pe12.4," dvf/dt= ",1pe12.4, &
            &" [1/s]")') &
        maxval(this%max_volfrac_change), &
        maxval(this%max_volfrac_change)/this%option%tran_dt
    endif
  endif
  if (this%option%print_file_flag) then  
    write(this%option%fid_out,&
            '("  --> max chng: dcmx= ",1pe12.4,"  dc/dt= ",1pe12.4, &
            &" [mol/s]")') &
      maxval(this%max_concentration_change), &
      maxval(this%max_concentration_change)/this%option%tran_dt
    if (this%realization%reaction%mineral%nkinmnrl > 0) then
      write(this%option%fid_out, &
        '("               dvfmx= ",1pe12.4," dvf/dt= ",1pe12.4," [1/s]")') &
        maxval(this%max_volfrac_change), &
        maxval(this%max_volfrac_change)/this%option%tran_dt
    endif
  endif
#endif
  
end subroutine PMOSRTFinalizeTimestep

! ************************************************************************** !

function PMOSRTAcceptSolution(this)
  ! 
  ! PMRichardsAcceptSolution:
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 

  implicit none
  
  class(pm_osrt_type) :: this
  
  PetscBool :: PMOSRTAcceptSolution
  
  PMOSRTAcceptSolution = PETSC_TRUE
  
end function PMOSRTAcceptSolution

! ************************************************************************** !

subroutine PMOSRTUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                num_newton_iterations,tfac, &
                              time_step_max_growth_factor)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 
  use Realization_Subsurface_class

  implicit none
  
  class(pm_osrt_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  PetscReal, parameter :: pert = 1.d-20
  PetscReal :: dtt
  PetscReal :: uvf

  dtt = 1.d20
  if (this%volfrac_change_governor < 1.d0) then
    uvf= this%volfrac_change_governor/(maxval(this%max_volfrac_change)+pert)
    dtt = 0.5d0 * dt * (1.d0 + uvf)
  endif

  dtt = min(time_step_max_growth_factor*dt,dtt)
  if (dtt > dt_max) dtt = dt_max
  ! geh: see comment above under flow stepper
  dtt = max(dtt,dt_min)
  dt = dtt
  
  call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)
  
end subroutine PMOSRTUpdateTimestep

! ************************************************************************** !

recursive subroutine PMOSRTFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 
  use Option_module

  implicit none
  
  class(pm_osrt_type) :: this


  if (OptionPrintToScreen(this%option)) then
    write(*,'(/," Transport Time: ", es12.4, " [sec]",/,&
               &"  Reaction Time: ", es12.4, " [sec]")') &
            this%cumulative_transport_time, & 
            this%cumulative_reaction_time
  endif
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMOSRTFinalizeRun

! ************************************************************************** !

subroutine PMOSRTStrip(this)
  ! 
  ! Strips members of OSRT process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 
  implicit none

  class(pm_osrt_type) :: this

  PetscErrorCode :: ierr

  call PMRTStrip(this)

  if (this%fixed_accum /= PETSC_NULL_VEC) then
    call VecDestroy(this%fixed_accum,ierr);CHKERRQ(ierr)
  endif
  if (this%rhs /= PETSC_NULL_VEC) then
    call VecDestroy(this%rhs,ierr);CHKERRQ(ierr)
  endif

end subroutine PMOSRTStrip

! ************************************************************************** !

subroutine PMOSRTDestroy(this)
  ! 
  ! Destroys OSRT process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/09/19
  ! 
  class(pm_osrt_type) :: this

  call PMOSRTStrip(this)

end subroutine PMOSRTDestroy
  
end module PM_OSRT_class
