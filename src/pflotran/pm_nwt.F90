module PM_NWT_class
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class 
  use Realization_Subsurface_class
  use Communicator_Base_module  
  use Option_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
    
  type, public :: pm_nwt_controls_type
    PetscReal, pointer :: max_concentration_change(:)
    PetscReal, pointer :: max_volfrac_change(:)
    PetscReal :: volfrac_change_governor
    PetscReal :: cfl_governor
    PetscReal :: newton_inf_rel_update_tol
    PetscReal :: newton_inf_scaled_res_tol
    PetscBool :: check_post_converged
    PetscBool :: check_post_convergence
    PetscBool :: check_update
#ifdef OS_STATISTICS
! use PetscReal for large counts
    PetscInt :: newton_call_count
    PetscReal :: sum_newton_call_count
    PetscInt :: newton_iterations
    PetscReal :: sum_newton_iterations
    PetscInt :: max_newton_iterations
    PetscInt :: overall_max_newton_iterations
#endif    
  end type pm_nwt_controls_type
  
  type, public :: pm_nwt_params_type
    PetscInt :: nphase
    PetscInt :: ncomp
    PetscInt :: nsorb
    PetscInt :: nmnrl
    PetscInt :: nauxiliary
    PetscBool :: calculate_transverse_dispersion
    PetscBool :: temperature_dependent_diffusion
  end type pm_nwt_params_type
  
  type, public, extends(pm_base_type) :: pm_nwt_type
  ! realization_base_type has the nwt object (equivalent to reaction)
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    type(pm_nwt_controls_type), pointer :: controls
    type(pm_nwt_params_type), pointer :: params
  contains
    procedure, public :: Setup => PMNWTSetup 
    procedure, public :: ReadSimulationBlock => PMNWTReadSimulationBlock
    procedure, public :: SetRealization => PMNWTSetRealization 
    procedure, public :: InitializeRun => PMNWTInitializeRun  
  end type pm_nwt_type
  
  public :: PMNWTCreate, PMNWTSetPlotVariables
  
  
contains

! ************************************************************************** !

function PMNWTCreate()
  ! 
  ! Creates the nuclear waste transport process model shell.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  ! 

  implicit none
  
  class(pm_nwt_type), pointer :: PMNWTCreate

  class(pm_nwt_type), pointer :: nwt_pm
    
  allocate(nwt_pm)
  nullify(nwt_pm%option)
  nullify(nwt_pm%output_option)
  nullify(nwt_pm%realization)
  nullify(nwt_pm%comm1)
    
  allocate(nwt_pm%controls)
  nullify(nwt_pm%controls%max_concentration_change)
  nullify(nwt_pm%controls%max_volfrac_change)
  nwt_pm%controls%volfrac_change_governor = 1.d0
  nwt_pm%controls%cfl_governor = UNINITIALIZED_DOUBLE
  nwt_pm%controls%newton_inf_rel_update_tol = UNINITIALIZED_DOUBLE
  nwt_pm%controls%newton_inf_scaled_res_tol = UNINITIALIZED_DOUBLE
  nwt_pm%controls%check_post_converged = PETSC_FALSE
  nwt_pm%controls%check_post_convergence = PETSC_FALSE
  nwt_pm%controls%check_update = PETSC_FALSE
#ifdef OS_STATISTICS
  nwt_pm%controls%newton_call_count = 0
  nwt_pm%controls%sum_newton_call_count = 0.d0
  nwt_pm%controls%newton_iterations = 0
  nwt_pm%controls%sum_newton_iterations = 0.d0
  nwt_pm%controls%max_newton_iterations = 0
  nwt_pm%controls%overall_max_newton_iterations = 0
#endif
  
  allocate(nwt_pm%params)
  nwt_pm%params%ncomp = 0
  nwt_pm%params%nphase = 1  ! For WIPP, we always assume liquid phase only
  nwt_pm%params%nsorb = 0
  nwt_pm%params%nmnrl = 0
  nwt_pm%params%nauxiliary = 0
  nwt_pm%params%calculate_transverse_dispersion = PETSC_FALSE
  nwt_pm%params%temperature_dependent_diffusion = PETSC_FALSE
  

  call PMBaseInit(nwt_pm)
  nwt_pm%name = 'Nuclear Waste Transport'
  nwt_pm%header = 'NUCLEAR WASTE TRANSPORT'
  
  PMNWTCreate => nwt_pm
  
end function PMNWTCreate

! ************************************************************************** !

subroutine PMNWTReadSimulationBlock(this,input)
  ! 
  ! Reads input file parameters associated with the nuclear waste transport 
  ! process model in the SIMULATION block.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use NW_Transport_Aux_module
 
  implicit none
  
  class(pm_nwt_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option
  
  error_string = 'NUCLEAR_WASTE_TRANSPORT OPTIONS'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    found = PETSC_FALSE
    call PMBaseReadSelectCase(this,input,keyword,found,error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('GLOBAL_IMPLICIT')
        option%transport%nw_transport_coupling = GLOBAL_IMPLICIT
      case('OPERATOR_SPLIT','OPERATOR_SPLITTING')
      case('MAX_VOLUME_FRACTION_CHANGE')
        call InputReadDouble(input,option,this%controls%volfrac_change_governor)
        call InputErrorMsg(input,option,'MAX_VOLUME_FRACTION_CHANGE', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
      case('ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,nwt_itol_rel_update)
        call InputErrorMsg(input,option,'ITOL_RELATIVE_UPDATE', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
        this%controls%check_post_convergence = PETSC_TRUE
      case('MINIMUM_SATURATION')
        call InputReadDouble(input,option,nwt_min_saturation)
        call InputErrorMsg(input,option,'MINIMUM_SATURATION', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
      case('NUMERICAL_JACOBIAN')
        option%transport%numerical_derivatives = PETSC_TRUE
      ! jenn:todo Why is temperature_dependent_diffusion in the SIMULATION block read?
      case('TEMPERATURE_DEPENDENT_DIFFUSION')
        this%params%temperature_dependent_diffusion = PETSC_TRUE
      case('MAX_CFL')
        call InputReadDouble(input,option,this%controls%cfl_governor)
        call InputErrorMsg(input,option,'MAX_CFL', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
      case('MULTIPLE_CONTINUUM')
        option%use_mc = PETSC_TRUE          
      case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  enddo
  
end subroutine PMNWTReadSimulationBlock

! ************************************************************************** !
  
subroutine PMNWTSetup(this)
  ! 
  ! Initializes variables associated with nuclear waste transport.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  ! 
  
  use NW_Transport_Aux_module

  implicit none
  
  class(pm_nwt_type) :: this
    
  type(nw_trans_realization_type), pointer :: nw_trans
  
  nw_trans => this%realization%nw_trans
  
  this%params%nphase = nw_trans%params%nphase
  this%params%ncomp = nw_trans%params%ncomp
  this%params%nsorb = nw_trans%params%nsorb
  this%params%nmnrl = nw_trans%params%nmnrl
  this%params%nauxiliary = nw_trans%params%nauxiliary
  this%params%calculate_transverse_dispersion = &
                               nw_trans%params%calculate_transverse_dispersion
  this%params%temperature_dependent_diffusion = &
                               nw_trans%params%temperature_dependent_diffusion
        
  ! set the communicator
  this%comm1 => this%realization%comm1
  
  allocate(this%controls%max_concentration_change(this%params%ncomp))
  allocate(this%controls%max_volfrac_change(this%params%ncomp))

end subroutine PMNWTSetup

! ************************************************************************** !

subroutine PMNWTSetRealization(this,realization)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  ! 

  use Realization_Subsurface_class  

  implicit none
  
  class(pm_nwt_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization
  
  if (this%realization%nw_trans%use_log_formulation) then
    this%solution_vec = realization%field%tran_log_xx
  else
    this%solution_vec = realization%field%tran_xx
  endif
  this%residual_vec = realization%field%tran_r
  
end subroutine PMNWTSetRealization

! ************************************************************************** !

subroutine PMNWTInitializeRun(this)
  ! 
  ! Initializes process model time stepping.
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/02/2019
  ! 
  
  use NW_Transport_module
  
  implicit none
  
  class(pm_nwt_type) :: this
  
  ! check for uninitialized flow variables
  call RealizUnInitializedVarsTran(this%realization)
  
  ! update boundary conditions (don't know if this is necessary)
  ! jenn:todo Is updating BCs necessary in PMNWTInitializeRun()?
  call NWTUpdateAuxVars(this%realization,PETSC_FALSE,PETSC_TRUE)
  
end subroutine PMNWTInitializeRun

! ************************************************************************** !

subroutine PMNWTSetPlotVariables(list,nw_trans,option,time_unit)
  ! 
  ! Adds variables to be printed for plotting.
  !
  ! Author: Jenn Frederick
  ! Date: 03/28/2019
  !
  
  use Option_module
  use Output_Aux_module
  use Variables_module
  use NW_Transport_Aux_module
    
  implicit none
  
  type(output_variable_list_type), pointer :: list
  type(nw_trans_realization_type), pointer :: nw_trans
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: time_unit
  
  character(len=MAXWORDLENGTH) :: name,  units
  PetscInt :: i
  
  ! jenn:todo Where should I be setting nw_trans%print ? In an equivalent
  ! block like CHEMISTRY,OUTPUT?
  if (nw_trans%print%molality) then
    do i=1,nw_trans%params%ncomp
      if (nw_trans%species_print(i)) then
        name = 'Molality ' // trim(nw_trans%species_names(i))
        units = 'm'
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     PRIMARY_MOLALITY,i)
      endif
    enddo
  endif 
  
end subroutine PMNWTSetPlotVariables

! ************************************************************************** !
  
! jenn:todo Remember to deallocate/destroy all PMNWT pointers and arrays

end module PM_NWT_class
