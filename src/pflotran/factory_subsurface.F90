module Factory_Subsurface_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Simulation_Subsurface_class

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none

  private

  public :: SubsurfaceInitialize, &
            SubsurfaceInitializePostPETSc, &
            SubsurfaceJumpStart, &
            ! move to init_subsurface
            SubsurfaceReadFlowPM, &
            SubsurfaceReadRTPM, &
            SubsurfaceReadNWTPM, &
            SubsurfaceReadWasteFormPM, &
            SubsurfaceReadUFDDecayPM, &
            SubsurfaceReadUFDBiospherePM

contains

! ************************************************************************** !

subroutine SubsurfaceInitialize(simulation)
  !
  ! Sets up PFLOTRAN subsurface simulation
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  !

  use WIPP_module
  use Klinkenberg_module

  implicit none

  class(simulation_subsurface_type) :: simulation

  ! Modules that must be initialized
  call WIPPInit()
  call KlinkenbergInit()

  ! NOTE: PETSc must already have been initialized here!
  call SubsurfaceInitializePostPetsc(simulation)

end subroutine SubsurfaceInitialize

! ************************************************************************** !

subroutine SubsurfaceInitializePostPetsc(simulation)
  !
  ! Sets up PFLOTRAN subsurface simulation
  ! framework after to PETSc initialization
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/13
  !

  use Option_module
  use PM_Subsurface_Flow_class
  use PM_Base_class
  use PM_RT_class
  use PM_NWT_class
  use PM_Waste_Form_class
  use PM_UFD_Decay_class
  use PM_UFD_Biosphere_class
  use PM_Auxiliary_class
  use Timestepper_BE_class
  use Realization_Subsurface_class
  use Logging_module
  use Simulation_Subsurface_class
  use Solver_module
  use Waypoint_module
  use Init_Common_module
  use Init_Subsurface_module
  use Input_Aux_module
  use String_module
  use Checkpoint_module

  implicit none

  class(simulation_subsurface_type) :: simulation

  type(option_type), pointer :: option
  class(pm_subsurface_flow_type), pointer :: pm_flow
  class(pm_rt_type), pointer :: pm_rt
  class(pm_nwt_type), pointer :: pm_nwt
  class(pm_waste_form_type), pointer :: pm_waste_form
  class(pm_ufd_decay_type), pointer :: pm_ufd_decay
  class(pm_ufd_biosphere_type), pointer :: pm_ufd_biosphere
  class(pm_auxiliary_type), pointer :: pm_auxiliary
  class(realization_subsurface_type), pointer :: realization

  option => simulation%option

  ! process command line arguments specific to subsurface
  call SubsurfInitCommandLineSettings(option)

  call ExtractPMsFromPMList(simulation,pm_flow,pm_rt,pm_nwt,pm_waste_form,&
                            pm_ufd_decay,pm_ufd_biosphere,pm_auxiliary)

  call SubsurfaceSetFlowMode(pm_flow,option)

  realization => RealizationCreate(option)
  simulation%realization => realization
  realization%output_option => simulation%output_option

  simulation%waypoint_list_subsurface => WaypointListCreate()

  ! Setup linkages between PMCs
  call SetupPMCLinkages(simulation,pm_flow,pm_rt,pm_nwt,pm_waste_form,&
    pm_ufd_decay,pm_ufd_biosphere,pm_auxiliary,realization)
  
  ! SubsurfaceInitSimulation() must be called after pmc linkages are set above.
  call SubsurfaceInitSimulation(simulation)

  ! set first process model coupler as the master
  simulation%process_model_coupler_list%is_master = PETSC_TRUE

end subroutine SubsurfaceInitializePostPetsc

! ************************************************************************** !

subroutine ExtractPMsFromPMList(simulation,pm_flow,pm_rt,pm_nwt,pm_waste_form,&
                                pm_ufd_decay,pm_ufd_biosphere,pm_auxiliary)
  !
  ! Extracts all possible PMs from the PM list
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PM_Subsurface_Flow_class
  use PM_Base_class
  use PM_RT_class
  use PM_NWT_class
  use PM_Waste_Form_class
  use PM_UFD_Decay_class
  use PM_UFD_Biosphere_class
  use PM_Auxiliary_class
  use Option_module
  use Simulation_Subsurface_class

  implicit none

  class(simulation_subsurface_type) :: simulation

  type(option_type), pointer :: option
  class(pm_subsurface_flow_type), pointer :: pm_flow
  class(pm_rt_type), pointer :: pm_rt
  class(pm_nwt_type), pointer :: pm_nwt
  class(pm_waste_form_type), pointer :: pm_waste_form
  class(pm_ufd_decay_type), pointer :: pm_ufd_decay
  class(pm_ufd_biosphere_type), pointer :: pm_ufd_biosphere
  class(pm_auxiliary_type), pointer :: pm_auxiliary
  class(pm_base_type), pointer :: cur_pm, prev_pm

  option => simulation%option

  nullify(pm_flow)
  nullify(pm_rt)
  nullify(pm_nwt)
  nullify(pm_waste_form)
  nullify(pm_ufd_decay)
  nullify(pm_ufd_biosphere)
  nullify(pm_auxiliary)

  cur_pm => simulation%process_model_list
  do
    if (.not.associated(cur_pm)) exit
    select type(cur_pm)
      class is(pm_subsurface_flow_type)
        pm_flow => cur_pm
      class is(pm_rt_type)
        pm_rt => cur_pm
      class is(pm_nwt_type)
        pm_nwt => cur_pm
      class is(pm_waste_form_type)
        pm_waste_form => cur_pm
      class is(pm_ufd_decay_type)
        pm_ufd_decay => cur_pm
      class is(pm_ufd_biosphere_type)
        pm_ufd_biosphere => cur_pm
      class is(pm_auxiliary_type)
        pm_auxiliary => cur_pm
      class default
        option%io_buffer = &
         'PM Class unrecognized in SubsurfaceInitializePostPetsc.'
        call PrintErrMsg(option)
    end select

    prev_pm => cur_pm
    cur_pm => cur_pm%next

    ! we must destroy the linkage between pms so that they are in independent
    ! lists among pmcs
    nullify(prev_pm%next)

  enddo

end subroutine ExtractPMsFromPMList

! ************************************************************************** !

subroutine SetupPMCLinkages(simulation,pm_flow,pm_rt,pm_nwt,pm_waste_form,&
                            pm_ufd_decay,pm_ufd_biosphere,pm_auxiliary, &
                            realization)
  !
  ! Sets up all PMC linkages
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PM_Subsurface_Flow_class
  use PM_Base_class
  use PM_RT_class
  use PM_NWT_class
  use PM_Waste_Form_class
  use PM_UFD_Decay_class
  use PM_UFD_Biosphere_class
  use PM_Auxiliary_class
  use Realization_Subsurface_class
  use Option_module
  use Input_Aux_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_subsurface_flow_type), pointer :: pm_flow
  class(pm_rt_type), pointer :: pm_rt
  class(pm_nwt_type), pointer :: pm_nwt
  class(pm_waste_form_type), pointer :: pm_waste_form
  class(pm_ufd_decay_type), pointer :: pm_ufd_decay
  class(pm_ufd_biosphere_type), pointer :: pm_ufd_biosphere
  class(pm_auxiliary_type), pointer :: pm_auxiliary
  class(realization_subsurface_type), pointer :: realization

  type(option_type), pointer :: option
  type(input_type), pointer :: input

  option => simulation%option

  if (associated(pm_flow)) &
    call AddPMCSubsurfaceFlow(simulation,pm_flow,'PMCSubsurfaceFlow', &
                              realization,option)

  if (associated(pm_rt))   &
    call AddPMCSubsurfaceRT(simulation,pm_rt,'PMCSubsurfaceTransport', &
                            realization,option)
                            
  if (associated(pm_nwt))   &
    call AddPMCSubsurfaceNWT(simulation,pm_nwt,'PMCSubsurfaceNWTransport', &
                            realization,option)
                            
  input => InputCreate(IN_UNIT,option%input_filename,option)
  call SubsurfaceReadRequiredCards(simulation,input)
  call SubsurfaceReadInput(simulation,input)

  if (associated(pm_waste_form)) &
    call AddPMCWasteForm(simulation,pm_waste_form,'PMC3PWasteForm',&
                         associated(pm_ufd_decay),realization,input,option)

  if (associated(pm_ufd_decay)) &
    call AddPMCUDFDecay(simulation,pm_ufd_decay,'PMC3PUFDDecay',realization, &
                        input,option)

  if (associated(pm_ufd_biosphere)) &
    call AddPMCUDFBiosphere(simulation,pm_ufd_biosphere,'PMC3PUFDBiosphere',&
                            associated(pm_ufd_decay),realization,input,option)

  if (associated(pm_auxiliary)) &
    call AddPMCAuxiliary(simulation,pm_auxiliary,'SALINITY',realization,option)

  call InputDestroy(input)

end subroutine SetupPMCLinkages

! ************************************************************************** !

subroutine AddPMCSubsurfaceFlow(simulation,pm_flow,pmc_name,realization,option)

  !
  ! Adds a subsurface flow PMC
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PM_Subsurface_Flow_class
  use PMC_Subsurface_class
  use Realization_Subsurface_class
  use String_module
  use Option_module
  use Logging_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_subsurface_flow_type), pointer :: pm_flow
  character(len=*) :: pmc_name
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option

  class(pmc_subsurface_type), pointer :: pmc_subsurface
  character(len=MAXSTRINGLENGTH) :: string

  pmc_subsurface => PMCSubsurfaceCreate()

  call pmc_subsurface%SetName(pmc_name)
  call pmc_subsurface%SetOption(option)
  call pmc_subsurface%SetCheckpointOption(simulation%checkpoint_option)
  call pmc_subsurface%SetWaypointList(simulation%waypoint_list_subsurface)

  pmc_subsurface%pm_list => pm_flow
  pmc_subsurface%pm_ptr%pm => pm_flow
  pmc_subsurface%realization => realization

  ! set up logging stage
  string = trim(pm_flow%name)
  call LoggingCreateStage(string,pmc_subsurface%stage)
  simulation%flow_process_model_coupler => pmc_subsurface
  simulation%process_model_coupler_list => simulation%flow_process_model_coupler

end subroutine AddPMCSubsurfaceFlow

! ************************************************************************** !

subroutine AddPMCSubsurfaceRT(simulation,pm_rt,pmc_name,realization,option)

  !
  ! Adds a subsurface reactive transport PMC
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PMC_Base_class
  use PM_RT_class
  use PMC_Subsurface_class
  use Realization_Subsurface_class
  use String_module
  use Option_module
  use Logging_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_rt_type), pointer :: pm_rt
  character(len=*) :: pmc_name
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option

  class(pmc_subsurface_type), pointer :: pmc_subsurface
  character(len=MAXSTRINGLENGTH) :: string
  class(pmc_base_type), pointer :: pmc_dummy

  nullify(pmc_dummy)

  pmc_subsurface => PMCSubsurfaceCreate()
  call pmc_subsurface%SetName(pmc_name)
  call pmc_subsurface%SetOption(option)
  call pmc_subsurface%SetCheckpointOption(simulation%checkpoint_option)
  call pmc_subsurface%SetWaypointList(simulation%waypoint_list_subsurface)
  pmc_subsurface%pm_list => pm_rt
  pmc_subsurface%pm_ptr%pm => pm_rt
  pmc_subsurface%realization => realization

  ! set up logging stage
  string = trim(pm_rt%name)
  call LoggingCreateStage(string,pmc_subsurface%stage)
  simulation%rt_process_model_coupler => pmc_subsurface

  if (.not.associated(simulation%process_model_coupler_list)) then
    simulation%process_model_coupler_list => pmc_subsurface
  else
    call PMCBaseSetChildPeerPtr(PMCCastToBase(pmc_subsurface),PM_CHILD, &
                      PMCCastToBase(simulation%flow_process_model_coupler), &
                      pmc_dummy,PM_INSERT)
  endif

end subroutine AddPMCSubsurfaceRT

! ************************************************************************** !

subroutine AddPMCSubsurfaceNWT(simulation,pm_nwt,pmc_name,realization,option)

  !
  ! Adds a subsurface nuclear waste reactive transport PMC.
  !
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  !

  use PMC_Base_class
  use PM_NWT_class
  use PMC_Subsurface_class
  use Realization_Subsurface_class
  use String_module
  use Option_module
  use Logging_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_nwt_type), pointer :: pm_nwt
  character(len=*) :: pmc_name
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option

  class(pmc_subsurface_type), pointer :: pmc_subsurface
  character(len=MAXSTRINGLENGTH) :: string
  class(pmc_base_type), pointer :: pmc_dummy

  nullify(pmc_dummy)

  pmc_subsurface => PMCSubsurfaceCreate()
  call pmc_subsurface%SetName(pmc_name)
  call pmc_subsurface%SetOption(option)
  call pmc_subsurface%SetCheckpointOption(simulation%checkpoint_option)
  call pmc_subsurface%SetWaypointList(simulation%waypoint_list_subsurface)
  pmc_subsurface%pm_list => pm_nwt
  pmc_subsurface%pm_ptr%pm => pm_nwt
  pmc_subsurface%realization => realization

  ! set up logging stage
  string = trim(pm_nwt%name)
  call LoggingCreateStage(string,pmc_subsurface%stage)
  simulation%nwt_process_model_coupler => pmc_subsurface

  if (.not.associated(simulation%process_model_coupler_list)) then
    simulation%process_model_coupler_list => pmc_subsurface
  else
    call PMCBaseSetChildPeerPtr(PMCCastToBase(pmc_subsurface),PM_CHILD, &
                      PMCCastToBase(simulation%flow_process_model_coupler), &
                      pmc_dummy,PM_INSERT)
  endif

end subroutine AddPMCSubsurfaceNWT

! ************************************************************************** !

subroutine AddPMCWasteForm(simulation,pm_waste_form,pmc_name,&
                           pm_ufd_decay_present,realization,input,option)

  !
  ! Adds a waste form PMC
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PMC_Base_class
  use PMC_Third_Party_class
  use PM_Waste_Form_class
  use Realization_Subsurface_class
  use String_module
  use Option_module
  use Logging_module
  use Input_Aux_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_waste_form_type), pointer :: pm_waste_form
  character(len=*) :: pmc_name
  logical :: pm_ufd_decay_present
  class(realization_subsurface_type), pointer :: realization
  type(input_type), pointer :: input
  type(option_type), pointer :: option

  class(pmc_third_party_type), pointer :: pmc_waste_form
  class(wf_mechanism_base_type), pointer :: cur_mechanism
  character(len=MAXSTRINGLENGTH) :: string
  class(pmc_base_type), pointer :: pmc_dummy

  nullify(pmc_dummy)

  string = 'WASTE_FORM_GENERAL'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  call pm_waste_form%ReadPMBlock(input)

  if (.not.associated(simulation%rt_process_model_coupler)) then
     option%io_buffer = 'The Waste Form process model requires &
          &reactive transport.'
     call PrintErrMsg(option)
  endif
  cur_mechanism => pm_waste_form%mechanism_list
  do
     if (.not.associated(cur_mechanism)) exit
     select type(cur_mechanism)
     type is(wf_mechanism_wipp_type)
        if (.not.pm_ufd_decay_present) then
           option%io_buffer = 'The WIPP type waste form mechanism requires &
                &the UFD_DECAY process model.'
           call PrintErrMsg(option)
        endif
     end select
     cur_mechanism => cur_mechanism%next
  enddo

  pmc_waste_form => PMCThirdPartyCreate()
  call pmc_waste_form%SetName(pmc_name)
  call pmc_waste_form%SetOption(option)
  call pmc_waste_form%SetCheckpointOption(simulation%checkpoint_option)
  pmc_waste_form%pm_list => pm_waste_form
  pmc_waste_form%pm_ptr%pm => pm_waste_form
  pmc_waste_form%realization => realization

  ! set up logging stage
  string = 'WASTE_FORM_GENERAL'
  call LoggingCreateStage(string,pmc_waste_form%stage)
  call PMCBaseSetChildPeerPtr(PMCCastToBase(pmc_waste_form),PM_CHILD, &
         PMCCastToBase(simulation%rt_process_model_coupler), &
         pmc_dummy,PM_APPEND)

end subroutine AddPMCWasteForm

! ************************************************************************** !

subroutine AddPMCUDFDecay(simulation,pm_ufd_decay,pmc_name,&
                          realization,input,option)

  !
  ! Adds a UFD decay PMC
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PMC_Base_class
  use PMC_Third_Party_class
  use PM_UFD_Decay_class
  use Realization_Subsurface_class
  use String_module
  use Option_module
  use Logging_module
  use Input_Aux_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_ufd_decay_type), pointer :: pm_ufd_decay
  character(len=*) :: pmc_name
  class(realization_subsurface_type), pointer :: realization
  type(input_type), pointer :: input
  type(option_type), pointer :: option

  class(pmc_third_party_type), pointer :: pmc_ufd_decay
  character(len=MAXSTRINGLENGTH) :: string
  class(pmc_base_type), pointer :: pmc_dummy

  nullify(pmc_dummy)

  string = 'UFD_DECAY'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  call pm_ufd_decay%ReadPMBlock(input)

  if (.not.associated(simulation%rt_process_model_coupler)) then
     option%io_buffer = 'The UFD_DECAY process model requires reactive &
          &transport.'
     call PrintErrMsg(option)
  endif

  pmc_ufd_decay => PMCThirdPartyCreate()
  call pmc_ufd_decay%SetName(pmc_name)
  call pmc_ufd_decay%SetOption(option)
  call pmc_ufd_decay%SetCheckpointOption(simulation%checkpoint_option)
  call pmc_ufd_decay%SetWaypointList(simulation%waypoint_list_subsurface)
  pmc_ufd_decay%pm_list => pm_ufd_decay
  pmc_ufd_decay%pm_ptr%pm => pm_ufd_decay
  pmc_ufd_decay%realization => realization

  ! set up logging stage
  string = 'UFD_DECAY'
  call LoggingCreateStage(string,pmc_ufd_decay%stage)
  call PMCBaseSetChildPeerPtr(PMCCastToBase(pmc_ufd_decay),PM_CHILD, &
         PMCCastToBase(simulation%rt_process_model_coupler), &
         pmc_dummy,PM_APPEND)

end subroutine AddPMCUDFDecay

! ************************************************************************** !

subroutine AddPMCUDFBiosphere(simulation,pm_ufd_biosphere,pmc_name,&
                              pm_ufd_decay_present,realization,input,option)

  !
  ! Adds a UFD biosphere PMC
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PMC_Base_class
  use PMC_Third_Party_class
  use PM_UFD_Biosphere_class
  use Realization_Subsurface_class
  use String_module
  use Option_module
  use Logging_module
  use Input_Aux_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_ufd_biosphere_type), pointer :: pm_ufd_biosphere
  character(len=*) :: pmc_name
  logical :: pm_ufd_decay_present
  class(realization_subsurface_type), pointer :: realization
  type(input_type), pointer :: input
  type(option_type), pointer :: option

  class(pmc_third_party_type), pointer :: pmc_ufd_biosphere
  character(len=MAXSTRINGLENGTH) :: string
  class(pmc_base_type), pointer :: pmc_dummy

  nullify(pmc_dummy)

  string = 'UFD_BIOSPHERE'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  call pm_ufd_biosphere%ReadPMBlock(input)
  if (.not.associated(simulation%rt_process_model_coupler)) then
     option%io_buffer = 'The UFD_BIOSPHERE process model requires reactive &
          &transport.'
     call PrintErrMsg(option)
  endif
  if (.not.pm_ufd_decay_present) then
     option%io_buffer = 'The UFD_BIOSPHERE process model requires the &
          &UFD_DECAY process model.'
     call PrintErrMsg(option)
  endif

  pmc_ufd_biosphere => PMCThirdPartyCreate()
  call pmc_ufd_biosphere%SetName(pmc_name)
  call pmc_ufd_biosphere%SetOption(option)
  call pmc_ufd_biosphere%SetCheckpointOption(simulation%checkpoint_option)
  call pmc_ufd_biosphere%SetWaypointList(simulation%waypoint_list_subsurface)
  pmc_ufd_biosphere%pm_list => pm_ufd_biosphere
  pmc_ufd_biosphere%pm_ptr%pm => pm_ufd_biosphere
  pmc_ufd_biosphere%realization => realization

  ! set up logging stage
  string = 'UFD_BIOSPHERE'
  call LoggingCreateStage(string,pmc_ufd_biosphere%stage)
  call PMCBaseSetChildPeerPtr(PMCCastToBase(pmc_ufd_biosphere),PM_CHILD, &
         PMCCastToBase(simulation%rt_process_model_coupler), &
         pmc_dummy,PM_APPEND)

end subroutine AddPMCUDFBiosphere

! ************************************************************************** !

subroutine AddPMCAuxiliary(simulation,pm_auxiliary,pmc_name, &
                           realization,option)

  !
  ! Adds an auxiliary PMC
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PMC_Base_class
  use PM_Auxiliary_class
  use PMC_Auxiliary_class
  use PMC_Subsurface_class
  use Realization_Subsurface_class
  use String_module
  use Option_module
  use Logging_module
  use Input_Aux_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_auxiliary_type), pointer :: pm_auxiliary
  character(len=*) :: pmc_name
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option

  class(pmc_auxiliary_type), pointer :: pmc_auxiliary
  character(len=MAXSTRINGLENGTH) :: string
  class(pmc_base_type), pointer :: pmc_dummy

  nullify(pmc_dummy)

  string = 'salinity'
  if (StringCompareIgnoreCase(pm_auxiliary%ctype,string)) then
    if (associated(simulation%rt_process_model_coupler)) then
      pmc_auxiliary => PMCAuxiliaryCreate()
      call pmc_auxiliary%SetName(pmc_name)
      call pmc_auxiliary%SetOption(option)
      pm_auxiliary%realization => realization
      pmc_auxiliary%pm_list => pm_auxiliary
      pmc_auxiliary%pm_aux => pm_auxiliary
      call PMCBaseSetChildPeerPtr(PMCCastToBase(pmc_auxiliary),PM_PEER, &
             PMCCastToBase(simulation%rt_process_model_coupler), &
             pmc_dummy,PM_APPEND)
    else
      option%io_buffer = 'Reactive transport must be included in the &
           &SIMULATION block in order to use the SALINITY process model.'
      call PrintErrMsg(option)
    endif
  endif

  call LoggingCreateStage(string,pmc_auxiliary%stage)

end subroutine AddPMCAuxiliary

! ************************************************************************** !

subroutine SubsurfInitCommandLineSettings(option)
  !
  ! Initializes PFLTORAN subsurface output
  ! filenames, etc.
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/13
  !

  use Option_module
  use Input_Aux_module

  implicit none

  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag

  string = '-multisimulation'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = MULTISIMULATION_SIM_TYPE
  endif

  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = STOCHASTIC_SIM_TYPE
  endif

end subroutine SubsurfInitCommandLineSettings

! ************************************************************************** !

subroutine SubsurfaceSetFlowMode(pm_flow,option)
  !
  ! Sets the flow mode (richards, vadose, mph, etc.)
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  use Option_module
  use PM_Subsurface_Flow_class
  use PM_Base_class
  use PM_Flash2_class
  use PM_General_class
  use PM_Hydrate_class
  use PM_WIPP_Flow_class
  use PM_Immis_class
  use PM_Miscible_class
  use PM_Mphase_class
  use PM_Richards_class
  use PM_TH_class
  use PM_TOilIms_class
  use PM_TOWG_class
  use PM_TOWG_Aux_module
  use PM_Richards_TS_class
  use PM_TH_TS_class
  use General_Aux_module

  implicit none

  type(option_type) :: option
  class(pm_subsurface_flow_type), pointer :: pm_flow

  if (.not.associated(pm_flow)) then
    option%nphase = 1
    option%liquid_phase = 1
    option%gas_phase = 2 ! still set gas phase to 2 for transport
    ! assume default isothermal when only transport
    option%use_isothermal = PETSC_TRUE
    return
  endif

  select type(pm_flow)
    class is (pm_flash2_type)
      option%iflowmode = FLASH2_MODE
      option%nphase = 2
      option%liquid_phase = 1
      option%gas_phase = 2
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
      option%use_isothermal = PETSC_FALSE
      option%water_id = 1
      option%air_id = 2
    class is (pm_wippflo_type)
      option%iflowmode = WF_MODE
      option%nphase = 2
      option%liquid_phase = 1  ! liquid_pressure
      option%gas_phase = 2     ! gas_pressure

      option%capillary_pressure_id = 3
      option%saturation_pressure_id = 4

      option%water_id = 1
      option%air_id = 2

      option%nflowdof = 2
      option%nflowspec = 2
    class is (pm_general_type)
      option%iflowmode = G_MODE
      option%nphase = 2
      option%liquid_phase = 1  ! liquid_pressure
      option%gas_phase = 2     ! gas_pressure

      option%air_pressure_id = 3
      option%capillary_pressure_id = 4
      option%vapor_pressure_id = 5
      option%saturation_pressure_id = 6

      option%water_id = 1
      option%air_id = 2
      option%energy_id = 3

      option%nflowdof = 3
      option%nflowspec = 2
      option%use_isothermal = PETSC_FALSE
    class is (pm_hydrate_type)
      call PMHydrateSetFlowMode(option)
    class is (pm_toil_ims_type)
      option%iflowmode = TOIL_IMS_MODE
      option%nphase = 2
      option%liquid_phase = 1           ! liquid_pressure
      option%oil_phase = 2              ! oil_pressure

      option%phase_map(1) = LIQUID_PHASE
      option%phase_map(2) = OIL_PHASE

! Check that MAX_PHASE is sufficiently large

      if( option%nphase > MAX_PHASE ) then
        option%io_buffer = 'ERROR: number of phases has exceeded MAX_PHASE'
        call PrintMsg(option)
      endif

      !option%capillary_pressure_id = 3  ! capillary pressure

      option%nflowdof = 3
      !two species (H2O,OIL): each present only in its own rich phase
      option%nflowspec = 2

      option%water_id = 1
      option%oil_id = 2
      option%energy_id = 3

      option%use_isothermal = PETSC_FALSE
    class is (pm_towg_type)
      option%iflowmode = TOWG_MODE

! Basic oil and gas

      option%nphase = 3
      option%liquid_phase = 1           ! liquid_pressure
      option%oil_phase = 2              ! oil_pressure
      option%gas_phase = 3              ! gas_pressure

!  Add solvent phase

      if (towg_miscibility_model == TOWG_SOLVENT_TL ) then
        option%nphase = 4
        option%solvent_phase = 4         ! solvent saturation
      endif

!  Phase maps

      option%phase_map(1) = LIQUID_PHASE
      option%phase_map(2) = OIL_PHASE
      option%phase_map(3) = GAS_PHASE

! Add solvent phase

      if ( towg_miscibility_model == TOWG_SOLVENT_TL ) then
        option%phase_map(4) = SOLVENT_PHASE
      endif

! Check that MAX_PHASE is sufficiently large

      if( option%nphase > MAX_PHASE ) then
        option%io_buffer = 'ERROR: number of phases has exceeded MAX_PHASE'
        call PrintMsg(option)
      endif

      option%energy_id = towg_energy_eq_idx
      select case (towg_miscibility_model)
        case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF,TOWG_BLACK_OIL)
          option%nflowdof = 4
          option%nflowspec = 3 !H20, Oil, Gas
        case(TOWG_SOLVENT_TL)
          option%nphase = 4
          option%nflowdof = 5
          option%nflowspec = 4 !H20, Oil, Gas, Solvent
        case default
          !option%io_buffer = 'SubsurfaceSetFlowMode: ' //
          !  'towg_miscibility_model must be intiialized'
          option%io_buffer = 'only immiscible TOWG currently implemented'
          call PrintErrMsg(option)
      end select
      option%use_isothermal = PETSC_FALSE
    class is (pm_immis_type)
      option%iflowmode = IMS_MODE
      option%nphase = 2
      option%liquid_phase = 1
      option%gas_phase = 2
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
      option%io_buffer = 'Material AuxVars must be refactored for IMMIS.'
      call PrintErrMsg(option)
    class is (pm_miscible_type)
      option%iflowmode = MIS_MODE
      option%nphase = 1
      option%liquid_phase = 1
      option%gas_phase = 2
      option%nflowdof = 2
      option%nflowspec = 2
      option%io_buffer = 'Material AuxVars must be refactored for MISCIBLE.'
      call PrintErrMsg(option)
    class is (pm_mphase_type)
      option%iflowmode = MPH_MODE
      option%nphase = 2
      option%liquid_phase = 1
      option%gas_phase = 2
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2 ! read CO2DATA0.dat
!     option%itable = 1 ! create CO2 database: co2data.dat
      option%use_isothermal = PETSC_FALSE
      option%water_id = 1
      option%air_id = 2
    class is (pm_richards_type)
      option%iflowmode = RICHARDS_MODE
      option%nphase = 1
      option%liquid_phase = 1
      option%nflowdof = 1
      option%nflowspec = 1
      option%use_isothermal = PETSC_TRUE
    class is (pm_th_type)
      option%iflowmode = TH_MODE
      option%nphase = 1
      option%liquid_phase = 1
      option%gas_phase = 2
      option%nflowdof = 2
      option%nflowspec = 1
      option%use_isothermal = PETSC_FALSE
      option%flow%store_fluxes = PETSC_TRUE
    class is (pm_richards_ts_type)
      option%iflowmode = RICHARDS_TS_MODE
      option%nphase = 1
      option%liquid_phase = 1
      option%nflowdof = 1
      option%nflowspec = 1
      option%use_isothermal = PETSC_TRUE
    class is (pm_th_ts_type)
      option%iflowmode = TH_TS_MODE
      option%nphase = 1
      option%liquid_phase = 1
      option%gas_phase = 2
      option%nflowdof = 2
      option%nflowspec = 1
      option%use_isothermal = PETSC_FALSE
      option%flow%store_fluxes = PETSC_TRUE
    class default
      option%io_buffer = ''
      call PrintErrMsg(option)

  end select

end subroutine SubsurfaceSetFlowMode

! ************************************************************************** !

subroutine SubsurfaceReadFlowPM(input,option,pm)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PMC_Base_class

  use PM_Base_class
  use PM_Flash2_class
  use PM_General_class
  use PM_Hydrate_class
  use PM_WIPP_Flow_class
  use PM_Immis_class
  use PM_Miscible_class
  use PM_Mphase_class
  use PM_Richards_class
  use PM_TH_class
  use PM_TOilIms_class
  use PM_TOWG_class
  use PM_Richards_TS_class
  use PM_TH_TS_class
  use Init_Common_module
  use General_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,SUBSURFACE_FLOW'

  nullify(pm)
  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('MODE')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'mode',error_string)
        call StringToUpper(word)
        select case(word)
          case('GENERAL','HYDRATE','TOIL_IMS','TOWG_IMMISCIBLE', &
               'TODD_LONGSTAFF','TOWG_MISCIBLE','BLACK_OIL', &
               'SOLVENT_TL','WIPP_FLOW')
          ! In OptionFlowInitRealization(), numerical_derivatives is set to
          ! PETSC_FALSE, but the default for GENERAL needs to be PETSC_TRUE.
          ! This is will eventually affect all flow modes with numerical
          ! derivatives as default if analytical derivatives are available
          ! and we are keying off this flag.
          option%flow%numerical_derivatives = PETSC_TRUE
        end select
        select case(word)
          case('GENERAL')
            pm => PMGeneralCreate()
          case('HYDRATE')
            pm => PMHydrateCreate()
          case('WIPP_FLOW')
            pm => PMWIPPFloCreate()
          case('BRAGFLO')
            option%io_buffer = 'BRAGFLO mode has been merged with WIPP_FLOW. &
              &Please use WIPP_FLOW instead.'
            call PrintErrMsg(option)
          case('MPHASE')
            pm => PMMphaseCreate()
          case('FLASH2')
            pm => PMFlash2Create()
          case('IMS','IMMIS','THS')
            pm => PMImmisCreate()
          case('MIS','MISCIBLE')
            pm => PMMiscibleCreate()
          case('RICHARDS')
            pm => PMRichardsCreate()
          case('TH')
            pm => PMTHCreate()
          case('TOIL_IMS')
            pm => PMTOilImsCreate()
          !case('TOWG')
          case('TOWG_IMMISCIBLE','TODD_LONGSTAFF','TOWG_MISCIBLE', &
               'BLACK_OIL','SOLVENT_TL')
            pm => PMTOWGCreate(input,word,option)
          case ('RICHARDS_TS')
            pm => PMRichardsTSCreate()
          case ('TH_TS')
            pm => PMTHTSCreate()
          case default
            error_string = trim(error_string) // ',MODE'
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
        pm%option => option
      case('OPTIONS')
        if (.not.associated(pm)) then
          option%io_buffer = 'MODE keyword must be read first under ' // &
                             trim(error_string)
          call PrintErrMsg(option)
        endif
        call pm%ReadSimulationBlock(input)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (.not.associated(pm)) then
    option%io_buffer = 'A flow MODE (card) must be included in the &
      &SUBSURFACE_FLOW block in ' // trim(error_string) // '.'
    call PrintErrMsg(option)
  endif

end subroutine SubsurfaceReadFlowPM

! ************************************************************************** !

subroutine SubsurfaceReadRTPM(input,option,pm)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PMC_Base_class
  use PM_Base_class
  use PM_RT_class

  use Init_Common_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,SUBSURFACE_TRANSPORT'

  pm => PMRTCreate()
  pm%option => option
  option%itranmode = RT_MODE

  call pm%ReadSimulationBlock(input)

end subroutine SubsurfaceReadRTPM

! ************************************************************************** !

subroutine SubsurfaceReadNWTPM(input,option,pm)
  !
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PMC_Base_class
  use PM_Base_class
  use PM_NWT_class

  use Init_Common_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,NUCLEAR_WASTE_TRANSPORT'

  pm => PMNWTCreate()
  pm%option => option
  option%itranmode = NWT_MODE

  call pm%ReadSimulationBlock(input)

end subroutine SubsurfaceReadNWTPM

! ************************************************************************** !

subroutine SubsurfaceReadWasteFormPM(input,option,pm)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_Waste_Form_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'SIMULATION,PROCESS_MODELS,WASTE_FORM'

  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)

    found = PETSC_FALSE
    call PMBaseReadSelectCase(pm,input,word,found,error_string,option)
    if (found) cycle

    select case(word)
      case('TYPE')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'mode',error_string)
        call StringToUpper(word)
        select case(word)
          case('GENERAL')
            pm => PMWFCreate()
          case default
            option%io_buffer = 'WASTE FORM type ' // trim(word) // &
              ' not recognized. Only TYPE GENERAL currently supported. &
              & TYPE GLASS or TYPE FMDM no longer supported.'
            call PrintErrMsg(option)
        end select
      case default
        option%io_buffer = 'Keyword ' // trim(word) // &
              ' not recognized for the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (.not.associated(pm)) then
    option%io_buffer = 'TYPE card missing in ' // trim(error_string)
    call PrintErrMsg(option)
  endif

  pm%option => option

end subroutine SubsurfaceReadWasteFormPM

! ************************************************************************** !

subroutine SubsurfaceReadUFDDecayPM(input,option,pm)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_UFD_Decay_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'SIMULATION,PROCESS_MODELS,UFD_DECAY'

  pm => PMUFDDecayCreate()
  pm%option => option

  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)

    found = PETSC_FALSE
    call PMBaseReadSelectCase(pm,input,word,found,error_string,option)
    if (found) cycle

    select case(word)
      case default
        option%io_buffer = 'Keyword ' // trim(word) // &
              ' not recognized for the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine SubsurfaceReadUFDDecayPM

! ************************************************************************** !

subroutine SubsurfaceReadUFDBiospherePM(input,option,pm)
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_UFD_Biosphere_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,UFD_BIOSPHERE'

  pm => PMUFDBCreate()
  pm%option => option

  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case default
        option%io_buffer = 'Keyword ' // trim(word) // &
              ' not recognized for the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine SubsurfaceReadUFDBiospherePM

! ************************************************************************** !

subroutine SubsurfaceInitSimulation(simulation)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Realization_Subsurface_class
  use Realization_Base_class
  use Discretization_module
  use Option_module
  use Output_module, only : Output
  use Output_Aux_module
  use Global_module
  use Init_Subsurface_module
  use Init_Subsurface_Flow_module
  use Init_Subsurface_Tran_module
  use Init_Common_module
  use Waypoint_module
  use Strata_module
  use Regression_module
  use PMC_Subsurface_class
  use PMC_Auxiliary_class
  use PMC_Base_class
  use PM_Base_class
  use PM_Base_Pointer_module
  use PM_Subsurface_Flow_class
  use PM_Auxiliary_class
  use Timestepper_BE_class
  use Waypoint_module

  implicit none

  class(simulation_subsurface_type) :: simulation

  class(pmc_subsurface_type), pointer :: flow_process_model_coupler
  class(pmc_subsurface_type), pointer :: tran_process_model_coupler
  class(pmc_auxiliary_type), pointer :: auxiliary_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pm_base_type), pointer :: cur_process_model
  class(pm_auxiliary_type), pointer :: pm_aux

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(waypoint_list_type), pointer :: sync_waypoint_list
  character(len=MAXSTRINGLENGTH) :: string
  SNESLineSearch :: linesearch
  PetscInt :: ndof
  PetscBool, allocatable :: dof_is_active(:)
  PetscBool :: failure
  PetscErrorCode :: ierr

  realization => simulation%realization
  option => realization%option

! begin from old Init()
  call SubsurfaceSetupRealization(simulation)
  call InitCommonAddOutputWaypoints(option,simulation%output_option, &
                                    simulation%waypoint_list_subsurface)

  !TODO(geh): refactor
  if (associated(simulation%flow_process_model_coupler)) then
    if (associated(simulation%flow_process_model_coupler%timestepper)) then
      simulation%flow_process_model_coupler%timestepper%cur_waypoint => &
        simulation%waypoint_list_subsurface%first
    endif
  endif
  if (associated(simulation%rt_process_model_coupler)) then
    if (associated(simulation%rt_process_model_coupler%timestepper)) then
      simulation%rt_process_model_coupler%timestepper%cur_waypoint => &
        simulation%waypoint_list_subsurface%first
    endif
  endif
  if (associated(simulation%nwt_process_model_coupler)) then
    if (associated(simulation%nwt_process_model_coupler%timestepper)) then
      simulation%nwt_process_model_coupler%timestepper%cur_waypoint => &
        simulation%waypoint_list_subsurface%first
    endif
  endif

  !TODO(geh): refactor
  ! initialize global auxiliary variable object
  call GlobalSetup(realization)

  ! always call the flow side since a velocity field still has to be
  ! set if no flow exists
  call InitSubsurfFlowSetupRealization(realization)
  if (option%ntrandof > 0) then
    call InitSubsurfTranSetupRealization(realization)
  endif
  ! InitSubsurfaceSetupZeroArray must come after InitSubsurfaceXXXRealization
  call InitSubsurfaceSetupZeroArrays(realization)
  call OutputVariableAppendDefaults(realization%output_option% &
                                      output_snap_variable_list,option)

  call RegressionCreateMapping(simulation%regression,realization)
! end from old Init()

  call DiscretizationPrintInfo(realization%discretization, &
                               realization%patch%grid,option)

  !----------------------------------------------------------------------------!
  ! This section for setting up new process model approach
  !----------------------------------------------------------------------------!

  if (StrataEvolves(realization%patch%strata_list)) then
    allocate(pm_aux)
    call PMAuxiliaryInit(pm_aux)
    string = 'EVOLVING_STRATA'
    call PMAuxiliarySetFunctionPointer(pm_aux,string)
    pm_aux%realization => realization
    pm_aux%option => option

    auxiliary_process_model_coupler => PMCAuxiliaryCreate()
    auxiliary_process_model_coupler%pm_list => pm_aux
    auxiliary_process_model_coupler%pm_aux => pm_aux
    auxiliary_process_model_coupler%option => option
    ! place the material process model as %peer for the top pmc
    simulation%process_model_coupler_list%peer => &
      auxiliary_process_model_coupler
  endif

  ! For each ProcessModel, set:
  ! - realization (subsurface or surface),
  ! - stepper (flow/trans/surf_flow),
  ! For each ProcessModelCoupler, set:
  ! - SNES functions (Residual/Jacobain), or TS function (RHSFunction)

  cur_process_model_coupler_top => simulation%process_model_coupler_list

  ! the following recursive subroutine will also call each pmc child
  ! and each pms's peers
  if (associated(cur_process_model_coupler_top)) then
    call SetUpPMApproach(cur_process_model_coupler_top,simulation)
  endif

  ! point the top process model coupler to Output
  simulation%process_model_coupler_list%Output => Output

  ! setup the waypoint list
  call SetupWaypointList(simulation)

  if (realization%debug%print_couplers) then
    call InitCommonVerifyAllCouplers(realization)
  endif

  call SubsurfaceJumpStart(simulation)

end subroutine SubsurfaceInitSimulation

! ************************************************************************** !

recursive subroutine SetUpPMApproach(pmc,simulation)
!
! Loops through all of the PMC's recursively and sets their realization,
! timestepper, and solver.
!
! Author: Jenn Frederick, SNL
! Date: 04/04/2016
!
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PMC_Base_class
  use PMC_Subsurface_class
  use PM_Base_Pointer_module
  use PM_Base_class
  use PM_Subsurface_Flow_class
  !TODO(geh): are these needed
  use PM_General_class
  use PM_Hydrate_class
  use PM_WIPP_Flow_class
  use PM_Richards_class
  use PM_TH_class
  use PM_RT_class
  use PM_NWT_class
  use PM_Waste_Form_class
  use PM_WIPP_SrcSink_class
  use PM_UFD_Decay_class
  use PM_UFD_Biosphere_class
  use PM_TOilIms_class
  use Option_module
  use Simulation_Subsurface_class
  use Realization_Subsurface_class
  use Timestepper_BE_class

  implicit none


  class(pmc_base_type), pointer :: pmc
  class(simulation_subsurface_type) :: simulation

  class(realization_subsurface_type), pointer :: realization
  class(pm_base_type), pointer :: cur_pm
  type(option_type), pointer :: option
  SNESLineSearch :: linesearch
  PetscErrorCode :: ierr

  realization => simulation%realization
  option => realization%option

  if (.not.associated(pmc)) return

  pmc%waypoint_list => simulation%waypoint_list_subsurface

  ! loop through this pmc's process models:
  cur_pm => pmc%pm_list
  do
    if (.not.associated(cur_pm)) exit
    ! set realization
    select type(cur_pm)
      class is(pm_rt_type)
        if (.not.associated(realization%reaction)) then
          option%io_buffer = 'SUBSURFACE_TRANSPORT specified as a &
            &process model without a corresponding CHEMISTRY block.'
          call PrintErrMsg(option)
        endif
        call cur_pm%SetRealization(realization)
        
      class is(pm_nwt_type)
        if (.not.associated(realization%reaction_nw)) then
          option%io_buffer = 'NUCLEAR_WASTE_TRANSPORT is specified as a &
            &process model in the SIMULATION block without the corresponding &
            &NUCLEAR_WASTE_CHEMISTRY block within the SUBSURFACE block.'
          call PrintErrMsg(option)
        endif
        call cur_pm%SetRealization(realization)

      class is(pm_subsurface_flow_type)
        call cur_pm%SetRealization(realization)

      class is(pm_waste_form_type)
        call cur_pm%SetRealization(realization)

      class is(pm_ufd_decay_type)
        call cur_pm%SetRealization(realization)

      class is(pm_ufd_biosphere_type)
        call cur_pm%SetRealization(realization)

    end select

    ! set time stepper
    select type(cur_pm)
      class is(pm_subsurface_flow_type)
        pmc%timestepper%dt = option%flow_dt

      class is(pm_rt_type)
        pmc%timestepper%dt = option%tran_dt
        
      class is(pm_nwt_type)
        pmc%timestepper%dt = option%tran_dt

    end select
    cur_pm%output_option => simulation%output_option
    call cur_pm%Setup()
    cur_pm => cur_pm%next
  enddo

  call pmc%SetupSolvers()

  ! call this function for this pmc's child
  if (associated(pmc%child)) then
    call SetUpPMApproach(pmc%child,simulation)
  endif

  ! call this function for this pmc's peer
  if (associated(pmc%peer)) then
    call SetUpPMApproach(pmc%peer,simulation)
  endif


end subroutine SetUpPMApproach

! ************************************************************************** !

subroutine SubsurfaceSetupRealization(simulation)
  !
  ! Initializes material property data structres and assign them to the domain.
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  !
  use Init_Subsurface_module
  use Simulation_Subsurface_class
  use Realization_Subsurface_class
  use Option_module
  use Logging_module
  use Waypoint_module
  use Init_Common_module
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  use Reaction_Database_module
  use EOS_module
  use Dataset_module
  use Patch_module
  use EOS_module !to be removed as already present above

  implicit none

  class(simulation_subsurface_type) :: simulation

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  PetscErrorCode :: ierr

  realization => simulation%realization
  option => realization%option

  call PetscLogEventBegin(logging%event_setup,ierr);CHKERRQ(ierr)

  ! process eos tables ready for evaluation
  call EOSProcess(option)

  ! set reference densities if not specified in input file.
  call EOSReferenceDensity(option)


  ! read reaction database
  if (associated(realization%reaction)) then
    if (realization%reaction%use_full_geochemistry) then
        call DatabaseRead(realization%reaction,option)
        call BasisInit(realization%reaction,option)
    else
      ! turn off activity coefficients since the database has not been read
      realization%reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
      ! jenn:todo Should I turn on print here too?
      allocate(realization%reaction%primary_species_print(option%ntrandof))
      realization%reaction%primary_species_print = PETSC_TRUE
    endif
  endif

  ! SK 09/30/13, Added to check if Mphase is called with OS
  if (option%transport%reactive_transport_coupling == OPERATOR_SPLIT .and. &
      option%iflowmode == MPH_MODE) then
    option%io_buffer = 'Operator split not implemented with MPHASE. &
                       &Switching to Global Implicit.'
    !geh: We should force the user to switch without automatically switching
!    call PrintWrnMsg(option)
    call PrintErrMsg(option)
    option%transport%reactive_transport_coupling = GLOBAL_IMPLICIT
  endif

  ! create grid and allocate vectors
  call RealizationCreateDiscretization(realization)

  ! read any regions provided in external files
  call InitCommonReadRegionFiles(realization)
  ! clip regions and set up boundary connectivity, distance
  call RealizationLocalizeRegions(realization)
  call RealizationPassPtrsToPatches(realization)
  call RealizationProcessDatasets(realization)
  if (realization%output_option%mass_balance_region_flag) then
    call PatchGetCompMassInRegionAssign(realization%patch%region_list, &
         realization%output_option%mass_balance_region_list,option)
  endif
  ! link conditions with regions through couplers and generate connectivity
  call RealProcessMatPropAndSatFunc(realization)
  ! must process conditions before couplers in order to determine dataset types
  call RealizationProcessConditions(realization)
  call RealizationProcessCouplers(realization)
  call SubsurfSandboxesSetup(realization)
  call RealProcessFluidProperties(realization)
  call SubsurfInitMaterialProperties(realization)
  ! assignVolumesToMaterialAuxVars() must be called after
  ! RealizInitMaterialProperties() where the Material object is created
  call SubsurfAssignVolsToMatAuxVars(realization)
  call RealizationInitAllCouplerAuxVars(realization)
  if (option%ntrandof > 0) then
    call PrintMsg(option,"  Setting up TRAN Realization ")
    call RealizationInitConstraints(realization)
    call PrintMsg(option,"  Finished setting up TRAN Realization ")
  endif
  call RealizationPrintCouplers(realization)
  if (.not.option%steady_state) then
    ! add waypoints associated with boundary conditions, source/sinks etc. to list
    call RealizationAddWaypointsToList(realization, &
                                       simulation%waypoint_list_subsurface)
    ! fill in holes in waypoint data
  endif
  call PetscLogEventEnd(logging%event_setup,ierr);CHKERRQ(ierr)

#ifdef OS_STATISTICS
  call RealizationPrintGridStatistics(realization)
#endif

#if !defined(HDF5_BROADCAST)
  call PrintMsg(option,"Default HDF5 method is used in Initialization")
#else
  call PrintMsg(option,"Glenn's HDF5 broadcast method is used in Initialization")
#endif

end subroutine SubsurfaceSetupRealization

! ************************************************************************** !

subroutine SetupWaypointList(simulation)
  ! 
  ! Sets up waypoint list
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  ! 

  use Checkpoint_module
  use Realization_Subsurface_class
  use Option_module
  use Waypoint_module

  implicit none

  class(simulation_subsurface_type) :: simulation

  class(realization_subsurface_type), pointer :: realization
  type(waypoint_list_type), pointer :: sync_waypoint_list
  type(option_type), pointer :: option

  realization => simulation%realization
  option => realization%option

  ! create sync waypoint list to be used a few lines below
  sync_waypoint_list => &
    WaypointCreateSyncWaypointList(simulation%waypoint_list_subsurface)

  ! merge in outer waypoints (e.g. checkpoint times)
  call WaypointListCopyAndMerge(simulation%waypoint_list_subsurface, &
                                simulation%waypoint_list_outer,option)

  ! add sync waypoints into outer list
  call WaypointListMerge(simulation%waypoint_list_outer,sync_waypoint_list, &
                         option)

  ! add in periodic time waypoints for checkpointing. these will not appear
  ! in the outer list
  call CheckpointPeriodicTimeWaypoints(simulation%checkpoint_option, &
                                       simulation%waypoint_list_subsurface, &
                                       option)

  ! clean up waypoints
  if (.not.option%steady_state) then
   ! fill in holes in waypoint data
    call WaypointListFillIn(simulation%waypoint_list_subsurface,option)
    call WaypointListRemoveExtraWaypnts(simulation%waypoint_list_subsurface, &
                                        option)
  endif

  ! debugging output
  if (realization%debug%print_waypoints) then
    call WaypointListPrint(simulation%waypoint_list_subsurface,option, &
                           realization%output_option)
  endif

end subroutine SetupWaypointList

! ************************************************************************** !

subroutine SubsurfaceJumpStart(simulation)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !

  use Realization_Subsurface_class
  use Option_module
  use Reactive_Transport_module, only : RTJumpStartKineticSorption

  implicit none

  type(simulation_subsurface_type) :: simulation

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: failure
  PetscErrorCode :: ierr

  realization => simulation%realization
  option => realization%option

  call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                           PETSC_NULL_CHARACTER, "-vecload_block_size", &
                           failure, ierr);CHKERRQ(ierr)

  if (option%transport%jumpstart_kinetic_sorption .and. &
      option%time < 1.d-40) then
    ! only user jumpstart for a restarted simulation
    if (.not. option%restart_flag) then
      option%io_buffer = 'Only use JUMPSTART_KINETIC_SORPTION on a &
        &restarted simulation.  ReactionEquilibrateConstraint() will &
        &appropriately set sorbed initial concentrations for a normal &
        &(non-restarted) simulation.'
      call PrintErrMsg(option)
    endif
    call RTJumpStartKineticSorption(realization)
  endif

end subroutine SubsurfaceJumpStart

! ************************************************************************** !

subroutine SubsurfaceReadRequiredCards(simulation,input)
  !
  ! Reads required cards from input file
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07, refactored 08/20/14, refactored 12/10/14
  !

  use Option_module
  use Discretization_module
  use Grid_module
  use Input_Aux_module
  use String_module
  use Patch_module
  use Realization_Subsurface_class
  use HDF5_Aux_module

  use Simulation_Subsurface_class
  use General_module
  use Reaction_module
  use Reaction_Aux_module
  use NW_Transport_Aux_module
  use Init_Common_module

  use Grid_Grdecl_module, only : SetUGrdEclCmplLocation

  implicit none

  class(simulation_subsurface_type) :: simulation

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  type(patch_type), pointer :: patch, patch2
  type(grid_type), pointer :: grid
  class(realization_subsurface_type), pointer :: realization
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  PetscInt :: ci,cj,ck,ckl,cku,ckll,ckuu
  PetscBool :: found
  PetscBool,parameter::cijk_d_true =PETSC_TRUE
  PetscBool,parameter::cijk_d_false=PETSC_FALSE
  PetscBool :: qerr

  character(len = MAXSTRINGLENGTH) :: wname

  realization => simulation%realization
  patch => realization%patch
  option => realization%option
  discretization => realization%discretization

  qerr  = PETSC_FALSE
  wname = '<missing>'
  found = PETSC_FALSE

  call InputPushBlock(input,'SUBSURFACE',option)

! Read in select required cards
!.........................................................................

!  Search for 'WELL_DATA' section and read well locations if found

  string = "WELL_DATA"
  call InputFindStringInFile(input,option,string,PETSC_FALSE,found)
  if( found ) then
    call InputPushBlock(input,'WELL_DATA',option)

! Read the WELL_DATA information

    call InputReadWord(input,option,wname,PETSC_TRUE)

! Search for completion locations

    call InputPushBlock(input,option)
    do
      call InputReadPflotranString(input,option)
      if (InputError(input)) exit

      call InputReadCard(input,option,word,PETSC_FALSE)
      call StringToUpper(word)
      card = trim(word)

      select case(trim(card))
        case('CIJK','CIJK_Z')
          ci  = 1
          cj  = 1
          ckl = 1
          cku = 1
          call InputReadInt(input,option,ci)
          call InputErrorMsg(input,option,'cijk I','WELL_DATA')
          call InputReadInt(input,option,cj)
          call InputErrorMsg(input,option,'cijk I','WELL_DATA')
          call InputReadInt(input,option,ckl)
          call InputErrorMsg(input,option,'cijk KL','WELL_DATA')
          call InputReadInt(input,option,cku)
          call InputErrorMsg(input,option,'cijk KU','WELL_DATA')
          ckll=min(ckl,cku)
          ckuu=max(ckl,cku)
          do ck=ckll,ckuu
            call SetUGrdEclCmplLocation(wname,ci,cj,ck,cijk_d_false,qerr)
          enddo
          if( qerr ) then
            input%ierr = 1
            call InputErrorMsg(input,option,'cijk','same well more than once')
          endif
        case('CIJK_D')
          ci  = 1
          cj  = 1
          ckl = 1
          cku = 1
          call InputReadInt(input,option,ci)
          call InputErrorMsg(input,option,'cijk_d I','WELL_DATA')
          call InputReadInt(input,option,cj)
          call InputErrorMsg(input,option,'cijk_d I','WELL_DATA')
          call InputReadInt(input,option,ckl)
          call InputErrorMsg(input,option,'cijk_d KL','WELL_DATA')
          call InputReadInt(input,option,cku)
          call InputErrorMsg(input,option,'cijk_d KU','WELL_DATA')
          ckll=min(ckl,cku)
          ckuu=max(ckl,cku)
          do ck=ckll,ckuu
            call SetUGrdEclCmplLocation(wname,ci,cj,ck,cijk_d_true,qerr)
          enddo
          if( qerr ) then
            input%ierr = 1
            call InputErrorMsg(input,option,'cijk_d','same well more than once')
          endif
        case('WELL_DATA') ! May be more than one
          call InputReadWord(input,option,wname,PETSC_TRUE)
      end select
    enddo
    call InputPopBlock(input,option)
    call InputPopBlock(input,option) ! for WELL_DATA
  endif

  ! GRID information - GRID is a required card for every simulation
  string = "GRID"
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call InputPushBlock(input,'GRID',option)
  call DiscretizationReadRequiredCards(discretization,input,option)

  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      patch => PatchCreate()
      patch%grid => discretization%grid
      if (.not.associated(realization%patch_list)) then
        realization%patch_list => PatchCreateList()
      endif
      call PatchAddToList(patch,realization%patch_list)
      realization%patch => patch
  end select
  call InputPopBlock(input,option)

  ! optional required cards - yes, an oxymoron, but we need to know if
  ! these exist before we can go any further.
  call InputRewind(input)
  call InputPushBlock(input,'REQUIRED_CARDS',option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    ! do not use InputReadCard here as this is a search operation
    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    select case(trim(card))

!....................
      case('DBASE_FILENAME')
        call InputPushCard(input,card,option)
        call InputReadFilename(input,option,string)
        call InputErrorMsg(input,option,'filename','DBASE_FILENAME')
        if (index(string,'.h5') > 0) then
          call HDF5ReadDbase(string,option)
        else
          call InputReadASCIIDbase(string,option)
        endif

!....................
      case('HDF5_WRITE_GROUP_SIZE')
        call InputPushCard(input,card,option)
        call InputReadInt(input,option,option%hdf5_write_group_size)
        call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')
        call InputSkipToEnd(input,option,'HDF5_WRITE_GROUP_SIZE')

      case('HDF5_READ_GROUP_SIZE')
        call InputPushCard(input,card,option)
        call InputReadInt(input,option,option%hdf5_read_group_size)
        call InputErrorMsg(input,option,'HDF5_READ_GROUP_SIZE','Group size')

!....................
      case('PROC')
        call InputPushCard(input,card,option)
        ! processor decomposition
        if (realization%discretization%itype == STRUCTURED_GRID) then
          grid => realization%patch%grid
          ! strip card from front of string
          call InputReadInt(input,option,grid%structured_grid%npx)
          call InputDefaultMsg(input,option,'npx')
          call InputReadInt(input,option,grid%structured_grid%npy)
          call InputDefaultMsg(input,option,'npy')
          call InputReadInt(input,option,grid%structured_grid%npz)
          call InputDefaultMsg(input,option,'npz')

          if (option%myrank == option%io_rank .and. &
              option%print_to_screen) then
            option%io_buffer = ' Processor Decomposition:'
            call PrintMsg(option)
            write(option%io_buffer,'("  npx   = ",3x,i4)') &
              grid%structured_grid%npx
            call PrintMsg(option)
            write(option%io_buffer,'("  npy   = ",3x,i4)') &
              grid%structured_grid%npy
            call PrintMsg(option)
            write(option%io_buffer,'("  npz   = ",3x,i4)') &
              grid%structured_grid%npz
            call PrintMsg(option)
          endif

          if (option%mycommsize /= grid%structured_grid%npx * &
                                 grid%structured_grid%npy * &
                                 grid%structured_grid%npz) then
            write(option%io_buffer,*) 'Incorrect number of processors &
              &specified: ',grid%structured_grid%npx*grid%structured_grid%npy* &
              grid%structured_grid%npz,' commsize = ',option%mycommsize
            call PrintErrMsg(option)
          endif
        endif
        
!....................
      case('CHEMISTRY')
        call InputPushCard(input,card,option)
        if (.not.associated(simulation%rt_process_model_coupler)) then
          option%io_buffer = 'CHEMISTRY card included when no &
            &SUBSURFACE_TRANSPORT process model included in SIMULATION block.'
          call PrintErrMsg(option)
        endif
        !geh: for some reason, we need this with CHEMISTRY read for
        !     multicontinuum
 !       option%use_mc = PETSC_TRUE
        call ReactionInit(realization%reaction,input,option)  
        realization%reaction_base => realization%reaction
        
!....................
      case('NUCLEAR_WASTE_CHEMISTRY')
        call InputPushCard(input,card,option)
        if (.not.associated(simulation%nwt_process_model_coupler)) then
          option%io_buffer = 'NUCLEAR_WASTE_CHEMISTRY card is &
            &included, but no NUCLEAR_WASTE_TRANSPORT process model found &
            &in the SIMULATION block.'
          call PrintErrMsg(option)
        endif     
        realization%reaction_nw => NWTReactionCreate()
        realization%reaction_base => realization%reaction_nw
        call NWTRead(realization%reaction_nw,input,option)
        
    end select
  enddo
  call InputPopBlock(input,option) ! REQUIRED_CARDS
  call InputPopBlock(input,option) ! SUBSURFACE

end subroutine SubsurfaceReadRequiredCards

! ************************************************************************** !

subroutine SubsurfaceReadInput(simulation,input)
  !
  ! Reads pflow input file
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !

  use Option_module
  use Field_module
  use Grid_module
  use Grid_Unstructured_Aux_module
  use Grid_Structured_module
  use Solver_module
  use Material_module
  use Saturation_Function_module
  use Characteristic_Curves_module
  use Creep_Closure_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Dataset_module
  use Dataset_Common_HDF5_class
  use Fluid_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Region_module
  use Condition_module
  use Transport_Constraint_Base_module
  use Transport_Constraint_RT_module
  use Transport_Constraint_NWT_module
  use Transport_Constraint_module
  use Coupler_module
  use Strata_module
  use Observation_module
  use Integral_Flux_module
  use Waypoint_module
  use Debug_module
  use Patch_module
  use Reaction_module
  use Reaction_Aux_module
  use NW_Transport_module
  use NW_Transport_Aux_module
  use Discretization_module
  use Input_Aux_module
  use String_module
  use Units_module
  use Reaction_Mineral_module
  use Regression_module
  use Output_Aux_module
  use Output_module
  use Output_Tecplot_module
  use Data_Mediator_Dataset_class
  use EOS_module
  use EOS_Water_module
  use SrcSink_Sandbox_module
  use Klinkenberg_module
  use WIPP_module
  use Utility_module
  use Checkpoint_module
  use Simulation_Subsurface_class
  use PMC_Subsurface_class
  use Timestepper_BE_class
  use Timestepper_Steady_class
  use Timestepper_TS_class
  use Well_Data_class
  use Hydrate_module

  use TH_Aux_module

#ifdef SOLID_SOLUTION
  use Reaction_Solid_Solution_module, only : SolidSolutionReadFromInputFile
#endif

  implicit none

  class(simulation_subsurface_type) :: simulation

  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXSTRINGLENGTH) :: string, temp_string
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: error_string

  character(len=1) :: backslash
  PetscReal :: temp_real, temp_real2
  PetscReal, pointer :: temp_real_array(:)
  PetscInt :: temp_int
  PetscInt :: id

  PetscBool :: vel_cent
  PetscBool :: vel_face
  PetscBool :: fluxes
  PetscBool :: mass_flowrate
  PetscBool :: energy_flowrate
  PetscBool :: aveg_mass_flowrate
  PetscBool :: aveg_energy_flowrate
  PetscBool :: bool_flag,unsupported_output

  PetscInt :: flag1, flag2

  type(region_type), pointer :: region
  type(flow_condition_type), pointer :: flow_condition
  class(well_data_type), pointer :: well_data
  type(tran_condition_type), pointer :: tran_condition
  class(tran_constraint_base_type), pointer :: tran_constraint
  class(tran_constraint_rt_type), pointer :: sec_tran_constraint
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  type(observation_type), pointer :: observation
  type(integral_flux_type), pointer :: integral_flux

  type(waypoint_type), pointer :: waypoint

  type(material_property_type), pointer :: material_property
  type(fluid_property_type), pointer :: fluid_property
  type(saturation_function_type), pointer :: saturation_function
  class(characteristic_curves_type), pointer :: characteristic_curves
  class(creep_closure_type), pointer :: creep_closure

  class(realization_subsurface_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  class(reaction_rt_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  class(dataset_base_type), pointer :: dataset
  class(dataset_ascii_type), pointer :: dataset_ascii
  class(data_mediator_dataset_type), pointer :: flow_data_mediator
  class(data_mediator_dataset_type), pointer :: rt_data_mediator
  type(waypoint_list_type), pointer :: waypoint_list
  type(input_type), pointer :: input, input_parent

  PetscReal :: dt_init
  PetscReal :: dt_min
  PetscReal :: units_conversion

  class(timestepper_base_type), pointer :: temp_timestepper
  class(timestepper_base_type), pointer :: flow_timestepper
  class(timestepper_BE_type), pointer :: tran_timestepper

  PetscInt::iwaytime,nwaytime,mwaytime
  PetscReal,dimension(:),pointer :: waytime
  PetscReal :: wtime, msfsalt, msfwatr, mlfsalt, mlfwatr

  internal_units = 'not_assigned'

  realization => simulation%realization
  output_option => simulation%output_option
  waypoint_list => simulation%waypoint_list_subsurface
  patch => realization%patch

  if (associated(patch)) grid => patch%grid

  option => realization%option
  field => realization%field
  reaction => realization%reaction

  if ((option%iflowmode == RICHARDS_TS_MODE) .or. &
      (option%iflowmode == TH_TS_MODE)) then
    flow_timestepper => TimestepperTSCreate()
  else
    flow_timestepper => TimestepperBECreate()
  endif
  flow_timestepper%solver%itype = FLOW_CLASS

  tran_timestepper => TimestepperBECreate()
  tran_timestepper%solver%itype = TRANSPORT_CLASS

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
                          
  call InputRewind(input)
  string = 'SUBSURFACE'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call InputPushBlock(input,'SUBSURFACE',option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    option%io_buffer = 'pflotran card:: ' // trim(card)
    call PrintMsg(option)

    select case(trim(card))

!....................
      case ('GRID')
        call DiscretizationRead(realization%discretization,input,option)

!....................
      case ('CHEMISTRY')
        call ReactionReadPass2(reaction,input,option)
        
!....................
      case('NUCLEAR_WASTE_CHEMISTRY')
        if (.not.associated(simulation%nwt_process_model_coupler)) then
          option%io_buffer = 'NUCLEAR_WASTE_CHEMISTRY card is included, but no &
            &NUCLEAR_WASTE_TRANSPORT process model found in SIMULATION block.'
          call PrintErrMsg(option)
        endif
        call NWTReadPass2(realization%reaction_nw,input,option)

!....................
      case ('SPECIFIED_VELOCITY')
        if (option%nflowdof > 0) then
          option%io_buffer = 'SPECIFIED_VELOCITY fields may not be used &
            &with a SUBSURFACE_FLOW mode.'
          call PrintErrMsg(option)
        endif
        internal_units = 'm/sec'
        flag1 = UNINITIALIZED_INTEGER ! uniform?
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','SPECIFIED_VELOCITY')
          call StringToUpper(word)
          select case(trim(word))
            case('UNIFORM?')
              flag1 = StringYesNoOther(input%buf)
            case('DATASET')
              if (flag1 == STRING_OTHER) then
                option%io_buffer = 'SPECIFIED_VELOCITY card "UNIFORM?" &
                  &must be answered with "YES"/"NO" before velocity data &
                  &can can be read.'
                call PrintErrMsg(option)
              endif
              if (flag1 == STRING_YES) then
                error_string = 'SPECIFIED_VELOCITY,UNIFORM,DATASET'
                dataset_ascii => DatasetAsciiCreate()
                dataset_ascii%data_type = DATASET_REAL
                dataset_ascii%array_width = 3 * &
                  max(option%nphase,option%transport%nphase)
                realization%uniform_velocity_dataset => dataset_ascii

                string = input%buf
                call InputReadDouble(input,option,temp_real)
                if (.not.InputError(input)) then
                  error_string = trim(error_string) // ',SINGLE'
                  input%buf = string
                  call DatasetAsciiReadSingle(dataset_ascii,input, &
                                              temp_string,internal_units, &
                                              error_string,option)
                else
                  input%buf = string
                  input%ierr = 0
                  call InputReadCard(input,option,word)
                  call InputErrorMsg(input,option,'keyword',error_string)
                  call StringToUpper(word)
                  select case(word)
                    case('FILE')
                      error_string = trim(error_string) // ',FILE'
                      call InputReadNChars(input,option,string, &
                                           MAXSTRINGLENGTH,PETSC_TRUE)
                      call InputErrorMsg(input,option,'filename',error_string)
                      call DatasetAsciiReadFile(dataset_ascii,string, &
                                                temp_string,internal_units, &
                                                error_string,option)
                    case('LIST')
                      error_string = trim(error_string) // ',LIST'
                      call DatasetAsciiReadList(dataset_ascii,input, &
                                                temp_string,internal_units, &
                                                error_string,option)
                    case default
                      call InputKeywordUnrecognized(input,word, &
                                                    error_string,option)
                  end select
                  if (dataset_ascii%time_storage%time_interpolation_method == &
                      INTERPOLATION_NULL) then
                    option%io_buffer = 'An INTERPOLATION method (LINEAR or &
                      &STEP) must be specified for: ' // trim(error_string)
                    call PrintErrMsg(option)
                  endif
                endif
                bool_flag = PETSC_FALSE
                call DatasetAsciiVerify(dataset_ascii,bool_flag,option)
                if (bool_flag) then
                  option%io_buffer = 'Error verifying ' // &
                    trim(error_string) // '.'
                  call PrintErrMsg(option)
                endif
              else
! Add interface for non-uniform dataset
                call InputReadFilename(input,option, &
                                       realization%nonuniform_velocity_filename)
                call InputErrorMsg(input,option,'filename', &
                                   'SPECIFIED_VELOCITY,NONUNIFORM,DATASET')
              endif
          end select
        enddo
        call InputPopBlock(input,option)
      case ('NONUNIFORM_VELOCITY')
        option%io_buffer = 'The NONUNIFORM_VELOCITY card within SUBSURFACE &
          &block has been deprecated. Use the SPECIFIED_VELOCITY block.'
        call PrintErrMsg(option)
      case ('UNIFORM_VELOCITY')
        option%io_buffer = 'The UNIFORM_VELOCITY card within SUBSURFACE &
          &block has been deprecated. Use the SPECIFIED_VELOCITY block.'
        call PrintErrMsg(option)
      case ('VELOCITY_DATASET')
        option%io_buffer = 'The VELOCITY_DATASET card within SUBSURFACE &
          &block has been deprecated. Use the SPECIFIED_VELOCITY block.'
        call PrintErrMsg(option)

!....................
      case ('DEBUG')
        call DebugRead(realization%debug,input,option)

!....................
      case ('PRINT_PRIMAL_GRID')
        !option%print_explicit_primal_grid = PETSC_TRUE
        option%io_buffer = 'PRINT_PRIMAL_GRID must now be entered under &
                            &OUTPUT card.'
        call PrintErrMsg(option)

!....................
      case ('PRINT_DUAL_GRID')
        !option%print_explicit_dual_grid = PETSC_TRUE
        option%io_buffer = 'PRINT_DUAL_GRID must now be entered under &
                            &OUTPUT card.'
        call PrintErrMsg(option)

!....................
      case ('PROC')

!....................
      case ('REGION')
        region => RegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','REGION')
        call PrintMsg(option,region%name)
        call RegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call RegionAddToList(region,realization%region_list)
        nullify(region)

!....................
      case ('FLOW_CONDITION')
        flow_condition => FlowConditionCreate(option)
        call InputReadWord(input,option,flow_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'FLOW_CONDITION','name')
        call PrintMsg(option,flow_condition%name)
        select case(option%iflowmode)
          case(G_MODE,WF_MODE)
            call FlowConditionGeneralRead(flow_condition,input,option)
          case(H_MODE)
            call FlowConditionHydrateRead(flow_condition,input,option)
          case(TOIL_IMS_MODE)
            call FlowConditionTOilImsRead(flow_condition,input,option)
          case(TOWG_MODE)
            call FlowConditionTOWGRead(flow_condition,input,option)
          case default
            call FlowConditionRead(flow_condition,input,option)
        end select
        call FlowConditionAddToList(flow_condition,realization%flow_conditions)
        nullify(flow_condition)

      case ('WELL_DATA')
        call WellDataSetFlag()
        well_data => WellDataCreate()
        call InputReadWord(input,option,well_data%w_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'WELL_DATA','name')
        call PrintMsg(option,well_data%w_name)
        nwaytime = 0
        mwaytime = 1
        allocate(waytime(mwaytime))
        call well_data%Read(input,option,waytime,nwaytime,mwaytime)
        do iwaytime=1,nwaytime
           wtime=waytime(iwaytime)
           waypoint => WaypointCreate()
           waypoint%time = wtime
           waypoint%update_conditions = PETSC_TRUE
           call WaypointInsertInList(waypoint,waypoint_list)
        enddo
        deallocate(waytime)
        call WellDataAddToList(well_data,realization%well_data)
        nullify(well_data)

!....................
      case ('TRANSPORT_CONDITION')
        if (option%itranmode == NULL_MODE) then
          option%io_buffer = 'TRANSPORT_CONDITIONs not supported without &
                             &CHEMISTRY or SUBSURFACE_NUCLEAR_WASTE_TRANSPORT.'
          call PrintErrMsg(option)
        endif
        tran_condition => TranConditionCreate(option)
        call InputReadWord(input,option,tran_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'TRANSPORT_CONDITION','name')
        call PrintMsg(option,tran_condition%name)
        call TranConditionRead(tran_condition, &
                               realization%transport_constraints, &
                               realization%reaction_base,input,option)
        call TranConditionAddToList(tran_condition, &
                                    realization%transport_conditions)
        nullify(tran_condition)

!....................
      case('CONSTRAINT')
        select case(option%itranmode)
          case(RT_MODE)
            tran_constraint => TranConstraintRTCreate(option)
          case(NWT_MODE)
            tran_constraint => TranConstraintNWTCreate(option)
          case default
            option%io_buffer = 'CONSTRAINTs not supported without CHEMISTRY &
                               &or SUBSURFACE_NUCLEAR_WASTE_TRANSPORT.'
          call PrintErrMsg(option)
        end select
        call InputReadWord(input,option,tran_constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'constraint','name')
        call PrintMsg(option,tran_constraint%name)
        select type(tc=>tran_constraint)
          class is(tran_constraint_rt_type)
            call TranConstraintRTRead(tc,reaction,input,option)
          class is(tran_constraint_nwt_type)
            call TranConstraintNWTRead(tc,realization%reaction_nw,input,option)
        end select
        call TranConstraintAddToList(tran_constraint, &
                                     realization%transport_constraints)
        nullify(tran_constraint)

!....................
      case ('BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Boundary Condition name')
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)

!....................
      case ('INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Initial Condition name')
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)

!....................
      case ('SOURCE_SINK')
        coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name')
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)

!....................
      case ('SOURCE_SINK_SANDBOX')
        call SSSandboxInit(option)
        call SSSandboxRead(input,option)

!....................
      case ('FLOW_MASS_TRANSFER')
        flow_data_mediator => DataMediatorDatasetCreate()
        call InputReadWord(input,option,flow_data_mediator%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Flow Mass Transfer name')
        call DataMediatorDatasetRead(flow_data_mediator,input,option)
        call flow_data_mediator%AddToList(realization%flow_data_mediator_list)
        nullify(flow_data_mediator)

!....................
      case ('RT_MASS_TRANSFER')
        rt_data_mediator => DataMediatorDatasetCreate()
        call InputReadWord(input,option,rt_data_mediator%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'RT Mass Transfer name')
        call DataMediatorDatasetRead(rt_data_mediator,input,option)
        call rt_data_mediator%AddToList(realization%tran_data_mediator_list)
        nullify(rt_data_mediator)

!....................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,input,option)
        call RealizationAddStrata(realization,strata)
        nullify(strata)

!.....................
      case ('DATASET')
        nullify(dataset)
        call DatasetRead(input,dataset,option)
        call DatasetBaseAddToList(dataset,realization%datasets)
        nullify(dataset)

!....................

      case('REFERENCE_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_pressure)
        call InputErrorMsg(input,option,'Reference Pressure','value')
        call InputReadAndConvertUnits(input,option%reference_pressure, &
                                      'Pa','Reference Pressure',option)
!....................

      case('REFERENCE_LIQUID_DENSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option, &
                             option%reference_density(option%liquid_phase))
        call InputErrorMsg(input,option,'Reference Liquid Density','value')
        call InputReadAndConvertUnits(input, &
                              option%reference_density(option%liquid_phase), &
                              'kg/m^3','Reference Density',option)
!....................

      case('REFERENCE_GAS_DENSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option, &
                             option%reference_density(option%gas_phase))
        call InputErrorMsg(input,option,'Reference Gas Density','value')
        call InputReadAndConvertUnits(input, &
                              option%reference_density(option%gas_phase), &
                              'kg/m^3','Reference Density',option)
!....................

      case('MINIMUM_HYDROSTATIC_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%minimum_hydrostatic_pressure)
        call InputErrorMsg(input,option,'Minimum Hydrostatic Pressure','value')
        call InputReadAndConvertUnits(input, &
                                      option%minimum_hydrostatic_pressure, &
                                    'Pa','Minimum Hydrostatic Pressure',option)
!......................

      case('REFERENCE_TEMPERATURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_temperature)
        call InputErrorMsg(input,option,'Reference Temperature','value')

!......................

      case('REFERENCE_POROSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_porosity)
        call InputErrorMsg(input,option,'Reference Porosity','value')

!......................

      case('REFERENCE_SATURATION')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_saturation)
        call InputErrorMsg(input,option,'Reference Saturation','value')

!......................

      case('NONISOTHERMAL')
        option%use_isothermal = PETSC_FALSE

!......................

      case('ISOTHERMAL')
        option%use_isothermal = PETSC_TRUE

!......................

      case('UPDATE_FLOW_PERMEABILITY')
        option%update_flow_perm = PETSC_TRUE

!......................

      case('DFN')
        grid%unstructured_grid%grid_type = TWO_DIM_GRID

!......................

      case("MULTIPLE_CONTINUUM")
        option%io_buffer = 'MULTIPLE_CONTINUUM must be entered under the &
          &SUBSURFACE_TRANSPORT block within the SIMULATION block.'
        call PrintErrMsg(option)

!......................

      case('SECONDARY_CONTINUUM_SOLVER')
        if (.not.option%use_mc) then
          option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can only be used &
                             &with MULTIPLE_CONTINUUM keyword.'
          call PrintErrMsg(option)
        endif
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('KEARST')
            option%secondary_continuum_solver = 1
          case('HINDMARSH')
            option%secondary_continuum_solver = 2
          case('THOMAS')
            option%secondary_continuum_solver = 3
          case default
            option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only &
                               &HINDMARSH or KEARST. For single component &
                               &chemistry THOMAS can be used.'
          call PrintErrMsg(option)
        end select
!....................

      case('SECONDARY_CONSTRAINT')
        if (.not.option%use_mc) then
          option%io_buffer = 'SECONDARY_CONSTRAINT can only be used with &
                             &MULTIPLE_CONTINUUM keyword.'
          call PrintErrMsg(option)
        endif
        if (.not.associated(reaction)) then
          option%io_buffer = 'SECONDARY_CONSTRAINT not supported without &
                             &CHEMISTRY.'
          call PrintErrMsg(option)
        endif
        sec_tran_constraint => TranConstraintRTCreate(option)
        call InputReadWord(input,option,sec_tran_constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'secondary constraint','name')
        call PrintMsg(option,sec_tran_constraint%name)
        call TranConstraintRTRead(sec_tran_constraint,reaction,input,option)
        realization%sec_transport_constraint => sec_tran_constraint
        nullify(sec_tran_constraint)

!......................

      case('BRIN','BRINE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%m_nacl)
        call InputDefaultMsg(input,option,'NaCl Concentration')

        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word(1:len_trim(word)))
          case('MOLAL')
          case('MASS')
            msfsalt = option%m_nacl
            msfwatr = 1.0 - msfsalt
            ! Convert mass salt (kg)/mass water(kg) to Kg-moles = kmol
            ! by dividing mass of salt by molecular weight salt
            ! Then factor of 1000 to convert Kg-mol to mol
            option%m_nacl = (1000.0*msfsalt/FMWNACL)/msfwatr
          case('MOLE')
            ! Convert Kg-mole salt (kg)/Kg mole water(kg) to Kg-mol/Kg water
            ! by multiplying mass of water by molecular weight water
            ! Then factor of 1000 to convert Kg-mol to mol
            mlfsalt = option%m_nacl
            mlfwatr = 1.0 - mlfsalt
            option%m_nacl = 1000.0*mlfsalt/(mlfwatr*FMWH2O)
          case default
            print *, 'Wrong unit: ', word(1:len_trim(word))
            stop
         end select

         ! Report to screen with units
         if (OptionPrintToScreen(option)) then
           print *,'Molality is ',option%m_nacl,' mol/Kg'
         endif

         ! Saturated molality is ~6.16, so above 10 is suspicious
         ! May be units error, so warn user
         if ( option%m_nacl > 10.0 ) then
           option%io_buffer = &
           'More that 10 mols/Kg ~ 584 gms/Kg '// &
           'is an unusually high brine molality'
           call PrintWrnMsg(option)
         endif

!......................

      case ('RESTART')
        option%io_buffer = 'The RESTART card within SUBSURFACE block has &
                           &been deprecated.'
        call PrintErrMsg(option)

!......................

      case ('CHECKPOINT')
        option%io_buffer = 'The CHECKPOINT card within SUBSURFACE block must &
                           &be moved to the SIMULATION block.'
        call PrintErrMsg(option)

!......................

      case ('NUMERICAL_JACOBIAN_FLOW')
        option%io_buffer = 'The NUMERICAL_JACOBIAN_FLOW card within &
          &SUBSURFACE block must be listed under the SIMULATION/&
          &PROCESS_MODELS/SUBSURFACE_FLOW/OPTIONS block as NUMERICAL_JACOBIAN.'
        call PrintErrMsg(option)

!......................

      case ('NUMERICAL_JACOBIAN_RXN')
        option%io_buffer = 'The NUMERICAL_JACOBIAN_RXN card within &
          &SUBSURFACE block must be listed under the SIMULATION/&
          &PROCESS_MODELS/SUBSURFACE_TRANSPORT block as &
          &NUMERICAL_JACOBIAN.'
        call PrintErrMsg(option)

!......................

      case ('NUMERICAL_JACOBIAN_MULTI_COUPLE')
        option%numerical_derivatives_multi_coupling = PETSC_TRUE

!......................

      case ('COMPUTE_STATISTICS')
        option%compute_statistics = PETSC_TRUE

!....................

      case ('CO2_DATABASE')
        call InputReadFilename(input,option,option%co2_database_filename)
        call InputErrorMsg(input,option,'CO2_DATABASE','filename')

!....................

      case ('TIMESTEPPER')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            call flow_timestepper%ReadInput(input,option)
            if (option%flow%resdef) then
              option%flow%flowTimestepperDone = PETSC_TRUE
            endif
          case('TRAN','TRANSPORT')
            call tran_timestepper%ReadInput(input,option)
          case default
            option%io_buffer = 'TIMESTEPPER must specify FLOW or TRANSPORT.'
            call PrintErrMsg(option)
        end select

!....................

      case ('LINEAR_SOLVER')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            call SolverReadLinear(flow_timestepper%solver,input,option)
            if (option%flow%resdef) then
              option%flow%flowSolverLinearDone = PETSC_TRUE
            endif
          case('TRAN','TRANSPORT')
            call SolverReadLinear(tran_timestepper%solver,input,option)
          case default
            option%io_buffer = 'LINEAR_SOLVER must specify FLOW or TRANSPORT.'
            call PrintErrMsg(option)
        end select

!....................

      case ('NEWTON_SOLVER')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            call SolverReadNewton(flow_timestepper%solver,input,option)
          case('TRAN','TRANSPORT')
            call SolverReadNewton(tran_timestepper%solver,input,option)
          case default
            option%io_buffer = 'NEWTON_SOLVER must specify FLOW or TRANSPORT.'
            call PrintErrMsg(option)
        end select
!....................

      case ('FLUID_PROPERTY')

        fluid_property => FluidPropertyCreate()
        call FluidPropertyRead(fluid_property,input,option)
        call FluidPropertyAddToList(fluid_property,realization%fluid_properties)
        nullify(fluid_property)

!....................

      case ('EOS')
        call EOSRead(input,option)

!....................

      case ('SATURATION_FUNCTION')
        if (option%iflowmode == RICHARDS_MODE .or. &
            option%iflowmode == TOIL_IMS_MODE .or. &
            option%iflowmode == TOWG_MODE .or. &
            option%iflowmode == G_MODE .or. &
            option%iflowmode == H_MODE .or. &
            (option%iflowmode == TH_MODE .and. &
             .not. th_use_freezing) .or. &
            option%iflowmode == TH_TS_MODE .or. &
            option%iflowmode == WF_MODE) then
          option%io_buffer = &
            'Must compile with legacy_saturation_function=1 to use the &
            &SATURATION_FUNCTION keyword.  Otherwise, use &
            &CHARACTERISTIC_CURVES.'
          call PrintErrMsg(option)
        endif
        saturation_function => SaturationFunctionCreate(option)
        call InputReadWord(input,option,saturation_function%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','SATURATION_FUNCTION')
        call SaturationFunctionRead(saturation_function,input,option)
        call SatFunctionComputePolynomial(option,saturation_function)
        call PermFunctionComputePolynomial(option,saturation_function)
        call SaturationFunctionAddToList(saturation_function, &
                                         realization%saturation_functions)
        nullify(saturation_function)

!....................

      case ('CHARACTERISTIC_CURVES')
        if (.not.(option%iflowmode == NULL_MODE .or. &
                  option%iflowmode == RICHARDS_MODE .or. &
                  option%iflowmode == RICHARDS_TS_MODE .or. &
                  option%iflowmode == TOIL_IMS_MODE .or. &
                  option%iflowmode == TOWG_MODE .or. &
                  option%iflowmode == G_MODE .or. &
                  option%iflowmode == H_MODE .or. &
                  option%iflowmode == TH_TS_MODE .or. &
                  (option%iflowmode == TH_MODE .and. &
                    .not. th_use_freezing) .or. &
                  option%iflowmode == WF_MODE)) then
          option%io_buffer = 'CHARACTERISTIC_CURVES not supported in flow &
            &modes other than RICHARDS, RICHARDS_TS, TOIL_IMS, WIPP_FLOW, TH, or GENERAL. &
            &Use SATURATION_FUNCTION.'
          call PrintErrMsg(option)
        endif
        characteristic_curves => CharacteristicCurvesCreate()
        call InputReadWord(input,option,characteristic_curves%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','CHARACTERISTIC_CURVES')
        option%io_buffer = '  Name :: ' // &
          trim(characteristic_curves%name)
        call PrintMsg(option)
        call CharacteristicCurvesRead(characteristic_curves,input,option)
!        call SatFunctionComputePolynomial(option,saturation_function)
!        call PermFunctionComputePolynomial(option,saturation_function)
        call CharacteristicCurvesAddToList(characteristic_curves, &
                                          realization%characteristic_curves)
        nullify(characteristic_curves)

!....................
      case('CREEP_CLOSURE_TABLE')
        wipp => WIPPGetPtr()
        option%flow%transient_porosity = PETSC_TRUE
        option%flow%creep_closure_on = PETSC_TRUE
        creep_closure => CreepClosureCreate()
        call InputReadWord(input,option,creep_closure%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','CREEP_CLOSURE_TABLE')
        option%io_buffer = '  Name :: ' // trim(creep_closure%name)
        call PrintMsg(option)
        call creep_closure%Read(input,option)
        call CreepClosureAddToList(creep_closure, &
             wipp%creep_closure_tables)
        nullify(creep_closure)

!....................

      case ('MATERIAL_PROPERTY')

        material_property => MaterialPropertyCreate()
        call InputReadWord(input,option,material_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')
        option%io_buffer = '  Name :: ' // trim(material_property%name)
        call PrintMsg(option)
        call MaterialPropertyRead(material_property,input,option)
        call MaterialPropertyAddToList(material_property, &
             realization%material_properties)
        nullify(material_property)

!....................

      case ('USE_TOUCH_OPTIONS')
        option%use_touch_options = PETSC_TRUE

      case ('MPI_IO')
!        call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
!                                       '-viewer_binary_mpiio')

      case ('HANDSHAKE_IO')
        call InputReadInt(input,option,option%io_handshake_buffer_size)
        call InputErrorMsg(input,option,'io_handshake_buffer_size', &
                           'HANDSHAKE_IO')

      case ('OVERWRITE_RESTART_TRANSPORT')
        option%io_buffer = 'OVERWRITE_RESTART_TRANSPORT no longer &
          &supported. Please use SKIP_RESTART in the SUBSURFACE_TRANSPORT &
          &process model options block.'
        call PrintErrMsg(option)

      case ('OVERWRITE_RESTART_FLOW_PARAMS')
        option%io_buffer = 'OVERWRITE_RESTART_FLOW_PARAMS no longer &
          &supported. Please use REVERT_PARAMETERS_ON_RESTART in &
          &the SUBSURFACE_FLOW process model options block.'
        call PrintErrMsg(option)

      case ('INITIALIZE_FLOW_FROM_FILE')
        call InputReadFilename(input,option,option%initialize_flow_filename)
        call InputErrorMsg(input,option,'filename','INITIALIZE_FLOW_FROM_FILE')

      case ('INITIALIZE_TRANSPORT_FROM_FILE')
        call InputReadFilename(input,option, &
                               option%initialize_transport_filename)
        call InputErrorMsg(input,option,'filename', &
                           'INITIALIZE_TRANSPORT_FROM_FILE')

      case ('CENTRAL_DIFFERENCE')
        option%use_upwinding = PETSC_FALSE

!....................
      case ('OBSERVATION')
        observation => ObservationCreate()
        call ObservationRead(observation,input,option)
        call ObservationAddToList(observation, &
                                  realization%patch%observation_list)
        nullify(observation)

!....................
      case ('INTEGRAL_FLUX')
        integral_flux => IntegralFluxCreate()
        call InputReadWord(input,option,integral_flux%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Integral Flux name')
        call IntegralFluxRead(integral_flux,input,option)
        call IntegralFluxAddToList(integral_flux, &
                                   realization%patch%integral_flux_list)
        nullify(integral_flux)

!.....................
      case ('WALLCLOCK_STOP')
        option%wallclock_stop_flag = PETSC_TRUE
        call InputReadDouble(input,option,option%wallclock_stop_time)
        call InputErrorMsg(input,option,'stop time','WALLCLOCK_STOP')

        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr /= 0) word = 'h'
        call InputDefaultMsg(input,option,'WALLCLOCK_STOP time units')
        internal_units = 'sec'
        units_conversion = UnitsConvertToInternal(word,internal_units,option)
        ! convert from hrs to seconds and add to start_time
        option%wallclock_stop_time = option%start_time + &
                                     option%wallclock_stop_time* &
                                     units_conversion

!....................
      case ('OUTPUT')
        vel_cent = PETSC_FALSE
        vel_face = PETSC_FALSE
        fluxes = PETSC_FALSE
        mass_flowrate = PETSC_FALSE
        energy_flowrate = PETSC_FALSE
        aveg_mass_flowrate = PETSC_FALSE
        aveg_energy_flowrate = PETSC_FALSE
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','OUTPUT')
          call StringToUpper(word)
        !----------------------------------------------------------------------
        !----- NEW INPUT FORMAT: ----------------------------------------------
        !----------------------------------------------------------------------
          select case(trim(word))
            case('OBSERVATION_FILE')
              call OutputFileRead(input,realization,output_option, &
                                  waypoint_list,trim(word))
            case('SNAPSHOT_FILE')
              call OutputFileRead(input,realization,output_option, &
                                  waypoint_list,trim(word))
            case('MASS_BALANCE_FILE')
              call OutputFileRead(input,realization,output_option, &
                                  waypoint_list,trim(word))
            case('ECLIPSE_FILE')
              call CreateOutputOptionEclipse(output_option)
              call OutputFileRead(input,realization,output_option, &
                                  waypoint_list,trim(word))
            case('LINEREPT')
               option%linerept = PETSC_TRUE
            case('TIME_UNITS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Output Time Units','OUTPUT')
              output_option%tunit = trim(word)
              internal_units = 'sec'
              output_option%tconv = &
                UnitsConvertToInternal(word,internal_units,option)
            case('VARIABLES')
              call OutputVariableRead(input,option, &
                                      output_option%output_variable_list)
            case('AVERAGE_VARIABLES')
              call OutputVariableRead(input,option, &
                                      output_option%aveg_output_variable_list)
            case('UNFILTER_NON_STATE_VARIABLES')
              output_option%filter_non_state_variables = PETSC_FALSE


        !----------------------------------------------------------------------
        !----- SUPPORT FOR OLD INPUT FORMAT: ----------------------------------
        !----------------------------------------------------------------------
            case('NO_FINAL','NO_PRINT_FINAL')
              output_option%print_final_obs = PETSC_FALSE
              output_option%print_final_snap = PETSC_FALSE
              output_option%print_final_massbal = PETSC_FALSE
            case('NO_INITIAL','NO_PRINT_INITIAL')
              output_option%print_initial_obs = PETSC_FALSE
              output_option%print_initial_snap = PETSC_FALSE
              output_option%print_initial_massbal = PETSC_FALSE
            case('PROCESSOR_ID')
              option%io_buffer = 'PROCESSOR_ID output must now be entered &
                                 &under OUTPUT/VARIABLES card as PROCESS_ID.'
              call PrintErrMsg(option)
!              output_option%print_iproc = PETSC_TRUE
            case('PERMEABILITY')
              option%io_buffer = 'PERMEABILITY output must now be entered &
                                 &under OUTPUT/VARIABLES card.'
              call PrintErrMsg(option)
!              output_option%print_permeability = PETSC_TRUE
            case('POROSITY')
              option%io_buffer = 'POROSITY output must now be entered under &
                                 &OUTPUT/VARIABLES card.'
              call PrintErrMsg(option)
!              output_option%print_porosity = PETSC_TRUE
            case('TORTUOSITY')
              option%io_buffer = 'TORTUOSITY output must now be entered under &
                                 &OUTPUT/VARIABLES card.'
              call PrintErrMsg(option)
!              output_option%print_tortuosity = PETSC_TRUE
            case('VOLUME')
              option%io_buffer = 'VOLUME output must now be entered under &
                                 &OUTPUT/VARIABLES card.'
              call PrintErrMsg(option)
!              output_option%print_volume = PETSC_TRUE
            case('MASS_BALANCE')
              option%compute_mass_balance_new = PETSC_TRUE
              output_option%periodic_msbl_output_ts_imod = 1
              call InputReadCard(input,option,word)
              call InputDefaultMsg(input,option, &
                                   'OUTPUT,MASS_BALANCE,DETAILED')
              if (len_trim(word) > 0) then
                call StringToUpper(word)
                select case(trim(word))
                  case('DETAILED')
                    option%mass_bal_detailed = PETSC_TRUE
                  case default
                    call InputKeywordUnrecognized(input,word, &
                           'OUTPUT,MASS_BALANCE',option)
                end select
              endif
            case('PRINT_COLUMN_IDS')
              output_option%print_column_ids = PETSC_TRUE

           case ('PRINT_PRIMAL_GRID')
             output_option%print_explicit_primal_grid = PETSC_TRUE

           !out_mesh_type defaults for primal_explicit grid is vetex_centered
           case ('EXPLICIT_GRID_PRIMAL_GRID_TYPE')
             if (associated(grid%unstructured_grid)) then
               if (associated(grid%unstructured_grid%explicit_grid)) then
                 call InputReadCard(input,option,word)
                 call InputErrorMsg(input,option,word, &
                       'EXPLICIT_GRID_PRIMAL_GRID_TYPE')
                 call PrintMsg(option,word)
                 call StringToUpper(word)
                   select case (trim(word))
                     case ('VERTEX_CENTERED')
                       grid%unstructured_grid%explicit_grid% &
                          output_mesh_type = VERTEX_CENTERED_OUTPUT_MESH
                     case ('CELL_CENTERED')
                       grid%unstructured_grid%explicit_grid% &
                          output_mesh_type = CELL_CENTERED_OUTPUT_MESH
                       call OptionSetBlocking(option,PETSC_FALSE)
                       if ( option%myrank == option%io_rank ) then
                         if (grid%unstructured_grid% &
                             explicit_grid%num_elems /= &
                             grid%unstructured_grid% &
                             explicit_grid%num_cells_global &
                            ) then
                           option%io_buffer = &
                             'EXPLICIT_GRID_PRIMAL_GRID_TYPE &
                             &if CELL_CENTERED option, the number of cells &
                             &of the grid to print and those &
                             &of the computational grid must be equal.'
                           call PrintErrMsg(option)
                         end if
                       end if
                       call OptionSetBlocking(option,PETSC_TRUE)
                       call OptionCheckNonBlockingError(option)
                     case default
                       option%io_buffer ='EXPLICIT_GRID_PRIMAL_GRID_TYPE &
                                  &only VERTEX_CENTERED and CELL_CENTERED &
                                  &are supported.'
                       call PrintErrMsg(option)
                   end select
               endif
             endif

           case ('PRINT_DUAL_GRID')
             output_option%print_explicit_dual_grid = PETSC_TRUE

            case('TIMES')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'units','OUTPUT,TIMES')
              internal_units = 'sec'
              units_conversion = &
                UnitsConvertToInternal(word,internal_units,option)
              string = 'OUTPUT,TIMES'
              nullify(temp_real_array)
              call UtilityReadArray(temp_real_array,NEG_ONE_INTEGER, &
                                    string,input,option)
              do temp_int = 1, size(temp_real_array)
                waypoint => WaypointCreate()
                waypoint%time = temp_real_array(temp_int)*units_conversion
                waypoint%print_snap_output = PETSC_TRUE
                call WaypointInsertInList(waypoint,waypoint_list)
              enddo
              call DeallocateArray(temp_real_array)
            case('OUTPUT_FILE')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'time increment', &
                                 'OUTPUT,OUTPUT_FILE')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_file = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%output_file_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,OUTPUT_FILE,PERIODIC')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'OUTPUT,OUTPUT_FILE',option)
              end select
            case('SCREEN')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'time increment','OUTPUT,SCREEN')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_screen = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%screen_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,SCREEN')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'OUTPUT,SCREEN',option)
              end select
            case('PERIODIC')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'time increment', &
                                 'OUTPUT,PERIODIC')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  internal_units = 'sec'
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'OUTPUT,PERIODIC,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'OUTPUT,PERIODIC,TIME')
                  units_conversion = UnitsConvertToInternal(word, &
                                     internal_units,option)
                  output_option%periodic_snap_output_time_incr = temp_real* &
                                                            units_conversion
                  call InputReadCard(input,option,word)
                  if (input%ierr == 0) then
                    if (StringCompareIgnoreCase(word,'between')) then
                      call InputReadDouble(input,option,temp_real)
                      call InputErrorMsg(input,option,'start time', &
                                         'OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'start time units', &
                                         'OUTPUT,PERIODIC,TIME')
                      units_conversion = UnitsConvertToInternal(word, &
                                         internal_units,option)
                      temp_real = temp_real * units_conversion
                      call InputReadCard(input,option,word)
                      if (.not.StringCompareIgnoreCase(word,'and')) then
                        input%ierr = 1
                      endif
                      call InputErrorMsg(input,option,'and', &
                                          'OUTPUT,PERIODIC,TIME"')
                      call InputReadDouble(input,option,temp_real2)
                      call InputErrorMsg(input,option,'end time', &
                                         'OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'end time units', &
                                         'OUTPUT,PERIODIC,TIME')
                      units_conversion = UnitsConvertToInternal(word, &
                                         internal_units,option)
                      temp_real2 = temp_real2 * units_conversion
                      do
                        waypoint => WaypointCreate()
                        waypoint%time = temp_real
                        waypoint%print_snap_output = PETSC_TRUE
                        call WaypointInsertInList(waypoint,waypoint_list)
                        temp_real = temp_real + &
                          output_option%periodic_snap_output_time_incr
                        if (temp_real > temp_real2) exit
                      enddo
                      output_option%periodic_snap_output_time_incr = 0.d0
                    else
                      input%ierr = 1
                      call InputErrorMsg(input,option,'between', &
                                          'OUTPUT,PERIODIC,TIME')
                    endif
                  endif
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_snap_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,TIMESTEP')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'OUTPUT,PERIODIC',option)
              end select
            case('OBSERVATION_TIMES')
              output_option%print_observation = PETSC_TRUE
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time units', &
                   'OUTPUT,OBSERVATION_TIMES')
              internal_units = 'sec'
              units_conversion = &
                UnitsConvertToInternal(word,internal_units,option)
              string = 'OBSERVATION_TIMES,TIMES'
              nullify(temp_real_array)
              call UtilityReadArray(temp_real_array,NEG_ONE_INTEGER, &
                                    string,input,option)
              do temp_int = 1, size(temp_real_array)
                waypoint => WaypointCreate()
                waypoint%time = temp_real_array(temp_int)*units_conversion
                waypoint%print_obs_output = PETSC_TRUE
                call WaypointInsertInList(waypoint,waypoint_list)
                waypoint => WaypointCreate()
                waypoint%time = temp_real_array(temp_int)*units_conversion
                waypoint%print_msbl_output = PETSC_TRUE
                call WaypointInsertInList(waypoint,waypoint_list)
              enddo
              call DeallocateArray(temp_real_array)
            case('PERIODIC_OBSERVATION')
              output_option%print_observation = PETSC_TRUE
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'time increment', &
                'OUTPUT, PERIODIC_OBSERVATION')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIME')
                  internal_units = 'sec'
                  units_conversion = UnitsConvertToInternal(word, &
                                     internal_units,option)
                  output_option%periodic_obs_output_time_incr = temp_real* &
                                                               units_conversion
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_obs_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIMESTEP')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'OUTPUT,PERIODIC_OBSERVATION',option)
              end select
            case('FORMAT')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'keyword','OUTPUT,FORMAT')
              call StringToUpper(word)
              select case(trim(word))
                case ('HDF5')
                  output_option%print_hdf5 = PETSC_TRUE
                  call InputReadCard(input,option,word)
                  call InputDefaultMsg(input,option, &
                                       'OUTPUT,FORMAT,HDF5,# FILES')
                  if (len_trim(word) > 0) then
                    call StringToUpper(word)
                    select case(trim(word))
                      case('SINGLE_FILE')
                        output_option%print_single_h5_file = PETSC_TRUE
                      case('MULTIPLE_FILES')
                        output_option%print_single_h5_file = PETSC_FALSE
                        output_option%times_per_h5_file = 1
                        call InputReadCard(input,option,word)
                        if (len_trim(word)>0) then
                          select case(trim(word))
                            case('TIMES_PER_FILE')
                              call InputReadInt(input,option, &
                                              output_option%times_per_h5_file)
                              call InputErrorMsg(input,option, &
                                'timestep increment', &
                                'OUTPUT,FORMAT,HDF5,MULTIPLE_FILES,&
                                &TIMES_PER_FILE')
                            case default
                              call InputKeywordUnrecognized(input,word, &
                                    'OUTPUT,FORMAT,HDF5,MULTIPLE_FILES',option)
                          end select
                        endif
                      case default
                        call InputKeywordUnrecognized(input,word, &
                               'OUTPUT,FORMAT,HDF5',option)
                    end select
                  endif
                case ('TECPLOT')
                  output_option%print_tecplot = PETSC_TRUE
                  call InputReadCard(input,option,word)
                  call InputErrorMsg(input,option,'TECPLOT','OUTPUT,FORMAT')
                  call StringToUpper(word)
                  select case(trim(word))
                    case('POINT')
                      output_option%tecplot_format = TECPLOT_POINT_FORMAT
                    case('BLOCK')
                      output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
                    case('FEBRICK')
                      output_option%tecplot_format = TECPLOT_FEBRICK_FORMAT
                    case default
                      call InputKeywordUnrecognized(input,word, &
                               'OUTPUT,FORMAT,TECPLOT',option)
                  end select
                  if (output_option%tecplot_format == TECPLOT_POINT_FORMAT &
                      .and. option%mycommsize > 1) then
                    output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
                  endif
                  if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
                    output_option%tecplot_format = TECPLOT_FEBRICK_FORMAT
                  endif
                case ('VTK')
                  output_option%print_vtk = PETSC_TRUE
                case default
                  call InputKeywordUnrecognized(input,word, &
                                                'OUTPUT,FORMAT',option)
              end select
            case('VELOCITY_AT_CENTER')
              vel_cent = PETSC_TRUE
            case('VELOCITY_AT_FACE')
              vel_face = PETSC_TRUE
            case('FLUXES')
              fluxes = PETSC_TRUE
            case('FLOWRATES','FLOWRATE')
              mass_flowrate = PETSC_TRUE
              energy_flowrate = PETSC_TRUE
            case('MASS_FLOWRATE')
              mass_flowrate = PETSC_TRUE
            case('ENERGY_FLOWRATE')
              energy_flowrate = PETSC_TRUE
            case('AVERAGE_FLOWRATES','AVERAGE_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
              aveg_energy_flowrate = PETSC_TRUE
            case('AVERAGE_MASS_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
            case('AVERAGE_ENERGY_FLOWRATE')
              aveg_energy_flowrate = PETSC_TRUE
            case ('HDF5_WRITE_GROUP_SIZE')
              call InputReadInt(input,option,option%hdf5_write_group_size)
              call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE', &
                                 'Group size')
            case('EXTEND_HDF5_TIME_FORMAT')
              output_option%extend_hdf5_time_format = PETSC_TRUE
            case default
              call InputKeywordUnrecognized(input,word,'OUTPUT',option)
          end select

        enddo
        call InputPopBlock(input,option)

  ! If VARIABLES were not specified within the *_FILE blocks, point their
  ! variable lists to the master variable list, which can be specified within
  ! the OUTPUT block. If no VARIABLES are specified for the master list, the
  ! defaults will be populated.
          if (.not.associated(output_option%output_snap_variable_list%first) &
              .and.(output_option%output_snap_variable_list%flow_vars .and. &
                    output_option%output_snap_variable_list%energy_vars)) then
            call OutputVariableListDestroy( &
                 output_option%output_snap_variable_list)
            output_option%output_snap_variable_list => &
                 output_option%output_variable_list
          endif
          if (.not.associated(output_option%output_obs_variable_list%first) &
              .and.(output_option%output_obs_variable_list%flow_vars .and. &
                    output_option%output_obs_variable_list%energy_vars)) then
            call OutputVariableListDestroy( &
                 output_option%output_obs_variable_list)
            output_option%output_obs_variable_list => &
                output_option%output_variable_list
          endif

        if (vel_cent) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_vel_cent = PETSC_TRUE
          if (output_option%print_hdf5) &
            output_option%print_hdf5_vel_cent = PETSC_TRUE
          if (output_option%print_vtk) &
            output_option%print_vtk_vel_cent = PETSC_TRUE
        endif
        if (vel_face) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_vel_face = PETSC_TRUE
          if (output_option%print_hdf5) &
           output_option%print_hdf5_vel_face = PETSC_TRUE
        endif
        if (fluxes) then
          output_option%print_fluxes = PETSC_TRUE
        endif
        if(output_option%aveg_output_variable_list%nvars>0) then
          if(Equal(output_option%periodic_snap_output_time_incr,0.d0)) then
            option%io_buffer = 'Keyword: AVERAGE_VARIABLES defined without &
                               &PERIODIC TIME being set.'
            call PrintErrMsg(option)
          endif
          if(.not.output_option%print_hdf5) then
            option%io_buffer = 'Keyword: AVERAGE_VARIABLES only defined for FORMAT HDF5'
            call PrintErrMsg(option)
          endif
        endif
        if (mass_flowrate .or. energy_flowrate .or. aveg_mass_flowrate .or. &
            aveg_energy_flowrate .or. fluxes) then
          if (output_option%print_hdf5) then
            output_option%print_hdf5_mass_flowrate = mass_flowrate
            output_option%print_hdf5_energy_flowrate = energy_flowrate
            output_option%print_hdf5_aveg_mass_flowrate = aveg_mass_flowrate
            output_option%print_hdf5_aveg_energy_flowrate = aveg_energy_flowrate
            if(aveg_mass_flowrate.or.aveg_energy_flowrate) then
              if(Equal(output_option%periodic_snap_output_time_incr,0.d0)) then
                option%io_buffer = 'Keyword: AVEGRAGE_FLOWRATES/ &
                  &AVEGRAGE_MASS_FLOWRATE/ENERGY_FLOWRATE defined without &
                  &PERIODIC TIME being set.'
                call PrintErrMsg(option)
              endif
            endif
          endif
          option%flow%store_fluxes = PETSC_TRUE
          if (associated(grid%unstructured_grid)) then
            if (associated(grid%unstructured_grid%explicit_grid)) then
              option%flow%store_fluxes = PETSC_TRUE
              output_option%print_explicit_flowrate = mass_flowrate
            endif
          endif
        endif
        if (associated(grid%unstructured_grid)) then
          if (associated(grid%unstructured_grid%explicit_grid)) then
            if( output_option%write_ecl .or. option%linerept ) then
              unsupported_output = output_option%print_tecplot &
                                    .or. output_option%print_vtk
            else
              unsupported_output = .not.output_option%print_hdf5
            endif

            if ( unsupported_output .and.  &
                 (grid%unstructured_grid%explicit_grid%output_mesh_type == &
                 CELL_CENTERED_OUTPUT_MESH) &
               ) then
                option%io_buffer = 'unstructured explicit grid &
                  &output_mesh_type = CELL_CENTERED supported for hdf5 only'
                call PrintErrMsg(option)
            end if
          end if
        end if

!.....................
      case ('REGRESSION')
        call RegressionRead(simulation%regression,input,option)

!.....................
      case ('TIME')
!        dt_init = UNINITIALIZED_DOUBLE
        dt_init = 1.d0
        dt_min = UNINITIALIZED_DOUBLE
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'word','TIME')
          select case(trim(word))
            case('SCREEN_UNITS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Screen Units','TIME')
              internal_units = 'sec'
              temp_real2 = UnitsConvertToInternal(word,internal_units,option)
              output_option%tunit = trim(word)
              output_option%tconv = temp_real2
            case('STEADY_STATE')
              option%steady_state = PETSC_TRUE
            case('FINAL_TIME')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Final Time','TIME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Final Time Units','TIME')
              internal_units = 'sec'
              temp_real2 = UnitsConvertToInternal(word,internal_units,option)
              if (len_trim(output_option%tunit) == 0) then
                output_option%tunit = trim(word)
                output_option%tconv = temp_real2
              endif
              waypoint => WaypointCreate()
              waypoint%final = PETSC_TRUE
              waypoint%time = temp_real*temp_real2
              waypoint%print_snap_output = PETSC_TRUE
              call WaypointInsertInList(waypoint,waypoint_list)
            case('INITIAL_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Initial Timestep Size','TIME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Initial Timestep Size Time &
                                              &Units','TIME')
              internal_units = 'sec'
              dt_init = temp_real*UnitsConvertToInternal(word, &
                                                         internal_units,option)
            case('MINIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Minimum Timestep Size','TIME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Minimum Timestep Size Time &
                                              &Units','TIME')
              internal_units = 'sec'
              dt_min = temp_real*UnitsConvertToInternal(word, &
                                                        internal_units,option)
            case('MAXIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Maximum Timestep Size','TIME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Maximum Timestep Size Time &
                                              &Units','TIME')
              waypoint => WaypointCreate()
              internal_units = 'sec'
              waypoint%dt_max = temp_real*UnitsConvertToInternal(word, &
                                          internal_units,option)
              call InputReadCard(input,option,word)
              if (input%ierr == 0) then
                call StringToUpper(word)
                if (StringCompare(word,'AT',TWO_INTEGER)) then
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'Maximum Timestep Size &
                                                  &Update Time','TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'Maximum Timestep Size &
                                                  &Update Time Units','TIME')
                  internal_units = 'sec'
                  waypoint%time = temp_real*UnitsConvertToInternal(word, &
                                            internal_units,option)
                else
                  option%io_buffer = 'Keyword under "MAXIMUM_TIMESTEP_SIZE" &
                                     &after maximum timestep size should &
                                     &be "at".'
                  call PrintErrMsg(option)
                endif
              else
                waypoint%time = 0.d0
              endif
              call WaypointInsertInList(waypoint,waypoint_list)
            case default
              call InputKeywordUnrecognized(input,word,'TIME',option)
          end select
        enddo
        call InputPopBlock(input,option)
        if (Initialized(dt_init)) then
          if (associated(flow_timestepper)) then
            flow_timestepper%dt_init = dt_init
            option%flow_dt = dt_init
          endif
          if (associated(tran_timestepper)) then
            tran_timestepper%dt_init = dt_init
            option%tran_dt = dt_init
          endif
        endif
        if (Initialized(dt_min)) then
          if (associated(flow_timestepper)) then
            flow_timestepper%dt_min = dt_min
          endif
          if (associated(tran_timestepper)) then
            tran_timestepper%dt_min = dt_min
          endif
        endif

!......................
      case ('HDF5_READ_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_read_group_size)
        call InputErrorMsg(input,option,'HDF5_READ_GROUP_SIZE','Group size')

!......................
      case ('HDF5_WRITE_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_write_group_size)
        call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')

!....................
      case('WIPP')
        wipp => WIPPGetPtr()
        call WIPPRead(input,option)

!....................
      case('KLINKENBERG_EFFECT')
        wipp => WIPPGetPtr()
        call KlinkenbergInit()
        klinkenberg => KlinkenbergCreate()
        call Klinkenberg%Read(input,option)

!....................
      case ('ONLY_VERTICAL_FLOW')
        option%flow%only_vertical_flow = PETSC_TRUE
        if (option%iflowmode /= TH_MODE .and. &
            option%iflowmode /= TH_TS_MODE .and. &
            option%iflowmode /= RICHARDS_MODE .and. &
            option%iflowmode /= RICHARDS_TS_MODE) then
          option%io_buffer = 'ONLY_VERTICAL_FLOW implemented in RICHARDS, &
                              &RICHARDS_TS and TH, TH_TS modes.'
          call PrintErrMsg(option)
        endif

!....................
      case ('QUASI_3D')
        option%flow%quasi_3d = PETSC_TRUE
        option%flow%only_vertical_flow = PETSC_TRUE
        if (option%iflowmode /= RICHARDS_MODE .and. &
            option%iflowmode /= RICHARDS_TS_MODE) then
          option%io_buffer = 'QUASI_3D implemented in RICHARDS and &
                              &RICHARDS_TS modes.'
          call PrintErrMsg(option)
        endif

!....................
      case ('ONLY_ENERGY_EQ')
        option%flow%only_energy_eq = PETSC_TRUE
        if (option%iflowmode /= TH_MODE .and. &
            option%iflowmode /= TH_TS_MODE) then
          option%io_buffer = 'ONLY_ENERGY_EQ applicable only in TH and &
                              &TH_TS modes.'
          call PrintErrMsg(option)
        endif

!....................
      case ('RELATIVE_PERMEABILITY_AVERAGE')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (trim(word))
          case ('UPWIND')
            option%rel_perm_aveg = UPWIND
          case ('HARMONIC')
            option%rel_perm_aveg = HARMONIC
          case ('DYNAMIC_HARMONIC')
            option%rel_perm_aveg = DYNAMIC_HARMONIC
          case default
            option%io_buffer = 'Cannot identify the specificed &
              &RELATIVE_PERMEABILITY_AVERAGE.'
            call PrintErrMsg(option)
          end select

!....................
      case ('DBASE_FILENAME')

!....................
      case ('END_SUBSURFACE')

        if (option%flow%resdef) then
          ! it is possible for a card to be skipped when we would like to set defaults for it,
          ! so we check that here
          if (.NOT. option%flow%flowSolverLinearDone) then
              call SolverReadLinear(flow_timestepper%solver,input,option)
              option%flow%flowSolverLinearDone = PETSC_TRUE
          endif
          if (.NOT. option%flow%flowTimestepperDone) then
            call flow_timestepper%ReadInput(input,option)
            option%flow%flowTimestepperDone = PETSC_TRUE
          endif
        endif

        exit

!....................
      case ('MIN_ALLOWABLE_SCALE')
        call InputReadDouble(input,option,option%min_allowable_scale)
        call InputErrorMsg(input,option,'minimium allowable scaling factor', &
                           'InitSubsurface')

!....................
      case ('HYDRATE')
        call HydrateRead(input,patch%methanogenesis,option)
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'SubsurfaceReadInput()',option)
    end select

  enddo
  call InputPopBlock(input,option) ! SUBSURFACE

  if (associated(simulation%flow_process_model_coupler)) then
    if (option%iflowmode == RICHARDS_TS_MODE) then
       if (option%steady_state) then
         option%io_buffer = 'Steady state not supported with PETSC TS'
         call PrintErrMsg(option)
       endif
    endif
    flow_timestepper%name = 'FLOW'
    if (option%steady_state) then
      !geh: This is a workaround for the Intel compiler which thinks that
      !     fts and/or flow_timestepper is not a pointer when passed within
      !     the select type statment.  This is too convoluted; it needs to
      !     be refactored!
      select type(fts=>flow_timestepper)
        class is(timestepper_BE_type)
          temp_timestepper => TimestepperSteadyCreateFromBE(fts)
        class is(timestepper_base_type)
          temp_timestepper =>  TimestepperSteadyCreateFromBase(fts)
      end select
      call flow_timestepper%Destroy()
      deallocate(flow_timestepper)
      nullify(flow_timestepper)
      flow_timestepper => temp_timestepper
      nullify(temp_timestepper)
    endif
    simulation%flow_process_model_coupler%timestepper => flow_timestepper
  else
    call flow_timestepper%Destroy()
    deallocate(flow_timestepper)
    nullify(flow_timestepper)
  endif
  if (associated(simulation%rt_process_model_coupler) .or. &
      associated(simulation%nwt_process_model_coupler)) then
    tran_timestepper%name = 'TRAN'
    if (option%steady_state) then
      ! transport is currently always BE
      temp_timestepper => TimestepperSteadyCreateFromBE(tran_timestepper)
      call tran_timestepper%Destroy()
      deallocate(tran_timestepper)
      nullify(tran_timestepper)
      flow_timestepper => temp_timestepper
      nullify(temp_timestepper)
    endif
    if (associated(simulation%rt_process_model_coupler)) then
      simulation%rt_process_model_coupler%timestepper => tran_timestepper
    endif
    if (associated(simulation%nwt_process_model_coupler)) then
      simulation%nwt_process_model_coupler%timestepper => tran_timestepper
    endif
  else
    call tran_timestepper%Destroy()
    deallocate(tran_timestepper)
    nullify(tran_timestepper)
  endif

end subroutine SubsurfaceReadInput

end module Factory_Subsurface_module
