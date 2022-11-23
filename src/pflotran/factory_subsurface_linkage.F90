module Factory_Subsurface_Linkage_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Simulation_Subsurface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: &
            SetupPMApproach, &
            ExtractPMsFromPMList, &
            SetupPMCLinkages, &
            FactorySubsurfSetPMCWaypointPtrs

contains

! ************************************************************************** !

recursive subroutine SetupPMApproach(pmc,simulation)
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
  use PMC_Geophysics_class
  use PM_Base_Pointer_module
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_RT_class
  use PM_NWT_class
  use PM_Waste_Form_class
  use PM_WIPP_SrcSink_class
  use PM_UFD_Decay_class
  use PM_UFD_Biosphere_class
  use PM_ERT_class
  use PM_Well_class
  use PM_Material_Transform_class
  use Option_module
  use Simulation_Subsurface_class
  use Realization_Subsurface_class

  implicit none


  class(pmc_base_type), pointer :: pmc
  class(simulation_subsurface_type) :: simulation

  class(realization_subsurface_type), pointer :: realization
  class(pm_base_type), pointer :: cur_pm
  type(option_type), pointer :: option

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
          option%io_buffer = 'SUBSURFACE_TRANSPORT MODE GIRT/OSRT is &
            &specified in the SIMULATION block without the corresponding &
            &process model without a corresponding CHEMISTRY block within &
            &the SUBSURFACE block.'
          call PrintErrMsg(option)
        endif
        call cur_pm%SetRealization(realization)

      class is(pm_nwt_type)
        if (.not.associated(realization%reaction_nw)) then
          option%io_buffer = 'SUBSURFACE_TRANSPORT MODE NWT is specified &
            &in the SIMULATION block without the corresponding &
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

      class is(pm_ert_type)
        call cur_pm%SetRealization(realization)

      class is(pm_well_type)
        call cur_pm%SetRealization(realization)

      class is(pm_material_transform_type)
        call cur_pm%SetRealization(realization)

    end select

    cur_pm%output_option => simulation%output_option
    call cur_pm%Setup()
    cur_pm => cur_pm%next
  enddo

  call pmc%SetupSolvers()

  ! call this function for this pmc's child
  if (associated(pmc%child)) then
    call SetupPMApproach(pmc%child,simulation)
  endif

  ! call this function for this pmc's peer
  if (associated(pmc%peer)) then
    call SetupPMApproach(pmc%peer,simulation)
  endif


end subroutine SetupPMApproach

! ************************************************************************** !

subroutine ExtractPMsFromPMList(simulation,pm_flow,pm_tran,pm_waste_form, &
                                pm_ufd_decay,pm_ufd_biosphere,pm_geop, &
                                pm_auxiliary,pm_well,pm_material_transform)
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
  use PM_ERT_class
  use PM_Auxiliary_class
  use PM_Well_class
  use PM_Material_Transform_class
  use Option_module
  use Simulation_Subsurface_class

  implicit none

  class(simulation_subsurface_type) :: simulation

  type(option_type), pointer :: option
  class(pm_subsurface_flow_type), pointer :: pm_flow
  class(pm_base_type), pointer :: pm_tran
  class(pm_waste_form_type), pointer :: pm_waste_form
  class(pm_ufd_decay_type), pointer :: pm_ufd_decay
  class(pm_ufd_biosphere_type), pointer :: pm_ufd_biosphere
  class(pm_base_type), pointer :: pm_geop
  class(pm_auxiliary_type), pointer :: pm_auxiliary
  class(pm_well_type), pointer :: pm_well
  class(pm_material_transform_type), pointer :: pm_material_transform
  class(pm_base_type), pointer :: cur_pm, prev_pm

  option => simulation%option

  nullify(pm_flow)
  nullify(pm_tran)
  nullify(pm_waste_form)
  nullify(pm_ufd_decay)
  nullify(pm_ufd_biosphere)
  nullify(pm_auxiliary)
  nullify(pm_well)
  nullify(pm_material_transform)

  cur_pm => simulation%process_model_list
  do
    if (.not.associated(cur_pm)) exit
    select type(cur_pm)
      class is(pm_subsurface_flow_type)
        pm_flow => cur_pm
      class is(pm_rt_type)
        pm_tran => cur_pm
      class is(pm_nwt_type)
        pm_tran => cur_pm
      class is(pm_waste_form_type)
        pm_waste_form => cur_pm
      class is(pm_ufd_decay_type)
        pm_ufd_decay => cur_pm
      class is(pm_ufd_biosphere_type)
        pm_ufd_biosphere => cur_pm
      class is(pm_ert_type)
        pm_geop => cur_pm
      class is(pm_auxiliary_type)
        pm_auxiliary => cur_pm
      class is(pm_well_type)
        pm_well => cur_pm
      class is(pm_material_transform_type)
        pm_material_transform => cur_pm
      class default
        option%io_buffer = &
         'PM Class unrecognized in FactorySubsurfaceInitPostPetsc.'
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

subroutine SetupPMCLinkages(simulation,pm_flow,pm_tran,pm_waste_form, &
                            pm_ufd_decay,pm_ufd_biosphere,pm_geop, &
                            pm_auxiliary,pm_well,pm_material_transform, &
                            realization)
  !
  ! Sets up all PMC linkages
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_Waste_Form_class
  use PM_UFD_Decay_class
  use PM_UFD_Biosphere_class
  use PM_Auxiliary_class
  use PM_Well_class
  use PM_Material_Transform_class
  use PM_WIPP_Flow_class
  use Factory_Subsurface_Read_module
  use Realization_Subsurface_class
  use Option_module
  use Input_Aux_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_subsurface_flow_type), pointer :: pm_flow
  class(pm_base_type), pointer :: pm_tran
  class(pm_waste_form_type), pointer :: pm_waste_form
  class(pm_ufd_decay_type), pointer :: pm_ufd_decay
  class(pm_ufd_biosphere_type), pointer :: pm_ufd_biosphere
  class(pm_base_type), pointer :: pm_geop
  class(pm_auxiliary_type), pointer :: pm_auxiliary
  class(pm_well_type), pointer :: pm_well
  class(pm_material_transform_type), pointer :: pm_material_transform
  class(realization_subsurface_type), pointer :: realization

  type(option_type), pointer :: option
  type(input_type), pointer :: input

  option => simulation%option

  if (associated(pm_flow)) &
    call AddPMCSubsurfaceFlow(simulation,pm_flow,'PMCSubsurfaceFlow', &
                              realization,option)

  if (associated(pm_tran))   &
    call AddPMCSubsurfaceTransport(simulation,pm_tran, &
                                   'PMCSubsurfaceTransport', &
                                   realization,option)

  if (associated(pm_geop)) &
    call AddPMCSubsurfaceGeophysics(simulation,pm_geop, &
                                    'PMCSubsurfaceGeophysics', &
                                    realization,option)

  input => InputCreate(IN_UNIT,option%input_filename,option)
  call SubsurfaceReadRequiredCards(simulation,input)
  call SubsurfaceReadInput(simulation,input)

  if (associated(pm_waste_form)) &
    call AddPMCWasteForm(simulation,pm_waste_form,'PMC3PWasteForm',&
                         associated(pm_ufd_decay),realization,input,option)

  if (associated(pm_ufd_decay)) &
    call AddPMCUFDDecay(simulation,pm_ufd_decay,'PMC3PUFDDecay',realization, &
                        input,option)

  if (associated(pm_ufd_biosphere)) &
    call AddPMCUDFBiosphere(simulation,pm_ufd_biosphere,'PMC3PUFDBiosphere',&
                            associated(pm_ufd_decay),realization,input,option)

  if (associated(pm_auxiliary)) &
    call AddPMCGeneral(simulation,pm_auxiliary,'SALINITY',realization,option)

  if (associated(pm_well)) &
    call AddPMCWell(simulation,pm_well,'PMCWell',realization,input,option)
  if (associated(pm_material_transform)) &
    call AddPMCMaterialTransform(simulation,pm_material_transform, &
                                 'PMC3MaterialTransform',realization,input, &
                                 option)

  if (associated(pm_flow)) then
    select type(pm_flow)
      class is (pm_wippflo_type)
        call AddPMWippSrcSink(realization,pm_flow,input)
    end select
  endif

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
  use PMC_Linear_class
  use Timestepper_TS_class
  use Timestepper_SNES_class
  use Timestepper_KSP_class
  use Timestepper_Steady_class
  use PM_TH_TS_class
  use PM_Richards_TS_class
  use PM_PNF_class
  use Realization_Subsurface_class
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

  select type(pm_flow)
    class is(pm_pnf_type)
      pmc_subsurface => PMCLinearCreate()
    class default
      pmc_subsurface => PMCSubsurfaceCreate()
  end select

  call pmc_subsurface%SetName(pmc_name)
  call pmc_subsurface%SetOption(option)
  call pmc_subsurface%SetWaypointList(simulation%waypoint_list_subsurface)

  pmc_subsurface%pm_list => pm_flow
  pmc_subsurface%pm_ptr%pm => pm_flow
  pmc_subsurface%realization => realization

  ! add time integrator
  if (pm_flow%steady_state) then
    option%flow%steady_state = PETSC_TRUE
    pmc_subsurface%timestepper => TimestepperSteadyCreate()
  else
    select type(pm_flow)
      class is(pm_richards_ts_type)
        pmc_subsurface%timestepper => TimestepperTSCreate()
      class is(pm_th_ts_type)
        pmc_subsurface%timestepper => TimestepperTSCreate()
      class is(pm_pnf_type)
        pmc_subsurface%timestepper => TimestepperKSPCreate()
      class default
        pmc_subsurface%timestepper => TimestepperSNESCreate()
    end select
  endif
  pmc_subsurface%timestepper%name = 'FLOW'

  ! add solver
  call pmc_subsurface%pm_list%InitializeSolver()
  pmc_subsurface%timestepper%solver => pmc_subsurface%pm_list%solver
  pmc_subsurface%timestepper%solver%itype = FLOW_CLASS

  ! set up logging stage
  string = trim(pm_flow%name)
  call LoggingCreateStage(string,pmc_subsurface%stage)
  simulation%flow_process_model_coupler => pmc_subsurface
  simulation%process_model_coupler_list => &
    simulation%flow_process_model_coupler

end subroutine AddPMCSubsurfaceFlow

! ************************************************************************** !

subroutine AddPMCSubsurfaceTransport(simulation,pm_base,pmc_name, &
                                     realization,option)
  !
  ! Adds a subsurface transport PMC
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/19
  !

  use PMC_Base_class
  use PM_Base_class
  use PM_RT_class
  use PM_NWT_class
  use PMC_Subsurface_class
  use PMC_Subsurface_OSRT_class
  use Timestepper_SNES_class
  use Timestepper_KSP_class
  use Timestepper_Steady_class
  use Realization_Subsurface_class
  use Option_module
  use Logging_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_base_type), pointer :: pm_base
  character(len=*) :: pmc_name
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option

  class(pmc_subsurface_type), pointer :: pmc_subsurface
  character(len=MAXSTRINGLENGTH) :: string
  class(pmc_base_type), pointer :: pmc_dummy

  nullify(pmc_dummy)

  select type(pm=>pm_base)
    class is(pm_rt_type)
      if (pm%operator_split) then
        pmc_subsurface => PMCSubsurfaceOSRTCreate()
      else
        pmc_subsurface => PMCSubsurfaceCreate()
      endif
    class is(pm_nwt_type)
      pmc_subsurface => PMCSubsurfaceCreate()
  end select
  call pmc_subsurface%SetName(pmc_name)
  call pmc_subsurface%SetOption(option)
  call pmc_subsurface%SetWaypointList(simulation%waypoint_list_subsurface)
  pmc_subsurface%pm_list => pm_base
  pmc_subsurface%pm_ptr%pm => pm_base
  pmc_subsurface%realization => realization

  ! add time integrator
  if (pm_base%steady_state) then
    option%transport%steady_state = PETSC_TRUE
    pmc_subsurface%timestepper => TimestepperSteadyCreate()
  else
    select type(pm=>pm_base)
      class is(pm_rt_type)
        if (pm%operator_split) then
          pmc_subsurface%timestepper => TimestepperKSPCreate()
        else
          pmc_subsurface%timestepper => TimestepperSNESCreate()
        endif
      class is(pm_nwt_type)
        pmc_subsurface%timestepper => TimestepperSNESCreate()
    end select
  endif
  pmc_subsurface%timestepper%name = 'TRAN'

  ! add solver
  call pmc_subsurface%pm_list%InitializeSolver()
  pmc_subsurface%timestepper%solver => pmc_subsurface%pm_list%solver
  pmc_subsurface%timestepper%solver%itype = TRANSPORT_CLASS

  ! set up logging stage
  string = trim(pm_base%name)
  call LoggingCreateStage(string,pmc_subsurface%stage)
  simulation%tran_process_model_coupler => pmc_subsurface

  if (.not.associated(simulation%process_model_coupler_list)) then
    simulation%process_model_coupler_list => pmc_subsurface
  else
    call PMCBaseSetChildPeerPtr(pmc_subsurface%CastToBase(),PM_CHILD, &
                      simulation%flow_process_model_coupler%CastToBase(), &
                      pmc_dummy,PM_INSERT)
  endif

end subroutine AddPMCSubsurfaceTransport

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

  if (option%itranmode /= RT_MODE .and. option%itranmode /= NWT_MODE) then
     option%io_buffer = 'The Waste Form process model requires &
          &a transport process model (GIRT/OSRT or NWT).'
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
  pmc_waste_form%pm_list => pm_waste_form
  pmc_waste_form%pm_ptr%pm => pm_waste_form
  pmc_waste_form%realization => realization

  ! set up logging stage
  string = 'WASTE_FORM_GENERAL'
  call LoggingCreateStage(string,pmc_waste_form%stage)
  call PMCBaseSetChildPeerPtr(pmc_waste_form%CastToBase(),PM_CHILD, &
         simulation%tran_process_model_coupler%CastToBase(), &
         pmc_dummy,PM_APPEND)

end subroutine AddPMCWasteForm

! ************************************************************************** !

subroutine AddPMCUFDDecay(simulation,pm_ufd_decay,pmc_name,&
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

  if (option%itranmode /= RT_MODE .and. option%itranmode /= NWT_MODE) then
     option%io_buffer = 'The UFD_DECAY process model requires a transport &
          &process model (GIRT/OSRT or NWT).'
     call PrintErrMsg(option)
  endif

  pmc_ufd_decay => PMCThirdPartyCreate()
  call pmc_ufd_decay%SetName(pmc_name)
  call pmc_ufd_decay%SetOption(option)
  call pmc_ufd_decay%SetWaypointList(simulation%waypoint_list_subsurface)
  pmc_ufd_decay%pm_list => pm_ufd_decay
  pmc_ufd_decay%pm_ptr%pm => pm_ufd_decay
  pmc_ufd_decay%realization => realization

  ! set up logging stage
  string = 'UFD_DECAY'
  call LoggingCreateStage(string,pmc_ufd_decay%stage)
  call PMCBaseSetChildPeerPtr(pmc_ufd_decay%CastToBase(),PM_CHILD, &
         simulation%tran_process_model_coupler%CastToBase(), &
         pmc_dummy,PM_APPEND)

end subroutine AddPMCUFDDecay

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
  if (option%itranmode /= RT_MODE) then
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
  call pmc_ufd_biosphere%SetWaypointList(simulation%waypoint_list_subsurface)
  pmc_ufd_biosphere%pm_list => pm_ufd_biosphere
  pmc_ufd_biosphere%pm_ptr%pm => pm_ufd_biosphere
  pmc_ufd_biosphere%realization => realization

  ! set up logging stage
  string = 'UFD_BIOSPHERE'
  call LoggingCreateStage(string,pmc_ufd_biosphere%stage)
  call PMCBaseSetChildPeerPtr(pmc_ufd_biosphere%CastToBase(),PM_CHILD, &
         simulation%tran_process_model_coupler%CastToBase(), &
         pmc_dummy,PM_APPEND)

end subroutine AddPMCUDFBiosphere

! ************************************************************************** !

subroutine AddPMCSubsurfaceGeophysics(simulation,pm_base,pmc_name, &
                                      realization,option)
  ! Adds a Geophysics PMC
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/25/21
  !

  use PM_Base_class
  use PM_ERT_class
  use PMC_Base_class
  use PMC_Geophysics_class
  use Timestepper_Steady_class
  use Realization_Subsurface_class
  use Option_module
  use Logging_module
  use Waypoint_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_base_type), pointer :: pm_base
  character(len=*) :: pmc_name
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option

  class(pmc_geophysics_type), pointer :: pmc_geophysics
  character(len=MAXSTRINGLENGTH) :: string

  pmc_geophysics => PMCGeophysicsCreate()

  call pmc_geophysics%SetName(pmc_name)
  call pmc_geophysics%SetOption(option)

  pmc_geophysics%pm_list => pm_base
  pmc_geophysics%pm_ptr%pm => pm_base
  pmc_geophysics%realization => realization

  ! add time integrator
  select type(pm=>pm_base)
    class is(pm_ert_type)
      pmc_geophysics%timestepper => TimestepperSteadyCreate()
      call WaypointListCopyAndMerge(simulation%waypoint_list_subsurface, &
                                    pm%waypoint_list,option)
      call pmc_geophysics%SetWaypointList(pm%waypoint_list)
    class default
      pmc_geophysics%timestepper => TimestepperSteadyCreate()
  end select
  pmc_geophysics%timestepper%name = 'GEOP'

  ! add solver
  call pmc_geophysics%pm_list%InitializeSolver()
  pmc_geophysics%timestepper%solver => pmc_geophysics%pm_list%solver
  pmc_geophysics%timestepper%solver%itype = GEOPHYSICS_CLASS

  ! set up logging stage
  string = trim(pm_base%name)
  call LoggingCreateStage(string,pmc_geophysics%stage)
  simulation%geop_process_model_coupler => pmc_geophysics
  if (associated(simulation%process_model_coupler_list)) then
    simulation%process_model_coupler_list%peer => pmc_geophysics
  else
    simulation%process_model_coupler_list => pmc_geophysics
  endif

end subroutine AddPMCSubsurfaceGeophysics

! ************************************************************************** !

subroutine AddPMCGeneral(simulation,pm_auxiliary,pmc_name, &
                         realization,option)

  !
  ! Adds an auxiliary PMC
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PMC_Base_class
  use PM_Auxiliary_class
  use PMC_General_class
  use PMC_Subsurface_class
  use Realization_Subsurface_class
  use Option_module
  use String_module
  use Logging_module
  use Input_Aux_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_auxiliary_type), pointer :: pm_auxiliary
  character(len=*) :: pmc_name
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option

  class(pmc_general_type), pointer :: pmc_general
  character(len=MAXSTRINGLENGTH) :: string
  class(pmc_base_type), pointer :: pmc_dummy

  nullify(pmc_dummy)

  pm_auxiliary%realization => realization

  string = 'salinity'
  if (StringCompareIgnoreCase(pm_auxiliary%ctype,string)) then
    if (option%itranmode == RT_MODE) then
      pmc_general => PMCGeneralCreate(pmc_name,pm_auxiliary%CastToBase())
      call PMCBaseSetChildPeerPtr(pmc_general%CastToBase(),PM_PEER, &
             simulation%tran_process_model_coupler%CastToBase(), &
             pmc_dummy,PM_APPEND)
    else
      option%io_buffer = 'Reactive transport must be included in the &
           &SIMULATION block in order to use the SALINITY process model.'
      call PrintErrMsg(option)
    endif
  endif

  call LoggingCreateStage(string,pmc_general%stage)

end subroutine AddPMCGeneral

! ************************************************************************** !

subroutine AddPMCMaterialTransform(simulation, pm_material_transform, pmc_name,&
                                   realization, input, option)
  !
  ! Adds a material transform PMC
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022
  !

  use PMC_Base_class
  use PMC_Third_Party_class
  use PM_Material_Transform_class
  use Realization_Subsurface_class
  use Option_module
  use Logging_module
  use Input_Aux_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_material_transform_type), pointer :: pm_material_transform
  character(len=*) :: pmc_name
  class(realization_subsurface_type), pointer :: realization
  type(input_type), pointer :: input
  type(option_type), pointer :: option

  class(pmc_third_party_type), pointer :: pmc_material_transform
  character(len=MAXSTRINGLENGTH) :: string
  class(pmc_base_type), pointer :: pmc_dummy

  nullify(pmc_dummy)

  string = 'MATERIAL_TRANSFORM_GENERAL'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  call pm_material_transform%ReadPMBlock(input)

  pmc_material_transform => PMCThirdPartyCreate()
  call pmc_material_transform%SetName(pmc_name)
  call pmc_material_transform%SetOption(option)
  call pmc_material_transform%SetWaypointList(simulation&
                                                %waypoint_list_subsurface)
  pmc_material_transform%pm_list => pm_material_transform
  pmc_material_transform%pm_ptr%pm => pm_material_transform
  pmc_material_transform%realization => realization

  ! set up logging stage
  string = 'MATERIAL_TRANSFORM_GENERAL'
  call LoggingCreateStage(string,pmc_material_transform%stage)

  ! Material transform is child of flow and peer of transport
  if (associated(simulation%tran_process_model_coupler) .and. &
      associated(simulation%flow_process_model_coupler)) then
    call PMCBaseSetChildPeerPtr(pmc_material_transform%CastToBase(), &
           PM_CHILD,simulation%flow_process_model_coupler%CastToBase(), &
           simulation%tran_process_model_coupler%CastToBase(),PM_INSERT)
  elseif(associated(simulation%flow_process_model_coupler)) then
    call PMCBaseSetChildPeerPtr(pmc_material_transform%CastToBase(), &
           PM_CHILD,simulation%flow_process_model_coupler%CastToBase(), &
           pmc_dummy,PM_INSERT)
  elseif(associated(simulation%tran_process_model_coupler)) then
    call PMCBaseSetChildPeerPtr(pmc_material_transform%CastToBase(), &
           PM_PEER,simulation%tran_process_model_coupler%CastToBase(), &
           pmc_dummy,PM_APPEND)
  endif

end subroutine AddPMCMaterialTransform

! ************************************************************************** !

subroutine AddPMCWell(simulation,pm_well,pmc_name,realization,input, &
                      option)

  !
  ! Adds a well PMC
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  !

  use PMC_Base_class
  use PMC_Third_Party_class
  use PM_Well_class
  use Realization_Subsurface_class
  use Option_module
  use Logging_module
  use Input_Aux_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_well_type), pointer :: pm_well
  character(len=*) :: pmc_name
  class(realization_subsurface_type), pointer :: realization
  type(input_type), pointer :: input
  type(option_type), pointer :: option

  class(pmc_third_party_type), pointer :: pmc_well
  character(len=MAXSTRINGLENGTH) :: string
  class(pmc_base_type), pointer :: pmc_dummy

  nullify(pmc_dummy)

  string = 'WELLBORE_MODEL'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  call pm_well%ReadPMBlock(input)

  if (option%iflowmode /= WF_MODE) then
     option%io_buffer = 'The WELLBORE_MODEL process model can only be &
                        &used with WIPP_FLOW mode at the moment.'
     call PrintErrMsg(option)
  endif
  if ( (option%itranmode /= NULL_MODE) .and. &
       (option%itranmode /= NWT_MODE) ) then
       option%io_buffer = 'The WELLBORE_MODEL process model can only be &
                        &used with NWT mode at the moment.'
     call PrintErrMsg(option)
  endif

  pmc_well => PMCThirdPartyCreate()
  call pmc_well%SetName(pmc_name)
  call pmc_well%SetOption(option)
  call pmc_well%SetWaypointList(simulation%waypoint_list_subsurface)
  pmc_well%pm_list => pm_well
  pmc_well%pm_ptr%pm => pm_well
  pmc_well%realization => realization

  ! set up logging stage
  string = 'WELLBORE_MODEL'
  call LoggingCreateStage(string,pmc_well%stage)

  if ( (option%itranmode /= NULL_MODE) .and. &
       (option%itranmode == NWT_MODE) ) then
    call PMCBaseSetChildPeerPtr(pmc_well%CastToBase(),PM_CHILD, &
         simulation%tran_process_model_coupler%CastToBase(), &
         pmc_dummy,PM_APPEND)
  else
    call PMCBaseSetChildPeerPtr(pmc_well%CastToBase(),PM_CHILD, &
         simulation%flow_process_model_coupler%CastToBase(), &
         pmc_dummy,PM_APPEND)
  endif

end subroutine AddPMCWell

! ************************************************************************** !

subroutine AddPMWippSrcSink(realization,pm_wippflo,input)

  use Input_Aux_module
  use Option_module
  use Realization_Subsurface_class
  use WIPP_Flow_Aux_module
  use PM_WIPP_Flow_class
  use PM_WIPP_SrcSink_class

  implicit none

  class(realization_subsurface_type), pointer :: realization
  class(pm_wippflo_type) :: pm_wippflo
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: block_string

  option => realization%option

  block_string = 'WIPP_SOURCE_SINK'
  call InputFindStringInFile(input,option,block_string)
  if (input%ierr == 0 .and. wippflo_use_gas_generation) then
    pm_wippflo%pmwss_ptr => PMWSSCreate()
    pm_wippflo%pmwss_ptr%option => option
    call pm_wippflo%pmwss_ptr%ReadPMBlock(input)
    call PMWSSSetRealization(pm_wippflo%pmwss_ptr,realization)
  endif

end subroutine AddPMWippSrcSink

! ************************************************************************** !

subroutine FactorySubsurfSetPMCWaypointPtrs(simulation)
  !
  ! Sets the process model coupler waypoint pointers to the first waypoint
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/22

  implicit none

  class(simulation_subsurface_type) :: simulation

  if (associated(simulation%flow_process_model_coupler)) then
    call simulation%flow_process_model_coupler% &
           SetWaypointPtr(simulation%waypoint_list_subsurface)
  endif
  if (associated(simulation%tran_process_model_coupler)) then
    call simulation%tran_process_model_coupler% &
           SetWaypointPtr(simulation%waypoint_list_subsurface)
  endif
  if (associated(simulation%geop_process_model_coupler)) then
    call simulation%geop_process_model_coupler% &
           SetWaypointPtr(simulation%waypoint_list_subsurface)
  endif

end subroutine FactorySubsurfSetPMCWaypointPtrs

end module Factory_Subsurface_Linkage_module
