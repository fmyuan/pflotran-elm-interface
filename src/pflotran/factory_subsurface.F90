module Factory_Subsurface_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Simulation_Subsurface_class

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal

  implicit none

  private

  public :: FactorySubsurfaceInitialize, &
            FactorySubsurfaceInitPostPetsc, &
            FactorySubsurfaceJumpStart

contains

! ************************************************************************** !

subroutine FactorySubsurfaceInitialize(simulation)
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
  call FactorySubsurfaceInitPostPetsc(simulation)

end subroutine FactorySubsurfaceInitialize

! ************************************************************************** !

subroutine FactorySubsurfaceInitPostPetsc(simulation)
  !
  ! Sets up PFLOTRAN subsurface simulation
  ! framework after to PETSc initialization
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/13
  !

  use Option_module
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_Waste_Form_class
  use PM_UFD_Decay_class
  use PM_UFD_Biosphere_class
  use PM_Auxiliary_class
  use PM_Well_class
  use PM_Material_Transform_class
  use Factory_Subsurface_Linkage_module
  use Realization_Subsurface_class
  use Simulation_Subsurface_class
  use Waypoint_module

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
  class(realization_subsurface_type), pointer :: realization

  option => simulation%option

  nullify(pm_flow)
  nullify(pm_tran)
  nullify(pm_waste_form)
  nullify(pm_ufd_decay)
  nullify(pm_ufd_biosphere)
  nullify(pm_geop)
  nullify(pm_auxiliary)
  nullify(pm_well)

  ! process command line arguments specific to subsurface
  call FactSubInitCommandLineSettings(option)

  call FactSubLinkExtractPMsFromPMList(simulation,pm_flow,pm_tran, &
                                       pm_waste_form,pm_ufd_decay, &
                                       pm_ufd_biosphere,pm_geop, &
                                       pm_auxiliary,pm_well, &
                                       pm_material_transform)

  call FactorySubsurfaceSetFlowMode(pm_flow,option)
  call FactorySubsurfaceSetGeopMode(pm_geop,option)

  realization => RealizationCreate(option)
  simulation%realization => realization
  realization%output_option => simulation%output_option

  ! Setup linkages between PMCs
  call FactSubLinkSetupPMCLinkages(simulation,pm_flow,pm_tran, &
                                       pm_waste_form,pm_ufd_decay, &
                                       pm_ufd_biosphere,pm_geop, &
                                       pm_auxiliary,pm_well, &
                                       pm_material_transform)

  ! SubsurfaceInitSimulation() must be called after pmc linkages are set above.
  call FactorySubsurfaceInitSimulation(simulation)

  ! set first process model coupler as the master
  simulation%process_model_coupler_list%is_master = PETSC_TRUE

end subroutine FactorySubsurfaceInitPostPetsc

! ************************************************************************** !

subroutine FactSubInitCommandLineSettings(option)
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

end subroutine FactSubInitCommandLineSettings

! ************************************************************************** !

subroutine FactorySubsurfaceSetFlowMode(pm_flow,option)
  !
  ! Sets the flow mode (richards, vadose, mph, etc.)
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  use Option_module
  use PM_Subsurface_Flow_class
  use PM_Base_class
  use PM_General_class
  use PM_Hydrate_class
  use PM_WIPP_Flow_class
  use PM_Mphase_class
  use PM_Richards_class
  use PM_TH_class
  use PM_Richards_TS_class
  use PM_TH_TS_class
  use PM_ZFlow_class
  use ZFlow_Aux_module
  use PM_PNF_class
  use General_Aux_module

  implicit none

  type(option_type) :: option
  class(pm_subsurface_flow_type), pointer :: pm_flow

  option%liquid_phase = 1
  option%gas_phase = 2 ! always set gas phase to 2 for transport

  if (.not.associated(pm_flow)) then
    option%nphase = 1
    ! assume default isothermal when only transport
    option%use_isothermal = PETSC_TRUE
    return
  endif

  select type(pm_flow)
    class is (pm_wippflo_type)
      option%iflowmode = WF_MODE
      option%nphase = 2
      option%capillary_pressure_id = 3
      option%saturation_pressure_id = 4
      option%water_id = 1
      option%air_id = 2
      option%nflowdof = 2
      option%nflowspec = 2
    class is (pm_general_type)
      call PMGeneralSetFlowMode(pm_flow,option)
    class is (pm_hydrate_type)
      call PMHydrateSetFlowMode(option)
    class is (pm_mphase_type)
      option%iflowmode = MPH_MODE
      option%nphase = 2
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
      option%nflowdof = 1
      option%nflowspec = 1
      option%use_isothermal = PETSC_TRUE
    class is (pm_zflow_type)
      option%iflowmode = ZFLOW_MODE
      option%nphase = 1
      option%nflowdof = 0
      option%nflowspec = 0
      if (Initialized(zflow_liq_flow_eq)) then
        option%nflowdof = option%nflowdof + 1
        option%nflowspec = option%nflowspec + 1
      endif
      if (Initialized(zflow_heat_tran_eq)) then
        option%nflowdof = option%nflowdof + 1
      else
        option%use_isothermal = PETSC_TRUE
      endif
      if (Initialized(zflow_sol_tran_eq)) then
        option%nflowdof = option%nflowdof + 1
        option%nflowspec = 1
      endif
      if (option%nflowdof == 0) then
        option%io_buffer=  'A process must be specified under ZFLOW,&
          &OPTIONS,PROCESSES.'
        call PrintErrMsg(option)
      endif
    class is (pm_pnf_type)
      option%iflowmode = PNF_MODE
      option%nphase = 1
      option%nflowdof = 1
      option%nflowspec = 1
      option%use_isothermal = PETSC_TRUE
    class is (pm_th_type)
      option%iflowmode = TH_MODE
      option%nphase = 1
      option%nflowdof = 2
      option%nflowspec = 1
      option%use_isothermal = PETSC_FALSE
      option%flow%store_fluxes = PETSC_TRUE
    class is (pm_richards_ts_type)
      option%iflowmode = RICHARDS_TS_MODE
      option%nphase = 1
      option%nflowdof = 1
      option%nflowspec = 1
      option%use_isothermal = PETSC_TRUE
    class is (pm_th_ts_type)
      option%iflowmode = TH_TS_MODE
      option%nphase = 1
      option%nflowdof = 2
      option%nflowspec = 1
      option%use_isothermal = PETSC_FALSE
      option%flow%store_fluxes = PETSC_TRUE
    class default
      option%io_buffer = ''
      call PrintErrMsg(option)

  end select

  if (option%nflowdof == 0) then
    option%io_buffer = 'Number of flow degrees of freedom is zero.'
    call PrintErrMsg(option)
  endif
  if (option%nphase == 0) then
    option%io_buffer = 'Number of flow phases is zero.'
    call PrintErrMsg(option)
  endif
  if (option%nflowspec == 0) then
    option%io_buffer = 'Number of flow species is zero.'
    call PrintErrMsg(option)
  endif

end subroutine FactorySubsurfaceSetFlowMode

! ************************************************************************** !

subroutine FactorySubsurfaceSetGeopMode(pm_geop,option)
  !
  ! Sets the geophysics mode (ert, sip, etc.)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/26/21
  !

  use Option_module
  use PM_Base_class
  use PM_ERT_class
  !use General_Aux_module

  implicit none

  type(option_type) :: option
  class(pm_base_type), pointer :: pm_geop

  if (.not.associated(pm_geop)) then
    return
  endif

  select type(pm_geop)
    class is (pm_ert_type)
      option%igeopmode = ERT_MODE
      option%geopmode = "ERT"
      option%ngeopdof = 1
    class default
      option%io_buffer = ''
      call PrintErrMsg(option)
  end select

end subroutine FactorySubsurfaceSetGeopMode

! ************************************************************************** !

subroutine FactorySubsurfaceInitSimulation(simulation)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Realization_Subsurface_class
  use Realization_Base_class
  use Discretization_module
  use Option_module
  use Output_module, only : Output
  use Output_Aux_module
  use Global_module
  use Factory_Subsurface_Linkage_module
  use Init_Subsurface_module
  use Init_Subsurface_Flow_module
  use Init_Subsurface_Tran_module
  use Init_Subsurface_Geop_module
  use Init_Common_module
  use Waypoint_module
  use Strata_module
  use Regression_module
  use PMC_Subsurface_class
  use PMC_General_class
  use PMC_Base_class
  use PM_Base_Pointer_module
  use PM_Inversion_class
  use PM_Subsurface_Flow_class
  use Timestepper_SNES_class
  use Waypoint_module

  implicit none

  class(simulation_subsurface_type) :: simulation

  class(pmc_base_type), pointer :: cur_process_model_coupler_top

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option

  realization => simulation%realization
  option => realization%option

! begin from old Init()
  call FactorySubsurfSetupRealization(simulation)

  call InitCommonAddOutputWaypoints(option,simulation%output_option, &
                                    simulation%waypoint_list_subsurface)

  ! initialize global auxiliary variable object
  call GlobalSetup(realization)

  ! always call the flow side since a velocity field still has to be
  ! set if no flow exists
  call InitSubsurfFlowSetupRealization(simulation)
  if (option%ntrandof > 0) then
    call InitSubsurfTranSetupRealization(realization)
  endif
  if (option%ngeopdof > 0) then
    call InitSubsurfGeopSetupRealization(realization)
  endif

  ! InitSubsurfaceSetupZeroArray must come after InitSubsurfaceXXXRealization
  call InitSubsurfaceSetupZeroArrays(realization)
  call OutputVariableAppendDefaults(realization%output_option% &
                                      output_snap_variable_list,option)

  call RegressionSetup(simulation%regression,realization)
! end from old Init()

  call DiscretizationPrintInfo(realization%discretization, &
                               realization%patch%grid,option)

  !----------------------------------------------------------------------------!
  ! This section for setting up new process model approach
  !----------------------------------------------------------------------------!
  call FactSubLinkAddPMCEvolvingStrata(simulation)
  call FactSubLinkAddPMCInversion(simulation)

  ! For each ProcessModel, set:
  ! - realization (subsurface or surface),
  ! - stepper (flow/trans/surf_flow),
  ! For each ProcessModelCoupler, set:
  ! - SNES functions (Residual/Jacobian), or TS function (RHSFunction)

  cur_process_model_coupler_top => simulation%process_model_coupler_list

  ! the following recursive subroutine will also call each pmc child
  ! and each pms's peers
  if (associated(cur_process_model_coupler_top)) then
    call FactSubLinkSetupPMApproach(cur_process_model_coupler_top, &
                                    simulation)
  endif

  ! point the top process model coupler to Output
  simulation%process_model_coupler_list%Output => Output

  ! setup the outer waypoint lists
  call FactorySubsurfSetupWaypointList(simulation)
  call FactSubLinkSetPMCWaypointPtrs(simulation)

  if (realization%debug%print_couplers) then
    call InitCommonVerifyAllCouplers(realization)
  endif

  call FactorySubsurfaceJumpStart(simulation)

end subroutine FactorySubsurfaceInitSimulation

! ************************************************************************** !

subroutine FactorySubsurfSetupRealization(simulation)
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

  select case(option%itranmode)
    case(RT_MODE)
      ! read reaction database
      if (realization%reaction%use_full_geochemistry) then
        call DatabaseRead(realization%reaction,option)
        call BasisInit(realization%reaction,option)
      else
        ! turn off activity coefficients since the database has not been read
        realization%reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
        !TODO(jenn) Should I turn on print here too?
        allocate(realization%reaction%primary_species_print(option%ntrandof))
        realization%reaction%primary_species_print = PETSC_TRUE
      endif

      ! SK 09/30/13, Added to check if Mphase is called with OS
      if (option%transport%reactive_transport_coupling == OPERATOR_SPLIT .and. &
          option%iflowmode == MPH_MODE) then
        option%io_buffer = 'Operator splitting currently not implemented with &
                   &MPHASE. Please switch reactive transport to MODE GIRT.'
        call PrintErrMsg(option)
        option%transport%reactive_transport_coupling = GLOBAL_IMPLICIT
      endif
    case(NWT_MODE)
  end select

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
  ! add waypoints associated with boundary conditions, source/sinks etc. to list
  call RealizationAddWaypointsToList(realization, &
                                     simulation%waypoint_list_subsurface)
  ! fill in holes in waypoint data
  if (option%ngeopdof > 0) then
    ! Read geophysics survey file
    call RealizationReadGeopSurveyFile(realization)
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

end subroutine FactorySubsurfSetupRealization

! ************************************************************************** !

subroutine FactorySubsurfSetupWaypointList(simulation)
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
  ! creates a copy of outer and merges to subsurface
  call WaypointListCopyAndMerge(simulation%waypoint_list_subsurface, &
                                simulation%waypoint_list_outer,option)

  ! add sync waypoints into outer list
  call WaypointListMerge(simulation%waypoint_list_outer,sync_waypoint_list, &
                         option)

  ! add in periodic time waypoints for checkpointing. these will not appear
  ! in the outer list
  call CheckpointPeriodicTimeWaypoints(simulation%waypoint_list_subsurface, &
                                       option)
 ! fill in holes in waypoint data
  call WaypointListFillIn(simulation%waypoint_list_subsurface,option)
  call WaypointListRemoveExtraWaypnts(simulation%waypoint_list_subsurface, &
                                      option)
  call WaypointListFindDuplicateTimes(simulation%waypoint_list_subsurface, &
                                      option)

  ! debugging output
  if (realization%debug%print_waypoints) then
    call WaypointListPrint(simulation%waypoint_list_subsurface,option, &
                           realization%output_option)
  endif

end subroutine FactorySubsurfSetupWaypointList

! ************************************************************************** !

subroutine FactorySubsurfaceJumpStart(simulation)
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

  PetscBool :: failure
  PetscErrorCode :: ierr

  realization => simulation%realization
  option => realization%option

  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           "-vecload_block_size",failure,ierr);CHKERRQ(ierr)

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

end subroutine FactorySubsurfaceJumpStart

end module Factory_Subsurface_module
