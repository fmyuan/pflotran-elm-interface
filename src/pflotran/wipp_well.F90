module WIPP_Well_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Well_class
  use PM_Base_class
  use Option_module
  use Geometry_module
  use Strata_module
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module
  use WIPP_Flow_Aux_module
  use NW_Transport_Aux_module
  use Well_Grid_module
  use Condition_module

  implicit none

  private

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !    TOP
  !   -------- ql_bc(2)
  !   | i=n|
  !   -------- liq%ql(k=n-1)
  !   |    |
  !   ------
  !   |    |
  !   -------- liq%ql(k=3)
  !   | i=3|
  !   -------- liq%ql(k=2)
  !   | i=2|
  !   -------- liq%ql(k=1)
  !   | i=1|
  !   -------- ql_bc(1)
  !    BOT
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PetscBool :: initialize_well_tran = PETSC_TRUE

  ! srcsink vector indexing
  PetscInt, parameter :: UNPERT = 0
  PetscInt, parameter :: PERT_WRT_PL = 1
  PetscInt, parameter :: PERT_WRT_SG = 2

  ! WIPP Quasi-Implicit Well Model
  type, public, extends(pm_well_qi_type) :: pm_well_wipp_qi_type
  contains
    procedure, public :: ReadPMBlock => PMWellReadPMBlockWIPPQI
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMWellReadSimOptionsBlockWIPPQI
    procedure, public :: Setup => PMWellSetupWIPPQI
    procedure, public :: Solve => PMWellSolveWIPPQI
    procedure, public :: SolveFlow => PMWellSolveFlowWIPPQI
    procedure, public :: ModifyFlowResidual => PMWellModifyFlowResWIPPQI
    procedure, public :: ModifyFlowJacobian => PMWellModifyFlowJacWIPPQI
    procedure, public :: UpdateFlowRates => PMWellUpdateFlowRatesWIPPQI
    procedure, public :: InitializeTimestep => PMWellInitializeTimestepWIPPQI
    procedure, public :: UpdateFlowProperties => UpdateFlowPropertiesWIPPQI
  end type pm_well_wipp_qi_type

  public :: PMWellWIPPQICreate, &
            PMWellQISolveTran

  contains

! ************************************************************************** !

function PMWellWIPPQICreate()
  !
  ! Creates the WIPP quasi-implicit well process model.
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  implicit none

  class(pm_well_wipp_qi_type), pointer :: PMWellWIPPQICreate
  class(pm_well_wipp_qi_type), pointer :: pm_well

  allocate(pm_well)
  call PMWellQIInit(pm_well)

  pm_well%intrusion_time_start = UNINITIALIZED_DOUBLE
  pm_well%bh_zero_value = 1.d-20

  pm_well%well%well_model_type = WELL_MODEL_WIPP_QI
  pm_well%flow_coupling = QUASI_IMPLICIT_WELL

  nullify(pm_well%next_well)

  PMWellWIPPQICreate => pm_well

end function PMWellWIPPQICreate

! ************************************************************************** !

subroutine PMWellSetupWIPP(pm_well)

  use Grid_module
  use Utility_module
  use String_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Input_Aux_module
  use Dataset_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Transport_Constraint_NWT_module
  use NW_Transport_Aux_module
  use Global_Aux_module
  use Reactive_Transport_Aux_module

  implicit none

  class(pm_well_sequential_type) :: pm_well

  type(option_type), pointer :: option
  type(grid_type), pointer :: res_grid
  type(well_grid_type), pointer :: well_grid
  type(coupler_type), pointer :: source_sink
  type(input_type) :: input_dummy
  class(realization_subsurface_type), pointer :: realization
  class(tran_constraint_coupler_nwt_type), pointer ::tran_constraint_coupler_nwt
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  PetscInt :: num_new_source_sinks, offset

  call PMWellSetupBase(pm_well)
  option => pm_well%option
  realization => pm_well%realization
  res_grid => realization%patch%grid
  well_grid => pm_well%well_grid

  pm_well%flow_soln%ndof = pm_well%nphase

  if (option%itranmode /= NULL_MODE) then
    pm_well%transport = PETSC_TRUE
    if (option%itranmode /= NWT_MODE) then
      option%io_buffer ='The only transport mode allowed with the &
      &WIPP WELLBORE_MODEL is NWT_MODE.'
      call PrintErrMsg(option)
    endif
    pm_well%nspecies = realization%reaction_nw%params%nspecies
    pm_well%tran_soln%ndof = pm_well%nspecies
  endif

  ! add a reservoir src/sink coupler for each well segment
  num_new_source_sinks = 0
  if (associated(realization%patch%aux%Global%auxvars_ss)) then
    offset = realization%patch%aux%Global%num_aux_ss
  else
    offset = 0
  endif
  do k = 1,well_grid%nsegments

    source_sink => CouplerCreate()
    source_sink%itype = SRC_SINK_COUPLER_TYPE
    source_sink%name = trim(pm_well%name) // '_well_segment_' // StringWrite(k)

    ! ----- flow ------------------
    source_sink%flow_condition_name = trim(pm_well%name) // '_well_segment_' // &
                                      StringWrite(k) // '_flow_srcsink'
    source_sink%flow_condition => FlowConditionCreate(option)
    source_sink%flow_condition%name = source_sink%flow_condition_name
    if (well_grid%h_rank_id(k) /= option%myrank) cycle

    source_sink%flow_condition%general => FlowGeneralConditionCreate(option)
    string = 'RATE'
    source_sink%flow_condition%general%rate => FlowGeneralSubConditionPtr( &
      input_dummy,string,source_sink%flow_condition%general,option)
    source_sink%flow_condition%general%rate%itype = SCALED_MASS_RATE_SS
    source_sink%flow_condition%general%liquid_pressure => &
          FlowGeneralSubConditionPtr(input_dummy,string,source_sink% &
                                      flow_condition%general,option)
    source_sink%flow_condition%general%gas_pressure => &
          FlowGeneralSubConditionPtr(input_dummy,string,source_sink% &
                                    flow_condition%general,option)
    allocate(source_sink%flow_condition%general%rate%dataset%rarray(4))
    source_sink%flow_condition%general%rate%dataset%rarray(:) = 0.d0

    source_sink%flow_condition%well => FlowSubConditionCreate(ONE_INTEGER)

    ! ----- transport -------------
    if (pm_well%transport) then
      source_sink%tran_condition_name = trim(pm_well%name) // &
                          '_well_segment_' // StringWrite(k) // '_tran_srcsink'
      source_sink%tran_condition => TranConditionCreate(option)
      source_sink%tran_condition%name = source_sink%tran_condition_name
      source_sink%tran_condition%itype = WELL_SS
      tran_constraint_coupler_nwt => TranConstraintCouplerNWTCreate(option)
      allocate(tran_constraint_coupler_nwt%nwt_auxvar)
      call NWTAuxVarInit(tran_constraint_coupler_nwt%nwt_auxvar,&
                         pm_well%realization%reaction_nw,option)
      tran_constraint_coupler_nwt%constraint_name = trim(pm_well%name) //  &
                                            '_well_segment_' // StringWrite(k)
      source_sink%tran_condition%cur_constraint_coupler => &
                                                   tran_constraint_coupler_nwt
    endif

    source_sink%connection_set => &
      ConnectionCreate(1,GENERIC_CONNECTION_TYPE,res_grid%itype)
    source_sink%connection_set%id_dn = well_grid%h_local_id(k)
    if (well_grid%h_local_id(k) < 0) then
      source_sink%connection_set%num_connections = 0
    else
      num_new_source_sinks = num_new_source_sinks + 1
      source_sink%connection_set%offset = num_new_source_sinks + offset - 1
    endif

    call CouplerAddToList(source_sink,pm_well%realization%patch%source_sink_list)
    nullify(source_sink)
  enddo

end subroutine PMWellSetupWIPP

! ************************************************************************** !

subroutine PMWellSetupWIPPQI(this)
  !
  ! Initializes variables associated with the base well process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_wipp_qi_type) :: this

  call PMWellSetupWIPP(this)

end subroutine PMWellSetupWIPPQI

! ************************************************************************** !

subroutine PMWellReadSimOptionsBlockWIPPQI(this,input)
  !
  ! Author: Michael Nole
  ! Date: 03/08/24
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  class(pm_well_wipp_qi_type) :: this

  call PMWellReadSimOptionsBlockBase(this,input)

end subroutine PMWellReadSimOptionsBlockWIPPQI

! ************************************************************************** !

subroutine PMWellReadPMBlockWIPP(pm_well,input)
  !
  ! Reads input file parameters associated with the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  class(pm_well_sequential_type) :: pm_well

  type(option_type), pointer :: option
  type(well_grid_type), pointer :: well_grid
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt, pointer :: temp_seg_nums(:)
  PetscBool :: found
  PetscInt :: k, num_read
  PetscInt :: read_max = 10000

  option => pm_well%option
  well_grid => pm_well%well_grid
  input%ierr = INPUT_ERROR_NONE
  error_string = 'WELLBORE_MODEL'

  option%io_buffer = 'pflotran card:: WELLBORE_MODEL'
  call PrintMsg(option)

  allocate(temp_seg_nums(read_max))

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    found = PETSC_FALSE

    ! Read keywords within WELLBORE_MODEL block:
    select case(trim(word))
    !-------------------------------------
      case('SKIP_RESTART')
        pm_well%skip_restart = PETSC_TRUE
        cycle
    !-------------------------------------
      case('SINGLE_PHASE')
        pm_well%nphase = 1
        cycle
    !-------------------------------------
      case('TWO_PHASE')
        pm_well%nphase = 2
        cycle
    !-------------------------------------
      case('PRINT_WELL_FILE')
        pm_well%print_well = PETSC_TRUE
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) cycle
        call StringToUpper(word)
        if (StringCompare(word,'SEGMENTS')) then
          ! count the segment numbers
          num_read = 0
          do k = 1,read_max
            call InputReadInt(input,option,temp_seg_nums(k))
            if (InputError(input)) exit
            if (temp_seg_nums(k) <= 0) then
              option%io_buffer = 'A value provided for SEGMENTS &
                &after the ' // trim(error_string) // ', PRINT_WELL_FILE &
                &keyword was 0 or negative. Only positive integers are valid.'
              call PrintErrMsg(option)
            endif
            num_read = num_read + 1
          enddo
          if (num_read == 0) then
            option%io_buffer = 'At least one value for SEGMENTS &
              &must be provided after the ' // trim(error_string) // ', &
              &PRINT_WELL_FILE keyword.'
            call PrintErrMsg(option)
          endif
        else
          option%io_buffer = 'Keyword ' // trim(word) // ' not recognized &
            &after the ' // trim(error_string) // ', PRINT_WELL_FILE keyword&
            &. Did you mean "SEGMENTS"?'
          call PrintErrMsg(option)
        endif
        allocate(pm_well%well%segments_for_output(num_read))
        pm_well%well%segments_for_output(1:num_read) = temp_seg_nums(1:num_read)
        cycle
    !-------------------------------------
      case('SPLIT_WELL_FILE')
        pm_well%split_output_file = PETSC_TRUE
        cycle
    !-------------------------------------
      case('CHECK_FOR_SS')
        pm_well%ss_check = PETSC_TRUE
        cycle
    !-------------------------------------
      case('WIPP_INTRUSION_START_TIME')
        call InputReadDouble(input,option,pm_well%intrusion_time_start)
        call InputErrorMsg(input,option,'WIPP_INTRUSION_START_TIME', &
                           error_string)
        call InputReadAndConvertUnits(input,pm_well%intrusion_time_start,'sec', &
                           'WELLBORE_MODEL, WIPP_INTRUSION_START_TIME',option)
        pm_well%well_on = PETSC_FALSE
        cycle
    !-------------------------------------
      case('WIPP_INTRUSION_ZERO_VALUE')  ! [mol/m3-bulk]
        call InputReadDouble(input,option,pm_well%bh_zero_value)
        call InputErrorMsg(input,option,'WIPP_INTRUSION_ZERO_VALUE', &
                           error_string)
        cycle
    !-------------------------------------
    end select

    ! Read sub-blocks within WELLBORE_MODEL block:
    error_string = 'WELLBORE_MODEL'
    call PMWellReadGrid(well_grid,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadWell(pm_well,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadWellBCs(pm_well,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadFlowSolver(pm_well,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadTranSolver(pm_well,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadWellConstraintType(pm_well,input,option,word,error_string, &
                                      found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadWellOutput(pm_well,input,option,word,error_string,found)
    if (found) cycle

    if (.not. found) then
      option%io_buffer = 'Keyword "' // trim(word) // &
                         '" does not exist for WELLBORE_MODEL.'
      call PrintErrMsg(option)
    endif

  enddo
  call InputPopBlock(input,option)

  deallocate(temp_seg_nums)

  if (Initialized(pm_well%well%bh_p)) pm_well%flow_soln%bh_p = PETSC_TRUE
  if (pm_well%well%bh_p_set_by_reservoir) pm_well%flow_soln%bh_p = PETSC_TRUE
  if (Initialized(pm_well%well%th_p)) pm_well%flow_soln%th_p = PETSC_TRUE

  if (Initialized(pm_well%well%bh_ql)) pm_well%flow_soln%bh_q = PETSC_TRUE
  if (Initialized(pm_well%well%th_ql)) pm_well%flow_soln%th_q = PETSC_TRUE

end subroutine PMWellReadPMBlockWIPP

! ************************************************************************** !

subroutine PMWellReadPMBlockWIPPQI(this,input)

  use Input_Aux_module

  implicit none

  class(pm_well_wipp_qi_type) :: this
  type(input_type), pointer :: input

  call PMWellReadPMBlockWIPP(this,input)

end subroutine PMWellReadPMBlockWIPPQI

! ************************************************************************** !

subroutine PMWellInitializeTimestepWIPP(pm_well)
  !
  ! Initializes and takes the time step for the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscReal :: curr_time

  curr_time = pm_well%option%time - pm_well%option%flow_dt
  pm_well%dt_flow = pm_well%realization%option%flow_dt

  if (any(pm_well%well_grid%h_rank_id == pm_well%option%myrank)) then
    call PMWellInitializeTimestepFlow(pm_well,curr_time)
  endif

end subroutine PMWellInitializeTimestepWIPP

! ************************************************************************** !

subroutine PMWellInitializeTimestepWIPPQI(this)

  implicit none

  class(pm_well_wipp_qi_type) :: this

  call PMWellInitializeTimestepWIPP(this)

end subroutine PMWellInitializeTimestepWIPPQI

! ************************************************************************** !

subroutine PMWellInitializeTimestepFlow(pm_well,curr_time)
  !
  ! Initializes and takes the time step for the well process model - flow.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 05/11/2023

  use Option_module

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: curr_time
  PetscInt :: i

  ! update the reservoir object with current reservoir properties
  call PMWellCopyReservoir(pm_well%well%reservoir,pm_well%well%reservoir_save,&
                              pm_well%transport)

  call PMWellUpdateReservoirWIPP(pm_well,UNINITIALIZED_INTEGER)

  call PMWellComputeWellIndex(pm_well)

  call PMWellUpdateStrata(pm_well,curr_time)

  if (initialize_well_flow) then
      ! enter here if its the very first timestep
      call PMWellInitializeWellFlow(pm_well)
  endif

  do i = 1,pm_well%option%nflowdof
    call PMWellCopyWell(pm_well%well,pm_well%well_pert(i),pm_well%transport)
  enddo

  if (pm_well%well%bh_p_set_by_reservoir) then
      pm_well%well%bh_p = pm_well%well%reservoir%p_l(1)
  endif
  if (pm_well%well%bh_sg_set_by_reservoir) then
      pm_well%well%bh_sg = pm_well%well%reservoir%s_g(1)
  endif

  call pm_well%UpdateFlowProperties(PETSC_FALSE,ZERO_INTEGER)

end subroutine PMWellInitializeTimestepFlow

! ************************************************************************** !

subroutine PMWellInitializeWellTran(pm_well)
  !
  ! Initializes the well for the first time step for transport.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscInt :: k

  ! set initial transport parameters to the reservoir transport parameters
  if (pm_well%transport) then
    if (Initialized(pm_well%intrusion_time_start)) then
      ! set the borehole concentrations to the borehole zero value now
      do k = 1,pm_well%well_grid%nsegments
        pm_well%well%aqueous_mass(:,k) = &
          pm_well%bh_zero_value * pm_well%well%volume(k)
        pm_well%well%aqueous_conc(:,k) = &
          pm_well%well%aqueous_mass(:,k) / &                           ! [mol]
          (pm_well%well%phi(k)*pm_well%well%volume(k)*pm_well%well% &
          liq%s(k)) ! [m3-liq]
      enddo
    else
      ! set the wellbore concentrations to the reservoir values
      pm_well%well%aqueous_conc = pm_well%well%reservoir%aqueous_conc
      do k = 1,pm_well%well_grid%nsegments
        pm_well%well%aqueous_mass(:,k) = pm_well%well%aqueous_conc(:,k) * &
                pm_well%well%phi(k) * pm_well%well%volume(k) * &
                pm_well%well%liq%s(k)
      enddo
    endif
    pm_well%tran_soln%prev_soln%aqueous_conc = pm_well%well%aqueous_conc
    pm_well%tran_soln%prev_soln%aqueous_mass = pm_well%well%aqueous_mass
    pm_well%tran_soln%prev_soln%resr_aqueous_conc = &
                                     pm_well%well%reservoir%aqueous_conc
  endif

  initialize_well_tran = PETSC_FALSE

end subroutine PMWellInitializeWellTran

! ************************************************************************** !

subroutine PMWellUpdateReservoirWIPP(pm_well,wippflo_update_index)
  !
  ! Updates the reservoir properties for the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

  use WIPP_Flow_Aux_module
  use Material_Aux_module
  use NW_Transport_Aux_module
  use Grid_module

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: wippflo_update_index

  type(well_reservoir_type), pointer :: reservoir
  type(wippflo_auxvar_type), pointer :: wippflo_auxvar
  type(nw_transport_auxvar_type), pointer :: nwt_auxvar
  type(material_auxvar_type), pointer :: material_auxvar
  type(grid_type), pointer :: res_grid
  type(option_type), pointer :: option
  type(well_comm_type), pointer :: well_comm
  PetscInt :: k, indx, vec_size
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  option => pm_well%option
  well_comm => pm_well%well_comm
  reservoir => pm_well%well%reservoir

  res_grid => pm_well%realization%patch%grid

  reservoir%tmp_flow = -MAX_DOUBLE
  reservoir%tmp_tran = -MAX_DOUBLE

  if (wippflo_update_index < ZERO_INTEGER) then
    indx = ZERO_INTEGER
  else
    indx = wippflo_update_index
  endif

  do k = 1,pm_well%well_grid%nsegments
    if (pm_well%well_grid%h_rank_id(k) /= option%myrank) cycle

    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    wippflo_auxvar => &
      pm_well%realization%patch%aux%wippflo%auxvars(indx,ghosted_id)
    if (pm_well%transport) then
      nwt_auxvar => &
        pm_well%realization%patch%aux%nwt%auxvars(ghosted_id)
    endif
    material_auxvar => &
      pm_well%realization%patch%aux%material%auxvars(ghosted_id)

    reservoir%tmp_flow(k,1) = wippflo_auxvar%pres(option%liquid_phase) ! p_l
    reservoir%tmp_flow(k,2) = wippflo_auxvar%pres(option%gas_phase) ! p_g
    reservoir%tmp_flow(k,3) = wippflo_auxvar%sat(option%liquid_phase) ! s_l
    reservoir%tmp_flow(k,4) = wippflo_auxvar%sat(option%gas_phase) ! s_g
    reservoir%tmp_flow(k,5) = wippflo_auxvar%mobility(option%liquid_phase) ! mob_l
    reservoir%tmp_flow(k,6) = wippflo_auxvar%mobility(option%gas_phase) ! mob_g
    reservoir%tmp_flow(k,7) = wippflo_auxvar%kr(option%liquid_phase) ! kr_l
    reservoir%tmp_flow(k,8) = wippflo_auxvar%kr(option%gas_phase) ! kr_g
    reservoir%tmp_flow(k,9) = wippflo_auxvar%den_kg(option%liquid_phase) ! den_l
    reservoir%tmp_flow(k,10) = wippflo_auxvar%den_kg(option%gas_phase) ! den_g
    reservoir%tmp_flow(k,11) = wippflo_auxvar%mu(option%liquid_phase) ! vis_l
    reservoir%tmp_flow(k,12) = wippflo_auxvar%mu(option%gas_phase) ! vis_g
    reservoir%tmp_flow(k,13) = wippflo_auxvar%effective_porosity ! e_por

    reservoir%tmp_flow(k,14) = material_auxvar%permeability(1) ! kx
    reservoir%tmp_flow(k,15) = material_auxvar%permeability(2) ! ky
    reservoir%tmp_flow(k,16) = material_auxvar%permeability(3) ! kz
    reservoir%tmp_flow(k,17) = material_auxvar%volume

    if (res_grid%itype == STRUCTURED_GRID) then
      reservoir%tmp_flow(k,18) = res_grid%structured_grid%dx(ghosted_id) ! dx
      reservoir%tmp_flow(k,19) = res_grid%structured_grid%dy(ghosted_id) ! dy
      reservoir%tmp_flow(k,20) = res_grid%structured_grid%dz(ghosted_id) ! dz
    else
      reservoir%tmp_flow(k,20) = pm_well%well_grid%res_dz(k) ! dz
      reservoir%tmp_flow(k,18) = sqrt(material_auxvar%volume/ &
                                  reservoir%tmp_flow(k,20)) ! dx
      reservoir%tmp_flow(k,19) = reservoir%tmp_flow(k,18) ! dy
    endif

    if (pm_well%transport) then
      !aqueous_conc
      reservoir%tmp_tran(:,k,1) = nwt_auxvar%aqueous_eq_conc(:)
      !aqueous_mass = aq_conc * e_por * volume * s_l
      reservoir%tmp_tran(:,k,2) = &
            reservoir%tmp_tran(:,k,1) * reservoir%tmp_flow(k,13) * &
            reservoir%tmp_flow(k,17) * reservoir%tmp_flow(k,3)
    endif
  enddo

  vec_size = pm_well%well_grid%nsegments
  if (well_comm%commsize > 1) then
    ! Updates reservoir property vector in place using the maximum value.
    ! The rank-updated value will be larger then the initialized value (-MAX_DOUBLE)
    call MPI_Allreduce(MPI_IN_PLACE,reservoir%tmp_flow,vec_size*20,&
                       MPI_DOUBLE_PRECISION,MPI_MAX,pm_well%well_comm%comm,ierr)
    if (pm_well%transport) then
      call MPI_Allreduce(MPI_IN_PLACE,reservoir%tmp_tran,&
                         vec_size*pm_well%nspecies*2,&
                         MPI_DOUBLE_PRECISION,MPI_MAX,pm_well%well_comm%comm,ierr)
    endif
  endif

  do k = 1,pm_well%well_grid%nsegments
    reservoir%p_l(k) = reservoir%tmp_flow(k,1)
    reservoir%p_g(k) = reservoir%tmp_flow(k,2)
    reservoir%s_l(k) = reservoir%tmp_flow(k,3)
    reservoir%s_g(k) = reservoir%tmp_flow(k,4)
    reservoir%mobility_l(k) = reservoir%tmp_flow(k,5)
    reservoir%mobility_g(k) = reservoir%tmp_flow(k,6)
    reservoir%kr_l(k) = reservoir%tmp_flow(k,7)
    reservoir%kr_g(k) = reservoir%tmp_flow(k,8)
    reservoir%den_l(k) = reservoir%tmp_flow(k,9)
    reservoir%den_g(k) = reservoir%tmp_flow(k,10)
    reservoir%visc_l(k) = reservoir%tmp_flow(k,11)
    reservoir%visc_g(k) = reservoir%tmp_flow(k,12)
    reservoir%e_por(k) = reservoir%tmp_flow(k,13)
    reservoir%kx(k) = reservoir%tmp_flow(k,14)
    reservoir%ky(k) = reservoir%tmp_flow(k,15)
    reservoir%kz(k) = reservoir%tmp_flow(k,16)
    reservoir%volume(k) = reservoir%tmp_flow(k,17)
    reservoir%dx(k) = reservoir%tmp_flow(k,18)
    reservoir%dy(k) = reservoir%tmp_flow(k,19)
    reservoir%dz(k) = reservoir%tmp_flow(k,20)
  enddo

  if (pm_well%transport) then
    do k = 1,pm_well%well_grid%nsegments
      reservoir%aqueous_conc(:,k) = reservoir%tmp_tran(:,k,1)
      reservoir%aqueous_mass(:,k) = reservoir%tmp_tran(:,k,2)
    enddo
  endif

end subroutine PMWellUpdateReservoirWIPP

! ************************************************************************** !

subroutine PMWellUpdateFlowRatesWIPPQI(this,well_pert,res_pert, &
                                       segment_index,ierr)
  !
  ! This subroutine performs the well rate computation when called from the
  ! fully- or quasi-coupled source/sink update.
  !
  ! Author: Michael Nole
  ! Date: 01/16/2023
  !

  Use Option_module

  implicit none

  class(pm_well_wipp_qi_type) :: this
  PetscErrorCode :: ierr
  PetscInt :: well_pert
  PetscInt :: res_pert
  PetscInt :: segment_index

  type(option_type), pointer :: option
  PetscReal :: time
  PetscInt :: k


  option => this%option
  time = this%realization%option%time
  !  if (this%flow_soln%n_steps < 1) return

  ! Need to limit well model timestepping
  this%min_dt_flow = option%flow_dt * min_flow_dt_scale

  if (.not. this%well_on .and. Initialized(this%intrusion_time_start) .and. &
      time < this%intrusion_time_start) then
      this%srcsink_water(well_pert,:) = 0.d0
      this%srcsink_gas(well_pert,:) = 0.d0
    return
  elseif (.not. this%well_on) then
    this%well_on = PETSC_TRUE
  endif

  this%print_output = PETSC_FALSE

  call PMWellUpdateReservoirWIPP(this,res_pert)

  if (initialize_well_flow) then
    call PMWellInitializeWellFlow(this)
    this%flow_soln%prev_soln%pl = this%well%pl
    this%flow_soln%prev_soln%sg = this%well%gas%s
    this%flow_soln%soln_save%pl = this%well%pl
    this%flow_soln%soln_save%sg = this%well%gas%s
    this%flow_soln%prev_soln%bh_p = this%well%bh_p
    this%flow_soln%soln_save%bh_p = this%well%bh_p
    do k = 1,option%nflowdof
      call PMWellCopyWell(this%well,this%well_pert(k),this%transport)
    enddo
  else
    this%well%pl = this%flow_soln%soln_save%pl
    this%well%gas%s = this%flow_soln%soln_save%sg
    this%well%bh_p = this%flow_soln%soln_save%bh_p
  endif

  ! Quasi-implicit
  do k = 1,option%nflowdof
    call PMWellCopyWell(this%well,this%well_pert(k),this%transport)
  enddo

  call this%UpdateFlowProperties(PETSC_FALSE,ZERO_INTEGER)

  this%dt_flow = this%realization%option%flow_dt
  call this%SolveFlow(well_pert,ierr)
  this%print_output = PETSC_TRUE

  this%srcsink_water(well_pert,:) = -1.d0 * &
                                          this%well%liq%Q(:)! [kmol/s]
  this%srcsink_gas(well_pert,:)   = -1.d0 * &
                                          this%well%gas%Q(:)! [kmol/s]


end subroutine PMWellUpdateFlowRatesWIPPQI

! ************************************************************************** !

subroutine PMWellModifyFlowResWIPPQI(this,residual)
  !
  ! This subroutine computes the well contribution to the reservoir residual
  ! when called from the quasi-coupled source/sink update in WIPP
  ! FLOW mode.
  !
  ! Author: Michael Nole
  ! Date: 01/16/2023
  !

  implicit none

  class(pm_well_wipp_qi_type) :: this
  PetscReal, pointer :: residual(:)

  PetscReal :: Res(this%flow_soln%ndof)
  PetscInt :: air_comp_id, wat_comp_id
  PetscInt :: local_id, local_start, local_end
  PetscInt :: ghosted_id
  PetscInt :: k

  air_comp_id = this%option%air_id
  wat_comp_id = this%option%water_id

  Res(:) = 0.d0

  if (this%well_comm%comm == MPI_COMM_NULL) return
  do k = 1,this%well_grid%nsegments
    if (this%well_grid%h_rank_id(k) /= this%option%myrank) cycle

    ghosted_id = this%well_grid%h_ghosted_id(k)
    local_id = this%realization%patch%grid%nG2L(ghosted_id)
    local_end = local_id * this%flow_soln%ndof
    local_start = local_end - this%flow_soln%ndof + 1

    ! kmol/sec
    Res(wat_comp_id) = this%srcsink_water(UNPERT,k)
    Res(air_comp_id) = this%srcsink_gas(UNPERT,k)

    call WIPPFloConvertUnitsToBRAGFlo(Res,this%realization%patch% &
                                      aux%Material% &
                                      auxvars(ghosted_id),this%option)
    residual(local_start:local_end) = residual(local_start:local_end) - &
                                      Res(:)

  enddo

end subroutine PMWellModifyFlowResWIPPQI

! ************************************************************************** !

subroutine PMWellModifyFlowJacWIPPQI(this,Jac,ierr)
  !
  ! This subroutine computes the well contribution to the reservoir
  ! Jacobian when called from the fully- or quasi-coupled source/sink update.
  !
  ! Author: Michael Nole
  ! Date: 01/16/2023
  !

  use Option_module

  implicit none

  class(pm_well_wipp_qi_type) :: this
  Mat :: Jac
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  PetscReal :: pert_pl, pert_sg
  PetscReal :: res, res_pert_pl, res_pert_sg
  PetscInt :: ghosted_id
  PetscReal, allocatable :: J_block(:,:)
  PetscInt :: air_comp_id, wat_comp_id
  PetscInt :: k

  option => this%option

  air_comp_id = this%option%air_id
  wat_comp_id = this%option%water_id

  if (this%well_comm%comm == MPI_COMM_NULL) return
  allocate(J_block(this%flow_soln%ndof,this%flow_soln%ndof))
  J_block = 0.d0

  do k = 1,this%well_grid%nsegments
    if (this%well_grid%h_rank_id(k) /= this%option%myrank) cycle

    J_block = 0.d0

    ghosted_id = this%well_grid%h_ghosted_id(k)
    pert_pl = this%realization%patch%aux%WIPPFlo% &
                    auxvars(WIPPFLO_LIQUID_PRESSURE_DOF,ghosted_id)%pert
    pert_sg = this%realization%patch%aux%WIPPFlo% &
                    auxvars(WIPPFLO_GAS_SATURATION_DOF,ghosted_id)%pert

    ! Liquid portion
    res = this%srcsink_water(UNPERT,k)
    res_pert_pl = this%srcsink_water(PERT_WRT_PL,k)
    res_pert_sg = this%srcsink_water(PERT_WRT_SG,k)
    J_block(wat_comp_id,WIPPFLO_LIQUID_PRESSURE_DOF) = &
            (res_pert_pl - res)/pert_pl
    J_block(wat_comp_id,WIPPFLO_GAS_SATURATION_DOF) = &
            (res_pert_sg - res)/pert_sg

    ! Gas portion
    res = this%srcsink_gas(UNPERT,k)
    res_pert_pl = this%srcsink_gas(PERT_WRT_PL,k)
    res_pert_sg = this%srcsink_gas(PERT_WRT_SG,k)
    J_block(air_comp_id,WIPPFLO_LIQUID_PRESSURE_DOF) = &
            (res_pert_pl - res)/pert_pl
    J_block(air_comp_id,WIPPFLO_GAS_SATURATION_DOF) = &
            (res_pert_sg - res)/pert_sg

    J_block = -1.d0*J_block
    call WIPPFloConvertUnitsToBRAGFlo(J_block,this%realization%patch%aux% &
                                  Material%auxvars(ghosted_id),this%option)

    call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1, &
                                  J_block,ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo

  if (allocated(J_block)) deallocate(J_block)

end subroutine PMWellModifyFlowJacWIPPQI

! ************************************************************************** !

subroutine PMWellSolveWIPP(pm_well,time,qi_coupling,ierr)

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: time
  PetscBool :: qi_coupling
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: curr_time, curr_time_converted

  curr_time = pm_well%option%time + pm_well%cumulative_dt_tran
  curr_time_converted = curr_time/pm_well%output_option%tconv

  ierr = 0 ! If pm_well is not set to zero, TS_STOP_FAILURE occurs if the solve
           ! routines are not entered, either due to an inactive well or due
           ! to being on a process that doesn't contain a well segment.


  if (Initialized(pm_well%intrusion_time_start) .and. &
      (curr_time < pm_well%intrusion_time_start)) then
    write(out_string,'(" Inactive.    Time =",1pe12.5," ",a4)') &
          curr_time_converted,pm_well%output_option%tunit
    call PrintMsg(pm_well%option,out_string)
    return
  endif

  if (qi_coupling) then
    write(out_string,'(" FLOW Step          Quasi-implicit wellbore flow &
                      &coupling is being used.")')
    call PrintMsg(pm_well%option,out_string)
    qi_coupling = PETSC_FALSE
  else
    call pm_well%SolveFlow(UNINITIALIZED_INTEGER,ierr)
  endif

  call PMWellCalcCumulativeQFlux(pm_well)

  !Debugging
  !call MPI_Barrier(pm_well%option%comm%communicator,ierr);CHKERRQ(ierr)
  if (pm_well%transport) then
    write(out_string,'(" TRAN Step          Quasi-implicit wellbore &
                     &transport coupling is being used.")')
    call PrintMsg(pm_well%option,out_string)

    ! must call prior to updating the prev_soln vectors
    call PMWellCalcCumulativeTranFlux(pm_well)

    pm_well%tran_soln%prev_soln%aqueous_conc = pm_well%well%aqueous_conc
    pm_well%tran_soln%prev_soln%aqueous_mass = pm_well%well%aqueous_mass
    pm_well%tran_soln%prev_soln%resr_aqueous_conc = &
                                  pm_well%well%reservoir%aqueous_conc
  endif

end subroutine PMWellSolveWIPP

! ************************************************************************** !

subroutine PMWellSolveWIPPQI(this,time,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  implicit none

  class(pm_well_wipp_qi_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  call PMWellSolveWIPP(this,time,this%update_for_flow_qi_coupling,ierr)

  if (this%transport .and. this%dt_tran < this%dt_flow) then
    this%update_for_flow_qi_coupling = PETSC_TRUE
  endif

end subroutine PMWellSolveWIPPQI

! ************************************************************************** !

subroutine WIPPWellSolveFlowSequential(this,perturbation_index,ierr)
  !
  ! Author: Michael Nole
  ! Date: 12/01/2021
  !

  use Option_module
  use Grid_module
  use EOS_Water_module
  use EOS_Gas_module
  use SCO2_Aux_module, only : fmw_comp, SCO2BrineDensity

  implicit none

  class(pm_well_sequential_type) :: this
  PetscInt :: perturbation_index
  PetscErrorCode :: ierr

  type(grid_type), pointer :: reservoir_grid
  type(option_type), pointer :: option
  type(well_soln_flow_type), pointer :: flow_soln
  character(len=MAXSTRINGLENGTH) :: out_string
  PetscInt :: n_iter,ts_cut,easy_converge_count
  PetscInt :: istart, iend
  PetscReal :: res(this%flow_soln%ndof)
  PetscReal :: res_fixed(this%flow_soln%ndof*this%well_grid%nsegments)
  PetscBool :: at_steady_state
  PetscReal :: ss_check_p(this%well_grid%nsegments,2), &
               ss_check_s(this%well_grid%nsegments,2)
  PetscInt :: i
  PetscInt :: ss_step_count, steps_to_declare_ss
  PetscReal :: gravity
  PetscReal, parameter :: threshold_p = 0.d0
  PetscReal, parameter :: epsilon = 1.d-14

  option => this%realization%option
  reservoir_grid => this%realization%patch%grid

  if (this%well_comm%comm == MPI_COMM_NULL) return

  flow_soln => this%flow_soln

  ierr = 0

  ts_cut = 0
  easy_converge_count = 0

  this%cumulative_dt_flow = 0.d0
  flow_soln%converged = PETSC_FALSE
  flow_soln%not_converged = PETSC_TRUE

  ss_check_p(:,1) = this%well%pl(:)
  ss_check_s(:,1) = this%well%gas%s(:)
  at_steady_state = PETSC_FALSE
  ss_step_count = 0
  steps_to_declare_ss = 10

  gravity = option%gravity(Z_DIRECTION)

  ! update well index
  call PMWellComputeWellIndex(this)

  do while (this%cumulative_dt_flow < this%realization%option%flow_dt)

    ! update the well src/sink Q vector at start of time step
    call PMWellUpdateWellQ(this%well,this%well%reservoir)

    call PMWellPreSolveFlow(this)

    ! Fixed accumulation term
    res_fixed = 0.d0
    res = 0.d0
    do i = 1,this%well_grid%nsegments
      call PMWellAccumulationFlow(this,this%well,i,res)
      istart = flow_soln%ndof*(i-1)+1
      iend = flow_soln%ndof*i
      res_fixed(istart:iend) = -1.d0 * res * this%dt_flow
    enddo

    n_iter = 0

    do while (flow_soln%not_converged)

      if (n_iter > (flow_soln%max_iter-1)) then
        flow_soln%cut_timestep = PETSC_TRUE
        if (this%print_output) then
          out_string = ' Maximum number of FLOW Newton iterations reached. &
                        &Cutting timestep!'
          call PrintMsg(this%option,out_string)
        endif
        call PMWellCutTimestepFlow(this)
        n_iter = 0
        ts_cut = ts_cut + 1
        easy_converge_count = 0

        if (ss_step_count > 2 .and. this%ss_check) then
          at_steady_state = PETSC_TRUE
          this%cumulative_dt_flow = this%realization%option%flow_dt
          WRITE(out_string,'(" PM Well FLOW convergence declared due to &
            &automatic time step control criterion. ")')
          call PrintMsg(this%option,out_string)
        endif


        exit
      endif
      if (ts_cut > flow_soln%max_ts_cut) then
        this%realization%option%io_buffer = &
          ' Maximum timestep cuts reached in PM Well FLOW. Solution has not &
           &converged. Exiting.'
        if (this%print_well) then
          call PMWellOutputSequential(this)
        endif
        call PrintErrMsg(this%realization%option)
      endif

      if (this%dt_flow <= this%min_dt_flow) then
        this%well_force_ts_cut = 1
        call PMWellCopyReservoir(this%well%reservoir_save, &
                                 this%well%reservoir,&
                                 this%transport)
        return
      endif

      flow_soln%residual = 0.d0
      flow_soln%residual = res_fixed / this%dt_flow

      easy_converge_count = easy_converge_count + 1

      call PMWellNewtonFlow(this)

      if (this%well_force_ts_cut > 0) return

      call PMWellCheckConvergenceFlow(this,n_iter,res_fixed)

    enddo

    if (easy_converge_count > 4 ) then
      if (this%cumulative_dt_flow + this%dt_flow * &
          flow_soln%ts_cut_factor < this%realization%option%flow_dt) then
        this%dt_flow = this%dt_flow * flow_soln%ts_cut_factor
        flow_soln%cut_timestep = PETSC_FALSE
      endif
    endif

    if (this%cumulative_dt_flow + this%dt_flow > &
        this%realization%option%flow_dt) then
      this%dt_flow = this%realization%option%flow_dt - &
                        this%cumulative_dt_flow
    endif

    ts_cut = 0
    flow_soln%n_steps = flow_soln%n_steps + 1

    ! Check if we're at steady-state

    ! TOUGH way:
    if (n_iter == 1 .and. flow_soln%converged) then
      ss_step_count = ss_step_count+1
    else
      ss_step_count = 0
    endif

    if (this%ss_check) then
      if (ss_step_count >= steps_to_declare_ss) then
        at_steady_state = PETSC_TRUE
        this%cumulative_dt_flow = this%realization%option%flow_dt
      endif
    endif

    ! Other way:
    !if (this%ss_check .and. flow_soln%converged) then
    !  ss_check_p(:,2) = this%well%pl(:)
    !  ss_check_s(:,2) = this%well%gas%s(:)

    !  dpdt = (ss_check_p(:,2) - ss_check_p(:,1)) / this%dt_flow
    !  dsdt = (ss_check_s(:,2) - ss_check_s(:,1)) / this%dt_flow

    !  if (maxval(abs(dpdt)) < eps_p) then
    !    if (maxval(abs(dsdt)) < eps_s) then
    !      ss_step_count = ss_step_count + 1
    !      if (ss_step_count > steps_to_declare_ss) at_steady_state=PETSC_TRUE
    !    else
    !      ss_step_count = 0
    !    endif
    !  endif
    !  if (at_steady_state) then
    !    this%cumulative_dt_flow = this%realization%option%flow_dt
    !  endif
    !  ss_check_p(:,1) = this%well%pl(:)
    !  ss_check_s(:,1) = this%well%gas%s(:)
    !endif

    call PMWellPostSolveFlow(this)

  enddo

end subroutine WIPPWellSolveFlowSequential

! ************************************************************************** !

subroutine PMWellPostSolveFlow(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: cur_time_converted

  cur_time_converted = pm_well%option%time/pm_well%output_option%tconv

  WRITE(out_string,'(" PM Well FLOW Step Complete!    Time=",1pe12.5," &
                    &",a4,"Total Newton Its =",i8)') &
                    cur_time_converted,pm_well%output_option%tunit, &
                    pm_well%flow_soln%n_newton
  call PrintMsg(pm_well%option,out_string)

end subroutine PMWellPostSolveFlow

! ************************************************************************** !

subroutine PMWellPostSolveTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/22/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscReal :: cur_time, cur_time_converted

  cur_time = pm_well%option%time + pm_well%option%tran_dt
  pm_well%tran_soln%tran_time = cur_time
  cur_time_converted = cur_time/pm_well%output_option%tconv

end subroutine PMWellPostSolveTran

! ************************************************************************** !

subroutine PMWellSolveFlowWIPPQI(this,perturbation_index,ierr)
  !
  ! Author: Michael Nole
  ! Date: 12/01/2021
  !

  implicit none

  class(pm_well_wipp_qi_type) :: this
  PetscInt :: perturbation_index
  PetscErrorCode :: ierr

  call WIPPWellSolveFlowSequential(this,perturbation_index,ierr)

end subroutine PMWellSolveFlowWIPPQI

! ************************************************************************** !

subroutine UpdateFlowPropertiesWIPP(pm_well,pert,index)
  !
  ! Updates flow well object properties, when WIPP_FLOW is the flow mode.
  !
  ! Author: Michael Nole
  ! Date: 01/06/2022
  !

  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_WIPP_module

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscBool :: pert
  PetscInt :: index

  type(well_type), pointer :: well
  type(characteristic_curves_ptr_type), pointer ::characteristic_curves_array(:)
  type(option_type), pointer :: option
  class(characteristic_curves_type), pointer :: characteristic_curves
  class(sat_func_base_type), pointer :: saturation_function
  type(strata_type), pointer :: strata
  PetscInt :: i,nsegments
  PetscReal :: T,dw,dg,dwmol,dwp,dwt,Psat,visl,visg
  PetscReal :: Pc,dpc_dsatl,krl,dkrl_dsatl,krg,dkrg_dsatl
  PetscErrorCode :: ierr

  if (pert) then
    well => pm_well%well_pert(index)
  else
    well => pm_well%well
  endif
  characteristic_curves_array => pm_well%realization%patch%characteristic_curves_array
  option => pm_well%realization%option

  T = option%flow%reference_temperature

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  nsegments =pm_well%well_grid%nsegments

  do i = 1,nsegments
    ! Material Properties
    strata => pm_well%strata_list%first
    do
      if (.not.associated(strata)) exit
      if (strata%id == pm_well%well_grid%strata_id(i)) then
        well%ccid(i) = strata%material_property%saturation_function_id
        well%permeability(i) = strata%material_property%permeability(3,3)
        well%phi(i) = strata%material_property%porosity
        exit
      endif
      strata => strata%next
    enddo

    ! Saturations
    well%gas%s(i) = max(well%gas%s(i),0.d0)
    well%gas%s(i) = min(well%gas%s(i),1.d0)
    well%liq%s(i) = 1.d0 - well%gas%s(i)

    ! Capillary Pressure
    characteristic_curves => characteristic_curves_array(well%ccid(i))%ptr
    saturation_function => characteristic_curves%saturation_function
    select type(sat_func => saturation_function)
      class is (sat_func_krp3_type)
        if (.not. option%flow%pct_updated) then
          sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
                   CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
        else
          call sat_func% &
                   CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
        endif
      class is (sat_func_krp4_type)
        if (.not. option%flow%pct_updated) then
          sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
               CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
        else
          call sat_func% &
               CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
        endif
      class is (sat_func_krp5_type)
            if (.not. option%flow%pct_updated) then
              sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
              option%flow%pct_updated = PETSC_TRUE
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            else
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            endif
      class default
        call sat_func% &
               CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
    end select
    well%pg(i) = well%pl(i) + Pc

    ! Relative Permeabilities
    call characteristic_curves%liq_rel_perm_function% &
               RelativePermeability(well%liq%s(i),krl,dkrl_dsatl,option)
    well%liq%kr(i) = krl

    call characteristic_curves%gas_rel_perm_function% &
               RelativePermeability(well%liq%s(i),krg,dkrg_dsatl,option)
    well%gas%kr(i) = krg

    !Density
    call EOSWaterDensityBRAGFLO(T,well%pl(i),PETSC_FALSE, &
                                dw,dwmol,dwp,dwt,ierr)
    call EOSGasDensity(T,well%pg(i),dg,ierr)

    well%liq%den(i) = dw
    !No water vapor in WIPP_Darcy mode
    well%gas%den(i) = dg

    !Viscosity
    call EOSWaterSaturationPressure(T,Psat,ierr)
    call EOSWaterViscosity(T,well%pl(i),Psat,visl,ierr)
    call EOSGasViscosity(T,well%pg(i),well%pg(i),dg,visg,ierr)

    well%liq%visc(i) = visl
    well%gas%visc(i) = visg

  enddo
end subroutine UpdateFlowPropertiesWIPP

! ************************************************************************** !

subroutine UpdateFlowPropertiesWIPPQI(this,pert,index)
  !
  ! Updates flow well object properties, when WIPP_FLOW is the flow mode.
  !
  ! Author: Michael Nole
  ! Date: 01/06/2022
  !

  implicit none

  class(pm_well_wipp_qi_type) :: this
  PetscBool :: pert
  PetscInt :: index

  call UpdateFlowPropertiesWIPP(this,pert,index)

end subroutine UpdateFlowPropertiesWIPPQI

! ************************************************************************** !

subroutine PMWellUpdatePropertiesTran(pm_well)
  !
  ! Updates transport related well object properties for each time step.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/13/2022
  !

  use Transport_Constraint_Base_module
  use Transport_Constraint_NWT_module
  use Transport_Constraint_module
  use Condition_module

  implicit none

  class(pm_well_type) :: pm_well

  type(tran_condition_type), pointer :: tran_condition
  class(tran_constraint_base_type), pointer :: cur_constraint

  ! update the top of hole boundary condition with current constraint
  tran_condition => pm_well%realization%transport_conditions%first
  do
    if (.not.associated(tran_condition)) exit
      if (trim(tran_condition%name) == &
          trim(pm_well%well%tran_condition_name)) exit
    tran_condition => tran_condition%next
  enddo
  cur_constraint => tran_condition%cur_constraint_coupler%constraint
  select type(constraint=>cur_constraint)
    class is (tran_constraint_nwt_type)
      if (any(constraint%nwt_species%constraint_type /=  &
          CONSTRAINT_AQ_EQUILIBRIUM)) then
        pm_well%option%io_buffer = 'TRANSPORT_CONDITION ' // &
          trim(pm_well%well%tran_condition_name) // ' for WELLBORE_MODEL,&
          &WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE CONSTRAINT must be of &
          &type "AQ".'
        call PrintErrMsg(pm_well%option)
      endif
      pm_well%well%aqueous_conc_th = constraint%nwt_species%constraint_conc
  end select

end subroutine PMWellUpdatePropertiesTran

! ************************************************************************** !

subroutine PMWellNewtonFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 01/20/2022
  !

  use Utility_module

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscReal :: identity(pm_well%nphase*pm_well%well_grid%nsegments,&
                        pm_well%nphase*pm_well%well_grid%nsegments)
  PetscReal :: new_dx(pm_well%nphase*pm_well%well_grid%nsegments)
  PetscInt :: indx(pm_well%nphase*pm_well%well_grid%nsegments)
  PetscInt :: i,j
  PetscInt :: d

  call PMWellUpdateWellQ(pm_well%well,pm_well%well%reservoir)

  call PMWellPerturb(pm_well)

  call PMWellResidualFlow(pm_well)

  call PMWellJacobianFlow(pm_well)

  do i = 1,pm_well%nphase*pm_well%well_grid%nsegments
    do j = 1,pm_well%nphase*pm_well%well_grid%nsegments
      if (i==j) then
        identity(i,j) = 1.d0
      else
        identity(i,j) = 0.d0
      endif
    enddo
  enddo
  call LUDecomposition(pm_well%flow_soln%Jacobian,pm_well%nphase*pm_well% &
                        well_grid%nsegments,indx,d)
  call LUBackSubstitution(pm_well%flow_soln%Jacobian, &
                          pm_well%nphase*pm_well%well_grid%nsegments,&
                          indx,pm_well%flow_soln%residual)
  new_dx = -1.d0 * pm_well%flow_soln%residual


  do i = 1,pm_well%well_grid%nsegments
    if (dabs(new_dx(i)) > 1.d15) then
      pm_well%well_force_ts_cut = 1
      return
    endif
    if (isnan(new_dx(i))) then
      pm_well%well_force_ts_cut = 1
      return
    endif
  enddo

  pm_well%flow_soln%update = new_dx

  call PMWellUpdateSolutionFlow(pm_well)

end subroutine PMWellNewtonFlow

! ************************************************************************** !

subroutine PMWellNewtonTran(pm_well,n_iter)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  use Utility_module

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: n_iter

  PetscInt :: nm, dummy
  PetscInt :: indx(pm_well%nspecies*pm_well%well_grid%nsegments)

  if (.not. any(pm_well%option%myrank == pm_well%well_grid%h_rank_id)) return

  nm = pm_well%nspecies * pm_well%well_grid%nsegments

  ! at this time, the tran_soln%residual vector has been zero'd out and has
  ! been loaded with the fixed accumulation divided by current dt

  call PMWellResidualTran(pm_well)

  call PMWellJacobianTran(pm_well)

  ! J dx = -R     => dx = J^(-1)(-R)
  ! [m3-bulk/sec] dx = -[mol/sec]
  ! dx in [mol/m3-bulk]

  call LUDecomposition(pm_well%tran_soln%Jacobian,nm,indx,dummy)

  call LUBackSubstitution(pm_well%tran_soln%Jacobian,nm,indx, &
                          pm_well%tran_soln%residual)

  pm_well%tran_soln%update = +1.0d0 * pm_well%tran_soln%residual ! [mol/m3-bulk]

  call PMWellUpdateSolutionTran(pm_well)

end subroutine PMWellNewtonTran

! ************************************************************************** !

subroutine PMWellUpdateWellQ(well,reservoir)
  !
  ! Updates the src/sink vector for the fluid object.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021
  !

  implicit none

  type(well_type) :: well
  type(well_reservoir_type), pointer :: reservoir

  type(well_fluid_type), pointer :: liq
  type(well_fluid_type), pointer :: gas

  PetscReal, parameter :: threshold_p = 0.d0 !1.d-2 !1.d-1
  PetscReal :: mobility, den_ave
  PetscBool :: upwind
  PetscInt :: i, nsegments

  liq => well%liq
  gas => well%gas

  nsegments = size(well%liq%Q)

  ! + Q goes out of well to reservoir
  ! - Q goes into well from reservoir

  do i = 1,nsegments
    if (dabs((reservoir%p_l(i)-well%pl(i)))/well%pl(i) > threshold_p) then
      upwind = reservoir%p_l(i) > well%pl(i)
      if (upwind) then
        mobility = reservoir%kr_l(i)/reservoir%visc_l(i)
      else
        mobility = liq%kr(i)/liq%visc(i)
      endif
      den_ave = 0.5d0 * (liq%den(i) + reservoir%den_l(i)) / FMWH2O
      ! Flowrate in kmol/s
      liq%Q(i) = den_ave*mobility*well%WI(i)* &
                  (reservoir%p_l(i)-well%pl(i))

      upwind = reservoir%p_g(i) > well%pg(i)
      if (upwind) then
        mobility = reservoir%kr_g(i)/reservoir%visc_g(i)
      else
        mobility = gas%kr(i)/gas%visc(i)
      endif
      den_ave = 0.5d0 * (gas%den(i) + reservoir%den_g(i)) / &
                          fmw_comp(TWO_INTEGER)

      ! Flowrate in kmol/s
      gas%Q(i) = den_ave*mobility*well%WI(i)* &
                  (reservoir%p_g(i)-well%pg(i))
    else
      liq%Q(i) = 0.d0
      gas%Q(i) = 0.d0
    endif
  enddo

end subroutine PMWellUpdateWellQ

! ************************************************************************** !

subroutine PMWellAccumulationFlow(pm_well,well,id,Res)
  !
  ! Computes the accumulation term for the flow residual based on the
  ! chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 12/23/2021
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  type(well_type) :: well
  PetscInt :: id
  PetscReal :: Res(pm_well%nphase)

  Res = 0.d0
  ! liquid accumulation term
  Res(1) = Res(1) + well%liq%s(id) * well%liq%den(id) / FMWH2O * &
            well%phi(id) * well%volume(id) / pm_well%dt_flow
  ! gas accumulation term
  Res(2) = Res(2) + well%gas%s(id) * well%gas%den(id) / &
            fmw_comp(TWO_INTEGER) * well%phi(id) * well%volume(id) / &
            pm_well%dt_flow

end subroutine PMWellAccumulationFlow

! ************************************************************************** !

subroutine PMWellSrcSink(pm_well,well,id,Res)
  !
  ! Computes the source/sink term for the flow residual based on the
  ! chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 02/26/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  type(well_type) :: well
  PetscInt :: id
  PetscReal :: Res(pm_well%nphase)

  Res = 0.d0

  call PMWellUpdateWellQ(pm_well%well,pm_well%well%reservoir)

  ! kmol/s
  Res(1) = Res(1) - well%liq%Q(id)
  ! kmol/s
  Res(2) = Res(2) - well%gas%Q(id)


end subroutine PMWellSrcSink

! ************************************************************************** !

subroutine PMWellAccumulationTran(pm_well,isegment,Res)
  !
  ! Computes the fixed accumulation term for the transport residual.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: isegment
  PetscReal :: Res(pm_well%nspecies)

  PetscInt :: ispecies, k

  Res = 0.d0

  ! porosity in [m^3-void/m^3-bulk]
  ! saturation in [m^3-liq/m^3-void]
  ! volume in [m^3-bulk]
  ! aqueous conc in [mol-species/m^3-liq]
  ! Res(:) in [mol-species]

  ! NOTE: division by dt occurs later, when the fixed accumulation is needed

  ! NOTE: this calculation uses the converged solution because it is using
  !       this%well%liq%s() and this%well%prev_soln%aqueous_conc() values

  do ispecies = 1,pm_well%nspecies
    k = ispecies
    Res(k) = pm_well%well%volume(isegment) * pm_well%well%phi(isegment) * &
             pm_well%well%liq%s(isegment) * &
             pm_well%tran_soln%prev_soln%aqueous_conc(ispecies,isegment)
  enddo

end subroutine PMWellAccumulationTran

! ************************************************************************** !

subroutine PMWellAccumDerivative(pm_well,local_id,Jac)
  !
  ! Computes the derivative of the accumulation term for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 01/05/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: local_id
  PetscReal :: Jac(pm_well%nphase,pm_well%nphase)

  PetscInt :: idof, irow
  PetscReal :: res(pm_well%nphase),res_pert(pm_well%nphase)

  call PMWellAccumulationFlow(pm_well,pm_well%well,local_id,res)

  do idof = 1, pm_well%nphase
    call PMWellAccumulationFlow(pm_well,pm_well%well_pert(idof),local_id, &
                                res_pert)
    do irow = 1, pm_well%nphase
      Jac(irow,idof) = (res_pert(irow)-res(irow))/pm_well%pert(local_id,idof)
    enddo !irow
  enddo ! idof

end subroutine PMWellAccumDerivative
! ************************************************************************** !

subroutine PMWellSrcSinkDerivative(pm_well,local_id,Jac)
  !
  ! Computes the derivative of the source/sink term for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 02/26/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: local_id
  PetscReal :: Jac(pm_well%nphase,pm_well%nphase)

  PetscInt :: idof, irow
  PetscReal :: res(pm_well%nphase),res_pert(pm_well%nphase)

  call PMWellSrcSink(pm_well,pm_well%well,local_id,res)

  do idof = 1, pm_well%nphase
    call PMWellSrcSink(pm_well,pm_well%well_pert(idof),local_id,res_pert)
    do irow = 1, pm_well%nphase
      Jac(irow,idof) = (res_pert(irow)-res(irow))/pm_well%pert(local_id,idof)
    enddo !irow
  enddo ! idof

end subroutine PMWellSrcSinkDerivative

! ************************************************************************** !

subroutine PMWellFlux(pm_well,well_up,well_dn,iup,idn,Res,save_flux)
  !
  ! Computes the internal flux terms for the residual based on
  ! the chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 12/23/2021
  !

  !MAN: will removing this break regression tests?
  ! use SCO2_Aux_module, only: sco2_fmw => fmw_comp

  implicit none

  class(pm_well_sequential_type) :: pm_well
  type(well_type) :: well_up, well_dn
  PetscInt :: iup, idn
  PetscReal :: Res(pm_well%nphase)
  PetscBool :: save_flux

  type(well_grid_type), pointer :: well_grid

  PetscReal :: perm_den_mu_area_ave_over_dist(2), perm_den_mu_area_up(2), &
               perm_den_mu_area_dn(2)
  PetscReal :: perm_up, perm_dn, dist_up, dist_dn, density_kg_ave, rel_perm
  PetscReal :: gravity_term, delta_pressure
  PetscReal :: v_darcy
  PetscReal :: density_ave_kmol, tot_mole_flux
  PetscReal :: up_scale, dn_scale
  PetscBool :: upwind
  PetscReal :: gravity


  well_grid => pm_well%well_grid

  gravity = pm_well%option%gravity(Z_DIRECTION)

  Res(:) = 0.d0

  ! Vertical well, no Klinkenberg

    perm_up = well_up%permeability(iup)
    perm_dn = well_dn%permeability(idn)
    dist_up = well_grid%dh(iup)/2.d0
    dist_dn = well_grid%dh(idn)/2.d0

    perm_den_mu_area_up(1) = perm_up * well_up%liq%den(iup) / &
                          FMWH2O / well_up%liq%visc(iup) * &
                          PI * (well_up%diameter(iup)/2.d0)**2
    perm_den_mu_area_up(2) = perm_up * well_up%gas%den(iup) / &
                          fmw_comp(TWO_INTEGER) / well_up%gas%visc(iup) * &
                          PI * (well_up%diameter(iup)/2.d0)**2
    perm_den_mu_area_dn(1) = perm_dn * well_dn%liq%den(idn) / &
                          FMWH2O / well_dn%liq%visc(idn) * &
                          PI * (well_dn%diameter(idn)/2.d0)**2
    perm_den_mu_area_dn(2) = perm_dn * well_dn%gas%den(idn) / &
                          fmw_comp(TWO_INTEGER) / well_dn%gas%visc(idn) * &
                          PI * (well_dn%diameter(idn)/2.d0)**2

    perm_den_mu_area_ave_over_dist(1) = &
            (perm_den_mu_area_up(1) * perm_den_mu_area_dn(1)) / &
            (dist_up*perm_den_mu_area_dn(1) + dist_dn*perm_den_mu_area_up(1))

    perm_den_mu_area_ave_over_dist(2) = &
            (perm_den_mu_area_up(2) * perm_den_mu_area_dn(2)) / &
            (dist_up*perm_den_mu_area_dn(2) + dist_dn*perm_den_mu_area_up(2))

    ! Liquid flux
    density_kg_ave = 0.5d0*(well_up%liq%den(iup)+well_dn%liq%den(idn))
    ! Assuming the well is always vertical and gravity is in the
    ! (-) direction
    gravity_term = density_kg_ave * gravity * &
                    0.5d0*(well_grid%dh(iup)+well_grid%dh(idn))
    delta_pressure = well_up%pl(iup) - well_dn%pl(idn) + &
                      gravity_term
    up_scale = 0.d0
    dn_scale = 0.d0

    upwind = delta_pressure > 0.d0

    ! Only upwinding the perm here, not mobility ratio
    ! following the LumpedHarmonic BRAGFLO formulation
    if (upwind) then
      up_scale = 1.d0
      rel_perm = well_up%liq%kr(iup)
    else
      dn_scale = 1.d0
      rel_perm = well_dn%liq%kr(idn)
    endif

    !kmol/sec
    tot_mole_flux = perm_den_mu_area_ave_over_dist(1) * rel_perm * &
              delta_pressure
    density_ave_kmol = density_kg_ave / fmw_comp(ONE_INTEGER)
    ! v_darcy = kmol/sec / kmol/m^3 / area[m^2]
    v_darcy = tot_mole_flux/density_ave_kmol/(5.d-1*(well_up%area(iup)+ &
              well_dn%area(idn)))
    ! Store flux calculation for consistency with transport
    if (save_flux) then
      well_up%ql(iup) = v_darcy
      well_up%ql_kmol(iup) = tot_mole_flux
    endif

    Res(1) = Res(1) + tot_mole_flux

    ! Gas flux
    density_kg_ave = 0.5d0*(well_up%gas%den(iup)+well_dn%gas%den(idn))
    ! Assuming the well is always vertical and gravity is in the
    ! (-) direction
    gravity_term = density_kg_ave * gravity * &
                    0.5d0*(well_grid%dh(iup)+well_grid%dh(idn))
    delta_pressure = well_up%pg(iup) - well_dn%pg(idn) + &
                      gravity_term

    up_scale = 0.d0
    dn_scale = 0.d0

    upwind = delta_pressure > 0.d0

    ! Only upwinding the perm here, not mobility ratio
    ! following the LumpedHarmonic BRAGFLO formulation
    if (upwind) then
      up_scale = 1.d0
      rel_perm = well_up%gas%kr(iup)
    else
      dn_scale = 1.d0
      rel_perm = well_dn%gas%kr(idn)
    endif

    !kmol/sec
    tot_mole_flux = perm_den_mu_area_ave_over_dist(2) * rel_perm * &
                    delta_pressure
    density_ave_kmol = density_kg_ave / fmw_comp(TWO_INTEGER)
    ! v_darcy [m/sec] = mole flux [kmol/sec] / den [kmol/m^3] / area[m^2]
    v_darcy = tot_mole_flux/density_ave_kmol/(5.d-1*(well_up%area(iup) + &
                        well_dn%area(idn)))
    ! Store flux calculation for consistency with transport
    if (save_flux) then
      well_up%qg(iup) = v_darcy
      well_up%qg_kmol(iup) = tot_mole_flux
    endif

    Res(2) = Res(2) + tot_mole_flux

end subroutine PMWellFlux

! ************************************************************************** !

subroutine PMWellFluxDerivative(pm_well,iup,idn,Jup,Jdn)
  !
  ! Computes the derivative of the internal flux terms for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 01/05/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: iup, idn, iphase, irow
  PetscReal :: Jup(pm_well%nphase,pm_well%nphase), &
               Jdn(pm_well%nphase,pm_well%nphase)

  PetscReal :: res_up(pm_well%nphase),res_dn(pm_well%nphase), &
               res_pert(pm_well%nphase)

  call PMWellFlux(pm_well,pm_well%well,pm_well%well,iup,idn,res_up,PETSC_FALSE)

  res_dn = res_up

  ! upgradient derivatives
  do iphase = 1,pm_well%nphase
    call PMWellFlux(pm_well,pm_well%well_pert(iphase),pm_well%well,iup,idn, &
                    res_pert,PETSC_FALSE)
    do irow = 1, pm_well%nphase
      Jup(irow,iphase) = (res_pert(irow)-res_up(irow)) / &
                         pm_well%pert(iup,iphase)
    enddo !irow
  enddo

  ! downgradient derivatives
  do iphase = 1,pm_well%nphase
    call PMWellFlux(pm_well,pm_well%well,pm_well%well_pert(iphase),iup,idn, &
                    res_pert,PETSC_FALSE)
    do irow = 1, pm_well%nphase
      Jdn(irow,iphase) = (res_pert(irow)-res_dn(irow)) / &
                         pm_well%pert(idn,iphase)
    enddo !irow
  enddo

end subroutine PMWellFluxDerivative

! ************************************************************************** !

subroutine PMWellBCFlux(pm_well,well,Res,save_flux)
  !
  ! Computes the boundary flux terms for the residual based on the
  ! chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 12/24/2021
  !

  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_WIPP_module

  implicit none

  class(pm_well_sequential_type) :: pm_well
  type(well_type) :: well
  PetscReal :: Res(2*pm_well%nphase)
  PetscBool :: save_flux

  type(option_type), pointer :: option
  type(well_grid_type), pointer :: well_grid
  type(well_reservoir_type), pointer :: reservoir
  class(characteristic_curves_type), pointer :: characteristic_curves
  class(sat_func_base_type), pointer :: saturation_function

  !MAN: clean these up
  PetscReal :: perm_ave_over_dist
  PetscReal :: gravity_term, delta_pressure
  PetscReal :: density_ave, tot_mole_flux
  PetscReal :: boundary_pressure, boundary_den
  PetscReal :: boundary_pg, boundary_krg, dn_scale
  PetscReal :: t,dwmol,dwp,dwt,Psat,visl,visg
  PetscReal :: Pc,dpc_dsatl,dkrl_dsatl,dkrg_dsatl
  PetscReal :: v_darcy,q,rel_perm
  PetscBool :: upwind
  PetscInt :: itop
  PetscErrorCode :: ierr
  PetscReal :: gravity

  option => pm_well%option

  well_grid => pm_well%well_grid
  reservoir => pm_well%well%reservoir

  gravity = option%gravity(Z_DIRECTION)
  t = 25.d0 !Constant temperature

  Res(:) = 0.d0
  ! Vertical well, no Klinkenberg.
  itop = pm_well%well_grid%nsegments
  if (pm_well%flow_soln%bh_p) then
    !Dirichlet pressure and saturation at the bottom

    characteristic_curves => pm_well%realization%patch% &
                            characteristic_curves_array(well%ccid(1))%ptr
    saturation_function => characteristic_curves%saturation_function

    ! Water Residual

    perm_ave_over_dist = well%permeability(1) / (well_grid%dh(1)/2.d0)
    boundary_pressure = well%bh_p
    gravity_term = well%liq%den(1) * gravity * &
                    well_grid%dh(1)/2.d0
    delta_pressure = boundary_pressure - well%pl(1) + gravity_term

    call EOSWaterSaturationPressure(t,Psat,ierr)
    call EOSWaterDensityBRAGFLO(t,boundary_pressure,PETSC_FALSE, &
                            boundary_den,dwmol,dwp,dwt,ierr)

    upwind = delta_pressure > 0.d0
    if (upwind) then
      call characteristic_curves%liq_rel_perm_function% &
            RelativePermeability(1.d0-well%bh_sg,rel_perm,dkrl_dsatl,option)
      call EOSWaterViscosity(t,boundary_pressure,Psat,visl,ierr)
    else
      dn_scale = 1.d0
      rel_perm = well%liq%kr(1)
      visl = well%liq%visc(1)
    endif

    ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
    !                    dP[Pa]]
    v_darcy = perm_ave_over_dist * rel_perm / visl * &
              delta_pressure
    if (upwind) then
      density_ave = boundary_den
    else
      density_ave = (well%liq%den(1)+boundary_den) / &
                    (2.d0 * fmw_comp(ONE_INTEGER))
    endif
    q = v_darcy * well%area(1)
    tot_mole_flux = q * density_ave
    ! Store boundary flux for consistency with transport
    if (save_flux) then
      well%ql_bc(1) = v_darcy
      well%ql_kmol_bc(1) = tot_mole_flux
    endif
    Res(1) = Res(1) + tot_mole_flux

    ! Gas Residual

    ! Capillary Pressure
    select type(sat_func => saturation_function)
      class is (sat_func_krp3_type)
        if (.not. option%flow%pct_updated) then
          sat_func%pct = sat_func%pct_a * well% permeability(1) ** &
                          sat_func%pct_exp
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
                CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
        else
          call sat_func% &
                CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
        endif
      class is (sat_func_krp4_type)
        if (.not. option%flow%pct_updated) then
          sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                          sat_func%pct_exp
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
                CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
        else
          call sat_func% &
                CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
        endif
      class is (sat_func_krp5_type)
        if (.not. option%flow%pct_updated) then
          sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                          sat_func%pct_exp
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
                CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
        else
          call sat_func% &
                CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
        endif
      class default
        call sat_func% &
            CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
    end select
    boundary_pg = boundary_pressure + Pc

    call characteristic_curves%gas_rel_perm_function% &
            RelativePermeability(1.d0-well%bh_sg,boundary_krg, &
            dkrg_dsatl,option)

    gravity_term = well%gas%den(1) * gravity * &
                    well_grid%dh(1)/2.d0
    delta_pressure = boundary_pg - well%pg(1) + gravity_term

    call EOSGasDensity(t,boundary_pg,boundary_den,ierr)

    upwind = delta_pressure > 0.d0
    if (upwind) then
      call characteristic_curves%gas_rel_perm_function% &
            RelativePermeability(1.d0-well%bh_sg,rel_perm,dkrl_dsatl,option)
      call EOSGasViscosity(t,boundary_pg,boundary_pg,boundary_den,visg,ierr)
    else
      dn_scale = 1.d0
      rel_perm = well%gas%kr(1)
      visg = well%gas%visc(1)
    endif

    v_darcy = perm_ave_over_dist * rel_perm/visg * &
              delta_pressure
    if (upwind) then
      density_ave = boundary_den
    else
      density_ave = (well%gas%den(1)+boundary_den) / &
                    (2.d0 *fmw_comp(TWO_INTEGER))
    endif
    q = v_darcy * well%area(1)
    tot_mole_flux = q * density_ave
    ! Store boundary flux for consistency with transport
    if (save_flux) then
      well%qg_bc(1) = v_darcy
      well%qg_kmol_bc(1) = tot_mole_flux
    endif
    Res(2) = Res(2) + tot_mole_flux

  else if (pm_well%flow_soln%bh_q) then
    !Neumann flux at the bottom
    v_darcy = well%bh_ql

    if (v_darcy > 0.d0) then
      density_ave = reservoir%den_l(1) / fmw_comp(ONE_INTEGER)
    else
      density_ave = well%liq%den(1) / fmw_comp(ONE_INTEGER)
    endif
    q = v_darcy * well%area(1)
    tot_mole_flux = q * density_ave
    ! Store boundary flux for consistency with transport
    if (save_flux) then
      well%ql_bc(1) = v_darcy
      well%ql_kmol_bc(1) = tot_mole_flux
    endif
    Res(1) = Res(1) + tot_mole_flux

    v_darcy = well%bh_qg

    if (v_darcy > 0.d0) then
      density_ave = reservoir%den_g(1) / fmw_comp(TWO_INTEGER)
    else
      density_ave = well%gas%den(1) / fmw_comp(TWO_INTEGER)
    endif
    q = v_darcy * well%area(1)
    tot_mole_flux = q * density_ave
    ! Store boundary flux for consistency with transport
    if (save_flux) then
      well%qg_bc(1) = v_darcy
      well%qg_kmol_bc(1) = tot_mole_flux
    endif
    Res(2) = Res(2) + tot_mole_flux
  else
    ! this should not happen once error messaging is updated
  endif

  if (pm_well%flow_soln%th_p) then
    !Dirichlet pressure and saturation at the top
    characteristic_curves => pm_well%realization%patch% &
                            characteristic_curves_array(well%ccid(itop))%ptr
    saturation_function => characteristic_curves%saturation_function

    ! Water Residual

    perm_ave_over_dist = well%permeability(itop) / (well_grid%dh(itop)/2.d0)
    boundary_pressure = well%th_p
    gravity_term = well%liq%den(itop) * gravity * &
                    well_grid%dh(itop)/2.d0
    delta_pressure = well%pl(itop) - boundary_pressure + gravity_term

    call EOSWaterSaturationPressure(t,Psat,ierr)
    call EOSWaterDensityBRAGFLO(t,boundary_pressure,PETSC_FALSE, &
                            boundary_den,dwmol,dwp,dwt,ierr)

    upwind = delta_pressure < 0.d0
    if (upwind) then
      call characteristic_curves%liq_rel_perm_function% &
            RelativePermeability(1.d0-well%th_sg,rel_perm,dkrl_dsatl,option)
      call EOSWaterViscosity(t,boundary_pressure,Psat,visl,ierr)
    else
      dn_scale = 1.d0
      rel_perm = well%liq%kr(itop)
      visl = well%liq%visc(itop)
    endif

    ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
    !                    dP[Pa]]
    v_darcy = perm_ave_over_dist * rel_perm / visl * &
              delta_pressure

    density_ave = (well%liq%den(itop)+boundary_den) / &
                  (2.d0 * fmw_comp(ONE_INTEGER))
    q = v_darcy * well%area(itop)
    tot_mole_flux = q * density_ave
    ! Store boundary flux for consistency with transport
    if (save_flux) then
      well%ql_bc(2) = v_darcy
      well%ql_kmol_bc(2) = tot_mole_flux
    endif
    Res(3) = Res(3) - tot_mole_flux

    ! Gas Residual

    ! Capillary Pressure
    select type(sat_func => saturation_function)
      class is (sat_func_krp3_type)
        if (.not. option%flow%pct_updated) then
          sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                          sat_func%pct_exp
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
                CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
        else
          call sat_func% &
                CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
        endif
      class is (sat_func_krp4_type)
        if (.not. option%flow%pct_updated) then
          sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                          sat_func%pct_exp
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
                CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
        else
          call sat_func% &
                CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
        endif
      class is (sat_func_krp5_type)
        if (.not. option%flow%pct_updated) then
          sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                          sat_func%pct_exp
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
                CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
        else
          call sat_func% &
                CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
        endif
      class default
        call sat_func% &
            CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
    end select
    boundary_pg = boundary_pressure + Pc

    call characteristic_curves%gas_rel_perm_function% &
            RelativePermeability(1.d0-well%th_sg,boundary_krg, &
            dkrg_dsatl,option)

    gravity_term = well%gas%den(itop) * gravity * &
                    well_grid%dh(itop)/2.d0
    delta_pressure = well%pg(itop) - boundary_pg + gravity_term

    call EOSGasDensity(t,boundary_pg,boundary_den,ierr)

    upwind = delta_pressure < 0.d0
    if (upwind) then
      call characteristic_curves%gas_rel_perm_function% &
            RelativePermeability(1.d0-well%th_sg,rel_perm,dkrl_dsatl,option)
      call EOSGasViscosity(t,boundary_pg,boundary_pg,boundary_den,visg,ierr)
    else
      dn_scale = 1.d0
      rel_perm = well%gas%kr(itop)
      visg = well%gas%visc(itop)
    endif

    v_darcy = perm_ave_over_dist * rel_perm/visg * &
              delta_pressure

    density_ave = (well%gas%den(itop)+boundary_den) / &
                  (2.d0 *fmw_comp(TWO_INTEGER))
    q = v_darcy * well%area(itop)
    tot_mole_flux = q * density_ave
    ! Store boundary flux for consistency with transport
    if (save_flux) then
      well%qg_bc(2) = v_darcy
      well%qg_kmol_bc(2) = tot_mole_flux
    endif
    Res(4) = Res(4) - tot_mole_flux
  else
    !Neumann flux at the top
    v_darcy = -well%th_ql

    ! Always take well density with tophole flux bc
    density_ave = well%liq%den(itop) / fmw_comp(ONE_INTEGER)
    q = v_darcy * well%area(itop)
    tot_mole_flux = q * density_ave
    ! Store boundary flux for consistency with transport
    if (save_flux) then
      well%ql_bc(2) = v_darcy
      well%ql_kmol_bc(2) = tot_mole_flux
    endif
    Res(3) = Res(3) - tot_mole_flux

    v_darcy = -well%th_qg

    ! Always take well density with tophole flux bc
    density_ave = well%gas%den(itop) / fmw_comp(TWO_INTEGER)
    q = v_darcy * well%area(itop)
    tot_mole_flux = q * density_ave
    ! Store boundary flux for consistency with transport
    if (save_flux) then
      well%qg_bc(2) = v_darcy
      well%qg_kmol_bc(2) = tot_mole_flux
    endif
    Res(4) = Res(4) - tot_mole_flux
  endif

end subroutine PMWellBCFlux

! ************************************************************************** !

subroutine PMWellBCFluxDerivative(pm_well,Jtop,Jbtm)
  !
  ! Computes the derivative of the boundary flux terms for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 01/05/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: Jtop(pm_well%flow_soln%ndof,pm_well%flow_soln%ndof), &
               Jbtm(pm_well%flow_soln%ndof,pm_well%flow_soln%ndof)

  PetscInt :: idof, irow
  PetscReal :: res(2*pm_well%flow_soln%ndof),res_pert(2*pm_well%flow_soln%ndof)

  call PMWellBCFlux(pm_well,pm_well%well,res,PETSC_FALSE)

  Jtop = 0.d0
  Jbtm = 0.d0

  ! downgradient derivatives
  do idof = 1,pm_well%nphase
    call PMWellBCFlux(pm_well,pm_well%well_pert(idof),res_pert,PETSC_FALSE)
    do irow = 1, pm_well%nphase
      Jbtm(irow,idof) = (res_pert(irow)-res(irow)) / &
                        pm_well%pert(1,idof)
    enddo
    do irow = 1, pm_well%nphase
      Jtop(irow,idof) = (res_pert(irow + pm_well%flow_soln%ndof)- &
                         res(irow + pm_well%flow_soln%ndof)) / &
                         pm_well%pert(pm_well%well_grid%nsegments,idof)
    enddo
  enddo


end subroutine PMWellBCFluxDerivative

! ************************************************************************** !

subroutine PMWellPerturb(pm_well)
  !
  ! Calculates the state variables for the perturbed well system.
  !
  ! Author: Michael Nole
  ! Date: 01/06/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscReal :: x(pm_well%well_grid%nsegments,pm_well%nphase), &
               pert(pm_well%well_grid%nsegments,pm_well%nphase)
  PetscInt :: i

  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-10

  x(:,ONE_INTEGER) = pm_well%well%pl
  x(:,TWO_INTEGER) = pm_well%well%gas%s

  ! Re-initialize perturbations
  pm_well%well_pert(ONE_INTEGER)%pl = pm_well%well%pl
  pm_well%well_pert(TWO_INTEGER)%pl = pm_well%well%pl
  pm_well%well_pert(ONE_INTEGER)%gas%s = pm_well%well%gas%s
  pm_well%well_pert(TWO_INTEGER)%gas%s = pm_well%well%gas%s

  pert(:,ONE_INTEGER) = perturbation_tolerance*x(:,ONE_INTEGER) + &
                        min_perturbation
  do i = 1,pm_well%well_grid%nsegments
    if (x(i,TWO_INTEGER) > 0.5d0) then
      pert(i,TWO_INTEGER) = -1.d0 * perturbation_tolerance
    else
      pert(i,TWO_INTEGER) = perturbation_tolerance
    endif
  enddo

  pm_well%well_pert(ONE_INTEGER)%pl = x(:,ONE_INTEGER) + pert(:,ONE_INTEGER)
  pm_well%well_pert(TWO_INTEGER)%gas%s = x(:,TWO_INTEGER) + pert(:,TWO_INTEGER)

  ! Update perturbed well properties
  call pm_well%UpdateFlowProperties(PETSC_TRUE, ONE_INTEGER)
  call pm_well%UpdateFlowProperties(PETSC_TRUE, TWO_INTEGER)
  ! Update perturbed source/sink term from the reservoir
  call PMWellUpdateWellQ(pm_well%well_pert(ONE_INTEGER), &
                         pm_well%well_pert(ONE_INTEGER)%reservoir)
  call PMWellUpdateWellQ(pm_well%well_pert(TWO_INTEGER), &
                         pm_well%well_pert(TWO_INTEGER)%reservoir)

  pm_well%pert = pert

end subroutine PMWellPerturb

! ************************************************************************** !


subroutine PMWellCheckConvergenceFlow(pm_well,n_iter,fixed_accum)
  !
  ! Checks flow solution convergence against prescribed tolerances.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 01/20/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: n_iter
  PetscReal :: fixed_accum(pm_well%flow_soln%ndof*pm_well%well_grid%nsegments)

  PetscReal, parameter :: zero_saturation = 1.d-15
  PetscReal, parameter :: zero_accumulation = 1.d-15

  type(well_soln_flow_type), pointer :: flow_soln
  character(len=MAXSTRINGLENGTH) :: out_string
  character(len=MAXSTRINGLENGTH) :: rsn_string
  PetscBool :: cnvgd_due_to_residual(pm_well%well_grid%nsegments* &
                                     pm_well%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_abs_res(pm_well%well_grid%nsegments* &
                                    pm_well%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_scaled_res(pm_well%well_grid%nsegments* &
                                       pm_well%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_update(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update_p(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update_s(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update_p(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update_s(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_on_pressure(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_on_saturation(pm_well%well_grid%nsegments)
  PetscReal :: update_p(pm_well%well_grid%nsegments) ! liquid pressure
  PetscReal :: update_s(pm_well%well_grid%nsegments) ! gas saturation
  PetscReal :: temp_real
  PetscReal :: max_scaled_residual,max_absolute_residual
  PetscReal :: max_relative_update_p,max_relative_update_s
  PetscReal :: max_absolute_update_p,max_absolute_update_s
  PetscInt :: loc_max_scaled_residual,loc_max_abs_residual
  PetscInt :: loc_max_rel_update_p,loc_max_rel_update_s
  PetscInt :: loc_max_abs_update_p,loc_max_abs_update_s
  PetscInt :: idof
  PetscInt :: k

  flow_soln => pm_well%flow_soln

  n_iter = n_iter + 1
  flow_soln%n_newton = flow_soln%n_newton + 1

  cnvgd_due_to_residual = PETSC_FALSE
  cnvgd_due_to_abs_res = PETSC_FALSE
  cnvgd_due_to_scaled_res = PETSC_FALSE
  cnvgd_due_to_update = PETSC_FALSE
  cnvgd_due_to_abs_update_p = PETSC_FALSE
  cnvgd_due_to_abs_update_s = PETSC_FALSE
  cnvgd_due_to_abs_update = PETSC_FALSE
  cnvgd_due_to_rel_update_p = PETSC_FALSE
  cnvgd_due_to_rel_update_s = PETSC_FALSE
  cnvgd_due_to_rel_update = PETSC_FALSE
  cnvgd_on_pressure = PETSC_FALSE
  cnvgd_on_saturation = PETSC_FALSE
  update_p = UNINITIALIZED_DOUBLE
  update_s = UNINITIALIZED_DOUBLE
  rsn_string = ''

  ! Update the residual
  flow_soln%residual = fixed_accum / pm_well%dt_flow
  call PMWellResidualFlow(pm_well)

  ! Update mass balance
  call PMWellMassBalance(pm_well)

  do k = 1,pm_well%well_grid%nsegments
    idof = flow_soln%ndof*(k-1)+1
    update_p(k) = flow_soln%update(idof)
    update_s(k) = flow_soln%update(idof+1)

    ! Absolute Solution Updates
    temp_real = dabs(update_p(k))
    if (temp_real < flow_soln%itol_abs_update_p) then
      cnvgd_due_to_abs_update_p(k) = PETSC_TRUE
    endif

    temp_real = dabs(update_s(k))
    if (temp_real > 0.d0) then
      if ((-1.d0*log10(temp_real)) >= &
           (-1.d0*log10(flow_soln%itol_abs_update_s))) then
        cnvgd_due_to_abs_update_s(k) = PETSC_TRUE
      endif
    else
      cnvgd_due_to_abs_update_s(k) = PETSC_TRUE
    endif

    ! Relative Solution Updates
    temp_real = dabs(update_p(k)/pm_well%well%pl(k))
    if (temp_real < flow_soln%itol_rel_update_p) then
      cnvgd_due_to_rel_update_p(k) = PETSC_TRUE
    endif
    temp_real = dabs(update_s(k)/pm_well%well%gas%s(k))
    if (temp_real < flow_soln%itol_rel_update_s) then
      cnvgd_due_to_rel_update_s(k) = PETSC_TRUE
    endif

    ! Liquid (water) Component
    if (dabs(fixed_accum(idof)) > zero_accumulation) then
      ! Absolute Residual
      temp_real = dabs(flow_soln%residual(idof))
      if (temp_real <= flow_soln%itol_abs_res) then
        cnvgd_due_to_abs_res(idof) = PETSC_TRUE
      endif

      ! Scaled Residual
      temp_real = dabs(flow_soln%residual(idof) / &
                       (fixed_accum(idof)/pm_well%dt_flow))
      if (temp_real <= flow_soln%itol_scaled_res) then
        cnvgd_due_to_scaled_res(idof) = PETSC_TRUE
      endif
    else
      cnvgd_due_to_abs_res(idof) = PETSC_TRUE
      cnvgd_due_to_scaled_res(idof) = PETSC_TRUE
    endif

    ! Gas (air) Component
    if (dabs(fixed_accum(idof+1)) > zero_accumulation) then
      ! Absolute Residual
      temp_real = dabs(flow_soln%residual(idof+1))
      if (temp_real <= flow_soln%itol_abs_res) then
        cnvgd_due_to_abs_res(idof+1) = PETSC_TRUE
      endif

      ! Scaled Residual
      temp_real = dabs(flow_soln%residual(idof+1) / &
                       (fixed_accum(idof+1)/pm_well%dt_flow))
      if (temp_real <= flow_soln%itol_scaled_res) then
        cnvgd_due_to_scaled_res(idof+1) = PETSC_TRUE
      endif
    else
      cnvgd_due_to_abs_res(idof+1) = PETSC_TRUE
      cnvgd_due_to_scaled_res(idof+1) = PETSC_TRUE
    endif
  enddo

  max_absolute_residual = maxval(dabs(flow_soln%residual))
  loc_max_abs_residual = maxloc(dabs(flow_soln%residual),1)

  max_scaled_residual = maxval(dabs(flow_soln%residual/ &
                                    (fixed_accum/pm_well%dt_flow)))
  loc_max_scaled_residual = maxloc(dabs(flow_soln%residual/ &
                                        (fixed_accum/pm_well%dt_flow)),1)

  max_absolute_update_p = maxval(dabs(update_p))
  loc_max_abs_update_p = maxloc(dabs(update_p),1)

  max_absolute_update_s = maxval(dabs(update_s))
  loc_max_abs_update_s = maxloc(dabs(update_s),1)

  max_relative_update_p = maxval(dabs(update_p/pm_well%well%pl))
  loc_max_rel_update_p = maxloc(dabs(update_p/pm_well%well%pl),1)

  max_relative_update_s = maxval(dabs(update_s/pm_well%well%gas%s))
  loc_max_rel_update_s = maxloc(dabs(update_s/pm_well%well%gas%s),1)

  do k = 1,pm_well%well_grid%nsegments*flow_soln%ndof
    if (cnvgd_due_to_scaled_res(k) .or. cnvgd_due_to_abs_res(k)) then
      cnvgd_due_to_residual(k) = PETSC_TRUE
    endif
  enddo
  do k = 1,pm_well%well_grid%nsegments
    if (cnvgd_due_to_abs_update_p(k) .or. &
        cnvgd_due_to_rel_update_p(k)) then
      cnvgd_on_pressure(k) = PETSC_TRUE
    endif
    if (cnvgd_due_to_abs_update_s(k) .or. &
        cnvgd_due_to_rel_update_s(k)) then
      cnvgd_on_saturation(k) = PETSC_TRUE
    endif
    if (cnvgd_on_pressure(k) .and. cnvgd_on_saturation(k)) then
      cnvgd_due_to_update(k) = PETSC_TRUE
    endif
  enddo

  if (all(cnvgd_due_to_abs_res)) then
    rsn_string = trim(rsn_string) // ' R '
  endif
  if (all(cnvgd_due_to_scaled_res)) then
    rsn_string = trim(rsn_string) // ' sR '
  endif
  if (all(cnvgd_due_to_abs_update)) then
    rsn_string = trim(rsn_string) // ' uP&uS '
  endif
  if (all(cnvgd_due_to_rel_update)) then
    rsn_string = trim(rsn_string) // ' ruP&ruS '
  endif

  if (all(cnvgd_due_to_residual) .and. all(cnvgd_due_to_update)) then
    flow_soln%converged = PETSC_TRUE
    flow_soln%not_converged = PETSC_FALSE
    pm_well%cumulative_dt_flow = pm_well%cumulative_dt_flow + pm_well%dt_flow
    pm_well%flow_soln%prev_soln%pl = pm_well%well%pl
    pm_well%flow_soln%prev_soln%sg = pm_well%well%gas%s
    !call PMWellUpdateSolutionFlow(pm_well)
  else
    flow_soln%converged = PETSC_FALSE
    flow_soln%not_converged = PETSC_TRUE
  endif

  if (pm_well%print_output) then
    write(out_string,'(i2," aR:",es10.2,"  sR:",es10.2,"  uP:", es10.2," &
          &  uS:",es10.2,"  ruP:",es10.2,"  ruS:",es10.2)') &
          n_iter,max_absolute_residual,max_scaled_residual, &
          max_absolute_update_p,max_absolute_update_s, &
          max_relative_update_p,max_relative_update_s
    call PrintMsg(pm_well%option,out_string)
    if (flow_soln%converged) then
      out_string = ' WELL FLOW Solution converged!  ---> ' // trim(rsn_string)
      call PrintMsg(pm_well%option,out_string)
    endif
  endif

end subroutine PMWellCheckConvergenceFlow

! ************************************************************************** !

subroutine PMWellCheckConvergenceTran(pm_well,n_iter,fixed_accum)
  !
  ! Checks transport solution convergence against prescribed tolerances.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: n_iter
  PetscReal :: fixed_accum(pm_well%tran_soln%ndof*pm_well%well_grid%nsegments)

  type(well_soln_tran_type), pointer :: soln
  character(len=MAXSTRINGLENGTH) :: out_string
  character(len=MAXSTRINGLENGTH) :: rsn_string
  PetscReal :: temp_real
  PetscBool :: cnvgd_due_to_residual(pm_well%well_grid%nsegments* &
                                     pm_well%tran_soln%ndof)
  PetscBool :: cnvgd_due_to_abs_res(pm_well%well_grid%nsegments* &
                                    pm_well%tran_soln%ndof)
  PetscBool :: cnvgd_due_to_scaled_res(pm_well%well_grid%nsegments* &
                                       pm_well%tran_soln%ndof)
  PetscBool :: cnvgd_due_to_update(pm_well%well_grid%nsegments*pm_well% &
                                   tran_soln%ndof)
  PetscReal :: vol_vec(pm_well%well_grid%nsegments*pm_well%tran_soln%ndof)
  PetscReal :: aq_mass_vec(pm_well%well_grid%nsegments*pm_well%tran_soln%ndof)
  PetscReal :: max_scaled_residual,max_absolute_residual
  PetscReal :: max_update
  PetscInt :: loc_max_scaled_residual,loc_max_abs_residual
  PetscInt :: loc_max_update
  PetscInt :: k,n,j,S,tag,last_rank,first_rank
  PetscInt :: isegment, ispecies
  PetscErrorCode :: ierr

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  soln => pm_well%tran_soln

  n_iter = n_iter + 1
  soln%n_newton = soln%n_newton + 1
  S = pm_well%well_grid%nsegments*pm_well%tran_soln%ndof

  cnvgd_due_to_residual = PETSC_FALSE
  cnvgd_due_to_abs_res = PETSC_FALSE
  cnvgd_due_to_scaled_res = PETSC_FALSE
  cnvgd_due_to_update = PETSC_FALSE
  rsn_string = ''

  ! Update the residual
  soln%residual = 0.d0
  if (any(pm_well%option%myrank == pm_well%well_grid%h_rank_id)) then
    soln%residual = fixed_accum/pm_well%dt_tran
    call PMWellResidualTran(pm_well)

    do k = 1,(pm_well%well_grid%nsegments*soln%ndof)
      ! Absolute Residual
      temp_real = dabs(soln%residual(k))
      if (temp_real < soln%itol_abs_res) then
        cnvgd_due_to_abs_res(k) = PETSC_TRUE
      endif
      ! Scaled Residual
      temp_real = dabs(soln%residual(k)/(fixed_accum(k)/pm_well%dt_tran))
      if (temp_real < soln%itol_scaled_res) then
        cnvgd_due_to_scaled_res(k) = PETSC_TRUE
      endif
    enddo

    ! Relative Update
    do n = 1,pm_well%well_grid%nsegments
      isegment = n
      do k = 1, soln%ndof
        ispecies = k
        j = ((isegment-1)*soln%ndof) + ispecies
        vol_vec(j) = pm_well%well%volume(isegment)
        aq_mass_vec(j) = pm_well%well%aqueous_mass(ispecies,isegment)
        temp_real = dabs(soln%update(j)*vol_vec(j)/aq_mass_vec(j))
        if (temp_real < soln%itol_rel_update) then
          cnvgd_due_to_update(j) = PETSC_TRUE
        endif
      enddo
    enddo

    max_absolute_residual = maxval(dabs(soln%residual))
    loc_max_abs_residual = maxloc(dabs(soln%residual),1)

    max_scaled_residual = maxval(dabs(soln%residual/ &
                                      (fixed_accum/pm_well%dt_tran)))
    loc_max_scaled_residual = maxloc(dabs(soln%residual/ &
                                          (fixed_accum/pm_well%dt_tran)),1)

    max_update = maxval(dabs(soln%update*vol_vec/aq_mass_vec))
    loc_max_update = maxloc(dabs(soln%update*vol_vec/aq_mass_vec),1)

    do k = 1,(pm_well%well_grid%nsegments*soln%ndof)
      if (cnvgd_due_to_scaled_res(k) .or. cnvgd_due_to_abs_res(k)) then
        cnvgd_due_to_residual(k) = PETSC_TRUE
      endif
    enddo
  endif

  call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
  first_rank = pm_well%well_comm%well_rank_list(1)
  last_rank = pm_well%well_comm%well_rank_list(pm_well%well_comm%commsize)
  if (pm_well%well_comm%commsize > 1) then
    tag = 0
    if (pm_well%well_comm%rank == last_rank) then
      call MPI_Send(cnvgd_due_to_abs_res,S,MPI_LOGICAL,0,tag, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call MPI_Send(cnvgd_due_to_scaled_res,S,MPI_LOGICAL,0,tag+1, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call MPI_Send(cnvgd_due_to_update,S,MPI_LOGICAL,0,tag+2, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
    endif
    if (pm_well%well_comm%rank == first_rank) then
      call MPI_Recv(cnvgd_due_to_abs_res,S,MPI_LOGICAL, &
                    last_rank,tag,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
      call MPI_Recv(cnvgd_due_to_scaled_res,S,MPI_LOGICAL, &
                    last_rank,tag+1,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
      call MPI_Recv(cnvgd_due_to_update,S,MPI_LOGICAL, &
                    last_rank,tag+2,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
    endif
  endif

  call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
  if (all(cnvgd_due_to_abs_res)) then
    rsn_string = trim(rsn_string) // ' aR '
  endif
  if (all(cnvgd_due_to_scaled_res)) then
    rsn_string = trim(rsn_string) // ' sR '
  endif
  if (all(cnvgd_due_to_update)) then
    rsn_string = trim(rsn_string) // ' rU '
  endif

  if (pm_well%well_comm%commsize > 1) then
    tag = 0
    if (pm_well%well_comm%rank == last_rank) then
      call MPI_Send(max_update,1,MPI_DOUBLE_PRECISION,0,tag, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call MPI_Send(max_absolute_residual,1,MPI_DOUBLE_PRECISION,0,tag+1, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call MPI_Send(max_scaled_residual,1,MPI_DOUBLE_PRECISION,0,tag+2, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
    endif
    if (pm_well%well_comm%rank == first_rank) then
      call MPI_Recv(max_update,1,MPI_DOUBLE_PRECISION, &
                    last_rank,tag,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
      call MPI_Recv(max_absolute_residual,1,MPI_DOUBLE_PRECISION, &
                    last_rank,tag+1,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
      call MPI_Recv(max_scaled_residual,1,MPI_DOUBLE_PRECISION, &
                    last_rank,tag+2,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
    endif
  endif

  call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
  write(out_string,'(i4,"    aR:",es10.3,"    sR:",es10.3,"    rU:", es10.3)')&
        n_iter,max_absolute_residual,max_scaled_residual, &
        max_update
  call PrintMsg(pm_well%option,out_string)

  if (pm_well%well_comm%commsize > 1) then
    tag = 0
    if (pm_well%well_comm%rank == last_rank) then
      call MPI_Send(cnvgd_due_to_residual,S,MPI_LOGICAL,0,tag, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call MPI_Send(cnvgd_due_to_update,S,MPI_LOGICAL,0,tag+1, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
    endif
    if (pm_well%well_comm%rank == first_rank) then
      call MPI_Recv(cnvgd_due_to_residual,S,MPI_LOGICAL, &
                    last_rank,tag,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
      call MPI_Recv(cnvgd_due_to_update,S,MPI_LOGICAL, &
                    last_rank,tag+1,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
    endif
  endif

  call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
  if (all(cnvgd_due_to_residual) .and. all(cnvgd_due_to_update)) then
    soln%converged = PETSC_TRUE
    soln%not_converged = PETSC_FALSE
    out_string = ' WELL TRAN Solution converged!  ---> ' // trim(rsn_string)
    call PrintMsg(pm_well%option,out_string)
    pm_well%cumulative_dt_tran = pm_well%cumulative_dt_tran + pm_well%dt_tran
  else
    soln%converged = PETSC_FALSE
    soln%not_converged = PETSC_TRUE
  endif

  pm_well%tran_soln%cut_ts_flag = PETSC_FALSE

end subroutine PMWellCheckConvergenceTran

! ************************************************************************** !

subroutine PMWellResidualFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscInt :: i, iup, idn
  PetscReal :: res_accum(pm_well%nphase)
  PetscReal :: res_src_sink(pm_well%nphase)
  PetscReal :: res_flux(pm_well%nphase)
  PetscReal :: res_flux_bc(2*pm_well%nphase)
  PetscReal :: res_temp(pm_well%flow_soln%ndof*pm_well%well_grid%nsegments)

  res_accum = 0.d0
  res_src_sink = 0.d0
  res_flux = 0.d0
  res_flux_bc = 0.d0

  res_temp(:) = pm_well%flow_soln%residual(:)

  call PMWellBCFlux(pm_well,pm_well%well,res_flux_bc,PETSC_TRUE)

  do i = 1,pm_well%well_grid%nsegments
    iup = i
    idn = i + 1

    ! Accumulation Term
    call PMWellAccumulationFlow(pm_well,pm_well%well,i,res_accum)

    res_temp(pm_well%flow_soln%ndof*(i-1)+1) = &
            res_temp(pm_well%flow_soln%ndof*(i-1)+1) + &
            res_accum(ONE_INTEGER)

    ! Source/Sink Term
    call PMWellSrcSink(pm_well,pm_well%well,i,res_src_sink)

    res_temp(pm_well%flow_soln%ndof*(i-1)+1) = &
          res_temp(pm_well%flow_soln%ndof*(i-1)+1) + &
          res_src_sink(ONE_INTEGER)

    ! Flux Term
    if (i < pm_well%well_grid%nsegments) then
      call PMWellFlux(pm_well,pm_well%well,pm_well%well,iup,idn,res_flux, &
                      PETSC_TRUE)
    endif

    if (i == 1) then
      ! Water mass residual in cell i+1: Subtract flux to i+1 cell
      res_temp(pm_well%flow_soln%ndof*i+1) = &
            res_temp(pm_well%flow_soln%ndof*i+1) &
            - res_flux(1)

      ! Water mass residual in cell i: Add flux in from BC,
      ! add flux to i+1 cell
      res_temp(pm_well%flow_soln%ndof*(i-1)+1) = &
            res_temp(pm_well%flow_soln%ndof*(i-1)+1) &
            - res_flux_bc(1) + res_flux(1)

    elseif (i < pm_well%well_grid%nsegments) then
      ! Water mass residual in cell i: Subtract flux to i+1 cell
      res_temp(pm_well%flow_soln%ndof*(i-1)+1) = &
            res_temp(pm_well%flow_soln%ndof*(i-1)+1) &
            + res_flux(1)
      ! Water mass residual in cell i+1: Add flux to i+1 cell
      res_temp(pm_well%flow_soln%ndof*i+1) = &
            res_temp(pm_well%flow_soln%ndof*i+1) &
            - res_flux(1)
    else
      ! Water mass residual in cell i: Subtract flux to BC
      res_temp(pm_well%flow_soln%ndof*(i-1)+1) = &
            res_temp(pm_well%flow_soln%ndof*(i-1)+1) &
            - res_flux_bc(3)
    endif

    if (pm_well%nphase == 2) then

      res_temp(pm_well%flow_soln%ndof*(i-1)+2) = &
            res_temp(pm_well%flow_soln%ndof*(i-1)+2) + &
            res_accum(TWO_INTEGER)

      res_temp(pm_well%flow_soln%ndof*(i-1)+2) = &
            res_temp(pm_well%flow_soln%ndof*(i-1)+2) + &
            res_src_sink(TWO_INTEGER)

      if (i == 1) then
        ! Air mass residual in cell i+1: Subtract flux to i+1 cell
        res_temp(pm_well%flow_soln%ndof*i+2) = &
              res_temp(pm_well%flow_soln%ndof*i+2) &
              - res_flux(2)
        ! Air mass residual in cell i: Subtract flux in from BC,
        ! add flux to i+1 cell
        res_temp(pm_well%flow_soln%ndof*(i-1)+2) = &
              res_temp(pm_well%flow_soln%ndof*(i-1)+2) &
              - res_flux_bc(2) + res_flux(2)
      elseif (i < pm_well%well_grid%nsegments) then
        ! Air mass residual in cell i: Subtract flux to i+1 cell
        res_temp(pm_well%flow_soln%ndof*(i-1)+2) = &
              res_temp(pm_well%flow_soln%ndof*(i-1)+2) &
              + res_flux(2)
        ! Air mass residual in cell i+1: Add flux to i+1 cell
        res_temp(pm_well%flow_soln%ndof*i+2) = &
              res_temp(pm_well%flow_soln%ndof*i+2) &
              - res_flux(2)
      else
        ! Air mass residual in cell i: Subtract flux to BC
        res_temp(pm_well%flow_soln%ndof*(i-1)+2) = &
              res_temp(pm_well%flow_soln%ndof*(i-1)+2) &
              - res_flux_bc(4)
      endif
    endif
  enddo
  pm_well%flow_soln%residual(:) = res_temp(:)

end subroutine PMWellResidualFlow

! ************************************************************************** !

subroutine PMWellQISolveTran(pm_well)
  !
  ! This routine is called by NWTResidual() in nw_transport.F90
  !
  ! Author: Jennifer M. Frederick
  ! Date: 03/13/2023
  !

  implicit none

  class(pm_well_qi_type) :: pm_well

  PetscReal :: curr_time
  PetscErrorCode :: ierr

  ierr = 0
  pm_well%tran_soln%cut_ts_flag = PETSC_FALSE
  curr_time = pm_well%option%time

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return


  call PMWellUpdateReservoirConcTran(pm_well)
  if (initialize_well_tran) then
    call PMWellInitializeWellTran(pm_well)
  endif

  call PMWellUpdatePropertiesTran(pm_well)
  pm_well%dt_tran = pm_well%option%tran_dt

  call PMWellSolveTran(pm_well,ierr)
  if (pm_well%tran_soln%cut_ts_flag) return

  call PMWellUpdateReservoirSrcSinkTran(pm_well)
  call PMWellUpdateReservoirConcTran(pm_well)

end subroutine PMWellQISolveTran

! ************************************************************************** !

subroutine PMWellSolveTran(pm_well,ierr)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/22/2022
  !

  implicit none

  class(pm_well_qi_type) :: pm_well
  PetscErrorCode :: ierr

  type(well_soln_tran_type), pointer :: soln
  character(len=MAXSTRINGLENGTH) :: out_string
  PetscLogDouble :: log_start_time, log_end_time
  PetscReal :: res_fixed(pm_well%tran_soln%ndof*pm_well%well_grid%nsegments)
  PetscReal :: master_dt
  PetscInt :: n_iter, ts_cut
  PetscInt :: istart, iend
  PetscInt :: k

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  soln => pm_well%tran_soln

  ierr = 0
  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

  ts_cut = 0

  pm_well%cumulative_dt_tran = 0.d0
  soln%converged = PETSC_FALSE
  soln%not_converged = PETSC_TRUE

  master_dt = pm_well%option%tran_dt

  do while (pm_well%cumulative_dt_tran < master_dt)

    call PMWellPreSolveTran(pm_well,master_dt)

    n_iter = 0

    do while (soln%not_converged)
      if (n_iter > (soln%max_iter-1)) then
        soln%cut_timestep = PETSC_TRUE
        soln%cut_ts_flag = PETSC_TRUE
        out_string = ' Maximum number of TRAN Newton iterations reached. &
                      &Cutting timestep!'
        call PrintMsg(pm_well%option,out_string)
        call PMWellCutTimestepTran(pm_well)
        ! make sure well-flow doesn't get re-solved:
        pm_well%update_for_flow_qi_coupling = PETSC_TRUE
        return
      endif
      if (ts_cut > soln%max_ts_cut) then
        pm_well%realization%option%io_buffer = &
          ' Maximum timestep cuts reached in PM Well TRAN. Solution has not &
           &converged. Exiting.'
        if (pm_well%print_well) then
          call PMWellOutputSequential(pm_well)
        endif
        call PrintErrMsg(pm_well%realization%option)
      endif

      soln%residual = 0.d0
      if (any(pm_well%option%myrank == pm_well%well_grid%h_rank_id)) then
        ! Get fixed accumulation term (not yet divided by dt)
        do k = 1,pm_well%well_grid%nsegments
          istart = soln%ndof*(k-1)+1
          iend = soln%ndof*k
          call PMWellAccumulationTran(pm_well,k,res_fixed(istart:iend))
        enddo
        soln%residual = res_fixed / pm_well%dt_tran
      endif
      call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call PMWellNewtonTran(pm_well,n_iter)
      call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call PMWellCheckConvergenceTran(pm_well,n_iter,res_fixed)

    enddo

    ! try to increase the time step, if possible
    if (soln%converged) then
      pm_well%dt_tran = soln%ts_ramp_factor * pm_well%dt_tran
    endif
    if (pm_well%dt_tran > master_dt) then
      pm_well%dt_tran = master_dt
    endif

    ! if pm_well next time step will overstep master_dt, then correct it
    if (pm_well%cumulative_dt_tran + pm_well%dt_tran > master_dt) then
      pm_well%dt_tran = master_dt - pm_well%cumulative_dt_tran
    endif

    soln%n_steps = soln%n_steps + 1

  enddo

  call PMWellPostSolveTran(pm_well)

  call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

end subroutine PMWellSolveTran

! ************************************************************************** !


subroutine PMWellResidualTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  ! at this time, the tran_soln%residual vector has been zero'd out and has
  ! been loaded with the fixed accumulation divided by current dt

  ! update the auxiliary variables (runs through the equilibrium
  ! diss/precip/sorb routine) - we only need to update the aqueous_conc from
  ! the aqueous_mass value, or vice versa?

  call PMWellResidualTranAccum(pm_well)

  ! calculate the source/sink terms (in/out of well segments)
  call PMWellResidualTranSrcSink(pm_well)

  ! calculate the rxn terms (decay/ingrowth)
  call PMWellResidualTranRxn(pm_well)

  ! calculate the flux terms
  call PMWellResidualTranFlux(pm_well)

end subroutine PMWellResidualTran

! ************************************************************************** !

subroutine PMWellResidualTranAccum(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend
  PetscReal :: Res(pm_well%nspecies)

  ! porosity in [m^3-void/m^3-bulk]
  ! saturation in [m^3-liq/m^3-void]
  ! volume in [m^3-bulk]
  ! aqueous conc in [mol-species/m^3-liq]
  ! residual in [mol-species/sec]

  ! calculate the accumulation term as:
  ! residual = (Res_accum/tran_dt)
  !            - residual(which is already = fixed_accum/dt)

  do isegment = 1,pm_well%well_grid%nsegments

    offset = (isegment-1)*pm_well%nspecies
    istart = offset + 1
    iend = offset + pm_well%nspecies

    do ispecies = 1,pm_well%nspecies
      k = ispecies
      Res(k) = pm_well%well%volume(isegment) * pm_well%well%phi(isegment) * &
               pm_well%well%liq%s(isegment) * &
               pm_well%well%aqueous_conc(ispecies,isegment)
      Res(k) = Res(k) / pm_well%dt_tran
    enddo

    pm_well%tran_soln%residual(istart:iend) = Res(:) - &
                                        pm_well%tran_soln%residual(istart:iend)
  enddo

end subroutine PMWellResidualTranAccum

! ************************************************************************** !

subroutine PMWellResidualTranSrcSink(pm_well)
  !
  ! Calculates the source sink terms (Q in/out of well) for the transport
  ! residual equation.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/06/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend
  PetscReal :: Res(pm_well%nspecies)

  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: resr
  PetscReal :: coef_Qin, coef_Qout ! into well, out of well
  PetscReal :: Qin, Qout
  PetscReal :: den_avg

  ! Q src/sink is in [kmol-liq/sec]
  ! FMWH2O is in [kg-liq/kmol-liq] where liq = water
  ! density is in [kg-liq/m^3-liq] where liq = water
  ! aqueous conc in [mol-species/m^3-liq]
  ! residual in [mol-species/sec]

  ! From the flow solution:
  ! + Q goes into well from reservoir
  ! - Q goes out of well into reservoir

  well => pm_well%well
  resr => pm_well%well%reservoir

  do isegment = 1,pm_well%well_grid%nsegments

    den_avg = 0.5d0*(well%liq%den(isegment)+resr%den_l(isegment))
    ! units of coef = [m^3-liq/sec]
    if (well%liq%Q(isegment) < 0.d0) then ! Q out of well
      coef_Qin = 0.d0
      coef_Qout = well%liq%Q(isegment)*FMWH2O/den_avg
    else ! Q into well
    !            [kmol-liq/sec]*[kg-liq/kmol-liq]/[kg-liq/m^3-liq]
      coef_Qin = well%liq%Q(isegment)*FMWH2O/den_avg
      coef_Qout = 0.d0
    endif

    offset = (isegment-1)*pm_well%nspecies
    istart = offset + 1
    iend = offset + pm_well%nspecies

    do ispecies = 1,pm_well%nspecies
      k = ispecies
      Qin = coef_Qin*resr%aqueous_conc(ispecies,isegment)
      Qout = coef_Qout*well%aqueous_conc(ispecies,isegment)
      Res(k) = Qin + Qout
    enddo

    pm_well%tran_soln%residual(istart:iend) = &
                          pm_well%tran_soln%residual(istart:iend) + Res(:)
  enddo

end subroutine PMWellResidualTranSrcSink

! ************************************************************************** !

subroutine PMWellResidualTranRxn(pm_well)
  !
  ! Calculates the decay/ingrowth of radioactive species for the transport
  ! residual equation.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/06/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend, parent_id
  PetscReal :: Res(pm_well%nspecies)

  ! decay_rate in [1/sec]
  ! aqueous mass in [mol-species]
  ! residual in [mol-species/sec]

  do isegment = 1,pm_well%well_grid%nsegments

    offset = (isegment-1)*pm_well%nspecies
    istart = offset + 1
    iend = offset + pm_well%nspecies

    do ispecies = 1,pm_well%nspecies
      k = ispecies
      ! Add in species decay
      Res(k) = -(pm_well%well%species_decay_rate(k)* &
                          pm_well%well%aqueous_mass(ispecies,isegment))
      ! Add in contribution from parent (if exists)
      parent_id = pm_well%well%species_parent_id(ispecies)
      if (parent_id > 0) then
        Res(k) = Res(k) + (pm_well%well%species_parent_decay_rate(k)* &
                 pm_well%well%aqueous_mass(parent_id,isegment))
      endif
    enddo

    pm_well%tran_soln%residual(istart:iend) = &
                          pm_well%tran_soln%residual(istart:iend) + Res(:)
  enddo

end subroutine PMWellResidualTranRxn

! ************************************************************************** !

subroutine PMWellResidualTranFlux(pm_well)
  !
  ! Calculates the interior and BC flux terms for the transport residual
  ! equation.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/24/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend
  PetscInt :: n_up, n_dn
  PetscReal :: area_up, area_dn
  PetscReal :: q_up, q_dn
  PetscReal :: conc
  PetscReal :: diffusion
  PetscReal :: Res(pm_well%nspecies)
  PetscReal :: Res_up(pm_well%nspecies), Res_dn(pm_well%nspecies)

  ! residual in [mol-species/sec]
  ! area in [m2-bulk]
  ! q_up, d_dn in [m3-liq/m2-bulk-sec]
  ! conc in [mol-species/m3-liq]

  ! NOTE: The up direction is towards well top, and the dn direction is
  !       towards the well bottom.
  !       +q flows up the well
  !       -q flows down the well

  n_dn = +1
  n_up = -1

  diffusion = 0.d0 ! for now, since WIPP has no diffusion

  ! ----------------------------------------INTERIOR-FLUXES------------------

  do isegment = 2,(pm_well%well_grid%nsegments-1)

    Res(:) = 0.d0; Res_up(:) = 0.d0; Res_dn(:) = 0.d0

    area_up = 0.5d0 * (pm_well%well%area(isegment) + &
              pm_well%well%area(isegment+1))
    area_dn = 0.5d0 * (pm_well%well%area(isegment) + &
              pm_well%well%area(isegment-1))

    q_up = pm_well%well%ql(isegment)
    q_dn = pm_well%well%ql(isegment-1)

    offset = (isegment-1)*pm_well%nspecies
    istart = offset + 1
    iend = offset + pm_well%nspecies

    do ispecies = 1,pm_well%nspecies
      k = ispecies

      ! north surface:
      if (q_up < 0.d0) then ! flow is down well
        conc = pm_well%well%aqueous_conc(k,isegment+1)
        Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
      elseif (q_up > 0.d0) then ! flow is up well
        conc = pm_well%well%aqueous_conc(k,isegment)
        Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
      else ! q_up = 0
        Res_up(k) = (n_up*area_up)*(0.d0 - diffusion)
      endif

      ! south surface:
      if (q_dn < 0.d0) then ! flow is down well
        conc = pm_well%well%aqueous_conc(k,isegment)
        Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
      elseif (q_dn > 0.d0) then ! flow up well
        conc = pm_well%well%aqueous_conc(k,isegment-1)
        Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
      else ! q_dn = 0
        Res_dn(k) = (n_dn*area_dn)*(0.d0 - diffusion)
      endif

      Res(k) = Res_up(k) + Res_dn(k)
    enddo

    pm_well%tran_soln%residual(istart:iend) = &
                      pm_well%tran_soln%residual(istart:iend) + Res(:)
  enddo

  ! ----------------------------------------BOUNDARY-FLUXES------------------

  ! ----- bottom of well -----
  isegment = 1
  Res(:) = 0.d0; Res_up(:) = 0.d0; Res_dn(:) = 0.d0

  area_up = 0.5d0 * (pm_well%well%area(isegment) + &
            pm_well%well%area(isegment+1))
  area_dn = pm_well%well%area(isegment)

  q_up = pm_well%well%ql(isegment)
  q_dn = pm_well%well%ql_bc(1) ! bottom of hole ql

  offset = (isegment-1)*pm_well%nspecies ! = 0
  istart = offset + 1
  iend = offset + pm_well%nspecies

  do ispecies = 1,pm_well%nspecies
    k = ispecies

    ! north surface:
    if (q_up < 0.d0) then ! flow is down the well
      conc = pm_well%well%aqueous_conc(k,isegment+1)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    elseif (q_up > 0.d0) then ! flow is up the well
      conc = pm_well%well%aqueous_conc(k,isegment)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    else ! q_up = 0
      Res_up(k) = (n_up*area_up)*(0.d0 - diffusion)
    endif

    ! south surface:
    if (q_dn < 0.d0) then ! flow is down the well
      conc = pm_well%well%aqueous_conc(k,isegment)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    elseif (q_dn > 0.d0) then ! flow is up the well
      conc = pm_well%well%reservoir%aqueous_conc(k,isegment)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    else ! q_dn = 0
      Res_dn(k) = (n_dn*area_dn)*(0.d0 - diffusion)
    endif

    Res(k) = Res_up(k) + Res_dn(k)
  enddo

  pm_well%tran_soln%residual(istart:iend) = &
                      pm_well%tran_soln%residual(istart:iend) + Res(:)


  ! ----- top of well -----
  isegment = pm_well%well_grid%nsegments
  Res(:) = 0.d0; Res_up(:) = 0.d0; Res_dn(:) = 0.d0

  area_up = pm_well%well%area(isegment)
  area_dn = 0.5d0 * (pm_well%well%area(isegment) + &
            pm_well%well%area(isegment-1))

  q_up = pm_well%well%ql_bc(2) ! top of hole ql
  q_dn = pm_well%well%ql(isegment-1)

  offset = (isegment-1)*pm_well%nspecies
  istart = offset + 1
  iend = offset + pm_well%nspecies

  do ispecies = 1,pm_well%nspecies
    k = ispecies

    ! north surface:
    if (q_up < 0.d0) then ! flow is down the well
      conc = pm_well%well%aqueous_conc_th(ispecies)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    elseif (q_up > 0.d0) then ! flow is up the well
      conc = pm_well%well%aqueous_conc(k,isegment)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    else ! q_up = 0
      Res_up(k) = (n_up*area_up)*(0.d0 - diffusion)
    endif

    ! south surface:
    if (q_dn < 0.d0) then ! flow is down the well
      conc = pm_well%well%aqueous_conc(k,isegment)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    elseif (q_dn > 0.d0) then ! flow is up the well
      conc = pm_well%well%aqueous_conc(k,isegment-1)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    else ! q_dn = 0
      Res_dn(k) = (n_dn*area_dn)*(0.d0 - diffusion)
    endif

    Res(k) = Res_up(k) + Res_dn(k)
  enddo

  pm_well%tran_soln%residual(istart:iend) = &
                      pm_well%tran_soln%residual(istart:iend) + Res(:)

end subroutine PMWellResidualTranFlux

! ************************************************************************** !

subroutine PMWellJacobianFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscInt :: local_id
  PetscInt :: iconn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: i,k
  Vec, parameter :: null_vec = tVec(0)
  PetscReal :: Jup(pm_well%nphase,pm_well%nphase), &
               Jdn(pm_well%nphase,pm_well%nphase), &
               Jtop(pm_well%nphase,pm_well%nphase), &
               Jbtm(pm_well%nphase,pm_well%nphase), &
               Jtmp(pm_well%nphase,pm_well%nphase), &
               Jac(pm_well%nphase*pm_well%well_grid%nsegments, &
                   pm_well%nphase*pm_well%well_grid%nsegments)

  pm_well%flow_soln%Jacobian = 0.d0
  Jac = 0.d0
  Jup = 0.d0
  Jtop = 0.d0
  Jbtm = 0.d0
  Jtmp = 0.d0

  ! Accumulation Term ------------------------------------
  do local_id = 1,pm_well%well_grid%nsegments
    call PMWellAccumDerivative(pm_well,local_id,Jup)
    call PMWellFillJacFlow(pm_well,Jac,Jup,local_id,local_id)
  enddo

  ! Source/Sink Term
  do local_id = 1,pm_well%well_grid%nsegments
    call PMWellSrcSinkDerivative(pm_well,local_id,Jup)
    call PMWellFillJacFlow(pm_well,Jac,Jup,local_id,local_id)
  enddo

  ! Interior Flux Terms -----------------------------------
  do iconn = 1,pm_well%well_grid%nconnections

    local_id_up = iconn
    local_id_dn = iconn+1

    call PMWellFluxDerivative(pm_well,local_id_up,local_id_dn,Jup,Jdn)

    Jtmp = Jup
    call PMWellFillJacFlow(pm_well,Jac,Jtmp,local_id_up,local_id_up)

    Jtmp = Jdn
    call PMWellFillJacFlow(pm_well,Jac,Jtmp,local_id_up,local_id_dn)

    Jup = -Jup
    Jdn = -Jdn
    Jtmp = Jdn
    call PMWellFillJacFlow(pm_well,Jac,Jtmp,local_id_dn,local_id_dn)

    Jtmp = Jup
    call PMWellFillJacFlow(pm_well,Jac,Jtmp,local_id_dn,local_id_up)

  enddo

  ! Boundary Flux Terms -----------------------------------
  local_id = 1
  call PMWellBCFluxDerivative(pm_well,Jtop,Jbtm)
  Jbtm = -Jbtm
  call PMWellFillJacFlow(pm_well,Jac,Jbtm,local_id,local_id)

  local_id = pm_well%well_grid%nsegments
  Jtop = -Jtop
  call PMWellFillJacFlow(pm_well,Jac,Jtop,local_id,local_id)

  !pm_well_ni_count = pm_well_ni_count + 1

  do i = 1,pm_well%nphase*pm_well%well_grid%nsegments
    do k = 1,pm_well%nphase*pm_well%well_grid%nsegments
      pm_well%flow_soln%Jacobian(i,k) = Jac(i,k)
    enddo
  enddo

  !pm_well%flow_soln%Jacobian = Jac

end subroutine PMWellJacobianFlow

! ************************************************************************** !
subroutine PMWellFillJacFlow(pm_well,Jac,Jtmp,id1,id2)
  !
  ! Author: Michael Nole
  ! Date: 01/10/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: Jtmp(pm_well%nphase,pm_well%nphase), &
               Jac(pm_well%nphase*pm_well%well_grid%nsegments, &
                   pm_well%nphase*pm_well%well_grid%nsegments)
  PetscInt :: id1,id2

  PetscInt :: i,j

  do i = 1,pm_well%nphase
    do j = 1,pm_well%nphase
      Jac((id1-1)*pm_well%nphase+i,(id2-1)*pm_well%nphase+j) = &
      Jac((id1-1)*pm_well%nphase+i,(id2-1)*pm_well%nphase+j) + Jtmp(i,j)
    enddo
  enddo

end subroutine PMWellFillJacFlow

! ************************************************************************** !

subroutine PMWellJacobianTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/14/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscReal :: Jblock(pm_well%nspecies,pm_well%nspecies)
  PetscInt :: k, nspecies
  PetscInt :: jstart, jend

  nspecies = pm_well%nspecies
  pm_well%tran_soln%Jacobian(:,:) = 0.d0

  do k = 1,pm_well%well_grid%nsegments

    Jblock(:,:) = 0.d0

    call PMWellJacTranAccum(pm_well,Jblock,k)

    call PMWellJacTranSrcSink(pm_well,Jblock,k)

    call PMWellJacTranFlux(pm_well,Jblock,k)

    call PMWellJacTranRxn(pm_well,Jblock,k)

    ! place JBlock into full Jac based on isegment
    jstart = (k-1)*nspecies + 1
    jend = jstart + nspecies - 1
    pm_well%tran_soln%Jacobian(jstart:jend,jstart:jend) = Jblock

  enddo

end subroutine PMWellJacobianTran

! ************************************************************************** !

subroutine PMWellJacTranAccum(pm_well,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/14/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: Jblock(pm_well%nspecies,pm_well%nspecies)
  PetscInt :: isegment

  PetscReal :: vol_dt
  PetscInt :: istart, iend, ispecies

  ! units of Jac = [m^3-bulk/sec]
  ! units of volume = [m^3-bulk]
  ! units of tran_dt = [sec]

  vol_dt = pm_well%well%volume(isegment)/pm_well%dt_tran

  istart = 1
  iend = pm_well%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + vol_dt

  enddo

end subroutine PMWellJacTranAccum

! ************************************************************************** !

subroutine PMWellJacTranSrcSink(pm_well,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: Jblock(pm_well%nspecies,pm_well%nspecies)
  PetscInt :: isegment

  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: resr
  PetscInt :: istart, iend, ispecies
  PetscReal :: Qin, Qout
  PetscReal :: SSin, SSout, SS
  PetscReal :: vol, den_avg

  ! units of Jac = [m^3-bulk/sec]
  ! units of volume = [m^3-bulk]
  ! units of liq%Q = [kmol-liq/sec]
  ! units of FMWH2O = [kg-liq/kmol-liq]
  ! units of density = [kg-liq/m^3-liq]
  ! units of Qin = [m^3-liq/sec]
  ! units of SS = [m^3-bulk/sec]

  well => pm_well%well
  resr => pm_well%well%reservoir

  ! From the flow solution:
  ! + Q goes into well from reservoir
  ! - Q goes out of well into reservoir

  vol = pm_well%well%volume(isegment)
  den_avg = 0.5d0*(well%liq%den(isegment)+resr%den_l(isegment))

  ! units of Qin/out = [m^3-liq/sec]
  if (well%liq%Q(isegment) < 0.d0) then ! Q out of well
    Qin = 0.d0
    Qout = well%liq%Q(isegment)*FMWH2O/den_avg
    if (well%liq%s(isegment) < 1.d-40) then
      pm_well%option%io_buffer = 'HINT: The liquid saturation is zero. &
        &Division by zero will occur in PMWellJacTranSrcSink().'
      call PrintMsg(pm_well%option)
    endif
  else ! Q into well
    Qin = 0.d0  !well%liq%Q(isegment)*FMWH2O/den_avg
    Qout = 0.d0
    if (resr%s_l(isegment) < 1.d-40) then
      pm_well%option%io_buffer = 'HINT: The liquid saturation is zero. &
        &Division by zero will occur in PMWellJacTranSrcSink().'
      call PrintMsg(pm_well%option)
    endif
  endif

  SSin = Qin / (resr%e_por(isegment)*resr%s_l(isegment))    ! [m3-bulk/sec]
  SSout = Qout / (well%phi(isegment)*well%liq%s(isegment))  ! [m3-bulk/sec]
  SS = SSin + SSout

  istart = 1
  iend = pm_well%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + vol*(SS/vol)

  enddo

end subroutine PMWellJacTranSrcSink

! ************************************************************************** !

subroutine PMWellJacTranFlux(pm_well,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: Jblock(pm_well%nspecies,pm_well%nspecies)
  PetscInt :: isegment

  type(well_type), pointer :: well
  PetscInt :: istart, iend, ispecies
  PetscInt :: n_up, n_dn
  PetscReal :: d_diffusion_dM
  PetscReal :: J_up, J_dn
  PetscReal :: area_up, area_dn
  PetscReal :: sat_up, sat_dn, por_up, por_dn
  PetscReal :: u_up, u_dn

  ! units of Jac = [m^3-bulk/sec]
  ! area in [m2-bulk]
  ! q in [m3-liq/m2-bulk-sec]
  ! u in [m-liq/sec]
  ! sat in [m2-liq/m2-void]

  well => pm_well%well

  ! NOTE: The up direction is towards well top, and the dn direction is
  !       towards the well bottom.
  !       +q flows up the well
  !       -q flows down the well

  n_dn = +1
  n_up = -1

  d_diffusion_dM = 0.d0 ! for now, since WIPP has no diffusion

  if ((isegment > 1) .and. (isegment < pm_well%well_grid%nsegments)) then
  ! ----------------------------------------INTERIOR-FLUXES------------------

    ! define face values with arithmetic averages:
    area_up = 0.5d0 * (well%area(isegment) + well%area(isegment+1))
    area_dn = 0.5d0 * (well%area(isegment) + well%area(isegment-1))
    por_up = 0.5d0 * (well%phi(isegment) + well%phi(isegment+1))
    por_dn = 0.5d0 * (well%phi(isegment) + well%phi(isegment-1))
    sat_up = 0.5d0 * (well%liq%s(isegment) + well%liq%s(isegment+1))
    sat_dn = 0.5d0 * (well%liq%s(isegment) + well%liq%s(isegment-1))

    u_up = well%ql(isegment)/(sat_up*por_up)
    u_dn = well%ql(isegment-1)/(sat_dn*por_dn)

  ! ----------------------------------------BOUNDARY-FLUXES------------------
  else if (isegment == 1) then
    ! ----- bottom of well -----

    ! define face values with arithmetic averages:
    area_up = 0.5d0 * (well%area(isegment) + well%area(isegment+1))
    area_dn = well%area(isegment)
    por_up = 0.5d0 * (well%phi(isegment) + well%phi(isegment+1))
    por_dn = well%phi(isegment)
    sat_up = 0.5d0 * (well%liq%s(isegment) + well%liq%s(isegment+1))
    sat_dn = well%liq%s(isegment)

    u_up = well%ql(isegment)/(sat_up*por_up)
    u_dn = well%ql_bc(1)/(sat_dn*por_dn)        ! bottom of hole ql = ql_bc(1)

  else if (isegment == pm_well%well_grid%nsegments) then
    ! ----- top of well -----

    ! define face values with arithmetic averages:
    area_up = well%area(isegment)
    area_dn = 0.5d0 * (well%area(isegment) + well%area(isegment-1))
    por_up = well%phi(isegment)
    por_dn = 0.5d0 * (well%phi(isegment) + well%phi(isegment-1))
    sat_up = well%liq%s(isegment)
    sat_dn = 0.5d0 * (well%liq%s(isegment) + well%liq%s(isegment-1))

    u_up = well%ql_bc(2)/(sat_up*por_up)      ! top of hole ql = ql_bc(2)
    u_dn = well%ql(isegment-1)/(sat_dn*por_dn)

  endif

  ! north surface:
  J_up = (n_up*area_up)*(u_up - d_diffusion_dM)

  ! south surface:
  J_dn = (n_dn*area_dn)*(u_dn - d_diffusion_dM)

  istart = 1
  iend = pm_well%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + (J_up + J_dn)

  enddo

end subroutine PMWellJacTranFlux

! ************************************************************************** !

subroutine PMWellJacTranRxn(pm_well,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: Jblock(pm_well%nspecies,pm_well%nspecies)
  PetscInt :: isegment

  PetscInt :: istart, iend, ispecies, parent_id
  PetscReal :: vol

  ! units of Jac = [m^3-bulk/sec]
  ! units of volume = [m^3-bulk]
  ! decay_rate in [1/sec]

  vol = pm_well%well%volume(isegment)

  istart = 1
  iend = pm_well%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + &
                               (vol*pm_well%well%species_decay_rate(ispecies))

    parent_id = pm_well%well%species_parent_id(ispecies)
    if (parent_id > 0) then
      Jblock(ispecies,parent_id) = Jblock(ispecies,parent_id) - &
                        (vol*pm_well%well%species_parent_decay_rate(ispecies))
    endif

  enddo

end subroutine PMWellJacTranRxn

! ************************************************************************** !

subroutine PMWellPreSolveFlow(pm_well)
    !
    ! Author: Jennifer M. Frederick
    ! Date: 08/04/2021
    !

    implicit none

    class(pm_well_sequential_type) :: pm_well

    character(len=MAXSTRINGLENGTH) :: out_string
    PetscReal :: cur_time, cur_time_converted
    PetscReal :: dt_converted

    pm_well%flow_soln%not_converged = PETSC_TRUE
    pm_well%flow_soln%converged = PETSC_FALSE

    cur_time = pm_well%option%time + pm_well%cumulative_dt_flow
    cur_time_converted = cur_time/pm_well%output_option%tconv
    dt_converted = pm_well%dt_flow/pm_well%output_option%tconv

    if (pm_well%print_output) then
      write(out_string,'(" FLOW Step ",i6,"   Time =",1pe12.5,"   Dt =", &
                         1pe12.5," ",a4)') &
                       (pm_well%flow_soln%n_steps+1),cur_time_converted, &
                       dt_converted,pm_well%output_option%tunit
      call PrintMsg(pm_well%option,out_string)
    endif

  end subroutine PMWellPreSolveFlow

  ! ************************************************************************** !

  subroutine PMWellPreSolveTran(pm_well,master_dt)
    !
    ! Author: Jennifer M. Frederick
    ! Date: 02/22/2022
    !

    implicit none

    class(pm_well_sequential_type) :: pm_well
    PetscReal :: master_dt

    character(len=MAXSTRINGLENGTH) :: out_string
    PetscReal :: cur_time, cur_time_converted
    PetscReal :: dt_converted

    pm_well%tran_soln%not_converged = PETSC_TRUE
    pm_well%tran_soln%converged = PETSC_FALSE

    cur_time = pm_well%option%time + pm_well%cumulative_dt_tran
    pm_well%tran_soln%tran_time = cur_time

    cur_time_converted = cur_time/pm_well%output_option%tconv
    dt_converted = pm_well%dt_tran/pm_well%output_option%tconv

    write(out_string,'(" WELL TRAN Step ",i6,"   Time =",1pe12.5,"   Dt =", &
                     1pe12.5," ",a4)') &
                     (pm_well%tran_soln%n_steps+1),cur_time_converted, &
                     dt_converted,pm_well%output_option%tunit

    call PrintMsg(pm_well%option,out_string)

  end subroutine PMWellPreSolveTran

! ************************************************************************** !

subroutine PMWellUpdateSolutionFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 01/21/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscInt :: i
  PetscInt :: idof
  PetscReal, parameter :: MIN_SAT = 1.d-6
  PetscReal, parameter :: MAX_SAT = 0.99999
  do i = 1,pm_well%well_grid%nsegments
    idof = pm_well%flow_soln%ndof*(i-1)
    pm_well%well%pl(i) = pm_well%well%pl(i) +  &
                          pm_well%flow_soln%update(idof+1)
    idof = pm_well%flow_soln%ndof*i
    if (pm_well%well%gas%s(i) + pm_well%flow_soln%update(idof) < &
        MIN_SAT) then
      pm_well%flow_soln%update(idof) = MIN_SAT - pm_well%well%gas%s(i)
      pm_well%well%gas%s(i) = MIN_SAT
    elseif (pm_well%well%gas%s(i) + pm_well%flow_soln%update(idof) > &
            MAX_SAT) then
      pm_well%flow_soln%update(idof) = MAX_SAT - pm_well%well%gas%s(i)
      pm_well%well%gas%s(i) = MAX_SAT
    else
      pm_well%well%gas%s(i) = pm_well%well%gas%s(i) + &
                              pm_well%flow_soln%update(idof)
    endif

  enddo

  call pm_well%UpdateFlowProperties(PETSC_FALSE, ZERO_INTEGER)

end subroutine PMWellUpdateSolutionFlow

! ************************************************************************** !


subroutine PMWellUpdateSolutionTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/22
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscInt :: isegment
  PetscInt :: offset, istart, iend
  PetscInt :: nspecies
  PetscReal :: vol

  ! update in [mol/m3-bulk]
  ! volume in [m3-bulk]

  nspecies = pm_well%nspecies

  do isegment = 1,pm_well%well_grid%nsegments

    offset = (isegment-1)*nspecies
    istart = offset + 1
    iend = offset + nspecies

    vol = pm_well%well%volume(isegment)

    pm_well%well%aqueous_mass(1:nspecies,isegment) = &              ! [mol]
                           pm_well%well%aqueous_mass(1:nspecies,isegment) + &
                           (pm_well%tran_soln%update(istart:iend) * vol)
    pm_well%well%aqueous_conc(1:nspecies,isegment) = &
          pm_well%well%aqueous_mass(1:nspecies,isegment) / &        ! [mol]
          (pm_well%well%phi(isegment)*vol*pm_well%well%liq%s(isegment)) ! [m3-liq]

  enddo

end subroutine PMWellUpdateSolutionTran

! ************************************************************************** !

subroutine PMWellCutTimestepFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 01/24/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  ! could make this smarter or call smarter timestepping routines
  pm_well%dt_flow = pm_well%dt_flow / pm_well%flow_soln%ts_cut_factor
  pm_well%dt_flow = max(pm_well%dt_flow,pm_well%min_dt_flow)
  pm_well%well%pl = pm_well%flow_soln%prev_soln%pl
  pm_well%well%gas%s = pm_well%flow_soln%prev_soln%sg
  call pm_well%UpdateFlowProperties(PETSC_FALSE,ZERO_INTEGER)

end subroutine PMWellCutTimestepFlow

! ************************************************************************** !

subroutine PMWellCutTimestepTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  implicit none

  class(pm_well_sequential_type) :: pm_well

  pm_well%well%aqueous_mass = pm_well%tran_soln%prev_soln%aqueous_mass
  pm_well%well%aqueous_conc = pm_well%tran_soln%prev_soln%aqueous_conc
  call PMWellUpdatePropertiesTran(pm_well)

end subroutine PMWellCutTimestepTran

! ************************************************************************** !
end module WIPP_Well_class
