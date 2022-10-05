module PM_ERT_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_class
  use Option_module
  use ERT_Aux_module
  use Survey_module
  use Waypoint_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_ert_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    type(survey_type), pointer :: survey
    ! this waypoint list should only hold the survey times
    type(waypoint_list_type), pointer :: waypoint_list
    PetscInt :: linear_iterations_in_step
    PetscLogDouble :: ksp_time
    Vec :: rhs
    Vec :: dconductivity_dsaturation
    Vec :: dconductivity_dconcentration
    ! EMPIRICAL Archie and Waxman-Smits options
    PetscInt :: conductivity_mapping_law
    PetscReal :: tortuosity_constant   ! a
    PetscReal :: cementation_exponent  ! m
    PetscReal :: saturation_exponent   ! n
    PetscReal :: water_conductivity
    PetscReal :: surface_conductivity
    PetscReal :: tracer_conductivity
    PetscReal :: clay_conductivity
    PetscReal :: clay_volume_factor
    PetscReal :: max_tracer_conc
    PetscReal, pointer :: species_conductivity_coef(:)
    character(len=MAXSTRINGLENGTH) :: mobility_database
    ! Starting sulution/potential
    PetscBool :: analytical_potential
    PetscBool :: coupled_ert_flow_jacobian
  contains
    procedure, public :: Setup => PMERTSetup
    procedure, public :: ReadSimulationOptionsBlock => PMERTReadSimOptionsBlock
    procedure, public :: SetRealization => PMERTSetRealization
    procedure, public :: InitializeRun => PMERTInitializeRun
    procedure, public :: SetupSolvers => PMERTSetupSolvers
    procedure, public :: PreSolve => PMERTPreSolve
    procedure, public :: Solve => PMERTSolve
    procedure, public :: FinalizeRun => PMERTFinalizeRun
    procedure, public :: AcceptSolution => PMERTAcceptSolution
    procedure, public :: UpdateSolution => PMERTUpdateSolution
    procedure, public :: UpdateAuxVars => PMERTUpdateAuxVars
    procedure, public :: InputRecord => PMERTInputRecord
    procedure, public :: Destroy => PMERTDestroy
  end type pm_ert_type

  public :: PMERTCreate, &
            PMERTInit, &
            PMERTCast, &
            PMERTStrip

contains

! ************************************************************************** !

function PMERTCreate()
  !
  ! Creates ert process model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !
  implicit none

  class(pm_ert_type), pointer :: PMERTCreate

  class(pm_ert_type), pointer :: pm_ert

  allocate(pm_ert)
  call PMERTInit(pm_ert)
  pm_ert%name = 'Electrical Resist. Tomography'
  pm_ert%header = 'ERT GEOPHYSICS'

  PMERTCreate => pm_ert

end function PMERTCreate

! ************************************************************************** !

subroutine PMERTInit(pm_ert)
  !
  ! Initializes ert process model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !
  implicit none

  class(pm_ert_type) :: pm_ert

  call PMBaseInit(pm_ert)
  nullify(pm_ert%realization)
  nullify(pm_ert%comm1)
  nullify(pm_ert%survey)
  pm_ert%waypoint_list => WaypointListCreate()

  pm_ert%linear_iterations_in_step = 0
  pm_ert%ksp_time = 0.d0
  pm_ert%rhs = PETSC_NULL_VEC
  pm_ert%dconductivity_dsaturation = PETSC_NULL_VEC
  pm_ert%dconductivity_dconcentration = PETSC_NULL_VEC

  ! Archie and Waxman-Smits default values
  pm_ert%conductivity_mapping_law = ARCHIE
  pm_ert%tortuosity_constant = 1.d0
  pm_ert%cementation_exponent = 1.9d0
  pm_ert%saturation_exponent = 2.d0
  pm_ert%water_conductivity = 0.01d0
  pm_ert%surface_conductivity = 0.d0 ! to modifie Archie's equation
  pm_ert%tracer_conductivity = 0.d0
  pm_ert%clay_conductivity = 0.03d0
  pm_ert%clay_volume_factor = 0.0d0  ! No clay -> clean sand
  pm_ert%max_tracer_conc = UNINITIALIZED_DOUBLE

  pm_ert%analytical_potential = PETSC_TRUE
  pm_ert%coupled_ert_flow_jacobian = PETSC_FALSE

  nullify(pm_ert%species_conductivity_coef)
  pm_ert%mobility_database = ''

end subroutine PMERTInit

! ************************************************************************** !

subroutine PMERTReadSimOptionsBlock(this,input)
  !
  ! Reads input file parameters associated with the ert process model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !
  use Input_Aux_module
  use String_module
  use Units_module
  use Utility_module

  implicit none

  class(pm_ert_type) :: this
  type(input_type), pointer :: input

  type(option_type), pointer :: option
  type(waypoint_type), pointer :: waypoint
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: internal_units
  PetscReal :: units_conversion
  PetscReal, pointer :: temp_real_array(:)
  PetscInt :: temp_int
  PetscBool :: found
  PetscBool :: output_all_surveys

  option => this%option

  error_string = 'ERT Options'

  output_all_surveys = PETSC_FALSE
  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMBaseReadSimOptionsSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))
      ! Add various options for ERT if needed here
      case('COMPUTE_JACOBIAN')
        option%geophysics%compute_jacobian = PETSC_TRUE
      case('NO_ANALYTICAL_POTENTIAL')
        this%analytical_potential = PETSC_FALSE
      case('CONDUCTIVITY_MAPPING_LAW')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringToUpper(word)
        select case(trim(word))
          case('ARCHIE')
            this%conductivity_mapping_law = ARCHIE
          case('WAXMAN_SMITS')
            this%conductivity_mapping_law = WAXMAN_SMITS
          case default
            option%io_buffer  = 'CONDUCTIVITY_MAPPING_LAW: ' &
                                      // trim(word) // ' unknown.'
            call PrintErrMsg(option)
        end select
      case('TORTUOSITY_CONSTANT')
        call InputReadDouble(input,option,this%tortuosity_constant)
        call InputErrorMsg(input,option,keyword,error_string)
      case('CEMENTATION_EXPONENT')
        call InputReadDouble(input,option,this%cementation_exponent)
        call InputErrorMsg(input,option,keyword,error_string)
      case('SATURATION_EXPONENT')
        call InputReadDouble(input,option,this%saturation_exponent)
        call InputErrorMsg(input,option,keyword,error_string)
      case('WATER_CONDUCTIVITY')
        call InputReadDouble(input,option,this%water_conductivity)
        call InputErrorMsg(input,option,keyword,error_string)
      case('SURFACE_CONDUCTIVITY')
        call InputReadDouble(input,option,this%surface_conductivity)
        call InputErrorMsg(input,option,keyword,error_string)
      case('TRACER_CONDUCTIVITY')
        call InputReadDouble(input,option,this%tracer_conductivity)
        call InputErrorMsg(input,option,keyword,error_string)
      case('CLAY_CONDUCTIVITY')
        call InputReadDouble(input,option,this%clay_conductivity)
        call InputErrorMsg(input,option,keyword,error_string)
      case('CLAY_VOLUME_FACTOR','SHALE_VOLUME_FACTOR')
        call InputReadDouble(input,option,this%clay_volume_factor)
        call InputErrorMsg(input,option,keyword,error_string)
      case('SURVEY_TIMES')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'units','OUTPUT,TIMES')
        internal_units = 'sec'
        units_conversion = &
          UnitsConvertToInternal(word,internal_units,option)
        string = trim(error_string) // 'SURVEY_TIMES'
        nullify(temp_real_array)
        call UtilityReadArray(temp_real_array,NEG_ONE_INTEGER, &
                              string,input,option)
        do temp_int = 1, size(temp_real_array)
          waypoint => WaypointCreate()
          waypoint%time = temp_real_array(temp_int)*units_conversion
          waypoint%sync = PETSC_TRUE
          call WaypointInsertInList(waypoint,this%waypoint_list)
        enddo
        call DeallocateArray(temp_real_array)
      case('MOBILITY_DATABASE')
        call InputReadFilename(input,option,this%mobility_database)
        call InputErrorMsg(input,option,keyword,error_string)
      case('OUTPUT_ALL_SURVEYS')
        output_all_surveys = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (output_all_surveys) then
    waypoint => this%waypoint_list%first
    do
      if (.not.associated(waypoint)) exit
      waypoint%print_snap_output = PETSC_TRUE
      waypoint => waypoint%next
    enddo
  endif

end subroutine PMERTReadSimOptionsBlock

! ************************************************************************** !

subroutine PMERTSetup(this)
  !
  ! Initializes variables associated with ert
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !
  implicit none

  class(pm_ert_type) :: this

  ! set the communicator
  this%comm1 => this%realization%comm1
  ! setup survey
  this%survey => this%realization%survey

end subroutine PMERTSetup

! ************************************************************************** !

function PMERTCast(this)
  !
  ! Initializes a base process model to ert
  !
  ! Author: Glenn Hammond
  ! Date: 06/08/22

  use Option_module

  implicit none

  class(pm_base_type), pointer :: this

  class(pm_ert_type), pointer :: PMERTCast

  nullify(PMERTCast)
  if (.not.associated(this)) return
  select type (this)
    class is (pm_ert_type)
      PMERTCast => this
    class default
      this%option%io_buffer = 'Cannot cast pm_base_type to pm_ert_type &
        &in PMERTCast.'
      call PrintErrMsg(this%option)
  end select

end function PMERTCast

! ************************************************************************** !

subroutine PMERTSetRealization(this,realization)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !

  use Realization_Subsurface_class

  implicit none

  class(pm_ert_type) :: this
  class(realization_subsurface_type), pointer :: realization

  this%realization => realization
  this%realization_base => realization

end subroutine PMERTSetRealization

! ************************************************************************** !

recursive subroutine PMERTInitializeRun(this)
  !
  ! Initializes the time stepping
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !
  use Connection_module
  use Coupler_module
  use Discretization_module
  use Grid_module
  use Input_Aux_module
  use Patch_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use String_module
  use Transport_Constraint_RT_module
  use ZFlow_Aux_module

  implicit none

  class(pm_ert_type) :: this

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  class(reaction_rt_type), pointer :: reaction
  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_bc(:)
  type(zflow_auxvar_type), pointer :: zflow_auxvars(:,:)
  type(zflow_auxvar_type), pointer :: zflow_auxvars_bcss(:)
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(patch_type), pointer :: patch
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: ispecies
  PetscBool :: flag
  PetscReal :: tempreal
  PetscReal, parameter :: ELEMENTARY_CHARGE = 1.6022d-19 ! C
  PetscReal, parameter :: AVOGADRO_NUMBER = 6.02214d23 ! atoms per mol
  PetscInt :: iconn, sum_connection
  PetscInt :: local_id, ghosted_id
  PetscInt :: i
  PetscErrorCode :: ierr


  patch => this%realization%patch
  reaction => patch%reaction
  grid => patch%grid
  option => this%option

  call DiscretizationDuplicateVector(this%realization%discretization, &
                                     this%realization%field%work,this%rhs)

  ! Initialize to zeros
  call VecZeroEntries(this%rhs,ierr);CHKERRQ(ierr)

  if (associated(option%inversion)) then
    if (option%inversion%coupled_flow_ert .and. &
        option%inversion%calculate_ert_jacobian) then
      this%coupled_ert_flow_jacobian = PETSC_TRUE
      if (option%iflowmode == ZFLOW_MODE) then
        call DiscretizationDuplicateVector(this%realization%discretization, &
                                           this%realization%field%work, &
                                           this%dconductivity_dsaturation)
        call VecZeroEntries(this%dconductivity_dsaturation,ierr);CHKERRQ(ierr)
        if (zflow_sol_tran_eq > 0 .or. option%itranmode == RT_MODE) then
          call DiscretizationDuplicateVector(this%realization%discretization, &
                                             this%realization%field%work, &
                                             this%dconductivity_dconcentration)
          call VecZeroEntries(this%dconductivity_dconcentration, &
                              ierr);CHKERRQ(ierr)
        endif
      endif
    endif
  endif

  if (option%iflowmode == ZFLOW_MODE .and. zflow_sol_tran_eq > 0) then
    zflow_auxvars => patch%aux%ZFlow%auxvars
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      this%max_tracer_conc = max(this%max_tracer_conc, &
                                 zflow_auxvars(ZERO_INTEGER,local_id)%conc)
    enddo
    zflow_auxvars_bcss => patch%aux%ZFlow%auxvars_bc
    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0
    do
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        this%max_tracer_conc = &
          max(this%max_tracer_conc, &
              zflow_auxvars_bcss(sum_connection)%conc)
      enddo
      boundary_condition => boundary_condition%next
    enddo
    zflow_auxvars_bcss => patch%aux%ZFlow%auxvars_ss
    source_sink => patch%source_sink_list%first
    sum_connection = 0
    do
      if (.not.associated(source_sink)) exit
      cur_connection_set => source_sink%connection_set
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        this%max_tracer_conc = &
          max(this%max_tracer_conc, &
              zflow_auxvars_bcss(sum_connection)%conc)
      enddo
      source_sink => source_sink%next
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,this%max_tracer_conc,ONE_INTEGER, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                       ierr);CHKERRQ(ierr)
  endif

  if (option%itranmode == RT_MODE) then
    ! calculate species conductivity coefficients if defined
    if (len_trim(this%mobility_database) > 0) then
      if (.not.associated(reaction%primary_spec_Z)) then
        option%io_buffer = 'The CHEMISTRY block must be include a DATABASE to &
          &calculate fluid conductivity as a function of species mobilities.'
        call PrintErrMsg(option)
      endif
      input => InputCreate(IUNIT_TEMP,this%mobility_database,option)
      allocate(this%species_conductivity_coef(reaction%naqcomp))
      this%species_conductivity_coef = UNINITIALIZED_DOUBLE
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        if (len_trim(input%buf) == 0) cycle
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'MOBILITY SPECIES NAME', &
                          'MOBILITY_DATABASE')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,'MOBILITY VALUE','MOBILITY_DATABASE')
        ispecies = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE, &
                                              this%option)
        if (Initialized(ispecies)) then
          this%species_conductivity_coef(ispecies) = & ! [m^2-charge-A/V-mol]
            tempreal * &                               ! mobility [m^2/V-s]
            abs(reaction%primary_spec_Z(ispecies)) * & ! [charge/atom]
            AVOGADRO_NUMBER * &                        ! [atom/mol]
            ELEMENTARY_CHARGE                          ! [A-s] or [C]
        endif
      enddo
      flag = PETSC_FALSE
      do ispecies = 1, reaction%naqcomp
        if (Uninitialized(this%species_conductivity_coef(ispecies))) then
          if (.not.flag) then
            flag = PETSC_TRUE
            option%io_buffer = ''
            call PrintMsg(option)
          endif
          option%io_buffer = reaction%primary_species_names(ispecies)
          call PrintMsg(option)
        endif
      enddo
      if (flag) then
        option%io_buffer = 'Electrical mobilities for the species above not &
          &defined in mobility database: ' // trim(this%mobility_database)
        call PrintErrMsg(option)
      endif
      call InputDestroy(input)
    else if (associated(patch%aux%RT)) then
      rt_auxvars => patch%aux%RT%auxvars
      rt_auxvars_bc => patch%aux%RT%auxvars_bc
      ispecies = 1
      do local_id = 1, grid%nlmax  ! For each local node do...
        ghosted_id = grid%nL2G(local_id)
        if (patch%imat(ghosted_id) <= 0) cycle
        this%max_tracer_conc = max(this%max_tracer_conc, &
                                  rt_auxvars(local_id)%total(ispecies,1))
      enddo
      boundary_condition => patch%boundary_condition_list%first
      sum_connection = 0
      do
        if (.not.associated(boundary_condition)) exit
        cur_connection_set => boundary_condition%connection_set
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1
          this%max_tracer_conc = &
            max(this%max_tracer_conc, &
                rt_auxvars_bc(sum_connection)%total(ispecies,1))
        enddo
        boundary_condition => boundary_condition%next
      enddo
      source_sink => patch%source_sink_list%first
      do
        if (.not.associated(source_sink)) exit
        cur_connection_set => source_sink%connection_set
        rt_auxvar => TranConstraintRTGetAuxVar(source_sink%tran_condition% &
                                              cur_constraint_coupler)
        this%max_tracer_conc = max(this%max_tracer_conc, &
                                  rt_auxvar%total(ispecies,1))
        source_sink => source_sink%next
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,this%max_tracer_conc,ONE_INTEGER, &
                         MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                         ierr);CHKERRQ(ierr)
    endif
  endif

  ! ensure that electrodes are not placed in inactive cells
  flag = PETSC_FALSE
  do i = 1, size(this%survey%ipos_electrode)
    local_id = this%survey%ipos_electrode(i)
    if (local_id <= 0) cycle ! not on process
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      option%io_buffer = 'Electrode in inactive grid cell: ' // &
        trim(StringWrite(grid%nG2A(ghosted_id)))
      call PrintErrMsgNoStopByRank(option)
      flag = PETSC_TRUE
    endif
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,flag,ONE_INTEGER_MPI,MPI_LOGICAL,MPI_LOR, &
                     option%mycomm,ierr);CHKERRQ(ierr)
  if (flag) then
    option%io_buffer = 'Electrodes in inactive cells (see above).'
    call PrintErrMsg(option)
  endif

  ! calculate initial electrical conductivity
  call PMERTPreSolve(this)

end subroutine PMERTInitializeRun

! ************************************************************************** !

subroutine PMERTDummyExtension(this)
  !
  ! Dummy routine to supercede requirement for extension from base
  !
  ! Author: Glenn Hammond
  ! Date: 03/26/21
  !
  implicit none

  class(pm_ert_type) :: this

end subroutine PMERTDummyExtension

! ************************************************************************** !

subroutine PMERTSetupSolvers(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/29/21
  !

  use Solver_module
  use Discretization_module

  implicit none

  class(pm_ert_type) :: this

  type(option_type), pointer :: option
  type(solver_type), pointer :: solver

  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string

  option => this%option
  solver => this%solver

  call SolverCreateKSP(solver,option%mycomm)

  call PrintMsg(option,"  Beginning setup of ERT KSP")
  ! TODO(pj): set ert as prefix
  call KSPSetOptionsPrefix(solver%ksp,"geop_",ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)

  solver%M_mat_type = MATAIJ
  solver%Mpre_mat_type = MATAIJ
  call DiscretizationCreateMatrix(this%realization%discretization, &
                                  ONEDOF, &
                                  solver%Mpre_mat_type, &
                                  solver%Mpre,option)

  call MatSetOptionsPrefix(solver%Mpre,"geop_",ierr);CHKERRQ(ierr)
  solver%M = solver%Mpre

  ! Have PETSc do a KSP_View() at the end of each solve if
  ! verbosity > 0.
  if (option%verbosity >= 2) then
    string = '-geop_ksp_view'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS,string, &
                                  ierr);CHKERRQ(ierr)
    string = '-geop_ksp_monitor'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS,string, &
                                  ierr);CHKERRQ(ierr)
  endif

  call PrintMsg(option,"  Finished setting up ERT KSP")

  ! TODO(pj): Whay do I need the follwing call as other pmc don't need?
  call KSPSetOperators(solver%ksp,solver%M,solver%Mpre,ierr);CHKERRQ(ierr)

  call SolverSetKSPOptions(solver,option)

end subroutine PMERTSetupSolvers

! ************************************************************************** !

subroutine PMERTPreSolve(this)

  ! Update flow and transport dependent variables (e.g. bulk conductivity)
  ! prior to solve.
  !
  ! Author: Glenn Hammond
  ! Date: 03/29/21
  !
  use Field_module
  use Global_Aux_module
  use Grid_module
  use Material_Aux_module
  use Option_module
  use Patch_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Realization_Base_class
  use Variables_module
  use ERT_module
  use ZFlow_Aux_module

  implicit none

  class(pm_ert_type) :: this

  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(zflow_auxvar_type), pointer :: zflow_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  class(reaction_rt_type), pointer :: reaction
  PetscInt :: ghosted_id,local_id
  PetscInt :: species_id
  PetscInt :: empirical_law
  PetscReal :: a,m,n,cond_w,cond_s,cond_c,Vc,cond  ! variables for Archie's law
  PetscReal :: por,sat
  PetscReal :: cond_sp,cond_w0
  PetscReal :: dcond_dsat,dcond_dconc
  PetscReal :: tracer_scale
  PetscReal, pointer :: dcond_dsat_vec_ptr(:),dcond_dconc_vec_ptr(:)
  PetscErrorCode :: ierr

  option => this%option
  if (option%iflowmode == NULL_MODE .and. option%itranmode == NULL_MODE) return

  patch => this%realization%patch
  grid => patch%grid
  reaction => patch%reaction

  empirical_law = this%conductivity_mapping_law
  a = this%tortuosity_constant
  m = this%cementation_exponent
  n = this%saturation_exponent
  Vc = this%clay_volume_factor
  cond_w = this%water_conductivity
  cond_s = this%surface_conductivity
  cond_c = this%clay_conductivity
  cond_w0 = cond_w

  global_auxvars => patch%aux%Global%auxvars
  nullify(rt_auxvars)
  nullify(zflow_auxvars)
  if (option%itranmode == RT_MODE) then
    rt_auxvars => patch%aux%RT%auxvars
  endif
  if (zflow_sol_tran_eq > 0) then
    zflow_auxvars => patch%aux%ZFlow%auxvars
  endif
  if (Initialized(this%max_tracer_conc)) then
    tracer_scale = this%tracer_conductivity/this%max_tracer_conc
  endif

  if (this%coupled_ert_flow_jacobian) then
    call VecGetArrayF90(this%dconductivity_dsaturation, &
                        dcond_dsat_vec_ptr,ierr);CHKERRQ(ierr)
    if (associated(rt_auxvars) .or. associated(zflow_auxvars)) then
      call VecGetArrayF90(this%dconductivity_dconcentration, &
                          dcond_dconc_vec_ptr,ierr);CHKERRQ(ierr)
    endif
  endif

  material_auxvars => patch%aux%Material%auxvars
  do ghosted_id = 1, grid%ngmax
    local_id = grid%nG2L(ghosted_id)
    dcond_dsat = 0.d0
    dcond_dconc = 0.d0
    if (patch%imat(ghosted_id) <= 0) cycle
    por = material_auxvars(ghosted_id)%porosity
    sat = global_auxvars(ghosted_id)%sat(1)
    if (associated(rt_auxvars)) then
      if (associated(this%species_conductivity_coef)) then
        ! assuming that we sum conductivity across species
        cond_sp = 0.d0
        do species_id = 1, reaction%naqcomp
          cond_sp = cond_sp + &                          ! S/m
            this%species_conductivity_coef(species_id)  * &![m^2-charge-A/V-mol]
            rt_auxvars(ghosted_id)%pri_molal(species_id)* &![mol/kg water]
            global_auxvars(ghosted_id)%den_kg(1)           ![kg water/m^3]
        enddo
        ! modify fluid conductivity for species contribution
        cond_w = cond_w0 + cond_sp
      else
        species_id = 1
        cond_sp = tracer_scale * rt_auxvars(ghosted_id)%total(species_id,1)
        cond_w = cond_w0 + cond_sp
        dcond_dconc = tracer_scale
      endif
    endif
    if (associated(zflow_auxvars)) then
      cond_sp = tracer_scale * zflow_auxvars(ZERO_INTEGER,ghosted_id)%conc
      cond_w = cond_w0 + cond_sp
      dcond_dconc = tracer_scale
    endif
    ! compute conductivity
    call ERTConductivityFromEmpiricalEqs(por,sat,a,m,n,Vc,cond_w,cond_s, &
                                         cond_c,empirical_law,cond,dcond_dsat)
    material_auxvars(ghosted_id)%electrical_conductivity(1) = cond
    if (this%coupled_ert_flow_jacobian) then
      if (local_id > 0) dcond_dsat_vec_ptr(local_id) = dcond_dsat
      if (associated(rt_auxvars) .or. associated(zflow_auxvars)) then
        if (local_id > 0) dcond_dconc_vec_ptr(local_id) = dcond_dconc
      endif
    endif
  enddo

  if (this%coupled_ert_flow_jacobian) then
    call VecRestoreArrayF90(this%dconductivity_dsaturation, &
                            dcond_dsat_vec_ptr,ierr);CHKERRQ(ierr)
    if (associated(rt_auxvars) .or. associated(zflow_auxvars)) then
      call VecRestoreArrayF90(this%dconductivity_dconcentration, &
                             dcond_dconc_vec_ptr,ierr);CHKERRQ(ierr)
    endif
  endif

end subroutine PMERTPreSolve

! ************************************************************************** !

subroutine PMERTSolve(this,time,ierr)

  ! Solves the linear systsem for ERT for all electrodes
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/27/21
  !
  use Patch_module
  use Debug_module
  use Grid_module
  use Solver_module
  use Field_module
  use Discretization_module
  use ERT_module
  use String_module

  implicit none

  class(pm_ert_type) :: this
  PetscReal :: time
  PetscErrorCode :: solve_ierr

  class(realization_subsurface_type), pointer :: realization
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(solver_type), pointer :: solver
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(survey_type), pointer :: survey
  type(ert_auxvar_type), pointer :: ert_auxvars(:)

  PetscInt :: ielec,nelec
  PetscInt :: elec_id
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: num_linear_iterations
  PetscReal :: val
  PetscReal :: average_cond
  PetscReal, pointer :: vec_ptr(:)
  character(len=MAXSTRINGLENGTH) :: string

  PetscViewer :: viewer
  !PetscLogDouble :: log_start_time
  PetscLogDouble :: log_ksp_start_time
  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscErrorCode :: ierr
  KSPConvergedReason :: ksp_reason

  solve_ierr = 0
  ! Forward solve start
  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

  solver => this%solver
  survey => this%survey
  realization => this%realization
  field => realization%field
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid

  ert_auxvars => patch%aux%ERT%auxvars

  ! Build System matrix
  call ERTCalculateMatrix(realization,solver%M, &
                          this%option%geophysics%compute_jacobian)
  call KSPSetOperators(solver%ksp,solver%M,solver%M,ierr);CHKERRQ(ierr)
  !call MatView(solver%M,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)

  ! Get Average Conductivity for a 3D model
  call ERTCalculateAverageConductivity(realization)
  average_cond = survey%average_conductivity

  nelec = survey%num_electrode
  this%linear_iterations_in_step = 0

  this%option%io_buffer = '  Solving for electrode:'
  call PrintMsgNoAdvance(this%option)
  do ielec=1,nelec
    write(this%option%io_buffer,'(x,a)') trim(StringWrite(ielec))
    call PrintMsgNoAdvance(this%option)
    if (this%analytical_potential) then
      ! Initial Solution -> analytic sol for a half-space
      ! Get Analytical potential for a half-space
      call ERTCalculateAnalyticPotential(realization,ielec,average_cond)
      ! assign analytic potential as initial solution
      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (patch%imat(ghosted_id) <= 0) cycle
        vec_ptr(local_id) = ert_auxvars(ghosted_id)%potential(ielec)
      enddo
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      call KSPSetInitialGuessNonzero(solver%ksp,PETSC_TRUE,ierr);CHKERRQ(ierr)
    else
      ! zero initial solution
      call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
    endif

    ! RHS
    call VecZeroEntries(this%rhs,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%rhs,vec_ptr,ierr);CHKERRQ(ierr)

    ! Get the local-id of ielec
    elec_id = survey%ipos_electrode(ielec)

    if (elec_id > 0) then
      ! DBG
      !print*,'Source cell-id: ',elec_id,grid%x(grid%nL2G(elec_id)), &
      !         grid%y(grid%nL2G(elec_id)),grid%z(grid%nL2G(elec_id))
      ! it should qualify on only one proc
      val = -1.0
      vec_ptr(elec_id) = val
    endif
    call VecRestoreArrayF90(this%rhs,vec_ptr,ierr);CHKERRQ(ierr)

    if (realization%debug%vecview_residual) then
      string = 'ERTrhs_' // trim(adjustl(StringWrite(elec_id)))
      call DebugCreateViewer(realization%debug,string,this%option,viewer)
      call VecView(this%rhs,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

    ! Solve system
    call PetscTime(log_ksp_start_time,ierr);CHKERRQ(ierr)
    call KSPSolve(solver%ksp,this%rhs,field%work,ierr);CHKERRQ(ierr)
    call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
    this%ksp_time = this%ksp_time + (log_end_time - log_ksp_start_time)

    if (realization%debug%vecview_solution) then
      string = 'ERTsolution_' // trim(adjustl(StringWrite(elec_id)))
      call DebugCreateViewer(realization%debug,string,this%option,viewer)
      call VecView(field%work,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

    call DiscretizationGlobalToLocal(discretization,field%work, &
                                     field%work_loc,ONEDOF)
    call VecGetArrayF90(field%work_loc,vec_ptr,ierr);CHKERRQ(ierr)
    ! store potentials for each electrode
    do ghosted_id=1,grid%ngmax
      if (patch%imat(ghosted_id) <= 0) cycle
      ert_auxvars(ghosted_id)%potential(ielec) = vec_ptr(ghosted_id)
    enddo
    call VecRestoreArrayF90(field%work_loc,vec_ptr,ierr);CHKERRQ(ierr)

    call KSPGetIterationNumber(solver%ksp,num_linear_iterations, &
                               ierr);CHKERRQ(ierr)
    call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr);CHKERRQ(ierr)
    this%linear_iterations_in_step = this%linear_iterations_in_step + &
                                     num_linear_iterations
  enddo
  call PrintMsg(this%option,'')

  ! Assemble solutions
  call PMERTAssembleSimulatedData(this,time)

  ! Build Jacobian
  if (this%option%geophysics%compute_jacobian) then
    if (associated(this%option%inversion)) then
      if (.not.this%option%inversion%coupled_flow_ert) then
        call PMERTBuildJacobian(this)
      elseif (this%option%inversion%coupled_flow_ert .and. &
              this%option%inversion%calculate_ert_jacobian) then
        call PMERTBuildCoupledJacobian(this)
      endif
    else
      call PMERTBuildJacobian(this)
    endif
  endif

end subroutine PMERTSolve

! ************************************************************************** !

subroutine PMERTAssembleSimulatedData(this,time)
  !
  ! Assembles ERT simulated data for each measurement
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/10/21
  !

  use Patch_module
  use Grid_module

  implicit none

  class(pm_ert_type) :: this
  PetscReal :: time

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(survey_type), pointer :: survey
  type(ert_auxvar_type), pointer :: ert_auxvars(:)


  PetscInt :: idata
  PetscInt :: ia,ib,im,in
  PetscInt :: local_id_m,local_id_n
  PetscInt :: ghosted_id_m,ghosted_id_n
  character(len=MAXWORDLENGTH) :: time_suffix
  PetscErrorCode :: ierr

  option => this%option
  patch => this%realization%patch
  grid => patch%grid
  survey => this%survey

  ert_auxvars => patch%aux%ERT%auxvars

  survey%dsim = 0.

  do idata=1,survey%num_measurement
    ! for A and B electrodes
    ia = survey%config(1,idata)
    ib = survey%config(2,idata)
    im = survey%config(3,idata)
    in = survey%config(4,idata)

    if (im /= 0) then
      local_id_m = survey%ipos_electrode(im)
      if (local_id_m > 0) then
        ghosted_id_m = grid%nL2G(local_id_m)
        ! Due to source A at +ve M
        if (ia /= 0) survey%dsim(idata) = survey%dsim(idata) + &
                             ert_auxvars(ghosted_id_m)%potential(ia)
        ! Due to sink B at +ve M
        if (ib /= 0) survey%dsim(idata) = survey%dsim(idata) - &
                             ert_auxvars(ghosted_id_m)%potential(ib)
      endif
    endif

    if (in /= 0) then
      local_id_n = survey%ipos_electrode(in)
      if (local_id_n > 0) then
        ghosted_id_n = grid%nL2G(local_id_n)
        ! Due to source at A at -ve N
        if (ia /= 0) survey%dsim(idata) = survey%dsim(idata) - &
                             ert_auxvars(ghosted_id_n)%potential(ia)
        ! Due to sink at B at -ve N
        if (ib /= 0) survey%dsim(idata) = survey%dsim(idata) + &
                                    ert_auxvars(ghosted_id_n)%potential(ib)
      endif
    endif

    ! Reduce/allreduce?
    call MPI_Allreduce(MPI_IN_PLACE,survey%dsim(idata),ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                       ierr);CHKERRQ(ierr)

  enddo

  ! write simuated data in a E4D .srv file
  if (OptionIsIORank(option)) then
    write(time_suffix,'(f15.4)') time/this%output_option%tconv
    time_suffix = trim(adjustl(time_suffix)) // &
                  trim(adjustl(this%output_option%tunit))
    call SurveyWriteERT(this%survey,time_suffix,option)
  endif

end subroutine PMERTAssembleSimulatedData

! ************************************************************************** !

subroutine PMERTBuildJacobian(this)
  !
  ! Builds ERT Jacobian Matrix distributed across processors
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 03/15/21
  !

  use Patch_module
  use Grid_module
  use Inversion_Coupled_Aux_module
  use Inversion_Parameter_module
  use Inversion_TS_Aux_module
  use Material_Aux_module
  use String_module
  use Timer_class
  use Utility_module

  implicit none

  class(pm_ert_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(survey_type), pointer :: survey
  type(ert_auxvar_type), pointer :: ert_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  class(timer_type), pointer ::timer
  type(inversion_coupled_soln_type), pointer :: solutions(:)
  type(inversion_parameter_type), pointer :: parameters(:)

  PetscInt, pointer :: cell_neighbors(:,:)
  PetscReal, allocatable :: phi_sor(:), phi_rec(:)
  PetscReal :: jacob
  PetscReal :: cond,wd,wd_cull
  PetscInt :: idata
  PetscInt :: ia,ib,im,in
  PetscInt :: local_id,ghosted_id
  PetscInt :: inbr,num_neighbors
  PetscErrorCode :: ierr

  option => this%option
  patch => this%realization%patch
  grid => patch%grid
  survey => this%survey

  ert_auxvars => patch%aux%ERT%auxvars
  material_auxvars => patch%aux%Material%auxvars
  cell_neighbors => grid%cell_neighbors_local_ghosted

  call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
  timer => TimerCreate()
  call timer%Start()

  if (OptionPrintToScreen(this%option)) then
    write(*,'(/," --> Building ERT Jacobian matrix:")')
  endif

  do idata=1,survey%num_measurement

    ! for A and B electrodes
    ia = survey%config(1,idata)
    ib = survey%config(2,idata)
    im = survey%config(3,idata)
    in = survey%config(4,idata)

    do local_id=1,grid%nlmax

      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      num_neighbors = cell_neighbors(0,local_id)
      allocate(phi_sor(num_neighbors+1), phi_rec(num_neighbors+1))
      phi_sor = 0.d0
      phi_rec = 0.d0

      ! Source electrode +A
      if(ia/=0) phi_sor(1) = phi_sor(1) + ert_auxvars(ghosted_id)%potential(ia)
      ! Source electrode -B
      if(ib/=0) phi_sor(1) = phi_sor(1) - ert_auxvars(ghosted_id)%potential(ib)
      ! Receiver electrode +M
      if(im/=0) phi_rec(1) = phi_rec(1) + ert_auxvars(ghosted_id)%potential(im)
      ! Receiver electrode -N
      if(in/=0) phi_rec(1) = phi_rec(1) - ert_auxvars(ghosted_id)%potential(in)

      jacob = phi_sor(1) * ert_auxvars(ghosted_id)%delM(1) * phi_rec(1)

      do inbr = 1,num_neighbors
        ! Source electrode +A
        if (ia/=0) then
          phi_sor(inbr+1) = phi_sor(inbr+1) +                               &
               ert_auxvars(abs(cell_neighbors(inbr,local_id)))%potential(ia)
        endif

        ! Source electrode -B
        if (ib/=0) then
          phi_sor(inbr+1) = phi_sor(inbr+1) -                               &
               ert_auxvars(abs(cell_neighbors(inbr,local_id)))%potential(ib)
        endif

        ! Receiver electrode +M
        if (im/=0) then
          phi_rec(inbr+1) = phi_rec(inbr+1) +                               &
               ert_auxvars(abs(cell_neighbors(inbr,local_id)))%potential(im)
        endif

        ! Receiver electrode -N
        if (in/=0) then
          phi_rec(inbr+1) = phi_rec(inbr+1) -                               &
               ert_auxvars(abs(cell_neighbors(inbr,local_id)))%potential(in)
        endif

        jacob = jacob +                                                     &
          phi_sor(1)*ert_auxvars(ghosted_id)%delM(1+inbr)*phi_rec(inbr+1) + &
          phi_sor(1+inbr) * (                                               &
                      ert_auxvars(ghosted_id)%delM(1+inbr)*phi_rec(1) -     &
                      ert_auxvars(ghosted_id)%delM(1+inbr)*phi_rec(1+inbr) )

      enddo

      cond = material_auxvars(ghosted_id)%electrical_conductivity(1)

      ! As phi_rec is due to -ve unit source but A^-1(p) gives field due to
      ! +ve unit source so phi_rec -> - phi_rec
      ! thus => jacob = phi_s * (dM/dcond) * phi_r
      ! wrt m=ln(cond) -> dV/dm = cond * dV/dcond
      wd = survey%Wd(idata)
      wd_cull = survey%Wd_cull(idata)
      ert_auxvars(ghosted_id)%jacobian(idata) = jacob * cond * wd * wd_cull

      deallocate(phi_sor, phi_rec)
    enddo
  enddo

  ! I can now deallocate potential and delM (and M just after solving)
  ! But what about potential field output?

  call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime())) &
    // ' seconds to build Jacobian.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine PMERTBuildJacobian

! ************************************************************************** !

subroutine PMERTBuildCoupledJacobian(this)
  !
  ! Builds ERT FLOW Coupled Jacobian Matrix distributed across processors
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/03/22
  !

  use Patch_module
  use Grid_module
  use Inversion_Coupled_Aux_module
  use Inversion_Parameter_module
  use Inversion_TS_Aux_module
  use Material_Aux_module
  use String_module
  use Timer_class
  use Utility_module
  use ZFlow_Aux_module

  implicit none

  class(pm_ert_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(survey_type), pointer :: survey
  type(ert_auxvar_type), pointer :: ert_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  class(timer_type), pointer ::timer
  type(inversion_coupled_soln_type), pointer :: solutions(:)
  type(inversion_parameter_type), pointer :: parameters(:)

  PetscInt, pointer :: cell_neighbors(:,:)
  PetscReal, allocatable :: phi_sor(:), phi_rec(:)
  PetscReal, pointer :: dsat_dparam_ptr(:),dconc_dparam_ptr(:)
  PetscReal, pointer :: dcond_dsat_vec_ptr(:),dcond_dconc_vec_ptr(:)
  PetscReal :: jacob,coupled_jacob
  PetscReal :: cond,wd,wd_cull
  PetscInt :: idata,ndata,imeasurement
  PetscInt :: ielec
  PetscInt :: ia,ib,im,in,isurvey,iparam
  PetscInt :: local_id,ghosted_id
  PetscInt :: local_id_a,local_id_b
  PetscInt :: ghosted_id_a,ghosted_id_b
  PetscInt :: local_id_m,local_id_n
  PetscInt :: ghosted_id_m,ghosted_id_n
  PetscInt :: inbr,num_neighbors
  PetscErrorCode :: ierr

  option => this%option
  patch => this%realization%patch
  grid => patch%grid
  survey => this%survey

  ert_auxvars => patch%aux%ERT%auxvars
  material_auxvars => patch%aux%Material%auxvars
  cell_neighbors => grid%cell_neighbors_local_ghosted

  call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
  timer => TimerCreate()
  call timer%Start()

  if (OptionPrintToScreen(this%option)) then
    write(*,'(/," --> Building partial ERT FLOW Coupled Jacobian matrix:")')
  endif

  ndata = survey%num_measurement

  if (this%coupled_ert_flow_jacobian) then
    solutions => &
      patch%aux%inversion_forward_aux%inversion_coupled_aux%solutions
    parameters => &
      patch%aux%inversion_forward_aux%inversion_coupled_aux%parameters

    call VecGetArrayF90(this%dconductivity_dsaturation, &
                        dcond_dsat_vec_ptr,ierr);CHKERRQ(ierr)
    if (Initialized(zflow_sol_tran_eq)) then
      call VecGetArrayF90(this%dconductivity_dconcentration, &
                          dcond_dconc_vec_ptr,ierr);CHKERRQ(ierr)
    endif

    do isurvey = 1, size(solutions)
      if (.not.Equal(solutions(isurvey)%time,option%time)) cycle
      do iparam = 1, size(parameters)
        call VecGetArrayF90(solutions(isurvey)% &
                              dsaturation_dparameter(iparam), &
                            dsat_dparam_ptr,ierr);CHKERRQ(ierr)
        if (Initialized(zflow_sol_tran_eq)) then
          call VecGetArrayF90(solutions(isurvey)% &
                                dsolute_dparameter(iparam), &
                              dconc_dparam_ptr,ierr);CHKERRQ(ierr)
        endif

        do idata=1,ndata

          ! for A and B electrodes
          ia = survey%config(1,idata)
          ib = survey%config(2,idata)
          im = survey%config(3,idata)
          in = survey%config(4,idata)

          imeasurement = idata + (isurvey - 1) * ndata

          do local_id=1,grid%nlmax

            ghosted_id = grid%nL2G(local_id)
            if (patch%imat(ghosted_id) <= 0) cycle

            num_neighbors = cell_neighbors(0,local_id)
            allocate(phi_sor(num_neighbors+1), phi_rec(num_neighbors+1))
            phi_sor = 0.d0
            phi_rec = 0.d0

            ! Source electrode +A
            if(ia/=0) phi_sor(1) = phi_sor(1) + &
                                   ert_auxvars(ghosted_id)%potential(ia)
            ! Source electrode -B
            if(ib/=0) phi_sor(1) = phi_sor(1) - &
                                   ert_auxvars(ghosted_id)%potential(ib)
            ! Receiver electrode +M
            if(im/=0) phi_rec(1) = phi_rec(1) + &
                                   ert_auxvars(ghosted_id)%potential(im)
            ! Receiver electrode -N
            if(in/=0) phi_rec(1) = phi_rec(1) - &
                                   ert_auxvars(ghosted_id)%potential(in)

            jacob = phi_sor(1) * ert_auxvars(ghosted_id)%delM(1) * phi_rec(1)

            do inbr = 1,num_neighbors
            ! Source electrode +A
              if (ia/=0) then
                phi_sor(inbr+1) = phi_sor(inbr+1) + &
                  ert_auxvars(abs(cell_neighbors(inbr,local_id)))%potential(ia)
              endif

              ! Source electrode -B
              if (ib/=0) then
                phi_sor(inbr+1) = phi_sor(inbr+1) - &
                  ert_auxvars(abs(cell_neighbors(inbr,local_id)))%potential(ib)
              endif

              ! Receiver electrode +M
              if (im/=0) then
                phi_rec(inbr+1) = phi_rec(inbr+1) + &
                  ert_auxvars(abs(cell_neighbors(inbr,local_id)))%potential(im)
              endif

              ! Receiver electrode -N
              if (in/=0) then
                phi_rec(inbr+1) = phi_rec(inbr+1) - &
                  ert_auxvars(abs(cell_neighbors(inbr,local_id)))%potential(in)
              endif

              jacob = jacob + &
                      phi_sor(1) * ert_auxvars(ghosted_id)%delM(1+inbr) * &
                        phi_rec(inbr+1) + &
                      phi_sor(1+inbr) * ( &
                        ert_auxvars(ghosted_id)%delM(1+inbr)*phi_rec(1) - &
                        ert_auxvars(ghosted_id)%delM(1+inbr)*phi_rec(1+inbr) )

            enddo

            ! dERT/dparam = dERT/dCond * dCond/dSat  * dSat/dParam + &
            !               dERT/dCond * dCond/dConc * dConc/dParam
            coupled_jacob = jacob * dsat_dparam_ptr(local_id) * &
                            dcond_dsat_vec_ptr(local_id)
            if (Initialized(zflow_sol_tran_eq)) then
              coupled_jacob = coupled_jacob + &
                              jacob * dconc_dparam_ptr(local_id) * &
                              dcond_dconc_vec_ptr(local_id)
            endif

            call MatSetValue(patch%aux%inversion_forward_aux% &
                               JsensitivityT_ptr, &
                             iparam-1,imeasurement-1,coupled_jacob, &
                             ADD_VALUES,ierr);CHKERRQ(ierr)

            deallocate(phi_sor, phi_rec)

          enddo
        enddo
        call VecRestoreArrayF90(solutions(isurvey)% &
                              dsaturation_dparameter(iparam), &
                            dsat_dparam_ptr,ierr);CHKERRQ(ierr)
        if (Initialized(zflow_sol_tran_eq)) then
          call VecRestoreArrayF90(solutions(isurvey)% &
                              dsolute_dparameter(iparam), &
                            dconc_dparam_ptr,ierr);CHKERRQ(ierr)
        endif
      enddo
    enddo

    call VecRestoreArrayF90(this%dconductivity_dsaturation, &
                            dcond_dsat_vec_ptr,ierr);CHKERRQ(ierr)
    if (Initialized(zflow_sol_tran_eq)) then
      call VecRestoreArrayF90(this%dconductivity_dconcentration, &
                              dcond_dconc_vec_ptr,ierr);CHKERRQ(ierr)
    endif

  endif

  call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime())) &
    // ' seconds to build partial coupled Jacobian.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine PMERTBuildCoupledJacobian

! ************************************************************************** !

function PMERTAcceptSolution(this)
  !
  ! PMERTAcceptSolution:
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !

  implicit none

  class(pm_ert_type) :: this

  PetscBool :: PMERTAcceptSolution

  ! do nothing
  PMERTAcceptSolution = PETSC_TRUE

end function PMERTAcceptSolution

! ************************************************************************** !

recursive subroutine PMERTFinalizeRun(this)
  !
  ! Finalizes the time stepping
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !

  implicit none

  class(pm_ert_type) :: this

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMERTFinalizeRun

! ************************************************************************** !

subroutine PMERTUpdateSolution(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !

  implicit none

  class(pm_ert_type) :: this

end subroutine PMERTUpdateSolution

! ************************************************************************** !

subroutine PMERTUpdateAuxVars(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21


  implicit none

  class(pm_ert_type) :: this

end subroutine PMERTUpdateAuxVars

! ************************************************************************** !

subroutine PMERTInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !

  implicit none

  class(pm_ert_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMERTInputRecord

! ************************************************************************** !

subroutine PMERTStrip(this)
  !
  ! Strips members of ERT process model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !
  use Utility_module

  implicit none

  class(pm_ert_type) :: this

  PetscErrorCode :: ierr

  call PMBaseDestroy(this)
  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%survey)
  call WaypointListDestroy(this%waypoint_list)

  call DeallocateArray(this%species_conductivity_coef)
  if (this%rhs /= PETSC_NULL_VEC) then
    call VecDestroy(this%rhs,ierr);CHKERRQ(ierr)
  endif
!  if (this%dconductivity_dsaturation /= PETSC_NULL_VEC) then
!    call VecDestroy(this%dconductivity_dsaturation,ierr);CHKERRQ(ierr)
!  endif
  if (this%dconductivity_dconcentration /= PETSC_NULL_VEC) then
   call VecDestroy(this%dconductivity_dconcentration,ierr);CHKERRQ(ierr)
  endif

end subroutine PMERTStrip

! ************************************************************************** !

subroutine PMERTDestroy(this)
  !
  ! Destroys ERT process model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !
  implicit none

  class(pm_ert_type) :: this

  call PMERTStrip(this)

end subroutine PMERTDestroy

end module PM_ERT_class
