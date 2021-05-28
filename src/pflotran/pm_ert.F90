module PM_ERT_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_module
  use Option_module
  use ERT_Aux_module
  use Survey_module
  use Inversion_module
  use Waypoint_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_ert_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    type(survey_type), pointer :: survey
    type(inversion_type), pointer :: inversion
    ! this waypoint list should only hold the survey times
    type(waypoint_list_type), pointer :: waypoint_list
    PetscInt :: linear_iterations_in_step
    PetscLogDouble :: ksp_time
    Vec :: rhs
    ! EMPIRICAL Archie and Waxman-Smits options
    PetscReal :: tortuosity_constant   ! a
    PetscReal :: cementation_exponent  ! m
    PetscReal :: saturation_exponent   ! n
    PetscReal :: water_conductivity
    PetscReal :: tracer_conductivity
    PetscReal :: clay_conductivity
    PetscReal :: clay_volume_factor
    PetscReal :: max_tracer_conc
    PetscReal, pointer :: species_conductivity_coef(:)
    character(len=MAXSTRINGLENGTH) :: mobility_database
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
  pm_ert%name = 'Electrical Resistivity Tomography'
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
  nullify(pm_ert%inversion)
  pm_ert%waypoint_list => WaypointListCreate()

  pm_ert%linear_iterations_in_step = 0
  pm_ert%ksp_time = 0.d0
  pm_ert%rhs = PETSC_NULL_VEC

  ! Archie and Waxman-Smits default values
  pm_ert%tortuosity_constant = 1.d0
  pm_ert%cementation_exponent = 1.9d0
  pm_ert%saturation_exponent = 2.d0
  pm_ert%water_conductivity = 0.01d0
  pm_ert%tracer_conductivity = 0.d0
  pm_ert%clay_conductivity = 0.03d0
  pm_ert%clay_volume_factor = 0.0d0  ! No clay -> clean sand
  pm_ert%max_tracer_conc = UNINITIALIZED_DOUBLE

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
      case('INVERSION')
        option%geophysics%inversion = PETSC_TRUE
        option%geophysics%compute_jacobian = PETSC_TRUE
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

  type(ert_type), pointer :: ert

  ! set the communicator
  this%comm1 => this%realization%comm1
  ! setup survey
  this%survey => this%realization%survey
  ! setup inversion
  if (associated(this%realization%inversion)) then
    this%inversion => this%realization%inversion
  endif

end subroutine PMERTSetup

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
  use Transport_Constraint_RT_module

  implicit none

  class(pm_ert_type) :: this

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  class(reaction_rt_type), pointer :: reaction
  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_bc(:)
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(patch_type), pointer :: patch
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: ispecies
  PetscBool :: flag
  PetscReal :: tempreal
  PetscReal, parameter :: ELEMENTARY_CHARGE = 1.6022d-19 ! C
  PetscReal, parameter :: AVOGADRO_NUMBER = 6.02214d23 ! atoms per mol
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: iconn, sum_connection
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr


  patch => this%realization%patch
  reaction => patch%reaction
  grid => patch%grid
  option => this%option

  call DiscretizationDuplicateVector(this%realization%discretization, &
                                     this%realization%field%work,this%rhs)

  ! Initialize to zeros
  call VecZeroEntries(this%rhs,ierr);CHKERRQ(ierr)

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
  endif

  if (option%geophysics%inversion) then
    if (.not.associated(this%inversion)) then
      option%io_buffer = 'There should be an INVERSION card in input file &
                          &if process model ERT has INVERSION option.'
      call PrintErrMsg(option)
    else
      call InversionConstrainedArraysFromList(this%inversion,patch,option)
    endif
  endif

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
  call KSPSetOptionsPrefix(solver%ksp, "geop_",ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)

  solver%M_mat_type = MATAIJ
  solver%Mpre_mat_type = MATAIJ
  !TODO(geh): XXXCreateJacobian -> XXXCreateMatrix
  call DiscretizationCreateJacobian(this%realization%discretization, &
                                    ONEDOF, &
                                    solver%Mpre_mat_type, &
                                    solver%Mpre,option)

  call MatSetOptionsPrefix(solver%Mpre,"geop_",ierr);CHKERRQ(ierr)
  solver%M = solver%Mpre

  ! Have PETSc do a KSP_View() at the end of each solve if
  ! verbosity > 0.
  if (option%verbosity >= 2) then
    string = '-geop_ksp_view'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                  string, ierr);CHKERRQ(ierr)
    string = '-geop_ksp_monitor'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                  string, ierr);CHKERRQ(ierr)
  endif

  call PrintMsg(option,"  Finished setting up ERT KSP")

  ! TODO(pj): Whay do I need the follwing call as other pmc don't need?
  call  KSPSetOperators(solver%ksp,solver%M,solver%Mpre, ierr)

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
  use Material_Aux_class
  use Option_module
  use Patch_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Realization_Base_class
  use Variables_module
  use ERT_module

  implicit none

  class(pm_ert_type) :: this

  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  class(reaction_rt_type), pointer :: reaction
  PetscInt :: ghosted_id
  PetscInt :: species_id
  PetscReal :: a,m,n,cond_w,cond_c,Vc,cond  ! variables for Archie's law
  PetscReal :: por,sat
  PetscReal :: cond_sp,cond_w0
  PetscReal :: tracer_scale
  PetscErrorCode :: ierr

  option => this%option
  if (option%iflowmode == NULL_MODE .and. option%itranmode == NULL_MODE) return

  patch => this%realization%patch
  grid => patch%grid
  reaction => patch%reaction

  a = this%tortuosity_constant
  m = this%cementation_exponent
  n = this%saturation_exponent
  Vc = this%clay_volume_factor
  cond_w = this%water_conductivity
  cond_c = this%clay_conductivity
  cond_w0 = cond_w

  global_auxvars => patch%aux%Global%auxvars
  nullify(rt_auxvars)
  if (associated(patch%aux%RT)) then
    rt_auxvars => patch%aux%RT%auxvars
    if (Initialized(this%max_tracer_conc)) then
      tracer_scale = this%tracer_conductivity/this%max_tracer_conc
    endif
  endif

  material_auxvars => patch%aux%Material%auxvars
  do ghosted_id = 1, grid%ngmax
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
        cond_w = cond_w0 + tracer_scale * &
                           rt_auxvars(ghosted_id)%total(species_id,1)
      endif
    endif
    ! compute conductivity
    call ERTConductivityFromEmpiricalEqs(por,sat,a,m,n,Vc,cond_w,cond_c,cond)
    material_auxvars(ghosted_id)%electrical_conductivity(1) = cond
  enddo
 
end subroutine PMERTPreSolve

! ************************************************************************** !

subroutine PMERTSolve(this,time,ierr)

  ! Solves the linear systsem for ERT for all electrodes
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/27/21
  !
  use Patch_module
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
  PetscInt :: elec_id, local_elec_id
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: num_linear_iterations
  PetscReal :: val
  PetscReal :: average_cond
  PetscReal, pointer :: vec_ptr(:)

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

  if (OptionPrintToScreen(this%option)) then
    write(*,'(" Solving for electrode:")',advance='no')
  endif
  do ielec=1,nelec
    if (OptionPrintToScreen(this%option)) then
      write(*,'(x,a)',advance='no') trim(StringWrite(ielec))
    endif

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

    ! NB. solution is stored in field%work -> this can be an initial guess
    !call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)

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
    !call VecView(this%rhs,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)

    ! Solve system
    call PetscTime(log_ksp_start_time,ierr); CHKERRQ(ierr)
    call KSPSolve(solver%ksp,this%rhs,field%work,ierr);CHKERRQ(ierr)
    call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
    this%ksp_time = this%ksp_time + (log_end_time - log_ksp_start_time)
    !call VecView(field%work,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)

    call DiscretizationGlobalToLocal(discretization,field%work, &
                                     field%work_loc,ONEDOF)
    call VecGetArrayF90(field%work_loc,vec_ptr,ierr);CHKERRQ(ierr)
    ! store potentials for each electrode
    do ghosted_id=1,grid%ngmax
      if (patch%imat(ghosted_id) <= 0) cycle
      ert_auxvars(ghosted_id)%potential(ielec) = vec_ptr(ghosted_id)
    enddo
    call VecRestoreArrayF90(field%work_loc,vec_ptr,ierr);CHKERRQ(ierr)

    call KSPGetIterationNumber(solver%ksp,num_linear_iterations,ierr)
    call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr)
    this%linear_iterations_in_step = this%linear_iterations_in_step + &
                                     num_linear_iterations
  enddo

  ! Assemble solutions
  call PMERTAssembleSimulatedData(this,time)

  ! Build Jacobian
  if (this%option%geophysics%compute_jacobian) call PMERTBuildJacobian(this)

  ! For inversion
  if (this%option%geophysics%inversion) &
                           call PMERTUpdateElectricalConductivity(this)

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
  PetscInt :: ielec
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
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

  enddo

  ! write simuated data in a E4D .srv file
  if (option%myrank == option%io_rank) then
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
  use Material_Aux_class

  implicit none 

  class(pm_ert_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(survey_type), pointer :: survey
  type(ert_auxvar_type), pointer :: ert_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt, pointer :: cell_neighbors(:,:)
  PetscReal, allocatable :: phi_sor(:), phi_rec(:)
  PetscReal :: jacob
  PetscReal :: cond
  PetscInt :: idata
  PetscInt :: ielec
  PetscInt :: ia,ib,im,in
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
      ert_auxvars(ghosted_id)%jacobian(idata) = jacob * cond

      deallocate(phi_sor, phi_rec)
    enddo    
  enddo  

  ! I can now deallocate potential and delM (and M just after solving)
  ! But what about potential field output?

end subroutine PMERTBuildJacobian

! ************************************************************************** !

subroutine PMERTUpdateElectricalConductivity(this)
  !
  ! Computes conducivity update del_cond and use it update
  ! previous conductivity as cond_new = cond_prev + del_cond
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/11/21
  !

  use Patch_module
  use Grid_module

  implicit none

  class(pm_ert_type) :: this
  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(survey_type), pointer :: survey  
  type(inversion_type), pointer :: inversion
  
  patch => this%realization%patch
  grid => patch%grid

  survey => this%survey
  inversion => this%inversion
  
  ! Build Wm matrix
  call PMERTBuildWm(this)

  call InversionAllocateWorkArrays(inversion,survey,grid)

  ! get inversion%del_cond
  call PMERTCGLSSolve(this)

  ! TODO: Update material_auxvars(:)%electrical_conductivity(1)

  call InversionDeallocateWorkArrays(inversion)

end subroutine PMERTUpdateElectricalConductivity

! ************************************************************************** !

subroutine PMERTCGLSSolve(this)
  !
  ! Implements CGLS solver for least sqaure equivalent 
  !            of the normal equations  
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/07/21
  !
 
  implicit none

  class(pm_ert_type) :: this

  type(option_type), pointer :: option
  type(survey_type), pointer :: survey
  type(inversion_type), pointer :: inversion

  PetscInt :: i,nm,ncons
  PetscReal :: alpha,gbeta,gamma,gamma1,delta1,delta2,delta
  PetscReal :: norms0,norms,normx,xmax
  PetscReal :: resNE,resNE_old
  PetscBool :: exit_info,indefinite
  PetscErrorCode :: ierr

  PetscReal, parameter :: delta_initer = 1e-3
  PetscReal, parameter :: initer_conv  = 1e-4

  option => this%option
  survey => this%survey
  inversion => this%inversion

  inversion%del_cond = 0.0d0

  if (OptionPrintToScreen(this%option)) then
    write(*,'(" --> Solving normal equation using CGLS solver:")') 
  endif

  nm = survey%num_measurement
  ncons = inversion%num_constraints_local

  ! Get RHS vector inversion%b
  call PMERTCGLSRhs(this)

  inversion%r = inversion%b

  ! get inversion%s = J^tr
  call PMERTComputeMatVecProductJtr(this)
  inversion%p = inversion%s

  gamma = dot_product(inversion%s,inversion%s)
  call MPI_Allreduce(MPI_IN_PLACE,gamma,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

  norms0 = sqrt(gamma)
  xmax = 0.d0
  normx = 0.d0
  resNE = 0.d0
  exit_info = PETSC_FALSE
  indefinite = PETSC_FALSE

  do i=1,inversion%maxiter

    if (exit_info) exit

    ! get inversion%q = Jp
    call PMERTComputeMatVecProductJp(this)

    delta1 = dot_product(inversion%q(1:nm),inversion%q(1:nm))
    delta2 = dot_product(inversion%q(nm+1:nm+ncons),inversion%q(nm+1:nm+ncons))
    call MPI_Allreduce(MPI_IN_PLACE,delta2,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
    delta = delta1 + delta2

    if (delta < 0) indefinite = PETSC_TRUE
    if (delta == 0) delta = epsilon(delta)

    alpha = gamma / delta

    inversion%del_cond = inversion%del_cond + alpha * inversion%p
    inversion%r = inversion%r - alpha * inversion%q

    ! get inversion%s = J^tr
    call PMERTComputeMatVecProductJtr(this)
 
    gamma1 = gamma
    gamma = dot_product(inversion%s,inversion%s)
    call MPI_Allreduce(MPI_IN_PLACE,gamma,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

    norms = sqrt(gamma)
    gbeta = gamma / gamma1
    inversion%p = inversion%s + gbeta * inversion%p

    normx = dot_product(inversion%del_cond,inversion%del_cond)
    call MPI_Allreduce(MPI_IN_PLACE,normx,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
    normx = sqrt(normx)
    if (xmax < normx) xmax = normx
    if ( (norms <= norms0 * initer_conv) .or. (normx * initer_conv >= 1)) &
                               exit_info = PETSC_TRUE

    resNE_old = resNE
    resNE = norms / norms0

    if( abs((resNE_old - resNe) /resNE_old) < delta_initer .and. &
        i > inversion%miniter) exit_info = PETSC_TRUE
  enddo

end subroutine PMERTCGLSSolve

! ************************************************************************** !

subroutine PMERTCGLSRhs(this)
  !
  ! Builds RHS for least-square equation for CGLS solver  
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/11/21
  !
 
  use Patch_module
  use Material_Aux_class

  implicit none

  class(pm_ert_type) :: this

  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(survey_type), pointer :: survey
  type(inversion_type), pointer :: inversion
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: idata,iconst,irb,num_measurement
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: cond_ce,cond_nb,x     ! cell's and neighbor's
  PetscReal :: wm,beta

  option => this%option
  patch => this%realization%patch
  material_auxvars => patch%aux%Material%auxvars

  survey => this%survey
  inversion => this%inversion
  constrained_block => inversion%constrained_block
  rblock => inversion%rblock

  inversion%b = 0.0d0

  num_measurement = survey%num_measurement

  do idata=1,num_measurement
    inversion%b(idata) = survey%Wd(idata) * survey%Wd_cull(idata) * &
                         ( survey%dobs(idata) - survey%dsim(idata) )
  enddo

  beta = inversion%beta

  do iconst=1,inversion%num_constraints_local
    if (inversion%Wm(iconst) == 0) cycle

    wm = inversion%Wm(iconst)

    cond_ce = material_auxvars(rblock(iconst,1))%electrical_conductivity(1)
    irb = rblock(iconst,3)

    select case(constrained_block%structure_metric(irb))
    case(1)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = log(cond_ce) - log(cond_nb)
    case(2)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = log(cond_ce) - log(cond_nb)
    case(3)
      x = log(cond_ce) - log(constrained_block%reference_conductivity(irb))
    case(4)
      x = log(cond_ce) - log(constrained_block%reference_conductivity(irb))
    case(5)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = log(cond_ce) - log(cond_nb)
      ! TODO: compute rx,ry, and rz
    case(6)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = log(cond_ce) - log(cond_nb)
      ! TODO: compute rx,ry, and rz
    case(7)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = (log(cond_ce) - log(constrained_block%reference_conductivity(irb))) &
         -(log(cond_nb) - log(constrained_block%reference_conductivity(irb)))
    case(8)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = (log(cond_ce) - log(constrained_block%reference_conductivity(irb))) &
         -(log(cond_nb) - log(constrained_block%reference_conductivity(irb)))
    case(9)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = log(cond_ce) - log(cond_nb)
    case(10)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = log(cond_ce) - log(cond_nb)
    case default
      option%io_buffer = 'Supported STRUCTURE_METRIC in INVERSION, &
                          &CONSTRAINED_BLOCKS is between 1 to 10'
      call PrintErrMsg(option)
    end select

    inversion%b(num_measurement + iconst) = - sqrt(beta) * wm * x

  enddo

end subroutine PMERTCGLSRhs

! ************************************************************************** !

subroutine PMERTBuildWm(this)
  !
  ! Builds model regularization matrix: Wm
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/28/21
  !

  use Patch_module
  use Material_Aux_class
  use Inversion_module

  implicit none

  class(pm_ert_type) :: this

  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(inversion_type), pointer :: inversion
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: i

  option => this%option
  patch => this%realization%patch
  material_auxvars => patch%aux%Material%auxvars

  inversion => this%inversion
  constrained_block => inversion%constrained_block

  if (.not.associated(inversion%Wm)) call PMERTGetInfoAllocateWm(this)

  do i=1,inversion%num_constraints_local
    call ComputeWm(i,inversion%Wm(i))
  enddo

contains
  subroutine ComputeWm(iconst,wm)
    ! computes an element of Wm matrix
    !
    ! Author: Piyoosh Jaysaval
    ! Date: 05/28/21

    implicit none

    PetscInt :: iconst
    PetscReal :: wm

    PetscInt :: irb
    PetscReal :: x,awx,awy,awz
    PetscReal :: cond_ce,cond_nb     ! cell's and neighbor's
    PetscReal :: mn,sd
    PetscInt, pointer :: rblock(:,:)

    rblock => inversion%rblock

    ! get cond & block of the ith constrained eq.
    cond_ce = material_auxvars(rblock(iconst,1))%electrical_conductivity(1)
    irb = rblock(iconst,3)

    select case(constrained_block%structure_metric(irb))
    case(1)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = log(cond_ce) - log(cond_nb)
    case(2)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = abs(log(cond_ce) - log(cond_nb))
    case(3)
      x = log(cond_ce) - log(constrained_block%reference_conductivity(irb))
    case(4)
      x = abs(log(cond_ce) - &
              log(constrained_block%reference_conductivity(irb)))
    case(5)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = log(cond_ce) - log(cond_nb)
      ! TODO: compute rx,ry, and rz
    case(6)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = abs(log(cond_ce) - log(cond_nb))
      ! TODO: compute rx,ry, and rz
    case(7)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = (log(cond_ce) - log(constrained_block%reference_conductivity(irb))) &
         -(log(cond_nb) - log(constrained_block%reference_conductivity(irb)))
    case(8)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = abs( &
          (log(cond_ce) - log(constrained_block%reference_conductivity(irb))) &
         -(log(cond_nb) - log(constrained_block%reference_conductivity(irb))) )
    case(9)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = log(cond_ce) - log(cond_nb)
    case(10)
      cond_nb = material_auxvars(rblock(iconst,2))%electrical_conductivity(1)
      x = abs(log(cond_ce) - log(cond_nb))
    case default
      option%io_buffer = 'Supported STRUCTURE_METRIC in INVERSION, &
                          &CONSTRAINED_BLOCKS is between 1 to 10'
      call PrintErrMsg(option)
    end select

    mn = constrained_block%wf_mean(irb)
    sd = constrained_block%wf_sdev(irb)
    ! Get weight w
    select case(constrained_block%wf_type(irb))
    case(1)
      wm = 0.5 * (1 - erf( (x-mn)/sqrt(2*sd*sd) ))
    case(2)
      wm = 0.5 * (1 + erf( (x-mn)/sqrt(2*sd*sd) ))
    case(3)
      wm = 1 - exp(-((x-mn)*(x-mn)) / (2*sd*sd))
    case(4)
      wm = exp(-((x-mn)*(x-mn)) / (2*sd*sd))
    case(5)
      if((x-mn) < 0) then
         wm = 1 / (sd*sd)
      else
         wm = sd*sd / (((x-mn)*(x-mn) + sd*sd)*((x-mn)*(x-mn) + sd*sd))
      end if
    case(6)
      if((x-mn) > 0) then
         wm = 1 / (sd*sd)
      else
        wm = sd*sd / (((x-mn)*(x-mn) + sd*sd)*((x-mn)*(x-mn) + sd*sd))
      end if
    case default
      option%io_buffer = 'Supported WEIGHING_FUNCTION in INVERSION, &
                          &CONSTRAINED_BLOCKS is between 1 to 6'
      call PrintErrMsg(option)
    end select

    if(constrained_block%structure_metric(irb) == 5 .or. &
       constrained_block%structure_metric(irb) == 6) then
      awx = constrained_block%aniso_weight(irb,1)
      awy = constrained_block%aniso_weight(irb,2)
      awz = constrained_block%aniso_weight(irb,3)
      ! TODO: compute rx,ry, and rz before
      !wm = (1 - abs( awx*rx + awy*ry + awz*rz))**2
    end if

    wm = constrained_block%relative_weight(irb) * wm

  end subroutine ComputeWm

  end subroutine PMERTBuildWm

! ************************************************************************** !

subroutine PMERTGetInfoAllocateWm(this)
  !
  ! Allocate and get info on Wm and rblock
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/27/21
  !

  use Patch_module
  use Grid_module
  use Inversion_module

  implicit none

  class(pm_ert_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(inversion_type), pointer :: inversion
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: local_id,ghosted_id,ghosted_id_nbr
  PetscInt :: iconblock,inbr,ilink
  PetscInt :: num_constraints
  PetscInt :: num_neighbor

  patch => this%realization%patch
  grid => patch%grid

  inversion => this%inversion
  constrained_block => inversion%constrained_block

  num_constraints = 0
  do local_id=1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    do iconblock=1,constrained_block%num_constrained_block
      if (constrained_block%structure_metric(iconblock) > 0) then
        if (constrained_block%material_id(iconblock) == &
            patch%imat(ghosted_id)) then
          if (constrained_block%structure_metric(iconblock) == 3 .or. &
              constrained_block%structure_metric(iconblock) == 4) then
            num_constraints = num_constraints + 1
          else
            num_neighbor = grid%cell_neighbors_local_ghosted(0,local_id)
            do inbr=1,num_neighbor
              ghosted_id_nbr = abs( &
                            grid%cell_neighbors_local_ghosted(inbr,local_id))
              if (patch%imat(ghosted_id_nbr) /= patch%imat(ghosted_id)) then
                do ilink=1,constrained_block%block_link(iconblock,1)
                  if (constrained_block%block_link(iconblock,ilink+1) == &
                      patch%imat(ghosted_id_nbr)) then
                    num_constraints = num_constraints + 1
                  endif
                enddo
              else
                if (constrained_block%structure_metric(iconblock) < 9 .or. &
                    constrained_block%structure_metric(iconblock) > 10) then
                  num_constraints = num_constraints + 1
                endif
              endif
            enddo
          endif
        endif
      endif
    enddo
  enddo

  inversion%num_constraints_local = num_constraints
  allocate(inversion%Wm(num_constraints))
  allocate(inversion%rblock(num_constraints,THREE_INTEGER))
  inversion%Wm = 0.d0
  inversion%rblock = 0

  ! repeat once num_constraints is known
  num_constraints = 0
  do local_id=1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    do iconblock=1,constrained_block%num_constrained_block
      if (constrained_block%structure_metric(iconblock) > 0) then
        if (constrained_block%material_id(iconblock) == &
            patch%imat(ghosted_id)) then
          if (constrained_block%structure_metric(iconblock) == 3 .or. &
              constrained_block%structure_metric(iconblock) == 4) then
            num_constraints = num_constraints + 1
            inversion%rblock(num_constraints,1) = ghosted_id
            inversion%rblock(num_constraints,3) = iconblock
          else
            num_neighbor = grid%cell_neighbors_local_ghosted(0,local_id)
            do inbr=1,num_neighbor
              ghosted_id_nbr = abs( &
                            grid%cell_neighbors_local_ghosted(inbr,local_id))
              if (patch%imat(ghosted_id_nbr) /= patch%imat(ghosted_id)) then
                do ilink=1,constrained_block%block_link(iconblock,1)
                  if (constrained_block%block_link(iconblock,ilink+1) == &
                      patch%imat(ghosted_id_nbr)) then
                    num_constraints = num_constraints + 1
                    inversion%rblock(num_constraints,1) = ghosted_id
                    inversion%rblock(num_constraints,2) = ghosted_id_nbr
                    inversion%rblock(num_constraints,3) = iconblock
                  endif
                enddo
              else
                if (constrained_block%structure_metric(iconblock) < 9 .or. &
                    constrained_block%structure_metric(iconblock) > 10) then
                  num_constraints = num_constraints + 1
                  inversion%rblock(num_constraints,1) = ghosted_id
                  inversion%rblock(num_constraints,2) = ghosted_id_nbr
                  inversion%rblock(num_constraints,3) = iconblock
                endif
              endif
            enddo
          endif
        endif
      endif
    enddo
  enddo

end subroutine PMERTGetInfoAllocateWm

! ************************************************************************** !

subroutine PMERTComputeMatVecProductJp(this)
  !
  ! Computes product of Jacobian J with a vector p = Jp 
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/06/21
  !

  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module

  implicit none

  class(pm_ert_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(survey_type), pointer :: survey
  type(inversion_type), pointer :: inversion
  type(constrained_block_type), pointer :: constrained_block
  type(ert_auxvar_type), pointer :: ert_auxvars(:)

  PetscInt :: idata,iconst,irb,num_measurement
  PetscInt :: local_id,ghosted_id
  PetscInt :: local_id_nb,ghosted_id_nb
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: beta,wm
  PetscReal, pointer :: pvec_ptr(:)
  PetscErrorCode :: ierr

  option => this%option
  field => this%realization%field
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  ert_auxvars => patch%aux%ERT%auxvars

  survey => this%survey
  inversion => this%inversion
  constrained_block => inversion%constrained_block
  rblock => inversion%rblock

  inversion%q = 0.d0

  ! Data part
  call VecGetArrayF90(field%work,pvec_ptr,ierr);CHKERRQ(ierr)
  pvec_ptr = 0.d0

  do idata=1,survey%num_measurement
    do local_id=1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)         
      if (patch%imat(ghosted_id) <= 0) cycle
      inversion%q(idata) = inversion%q(idata) + &
                           ert_auxvars(ghosted_id)%jacobian(idata) * &
                           inversion%p(local_id)
      if (idata == 1) pvec_ptr(local_id) = inversion%p(local_id)
    enddo
    
    call MPI_Allreduce(MPI_IN_PLACE,inversion%q(idata),ONE_INTEGER_MPI, &  
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
  enddo

  call VecRestoreArrayF90(field%work,pvec_ptr,ierr);CHKERRQ(ierr)

  ! Model part
  ! Get local inversion%p to ghosted in pvec_ptr
  call DiscretizationGlobalToLocal(discretization,field%work, &
                                   field%work_loc,ONEDOF)
  call VecGetArrayF90(field%work_loc,pvec_ptr,ierr);CHKERRQ(ierr)

  num_measurement = survey%num_measurement
  beta = inversion%beta

  do iconst=1,inversion%num_constraints_local
    if (inversion%Wm(iconst) == 0) cycle

    wm = inversion%Wm(iconst)
    irb = rblock(iconst,3)
    ghosted_id = rblock(iconst,1)

    if (constrained_block%structure_metric(irb) == 3 .or. &
        constrained_block%structure_metric(irb) == 4) then
          inversion%q(num_measurement + iconst) = &
                sqrt(beta) * wm * pvec_ptr(ghosted_id)
    else
      ghosted_id_nb = rblock(iconst,2)
      inversion%q(num_measurement + iconst) = &
          sqrt(beta) * wm * (pvec_ptr(ghosted_id) - pvec_ptr(ghosted_id_nb))
    endif
  enddo

  call VecRestoreArrayF90(field%work_loc,pvec_ptr,ierr);CHKERRQ(ierr)

end subroutine PMERTComputeMatVecProductJp

! ************************************************************************** !

subroutine PMERTComputeMatVecProductJtr(this)
  !
  ! Computes product of Jacobian J transpose with a vector r = J^t x r  
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/06/21
  !

  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module

  implicit none

  class(pm_ert_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(survey_type), pointer :: survey
  type(inversion_type), pointer :: inversion
  type(constrained_block_type), pointer :: constrained_block
  type(ert_auxvar_type), pointer :: ert_auxvars(:)

  PetscInt :: idata,iconst,irb,num_measurement
  PetscInt :: local_id,ghosted_id
  PetscInt :: local_id_nb,ghosted_id_nb
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: beta,wm
  PetscReal, pointer :: svec_ptr(:)
  PetscErrorCode :: ierr

  field => this%realization%field
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  ert_auxvars => patch%aux%ERT%auxvars

  survey => this%survey
  inversion => this%inversion
  constrained_block => inversion%constrained_block
  rblock => inversion%rblock

  inversion%s = 0.0d0

  ! Model part
  call VecGetArrayF90(field%work_loc,svec_ptr,ierr);CHKERRQ(ierr)
  svec_ptr = 0.d0

  num_measurement = survey%num_measurement
  beta = inversion%beta

  do iconst=1,inversion%num_constraints_local
    if (inversion%Wm(iconst) == 0) cycle

    wm = inversion%Wm(iconst)
    irb = rblock(iconst,3)
    ghosted_id = rblock(iconst,1)

    if (constrained_block%structure_metric(irb) == 3 .or. &
        constrained_block%structure_metric(irb) == 4) then
          svec_ptr(ghosted_id) = svec_ptr(ghosted_id) + &
                sqrt(beta) * wm * inversion%r(num_measurement + iconst)
    else
      ghosted_id_nb = rblock(iconst,2)
      svec_ptr(ghosted_id) = svec_ptr(ghosted_id) + &
              sqrt(beta) * wm * inversion%r(num_measurement + iconst)
      svec_ptr(ghosted_id_nb) = svec_ptr(ghosted_id_nb) - &
              sqrt(beta) * wm * inversion%r(num_measurement + iconst)
    endif
  enddo

  call VecRestoreArrayF90(field%work_loc,svec_ptr,ierr);CHKERRQ(ierr)

  ! Data part
  call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
  call DiscretizationLocalToGlobalAdd(discretization,field%work_loc, &
                                   field%work,ONEDOF)

  call VecGetArrayF90(field%work,svec_ptr,ierr);CHKERRQ(ierr)

  do local_id=1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    do idata=1,survey%num_measurement
      svec_ptr(local_id) = svec_ptr(local_id) + &
                             ert_auxvars(ghosted_id)%jacobian(idata) * &
                             inversion%r(idata)
    enddo
    inversion%s(local_id) = svec_ptr(local_id)
  enddo

  call VecRestoreArrayF90(field%work,svec_ptr,ierr);CHKERRQ(ierr)

end subroutine PMERTComputeMatVecProductJtr

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

  character(len=MAXWORDLENGTH) :: word
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
  nullify(this%inversion)
  call WaypointListDestroy(this%waypoint_list)

  call DeallocateArray(this%species_conductivity_coef)
  if (this%rhs /= PETSC_NULL_VEC) then
    call VecDestroy(this%rhs,ierr);CHKERRQ(ierr)
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
