module PM_ERT_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_module
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
    type(waypoint_list_type), pointer :: waypoint_list
    PetscInt :: linear_iterations_in_step
    PetscLogDouble :: ksp_time
    Vec :: rhs
  contains
    procedure, public :: Setup => PMERTSetup
    procedure, public :: ReadSimulationOptionsBlock => PMERTReadSimOptionsBlock
    procedure, public :: SetRealization => PMERTSetRealization
    procedure, public :: InitializeRun => PMERTInitializeRun
    procedure, public :: SetupSolvers => PMERTSetupSolvers
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
  pm_ert%waypoint_list => WaypointListCreate()

  pm_ert%linear_iterations_in_step = 0
  pm_ert%ksp_time = 0.d0
  pm_ert%rhs = PETSC_NULL_VEC

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

  option => this%option

  error_string = 'ERT Options'

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
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

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
  use Discretization_module

  implicit none

  class(pm_ert_type) :: this
  PetscErrorCode :: ierr

  PetscReal, pointer :: vec_ptr(:)

  call DiscretizationDuplicateVector(this%realization%discretization, &
                                     this%realization%field%work,this%rhs)

  ! Initialize to zeros
  call VecZeroEntries(this%rhs,ierr);CHKERRQ(ierr)

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
  use ERT_module
  use String_module
  use Survey_module

  implicit none

  class(pm_ert_type) :: this
  PetscReal :: time
  PetscErrorCode :: solve_ierr

  class(realization_subsurface_type), pointer :: realization
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(solver_type), pointer :: solver
  type(field_type), pointer :: field
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

  call PMBasePrintHeader(this)
  solve_ierr = 0
  ! Forward solve start
  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

  solver => this%solver
  survey => this%survey
  realization => this%realization
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  ert_auxvars => patch%aux%ERT%auxvars

  ! Build System matrix
  call ERTCalculateMatrix(realization,solver%M)
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

    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

    ! store potentials for each electrode
    do local_id=1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      ert_auxvars(ghosted_id)%potential(ielec) = vec_ptr(local_id)
    enddo

    call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call KSPGetIterationNumber(solver%ksp,num_linear_iterations,ierr)
    call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr)
    this%linear_iterations_in_step = this%linear_iterations_in_step + &
                                     num_linear_iterations
  enddo

  ! Assemble solutions
  call PMERTAssembleSimulatedData(this,time)
  !print*,survey%dsim

end subroutine PMERTSolve

! ************************************************************************** !

subroutine PMERTAssembleSimulatedData(this,time)

  use Patch_module
  use Grid_module
  use Survey_module

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

  implicit none

  class(pm_ert_type) :: this

  PetscErrorCode :: ierr

  call PMBaseDestroy(this)
  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%survey)
  call WaypointListDestroy(this%waypoint_list)

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
