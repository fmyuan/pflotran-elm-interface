module Inversion_Subsurface_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Aux_module
  use Inversion_TS_Aux_module
  use Inversion_Base_class
  use Realization_Subsurface_class
  use Simulation_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_base_type) :: inversion_subsurface_type
    character(len=MAXSTRINGLENGTH) :: forward_simulation_filename
    class(simulation_subsurface_type), pointer :: forward_simulation
    class(realization_subsurface_type), pointer :: realization
    type(inversion_aux_type), pointer :: inversion_aux
    PetscReal, pointer :: measurement(:)
    PetscInt, pointer :: imeasurement(:)
    PetscInt :: measurement_offset
    PetscInt :: iqoi(2)
    PetscInt :: iobsfunc
    PetscInt :: n_qoi_per_cell
    Vec :: quantity_of_interest
    Vec :: ref_quantity_of_interest
    Vec :: measurement_vec
    character(len=MAXWORDLENGTH) :: ref_qoi_dataset_name
    PetscBool :: print_sensitivity_jacobian
    PetscBool :: debug_adjoint
    PetscInt :: debug_verbosity
    PetscBool :: local_adjoint
    VecScatter :: scatter_global_to_measurement
  contains
    procedure, public :: Init => InversionSubsurfaceInit
    procedure, public :: ReadBlock => InversionSubsurfReadBlock
    procedure, public :: Initialize => InversionSubsurfInitialize
    procedure, public :: Step => InversionSubsurfaceStep
    procedure, public :: ConnectToForwardRun => InvSubsurfConnectToForwardRun
    procedure, public :: CalculateSensitivity => InvSubsurfCalculateSensitivity
    procedure, public :: OutputSensitivity => InvSubsurfOutputSensitivity
    procedure, public :: Invert => InversionSubsurfaceInvert
    procedure, public :: Strip => InversionSubsurfaceStrip
  end type inversion_subsurface_type

  public :: InversionSubsurfaceCreate, &
            InversionSubsurfaceInit, &
            InversionSubsurfReadSelectCase, &
            InversionSubsurfInitialize, &
            InvSubsurfConnectToForwardRun, &
            InvSubsurfOutputSensitivity, &
            InversionSubsurfaceStrip

contains

! ************************************************************************** !

function InversionSubsurfaceCreate(driver)
  !
  ! Allocates and initializes a new subsurface inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/21
  !
  use Driver_module

  class(driver_type), pointer :: driver

  class(inversion_subsurface_type), pointer :: InversionSubsurfaceCreate

  allocate(InversionSubsurfaceCreate)
  call InversionSubsurfaceCreate%Init(driver)

end function InversionSubsurfaceCreate

! ************************************************************************** !

subroutine InversionSubsurfaceInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  use Driver_module
  use ZFlow_Aux_module, only : zflow_calc_adjoint

  class(inversion_subsurface_type) :: this
  class(driver_type), pointer :: driver

  call InversionBaseInit(this,driver)

  this%quantity_of_interest = PETSC_NULL_VEC
  this%iqoi = UNINITIALIZED_INTEGER
  this%iobsfunc = UNINITIALIZED_INTEGER
  this%n_qoi_per_cell = UNINITIALIZED_INTEGER
  this%measurement_vec = PETSC_NULL_VEC
  this%ref_quantity_of_interest = PETSC_NULL_VEC
  this%ref_qoi_dataset_name = ''
  this%forward_simulation_filename = ''
  this%print_sensitivity_jacobian = PETSC_FALSE
  this%debug_adjoint = PETSC_FALSE
  this%local_adjoint = PETSC_FALSE
  this%debug_verbosity = UNINITIALIZED_INTEGER
  this%scatter_global_to_measurement = PETSC_NULL_VECSCATTER
  this%measurement_offset = UNINITIALIZED_INTEGER

  nullify(this%measurement)
  nullify(this%imeasurement)

  nullify(this%forward_simulation)
  nullify(this%realization)
  nullify(this%inversion_aux)

  zflow_calc_adjoint = PETSC_TRUE

end subroutine InversionSubsurfaceInit

! ************************************************************************** !

subroutine InversionSubsurfReadBlock(this,input,option)

  use Input_Aux_module
  use Option_module
  use String_module

  class(inversion_subsurface_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'Test Inversion'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_TRUE
    call InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine InversionSubsurfReadBlock

! ************************************************************************** !

subroutine InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                          error_string,option)

  use Input_Aux_module
  use Option_module
  use String_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, POROSITY, &
                               LIQUID_PRESSURE, LIQUID_SATURATION
  use Material_Aux_class, only : POROSITY_BASE
  use Utility_module

  class(inversion_subsurface_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, pointer :: tempint(:)
  PetscReal, pointer :: tempreal(:)

  found = PETSC_TRUE
  call InversionBaseReadSelectCase(this,input,keyword,found, &
                                   error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('FORWARD_SIMULATION_FILENAME')
      call InputReadFilename(input,option,this%forward_simulation_filename)
      call InputErrorMsg(input,option,keyword,error_string)
    case('QUANTITY_OF_INTEREST')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,keyword,error_string)
      call StringToUpper(word)
      select case(word)
        case('ELECTRICAL_CONDUCTIVITY')
          this%iqoi(1) = ELECTRICAL_CONDUCTIVITY
        case('PERMEABILITY')
          this%iqoi(1) = PERMEABILITY
        case('POROSITY')
          this%iqoi(1) = POROSITY
          this%iqoi(2) = POROSITY_BASE
        case default
          call InputKeywordUnrecognized(input,word,trim(error_string)// &
                                        & ','//trim(keyword),option)
      end select
    case('REFERENCE_QUANTITY_OF_INTEREST')
      call InputReadNChars(input,option,this%ref_qoi_dataset_name, &
                           MAXWORDLENGTH,PETSC_TRUE)
      call InputErrorMsg(input,option,'DATASET NAME', &
                         keyword)
    case('OBSERVATION_FUNCTION')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,keyword,error_string)
      call StringToUpper(word)
      select case(word)
        case('LIQUID_PRESSURE')
          this%iobsfunc = LIQUID_PRESSURE
        case('LIQUID_SATURATION')
          this%iobsfunc = LIQUID_SATURATION
        case default
          call InputKeywordUnrecognized(input,word,trim(error_string)// &
                                        & ','//trim(keyword),option)
      end select
    case('MEASUREMENTS')
      string = trim(error_string)//keyword
      i = 10
      allocate(tempint(i))
      tempint = UNINITIALIZED_INTEGER
      allocate(tempreal(i))
      tempreal = UNINITIALIZED_DOUBLE
      i = 0
      do
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,error_string)
        if (InputCheckExit(input,option)) exit
        i = i + 1
        if (i > size(tempint)) then
          call ReallocateArray(tempint)
          call ReallocateArray(tempreal)
        endif
        call InputReadInt(input,option,tempint(i))
        call InputErrorMsg(input,option,'cell id',string)
        call InputReadDouble(input,option,tempreal(i))
        call InputErrorMsg(input,option,'measurement',string)
      enddo
      allocate(this%imeasurement(i))
      this%imeasurement(:) = tempint(1:i)
      allocate(this%measurement(i))
      this%measurement(:) = tempreal(1:i)
      call DeallocateArray(tempint)
      call DeallocateArray(tempreal)
    case('PRINT_SENSITIVITY_JACOBIAN')
      this%print_sensitivity_jacobian = PETSC_TRUE
    case('DEBUG_ADJOINT')
      this%debug_adjoint = PETSC_TRUE
      call InputReadInt(input,option,i)
      if (input%ierr == 0) then
        this%debug_verbosity = i
      endif
    case('LOCAL_ADJOINT')
      this%local_adjoint = PETSC_TRUE
    case default
      found = PETSC_FALSE
  end select

end subroutine InversionSubsurfReadSelectCase

! ************************************************************************** !

subroutine InversionSubsurfInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Connection_module
  use Discretization_module
  use Coupler_module
  use Grid_module
  use Patch_module

  class(inversion_subsurface_type) :: this

  type(patch_type), pointer :: patch
  type(inversion_ts_aux_type), pointer :: inversion_ts_aux
  PetscInt :: i
  PetscInt :: sum_connection
  PetscInt :: num_measurements, num_measurements_local
  PetscInt :: num_parameters_local, num_parameters_global
  PetscInt, allocatable :: int_array(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal :: tempreal
  Vec :: v
  IS :: is_petsc
  PetscErrorCode :: ierr

  nullify(vec_ptr)

  patch => this%realization%patch
  if (.not.associated(this%inversion_aux)) then
    this%n_qoi_per_cell = 1 ! 1 perm per cell

    this%inversion_aux => InversionAuxCreate()
    num_measurements = size(this%imeasurement)
    num_parameters_local = patch%grid%nlmax*this%n_qoi_per_cell
    num_parameters_global = patch%grid%nmax*this%n_qoi_per_cell
    ! JsensitivityT is the transpose of the sensitivity Jacobian
    ! with num measurement columns and num parameter rows
    call MatCreateDense(this%driver%comm%mycomm, &
                        num_parameters_local,PETSC_DECIDE, &
                        num_parameters_global,num_measurements, &
                        PETSC_NULL_SCALAR, &
                        this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
    call MatZeroEntries(this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
    ! cannot pass in this%measurement_vec as it is initialized to
    ! PETSC_NULL_VEC and MatCreateVecs keys off that input
    call MatCreateVecs(this%inversion_aux%JsensitivityT,v,PETSC_NULL_VEC, &
                       ierr);CHKERRQ(ierr)
    this%measurement_vec = v
    call MatGetLocalSize(this%inversion_aux%JsensitivityT, &
                         PETSC_NULL_INTEGER, &
                         num_measurements_local,ierr);CHKERRQ(ierr)
    ! must initialize to zero...must be a bug in MPICh
    this%measurement_offset = 0
    call MPI_Exscan(num_measurements_local,this%measurement_offset, &
                    ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                    this%driver%comm%mycomm,ierr);CHKERRQ(ierr)
    !TODO(geh): compress the mappings
    !call GridMapCellsToConnections(patch%grid, &
    !                           this%inversion_aux%cell_to_internal_connection)

#if 0
    if (this%local_adjoint) then
      sum_connection = &
        ConnectionGetNumberInList(patch%grid%internal_connection_set_list)
      sum_connection2 = &
        CouplerGetNumConnectionsInList(patch%boundary_condition_list)
      ! set up pointer to solution vec
      inversion_ts_aux%mat_vec_solution_ptr%solution = &
        this%realization%field%flow_xx
! the old flux approach
      call InvTSAuxAllocateFluxCoefArrays(inversion_ts_aux,sum_connection, &
                                          sum_connection2)
      allocate(int_array(patch%grid%nlmax))
      int_array = 0
      boundary_condition => patch%boundary_condition_list%first
      sum_connection = 0
      do
        if (.not.associated(boundary_condition)) exit
        cur_connection_set => boundary_condition%connection_set
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1
          local_id = cur_connection_set%id_dn(iconn)
          int_array(local_id) = int_array(local_id) + 1
        enddo
        boundary_condition => boundary_condition%next
      enddo
      iconn = maxval(int_array)
      deallocate(int_array)
      allocate(this%inversion_aux% &
                cell_to_bc_connection(0:iconn,patch%grid%nlmax))
      this%inversion_aux%cell_to_bc_connection = 0
      boundary_condition => patch%boundary_condition_list%first
      sum_connection = 0
      do
        if (.not.associated(boundary_condition)) exit
        cur_connection_set => boundary_condition%connection_set
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1
          local_id = cur_connection_set%id_dn(iconn)
          this%inversion_aux%cell_to_bc_connection(0,local_id) = &
            this%inversion_aux%cell_to_bc_connection(0,local_id) + 1
            this%inversion_aux%cell_to_bc_connection( &
              this%inversion_aux%cell_to_bc_connection(0,local_id), &
              local_id) = sum_connection
        enddo
        boundary_condition => boundary_condition%next
      enddo
    endif
#endif

    ! map measurement vec to the solution vector
    if (this%driver%comm%myrank == 0) then
      do i = 1, num_measurements
        tempreal = dble(this%imeasurement(i))
        call VecSetValue(this%measurement_vec,i-1,tempreal, &
                         INSERT_VALUES,ierr);CHKERRQ(ierr)
      enddo
    endif
    call VecAssemblyBegin(this%measurement_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(this%measurement_vec,ierr);CHKERRQ(ierr)
    allocate(int_array(num_measurements_local))
    call VecGetArrayF90(this%measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)
    do i = 1, num_measurements_local
      int_array(i) = int(vec_ptr(i)+1.d-5)
    enddo
    int_array = int_array - 1
    call VecRestoreArrayF90(this%measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)
    call DiscretAOApplicationToPetsc(this%realization%discretization, &
                                     int_array)
    call ISCreateGeneral(this%driver%comm%mycomm,num_measurements_local, &
                         int_array, &
                         PETSC_COPY_VALUES,is_petsc,ierr);CHKERRQ(ierr)
    call VecScatterCreate(this%realization%field%work,is_petsc, &
                          this%measurement_vec,PETSC_NULL_IS, &
                          this%scatter_global_to_measurement, &
                          ierr);CHKERRQ(ierr)
    call ISDestroy(is_petsc,ierr);
    deallocate(int_array)
  endif

  ! create inversion_ts_aux for first time step
  nullify(inversion_ts_aux) ! must pass in null object
  inversion_ts_aux => InversionTSAuxCreate(inversion_ts_aux)
  ! set up pointer to M matrix
  inversion_ts_aux%mat_vec_solution_ptr%M = &
    this%forward_simulation%flow_process_model_coupler%timestepper%solver%M

  this%inversion_aux%inversion_ts_aux_list => inversion_ts_aux
  call InvTSAuxAllocate(inversion_ts_aux,patch%grid%nlmax)


end subroutine InversionSubsurfInitialize

! ************************************************************************** !

subroutine InversionSubsurfaceStep(this)
  !
  ! Execute a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Option_module
  use Factory_Forward_module

  class(inversion_subsurface_type) :: this

  type(option_type), pointer :: option

  option => OptionCreate()
  write(option%group_prefix,'(i6)') this%iteration
  option%group_prefix = 'Run' // trim(adjustl(option%group_prefix))
  call OptionSetDriver(option,this%driver)
  call OptionSetInversionOption(option,this%inversion_option)
  call FactoryForwardInitialize(this%forward_simulation, &
                                this%forward_simulation_filename,option)
  this%realization => this%forward_simulation%realization
  call this%Initialize()
  call this%forward_simulation%InitializeRun()
  call this%ConnectToForwardRun()
  if (option%status == PROCEED) then
    call this%forward_simulation%ExecuteRun()
  endif
  call this%CalculateSensitivity()
  call this%OutputSensitivity('')
  nullify(this%realization)
  call this%forward_simulation%FinalizeRun()
  call this%forward_simulation%Strip()
  deallocate(this%forward_simulation)
  nullify(this%forward_simulation)
  call this%Invert()

  this%converg_flag = PETSC_FALSE
  if (this%iteration > this%maximum_iteration) this%converg_flag = PETSC_TRUE

end subroutine InversionSubsurfaceStep

! ************************************************************************** !

subroutine InvSubsurfConnectToForwardRun(this)
  !
  ! Sets up the interface between a single forward run and the outer
  ! inversion wrapper
  !
  ! Author: Glenn Hammond
  ! Date: 10/18/21
  !
  use Discretization_module
  use Material_module

  class(inversion_subsurface_type) :: this

  PetscErrorCode :: ierr

  this%realization%patch%aux%inversion_ts_aux => &
    this%inversion_aux%inversion_ts_aux_list

  ! on first pass, store and set thereafter
  if (this%quantity_of_interest == PETSC_NULL_VEC) then
    ! can't move VecDuplicate earlier as it will not be null for the
    ! conditional above
    call VecDuplicate(this%realization%field%work, &
                      this%quantity_of_interest,ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi(1),this%iqoi(2))
    call DiscretizationLocalToGlobal(this%realization%discretization, &
                                     this%realization%field%work_loc, &
                                     this%quantity_of_interest,ONEDOF)
  endif
  call DiscretizationGlobalToLocal(this%realization%discretization, &
                                   this%quantity_of_interest, &
                                   this%realization%field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                               this%realization%field%work_loc, &
                               this%iqoi(1),this%iqoi(2))

end subroutine InvSubsurfConnectToForwardRun

! ************************************************************************** !

subroutine InvSubsurfCalculateSensitivity(this)
  !
  ! Calculates sensitivity matrix Jsensitivity
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Option_module
  use String_module
  use Timer_class
  use Utility_module

  class(inversion_subsurface_type) :: this

  type(inversion_aux_type), pointer :: inversion_aux
  type(inversion_ts_aux_type), pointer :: cur_inversion_ts_aux
  type(inversion_ts_aux_type), pointer :: prev_inversion_ts_aux
  type(option_type), pointer :: option
  class(timer_type), pointer :: timer
  PetscErrorCode :: ierr

  option => this%realization%option
  inversion_aux => this%inversion_aux

  timer => TimerCreate()

  call timer%Start()

  call PrintHeader('SENSITIVITY JACOBIAN',option)

  ! go to end of list
  cur_inversion_ts_aux => inversion_aux%inversion_ts_aux_list
  if (.not.associated(cur_inversion_ts_aux)) then
    option%io_buffer = 'Inversion timestep auxiliary list is NULL.'
    call PrintErrMsg(option)
  endif
  do
    if (.not.associated(cur_inversion_ts_aux%next)) exit
    cur_inversion_ts_aux => cur_inversion_ts_aux%next
  enddo

  ! the last link should be allocated, but not populated. this is by design
  if (cur_inversion_ts_aux%dResdu /= PETSC_NULL_MAT) then
    option%io_buffer = 'Last link in Inversion timestep auxiliary list &
      &is not NULL.'
    call PrintErrMsg(option)
  else
    ! remove the last link
    prev_inversion_ts_aux => cur_inversion_ts_aux%prev
    if (.not.associated(prev_inversion_ts_aux)) then
      option%io_buffer = 'Next to last link in Inversion timestep &
        &auxiliary list is NULL.'
      call PrintErrMsg(option)
    endif
    nullify(prev_inversion_ts_aux%next)
    call InversionTSAuxDestroy(cur_inversion_ts_aux)
    ! point cur_inversion_ts_aux to the end of the list
    cur_inversion_ts_aux => prev_inversion_ts_aux
    this%inversion_aux%max_ts = cur_inversion_ts_aux%timestep
  endif

  ! work backward through list
  call OptionPrint(' Working backward through inversion_ts_aux list &
                   &calculating lambdas.',option)
  do
    if (.not.associated(cur_inversion_ts_aux)) exit
    print *, '  call InvSubsurfCalcLambda: ', cur_inversion_ts_aux%timestep
    call InvSubsurfCalcLambda(this,cur_inversion_ts_aux)
    cur_inversion_ts_aux => cur_inversion_ts_aux%prev
  enddo

  ! work forward through list
  call OptionPrint(' Working forward through inversion_ts_aux list &
                   &calculating sensitivity coefficients.',option)
  call MatZeroEntries(inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
  cur_inversion_ts_aux => inversion_aux%inversion_ts_aux_list
  do
    if (.not.associated(cur_inversion_ts_aux)) exit
    print *, '  call InvSubsurfAddSensitivity: ', cur_inversion_ts_aux%timestep
    call InvSubsurfAddSensitivity(this,cur_inversion_ts_aux)
    cur_inversion_ts_aux => cur_inversion_ts_aux%next
  enddo
  call MatAssemblyBegin(inversion_aux%JsensitivityT, &
                        MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(inversion_aux%JsensitivityT, &
                      MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  call InversionTSAuxListDestroy(inversion_aux%inversion_ts_aux_list,PETSC_TRUE)

  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime())) // &
    ' seconds to build all sensitivities.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InvSubsurfCalculateSensitivity

! ************************************************************************** !

subroutine InvSubsurfCalcLambda(this,inversion_ts_aux)
  !
  ! Calculates sensitivity matrix Jsensitivity
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Connection_module
  use Debug_module
  use Discretization_module
  use Grid_module
  use Option_module
  use Patch_module
  use Realization_Base_class
  use Solver_module
  use String_module
  use Timer_class
  use ZFlow_Aux_module
  use Variables_module

  use PM_Base_class
  use PM_ZFlow_class

  class(inversion_subsurface_type) :: this

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(inversion_aux_type), pointer :: inversion_aux
  type(inversion_ts_aux_type), pointer :: inversion_ts_aux
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(solver_type), pointer :: solver
  type(zflow_auxvar_type), pointer :: zflow_auxvars(:,:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal :: tempreal
  PetscInt :: imeasurement
  PetscInt :: icell_measurement
  Vec :: work
  Vec :: p
  Vec :: rhs
  ! derivative of residual at k+1 time level wrt unknown at k time level
  ! times lambda at k time level
  Vec :: dReskp1_duk_lambdak
  Vec :: natural_vec
  class(timer_type), pointer :: timer
  PetscErrorCode :: ierr

  nullify(vec_ptr)

  solver => this%forward_simulation%flow_process_model_coupler% &
              timestepper%solver
  option => this%realization%option
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  inversion_aux => this%inversion_aux
  zflow_auxvars => patch%aux%ZFlow%auxvars

  timer => TimerCreate()

  call timer%Start()

  work = this%realization%field%work ! DO NOT DESTROY!
  call VecDuplicate(work,p,ierr);CHKERRQ(ierr)
  call VecDuplicate(work,rhs,ierr);CHKERRQ(ierr)
  call VecDuplicate(work,dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)

  if (this%debug_verbosity > 2) then
    if (OptionPrintToScreen(option)) print *, 'M'
    call MatView(inversion_ts_aux%dResdu,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
  endif
  call VecDuplicateVecsF90(work,size(this%imeasurement), &
                           inversion_ts_aux%lambda,ierr);CHKERRQ(ierr)
  call KSPSetOperators(solver%ksp,inversion_ts_aux%dResdu, &
                       inversion_ts_aux%dResdu,ierr);CHKERRQ(ierr)
  do imeasurement = 1, size(this%imeasurement)
    if (inversion_aux%max_ts == inversion_ts_aux%timestep) then
      call VecZeroEntries(natural_vec,ierr);CHKERRQ(ierr)
      if (option%myrank == 0) then
        icell_measurement = this%imeasurement(imeasurement)
        tempreal = -1.d0
        call VecSetValue(natural_vec,icell_measurement-1,tempreal, &
                        INSERT_VALUES,ierr);CHKERRQ(ierr)
      endif
      call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
      call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
      call DiscretizationNaturalToGlobal(discretization,natural_vec,p,ONEDOF)
      select case(this%iobsfunc)
        case(LIQUID_PRESSURE)
        case(LIQUID_SATURATION)
          call RealizationGetVariable(this%realization,work,DERIVATIVE, &
                                      ZFLOW_LIQ_SAT_WRT_LIQ_PRES)
          call VecPointwiseMult(p,p,work,ierr);CHKERRQ(ierr)
      end select
      if (this%debug_verbosity > 2) then
        if (OptionPrintToScreen(option)) print *, 'p'
        call VecView(p,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      endif
      call VecZeroEntries(dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
    else
      call VecZeroEntries(p,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(dReskp1_duk_lambdak,vec_ptr,ierr);CHKERRQ(ierr)
      vec_ptr(:) = inversion_ts_aux%next%dRes_du_k(:)
      call VecRestoreArrayF90(dReskp1_duk_lambdak,vec_ptr,ierr);CHKERRQ(ierr)
      call VecPointwiseMult(dReskp1_duk_lambdak,dReskp1_duk_lambdak, &
                            inversion_ts_aux%next%lambda(imeasurement), &
                            ierr);CHKERRQ(ierr)
    endif
    call VecWAXPY(rhs,-1.d0,dReskp1_duk_lambdak,p,ierr);CHKERRQ(ierr)
    call KSPSolveTranspose(solver%ksp,rhs, &
                           inversion_ts_aux%lambda(imeasurement), &
                           ierr);CHKERRQ(ierr)
    if (this%debug_verbosity > 2) then
      if (OptionPrintToScreen(option)) print *, 'lambda'
      call VecView(inversion_ts_aux%lambda(imeasurement), &
                   PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      if (option%comm%mycommsize == 1) then
        call VecGetArrayF90(inversion_ts_aux%lambda(imeasurement), &
                            vec_ptr,ierr);CHKERRQ(ierr)
        print *, vec_ptr(:)
        call VecRestoreArrayF90(inversion_ts_aux%lambda(imeasurement), &
                                vec_ptr,ierr);CHKERRQ(ierr)
      endif
    endif
  enddo
  call VecDestroy(p,ierr);CHKERRQ(ierr)
  call VecDestroy(rhs,ierr);CHKERRQ(ierr)
  call VecDestroy(dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)

  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime()))
  option%io_buffer = trim(option%io_buffer) // &
    ' seconds to calculate lambdas.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InvSubsurfCalcLambda

! ************************************************************************** !

subroutine InvSubsurfAddSensitivity(this,inversion_ts_aux)
  !
  ! Calculates sensitivity matrix Jsensitivity
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Connection_module
  use Debug_module
  use Discretization_module
  use Grid_module
  use Option_module
  use Patch_module
  use Solver_module
  use String_module
  use Timer_class

  use PM_Base_class
  use PM_ZFlow_class

  class(inversion_subsurface_type) :: this

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(inversion_aux_type), pointer :: inversion_aux
  type(inversion_ts_aux_type), pointer :: inversion_ts_aux
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(solver_type), pointer :: solver
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscReal :: tempreal
  PetscInt :: iparameter, imeasurement
  PetscInt :: natural_id
  Vec :: work
  Vec :: work_loc, lambda_loc
  Vec :: dResdKLambda
  PetscViewer :: viewer
  class(timer_type), pointer :: timer
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  nullify(vec_ptr)
  nullify(vec_ptr2)

  lambda_loc = PETSC_NULL_VEC
  dResdKLambda = PETSC_NULL_VEC

  solver => this%forward_simulation%flow_process_model_coupler% &
              timestepper%solver
  option => this%realization%option
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  inversion_aux => this%inversion_aux

  timer => TimerCreate()

  call timer%Start()

  work = this%realization%field%work ! DO NOT DESTROY!

  if (this%debug_adjoint) then
    string = 'dResdK_ts'//trim(StringWrite(inversion_ts_aux%timestep))//'.txt'
    call PetscViewerASCIIOpen(option%mycomm,string, &
                              viewer,ierr);CHKERRQ(ierr)
    call MatView(inversion_ts_aux%dResdK,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  do imeasurement = 1, size(this%imeasurement)
    if (this%debug_adjoint) then
      string = 'lambda_ts'//trim(StringWrite(inversion_ts_aux%timestep)) // &
                '_' // trim(StringWrite(this%imeasurement(imeasurement))) // &
                '.txt'
      call PetscViewerASCIIOpen(option%mycomm,string, &
                                viewer,ierr);CHKERRQ(ierr)
      call VecView(inversion_ts_aux%lambda(imeasurement),viewer, &
                    ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif
    call VecDuplicate(work,dResdKLambda,ierr);CHKERRQ(ierr)
    call MatMultTranspose(inversion_ts_aux%dResdK, &
                          inversion_ts_aux%lambda(imeasurement), &
                          dResdKLambda,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(dResdKLambda,vec_ptr,ierr);CHKERRQ(ierr)
    do iparameter = 1, grid%nlmax
      natural_id = grid%nG2A(grid%nL2G(iparameter))
      call MatSetValue(inversion_aux%JsensitivityT,natural_id-1,imeasurement-1, &
                      vec_ptr(iparameter),ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
    call VecRestoreArrayF90(dResdKLambda,vec_ptr,ierr);CHKERRQ(ierr)
    call VecDestroy(dResdKLambda,ierr);CHKERRQ(ierr)
  enddo

  if (lambda_loc /= PETSC_NULL_VEC) then
    call VecDestroy(lambda_loc,ierr);CHKERRQ(ierr)
  endif
  if (dResdKLambda /= PETSC_NULL_VEC) then
    call VecDestroy(dResdKLambda,ierr);CHKERRQ(ierr)
  endif

  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime()))
  option%io_buffer = trim(option%io_buffer) // &
    ' seconds to add contributions to Jsensitivity.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InvSubsurfAddSensitivity

! ************************************************************************** !

subroutine InvSubsrfBMinusSM(this,local_id,lambda_ptr,solution,value_)
  !
  !
  !
  ! Author: Glenn Hammond
  ! Date: 11/15/21
  !
  use Connection_module
  use Grid_module
  use Option_module

  class(inversion_subsurface_type) :: this
  PetscInt :: local_id
  PetscReal, pointer :: lambda_ptr(:)
  PetscReal, pointer :: solution(:)
  PetscReal :: value_

  type(connection_set_type), pointer :: connection_set
  type(grid_type), pointer :: grid
  type(inversion_aux_type), pointer :: inversion_aux
  type(inversion_ts_aux_type), pointer :: inversion_ts_aux
  PetscReal :: Mlambda_up
  PetscReal :: Mlambda_dn
  PetscReal :: rhs_up
  PetscReal :: rhs_dn
  PetscReal :: Mlambda(this%realization%patch%grid%ngmax)
  PetscReal :: rhs(this%realization%patch%grid%ngmax)
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: ghosted_id
  PetscInt :: iconn
  PetscInt :: i

  grid => this%realization%patch%grid
  connection_set => grid%internal_connection_set_list%first
  inversion_aux => this%inversion_aux
  inversion_ts_aux => inversion_aux%inversion_ts_aux_list

  Mlambda = 0.d0
  rhs = 0.d0
  ghosted_id = grid%nL2G(local_id)
  do i = 1, inversion_aux%cell_to_internal_connection(0,local_id)
    iconn = inversion_aux%cell_to_internal_connection(i,local_id)
    ghosted_id_up = connection_set%id_up(iconn)
    ghosted_id_dn = connection_set%id_dn(iconn)
    Mlambda_up = 0.d0
    Mlambda_dn = 0.d0
    rhs_up = 0.d0
    rhs_dn = 0.d0
    if (ghosted_id == ghosted_id_up) then
      Mlambda_up = Mlambda_up + &
        inversion_ts_aux%dFluxdIntConn(1,iconn)*lambda_ptr(ghosted_id_up)
      Mlambda_dn = Mlambda_dn + &
        inversion_ts_aux%dFluxdIntConn(3,iconn)*lambda_ptr(ghosted_id_up)
      Mlambda_dn = Mlambda_dn - &
        inversion_ts_aux%dFluxdIntConn(3,iconn)*lambda_ptr(ghosted_id_dn)
      Mlambda_up = Mlambda_up - &
        inversion_ts_aux%dFluxdIntConn(1,iconn)*lambda_ptr(ghosted_id_dn)
      rhs_up = rhs_up - inversion_ts_aux%dFluxdIntConn(5,iconn)
      rhs_dn = rhs_dn + inversion_ts_aux%dFluxdIntConn(5,iconn)
    elseif (ghosted_id == ghosted_id_dn) then
      Mlambda_up = Mlambda_up + &
        inversion_ts_aux%dFluxdIntConn(2,iconn)*lambda_ptr(ghosted_id_up)
      Mlambda_dn = Mlambda_dn + &
        inversion_ts_aux%dFluxdIntConn(4,iconn)*lambda_ptr(ghosted_id_up)
      Mlambda_dn = Mlambda_dn - &
        inversion_ts_aux%dFluxdIntConn(4,iconn)*lambda_ptr(ghosted_id_dn)
      Mlambda_up = Mlambda_up - &
        inversion_ts_aux%dFluxdIntConn(2,iconn)*lambda_ptr(ghosted_id_dn)
      rhs_up = rhs_up - inversion_ts_aux%dFluxdIntConn(6,iconn)
      rhs_dn = rhs_dn + inversion_ts_aux%dFluxdIntConn(6,iconn)
    else
      this%realization%option%io_buffer = 'Incorrect mapping of connection'
      call PrintErrMsg(this%realization%option)
    endif
    Mlambda(ghosted_id_up) = Mlambda(ghosted_id_up) + Mlambda_up
    Mlambda(ghosted_id_dn) = Mlambda(ghosted_id_dn) + Mlambda_dn
    rhs(ghosted_id_up) = rhs(ghosted_id_up) + rhs_up
    rhs(ghosted_id_dn) = rhs(ghosted_id_dn) + rhs_dn
  enddo

  do i = 1, inversion_aux%cell_to_bc_connection(0,local_id)
    iconn = inversion_aux%cell_to_bc_connection(i,local_id)
    Mlambda(ghosted_id) = Mlambda(ghosted_id) + &
      inversion_ts_aux%dFluxdBCConn(1,iconn)*lambda_ptr(ghosted_id)
    rhs(ghosted_id) = rhs(ghosted_id) + &
      inversion_ts_aux%dFluxdBCConn(2,iconn)
  enddo

  value_ = dot_product(rhs,lambda_ptr) - dot_product(solution,Mlambda)

end subroutine InvSubsrfBMinusSM

! ************************************************************************** !

subroutine InvSubsurfScaleSensitivity(this,JsensitivityT)
  !
  ! Writes sensitivity Jacobian to an ASCII output file
  !
  ! Author: Glenn hammond
  ! Date: 10/11/21
  !
  use Realization_Base_class
  use Variables_module, only : PERMEABILITY

  class(inversion_subsurface_type) :: this
  Mat :: JsensitivityT

  PetscErrorCode :: ierr

  call RealizationGetVariable(this%realization, &
                              this%realization%field%work, &
                              PERMEABILITY,ZERO_INTEGER)
  call MatDiagonalScale(JsensitivityT, &
                        this%realization%field%work, & ! scales rows
                        PETSC_NULL_VEC, &  ! scales columns
                        ierr);CHKERRQ(ierr)

end subroutine InvSubsurfScaleSensitivity

! ************************************************************************** !

subroutine InvSubsurfOutputSensitivity(this,suffix)
  !
  ! Writes sensitivity Jacobian
  !
  ! Author: Glenn hammond
  ! Date: 10/11/21
  !

  class(inversion_subsurface_type) :: this
  character(len=*) :: suffix

  character(len=MAXSTRINGLENGTH) :: filename_prefix

  filename_prefix = trim(this%driver%global_prefix) // '_Jsense'
  if (len_trim(suffix) > 0) filename_prefix = trim(filename_prefix) // '_' // &
                            suffix
  call InvSubsurfOutputSensitivityASCII(this,this%inversion_aux%JsensitivityT, &
                                        filename_prefix)
  call InvSubsurfOutputSensitivityHDF5(this,this%inversion_aux%JsensitivityT, &
                                       filename_prefix)

end subroutine InvSubsurfOutputSensitivity

! ************************************************************************** !

subroutine InvSubsurfOutputSensitivityASCII(this,JsensitivityT,filename_prefix)
  !
  ! Writes sensitivity Jacobian to an ASCII output file
  !
  ! Author: Glenn hammond
  ! Date: 10/11/21
  !
  use Realization_Subsurface_class
  use Variables_module, only : PERMEABILITY

  class(inversion_subsurface_type) :: this
  Mat :: JsensitivityT
  character(len=*) :: filename_prefix

  character(len=MAXSTRINGLENGTH) :: string
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  if (.not.associated(this%realization)) then
    call this%driver%PrintErrMsg('InvSubsurfOutputSensitivityASCII must be &
           &called before the forward simulation is destroyed.')
  endif

  string = trim(filename_prefix) // '.txt'
  call PetscViewerASCIIOpen(this%realization%option%mycomm,string,viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(JsensitivityT,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

end subroutine InvSubsurfOutputSensitivityASCII

! ************************************************************************** !

subroutine InvSubsurfOutputSensitivityHDF5(this,JsensitivityT,filename_prefix)
  !
  ! Writes sensitivity Jacobian to an HDF5 output file
  !
  ! Author: Glenn hammond
  ! Date: 10/18/21
  !
  use hdf5
  use HDF5_module
  use Output_HDF5_module
  use String_module

  class(inversion_subsurface_type) :: this
  Mat :: JsensitivityT
  character(len=*) :: filename_prefix

  Vec :: row_vec
  PetscReal, pointer :: row_ptr(:)
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscInt :: imeasurement
  PetscInt :: num_measurement
  PetscErrorCode :: ierr

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
  character(len=MAXSTRINGLENGTH) :: string
  integer :: hdf5_err

  if (.not.associated(this%realization)) then
    call this%driver%PrintErrMsg('InvSubsurfOutputSensitivityHDF5 must be &
           &called before the forward simulation is destroyed.')
  endif

  !HDF5 formatted output
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,this%realization%option%mycomm, &
                          MPI_INFO_NULL,hdf5_err)
#endif
  string = trim(filename_prefix) // '.h5'
  call h5fcreate_f(string,H5F_ACC_TRUNC_F,file_id,hdf5_err, &
                   H5P_DEFAULT_F,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  call OutputHDF5WriteStructCoordGroup(file_id, &
                                       this%realization%discretization, &
                                       this%realization%patch%grid, &
                                       this%realization%option)
  ! create a group for the data set
  this%realization%option%time = 0.d0
  write(string,'(''Time:'',es13.5,x,a1)') &
        this%realization%option%time/this%realization%output_option%tconv, &
        this%realization%output_option%tunit
  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)

  num_measurement = size(this%imeasurement)
  call VecCreateMPI(this%realization%option%mycomm,PETSC_DECIDE, &
                    num_measurement, &
                    row_vec,ierr);CHKERRQ(ierr)
  do imeasurement = 1, num_measurement
    call VecZeroEntries(row_vec,ierr);CHKERRQ(ierr)
    if (this%realization%option%myrank == 0) then
      call VecSetValue(row_vec,imeasurement-1,1.d0,INSERT_VALUES, &
                       ierr);CHKERRQ(ierr)
    endif
    call VecAssemblyBegin(row_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(row_vec,ierr);CHKERRQ(ierr)
    call MatMult(JsensitivityT,row_vec, &
                 this%realization%field%work, &
                 ierr);CHKERRQ(ierr)
    string = 'Measurement ' // StringWrite(this%imeasurement(imeasurement))
    call HDF5WriteStructDataSetFromVec(string,this%realization, &
                                       this%realization%field%work,grp_id, &
                                       H5T_NATIVE_DOUBLE)
  enddo
  call VecDestroy(row_vec,ierr);CHKERRQ(ierr)
  call h5gclose_f(grp_id,hdf5_err)
  call OutputHDF5CloseFile(this%realization%option,file_id)

end subroutine InvSubsurfOutputSensitivityHDF5

! ************************************************************************** !

subroutine InversionSubsurfaceInvert(this)
  !
  ! Inverts for a parameter update
  !
  ! Author: Glenn hammond
  ! Date: 09/20/21
  !

  class(inversion_subsurface_type) :: this

end subroutine InversionSubsurfaceInvert

! ************************************************************************** !

subroutine InversionSubsurfaceStrip(this)
  !
  ! Deallocates members of inversion Subsurface
  !
  ! Author: Glenn hammond
  ! Date: 09/20/21
  !
  use Utility_module

  class(inversion_subsurface_type) :: this

  PetscErrorCode :: ierr

  call InversionBaseStrip(this)

  nullify(this%realization)
  call InversionAuxDestroy(this%inversion_aux)
  if (associated(this%forward_simulation)) then
    print *, 'Why is forward simulation still associated in &
             &InversionSubSurfStrip?'
    stop
  endif
  nullify(this%forward_simulation)
  if (this%quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%quantity_of_interest,ierr);CHKERRQ(ierr)
  endif
  if (this%ref_quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%ref_quantity_of_interest,ierr);CHKERRQ(ierr)
  endif
  if (this%measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%measurement_vec,ierr);CHKERRQ(ierr)
  endif
  if (this%scatter_global_to_measurement /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(this%scatter_global_to_measurement, &
                           ierr);CHKERRQ(ierr)
  endif

  call DeallocateArray(this%imeasurement)
  call DeallocateArray(this%measurement)

end subroutine InversionSubsurfaceStrip

end module Inversion_Subsurface_class
