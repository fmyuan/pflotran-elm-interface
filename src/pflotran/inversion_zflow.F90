module Inversion_ZFlow_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Aux_module
  use Inversion_Base_class
  use Inversion_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_subsurface_type) :: inversion_zflow_type
    PetscInt :: info_format
    PetscInt :: start_iteration          ! Starting iteration number
    PetscInt :: miniter,maxiter          ! min/max CGLS iterations
    PetscBool :: string_color

    PetscReal :: beta                    ! regularization parameter
    PetscReal :: beta_red_factor         ! beta reduction factor
    PetscReal :: minparam,maxparam       ! min/max paramter value
    PetscReal :: target_chi2             ! target CHI^2 norm
    PetscReal :: current_chi2

    ! For joint inversion
    PetscReal :: alpha_liquid_pressure      ! weight to liquid pressure cost
    PetscReal :: alpha_liquid_saturation    ! weight to saturation cost
    PetscReal :: alpha_solute_concentration ! weight to concentration cost
    PetscReal :: alpha_ert_measurement      ! weight to ERT cost

    ! Cost/objective functions
    PetscReal :: min_phi_red             ! min change in cost function
    PetscReal :: phi_total_0,phi_total
    PetscReal :: phi_data_0,phi_data
    PetscReal :: phi_model_0,phi_model

    ! to check divergence of gamma in CGLS
    PetscBool :: check_gamma_divergence

    ! arrays for CGLS algorithm
    PetscReal, pointer :: b(:)           ! vector for CGLS RHS
    PetscReal, pointer :: p(:)           ! vector of dim -> num of inv cells
    PetscReal, pointer :: q(:)           ! product of Jacobian with p = Jp
    PetscReal, pointer :: r(:)           ! vector of dim -> num of measur
    PetscReal, pointer :: s(:)           ! product Jacobian transpose with r
    PetscReal, pointer :: del_param(:)   ! parameter update vector

    ! For Wm
    PetscInt :: num_constraints_local    ! Number of constraints
    PetscInt :: num_constraints_total    ! Total number of constraints
    PetscInt, pointer :: rblock(:,:)     ! array stores info about reg.
    PetscReal, pointer :: Wm(:)          ! Regularization matrix
    Vec :: parameter_tmp_vec
    Vec :: dist_parameter_tmp_vec

    type(constrained_block_type), pointer :: constrained_block

  contains
    procedure, public :: Init => InversionZFlowInit
    procedure, public :: ReadBlock => InversionZFlowReadBlock
    procedure, public :: SetupForwardRunLinkage => &
                           InvZFlowSetupForwardRunLinkage
    procedure, public :: EvaluateCostFunction => InvZFlowEvaluateCostFunction
    procedure, public :: Checkpoint => InversionZFlowCheckpoint
    procedure, public :: RestartReadData => InversionZFlowRestartReadData
    procedure, public :: CheckConvergence => InversionZFlowCheckConvergence
    procedure, public :: WriteIterationInfo => InvZFlowWriteIterationInfo
    procedure, public :: ScaleSensitivity => InversionZFlowScaleSensitivity
    procedure, public :: CalculateUpdate => InversionZFlowCalculateUpdate
    procedure, public :: UpdateRegularizationParameters => &
                           InvZFlowUpdateRegularizParams
    procedure, public :: Finalize => InversionZFlowFinalize
    procedure, public :: Strip => InversionZFlowStrip
  end type inversion_zflow_type

  type, public :: constrained_block_type
    PetscInt :: num_constrained_block
    PetscInt :: max_num_block_link
    type(constrained_block_par_type), pointer :: constrained_block_list

    ! arrays from the linked list
    character(len=MAXWORDLENGTH), pointer :: material_name(:)
    PetscInt, pointer :: material_id(:)
    PetscInt, pointer :: structure_metric(:)
    PetscInt, pointer :: wf_type(:)
    PetscInt, pointer :: block_link(:,:)
    PetscReal, pointer :: wf_mean(:)
    PetscReal, pointer :: wf_sdev(:)
    PetscReal, pointer :: relative_weight(:)
    PetscReal, pointer :: aniso_weight(:,:)
    PetscReal, pointer :: reference_parameter(:)
  end type constrained_block_type

  type constrained_block_par_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: structure_metric
    PetscInt :: weighing_function
    PetscInt :: num_block_link
    character(len=MAXWORDLENGTH), pointer :: block_link(:)

    PetscReal :: aniso_weight(3)
    PetscReal :: relative_weight
    PetscReal :: weighing_function_mean
    PetscReal :: weighing_function_std_dev
    PetscReal :: reference_parameter
    type(constrained_block_par_type), pointer :: next
  end type constrained_block_par_type

  public :: InversionZFlowCreate, &
            InversionZFlowFinalize, &
            InversionZFlowStrip, &
            InversionZFlowDestroy

contains

! ************************************************************************** !

function InversionZFlowCreate(driver)
  !
  ! Allocates and initializes a new zflow inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/05/22
  !

  use Driver_class

  class(driver_type), pointer :: driver

  class(inversion_zflow_type), pointer :: InversionZFlowCreate

  allocate(InversionZFlowCreate)
  call InversionZFlowCreate%Init(driver)

end function InversionZFlowCreate

! ************************************************************************** !

subroutine InversionZFlowInit(this,driver)
  !
  ! Initializes a new zflow inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/05/22
  !
  use Driver_class

  implicit none

  class(inversion_zflow_type) :: this
  class(driver_type), pointer :: driver

  call InversionSubsurfaceInit(this,driver)

  ! Default inversion parameters
  this%info_format = 0
  this%miniter = 10
  this%maxiter = 50
  this%string_color = PETSC_TRUE

  this%beta = 100.d0
  this%beta_red_factor = 0.5d0
  this%minparam = 1d-17
  this%maxparam = 1d-07
  this%target_chi2 = 1.d0
  this%min_phi_red = 0.2d0

  this%alpha_liquid_pressure = 1.d0
  this%alpha_liquid_saturation = 1.d0
  this%alpha_solute_concentration = 1.d0
  this%alpha_ert_measurement = 1.d0

  this%start_iteration = 1
  this%maximum_iteration = 20
  this%num_constraints_local = 0
  this%num_constraints_total = 0
  this%current_chi2 = UNINITIALIZED_DOUBLE
  this%phi_total_0 = UNINITIALIZED_DOUBLE
  this%phi_data_0 = UNINITIALIZED_DOUBLE
  this%phi_model_0 = UNINITIALIZED_DOUBLE
  this%phi_total = UNINITIALIZED_DOUBLE
  this%phi_data = UNINITIALIZED_DOUBLE
  this%phi_model = UNINITIALIZED_DOUBLE

  this%check_gamma_divergence = PETSC_FALSE

  this%parameter_tmp_vec = PETSC_NULL_VEC
  this%dist_parameter_tmp_vec = PETSC_NULL_VEC

  nullify(this%b)
  nullify(this%p)
  nullify(this%q)
  nullify(this%r)
  nullify(this%s)
  nullify(this%del_param)
  nullify(this%Wm)
  nullify(this%rblock)

  this%constrained_block => ConstrainedBlockCreate()

end subroutine InversionZFlowInit

! ************************************************************************** !

function ConstrainedBlockCreate()
  !
  ! Creates Constrained Block type
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/05/22
  !

  implicit none

  type(constrained_block_type), pointer :: ConstrainedBlockCreate
  type(constrained_block_type), pointer :: constrained_block

  allocate(constrained_block)

  constrained_block%num_constrained_block = 0
  nullify(constrained_block%constrained_block_list)

  nullify(constrained_block%material_name)
  nullify(constrained_block%material_id)
  nullify(constrained_block%structure_metric)
  nullify(constrained_block%wf_type)
  nullify(constrained_block%block_link)
  nullify(constrained_block%wf_mean)
  nullify(constrained_block%wf_sdev)
  nullify(constrained_block%relative_weight)
  nullify(constrained_block%aniso_weight)
  nullify(constrained_block%reference_parameter)

  ConstrainedBlockCreate => constrained_block

end function ConstrainedBlockCreate

! ************************************************************************** !

function ConstrainedBlockParCreate()
  !
  ! Creates Constrained Block Par type
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/05/22
  !

  implicit none

  type(constrained_block_par_type), pointer :: ConstrainedBlockParCreate
  type(constrained_block_par_type), pointer :: constrained_block

  allocate(constrained_block)
  constrained_block%id = 0
  constrained_block%name = ''

  ! Deafult is set the smoothness constraint
  constrained_block%structure_metric = 2
  constrained_block%weighing_function = 1
  constrained_block%num_block_link = 0
  nullify(constrained_block%block_link)

  constrained_block%aniso_weight = 1.d0
  constrained_block%relative_weight = 1.d0
  constrained_block%weighing_function_mean = 10.d0
  constrained_block%weighing_function_std_dev = 0.001d0
  constrained_block%reference_parameter = 0.d0
  nullify(constrained_block%next)

  ConstrainedBlockParCreate => constrained_block

end function ConstrainedBlockParCreate

! ************************************************************************** !

subroutine InversionZFlowAllocateWorkArrays(this)
  !
  ! Initialize inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/05/22
  !

  use Grid_module

  implicit none

  class(inversion_zflow_type) :: this

  PetscInt :: num_measurement
  PetscInt :: num_constraints

  num_measurement = size(this%inversion_aux%measurements)
  num_constraints = this%num_constraints_local

  allocate(this%b(num_measurement + num_constraints))
  allocate(this%p(this%num_parameters_local))
  allocate(this%q(num_measurement + num_constraints))
  allocate(this%r(num_measurement + num_constraints))
  allocate(this%s(this%num_parameters_local))
  allocate(this%del_param(this%num_parameters_local))

  this%b = 0.d0
  this%p = 0.d0
  this%q = 0.d0
  this%r = 0.d0
  this%s = 0.d0
  this%del_param = 0.d0

end subroutine InversionZFlowAllocateWorkArrays

! ************************************************************************** !

subroutine InversionZFlowDeallocateWorkArrays(this)
  !
  ! Initialize inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/05/22
  !

  use Utility_module, only : DeallocateArray

  implicit none

  class(inversion_zflow_type) :: this

  call DeallocateArray(this%b)
  call DeallocateArray(this%p)
  call DeallocateArray(this%q)
  call DeallocateArray(this%r)
  call DeallocateArray(this%s)
  call DeallocateArray(this%del_param)

end subroutine InversionZFlowDeallocateWorkArrays

! ************************************************************************** !

subroutine InversionZFlowConstrainedArraysFromList(this)
  !
  ! Gets Constrained Block parameter arrays from linked list
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/05/22
  !

  use Option_module
  use Material_module
  use Patch_module

  implicit none

  class(inversion_zflow_type) :: this

  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(material_property_type), pointer :: material_property
  type(constrained_block_type), pointer :: constrained_block
  type(constrained_block_par_type), pointer :: cur_constrained_block

  PetscInt :: i,iconblock
  PetscInt :: nconblock

  patch => this%realization%patch
  option => this%realization%option

  constrained_block => this%constrained_block
  nconblock = constrained_block%num_constrained_block

  if (nconblock > 0) then
    allocate(constrained_block%material_name(nconblock))
    constrained_block%material_name = ''
    allocate(constrained_block%material_id(nconblock))
    constrained_block%material_id = 0
    allocate(constrained_block%structure_metric(nconblock))
    constrained_block%structure_metric = 2
    allocate(constrained_block%wf_type(nconblock))
    constrained_block%wf_type = 1
    allocate(constrained_block%wf_mean(nconblock))
    constrained_block%wf_mean = 10.d0
    allocate(constrained_block%wf_sdev(nconblock))
    constrained_block%wf_sdev = 0.001d0
    allocate(constrained_block%relative_weight(nconblock))
    constrained_block%relative_weight = 1.d0
    allocate(constrained_block%aniso_weight(nconblock,THREE_INTEGER))
    constrained_block%aniso_weight = 1.d0
    allocate(constrained_block%reference_parameter(nconblock))
    constrained_block%reference_parameter = 0.d0
    allocate(constrained_block%block_link(nconblock,&
             constrained_block%max_num_block_link+1))
    constrained_block%block_link = 0

    cur_constrained_block => constrained_block%constrained_block_list
    iconblock = 1
    do
      if (.not. associated(cur_constrained_block)) exit
      constrained_block%material_name(iconblock) = cur_constrained_block%name
      material_property => &
          MaterialPropGetPtrFromArray(cur_constrained_block%name, &
                                      patch%material_property_array)
      if (.not.associated(material_property)) then
        option%io_buffer = 'Constrained block " &
                           &' // trim(cur_constrained_block%name) // &
                           &'" not found in material list'
        call PrintErrMsg(option)
      endif

      constrained_block%material_id(iconblock) = material_property%internal_id
      constrained_block%structure_metric(iconblock) = &
                        cur_constrained_block%structure_metric
      constrained_block%wf_type(iconblock) = &
                        cur_constrained_block%weighing_function
      constrained_block%wf_mean(iconblock) = &
                        cur_constrained_block%weighing_function_mean
      constrained_block%wf_sdev(iconblock) = &
                        cur_constrained_block%weighing_function_std_dev
      constrained_block%relative_weight(iconblock) = &
                        cur_constrained_block%relative_weight
      constrained_block%aniso_weight(iconblock,:) = &
                        cur_constrained_block%aniso_weight
      constrained_block%reference_parameter(iconblock) = &
                        cur_constrained_block%reference_parameter
      constrained_block%block_link(iconblock,ONE_INTEGER) = &
                        cur_constrained_block%num_block_link
      do i=1,cur_constrained_block%num_block_link
        material_property => &
            MaterialPropGetPtrFromArray(cur_constrained_block%block_link(i), &
                                      patch%material_property_array)
        if (.not.associated(material_property)) then
          option%io_buffer = 'Linked block "&
                             &'//trim(cur_constrained_block%block_link(i)) // &
                             &'" in constrained block "&
                             &'//trim(cur_constrained_block%name) // &
                             &'" not found in material list'
          call PrintErrMsg(option)
        endif
        constrained_block%block_link(iconblock,i+1) = &
                                      material_property%internal_id
      enddo

      cur_constrained_block => cur_constrained_block%next
      iconblock = iconblock + 1

    enddo
  endif

end subroutine InversionZFlowConstrainedArraysFromList

! ************************************************************************** !

subroutine InversionZFlowReadBlock(this,input,option)
  !
  ! Reads input file parameters associated an ZFlow inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/06/22
  !

  use Input_Aux_module
  use Option_module
  use String_module

  class(inversion_zflow_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'ZFlow Inversion'

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
      case('MIN_PERMEABILITY','MIN_PARAMETER')
        call InputReadDouble(input,option,this%minparam)
        call InputErrorMsg(input,option,'MIN_PARAMETER', &
                           error_string)
      case('MAX_PERMEABILITY','MAX_PARAMETER')
        call InputReadDouble(input,option,this%maxparam)
        call InputErrorMsg(input,option,'MAX_PARAMETER', &
                           error_string)
      case('MIN_CGLS_ITERATION')
        call InputReadInt(input,option,this%miniter)
        call InputErrorMsg(input,option,'MIN_CGLS_ITERATION',error_string)
      case('MAX_CGLS_ITERATION')
        call InputReadInt(input,option,this%maxiter)
        call InputErrorMsg(input,option,'MAX_CGLS_ITERATION',error_string)
      case('CHECK_CGLS_GAMMA_DIVERGENCE')
        this%check_gamma_divergence = PETSC_TRUE
      case('BETA')
        call InputReadDouble(input,option,this%beta)
        call InputErrorMsg(input,option,'BETA',error_string)
      case('BETA_REDUCTION_FACTOR')
        call InputReadDouble(input,option,this%beta_red_factor)
        call InputErrorMsg(input,option,'BETA_REDUCTION_FACTOR',error_string)
      case('ALPHA_LIQUID_PRESSURE')
        call InputReadDouble(input,option,this%alpha_liquid_pressure)
        call InputErrorMsg(input,option,'ALPHA_LIQUID_PRESSURE',error_string)
      case('ALPHA_LIQUID_SATURATION')
        call InputReadDouble(input,option,this%alpha_liquid_saturation)
        call InputErrorMsg(input,option,'ALPHA_LIQUID_SATURATION',error_string)
      case('ALPHA_SOLUTE_CONCENTRATION')
        call InputReadDouble(input,option,this%alpha_solute_concentration)
        call InputErrorMsg(input,option,'ALPHA_SOLUTE_CONCENTRATION', &
                           error_string)
      case('ALPHA_ERT_MEASUREMENT')
        call InputReadDouble(input,option,this%alpha_ert_measurement)
        call InputErrorMsg(input,option,'ALPHA_ERT_MEASUREMENT',error_string)
      case('TARGET_CHI2')
        call InputReadDouble(input,option,this%target_chi2)
        call InputErrorMsg(input,option,'TARGET_CHI2',error_string)
      case('MIN_COST_REDUCTION')
        call InputReadDouble(input,option,this%min_phi_red)
        call InputErrorMsg(input,option,'MIN_COST_REDUCTION',error_string)
      case('CONSTRAINED_BLOCKS')
        call ConstrainedBlockRead(this%constrained_block,input,option)
      case('NO_STRING_COLOR')
        this%string_color = PETSC_FALSE
      case('INFO_FORMAT_GLENN')
        this%info_format = 1
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine InversionZFlowReadBlock

! ************************************************************************** !

subroutine ConstrainedBlockRead(constrained_block,input,option)
  !
  ! Read constrained blocks options
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/06/22

  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(constrained_block_type) :: constrained_block
  type(input_type), pointer :: input
  type(option_type) :: option

  type(constrained_block_par_type), pointer :: cur_constrained_block
  type(constrained_block_par_type), pointer :: prev_constrained_block

  character(len=MAXWORDLENGTH) :: error_string

  error_string = 'INVERSION,CONSTRAINED_BLOCKS'

  nullify(prev_constrained_block)
  call InputPushBlock(input,option)
  constrained_block%max_num_block_link = 0
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    constrained_block%num_constrained_block = &
                          constrained_block%num_constrained_block + 1

    cur_constrained_block => ConstrainedBlockParCreate()
    call InputReadCard(input,option,cur_constrained_block%name)
    call InputErrorMsg(input,option,'keyword',error_string)

    call ConstrainedBlockParRead(cur_constrained_block,input,option)
    if (constrained_block%max_num_block_link < &
        cur_constrained_block%num_block_link) &
      constrained_block%max_num_block_link = &
        cur_constrained_block%num_block_link

    if (.not.associated(constrained_block%constrained_block_list)) then
      constrained_block%constrained_block_list => cur_constrained_block
      cur_constrained_block%id = 1
    endif
    if (associated(prev_constrained_block)) then
      prev_constrained_block%next => cur_constrained_block
      cur_constrained_block%id = prev_constrained_block%id + 1
    endif
    prev_constrained_block => cur_constrained_block
    nullify(cur_constrained_block)
  enddo
  call InputPopBlock(input,option)

end subroutine ConstrainedBlockRead

! ************************************************************************** !

subroutine ConstrainedBlockParRead(constrained_block,input,option)
  !
  ! Read constrained blocks parameters options
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/06/22

  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(constrained_block_par_type) :: constrained_block
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i,num_block_link
  PetscReal :: norm_factor
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: error_string

  error_string = 'INVERSION,CONSTRAINED_BLOCKS'

  word = ''

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(trim(word))
      case('STRUCTURE_METRIC')
        call InputReadInt(input,option,constrained_block%structure_metric)
        call InputErrorMsg(input,option,'STRUCTURE_METRIC',error_string)
      case('WEIGHING_FUNCTION')
        call InputReadInt(input,option,constrained_block%weighing_function)
        call InputErrorMsg(input,option,'WEIGHING_FUNCTION',error_string)
      case('WEIGHING_FUNCTION_MEAN')
        call InputReadDouble(input,option, &
                            constrained_block%weighing_function_mean)
        call InputErrorMsg(input,option,'WEIGHING_FUNCTION_MEAN',error_string)
      case('WEIGHING_FUNCTION_STD_DEVIATION')
        call InputReadDouble(input,option, &
                            constrained_block%weighing_function_std_dev)
        call InputErrorMsg(input,option,'WEIGHING_FUNCTION_STD_DEVIATION', &
                          error_string)
      case('BLOCK_LINKS')
        call InputReadInt(input,option,num_block_link)
        call InputErrorMsg(input,option,'BLOCK_LINKS',error_string)
        constrained_block%num_block_link = num_block_link
        allocate(constrained_block%block_link(num_block_link))
        do i=1,num_block_link
          call InputReadCard(input,option,constrained_block%block_link(i))
          call InputErrorMsg(input,option,'BLOCK_LINKS',error_string)
        enddo
      case('ANISOTROPIC_WEIGHTS')
        do i=1,THREE_INTEGER
          call InputReadDouble(input,option,constrained_block%aniso_weight(i))
          call InputErrorMsg(input,option,'ANISOTROPY_WEIGHTS',error_string)
        enddo
        norm_factor = norm2(constrained_block%aniso_weight)
        if (norm_factor > 0.) constrained_block%aniso_weight = &
                                constrained_block%aniso_weight / norm_factor
      case('RELATIVE_WEIGHT')
        call InputReadDouble(input,option,constrained_block%relative_weight)
        call InputErrorMsg(input,option,'RELATIVE_WEIGHT',error_string)
      case('REFERENCE_PERMEABILITY','REFERENCE_PARAMETER')
        call InputReadDouble(input,option, &
                            constrained_block%reference_parameter)
        call InputErrorMsg(input,option,'REFERENCE_PARAMETER',error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine ConstrainedBlockParRead

! ************************************************************************** !

function ConstrainedBlockGetBlockIDFromMatID(constrained_block,material_id, &
                                             option)
  !
  ! Read constrained blocks parameters options
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 07/18/22

  use Option_module

  type(constrained_block_type) :: constrained_block
  PetscInt :: material_id
  PetscInt :: ConstrainedBlockGetBlockIDFromMatID

  type(option_type), pointer :: option

  PetscInt :: iblock
  PetscBool :: found

  found = .false.

  do iblock = 1,constrained_block%num_constrained_block
    if (material_id == constrained_block%material_id(iblock)) then
      ConstrainedBlockGetBlockIDFromMatID = iblock
      found = .true.
      return
    endif
  enddo

  if (.not.found) then
      option%io_buffer = 'ConstrainedBlockGetBlockIDFromMatID failed to find &
                         &Block ID as specified Material ID is not &
                         &associated with any constrained blocks.'
      call PrintErrMsg(option)
    endif

end function ConstrainedBlockGetBlockIDFromMatID

! ************************************************************************** !

subroutine InvZFlowSetupForwardRunLinkage(this)
  !
  ! Initializes inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/06/22
  !
  use Discretization_module
  use Inversion_TS_Aux_module
  use Inversion_Measurement_Aux_module
  use Inversion_Parameter_module
  use Option_module
  use Variables_module, only : PERMEABILITY,POROSITY,ELECTRICAL_CONDUCTIVITY, &
                               VG_ALPHA,VG_SR,VG_M

  implicit none

  class(inversion_zflow_type) :: this

  PetscBool :: exists
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: iqoi(2)
  PetscInt :: i
  PetscErrorCode :: ierr

  call InvSubsurfSetupForwardRunLinkage(this)

  call VecDuplicate(this%inversion_aux%dist_parameter_vec,this%dist_parameter_tmp_vec, &
                    ierr);CHKERRQ(ierr)

  ! check to ensure that quantity of interest exists
  exists = PETSC_FALSE
  iqoi = InversionParameterIntToQOIArray(this%inversion_aux%parameters(1))
  select case(iqoi(1))
    case(PERMEABILITY)
      if (this%realization%option%iflowmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'PERMEABILITY'
    case(POROSITY)
      if (this%realization%option%iflowmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'POROSITY'
    case(ELECTRICAL_CONDUCTIVITY)
      if (this%realization%option%igeopmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'ELECTRICAL_CONDUCTIVITY'
    case(VG_ALPHA)
      if (this%realization%option%iflowmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'VG_ALPHA'
    case(VG_SR)
      if (this%realization%option%iflowmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'VG_SR'
    case(VG_M)
      if (this%realization%option%iflowmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'VG_M'
    case default
      word = 'unknown_parameter'
  end select
  if (.not.exists) then
    this%realization%option%io_buffer = 'Inversion for ' // trim(word) // &
      &' cannot be performed with the specified process models.'
    call PrintErrMsg(this%realization%option)
  endif

  call InversionZFlowConstrainedArraysFromList(this)

  ! scale data weight by a scalar weight for joint inversion
  if (this%iteration==1) then
    do i = 1, size(this%inversion_aux%measurements)
      select case(this%inversion_aux%measurements(i)%iobs_var)
        case(OBS_LIQUID_PRESSURE)
          this%inversion_aux%measurements(i)%weight = &
            this%alpha_liquid_pressure * &
            this%inversion_aux%measurements(i)%weight
        case(OBS_LIQUID_SATURATION)
          this%inversion_aux%measurements(i)%weight = &
            this%alpha_liquid_saturation * &
            this%inversion_aux%measurements(i)%weight
        case(OBS_SOLUTE_CONCENTRATION)
          this%inversion_aux%measurements(i)%weight = &
            this%alpha_solute_concentration * &
            this%inversion_aux%measurements(i)%weight
        case(OBS_ERT_MEASUREMENT)
          this%inversion_aux%measurements(i)%weight = &
            this%alpha_ert_measurement * &
            this%inversion_aux%measurements(i)%weight
      end select
    enddo
  endif

  ! Build Wm matrix
  call InversionZFlowBuildWm(this)

  if (.not.this%inversion_aux%qoi_is_full_vector) then
    call VecDuplicate(this%inversion_aux%parameter_vec, &
                      this%parameter_tmp_vec, &
                      ierr);CHKERRQ(ierr)
  endif

end subroutine InvZFlowSetupForwardRunLinkage

! ************************************************************************** !

subroutine InversionZFlowCheckConvergence(this)
  !
  ! Check Inversion convergence
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/15/21

  implicit none

  class(inversion_zflow_type) :: this

  this%converged = PETSC_FALSE
  call this%EvaluateCostFunction()
  if ((this%current_chi2 <= this%target_chi2) .or. &
      (this%iteration > this%maximum_iteration)) this%converged = PETSC_TRUE

end subroutine InversionZFlowCheckConvergence

! ************************************************************************** !

subroutine InvZFlowEvaluateCostFunction(this)
  !
  ! Evaluates cost functions for inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/15/21

  use Realization_Base_class
  use Option_module
  use Patch_module
  use Material_Aux_module
  use String_module
  use Variables_module, only : PERMEABILITY,POROSITY,ELECTRICAL_CONDUCTIVITY

  implicit none

  class(inversion_zflow_type) :: this

  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: idata,num_measurement
  PetscInt :: iconst,num_constraints
  PetscInt :: irb,ghosted_id,ghosted_id_nb
  PetscInt :: imat_id,imat_id_nb
  PetscInt :: iparameter
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: wd,tempreal
  PetscReal :: param_ce,param_nb              ! cell's and neighbor's
  PetscReal :: wm,x
  PetscReal, allocatable :: model_vector(:)
  PetscBool :: use_neighbor
  PetscErrorCode :: ierr

  option => this%realization%option
  patch => this%realization%patch
  material_auxvars => patch%aux%Material%auxvars

  constrained_block => this%constrained_block
  rblock => this%rblock

  num_measurement = size(this%inversion_aux%measurements)

  ! Data part
  this%phi_data = 0.d0
  do idata=1,num_measurement
    wd = this%inversion_aux%measurements(idata)%weight
    tempreal = wd * (this%inversion_aux%measurements(idata)%value - &
                     this%inversion_aux%measurements(idata)%simulated_value)
    this%phi_data = this%phi_data + tempreal * tempreal

  enddo

  this%current_chi2 = this%phi_data / num_measurement

  ! model cost function
  this%phi_model = 0.d0

  if (this%inversion_aux%qoi_is_full_vector) then

    num_constraints = this%num_constraints_local
    ! allocate to at least size 1 to allow for inner product
    allocate(model_vector(max(num_constraints,1)))
    model_vector = 0.d0

    do iconst=1,num_constraints
      if (this%Wm(iconst) == 0) cycle

      wm = this%Wm(iconst)

      use_neighbor = PETSC_FALSE
      ghosted_id = rblock(iconst,1)
      if (patch%imat(ghosted_id) <= 0) cycle
      irb = rblock(iconst,3)
      select case(constrained_block%structure_metric(irb))
        case(1:2,5:10)
          ghosted_id_nb = rblock(iconst,2)
          if (patch%imat(ghosted_id_nb) <=0 ) cycle
          use_neighbor = PETSC_TRUE
        case default
      end select

      iparameter = this%inversion_aux%parameters(1)%iparameter

      call InvAuxGetParamValueByCell(this%inversion_aux,param_ce, &
                                     iparameter, &
                                     patch%imat(ghosted_id), &
                                     material_auxvars(ghosted_id))
      if (use_neighbor) then
        call InvAuxGetParamValueByCell(this%inversion_aux,param_nb, &
                                       iparameter, &
                                       patch%imat(ghosted_id_nb), &
                                       material_auxvars(ghosted_id_nb))
      endif

      x = 0.d0

      select case(constrained_block%structure_metric(irb))
        case(1)
          x = log(param_ce) - log(param_nb)
        case(2)
          x = abs(log(param_ce) - log(param_nb))
        case(3)
          x = log(param_ce) - log(constrained_block%reference_parameter(irb))
        case(4)
          x = abs(log(param_ce) - &
                  log(constrained_block%reference_parameter(irb)))
        case(5)
          x = log(param_ce) - log(param_nb)
        case(6)
          x = abs(log(param_ce) - log(param_nb))
        case(7)
          x = (log(param_ce) - log(constrained_block%reference_parameter(irb)))&
             -(log(param_nb) - log(constrained_block%reference_parameter(irb)))
        case(8)
          x = abs( &
              (log(param_ce) - log(constrained_block%reference_parameter(irb)))&
             -(log(param_nb) - log(constrained_block%reference_parameter(irb))))
        case(9)
          x = log(param_ce) - log(param_nb)
        case(10)
          x = abs(log(param_ce) - log(param_nb))
        case default

      end select

      model_vector(iconst) = wm * x

    enddo

    this%phi_model = this%beta * dot_product(model_vector,model_vector)
    call MPI_Allreduce(MPI_IN_PLACE,this%phi_model,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                       ierr);CHKERRQ(ierr)
    deallocate(model_vector)

  else

    num_constraints = this%num_constraints_local
    ! allocate to at least size 1 to allow for inner product
    allocate(model_vector(max(num_constraints,1)))
    model_vector = 0.d0

    do iconst=1,num_constraints
      if (this%Wm(iconst) == 0) cycle

      wm = this%Wm(iconst)

      use_neighbor = PETSC_FALSE
      imat_id = rblock(iconst,1)
      if (imat_id <= 0) cycle
      irb = rblock(iconst,3)
      select case(constrained_block%structure_metric(irb))
        case(1:2,7:10)
          imat_id_nb = rblock(iconst,2)
          if (imat_id_nb <=0 ) cycle
          use_neighbor = PETSC_TRUE
        case(5:6)
          option%io_buffer = 'Supported STRUCTURE_METRIC in INVERSION, &
                             &CONSTRAINED_BLOCKS is between 1-4 and 7-10 for &
                             &parameter based approach.'
          call PrintErrMsg(option)
        case default
      end select

      iparameter = this%inversion_aux%parameters(1)%iparameter

      call InvAuxGetSetParamValueByMat(this%inversion_aux,param_ce, &
                                       iparameter, &
                                       imat_id,INVAUX_GET_MATERIAL_VALUE)
      if (use_neighbor) then
        call InvAuxGetSetParamValueByMat(this%inversion_aux,param_nb, &
                                         iparameter, &
                                         imat_id_nb, &
                                         INVAUX_GET_MATERIAL_VALUE)
      endif

      x = 0.d0

      select case(constrained_block%structure_metric(irb))
        case(1)
          x = log(param_ce) - log(param_nb)
        case(2)
          x = abs(log(param_ce) - log(param_nb))
        case(3)
          x = log(param_ce) - log(constrained_block%reference_parameter(irb))
        case(4)
          x = abs(log(param_ce) - &
                  log(constrained_block%reference_parameter(irb)))
        case(5)
          x = log(param_ce) - log(param_nb)
        case(6)
          x = abs(log(param_ce) - log(param_nb))
        case(7)
          x = (log(param_ce) - log(constrained_block%reference_parameter(irb)))&
             -(log(param_nb) - log(constrained_block%reference_parameter(irb)))
        case(8)
          x = abs( &
              (log(param_ce) - log(constrained_block%reference_parameter(irb)))&
             -(log(param_nb) - log(constrained_block%reference_parameter(irb))))
        case(9)
          x = log(param_ce) - log(param_nb)
        case(10)
          x = abs(log(param_ce) - log(param_nb))
        case default

      end select

      model_vector(iconst) = wm * x

    enddo

    this%phi_model = this%beta * dot_product(model_vector,model_vector)
    deallocate(model_vector)

  endif

  this%phi_total = this%phi_data + this%phi_model

  if (this%iteration == this%start_iteration) then
    this%phi_data_0 = this%phi_data
    this%phi_model_0 = this%phi_model
    this%phi_total_0 = this%phi_total
  endif

end subroutine InvZFlowEvaluateCostFunction

! ************************************************************************** !

subroutine InvZFlowUpdateRegularizParams(this)
  !
  ! Check Beta if it needs cooling/reduction
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/18/21

  implicit none

  class(inversion_zflow_type) :: this

  if (this%iteration == this%start_iteration) return

  if ( (this%phi_total_0 - this%phi_total)/this%phi_total_0 <= &
                                                      this%min_phi_red ) then
    this%beta = this%beta * this%beta_red_factor
    this%phi_model = this%beta_red_factor * this%phi_model
  endif

  ! update the cost functions
  this%phi_data_0 = this%phi_data
  this%phi_model_0 = this%phi_model
  this%phi_total_0 = this%phi_total

end subroutine InvZFlowUpdateRegularizParams

! ************************************************************************** !

subroutine InversionZFlowCalculateUpdate(this)
  !
  ! Calculates updated model parameters
  ! using m_new = m_old + del_m
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/13/21
  !
  use Discretization_module
  use Patch_module
  use Grid_module

  implicit none

  class(inversion_zflow_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid

  PetscInt :: iparameter, ghosted_id
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec2_ptr(:)
  PetscReal :: new_value
  Vec :: work_dup
  Vec :: dist_del_param_vec
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid

  ! simply setting a local pointer for clarity
  dist_del_param_vec = this%dist_parameter_tmp_vec

  call InversionZFlowAllocateWorkArrays(this)

  ! get inversion%del_param
  call InversionZFlowCGLSSolve(this)

  call VecGetArrayF90(dist_del_param_vec,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(:) = this%del_param(:)
  call VecRestoreArrayF90(dist_del_param_vec,vec_ptr,ierr);CHKERRQ(ierr)

  if (this%inversion_aux%qoi_is_full_vector) then
    ! have to copy values to global work vecs in order to loop over
    ! ghosted ids

    ! dist_del_param_vec holds the update
    call InvAuxScatGlobalToDistParam(this%inversion_aux, &
                                     this%realization%field%work, &
                                     dist_del_param_vec, &
                                     INVAUX_SCATREVERSE)
    call VecDuplicate(this%realization%field%work,work_dup, &
                      ierr);CHKERRQ(ierr)
    ! dist_parameter_vec holds the original value
    call InvAuxScatGlobalToDistParam(this%inversion_aux, &
                                     work_dup, &
                                     this%inversion_aux%dist_parameter_vec, &
                                     INVAUX_SCATREVERSE)

    ! Get updated parameter as m_new = m_old + del_m (where m = log(param))
    call VecGetArrayF90(work_dup,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%realization%field%work,vec2_ptr, &
                        ierr);CHKERRQ(ierr)
    do iparameter = 1, this%num_parameters_local
      if (this%inversion_aux%qoi_is_full_vector) then
        ghosted_id = grid%nL2G(iparameter)
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      vec_ptr(iparameter) = exp(log(vec_ptr(iparameter)) + vec2_ptr(iparameter))
      if (vec_ptr(iparameter) > this%maxparam) &
        vec_ptr(iparameter) = this%maxparam
      if (vec_ptr(iparameter) < this%minparam) &
        vec_ptr(iparameter) = this%minparam
    enddo
    call VecRestoreArrayF90(work_dup,vec_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(this%realization%field%work,vec2_ptr, &
                            ierr);CHKERRQ(ierr)
    call InversionZFlowDeallocateWorkArrays(this)

    ! copy back to dist_parameter_vec
    call InvAuxScatGlobalToDistParam(this%inversion_aux, &
                                     work_dup, &
                                     this%inversion_aux%dist_parameter_vec, &
                                     INVAUX_SCATFORWARD)
    call VecDestroy(work_dup,ierr);CHKERRQ(ierr)
  else
    call VecGetArrayF90(this%inversion_aux%dist_parameter_vec, &
                        vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(dist_del_param_vec,vec2_ptr,ierr);CHKERRQ(ierr)
    do iparameter = 1, this%num_parameters_local
#if 0
      vec_ptr(iparameter) = exp(log(vec_ptr(iparameter)) + vec2_ptr(iparameter))
      if (vec_ptr(iparameter) > this%maxparam) &
        vec_ptr(iparameter) = this%maxparam
      if (vec_ptr(iparameter) < this%minparam) &
        vec_ptr(iparameter) = this%minparam
#else
      new_value = exp(log(vec_ptr(iparameter)) + vec2_ptr(iparameter))
      if (new_value > this%maxparam) then
        new_value = this%maxparam
      else if (new_value < this%minparam) then
        new_value = this%minparam
      endif
      vec2_ptr(iparameter) = new_value - vec_ptr(iparameter)
      vec_ptr(iparameter) = new_value
#endif
    enddo
    call VecRestoreArrayF90(this%inversion_aux%dist_parameter_vec,vec_ptr, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(dist_del_param_vec,vec2_ptr,ierr);CHKERRQ(ierr)

    call InvAuxScatParamToDistParam(this%inversion_aux, &
                                    this%inversion_aux%parameter_vec, &
                                    this%inversion_aux%dist_parameter_vec, &
                                    INVAUX_SCATREVERSE)
    call InvAuxScatParamToDistParam(this%inversion_aux, &
                                    this%inversion_aux%del_parameter_vec, &
                                    dist_del_param_vec, &
                                    INVAUX_SCATREVERSE)
    call InvSubsurfPrintCurParamUpdate(this)
  endif

end subroutine InversionZFlowCalculateUpdate

! ************************************************************************** !

subroutine InversionZFlowCGLSSolve(this)
  !
  ! Implements CGLS solver for least sqaure equivalent
  !            of the ZFlow normal equations
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/14/21
  !

  use Option_module
  use Timer_class
  use String_module

  implicit none

  class(inversion_zflow_type) :: this

  type(option_type), pointer :: option
  class(timer_type), pointer :: timer

  PetscInt :: i,nm,ncons
  PetscReal :: alpha,gbeta,gamma,gamma1,delta1,delta2,delta
  PetscReal :: norms0,norms,normx,xmax
  PetscReal :: resNE,resNE_old
  PetscBool :: exit_info,indefinite
  PetscErrorCode :: ierr

  PetscReal, parameter :: delta_initer = 1e-23
  PetscReal, parameter :: initer_conv  = 1e-24

  option => this%realization%option

  this%del_param = 0.0d0

  alpha = 0.d0
  gbeta = 0.d0
  gamma = 0.d0
  delta = 0.d0
  gamma1 = 0.d0
  delta1 = 0.d0
  delta2 = 0.d0

  timer => TimerCreate()
  call timer%Start()

  if (OptionPrintToScreen(option)) then
    write(*,'(" --> Solving ZFlow normal equation using CGLS solver:")')
  endif

  nm = size(this%inversion_aux%measurements)
  ncons = this%num_constraints_local

  ! Get RHS vector this%b
  call InversionZFlowCGLSRhs(this)

  this%r = this%b

  ! get this%s = J^tr
  call InversionZFlowComputeMatVecProductJtr(this)
  this%p = this%s

  gamma = dot_product(this%s,this%s)
  call MPI_Allreduce(MPI_IN_PLACE,gamma,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)

  norms0 = sqrt(gamma)
  xmax = 0.d0
  normx = 0.d0
  resNE = 0.d0
  exit_info = PETSC_FALSE
  indefinite = PETSC_FALSE

  do i=1,this%maxiter

    if (exit_info) exit

    ! get this%q = Jp
    call InversionZFlowComputeMatVecProductJp(this)

    delta1 = dot_product(this%q(1:nm),this%q(1:nm))
    if (ncons > 0) &
      delta2 = dot_product(this%q(nm+1:nm+ncons),this%q(nm+1:nm+ncons))

    if (this%inversion_aux%qoi_is_full_vector) &
      call MPI_Allreduce(MPI_IN_PLACE,delta2,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                         ierr);CHKERRQ(ierr)
    delta = delta1 + delta2

    if (delta < 0) indefinite = PETSC_TRUE
    if (delta == 0) delta = epsilon(delta)

    alpha = gamma / delta

    this%del_param = this%del_param + alpha * this%p
    this%r = this%r - alpha * this%q

    ! get this%s = J^tr
    call InversionZFlowComputeMatVecProductJtr(this)

    gamma1 = gamma
    gamma = dot_product(this%s,this%s)
    call MPI_Allreduce(MPI_IN_PLACE,gamma,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                       ierr);CHKERRQ(ierr)

    norms = sqrt(gamma)
    gbeta = gamma / gamma1
    this%p = this%s + gbeta * this%p

    normx = dot_product(this%del_param,this%del_param)
    call MPI_Allreduce(MPI_IN_PLACE,normx,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                       ierr);CHKERRQ(ierr)
    normx = sqrt(normx)
    if (xmax < normx) xmax = normx
    if ( (norms <= norms0 * initer_conv) .or. (normx * initer_conv >= 1)) &
                               exit_info = PETSC_TRUE

    resNE_old = resNE
    resNE = norms / norms0

    if ( abs((resNE_old - resNe) /resNE_old) < delta_initer .and. &
        i > this%miniter) exit_info = PETSC_TRUE

    ! PJ: for coupled flow, transport, and ert we need following condition --> ?
    if (this%check_gamma_divergence .or. &
        this%inversion_option%coupled_flow_ert) then
      if (gamma > gamma1) exit_info = PETSC_TRUE
    endif

  enddo

  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime())) &
    // ' seconds and ' // trim(StringWrite(i)) // &
    ' iterations to solve ZFlow normal equation.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InversionZFlowCGLSSolve

! ************************************************************************** !

subroutine InversionZFlowCGLSRhs(this)
  !
  ! Builds RHS for least-square equation for CGLS solver
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/23/21
  !

  use Realization_Base_class
  use Patch_module
  use Material_Aux_module
  use Option_module
  use String_module
  use Variables_module, only : PERMEABILITY,POROSITY,ELECTRICAL_CONDUCTIVITY

  implicit none

  class(inversion_zflow_type) :: this

  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: idata,iconst,irb,num_measurement
  PetscInt :: ghosted_id,ghosted_id_nb
  PetscInt :: imat_id,imat_id_nb
  PetscInt :: iparameter
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: param_ce,param_nb,x     ! cell's and neighbor's
  PetscReal :: wm,beta
  PetscReal :: wd
  PetscBool :: use_neighbor

  option => this%realization%option
  patch => this%realization%patch
  material_auxvars => patch%aux%Material%auxvars

  constrained_block => this%constrained_block
  rblock => this%rblock

  this%b = 0.0d0

  num_measurement = size(this%inversion_aux%measurements)

  ! Data part
  do idata=1,num_measurement
    wd = this%inversion_aux%measurements(idata)%weight
    this%b(idata) = wd * (this%inversion_aux%measurements(idata)%value - &
                        this%inversion_aux%measurements(idata)%simulated_value)
  enddo

  ! Model part
  beta = this%beta

  if (this%inversion_aux%qoi_is_full_vector) then

    do iconst=1,this%num_constraints_local
      if (this%Wm(iconst) == 0) cycle

      wm = this%Wm(iconst)

      use_neighbor = PETSC_FALSE
      ghosted_id = rblock(iconst,1)
      if (patch%imat(ghosted_id) <= 0) cycle
      irb = rblock(iconst,3)
      select case(constrained_block%structure_metric(irb))
        case(1:2,5:10)
          ghosted_id_nb = rblock(iconst,2)
          if (patch%imat(ghosted_id_nb) <=0 ) cycle
          use_neighbor = PETSC_TRUE
        case default
      end select

      iparameter = this%inversion_aux%parameters(1)%iparameter

      call InvAuxGetParamValueByCell(this%inversion_aux,param_ce, &
                                     iparameter, &
                                     patch%imat(ghosted_id), &
                                     material_auxvars(ghosted_id))
      if (use_neighbor) then
        call InvAuxGetParamValueByCell(this%inversion_aux,param_nb, &
                                       iparameter, &
                                       patch%imat(ghosted_id_nb), &
                                       material_auxvars(ghosted_id_nb))
      endif

      x = 0.0d0

      select case(constrained_block%structure_metric(irb))
        case(1)
          x = log(param_ce) - log(param_nb)
        case(2)
          x = log(param_ce) - log(param_nb)
        case(3)
          x = log(param_ce) - log(constrained_block%reference_parameter(irb))
        case(4)
          x = log(param_ce) - log(constrained_block%reference_parameter(irb))
        case(5)
          x = log(param_ce) - log(param_nb)
          ! TODO: compute rx,ry, and rz
        case(6)
          x = log(param_ce) - log(param_nb)
          ! TODO: compute rx,ry, and rz
        case(7)
          x = (log(param_ce) - log(constrained_block%reference_parameter(irb)))&
             -(log(param_nb) - log(constrained_block%reference_parameter(irb)))
        case(8)
          x = (log(param_ce) - log(constrained_block%reference_parameter(irb)))&
             -(log(param_nb) - log(constrained_block%reference_parameter(irb)))
        case(9)
          x = log(param_ce) - log(param_nb)
        case(10)
          x = log(param_ce) - log(param_nb)
        case default
          option%io_buffer = 'Supported STRUCTURE_METRIC in INVERSION, &
                              &CONSTRAINED_BLOCKS is between 1 to 10'
          call PrintErrMsg(option)
      end select

      this%b(num_measurement + iconst) = - sqrt(beta) * wm * x

    enddo

  else

    do iconst=1,this%num_constraints_local
      if (this%Wm(iconst) == 0) cycle

      wm = this%Wm(iconst)

      use_neighbor = PETSC_FALSE
      imat_id = rblock(iconst,1)
      if (imat_id <= 0) cycle
      irb = rblock(iconst,3)
      select case(constrained_block%structure_metric(irb))
        case(1:2,7:10)
          imat_id_nb = rblock(iconst,2)
          if (imat_id_nb <=0 ) cycle
          use_neighbor = PETSC_TRUE
        case(5:6)
          option%io_buffer = 'Supported STRUCTURE_METRIC in INVERSION, &
                             &CONSTRAINED_BLOCKS is between 1-4 and 7-10 for &
                             &parameter based approach.'
          call PrintErrMsg(option)
        case default
      end select

      iparameter = this%inversion_aux%parameters(1)%iparameter

      call InvAuxGetSetParamValueByMat(this%inversion_aux,param_ce, &
                                       iparameter, &
                                       imat_id,INVAUX_GET_MATERIAL_VALUE)
      if (use_neighbor) then
        call InvAuxGetSetParamValueByMat(this%inversion_aux,param_nb, &
                                         iparameter,imat_id_nb, &
                                         INVAUX_GET_MATERIAL_VALUE)
      endif

      x = 0.0d0

      select case(constrained_block%structure_metric(irb))
        case(1)
          x = log(param_ce) - log(param_nb)
        case(2)
          x = log(param_ce) - log(param_nb)
        case(3)
          x = log(param_ce) - log(constrained_block%reference_parameter(irb))
        case(4)
          x = log(param_ce) - log(constrained_block%reference_parameter(irb))
        case(5)
          x = log(param_ce) - log(param_nb)
          ! TODO: compute rx,ry, and rz
        case(6)
          x = log(param_ce) - log(param_nb)
          ! TODO: compute rx,ry, and rz
        case(7)
          x = (log(param_ce) - log(constrained_block%reference_parameter(irb)))&
             -(log(param_nb) - log(constrained_block%reference_parameter(irb)))
        case(8)
          x = (log(param_ce) - log(constrained_block%reference_parameter(irb)))&
             -(log(param_nb) - log(constrained_block%reference_parameter(irb)))
        case(9)
          x = log(param_ce) - log(param_nb)
        case(10)
          x = log(param_ce) - log(param_nb)
        case default
          option%io_buffer = 'Supported STRUCTURE_METRIC in INVERSION, &
                              &CONSTRAINED_BLOCKS is between 1 to 10'
          call PrintErrMsg(option)
      end select

      this%b(num_measurement + iconst) = - sqrt(beta) * wm * x

    enddo

  endif

end subroutine InversionZFlowCGLSRhs

! ************************************************************************** !

subroutine InversionZFlowBuildWm(this)
  !
  ! Builds model regularization matrix: Wm
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/23/21
  !

  use Patch_module
  use Grid_module
  use Material_Aux_module
  use Option_module
  use String_module
  use Variables_module, only : PERMEABILITY,POROSITY,ELECTRICAL_CONDUCTIVITY

  implicit none

  class(inversion_zflow_type) :: this

  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: i

  option => this%realization%option
  patch => this%realization%patch
  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars

  constrained_block => this%constrained_block

  if ((.not.associated(this%Wm)).or.(.not.associated(this%rblock))) &
                                          call InversionZFlowAllocateWm(this)

  do i=1,this%num_constraints_local
    call ComputeWm(i,this%Wm(i))
  enddo

contains

  subroutine ComputeWm(iconst,wm)
    ! computes an element of Wm matrix
    !
    ! Author: Piyoosh Jaysaval
    ! Date: 06/18/21

    implicit none

    PetscInt :: iconst
    PetscReal :: wm

    PetscInt :: irb,ghosted_id,ghosted_id_nb
    PetscInt :: imat_id,imat_id_nb
    PetscInt :: iparameter
    PetscReal :: x,awx,awy,awz
    PetscReal :: param_ce,param_nb     ! cell's and neighbor's
    PetscReal :: mn,sd
    PetscReal :: rx,ry,rz,r
    PetscInt, pointer :: rblock(:,:)
    PetscBool :: use_neighbor

    rblock => this%rblock

    if (this%inversion_aux%qoi_is_full_vector) then
      ! get param & block of the ith constrained eq.
      use_neighbor = PETSC_FALSE
      ghosted_id = rblock(iconst,1)
      if (patch%imat(ghosted_id) <= 0) return
      irb = rblock(iconst,3)
      select case(constrained_block%structure_metric(irb))
        case(1:2,5:10)
          ghosted_id_nb = rblock(iconst,2)
          if (patch%imat(ghosted_id_nb) <=0 ) return
          use_neighbor = PETSC_TRUE
        case default
      end select

      iparameter = this%inversion_aux%parameters(1)%iparameter

      call InvAuxGetParamValueByCell(this%inversion_aux,param_ce, &
                                     iparameter, &
                                     patch%imat(ghosted_id), &
                                     material_auxvars(ghosted_id))
      if (use_neighbor) then
        call InvAuxGetParamValueByCell(this%inversion_aux,param_nb, &
                                       iparameter, &
                                       patch%imat(ghosted_id_nb), &
                                       material_auxvars(ghosted_id_nb))
      endif

    else
      ! get material param & block of the ith constrained eq.
      use_neighbor = PETSC_FALSE
      imat_id = rblock(iconst,1)
      if (imat_id <= 0) return
      irb = rblock(iconst,3)
      select case(constrained_block%structure_metric(irb))
        case(1:2,7:10)
          imat_id_nb = rblock(iconst,2)
          if (imat_id_nb <= 0) return
          use_neighbor = PETSC_TRUE
        case(5:6)
          option%io_buffer = 'Supported STRUCTURE_METRIC in INVERSION, &
                             &CONSTRAINED_BLOCKS is between 1-4 and 7-10 for &
                             &parameter based approach.'
        call PrintErrMsg(option)
        case default
      end select

      iparameter = this%inversion_aux%parameters(1)%iparameter

      call InvAuxGetSetParamValueByMat(this%inversion_aux,param_ce, &
                                       iparameter, &
                                       imat_id,INVAUX_GET_MATERIAL_VALUE)
      if (use_neighbor) then
        call InvAuxGetSetParamValueByMat(this%inversion_aux,param_nb, &
                                         iparameter, &
                                         imat_id_nb, &
                                         INVAUX_GET_MATERIAL_VALUE)
      endif

    endif

    x = 0.d0

    select case(constrained_block%structure_metric(irb))
      case(1)
        x = log(param_ce) - log(param_nb)
      case(2)
        x = abs(log(param_ce) - log(param_nb))
      case(3)
        x = log(param_ce) - log(constrained_block%reference_parameter(irb))
      case(4)
        x = abs(log(param_ce) - log(constrained_block%reference_parameter(irb)))
      case(5)
        !x = log(param_ce) - log(param_nb)

        ! compute unit vectors: rx,ry, and rz
        rx = grid%x(ghosted_id) - grid%x(ghosted_id_nb)
        ry = grid%y(ghosted_id) - grid%y(ghosted_id_nb)
        rz = grid%z(ghosted_id) - grid%z(ghosted_id_nb)
        r = sqrt(rx*rx + ry*ry + rz*rz)
        rx = rx / r
        ry = ry / r
        rz = rz / r
      case(6)
        !x = abs(log(param_ce) - log(param_nb))

        ! compute unit vectors: rx,ry, and rz
        rx = abs(grid%x(ghosted_id) - grid%x(ghosted_id_nb))
        ry = abs(grid%y(ghosted_id) - grid%y(ghosted_id_nb))
        rz = abs(grid%z(ghosted_id) - grid%z(ghosted_id_nb))
        r = sqrt(rx*rx + ry*ry + rz*rz)
        rx = rx / r
        ry = ry / r
        rz = rz / r
      case(7)
        x = (log(param_ce) - log(constrained_block%reference_parameter(irb))) &
           -(log(param_nb) - log(constrained_block%reference_parameter(irb)))
      case(8)
        x = abs( &
            (log(param_ce) - log(constrained_block%reference_parameter(irb))) &
           -(log(param_nb) - log(constrained_block%reference_parameter(irb))))
      case(9)
        x = log(param_ce) - log(param_nb)
      case(10)
        x = abs(log(param_ce) - log(param_nb))
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
      if ((x-mn) < 0) then
        wm = 1 / (sd*sd)
      else
        wm = sd*sd / (((x-mn)*(x-mn) + sd*sd)*((x-mn)*(x-mn) + sd*sd))
      end if
    case(6)
      if ((x-mn) > 0) then
        wm = 1 / (sd*sd)
      else
        wm = sd*sd / (((x-mn)*(x-mn) + sd*sd)*((x-mn)*(x-mn) + sd*sd))
      end if
    case default
      option%io_buffer = 'Supported WEIGHING_FUNCTION in INVERSION, &
                          &CONSTRAINED_BLOCKS is between 1 to 6'
      call PrintErrMsg(option)
    end select

    if (constrained_block%structure_metric(irb) == 5 .or. &
        constrained_block%structure_metric(irb) == 6) then
      awx = constrained_block%aniso_weight(irb,1)
      awy = constrained_block%aniso_weight(irb,2)
      awz = constrained_block%aniso_weight(irb,3)
      wm = (1 - abs( awx*rx + awy*ry + awz*rz))**2
    end if

    wm = constrained_block%relative_weight(irb) * wm

  end subroutine ComputeWm

end subroutine InversionZFlowBuildWm

! ************************************************************************** !

subroutine InversionZFlowAllocateWm(this)
  !
  ! Allocate and get info on Wm and rblock
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/18/21
  !

  use Patch_module
  use Grid_module
  use Option_module

  implicit none

  class(inversion_zflow_type) :: this

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: local_id,ghosted_id,ghosted_id_nbr
  PetscInt :: imat_id,imat_id_nb
  PetscInt :: iconblock,inbr,ilink
  PetscInt :: num_constraints
  PetscInt :: num_neighbor
  PetscErrorCode :: ierr

  option => this%realization%option
  patch => this%realization%patch
  grid => patch%grid

  constrained_block => this%constrained_block

  if (this%inversion_aux%qoi_is_full_vector) then
    num_constraints = 0
    do local_id=1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
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
                if (patch%imat(ghosted_id_nbr) <= 0) cycle
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

    this%num_constraints_local = num_constraints
    call MPI_Allreduce(num_constraints,this%num_constraints_total, &
                       ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm, &
                       ierr);CHKERRQ(ierr)
    allocate(this%Wm(num_constraints))
    allocate(this%rblock(num_constraints,THREE_INTEGER))
    this%Wm = 0.d0
    this%rblock = 0

    ! repeat once num_constraints is known
    num_constraints = 0
    do local_id=1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      do iconblock=1,constrained_block%num_constrained_block
        if (constrained_block%structure_metric(iconblock) > 0) then
          if (constrained_block%material_id(iconblock) == &
              patch%imat(ghosted_id)) then
            if (constrained_block%structure_metric(iconblock) == 3 .or. &
                constrained_block%structure_metric(iconblock) == 4) then
              num_constraints = num_constraints + 1
              this%rblock(num_constraints,1) = ghosted_id
              this%rblock(num_constraints,3) = iconblock
            else
              num_neighbor = grid%cell_neighbors_local_ghosted(0,local_id)
              do inbr=1,num_neighbor
                ghosted_id_nbr = abs( &
                              grid%cell_neighbors_local_ghosted(inbr,local_id))
                if (patch%imat(ghosted_id_nbr) <= 0) cycle
                if (patch%imat(ghosted_id_nbr) /= patch%imat(ghosted_id)) then
                  do ilink=1,constrained_block%block_link(iconblock,1)
                    if (constrained_block%block_link(iconblock,ilink+1) == &
                        patch%imat(ghosted_id_nbr)) then
                      num_constraints = num_constraints + 1
                      this%rblock(num_constraints,1) = ghosted_id
                      this%rblock(num_constraints,2) = ghosted_id_nbr
                      this%rblock(num_constraints,3) = iconblock
                    endif
                  enddo
                else
                  if (constrained_block%structure_metric(iconblock) < 9 .or. &
                      constrained_block%structure_metric(iconblock) > 10) then
                    num_constraints = num_constraints + 1
                    this%rblock(num_constraints,1) = ghosted_id
                    this%rblock(num_constraints,2) = ghosted_id_nbr
                    this%rblock(num_constraints,3) = iconblock
                  endif
                endif
              enddo
            endif
          endif
        endif
      enddo
    enddo
  else
    num_constraints = 0
    do iconblock=1,constrained_block%num_constrained_block
      imat_id = constrained_block%material_id(iconblock)
      if (constrained_block%structure_metric(iconblock) > 0 .and. &
          imat_id > 0) then
        if (constrained_block%structure_metric(iconblock) == 3 .or. &
            constrained_block%structure_metric(iconblock) == 4) then
          num_constraints = num_constraints + 1
        else
          do ilink=1,constrained_block%block_link(iconblock,1)
            num_constraints = num_constraints + 1
          enddo
        endif
      endif
    enddo

    this%num_constraints_local = num_constraints
    this%num_constraints_total = this%num_constraints_local ! both are same
    allocate(this%Wm(num_constraints))
    allocate(this%rblock(num_constraints,THREE_INTEGER))
    this%Wm = 0.d0
    this%rblock = 0

    ! repeat once num_constraints is known
    num_constraints = 0
    do iconblock=1,constrained_block%num_constrained_block
      imat_id = constrained_block%material_id(iconblock)
      if (constrained_block%structure_metric(iconblock) > 0 .and. &
          imat_id > 0) then
        if (constrained_block%structure_metric(iconblock) == 3 .or. &
            constrained_block%structure_metric(iconblock) == 4) then
          num_constraints = num_constraints + 1
          this%rblock(num_constraints,1) = imat_id
          this%rblock(num_constraints,3) = iconblock
        else
          do ilink=1,constrained_block%block_link(iconblock,1)
            imat_id_nb = constrained_block%block_link(iconblock,1+ilink)
            num_constraints = num_constraints + 1
            this%rblock(num_constraints,1) = imat_id
            this%rblock(num_constraints,2) = imat_id_nb
            this%rblock(num_constraints,3) = iconblock
          enddo
        endif
      endif
    enddo

  endif

end subroutine InversionZFlowAllocateWm

! ************************************************************************** !

subroutine InversionZFlowComputeMatVecProductJp(this)
  !
  ! Computes product of Jacobian J with a vector p = Jp
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/14/21
  !

  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module
  use Option_module

  implicit none

  class(inversion_zflow_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(constrained_block_type), pointer :: constrained_block
  type(inversion_aux_type), pointer :: inversion_aux

  PetscInt :: iconst,irb,num_measurement
  PetscInt :: ghosted_id,ghosted_id_nb
  PetscInt :: imat_id,imat_id_nb
  PetscInt :: iblock_id,iblock_id_nb
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: beta,wm
  PetscReal, pointer :: pvec_ptr(:)
  PetscReal, pointer :: q1vec_ptr(:)
  PetscErrorCode :: ierr

  Vec :: p1
  Vec :: q1
  Vec :: q1_dist

  option => this%realization%option
  field => this%realization%field
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  inversion_aux => this%inversion_aux

  constrained_block => this%constrained_block
  rblock => this%rblock

  this%q = 0.d0

  num_measurement = size(this%inversion_aux%measurements)

  ! Data part
  call VecDuplicate(this%inversion_aux%dist_parameter_vec,p1,ierr);CHKERRQ(ierr)
  call VecDuplicate(this%inversion_aux%dist_measurement_vec,q1_dist, &
                    ierr);CHKERRQ(ierr)
  call VecDuplicate(this%inversion_aux%measurement_vec,q1,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(p1,pvec_ptr,ierr);CHKERRQ(ierr)
  pvec_ptr = this%p
  call VecRestoreArrayF90(p1,pvec_ptr,ierr);CHKERRQ(ierr)

  ! q = Jp -> data part
  call MatMultTranspose(inversion_aux%JsensitivityT,p1,q1_dist, &
                        ierr);CHKERRQ(ierr)

  call InvAuxScatMeasToDistMeas(this%inversion_aux, &
                                q1, &
                                q1_dist, &
                                INVAUX_SCATREVERSE)

  call VecGetArrayF90(q1,q1vec_ptr,ierr);CHKERRQ(ierr)
  this%q(1:num_measurement) = q1vec_ptr
  call VecRestoreArrayF90(q1,q1vec_ptr,ierr);CHKERRQ(ierr)

  ! Model part -> q2
  beta = this%beta

  if (this%inversion_aux%qoi_is_full_vector) then
    ! Get local this%p to ghosted in pvec_ptr
    call InvAuxScatGlobalToDistParam(this%inversion_aux, &
                                     this%realization%field%work, &
                                     p1, &
                                     INVAUX_SCATREVERSE)
    call DiscretizationGlobalToLocal(discretization,field%work, &
                                    field%work_loc,ONEDOF)
    call VecGetArrayF90(field%work_loc,pvec_ptr,ierr);CHKERRQ(ierr)

    do iconst=1,this%num_constraints_local
      if (this%Wm(iconst) == 0) cycle

      wm = this%Wm(iconst)
      irb = rblock(iconst,3)
      ghosted_id = rblock(iconst,1)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (constrained_block%structure_metric(irb) == 3 .or. &
          constrained_block%structure_metric(irb) == 4) then
        this%q(num_measurement + iconst) = &
          sqrt(beta) * wm * pvec_ptr(ghosted_id)
      else
        ghosted_id_nb = rblock(iconst,2)
        if (patch%imat(ghosted_id_nb) <= 0) cycle
        this%q(num_measurement + iconst) = &
            sqrt(beta) * wm * (pvec_ptr(ghosted_id) - pvec_ptr(ghosted_id_nb))
      endif
    enddo

    call VecRestoreArrayF90(field%work_loc,pvec_ptr,ierr);CHKERRQ(ierr)

  else

    call VecZeroEntries(this%parameter_tmp_vec,ierr);CHKERRQ(ierr)
    call InvAuxScatParamToDistParam(this%inversion_aux, &
                                    this%parameter_tmp_vec, &
                                    p1, &
                                    INVAUX_SCATREVERSE)

    call VecGetArrayF90(this%parameter_tmp_vec,pvec_ptr,ierr);CHKERRQ(ierr)

    do iconst=1,this%num_constraints_local
      if (this%Wm(iconst) == 0) cycle

      wm = this%Wm(iconst)
      irb = rblock(iconst,3)
      imat_id = rblock(iconst,1)
      if (imat_id <= 0) cycle
      iblock_id = ConstrainedBlockGetBlockIDFromMatID(constrained_block, &
                                            imat_id,this%realization%option)
      if (constrained_block%structure_metric(irb) == 3 .or. &
          constrained_block%structure_metric(irb) == 4) then
        this%q(num_measurement + iconst) = &
          sqrt(beta) * wm * pvec_ptr(iblock_id)
      else
        imat_id_nb = rblock(iconst,2)
        if (imat_id_nb <= 0) cycle
        iblock_id_nb = ConstrainedBlockGetBlockIDFromMatID(constrained_block, &
                                         imat_id_nb,this%realization%option)
        this%q(num_measurement + iconst) = &
            sqrt(beta) * wm * (pvec_ptr(iblock_id) - pvec_ptr(iblock_id_nb))
      endif
    enddo

    call VecRestoreArrayF90(this%parameter_tmp_vec,pvec_ptr,ierr);CHKERRQ(ierr)

  endif

  call VecDestroy(p1,ierr);CHKERRQ(ierr)
  call VecDestroy(q1,ierr);CHKERRQ(ierr)
  call VecDestroy(q1_dist,ierr);CHKERRQ(ierr)

end subroutine InversionZFlowComputeMatVecProductJp

!************************************************************************** !

subroutine InversionZFlowComputeMatVecProductJtr(this)
  !
  ! Computes product of Jacobian J transpose with a vector r = J^t x r
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/14/21
  !

  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module

  implicit none

  class(inversion_zflow_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(constrained_block_type), pointer :: constrained_block
  type(inversion_aux_type), pointer :: inversion_aux

  PetscInt :: iconst,irb,num_measurement
  PetscInt :: ghosted_id,ghosted_id_nb
  PetscInt :: imat_id,imat_id_nb
  PetscInt :: iblock_id,iblock_id_nb
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: beta,wm
  PetscReal, pointer :: r1vec_ptr(:)
  PetscReal, pointer :: s1vec_ptr(:)
  PetscReal, pointer :: s2vec_ptr(:)
  PetscErrorCode :: ierr

  Vec :: r1
  Vec :: s1

  field => this%realization%field
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  inversion_aux => this%inversion_aux

  constrained_block => this%constrained_block
  rblock => this%rblock

  this%s = 0.0d0

  num_measurement = size(this%inversion_aux%measurements)

  call VecZeroEntries(this%dist_parameter_tmp_vec,ierr);CHKERRQ(ierr)

  ! Model part -> s2

  beta = this%beta

  if (this%inversion_aux%qoi_is_full_vector) then

    call VecGetArrayF90(field%work_loc,s2vec_ptr,ierr);CHKERRQ(ierr)
    s2vec_ptr = 0.d0

    do iconst=1,this%num_constraints_local
      if (this%Wm(iconst) == 0) cycle

      wm = this%Wm(iconst)
      irb = rblock(iconst,3)
      ghosted_id = rblock(iconst,1)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (constrained_block%structure_metric(irb) == 3 .or. &
          constrained_block%structure_metric(irb) == 4) then
        s2vec_ptr(ghosted_id) = s2vec_ptr(ghosted_id) + &
          sqrt(beta) * wm * this%r(num_measurement + iconst)
      else
        ghosted_id_nb = rblock(iconst,2)
        if (patch%imat(ghosted_id_nb) <= 0) cycle
        s2vec_ptr(ghosted_id) = s2vec_ptr(ghosted_id) + &
                sqrt(beta) * wm * this%r(num_measurement + iconst)
        s2vec_ptr(ghosted_id_nb) = s2vec_ptr(ghosted_id_nb) - &
                sqrt(beta) * wm * this%r(num_measurement + iconst)
      endif
    enddo

    call VecRestoreArrayF90(field%work_loc,s2vec_ptr,ierr);CHKERRQ(ierr)

    ! s2 in field%work
    call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
    call DiscretizationLocalToGlobalAdd(discretization,field%work_loc, &
                                        field%work,ONEDOF)
    call InvAuxScatGlobalToDistParam(this%inversion_aux, &
                                     this%realization%field%work, &
                                     this%dist_parameter_tmp_vec, &
                                     INVAUX_SCATFORWARD)

  else

    call VecZeroEntries(this%parameter_tmp_vec,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%parameter_tmp_vec,s2vec_ptr,ierr);CHKERRQ(ierr)
    s2vec_ptr = 0.d0

    do iconst=1,this%num_constraints_local
      if (this%Wm(iconst) == 0) cycle

      wm = this%Wm(iconst)
      irb = rblock(iconst,3)
      imat_id = rblock(iconst,1)
      if (imat_id <= 0) cycle
      iblock_id = ConstrainedBlockGetBlockIDFromMatID(constrained_block, &
                                            imat_id,this%realization%option)
      if (constrained_block%structure_metric(irb) == 3 .or. &
          constrained_block%structure_metric(irb) == 4) then
        s2vec_ptr(iblock_id) = s2vec_ptr(iblock_id) + &
                        sqrt(beta) * wm * this%r(num_measurement + iconst)
      else
        imat_id_nb = rblock(iconst,2)
        if (imat_id_nb <= 0) cycle
        iblock_id_nb = ConstrainedBlockGetBlockIDFromMatID(constrained_block, &
                                         imat_id_nb,this%realization%option)
        s2vec_ptr(iblock_id) = s2vec_ptr(iblock_id) + &
                        sqrt(beta) * wm * this%r(num_measurement + iconst)
        s2vec_ptr(iblock_id_nb) = s2vec_ptr(iblock_id_nb) - &
                        sqrt(beta) * wm * this%r(num_measurement + iconst)
      endif
    enddo

    call VecRestoreArrayF90(this%parameter_tmp_vec,s2vec_ptr,ierr);CHKERRQ(ierr)

    call InvAuxScatParamToDistParam(this%inversion_aux, &
                                    this%parameter_tmp_vec, &
                                    this%dist_parameter_tmp_vec, &
                                    INVAUX_SCATFORWARD)
  endif

  ! Data part
  call VecDuplicate(this%inversion_aux%measurement_vec,r1,ierr);CHKERRQ(ierr)
  call VecDuplicate(this%dist_parameter_tmp_vec,s1,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(r1,r1vec_ptr,ierr);CHKERRQ(ierr)
  r1vec_ptr = this%r(1:num_measurement)
  call VecRestoreArrayF90(r1,r1vec_ptr,ierr);CHKERRQ(ierr)
  call InvAuxScatMeasToDistMeas(this%inversion_aux, &
                                r1, &
                                this%inversion_aux%dist_measurement_vec, &
                                INVAUX_SCATFORWARD)
  call VecGetArrayF90(this%inversion_aux%dist_measurement_vec,r1vec_ptr, &
                      ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(this%inversion_aux%dist_measurement_vec,r1vec_ptr, &
                          ierr);CHKERRQ(ierr)

  ! s = J^T*r -> data part
  call MatMult(inversion_aux%JsensitivityT, &
               this%inversion_aux%dist_measurement_vec,s1,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(s1,s1vec_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(this%dist_parameter_tmp_vec,s2vec_ptr, &
                      ierr);CHKERRQ(ierr)
  this%s = s1vec_ptr + s2vec_ptr
  call VecRestoreArrayF90(s1,s1vec_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(this%dist_parameter_tmp_vec,s2vec_ptr, &
                          ierr);CHKERRQ(ierr)

  call VecDestroy(r1,ierr);CHKERRQ(ierr)
  call VecDestroy(s1,ierr);CHKERRQ(ierr)

end subroutine InversionZFlowComputeMatVecProductJtr

! ************************************************************************** !

subroutine InvZFlowWriteIterationInfo(this)

  class(inversion_zflow_type) :: this

  if (this%info_format == 1) then
    call InvZFlowWriteIterationInfo2(this)
  else
    call InvZFlowWriteIterationInfo1(this)
  endif

end subroutine InvZFlowWriteIterationInfo

! ************************************************************************** !

subroutine InvZFlowWriteIterationInfo1(this)
  !
  ! Writes inversion run info
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/20/21
  !

  use String_module

  implicit none

  class(inversion_zflow_type) :: this

  PetscInt :: fid
  PetscInt, parameter :: zeronum = 0

  call InvSubsurfWriteIterationInfo(this)
  if (this%driver%PrintToScreen()) then
    write(*,*)
    write(*,98)
    if (this%string_color) then
      if (this%iteration == this%start_iteration) then
        write(*,'(/,2x,a,i6.4,/)') StringColor("CONVERGENCE STATISTICS AT &
                                &STARTING ITERATION:",C_RED), &
                                this%start_iteration
      else
        write(*,'(/,2x,a,i6.4,/)') StringColor("CONVERGENCE STATISTICS AFTER &
                                  &ITERATION:",C_RED),this%iteration
      endif
    else
      if (this%iteration == this%start_iteration) then
        write(*,'(/,2x,a,i6.4,/)') "CONVERGENCE STATISTICS AT &
                                &STARTING ITERATION:", &
                                this%start_iteration
      else
        write(*,'(/,2x,a,i6.4,/)') "CONVERGENCE STATISTICS AFTER &
                                  &ITERATION:", this%iteration
      endif
    endif
    write(*,99)
    if (this%string_color) then
      write(*,*) StringColor("  Phi_data   ",C_GREEN), &
                StringColor("   Phi_Model  ",C_BLUE), &
                StringColor("  Phi_Model/Beta",C_MAGENTA), &
                StringColor("    Phi_Total   ",C_CYAN)
    else
      write(*,*) "  Phi_data   ", &
                 "   Phi_Model  ", &
                 "  Phi_Model/Beta", &
                 "    Phi_Total   "
    endif
    write(*,102) this%phi_data,this%phi_model,this%phi_model/this%beta, &
                 this%phi_total
    write(*,*)
    if (this%num_constraints_total >= 0) then
      write(*,103) this%num_constraints_total
    else
      write(*,103) zeronum
    endif
    write(*,104) this%current_chi2
    write(*,105) this%target_chi2
    write(*,106) sqrt(this%current_chi2)
    write(*,107) this%beta
    write(*,108) this%beta_red_factor
    write(*,109) 100.d0*(this%phi_total_0 - this%phi_total)/this%phi_total_0
    write(*,110) 100.d0*this%min_phi_red
    write(*,99)
    write(*,*)
    write(*,98)
  endif

  if (this%driver%PrintToFile()) then
    fid = this%driver%fid_out
    write(fid,*)
    write(fid,98)
    if (this%iteration == this%start_iteration) then
      write(fid,'(/,2x,a,i6.4,/)') "CONVERGENCE STATISTICS AT STARTING &
                                   &ITERATION:", this%start_iteration
    else
      write(fid,'(/,2x,a,i6.4,/)') "CONVERGENCE STATISTICS AFTER ITERATION:", &
                                    this%iteration
    endif
    write(fid,99)
    write(fid,101) "  Phi_data   ","   Phi_Model   "," Phi_Model/Beta", &
                  &"   Phi_Total   "
    write(fid,102) this%phi_data,this%phi_model,this%phi_model/this%beta, &
                   this%phi_total
    write(fid,*)
    if (this%num_constraints_total >= 0) then
      write(fid,103) this%num_constraints_total
    else
      write(fid,103) zeronum
    endif
    write(fid,104) this%current_chi2
    write(fid,105) this%target_chi2
    write(fid,106) sqrt(this%current_chi2)
    write(fid,107) this%beta
    write(fid,108) this%beta_red_factor
    write(fid,109) 100.d0*(this%phi_total_0 - this%phi_total)/this%phi_total_0
    write(fid,110) 100.d0*this%min_phi_red
    write(fid,99)
    write(fid,*)
    write(fid,98)
    flush(fid)
  endif

98 format(40('=+'))
99 format(80('~'))
101 format(4a15)
102 format(4g15.7)

103 format(4x,'Number of Constraint Eqs:      ',2x,i15.10)
104 format(4x,'Current Chi2:                  ',2x,f15.4)
105 format(4x,'Target Ch2:                    ',2x,f15.4)
106 format(4x,'RMS error:                     ',2x,f15.4)
107 format(4x,'Beta:                          ',2x,f15.4)
108 format(4x,'Beta reduction factor:         ',2x,f15.4)
109 format(4x,'Reduction in Phi_Total:        ',2x,f15.4," %")
110 format(4x,'Minimum reduction in Phi_Total ' /,8x, &
                 &'before Beta reduction:     ',2x,f15.4," %")

end subroutine InvZFlowWriteIterationInfo1

! ************************************************************************** !

subroutine InvZFlowWriteIterationInfo2(this)
  !
  ! Writes inversion run info
  !
  ! Author: Glenn Hammond
  ! Date: 12/16/22
  !

  use String_module

  implicit none

  class(inversion_zflow_type) :: this

  character(len=:), allocatable :: string
  character(len=:), allocatable :: nl
  character(len=80) :: divider

  nl = new_line('a')
  write(divider,'(40("=+"))')
  string = nl // trim(divider)
  call this%driver%PrintMsg(string)
  string = nl // ' Iteration ' // &
           StringWrite(this%iteration) // nl

  call InvSubsurfWriteIterationInfo(this)
  string = nl // ' Convergence statistics' // nl
  call this%driver%PrintMsg(string)
  string = Helper1('Phi_Data') // &
             Helper2(StringWriteF(this%phi_data)) // nl // &
           Helper1('Phi_Model') // &
             Helper2(StringWriteF(this%phi_model)) // nl // &
           Helper1('Phi_Model/Beta') // &
             Helper2(StringWriteF(this%phi_model/this%beta)) // nl // &
           Helper1('Phi_Total') // &
             Helper2(StringWriteF(this%phi_total)) // nl
  call this%driver%PrintMsg(string)
  string = Helper1('Number of Constraint Eqs') // &
             Helper2(StringWrite(this%num_constraints_total)) // nl // &
           Helper1('Current Chi2') // &
             Helper2(StringWriteF(this%current_chi2)) // nl // &
           Helper1('Target Chi2') // &
             Helper2(StringWriteF(this%target_chi2)) // nl // &
           Helper1('RMS error') // &
             Helper2(StringWriteF(sqrt(this%current_chi2))) // nl // &
           Helper1('Beta') // &
             Helper2(StringWriteF(this%beta)) // nl // &
           Helper1('Beta reduction factor') // &
             Helper2(StringWriteF(this%beta_red_factor)) // nl // &
           Helper1('Reduction in Phi_Total') // &
             Helper2(StringWriteF(100.d0*(this%phi_total_0 - &
                                 this%phi_total)/this%phi_total_0)) // &
             ' %' // nl // &
           '  Minimum reduction in Phi_Total' // nl // &
           Helper1('before Beta reduction') // &
             Helper2(StringWriteF(100.d0*this%min_phi_red)) // ' %'

  call this%driver%PrintMsg(string)

  string = nl // divider // nl
  call this%driver%PrintMsg(string)

contains

function Helper1(str)

  use String_module

  character(len=*) :: str

  character(len=35) :: Helper1

  Helper1 = trim(str) // ' :'
  Helper1 = adjustr(Helper1)

end function Helper1

function Helper2(str)

  use String_module

  character(len=*) :: str

  character(len=20) :: Helper2

  Helper2 = trim(str)
  Helper2 = adjustr(Helper2)

end function Helper2

end subroutine InvZFlowWriteIterationInfo2

! ************************************************************************** !

subroutine InversionZFlowScaleSensitivity(this)
  !
  ! Scales sensitivity Jacobian for ln(K)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/11/22
  !
  use Discretization_module
  use Realization_Base_class

  class(inversion_zflow_type) :: this

  Vec :: wd_vec
  PetscInt :: idata,num_measurement
  PetscReal, pointer :: wdvec_ptr(:)
  PetscErrorCode :: ierr

  num_measurement = size(this%inversion_aux%measurements)
  call VecDuplicate(this%inversion_aux%measurement_vec,wd_vec, &
                    ierr);CHKERRQ(ierr)
  call VecZeroEntries(wd_vec,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(wd_vec,wdvec_ptr,ierr);CHKERRQ(ierr)
  do idata = 1, num_measurement
    wdvec_ptr(idata) = this%inversion_aux%measurements(idata)%weight
  enddo
  call VecRestoreArrayF90(wd_vec,wdvec_ptr,ierr);CHKERRQ(ierr)
  call InvAuxScatMeasToDistMeas(this%inversion_aux, &
                                wd_vec, &
                                this%inversion_aux%dist_measurement_vec, &
                                INVAUX_SCATFORWARD)

  ! Column Scale with wd
  call MatDiagonalScale(this%inversion_aux%JsensitivityT,PETSC_NULL_VEC, &
                        this%inversion_aux%dist_measurement_vec, &
                        ierr);CHKERRQ(ierr)
  ! Row scale with parameter
  call MatDiagonalScale(this%inversion_aux%JsensitivityT, &
                        this%inversion_aux%dist_parameter_vec,PETSC_NULL_VEC, &
                        ierr);CHKERRQ(ierr)
  call VecDestroy(wd_vec,ierr);CHKERRQ(ierr)

end subroutine InversionZFlowScaleSensitivity

! ************************************************************************** !

subroutine InversionZFlowCheckpoint(this)
  !
  ! Checkpoints the values of parameters and inversion settings for the
  ! current iterate
  !
  ! Author: Glenn Hammond
  ! Date: 12/09/22
  !
  use hdf5
  use Driver_class
  use HDF5_Aux_module
  use String_module

  class(inversion_zflow_type) :: this

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  character(len=MAXSTRINGLENGTH) :: string
  integer :: hdf5_err
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  if (len_trim(this%checkpoint_filename) == 0) return

  call this%driver%PrintMsg('Checkpointing inversion iteration ' // &
                            trim(StringWrite(this%iteration)) // '.')
  call HDF5FileOpen(this%checkpoint_filename,file_id,(this%iteration==1), &
                    this%driver)
  call HDF5AttributeWrite(file_id,H5T_NATIVE_INTEGER,'Last Iteration', &
                          this%iteration,this%driver)
  string = 'Iteration ' // trim(StringWrite(this%iteration))
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  call HDF5AttributeWrite(grp_id,H5T_NATIVE_DOUBLE,'Phi Total', &
                          this%phi_total,this%driver)
  call HDF5AttributeWrite(grp_id,H5T_NATIVE_DOUBLE,'Phi Data', &
                          this%phi_data,this%driver)
  call HDF5AttributeWrite(grp_id,H5T_NATIVE_DOUBLE,'Phi Model', &
                          this%phi_model,this%driver)
  call VecGetArrayReadF90(this%inversion_aux%parameter_vec,vec_ptr, &
                          ierr);CHKERRQ(ierr)
  call HDF5DatasetWrite(grp_id,'Parameter Values',vec_ptr,this%driver)
  call VecRestoreArrayReadF90(this%inversion_aux%parameter_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(this%inversion_aux%measurement_vec,vec_ptr, &
                          ierr);CHKERRQ(ierr)
  call HDF5DatasetWrite(grp_id,'Measurement Values',vec_ptr,this%driver)
  call VecRestoreArrayReadF90(this%inversion_aux%measurement_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
  call h5gclose_f(grp_id,hdf5_err)
  call HDF5FileClose(file_id)

end subroutine InversionZFlowCheckpoint

! ************************************************************************** !

subroutine InversionZFlowRestartReadData(this)
  !
  ! Reads inversion parameters for a specific iteration from the inversion
  ! checkpoint file
  !
  ! Author: Glenn Hammond
  ! Date: 12/09/22
  !
  use hdf5
  use Driver_class
  use HDF5_Aux_module
  use String_module

  class(inversion_zflow_type) :: this

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  character(len=MAXSTRINGLENGTH) :: string
  integer :: hdf5_err
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  if (.not.this%inversion_aux%startup_phase .or. &
      Uninitialized(this%restart_iteration)) return

  call this%driver%PrintMsg('Reading inversion parameters from inversion &
                            &checkpoint file "' // &
                            trim(this%restart_filename) // '".')
  call HDF5FileOpen(this%restart_filename,file_id,PETSC_FALSE,this%driver)
  string = 'Iteration 1'
  call HDF5GroupOpen(file_id,string,grp_id,this%driver)
  call HDF5AttributeRead(grp_id,H5T_NATIVE_DOUBLE,'Phi Total', &
                         this%phi_total_0,this%driver)
  call HDF5AttributeRead(grp_id,H5T_NATIVE_DOUBLE,'Phi Data', &
                         this%phi_data_0,this%driver)
  call HDF5AttributeRead(grp_id,H5T_NATIVE_DOUBLE,'Phi Model', &
                         this%phi_model_0,this%driver)
  call h5gclose_f(grp_id,hdf5_err)
  string = 'Iteration ' // trim(StringWrite(this%restart_iteration))
  call HDF5GroupOpen(file_id,string,grp_id,this%driver)
  call HDF5AttributeRead(grp_id,H5T_NATIVE_DOUBLE,'Phi Total', &
                         this%phi_total,this%driver)
  call HDF5AttributeRead(grp_id,H5T_NATIVE_DOUBLE,'Phi Data', &
                         this%phi_data,this%driver)
  call HDF5AttributeRead(grp_id,H5T_NATIVE_DOUBLE,'Phi Model', &
                         this%phi_model,this%driver)
  call VecGetArrayReadF90(this%inversion_aux%parameter_vec,vec_ptr, &
                          ierr);CHKERRQ(ierr)
  call HDF5DatasetRead(grp_id,'Parameter Values',vec_ptr,this%driver)
  call VecRestoreArrayReadF90(this%inversion_aux%parameter_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
  call h5gclose_f(grp_id,hdf5_err)
  call HDF5FileClose(file_id)

  call InvAuxCopyParamToFromParamVec(this%inversion_aux, &
                                     INVAUX_PARAMETER_VALUE, &
                                     INVAUX_COPY_FROM_VEC)
  call InvAuxScatParamToDistParam(this%inversion_aux, &
                                  this%inversion_aux%parameter_vec, &
                                  this%inversion_aux%dist_parameter_vec, &
                                  INVAUX_SCATFORWARD)

end subroutine InversionZFlowRestartReadData

! ************************************************************************** !

subroutine InversionZFlowFinalize(this)
  !
  ! Finalizes inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/06/22
  !

  implicit none

  class(inversion_zflow_type) :: this

  call InversionBaseFinalize(this)

end subroutine InversionZFlowFinalize

! ************************************************************************** !

subroutine InversionZFlowStrip(this)
  !
  ! Deallocates members of inversion zflow
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/06/22
  !

  implicit none

  class(inversion_zflow_type) :: this

   call InversionSubsurfaceStrip(this)

end subroutine InversionZFlowStrip

! ************************************************************************** !

subroutine ConstrainedBlockParDestroy(constrained_block)
  !
  ! ConstrainedBlockParDestroy: Deallocates a constrained block par object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/06/22
  !

  use Utility_module, only : DeallocateArray

  implicit none

  type(constrained_block_par_type), pointer :: constrained_block

  if (.not.associated(constrained_block)) return

  call DeallocateArray(constrained_block%block_link)

  deallocate(constrained_block)
  nullify(constrained_block)

end subroutine ConstrainedBlockParDestroy

! ************************************************************************** !

subroutine ConstrainedBlockDestroy(constrained_block)
  !
  ! ConstrainedBlockParDestroy: Deallocates a constrained block par object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/06/22
  !

  use Utility_module, only : DeallocateArray

  implicit none

  type(constrained_block_type), pointer :: constrained_block

  type(constrained_block_par_type), pointer :: cur_constrained_block
  type(constrained_block_par_type), pointer :: prev_constrained_block

  if (.not.associated(constrained_block)) return

  cur_constrained_block => constrained_block%constrained_block_list
  do
    if (.not.associated(cur_constrained_block)) exit
    prev_constrained_block => cur_constrained_block
    cur_constrained_block => cur_constrained_block%next
    call ConstrainedBlockParDestroy(prev_constrained_block)
  enddo
  nullify(constrained_block%constrained_block_list)

  call DeallocateArray(constrained_block%material_name)
  call DeallocateArray(constrained_block%material_id)
  call DeallocateArray(constrained_block%structure_metric)
  call DeallocateArray(constrained_block%wf_type)
  call DeallocateArray(constrained_block%block_link)
  call DeallocateArray(constrained_block%wf_mean)
  call DeallocateArray(constrained_block%wf_sdev)
  call DeallocateArray(constrained_block%relative_weight)
  call DeallocateArray(constrained_block%aniso_weight)
  call DeallocateArray(constrained_block%reference_parameter)

  deallocate(constrained_block)
  nullify(constrained_block)

end subroutine ConstrainedBlockDestroy

! ************************************************************************** !

subroutine InversionZFlowDestroy(this)
  !
  ! Deallocates a inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !

  use Utility_module, only : DeallocateArray

  implicit none

  class(inversion_zflow_type), pointer :: this

  PetscErrorCode :: ierr

  if (.not.associated(this)) return

  call DeallocateArray(this%b)
  call DeallocateArray(this%p)
  call DeallocateArray(this%q)
  call DeallocateArray(this%r)
  call DeallocateArray(this%s)
  call DeallocateArray(this%del_param)
  call DeallocateArray(this%Wm)
  call DeallocateArray(this%rblock)

  if (this%dist_parameter_tmp_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%dist_parameter_tmp_vec,ierr);CHKERRQ(ierr)
  endif

  if (this%parameter_tmp_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%parameter_tmp_vec,ierr);CHKERRQ(ierr)
  endif

  call ConstrainedBlockDestroy(this%constrained_block)

  call this%Strip()
  deallocate(this)
  nullify(this)

end subroutine InversionZFlowDestroy

end module Inversion_ZFlow_class
