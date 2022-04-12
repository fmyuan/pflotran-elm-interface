module Inversion_ZFlow_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Inversion_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_subsurface_type) :: inversion_zflow_type
    PetscInt :: start_iteration          ! Starting iteration number
    PetscInt :: miniter,maxiter          ! min/max CGLS iterations

    PetscReal :: beta                    ! regularization parameter
    PetscReal :: beta_red_factor         ! beta reduction factor
    PetscReal :: minperm,maxperm         ! min/max permeability
    PetscReal :: target_chi2             ! target CHI^2 norm
    PetscReal :: current_chi2

    ! Cost/objective functions
    PetscReal :: min_phi_red             ! min change in cost function
    PetscReal :: phi_total_0,phi_total
    PetscReal :: phi_data_0,phi_data
    PetscReal :: phi_model_0,phi_model

    ! arrays for CGLS algorithm
    PetscReal, pointer :: b(:)           ! vector for CGLS RHS
    PetscReal, pointer :: p(:)           ! vector of dim -> num of inv cells
    PetscReal, pointer :: q(:)           ! product of Jacobian with p = Jp
    PetscReal, pointer :: r(:)           ! vector of dim -> num of measur
    PetscReal, pointer :: s(:)           ! product Jacobian transpose with r
    PetscReal, pointer :: del_perm(:)    ! permeability update vector

    ! For Wm
    PetscInt :: num_constraints_local    ! Number of constraints
    PetscInt :: num_constraints_total    ! Total number of constraints
    PetscInt, pointer :: rblock(:,:)     ! array stores info about reg.
    PetscReal, pointer :: Wm(:)          ! Regularization matrix
    Vec :: dist_parameter_tmp_vec

    type(constrained_block_type), pointer :: constrained_block

  contains
    procedure, public :: Init => InversionZFlowInit
    procedure, public :: ReadBlock => InversionZFlowReadBlock
    procedure, public :: Initialize => InversionZFlowInitialize
    procedure, public :: EvaluateCostFunction => InvZFlowEvaluateCostFunction
    procedure, public :: CheckConvergence => InversionZFlowCheckConvergence
    procedure, public :: WriteIterationInfo => InversionZFlowWriteIterationInfo
    procedure, public :: ScaleSensitivity => InversionZFlowScaleSensitivity
    procedure, public :: CalculateUpdate => InversionZFlowCalculateUpdate
    procedure, public :: UpdateParameters => InversionZFlowUpdateParameters
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
    PetscReal, pointer :: reference_permeability(:)
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
    PetscReal :: reference_permeability
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

  use Driver_module

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
  use Driver_module

  implicit none

  class(inversion_zflow_type) :: this
  class(driver_type), pointer :: driver

  call InversionSubsurfaceInit(this,driver)

  ! Default inversion parameters
  this%miniter = 10
  this%maxiter = 50

  this%beta = 100.d0
  this%beta_red_factor = 0.5d0
  this%minperm = 1d-17
  this%maxperm = 1d-07
  this%target_chi2 = 1.d0
  this%min_phi_red = 0.2d0

  this%start_iteration = 1
  this%maximum_iteration = 20
  this%num_constraints_local = UNINITIALIZED_INTEGER
  this%num_constraints_total = UNINITIALIZED_INTEGER
  this%current_chi2 = UNINITIALIZED_DOUBLE
  this%phi_total_0 = UNINITIALIZED_DOUBLE
  this%phi_data_0 = UNINITIALIZED_DOUBLE
  this%phi_model_0 = UNINITIALIZED_DOUBLE
  this%phi_total = UNINITIALIZED_DOUBLE
  this%phi_data = UNINITIALIZED_DOUBLE
  this%phi_model = UNINITIALIZED_DOUBLE

  this%dist_parameter_tmp_vec = PETSC_NULL_VEC

  nullify(this%b)
  nullify(this%p)
  nullify(this%q)
  nullify(this%r)
  nullify(this%s)
  nullify(this%del_perm)
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
  nullify(constrained_block%reference_permeability)

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
  constrained_block%reference_permeability = 0.d0
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

  type(grid_type), pointer :: grid

  PetscInt :: num_measurement
  PetscInt :: num_constraints

  grid => this%realization%patch%grid

  num_measurement = size(this%measurements)
  num_constraints = this%num_constraints_local

  allocate(this%b(num_measurement + num_constraints))
  allocate(this%p(this%num_parameters_local))
  allocate(this%q(num_measurement + num_constraints))
  allocate(this%r(num_measurement + num_constraints))
  allocate(this%s(this%num_parameters_local))
  allocate(this%del_perm(this%num_parameters_local))

  this%b = 0.d0
  this%p = 0.d0
  this%q = 0.d0
  this%r = 0.d0
  this%s = 0.d0
  this%del_perm = 0.d0

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
  call DeallocateArray(this%del_perm)

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
    allocate(constrained_block%reference_permeability(nconblock))
    constrained_block%reference_permeability = 0.d0
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
        option%io_buffer = 'Contrained block " &
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
      constrained_block%reference_permeability(iconblock) = &
                        cur_constrained_block%reference_permeability
      constrained_block%block_link(iconblock,ONE_INTEGER) = &
                        cur_constrained_block%num_block_link
      do i=1,cur_constrained_block%num_block_link
        material_property => &
            MaterialPropGetPtrFromArray(cur_constrained_block%block_link(i), &
                                      patch%material_property_array)
        if (.not.associated(material_property)) then
          option%io_buffer = 'Linked block "&
                             &'//trim(cur_constrained_block%block_link(i)) // &
                             &'" in contrained block "&
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
      case('MIN_PERMEABILITY')
        call InputReadDouble(input,option,this%minperm)
        call InputErrorMsg(input,option,'MIN_PERMEABILITY', &
                           error_string)
      case('MAX_PERMEABILITY')
        call InputReadDouble(input,option,this%maxperm)
        call InputErrorMsg(input,option,'MAX_PERMEABILITY', &
                           error_string)
      case('MIN_CGLS_ITERATION')
        call InputReadInt(input,option,this%miniter)
        call InputErrorMsg(input,option,'MIN_CGLS_ITERATION',error_string)
      case('MAX_CGLS_ITERATION')
        call InputReadInt(input,option,this%maxiter)
        call InputErrorMsg(input,option,'MAX_CGLS_ITERATION',error_string)
      case('BETA')
        call InputReadDouble(input,option,this%beta)
        call InputErrorMsg(input,option,'BETA',error_string)
      case('BETA_REDUCTION_FACTOR')
        call InputReadDouble(input,option,this%beta_red_factor)
        call InputErrorMsg(input,option,'BETA_REDUCTION_FACTOR',error_string)
      case('TARGET_CHI2')
        call InputReadDouble(input,option,this%target_chi2)
        call InputErrorMsg(input,option,'TARGET_CHI2',error_string)
      case('MIN_COST_REDUCTION')
        call InputReadDouble(input,option,this%min_phi_red)
        call InputErrorMsg(input,option,'MIN_COST_REDUCTION',error_string)
      case('CONSTRAINED_BLOCKS')
        call ConstrainedBlockRead(this%constrained_block,input,option)
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
      case('REFERENCE_PERMEABILITY')
        call InputReadDouble(input,option, &
                            constrained_block%reference_permeability)
        call InputErrorMsg(input,option,'REFERENCE_PERMEABILITY',error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine ConstrainedBlockParRead

! ************************************************************************** !

subroutine InversionZFlowInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/06/22
  !
  use Discretization_module
  use Inversion_TS_Aux_module
  use Inversion_Parameter_module
  use Option_module
  use Variables_module, only : PERMEABILITY

  implicit none

  class(inversion_zflow_type) :: this

  PetscBool :: exists
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: iqoi(2)
  PetscErrorCode :: ierr

  call InversionSubsurfInitialize(this)

  call VecDuplicate(this%dist_parameter_vec,this%dist_parameter_tmp_vec, &
                    ierr);CHKERRQ(ierr)

  ! check to ensure that quantity of interest exists
  exists = PETSC_FALSE
  iqoi = InversionParameterIntToQOIArray(this%parameters(1))
  select case(iqoi(1))
    case(PERMEABILITY)
      if (this%realization%option%iflowmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'PERMEABILITY'
    case default
  end select
  if (.not.exists) then
    this%realization%option%io_buffer = 'Inversion for ' // trim(word) // &
      &' cannot be performed with the specified process models.'
    call PrintErrMsg(this%realization%option)
  endif

  call InversionZFlowConstrainedArraysFromList(this)

  ! Build Wm matrix
  call InversionZFlowBuildWm(this)

end subroutine InversionZFlowInitialize

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

  implicit none

  class(inversion_zflow_type) :: this

  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: idata,num_measurement
  PetscInt :: iconst,num_constraints
  PetscInt :: irb,ghosted_id,ghosted_id_nb
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: wd,tempreal
  PetscReal, pointer :: vec_ptr(:)
  PetscReal :: perm_ce,perm_nb              ! cell's and neighbor's
  PetscReal :: wm,x
  PetscReal, allocatable :: model_vector(:)
  PetscErrorCode :: ierr

  option => this%realization%option
  patch => this%realization%patch
  material_auxvars => patch%aux%Material%auxvars

  constrained_block => this%constrained_block
  rblock => this%rblock

  num_measurement = size(this%measurements)

  ! Data part
  this%phi_data = 0.d0
  do idata=1,num_measurement

    wd = 0.05 * this%measurements(idata)%value
    wd = 1/wd

    tempreal = wd * (this%measurements(idata)%value - &
                     this%measurements(idata)%simulated_value)
    this%phi_data = this%phi_data + tempreal * tempreal

  enddo

  this%current_chi2 = this%phi_data / num_measurement

  ! model cost function
  this%phi_model = 0.d0
  num_constraints = this%num_constraints_local
  ! allocate to at least size 1 to allow for inner product
  allocate(model_vector(max(num_constraints,1)))
  model_vector = 0.d0

  do iconst=1,num_constraints
    if (this%Wm(iconst) == 0) cycle

    wm = this%Wm(iconst)

    ! get perm & block of the ith constrained eq.
    ghosted_id = rblock(iconst,1)
    ghosted_id_nb = rblock(iconst,2)
    if ((patch%imat(ghosted_id) <= 0) .or. &
        (patch%imat(ghosted_id_nb) <= 0)) cycle
    irb = rblock(iconst,3)
    perm_ce = material_auxvars(ghosted_id)%permeability(perm_xx_index)
    x = 0.d0

    select case(constrained_block%structure_metric(irb))
      case(1)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
      case(2)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = abs(log(perm_ce) - log(perm_nb))
      case(3)
        x = log(perm_ce) - log(constrained_block%reference_permeability(irb))
      case(4)
        x = abs(log(perm_ce) - &
                log(constrained_block%reference_permeability(irb)))
      case(5)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
      case(6)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = abs(log(perm_ce) - log(perm_nb))
      case(7)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = (log(perm_ce)- log(constrained_block%reference_permeability(irb))) &
          -(log(perm_nb) - log(constrained_block%reference_permeability(irb)))
      case(8)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = abs( &
            (log(perm_ce)- log(constrained_block%reference_permeability(irb))) &
          -(log(perm_nb) - log(constrained_block%reference_permeability(irb))) )
      case(9)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
      case(10)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = abs(log(perm_ce) - log(perm_nb))
      case default

    end select

    model_vector(iconst) = wm * x

  enddo

  this%phi_model = this%beta * dot_product(model_vector,model_vector)
  call MPI_Allreduce(MPI_IN_PLACE,this%phi_model,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
  deallocate(model_vector)

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

subroutine InversionZFlowUpdateParameters(this)
  !
  ! Updates input parameters
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/14/21
  !

  use Material_module
  use Discretization_module
  use Inversion_Parameter_module
  use Field_module

  class(inversion_zflow_type) :: this

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization

  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: iqoi(2)
  PetscErrorCode :: ierr

  field => this%realization%field
  discretization => this%realization%discretization

  iqoi = InversionParameterIntToQOIArray(this%parameters(1))
  call DiscretizationGlobalToLocal(discretization,this%dist_parameter_vec, &
                                   field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                               field%work_loc,iqoi(1),iqoi(2))

end subroutine InversionZFlowUpdateParameters

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
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid

  if (this%dist_parameter_vec /= PETSC_NULL_VEC) then

    call InversionZFlowAllocateWorkArrays(this)

    ! get inversion%del_perm
    call InversionZFlowCGLSSolve(this)

    call VecGetArrayF90(this%dist_parameter_tmp_vec,vec_ptr, &
                        ierr);CHKERRQ(ierr)
    vec_ptr(:) = this%del_perm(:)
    call VecRestoreArrayF90(this%dist_parameter_tmp_vec,vec_ptr, &
                            ierr);CHKERRQ(ierr)
    call InvSubsurfScatGlobalToDistParam(this, &
                                         this%realization%field%work, &
                                         this%dist_parameter_tmp_vec, &
                                         INVSUBSCATREVERSE)

    ! Get updated permeability as m_new = m_old + del_m (where m = log(perm))
    call VecGetArrayF90(this%dist_parameter_vec,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%realization%field%work,vec2_ptr,ierr);CHKERRQ(ierr)
    do iparameter = 1, this%num_parameters_local
      if (this%qoi_is_full_vector) then
        ghosted_id = grid%nL2G(iparameter)
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      vec_ptr(iparameter) = exp(log(vec_ptr(iparameter)) + vec2_ptr(iparameter))
      if (vec_ptr(iparameter) > this%maxperm) vec_ptr(iparameter) = this%maxperm
      if (vec_ptr(iparameter) < this%minperm) vec_ptr(iparameter) = this%minperm
    enddo
    call VecRestoreArrayF90(this%dist_parameter_vec,vec_ptr, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(this%realization%field%work,vec2_ptr, &
                            ierr);CHKERRQ(ierr)
    call InversionZFlowDeallocateWorkArrays(this)

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
  PetscBool :: lprint, l2print
  PetscErrorCode :: ierr

  PetscReal, parameter :: delta_initer = 1e-23
  PetscReal, parameter :: initer_conv  = 1e-24

  option => this%realization%option

  this%del_perm = 0.0d0

  timer => TimerCreate()
  call timer%Start()

  if (OptionPrintToScreen(option)) then
    write(*,'(" --> Solving ZFlow normal equation using CGLS solver:")')
  endif

  nm = size(this%measurements)
  ncons = this%num_constraints_local

  ! Get RHS vector this%b
  call InversionZFlowCGLSRhs(this)

  this%r = this%b

  ! get this%s = J^tr
  call InversionZFlowComputeMatVecProductJtr(this)
  this%p = this%s

  gamma = dot_product(this%s,this%s)
  call MPI_Allreduce(MPI_IN_PLACE,gamma,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

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
    delta2 = dot_product(this%q(nm+1:nm+ncons),this%q(nm+1:nm+ncons))
    call MPI_Allreduce(MPI_IN_PLACE,delta2,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
    delta = delta1 + delta2

    if (delta < 0) indefinite = PETSC_TRUE
    if (delta == 0) delta = epsilon(delta)

    alpha = gamma / delta

    this%del_perm = this%del_perm + alpha * this%p
    this%r = this%r - alpha * this%q

    ! get this%s = J^tr
    call InversionZFlowComputeMatVecProductJtr(this)

    gamma1 = gamma
    gamma = dot_product(this%s,this%s)
    call MPI_Allreduce(MPI_IN_PLACE,gamma,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

    norms = sqrt(gamma)
    gbeta = gamma / gamma1
    this%p = this%s + gbeta * this%p

    normx = dot_product(this%del_perm,this%del_perm)
    call MPI_Allreduce(MPI_IN_PLACE,normx,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
    normx = sqrt(normx)
    if (xmax < normx) xmax = normx
    if ( (norms <= norms0 * initer_conv) .or. (normx * initer_conv >= 1)) &
                               exit_info = PETSC_TRUE

    resNE_old = resNE
    resNE = norms / norms0

    if ( abs((resNE_old - resNe) /resNE_old) < delta_initer .and. &
        i > this%miniter) exit_info = PETSC_TRUE

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

  implicit none

  class(inversion_zflow_type) :: this

  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: idata,iconst,irb,num_measurement
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: perm_ce,perm_nb,x     ! cell's and neighbor's
  PetscReal :: wm,beta
  PetscReal :: wd
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  option => this%realization%option
  patch => this%realization%patch
  material_auxvars => patch%aux%Material%auxvars

  constrained_block => this%constrained_block
  rblock => this%rblock

  this%b = 0.0d0

  num_measurement = size(this%measurements)

  ! Data part
  do idata=1,num_measurement

    wd = 0.05 * this%measurements(idata)%value
    wd = 1/wd

    this%b(idata) = wd * (this%measurements(idata)%value - &
                          this%measurements(idata)%simulated_value)
  enddo

  ! Model part
  beta = this%beta

  do iconst=1,this%num_constraints_local
    if (this%Wm(iconst) == 0) cycle

    wm = this%Wm(iconst)

    perm_ce = material_auxvars(rblock(iconst,1))%permeability(perm_xx_index)
    irb = rblock(iconst,3)

    select case(constrained_block%structure_metric(irb))
      case(1)
        perm_nb = material_auxvars(rblock(iconst,2))%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
      case(2)
        perm_nb = material_auxvars(rblock(iconst,2))%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
      case(3)
        x = log(perm_ce) - log(constrained_block%reference_permeability(irb))
      case(4)
        x = log(perm_ce) - log(constrained_block%reference_permeability(irb))
      case(5)
        perm_nb = material_auxvars(rblock(iconst,2))%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
        ! TODO: compute rx,ry, and rz
      case(6)
        perm_nb = material_auxvars(rblock(iconst,2))%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
        ! TODO: compute rx,ry, and rz
      case(7)
        perm_nb = material_auxvars(rblock(iconst,2))%permeability(perm_xx_index)
        x = (log(perm_ce)- log(constrained_block%reference_permeability(irb))) &
          -(log(perm_nb) - log(constrained_block%reference_permeability(irb)))
      case(8)
        perm_nb = material_auxvars(rblock(iconst,2))%permeability(perm_xx_index)
        x = (log(perm_ce)- log(constrained_block%reference_permeability(irb))) &
          -(log(perm_nb) - log(constrained_block%reference_permeability(irb)))
      case(9)
        perm_nb = material_auxvars(rblock(iconst,2))%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
      case(10)
        perm_nb = material_auxvars(rblock(iconst,2))%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
      case default
        option%io_buffer = 'Supported STRUCTURE_METRIC in INVERSION, &
                            &CONSTRAINED_BLOCKS is between 1 to 10'
        call PrintErrMsg(option)
    end select

    this%b(num_measurement + iconst) = - sqrt(beta) * wm * x

  enddo

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
    PetscReal :: x,awx,awy,awz
    PetscReal :: perm_ce,perm_nb     ! cell's and neighbor's
    PetscReal :: mn,sd
    PetscReal :: rx,ry,rz,r
    PetscInt, pointer :: rblock(:,:)

    rblock => this%rblock

    ! get perm & block of the ith constrained eq.
    ghosted_id = rblock(iconst,1)
    ghosted_id_nb = rblock(iconst,2)
    if (patch%imat(ghosted_id) <= 0 .or.   &
        patch%imat(ghosted_id_nb) <=0 ) return
    irb = rblock(iconst,3)
    perm_ce = material_auxvars(ghosted_id)%permeability(perm_xx_index)
    x = 0.d0

    select case(constrained_block%structure_metric(irb))
      case(1)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
      case(2)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = abs(log(perm_ce) - log(perm_nb))
      case(3)
        x = log(perm_ce) - log(constrained_block%reference_permeability(irb))
      case(4)
        x = abs(log(perm_ce) - &
                log(constrained_block%reference_permeability(irb)))
      case(5)
        !perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        !x = log(perm_ce) - log(perm_nb)

        ! compute unit vectors: rx,ry, and rz
        rx = grid%x(ghosted_id) - grid%x(ghosted_id_nb)
        ry = grid%y(ghosted_id) - grid%y(ghosted_id_nb)
        rz = grid%z(ghosted_id) - grid%z(ghosted_id_nb)
        r = sqrt(rx*rx + ry*ry + rz*rz)
        rx = rx / r
        ry = ry / r
        rz = rz / r
      case(6)
        !perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        !x = abs(log(perm_ce) - log(perm_nb))

        ! compute unit vectors: rx,ry, and rz
        rx = abs(grid%x(ghosted_id) - grid%x(ghosted_id_nb))
        ry = abs(grid%y(ghosted_id) - grid%y(ghosted_id_nb))
        rz = abs(grid%z(ghosted_id) - grid%z(ghosted_id_nb))
        r = sqrt(rx*rx + ry*ry + rz*rz)
        rx = rx / r
        ry = ry / r
        rz = rz / r
      case(7)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = (log(perm_ce)- log(constrained_block%reference_permeability(irb))) &
          -(log(perm_nb) - log(constrained_block%reference_permeability(irb)))
      case(8)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = abs( &
            (log(perm_ce)- log(constrained_block%reference_permeability(irb))) &
          -(log(perm_nb) - log(constrained_block%reference_permeability(irb))) )
      case(9)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = log(perm_ce) - log(perm_nb)
      case(10)
        perm_nb = material_auxvars(ghosted_id_nb)%permeability(perm_xx_index)
        x = abs(log(perm_ce) - log(perm_nb))
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
  PetscInt :: iconblock,inbr,ilink
  PetscInt :: num_constraints
  PetscInt :: num_neighbor
  PetscErrorCode :: ierr

  option => this%realization%option
  patch => this%realization%patch
  grid => patch%grid

  constrained_block => this%constrained_block

  if (this%qoi_is_full_vector) then
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
  endif

  this%num_constraints_local = num_constraints
  call MPI_Allreduce(num_constraints,this%num_constraints_total, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  allocate(this%Wm(num_constraints))
  allocate(this%rblock(num_constraints,THREE_INTEGER))
  this%Wm = 0.d0
  this%rblock = 0

  if (this%qoi_is_full_vector) then

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
  use Inversion_Aux_module

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

  num_measurement = size(this%measurements)

  ! Data part
  call VecDuplicate(this%dist_parameter_vec,p1,ierr);CHKERRQ(ierr)
  call VecDuplicate(this%dist_measurement_vec,q1_dist,ierr);CHKERRQ(ierr)
  call VecDuplicate(this%measurement_vec,q1,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(p1,pvec_ptr,ierr);CHKERRQ(ierr)
  pvec_ptr = this%p
  call VecRestoreArrayF90(p1,pvec_ptr,ierr);CHKERRQ(ierr)

  ! q = Jp -> data part
  call MatMultTranspose(inversion_aux%JsensitivityT,p1,q1_dist, &
                        ierr);CHKERRQ(ierr)

  call InvSubsurfScatMeasToDistMeas(this, &
                                    q1, &
                                    q1_dist, &
                                    INVSUBSCATREVERSE)

  call VecGetArrayF90(q1,q1vec_ptr,ierr);CHKERRQ(ierr)
  this%q(1:num_measurement) = q1vec_ptr
  call VecRestoreArrayF90(q1,q1vec_ptr,ierr);CHKERRQ(ierr)

  ! Model part -> q2
  ! Get local this%p to ghosted in pvec_ptr
  call InvSubsurfScatGlobalToDistParam(this, &
                                       this%realization%field%work, &
                                       p1, &
                                       INVSUBSCATREVERSE)
  call DiscretizationGlobalToLocal(discretization,field%work, &
                                   field%work_loc,ONEDOF)
  call VecGetArrayF90(field%work_loc,pvec_ptr,ierr);CHKERRQ(ierr)

  beta = this%beta

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
  use Inversion_Aux_module

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

  num_measurement = size(this%measurements)

  ! Model part -> s2
  call VecGetArrayF90(field%work_loc,s2vec_ptr,ierr);CHKERRQ(ierr)
  s2vec_ptr = 0.d0

  beta = this%beta

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
  call InvSubsurfScatGlobalToDistParam(this, &
                                       this%realization%field%work, &
                                       this%dist_parameter_tmp_vec, &
                                       INVSUBSCATFORWARD)

  ! Data part
  call VecDuplicate(this%measurement_vec,r1,ierr);CHKERRQ(ierr)
  call VecDuplicate(this%dist_parameter_tmp_vec,s1,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(r1,r1vec_ptr,ierr);CHKERRQ(ierr)
  r1vec_ptr = this%r(1:num_measurement)
  call VecRestoreArrayF90(r1,r1vec_ptr,ierr);CHKERRQ(ierr)
  call InvSubsurfScatMeasToDistMeas(this, &
                                    r1, &
                                    this%dist_measurement_vec, &
                                    INVSUBSCATFORWARD)
  call VecGetArrayF90(this%dist_measurement_vec,r1vec_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(this%dist_measurement_vec,r1vec_ptr,ierr);CHKERRQ(ierr)

  ! s = J^T*r -> data part
  call MatMult(inversion_aux%JsensitivityT,this%dist_measurement_vec, &
               s1,ierr);CHKERRQ(ierr)

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

subroutine InversionZFlowWriteIterationInfo(this)
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
  character(len=MAXWORDLENGTH) :: string

  if (this%driver%PrintToScreen()) then
    write(*,*)
    write(*,98)
    if (this%iteration == this%start_iteration) then
      write(*,'(/,2x,a,i6.4,/)') StringColor("CONVERGENCE STATISTICS AT &
                              &STARTING ITERATION:",C_RED), this%start_iteration
    else
      write(*,'(/,2x,a,i6.4,/)') StringColor("CONVERGENCE STATISTICS AFTER &
                                 &ITERATION:",C_RED),this%iteration
    endif
    write(*,99)
    write(*,*) StringColor("  Phi_data   ",C_GREEN), &
               StringColor("   Phi_Model  ",C_BLUE), &
               StringColor("  Phi_Model/Beta",C_MAGENTA), &
               StringColor("    Phi_Total   ",C_CYAN)
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

end subroutine InversionZFlowWriteIterationInfo

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
  use Variables_module, only : PERMEABILITY

  class(inversion_zflow_type) :: this

  Vec :: wd_vec
  PetscInt :: idata,num_measurement
  PetscReal :: wd
  PetscReal, pointer :: wdvec_ptr(:)
  PetscErrorCode :: ierr

  num_measurement = size(this%measurements)
  call VecDuplicate(this%measurement_vec,wd_vec,ierr);CHKERRQ(ierr)
  call VecZeroEntries(wd_vec,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(wd_vec,wdvec_ptr,ierr);CHKERRQ(ierr)
  do idata = 1, num_measurement
    wd = 0.05 * this%measurements(idata)%value
    wd = 1/wd
    wdvec_ptr(idata) = wd
  enddo
  call VecRestoreArrayF90(wd_vec,wdvec_ptr,ierr);CHKERRQ(ierr)
  call InvSubsurfScatMeasToDistMeas(this, &
                                    wd_vec, &
                                    this%dist_measurement_vec, &
                                    INVSUBSCATFORWARD)

  ! Column Scale with wd
  call MatDiagonalScale(this%inversion_aux%JsensitivityT, &
                        PETSC_NULL_VEC, & ! scales rows
                        this%dist_measurement_vec, &  ! scales columns
                        ierr);CHKERRQ(ierr)
  ! Row scale with perm
  call MatDiagonalScale(this%inversion_aux%JsensitivityT, &
                        this%dist_parameter_vec, & ! scales rows
                        PETSC_NULL_VEC, &  ! scales columns
                        ierr);CHKERRQ(ierr)
  call VecDestroy(wd_vec,ierr);CHKERRQ(ierr)

end subroutine InversionZFlowScaleSensitivity

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
  call DeallocateArray(constrained_block%reference_permeability)

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
  call DeallocateArray(this%del_perm)
  call DeallocateArray(this%Wm)
  call DeallocateArray(this%rblock)

  if (this%dist_parameter_tmp_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%dist_parameter_tmp_vec,ierr);CHKERRQ(ierr)
  endif

  call ConstrainedBlockDestroy(this%constrained_block)

  call this%Strip()
  deallocate(this)
  nullify(this)

end subroutine InversionZFlowDestroy

end module Inversion_ZFlow_class
