module Inversion_ERT_class

#include "petsc/finclude/petscvec.h"
  use petscvec

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Inversion_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_subsurface_type) :: inversion_ert_type
    PetscInt :: start_iteration          ! Starting iteration number
    PetscInt :: miniter,maxiter          ! min/max CGLS iterations

    PetscReal :: beta                    ! regularization parameter
    PetscReal :: beta_red_factor         ! beta reduction factor
    PetscReal :: mincond,maxcond         ! min/max conductivity
    PetscReal :: target_chi2             ! target CHI^2 norm
    PetscReal :: current_chi2

    ! Cost/objective functions
    PetscReal :: min_phi_red             ! min change in cost function
    PetscReal :: phi_total_0,phi_total
    PetscReal :: phi_data_0,phi_data
    PetscReal :: phi_model_0,phi_model

    PetscBool :: cull_flag               ! flag to ignore data outliers
    PetscReal :: cull_dev                ! data cull cutoff (std deviation)

    PetscBool :: app_cond_start_model    ! apparent cond as start model

    ! arrays for CGLS algorithm
    PetscReal, pointer :: b(:)           ! vector for CGLS RHS
    PetscReal, pointer :: p(:)           ! vector of dim -> num of inv cells
    PetscReal, pointer :: q(:)           ! product of Jacobian with p = Jp
    PetscReal, pointer :: r(:)           ! vector of dim -> num of measur
    PetscReal, pointer :: s(:)           ! product Jacobian transpose with r
    PetscReal, pointer :: del_cond(:)    ! conductivity update vector

    ! For Wm
    PetscInt :: num_constraints_local    ! Number of constraints
    PetscInt :: num_constraints_total    ! Total number of constraints
    PetscInt, pointer :: rblock(:,:)     ! array stores info about reg.
    PetscReal, pointer :: Wm(:)          ! Regularization matrix

    type(constrained_block_type), pointer :: constrained_block

  contains
    procedure, public :: Init => InversionERTInit
    procedure, public :: Initialize => InversionERTInitialize
    procedure, public :: ReadBlock => InversionERTReadBlock
    procedure, public :: InitializeIterationNumber => &
                           InversionERTInitIterationNumber
    procedure, public :: Step => InversionERTStep
    procedure, public :: UpdateParameters => InversionERTUpdateParameters
    procedure, public :: CalculateUpdate => InversionERTCalculateUpdate
    procedure, public :: CheckConvergence => InversionERTCheckConvergence
    procedure, public :: EvaluateCostFunction => InvERTEvaluateCostFunction
    procedure, public :: UpdateRegularizeParameters => &
                           InvERTUpdateRegularizParams
    procedure, public :: WriteIterationInfo => InversionERTWriteIterationInfo
    procedure, public :: Finalize => InversionERTFinalize
    procedure, public :: Strip => InversionERTStrip
  end type inversion_ert_type

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
    PetscReal, pointer :: reference_conductivity(:)
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
    PetscReal :: reference_conductivity
    type(constrained_block_par_type), pointer :: next
  end type constrained_block_par_type

  public :: InversionERTCreate, &
            InversionERTStrip,  &
            InversionERTDestroy

contains

! ************************************************************************** !

function InversionERTCreate(driver)
  !
  ! Allocates and initializes a new inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  use Driver_module

  class(driver_type), pointer :: driver

  class(inversion_ert_type), pointer :: InversionERTCreate

  allocate(InversionERTCreate)
  call InversionERTCreate%Init(driver)

end function InversionERTCreate

! ************************************************************************** !

subroutine InversionERTInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY
  use Driver_module

  class(inversion_ert_type) :: this
  class(driver_type), pointer :: driver

  call InversionSubsurfaceInit(this,driver)
  ! override default set in InversionSubsurfaceInit
  this%iqoi = ELECTRICAL_CONDUCTIVITY

  ! Default inversion parameters
  this%miniter = 10
  this%maxiter = 50

  this%beta = 100.d0
  this%beta_red_factor = 0.5d0
  this%mincond = 0.00001d0
  this%maxcond = 10.d0
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

  this%cull_flag = PETSC_FALSE
  this%cull_dev = UNINITIALIZED_DOUBLE

  this%app_cond_start_model = PETSC_TRUE

  nullify(this%b)
  nullify(this%p)
  nullify(this%q)
  nullify(this%r)
  nullify(this%s)
  nullify(this%del_cond)
  nullify(this%Wm)
  nullify(this%rblock)

  this%constrained_block => ConstrainedBlockCreate()

  nullify(this%realization)

end subroutine InversionERTInit

! ************************************************************************** !

function ConstrainedBlockCreate()
  !
  ! Creates Constrained Block type
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
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
  nullify(constrained_block%reference_conductivity)

  ConstrainedBlockCreate => constrained_block

end function ConstrainedBlockCreate

! ************************************************************************** !

function ConstrainedBlockParCreate()
  !
  ! Creates Constrained Block Par type
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
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
  constrained_block%reference_conductivity = 0.d0
  nullify(constrained_block%next)

  ConstrainedBlockParCreate => constrained_block

end function ConstrainedBlockParCreate

! ************************************************************************** !

subroutine InversionERTAllocateWorkArrays(this)
  !
  ! Initialize inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/17/21
  !

  use Grid_module
  use Survey_module

  implicit none

  class(inversion_ert_type) :: this

  type(survey_type), pointer :: survey
  type(grid_type), pointer :: grid

  PetscInt :: num_constraints

  survey => this%realization%survey
  grid => this%realization%patch%grid

  num_constraints = this%num_constraints_local

  allocate(this%b(survey%num_measurement + num_constraints))
  allocate(this%p(grid%nlmax))
  allocate(this%q(survey%num_measurement + num_constraints))
  allocate(this%r(survey%num_measurement + num_constraints))
  allocate(this%s(grid%nlmax))
  allocate(this%del_cond(grid%nlmax))

  this%b = 0.d0
  this%p = 0.d0
  this%q = 0.d0
  this%r = 0.d0
  this%s = 0.d0
  this%del_cond = 0.d0

end subroutine InversionERTAllocateWorkArrays

! ************************************************************************** !

subroutine InversionERTDeallocateWorkArrays(this)
  !
  ! Initialize inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/17/21
  !

  use Utility_module, only : DeallocateArray

  implicit none

  class(inversion_ert_type) :: this

  call DeallocateArray(this%b)
  call DeallocateArray(this%p)
  call DeallocateArray(this%q)
  call DeallocateArray(this%r)
  call DeallocateArray(this%s)
  call DeallocateArray(this%del_cond)

end subroutine InversionERTDeallocateWorkArrays

! ************************************************************************** !

subroutine InversionERTConstrainedArraysFromList(this)
  !
  ! Gets Constrained Block parameter arrays from linked list
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !

  use Option_module
  use Material_module
  use Patch_module

  implicit none

  class(inversion_ert_type) :: this

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
    allocate(constrained_block%reference_conductivity(nconblock))
    constrained_block%reference_conductivity = 0.d0
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
      constrained_block%reference_conductivity(iconblock) = &
                        cur_constrained_block%reference_conductivity
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

end subroutine InversionERTConstrainedArraysFromList

! ************************************************************************** !

subroutine InversionERTReadBlock(this,input,option)
  !
  ! Reads input file parameters associated an ERT inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use Variables_module, only : PERMEABILITY, ELECTRICAL_CONDUCTIVITY

  class(inversion_ert_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'INVERSION'

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
    call InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('MIN_ELECTRICAL_CONDUCTIVITY')
        call InputReadDouble(input,option,this%mincond)
        call InputErrorMsg(input,option,'MIN_ELECTRICAL_CONDUCTIVITY', &
                           error_string)
      case('MAX_ELECTRICAL_CONDUCTIVITY')
        call InputReadDouble(input,option,this%maxcond)
        call InputErrorMsg(input,option,'MAX_ELECTRICAL_CONDUCTIVITY', &
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
      case('STARTING_MODEL')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringToUpper(word)
        select case(word)
          case('APPARENT_ELECTRICAL_CONDUCTIVITY')
            this%app_cond_start_model = PETSC_TRUE
          case('INPUT_ELECTRICAL_CONDUCTIVITY')
            this%app_cond_start_model = PETSC_FALSE
          case default
            call InputKeywordUnrecognized(input,word,trim(error_string)// &
                             & ',STARTING_MODEL',option)
        end select
      case('START_INVERSION_ITERATION')
        call InputReadInt(input,option,this%start_iteration)
        call InputErrorMsg(input,option,'START_INVERSION_ITERATION', &
                           error_string)
        this%iteration = this%start_iteration
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine InversionERTReadBlock

! ************************************************************************** !

subroutine ConstrainedBlockRead(constrained_block,input,option)
  !
  ! Read constrained blocks options
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21

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
  ! Date: 06/14/21

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
      case('REFERENCE_CONDUCTIVITY')
        call InputReadDouble(input,option, &
                            constrained_block%reference_conductivity)
        call InputErrorMsg(input,option,'REFERENCE_CONDUCTIVITY',error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine ConstrainedBlockParRead

! ************************************************************************** !

subroutine InversionERTInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  use Dataset_Base_class
  use Dataset_module
  use Discretization_module
  use Init_Subsurface_module
  use Material_module
  use Option_module
  use Variables_module, only : PERMEABILITY, ELECTRICAL_CONDUCTIVITY

  class(inversion_ert_type) :: this

  PetscBool :: exists
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string
  class(dataset_base_type), pointer :: dataset
  PetscErrorCode :: ierr

  ! theck to ensure that quantity of interest exists
  exists = PETSC_FALSE
  select case(this%iqoi)
    case(ELECTRICAL_CONDUCTIVITY)
      if (this%realization%option%igeopmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'ELECTRICAL_CONDUCTIVITY'
    case default
  end select
  if (.not.exists) then
    this%realization%option%io_buffer = 'Inversion for ' // trim(word) // &
      &' cannot be performed with the specified process models.'
    call PrintErrMsg(this%realization%option)
  endif

  if (this%app_cond_start_model) then
    ! non-ghosted Vec
    call VecDuplicate(this%realization%field%work, &
                      this%quantity_of_interest,ierr);CHKERRQ(ierr)
    call VecSet(this%quantity_of_interest, &
            this%realization%survey%apparent_conductivity,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(this%realization%discretization, &
                                     this%quantity_of_interest, &
                                     this%realization%field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi,ZERO_INTEGER)
  else
    ! non-ghosted Vec
    call VecDuplicate(this%realization%field%work, &
                      this%quantity_of_interest,ierr);CHKERRQ(ierr)
    ! ghosted Vec
    call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(this%realization%discretization, &
                                     this%realization%field%work_loc, &
                                     this%quantity_of_interest,ONEDOF)
  endif

  call InversionERTConstrainedArraysFromList(this)

  if (len_trim(this%ref_qoi_dataset_name) > 0) then
    call VecDuplicate(this%quantity_of_interest, &
                      this%ref_quantity_of_interest,ierr);CHKERRQ(ierr)
    string = 'Reference QOI dataset'
    dataset => DatasetBaseGetPointer(this%realization%datasets, &
                                     this%ref_qoi_dataset_name, &
                                     string,this%realization%option)
    call SubsurfReadDatasetToVecWithMask(this%realization,dataset, &
                                         ZERO_INTEGER,PETSC_TRUE, &
                                         this%ref_quantity_of_interest)
    ! do not destroy the dataset as the pointer is owned by the list; just
    ! strip to free up memory
    call DatasetStrip(dataset)
  endif

end subroutine InversionERTInitialize

! ************************************************************************** !

subroutine InversionERTStep(this)
  !
  ! Execute a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Option_module
  use Factory_Forward_module

  class(inversion_ert_type) :: this

  type(option_type), pointer :: option

  option => OptionCreate()
  write(option%group_prefix,'(i6)') this%iteration+1
  option%group_prefix = 'Run' // trim(adjustl(option%group_prefix))
  call OptionSetDriver(option,this%driver)
  call FactoryForwardInitialize(this%forward_simulation, &
                                this%forward_simulation_filename,option)
  this%realization => this%forward_simulation%realization
  call this%UpdateParameters()
  call this%forward_simulation%InitializeRun()
  if (option%status == PROCEED) then
    call this%forward_simulation%ExecuteRun()
  endif
  call this%CheckConvergence()
  call this%CalculateUpdate()
  call this%WriteIterationInfo()
  call this%UpdateRegularizationParameters()
  call this%forward_simulation%FinalizeRun()
  call this%forward_simulation%Strip()
  deallocate(this%forward_simulation)
  nullify(this%forward_simulation)

end subroutine InversionERTStep

! ************************************************************************** !

subroutine InversionERTCheckConvergence(this)
  !
  ! Check Inversion convergence
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/15/21

  use Survey_module

  implicit none

  class(inversion_ert_type) :: this

  type(survey_type), pointer :: survey

  survey => this%realization%survey

  this%converg_flag = PETSC_FALSE
  call this%EvaluateCostFunction()
  if ((this%current_chi2 <= this%target_chi2) .or. &
      (this%iteration > this%maximum_iteration)) this%converg_flag = PETSC_TRUE

end subroutine InversionERTCheckConvergence

! ************************************************************************** !

subroutine InvERTEvaluateCostFunction(this)
  !
  ! Evaluates cost functions for inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/16/21

  use Option_module
  use Patch_module
  use Material_Aux_class
  use Survey_module

  implicit none

  class(inversion_ert_type) :: this

  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(survey_type), pointer :: survey
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: idata,ndata,ncull
  PetscInt :: iconst,num_constraints
  PetscInt :: irb,ghosted_id,ghosted_id_nb
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: err_mean,err_sdev
  PetscReal :: cond_ce,cond_nb              ! cell's and neighbor's
  PetscReal :: wm,x
  PetscReal, allocatable :: data_vector(:)
  PetscReal, allocatable :: model_vector(:)
  PetscErrorCode :: ierr

  option => this%realization%option
  patch => this%realization%patch
  material_auxvars => patch%aux%Material%auxvars

  survey => this%realization%survey
  constrained_block => this%constrained_block
  rblock => this%rblock

  ndata = survey%num_measurement
  allocate(data_vector(ndata))
  data_vector = 0.d0

  data_vector = survey%Wd * (survey%dobs - survey%dsim)

  ncull = 0
  if (this%cull_flag) then
    err_mean = sum(data_vector) / ndata
    err_sdev = sqrt( dot_product( (data_vector - err_mean), &
                                  (data_vector - err_mean) ) / ndata )
    do idata=1,ndata
      if (data_vector(idata) < (err_mean - err_sdev*this%cull_dev) .or. &
          data_vector(idata) > (err_mean + err_sdev*this%cull_dev)) then
        survey%Wd_cull(idata) = 0
        ncull = ncull + 1
      endif
    enddo
  endif

  data_vector = survey%Wd_cull * data_vector

  this%phi_data = dot_product(data_vector,data_vector)
  this%current_chi2 = this%phi_data / (ndata - ncull)

  deallocate(data_vector)

  ! model cost function
  num_constraints = this%num_constraints_local
  allocate(model_vector(num_constraints))
  model_vector = 0.d0

  do iconst=1,num_constraints
    if (this%Wm(iconst) == 0) cycle

    wm = this%Wm(iconst)

    ! get cond & block of the ith constrained eq.
    ghosted_id = rblock(iconst,1)
    ghosted_id_nb = rblock(iconst,2)
    irb = rblock(iconst,3)
    cond_ce = material_auxvars(ghosted_id)%electrical_conductivity(1)
    x = 0.d0

    select case(constrained_block%structure_metric(irb))
      case(1)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = log(cond_ce) - log(cond_nb)
      case(2)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = abs(log(cond_ce) - log(cond_nb))
      case(3)
        x = log(cond_ce) - log(constrained_block%reference_conductivity(irb))
      case(4)
        x = abs(log(cond_ce) - &
                log(constrained_block%reference_conductivity(irb)))
      case(5)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = log(cond_ce) - log(cond_nb)
      case(6)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = abs(log(cond_ce) - log(cond_nb))
      case(7)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = (log(cond_ce) - log(constrained_block%reference_conductivity(irb))) &
          -(log(cond_nb) - log(constrained_block%reference_conductivity(irb)))
      case(8)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = abs( &
            (log(cond_ce) - log(constrained_block%reference_conductivity(irb))) &
          -(log(cond_nb) - log(constrained_block%reference_conductivity(irb))) )
      case(9)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = log(cond_ce) - log(cond_nb)
      case(10)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = abs(log(cond_ce) - log(cond_nb))
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

end subroutine InvERTEvaluateCostFunction

! ************************************************************************** !

subroutine InvERTUpdateRegularizParams(this)
  !
  ! Check Beta if it needs cooling/reduction
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/15/21

  implicit none

  class(inversion_ert_type) :: this

  ! update iteration number
  this%iteration = this%iteration + 1

  if (this%iteration - 1 == this%start_iteration) return

  if ( (this%phi_total_0 - this%phi_total)/this%phi_total_0 <= &
                                                      this%min_phi_red ) then
    this%beta = this%beta * this%beta_red_factor
    this%phi_model = this%beta_red_factor * this%phi_model
  endif

  ! update the cost functions
  this%phi_data_0 = this%phi_data
  this%phi_model_0 = this%phi_model
  this%phi_total_0 = this%phi_total

end subroutine InvERTUpdateRegularizParams

! ************************************************************************** !

subroutine InversionERTInitIterationNumber(this)
  !
  ! Sets starting iteration number
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 07/09/21

  class(inversion_ert_type) :: this

  this%iteration = this%start_iteration - 1

end subroutine InversionERTInitIterationNumber

! ************************************************************************** !

subroutine InversionERTUpdateParameters(this)
  !
  ! Updates input parameters
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !

  use Material_module
  use Discretization_module
  use Field_module

  class(inversion_ert_type) :: this

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization

  PetscInt :: local_id,ghosted_id
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  field => this%realization%field
  discretization => this%realization%discretization

  if (this%quantity_of_interest == PETSC_NULL_VEC) then
    call this%Initialize()
  else
    call DiscretizationGlobalToLocal(discretization, &
                                     this%quantity_of_interest, &
                                     field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 field%work_loc,this%iqoi,ZERO_INTEGER)
  endif

  ! Build Wm matrix
  call InversionERTBuildWm(this)

end subroutine InversionERTUpdateParameters

! ************************************************************************** !

subroutine InversionERTCalculateUpdate(this)
  !
  ! Calculates updated model parameters
  ! using m_new = m_old + del_m
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !

  use Patch_module
  use Grid_module

  implicit none

  class(inversion_ert_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid

  PetscInt :: local_id
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid

  if (this%quantity_of_interest /= PETSC_NULL_VEC) then

    call InversionERTAllocateWorkArrays(this)

    ! get inversion%del_cond
    call InversionERTCGLSSolve(this)

    ! Get updated conductivity as m_new = m_old + del_m (where m = log(sigma))
    call VecGetArrayF90(this%quantity_of_interest,vec_ptr,ierr);CHKERRQ(ierr)
    do local_id=1,grid%nlmax
      vec_ptr(local_id) = exp(log(vec_ptr(local_id)) + this%del_cond(local_id))
      if (vec_ptr(local_id) > this%maxcond) vec_ptr(local_id) = this%maxcond
      if (vec_ptr(local_id) < this%mincond) vec_ptr(local_id) = this%mincond
    enddo
    call VecRestoreArrayF90(this%quantity_of_interest,vec_ptr, &
                                                          ierr);CHKERRQ(ierr)
    call InversionERTDeallocateWorkArrays(this)
  endif

end subroutine InversionERTCalculateUpdate

! ************************************************************************** !

subroutine InversionERTCGLSSolve(this)
  !
  ! Implements CGLS solver for least sqaure equivalent
  !            of the normal equations
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/17/21
  !

  use Option_module
  use Survey_module
  use Timer_class
  use String_module

  implicit none

  class(inversion_ert_type) :: this

  type(option_type), pointer :: option
  type(survey_type), pointer :: survey
  class(timer_type), pointer ::timer

  PetscInt :: i,nm,ncons
  PetscReal :: alpha,gbeta,gamma,gamma1,delta1,delta2,delta
  PetscReal :: norms0,norms,normx,xmax
  PetscReal :: resNE,resNE_old
  PetscBool :: exit_info,indefinite
  PetscErrorCode :: ierr

  PetscReal, parameter :: delta_initer = 1e-3
  PetscReal, parameter :: initer_conv  = 1e-4

  option => this%realization%option
  survey => this%realization%survey

  this%del_cond = 0.0d0

  timer => TimerCreate()
  call timer%Start()

  if (OptionPrintToScreen(option)) then
    write(*,'(" --> Solving normal equation using CGLS solver:")')
  endif

  nm = survey%num_measurement
  ncons = this%num_constraints_local

  ! Get RHS vector this%b
  call InversionERTCGLSRhs(this)

  this%r = this%b

  ! get this%s = J^tr
  call InversionERTComputeMatVecProductJtr(this)
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
    call InversionERTComputeMatVecProductJp(this)

    delta1 = dot_product(this%q(1:nm),this%q(1:nm))
    delta2 = dot_product(this%q(nm+1:nm+ncons),this%q(nm+1:nm+ncons))
    call MPI_Allreduce(MPI_IN_PLACE,delta2,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
    delta = delta1 + delta2

    if (delta < 0) indefinite = PETSC_TRUE
    if (delta == 0) delta = epsilon(delta)

    alpha = gamma / delta

    this%del_cond = this%del_cond + alpha * this%p
    this%r = this%r - alpha * this%q

    ! get this%s = J^tr
    call InversionERTComputeMatVecProductJtr(this)

    gamma1 = gamma
    gamma = dot_product(this%s,this%s)
    call MPI_Allreduce(MPI_IN_PLACE,gamma,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

    norms = sqrt(gamma)
    gbeta = gamma / gamma1
    this%p = this%s + gbeta * this%p

    normx = dot_product(this%del_cond,this%del_cond)
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
    ' iterations to solve normal equation.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InversionERTCGLSSolve

! ************************************************************************** !

subroutine InversionERTCGLSRhs(this)
  !
  ! Builds RHS for least-square equation for CGLS solver
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/17/21
  !

  use Patch_module
  use Material_Aux_class
  use Option_module
  use Survey_module

  implicit none

  class(inversion_ert_type) :: this

  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(survey_type), pointer :: survey
  type(constrained_block_type), pointer :: constrained_block

  PetscInt :: idata,iconst,irb,num_measurement
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: cond_ce,cond_nb,x     ! cell's and neighbor's
  PetscReal :: wm,beta

  option => this%realization%option
  patch => this%realization%patch
  material_auxvars => patch%aux%Material%auxvars

  survey => this%realization%survey
  constrained_block => this%constrained_block
  rblock => this%rblock

  this%b = 0.0d0

  num_measurement = survey%num_measurement

  do idata=1,num_measurement
    this%b(idata) = survey%Wd(idata) * survey%Wd_cull(idata) * &
                         ( survey%dobs(idata) - survey%dsim(idata) )
  enddo

  beta = this%beta

  do iconst=1,this%num_constraints_local
    if (this%Wm(iconst) == 0) cycle

    wm = this%Wm(iconst)

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

    this%b(num_measurement + iconst) = - sqrt(beta) * wm * x

  enddo

end subroutine InversionERTCGLSRhs

! ************************************************************************** !

subroutine InversionERTBuildWm(this)
  !
  ! Builds model regularization matrix: Wm
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/18/21
  !

  use Patch_module
  use Grid_module
  use Material_Aux_class
  use Option_module

  implicit none

  class(inversion_ert_type) :: this

  class(material_auxvar_type), pointer :: material_auxvars(:)
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
                                          call InversionERTAllocateWm(this)

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
    PetscReal :: cond_ce,cond_nb     ! cell's and neighbor's
    PetscReal :: mn,sd
    PetscReal :: rx,ry,rz,r
    PetscInt, pointer :: rblock(:,:)

    rblock => this%rblock

    ! get cond & block of the ith constrained eq.
    ghosted_id = rblock(iconst,1)
    ghosted_id_nb = rblock(iconst,2)
    irb = rblock(iconst,3)
    cond_ce = material_auxvars(ghosted_id)%electrical_conductivity(1)
    x = 0.d0

    select case(constrained_block%structure_metric(irb))
      case(1)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = log(cond_ce) - log(cond_nb)
      case(2)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = abs(log(cond_ce) - log(cond_nb))
      case(3)
        x = log(cond_ce) - log(constrained_block%reference_conductivity(irb))
      case(4)
        x = abs(log(cond_ce) - &
                log(constrained_block%reference_conductivity(irb)))
      case(5)
        !cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        !x = log(cond_ce) - log(cond_nb)

        ! compute unit vectors: rx,ry, and rz
        rx = grid%x(ghosted_id) - grid%x(ghosted_id_nb)
        ry = grid%y(ghosted_id) - grid%y(ghosted_id_nb)
        rz = grid%z(ghosted_id) - grid%z(ghosted_id_nb)
        r = sqrt(rx*rx + ry*ry + rz*rz)
        rx = rx / r
        ry = ry / r
        rz = rz / r
      case(6)
        !cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        !x = abs(log(cond_ce) - log(cond_nb))

        ! compute unit vectors: rx,ry, and rz
        rx = abs(grid%x(ghosted_id) - grid%x(ghosted_id_nb))
        ry = abs(grid%y(ghosted_id) - grid%y(ghosted_id_nb))
        rz = abs(grid%z(ghosted_id) - grid%z(ghosted_id_nb))
        r = sqrt(rx*rx + ry*ry + rz*rz)
        rx = rx / r
        ry = ry / r
        rz = rz / r
      case(7)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = (log(cond_ce) - log(constrained_block%reference_conductivity(irb))) &
          -(log(cond_nb) - log(constrained_block%reference_conductivity(irb)))
      case(8)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = abs( &
            (log(cond_ce) - log(constrained_block%reference_conductivity(irb))) &
          -(log(cond_nb) - log(constrained_block%reference_conductivity(irb))) )
      case(9)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
        x = log(cond_ce) - log(cond_nb)
      case(10)
        cond_nb = material_auxvars(ghosted_id_nb)%electrical_conductivity(1)
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

end subroutine InversionERTBuildWm

! ************************************************************************** !

subroutine InversionERTAllocateWm(this)
  !
  ! Allocate and get info on Wm and rblock
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/17/21
  !

  use Patch_module
  use Grid_module
  use Option_module

  implicit none

  class(inversion_ert_type) :: this

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

  this%num_constraints_local = num_constraints
  call MPI_Allreduce(num_constraints,this%num_constraints_total, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  allocate(this%Wm(num_constraints))
  allocate(this%rblock(num_constraints,THREE_INTEGER))
  this%Wm = 0.d0
  this%rblock = 0

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
            this%rblock(num_constraints,1) = ghosted_id
            this%rblock(num_constraints,3) = iconblock
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

end subroutine InversionERTAllocateWm

! ************************************************************************** !

subroutine InversionERTComputeMatVecProductJp(this)
  !
  ! Computes product of Jacobian J with a vector p = Jp
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/15/21
  !

  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module
  use Option_module
  use Survey_module
  use ERT_Aux_module

  implicit none

  class(inversion_ert_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(survey_type), pointer :: survey
  type(constrained_block_type), pointer :: constrained_block
  type(ert_auxvar_type), pointer :: ert_auxvars(:)

  PetscInt :: idata,iconst,irb,num_measurement
  PetscInt :: local_id,ghosted_id
  PetscInt :: ghosted_id_nb
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: beta,wm
  PetscReal, pointer :: pvec_ptr(:)
  PetscErrorCode :: ierr

  option => this%realization%option
  field => this%realization%field
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  ert_auxvars => patch%aux%ERT%auxvars

  survey => this%realization%survey
  constrained_block => this%constrained_block
  rblock => this%rblock

  this%q = 0.d0

  ! Data part
  call VecGetArrayF90(field%work,pvec_ptr,ierr);CHKERRQ(ierr)
  pvec_ptr = 0.d0

  do idata=1,survey%num_measurement
    do local_id=1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      this%q(idata) = this%q(idata) + &
                           ert_auxvars(ghosted_id)%jacobian(idata) * &
                           this%p(local_id)
      if (idata == 1) pvec_ptr(local_id) = this%p(local_id)
    enddo

    call MPI_Allreduce(MPI_IN_PLACE,this%q(idata),ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
  enddo

  call VecRestoreArrayF90(field%work,pvec_ptr,ierr);CHKERRQ(ierr)

  ! Model part
  ! Get local this%p to ghosted in pvec_ptr
  call DiscretizationGlobalToLocal(discretization,field%work, &
                                   field%work_loc,ONEDOF)
  call VecGetArrayF90(field%work_loc,pvec_ptr,ierr);CHKERRQ(ierr)

  num_measurement = survey%num_measurement
  beta = this%beta

  do iconst=1,this%num_constraints_local
    if (this%Wm(iconst) == 0) cycle

    wm = this%Wm(iconst)
    irb = rblock(iconst,3)
    ghosted_id = rblock(iconst,1)

    if (constrained_block%structure_metric(irb) == 3 .or. &
        constrained_block%structure_metric(irb) == 4) then
      this%q(num_measurement + iconst) = &
        sqrt(beta) * wm * pvec_ptr(ghosted_id)
    else
      ghosted_id_nb = rblock(iconst,2)
      this%q(num_measurement + iconst) = &
          sqrt(beta) * wm * (pvec_ptr(ghosted_id) - pvec_ptr(ghosted_id_nb))
    endif
  enddo

  call VecRestoreArrayF90(field%work_loc,pvec_ptr,ierr);CHKERRQ(ierr)

end subroutine InversionERTComputeMatVecProductJp

! ************************************************************************** !

subroutine InversionERTComputeMatVecProductJtr(this)
  !
  ! Computes product of Jacobian J transpose with a vector r = J^t x r
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/16/21
  !

  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module
  use Survey_module
  use ERT_Aux_module

  implicit none

  class(inversion_ert_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(survey_type), pointer :: survey
  type(constrained_block_type), pointer :: constrained_block
  type(ert_auxvar_type), pointer :: ert_auxvars(:)

  PetscInt :: idata,iconst,irb,num_measurement
  PetscInt :: local_id,ghosted_id
  PetscInt :: ghosted_id_nb
  PetscInt, pointer :: rblock(:,:)
  PetscReal :: beta,wm
  PetscReal, pointer :: svec_ptr(:)
  PetscErrorCode :: ierr

  field => this%realization%field
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  ert_auxvars => patch%aux%ERT%auxvars

  survey => this%realization%survey
  constrained_block => this%constrained_block
  rblock => this%rblock

  this%s = 0.0d0

  ! Model part
  call VecGetArrayF90(field%work_loc,svec_ptr,ierr);CHKERRQ(ierr)
  svec_ptr = 0.d0

  num_measurement = survey%num_measurement
  beta = this%beta

  do iconst=1,this%num_constraints_local
    if (this%Wm(iconst) == 0) cycle

    wm = this%Wm(iconst)
    irb = rblock(iconst,3)
    ghosted_id = rblock(iconst,1)

    if (constrained_block%structure_metric(irb) == 3 .or. &
        constrained_block%structure_metric(irb) == 4) then
      svec_ptr(ghosted_id) = svec_ptr(ghosted_id) + &
        sqrt(beta) * wm * this%r(num_measurement + iconst)
    else
      ghosted_id_nb = rblock(iconst,2)
      svec_ptr(ghosted_id) = svec_ptr(ghosted_id) + &
              sqrt(beta) * wm * this%r(num_measurement + iconst)
      svec_ptr(ghosted_id_nb) = svec_ptr(ghosted_id_nb) - &
              sqrt(beta) * wm * this%r(num_measurement + iconst)
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
                             this%r(idata)
    enddo
    this%s(local_id) = svec_ptr(local_id)
  enddo

  call VecRestoreArrayF90(field%work,svec_ptr,ierr);CHKERRQ(ierr)

end subroutine InversionERTComputeMatVecProductJtr

! ************************************************************************** !

subroutine InversionERTWriteIterationInfo(this)
  !
  ! Writes inversion run info
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 07/02/21
  !

  use String_module

  implicit none

  class(inversion_ert_type) :: this

  PetscInt :: fid
  character(len=MAXWORDLENGTH) :: string

  if (this%driver%PrintToScreen()) then
    write(*,*)
    write(*,98)
    if (this%iteration == this%start_iteration) then
      write(*,'(/,2x,a,i6.4,/)') StringColor("CONVERGENCE STATISTICS AT STARTING &
                                 &ITERATION:",C_RED), this%start_iteration
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
    write(*,103) this%num_constraints_total
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
    write(fid,103) this%num_constraints_total
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

end subroutine InversionERTWriteIterationInfo

! ************************************************************************** !

subroutine InversionERTFinalize(this)
  !
  ! Finalizes inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  class(inversion_ert_type) :: this

  call InversionBaseFinalize(this)

end subroutine InversionERTFinalize

! ************************************************************************** !

subroutine InversionERTStrip(this)
  !
  ! Deallocates members of inversion ERT
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  class(inversion_ert_type) :: this

  PetscErrorCode :: ierr

  call InversionSubsurfaceStrip(this)

  nullify(this%realization)
  if (this%quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%quantity_of_interest,ierr);CHKERRQ(ierr)
  endif
  if (this%ref_quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%ref_quantity_of_interest,ierr);CHKERRQ(ierr)
  endif

end subroutine InversionERTStrip

! ************************************************************************** !

subroutine ConstrainedBlockParDestroy(constrained_block)
  !
  ! ConstrainedBlockParDestroy: Deallocates a constrained block par object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/16/21
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
  ! Date: 06/16/21
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
  call DeallocateArray(constrained_block%reference_conductivity)

  deallocate(constrained_block)
  nullify(constrained_block)

end subroutine ConstrainedBlockDestroy

! ************************************************************************** !

subroutine InversionERTDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !

  use Utility_module, only : DeallocateArray

  implicit none

  class(inversion_ert_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call DeallocateArray(inversion%b)
  call DeallocateArray(inversion%p)
  call DeallocateArray(inversion%q)
  call DeallocateArray(inversion%r)
  call DeallocateArray(inversion%s)
  call DeallocateArray(inversion%del_cond)
  call DeallocateArray(inversion%Wm)
  call DeallocateArray(inversion%rblock)

  call ConstrainedBlockDestroy(inversion%constrained_block)

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionERTDestroy

end module Inversion_ERT_class
