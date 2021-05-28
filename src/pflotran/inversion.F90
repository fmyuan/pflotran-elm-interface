module Inversion_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Survey_module

  implicit none

  private
    
  type, public :: inversion_type
    PetscInt :: iteration                ! iteration number
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

    PetscBool :: converg_flag            ! convergence flag

    ! arrays for CGLS algorithm
    PetscReal, pointer :: b(:)           ! vector for CGLS RHS
    PetscReal, pointer :: p(:)           ! vector of dim -> num of inv cells
    PetscReal, pointer :: q(:)           ! product of Jacobian with p = Jp
    PetscReal, pointer :: r(:)           ! vector of dim -> num of measur
    PetscReal, pointer :: s(:)           ! product Jacobian transpose with r
    PetscReal, pointer :: del_cond(:)    ! conductivity update vector

    ! For Wm
    PetscInt :: num_constraints_local    ! Number of constraints
    PetscInt, pointer :: rblock(:,:)     ! array stores info about reg.
    PetscReal, pointer :: Wm(:)          ! Regularization matrix
  
    type(constrained_block_type), pointer :: constrained_block
  
    end type inversion_type

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

  public :: InversionCreate, &
            InversionRead, &
            InversionConstrainedArraysFromList, &
            InversionAllocateWorkArrays, &
            InversionDeallocateWorkArrays, &
            InversionDestroy

contains

! ************************************************************************** !

function InversionCreate()
  !
  ! Creates inversion type
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/03/21
  !

  implicit none
      
  type(inversion_type), pointer :: InversionCreate
  type(inversion_type), pointer :: inversion

  allocate(inversion)

  ! Default inversion parameters
  inversion%miniter = 10
  inversion%maxiter = 50

  inversion%beta = 100.d0
  inversion%beta_red_factor = 0.5d0
  inversion%mincond = 0.00001d0
  inversion%maxcond = 10.d0
  inversion%target_chi2 = 1.d0
  inversion%min_phi_red = 0.2d0

  inversion%iteration = UNINITIALIZED_INTEGER
  inversion%num_constraints_local = UNINITIALIZED_INTEGER
  inversion%current_chi2 = UNINITIALIZED_DOUBLE
  inversion%phi_total_0 = UNINITIALIZED_DOUBLE
  inversion%phi_data_0 = UNINITIALIZED_DOUBLE
  inversion%phi_model_0 = UNINITIALIZED_DOUBLE
  inversion%phi_total = UNINITIALIZED_DOUBLE
  inversion%phi_data = UNINITIALIZED_DOUBLE
  inversion%phi_model = UNINITIALIZED_DOUBLE

  inversion%cull_flag = PETSC_FALSE
  inversion%cull_dev = UNINITIALIZED_DOUBLE

  inversion%converg_flag = PETSC_FALSE

  nullify(inversion%b)
  nullify(inversion%p)
  nullify(inversion%q)
  nullify(inversion%r)
  nullify(inversion%s)
  nullify(inversion%del_cond)
  nullify(inversion%Wm)
  nullify(inversion%rblock)

  inversion%constrained_block => ConstrainedBlockCreate()

  InversionCreate => inversion

end function InversionCreate

! ************************************************************************** !

function ConstrainedBlockCreate()
  !
  ! Creates Constrained Block type
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/20/21
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
  ! Date: 05/20/21
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

subroutine InversionConstrainedArraysFromList(inversion,patch,option)
  !
  ! Gets Constrained Block parameter arrays from linked list
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/24/21
  !

  use Option_module
  use Material_module
  use Patch_module

  implicit none

  type(inversion_type) :: inversion
  type(patch_type) :: patch
  type(option_type) :: option

  type(material_property_type), pointer :: material_property
  type(constrained_block_type), pointer :: constrained_block
  type(constrained_block_par_type), pointer :: cur_constrained_block

  PetscInt :: i,iconblock
  PetscInt :: nconblock

  constrained_block => inversion%constrained_block
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
                           ' // trim(cur_constrained_block%name) // &
                           '" not found in material list'
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
                             '//trim(cur_constrained_block%block_link(i)) // &
                             '" in contrained block "&
                             '//trim(cur_constrained_block%name) // &
                             '" not found in material list'
          call PrintErrMsg(option)
        endif
        constrained_block%block_link(iconblock,i+1) = &
                                      material_property%internal_id
      enddo

      cur_constrained_block => cur_constrained_block%next
      iconblock = iconblock + 1

    enddo
  endif

end subroutine InversionConstrainedArraysFromList

! ************************************************************************** !

subroutine InversionAllocateWorkArrays(inversion,survey,grid)
  ! 
  ! Initialize inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/03/21
  !

  use Grid_module
       
  implicit none
      
  type(inversion_type) :: inversion
  type(survey_type) :: survey
  type(grid_type), pointer :: grid

  PetscInt :: num_constraints

  num_constraints = inversion%num_constraints_local

  allocate(inversion%b(survey%num_measurement + num_constraints))
  allocate(inversion%p(grid%nlmax))
  allocate(inversion%q(survey%num_measurement + num_constraints))
  allocate(inversion%r(survey%num_measurement + num_constraints))
  allocate(inversion%s(grid%nlmax))
  allocate(inversion%del_cond(grid%nlmax))

  inversion%b = 0.d0
  inversion%p = 0.d0
  inversion%q = 0.d0
  inversion%r = 0.d0
  inversion%s = 0.d0
  inversion%del_cond = 0.d0

end subroutine InversionAllocateWorkArrays

! ************************************************************************** !

subroutine InversionDeallocateWorkArrays(inversion)
  ! 
  ! Initialize inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/11/21
  !
  
  use Utility_module, only : DeallocateArray

  implicit none
      
  type(inversion_type) :: inversion

  call DeallocateArray(inversion%b)
  call DeallocateArray(inversion%p)
  call DeallocateArray(inversion%q)
  call DeallocateArray(inversion%r)
  call DeallocateArray(inversion%s)
  call DeallocateArray(inversion%del_cond)

end subroutine InversionDeallocateWorkArrays

! ************************************************************************** !

subroutine InversionRead(inversion,input,option)
  !
  ! Read inversion options
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/04/21
  
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none
      
  type(inversion_type) :: inversion
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  ! we initialize the word to blanks to avoid error reported by valgrind
  word = ''
  
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (input%ierr /= 0) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','INVERSION')
    call StringToUpper(word)
    select case(trim(word))
    case('MIN_ELECTRICAL_CONDUCTIVITY')
      call InputReadDouble(input,option,inversion%mincond)
      call InputErrorMsg(input,option,'MIN_ELECTRICAL_CONDUCTIVITY', &
                           'INVERSION')
    case('MAX_ELECTRICAL_CONDUCTIVITY')
      call InputReadDouble(input,option,inversion%maxcond)
      call InputErrorMsg(input,option,'MAX_ELECTRICAL_CONDUCTIVITY', &
                           'INVERSION')
    case('MIN_CGLS_ITERATION')
      call InputReadInt(input,option,inversion%miniter)
      call InputErrorMsg(input,option,'MIN_CGLS_ITERATION','INVERSION')
    case('MAX_CGLS_ITERATION')
      call InputReadInt(input,option,inversion%maxiter)
      call InputErrorMsg(input,option,'MAX_CGLS_ITERATION','INVERSION')
    case('BETA')
      call InputReadDouble(input,option,inversion%beta)
      call InputErrorMsg(input,option,'BETA','INVERSION')
    case('BETA_REDUCTION_FACTOR')
      call InputReadDouble(input,option,inversion%beta_red_factor)
      call InputErrorMsg(input,option,'BETA_REDUCTION_FACTOR','INVERSION')
    case('TARGET_CHI2')
      call InputReadDouble(input,option,inversion%target_chi2)
      call InputErrorMsg(input,option,'TARGET_CHI2','INVERSION')
    case('MIN_COST_REDUCTION')
      call InputReadDouble(input,option,inversion%min_phi_red)
      call InputErrorMsg(input,option,'MIN_COST_REDUCTION','INVERSION')
    case('CONSTRAINED_BLOCKS')
      call ConstrainedBlockRead(inversion%constrained_block,input,option)
    case default
      call InputKeywordUnrecognized(input,word,'INVERSION',option)
    end select
  enddo
  call InputPopBlock(input,option)
        
end subroutine InversionRead

! ************************************************************************** !

subroutine ConstrainedBlockRead(constrained_block,input,option)
  !
  ! Read constrained blocks options
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/21/21
  
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(constrained_block_type) :: constrained_block
  type(input_type), pointer :: input
  type(option_type) :: option  

  type(constrained_block_par_type), pointer :: cur_constrained_block
  type(constrained_block_par_type), pointer :: prev_constrained_block

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
    call InputErrorMsg(input,option,'keyword','INVERSION,CONSTRAINED_BLOCKS')
   
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
  ! Date: 05/21/21
  
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(constrained_block_par_type) :: constrained_block
  type(input_type), pointer :: input
  type(option_type) :: option  

  PetscInt :: i,num_block_link
  character(len=MAXWORDLENGTH) :: word

  word = ''

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','INVERSION,CONSTRAINED_BLOCKS')
    call StringToUpper(word)
    select case(trim(word))
    case('STRUCTURE_METRIC')
      call InputReadInt(input,option,constrained_block%structure_metric)
      call InputErrorMsg(input,option,'STRUCTURE_METRIC', &
                         'INVERSION,CONSTRAINED_BLOCKS')
    case('WEIGHING_FUNCTION')
      call InputReadInt(input,option,constrained_block%weighing_function)
      call InputErrorMsg(input,option,'WEIGHING_FUNCTION', &
                         'INVERSION,CONSTRAINED_BLOCKS')
    case('WEIGHING_FUNCTION_MEAN')
      call InputReadDouble(input,option,constrained_block%weighing_function_mean)
      call InputErrorMsg(input,option,'WEIGHING_FUNCTION_MEAN', &
                         'INVERSION,CONSTRAINED_BLOCKS')
    case('WEIGHING_FUNCTION_STD_DEVIATION')
      call InputReadDouble(input,option,constrained_block%weighing_function_std_dev)
      call InputErrorMsg(input,option,'WEIGHING_FUNCTION_STD_DEVIATION', &
                         'INVERSION,CONSTRAINED_BLOCKS')
    case('BLOCK_LINKS')
      call InputReadInt(input,option,num_block_link)
      call InputErrorMsg(input,option,'BLOCK_LINKS', &
                         'INVERSION,CONSTRAINED_BLOCKS')
      constrained_block%num_block_link = num_block_link
      allocate(constrained_block%block_link(num_block_link))
      do i=1,num_block_link                   
        call InputReadCard(input,option,constrained_block%block_link(i))
        call InputErrorMsg(input,option,'BLOCK_LINKS', &
                         'INVERSION,CONSTRAINED_BLOCKS')
      enddo
    case('ANISOTROPIC_WEIGHTS')
      do i=1,THREE_INTEGER
        call InputReadDouble(input,option,constrained_block%aniso_weight(i))
        call InputErrorMsg(input,option,'ANISOTROPY_WEIGHTS', &
                         'INVERSION,CONSTRAINED_BLOCKS')
      enddo
    case('RELATIVE_WEIGHT')
      call InputReadDouble(input,option,constrained_block%relative_weight)
      call InputErrorMsg(input,option,'RELATIVE_WEIGHT', &
                         'INVERSION,CONSTRAINED_BLOCKS')
    case('REFERENCE_CONDUCTIVITY')
      call InputReadDouble(input,option,constrained_block%reference_conductivity)
      call InputErrorMsg(input,option,'REFERENCE_CONDUCTIVITY', &
                         'INVERSION,CONSTRAINED_BLOCKS')
    case default
      call InputKeywordUnrecognized(input,word, &
                                    'INVERSION,CONSTRAINED_BLOCKS',option)
    end select   
  enddo
  call InputPopBlock(input,option)

end subroutine ConstrainedBlockParRead

! ************************************************************************** !

subroutine InversionCheckConvergence(inversion,survey)
  !
  ! Check Inversion convergence
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/05/21
      
  implicit none
  
  type(inversion_type) :: inversion
  type(survey_type) :: survey

  inversion%converg_flag = PETSC_FALSE
  call InversionEvaluateCostFunctions(inversion,survey)
  if (inversion%current_chi2 <= inversion%target_chi2) &
                           inversion%converg_flag = PETSC_TRUE
      
end subroutine InversionCheckConvergence

! ************************************************************************** !

subroutine InversionEvaluateCostFunctions(inversion,survey)
  !
  ! Evaluates cost functions for inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/05/21
        
  implicit none
        
  type(inversion_type) :: inversion
  type(survey_type) :: survey
  
  PetscInt :: idata,ndata,ncull
  PetscReal :: err_mean,err_sdev
  PetscReal, allocatable :: data_vector(:)

  ndata = survey%num_measurement
  allocate(data_vector(ndata))
  data_vector = 0.d0

  data_vector = survey%Wd * (survey%dobs - survey%dsim)
    
  ncull = 0
  if (inversion%cull_flag) then
    err_mean = sum(data_vector) / ndata
    err_sdev = sqrt( dot_product( (data_vector - err_mean), &
                                  (data_vector - err_mean) ) / ndata )
    do idata=1,ndata
      if (data_vector(idata) < (err_mean - err_sdev*inversion%cull_dev) .or. &
          data_vector(idata) > (err_mean + err_sdev*inversion%cull_dev)) then
        survey%Wd_cull(idata) = 0
        ncull = ncull + 1
      endif
    enddo    
  endif

  data_vector = survey%Wd_cull * data_vector
    
  inversion%phi_data = dot_product(data_vector,data_vector)
  inversion%current_chi2 = inversion%phi_data / (ndata - ncull)
    
  deallocate(data_vector)

  ! TODO: compute phi_model
  inversion%phi_total = inversion%phi_data + inversion%phi_model

  if (inversion%iteration == 1) then
    inversion%phi_data_0 = inversion%phi_data
    inversion%phi_model_0 = inversion%phi_model
    inversion%phi_total_0 = inversion%phi_total
  endif  

end subroutine InversionEvaluateCostFunctions
  
! ************************************************************************** !

subroutine InversionCheckBeta(inversion)
  !
  ! Check Beta if it needs cooling/reduction
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/05/21
        
  implicit none
        
  type(inversion_type) :: inversion

  if ( abs(inversion%phi_total_0 - inversion%phi_total) <= &
       inversion%min_phi_red ) then
    inversion%beta = inversion%beta * inversion%beta_red_factor
    inversion%phi_model = inversion%beta_red_factor * inversion%phi_model
  endif

  ! update the cost functions
  inversion%phi_data_0 = inversion%phi_data
  inversion%phi_model_0 = inversion%phi_model
  inversion%phi_total_0 = inversion%phi_total

end subroutine InversionCheckBeta

! ************************************************************************** !

subroutine ConstrainedBlockParDestroy(constrained_block)
  ! 
  ! ConstrainedBlockParDestroy: Deallocates a constrained block par object
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 05/20/21
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
  ! Date: 05/20/21
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

subroutine InversionDestroy(inversion)
  !
  ! Deallocates inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 05/03/21
  !

  use Utility_module, only : DeallocateArray
      
  implicit none
      
  type(inversion_type), pointer :: inversion
      
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

  deallocate(inversion)
  nullify(inversion)

end subroutine InversionDestroy

end module Inversion_module
