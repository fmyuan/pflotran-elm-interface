module Dataset_Global_Hdf5_list_class
! A generalized data_mediator_dataset class, i.e. not specifically for MASS TRANSFER

#include "petsc/finclude/petscvec.h"
  use petscvec

  use PFLOTRAN_Constants_module
  use Data_Mediator_Base_class
  use Dataset_Global_HDF5_class
  
  implicit none

  private

  type, public, extends(data_mediator_base_type) :: dataset_global_hdf5_list_type
    PetscInt :: idx  ! index (order) in the list
    class(dataset_global_hdf5_type), pointer :: dataset
  contains
    procedure, public :: Update => Update
    procedure, public :: Strip => DatasetGlobalHd5ListStrip
  end type data_mediator_dataset_type
  
  public :: DatasetGlobalHd5ListCreate, &
            DatasetGlobalHd5ListRead, &
            DatasetGlobalHd5ListInit

contains

! ************************************************************************** !

function DatasetGlobalHd5ListCreate()
  ! 
  ! Creates a data mediator object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 
  
  implicit none

  class(data_mediator_dataset_type), pointer :: DatasetGlobalHd5ListCreate
  
  class(data_mediator_dataset_type), pointer :: data_mediator
  
  allocate(data_mediator)
  call DataMediatorBaseCreate(data_mediator)
  data_mediator%idof = 0
  nullify(data_mediator%dataset)
  DatasetGlobalHd5ListCreate => data_mediator

end function DatasetGlobalHd5ListCreate

! ************************************************************************** !

subroutine DatasetGlobalHd5ListRead(data_mediator,input,option)
  ! 
  ! Reads in contents of a data mediator card
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(data_mediator_dataset_type) :: data_mediator
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word

  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','HDF5_LIST')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('IDX')
        call InputReadInt(input,option,data_mediator%idx)
        call InputErrorMsg(input,option,'idx','HDF5_LIST')
      case('DATASET')
        data_mediator%dataset => DatasetGlobalHDF5Create()
        call InputReadNChars(input,option, &
                             data_mediator%dataset%name,&
                             MAXWORDLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'DATASET,NAME','HDF5_LIST')
      case default
        call InputKeywordUnrecognized(input,keyword,'HDF5_LIST',option)
    end select
    
  enddo  
  call InputPopBlock(input,option)

end subroutine DatasetGlobalHd5ListRead

! ************************************************************************** !

subroutine DatasetGlobalHd5ListInit(data_mediator, discretization, &
                                   available_datasets, option)
  ! 
  ! Initializes data mediator object opening dataset to
  ! set up times, vectors, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/09/13
  ! 
  use Discretization_module
  use Dataset_Base_class
  use Dataset_Common_HDF5_class
  use Option_module
  use String_module

  implicit none
  
  class(data_mediator_dataset_type) :: data_mediator
  type(discretization_type) :: discretization
  class(dataset_base_type), pointer :: available_datasets
  type(option_type) :: option
  
  class(dataset_base_type), pointer :: dataset_base_ptr
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  if (.not.associated(data_mediator%dataset)) then
    option%io_buffer = 'A "global" DATASET does not exist for ' // &
      'MASS_TRANSFER object "' // trim(data_mediator%name) // '".'
    call PrintErrMsg(option)
  endif

  string = 'Data Mediator ' // trim(data_mediator%name)
  dataset_base_ptr => &
    DatasetBaseGetPointer(available_datasets,data_mediator%dataset%name, &
                          string,option)
  call DatasetGlobalHDF5Destroy(data_mediator%dataset)
  select type(dataset => dataset_base_ptr)
    class is(dataset_global_hdf5_type)
      data_mediator%dataset => dataset
    class default
      option%io_buffer = 'DATASET ' // trim(dataset%name) // 'is not of ' // &
        'GLOBAL type, which is necessary for all HDF5_LIST objects.'
      call PrintErrMsg(option)
  end select
  ! dm_wrapper is solely a pointer; it should not be allocated
  data_mediator%dataset%dm_wrapper => discretization%dm_1dof
  data_mediator%dataset%local_size = discretization%grid%nlmax
  data_mediator%dataset%global_size = discretization%grid%nmax


#ifdef ELM_PFLOTRAN
  ! by-passing reading dataset%time_storage from the dummy hdf5 file ('dummy.h5'),
  ! which implies that CLM directly updates 'this%dataset%rarray' at each clm-timestep
  if (.not.(StringEndswith(trim(data_mediator%dataset%filename), "dummy.h5")) ) then
      option%io_buffer = 'DATASET filename ' // trim(data_mediator%dataset%filename) // &
        'should be named as dummy.h5 (none) for coupling with CLM '
      call printErrMsg(option)
  endif
#else
  
  if (.not.associated(data_mediator%dataset%time_storage)) then
    call DatasetCommonHDF5ReadTimes(data_mediator%dataset%filename, &
                                    data_mediator%dataset%hdf5_dataset_name, &
                                    data_mediator%dataset%time_storage,option)
    ! if time interpolation methods not set in hdf5 file, set to default of STEP
    if (data_mediator%dataset%time_storage%time_interpolation_method == &
        INTERPOLATION_NULL) then
      data_mediator%dataset%time_storage%time_interpolation_method = &
        INTERPOLATION_STEP
    endif
  endif 

#endif
  
end subroutine DatasetGlobalHd5ListInit

! ************************************************************************** !

recursive subroutine DatasetGlobalHd5ListUpdate(this,data_mediator_vec,option)
  ! 
  ! Updates a data mediator object transfering data from
  ! the buffer into the PETSc Vec
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 
  use Option_module
  use String_module
  
  implicit none
  
  class(data_mediator_dataset_type) :: this
  Vec :: data_mediator_vec
  type(option_type) :: option  
  
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: ndof_per_cell
  PetscInt :: mdof_local_size
  PetscInt :: offset
  PetscInt :: i
  PetscErrorCode :: ierr
  
#ifdef ELM_PFLOTRAN
  ! by-passing reading dataset from the dummy hdf5 file,
  ! which implies that CLM directly updates 'this%dataset%rarray' instead of reading
  if (.not.(StringEndswith(trim(this%dataset%filename), "dummy.h5")) ) then
    option%io_buffer = 'DATASET filename ' // trim(this%dataset%filename) // &
      'must be named as dummy.h5 (none) for coupling with CLM '
    call printErrMsg(option)

  else
    ! need to initialize this%dataset%rarray, when first-time uses it (and skip reading from h5 file)
    if (.not.associated(this%dataset%rarray)) then
      if (this%dataset%local_size == 0) then
        option%io_buffer = 'Local size of Global Dataset has not been set.'
        call printErrMsg(option)
      endif
      allocate(this%dataset%rarray(this%dataset%local_size))
      this%dataset%rarray = 0.d0
    endif

  endif
#else

  call DatasetGlobalHDF5Load(this%dataset,option)
#endif

  call VecGetLocalSize(data_mediator_vec,mdof_local_size,ierr);CHKERRQ(ierr)
  ndof_per_cell = mdof_local_size / this%dataset%local_size
  if (mod(mdof_local_size,this%dataset%local_size) > 0) then
    option%io_buffer = 'Mismatched vector size in MassTransferUpdate.'
    call PrintErrMsg(option)
  endif
  call VecGetArrayF90(data_mediator_vec,vec_ptr,ierr);CHKERRQ(ierr)
  offset = this%idx
  do i = 1, this%dataset%local_size
    vec_ptr(offset) = vec_ptr(offset) + this%dataset%rarray(i)
    offset = offset + ndof_per_cell
  enddo
  call VecRestoreArrayF90(data_mediator_vec,vec_ptr,ierr);CHKERRQ(ierr)
  
end subroutine DatasetGlobalHd5ListUpdate

! ************************************************************************** !

recursive subroutine DatasetGlobalHd5ListStrip(this)
  ! 
  ! Destroys a data mediator object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/01/13
  ! 

  implicit none
  
  class(data_mediator_dataset_type) :: this
  
  PetscErrorCode :: ierr
  
  ! update the next one
  if (associated(this%next)) then
    call this%next%Strip()
    deallocate(this%next)
    nullify(this%next)
  endif 
  
  ! Simply nullify the pointer as the dataset resides in a list to be
  ! destroyed separately.
  nullify(this%dataset)

end subroutine DatasetGlobalHd5ListStrip

end module Dataset_Global_Hdf5_List_class
