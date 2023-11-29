module HDF5_Aux_module

#include "petsc/finclude/petscsys.h"
  use, intrinsic :: iso_c_binding
  use petscsys
  use hdf5
  use h5lt

  use Driver_class
  use Logging_module
  use Option_module
  use PFLOTRAN_Constants_module
  use Print_module

  implicit none

  private

  PetscInt, parameter, public :: HDF5_READ_BUFFER_SIZE = 1000000
  PetscMPIInt, parameter :: ON=1, OFF=0
!#define HDF5_BROADCAST

  PetscErrorCode :: ierr

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  interface HDF5AttributeRead
    module procedure :: HDF5AttributeReadBase
    module procedure :: HDF5AttributeReadDouble1
    module procedure :: HDF5AttributeReadDouble2
    module procedure :: HDF5AttributeReadInteger1
    module procedure :: HDF5AttributeReadInteger2
    module procedure :: HDF5AttributeReadIntegerArray1D1
    module procedure :: HDF5AttributeReadIntegerArray1D2
  end interface

  interface HDF5AttributeWrite
    module procedure :: HDF5AttributeWriteBase
    module procedure :: HDF5AttributeWriteDouble1
    module procedure :: HDF5AttributeWriteDouble2
    module procedure :: HDF5AttributeWriteInteger1
    module procedure :: HDF5AttributeWriteInteger2
    module procedure :: HDF5AttributeWriteIntArray1D1
    module procedure :: HDF5AttributeWriteIntArray1D2
  end interface

  interface HDF5DatasetRead
    module procedure :: HDF5DatasetReadInteger1
    module procedure :: HDF5DatasetReadInteger2
    module procedure :: HDF5DatasetReadDoubleArray1D1
    module procedure :: HDF5DatasetReadDoubleArray1D2
  end interface

  interface HDF5DatasetWrite
    module procedure :: HDF5DatasetWriteInteger1
    module procedure :: HDF5DatasetWriteInteger2
    module procedure :: HDF5DatasetWriteDoubleArray1D1
    module procedure :: HDF5DatasetWriteDoubleArray1D2
  end interface

  interface HDF5FileOpen
    module procedure :: HDF5FileOpen1
    module procedure :: HDF5FileOpen2
    module procedure :: HDF5FileOpen3
    module procedure :: HDF5FileOpen4
    module procedure :: HDF5FileOpen5
  end interface

  interface HDF5FileOpenReadOnly
    module procedure :: HDF5FileOpenReadOnly1
    module procedure :: HDF5FileOpenReadOnly2
    module procedure :: HDF5FileOpenReadOnly3
  end interface

  interface HDF5FileClose
    module procedure :: HDF5FileClose1
    module procedure :: HDF5FileClose2
    module procedure :: HDF5FileClose3
  end interface

  interface HDF5GroupCreate
    module procedure :: HDF5GroupCreate1
    module procedure :: HDF5GroupCreate2
    module procedure :: HDF5GroupCreate3
  end interface

  interface HDF5GroupOpen
    module procedure :: HDF5GroupOpen1
    module procedure :: HDF5GroupOpen2
    module procedure :: HDF5GroupOpen3
  end interface

  interface HDF5GroupOpenOrCreate
    module procedure :: HDF5GroupOpenOrCreate1
    module procedure :: HDF5GroupOpenOrCreate2
    module procedure :: HDF5GroupOpenOrCreate3
  end interface

  interface HDF5GroupClose
    module procedure :: HDF5GroupClose1
    module procedure :: HDF5GroupClose2
    module procedure :: HDF5GroupClose3
  end interface

  interface HDF5DatasetOpen
    module procedure :: HDF5DatasetOpen1
    module procedure :: HDF5DatasetOpen2
    module procedure :: HDF5DatasetOpen3
  end interface

  interface HDF5DatasetClose
    module procedure :: HDF5DatasetClose1
    module procedure :: HDF5DatasetClose2
    module procedure :: HDF5DatasetClose3
  end interface

  interface HDF5GroupExists
    module procedure :: HDF5GroupExists1
    module procedure :: HDF5GroupExists2
    module procedure :: HDF5GroupExists3
  end interface

  interface HDF5DatasetExists
    module procedure :: HDF5DatasetExists1
    module procedure :: HDF5DatasetExists2
    module procedure :: HDF5DatasetExists3
  end interface


  public :: HDF5ReadNDimRealArray, &
            HDF5GroupExists, &
            HDF5DatasetExists, &
            HDF5MakeStringCompatible, &
            HDF5ReadDbase, &
            HDF5AttributeRead, &
            HDF5AttributeWrite, &
            HDF5DatasetRead, &
            HDF5DatasetWrite, &
            HDF5FileOpen, &
            HDF5FileTryOpen, &
            HDF5FileOpenReadOnly, &
            HDF5FileClose, &
            HDF5GroupOpen, &
            HDF5GroupCreate, &
            HDF5GroupOpenOrCreate, &
            HDF5GroupClose, &
            HDF5DatasetOpen, &
            HDF5DatasetClose, &
            HDF5Init, &
            HDF5Finalize

contains

! ************************************************************************** !

subroutine HDF5Init()
  !
  ! From the HDF5 library documentation:
  !
  ! When the HDF5 Library is employed in a Fortran90 application, h5open_f
  ! initializes global variables (for example, predefined types) and performs
  ! other tasks required to initialize the HDF5 Fortran Library. h5open_f and
  ! h5close_f are required calls in HDF5 Fortran applications.
  !
  ! It only needs to be called once, but no damage if more than once
  !
  ! Author: Glenn Hammond
  ! Date: 07/06/20
  !

  integer :: hdf5_err

  call h5open_f(hdf5_err)

end subroutine HDF5Init

! ************************************************************************** !

subroutine HDF5ReadNDimRealArray(option,file_id,dataset_name,ndims,dims, &
                                 real_array)
  !
  ! Read in an n-dimensional array from an hdf5 file
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/10
  !
  implicit none

  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  integer(HID_T) :: file_id
  PetscInt :: ndims
  PetscInt, pointer :: dims(:)
  PetscReal, pointer :: real_array(:)

  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: offset(1), length(1), stride(1)
  PetscMPIInt :: rank_mpi
  integer(HSIZE_T) :: num_reals_in_dataset
  PetscInt :: temp_int, i
  integer :: ndims_hdf5
  integer :: hdf5_err

  call PetscLogEventBegin(logging%event_read_ndim_real_array_hdf5, &
                          ierr);CHKERRQ(ierr)

  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space

  call h5sget_simple_extent_ndims_f(file_space_id,ndims_hdf5,hdf5_err)
  ndims = ndims_hdf5
  allocate(dims_h5(ndims))
  allocate(max_dims_h5(ndims))
  allocate(dims(ndims))
  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
  dims = int(dims_h5)
  call h5sget_simple_extent_npoints_f(file_space_id,num_reals_in_dataset, &
                                      hdf5_err)
  temp_int = dims(1)
  do i = 2, ndims
    temp_int = temp_int * dims(i)
  enddo

  rank_mpi = 1
  offset = 0
  length = num_reals_in_dataset
  stride = 1

  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err,length)

  allocate(real_array(num_reals_in_dataset))
  real_array = 0.d0
#ifdef HDF5_BROADCAST
  if (OptionIsIORank(option)) then
#endif
    call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
    call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_array,length, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
#ifdef HDF5_BROADCAST
  endif
  if (option%comm%size > 1) then
    int_mpi = num_reals_in_dataset
    call MPI_Bcast(real_array,int_mpi,MPI_DOUBLE_PRECISION,option%io_rank, &
                   option%mycomm,ierr);CHKERRQ(ierr)
  endif
#endif

  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call HDF5DatasetClose(data_set_id,option)

  deallocate(dims_h5)
  deallocate(max_dims_h5)

  call PetscLogEventEnd(logging%event_read_ndim_real_array_hdf5, &
                        ierr);CHKERRQ(ierr)

end subroutine HDF5ReadNDimRealArray

#if defined(PARALLELIO_LIB)

! ************************************************************************** !

subroutine HDF5ReadDatasetReal1D(filename,dataset_name,read_option,option, &
                                 data,data_dims,dataset_dims)
  !
  ! Author: Gautam Bisht
  ! Date: 05/13/2010
  !

  implicit none

#if defined(PARALLELIO_LIB)
  include "piof.h"
#endif

  ! in
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: dataset_name
  integer                        :: read_option
  type(option_type)              :: option

  ! out
  PetscReal,pointer              :: data(:)
  PetscInt                       :: data_dims(1)
  PetscInt                       :: dataset_dims(1)

  ! local
  integer :: file_id
  integer :: ndims
  PetscInt :: ii, remainder

  PetscErrorCode :: ierr

  ! Open file collectively
  filename = trim(filename) // CHAR(0)
  call parallelIO_open_file(filename, option%ioread_group_id, &
                            FILE_READONLY, file_id, ierr)

  ! Get dataset dimnesions
  call parallelIO_get_dataset_ndims(ndims, file_id, dataset_name, &
                                    option%ioread_group_id, ierr)
  if (ndims.ne.1) then
    option%io_buffer='Dimension of ' // dataset_name // ' dataset in ' // &
      filename // ' is not equal to 1.'
    call PrintErrMsg(option)
  endif

  ! Get size of each dimension
  call parallelIO_get_dataset_dims(dataset_dims, file_id, dataset_name, &
                                   option%ioread_group_id, ierr)

  data_dims(1) = dataset_dims(1)/option%comm%size

  remainder = dataset_dims(1) - data_dims(1)*option%comm%size
  if (option%myrank < remainder) data_dims(1) = data_dims(1) + 1


  allocate(data(data_dims(1)))

  !call parallelIO_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)

  ! Read the dataset collectively
  call parallelIO_read_dataset( data, PIO_DOUBLE, ndims, dataset_dims, &
                               data_dims, file_id, dataset_name, &
                               option%ioread_group_id, &
                               NONUNIFORM_CONTIGUOUS_READ, ierr)

  !data_dims(1) = data_dims(1) + data_dims(2)
  !data_dims(2) = data_dims(1) - data_dims(2)
  !data_dims(1) = data_dims(1) - data_dims(2)

  !dataset_dims(1) = dataset_dims(1) + dataset_dims(2)
  !dataset_dims(2) = dataset_dims(1) - dataset_dims(2)
  !dataset_dims(1) = dataset_dims(1) - dataset_dims(2)

  ! Close file
  call parallelIO_close_file( file_id, option%ioread_group_id, ierr)

end subroutine HDF5ReadDatasetReal1D

#endif
! PARALLELIO_LIB

! ************************************************************************** !

subroutine HDF5GroupOpen1(parent_id,group_name,group_id,option)
  !
  ! Opens an HDF5 group with proper error messaging when not found.
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/23
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: group_name
  integer(HID_T) :: group_id
  type(option_type) :: option

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5GroupOpen(parent_id,group_name,group_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5GroupOpen1

! ************************************************************************** !

subroutine HDF5GroupOpen2(parent_id,group_name,group_id,driver)
  !
  ! Opens an HDF5 group with proper error messaging when not found.
  !
  ! Author: Glenn Hammond
  ! Date: 06/28/18
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: group_name
  integer(HID_T) :: group_id
  class(driver_type) :: driver

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5GroupOpen(parent_id,group_name,group_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5GroupOpen2

! ************************************************************************** !

subroutine HDF5GroupOpen3(parent_id,group_name,group_id,print_handler)
  !
  ! Opens an HDF5 group with proper error messaging when not found.
  !
  ! Author: Glenn Hammond
  ! Date: 06/28/18
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: group_name
  integer(HID_T) :: group_id
  type(print_handler_type) print_handler

  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt :: hdf5_err

  string = trim(group_name)
  call h5gopen_f(parent_id,string,group_id,hdf5_err)
  if (hdf5_err < 0) then
    call PrintErrorMessage(print_handler, &
                            'HDF5 Group "' // trim(string) // '" not found.')
  endif

end subroutine HDF5GroupOpen3

! ************************************************************************** !

subroutine HDF5GroupCreate1(parent_id,group_name,group_id,option)
  !
  ! Creates an HDF5 group with proper error messaging when failing
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: group_name
  integer(HID_T) :: group_id
  type(option_type) :: option

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5GroupCreate(parent_id,group_name,group_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5GroupCreate1

! ************************************************************************** !

subroutine HDF5GroupCreate2(parent_id,group_name,group_id,driver)
  !
  ! Creates an HDF5 group with proper error messaging when failing
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: group_name
  integer(HID_T) :: group_id
  class(driver_type) :: driver

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5GroupCreate(parent_id,group_name,group_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5GroupCreate2

! ************************************************************************** !

subroutine HDF5GroupCreate3(parent_id,group_name,group_id,print_handler)
  !
  ! Creates an HDF5 group with proper error messaging when failing
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: group_name
  integer(HID_T) :: group_id
  type(print_handler_type) :: print_handler

  character(len=MAXSTRINGLENGTH) :: string
  integer :: hdf5_err

  string = trim(group_name)
  call h5gcreate_f(parent_id,string,group_id,hdf5_err, &
                   OBJECT_NAMELEN_DEFAULT_F)
  if (hdf5_err < 0) then
    if (len_trim(string) == 0) then
      string = 'An HDF5 Group could not be created &
                &due to the name being blank.'
    else
      string = 'HDF5 Group "' // trim(string) // &
        '" could not be created.'
    endif
    call PrintErrorMessage(print_handler,string)
  endif

end subroutine HDF5GroupCreate3

! ************************************************************************** !

subroutine HDF5GroupOpenOrCreate1(parent_id,group_name,group_id,option)
  !
  ! Opens an HDF5 group or creates it if it does not exist
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: group_name
  integer(HID_T) :: group_id
  type(option_type) :: option

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5GroupOpenOrCreate(parent_id,group_name,group_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5GroupOpenOrCreate1

! ************************************************************************** !

subroutine HDF5GroupOpenOrCreate2(parent_id,group_name,group_id,driver)
  !
  ! Opens an HDF5 group or creates it if it does not exist
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: group_name
  integer(HID_T) :: group_id
  class(driver_type) :: driver

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5GroupOpenOrCreate(parent_id,group_name,group_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5GroupOpenOrCreate2

! ************************************************************************** !

subroutine HDF5GroupOpenOrCreate3(parent_id,group_name,group_id,print_handler)
  !
  ! Opens an HDF5 group or creates it if it does not exist
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: group_name
  integer(HID_T) :: group_id
  type(print_handler_type) :: print_handler

  character(len=MAXSTRINGLENGTH) :: string
  integer :: hdf5_err, hdf5_err2

  string = trim(group_name)
  call h5eset_auto_f(OFF,hdf5_err2)
  call h5gopen_f(parent_id,string,group_id,hdf5_err)
  call h5eset_auto_f(ON,hdf5_err2)
  if (hdf5_err /= 0) then
    call HDF5GroupCreate(parent_id,group_name,group_id,print_handler)
  endif

end subroutine HDF5GroupOpenOrCreate3

! ************************************************************************** !

function HDF5GroupExists1(filename,group_name,option)
  !
  ! Returns true if a group exists
  !
  ! Author: Glenn Hammond
  ! Date: 03/26/2012
  !
  implicit none

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: group_name
  type(option_type) :: option

  PetscBool :: HDF5GroupExists1

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  HDF5GroupExists1 = HDF5GroupExists(filename,group_name,print_handler)
  call PrintDestroyHandler(print_handler)

end function HDF5GroupExists1

! ************************************************************************** !

function HDF5GroupExists2(filename,group_name,driver)
  !
  ! Returns true if a group exists
  !
  ! Author: Glenn Hammond
  ! Date: 03/26/2012
  !
  implicit none

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: group_name
  class(driver_type) :: driver

  PetscBool :: HDF5GroupExists2

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  HDF5GroupExists2 = HDF5GroupExists(filename,group_name,print_handler)
  call PrintDestroyHandler(print_handler)

end function HDF5GroupExists2

! ************************************************************************** !

function HDF5GroupExists3(filename,group_name,print_handler)
  !
  ! Returns true if a group exists
  !
  ! Author: Glenn Hammond
  ! Date: 03/26/2012
  !
  use Communicator_Aux_module

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: group_name
  type(print_handler_type) :: print_handler

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscBool :: group_exists
  integer :: hdf5_err, hdf5_err2

  PetscBool :: HDF5GroupExists3

  ! open the file
  call HDF5FileOpenReadOnly(filename,file_id,PETSC_TRUE,'', &
                            print_handler,print_handler%comm)

  call PrintMessage(print_handler,'Testing group: ' // trim(group_name))
  ! I turn off error messaging since if the group does not exist, an error
  ! will be printed, but the user does not need to see this.
  call h5eset_auto_f(OFF,hdf5_err2)
  call h5gopen_f(file_id,group_name,grp_id,hdf5_err)
  group_exists = .not.(hdf5_err < 0)
  call h5eset_auto_f(ON,hdf5_err2)

  if (group_exists) then
    HDF5GroupExists3 = PETSC_TRUE
    call HDF5GroupClose(grp_id,print_handler)
    string = 'Group "' // trim(group_name) // '" in HDF5 file "' // &
      trim(filename) // '" found in file.'
  else
    HDF5GroupExists3 = PETSC_FALSE
    string = 'Group "' // trim(group_name) // '" in HDF5 file "' // &
      trim(filename) // '" not found in file.  Therefore, assuming a ' // &
      'cell-indexed dataset.'
  endif
  call PrintMessage(print_handler,string)

  call HDF5FileClose(file_id,print_handler)

end function HDF5GroupExists3

! ************************************************************************** !

subroutine HDF5GroupClose1(group_id,option)
  !
  ! Closes an HDF5 group with proper error messaging when it fails.
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: group_id
  type(option_type) :: option

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5GroupClose(group_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5GroupClose1

! ************************************************************************** !

subroutine HDF5GroupClose2(group_id,driver)
  !
  ! Closes an HDF5 group with proper error messaging when it fails.
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: group_id
  class(driver_type) :: driver

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5GroupClose(group_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5GroupClose2

! ************************************************************************** !

subroutine HDF5GroupClose3(group_id,print_handler)
  !
  ! Closes an HDF5 group with proper error messaging when it fails.
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: group_id
  type(print_handler_type) :: print_handler

  integer :: hdf5_err

  call h5gclose_f(group_id,hdf5_err)
  call HDF5CloseCheckError(hdf5_err,group_id,print_handler)

end subroutine HDF5GroupClose3

! ************************************************************************** !

function HDF5DatasetExists1(filename,group_name,dataset_name,option)
  !
  ! Returns true if a dataset exists
  !
  ! Author: Gautam Bisht
  ! Date: 04/30/2015
  !
  implicit none

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: group_name
  character(len=MAXWORDLENGTH) :: dataset_name
  type(option_type) :: option

  PetscBool :: HDF5DatasetExists1

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  HDF5DatasetExists1 = HDF5DatasetExists(filename,group_name,dataset_name, &
                                         print_handler,option%comm)
  call PrintDestroyHandler(print_handler)

end function HDF5DatasetExists1

! ************************************************************************** !

function HDF5DatasetExists2(filename,group_name,dataset_name,driver)
  !
  ! Returns true if a dataset exists
  !
  ! Author: Gautam Bisht
  ! Date: 04/30/2015
  !
  implicit none

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: group_name
  character(len=MAXWORDLENGTH) :: dataset_name
  class(driver_type) :: driver

  PetscBool :: HDF5DatasetExists2

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  HDF5DatasetExists2 = HDF5DatasetExists(filename,group_name,dataset_name, &
                                         print_handler,driver%comm)
  call PrintDestroyHandler(print_handler)

end function HDF5DatasetExists2

! ************************************************************************** !

function HDF5DatasetExists3(filename,group_name,dataset_name, &
                            print_handler,comm)
  !
  ! Returns true if a dataset exists
  !
  ! Author: Gautam Bisht
  ! Date: 04/30/2015
  !
  use Communicator_Aux_module

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: group_name
  character(len=MAXWORDLENGTH) :: dataset_name
  type(print_handler_type) :: print_handler
  type(comm_type) :: comm

  character(len=MAXWORDLENGTH) :: group_name_local
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: dataset_id
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscBool :: group_exists
  PetscBool :: dataset_exists
  integer :: hdf5_err, hdf5_err2

  PetscBool :: HDF5DatasetExists3

  if (len_trim(group_name) == 0) then
    group_name_local = "/" // CHAR(0)
  else
    group_name_local = trim(group_name) // "/" // CHAR(0)
  endif

  call HDF5FileOpenReadOnly(filename,file_id,PETSC_TRUE,'', &
                            print_handler,comm)

  ! I turn off error messaging since if the group does not exist, an error
  ! will be printed, but the user does not need to see this.
  call h5eset_auto_f(OFF,hdf5_err2)
  call h5gopen_f(file_id,group_name_local,grp_id,hdf5_err)
  group_exists = .not.(hdf5_err < 0)
  call h5eset_auto_f(ON,hdf5_err2)

  if (.not.group_exists) then
    HDF5DatasetExists3 = PETSC_FALSE
  endif

  call HDF5DatasetOpen(grp_id,dataset_name,dataset_id,print_handler)
  dataset_exists = .not.(hdf5_err < 0)

  if (.not.dataset_exists) then
    HDF5DatasetExists3 = PETSC_FALSE
  else
    HDF5DatasetExists3 = PETSC_TRUE
    call HDF5DatasetClose(dataset_id,print_handler)
  endif

  if (group_exists) call HDF5GroupClose(grp_id,print_handler)

  call HDF5FileClose(file_id,print_handler)

end function HDF5DatasetExists3

! ************************************************************************** !

subroutine HDF5DatasetOpen1(parent_id,dataset_name,dataset_id,option)
  !
  ! Opens an HDF5 dataset
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: dataset_name
  integer(HID_T) :: dataset_id
  type(option_type) :: option

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5DatasetOpen(parent_id,dataset_name,dataset_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetOpen1

! ************************************************************************** !

subroutine HDF5DatasetOpen2(parent_id,dataset_name,dataset_id,driver)
  !
  ! Opens an HDF5 dataset
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: dataset_name
  integer(HID_T) :: dataset_id
  class(driver_type) :: driver

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5DatasetOpen(parent_id,dataset_name,dataset_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetOpen2

! ************************************************************************** !

subroutine HDF5DatasetOpen3(parent_id,dataset_name,dataset_id,print_handler)
  !
  ! Opens an HDF5 dataset
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: parent_id
  character(len=*) :: dataset_name
  integer(HID_T) :: dataset_id
  type(print_handler_type) :: print_handler

  character(len=MAXSTRINGLENGTH) :: string
  integer :: hdf5_err

  string = trim(dataset_name)
  call h5dopen_f(parent_id,string,dataset_id,hdf5_err)
  if (hdf5_err < 0) then
    call PrintErrorMessage(print_handler, &
      'HDF5 Dataset ' // trim(dataset_name) // &
      ' could not be opened within ' // &
      trim(HDF5ObjectGetNameTypeString(parent_id,print_handler))  // '.')
  endif

end subroutine HDF5DatasetOpen3

! ************************************************************************** !

subroutine HDF5DatasetClose1(dataset_id,option)
  !
  ! Closes an HDF5 dataset
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: dataset_id
  type(option_type) :: option

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5DatasetClose3(dataset_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetClose1

! ************************************************************************** !

subroutine HDF5DatasetClose2(dataset_id,driver)
  !
  ! Closes an HDF5 dataset
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: dataset_id
  class(driver_type) :: driver

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5DatasetClose3(dataset_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetClose2

! ************************************************************************** !

subroutine HDF5DatasetClose3(dataset_id,print_handler)
  !
  ! Closes an HDF5 dataset
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  implicit none

  integer(HID_T) :: dataset_id
  type(print_handler_type) :: print_handler

  integer :: hdf5_err

  call h5dclose_f(dataset_id,hdf5_err)
  call HDF5CloseCheckError(hdf5_err,dataset_id,print_handler)

end subroutine HDF5DatasetClose3

! ************************************************************************** !

subroutine HDF5MakeStringCompatible(name)
  !
  ! Replaces '/' in string with '_' for hdf5 names
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

  implicit none

  character(len=*) :: name

  PetscInt :: len, ichar

  len = len_trim(name)
  do ichar = 1, len
    if (name(ichar:ichar) == '/') then
      name(ichar:ichar) = '_'
    endif
  enddo

  name = trim(name)

end subroutine HDF5MakeStringCompatible

! ************************************************************************** !

subroutine HDF5ReadDbase(filename,option)
  !
  ! Read in an ASCII database
  !
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  !
  use String_module
  use Input_Aux_module, only : dbase
  use h5lt

  implicit none

  character(len=*) :: filename
  type(option_type) :: option

  character(len=MAXWORDLENGTH), allocatable :: wbuffer(:)
  character(len=MAXWORDLENGTH) :: wbuffer_word
  PetscReal, allocatable :: rbuffer(:)
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, allocatable :: ibuffer(:)
  PetscInt :: value_index
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: object_name
  character(len=MAXWORDLENGTH) :: word
  integer(HID_T) :: file_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: datatype_id
  integer(HID_T) :: datatype_id2
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HSIZE_T) :: num_values_in_dataset
  integer(SIZE_T) size_t_int
  integer(HSIZE_T) :: offset(1), length(1), stride(1)
  PetscMPIInt :: rank_mpi
  integer :: num_objects
  integer :: i_object
  integer :: object_type
  integer :: class_id
  integer :: hdf5_err
  PetscInt :: num_ints
  PetscInt :: num_reals
  PetscInt :: num_words
  PetscInt :: temp_int

  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call PrintMsg(option)
  call HDF5FileOpenReadOnly(filename,file_id,PETSC_TRUE,'',option)
  call h5gn_members_f(file_id, '.',num_objects, hdf5_err)
  num_ints = 0
  num_reals = 0
  num_words = 0
  ! index is zero-based
  do i_object = 0, num_objects-1
    ! read in to string in case the name is too large.
    call h5gget_obj_info_idx_f(file_id,'.',i_object,string, &
                               object_type,hdf5_err)
    if (len_trim(string) > MAXWORDLENGTH) then
      option%io_buffer = 'HDF5 DBASE object names must be shorter than &
        &32 characters: ' // trim(string)
      call PrintErrMsg(option)
    endif
    object_name = trim(string)
    if (object_type == H5G_DATASET_F) then
      call h5dopen_f(file_id,object_name,dataset_id,hdf5_err)
      call h5dget_type_f(dataset_id, datatype_id, hdf5_err)
      call h5tget_class_f(datatype_id, class_id, hdf5_err)
      ! cannot use a select case statement since the H5T definitions are not
      ! guaranteed to be constant.  the preprocessor throws an error
      if (class_id == H5T_INTEGER_F) then
        num_ints = num_ints + 1
      else if (class_id == H5T_FLOAT_F) then
        num_reals = num_reals + 1
      else if (class_id == H5T_STRING_F) then
        num_words = num_words + 1
      else
        option%io_buffer = 'Unrecognized HDF5 datatype in Dbase: ' // &
          trim(object_name)
        call PrintErrMsg(option)
      endif
      call h5tclose_f(datatype_id, hdf5_err)
      call HDF5DatasetClose(dataset_id,option)
    endif
  enddo
  allocate(dbase)
  nullify(dbase%icard)
  nullify(dbase%rcard)
  nullify(dbase%ccard)
  nullify(dbase%ivalue)
  nullify(dbase%rvalue)
  nullify(dbase%cvalue)
  if (num_ints > 0) then
    allocate(dbase%icard(num_ints))
    dbase%icard = ''
    allocate(dbase%ivalue(num_ints))
    dbase%ivalue = UNINITIALIZED_INTEGER
  endif
  if (num_reals > 0) then
    allocate(dbase%rcard(num_reals))
    dbase%rcard = ''
    allocate(dbase%rvalue(num_reals))
    dbase%rvalue = UNINITIALIZED_DOUBLE
  endif
  if (num_words > 0) then
    allocate(dbase%ccard(num_words))
    dbase%ccard = ''
    allocate(dbase%cvalue(num_words))
    dbase%cvalue = '-999'
  endif
  value_index = 1
  if (option%id > 0) then
    value_index = option%id
  endif
  num_ints = 0
  num_reals = 0
  num_words = 0
  do i_object = 0, num_objects-1
    call h5gget_obj_info_idx_f(file_id,'.',i_object,object_name, &
                               object_type,hdf5_err)
    if (object_type == H5G_DATASET_F) then
! use once HDF5 lite is linked in PETSc
!      call h5ltget_dataset_info_f(file_id,object_name,dims,dummy_int, &
!                                  type_size,hdf5_err)
!      allocate(buffer(dims(1)))
!      buffer = 0.d0
!      call h5ltread_dataset_double_f(file_id,object_name,buffer, &
!                                     dims,hdf5_err)
!      dbase%card(icount) = trim(object_name)
!      if (option%id > 0) then
!        if (option%id > dims(1)) then
!          write(word,*) dims(1)
!          option%io_buffer = 'DBASE dataset "' // trim(object_name) // &
!            '" is too small (' // trim(adjustl(word)) // &
!            ') for number of realizations.'
!          call PrintErrMsg(option)
!        endif
!        dbase%value(icount) = buffer(option%id)
!      else
!        dbase%value(icount) = buffer(1)
!      endif
!      deallocate(buffer)

      call h5dopen_f(file_id,object_name,dataset_id,hdf5_err)
      call h5dget_space_f(dataset_id,file_space_id,hdf5_err)
      ! should be a rank=1 data space
      call h5sget_simple_extent_npoints_f(file_space_id, &
                                          num_values_in_dataset,hdf5_err)
      if (option%id > 0) then
        if (option%id > num_values_in_dataset) then
          temp_int = int(num_values_in_dataset)
          option%io_buffer = 'Data in DBASE_FILENAME "' // &
            trim(filename) // '" object "' // &
            trim(object_name) // &
            '" is too small (' // StringWrite(temp_int) // &
            ') for number of realizations.'
          call PrintErrMsg(option)
        endif
      else if (value_index > num_values_in_dataset) then
        write(word,*) num_values_in_dataset
        option%io_buffer = 'Data in DBASE_FILENAME "' // &
          trim(filename) // '" object "' // &
          trim(object_name) // &
          '" must have at least one value in each dataset.'
        call PrintErrMsg(option)
      endif
      rank_mpi = 1
      offset = 0
      length = num_values_in_dataset
      stride = 1
      call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
      call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
      call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err, &
                              length)

      call h5dget_type_f(dataset_id, datatype_id, hdf5_err)
      call h5tget_class_f(datatype_id, class_id, hdf5_err)
      call h5tclose_f(datatype_id, hdf5_err)
#ifdef HDF5_BROADCAST
      if (OptionIsIORank(option)) then
#endif
      call PetscLogEventBegin(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
      if (class_id == H5T_INTEGER_F) then
        allocate(ibuffer(num_values_in_dataset))
        ibuffer = 0
        call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,ibuffer,length, &
                       hdf5_err,memory_space_id,file_space_id,prop_id)
      else if (class_id == H5T_FLOAT_F) then
        allocate(rbuffer(num_values_in_dataset))
        rbuffer = 0.d0
        call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,rbuffer,length, &
                       hdf5_err,memory_space_id,file_space_id,prop_id)
      else if (class_id == H5T_STRING_F) then
        call h5tcopy_f(H5T_NATIVE_CHARACTER,datatype_id2,hdf5_err)
        size_t_int = MAXWORDLENGTH
        call h5tset_size_f(datatype_id2,size_t_int,hdf5_err)
        allocate(wbuffer(num_values_in_dataset))
        wbuffer = ''
        call h5dread_f(dataset_id,datatype_id2,wbuffer,length, &
                       hdf5_err,memory_space_id,file_space_id,prop_id)
        wbuffer_word = wbuffer(value_index)
        deallocate(wbuffer)
        call h5tclose_f(datatype_id2,hdf5_err)
      endif
      call PetscLogEventEnd(logging%event_h5dread_f,ierr);CHKERRQ(ierr)
#ifdef HDF5_BROADCAST
      endif
      if (option%comm%size > 1) then
        int_mpi = num_values_in_dataset
        if (class_id == H5T_INTEGER_F) then
          call MPI_Bcast(ibuffer,int_mpi,MPI_INTEGER,option%io_rank, &
                         option%mycomm,ierr);CHKERRQ(ierr)
        else if (class_id == H5T_FLOAT_F) then
          call MPI_Bcast(rbuffer,int_mpi,MPI_DOUBLE_PRECISION,option%io_rank, &
                         option%mycomm,ierr);CHKERRQ(ierr)
        else if (class_id == H5T_STRING_F) then
          int_mpi = MAXWORDLENGTH
          call MPI_Bcast(wbuffer_word,int_mpi,MPI_CHARACTER,option%io_rank, &
                         option%mycomm,ierr);CHKERRQ(ierr)
        endif
      endif
#endif
      call h5pclose_f(prop_id,hdf5_err)
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      call h5sclose_f(file_space_id,hdf5_err)
      call HDF5DatasetClose(dataset_id,option)
      call StringToUpper(object_name)
      ! these conditionals must come after the bcasts above!!!
      if (class_id == H5T_INTEGER_F) then
        num_ints = num_ints + 1
        dbase%icard(num_ints) = trim(object_name)
        dbase%ivalue(num_ints) = ibuffer(value_index)
        deallocate(ibuffer)
      else if (class_id == H5T_FLOAT_F) then
        num_reals = num_reals + 1
        dbase%rcard(num_reals) = trim(object_name)
        dbase%rvalue(num_reals) = rbuffer(value_index)
        deallocate(rbuffer)
      else if (class_id == H5T_STRING_F) then
        num_words = num_words + 1
        dbase%ccard(num_words) = trim(object_name)
        dbase%cvalue(num_words) = wbuffer_word
      endif
    endif
  enddo
  call HDF5FileClose(file_id,option)

end subroutine HDF5ReadDbase

! ************************************************************************** !

  subroutine HDF5AttributeReadDouble1(parent_id,attr_type,attr_name, &
                                      buf,option)
  !
  ! Reads a single PetscReal dataset from an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscReal, target, intent(out) :: buf
  type(option_type), intent(in)  :: option

  type(print_handler_type), pointer :: print_handler
  type(c_ptr) :: ptr

  print_handler => OptionCreatePrintHandler(option)
  ptr = c_loc(buf)
  call HDF5AttributeReadBase(parent_id,attr_type,attr_name, &
                             shape(buf),ptr,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeReadDouble1

! ************************************************************************** !

  subroutine HDF5AttributeReadDouble2(parent_id,attr_type,attr_name, &
                                      buf,driver)
  !
  ! Reads a single PetscReal dataset from an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscReal, target, intent(out) :: buf
  class(driver_type), intent(in) :: driver

  type(print_handler_type), pointer :: print_handler
  type(c_ptr) :: ptr

  print_handler => driver%CreatePrintHandler()
  ptr = c_loc(buf)
  call HDF5AttributeReadBase(parent_id,attr_type,attr_name, &
                             shape(buf),ptr,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeReadDouble2

! ************************************************************************** !

subroutine HDF5AttributeReadInteger1(parent_id,attr_type,attr_name, &
                                     buf,option)
  !
  ! Reads a single PetscInt dataset from an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscInt, target, intent(out) :: buf
  type(option_type), intent(in)  :: option

  type(print_handler_type), pointer :: print_handler
  type(c_ptr) :: ptr

  print_handler => OptionCreatePrintHandler(option)
  ptr = c_loc(buf)
  call HDF5AttributeReadBase(parent_id,attr_type,attr_name, &
                             shape(buf),ptr,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeReadInteger1

! ************************************************************************** !

subroutine HDF5AttributeReadInteger2(parent_id,attr_type,attr_name, &
                                     buf,driver)
  !
  ! Reads a single PetscInt dataset from an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscInt, target, intent(out) :: buf
  class(driver_type), intent(in) :: driver

  type(print_handler_type), pointer :: print_handler
  type(c_ptr) :: ptr

  print_handler => driver%CreatePrintHandler()
  ptr = c_loc(buf)
  call HDF5AttributeReadBase(parent_id,attr_type,attr_name, &
                             shape(buf),ptr,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeReadInteger2

! ************************************************************************** !

subroutine HDF5AttributeReadIntegerArray1D1(parent_id,attr_type,attr_name, &
                                            buf,option)
  !
  ! Reads a 1D PetscInt array dataset from an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscInt, target, intent(inout) :: buf(:)
  type(option_type), intent(in)  :: option

  type(print_handler_type), pointer :: print_handler
  type(c_ptr) :: ptr

  print_handler => OptionCreatePrintHandler(option)
  ptr = c_loc(buf)
  call HDF5AttributeReadBase(parent_id,attr_type,attr_name, &
                             shape(buf),ptr,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeReadIntegerArray1D1

! ************************************************************************** !

subroutine HDF5AttributeReadIntegerArray1D2(parent_id,attr_type,attr_name, &
                                            buf,driver)
  !
  ! Reads a 1D PetscInt array dataset from an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscInt, target, intent(inout) :: buf(:)
  class(driver_type), intent(in) :: driver

  type(print_handler_type), pointer :: print_handler
  type(c_ptr) :: ptr

  print_handler => driver%CreatePrintHandler()
  ptr = c_loc(buf)
  call HDF5AttributeReadBase(parent_id,attr_type,attr_name, &
                             shape(buf),ptr,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeReadIntegerArray1D2

! ************************************************************************** !

subroutine HDF5AttributeReadBase(parent_id,attr_type,attr_name, &
                                 attr_shape,buf,print_handler)
  !
  ! Reads a dataset from an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  use String_module

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscInt, intent(in) :: attr_shape(:)
  type(c_ptr), intent(inout) :: buf
  type(print_handler_type), intent(in)  :: print_handler

  integer(HID_T) :: attr_id
  integer(HID_T) :: attr_type2
  integer(HSIZE_T) :: dims(10)
  integer(HSIZE_T) :: size_
  integer(HSIZE_T) :: size2
  integer :: ndims
  integer :: i
  integer :: hdf5_err
  PetscBool :: is_equal

  dims = 0
  call h5aopen_f(parent_id,attr_name,attr_id,hdf5_err)
  if (hdf5_err /= 0) then
    call PrintErrorMessage(print_handler, &
                           'Opening attribute "' // trim(attr_name) // '".')
  endif
  ndims = size(attr_shape)
  call h5aget_type_f(attr_id,attr_type2,hdf5_err)
  call h5tequal_f(attr_type,attr_type2,is_equal,hdf5_err)
  if (.not.is_equal) then
    call PrintErrorMessage(print_handler, &
                           'Mismatch in attribute type for "' // &
                           trim(attr_name) // '" read.')
  endif
  call h5aget_storage_size_f(attr_id,size_,hdf5_err)
  call h5tget_size_f(attr_type,size2,hdf5_err)
  do i = 1, ndims
    size2 = size2 * attr_shape(i)
  enddo
  if (size_ /= size2) then
    call PrintErrorMessage(print_handler, &
                           'Mismatch in attribute size for "' // &
                           trim(attr_name) // '" read: ' // &
                           trim(StringWrite(int(size_))) // ' vs ' // &
                           trim(StringWrite(int(size2))))
  endif
  call h5aread_f(attr_id,attr_type,buf,hdf5_err)
  if (hdf5_err /= 0) then
    call PrintErrorMessage(print_handler, &
                           'Error reading attribute "' // &
                           trim(attr_name) // '".')
  endif
  call h5aclose_f(attr_id,hdf5_err)

end subroutine HDF5AttributeReadBase

! ************************************************************************** !

subroutine HDF5AttributeWriteDouble1(parent_id,attr_type,attr_name, &
                                     buf,option)
  !
  ! Writes a single PetscInt dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscReal, target, intent(in) :: buf
  type(option_type), intent(in)  :: option

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5AttributeWrite(parent_id,attr_type,attr_name, &
                          shape(buf),c_loc(buf),print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeWriteDouble1

! ************************************************************************** !

subroutine HDF5AttributeWriteDouble2(parent_id,attr_type,attr_name, &
                                     buf,driver)
  !
  ! Writes a single PetscInt dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscReal, target, intent(in) :: buf
  class(driver_type), intent(in) :: driver

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5AttributeWrite(parent_id,attr_type,attr_name, &
                          shape(buf),c_loc(buf),print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeWriteDouble2

! ************************************************************************** !

subroutine HDF5AttributeWriteInteger1(parent_id,attr_type,attr_name, &
                                      buf,option)
  !
  ! Writes a single PetscInt dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscInt, target, intent(in) :: buf
  type(option_type), intent(in)  :: option

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5AttributeWrite(parent_id,attr_type,attr_name, &
                          shape(buf),c_loc(buf),print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeWriteInteger1

! ************************************************************************** !

subroutine HDF5AttributeWriteInteger2(parent_id,attr_type,attr_name, &
                                      buf,driver)
  !
  ! Writes a single PetscInt dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscInt, target, intent(in) :: buf
  class(driver_type), intent(in) :: driver

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5AttributeWrite(parent_id,attr_type,attr_name, &
                          shape(buf),c_loc(buf),print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeWriteInteger2

! ************************************************************************** !

subroutine HDF5AttributeWriteIntArray1D1(parent_id,attr_type,attr_name, &
                                         buf,option)
  !
  ! Writes a 1D PetscInt array dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscInt, target, intent(in) :: buf(:)
  type(option_type), intent(in)  :: option

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5AttributeWrite(parent_id,attr_type,attr_name, &
                          shape(buf),c_loc(buf),print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeWriteIntArray1D1

! ************************************************************************** !

subroutine HDF5AttributeWriteIntArray1D2(parent_id,attr_type,attr_name, &
                                         buf,driver)
  !
  ! Writes a 1D PetscInt array dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscInt, target, intent(in) :: buf(:)
  class(driver_type), intent(in) :: driver

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5AttributeWrite(parent_id,attr_type,attr_name, &
                          shape(buf),c_loc(buf),print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5AttributeWriteIntArray1D2

! ************************************************************************** !

subroutine HDF5AttributeWriteBase(parent_id,attr_type,attr_name, &
                                  attr_shape,buf,print_handler)
  !
  ! Writes a dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  use String_module

  integer(HID_T), intent(in) :: parent_id
  integer(HID_T), intent(in) :: attr_type
  character(len=*), intent(in) :: attr_name
  PetscInt, intent(in) :: attr_shape(:)
  type(c_ptr), intent(in) :: buf
  type(print_handler_type), intent(in)  :: print_handler

  integer(HID_T) :: attr_id
  integer(HID_T) :: attr_type2
  integer(HID_T) :: dspace_id
  integer(HSIZE_T) :: dims(10)
  integer(HSIZE_T) :: size_
  integer(HSIZE_T) :: size2
  integer :: ndims
  integer :: i
  integer :: hdf5_err
  integer :: hdf5_err2
  PetscBool :: is_equal

  dims = 0
  call h5eset_auto_f(OFF,hdf5_err2)
  call h5aopen_f(parent_id,attr_name,attr_id,hdf5_err)
  call h5eset_auto_f(ON,hdf5_err2)
  ndims = size(attr_shape)
  if (hdf5_err == 0) then
    ! if the attribute exists, overwrite it
    call h5aget_type_f(attr_id,attr_type2,hdf5_err)
    call h5tequal_f(attr_type,attr_type2,is_equal,hdf5_err)
    if (.not.is_equal) then
      call PrintErrorMessage(print_handler, &
                             'Mismatch in attribute type for "' // &
                             trim(attr_name) // '" overwrite.')
    endif
    call h5aget_storage_size_f(attr_id,size_,hdf5_err)
    call h5tget_size_f(attr_type,size2,hdf5_err)
    do i = 1, ndims
      size2 = size2 * attr_shape(i)
    enddo
    if (size_ /= size2) then
      call PrintErrorMessage(print_handler, &
                             'Mismatch in attribute size for "' // &
                             trim(attr_name) // '" overwrite: ' // &
                             trim(StringWrite(int(size_))) // ' vs ' // &
                             trim(StringWrite(int(size2))))
    endif
  else
    ! otherwise, create it
    if (ndims == 0) then
      call h5screate_f(H5S_SCALAR_F,dspace_id,hdf5_err)
    else
      dims(1:ndims) = attr_shape
      call h5screate_simple_f(ndims,dims,dspace_id,hdf5_err)
    endif
    call h5acreate_f(parent_id,attr_name,attr_type,dspace_id, &
                     attr_id,hdf5_err)
    call h5sclose_f(dspace_id,hdf5_err)
  endif
  call h5awrite_f(attr_id,attr_type,buf,hdf5_err)
  if (hdf5_err /= 0) then
    call PrintErrorMessage(print_handler, &
                           'Error writing attribute "' // &
                           trim(attr_name) // '".')
  endif
  call h5aclose_f(attr_id,hdf5_err)

end subroutine HDF5AttributeWriteBase

! ************************************************************************** !

subroutine HDF5DatasetReadDoubleArray1D1(loc_id,dset_name,buf,option)
  !
  ! Writes a 1D PetscReal dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: dset_name  ! must be of variable length
  PetscReal, target, intent(inout) :: buf(:)
  type(option_type), intent(inout) :: option

  type(print_handler_type), pointer :: print_handler
  integer(HSIZE_T) :: dims(1)
  type(c_ptr) :: ptr

  print_handler => OptionCreatePrintHandler(option)
  dims(1) = size(buf,1)
  ptr = c_loc(buf)
  call HDF5DatasetReadBase(loc_id,dset_name,H5T_NATIVE_DOUBLE, &
                           dims,ptr,'1D PetscReal',print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetReadDoubleArray1D1

! ************************************************************************** !

subroutine HDF5DatasetReadDoubleArray1D2(loc_id,dset_name,buf,driver)
  !
  ! Writes a 1D PetscReal dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: dset_name  ! must be of variable length
  PetscReal, target, intent(inout) :: buf(:)
  class(driver_type), intent(in) :: driver

  type(print_handler_type), pointer :: print_handler
  integer(HSIZE_T) :: dims(1)
  type(c_ptr) :: ptr

  print_handler => driver%CreatePrintHandler()
  dims(1) = size(buf,1)
  ptr = c_loc(buf)
  call HDF5DatasetReadBase(loc_id,dset_name,H5T_NATIVE_DOUBLE, &
                           dims,ptr,'1D PetscReal',print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetReadDoubleArray1D2

! ************************************************************************** !

subroutine HDF5DatasetReadInteger1(loc_id,dset_name,buf,option)
  !
  ! Writes a single PetscInt dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: dset_name  ! must be of variable length
  PetscInt, target, intent(in) :: buf
  type(option_type), intent(in)  :: option

  type(print_handler_type), pointer :: print_handler
  integer(HSIZE_T) :: dims(1)
  type(c_ptr) :: ptr

  print_handler => OptionCreatePrintHandler(option)
  dims(1) = 1
  ptr = c_loc(buf)
  call HDF5DatasetReadBase(loc_id,dset_name,H5T_NATIVE_INTEGER, &
                           dims,ptr,'single PetscInt',print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetReadInteger1

! ************************************************************************** !

subroutine HDF5DatasetReadInteger2(loc_id,dset_name,buf,driver)
  !
  ! Writes a single PetscInt dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: dset_name  ! must be of variable length
  PetscInt, target, intent(in) :: buf
  class(driver_type), intent(in) :: driver

  type(print_handler_type), pointer :: print_handler
  integer(HSIZE_T) :: dims(1)
  type(c_ptr) :: ptr

  print_handler => driver%CreatePrintHandler()
  dims(1) = 1
  ptr = c_loc(buf)
  call HDF5DatasetReadBase(loc_id,dset_name,H5T_NATIVE_INTEGER, &
                           dims,ptr,'single PetscInt',print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetReadInteger2

! ************************************************************************** !

subroutine HDF5DatasetReadBase(parent_id,dset_name,dset_type_expected, &
                               dset_shape_expected,buf,dset_err_str, &
                               print_handler)
  !
  ! Reads a dataset from an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  use String_module

  integer(HID_T), intent(in) :: parent_id
  character(len=*), intent(in) :: dset_name
  integer(HID_T), intent(in) :: dset_type_expected
  integer(HSIZE_T) :: dset_shape_expected(:)
  type(c_ptr), intent(inout) :: buf
  character(len=*), intent(in) :: dset_err_str
  type(print_handler_type), intent(in)  :: print_handler

  integer(HID_T) :: dset_id
  integer(HID_T) :: dset_type_in_file
  integer(HSIZE_T) :: size_in_file
  integer(HSIZE_T) :: size_expected
  integer :: ndims_in_file
  integer :: ndims_expected
  integer :: i
  integer :: hdf5_err
  PetscBool :: is_equal

  call h5dopen_f(parent_id,dset_name,dset_id,hdf5_err)
  if (hdf5_err /= 0) then
    call PrintErrorMessage(print_handler, &
                           'Opening dataset "' // trim(dset_name) // '".')
  endif
  call h5ltget_dataset_ndims_f(parent_id,dset_name,ndims_in_file,hdf5_err)
  if (hdf5_err /= 0) then
    call PrintErrorMessage(print_handler, &
                           'Reading ' // trim(dset_err_str) // &
                           ' dataset "' // &
                           trim(dset_name) // '" number of dimensions.')
  endif
  ndims_expected = size(dset_shape_expected)
  if (ndims_in_file /= ndims_expected) then
    call PrintErrorMessage(print_handler, &
                           trim(dset_err_str) // ' dataset "' // &
                           trim(dset_name) // &
                           '" has incorrect dimensions: ' // &
                           trim(StringWrite(ndims_in_file)) // &
                           ' in file versus ' // &
                           trim(StringWrite(ndims_expected)) // &
                           ' expected.')
  endif
  call h5dget_type_f(dset_id,dset_type_in_file,hdf5_err)
  call h5tequal_f(dset_type_expected,dset_type_in_file,is_equal,hdf5_err)
  if (.not.is_equal) then
    call PrintErrorMessage(print_handler, &
                           trim(dset_err_str) // &
                           ' has a mismatch in dataset type for "' // &
                           trim(dset_name) // '" read.')
  endif
  call h5dget_storage_size_f(dset_id,size_in_file,hdf5_err)
  call h5tget_size_f(dset_type_expected,size_expected,hdf5_err)
  do i = 1, ndims_expected
    size_expected = size_expected * dset_shape_expected(i)
  enddo
  if (size_in_file /= size_expected) then
    call PrintErrorMessage(print_handler, &
                           trim(dset_err_str) // ' dataset "' // &
                           trim(dset_name) // &
                           '" has a mismatch in size: ' // &
                           trim(StringWrite(int(size_in_file))) // &
                           ' in file versus ' // &
                           trim(StringWrite(int(size_expected))) // &
                           'expected.')
  endif
  call h5dread_f(dset_id,dset_type_expected,buf,hdf5_err)
  if (hdf5_err /= 0) then
    call PrintErrorMessage(print_handler, &
                           trim(dset_err_str) // ' dataset "' // &
                           trim(dset_name) // &
                           '" has an error while reading data.')
  endif
  call HDF5DatasetClose(dset_id,print_handler)

end subroutine HDF5DatasetReadBase

! ************************************************************************** !

subroutine HDF5DatasetWriteInteger1(loc_id,dset_name,buf,option)
  !
  ! Writes a single PetscInt dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: dset_name  ! must be of variable length
  PetscInt, target, intent(inout) :: buf
  type(option_type), intent(in)  :: option

  type(print_handler_type), pointer :: print_handler
  integer(HSIZE_T) :: dims(1)
  type(c_ptr) :: ptr

  print_handler => OptionCreatePrintHandler(option)
  dims(1) = 1
  ptr = c_loc(buf)
  call HDF5DatasetWriteBase(loc_id,dset_name,H5T_NATIVE_INTEGER,dims, &
                            ptr,'single PetscInt',print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetWriteInteger1

! ************************************************************************** !

subroutine HDF5DatasetWriteInteger2(loc_id,dset_name,buf,driver)
  !
  ! Writes a single PetscInt dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: dset_name  ! must be of variable length
  PetscInt, target, intent(inout) :: buf
  class(driver_type), intent(in) :: driver

  type(print_handler_type), pointer :: print_handler
  integer(HSIZE_T) :: dims(1)
  type(c_ptr) :: ptr

  print_handler => driver%CreatePrintHandler()
  dims(1) = 1
  ptr = c_loc(buf)
  call HDF5DatasetWriteBase(loc_id,dset_name,H5T_NATIVE_INTEGER,dims, &
                            ptr,'single PetscInt',print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetWriteInteger2

! ************************************************************************** !

function HDF5ObjectGetNameTypeString(loc_id,print_handler)
  !
  ! Returns the type and name of an object
  !
  ! Author: Glenn Hammond
  ! Date: 03/03/23
  !
  integer(HID_T), intent(in) :: loc_id
  type(print_handler_type) :: print_handler

  character(len=MAXSTRINGLENGTH) :: HDF5ObjectGetNameTypeString

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  integer :: itype
  integer(SIZE_T) :: buf_size, name_size
  integer :: hdf5_err

  buf_size = len(string)-1
  call h5iget_name_f(loc_id,string,buf_size,name_size,hdf5_err)
  if (hdf5_err /= 0) then
    call PrintErrorMessage(print_handler,'Error getting name of HDF5 object')
  endif
  call h5iget_type_f(loc_id,itype,hdf5_err)
  if (hdf5_err /= 0) then
    call PrintErrorMessage(print_handler,'Error getting type of HDF5 object')
  endif
  ! cannot use a select case here
  if (itype == H5I_FILE_F) then
    word = 'File'
  else if (itype == H5I_GROUP_F) then
    word = 'Group'
  else if (itype == H5I_DATATYPE_F) then
    word = 'Datatype'
  else if (itype == H5I_DATASPACE_F) then
    word = 'Dataspace'
  else if (itype == H5I_DATASET_F) then
    word = 'Dataset'
  else if (itype == H5I_ATTR_F) then
    word = 'Attribute'
  else
    word = 'Unknown object'
  endif

  HDF5ObjectGetNameTypeString = trim(word) // ' "' // trim(string)
  HDF5ObjectGetNameTypeString = trim(HDF5ObjectGetNameTypeString) // '"'

end function HDF5ObjectGetNameTypeString

! ************************************************************************** !

subroutine HDF5CloseCheckError(hdf5_err,obj_id,print_handler)
  !
  ! Checks status of error flag and prints an message if an error
  !
  ! Author: Glenn Hammond
  ! Date: 03/03/23
  !
  integer :: hdf5_err
  integer(HID_T), intent(in) :: obj_id
  type(print_handler_type) :: print_handler

  character(len=MAXSTRINGLENGTH) :: string

  if (hdf5_err < 0) then
    string = 'HDF5 ' // &
             trim(HDF5ObjectGetNameTypeString(obj_id,print_handler)) // &
             ' could not be closed.'
    call PrintErrorMessage(print_handler,string)
  endif

end subroutine HDF5CloseCheckError

! ************************************************************************** !

subroutine HDF5DatasetWriteDoubleArray1D1(loc_id,dset_name,buf,option)
  !
  ! Writes a 1D PetscReal dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: dset_name  ! must be of variable length
  PetscReal, target, intent(in) :: buf(:)
  type(option_type), intent(in) :: option

  type(print_handler_type), pointer :: print_handler
  integer(HSIZE_T) :: dims(1)
  type(c_ptr) :: ptr

  print_handler => OptionCreatePrintHandler(option)
  dims(1) = size(buf)
  ptr = c_loc(buf)
  call HDF5DatasetWriteBase(loc_id,dset_name,H5T_NATIVE_DOUBLE,dims, &
                            ptr,'1D PetscReal',print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetWriteDoubleArray1D1

! ************************************************************************** !

subroutine HDF5DatasetWriteDoubleArray1D2(loc_id,dset_name,buf,driver)
  !
  ! Writes a 1D PetscReal dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: dset_name  ! must be of variable length
  PetscReal, target, intent(in) :: buf(:)
  class(driver_type), intent(in) :: driver

  type(print_handler_type), pointer :: print_handler
  integer(HSIZE_T) :: dims(1)
  type(c_ptr) :: ptr

  print_handler => driver%CreatePrintHandler()
  dims(1) = size(buf)
  ptr = c_loc(buf)
  call HDF5DatasetWriteBase(loc_id,dset_name,H5T_NATIVE_DOUBLE,dims, &
                            ptr,'1D PetscReal',print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5DatasetWriteDoubleArray1D2

! ************************************************************************** !

subroutine HDF5DatasetWriteBase(loc_id,dset_name,dset_type,dset_dims, &
                                buf,dset_err_str,print_handler)
  !
  ! Writes a single PetscInt dataset to an object
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !

  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: dset_name  ! must be of variable length
  integer(HID_T), intent(in) :: dset_type
  integer(HSIZE_T) :: dset_dims(:)
  type(c_ptr) :: buf
  character(len=*), intent(in) :: dset_err_str
  type(print_handler_type) :: print_handler

  integer :: ndims
  integer :: hdf5_err

  ndims = size(dset_dims)
  call h5ltmake_dataset_f(loc_id,dset_name,ndims,dset_dims, &
                          dset_type,buf,hdf5_err)
  if (hdf5_err /= 0) then
    call PrintErrorMessage(print_handler, &
                           'Error writing ' // trim(dset_err_str) // &
                           'dataset "' // trim(dset_name) // '".')
  endif

end subroutine HDF5DatasetWriteBase

! ************************************************************************** !

subroutine HDF5FileTryOpen(filename,file_id,failed,comm)
  !
  ! Attempts to open an HDF5 file. If it fails, it sets the failed flag
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  use Communicator_Aux_module

  character(len=*) :: filename  ! must be of variable length
  integer(HID_T) :: file_id
  PetscBool :: failed
  type(comm_type) :: comm

  character(len=MAXSTRINGLENGTH) :: string
  integer(HID_T) :: prop_id
  PetscMPIInt, parameter :: ON=1, OFF=0
  integer :: hdf5_err, hdf5_err2

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,comm%communicator, &
                          MPI_INFO_NULL,hdf5_err)
#endif
  hdf5_err = 0
  call h5eset_auto_f(OFF,hdf5_err2)
  string = trim(filename)
  call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf5_err,prop_id)
  if (hdf5_err /= 0) failed = PETSC_TRUE
  call h5eset_auto_f(ON,hdf5_err2)
  call h5pclose_f(prop_id,hdf5_err)

end subroutine HDF5FileTryOpen

! ************************************************************************** !

subroutine HDF5FileOpen1(filename,file_id,create,option)
  !
  ! Opens an HDF5 file.  This wrapper provides error messaging if the file
  ! does not exist.
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  use Communicator_Aux_module

  character(len=*) :: filename
  integer(HID_T) :: file_id
  PetscBool :: create
  type(option_type) :: option

  type(print_handler_type), pointer :: print_handler
  type(comm_type), pointer :: comm

  print_handler => OptionCreatePrintHandler(option)
  comm => option%comm
  call HDF5FileOpen(filename,file_id,create,PETSC_TRUE,print_handler,comm)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5FileOpen1

! ************************************************************************** !

subroutine HDF5FileOpen2(filename,file_id,create,collective,option)
  !
  ! Opens an HDF5 file.  This wrapper provides error messaging if the file
  ! does not exist.
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  use Communicator_Aux_module

  character(len=*) :: filename
  integer(HID_T) :: file_id
  PetscBool :: create
  PetscBool :: collective
  type(option_type) :: option

  type(print_handler_type), pointer :: print_handler
  type(comm_type), pointer :: comm

  print_handler => OptionCreatePrintHandler(option)
  comm => option%comm
  call HDF5FileOpen(filename,file_id,create,collective,print_handler,comm)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5FileOpen2

! ************************************************************************** !

subroutine HDF5FileOpen3(filename,file_id,create,driver)
  !
  ! Opens an HDF5 file.  This wrapper provides error messaging if the file
  ! does not exist.
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  use Communicator_Aux_module

  character(len=*) :: filename
  integer(HID_T) :: file_id
  PetscBool :: create
  class(driver_type) :: driver

  type(print_handler_type), pointer :: print_handler
  type(comm_type), pointer :: comm

  print_handler => driver%CreatePrintHandler()
  comm => driver%comm
  call HDF5FileOpen(filename,file_id,create,PETSC_TRUE,print_handler,comm)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5FileOpen3

! ************************************************************************** !

subroutine HDF5FileOpen4(filename,file_id,create,collective,driver)
  !
  ! Opens an HDF5 file.  This wrapper provides error messaging if the file
  ! does not exist.
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/23
  !
  use Communicator_Aux_module

  character(len=*) :: filename
  integer(HID_T) :: file_id
  PetscBool :: create
  PetscBool :: collective
  class(driver_type) :: driver

  type(print_handler_type), pointer :: print_handler
  type(comm_type), pointer :: comm

  print_handler => driver%CreatePrintHandler()
  comm => driver%comm
  call HDF5FileOpen(filename,file_id,create,collective,print_handler,comm)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5FileOpen4

! ************************************************************************** !

subroutine HDF5FileOpen5(filename,file_id,create,collective, &
                         print_handler,comm)
  !
  ! Opens an HDF5 file.  This wrapper provides error messaging if the file
  ! does not exist.
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  use Communicator_Aux_module

  character(len=*) :: filename  ! must be of variable length
  integer(HID_T) :: file_id
  PetscBool :: create
  PetscBool :: collective
  type(print_handler_type) :: print_handler
  type(comm_type) :: comm

  integer(HID_T) :: prop_id
  PetscMPIInt, parameter :: ON=1, OFF=0
  integer :: hdf5_err, hdf5_err2
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (collective) then
    call h5pset_fapl_mpio_f(prop_id,comm%communicator, &
                            MPI_INFO_NULL,hdf5_err)
  endif
#endif
  hdf5_err = 0
  string = trim(filename)
  if (create) then
    call h5eset_auto_f(OFF, hdf5_err2)
    call h5fcreate_f(string,H5F_ACC_TRUNC_F,file_id,hdf5_err, &
                     H5P_DEFAULT_F,prop_id)
    call h5eset_auto_f(ON, hdf5_err2)
  else
    call h5eset_auto_f(OFF, hdf5_err2)
    call h5fopen_f(string,H5F_ACC_RDWR_F,file_id,hdf5_err,prop_id)
    call h5eset_auto_f(ON, hdf5_err2)
  endif
  if (hdf5_err /= 0) then
    if (create) then
      word = 'creating'
    else
      word = 'opening'
    endif
    ! the last arg is whether to print "byrank" which is false if collective
    call PrintErrorMessage(print_handler%print_flags, &
                           comm,print_handler%fid, &
                           'Error ' // trim(word) // ' HDF5 file "' // &
                              trim(filename) // '".', &
                           print_handler%exit_code, &
                           PETSC_FALSE,.not.collective)
  endif
  call h5pclose_f(prop_id,hdf5_err)

end subroutine HDF5FileOpen5

! ************************************************************************** !

subroutine HDF5FileOpenReadOnly1(filename,file_id,is_collective, &
                                 error_string,option,hdfopen_err)
  !
  ! Opens an HDF5 file.  This wrapper provides error messaging if the file
  ! does not exist.
  !
  ! Author: Glenn Hammond
  ! Date: 06/22/15
  !
  character(len=*) :: filename  ! must be of variable length
  integer(HID_T) :: file_id
  PetscBool :: is_collective
  character(len=*) :: error_string
  type(option_type) :: option
  integer, optional :: hdfopen_err

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5FileOpenReadOnly(filename,file_id,is_collective, &
                            error_string,print_handler,option%comm, &
                            hdfopen_err)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5FileOpenReadOnly1

! ************************************************************************** !

subroutine HDF5FileOpenReadOnly2(filename,file_id,is_collective, &
                                 error_string,driver,hdfopen_err)
  !
  ! Opens an HDF5 file.  This wrapper provides error messaging if the file
  ! does not exist.
  !
  ! Author: Glenn Hammond
  ! Date: 06/22/15
  !
  character(len=*) :: filename  ! must be of variable length
  integer(HID_T) :: file_id
  PetscBool :: is_collective
  character(len=*) :: error_string
  class(driver_type) :: driver
  integer, optional :: hdfopen_err

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5FileOpenReadOnly(filename,file_id,is_collective, &
                            error_string,print_handler,driver%comm, &
                            hdfopen_err)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5FileOpenReadOnly2

! ************************************************************************** !

subroutine HDF5FileOpenReadOnly3(filename,file_id,is_collective, &
                                 error_string,print_handler,comm, &
                                 hdfopen_err)
  !
  ! Opens an HDF5 file.  This wrapper provides error messaging if the file
  ! does not exist.
  !
  ! Author: Glenn Hammond
  ! Date: 06/22/15
  !
  use Communicator_Aux_module

  character(len=*) :: filename  ! must be of variable length
  integer(HID_T) :: file_id
  PetscBool :: is_collective
  character(len=*) :: error_string
  type(print_handler_type) :: print_handler
  type(comm_type) :: comm
  integer, optional :: hdfopen_err

  integer(HID_T) :: prop_id
  PetscMPIInt, parameter :: ON=1, OFF=0
  integer :: hdf5_err, hdf5_err2
  character(len=MAXSTRINGLENGTH) :: string

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (is_collective) then
    call h5pset_fapl_mpio_f(prop_id,comm%communicator, &
                            MPI_INFO_NULL,hdf5_err)
  endif
#endif
  call h5eset_auto_f(OFF, hdf5_err2)
  string = trim(filename)
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  ! for non-collective operations, error handing may need to be
  ! handled externally
  if (.not.is_collective .and. present(hdfopen_err)) then
    hdfopen_err = hdf5_err
  else if (hdf5_err /= 0) then
    if (len_trim(error_string) > 0) then
      string = trim(error_string)
    else
      string = 'HDF5 file "' // trim(filename) // '" not found.'
    endif
    call PrintErrorMessage(print_handler,string)
  endif
  call h5eset_auto_f(ON, hdf5_err2)
  call h5pclose_f(prop_id,hdf5_err)

end subroutine HDF5FileOpenReadOnly3

! ************************************************************************** !

subroutine HDF5FileClose1(file_id,option)
  !
  ! Closes an HDF5 file
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  integer(HID_T), intent(in) :: file_id
  type(option_type) :: option

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  call HDF5FileClose(file_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5FileClose1

! ************************************************************************** !

subroutine HDF5FileClose2(file_id,driver)
  !
  ! Closes an HDF5 file
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  integer(HID_T), intent(in) :: file_id
  class(driver_type) :: driver

  type(print_handler_type), pointer :: print_handler

  print_handler => driver%CreatePrintHandler()
  call HDF5FileClose(file_id,print_handler)
  call PrintDestroyHandler(print_handler)

end subroutine HDF5FileClose2

! ************************************************************************** !

subroutine HDF5FileClose3(file_id,print_handler)
  !
  ! Closes an HDF5 file
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/22
  !
  integer(HID_T), intent(in) :: file_id
  type(print_handler_type) :: print_handler

  integer :: hdf5_err

  call h5fclose_f(file_id,hdf5_err)
  call HDF5CloseCheckError(hdf5_err,file_id,print_handler)

end subroutine HDF5FileClose3

! ************************************************************************** !

subroutine HDF5Finalize()
  !
  ! Closes the HDF5 library interface for Fortran.
  !
  ! Author: Glenn Hammond
  ! Date: 07/06/20
  !

  integer :: hdf5_err

  call h5close_f(hdf5_err)

end subroutine HDF5Finalize

end module HDF5_Aux_module
