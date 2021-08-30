module Debug_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter, public :: DEBUG_ASCII_FORMAT = 1
  PetscInt, parameter, public :: DEBUG_BINARY_FORMAT = 2
  PetscInt, parameter, public :: DEBUG_MATLAB_FORMAT = 3
  PetscInt, parameter, public :: DEBUG_NATIVE_FORMAT = 4

  type, public :: debug_type
    PetscBool :: vecview_residual
    PetscBool :: vecview_solution
    PetscBool :: matview_Matrix
    PetscBool :: matview_Matrix_detailed
    PetscBool :: norm_Matrix

    PetscInt  :: output_format
    PetscBool :: verbose_filename

    PetscBool :: print_couplers
    PetscBool :: print_regions
    character(len=MAXSTRINGLENGTH) :: coupler_string
    PetscBool :: print_waypoints
  end type debug_type

  public :: DebugCreate, &
            DebugRead, &
            DebugCreateViewer, &
            DebugWriteFilename, &
            DebugViewerDestroy, &
            DebugDestroy

contains

! ************************************************************************** !

function DebugCreate()
  !
  ! Create object that stores debugging options for PFLOW
  !
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  !

  implicit none

  type(debug_type), pointer :: DebugCreate

  type(debug_type), pointer :: debug

  allocate(debug)

  debug%vecview_residual = PETSC_FALSE
  debug%vecview_solution = PETSC_FALSE
  debug%matview_Matrix = PETSC_FALSE
  debug%matview_Matrix_detailed = PETSC_FALSE
  debug%norm_Matrix = PETSC_FALSE

  debug%output_format = DEBUG_ASCII_FORMAT
  debug%verbose_filename = PETSC_FALSE

  debug%print_couplers = PETSC_FALSE
  debug%print_regions = PETSC_FALSE
  debug%coupler_string = ''
  debug%print_waypoints = PETSC_FALSE

  DebugCreate => debug

end function DebugCreate

! ************************************************************************** !

subroutine DebugRead(debug,input,option)
  !
  ! Reads debugging data from the input file
  !
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  !

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(debug_type) :: debug
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword','DEBUG')
    call StringToUpper(keyword)

    select case(trim(keyword))

      case('PRINT_SOLUTION','VECVIEW_SOLUTION','VIEW_SOLUTION')
        debug%vecview_solution = PETSC_TRUE
      case('PRINT_RESIDUAL','VECVIEW_RESIDUAL','VIEW_RESIDUAL')
        debug%vecview_residual = PETSC_TRUE
      case('PRINT_JACOBIAN','matview_Matrix','VIEW_JACOBIAN')
        debug%matview_Matrix = PETSC_TRUE
      case('PRINT_JACOBIAN_NORM','norm_Matrix')
        debug%norm_Matrix = PETSC_TRUE
      case('PRINT_MATRIX','MATVIEW_MATRIX','VIEW_MATRIX')
        debug%matview_Matrix = PETSC_TRUE
      case('PRINT_REGIONS')
        debug%print_regions = PETSC_TRUE
      case('PRINT_COUPLERS','PRINT_COUPLER')
        debug%print_couplers = PETSC_TRUE
        debug%coupler_string = trim(adjustl(input%buf))
      case('PRINT_JACOBIAN_DETAILED','matview_Matrix_DETAILED', &
           'VIEW_JACOBIAN_DETAILED')
        debug%matview_Matrix_detailed = PETSC_TRUE
      case('PRINT_WAYPOINTS')
        debug%print_waypoints = PETSC_TRUE
      case('APPEND_COUNTS_TO_FILENAME','APPEND_COUNTS_TO_FILENAMES')
        debug%verbose_filename = PETSC_TRUE
      case('FORMAT')
        call InputReadCard(input,option,keyword)
        call InputErrorMsg(input,option,'keyword','DEBUG,FORMAT')
        call StringToUpper(keyword)
        select case(keyword)
          case('ASCII')
            debug%output_format = DEBUG_ASCII_FORMAT
          case('BINARY')
            debug%output_format = DEBUG_BINARY_FORMAT
          case('MATLAB')
            debug%output_format = DEBUG_MATLAB_FORMAT
          case('NATIVE','PARALLEL')
            debug%output_format = DEBUG_NATIVE_FORMAT
        end select
      case default
        call InputKeywordUnrecognized(input,keyword,'DEBUG',option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine DebugRead

! ************************************************************************** !

subroutine DebugCreateViewer(debug,viewer_name_prefix,option,viewer)
  !
  ! Creates a PETSc viewer for saving PETSc vector or matrix in ASCII or
  ! binary format
  !
  ! Author: Gautam Bisht
  ! Date: 09/23/14
  !

  use Option_module

  implicit none

  type(debug_type), pointer :: debug
  character(len=MAXSTRINGLENGTH), intent(in) :: viewer_name_prefix
  type(option_type) :: option
  PetscViewer, intent (inout) :: viewer

  character(len=MAXWORDLENGTH) :: viewer_name
  PetscErrorCode :: ierr


  select case(debug%output_format)
    case(DEBUG_ASCII_FORMAT)
      viewer_name = trim(viewer_name_prefix) // '.out'
      call PetscViewerASCIIOpen(option%mycomm,viewer_name,viewer, &
                                ierr);CHKERRQ(ierr)
    case(DEBUG_BINARY_FORMAT)
      viewer_name = trim(adjustl(viewer_name_prefix)) // '.bin'
      call PetscViewerBinaryOpen(option%mycomm,viewer_name, &
                                 FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
    case(DEBUG_MATLAB_FORMAT)
      viewer_name = trim(viewer_name_prefix) // '.mat'
      call PetscViewerASCIIOpen(option%mycomm,viewer_name, &
                                 viewer,ierr);CHKERRQ(ierr)
      call PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB, &
                                 ierr);CHKERRQ(ierr)
    case(DEBUG_NATIVE_FORMAT)
      viewer_name = trim(viewer_name_prefix) // '.bin'
      call PetscViewerBinaryOpen(option%mycomm,viewer_name, &
                                 FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerPushFormat(viewer,PETSC_VIEWER_NATIVE, &
                                 ierr);CHKERRQ(ierr)
  end select

end subroutine DebugCreateViewer

! ************************************************************************** !

subroutine DebugWriteFilename(debug,filename,prefix,suffix,ts,ts_cut,ni)
  !
  ! Appends timestep, timestep cut, and Newton iteration counts to a
  ! filename.
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/18
  !

  implicit none

  type(debug_type) :: debug
  character(len=*) :: filename
  character(len=*) :: prefix
  character(len=*) :: suffix
  PetscInt :: ts
  PetscInt :: ts_cut
  PetscInt :: ni

  character(len=MAXWORDLENGTH) :: word

  filename = adjustl(prefix)
  if (debug%verbose_filename) then
    write(word,*) ts
    filename = trim(filename) // '_ts' // adjustl(word)
    write(word,*) ts_cut
    filename = trim(filename) // '_tc' // adjustl(word)
    write(word,*) ni
    filename = trim(filename) // '_ni' // adjustl(word)
  endif
  if (len_trim(suffix) > 0) then
    filename = trim(filename) // '.' // adjustl(suffix)
  endif

end subroutine DebugWriteFilename

! ************************************************************************** !

subroutine DebugViewerDestroy(debug,viewer)
  !
  ! Deallocates PETSc Viewer
  !
  ! Author: Heeho Park
  ! Date: 11/08/18
  !
  implicit none

  type(debug_type) :: debug
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  !geh: must use an 'or' operation since new formats greater than
  !     DEBUG_NATIVE_FORMAT may not require PopFormat()
  if (debug%output_format == DEBUG_MATLAB_FORMAT .or. &
      debug%output_format == DEBUG_NATIVE_FORMAT) then
  !  DEBUG_ASCII_FORMAT = 1
  !  DEBUG_BINARY_FORMAT = 2
  !  DEBUG_MATLAB_FORMAT = 3  popformat required
  !  DEBUG_NATIVE_FORMAT = 4  popformat required
    call PetscViewerPopFormat(viewer,ierr);CHKERRQ(ierr)
  endif
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)


end subroutine DebugViewerDestroy

! ************************************************************************** !

subroutine DebugDestroy(debug)
  !
  ! Deallocates memory associated with debug object
  !
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  !
  implicit none

  type(debug_type), pointer :: debug

  if (.not.associated(debug)) return

  deallocate(debug)
  nullify(debug)

end subroutine DebugDestroy

end module Debug_module
