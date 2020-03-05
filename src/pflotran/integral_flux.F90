module Integral_Flux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Geometry_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  

  PetscInt, parameter, public :: INTEGRATE_FLOW = 1
  PetscInt, parameter, public :: INTEGRATE_TRANSPORT = 2
  
  PetscInt, parameter, public :: SIGNED_FLUXES = 0
  PetscInt, parameter, public :: POSITIVE_FLUXES_ONLY = 1
  PetscInt, parameter, public :: ABSOLUTE_FLUXES = 2

  type, public :: integral_flux_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    type(point3d_type), pointer :: polygon(:)
    type(plane_type), pointer :: plane
    PetscReal, pointer :: coordinates_and_directions(:,:)
    PetscInt, pointer :: vertices(:,:)
    PetscInt, pointer :: cell_ids(:,:)
    PetscBool :: invert_direction
    PetscInt :: flux_calculation_option !0=signed, 1=positive_only, 2=absolute
    PetscInt, pointer :: internal_connections(:)
    PetscInt, pointer :: boundary_connections(:)
    PetscReal, pointer :: integral_value(:)
    type(integral_flux_type), pointer :: next
  end type integral_flux_type
  
  type, public :: integral_flux_list_type
    PetscInt :: num_integral_fluxes
    type(integral_flux_type), pointer :: first
    type(integral_flux_type), pointer :: last
    type(integral_flux_type), pointer :: array(:)
  end type integral_flux_list_type

  public :: IntegralFluxCreate, &
            IntegralFluxDestroy, &
            IntegralFluxRead, &
            IntegralFluxAddToList, &
            IntegralFluxInitList, &
            IntegralFluxDestroyList, &
            IntegralFluxGetPtrFromList, &
            IntegralFluxSizeStorage, &
            IntegralFluxUpdate, &
            IntegralFluxGetInstantaneous

contains

! ************************************************************************** !

function IntegralFluxCreate()
  ! 
  ! Create object that stores integral flux information
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_type), pointer :: IntegralFluxCreate
  
  type(integral_flux_type), pointer :: integral_flux
  
  allocate(integral_flux)
  
  integral_flux%name = ''
  integral_flux%id = 0
  integral_flux%invert_direction = PETSC_FALSE
  integral_flux%flux_calculation_option = SIGNED_FLUXES !signed
  nullify(integral_flux%polygon)
  nullify(integral_flux%plane)
  nullify(integral_flux%coordinates_and_directions)
  nullify(integral_flux%vertices)
  nullify(integral_flux%cell_ids)
  nullify(integral_flux%internal_connections)
  nullify(integral_flux%boundary_connections)
  nullify(integral_flux%integral_value)
  nullify(integral_flux%next)
  
  IntegralFluxCreate => integral_flux

end function IntegralFluxCreate

! ************************************************************************** !

subroutine IntegralFluxRead(integral_flux,input,option)
  ! 
  ! Reads integral flux data from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  use Utility_module, only : ReallocateArray, DeallocateArray
  
  implicit none
  
  type(integral_flux_type) :: integral_flux
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: icount
  PetscReal :: x(3), y(3), z(3)
  PetscInt, pointer :: int_array(:,:)
  PetscReal, pointer :: real_array(:,:)
  
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)
    
    if (InputCheckExit(input,option)) exit  

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword','integral_flux')   
      
    select case(trim(keyword))
      case('NAME')
        call InputReadWord(input,option,integral_flux%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','INTEGRAL_FLUX')    
      case('INVERT_DIRECTION')
        integral_flux%invert_direction = PETSC_TRUE
      case('FLUXES_OPTION')
        call InputReadCard(input,option,keyword)
        call StringToUpper(keyword)
        select case(trim(keyword))
          case('SIGNED_FLUXES')
            integral_flux%flux_calculation_option = SIGNED_FLUXES
          case('POSITIVE_FLUXES_ONLY')
            integral_flux%flux_calculation_option = POSITIVE_FLUXES_ONLY
          case('ABSOLUTE_FLUXES')
            integral_flux%flux_calculation_option = ABSOLUTE_FLUXES
          case default
            call InputKeywordUnrecognized(input,keyword,'INTEGRAL_FLUX,FLUXES_OPTION',option)
        end select
      case('POLYGON')
        call GeometryReadCoordinates(input,option,integral_flux%name, &
                                     integral_flux%polygon)
      case('PLANE')
        error_string = 'INTEGRAL_FLUX,PLANE'
        icount = 0
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit  
          icount = icount + 1
          ! only use first three points
          if (icount < 4) then
            call InputReadDouble(input,option,x(icount))
            call InputErrorMsg(input,option,'x-coordinate',error_string)
            call InputReadDouble(input,option,y(icount))
            call InputErrorMsg(input,option,'y-coordinate',error_string)
            call InputReadDouble(input,option,z(icount))
            call InputErrorMsg(input,option,'z-coordinate',error_string)
          endif
        enddo
        allocate(integral_flux%plane)
        call GeometryComputePlaneWithPoints(integral_flux%plane, &
                                            x(1),y(1),z(1), &
                                            x(2),y(2),z(2), &
                                            x(3),y(3),z(3))
      case('COORDINATES')
        option%io_buffer = "COORDINATES has been deprecated within the & 
          INTEGRAL_FLUX block in favor of COORDINATES_AND_DIRECTIONS. &
          Please see the user guide for instructions on the new card's use."
        call PrintErrMsg(option)
      case('COORDINATES_AND_DIRECTIONS')
        allocate(real_array(6,100))
        real_array = UNINITIALIZED_DOUBLE
        error_string = 'INTEGRAL_FLUX,COORDINATES_AND_DIRECTIONS'
        icount = 0
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit  
          icount = icount + 1
          if (icount > size(real_array,2)) then
            call ReallocateArray(real_array)
          endif
          call InputReadDouble(input,option,real_array(1,icount))
          call InputErrorMsg(input,option,'x-coordinate',error_string)
          call InputReadDouble(input,option,real_array(2,icount))
          call InputErrorMsg(input,option,'y-coordinate',error_string)
          call InputReadDouble(input,option,real_array(3,icount))
          call InputErrorMsg(input,option,'z-coordinate',error_string)
          call InputReadDouble(input,option,real_array(4,icount))
          call InputErrorMsg(input,option,'x-direction',error_string)
          call InputReadDouble(input,option,real_array(5,icount))
          call InputErrorMsg(input,option,'y-direction',error_string)
          call InputReadDouble(input,option,real_array(6,icount))
          call InputErrorMsg(input,option,'z-direction',error_string)
        enddo
        allocate(integral_flux%coordinates_and_directions(6,icount))
        integral_flux%coordinates_and_directions = real_array(:,1:icount)
        call DeallocateArray(real_array)
      case('VERTICES')
        allocate(int_array(4,100))
        int_array = UNINITIALIZED_INTEGER
        error_string = 'INTEGRAL_FLUX,VERTICES'
        icount = 0
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit  
          icount = icount + 1
          if (icount > size(int_array,2)) then
            call ReallocateArray(int_array)
          endif
          call InputReadInt(input,option,int_array(1,icount))
          call InputErrorMsg(input,option,'vertex 1',error_string)
          call InputReadInt(input,option,int_array(2,icount))
          call InputErrorMsg(input,option,'vertex 2',error_string)
          call InputReadInt(input,option,int_array(3,icount))
          call InputErrorMsg(input,option,'vertex 3',error_string)
          call InputReadInt(input,option,int_array(4,icount))
          ! fourth value is optional, no error message
        enddo
        allocate(integral_flux%vertices(4,icount))
        integral_flux%vertices = int_array(:,1:icount)
        call DeallocateArray(int_array)
      case('CELL_IDS')
        allocate(int_array(2,100))
        int_array = UNINITIALIZED_INTEGER
        error_string = 'INTEGRAL_FLUX,CELL_IDS'
        icount = 0
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit  
          icount = icount + 1
          if (icount > size(int_array,2)) then
            call ReallocateArray(int_array)
          endif
          call InputReadInt(input,option,int_array(1,icount))
          call InputErrorMsg(input,option,'cell id 1',error_string)
          call InputReadInt(input,option,int_array(2,icount))
          call InputErrorMsg(input,option,'cell id 2',error_string)
        enddo
        allocate(integral_flux%cell_ids(2,icount))
        integral_flux%cell_ids = int_array(:,1:icount)
        call DeallocateArray(int_array)
      case default
        call InputKeywordUnrecognized(input,keyword,'INTEGRAL_FLUX',option)
    end select 
  
  enddo  
  call InputPopBlock(input,option)

  if (len_trim(integral_flux%name) < 1) then
    option%io_buffer = 'All INTEGRAL_FLUXes must have a name.'
    call PrintErrMsg(option)
  endif

end subroutine IntegralFluxRead

! ************************************************************************** !

subroutine IntegralFluxSizeStorage(integral_flux,option)
  ! 
  ! Sizes the arrays that store the integrated flux
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use Option_module
  
  implicit none
  
  type(integral_flux_type) :: integral_flux
  type(option_type) :: option
  
  allocate(integral_flux%integral_value(option%nflowdof+option%ntrandof))
  integral_flux%integral_value = 0.d0

end subroutine IntegralFluxSizeStorage

! ************************************************************************** !

subroutine IntegralFluxUpdate(integral_flux_list,internal_fluxes, &
                              boundary_fluxes,iflag,option)
  ! 
  ! Updates the stored integrated value of each integral flux measurement
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use Option_module
  
  implicit none
  
  type(integral_flux_list_type) :: integral_flux_list
  PetscReal :: internal_fluxes(:,:)
  PetscReal :: boundary_fluxes(:,:)
  PetscInt :: iflag
  type(option_type) :: option
  
  type(integral_flux_type), pointer :: integral_flux
  PetscReal, allocatable :: sum_array(:)
  PetscInt :: offset
  PetscInt :: num_values
  PetscReal :: dt
 
  if (.not.associated(integral_flux_list%first)) return

  select case(iflag)
    case(INTEGRATE_FLOW)
      offset = 0
      num_values = option%nflowdof
      dt = option%flow_dt
    case(INTEGRATE_TRANSPORT)
      offset = option%nflowdof
      num_values = option%ntrandof
      dt = option%tran_dt
    case default
      offset = -1 ! to catch bugs
  end select
  
  allocate(sum_array(num_values))
  integral_flux => integral_flux_list%first
  do
    if (.not.associated(integral_flux)) exit
    sum_array = 0.d0
    call IntegralFluxGetInstantaneous(integral_flux, internal_fluxes, &
                                      boundary_fluxes,num_values, &
                                      sum_array,option)
    integral_flux%integral_value(offset+1:offset+num_values) = &
      integral_flux%integral_value(offset+1:offset+num_values) + &
      sum_array(1:num_values)*dt
    integral_flux => integral_flux%next
  enddo
  deallocate(sum_array)

end subroutine IntegralFluxUpdate


! ************************************************************************** !

subroutine IntegralFluxGetInstantaneous(integral_flux, internal_fluxes, &
                                        boundary_fluxes,num_values, &
                                        sum_array,option)
  ! 
  ! Returns the instantaneous mole flux for an integral flux object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use Option_module
  
  implicit none
  
  type(integral_flux_type) :: integral_flux
  PetscReal :: internal_fluxes(:,:)
  PetscReal :: boundary_fluxes(:,:)
  PetscInt :: num_values
  PetscReal :: sum_array(:)
  type(option_type) :: option
  
  PetscInt :: i
  PetscInt :: j
  PetscInt :: iconn
  
  sum_array = 0.d0
  
  if (associated(integral_flux%internal_connections)) then
    do i = 1, size(integral_flux%internal_connections)
      iconn = integral_flux%internal_connections(i)
      select case(integral_flux%flux_calculation_option)
        case(POSITIVE_FLUXES_ONLY)
          if (integral_flux%invert_direction) then
            if (iconn > 0) then
              do j=1, num_values
                sum_array(j) = sum_array(j) + min(internal_fluxes(j,iconn),0.)
              enddo
            else
              do j=1, num_values
                sum_array(j) = sum_array(j) + min(-internal_fluxes(j,-iconn),0.)
              enddo
            endif
          else
            if (iconn > 0) then
              do j=1, num_values
                sum_array(j) = sum_array(j) + max(internal_fluxes(j,iconn),0.)
              enddo
            else
              ! negative connection ids indicate inversion of flux
              do j=1, num_values
                sum_array(j) = sum_array(j) + max(-internal_fluxes(j,-iconn),0.)
              enddo
            endif
          endif
        case(ABSOLUTE_FLUXES)
          if (iconn > 0) then
            sum_array(1:num_values) = sum_array(1:num_values) + &
                                       abs(internal_fluxes(1:num_values,iconn))
          else
            sum_array(1:num_values) = sum_array(1:num_values) + &
                                       abs(-internal_fluxes(1:num_values,-iconn))
          endif
        case(SIGNED_FLUXES)
          if (iconn > 0) then
            sum_array(1:num_values) = sum_array(1:num_values) + &
                                      internal_fluxes(1:num_values,iconn)
          else
            ! negative connection ids indicate inversion of flux
            sum_array(1:num_values) = sum_array(1:num_values) - &
                                      internal_fluxes(1:num_values,-iconn)
          endif
      end select
    enddo
  endif
  
  if (associated(integral_flux%boundary_connections)) then
    do i = 1, size(integral_flux%boundary_connections)
      iconn = integral_flux%boundary_connections(i)
      select case(integral_flux%flux_calculation_option)
        case(POSITIVE_FLUXES_ONLY)
          if (integral_flux%invert_direction) then
            if (iconn > 0) then
              do j=1, num_values
                sum_array(j) = sum_array(j) + min(boundary_fluxes(j,iconn),0.d0)
              enddo
            else
              do j=1, num_values
                sum_array(j) = sum_array(j) + min(-boundary_fluxes(j,-iconn),0.d0)
              enddo
            endif
          else
            if (iconn > 0) then
              do j=1, num_values
                sum_array(j) = sum_array(j) + max(boundary_fluxes(j,iconn),0.d0)
              enddo
            else
              ! negative connection ids indicate inversion of flux
              do j=1, num_values
                sum_array(j) = sum_array(j) + max(-boundary_fluxes(j,-iconn),0.d0)
              enddo
            endif
          endif
        case(ABSOLUTE_FLUXES)
          if (iconn > 0) then
            sum_array(1:num_values) = sum_array(1:num_values) + &
                                       abs(boundary_fluxes(1:num_values,iconn))
          else
            sum_array(1:num_values) = sum_array(1:num_values) + &
                                       abs(-boundary_fluxes(1:num_values,-iconn))
          endif
        case(SIGNED_FLUXES) !default
          if (iconn > 0) then
            sum_array(1:num_values) = sum_array(1:num_values) + &
                                       boundary_fluxes(1:num_values,iconn)
          else
            sum_array(1:num_values) = sum_array(1:num_values) - &
                                       boundary_fluxes(1:num_values,-iconn)
          endif
      end select
    enddo
  endif
  
  if (integral_flux%invert_direction .and. &
      integral_flux%flux_calculation_option == SIGNED_FLUXES) then
    sum_array = -1.d0 * sum_array
  endif

end subroutine IntegralFluxGetInstantaneous

! ************************************************************************** !

subroutine IntegralFluxInitList(list)
  ! 
  ! Initializes a integral flux list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none

  type(integral_flux_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_integral_fluxes = 0

end subroutine IntegralFluxInitList

! ************************************************************************** !

subroutine IntegralFluxAddToList(new_integral_flux,list)
  ! 
  ! Adds a new integral_flux to a integral flux list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_type), pointer :: new_integral_flux
  type(integral_flux_list_type) :: list
  
  list%num_integral_fluxes = list%num_integral_fluxes + 1
  new_integral_flux%id = list%num_integral_fluxes
  if (.not.associated(list%first)) list%first => new_integral_flux
  if (associated(list%last)) list%last%next => new_integral_flux
  list%last => new_integral_flux
  
end subroutine IntegralFluxAddToList

! ************************************************************************** !

function IntegralFluxGetPtrFromList(integral_flux_name,integral_flux_list)
  ! 
  ! Returns a pointer to the integral flux matching integral_flux_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use String_module

  implicit none
  
  type(integral_flux_type), pointer :: IntegralFluxGetPtrFromList
  character(len=MAXWORDLENGTH) :: integral_flux_name
  type(integral_flux_list_type) :: integral_flux_list
 
  PetscInt :: length
  type(integral_flux_type), pointer :: integral_flux
    
  nullify(IntegralFluxGetPtrFromList)
  integral_flux => integral_flux_list%first
  
  do 
    if (.not.associated(integral_flux)) exit
    length = len_trim(integral_flux_name)
    if (length == len_trim(integral_flux%name) .and. &
        StringCompare(integral_flux%name,integral_flux_name, &
                        length)) then
      IntegralFluxGetPtrFromList => integral_flux
      return
    endif
    integral_flux => integral_flux%next
  enddo
  
end function IntegralFluxGetPtrFromList

! ************************************************************************** !

subroutine IntegralFluxDestroyList(integral_flux_list)
  ! 
  ! Deallocates a list of integral fluxes
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_list_type), pointer :: integral_flux_list
  
  type(integral_flux_type), pointer :: integral_flux, prev_integral_flux
  
  if (.not.associated(integral_flux_list)) return
  
  integral_flux => integral_flux_list%first
  do 
    if (.not.associated(integral_flux)) exit
    prev_integral_flux => integral_flux
    integral_flux => integral_flux%next
    call IntegralFluxDestroy(prev_integral_flux)
  enddo
  
  integral_flux_list%num_integral_fluxes = 0
  nullify(integral_flux_list%first)
  nullify(integral_flux_list%last)
  if (associated(integral_flux_list%array)) deallocate(integral_flux_list%array)
  nullify(integral_flux_list%array)
  
  deallocate(integral_flux_list)
  nullify(integral_flux_list)

end subroutine IntegralFluxDestroyList

! ************************************************************************** !

subroutine IntegralFluxDestroy(integral_flux)
  ! 
  ! Deallocates a integral flux
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 
  use Utility_module
  
  implicit none
  
  type(integral_flux_type), pointer :: integral_flux
  
  PetscInt :: i
  
  if (.not.associated(integral_flux)) return
  
  if (associated(integral_flux%polygon)) &
    deallocate(integral_flux%polygon)
  nullify(integral_flux%polygon)
  if (associated(integral_flux%plane)) &
    deallocate(integral_flux%plane)
  nullify(integral_flux%plane) 
  call DeallocateArray(integral_flux%coordinates_and_directions)
  call DeallocateArray(integral_flux%vertices)
  call DeallocateArray(integral_flux%cell_ids)
  call DeallocateArray(integral_flux%internal_connections)
  call DeallocateArray(integral_flux%boundary_connections)
  call DeallocateArray(integral_flux%integral_value)
  deallocate(integral_flux)
  nullify(integral_flux)

end subroutine IntegralFluxDestroy

end module Integral_Flux_module
