module PM_Well_class

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscsnes.h"
  use petscsys
  use petscsnes
  use PM_Base_class
  use Option_module
  use Realization_Subsurface_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type :: well_grid_type
    ! number of well segments
    PetscInt :: nsegments      
    ! delta h discretization of each segment center       
    PetscReal, pointer :: dh(:)      
    ! h coordinate of each segment center (seg#,x-y-z)
    PetscReal, pointer :: h(:,:)    
    ! coordinate of the top/bottom of the well 
    PetscReal :: tophole(3)
    PetscReal :: bottomhole(3)     
    ! gravity vector magnitude
    PetscReal :: g
  end type well_grid_type

  type :: well_type
    ! cross-sectional area of each well segment [calc'd]
    PetscReal, pointer :: area(:)         
    ! diameter of each well segment       
    PetscReal, pointer :: diameter(:) 
    ! volume of each well segment [calc'd]
    PetscReal, pointer :: volume(:) 
    ! friction ceofficient of each well segment        
    PetscReal, pointer :: f(:)      
    ! well index of each well segment [0,1]  0 = cased; 1 = open     
    PetscReal, pointer :: WI(:)    
  end type well_type

  type, public, extends(pm_base_type) :: pm_well_type
    class(realization_subsurface_type), pointer :: realization
    type(well_grid_type), pointer :: grid
    type(well_type), pointer :: well
  contains
    procedure, public :: Setup => PMWellSetup
    procedure, public :: ReadPMBlock => PMWellReadPMBlock
    procedure, public :: SetRealization => PMWellSetRealization
    procedure, public :: InitializeRun => PMWellInitializeRun
    procedure, public :: FinalizeRun => PMWellFinalizeRun
    procedure, public :: InitializeTimestep => PMWellInitializeTimestep
    procedure, public :: UpdateTimestep => PMWellUpdateTimestep
    procedure, public :: FinalizeTimestep => PMWellFinalizeTimestep
    procedure, public :: Residual => PMWellResidual
    procedure, public :: Jacobian => PMWellJacobian
    procedure, public :: PreSolve => PMWellPreSolve
    procedure, public :: PostSolve => PMWellPostSolve
    procedure, public :: Destroy => PMWellDestroy
  end type pm_well_type

  public :: PMWellCreate, PMWellReadPMBlock, PMWellReadPass2

  contains

! ************************************************************************** !

function PMWellCreate()
  ! 
  ! Creates the well process model.
  ! 
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type), pointer :: PMWellCreate
  
  allocate(PMWellCreate)
  call PMBaseInit(PMWellCreate)

  PMWellCreate%header = 'WELLBORE_MODEL'

  nullify(PMWellCreate%realization)

  ! create the grid object:
  allocate(PMWellCreate%grid)
  PMWellCreate%grid%nsegments = UNINITIALIZED_INTEGER
  nullify(PMWellCreate%grid%dh)
  nullify(PMWellCreate%grid%h)
  PMWellCreate%grid%tophole(:) = UNINITIALIZED_DOUBLE
  PMWellCreate%grid%bottomhole(:) = UNINITIALIZED_DOUBLE
  PMWellCreate%grid%g = 9.81

  ! create the well object:
  allocate(PMWellCreate%well)
  nullify(PMWellCreate%well%area)
  nullify(PMWellCreate%well%diameter)
  nullify(PMWellCreate%well%volume)
  nullify(PMWellCreate%well%f)
  nullify(PMWellCreate%well%WI)


end function PMWellCreate

! ************************************************************************** !
  
subroutine PMWellSetup(this)
  ! 
  ! Initializes variables associated with the well process model.
  ! 
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  ! 

  implicit none
  
  class(pm_well_type) :: this
  
  PetscReal :: diff_x,diff_y,diff_z
  PetscReal :: dh_x,dh_y,dh_z
  PetscReal :: total_length
  PetscInt :: k

  allocate(this%grid%dh(this%grid%nsegments))
  allocate(this%grid%h(this%grid%nsegments,3))

  diff_x = this%grid%tophole(1)-this%grid%bottomhole(1)
  diff_y = this%grid%tophole(2)-this%grid%bottomhole(2)
  diff_z = this%grid%tophole(3)-this%grid%bottomhole(3)

  dh_x = diff_x/this%grid%nsegments
  dh_y = diff_y/this%grid%nsegments
  dh_z = diff_z/this%grid%nsegments

  diff_x = diff_x*diff_x
  diff_y = diff_y*diff_y
  diff_z = diff_z*diff_z

  total_length = sqrt(diff_x+diff_y+diff_z)

  do k = 1,this%grid%nsegments
    this%grid%h(k,1) = this%grid%bottomhole(1)+(dh_x*(k-0.5))
    this%grid%h(k,2) = this%grid%bottomhole(2)+(dh_y*(k-0.5))
    this%grid%h(k,3) = this%grid%bottomhole(3)+(dh_z*(k-0.5))
  enddo

  this%grid%dh(:) = total_length/this%grid%nsegments

  allocate(this%well%area(this%grid%nsegments))
  allocate(this%well%volume(this%grid%nsegments))

end subroutine PMWellSetup

! ************************************************************************** !

subroutine PMWellReadPMBlock(this,input)
  ! 
  ! Reads input file parameters associated with the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  use Input_Aux_module
  use String_module
  use Utility_module
  
  implicit none

  class(pm_well_type) :: this
  type(input_type), pointer :: input

  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option
  input%ierr = 0
  error_string = 'WELLBORE_MODEL'
  
  option%io_buffer = 'pflotran card:: WELLBORE_MODEL'
  call PrintMsg(option)

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    found = PETSC_FALSE

    ! Read keywords within WELLBORE_MODEL block:
    select case(trim(word))
    !-------------------------------------
      case('ACTION1')
        ! do some action or assignment
        cycle
    !-------------------------------------
      case('ACTION2')
        ! do some action or assignment
        cycle
    !-------------------------------------
    end select

    ! Read sub-blocks within WELLBORE_MODEL block:
    error_string = 'WELLBORE_MODEL'
    call PMWellReadGrid(this,input,option,word,error_string,found)
    if (found) cycle
    
    error_string = 'WELLBORE_MODEL'
    call PMWellReadWell(this,input,option,word,error_string,found)
    if (found) cycle

    if (.not. found) then
      option%io_buffer = 'Keyword "' // trim(word) // &
                         '" does not exist for WELLBORE_MODEL.'
      call PrintErrMsg(option)
    endif

  enddo
  call InputPopBlock(input,option)

  ! error checking - did the user provide the right amount of vector values
  ! that matches nsegments?

end subroutine PMWellReadPMBlock

! ************************************************************************** !

subroutine PMWellReadGrid(this,input,option,keyword,error_string,found)
  ! 
  ! Reads input file parameters associated with the well model grid.
  ! 
  ! Author: Jenn Frederick
  ! Date: 11/03/2021
  !
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  class(pm_well_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscInt :: num_errors
  character(len=MAXWORDLENGTH) :: word

  error_string = trim(error_string) // ',GRID'
  found = PETSC_TRUE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('GRID')
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !-----------------------------
          case('NUMBER_OF_SEGMENTS')
            call InputReadInt(input,option,this%grid%nsegments)
            call InputErrorMsg(input,option,'NUMBER_OF_SEGMENTS',error_string)
        !-----------------------------
          case('TOP_OF_HOLE')
            call InputReadDouble(input,option,this%grid%tophole(1))
            call InputErrorMsg(input,option,'TOP_OF_HOLE x-coordinate', &
                               error_string)
            call InputReadDouble(input,option,this%grid%tophole(2))
            call InputErrorMsg(input,option,'TOP_OF_HOLE y-coordinate', &
                               error_string)
            call InputReadDouble(input,option,this%grid%tophole(3))
            call InputErrorMsg(input,option,'TOP_OF_HOLE z-coordinate', &
                               error_string)
        !-----------------------------
          case('BOTTOM_OF_HOLE')
            call InputReadDouble(input,option,this%grid%bottomhole(1))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE x-coordinate', &
                               error_string)
            call InputReadDouble(input,option,this%grid%bottomhole(2))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE y-coordinate', &
                               error_string)
            call InputReadDouble(input,option,this%grid%bottomhole(3))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE z-coordinate', &
                               error_string)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      ! ----------------- error messaging -------------------------------------
      if (Uninitialized(this%grid%nsegments)) then
        option%io_buffer = 'ERROR: NUMBER_OF_SEGMENTS must be specified &
                           &in the ' // trim(error_string) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Uninitialized(this%grid%tophole(1)) .or. &
          Uninitialized(this%grid%tophole(2)) .or. &
          Uninitialized(this%grid%tophole(3))) then
        option%io_buffer = 'ERROR: TOP_OF_HOLE must be specified &
                           &in the ' // trim(error_string) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Uninitialized(this%grid%bottomhole(1)) .or. &
          Uninitialized(this%grid%bottomhole(2)) .or. &
          Uninitialized(this%grid%bottomhole(3))) then
        option%io_buffer = 'ERROR: BOTTOM_OF_HOLE must be specified &
                           &in the ' // trim(error_string) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the WELLBORE_MODEL,GRID block. See above error messages.'
    call PrintErrMsg(option)
  endif

  end subroutine PMWellReadGrid

! ************************************************************************** !

subroutine PMWellReadWell(this,input,option,keyword,error_string,found)
  ! 
  ! Reads input file parameters associated with the well model 
  ! well properties.
  ! 
  ! Author: Jenn Frederick
  ! Date: 11/03/2021
  !
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  class(pm_well_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscInt :: k
  PetscInt :: num_errors,read_max,num_read
  character(len=MAXWORDLENGTH) :: word
  PetscReal, pointer :: temp_diameter(:)
  PetscReal, pointer :: temp_friction(:)
  PetscReal, pointer :: temp_well_index(:)

  error_string = trim(error_string) // ',WELL'
  found = PETSC_TRUE
  num_errors = 0

  read_max = 100
  allocate(temp_diameter(read_max))
  allocate(temp_friction(read_max))
  allocate(temp_well_index(read_max))
  temp_diameter(:) = UNINITIALIZED_DOUBLE
  temp_friction(:) = UNINITIALIZED_DOUBLE
  temp_well_index(:) = UNINITIALIZED_DOUBLE

  select case(trim(keyword))
  !-------------------------------------
    case('WELL')
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        num_read = 0
        select case(trim(word))
        !-----------------------------
          case('AREA')
        !-----------------------------
          case('DIAMETER')
            do k = 1,read_max
              call InputReadDouble(input,option,temp_diameter(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for DIAMETER must be &
                &provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(this%well%diameter(num_read))
            this%well%diameter(1:num_read) = temp_diameter(1:num_read)
        !-----------------------------
          case('FRICTION_COEFFICIENT')
            do k = 1,read_max
              call InputReadDouble(input,option,temp_friction(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for FRICTION_COEFFICIENT &
                &must be provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(this%well%f(num_read))
            this%well%f(1:num_read) = temp_friction(1:num_read)
        !-----------------------------
          case('WELL_INDEX')
            do k = 1,read_max
              call InputReadDouble(input,option,temp_well_index(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for WELL_INDEX &
                &must be provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(this%well%WI(num_read))
            this%well%WI(1:num_read) = temp_well_index(1:num_read)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      ! ----------------- error messaging -------------------------------------
      

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                  &the WELLBORE_MODEL,WELL block. See above error messages.'
    call PrintErrMsg(option)
  endif

  deallocate(temp_well_index)
  deallocate(temp_friction)
  deallocate(temp_diameter)

  end subroutine PMWellReadWell

! ************************************************************************** !

subroutine PMWellReadPass2(input,option)
  ! 
  ! Reads input file parameters associated with the well process model.
  ! Second pass dummy read.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 11/03/2021

  use Input_Aux_module
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'SUBSURFACE,WELLBORE_MODEL'

  input%ierr = 0
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
    !--------------------
    case('GRID')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
        enddo
    !--------------------
    case('WELL')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
        enddo
    !--------------------
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine PMWellReadPass2

! ************************************************************************** !

subroutine PMWellSetRealization(this,realization)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  use Realization_Subsurface_class

  implicit none
  
  class(pm_well_type) :: this
  class(realization_subsurface_type), pointer :: realization

  ! There will probably need to be a realization_well created because we need
  ! our own grid and solution vectors and all that. But we will still need to
  ! be connected to the subsurface realization too.
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMWellSetRealization

! ************************************************************************** !

recursive subroutine PMWellInitializeRun(this)
  ! 
  ! Initializes the simulation run for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  
  implicit none

  class(pm_well_type) :: this
  
  ! placeholder
  
end subroutine PMWellInitializeRun

! ************************************************************************** !

recursive subroutine PMWellFinalizeRun(this)
  ! 
  ! Finalizes the simulation run for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  
  ! placeholder 
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMWellFinalizeRun

! ************************************************************************** !

subroutine PMWellInitializeTimestep(this)
  !
  ! Initializes the time step for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  
  implicit none
  
  class(pm_well_type) :: this

  ! placeholder

end subroutine PMWellInitializeTimestep

! ************************************************************************** !

subroutine PMWellUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                               num_newton_iterations,tfac, &
                               time_step_max_growth_factor)
  !
  ! Updates the time step for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  
  implicit none
  
  class(pm_well_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  ! placeholder

end subroutine PMWellUpdateTimestep

! ************************************************************************** !

subroutine PMWellFinalizeTimestep(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this

  ! placeholder
  
end subroutine PMWellFinalizeTimestep

! ************************************************************************** !

subroutine PMWellResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
    
  ! placeholder

end subroutine PMWellResidual

! ************************************************************************** !

subroutine PMWellJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  ! placeholder

end subroutine PMWellJacobian

! ************************************************************************** !

subroutine PMWellPreSolve(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  
  ! placeholder
  
end subroutine PMWellPreSolve

! ************************************************************************** !

subroutine PMWellPostSolve(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  
  ! placeholder
  
end subroutine PMWellPostSolve

! ************************************************************************** !

subroutine PMWellDestroy(this)
  ! 
  ! Destroys the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !
  use Utility_module, only : DeallocateArray
  use Option_module

  implicit none
  
  class(pm_well_type) :: this
    
  call PMBaseDestroy(this)

  call DeallocateArray(this%grid%h)
  
end subroutine PMWellDestroy

! ************************************************************************** !

end module PM_Well_class