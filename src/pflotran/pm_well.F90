module PM_Well_class

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscsnes.h"
  use petscsys
  use petscsnes
  use PM_Base_class
  use Option_module
  use Geometry_module
  use Realization_Subsurface_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type :: well_grid_type
    ! number of well segments
    PetscInt :: nsegments      
    ! delta h discretization of each segment center [m]       
    PetscReal, pointer :: dh(:)      
    ! h coordinate of each segment center [m]
    type(point3d_type), pointer :: h(:)  
    ! the local id of the reservoir grid cell within which each segment 
    ! center resides 
    PetscInt, pointer :: h_local_id(:) 
    ! the ghosted id of the reservoir grid cell within which each segment 
    ! center resides 
    PetscInt, pointer :: h_ghosted_id(:)
    ! coordinate of the top/bottom of the well [m]
    PetscReal :: tophole(3)
    PetscReal :: bottomhole(3)     
    ! gravity vector magnitude [m/s2]
    PetscReal :: g
  end type well_grid_type

  type :: well_reservoir_type
    ! reservoir pressure [Pa]    
    PetscReal, pointer :: p(:) 
    ! reservoir liquid saturation
    PetscReal, pointer :: s_l(:)
    ! reservoir gas saturation
    PetscReal, pointer :: s_g(:)
    ! reservoir liquid mobility
    PetscReal, pointer :: mobility_l(:)
    ! reservoir gas mobility
    PetscReal, pointer :: mobility_g(:)
  end type

  type :: well_type
    ! cross-sectional area of each well segment [calc'd] [m2]
    PetscReal, pointer :: area(:)         
    ! diameter of each well segment [m]      
    PetscReal, pointer :: diameter(:) 
    ! volume of each well segment [calc'd] [m3]
    PetscReal, pointer :: volume(:) 
    ! friction ceofficient of each well segment        
    PetscReal, pointer :: f(:)      
    ! well index of each well segment [0,1]  0 = cased; 1 = open     
    PetscReal, pointer :: WI(:)    
    ! well mixture density [kg/m3]    
    PetscReal, pointer :: mixrho(:)
    ! well pressure [Pa]     
    PetscReal, pointer :: p(:) 
    ! well mixture velocity [m/s]   
    PetscReal, pointer :: vm(:)
    ! well bottom of hole pressure BC flag
    PetscBool :: bh_p_set_by_reservoir
    ! well bottom of hole pressure BC [Pa]
    PetscReal :: bh_p 
    ! well top of hole pressure BC [Pa]
    PetscReal :: th_p
  end type well_type

  type :: well_fluid_type
    ! fluid phase ID (liq/gas)
    PetscInt :: ifluid
    ! fluid mobility (phase_permeability/viscosity)
    PetscReal :: mobility
    ! fluid reference density [kg/m3]
    PetscReal :: rho0
    ! fluid density [kg/m3]
    PetscReal, pointer :: rho(:)
    ! fluid saturation
    PetscReal, pointer :: s(:)
    ! fluid source/sink in/out of well [kg/s]??   
    PetscReal, pointer :: Q(:)
    ! equation of state for density
    !procedure(rho_interface), pointer, nopass :: update_rho_ptr => null()
  end type well_fluid_type

  !interface
  !  subroutine rho_interface(fluid)
  !    type(well_fluid_type) :: fluid
  !  end subroutine rho_interface
  !end interface

  type :: well_soln_type
    ! number of primary variables
    PetscInt :: ndof
    ! residual vector
    PetscReal, pointer :: residual(:)
    ! residual vector component - pressure (P)
    PetscReal, pointer :: res_p(:)
    ! residual vector component - mixture velocity (vm)
    PetscReal, pointer :: res_vm(:)
    ! residual vector component - gas saturation (sg)
    PetscReal, pointer :: res_sg(:)
    ! flag for bottom of hole pressure BC
    PetscBool :: bh_p
    ! flag for top of hole pressure BC
    PetscBool :: th_p
  end type well_soln_type

  type, public, extends(pm_base_type) :: pm_well_type
    class(realization_subsurface_type), pointer :: realization
    type(well_grid_type), pointer :: grid
    type(well_type), pointer :: well
    type(well_reservoir_type), pointer :: reservoir 
    type(well_fluid_type), pointer :: liq
    type(well_fluid_type), pointer :: gas
    type(well_soln_type), pointer :: soln
    PetscInt :: nphase 
  contains
    procedure, public :: Setup => PMWellSetup
    procedure, public :: ReadPMBlock => PMWellReadPMBlock
    procedure, public :: SetRealization => PMWellSetRealization
    procedure, public :: InitializeRun => PMWellInitializeRun
    procedure, public :: FinalizeRun => PMWellFinalizeRun
    procedure, public :: InitializeTimestep => PMWellInitializeTimestep
    procedure, public :: UpdateTimestep => PMWellUpdateTimestep
    procedure, public :: FinalizeTimestep => PMWellFinalizeTimestep
    !procedure, public :: Residual => PMWellResidual
    !procedure, public :: Jacobian => PMWellJacobian
    procedure, public :: PreSolve => PMWellPreSolve
    procedure, public :: Solve => PMWellSolve
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

  PMWellCreate%header = 'WELLBORE MODEL'

  nullify(PMWellCreate%realization)
  PMWellCreate%nphase = 0

  ! create the well grid object:
  allocate(PMWellCreate%grid)
  PMWellCreate%grid%nsegments = UNINITIALIZED_INTEGER
  nullify(PMWellCreate%grid%dh)
  nullify(PMWellCreate%grid%h)
  nullify(PMWellCreate%grid%h_local_id)
  nullify(PMWellCreate%grid%h_ghosted_id)
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
  nullify(PMWellCreate%well%mixrho)
  nullify(PMWellCreate%well%p)
  nullify(PMWellCreate%well%vm)
  PMWellCreate%well%bh_p_set_by_reservoir = PETSC_FALSE
  PMWellCreate%well%bh_p = UNINITIALIZED_DOUBLE
  PMWellCreate%well%th_p = UNINITIALIZED_DOUBLE

  ! create the reservoir object:
  allocate(PMWellCreate%reservoir)
  nullify(PMWellCreate%reservoir%p)
  nullify(PMWellCreate%reservoir%s_l)
  nullify(PMWellCreate%reservoir%s_g)
  nullify(PMWellCreate%reservoir%mobility_l)
  nullify(PMWellCreate%reservoir%mobility_g)

  ! create the fluid/liq objects:
  allocate(PMWellCreate%liq)
  PMWellCreate%liq%ifluid = 1 
  PMWellCreate%liq%mobility = UNINITIALIZED_DOUBLE
  PMWellCreate%liq%rho0 = UNINITIALIZED_DOUBLE
  nullify(PMWellCreate%liq%rho)
  nullify(PMWellCreate%liq%s)
  nullify(PMWellCreate%liq%Q)
  !PMWellCreate%liq%update_rho_ptr => PMWellRhoIncompress

  ! create the fluid/gas objects:
  allocate(PMWellCreate%gas)
  PMWellCreate%gas%ifluid = 2 
  PMWellCreate%gas%mobility = UNINITIALIZED_DOUBLE
  PMWellCreate%gas%rho0 = UNINITIALIZED_DOUBLE
  nullify(PMWellCreate%gas%rho)
  nullify(PMWellCreate%gas%s)
  nullify(PMWellCreate%gas%Q)
  !PMWellCreate%gas%update_rho_ptr => PMWellRhoIncompress

  ! create the well solution object:
  allocate(PMWellCreate%soln)
  nullify(PMWellCreate%soln%residual)
  nullify(PMWellCreate%soln%res_p)
  nullify(PMWellCreate%soln%res_vm)
  nullify(PMWellCreate%soln%res_sg)
  PMWellCreate%soln%ndof = UNINITIALIZED_INTEGER
  PMWellCreate%soln%bh_p = PETSC_FALSE
  PMWellCreate%soln%th_p = PETSC_FALSE


end function PMWellCreate

! ************************************************************************** !
  
subroutine PMWellSetup(this)
  ! 
  ! Initializes variables associated with the well process model.
  ! 
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  ! 

  use Option_module
  use Grid_module

  implicit none
  
  class(pm_well_type) :: this
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: res_grid
  PetscReal :: diff_x,diff_y,diff_z
  PetscReal :: dh_x,dh_y,dh_z
  PetscReal :: total_length
  PetscReal :: temp_real
  PetscInt :: local_id
  PetscInt :: nsegments
  PetscInt :: k

  option => this%option
  res_grid => this%realization%patch%grid
  nsegments = this%grid%nsegments

  allocate(this%grid%dh(nsegments))
  allocate(this%grid%h(nsegments))
  allocate(this%grid%h_local_id(nsegments))
  allocate(this%grid%h_ghosted_id(nsegments))

  diff_x = this%grid%tophole(1)-this%grid%bottomhole(1)
  diff_y = this%grid%tophole(2)-this%grid%bottomhole(2)
  diff_z = this%grid%tophole(3)-this%grid%bottomhole(3)

  dh_x = diff_x/nsegments
  dh_y = diff_y/nsegments
  dh_z = diff_z/nsegments

  diff_x = diff_x*diff_x
  diff_y = diff_y*diff_y
  diff_z = diff_z*diff_z

  total_length = sqrt(diff_x+diff_y+diff_z)

  do k = 1,this%grid%nsegments
    this%grid%h(k)%id = k
    this%grid%h(k)%x = this%grid%bottomhole(1)+(dh_x*(k-0.5))
    this%grid%h(k)%y = this%grid%bottomhole(2)+(dh_y*(k-0.5))
    this%grid%h(k)%z = this%grid%bottomhole(3)+(dh_z*(k-0.5))
  enddo

  this%grid%dh(:) = total_length/nsegments

  ! Get the local_id for each well segment center from the reservoir grid
  this%grid%h_local_id(:) = -1  
  do k = 1,this%grid%nsegments
    call GridGetLocalIDFromCoordinate(res_grid,this%grid%h(k), &
                                      option,local_id)
    this%grid%h_local_id(k) = local_id
    this%grid%h_ghosted_id(k) = res_grid%nL2G(local_id)
  enddo

  if (size(this%well%diameter) /= nsegments) then
    if (size(this%well%diameter) == 1) then
      temp_real = this%well%diameter(1)
      deallocate(this%well%diameter)
      allocate(this%well%diameter(nsegments))
      this%well%diameter(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,WELL,DIAMETER must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,WELL,DIAMETER, &
        &it will be applied to all segments of the well.'
      call PrintErrMsg(option)
    endif
  endif
  if (size(this%well%WI) /= nsegments) then
    if (size(this%well%WI) == 1) then
      temp_real = this%well%WI(1)
      deallocate(this%well%WI)
      allocate(this%well%WI(nsegments))
      this%well%WI(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,WELL,WELL_INDEX must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,WELL,WELL_INDEX, &
        &it will be applied to all segments of the well.'
      call PrintErrMsg(option)
    endif
  endif
  if (size(this%well%f) /= nsegments) then
    if (size(this%well%f) == 1) then
      temp_real = this%well%f(1)
      deallocate(this%well%f)
      allocate(this%well%f(nsegments))
      this%well%f(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,WELL,FRICTION_COEFFICIENT must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,WELL,FRICTION_COEFFICIENT, &
        &it will be applied to all segments of the well.'
      call PrintErrMsg(option)
    endif
  endif

  allocate(this%well%area(nsegments))
  this%well%area = 3.14159*(this%well%diameter/2.0)*(this%well%diameter/2.0)

  allocate(this%well%volume(nsegments))
  this%well%volume = this%well%area*this%grid%dh

  this%soln%ndof = this%nphase + 1

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
      case('SINGLE_PHASE')
        this%nphase = 1
        cycle
    !-------------------------------------
      case('TWO_PHASE')
        this%nphase = 2
        cycle
    !-------------------------------------
      case('LIQUID_MOBILITY')
        call InputReadDouble(input,option,this%liq%mobility)
        call InputErrorMsg(input,option,'LIQUID_MOBILITY value', &
                           error_string)
        cycle
    !-------------------------------------
      case('GAS_MOBILITY')
        call InputReadDouble(input,option,this%gas%mobility)
        call InputErrorMsg(input,option,'GAS_MOBILITY value', &
                           error_string)
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

    error_string = 'WELLBORE_MODEL'
    call PMWellReadWellBCs(this,input,option,word,error_string,found)
    if (found) cycle

    if (.not. found) then
      option%io_buffer = 'Keyword "' // trim(word) // &
                         '" does not exist for WELLBORE_MODEL.'
      call PrintErrMsg(option)
    endif

  enddo
  call InputPopBlock(input,option)

  if (this%nphase == 0) then
    option%io_buffer = 'The number of fluid phases must be indicated &
      &in the WELLBORE_MODEL block using one of these keywords: &
      &SINGLE_PHASE, TWO_PHASE.'
    call PrintErrMsg(option)
  endif

  if (Uninitialized(this%liq%mobility)) then
    option%io_buffer = 'LIQUID_MOBILITY must be provided in the &
                       &WELLBORE_MODEL block.'
    call PrintErrMsg(option)
  endif

  if ((this%nphase == 2) .and. (Uninitialized(this%gas%mobility))) then
    option%io_buffer = 'GAS_MOBILITY must be provided in the &
                       &WELLBORE_MODEL block.'
    call PrintErrMsg(option)
  endif
    

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

  read_max = 200
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
      if (.not.associated(this%well%WI)) then
        option%io_buffer = 'Keyword WELL_INDEX must be provided in &
                           &the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
      endif
      if (.not.associated(this%well%f)) then
        option%io_buffer = 'Keyword FRICTION_COEFFICIENT must be provided in &
                           &the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
      endif
      if (.not.associated(this%well%diameter)) then
        option%io_buffer = 'Keyword DIAMETER must be provided in &
                           &the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
      endif    

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

subroutine PMWellReadWellBCs(this,input,option,keyword,error_string,found)
  ! 
  ! Reads input file parameters associated with the well model 
  ! boundary conditions.
  ! 
  ! Author: Jenn Frederick
  ! Date: 12/01/2021
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

  character(len=MAXWORDLENGTH) :: word

  error_string = trim(error_string) // ',WELL_BOUNDARY_CONDITIONS'
  found = PETSC_TRUE

  select case(trim(keyword))
  !-------------------------------------
    case('WELL_BOUNDARY_CONDITIONS')
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
          case('BOTTOM_OF_HOLE')
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
                    case('PRESSURE')
                      call InputReadDouble(input,option,this%well%bh_p)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('PRESSURE_SET_BY_RESERVOIR')
                      this%well%bh_p_set_by_reservoir = PETSC_TRUE
                  !-----------------------------
                    case('FLUX')
                  !-----------------------------
                    case default
                      call InputKeywordUnrecognized(input,word, &
                                                    error_string,option)
                  !-----------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case('TOP_OF_HOLE')
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
                    case('PRESSURE')
                      call InputReadDouble(input,option,this%well%th_p)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('FLUX')
                  !-----------------------------
                    case default
                      call InputKeywordUnrecognized(input,word, &
                                                    error_string,option)
                  !-----------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      ! ----------------- error messaging -------------------------------------
      if (Initialized(this%well%bh_p) .and. &
          this%well%bh_p_set_by_reservoir) then
        option%io_buffer = 'Either keyword BOTTOM_OF_HOLE,PRESSURE or keyword &
          &BOTTOM_OF_HOLE,PRESSURE_SET_BY_RESERVOIR must be provided in the ' &
          // trim(error_string) // ' block, but NOT BOTH.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%well%th_p)) then
        if (Uninitialized(this%well%bh_p) .and. &
            .not.this%well%bh_p_set_by_reservoir)  then
          option%io_buffer = 'Keyword BOTTOM_OF_HOLE,PRESSURE/PRESSURE_SET_&
          &BY_RESERVOIR or keyword TOP_OF_HOLE,PRESSURE must be provided in &
          &the ' // trim(error_string) // ' block.'
          call PrintErrMsg(option) 
        endif
      endif
      if (Initialized(this%well%bh_p) .and. Initialized(this%well%th_p)) then
        option%io_buffer = 'Either keyword BOTTOM_OF_HOLE,PRESSURE or keyword &
          &TOP_OF_HOLE,PRESSURE must be provided in the ' &
          // trim(error_string) // ' block, but NOT BOTH.'
        call PrintErrMsg(option)
      endif
      if (this%well%bh_p_set_by_reservoir .and. &
          Initialized(this%well%th_p)) then
        option%io_buffer = 'Either keyword BOTTOM_OF_HOLE,PRESSURE_SET_BY_&
          &RESERVOIR or keyword TOP_OF_HOLE,PRESSURE must be provided in &
          &the ' // trim(error_string) // ' block, but NOT BOTH.'
        call PrintErrMsg(option)
      endif

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (Initialized(this%well%bh_p)) this%soln%bh_p = PETSC_TRUE
  if (this%well%bh_p_set_by_reservoir) this%soln%bh_p = PETSC_TRUE
  if (Initialized(this%well%th_p)) this%soln%th_p = PETSC_TRUE

  end subroutine PMWellReadWellBCs

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
  character(len=MAXWORDLENGTH) :: card
  
  error_string = 'SUBSURFACE,WELLBORE_MODEL'

  input%ierr = 0
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,card)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
    !--------------------
    case('GRID','WELL')
      call InputSkipToEND(input,option,card)
    !--------------------
    case('WELL_BOUNDARY_CONDITIONS')
      call InputSkipToEND(input,option,card) 
      call InputSkipToEND(input,option,card)
    !--------------------
    end select

  enddo
  call InputPopBlock(input,option)
  call InputSkipToEND(input,option,card)

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
  
  PetscInt :: nsegments 

  nsegments = this%grid%nsegments 

  allocate(this%soln%residual(nsegments*this%soln%ndof))
  this%soln%residual(:) = UNINITIALIZED_DOUBLE

  allocate(this%soln%res_p(nsegments))
  allocate(this%soln%res_vm(nsegments))
  this%soln%res_p(:) = UNINITIALIZED_DOUBLE
  this%soln%res_vm(:) = UNINITIALIZED_DOUBLE

  if (this%nphase == 2) then
    allocate(this%soln%res_sg(nsegments))
    this%soln%res_sg(:) = UNINITIALIZED_DOUBLE
  endif

  allocate(this%well%mixrho(nsegments))
  allocate(this%well%p(nsegments))
  allocate(this%well%vm(nsegments))

  allocate(this%liq%s(nsegments))
  this%liq%rho0 = this%option%flow%reference_density(1)
  allocate(this%liq%rho(nsegments))
  this%liq%rho(:) = this%liq%rho0
  allocate(this%liq%Q(nsegments))
  if (this%nphase == 2) then
    allocate(this%gas%s(nsegments))
    this%gas%rho0 = this%option%flow%reference_density(2)
    allocate(this%gas%rho(nsegments))
    this%gas%rho(:) = this%gas%rho0
    allocate(this%gas%Q(nsegments))
  endif

  allocate(this%reservoir%p(nsegments))
  allocate(this%reservoir%s_l(nsegments))
  allocate(this%reservoir%s_g(nsegments))
  allocate(this%reservoir%mobility_l(nsegments))
  allocate(this%reservoir%mobility_g(nsegments))
  
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
  ! Initializes and takes the time step for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  
  implicit none
  
  class(pm_well_type) :: this

  ! update the reservoir object's pressure and saturations
  call PMWellUpdateReservoir(this)

end subroutine PMWellInitializeTimestep

! ************************************************************************** !

subroutine PMWellUpdateReservoir(this)
  !
  ! Updates the reservoir properties for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

  use WIPP_Flow_Aux_module
  
  implicit none
  
  class(pm_well_type) :: this

  type(wippflo_auxvar_type), pointer :: wippflo_auxvar
  type(option_type), pointer :: option
  PetscInt :: k
  PetscInt :: ghosted_id

  option => this%option 

  do k=1,size(this%grid%h_ghosted_id)
    ghosted_id = this%grid%h_ghosted_id(k)

    wippflo_auxvar => &
      this%realization%patch%aux%wippflo%auxvars(0,ghosted_id)

    this%reservoir%p(k) = wippflo_auxvar%pres(option%liquid_phase)
    this%reservoir%s_l(k) = wippflo_auxvar%sat(option%liquid_phase)
    this%reservoir%s_g(k) = wippflo_auxvar%sat(option%gas_phase)
    this%reservoir%mobility_l(k) = &
      wippflo_auxvar%mobility(option%liquid_phase)
    this%reservoir%mobility_g(k) = wippflo_auxvar%mobility(option%gas_phase)
    ! note: the following other things available in wippflo_auxvar object:
    ! den, den_kg, xmol, kr, mu, effective_porosity, alpha, elevation,
    ! fracture_perm_scaling_factor, klinkenberg 
    if ((k == 1) .and. this%well%bh_p_set_by_reservoir) then
      this%well%bh_p = this%reservoir%p(k)
    endif
  enddo

end subroutine PMWellUpdateReservoir

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

subroutine PMWellRhoIncompress(fluid)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 11/09/2021

  implicit none
  
  type(well_fluid_type) :: fluid

  fluid%rho(:) = fluid%rho0

end subroutine PMWellRhoIncompress

! ************************************************************************** !

subroutine PMWellResidual(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  !SNES :: snes
  !Vec :: xx
  !Vec :: r
  !PetscErrorCode :: ierr
  PetscInt :: i, k 
    
  call PMWellResidualP(this)
  
  call PMWellResidualVm(this)

  if (this%nphase == 2) then
    call PMWellResidualSg(this)
  endif

  ! form the full residual vector
  ! -------------------------------------
  ! |1,2,3, |4,5,6, |7,8,9, |...|
  ! |p,vm,sg|p,vm,sg|p,vm,sg|...|p,vm,sg|
  ! | seg 1 | seg 2 | seg 3 |...| seg n |
  ! -------------------------------------
  ! |bottom | inter | inter |...| top   |
  i = 1
  do k = 1,(this%grid%nsegments*this%soln%ndof)
    if (mod((k-1),this%soln%ndof) == 0) then
      this%soln%residual(k) = this%soln%res_p(i)
      this%soln%residual(k+1) = this%soln%res_vm(i)
      if (this%nphase == 2) then
        this%soln%residual(k+2) = this%soln%res_sg(i) 
      endif 
      i = i + 1
    endif 
  enddo

end subroutine PMWellResidual

! ************************************************************************** !

subroutine PMWellResidualP(this)
  ! 
  ! Calculates the residual for the well pressure. 
  !
  ! Author: Jennifer M. Frederick
  ! Date: 11/09/2021

  implicit none
  
  class(pm_well_type) :: this

  type(well_grid_type), pointer :: grid
  type(well_type), pointer :: well 
  type(well_soln_type), pointer :: soln
  PetscReal :: dpress(this%grid%nsegments)
  PetscReal :: dpress_a(this%grid%nsegments)
  PetscReal :: dpress_f(this%grid%nsegments)
  PetscReal :: dpress_h(this%grid%nsegments)
  PetscInt :: k

  grid => this%grid
  well => this%well
  soln => this%soln

  ! note: booleans soln%bh_p and soln%th_p indicate which pressure boundary 
  !       condition was given (only one could be specified)

  do k = 1,grid%nsegments
    dpress_h(k) = well%mixrho(k)*grid%g*grid%dh(k)
    dpress_f(k) = (well%f(k)*well%mixrho(k)*well%vm(k)* &
                  abs(well%vm(k))*grid%dh(k))/(2.d0*well%diameter(k))
    if (k == 1) then ! top segment
      dpress(k) = 0.0d0 ! hold
      dpress_a(k) = 0.0d0 ! hold
    else if (k == grid%nsegments) then ! bottom segment
      dpress(k) = 0.0d0 ! hold
      dpress_a(k) = 0.0d0 ! hold
    else ! interior segment(s)
      dpress(k) = well%p(k+1)-well%p(k)
      dpress_a(k) = 0.0d0
    endif
  enddo
  soln%res_p = dpress - dpress_h - dpress_f - dpress_a

end subroutine PMWellResidualP

! ************************************************************************** !

subroutine PMWellResidualVm(this)
  ! 
  ! Calculates the residual for the well mixture velocity. 
  !
  ! Author: Jennifer M. Frederick
  ! Date: 11/09/2021

  implicit none
  
  class(pm_well_type) :: this

  type(well_grid_type), pointer :: grid
  type(well_type), pointer :: well 
  type(well_soln_type), pointer :: soln
  PetscReal :: accum(this%grid%nsegments)
  PetscReal :: flux(this%grid%nsegments)
  PetscReal :: srcsink(this%grid%nsegments)
  PetscInt :: k

  grid => this%grid
  well => this%well
  soln => this%soln

  do k = 1,this%grid%nsegments
    if (k == 1) then ! top segment
      accum(k) = 1.0d0
      flux(k) = 1.0d0
      srcsink(k) = 1.0d0
    else if (k == this%grid%nsegments) then ! bottom segment
      accum(k) = 1.0d0
      flux(k) = 1.0d0
      srcsink(k) = 1.0d0
    else ! interior segment(s)
      accum(k) = 1.0d0
      flux(k) = 1.0d0
      srcsink(k) = 1.0d0      
    endif
  enddo
  soln%res_vm = accum + flux - srcsink

end subroutine PMWellResidualVm

! ************************************************************************** !

subroutine PMWellResidualSg(this)
  ! 
  ! Calculates the residual for the well mixture velocity. 
  !
  ! Author: Jennifer M. Frederick
  ! Date: 11/09/2021

  implicit none
  
  class(pm_well_type) :: this

  PetscReal :: accum(this%grid%nsegments)
  PetscReal :: flux(this%grid%nsegments)
  PetscReal :: srcsink(this%grid%nsegments)
  PetscInt :: k

  do k = 1,this%grid%nsegments
    if (k == 1) then ! top segment
      accum(k) = 1.5d0
      flux(k) = 1.5d0
      srcsink(k) = 1.5d0
    else if (k == this%grid%nsegments) then ! bottom segment
      accum(k) = 1.5d0
      flux(k) = 1.5d0
      srcsink(k) = 1.5d0
    else ! interior segment(s)
      accum(k) = 1.5d0
      flux(k) = 1.5d0
      srcsink(k) = 1.5d0
    endif
  enddo
  this%soln%res_sg = accum + flux - srcsink

end subroutine PMWellResidualSg

! ************************************************************************** !

subroutine PMWellJacobian(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  
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

subroutine PMWellSolve(this,time,ierr)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

  implicit none
  
  class(pm_well_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  PetscLogDouble :: log_start_time, log_end_time

  ierr = 0
  call PetscTime(log_start_time, ierr);CHKERRQ(ierr)
  
  call PMWellResidual(this)
  ! the residual equations will include the src/sink term which will be
  ! defined by the difference between p_res and p_well.

  call PMWellJacobian(this)

  ! use Newton's method to solve for the well pressure
  !call PMWellNewton(this)

  ! update the well src/sink Q vector
  call PMWellUpdateWellQ(this)

  ! calculate the phase Darcy velocity
  call PMWellCalcVelocity(this)

  call PetscTime(log_end_time, ierr);CHKERRQ(ierr)
  
end subroutine PMWellSolve

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

subroutine PMWellUpdateWellQ(this)
  !
  ! Updates the src/sink vector for the fluid object.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

  implicit none
  
  class(pm_well_type) :: this

  type(well_fluid_type), pointer :: liq
  type(well_fluid_type), pointer :: gas

  liq => this%liq
  gas => this%gas
  
  liq%Q = liq%rho*liq%mobility*this%well%WI*(this%reservoir%p-this%well%p)
  gas%Q = gas%rho*gas%mobility*this%well%WI*(this%reservoir%p-this%well%p)
  
end subroutine PMWellUpdateWellQ

! ************************************************************************** !

subroutine PMWellCalcVelocity(this)
  !
  ! Calculates the Darcy velocity in the well given the well pressure.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

  implicit none
  
  class(pm_well_type) :: this
  
  ! placeholder
  
end subroutine PMWellCalcVelocity

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

  call DeallocateArray(this%grid%h_local_id)
  call DeallocateArray(this%grid%h_ghosted_id)
  call DeallocateArray(this%grid%dh)
  nullify(this%grid)

  call DeallocateArray(this%reservoir%p)
  call DeallocateArray(this%reservoir%s_l)
  call DeallocateArray(this%reservoir%s_g)
  call DeallocateArray(this%reservoir%mobility_l)
  call DeallocateArray(this%reservoir%mobility_g)
  nullify(this%reservoir)

  call DeallocateArray(this%well%area)
  call DeallocateArray(this%well%diameter)
  call DeallocateArray(this%well%volume)
  call DeallocateArray(this%well%f)
  call DeallocateArray(this%well%WI)
  call DeallocateArray(this%well%mixrho)
  call DeallocateArray(this%well%p)
  call DeallocateArray(this%well%vm)
  nullify(this%well)

  call DeallocateArray(this%liq%rho)
  call DeallocateArray(this%liq%s)
  call DeallocateArray(this%liq%Q)
  call DeallocateArray(this%gas%rho)
  call DeallocateArray(this%gas%s)
  call DeallocateArray(this%gas%Q)
  nullify(this%liq)
  nullify(this%gas)

  call DeallocateArray(this%soln%residual)
  call DeallocateArray(this%soln%res_p)
  call DeallocateArray(this%soln%res_vm)
  call DeallocateArray(this%soln%res_sg)
  nullify(this%soln)
  
end subroutine PMWellDestroy

! ************************************************************************** !

end module PM_Well_class