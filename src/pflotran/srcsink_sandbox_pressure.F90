module SrcSink_Sandbox_Pressure_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use SrcSink_Sandbox_Base_class

  implicit none

  private

  PetscInt, parameter, public :: SS_PRES_LIQUID_PRESSURE_OFFSET = 1
  PetscInt, parameter, public :: SS_PRES_GAS_PRESSURE_OFFSET = 2
  PetscInt, parameter, public :: SS_PRES_LIQUID_DENSITY_OFFSET = 3
  PetscInt, parameter, public :: SS_PRES_GAS_DENSITY_OFFSET = 4

  type, public, &
    extends(srcsink_sandbox_base_type) :: srcsink_sandbox_pressure_type
    PetscBool :: scale_maximum_mass_rate
    PetscBool :: inhibit_flow_above_pressure
    PetscInt :: iphase
    PetscReal :: pressure
    PetscReal :: max_mass_rate
    PetscReal :: pressure_span
    PetscReal :: sum_volume
  contains
    procedure, public :: ReadInput => PressureRead
    procedure, public :: Setup => PressureSetup
    procedure, public :: Evaluate => PressureSrcSink
    procedure, public :: Destroy => PressureDestroy
  end type srcsink_sandbox_pressure_type

  public :: PressureCreate

contains

! ************************************************************************** !

function PressureCreate()
  !
  ! Allocates pressure src/sink object.
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/23

  implicit none

  class(srcsink_sandbox_pressure_type), pointer :: PressureCreate

  allocate(PressureCreate)
  call SSSandboxBaseInit(PressureCreate)
  PressureCreate%scale_maximum_mass_rate = PETSC_FALSE
  PressureCreate%inhibit_flow_above_pressure = PETSC_TRUE
  PressureCreate%iphase = UNINITIALIZED_INTEGER
  PressureCreate%max_mass_rate = UNINITIALIZED_DOUBLE
  PressureCreate%pressure = UNINITIALIZED_DOUBLE
  PressureCreate%pressure_span = 1.d4
  PressureCreate%sum_volume = 0.d0

end function PressureCreate

! ************************************************************************** !

subroutine PressureRead(this,input,option)
  !
  ! Reads input deck for pressure src/sink parameters
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/23
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(srcsink_sandbox_pressure_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units, word2
  PetscBool :: found
  PetscBool :: inhibit_term_read

  inhibit_term_read = PETSC_FALSE
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'SOURCE_SINK_SANDBOX,PRESSURE')
    call StringToUpper(word)

    call SSSandboxBaseSelectCase(this,input,option,word,found)
    if (found) cycle

    select case(trim(word))

      case('PRESSURE_SPAN')
        call InputReadDouble(input,option,this%pressure_span)
        call InputErrorMsg(input,option,word, &
                           'SOURCE_SINK_SANDBOX,PRESSURE,'//trim(word))
      case('PHASE')
        call InputReadWord(input,option,word2,PETSC_TRUE)
        call InputErrorMsg(input,option,word2, &
                           'SOURCE_SINK_SANDBOX,PRESSURE,'//trim(word))
        call StringToUpper(word2)
        select case(word2)
          case('LIQUID')
            this%iphase = option%liquid_phase
          case('GAS')
            this%iphase = option%gas_phase
          case default
            call InputKeywordUnrecognized(input,word2, &
                                'SRCSINK_SANDBOX,PRESSURE,PHASE',option)
        end select
      case('SCALE_MAXIMUM_MASS_RATE')
        this%scale_maximum_mass_rate = PETSC_TRUE
      case('INHIBIT_FLOW_ABOVE_PRESSURE')
        this%inhibit_flow_above_pressure = PETSC_TRUE
        inhibit_term_read = PETSC_TRUE
      case('INHIBIT_FLOW_BELOW_PRESSURE')
        this%inhibit_flow_above_pressure = PETSC_FALSE
        inhibit_term_read = PETSC_TRUE
      case('MAXIMUM_MASS_RATE')
        call InputReadDouble(input,option,this%max_mass_rate)
        call InputErrorMsg(input,option,word, &
                           'SOURCE_SINK_SANDBOX,PRESSURE,'//trim(word))
        internal_units = 'kg/sec'
        call InputReadAndConvertUnits(input,this%max_mass_rate, &
                                      internal_units, &
                                      'SOURCE_SINK_SANDBOX,PRESSURE,' // &
                                      trim(word),option)
      case('PRESSURE')
        internal_units = 'Pa'
        call InputReadDouble(input,option,this%pressure)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
        call InputReadAndConvertUnits(input,this%pressure,internal_units, &
                                'SOURCE_SINK_SANDBOX,PRESSURE,'//trim(word), &
                                      option)
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'SRCSINK_SANDBOX,PRESSURE',option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (.not.Initialized(this%iphase)) then
    option%io_buffer = 'PHASE LIQUID/GAS must be assigned in a PRESSURE &
      &SOURCE_SINK_SANDBOX.'
    call PrintErrMsg(option)
  endif
  if (.not.Initialized(this%max_mass_rate)) then
    option%io_buffer = 'A MAXIMUM_MASS_RATE must be assigned in a PRESSURE &
      &SOURCE_SINK_SANDBOX.'
    call PrintErrMsg(option)
  endif
  if (.not.Initialized(this%pressure)) then
    option%io_buffer = 'A PRESSURE must be assigned for each phase in a &
      &PRESSURE SOURCE_SINK_SANDBOX.'
    call PrintErrMsg(option)
  endif
  if (.not.inhibit_term_read) then
    option%io_buffer = 'INHIBIT_FLOW_ABOVE_PRESSURE or &
      &INHIBIT_FLOW_BELOW_PRESSURE must be specified in a PRESSURE &
      &SOURCE_SINK_SANDBOX.'
    call PrintErrMsg(option)
  endif
  if (this%inhibit_flow_above_pressure .and. this%max_mass_rate < 0.d0) then
    option%io_buffer = 'A negative maximum mass rate indicating extraction &
      &has been assigned to a PRESSURE SOURCE_SINK_SANDBOX with the flag &
      &INHIBIT_FLOW_ABOVE_PRESSURE. Mass will be continue to be extracted &
      &below the specified PRESSURE.'
    call PrintErrMsg(option)
  endif
  if (.not.this%inhibit_flow_above_pressure .and. &
      this%max_mass_rate > 0.d0) then
    option%io_buffer = 'A postive maximum mass rate indicating injection &
      &has been assigned to a PRESSURE SOURCE_SINK_SANDBOX with the flag &
      &INHIBIT_FLOW_BELOW_PRESSURE. Mass will be continue to be injected &
      &above the specified PRESSURE.'
    call PrintErrMsg(option)
  endif

end subroutine PressureRead

! ************************************************************************** !

subroutine PressureSetup(this,grid,region_list,material_auxvars,option)
  !
  ! Sets up the pressure src/sink
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/23

  use Option_module
  use Grid_module
  use Material_Aux_module, only: material_auxvar_type
  use Region_module

  implicit none

  class(srcsink_sandbox_pressure_type) :: this
  type(grid_type) :: grid
  type(region_list_type) :: region_list
  type(material_auxvar_type) :: material_auxvars(:)
  type(option_type) :: option

  PetscInt :: icell
  PetscErrorCode :: ierr

  call SSSandboxBaseSetup(this,grid,region_list,material_auxvars,option)
  if (this%scale_maximum_mass_rate) then
    this%sum_volume = 0.d0
    if (associated(this%local_cell_ids)) then
      do icell = 1, size(this%local_cell_ids)
        this%sum_volume = this%sum_volume + &
          material_auxvars(grid%nL2G(this%local_cell_ids(icell)))%volume
      enddo
    endif
    call MPI_Allreduce(MPI_IN_PLACE,this%sum_volume,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                       ierr);CHKERRQ(ierr)
  endif

end subroutine PressureSetup

! ************************************************************************** !

subroutine PressureSrcSink(this,Residual,Jacobian,compute_derivative, &
                           material_auxvar,aux_real,option)
  !
  ! Evaluates src/sink storing residual and/or Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/23

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module
  use String_module
  use Utility_module

  implicit none

  class(srcsink_sandbox_pressure_type) :: this
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: aux_real(:)

  PetscReal :: max_rate_kmol
  PetscReal :: smoothstep_scale
  PetscReal :: derivative
  PetscReal :: volume_scale
  PetscReal :: xmin, xmax

  if (this%inhibit_flow_above_pressure) then
    xmin = this%pressure-this%pressure_span
    xmax = this%pressure
  else
    xmin = this%pressure
    xmax = this%pressure+this%pressure_span
  endif

  max_rate_kmol = this%max_mass_rate/FMWH2O
  volume_scale = 1.d0
  if (this%scale_maximum_mass_rate) then
    volume_scale = material_auxvar%volume / this%sum_volume
  endif
  call Smoothstep(aux_real(this%iphase),xmin,xmax,smoothstep_scale,derivative)
  if (this%inhibit_flow_above_pressure) then
    smoothstep_scale = 1.d0-smoothstep_scale
    derivative = -1.d0 * derivative
  endif
  Residual(this%iphase) = max_rate_kmol*volume_scale*smoothstep_scale
!  print *, aux_real(iphase), scale, derivative
!  print *, 'pumping rate: ' // &
!    StringWrite(Residual(1)/aux_real(SS_PRES_LIQUID_DENSITY_OFFSET)*3600.)
!  print *, material_auxvar%id, Residual(this%iphase)

  if (compute_derivative) then
    Jacobian(this%iphase,this%iphase) = -max_rate_kmol*volume_scale*derivative
  endif

end subroutine PressureSrcSink

! ************************************************************************** !

subroutine PressureDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/23

  use Utility_module

  implicit none

  class(srcsink_sandbox_pressure_type) :: this

  call SSSandboxBaseDestroy(this)

end subroutine PressureDestroy

end module SrcSink_Sandbox_Pressure_class
