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
    PetscInt :: iphase
    PetscReal :: pressure
    PetscReal :: max_mass_rate
    PetscReal :: pressure_span
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
  PressureCreate%iphase = UNINITIALIZED_INTEGER
  PressureCreate%max_mass_rate = UNINITIALIZED_DOUBLE
  PressureCreate%pressure = UNINITIALIZED_DOUBLE
  PressureCreate%pressure_span = 1.d4

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

end subroutine PressureRead

! ************************************************************************** !

subroutine PressureSetup(this,grid,option)
  !
  ! Sets up the pressure src/sink
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/23

  use Option_module
  use Grid_module

  implicit none

  class(srcsink_sandbox_pressure_type) :: this
  type(grid_type) :: grid
  type(option_type) :: option

  call SSSandboxBaseSetup(this,grid,option)

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
  PetscReal :: scale, derivative
  PetscReal :: xmin, xmax

  xmin = this%pressure-this%pressure_span
  xmax = this%pressure

  max_rate_kmol = this%max_mass_rate/FMWH2O
  call Smoothstep(aux_real(this%iphase),xmin,xmax,scale,derivative)
  scale = 1.d0-scale
  derivative = -1.d0 * derivative
  Residual(this%iphase) = max_rate_kmol*scale
!  print *, aux_real(iphase), scale, derivative
!  print *, 'pumping rate: ' // &
!    StringWrite(Residual(1)/aux_real(SS_PRES_LIQUID_DENSITY_OFFSET)*3600.)

  if (compute_derivative) then
    Jacobian(this%iphase,this%iphase) = -max_rate_kmol*derivative
  endif

end subroutine PressureSrcSink

! ************************************************************************** !

subroutine PressureDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: Glenn Hammond
  ! Date: 12/08/23

  implicit none

  class(srcsink_sandbox_pressure_type) :: this

  call SSSandboxBaseDestroy(this)

end subroutine PressureDestroy

end module SrcSink_Sandbox_Pressure_class
