module SrcSink_Sandbox_Mass_Rate_class

! Sandbox srcsink for mass rate source terms.  This source/sink is identical
! to the standard mass rate source/sink in PFLOTRAN, but is more for
! illustrative purposes

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use SrcSink_Sandbox_Base_class

  implicit none

  private

  type, public, &
    extends(srcsink_sandbox_base_type) :: srcsink_sandbox_mass_rate_type
    PetscReal, pointer :: rate(:)
  contains
    procedure, public :: ReadInput => MassRateRead
    procedure, public :: Setup => MassRateSetup
    procedure, public :: Evaluate => MassRateSrcSink
    procedure, public :: Destroy => MassRateDestroy
  end type srcsink_sandbox_mass_rate_type

  public :: MassRateCreate

contains

! ************************************************************************** !

function MassRateCreate()
  !
  ! Allocates mass rate src/sink object.
  !
  ! Author: Glenn Hammond
  ! Date: 05/06/14

  implicit none

  class(srcsink_sandbox_mass_rate_type), pointer :: MassRateCreate

  allocate(MassRateCreate)
  call SSSandboxBaseInit(MassRateCreate)
  nullify(MassRateCreate%rate)


end function MassRateCreate

! ************************************************************************** !

subroutine MassRateRead(this,input,option)
  !
  ! Reads input deck for mass rate src/sink parameters
  !
  ! Author: Glenn Hammond
  ! Date: 05/06/14
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(srcsink_sandbox_mass_rate_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units
  PetscBool :: found

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'SOURCE_SINK_SANDBOX,MASS_RATE')
    call StringToUpper(word)

    call SSSandboxBaseSelectCase(this,input,option,word,found)
    if (found) cycle

    select case(trim(word))

      case('RATE')
        allocate(this%rate(option%nflowdof))
        do i = 1, option%nflowdof
          word = ''
          select case(i)
            case(ONE_INTEGER)
              word = 'Liquid Component Rate'
              internal_units = 'kg/sec'
            case(TWO_INTEGER)
              word = 'Gas Component Rate'
              internal_units = 'kg/sec'
            case(THREE_INTEGER)
              word = 'Energy Rate'
              internal_units = 'MW|MJ/sec'
            case default
              write(word,*) i
              option%io_buffer = 'Unknown dof #' // trim(adjustl(word)) // &
                                 ' in MassRateRead.'
              call PrintErrMsg(option)
          end select
          call InputReadDouble(input,option,this%rate(i))
          call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,XXX_RATE')
          call InputReadAndConvertUnits(input,this%rate(i),internal_units, &
                                        'SOURCE_SINK_SANDBOX,'//trim(word), &
                                        option)
        enddo
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'SRCSINK_SANDBOX,MASS_RATE',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine MassRateRead

! ************************************************************************** !

subroutine MassRateSetup(this,grid,region_list,material_auxvars,option)
  !
  ! Sets up the mass rate src/sink
  !
  ! Author: Glenn Hammond
  ! Date: 05/06/14

  use Option_module
  use Grid_module
  use General_Aux_module, only : general_fmw => fmw_comp
  use WIPP_Flow_Aux_module, only : wipp_flow_fmw => fmw_comp
  use Material_Aux_module, only: material_auxvar_type
  use Region_module

  implicit none

  class(srcsink_sandbox_mass_rate_type) :: this
  type(grid_type) :: grid
  type(region_list_type) :: region_list
  type(material_auxvar_type) :: material_auxvars(:)
  type(option_type) :: option

  call SSSandboxBaseSetup(this,grid,region_list,material_auxvars,option)
  ! convert rate from kg/s to mol/s
  select case(option%iflowmode)
    case(RICHARDS_MODE)
      this%rate(1) = this%rate(1) / FMWH2O
    case(G_MODE,H_MODE)
      this%rate(:) = this%rate(:) / general_fmw(:)
    case(WF_MODE)
      this%rate(:) = this%rate(:) / wipp_flow_fmw(:)
    case default
      option%io_buffer = 'Rate conversion not set up for flow mode in ' // &
                         'MassRateSetup'
      call PrintErrMsg(option)
  end select

end subroutine MassRateSetup

! ************************************************************************** !

subroutine MassRateSrcSink(this,Residual,Jacobian,compute_derivative, &
                           material_auxvar,aux_real,option)
  !
  ! Evaluates src/sink storing residual and/or Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 05/06/14

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(srcsink_sandbox_mass_rate_type) :: this
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: aux_real(:)

  PetscInt :: idof

  do idof = 1, option%nflowdof
    Residual(idof) = this%rate(idof)
  enddo

  if (compute_derivative) then

    ! since the rates are constant, there is no derivative

  endif

end subroutine MassRateSrcSink

! ************************************************************************** !

subroutine MassRateDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: Glenn Hammond
  ! Date: 05/06/14

  implicit none

  class(srcsink_sandbox_mass_rate_type) :: this

  call SSSandboxBaseDestroy(this)
  deallocate(this%rate)
  nullify(this%rate)

end subroutine MassRateDestroy

end module SrcSink_Sandbox_Mass_Rate_class
