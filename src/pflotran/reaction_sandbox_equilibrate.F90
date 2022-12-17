module Reaction_Sand_Equilibrate_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_equilibrate_type
    character(len=MAXWORDLENGTH) :: species_name
    PetscInt :: species_id
    PetscReal :: equilibrium_concentration
    PetscReal :: rate_constant
  contains
    procedure, public :: ReadInput => EquilibrateRead
    procedure, public :: Setup => EquilibrateSetup
    procedure, public :: Evaluate => EquilibrateEvaluate
    procedure, public :: Destroy => EquilibrateDestroy
  end type reaction_sandbox_equilibrate_type

  public :: EquilibrateCreate

contains

! ************************************************************************** !

function EquilibrateCreate()
  !
  ! Allocates equilibrate reaction object.
  !
  ! Author: Glenn Hammond
  ! Date: 12/07/22
  !
  implicit none

  class(reaction_sandbox_equilibrate_type), pointer :: EquilibrateCreate

  allocate(EquilibrateCreate)
  EquilibrateCreate%species_name = ''
  EquilibrateCreate%species_id = 0
  EquilibrateCreate%equilibrium_concentration = UNINITIALIZED_DOUBLE
  EquilibrateCreate%rate_constant = UNINITIALIZED_DOUBLE
  nullify(EquilibrateCreate%next)

end function EquilibrateCreate

! ************************************************************************** !

subroutine EquilibrateRead(this,input,option)
  !
  ! Reads input deck for equilibrate reaction parameters
  !
  ! Author: Glenn Hammond
  ! Date: 12/07/22
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_equilibrate_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units
  PetscReal :: half_life

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,EQUILIBRATE')
    call StringToUpper(word)

    select case(trim(word))
      case('SPECIES_NAME')
        call InputReadWord(input,option,this%species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,EQUILIBRATE')
      case('EQUILIBRIUM_CONCENTRATION')
        call InputReadDouble(input,option,this%equilibrium_concentration)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,EQUILIBRATE')
      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,EQUILIBRATE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) then
          input%err_buf = 'REACTION_SANDBOX,EQUILIBRATE,RATE_CONSTANT UNITS'
          call InputDefaultMsg(input,option)
        else
          ! If units exist, convert to internal units of 1/s.
          internal_units = 'unitless/sec'
          this%rate_constant = this%rate_constant * &
            UnitsConvertToInternal(word,internal_units,option)
        endif
      case('HALF_LIFE')
        call InputReadDouble(input,option,half_life)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,EQUILIBRATE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) then
          input%err_buf = 'REACTION_SANDBOX,EQUILIBRATE,HALF_LIFE UNITS'
          call InputDefaultMsg(input,option)
        else
          ! If units exist, convert to internal units of sec.
          internal_units = 'sec'
          half_life = half_life * &
            UnitsConvertToInternal(word,internal_units,option)
        endif
        this%rate_constant = -1.d0*log(0.5d0)/half_life
      case default
        call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,EQUILIBRATE',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine EquilibrateRead

! ************************************************************************** !

subroutine EquilibrateSetup(this,reaction,option)
  !
  ! Sets up the equilibrate reaction with parameters either read from the
  ! input deck or hardwired.
  !
  ! Author: Glenn Hammond
  ! Date: 12/07/22
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_equilibrate_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  this%species_id = &
    GetPrimarySpeciesIDFromName(this%species_name,reaction,option)

  if (Uninitialized(this%equilibrium_concentration)) then
    option%io_buffer = 'An EQUILIBRIUM_CONCENTRATION must be specified for &
      &species "' // trim(this%species_name) // &
      '" in REACTION_SANDBOX EQUILIBRATE.'
    call PrintErrMsg(option)
  endif

  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'A RATE_CONSTANT or HALF_LIFE must be specified for &
      &species "' // trim(this%species_name) // &
      '" in REACTION_SANDBOX EQUILIBRATE. If you are unsure what &
      &to specify, use the HALF_LIFE keyword and begin by setting &
      &the value to 0.01 times the time step size (0.01*dt).'
    call PrintErrMsg(option)
  endif

end subroutine EquilibrateSetup

! ************************************************************************** !

subroutine EquilibrateEvaluate(this,Residual,Jacobian,compute_derivative, &
                               rt_auxvar,global_auxvar,material_auxvar, &
                               reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 12/07/22
  !
  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_equilibrate_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water

  ! Units of the Residual must be in moles/second.
  ! 1.d3 converts m^3 water -> L water
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  Residual(this%species_id) = Residual(this%species_id) - &
    (-1.d0) * & ! negative stoichiometry
    this%rate_constant * &  ! 1/sec
    L_water * & ! L water
    (rt_auxvar%total(this%species_id,iphase) - & ! mol/L water
     this%equilibrium_concentration)

  if (compute_derivative) then

    ! Units = (mol/sec)*(kg water/mol) = kg water/sec
    Jacobian(this%species_id,this%species_id) = &
    Jacobian(this%species_id,this%species_id) - &
      (-1.d0) * & ! negative stoichiometry
      this%rate_constant * & ! 1/sec
      L_water * & ! L water
      ! kg water/L water
      rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase)

  endif

end subroutine EquilibrateEvaluate

! ************************************************************************** !

subroutine EquilibrateDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: Glenn Hammond
  ! Date: 12/07/22
  !
  implicit none

  class(reaction_sandbox_equilibrate_type) :: this

end subroutine EquilibrateDestroy

end module Reaction_Sand_Equilibrate_class
