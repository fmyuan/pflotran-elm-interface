module Reaction_Sandbox_Arctic_GHG_class

#include "petsc/finclude/petscsys.h"
  use petscsys

! 1. Change all references to "Example" as desired to rename the module and
!    and subroutines within the module.

  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

! 2. Add module variables here.  Note that one must use the PETSc data types
!    PetscInt, PetscReal, PetscBool to declare variables of type integer
!    float/real*8, and logical respectively.  E.g.,
!
! PetscReal, parameter :: formula_weight_of_water = 18.01534d0

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_arctic_GHG_type
! 3. Add variables/arrays associated with new reaction.
    character(len=MAXWORDLENGTH) :: species_name
    PetscInt :: species_id
    PetscReal :: rate_constant
  contains
    procedure, public :: ReadInput => ArcticGHGRead
    procedure, public :: Setup => ArcticGHGSetup
    procedure, public :: Evaluate => ArcticGHGEvaluate
    procedure, public :: Destroy => ArcticGHGDestroy
  end type reaction_sandbox_example_type

  public :: ArcticGHGCreate

contains

! ************************************************************************** !

function ArcticGHGCreate()
  !
  ! Allocates arctic GHG reaction object.
  !
  ! Author: David Fukuyama
  ! Date: 05/09/23
  !

  implicit none

  class(reaction_sandbox_example_type), pointer :: ArcticGHGCreate

! 4. Add code to allocate the object, initialize all variables to zero and
!    nullify all pointers. E.g.,
  allocate(ArcticGHGCreate)
  ArcticGHGCreate%species_name = ''
  ArcticGHGCreate%species_id = 0
  ArcticGHGCreate%rate_constant = 0.d0
  nullify(ArcticGHGCreate%next)

end function ArcticGHGCreate

! ************************************************************************** !

subroutine ArcticGHGRead(this,input,option)
  !
  ! Reads input deck for example reaction parameters (if any)
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_example_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,EXAMPLE')
    call StringToUpper(word)

    select case(trim(word))

      ! Example Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !     # begin user-defined input
      !     EXAMPLE
      !       EXAMPLE_INTEGER 1
      !       EXAMPLE_INTEGER_ARRAY 2 3 4
      !     END
      !     # end user defined input
      !   END
      !   ...
      ! END

! 5. Add a case statement for reading variables.
      case('SPECIES_NAME')
! 6. Read the variable.
        ! Read the character string indicating which of the primary species
        ! is being decayed.
        call InputReadWord(input,option,this%species_name,PETSC_TRUE)
! 7. Inform the user of any errors if not read correctly.
        call InputErrorMsg(input,option,'SPECIES_NAME', &
                           'CHEMISTRY,REACTION_SANDBOX,EXAMPLE')
! 8. Repeat for other variables.
      case('RATE_CONSTANT')
        ! Read the double precision rate constant.
        call InputReadDouble(input,option,this%rate_constant)
        ! Note the use of character variable 'word' instead of 'RATE_CONSTANT'
        ! in the error message, as they are identical.
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,EXAMPLE')
        ! Read the optional units and convert to internal
        ! units of 1/s.
        internal_units = 'unitless/sec'
        call InputReadAndConvertUnits(input,this%rate_constant, &
                                internal_units,'CHEMISTRY,REACTION_SANDBOX,&
                                &EXAMPLE,RATE_CONSTANT',option)
      case default
        call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,EXAMPLE',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine ArcticGHGRead

! ************************************************************************** !

subroutine ArcticGHGSetup(this,reaction,option)
  !
  ! Sets up the example reaction with parameters either read from the
  ! input deck or hardwired.
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_example_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize.
  this%species_id = &
    GetPrimarySpeciesIDFromName(this%species_name,reaction,option)

end subroutine ArcticGHGSetup

! ************************************************************************** !

subroutine ArcticGHGEvaluate(this,Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_example_type) :: this
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

  ! Description of subroutine arguments:

  ! Residual - 1D array storing Residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entries in units kg water/sec
  !
  !  Jacobian [kg water/sec] * dc [mol/kg water] = -Res [mol/sec]
  !
  ! compute_derivative - Flag indicating whether analytical derivatives will
  !   be calculated.  The user must provide either the analytical derivatives
  !   or a numerical approximation unless always running with
  !   NUMERICAL_JACOBIAN defined within the NUMERICAL_METHODS TRANSPORT,
  !   NEWTON_SOLVER block of the input deck.  If the use of
  !   NUMERICAL_JACOBIAN is assumed, the user should provide an error
  !   message when compute_derivative is true.  E.g.,
  !
  !   if (compute_derivative) then
  !     option%io_buffer = 'NUMERICAL_JACOBIAN must be specified within &
  !       &the NEWTON_SOLVER block of NUMERICAL_METHODS TRANSPORT due to &
  !       &assumptions made in ExampleEvaluate.'
  !     call PrintErrMsg(option)
  !   endif
  !
  ! rt_auxvar - Object holding chemistry information (e.g., concentrations,
  !   activity coefficients, mineral volume fractions, etc.).  See
  !   reactive_transport_aux.F90.
  !
  !   Useful variables:
  !     rt_auxvar%total(:,iphase) - total component concentrations
  !                                 [mol/L water] for phase
  !     rt_auxvar%pri_molal(:) - free ion concentrations [mol/kg water]
  !     rt_auxvar%pri_act_coef(:) - activity coefficients for primary species
  !     rt_auxvar%aqueous%dtotal(:,iphase) - derivative of total component
  !                 concentration with respect to free ion [kg water/L water]
  !
  ! global_auxvar - Object holding information on flow (e.g., saturation,
  !   density, viscosity, temperature, etc)
  !
  !   Useful variables:
  !     global_auxvar%den(iphase) - liquid density [mol/m^3]
  !     global_auxvar%den_kg(iphase) - liquid density [kg/m^3]
  !     global_auxvar%sat(iphase) - liquid saturation [m^3 water/m^3 pore]
  !     global_auxvar%temp - temperature [C]
  !
  ! porosity - effective porosity of grid cell [m^3 pore/m^3 bulk]
  ! volume - volume of grid cell [m^3]
  ! reaction - Provides access to variable describing chemistry.  E.g.,
  !   reaction%ncomp - # chemical degrees of freedom (mobile and immobile)
  !   reaction%naqcomp - # chemical degrees of freedom on water
  !   reaction%primary_species_names(:) - names of primary species
  !
  ! option - Provides handle for controlling simulation, catching and
  !          reporting errors.

! 10. Add code for the Residual evaluation.

  ! Units of the Residual must be in moles/second.
  ! 1.d3 converts m^3 water -> L water
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  ! Always "subtract" the contribution from the Residual.
  Residual(this%species_id) = Residual(this%species_id) - &
    (-1.d0) * & ! negative stoichiometry
    this%rate_constant * &  ! 1/sec
    L_water * & ! L water
    rt_auxvar%total(this%species_id,iphase) ! mol/L water

  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for the Jacobian evaluation.

    ! Always "add" the contribution to the Jacobian.
    ! Units = (mol/sec)*(kg water/mol) = kg water/sec
    Jacobian(this%species_id,this%species_id) = &
    Jacobian(this%species_id,this%species_id) - &
      (-1.d0) * & ! negative stoichiometry
      this%rate_constant * & ! 1/sec
      L_water * & ! L water
      ! rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) =
      !   derivative of total component concentration with respect to the
      !   free ion concentration of the same species.
      ! kg water/L water
      rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase)

  endif

end subroutine ArcticGHGEvaluate

! ************************************************************************** !

subroutine ArcticGHGDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  implicit none

  class(reaction_sandbox_example_type) :: this

! 12. Add code to deallocate dynamic members of reaction_sandbox_example_type.

end subroutine ArcticGHGDestroy

end module Reaction_Sandbox_ArcticGHG_class
