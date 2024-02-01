module Carbon_Sandbox_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: carbon_sandbox_base_type
    PetscInt, pointer :: specid(:)
    PetscReal, pointer :: stoich(:)
    PetscReal :: rate_constant
    procedure(CarbonEvaluate), pointer :: Evaluate => null()
    class(carbon_sandbox_base_type), pointer :: next
  contains
    procedure, public :: ReadInput => CarbonBaseReadInput
    procedure, public :: Setup => CarbonBaseSetup
    procedure, public :: Destroy => CarbonBaseDestroy
  end type carbon_sandbox_base_type

  interface
    subroutine CarbonEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar, &
                          reaction,option)
      use Global_Aux_module
      use Material_Aux_module
      use Option_module
      use Reaction_Aux_module
      use Reactive_Transport_Aux_module
      import :: carbon_sandbox_base_type
      implicit none
      class(carbon_sandbox_base_type) :: this
      PetscReal :: Residual(:)
      PetscReal :: Jacobian(:,:)
      PetscBool :: compute_derivative
      type(reactive_transport_auxvar_type) :: rt_auxvar
      type(global_auxvar_type) :: global_auxvar
      type(material_auxvar_type) :: material_auxvar
      class(reaction_rt_type) :: reaction
      type(option_type) :: option
    end subroutine CarbonEvaluate
  end interface

  public :: CarbonBaseCreate

contains

! ************************************************************************** !

function CarbonBaseCreate()
  !
  ! Allocates and initializes the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  class(carbon_sandbox_base_type), pointer :: CarbonBaseCreate

  class(carbon_sandbox_base_type), pointer :: this

  allocate(this)
  nullify(this%specid)
  nullify(this%stoich)
  this%rate_constant = UNINITIALIZED_DOUBLE
  this%Evaluate => null()
  nullify(this%next)

  CarbonBaseCreate => this

end function CarbonBaseCreate

! ************************************************************************** !

subroutine CarbonBaseSetup(this,reaction,option)
  !
  ! Configures the reaction and associated data structures
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Option_module
  use Reaction_Aux_module

  class(carbon_sandbox_base_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

end subroutine CarbonBaseSetup

! ************************************************************************** !

subroutine CarbonBaseReadInput(this,input,option)
  !
  ! Reads parameters from input deck
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Option_module
  use Input_Aux_module

  class(carbon_sandbox_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

end subroutine CarbonBaseReadInput

! ************************************************************************** !

subroutine CarbonBaseEvaluate(this,Residual,Jacobian,compute_derivative, &
                              rt_auxvar,global_auxvar,material_auxvar, &
                              reaction,option)
  !
  ! Evaluates the rate expression
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Global_Aux_module
  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Reaction_Inhibition_Aux_module

  class(carbon_sandbox_base_type) :: this
  class(reaction_rt_type) :: reaction
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscBool :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option

  PetscReal :: effective_rate
  PetscReal :: mol_spec(reaction%ncomp)
  PetscReal :: ln_mol_spec(reaction%ncomp)
  PetscReal :: conc,inhibition,inhibition_factor
  PetscReal :: liter_water
  PetscReal :: dummy
  PetscInt :: i, icomp, ncomp

  ncomp = this%specid(0)

  liter_water = material_auxvar%volume* &
                material_auxvar%porosity* &
                global_auxvar%sat(1)*1.d3

  do i = 1, reaction%naqcomp
    mol_spec(i) = rt_auxvar%pri_molal(i)*liter_water
  enddo
  do i = 1, reaction%offset_immobile+reaction%immobile%nimmobile
    mol_spec(reaction%offset_immobile+i) = &
      rt_auxvar%immobile(i)*material_auxvar%volume
  enddo
  ln_mol_spec = log(mol_spec)

  effective_rate = this%rate_constant
  inhibition = 0.d0
  do i = 1, ncomp
    if (this%stoich(i) > 0.d0) cycle
    icomp = this%specid(i)
    effective_rate = effective_rate + this%stoich(i) * ln_mol_spec(icomp)
    conc = mol_spec(icomp) / liter_water
    call ReactionInhibitionSmoothStep(conc,1.d-20,inhibition_factor,dummy)
    inhibition = inhibition + log(inhibition_factor)
  enddo
  effective_rate = exp(effective_rate+inhibition)
  do i = 1, ncomp
    icomp = this%specid(i)
    Residual(icomp) = Residual(icomp) - this%stoich(i)*effective_rate
  enddo

end subroutine CarbonBaseEvaluate

! ************************************************************************** !

subroutine CarbonBaseDestroy(this)
  !
  ! Destroys the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  class(carbon_sandbox_base_type) :: this

end subroutine CarbonBaseDestroy

end module Carbon_Sandbox_Base_class
