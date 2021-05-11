module Reaction_Sandbox_Biodeg_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_biodeg_type
    ! Aqueous species
    PetscInt :: species_Aaq_id
    PetscInt :: species_Baq_id
    PetscInt :: species_Caq_id
    PetscInt :: species_Daq_id
    ! Immobile species (e.g. biomass)
    PetscInt :: species_Xim_id
  contains
    procedure, public :: Setup => BiodegSetup
    procedure, public :: Evaluate => BiodegEvaluate
  end type reaction_sandbox_biodeg_type

  public :: BiodegCreate, &
            BiodegSetup

contains

! ************************************************************************** !

function BiodegCreate()
  ! 
  ! Allocates biodegradtion reaction object.
  ! 
  implicit none
  
  class(reaction_sandbox_biodeg_type), pointer :: BiodegCreate

  allocate(BiodegCreate)
  nullify(BiodegCreate%next)  
      
end function BiodegCreate

! ************************************************************************** !

subroutine BiodegSetup(this,reaction,option)
  ! 
  ! Sets up the biodegradation reaction with hardwired parameters
  ! 
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_biodeg_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word

  ! Aqueous species
  word = 'Aaq'
  this%species_Aaq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Baq'
  this%species_Baq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Caq'
  this%species_Caq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Daq'
  this%species_Daq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)

  ! Immobile species
  word = 'Xim'
  this%species_Xim_id = &
    GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
      
end subroutine BiodegSetup

! ************************************************************************** !

subroutine BiodegEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar,reaction, &
                          option)
  ! 
  ! Evaluates biodegradation reaction storing residual but no Jacobian
  ! 
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(reaction_sandbox_biodeg_type) :: this  
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp) ! [mole / sec]
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: volume               ! [m^3 bulk volume]
  PetscReal :: molality_to_molarity ! [kg water / L water]
  PetscReal :: Aaq, Baq, Caq, Daq   ! [mole / L water]
  PetscReal :: Xim                  ! [mole biomass / m^3 bulk volume]
  PetscReal :: k_max                ! [mole rxn / mole biomass]
  PetscReal :: k_decay              ! [1 / sec]
  PetscReal :: K_Aaq, K_Baq, I_Caq  ! [mole / L water]
  PetscReal :: yield                ! [mole biomass / mole rxn]
  PetscReal :: n                    ! [-] Hill constant
  PetscReal :: I_r                  ! [mole rxn / m^3 bulk volume-sec]
  PetscReal :: I                    ! [mole rxn / sec]
  PetscReal :: stoichA, stoichB, stoichC, stoichD ! [mole / mole rxn]
  PetscReal :: RateA, RateB, RateC, RateD, RateX  ! [mole / sec]
  
  volume = material_auxvar%volume        ! den_kg [kg fluid / m^3 fluid]
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
  
  Aaq = rt_auxvar%pri_molal(this%species_Aaq_id)*molality_to_molarity
  Baq = rt_auxvar%pri_molal(this%species_Baq_id)*molality_to_molarity
  Caq = rt_auxvar%pri_molal(this%species_Caq_id)*molality_to_molarity
  Daq = rt_auxvar%pri_molal(this%species_Daq_id)*molality_to_molarity
  Xim = rt_auxvar%immobile(this%species_Xim_id)

  k_max = 9.d-2
  k_decay = 1.d-6
  K_Aaq = 2.d-4
  K_Baq = 1.25d-5
  I_Caq = 2.5d-4

  yield = 1.d-4
  n = 1.d0

  stoichA = -1.d0
  stoichB = -0.25d0
  stoichC = 0.33d0
  stoichD = 1.d0

  I_r = k_max * Xim * Aaq**n / (K_Aaq**n + Aaq**n) * &
                      Baq / (K_Baq + Baq) * &
                      I_Caq / (I_Caq + Caq)

  I = I_r * volume
  RateA = stoichA * I
  RateB = stoichB * I
  RateC = stoichC * I
  RateD = stoichD * I
  RateX = yield * I - k_decay * Xim * volume
  
  ! Note: Always subtract contribution from residual
  Residual(this%species_Aaq_id) = Residual(this%species_Aaq_id) - RateA
  Residual(this%species_Baq_id) = Residual(this%species_Baq_id) - RateB
  Residual(this%species_Caq_id) = Residual(this%species_Caq_id) - RateC
  Residual(this%species_Daq_id) = Residual(this%species_Daq_id) - RateD
  Residual(this%species_Xim_id + reaction%offset_immobile) = &
    Residual(this%species_Xim_id + reaction%offset_immobile) - RateX
  
end subroutine BiodegEvaluate

end module Reaction_Sandbox_Biodeg_class
