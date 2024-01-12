module Reaction_Sand_FlexBioHill_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_BioHill_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, &
    extends(reaction_sandbox_biohill_type) :: reaction_sandbox_flexbiohill_type
    PetscReal :: k_max
    PetscReal :: K_Aaq_n
    PetscReal :: K_Baq
    PetscReal :: I_Caq
    PetscReal :: yield
    PetscReal :: k_decay
    PetscReal :: n
    PetscBool :: molarity_units
    PetscReal, pointer :: stoich(:)
  contains
    procedure, public :: ReadInput => FlexBioHillReadInput
    procedure, public :: Setup => FlexBioHillSetup
    procedure, public :: Evaluate => FlexBioHillEvaluate
    procedure, public :: Destroy => FlexBioHillDestroy
  end type reaction_sandbox_flexbiohill_type

  public :: FlexBioHillCreate

contains

! ************************************************************************** !

function FlexBioHillCreate()
  !
  ! Allocates flexible biodegradation reaction object.
  !
  implicit none

  class(reaction_sandbox_flexbiohill_type), pointer :: FlexBioHillCreate

  allocate(FlexBioHillCreate)
  FlexBioHillCreate%k_max = UNINITIALIZED_DOUBLE
  FlexBioHillCreate%K_Aaq_n = UNINITIALIZED_DOUBLE
  FlexBioHillCreate%K_Baq = UNINITIALIZED_DOUBLE
  FlexBioHillCreate%I_Caq = UNINITIALIZED_DOUBLE
  FlexBioHillCreate%yield = UNINITIALIZED_DOUBLE
  FlexBioHillCreate%k_decay = UNINITIALIZED_DOUBLE
  FlexBioHillCreate%n = 1.d0
  FlexBioHillCreate%molarity_units = PETSC_TRUE
  nullify(FlexBioHillCreate%stoich)
  nullify(FlexBioHillCreate%next)

end function FlexBioHillCreate

! ************************************************************************** !

subroutine FlexBioHillReadInput(this,input,option)
  !
  ! Reads flexible biodegradation reaction parameters
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(reaction_sandbox_flexbiohill_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscReal :: K_Aaq

  K_Aaq = UNINITIALIZED_DOUBLE
  error_string = 'CHEMISTRY,REACTION_SANDBOX,FLEXIBLE_BIODEGRADATION'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    select case(word)
      case('AQUEOUS_CONCENTRATION_UNITS')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case(word)
          case('MOLARITY')
            this%molarity_units = PETSC_TRUE
          case('MOLALITY')
            this%molarity_units = PETSC_FALSE
          case default
            call InputKeywordUnrecognized(input,word, &
                         trim(error_string)//&
                         'AQUEOUS_CONCENTRATION_UNITS',option)
        end select
      case('MAX_SPECIFIC_UTILIZATION_RATE')
        call InputReadDouble(input,option,this%k_max)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%k_max,'1/sec|mol/mol-sec', &
                                      trim(error_string)//','//word,option)
      case('AAQ_HALF_SATURATION_CONSTANT')
        call InputReadDouble(input,option,K_Aaq)
        call InputErrorMsg(input,option,word,error_string)
      case('BAQ_HALF_SATURATION_CONSTANT')
        call InputReadDouble(input,option,this%K_Baq)
        call InputErrorMsg(input,option,word,error_string)
      case('CAQ_MONOD_INHIBITION_CONSTANT')
        call InputReadDouble(input,option,this%I_Caq)
        call InputErrorMsg(input,option,word,error_string)
      case('YIELD')
        call InputReadDouble(input,option,this%yield)
        call InputErrorMsg(input,option,word,error_string)
      case('BIOMASS_DECAY_RATE_CONSTANT')
        call InputReadDouble(input,option,this%k_decay)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%k_decay,'1/sec', &
                                      trim(error_string)//','//word,option)
      case('HILL_EXPONENT')
        call InputReadDouble(input,option,this%n)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  error_string = ''
  if (Uninitialized(this%k_max)) then
    error_string = 'MAX_SPECIFIC_UTILIZATION_RATE,'
  endif
  if (Uninitialized(this%k_decay)) then
    error_string = trim(error_string) // 'BIOMASS_DECAY_RATE_CONSTANT,'
  endif
  if (Uninitialized(K_Aaq)) then
    error_string = trim(error_string) // 'AAQ_HALF_SATURATION_CONSTANT,'
  endif
  if (Uninitialized(this%K_Baq)) then
    error_string = trim(error_string) // 'BAQ_HALF_SATURATION_CONSTANT,'
  endif
  if (Uninitialized(this%I_Caq)) then
    error_string = trim(error_string) // 'CAQ_MONOD_INHIBITION_CONSTANT,'
  endif
  if (Uninitialized(this%yield)) then
    error_string = trim(error_string) // 'YIELD,'
  endif

  if (len_trim(error_string) > 0) then
    option%io_buffer = 'Reaction Sandbox FLEXIBLE_BIODEGRADATION has &
      &uninitialized parameters: ' &
      // error_string(1:len_trim(error_string)-1)
    call PrintErrMsg(option)
  endif

  this%K_Aaq_n = K_Aaq**this%n

end subroutine FlexBioHillReadInput

! ************************************************************************** !

subroutine FlexBioHillSetup(this,reaction,option)
  !
  ! Sets up the flexible biodegradation reaction with hardwired parameters
  !
  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module
  use Option_module

  implicit none

  class(reaction_sandbox_flexbiohill_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  call BioHillSetup(this,reaction,option)
  allocate(this%stoich(reaction%ncomp))
  this%stoich = 0.d0
  this%stoich(this%species_Aaq_id) = -1.d0
  this%stoich(this%species_Baq_id) = -0.25d0
  this%stoich(this%species_Caq_id) = 0.33d0
  this%stoich(this%species_Daq_id) = 1.d0
  this%stoich(this%species_Xim_id+reaction%offset_immobile) = this%yield

end subroutine FlexBioHillSetup

! ************************************************************************** !

subroutine FlexBioHillEvaluate(this,Residual,Jacobian,compute_derivative, &
                               rt_auxvar,global_auxvar,material_auxvar, &
                               reaction,option)
  !
  ! Evaluates flexible biodegradation reaction storing residual and Jacobian
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_flexbiohill_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp) ! [mole / sec]
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp) ! [kg water / sec]
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: volume               ! [m^3 bulk volume]
  PetscReal :: molality_to_molarity ! [kg water / L water]
  PetscReal :: Aaq, Baq, Caq, Daq   ! [mole / L water] or [mole / kg water]
  PetscReal :: Xim                  ! [mole biomass / m^3 bulk volume]
  PetscReal :: I_r                  ! [mole rxn / m^3 bulk volume-sec]
  PetscReal :: I                    ! [mole rxn / sec]
  PetscReal :: dIdx(reaction%ncomp) ! mixed units
  PetscInt :: icomp, jcomp, Xim_offset

  Xim_offset = this%species_Xim_id + reaction%offset_immobile
  volume = material_auxvar%volume
  if (this%molarity_units) then
    molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
  else
    molality_to_molarity = 1.d0
  endif

  Aaq = rt_auxvar%pri_molal(this%species_Aaq_id)*molality_to_molarity
  Baq = rt_auxvar%pri_molal(this%species_Baq_id)*molality_to_molarity
  Caq = rt_auxvar%pri_molal(this%species_Caq_id)*molality_to_molarity
  Daq = rt_auxvar%pri_molal(this%species_Daq_id)*molality_to_molarity
  Xim = rt_auxvar%immobile(this%species_Xim_id)

  I_r = this%k_max * Xim * Aaq**this%n / (this%K_Aaq_n + Aaq**this%n) * &
                           Baq / (this%K_Baq + Baq) * &
                           this%I_Caq / (this%I_Caq + Caq)
  I = I_r * volume

  ! Note: Always subtract contribution from residual
  do icomp = 1, reaction%ncomp
    Residual(icomp) = Residual(icomp) - this%stoich(icomp) * I
  enddo
  ! Add due to negative stoichiometry for decay
  Residual(Xim_offset) = Residual(Xim_offset) + this%k_decay * Xim * volume

  if (compute_derivative) then
    ! Note: Always subtract contribution from Jacobian
    dIdx = 0.d0
    ! [mole rxn - kg water / sec] for Aaq, Baq, Caq
    dIdx(this%species_Aaq_id) = &
      this%n * I / rt_auxvar%pri_molal(this%species_Aaq_id) - &
      I / (this%K_Aaq_n + Aaq**this%n) * &
        this%n * Aaq**this%n / rt_auxvar%pri_molal(this%species_Aaq_id)
    dIdx(this%species_Baq_id) = &
      I / rt_auxvar%pri_molal(this%species_Baq_id) - &
      I / (this%K_Baq + Baq) * &
        Baq / rt_auxvar%pri_molal(this%species_Baq_id)
    dIdx(this%species_Caq_id) = &
      -1.d0 * I / (this%I_Caq + Caq) * &
        Caq / rt_auxvar%pri_molal(this%species_Caq_id)
    ! [mole rxn - m^3 bulk / sec] for Xim
    dIdx(Xim_offset) = I / Xim
    do icomp = 1, reaction%ncomp
      do jcomp = 1, reaction%ncomp
        Jacobian(icomp,jcomp) = Jacobian(icomp,jcomp) - &
          this%stoich(icomp) * dIdx(jcomp)
      enddo
    enddo
    Jacobian(Xim_offset,Xim_offset) = &
      Jacobian(Xim_offset,Xim_offset) + this%k_decay * volume
  endif

end subroutine FlexBioHillEvaluate

! ************************************************************************** !

subroutine FlexBioHillDestroy(this)
  !
  ! Deallocates dynamic memory
  !
  use Utility_module, only : DeallocateArray

  implicit none

  class(reaction_sandbox_flexbiohill_type) :: this

  call DeallocateArray(this%stoich)

end subroutine FlexBioHillDestroy

end module Reaction_Sand_FlexBioHill_class
