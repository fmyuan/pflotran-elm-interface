module Reaction_Sandbox_Gas_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, &
       extends(reaction_sandbox_base_type) :: reaction_sandbox_gas_type
    PetscInt :: nspecies
    ! aqueous species stored in a vector
    PetscInt, allocatable :: aq_vec(:)
    ! immobile (sorbed) species stored in a vector
    PetscInt, allocatable :: im_vec(:)
    ! gas species stored in a vector
    PetscInt, allocatable :: gas_vec(:)
    PetscInt :: material_id_skip

    ! create an allocatable array of word-length strings
    character(len=MAXWORDLENGTH), allocatable :: name_vec(:)

    PetscReal, allocatable :: k(:) ! rate constant [mol/m^3 gas]
    PetscReal, allocatable :: Keq(:) ! equilibrium constant [L_gas/m^3 bulk]

  contains
    procedure, public :: ReadInput => GasRead
    procedure, public :: Setup => GasSetup
    procedure, public :: Evaluate => GasReact
    procedure, public :: Destroy => GasDestroy
  end type reaction_sandbox_gas_type

  public :: GasCreate

contains

! ************************************************************************** !

function GasCreate()
  !
  ! Allocates gas reaction object.
  !
  ! Author: Kris Kuhlman
  ! Date: July 2018
  !

  implicit none

  class(reaction_sandbox_gas_type), pointer :: GasCreate
  PetscInt, parameter :: ns = 0

  allocate(GasCreate)
  GasCreate%nspecies = ns
  allocate(GasCreate%aq_vec(ns))
  allocate(GasCreate%im_vec(ns))
  allocate(GasCreate%gas_vec(ns))
  GasCreate%aq_vec = -999
  GasCreate%im_vec = -999
  GasCreate%gas_vec = -999
  allocate(GasCreate%k(ns))
  allocate(GasCreate%Keq(ns))
  GasCreate%k = 0.0d0
  GasCreate%Keq = 0.0d0
  GasCreate%material_id_skip = 0
  allocate(GasCreate%name_vec(ns))
  nullify(GasCreate%next)

end function GasCreate

! ************************************************************************** !

subroutine GasRead(this,input,option)
  !
  ! Reads input deck for gas reaction parameters
  !
  ! Author: Kris Kuhlman
  ! Date: July 2018
  !
  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  class(reaction_sandbox_gas_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=21) :: buffer
  PetscInt :: i, previous_ns
  character(len=MAXWORDLENGTH) :: word

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
         'CHEMISTRY,REACTION_SANDBOX,GAS')

    call StringToUpper(word)
    select case(trim(word))
    case('NUM_GAS_SPECIES')
      ! need to call this first (or uses default)
      ! if this comes after reading other components, they will be blown away
      previous_ns = this%nspecies
      call InputReadInt(input,option,this%nspecies)
      call InputErrorMsg(input,option,'num_gas_species', &
           'CHEMISTRY,REACTION_SANDBOX,GAS')
      if (previous_ns /= this%nspecies) then
        deallocate(this%aq_vec);allocate(this%aq_vec(this%nspecies))
        deallocate(this%im_vec);allocate(this%im_vec(this%nspecies))
        deallocate(this%gas_vec);allocate(this%gas_vec(this%nspecies))
        this%aq_vec = -999
        this%im_vec = -999
        this%gas_vec = -999
        deallocate(this%k);allocate(this%k(this%nspecies))
        deallocate(this%Keq);allocate(this%Keq(this%nspecies))
        this%k = 0.0d0
        this%Keq = 0.0d0
        deallocate(this%name_vec);allocate(this%name_vec(this%nspecies))
      end if

    case('GAS_SPECIES_NAMES')
      buffer = 'gas_species_names_XXX'
      do i=1, this%nspecies
        call InputReadWord(input,option,word,PETSC_TRUE)
        this%name_vec(i) = trim(word)
        write(buffer(19:21),'(I3.3)') i
        call InputErrorMsg(input,option,buffer, &
             'CHEMISTRY,REACTION_SANDBOX,GAS')
      end do

    case('RATE_CONSTANTS')
      ! units of: moles / (m^3 gas * sec)
      call InputReadNDoubles(input,option,this%k,this%nspecies)
      call InputErrorMsg(input,option,'rate_constants', &
           'CHEMISTRY,REACTION_SANDBOX,GAS')

    case('EQUILIBRIUM_CONSTANTS')
      ! units of: liters gas / (m^3 bulk )
      call InputReadNDoubles(input,option,this%Keq,this%nspecies)
      call InputErrorMsg(input,option,'equilibrium_constants', &
           'CHEMISTRY,REACTION_SANDBOX,GAS')

    case('MATERIAL_ID_SKIP')
      ! read in material ID to skip gas sorption reaction calc
      call InputReadInt(input,option,this%material_id_skip)
      call InputErrorMsg(input,option,'material_id_skip', &
           'CHEMISTRY,REACTION_SANDBOX,GAS')

    case default
      call InputKeywordUnrecognized(input,word, &
           'CHEMISTRY,REACTION_SANDBOX,GAS',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine GasRead

! ************************************************************************** !

subroutine GasSetup(this,reaction,option)
  !
  ! Sets up the gas reaction either with parameters either
  ! read from the input deck or hardwired.
  !
  ! Author: Kris Kuhlman
  ! Date: July 2018
  !

  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module
  use Reaction_Gas_Aux_module
  use Option_module

  implicit none

  class(reaction_sandbox_gas_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word
  character(len=4), parameter :: aq = '(aq)', im = '(im)'
  character(len=3) :: g = '(g)'

  ! species read from file must exist in database in both (aq) and (g) form

  do i = 1, this%nspecies
    word = trim(this%name_vec(i))//aq
    this%aq_vec(i) = ReactionGetPriSpeciesIDFromName(word,reaction,option)
    word = trim(this%name_vec(i))//im
    this%im_vec(i) = ImmobileGetSpeciesIDFromName(word,reaction%immobile,option)
    word = trim(this%name_vec(i))//g
    this%gas_vec(i) = GasGetIDFromName(reaction%gas,word)
  end do

end subroutine GasSetup

! ************************************************************************** !

subroutine GasReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  !
  ! Evaluates reaction storing residual and/or Jacobian
  !
  ! Author: Kris Kuhlman
  ! Date: July 2018
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_gas_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative

  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt :: i
  PetscInt, parameter :: iphase = 2

  PetscReal :: im_rate  ! specified rate (mol/sec)
  PetscReal :: g_rate  ! specified rate (mol/sec)
  PetscReal :: Q  ! computed/observed Cs/Cg ratio
  PetscReal :: Cg ! concentration in gas: mol/(L gas)
  !PetscReal :: Cs ! concentration in solid: mol/(m^3 bulk)
  PetscReal :: RT

  if (material_auxvar%id /= this%material_id_skip) then

    ! converts pressure in Pa to moles per m^3 gas
    RT = IDEAL_GAS_CONSTANT*(global_auxvar%temp + 273.15d0)*1.0d3

    do i = 1, this%nspecies

      ! rate expression
      ! rate = -k*(1.d0-Q/K)
      !   k = rate constant [mol/m^3 gas-sec]
      !   Keq = equilibrium constant [L gas/m^3 bulk]
      !   Q = ratio of gas concentration [Cg] to sorbed concentration [Cs]
      !   Cg = gas concentration [mol/L gas]
      !   Cs = gas concentration [mol/m^3 bulk]
      !
      !   Cs = Keq * Cg, therefore, Keq = Cs/Cg at equilibrium
      !   Q = Cs/Cg as measured

      ! gas concentration from partial pressure
      Cg = rt_auxvar%gas_pp(this%gas_vec(i)) * 1.0d5 / RT

      Q = rt_auxvar%immobile(this%im_vec(i))/Cg   ! mol/m^3 bulk

      ! gas rate (scaled by cell volume, gas saturation and porosity)
      g_rate = -this%k(i) * (1.d0 - Q/this%Keq(i)) * &
            material_auxvar%volume * global_auxvar%sat(iphase) * &
            material_auxvar%porosity

      ! immobile rate (1:1 stoichiometry)
      im_rate = -g_rate

      Residual(this%aq_vec(i)) = Residual(this%aq_vec(i)) - g_rate

      Residual(this%im_vec(i) + reaction%offset_immobile) = &
           Residual(this%im_vec(i) + reaction%offset_immobile) - im_rate

    end do

    if (compute_derivative) then
      option%io_buffer = 'NUMERICAL_JACOBIAN must always be used in SUBSURFACE_TRANSPORT' // &
           ' process model due to assumptions in gas reaction sandbox'
      call PrintErrMsg(option)
    endif
  end if

end subroutine GasReact

! ************************************************************************** !

subroutine GasDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Kris Kuhlman
  ! Date: July 2018
  !

  implicit none

  class(reaction_sandbox_gas_type) :: this

  deallocate(this%aq_vec)
  deallocate(this%im_vec)
  deallocate(this%gas_vec)
  deallocate(this%k)
  deallocate(this%Keq)
  deallocate(this%name_vec)

end subroutine GasDestroy

end module Reaction_Sandbox_Gas_class
