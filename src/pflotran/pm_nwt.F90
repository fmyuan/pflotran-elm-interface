module PM_NWT_class
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class 
  use Realization_Subsurface_class
  use Communicator_Base_module  
  use Option_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
  type, public :: species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: a0  ! don't know what this is
    PetscReal :: molar_weight
    PetscReal :: Z  ! don't know what this is
    PetscBool :: print_me
    type(species_type), pointer :: next
  end type species_type
  
  type, public :: radioactive_decay_rxn_type
    PetscInt :: id
    character(len=MAXSTRINGLENGTH) :: reaction
    PetscReal :: rate_constant
    PetscReal :: half_life
    PetscBool :: print_me
    type(radioactive_decay_rxn_type), pointer :: next
  end type radioactive_decay_rxn_type
  
  type, public :: nw_transport_rxn_type
    PetscInt :: offset_auxiliary
    PetscBool :: use_log_formulation
    PetscReal, pointer :: diffusion_coefficient(:,:)
    PetscReal, pointer :: diffusion_activation_energy(:,:)
    character(len=MAXWORDLENGTH), pointer :: species_names(:)
    type(species_type), pointer :: species_list
    PetscBool, pointer :: species_print(:)
    type(radioactive_decay_rxn_type), pointer :: rad_decay_rxn_list
  end type nw_transport_rxn_type
  
  type, public :: pm_nwt_controls_type
    ! governs the size of subsequent time steps
    PetscReal, pointer :: max_concentration_change(:)
    PetscReal, pointer :: max_volfrac_change(:)
    PetscReal :: volfrac_change_governor
    PetscReal :: cfl_governor
  end type pm_nwt_controls_type
  
  type, public :: pm_nwt_params_type
    PetscInt :: nphase
    PetscInt :: ncomp
    PetscInt :: nsorb
    PetscInt :: nmnrl
    PetscInt :: nauxiliary
#ifdef OS_STATISTICS
! use PetscReal for large counts
    PetscInt :: newton_call_count
    PetscReal :: sum_newton_call_count
    PetscInt :: newton_iterations
    PetscReal :: sum_newton_iterations
    PetscInt :: max_newton_iterations
    PetscInt :: overall_max_newton_iterations
#endif    
    PetscReal :: newton_inf_rel_update_tol
    PetscReal :: newton_inf_scaled_res_tol
    PetscBool :: check_post_converged
    PetscBool :: calculate_transverse_dispersion
    PetscBool :: temperature_dependent_diffusion
  end type pm_nwt_params_type
  
  type, public, extends(pm_base_type) :: pm_nwt_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    ! local variables
    PetscBool :: check_post_convergence
    PetscBool :: temperature_dependent_diffusion
    type(pm_nwt_controls_type), pointer :: controls
    type(pm_nwt_params_type), pointer :: params
    type(nw_transport_rxn_type), pointer :: rxn 
  contains
    procedure, public :: Setup => PMNWTSetup 
    procedure, public :: ReadSimulationBlock => PMNWTReadSimulationBlock
    procedure, public :: ReadPMBlock => PMNWTReadPMBlock
    procedure, public :: SetRealization => PMNWTSetRealization   
  end type pm_nwt_type
  
  public :: PMNWTCreate, PMNWTSpeciesCreate, PMNWTRadDecayRxnCreate
  
  
contains

! ************************************************************************** !

function PMNWTCreate()
  ! 
  ! Creates the nuclear waste transport process model shell.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  ! 

  implicit none
  
  class(pm_nwt_type), pointer :: PMNWTCreate

  class(pm_nwt_type), pointer :: nwt_pm
    
  allocate(nwt_pm)
  nullify(nwt_pm%option)
  nullify(nwt_pm%output_option)
  nullify(nwt_pm%realization)
  nullify(nwt_pm%comm1)
  
  ! local variables
  nwt_pm%check_post_convergence = PETSC_FALSE
  nwt_pm%temperature_dependent_diffusion = PETSC_FALSE
  
  allocate(nwt_pm%controls)
  nullify(nwt_pm%controls%max_concentration_change)
  nullify(nwt_pm%controls%max_volfrac_change)
  nwt_pm%controls%volfrac_change_governor = 1.d0
  nwt_pm%controls%cfl_governor = UNINITIALIZED_DOUBLE
  
  allocate(nwt_pm%rxn)
  nwt_pm%rxn%offset_auxiliary = 0
  nwt_pm%rxn%use_log_formulation = PETSC_FALSE
  nullify(nwt_pm%rxn%diffusion_coefficient)
  ! later on in setup:
  ! allocate(nwt_rxn%diffusion_coefficient(ncomp,nphase))
  ! nwt_rxn%diffusion_coefficient = 1.d-9
  nullify(nwt_pm%rxn%diffusion_activation_energy)
  ! later on in setup:
  ! allocate(nwt_rxn%diffusion_activation_energy(ncomp,nphase))
  ! nwt_rxn%diffusion_activation_energy = 0.d0
  nullify(nwt_pm%rxn%species_names)
  ! later on in setup:
  ! allocate(nwt_rxn%species_names(ncomp))
  ! nwt_rxn%species_names = ''
  nullify(nwt_pm%rxn%species_list)
  nullify(nwt_pm%rxn%species_print)
  ! later on in setup:
  ! allocate(nwt_rxn%species_print(ncomp))
  ! nwt_rxn%species_print = PETSC_FALSE
  nullify(nwt_pm%rxn%rad_decay_rxn_list)
  
  allocate(nwt_pm%params)
  nwt_pm%params%ncomp = 0
  nwt_pm%params%nphase = 1  ! For WIPP, we always assume liquid phase only
  nwt_pm%params%nsorb = 0
  nwt_pm%params%nmnrl = 0
  nwt_pm%params%nauxiliary = 0
  nwt_pm%params%calculate_transverse_dispersion = PETSC_FALSE
  nwt_pm%params%temperature_dependent_diffusion = PETSC_FALSE
  nwt_pm%params%check_post_converged = PETSC_FALSE
  nwt_pm%params%newton_inf_rel_update_tol = UNINITIALIZED_DOUBLE
  nwt_pm%params%newton_inf_scaled_res_tol = UNINITIALIZED_DOUBLE
#ifdef OS_STATISTICS
  nwt_pm%params%newton_call_count = 0
  nwt_pm%params%sum_newton_call_count = 0.d0
  nwt_pm%params%newton_iterations = 0
  nwt_pm%params%sum_newton_iterations = 0.d0
  nwt_pm%params%max_newton_iterations = 0
  nwt_pm%params%overall_max_newton_iterations = 0
#endif

  call PMBaseInit(nwt_pm)
  nwt_pm%name = 'Nuclear Waste Transport'
  nwt_pm%header = 'NUCLEAR WASTE TRANSPORT'
  
  PMNWTCreate => nwt_pm
  
end function PMNWTCreate

! ************************************************************************** !

function PMNWTSpeciesCreate()
  ! 
  ! Allocate and initialize a species object.
  ! 
  ! Author: Jenn Frederick      
  ! Date: 03/11/2019
  ! 
    
  implicit none
  
  type(species_type), pointer :: PMNWTSpeciesCreate
  
  type(species_type), pointer :: species

  allocate(species) 
  species%id = 0 
  species%name = ''
  species%a0 = 0.d0
  species%molar_weight = 0.d0
  species%Z = 0.d0
  species%print_me = PETSC_FALSE
  nullify(species%next)

  PMNWTSpeciesCreate => species
  
end function PMNWTSpeciesCreate

! ************************************************************************** !

function PMNWTRadDecayRxnCreate()
  ! 
  ! Allocate and initialize a radioactive decay reaction object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 
  
  implicit none
    
  type(radioactive_decay_rxn_type), pointer :: PMNWTRadDecayRxnCreate

  type(radioactive_decay_rxn_type), pointer :: rxn
  
  allocate(rxn)
  rxn%id = 0
  rxn%reaction = ''
  rxn%rate_constant = 0.d0
  rxn%half_life = 0.d0
  rxn%print_me = PETSC_FALSE
  nullify(rxn%next)
  
  PMNWTRadDecayRxnCreate => rxn
  
end function PMNWTRadDecayRxnCreate

! ************************************************************************** !
  
subroutine PMNWTSetup(this)
  ! 
  ! Initializes variables associated with nuclear waste transport.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  ! 
  
  use NW_Transport_Aux_module

  implicit none
  
  class(pm_nwt_type) :: this
  
  type(nw_transport_param_type), pointer :: nwt_parameter
  
  ! jenn:todo get NWT in patch%aux pointed to correct place?
  nwt_parameter => this%realization%patch%aux%NWT%nwt_parameter
  
  ! pass down flags from PMNWT class
  nwt_parameter%temperature_dependent_diffusion = &
    this%temperature_dependent_diffusion
  
  ! set the communicator
  this%comm1 => this%realization%comm1
  
  allocate(this%controls%max_concentration_change(nwt_parameter%ncomp))
  allocate(this%controls%max_volfrac_change(nwt_parameter%ncomp))

end subroutine PMNWTSetup

! ************************************************************************** !

subroutine PMNWTReadSimulationBlock(this,input)
  ! 
  ! Reads input file parameters associated with the nuclear waste transport 
  ! process model in the SIMULATION block.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use NW_Transport_Aux_module
 
  implicit none
  
  class(pm_nwt_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option
  
  error_string = 'NUCLEAR_WASTE_TRANSPORT OPTIONS'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    found = PETSC_FALSE
    call PMBaseReadSelectCase(this,input,keyword,found,error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('GLOBAL_IMPLICIT','OPERATOR_SPLIT','OPERATOR_SPLITTING')
      case('MAX_VOLUME_FRACTION_CHANGE')
        call InputReadDouble(input,option,this%controls%volfrac_change_governor)
        call InputErrorMsg(input,option,'MAX_VOLUME_FRACTION_CHANGE', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
      case('ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,nwt_itol_rel_update)
        call InputErrorMsg(input,option,'ITOL_RELATIVE_UPDATE', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
        this%check_post_convergence = PETSC_TRUE
      case('MINIMUM_SATURATION')
        call InputReadDouble(input,option,nwt_min_saturation)
        call InputErrorMsg(input,option,'MINIMUM_SATURATION', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
      case('NUMERICAL_JACOBIAN')
        option%transport%numerical_derivatives = PETSC_TRUE
      case('TEMPERATURE_DEPENDENT_DIFFUSION')
        this%temperature_dependent_diffusion = PETSC_TRUE
      case('MAX_CFL')
        call InputReadDouble(input,option,this%controls%cfl_governor)
        call InputErrorMsg(input,option,'MAX_CFL', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
      case('MULTIPLE_CONTINUUM')
        option%use_mc = PETSC_TRUE          
      case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  enddo
  
end subroutine PMNWTReadSimulationBlock

! ************************************************************************** !

subroutine PMNWTReadPMBlock(this,input)
  ! 
  ! Reads input file parameters associated with the nuclear waste transport 
  ! process model within the SUBSURFACE block.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use NW_Transport_Aux_module
 
  implicit none
  
  class(pm_nwt_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscInt :: k
  character(len=MAXWORDLENGTH), pointer :: temp_species_names(:)

  option => this%option
  
  error_string = 'SUBSURFACE,NUCLEAR_WASTE_TRANSPORT'
  allocate(temp_species_names(50))
  temp_species_names = ''
  k = 0
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    select case(trim(keyword))
      case('SPECIES')
        error_string = trim(error_string) // ',SPECIES'
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          this%params%ncomp = this%params%ncomp + 1
          option%ntrandof = this%params%ncomp
          k = k + 1
          if (k > 50) then
            option%io_buffer = 'More than 50 species are provided in the ' &
                               // trim(error_string) // ', SPECIES block.'
            call PrintErrMsgToDev('if reducing to less than 50 is not &
                                  &an option.',option)
          endif
          
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'species name',error_string)
          call StringToUpper(word)
          temp_species_names(k) = trim(word)
        enddo
        if (k == 0) then
          option%io_buffer = 'ERROR: At least one radionuclide species &
                              &must be provided in the ' // &
                              trim(error_string) // ' block.'
          call printErrMsg(option)
        endif
        allocate(this%rxn%species_names(k))
        this%rxn%species_names(1:k) = temp_species_names(1:k)
        deallocate(temp_species_names)
      case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
    
  enddo
  
end subroutine PMNWTReadPMBlock

! ************************************************************************** !

subroutine PMNWTSetRealization(this,realization)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  ! 

  use Realization_Subsurface_class  

  implicit none
  
  class(pm_nwt_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization
  
  if (this%rxn%use_log_formulation) then
    this%solution_vec = realization%field%tran_log_xx
  else
    this%solution_vec = realization%field%tran_xx
  endif
  this%residual_vec = realization%field%tran_r
  
end subroutine PMNWTSetRealization

! ************************************************************************** !
  
! jenn:todo Remember to deallocate/destroy all PMNWT pointers and arrays

end module PM_NWT_class