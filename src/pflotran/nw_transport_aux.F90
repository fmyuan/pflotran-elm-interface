module NW_Transport_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Reaction_Base_module
  use Matrix_Zeroing_module

  implicit none
  
  private 


  PetscReal, public :: nwt_itol_scaled_res = UNINITIALIZED_DOUBLE
  PetscReal, public :: nwt_itol_rel_update = UNINITIALIZED_DOUBLE
  PetscReal, public :: MIN_LIQ_SAT = 1.0d-20
    
  type, public :: nw_transport_auxvar_type
    ! total mass as bulk concentration
    PetscReal, pointer :: total_bulk_conc(:)   ! mol-species/m^3-bulk
    ! dissolved 
    PetscReal, pointer :: aqueous_eq_conc(:)   ! mol-species/m^3-liq
    ! sorbed 
    PetscReal, pointer :: sorb_eq_conc(:)      ! mol-species/m^3-bulk
    ! precipitated 
    PetscReal, pointer :: mnrl_eq_conc(:)      ! mol-species/m^3-bulk 
    PetscReal, pointer :: mnrl_vol_frac(:)     ! m^3-mnrl/m^3-void 
    ! auxiliary array to store miscellaneous data 
    PetscReal, pointer :: auxiliary_data(:)
    PetscReal, pointer :: mass_balance(:,:)
    PetscReal, pointer :: mass_balance_delta(:,:)
  end type nw_transport_auxvar_type
    
  type, public :: nw_transport_type
    ! number of nwt_auxvars objects for local and ghosted cells
    PetscInt :: num_aux  
    ! number of nwt_auxvars objects for boundary connections
    PetscInt :: num_aux_bc  
    ! number of nwt_auxvars objects for source/sinks
    PetscInt :: num_aux_ss  
    PetscBool :: inactive_cells_exist
    PetscBool :: truncate_output
    ! nwt_auxvars for local and ghosted grid cells
    type(nw_transport_auxvar_type), pointer :: auxvars(:)
    ! nwt_auxvars for boundary connections
    type(nw_transport_auxvar_type), pointer :: auxvars_bc(:)
    ! nwt_auxvars for source/sinks
    type(nw_transport_auxvar_type), pointer :: auxvars_ss(:)
    ! matrix zeroing handling inactive cells
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type nw_transport_type
  
  type, public :: nwt_params_type
    ! number of fluid phases (1=liquid, 2=gas)
    PetscInt :: nphase
    ! number of species
    PetscInt :: nspecies
    ! number of mass components (1=aqueous, 2=precip, 3=sorbed)
    PetscInt :: ncomp
    PetscInt :: nauxiliary
    PetscBool :: calculate_transverse_dispersion
    PetscBool :: temperature_dependent_diffusion
    PetscReal :: truncated_concentration
  end type nwt_params_type
  
  type, public :: nwt_print_type
    PetscBool :: aqueous_eq_conc
    PetscBool :: mnrl_eq_conc
    PetscBool :: mnrl_vol_frac
    PetscBool :: sorb_eq_conc
    PetscBool :: total_bulk_conc
    PetscBool :: all_species
    PetscBool :: all_concs
  end type nwt_print_type
  
  type, public :: species_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: id
    PetscReal :: molar_weight
    PetscReal :: mnrl_molar_density  ! [mol/m^3-mnrl]
    PetscReal :: solubility_limit    ! [mol/m^3-liq]
    PetscReal :: ele_kd              ! [m^3-water/m^3-bulk
    PetscBool :: radioactive
    PetscBool :: print_me
    type(species_type), pointer :: next
  end type species_type
  
  type, public :: radioactive_decay_rxn_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: daughter_name
    character(len=MAXWORDLENGTH) :: parent_name
    PetscInt :: species_id
    PetscInt :: daughter_id
    PetscInt :: parent_id
    PetscReal :: rate_constant
    PetscReal :: rate_constant_parent
    PetscReal :: half_life
    PetscBool :: print_me
    type(radioactive_decay_rxn_type), pointer :: next
  end type radioactive_decay_rxn_type

  type, public :: nwt_species_constraint_type
    ! Any changes here must be incorporated within NWTProcessConstraint(),
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    PetscInt, pointer :: constraint_type(:)
  end type nwt_species_constraint_type
  
  ! this is the equivalent to reaction_rt_type as in realization%reaction
  !TODO(jenn): need to place transport-related variables in nwt_params_type
  !            only reaction-related variables should be in reaction_nw_type
  type, public, extends(reaction_base_type) :: reaction_nw_type
    PetscInt :: offset_auxiliary ! even used?
    PetscReal, pointer :: diffusion_coefficient(:,:) !TODO(jenn): move to nwt_params_type
    PetscReal, pointer :: diffusion_activation_energy(:,:) !TODO(jenn): move to nwt_params_type
    character(len=MAXWORDLENGTH), pointer :: species_names(:)
    type(species_type), pointer :: species_list
    PetscBool, pointer :: species_print(:)
    type(radioactive_decay_rxn_type), pointer :: rad_decay_rxn_list
    type(nwt_params_type), pointer :: params !TODO(jenn): move to nw_transport_type
    type(nwt_print_type), pointer :: print_what 
    PetscBool :: reaction_nw_on
    PetscBool :: truncate_output
  end type reaction_nw_type

  interface NWTAuxVarDestroy
    module procedure NWTAuxVarSingleDestroy
    module procedure NWTAuxVarArrayDestroy
  end interface NWTAuxVarDestroy
  
  public :: NWTAuxCreate, &
            NWTSpeciesCreate, &
            NWTSpeciesConstraintCreate, &
            NWTRadDecayRxnCreate, &
            NWTReactionCreate, &
            NWTReactionCast, &
            NWTRead, &
            NWTReadPass2, &
            NWTSetPlotVariables, &
            NWTAuxVarInit, &
            NWTAuxVarCopy, &
            NWTAuxVarCopyInitialGuess, &
            NWTAuxDestroy, &
            NWTAuxVarDestroy, &
            NWTAuxVarStrip, &
            NWTReactionDestroy
            
            
contains

! ************************************************************************** !

function NWTAuxCreate()
  ! 
  ! Allocate and initialize nuclear waste transport auxiliary objects.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 
  implicit none
  
  type(nw_transport_type), pointer :: NWTAuxCreate
  
  type(nw_transport_type), pointer :: aux

  allocate(aux)  
  aux%num_aux = 0      
  aux%num_aux_bc = 0   
  aux%num_aux_ss = 0     
  nullify(aux%auxvars)      
  nullify(aux%auxvars_bc)   
  nullify(aux%auxvars_ss)   
  nullify(aux%matrix_zeroing)
  aux%inactive_cells_exist = PETSC_FALSE
  aux%truncate_output = PETSC_FALSE

  NWTAuxCreate => aux
  
end function NWTAuxCreate

! ************************************************************************** !

subroutine NWTAuxVarInit(auxvar,reaction_nw,option)
  ! 
  ! Initializes the nuclear waste auxiliary object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  !
  
  use Option_module

  implicit none
  
  type(nw_transport_auxvar_type) :: auxvar
  class(reaction_nw_type) :: reaction_nw
  type(option_type) :: option
  
  PetscInt :: nspecies, nauxiliary, nphase
  
  nspecies = reaction_nw%params%nspecies
  nauxiliary = reaction_nw%params%nauxiliary
  nphase = reaction_nw%params%nphase
    
  allocate(auxvar%total_bulk_conc(nspecies))
  auxvar%total_bulk_conc = 0.d0
  allocate(auxvar%aqueous_eq_conc(nspecies))
  auxvar%aqueous_eq_conc = 0.d0
  allocate(auxvar%sorb_eq_conc(nspecies))
  auxvar%sorb_eq_conc = 0.d0
  allocate(auxvar%mnrl_eq_conc(nspecies))
  auxvar%mnrl_eq_conc = 0.d0
  allocate(auxvar%mnrl_vol_frac(nspecies))
  auxvar%mnrl_vol_frac = 0.d0

  if (nauxiliary > 0) then
    allocate(auxvar%auxiliary_data(nauxiliary))
    auxvar%auxiliary_data = 0.d0
  else
    nullify(auxvar%auxiliary_data)
  endif
  
  if (option%iflag /= 0 .and. option%compute_mass_balance_new) then
    allocate(auxvar%mass_balance(nspecies,nphase))
    auxvar%mass_balance = 0.d0
    allocate(auxvar%mass_balance_delta(nspecies,nphase))
    auxvar%mass_balance_delta = 0.d0
  else
    nullify(auxvar%mass_balance)
    nullify(auxvar%mass_balance_delta)
  endif
  
end subroutine NWTAuxVarInit

! ************************************************************************** !

function NWTReactionCreate()
  ! 
  ! Allocates and initializes a NWT realization object
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/19/2019
  ! 
  implicit none
  
  class(reaction_nw_type), pointer :: NWTReactionCreate
  
  class(reaction_nw_type), pointer :: reaction_nw 

  allocate(reaction_nw)
  call ReactionBaseInit(reaction_nw)

  reaction_nw%offset_auxiliary = 0
  nullify(reaction_nw%diffusion_coefficient)
  nullify(reaction_nw%diffusion_activation_energy)
  nullify(reaction_nw%species_names)
  nullify(reaction_nw%species_list)
  nullify(reaction_nw%species_print)
  nullify(reaction_nw%rad_decay_rxn_list)
  reaction_nw%reaction_nw_on = PETSC_TRUE
  reaction_nw%truncate_output = PETSC_FALSE
  
  nullify(reaction_nw%params)
  allocate(reaction_nw%params)
  reaction_nw%params%ncomp = 0
  reaction_nw%params%nphase = 0
  reaction_nw%params%nspecies = 0
  reaction_nw%params%nauxiliary = 0
  reaction_nw%params%calculate_transverse_dispersion = PETSC_FALSE
  reaction_nw%params%temperature_dependent_diffusion = PETSC_FALSE
  reaction_nw%params%truncated_concentration = UNINITIALIZED_DOUBLE
  
  nullify(reaction_nw%print_what)
  allocate(reaction_nw%print_what)
  reaction_nw%print_what%aqueous_eq_conc = PETSC_FALSE
  reaction_nw%print_what%mnrl_eq_conc = PETSC_FALSE
  reaction_nw%print_what%mnrl_vol_frac = PETSC_FALSE
  reaction_nw%print_what%sorb_eq_conc = PETSC_FALSE
  reaction_nw%print_what%total_bulk_conc = PETSC_FALSE
  reaction_nw%print_what%all_species = PETSC_FALSE
  reaction_nw%print_what%all_concs = PETSC_FALSE
  
  NWTReactionCreate => reaction_nw

end function NWTReactionCreate

! ************************************************************************** !

function NWTReactionCast(reaction_base)
  ! 
  ! Casts a reaction_base type to reaction_nw type if applicable.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/19
  ! 
  implicit none
  
  class(reaction_base_type), pointer :: reaction_base

  class(reaction_nw_type), pointer :: NWTReactionCast

  nullify(NWTReactionCast)
  select type(r=>reaction_base)
    class is(reaction_nw_type)
      NWTReactionCast => r
  end select

end function NWTReactionCast
  
! ************************************************************************** !

subroutine NWTRead(reaction_nw,input,option)
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
 
  implicit none
  
  class(reaction_nw_type), pointer :: reaction_nw
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word, parent_name_hold
  character(len=MAXSTRINGLENGTH) :: error_string_base, error_string
  PetscInt :: k, j
  type(species_type), pointer :: new_species, prev_species
  character(len=MAXWORDLENGTH), pointer :: temp_species_names(:)
  character(len=MAXWORDLENGTH), pointer :: temp_species_parents(:)
  type(radioactive_decay_rxn_type), pointer :: new_rad_rxn, prev_rad_rxn
  
  error_string_base = 'SUBSURFACE,NUCLEAR_WASTE_CHEMISTRY'
  nullify(prev_species)
  allocate(temp_species_names(50))
  allocate(temp_species_parents(50))
  temp_species_names = ''
  temp_species_parents = ''
  nullify(prev_rad_rxn)
  nullify(prev_species)
  k = 0
  
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    select case(trim(keyword))
      case('SPECIES')
        error_string = trim(error_string_base) // ',SPECIES'
        
        reaction_nw%params%nspecies = reaction_nw%params%nspecies + 1
        option%ntrandof = reaction_nw%params%nspecies
        k = k + 1
        if (k > 50) then
          option%io_buffer = 'More than 50 species are provided using ' &
                             // trim(error_string) // ' blocks.'
          call PrintErrMsgToDev(option, 'if reducing to less than 50 is not &
                                &an option.')
        endif
        new_species => NWTSpeciesCreate()
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          call InputReadCard(input,option,keyword)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(keyword)
          
          select case(trim(keyword))
            case('NAME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'species name',error_string)
              call StringToUpper(word)
              temp_species_names(k) = trim(word)
              new_species%name = trim(word)
            case('SOLUBILITY')
              call InputReadDouble(input,option,new_species%solubility_limit)
              call InputErrorMsg(input,option,'species solubility',error_string)
            case('PRECIP_MOLAR_DENSITY','PRECIPITATE_MOLAR_DENSITY')
              call InputReadDouble(input,option,new_species%mnrl_molar_density)
              call InputErrorMsg(input,option,'species mineral molar density', &
                                 error_string)
            case('ELEMENTAL_KD')
              call InputReadDouble(input,option,new_species%ele_kd)
              call InputErrorMsg(input,option,'species elemental Kd', &
                                 error_string)
            case default
              call InputKeywordUnrecognized(input,keyword,error_string,option)
          end select
        enddo
        call InputPopBlock(input,option)
        
        if (new_species%name == '') then
          option%io_buffer = 'NAME not provided in ' // trim(error_string) // &
                             ' block.'
          call PrintErrMsg(option)
        endif
        if (Uninitialized(new_species%solubility_limit)) then
          option%io_buffer = 'SOLUBILITY not provided in ' // &
                             trim(error_string) // ' block for SPECIES ' // &
                             trim(new_species%name) // '.'
          call PrintErrMsg(option)
        endif
        if (Uninitialized(new_species%mnrl_molar_density)) then
          option%io_buffer = 'PRECIPITATE_MOLAR_DENSITY not provided in ' // &
                             trim(error_string) // ' block for SPECIES ' // &
                             trim(new_species%name) // '.'
          call PrintErrMsg(option)
        endif
        if (Uninitialized(new_species%ele_kd)) then
          option%io_buffer = 'ELEMENTAL_KD not provided in ' // &
                             trim(error_string) // ' block for SPECIES ' // &
                             trim(new_species%name) // '.'
          call PrintErrMsg(option)
        endif

        if (.not.associated(reaction_nw%species_list)) then
          reaction_nw%species_list => new_species
          new_species%id = 1
        endif
        if (associated(prev_species)) then
          prev_species%next => new_species
          new_species%id = prev_species%id + 1
        endif
        prev_species => new_species
        nullify(new_species)
        
      case('RADIOACTIVE_DECAY')
        nullify(prev_rad_rxn)
        error_string = trim(error_string_base) // ',RADIOACTIVE_DECAY'
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          new_rad_rxn => NWTRadDecayRxnCreate()
          call InputReadDouble(input,option,new_rad_rxn%rate_constant)
          call InputErrorMsg(input,option,'radioactive species decay rate &
                             &constant',error_string)
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'radioactive species name', &
                             error_string)
          call StringToUpper(word)
          new_rad_rxn%name = trim(word)
          parent_name_hold = trim(word)
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (input%ierr == 0) then ! '->' was read (or anything)
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'radioactive species daughter &
                               &name',error_string)
            call StringToUpper(word)
            new_rad_rxn%daughter_name = trim(word)
            j = 0
            ! record which species was the parent
            do while (j < 50)
              j = j + 1
              if (trim(temp_species_names(j)) == &
                  trim(new_rad_rxn%daughter_name)) then
                temp_species_parents(j) = trim(parent_name_hold)
                exit
              endif
            enddo
          endif

          if (.not.associated(reaction_nw%rad_decay_rxn_list)) then
            reaction_nw%rad_decay_rxn_list => new_rad_rxn
          endif
          if (associated(prev_rad_rxn)) then
            prev_rad_rxn%next => new_rad_rxn
          endif
          prev_rad_rxn => new_rad_rxn
          nullify(new_rad_rxn)
        enddo
      case('LOG_FORMULATION')
        reaction_nw%use_log_formulation = PETSC_TRUE
      case('TRUNCATE_CONCENTRATION')
        error_string = trim(error_string_base) // ',TRUNCATE_CONCENTRATION'
        call InputReadDouble(input,option, &
                             reaction_nw%params%truncated_concentration)
        call InputErrorMsg(input,option,'concentration value',error_string)
      case('OUTPUT')
        call NWTReadOutput(reaction_nw,input,option)
      case default
        call InputKeywordUnrecognized(input,keyword,error_string_base,option)
    end select
            
  enddo
  call InputPopBlock(input,option)
  
  if (k == 0) then
     option%io_buffer = 'ERROR: At least one species must be provided &
                        &in a ' // trim(error_string_base) // ',SPECIES block.'
     call PrintErrMsg(option)
  endif
  allocate(reaction_nw%species_names(k))
  reaction_nw%species_names(1:k) = temp_species_names(1:k)
  
  ! assign species_id, parent_id to the rad_rxn objects
  ! check that all radioactive species were listed in the SPECIES block
  call NWTVerifySpecies(reaction_nw%species_list,reaction_nw%rad_decay_rxn_list, &
                        temp_species_names,temp_species_parents,option)
                        
   deallocate(temp_species_names)
   deallocate(temp_species_parents)
  
end subroutine NWTRead

! ************************************************************************** !

subroutine NWTReadOutput(reaction_nw,input,option)
  ! 
  ! Reads species and concentration types to be printed in output
  ! 
  ! Author: Jenn Frederick
  ! Date: 06/26/2019
  ! 

  use Input_Aux_module
  use String_module  
  use Option_module
  
  implicit none
  
  class(reaction_nw_type) :: reaction_nw
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: word
  
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    error_string = 'NUCLEAR_WASTE_CHEMISTRY,OUTPUT'
    call InputReadCard(input,option,word)  
    call InputErrorMsg(input,option,'keyword',error_string)
                     
    call StringToUpper(word)
    select case(word)
      case('ALL_SPECIES')
        reaction_nw%print_what%all_species = PETSC_TRUE
      case('ALL_CONCENTRATIONS')
        reaction_nw%print_what%all_concs = PETSC_TRUE
      case('TOTAL_BULK_CONCENTRATION')
        reaction_nw%print_what%total_bulk_conc = PETSC_TRUE
      case('AQUEOUS_CONCENTRATION')
        reaction_nw%print_what%aqueous_eq_conc = PETSC_TRUE
      case('MINERAL_CONCENTRATION')
        reaction_nw%print_what%mnrl_eq_conc = PETSC_TRUE
      case('SORBED_CONCENTRATION')
        reaction_nw%print_what%sorb_eq_conc = PETSC_TRUE
      case('MINERAL_VOLUME_FRACTION')
        reaction_nw%print_what%mnrl_vol_frac = PETSC_TRUE
      case('TRUNCATE_OUTPUT')
        reaction_nw%truncate_output = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
    
  enddo
  call InputPopBlock(input,option)
  
  if (reaction_nw%print_what%all_concs) then
    reaction_nw%print_what%total_bulk_conc= PETSC_TRUE
    reaction_nw%print_what%aqueous_eq_conc= PETSC_TRUE
    reaction_nw%print_what%mnrl_eq_conc= PETSC_TRUE
    reaction_nw%print_what%sorb_eq_conc= PETSC_TRUE
  endif

end subroutine NWTReadOutput

! ************************************************************************** !

subroutine NWTReadPass2(reaction_nw,input,option)
  ! 
  ! Reads input file parameters associated with the nuclear waste transport 
  ! process model within the SUBSURFACE block for a second pass.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/20/2019
  !
  use Input_Aux_module
  use String_module
  use Option_module
 
  implicit none
  
  class(reaction_nw_type), pointer :: reaction_nw
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'SUBSURFACE,NUCLEAR_WASTE_CHEMISTRY'
  
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    select case(trim(keyword))
      case('SPECIES')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
        enddo
      case('RADIOACTIVE_DECAY')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
        enddo
      case('LOG_FORMULATION')
      case('TRUNCATE_CONCENTRATION')
      case('OUTPUT')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
        enddo
      !case default
      !  call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
    
  enddo
  call InputPopBlock(input,option)
  
end subroutine NWTReadPass2

! ************************************************************************** !

subroutine NWTSetPlotVariables(list,reaction_nw,option,time_unit)
  ! 
  ! Adds variables to be printed for plotting.
  !
  ! Author: Jenn Frederick
  ! Date: 03/28/2019
  !
  
  use Output_Aux_module
  use Variables_module
  use Option_module
    
  implicit none
  
  type(output_variable_list_type), pointer :: list
  class(reaction_nw_type), pointer :: reaction_nw
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: time_unit
  
  character(len=MAXWORDLENGTH) :: name,  units
  PetscInt :: i
  
  !TODO(jenn) Right now, this assumes ALL_SPECIES are printed by default.
  do i=1,reaction_nw%params%nspecies
    reaction_nw%species_print(i) = PETSC_TRUE
  enddo
  
  if (reaction_nw%print_what%total_bulk_conc) then
    do i=1,reaction_nw%params%nspecies
      if (reaction_nw%species_print(i)) then
        name = 'Total Bulk Conc. ' // trim(reaction_nw%species_names(i))
        units = 'mol/m^3-bulk'
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     TOTAL_BULK_CONC,i)
      endif
    enddo
  endif 
    
  if (reaction_nw%print_what%aqueous_eq_conc) then
    do i=1,reaction_nw%params%nspecies
      if (reaction_nw%species_print(i)) then
        name = 'Aq. Conc. ' // trim(reaction_nw%species_names(i))
        units = 'mol/m^3-liq'
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     AQUEOUS_EQ_CONC,i)
      endif
    enddo
  endif 
  
  if (reaction_nw%print_what%mnrl_eq_conc) then
    do i=1,reaction_nw%params%nspecies
      if (reaction_nw%species_print(i)) then
        name = 'Mnrl. Conc. ' // trim(reaction_nw%species_names(i))
        units = 'mol/m^3-bulk'
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     MNRL_EQ_CONC,i)
      endif
    enddo
  endif 
  
  if (reaction_nw%print_what%sorb_eq_conc) then
    do i=1,reaction_nw%params%nspecies
      if (reaction_nw%species_print(i)) then
        name = 'Sorb. Conc. ' // trim(reaction_nw%species_names(i))
        units = 'mol/m^3-bulk'
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     SORB_EQ_CONC,i)
      endif
    enddo
  endif
  
  if (reaction_nw%print_what%mnrl_vol_frac) then
    do i=1,reaction_nw%params%nspecies
      if (reaction_nw%species_print(i)) then
        name = 'Mnrl. Vol. Frac. ' // trim(reaction_nw%species_names(i))
        units = 'm^3-mnrl/m^3-void'
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     MNRL_VOLUME_FRACTION,i)
      endif
    enddo
  endif
  
end subroutine NWTSetPlotVariables

! ************************************************************************** !

function NWTSpeciesCreate()
  ! 
  ! Allocate and initialize a species object.
  ! 
  ! Author: Jenn Frederick      
  ! Date: 03/11/2019
  ! 
    
  implicit none
  
  type(species_type), pointer :: NWTSpeciesCreate
  
  type(species_type), pointer :: species

  allocate(species) 
  species%id = 0 
  species%name = ''
  species%molar_weight = UNINITIALIZED_DOUBLE
  species%mnrl_molar_density = UNINITIALIZED_DOUBLE
  species%solubility_limit = UNINITIALIZED_DOUBLE 
  species%ele_kd = UNINITIALIZED_DOUBLE 
  species%radioactive = PETSC_FALSE
  species%print_me = PETSC_FALSE
  nullify(species%next)

  NWTSpeciesCreate => species
  
end function NWTSpeciesCreate

! ************************************************************************** !

function NWTRadDecayRxnCreate()
  ! 
  ! Allocate and initialize a radioactive decay reaction object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 
  
  implicit none
    
  type(radioactive_decay_rxn_type), pointer :: NWTRadDecayRxnCreate

  type(radioactive_decay_rxn_type), pointer :: rxn
  
  allocate(rxn)
  rxn%daughter_id = 0
  rxn%species_id = 0
  rxn%parent_id = 0
  rxn%name = ''
  rxn%daughter_name = ''
  rxn%parent_name = ''
  rxn%rate_constant = 0.d0
  rxn%rate_constant_parent = 0.d0
  rxn%half_life = 0.d0
  rxn%print_me = PETSC_FALSE
  nullify(rxn%next)
  
  NWTRadDecayRxnCreate => rxn
  
end function NWTRadDecayRxnCreate

! ************************************************************************** !

function NWTSpeciesConstraintCreate(reaction_nw,option)
  ! 
  ! Creates a nuclear waste transport species constraint object
  ! 
  ! Author: Jenn Frederick      
  ! Date: 03/21/2019
  ! 
  use Option_module

  implicit none

  class(reaction_nw_type) :: reaction_nw
  type(option_type) :: option
  type(nwt_species_constraint_type), pointer :: NWTSpeciesConstraintCreate

  type(nwt_species_constraint_type), pointer :: constraint

  allocate(constraint)
  allocate(constraint%names(reaction_nw%params%nspecies))
  constraint%names = ''
  allocate(constraint%constraint_conc(reaction_nw%params%nspecies))
  constraint%constraint_conc = 0.d0
  allocate(constraint%constraint_type(reaction_nw%params%nspecies))
  constraint%constraint_type = 0

  NWTSpeciesConstraintCreate => constraint

end function NWTSpeciesConstraintCreate

! ************************************************************************** !

subroutine NWTVerifySpecies(species_list,rad_decay_rxn_list,species_names, &
                            parent_names,option)
  ! 
  ! Assigns species_id, parent_id to the rad_rxn objects after checking that
  ! all radioactive species were listed in the SPECIES block.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/09/2019
  !
  
  use Option_module

  implicit none
  
  type(species_type), pointer :: species_list
  type(radioactive_decay_rxn_type), pointer :: rad_decay_rxn_list
  character(len=MAXWORDLENGTH), pointer :: species_names(:)
  character(len=MAXWORDLENGTH), pointer :: parent_names(:) 
  type(option_type) :: option
  
  type(species_type), pointer :: species
  type(radioactive_decay_rxn_type), pointer :: rad_rxn, cur_rad_rxn
  PetscInt :: k 
  
  PetscBool :: parent_found, daughter_found
  
  rad_rxn => rad_decay_rxn_list
  do
    if (.not.associated(rad_rxn)) exit
    ! Check if the parent and daughter radioactive species are listed within 
    ! the SPECIES block. 
    species => species_list
    do
      if (.not.associated(species)) exit
      if (trim(rad_rxn%name) == trim(species%name)) then 
        parent_found = PETSC_TRUE
        rad_rxn%species_id = species%id
        species%radioactive = PETSC_TRUE
      endif
      if (len(trim(rad_rxn%daughter_name)) > 0) then
        if (trim(rad_rxn%daughter_name) == trim(species%name)) then 
          daughter_found = PETSC_TRUE
          rad_rxn%daughter_id = species%id
        endif
      else
        daughter_found = PETSC_TRUE
      endif
                  
      species => species%next
    enddo
    if (.not.parent_found) then
      option%io_buffer = 'Radioactive species ' // trim(rad_rxn%name) &
                         // ' must also be included using a SPECIES block.'
      call PrintErrMsg(option)
    endif
    if (.not.daughter_found) then
      option%io_buffer = 'Radioactive species ' // trim(rad_rxn%daughter_name) &
                         // ' must also be included using a SPECIES block.'
      call PrintErrMsg(option)
    endif
    
    ! assign the parent information
    k = 0
    do
      k = k + 1
      if (trim(species_names(k)) == trim(rad_rxn%name)) then
        rad_rxn%parent_name = trim(parent_names(k))
        exit
      endif
    enddo
    species => species_list
    do
      if (.not.associated(species)) exit
      if (trim(rad_rxn%parent_name) == trim(species%name)) then
        rad_rxn%parent_id = species%id
        exit
      endif
      
      species => species%next
    enddo
    
    daughter_found = PETSC_FALSE
    parent_found = PETSC_FALSE
    
    rad_rxn => rad_rxn%next
  enddo
  
  cur_rad_rxn => rad_decay_rxn_list
  do
    if (.not.associated(cur_rad_rxn)) exit
    ! check if current radioactive species has a parent
    if (cur_rad_rxn%parent_id > 0) then
      rad_rxn => rad_decay_rxn_list
      do
        if (.not.associated(rad_rxn)) exit
        ! assign the parent's decay rate to rate_constant_parent
        if (rad_rxn%species_id == cur_rad_rxn%parent_id) then
          cur_rad_rxn%rate_constant_parent = rad_rxn%rate_constant
          exit
        endif
        rad_rxn => rad_rxn%next
      enddo
    endif
    cur_rad_rxn => cur_rad_rxn%next
  enddo
    
end subroutine NWTVerifySpecies

! ************************************************************************** !

subroutine NWTAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies the nuclear waste transport auxiliary object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 
  
  use Option_module
  
  implicit none
  
  type(nw_transport_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option
  
end subroutine NWTAuxVarCopy

! ************************************************************************** !

subroutine NWTAuxVarCopyInitialGuess(auxvar,auxvar2,option)
  ! 
  ! Copies the molality in nwt_auxvar to serve as an initial guess when
  ! equilibrating constraints.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 

  use Option_module

  implicit none
  
  type(nw_transport_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option  
  
  auxvar%total_bulk_conc = auxvar2%total_bulk_conc
  
end subroutine NWTAuxVarCopyInitialGuess

! ************************************************************************** !

subroutine NWTAuxVarSingleDestroy(auxvar)
  ! 
  ! Deallocates the nuclear waste transport auxiliary object (single)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 

  implicit none

  type(nw_transport_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call NWTAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)  

end subroutine NWTAuxVarSingleDestroy

! ************************************************************************** !

subroutine NWTAuxVarArrayDestroy(auxvars)
  ! 
  ! Deallocates the nuclear waste transport auxiliary object (array)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  !  

  implicit none

  type(nw_transport_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call NWTAuxVarStrip(auxvars(iaux))
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine NWTAuxVarArrayDestroy

! ************************************************************************** !

subroutine NWTAuxVarStrip(auxvar)
  ! 
  ! Deallocates all members of a single nuclear waste transport aux object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 

  use Utility_module, only: DeallocateArray

  implicit none

  type(nw_transport_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%total_bulk_conc)
  call DeallocateArray(auxvar%aqueous_eq_conc)
  call DeallocateArray(auxvar%sorb_eq_conc)
  call DeallocateArray(auxvar%mnrl_eq_conc)
  call DeallocateArray(auxvar%mnrl_vol_frac)
  call DeallocateArray(auxvar%auxiliary_data)
  
end subroutine NWTAuxVarStrip

! ************************************************************************** !

subroutine NWTAuxDestroy(aux)
  ! 
  ! Deallocates a nuclear waste transport auxiliary object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 

  use Utility_module, only: DeallocateArray
  
  implicit none

  type(nw_transport_type), pointer :: aux
  
  if (.not.associated(aux)) return
  
  call NWTAuxVarDestroy(aux%auxvars)
  call NWTAuxVarDestroy(aux%auxvars_bc)
  call NWTAuxVarDestroy(aux%auxvars_ss)
  call MatrixZeroingDestroy(aux%matrix_zeroing)

  deallocate(aux)
  nullify(aux)

end subroutine NWTAuxDestroy

! ************************************************************************** !

subroutine NWTReactionDestroy(reaction_nw,option)
  ! 
  ! Deallocates a nuclear waste transport realization object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  ! 

  use Utility_module, only: DeallocateArray
  use Option_module
  
  implicit none
  
  class(reaction_nw_type), pointer :: reaction_nw
  type(option_type) :: option
  
  type(radioactive_decay_rxn_type), pointer :: rad_decay_rxn,prev_rad_decay_rxn
  type(species_type), pointer :: species, prev_species
  
  if (.not.associated(reaction_nw)) return

  call ReactionBaseStrip(reaction_nw)
  
  call DeallocateArray(reaction_nw%diffusion_coefficient)
  call DeallocateArray(reaction_nw%diffusion_activation_energy)
  call DeallocateArray(reaction_nw%species_names)
  call DeallocateArray(reaction_nw%species_print)
  
  nullify(reaction_nw%params)
  nullify(reaction_nw%print_what)
  
  ! radioactive decay reactions
  rad_decay_rxn => reaction_nw%rad_decay_rxn_list
  do
    if (.not.associated(rad_decay_rxn)) exit
    prev_rad_decay_rxn => rad_decay_rxn
    rad_decay_rxn => rad_decay_rxn%next
    nullify(prev_rad_decay_rxn%next)
    deallocate(prev_rad_decay_rxn)  
    nullify(prev_rad_decay_rxn)
  enddo    
  nullify(reaction_nw%rad_decay_rxn_list)
  
  ! species
  species => reaction_nw%species_list
  do
    if (.not.associated(species)) exit
    prev_species => species
    species => species%next
    nullify(prev_species%next)
    deallocate(prev_species)  
    nullify(prev_species)
  enddo    
  nullify(reaction_nw%species_list)
  
  deallocate(reaction_nw)
  nullify(reaction_nw)

end subroutine NWTReactionDestroy

! ************************************************************************** !

end module NW_Transport_Aux_module
