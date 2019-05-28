module NW_Transport_Aux_module

  ! this module cannot depend on any other modules besides Option_module
  ! and Matrix_Block_Aux_module
  ! At this point I am violating this rule.
  use Matrix_Block_Aux_module

  use PFLOTRAN_Constants_module
  use petscsys

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscReal, public :: nwt_itol_scaled_res = UNINITIALIZED_DOUBLE
  PetscReal, public :: nwt_itol_rel_update = UNINITIALIZED_DOUBLE
  PetscReal, public :: nwt_min_saturation = 0.d0
  
  type, public :: nw_transport_auxvar_type
    ! total mass as bulk concentration
    PetscReal, pointer :: total_bulk_conc(:)   ! mol-species/m^3-bulk
    ! dissolved 
    PetscReal, pointer :: aqueous_eq_conc(:)   ! mol-species/m^3-liq
    ! sorbed 
    PetscReal, pointer :: sorb_eq_conc(:)      ! mol-species/m^3-sorb
    ! precipitated 
    PetscReal, pointer :: mnrl_eq_conc(:)      ! mol-species/m^3-mnrl 
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
    ! ids of zero rows in local, non-ghosted numbering
    PetscInt, pointer :: zero_rows_local(:) 
    ! ids of zero rows in ghosted numbering
    PetscInt, pointer :: zero_rows_local_ghosted(:)
    ! number of zeroed rows in Jacobian for inactive cells
    PetscInt :: n_zero_rows
    PetscBool :: inactive_cells_exist
    ! nwt_auxvars for local and ghosted grid cells
    type(nw_transport_auxvar_type), pointer :: auxvars(:)
    ! nwt_auxvars for boundary connections
    type(nw_transport_auxvar_type), pointer :: auxvars_bc(:)
    ! nwt_auxvars for source/sinks
    type(nw_transport_auxvar_type), pointer :: auxvars_ss(:)
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
  end type nwt_params_type
  
  type, public :: nwt_print_type
    PetscBool :: molality
    PetscBool :: thingB
  end type nwt_print_type
  
  type, public :: species_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: id
    PetscReal :: molar_weight
    PetscBool :: radioactive
    PetscBool :: print_me
    type(species_type), pointer :: next
  end type species_type
  
  type, public :: nwt_species_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered??????????
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    PetscInt, pointer :: constraint_spec_id(:)
  end type nwt_species_constraint_type
  
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
  
  ! this is the equivalent to reaction_type as in realization%reaction
  type, public :: nw_trans_realization_type
    PetscInt :: offset_auxiliary
    PetscBool :: use_log_formulation
    PetscReal, pointer :: diffusion_coefficient(:,:)
    PetscReal, pointer :: diffusion_activation_energy(:,:)
    character(len=MAXWORDLENGTH), pointer :: species_names(:)
    type(species_type), pointer :: species_list
    PetscBool, pointer :: species_print(:)
    type(radioactive_decay_rxn_type), pointer :: rad_decay_rxn_list
    type(nwt_params_type), pointer :: params
    type(nwt_print_type), pointer :: print 
    PetscBool :: nw_trans_on
  end type nw_trans_realization_type

  interface NWTAuxVarDestroy
    module procedure NWTAuxVarSingleDestroy
    module procedure NWTAuxVarArrayDestroy
  end interface NWTAuxVarDestroy
  
  public :: NWTAuxCreate, NWTSpeciesCreate, NWTRadDecayRxnCreate, &
            NWTRealizCreate, NWTSpeciesConstraintCreate, &
            NWTRead, NWTReadPass2, NWTProcessConstraint, &
            NWTAuxVarInit, NWTAuxVarCopy, NWTAuxVarCopyInitialGuess, &
            NWTAuxDestroy, NWTAuxVarDestroy, NWTAuxVarStrip, &
            NWTSpeciesConstraintDestroy, NWTransDestroy
            
contains

! ************************************************************************** !

function NWTAuxCreate()
  ! 
  ! Allocate and initialize nuclear waste transport auxiliary objects.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys

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
  aux%n_zero_rows = 0    
  nullify(aux%zero_rows_local)  
  nullify(aux%zero_rows_local_ghosted) 
  aux%inactive_cells_exist = PETSC_FALSE

  NWTAuxCreate => aux
  
end function NWTAuxCreate

! ************************************************************************** !

subroutine NWTAuxVarInit(auxvar,nw_trans,option)
  ! 
  ! Initializes the nuclear waste auxiliary object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  !
  
  use Option_module

  implicit none
  
  type(nw_transport_auxvar_type) :: auxvar
  type(nw_trans_realization_type) :: nw_trans
  type(option_type) :: option
  
  PetscInt :: nspecies, nauxiliary, nphase
  
  nspecies = nw_trans%params%nspecies
  nauxiliary = nw_trans%params%nauxiliary
  nphase = nw_trans%params%nphase
  
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

function NWTRealizCreate()
  ! 
  ! Allocates and initializes a NWT realization object
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/19/2019
  ! 
  implicit none
  
  type(nw_trans_realization_type), pointer :: NWTRealizCreate
  
  type(nw_trans_realization_type), pointer :: nwtr 

  allocate(nwtr)
  nwtr%offset_auxiliary = 0
  nwtr%use_log_formulation = PETSC_FALSE
  nullify(nwtr%diffusion_coefficient)
  nullify(nwtr%diffusion_activation_energy)
  nullify(nwtr%species_names)
  nullify(nwtr%species_list)
  nullify(nwtr%species_print)
  nullify(nwtr%rad_decay_rxn_list)
  nwtr%nw_trans_on = PETSC_TRUE
  
  nullify(nwtr%params)
  allocate(nwtr%params)
  nwtr%params%ncomp = 0
  nwtr%params%nphase = 0
  nwtr%params%nspecies = 0
  nwtr%params%nauxiliary = 0
  nwtr%params%calculate_transverse_dispersion = PETSC_FALSE
  nwtr%params%temperature_dependent_diffusion = PETSC_FALSE
  
  nullify(nwtr%print)
  allocate(nwtr%print)
  nwtr%print%molality = PETSC_FALSE
  nwtr%print%thingB = PETSC_FALSE
  
  NWTRealizCreate => nwtr

end function NWTRealizCreate

! ************************************************************************** !

subroutine NWTRead(nw_trans,input,option)
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
  
  type(nw_trans_realization_type), pointer :: nw_trans
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word, parent_name_hold
  character(len=MAXSTRINGLENGTH) :: error_string_base, error_string
  PetscInt :: k, j
  type(species_type), pointer :: new_species, prev_species
  character(len=MAXWORDLENGTH), pointer :: temp_species_names(:)
  character(len=MAXWORDLENGTH), pointer :: temp_species_parents(:)
  type(radioactive_decay_rxn_type), pointer :: new_rad_rxn, prev_rad_rxn
  
  error_string_base = 'SUBSURFACE,NUCLEAR_WASTE_TRANSPORT'
  nullify(prev_species)
  allocate(temp_species_names(50))
  allocate(temp_species_parents(50))
  temp_species_names = ''
  temp_species_parents = ''
  nullify(prev_rad_rxn)
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
        nullify(prev_species)
        error_string = trim(error_string_base) // ',SPECIES'
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          nw_trans%params%nspecies = nw_trans%params%nspecies + 1
          option%ntrandof = nw_trans%params%nspecies
          k = k + 1
          if (k > 50) then
            option%io_buffer = 'More than 50 species are provided in the ' &
                               // trim(error_string) // ', SPECIES block.'
            call PrintErrMsgToDev('if reducing to less than 50 is not &
                                  &an option.',option)
          endif
          new_species => NWTSpeciesCreate()
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'species name',error_string)
          call StringToUpper(word)
          temp_species_names(k) = trim(word)
          new_species%name = trim(word)
          if (.not.associated(nw_trans%species_list)) then
            nw_trans%species_list => new_species
            new_species%id = 1
          endif
          if (associated(prev_species)) then
            prev_species%next => new_species
            new_species%id = prev_species%id + 1
          endif
          prev_species => new_species
          nullify(new_species)
        enddo
        if (k == 0) then
          option%io_buffer = 'ERROR: At least one species &
                              &must be provided in the ' // &
                              trim(error_string) // ' block.'
          call printErrMsg(option)
        endif
        allocate(nw_trans%species_names(k))
        nw_trans%species_names(1:k) = temp_species_names(1:k)
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
            do
              j = j + 1
              if (trim(temp_species_names(j)) == &
                  trim(new_rad_rxn%daughter_name)) then
                temp_species_parents(j) = trim(parent_name_hold)
                exit
              endif
            enddo
          endif

          if (.not.associated(nw_trans%rad_decay_rxn_list)) then
            nw_trans%rad_decay_rxn_list => new_rad_rxn
          endif
          if (associated(prev_rad_rxn)) then
            prev_rad_rxn%next => new_rad_rxn
          endif
          prev_rad_rxn => new_rad_rxn
          nullify(new_rad_rxn)
        enddo
      case('LOG_FORMULATION')
        nw_trans%use_log_formulation = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
            
  enddo
  
  ! assign species_id, parent_id to the rad_rxn objects
  ! check that all radioactive species were listed in the SPECIES block
  call NWTVerifySpecies(nw_trans%species_list,nw_trans%rad_decay_rxn_list, &
                        temp_species_names,temp_species_parents,option)
                        
   deallocate(temp_species_names)
   deallocate(temp_species_parents)
  
end subroutine NWTRead

! ************************************************************************** !

subroutine NWTReadPass2(nw_trans,input,option)
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
  
  type(nw_trans_realization_type), pointer :: nw_trans
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'SUBSURFACE,SUBSURFACE_NUCLEAR_WASTE_TRANSPORT'
  
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
      case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
    
  enddo
  
end subroutine NWTReadPass2

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
  species%molar_weight = 0.d0
  species%radioactive = PETSC_FALSE
  species%print_me = PETSC_FALSE
  nullify(species%next)

  NWTSpeciesCreate => species
  
end function NWTSpeciesCreate

! ************************************************************************** !

function NWTSpeciesConstraintCreate(nw_trans,option)
  ! 
  ! Creates a nuclear waste transport species constraint object
  ! 
  ! Author: Jenn Frederick      
  ! Date: 03/21/2019
  ! 
  use Option_module
  
  implicit none
  
  type(nw_trans_realization_type) :: nw_trans
  type(option_type) :: option
  type(nwt_species_constraint_type), pointer :: NWTSpeciesConstraintCreate

  type(nwt_species_constraint_type), pointer :: constraint
  
  allocate(constraint)
  allocate(constraint%names(nw_trans%params%nspecies))
  constraint%names = ''
  allocate(constraint%constraint_conc(nw_trans%params%nspecies))
  constraint%constraint_conc = 0.d0
  allocate(constraint%constraint_spec_id(nw_trans%params%nspecies))
  constraint%constraint_spec_id = 0

  NWTSpeciesConstraintCreate => constraint

end function NWTSpeciesConstraintCreate

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
      option%io_buffer = 'ERROR: Radioactive species ' // trim(rad_rxn%name) &
                         // ' must also be included in the SPECIES block.'
      call printErrMsg(option)
    endif
    if (.not.daughter_found) then
      option%io_buffer = 'ERROR: Radioactive species ' // &
                         trim(rad_rxn%daughter_name) &
                         // ' must also be included in the SPECIES block.'
      call printErrMsg(option)
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

subroutine NWTProcessConstraint(nw_trans,constraint_name, &
                                nwt_species_constraint,option)
  ! 
  ! Ensures ordering of species is consistant between the nw_trans object
  ! and the constraint object. I don't know what the point of this is.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/22/2019
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use String_module
  use Utility_module
  
  implicit none
  
  type(nw_trans_realization_type), pointer :: nw_trans
  character(len=MAXWORDLENGTH) :: constraint_name
  type(nwt_species_constraint_type), pointer :: nwt_species_constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: ispecies, jspecies
  PetscReal :: constraint_conc(nw_trans%params%nspecies)
  character(len=MAXWORDLENGTH) :: constraint_species_names( &
                                                     nw_trans%params%nspecies)
  
  constraint_conc = 0.d0
  constraint_species_names = ''
  
  do ispecies = 1, nw_trans%params%nspecies
    found = PETSC_FALSE
    do jspecies = 1, nw_trans%params%nspecies
      if (StringCompare(nwt_species_constraint%names(ispecies), &
                        nw_trans%species_names(jspecies),MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
               'Species ' // trim(nwt_species_constraint%names(ispecies)) // &
               ' from CONSTRAINT ' // trim(constraint_name) // &
               ' not found among species.'
      call printErrMsg(option)
    else
      constraint_conc(jspecies) = &
                               nwt_species_constraint%constraint_conc(ispecies)
      constraint_species_names(jspecies) = &
                                         nwt_species_constraint%names(ispecies)
    endif
  enddo
  
  ! place ordered constraint parameters back in original arrays
  nwt_species_constraint%constraint_conc = constraint_conc
  nwt_species_constraint%names = constraint_species_names
  
end subroutine NWTProcessConstraint

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
  call DeallocateArray(aux%zero_rows_local)
  call DeallocateArray(aux%zero_rows_local_ghosted)

  deallocate(aux)
  nullify(aux)

end subroutine NWTAuxDestroy

! ************************************************************************** !

subroutine NWTSpeciesConstraintDestroy(constraint)
  ! 
  ! Deallocates a nuclear waste transport species constraint object
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/21/2019
  ! 

  use Utility_module, only: DeallocateArray
  
  implicit none
  
  type(nwt_species_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_conc)
  call DeallocateArray(constraint%constraint_spec_id)

  deallocate(constraint)
  nullify(constraint)

end subroutine NWTSpeciesConstraintDestroy

! ************************************************************************** !

subroutine NWTransDestroy(nw_trans,option)
  ! 
  ! Deallocates a nuclear waste transport realization object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  ! 

  use Utility_module, only: DeallocateArray
  use Option_module
  
  implicit none
  
  type(nw_trans_realization_type), pointer :: nw_trans
  type(option_type) :: option
  
  type(radioactive_decay_rxn_type), pointer :: rad_decay_rxn,prev_rad_decay_rxn
  type(species_type), pointer :: species, prev_species
  
  if (.not.associated(nw_trans)) return
  
  call DeallocateArray(nw_trans%diffusion_coefficient)
  call DeallocateArray(nw_trans%diffusion_activation_energy)
  call DeallocateArray(nw_trans%species_names)
  call DeallocateArray(nw_trans%species_print)
  
  nullify(nw_trans%params)
  nullify(nw_trans%print)
  
  ! radioactive decay reactions
  rad_decay_rxn => nw_trans%rad_decay_rxn_list
  do
    if (.not.associated(rad_decay_rxn)) exit
    prev_rad_decay_rxn => rad_decay_rxn
    rad_decay_rxn => rad_decay_rxn%next
    nullify(prev_rad_decay_rxn%next)
    deallocate(prev_rad_decay_rxn)  
    nullify(prev_rad_decay_rxn)
  enddo    
  nullify(nw_trans%rad_decay_rxn_list)
  
  ! species
  species => nw_trans%species_list
  do
    if (.not.associated(species)) exit
    prev_species => species
    species => species%next
    nullify(prev_species%next)
    deallocate(prev_species)  
    nullify(prev_species)
  enddo    
  nullify(nw_trans%species_list)
  
  deallocate(nw_trans)
  nullify(nw_trans)

end subroutine NWTransDestroy

! ************************************************************************** !

end module NW_Transport_Aux_module
