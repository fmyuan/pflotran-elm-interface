module NW_Transport_Aux_module

  ! this module cannot depend on any other modules besides Option_module
  ! and Matrix_Block_Aux_module
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
    ! molality
    PetscReal, pointer :: molality(:)     ! mol/kg water
    ! sorbed totals
    PetscReal, pointer :: total_sorb_eq(:)    ! mol/m^3 bulk
    PetscReal, pointer :: dtotal_sorb_eq(:,:) ! kg water/m^3 bulk
    ! mineral reactions
    PetscReal, pointer :: mnrl_volfrac0(:)
    PetscReal, pointer :: mnrl_volfrac(:)
    PetscReal, pointer :: mnrl_area0(:)
    PetscReal, pointer :: mnrl_area(:)
    PetscReal, pointer :: mnrl_rate(:)
    ! auxiliary array to store miscellaneous data (e.g. reaction rates, 
    ! cumulative mass, etc.
    PetscReal, pointer :: auxiliary_data(:)
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
    !type(nw_transport_param_type), pointer :: nwt_parameter
    ! nwt_auxvars for local and ghosted grid cells
    type(nw_transport_auxvar_type), pointer :: auxvars(:)
    ! nwt_auxvars for boundary connections
    type(nw_transport_auxvar_type), pointer :: auxvars_bc(:)
    ! nwt_auxvars for source/sinks
    type(nw_transport_auxvar_type), pointer :: auxvars_ss(:)
  end type nw_transport_type
  
  type, public :: nwt_params_type
    PetscInt :: nphase
    PetscInt :: ncomp
    PetscInt :: nsorb
    PetscInt :: nmnrl
    PetscInt :: nauxiliary
    PetscBool :: calculate_transverse_dispersion
    PetscBool :: temperature_dependent_diffusion
  end type nwt_params_type
  
  type, public :: species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: a0  ! don't know what this is
    PetscReal :: molar_weight
    PetscReal :: Z  ! don't know what this is
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
    PetscInt :: id
    character(len=MAXSTRINGLENGTH) :: reaction
    PetscReal :: rate_constant
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
    PetscBool :: nw_trans_on
  end type nw_trans_realization_type

  interface NWTAuxVarDestroy
    module procedure NWTAuxVarSingleDestroy
    module procedure NWTAuxVarArrayDestroy
  end interface NWTAuxVarDestroy
  
  public :: NWTAuxCreate, NWTSpeciesCreate, NWTRadDecayRxnCreate, &
            NWTRealizCreate, NWTSpeciesConstraintCreate, &
            NWTRead, NWTReadPass2, &
            NWTAuxVarInit, NWTAuxVarCopy, NWTAuxVarCopyInitialGuess, &
            NWTAuxDestroy, NWTAuxVarDestroy, NWTAuxVarStrip, &
            NWTSpeciesConstraintDestroy
            
contains

! ************************************************************************** !

function NWTAuxCreate(ncomp,nphase)
  ! 
  ! Allocate and initialize nuclear waste transport auxiliary objects.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none
  
  ! jenn:todo Why do we need these here? ncomp, nphase
  PetscInt :: ncomp
  PetscInt :: nphase
  
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
  
  PetscInt :: ncomp, nmnrl, nsorb, nauxiliary
  
  ncomp = nw_trans%params%ncomp
  nmnrl = nw_trans%params%nmnrl
  nsorb = nw_trans%params%nsorb
  nauxiliary = nw_trans%params%nauxiliary
  
  allocate(auxvar%molality(ncomp))
  auxvar%molality = 0.d0
  
  if (nsorb > 0) then
    allocate(auxvar%total_sorb_eq(ncomp))
    auxvar%total_sorb_eq = 0.d0
    allocate(auxvar%dtotal_sorb_eq(ncomp,ncomp))
    auxvar%dtotal_sorb_eq = 0.d0
  else
    nullify(auxvar%total_sorb_eq)
    nullify(auxvar%dtotal_sorb_eq)
  endif
  
  if (nmnrl > 0) then
    allocate(auxvar%mnrl_volfrac0(nmnrl))
    auxvar%mnrl_volfrac0 = 0.d0
    allocate(auxvar%mnrl_volfrac(nmnrl))
    auxvar%mnrl_volfrac = 0.d0
    allocate(auxvar%mnrl_area0(nmnrl))
    auxvar%mnrl_area0 = 0.d0
    allocate(auxvar%mnrl_area(nmnrl))
    auxvar%mnrl_area = 0.d0
    allocate(auxvar%mnrl_rate(nmnrl))
    auxvar%mnrl_rate = 0.d0
  else
    nullify(auxvar%mnrl_volfrac0)
    nullify(auxvar%mnrl_volfrac)
    nullify(auxvar%mnrl_area0)
    nullify(auxvar%mnrl_area)
    nullify(auxvar%mnrl_rate)
  endif
  
  if (nauxiliary > 0) then
    allocate(auxvar%auxiliary_data(nauxiliary))
    auxvar%auxiliary_data = 0.d0
  else
    nullify(auxvar%auxiliary_data)
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
  nwtr%params%nphase = 1  ! For WIPP, we always assume liquid phase only
  nwtr%params%nsorb = 0
  nwtr%params%nmnrl = 0
  nwtr%params%nauxiliary = 0
  nwtr%params%calculate_transverse_dispersion = PETSC_FALSE
  nwtr%params%temperature_dependent_diffusion = PETSC_FALSE
  
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
  
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: k
  type(species_type), pointer :: new_species, prev_species
  character(len=MAXWORDLENGTH), pointer :: temp_species_names(:)
  
  error_string = 'SUBSURFACE,NUCLEAR_WASTE_TRANSPORT'
  nullify(prev_species)
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
        nullify(prev_species)
        error_string = trim(error_string) // ',SPECIES'
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          nw_trans%params%ncomp = nw_trans%params%ncomp + 1
          option%ntrandof = nw_trans%params%ncomp
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
          option%io_buffer = 'ERROR: At least one radionuclide species &
                              &must be provided in the ' // &
                              trim(error_string) // ' block.'
          call printErrMsg(option)
        endif
        allocate(nw_trans%species_names(k))
        nw_trans%species_names(1:k) = temp_species_names(1:k)
        deallocate(temp_species_names)
      case('LOG_FORMULATION')
        nw_trans%use_log_formulation = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
    
  enddo
  
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
  species%a0 = 0.d0
  species%molar_weight = 0.d0
  species%Z = 0.d0
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
  allocate(constraint%names(nw_trans%params%ncomp))
  constraint%names = ''
  allocate(constraint%constraint_conc(nw_trans%params%ncomp))
  constraint%constraint_conc = 0.d0
  allocate(constraint%constraint_spec_id(nw_trans%params%ncomp))
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
  rxn%id = 0
  rxn%reaction = ''
  rxn%rate_constant = 0.d0
  rxn%half_life = 0.d0
  rxn%print_me = PETSC_FALSE
  nullify(rxn%next)
  
  NWTRadDecayRxnCreate => rxn
  
end function NWTRadDecayRxnCreate

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
  
  auxvar%molality = auxvar2%molality
  
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
  
  call DeallocateArray(auxvar%molality)
  call DeallocateArray(auxvar%total_sorb_eq)
  call DeallocateArray(auxvar%dtotal_sorb_eq)
    
  call DeallocateArray(auxvar%mnrl_volfrac0)
  call DeallocateArray(auxvar%mnrl_volfrac)
  call DeallocateArray(auxvar%mnrl_area0)
  call DeallocateArray(auxvar%mnrl_area)
  call DeallocateArray(auxvar%mnrl_rate)
    
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

end module NW_Transport_Aux_module