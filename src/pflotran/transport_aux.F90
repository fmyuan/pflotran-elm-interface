module Transport_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  ! this module cannot depend on any other modules besides Option_module
  ! and Matrix_Block_Aux_module
  use Matrix_Block_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  PetscReal, public :: tran_itol_scaled_res = UNINITIALIZED_DOUBLE
  PetscReal, public :: tran_itol_rel_update = UNINITIALIZED_DOUBLE
  PetscReal, public :: tran_min_saturation = 0.d0

  type, public :: tran_species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscInt  :: nphase   ! default 2, ie gaseous and aqueous
    PetscReal :: molar_weight
    PetscReal :: molar_gas_volume
    PetscBool :: is_reactive
    PetscReal, pointer :: diffusion_coefficient(:)
    PetscReal, pointer :: diffusion_activation_energy(:)
    PetscReal, pointer :: fugacoeff(:)

    type(tran_species_type), pointer :: next
  end type tran_species_type

  type, public :: tran_species_constraint_type
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    PetscReal, pointer :: molarity(:)
    PetscInt, pointer :: constraint_type(:)
    PetscInt, pointer :: constraint_specid(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_aux_string(:)
    PetscBool, pointer :: external_dataset(:)
  end type tran_species_constraint_type

  type, public :: transport_auxvar_type
    type(tran_species_type), pointer :: tran_species(:)

    ! molality in solution (aq.)
    PetscReal, pointer :: aq_molal(:)         ! mol/kg water
    ! aq-phase dependent totals
    PetscReal, pointer :: aq_total(:,:)       ! mol solute/L water
    type(matrix_block_auxvar_type), pointer :: aqueous

    ! gases - assuming same species names for both aq. and gas phase
    PetscReal, pointer :: gas_pp(:) ! gas partial pressure in bars in air
    type(matrix_block_auxvar_type), pointer :: gaseous
    
    PetscReal, pointer :: mass_balance(:,:)
    PetscReal, pointer :: mass_balance_delta(:,:)
    
    ! auxiliary array to store miscellaneous data (e.g. cumulative mass, etc.)
    PetscReal, pointer :: auxiliary_data(:)
    
  end type transport_auxvar_type

  type, public :: transport_param_type
    PetscInt :: ncomp
    PetscInt :: nphase
    PetscReal, pointer :: diffusion_coefficient(:,:)         ! (ncomp,nphase)
    PetscReal, pointer :: diffusion_activation_energy(:,:)   ! (ncomp,nphase)
    PetscReal :: newton_inf_rel_update_tol
    PetscReal :: newton_inf_scaled_res_tol
    PetscBool :: check_post_converged
    PetscBool :: calculate_transverse_dispersion
    PetscBool :: temperature_dependent_diffusion
  end type transport_param_type

  type, public :: transport_type
    ! compressed arrays for efficient computation
    PetscInt :: ncomp
    PetscInt :: naqcomp
    PetscInt :: ngascomp

    ! offsets
    PetscInt :: offset_aqueous
    PetscInt :: offset_gaseous

    type(tran_species_type), pointer :: tran_species_list

    !
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscInt :: n_zero_rows
    PetscBool :: inactive_cells_exist
    type(transport_param_type), pointer :: tran_parameter
    type(transport_auxvar_type), pointer :: auxvars(:)
    type(transport_auxvar_type), pointer :: auxvars_bc(:)
    type(transport_auxvar_type), pointer :: auxvars_ss(:)
  end type transport_type

  interface GetTranSpeciesIDFromName
    module procedure GetTranSpeciesIDFromName1
    module procedure GetTranSpeciesIDFromName2
  end interface

  interface TranAuxVarDestroy
    module procedure TranAuxVarSingleDestroy
    module procedure TranAuxVarArrayDestroy
  end interface TranAuxVarDestroy
  
  public :: TranAuxCreate, TranAuxDestroy, &
            TranAuxVarInit, TranAuxVarCopy, TranAuxVarDestroy, &
            TranAuxVarStrip, TranAuxVarCopyInitialGuess

  public :: TransportCreate, &
            GetTranSpeciesCount, &
            GetTranSpeciesNames, &
            GetTranSpeciesIDFromName, &
            TranSpeciesCreate, &
            TranSpeciesDestroy, &
            TranSpeciesConstraintCreate, &
            TranSpeciesConstraintDestroy

contains

! ************************************************************************** !

function TranAuxCreate(ncomp,nphase)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 
  use Option_module

  implicit none
  
  PetscInt :: ncomp
  PetscInt :: nphase
  type(transport_type), pointer :: TranAuxCreate
  type(transport_type), pointer :: aux

  allocate(aux)  

  !species-wised
  aux%ncomp = ncomp
  aux%naqcomp = ncomp
  aux%ngascomp = ncomp

  aux%offset_aqueous = 0
  aux%offset_gaseous = 0

  nullify(aux%tran_species_list)

  ! grid-wised
  aux%num_aux = 0      ! number of tran_auxvars objects for local and ghosted cells
  aux%num_aux_bc = 0   ! number of tran_auxvars objects for boundary connections
  aux%num_aux_ss = 0   ! number of tran_auxvars objects for source/sinks
  nullify(aux%auxvars)      ! tran_auxvars for local and ghosted grid cells
  nullify(aux%auxvars_bc)   ! tran_auxvars for boundary connections
  nullify(aux%auxvars_ss)   ! tran_auxvars for source/sinks
  aux%n_zero_rows = 0       ! number of zeroed rows in Jacobian for inactive cells
  nullify(aux%zero_rows_local)  ! ids of zero rows in local, non-ghosted numbering
  nullify(aux%zero_rows_local_ghosted) ! ids of zero rows in ghosted numbering
  aux%inactive_cells_exist = PETSC_FALSE

  allocate(aux%tran_parameter)
  allocate(aux%tran_parameter%diffusion_coefficient(ncomp,nphase))
  allocate(aux%tran_parameter%diffusion_activation_energy(ncomp,nphase))
  aux%tran_parameter%ncomp = ncomp
  aux%tran_parameter%nphase = nphase
  aux%tran_parameter%diffusion_coefficient = 1.d-9
  aux%tran_parameter%diffusion_activation_energy = 0.d0
  aux%tran_parameter%calculate_transverse_dispersion = PETSC_FALSE
  aux%tran_parameter%temperature_dependent_diffusion = PETSC_FALSE

  TranAuxCreate => aux
  
end function TranAuxCreate

! ************************************************************************** !
subroutine TranAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(transport_type) :: auxvar
  type(option_type) :: option
  
#if 0
! NOT YET (TODO)

  allocate(auxvar%molal(transport%naqcomp))
  auxvar%molal = 0.d0

  allocate(auxvar%total(transport%naqcomp,transport%nphase))
  auxvar%total = 0.d0
  auxvar%aqueous => MatrixBlockAuxVarCreate(option)
  call MatrixBlockAuxVarInit(auxvar%aqueous,transport%naqcomp, &
                             transport%naqcomp,transport%nphase,option)
  

  if (reaction%gas%nactive_gas > 0) then
    allocate(auxvar%gas_pp(reaction%gas%nactive_gas))
    auxvar%gas_pp = 0.d0
  else
    nullify(auxvar%gas_pp)
  endif

  allocate(auxvar%act_coef(transport%naqcomp))
  auxvar%act_coef = 1.d0

! initialize ln activity H2O
  auxvar%ln_act_h2o = 0.d0
  
  if (option%iflag /= 0 .and. option%compute_mass_balance_new) then
    allocate(auxvar%mass_balance(transport%ncomp,transport%nphase))
    auxvar%mass_balance = 0.d0
    allocate(auxvar%mass_balance_delta(transport%ncomp,transport%nphase))
    auxvar%mass_balance_delta = 0.d0
  else
    nullify(auxvar%mass_balance)
    nullify(auxvar%mass_balance_delta)
  endif
  
  if (transport%nauxiliary > 0) then
    allocate(auxvar%auxiliary_data(transportn%nauxiliary))
    auxvar%auxiliary_data = 0.d0
  else
    nullify(auxvar%auxiliary_data)
  endif

#endif

end subroutine TranAuxVarInit

! ************************************************************************** !

subroutine TranAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copys an auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/05/08
  ! 

  use Option_module

  implicit none
  
  type(transport_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option  
  
  auxvar%aq_molal = auxvar2%aq_molal
  auxvar%aq_total = auxvar2%aq_total
  call MatrixBlockAuxVarCopy(auxvar%aqueous,auxvar2%aqueous,option)
  
  if (associated(auxvar%gas_pp)) then
    auxvar%gas_pp = auxvar2%gas_pp
    call MatrixBlockAuxVarCopy(auxvar%aqueous,auxvar2%aqueous,option)
  endif

  if (associated(auxvar%mass_balance)) then
    auxvar%mass_balance = auxvar2%mass_balance
    auxvar%mass_balance_delta = auxvar2%mass_balance_delta
  endif

  if (associated(auxvar%auxiliary_data)) then
    auxvar%auxiliary_data = auxvar2%auxiliary_data
  endif
  
end subroutine TranAuxVarCopy

! ************************************************************************** !

subroutine TranAuxVarCopyInitialGuess(auxvar,auxvar2,option)
  ! 
  ! Copies select entries in rt_auxvar that serve as an initial guess when
  ! equilibrating constraints
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/18/17
  ! 

  use Option_module

  implicit none
  
  type(transport_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option  
  
  auxvar%aq_molal = auxvar2%aq_molal
  if (associated(auxvar%gas_pp)) then
    auxvar%gas_pp = auxvar2%gas_pp
  endif

  
end subroutine TranAuxVarCopyInitialGuess

! ************************************************************************** !

subroutine TranAuxVarSingleDestroy(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(transport_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call TranAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)  

end subroutine TranAuxVarSingleDestroy

! ************************************************************************** !

subroutine TranAuxVarArrayDestroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(transport_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call TranAuxVarStrip(auxvars(iaux))
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine TranAuxVarArrayDestroy

! ************************************************************************** !

subroutine TranAuxVarStrip(auxvar)
  ! 
  ! Deallocates all members of single auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Utility_module, only: DeallocateArray

  implicit none

  type(transport_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%aq_molal)
  call DeallocateArray(auxvar%aq_total)
  call MatrixBlockAuxVarDestroy(auxvar%aqueous)

  call DeallocateArray(auxvar%gas_pp)
  call MatrixBlockAuxVarDestroy(auxvar%gaseous)

  call DeallocateArray(auxvar%mass_balance)
  call DeallocateArray(auxvar%mass_balance_delta)
  
  
  call DeallocateArray(auxvar%auxiliary_data)
  
end subroutine TranAuxVarStrip

! ************************************************************************** !

subroutine TranAuxDestroy(aux)
  ! 
  ! Deallocates a reactive transport auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Utility_module, only: DeallocateArray
  
  implicit none

  type(transport_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  call TranAuxVarDestroy(aux%auxvars)
  call TranAuxVarDestroy(aux%auxvars_bc)
  call TranAuxVarDestroy(aux%auxvars_ss)
  call DeallocateArray(aux%zero_rows_local)
  call DeallocateArray(aux%zero_rows_local_ghosted)

  if (associated(aux%tran_parameter)) then
    call DeallocateArray(aux%tran_parameter%diffusion_coefficient)
    call DeallocateArray(aux%tran_parameter%diffusion_activation_energy)
    deallocate(aux%tran_parameter)
  endif
  nullify(aux%tran_parameter)

  deallocate(aux)
  nullify(aux)

  end subroutine TranAuxDestroy


! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !

function TransportCreate()
  !
  ! Allocate and initialize transport object
  !
  ! Author: Glenn Hammond
  ! Date: 05/02/08
  !
  implicit none

  class(transport_type), pointer :: TransportCreate

  class(transport_type), pointer :: transport

  allocate(transport)

  nullify(transport%tran_species_list)

  transport%ncomp = 0
  transport%naqcomp = 0
  transport%ngascomp = 0
  transport%offset_aqueous = 0
  transport%offset_gaseous = 0

  call TranInit(transport)

  TransportCreate => transport

end function TransportCreate

! ************************************************************************** !

function TranSpeciesCreate()
  !
  ! Allocate and initialize a transport species object
  !
  ! Author: Glenn Hammond
  ! Date: 05/02/08
  !
  implicit none

  type(tran_species_type), pointer :: TranSpeciesCreate

  type(tran_species_type), pointer :: species

  allocate(species)
  species%id = 0
  species%name = ''
  species%molar_weight = 0.d0
  species%is_reactive = PETSC_FALSE
  nullify(species%next)

  TranSpeciesCreate => species

end function TranSpeciesCreate

! ************************************************************************** !

function TranSpeciesConstraintCreate(transport,option)
  !
  ! Creates a transport species constraint
  ! object
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  !
  use Option_module

  implicit none

  class(transport_type) :: transport
  type(option_type) :: option
  type(tran_species_constraint_type), pointer :: TranSpeciesConstraintCreate

  type(tran_species_constraint_type), pointer :: constraint

  allocate(constraint)
  allocate(constraint%names(transport%ncomp))
  constraint%names = ''
  allocate(constraint%constraint_conc(transport%ncomp))
  constraint%constraint_conc = 0.d0
  allocate(constraint%molarity(transport%ncomp))
  constraint%molarity = 0.d0
  allocate(constraint%constraint_specid(transport%ncomp))
  constraint%constraint_specid = 0
  allocate(constraint%constraint_type(transport%ncomp))
  constraint%constraint_type = 0
  allocate(constraint%constraint_aux_string(transport%ncomp))
  constraint%constraint_aux_string = ''
  allocate(constraint%external_dataset(transport%ncomp))
  constraint%external_dataset = PETSC_FALSE

  TranSpeciesConstraintCreate => constraint

end function TranSpeciesConstraintCreate
! ************************************************************************** !

function GetTranSpeciesNames(transport)
  !
  ! Returns the names of transport species in an array
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/08
  !

  implicit none

  character(len=MAXWORDLENGTH), pointer :: GetTranSpeciesNames(:)
  class(transport_type) :: transport

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(tran_species_type), pointer :: species

  count = GetTranSpeciesCount(transport)
  allocate(names(count))

  count = 1
  species => transport%tran_species_list
  do
    if (.not.associated(species)) exit
    names(count) = species%name
    count = count + 1
    species => species%next
  enddo

  GetTranSpeciesNames => names

end function GetTranSpeciesNames

! ************************************************************************** !

function GetTranSpeciesCount(transport)
  !
  ! Returns the number of transport species
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/08
  !

  implicit none

  PetscInt :: GetTranSpeciesCount
  class(transport_type) :: transport

  type(tran_species_type), pointer :: species

  GetTranSpeciesCount = 0
  species => transport%tran_species_list
  do
    if (.not.associated(species)) exit
    GetTranSpeciesCount = GetTranSpeciesCount + 1
    species => species%next
  enddo

end function GetTranSpeciesCount

! ************************************************************************** !

function GetTranSpeciesIDFromName1(name,transport,option)
  !
  ! Returns the id of named transport species
  !
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  !
  use Option_module
  use String_module

  implicit none

  character(len=MAXWORDLENGTH) :: name
  class(transport_type) :: transport
  type(option_type) :: option

  PetscInt :: GetTranSpeciesIDFromName1

  GetTranSpeciesIDFromName1 = GetTranSpeciesIDFromName2(name,transport, &
          PETSC_TRUE, option)

end function GetTranSpeciesIDFromName1

! ************************************************************************** !

function GetTranSpeciesIDFromName2(name,transport,return_error,option)
  !
  ! Returns the id of named primary species
  !
  ! Author: Glenn Hammond
  ! Date: 10/30/12


  use Option_module
  use String_module

  implicit none

  character(len=MAXWORDLENGTH) :: name
  class(transport_type) :: transport
  type(option_type) :: option

  PetscInt :: GetTranSpeciesIDFromName2

  type(tran_species_type), pointer :: species
  PetscInt :: i
  PetscBool :: return_error

  GetTranSpeciesIDFromName2 = UNINITIALIZED_INTEGER

  species => transport%tran_species_list
  i = 0
  do
    if (.not.associated(species)) exit
    i = i + 1
    if (StringCompare(name,species%name,MAXWORDLENGTH)) then
      GetTranSpeciesIDFromName2 = i
      exit
    endif
    species => species%next
  enddo

  if (return_error .and. GetTranSpeciesIDFromName2 <= 0) then
    option%io_buffer = 'Species "' // trim(name) // &
      '" not found among transport species in GetTranSpeciesIDFromName().'
    call PrintErrMsg(option)
  endif

end function GetTranSpeciesIDFromName2

! ************************************************************************** !

subroutine TranSpeciesDestroy(species)
  !
  ! Deallocates an aqueous species
  !
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  !

  implicit none

  type(tran_species_type), pointer :: species

  deallocate(species)
  nullify(species)

end subroutine TranSpeciesDestroy

! ************************************************************************** !

subroutine TranSpeciesListDestroy(species_list)
  !
  ! Deallocates a transport species
  !
  ! Author: Glenn Hammond
  ! Date: 09/03/10
  !

  !TODO(geh): make these destructors recursive
  implicit none

  type(tran_species_type), pointer :: species_list

  type(tran_species_type), pointer :: species, prev_species

  species => species_list
  do
    if (.not.associated(species)) exit
    prev_species => species
    species => species%next
    call TranSpeciesDestroy(prev_species)
  enddo
  nullify(species_list)

end subroutine TranSpeciesListDestroy

! ************************************************************************** !

subroutine TranSpeciesConstraintDestroy(constraint)
  !
  ! Destroys an aqueous species constraint
  ! object
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  !

  use Utility_module, only: DeallocateArray

  implicit none

  type(tran_species_constraint_type), pointer :: constraint

  if (.not.associated(constraint)) return

  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_conc)
  call DeallocateArray(constraint%molarity)
  call DeallocateArray(constraint%constraint_type)
  call DeallocateArray(constraint%constraint_specid)
  call DeallocateArray(constraint%constraint_aux_string)
  call DeallocateArray(constraint%external_dataset)

  deallocate(constraint)
  nullify(constraint)

end subroutine TranSpeciesConstraintDestroy

! ************************************************************************** !

subroutine TransportDestroy(transport,option)
  !
  ! Deallocates a transport object
  !
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  !

  use Utility_module, only: DeallocateArray
  use Option_module

  implicit none

  class(transport_type), pointer :: transport

  type(tran_species_type), pointer :: tran_species, prev_tran_species
  type(option_type) :: option

  if (.not.associated(transport)) return

  ! all transport species
  if (associated(transport%tran_species_list)) &
    call TranSpeciesListDestroy(transport%tran_species_list)
  nullify(transport%tran_species_list)

  deallocate(transport)
  nullify(transport)

end subroutine TransportDestroy

! ************************************************************************** !


! ************************************************************************** !




! ************************************************************************** !
! ************************************************************************** !


end module Transport_Aux_module
