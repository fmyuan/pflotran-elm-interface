module NW_Transport_Aux_module

  ! this module cannot depend on any other modules besides Option_module
  ! and Matrix_Block_Aux_module
  use Matrix_Block_Aux_module

  use PFLOTRAN_Constants_module

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
  
  type, public :: nw_transport_param_type
    PetscInt :: nphase
    PetscInt :: ncomp
    PetscInt :: nsorb
    PetscInt :: nmnrl
    PetscInt :: nauxiliary
    PetscInt :: offset_auxiliary
    PetscReal, pointer :: diffusion_coefficient(:,:)
    PetscReal, pointer :: diffusion_activation_energy(:,:)
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
  end type nw_transport_param_type
  
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
    type(nw_transport_param_type), pointer :: nwt_parameter
    ! nwt_auxvars for local and ghosted grid cells
    type(nw_transport_auxvar_type), pointer :: auxvars(:)
    ! nwt_auxvars for boundary connections
    type(nw_transport_auxvar_type), pointer :: auxvars_bc(:)
    ! nwt_auxvars for source/sinks
    type(nw_transport_auxvar_type), pointer :: auxvars_ss(:)
  end type nw_transport_type

  interface NWTAuxVarDestroy
    module procedure NWTAuxVarSingleDestroy
    module procedure NWTAuxVarArrayDestroy
  end interface NWTAuxVarDestroy
  
  public :: NWTAuxCreate, NWTAuxDestroy, &
            NWTAuxVarInit, NWTAuxVarCopy, NWTAuxVarDestroy, &
            NWTAuxVarStrip, NWTAuxVarCopyInitialGuess
            
contains

! ************************************************************************** !

function NWTAuxCreate(ncomp,nphase)
  ! 
  ! Allocate and initialize nuclear waste transport auxiliary object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module

  implicit none
  
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

  allocate(aux%nwt_parameter)
  aux%nwt_parameter%ncomp = ncomp
  aux%nwt_parameter%nphase = nphase
  aux%nwt_parameter%nsorb = 0
  aux%nwt_parameter%nmnrl = 0
  aux%nwt_parameter%nauxiliary = 0
  allocate(aux%nwt_parameter%diffusion_coefficient(ncomp,nphase))
  allocate(aux%nwt_parameter%diffusion_activation_energy(ncomp,nphase))
  aux%nwt_parameter%diffusion_coefficient = 1.d-9
  aux%nwt_parameter%diffusion_activation_energy = 0.d0
  aux%nwt_parameter%offset_auxiliary = 0
  aux%nwt_parameter%calculate_transverse_dispersion = PETSC_FALSE
  aux%nwt_parameter%temperature_dependent_diffusion = PETSC_FALSE
#ifdef OS_STATISTICS
  aux%nwt_parameter%newton_call_count = 0
  aux%nwt_parameter%sum_newton_call_count = 0.d0
  aux%nwt_parameter%newton_iterations = 0
  aux%nwt_parameter%sum_newton_iterations = 0.d0
  aux%nwt_parameter%max_newton_iterations = 0
  aux%nwt_parameter%overall_max_newton_iterations = 0
#endif  

  NWTAuxCreate => aux
  
end function NWTAuxCreate

! ************************************************************************** !

subroutine NWTAuxVarInit(auxvar,nwt_parameter,option)
  ! 
  ! Initializes the nuclear waste auxiliary object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/11/2019
  !
  
  use Option_module

  implicit none
  
  type(nw_transport_auxvar_type) :: auxvar
  type(nw_transport_param_type) :: nwt_parameter
  type(option_type) :: option
  
  PetscInt :: ncomp, nmnrl
  
  ncomp = nwt_parameter%ncomp
  nmnrl = nwt_parameter%nmnrl
  
  allocate(auxvar%molality(ncomp))
  auxvar%molality = 0.d0
  
  if (nwt_parameter%nsorb > 0) then
    allocate(auxvar%total_sorb_eq(ncomp))
    auxvar%total_sorb_eq = 0.d0
    allocate(auxvar%dtotal_sorb_eq(ncomp,ncomp))
    auxvar%dtotal_sorb_eq = 0.d0
  else
    nullify(auxvar%total_sorb_eq)
    nullify(auxvar%dtotal_sorb_eq)
  endif
  
  if (nwt_parameter%nmnrl > 0) then
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
  
  if (nwt_parameter%nauxiliary > 0) then
    allocate(auxvar%auxiliary_data(nwt_parameter%nauxiliary))
    auxvar%auxiliary_data = 0.d0
  else
    nullify(auxvar%auxiliary_data)
  endif
  
end subroutine NWTAuxVarInit

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

  if (associated(aux%nwt_parameter)) then
    call DeallocateArray(aux%nwt_parameter%diffusion_coefficient)
    call DeallocateArray(aux%nwt_parameter%diffusion_activation_energy)
    deallocate(aux%nwt_parameter)
  endif
  nullify(aux%nwt_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine NWTAuxDestroy

! ************************************************************************** !

end module NW_Transport_Aux_module