module Reactive_Transport_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  ! this module cannot depend on any other modules besides Option_module
  ! and Matrix_Block_Aux_module
  use Matrix_Block_Aux_module
  use Matrix_Zeroing_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscReal, public :: rt_itol_scaled_res = UNINITIALIZED_DOUBLE
  PetscReal, public :: rt_itol_rel_update = UNINITIALIZED_DOUBLE
  PetscReal, public :: rt_min_saturation = 1.d-40

  type, public :: reactive_transport_auxvar_type
    ! molality
    PetscReal, pointer :: pri_molal(:)     ! mol/kg water

    ! phase dependent totals
    PetscReal, pointer :: total(:,:)       ! mol solute/L water
    type(matrix_block_auxvar_type), pointer :: aqueous

    ! sorbed totals
    PetscReal, pointer :: total_sorb_eq(:)    ! mol/m^3 bulk
    PetscReal, pointer :: dtotal_sorb_eq(:,:) ! kg water/m^3 bulk

    ! aqueous complexes
    PetscReal, pointer :: sec_molal(:)

    ! gases
    PetscReal, pointer :: gas_pp(:) ! gas partial pressure in bars

    ! sorption reactions
    PetscReal, pointer :: srfcplxrxn_free_site_conc(:)
    PetscReal, pointer :: kinsrfcplx_conc(:,:) ! S_{i\alpha}^k
    PetscReal, pointer :: kinsrfcplx_conc_kp1(:,:) ! S_{i\alpha}^k+1
    PetscReal, pointer :: kinsrfcplx_free_site_conc(:)  ! S_\alpha
    PetscReal, pointer :: eqsrfcplx_conc(:)
    PetscReal, pointer :: eqionx_ref_cation_sorbed_conc(:)
    PetscReal, pointer :: eqionx_conc(:,:)

    ! mineral reactions
    PetscReal, pointer :: mnrl_volfrac0(:)
    PetscReal, pointer :: mnrl_volfrac(:)
    PetscReal, pointer :: mnrl_area0(:)
    PetscReal, pointer :: mnrl_area(:)
    PetscReal, pointer :: mnrl_rate(:)

    ! activity coefficients
!   PetscReal :: act_h2o
    PetscReal, pointer :: pri_act_coef(:)
    PetscReal, pointer :: sec_act_coef(:)

    PetscReal :: ln_act_h2o

    PetscReal, pointer :: mass_balance(:,:)
    PetscReal, pointer :: mass_balance_delta(:,:)

    PetscReal, pointer :: kinmr_total_sorb(:,:,:)

    ! immobile species such as biomass
    PetscReal, pointer :: immobile(:) ! mol/m^3 bulk

    ! auxiliary array to store miscellaneous data (e.g. reaction rates,
    ! cumulative mass, etc.
    PetscReal, pointer :: auxiliary_data(:)

  end type reactive_transport_auxvar_type

  type, public :: reactive_transport_param_type
    PetscInt :: nphase
    PetscInt :: ncomp
    PetscInt :: naqcomp
    PetscInt :: nimcomp
    PetscInt :: ngas
    PetscInt :: offset_aqueous
    PetscInt :: offset_immobile
    PetscInt :: offset_auxiliary
    PetscBool :: species_dependent_diffusion
    PetscBool :: millington_quirk_tortuosity
    PetscReal, pointer :: diffusion_coefficient(:,:)
    PetscReal, pointer :: diffusion_activation_energy(:,:)
    PetscReal, pointer :: pri_spec_diff_coef(:)
    PetscReal, pointer :: sec_spec_diff_coef(:)
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
  end type reactive_transport_param_type

  type, public :: reactive_transport_type
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    PetscBool :: inactive_cells_exist
    type(reactive_transport_param_type), pointer :: rt_parameter
    type(reactive_transport_auxvar_type), pointer :: auxvars(:)
    type(reactive_transport_auxvar_type), pointer :: auxvars_bc(:)
    type(reactive_transport_auxvar_type), pointer :: auxvars_ss(:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type reactive_transport_type

  interface RTAuxVarDestroy
    module procedure RTAuxVarSingleDestroy
    module procedure RTAuxVarArrayDestroy
  end interface RTAuxVarDestroy

  public :: RTAuxCreate, RTAuxDestroy, &
            RTAuxVarInit, RTAuxVarCopy, RTAuxVarDestroy, &
            RTAuxVarStrip, RTAuxVarCopyInitialGuess

contains

! ************************************************************************** !

function RTAuxCreate(naqcomp,nphase)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  !
  use Option_module

  implicit none

  PetscInt :: naqcomp
  PetscInt :: nphase
  type(reactive_transport_type), pointer :: RTAuxCreate

  type(reactive_transport_type), pointer :: aux

  allocate(aux)
  aux%num_aux = 0      ! number of rt_auxvars objects for local and ghosted cells
  aux%num_aux_bc = 0   ! number of rt_auxvars objects for boundary connections
  aux%num_aux_ss = 0   ! number of rt_auxvars objects for source/sinks
  nullify(aux%auxvars)      ! rt_auxvars for local and ghosted grid cells
  nullify(aux%auxvars_bc)   ! rt_auxvars for boundary connections
  nullify(aux%auxvars_ss)   ! rt_auxvars for source/sinks
  nullify(aux%matrix_zeroing)
  aux%inactive_cells_exist = PETSC_FALSE

  allocate(aux%rt_parameter)
  aux%rt_parameter%naqcomp = naqcomp
  aux%rt_parameter%nphase = nphase
  allocate(aux%rt_parameter%diffusion_coefficient(naqcomp,nphase))
  allocate(aux%rt_parameter%diffusion_activation_energy(naqcomp,nphase))
  aux%rt_parameter%diffusion_coefficient = 1.d-9
  aux%rt_parameter%diffusion_activation_energy = 0.d0
  nullify(aux%rt_parameter%pri_spec_diff_coef)
  nullify(aux%rt_parameter%sec_spec_diff_coef)
  aux%rt_parameter%ncomp = 0
  aux%rt_parameter%nimcomp = 0
  aux%rt_parameter%ngas = 0
  aux%rt_parameter%offset_aqueous = 0
  aux%rt_parameter%offset_immobile = 0
  aux%rt_parameter%offset_auxiliary = 0
  aux%rt_parameter%species_dependent_diffusion = PETSC_FALSE
  aux%rt_parameter%millington_quirk_tortuosity = PETSC_FALSE
  aux%rt_parameter%calculate_transverse_dispersion = PETSC_FALSE
  aux%rt_parameter%temperature_dependent_diffusion = PETSC_FALSE
#ifdef OS_STATISTICS
  aux%rt_parameter%newton_call_count = 0
  aux%rt_parameter%sum_newton_call_count = 0.d0
  aux%rt_parameter%newton_iterations = 0
  aux%rt_parameter%sum_newton_iterations = 0.d0
  aux%rt_parameter%max_newton_iterations = 0
  aux%rt_parameter%overall_max_newton_iterations = 0
#endif

  RTAuxCreate => aux

end function RTAuxCreate

! ************************************************************************** !

subroutine RTAuxVarInit(auxvar,reaction,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  !

  use Option_module
  use Reaction_Aux_module
  use Reaction_Surface_Complexation_Aux_module

  implicit none

  type(reactive_transport_auxvar_type) :: auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  type(surface_complexation_type), pointer :: surface_complexation

  surface_complexation => reaction%surface_complexation

  allocate(auxvar%pri_molal(reaction%naqcomp))
  auxvar%pri_molal = 0.d0

  allocate(auxvar%total(reaction%naqcomp,reaction%nphase))
  auxvar%total = 0.d0
  auxvar%aqueous => MatrixBlockAuxVarCreate(option)
  call MatrixBlockAuxVarInit(auxvar%aqueous,reaction%naqcomp, &
                             reaction%naqcomp,reaction%nphase,option)

  if (reaction%neqcplx > 0) then
    allocate(auxvar%sec_molal(reaction%neqcplx))
    auxvar%sec_molal = 0.d0
  else
    nullify(auxvar%sec_molal)
  endif

  if (reaction%gas%nactive_gas > 0) then
    allocate(auxvar%gas_pp(reaction%gas%nactive_gas))
    auxvar%gas_pp = 0.d0
  else
    nullify(auxvar%gas_pp)
  endif

  if (reaction%neqsorb > 0) then
    allocate(auxvar%total_sorb_eq(reaction%naqcomp))
    auxvar%total_sorb_eq = 0.d0
    allocate(auxvar%dtotal_sorb_eq(reaction%naqcomp,reaction%naqcomp))
    auxvar%dtotal_sorb_eq = 0.d0
  else
    nullify(auxvar%total_sorb_eq)
    nullify(auxvar%dtotal_sorb_eq)
  endif

  ! surface complexation
  nullify(auxvar%eqsrfcplx_conc)
  nullify(auxvar%srfcplxrxn_free_site_conc)
  nullify(auxvar%kinsrfcplx_conc)
  nullify(auxvar%kinsrfcplx_conc_kp1)
  nullify(auxvar%kinsrfcplx_free_site_conc)
  nullify(auxvar%kinmr_total_sorb)
  if (surface_complexation%nsrfcplxrxn > 0) then
    allocate(auxvar%srfcplxrxn_free_site_conc(surface_complexation%nsrfcplxrxn))
    auxvar%srfcplxrxn_free_site_conc = 1.d-9 ! initialize to guess
    if (surface_complexation%neqsrfcplxrxn > 0) then
      allocate(auxvar%eqsrfcplx_conc(surface_complexation%nsrfcplx))
      auxvar%eqsrfcplx_conc = 0.d0
    endif
    if (surface_complexation%nkinsrfcplxrxn > 0) then
      !geh: currently hardwired to only 1 reaction
      allocate(auxvar%kinsrfcplx_conc(surface_complexation%nkinsrfcplx,1))
      auxvar%kinsrfcplx_conc = 0.d0

      allocate(auxvar%kinsrfcplx_conc_kp1(surface_complexation%nkinsrfcplx,1))
      auxvar%kinsrfcplx_conc_kp1 = 0.d0
    endif
    if (surface_complexation%nkinmrsrfcplxrxn > 0) then
      ! the zeroth entry here stores the equilibrium concentration used in the
      ! update
      ! the zeroth entry of kinmr_nrate holds the maximum number of rates
      ! prescribed in a multirate reaction...required for appropriate sizing
      allocate(auxvar%kinmr_total_sorb(reaction%naqcomp, &
                                        0:surface_complexation%kinmr_nrate(0), &
                                        surface_complexation%nkinmrsrfcplxrxn))
      auxvar%kinmr_total_sorb = 0.d0
    endif
  endif

  if (reaction%neqionxrxn > 0) then
    allocate(auxvar%eqionx_ref_cation_sorbed_conc(reaction%neqionxrxn))
    auxvar%eqionx_ref_cation_sorbed_conc = 1.d-9 ! initialize to guess

    allocate(auxvar%eqionx_conc(reaction%neqionxcation,reaction%neqionxrxn))
    auxvar%eqionx_conc = 1.d-9

!   allocate(auxvar%eqionx_cec(reaction%neqionxcation))
!   auxvar%eqionx_cec = 0.d0
  else
    nullify(auxvar%eqionx_ref_cation_sorbed_conc)
    nullify(auxvar%eqionx_conc)
!   nullify(auxvar%eqionx_cec)
  endif

  if (associated(reaction%mineral)) then
    if (reaction%mineral%nkinmnrl > 0) then
      allocate(auxvar%mnrl_volfrac0(reaction%mineral%nkinmnrl))
      auxvar%mnrl_volfrac0 = 0.d0
      allocate(auxvar%mnrl_volfrac(reaction%mineral%nkinmnrl))
      auxvar%mnrl_volfrac = 0.d0
      allocate(auxvar%mnrl_area0(reaction%mineral%nkinmnrl))
      auxvar%mnrl_area0 = 0.d0
      allocate(auxvar%mnrl_area(reaction%mineral%nkinmnrl))
      auxvar%mnrl_area = 0.d0
      allocate(auxvar%mnrl_rate(reaction%mineral%nkinmnrl))
      auxvar%mnrl_rate = 0.d0
    else
      nullify(auxvar%mnrl_volfrac0)
      nullify(auxvar%mnrl_volfrac)
      nullify(auxvar%mnrl_area0)
      nullify(auxvar%mnrl_area)
      nullify(auxvar%mnrl_rate)
    endif
  endif

  allocate(auxvar%pri_act_coef(reaction%naqcomp))
  auxvar%pri_act_coef = 1.d0
  if (reaction%neqcplx > 0) then
    allocate(auxvar%sec_act_coef(reaction%neqcplx))
    auxvar%sec_act_coef = 1.d0
  else
    nullify(auxvar%sec_act_coef)
  endif

! initialize ln activity H2O
  auxvar%ln_act_h2o = 0.d0

  if (option%iflag /= 0 .and. option%compute_mass_balance_new) then
    allocate(auxvar%mass_balance(reaction%ncomp,reaction%nphase))
    auxvar%mass_balance = 0.d0
    allocate(auxvar%mass_balance_delta(reaction%ncomp,reaction%nphase))
    auxvar%mass_balance_delta = 0.d0
  else
    nullify(auxvar%mass_balance)
    nullify(auxvar%mass_balance_delta)
  endif

  if (reaction%immobile%nimmobile > 0) then
    allocate(auxvar%immobile(reaction%immobile%nimmobile))
    auxvar%immobile = 0.d0
  else
    nullify(auxvar%immobile)
  endif

  if (reaction%nauxiliary > 0) then
    allocate(auxvar%auxiliary_data(reaction%nauxiliary))
    auxvar%auxiliary_data = 0.d0
  else
    nullify(auxvar%auxiliary_data)
  endif

end subroutine RTAuxVarInit

! ************************************************************************** !

subroutine RTAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copys an auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/05/08
  !

  use Option_module

  implicit none

  type(reactive_transport_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pri_molal = auxvar%pri_molal

  auxvar2%total = auxvar%total

  call MatrixBlockAuxVarCopy(auxvar%aqueous,auxvar2%aqueous,option)

  if (associated(auxvar%sec_molal)) &
    auxvar2%sec_molal = auxvar%sec_molal
  if (associated(auxvar%total_sorb_eq)) then
    auxvar2%total_sorb_eq = auxvar%total_sorb_eq
  endif
  if (associated(auxvar%dtotal_sorb_eq)) then
    auxvar2%dtotal_sorb_eq = auxvar%dtotal_sorb_eq
  endif

  if (associated(auxvar%gas_pp)) &
    auxvar2%gas_pp = auxvar%gas_pp

  if (associated(auxvar%srfcplxrxn_free_site_conc)) then
    auxvar2%srfcplxrxn_free_site_conc = auxvar%srfcplxrxn_free_site_conc
  endif

  if (associated(auxvar%eqsrfcplx_conc)) then
    auxvar2%eqsrfcplx_conc = auxvar%eqsrfcplx_conc
  endif

  if (associated(auxvar%kinsrfcplx_conc)) then
    auxvar2%kinsrfcplx_conc = auxvar%kinsrfcplx_conc
    auxvar2%kinsrfcplx_conc_kp1 = auxvar%kinsrfcplx_conc_kp1
    auxvar2%kinsrfcplx_free_site_conc = auxvar%kinsrfcplx_free_site_conc
  endif

  if (associated(auxvar%eqionx_ref_cation_sorbed_conc)) then
    auxvar2%eqionx_ref_cation_sorbed_conc = &
      auxvar%eqionx_ref_cation_sorbed_conc
    auxvar2%eqionx_conc = auxvar%eqionx_conc
  endif

  if (associated(auxvar%mnrl_volfrac)) then
    auxvar2%mnrl_volfrac0 = auxvar%mnrl_volfrac0
    auxvar2%mnrl_volfrac = auxvar%mnrl_volfrac
    auxvar2%mnrl_area0 = auxvar%mnrl_area0
    auxvar2%mnrl_area = auxvar%mnrl_area
    auxvar2%mnrl_rate = auxvar%mnrl_rate
  endif

  auxvar2%ln_act_h2o = auxvar%ln_act_h2o

  auxvar2%pri_act_coef = auxvar%pri_act_coef
  if (associated(auxvar%sec_act_coef)) &
    auxvar2%sec_act_coef = auxvar%sec_act_coef

  if (associated(auxvar%mass_balance)) then
    auxvar2%mass_balance = auxvar%mass_balance
    auxvar2%mass_balance_delta = auxvar%mass_balance_delta
  endif

  if (associated(auxvar%kinmr_total_sorb)) then
    auxvar2%kinmr_total_sorb = auxvar%kinmr_total_sorb
  endif

  if (associated(auxvar%immobile)) then
    auxvar2%immobile = auxvar%immobile
  endif

  if (associated(auxvar%auxiliary_data)) then
    auxvar2%auxiliary_data = auxvar%auxiliary_data
  endif

end subroutine RTAuxVarCopy

! ************************************************************************** !

subroutine RTAuxVarCopyInitialGuess(auxvar,auxvar2,option)
  !
  ! Copies select entries in rt_auxvar that serve as an initial guess when
  ! equilibrating constraints
  !
  ! Author: Glenn Hammond
  ! Date: 07/18/17
  !

  use Option_module

  implicit none

  type(reactive_transport_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pri_molal = auxvar%pri_molal

#if 0
  if (associated(auxvar%srfcplxrxn_free_site_conc)) then
    auxvar2%srfcplxrxn_free_site_conc = auxvar%srfcplxrxn_free_site_conc
  endif

  auxvar2%pri_act_coef = auxvar%pri_act_coef
  if (associated(auxvar%sec_act_coef)) &
    auxvar2%sec_act_coef = auxvar%sec_act_coef

  if (associated(auxvar%immobile)) then
    auxvar2%immobile = auxvar%immobile
  endif
#endif

end subroutine RTAuxVarCopyInitialGuess

! ************************************************************************** !

subroutine RTAuxVarSingleDestroy(auxvar)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  !

  implicit none

  type(reactive_transport_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call RTAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine RTAuxVarSingleDestroy

! ************************************************************************** !

subroutine RTAuxVarArrayDestroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  !

  implicit none

  type(reactive_transport_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call RTAuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine RTAuxVarArrayDestroy

! ************************************************************************** !

subroutine RTAuxVarStrip(auxvar)
  !
  ! Deallocates all members of single auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  !

  use Utility_module, only: DeallocateArray

  implicit none

  type(reactive_transport_auxvar_type) :: auxvar

  call DeallocateArray(auxvar%pri_molal)
  call DeallocateArray(auxvar%total)

  call MatrixBlockAuxVarDestroy(auxvar%aqueous)

  call DeallocateArray(auxvar%sec_molal)
  call DeallocateArray(auxvar%gas_pp)
  call DeallocateArray(auxvar%total_sorb_eq)
  call DeallocateArray(auxvar%dtotal_sorb_eq)

  call DeallocateArray(auxvar%eqsrfcplx_conc)
  call DeallocateArray(auxvar%srfcplxrxn_free_site_conc)
  call DeallocateArray(auxvar%kinsrfcplx_conc)
  call DeallocateArray(auxvar%kinsrfcplx_conc_kp1)
  call DeallocateArray(auxvar%kinsrfcplx_free_site_conc)

  call DeallocateArray(auxvar%eqionx_ref_cation_sorbed_conc)
  call DeallocateArray(auxvar%eqionx_conc)

  call DeallocateArray(auxvar%mnrl_volfrac0)
  call DeallocateArray(auxvar%mnrl_volfrac)
  call DeallocateArray(auxvar%mnrl_area0)
  call DeallocateArray(auxvar%mnrl_area)
  call DeallocateArray(auxvar%mnrl_rate)

  call DeallocateArray(auxvar%pri_act_coef)
  call DeallocateArray(auxvar%sec_act_coef)

  call DeallocateArray(auxvar%mass_balance)
  call DeallocateArray(auxvar%mass_balance_delta)

  call DeallocateArray(auxvar%kinmr_total_sorb)

  call DeallocateArray(auxvar%immobile)
  call DeallocateArray(auxvar%auxiliary_data)

end subroutine RTAuxVarStrip

! ************************************************************************** !

subroutine RTAuxDestroy(aux)
  !
  ! Deallocates a reactive transport auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  !

  use Utility_module, only: DeallocateArray

  implicit none

  type(reactive_transport_type), pointer :: aux

  if (.not.associated(aux)) return

  call RTAuxVarDestroy(aux%auxvars)
  call RTAuxVarDestroy(aux%auxvars_bc)
  call RTAuxVarDestroy(aux%auxvars_ss)
  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%rt_parameter)) then
    call DeallocateArray(aux%rt_parameter%diffusion_coefficient)
    call DeallocateArray(aux%rt_parameter%diffusion_activation_energy)
    call DeallocateArray(aux%rt_parameter%pri_spec_diff_coef)
    call DeallocateArray(aux%rt_parameter%pri_spec_diff_coef)
    deallocate(aux%rt_parameter)
  endif
  nullify(aux%rt_parameter)

  deallocate(aux)
  nullify(aux)

  end subroutine RTAuxDestroy

end module Reactive_Transport_Aux_module
