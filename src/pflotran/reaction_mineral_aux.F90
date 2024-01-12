module Reaction_Mineral_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Database_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  ! mineral types
  PetscInt, parameter, public :: MINERAL_REFERENCE = 1
  PetscInt, parameter, public :: MINERAL_KINETIC = 2
  PetscInt, parameter, public :: MINERAL_EQUILIBRIUM = 3

  type, public :: mineral_rxn_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: molar_volume
    PetscReal :: molar_weight
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn
    type(transition_state_rxn_type), pointer :: tstrxn
    type(mineral_rxn_type), pointer :: next
  end type mineral_rxn_type

  type, public :: transition_state_rxn_type
    PetscReal :: min_scale_factor
    PetscReal :: affinity_factor_sigma
    PetscReal :: affinity_factor_beta
    PetscReal :: affinity_threshold
    PetscReal :: rate_limiter
    PetscReal :: surf_area_vol_frac_pwr
    PetscReal :: surf_area_porosity_pwr
    PetscInt :: irreversible
    PetscReal :: rate
    PetscReal :: activation_energy
    character(len=MAXWORDLENGTH) :: armor_min_name
    PetscReal :: armor_pwr
    PetscReal :: armor_crit_vol_frac
    PetscReal :: surf_area_epsilon
    PetscReal :: vol_frac_epsilon
    type(transition_state_prefactor_type), pointer :: prefactor
    type(transition_state_rxn_type), pointer :: next
  end type transition_state_rxn_type

  type, public :: transition_state_prefactor_type
    type(ts_prefactor_species_type), pointer :: species
    ! these supercede the those above in transition_state_rxn_type
    PetscReal :: rate
    PetscReal :: activation_energy
    type(transition_state_prefactor_type), pointer :: next
  end type transition_state_prefactor_type

  type, public :: ts_prefactor_species_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: id
    PetscReal :: alpha
    PetscReal :: beta
    PetscReal :: attenuation_coef
    type(ts_prefactor_species_type), pointer :: next
  end type ts_prefactor_species_type

  type, public :: mineral_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_vol_frac(:)
    PetscReal, pointer :: constraint_area(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_vol_frac_string(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_area_string(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_area_units(:)
    PetscReal, pointer :: constraint_area_conv_factor(:)
    PetscBool, pointer :: area_per_unit_mass(:)
    PetscBool, pointer :: external_vol_frac_dataset(:)
    PetscBool, pointer :: external_area_dataset(:)
  end type mineral_constraint_type

  type, public :: mineral_type

    PetscInt :: nmnrl
    PetscBool :: print_all
    PetscBool :: print_saturation_index
    PetscBool :: print_rate
    PetscBool :: print_volume_fraction
    PetscBool :: print_surface_area
    character(len=MAXWORDLENGTH), pointer :: mineral_names(:)

    type(mineral_rxn_type), pointer :: mineral_list

    ! for saturation states
    PetscInt, pointer :: mnrlspecid(:,:)
    PetscReal, pointer :: mnrlstoich(:,:)
    PetscInt, pointer :: mnrlh2oid(:)
    PetscReal, pointer :: mnrlh2ostoich(:)
    PetscReal, pointer :: mnrl_logK(:)
    PetscReal, pointer :: mnrl_logKcoef(:,:)
    PetscBool, pointer :: mnrl_print(:)

    ! for kinetic reactions
    PetscInt :: nkinmnrl
    character(len=MAXWORDLENGTH), pointer :: kinmnrl_names(:)
    character(len=MAXWORDLENGTH), pointer :: kinmnrl_armor_min_names(:)
    PetscBool, pointer :: kinmnrl_print(:)
    PetscInt, pointer :: kinmnrlspecid(:,:)
    PetscReal, pointer :: kinmnrlstoich(:,:)
    PetscInt, pointer :: kinmnrlh2oid(:)
    PetscReal, pointer :: kinmnrlh2ostoich(:)
    PetscReal, pointer :: kinmnrl_logK(:)
    PetscReal, pointer :: kinmnrl_logKcoef(:,:)
    PetscReal, pointer :: kinmnrl_rate_constant(:)
    PetscReal, pointer :: kinmnrl_activation_energy(:)
    PetscReal, pointer :: kinmnrl_molar_vol(:)
    PetscReal, pointer :: kinmnrl_molar_wt(:)
    PetscInt, pointer :: kinmnrl_num_prefactors(:)
    PetscInt, pointer :: kinmnrl_prefactor_id(:,:,:)
    PetscReal, pointer :: kinmnrl_pref_alpha(:,:,:)
    PetscReal, pointer :: kinmnrl_pref_beta(:,:,:)
    PetscReal, pointer :: kinmnrl_pref_atten_coef(:,:,:)
    PetscReal, pointer :: kinmnrl_pref_rate(:,:)
    PetscReal, pointer :: kinmnrl_pref_activation_energy(:,:)
    PetscReal, pointer :: kinmnrl_min_scale_factor(:)
    PetscReal, pointer :: kinmnrl_Temkin_const(:)
    PetscReal, pointer :: kinmnrl_affinity_power(:)
    PetscReal, pointer :: kinmnrl_affinity_threshold(:)
    PetscReal, pointer :: kinmnrl_rate_limiter(:)
    PetscReal, pointer :: kinmnrl_surf_area_vol_frac_pwr(:)
    PetscReal, pointer :: kinmnrl_surf_area_porosity_pwr(:)
    PetscReal, pointer :: kinmnrl_armor_crit_vol_frac(:)
    PetscReal, pointer :: kinmnrl_armor_pwr(:)
    PetscReal, pointer :: kinmnrl_surf_area_epsilon(:)
    PetscReal, pointer :: kinmnrl_vol_frac_epsilon(:)
    PetscInt, pointer :: kinmnrl_irreversible(:)

  end type mineral_type

  interface ReactionMnrlGetMnrlIDFromName
    module procedure ReactionMnrlGetMnrlIDFromName1
    module procedure ReactionMnrlGetMnrlIDFromName2
  end interface

  public :: ReactionMnrlCreateMineralObject, &
            ReactionMnrlGetMineralCount, &
            ReactionMnrlGetMineralNames, &
            ReactionMnrlGetMnrlIDFromName, &
            ReactionMnrlGetKinMnrlIDFromName, &
            ReactionMnrlGetMnrlPtrFromName, &
            ReactionMnrlCreateTSTRxn, &
            ReactionMnrlCreateTSTPrefactor, &
            ReactionMnrlCreateTSTPrefSpec, &
            ReactionMnrlDestroyTSTRxn, &
            ReactionMnrlCreateMineralRxn, &
            ReactionMnrlDestroyMineralRxn, &
            ReactionMnrlCreateMnrlConstraint, &
            ReactionMnrlDestMnrlConstraint, &
            ReactionMnrlDestoyMineral

contains

! ************************************************************************** !

function ReactionMnrlCreateMineralObject()
  !
  ! Allocate and initialize mineral reaction object
  !
  ! Author: Glenn Hammond
  ! Date: 08/16/12
  !
  implicit none

  type(mineral_type), pointer :: ReactionMnrlCreateMineralObject

  type(mineral_type), pointer :: mineral

  allocate(mineral)

  nullify(mineral%mineral_list)
  mineral%print_all = PETSC_FALSE
  mineral%print_volume_fraction = PETSC_TRUE
  mineral%print_rate = PETSC_TRUE
  mineral%print_saturation_index = PETSC_FALSE
  mineral%print_surface_area = PETSC_FALSE

  ! for saturation states
  mineral%nmnrl = 0
  nullify(mineral%mineral_names)
  nullify(mineral%mnrl_print)
  nullify(mineral%mnrlspecid)
  nullify(mineral%mnrlh2oid)
  nullify(mineral%mnrlstoich)
  nullify(mineral%mnrlh2ostoich)
  nullify(mineral%mnrl_logK)
  nullify(mineral%mnrl_logKcoef)

  ! for kinetic mineral reactions
  mineral%nkinmnrl = 0
  nullify(mineral%kinmnrl_names)
  nullify(mineral%kinmnrl_print)
  nullify(mineral%kinmnrlspecid)
  nullify(mineral%kinmnrlstoich)
  nullify(mineral%kinmnrlh2oid)
  nullify(mineral%kinmnrlh2ostoich)
  nullify(mineral%kinmnrl_logK)
  nullify(mineral%kinmnrl_logKcoef)
  nullify(mineral%kinmnrl_rate_constant)
  nullify(mineral%kinmnrl_activation_energy)
  nullify(mineral%kinmnrl_molar_vol)
  nullify(mineral%kinmnrl_molar_wt)

  nullify(mineral%kinmnrl_num_prefactors)
  nullify(mineral%kinmnrl_prefactor_id)
  nullify(mineral%kinmnrl_pref_alpha)
  nullify(mineral%kinmnrl_pref_beta)
  nullify(mineral%kinmnrl_pref_atten_coef)
  nullify(mineral%kinmnrl_pref_rate)
  nullify(mineral%kinmnrl_pref_activation_energy)

  nullify(mineral%kinmnrl_min_scale_factor)
  nullify(mineral%kinmnrl_Temkin_const)
  nullify(mineral%kinmnrl_affinity_power)
  nullify(mineral%kinmnrl_affinity_threshold)
  nullify(mineral%kinmnrl_irreversible)
  nullify(mineral%kinmnrl_rate_limiter)
  nullify(mineral%kinmnrl_surf_area_vol_frac_pwr)
  nullify(mineral%kinmnrl_surf_area_porosity_pwr)

  nullify(mineral%kinmnrl_armor_min_names)
  nullify(mineral%kinmnrl_armor_crit_vol_frac)
  nullify(mineral%kinmnrl_armor_pwr)

  nullify(mineral%kinmnrl_surf_area_epsilon)
  nullify(mineral%kinmnrl_vol_frac_epsilon)

  ReactionMnrlCreateMineralObject => mineral

end function ReactionMnrlCreateMineralObject

! ************************************************************************** !

function ReactionMnrlCreateMineralRxn()
  !
  ! Allocate and initialize a mineral object
  !
  ! Author: Glenn Hammond
  ! Date: 05/02/08
  !
  implicit none

  type(mineral_rxn_type), pointer :: ReactionMnrlCreateMineralRxn

  type(mineral_rxn_type), pointer :: mineral

  allocate(mineral)
  mineral%id = 0
  mineral%itype = 0
  mineral%name = ''
  mineral%molar_volume = 0.d0
  mineral%molar_weight = 0.d0
  mineral%print_me = PETSC_FALSE
  nullify(mineral%tstrxn)
  nullify(mineral%next)

  ReactionMnrlCreateMineralRxn => mineral

end function ReactionMnrlCreateMineralRxn

! ************************************************************************** !

function ReactionMnrlCreateTSTRxn()
  !
  ! Allocate and initialize a transition state
  ! theory reaction
  !
  ! Author: Glenn Hammond
  ! Date: 09/01/08
  !

  implicit none

  type(transition_state_rxn_type), pointer :: ReactionMnrlCreateTSTRxn

  type(transition_state_rxn_type), pointer :: tstrxn

  allocate(tstrxn)
  tstrxn%min_scale_factor = UNINITIALIZED_DOUBLE
  tstrxn%affinity_factor_sigma = UNINITIALIZED_DOUBLE
  tstrxn%affinity_factor_beta = UNINITIALIZED_DOUBLE
  tstrxn%affinity_threshold = 0.d0
  tstrxn%surf_area_vol_frac_pwr = 0.d0
  tstrxn%surf_area_porosity_pwr = 0.d0
  tstrxn%rate_limiter = 0.d0
  tstrxn%irreversible = 0
  tstrxn%activation_energy = 0.d0
  tstrxn%armor_min_name = ''
  tstrxn%armor_pwr = 0.d0
  tstrxn%armor_crit_vol_frac = 0.d0
  tstrxn%surf_area_epsilon = 0.d0
  tstrxn%vol_frac_epsilon = 0.d0
  tstrxn%rate = 0.d0
  nullify(tstrxn%prefactor)
  nullify(tstrxn%next)

  ReactionMnrlCreateTSTRxn => tstrxn

end function ReactionMnrlCreateTSTRxn

! ************************************************************************** !

function ReactionMnrlCreateTSTPrefactor()
  !
  ! Allocate and initialize a transition state
  ! theory prefactor
  !
  ! Author: Glenn Hammond
  ! Date: 07/29/11
  !

  implicit none

  type(transition_state_prefactor_type), pointer :: &
    ReactionMnrlCreateTSTPrefactor

  type(transition_state_prefactor_type), pointer :: prefactor

  allocate(prefactor)
  prefactor%rate = 0.d0
  prefactor%activation_energy = 0.d0
  nullify(prefactor%species)
  nullify(prefactor%next)

  ReactionMnrlCreateTSTPrefactor => prefactor

end function ReactionMnrlCreateTSTPrefactor

! ************************************************************************** !

function ReactionMnrlCreateTSTPrefSpec()
  !
  ! Allocate and initialize a transition state
  ! theory prefactor species
  !
  ! Author: Glenn Hammond
  ! Date: 08/01/11
  !

  implicit none

  type(ts_prefactor_species_type), pointer :: ReactionMnrlCreateTSTPrefSpec

  type(ts_prefactor_species_type), pointer :: species

  allocate(species)
  species%name = ''
  species%id = 0
  species%alpha = 0.d0
  species%beta = 0.d0
  species%attenuation_coef = 0.d0
  nullify(species%next)

  ReactionMnrlCreateTSTPrefSpec => species

end function ReactionMnrlCreateTSTPrefSpec

! ************************************************************************** !

function ReactionMnrlCreateMnrlConstraint(mineral,option)
  !
  ! Creates a mineral constraint object
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  !
  use Option_module

  implicit none

  type(mineral_type) :: mineral
  type(option_type) :: option
  type(mineral_constraint_type), pointer :: ReactionMnrlCreateMnrlConstraint

  type(mineral_constraint_type), pointer :: constraint

  allocate(constraint)
  allocate(constraint%names(mineral%nkinmnrl))
  constraint%names = ''
  allocate(constraint%constraint_vol_frac(mineral%nkinmnrl))
  constraint%constraint_vol_frac = UNINITIALIZED_DOUBLE
  allocate(constraint%constraint_area(mineral%nkinmnrl))
  constraint%constraint_area = UNINITIALIZED_DOUBLE
  allocate(constraint%constraint_vol_frac_string(mineral%nkinmnrl))
  constraint%constraint_vol_frac_string = ''
  allocate(constraint%constraint_area_string(mineral%nkinmnrl))
  constraint%constraint_area_string = ''
  allocate(constraint%constraint_area_units(mineral%nkinmnrl))
  constraint%constraint_area_units = ''
  allocate(constraint%constraint_area_conv_factor(mineral%nkinmnrl))
  constraint%constraint_area_conv_factor = UNINITIALIZED_DOUBLE
  allocate(constraint%area_per_unit_mass(mineral%nkinmnrl))
  constraint%area_per_unit_mass = PETSC_FALSE
  allocate(constraint%external_vol_frac_dataset(mineral%nkinmnrl))
  constraint%external_vol_frac_dataset = PETSC_FALSE
  allocate(constraint%external_area_dataset(mineral%nkinmnrl))
  constraint%external_area_dataset = PETSC_FALSE

  ReactionMnrlCreateMnrlConstraint => constraint

end function ReactionMnrlCreateMnrlConstraint

! ************************************************************************** !

function ReactionMnrlGetMnrlPtrFromName(name,mineral)
  !
  ! Returns the mineral corresponding to the name
  !
  ! Author: Glenn Hammond
  ! Date: 11/22/17
  !
  use String_module

  implicit none

  character(len=MAXWORDLENGTH) :: name
  type(mineral_type) :: mineral

  type(mineral_rxn_type), pointer :: ReactionMnrlGetMnrlPtrFromName

  ReactionMnrlGetMnrlPtrFromName => mineral%mineral_list
  do
    if (.not.associated(ReactionMnrlGetMnrlPtrFromName)) exit
    if (StringCompare(name,ReactionMnrlGetMnrlPtrFromName%name, &
                      MAXWORDLENGTH)) return
    ReactionMnrlGetMnrlPtrFromName => ReactionMnrlGetMnrlPtrFromName%next
  enddo

end function ReactionMnrlGetMnrlPtrFromName

! ************************************************************************** !

function ReactionMnrlGetMineralNames(mineral)
  !
  ! Returns the names of minerals in an array
  !
  ! Author: Glenn Hammond
  ! Date: 09/04/08
  !

  implicit none

  character(len=MAXWORDLENGTH), pointer :: ReactionMnrlGetMineralNames(:)
  type(mineral_type) :: mineral

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(mineral_rxn_type), pointer :: cur_mineral

  count = ReactionMnrlGetMineralCount(mineral)
  allocate(names(count))

  count = 1
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    names(count) = cur_mineral%name
    count = count + 1
    cur_mineral => cur_mineral%next
  enddo

  ReactionMnrlGetMineralNames => names

end function ReactionMnrlGetMineralNames

! ************************************************************************** !

function ReactionMnrlGetMineralCount(mineral)
  !
  ! Returns the number of minerals
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/08
  !

  implicit none

  PetscInt :: ReactionMnrlGetMineralCount
  type(mineral_type) :: mineral

  type(mineral_rxn_type), pointer :: cur_mineral

  ReactionMnrlGetMineralCount = 0
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    ReactionMnrlGetMineralCount = ReactionMnrlGetMineralCount + 1
    cur_mineral => cur_mineral%next
  enddo

end function ReactionMnrlGetMineralCount

! ************************************************************************** !

function ReactionMnrlGetMnrlIDFromName1(name,mineral,option)
  !
  ! Returns the id of mineral with the corresponding name
  !
  ! Author: Glenn Hammond
  ! Date: 09/04/08
  !
  use Option_module
  use String_module

  implicit none

  type(mineral_type) :: mineral
  character(len=MAXWORDLENGTH) :: name
  type(option_type) :: option

  PetscInt :: ReactionMnrlGetMnrlIDFromName1

  ReactionMnrlGetMnrlIDFromName1 = &
    ReactionMnrlGetMnrlIDFromName(name,mineral,PETSC_FALSE,PETSC_TRUE,option)

end function ReactionMnrlGetMnrlIDFromName1

! ************************************************************************** !

function ReactionMnrlGetMnrlIDFromName2(name,mineral,must_be_kinetic, &
                                        throw_error,option)
  !
  ! Returns the id of mineral with the corresponding name
  !
  ! Author: Glenn Hammond
  ! Date: 09/04/08
  !
  use Option_module
  use String_module

  implicit none

  type(mineral_type) :: mineral
  character(len=MAXWORDLENGTH) :: name
  PetscBool :: must_be_kinetic
  PetscBool :: throw_error
  type(option_type) :: option

  PetscInt :: ReactionMnrlGetMnrlIDFromName2
  type(mineral_rxn_type), pointer :: cur_mineral
  PetscInt :: ikinmnrl

  ReactionMnrlGetMnrlIDFromName2 = -1

  cur_mineral => mineral%mineral_list
  ikinmnrl = 0
  do
    if (.not.associated(cur_mineral)) exit
    if (cur_mineral%itype == MINERAL_KINETIC) then
      ikinmnrl = ikinmnrl + 1
    endif
    if (StringCompare(name,cur_mineral%name,MAXWORDLENGTH)) then
      if (must_be_kinetic .and. cur_mineral%itype /= MINERAL_KINETIC) exit
      if (must_be_kinetic) then
        ReactionMnrlGetMnrlIDFromName2 = ikinmnrl
        exit
      endif
      ReactionMnrlGetMnrlIDFromName2 = cur_mineral%id
      exit
    endif
    cur_mineral => cur_mineral%next
  enddo

  if (throw_error .and. ReactionMnrlGetMnrlIDFromName2 <= 0) then
    option%io_buffer = 'Mineral "' // trim(name) // &
      '" not found among minerals in ReactionMnrlGetMnrlIDFromName().'
    call PrintErrMsg(option)
  endif

end function ReactionMnrlGetMnrlIDFromName2

! ************************************************************************** !

function ReactionMnrlGetKinMnrlIDFromName(name,mineral,option)
  !
  ! Returns the id of kinetic mineral with the
  ! corresponding name
  !
  ! Author: Glenn Hammond
  ! Date: 03/11/13
  !
  use Option_module
  use String_module

  implicit none

  character(len=MAXWORDLENGTH) :: name
  type(mineral_type) :: mineral
  type(option_type) :: option

  PetscInt :: ReactionMnrlGetKinMnrlIDFromName

  ReactionMnrlGetKinMnrlIDFromName = &
    ReactionMnrlGetMnrlIDFromName(name,mineral,PETSC_TRUE,PETSC_FALSE,option)

  if (ReactionMnrlGetKinMnrlIDFromName <= 0) then
    option%io_buffer = 'Mineral "' // trim(name) // &
      '" not found among kinetic minerals in &
      &ReactionMnrlGetKinMnrlIDFromName().'
    call PrintErrMsg(option)
  endif

end function ReactionMnrlGetKinMnrlIDFromName

! ************************************************************************** !

subroutine ReactionMnrlDestroyMineralRxn(mineral)
  !
  ! ReactionMnrlDestoyMineral: Deallocates a mineral rxn object
  !
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  !

  implicit none

  type(mineral_rxn_type), pointer :: mineral

  if (associated(mineral%dbaserxn)) &
    call DatabaseRxnDestroy(mineral%dbaserxn)
  if (associated(mineral%tstrxn)) &
    call ReactionMnrlDestroyTSTRxn(mineral%tstrxn)

  deallocate(mineral)
  nullify(mineral)

end subroutine ReactionMnrlDestroyMineralRxn

! ************************************************************************** !

recursive subroutine ReactionMnrlDestroyTSTRxn(tstrxn)
  !
  ! Deallocates a transition state reaction
  !
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  !

  implicit none

  type(transition_state_rxn_type), pointer :: tstrxn

  if (.not.associated(tstrxn)) return

  call ReactionMnrlDestroyTSTRxn(tstrxn%next)
  call TransitionStatePrefactorDestroy(tstrxn%prefactor)

  deallocate(tstrxn)
  nullify(tstrxn)

end subroutine ReactionMnrlDestroyTSTRxn

! ************************************************************************** !

recursive subroutine TransitionStatePrefactorDestroy(prefactor)
  !
  ! Deallocates a transition state prefactor
  !
  ! Author: Glenn Hammond
  ! Date: 07/29/11
  !

  implicit none

  type(transition_state_prefactor_type), pointer :: prefactor

  if (.not.associated(prefactor)) return

  call TransitionStatePrefactorDestroy(prefactor%next)
  call TSPrefactorSpeciesDestroy(prefactor%species)

  deallocate(prefactor)
  nullify(prefactor)

end subroutine TransitionStatePrefactorDestroy

! ************************************************************************** !

recursive subroutine TSPrefactorSpeciesDestroy(species)
  !
  ! Deallocates a transition state prefactor
  !
  ! Author: Glenn Hammond
  ! Date: 08/01/11
  !

  implicit none

  type(ts_prefactor_species_type), pointer :: species

  if (.not.associated(species)) return

  call TSPrefactorSpeciesDestroy(species%next)

  deallocate(species)
  nullify(species)

end subroutine TSPrefactorSpeciesDestroy

! ************************************************************************** !

subroutine ReactionMnrlDestMnrlConstraint(constraint)
  !
  ! Destroys a mineral constraint object
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  !

  use Utility_module, only: DeallocateArray

  implicit none

  type(mineral_constraint_type), pointer :: constraint

  if (.not.associated(constraint)) return

  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_vol_frac)
  call DeallocateArray(constraint%constraint_area)
  call DeallocateArray(constraint%constraint_vol_frac_string)
  call DeallocateArray(constraint%constraint_area_string)
  call DeallocateArray(constraint%constraint_area_units)
  call DeallocateArray(constraint%constraint_area_conv_factor)
  call DeallocateArray(constraint%area_per_unit_mass)
  call DeallocateArray(constraint%external_area_dataset)
  call DeallocateArray(constraint%external_vol_frac_dataset)

  deallocate(constraint)
  nullify(constraint)

end subroutine ReactionMnrlDestMnrlConstraint

! ************************************************************************** !

subroutine ReactionMnrlDestoyMineral(mineral)
  !
  ! Deallocates a mineral object
  !
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  !

  use Utility_module, only: DeallocateArray

  implicit none

  type(mineral_type), pointer :: mineral

  type(mineral_rxn_type), pointer :: cur_mineral, prev_mineral

  if (.not.associated(mineral)) return

  ! mineral species
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    prev_mineral => cur_mineral
    cur_mineral => cur_mineral%next
    call ReactionMnrlDestroyMineralRxn(prev_mineral)
  enddo
  nullify(mineral%mineral_list)

  call DeallocateArray(mineral%mineral_names)
  call DeallocateArray(mineral%kinmnrl_names)
  call DeallocateArray(mineral%mnrl_print)
  call DeallocateArray(mineral%kinmnrl_print)
  call DeallocateArray(mineral%mnrlspecid)
  call DeallocateArray(mineral%mnrlstoich)
  call DeallocateArray(mineral%mnrlh2oid)
  call DeallocateArray(mineral%mnrlh2ostoich)
  call DeallocateArray(mineral%mnrl_logK)
  call DeallocateArray(mineral%mnrl_logKcoef)

  call DeallocateArray(mineral%kinmnrlspecid)
  call DeallocateArray(mineral%kinmnrlstoich)
  call DeallocateArray(mineral%kinmnrlh2oid)
  call DeallocateArray(mineral%kinmnrlh2ostoich)
  call DeallocateArray(mineral%kinmnrl_logK)
  call DeallocateArray(mineral%kinmnrl_logKcoef)
  call DeallocateArray(mineral%kinmnrl_rate_constant)
  call DeallocateArray(mineral%kinmnrl_molar_vol)
  call DeallocateArray(mineral%kinmnrl_molar_wt)

  call DeallocateArray(mineral%kinmnrl_num_prefactors)
  call DeallocateArray(mineral%kinmnrl_prefactor_id)
  call DeallocateArray(mineral%kinmnrl_pref_alpha)
  call DeallocateArray(mineral%kinmnrl_pref_beta)
  call DeallocateArray(mineral%kinmnrl_pref_atten_coef)
  call DeallocateArray(mineral%kinmnrl_pref_rate)
  call DeallocateArray(mineral%kinmnrl_pref_activation_energy)

  call DeallocateArray(mineral%kinmnrl_min_scale_factor)
  call DeallocateArray(mineral%kinmnrl_Temkin_const)
  call DeallocateArray(mineral%kinmnrl_affinity_power)
  call DeallocateArray(mineral%kinmnrl_affinity_threshold)
  call DeallocateArray(mineral%kinmnrl_activation_energy)
  call DeallocateArray(mineral%kinmnrl_rate_limiter)
  call DeallocateArray(mineral%kinmnrl_surf_area_vol_frac_pwr)
  call DeallocateArray(mineral%kinmnrl_surf_area_porosity_pwr)
  call DeallocateArray(mineral%kinmnrl_armor_min_names)
  call DeallocateArray(mineral%kinmnrl_armor_pwr)
  call DeallocateArray(mineral%kinmnrl_armor_crit_vol_frac)
  call DeallocateArray(mineral%kinmnrl_irreversible)

  call DeallocateArray(mineral%kinmnrl_surf_area_epsilon)
  call DeallocateArray(mineral%kinmnrl_vol_frac_epsilon)

  deallocate(mineral)
  nullify(mineral)

end subroutine ReactionMnrlDestoyMineral

end module Reaction_Mineral_Aux_module
