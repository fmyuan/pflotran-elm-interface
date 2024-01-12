module Reaction_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Generic_module
  use Reaction_Base_module
  use Reaction_Database_Aux_module
  use Reaction_Gas_Aux_module
  use Reaction_Immobile_Aux_module
  use Reaction_Isotherm_Aux_module
  use Reaction_Microbial_Aux_module
  use Reaction_Mineral_Aux_module
  use Reaction_Surface_Complexation_Aux_module

#ifdef SOLID_SOLUTION
  use Reaction_Solid_Soln_Aux_module
#endif

  implicit none

  private

  ! activity coefficients
  PetscInt, parameter, public :: ACT_COEF_FREQUENCY_OFF = 0
  PetscInt, parameter, public :: ACT_COEF_FREQUENCY_TIMESTEP = 1
  PetscInt, parameter, public :: ACT_COEF_FREQUENCY_NEWTON_ITER = 2
  PetscInt, parameter, public :: ACT_COEF_ALGORITHM_LAG = 3
  PetscInt, parameter, public :: ACT_COEF_ALGORITHM_NEWTON = 4
  PetscInt, parameter, public :: NO_BDOT = 5

  type, public :: species_idx_type
    PetscInt :: h2o_aq_id
    PetscInt :: h_ion_id
    PetscInt :: na_ion_id
    PetscInt :: cl_ion_id
    PetscInt :: co2_aq_id
    PetscInt :: tracer_aq_id
    PetscInt :: co2_gas_id
    PetscInt :: o2_gas_id
    PetscInt :: water_age_id
    PetscInt :: tracer_age_id
  end type species_idx_type

  type, public :: aq_species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: a0
    PetscReal :: molar_weight
    PetscReal :: Z
    PetscBool :: print_me
    PetscBool :: is_redox
    type(database_rxn_type), pointer :: dbaserxn
    type(aq_species_type), pointer :: next
  end type aq_species_type

  type, public :: ion_exchange_rxn_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: mineral_name
    type(ion_exchange_cation_type), pointer :: cation_list
    PetscReal :: CEC
    type(ion_exchange_rxn_type), pointer :: next
  end type ion_exchange_rxn_type

  type, public :: ion_exchange_cation_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: k
    type(ion_exchange_cation_type), pointer :: next
  end type ion_exchange_cation_type

  type, public :: dynamic_kd_rxn_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: kd_species_name
    character(len=MAXWORDLENGTH) :: ref_species_name
    PetscReal :: ref_species_high
    PetscReal :: KD_high
    PetscReal :: KD_low
    PetscReal :: KD_power
    type(dynamic_kd_rxn_type), pointer :: next
  end type dynamic_kd_rxn_type

  type, public :: radioactive_decay_rxn_type
    PetscInt :: id
    character(len=MAXSTRINGLENGTH) :: reaction
    PetscReal :: rate_constant
    PetscReal :: half_life
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn
    type(radioactive_decay_rxn_type), pointer :: next
  end type radioactive_decay_rxn_type

  type, public :: general_rxn_type
    PetscInt :: id
    character(len=MAXSTRINGLENGTH) :: reaction
    PetscReal :: forward_rate
    PetscReal :: backward_rate
    PetscBool :: print_me
    type(database_rxn_type), pointer :: dbaserxn
    type(general_rxn_type), pointer :: next
  end type general_rxn_type

  type, public :: aq_species_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    PetscReal, pointer :: basis_molarity(:)
    PetscInt, pointer :: constraint_type(:)
    PetscInt, pointer :: constraint_spec_id(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_aux_string(:)
    PetscBool, pointer :: external_dataset(:)
  end type aq_species_constraint_type

  type, public :: guess_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: conc(:)
  end type guess_constraint_type

  type, public, extends(reaction_base_type) :: reaction_rt_type
    character(len=MAXSTRINGLENGTH) :: database_filename
    PetscBool :: use_full_geochemistry
    PetscReal :: truncated_concentration
    PetscBool :: check_update
    PetscBool :: print_all_species
    PetscBool :: print_all_primary_species
    PetscBool :: print_all_secondary_species
    PetscBool :: print_pH
    PetscBool :: print_Eh
    PetscBool :: print_pe
    PetscBool :: print_O2
    PetscBool :: print_kd
    PetscBool :: print_total_sorb
    PetscBool :: print_total_sorb_mobile
    PetscBool :: print_act_coefs
    PetscBool :: print_total_component
    PetscBool :: print_free_ion
    PetscBool :: print_total_bulk ! total in aq and sorbed phases
    PetscBool :: initialize_with_molality
    PetscBool :: print_age
    PetscBool :: print_auxiliary
    PetscBool :: print_verbose_constraints
    PetscBool :: use_geothermal_hpt
    PetscInt :: print_free_conc_type
    PetscInt :: print_tot_conc_type
    PetscInt :: print_secondary_conc_type
    PetscBool :: print_total_mass_kg
    PetscInt :: num_dbase_temperatures
    PetscInt :: num_dbase_parameters

    PetscInt :: logging_verbosity
    PetscInt :: maximum_reaction_cuts
    PetscInt :: maximum_reaction_iterations
    PetscInt :: io_rank
    PetscBool :: stop_on_rreact_failure
    PetscBool :: use_total_as_guess

    PetscReal, pointer :: dbase_temperatures(:)
    type(species_idx_type), pointer :: species_idx

    type(aq_species_type), pointer :: primary_species_list
    type(aq_species_type), pointer :: secondary_species_list
    type(ion_exchange_rxn_type), pointer :: ion_exchange_rxn_list
    type(general_rxn_type), pointer :: general_rxn_list
    type(radioactive_decay_rxn_type), pointer :: radioactive_decay_rxn_list
    type(dynamic_kd_rxn_type), pointer :: dynamic_kd_rxn_list
    type(aq_species_type), pointer :: redox_species_list
    type(generic_parameter_type), pointer :: aq_diffusion_coefficients
    type(generic_parameter_type), pointer :: gas_diffusion_coefficients
    PetscInt :: act_coef_update_frequency
    PetscInt :: act_coef_update_algorithm
    PetscBool :: checkpoint_activity_coefs
    PetscBool :: act_coef_use_bdot
    PetscBool :: use_activity_h2o
    PetscBool :: calculate_water_age
    PetscBool :: calculate_tracer_age

    ! new reaction objects
    type(surface_complexation_type), pointer :: surface_complexation
    type(mineral_type), pointer :: mineral
    type(microbial_type), pointer :: microbial
    type(immobile_type), pointer :: immobile
    type(gas_type), pointer :: gas
    type(isotherm_type), pointer :: isotherm

#ifdef SOLID_SOLUTION
    type(solid_solution_type), pointer :: solid_solution_list
#endif

    ! phases
    PetscInt :: nphase

    ! compressed arrays for efficient computation
    ! primary aqueous complexes
    PetscInt :: ncomp
    PetscInt :: naqcomp
    PetscInt :: ncollcomp
    PetscInt :: nimcomp

    ! offsets
    PetscInt :: offset_aqueous
    PetscInt :: offset_immobile

    character(len=MAXWORDLENGTH), pointer :: primary_species_names(:)
    PetscBool, pointer :: primary_species_print(:)
    PetscReal, pointer :: primary_spec_a0(:)
    PetscReal, pointer :: primary_spec_Z(:)
    PetscReal, pointer :: primary_spec_molar_wt(:)

    ! aqueous complexes
    PetscInt :: neqcplx
    character(len=MAXWORDLENGTH), pointer :: secondary_species_names(:)
    PetscBool, pointer :: secondary_species_print(:)
    character(len=MAXWORDLENGTH), pointer :: eqcplx_basis_names(:,:)
    PetscBool, pointer :: eqcplx_basis_print(:)
    PetscInt, pointer :: eqcplxspecid(:,:)   ! (0:ncomp in rxn)
    PetscReal, pointer :: eqcplxstoich(:,:)
    PetscInt, pointer :: eqcplxh2oid(:)       ! id of water, if present
    PetscReal, pointer :: eqcplxh2ostoich(:)  ! stoichiometry of water, if present
    PetscReal, pointer :: eqcplx_a0(:)  ! Debye-Huckel constant
    PetscReal, pointer :: eqcplx_Z(:)
    PetscReal, pointer :: eqcplx_molar_wt(:)
    PetscReal, pointer :: eqcplx_logK(:)
    PetscReal, pointer :: eqcplx_logKcoef(:,:)
    ! Debye-Huckel
    PetscReal :: debyeA  ! Debye-Huckel A coefficient
    PetscReal :: debyeB  ! Debye-Huckel B coefficient
    PetscReal :: debyeBdot  ! Debye-Huckel Bdot coefficient

    PetscInt :: nsorb
    PetscInt :: neqsorb
    PetscBool, pointer :: kd_print(:)
    PetscBool, pointer :: total_sorb_print(:)

    ! ionx exchange reactions
    PetscInt :: neqionxrxn
    PetscInt :: neqionxcation
    PetscBool, pointer :: eqionx_rxn_Z_flag(:)
    PetscInt, pointer :: eqionx_rxn_cation_X_offset(:)
    PetscReal, pointer :: eqionx_rxn_CEC(:)
    PetscInt, pointer :: eqionx_rxn_to_surf(:)
    PetscReal, pointer :: eqionx_rxn_k(:,:)
    PetscInt, pointer :: eqionx_rxn_cationid(:,:)
#if 0
    PetscReal, pointer :: kinionx_rxn_CEC(:)
    PetscReal, pointer :: kinionx_rxn_k(:,:)
    PetscInt, pointer :: kinionx_rxn_cationid(:)
#endif

    ! radioactive decay rxn
    PetscInt :: nradiodecay_rxn
    ! ids and stoichiometries for species involved in reaction
    PetscInt, pointer :: radiodecayspecid(:,:)
    PetscReal, pointer :: radiodecaystoich(:,:)
    ! index of radiodecayspecid for species in forward
    ! reaction equation
    PetscInt, pointer :: radiodecayforwardspecid(:)
    PetscReal, pointer :: radiodecay_kf(:)

    ! general rxn
    PetscInt :: ngeneral_rxn
    ! ids and stoichiometries for species involved in reaction
    PetscInt, pointer :: generalspecid(:,:)
    PetscReal, pointer :: generalstoich(:,:)
    ! index of generalspecid & generalstoich for species in forward
    ! reaction equation
    PetscInt, pointer :: generalforwardspecid(:,:)
    PetscReal, pointer :: generalforwardstoich(:,:)
    ! index of generalspecid & generalstoich for species in backward
    ! reaction equation
    PetscInt, pointer :: generalbackwardspecid(:,:)
    PetscReal, pointer :: generalbackwardstoich(:,:)
    PetscInt, pointer :: generalh2oid(:)
    PetscReal, pointer :: generalh2ostoich(:)
    PetscReal, pointer :: general_kf(:)
    PetscReal, pointer :: general_kr(:)

    ! dynamic kd rxn
    PetscInt :: neqdynamickdrxn
    PetscInt, pointer :: eqdynamickdspecid(:)
    PetscInt, pointer :: eqdynamickdrefspecid(:)
    PetscReal, pointer :: eqdynamickdrefspechigh(:)
    PetscReal, pointer :: eqdynamickdlow(:)
    PetscReal, pointer :: eqdynamickdhigh(:)
    PetscReal, pointer :: eqdynamickdpower(:)

    PetscReal :: max_dlnC
    PetscReal :: max_dlnC_rreact
    PetscReal :: max_relative_change_tolerance
    PetscReal :: max_residual_tolerance
    PetscReal :: max_rel_residual_tolerance

    PetscBool :: update_permeability
    PetscBool :: update_tortuosity
    PetscBool :: update_porosity
    PetscBool :: calculate_initial_porosity
    PetscReal :: minimum_porosity
    PetscBool :: update_mineral_surface_area
    PetscBool :: update_mnrl_surf_with_porosity

    PetscBool :: update_armor_mineral_surface
    PetscInt :: update_armor_mineral_surface_flag

    PetscBool :: use_sandbox
    PetscInt :: nauxiliary
    PetscInt :: mc_flag

 end type reaction_rt_type


  interface ReactionGetPriSpeciesIDFromName
    module procedure ReactionGetPriSpeciesIDFromName1
    module procedure ReactionGetPriSpeciesIDFromName2
  end interface

  interface ReactionGetSecSpeciesIDFromName
    module procedure ReactionGetSecSpeciesIDFromName1
    module procedure ReactionGetSecSpeciesIDFromName2
  end interface

  public :: ReactionCreate, &
            ReactionCast, &
            SpeciesIndexCreate, &
            GasSpeciesCreate, &
            ReactionGetPriSpeciesCount, &
            ReactionGetPriSpeciesNames, &
            ReactionGetPriSpeciesIDFromName, &
            ReactionGetSecSpeciesCount, &
            ReactionGetSecSpeciesNames, &
            ReactionGetSecSpeciesIDFromName, &
            ReactionGetImmobileCount, &
            ReactionFitLogKCoef, &
            ReactionInitializeLogK, &
            ReactionInterpolateLogK, &
            ReactionInitializeLogK_hpt, &
            ReactionInterpolateLogK_hpt, &
            TransitionStateTheoryRxnCreate, &
            TransitionStatePrefactorCreate, &
            TSPrefactorSpeciesCreate, &
            TransitionStateTheoryRxnDestroy, &
            AqueousSpeciesCreate, &
            AqueousSpeciesDestroy, &
            AqueousSpeciesConstraintCreate, &
            AqueousSpeciesConstraintDestroy, &
            GuessConstraintCreate, &
            GuessConstraintDestroy, &
            MineralConstraintCreate, &
            MineralConstraintDestroy, &
            RadioactiveDecayRxnCreate, &
            RadioactiveDecayRxnDestroy, &
            GeneralRxnCreate, &
            GeneralRxnDestroy, &
            DynamicKDRxnCreate, &
            DynamicKDRxnDestroy, &
            IonExchangeRxnCreate, &
            IonExchangeCationCreate, &
            ReactionInputRecord, &
            ReactionNetworkToStoich, &
            ReactionDestroy

contains

! ************************************************************************** !

function ReactionCreate()
  !
  ! Allocate and initialize reaction object
  !
  ! Author: Glenn Hammond
  ! Date: 05/02/08
  !
  implicit none

  class(reaction_rt_type), pointer :: ReactionCreate

  class(reaction_rt_type), pointer :: reaction

  allocate(reaction)
  call ReactionBaseInit(reaction)

  reaction%database_filename = ''
  reaction%num_dbase_temperatures = 0
  nullify(reaction%dbase_temperatures)

  reaction%act_coef_use_bdot = PETSC_TRUE
  reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
  reaction%act_coef_update_algorithm = ACT_COEF_ALGORITHM_LAG
  reaction%checkpoint_activity_coefs = PETSC_TRUE
  reaction%print_all_species = PETSC_FALSE
  reaction%print_all_primary_species = PETSC_FALSE
  reaction%print_all_secondary_species = PETSC_FALSE
  reaction%print_pH = PETSC_FALSE
  reaction%print_Eh = PETSC_FALSE
  reaction%print_pe = PETSC_FALSE
  reaction%print_O2 = PETSC_FALSE
  reaction%print_kd = PETSC_FALSE
  reaction%print_total_sorb = PETSC_FALSE
  reaction%print_total_sorb_mobile = PETSC_FALSE
  reaction%print_act_coefs = PETSC_FALSE
  reaction%truncated_concentration = UNINITIALIZED_DOUBLE
  reaction%check_update = PETSC_TRUE
  reaction%use_full_geochemistry = PETSC_FALSE
  reaction%use_activity_h2o = PETSC_FALSE
  reaction%calculate_tracer_age = PETSC_FALSE
  reaction%calculate_water_age = PETSC_FALSE
  reaction%print_age = PETSC_FALSE
  reaction%print_auxiliary = PETSC_FALSE
  reaction%print_verbose_constraints = PETSC_FALSE
  reaction%print_total_component = PETSC_FALSE
  reaction%print_free_ion = PETSC_FALSE
  reaction%print_total_bulk = PETSC_FALSE
  reaction%use_geothermal_hpt = PETSC_FALSE
  reaction%print_total_mass_kg = PETSC_FALSE

  reaction%initialize_with_molality = PETSC_FALSE
  reaction%print_free_conc_type = 0
  reaction%print_tot_conc_type = 0
  reaction%print_secondary_conc_type = 0

  reaction%logging_verbosity = 0
  reaction%maximum_reaction_cuts = 10
  reaction%maximum_reaction_iterations = 20
  reaction%io_rank = UNINITIALIZED_INTEGER
  reaction%stop_on_rreact_failure = PETSC_TRUE
  reaction%use_total_as_guess = PETSC_FALSE

  nullify(reaction%species_idx)

  nullify(reaction%primary_species_list)
  nullify(reaction%secondary_species_list)
  nullify(reaction%ion_exchange_rxn_list)
  nullify(reaction%radioactive_decay_rxn_list)
  nullify(reaction%general_rxn_list)
  nullify(reaction%dynamic_kd_rxn_list)
  nullify(reaction%redox_species_list)
  nullify(reaction%aq_diffusion_coefficients)
  nullify(reaction%gas_diffusion_coefficients)

  ! new reaction objects
  reaction%surface_complexation => SurfaceComplexationCreate()
  reaction%mineral => MineralCreate()
  reaction%microbial => MicrobialCreate()
  reaction%immobile => ReactionImCreate()
  reaction%gas => GasCreate()
  reaction%isotherm => IsothermCreate()
#ifdef SOLID_SOLUTION
  nullify(reaction%solid_solution_list)
#endif

  nullify(reaction%primary_species_names)
  nullify(reaction%secondary_species_names)
  nullify(reaction%eqcplx_basis_names)

  nullify(reaction%primary_species_print)
  nullify(reaction%secondary_species_print)
  nullify(reaction%eqcplx_basis_print)
  nullify(reaction%kd_print)
  nullify(reaction%total_sorb_print)

  reaction%nphase = 0
  reaction%ncomp = 0
  reaction%naqcomp = 0
  reaction%offset_aqueous = 0
  reaction%offset_immobile = 0
  nullify(reaction%primary_spec_a0)
  nullify(reaction%primary_spec_Z)
  nullify(reaction%primary_spec_molar_wt)

  reaction%neqcplx = 0
  nullify(reaction%eqcplxspecid)
  nullify(reaction%eqcplxstoich)
  nullify(reaction%eqcplxh2oid)
  nullify(reaction%eqcplxh2ostoich)
  nullify(reaction%eqcplx_a0)
  nullify(reaction%eqcplx_Z)
  nullify(reaction%eqcplx_molar_wt)
  nullify(reaction%eqcplx_logK)
  nullify(reaction%eqcplx_logKcoef)

  reaction%debyeA = 0.5114d0
  reaction%debyeB = 0.3288d0
  reaction%debyeBdot = 0.0410d0

  reaction%nsorb = 0
  reaction%neqsorb = 0

  reaction%neqionxrxn = 0
  reaction%neqionxcation = 0
  nullify(reaction%eqionx_rxn_Z_flag)
  nullify(reaction%eqionx_rxn_cation_X_offset)
  nullify(reaction%eqionx_rxn_CEC)
  nullify(reaction%eqionx_rxn_to_surf)
  nullify(reaction%eqionx_rxn_k)
  nullify(reaction%eqionx_rxn_cationid)
#if 0
  nullify(reaction%kinionx_CEC)
  nullify(reaction%kinionx_k)
  nullify(reaction%kinionx_cationid)
#endif

  reaction%ngeneral_rxn = 0
  nullify(reaction%generalspecid)
  nullify(reaction%generalstoich)
  nullify(reaction%generalforwardspecid)
  nullify(reaction%generalforwardstoich)
  nullify(reaction%generalbackwardspecid)
  nullify(reaction%generalbackwardstoich)
  nullify(reaction%generalh2oid)
  nullify(reaction%generalh2ostoich)
  nullify(reaction%general_kf)
  nullify(reaction%general_kr)

  reaction%nradiodecay_rxn = 0
  nullify(reaction%radiodecayspecid)
  nullify(reaction%radiodecaystoich)
  nullify(reaction%radiodecayforwardspecid)
  nullify(reaction%radiodecay_kf)

  reaction%neqdynamickdrxn = 0
  nullify(reaction%eqdynamickdspecid)
  nullify(reaction%eqdynamickdrefspecid)
  nullify(reaction%eqdynamickdrefspechigh)
  nullify(reaction%eqdynamickdlow)
  nullify(reaction%eqdynamickdhigh)
  nullify(reaction%eqdynamickdpower)

  reaction%max_dlnC = 5.d0
  reaction%max_dlnC_rreact = 5.d0
  reaction%max_relative_change_tolerance = 1.d-6
  reaction%max_residual_tolerance = 1.d-12
  reaction%max_rel_residual_tolerance = 1.d-8

  reaction%update_permeability = PETSC_FALSE
  reaction%update_tortuosity = PETSC_FALSE
  reaction%update_porosity = PETSC_FALSE
  reaction%calculate_initial_porosity = PETSC_FALSE
  reaction%minimum_porosity = 0.d0
  reaction%update_mineral_surface_area = PETSC_FALSE
  reaction%update_mnrl_surf_with_porosity = PETSC_FALSE

  reaction%update_armor_mineral_surface = PETSC_FALSE
  reaction%update_armor_mineral_surface_flag = 0

  reaction%use_sandbox = PETSC_FALSE
  reaction%nauxiliary = 0
  reaction%mc_flag = 0

  ReactionCreate => reaction

end function ReactionCreate

! ************************************************************************** !

function ReactionCast(reaction_base)
  !
  ! Casts a reaction_base type to reaction_nw type if applicable.
  !
  ! Author: Glenn Hammond
  ! Date: 10/21/19
  !
  implicit none

  class(reaction_base_type), pointer :: reaction_base

  class(reaction_rt_type), pointer :: ReactionCast

  nullify(ReactionCast)
  if (.not.associated(reaction_base)) return
  select type(r=>reaction_base)
    class is(reaction_rt_type)
      ReactionCast => r
  end select

end function ReactionCast

! ************************************************************************** !

function SpeciesIndexCreate()
  !
  ! Allocate and initialize a species index object
  !
  ! Author: Peter Lichtner
  ! Date: 01/29/10
  !

  implicit none

  type(species_idx_type), pointer :: SpeciesIndexCreate

  type(species_idx_type), pointer :: species_idx

  allocate(species_idx)

  species_idx%h2o_aq_id = 0
  species_idx%h_ion_id = 0
  species_idx%na_ion_id = 0
  species_idx%cl_ion_id = 0
  species_idx%co2_aq_id = 0
  species_idx%tracer_aq_id = 0
  species_idx%co2_gas_id = 0
  species_idx%o2_gas_id = 0
  species_idx%tracer_age_id = 0
  species_idx%water_age_id = 0

  SpeciesIndexCreate => species_idx

end function SpeciesIndexCreate

! ************************************************************************** !

function AqueousSpeciesCreate()
  !
  ! Allocate and initialize an aqueous species object
  !
  ! Author: Glenn Hammond
  ! Date: 05/02/08
  !
  implicit none

  type(aq_species_type), pointer :: AqueousSpeciesCreate

  type(aq_species_type), pointer :: species

  allocate(species)
  species%id = 0
  species%name = ''
  species%a0 = 0.d0
  species%molar_weight = 0.d0
  species%Z = 0.d0
  species%print_me = PETSC_FALSE
  species%is_redox = PETSC_FALSE
  nullify(species%dbaserxn)
  nullify(species%next)

  AqueousSpeciesCreate => species

end function AqueousSpeciesCreate

! ************************************************************************** !

function IonExchangeRxnCreate()
  !
  ! Allocate and initialize an ion exchange reaction
  !
  ! Author: Peter Lichtner
  ! Date: 10/24/08
  !

  implicit none

  type(ion_exchange_rxn_type), pointer :: IonExchangeRxnCreate

  type(ion_exchange_rxn_type), pointer :: ionxrxn

  allocate(ionxrxn)
  ionxrxn%id = 0
  ionxrxn%mineral_name = ''
  ionxrxn%CEC = UNINITIALIZED_DOUBLE
  nullify(ionxrxn%cation_list)
  nullify(ionxrxn%next)

  IonExchangeRxnCreate => ionxrxn

end function IonExchangeRxnCreate

! ************************************************************************** !

function IonExchangeCationCreate()
  !
  ! Allocate and initialize a cation associated with
  ! an ion exchange reaction
  !
  ! Author: Peter Lichtner
  ! Date: 10/24/08
  !

  implicit none

  type(ion_exchange_cation_type), pointer :: IonExchangeCationCreate

  type(ion_exchange_cation_type), pointer :: cation

  allocate(cation)
  cation%name = ''
  cation%k = 0.d0
  nullify(cation%next)

  IonExchangeCationCreate => cation

end function IonExchangeCationCreate

! ************************************************************************** !

function RadioactiveDecayRxnCreate()
  !
  ! Allocate and initialize a radioactive decay
  ! reaction
  !
  ! Author: Glenn Hammond
  ! Date: 01/07/14
  !
  implicit none

  type(radioactive_decay_rxn_type), pointer :: RadioactiveDecayRxnCreate

  type(radioactive_decay_rxn_type), pointer :: rxn

  allocate(rxn)
  rxn%id = 0
  rxn%reaction = ''
  rxn%rate_constant = 0.d0
  rxn%half_life = 0.d0
  rxn%print_me = PETSC_FALSE
  nullify(rxn%dbaserxn)
  nullify(rxn%next)

  RadioactiveDecayRxnCreate => rxn

end function RadioactiveDecayRxnCreate

! ************************************************************************** !

function GeneralRxnCreate()
  !
  ! Allocate and initialize a general reaction
  !
  ! Author: Glenn Hammond
  ! Date: 09/03/10
  !
  implicit none

  type(general_rxn_type), pointer :: GeneralRxnCreate

  type(general_rxn_type), pointer :: rxn

  allocate(rxn)
  rxn%id = 0
  rxn%reaction = ''
  rxn%forward_rate = 0.d0
  rxn%backward_rate = 0.d0
  rxn%print_me = PETSC_FALSE
  nullify(rxn%dbaserxn)
  nullify(rxn%next)

  GeneralRxnCreate => rxn

end function GeneralRxnCreate

! ************************************************************************** !

function DynamicKDRxnCreate()
  !
  ! Allocate and initialize a dynamic KD sorption reaction
  !
  ! Author: Glenn Hammond
  ! Date: 12/21/19
  !

  implicit none

  type(dynamic_kd_rxn_type), pointer :: DynamicKDRxnCreate

  type(dynamic_kd_rxn_type), pointer :: rxn

  allocate(rxn)
  rxn%id = 0
  rxn%kd_species_name = ''
  rxn%ref_species_name = ''
  rxn%ref_species_high = UNINITIALIZED_DOUBLE
  rxn%KD_low = UNINITIALIZED_DOUBLE
  rxn%KD_high = UNINITIALIZED_DOUBLE
  rxn%KD_power = UNINITIALIZED_DOUBLE
  nullify(rxn%next)

  DynamicKDRxnCreate => rxn

end function DynamicKDRxnCreate

! ************************************************************************** !

function AqueousSpeciesConstraintCreate(reaction,option)
  !
  ! Creates an aqueous species constraint
  ! object
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  !
  use Option_module

  implicit none

  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  type(aq_species_constraint_type), pointer :: AqueousSpeciesConstraintCreate

  type(aq_species_constraint_type), pointer :: constraint

  allocate(constraint)
  allocate(constraint%names(reaction%naqcomp))
  constraint%names = ''
  allocate(constraint%constraint_conc(reaction%naqcomp))
  constraint%constraint_conc = 0.d0
  allocate(constraint%basis_molarity(reaction%naqcomp))
  constraint%basis_molarity = 0.d0
  allocate(constraint%constraint_spec_id(reaction%naqcomp))
  constraint%constraint_spec_id = 0
  allocate(constraint%constraint_type(reaction%naqcomp))
  constraint%constraint_type = 0
  allocate(constraint%constraint_aux_string(reaction%naqcomp))
  constraint%constraint_aux_string = ''
  allocate(constraint%external_dataset(reaction%naqcomp))
  constraint%external_dataset = PETSC_FALSE

  AqueousSpeciesConstraintCreate => constraint

end function AqueousSpeciesConstraintCreate

! ************************************************************************** !

function GuessConstraintCreate(reaction,option)
  !
  ! Creates an aqueous species constraint
  ! object
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  !

  use Option_module

  implicit none

  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  type(guess_constraint_type), pointer :: GuessConstraintCreate

  type(guess_constraint_type), pointer :: constraint

  allocate(constraint)
  allocate(constraint%names(reaction%naqcomp))
  constraint%names = ''
  allocate(constraint%conc(reaction%naqcomp))
  constraint%conc = 0.d0

  GuessConstraintCreate => constraint

end function GuessConstraintCreate

! ************************************************************************** !

function ReactionGetPriSpeciesNames(reaction)
  !
  ! Returns the names of primary species in an array
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/08
  !

  implicit none

  character(len=MAXWORDLENGTH), pointer :: ReactionGetPriSpeciesNames(:)
  class(reaction_rt_type) :: reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(aq_species_type), pointer :: species

  count = ReactionGetPriSpeciesCount(reaction)
  allocate(names(count))

  count = 1
  species => reaction%primary_species_list
  do
    if (.not.associated(species)) exit
    names(count) = species%name
    count = count + 1
    species => species%next
  enddo

  ReactionGetPriSpeciesNames => names

end function ReactionGetPriSpeciesNames

! ************************************************************************** !

function ReactionGetPriSpeciesCount(reaction)
  !
  ! Returns the number of primary species
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/08
  !

  implicit none

  PetscInt :: ReactionGetPriSpeciesCount
  class(reaction_rt_type) :: reaction

  type(aq_species_type), pointer :: species

  ReactionGetPriSpeciesCount = 0
  species => reaction%primary_species_list
  do
    if (.not.associated(species)) exit
    ReactionGetPriSpeciesCount = ReactionGetPriSpeciesCount + 1
    species => species%next
  enddo

end function ReactionGetPriSpeciesCount

! ************************************************************************** !

function ReactionGetPriSpeciesIDFromName1(name,reaction,option)
  !
  ! Returns the id of named primary species
  !
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  !
  use Option_module
  use String_module

  implicit none

  character(len=MAXWORDLENGTH) :: name
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscInt :: ReactionGetPriSpeciesIDFromName1

  ReactionGetPriSpeciesIDFromName1 = &
    ReactionGetPriSpeciesIDFromName2(name,reaction,PETSC_TRUE, option)

end function ReactionGetPriSpeciesIDFromName1

! ************************************************************************** !

function ReactionGetPriSpeciesIDFromName2(name,reaction,return_error,option)
  !
  ! Returns the id of named primary species
  !
  ! Author: Glenn Hammond
  ! Date: 10/30/12


  use Option_module
  use String_module

  implicit none

  character(len=MAXWORDLENGTH) :: name
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscInt :: ReactionGetPriSpeciesIDFromName2

  type(aq_species_type), pointer :: species
  PetscInt :: i
  PetscBool :: return_error

  ReactionGetPriSpeciesIDFromName2 = UNINITIALIZED_INTEGER

  ! if the primary species name list exists
  if (associated(reaction%primary_species_names)) then
    do i = 1, size(reaction%primary_species_names)
      if (StringCompare(name,reaction%primary_species_names(i), &
                        MAXWORDLENGTH)) then
        ReactionGetPriSpeciesIDFromName2 = i
        exit
      endif
    enddo
  else
    species => reaction%primary_species_list
    i = 0
    do
      if (.not.associated(species)) exit
      i = i + 1
      if (StringCompare(name,species%name,MAXWORDLENGTH)) then
        ReactionGetPriSpeciesIDFromName2 = i
        exit
      endif
      species => species%next
    enddo
  endif

  if (return_error .and. ReactionGetPriSpeciesIDFromName2 <= 0) then
    option%io_buffer = 'Species "' // trim(name) // &
      '" not found among primary species in ReactionGetPriSpeciesIDFromName2().'
    call PrintErrMsg(option)
  endif

end function ReactionGetPriSpeciesIDFromName2

! ************************************************************************** !

function ReactionGetSecSpeciesNames(reaction)
  !
  ! Returns the names of secondary species in an array
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/08
  !

  implicit none

  character(len=MAXWORDLENGTH), pointer :: ReactionGetSecSpeciesNames(:)
  class(reaction_rt_type) :: reaction

  PetscInt :: count
  character(len=MAXWORDLENGTH), pointer :: names(:)
  type(aq_species_type), pointer :: species

  count = ReactionGetSecSpeciesCount(reaction)
  allocate(names(count))

  count = 1
  species => reaction%secondary_species_list
  do
    if (.not.associated(species)) exit
    names(count) = species%name
    count = count + 1
    species => species%next
  enddo

  ReactionGetSecSpeciesNames => names

end function ReactionGetSecSpeciesNames

! ************************************************************************** !

function ReactionGetSecSpeciesCount(reaction)
  !
  ! Returns the number of secondary species
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/08
  !

  implicit none

  PetscInt :: ReactionGetSecSpeciesCount
  class(reaction_rt_type) :: reaction

  type(aq_species_type), pointer :: species

  ReactionGetSecSpeciesCount = 0
  species => reaction%secondary_species_list
  do
    if (.not.associated(species)) exit
    ReactionGetSecSpeciesCount = ReactionGetSecSpeciesCount + 1
    species => species%next
  enddo

end function ReactionGetSecSpeciesCount

! ************************************************************************** !

function ReactionGetSecSpeciesIDFromName1(name,reaction,option)
  !
  ! Returns the id of named secondary species
  !
  ! Author: Peter Rieke
  ! Date: 09/16/2016
  !
  use Option_module
  use String_module
  implicit none
  character(len=MAXWORDLENGTH) :: name
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscInt :: ReactionGetSecSpeciesIDFromName1
  ReactionGetSecSpeciesIDFromName1 = &
    ReactionGetSecSpeciesIDFromName2(name,reaction, PETSC_TRUE, option)

end function ReactionGetSecSpeciesIDFromName1

! ************************************************************************** !

function ReactionGetSecSpeciesIDFromName2(name,reaction,return_error,option)
  !
  ! Returns the id of named secondary species
  !
  ! Author: Peter Rieke
  ! Date: 09/16/2016
  !
  use Option_module
  use String_module
  implicit none
  character(len=MAXWORDLENGTH) :: name
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  PetscInt :: ReactionGetSecSpeciesIDFromName2
  type(aq_species_type), pointer :: species
  PetscInt :: i
  PetscBool :: return_error

  ReactionGetSecSpeciesIDFromName2 = UNINITIALIZED_INTEGER

  ! if the Secondary species name list exists
  if (associated(reaction%Secondary_species_names)) then
    do i = 1, size(reaction%Secondary_species_names)
      if (StringCompare(name,reaction%Secondary_species_names(i), &
                        MAXWORDLENGTH)) then
        ReactionGetSecSpeciesIDFromName2 = i
        exit
      endif
    enddo
  else
    species => reaction%Secondary_species_list
    i = 0
    do
      if (.not.associated(species)) exit
      i = i + 1
      if (StringCompare(name,species%name,MAXWORDLENGTH)) then
        ReactionGetSecSpeciesIDFromName2 = i
        exit
      endif
      species => species%next
    enddo
  endif

  if (return_error .and. ReactionGetSecSpeciesIDFromName2 <= 0) then
    option%io_buffer = 'Species "' // trim(name) // '" not found among &
      &secondary species in ReactionGetSecSpeciesIDFromName().'
    call PrintErrMsg(option)
  endif

end function ReactionGetSecSpeciesIDFromName2

! ************************************************************************** !

function ReactionGetImmobileCount(reaction)
  !
  ! Returns the number of immobile species
  !
  ! Author: Glenn Hammond
  ! Date: 01/02/13
  !

  implicit none

  PetscInt :: ReactionGetImmobileCount
  class(reaction_rt_type) :: reaction

  ReactionGetImmobileCount = ReactionImGetCount(reaction%immobile)

end function ReactionGetImmobileCount

! ************************************************************************** !

subroutine ReactionFitLogKCoef(coefs,logK,name,option,reaction)
  !
  ! Least squares fit to log K over database temperature
  ! range
  !
  ! Author: P.C. Lichtner
  ! Date: 02/13/09
  !

  use Option_module
  use Utility_module

  implicit none

  class(reaction_rt_type) :: reaction
  PetscReal :: coefs(FIVE_INTEGER)
  character(len=MAXWORDLENGTH) :: name
  PetscReal :: logK(reaction%num_dbase_temperatures)
  type(option_type) :: option

  PetscInt :: temp_int(reaction%num_dbase_temperatures), &
              indx(reaction%num_dbase_temperatures)
  PetscReal :: a(FIVE_INTEGER,FIVE_INTEGER), &
               vec(FIVE_INTEGER,reaction%num_dbase_temperatures), &
               temperature_kelvin

  PetscInt :: i, j, k, iflag

  ! need to fill in vec with equations for temperatures vs coefs.

  do i = 1, reaction%num_dbase_temperatures
    temperature_kelvin = reaction%dbase_temperatures(i) + 273.15d0
    vec(1,i) = log(temperature_kelvin)
    vec(2,i) = 1.d0
    vec(3,i) = temperature_kelvin
    vec(4,i) = 1.d0/temperature_kelvin
    vec(5,i) = 1.d0/(temperature_kelvin*temperature_kelvin)
  enddo

  iflag = 0
  do j = 1, FIVE_INTEGER
    coefs(j) = 0.d0
    do i = 1, reaction%num_dbase_temperatures
      if (dabs(logK(i) - 500.) < 1.d-10) then
        iflag = 1
        temp_int(i) = ZERO_INTEGER
      else
        coefs(j) = coefs(j) + vec(j,i)*logK(i)
        temp_int(i) = ONE_INTEGER
      endif
    enddo
  enddo

  if (iflag == 1) then
    option%io_buffer = 'In ReactionFitLogKCoef: log K = 500 for ' // trim(name)
    call PrintWrnMsg(option)
  endif

  do j = 1, FIVE_INTEGER
    do k = j, FIVE_INTEGER
      a(j,k) = 0.d0
      do i = 1, reaction%num_dbase_temperatures
        if (temp_int(i) == 1) then
          a(j,k) = a(j,k) + vec(j,i)*vec(k,i)
        endif
      enddo
      if (j .ne. k) a(k,j) = a(j,k)
    enddo
  enddo

  call LUDecomposition(a,FIVE_INTEGER,indx,i)
  call LUBackSubstitution(a,FIVE_INTEGER,indx,coefs)

end subroutine ReactionFitLogKCoef

! ************************************************************************** !

subroutine ReactionInitializeLogK(logKcoef,logKs,logK,option,reaction)
  !
  ! Least squares fit to log K over database temperature range
  !
  ! Author: P.C. Lichtner
  ! Date: 02/13/09
  !

  use Option_module

  implicit none

  class(reaction_rt_type) :: reaction
  PetscReal :: logKcoef(FIVE_INTEGER)
  PetscReal :: logKs(reaction%num_dbase_temperatures)
  PetscReal :: logK, logK_1D_Array(ONE_INTEGER)
  type(option_type) :: option

  PetscReal :: coefs(FIVE_INTEGER,ONE_INTEGER)
  PetscReal :: temperature
  PetscInt :: itemperature
  PetscInt :: i

  ! we always initialize on reference temperature
  temperature = option%flow%reference_temperature

  itemperature = 0
  if (option%use_isothermal) then ! find database temperature if relevant
    do i = 1, reaction%num_dbase_temperatures
      if (dabs(option%flow%reference_temperature - &
               reaction%dbase_temperatures(i)) < 1.d-10) then
        itemperature = i
        exit
      endif
    enddo
  endif

  if (itemperature > 0) then ! use database temperature
    logK = logKs(itemperature)
  else                       ! interpolate
    coefs(:,ONE_INTEGER) = logKcoef(:)
    call ReactionInterpolateLogK(coefs,logK_1D_Array,temperature,ONE_INTEGER)
    logK = logK_1D_Array(ONE_INTEGER)
  endif

end subroutine ReactionInitializeLogK

! ************************************************************************** !

subroutine ReactionInterpolateLogK(coefs,logKs,temp,n)
  !
  ! Interpolation log K function: temp - temperature [C]
  ! b - fit coefficients determined from fit(...)
  !
  ! Author: P.C. Lichtner
  ! Date: 02/13/09
  !

  implicit none

  PetscInt :: n
  PetscReal :: coefs(5,n), logKs(n), temp

  PetscInt :: i
  PetscReal :: temp_kelvin

  temp_kelvin = temp + 273.15d0

  do i = 1, n
    logKs(i) = coefs(1,i)*log(temp_kelvin) &
             + coefs(2,i)           &
             + coefs(3,i)*temp_kelvin      &
             + coefs(4,i)/temp_kelvin      &
             + coefs(5,i)/(temp_kelvin*temp_kelvin)
  enddo

end subroutine ReactionInterpolateLogK

! ************************************************************************** !

subroutine ReactionInitializeLogK_hpt(logKcoef,logK,option,reaction)
  !
  ! ReactionInitializeLogK: Least squares fit to log K over database temperature range
  !
  ! Author: Chuan Lu
  ! Date: 12/29/11
  !

  use Option_module

  implicit none

  class(reaction_rt_type) :: reaction
  PetscReal :: logKcoef(17)
  PetscReal :: logK, logK_1D_Array(ONE_INTEGER)
  type(option_type) :: option

  PetscReal :: coefs(17,ONE_INTEGER)
  PetscReal :: temperature, pressure

  ! we always initialize on reference temperature
  temperature = option%flow%reference_temperature
  pressure = option%flow%reference_pressure


  coefs(:,ONE_INTEGER) = logKcoef(:)
  call ReactionInterpolateLogK_hpt(coefs,logK_1D_Array,temperature,pressure, &
                                   ONE_INTEGER)
  logK = logK_1D_Array(ONE_INTEGER)
!   print *,'ReactionInitializeLogK_hpt: ', pressure,temperature, logK

end subroutine ReactionInitializeLogK_hpt

! ************************************************************************** !

subroutine ReactionInterpolateLogK_hpt(coefs,logKs,temp,pres,n)
  !
  ! ReactionInterpolateLogK: Interpolation log K function: temp - temperature [C]
  ! b - fit coefficients determined from fit(...)
  !
  ! Author: P.C. Lichtner
  ! Date: 02/13/09
  !

  implicit none

  PetscInt :: n
  PetscReal :: coefs(17,n), logKs(n), temp, pres

  PetscInt :: i
  PetscReal :: temp_kelvin, tr, pr, logtr

  temp_kelvin = temp + 273.15d0
  tr = temp_kelvin/273.15d0
  pr = pres/1.d7
  logtr = log(tr)/log(10.d0)

  do i = 1, n
    logKs(i) = coefs(1,i)                 &
             + coefs(2,i) * tr            &
             + coefs(3,i) / tr            &
             + coefs(4,i) * logtr         &
             + coefs(5,i) * tr * tr       &
             + coefs(6,i) / tr / tr       &
             + coefs(7,i) * sqrt(tr)      &
             + coefs(8,i) * pr            &
             + coefs(9,i) * pr * tr       &
             + coefs(10,i) * pr / tr      &
             + coefs(11,i) * pr * logtr   &
             + coefs(12,i) / pr           &
             + coefs(13,i) / pr * tr      &
             + coefs(14,i) / pr / tr      &
             + coefs(15,i) * pr * pr      &
             + coefs(16,i) * pr * pr * tr &
             + coefs(17,i) * pr * pr / tr
  enddo
 ! print *,'ReactionInterpolateLogK_hpt: ', pres,temp, logKs, coefs
end subroutine ReactionInterpolateLogK_hpt

! ************************************************************************** !

subroutine ReactionInputRecord(rxn)
  !
  ! Prints ingested chemistry and reactive transport information to the input
  ! record file.
  !
  ! Author: Jenn Frederick
  ! Date: 04/12/2016
  !
  use Reaction_Immobile_Aux_module

  implicit none

  class(reaction_rt_type), pointer :: rxn

  type(aq_species_type), pointer :: cur_aq_species
  type(gas_species_type), pointer :: cur_gas_species
  type(immobile_species_type), pointer :: cur_imm_species
  type(radioactive_decay_rxn_type), pointer :: cur_rad_decay_rxn
  type(isotherm_link_type), pointer :: cur_isotherm_rxn
  character(len=MAXWORDLENGTH) :: word1, word2
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
       &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'CHEMISTRY'

  if (.not.associated(rxn)) return

! --------- primary species list ---------------------------------------------
  if (associated(rxn%primary_species_list)) then
    write(id,'(a29)',advance='no') 'primary species list: '
    cur_aq_species => rxn%primary_species_list
    write(id,'(a)') trim(cur_aq_species%name)
    cur_aq_species => cur_aq_species%next
    do
      if (.not.associated(cur_aq_species)) exit
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') trim(cur_aq_species%name)
      cur_aq_species => cur_aq_species%next
    enddo
    write(id,'(a29)') '---------------------------: '
  endif
! --------- secondary species list -------------------------------------------
  if (associated(rxn%secondary_species_list)) then
    write(id,'(a29)',advance='no') 'secondary species list: '
    cur_aq_species => rxn%secondary_species_list
    write(id,'(a)') trim(cur_aq_species%name)
    cur_aq_species => cur_aq_species%next
    do
      if (.not.associated(cur_aq_species)) exit
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') trim(cur_aq_species%name)
      cur_aq_species => cur_aq_species%next
    enddo
    write(id,'(a29)') '---------------------------: '
  endif
! --------- active gas species list ------------------------------------------
  if (associated(rxn%gas%list)) then
    write(id,'(a29)',advance='no') 'active gas species list: '
    cur_gas_species => rxn%gas%list
    write(id,'(a)') trim(cur_gas_species%name)
    cur_gas_species => cur_gas_species%next
    do
      if (.not.associated(cur_gas_species)) exit
      if (cur_gas_species%itype == ACTIVE_GAS .or. &
          cur_gas_species%itype == ACTIVE_AND_PASSIVE_GAS) then
        write(id,'(a29)',advance='no') ' '
        write(id,'(a)') trim(cur_gas_species%name)
      endif
      cur_gas_species => cur_gas_species%next
    enddo
    write(id,'(a29)') '---------------------------: '
  endif
! --------- passive gas species list -----------------------------------------
  if (associated(rxn%gas%list)) then
    write(id,'(a29)',advance='no') 'passive gas species list: '
    cur_gas_species => rxn%gas%list
    write(id,'(a)') trim(cur_gas_species%name)
    cur_gas_species => cur_gas_species%next
    do
      if (.not.associated(cur_gas_species)) exit
      if (cur_gas_species%itype == PASSIVE_GAS .or. &
          cur_gas_species%itype == ACTIVE_AND_PASSIVE_GAS) then
        write(id,'(a29)',advance='no') ' '
        write(id,'(a)') trim(cur_gas_species%name)
      endif
      cur_gas_species => cur_gas_species%next
    enddo
    write(id,'(a29)') '---------------------------: '
  endif
! --------- immobile species list --------------------------------------------
  if (associated(rxn%immobile%list)) then
    write(id,'(a29)',advance='no') 'immobile species list: '
    cur_imm_species => rxn%immobile%list
    write(id,'(a)') trim(cur_imm_species%name)
    cur_imm_species => cur_imm_species%next
    do
      if (.not.associated(cur_imm_species)) exit
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') trim(cur_imm_species%name)
      cur_imm_species => cur_imm_species%next
    enddo
    write(id,'(a29)') '---------------------------: '
  endif

! --------- radioactive decay reaction list ----------------------------------
  if (associated(rxn%radioactive_decay_rxn_list)) then
    cur_rad_decay_rxn => rxn%radioactive_decay_rxn_list
    do
      if (.not.associated(cur_rad_decay_rxn)) exit
      write(id,'(a29)',advance='no') 'radioactive decay reaction: '
      write(id,'(a)') adjustl(trim(cur_rad_decay_rxn%reaction))
      write(id,'(a29)',advance='no') 'decay rate: '
      write(word1,*) cur_rad_decay_rxn%rate_constant
      write(id,'(a)') adjustl(trim(word1)) // ' 1/sec'

      write(id,'(a29)') '---------------------------: '
      cur_rad_decay_rxn => cur_rad_decay_rxn%next
    enddo
  endif

! --------- sorption isotherm reaction list ----------------------------------
  if (associated(rxn%isotherm%isotherm_list)) then
    cur_isotherm_rxn => rxn%isotherm%isotherm_list
    do
      if (.not.associated(cur_isotherm_rxn)) exit
      write(id,'(a29)',advance='no') 'sorption, isotherm reaction: '
      write(id,'(a)') adjustl(trim(cur_isotherm_rxn%species_name))
      write(id,'(a29)',advance='no') 'type: '
      select case (cur_isotherm_rxn%itype)
        case (SORPTION_LINEAR)
          write(id,'(a)') 'linear sorption'
        case (SORPTION_LANGMUIR)
          write(id,'(a)') 'langmuir sorption'
          write(id,'(a29)',advance='no') 'langmuir b: '
          write(word1,*) cur_isotherm_rxn%Langmuir_B
          write(id,'(a)') adjustl(trim(word1))
        case (SORPTION_FREUNDLICH)
          write(id,'(a)') 'freundlich sorption'
          write(id,'(a29)',advance='no') 'freundlich n: '
          write(word1,*) cur_isotherm_rxn%Freundlich_N
          write(id,'(a)') adjustl(trim(word1))
      end select
      if (len_trim(cur_isotherm_rxn%kd_mineral_name) > 0) then  !UPDATE
        write(id,'(a29)',advance='no') 'Kd mineral name: '
        write(id,'(a)') adjustl(trim(cur_isotherm_rxn%kd_mineral_name))
        word2 = ' L/kg'
      else
        word2 = ' kg/m^3'
      endif
      write(id,'(a29)',advance='no') 'distribution coeff. / Kd: '
      write(word1,*) cur_isotherm_rxn%Kd
      write(id,'(a)') adjustl(trim(word1)) // adjustl(trim(word2))

      write(id,'(a29)') '---------------------------: '
      cur_isotherm_rxn => cur_isotherm_rxn%next
    enddo
  endif

end subroutine ReactionInputRecord

! ************************************************************************** !

subroutine ReactionNetworkToStoich(reaction,filename,spec_ids,stoich,option)

  ! Reads in a reaction network and parses the stoichiometries and species ids

  ! Authors: Glenn Hammond
  ! Date: 05/24/22

  use Input_Aux_module
  use Option_module
  use Reaction_Database_Aux_module

  implicit none

  class(reaction_rt_type) :: reaction
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt, pointer :: spec_ids(:,:)
  PetscReal, pointer :: stoich(:,:)
  type(option_type) :: option

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: irxn, nrxn
  type(database_rxn_ptr_type), pointer :: cur_rxn, rxn_list, last_rxn

  input => InputCreate(IUNIT_TEMP,filename,option)
  input%ierr = 0
  call InputPushBlock(input,option)
  nullify(rxn_list)
  nrxn = 0
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    cur_rxn => DatabaseRxnPtrCreate()
    string = input%buf
    cur_rxn%dbaserxn => &
      DatabaseRxnCreateFromRxnString(string, &
                                     reaction%naqcomp, &
                                     reaction%offset_aqueous, &
                                     reaction%primary_species_names, &
                                     reaction%nimcomp, &
                                     reaction%offset_immobile, &
                                     reaction%immobile%names, &
                                     PETSC_FALSE,option)
    if (.not.associated(rxn_list)) then
      rxn_list => cur_rxn
    else
      last_rxn%next => cur_rxn
    endif
    last_rxn => cur_rxn
    nrxn = nrxn + 1
  enddo
  call InputPopBlock(input,option)
  call InputDestroy(input)

  allocate(spec_ids(0:reaction%naqcomp,nrxn))
  allocate(stoich(reaction%naqcomp,nrxn))
  spec_ids = 0
  stoich = 0.d0

  cur_rxn => rxn_list
  irxn = 0
  do
    if (.not.associated(cur_rxn)) exit
    irxn = irxn + 1
    spec_ids(0,irxn) = cur_rxn%dbaserxn%nspec
    spec_ids(1:cur_rxn%dbaserxn%nspec,irxn) = &
      cur_rxn%dbaserxn%spec_ids(:)
    stoich(1:cur_rxn%dbaserxn%nspec,irxn) = cur_rxn%dbaserxn%stoich(:)
    cur_rxn => cur_rxn%next
  enddo

  call DatabaseRxnPtrDestroy(rxn_list)

end subroutine ReactionNetworkToStoich

! ************************************************************************** !

subroutine SpeciesIndexDestroy(species_idx)
  !
  ! Deallocates a species index object
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/10
  !

  implicit none

  type(species_idx_type), pointer :: species_idx

  if (associated(species_idx)) deallocate(species_idx)
  nullify(species_idx)

end subroutine SpeciesIndexDestroy

! ************************************************************************** !

subroutine AqueousSpeciesDestroy(species)
  !
  ! Deallocates an aqueous species
  !
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  !

  implicit none

  type(aq_species_type), pointer :: species

  if (associated(species%dbaserxn)) call DatabaseRxnDestroy(species%dbaserxn)
  deallocate(species)
  nullify(species)

end subroutine AqueousSpeciesDestroy

! ************************************************************************** !

subroutine AqueousSpeciesListDestroy(aq_species_list)
  !
  ! Deallocates an aqueous species
  !
  ! Author: Glenn Hammond
  ! Date: 09/03/10
  !

  !TODO(geh): make these destructors recursive
  implicit none

  type(aq_species_type), pointer :: aq_species_list

  type(aq_species_type), pointer :: species, prev_species

  species => aq_species_list
  do
    if (.not.associated(species)) exit
    prev_species => species
    species => species%next
    call AqueousSpeciesDestroy(prev_species)
  enddo
  nullify(aq_species_list)

end subroutine AqueousSpeciesListDestroy

! ************************************************************************** !

subroutine IonExchangeRxnDestroy(ionxrxn)
  !
  ! Deallocates an ion exchange reaction
  !
  ! Author: Glenn Hammond
  ! Date: 10/24/08
  !

  implicit none

  type(ion_exchange_rxn_type), pointer :: ionxrxn

  type(ion_exchange_cation_type), pointer :: cur_cation, prev_cation

  if (.not.associated(ionxrxn)) return

  cur_cation => ionxrxn%cation_list
  do
    if (.not.associated(cur_cation)) exit
    prev_cation => cur_cation
    cur_cation => cur_cation%next
    deallocate(prev_cation)
    nullify(prev_cation)
  enddo

  nullify(ionxrxn%next)

  deallocate(ionxrxn)
  nullify(ionxrxn)

end subroutine IonExchangeRxnDestroy

! ************************************************************************** !

subroutine RadioactiveDecayRxnDestroy(rxn)
  !
  ! Deallocates a general reaction
  !
  ! Author: Glenn Hammond
  ! Date: 01/07/14
  !

  implicit none

  type(radioactive_decay_rxn_type), pointer :: rxn

  if (.not.associated(rxn)) return

  if (associated(rxn%dbaserxn)) &
    call DatabaseRxnDestroy(rxn%dbaserxn)
  nullify(rxn%dbaserxn)
  nullify(rxn%next)

  deallocate(rxn)
  nullify(rxn)

end subroutine RadioactiveDecayRxnDestroy

! ************************************************************************** !

subroutine GeneralRxnDestroy(rxn)
  !
  ! Deallocates a general reaction
  !
  ! Author: Glenn Hammond
  ! Date: 09/03/10
  !

  implicit none

  type(general_rxn_type), pointer :: rxn

  if (.not.associated(rxn)) return

  if (associated(rxn%dbaserxn)) &
    call DatabaseRxnDestroy(rxn%dbaserxn)
  nullify(rxn%dbaserxn)
  nullify(rxn%next)

  deallocate(rxn)
  nullify(rxn)

end subroutine GeneralRxnDestroy

! ************************************************************************** !

subroutine DynamicKDRxnDestroy(rxn)
  !
  ! Deallocates a dynamic KD reaction
  !
  ! Author: Glenn Hammond
  ! Date: 12/21/19
  !

  implicit none

  type(dynamic_kd_rxn_type), pointer :: rxn

  if (.not.associated(rxn)) return

  deallocate(rxn)
  nullify(rxn)

end subroutine DynamicKDRxnDestroy

! ************************************************************************** !

subroutine AqueousSpeciesConstraintDestroy(constraint)
  !
  ! Destroys an aqueous species constraint
  ! object
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  !

  use Utility_module, only: DeallocateArray

  implicit none

  type(aq_species_constraint_type), pointer :: constraint

  if (.not.associated(constraint)) return

  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_conc)
  call DeallocateArray(constraint%basis_molarity)
  call DeallocateArray(constraint%constraint_type)
  call DeallocateArray(constraint%constraint_spec_id)
  call DeallocateArray(constraint%constraint_aux_string)
  call DeallocateArray(constraint%external_dataset)

  deallocate(constraint)
  nullify(constraint)

end subroutine AqueousSpeciesConstraintDestroy

! ************************************************************************** !

subroutine GuessConstraintDestroy(constraint)
  !
  ! Destroys an aqueous species constraint
  ! object
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  !

  use Utility_module, only: DeallocateArray

  implicit none

  type(guess_constraint_type), pointer :: constraint

  if (.not.associated(constraint)) return

  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%conc)

  deallocate(constraint)
  nullify(constraint)

end subroutine GuessConstraintDestroy

! ************************************************************************** !

subroutine ReactionDestroy(reaction,option)
  !
  ! Deallocates a reaction object
  !
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  !

  use Utility_module, only: DeallocateArray
  use Option_module

  implicit none

  class(reaction_rt_type), pointer :: reaction

  type(ion_exchange_rxn_type), pointer :: ionxrxn, prev_ionxrxn
  type(general_rxn_type), pointer :: general_rxn, prev_general_rxn
  type(radioactive_decay_rxn_type), pointer :: radioactive_decay_rxn, &
                                               prev_radioactive_decay_rxn
  type(dynamic_kd_rxn_type), pointer :: dynamic_kd_rxn, prev_dynamic_kd_rxn
  type(option_type) :: option

  if (.not.associated(reaction)) return

  call ReactionBaseStrip(reaction)

  !species index
  call SpeciesIndexDestroy(reaction%species_idx)

  ! primary species
  if (associated(reaction%primary_species_list)) &
    call AqueousSpeciesListDestroy(reaction%primary_species_list)
  nullify(reaction%primary_species_list)

  ! secondary species
  if (associated(reaction%secondary_species_list)) &
    call AqueousSpeciesListDestroy(reaction%secondary_species_list)
  nullify(reaction%secondary_species_list)

  ! ionx exchange reactions
  ionxrxn => reaction%ion_exchange_rxn_list
  do
    if (.not.associated(ionxrxn)) exit
    prev_ionxrxn => ionxrxn
    ionxrxn => ionxrxn%next
    call IonExchangeRxnDestroy(prev_ionxrxn)
  enddo
  nullify(reaction%ion_exchange_rxn_list)

  ! radioactive decay reactions
  radioactive_decay_rxn => reaction%radioactive_decay_rxn_list
  do
    if (.not.associated(radioactive_decay_rxn)) exit
    prev_radioactive_decay_rxn => radioactive_decay_rxn
    radioactive_decay_rxn => radioactive_decay_rxn%next
    call RadioactiveDecayRxnDestroy(prev_radioactive_decay_rxn)
  enddo
  nullify(reaction%radioactive_decay_rxn_list)

  ! general reactions
  general_rxn => reaction%general_rxn_list
  do
    if (.not.associated(general_rxn)) exit
    prev_general_rxn => general_rxn
    general_rxn => general_rxn%next
    call GeneralRxnDestroy(prev_general_rxn)
  enddo
  nullify(reaction%general_rxn_list)

  ! dynamic kd reactions
  dynamic_kd_rxn => reaction%dynamic_kd_rxn_list
  do
    if (.not.associated(dynamic_kd_rxn)) exit
    prev_dynamic_kd_rxn => dynamic_kd_rxn
    dynamic_kd_rxn => dynamic_kd_rxn%next
    call DynamicKDRxnDestroy(prev_dynamic_kd_rxn)
  enddo
  nullify(reaction%dynamic_kd_rxn_list)

  call SurfaceComplexationDestroy(reaction%surface_complexation)
  call MineralDestroy(reaction%mineral)
  call MicrobialDestroy(reaction%microbial)
  call ReactionImDestroy(reaction%immobile)
  call GasDestroy(reaction%gas)
  call IsothermDestroy(reaction%isotherm,option)
#ifdef SOLID_SOLUTION
  call SolidSolutionDestroy(reaction%solid_solution_list)
#endif

  if (associated(reaction%dbase_temperatures)) &
    deallocate(reaction%dbase_temperatures)
  nullify(reaction%dbase_temperatures)

  ! redox species
  if (associated(reaction%redox_species_list)) &
    call AqueousSpeciesListDestroy(reaction%redox_species_list)
  nullify(reaction%redox_species_list)

  call GenericParameterDestroy(reaction%aq_diffusion_coefficients)
  call GenericParameterDestroy(reaction%gas_diffusion_coefficients)

  call DeallocateArray(reaction%primary_species_names)
  call DeallocateArray(reaction%secondary_species_names)
  call DeallocateArray(reaction%eqcplx_basis_names)

  call DeallocateArray(reaction%primary_species_print)
  call DeallocateArray(reaction%secondary_species_print)
  call DeallocateArray(reaction%eqcplx_basis_print)
  call DeallocateArray(reaction%kd_print)
  call DeallocateArray(reaction%total_sorb_print)

  call DeallocateArray(reaction%primary_spec_a0)
  call DeallocateArray(reaction%primary_spec_Z)
  call DeallocateArray(reaction%primary_spec_molar_wt)

  call DeallocateArray(reaction%eqcplxspecid)
  call DeallocateArray(reaction%eqcplxstoich)
  call DeallocateArray(reaction%eqcplxh2oid)
  call DeallocateArray(reaction%eqcplxh2ostoich)
  call DeallocateArray(reaction%eqcplx_a0)
  call DeallocateArray(reaction%eqcplx_Z)
  call DeallocateArray(reaction%eqcplx_molar_wt)
  call DeallocateArray(reaction%eqcplx_logK)
  call DeallocateArray(reaction%eqcplx_logKcoef)

  call DeallocateArray(reaction%eqionx_rxn_Z_flag)
  call DeallocateArray(reaction%eqionx_rxn_cation_X_offset)
  call DeallocateArray(reaction%eqionx_rxn_to_surf)
  call DeallocateArray(reaction%eqionx_rxn_CEC)
  call DeallocateArray(reaction%eqionx_rxn_k)
  call DeallocateArray(reaction%eqionx_rxn_cationid)

#if 0
  call DeallocateArray(reaction%kinionx_CEC)
  call DeallocateArray(reaction%kinionx_k)
  call DeallocateArray(reaction%kinionx_cationid)
#endif

  call DeallocateArray(reaction%radiodecayspecid)
  call DeallocateArray(reaction%radiodecaystoich)
  call DeallocateArray(reaction%radiodecayforwardspecid)
  call DeallocateArray(reaction%radiodecay_kf)

  call DeallocateArray(reaction%generalspecid)
  call DeallocateArray(reaction%generalstoich)
  call DeallocateArray(reaction%generalforwardspecid)
  call DeallocateArray(reaction%generalforwardstoich)
  call DeallocateArray(reaction%generalbackwardspecid)
  call DeallocateArray(reaction%generalbackwardstoich)
  call DeallocateArray(reaction%generalh2oid)
  call DeallocateArray(reaction%generalh2ostoich)
  call DeallocateArray(reaction%general_kf)
  call DeallocateArray(reaction%general_kr)

  call DeallocateArray(reaction%eqdynamickdspecid)
  call DeallocateArray(reaction%eqdynamickdrefspecid)
  call DeallocateArray(reaction%eqdynamickdrefspechigh)
  call DeallocateArray(reaction%eqdynamickdlow)
  call DeallocateArray(reaction%eqdynamickdhigh)
  call DeallocateArray(reaction%eqdynamickdpower)

  deallocate(reaction)
  nullify(reaction)

end subroutine ReactionDestroy

end module Reaction_Aux_module
