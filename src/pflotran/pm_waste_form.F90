module PM_Waste_Form_class

! MODULE DESCRIPTION:
! ===========================================================================
! This process model calculates the radionuclide source term due to
! nuclear waste form dissolution. The process model includes a waste
! package degradation model (including an instant release fraction), and
! several waste form dissolution models. Radionuclide decay is calculated
! to adjust radionuclide mass fractions (concentrations) within the waste
! form before and after waste package breach.
! ===========================================================================

#include "petsc/finclude/petscvec.h"
  use petscvec

  use PM_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Geometry_module
  use Data_Mediator_Vec_class
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Region_module
  use Checkpoint_module
  use Kdtree_module

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  use Lookup_Table_module

  implicit none

  private

  PetscBool, public :: bypass_warning_message = PETSC_FALSE
  PetscBool, public :: FMDM_surrogate_knnr = PETSC_FALSE

! OBJECT rad_species_type:
! ========================
! ---------------------------------------------------------------------------
! Description:  This object describes a radionuclide (RN) inside a waste
! form. A linked list of these objects is a member of the base waste form
! mechanism object.
! ---------------------------------------------------------------------------
! formula_weight: [g-RN/mol] molar mass of the radionuclide (RN)
! decay_constant: [1/sec] decay rate constant of the radionuclide
! mass_fraction: [g-RN/g-bulk] mass fraction of radionuclide (RN) as defined
!    by the mass of RN over the mass of the bulk waste form
! inst_release_fraction: [-] the fraction of the radionuclide mass that is
!    instantly released from the waste package upon breach
! daugh_id: [-] daughter radionuclide id number
! daughter: daughter name string
! ispecies: [-] primary species id number of radionuclide in reactive
!    transport process model
! name: name string for radionuclide
! -----------------------------------------
  type, public :: rad_species_type
   PetscReal :: formula_weight
   PetscReal :: decay_constant
   PetscReal :: mass_fraction
   PetscReal :: inst_release_fraction
   PetscInt :: daugh_id
   character(len=MAXWORDLENGTH) :: daughter
   PetscInt :: ispecies
   character(len=MAXWORDLENGTH) :: name
  end type rad_species_type
! -----------------------------------------

! OBJECT wf_mechanism_base_type:
! ==============================
! ---------------------------------------------------------------------------
! Description:  This object describes/defines the behavior of the waste form
! object, and contains a linked list of radionuclides, parameters defining
! properties of the waste form object important to its dissolution behavior,
! and waste package degradation model parameters. This is the base mechanism
! object, and therefore the member procedure "Dissolution" must be extended.
! ---------------------------------------------------------------------------
! rad_species_list(:): pointer to a linked list of radionuclide objects
! num_species: [-] number of radionuclides in the waste form inventory
! seed: [-] integer that provides a seed for the random number
!    generator for normal distribution
! canister_degradation_model: Boolean which indicates if the waste package
!    degradation model is on or off
! vitality_rate_mean: [log10/yr] mean waste package degradation rate, used
!    to calculate a degradation rate from a normal distribution
! vitality_rate_stdev: [log10/yr] standard deviation of the waste package
!    degradation rate, used to calculate a degradation rate from a normal
!    distribution
! vitality_rate_trunc: [log10/yr] waste package degradation rate truncation
!    value, which is used to truncate the normal distribution
! canister_material_constant: [-] waste package material constant
! matrix_density: [kg/m3] waste form bulk matrix density
! specific_surface_area: [m2/kg] waste form surface area per waste form mass
! name: name string of the mechanism object
! next: pointer to next mechanism object in linked list
! Dissolution (procedure): must be extended; defines the dissolution
!    behavior of the waste form after breach occurs
! -----------------------------------------------------------
  type, public :: wf_mechanism_base_type
    type(rad_species_type), pointer :: rad_species_list(:)
    PetscInt :: num_species
    PetscInt :: seed
    PetscBool :: canister_degradation_model
    PetscReal :: vitality_rate_mean
    PetscReal :: vitality_rate_stdev
    PetscReal :: vitality_rate_trunc
    PetscReal :: canister_material_constant
    PetscReal :: matrix_density
    PetscReal :: specific_surface_area
    character(len=MAXWORDLENGTH) :: name
    class(wf_mechanism_base_type), pointer :: next
  contains
    procedure, public :: Dissolution => WFMechBaseDissolution
  end type wf_mechanism_base_type
! -----------------------------------------------------------

! OBJECT wf_mechanism_glass_type:
! ===============================
! ---------------------------------------------------------------------------
! Description:  Defines the dissolution behavior of a glass log type of
! waste form containing high level nuclear waste. This object extends the
! base mechanism object.
! ---------------------------------------------------------------------------
! dissolution_rate: [kg-glass/m2/sec] glass waste form dissolution rate
! k0: [kg-glass/m2/sec] base glass waste form dissolution rate
! k_long: [kg-glass/m2/sec] long-term glass waste form dissolution rate
! nu: [-] pH dependence parameter
! Ea: [J/mol] effective activation energy
! Q: [-] ion activity product of H4SiO4
! K: [-] equilibrium constant for rate limiting step, which is the activity
!    of H4SiO4 at saturation with glass
! v: [-] affinity term exponent
! pH: [-] the pH in the waste form region
! use_pH: Boolean that indicates if pH should be calculated from simulation
! use_Q: Boolean that indicates if Q should be calulated from simulation
! h_ion_id: [-] species id of H+ ion in reactive transport process model
! SiO2_id: [-] species id of SiO2 in reactive transport process model
! Dissolution (procedure): calculates the glass dissolution rate
! ------------------------------------------------------------------------
  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_glass_type
    PetscReal :: dissolution_rate
    PetscReal :: k0
    PetscReal :: k_long
    PetscReal :: nu
    PetscReal :: Ea
    PetscReal :: Q
    PetscReal :: K
    PetscReal :: v
    PetscReal :: pH
    PetscBool :: use_pH
    PetscBool :: use_Q
    PetscInt :: h_ion_id
    PetscInt :: SiO2_id
  contains
    procedure, public :: Dissolution => WFMechGlassDissolution
  end type wf_mechanism_glass_type
! ------------------------------------------------------------------------

! OBJECT wf_mechanism_dsnf_type:
! ==============================
! ---------------------------------------------------------------------------
! Description:  Defines the dissolution behavior of defense-related spent
! nuclear fuel type of waste form containing high level nuclear waste. This
! object extends the base mechanism object.
! ---------------------------------------------------------------------------
! frac_dissolution_rate: [1/sec] fractional dissolution rate of the waste
!    form
! Dissolution (procedure): calculates the DSNF dissolution rate
! -----------------------------------------------------------------------
  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_dsnf_type
    PetscReal :: frac_dissolution_rate    ! 1/sec
  contains
    procedure, public :: Dissolution => WFMechDSNFDissolution
  end type wf_mechanism_dsnf_type
! -----------------------------------------------------------------------

! OBJECT wf_mechanism_wipp_type:
! ==============================
! ---------------------------------------------------------------------------
! Description:  Defines the dissolution behavior of transuranic waste at the
! Waste Isolation Pilot Plant (WIPP). This object extends the DSNF mechanism
! object. When using the WIPP waste form mechanism, the UFD_DECAY process
! model must also be used.
! Note: when selecting for DSNF and WIPP together, class is() can be used,
! but when selecting for either DSNF or WIPP, type is() should be used
! ---------------------------------------------------------------------------
! All member variables and procedures are defined by the DSNF mechanism
!    object.
! -----------------------------------------------------------------------
  type, public, extends(wf_mechanism_dsnf_type) :: wf_mechanism_wipp_type
  end type wf_mechanism_wipp_type
! -----------------------------------------------------------------------

! OBJECT wf_mechanism_fmdm_type:
! ==============================
! ---------------------------------------------------------------------------
! Description:  Defines the dissolution behavior of uranium dioxide high
! level nuclear waste through coupling to an external model called the
! Fuel Matrix Degradation Model (FMDM). This object extends the base
! mechanism object.
! ---------------------------------------------------------------------------
! dissolution_rate: [kg-bulk/m2/sec] bulk dissolution rate of the waste form
! frac_dissolution_rate: [1/sec] fractional dissolution rate of the waste form
! burnup: [GWd/MTHM] waste form burnup if the FMDM is linked
! burnup: [kg-bulk/m2/sec] used as bulk dissolution rate of the waste form
!    if the FMDM is not linked
! num_grid_cells_in_waste_form: [-] number of grid cells in the 1D
!    calculations within the FMDM (currently hardwired to 40)
! mapping_fmdm(:): [-] mapping of fmdm species into fmdm concentration array
! mapping_fmdm_to_pflotran(:): [-] mapping of species in fmdm concentration
!    array to pflotran
! concentration(:,:): [mol/L] concentrations of chemical species relevant to
!    the calculation of dissolution rate in the FMDM, sized by
!    (num_concentrations,num_grid_cells_in_waste_form)
! num_concentrations: [-] number of chemical species (currently hardwired
!    to 11)
! i*: [-] species id number
! Dissolution (procedure): calculates the FMDM dissolution rate
! -----------------------------------------------------------------------
  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_fmdm_type
    PetscReal :: dissolution_rate
    PetscReal :: frac_dissolution_rate
    PetscReal :: burnup
    PetscInt :: num_grid_cells_in_waste_form
    PetscInt, pointer :: mapping_fmdm(:)
    PetscInt, pointer :: mapping_fmdm_to_pflotran(:)
    PetscReal, pointer :: concentration(:,:)
    PetscInt :: num_concentrations
    PetscInt :: iUO2_2p
    PetscInt :: iUCO3_2n
    PetscInt :: iUO2
    PetscInt :: iCO3_2n
    PetscInt :: iO2
    PetscInt :: iH2O2
    PetscInt :: iFe_2p
    PetscInt :: iH2
    PetscInt :: iUO2_sld
    PetscInt :: iUO3_sld
    PetscInt :: iUO4_sld
  contains
    procedure, public :: Dissolution => WFMechFMDMDissolution
  end type wf_mechanism_fmdm_type
! -----------------------------------------------------------------------

! OBJECT wf_mechanism_fmdm_surrogate_type:
! ========================================
! ---------------------------------------------------------------------------
! Description:  Defines the dissolution behavior of uranium dioxide high
! level nuclear waste through coupling to a single-layer feed-forward
! artificial neural network or k nearest neighbors SURROGATE APPROXIMATION
! of an external model called the Fuel Matrix Degradation Model (FMDM). This
! object extends the base mechanism object.
! ---------------------------------------------------------------------------
! dissolution_rate: [kg-bulk/m2/sec] bulk dissolution rate of the waste form
! frac_dissolution_rate: [1/sec] fractional dissolution rate of the waste form
! burnup: [GWd/MTHM] waste form burnup if the FMDM is linked
! mapping_fmdm(:): [-] mapping of fmdm species into fmdm concentration array
! mapping_fmdm_to_pflotran(:): [-] mapping of species in fmdm concentration
!    array to pflotran
! num_concentrations: [-] number of chemical species (currently hardwired
!    to 4 environment concentrations)
! i*: [-] species id number
! Dissolution (procedure): calculates the FMDM dissolution rate
! -----------------------------------------------------------------------
  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_fmdm_surrogate_type
    PetscReal :: dissolution_rate
    PetscReal :: frac_dissolution_rate
    PetscReal :: burnup
    PetscReal :: decay_time ! specific to this waste form
    PetscInt, pointer :: mapping_surrfmdm(:)
    PetscInt, pointer :: mapping_surrfmdm_to_pflotran(:)
    PetscReal, pointer :: concentration(:)
    PetscInt :: num_concentrations
    PetscInt :: iCO3_2n
    PetscInt :: iO2
    PetscInt :: iFe_2p
    PetscInt :: iH2
    ! ANN parameters
    PetscReal :: input_hidden1_weights(6,64)
    PetscReal :: input_hidden1_bias(64)
    PetscReal :: hidden1_hidden2_weights(64,64)
    PetscReal :: hidden1_hidden2_bias(64)
    PetscReal :: hidden2_output_weights(64)
    PetscReal :: hidden2_output_bias
    PetscReal :: scaler_offsets(6)
    PetscReal :: scaler_scales(6)
    ! kNNr variables
    PetscInt :: num_nearest_neighbor
    type(kdtree), pointer :: tree
    PetscReal, pointer :: knnr_array(:,:)
    PetscInt :: num_qoi
    PetscReal, pointer :: table_data(:,:)
    PetscReal :: knnr_eps
  contains
    procedure, public :: Dissolution => WFMechFMDMSurrogateDissolution
  end type wf_mechanism_fmdm_surrogate_type
! -----------------------------------------------------------------------

! OBJECT wf_mechanism_custom_type:
! ================================
! ---------------------------------------------------------------------------
! Description:  Defines the dissolution behavior of a custom type of waste
! form containing high level nuclear waste. This object extends the base
! mechanism object.
! ---------------------------------------------------------------------------
! dissolution_rate: [kg-bulk/m2/sec] bulk dissolution rate of the waste form
! frac_dissolution_rate: [1/sec] fractional dissolution rate of the waste form
! Dissolution (procedure): calculates the CUSTOM dissolution rate
! -------------------------------------------------------------------------
  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_custom_type
    PetscReal :: dissolution_rate
    PetscReal :: frac_dissolution_rate
    PetscBool :: frac_diss_vol_init
  contains
    procedure, public :: Dissolution => WFMechCustomDissolution
  end type wf_mechanism_custom_type
! -------------------------------------------------------------------------

! OBJECT waste_form_base_type:
! ============================
! ---------------------------------------------------------------------------
! Description:
! ---------------------------------------------------------------------------
! id: [-] waste form id number
! rank_list(:): [-] array of 1's and 0's used to determine local waste forms
! coordinate: pointer to coordinate object that stores waste form's location
!    if a coordinate point was given
! region_name: name string for waste form's region
! region: pointer to the waste form's region object
! scaling_factor(:): [-] array of volume scaling factor for each grid cell
!    in the waste form's region object
! init_volume: [m3] initial waste form volume
! volume: [m3] current waste form volume
! exposure_factor: [-] multiplying factor to the waste form dissolution rate,
!    by default the value is 1.d0
! eff_dissolution_rate: [kg-bulk/sec] effective waste form dissolution rate
!    which takes into account the specific surface area, matrix density,
!    volume, and exposure factor
! instantaneous_mass_rate(:): [mol/sec] radionuclide source term
! cumulative_mass(:): [mol] cumulative mass of radionuclide released
! rad_mass_fraction(:): [g-RN/g-bulk] current radionuclide (RN) mass fraction
!    in the waste form
! rad_concentration(:): [mol-RN/g-bulk] current radionuclide (RN)
!    concentration in the waste form
! inst_release_amount(:): [mol-RN/g-bulk] fraction of current radionuclide
!    (RN) concentration that is instantly released upon waste package breach
! canister_degradation_flag: Boolean that indicates if the waste package
!    degradation model is on of off
! canister_vitality: [%] current waste package vitality (between 0 and 1)
! canister_vitality_rate: [%/sec] base rate of vitality degradation
! eff_canister_vit_rate: [%/sec] effective rate of vitality degradation
!    after effects of temperature and canister material constant
! breach_time: [sec] time of waste package breach
! breached: Boolean indicating if waste package has breached
! decay_start_time: [sec] time in simuation when radionuclides in the waste
!    form should start decaying, default time is 0.d0 sec
! mech_name: name string for the waste form mechanism object
! mechanism: pointer to waste form's mechanism object
! crit: pointer to the criticality object
! next: pointer to next waste form object in linked list
! -----------------------------------------------------
  type :: waste_form_base_type
    PetscInt :: id
    PetscInt, pointer :: rank_list(:)
    type(point3d_type) :: coordinate
    character(len=MAXWORDLENGTH) :: region_name
    type(region_type), pointer :: region
    PetscReal, pointer :: scaling_factor(:)
    PetscReal :: init_volume
    PetscReal :: volume
    PetscReal :: exposure_factor
    PetscReal :: eff_dissolution_rate
    PetscReal, pointer :: instantaneous_mass_rate(:)
    PetscReal, pointer :: cumulative_mass(:)
    PetscReal, pointer :: rad_mass_fraction(:)
    PetscReal, pointer :: rad_concentration(:)
    PetscReal, pointer :: inst_release_amount(:)
    PetscBool :: canister_degradation_flag
    PetscBool :: spacer_degradation_flag
    PetscReal :: canister_vitality
    PetscReal :: canister_vitality_rate
    PetscReal :: eff_canister_vit_rate
    PetscReal :: spacer_vitality
    PetscReal :: spacer_vitality_rate
    PetscReal :: breach_time
    PetscBool :: breached
    PetscReal :: decay_start_time
    character(len=MAXWORDLENGTH) :: mech_name
    character(len=MAXWORDLENGTH) :: spacer_mech_name
    character(len=MAXWORDLENGTH) :: criticality_mech_name
    class(wf_mechanism_base_type), pointer :: mechanism
    class(spacer_mechanism_base_type), pointer :: spacer_mechanism
    class(crit_mechanism_base_type), pointer :: criticality_mechanism
    class(waste_form_base_type), pointer :: next

  end type waste_form_base_type
! -----------------------------------------------------

! OBJECT pm_waste_form_type:
! ==========================
! ---------------------------------------------------------------------------
! Description:  This is the waste form process model object. It has a list of
! waste forms, mechanisms, and a data mediator vector. Several procedures
! allow interfacing with the process model structure and extend the
! pm_base_type procedures. This is the highest level object in this module.
! ---------------------------------------------------------------------------
! realization: pointer to the realization object
! data_mediator: pointer to the data mediator array which stores the values
!    of the radionuclide source terms
! waste_form_list: pointer to the linked list of waste form objects
! mechanism_list: pointer to the linked list of mechanism objects
! print_mass_balance: Boolean that indicates if the *.wf file is generated
! implicit_solution: Boolean that indicates if the explicit solution method
!    is used for calculation of isotope decay and ingrowth
! -------------------------------------------------------------------
  type, public, extends(pm_base_type) :: pm_waste_form_type
    class(realization_subsurface_type), pointer :: realization
    class(data_mediator_vec_type), pointer :: data_mediator
    class(waste_form_base_type), pointer :: waste_form_list
    class(wf_mechanism_base_type), pointer :: mechanism_list
    class(spacer_mechanism_base_type), pointer :: spacer_mech_list
    type(criticality_mediator_type), pointer :: criticality_mediator
    PetscBool :: print_mass_balance
    PetscBool :: implicit_solution
    PetscLogDouble :: cumulative_time
  contains
    procedure, public :: SetRealization => PMWFSetRealization
    procedure, public :: Setup => PMWFSetup
    procedure, public :: ReadPMBlock => PMWFReadPMBlock
    procedure, public :: InitializeRun => PMWFInitializeRun
    procedure, public :: InitializeTimestep => PMWFInitializeTimestep
    procedure, public :: FinalizeTimestep => PMWFFinalizeTimestep
    procedure, public :: UpdateSolution => PMWFUpdateSolution
    procedure, public :: Solve => PMWFSolve
    procedure, public :: CheckpointHDF5 => PMWFCheckpointHDF5
    procedure, public :: CheckpointBinary => PMWFCheckpointBinary
    procedure, public :: RestartHDF5 => PMWFRestartHDF5
    procedure, public :: RestartBinary => PMWFRestartBinary
    procedure, public :: InputRecord => PMWFInputRecord
    procedure, public :: Destroy => PMWFDestroy
  end type pm_waste_form_type

! -----------------------------------------------------

  ! Stores information regarding the criticality event
  type, public :: criticality_event_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: steady_state
    PetscReal :: crit_start
    PetscReal :: crit_end
    PetscBool :: crit_flag
  end type criticality_event_type

! -------------------------------------------------------------------

  ! Stores variables for the criticality heat emission lookup table
  type, public :: crit_heat_type
    character(len=MAXSTRINGLENGTH) :: file_name
    PetscInt :: num_start_times
    PetscInt :: num_values_per_start_time
    PetscReal :: start_time_datamax
    PetscReal :: temp_datamax
    PetscReal :: power_datamax
    class(lookup_table_general_type), pointer :: lookup_table
    class(crit_heat_type), pointer :: next
  contains
    procedure, public :: Read => CritHeatRead
    procedure, public :: Evaluate => CritHeatEvaluate
  end type crit_heat_type

! -------------------------------------------------------------------

  ! Provide one lookup table per nuclide
  type, public :: crit_inventory_lookup_type
    class(lookup_table_general_type), pointer :: lookup
    type(crit_inventory_lookup_type), pointer :: next
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: use_log10_time
    PetscReal :: log10_time_zero_sub
  contains
    procedure, public :: Evaluate => CritInventoryEvaluate
  end type crit_inventory_lookup_type

  ! Stores variables for the criticality inventory lookup table
  type, public :: crit_inventory_type
    character(len=MAXSTRINGLENGTH) :: file_name
    PetscInt :: total_points
    PetscInt :: num_start_times
    PetscInt :: num_powers
    PetscInt :: num_real_times
    PetscInt :: num_species
    PetscReal :: start_time_datamax
    PetscReal :: power_datamax
    PetscReal :: real_time_datamax
    PetscBool :: switch_implicit
    PetscBool :: allow_implicit
    PetscBool :: allow_extrap
    PetscBool :: continue_lookup
    type(crit_inventory_lookup_type), pointer :: radionuclide_table
    class(crit_inventory_type), pointer :: next
  contains
    procedure, public :: Read => CritInventoryRead
  end type crit_inventory_type

! -------------------------------------------------------------------

! Stores variables relevant to criticality calculations
  type, public :: crit_mechanism_base_type
    character(len=MAXWORDLENGTH) :: mech_name
    character(len=MAXSTRINGLENGTH) :: rad_dataset_name
    character(len=MAXSTRINGLENGTH) :: heat_dataset_name
    PetscReal :: decay_heat
    PetscReal :: crit_heat
    PetscReal :: sw
    PetscReal :: rho_w
    PetscReal :: temperature
    PetscReal :: k_effective
    PetscInt :: heat_source_cond
    type(criticality_event_type), pointer :: crit_event
    class(dataset_ascii_type), pointer :: rad_dataset
    class(dataset_ascii_type), pointer :: heat_dataset
    class(crit_heat_type), pointer :: crit_heat_dataset
    class(crit_inventory_type), pointer :: inventory_dataset
    class(crit_mechanism_base_type), pointer :: next
  end type crit_mechanism_base_type

! -----------------------------------------------------

  type, public :: criticality_mediator_type
    class(crit_mechanism_base_type), pointer :: crit_mech_list
    class(data_mediator_vec_type), pointer :: data_mediator
    PetscInt :: total_num_cells
  end type criticality_mediator_type

! -----------------------------------------------------

  type, public :: spacer_mechanism_base_type
    character(len=MAXWORDLENGTH) :: mech_name
    PetscReal :: threshold_sat ! threshold saturation for asm. exposure to water
    PetscReal :: alteration_rate  ! saturation-based factor for altering rate
    PetscReal :: spacer_mass  ! total mass of grid spacers
    PetscReal :: spacer_surface_area ! total surface area of grid spacers
    PetscReal :: spacer_coeff  ! empirical coefficient of Arrhenius term
    PetscReal :: spacer_activation_energy  ! activation energy
    class(spacer_mechanism_base_type), pointer :: next
  contains
    procedure, public :: Degradation => SpacerMechBaseDegradation
  end type spacer_mechanism_base_type

! -----------------------------------------------------

  public :: PMWFCreate, &
            PMWFSetup, &
            PMWFMechanismGlassCreate, &
            PMWFMechanismDSNFCreate, &
            PMWFMechanismWIPPCreate, &
            PMWFMechanismCustomCreate, &
            PMWFMechanismFMDMCreate, &
            PMWFMechanismFMDMSurrogateCreate

contains

! ************************************************************************** !

subroutine PMWFMechanismInit(this)
  !
  ! Initializes the base waste form mechanism
  !
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): base mechanism object
! -------------------------------------
  class(wf_mechanism_base_type) :: this
! -------------------------------------

  nullify(this%next)
  nullify(this%rad_species_list)
  this%num_species = 0
  this%seed = 1
  this%matrix_density = UNINITIALIZED_DOUBLE
  this%specific_surface_area = UNINITIALIZED_DOUBLE
  this%name = ''
 !---- canister degradation model ----------------------
  this%canister_degradation_model = PETSC_FALSE
  this%vitality_rate_mean = UNINITIALIZED_DOUBLE
  this%vitality_rate_stdev = UNINITIALIZED_DOUBLE
  this%vitality_rate_trunc = UNINITIALIZED_DOUBLE
  this%canister_material_constant = UNINITIALIZED_DOUBLE
 !------------------------------------------------------

end subroutine PMWFMechanismInit

! ************************************************************************** !

function PMWFMechanismGlassCreate()
  !
  ! Creates the glass waste form mechanism
  !
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none

! LOCAL VARIABLES:
! ================
! PMWFMechanismGlassCreate (output): new GLASS mechanism object
! glass: new GLASS mechanism object with shorter name
! -------------------------------------------------------------------
  class(wf_mechanism_glass_type), pointer :: PMWFMechanismGlassCreate
  class(wf_mechanism_glass_type), pointer :: glass
! -------------------------------------------------------------------

  allocate(glass)
  call PMWFMechanismInit(glass)
  glass%dissolution_rate = 0.d0        ! [kg/m^2/sec]
  glass%k0 = UNINITIALIZED_DOUBLE      ! [kg/m^2/sec]
  glass%k_long = UNINITIALIZED_DOUBLE  ! [kg/m^2/sec]
  glass%nu = UNINITIALIZED_DOUBLE      ! [-]
  glass%Ea = UNINITIALIZED_DOUBLE      ! [J/mol]
  glass%Q = UNINITIALIZED_DOUBLE       ! [-]
  glass%K = UNINITIALIZED_DOUBLE       ! [-]
  glass%v = UNINITIALIZED_DOUBLE       ! [-]
  glass%pH = UNINITIALIZED_DOUBLE      ! [-]
  glass%use_pH = PETSC_FALSE
  glass%use_Q = PETSC_FALSE
  glass%h_ion_id = 0
  glass%SiO2_id = 0

  PMWFMechanismGlassCreate => glass

end function PMWFMechanismGlassCreate

! ************************************************************************** !

function PMWFMechanismDSNFCreate()
  !
  ! Creates the DSNF (Defense Spent Nuclear Fuel) waste form mechanism
  !
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none

! LOCAL VARIABLES:
! ================
! PMWFMechanismDSNFCreate (output): new DSNF mechanism object
! dsnf: new DSNF mechanism object with shorter name
! -----------------------------------------------------------------
  class(wf_mechanism_dsnf_type), pointer :: PMWFMechanismDSNFCreate
  class(wf_mechanism_dsnf_type), pointer :: dsnf
! -----------------------------------------------------------------

  allocate(dsnf)
  call PMWFMechanismInit(dsnf)
  dsnf%frac_dissolution_rate = UNINITIALIZED_DOUBLE  ! 1/sec

  PMWFMechanismDSNFCreate => dsnf

end function PMWFMechanismDSNFCreate

! ************************************************************************** !

function PMWFMechanismWIPPCreate()
  !
  ! Creates the WIPP (Waste Isolation Pilot Plant) waste form mechanism
  !
  ! Author: Jenn Frederick
  ! Date: 012/7/2016

  implicit none

! LOCAL VARIABLES:
! ================
! PMWFMechanismWIPPCreate (output): new WIPP mechanism object
! wipp: new WIPP mechanism object with shorter name
! -----------------------------------------------------------------
  class(wf_mechanism_wipp_type), pointer :: PMWFMechanismWIPPCreate
  class(wf_mechanism_wipp_type), pointer :: wipp
! -----------------------------------------------------------------

  allocate(wipp)
  call PMWFMechanismInit(wipp)
  wipp%frac_dissolution_rate = UNINITIALIZED_DOUBLE  ! 1/sec

  PMWFMechanismWIPPCreate => wipp

end function PMWFMechanismWIPPCreate

! ************************************************************************** !

function PMWFMechanismFMDMCreate()
  !
  ! Creates the FMDM waste form mechanism package
  !
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none

! LOCAL VARIABLES:
! ================
! PMWFMechanismFMDMCreate (output): new FMDM mechanism object
! fmdm: new FMDM mechanism object with shorter name
! -----------------------------------------------------------------
  class(wf_mechanism_fmdm_type), pointer :: PMWFMechanismFMDMCreate
  class(wf_mechanism_fmdm_type), pointer :: fmdm
! -----------------------------------------------------------------

  allocate(fmdm)
  call PMWFMechanismInit(fmdm)

  fmdm%dissolution_rate = UNINITIALIZED_DOUBLE       ! kg/m^2/sec
  fmdm%frac_dissolution_rate = UNINITIALIZED_DOUBLE  ! 1/day
  fmdm%burnup = UNINITIALIZED_DOUBLE                 ! GWd/MTHM or (kg/m^2/sec)

  fmdm%num_grid_cells_in_waste_form = 40  ! hardwired

  nullify(fmdm%concentration)
  fmdm%num_concentrations = 11  ! hardwired
  fmdm%iUO2_2p = 1
  fmdm%iUCO3_2n = 2
  fmdm%iUO2 = 3
  fmdm%iCO3_2n = 4
  fmdm%iO2 = 5
  fmdm%iH2O2 = 6
  fmdm%iFe_2p = 7
  fmdm%iH2 = 8
  fmdm%iUO2_sld = 9
  fmdm%iUO3_sld = 10
  fmdm%iUO4_sld = 11

  allocate(fmdm%mapping_fmdm_to_pflotran(fmdm%num_concentrations))
  fmdm%mapping_fmdm_to_pflotran = UNINITIALIZED_INTEGER

  ! concentration can be allocated here because we hardwired
  ! the num_grid_cells_in_waste_form value, but if it becomes
  ! user defined, then allocation must be delayed until PMWFSetup
  allocate(fmdm%concentration(fmdm%num_concentrations, &
                              fmdm%num_grid_cells_in_waste_form))
  fmdm%concentration = 1.d-13

  allocate(fmdm%mapping_fmdm(4))
  fmdm%mapping_fmdm = [fmdm%iO2,fmdm%iCO3_2n,fmdm%iH2,fmdm%iFe_2p]

  PMWFMechanismFMDMCreate => fmdm

end function PMWFMechanismFMDMCreate

! ************************************************************************** !

function PMWFMechanismFMDMSurrogateCreate(option)
  !
  ! Creates the FMDM surrogate waste form mechanism package
  !
  ! Author: Tom Seidl
  ! Date: 03/05/2019

  implicit none

  type(option_type) :: option

! LOCAL VARIABLES:
! ================
! PMWFMechanismFMDMSurrogateCreate (output): new FMDM surrogate
! mechanism object
! surrfmdm: new FMDM mechanism object with shorter name
! -----------------------------------------------------------------
  class(wf_mechanism_fmdm_surrogate_type), pointer :: PMWFMechanismFMDMSurrogateCreate
  class(wf_mechanism_fmdm_surrogate_type), pointer :: surrfmdm
! -----------------------------------------------------------------

  allocate(surrfmdm)
  call PMWFMechanismInit(surrfmdm)

  surrfmdm%dissolution_rate = UNINITIALIZED_DOUBLE       ! kg/m^2/sec
  surrfmdm%frac_dissolution_rate = UNINITIALIZED_DOUBLE  ! 1/day
  surrfmdm%burnup = UNINITIALIZED_DOUBLE                 ! GWd/MTHM
  surrfmdm%decay_time = UNINITIALIZED_DOUBLE    ! K
  ! TS - may want to change to C later

  surrfmdm%num_concentrations = 4  ! hardwired
  surrfmdm%iCO3_2n = 1
  surrfmdm%iO2 = 2
  surrfmdm%iFe_2p = 3
  surrfmdm%iH2 = 4

  allocate(surrfmdm%mapping_surrfmdm_to_pflotran(surrfmdm%num_concentrations))
  surrfmdm%mapping_surrfmdm_to_pflotran = UNINITIALIZED_INTEGER

  allocate(surrfmdm%concentration(surrfmdm%num_concentrations))
  surrfmdm%concentration = 1.d-13

  allocate(surrfmdm%mapping_surrfmdm(4))
  surrfmdm%mapping_surrfmdm = [surrfmdm%iO2,surrfmdm%iCO3_2n, &
                               surrfmdm%iH2,surrfmdm%iFe_2p]

  if (FMDM_surrogate_knnr) then
    surrfmdm%knnr_eps = tiny (0.0d0)
    call KnnrInit(surrfmdm,option)
  else
    call ANNReadH5File(surrfmdm,option)
  endif

  PMWFMechanismFMDMSurrogateCreate => surrfmdm

end function PMWFMechanismFMDMSurrogateCreate

! ************************************************************************** !

function PMWFMechanismCustomCreate()
  !
  ! Creates the 'custom' waste form mechanism package
  !
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none

! LOCAL VARIABLES:
! ================
! PMWFMechanismCustomCreate (output): new CUSTOM mechanism object
! custom: new CUSTOM mechanism object with shorter name
! ---------------------------------------------------------------------
  class(wf_mechanism_custom_type), pointer :: PMWFMechanismCustomCreate
  class(wf_mechanism_custom_type), pointer :: custom
! ---------------------------------------------------------------------

  allocate(custom)
  call PMWFMechanismInit(custom)
  custom%dissolution_rate = UNINITIALIZED_DOUBLE      ! kg/m^2/sec
  custom%frac_dissolution_rate = UNINITIALIZED_DOUBLE ! 1/sec
  custom%frac_diss_vol_init = PETSC_FALSE

  PMWFMechanismCustomCreate => custom

end function PMWFMechanismCustomCreate


! ************************************************************************** !

subroutine PMWFRadSpeciesInit(rad_species)
  !
  ! Initializes a radioactive species in the waste form mechanism package
  !
  ! Author: Jenn Frederick
  ! Date: 03/09/16

  implicit none

! ----------------------------------------------
  type(rad_species_type) :: rad_species
! ----------------------------------------------

  rad_species%name = ''
  rad_species%daughter = ''
  rad_species%daugh_id = UNINITIALIZED_INTEGER
  rad_species%formula_weight = UNINITIALIZED_DOUBLE
  rad_species%decay_constant = UNINITIALIZED_DOUBLE
  rad_species%mass_fraction = UNINITIALIZED_DOUBLE
  rad_species%inst_release_fraction = UNINITIALIZED_DOUBLE
  rad_species%ispecies = UNINITIALIZED_INTEGER

end subroutine PMWFRadSpeciesInit

! ************************************************************************** !

function PMWFSpacerMechCreate()
!
! Creates a spacer grid degradation model in the waste form
!
! Author: Alex Salazar III
! Date: 05/06/2021

  implicit none

! LOCAL VARIABLES:
! ================
! ----------------------------------------------
  class(spacer_mechanism_base_type), pointer :: PMWFSpacerMechCreate
  class(spacer_mechanism_base_type), pointer :: spc
! ----------------------------------------------

  allocate(spc)
  nullify(spc%next)

  spc%mech_name = ''
  spc%threshold_sat = 0.0d0
  spc%alteration_rate = UNINITIALIZED_DOUBLE
  spc%spacer_mass = UNINITIALIZED_DOUBLE
  spc%spacer_surface_area = UNINITIALIZED_DOUBLE
  spc%spacer_coeff = UNINITIALIZED_DOUBLE
  spc%spacer_activation_energy = UNINITIALIZED_DOUBLE

  PMWFSpacerMechCreate => spc

end function PMWFSpacerMechCreate

! ************************************************************************** !

function PMWFWasteFormCreate()
  !
  ! Creates a waste form and initializes all parameters
  !
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none

! LOCAL VARIABLES:
! ================
! PMWFWasteFormCreate (output): new waste form object
! wf: new waste form object with a shorter name
! ----------------------------------------------------------
  type(waste_form_base_type), pointer :: PMWFWasteFormCreate
  type(waste_form_base_type), pointer :: wf
! ----------------------------------------------------------

  allocate(wf)
  wf%id = UNINITIALIZED_INTEGER
  nullify(wf%rank_list)
  wf%coordinate%x = UNINITIALIZED_DOUBLE
  wf%coordinate%y = UNINITIALIZED_DOUBLE
  wf%coordinate%z = UNINITIALIZED_DOUBLE
  nullify(wf%region)
  wf%region_name = ''
  nullify(wf%scaling_factor)            ! [-]
  wf%init_volume = UNINITIALIZED_DOUBLE ! m3
  wf%volume = UNINITIALIZED_DOUBLE      ! m3
  wf%exposure_factor = 1.0d0            ! [-]
  wf%eff_dissolution_rate = UNINITIALIZED_DOUBLE ! kg/sec
  wf%mech_name = ''
  wf%spacer_mech_name = ''
  wf%criticality_mech_name = ''
  wf%decay_start_time = 0.d0          ! sec (default value)
  nullify(wf%instantaneous_mass_rate) ! mol-rad/sec
  nullify(wf%cumulative_mass)         ! mol-rad
  nullify(wf%rad_mass_fraction)       ! g-rad/g-matrix
  nullify(wf%rad_concentration)       ! mol-rad/g-matrix
  nullify(wf%inst_release_amount)     ! mol-rad/g-matrix
  nullify(wf%mechanism)
  nullify(wf%spacer_mechanism)
  nullify(wf%criticality_mechanism)
  nullify(wf%next)
 !------- canister degradation model -----------------
  wf%canister_degradation_flag = PETSC_FALSE
  wf%breached = PETSC_FALSE
  wf%breach_time = UNINITIALIZED_DOUBLE
  wf%canister_vitality = 0.d0
  wf%canister_vitality_rate = UNINITIALIZED_DOUBLE
  wf%eff_canister_vit_rate = UNINITIALIZED_DOUBLE
  !------- spacer degradation model ------------------
  wf%spacer_degradation_flag = PETSC_FALSE
  wf%spacer_vitality = 1.d0
  wf%spacer_vitality_rate = UNINITIALIZED_DOUBLE
 !----------------------------------------------------

 PMWFWasteFormCreate => wf

end function PMWFWasteFormCreate

! ************************************************************************** !

function PMWFCreate()
  !
  ! Creates and initializes the waste form process model
  !
  ! Author: Glenn Hammond
  ! Date: 01/15/15, 07/20/15
  ! Notes: Modified by Jenn Frederick 03/24/2016

  implicit none

! LOCAL VARIABLES:
! ================
! PMWFCreate (output): new waste form process model object
! ------------------------------------------------
  class(pm_waste_form_type), pointer :: PMWFCreate
! ------------------------------------------------

  allocate(PMWFCreate)
  call PMBaseInit(PMWFCreate)
  nullify(PMWFCreate%realization)
  nullify(PMWFCreate%data_mediator)
  nullify(PMWFCreate%waste_form_list)
  nullify(PMWFCreate%mechanism_list)
  nullify(PMWFCreate%spacer_mech_list)
  nullify(PMWFCreate%criticality_mediator)
  PMWFCreate%print_mass_balance = PETSC_FALSE
  PMWFCreate%implicit_solution = PETSC_FALSE
  PMWFCreate%cumulative_time = 0.d0
  PMWFCreate%name = 'waste form general'
  PMWFCreate%header = 'WASTE FORM (GENERAL)'

end function PMWFCreate

! ************************************************************************** !

subroutine PMWFReadPMBlock(this,input)
  !
  ! Reads input file parameters associated with the waste form process model
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Notes: Modified by Jenn Frederick, 03/24/2016

  use Input_Aux_module
  use Option_module
  use String_module
  use Region_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (intput/output): waste form process model object
! input (input/output): pointer to input object
! ----------------------------------
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
!   class(simulation_subsurface_type) :: simulation
! ----------------------------------

! LOCAL VARIABLES:
! ================
! cur_waste_form: pointer to current waste form object
! cur_mechanism: pointer to current mechanism object
! option: pointer to option object
! word: temporary string
! error_string: error message string
! found: Boolean helper
! -------------------------------------------------------
  class(waste_form_base_type), pointer :: cur_waste_form
  class(wf_mechanism_base_type), pointer :: cur_mechanism
  class(spacer_mechanism_base_type), pointer :: cur_sp_mech
  class(crit_mechanism_base_type), pointer :: cur_crit_mech
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: assigned
  PetscInt  :: id = INPUT_RECORD_UNIT
  PetscBool :: is_open
! -------------------------------------------------------

  option => this%option
  input%ierr = 0
  error_string = 'WASTE_FORM_GENERAL'

  option%io_buffer = 'pflotran card:: ' // trim(error_string)
  call PrintMsg(option)

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    found = PETSC_FALSE

    select case(trim(word))
    !-------------------------------------
      case('PRINT_MASS_BALANCE')
        this%print_mass_balance = PETSC_TRUE
        cycle
    !-------------------------------------
      case('IMPLICIT_SOLUTION')
        this%implicit_solution = PETSC_TRUE
        cycle
    end select

    error_string = 'WASTE_FORM_GENERAL'
    call PMWFReadMechanism(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WASTE_FORM_GENERAL'
    call PMWFReadWasteForm(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WASTE_FORM_GENERAL'
    call PMWFReadSpacerMech(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WASTE_FORM_GENERAL'
    call ReadCriticalityMech(this,input,option,word,error_string,found)
    if (found) cycle

    if (.not. found) then
      option%io_buffer = 'Keyword "' // trim(word) // &
                         '" not applicable for the waste form process model.'
      call PrintErrMsg(option)
    endif

  enddo
  call InputPopBlock(input,option)

  ! Assign chosen mechanism to each waste form
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    cur_mechanism => this%mechanism_list
    do
      if (.not.associated(cur_mechanism)) exit
      if (StringCompare(cur_waste_form%mech_name,cur_mechanism%name)) then
        cur_waste_form%mechanism => cur_mechanism
        exit
      endif
      cur_mechanism => cur_mechanism%next
    enddo

    ! Assign chosen spacer grid degradation mechanism to each waste form object
    if (associated(this%spacer_mech_list)) then
      cur_sp_mech => this%spacer_mech_list
      do
        if (.not. associated(cur_sp_mech)) exit
        assigned = PETSC_FALSE
        if (StringCompare(cur_waste_form%spacer_mech_name, &
                          cur_sp_mech%mech_name)) then
          cur_waste_form%spacer_mechanism => cur_sp_mech
          assigned = PETSC_TRUE
          exit
        endif
        cur_sp_mech => cur_sp_mech%next
      enddo
      if (.not. assigned .and. &
          len_trim(cur_waste_form%spacer_mech_name) > 1) then
        option%io_buffer = 'Spacer degradation mechanism "' &
        // trim(cur_waste_form%spacer_mech_name) &
        //'" specified for waste form not found among the ' &
        //'available options.'
        call PrintErrMsg(option)
      endif
    endif

    ! Assign chosen criticality mechanism to each waste form object
    if (associated(this%criticality_mediator)) then
      cur_crit_mech => this%criticality_mediator%crit_mech_list
      do
        if (.not. associated(cur_crit_mech)) exit
        assigned = PETSC_FALSE
        if (StringCompare(cur_waste_form%criticality_mech_name, &
                          cur_crit_mech%mech_name)) then
          cur_waste_form%criticality_mechanism => cur_crit_mech
          assigned = PETSC_TRUE
          exit
        endif
        cur_crit_mech => cur_crit_mech%next
      enddo
      if (.not. assigned .and. &
          len_trim(cur_waste_form%criticality_mech_name) > 1) then
        option%io_buffer = 'Criticality mechanism "' &
        // trim(cur_waste_form%criticality_mech_name) &
        //'" specified for waste form not found among the ' &
        //'available options.'
        call PrintErrMsg(option)
      endif
    endif

    ! error messaging: ----------------------------------------------
    if (.not.associated(cur_waste_form%mechanism)) then
      option%io_buffer = 'WASTE_FORM MECHANISM ' // &
                         trim(cur_waste_form%mech_name) // &
                         ' not found among given mechanism names.'
      call PrintErrMsg(option)
    endif

    if (.not.cur_waste_form%mechanism%canister_degradation_model) then
      ! canister vitality specified, but can.deg. model is off:
      if (initialized(cur_waste_form%canister_vitality_rate)) then
        option%io_buffer = 'WASTE_FORM MECHANISM ' // &
          trim(cur_waste_form%mech_name) // ' does not have the canister &
          &degradation model turned on, but at least one of the waste forms &
          &assigned to this mechanism specifies a canister vitality rate.'
        call PrintErrMsg(option)
      endif
      ! canister breach time specified, but can.deg. model is off:
      if (initialized(cur_waste_form%breach_time)) then
        option%io_buffer = 'WASTE_FORM MECHANISM ' // &
          trim(cur_waste_form%mech_name) // ' does not have the canister &
          &degradation model turned on, but at least one of the waste forms &
          &assigned to this mechanism specifies a canister breach time.'
        call PrintErrMsg(option)
      endif
    endif

    ! both waste form and mechanism canister vitality rate parameters
    ! are specified:
    if (initialized(cur_waste_form%canister_vitality_rate) .and. &
        ( initialized(cur_waste_form%mechanism%vitality_rate_mean) .or. &
          initialized(cur_waste_form%mechanism%vitality_rate_stdev) .or. &
          initialized(cur_waste_form%mechanism%vitality_rate_trunc) )) then
      option%io_buffer = 'Either CANISTER_VITALITY_RATE within the &
        &WASTE_FORM blocks -or- the VITALITY_LOG10_MEAN, &
        &VITALITY_LOG10_STDEV, and VITALITY_UPPER_TRUNCATION within &
        &the WASTE_FORM MECHANISM ' // trim(cur_waste_form%mechanism%name) &
        // ' block should be specified, but not both.'
      call PrintErrMsg(option)
    endif

    ! the canister degradation model is on, but there are problems with
    ! the parameters provided:
    if (cur_waste_form%mechanism%canister_degradation_model) then
      ! all parameters are missing:
      if ( (Uninitialized(cur_waste_form%mechanism%vitality_rate_mean) .or. &
            Uninitialized(cur_waste_form%mechanism%vitality_rate_stdev) .or. &
            Uninitialized(cur_waste_form%mechanism%vitality_rate_trunc) ) .and. &
          Uninitialized(cur_waste_form%canister_vitality_rate) .and. &
          Uninitialized(cur_waste_form%breach_time)                 )  then
        option%io_buffer = 'CANISTER_VITALITY_RATE within the WASTE_FORM &
          &blocks -or- CANISTER_BREACH_TIME within the WASTE_FORM blocks &
          &-or- the VITALITY_LOG10_MEAN, VITALITY_LOG10_STDEV, and &
          &VITALITY_UPPER_TRUNCATION within the WASTE_FORM with MECHANISM ' // &
          trim(cur_waste_form%mechanism%name) // ' is missing.'
        call PrintErrMsg(option)
      endif
      ! all parameters are given:
      if ( (initialized(cur_waste_form%mechanism%vitality_rate_mean) .or. &
            initialized(cur_waste_form%mechanism%vitality_rate_stdev) .or. &
            initialized(cur_waste_form%mechanism%vitality_rate_trunc) ) .and. &
          initialized(cur_waste_form%canister_vitality_rate) .and. &
          initialized(cur_waste_form%breach_time)                 )  then
        option%io_buffer = 'CANISTER_VITALITY_RATE within the WASTE_FORM &
          &blocks -or- CANISTER_BREACH_TIME within the WASTE_FORM blocks &
          &-or- the VITALITY_LOG10_MEAN, VITALITY_LOG10_STDEV, and &
          &VITALITY_UPPER_TRUNCATION within the WASTE_FORM with MECHANISM ' // &
          trim(cur_waste_form%mechanism%name) // ' should be specified, &
          &but not all.'
        call PrintErrMsg(option)
      endif
      ! both breach time and can. deg. rate were given
      if (initialized(cur_waste_form%canister_vitality_rate) .and. &
          initialized(cur_waste_form%breach_time)) then
        option%io_buffer = 'Either CANISTER_VITALITY_RATE -or- &
          &CANISTER_BREACH_TIME within the WASTE_FORM block with &
          &WASTE_FORM MECHANISM ' // trim(cur_waste_form%mechanism%name) &
          // ' should be specified, but not both.'
        call PrintErrMsg(option)
      endif
      ! both breach time and can. deg. distribution were given
      if ((initialized(cur_waste_form%mechanism%vitality_rate_mean) .or. &
           initialized(cur_waste_form%mechanism%vitality_rate_stdev) .or. &
           initialized(cur_waste_form%mechanism%vitality_rate_trunc)) .and. &
          initialized(cur_waste_form%breach_time)) then
        option%io_buffer = 'Either CANISTER_BREACH_TIME within the &
          &WASTE_FORM block with WASTE_FORM MECHANISM ' &
          // trim(cur_waste_form%mechanism%name) // ' -or- the &
          &VITALITY_LOG10_MEAN, VITALITY_LOG10_STDEV, and &
          &VITALITY_UPPER_TRUNCATION within the WASTE_FORM with MECHANISM ' // &
          trim(cur_waste_form%mechanism%name) // ' should be specified, &
          &but not both.'
        call PrintErrMsg(option)
      endif
    endif

    if (cur_waste_form%spacer_degradation_flag .and. &
        .not. associated(this%spacer_mech_list)) then
      option%io_buffer = 'SPACER_MECHANISM_NAME "' &
                       // trim(cur_waste_form%spacer_mech_name) &
                       //'" was specified for waste form "' &
                       // trim(cur_waste_form%mechanism%name) &
                       //'" but no SPACER_DEGRADATION_MECHANISM block was ' &
                       //'actually provided.'
      call PrintErrMsg(option)
    endif

    cur_waste_form => cur_waste_form%next
  enddo

  inquire(id, OPENED=is_open)
  if (is_open .and. OptionPrintToFile(option)) then
    if (associated(this%waste_form_list)) then
      call WasteFormInputRecord(this%waste_form_list)
    endif
  endif

end subroutine PMWFReadPMBlock

! ************************************************************************** !

subroutine PMWFReadMechanism(this,input,option,keyword,error_string,found)
  !
  ! Reads input file parameters associated with the waste form mechanism
  !
  ! Author: Jenn Frederick
  ! Date: 03/24/2016
  !
  use Input_Aux_module
  use Reaction_Aux_module, only: GetPrimarySpeciesIDFromName
  use Option_module
  use Condition_module, only : ConditionReadValues
  use String_module
  use Units_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! input (input/output): pointer to input object
! option (input/output): pointer to option object
! keyword (input): keyword string
! error_string (input/output): error message string
! found (input/output): Boolean helper
! ----------------------------------------------
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
! ----------------------------------------------

! LOCAL VARIABLES:
! ================
! added: Boolean helper
! num_errors: [-] number of errors that occured during read
! word: temporary string
! temp_buf: temporary buffer string
! temp_species_array(:): temporary array of radionuclide species objects
! new_mechanism: pointer to new mechanism object
! cur_mechanism: pointer to current mechanism object
! k, j: [-] looping index integers
! double: temporary double precision number
! ----------------------------------------------------------------------
  PetscBool :: added
  PetscInt :: num_errors
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: temp_buf
  type(rad_species_type), allocatable :: temp_species_array(:)
  class(wf_mechanism_base_type), pointer :: new_mechanism, cur_mechanism
  PetscInt :: k, j
  PetscReal :: double
  PetscInt :: integer
  PetscErrorCode :: ierr
  PetscLogDouble :: log_start_time, log_end_time
! ----------------------------------------------------------------------

  error_string = trim(error_string) // ',MECHANISM'
  found = PETSC_TRUE
  added = PETSC_FALSE
  input%ierr = 0
  k = 0
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('MECHANISM')
      call InputReadCard(input,option,word)
      call InputErrorMsg(input,option,'mechanism type',error_string)
      num_errors = 0
      call StringToUpper(word)
      select case(trim(word))
      !---------------------------------
        case('GLASS')
          error_string = trim(error_string) // ' GLASS'
          allocate(new_mechanism)
          new_mechanism => PMWFMechanismGlassCreate()
      !---------------------------------
        case('DSNF')
          error_string = trim(error_string) // ' DSNF'
          allocate(new_mechanism)
          new_mechanism => PMWFMechanismDSNFCreate()
      !---------------------------------
        case('WIPP')
          error_string = trim(error_string) // ' WIPP'
          allocate(new_mechanism)
          new_mechanism => PMWFMechanismWIPPCreate()
      !---------------------------------
        case('FMDM')
          ! for now, set bypass_warning_message = TRUE so we can run
          ! the fmdm model even though its not included/linked
          bypass_warning_message = PETSC_TRUE
#ifndef FMDM_MODEL
          this%option%io_buffer = 'Preprocessing statement FMDM_MODEL must &
            &be defined and the ANL FMDM library must be linked to PFLOTRAN &
            &to employ the fuel matrix degradation model.'
          if (.not.bypass_warning_message) then
            call PrintErrMsg(this%option)
          endif
#endif
          error_string = trim(error_string) // ' FMDM'
          allocate(new_mechanism)
          new_mechanism => PMWFMechanismFMDMCreate()
      !---------------------------------
        case('FMDM_SURROGATE')
          error_string = trim(error_string) // ' FMDM_SURROGATE'
          allocate(new_mechanism)
          call PetscTime(log_start_time,ierr);CHKERRQ(ierr)
          new_mechanism => PMWFMechanismFMDMSurrogateCreate(option)
          call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
          this%cumulative_time = this%cumulative_time + (log_end_time - log_start_time)
      !---------------------------------
        case('FMDM_SURROGATE_KNNR')
          FMDM_surrogate_knnr = PETSC_TRUE
          error_string = trim(error_string) // ' FMDM_SURROGATE_KNNR'
          allocate(new_mechanism)
          call PetscTime(log_start_time,ierr);CHKERRQ(ierr)
          new_mechanism => PMWFMechanismFMDMSurrogateCreate(option)
          call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
          this%cumulative_time = this%cumulative_time + (log_end_time - log_start_time)
      !---------------------------------
        case('CUSTOM')
          error_string = trim(error_string) // ' CUSTOM'
          new_mechanism => PMWFMechanismCustomCreate()
      !---------------------------------
        case default
          option%io_buffer = 'Unrecognized mechanism type &
                             &in the ' // trim(error_string) // ' block.'
          call PrintErrMsg(option)
      !---------------------------------
      end select

      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !--------------------------
          case('NAME')
            call InputReadCard(input,option,word)
            call InputErrorMsg(input,option,'mechanism name',error_string)
            call StringToUpper(word)
            new_mechanism%name = trim(word)
        !--------------------------
          case('SPECIFIC_SURFACE_AREA')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'specific surface area', &
                               error_string)
            call InputReadAndConvertUnits(input,double,'m^2/kg', &
                            trim(error_string)//',specific surface area',option)
            select type(new_mechanism)
              class is(wf_mechanism_dsnf_type)
                ! applies to dsnf & wipp types
                option%io_buffer = 'ERROR: SPECIFIC_SURFACE_AREA cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
              class default
                new_mechanism%specific_surface_area = double
            end select
        !--------------------------
          case('MATRIX_DENSITY')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'matrix density',error_string)
            call InputReadAndConvertUnits(input,double,'kg/m^3', &
                                   trim(error_string)//',matrix density',option)
            select type(new_mechanism)
              class default
                new_mechanism%matrix_density = double
            end select
        !--------------------------
          case('SEED')
            call InputReadInt(input,option,integer)
            call InputErrorMsg(input,option,'seed',error_string)
            select type(new_mechanism)
              class default
                new_mechanism%seed = integer
            end select
        !--------------------------
          case('FRACTIONAL_DISSOLUTION_RATE_VI')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'fractional dissolution rate vi', &
                               error_string)
            call InputReadAndConvertUnits(input,double,'unitless/sec', &
                   trim(error_string)//',fractional dissolution rate vi',option)
            select type(new_mechanism)
              type is(wf_mechanism_custom_type)
                new_mechanism%frac_dissolution_rate = double
                new_mechanism%frac_diss_vol_init = PETSC_TRUE
              class default
                option%io_buffer = 'ERROR: FRACTIONAL_DISSOLUTION_RATE_VI &
                           &cannot be specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('FRACTIONAL_DISSOLUTION_RATE')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'fractional dissolution rate', &
                               error_string)
            call InputReadAndConvertUnits(input,double,'unitless/sec', &
                      trim(error_string)//',fractional dissolution rate',option)
            select type(new_mechanism)
              type is(wf_mechanism_custom_type)
                new_mechanism%frac_dissolution_rate = double
              class default
                option%io_buffer = 'ERROR: FRACTIONAL_DISSOLUTION_RATE cannot &
                                   &be specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('DISSOLUTION_RATE')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'dissolution rate',error_string)
            call InputReadAndConvertUnits(input,double,'kg/m^2-sec', &
                                 trim(error_string)//',dissolution_rate',option)
            select type(new_mechanism)
              type is(wf_mechanism_custom_type)
                new_mechanism%dissolution_rate = double
              class default
                option%io_buffer = 'ERROR: DISSOLUTION_RATE cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('K0')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'K0 (intrinsic dissolution rate)', &
                               error_string)
            call InputReadAndConvertUnits(input,double,'kg/m^2-sec', &
                  trim(error_string)//',K0 (intrinsic dissolution rate)',option)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%k0 = double
              class default
                option%io_buffer = 'ERROR: K0 (intrinsic dissolution rate) &
                                &cannot be specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('K_LONG')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'K_LONG (dissolution rate)', &
                               error_string)
            call InputReadAndConvertUnits(input,double,'kg/m^2-sec', &
                    trim(error_string)//',K_LONG (dissolution rate)',option)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%k_long = double
              class default
                option%io_buffer = 'ERROR: K_LONG (dissolution rate) cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('NU')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'NU (pH dependence parameter)', &
                               error_string)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%nu = double
              class default
                option%io_buffer = 'ERROR: NU (pH dependence parameter) cannot &
                                   &be specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('EA')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'EA (effective activation energy)',&
                               error_string)
            call InputReadAndConvertUnits(input,double,'J/mol', &
                 trim(error_string)//',EA (effective activation energy)',option)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%Ea = double
              class default
                option%io_buffer = 'ERROR: EA (effective activation energy) &
                                &cannot be specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('Q')
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                temp_buf = input%buf
                call InputReadDouble(input,option,double)
                if (InputError(input)) then
                  word = adjustl(trim(temp_buf))
                  call StringToUpper(word)
                  if (trim(word) == 'AS_CALCULATED') then
                    new_mechanism%use_Q = PETSC_TRUE
                  else
                    option%io_buffer = 'ERROR: Q value (ion activity product) &
                                     &was not provided, or Q instructions not &
                                     &understood for ' // trim(error_string)
                    call PrintMsg(option)
                    num_errors = num_errors + 1
                  endif
                endif
                if (new_mechanism%use_Q) then
                  new_mechanism%Q = 0.d0  ! initializes to value
                else
                  new_mechanism%Q = double
                endif
              class default
                option%io_buffer = 'ERROR: Q (ion activity product) cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('K')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'K (equilibrium constant)',&
                               error_string)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%K = double
              class default
                option%io_buffer = 'ERROR: K (equilibrium constant) cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('V')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'V (exponent parameter)',&
                               error_string)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%v = double
              class default
                option%io_buffer = 'ERROR: V (exponent parameter) cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('PH')
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                temp_buf = input%buf
                call InputReadDouble(input,option,double)
                if (InputError(input)) then
                  word = adjustl(trim(temp_buf))
                  call StringToUpper(word)
                  if (trim(word) == 'AS_CALCULATED') then
                    new_mechanism%use_pH = PETSC_TRUE
                  else
                    option%io_buffer = 'ERROR: PH value was not provided, or &
                                       &PH instructions not understood for ' &
                                       // trim(error_string) // '.'
                    call PrintMsg(option)
                    num_errors = num_errors + 1
                  endif
                endif
                if (new_mechanism%use_pH) then
                  new_mechanism%pH = 7.d0  ! initializes to neutral value
                else
                  new_mechanism%pH = double
                endif
              class default
                option%io_buffer = 'ERROR: PH cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('KIENZLER_DISSOLUTION')
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%k0 = 560.d0/(24.d0*3600.d0)  ! kg/m^2-sec
                new_mechanism%k_long = 0.d0
                new_mechanism%nu = 0.d0
                new_mechanism%Ea = 7397.d0*8.314d0
                new_mechanism%Q = 0.d0
                new_mechanism%K = 1.d0     ! This value doesn't matter since Q=0
                new_mechanism%v = 1.d0
                new_mechanism%pH = 0.d0
              class default
                option%io_buffer = 'ERROR: KIENZLER_DISSOLUTION cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('BURNUP')
            select type(new_mechanism)
              type is(wf_mechanism_fmdm_type)
                call InputReadDouble(input,option,new_mechanism%burnup)
                call InputErrorMsg(input,option,'burnup',error_string)
#ifndef FMDM_MODEL
                ! if fmdm model is not on, then burnup is dissolution rate
                call InputReadAndConvertUnits(input,new_mechanism%burnup, &
                                'kg/m^2-sec',trim(error_string)//',burnup', &
                                option)
                option%io_buffer = 'Warning: FMDM is not linked, but an &
                                   &FMDM mechanism was defined. BURNUP &
                                   &will be used for fuel dissolution rate.'
                call PrintMsg(option)
#else
                option%io_buffer = 'FMDM is linked.'
                call PrintMsg(option)
#endif
              type is(wf_mechanism_fmdm_surrogate_type)
                call InputReadDouble(input,option,new_mechanism%burnup)
                call InputErrorMsg(input,option,'burnup',error_string)
                option%io_buffer = 'FMDM surrogate selected.'
                call PrintMsg(option)
              class default
                option%io_buffer = 'ERROR: BURNUP cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('DECAY_TIME')
            select type(new_mechanism)
              type is(wf_mechanism_fmdm_surrogate_type)
                call InputReadDouble(input,option, &
                                     new_mechanism%decay_time)
                call InputErrorMsg(input,option,'DECAY_TIME', &
                                   error_string)
                call InputReadAndConvertUnits(input, &
                  new_mechanism%decay_time,'year','DECAY_TIME', &
                  option)
              class default
                option%io_buffer = 'ERROR: DECAY_TIME cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('NEAREST_NEIGHBOR')
            select type(new_mechanism)
              type is(wf_mechanism_fmdm_surrogate_type)
               call InputReadInt(input,option, &
                     new_mechanism%num_nearest_neighbor)
                call InputErrorMsg(input,option,'NEAREST_NEIGHBOR', &
                                   error_string)
              class default
                option%io_buffer = 'ERROR: NEAREST_NEIGHBOR cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
             end select
        !--------------------------
          case('KNNR_EPS')
            select type(new_mechanism)
              type is(wf_mechanism_fmdm_surrogate_type)
               call InputReadDouble(input,option, &
                     new_mechanism%knnr_eps)
                call InputErrorMsg(input,option,'KNNR_EPS', &
                                   error_string)
              class default
                option%io_buffer = 'ERROR: KNNR_EPS cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
             end select
        !--------------------------
          case('SPECIES')
            allocate(temp_species_array(50))
            do
              call InputReadPflotranString(input,option)
              if (InputCheckExit(input,option)) exit
              k = k + 1
              if (k > 50) then
                option%io_buffer = 'More than 50 radionuclide species are &
                  &provided in the ' // trim(error_string) // &
                  ', SPECIES block.'
                call PrintErrMsgToDev(option, &
                                       'if reducing to less than 50 is not &
                                       &an option.')
              endif
              call PMWFRadSpeciesInit(temp_species_array(k))
              ! read species name
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'SPECIES name',error_string)
              temp_species_array(k)%name = trim(word)
              ! read species formula weight [g-species/mol-species]
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES formula weight', &
                                 error_string)
              temp_species_array(k)%formula_weight = double
              ! read species decay constant [1/sec]
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES decay rate constant', &
                                 error_string)
              temp_species_array(k)%decay_constant = double
              ! read species initial mass fraction [g-species/g-bulk]
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES initial mass fraction', &
                                 error_string)
              temp_species_array(k)%mass_fraction = double
              ! read species instant release fraction [-]
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES instant release &
                                 &fraction',error_string)
              temp_species_array(k)%inst_release_fraction = double
              ! read species daughter
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                temp_species_array(k)%daughter = trim(word)
              else
                temp_species_array(k)%daughter = 'no_daughter'
              endif
              new_mechanism%num_species = k
            enddo
            if (k == 0) then
              option%io_buffer = 'ERROR: At least one radionuclide species &
                                 &must be provided in the ' // &
                                 trim(error_string) // ', SPECIES block.'
              call PrintMsg(option)
              num_errors = num_errors + 1
            endif
            allocate(new_mechanism%rad_species_list(k))
            new_mechanism%rad_species_list(1:k) = temp_species_array(1:k)
            deallocate(temp_species_array)
            k = 0
            do while (k < new_mechanism%num_species)
              k = k + 1
              if (trim(new_mechanism%rad_species_list(k)%daughter) == &
                  'no_daughter') then
                new_mechanism%rad_species_list(k)%daugh_id = 0
              else
                j = 0
                do while (j < new_mechanism%num_species)
                  j = j + 1
                  if (trim(new_mechanism%rad_species_list(k)%daughter) == &
                       trim(new_mechanism%rad_species_list(j)%name)) then
                    new_mechanism%rad_species_list(k)%daugh_id = j
                    exit
                  endif
                enddo
              endif
            enddo
        !--------------------------
          case('CANISTER_DEGRADATION_MODEL')
            new_mechanism%canister_degradation_model = PETSC_TRUE
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputCheckExit(input,option)) exit
              call InputReadCard(input,option,word)
              call StringToUpper(word)
              select case(trim(word))
              case('VITALITY_LOG10_MEAN')
                call InputReadDouble(input,option, &
                                     new_mechanism%vitality_rate_mean)
                call InputErrorMsg(input,option,'canister vitality log-10 &
                                   &mean value',error_string)
              case('VITALITY_LOG10_STDEV')
                call InputReadDouble(input,option, &
                                     new_mechanism%vitality_rate_stdev)
                call InputErrorMsg(input,option,'canister vitality log-10 &
                                   &st. dev. value',error_string)
              case('VITALITY_UPPER_TRUNCATION')
                call InputReadDouble(input,option, &
                                     new_mechanism%vitality_rate_trunc)
                call InputErrorMsg(input,option,'canister vitality log-10 &
                                   &upper truncation value',error_string)
              case('CANISTER_MATERIAL_CONSTANT')
                call InputReadDouble(input,option, &
                                     new_mechanism%canister_material_constant)
                call InputErrorMsg(input,option,'canister material constant', &
                                   error_string)
              case default
                option%io_buffer = 'Keyword ' // trim(word) // &
                                   ' not recognized in the ' // &
                                   trim(error_string) // &
                                   ' CANISTER_DEGRADATION_MODEL block.'
                call PrintErrMsg(option)
              end select
            enddo
            call InputPopBlock(input,option)
        !--------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !--------------------------
        end select
      enddo
      call InputPopBlock(input,option)

     !----------- error messaging ----------------------------------------------
      if (new_mechanism%name == '') then
        option%io_buffer = 'ERROR: NAME must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      select type(new_mechanism)
        type is(wf_mechanism_glass_type)
          if (Uninitialized(new_mechanism%specific_surface_area)) then
            option%io_buffer = 'ERROR: SPECIFIC_SURFACE_AREA must be specified &
                               &in ' // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if (Uninitialized(new_mechanism%k0)) then
            option%io_buffer = 'ERROR: K0 must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if (Uninitialized(new_mechanism%k_long)) then
            option%io_buffer = 'ERROR: K_LONG must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if (Uninitialized(new_mechanism%nu)) then
            option%io_buffer = 'ERROR: NU must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if (Uninitialized(new_mechanism%Ea)) then
            option%io_buffer = 'ERROR: EA must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if (Uninitialized(new_mechanism%Q)) then
            option%io_buffer = 'ERROR: Q must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if (Uninitialized(new_mechanism%K)) then
            option%io_buffer = 'ERROR: K must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if (Uninitialized(new_mechanism%pH)) then
            option%io_buffer = 'ERROR: PH must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if (Uninitialized(new_mechanism%v)) then
            option%io_buffer = 'ERROR: V must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
        type is(wf_mechanism_custom_type)
          if (Uninitialized(new_mechanism%specific_surface_area) .and. &
              Uninitialized(new_mechanism%dissolution_rate) .and. &
              Uninitialized(new_mechanism%frac_dissolution_rate)) then
            option%io_buffer = 'ERROR: FRACTIONAL_DISSOLUTION_RATE, &
                               &FRACTIONAL_DISSOLUTION_RATE_VI, or &
                               &DISSOLUTION_RATE with SPECIFIC_SURFACE_AREA &
                               &must be specified in ' // trim(error_string) &
                               // ' ' // trim(new_mechanism%name) // ' block.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if ( (Initialized(new_mechanism%frac_dissolution_rate) .and. &
                Initialized(new_mechanism%dissolution_rate)    ) .or. &
               (Uninitialized(new_mechanism%frac_dissolution_rate) .and. &
                Uninitialized(new_mechanism%dissolution_rate)    ) ) then
            option%io_buffer = 'ERROR: Either FRACTIONAL_DISSOLUTION_RATE(_VI) &
                               &or DISSOLUTION_RATE with SPECIFIC_SURFACE_AREA &
                               &must be specified in ' // trim(error_string) &
                               // ' ' // trim(new_mechanism%name) // ' block. &
                               &Both types of dissolution rates cannot be &
                               &specified.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if ( (Initialized(new_mechanism%specific_surface_area) .and. &
                Uninitialized(new_mechanism%dissolution_rate)  ) .or. &
               (Uninitialized(new_mechanism%specific_surface_area) .and. &
                initialized(new_mechanism%dissolution_rate)      ) ) then
            option%io_buffer = 'ERROR: FRACTIONAL_DISSOLUTION_RATE(_VI) or &
                               &DISSOLUTION_RATE with SPECIFIC_SURFACE_AREA &
                               &must be specified in ' // trim(error_string) &
                               // ' ' // trim(new_mechanism%name) // ' block.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
        type is(wf_mechanism_fmdm_type)
          if (Uninitialized(new_mechanism%burnup)) then
            option%io_buffer = 'ERROR: BURNUP must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
          if (Uninitialized(new_mechanism%specific_surface_area)) then
            option%io_buffer = 'ERROR: SPECIFIC_SURFACE_AREA must be specified &
                               &in ' // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
      end select
      if (Uninitialized(new_mechanism%matrix_density)) then
        option%io_buffer = 'ERROR: MATRIX_DENSITY must be specified in ' // &
                           trim(error_string) // ' ' // &
                           trim(new_mechanism%name) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif

      if (new_mechanism%canister_degradation_model .and. &
          Uninitialized(new_mechanism%canister_material_constant)) then
        option%io_buffer = 'ERROR: CANISTER_MATERIAL_CONSTANT must be given in &
                           &the ' // trim(error_string) // ' ' // &
                           trim(new_mechanism%name) // &
                           ', CANISTER_DEGRADATION_MODEL block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif

      if (.not.associated(new_mechanism%rad_species_list)) then
        option%io_buffer = 'ERROR: At least one SPECIES must be specified in &
                           &the ' // trim(error_string) // ' ' // &
                           trim(new_mechanism%name) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif

      if (.not.associated(this%mechanism_list)) then
        this%mechanism_list => new_mechanism
      else
        cur_mechanism => this%mechanism_list
        do
          if (.not.associated(cur_mechanism)) exit
          if (.not.associated(cur_mechanism%next)) then
            cur_mechanism%next => new_mechanism
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_mechanism => cur_mechanism%next
        enddo
      endif
      nullify(new_mechanism)
  !-------------------------------------
    case default !(MECHANISM keyword not found)
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the MECHANISM block. See above.'
    call PrintErrMsg(option)
  endif

end subroutine PMWFReadMechanism

! ************************************************************************** !

subroutine PMWFReadWasteForm(this,input,option,keyword,error_string,found)
  !
  ! Reads input file parameters associated with the waste form
  !
  ! Author: Jenn Frederick
  ! Date: 03/24/2016
  !
  use Input_Aux_module
  use Reaction_Aux_module, only: GetPrimarySpeciesIDFromName
  use Option_module
  use Condition_module, only : ConditionReadValues
  use Dataset_Ascii_class
  use String_module
  use Units_module
  use Region_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! input (input/output): pointer to input object
! option (input/output): pointer to option object
! keyword (input): keyword string
! error_string (input/output): error message string
! found (input/output): Boolean helper
! ----------------------------------------------
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
! ----------------------------------------------

! LOCAL VARIABLES:
! ================
! added: Boolean helper
! num_errors: [-] number of errors during read
! word: temporary string
! new_waste_form: pointer to new waste form object
! cur_waste_form: pointer to current waste form object
! ----------------------------------------------------------------------
  PetscBool :: added
  PetscInt :: num_errors
  character(len=MAXWORDLENGTH) :: word
  class(waste_form_base_type), pointer :: new_waste_form, cur_waste_form
! ----------------------------------------------------------------------

  error_string = trim(error_string) // ',WASTE_FORM'
  found = PETSC_TRUE
  added = PETSC_FALSE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('WASTE_FORM')
      new_waste_form => PMWFWasteFormCreate()
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !-----------------------------
          case('EXPOSURE_FACTOR')
            call InputReadDouble(input,option,new_waste_form%exposure_factor)
            call InputErrorMsg(input,option,'exposure factor',error_string)
        !-----------------------------
          case('VOLUME')
            call InputReadDouble(input,option,new_waste_form%volume)
            call InputErrorMsg(input,option,'volume',error_string)
            call InputReadAndConvertUnits(input,new_waste_form%volume, &
                                          'm^3',trim(error_string)//',volume', &
                                          option)
            new_waste_form%init_volume = new_waste_form%volume
        !-----------------------------
          case('COORDINATE')
            call GeometryReadCoordinate(input,option, &
                                        new_waste_form%coordinate,error_string)
        !-----------------------------
          case('REGION')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'region assignment',error_string)
            new_waste_form%region_name = trim(word)
        !-----------------------------
          case('MECHANISM_NAME')
            call InputReadCard(input,option,word)
            call InputErrorMsg(input,option,'mechanism assignment',error_string)
            call StringToUpper(word)
            new_waste_form%mech_name = trim(word)
        !-----------------------------
          case('CANISTER_VITALITY_RATE')
            call InputReadDouble(input,option, &
                                 new_waste_form%canister_vitality_rate)
            call InputErrorMsg(input,option,'canister vitality rate', &
                               error_string)
            call InputReadAndConvertUnits(input, &
                        new_waste_form%canister_vitality_rate,'unitless/sec', &
                        trim(error_string)//'canister vitality rate',option)
        !-----------------------------
          case('CANISTER_BREACH_TIME')
            call InputReadDouble(input,option, &
                                 new_waste_form%breach_time)
            call InputErrorMsg(input,option,'CANISTER_BREACH_TIME',error_string)
            call InputReadAndConvertUnits(input,new_waste_form%breach_time, &
                           'sec',trim(error_string)//',CANISTER_BREACH_TIME', &
                           option)
        !-----------------------------
          case('DECAY_START_TIME')
            call InputReadDouble(input,option, &
                                 new_waste_form%decay_start_time)
            call InputErrorMsg(input,option,'DECAY_START_TIME',error_string)
            call InputReadAndConvertUnits(input, &
                 new_waste_form%decay_start_time,'sec',trim(error_string)// &
                 ',DECAY_START_TIME',option)
        !-----------------------------
          case('CRITICALITY_MECHANISM_NAME')
            call InputReadCard(input,option,word)
            call InputErrorMsg(input,option,'criticality mechanism ' &
                             //'assignment',error_string)
            call StringToUpper(word)
            new_waste_form%criticality_mech_name = trim(word)
        !-----------------------------
          case('SPACER_MECHANISM_NAME')
            call InputReadCard(input,option,word)
            call InputErrorMsg(input,option,'spacer grid degradation ' &
                            //'mechanism assignment',error_string)
            call StringToUpper(word)
            new_waste_form%spacer_mech_name = trim(word)
            new_waste_form%spacer_degradation_flag = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select

      enddo
      call InputPopBlock(input,option)


     ! ----------------- error messaging -------------------------------------
      if (Uninitialized(new_waste_form%volume)) then
        option%io_buffer = 'ERROR: VOLUME must be specified for all &
                           &waste forms.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Uninitialized(new_waste_form%coordinate%z) .and. &
          (len_trim(new_waste_form%region_name) == 0)) then
        option%io_buffer = 'ERROR: Either COORDINATE or REGION must be &
                           &specified for all waste forms.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Initialized(new_waste_form%coordinate%z) .and. &
          (len_trim(new_waste_form%region_name) > 0)) then
        option%io_buffer = 'ERROR: Either COORDINATE or REGION must be &
                           &specified for all waste forms, but not both.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (new_waste_form%mech_name == '') then
        option%io_buffer = 'ERROR: MECHANISM_NAME must be specified for &
                           &all waste forms.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      !note: do not throw error if EXPOSURE_FACTOR isn't specified (default = 1)

      if (.not.associated(this%waste_form_list)) then
        this%waste_form_list => new_waste_form
      else
        cur_waste_form => this%waste_form_list
        do
          if (.not.associated(cur_waste_form)) exit
          if (.not.associated(cur_waste_form%next)) then
            cur_waste_form%next => new_waste_form
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_waste_form => cur_waste_form%next
        enddo
      endif
      nullify(new_waste_form)
  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (.not.associated(this%waste_form_list)) then
    option%io_buffer = 'ERROR: At least one WASTE_FORM must be specified &
                       &in the WASTE_FORM_GENERAL block.'
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif

  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the WASTE_FORM_GENERAL,WASTE_FORM block(s). See above.'
    call PrintErrMsg(option)
  endif

end subroutine PMWFReadWasteForm

! ************************************************************************** !

subroutine PMWFReadSpacerMech(this,input,option,keyword,error_string,found)
  !
  ! Reads input file parameters associated with the spacer grid
  !   degradation model
  !
  ! Author: Alex Salazar III
  ! Date: 05/04/2021
  !
  use Input_Aux_module
  use Option_module
  use Condition_module, only : ConditionReadValues
  use Dataset_Ascii_class
  use String_module
  use Units_module
  use Region_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! input (input/output): pointer to input object
! option (input/output): pointer to option object
! keyword (input): keyword string
! error_string (input/output): error message string
! found (input/output): Boolean helper
! ----------------------------------------------
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
! ----------------------------------------------

! LOCAL VARIABLES:
! ================
! added: Boolean helper
! num_errors: [-] number of errors during read
! word: temporary string
! new_waste_form: pointer to new waste form object
! cur_waste_form: pointer to current waste form object
! ----------------------------------------------------------------------
  PetscBool :: added
  PetscInt :: num_errors
  character(len=MAXWORDLENGTH) :: word
  class(spacer_mechanism_base_type), pointer :: new_sp_mech, cur_sp_mech
! ----------------------------------------------------------------------

  error_string = trim(error_string) // ',SPACER_DEGRADATION_MECHANISM'
  found = PETSC_TRUE
  added = PETSC_FALSE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('SPACER_DEGRADATION_MECHANISM')
      allocate(new_sp_mech)
      new_sp_mech => PMWFSpacerMechCreate()
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !-----------------------------
         case('NAME')
           call InputReadCard(input,option,word)
           call InputErrorMsg(input,option, &
                              'spacer grid degradation mechanism assignment', &
                              error_string)
           call StringToUpper(word)
           new_sp_mech%mech_name = trim(word)
        !-----------------------------
          case('EXPOSURE_LEVEL')
            call InputReadDouble(input,option, &
                                 new_sp_mech%threshold_sat)
            call InputErrorMsg(input,option,'grid spacer exposure &
                               &saturation limit',error_string)
        !-----------------------------
          case('MASS')
            call InputReadDouble(input,option, &
                                 new_sp_mech%spacer_mass)
            call InputErrorMsg(input,option,'grid spacer total mass', &
                               error_string)
            call InputReadAndConvertUnits(input, new_sp_mech%spacer_mass, &
                                          'kg','MASS',option)
        !-----------------------------
          case('SURFACE_AREA')
            call InputReadDouble(input,option, &
                                 new_sp_mech%spacer_surface_area)
            call InputErrorMsg(input,option,'grid spacer total surface &
                               &area',error_string)
            call InputReadAndConvertUnits(input, &
                                          new_sp_mech%spacer_surface_area, &
                                          'm^2','SURFACE_AREA',option)
        !-----------------------------
          case('C')
            call InputReadDouble(input,option, &
                                 new_sp_mech%spacer_coeff)
            call InputErrorMsg(input,option,'grid spacer degradation model &
                               &empirical constant',error_string)
            call InputReadAndConvertUnits(input, &
                                          new_sp_mech%spacer_coeff, &
                                          'kg/m^2-s','C',option)
        !-----------------------------
          case('Q')
            call InputReadDouble(input,option, &
                                 new_sp_mech%spacer_activation_energy)
            call InputErrorMsg(input,option,'grid spacer degradation model &
                               activation energy',error_string)
            call InputReadAndConvertUnits(input, &
                                          new_sp_mech%spacer_activation_energy,&
                                          'J/mol','Q',option)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select

      enddo
      call InputPopBlock(input,option)


      ! --------------------------- error messaging ---------------------------
      if (len_trim(new_sp_mech%mech_name) < 1) then
        option%io_buffer = 'Name must be specified for spacer grid ' &
                         //'degradation mechanism in order to be associated ' &
                         //'with a waste form.'
        call PrintWrnMsg(option)
      endif
      if (Uninitialized(new_sp_mech%spacer_mass)) then
        option%io_buffer = 'ERROR: Total spacer grid mass must be specified &
                           &for spacer grid degradation mechanism.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Uninitialized(new_sp_mech%spacer_activation_energy)) then
        option%io_buffer = 'ERROR: Activation energy must be specified &
                           &for spacer grid degradation mechanism.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Uninitialized(new_sp_mech%spacer_coeff)) then
        option%io_buffer = 'ERROR: Scaling coefficient must be specified &
                           &for spacer grid degradation mechanism.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Uninitialized(new_sp_mech%spacer_surface_area)) then
        option%io_buffer = 'ERROR: Total spacer grid surface area must be &
                           &specified for spacer grid degradation mechanism.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (new_sp_mech%threshold_sat > 1.d0 .or. &
          new_sp_mech%threshold_sat < 0.d0) then
        option%io_buffer = 'ERROR: Saturation limit for exposure must be ' &
                         //'between 0 and 1 in the spacer grid degradation ' &
                         //'mechanism.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif

      if (.not.associated(this%spacer_mech_list)) then
        this%spacer_mech_list => new_sp_mech
      else
        cur_sp_mech => this%spacer_mech_list
        do
          if (.not.associated(cur_sp_mech)) exit
          if (.not.associated(cur_sp_mech%next)) then
            cur_sp_mech%next => new_sp_mech
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_sp_mech => cur_sp_mech%next
        enddo
      endif
      nullify(new_sp_mech)

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in ' &
                    //'the WASTE_FORM_GENERAL,SPACER_DEGRADATION_MECHANISM ' &
                    //'block(s). See above.'
    call PrintErrMsg(option)
  endif

end subroutine PMWFReadSpacerMech

! ************************************************************************** !

subroutine PMWFAssociateRegion(this,region_list)
  !
  ! Associates the waste form to its assigned region via the REGION keyword
  ! or the COORDINATE keyword.
  !
  ! Author: Jenn Frederick
  ! Date: 10/24/2016
  !

  use Region_module
  use Option_module
  use String_module
  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form object
! region_list (input): pointer to linked list of region objects
! ----------------------------------------------
  class(pm_waste_form_type) :: this
  type(region_list_type), pointer :: region_list
! ----------------------------------------------

! LOCAL VARIABLES:
! ================
! cur_region: pointer to current region object
! new_region: pointer to new region object
! cur_waste_form: pointer to current waste form object
! option: pointer to option object
! grid: pointer to grid object
! word*: temporary word string
! x, y, z: [m] coordinates
! i, j, k: indexing for structure grid
! local_id: [-] grid cell id number
! coordinate_counter: [-] counting integer helper for naming regions
! ------------------------------------------------------
  type(region_type), pointer :: cur_region
  type(region_type), pointer :: new_region
  class(waste_form_base_type), pointer :: cur_waste_form
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  character(len=MAXWORDLENGTH) :: word1, word2
  PetscReal :: x, y, z
  PetscInt :: i, j, k
  PetscInt :: local_id(1)
  PetscInt :: coordinate_counter
! ------------------------------------------------------

  option => this%option
  grid => this%realization%patch%grid
  coordinate_counter = 0

  cur_waste_form => this%waste_form_list

  do
    if (.not.associated(cur_waste_form)) exit
    ! if COORDINATE was given, auto-create a region for it
    if (Initialized(cur_waste_form%coordinate%z)) then
      local_id(1) = -1
      coordinate_counter = coordinate_counter + 1
      x = cur_waste_form%coordinate%x
      y = cur_waste_form%coordinate%y
      z = cur_waste_form%coordinate%z
      select case(grid%itype)
        case(STRUCTURED_GRID)
          call StructGridGetIJKFromCoordinate(grid%structured_grid,x,y,z, &
                                              i,j,k)
          if (i > 0 .and. j > 0 .and. k > 0) then
            local_id(1) = i + (j-1)*grid%structured_grid%nlx + &
                        (k-1)*grid%structured_grid%nlxy
          endif
        case(IMPLICIT_UNSTRUCTURED_GRID)
          call UGridGetCellFromPoint(x,y,z, &
                                     grid%unstructured_grid,option,local_id(1))
        case default
          option%io_buffer = 'Only STRUCTURED_GRID and &
                    &IMPLICIT_UNSTRUCTURED_GRID types supported in PMWasteForm.'
          call PrintErrMsg(option)
      end select
      ! create the region only if the current process owns the waste form
      if (local_id(1) > 0) then
        new_region => RegionCreate(local_id)
        write(word1,'(i6)') coordinate_counter
        write(word2,'(i6)') option%myrank
        new_region%name = 'WF_COORDINATE_' // trim(adjustl(word1)) // '_p' //  &
                          trim(adjustl(word2))
        cur_waste_form%region => new_region

        allocate(cur_waste_form%scaling_factor(1))
        cur_waste_form%scaling_factor(1) = 1.d0
      endif
    else
      cur_region => region_list%first
      do
        if (.not.associated(cur_region)) exit
        if (StringCompare(cur_region%name,cur_waste_form%region_name)) then
          cur_waste_form%region => cur_region

          exit
        endif
        cur_region => cur_region%next
      enddo
      if (.not.associated(cur_waste_form%region)) then
        option%io_buffer = 'WASTE_FORM REGION ' // &
                           trim(cur_waste_form%region_name) // ' not found.'
        call PrintErrMsg(option)
      endif
    endif
    !
    cur_waste_form => cur_waste_form%next
  enddo

end subroutine PMWFAssociateRegion

! ************************************************************************** !

subroutine PMWFSetRegionScaling(this,waste_form)
  !
  ! Calculates and sets the scaling factor vector for each of the waste forms
  ! that have assigned regions. This function is called only if a region was
  ! just associated with it. It assumes the volume of the cells that make up
  ! the region do not change over the course of the simulation.
  !
  ! Author: Jenn Frederick
  ! Date: 10/21/2016
  !

  use Material_Aux_module
  use Grid_module
  use Utility_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! waste_form (input/output): pointer to waste form object
! --------------------------------------------------
  class(pm_waste_form_type) :: this
  class(waste_form_base_type), pointer :: waste_form
! --------------------------------------------------

! LOCAL VARIABLES:
! ================
! material_auxvars(:): pointer to material auxvar object, which stores
!    the grid cell volume and indexed by the ghosted grid cell id
! grid: pointer to the grid object
! k: [-] looping index integer
! local_id: [-] local grid cell id
! ghosted_id: [-] ghosted grid cell id
! total_volume_local: [m3] total local waste form region volume
! total_volume_global: [m3] total global waste form region volume
! -----------------------------------------------------------
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type), pointer :: grid
  PetscInt :: k
  PetscInt :: local_id, ghosted_id
  PetscReal :: total_volume_local, total_volume_global
! -----------------------------------------------------------

  material_auxvars => this%realization%patch%aux%Material%auxvars
  grid => this%realization%patch%grid
  allocate(waste_form%scaling_factor(waste_form%region%num_cells))
  total_volume_local = 0.d0
  total_volume_global = 0.d0

  ! scale by cell volume
  do k = 1,waste_form%region%num_cells
    local_id = waste_form%region%cell_ids(k)
    ghosted_id = grid%nL2G(local_id)
    waste_form%scaling_factor(k) = material_auxvars(ghosted_id)%volume ! [m^3]
    total_volume_local = total_volume_local &
                         + material_auxvars(ghosted_id)%volume  ! [m^3]
  enddo
  call CalcParallelSUM(this%option,waste_form%rank_list, &
                       total_volume_local,total_volume_global)
  waste_form%scaling_factor = waste_form%scaling_factor/total_volume_global

end subroutine PMWFSetRegionScaling

! ************************************************************************** !

subroutine PMWFSetRealization(this,realization)
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Realization_Subsurface_class

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! realization (input): subsurface realization object
! ----------------------------------------------------------
  class(pm_waste_form_type) :: this
  class(realization_subsurface_type), pointer :: realization
! ----------------------------------------------------------

  this%realization => realization
  this%realization_base => realization

end subroutine PMWFSetRealization

! ************************************************************************** !

subroutine PMWFSetup(this)
  !
  ! Associates the waste forms to their regions and sets the waste form id.
  ! Creates an MPI group/communicator for processes that own a waste form.
  ! Throws out waste forms on processes that do not own the waste form region.
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  !
  ! Notes: Updated/modified by Jennifer Frederick, 10/24/2016.

  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Option_module
  use Reaction_Aux_module
  use NW_Transport_Aux_module
  use Utility_module, only : GetRndNumFromNormalDist
  use String_module

  implicit none

! INPUT ARGUEMENTS:
! =================
! this (input/output): waste form process model object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

! LOCAL VARIABLES:
! ================
! option: pointer to option object
! reaction: pointer to reaction object
! species_name: species name string
! names(:): array of name strings
! cur_waste_form: pointer to current waste form object
! prev_waste_form: pointer to previous waste form object
! next_waste_form: pointer to next waste form object
! cur_mechanism: pointer to current mechanism object
! waste_form_id: [-] waste form id number
! i, j,: [-] looping index integers
! local: Boolean helper indicating if waste form is local
! found: Boolean helper
! ierr: [-] PETSc error integer
! newcomm_size: [-] size of new MPI communicator number
! ranks(:): array of size(mycommsize) used to find local waste form objects
! -------------------------------------------------------
  type(option_type), pointer :: option
  class(reaction_rt_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: species_name
  character(len=MAXWORDLENGTH), pointer :: names(:)
  class(waste_form_base_type), pointer :: cur_waste_form
  class(waste_form_base_type), pointer :: prev_waste_form
  class(waste_form_base_type), pointer :: next_waste_form
  class(wf_mechanism_base_type), pointer :: cur_mechanism
  PetscInt :: waste_form_id, num_species
  PetscInt :: i, j
  PetscBool :: local, found
  PetscErrorCode :: ierr
  PetscMPIInt :: newcomm_size
  PetscInt, pointer :: ranks(:)
! -------------------------------------------------------

  option => this%realization%option
  reaction => this%realization%reaction

  ! point the waste form region to the desired region
  call PMWFAssociateRegion(this,this%realization%patch%region_list)

  allocate(ranks(option%comm%mycommsize))

  waste_form_id = 0
  nullify(prev_waste_form)
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    waste_form_id = waste_form_id + 1
    local = PETSC_FALSE
   !--------- canister degradation model --------------------
    if (cur_waste_form%mechanism%canister_degradation_model) then
      cur_waste_form%canister_degradation_flag = PETSC_TRUE
      cur_waste_form%canister_vitality = 1.d0
      ! waste form breach time specified:
      if (initialized(cur_waste_form%breach_time) .and. &
          Uninitialized(cur_waste_form%canister_vitality_rate)) then
        cur_waste_form%eff_canister_vit_rate = &
          (1.d0/cur_waste_form%breach_time)
      ! distribution for canister degradation rate specified:
      elseif (Uninitialized(cur_waste_form%canister_vitality_rate) .and. &
              Uninitialized(cur_waste_form%breach_time)) then
        ! call to random number generator must be done while each processor
        ! knows about every other processor's waste forms, otherwise the
        ! memory of the random number generator will not be global
        call GetRndNumFromNormalDist( &
             cur_waste_form%mechanism%vitality_rate_mean, &
             cur_waste_form%mechanism%vitality_rate_stdev, &
             cur_waste_form%canister_vitality_rate, &
             cur_waste_form%mechanism%seed)
        write(*,*) 'cur_waste_form%canister_vitality_rate = '
        write(*,*) cur_waste_form%canister_vitality_rate
        if (cur_waste_form%canister_vitality_rate > &
            cur_waste_form%mechanism%vitality_rate_trunc) then
          cur_waste_form%canister_vitality_rate = &
            cur_waste_form%mechanism%vitality_rate_trunc
        endif
        ! Given rates are in units of log-10/yr, so convert to 1/yr:
        cur_waste_form%canister_vitality_rate = &
          10.d0**(cur_waste_form%canister_vitality_rate)
        ! Convert rates from 1/yr to internal units of 1/sec
        cur_waste_form%canister_vitality_rate = &
          cur_waste_form%canister_vitality_rate * &
          (1.d0/DAYS_PER_YEAR/24.d0/3600.d0)
      endif
    endif
   !----------------------------------------------------------
    if (associated(cur_waste_form%region)) then
      if (cur_waste_form%region%num_cells > 0) then
          local = PETSC_TRUE
      endif
    endif
    ranks(:) = 0
    newcomm_size = 0
    if (local) then
      cur_waste_form%id = waste_form_id
      ranks(option%myrank+1) = 1
    else
      cur_waste_form%id = 0
      ranks(option%myrank+1) = 0
    endif
    call MPI_Allreduce(MPI_IN_PLACE,ranks,option%comm%mycommsize,MPI_INTEGER, &
                       MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
    newcomm_size = sum(ranks)
    allocate(cur_waste_form%rank_list(newcomm_size))
    j = 0
    do i = 1,option%comm%mycommsize
      if (ranks(i) == 1) then
        j = j + 1
        cur_waste_form%rank_list(j) = (i - 1)  ! (world ranks)
      endif
    enddo
    if (local) then
      call PMWFSetRegionScaling(this,cur_waste_form)
      prev_waste_form => cur_waste_form
      cur_waste_form => cur_waste_form%next
    else
      ! remove waste form because it is not local
      next_waste_form => cur_waste_form%next
      if (associated(prev_waste_form)) then
        prev_waste_form%next => next_waste_form
      else
        this%waste_form_list => next_waste_form
      endif
      call PMWFDestroyWasteForm(cur_waste_form)
      cur_waste_form => next_waste_form
    endif
  enddo

  deallocate(ranks)

  ! check if the mechanism list includes fmdm or glass mechanisms:
  cur_mechanism => this%mechanism_list
  do
    if (.not.associated(cur_mechanism)) exit
    select type(cur_mechanism)
      ! set up indexing of solute concentrations for fmdm model:
      type is(wf_mechanism_fmdm_type)
        species_name = 'O2(aq)'
        cur_mechanism%mapping_fmdm_to_pflotran(cur_mechanism%iO2) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
        species_name = 'CO3--'
        cur_mechanism%mapping_fmdm_to_pflotran(cur_mechanism%iCO3_2n) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
        species_name = 'H2(aq)'
        cur_mechanism%mapping_fmdm_to_pflotran(cur_mechanism%iH2) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
        species_name = 'Fe++'
        cur_mechanism%mapping_fmdm_to_pflotran(cur_mechanism%iFe_2p) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
      type is(wf_mechanism_fmdm_surrogate_type)
        species_name = 'O2(aq)'
        cur_mechanism%mapping_surrfmdm_to_pflotran(cur_mechanism%iO2) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
        species_name = 'CO3--'
        cur_mechanism%mapping_surrfmdm_to_pflotran(cur_mechanism%iCO3_2n) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
        species_name = 'H2(aq)'
        cur_mechanism%mapping_surrfmdm_to_pflotran(cur_mechanism%iH2) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
        species_name = 'Fe++'
        cur_mechanism%mapping_surrfmdm_to_pflotran(cur_mechanism%iFe_2p) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
      type is(wf_mechanism_glass_type)
        if (cur_mechanism%use_pH) then
          if ((associated(this%realization%reaction%species_idx) .and. &
              (this%realization%reaction%species_idx%h_ion_id == 0)) .or. &
              (.not.associated(this%realization%reaction%species_idx))) then
            option%io_buffer = 'pH may not be calculated when H+ is not &
                               &defined as a species and/or the card&
                               &USE_FULL_GEOCHEMISTRY is not specified in the &
                               &CHEMISTRY block - (MECHANISM GLASS).'
            call PrintErrMsg(option)
          else
            cur_mechanism%h_ion_id = &
                    this%realization%reaction%species_idx%h_ion_id
          endif
        endif
        if (cur_mechanism%use_Q) then
          species_name = 'SiO2(aq)'
          if (associated(this%realization%reaction)) then
            ! search through the species names so that the generic error
            ! message from GetPrimarySpeciesIDFromName is not thrown first
            ! when SiO2 is missing
            allocate(names(GetPrimarySpeciesCount(this%realization%reaction)))
            names => GetPrimarySpeciesNames(this%realization%reaction)
            i = 0
            found = PETSC_FALSE
            do while (i < len(names))
              i = i + 1
              if (adjustl(trim(species_name)) == adjustl(trim(names(i)))) then
                cur_mechanism%SiO2_id = &
                                     GetPrimarySpeciesIDFromName(species_name, &
                                     this%realization%reaction,option)
                found = PETSC_TRUE
              endif
              if (found) exit
              if ((.not.found) .and. (i == len(names))) then
                deallocate(names)
                allocate(names(GetSecondarySpeciesCount(this%realization%reaction)))
                names => GetSecondarySpeciesNames(this%realization%reaction)
                i = 0
                do while (i < len(names))
                  i = i + 1
                  if (adjustl(trim(species_name)) == &
                      adjustl(trim(names(i)))) then
                    cur_mechanism%SiO2_id = &
                                   GetSecondarySpeciesIDFromName(species_name, &
                                   this%realization%reaction,option)
                    found = PETSC_TRUE
                  endif
                  if (found) exit
                  if ((.not.found) .and. (i == len(names))) then
                    option%io_buffer = 'Q may not be calculated when SiO2(aq) &
                              &is not defined as a primary or secondary &
                              &species - (MECHANISM GLASS).'
                    call PrintErrMsg(option)
                  endif
                enddo
              endif
            enddo
            deallocate(names)
          else
            option%io_buffer = 'Q may not be calculated when SiO2(aq) is not &
                               &defined as a species and/or the card &
                               &USE_FULL_GEOCHEMISTRY is not specified in the &
                               &CHEMISTRY block - (MECHANISM GLASS).'
            call PrintErrMsg(option)
          endif
        endif
    end select
    cur_mechanism => cur_mechanism%next
  enddo

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    num_species = cur_waste_form%mechanism%num_species

    allocate(cur_waste_form%cumulative_mass(num_species))
    cur_waste_form%cumulative_mass = 0.d0
    allocate(cur_waste_form%rad_mass_fraction(num_species))
    cur_waste_form%rad_mass_fraction = &
      cur_waste_form%mechanism%rad_species_list%mass_fraction
    allocate(cur_waste_form%rad_concentration(num_species))
    cur_waste_form%rad_concentration = 0.d0
    allocate(cur_waste_form%instantaneous_mass_rate(num_species))
    cur_waste_form%instantaneous_mass_rate = 0.d0
    allocate(cur_waste_form%inst_release_amount(num_species))
    cur_waste_form%inst_release_amount = 0.d0
    do j = 1, num_species
      if (associated(this%realization%reaction_nw)) then
        cur_waste_form%mechanism%rad_species_list(j)%ispecies = &
          NWTGetSpeciesIDFromName( &
          cur_waste_form%mechanism%rad_species_list(j)%name, &
          this%realization%reaction_nw,this%option)
      elseif (associated(this%realization%reaction)) then
        cur_waste_form%mechanism%rad_species_list(j)%ispecies = &
          GetPrimarySpeciesIDFromName( &
          cur_waste_form%mechanism%rad_species_list(j)%name, &
          this%realization%reaction,this%option)
      else
        option%io_buffer = 'Currently, a transport process model (GIRT/OSRT ' &
                         //'or NWT) is required to use the WASTE_FORM_GENERAL '&
                         //'process model.'
        call PrintErrMsg(option)
      endif
    enddo
    cur_waste_form => cur_waste_form%next
  enddo


end subroutine PMWFSetup

! ************************************************************************** !

 subroutine PMWFInitializeRun(this)
  !
  ! Initializes the process model for the simulation
  !
  ! Author: Glenn Hammond
  ! Date: 08/25/15
  use Reaction_Aux_module
  use Reaction_Gas_Aux_module
  use Reaction_Sandbox_module
  use Realization_Base_class

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

! LOCAL VARIABLES:
! ================
! is: PETSc index set object
! cur_waste_form: pointer to current waste form object
! num_waste_form_cells: [-] number of cells in waste form region
! num_species: [-] number of readionuclide species objects
! size_of_vec: [-] size of data mediator array
! i, j, k: [-] looping index integers
! species_indices_in_residual(:): [-] array of indexes in the residual
! ierr: [-] PETSc error integer
! -------------------------------------------------------
  IS :: is
  class(waste_form_base_type), pointer :: cur_waste_form
  class(reaction_rt_type), pointer :: reaction
  PetscInt :: num_waste_form_cells
  PetscInt :: size_of_vec
  PetscInt :: i, j, k
  PetscInt, allocatable :: species_indices_in_residual(:), &
                           energy_indices_in_residual(:)
  PetscErrorCode :: ierr
! -------------------------------------------------------

  reaction => this%realization%reaction

  ! ensure that waste form is not being used with other reactive transport
  ! capability
  if (associated(reaction)) then
    if (reaction%neqcplx + reaction%nsorb + reaction%ngeneral_rxn + &
        reaction%microbial%nrxn + reaction%nradiodecay_rxn + &
        reaction%immobile%nimmobile > 0 .or. &
        GasGetCount(reaction%gas,ACTIVE_AND_PASSIVE_GAS) > 0 .or. &
        associated(rxn_sandbox_list)) then
      this%option%io_buffer = 'The UFD_DECAY process model may not be used &
        &with other reactive transport capability within PFLOTRAN. &
        &Minerals (and mineral kinetics) are used because their data &
        &structures are leveraged by UFD_DECAY, but no "conventional" &
        &mineral precipitation-dissolution capability is used.'
      call PrintErrMsg(this%option)
    endif
  endif

  if (this%print_mass_balance) then
    call PMWFOutputHeader(this)
    call PMWFOutput(this)
  endif

  !---------------- set up heat transfer (criticality) ---------------!

  if (associated(this%criticality_mediator)) then
    call RealizCreateFlowMassTransferVec(this%realization)
    this%criticality_mediator%data_mediator => DataMediatorVecCreate()
    call this%criticality_mediator%data_mediator%AddToList(this%realization% &
                                                        flow_data_mediator_list)
    size_of_vec = 0
    cur_waste_form => this%waste_form_list
    do
      if (.not.associated(cur_waste_form)) exit
      if (associated(cur_waste_form%criticality_mechanism)) then
        size_of_vec = size_of_vec + cur_waste_form%region%num_cells
      endif
      cur_waste_form => cur_waste_form%next
    enddo

    call VecCreateSeq(PETSC_COMM_SELF,size_of_vec, &
                      this%criticality_mediator%data_mediator%vec, &
                      ierr);CHKERRQ(ierr)
    call VecSetFromOptions(this%criticality_mediator%data_mediator%vec, &
                           ierr);CHKERRQ(ierr)

    cur_waste_form => this%waste_form_list
    allocate(energy_indices_in_residual(size_of_vec))
    j = 0
    do
      if (.not. associated(cur_waste_form)) exit
      if (associated(cur_waste_form%criticality_mechanism)) then
        do i = 1, cur_waste_form%region%num_cells
          j = j + 1
          energy_indices_in_residual(j) = cur_waste_form%region% &
                                          cell_ids(i) * this%option%nflowdof - 1
        enddo
      endif
      cur_waste_form => cur_waste_form%next
    enddo
    energy_indices_in_residual(:) = energy_indices_in_residual(:) + &
            this%realization%patch%grid%global_offset*this%option%nflowdof

    this%criticality_mediator%total_num_cells = j

    call ISCreateGeneral(this%option%mycomm,size_of_vec, &
                         energy_indices_in_residual,PETSC_COPY_VALUES,is, &
                         ierr);CHKERRQ(ierr)
    call VecScatterCreate(this%criticality_mediator%data_mediator%vec, &
                          PETSC_NULL_IS,this%realization%field%flow_r,is, &
                          this%criticality_mediator%data_mediator% &
                            scatter_ctx, &
                          ierr);CHKERRQ(ierr)
    if (allocated(energy_indices_in_residual)) then
        deallocate(energy_indices_in_residual)
    endif

    call ISDestroy(is,ierr);CHKERRQ(ierr)
  endif
  !------------- set up mass transfer ----------------------!

  call RealizCreateTranMassTransferVec(this%realization)
  this%data_mediator => DataMediatorVecCreate()
  call this%data_mediator%AddToList(this%realization%tran_data_mediator_list)
  ! create a Vec sized by # waste packages * # waste package cells in region *
  ! # primary dofs influenced by waste package
  ! count of waste form cells
  cur_waste_form => this%waste_form_list
  num_waste_form_cells = 0
  size_of_vec = 0
  do
    if (.not.associated(cur_waste_form)) exit
    size_of_vec = size_of_vec + (cur_waste_form%mechanism%num_species * &
                                 cur_waste_form%region%num_cells)
    num_waste_form_cells = num_waste_form_cells + 1
    cur_waste_form => cur_waste_form%next
  enddo
  call VecCreateSeq(PETSC_COMM_SELF,size_of_vec,this%data_mediator%vec, &
                    ierr);CHKERRQ(ierr)
  call VecSetFromOptions(this%data_mediator%vec,ierr);CHKERRQ(ierr)

  if (num_waste_form_cells > 0) then
    allocate(species_indices_in_residual(size_of_vec))
    species_indices_in_residual = 0
    cur_waste_form => this%waste_form_list
    i = 0
    do
      if (.not.associated(cur_waste_form)) exit
      do k = 1,cur_waste_form%region%num_cells
        do j = 1,cur_waste_form%mechanism%num_species
          i = i + 1
          species_indices_in_residual(i) = &
              (cur_waste_form%region%cell_ids(k)-1)*this%option%ntrandof + &
              cur_waste_form%mechanism%rad_species_list(j)%ispecies
        enddo
      enddo
      cur_waste_form => cur_waste_form%next
    enddo                             ! zero-based indexing
    species_indices_in_residual(:) = species_indices_in_residual(:) - 1
    ! set to global petsc index
    species_indices_in_residual(:) = species_indices_in_residual(:) + &
      this%realization%patch%grid%global_offset*this%option%ntrandof
  endif
  call ISCreateGeneral(this%option%mycomm,size_of_vec, &
                       species_indices_in_residual,PETSC_COPY_VALUES,is, &
                       ierr);CHKERRQ(ierr)
  if (allocated(species_indices_in_residual)) &
    deallocate(species_indices_in_residual)
  call VecScatterCreate(this%data_mediator%vec,PETSC_NULL_IS, &
                        this%realization%field%tran_r,is, &
                        this%data_mediator%scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)

  call PMWFSolve(this,0.d0,ierr)

end subroutine PMWFInitializeRun

! ************************************************************************** !

subroutine PMWFInitializeTimestep(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Notes: Modified by Jenn Frederick 03/28/2016

  use Utility_module
  use Global_Aux_module
  use Material_Aux_module
  use Field_module
  use Option_module
  use Grid_module
  use Patch_module
  use Utility_module
  use Dataset_Ascii_class
  use String_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

! LOCAL VARIABLES:
! ================
! cur_waste_form: pointer to current waste form object
! cwfm: pointer to current waste form mechanism object
! global_auxvars(:): pointer to global auxvar object, which stores the
!    temperature, liquid saturation, and liquid density, and is indexed
!    by the ghosted cell id
! material_auxvars(:): pointer to material auxvar object, which stores
!    the porosity and grid cell volume, and is indexed by the ghosted cell id
! field: pointer to field object
! option: pointer to option object
! grid: pointer to grid object
! dV: [m3] volume change
! dt: [sec] transport time step length
! avg_temp_local/global: [C] average temperature in waste form region
! local_id: [-] local grid cell id
! ghosted_id: [-] ghosted grid cell id
! idof: [-] degree of freedom number
! i, k, p, g, d, f: [-] looping index integers
! num_species: [-] number of readionuclide species objects
! ierr: [-] PETSc error integer
! Coeff(:): [mol-RN/g-bulk] coefficient in radionuclide (RN) decay equations
! concentration_old(:): [mol-RN/g-bulk] radionuclide concentration from
!    previous time step in decay equations
! inst_release_molality: [mol-RN/kg-water] instant release fraction of a
!    radionuclide (RN) in molality units
! conversion: [1/day] --> [1/sec]
! xx_p(:): [mol-RN/kg-water] pointer to solution vector for species
!    concentration
! norm: [-] norm calculation value
! residual(:): [mol/g-bulk/sec] residual array for implicit calculation
! solution(:): [mol/g-bulk] solution array for implicit calculation
! rhs(:): [mol/g-bulk/sec] right hand side array for implicit calculation
! indices(:): [-] array of indices
! Jacobian(:): [1/sec] Jacobian matrix for implicit calculation
! rate: [mol/g-bulk/sec] isotope mass decay rate
! rate_constant: [1/sec] isotope decay constant
! one_over_dt: [1/sec] helper variable to avoid dividing
! tolerance: [-] tolerance parameter for implicit calculation
! idaughter: [-] daughter integer number
! it: [-] iteration number for implicit calculation
! -----------------------------------------------------------
  class(waste_form_base_type), pointer :: cur_waste_form
  class(wf_mechanism_base_type), pointer :: cwfm
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  class(crit_mechanism_base_type), pointer :: cur_criticality
  PetscReal :: dV
  PetscReal :: dt
  PetscReal :: avg_temp_local, avg_temp_global
  PetscReal :: avg_sat_local, avg_sat_global
  PetscInt :: local_id, ghosted_id
  PetscInt :: idof
  PetscInt :: i, k, p, g, d, f, j
  PetscInt :: num_species
  PetscErrorCode :: ierr
  PetscReal, allocatable :: Coeff(:)
  PetscReal, allocatable :: concentration_old(:)
  PetscReal :: inst_release_molality
  PetscReal, pointer :: xx_p(:)
  ! implicit solution parameters
  PetscReal :: norm
  PetscReal, allocatable :: residual(:)
  PetscReal, allocatable :: solution(:)
  PetscReal, allocatable :: rhs(:)
  PetscInt, allocatable :: indices(:)
  PetscReal, allocatable :: Jacobian(:,:)
  PetscReal :: rate, rate_constant, one_over_dt
  PetscReal, parameter :: tolerance = 1.d-12
  PetscInt :: it, iiso, idaughter
! -----------------------------------------------------------
  PetscReal :: t_low, t_high
  PetscReal, pointer :: times(:)
  PetscBool :: dataset_solution
  PetscBool :: expanded_dataset_solution
  PetscBool :: switch_to_implicit ! warning message for switch to implicit soln
  class(dataset_ascii_type), pointer :: dataset
  class(crit_inventory_type), pointer :: crit_inventory
  type(crit_inventory_lookup_type), pointer :: inventory_table
  type(rad_species_type), pointer :: rad_species(:)

  avg_temp_global = UNINITIALIZED_DOUBLE
  avg_sat_global = UNINITIALIZED_DOUBLE

  global_auxvars => this%realization%patch%aux%Global%auxvars
  material_auxvars => this%realization%patch%aux%Material%auxvars
  field => this%realization%field
  option => this%option
  grid => this%realization%patch%grid
  dt = option%tran_dt
  ! zero entries from previous time step
  call VecZeroEntries(this%data_mediator%vec,ierr);CHKERRQ(ierr)

  if (associated(this%criticality_mediator)) then
    call VecZeroEntries(this%criticality_mediator%data_mediator%vec, &
                        ierr);CHKERRQ(ierr)
  endif

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    cwfm => cur_waste_form%mechanism
    num_species = cwfm%num_species
    allocate(Coeff(num_species))
    allocate(concentration_old(num_species))
    ! ------ update mass balances after transport step ---------------------
    select type(cwfm => cur_waste_form%mechanism)
      class is(wf_mechanism_dsnf_type)
        ! note: do nothing here because the cumulative mass update for dsnf
        ! and/or wipp mechanisms has already occured (if breached)
      class default
        cur_waste_form%cumulative_mass = cur_waste_form%cumulative_mass + &
          cur_waste_form%instantaneous_mass_rate*dt
    end select
    ! ------ update matrix volume ------------------------------------------
    select type(cwfm => cur_waste_form%mechanism)
      class is(wf_mechanism_dsnf_type)
        ! note: do nothing here because the volume update for dsnf and/or
        ! wipp mechanisms has already occured (if breached)
      class default
         dV = cur_waste_form%eff_dissolution_rate / &   ! kg-matrix/sec
           cwfm%matrix_density * &                      ! kg-matrix/m^3-matrix
           dt                                           ! sec
         cur_waste_form%volume = cur_waste_form%volume - dV
    end select
    if (cur_waste_form%volume <= 1.d-8) then
      cur_waste_form%volume = 0.d0
    endif

    ! ------ get species concentrations from mass fractions ----------------
    do k = 1,num_species
      if (cur_waste_form%volume <= 0.d0) then
        cur_waste_form%rad_concentration(k) = 0.d0
        cur_waste_form%rad_mass_fraction(k) = 0.d0
      else
        cur_waste_form%rad_concentration(k) = &
          cur_waste_form%rad_mass_fraction(k) / &
          cwfm%rad_species_list(k)%formula_weight
      endif
    enddo

    !---------------- vitality degradation function ------------------------
    if (cur_waste_form%canister_degradation_flag .and. &
        (cur_waste_form%canister_vitality > 1.d-3)) then
      if (.not.cur_waste_form%breached .and. &
          initialized(cur_waste_form%breach_time)) then
        ! do not modify eff_canister_vit_rate from what it was set to
        cur_waste_form%eff_canister_vit_rate = &
          cur_waste_form%eff_canister_vit_rate
      else
        avg_temp_local = 0.d0
        do i = 1,cur_waste_form%region%num_cells
          local_id = cur_waste_form%region%cell_ids(i)
          ghosted_id = grid%nL2G(local_id)
          avg_temp_local = avg_temp_local + global_auxvars(ghosted_id)%temp * &
                           cur_waste_form%scaling_factor(i)
        enddo
        call CalcParallelSUM(option,cur_waste_form%rank_list,avg_temp_local, &
                             avg_temp_global)
        avg_temp_global = avg_temp_global+273.15d0   ! Kelvin
        cur_waste_form%eff_canister_vit_rate = &
          cur_waste_form%canister_vitality_rate * &
          exp( cwfm%canister_material_constant * ( (1.d0/333.15d0) - &
          (1.d0/(avg_temp_global))) )
      endif
      cur_waste_form%canister_vitality = cur_waste_form%canister_vitality &
                                 - (cur_waste_form%eff_canister_vit_rate*dt)
      if (cur_waste_form%canister_vitality <= 1.d-3) then
        cur_waste_form%canister_vitality = 0.d0
        cur_waste_form%eff_canister_vit_rate = 0.d0
        cur_waste_form%canister_vitality_rate = 0.d0
      endif
    endif

    !----------------- spacer degradation function -------------------------
    if (cur_waste_form%spacer_degradation_flag .and. &
       (cur_waste_form%spacer_vitality > 1.d-2) .and. &
       cur_waste_form%breached) then

        ! Get average saturation
        if (.not. Initialized(avg_sat_global)) then
          avg_sat_local = 0.d0
          do i = 1,cur_waste_form%region%num_cells
            local_id = cur_waste_form%region%cell_ids(i)
            ghosted_id = grid%nL2G(local_id)
            avg_sat_local = avg_sat_local + &
                            global_auxvars(ghosted_id)%sat(LIQUID_PHASE) * &
                            cur_waste_form%scaling_factor(i)
          enddo
          call CalcParallelSUM(option,cur_waste_form%rank_list,avg_sat_local, &
                               avg_sat_global)
        endif

        if (avg_sat_global >= &
            cur_waste_form%spacer_mechanism%threshold_sat) then
          cur_waste_form%spacer_mechanism%alteration_rate = 1.0d0
        else
          cur_waste_form%spacer_mechanism%alteration_rate = avg_sat_global / &
            cur_waste_form%spacer_mechanism%threshold_sat
        endif

        ! Get average temperature
        if ( .not. Initialized(avg_temp_global)) then
          avg_temp_local = 0.d0
          do i = 1,cur_waste_form%region%num_cells
            local_id = cur_waste_form%region%cell_ids(i)
            ghosted_id = grid%nL2G(local_id)
            avg_temp_local = avg_temp_local + global_auxvars(ghosted_id)%temp* &
                             cur_waste_form%scaling_factor(i)
          enddo
          call CalcParallelSUM(option,cur_waste_form%rank_list,avg_temp_local, &
                               avg_temp_global)
          avg_temp_global = avg_temp_global + 273.15d0   ! Kelvin
        endif

      call cur_waste_form%spacer_mechanism%Degradation(cur_waste_form, this, &
                                                       avg_sat_global, &
                                                       avg_temp_global, &
                                                       dt, ierr)
    endif

    ! ------------------ criticality termination criterion -----------------
    if (cur_waste_form%spacer_degradation_flag .and. &
       (cur_waste_form%spacer_vitality <= 1.d-2)) then
      if (associated(this%criticality_mediator)) then
        cur_criticality => this%criticality_mediator%crit_mech_list
        do
          if (.not. associated(cur_criticality)) exit
          if (StringCompare(cur_criticality%mech_name, &
                            cur_waste_form%criticality_mech_name)) then
            if (option%time < cur_criticality%crit_event%crit_end) then
              cur_criticality%crit_event%crit_end = option%time
            endif
          endif
          cur_criticality => cur_criticality%next
        enddo
      endif
      cur_waste_form%spacer_degradation_flag = PETSC_FALSE
    endif

    !------- instantaneous release -----------------------------------------
    if ((.not.cur_waste_form%breached .and. &
         cur_waste_form%canister_vitality < 1.d-3) .or. &
        (.not.cur_waste_form%breached .and. &
         initialized(cur_waste_form%breach_time) .and. &
         option%time > cur_waste_form%breach_time)) then
      call VecGetArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
      do k = 1,num_species
        cur_waste_form%inst_release_amount(k) = &
           (cwfm%rad_species_list(k)%inst_release_fraction * &
            cur_waste_form%rad_concentration(k))

        ! Update cumulative release (mol) to include instantaneous release
        cur_waste_form%cumulative_mass(k) = cur_waste_form% &
                  cumulative_mass(k) + cur_waste_form%inst_release_amount(k) * &
                  cur_waste_form%volume * cwfm%matrix_density * 1.d3

        cur_waste_form%rad_concentration(k) = &
           cur_waste_form%rad_concentration(k) - &
           cur_waste_form%inst_release_amount(k)
        ! update mass fractions after instantaneous release
        cur_waste_form%rad_mass_fraction(k) = &
           cur_waste_form%rad_concentration(k) * &
           cwfm%rad_species_list(k)%formula_weight
        ! update transport solution vector with mass injection molality
        ! as an alternative to a source term (issue with tran_dt changing)
        do f = 1, cur_waste_form%region%num_cells
          local_id = cur_waste_form%region%cell_ids(f)
          ghosted_id = grid%nL2G(local_id)
          inst_release_molality = &                    ! [mol-rad/kg-water]
            ! [mol-rad]
            (cur_waste_form%inst_release_amount(k) * & ! [mol-rad/g-matrix]
             cur_waste_form%volume * &                 ! [m^3-matrix]
             cwfm%matrix_density * &                   ! [kg-matrix/m^3-matrix]
             1.d3) / &                               ! [kg-matrix] -> [g-matrix]
             ! [kg-water]
            (material_auxvars(ghosted_id)%porosity * &         ! [-]
             global_auxvars(ghosted_id)%sat(LIQUID_PHASE) * &  ! [-]
             material_auxvars(ghosted_id)%volume * &           ! [m^3]
             global_auxvars(ghosted_id)%den_kg(LIQUID_PHASE))  ! [kg/m^3-water]
          idof = cwfm%rad_species_list(k)%ispecies + &
                 ((local_id - 1) * option%ntrandof)
          xx_p(idof) = xx_p(idof) + &
                       (inst_release_molality*cur_waste_form%scaling_factor(f))
        enddo

      enddo

      cur_waste_form%breached = PETSC_TRUE
      cur_waste_form%breach_time = option%time
      call VecRestoreArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
    endif

    ! Save the concentration after inst. release for the decay step
    concentration_old = cur_waste_form%rad_concentration

    if ((cur_waste_form%volume >= 0.d0) .and. &
        (option%time >= cur_waste_form%decay_start_time)) then !--------------

    dataset_solution = PETSC_FALSE
    expanded_dataset_solution = PETSC_FALSE
    switch_to_implicit = PETSC_FALSE
    if (associated(cur_waste_form%criticality_mechanism)) then
      cur_criticality => cur_waste_form%criticality_mechanism
      if (associated(cur_criticality%rad_dataset)) then
        dataset_solution = PETSC_TRUE
        if (this%implicit_solution) this%implicit_solution = PETSC_FALSE
      endif
      if (associated(cur_criticality%inventory_dataset)) then
        crit_inventory => cur_criticality%inventory_dataset
        if (crit_inventory%switch_implicit) then
          ! extrapolation detected in the expanded dataset - switch to implicit
          this%implicit_solution = PETSC_TRUE
        elseif (.not. crit_inventory%continue_lookup .and. &
                option%time > cur_criticality%crit_event%crit_end) then
          ! use implicit solution if the criticality end time is exceeded
          this%implicit_solution = PETSC_TRUE
          crit_inventory%switch_implicit = PETSC_TRUE ! stop using lookup table
          switch_to_implicit =  PETSC_TRUE ! give warning message
        else
          ! use 3D lookup table (expanded dataset) for radionuclide inventories
          expanded_dataset_solution = PETSC_TRUE
          if (this%implicit_solution) this%implicit_solution = PETSC_FALSE
        endif
      endif
    endif

    if (expanded_dataset_solution) then
      ! interpolate radionuclide inventories using lookup tables from external
      !   neutronics calculations

      ! check number of species
      if (num_species > crit_inventory%num_species) then
        option%io_buffer = 'Number of species listed in the criticality ' &
                         //'inventory lookup table "' &
                         //trim(crit_inventory%file_name) &
                         //'" is less than the number specified ' &
                         //'in Waste Form Process Model.'
        call PrintErrMsg(option)
      endif

      rad_species => cwfm%rad_species_list
      do j = 1, num_species
        k = 0
        inventory_table => crit_inventory%radionuclide_table
        do
          if (.not. associated (inventory_table)) exit

          ! find index for cur_waste_form%rad_mass_fraction(:)
          if (rad_species(j)%name == inventory_table%name) then
            k = j

            cur_waste_form%rad_mass_fraction(k) = &
              inventory_table%Evaluate(cur_criticality%crit_event%crit_start, &
                                       cur_criticality%crit_heat, &
                                       option%time)

            if (inventory_table%lookup%axis3%extrapolate .and. &
                .not. crit_inventory%allow_extrap) then
             ! Fallback options if extrapolation is detected
             if (crit_inventory%allow_implicit) then
               ! Resort to implicit solution
               this%implicit_solution = PETSC_TRUE
               crit_inventory%switch_implicit = PETSC_TRUE
               switch_to_implicit =  PETSC_TRUE
             else
               ! Alert user that extrapolation has been attempted
               option%io_buffer = 'Extrapolation of inventory lookup table "' &
                                // trim(crit_inventory%file_name) // '" ' &
                                //'has been detected. Extrapolation can be ' &
                                //'enabled with keyword '&
                                //'USE_LOOKUP_AND_EXTRAPOLATION in the OPTION '&
                                //'sub-block.'
               call PrintErrMsg(option)
             endif
            else
              ! Modify the radionuclide concentration
              cur_waste_form%rad_concentration(k) = &
                cur_waste_form%rad_mass_fraction(k) / &
                cwfm%rad_species_list(k)%formula_weight
            endif
          endif

          inventory_table => inventory_table%next
        enddo
        if (k == 0) then
          option%io_buffer = 'Radionuclide "'// trim(rad_species(j)%name) &
                           //'" in the waste form species list was not found ' &
                           //'in criticality inventory lookup table "' &
                           //trim(crit_inventory%file_name) // '"'
          call PrintErrMsg(option)
        endif
      enddo
      nullify(inventory_table)

    elseif (dataset_solution) then
      !Import radionuclide inventory from external neutronics code calculations
      dataset => cur_waste_form%criticality_mechanism%rad_dataset
      times => dataset%time_storage%times
      if (num_species > dataset%dims(1)) then
        option%io_buffer = 'Number of species in dataset is less than ' // &
                     'the number specified in Waste Form Process Model.'
        call PrintErrMsg(option)
      elseif (num_species < dataset%dims(1)) then
        option%io_buffer = 'Number of species in dataset is greater than ' //&
                       'the number specified in Waste Form Process Model.'
        call PrintErrMsg(option)
      endif
      j=1
      t_low = times(j)
      t_high = t_low
      do
        if(option%time < times(j)) exit
        if(j == size(times)) exit
        t_low = times(j)
        j = j+1
        t_high = times(j)
      enddo

      do k = 1,num_species
        if (j == size(times) .and. option%time >= times(j)) then
          cur_waste_form%rad_mass_fraction(k) = dataset%rbuffer((j-1)* &
                  num_species+k)
        elseif (j==1) then
          cur_waste_form%rad_mass_fraction(k) = 1.d-20
        else
          cur_waste_form%rad_mass_fraction(k) = dataset%rbuffer((j-2)* &
                num_species+k) + (option%time-t_low)/(t_high-t_low)* &
                (dataset%rbuffer((j-2)*num_species+k+num_species)-dataset% &
                rbuffer((j-2)*num_species+k))
       endif
       cur_waste_form%rad_concentration(k) = &
          cur_waste_form%rad_mass_fraction(k) / &
          cwfm%rad_species_list(k)%formula_weight
      enddo

    elseif (.not.this%implicit_solution) then !-----------------------------------

    ! 3-generation analytical solution derived for multiple parents and
    ! grandparents and non-zero initial daughter concentrations (see Section
    ! 3.2.3 of Mariner et al. (2016), SAND2016-9610R), where the solution is
    ! obtained explicitly in time

      !------- decay the radionuclide species --------------------------------
      ! FIRST PASS =====================
      do d = 1,num_species
        ! Update the initial value of the species coefficient
        Coeff(d) = cur_waste_form%rad_concentration(d)
        do p = 1,num_species
          ! If the daughter has a parent(s):
          if (d == cwfm%rad_species_list(p)%daugh_id) then
            Coeff(d) = Coeff(d) - &
              (cwfm%rad_species_list(p)%decay_constant * &
               concentration_old(p)) / &
              (cwfm%rad_species_list(d)%decay_constant - &
               cwfm%rad_species_list(p)%decay_constant)
            do g = 1,num_species
              ! If the daughter has a grandparent(s):
              if (p == cwfm%rad_species_list(g)%daugh_id) then
                Coeff(d) = Coeff(d) - &
                  ((cwfm%rad_species_list(p)%decay_constant* &
                    cwfm%rad_species_list(g)%decay_constant* &
                    concentration_old(g)) / &
                  ((cwfm%rad_species_list(p)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant)* &
                   (cwfm%rad_species_list(d)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant))) + &
                  ((cwfm%rad_species_list(p)%decay_constant* &
                    cwfm%rad_species_list(g)%decay_constant* &
                    concentration_old(g)) / &
                  ((cwfm%rad_species_list(p)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant)* &
                   (cwfm%rad_species_list(d)%decay_constant - &
                    cwfm%rad_species_list(p)%decay_constant)))
              endif
            enddo ! grandparent loop
          endif
        enddo ! parent loop
      enddo
      ! SECOND PASS ====================
      do d = 1,num_species
        ! Decay the species
        cur_waste_form%rad_concentration(d) = Coeff(d) * exp(-1.d0 * &
          cwfm%rad_species_list(d)%decay_constant * dt)
        do p = 1,num_species
          ! If the daughter has a parent(s):
          if (d == cwfm%rad_species_list(p)%daugh_id) then
            cur_waste_form%rad_concentration(d) = &
              cur_waste_form%rad_concentration(d) + &
              (((cwfm%rad_species_list(p)%decay_constant* &
                 concentration_old(p)) / &
                (cwfm%rad_species_list(d)%decay_constant - &
                 cwfm%rad_species_list(p)%decay_constant)) * &
               exp(-1.d0 * cwfm%rad_species_list(p)%decay_constant * dt))
            do g = 1,num_species
              ! If the daughter has a grandparent(s):
              if (p == cwfm%rad_species_list(g)%daugh_id) then
                cur_waste_form%rad_concentration(d) = &
                  cur_waste_form%rad_concentration(d) - &
                  ((cwfm%rad_species_list(p)%decay_constant* &
                    cwfm%rad_species_list(g)%decay_constant* &
                    concentration_old(g)*exp(-1.d0* &
                    cwfm%rad_species_list(p)%decay_constant*dt)) / &
                  ((cwfm%rad_species_list(p)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant)* &
                   (cwfm%rad_species_list(d)%decay_constant - &
                    cwfm%rad_species_list(p)%decay_constant))) + &
                  ((cwfm%rad_species_list(p)%decay_constant* &
                    cwfm%rad_species_list(g)%decay_constant* &
                    concentration_old(g)*exp(-1.d0* &
                    cwfm%rad_species_list(g)%decay_constant*dt)) / &
                  ((cwfm%rad_species_list(p)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant)* &
                   (cwfm%rad_species_list(d)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant)))
              endif
            enddo ! grandparent loop
          endif
        enddo ! parent loop
      enddo

    endif !--------------------------------------------------------------------

    ! Warning message for switch to implicit solution from lookup table
    if (switch_to_implicit) then
      if (associated(cur_waste_form%criticality_mechanism%inventory_dataset)) &
        then
        crit_inventory => cur_waste_form%criticality_mechanism%inventory_dataset
        option%io_buffer = 'Switching from inventory lookup table "' &
                         // trim(crit_inventory%file_name) &
                         //'" to implicit solution.'
        call PrintWrnMsg(option)
      endif
    endif

    if (this%implicit_solution) then

    ! implicit solution based on Bateman's equations and Newton's method

    allocate(residual(num_species))
    allocate(solution(num_species))
    allocate(rhs(num_species))
    allocate(indices(num_species))
    allocate(Jacobian(num_species,num_species))

    residual = 1.d0 ! to start, must set bigger than tolerance
    ! to start, set solution to old concentration, but ensure it's not zero
    solution = max(concentration_old,1.d-40)
    one_over_dt = 1.d0/dt
    it = 0
    do ! nonlinear loop
      if (dot_product(residual,residual) < tolerance) exit ! 2-norm(residual)
      it = it + 1
      residual = 0.d0 ! set to zero because we are summing
      ! f(C^{k+1,p}) = (C^{k+1,p} - C^k)/dt -R(C^{k+1,p})
      Jacobian = 0.d0 ! set to zero because we are summing
      ! J_ij = del[f_i(C^{k+1,p})]/del[C_j^{k+1,p}]
      ! isotope loop
      do iiso = 1, num_species
        ! ----accumulation term for isotope------------------------!-units--
        ! dC/dt = (C^{k+1,p} - C^k)/dt
        residual(iiso) = residual(iiso) + &                        ! mol/g/sec
                     (solution(iiso) - concentration_old(iiso)) * &! mol/g
                     one_over_dt                                   ! 1/sec
        ! d[(C^{k+1,p} - C^k)/dt]/d[C^{k+1,p}] = 1/dt
        Jacobian(iiso,iiso) = Jacobian(iiso,iiso) + &              ! 1/sec
                              one_over_dt                          ! 1/sec
        ! ----source/sink term for isotope-------------------------!-units--
        ! -R(C^{k+1,p}) = -(-L*(C^{k+1,p}))    L=lambda
        rate_constant = cwfm%rad_species_list(iiso)%decay_constant ! 1/sec
        rate = rate_constant * solution(iiso)                      ! mol/g/sec
        residual(iiso) = residual(iiso) + rate                     ! mol/g/sec
        ! d[-(-L*(C^{k+1,p}))]/d[C^{k+1,p}] = L
        Jacobian(iiso,iiso) = Jacobian(iiso,iiso) + rate_constant  ! 1/sec
        ! daughter loop (there is only one daughter currently)
        !do i = 1, this%isotope_daughters(0,iiso)
        ! ----source/sink term for daughter----------------------!-units--
        idaughter = cwfm%rad_species_list(iiso)%daugh_id
        if (idaughter > 0) then
          ! -R(C^{k+1,p}) = -(L*S*(C^{k+1,p}))    L=lambda
          residual(idaughter) = residual(idaughter) - rate         ! mol/g/sec
          ! d[-(L*S*(C^{k+1,p}))]/d[C^{k+1,p}] = -L*S
          Jacobian(idaughter,iiso) = Jacobian(idaughter,iiso) - &  ! 1/sec
                                     rate_constant                 ! 1/sec
        endif
        !enddo
        ! k=time, p=iterate, C=isotope concentration
      enddo
      ! scale Jacobian
      do iiso = 1, num_species
        norm = max(1.d0,maxval(abs(Jacobian(iiso,:))))
        norm = 1.d0/norm
        rhs(iiso) = residual(iiso)*norm
        ! row scaling
        Jacobian(iiso,:) = Jacobian(iiso,:)*norm
      enddo
      ! log formulation for derivatives, column scaling
      do iiso = 1, num_species
        Jacobian(:,iiso) = Jacobian(:,iiso)*solution(iiso)
      enddo
      ! linear solve steps
      ! solve step 1/2: get LU decomposition
      call LUDecomposition(Jacobian,num_species,indices,i)
      ! solve step 2/2: LU back substitution linear solve
      call LUBackSubstitution(Jacobian,num_species,indices,rhs)
      rhs = dsign(1.d0,rhs)*min(dabs(rhs),10.d0)
      ! update the solution
      solution = solution*exp(-rhs)
    enddo
    ! if old_concentration was zero, set it back to zero here
    do iiso = 1,num_species
      if (solution(iiso) <= 1.d-40) solution(iiso) = 0.d0
    enddo
    cur_waste_form%rad_concentration = solution

    deallocate(residual)
    deallocate(solution)
    deallocate(rhs)
    deallocate(indices)
    deallocate(Jacobian)

    endif !-----implicit_solution---------------------------------------------

    ! ------ update species mass fractions ---------------------------------
    do k = 1,num_species
      cur_waste_form%rad_mass_fraction(k) = &       ! [g-rad/g-wf]
      cur_waste_form%rad_concentration(k) * &       ! [mol-rad/g-wf]
        cur_waste_form%mechanism%rad_species_list(k)%formula_weight
      ! to avoid errors in plotting data when conc is very very low:
      if (cur_waste_form%rad_mass_fraction(k) <= 1d-40) then
        cur_waste_form%rad_mass_fraction(k) = 0.d0
      endif
    enddo

    endif !-------------------------------------------------------------------

    deallocate(concentration_old)
    deallocate(Coeff)
    cur_waste_form => cur_waste_form%next
  enddo

  if (this%print_mass_balance) then
    call PMWFOutput(this)
  endif

end subroutine PMWFInitializeTimestep

! ************************************************************************** !

subroutine PMWFSolve(this,time,ierr)
  !
  ! Updates the source term based on the dissolution model chosen
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Updated/modified by Jenn Frederick 04/2016
  !
  ! Notes: The species loop must be the inner loop, while the grid cell loop
  ! must be the outer loop, in order for the vec_p(i) indexing to work.

  use Global_Aux_module
  use Material_Aux_module
  use Reactive_Transport_Aux_module, only : rt_min_saturation
  use Grid_module
  use Option_module
  use Utility_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! time (input): [sec] simulation time
! ierr (input/output): [-] PETSc error integer
! ---------------------------------
  class(pm_waste_form_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
! ---------------------------------

! LOCAL VARIABLES:
! ================
! cur_waste_form: pointer to current waste form object
! i, j, k: [-] looping index integers
! num_species: [-] number of radionuclide species objects
! local_id: [-] local grid cell id number
! ghosted_id: [-] ghosted grid cell id number
! idof: [-] degree of freedom index
! inst_diss_molality: [mol-RN/kg-water] molality of radionuclide (RN) that
!    is dissolved from waste form
! vec_p(:): [mol-RN/sec] pointer to data mediator array that stores the
!    radionuclide source term
! xx_p(:): [mol-RN/kg-water] pointer to solution vector for species
!    concentration
! fmdm_count_*: [-] number of time the FMDM is called
! word: temporary string
! global_auxvars(:): pointer to global auxvar object, which stores the
!    liquid saturation and liquid density at each grid cell, indexed by
!    the ghosted cell id
! material_auxvars(:): pointer to material auxvar object, which stores the
!    porosity and grid cell volume at each grid cell, indexed by the
!    ghosted cell id
! grid: pointer to the grid object
! -----------------------------------------------------------
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: i, j, k, m
  PetscInt :: num_species
  PetscInt :: local_id, ghosted_id
  PetscInt :: idof
  PetscReal :: inst_diss_molality
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: xx_p(:), heat_source(:)
  PetscReal :: avg_temp_local, avg_temp_global
  PetscReal :: avg_sw_local, avg_sw_global
  PetscReal :: avg_rho_w_local, avg_rho_w_global
  PetscInt :: fmdm_count_global, fmdm_count_local
  PetscLogDouble :: log_start_time, log_end_time
  character(len=MAXWORDLENGTH) :: word
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type), pointer :: grid
  class(crit_mechanism_base_type), pointer :: cur_criticality
  type(option_type), pointer :: option
! -----------------------------------------------------------

  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

  fmdm_count_global = 0
  fmdm_count_local = 0
  avg_temp_local = 0.d0
  avg_sw_local  = 0.d0
  avg_rho_w_local = 0.d0
  avg_temp_global = UNINITIALIZED_DOUBLE
  avg_sw_global  = UNINITIALIZED_DOUBLE
  avg_rho_w_global = UNINITIALIZED_DOUBLE
  global_auxvars => this%realization%patch%aux%Global%auxvars
  material_auxvars => this%realization%patch%aux%Material%auxvars
  grid => this%realization%patch%grid

  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(this%realization%field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  if (associated(this%criticality_mediator)) then
    call VecGetArrayF90(this%criticality_mediator%data_mediator%vec, &
                        heat_source,ierr);CHKERRQ(ierr)
  endif

  cur_waste_form => this%waste_form_list
  i = 0
  m = 0
  do
    if (.not.associated(cur_waste_form)) exit
    num_species = cur_waste_form%mechanism%num_species
    !
    if ((cur_waste_form%volume > 0.d0) .and. &
        (cur_waste_form%canister_vitality <= 1.d-40)) then
    !---------------------------------------------------------------------------
      ! calculate the mechanism-specific eff_dissolution_rate [kg-bulk/sec]:
      call cur_waste_form%mechanism%Dissolution(cur_waste_form,this,ierr)
      select type(cwfm => cur_waste_form%mechanism)
      !-----------------------------------------------------------------------
        ! ignore source term if dsnf/wipp type, and directly update the
        ! solution vector instead (see note in WFMech[DSNF/WIPP]Dissolution):
        class is(wf_mechanism_dsnf_type)
          do k = 1,cur_waste_form%region%num_cells
            local_id = cur_waste_form%region%cell_ids(k)
            ghosted_id = grid%nL2G(local_id)
            if (global_auxvars(ghosted_id)%sat(LIQUID_PHASE) < &
                rt_min_saturation) then
              cycle
            endif
            do j = 1,num_species
              i = i + 1
              cur_waste_form%instantaneous_mass_rate(j) = &   ! mol-rad/sec
                (cur_waste_form%eff_dissolution_rate / &      ! kg-bulk/sec
                 cur_waste_form%mechanism%rad_species_list(j)%formula_weight * &! kg-rad/kmol-rad
                 cur_waste_form%rad_mass_fraction(j) * &      ! kg-rad/kg-bulk
                 1.d3)
              inst_diss_molality = &                          ! mol-rad/kg-water
                cur_waste_form%instantaneous_mass_rate(j) * & ! mol-rad/sec
                this%realization%option%tran_dt / &           ! sec
                ! [kg-water]
                (material_auxvars(ghosted_id)%porosity * &        ! [-]
                 global_auxvars(ghosted_id)%sat(LIQUID_PHASE) * & ! [-]
                 material_auxvars(ghosted_id)%volume * &          ! [m^3]
                 global_auxvars(ghosted_id)%den_kg(LIQUID_PHASE)) ! [kg/m^3-water]
              idof = cwfm%rad_species_list(j)%ispecies + &
                     ((local_id - 1) * this%option%ntrandof)
              xx_p(idof) = xx_p(idof) + &                     ! mol-rad/kg-water
                           (inst_diss_molality*cur_waste_form%scaling_factor(k))
              vec_p(i) = 0.d0
              if (k == 1) then
                ! update the cumulative mass now, not at next timestep:
                cur_waste_form%cumulative_mass(j) = &            ! mol-rad
                  cur_waste_form%cumulative_mass(j) + &          ! mol-rad
                  cur_waste_form%instantaneous_mass_rate(j) * &  ! mol-rad/sec
                  this%realization%option%tran_dt                ! sec
                ! update the volume now, not at next timestep:
                cur_waste_form%volume = 0.d0
              endif
            enddo
          enddo
      !-----------------------------------------------------------------------
        ! for all other waste form types, load the source term, and update
        ! the cumulative mass and volume at next timestep:
        class default
          do k = 1,cur_waste_form%region%num_cells
            do j = 1,num_species
              i = i + 1
              cur_waste_form%instantaneous_mass_rate(j) = &   ! mol-red/sec
                (cur_waste_form%eff_dissolution_rate / &      ! kg-matrix/sec
                 cur_waste_form%mechanism%rad_species_list(j)%formula_weight * &! kg-rad/kmol-rad
                 cur_waste_form%rad_mass_fraction(j) * &      ! kg-rad/kg-matrix
                 1.d3)
              vec_p(i) = cur_waste_form%instantaneous_mass_rate(j) * &
                         cur_waste_form%scaling_factor(k)  ! mol-rad/sec * [-]
            enddo
          enddo
      !-------------------------------------------------------------------------
      end select
      ! count the number of times FMDM was called:
      select type(cwfm => cur_waste_form%mechanism)
        type is(wf_mechanism_fmdm_type)
          fmdm_count_local = fmdm_count_local + 1
        type is (wf_mechanism_fmdm_surrogate_type)
          fmdm_count_local = fmdm_count_local + 1
      end select
    !---------------------------------------------------------------------------
    else ! (canister not breached, or all waste form has dissolved already)
      vec_p((i+1):(i+num_species*cur_waste_form%region%num_cells)) = 0.d0
      i = i + num_species*cur_waste_form%region%num_cells
      cur_waste_form%eff_dissolution_rate = 0.d0
      cur_waste_form%instantaneous_mass_rate = 0.d0
    endif

    !---------- Criticality Calculation --------------------------------------
    if (associated(cur_waste_form%criticality_mechanism)) then
      cur_criticality => cur_waste_form%criticality_mechanism
      if (time >= cur_criticality%crit_event%crit_start .and. time < &
              cur_criticality%crit_event%crit_end .and. &
              cur_waste_form%spacer_vitality > 0.d0) then
        cur_criticality%crit_event%crit_flag = PETSC_TRUE
      else
        cur_criticality%crit_event%crit_flag = PETSC_FALSE
      endif

      call CriticalityCalc(cur_criticality,time,ierr)

      ! Use heat of criticality lookup table if applicable
      if (associated(cur_criticality%crit_heat_dataset) .and. &
          cur_criticality%crit_event%crit_flag) then
        do j = 1,cur_waste_form%region%num_cells
          ghosted_id = grid%nL2G(cur_waste_form%region%cell_ids(j))
          avg_temp_local = avg_temp_local + &  ! Celsius
            (global_auxvars(ghosted_id)%temp*cur_waste_form%scaling_factor(j))
        enddo
        call CalcParallelSUM(option,cur_waste_form%rank_list,avg_temp_local, &
                             avg_temp_global)
        avg_temp_global = avg_temp_global / size(cur_waste_form%rank_list)
        cur_criticality%temperature = avg_temp_global
        cur_criticality%crit_heat = cur_criticality%crit_heat_dataset% &
          Evaluate(cur_criticality%crit_event%crit_start, &
                   cur_criticality%temperature)
      endif

      ! Criticality termination - water saturation
      if (Initialized(cur_criticality%sw) .and. &
          cur_criticality%crit_event%crit_flag) then
        do j = 1,cur_waste_form%region%num_cells
          ghosted_id = grid%nL2G(cur_waste_form%region%cell_ids(j))
          avg_sw_local = avg_sw_local + &  ! Liquid Saturation
            (global_auxvars(ghosted_id)%sat(LIQUID_PHASE) * &
            cur_waste_form%scaling_factor(j))
        enddo
        call CalcParallelSUM(option,cur_waste_form%rank_list,avg_sw_local, &
                             avg_sw_global)
        avg_sw_global = avg_sw_global / size(cur_waste_form%rank_list)
        if (avg_sw_global < cur_criticality%sw) then
          cur_criticality%crit_event%crit_flag = PETSC_FALSE
        endif
      endif

      ! Criticality termination - water density
      if (Initialized(cur_criticality%rho_w) .and. &
          cur_criticality%crit_event%crit_flag) then
        do j = 1,cur_waste_form%region%num_cells
          ghosted_id = grid%nL2G(cur_waste_form%region%cell_ids(j))
          avg_rho_w_local = avg_rho_w_local + &  ! Liquid Density
            (global_auxvars(ghosted_id)%den_kg(LIQUID_PHASE) * &
            cur_waste_form%scaling_factor(j))
        enddo
        call CalcParallelSUM(option,cur_waste_form%rank_list,avg_rho_w_local, &
                             avg_rho_w_global)
        avg_rho_w_global = avg_rho_w_global / size(cur_waste_form%rank_list)
        if (avg_rho_w_global < cur_criticality%rho_w) then
          cur_criticality%crit_event%crit_flag = PETSC_FALSE
        endif
      endif

      ! Define heat source
      do j = 1, cur_waste_form%region%num_cells
        m = m + 1
        heat_source(m) = cur_criticality%decay_heat
        if (cur_criticality%crit_event%crit_flag) then
          heat_source(m) = heat_source(m) + cur_criticality%crit_heat
        endif
        ! Distribute heat source throughout all cells in a waste package
        heat_source(m) = heat_source(m) * cur_waste_form%scaling_factor(j)
      enddo

    endif

    cur_waste_form => cur_waste_form%next
  enddo

  ! ideally, this print statement would go inside the dissolution subroutine
  call MPI_Allreduce(fmdm_count_local,fmdm_count_global,ONE_INTEGER_MPI, &
                     MPI_INTEGER,MPI_SUM,this%realization%option%mycomm, &
                     ierr);CHKERRQ(ierr)
  if ((fmdm_count_global > 0) .and. &
      this%realization%option%print_screen_flag) then
    write(word,'(i5)') fmdm_count_global
  ! ** START (this can be removed after FMDM profiling is finished) **
    write(*,'(a)') '== ' // adjustl(trim(word)) // ' call(s) to FMDM.'
  ! ** END (this can be removed after FMDM profiling is finished) **
  endif

  if (associated(this%criticality_mediator)) then
    call VecRestoreArrayF90(this%criticality_mediator%data_mediator%vec, &
                            heat_source,ierr);CHKERRQ(ierr)
  endif
  call VecRestoreArrayF90(this%realization%field%tran_xx,xx_p, &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)

  call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

  this%cumulative_time = this%cumulative_time + (log_end_time - log_start_time)

end subroutine PMWFSolve

! ************************************************************************** !

subroutine WFMechBaseDissolution(this,waste_form,pm,ierr)
  !
  ! Calculates the waste form dissolution rate; must be extended
  !
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): base mechanism object
! waste_form (input/output): base waste form object
! pm (input/output): waste form process model
! ierr (input/output): [-] PETSc error integer
! -----------------------------------------
  class(wf_mechanism_base_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr
! -----------------------------------------

  ! This routine must be extended.
  print *, 'subroutine WFMechBaseDissolution must be extended!'
  stop

end subroutine WFMechBaseDissolution

! ************************************************************************** !

subroutine WFMechGlassDissolution(this,waste_form,pm,ierr)
  !
  ! Calculates the glass waste form dissolution rate
  !
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  use Grid_module
  use Utility_module
  use Global_Aux_module
  use String_module
  use Reactive_Transport_Aux_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): base mechanism object
! waste_form (input/output): base waste form object
! pm (input/output): waste form process model
! ierr (input/output): [-] PETSc error integer
! -----------------------------------------
  class(wf_mechanism_glass_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr
! -----------------------------------------

! LOCAL VARIABLES:
! ================
! grid: pointer to grid object
! global_auxvars(:): pointer to global auxvar object, which stores the
!    temperature at each grid cell, and is indexed by the ghosted cell id
! rt_auxvars(:): pointer to the reactive transport auxvar object, which
!    stores the molality of each species and the activity coefficients of
!    each species, and is indexed by the ghosted cell id
! time_conversion: [1/day] --> [1/sec]
! avg_temp_local/global: [C] average temperature in waste form region
! ph_local/global: [-] average pH in waste form region
! Q_local/global: [-] average activity coefficient for H4SiO4 in waste
!    form region
! ghosted_id: [-] ghosted grid cell id
! i: [-] looping index integer
! --------------------------------------------------------------
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscReal :: avg_temp_local, avg_temp_global
  PetscReal :: ph_local, ph_global
  PetscReal :: Q_local, Q_global
  PetscInt :: ghosted_id
  PetscInt :: i
! --------------------------------------------------------------

  grid => pm%realization%patch%grid
  global_auxvars => pm%realization%patch%aux%Global%auxvars
  rt_auxvars => pm%realization%patch%aux%RT%auxvars

  ierr = 0

  ! Glass dissolution equation: Kienzler et al. (2012) Eq. 6 pg. 17
  ! Kienzler, B., M. Altmaier, et al. (2012) Radionuclide Source Term form
  ! HLW Glass, Spent Nuclear Fuel, and Compacted Hulls and End Pieces
  ! (CSD-C Waste). KIT Scientific Reports 7624. Karlsruhe Institute of
  ! Technology, Baden-Wurttemberg, Germany.
  ! Generalized glass dissolution equation comes from Eq. 2.3.7-6 in
  ! Yucca Mountain Repository SAR, Section 2.3.7, DOE/RW-0573 Rev.0

  avg_temp_local = 0.d0
  do i = 1,waste_form%region%num_cells
    ghosted_id = grid%nL2G(waste_form%region%cell_ids(i))
    avg_temp_local = avg_temp_local + &  ! Celsius
               (global_auxvars(ghosted_id)%temp * waste_form%scaling_factor(i))
  enddo
  call CalcParallelSUM(pm%option,waste_form%rank_list,avg_temp_local, &
                       avg_temp_global)
  avg_temp_global = avg_temp_global+273.15d0   ! Kelvin

  if (this%use_pH) then  ! pH ------------------------------------------------
    if (this%h_ion_id > 0) then   ! primary species
      ph_local = 0.d0
      do i = 1,waste_form%region%num_cells
        ghosted_id = grid%nL2G(waste_form%region%cell_ids(i))
        ph_local = ph_local + &
               ( -log10((rt_auxvars(ghosted_id)%pri_molal(this%h_ion_id) * &
                         rt_auxvars(ghosted_id)%pri_act_coef(this%h_ion_id) * &
                         waste_form%scaling_factor(i))) )
      enddo
    elseif (this%h_ion_id < 0) then   ! secondary species
      ph_local = 0.d0
      do i = 1,waste_form%region%num_cells
        ghosted_id = grid%nL2G(waste_form%region%cell_ids(i))
        ph_local = ph_local + &
               ( -log10((rt_auxvars(ghosted_id)%sec_molal(this%h_ion_id) * &
                         rt_auxvars(ghosted_id)%sec_act_coef(this%h_ion_id) * &
                         waste_form%scaling_factor(i))) )
      enddo
    endif
    call CalcParallelSUM(pm%option,waste_form%rank_list,ph_local,ph_global)
    this%pH = ph_global
  endif ! pH -----------------------------------------------------------------

  if (this%use_Q) then  ! Q --------------------------------------------------
    if (this%SiO2_id > 0) then   ! primary species
      Q_local = 0.d0
      do i = 1,waste_form%region%num_cells
        ghosted_id = grid%nL2G(waste_form%region%cell_ids(i))
        Q_local = Q_local + &
                 (rt_auxvars(ghosted_id)%pri_molal(abs(this%SiO2_id)) * &
                  rt_auxvars(ghosted_id)%pri_act_coef(abs(this%SiO2_id)) * &
                  waste_form%scaling_factor(i))
      enddo
    elseif (this%SiO2_id < 0) then   ! secondary species
      Q_local = 0.d0
      do i = 1,waste_form%region%num_cells
        ghosted_id = grid%nL2G(waste_form%region%cell_ids(i))
        Q_local = Q_local + &
                 (rt_auxvars(ghosted_id)%sec_molal(abs(this%SiO2_id)) * &
                  rt_auxvars(ghosted_id)%sec_act_coef(abs(this%SiO2_id)) * &
                  waste_form%scaling_factor(i))
      enddo
    endif
    call CalcParallelSUM(pm%option,waste_form%rank_list,Q_local,Q_global)
    this%Q = Q_global
  endif  ! Q -----------------------------------------------------------------

  ! kg-glass/m^2/sec
  this%dissolution_rate = this%k0 * (10.d0**(this%nu*this%pH)) * &
                          exp(-this%Ea/(8.314d0*avg_temp_global)) * &
                          (1.d0 - (this%Q/this%K)**(1.d0/this%v)) + this%k_long

  ! kg-glass/sec
  waste_form%eff_dissolution_rate = &
    this%dissolution_rate * &          ! kg-glass/m^2/sec
    this%specific_surface_area * &     ! m^2/kg glass
    this%matrix_density * &            ! kg-glass/m^3-glass
    waste_form%volume * &              ! m^3-glass
    waste_form%exposure_factor         ! [-]

end subroutine WFMechGlassDissolution

! ************************************************************************** !

subroutine WFMechDSNFDissolution(this,waste_form,pm,ierr)
  !
  ! Calculates the DSNF waste form dissolution rate
  !
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): base mechanism object
! waste_form (input/output): base waste form object
! pm (input/output): waste form process model
! ierr (input/output): [-] PETSc error integer
! -----------------------------------------
  class(wf_mechanism_dsnf_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr
! -----------------------------------------

  ierr = 0

  ! Because the DSNF dissolution rate is instantaneous, the amount of
  ! released isotopes gets updated directly in the solution vector after
  ! this routine is called, within PMWFSolve.
  ! Doing the direct update to the solution vector resolves the potential
  ! error that may occur if the next timestep size is different from the
  ! current timestep size, when the dissolution rate would have been
  ! calculated. This potential error is greatly reduced in magnitude for
  ! the other dissolution models, so we only do the direct update for DSNF.

  ! the entire waste form dissolves in the current timestep:
  this%frac_dissolution_rate = 1.d0 / pm%realization%option%tran_dt

  ! kg-matrix/sec
  waste_form%eff_dissolution_rate = &
    this%frac_dissolution_rate * &           ! 1/sec
    this%matrix_density * &                  ! kg matrix/m^3 matrix
    waste_form%volume * &                    ! m^3 matrix
    waste_form%exposure_factor               ! [-]

end subroutine WFMechDSNFDissolution

! ************************************************************************** !

subroutine WFMechWIPPDissolution(this,waste_form,pm,ierr)
  !
  ! Calculates the WIPP waste form dissolution rate
  !
  ! Author: Jenn Frederick
  ! Date: 012/8/2016

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): base mechanism object
! waste_form (input/output): base waste form object
! pm (input/output): waste form process model
! ierr (input/output): [-] PETSc error integer
! -----------------------------------------
  class(wf_mechanism_wipp_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr
! -----------------------------------------

  ! This subroutine is only a placeholder.
  ! The WIPP waste form mechanism is an extension of the DSNF waste form
  ! mechanism and does not have its own dissolution routine, however, this
  ! placeholder exists in case a different one should be implemented.

end subroutine WFMechWIPPDissolution

! ************************************************************************** !

subroutine WFMechFMDMDissolution(this,waste_form,pm,ierr)
  !
  ! Calculates the FMDM waste form dissolution rate using the FMDM model
  !
  ! Author: Jenn Frederick (with old code by Glenn Hammond)
  ! Date: 05/05/2016

  use Grid_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Option_module
  use Utility_module

  implicit none

 ! FMDM model:
 !===================================================================
  interface
    subroutine AMP_step ( burnup, sTme, temperature_C, conc, &
                          initialRun, fuelDisRate, Usource, success )
      real ( kind = 8), intent( in ) :: burnup
      real ( kind = 8), intent( in ) :: sTme
      real ( kind = 8), intent( in ) :: temperature_C
      real ( kind = 8), intent( inout ),  dimension (:,:) :: conc
      logical ( kind = 4), intent( in ) :: initialRun
      ! sum of fluxes of 3 uranium compounds (UO2,2+;UCO3,2+;UO2)
      ! units: g/m^2/yr where g = sum of uranium compound mass
      real ( kind = 8), intent(out) :: fuelDisRate
      ! flux of just the uranium from the 3 uranium compounds
      ! units: g/m^2/yr where g = uranium mass
      real ( kind = 8), intent(out) :: Usource
      integer ( kind = 4), intent(out) :: success
    end subroutine
  end interface
 !===================================================================

! INPUT ARGUMENTS:
! ================
! this (input/output): base mechanism object
! waste_form (input/output): base waste form object
! pm (input/output): waste form process model
! ierr (input/output): [-] PETSc error integer
! -----------------------------------------
  class(wf_mechanism_fmdm_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr
! -----------------------------------------

! LOCAL VARIABLES:
! ================
! grid: pointer to grid object
! rt_auxvars(:): pointer to reactive transport auxvar object, which stores
!    the total component concentration, and is indexed by the ghosted cell id
! i, k: [-] looping index integers
! icomp_fmdm: [-] FMDM species component number
! icomp_pflotran: [-] species component number mapped from FMDM to PFLOTRAN
! ghosted_id: [-] ghosted grid cell id
! avg_temp_local: [C] average temperature in the waste form region
! --------------------------------------------------------------
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscInt :: i, k
  PetscInt :: icomp_fmdm
  PetscInt :: icomp_pflotran
  PetscInt :: ghosted_id
! --------------------------------------------------------------

 ! FMDM model:
 !=======================================================
  integer ( kind = 4) :: success
  logical ( kind = 4) :: initialRun
  PetscReal :: Usource
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(option_type), pointer :: option
 !========================================================

  grid => pm%realization%patch%grid
  rt_auxvars => pm%realization%patch%aux%RT%auxvars
  global_auxvars => pm%realization%patch%aux%Global%auxvars
  option => pm%realization%option

  ierr = 0

  do k = 1,waste_form%region%num_cells
    ghosted_id = grid%nL2G(waste_form%region%cell_ids(k))
    ! overwrite the components in mapping_pflotran array
    do i = 1, size(this%mapping_fmdm)
      icomp_fmdm = this%mapping_fmdm(i)
      icomp_pflotran = this%mapping_fmdm_to_pflotran(icomp_fmdm)
      ! original -- set ALL cells to the concentration values
      this%concentration(icomp_fmdm,:) = &
        rt_auxvars(ghosted_id)%total(icomp_pflotran,LIQUID_PHASE)
      ! It seems like only the boundary should be set by PFLOTRAN,
      ! i.e.
      !geh: the boundary should be at beginning of the array, not the end
      !this%concentration(icomp_fmdm,this%num_grid_cells_in_waste_form) = &
      !  rt_auxvars(ghosted_id)%total(icomp_pflotran,LIQUID_PHASE)
    enddo
  enddo

  ! convert total component concentration from mol/L to mol/m3 (*1.d3)
  this%concentration = this%concentration*1.d3

  if (.not. Equal(waste_form%volume,waste_form%init_volume)) then
    initialRun = PETSC_FALSE
  else
    initialRun = PETSC_TRUE
  endif

#ifdef FMDM_MODEL
 ! FMDM model calculates this%dissolution_rate and Usource [g/m^2/yr]:
 !====================================================================
  time = option%time

  avg_temp_local = 0.d0
  do i = 1,waste_form%region%num_cells
    ghosted_id = grid%nL2G(waste_form%region%cell_ids(i))
    avg_temp_local = avg_temp_local + &  ! Celsius
               global_auxvars(ghosted_id)%temp * waste_form%scaling_factor(i)
  enddo
  call CalcParallelSUM(option,waste_form%rank_list,avg_temp_local, &
                       avg_temp_global)
  call AMP_step(this%burnup, time, avg_temp_global, this%concentration, &
                initialRun, this%dissolution_rate, Usource, success)
  !write(*,*) this%dissolution_rate
  ! convert total component concentration from mol/m3 back to mol/L (/1.d3)
  this%concentration = this%concentration/1.d3
  ! convert this%dissolution_rate from fmdm to pflotran units:
  ! g/m^2/yr => kg/m^2/sec
  this%dissolution_rate = this%dissolution_rate / &
                          (1000.d0*24.d0*3600.d0*DAYS_PER_YEAR)
  Usource = Usource / (1000.d0*24.d0*3600.d0*DAYS_PER_YEAR)
 !====================================================================
#else
  ! if no FMDM model, use the burnup as this%dissolution_rate:
  ! if no FMDM model, the units of burnup should already be kg-matrix/m^2/sec:
  success = 1
  this%dissolution_rate = this%burnup
  Usource = this%burnup
#endif

  if (success == 0) then
    ierr = 1
    return
  endif

  !==================
  this%frac_dissolution_rate = &    ! 1/sec
    this%dissolution_rate * &       ! kg-matrix/m^2/sec
    this%specific_surface_area      ! m^2/kg-matrix
  !==================

  ! kg-matrix / sec
  waste_form%eff_dissolution_rate = &
     this%dissolution_rate * &         ! kg-matrix/m^2/sec
     this%specific_surface_area * &    ! m^2/kg-matrix
     this%matrix_density * &           ! kg-matrix/m^3-matrix
     waste_form%volume * &             ! m^3-matrix
     waste_form%exposure_factor        ! [-]

end subroutine WFMechFMDMDissolution

! ************************************************************************** !

subroutine WFMechFMDMSurrogateDissolution(this,waste_form,pm,ierr)
  !
  ! Calculates the FMDM waste form dissolution rate using a
  ! single-layer feed-forward artifical neural network
  ! SURROGATE APPROXIMATION of the FMDM model
  !
  ! Author: Tom Seidl (with old code by Jenn Frederick and Glenn Hammond)
  ! Date: 11/12/2019

  use Grid_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Option_module
  use Utility_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): base mechanism object
! waste_form (input/output): base waste form object
! pm (input/output): waste form process model
! ierr (input/output): [-] PETSc error integer
! -----------------------------------------
  class(wf_mechanism_fmdm_surrogate_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr
! -----------------------------------------

! LOCAL VARIABLES:
! ================
! grid: pointer to grid object
! rt_auxvars(:): pointer to reactive transport auxvar object, which stores
!    the total component concentration, and is indexed by the ghosted cell id
! i, k: [-] looping index integers
! icomp_fmdm: [-] FMDM species component number
! icomp_pflotran: [-] species component number mapped from FMDM to PFLOTRAN
! ghosted_id: [-] ghosted grid cell id
! avg_temp_local: [C] average temperature in the waste form region
! --------------------------------------------------------------
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscInt :: i, k
  PetscInt :: icomp_surrfmdm
  PetscInt :: icomp_pflotran
  PetscInt :: ghosted_id
  PetscReal :: avg_temp_local
! --------------------------------------------------------------

 ! FMDM surrogate model:
 !=======================================================
  PetscReal :: time
  PetscReal :: avg_temp_global
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(option_type), pointer :: option
 !========================================================

  grid => pm%realization%patch%grid
  rt_auxvars => pm%realization%patch%aux%RT%auxvars
  global_auxvars => pm%realization%patch%aux%Global%auxvars
  option => pm%realization%option

  do k = 1,waste_form%region%num_cells
    ghosted_id = grid%nL2G(waste_form%region%cell_ids(k))
    ! overwrite the components in mapping_pflotran array
    do i = 1, size(this%mapping_surrfmdm)
      icomp_surrfmdm = this%mapping_surrfmdm(i)
      icomp_pflotran = this%mapping_surrfmdm_to_pflotran(icomp_surrfmdm)
      this%concentration(icomp_surrfmdm) = &
        rt_auxvars(ghosted_id)%total(icomp_pflotran,LIQUID_PHASE)
    enddo
  enddo

  ! convert total component concentration from mol/L to mol/m3 (*1.d3)
  this%concentration = this%concentration*1.d3

 ! FMDM surrogate model calculates this%dissolution_rate [g/m^2/yr]:
 !====================================================================
  time = option%time

  avg_temp_local = 0.d0
  do i = 1,waste_form%region%num_cells
    ghosted_id = grid%nL2G(waste_form%region%cell_ids(i))
    avg_temp_local = avg_temp_local + &  ! Celsius
               global_auxvars(ghosted_id)%temp * waste_form%scaling_factor(i)
  enddo
  call CalcParallelSUM(option,waste_form%rank_list,avg_temp_local, &
                       avg_temp_global)

  if (FMDM_surrogate_knnr) then
    call KnnrQuery(this, time, avg_temp_global)
  else
    call AMP_ann_surrogate_step(this, time, avg_temp_global)
  endif

  ! convert total component concentration from mol/m3 back to mol/L (/1.d3)
  this%concentration = this%concentration/1.d3
  ! convert this%dissolution_rate from fmdm to pflotran units:
  ! g/m^2/yr => kg/m^2/sec
  this%dissolution_rate = this%dissolution_rate / &
                          (1000.d0*24.d0*3600.d0*DAYS_PER_YEAR)

  ierr = 0
  !==================
  this%frac_dissolution_rate = &    ! 1/sec
    this%dissolution_rate * &       ! kg-matrix/m^2/sec
    this%specific_surface_area      ! m^2/kg-matrix
  !==================

  ! kg-matrix / sec
  waste_form%eff_dissolution_rate = &
     this%dissolution_rate * &         ! kg-matrix/m^2/sec
     this%specific_surface_area * &    ! m^2/kg-matrix
     this%matrix_density * &           ! kg-matrix/m^3-matrix
     waste_form%volume * &             ! m^3-matrix
     waste_form%exposure_factor        ! [-]

end subroutine WFMechFMDMSurrogateDissolution

! ************************************************************************** !

subroutine WFMechCustomDissolution(this,waste_form,pm,ierr)
  !
  ! Calculates the "custom" waste form dissolution rate
  !
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  use Grid_module
  use Global_Aux_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): base mechanism object
! waste_form (input/output): base waste form object
! pm (input/output): waste form process model
! ierr (input/output): [-] PETSc error integer
! -----------------------------------------
  class(wf_mechanism_custom_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr
! -----------------------------------------

  ! Note: Units for dissolution rates have already been converted to
  ! internal units within the PMWFRead routine.

  ierr = 0

  if (Uninitialized(this%frac_dissolution_rate)) then
    ! kg-matrix / sec
    waste_form%eff_dissolution_rate = &
       this%dissolution_rate * &         ! kg-matrix/m^2/sec
       this%specific_surface_area * &    ! m^2/kg-matrix
       this%matrix_density * &           ! kg-matrix/m^3-matrix
       waste_form%volume * &             ! m^3-matrix
       waste_form%exposure_factor        ! [-]
  else if (this%frac_diss_vol_init) then
    ! kg-matrix / sec
    waste_form%eff_dissolution_rate = &
       this%frac_dissolution_rate * &     ! [-]/sec
       this%matrix_density * &            ! kg-matrix/m^3-matrix
       waste_form%init_volume * &         ! m^3 matrix
       waste_form%exposure_factor         ! [-]
  else
    ! kg-matrix / sec
    waste_form%eff_dissolution_rate = &
       this%frac_dissolution_rate * &     ! [-]/sec
       this%matrix_density * &            ! kg-matrix/m^3-matrix
       waste_form%volume * &              ! m^3 matrix
       waste_form%exposure_factor         ! [-]
  endif

end subroutine WFMechCustomDissolution

! ************************************************************************** !

subroutine PMWFFinalizeTimestep(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process mode object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

end subroutine PMWFFinalizeTimestep

! ************************************************************************** !

subroutine PMWFUpdateSolution(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process mode object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

end subroutine PMWFUpdateSolution

! ************************************************************************** !

recursive subroutine PMWFFinalizeRun(this)
  !
  ! Finalizes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process mode object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMWFFinalizeRun

! ************************************************************************** !

subroutine PMWFOutput(this)
  !
  ! Sets up output for a waste form process model
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Option_module
  use Output_Aux_module
  use Global_Aux_module
  use Grid_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process mode object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

! LOCAL VARIABLES:
! ================
! option: pointer to option object
! output_option: pointer to output option object
! cur_waste_form: pointer to curent waste form object
! grid: pointer to grid object
! filename: filename string
! fid: [-] file id number
! i: [-] looping index integer
! ------------------------------------------------------
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  class(waste_form_base_type), pointer :: cur_waste_form
  type(grid_type), pointer :: grid
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: i
! ------------------------------------------------------

  if (.not.associated(this%waste_form_list)) return

100 format(100es18.8)
101 format(1I6.1)

  option => this%realization%option
  output_option => this%realization%output_option
  grid => this%realization%patch%grid

  fid = 86
  filename = PMWFOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")

  ! this time is set at the end of the reactive transport step
  write(fid,100,advance="no") option%time / output_option%tconv

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    write(fid,101,advance="no") cur_waste_form%id
    do i = 1, cur_waste_form%mechanism%num_species
      write(fid,100,advance="no") cur_waste_form%cumulative_mass(i), &
                                  cur_waste_form%instantaneous_mass_rate(i), &
                                  cur_waste_form%rad_mass_fraction(i)
    enddo
    write(fid,100,advance="no") cur_waste_form%eff_dissolution_rate, &
                                cur_waste_form%volume, &
                                cur_waste_form%eff_canister_vit_rate, &
                                cur_waste_form%canister_vitality*100.d0

    if (associated(cur_waste_form%spacer_mechanism)) then
      write(fid,100,advance="no") cur_waste_form%spacer_vitality_rate, &
                                  cur_waste_form%spacer_vitality
    endif

    cur_waste_form => cur_waste_form%next
  enddo
  close(fid)

end subroutine PMWFOutput

! ************************************************************************** !

function PMWFOutputFilename(option)
  !
  ! Generates filename for waste form output
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Option_module

  implicit none

! INPUT ARGUMENTS:
! ================
! option (input): pointer to option object
! ------------------------------------
  type(option_type), pointer :: option
! ------------------------------------

! LOCAL VARIABLES:
! ================
! PMWFOutputFilename (output): filename string
! word: temporary string
! ----------------------------------------------------
  character(len=MAXSTRINGLENGTH) :: PMWFOutputFilename
  character(len=MAXWORDLENGTH) :: word
! ----------------------------------------------------

  write(word,'(i6)') option%myrank
  PMWFOutputFilename = trim(option%global_prefix) // &
                       trim(option%group_prefix) // &
                       '-' // trim(adjustl(word)) // '.wf'

end function PMWFOutputFilename

! ************************************************************************** !

subroutine PMWFOutputHeader(this)
  !
  ! Writes header for waste form output file
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Output_Aux_module
  use Grid_module
  use Utility_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

! LOCAL VARIABLES:
! ================
! output_option: pointer to output option object
! grid: pointer to grid object
! cur_waste_form: pointer to current waste form object
! cell_string: cell string
! x_string, y_string, z_string: coordinate strings
! units_string: units string
! variable_string: variable string
! filename: filename string
! fid: [-] file id number
! icolumn: [-] column number
! i: [-] looping index integer
! exist: Boolean helper to check is file exists
! -------------------------------------------------------------
  type(output_option_type), pointer :: output_option
  type(grid_type), pointer :: grid
  class(waste_form_base_type), pointer :: cur_waste_form
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: icolumn, i
  PetscBool :: exist
! -------------------------------------------------------------

  if (.not.associated(this%waste_form_list)) return

  output_option => this%realization%output_option
  grid => this%realization%patch%grid

  fid = 86
  filename = PMWFOutputFilename(this%option)
  exist = FileExists(trim(filename))
  if (this%option%restart_flag .and. exist) return
  open(unit=fid,file=filename,action="write",status="replace")

  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif

  write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    if (initialized(cur_waste_form%coordinate%z)) then
      ! cell natural id
      write(cell_string,*) grid%nG2A(grid%nL2G(cur_waste_form%region%cell_ids(1)))
      cell_string = ' (' // trim(adjustl(cell_string)) // ')'
      ! coordinate of waste form
      x_string = BestFloat(cur_waste_form%coordinate%x,1.d4,1.d-2)
      y_string = BestFloat(cur_waste_form%coordinate%y,1.d4,1.d-2)
      z_string = BestFloat(cur_waste_form%coordinate%z,1.d4,1.d-2)
      cell_string = trim(cell_string) // &
               ' (' // trim(adjustl(x_string)) // &
               ' ' // trim(adjustl(y_string)) // &
               ' ' // trim(adjustl(z_string)) // ')'
    else
      cell_string = trim(cur_waste_form%region_name)
    endif
    variable_string = 'WF ID#'
    units_string = ''
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    do i = 1, cur_waste_form%mechanism%num_species
      variable_string = trim(cur_waste_form%mechanism%rad_species_list(i)%name) &
                        // ' Cum. Release'
      ! cumulative
      units_string = 'mol'
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      variable_string = trim(cur_waste_form%mechanism%rad_species_list(i)%name) &
                        // ' Release Rate'
      ! instantaneous
      units_string = 'mol/s' !// trim(adjustl(output_option%tunit))
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      variable_string = trim(cur_waste_form%mechanism%rad_species_list(i)%name) &
                        // ' Mass Frac.'
      units_string = 'g-rad/g-matrix'
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
    enddo
    variable_string = 'WF Dissolution Rate'
    units_string = 'kg/s'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'WF Volume'
    units_string = 'm^3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Canister Vitality Deg. Rate'
    units_string = '1/yr'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Canister Vitality'
    units_string = '%'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)

    if (associated(cur_waste_form%spacer_mechanism)) then
      variable_string = 'Spacer Vitality Deg. Rate'
      units_string = '1/yr'
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)

      variable_string = 'Spacer Vitality'
      units_string = '%'
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
    endif

    cur_waste_form => cur_waste_form%next
  enddo

  close(fid)

end subroutine PMWFOutputHeader

! ***************************************************************************** !

subroutine PMWFCheckpointHDF5(this,pm_grp_id)
  !
  ! Checkpoints data associated with the waste form process model
  ! into the "canister_properties" dataset for a given wf pm:
  ! canister vitality (1), canister volume (2), breach time (3),
  ! spacer grid vitality (4),
  ! radionuclide mass fraction (5:2:end-1),
  ! cumulative mass released (6:2:end) .

  !
  ! Author: Michael Nole
  ! Date: 09/21/18
  !

  use Option_module
  use Realization_Subsurface_class
  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec

  implicit none

  ! Input Arguments
  class(pm_waste_form_type) :: this

  integer(HID_T) :: pm_grp_id

  ! Local Variables
  IS :: is
  VecScatter :: scatter_ctx
  Vec :: global_wf_vec, local_wf_vec

  character(len=MAXSTRINGLENGTH) :: dataset_name

  PetscErrorCode :: ierr
  PetscInt :: local_stride, n_wf_local, n_wf_global, n_check_vars, &
              n_vecs, local_stride_tmp, i, j, num_species, stride
  PetscInt, allocatable :: indices(:), int_array(:)
  PetscReal, allocatable :: check_vars(:)

  class(waste_form_base_type), pointer :: cur_waste_form

  cur_waste_form => this%waste_form_list

  local_stride=0
  local_stride_tmp=0
  n_wf_local=0
  n_wf_global=0

  n_check_vars=4 !number of scalar wf checkpoint variables
  n_vecs=2 !number of vector wf checkpoint variables (by species)

  do
    if (.not.associated(cur_waste_form)) exit
    n_wf_local=n_wf_local+1
    local_stride_tmp=local_stride_tmp+n_check_vars+&
                     cur_waste_form%mechanism%num_species*n_vecs
    cur_waste_form => cur_waste_form%next
    if (local_stride_tmp>local_stride) then
      local_stride=local_stride_tmp
    endif
    local_stride_tmp=0
  enddo

  cur_waste_form => this%waste_form_list
  allocate(int_array(n_wf_local))
  i=1
  do
    if (.not.associated(cur_waste_form)) exit
    int_array(i)=cur_waste_form%id-1
    i=i+1
    cur_waste_form => cur_waste_form%next
  enddo


  !Gather relevant information from all processes
  call MPI_Allreduce(local_stride,stride,ONE_INTEGER_MPI,MPI_INTEGER,MPI_MAX, &
                     this%option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(n_wf_local,n_wf_global,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_SUM,this%option%mycomm,ierr);CHKERRQ(ierr)

  !Create MPI vector and sequential vector for mapping
  call VecCreateMPI(this%option%mycomm,n_wf_local*stride,n_wf_global*stride, &
                    global_wf_vec,ierr);CHKERRQ(ierr)

  call VecCreateSeq(PETSC_COMM_SELF,n_wf_local*stride,local_wf_vec, &
                    ierr);CHKERRQ(ierr)

  call VecSetBlockSize(global_wf_vec,stride,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(local_wf_vec,stride,ierr);CHKERRQ(ierr)


  allocate(check_vars(stride))
  allocate(indices(stride))

  !Collect data for checkpointing
  j=1
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    num_species=cur_waste_form%mechanism%num_species
    check_vars(1)=cur_waste_form%canister_vitality
    check_vars(2)=cur_waste_form%volume
    check_vars(3)=cur_waste_form%breach_time
    check_vars(4)=cur_waste_form%spacer_vitality

    do i = 1,num_species
      check_vars(n_vecs*i-1+n_check_vars)=cur_waste_form%rad_mass_fraction(i)
      check_vars(n_vecs*i-1+n_check_vars+1)=cur_waste_form%cumulative_mass(i)
    enddo
    i=n_vecs*num_species+n_check_vars+1
    do
      if (i>stride) exit
      check_vars(i)=-9999
      i=i+1
    enddo

    do i = 1,stride
      indices(i)=(j-1)*stride+i-1
    enddo
    j=j+1

    call VecSetValues(local_wf_vec,stride,indices,check_vars,INSERT_VALUES, &
                      ierr);CHKERRQ(ierr)
    cur_waste_form => cur_waste_form%next

  enddo

  !Create map and add values from the sequential vector to the global
  call ISCreateBlock(this%option%mycomm,stride,n_wf_local,int_array, &
                     PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)

  call VecScatterCreate(local_wf_vec,PETSC_NULL_IS,global_wf_vec,is, &
                        scatter_ctx,ierr);CHKERRQ(ierr)

  call VecScatterBegin(scatter_ctx,local_wf_vec,global_wf_vec,INSERT_VALUES, &
                       SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx,local_wf_vec,global_wf_vec,INSERT_VALUES, &
                     SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  dataset_name='canister_properties'

  !Write the checkpoint file
  call HDF5WriteDataSetFromVec(dataset_name, this%option, global_wf_vec,&
           pm_grp_id, H5T_NATIVE_DOUBLE)

  call VecScatterDestroy(scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  call VecDestroy(global_wf_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(local_wf_vec,ierr);CHKERRQ(ierr)

end subroutine PMWFCheckpointHDF5

! ***************************************************************************** !


subroutine PMWFRestartHDF5(this,pm_grp_id)
  !
  ! Restarts data associated with waste form process model
  ! from the "canister_properties" dataset for a given wf pm:
  ! canister vitality (1), canister volume (2), breach time (3),
  ! spacer grid vitality (4),
  ! radionuclide mass fraction (5:2:end-1),
  ! cumulative mass released (6:2:end) .
  !
  ! Author: Michael Nole
  ! Date: 10/03/18

  use Option_module
  use Realization_Subsurface_class
  use hdf5
  use HDF5_module, only : HDF5ReadDataSetInVec

  implicit none

  ! Input Arguments
  class(pm_waste_form_type) :: this

  integer(HID_T) :: pm_grp_id

  ! Local Variables
  IS :: is
  VecScatter :: scatter_ctx
  character(len=MAXSTRINGLENGTH) :: dataset_name

  PetscErrorCode :: ierr
  PetscInt :: local_stride, n_wf_local, n_wf_global, n_check_vars, &
              n_vecs, local_stride_tmp, i, j, num_species, stride
  PetscInt, allocatable :: int_array(:)
  PetscReal, pointer :: local_wf_array(:)

  class(waste_form_base_type), pointer :: cur_waste_form

  Vec :: global_wf_vec, local_wf_vec

  cur_waste_form => this%waste_form_list

  local_stride=0
  local_stride_tmp=0
  n_wf_local=0
  n_wf_global=0

  n_check_vars=4 !number of scalar wf checkpoint variables
  n_vecs=2 !number of vector wf checkpoint variables (by species)

  do
    if (.not.associated(cur_waste_form)) exit
    n_wf_local=n_wf_local+1
    local_stride_tmp=local_stride_tmp+n_check_vars+&
                     cur_waste_form%mechanism%num_species*n_vecs
    cur_waste_form => cur_waste_form%next
    if (local_stride_tmp>local_stride) then
      local_stride=local_stride_tmp
    endif
    local_stride_tmp=0
  enddo

  cur_waste_form => this%waste_form_list
  allocate(int_array(n_wf_local))
  i=1
  do
    if (.not.associated(cur_waste_form)) exit
    int_array(i)=cur_waste_form%id-1
    i=i+1
    cur_waste_form => cur_waste_form%next
  enddo

  !Gather relevant information
  call MPI_Allreduce(local_stride,stride,ONE_INTEGER_MPI,MPI_INTEGER,MPI_MAX, &
                     this%option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(n_wf_local,n_wf_global,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_SUM,this%option%mycomm,ierr);CHKERRQ(ierr)

  !Create MPI vector into which HDF5 will read, and sequential vector
  !for wf information on a given process.
  call VecCreateMPI(this%option%mycomm,n_wf_local*stride,n_wf_global*stride, &
                    global_wf_vec,ierr);CHKERRQ(ierr)

  call VecCreateSeq(PETSC_COMM_SELF,n_wf_local*stride,local_wf_vec, &
                    ierr);CHKERRQ(ierr)

  call VecSetBlockSize(global_wf_vec,stride,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(local_wf_vec,stride,ierr);CHKERRQ(ierr)

  !Read the data
  dataset_name = 'canister_properties'
  call HDF5ReadDataSetInVec(dataset_name, this%option, global_wf_vec, &
                             pm_grp_id, H5T_NATIVE_DOUBLE)

  !Create map between MPI and sequential vectors
  call ISCreateBlock(this%option%mycomm,stride,n_wf_local,int_array, &
                     PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)

  call VecScatterCreate(global_wf_vec,is,local_wf_vec,PETSC_NULL_IS, &
                        scatter_ctx,ierr);CHKERRQ(ierr)

  !Get the data from the MPI vector
  call VecScatterBegin(scatter_ctx,global_wf_vec,local_wf_vec,INSERT_VALUES, &
                       SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx,global_wf_vec,local_wf_vec,INSERT_VALUES, &
                     SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  !Convert the data to a Fortran array
  call VecGetArrayF90(local_wf_vec,local_wf_array,ierr);CHKERRQ(ierr)

  !Assign checkpointed waste form attribute values
  i=1
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    num_species=cur_waste_form%mechanism%num_species

    allocate(cur_waste_form%rad_mass_fraction(num_species))
    allocate(cur_waste_form%cumulative_mass(num_species))

    cur_waste_form%canister_vitality=local_wf_array(i)
    cur_waste_form%volume=local_wf_array(i+1)
    cur_waste_form%breach_time=local_wf_array(i+2)
    cur_waste_form%spacer_vitality=local_wf_array(i+3)
    if (cur_waste_form%breach_time<0) then
      cur_waste_form%breached=PETSC_FALSE
    else
      cur_waste_form%breached=PETSC_TRUE
    endif

    do j = 1,num_species
      cur_waste_form%rad_mass_fraction(j)=local_wf_array(2*(j-1)+i+n_check_vars)
      cur_waste_form%cumulative_mass(j)=local_wf_array(2*j-1+i+n_check_vars)
    enddo
    cur_waste_form => cur_waste_form%next
    i=i+stride
  enddo

  call VecRestoreArrayF90(local_wf_vec,local_wf_array,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  call VecDestroy(global_wf_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(local_wf_vec,ierr);CHKERRQ(ierr)

end subroutine PMWFRestartHDF5

! ************************************************************************** !

subroutine PMWFCheckpointBinary(this, viewer)
  !
  ! Checkpoints data associated with the waste form process model
  ! into a checkpiont binary file for a given wf pm:
  ! canister vitality (1), canister volume (2), breach time (3),
  ! spacer grid vitality (4),
  ! radionuclide mass fraction (5:2:end-1),
  ! cumulative mass released (6:2:end)

  !
  ! Author: Michael Nole
  ! Date: 10/09/18
  !

  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Field_module
  use Discretization_module
  use Grid_module

  implicit none

  !Input Arguments
  PetscViewer :: viewer
  class(pm_waste_form_type) :: this

  ! Local Variables

  IS :: is
  VecScatter :: scatter_ctx
  Vec :: global_wf_vec, local_wf_vec

  PetscErrorCode :: ierr
  PetscInt :: local_stride, n_wf_local, n_wf_global, n_check_vars, &
              n_vecs, local_stride_tmp, i, j, num_species, stride
  PetscInt, allocatable :: indices(:), int_array(:)
  PetscReal, allocatable :: check_vars(:)

  class(waste_form_base_type), pointer :: cur_waste_form

  cur_waste_form => this%waste_form_list

  local_stride=0
  local_stride_tmp=0
  n_wf_local=0
  n_wf_global=0

  n_check_vars=4 !number of scalar wf checkpoint variables
  n_vecs=2 !number of vector wf checkpoint variables (by species)

  do
    if (.not.associated(cur_waste_form)) exit
    n_wf_local=n_wf_local+1
    local_stride_tmp=local_stride_tmp+n_check_vars+&
                     cur_waste_form%mechanism%num_species*n_vecs
    cur_waste_form => cur_waste_form%next
    if (local_stride_tmp>local_stride) then
      local_stride=local_stride_tmp
    endif
    local_stride_tmp=0
  enddo

  cur_waste_form => this%waste_form_list
  allocate(int_array(n_wf_local))
  i=1
  do
    if (.not.associated(cur_waste_form)) exit
    int_array(i)=cur_waste_form%id-1
    i=i+1
    cur_waste_form => cur_waste_form%next
  enddo


  !Gather relevant information from all processes
  call MPI_Allreduce(local_stride,stride,ONE_INTEGER_MPI,MPI_INTEGER,MPI_MAX, &
                     this%option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(n_wf_local,n_wf_global,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_SUM,this%option%mycomm,ierr);CHKERRQ(ierr)

  !Create MPI vector and sequential vector for mapping
  call VecCreateMPI(this%option%mycomm,n_wf_local*stride,n_wf_global*stride, &
                    global_wf_vec,ierr);CHKERRQ(ierr)

  call VecCreateSeq(PETSC_COMM_SELF,n_wf_local*stride,local_wf_vec, &
                    ierr);CHKERRQ(ierr)

  call VecSetBlockSize(global_wf_vec,stride,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(local_wf_vec,stride,ierr);CHKERRQ(ierr)


  allocate(check_vars(stride))
  allocate(indices(stride))

  !Collect data for checkpointing
  j=1
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    num_species=cur_waste_form%mechanism%num_species
    check_vars(1)=cur_waste_form%canister_vitality
    check_vars(2)=cur_waste_form%volume
    check_vars(3)=cur_waste_form%breach_time
    check_vars(4)=cur_waste_form%spacer_vitality

    do i = 1,num_species
      check_vars(n_vecs*i-1+n_check_vars)=cur_waste_form%rad_mass_fraction(i)
      check_vars(n_vecs*i-1+n_check_vars+1)=cur_waste_form%cumulative_mass(i)
    enddo
    i=n_vecs*num_species+n_check_vars+1
    do
      if (i>stride) exit
      check_vars(i)=-9999
      i=i+1
    enddo

    do i = 1,stride
      indices(i)=(j-1)*stride+i-1
    enddo
    j=j+1

    call VecSetValues(local_wf_vec,stride,indices,check_vars,INSERT_VALUES, &
                      ierr);CHKERRQ(ierr)
    cur_waste_form => cur_waste_form%next

  enddo

  !Create map and add values from the sequential vector to the global
  call ISCreateBlock(this%option%mycomm,stride,n_wf_local,int_array, &
                     PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)

  call VecScatterCreate(local_wf_vec,PETSC_NULL_IS,global_wf_vec,is, &
                        scatter_ctx,ierr);CHKERRQ(ierr)

  call VecScatterBegin(scatter_ctx,local_wf_vec,global_wf_vec,INSERT_VALUES, &
                       SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx,local_wf_vec,global_wf_vec,INSERT_VALUES, &
                     SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  !Write the checkpoint file

  call VecView(global_wf_vec,viewer,ierr);CHKERRQ(ierr)

  call VecScatterDestroy(scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  call VecDestroy(global_wf_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(local_wf_vec,ierr);CHKERRQ(ierr)



end subroutine PMWFCheckpointBinary

! ***************************************************************************** !

subroutine PMWFRestartBinary(this, viewer)

  !
  ! Restarts data associated with waste form process model
  ! from a checkpoint binary file for a given wf pm:
  ! canister vitality (1), canister volume (2), breach time (3),
  ! spacer grid vitality (4),
  ! radionuclide mass fraction (5:2:end-1),
  ! cumulative mass released (6:2:end) .
  !
  ! Author: Michael Nole
  ! Date: 10/09/18

  use Option_module
  use Realization_Subsurface_class
  use hdf5
  use HDF5_module, only : HDF5ReadDataSetInVec

  implicit none

  ! Input Arguments
  PetscViewer :: viewer
  class(pm_waste_form_type) :: this

  ! Local Variables
  IS :: is
  VecScatter :: scatter_ctx

  PetscErrorCode :: ierr
  PetscInt :: local_stride, n_wf_local, n_wf_global, n_check_vars, &
              n_vecs, local_stride_tmp, i, j, num_species, stride
  PetscInt, allocatable :: int_array(:)
  PetscReal, pointer :: local_wf_array(:)

  class(waste_form_base_type), pointer :: cur_waste_form

  Vec :: global_wf_vec, local_wf_vec

  cur_waste_form => this%waste_form_list

  local_stride=0
  local_stride_tmp=0
  n_wf_local=0
  n_wf_global=0

  n_check_vars=4 !number of scalar wf checkpoint variables
  n_vecs=2 !number of vector wf checkpoint variables (by species)

  do
    if (.not.associated(cur_waste_form)) exit
    n_wf_local=n_wf_local+1
    local_stride_tmp=local_stride_tmp+n_check_vars+&
                     cur_waste_form%mechanism%num_species*n_vecs
    cur_waste_form => cur_waste_form%next
    if (local_stride_tmp>local_stride) then
      local_stride=local_stride_tmp
    endif
    local_stride_tmp=0
  enddo

  cur_waste_form => this%waste_form_list
  allocate(int_array(n_wf_local))
  i=1
  do
    if (.not.associated(cur_waste_form)) exit
    int_array(i)=cur_waste_form%id-1
    i=i+1
    cur_waste_form => cur_waste_form%next
  enddo

  !Gather relevant information
  call MPI_Allreduce(local_stride,stride,ONE_INTEGER_MPI,MPI_INTEGER,MPI_MAX, &
                     this%option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(n_wf_local,n_wf_global,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_SUM,this%option%mycomm,ierr);CHKERRQ(ierr)

  !Create MPI vector into which HDF5 will read, and sequential vector
  !for wf information on a given process.
  call VecCreateMPI(this%option%mycomm,n_wf_local*stride,n_wf_global*stride, &
                    global_wf_vec,ierr);CHKERRQ(ierr)

  call VecCreateSeq(PETSC_COMM_SELF,n_wf_local*stride,local_wf_vec, &
                    ierr);CHKERRQ(ierr)

  call VecSetBlockSize(global_wf_vec,stride,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(local_wf_vec,stride,ierr);CHKERRQ(ierr)

  !Read the data
  call VecLoad(global_wf_vec,viewer,ierr);CHKERRQ(ierr)

  !Create map between MPI and sequential vectors
  call ISCreateBlock(this%option%mycomm,stride,n_wf_local,int_array, &
                     PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)

  call VecScatterCreate(global_wf_vec,is,local_wf_vec,PETSC_NULL_IS, &
                        scatter_ctx,ierr);CHKERRQ(ierr)

  !Get the data from the MPI vector
  call VecScatterBegin(scatter_ctx,global_wf_vec,local_wf_vec,INSERT_VALUES, &
                       SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx,global_wf_vec,local_wf_vec,INSERT_VALUES, &
                     SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  !Convert the data to a Fortran array
  call VecGetArrayF90(local_wf_vec,local_wf_array,ierr);CHKERRQ(ierr)

  !Assign checkpointed waste form attribute values
  i=1
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    num_species=cur_waste_form%mechanism%num_species

    allocate(cur_waste_form%rad_mass_fraction(num_species))
    allocate(cur_waste_form%cumulative_mass(num_species))

    cur_waste_form%canister_vitality=local_wf_array(i)
    cur_waste_form%volume=local_wf_array(i+1)
    cur_waste_form%breach_time=local_wf_array(i+2)
    cur_waste_form%spacer_vitality=local_wf_array(i+3)
    if (cur_waste_form%breach_time<0) then
      cur_waste_form%breached=PETSC_FALSE
    else
      cur_waste_form%breached=PETSC_TRUE
    endif

    do j = 1,num_species
      cur_waste_form%rad_mass_fraction(j)=local_wf_array(2*(j-1)+i+n_check_vars)
      cur_waste_form%cumulative_mass(j)=local_wf_array(2*j-1+i+n_check_vars)
    enddo
    cur_waste_form => cur_waste_form%next
    i=i+stride
  enddo

  call VecRestoreArrayF90(local_wf_vec,local_wf_array,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  call VecDestroy(global_wf_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(local_wf_vec,ierr);CHKERRQ(ierr)


end subroutine PMWFRestartBinary

! ***************************************************************************** !

subroutine PMWFInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  !
  ! Modified by Alex Salazar III
  ! Date: 05/13/2021

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

! LOCAL VARIABLES:
! ================
! id: [-] file id number
! --------------
  PetscInt :: id
! --------------

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  if (associated(this%mechanism_list)) then
    call MechanismInputRecord(this%mechanism_list);
  endif

  if (associated(this%spacer_mech_list)) then
    call SpacerMechInputRecord(this%spacer_mech_list)
  endif

  if (associated(this%criticality_mediator)) then
    call CritMechInputRecord(this%criticality_mediator);
  endif


end subroutine PMWFInputRecord

! ************************************************************************** !

subroutine WasteFormInputRecord(this)
  !
  ! Writes waste form information to the input record file.
  !
  ! Author: Alex Salazar, SNL
  ! Date: 10/22/2020
  !
  implicit none
  ! INPUT ARGUMENTS:
  ! ================
  class(waste_form_base_type), pointer :: this
  ! ---------------------------------
  ! LOCAL VARIABLES:
  ! ================
  class(waste_form_base_type), pointer :: cur_waste_form
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
  ! ---------------------------------

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'WASTE FORMS'

  cur_waste_form => this
  do
    if (.not. associated(cur_waste_form)) exit

    if (len_trim(adjustl(cur_waste_form%region_name)) > 0) then
      write(id,'(a29)',advance='no') 'region: '
      write(id,'(a)') cur_waste_form%region_name
    endif

    if (Initialized(cur_waste_form%coordinate%x)) then
      write(id,'(a29)',advance='no') 'x-coordinate: '
      write(word,'(es12.5)') cur_waste_form%coordinate%x
      write(id,'(a)') trim(adjustl(word))
    endif

    if (Initialized(cur_waste_form%coordinate%y)) then
      write(id,'(a29)',advance='no') 'y-coordinate: '
      write(word,'(es12.5)') cur_waste_form%coordinate%y
      write(id,'(a)') trim(adjustl(word))
    endif

    if (Initialized(cur_waste_form%coordinate%z)) then
      write(id,'(a29)',advance='no') 'z-coordinate: '
      write(word,'(es12.5)') cur_waste_form%coordinate%z
      write(id,'(a)') trim(adjustl(word))
    endif

    if (Initialized(cur_waste_form%exposure_factor)) then
      write(id,'(a29)',advance='no') 'exposure_factor: '
      write(word,'(es12.5)') cur_waste_form%exposure_factor
      write(id,'(a)') trim(adjustl(word))
    endif

    if (Initialized(cur_waste_form%volume)) then
      write(id,'(a29)',advance='no') 'volume: '
      write(word,'(es12.5)') cur_waste_form%volume
      write(id,'(a)') trim(adjustl(word)) // ' m^3'
    endif

    if (len_trim(adjustl(cur_waste_form%mech_name)) > 0) then
      write(id,'(a29)',advance='no') 'mechanism: '
      write(id,'(a)') cur_waste_form%mech_name
    endif

    if (len_trim(adjustl(cur_waste_form%spacer_mech_name)) > 0) then
      write(id,'(a29)',advance='no') 'spacer mechanism: '
      write(id,'(a)') cur_waste_form%spacer_mech_name
    endif

    if (len_trim(adjustl(cur_waste_form%criticality_mech_name)) > 0) then
      write(id,'(a29)',advance='no') 'criticality mechanism: '
      write(id,'(a)') cur_waste_form%criticality_mech_name
    endif

    if (Initialized(cur_waste_form%canister_vitality_rate)) then
      write(id,'(a29)',advance='no') 'canister vitality rate: '
      write(word,'(es12.5)') cur_waste_form%canister_vitality_rate
      write(id,'(a)') trim(adjustl(word)) // ' sec^-1'
    endif

    if (Initialized(cur_waste_form%breach_time)) then
      write(id,'(a29)',advance='no') 'canister breach time: '
      write(word,'(es12.5)') cur_waste_form%breach_time
      write(id,'(a)') trim(adjustl(word)) // ' sec'
    endif

    if (cur_waste_form%decay_start_time > 0.0d0) then
      write(id,'(a29)',advance='no') 'decay start time: '
      write(word,'(es12.5)') cur_waste_form%decay_start_time
      write(id,'(a)') trim(adjustl(word)) // ' sec'
    endif

    write(id,'(a29)') '---------------------------: '
    cur_waste_form => cur_waste_form%next
  enddo

end subroutine WasteFormInputRecord

! ************************************************************************** !

subroutine MechanismInputRecord(this)
  !
  ! Writes waste form mechanism information to the input record file.
  !
  ! Author: Alex Salazar, SNL
  ! Date: 10/07/2020
  !
  implicit none
  ! INPUT ARGUMENTS:
  ! ================
  class(wf_mechanism_base_type), pointer :: this
  ! ---------------------------------
  ! LOCAL VARIABLES:
  ! ================
  class(wf_mechanism_base_type), pointer :: cur_mech
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i
  PetscInt :: id = INPUT_RECORD_UNIT
  ! ---------------------------------

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'WASTE FORM MECHANISMS'

  cur_mech => this
  do
    if (.not. associated(cur_mech)) exit

    ! Base variables
    if (len_trim(adjustl(cur_mech%name)) > 0) then
      write(id,'(a29)',advance='no') 'mechanism name: '
      write(id,'(a)') trim(adjustl(cur_mech%name))
    endif

    if (Initialized(cur_mech%specific_surface_area)) then
      write(id,'(a29)',advance='no') 'specific surface area: '
      write(word,'(es12.5)') cur_mech%specific_surface_area
      write(id,'(a)') trim(adjustl(word)) // ' m^2/kg'
    endif

    if (Initialized(cur_mech%matrix_density)) then
      write(id,'(a29)',advance='no') 'matrix density: '
      write(word,'(es12.5)') cur_mech%matrix_density
      write(id,'(a)') trim(adjustl(word)) // ' kg/m^3'
    endif

    if (.not. cur_mech%seed == 1) then
      write(id,'(a29)',advance='no') 'random seed: '
      write(word,'(I12)') cur_mech%seed
      write(id,'(a)') trim(adjustl(word)) // ' kg/m^3'
    endif

    if (cur_mech%num_species > 0) then
      write(id,'(a29)') 'SPECIES: '
      do i = 1,cur_mech%num_species
        write(id,'(a29)',advance='no') ''
        write(id,'(a12,1X)',advance='no') &
          cur_mech%rad_species_list(i)%name
        write(id,'(es12.5,1X)',advance='no') &
          cur_mech%rad_species_list(i)%formula_weight
        write(id,'(es12.5,1X)',advance='no') &
          cur_mech%rad_species_list(i)%decay_constant
        write(id,'(es12.5,1X)',advance='no') &
          cur_mech%rad_species_list(i)%mass_fraction
        write(id,'(es12.5,1X)',advance='no') &
          cur_mech%rad_species_list(i)%inst_release_fraction
        write(id,'(a12)') cur_mech%rad_species_list(i)%daughter
      enddo
    endif

    if (cur_mech%canister_degradation_model) then

      if (Initialized(cur_mech%vitality_rate_mean)) then
        write(id,'(a29)',advance='no') 'mean degradation rate: '
        write(word,'(es12.5)') cur_mech%vitality_rate_mean
        write(id,'(a)') trim(adjustl(word)) // ' log10/yr'
      endif

      if (Initialized(cur_mech%vitality_rate_stdev)) then
        write(id,'(a29)',advance='no') 'stdev degradation rate: '
        write(word,'(es12.5)') cur_mech%vitality_rate_stdev
        write(id,'(a)') trim(adjustl(word)) // ' log10/yr'
      endif

      if (Initialized(cur_mech%vitality_rate_trunc)) then
        write(id,'(a29)',advance='no') 'degradation rate truncation: '
        write(word,'(es12.5)') cur_mech%vitality_rate_trunc
        write(id,'(a)') trim(adjustl(word)) // ' log10/yr'
      endif

      if (Initialized(cur_mech%canister_material_constant)) then
        write(id,'(a29)',advance='no') 'canister material constant: '
        write(word,'(es12.5)') cur_mech%canister_material_constant
        write(id,'(a)') trim(adjustl(word))
      endif

    endif

    ! Mechanism types
    select type(cm => cur_mech)
    class is(wf_mechanism_glass_type)

      if (cm%dissolution_rate > 0.0d0) then
        write(id,'(a29)',advance='no') 'dissolution rate: '
        write(word,'(es12.5)') cm%dissolution_rate
        write(id,'(a)') trim(adjustl(word)) // ' kg/m^2/sec'
      endif

      if (Initialized(cm%k0)) then
        write(id,'(a29)',advance='no') 'K_0 (int. dissolution rate): '
        write(word,'(es12.5)') cm%k0
        write(id,'(a)') trim(adjustl(word)) // ' kg/m^2/sec'
      endif

      if (Initialized(cm%k_long)) then
        write(id,'(a29)',advance='no') 'K_LONG (dissolution rate): '
        write(word,'(es12.5)') cm%k_long
        write(id,'(a)') trim(adjustl(word)) // ' kg/m^2/sec'
      endif

      if (Initialized(cm%nu)) then
        write(id,'(a29)',advance='no') 'nu (pH dependence): '
        write(word,'(es12.5)') cm%nu
        write(id,'(a)') trim(adjustl(word))
      endif

      if (Initialized(cm%ea)) then
        write(id,'(a29)',advance='no') 'effective activation energy: '
        write(word,'(es12.5)') cm%ea
        write(id,'(a)') trim(adjustl(word)) // ' J/mol'
      endif

      if (Initialized(cm%Q)) then
        write(id,'(a29)',advance='no') 'Q value: '
        write(word,'(es12.5)') cm%Q
        write(id,'(a)') trim(adjustl(word))
      endif

      if (Initialized(cm%K)) then
        write(id,'(a29)',advance='no') 'K (equilibrium constant): '
        write(word,'(es12.5)') cm%K
        write(id,'(a)') trim(adjustl(word))
      endif

      if (Initialized(cm%V)) then
        write(id,'(a29)',advance='no') 'V (exponent parameter): '
        write(word,'(es12.5)') cm%V
        write(id,'(a)') trim(adjustl(word))
      endif

      if (Initialized(cm%pH)) then
        write(id,'(a29)',advance='no') 'pH: '
        write(word,'(es12.5)') cm%pH
        write(id,'(a)') trim(adjustl(word))
      endif

    class is(wf_mechanism_dsnf_type)

      if (Initialized(cm%frac_dissolution_rate)) then
        write(id,'(a29)',advance='no') 'fractional dissolution rate: '
        write(word,'(es12.5)') cm%frac_dissolution_rate
        write(id,'(a)') trim(adjustl(word)) // ' sec^-1'
      endif

    class is(wf_mechanism_wipp_type)

      if (Initialized(cm%frac_dissolution_rate)) then
        write(id,'(a29)',advance='no') 'fractional dissolution rate: '
        write(word,'(es12.5)') cm%frac_dissolution_rate
        write(id,'(a)') trim(adjustl(word)) // ' sec^-1'
      endif

    class is(wf_mechanism_fmdm_type)

      if (Initialized(cm%dissolution_rate)) then
        write(id,'(a29)',advance='no') 'dissolution rate: '
        write(word,'(es12.5)') cm%dissolution_rate
        write(id,'(a)') trim(adjustl(word)) // ' kg/m^2/sec'
      endif

      if (Initialized(cm%frac_dissolution_rate)) then
        write(id,'(a29)',advance='no') 'fractional dissolution rate: '
        write(word,'(es12.5)') cm%frac_dissolution_rate
        write(id,'(a)') trim(adjustl(word)) // ' sec^-1'
      endif

      if (Initialized(cm%burnup)) then
        write(id,'(a29)',advance='no') 'burnup: '
        write(word,'(es12.5)') cm%burnup
        write(id,'(a)') trim(adjustl(word)) // ' GWd/MTHM'
      endif

    class is(wf_mechanism_fmdm_surrogate_type)

      if (Initialized(cm%dissolution_rate)) then
        write(id,'(a29)',advance='no') 'dissolution rate: '
        write(word,'(es12.5)') cm%dissolution_rate
        write(id,'(a)') trim(adjustl(word)) // ' kg/m^2/sec'
      endif

      if (Initialized(cm%frac_dissolution_rate)) then
        write(id,'(a29)',advance='no') 'fractional dissolution rate: '
        write(word,'(es12.5)') cm%frac_dissolution_rate
        write(id,'(a)') trim(adjustl(word)) // ' sec^-1'
      endif

      if (Initialized(cm%burnup)) then
        write(id,'(a29)',advance='no') 'burnup: '
        write(word,'(es12.5)') cm%burnup
        write(id,'(a)') trim(adjustl(word)) // ' GWd/MTHM'
      endif

      if (Initialized(cm%decay_time)) then
        write(id,'(a29)',advance='no') 'decay time: '
        write(word,'(es12.5)') cm%decay_time
        write(id,'(a)') trim(adjustl(word)) // ' sec'
      endif

      if (Initialized(cm%num_nearest_neighbor)) then
        write(id,'(a29)',advance='no') 'nearest neighbor: '
        write(word,'(I12)') cm%num_nearest_neighbor
        write(id,'(a)') trim(adjustl(word))
      endif

      if (Initialized(cm%knnr_eps)) then
        write(id,'(a29)',advance='no') 'KNNR EPS: '
        write(word,'(es12.5)') cm%knnr_eps
        write(id,'(a)') trim(adjustl(word))
      endif

    class is(wf_mechanism_custom_type)

      if (Initialized(cm%dissolution_rate)) then
        write(id,'(a29)',advance='no') 'dissolution rate: '
        write(word,'(es12.5)') cm%dissolution_rate
        write(id,'(a)') trim(adjustl(word)) // ' kg/m^2/sec'
      endif

      if (Initialized(cm%frac_dissolution_rate)) then
        write(id,'(a29)',advance='no') 'fractional dissolution rate: '
        write(word,'(es12.5)') cm%frac_dissolution_rate
        write(id,'(a)') trim(adjustl(word)) // ' sec^-1'
      endif

    end select

    write(id,'(a29)') '---------------------------: '
    cur_mech => cur_mech%next
  enddo

end subroutine MechanismInputRecord

! ************************************************************************** !

subroutine SpacerMechInputRecord(this)
  !
  ! Writes spacer grid degradation mechanism information to the
  !   input record file.
  !
  ! Author: Alex Salazar, SNL
  ! Date: 05/13/2021
  !
  implicit none
  ! INPUT ARGUMENTS:
  ! ================
  class(spacer_mechanism_base_type), pointer :: this
  ! ---------------------------------
  ! LOCAL VARIABLES:
  ! ================
  class(spacer_mechanism_base_type), pointer :: cur_sp_mech
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
  ! ---------------------------------

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'SPACER DEGRADATION MECHANISMS'

  cur_sp_mech => this
  do
    if (.not.associated(cur_sp_mech)) exit

    if (len_trim(adjustl(cur_sp_mech%mech_name)) > 0) then
      write(id,'(a29)',advance='no') 'spacer mechanism name: '
      write(id,'(a)') trim(adjustl(cur_sp_mech%mech_name))
    endif

    if (Initialized(cur_sp_mech%spacer_mass)) then
      write(id,'(a29)',advance='no') 'grid spc. tot. mass: '
      write(word,'(es12.5)') cur_sp_mech%spacer_mass
      write(id,'(a)') trim(adjustl(word)) // ' kg'
    endif

    if (Initialized(cur_sp_mech%spacer_surface_area)) then
      write(id,'(a29)',advance='no') 'grid spc. tot. surface area: '
      write(word,'(es12.5)') cur_sp_mech%spacer_surface_area
      write(id,'(a)') trim(adjustl(word)) // ' m^2'
    endif

    if (Initialized(cur_sp_mech%spacer_coeff)) then
      write(id,'(a29)',advance='no') 'grid spc. mech constant: '
      write(word,'(es12.5)') cur_sp_mech%spacer_coeff
      write(id,'(a)') trim(adjustl(word)) // ' kg/m^2-s'
    endif

    if (Initialized(cur_sp_mech%spacer_activation_energy)) then
      write(id,'(a29)',advance='no') 'grid spc. mech act. energy: '
      write(word,'(es12.5)') cur_sp_mech%spacer_activation_energy
      write(id,'(a)') trim(adjustl(word)) // ' J/mol'
    endif

    if (Initialized(cur_sp_mech%threshold_sat)) then
      write(id,'(a29)',advance='no') 'threshold saturation: '
      write(word,'(es12.5)') cur_sp_mech%threshold_sat
      write(id,'(a)') trim(adjustl(word))
    endif

    write(id,'(a29)') '---------------------------: '
    cur_sp_mech => cur_sp_mech%next
  enddo

end subroutine SpacerMechInputRecord

! ************************************************************************** !

subroutine CritMechInputRecord(this)
  !
  ! Writes criticality mechanism information to the input record file.
  !
  ! Author: Alex Salazar, SNL
  ! Date: 10/06/2020
  !
  implicit none
  ! INPUT ARGUMENTS:
  ! ================
  type(criticality_mediator_type), pointer :: this
  ! ---------------------------------
  ! LOCAL VARIABLES:
  ! ================
  class(crit_mechanism_base_type), pointer :: cur_crit_mech
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
  ! ---------------------------------

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'CRITICALITY MECHANISMS'

  cur_crit_mech => this%crit_mech_list
  do
    if (.not.associated(cur_crit_mech)) exit

    if (len_trim(adjustl(cur_crit_mech%mech_name)) > 0) then
      write(id,'(a29)',advance='no') 'criticality mechanism name: '
      write(id,'(a)') trim(adjustl(cur_crit_mech%mech_name))
    endif

    if (associated(cur_crit_mech%crit_event)) then
      if (Initialized(cur_crit_mech%crit_event%crit_start)) then
        write(id,'(a29)',advance='no') 'criticality start time: '
        write(word,'(es12.5)') cur_crit_mech%crit_event%crit_start
        write(id,'(a)') trim(adjustl(word)) // ' s'
      endif

      if (Initialized(cur_crit_mech%crit_event%crit_end)) then
        write(id,'(a29)',advance='no') 'criticality end time: '
        write(word,'(es12.5)') cur_crit_mech%crit_event%crit_end
        write(id,'(a)') trim(adjustl(word)) // ' s'
      endif
    endif

    if (associated(cur_crit_mech%inventory_dataset)) then
      write(id,'(a29)',advance='no') 'crit. inv. lookup table: '
      write(id,'(a)') trim(adjustl(cur_crit_mech%inventory_dataset%file_name))
    endif

    if (associated(cur_crit_mech%crit_heat_dataset)) then
      write(id,'(a29)',advance='no') 'crit. heat lookup table: '
      write(id,'(a)') trim(adjustl(cur_crit_mech%crit_heat_dataset%file_name))
    elseif (cur_crit_mech%crit_heat > 0.0d0) then
      write(id,'(a29)',advance='no') 'heat of criticality: '
      write(word,'(es12.5)') cur_crit_mech%crit_heat
      write(id,'(a)') trim(adjustl(word)) // ' MW'
    endif

    if (cur_crit_mech%sw > 0.0d0) then
      write(id,'(a29)',advance='no') 'critical water saturation: '
      write(word,'(es12.5)') cur_crit_mech%sw
      write(id,'(a)') trim(adjustl(word))
    endif

    if (cur_crit_mech%rho_w > 0.0d0) then
      write(id,'(a29)',advance='no') 'critical water density: '
      write(word,'(es12.5)') cur_crit_mech%rho_w
      write(id,'(a)') trim(adjustl(word)) // ' kg/m^3'
    endif

    if (len_trim(adjustl(cur_crit_mech%heat_dataset_name)) > 0) then
      write(id,'(a29)',advance='no') 'decay heat dataset: '
      write(id,'(a)') trim(adjustl(cur_crit_mech%heat_dataset_name))
    endif

    if (.not. cur_crit_mech%heat_source_cond == 0) then
      write(id,'(a29)',advance='no') 'decay heat source condition: '
      write(word,'(I1)') cur_crit_mech%heat_source_cond
      write(id,'(a)') trim(adjustl(word))
    endif

    if (len_trim(adjustl(cur_crit_mech%rad_dataset_name)) > 0) then
      write(id,'(a29)',advance='no') 'inventory dataset: '
      write(id,'(a)') trim(adjustl(cur_crit_mech%rad_dataset_name))
    endif

    write(id,'(a29)') '---------------------------: '
    cur_crit_mech => cur_crit_mech%next
  enddo

end subroutine CritMechInputRecord

! ************************************************************************** !

subroutine PMWFStrip(this)
  !
  ! Strips the waste form process model
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Notes: Modified by Jenn Frederick, 03/28/2016

  use Utility_module, only : DeallocateArray

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

! LOCAL VARIABLES:
! ================
! cur_waste_form: pointer to current waste form object
! prev_waste_form: pointer to previous waste form object
! -------------------------------------------------------
  class(waste_form_base_type), pointer :: cur_waste_form
  class(waste_form_base_type), pointer :: prev_waste_form
! -------------------------------------------------------

  nullify(this%realization)
  nullify(this%data_mediator)

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    prev_waste_form => cur_waste_form
    cur_waste_form => cur_waste_form%next
    call PMWFDestroyWasteForm(prev_waste_form)
  enddo
  nullify(this%waste_form_list)
  call PMWFMechanismStrip(this)
  call PMWFSpacerMechStrip(this)
  call CriticalityStrip(this%criticality_mediator)

end subroutine PMWFStrip

! ************************************************************************** !

subroutine PMWFMechanismStrip(this)
  !
  ! Strips the waste form mechanisms in the waste form process model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/28/2016
  !

  use Utility_module, only : DeallocateArray

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

! LOCAL VARIABLES:
! ================
! cur_mechanism: pointer to current mechanism object
! prev_mechanism: pointer to previous mechanism object
! --------------------------------------------------------
  class(wf_mechanism_base_type), pointer :: cur_mechanism
  class(wf_mechanism_base_type), pointer :: prev_mechanism
! --------------------------------------------------------

  cur_mechanism => this%mechanism_list
  do
    if (.not.associated(cur_mechanism)) exit
    prev_mechanism => cur_mechanism
    cur_mechanism => cur_mechanism%next
    deallocate(prev_mechanism%rad_species_list)
    nullify(prev_mechanism%rad_species_list)
    select type(prev_mechanism)
      type is(wf_mechanism_fmdm_type)
        call DeallocateArray(prev_mechanism%concentration)
        call DeallocateArray(prev_mechanism%mapping_fmdm)
        call DeallocateArray(prev_mechanism%mapping_fmdm_to_pflotran)
      type is (wf_mechanism_fmdm_surrogate_type)
        call DeallocateArray(prev_mechanism%concentration)
        call DeallocateArray(prev_mechanism%mapping_surrfmdm)
        call DeallocateArray(prev_mechanism%mapping_surrfmdm_to_pflotran)
        if (FMDM_surrogate_knnr) then
          call DeallocateArray(prev_mechanism%knnr_array)
          call DeallocateArray(prev_mechanism%table_data)
          call KdtreeDestroy(prev_mechanism%tree)
        endif
    end select
    deallocate(prev_mechanism)
    nullify(prev_mechanism)
  enddo
  nullify(this%mechanism_list)

end subroutine PMWFMechanismStrip

! ************************************************************************** !

subroutine PMWFSpacerMechStrip(this)
  !
  ! Strips the spacer grid degradation mechanisms in the waste form
  !   process model.
  !
  ! Author: Alex Salazar III
  ! Date: 05/10/2021
  !

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------

! LOCAL VARIABLES:
! ================
! cur_mechanism: pointer to current spacer grid degradation mechanism object
! prev_mechanism: pointer to previous spacer grid degradation mechanism object
! --------------------------------------------------------
  class(spacer_mechanism_base_type), pointer :: cur_mechanism
  class(spacer_mechanism_base_type), pointer :: prev_mechanism
! --------------------------------------------------------

  if (associated(this%spacer_mech_list)) then
    cur_mechanism => this%spacer_mech_list
    do
      if (.not.associated(cur_mechanism)) exit
      prev_mechanism => cur_mechanism
      cur_mechanism => cur_mechanism%next
      deallocate(prev_mechanism)
      nullify(prev_mechanism)
    enddo
    nullify(this%spacer_mech_list)
  endif

end subroutine PMWFSpacerMechStrip

! ************************************************************************** !

subroutine PMWFDestroyWasteForm(waste_form)
  !
  ! Destroys a waste form in the waste form process model
  !
  ! Author: Jenn Frederick
  ! Date: 03/28/17
  !

  use Utility_module, only : DeallocateArray

  implicit none

! INPUT ARGUMENTS:
! ================
! waste_form (input/output): waste form object to be destroyed
! --------------------------------------------------
  class(waste_form_base_type), pointer :: waste_form
! --------------------------------------------------

  call DeallocateArray(waste_form%rad_mass_fraction)
  call DeallocateArray(waste_form%rad_concentration)
  call DeallocateArray(waste_form%inst_release_amount)
  call DeallocateArray(waste_form%instantaneous_mass_rate)
  call DeallocateArray(waste_form%cumulative_mass)
  call DeallocateArray(waste_form%scaling_factor)
  call DeallocateArray(waste_form%rank_list)
  nullify(waste_form%mechanism)
  nullify(waste_form%region)
  deallocate(waste_form)
  nullify(waste_form)

end subroutine PMWFDestroyWasteForm

! ************************************************************************** !

subroutine PMWFDestroy(this)
  !
  ! Destroys the waste form process model
  !
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use String_module

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------
  character(len=MAXWORDLENGTH) :: word

  if (OptionPrintToScreen(this%option)) then
    word = StringWrite('(es12.4)',this%cumulative_time)
    write(*,'(/,a)') 'PM Waste Form time = ' // trim(adjustl(word)) // ' seconds'
  endif

  call PMBaseDestroy(this)
  call PMWFStrip(this)

end subroutine PMWFDestroy

! ************************************************************************** !

subroutine SpacerMechBaseDegradation(this,waste_form,pm,sat,temp,dt,ierr)
  !
  ! Computes spacer degradation using the base mechanism,
  ! an Arrhenius relationship.
  !

  implicit none

  class(spacer_mechanism_base_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscReal, intent(in) :: sat
  PetscReal, intent(in) :: temp ! Kelvin
  PetscReal, intent(in) :: dt
  PetscErrorCode :: ierr

  PetscReal :: dspv

  ! Spacer vitality rate - apply Arrhenius-type corrosion model [kg/m^2-s]
  waste_form%spacer_vitality_rate = this%spacer_coeff * exp(-1.0d0 * &
                                    this%spacer_activation_energy / &
                                    (IDEAL_GAS_CONSTANT * temp))

  ! Modify rate with total surface area and saturation factor [kg/s]
  waste_form%spacer_vitality_rate = waste_form%spacer_vitality_rate* &
                                    this%spacer_surface_area* &
                                    this%alteration_rate

  ! Change in spacer vitality [kg/kg]
  dspv = (1.0d0 / this%spacer_mass) * waste_form%spacer_vitality_rate * dt

  ! Spacer vitality [kg/kg]
  waste_form%spacer_vitality = waste_form%spacer_vitality - dspv

  ! Ensure value between 0 and 1
  if (waste_form%spacer_vitality > 1.0d0) then
    waste_form%spacer_vitality = 1.0d0
  elseif (waste_form%spacer_vitality < 0.0d0) then
    waste_form%spacer_vitality = 0.0d0
  endif

end subroutine SpacerMechBaseDegradation

! ************************************************************************** !

subroutine CriticalityMechInit(this)
  !
  ! Initializes the base criticality mechanism.
  !
  ! Author: Michael Nole
  ! Date: 11/01/18

  use Dataset_Ascii_class

  implicit none

  class(crit_mechanism_base_type), pointer :: this

  allocate(this)
  allocate(this%crit_event)
  nullify(this%next)

  this%mech_name = ''
  this%heat_dataset_name = ''
  this%rad_dataset_name = ''
  this%decay_heat = 0.d0
  this%crit_heat = 0.d0
  this%sw = UNINITIALIZED_DOUBLE
  this%rho_w = UNINITIALIZED_DOUBLE
  this%temperature = UNINITIALIZED_DOUBLE
  this%k_effective = 0.d0

  this%crit_event%name = ''
  this%crit_event%steady_state = PETSC_FALSE
  this%crit_event%crit_start = UNINITIALIZED_DOUBLE
  this%crit_event%crit_end = UNINITIALIZED_DOUBLE
  this%crit_event%crit_flag = PETSC_FALSE

  nullify(this%rad_dataset)
  nullify(this%heat_dataset)
  nullify(this%crit_heat_dataset)
  nullify(this%inventory_dataset)

end subroutine CriticalityMechInit

! ************************************************************************** !

subroutine CriticalityMediatorInit(this)
  !
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none

  type(criticality_mediator_type), pointer :: this

  nullify(this%data_mediator)
  nullify(this%crit_mech_list)
  this%total_num_cells = 0

end subroutine CriticalityMediatorInit

! ************************************************************************** !

function CriticalityMediatorCreate()

  !
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none

  type(criticality_mediator_type), pointer :: CriticalityMediatorCreate
  type(criticality_mediator_type), pointer :: crit

  allocate(crit)
  call CriticalityMediatorInit(crit)

  CriticalityMediatorCreate => crit

end function CriticalityMediatorCreate


! ************************************************************************** !

function CriticalityMechCreate()

  !
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none

  class(crit_mechanism_base_type), pointer :: CriticalityMechCreate
  class(crit_mechanism_base_type), pointer :: crit

  allocate(crit)
  call CriticalityMechInit(crit)

  CriticalityMechCreate => crit

end function CriticalityMechCreate

! ************************************************************************** !

subroutine ReadCriticalityMech(pmwf,input,option,keyword,error_string,found)

  !
  ! Author: Michael Nole
  ! Date: 11/01/18

  use Input_Aux_module
  use Option_module
  use String_module
  use Dataset_Ascii_class

  implicit none

  class(pm_waste_form_type) :: pmwf
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword, internal_units
  character(len=MAXSTRINGLENGTH) :: error_string,temp_string

  PetscBool :: found, added
  PetscInt :: num_errors
  PetscReal :: tempreal

  character(len=MAXWORDLENGTH) :: word
  class(crit_mechanism_base_type), pointer :: new_crit_mech, cur_crit_mech

  error_string = trim(error_string) // ',CRITICALITY_MECH'
  added = PETSC_FALSE
  found = PETSC_TRUE
  num_errors = 0
  tempreal = UNINITIALIZED_DOUBLE
  select case(trim(keyword))
  !-------------------------------------
    case('CRITICALITY_MECH')
      allocate(new_crit_mech)
      new_crit_mech => CriticalityMechCreate()
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input, option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case (trim(word))
        !-------------------------------------
          case('NAME')
            call InputReadCard(input,option,word)
            call InputErrorMsg(input,option, &
                  'criticality mechanism assignment',error_string)
            call StringToUpper(word)
            new_crit_mech%mech_name = trim(word)
        !-------------------------------------
          case('CRIT_START')
            call InputReadDouble(input,option,new_crit_mech% &
                                 crit_event%crit_start)
            call InputErrorMsg(input,option,'CRIT_START',error_string)
            call InputReadAndConvertUnits(input,new_crit_mech% &
                                 crit_event%crit_start,'sec', &
                                 trim(error_string)//',CRIT_START', &
                                 option)
        !-------------------------------------
          case('CRIT_END')
            call InputReadDouble(input,option,new_crit_mech% &
                                 crit_event%crit_end)
            call InputErrorMsg(input,option,'CRIT_END',error_string)
            call InputReadAndConvertUnits(input,new_crit_mech% &
                                          crit_event%crit_end,'sec', &
                                          trim(error_string)//',CRIT_END', &
                                          option)
        !-------------------------------------
          case('CRITICAL_WATER_SATURATION')
            call InputReadDouble(input,option,new_crit_mech%sw)
            call InputErrorMsg(input,option,'CRITICAL WATER SATURATION',&
              error_string)
        !-------------------------------------
          case('CRITICAL_WATER_DENSITY')
            call InputReadDouble(input,option,new_crit_mech%rho_w)
            call InputErrorMsg(input,option,'CRITICAL WATER DENSITY',&
              error_string)
            call InputReadAndConvertUnits(input,new_crit_mech%rho_w, &
              'kg/m^3',trim(error_string)//',CRITICAL WATER DENSITY',option)
        !-------------------------------------
          case('HEAT_OF_CRITICALITY')
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit
              call InputReadCard(input,option,word,PETSC_FALSE)
              select case(trim(word))
              !-------------------------------------
                case('CONSTANT_POWER')
                  internal_units = 'MW'
                  call InputReadDouble(input,option,new_crit_mech%crit_heat)
                  call InputErrorMsg(input,option,'HEAT_OF_CRITICALITY, ' &
                    //'CONSTANT_POWER',error_string)
                  call InputReadAndConvertUnits(input,new_crit_mech%crit_heat, &
                    internal_units,trim(error_string)//',HEAT_OF_CRITICALITY,' &
                    //' CONSTANT_POWER',option)
              !-------------------------------------
                case('DATASET')
                  allocate(new_crit_mech%crit_heat_dataset)
                  new_crit_mech%crit_heat_dataset => CritHeatCreate()
                  call InputReadFilename(input,option,new_crit_mech% &
                                         crit_heat_dataset%file_name)
                  call new_crit_mech%crit_heat_dataset%Read(new_crit_mech% &
                                                            crit_heat_dataset% &
                                                            file_name,option)
              !-----------------------------
                case default
                  call InputKeywordUnrecognized(input,word,error_string,option)
              !-----------------------------
              end select
            enddo
            call InputPopBlock(input,option)
        !-------------------------------------
          case('DECAY_HEAT')
            call InputReadCard(input,option,word)
            select case (trim(word))
            !-----------------------------
              case('TOTAL')
                new_crit_mech%heat_source_cond = 1
            !-----------------------------
              case('ADDITIONAL')
                new_crit_mech%heat_source_cond = 2
            !-----------------------------
              case('CYCLIC')
                new_crit_mech%heat_source_cond = 3
            !-----------------------------
            end select
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit
              call InputReadCard(input,option,word,PETSC_FALSE)
              select case(trim(word))
              !-----------------------------
                case('DATASET')
                  internal_units = 'MW'
                  new_crit_mech%heat_dataset => DatasetAsciiCreate()
                  call InputReadFilename(input,option,new_crit_mech% &
                          heat_dataset_name)
                  call DatasetAsciiReadFile(new_crit_mech%heat_dataset, &
                          new_crit_mech%heat_dataset_name,temp_string, &
                          internal_units,error_string,option)
                  new_crit_mech%heat_dataset%time_storage% &
                          time_interpolation_method = 2
              !-----------------------------
                case default
                  call InputKeywordUnrecognized(input,word,error_string,option)
              !-----------------------------
              end select
            enddo
            call InputPopBlock(input,option)
        !-------------------------------------
          case('INVENTORY')
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit
              call InputReadCard(input,option,word,PETSC_FALSE)
              select case(trim(word))
              !-----------------------------
                case('DATASET')
                  internal_units = 'g/g'
                  new_crit_mech%rad_dataset => DatasetAsciiCreate()
                  call InputReadFilename(input,option,new_crit_mech% &
                          rad_dataset_name)
                  call DatasetAsciiReadFile(new_crit_mech%rad_dataset, &
                          new_crit_mech%rad_dataset_name,temp_string, &
                          internal_units,error_string,option)
                  new_crit_mech%rad_dataset%time_storage% &
                          time_interpolation_method = 2
              !-----------------------------
                case('EXPANDED_DATASET')
                  allocate(new_crit_mech%inventory_dataset)
                  new_crit_mech%inventory_dataset => CritInventoryCreate()
                  call InputReadFilename(input,option,new_crit_mech% &
                         inventory_dataset%file_name)
                  call new_crit_mech%inventory_dataset%Read(new_crit_mech% &
                                                            inventory_dataset% &
                                                            file_name,option)
              !-----------------------------
                case('OPTION')
                  call InputPushBlock(input,option)
                  do
                    call InputReadPflotranString(input,option)
                    if (InputError(input)) exit
                    if (InputCheckExit(input,option)) exit
                    call InputReadCard(input,option,word,PETSC_FALSE)
                    select case(trim(word))
                    !-----------------------------
                      case('USE_LOOKUP_AND_IMPLICIT')
                        ! Implicit solution becomes fallback if extrapolation
                        !   is needed from real time array
                        if (associated(new_crit_mech%inventory_dataset)) then
                          new_crit_mech%inventory_dataset%allow_implicit = &
                            PETSC_TRUE
                          ! Detect conflicting options
                          if (new_crit_mech%inventory_dataset%allow_extrap) then
                            option%io_buffer = 'Option "' // trim(word) // &
                                               '" conflicts with others listed.'
                            call PrintErrMsg(option)
                          endif
                        else
                          option%io_buffer = 'Option "' // trim(word) // &
                                             '" requires specification of ' //&
                                             'EXPANDED_DATASET.'
                          call PrintErrMsg(option)
                        endif
                    !-----------------------------
                      case('USE_LOOKUP_AND_EXTRAPOLATION')
                        ! Allow extrapolation from lookup table
                        if (associated(new_crit_mech%inventory_dataset)) then
                          new_crit_mech%inventory_dataset%allow_extrap = &
                            PETSC_TRUE
                          ! Detect conflicting options
                          if (new_crit_mech%inventory_dataset% &
                              allow_implicit) then
                            option%io_buffer = 'Option "' // trim(word) // &
                                               '" conflicts with others listed.'
                            call PrintErrMsg(option)
                          endif
                        else
                          option%io_buffer = 'Option "' // trim(word) // &
                                             '" requires specification of ' //&
                                             'EXPANDED_DATASET.'
                          call PrintErrMsg(option)
                        endif
                    !-----------------------------
                      case('USE_LOOKUP_AFTER_CRITICALITY')
                        ! This allows the lookup table to be used after the
                        !   criticality end time. Otherwise, the implicit
                        !   solution is used after the criticality event.
                        if (associated(new_crit_mech%inventory_dataset)) then
                          new_crit_mech%inventory_dataset%continue_lookup = &
                            PETSC_TRUE
                        else
                          option%io_buffer = 'Option "' // trim(word) // &
                                             '" requires specification of ' //&
                                             'EXPANDED_DATASET.'
                          call PrintErrMsg(option)
                        endif
                    !-----------------------------
                      case('LOG10_TIME_INTERPOLATION')
                        ! Use log10 basis for time interpolation
                        if (.not. associated(new_crit_mech% &
                                             inventory_dataset)) then
                          option%io_buffer = 'Option "' // trim(word) // &
                                             '" requires specification of ' //&
                                             'EXPANDED_DATASET.'
                          call PrintErrMsg(option)
                        endif

                        ! User can specify substitute value for zero (optional)
                        internal_units = 's'
                        call InputReadDouble(input,option,tempreal)
                        if (Initialized(tempreal)) then
                          call InputReadAndConvertUnits(input,tempreal, &
                                 internal_units,trim(error_string) &
                                 //',EXPANDED_DATASET,OPTION' &
                                 //',LOG10_TIME_INTERPOLATION,' &
                                 //' zero substitute (' &
                                 //trim(new_crit_mech%inventory_dataset% &
                                 file_name)//')',option)
                          call CritInventoryCheckZeroSub(tempreal, &
                                                         new_crit_mech% &
                                                         inventory_dataset, &
                                                         option)
                        endif

                        ! Modify lookup tables
                        call CritInventoryUseLog10(new_crit_mech% &
                                                   inventory_dataset,option)
                    !-----------------------------
                      case default
                        call InputKeywordUnrecognized(input,word,error_string, &
                                                      option)
                    !-----------------------------
                    end select
                  enddo
                  call InputPopBlock(input,option)
              !-----------------------------
                case default
                  call InputKeywordUnrecognized(input,word,error_string,option)
              !-----------------------------
              end select
            enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-------------------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      ! --------------------------- error messaging ---------------------------
      if (len_trim(new_crit_mech%mech_name) < 1) then
        option%io_buffer = 'Name must be specified for criticality mechanism ' &
                         //'in order to be associated with a waste form.'
        call PrintWrnMsg(option)
      endif

      if (Uninitialized(new_crit_mech%crit_event%crit_start)) then
        option%io_buffer = 'ERROR: Criticality start time must be specified ' &
                         //'for criticality mechanism.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif

      if (Uninitialized(new_crit_mech%crit_event%crit_end)) then
        option%io_buffer = 'ERROR: Criticality end time must be specified ' &
                         //'for criticality mechanism.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif

      if (Initialized(new_crit_mech%sw)) then
        if (new_crit_mech%sw > 1.d0 .or. new_crit_mech%sw < 0.d0) then
          option%io_buffer = 'ERROR: Critical water saturation must be ' &
                           //'between 0 and 1 in the criticality mechanism.'
          call PrintMsg(option)
          num_errors = num_errors + 1
        endif
      endif

      if (associated(new_crit_mech%crit_heat_dataset)) then
        if (Initialized(new_crit_mech%crit_event%crit_start)) then
          if (new_crit_mech%crit_event%crit_start > &
              new_crit_mech%crit_heat_dataset%start_time_datamax) then
            option%io_buffer = 'ERROR: Criticality start time exceeds ' &
                             //'maximum value in heat of criticality lookup ' &
                             //'table.'
            call PrintMsg(option)
            num_errors = num_errors + 1
          endif
        endif
      endif

      if (.not.associated(pmwf%criticality_mediator)) then
        pmwf%criticality_mediator => CriticalityMediatorCreate()
      endif

      if (.not. associated(pmwf%criticality_mediator%crit_mech_list)) then
        pmwf%criticality_mediator%crit_mech_list => new_crit_mech
      else
        cur_crit_mech => pmwf%criticality_mediator%crit_mech_list
        do
          if (.not. associated(cur_crit_mech)) exit
          if (.not. associated(cur_crit_mech%next)) then
            cur_crit_mech%next => new_crit_mech
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_crit_mech => cur_crit_mech%next
        enddo
      endif
      nullify(new_crit_mech)
  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in ' &
                    //'the WASTE_FORM_GENERAL,CRITICALITY_MECH ' &
                    //'block(s). See above.'
    call PrintErrMsg(option)
  endif

end subroutine ReadCriticalityMech

! ************************************************************************** !

subroutine CriticalityCalc(this,time,ierr)

  ! Calculate mass and heat source terms as a function of time.
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none

  class(crit_mechanism_base_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  PetscReal :: t_low, t_high
  PetscReal, pointer :: times(:)
  class(dataset_ascii_type), pointer :: dataset
  PetscInt :: j

  dataset => this%heat_dataset
  if (associated(dataset%time_storage)) then
    times => dataset%time_storage%times
    j=1
    t_low = times(j)
    t_high = t_low
    do
      if(time < times(j)) exit
      if(j == size(times)) exit
      t_low = times(j)
      j = j+1
      t_high = times(j)
    enddo

    if (j == size(times) .and. time >= times(j)) then
      this%decay_heat = dataset%rbuffer(j)
    elseif (j==1) then
      this%decay_heat = 0.d0
    else
      this%decay_heat = dataset%rbuffer(j-1) + (time-t_low)/(t_high-t_low)* &
               (dataset%rbuffer(j)-dataset%rbuffer(j-1))
    endif
  else
    this%decay_heat = 0.d0
  endif
end subroutine CriticalityCalc

! ************************************************************************** !

subroutine CritReadValues(input, option, keyword, dataset_base, &
                          data_external_units, data_internal_units)
  use Input_Aux_module
  use String_module
  use Option_module
  use Logging_module
  use HDF5_Aux_module
  use Units_module
  use Dataset_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use hdf5

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  class(dataset_base_type), pointer :: dataset_base
  character(len=*) :: data_external_units
  character(len=*) :: data_internal_units

  character(len=MAXSTRINGLENGTH), pointer :: internal_unit_strings(:)
  class(dataset_ascii_type), pointer :: dataset_ascii
  character(len=MAXSTRINGLENGTH) :: string2, filename, hdf5_path
  character(len=MAXWORDLENGTH) :: word, realization_word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: length, i
  PetscErrorCode :: ierr

  ! dataset_base, though of type dataset_base_type, should always be created
  ! as dataset_ascii_type.
  dataset_ascii => DatasetAsciiCast(dataset_base)
  if (.not.associated(dataset_ascii)) then
    ! The dataset was not of type dataset_ascii and was likely set to a different
    ! type.  There is a bug in the input file.
    option%io_buffer = 'Dataset associated with ' // trim(keyword) // &
      ' in the input file is already associated with a different dataset &
      &type.  Check for duplicate definitions of ' // trim(keyword) // '.'
    call PrintErrMsg(option)
  endif

  filename = ''
  realization_word = ''
  hdf5_path = ''

  internal_unit_strings => StringSplit(data_internal_units,',')

  input%ierr = 0
  string2 = trim(input%buf)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'file or value','CONDITION')
  call StringToUpper(word)
  length = len_trim(word)
  if (StringStartsWithAlpha(word)) then
    call InputPushCard(input,word,option)
    if (length == FOUR_INTEGER .and. &
        StringCompare(word,'file',FOUR_INTEGER)) then
      input%err_buf2 = trim(keyword) // ', FILE'
      input%err_buf = 'keyword'
      call InputReadFilename(input,option,string2)
      if (input%ierr == 0) then
        filename = string2
      else
        option%io_buffer = 'The ability to read realization dependent &
          &datasets outside the DATASET block is no longer supported'
        call PrintErrMsg(option)
      endif

      if (len_trim(filename) < 2) then
        option%io_buffer = 'No filename listed under Flow_Condition: ' // &
                           trim(keyword)
        call PrintErrMsg(option)
      endif
      if (index(filename,'.h5') > 0) then
        write(option%io_buffer,'("Reading of HDF5 datasets for flow ", &
                                 &"conditions not currently supported.")')
        call PrintErrMsg(option)
      else
        i = index(filename,'.',PETSC_TRUE)
        if (i > 2) then
          filename = filename(1:i-1) // trim(realization_word) // filename(i:)
        else
          filename = trim(filename) // trim(realization_word)
        endif
        error_string = 'CONDITION,' // trim(keyword) // ',FILE'
        call DatasetAsciiReadFile(dataset_ascii,filename,data_external_units, &
                                  data_internal_units,error_string,option)
        dataset_ascii%filename = filename
      endif

    else if (StringCompare(word,'dataset')) then
      call InputReadWord(input,option,word,PETSC_TRUE)
      input%err_buf2 = trim(keyword) // ', DATASET'
      input%err_buf = 'dataset name'
      call InputErrorMsg(input,option)
      call DatasetDestroy(dataset_base)
      dataset_base => DatasetBaseCreate()
      dataset_base%name = word
    else if (length==FOUR_INTEGER .and. StringCompare(word,'list',length)) then
      error_string = 'CONDITION,' // trim(keyword) // ',LIST'
      call DatasetAsciiReadList(dataset_ascii,input,data_external_units, &
                                data_internal_units,error_string,option)
    else
      option%io_buffer = 'Keyword "' // trim(word) // &
        '" not recognized in when reading condition values for "' // &
        trim(keyword) // '".'
      call PrintErrMsg(option)
    endif
  else
    input%buf = trim(string2)
    error_string = 'CONDITION,' // trim(keyword) // ',SINGLE'
    call DatasetAsciiReadSingle(dataset_ascii,input,data_external_units, &
                                data_internal_units,error_string,option)
#if 0
    allocate(dataset_ascii%rarray(dataset_ascii%array_width))
    do icol=1,dataset_ascii%array_width
      call InputReadDouble(input,option,dataset_ascii%rarray(icol))
      write(input%err_buf,'(a,i2)') trim(keyword) // &
                                    ' dataset_values, icol = ', icol
      input%err_buf2 = 'CONDITION'
      call InputErrorMsg(input,option)
    enddo
    string2 = input%buf
    call InputReadWord(input,option,word,PETSC_TRUE)
    if (InputError(input)) then
      call InputCheckMandatoryUnits(input,option)
      word = trim(keyword) // ' UNITS'
      call InputDefaultMsg(input,option,word)
    else
      input%buf = string2
      units = ''
      do icol=1,dataset_ascii%array_width
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,'CONDITION')
        dataset_ascii%rarray(icol) = UnitsConvertToInternal(word, &
                                     internal_unit_strings(icol),option) * &
                                     dataset_ascii%rarray(icol)
        units = trim(units) // ' ' // trim(word)
      enddo
    endif
#endif
  endif

  deallocate(internal_unit_strings)
  nullify(internal_unit_strings)

  call PetscLogEventEnd(logging%event_flow_condition_read_values, &
                        ierr);CHKERRQ(ierr)

end subroutine CritReadValues

! ************************************************************************** !

subroutine CriticalityStrip(this)

  !
  ! Author: Michael Nole
  ! Date: 11/01/18
  implicit none

  type(criticality_mediator_type), pointer :: this
  nullify(this)


end subroutine CriticalityStrip

! ************************************************************************** !

subroutine CritHeatRead(this,filename,option)
  !
  ! Author: Alex Salazar III
  ! Date: 05/12/2021
  !
  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module
  use Units_module

  implicit none

  class(crit_heat_type) :: this
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: keyword, word, internal_units
  character(len=MAXSTRINGLENGTH) :: error_string
  type(input_type), pointer :: input2
  PetscInt :: temp_int
  PetscReal :: time_units_conversion
  PetscReal :: temp_units_conversion
  PetscReal :: power_units_conversion

  time_units_conversion = 1.d0
  temp_units_conversion = 1.d0
  power_units_conversion = 1.d0

  if (len_trim(filename) < 1) then
    option%io_buffer = 'Filename must be specified for heat of criticality ' &
                     //'lookup table.'
    call PrintErrMsg(option)
  endif

  this%lookup_table => LookupTableCreateGeneral(TWO_INTEGER)
  error_string = 'heat of criticality lookup table'
  input2 => InputCreate(IUNIT_TEMP,filename,option)
  input2%ierr = 0
  do
    call InputReadPflotranString(input2,option)
    if (InputError(input2)) exit

    call InputReadCard(input2,option,keyword)
    call InputErrorMsg(input2,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      case('NUM_START_TIMES')
        call InputReadInt(input2,option,this%num_start_times)
        call InputErrorMsg(input2,option,'number of start times',error_string)
      case('NUM_VALUES_PER_START_TIME')
        call InputReadInt(input2,option,this%num_values_per_start_time)
        call InputErrorMsg(input2,option,'number of values per start time', &
                           error_string)
      case('TIME_UNITS')
        internal_units = 'sec'
        call InputReadWord(input2,option,word,PETSC_TRUE)
        call InputErrorMsg(input2,option,'UNITS','CONDITION')
        time_units_conversion = UnitsConvertToInternal(word, &
                                internal_units,option)
      case('TEMPERATURE_UNITS')
        internal_units = 'C'
        call InputReadWord(input2,option,word,PETSC_TRUE)
        call InputErrorMsg(input2,option,'UNITS','CONDITION')
        call StringToUpper(word)
        temp_units_conversion = UnitsConvertToInternal(word, &
                                internal_units,option)
      case('POWER_UNITS')
        internal_units = 'MW'
        call InputReadWord(input2,option,word,PETSC_TRUE)
        call InputErrorMsg(input2,option,'UNITS','CONDITION')
        power_units_conversion = UnitsConvertToInternal(word, &
                                 internal_units,option)
      case('START_TIME')
        if (Uninitialized(this%num_start_times) .or. &
            Uninitialized(this%num_values_per_start_time)) then
          option%io_buffer = 'NUM_START_TIMES and NUM_VALUES_PER_START_TIME ' &
                           //'must be specified prior to reading the ' &
                           //'corresponding arrays.'
          call PrintErrMsg(option)
        endif
        this%lookup_table%dims(1) = this%num_start_times
        this%lookup_table%dims(2) = this%num_values_per_start_time
        temp_int = this%num_start_times*this%num_values_per_start_time
        allocate(this%lookup_table%axis1%values(this%num_start_times))
        allocate(this%lookup_table%axis2%values(temp_int))
        allocate(this%lookup_table%data(temp_int))
        string = 'START_TIME in heat of criticality lookup table'
        call UtilityReadArray(this%lookup_table%axis1%values, &
                              NEG_ONE_INTEGER,string, &
                              input2,option)
        this%lookup_table%axis1%values = this%lookup_table%axis1%values * &
          time_units_conversion
      case('TEMPERATURE')
        string = 'TEMPERATURE in heat of criticality lookup table'
        call UtilityReadArray(this%lookup_table%axis2%values, &
                              NEG_ONE_INTEGER, &
                              string,input2,option)
        this%lookup_table%axis2%values = this%lookup_table%axis2%values * &
          temp_units_conversion
      case('POWER')
        string = 'POWER in heat of criticality lookup table'
        call UtilityReadArray(this%lookup_table%data, &
                              NEG_ONE_INTEGER, &
                              string,input2,option)
        this%lookup_table%data = this%lookup_table%data * power_units_conversion
     case default
        error_string = trim(error_string) // ': ' // filename
        call InputKeywordUnrecognized(input2,keyword,error_string,option)
    end select
  enddo
  call InputDestroy(input2)

  if (size(this%lookup_table%axis1%values) /= this%num_start_times) then
    option%io_buffer = 'Number of start times does not match NUM_START_TIMES.'
    call PrintErrMsg(option)
  endif
  if (size(this%lookup_table%axis2%values) /= &
      this%num_start_times*this%num_values_per_start_time) then
    option%io_buffer = 'Number of temperatures does not match ' &
                     //'NUM_START_TIMES * NUM_VALUES_PER_START_TIME.'
    call PrintErrMsg(option)
  endif
  if (size(this%lookup_table%data) /= &
      this%num_start_times*this%num_values_per_start_time) then
    option%io_buffer = 'Number of criticality power outputs does not match ' &
                     //'NUM_START_TIMES * NUM_VALUES_PER_START_TIME.'
    call PrintErrMsg(option)
  endif

  ! set limits
  this%start_time_datamax = maxval(this%lookup_table%axis1%values)
  this%temp_datamax = maxval(this%lookup_table%axis2%values)
  this%power_datamax = maxval(this%lookup_table%data)

end subroutine CritHeatRead

! ************************************************************************** !

subroutine CritInventoryRead(this,filename,option)
  !
  ! Reads in tables of time-dependent radionuclide inventories parametrized by
  !   the criticality start time and power output.
  !
  ! Author: Alex Salazar III
  ! Date: 02/16/2022
  !
  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module
  use Units_module
  !
  implicit none
  ! ----------------------------------
  class(crit_inventory_type) :: this
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  ! ----------------------------------
  character(len=MAXSTRINGLENGTH) :: string, error_string, str1, str2, str3
  character(len=MAXWORDLENGTH) :: keyword, word, internal_units
  type(input_type), pointer :: input
  PetscInt :: i
  PetscInt :: ict
  PetscInt :: num_partitions
  PetscInt :: arr_size
  PetscReal :: time_units_conversion
  PetscReal :: power_units_conversion
  PetscReal :: data_units_conversion
  PetscReal, pointer :: tmpaxis1(:)
  PetscReal, pointer :: tmpaxis2(:)
  PetscReal, pointer :: tmpaxis3(:)
  PetscReal, pointer :: tmpdata(:) ! temp inventory data pointer from input read
  type(lookup_irregular_array_type), allocatable :: tmpdata1(:) ! saves data - old
  type(lookup_irregular_array_type), allocatable :: tmpdata2(:) ! saves data - new
  character(len=MAXWORDLENGTH), allocatable :: tmpname1(:)
  character(len=MAXWORDLENGTH), allocatable :: tmpname2(:)
  PetscInt :: tmpdatarank ! rank of array before new inventory is detected
  PetscInt :: mode
  PetscInt, parameter :: mode_lp = 1
  PetscInt, parameter :: mode_tl = 2
  type(crit_inventory_lookup_type), pointer :: cur_inventory, new_inventory
  ! ----------------------------------

  ict = 0
  num_partitions = 0
  arr_size = 0
  time_units_conversion  = 1.d0
  power_units_conversion = 1.d0
  data_units_conversion  = 1.d0
  mode = UNINITIALIZED_INTEGER

  tmpdatarank = 1
  allocate(tmpdata1(tmpdatarank))
  allocate(tmpdata2(tmpdatarank))
  allocate(tmpname1(tmpdatarank))
  allocate(tmpname2(tmpdatarank))

  if (len_trim(filename) < 1) then
    option%io_buffer = 'Filename must be specified for criticality inventory ' &
                     //'lookup table.'
    call PrintErrMsg(option)
  endif

  ! this%lookup_table => LookupTableCreateGeneral(THREE_INTEGER) !3D interpolation
  error_string = 'criticality inventory lookup table "' // trim(filename) // '"'
  input => InputCreate(IUNIT_TEMP,filename,option)
  input%ierr = 0
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
    !-------------------------------------
      case('MODE')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'interpolation mode',error_string)
        call StringToUpper(word)
        select case(word)
        !-------------------------------------
          case('POLYNOMIAL','LAGRANGE_POLYNOMIAL')
            mode = mode_lp
        !-------------------------------------
          case('LINEAR','TRILINEAR')
            mode = mode_tl
        !-------------------------------------
          case default
            call InputKeywordUnrecognized(input,word, &
                   'interpolation mode',option)
        end select
    !-------------------------------------
      case('NUM_START_TIMES')
        call InputReadInt(input,option,this%num_start_times)
        call InputErrorMsg(input,option,'number of start times',error_string)
    !-------------------------------------
      case('NUM_POWERS')
        call InputReadInt(input,option,this%num_powers)
        call InputErrorMsg(input,option,'number of power outputs',error_string)
    !-------------------------------------
      case('NUM_REAL_TIMES')
        call InputReadInt(input,option,this%num_real_times)
        call InputErrorMsg(input,option,'maximum length of inventory ' &
                                      //'evaluation times',error_string)
    !-------------------------------------
      case('TOTAL_POINTS')
        call InputReadInt(input,option,this%total_points)
        call InputErrorMsg(input,option,'total inventory evaluation points', &
                           error_string)
    !-------------------------------------
      case('NUM_SPECIES')
        call InputReadInt(input,option,this%num_species)
        call InputErrorMsg(input,option,'number of species in inventory', &
                           error_string)
    !-------------------------------------
      case('TIME_UNITS')
        internal_units = 'sec'
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'time units',error_string)
        time_units_conversion = UnitsConvertToInternal(word, &
                                internal_units,option)
    !-------------------------------------
      case('POWER_UNITS')
        internal_units = 'MW'
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'power units',error_string)
        power_units_conversion = UnitsConvertToInternal(word, &
                                 internal_units,option)
    !-------------------------------------
      case('DATA_UNITS')
        internal_units = 'g/g'
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'data units',error_string)
        data_units_conversion = UnitsConvertToInternal(word, &
                                 internal_units,option)
    !-------------------------------------
      case('START_TIME')
        string = 'START_TIME in criticality inventory lookup table "' &
                 // trim(filename) // '"'

        nullify(tmpaxis1)

        call UtilityReadArray(tmpaxis1, &
                              NEG_ONE_INTEGER,string, &
                              input,option)
    !-------------------------------------
      case('POWER')
        string = 'POWER in criticality inventory lookup table "' &
                 // trim(filename) // '"'

        nullify(tmpaxis2)

        call UtilityReadArray(tmpaxis2, &
                              NEG_ONE_INTEGER, &
                              string,input,option)
    !-------------------------------------
      case('REAL_TIME')
        string = 'REAL_TIME in criticality inventory lookup table "' &
                 // trim(filename) // '"'

        nullify(tmpaxis3)

        call UtilityReadArray(tmpaxis3, &
                              NEG_ONE_INTEGER, &
                              string,input,option)
    !-------------------------------------
      case('INVENTORY','INVENTORIES')
        ! NEST ORDER
        !
        ! I  start time i
        ! I  start time i+1
        ! I  ...
        ! I  start time imax
        !
        ! J    power output j
        ! J    power output j+1
        ! J    ...
        ! J    power output jmax
        !
        ! K      start i    power j   : times(i   ,j   )
        ! K      start i    power j+1 : times(i   ,j+1 )
        ! K      ...
        ! K      start i    power jmax: times(i  ,jmax )
        ! K      start i+1  power j   : times(i+1,j    )
        ! K      ...
        ! K      start imax power jmax: times(imax,jmax)
        !
        !  There is a lookup table for each nuclide z
        !      lookup_table(1)
        ! I,J,K    start i    power j    nuclide z   : inventory(t) for times(i   ,j   )
        ! I,J,K    start i    power j+1  nuclide z   : inventory(t) for times(i   ,j+1 )
        ! I,J,K    start i    power j+2  nuclide z   : inventory(t) for times(i   ,j+2 )
        ! I,J,K    ...
        ! I,J,K    start i    power jmax nuclide z   : inventory(t) for times(i   ,jmax)
        ! I,J,K    start i+1  power j    nuclide z   : inventory(t) for times(i+1 ,jmax)
        ! I,J,K    start i+1  power j+1  nuclide z   : inventory(t) for times(i+1 ,j+1 )
        ! I,J,K    ...
        ! I,J,K    start imax power jmax nuclide z   : inventory(t) for times(imax,jmax)
        !      lookup_table(2)
        ! I,J,K    start i    power j    nuclide z+1 : inventory(t) for times(i   ,j   )
        ! I,J,K    start i    power j+1  nuclide z+1 : inventory(t) for times(i   ,j+1 )
        ! I,J,K    ...
        ! I,J,K
        ! I,J,K    start imax power jmax nuclide zmax: inventory(t) for times(imax,jmax)

        ict = ict + 1

        write(str1,*) ict

        string = 'INVENTORY (#'// trim(adjustl(str1)) //')'

        ! Get species name
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,string //': species name',error_string)

        ! Re-allocate temporary arrays of data
        if (ict > tmpdatarank) then
          tmpdata1(1:tmpdatarank) = tmpdata2(1:tmpdatarank)
          if (allocated(tmpdata2)) deallocate(tmpdata2)
          allocate(tmpdata2(ict))
          tmpdata2(1:tmpdatarank) = tmpdata1(1:tmpdatarank)
          if (allocated(tmpdata1)) deallocate(tmpdata1)
          allocate(tmpdata1(ict))
          tmpdata1(1:tmpdatarank) = tmpdata2(1:tmpdatarank)

          tmpname1(1:tmpdatarank) = tmpname2(1:tmpdatarank)
          if (allocated(tmpname2)) deallocate(tmpname2)
          allocate(tmpname2(ict))
          tmpname2(1:tmpdatarank) = tmpname1(1:tmpdatarank)
          if (allocated(tmpname1)) deallocate(tmpname1)
          allocate(tmpname1(ict))
          tmpname1(1:tmpdatarank) = tmpname2(1:tmpdatarank)

          tmpdatarank = ict
        endif

        nullify(tmpdata)

        call UtilityReadArray(tmpdata, &
                              NEG_ONE_INTEGER, &
                              string // " : " // error_string,input,option)
        allocate(tmpdata2(ict)%data(1:size(tmpdata)))

        tmpname2(ict) = word ! save species name
        tmpdata2(ict)%data = tmpdata ! save data from inventory block
    !-------------------------------------
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputDestroy(input)

  ! Check for errors after input read
  if (Uninitialized(this%num_start_times)) then
    this%num_start_times = size(tmpaxis1)
  elseif (this%num_start_times /= size(tmpaxis1)) then
    write(str1,*) size(tmpaxis1)
    option%io_buffer = 'Number of start times listed in START_TIME ('&
                     // trim(adjustl(str1)) &
                     //') does not match NUM_START_TIMES in "' &
                     // trim(filename) // '".'
    call PrintErrMsg(option)
  endif

  if (Uninitialized(this%num_powers)) then
    this%num_powers = size(tmpaxis2)
  elseif (this%num_powers /= size(tmpaxis2)) then
    write(str1,*) size(tmpaxis2)
    option%io_buffer = 'Number of powers listed in POWER ('&
                     // trim(adjustl(str1)) &
                     //') does not match NUM_POWERS in "' &
                     // trim(filename) // '".'
    call PrintErrMsg(option)
  endif

  if (Uninitialized(this%total_points) .and. &
      Uninitialized(this%num_real_times)) then
    this%total_points = size(tmpaxis3)
  endif

  if (Initialized(this%total_points)) then
    arr_size = this%total_points
  else
    if (Uninitialized(this%num_real_times)) then
      option%io_buffer = 'NUM_REAL_TIMES must be specified in "' &
                       // trim(filename) // '" to provide axis3 dimension. '&
                       //'TOTAL_POINTS may also be used if the axis3 lengths '&
                       //'are non-uniform.'
      call PrintErrMsg(option)
    endif
    arr_size = this%num_start_times*this%num_powers*this%num_real_times
  endif

  if (size(tmpaxis3) /= arr_size) then
    write(str1,*) size(tmpaxis3)
    if (Initialized(this%total_points)) then
      option%io_buffer = 'Number of points listed for REAL_TIME ('&
                       // trim(adjustl(str1)) &
                       //') does not match TOTAL_POINTS in "' &
                       // trim(filename) // '".'
                       call PrintErrMsg(option)
    else
      option%io_buffer = 'Number of points listed for REAL_TIME ('&
                       // trim(adjustl(str1)) &
                       //') does not match ' &
                       //'NUM_START_TIMES * NUM_POWERS * NUM_REAL_TIMES in "' &
                       // trim(filename) // '".'
                       call PrintWrnMsg(option)
    endif
  endif

  if (Uninitialized(this%num_species)) then
    if (ict > 0) then
      this%num_species = ict
    else
      option%io_buffer = 'No INVENTORY blocks were detected in "' &
                       // trim(filename) // '".'
      call PrintErrMsg(option)
    endif
  endif

  if (ict /= this%num_species) then
    write(str1,*)ict
    write(str2,*)this%num_species
    option%io_buffer = 'Total number of inventories in data table (' &
                     // trim(adjustl(str1)) &
                     //') does not match NUM_SPECIES specified (' &
                     // trim(adjustl(str2)) //') in "'// trim(filename) // '".'
    call PrintErrMsg(option)
  endif

  ! Adjust units of axes if neccessary and record maxima
  tmpaxis1 = tmpaxis1 * time_units_conversion
  tmpaxis2 = tmpaxis2 * power_units_conversion
  tmpaxis3 = tmpaxis3 * time_units_conversion
  this%start_time_datamax = maxval(tmpaxis1)
  this%power_datamax = maxval(tmpaxis2)
  this%real_time_datamax = maxval(tmpaxis3)

  ! Setup nuclide lookup tables after file read
  num_partitions = (this%num_start_times*this%num_powers)

  cur_inventory => CritInventoryLookupCreate()
  do i = 1, this%num_species

    new_inventory  => CritInventoryLookupCreate()

    new_inventory%name = tmpname2(i)
    new_inventory%lookup => LookupTableCreateGeneral(THREE_INTEGER)
    new_inventory%lookup%dims(1) = this%num_start_times
    new_inventory%lookup%dims(2) = this%num_powers
    new_inventory%lookup%dims(3) = this%num_real_times
    new_inventory%lookup%axis3%num_partitions = num_partitions

    allocate(new_inventory%lookup%axis1%values(this%num_start_times))
    allocate(new_inventory%lookup%axis2%values(this%num_powers))
    allocate(new_inventory%lookup%axis3%values(arr_size))
    allocate(new_inventory%lookup%data(arr_size))

    if (Initialized(this%total_points)) then
      allocate(new_inventory%lookup%axis3%bounds(num_partitions))
      allocate(new_inventory%lookup%axis3%partition(num_partitions))
    endif

    if (Initialized(mode)) then
      new_inventory%lookup%mode = mode
    endif

    new_inventory%lookup%axis1%values => tmpaxis1

    new_inventory%lookup%axis2%values => tmpaxis2

    new_inventory%lookup%axis3%values => tmpaxis3

    if (Initialized(this%total_points)) then
      call CritInventoryRealTimeSections(new_inventory%lookup,filename,option)
    endif

    call CritInventoryCheckDuplicates(new_inventory%lookup,filename,option)

    if (size(tmpdata2(i)%data) /= arr_size) then
      write(str1,*) i
      write(str2,*) size(tmpdata2(i)%data)
      write(str3,*) arr_size
      option%io_buffer = 'Number of inventories listed for nuclide #' &
                       // trim(adjustl(str1)) // ' (' &
                       // trim(adjustl(str2)) // ') ' &
                       //'does not match length of REAL_TIME array (' &
                       // trim(adjustl(str3)) //') in "'// trim(filename) //'".'
      call PrintErrMsg(option)
    endif

    new_inventory%lookup%data = tmpdata2(i)%data * data_units_conversion

    if (Initialized(this%total_points)) then
      write(str1,*) i
      error_string = trim(filename) //' for nuclide #' &
        // trim(adjustl(str1))
      call CritInventoryDataSections(new_inventory%lookup, &
                                     error_string,option)
    endif

    ! Add lookup table to linked list
    if (associated(this%radionuclide_table)) then
      cur_inventory => this%radionuclide_table
      do
        if (.not. associated(cur_inventory%next)) exit
        cur_inventory => cur_inventory%next
      enddo
      cur_inventory%next => new_inventory
    else
      this%radionuclide_table => new_inventory
    endif

    nullify(new_inventory)

  enddo
  nullify(cur_inventory)

  if (associated(tmpaxis1)) nullify(tmpaxis1)
  if (associated(tmpaxis2)) nullify(tmpaxis2)
  if (associated(tmpaxis3)) nullify(tmpaxis3)
  if (associated(tmpdata))  nullify(tmpdata)
  if (allocated(tmpdata1)) deallocate(tmpdata1)
  if (allocated(tmpdata2)) deallocate(tmpdata2)

end subroutine CritInventoryRead

! ************************************************************************** !

subroutine CritInventoryRealTimeSections(this,string,option)
  !
  ! Partition of monontically increasing real time values (axis3) for
  !   criticality inventory lookup tables
  !
  ! Author: Alex Salazar III
  ! Date: 02/18/2022
  !
  use Option_module
  !
  implicit none
  ! ----------------------------------
  class(lookup_table_general_type), pointer :: this ! lookup table
  character(len=MAXSTRINGLENGTH), intent(in) :: string
  class (option_type), intent(inout) :: option
  ! ----------------------------------
  PetscReal, pointer :: array(:)    ! array of real times
  PetscInt  :: i, j, k, l
  PetscInt  :: sz
  PetscInt  :: bnd1, bnd2
  PetscReal :: tmp1, tmp2
  ! ----------------------------------

  ! This subroutine only operates upon axis3
  if (.not. associated(this%axis3)) return

  ! This subroutine is not needed if the axis3 does not need to be partitioned
  if (this%axis3%num_partitions <= 0) return

  ! Allocate axis3 objects
  if (.not. allocated(this%axis3%bounds)) then
    allocate(this%axis3%bounds(this%axis3%num_partitions))
  endif

  if (.not. allocated(this%axis3%partition)) then
    allocate(this%axis3%partition(this%axis3%num_partitions))
  endif

  if (associated(this%axis3%values)) then
    array => this%axis3%values
  else
    option%io_buffer = 'Values for REAL_TIME (axis3) were not associated in '&
                     // trim(string) // '.'
    call PrintErrMsg(option)
  endif

  ! Assuming monotonic real times, identify bounds for the different partitions
  j = 1
  tmp1 = 0.d0
  tmp2 = 0.d0
  do i = 1, size(array)
    tmp1 = array(i)
    if (i > 1) then
      if (tmp1 < tmp2) then
        this%axis3%bounds(j) = i - 1
        j = j + 1
    else
      if (j > this%axis3%num_partitions .and. i == size(array)) then
        option%io_buffer = 'Values for REAL_TIME must monotonically '&
                         //'increase for each inventory evaluation ' &
                         //'in the table for ' // trim(string) // '.'
        call PrintErrMsg(option)
      endif
    endif
    endif
    tmp2 = tmp1 ! store previous axis3 value
  enddo
  ! The last bound is the size of the unparitioned array
  this%axis3%bounds(j) = size(array)

  ! Partition the array
  j = 1
  sz = 0
  bnd1 = 0
  bnd2 = 0
  do i = 1, size(this%axis3%bounds)
    if (i == 1) then
      bnd1 = 1 ! start of unpartitioned array
    else
      bnd1 = this%axis3%bounds(i-1) + 1
    endif
    bnd2 = this%axis3%bounds(i)
    sz = abs(bnd2 - bnd1) + 1
    l = 1
    allocate(this%axis3%partition(j)%data(sz))
    do k = bnd1, bnd2
      this%axis3%partition(j)%data(l) = array(k)
      l = l + 1
    enddo
    j = j + 1
  enddo

  if (associated(array)) nullify(array)

end subroutine CritInventoryRealTimeSections

! ************************************************************************** !

subroutine CritInventoryDataSections(this,string,option)
  !
  ! Partition of nuclide inventory data for criticality inventory lookup tables
  !
  ! Author: Alex Salazar III
  ! Date: 02/18/2022
  !
  use Option_module
  !
  implicit none
  ! ----------------------------------
  class(lookup_table_general_type), pointer :: this ! lookup table
  character(len=MAXSTRINGLENGTH), intent(in) :: string
  class (option_type), intent(inout) :: option
  ! ----------------------------------
  PetscReal, pointer :: array(:)    ! array of real times
  PetscInt  :: i, j, k, l
  PetscInt  :: szlim, sz
  PetscInt  :: bnd1, bnd2
  ! ----------------------------------

  ! This subroutine is not needed if axis3 was not partitioned
  if (this%axis3%num_partitions <= 0) return

  ! Verify inventory data values are associated
  if (associated(this%data)) then
    array => this%data
    szlim = size(this%data)
  else
    option%io_buffer = 'Values for inventory data were not associated in ' &
                     // trim(string) // '.'
    call PrintErrMsg(option)
  endif

  ! Verify axis3 partition bounds exist
  if (.not. allocated(this%axis3%bounds)) then
    option%io_buffer = 'Boundaries for axis3 were not defined in ' &
                      // trim(string) // '.'
    call PrintErrMsg(option)
  endif

  ! Allocate data partition based on axis3 partitions if needed
  if (.not. allocated(this%partition)) then
    allocate(this%partition(this%axis3%num_partitions))
  endif

  if (szlim /= maxval(this%axis3%bounds)) then
    option%io_buffer = 'Array length mismatch between axis3 values and data ' &
                     //'in ' // trim(string) // '.'
    call PrintErrMsg(option)
  endif

  ! Partition the array
  j = 1
  sz = 0
  bnd1 = 0
  bnd2 = 0
  do i = 1, size(this%axis3%bounds)
    if (i == 1) then
      bnd1 = 1 ! start of unpartitioned array
    else
      bnd1 = this%axis3%bounds(i-1) + 1
    endif
    bnd2 = this%axis3%bounds(i)
    sz = abs(bnd2 - bnd1) + 1
    l = 1
    allocate(this%partition(j)%data(sz))
    do k = bnd1, bnd2
      this%partition(j)%data(l) = array(k)
      l = l + 1
    enddo
    j = j + 1
  enddo

  if (associated(array)) nullify(array)

end subroutine CritInventoryDataSections

! ************************************************************************** !

subroutine CritInventoryCheckDuplicates(this,string,option)
  !
  ! Checks for duplicate entries within the lookup table axes
  !
  ! Author: Alex Salazar III
  ! Date: 03/01/2022
  !
  use Option_module
  !
  implicit none
  ! ----------------------------------
  class(lookup_table_general_type), pointer :: this ! lookup table
  character(len=MAXSTRINGLENGTH), intent(in) :: string
  type(option_type), intent(inout) :: option
  ! ----------------------------------
  PetscReal, pointer :: values(:) ! values under inspection
  PetscInt :: i, j, k ! iterators
  PetscReal :: ref1, ref2 ! comparison
  PetscInt :: kstart, kend ! start and end of axis3 array to interpolate
  PetscInt :: nk ! number of lists expected in axis3
  character(len=MAXSTRINGLENGTH) :: sref1, sref2
  ! ----------------------------------

  ! Check axis1 for duplicates
  if (associated(this%axis1)) then
    values => this%axis1%values
    do i = 1, size(values)
      ref1 = values(i)
      do j = 1, size(values)
        if (i == j) cycle
        ref2 = values(j)
        if (ref1 == ref2) then
          write(sref1,'(es12.5)') ref1
          option%io_buffer = 'Duplicate entry (' // trim(adjustl(sref1)) &
                           //') detected in axis1 for "' &
                           // trim(string) // '".'
          call PrintErrMsg(option)
        end if
      enddo
    end do
  endif

  ! Check axis2 for duplicates
  if (associated(this%axis2)) then
    values => this%axis2%values
    do i = 1, size(values)
      ref1 = values(i)
      do j = 1, size(values)
        if (i == j) cycle
        ref2 = values(j)
        if (ref1 == ref2) then
          write(sref1,'(es12.5)') ref1
          option%io_buffer = 'Duplicate entry (' // trim(adjustl(sref1)) &
                           //') detected in axis2 for "' &
                           // trim(string) // '".'
          call PrintErrMsg(option)
        end if
      enddo
    end do
  endif

  ! Check axis3 for duplicates
  if (associated(this%axis3)) then

    if (allocated(this%axis3%partition)) then
      ! ---> axis3 is has defined partitions (non-rectangular)
      do k = 1, size(this%axis3%partition)
        values => this%axis3%partition(k)%data
        do i = 1, size(values)
          ref1 = values(i)
          do j = 1, size(values)
            if (i == j) cycle
            ref2 = values(j)
            if (ref1 == ref2) then
              write(sref1,'(es12.5)') ref1
              write(sref2,'(i3)') k
              option%io_buffer = 'Duplicate entry (' // trim(adjustl(sref1)) &
                               //') detected in partition ' &
                               // trim(adjustl(sref2)) //' of axis3 for "' &
                               // trim(string) // '".'
              call PrintErrMsg(option)
            end if
          enddo
        end do
      enddo

    else
      ! ---> axis3 is described by the dim(3) value (rectangular)
      nk = size(this%axis3%values)/this%dims(3)
      kstart = 0
      kend = 0
      do k = 1, nk
        kstart = (k - 1)*this%dims(3) + 1
        kend = k*this%dims(3)

        values => this%axis3%values(kstart:kend)
        do i = 1, size(values)
          ref1 = values(i)
          do j = 1, size(values)
            if (i == j) cycle
            ref2 = values(j)
            if (ref1 == ref2) then
              write(sref1,'(es12.5)') ref1
              write(sref2,'(i3)') k
              option%io_buffer = 'Duplicate entry (' // trim(adjustl(sref1)) &
                               //') detected in dataset ' &
                               // trim(adjustl(sref2)) //' of axis3 for "' &
                               // trim(string) // '".'
              call PrintErrMsg(option)
            end if
          enddo
        end do
      enddo
    endif
  endif

end subroutine CritInventoryCheckDuplicates

! ************************************************************************** !

subroutine CritInventoryUseLog10(inventory,option)
  !
  ! Transform radionuclide inventory lookup tables to use log10(time) basis
  !
  ! Author: Alex Salazar III
  ! Date: 08/17/2022
  !
  use Option_module
  !
  implicit none
  ! ----------------------------------
  class(crit_inventory_type), pointer :: inventory ! criticality inventory data
  type(option_type) :: option
  ! ----------------------------------
  type(crit_inventory_lookup_type), pointer :: inv ! radionuclide inventory
  PetscReal :: zero_sub ! substitute value for instances of zero
  PetscInt :: i, psize
  PetscBool :: modified
  ! ----------------------------------

  ! axis3 is modified (rectangular array)
  modified = PETSC_FALSE

  ! loop through linked list of lookup tables and modify axis 3 values
  inv => inventory%radionuclide_table
  do
    if (.not. associated(inv)) exit
    inv%use_log10_time = PETSC_TRUE
    zero_sub = inv%log10_time_zero_sub
    if (allocated(inv%lookup%axis3%partition)) then
      ! time axis has defined partitions (non-rectangular)
      psize = size(inv%lookup%axis3%partition)
      do i = 1, psize
        call CritInventoryLog10Array(inv%lookup%axis3%partition(i)%data, &
                                     zero_sub,option)
      enddo
    else
      ! time axis is described by the dim(3) value (rectangular)
      if (.not. modified) then
        call CritInventoryLog10Array(inv%lookup%axis3%values,zero_sub,option)
        modified = PETSC_TRUE
      endif
    endif
    inv => inv%next
  enddo
  nullify(inv)

end subroutine CritInventoryUseLog10

! ************************************************************************** !

subroutine CritInventoryLog10Array(array1,zero_sub,option)
  !
  ! Modify array to use log10 of original values
  !
  ! Author: Alex Salazar III
  ! Date: 08/16/2022
  !
  use Option_module
  !
  implicit none
  ! ----------------------------------
  PetscReal, pointer :: array1(:)
  type(option_type) :: option
  PetscReal :: zero_sub ! substitute value for instances of zero
  ! ----------------------------------
  PetscReal, pointer :: array2(:)
  PetscInt :: i, asize
  PetscReal :: tval ! values from array to be modified
  ! ----------------------------------

  ! allocate dummy array for log10 transformation
  asize = size(array1)
  allocate(array2(asize))

  ! populate dummy array with log10 values of original array
  do i = 1, asize
    tval = array1(i)
    ! replace zeros in original array with a substitute value
    if (tval <= 0) tval = zero_sub
    array2(i) = log10(tval)
  enddo

  ! replace original array with the log10 array
  array1(1:asize) = array2(1:asize)
  deallocate(array2)
  nullify(array2)

end subroutine CritInventoryLog10Array

! ************************************************************************** !

subroutine CritInventoryCheckZeroSub(zero_sub,inventory,option)
  !
  ! For log10 interpolation, make sure zero substitute is below minimum nonzero
  !   value in the time arrays
  !
  ! Author: Alex Salazar III
  ! Date: 08/18/2022
  !
  use Option_module
  !
  implicit none
  ! ----------------------------------
  PetscReal :: zero_sub ! substitute value for instances of zero
  class(crit_inventory_type), pointer :: inventory ! criticality inventory data
  type(option_type) :: option
  ! ----------------------------------
  type(crit_inventory_lookup_type), pointer :: inv ! radionuclide inventory
  PetscInt :: i, psize
  PetscReal :: min_nonzero ! smallest nonzero number in array
  character(len=MAXSTRINGLENGTH) :: word1, word2
  PetscBool :: checked
  ! ----------------------------------

  ! loop through linked list of lookup tables and modify axis 3 values
  checked = PETSC_FALSE
  inv => inventory%radionuclide_table
  do
    if (.not. associated(inv)) exit
    inv%log10_time_zero_sub = zero_sub
    if (allocated(inv%lookup%axis3%partition)) then
      ! time axis has defined partitions (non-rectangular)
      psize = size(inv%lookup%axis3%partition)
      do i = 1, psize
        min_nonzero = &
          CritInventoryMinNonzeroTime(inv%lookup%axis3%partition(i)%data,option)
        if (zero_sub >= min_nonzero) then
          write(word1,'(es11.5)') zero_sub
          write(word2,'(es11.5)') min_nonzero
          option%io_buffer = 'For log10 time interpolation, zero substitute (' &
                           // trim(word1) &
                           //') must remain below minimum nonzero value (' &
                           // trim(word2) //') in REAL_TIME array in file "' &
                           // trim(inventory%file_name) // '."'
          call PrintErrMsg(option)
        endif
      enddo
    elseif (.not. checked) then
      ! time axis is described by the dim(3) value (rectangular)
      min_nonzero = CritInventoryMinNonzeroTime(inv%lookup%axis3%values,option)
      if (zero_sub >= min_nonzero) then
        write(word1,'(es11.5)') zero_sub
        write(word2,'(es11.5)') min_nonzero
        option%io_buffer = 'For log10 time interpolation, zero substitute (' &
                         // trim(word1) &
                         //') must remain below minimum nonzero value (' &
                         // trim(word2) //') in REAL_TIME array in file "' &
                         // trim(inventory%file_name) // '."'
        call PrintErrMsg(option)
      endif
      checked = PETSC_TRUE
    endif
    inv => inv%next
  enddo
  nullify(inv)

end subroutine CritInventoryCheckZeroSub

! ************************************************************************** !

function CritInventoryMinNonzeroTime(array1,option)
  !
  ! Find minimum nonzero time
  !
  ! Author: Alex Salazar III
  ! Date: 08/18/2022
  !
  use Option_module
  !
  implicit none
  ! ----------------------------------
  PetscReal :: CritInventoryMinNonzeroTime
  PetscReal, pointer :: array1(:)
  type(option_type) :: option
  ! ----------------------------------
  PetscReal, allocatable :: array2(:)
  PetscInt :: i, j, asize1, asize2
  PetscReal :: val1
  ! ----------------------------------

  asize1 = size(array1)
  val1 = 0.0d0

  ! check number of non-zero entries
  asize2 = 0
  do i = 1, asize1
    val1 = array1(i)
    if (val1 > 0.0d0) asize2 = asize2 + 1
  enddo

  ! group nonzero values
  allocate(array2(asize2))
  j = 1
  do i = 1, asize1
    val1 = array1(i)
    if (val1 > 0.0d0) then
      array2(j) = val1
      j = j + 1
    endif
  enddo

  ! get smallest nonzero number
  CritInventoryMinNonzeroTime = minval(array2)
  deallocate(array2)

end function CritInventoryMinNonzeroTime

! ************************************************************************** !

function CritHeatEvaluate(this,start_time,temperature)
  !
  ! Author: Alex Salazar III
  ! Date: 05/12/2021
  !

  implicit none

  class(crit_heat_type) :: this
  PetscReal :: start_time
  PetscReal :: temperature

  PetscReal :: CritHeatEvaluate

  CritHeatEvaluate = this%lookup_table%Sample(start_time,temperature)

end function CritHeatEvaluate

! ************************************************************************** !

function CritInventoryEvaluate(this,start_time,power,time)
  !
  ! Evaluate radionuclide mass fraction from lookup table
  !
  ! Author: Alex Salazar III
  ! Date: 02/21/2022
  !
  implicit none
  ! ----------------------------------
  class(crit_inventory_lookup_type) :: this
  PetscReal :: start_time
  PetscReal :: power
  PetscReal :: time
  PetscReal :: CritInventoryEvaluate
  ! ----------------------------------
  PetscReal :: t ! time value passed to lookup table
  PetscReal :: zero_sub ! substitute time for t=0 (log10 interpolation)
  ! ----------------------------------

  ! time value passed to lookup table
  if (this%use_log10_time) then
    if (time <= 0.d0) then
      zero_sub = this%log10_time_zero_sub
      t = log10(zero_sub) ! substitute for zero
    else
      t = log10(time)
    endif
  else
    t = time
  endif

  CritInventoryEvaluate = this%lookup%Sample(start_time,power,t)

end function CritInventoryEvaluate

! ************************************************************************** !

function CritHeatCreate()
  !
  ! Author: Alex Salazar III
  ! Date: 05/12/2021
  !

  implicit none

  class(crit_heat_type), pointer :: CritHeatCreate
  class(crit_heat_type), pointer :: ch

  allocate(ch)
  nullify(ch%next)
  nullify(ch%lookup_table)

  ch%file_name = ''
  ch%num_start_times = UNINITIALIZED_INTEGER
  ch%num_values_per_start_time = UNINITIALIZED_INTEGER
  ch%start_time_datamax = UNINITIALIZED_DOUBLE
  ch%temp_datamax = UNINITIALIZED_DOUBLE
  ch%power_datamax = UNINITIALIZED_DOUBLE

  CritHeatCreate => ch

end function CritHeatCreate

! ************************************************************************** !

function CritInventoryCreate()
  !
  ! Author: Alex Salazar III
  ! Date: 02/16/2022
  !

  implicit none

  class(crit_inventory_type), pointer :: CritInventoryCreate
  class(crit_inventory_type), pointer :: ci

  allocate(ci)
  nullify(ci%next)
  nullify(ci%radionuclide_table)

  ci%file_name = ''
  ci%total_points    = UNINITIALIZED_INTEGER
  ci%num_start_times = UNINITIALIZED_INTEGER
  ci%num_powers      = UNINITIALIZED_INTEGER
  ci%num_real_times  = UNINITIALIZED_INTEGER
  ci%num_species     = UNINITIALIZED_INTEGER
  ci%start_time_datamax = UNINITIALIZED_DOUBLE
  ci%power_datamax      = UNINITIALIZED_DOUBLE
  ci%real_time_datamax  = UNINITIALIZED_DOUBLE
  ci%switch_implicit = PETSC_FALSE
  ci%allow_implicit = PETSC_FALSE
  ci%allow_extrap   = PETSC_FALSE
  ci%continue_lookup = PETSC_FALSE

  CritInventoryCreate => ci

end function CritInventoryCreate

! ************************************************************************** !

function CritInventoryLookupCreate()
  !
  ! Author: Alex Salazar III
  ! Date: 04/19/2022
  !

  implicit none

  class(crit_inventory_lookup_type), pointer :: CritInventoryLookupCreate
  class(crit_inventory_lookup_type), pointer :: cl

  allocate(cl)
  nullify(cl%next)
  nullify(cl%lookup)

  cl%name = ''
  cl%use_log10_time = PETSC_FALSE
  cl%log10_time_zero_sub = 1.0d-20

  CritInventoryLookupCreate => cl

end function CritInventoryLookupCreate

! ************************************************************************** !

function relu(x)
! Rectified linear unit
  implicit none
  PetscReal :: relu
  PetscReal, intent(in) :: x
  if (x >= 0.0d0) then
    relu = x
  else
    relu = 0.0d0
  endif
end function relu

! ************************************************************************** !

function dose_rate(years_time, decay_time, burnup)
! Computes the dose rate at the fuel surface
  implicit none
  PetscReal :: dose_rate
  PetscReal, intent(in) :: years_time, decay_time, burnup
  PetscReal :: f1, f2, f3, f4, f5
  PetscReal :: AOF, rad0a

  AOF = years_time + decay_time

  f2 = log(AOF)
  f1 = f2**2.0d0
  f3 = 1.0d0/f2
  f4 = f2/AOF
  f5 = exp(burnup/25.26892627636246d0)

  rad0a = -206.0634818750711d0   - 0.7631591788870090d0*f1 &
        + 20.97112373957833d0*f2 + 678.8463343193430d0*f3 &
        - 506.7149017370657d0*f4 + 0.1555448893425319d0*f5
  dose_rate = max(exp(rad0a),5.0d-3)
end function dose_rate

! ************************************************************************** !
subroutine AMP_ann_surrogate_step(this, sTme, current_temp_C)

  implicit none
  class(wf_mechanism_fmdm_surrogate_type) :: this
  PetscReal, intent(in) :: sTme
  PetscReal, intent(in) :: current_temp_C
  ! constants
  PetscInt, parameter :: num_features = 6 ! number of inputs to ANN
  PetscInt, parameter :: N = 64 ! number of nodes per hidden layer
  PetscReal, parameter :: UO2_molar_mass = 270.0d0 ! g/mol

  ! local variables
  PetscReal :: yTme
  PetscInt :: i
  ! features
  PetscReal :: f(6)
  ! hidden layer nodes values
  PetscReal :: h1(64), h2(64)

  yTme = sTme/60.0d0/60.0d0/24.0d0/DAYS_PER_YEAR

  ! features
  f(1) = current_temp_C + 273.15d0
  f(2) = log10(this%concentration(1)) ! Env_CO3_2n
  f(3) = log10(this%concentration(2)) ! Env_O2
  f(4) = log10(this%concentration(3)) ! Env_Fe_2p
  f(5) = log10(this%concentration(4)) ! Env_H2
  f(6) = log10(dose_rate(yTme,this%decay_time,this%burnup))

  ! standardize
  do i = 1,num_features
    f(i) = (f(i) - this%scaler_offsets(i))/this%scaler_scales(i)
  enddo

  ! Input - Hidden Layer 1
  do i = 1,N
    h1(i) = relu(dot_product(f, this%input_hidden1_weights(:,i)) &
          + this%input_hidden1_bias(i))
  enddo

  ! Hidden Layer 1 - Hidden Layer 2
  do i = 1,N
    h2(i) = relu(dot_product(h1, this%hidden1_hidden2_weights(:,i)) &
          + this%hidden1_hidden2_bias(i))
  enddo

  ! Hidden Layer 2 - Output
  this%dissolution_rate = 10**(dot_product(h2, this%hidden2_output_weights) &
                        + this%hidden2_output_bias)*UO2_molar_mass

end subroutine AMP_ann_surrogate_step

! ************************************************************************** !

subroutine ANNReadH5File(this, option)

  use hdf5
  use HDF5_Aux_module

  implicit none

  type(option_type) :: option
  class(wf_mechanism_fmdm_surrogate_type) :: this

  character(len=MAXSTRINGLENGTH) :: h5_name = 'fmdm_ann_coeffs.h5'
  character(len=MAXSTRINGLENGTH) :: group_name = '/'
  character(len=MAXSTRINGLENGTH) :: dataset_name

  integer(HID_T) :: prop_id
  integer(HID_T) :: file_id
  integer(HID_T) :: group_id
  integer(HID_T) :: dataset_id

  integer(HSIZE_T), allocatable :: dims_h5(:)

  PetscMPIInt :: hdf5_err

  call h5open_f(hdf5_err)
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
  call HDF5OpenFileReadOnly(h5_name,file_id,prop_id,'',option)
  call HDF5GroupOpen(file_id,group_name,group_id,option)

  dataset_name = 'input_hidden1_weights'
  call ANNGetH5DatasetInfo(group_id,option,h5_name,dataset_name,dataset_id, &
       dims_h5)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, this%input_hidden1_weights, dims_h5, &
                 hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)
  deallocate(dims_h5)

  dataset_name = 'input_hidden1_bias'
  call ANNGetH5DatasetInfo(group_id,option,h5_name,dataset_name,dataset_id, &
       dims_h5)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, this%input_hidden1_bias,dims_h5, &
                 hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)
  deallocate(dims_h5)

  dataset_name = 'hidden1_hidden2_weights'
  call ANNGetH5DatasetInfo(group_id,option,h5_name,dataset_name,dataset_id, &
       dims_h5)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, this%hidden1_hidden2_weights, dims_h5, &
                 hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)
  deallocate(dims_h5)

  dataset_name = 'hidden1_hidden2_bias'
  call ANNGetH5DatasetInfo(group_id,option,h5_name,dataset_name,dataset_id, &
       dims_h5)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, this%hidden1_hidden2_bias, dims_h5, &
                 hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)
  deallocate(dims_h5)

  dataset_name = 'hidden2_output_weights'
  call ANNGetH5DatasetInfo(group_id,option,h5_name,dataset_name,dataset_id, &
       dims_h5)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%hidden2_output_weights,dims_h5, &
                   hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)
  deallocate(dims_h5)

  dataset_name = 'hidden2_output_bias'
  call ANNGetH5DatasetInfo(group_id,option,h5_name,dataset_name,dataset_id, &
       dims_h5)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%hidden2_output_bias,dims_h5, &
                 hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)
  deallocate(dims_h5)

  dataset_name = 'scaler_offsets'
  call ANNGetH5DatasetInfo(group_id,option,h5_name,dataset_name,dataset_id, &
       dims_h5)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%scaler_offsets,dims_h5, &
                 hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)
  deallocate(dims_h5)

  dataset_name = 'scaler_scales'
  call ANNGetH5DatasetInfo(group_id,option,h5_name,dataset_name,dataset_id, &
       dims_h5)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%scaler_scales,dims_h5, &
                 hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)
  deallocate(dims_h5)

  call h5gclose_f(group_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5pclose_f(prop_id,hdf5_err)

end subroutine ANNReadH5File

! ************************************************************************** !

subroutine ANNGetH5DatasetInfo(group_id,option,h5_name,dataset_name,dataset_id,&
                               dims_h5)

  use hdf5

  implicit none

  type(option_type) :: option

  integer(HID_T) :: group_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: file_space_id

  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)

  PetscInt :: ndims_h5

  character(len=MAXSTRINGLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: h5_name
  PetscMPIInt :: hdf5_err

  call h5dopen_f(group_id,dataset_name,dataset_id,hdf5_err)

  if (hdf5_err < 0) then
    option%io_buffer = 'A dataset named "' // trim(dataset_name) // '" not found in HDF5 file "' // &
    trim(h5_name) // '".'
    call PrintErrMsg(option)
  endif

  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)
  call h5sget_simple_extent_ndims_f(file_space_id,ndims_h5,hdf5_err)

  allocate(dims_h5(ndims_h5))
  allocate(max_dims_h5(ndims_h5))

  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)

  deallocate(max_dims_h5)

end subroutine ANNGetH5DatasetInfo

! ************************************************************************** !

subroutine KnnrInit(this,option)

  implicit none

  type(option_type) :: option
  class(wf_mechanism_fmdm_surrogate_type) :: this
  PetscInt :: i_d, d
  PetscInt :: data_array_shape(2)

  call KnnrReadH5File(this, option)

  data_array_shape = shape(this%table_data)
  ! Quantities of Interest (QoI) is not part of the search query of the search query.
  this%num_qoi = data_array_shape(1)-1
  d = data_array_shape(2)

  allocate(this%knnr_array(this%num_qoi,d))

  do i_d = 1, this%num_qoi
    this%knnr_array(i_d,:) = this%table_data(i_d,:)
  end do

  this%tree => KdtreeCreate()

  call KdtreeConstruct(this%tree,this%knnr_array,sort=PETSC_FALSE,rearrange=PETSC_FALSE)

end subroutine KnnrInit

! ************************************************************************** !

subroutine KnnrQuery(this,sTme,current_temp_C)

  implicit none

  class(wf_mechanism_fmdm_surrogate_type) :: this

  PetscReal :: current_temp_C
  PetscReal :: decay_time
  PetscReal :: burnup
  PetscReal :: sTme
  PetscInt :: nn

  ! features
  PetscReal :: f(4)
  PetscReal :: yTme

  PetscReal :: qoi_ave
  PetscReal, parameter :: UO2_molar_mass = 270.0d0 !g/mol

  type(kdtree_result), allocatable :: knnr_results(:)

  decay_time = this%decay_time
  burnup = this%burnup
  nn = this%num_nearest_neighbor

  yTme = sTme/60.0d0/60.0d0/24.0d0/DAYS_PER_YEAR

  f(1) = log10(current_temp_C + 273.15d0)
  f(2) = log10(this%concentration(1)) ! Env_CO3_2n
  f(3) = log10(this%concentration(4)) ! Env_H2
  f(4) = log10(dose_rate(yTme,decay_time,burnup))

  allocate(knnr_results(nn))

  call kdtreeNNearest(tp=this%tree,qv=f,nn=nn,results=knnr_results)

  call KnnrInverseDistance(knnr_results,nn,this%table_data,this%num_qoi,this%knnr_eps,qoi_ave)

  this%dissolution_rate = (qoi_ave) * UO2_molar_mass !convert units

end subroutine KnnrQuery

! ************************************************************************** !

subroutine KnnrReadH5File(this, option)

  use hdf5
  use HDF5_Aux_module

  implicit none

  type(option_type) :: option
  class(wf_mechanism_fmdm_surrogate_type) :: this

  character(len=MAXSTRINGLENGTH) :: h5_name = 'FMDM_knnr_data.h5'
  character(len=MAXSTRINGLENGTH) :: group_name = '/'
  character(len=MAXSTRINGLENGTH) :: dataset_name

  integer(HID_T) :: prop_id
  integer(HID_T) :: file_id
  integer(HID_T) :: group_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: file_space_id

  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)

  PetscInt :: ndims_h5

  PetscMPIInt :: hdf5_err

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)

  call HDF5OpenFileReadOnly(h5_name,file_id,prop_id,'',option)

  call h5pclose_f(prop_id,hdf5_err)

  !hdf5groupopen
  call HDF5GroupOpen(file_id,group_name,group_id,option)

  !Get Nearest Neighbors
  call KnnrGetNearestNeighbors(this,group_id,h5_name,option)

  !Read features
  dataset_name = 'Temp'
  call h5dopen_f(group_id,dataset_name,dataset_id,hdf5_err)

  if (hdf5_err < 0) then
    option%io_buffer = 'A dataset named "' // trim(dataset_name) // '" not found in HDF5 file "' // &
    trim(h5_name) // '".'
    call PrintErrMsg(option)
  endif

  ! get dataspace ID
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  call h5sget_simple_extent_ndims_f(file_space_id,ndims_h5,hdf5_err)

  allocate(dims_h5(ndims_h5))
  allocate(max_dims_h5(ndims_h5))

  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)

  allocate(this%table_data(5,dims_h5(1)))


  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, this%table_data(1,:), dims_h5, &
                   hdf5_err)

  call h5dclose_f(dataset_id,hdf5_err)

  dataset_name = 'Env_CO3_2n'
  call KnnrReadH5Dataset(this,group_id,dims_h5,option,h5_name,dataset_name,2)
  dataset_name = 'Env_H2'
  call KnnrReadH5Dataset(this,group_id,dims_h5,option,h5_name,dataset_name,3)
  dataset_name = 'Dose Rate d0'
  call KnnrReadH5Dataset(this,group_id,dims_h5,option,h5_name,dataset_name,4)
  dataset_name = 'UO2 Surface Flux'
  call KnnrReadH5Dataset(this,group_id,dims_h5,option,h5_name,dataset_name,5)

  deallocate(dims_h5)
  deallocate(max_dims_h5)


  call h5gclose_f(group_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)

  this%table_data(1,:) = log10(this%table_data(1,:))
  this%table_data(2,:) = log10(this%table_data(2,:))
  this%table_data(3,:) = log10(this%table_data(3,:))
  this%table_data(4,:) = log10(this%table_data(4,:))


end subroutine KnnrReadH5File

! ************************************************************************** !

subroutine KnnrGetNearestNeighbors(this,group_id,h5_name,option)

  use hdf5

  implicit none

  type(option_type) :: option
  class(wf_mechanism_fmdm_surrogate_type) :: this

  integer(HID_T) :: group_id
  integer(HID_T) :: dataset_id

  integer(HSIZE_T) :: dims_h5(1) = 1

  character(len=MAXSTRINGLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: h5_name
  PetscMPIInt :: hdf5_err

  dataset_name = 'Nearest Neighbors Num'

  call h5dopen_f(group_id,dataset_name,dataset_id,hdf5_err)

  if (hdf5_err < 0) then
    option%io_buffer = 'A dataset named "' // trim(dataset_name) // '" not found in HDF5 file "' // &
    trim(h5_name) // '".'
    call PrintErrMsg(option)
  else
     call h5dread_f(dataset_id,H5T_NATIVE_INTEGER, this%num_nearest_neighbor, dims_h5, &
       hdf5_err)

     call h5dclose_f(dataset_id,hdf5_err)
  endif

end subroutine KnnrGetNearestNeighbors

! ************************************************************************** !

subroutine KnnrReadH5Dataset(this,group_id,dims_h5,option,h5_name,dataset_name,i)

  use hdf5

  implicit none

  type(option_type) :: option
  class(wf_mechanism_fmdm_surrogate_type) :: this

  integer(HID_T) :: group_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: file_space_id

  integer(HSIZE_T),allocatable :: dims_h5(:)

  PetscInt :: i

  character(len=MAXSTRINGLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: h5_name
  PetscMPIInt :: hdf5_err

  call h5dopen_f(group_id,dataset_name,dataset_id,hdf5_err)

  if (hdf5_err < 0) then
    option%io_buffer = 'A dataset named "' // trim(dataset_name) // '" not found in HDF5 file "' // &
    trim(h5_name) // '".'
    call PrintErrMsg(option)
  endif

  ! get dataspace ID
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, this%table_data(i,:), dims_h5, &
       hdf5_err)

  call h5dclose_f(dataset_id,hdf5_err)


end subroutine KnnrReadH5Dataset

! ************************************************************************** !

subroutine KnnrInverseDistance(knnr_results,nn,table_data,n,eps,qoi_ave)

  implicit none

  PetscReal :: qoi_i, qoi_sum, qoi_ave, qoi_weights, weight, dis

  type(kdtree_result), allocatable :: knnr_results(:)
  PetscReal :: eps, table_data(:,:)
  PetscInt :: n

  type(kdtree_result) :: knnr_qoi
  PetscInt :: i_d,nn

  qoi_weights = 0.0
  qoi_sum = 0.0

  do i_d = 1,nn
    knnr_qoi = knnr_results(i_d)

    qoi_i = log10(table_data(n+1,knnr_qoi%idx))

    dis = knnr_qoi%dis

    if (abs(dis) <= eps) then
      qoi_weights = 1.0
      qoi_sum = qoi_i

      exit
    elseif (KnnrIsInfinite(abs(1/dis))) then
      qoi_weights = 1.0
      qoi_sum = qoi_i

      exit
    else

       weight = 1 / dis

       qoi_sum = qoi_sum + qoi_i * weight

       qoi_weights = qoi_weights + weight

    endif

  enddo

  qoi_ave = qoi_sum/qoi_weights
  qoi_ave = 10**(qoi_ave)

end subroutine KnnrInverseDistance

! ************************************************************************** !

function KnnrIsInfinite(value1)

  implicit none

  PetscBool :: KnnrIsInfinite
  PetscReal :: value1
  PetscReal :: infinity

  KnnrIsInfinite = PETSC_FALSE

  infinity = huge(0.0d0)

  if (value1 >= infinity) then
    KnnrIsInfinite = PETSC_TRUE
  endif

end function KnnrIsInfinite

! ************************************************************************** !

end module PM_Waste_Form_class

