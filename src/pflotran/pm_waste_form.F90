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
 
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none

  private

  PetscBool, public :: bypass_warning_message = PETSC_FALSE

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
! artificial neural network SURROGATE APPROXIMATION
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
    PetscReal, pointer :: outer_weights(:)
    PetscReal, pointer :: inner_weights(:,:)
    PetscReal, pointer :: scaler_means(:)
    PetscReal, pointer :: scaler_variances(:)
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
    PetscReal :: canister_vitality                     
    PetscReal :: canister_vitality_rate
    PetscReal :: eff_canister_vit_rate
    PetscReal :: breach_time                           
    PetscBool :: breached
    PetscReal :: decay_start_time                      
    character(len=MAXWORDLENGTH) :: mech_name
    class(wf_mechanism_base_type), pointer :: mechanism
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
    type(criticality_mediator_type), pointer :: criticality_mediator
    PetscBool :: print_mass_balance
    PetscBool :: implicit_solution
    PetscBool :: dataset_solution
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

! -------------------------------------------------------------------

! Stores variables relevant to criticality calculations
  type, public :: criticality_mechanism_type
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
    class(dataset_ascii_type), pointer :: rad_dataset
    class(dataset_ascii_type), pointer :: heat_dataset
    type(criticality_mechanism_type), pointer :: next
  end type criticality_mechanism_type

! -----------------------------------------------------

  ! Stores information regarding the criticality event
  type, public :: criticality_event_type
    character(len=MAXWORDLENGTH) :: mech_name
    PetscBool :: steady_state
    PetscReal :: crit_start
    PetscReal :: crit_end
    PetscBool :: crit_flag
  end type criticality_event_type

! -----------------------------------------------------

  ! Criticality process model object. Includes relevant variables for 
  ! criticality consequence calculations.
  type, public :: criticality_type
    type(criticality_event_type), pointer :: crit_event
    type(criticality_mechanism_type), pointer :: crit_mech
    type(region_type), pointer :: region
    type(criticality_type), pointer :: next
  end type criticality_type

! -----------------------------------------------------

  type, public :: criticality_mediator_type
    type(criticality_type), pointer :: criticality_list
    type(criticality_mechanism_type), pointer :: crit_mech_list
    class(data_mediator_vec_type), pointer :: data_mediator
    PetscInt :: total_num_cells
  end type criticality_mediator_type

! -----------------------------------------------------
  
  public :: PMWFCreate, &
            PMWFSetup, &
            PMWFMechanismGlassCreate, &
            PMWFMechanismDSNFCreate, &
            PMWFMechanismWIPPCreate, &
            PMWFMechanismCustomCreate, &
            PMWFMechanismFMDMCreate, &
            PMWFMechanismFMDMSurrogateCreate, &
            PMWFRadSpeciesCreate, &
            CriticalityMediatorCreate, &
            CriticalityMechCreate, &
            CriticalityCreate, &
            ReadCriticalityMech, &
            CriticalityCalc, &
            CriticalityInitializeRun, &
            AssignCritMech, &
            CriticalitySolve
  
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

function PMWFMechanismFMDMSurrogateCreate()
  ! 
  ! Creates the FMDM surrogate waste form mechanism package
  ! 
  ! Author: Tom Seidl
  ! Date: 03/05/2019

  implicit none

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

  allocate(surrfmdm%outer_weights(101))
  open(IUNIT_TEMP,file="ann_surrogate/outer_weights.txt")
  read(IUNIT_TEMP,*) surrfmdm%outer_weights
  close(IUNIT_TEMP)

  allocate(surrfmdm%inner_weights(7,100))
  open(IUNIT_TEMP,file="ann_surrogate/inner_weights.txt")
  read(IUNIT_TEMP,*) surrfmdm%inner_weights
  close(IUNIT_TEMP)

  allocate(surrfmdm%scaler_means(6))
  open(IUNIT_TEMP,file="ann_surrogate/means.txt")
  read(IUNIT_TEMP,*) surrfmdm%scaler_means
  close(IUNIT_TEMP)

  allocate(surrfmdm%scaler_variances(6))
  open(IUNIT_TEMP,file="ann_surrogate/vars.txt")
  read(IUNIT_TEMP,*) surrfmdm%scaler_variances
  close(IUNIT_TEMP)

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

function PMWFRadSpeciesCreate()
  ! 
  ! Creates a radioactive species in the waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/09/16

  implicit none

! LOCAL VARIABLES:
! ================
! PMWFRadSpeciesCreate (output): new radionuclide species object
! ---------------------------------------------- 
  type(rad_species_type) :: PMWFRadSpeciesCreate
! ----------------------------------------------

  PMWFRadSpeciesCreate%name = ''
  PMWFRadSpeciesCreate%daughter = ''
  PMWFRadSpeciesCreate%daugh_id = UNINITIALIZED_INTEGER
  PMWFRadSpeciesCreate%formula_weight = UNINITIALIZED_DOUBLE
  PMWFRadSpeciesCreate%decay_constant = UNINITIALIZED_DOUBLE
  PMWFRadSpeciesCreate%mass_fraction = UNINITIALIZED_DOUBLE
  PMWFRadSpeciesCreate%inst_release_fraction = UNINITIALIZED_DOUBLE
  PMWFRadSpeciesCreate%ispecies = UNINITIALIZED_INTEGER

end function PMWFRadSpeciesCreate

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
  wf%decay_start_time = 0.d0          ! sec (default value)
  nullify(wf%instantaneous_mass_rate) ! mol-rad/sec
  nullify(wf%cumulative_mass)         ! mol-rad
  nullify(wf%rad_mass_fraction)       ! g-rad/g-matrix
  nullify(wf%rad_concentration)       ! mol-rad/g-matrix
  nullify(wf%inst_release_amount)     ! mol-rad/g-matrix
  nullify(wf%mechanism)
  nullify(wf%next)
 !------- canister degradation model -----------------
  wf%canister_degradation_flag = PETSC_FALSE
  wf%breached = PETSC_FALSE
  wf%breach_time = UNINITIALIZED_DOUBLE
  wf%canister_vitality = 0.d0
  wf%canister_vitality_rate = UNINITIALIZED_DOUBLE
  wf%eff_canister_vit_rate = UNINITIALIZED_DOUBLE
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
  nullify(PMWFCreate%criticality_mediator)
  PMWFCreate%print_mass_balance = PETSC_FALSE
  PMWFCreate%implicit_solution = PETSC_FALSE
  PMWFCreate%dataset_solution = PETSC_FALSE
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
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
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
    !-------------------------------------
      case('DATASET_SOLUTION')
        this%dataset_solution = PETSC_TRUE
        cycle
    end select

    error_string = 'WASTE_FORM_GENERAL'
    call PMWFReadMechanism(this,input,option,word,error_string,found)
    if (found) cycle
    
    error_string = 'WASTE_FORM_GENERAL'
    call PMWFReadWasteForm(this,input,option,word,error_string,found)
    if (found) cycle
    
    error_string = 'WASTE_FORM_GENERAL'
    call ReadCriticalityMech(this%criticality_mediator,input,option, &
                             word,error_string,found)
    if (found) cycle
   
  enddo
  call InputPopBlock(input,option)
  
  ! Assign chosen mechanism to each criticality object
  if (associated(this%criticality_mediator)) then
    call AssignCritMech(this%criticality_mediator)
  endif
  
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
    ! error messaging: ----------------------------------------------
    if (.not.associated(cur_waste_form%mechanism)) then
      option%io_buffer = 'WASTE_FORM MECHANISM ' // &
                         trim(cur_waste_form%mech_name) // &
                         ' not found amoung given mechanism names.'
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
    
    cur_waste_form => cur_waste_form%next
  enddo
    
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
  type(rad_species_type), pointer :: temp_species_array(:)
  class(wf_mechanism_base_type), pointer :: new_mechanism, cur_mechanism
  PetscInt :: k, j
  PetscReal :: double
  PetscInt :: integer
! ----------------------------------------------------------------------

  error_string = trim(error_string) // ',MECHANISM'
  found = PETSC_TRUE
  added = PETSC_FALSE
  input%ierr = 0
  allocate(temp_species_array(50))
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
          new_mechanism => PMWFMechanismFMDMSurrogateCreate()
      !---------------------------------
        case('CUSTOM')
          error_string = trim(error_string) // ' CUSTOM'
          allocate(new_mechanism)
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
                  new_mechanism%decay_time,'day','DECAY_TIME', &
                  option)
              class default
                option%io_buffer = 'ERROR: DECAY_TIME cannot be &
                                   &specified for ' // trim(error_string)
                call PrintMsg(option)
                num_errors = num_errors + 1
            end select
        !--------------------------
          case('SPECIES')
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
              temp_species_array(k) = PMWFRadSpeciesCreate() 
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
  type(criticality_type), pointer :: cur_criticality, new_criticality
! ----------------------------------------------------------------------

  error_string = trim(error_string) // ',WASTE_FORM'
  found = PETSC_TRUE
  added = PETSC_FALSE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('WASTE_FORM')
      allocate(new_waste_form)
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
          case('CRITICALITY')
            new_criticality => CriticalityCreate()
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input, option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'keyword',error_string)
              call StringToUpper(word)
              select case (trim(word))
                case('MECH_NAME')
                  call InputReadCard(input,option,word)
                  call InputErrorMsg(input,option,'criticality mechanism assignment', &
                                     error_string)
                  call StringToUpper(word)
                  new_criticality%crit_event%mech_name= trim(word)
                case('CRIT_START')
                  call InputReadDouble(input,option,new_criticality% &
                                       crit_event%crit_start)
                  call InputErrorMsg(input,option,'CRIT_START',error_string)
                  call InputReadAndConvertUnits(input,new_criticality% &
                           crit_event%crit_start,'sec', &
                           trim(error_string)//',CRIT_START', &
                           option)
                case('CRIT_END')
                  call InputReadDouble(input,option,new_criticality% &
                                       crit_event%crit_end)
                  call InputErrorMsg(input,option,'CRIT_END',error_string)
                  call InputReadAndConvertUnits(input,new_criticality% &
                           crit_event%crit_end,'sec', &
                           trim(error_string)//',CRIT_END', &
                           option)
                case default
                  call InputKeywordUnrecognized(input,word,error_string,option)
              end select
            enddo      
            call InputPopBlock(input,option)
            if (.not. associated(this%criticality_mediator)) then
              this%criticality_mediator => CriticalityMediatorCreate()
              this%criticality_mediator%criticality_list => new_criticality
            else
              cur_criticality => this%criticality_mediator%criticality_list
              do
                if (.not. associated(cur_criticality)) exit
                if (.not. associated(cur_criticality%next)) then
                  cur_criticality%next => new_criticality
                endif
                cur_criticality => cur_criticality%next
              enddo
            endif
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
          (len(trim(new_waste_form%region_name)) == 0)) then
        option%io_buffer = 'ERROR: Either COORDINATE or REGION must be &
                           &specified for all waste forms.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Initialized(new_waste_form%coordinate%z) .and. &
          (len(trim(new_waste_form%region_name)) > 0)) then
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
  type(criticality_type), pointer :: cur_criticality
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
  if (associated(this%criticality_mediator)) then
    cur_criticality => this%criticality_mediator%criticality_list
  else
    nullify(cur_criticality)
  endif
  
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
        
        if (associated(cur_criticality)) then
          cur_criticality%region => new_region
          cur_criticality => cur_criticality%next
        endif
        
        allocate(cur_waste_form%scaling_factor(1))
        cur_waste_form%scaling_factor(1) = 1.d0
      endif
    else
      cur_region => region_list%first
      do
        if (.not.associated(cur_region)) exit
        if (StringCompare(cur_region%name,cur_waste_form%region_name)) then
          cur_waste_form%region => cur_region
          
          if (associated(cur_criticality)) then
            cur_criticality%region => cur_region
            cur_criticality => cur_criticality%next
          endif
          
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
  
  if (associated(this%criticality_mediator)) then
    cur_criticality => this%criticality_mediator%criticality_list
    do
      if (.not.associated(cur_criticality)) exit
      if (.not.associated(cur_criticality%region)) then
        nullify(this%criticality_mediator%criticality_list)
        exit
      endif
      cur_criticality => cur_criticality%next
    enddo
  endif
  
  
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

  use Material_Aux_class
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
  class(material_auxvar_type), pointer :: material_auxvars(:)
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
  
  allocate(ranks(option%mycommsize))
  
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
    call MPI_Allreduce(MPI_IN_PLACE,ranks,option%mycommsize,MPI_INTEGER, &
                       MPI_SUM,option%mycomm,ierr)
    newcomm_size = sum(ranks)
    allocate(cur_waste_form%rank_list(newcomm_size))
    j = 0
    do i = 1,option%mycommsize
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
      cur_waste_form%mechanism%rad_species_list(j)%ispecies = &
        GetPrimarySpeciesIDFromName( &
        cur_waste_form%mechanism%rad_species_list(j)%name, &
        this%realization%reaction,this%option)
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
  PetscInt :: num_waste_form_cells
  PetscInt :: num_species
  PetscInt :: size_of_vec
  PetscInt :: i, j, k
  PetscInt, allocatable :: species_indices_in_residual(:)
  PetscErrorCode :: ierr
! -------------------------------------------------------
  
  if (this%option%restart_flag .and. .not.this%skip_restart) then
      ! need to verify whether this is needed anymore
      call PMWFSetup(this)
  endif
  
  if (this%print_mass_balance) then
    call PMWFOutputHeader(this)
    call PMWFOutput(this)
  endif

  ! set up heat transfer
  
  if (associated(this%criticality_mediator)) then
    call CriticalityInitializeRun(this%criticality_mediator,this%realization, &
                                  this%option)
  endif

  ! set up mass transfer
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
  call VecCreateSeq(PETSC_COMM_SELF,size_of_vec, &
                    this%data_mediator%vec,ierr);CHKERRQ(ierr)
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
                       species_indices_in_residual, &
                       PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)
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
  use Material_Aux_class
  use Field_module
  use Option_module
  use Grid_module
  use Patch_module
  use Utility_module
  use Dataset_Ascii_class
  
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
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(criticality_type), pointer :: cur_criticality
  PetscReal :: dV
  PetscReal :: dt
  PetscReal :: avg_temp_local, avg_temp_global
  PetscInt :: local_id, ghosted_id
  PetscInt :: idof
  PetscInt :: i, k, p, g, d, f, j
  PetscInt :: num_species
  PetscErrorCode :: ierr
  PetscReal, allocatable :: Coeff(:)
  PetscReal, allocatable :: concentration_old(:)
  PetscReal :: inst_release_molality
  PetscReal, parameter :: conversion = 1.d0/(24.d0*3600.d0)
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
  class(dataset_ascii_type), pointer :: dataset

  global_auxvars => this%realization%patch%aux%Global%auxvars
  material_auxvars => this%realization%patch%aux%Material%auxvars
  field => this%realization%field
  option => this%option
  grid => this%realization%patch%grid
  dt = option%tran_dt
  
  call PMBasePrintHeader(this)

  ! zero entries from previous time step
  call VecZeroEntries(this%data_mediator%vec,ierr);CHKERRQ(ierr)

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

    if (this%dataset_solution) then
      !Import radionuclide inventory from external neutronics code calculations
      if (.not.cur_waste_form%breached) then
        do k = 1,num_species
          cur_waste_form%rad_mass_fraction(k) = 1.d-20
          cur_waste_form%rad_concentration(k) = &
            cur_waste_form%rad_mass_fraction(k) / &
            cwfm%rad_species_list(k)%formula_weight
        enddo
      else
        cur_criticality => this%criticality_mediator%criticality_list 
        dataset => cur_criticality%crit_mech%rad_dataset
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
      endif
      
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

    else !--------------------------------------------------------------------

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
      call ludcmp(Jacobian,num_species,indices,i)
      ! solve step 2/2: LU back substitution linear solve
      call lubksb(Jacobian,num_species,indices,rhs)
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
  use Material_Aux_class
  use Reactive_Transport_Aux_module, only : rt_min_saturation
  use Grid_module
  
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
  PetscInt :: i, j, k
  PetscInt :: num_species
  PetscInt :: local_id, ghosted_id
  PetscInt :: idof
  PetscReal :: inst_diss_molality
  PetscReal, pointer :: vec_p(:)  
  PetscReal, pointer :: xx_p(:)
  PetscInt :: fmdm_count_global, fmdm_count_local
  PetscLogDouble :: log_start_time, log_end_time
  character(len=MAXWORDLENGTH) :: word
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type), pointer :: grid
! -----------------------------------------------------------

  call PetscTime(log_start_time, ierr);CHKERRQ(ierr)
  
  fmdm_count_global = 0
  fmdm_count_local = 0
  global_auxvars => this%realization%patch%aux%Global%auxvars
  material_auxvars => this%realization%patch%aux%Material%auxvars
  grid => this%realization%patch%grid

  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(this%realization%field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  
  cur_waste_form => this%waste_form_list
  i = 0
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
    !
    if (associated(this%criticality_mediator)) then
      call CriticalitySolve(this%criticality_mediator,this%realization,time, &
                            cur_waste_form%scaling_factor,ierr)
    endif
    cur_waste_form => cur_waste_form%next
  enddo
 
  ! ideally, this print statement would go inside the dissolution subroutine
  call MPI_Allreduce(fmdm_count_local,fmdm_count_global,ONE_INTEGER_MPI, &
                     MPI_INTEGER,MPI_SUM,this%realization%option%mycomm,ierr)
  if ((fmdm_count_global > 0) .and. &
      this%realization%option%print_screen_flag) then
    write(word,'(i5)') fmdm_count_global
  ! ** START (this can be removed after FMDM profiling is finished) **
    write(*,'(a)') '== ' // adjustl(trim(word)) // ' call(s) to FMDM.'
  ! ** END (this can be removed after FMDM profiling is finished) **
  endif
  
  call VecRestoreArrayF90(this%realization%field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
  call PetscTime(log_end_time, ierr);CHKERRQ(ierr)

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
  PetscReal, parameter :: time_conversion = 1.d0/(24.d0*3600.d0)
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
    avg_temp_local = avg_temp_local + &  ! Celcius
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
  PetscReal :: avg_temp_local
! --------------------------------------------------------------
  
 ! FMDM model: 
 !=======================================================
  integer ( kind = 4) :: success
  logical ( kind = 4) :: initialRun
  PetscReal :: time
  PetscReal :: Usource
  PetscReal :: avg_temp_global
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
    avg_temp_local = avg_temp_local + &  ! Celcius
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
    avg_temp_local = avg_temp_local + &  ! Celcius
               global_auxvars(ghosted_id)%temp * waste_form%scaling_factor(i)
  enddo
  call CalcParallelSUM(option,waste_form%rank_list,avg_temp_local, &
                       avg_temp_global)

  call AMP_surrogate_step(this%burnup, time, avg_temp_global, &
                          this%concentration, this%decay_time, &
                          this%outer_weights, &
                          this%inner_weights, this%scaler_means, &
                          this%scaler_variances, this%dissolution_rate)

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
  ! radionuclide mass fraction (4:2:end-1), 
  ! cumulative mass released (5:2:end) .

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
  
  n_check_vars=3 !number of scalar wf checkpoint variables
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
  call MPI_Allreduce(local_stride,stride,ONE_INTEGER_MPI, &
                  MPI_INTEGER,MPI_MAX,this%option%mycomm,ierr)
  call MPI_Allreduce(n_wf_local,n_wf_global,ONE_INTEGER_MPI, &
                     MPI_INTEGER,MPI_SUM,this%option%mycomm,ierr)
  
  !Create MPI vector and sequential vector for mapping
  call VecCreateMPI(this%option%mycomm,n_wf_local*stride,&
                    n_wf_global*stride,& 
                    global_wf_vec,ierr);CHKERRQ(ierr)
                    
  call VecCreateSeq(PETSC_COMM_SELF, n_wf_local*stride,local_wf_vec, &
                    ierr); CHKERRQ(ierr)
  
  call VecSetBlockSize(global_wf_vec, stride, ierr);CHKERRQ(ierr)
  call VecSetBlockSize(local_wf_vec, stride, ierr);CHKERRQ(ierr)
                                    
                    
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
    
    call VecSetValues(local_wf_vec,stride,indices,check_vars, &
                     INSERT_VALUES,ierr);CHKERRQ(ierr)
    cur_waste_form => cur_waste_form%next

  enddo
  
  !Create map and add values from the sequential vector to the global 
  call ISCreateBlock(this%option%mycomm,stride,n_wf_local,int_array, &
                     PETSC_COPY_VALUES,is, ierr); CHKERRQ(ierr)
  
  call VecScatterCreate(local_wf_vec,PETSC_NULL_IS,global_wf_vec, &
                        is,scatter_ctx, ierr);CHKERRQ(ierr)
                        
  call VecScatterBegin(scatter_ctx, local_wf_vec, global_wf_vec, &  
                       INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx, local_wf_vec, global_wf_vec, & 
                     INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  
  dataset_name='canister_properties'
  
  !Write the checkpoint file
  call HDF5WriteDataSetFromVec(dataset_name, this%option, global_wf_vec,&
           pm_grp_id, H5T_NATIVE_DOUBLE)
           
  call VecScatterDestroy(scatter_ctx, ierr);CHKERRQ(ierr)
  call ISDestroy(is, ierr);CHKERRQ(ierr)
  call VecDestroy(global_wf_vec, ierr);CHKERRQ(ierr)
  call VecDestroy(local_wf_vec, ierr);CHKERRQ(ierr)

end subroutine PMWFCheckpointHDF5

! ***************************************************************************** !


subroutine PMWFRestartHDF5(this,pm_grp_id)
  ! 
  ! Restarts data associated with waste form process model
  ! from the "canister_properties" dataset for a given wf pm: 
  ! canister vitality (1), canister volume (2), breach time (3),
  ! radionuclide mass fraction (4:2:end-1), 
  ! cumulative mass released (5:2:end) .
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
  
  n_check_vars=3 !number of scalar wf checkpoint variables
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
  call MPI_Allreduce(local_stride,stride,ONE_INTEGER_MPI, &
                  MPI_INTEGER,MPI_MAX,this%option%mycomm,ierr)
  call MPI_Allreduce(n_wf_local,n_wf_global,ONE_INTEGER_MPI, &
                     MPI_INTEGER,MPI_SUM,this%option%mycomm,ierr)
                     
  !Create MPI vector into which HDF5 will read, and sequential vector
  !for wf information on a given process.
  call VecCreateMPI(this%option%mycomm,n_wf_local*stride,&
                    n_wf_global*stride,& 
                    global_wf_vec,ierr);CHKERRQ(ierr)
  
  call VecCreateSeq(PETSC_COMM_SELF, n_wf_local*stride,local_wf_vec, &
                    ierr); CHKERRQ(ierr)
  
  call VecSetBlockSize(global_wf_vec, stride, ierr);CHKERRQ(ierr)
  call VecSetBlockSize(local_wf_vec, stride, ierr);CHKERRQ(ierr)
  
  !Read the data
  dataset_name = 'canister_properties'
  call HDF5ReadDataSetInVec(dataset_name, this%option, global_wf_vec, &
                             pm_grp_id, H5T_NATIVE_DOUBLE)
  
  !Create map between MPI and sequential vectors
  call ISCreateBlock(this%option%mycomm,stride,n_wf_local,int_array, &
                     PETSC_COPY_VALUES,is, ierr); CHKERRQ(ierr)
                     
  call VecScatterCreate(global_wf_vec,is,local_wf_vec, &
                        PETSC_NULL_IS,scatter_ctx, ierr);CHKERRQ(ierr)
  
  !Get the data from the MPI vector
  call VecScatterBegin(scatter_ctx, global_wf_vec, local_wf_vec, &  
                       INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx, global_wf_vec, local_wf_vec, & 
                     INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  
  !Convert the data to a Fortran array
  call VecGetArrayF90(local_wf_vec, local_wf_array, ierr); CHKERRQ(ierr)

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
  
  call VecRestoreArrayF90(local_wf_vec, local_wf_array, ierr);CHKERRQ(ierr)
  call VecScatterDestroy(scatter_ctx, ierr);CHKERRQ(ierr)
  call ISDestroy(is, ierr);CHKERRQ(ierr)
  call VecDestroy(global_wf_vec, ierr);CHKERRQ(ierr)
  call VecDestroy(local_wf_vec, ierr);CHKERRQ(ierr)

end subroutine PMWFRestartHDF5

! ************************************************************************** !

subroutine PMWFCheckpointBinary(this, viewer)
  ! 
  ! Checkpoints data associated with the waste form process model
  ! into a checkpiont binary file for a given wf pm: 
  ! canister vitality (1), canister volume (2), breach time (3),
  ! radionuclide mass fraction (4:2:end-1), 
  ! cumulative mass released (5:2:end) .

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
  
  n_check_vars=3 !number of scalar wf checkpoint variables
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
  call MPI_Allreduce(local_stride,stride,ONE_INTEGER_MPI, &
                  MPI_INTEGER,MPI_MAX,this%option%mycomm,ierr)
  call MPI_Allreduce(n_wf_local,n_wf_global,ONE_INTEGER_MPI, &
                     MPI_INTEGER,MPI_SUM,this%option%mycomm,ierr)                  
                     
  !Create MPI vector and sequential vector for mapping
  call VecCreateMPI(this%option%mycomm,n_wf_local*stride,&
                    n_wf_global*stride,& 
                    global_wf_vec,ierr);CHKERRQ(ierr)
                    
  call VecCreateSeq(PETSC_COMM_SELF, n_wf_local*stride,local_wf_vec, &
                    ierr); CHKERRQ(ierr)
  
  call VecSetBlockSize(global_wf_vec, stride, ierr);CHKERRQ(ierr)
  call VecSetBlockSize(local_wf_vec, stride, ierr);CHKERRQ(ierr)
                                    
                    
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
    
    call VecSetValues(local_wf_vec,stride,indices,check_vars, &
                     INSERT_VALUES,ierr);CHKERRQ(ierr)
    cur_waste_form => cur_waste_form%next

  enddo
  
  !Create map and add values from the sequential vector to the global 
  call ISCreateBlock(this%option%mycomm,stride,n_wf_local,int_array, &
                     PETSC_COPY_VALUES,is, ierr); CHKERRQ(ierr)
  
  call VecScatterCreate(local_wf_vec,PETSC_NULL_IS,global_wf_vec, &
                        is,scatter_ctx, ierr);CHKERRQ(ierr)
                        
  call VecScatterBegin(scatter_ctx, local_wf_vec, global_wf_vec, &  
                       INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx, local_wf_vec, global_wf_vec, & 
                     INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  
  !Write the checkpoint file

  call VecView(global_wf_vec,viewer,ierr);CHKERRQ(ierr)         
           
  call VecScatterDestroy(scatter_ctx, ierr);CHKERRQ(ierr)
  call ISDestroy(is, ierr);CHKERRQ(ierr)
  call VecDestroy(global_wf_vec, ierr);CHKERRQ(ierr)
  call VecDestroy(local_wf_vec, ierr);CHKERRQ(ierr)
  
  

end subroutine PMWFCheckpointBinary

! ***************************************************************************** !

subroutine PMWFRestartBinary(this, viewer)

  ! 
  ! Restarts data associated with waste form process model
  ! from a checkpoint binary file for a given wf pm: 
  ! canister vitality (1), canister volume (2), breach time (3),
  ! radionuclide mass fraction (4:2:end-1), 
  ! cumulative mass released (5:2:end) .
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
  
  n_check_vars=3 !number of scalar wf checkpoint variables
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
  call MPI_Allreduce(local_stride,stride,ONE_INTEGER_MPI, &
                  MPI_INTEGER,MPI_MAX,this%option%mycomm,ierr)
  call MPI_Allreduce(n_wf_local,n_wf_global,ONE_INTEGER_MPI, &
                     MPI_INTEGER,MPI_SUM,this%option%mycomm,ierr)
                     
  !Create MPI vector into which HDF5 will read, and sequential vector
  !for wf information on a given process.
  call VecCreateMPI(this%option%mycomm,n_wf_local*stride,&
                    n_wf_global*stride,& 
                    global_wf_vec,ierr);CHKERRQ(ierr)
  
  call VecCreateSeq(PETSC_COMM_SELF, n_wf_local*stride,local_wf_vec, &
                    ierr); CHKERRQ(ierr)
  
  call VecSetBlockSize(global_wf_vec, stride, ierr);CHKERRQ(ierr)
  call VecSetBlockSize(local_wf_vec, stride, ierr);CHKERRQ(ierr)
  
  !Read the data
  call VecLoad(global_wf_vec,viewer, ierr);CHKERRQ(ierr)
  
  !Create map between MPI and sequential vectors
  call ISCreateBlock(this%option%mycomm,stride,n_wf_local,int_array, &
                     PETSC_COPY_VALUES,is, ierr); CHKERRQ(ierr)
                     
  call VecScatterCreate(global_wf_vec,is,local_wf_vec, &
                        PETSC_NULL_IS,scatter_ctx, ierr);CHKERRQ(ierr)
  
  !Get the data from the MPI vector
  call VecScatterBegin(scatter_ctx, global_wf_vec, local_wf_vec, &  
                       INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx, global_wf_vec, local_wf_vec, & 
                     INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  
  !Convert the data to a Fortran array
  call VecGetArrayF90(local_wf_vec, local_wf_array, ierr); CHKERRQ(ierr)

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
  
  call VecRestoreArrayF90(local_wf_vec, local_wf_array, ierr);CHKERRQ(ierr)
  call VecScatterDestroy(scatter_ctx, ierr);CHKERRQ(ierr)
  call ISDestroy(is, ierr);CHKERRQ(ierr)
  call VecDestroy(global_wf_vec, ierr);CHKERRQ(ierr)
  call VecDestroy(local_wf_vec, ierr);CHKERRQ(ierr)


end subroutine PMWFRestartBinary

! ***************************************************************************** !

subroutine PMWFInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
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
! id: [-] file id number
! --------------
  PetscInt :: id
! --------------

  id = INPUT_RECORD_UNIT
  
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  
end subroutine PMWFInputRecord

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
!   call CriticalityStrip(this%criticality_mediator)

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
        call DeallocateArray(prev_mechanism%outer_weights)
        call DeallocateArray(prev_mechanism%inner_weights)
        call DeallocateArray(prev_mechanism%scaler_means)
        call DeallocateArray(prev_mechanism%scaler_variances)
    end select
    deallocate(prev_mechanism)
    nullify(prev_mechanism)
  enddo
  nullify(this%mechanism_list)

end subroutine PMWFMechanismStrip

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

  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): waste form process model object
! ---------------------------------
  class(pm_waste_form_type) :: this
! ---------------------------------
  character(len=MAXWORDLENGTH) :: word

  write(word,'(f12.1)') this%cumulative_time
  
  write(*,'(/,a)') 'PM Waste Form time = ' // trim(adjustl(word)) // ' seconds'
  call PMBaseDestroy(this)
  call PMWFStrip(this)
  
end subroutine PMWFDestroy

! ************************************************************************** !

subroutine CriticalityInit(this)
  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none

  type(criticality_type), pointer :: this

  allocate(this)
  allocate(this%crit_event)
  allocate(this%crit_mech)
  nullify(this%region)
  nullify(this%next)

  this%crit_event%steady_state = PETSC_FALSE
  this%crit_event%crit_flag = PETSC_FALSE
  this%crit_event%crit_start = 0.d0
  this%crit_event%crit_end = 0.d0



end subroutine CriticalityInit

! ************************************************************************** !

subroutine CriticalityMechInit(this)
  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18

  use Dataset_Ascii_class

  implicit none

  type(criticality_mechanism_type), pointer :: this

  allocate(this)
  nullify(this%next)

  this%decay_heat = 0.d0
  this%crit_heat = 0.d0
  this%sw = 0.d0
  this%rho_w = 0.d0
  this%temperature = 0.d0
  this%k_effective = 0.d0

  this%rad_dataset => DatasetAsciiCreate()
  this%heat_dataset => DatasetAsciiCreate()

end subroutine CriticalityMechInit

! ************************************************************************** !

subroutine CriticalityMediatorInit(this)
  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none

  type(criticality_mediator_type), pointer :: this

  nullify(this%data_mediator)
  nullify(this%criticality_list)
  nullify(this%crit_mech_list)


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
  
function CriticalityCreate()

  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none

  type(criticality_type), pointer :: CriticalityCreate
  type(criticality_type), pointer :: crit

  allocate(crit)
  call CriticalityInit(crit)

  CriticalityCreate => crit

end function CriticalityCreate

! ************************************************************************** !

function CriticalityMechCreate()

  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none

  type(criticality_mechanism_type), pointer :: CriticalityMechCreate
  type(criticality_mechanism_type), pointer :: crit

  allocate(crit)
  call CriticalityMechInit(crit)

  CriticalityMechCreate => crit

end function CriticalityMechCreate

! ************************************************************************** !

subroutine ReadCriticalityMech(this,input,option,keyword,error_string,found)

  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18

  use Input_Aux_module
  use Option_module
  use String_module
  use Dataset_Ascii_class

  implicit none

  type(criticality_mediator_type), pointer :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword, internal_units
  character(len=MAXSTRINGLENGTH) :: error_string,temp_string

  PetscBool :: found

  character(len=MAXWORDLENGTH) :: word
  type(criticality_mechanism_type), pointer :: new_crit_mech, cur_crit_mech

  error_string = trim(error_string) // ',CRITICALITY'
  found = PETSC_TRUE
  select case(trim(keyword))
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
          case('NAME')
            call InputReadCard(input,option,word)
            call InputErrorMsg(input,option, &
                  'criticality mechanism assignment',error_string)
            call StringToUpper(word)
            new_crit_mech%mech_name = trim(word)
          case('HEAT_OF_CRITICALITY')
            call InputReadDouble(input,option,new_crit_mech%crit_heat)
            call InputErrorMsg(input,option,'HEAT_OF_CRITICALITY',error_string)
          case('DECAY_HEAT')
            call InputReadCard(input,option,word)
            select case (trim(word))
              case('TOTAL')
                new_crit_mech%heat_source_cond = 1
              case('ADDITIONAL')
                new_crit_mech%heat_source_cond = 2
              case('CYCLIC')
                new_crit_mech%heat_source_cond = 3
            end select
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit
              call InputReadCard(input,option,word,PETSC_FALSE)
              select case(trim(word))
                case('DATASET')
                  internal_units = 'MW'
                  call InputReadFilename(input,option,new_crit_mech% &
                          heat_dataset_name)
                  call DatasetAsciiReadFile(new_crit_mech%heat_dataset, &
                          new_crit_mech%heat_dataset_name,temp_string, &
                          internal_units,error_string,option)
                  new_crit_mech%heat_dataset%time_storage% &
                          time_interpolation_method = 2
              end select
            enddo
            call InputPopBlock(input,option)
          case('INVENTORY')
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit
              call InputReadCard(input,option,word,PETSC_FALSE)
              select case(trim(word))
                case('DATASET')
                  internal_units = 'g/g'
                  call InputReadFilename(input,option,new_crit_mech% &
                          rad_dataset_name)
                  call DatasetAsciiReadFile(new_crit_mech%rad_dataset, &
                          new_crit_mech%rad_dataset_name,temp_string, &
                          internal_units,error_string,option)
                  new_crit_mech%rad_dataset%time_storage% &
                          time_interpolation_method = 2
              end select
            enddo
            call InputPopBlock(input,option)
        end select
      enddo
      call InputPopBlock(input,option)
      if (.not. associated(this%crit_mech_list)) then
        this%crit_mech_list => new_crit_mech
      else
        cur_crit_mech => this%crit_mech_list
        do
          if (.not. associated(cur_crit_mech)) exit
          if (.not. associated(cur_crit_mech%next)) then
            cur_crit_mech => new_crit_mech
          endif
          cur_crit_mech => cur_crit_mech%next
        enddo
      endif
      nullify(cur_crit_mech)
    case default
      found = PETSC_FALSE
  end select

end subroutine ReadCriticalityMech

! ************************************************************************** !

subroutine CriticalityCalc(this,time,ierr)

  ! Calculate mass and heat source terms as a function of time.
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none

  type(criticality_mechanism_type) :: this
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

subroutine CriticalityInitializeRun(this, realization, option)

  ! Author: Michael Nole
  ! Date: 11/01/18

  use Realization_Subsurface_class
  use Realization_Base_class
  use Option_module

  implicit none

  type(criticality_mediator_type), pointer :: this
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option

  type(criticality_type), pointer :: cur_criticality
  PetscInt :: vec_size, i, j
  PetscInt, allocatable :: energy_indices_in_residual(:)
  PetscErrorCode :: ierr
  IS :: is

  call RealizCreateFlowMassTransferVec(realization)
  this%data_mediator => DataMediatorVecCreate()
  call this%data_mediator%AddToList(realization%flow_data_mediator_list)

  cur_criticality => this%criticality_list
  vec_size = 0

  do
    if (.not. associated(cur_criticality)) exit
    vec_size = vec_size + cur_criticality%region%num_cells
    cur_criticality => cur_criticality%next
  enddo

  call VecCreateSeq(PETSC_COMM_SELF, vec_size,this%data_mediator%vec, &
                    ierr);CHKERRQ(ierr)
  call VecSetFromOptions(this%data_mediator%vec,ierr); CHKERRQ(ierr)

  cur_criticality => this%criticality_list
  allocate(energy_indices_in_residual(vec_size))
  j = 0
  do
    if (.not. associated(cur_criticality)) exit
      do i = 1, cur_criticality%region%num_cells
        j = j + 1
        energy_indices_in_residual(j) = cur_criticality%region% &
                                        cell_ids(i) * option%nflowdof - 1
      enddo
    cur_criticality => cur_criticality%next
  enddo
  energy_indices_in_residual(:) = energy_indices_in_residual(:) + &
      realization%patch%grid%global_offset*option%nflowdof

  this%total_num_cells = j

  call ISCreateGeneral(option%mycomm,vec_size, &
                       energy_indices_in_residual, &
                       PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)
  call VecScatterCreate(this%data_mediator%vec,PETSC_NULL_IS, &
                        realization%field%flow_r, is, &
                        this%data_mediator%scatter_ctx, ierr); CHKERRQ(ierr)
  if (allocated(energy_indices_in_residual)) then
      deallocate(energy_indices_in_residual)
  endif

  call ISDestroy(is,ierr);CHKERRQ(ierr)

end subroutine CriticalityInitializeRun

! ************************************************************************** !

subroutine AssignCritMech(this)

  use String_module

  implicit none

  type(criticality_mediator_type), pointer :: this

  type(criticality_mechanism_type), pointer :: cur_mechanism
  type(criticality_type), pointer :: cur_criticality

  cur_criticality => this%criticality_list
  do
    if(.not. associated(cur_criticality)) exit
    cur_mechanism => this%crit_mech_list
    do
      if (.not. associated(cur_mechanism)) exit
      if (StringCompare(cur_criticality%crit_event%mech_name, &
          cur_mechanism%mech_name)) then
        cur_criticality%crit_mech => cur_mechanism
        exit
      endif
    enddo
    cur_criticality => cur_criticality%next
  enddo

end subroutine AssignCritMech

! ************************************************************************** !

subroutine CriticalitySolve(this,realization,time,scaling_factor,ierr)
  !
  !Author: Michael Nole
  !Date: 11/05/18
  !
  use Realization_Subsurface_class

  implicit none

  type(criticality_mediator_type), pointer :: this
  class(realization_subsurface_type), pointer :: realization
  PetscReal :: time
  PetscReal, pointer :: scaling_factor(:)
  PetscErrorCode :: ierr

  PetscInt :: i,j
  type(criticality_type), pointer :: cur_criticality
  PetscReal, pointer :: heat_source(:)

  call VecGetArrayF90(this%data_mediator%vec,heat_source, &
                      ierr);CHKERRQ(ierr)

  cur_criticality => this%criticality_list
  j = 0
  do
    if (.not. associated(cur_criticality)) exit

    if (time >= cur_criticality%crit_event%crit_start .and. time < &
            cur_criticality%crit_event%crit_end) then
      cur_criticality%crit_event%crit_flag = PETSC_TRUE
    else
      cur_criticality%crit_event%crit_flag = PETSC_FALSE
    endif

    call CriticalityCalc(cur_criticality%crit_mech,time,ierr)

    do i = 1, cur_criticality%region%num_cells
      j = j + 1
      heat_source(j) = cur_criticality%crit_mech%decay_heat

      if (cur_criticality%crit_event%crit_flag) then
        heat_source(j) = heat_source(j) + cur_criticality%crit_mech%crit_heat
      endif

      ! Distribute heat source throughout all cells in a waste package
      heat_source(j) = heat_source(j) * scaling_factor(j)
    enddo

    cur_criticality => cur_criticality%next
  enddo

  call VecRestoreArrayF90(this%data_mediator%vec,heat_source, &
                          ierr);CHKERRQ(ierr)

end subroutine CriticalitySolve

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
  PetscInt :: length, i, icount
  PetscInt :: icol
  PetscInt :: ndims
  PetscInt, pointer :: dims(:)
  PetscReal, pointer :: real_buffer(:)
  PetscErrorCode :: ierr

  integer(HID_T) :: file_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: hdf5_err
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

subroutine AMP_surrogate_step (burnup, sTme, current_temp_C, &
                               conc, decay_time, &
                               outer_weights, inner_weights, &
                               scaler_means, scaler_variances, &
                               fuelDisRate)
  implicit none
  PetscReal, intent(in) :: sTme
  PetscReal, intent(in) :: burnup
  PetscReal, intent(in) :: current_temp_C 
  PetscReal, intent(in) :: decay_time
  ! four environmental concentrations
  PetscReal, intent(in) :: conc(:)
  ! ANN weights
  PetscReal, intent(in) :: outer_weights(:)
  PetscReal, intent(in) :: inner_weights(:,:)
  ! standardization scaler parameters
  PetscReal, intent(in) :: scaler_means(:)
  PetscReal, intent(in) :: scaler_variances(:)
  ! output
  PetscReal, intent(out) :: fuelDisRate ! g/m2/yr
  
  ! local variables
  PetscReal :: yTme
  PetscReal :: f1, f2, f3, f4, f5
  PetscReal :: AOF, rad0a, rad0
  PetscReal :: node_sum
  PetscInt :: i
  PetscInt :: N ! number of nodes in ANN
  ! features
  PetscReal, dimension(6) :: f
  ! hidden layer nodes values
  PetscReal, dimension(100) :: h
  ! constants
  PetscReal, parameter :: UO2_molar_mass = 270.0d0 ! g/mol

  yTme = sTme/60.0d0/60.0d0/24.0d0/DAYS_PER_YEAR

  ! calculate dose rate at the fuel surface (rad0)
  AOF = yTme + decay_time

  f2 = log(AOF)
  f1 = f2**2.0d0
  f3 = 1.0d0/f2
  f4 = f2/AOF
  f5 = exp(burnup/25.26892627636246d0)

  rad0a = -206.0634818750711d0   - 0.7631591788870090d0*f1 &
        + 20.97112373957833d0*f2 + 678.8463343193430d0*f3 &
        - 506.7149017370657d0*f4 + 0.1555448893425319d0*f5
  rad0 = max(exp(rad0a),5.0d-3)

  ! features
  f(1) = current_temp_C + 273.15d0
  f(2) = log10(conc(1)) ! Env_CO3_2n
  f(3) = log10(conc(2)) ! Env_O2
  f(4) = log10(conc(3)) ! Env_Fe_2p
  f(5) = log10(conc(4)) ! Env_H2
  f(6) = log10(rad0)

  ! standardize
  do i = 1,6
    f(i) = (f(i) - scaler_means(i))/sqrt(scaler_variances(i))
  enddo

  N = 100 ! 100 nodes

  do i = 1,N
    node_sum = dot_product(f, inner_weights(1:6,i)) + inner_weights(7,i)
    if (node_sum < 0.0d0) then
      node_sum = 0.0d0
    endif
    h(i) = node_sum
  enddo

  fuelDisRate= 10**(dot_product(h,outer_weights(1:N)) + &
                          outer_weights(101))*UO2_molar_mass

end subroutine AMP_surrogate_step

! ************************************************************************** !

end module PM_Waste_Form_class

