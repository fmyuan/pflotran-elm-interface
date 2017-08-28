module PM_WIPP_SrcSink_class

! MODULE DESCRIPTION:
! ============================================================================
! This process model tracks chemical species in waste panel inventories 
! involved in the generation/uptake of H2 gas and brine (water).
! The chemical species are tracked for accounting purposes and also for the
! calculation of gas and brine generation rates as determined by reaction
! rate constants of chemical reactions between the chemical species, and the 
! availability of limiting chemical species.
! The calculated gas and brine generation rates are assigned to fluid and gas
! source terms, respectively, in the flow process model.
!=============================================================================

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PM_Base_class
  use Region_module
  use PFLOTRAN_Constants_module
  use Realization_Subsurface_class
  use Data_Mediator_Vec_class

  implicit none

  private
  
  ! molar mass parameters [kg/mol]
  PetscReal, parameter :: MW_FE = 0.055847d0
  PetscReal, parameter :: MW_CELL = 0.027023d0
  PetscReal, parameter :: MW_FEOH2 = 0.08986d0
  PetscReal, parameter :: MW_CO2 = 0.0440098d0
  PetscReal, parameter :: MW_N2 = 0.02801348d0
  PetscReal, parameter :: MW_H2S = 0.03408188d0
  PetscReal, parameter :: MW_FES = 0.087911d0
  PetscReal, parameter :: MW_MGO = 0.040304d0
  PetscReal, parameter :: MW_MGOH2 = 0.05832d0
  PetscReal, parameter :: MW_MGCO3 = 0.084314d0
  PetscReal, parameter :: MW_HYDRO = 0.467636d0
  PetscReal, parameter :: MW_NO3 = (14.0067d-3 + 3.d0*15.9994d-3) ! not in database
  PetscReal, parameter :: MW_SO4 = (32.065d-3 + 4.d0*15.9994d-3) ! not in database
  PetscReal, parameter :: MW_H2 = 2.015880d-3
  PetscReal, parameter :: MW_H2O = 1.801528d-2
  
  ! density parameters [kg/m3]
  PetscReal, parameter :: DN_FE = 7870.d0
  PetscReal, parameter :: DN_FEOH2 = 3400.d0
  PetscReal, parameter :: DN_CELL = 1100.d0
  PetscReal, parameter :: DN_FES = 4700.d0
  PetscReal, parameter :: DN_HYDRO = 2300.d0
  PetscReal, parameter :: DN_MGCO3 = 3050.d0
  PetscReal, parameter :: DN_MGO = 3600.d0
  PetscReal, parameter :: DN_MGOH2 = 2370.d0
  PetscReal, parameter :: DN_SALT = 2180.d0

! OBJECT chem_species_type:
! =========================
! ---------------------------------------------------------------------------
! Description:  This object describes a chemical species in the waste panel
! that is tracked for accounting purposes, or influences a reaction rate
! constant for gas/brine generation. It is a member object of the
! inventory_type object.
! ---------------------------------------------------------------------------
! initial_conc_mol(:): [mol-species/m3-bulk] initial molar concentration of
!    a chemical species in the waste panel, and indexed by each cell in the
!    waste panel region
! inst_rate(:): [mol/m3-bulk/sec] instantaneous reaction rate constant of a
!    chemical species in the waste panel, and indexed by each cell in the
!    waste panel region
! current_conc_mol(:): [mol-species/m3-bulk] current molar concentration of
!    a chemical species in the waste panel, and indexed by each cell in the
!    waste panel region
! current_conc_kg(:): [kg-species/m3-bulk] current concentration of a
!    chemical species in the waste panel (in BRAGFLO units), and indexed by 
!    each cell in the waste panel region
! molar_mass: [kg/mol] molar mass of the chemical species
! tot_mass_in_panel: [kg] total mass of the chemical species in the entire
!    panel volume
! -------------------------------------------
  type, public :: chem_species_type
    PetscReal, pointer :: initial_conc_mol(:)
    PetscReal, pointer :: inst_rate(:)  
    PetscReal, pointer :: current_conc_mol(:)
    PetscReal, pointer :: current_conc_kg(:)
    PetscReal :: molar_mass              
    PetscReal :: tot_mass_in_panel           
  end type chem_species_type
! -------------------------------------------
  
! OBJECT inventory_type:
! ======================
! ---------------------------------------------------------------------------
! Description:  This object describes a waste panel's current inventory, and  
! is made up of several chem_species_type objects. It also stores the initial 
! values of several tracked species as well as the steel drum content. A 
! pointer to a pre-inventory object aids in this object's initialization.
! ---------------------------------------------------------------------------
! name: name string
! Fe_s: solid iron species object
! FeOH2_s: solid corrosion product iron hydroxide species object
! BioDegs_s: solid biodegradables (cellulosics+rubber+1.7*plastics) species object
! FeS_s: solid iron sulfide species object
! MgO_s: solid magnesium oxide species object
! MgOH2_s: solid magnesium hydroxide (brucite) species object
! Mg5CO34OH24H2_s: solid hydromagnesite species object
! MgCO3_s: solid magnesium carbonate species object
! *_in_panel: [kg] total initial mass of species/material in waste panel
! drum_conc: [-/m3] number of steel drums per unit volume of waste panel
!    note: this parameter is defined by the parameter of the same name in 
!    the preinventory: pointer to the pre-inventory object, which stores the
!    initial inventory values
! ---------------------------------------------------
  type, public :: inventory_type
    character(len=MAXWORDLENGTH) :: name
    type(chem_species_type) :: Fe_s 
    type(chem_species_type) :: FeOH2_s
    type(chem_species_type) :: BioDegs_s
    type(chem_species_type) :: FeS_s
    type(chem_species_type) :: MgO_s
    type(chem_species_type) :: MgOH2_s
    type(chem_species_type) :: Mg5CO34OH24H2_s
    type(chem_species_type) :: MgCO3_s
  ! Initial Values:
    PetscReal :: Fe_in_panel         
    PetscReal :: MgO_in_panel        
    PetscReal :: Cellulose_in_panel  
    PetscReal :: RubberPlas_in_panel 
    PetscReal :: Biodegs_in_panel    
    PetscReal :: Nitrate_in_panel    
    PetscReal :: Sulfate_in_panel    
    PetscReal :: drum_conc    
    type(pre_inventory_type), pointer :: preinventory
  end type inventory_type
! ---------------------------------------------------
  
! OBJECT pre_inventory_type:
! ==========================
! ---------------------------------------------------------------------------
! Description:  This object describes the initial waste inventory emplaced in
! a waste panel. Several of the member parameters are ALGEBRA input 
! parameters, preprocessed by ALGEBRA for BRAGFLO input. A pre-inventory object 
! can be defined for each waste panel. Alternatively, a single pre-inventory 
! object can be defined for an entire repository that is made up of several 
! waste panels. In this latter case, "vrepos" should be specifed and each 
! waste panel will get assigned a portion of the single pre-inventory values 
! which are scaled by the volume ratio vol-panel/vrepos.
! Note: RH-remotely handled; CH-contact handled
! ---------------------------------------------------------------------------
! name: name string
! ironchw: [kg] mass of Fe-based material in CH waste
! ironrhw: [kg] mass of Fe-based material in RH waste
! irncchw: [kg] mass of Fe containers for CH waste
! irncrhw: [kg] mass of Fe containers for RH waste
! cellchw: [kg] mass of cellulosics in CH waste
! cellrhw: [kg] mass of cellulosics in RH waste
! celcchw: [kg] mass of cellulosics in container materials for CH waste
! celcrhw: [kg] mass of cellulosics in container materials for RH waste
! celechw: [kg] mass of cellulosics in emplacement materials for CH waste
! celerhw: [kg] mass of cellulosics in emplacement materials for RH waste
! rubbchw: [kg] mass of rubber in CH waste
! rubbrhw: [kg] mass of rubber in RH waste
! rubcchw: [kg] mass of rubber in container materials for CH waste
! rubcrhw: [kg] mass of rubber in container materials for RH waste
! rubechw: [kg] mass of rubber in emplacement materials for CH waste
! ruberhw: [kg] mass of rubber in emplacement materials for RH waste  
! plaschw: [kg] mass of plastics in CH waste
! plasrhw: [kg] mass of plastics in RH waste
! plscchw: [kg] mass of plastics in container materials for CH waste
! plscrhw: [kg] mass of plastics in container materials for RH waste
! plsechw: [kg] mass of plastics in emplacement materials for CH waste
! plserhw: [kg] mass of plastics in emplacement materials for RH waste
! plasfac: [-] mass ratio of plastics to equivalent carbon
! mgo_ef: [-] MgO excess factor; ratio mol-MgO/mol-organic-carbon
! vrepos: [m3] volume of total repository (optional) note: this parameter
!    should be defined if the waste panel inventories are being distributed 
!    from a single pre-inventory, where the amount of pre-inventory going to 
!    each waste panel inventory is scaled by the volume ratio vol-panel/vrepos
! drum_conc: [-/m3] number of steel drums per unit volume of inventory
!    volume (which can be a waste panel volume or an entire repository volume)
!    note: this defines the parameter of the same name in the inventory object
! nitrate: [mol] initial amount of aqueous nitrate in inventory
! sulfate: [mol] initial amount of aqueous sulfate in inventory
! next: pointer to the next pre-inventory object in a linked list
! -------------------------------------------
  type, public :: pre_inventory_type
    character(len=MAXWORDLENGTH) :: name
  ! ALGEBRA parameters:
    PetscReal :: ironchw 
    PetscReal :: ironrhw
    PetscReal :: irncchw 
    PetscReal :: irncrhw
    PetscReal :: cellchw 
    PetscReal :: cellrhw
    PetscReal :: celcchw 
    PetscReal :: celcrhw  
    PetscReal :: celechw
    PetscReal :: celerhw  
    PetscReal :: rubbchw
    PetscReal :: rubbrhw 
    PetscReal :: rubcchw 
    PetscReal :: rubcrhw 
    PetscReal :: rubechw  
    PetscReal :: ruberhw    
    PetscReal :: plaschw  
    PetscReal :: plasrhw  
    PetscReal :: plscchw  
    PetscReal :: plscrhw  
    PetscReal :: plsechw  
    PetscReal :: plserhw   
    PetscReal :: plasfac   
    PetscReal :: mgo_ef    
    PetscReal :: vrepos    
    PetscReal :: drum_conc 
    PetscReal :: nitrate   
    PetscReal :: sulfate   
    type(pre_inventory_type), pointer :: next
  end type pre_inventory_type
! -------------------------------------------
  
! OBJECT srcsink_panel_type:
! ==========================
! ---------------------------------------------------------------------------
! Description:  This object represents a waste panel within a repository. It
! contains member pointers to it's inventory and region objects. It also
! stores the calculated gas/brine generation rates, indexed by each grid cell
! within the waste panel region. The pre-processed ALGEBRA parameters are
! also stored per waste panel because some depend on waste panel volume.
! MPI information is stored for efficient parallel processing algorithms used
! at the waste panel level.
! ---------------------------------------------------------------------------
! name: name string of the waste panel
! region: pointer to the waste panel's region object
! region_name: name string of the waste panel's region object
! inventory: pointer to the waste panel's inventory object
! inventory_name: name string of the waste panel's inventory object
! scaling_factor(:): [-] array of volume scaling factors for each grid cell
!    in the waste panel region
! gas_generation_rate(:): [mol-H2/m3-bulk/sec] array of the current gas 
!    generation rate for each grid cell in the waste panel region
! brine_generation_rate(:): [mol-H2O/m3-bulk/sec] array of the current brine 
!    generation rate for each grid cell in the waste panel region
! rxnrate_Fe_corrosion(:): [mol-Fe/m3-bulk/sec] array of the current anoxic iron 
!    corrosion rate constant
! rxnrate_cell_biodeg(:): [mol-cell/m3-bulk/sec] array of the current
!     biodegradation rate constant
! rxnrate_Fe_sulf(:): [mol-Fe/m3-bulk/sec] array of the current (metallic) iron 
!    sulfidation to FeS rate constant
! rxnrate_FeOH2_sulf(:): [mol-FeOH2/m3-bulk/sec] array of the current iron 
!    corrosion product iron hydroxide (Fe(OH)2) sulfidation to FeS rate constant
! rxnrate_MgO_hyd(:): [mol-MgO/m3-bulk/sec] array of the current MgO hydration to 
!    brucite (Mg(OH)2) rate constant
! rxnrate_MgOH2_carb(:): [mol-hydromagnesite/m3-bulk/sec] array of the current 
!    brucite (Mg(OH)2) carbonation to hydromagnesite (Mg5(CO3)4(OH)2:4H2O) rate constant
! rxnrate_MgO_carb(:): [mol-MgO/m3-bulk/sec] array of the current 
!    MgO carbonation to magnesite (MgCO3) rate constant
! rxnrate_hydromag_conv(:): [mol-hydromagnesite/m3-bulk/sec] array of the current 
!    hydromagnesite (Mg5(CO3)4(OH)2:4H2O) conversion to 
!    brucite (Mg(OH)2) and magnesite (MgCO3) rate constant
! inundated_corrosion_rate: [mol-Fe/m3-bulk/sec] corrosion rate of iron
!    when inundated in brine
! humid_corrosion_rate: [mol-Fe/m3-bulk/sec], [-] corrosion rate of iron
!    in a humid environment
! inundated_biodeg_rate: [mol-cellulosics/m3-bulk/sec] biodegradation rate
!    when inundated in brine
! humid_biodeg_rate: [mol-cellulosics/m3-bulk/sec] biodegradation rate
!    in a humid environment
! inundated_brucite_rate: [mol-MgOH2/m3-bulk/sec] rate of MgO hydration
!    when inundated in brine
! humid_brucite_rate: [mol-MgOH2/m3-bulk/sec] rate of MgO hydration
!    in a humid environment
! solids_production(:): [-] normalized volume of solids produced in the panel
! F_NO3: [-] fraction of carbon consumed through denitrification process; 
!    eq. PA.85, section PA-4.2.5
! F_SO4: [-] fraction of carbon consumed through sulfate reduction process; 
!    eq. PA.85, section PA-4.2.5
! RXH2S_factor: [-] mol H2S / mol Carbon (cellulose) consumed by biodegradation
! RXCO2_factor: [-] mol CO2 / mol Carbon (cellulose) consumed by biodegradation
! RXH2_factor:  [-] mol H2  / mol Carbon (cellulose) consumed by biodegradation
! RXH2O_factor: [-] mol H2O / mol Carbon (cellulose) consumed by biodegradation
! volume: [m3] waste panel volume
! scale_by_volume: Boolean flag to scale given inventory to waste panel volume
! prev_dt: [sec] value of the previous time step
! id: [-] waste panel id number
! myMPIcomm: [-] MPI communicator number object
! myMPIgroup: [-] MPI group number object
! rank_list(:): [-] array of 1's and 0's used to determine local waste panels
! next: pointer to the next waste panel object in linked list
! ------------------------------------------------
  type, public :: srcsink_panel_type
    character(len=MAXWORDLENGTH) :: name
    type(region_type), pointer :: region
    character(len=MAXWORDLENGTH) :: region_name
    type(inventory_type) :: inventory
    character(len=MAXWORDLENGTH) :: inventory_name
    PetscReal, pointer :: scaling_factor(:)        

    PetscReal, pointer :: solids_production(:)  
    PetscReal, pointer :: gas_generation_rate(:) 
    PetscReal, pointer :: brine_generation_rate(:)
    
    PetscReal, pointer :: rxnrate_Fe_corrosion_inund(:)
    PetscReal, pointer :: rxnrate_Fe_corrosion_humid(:)
    PetscReal, pointer :: rxnrate_cell_biodeg_inund(:)
    PetscReal, pointer :: rxnrate_cell_biodeg_humid(:)
    PetscReal, pointer :: rxnrate_MgO_hyd_inund(:)
    PetscReal, pointer :: rxnrate_MgO_hyd_humid(:)
    
    PetscReal, pointer :: rxnrate_Fe_corrosion(:)
    PetscReal, pointer :: rxnrate_cell_biodeg(:)
    PetscReal, pointer :: rxnrate_FeOH2_sulf(:)
    PetscReal, pointer :: rxnrate_Fe_sulf(:)
    PetscReal, pointer :: rxnrate_MgO_hyd(:)
    PetscReal, pointer :: rxnrate_MgOH2_carb(:)
    PetscReal, pointer :: rxnrate_MgO_carb(:)
    PetscReal, pointer :: rxnrate_hydromag_conv(:)
    
    PetscReal :: inundated_corrosion_rate    
    PetscReal :: humid_corrosion_rate        
    PetscReal :: inundated_biodeg_rate       
    PetscReal :: humid_biodeg_rate           
    PetscReal :: inundated_brucite_rate      
    PetscReal :: humid_brucite_rate    
    
    PetscReal :: F_NO3
    PetscReal :: F_SO4
    PetscReal :: RXH2S_factor       
    PetscReal :: RXCO2_factor       
    PetscReal :: RXH2_factor       
    PetscReal :: RXH2O_factor       
    PetscReal :: volume           
    PetscBool :: scale_by_volume  
    PetscReal :: prev_dt
    PetscInt :: id
    PetscMPIInt :: myMPIcomm
    PetscMPIInt :: myMPIgroup
    PetscInt, pointer :: rank_list(:)
    type(srcsink_panel_type), pointer :: next
  end type srcsink_panel_type
! ------------------------------------------------

! OBJECT pm_wipp_srcsink_type:
! ============================
! ---------------------------------------------------------------------------
! Description:  This is the wipp-srcsink process model object. It contains
! several ALGEBRA parameters. It has a list of waste panels, and a list of 
! pre-inventory member objects. Several procedures, allow interfacing with
! the process model structure and extend the pm_base_type procedures.
! This is the highest level object in this module.
! ---------------------------------------------------------------------------
! alpharxn: [-] smoothing parameter used in s_eff calculation
! smin: [-] minimum brine saturation where a grid cell is considered dry,
!    note: this is not brine residual saturation, but much smaller
! satwick: [-] wicking saturation
! corrmco2: [m/s] iron corrosion rate in inundated conditions
! humcorr: [m/s] iron corrosion rate in humid conditions
! gratmici: [mol-cellulosics/kg-cellulosics/sec] biodegradation rate in inundated 
!    conditions 
! gratmich: [mol-cellulosics/kg-cellulosics/sec] biodegradation rate in humid 
!    conditions 
! brucitei: [mol-MgOH2/kg-MgO/sec] MgO hydration rate in inundated conditions
! bruciteh: [mol-MgOH2/kg-MgO/sec] MgO hydration rate in humid conditions
! hymagcon_rate: [mol-hydromagnesite/kg-hydromagnesite/sec] 
!                rate of hydromagnesite conversion
! drum_surface_area: [m2/drum] surface area of steel drums
! biogenfc: [-] parameter uniformly sampled between 0 and 1, used to account
!    for the uncertainty in whether microbial gas generation could be realized 
!    in the WIPP at experimentally measured rates
! probdeg: [-] flag value of 0, 1, or 2; indicates whether gas
!    generation is produced by biodegradation, and/or rubbers and plastics,
!    in addition to iron corrosion
! bioidx: [-] flag value of 0 or 1; indicates whether gas generation is
!    produced by biodegradation
! plasidx: [-] flag value of 0 or 1; indicates whether gas generation is
!    produced by rubbers and plastics
! stoic_mat(8,10): [mol/mol] stoichiometry matrix for the 8 reactions, and for 
!    10 reactant species; the numbering for the rows is the same as the first
!    number of the STCO_## parameters in the PFD analysis of the database,
!    but the column number is one off due to >0 indexing, e.g., STCO_34
!    is contained in stoic_mat(3,5)
! output_start_time: [sec] the time when output *.pnl files are generated,
!    with the default value set at 0.d0 sec
! waste_panel_list: linked list of waste panel objects that make up a
!    repository setting
! pre_inventory_list: linked list of pre-inventory objects
! data_mediator: pointer to data mediator object which stores the gas and
!    brine generation source terms 
! realization: pointer to the realization object
! --------------------------------------------------------------------
  type, public, extends(pm_base_type) :: pm_wipp_srcsink_type
    PetscReal :: alpharxn        
    PetscReal :: smin           
    PetscReal :: satwick        
    PetscReal :: corrmco2         
    PetscReal :: humcorr          
    PetscReal :: gratmici         
    PetscReal :: gratmich         
    PetscReal :: brucitei         
    PetscReal :: bruciteh         
    PetscReal :: hymagcon_rate      
    PetscReal :: drum_surface_area  
    PetscReal :: biogenfc           
    PetscInt :: probdeg             
    PetscInt :: bioidx              
    PetscInt :: plasidx             
    PetscReal :: output_start_time
    PetscReal :: stoic_mat(8,10)
    type(srcsink_panel_type), pointer :: waste_panel_list
    type(pre_inventory_type), pointer :: pre_inventory_list
    class(data_mediator_vec_type), pointer :: data_mediator
    class(realization_subsurface_type), pointer :: realization
  contains
    procedure, public :: PMWSSSetRealization
    procedure, public :: Setup => PMWSSSetup
    procedure, public :: Read => PMWSSRead
    procedure, public :: InitializeRun => PMWSSInitializeRun
    procedure, public :: InitializeTimestep => PMWSSInitializeTimestep
    procedure, public :: FinalizeTimestep => PMWSSFinalizeTimestep
    procedure, public :: Solve => PMWSSSolve
    procedure, public :: Output => PMWSSOutput
    procedure, public :: InputRecord => PMWSSInputRecord
    procedure, public :: Destroy => PMWSSDestroy
  end type pm_wipp_srcsink_type
! --------------------------------------------------------------------
  
  interface PMWSSTaperRxnrate
    module procedure PMWSSTaperRxnrate1
    module procedure PMWSSTaperRxnrate2
  end interface

  public :: PMWSSCreate, &
            PMWSSWastePanelCreate, &
            PMWSSPreInventoryCreate

contains

! *************************************************************************** !

function PMWSSCreate()
  !
  ! Creates and initializes the WIPP source/sink process model.
  !
  ! Author: Jenn Frederick
  ! Date: 1/31/2017
  !
  
  implicit none
  
! LOCAL VARIABLES:
! ================
! PMWSSCreate (output): pointer to new wipp-srcsink process model object
! pm: pointer to new wipp-srcsink process model object with shorter name
! ---------------------------------------------------
  class(pm_wipp_srcsink_type), pointer :: PMWSSCreate
  class(pm_wipp_srcsink_type), pointer :: pm
! ---------------------------------------------------
  
  allocate(pm)
  nullify(pm%waste_panel_list)
  nullify(pm%pre_inventory_list)
  nullify(pm%data_mediator)
  nullify(pm%realization)
  pm%name = 'wipp source sink'
  pm%alpharxn = UNINITIALIZED_DOUBLE
  pm%smin = UNINITIALIZED_DOUBLE
  pm%satwick = UNINITIALIZED_DOUBLE
  pm%corrmco2 = UNINITIALIZED_DOUBLE
  pm%humcorr = UNINITIALIZED_DOUBLE
  pm%gratmici = UNINITIALIZED_DOUBLE
  pm%gratmich = UNINITIALIZED_DOUBLE
  pm%brucitei = UNINITIALIZED_DOUBLE
  pm%bruciteh = UNINITIALIZED_DOUBLE
  pm%hymagcon_rate = UNINITIALIZED_DOUBLE
  pm%drum_surface_area = UNINITIALIZED_DOUBLE
  pm%biogenfc = UNINITIALIZED_DOUBLE
  pm%probdeg = UNINITIALIZED_INTEGER
  pm%bioidx = UNINITIALIZED_INTEGER
  pm%plasidx = UNINITIALIZED_INTEGER
  pm%output_start_time = 0.d0  ! [sec] default value
  pm%stoic_mat = UNINITIALIZED_DOUBLE
  
  call PMBaseInit(pm)
  
  PMWSSCreate => pm
  
end function PMWSSCreate

! *************************************************************************** !

function PMWSSWastePanelCreate()
  !
  ! Creates and initializes a waste panel type.
  !
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !
  
  implicit none
  
! LOCAL VARIABLES:
! ================
! PMWSSWastePanelCreate (output): pointer to new waste panel object
! panel: pointer to new waste panel object with shorter name
! ----------------------------------------------------------
  type(srcsink_panel_type), pointer :: PMWSSWastePanelCreate
  type(srcsink_panel_type), pointer :: panel
! ----------------------------------------------------------
  
  allocate(panel)
  
  nullify(panel%next)
  nullify(panel%region)
  nullify(panel%scaling_factor)
  nullify(panel%solids_production)
  nullify(panel%gas_generation_rate)
  nullify(panel%brine_generation_rate)
  nullify(panel%rxnrate_Fe_corrosion_inund)
  nullify(panel%rxnrate_Fe_corrosion_humid)
  nullify(panel%rxnrate_cell_biodeg_inund)
  nullify(panel%rxnrate_cell_biodeg_humid)
  nullify(panel%rxnrate_MgO_hyd_inund)
  nullify(panel%rxnrate_MgO_hyd_humid)
  nullify(panel%rxnrate_Fe_corrosion)
  nullify(panel%rxnrate_cell_biodeg)
  nullify(panel%rxnrate_FeOH2_sulf)
  nullify(panel%rxnrate_Fe_sulf)
  nullify(panel%rxnrate_MgO_hyd)
  nullify(panel%rxnrate_MgOH2_carb)
  nullify(panel%rxnrate_MgO_carb)
  nullify(panel%rxnrate_hydromag_conv)
  
  nullify(panel%rank_list)
  call PMWSSInventoryInit(panel%inventory)
  panel%name = ''
  panel%region_name = ''
  panel%inventory_name = ''
  panel%volume = UNINITIALIZED_DOUBLE
  panel%scale_by_volume = PETSC_FALSE
  panel%inundated_corrosion_rate = UNINITIALIZED_DOUBLE
  panel%humid_corrosion_rate = UNINITIALIZED_DOUBLE
  panel%inundated_biodeg_rate = UNINITIALIZED_DOUBLE
  panel%humid_biodeg_rate = UNINITIALIZED_DOUBLE
  panel%inundated_brucite_rate = UNINITIALIZED_DOUBLE
  panel%humid_brucite_rate = UNINITIALIZED_DOUBLE
  panel%F_NO3 = UNINITIALIZED_DOUBLE
  panel%F_SO4 = UNINITIALIZED_DOUBLE
  panel%RXH2S_factor = UNINITIALIZED_DOUBLE
  panel%RXCO2_factor = UNINITIALIZED_DOUBLE
  panel%RXH2_factor = UNINITIALIZED_DOUBLE
  panel%RXH2O_factor = UNINITIALIZED_DOUBLE
  panel%prev_dt = 0.d0
  panel%id = 0
  panel%myMPIgroup = 0
  panel%myMPIcomm = 0
  
  PMWSSWastePanelCreate => panel

end function PMWSSWastePanelCreate

! *************************************************************************** !

function PMWSSPreInventoryCreate()
  !
  ! Creates and initializes a waste panel pre-inventory.
  !
  ! Author: Jenn Frederick
  ! Date: 3/01/2017
  !
  
  implicit none
  
! LOCAL VARIABLES:
! ================
! PMWSSPreInventoryCreate (output): pointer to new pre-inventory object
! preinv: pointer to new pre-inventory object with shorter name
! ------------------------------------------------------------
  type(pre_inventory_type), pointer :: PMWSSPreInventoryCreate
  type(pre_inventory_type), pointer :: preinv
! ------------------------------------------------------------
  
  allocate(preinv)
  nullify(preinv%next)

  preinv%name = ''
  ! ALGEBRA parameters:
  preinv%ironchw = UNINITIALIZED_DOUBLE
  preinv%ironrhw = UNINITIALIZED_DOUBLE
  preinv%irncchw = UNINITIALIZED_DOUBLE
  preinv%irncrhw = UNINITIALIZED_DOUBLE
  preinv%cellchw = UNINITIALIZED_DOUBLE
  preinv%cellrhw = UNINITIALIZED_DOUBLE
  preinv%celcchw = UNINITIALIZED_DOUBLE
  preinv%celcrhw = UNINITIALIZED_DOUBLE
  preinv%celechw = UNINITIALIZED_DOUBLE
  preinv%celerhw = UNINITIALIZED_DOUBLE
  preinv%rubbchw = UNINITIALIZED_DOUBLE
  preinv%rubbrhw = UNINITIALIZED_DOUBLE
  preinv%rubcchw = UNINITIALIZED_DOUBLE
  preinv%rubcrhw = UNINITIALIZED_DOUBLE
  preinv%rubechw = UNINITIALIZED_DOUBLE
  preinv%ruberhw = UNINITIALIZED_DOUBLE
  preinv%plaschw = UNINITIALIZED_DOUBLE
  preinv%plasrhw = UNINITIALIZED_DOUBLE
  preinv%plscchw = UNINITIALIZED_DOUBLE
  preinv%plscrhw = UNINITIALIZED_DOUBLE
  preinv%plsechw = UNINITIALIZED_DOUBLE
  preinv%plserhw = UNINITIALIZED_DOUBLE
  preinv%plasfac = UNINITIALIZED_DOUBLE
  preinv%mgo_ef = UNINITIALIZED_DOUBLE
  preinv%vrepos = UNINITIALIZED_DOUBLE
  preinv%drum_conc = UNINITIALIZED_DOUBLE
  preinv%nitrate = UNINITIALIZED_DOUBLE
  preinv%sulfate = UNINITIALIZED_DOUBLE
  
  PMWSSPreInventoryCreate => preinv

end function PMWSSPreInventoryCreate

! *************************************************************************** !

subroutine PMWSSInventoryInit(inventory)
  !
  ! Initializes a waste panel inventory object.
  !
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !
  
  implicit none

! INPUT ARGUMENTS:
! ================
! inventory (input/output): waste panel inventory object
! ---------------------------------
  type(inventory_type) :: inventory
! ---------------------------------
  
! LOCAL VARIABLES:
! ================
! molar_mass: [kg/mol] molar mass of a chemical species object
! -----------------------
  PetscReal :: molar_mass
! -----------------------
  
  nullify(inventory%preinventory)
  inventory%name = ''
  
  molar_mass = MW_FE ! iron
  call PMWSSInitChemSpecies(inventory%Fe_s,molar_mass) 
  
  molar_mass = MW_FEOH2 ! iron hydroxide
  call PMWSSInitChemSpecies(inventory%FeOH2_s,molar_mass) 
  
  molar_mass = MW_CELL ! cellulosics/rubbers/plastics
  call PMWSSInitChemSpecies(inventory%BioDegs_s,molar_mass)  
  
  molar_mass = MW_FES ! iron sulfide
  call PMWSSInitChemSpecies(inventory%FeS_s,molar_mass)   
  
  molar_mass = MW_MGO ! magnesium oxide
  call PMWSSInitChemSpecies(inventory%MgO_s,molar_mass)  
  
  molar_mass = MW_MGOH2 ! magnesium hydroxide
  call PMWSSInitChemSpecies(inventory%MgOH2_s,molar_mass)                
  
  molar_mass = MW_HYDRO ! hydromagnesite
  call PMWSSInitChemSpecies(inventory%Mg5CO34OH24H2_s,molar_mass)  
  
  molar_mass = MW_MGCO3 ! magnesium carbonate
  call PMWSSInitChemSpecies(inventory%MgCO3_s,molar_mass)  
  
  inventory%Fe_in_panel = UNINITIALIZED_DOUBLE
  inventory%MgO_in_panel = UNINITIALIZED_DOUBLE
  inventory%Cellulose_in_panel = UNINITIALIZED_DOUBLE
  inventory%RubberPlas_in_panel = UNINITIALIZED_DOUBLE
  inventory%Biodegs_in_panel = UNINITIALIZED_DOUBLE
  inventory%Nitrate_in_panel = UNINITIALIZED_DOUBLE
  inventory%Sulfate_in_panel = UNINITIALIZED_DOUBLE
  inventory%drum_conc = UNINITIALIZED_DOUBLE

end subroutine PMWSSInventoryInit

! *************************************************************************** !

subroutine PMWSSInitChemSpecies(chem_species,molar_mass)
  !
  ! Initializes a waste panel inventory's chemical species.
  !
  ! Author: Jenn Frederick
  ! Date: 2/23/2017
  !
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! chem_species (input/output): chemical species object
! molar_mass (input): [kg/mol] molar mass of the chemical species object
! ---------------------------------------
  type(chem_species_type) :: chem_species
  PetscReal :: molar_mass
! ---------------------------------------
  
  nullify(chem_species%current_conc_mol)        ! [mol/m3]
  nullify(chem_species%current_conc_kg)         ! [kg/m3]
  nullify(chem_species%initial_conc_mol)        ! [mol/m3]
  nullify(chem_species%inst_rate)               ! [mol/m3/sec]
  chem_species%molar_mass = molar_mass          ! [kg/mol]
  chem_species%tot_mass_in_panel = 0.d0         ! [kg/panel-volume]
  
end subroutine PMWSSInitChemSpecies

! *************************************************************************** !

subroutine PMWSSSetRealization(this,realization)
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/02/2017
  !

  use Realization_Subsurface_class

  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! realization (input): pointer to subsurface realization object
! ----------------------------------------------------------
  class(pm_wipp_srcsink_type) :: this
  class(realization_subsurface_type), pointer :: realization
! ----------------------------------------------------------
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMWSSSetRealization

! *************************************************************************** !

subroutine PMWSSAssociateRegion(this,region_list)
  ! 
  ! Associates the waste panel to its assigned region via the REGION keyword.
  ! 
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !

  use Region_module
  use Option_module
  use String_module
  
  implicit none
  
! INPUT ARGUMENTS:
! =================
! this (input/output): wipp-srcsink process model object
! region_list (input): pointer to list object of region objects in the 
!    simulation
! ----------------------------------------------
  class(pm_wipp_srcsink_type) :: this
  type(region_list_type), pointer :: region_list
! ----------------------------------------------
  
! LOCAL VARIABLES:
! ================
! cur_region: pointer to current region object
! cur_waste_panel: pointer to current waste panel object
! option: pointer to option object
! -----------------------------------------------------
  type(region_type), pointer :: cur_region
  class(srcsink_panel_type), pointer :: cur_waste_panel
  type(option_type), pointer :: option
! -----------------------------------------------------
  
  option => this%option
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
      cur_region => region_list%first     
      do
        if (.not.associated(cur_region)) exit
        if (StringCompare(cur_region%name,cur_waste_panel%region_name)) then
          cur_waste_panel%region => cur_region
          exit
        endif
        cur_region => cur_region%next
      enddo      
      if (.not.associated(cur_waste_panel%region)) then
        option%io_buffer = 'WASTE_PANEL REGION ' // &
                           trim(cur_waste_panel%region_name) // ' not found.'
        call printErrMsg(option)
      endif
    cur_waste_panel => cur_waste_panel%next
  enddo
  
end subroutine PMWSSAssociateRegion

! *************************************************************************** !

subroutine PMWSSAssociateInventory(this)
  ! 
  ! Associates the waste panel to its assigned inventory via INVENTORY keyword.
  ! 
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !

  use Option_module
  use String_module
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
! -----------------------------------
  
! LOCAL VARIABLES:
! ================
! cur_preinventory: pointer to current preinventory object
! cur_waste_panel: pointer to current waste panel object
! option: pointer to option object
! -----------------------------------------------------
  type(pre_inventory_type), pointer :: cur_preinventory
  class(srcsink_panel_type), pointer :: cur_waste_panel
  type(option_type), pointer :: option
! -----------------------------------------------------
  
  option => this%option
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    cur_preinventory => this%pre_inventory_list     
    do
      if (.not.associated(cur_preinventory)) exit
      if (StringCompare(cur_preinventory%name, &
                        cur_waste_panel%inventory_name)) then
        call PMWSSCopyPreInvToInv(cur_preinventory,cur_waste_panel%inventory)
        exit
      endif
      cur_preinventory => cur_preinventory%next
    enddo      
    if (.not.associated(cur_waste_panel%inventory%preinventory)) then
      option%io_buffer = 'WASTE_PANEL INVENTORY ' // &
                         trim(cur_waste_panel%inventory_name) // ' not found.'
      call printErrMsg(option)
    endif
    if (cur_waste_panel%scale_by_volume .and. &
        Uninitialized(cur_waste_panel%inventory%preinventory%vrepos)) then
      option%io_buffer = 'ERROR: WASTE_PANEL ' // trim(cur_waste_panel%name) &
                        // ' indicated SCALE_BY VOLUME = YES, but keyword &
                        &VREPOS was not given in INVENTORY ' // &
                        trim(cur_waste_panel%inventory%preinventory%name) // '.'
      call printErrMsg(option)
    endif
    cur_waste_panel => cur_waste_panel%next
  enddo
  
end subroutine PMWSSAssociateInventory

! *************************************************************************** !

subroutine PMWSSCopyPreInvToInv(preinventory,inventory)
  ! 
  ! Copies information from a pre-inventory to an inventory object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! preinventory (input/output): pointer to pre-inventory object
! inventory (input): inventory object
! -------------------------------------------------
  type(pre_inventory_type), pointer :: preinventory
  type(inventory_type) :: inventory
! -------------------------------------------------
  
  inventory%name = preinventory%name
  inventory%drum_conc = preinventory%drum_conc
  inventory%preinventory => preinventory
  
end subroutine PMWSSCopyPreInvToInv

! *************************************************************************** !

subroutine PMWSSSetRegionScaling(this,waste_panel)
  ! 
  ! Calculates and sets the scaling factor vector for each of the waste panels
  ! that have assigned regions. It assumes the volume of the cells that make up
  ! the region do not change over the course of the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !

  use Material_Aux_class
  use Grid_module

  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! waste_panel (input/output): pointer to waste panel object
! ------------------------------------------------
  class(pm_wipp_srcsink_type) :: this
  type(srcsink_panel_type), pointer :: waste_panel
! ------------------------------------------------
  
! LOCAL VARIABLES:
! ================
! material_auxvars(:): pointer to material auxvar object which stores 
!    grid cell volume indexed by the ghosted cell id
! grid: pointer to grid object
! k: [-] looping index
! local_id: [-] local grid cell id
! ghosted_id: [-] ghosted grid cell id
! total_volume_local: [m3] total local volume of grid cells in a waste panel
! total_volume_global: [m3] total global volume of grid cells in a waste panel
! ierr: PETSc error integer
! -----------------------------------------------------------
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type), pointer :: grid
  PetscInt :: k
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: total_volume_local
  PetscReal :: total_volume_global
  PetscErrorCode :: ierr
! -----------------------------------------------------------
  
  material_auxvars => this%realization%patch%aux%Material%auxvars
  grid => this%realization%patch%grid
  allocate(waste_panel%scaling_factor(waste_panel%region%num_cells))
  total_volume_local = 0.d0
  total_volume_global = 0.d0
  
  ! scale by cell volume
  do k = 1,waste_panel%region%num_cells
    local_id = waste_panel%region%cell_ids(k)
    ghosted_id = grid%nL2G(local_id)
    waste_panel%scaling_factor(k) = material_auxvars(ghosted_id)%volume ! [m^3]
    total_volume_local = total_volume_local &
                         + material_auxvars(ghosted_id)%volume          ! [m^3]
  enddo
  call MPI_Allreduce(total_volume_local,total_volume_global,ONE_INTEGER_MPI, &
              MPI_DOUBLE_PRECISION,MPI_SUM,waste_panel%myMPIcomm,ierr)
  waste_panel%scaling_factor = waste_panel%scaling_factor/total_volume_global 
  waste_panel%volume = total_volume_global
  
end subroutine PMWSSSetRegionScaling

! *************************************************************************** !

subroutine PMWSSRead(this,input)
  !
  ! Reads input file parameters for the WIPP source/sink process model.
  !
  ! Author: Jenn Frederick
  ! Date: 1/31/2017
  !
  
  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! input (input/output): pointer to input object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
  type(input_type), pointer :: input
! -----------------------------------
  
! LOCAL VARIABLES:
! ================
! option: pointer to option object
! word, word2: temporary strings
! double: temporary double precision number
! error_string#: temporary string
! new_waste_panel: pointer to waste panel object for new allocation
! cur_waste_panel: pointer to current waste panel object
! new_inventory: pointer to inventory object for new allocation
! cur_preinventory: pointer to current pre-inventory object
! num_errors: [-] number of errors integer
! rxn_num: [-] looping index integer for the reaction number
! spec_num: [-] looping index integer for the species number
! added: temporary Boolean
! ----------------------------------------------------------------------------
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word, word2
  PetscReal :: double
  character(len=MAXSTRINGLENGTH) :: error_string, error_string2, error_string3
  type(srcsink_panel_type), pointer :: new_waste_panel
  type(srcsink_panel_type), pointer :: cur_waste_panel
  type(pre_inventory_type), pointer :: new_inventory
  type(pre_inventory_type), pointer :: cur_preinventory
  PetscInt :: num_errors
  PetscInt :: rxn_num, spec_num
  PetscBool :: added
! ----------------------------------------------------------------------------
  
  option => this%option
  input%ierr = 0
  option%io_buffer = 'pflotran card:: WIPP_SOURCE_SINK'
  call printMsg(option)
  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    num_errors = 0
    error_string = 'WIPP_SOURCE_SINK'
    call StringToUpper(word)
    select case(trim(word))
    !-----------------------------------------
      case('ALPHARXN')
        call InputReadDouble(input,option,this%alpharxn)
        call InputErrorMsg(input,option,'ALPHARXN',error_string)
      case('SOCMIN')
        call InputReadDouble(input,option,this%smin)
        call InputErrorMsg(input,option,'SOCMIN',error_string)
      case('SAT_WICK')
        call InputReadDouble(input,option,this%satwick)
        call InputErrorMsg(input,option,'wicking saturation parameter &
                           &(SAT_WICK)',error_string)
      case('CORRMCO2')
        call InputReadDouble(input,option,this%corrmco2)
        call InputErrorMsg(input,option,'inundated steel corrosion rate &
                           &(CORRMCO2)',error_string)
      case('HUMCORR')
        call InputReadDouble(input,option,this%humcorr)
        call InputErrorMsg(input,option,'humid steel corrosion rate &
                           &(HUMCORR)',error_string)
      case('GRATMICI')
        call InputReadDouble(input,option,this%gratmici)
        call InputErrorMsg(input,option,'inundated biodegradation rate for &
                           &cellulose (GRATMICI)',error_string)
      case('GRATMICH')
        call InputReadDouble(input,option,this%gratmich)
        call InputErrorMsg(input,option,'humid diodegradation rate for &
                           &cellulose (GRATMICH)',error_string)
      case('BRUCITEC','BRUCITES','BRUCITEI')
        call InputReadDouble(input,option,this%brucitei)
        call InputErrorMsg(input,option,'MgO inundated hydration rate in &
                           &Castile or Salado brine (BRUCITE[C/S])', &
                           error_string)
      case('BRUCITEH')
        call InputReadDouble(input,option,this%bruciteh)
        call InputErrorMsg(input,option,'MgO humid hydration rate (BRUCITEH)', &
                           error_string)
      case('HYMAGCON')
        call InputReadDouble(input,option,this%hymagcon_rate)
        call InputErrorMsg(input,option,'hydromagnesite to magnesite &
                           &conversion rate (HYMAGCON)',error_string)
      case('ASDRUM')
        call InputReadDouble(input,option,this%drum_surface_area)
        call InputErrorMsg(input,option,'surface area of corrodable metal &
                           &per drum (ASDRUM)',error_string)
      case('BIOGENFC')
        call InputReadDouble(input,option,this%biogenfc)
        call InputErrorMsg(input,option,'probability of attaining sampled &
                           &microbial gas generation rates (BIOGENFC)', &
                           error_string)
      case('PROBDEG')
        call InputReadInt(input,option,this%probdeg)
        call InputErrorMsg(input,option,'flag: (PROBDEG)',error_string)
      case('OUTPUT_START_TIME')
        call InputReadDouble(input,option,double)
        call InputErrorMsg(input,option,'OUTPUT_START_TIME',error_string)
        call InputReadAndConvertUnits(input,double,'sec',trim(error_string) &
                                      // ',OUTPUT_START_TIME units',option)
        this%output_start_time = double
    !-----------------------------------------
    !-----------------------------------------
      case('WASTE_PANEL')
        error_string = trim(error_string) // ',WASTE_PANEL'
        allocate(new_waste_panel)
        new_waste_panel => PMWSSWastePanelCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        new_waste_panel%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_waste_panel%name)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
          !-----------------------------------
            case('REGION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'region assignment',error_string)
              new_waste_panel%region_name = trim(word)
          !-----------------------------------
            case('INVENTORY')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'inventory assignment', &        
                                 error_string)
              new_waste_panel%inventory_name = trim(word)
          !-----------------------------------
            case('SCALE_BY_VOLUME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'SCALE_BY_VOLUME',error_string)
              call StringToUpper(word)
              select case(trim(word))
                case('YES')
                  new_waste_panel%scale_by_volume = PETSC_TRUE
                case('NO')
                  new_waste_panel%scale_by_volume = PETSC_FALSE
                case default
                  call InputKeywordUnrecognized(word,'SCALE_BY_VOLUME (must &
                  &be "YES" or "NO")',option)
              end select
          !-----------------------------------    
            case default
              call InputKeywordUnrecognized(word,error_string,option)
          !----------------------------------- 
          end select
        enddo
        ! error messages ---------------------
        if (new_waste_panel%region_name == '') then
          option%io_buffer = 'ERROR: REGION must be specified in the ' // &
                 trim(error_string) // ' block. WASTE_PANEL name "' // &
                 trim(new_waste_panel%name) // '".'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (new_waste_panel%inventory_name == '') then
          option%io_buffer = 'ERROR: INVENTORY must be specified in the ' // &
                 trim(error_string) // ' block. WASTE_PANEL name "' // &
                 trim(new_waste_panel%name) // '".'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (num_errors > 0) then
          write(option%io_buffer,*) num_errors
          option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                             &WIPP_SOURCE_SINK,WASTE_PANEL block. See above.'
          call printErrMsg(option)
        endif
        added = PETSC_FALSE
        if (.not.associated(this%waste_panel_list)) then
          this%waste_panel_list => new_waste_panel
        else
          cur_waste_panel => this%waste_panel_list
          do
            if (.not.associated(cur_waste_panel)) exit
            if (.not.associated(cur_waste_panel%next)) then
              cur_waste_panel%next => new_waste_panel
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_waste_panel => cur_waste_panel%next
          enddo
        endif
        nullify(new_waste_panel)
    !-----------------------------------------
    !-----------------------------------------
      case('STOICHIOMETRIC_MATRIX')
        error_string = trim(error_string) // ',STOICHIOMETRIC_MATRIX'
        do rxn_num = 1,8
          call InputReadPflotranString(input,option)
          do spec_num = 1,10
            write(word,'(i1)') rxn_num
            write(word2,'(i2)') spec_num
            call InputReadDouble(input,option,this%stoic_mat(rxn_num,spec_num))
            call InputErrorMsg(input,option,'ROW ' // trim(adjustl(word)) // &
                               ', COL ' // trim(adjustl(word2)),error_string)
          enddo
        enddo
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        !if (InputCheckExit(input,option)) exit
    !-----------------------------------------
    !-----------------------------------------
      case('INVENTORY')
        error_string = trim(error_string) // ',INVENTORY'
        allocate(new_inventory)
        new_inventory => PMWSSPreInventoryCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        new_inventory%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_inventory%name)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
          !-----------------------------------
            case('SOLIDS','SOLID')
              error_string2 = trim(error_string) // ',SOLIDS'
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'keyword',error_string2)
                call StringToUpper(word)
                select case(trim(word))
                !-----------------------------
                  case('IRONCHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'IRONCHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                          trim(error_string2) // ',IRONCHW',option)
                    new_inventory%ironchw = double
                !-----------------------------
                  case('IRONRHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'IRONRHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                          trim(error_string2) // ',IRONRHW',option)
                    new_inventory%ironrhw = double
                !-----------------------------
                  case('IRNCCHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'IRNCCHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                          trim(error_string2) // ',IRNCCHW',option)
                    new_inventory%irncchw = double
                !-----------------------------
                  case('IRNCRHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'IRNCRHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                          trim(error_string2) // ',IRNCRHW',option)
                    new_inventory%irncrhw = double
                !-----------------------------
                  case('MGO_EF')
                    call InputReadDouble(input,option,new_inventory%mgo_ef)
                    call InputErrorMsg(input,option,'MGO_EF',error_string2)
                !-----------------------------
                  case('CELLCHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'CELLCHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',CELLCHW',option)
                    new_inventory%cellchw = double
                !-----------------------------
                  case('CELLRHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'CELLRHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',CELLRHW',option)
                    new_inventory%cellrhw = double
                !-----------------------------
                  case('CELCCHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'CELCCHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',CELCCHW',option)
                    new_inventory%celcchw = double
                !-----------------------------
                  case('CELCRHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'CELCRHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',CELCRHW',option)
                    new_inventory%celcrhw = double
                !-----------------------------
                  case('CELECHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'CELECHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',CELECHW',option)
                    new_inventory%celechw = double
                !-----------------------------
                  case('CELERHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'CELERHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',CELERHW',option)
                    new_inventory%celerhw = double
                !-----------------------------
                  case('RUBBCHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'RUBBCHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',RUBBCHW',option)
                    new_inventory%rubbchw = double
                !-----------------------------
                  case('RUBBRHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'RUBBRHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',RUBBRHW',option)
                    new_inventory%rubbrhw = double
                !-----------------------------
                  case('RUBCCHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'RUBCCHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',RUBCCHW',option)
                    new_inventory%rubcchw = double
                !-----------------------------
                  case('RUBCRHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'RUBCRHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',RUBCRHW',option)
                    new_inventory%rubcrhw = double
                !-----------------------------
                  case('RUBECHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'RUBECHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',RUBECHW',option)
                    new_inventory%rubechw = double
                !-----------------------------
                  case('RUBERHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'RUBERHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',RUBERHW',option)
                    new_inventory%ruberhw = double
                !-----------------------------
                  case('PLASCHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'PLASCHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',PLASCHW',option)
                    new_inventory%plaschw = double
                !-----------------------------
                  case('PLASRHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'PLASRHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',PLASRHW',option)
                    new_inventory%plasrhw = double    
                !-----------------------------
                  case('PLSCCHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'PLSCCHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',PLSCCHW',option)
                    new_inventory%plscchw = double
                !-----------------------------
                  case('PLSCRHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'PLSCRHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',PLSCRHW',option)
                    new_inventory%plscrhw = double
                !-----------------------------
                  case('PLSECHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'PLSECHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',PLSECHW',option)
                    new_inventory%plsechw = double
                !-----------------------------
                  case('PLSERHW')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'PLSERHW',error_string2)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string2) // ',PLSERHW',option)
                    new_inventory%plserhw = double
                !-----------------------------
                  case('PLASFAC')
                    call InputReadDouble(input,option,new_inventory%plasfac)
                    call InputErrorMsg(input,option,'PLASFAC',error_string2)
                !-----------------------------------
                  case('DRMCONC')
                    call InputReadDouble(input,option,new_inventory%drum_conc)
                    call InputErrorMsg(input,option,'DRMCONC',error_string2)
                !-----------------------------
                  case default
                    call InputKeywordUnrecognized(word,error_string2,option)
                !-----------------------------
                end select
              enddo
          !-----------------------------------
            case('AQUEOUS')
              error_string3 = trim(error_string) // ',AQUEOUS'
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'keyword',error_string3)
                call StringToUpper(word)
                select case(trim(word))
                !-----------------------------
                  case('NITRATE')
                    call InputReadDouble(input,option,new_inventory%nitrate)
                    call InputErrorMsg(input,option,'initial nitrate moles &
                                       &(NITRATE)',error_string3)
                !-----------------------------
                  case('SULFATE')
                    call InputReadDouble(input,option,new_inventory%sulfate)
                    call InputErrorMsg(input,option,'initial sulfate moles &
                                       &(SULFATE)',error_string3)
                !-----------------------------
                  case default
                    call InputKeywordUnrecognized(word,error_string3,option)
                !-----------------------------
                end select
              enddo
          !-----------------------------------
            case('VREPOS')
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'VREPOS',error_string)
              call InputReadAndConvertUnits(input,double,'m^3', &
                      trim(error_string) // ',VREPOS volume units',option)
              new_inventory%vrepos = double
          !-----------------------------------    
            case default
              call InputKeywordUnrecognized(word,error_string,option)
          !----------------------------------- 
          end select
        enddo
        ! error messages ---------------------
        if (Uninitialized(new_inventory%drum_conc)) then
          option%io_buffer = 'ERROR: Number of metal drums per m3 of waste &
                        &area must be specified using the SOLIDS,DRMCONC &
                        &keyword in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
      !----- IRON -----!
        if (Uninitialized(new_inventory%ironchw)) then
          option%io_buffer = 'ERROR: Initial mass of Fe-based material in CH &
                        &waste must be specified using the SOLIDS,IRONCHW &
                        &keyword in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%ironrhw)) then
          option%io_buffer = 'ERROR: Initial mass of Fe-based material in RH &
                        &waste must be specified using the SOLIDS,IRONRHW &
                        &keyword in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%irncchw)) then
          option%io_buffer = 'ERROR: Initial mass of Fe containers for CH &
                        &waste must be specified using the SOLIDS,IRNCCHW &
                        &keyword in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%irncrhw)) then
          option%io_buffer = 'ERROR: Initial mass of Fe containers for RH &
                        &waste must be specified using the SOLIDS,IRNCRHW &
                        &keyword in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
      !----- MGO -----!
        if (Uninitialized(new_inventory%mgo_ef)) then
          option%io_buffer = 'ERROR: MgO excess factor must be &
                        &specified using the SOLIDS,MGO_EF keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
      !----- CELLULOSICS -----!
        if (Uninitialized(new_inventory%cellchw)) then
          option%io_buffer = 'ERROR: Initial cellulosics mass for CH waste &
                        &must be specified using the SOLIDS,CELLCHW keyword &
                        &in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%cellrhw)) then
          option%io_buffer = 'ERROR: Initial cellulosics mass for RH waste &
                        &must be specified using the SOLIDS,CELLRHW keyword &
                        &in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%celcchw)) then
          option%io_buffer = 'ERROR: Initial cellulosics mass in container &
                        &materials for CH waste must be specified using the &
                        &SOLIDS,CELCCHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%celcrhw)) then
          option%io_buffer = 'ERROR: Initial cellulosics mass in container &
                        &materials for RH waste must be specified using the &
                        &SOLIDS,CELCRHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%celechw)) then
          option%io_buffer = 'ERROR: Initial cellulosics mass in emplacement &
                        &materials for CH waste must be specified using the &
                        &SOLIDS,CELECHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%celerhw)) then
          option%io_buffer = 'ERROR: Initial cellulosics mass in emplacement &
                        &materials for RH waste must be specified using the &
                        &SOLIDS,CELERHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
      !----- RUBBER -----!
        if (Uninitialized(new_inventory%rubbchw)) then
          option%io_buffer = 'ERROR: Initial rubber mass for CH waste must be &
                        &specified using the SOLIDS,RUBBCHW keyword &
                        &in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%rubbrhw)) then
          option%io_buffer = 'ERROR: Initial rubber mass for CH waste must be &
                        &specified using the SOLIDS,RUBBRHW keyword &
                        &in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%rubcchw)) then
          option%io_buffer = 'ERROR: Initial rubber mass in container &
                        &materials for CH waste must be specified using the &
                        &SOLIDS,RUBCCHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%rubcrhw)) then
          option%io_buffer = 'ERROR: Initial rubber mass in container &
                        &materials for RH waste must be specified using the &
                        &SOLIDS,RUBCRHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%rubechw)) then
          option%io_buffer = 'ERROR: Initial rubber mass in emplacement &
                        &materials for CH waste must be specified using the &
                        &SOLIDS,RUBECHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%ruberhw)) then
          option%io_buffer = 'ERROR: Initial rubber mass in emplacement &
                        &materials for RH waste must be specified using the &
                        &SOLIDS,RUBERHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
      !----- PLASTICS -----!
        if (Uninitialized(new_inventory%plaschw)) then
          option%io_buffer = 'ERROR: Initial plastics mass for CH waste must &
                        &be specified using the SOLIDS,PLASCHW keyword &
                        &in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%plasrhw)) then
          option%io_buffer = 'ERROR: Initial plastics mass for CH waste must &
                        &be specified using the SOLIDS,PLASRHW keyword &
                        &in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%plscchw)) then
          option%io_buffer = 'ERROR: Initial plastics mass in container &
                        &materials for CH waste must be specified using the &
                        &SOLIDS,PLSCCHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%plscrhw)) then
          option%io_buffer = 'ERROR: Initial plastics mass in container &
                        &materials for RH waste must be specified using the &
                        &SOLIDS,PLSCRHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%plsechw)) then
          option%io_buffer = 'ERROR: Initial plastics mass in emplacement &
                        &materials for CH waste must be specified using the &
                        &SOLIDS,PLSECHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%plserhw)) then
          option%io_buffer = 'ERROR: Initial plastics mass in emplacement &
                        &materials for RH waste must be specified using the &
                        &SOLIDS,PLSERHW keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%plasfac)) then
          option%io_buffer = 'ERROR: Mass ratio of plastics to equivalent &
                        &carbon must be specified using the SOLIDS,PLASFAC &
                        &keyword in the WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
      !----- AQUEOUS -----!
        if (Uninitialized(new_inventory%nitrate)) then
          option%io_buffer = 'ERROR: Initial nitrate moles inventory must &
                        &be specified using the AQUEOUS,NITRATE keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_inventory%sulfate)) then
          option%io_buffer = 'ERROR: Initial sulfate moles inventory must &
                        &be specified using the AQUEOUS,SULFATE keyword in the &
                        &WIPP_SOURCE_SINK,INVENTORY ' // &
                        trim(new_inventory%name) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
      !----- END COUNT -----!
        if (num_errors > 0) then
          write(option%io_buffer,*) num_errors
          option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                             &WIPP_SOURCE_SINK,INVENTORY block. See above.'
          call printErrMsg(option)
        endif
        added = PETSC_FALSE
        if (.not.associated(this%pre_inventory_list)) then
          this%pre_inventory_list => new_inventory
        else
          cur_preinventory => this%pre_inventory_list
          do
            if (.not.associated(cur_preinventory)) exit
            if (.not.associated(cur_preinventory%next)) then
              cur_preinventory%next => new_inventory
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_preinventory => cur_preinventory%next
          enddo
        endif
        nullify(new_inventory)
    !-----------------------------------------
      case default
        call InputKeywordUnrecognized(word,'WIPP_SOURCE_SINK',option)
    !-----------------------------------------
    end select  
  enddo
  
  if (.not.associated(this%waste_panel_list)) then
    option%io_buffer = 'ERROR: At least one WASTE_PANEL must be specified &
                       &in the WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (.not.associated(this%pre_inventory_list)) then
    option%io_buffer = 'ERROR: At least one INVENTORY must be specified in &
                       &the WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%alpharxn)) then
    option%io_buffer = 'ERROR: ALPHARXN must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%smin)) then
    option%io_buffer = 'ERROR: SOCMIN must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%satwick)) then
    option%io_buffer = 'ERROR: SAT_WICK (wicking saturation parameter) must be &
                       &specified in the WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%gratmici)) then
    option%io_buffer = 'ERROR: GRATMICI (inundated biodegradation rate for &
                       &cellulose) must be specified in the WIPP_SOURCE_SINK &
                       &block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%brucitei)) then
    option%io_buffer = 'ERROR: BRUCITE[C/S] (MgO inundated hydration rate in &
                       &Castile or Salado brine) must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%corrmco2)) then
    option%io_buffer = 'ERROR: CORRMCO2 (inundated steel corrosion rate) must &
                       &be specified in the WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%gratmich)) then
    option%io_buffer = 'ERROR: GRATMICH (humid biodegradation rate for &
                       &cellulose) must be specified in the WIPP_SOURCE_SINK &
                       &block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%bruciteh)) then
    option%io_buffer = 'ERROR: BRUCITEH (MgO humid hydration rate) must be &
                       &specified in the WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%humcorr)) then
    option%io_buffer = 'ERROR: HUMCORR (humid steel corrosion rate) must be &
                       &specified in the WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%hymagcon_rate)) then
    option%io_buffer = 'ERROR: HYMAGCON (hydromagnesite to magnesite &
                       &conversion rate) must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%drum_surface_area)) then
    option%io_buffer = 'ERROR: ASDRUM (metal drum surface area) &
                       &must be specified in the WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%biogenfc)) then
    option%io_buffer = 'ERROR: BIOGENFC (microbial gas generation probability) &
                       &must be specified in the WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%probdeg)) then
    option%io_buffer = 'ERROR: PROBDEG (biodegradation and/or plastics &
                       &inclusion flag) must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  do rxn_num = 1,8
    do spec_num = 1,10
      if (Uninitialized(this%stoic_mat(rxn_num,spec_num))) then
        write(word,'(i1)') rxn_num
        write(word2,'(i2)') spec_num
        option%io_buffer = 'ERROR: There is a missing value in the &
                           &STOICHIOMETRIC_MATRIX at ROW ' // &
                           trim(adjustl(word)) // ', COL ' // &
                           trim(adjustl(word2)) // '. This may mean you did &
                           &not provide enough values for the matrix. The &
                           &matrix should be sized 8x10.'
      call printMsg(option)
      num_errors = num_errors + 1
      endif
    enddo
  enddo
  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the WIPP_SOURCE_SINK block. See above.'
    call printErrMsg(option)
  endif
  
  call PMWSSAssociateInventory(this)
  
end subroutine PMWSSRead

! *************************************************************************** !

subroutine PMWSSProcessAfterRead(this,waste_panel)
  !
  ! After reading input parameters, ALGEBRA processing is done to get final 
  ! input parameters required.
  !
  ! Author: Jenn Frederick
  ! Date: 3/10/2017
  !
  
  use Option_module

  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! waste_panel (input/output): pointer to waste panel object
! ------------------------------------------------
  class(pm_wipp_srcsink_type) :: this
  type(srcsink_panel_type), pointer :: waste_panel
! ------------------------------------------------
  
! LOCAL VARIABLES:
! ================
! preinventory: pointer to pre-inventory object
! inventory: pointer to inventory object
! vol_scaling_factor: [-] grid cell volume scaling factor
! MOL_NO3: [mol] moles of nitrate
! F_NO3: [-] fraction of carbon consumed through the denitrification reaction
! MAX_C, A1, A2: intermediate calculation parameters
! D_c: [kg/m3] mass concentration of biodegradable materials
! D_m: [kg/m3] mass concentration of MgO
! D_s: [m2/m3] surface area concentration of iron steel
! -------------------------------------------------
  type(pre_inventory_type), pointer :: preinventory
  type(inventory_type), pointer :: inventory
  PetscReal :: vol_scaling_factor
  PetscReal :: MOL_NO3
  PetscReal :: F_NO3
  PetscReal :: MAX_C, A1, A2         
  PetscReal :: D_c 
  PetscReal :: D_m 
  PetscReal :: D_s
! -------------------------------------------------
  
  !-----PROBDEG-calculations-----------------------------------------------
  !-----(see Table PA-6, section PA-4.2.5)---------------------------------
  select case(this%probdeg)
  !--no-biodegradation--no-rubbers/plastics-degradation----
    case(0)
      this%bioidx = 0
      this%plasidx = 0
  !--yes-biodegradation--no-rubbers/plastics-degradation----
    case(1)
      this%bioidx = 1
      this%plasidx = 0
  !--yes-biodegradation--yes-rubbers/plastics-degradation----
    case(2)
      this%bioidx = 1
      this%plasidx = 1
    case default
      this%option%io_buffer = 'WIPP_SOURCE_SINK,PROBDEG values: 0,1,2 only.'
      call printErrMsg(this%option)
  end select
  
  preinventory => waste_panel%inventory%preinventory
  inventory => waste_panel%inventory
  !-----inventory-totals-------------------------------------units---------
  inventory%Fe_in_panel = &                                ! [kg]
        preinventory%ironchw + preinventory%ironrhw + &    ! [kg]
        preinventory%irncchw + preinventory%irncrhw        ! [kg]
  inventory%Cellulose_in_panel = &                         ! [kg]
        preinventory%cellchw + preinventory%cellrhw + &    ! [kg]
        preinventory%celcchw + preinventory%celcrhw + &    ! [kg]
        preinventory%celechw + preinventory%celerhw        ! [kg]
  inventory%RubberPlas_in_panel = &                        ! [kg]
       (preinventory%rubbchw + preinventory%rubbrhw + &    ! [kg]
        preinventory%rubcchw + preinventory%rubcrhw + &    ! [kg]
        preinventory%rubechw + preinventory%ruberhw) + &   ! [kg]
        preinventory%plasfac * &                           ! [-]
       (preinventory%plaschw + preinventory%plasrhw + &    ! [kg]
        preinventory%plscchw + preinventory%plscrhw + &    ! [kg]
        preinventory%plsechw + preinventory%plserhw)       ! [kg]
  inventory%Biodegs_in_panel = &                           ! [kg]
        inventory%Cellulose_in_panel + &                   ! [kg]
       (inventory%RubberPlas_in_panel*this%plasidx)        ! [kg]
  inventory%MgO_in_panel = &                               ! [kg]
       (inventory%Cellulose_in_panel + &                   ! [kg]
        inventory%RubberPlas_in_panel) * &                 ! [kg]
        preinventory%mgo_ef * MW_MGO / MW_CELL             ! [-]
  inventory%Nitrate_in_panel = &                           ! [mol]
        preinventory%nitrate                               ! [mol]
  inventory%Sulfate_in_panel = &                           ! [mol]
        preinventory%sulfate                               ! [mol]
  !-----scale-inventory-totals-------------------------------units---------
  if (waste_panel%scale_by_volume) then
    vol_scaling_factor = waste_panel%volume / &            ! [m3]
                         preinventory%vrepos               ! [m3]
    inventory%Fe_in_panel = &                              ! [kg] 
        inventory%Fe_in_panel*vol_scaling_factor           ! [kg]
    inventory%Cellulose_in_panel = &                       ! [kg] 
        inventory%Cellulose_in_panel*vol_scaling_factor    ! [kg]
    inventory%RubberPlas_in_panel = &                      ! [kg] 
        inventory%RubberPlas_in_panel*vol_scaling_factor   ! [kg]
    inventory%Biodegs_in_panel = &                         ! [kg] 
        inventory%Biodegs_in_panel*vol_scaling_factor      ! [kg]
    inventory%MgO_in_panel = &                             ! [kg] 
        inventory%MgO_in_panel*vol_scaling_factor          ! [kg]
    inventory%Nitrate_in_panel = &                         ! [kg] 
        inventory%Nitrate_in_panel*vol_scaling_factor      ! [kg]
    inventory%Sulfate_in_panel = &                         ! [kg] 
        inventory%Sulfate_in_panel*vol_scaling_factor      ! [kg]
  endif
  !-----(see equation PA.76, PA.75, PA.93, section PA-4.2.5)---------------
  !-----mass-concentrations----------------------------------units---------
  D_c = inventory%Biodegs_in_panel / &                     ! [kg]
        waste_panel%volume                                 ! [m3]
  D_s = this%drum_surface_area * &                         ! [m2]
        inventory%drum_conc                                ! [-/m3]
  D_m = inventory%MgO_in_panel / &                         ! [kg]
        waste_panel%volume                                 ! [m3]
  !-------------------------------------------------------------------------
  !-----(see equation PA.67, section PA-4.2.5)------------------------------
  !-----anoxic-iron-corrosion--------------------------------units----------
  waste_panel%inundated_corrosion_rate = &                 ! [mol-Fe/m3/sec]
        this%corrmco2 * &                                  ! [m/s]
        D_s * &                                            ! [m2/m3]
        DN_FE / &                                          ! [kg/m3]
        MW_FE                                              ! [kg/mol]
  waste_panel%humid_corrosion_rate = &                     ! [mol-Fe/m3/sec]
        this%humcorr * &                                   ! [m/s]
        D_s * &                                            ! [m2/m3]
        DN_FE / &                                          ! [kg/m3]
        MW_FE                                              ! [kg/mol]
  !-----(see equation PA.69, section PA-4.2.5)------------------------------
  !-----biodegradation------------------------------------units-------------
  waste_panel%inundated_biodeg_rate = &                 ! [mol-cell/m3/sec]
        this%gratmici * &                               ! [mol-cell/kg/sec]
        D_c * &                                         ! [kg/m3]
        this%biogenfc * &                               ! [-]
        this%bioidx                                     ! [-]
  waste_panel%humid_biodeg_rate = &                     ! [mol-cell/m3/sec]
        this%gratmich * &                               ! [mol-cell/kg/sec]
        D_c * &                                         ! [kg/m3]
        this%biogenfc * &                               ! [-]
        this%bioidx                                     ! [-]
  !-----(see equation PA.86, PA.87, PA.88, section PA-4.2.5)----------------
  !-----iron-sulfidation----------------------------------units-------------
  MOL_NO3 = inventory%Nitrate_in_panel                  ! [mol]
  A1 = inventory%Biodegs_in_panel / MW_CELL             ! [mol]
  A2 = this%gratmici * &                                ! [mol/kg/sec]
       (inventory%Biodegs_in_panel) * &                 ! [kg]
       (31556930.d0) * 10000.d0                         ! [sec/year]*[year]
  MAX_C = min(A1,A2)                                    ! [mol]
  F_NO3 = MOL_NO3 * (6.d0/4.8d0) / MAX_C                ! [-]
  F_NO3 = min(F_NO3,1.0)                                ! [-]
  waste_panel%F_NO3 = F_NO3                             ! [-]
  waste_panel%F_SO4 = 1.d0 - F_NO3                      ! [-]
  
  waste_panel%RXH2_factor = &                           ! [mol-H2/mol-cell]
        waste_panel%F_NO3*(2.4d0/6.0d0) + &
        waste_panel%F_SO4*(3.0d0/6.0d0)
  waste_panel%RXH2O_factor = &                          ! [mol-H2O/mol-cell]
        waste_panel%F_NO3*(7.4d0/6.0d0) + &
        waste_panel%F_SO4*(5.0d0/6.0d0)
  
  ! algebra/pre-brag/bragflo use the H2 (i.e. total gas) value for H2S
  waste_panel%RXH2S_factor = waste_panel%RXH2_factor   ! [mol-H2S/mol-cell]
  ! this is the correct, H2S-only, value
  !waste_panel%RXH2S_factor = waste_panel%F_SO4*(3.0d0/6.0d0) ! [mol-H2S/mol-cell]
  
  ! bragflo apparently sets RXH2O_factor to zero
  ! algebra sets STCO_22 as SMIC_H20, which is RXH2O_factor
  ! also, it includes salt weight in the brine (water) weight in the
  ! brine generation rate output
  waste_panel%RXH2O_factor = 0.0d0
  
  ! BRAGFLO User's Manual Eq. 155, based on Eqs. 145 & 146 stoichiometry 
  waste_panel%RXCO2_factor = 1.0d0                     ! [mol-CO2/mol-cell]
  
  !-----(see equation PA.73, section PA-4.2.5)------------------------------
  !-----MgO-hydration-------------------------------------units-------------
  waste_panel%inundated_brucite_rate = &                ! [mol-bruc/m3/sec]
        max(this%brucitei,this%bruciteh) * &            ! [mol-bruc/kg-MgO/sec]
        D_m                                             ! [kg-MgO/m3]
  waste_panel%humid_brucite_rate = &                    ! [mol-bruc/m3/sec]
        this%bruciteh * &                               ! [mol-bruc/kg-MgO/sec]
        D_m                                             ! [kg-MgO/m3]
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  
end subroutine PMWSSProcessAfterRead

! *************************************************************************** !

subroutine PMWSSSetup(this)
  !
  ! Associates the waste panels to their regions and sets the waste panel id.
  ! Creates an MPI group/communicator for processes that own a waste panel.
  ! Throws out waste panels on processes that do not own the waste panel region.
  !
  ! Author: Jenn Frederick
  ! Date: 2/06/2017
  !
  
  use Option_module

  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
! -----------------------------------

! LOCAL VARIABLES:
! ================
! option: pointer to option object  
! cur_waste_panel: pointer to current waste panel object
! prev_waste_panel: pointer to previous waste panel object
! next_waste_panel: pointer to next waste panel object
! waste_panel_id: [-] waste panel id number
! i, j: [-] looping index integers
! local: Boolean that indicates whether a waste panel object is local to
!    the current process
! ierr: PETSc error integer
! newcomm_size: [-] new MPI communicator size
! ranks(:): [-] pointer to array of size(ranks) used to find local waste 
!    panel objects
! -----------------------------------------------------
  type(option_type), pointer :: option
  type(srcsink_panel_type), pointer :: cur_waste_panel
  type(srcsink_panel_type), pointer :: prev_waste_panel
  type(srcsink_panel_type), pointer :: next_waste_panel
  PetscInt :: waste_panel_id
  PetscInt :: i, j
  PetscBool :: local
  PetscErrorCode :: ierr
  PetscMPIInt :: newcomm_size
  PetscInt, pointer :: ranks(:)
! -----------------------------------------------------
  
  option => this%realization%option
  
  ! point the waste panel region to the desired region 
  call PMWSSAssociateRegion(this,this%realization%patch%region_list)
  
  allocate(ranks(option%mycommsize))
  
  waste_panel_id = 0
  nullify(prev_waste_panel)
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    waste_panel_id = waste_panel_id + 1
    local = PETSC_FALSE
    if (associated(cur_waste_panel%region)) then
      if (cur_waste_panel%region%num_cells > 0) then
          local = PETSC_TRUE
      endif
    endif
    ranks(:) = 0
    newcomm_size = 0
    if (local) then
      cur_waste_panel%id = waste_panel_id
      ranks(option%myrank+1) = 1
    else
      cur_waste_panel%id = 0
      ranks(option%myrank+1) = 0
    endif
    ! count the number of processes that own the waste panel
    call MPI_Allreduce(MPI_IN_PLACE,ranks,option%mycommsize,MPI_INTEGER, &
                       MPI_SUM,option%mycomm,ierr)
    newcomm_size = sum(ranks)
    allocate(cur_waste_panel%rank_list(newcomm_size))
    j = 0
    do i = 1,option%mycommsize
      if (ranks(i) == 1) then
        j = j + 1
        cur_waste_panel%rank_list(j) = (i - 1)
      endif
    enddo
    ! create an MPI group and communicator for each waste panel
    call MPI_Group_incl(option%mygroup,newcomm_size,cur_waste_panel%rank_list, &
                        cur_waste_panel%myMPIgroup,ierr)
    call MPI_Comm_create(option%mycomm,cur_waste_panel%myMPIgroup, &
                         cur_waste_panel%myMPIcomm,ierr)
    if (local) then
      call PMWSSSetRegionScaling(this,cur_waste_panel)
      call PMWSSProcessAfterRead(this,cur_waste_panel)
      call PMWSSInventoryAllocate(cur_waste_panel%inventory, &
                             cur_waste_panel%region%num_cells, &
                             cur_waste_panel%volume)
      ! allocate all of the rate arrays here
      allocate(cur_waste_panel%gas_generation_rate( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%gas_generation_rate(:) = 0.d0
      allocate(cur_waste_panel%brine_generation_rate( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%brine_generation_rate(:) = 0.d0
      
      allocate(cur_waste_panel%rxnrate_Fe_corrosion_inund( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_Fe_corrosion_inund(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_Fe_corrosion_humid( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_Fe_corrosion_humid(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_cell_biodeg_inund( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_cell_biodeg_inund(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_cell_biodeg_humid( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_cell_biodeg_humid(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_MgO_hyd_inund( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_MgO_hyd_inund(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_MgO_hyd_humid( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_MgO_hyd_humid(:) = 0.d0
      
      allocate(cur_waste_panel%rxnrate_Fe_corrosion( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_Fe_corrosion(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_cell_biodeg( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_cell_biodeg(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_FeOH2_sulf( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_FeOH2_sulf(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_Fe_sulf( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_Fe_sulf(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_MgO_hyd( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_MgO_hyd(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_MgOH2_carb( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_MgOH2_carb(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_MgO_carb( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_MgO_carb(:) = 0.d0
      allocate(cur_waste_panel%rxnrate_hydromag_conv( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%rxnrate_hydromag_conv(:) = 0.d0
      
      allocate(cur_waste_panel%solids_production( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%solids_production(:) = 0.d0
      prev_waste_panel => cur_waste_panel
      cur_waste_panel => cur_waste_panel%next
    else 
      ! remove waste panel because it is not local
      next_waste_panel => cur_waste_panel%next
      if (associated(prev_waste_panel)) then
        prev_waste_panel%next => next_waste_panel
      else
        this%waste_panel_list => next_waste_panel
      endif
      call PMWSSDestroyPanel(cur_waste_panel)
      cur_waste_panel => next_waste_panel
    endif
  enddo
  
  deallocate(ranks)
  
end subroutine PMWSSSetup

! ************************************************************************** !

subroutine PMWSSInventoryAllocate(inventory,num_cells,volume)
  ! 
  ! Allocates the size of the chemical species arrays within the inventory to 
  ! the number of cells in the waste panel region, and assigns the initial
  ! values.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/23/2017
  !
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! inventory (input/output): inventory object
! num_cells (input): [-] number of cells in a waste panel region
! volume (input): [m3] volume of a waste panel region
! ---------------------------------
  type(inventory_type) :: inventory
  PetscInt :: num_cells
  PetscReal :: volume
! ---------------------------------
  
  !----species-with-initial-inventory-values------------------
  call PMWSSChemSpeciesAllocate(num_cells,inventory%Fe_s, &
                                inventory%Fe_in_panel,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%MgO_s, &
                                inventory%MgO_in_panel,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%BioDegs_s, &
                                inventory%Biodegs_in_panel,volume)
  !----species-without-initial-inventory-values---------------
  !----assign-0.d0-kg-as-initial-value------------------------
  call PMWSSChemSpeciesAllocate(num_cells,inventory%FeOH2_s,0.d0,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%FeS_s,0.d0,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%MgOH2_s,0.d0,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%Mg5CO34OH24H2_s,0.d0,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%MgCO3_s,0.d0,volume)
  
end subroutine PMWSSInventoryAllocate

! ************************************************************************** !

subroutine PMWSSChemSpeciesAllocate(num_cells,chem_species,initial_mass,volume)
  ! 
  ! Allocates the size of the chemical species arrays within the inventory to 
  ! the number of cells in the waste panel region.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/23/2017
  !
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! num_cells (input): [-] number of cells in a waste panel region
! chem_species (input/output): chemical species object
! initial_mass (input): [kg] initial mass of chemical species in waste panel
! volume (input): [m3] volume of waste panel object
! ---------------------------------------
  PetscInt :: num_cells
  type(chem_species_type) :: chem_species
  PetscReal :: initial_mass 
  PetscReal :: volume 
! ---------------------------------------
  
  allocate(chem_species%initial_conc_mol(num_cells))
  allocate(chem_species%current_conc_mol(num_cells))
  allocate(chem_species%current_conc_kg(num_cells))
  allocate(chem_species%inst_rate(num_cells))
  
  !----------------------------------------------------------!-[units]--------
  chem_species%current_conc_kg(:) = initial_mass/volume      ! [kg/m3]
  chem_species%initial_conc_mol(:) = initial_mass/volume/ &  ! [kg/m3]
                                     chem_species%molar_mass ! [kg/mol]
  chem_species%current_conc_mol(:) =  &
                            chem_species%initial_conc_mol(:) ! [mol/m3]
  chem_species%inst_rate(:) = 0.d0                           ! [mol/m3/sec]        
  chem_species%tot_mass_in_panel = initial_mass              ! [kg]
  
end subroutine PMWSSChemSpeciesAllocate

! ************************************************************************** !

subroutine PMWSSInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/14/2017
  !
  
#include "petsc/finclude/petscis.h"
  use petscis
  use Data_Mediator_Vec_class
  use Realization_Base_class
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/ouput): wipp-srcsink process model object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
! -----------------------------------
  
! LOCAL VARIABLES:
! ================
! is: PETSc index set object
! i, j, k: [-] looping index integers
! ierr: PETSc error integer
! size_of_vec: [-] size of array
! dofs_in_residual(:): [-] degrees of freedom in residual array
! cur_waste_panel: pointer to curent waste panel object
! ----------------------------------------------------
  IS :: is
  PetscInt :: i, j, k
  PetscErrorCode :: ierr
  PetscInt :: size_of_vec
  PetscInt, allocatable :: dofs_in_residual(:)
  type(srcsink_panel_type), pointer :: cur_waste_panel
! ----------------------------------------------------
  
  ierr = 0
  
  ! set up mass transfer vector
  call RealizCreateFlowMassTransferVec(this%realization)
  this%data_mediator => DataMediatorVecCreate()
  call this%data_mediator%AddToList(this%realization%flow_data_mediator_list)
  
  ! create a Vec sized by # waste panels * # waste panel cells in region *
  ! # src/sinks (water [mol/s], gas [mol/s], energy [MJ/s])
  cur_waste_panel => this%waste_panel_list
  size_of_vec = 0
  do
    if (.not.associated(cur_waste_panel)) exit
    size_of_vec = size_of_vec + (cur_waste_panel%region%num_cells * &
                                 this%option%nflowdof)
    cur_waste_panel => cur_waste_panel%next
  enddo
  call VecCreateSeq(PETSC_COMM_SELF,size_of_vec, &
                    this%data_mediator%vec,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(this%data_mediator%vec,ierr);CHKERRQ(ierr)
  
  if (size_of_vec > 0) then
    allocate(dofs_in_residual(size_of_vec))
    dofs_in_residual = 0
    i = 0
    cur_waste_panel => this%waste_panel_list
    do
      if (.not.associated(cur_waste_panel)) exit
      do k = 1,cur_waste_panel%region%num_cells
        do j = 1,this%option%nflowdof
          i = i + 1
          dofs_in_residual(i) = &
            (cur_waste_panel%region%cell_ids(k)-1)*this%option%nflowdof + j
        enddo
      enddo
      cur_waste_panel => cur_waste_panel%next
    enddo
    ! zero-based indexing
    dofs_in_residual(:) = dofs_in_residual(:) - 1
    ! index to global petsc ordering
    dofs_in_residual(:) = dofs_in_residual(:) + &
               this%realization%patch%grid%global_offset * this%option%nflowdof
  endif
  ! create the index set (IS)
  call ISCreateGeneral(this%option%mycomm,size_of_vec, &
                       dofs_in_residual, &
                       PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)
  if (allocated(dofs_in_residual)) deallocate(dofs_in_residual)
  ! load the data mediator vec scatter context with the IS
  call VecScatterCreate(this%data_mediator%vec,PETSC_NULL_VEC, &
                        this%realization%field%flow_r,is, &
                        this%data_mediator%scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)

  ! write header in the *.pnl files
  call PMWSSOutputHeader(this)
  call PMWSSOutput(this)

  ! solve for initial process model state
  call PMWSSSolve(this,0.d0,ierr)
  
end subroutine PMWSSInitializeRun

! *************************************************************************** !

subroutine PMWSSInitializeTimestep(this)
  ! 
  ! Initializes the process model to take a time step in the simulation:
  ! Updates the waste panel inventory.
  ! Writes data to output file from end of previous time step.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  !
  
  use Option_module
  
  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
! -----------------------------------
  
! LOCAL VARIABLES:
! ================
! cur_waste_panel: pointer to current waste panel object
! option: pointer to option object
! dt: [sec] flow time step value
! ----------------------------------------------------
  type(srcsink_panel_type), pointer :: cur_waste_panel
  type(option_type), pointer :: option
  PetscReal :: dt
! ----------------------------------------------------

  option => this%option
  dt = option%flow_dt  ! [sec]
  
  if (option%print_screen_flag) then
    write(*,'(/,2("=")," WIPP SRC/SINK PANEL MODEL ",51("="))')
  endif
  
  ! update the waste panel inventory
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    call PMWSSUpdateInventory(cur_waste_panel,option)
    ! set current dt after using the previous in the update above
    cur_waste_panel%prev_dt = dt
    cur_waste_panel => cur_waste_panel%next
  enddo
  
  ! write data to *.pnl output files from previous time step
  if (option%time >= this%output_start_time) then
    call PMWSSOutput(this)
  endif
  
end subroutine PMWSSInitializeTimestep

! *************************************************************************** !

subroutine PMWSSUpdateInventory(waste_panel,option)
  !
  ! Updates the waste panel tracked species inventory concentrations.
  !
  ! Author: Jenn Frederick
  ! Date: 02/10/2017
  !
  
  use Option_module
 
  implicit none
  
! INPUT ARGUMENTS:
! ================
! waste_panel (input/output): waste panel object
! option (input/output): pointer to option object
! ---------------------------------------
  type(srcsink_panel_type) :: waste_panel
  type(option_type), pointer :: option
! ---------------------------------------

! LOCAL VARIABLES:
! ================
! dt: [sec] flow time step value (flow_dt)
! ---------------
  PetscReal :: dt
! ---------------

  dt =  waste_panel%prev_dt
 
  call PMWSSUpdateChemSpecies(waste_panel%inventory%Fe_s,waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%FeOH2_s,waste_panel,dt, &
                              option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%BioDegs_s,waste_panel,dt, &
                              option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%FeS_s,waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%MgO_s,waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%MgOH2_s,waste_panel,dt, &
                              option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%Mg5CO34OH24H2_s, &
                              waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%MgCO3_s,waste_panel,dt, &
                              option)
               
  ! solids production calculation
  ! BRAGFLO User's Manual Section 4.13.7, eq. 159, eq. 160
  ! this sum does not include salt
  waste_panel%solids_production = 0.d0
  waste_panel%solids_production = &
    ((waste_panel%inventory%Fe_s%current_conc_kg - &
      (waste_panel%inventory%Fe_s%initial_conc_mol* &
       waste_panel%inventory%Fe_s%molar_mass))/DN_FE) + &
    ((waste_panel%inventory%FeOH2_s%current_conc_kg - &
      (waste_panel%inventory%FeOH2_s%initial_conc_mol* &
       waste_panel%inventory%FeOH2_s%molar_mass))/DN_FEOH2) + &
    ((waste_panel%inventory%FeS_s%current_conc_kg - &
      (waste_panel%inventory%FeS_s%initial_conc_mol* &
       waste_panel%inventory%FeS_s%molar_mass))/DN_FES) + &
    ((waste_panel%inventory%BioDegs_s%current_conc_kg - &
      (waste_panel%inventory%BioDegs_s%initial_conc_mol* &
       waste_panel%inventory%BioDegs_s%molar_mass))/DN_CELL) + &
    ((waste_panel%inventory%MgO_s%current_conc_kg - &
      (waste_panel%inventory%MgO_s%initial_conc_mol* &
       waste_panel%inventory%MgO_s%molar_mass))/DN_MGO) + &
    ((waste_panel%inventory%MgOH2_s%current_conc_kg - &
      (waste_panel%inventory%MgOH2_s%initial_conc_mol* &
       waste_panel%inventory%MgOH2_s%molar_mass))/DN_MGOH2) + &
    ((waste_panel%inventory%Mg5CO34OH24H2_s%current_conc_kg - &
      (waste_panel%inventory%Mg5CO34OH24H2_s%initial_conc_mol* &
       waste_panel%inventory%Mg5CO34OH24H2_s%molar_mass))/DN_HYDRO) + &
    ((waste_panel%inventory%MgCO3_s%current_conc_kg - &
      (waste_panel%inventory%MgCO3_s%initial_conc_mol* &
       waste_panel%inventory%MgCO3_s%molar_mass))/DN_MGCO3)
                                      
 end subroutine PMWSSUpdateInventory
 
! *************************************************************************** !

subroutine PMWSSUpdateChemSpecies(chem_species,waste_panel,dt,option)
  !
  ! Updates the waste panel tracked species inventory concentrations.
  !
  ! Author: Jenn Frederick
  ! Date: 2/23/2017
  !
  
  use Option_module
  use Utility_module
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! chem_species (input/output): chemical species object
! waste_panel (input): waste panel object
! dt (input): [sec] flow time step value
! option (input/output): pointer to option object
! ---------------------------------------
  type(chem_species_type) :: chem_species
  type(srcsink_panel_type) :: waste_panel
  PetscReal :: dt
  type(option_type), pointer :: option
! ---------------------------------------
  
! LOCAL VARIABLES:
! ================
! k: [-] looping index integer
! num_cells: [-] number of cells in a waste panel region
! local_conc_kg: [kg/m3] local concentration of a chemical species
! global_conc_kg: [kg/m3] global concentration of a chemical species
! ---------------------------
  PetscInt :: k
  PetscInt :: num_cells
  PetscReal :: local_conc_kg 
  PetscReal :: global_conc_kg
! ---------------------------
  
  num_cells = waste_panel%region%num_cells
  local_conc_kg = 0.d0  
  do k = 1,num_cells
    !--[mol/m3]-----------------------------------------!-units---------------
    chem_species%current_conc_mol(k) = &                ! [mol/m3]
                 chem_species%current_conc_mol(k) + &   ! [mol/m3]
                 chem_species%inst_rate(k) * &          ! [mol/m3/sec]
                 dt                                     ! [sec]
    chem_species%current_conc_mol = max(0.d0,chem_species%current_conc_mol)
    
    !--[kg/m3]------------------------------------------!-units---------------
    chem_species%current_conc_kg(k) = &                 ! [kg/m3]
                 chem_species%current_conc_mol(k) * &   ! [mol/m3]
                 chem_species%molar_mass                ! [kg/mol]
    chem_species%current_conc_kg = max(0.d0,chem_species%current_conc_kg)
    
    !--[kg/m3]------------------------------------------!-units---------------             
    local_conc_kg = local_conc_kg + &                   ! [kg/m3]
                    (chem_species%current_conc_kg(k) * &! [kg/m3]
                     waste_panel%scaling_factor(k))     ! [-]
  enddo
  call CalcParallelSUM(option,waste_panel%rank_list,local_conc_kg, &
                       global_conc_kg)
  !--[kg]-----------------------------------------------!-units---------------
  chem_species%tot_mass_in_panel = global_conc_kg * &   ! [kg/m3]
                                   waste_panel%volume   ! [m3]
                                   
end subroutine PMWSSUpdateChemSpecies

! *************************************************************************** !

 subroutine PMWSSSolve(this,time,ierr)
  ! 
  ! Calculates reaction rates, gas generation rate, and brine generation rate.
  ! Sets the fluid and energy source terms.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  !
  
  use Option_module
  use Grid_module
  use General_Aux_module
  use WIPP_Flow_Aux_module
  use Material_Aux_class
  use Global_Aux_module
  use EOS_Gas_module
  use EOS_Water_module
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! time (input): [sec] simulation time
! ierr (input/output): PETSc error integer
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
! -----------------------------------
  
! LOCAL VARIABLES:
! ================
! option: pointer to option object
! grid: pointer to grid object
! gen_auxvar(:,:): pointer to general mode auxvar object, which stores the 
!    liquid saturation, temperature, and pressure at each grid cell, and
!    indexed by the ghosted grid cell id
! global_auxvar(:): pointer to global auxvar object, which stores the phase
!    state of the system (gas, liquid, or two-phase state), and indexed by
!    the ghosted grid cell id
! material_auxvars(:): pointer to material auxvars object, which stores the
!    grid cell volume, and indexed by the ghosted grid cell id
! vec_p(:): pointer to the data mediator array, which stores the liquid
!    source term [kmol/sec], gas source term [kmol/sec], and energy source
!    term [MJ/sec], and is indexed by j
! cwp: pointer to current waste panel object
! i, j: [-] looping index integers
! local_id, ghosted_id: [-] local and ghosted grid cell id number
! num_cells: number of grid cells in a waste panel region
! SOCEXP: exponent value in effective brine saturation equation
! water_saturation: [-] liquid saturation from simulation
! s_eff: [-] effective brine saturation due to capillary action
!    in the waste materials
! -----------------------------------------------------------
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(general_auxvar_type), pointer :: general_auxvar(:,:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvar(:,:)
  type(global_auxvar_type), pointer :: global_auxvar(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, pointer :: vec_p(:)
  type(srcsink_panel_type), pointer :: cwp
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id
  PetscInt :: num_cells
  PetscReal :: SOCEXP
  PetscReal :: dt
  
  ! brine/gas generation variables
  PetscReal :: water_saturation
  PetscReal :: s_eff
  PetscReal :: sg_eff
  ! enthalpy calculation variables
  PetscReal :: temperature
  PetscReal :: pressure_liq
  PetscReal :: pressure_gas
  PetscReal :: H_liq
  PetscReal :: H_gas
  PetscReal :: U_gas
  PetscReal :: gas_energy
  PetscReal :: brine_energy
! -----------------------------------------------------------
  
  option => this%realization%option
  grid => this%realization%patch%grid
  nullify(general_auxvar)
  nullify(wippflo_auxvar)
  if (associated(this%realization%patch%aux%General)) then
    general_auxvar => this%realization%patch%aux%General%auxvars
  else
    wippflo_auxvar => this%realization%patch%aux%WIPPFlo%auxvars
  endif
  material_auxvars => this%realization%patch%aux%Material%auxvars
  global_auxvar => this%realization%patch%aux%Global%auxvars
  
  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
  dt = option%flow_dt  ! [sec]
  
  j = 0  ! (j indexes the data mediator vec)
  cwp => this%waste_panel_list
  do
    if (.not.associated(cwp)) exit
    num_cells = cwp%region%num_cells
    !-----zero-out-reaction-rates---------------------------------------------
    water_saturation = 0.d0
    s_eff = 0.d0
    sg_eff = 0.d0
    
    cwp%rxnrate_Fe_corrosion_inund(:) = 0.d0
    cwp%rxnrate_Fe_corrosion_humid(:) = 0.d0
    cwp%rxnrate_cell_biodeg_inund(:) = 0.d0
    cwp%rxnrate_cell_biodeg_humid(:) = 0.d0
    cwp%rxnrate_MgO_hyd_inund(:) = 0.d0
    cwp%rxnrate_MgO_hyd_humid(:) = 0.d0
    cwp%rxnrate_Fe_corrosion(:) = 0.d0
    cwp%rxnrate_cell_biodeg(:) = 0.d0
    cwp%rxnrate_FeOH2_sulf(:) = 0.d0
    cwp%rxnrate_Fe_sulf(:) = 0.d0
    cwp%rxnrate_MgO_hyd(:) = 0.d0
    cwp%rxnrate_MgOH2_carb(:) = 0.d0
    cwp%rxnrate_MgO_carb(:) = 0.d0
    cwp%rxnrate_hydromag_conv(:) = 0.d0
    

    do i = 1,num_cells
      local_id = cwp%region%cell_ids(i)
      ghosted_id = grid%nL2G(local_id)
    !-----effective-brine-saturation------------------------------------------
    !-----(see equation PA.99, section PA-4.2.6)------------------------------
      if (associated(general_auxvar)) then
        water_saturation = &
          general_auxvar(ZERO_INTEGER,ghosted_id)%sat(option%liquid_phase)
      else
        water_saturation = &
          wippflo_auxvar(ZERO_INTEGER,ghosted_id)%sat(option%liquid_phase)
      endif
      
      SOCEXP = 200.d0*(max((water_saturation-this%smin),0.d0))**2.d0
      s_eff = water_saturation-this%smin+(this%satwick * &
                                          (1.d0-exp(this%alpharxn*SOCEXP)))
      if (water_saturation <= this%smin) s_eff = 0.d0
      if (water_saturation > (1.d0-this%satwick+this%smin)) s_eff = 1.d0
      
      s_eff = min(s_eff,1.d0)
      s_eff = max(s_eff,0.d0)
      ! sg_eff is set to zero if sw_eff is also zero
      sg_eff = 0.d0
      if (s_eff > 1.d-16) then
        sg_eff = (1.d0-s_eff)
      endif
      
    ! note, bragflo applies a separate smoothing on the humid rates using s_eff 
    !
    ! bragflo calculates and reports inundated and humid rates separately
    ! inundated_rate_cons*s_eff and humid_rate_constant*sg_eff
    ! and sequentially - humid proceeds only if there is remaining solid after inundated
    
    
    !-----anoxic-iron-corrosion-[mol-Fe/m3/sec]-------------------------------
    !-----(see equation PA.67, PA.77, section PA-4.2.5)-----------------------
    
      ! CORSAT
      cwp%rxnrate_Fe_corrosion_inund(i) = cwp%inundated_corrosion_rate*s_eff
      call PMWSSSmoothRxnrate(cwp%rxnrate_Fe_corrosion_inund(i),i, &
                              cwp%inventory%Fe_s,this%alpharxn) 
      call PMWSSTaperRxnrate(cwp%rxnrate_Fe_corrosion_inund(i),i, &
                              cwp%inventory%Fe_s,1.d0,dt)
      
      ! CORHUM
      cwp%rxnrate_Fe_corrosion_humid(i) = cwp%humid_corrosion_rate*sg_eff
      call PMWSSSmoothRxnrate(cwp%rxnrate_Fe_corrosion_humid(i),i, &
                              cwp%inventory%Fe_s,this%alpharxn) 
      call PMWSSTaperRxnrate(cwp%rxnrate_Fe_corrosion_humid(i),i, &
                              cwp%inventory%Fe_s,1.d0,dt)
      ! the humid rate is also smoothed based on sg_eff; neglected here
      
      ! total corrosion
      cwp%rxnrate_Fe_corrosion(i) = & 
          cwp%rxnrate_Fe_corrosion_inund(i) + cwp%rxnrate_Fe_corrosion_humid(i)
      
      
    !-----biodegradation-[mol-cell/m3/sec]------------------------------------
    !-----(see equation PA.69, PA.82, PA.83, section PA-4.2.5)----------------
      
      ! BIOSAT
      cwp%rxnrate_cell_biodeg_inund(i) = cwp%inundated_biodeg_rate*s_eff
      call PMWSSSmoothRxnrate(cwp%rxnrate_cell_biodeg_inund(i),i, &
                              cwp%inventory%BioDegs_s,this%alpharxn) 
      call PMWSSTaperRxnrate(cwp%rxnrate_cell_biodeg_inund(i),i, &
                              cwp%inventory%BioDegs_s,1.d0,dt)
      
      ! BIOHUM
      cwp%rxnrate_cell_biodeg_humid(i) = cwp%humid_biodeg_rate*sg_eff
      call PMWSSSmoothRxnrate(cwp%rxnrate_cell_biodeg_humid(i),i, &
                              cwp%inventory%BioDegs_s,this%alpharxn) 
      call PMWSSTaperRxnrate(cwp%rxnrate_cell_biodeg_humid(i),i, &
                              cwp%inventory%BioDegs_s,1.d0,dt)
      ! the humid rate is also smoothed based on sg_eff; neglected here
      
      ! total microbial gas generation
      cwp%rxnrate_cell_biodeg(i) = & 
          cwp%rxnrate_cell_biodeg_inund(i) + cwp%rxnrate_cell_biodeg_humid(i)
      
      
    !-----iron-sulfidation-[mol-H2S/m3/sec]-----------------------------------
    !-----(see equation PA.68, PA.89, PA.90, section PA-4.2.5)----------------
      
      ! BIOFES
      ! FeOH2 sulfidation is assumed to kinetically dominate Fe sulfidation.
      ! The H2S generation rate is proportioned between FeOH and Fe.
      ! FeOH2 is sulifidized first, then Fe is sulifidized with remaining H2S
      cwp%rxnrate_FeOH2_sulf(i) = cwp%rxnrate_cell_biodeg(i)*cwp%RXH2S_factor
      ! bragflo uses Fe as initial concentration, not FeOH2, so this is correct
      
      ! for smoothing, the initial and current species are different
      call PMWSSSmoothRxnrate2(cwp%rxnrate_FeOH2_sulf(i), &
                               cwp%inventory%FeOH2_s%current_conc_mol(i), & 
                               cwp%inventory%Fe_s%initial_conc_mol(i), & 
                               this%alpharxn)
      
      ! here tapering actually truncates the rate, so this call is important
      call PMWSSTaperRxnrate(cwp%rxnrate_FeOH2_sulf(i),i,cwp%inventory%FeOH2_s,1.d0,dt)
      
      ! FeS rate is caculated based on the left over H2S rate
      cwp%rxnrate_Fe_sulf(i) = cwp%rxnrate_cell_biodeg(i)*cwp%RXH2S_factor - &
                              cwp%rxnrate_FeOH2_sulf(i)
      call PMWSSTaperRxnrate(cwp%rxnrate_Fe_sulf(i),i,cwp%inventory%Fe_s,1.d0,dt)
      
      
    !-----MgO-hydration-[mol-MgO/m3/sec]--------------------------------------
    !-----(see equation PA.73, PA.94, section PA-4.2.5)-----------------------
      
      ! CORMGO
      cwp%rxnrate_MgO_hyd_inund(i) = cwp%inundated_brucite_rate*s_eff
      cwp%rxnrate_MgO_hyd_humid(i) = cwp%humid_brucite_rate*sg_eff
      ! the humid rate is also smoothed based on sg_eff; neglected here
      
      ! unlike for corrorosion and biodegradation, bragflo doesn't
      ! separately store and report inundated and humid rates for MgO hydration
      cwp%rxnrate_MgO_hyd(i) = & 
          cwp%rxnrate_MgO_hyd_inund(i) + cwp%rxnrate_MgO_hyd_humid(i)
      
      call PMWSSSmoothRxnrate(cwp%rxnrate_MgO_hyd(i),i, &
                              cwp%inventory%MgO_s,this%alpharxn) 
      call PMWSSTaperRxnrate(cwp%rxnrate_MgO_hyd(i),i, &
                              cwp%inventory%MgO_s,1.d0,dt)
      
    !-----hydromagnesite-[mol-hydromagnesite/m3-bulk/sec]---------------------
    !-----(see equation PA.74, PA.96, section PA-4.2.5)-----------------------
      
      ! BIOMGO
      cwp%rxnrate_MgOH2_carb(i) = cwp%rxnrate_cell_biodeg(i)*cwp%RXCO2_factor
      
      ! for smoothing, the initial and current species are different
      call PMWSSSmoothRxnrate2(cwp%rxnrate_MgOH2_carb(i), &
                               cwp%inventory%MgOH2_s%current_conc_mol(i), & 
                               cwp%inventory%MgO_s%initial_conc_mol(i), & 
                               this%alpharxn)
      
      ! here tapering actually truncates the rate, so this call is important
      call PMWSSTaperRxnrate(cwp%rxnrate_MgOH2_carb(i),i,cwp%inventory%MgOH2_s,1.d0,dt)
      
      ! the remaining CO2 reacts with MgO
      cwp%rxnrate_MgO_carb(i) = cwp%rxnrate_cell_biodeg(i)*cwp%RXCO2_factor - &
                                cwp%rxnrate_MgOH2_carb(i)
      call PMWSSTaperRxnrate(cwp%rxnrate_MgO_carb(i),i,cwp%inventory%MgO_s,1.d0,dt)
      
    !-----hydromagnesite-conversion-[mol-hydromagnesite/m3-bulk/sec]----------
    !-----(see equation PA.74, PA.97, section PA-4.2.5)-----------------------
      
      ! HYDROCONV
      ! this is a first-order reaction, not zero-order; also rate is on /kg basis
      cwp%rxnrate_hydromag_conv(i) = this%hymagcon_rate* &
                                cwp%inventory%Mg5CO34OH24H2_s%current_conc_kg(i)
      ! for smoothing, the initial and current species are different
      call PMWSSSmoothRxnrate2(cwp%rxnrate_hydromag_conv(i), &
                               cwp%inventory%Mg5CO34OH24H2_s%current_conc_mol(i), & 
                               cwp%inventory%MgO_s%initial_conc_mol(i), & 
                               this%alpharxn)
      call PMWSSTaperRxnrate(cwp%rxnrate_hydromag_conv(i),i, &
                             cwp%inventory%Mg5CO34OH24H2_s,1.d0,dt)
      
    !-----tracked-species-[mol-species/m3-bulk/sec]---------------------------
    !-----note:column-id-is-shifted-by-+1-since-a-0-index-not-possible--------
      
      cwp%inventory%FeOH2_s%inst_rate(i) = &
                    this%stoic_mat(1,6)*cwp%rxnrate_Fe_corrosion(i) + &
                    this%stoic_mat(3,6)*cwp%rxnrate_FeOH2_sulf(i)
      cwp%inventory%Fe_s%inst_rate(i) = &
                    this%stoic_mat(1,4)*cwp%rxnrate_Fe_corrosion(i) + &
                    this%stoic_mat(4,4)*cwp%rxnrate_Fe_sulf(i) 
      cwp%inventory%FeS_s%inst_rate(i) = &
                    this%stoic_mat(4,7)*cwp%rxnrate_Fe_sulf(i) + &
                    this%stoic_mat(3,7)*cwp%rxnrate_FeOH2_sulf(i)
      cwp%inventory%BioDegs_s%inst_rate(i) = &
                    this%stoic_mat(2,5)*cwp%rxnrate_cell_biodeg(i)
      cwp%inventory%MgO_s%inst_rate(i) = &
                    this%stoic_mat(5,8)*cwp%rxnrate_MgO_hyd(i) + &
                    this%stoic_mat(7,8)*cwp%rxnrate_MgO_carb(i)
      cwp%inventory%MgOH2_s%inst_rate(i) = &
                    this%stoic_mat(5,9)*cwp%rxnrate_MgO_hyd(i) + & 
                    this%stoic_mat(6,9)*cwp%rxnrate_MgOH2_carb(i) + &
                    this%stoic_mat(8,9)*cwp%rxnrate_hydromag_conv(i)
      cwp%inventory%Mg5CO34OH24H2_s%inst_rate(i) = &
                    this%stoic_mat(6,1)*cwp%rxnrate_MgOH2_carb(i) + &
                    this%stoic_mat(8,1)*cwp%rxnrate_hydromag_conv(i) 
      cwp%inventory%MgCO3_s%inst_rate(i) = &
                    this%stoic_mat(7,10)*cwp%rxnrate_MgO_carb(i) + &
                    this%stoic_mat(8,10)*cwp%rxnrate_hydromag_conv(i)
      
    !-----gas-generation-[mol-H2/m3-bulk/sec]---------------------------------
    !-----(see equations PA.67-69, PA.77, PA.82-83, PA.86, PA.89, sec PA-4.2.5)
    
      cwp%gas_generation_rate(i) = &
                    this%stoic_mat(1,2)*cwp%rxnrate_Fe_corrosion(i) + &
                    this%stoic_mat(3,2)*cwp%rxnrate_FeOH2_sulf(i) + &
                    this%stoic_mat(4,2)*cwp%rxnrate_Fe_sulf(i) + &
                    this%stoic_mat(2,2)*cwp%rxnrate_cell_biodeg(i) + & ! zero
                    cwp%RXH2_factor*cwp%rxnrate_cell_biodeg(i)         ! H2
      
    !-----brine-generation-[mol-H2O/m3-bulk/sec]------------------------------
    !-----(see equations PA.77, PA.78, PA.82, PA.83, PA.90, PA.96, PA.97,)----
    !-----(section PA-4.2.5)--------------------------------------------------
      
      cwp%brine_generation_rate(i) = &
                    this%stoic_mat(1,3)*cwp%rxnrate_Fe_corrosion(i) + &
                    this%stoic_mat(3,3)*cwp%rxnrate_FeOH2_sulf(i) + &
                    this%stoic_mat(5,3)*cwp%rxnrate_MgO_hyd(i) + &
                    this%stoic_mat(8,3)*cwp%rxnrate_hydromag_conv(i) + &
                    this%stoic_mat(2,3)*cwp%rxnrate_cell_biodeg(i) + & ! STCO_22=SMIC_H20
                    cwp%RXH2O_factor*cwp%rxnrate_cell_biodeg(i)        ! zero
      
    !------source-term-calculation--------------------------------------------
      j = j + 1
      !---liquid-source-term-[kmol/sec]------------------------!-[units]------
      vec_p(j) = cwp%brine_generation_rate(i) * &              ! [mol/m3/sec]
                 material_auxvars(ghosted_id)%volume / &       ! [m3-bulk]
                 1.d3                                          ! [mol -> kmol]
      print *, 'vec_p brine [kmol/sec]'
      print *, vec_p(j)
      print *, 'brine_generation_rate [mol/m3/sec]'
      print *, cwp%brine_generation_rate(i)
      j = j + 1
      !---gas-source-term-[kmol/sec]---------------------------!-[units]------
      vec_p(j) = cwp%gas_generation_rate(i) * &                ! [mol/m3/sec]
                 material_auxvars(ghosted_id)%volume / &       ! [m3-bulk]
                 1.d3                                          ! [mol -> kmol]
      print *, 'vec_p gas [kmol/sec]'
      print *, vec_p(j)
      print *, 'gas_generation_rate [mol/m3/sec]'
      print *, cwp%gas_generation_rate(i)
      if (associated(general_auxvar)) then
        j = j + 1
        !---energy-source-term-[MJ/sec];-H-from-EOS-[J/kmol]-------------------
        brine_energy = 0.d0
        gas_energy = 0.d0
        temperature = general_auxvar(ZERO_INTEGER,ghosted_id)%temp
        select case(global_auxvar(ghosted_id)%istate)
          case(GAS_STATE) !----------------------------------------------------
            pressure_gas = general_auxvar(ZERO_INTEGER,ghosted_id)% &
                           pres(option%gas_phase)
            call EOSGasEnergy(temperature,pressure_gas,H_gas,U_gas,ierr)
            gas_energy = & !---[MJ/sec]-----------------------!-[units]--------
                cwp%gas_generation_rate(i) * &                ! [mol/m3/sec]
                material_auxvars(ghosted_id)%volume * &       ! [m3-bulk] 
                H_gas * 1.d-3 * 1.d-6                         ! [MJ/mol]
          case(LIQUID_STATE) !-------------------------------------------------
            pressure_liq = general_auxvar(ZERO_INTEGER,ghosted_id)% &
                           pres(option%liquid_phase)
            call EOSWaterEnthalpy(temperature,pressure_liq,H_liq,ierr)
            brine_energy = & !---[MJ/sec]---------------------!-[units]--------
                cwp%brine_generation_rate(i) * &              ! [mol/m3/sec]
                material_auxvars(ghosted_id)%volume * &       ! [m3-bulk] 
                H_liq * 1.d-3 * 1.d-6                         ! [MJ/mol]
          case(TWO_PHASE_STATE) !----------------------------------------------
            pressure_liq = general_auxvar(ZERO_INTEGER,ghosted_id)% &
                           pres(option%liquid_phase)
            pressure_gas = general_auxvar(ZERO_INTEGER,ghosted_id)% &
                           pres(option%gas_phase)
            call EOSWaterEnthalpy(temperature,pressure_liq,H_liq,ierr)
            call EOSGasEnergy(temperature,pressure_gas,H_gas,U_gas,ierr)
            brine_energy = & !---[MJ/sec]---------------------!-[units]--------
                cwp%brine_generation_rate(i) * &              ! [mol/m3/sec]
                material_auxvars(ghosted_id)%volume * &       ! [m3-bulk] 
                H_liq * 1.d-3 * 1.d-6                         ! [MJ/mol]
            gas_energy = & !---[MJ/sec]-----------------------!-[units]--------
                cwp%gas_generation_rate(i) * &                ! [mol/m3/sec]
                material_auxvars(ghosted_id)%volume * &       ! [m3-bulk] 
                H_gas * 1.d-3 * 1.d-6                         ! [MJ/mol]
        end select
        vec_p(j) = brine_energy + gas_energy
      endif
    enddo
    !-------------------------------------------------------------------------
    cwp => cwp%next
  enddo
  
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine PMWSSSolve

! ************************************************************************** !

subroutine PMWSSSmoothRxnrate(rxnrate,cell_num,limiting_species,alpharxn)
  !
  ! Smooths the reaction rate near the point where the reaction runs out of a
  ! limiting relevant reactant/species. This implements Eq. 158 in the BRAGFLO
  ! User's Manual.
  !
  ! Author: Jennifer Frederick
  ! Date: 03/27/3017
  !
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! rxnrate (input/output): [mol-species/m3-bulk/sec] reaction rate constant 
! cell_num (input): [-] grid cell numbering in waste panel region
! limiting_species (input): chemical species object that is the limiting
!    species of the reaction with reaction rate constant "rxnrate"
! alpharxn (input): [-] BRAGFLO parameter ALPHARXN
! -------------------------------------------
  PetscReal :: rxnrate
  PetscInt :: cell_num
  type(chem_species_type) :: limiting_species
  PetscReal :: alpharxn
! -------------------------------------------
  
! LOCAL VARIABLES:
! ================
! conc_ratio: [-] concentratio ratio of current molar concentration to
!    initial molar concentration of a chemical species
! -----------------------
  PetscReal :: conc_ratio
! -----------------------
  if ( limiting_species%initial_conc_mol(cell_num) > 0.0d0) then
    conc_ratio = ( limiting_species%current_conc_mol(cell_num) / &
                   limiting_species%initial_conc_mol(cell_num) ) 
    conc_ratio = min(1.d0,conc_ratio)
    ! K_smoothed = K * (1.0 - exp(A*C/Ci)  BRAGFLO User's Manual Eq. 158
    ! However, the above equation is an error. The correct equation is below:
    rxnrate = rxnrate * (1.d0 - exp(alpharxn*conc_ratio))
  else
    rxnrate = 0.0d0
  end if
  
end subroutine PMWSSSmoothRxnrate


! ************************************************************************** !

subroutine PMWSSSmoothRxnrate2(rxnrate,current_conc,initial_conc,alpharxn)
  !
  ! Smooths the reaction rate near the point where the reaction runs out of a
  ! limiting relevant reactant/species. This implements Eq. 158 in the BRAGFLO
  ! User's Manual.
  !
  ! Author: Jennifer Frederick
  ! Date: 03/27/3017
  !
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! rxnrate (input/output): [mol-species/m3-bulk/sec] reaction rate constant 
! cell_num (input): [-] grid cell numbering in waste panel region
! limiting_species (input): chemical species object that is the limiting
!    species of the reaction with reaction rate constant "rxnrate"
! alpharxn (input): [-] BRAGFLO parameter ALPHARXN
! -------------------------------------------
  PetscReal :: rxnrate
  PetscReal :: current_conc
  PetscReal :: initial_conc
  PetscReal :: alpharxn
! -------------------------------------------
  
! LOCAL VARIABLES:
! ================
! conc_ratio: [-] concentratio ratio of current molar concentration to
!    initial molar concentration of a chemical species
! -----------------------
  PetscReal :: conc_ratio
! -----------------------
  
  if (initial_conc > 0.0d0) then
    conc_ratio = current_conc/initial_conc
    conc_ratio = min(1.d0,conc_ratio)
    ! K_smoothed = K * (1.0 - exp(A*C/Ci)  BRAGFLO User's Manual Eq. 158
    ! However, the above equation is an error. The correct equation is below:
    rxnrate = rxnrate * (1.d0 - exp(alpharxn*conc_ratio))
  else
    rxnrate = 0.0d0
  end if
  
end subroutine PMWSSSmoothRxnrate2


! ************************************************************************** !

subroutine PMWSSTaperRxnrate1(rxnrate,cell_num,limiting_species1,stocoef,dt)
  !
  ! Tapers the reaction rate if the reaction runs out of a single
  ! limiting relevant reactant/species. The limiting reactant/species is
  ! chosen according to the equations in the BRAGFLO User's Manual, 
  ! Section 14.13.
  !
  ! Author: Jennifer Frederick
  ! Date: 03/27/3017
  !
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! rxnrate (input/output): [mol-species/m3-bulk/sec] reaction rate constant
! cell_num (input): [-] grid cell numbering in waste panel region
! limiting_species1 (input): chemical species object that is the limiting
!    species of the reaction with reaction rate constant "rxnrate"
! stocoef: [mol/mol] stoichiometric coefficient for limiting reactant
! dt: timestep size [sec]
! --------------------------------------------
  PetscReal :: rxnrate
  PetscInt :: cell_num
  type(chem_species_type) :: limiting_species1
  PetscReal :: stocoef
  PetscReal :: dt
  
  PetscReal :: available_concentration
  PetscReal :: reacted_concentration
! --------------------------------------------
  available_concentration = limiting_species1%current_conc_mol(cell_num)
  
  if (available_concentration <= 0.d0) then
    rxnrate = 0.d0
  else
    reacted_concentration = dt*stocoef*rxnrate
    if (reacted_concentration >= available_concentration) then
      rxnrate = available_concentration/(stocoef*dt)
    endif
  endif
  
end subroutine PMWSSTaperRxnrate1

! ************************************************************************** !

subroutine PMWSSTaperRxnrate2(rxnrate,cell_num,limiting_species1, &
                              limiting_species2)
  !
  ! Tapers the reaction rate if the reaction runs out of one of two possible
  ! limiting relevant reactants/species. The limiting reactants/species are
  ! chosen according to the equations in the BRAGFLO User's Manual, 
  ! Section 14.13.
  !
  ! Author: Jennifer Frederick
  ! Date: 03/27/3017
  !
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! rxnrate (input/output): [mol-species/m3-bulk/sec] reaction rate constant
! cell_num (input): [-] grid cell numbering in waste panel region
! limiting_species1 (input): first possible chemical species object that 
!    is the limiting species of the reaction with reaction rate 
!    constant "rxnrate"
! limiting_species2 (input): second possible chemical species object that 
!    is the limiting species of the reaction with reaction rate  
!    constant "rxnrate"
! --------------------------------------------
  PetscReal :: rxnrate
  PetscInt :: cell_num
  type(chem_species_type) :: limiting_species1
  type(chem_species_type) :: limiting_species2
! --------------------------------------------
  
! LOCAL_VARIABLES:
! ================
! rxnrate1, rxnrate2: [mol-species/m3-bulk/sec] reaction rate constant
! ---------------------
  PetscReal :: rxnrate1
  PetscReal :: rxnrate2
! ---------------------
  
  if (limiting_species1%current_conc_mol(cell_num) <= 0.d0) then
    rxnrate1 = 0.d0
  else
    rxnrate1 = rxnrate
  endif
  if (limiting_species2%current_conc_mol(cell_num) <= 0.d0) then
    rxnrate2 = 0.d0
  else
    rxnrate2 = rxnrate
  endif
  rxnrate = min(rxnrate1,rxnrate2)
  
end subroutine PMWSSTaperRxnrate2

! ************************************************************************** !

subroutine PMWSSFinalizeTimestep(this)
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/22/2017
  !

  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
! -----------------------------------
  
end subroutine PMWSSFinalizeTimestep

! *************************************************************************** !
 subroutine PMWSSOutput(this)
  ! 
  ! Sets up output for the process model in the *.pnl files.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/27/2017
  !
  
  use Option_module
  use Output_Aux_module
  use Utility_module
  
  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
! -----------------------------------
  
! LOCAL VARIABLES:
! ================
! option: pointer to option object
! output_option: pointer to output option object
! cwp: pointer to current waste panel object
! filename: filename string
! fid: [-] file id integer
! i: [-] looping index integer
! local_variables: [units vary] an array of local variable values for the 
!    parallel averaging calculation
! global_variables: [units vary] an array of global variable values for the 
!    parallel averaging calculation
! -----------------------------------------------------
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  class(srcsink_panel_type), pointer :: cwp
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: i
  PetscReal :: local_variables(21)
  PetscReal :: global_variables(21)
! -----------------------------------------------------
  
  if (.not.associated(this%waste_panel_list)) return
  
100 format(100es18.8)  ! double
101 format(1I6.1)      ! integer

  option => this%realization%option
  output_option => this%realization%output_option
  
  ! open file and set write action to append
  fid = 88
  filename = PMWSSOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")
       
  ! write time in user's time units
  write(fid,100,advance="no") option%time / output_option%tconv  ! [T units]
  
  cwp => this%waste_panel_list
  do
    if (.not.associated(cwp)) exit
  ! ----- pre-calculations for writing data ---------------------------------
    local_variables = 0.d0
    global_variables = 0.d0
    do i = 1,cwp%region%num_cells            
      ! H2RATE [mol-gas/m3-bulk/sec] index=1
      local_variables(1) = local_variables(1) + &    
        (cwp%gas_generation_rate(i)*cwp%scaling_factor(i)) 
      ! BRINRATE [mol-brine/m3-bulk/sec] index=2
      local_variables(2) = local_variables(2) + &  
        (cwp%brine_generation_rate(i)*cwp%scaling_factor(i)) 
      ! FERATE [mol-Fe/m3-bulk/sec] index=3
      local_variables(3) = local_variables(3) + &
        (cwp%inventory%Fe_s%inst_rate(i)*cwp%scaling_factor(i))
      ! CELLRATE [mol-Cell/m3-bulk/sec] index=4
      local_variables(4) = local_variables(4) + &
        (cwp%inventory%BioDegs_s%inst_rate(i)*cwp%scaling_factor(i))
      ! FEOH2R [mol-FeOH2/m3-bulk/sec] index=5
      local_variables(5) = local_variables(5) + &
        (cwp%inventory%FeOH2_s%inst_rate(i)*cwp%scaling_factor(i))
      ! FESR [mol-FeS/m3-bulk/sec] index=6
      local_variables(6) = local_variables(6) + &
        (cwp%inventory%FeS_s%inst_rate(i)*cwp%scaling_factor(i))
      ! MGOR [mol-MgO/m3-bulk/sec] index=7
      local_variables(7) = local_variables(7) + &
        (cwp%inventory%MgO_s%inst_rate(i)*cwp%scaling_factor(i))
      ! MGOH2R [mol-Mg(OH)2/m3-bulk/sec] index=8
      local_variables(8) = local_variables(8) + &
        (cwp%inventory%MgOH2_s%inst_rate(i)*cwp%scaling_factor(i))
      ! HYMAGR [mol-Hydromagnesite/m3-bulk/sec] index=9
      local_variables(9) = local_variables(9) + &
        (cwp%inventory%Mg5CO34OH24H2_s%inst_rate(i)*cwp%scaling_factor(i))    
      ! MGCO3R [mol-MgCO3/m3-bulk/sec] index=10
      local_variables(10) = local_variables(10) + &
        (cwp%inventory%MgCO3_s%inst_rate(i)*cwp%scaling_factor(i))  
      
      ! CORRATI [mol/m3-bulk/sec] index=11
      local_variables(11) = local_variables(11) + &
        (cwp%rxnrate_Fe_corrosion_inund(i)*cwp%scaling_factor(i))
      ! CORRATH [mol/m3-bulk/sec] index=12
      local_variables(12) = local_variables(12) + &
        (cwp%rxnrate_Fe_corrosion_humid(i)*cwp%scaling_factor(i))
      ! BIORATI [mol/m3-bulk/sec] index=13
      local_variables(13) = local_variables(13) + &
        (cwp%rxnrate_cell_biodeg_inund(i)*cwp%scaling_factor(i))
      ! BIORATH [mol/m3-bulk/sec] index=14
      local_variables(14) = local_variables(14) + &
        (cwp%rxnrate_cell_biodeg_humid(i)*cwp%scaling_factor(i))
      ! FEOH2_SR [mol/m3-bulk/sec] index=15
      local_variables(15) = local_variables(15) + &
        (cwp%rxnrate_FeOH2_sulf(i)*cwp%scaling_factor(i))
      ! FE_SR [mol/m3-bulk/sec] index=16
      local_variables(16) = local_variables(16) + &
        (cwp%rxnrate_Fe_sulf(i)*cwp%scaling_factor(i))
      ! MGO_HR [mol/m3-bulk/sec] index=17
      local_variables(17) = local_variables(17) + &
        (cwp%rxnrate_MgO_hyd(i)*cwp%scaling_factor(i))
      ! MGOH2_CR [mol/m3-bulk/sec] index=18
      local_variables(18) = local_variables(18) + &
        (cwp%rxnrate_MgOH2_carb(i)*cwp%scaling_factor(i))
      ! MGO_CR [mol/m3-bulk/sec] index=19
      local_variables(19) = local_variables(19) + &
        (cwp%rxnrate_MgO_carb(i)*cwp%scaling_factor(i))
      ! HYMAG_CR [mol/m3-bulk/sec] index=20
      local_variables(20) = local_variables(20) + &
        (cwp%rxnrate_hydromag_conv(i)*cwp%scaling_factor(i))
      
      ! PORSOLID [-] index=21
      local_variables(21) = local_variables(21) + &
        (cwp%solids_production(i)*cwp%scaling_factor(i))
    enddo
    call CalcParallelSUM(option,cwp%rank_list,local_variables,global_variables)
  ! ----- write data to file ------------------------------------------------
  ! ----- see table 7 in BRAGFLO's User Manual, pg. 86 ----------------------
  ! -------------------------------------------------------------------------
  ! PNL_ID [-]
    write(fid,101,advance="no") cwp%id 
  ! CORRATI [mol-Fe/sec] #39 inundated iron corrosion rate
    write(fid,100,advance="no") global_variables(11)*cwp%volume
  ! CORRATH [mol-Fe/sec] #40 humid iron corrosion rate
    write(fid,100,advance="no") global_variables(12)*cwp%volume
  ! BIORATI [mol-Cell/sec] #41 inundated biodegradation rate
    write(fid,100,advance="no") global_variables(13)*cwp%volume
  ! BIORATH [mol-Cell/sec] #42 humid biodegradation rate
    write(fid,100,advance="no") global_variables(14)*cwp%volume
  ! FEOH2_SR [mol/sec] #43 iron hydroxide sulfidation rate
    write(fid,100,advance="no") global_variables(15)*cwp%volume
  ! FE_SR [mol/sec] #44 iron sulfidation rate
    write(fid,100,advance="no") global_variables(16)*cwp%volume
  ! MGO_HR [mol/sec] #45 MgO hydration rate
    write(fid,100,advance="no") global_variables(17)*cwp%volume
  ! MGOH2_CR [mol/sec] #46 Mg(OH)2 carbonation rate
    write(fid,100,advance="no") global_variables(18)*cwp%volume
  ! MGO_CR [mol/sec] #47 MgO carbonation rate
    write(fid,100,advance="no") global_variables(19)*cwp%volume
  ! HYMAG_CR [mol/sec] #48 hydromagnesite conversion rate
    write(fid,100,advance="no") global_variables(20)*cwp%volume
    
  ! H2RATE [kg/m3-bulk/sec] #49 gas generation rate
    write(fid,100,advance="no") global_variables(1)*MW_H2
  ! BRINRATE [kg/m3-bulk/sec] #50 brine consumption rate
    ! actually, bragflo uses brine weight (with salt), not water weight
    write(fid,100,advance="no") global_variables(2)*MW_H2O
  ! FERATE [kg/m3-bulk/sec] #51 Fe consumption rate
    write(fid,100,advance="no") global_variables(3)*MW_FE
  ! CELLRATE [kg/m3-bulk/sec] #52 Biodegradables consumption rate
    write(fid,100,advance="no") global_variables(4)*MW_CELL
  ! FEOH2R [kg-FeOH2/m3-bulk/sec] #53 Fe(OH)2 consumption rate
    write(fid,100,advance="no") global_variables(5)*MW_FEOH2
  ! FESR [kg-FeS/m3-bulk/sec] #54 FeS consumption rate
    write(fid,100,advance="no") global_variables(6)*MW_FES
  ! MGOR [kg-MgO/m3-bulk/sec] #55 MgO consumption rate
    write(fid,100,advance="no") global_variables(7)*MW_MGO
  ! MGOH2R [kg-Mg(OH)2/m3-bulk/sec] #56 Mg(OH)2 consumption rate
    write(fid,100,advance="no") global_variables(8)*MW_MGOH2
  ! HYMAGR [kg-Hymag/m3-bulk/sec] #57 Hydromagnesite consumption rate
    write(fid,100,advance="no") global_variables(9)*MW_HYDRO
  ! MGCO3R [kg-MgCO3/m3-bulk/sec] #58 MgCO3 consumption rate
    write(fid,100,advance="no") global_variables(10)*MW_MGCO3
  ! FECONC [kg/m3-bulk] #59 Fe concentration
    write(fid,100,advance="no") &
      (cwp%inventory%Fe_s%tot_mass_in_panel)/cwp%volume
  ! CELLCONC [kg/m3-bulk] #60 Biodegradables concentration
    write(fid,100,advance="no") &
      (cwp%inventory%BioDegs_s%tot_mass_in_panel)/cwp%volume
  ! FEOH2C [kg/m3-bulk] #61 Fe(OH)2 concentration
    write(fid,100,advance="no") &
      (cwp%inventory%FeOH2_s%tot_mass_in_panel)/cwp%volume
  ! FESC [kg/m3-bulk] #62 FeS concentration
    write(fid,100,advance="no") &
      (cwp%inventory%FeS_s%tot_mass_in_panel)/cwp%volume
  ! MGOC [kg/m3-bulk] #63 MgO concentration
    write(fid,100,advance="no") &
      (cwp%inventory%MgO_s%tot_mass_in_panel)/cwp%volume
  ! MGOH2C [kg/m3-bulk] #64 Mg(OH)2 concentration
    write(fid,100,advance="no") &
      (cwp%inventory%MgOH2_s%tot_mass_in_panel)/cwp%volume
  ! HYMAGC [kg/m3-bulk] #65 Hydromagnesite concentration
    write(fid,100,advance="no") &
      (cwp%inventory%Mg5CO34OH24H2_s%tot_mass_in_panel)/cwp%volume
  ! MGCO3C [kg/m3-bulk] #66 MgCO3 concentration
    write(fid,100,advance="no") &
      (cwp%inventory%MgCO3_s%tot_mass_in_panel)/cwp%volume
      
  ! PORSOLID [-] #68 Solids production normalized by panel volume
    write(fid,100,advance="no") global_variables(21)
  !----------------------------------------
    cwp => cwp%next
  enddo
  
  close(fid)
  
end subroutine PMWSSOutput

! ************************************************************************** !

function PMWSSOutputFilename(option)
  ! 
  ! Generates filename for waste panel output.
  ! 
  ! Author: Jennifer Frederick
  ! Date: 03/27/17
  !
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
! PMWSSOutputFilename (output): filename string
! word: temporary string
! -----------------------------------------------------
  character(len=MAXSTRINGLENGTH) :: PMWSSOutputFilename
  character(len=MAXWORDLENGTH) :: word
! -----------------------------------------------------

  write(word,'(i6)') option%myrank
  PMWSSOutputFilename = trim(option%global_prefix) // &
                       trim(option%group_prefix) // &
                       '-' // trim(adjustl(word)) // '.pnl'
  
end function PMWSSOutputFilename  

! ************************************************************************** !

subroutine PMWSSOutputHeader(this)
  ! 
  ! Writes the header for a waste panel output file.
  ! 
  ! Author: Jennifer Frederick
  ! Date: 03/27/17

  use Output_Aux_module
  use Utility_module
    
  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
! -----------------------------------
  
! LOCAL VARIABLES:
! ================
! output_option: pointer to output option object
! cur_waste_panel: pointer to current waste panel object
! filename: filename string
! units_string: unit string
! variable_string: variable string
! cell_string: cell string
! fid: [-] file id integer
! icolumn: [-] column integer
! exist: Boolean to check if file exists
! -----------------------------------------------------
  type(output_option_type), pointer :: output_option
  class(srcsink_panel_type), pointer :: cur_waste_panel
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: units_string 
  character(len=MAXWORDLENGTH) :: variable_string
  character(len=MAXWORDLENGTH) :: cell_string
  PetscInt :: fid
  PetscInt :: icolumn
  PetscBool :: exist
! -----------------------------------------------------
  
  if (.not.associated(this%waste_panel_list)) return
  
  output_option => this%realization%output_option
  
  fid = 88
  filename = PMWSSOutputFilename(this%option)
  ! check if file exists
  exist = FileExists(trim(filename))
  if (this%option%restart_flag .and. exist) return
  open(unit=fid,file=filename,action="write",status="replace")
  
  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif
   
  ! write time in user's time units
  write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) &
                                // ']"'
  
  ! write the rest of the data
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    cell_string = trim(cur_waste_panel%region_name)
    variable_string = 'PNL_ID#'
    units_string = ''
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'CORRATI'  ! #39
    units_string = 'mol-Fe/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'CORRATH'  ! #40
    units_string = 'mol-Fe/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'BIORATI'  ! #41
    units_string = 'mol-Cell/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'BIORATH'  ! #42
    units_string = 'mol-Cell/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'FEOH2_SR'  ! #43
    units_string = 'mol/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'FE_SR'  ! #44
    units_string = 'mol/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MGO_HR'  ! #45
    units_string = 'mol/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MGOH2_CR'  ! #46
    units_string = 'mol/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MGO_CR'  ! #47
    units_string = 'mol/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'HYMAG_CR'  ! #48
    units_string = 'mol/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'H2RATE'  ! #49
    units_string = 'kg/m3/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'BRINRATE'  ! #50
    units_string = 'kg/m3/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'FERATE'  ! #51
    units_string = 'kg/m3/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'CELLRATE'  ! #52
    units_string = 'kg/m3/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'FEOH2R'  ! #53
    units_string = 'kg/m3/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'FESR'  ! #54
    units_string = 'kg/m3/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MGOR'  ! #55
    units_string = 'kg/m3/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MGOH2R'  ! #56
    units_string = 'kg/m3/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'HYMAGR'  ! #57
    units_string = 'kg/m3/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MGCO3R'  ! #58
    units_string = 'kg/m3/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)                         
    variable_string = 'FECONC'  ! #59
    units_string = 'kg/m3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'CELLCONC'  ! #60
    units_string = 'kg/m3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'FEOH2C'  ! #61
    units_string = 'kg/m3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'FESC'  ! #62
    units_string = 'kg/m3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MGOC'  ! #63
    units_string = 'kg/m3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MGOH2C'  ! #64
    units_string = 'kg/m3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'HYMAGC'  ! #65
    units_string = 'kg/m3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MGCO3C'  ! #66
    units_string = 'kg/m3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'PORSOLID'  ! #68
    units_string = '-'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
  !----------------------------------------
    cur_waste_panel => cur_waste_panel%next
  enddo
  
  close(fid)
  
end subroutine PMWSSOutputHeader

! *************************************************************************** !

subroutine PMWSSInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  ! 
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
! -----------------------------------

! LOCAL VARIABLES:
! ================
! id: filename id for the *.rec file
! --------------
  PetscInt :: id
! --------------

  id = INPUT_RECORD_UNIT
  
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMWSSInputRecord

! *************************************************************************** !

subroutine PMWSSDestroy(this)
  ! 
  ! Strips and destroys the WIPP source/sink process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  !

  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
! -----------------------------------
  
  call PMWSSStrip(this)
  
end subroutine PMWSSDestroy

! ************************************************************************** !

subroutine PMWSSStrip(this)
  ! 
  ! Strips the WIPP source/sink process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  !
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! this (input/output): wipp-srcsink process model object
! -----------------------------------
  class(pm_wipp_srcsink_type) :: this
! -----------------------------------
  
! LOCAL VARIABLES:
! ================
! cur_waste_panel: pointer to current waste panel object
! prev_waste_panel: pointer to previous waste panel object
! cur_preinventory: pointer to current pre-inventory object
! prev_preinventory: pointer to previous pre-inventory object
! ------------------------------------------------------------------------
  type(srcsink_panel_type), pointer :: cur_waste_panel, prev_waste_panel
  type(pre_inventory_type), pointer :: cur_preinventory, prev_preinventory
! ------------------------------------------------------------------------

  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    prev_waste_panel => cur_waste_panel
    cur_waste_panel => cur_waste_panel%next
    call PMWSSDestroyPanel(prev_waste_panel)
  enddo
  nullify(this%waste_panel_list)
  
  cur_preinventory => this%pre_inventory_list
  do
    if (.not.associated(cur_preinventory)) exit
    prev_preinventory => cur_preinventory
    cur_preinventory => cur_preinventory%next
    deallocate(prev_preinventory)
    nullify(prev_preinventory)
  enddo
  nullify(this%pre_inventory_list)

end subroutine PMWSSStrip

! ************************************************************************** !

subroutine PMWSSDestroyPanel(waste_panel)
  ! 
  ! Strips a waste panel in the WIPP source/sink process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2017
  !
  use Utility_module, only : DeallocateArray
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! waste_panel (input/output): pointer to waste panel object
! ------------------------------------------------
  type(srcsink_panel_type), pointer :: waste_panel
! ------------------------------------------------
  
  if (.not.associated(waste_panel)) return
  call PMWSSInventoryDestroy(waste_panel%inventory)
  call DeallocateArray(waste_panel%scaling_factor)
  call DeallocateArray(waste_panel%solids_production)
  call DeallocateArray(waste_panel%gas_generation_rate)
  call DeallocateArray(waste_panel%brine_generation_rate)
  call DeallocateArray(waste_panel%rxnrate_Fe_corrosion_inund)
  call DeallocateArray(waste_panel%rxnrate_Fe_corrosion_humid)
  call DeallocateArray(waste_panel%rxnrate_cell_biodeg_inund)
  call DeallocateArray(waste_panel%rxnrate_cell_biodeg_humid)
  call DeallocateArray(waste_panel%rxnrate_MgO_hyd_inund)
  call DeallocateArray(waste_panel%rxnrate_MgO_hyd_humid)
  call DeallocateArray(waste_panel%rxnrate_Fe_corrosion)
  call DeallocateArray(waste_panel%rxnrate_cell_biodeg)
  call DeallocateArray(waste_panel%rxnrate_FeOH2_sulf)
  call DeallocateArray(waste_panel%rxnrate_Fe_sulf)
  call DeallocateArray(waste_panel%rxnrate_MgO_hyd)
  call DeallocateArray(waste_panel%rxnrate_MgOH2_carb)
  call DeallocateArray(waste_panel%rxnrate_MgO_carb)
  call DeallocateArray(waste_panel%rxnrate_hydromag_conv)
  
  call DeallocateArray(waste_panel%rank_list)
  nullify(waste_panel%next)
  deallocate(waste_panel)
  nullify(waste_panel)
    
end subroutine PMWSSDestroyPanel

! ************************************************************************** !

subroutine PMWSSInventoryDestroy(inventory)
  ! 
  ! Deallocates the inventory's chemical species.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/01/2017
  !
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! inventory (input/output): inventory object
! ---------------------------------
  type(inventory_type) :: inventory
! ---------------------------------
  
  call PMWSSChemSpeciesDeallocate(inventory%Fe_s)
  call PMWSSChemSpeciesDeallocate(inventory%FeOH2_s)
  call PMWSSChemSpeciesDeallocate(inventory%BioDegs_s)
  call PMWSSChemSpeciesDeallocate(inventory%FeS_s)
  call PMWSSChemSpeciesDeallocate(inventory%MgO_s)
  call PMWSSChemSpeciesDeallocate(inventory%MgOH2_s)
  call PMWSSChemSpeciesDeallocate(inventory%Mg5CO34OH24H2_s)
  call PMWSSChemSpeciesDeallocate(inventory%MgCO3_s)
  nullify(inventory%preinventory)
  
end subroutine PMWSSInventoryDestroy

! ************************************************************************** !

subroutine PMWSSChemSpeciesDeallocate(chem_species)
  ! 
  ! Deallocates the inventory's chemical species.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/01/2017
  !
  use Utility_module, only : DeallocateArray
  
  implicit none
  
! INPUT ARGUMENTS:
! ================
! chem_species (input/output): chemical species object
! ---------------------------------------
  type(chem_species_type) :: chem_species
! ---------------------------------------
  
  call DeallocateArray(chem_species%initial_conc_mol)
  call DeallocateArray(chem_species%current_conc_mol)
  call DeallocateArray(chem_species%current_conc_kg)
  call DeallocateArray(chem_species%inst_rate)
  
end subroutine PMWSSChemSpeciesDeallocate

! *************************************************************************** !

end module PM_WIPP_SrcSink_class
