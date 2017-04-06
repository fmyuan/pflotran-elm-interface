module PM_WIPP_SrcSink_class

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PM_Base_class
  use Region_module
  use PFLOTRAN_Constants_module
  use Realization_Subsurface_class
  use Data_Mediator_Vec_class

  implicit none

  private

  type, public :: chem_species_type
    PetscReal, pointer :: initial_conc_mol(:)  ! [mol/m3-bulk]
    PetscReal, pointer :: inst_rate(:)         ! [mol/m3-bulk/sec]
    PetscReal, pointer :: current_conc_mol(:)  ! [mol/m3-bulk]
    PetscReal, pointer :: current_conc_kg(:)   ! [kg/m3-bulk] per BRAGFLO U.M.
    PetscReal :: molar_mass                    ! [kg/mol]
    PetscReal :: tot_mass_in_panel             ! [kg/panel-volume]
  end type chem_species_type
  
  type, public :: inventory_type
    character(len=MAXWORDLENGTH) :: name
    type(chem_species_type) :: Fe_s 
    type(chem_species_type) :: FeOH2_s
    type(chem_species_type) :: C6H10O5_s
    type(chem_species_type) :: RuPl_s
    type(chem_species_type) :: NO3_minus_aq
    type(chem_species_type) :: CO2_g
    type(chem_species_type) :: N2_g
    type(chem_species_type) :: SO42_minus_aq
    type(chem_species_type) :: H2S_g
    type(chem_species_type) :: FeS_s
    type(chem_species_type) :: MgO_s
    type(chem_species_type) :: MgOH2_s
    type(chem_species_type) :: Mg5CO34OH24H2_s
    type(chem_species_type) :: MgCO3_s
  ! Initial Values:
    PetscReal :: Fe_in_panel          ! [total initial kg in waste panel]
    PetscReal :: MgO_in_panel         ! [total initial kg in waste panel]
    PetscReal :: Cellulose_in_panel   ! [total initial kg in waste panel]
    PetscReal :: RubberPlas_in_panel  ! [total initial kg in waste panel]
    PetscReal :: Biodegs_in_panel     ! [total initial kg in waste panel]
    PetscReal :: Nitrate_in_panel     ! [total initial kg in waste panel]
    PetscReal :: Sulfate_in_panel     ! [total initial kg in waste panel]
    PetscReal :: drum_conc            ! [number of steel drums per m3 waste panel]
    type(pre_inventory_type), pointer :: preinventory
  end type inventory_type
  
  ! pre-inventory describes initial waste; what is emplaced in a waste panel
  type, public :: pre_inventory_type
    character(len=MAXWORDLENGTH) :: name
  ! ALGEBRA parameters:
    PetscReal :: ironchw    ! [kg] mass of Fe-based material in CH waste
    PetscReal :: ironrhw    ! [kg] mass of Fe-based material in RH waste
    PetscReal :: irncchw    ! [kg] mass of Fe containers for CH waste
    PetscReal :: irncrhw    ! [kg] mass of Fe containers for RH waste
    PetscReal :: cellchw    ! [kg] mass of cellulosics in CH waste
    PetscReal :: cellrhw    ! [kg] mass of cellulosics in RH waste
    PetscReal :: celcchw    ! [kg] mass of cellulosics in container materials for CH waste
    PetscReal :: celcrhw    ! [kg] mass of cellulosics in container materials for RH waste
    PetscReal :: celechw    ! [kg] mass of cellulosics in emplacement materials for CH waste
    PetscReal :: celerhw    ! [kg] mass of cellulosics in emplacement materials for RH waste
    PetscReal :: rubbchw    ! [kg] mass of rubber in CH waste
    PetscReal :: rubbrhw    ! [kg] mass of rubber in RH waste
    PetscReal :: rubcchw    ! [kg] mass of rubber in container materials for CH waste
    PetscReal :: rubcrhw    ! [kg] mass of rubber in container materials for RH waste
    PetscReal :: rubechw    ! [kg] mass of rubber in emplacement materials for CH waste
    PetscReal :: ruberhw    ! [kg] mass of rubber in emplacement materials for RH waste  
    PetscReal :: plaschw    ! [kg] mass of plastics in CH waste
    PetscReal :: plasrhw    ! [kg] mass of plastics in RH waste
    PetscReal :: plscchw    ! [kg] mass of plastics in container materials for CH waste
    PetscReal :: plscrhw    ! [kg] mass of plastics in container materials for RH waste
    PetscReal :: plsechw    ! [kg] mass of plastics in emplacement materials for CH waste
    PetscReal :: plserhw    ! [kg] mass of plastics in emplacement materials for RH waste 
    PetscReal :: plasfac    ! [-] mass ratio of plastics to equivalent carbon
    PetscReal :: mgo_ef     ! [-] MgO excess factor; ratio mol-MgO/mol-organic-C
    PetscReal :: vrepos     ! [m3] volume of total repository (not always needed)
    PetscReal :: drum_conc  ! [-/m3] number of steel drums per m3 waste panel
    PetscReal :: nitrate    ! [kg] mass of nitrate
    PetscReal :: sulfate    ! [kg] mass of sulfate
    type(pre_inventory_type), pointer :: next
  end type pre_inventory_type
  
  type, public :: srcsink_panel_type
    character(len=MAXWORDLENGTH) :: name
    type(region_type), pointer :: region
    character(len=MAXWORDLENGTH) :: region_name
    type(inventory_type) :: inventory
    character(len=MAXWORDLENGTH) :: inventory_name
    PetscReal, pointer :: scaling_factor(:)        ! [-]
    PetscReal, pointer :: gas_generation_rate(:)   ! [mol/m3-bulk/sec]
    PetscReal, pointer :: brine_generation_rate(:) ! [mol/m3-bulk/sec]
    PetscReal :: inundated_corrosion_rate          ! [mol/m3-bulk/sec]
    PetscReal :: humid_corrosion_rate              ! [mol/m3-bulk/sec], [-]
    PetscReal :: inundated_biodeg_rate             ! [mol/m3-bulk/sec]
    PetscReal :: humid_biodeg_rate                 ! [mol/m3-bulk/sec]
    PetscReal :: inundated_brucite_rate            ! [mol/m3-bulk/sec]
    PetscReal :: humid_brucite_rate                ! [mol/m3-bulk/sec]
    PetscReal :: RXH2S_factor                      ! [-]
    PetscReal :: volume                            ! [m3]
    PetscBool :: scale_by_volume                   ! flag to scale given inventory to waste panel volume
    PetscInt :: id
    PetscMPIInt :: myMPIcomm
    PetscMPIInt :: myMPIgroup
    PetscInt, pointer :: rank_list(:)
    type(srcsink_panel_type), pointer :: next
  end type srcsink_panel_type

  type, public, extends(pm_base_type) :: pm_wipp_srcsink_type
    PetscReal :: alpharxn           ! [-] 
    PetscReal :: smin               ! [-]
    PetscReal :: satwick            ! [-]
    PetscReal :: corrmco2           ! [m/s]
    PetscReal :: humcorr            ! [m/s]
    PetscReal :: gratmici           ! [mol/kg/sec]
    PetscReal :: gratmich           ! [mol/kg/sec]
    PetscReal :: brucites           ! [mol/kg/sec]
    PetscReal :: brucitec           ! [mol/kg/sec]
    PetscReal :: bruciteh           ! [mol/kg/sec]
    PetscReal :: RXCO2_factor       ! [-]
    PetscReal :: hymagcon_rate      ! [mol/kg/sec]
    PetscReal :: drum_surface_area  ! [m2/drum]
    PetscReal :: biogenfc           ! [-]
    PetscInt :: probdeg             ! [-]
    PetscInt :: bioidx              ! [-] flag
    PetscInt :: plasidx             ! [-] flag
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
  
  class(pm_wipp_srcsink_type), pointer :: PMWSSCreate
  class(pm_wipp_srcsink_type), pointer :: pm
  
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
  pm%brucites = UNINITIALIZED_DOUBLE
  pm%brucitec = UNINITIALIZED_DOUBLE
  pm%bruciteh = UNINITIALIZED_DOUBLE
  pm%RXCO2_factor = UNINITIALIZED_DOUBLE
  pm%hymagcon_rate = UNINITIALIZED_DOUBLE
  pm%drum_surface_area = UNINITIALIZED_DOUBLE
  pm%biogenfc = UNINITIALIZED_DOUBLE
  pm%probdeg = UNINITIALIZED_INTEGER
  pm%bioidx = UNINITIALIZED_INTEGER
  pm%plasidx = UNINITIALIZED_INTEGER
  
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
  
  type(srcsink_panel_type), pointer :: PMWSSWastePanelCreate
  type(srcsink_panel_type), pointer :: panel
  
  allocate(panel)
  
  nullify(panel%next)
  nullify(panel%region)
  nullify(panel%scaling_factor)
  nullify(panel%gas_generation_rate)
  nullify(panel%brine_generation_rate)
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
  panel%RXH2S_factor = UNINITIALIZED_DOUBLE
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
  
  type(pre_inventory_type), pointer :: PMWSSPreInventoryCreate
  type(pre_inventory_type), pointer :: preinv
  
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
  
  type(inventory_type) :: inventory
  
  PetscReal :: molar_mass
  
  nullify(inventory%preinventory)
  inventory%name = ''
  
  molar_mass = (0.055847d0) ! [kg/mol MW_FE]
  call PMWSSInitChemSpecies(inventory%Fe_s,molar_mass)  ! iron
  
  molar_mass = (0.08986d0) ! [kg/mol MW_FEOH2]
  call PMWSSInitChemSpecies(inventory%FeOH2_s,molar_mass) ! iron hydroxide
  
  molar_mass = (0.027023d0) ! [kg/mol MW_CELL]
  call PMWSSInitChemSpecies(inventory%C6H10O5_s,molar_mass)  ! cellulose
  
  molar_mass = (0.027023d0) ! [rubber/plastic kg/mol, same as MW_CELL]
  call PMWSSInitChemSpecies(inventory%RuPl_s,molar_mass)  ! rubber/plastics
  
  molar_mass = (14.0067d-3 + 3.d0*15.9994d-3) ! [NO3- kg/mol]
  call PMWSSInitChemSpecies(inventory%NO3_minus_aq,molar_mass)  ! nitrate
  
  molar_mass = (0.0440098d0) ! [kg/mol MW_CO2]
  call PMWSSInitChemSpecies(inventory%CO2_g,molar_mass)  ! carbon dioxide gas 
  
  molar_mass = (0.02801348d0) ! [kg/mol MW_N2]
  call PMWSSInitChemSpecies(inventory%N2_g,molar_mass)  ! nitrogen gas  
  
  molar_mass = (32.065d-3 + 4.d0*15.9994d-3) ! [SO42- kg/mol]
  call PMWSSInitChemSpecies(inventory%SO42_minus_aq,molar_mass)  ! sulfate
  
  molar_mass = (0.03408188d0) ! [kg/mol MW_H2S]
  call PMWSSInitChemSpecies(inventory%H2S_g,molar_mass)  
  
  molar_mass = (0.087911d0) ! [kg/mol MW_FES]
  call PMWSSInitChemSpecies(inventory%FeS_s,molar_mass)  ! iron sulfide 
  
  molar_mass = (0.040304d0) ! [kg/mol MW_MGO]
  call PMWSSInitChemSpecies(inventory%MgO_s,molar_mass)  ! magnesium oxide
  
  molar_mass = (0.05832d0) ! [kg/mol MW_MGOH2]
  call PMWSSInitChemSpecies(inventory%MgOH2_s,molar_mass)  ! magnesium hydroxide              
  
  molar_mass = (5.d0*24.305d-3 + 4.d0*12.0107d-3 + 4.d0*3.d0*15.9994d-3 + &
                2.d0*15.9994d-3 + 2.d0*1.01d-3 + 8.d0*1.01d-3 + &
                4.d0*15.9994d-3) ! [Mg5CO34OH24H2 kg/mol]
  call PMWSSInitChemSpecies(inventory%Mg5CO34OH24H2_s,molar_mass)  ! hydromagnesite
  
  molar_mass = (0.084314d0) ! [kg/mol MW_MGCO3]
  call PMWSSInitChemSpecies(inventory%MgCO3_s,molar_mass)  ! magnesium carbonate
  
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
  
  type(chem_species_type) :: chem_species
  PetscReal :: molar_mass
  
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
  
  class(pm_wipp_srcsink_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
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
  
  class(pm_wipp_srcsink_type) :: this
  type(region_list_type), pointer :: region_list
  
  type(region_type), pointer :: cur_region
  class(srcsink_panel_type), pointer :: cur_waste_panel
  type(option_type), pointer :: option
  
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
  
  class(pm_wipp_srcsink_type) :: this
  
  type(pre_inventory_type), pointer :: cur_preinventory
  class(srcsink_panel_type), pointer :: cur_waste_panel
  type(option_type), pointer :: option
  
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
  
  type(pre_inventory_type), pointer :: preinventory
  type(inventory_type) :: inventory
  
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
  
  class(pm_wipp_srcsink_type) :: this
  type(srcsink_panel_type), pointer :: waste_panel
  
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type), pointer :: grid
  PetscInt :: k
  PetscInt :: local_id, ghosted_id
  PetscReal :: total_volume_local, total_volume_global
  PetscErrorCode :: ierr
  
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
  
  class(pm_wipp_srcsink_type) :: this
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: double
  character(len=MAXSTRINGLENGTH) :: error_string, error_string2, error_string3
  type(srcsink_panel_type), pointer :: new_waste_panel
  type(srcsink_panel_type), pointer :: cur_waste_panel
  type(pre_inventory_type), pointer :: new_inventory
  type(pre_inventory_type), pointer :: cur_preinventory
  PetscInt :: num_errors
  PetscBool :: added
  
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
      case('BRUCITEC')
        call InputReadDouble(input,option,this%brucitec)
        call InputErrorMsg(input,option,'MgO inundated hydration rate in &
                           &Castile brine (BRUCITEI)',error_string)
      case('BRUCITES')
        call InputReadDouble(input,option,this%brucites)
        call InputErrorMsg(input,option,'MgO inundated hydration rate in &
                           &Salado brine (BRUCITEI)',error_string)
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
              call InputErrorMsg(input,option,'inventory assignment',error_string)
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
  if (Uninitialized(this%brucitec)) then
    option%io_buffer = 'ERROR: BRUCITEI (MgO inundated hydration rate in &
                       &Castile brine) must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%brucites)) then
    option%io_buffer = 'ERROR: BRUCITEI (MgO inundated hydration rate in &
                       &Salado brine) must be specified in the &
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
  
  class(pm_wipp_srcsink_type) :: this
  type(srcsink_panel_type), pointer :: waste_panel
  
  type(pre_inventory_type), pointer :: preinventory
  type(inventory_type), pointer :: inventory
  PetscReal :: vol_scaling_factor
  PetscReal :: MOL_NO3                        ! moles of nitrate
  PetscReal :: F_NO3
  PetscReal :: MAX_C, A1, A2                  ! intermediate parameters
  PetscReal, parameter :: DN_FE = 7870.d0     ! [kg/m3] density of iron
  PetscReal, parameter :: MW_FE = 5.5847d-2   ! [kg/mol] mol weight of iron
  PetscReal, parameter :: MW_MGO = 4.0304d-2  ! [kg/mol] mol weight of MgO
  PetscReal, parameter :: MW_C = 2.7023d-2    ! [kg/mol] mol weight of cellulosics
  PetscReal, parameter :: MW_NO3 = 6.20d-2    ! [kg/mol] mol weight of nitrate
  PetscReal, parameter :: MW_SO4 = 9.60d-2    ! [kg/mol] mol weight of sulfate
  PetscReal :: D_c                            ! [kg/m3] mass conc biodegradables
  PetscReal :: D_m                            ! [kg/m3] mass conc MgO
  PetscReal :: D_s                            ! [m2/m3] area conc iron steel
  
  !-----PROBDEG-calculations-----------------------------------------------
  select case(this%probdeg)
    case(0)
      this%bioidx = 0
      this%plasidx = 0
    case(1)
      this%bioidx = 1
      this%plasidx = 0
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
        preinventory%mgo_ef * MW_MGO / MW_C                ! [-]
  inventory%Nitrate_in_panel = &                           ! [kg]
        preinventory%nitrate * MW_NO3                      ! [kg]
  inventory%Sulfate_in_panel = &                           ! [kg]
        preinventory%sulfate * MW_SO4                      ! [kg]
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
  !-----mass-concentrations----------------------------------units---------
  D_c = inventory%Biodegs_in_panel / &                     ! [kg]
        waste_panel%volume                                 ! [m3]
  D_s = this%drum_surface_area * &                         ! [m2]
        inventory%drum_conc                                ! [-/m3]
  D_m = inventory%MgO_in_panel / &                         ! [kg]
        waste_panel%volume                                 ! [m3]
  !-------------------------------------------------------------------------
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
  waste_panel%humid_corrosion_rate = &                     ! [-]
        waste_panel%humid_corrosion_rate / &
        waste_panel%inundated_corrosion_rate
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
  waste_panel%humid_biodeg_rate = &                     ! [-]
        waste_panel%humid_biodeg_rate / &
        waste_panel%inundated_biodeg_rate
  !-----iron-sulfidation----------------------------------------------------
  MOL_NO3 = inventory%Nitrate_in_panel / MW_NO3         ! [mol]
  A1 = inventory%Biodegs_in_panel / MW_C                ! [mol]
  A2 = this%gratmici * &                                ! [mol/kg/sec]
       (inventory%Biodegs_in_panel) * &                 ! [kg]
       (31556930.d0) * 10000.d0                         ! [sec/year]*[year]
  MAX_C = min(A1,A2)                                    ! [mol]
  F_NO3 = MOL_NO3 * (6.d0/4.8d0) / MAX_C                ! [-]
  F_NO3 = min(F_NO3,1.0)                                ! [-]
  waste_panel%RXH2S_factor = 1.0 - F_NO3                ! [-]
  !-----MgO-hydration-------------------------------------units-------------
  waste_panel%inundated_brucite_rate = &                ! [mol-bruc/m3/sec]
        max(this%brucites,this%bruciteh) * &            ! [mol-bruc/kg/sec]
        D_m                                             ! [kg/m3]
  waste_panel%humid_brucite_rate = &                    ! [mol-bruc/m3/sec]
        this%bruciteh * &                               ! [mol-bruc/kg/sec]
        D_m                                             ! [kg/m3]
  !-------------------------------------------------------------------------
  this%RXCO2_factor = 1.d0  ! BRAGFLO User's Manual Eq. 155, based on
                            ! Eqs. 145 & 146 stoichiometry 
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
  
  class(pm_wipp_srcsink_type) :: this
  
  type(option_type), pointer :: option
  type(srcsink_panel_type), pointer :: cur_waste_panel, prev_waste_panel
  type(srcsink_panel_type), pointer :: next_waste_panel
  PetscInt :: waste_panel_id
  PetscInt :: i, j
  PetscBool :: local
  PetscErrorCode :: ierr
  PetscMPIInt :: newcomm_size
  PetscInt, pointer :: ranks(:)
  
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
      allocate(cur_waste_panel%gas_generation_rate( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%gas_generation_rate(:) = 0.d0
      allocate(cur_waste_panel%brine_generation_rate( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%brine_generation_rate(:) = 0.d0
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
  
  type(inventory_type) :: inventory
  PetscInt :: num_cells
  PetscReal :: volume
  
  call PMWSSChemSpeciesAllocate(num_cells,inventory%Fe_s, &
                                inventory%Fe_in_panel,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%MgO_s, &
                                inventory%MgO_in_panel,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%RuPl_s, &
                                inventory%RubberPlas_in_panel,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%C6H10O5_s, &
                                inventory%Cellulose_in_panel,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%SO42_minus_aq, &
                                inventory%Sulfate_in_panel,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%NO3_minus_aq, &
                                inventory%Nitrate_in_panel,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%FeOH2_s,0.d0,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%CO2_g,0.d0,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%N2_g,0.d0,volume)
  call PMWSSChemSpeciesAllocate(num_cells,inventory%H2S_g,0.d0,volume)
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
  
  PetscInt :: num_cells
  type(chem_species_type) :: chem_species
  PetscReal :: initial_mass     ! mass species in waste panel [kg]
  PetscReal :: volume           ! volume of waste panel [m3]
  
  allocate(chem_species%initial_conc_mol(num_cells))
  allocate(chem_species%current_conc_mol(num_cells))
  allocate(chem_species%current_conc_kg(num_cells))
  allocate(chem_species%inst_rate(num_cells))
  
  chem_species%current_conc_kg(:) = initial_mass / volume             ! [kg/m3]
  chem_species%initial_conc_mol(:) = initial_mass / volume / &
                                     chem_species%molar_mass
  chem_species%current_conc_mol(:) = chem_species%initial_conc_mol(:) ! [mol/m3]
  chem_species%inst_rate(:) = 0.d0        
  chem_species%tot_mass_in_panel = initial_mass
  
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
  
  class(pm_wipp_srcsink_type) :: this
  
  IS :: is
  PetscInt :: i, j, k
  PetscErrorCode :: ierr
  PetscInt :: size_of_vec
  PetscInt, allocatable :: dofs_in_residual(:)
  type(srcsink_panel_type), pointer :: cur_waste_panel
  
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

  call PMWSSOutputHeader(this)
  call PMWSSOutput(this)
  call PMWSSSolve(this,0.d0,ierr)
  
end subroutine PMWSSInitializeRun

! *************************************************************************** !

subroutine PMWSSInitializeTimestep(this)
  ! 
  ! Initializes the process model to take a time step in the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  !
  
  implicit none

  class(pm_wipp_srcsink_type) :: this
  
  type(srcsink_panel_type), pointer :: cur_waste_panel
  PetscReal :: dt

  dt = this%option%flow_dt
  
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," WIPP SRC/SINK PANEL MODEL ",51("="))')
  endif
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    call PMWSSUpdateInventory(cur_waste_panel,dt,this%option)
    cur_waste_panel => cur_waste_panel%next
  enddo
  
  call PMWSSOutput(this)
  
end subroutine PMWSSInitializeTimestep

! *************************************************************************** !

subroutine PMWSSUpdateInventory(waste_panel,dt,option)
  !
  ! Updates the waste panel tracked species inventory concentrations.
  !
  ! Author: Jenn Frederick
  ! Date: 02/10/2017
  !
  
  use Option_module
 
  implicit none
  
  type(srcsink_panel_type) :: waste_panel
  PetscReal :: dt ! [sec; flow_dt]
  type(option_type) :: option
 
  call PMWSSUpdateChemSpecies(waste_panel%inventory%Fe_s,waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%FeOH2_s,waste_panel,dt,&
                              option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%C6H10O5_s,waste_panel,dt, &
                              option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%RuPl_s,waste_panel,dt, &
                              option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%NO3_minus_aq,waste_panel, &
                              dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%N2_g,waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%CO2_g,waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%SO42_minus_aq,waste_panel, &
                              dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%H2S_g,waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%FeS_s,waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%MgO_s,waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%MgOH2_s,waste_panel,dt, &
                              option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%Mg5CO34OH24H2_s, &
                              waste_panel,dt,option)
  call PMWSSUpdateChemSpecies(waste_panel%inventory%MgCO3_s,waste_panel,dt, &
                              option)
                                      
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
  
  implicit none
  
  type(chem_species_type) :: chem_species
  type(srcsink_panel_type) :: waste_panel
  PetscReal :: dt       ! [sec; flow_dt]
  type(option_type) :: option
  
  PetscInt :: k
  PetscInt :: num_cells
  PetscReal :: local_conc_kg, global_conc_kg
  
  num_cells = waste_panel%region%num_cells
  local_conc_kg = 0.d0  
  do k = 1,num_cells
    ! [mol/m3]
    chem_species%current_conc_mol(k) = &
                 chem_species%current_conc_mol(k) + &   ! [mol/m3]
                 chem_species%inst_rate(k) * &          ! [mol/m3/sec]
                 dt                                     ! [sec]
    ! [kg/m3]
    chem_species%current_conc_kg(k) = &
                 chem_species%current_conc_mol(k) * &   ! [mol/m3]
                 chem_species%molar_mass                ! [kg/mol]
    ! [kg/m3]             
    local_conc_kg = local_conc_kg + &
                (chem_species%current_conc_kg(k)*waste_panel%scaling_factor(k))
  enddo
  ! [kg]
  call PMWSSCalcParallelSUM(option,waste_panel,local_conc_kg,global_conc_kg)
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
  use Material_Aux_class
  use Global_Aux_module
  use EOS_Gas_module
  use EOS_Water_module
  
  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(general_auxvar_type), pointer :: gen_auxvar(:,:)
  type(global_auxvar_type), pointer :: global_auxvar(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, pointer :: vec_p(:)
  type(srcsink_panel_type), pointer :: cur_waste_panel
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id
  PetscInt :: num_cells
  ! brine/gas generation variable
  PetscReal :: water_saturation
  PetscReal :: s_eff
  PetscReal :: rxnrate_corrosion
  PetscReal :: rxnrate_biodeg_nitrate
  PetscReal :: rxnrate_biodeg_sulfate
  PetscReal :: rxnrate_FeS_Fe
  PetscReal :: rxnrate_FeS_FeOH2
  PetscReal :: rxnrate_mgoh2
  PetscReal :: rxnrate_hydromag
  PetscReal :: rxnrate_hymagcon
  ! enthalpy calculation variables
  PetscReal :: temperature
  PetscReal :: pressure_liq
  PetscReal :: pressure_gas
  PetscReal :: H_liq
  PetscReal :: H_gas
  PetscReal :: U_gas
  PetscReal :: gas_energy
  PetscReal :: brine_energy
  
  !return
  
  option => this%realization%option
  grid => this%realization%patch%grid
  gen_auxvar => this%realization%patch%aux%General%auxvars
  material_auxvars => this%realization%patch%aux%Material%auxvars
  global_auxvar => this%realization%patch%aux%Global%auxvars
  
  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
  j = 0  ! (j indexes the data mediator vec)
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    num_cells = cur_waste_panel%region%num_cells
    !-----zero-out-reaction-rates---------------------------------------------
    water_saturation = 0.d0
    s_eff = 0.d0
    rxnrate_corrosion = 0.d0
    rxnrate_biodeg_nitrate = 0.d0
    rxnrate_biodeg_sulfate = 0.d0
    rxnrate_FeS_Fe = 0.d0
    rxnrate_FeS_FeOH2 = 0.d0
    rxnrate_mgoh2 = 0.d0
    rxnrate_hydromag = 0.d0
    rxnrate_hymagcon = 0.d0

    do i = 1,num_cells
      local_id = cur_waste_panel%region%cell_ids(i)
      ghosted_id = grid%nL2G(local_id)
    !-----effective-brine-saturation------------------------------------------
      water_saturation = &
        gen_auxvar(ZERO_INTEGER,ghosted_id)%sat(option%liquid_phase)
      s_eff = water_saturation - this%smin + this%satwick*(1.d0 - &
        exp(200.d0*this%alpharxn*(max((water_saturation-this%smin),0.d0))**2.d0))
      if (s_eff > 1.d0) s_eff = 1.d0
    !-----anoxic-iron-corrosion-[mol-Fe/m3/sec]-------------------------------
      rxnrate_corrosion = (cur_waste_panel%inundated_corrosion_rate*s_eff) + &
                          (cur_waste_panel%humid_corrosion_rate*(1.d0-s_eff))
      call PMWSSSmoothRxnrate(rxnrate_corrosion,i, &
                              cur_waste_panel%inventory%Fe_s,this%alpharxn) 
      call PMWSSTaperRxnrate(rxnrate_corrosion,i,cur_waste_panel%inventory%Fe_s)
    !-----biodegradation-[mol-cell/m3/sec]------------------------------------
      rxnrate_biodeg_nitrate = (cur_waste_panel%inundated_biodeg_rate*s_eff) + &
                               (cur_waste_panel%humid_biodeg_rate*(1.d0-s_eff))
      rxnrate_biodeg_sulfate = (cur_waste_panel%inundated_biodeg_rate*s_eff) + &
                               (cur_waste_panel%humid_biodeg_rate*(1.d0-s_eff))
      call PMWSSSmoothRxnrate(rxnrate_biodeg_nitrate,i, &
                              cur_waste_panel%inventory%C6H10O5_s,this%alpharxn)
      call PMWSSSmoothRxnrate(rxnrate_biodeg_sulfate,i, &
                              cur_waste_panel%inventory%C6H10O5_s,this%alpharxn)
      call PMWSSTaperRxnrate(rxnrate_biodeg_nitrate,i, &
                             cur_waste_panel%inventory%C6H10O5_s, &
                             cur_waste_panel%inventory%NO3_minus_aq)
      call PMWSSTaperRxnrate(rxnrate_biodeg_sulfate,i, &
                             cur_waste_panel%inventory%C6H10O5_s, &
                             cur_waste_panel%inventory%SO42_minus_aq)
    !-----iron-sulfidation-[mol-H2S/m3/sec]-----------------------------------
      rxnrate_FeS_Fe = rxnrate_biodeg_sulfate*cur_waste_panel%RXH2S_factor
      rxnrate_FeS_FeOH2 = rxnrate_biodeg_sulfate*cur_waste_panel%RXH2S_factor
      call PMWSSSmoothRxnrate(rxnrate_FeS_Fe,i,cur_waste_panel%inventory%Fe_s, &
                              this%alpharxn) 
      call PMWSSSmoothRxnrate(rxnrate_FeS_FeOH2,i, &
                              cur_waste_panel%inventory%Fe_s,this%alpharxn)
      call PMWSSTaperRxnrate(rxnrate_FeS_Fe,i, &
                             cur_waste_panel%inventory%Fe_s, &
                             cur_waste_panel%inventory%H2S_g)
      call PMWSSTaperRxnrate(rxnrate_FeS_FeOH2,i, &
                             cur_waste_panel%inventory%FeOH2_s, &
                             cur_waste_panel%inventory%H2S_g)
    !-----MgO-hydration-[mol-MgO/m3/sec]--------------------------------------
      rxnrate_mgoh2 = (cur_waste_panel%inundated_brucite_rate*s_eff) + &
                      ((cur_waste_panel%humid_brucite_rate)*(1.d0-s_eff))
      call PMWSSSmoothRxnrate(rxnrate_mgoh2,i,cur_waste_panel%inventory%MgO_s, &
                              this%alpharxn)
      call PMWSSTaperRxnrate(rxnrate_mgoh2,i,cur_waste_panel%inventory%MgO_s)
    !-----hydromagnesite-[mol/m3-bulk/sec]------------------------------------
      rxnrate_hydromag = max(rxnrate_biodeg_nitrate,rxnrate_biodeg_sulfate) * &
                         this%RXCO2_factor
      call PMWSSSmoothRxnrate(rxnrate_hydromag,i, &
                              cur_waste_panel%inventory%MgO_s,this%alpharxn)
      call PMWSSTaperRxnrate(rxnrate_hydromag,i, &
                             cur_waste_panel%inventory%MgOH2_s)
    !-----hydromagnesite-conversion-[mol/m3-bulk/sec]-------------------------
      rxnrate_hymagcon = this%hymagcon_rate* &
                  cur_waste_panel%inventory%Mg5CO34OH24H2_s%current_conc_kg(i)
      call PMWSSSmoothRxnrate(rxnrate_hymagcon,i, &
                              cur_waste_panel%inventory%MgO_s,this%alpharxn)
      call PMWSSTaperRxnrate(rxnrate_hydromag,i, &
                             cur_waste_panel%inventory%Mg5CO34OH24H2_s)
    !-----tracked-species-[mol-species/m3-bulk/sec]---------------------------
      cur_waste_panel%inventory%FeOH2_s%inst_rate(i) = &
          1.d0*rxnrate_corrosion + (-1.d0*rxnrate_FeS_FeOH2)
      cur_waste_panel%inventory%Fe_s%inst_rate(i) = &
          (-1.d0*rxnrate_corrosion) + (-1.d0*rxnrate_FeS_Fe) 
      cur_waste_panel%inventory%FeS_s%inst_rate(i) = &
          1.d0*rxnrate_FeS_Fe + 1.d0*rxnrate_FeS_FeOH2
      cur_waste_panel%inventory%C6H10O5_s%inst_rate(i) = &
          (-1.d0*rxnrate_biodeg_nitrate) + (-1.d0*rxnrate_biodeg_sulfate)
      cur_waste_panel%inventory%RuPl_s%inst_rate(i) = &
          (-1.d0*rxnrate_biodeg_nitrate) + (-1.d0*rxnrate_biodeg_sulfate)
      cur_waste_panel%inventory%NO3_minus_aq%inst_rate(i) = &
          (-4.8d0*rxnrate_biodeg_nitrate)
      cur_waste_panel%inventory%CO2_g%inst_rate(i) = &
          6.d0*rxnrate_biodeg_sulfate + 6.d0*rxnrate_biodeg_nitrate + &
          (-4.d0*rxnrate_hydromag)
      cur_waste_panel%inventory%N2_g%inst_rate(i) = &
          2.4d0*rxnrate_biodeg_nitrate
      cur_waste_panel%inventory%SO42_minus_aq%inst_rate(i) = &
          (-3.d0*rxnrate_biodeg_sulfate)
      cur_waste_panel%inventory%H2S_g%inst_rate(i) = &
          3.d0*rxnrate_biodeg_sulfate + (-1.d0*rxnrate_FeS_Fe) + &
          (-1.d0*rxnrate_FeS_FeOH2)
      cur_waste_panel%inventory%MgO_s%inst_rate(i) = &
          (-1.d0*rxnrate_mgoh2)
      cur_waste_panel%inventory%MgOH2_s%inst_rate(i) = &
          1.d0*rxnrate_mgoh2 + (-5.d0*rxnrate_hydromag) + 1.d0*rxnrate_hymagcon
      cur_waste_panel%inventory%Mg5CO34OH24H2_s%inst_rate(i) = &
          1.d0*rxnrate_hydromag + (-1.d0*rxnrate_hymagcon) 
      cur_waste_panel%inventory%MgCO3_s%inst_rate(i) = &
          4.d0*rxnrate_hymagcon 
    !-----gas-generation-[mol-H2/m3-bulk/sec]---------------------------------
      cur_waste_panel%gas_generation_rate(i) = &
          1.d0*rxnrate_corrosion + 1.d0*rxnrate_FeS_Fe + &
          2.4d0*rxnrate_biodeg_nitrate
    !-----brine-generation-[mol-H2O/m3-bulk/sec]------------------------------
      cur_waste_panel%brine_generation_rate(i) = &
          (-2.d0*rxnrate_corrosion) + 7.4d0*rxnrate_biodeg_nitrate + &
          5.d0*rxnrate_biodeg_sulfate + 2.d0*rxnrate_FeS_FeOH2 + &
          (-1.d0*rxnrate_mgoh2) + 4.d0*rxnrate_hymagcon
    !------source-term-calculation--------------------------------------------
      j = j + 1
      ! liquid source term [kmol/sec]
      vec_p(j) = cur_waste_panel%brine_generation_rate(i) * &  ! [mol/m3/sec]
                 material_auxvars(ghosted_id)%volume / &       ! [m3-bulk]
                 1.d3                                          ! [mol -> kmol]
      j = j + 1
      ! gas source term [kmol/sec]
      vec_p(j) = cur_waste_panel%gas_generation_rate(i) * &    ! [mol/m3/sec]
                 material_auxvars(ghosted_id)%volume / &       ! [m3-bulk]
                 1.d3                                          ! [mol -> kmol]
      j = j + 1
      ! energy source term [MJ/sec]; H from EOS [J/kmol]
      brine_energy = 0.d0
      gas_energy = 0.d0
      temperature = gen_auxvar(ZERO_INTEGER,ghosted_id)%temp
      select case(global_auxvar(ghosted_id)%istate)
        case(GAS_STATE) !------------------------------------------------------
          pressure_gas = gen_auxvar(ZERO_INTEGER,ghosted_id)% &
                         pres(option%gas_phase)
          call EOSGasEnergy(temperature,pressure_gas,H_gas,U_gas,ierr)
          gas_energy = &
              cur_waste_panel%gas_generation_rate(i) * &    ! [mol/m3/sec]
              material_auxvars(ghosted_id)%volume * &       ! [m3-bulk] 
              H_gas * 1.d-3 * 1.d-6                         ! [MJ/mol]
        case(LIQUID_STATE) !---------------------------------------------------
          pressure_liq = gen_auxvar(ZERO_INTEGER,ghosted_id)% &
                         pres(option%liquid_phase)
          call EOSWaterEnthalpy(temperature,pressure_liq,H_liq,ierr)
          brine_energy = &
              cur_waste_panel%brine_generation_rate(i) * &  ! [mol/m3/sec]
              material_auxvars(ghosted_id)%volume * &       ! [m3-bulk] 
              H_liq * 1.d-3 * 1.d-6                         ! [MJ/mol]
        case(TWO_PHASE_STATE) !------------------------------------------------
          pressure_liq = gen_auxvar(ZERO_INTEGER,ghosted_id)% &
                         pres(option%liquid_phase)
          pressure_gas = gen_auxvar(ZERO_INTEGER,ghosted_id)% &
                         pres(option%gas_phase)
          call EOSWaterEnthalpy(temperature,pressure_liq,H_liq,ierr)
          call EOSGasEnergy(temperature,pressure_gas,H_gas,U_gas,ierr)
          brine_energy = &
              cur_waste_panel%brine_generation_rate(i) * &  ! [mol/m3/sec]
              material_auxvars(ghosted_id)%volume * &       ! [m3-bulk] 
              H_liq * 1.d-3 * 1.d-6                         ! [MJ/mol]
          gas_energy = &
              cur_waste_panel%gas_generation_rate(i) * &    ! [mol/m3/sec]
              material_auxvars(ghosted_id)%volume * &       ! [m3-bulk] 
              H_gas * 1.d-3 * 1.d-6                         ! [MJ/mol]
      end select
      vec_p(j) = brine_energy + gas_energy
    enddo
    !-------------------------------------------------------------------------
    cur_waste_panel => cur_waste_panel%next
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
  
  PetscReal :: rxnrate
  PetscInt :: cell_num
  type(chem_species_type) :: limiting_species
  PetscReal :: alpharxn
  
  PetscReal :: conc_ratio
  
  conc_ratio = ( limiting_species%initial_conc_mol(cell_num) / &
                 limiting_species%current_conc_mol(cell_num) ) 
  ! K_smoothed = K * (1.0 - exp(A*C/Ci)  BRAGFLO User's Manual Eq. 158
  rxnrate = rxnrate * (1.d0 - exp(alpharxn*conc_ratio))
  
end subroutine PMWSSSmoothRxnrate

! ************************************************************************** !

subroutine PMWSSTaperRxnrate1(rxnrate,cell_num,limiting_species1)
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
  
  PetscReal :: rxnrate
  PetscInt :: cell_num
  type(chem_species_type) :: limiting_species1
  
  if (limiting_species1%current_conc_mol(cell_num) <= 0.d0) then
    rxnrate = 0.d0
  else
    rxnrate = rxnrate
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
  
  PetscReal :: rxnrate
  PetscInt :: cell_num
  type(chem_species_type) :: limiting_species1
  type(chem_species_type) :: limiting_species2
  
  PetscReal :: rxnrate1
  PetscReal :: rxnrate2
  
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
  
  class(pm_wipp_srcsink_type) :: this
  
end subroutine PMWSSFinalizeTimestep

! ************************************************************************** !

subroutine PMWSSCalcParallelSUM(option,waste_panel,local_val,global_sum)
  ! 
  ! Calculates global sum for a MPI_DOUBLE_PRECISION number over a
  ! waste panel region. This function uses only MPI_Send and MPI_Recv functions
  ! and does not need a communicator object. It reduces communication to the
  ! processes that are in the waste panel's rank_list object rather than using
  ! a call to MPI_Allreduce.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/23/17
  
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(srcsink_panel_type) :: waste_panel
  PetscReal :: local_val
  PetscReal :: global_sum

  PetscReal, pointer :: temp_array(:)
  PetscInt :: num_ranks
  PetscInt :: m
  PetscInt :: TAG
  PetscErrorCode :: ierr
  
  num_ranks = size(waste_panel%rank_list)
  allocate(temp_array(num_ranks))
  temp_array = 0.d0
  TAG = 0
  
  if (num_ranks > 1) then
  !------------------------------------------
    if (option%myrank .ne. waste_panel%rank_list(1)) then
      call MPI_Send(local_val,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                    waste_panel%rank_list(1),TAG,option%mycomm,ierr)
    else
      temp_array(1) = local_val
      do m = 2,num_ranks
        call MPI_Recv(local_val,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                      waste_panel%rank_list(m),TAG,option%mycomm, &
                      MPI_STATUS_IGNORE,ierr)
        temp_array(m) = local_val
      enddo
      global_sum = sum(temp_array)
    endif
    if (option%myrank == waste_panel%rank_list(1)) then
      do m = 2,num_ranks
        call MPI_Send(global_sum,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                      waste_panel%rank_list(m),TAG,option%mycomm,ierr)
      enddo
    else
      call MPI_Recv(global_sum,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                    waste_panel%rank_list(1),TAG,option%mycomm, &
                    MPI_STATUS_IGNORE,ierr)
    endif             
  !------------------------------------------        
  else 
    global_sum = local_val
  endif
  
  deallocate(temp_array)

end subroutine PMWSSCalcParallelSUM

! *************************************************************************** !

 subroutine PMWSSOutput(this)
  ! 
  ! Sets up output for the process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/27/2017
  !
  
  use Option_module
  use Output_Aux_module
  
  implicit none

  class(pm_wipp_srcsink_type) :: this
  
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  class(srcsink_panel_type), pointer :: cur_waste_panel
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: i
  PetscReal :: local_gas_rate, global_gas_rate
  PetscReal :: local_brine_rate, global_brine_rate
  
  if (.not.associated(this%waste_panel_list)) return
  
100 format(100es18.8)
101 format(1I6.1)

  option => this%realization%option
  output_option => this%realization%output_option
  
  fid = 88
  filename = PMWSSOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")
       
  write(fid,100,advance="no") option%time / output_option%tconv
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    ! pre-calculations
    local_gas_rate = 0.d0 
    local_brine_rate = 0.d0
    do i = 1,cur_waste_panel%region%num_cells            
      local_gas_rate = local_gas_rate + &    ! [mol/m3-bulk/sec]
                   (cur_waste_panel%gas_generation_rate(i) * &
                    cur_waste_panel%scaling_factor(i))
      local_brine_rate = local_brine_rate + &    ! [mol/m3-bulk/sec]
                   (cur_waste_panel%brine_generation_rate(i) * &
                    cur_waste_panel%scaling_factor(i))
    enddo
    call PMWSSCalcParallelSUM(option,cur_waste_panel,local_gas_rate, &
                              global_gas_rate)
    call PMWSSCalcParallelSUM(option,cur_waste_panel,local_brine_rate, &
                              global_brine_rate)
    write(fid,101,advance="no") cur_waste_panel%id
    write(fid,100,advance="no") &
      cur_waste_panel%inventory%Fe_s%tot_mass_in_panel, &
      cur_waste_panel%inventory%MgO_s%tot_mass_in_panel, &
      cur_waste_panel%inventory%C6H10O5_s%tot_mass_in_panel, &
      cur_waste_panel%inventory%RuPl_s%tot_mass_in_panel, &
      global_gas_rate, &
      global_brine_rate
  !----------------------------------------
    cur_waste_panel => cur_waste_panel%next
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

  use Option_module

  implicit none
  
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: PMWSSOutputFilename
  character(len=MAXWORDLENGTH) :: word

  write(word,'(i6)') option%myrank
  PMWSSOutputFilename = trim(option%global_prefix) // &
                       trim(option%group_prefix) // &
                       '-' // trim(adjustl(word)) // '.pnl'
  
end function PMWSSOutputFilename  

! ************************************************************************** !

subroutine PMWSSOutputHeader(this)
  ! 
  ! Writes header for waste panel output file.
  ! 
  ! Author: Jennifer Frederick
  ! Date: 03/27/17

  use Output_Aux_module
  use Utility_module
    
  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  
  type(output_option_type), pointer :: output_option
  
  class(srcsink_panel_type), pointer :: cur_waste_panel
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: units_string 
  character(len=MAXWORDLENGTH) :: variable_string
  character(len=MAXWORDLENGTH) :: cell_string
  PetscInt :: fid
  PetscInt :: icolumn
  PetscBool :: exist
  
  if (.not.associated(this%waste_panel_list)) return
  
  output_option => this%realization%output_option
  
  fid = 88
  filename = PMWSSOutputFilename(this%option)
  exist = FileExists(trim(filename))
  if (this%option%restart_flag .and. exist) return
  open(unit=fid,file=filename,action="write",status="replace")
  
  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif
  
  write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    cell_string = trim(cur_waste_panel%region_name)
    variable_string = 'WP ID#'
    units_string = ''
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Fe(s) mass'
    units_string = 'kg'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MgO(s) mass'
    units_string = 'kg'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Cellulosics(s) mass'
    units_string = 'kg'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Rubber/Plastics(s) mass'
    units_string = 'kg'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Avg. gas gen. rate'
    units_string = 'mol/m3-bulk/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Avg. brine gen. rate'
    units_string = 'mol/m3-bulk/sec'
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
  
  class(pm_wipp_srcsink_type) :: this

  PetscInt :: id

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
  
  class(pm_wipp_srcsink_type) :: this
  
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
  
  class(pm_wipp_srcsink_type) :: this
  
  type(srcsink_panel_type), pointer :: cur_waste_panel, prev_waste_panel
  type(pre_inventory_type), pointer :: cur_preinventory, prev_preinventory

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
  
  type(srcsink_panel_type), pointer :: waste_panel
  
  if (.not.associated(waste_panel)) return
  call PMWSSInventoryDestroy(waste_panel%inventory)
  call DeallocateArray(waste_panel%scaling_factor)
  call DeallocateArray(waste_panel%gas_generation_rate)
  call DeallocateArray(waste_panel%brine_generation_rate)
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
  
  type(inventory_type) :: inventory
  
  call PMWSSChemSpeciesDeallocate(inventory%Fe_s)
  call PMWSSChemSpeciesDeallocate(inventory%FeOH2_s)
  call PMWSSChemSpeciesDeallocate(inventory%C6H10O5_s)
  call PMWSSChemSpeciesDeallocate(inventory%RuPl_s)
  call PMWSSChemSpeciesDeallocate(inventory%NO3_minus_aq)
  call PMWSSChemSpeciesDeallocate(inventory%CO2_g)
  call PMWSSChemSpeciesDeallocate(inventory%N2_g)
  call PMWSSChemSpeciesDeallocate(inventory%SO42_minus_aq)
  call PMWSSChemSpeciesDeallocate(inventory%H2S_g)
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
  
  type(chem_species_type) :: chem_species
  
  call DeallocateArray(chem_species%initial_conc_mol)
  call DeallocateArray(chem_species%current_conc_mol)
  call DeallocateArray(chem_species%current_conc_kg)
  call DeallocateArray(chem_species%inst_rate)
  
end subroutine PMWSSChemSpeciesDeallocate

! *************************************************************************** !

end module PM_WIPP_SrcSink_class
