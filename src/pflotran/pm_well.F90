module PM_Well_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use Option_module
  use Geometry_module
  use Strata_module
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module
  use WIPP_Flow_Aux_module
  use NW_Transport_Aux_module
  use Well_Grid_module

  implicit none

  private

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !    TOP
  !   -------- ql_bc(2)
  !   | i=n|
  !   -------- liq%ql(k=n-1)
  !   |    |
  !   ------
  !   |    |
  !   -------- liq%ql(k=3)
  !   | i=3|
  !   -------- liq%ql(k=2)
  !   | i=2|
  !   -------- liq%ql(k=1)
  !   | i=1|
  !   -------- ql_bc(1)
  !    BOT
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  PetscBool :: initialize_well_flow = PETSC_TRUE
  PetscBool :: initialize_well_tran = PETSC_TRUE
  PetscReal :: min_flow_dt_scale = 1.d-3
  PetscReal, parameter :: gravity = -9.80665d0 ! [m/s2]

  PetscInt, parameter :: PEACEMAN_ISO = 1
  PetscInt, parameter :: PEACEMAN_2D = 2
  PetscInt, parameter :: PEACEMAN_3D = 3

  ! srcsink vector indexing
  PetscInt, parameter :: UNPERT = 0
  PetscInt, parameter :: PERT_WRT_PL = 1
  PetscInt, parameter :: PERT_WRT_SG = 2

  ! Well model types
  PetscInt, parameter, public :: WELL_MODEL_HYDROSTATIC = 1
  PetscInt, parameter, public :: WELL_MODEL_WIPP_DARCY = 2
  PetscInt, parameter, public :: WELL_MODEL_FULL_MOMENTUM = 3

  ! Well constraint types
  PetscInt, parameter :: CONSTANT_PRESSURE = 1
  PetscInt, parameter :: CONSTANT_PRESSURE_HYDROSTATIC = 2
  PetscInt, parameter :: CONSTANT_RATE = 3

  ! Well coupling options
  PetscInt, parameter, public :: FULLY_IMPLICIT_WELL = ONE_INTEGER
  PetscInt, parameter, public :: QUASI_IMPLICIT_WELL = TWO_INTEGER
  PetscInt, parameter, public :: SEQUENTIAL_WELL = THREE_INTEGER


  type :: well_reservoir_type
    ! reservoir liquid pressure [Pa]
    PetscReal, pointer :: p_l(:)
    ! reservoir gas pressure [Pa]
    PetscReal, pointer :: p_g(:)
    ! reservoir liquid saturation [-]
    PetscReal, pointer :: s_l(:)
    ! reservoir gas saturation [-]
    PetscReal, pointer :: s_g(:)
    ! reservoir temperature [C]
    PetscReal, pointer :: temp(:)
    ! reservoir liquid mobility
    PetscReal, pointer :: mobility_l(:)
    ! reservoir gas mobility
    PetscReal, pointer :: mobility_g(:)
    ! reservoir liquid relative permeability [-]
    PetscReal, pointer :: kr_l(:)
    ! reservoir gas relative permeability [-]
    PetscReal, pointer :: kr_g(:)
    ! reservoir liquid density [kg/m3]
    PetscReal, pointer :: den_l(:)
    ! reservoir gas density [kg/m3]
    PetscReal, pointer :: den_g(:)
    ! reservoir liquid dynamic viscosity [Pa-s]
    PetscReal, pointer :: visc_l(:)
    ! reservoir gas dynamic viscosity [Pa-s]
    PetscReal, pointer :: visc_g(:)
    ! reservoir effective porosity (factors in compressibility) [m3/m3]
    PetscReal, pointer :: e_por(:)
    ! reservoir species aqueous concentration [mol/m3-liq] (idof,conc@segment)
    PetscReal, pointer :: aqueous_conc(:,:)
    ! reservoir species aqueous mass [mol] (idof,mass@segment)
    PetscReal, pointer :: aqueous_mass(:,:)
    ! reservior fluid compositions (flow)
    PetscReal, pointer :: xmass_liq(:,:)
    PetscReal, pointer :: xmass_gas(:,:)
    ! reservoir permeabilities [m2]
    PetscReal, pointer :: kx(:)
    PetscReal, pointer :: ky(:)
    PetscReal, pointer :: kz(:)
    ! reservoir discretization [m]
    PetscReal, pointer :: dx(:)
    PetscReal, pointer :: dy(:)
    PetscReal, pointer :: dz(:)
    ! reservoir cell volume [m3]
    PetscReal, pointer :: volume(:)
  end type

  type :: well_type
    ! type of well model to use (flow mode)
    PetscInt :: well_model_type
    ! type of well model constraint (pressure control, rate control)
    PetscInt :: well_constraint_type
    ! reservoir variables associated with a particular well or well pert
    type(well_reservoir_type), pointer :: reservoir
    type(well_reservoir_type), pointer :: reservoir_save
    ! mass balance of liquid [kmol/s????????]
    PetscReal, pointer :: mass_balance_liq(:)
    ! well fluid properties
    type(well_fluid_type), pointer :: liq
    type(well_fluid_type), pointer :: gas
    ! cross-sectional area of each well segment [calc'd] [m2]
    PetscReal, pointer :: area(:)
    ! diameter of each well segment [m]
    PetscReal, pointer :: diameter(:)
    ! volume of each well segment [calc'd] [m3]
    PetscReal, pointer :: volume(:)
    ! friction ceofficient of each well segment
    PetscReal, pointer :: friction_factor(:)
    ! skin factor of each well segment
    PetscReal, pointer :: skin(:)
    ! well index of each well segment [0,1]  0 = cased; 1 = open
    PetscReal, pointer :: WI_base(:)
    ! total well index for each well segment (including reservoir effects)
    PetscReal, pointer :: WI(:)
    ! Peaceman equivalent radius
    PetscReal, pointer :: r0(:)
    ! well index model (probably has to get moved out of well_type)
    PetscInt :: WI_model = PEACEMAN_3D
    ! well liquid pressure [Pa]
    PetscReal, pointer :: pl(:)
    ! well gas pressure [Pa]
    PetscReal, pointer :: pg(:)
    ! well temperature [C]
    PetscReal, pointer :: temp(:)
    ! flag for output
    PetscBool :: output_pl 
    ! flag for output
    PetscBool :: output_pg
    ! flag for output
    PetscBool :: output_sl
    ! flag for output
    PetscBool :: output_sg 
    ! well liquid Darcy flux [m3/m2-s] of interior interfaces
    PetscReal, pointer :: ql(:)
    ! well gas Darcy flux [m3/m2-s] of interior interfaces
    PetscReal, pointer :: qg(:)
    ! well liquid Darcy flux [m3/m2-s] of top/bottom boundaries
    PetscReal, pointer :: ql_bc(:)
    ! well gas Darcy flux [m3/m2-s] of top/bottom boundaries
    PetscReal, pointer :: qg_bc(:)
    ! well liquid mass flux [kmol/s] of interior interfaces
    PetscReal, pointer :: ql_kmol(:)
    ! well gas mass flux [kmol/s] of interior interfaces
    PetscReal, pointer :: qg_kmol(:)
    ! well liquid mass flux [kmol/s] of top/bottom boundaries
    PetscReal, pointer :: ql_kmol_bc(:)
    ! well gas mass flux [kmol/s] of top/bottom boundaries
    PetscReal, pointer :: qg_kmol_bc(:)
    ! well liquid cumulative mass [kmol] in each segment
    PetscReal, pointer :: liq_cum_mass(:)
    ! well liquid instantaneous mass [kmol] in each segment
    PetscReal, pointer :: liq_mass(:)
    ! well species names
    character(len=MAXWORDLENGTH), pointer :: species_names(:)
    ! well species parent id number
    PetscInt, pointer :: species_parent_id(:)
    ! well species radioactive flag
    PetscBool, pointer :: species_radioactive(:)
    ! well species decay rate [1/sec]
    PetscReal, pointer :: species_decay_rate(:)
    ! well species' parent decay rate [1/sec]
    PetscReal, pointer :: species_parent_decay_rate(:)
    ! well species aqueous concentration [mol/m3-liq] (ispecies,conc@segment)
    PetscReal, pointer :: aqueous_conc(:,:)
    ! well species aqueous mass [mol] (ispecies,mass@segment)
    PetscReal, pointer :: aqueous_mass(:,:)
    ! flag for output
    PetscBool :: output_aqc 
    ! flag for output
    PetscBool :: output_aqm
    ! well species aqueous concentration top of hole BC [mol/m3-liq]
    PetscReal, pointer :: aqueous_conc_th(:)
    ! well bottom of hole pressure BC flag
    PetscBool :: bh_p_set_by_reservoir
    ! well bottom of hole Sg BC flag
    PetscBool :: bh_sg_set_by_reservoir
    ! well bottom of hole pressure BC [Pa]
    PetscReal :: bh_p
    ! well top of hole pressure BC [Pa]
    PetscReal :: th_p
    ! well bottom of hole saturation BC [Pa]
    PetscReal :: bh_sg
    ! well top of hole saturation BC [Pa]
    PetscReal :: th_sg
    ! well bottom of hole rate BC [kg/s]
    PetscReal :: bh_ql
    PetscReal :: bh_qg
    ! well top of hole rate BC [kg/s]
    PetscReal :: th_ql
    PetscReal :: th_qg
    ! well transport constraint name
    character(len=MAXWORDLENGTH) :: tran_condition_name
    ! Link to characteristic curves
    PetscInt, pointer :: ccid(:)
    ! permeability along the well [m2]
    PetscReal, pointer :: permeability(:)
    ! porosity
    PetscReal, pointer :: phi(:)
  end type well_type

  type :: well_fluid_type
    ! fluid phase ID (liq/gas)
    PetscInt :: ifluid
    ! fluid rel perm
    PetscReal, pointer :: kr(:)
    ! fluid reference density [kg/m3]
    PetscReal :: den0
    ! fluid density [kg/m3]
    PetscReal, pointer :: den(:)
    ! fluid viscosity [Pa-s]
    PetscReal, pointer :: visc(:)
    ! fluid saturation
    PetscReal, pointer :: s(:)
    ! fluid composition (flow)
    PetscReal, pointer :: xmass(:,:)
    ! fluid source/sink in/out of well [kmol/s]
    PetscReal, pointer :: Q(:)
    ! flag for output
    PetscBool :: output_Q
  end type well_fluid_type

  ! primary variables necessary to reset flow solution
  type :: well_flow_save_type
    PetscReal, pointer :: pl(:)
    PetscReal, pointer :: sg(:)
    PetscReal :: bh_p
  end type well_flow_save_type

  ! primary variables necessary to reset transport solution
  type :: well_tran_save_type
    PetscReal, pointer :: aqueous_conc(:,:)
    PetscReal, pointer :: aqueous_mass(:,:)
  end type well_tran_save_type

  type :: well_soln_base_type
    ! number of primary variables
    PetscInt :: ndof
    ! residual vector
    PetscReal, pointer :: residual(:)
    ! Jacobian matrix
    PetscReal, pointer :: Jacobian(:,:)
    ! Solution update
    PetscReal, pointer :: update(:)
    ! convergence flags
    PetscBool :: not_converged
    PetscBool :: converged
    PetscBool :: cut_timestep
    ! maximum number of iterations
    PetscInt :: max_iter
    ! maximum number of time step cuts
    PetscInt :: max_ts_cut
    ! time step cut factor (divides the time step)
    PetscReal :: ts_cut_factor
    ! time step ramp-up factor (multiplies the time step)
    PetscReal :: ts_ramp_factor
    ! convergence criterial
    PetscReal :: itol_abs_res
    PetscReal :: itol_scaled_res
    ! solver statistics
    PetscInt :: n_steps
    PetscInt :: n_newton
  end type well_soln_base_type

  type, extends(well_soln_base_type) :: well_soln_flow_type
    ! Previously converged flow solution
    type(well_flow_save_type) :: prev_soln
    ! flow solution at the end of a WIPP_FLOW timestep
    type(well_flow_save_type) :: soln_save
    ! flag for bottom of hole pressure BC
    PetscBool :: bh_p
    ! flag for top of hole pressure BC
    PetscBool :: th_p
    ! flag for bottom of hole rate BC
    PetscBool :: bh_q
    ! flag for top of hole rate BC
    PetscBool :: th_q
    ! convergence criterial
    PetscReal :: itol_abs_update_p
    PetscReal :: itol_abs_update_s
    PetscReal :: itol_rel_update_p
    PetscReal :: itol_rel_update_s
  end type well_soln_flow_type

  type, extends(well_soln_base_type) :: well_soln_tran_type
    ! Previously converged transport solution
    type(well_tran_save_type) :: prev_soln
    ! convergence criterial
    PetscReal :: itol_abs_update
    PetscReal :: itol_rel_update
    ! time tracking
    PetscReal :: tran_time
    ! time step cut flag for quasi-implicit coupling
    PetscBool :: cut_ts_flag
  end type well_soln_tran_type

  type :: well_comm_type
    ! rank in PETSC_COMM_WORLD (from option%myrank)
    PetscMPIInt :: petsc_rank 
    PetscMPIInt, pointer :: petsc_rank_list(:) 
    PetscMPIInt, pointer :: well_rank_list(:) 
    ! WELL_COMM_WORLD
    PetscMPIInt :: comm      
    ! group in WELL_COMM_WORLD
    PetscMPIInt :: group      
    ! rank in WELL_COMM_WORLD
    PetscMPIInt :: rank        
    ! size of WELL_COMM_WORLD
    PetscMPIInt :: commsize
  end type well_comm_type

  type, public, extends(pm_base_type) :: pm_well_type
    class(realization_subsurface_type), pointer :: realization
    type(well_grid_type), pointer :: well_grid
    type(well_type), pointer :: well
    type(well_type), pointer :: well_pert(:)
    type(well_soln_flow_type), pointer :: flow_soln
    type(well_soln_tran_type), pointer :: tran_soln
    type(strata_list_type), pointer :: strata_list
    type(well_comm_type), pointer :: well_comm
    PetscReal, pointer :: pert(:,:)
    PetscReal, pointer :: srcsink_water(:,:)
    PetscReal, pointer :: srcsink_gas(:,:)
    PetscInt :: nphase
    PetscInt :: nspecies
    PetscReal :: intrusion_time_start           ! [sec]
    PetscReal :: bh_zero_value                  ! [mol/m3-bulk]
    PetscReal :: dt_flow, dt_tran               ! [sec]
    PetscReal :: min_dt_flow, min_dt_tran       ! [sec]
    PetscReal :: cumulative_dt_flow             ! [sec]
    PetscReal :: cumulative_dt_tran             ! [sec]
    PetscBool :: transport
    PetscBool :: ss_check
    PetscBool :: well_on    !Turns the well on, regardless of other checks
    PetscInt :: well_force_ts_cut
    PetscBool :: update_for_wippflo_qi_coupling
    PetscInt :: flow_coupling
    PetscBool :: print_well
    PetscBool :: print_output
  contains
    procedure, public :: Setup => PMWellSetup
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMWellReadSimOptionsBlock
    procedure, public :: ReadPMBlock => PMWellReadPMBlock
    procedure, public :: SetRealization => PMWellSetRealization
    procedure, public :: InitializeRun => PMWellInitializeRun
    procedure, public :: FinalizeRun => PMWellFinalizeRun
    procedure, public :: InitializeTimestep => PMWellInitializeTimestep
    procedure, public :: UpdateTimestep => PMWellUpdateTimestep
    procedure, public :: FinalizeTimestep => PMWellFinalizeTimestep
    procedure, public :: PreSolve => PMWellPreSolve
    procedure, public :: Solve => PMWellSolve
    procedure, public :: PostSolve => PMWellPostSolve
    procedure, public :: InputRecord => PMWellInputRecord
    procedure, public :: Destroy => PMWellDestroy
    procedure, public :: Perturb => PMWellBasePerturb
    procedure, public :: ModifyFlowJacobian => PMWellModifyFlowJacobian
    procedure, public :: ModifyFlowResidual => PMWellModifyFlowResidual
    procedure, public :: UpdateFlowRates => PMWellUpdateFlowRates
  end type pm_well_type

  public :: PMWellCreate, &
            PMWellSetupGrid, &
            PMWellReadPMBlock, &
            PMWellReadGrid, &
            PMWellReadPass2, &
            PMWellQISolveTran, &
            PMWellUpdateReservoirSrcSinkFlow

  contains

! ************************************************************************** !

function PMWellCreate()
  !
  ! Creates the well process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type), pointer :: PMWellCreate
  class(pm_well_type), pointer :: this

  allocate(this)
  call PMBaseInit(this)

  this%header = 'WELLBORE MODEL'

  nullify(this%realization)
  this%well_grid => WellGridCreate()
  allocate(this%well)
  call PMWellVarCreate(this%well)
  call PMWellFlowCreate(this)
  call PMWellTranCreate(this)

  ! strata list specific to well
  allocate(this%strata_list)
  nullify(this%strata_list%first)
  nullify(this%strata_list%last)
  nullify(this%strata_list%array)
  this%strata_list%num_strata = 0

  allocate(this%well_comm)
  nullify(this%well_comm%petsc_rank_list)
  nullify(this%well_comm%well_rank_list)
  this%well_comm%petsc_rank = 0
  this%well_comm%comm = 0
  this%well_comm%group = 0
  this%well_comm%rank = 0
  this%well_comm%commsize = 0

  nullify(this%pert)
  nullify(this%srcsink_water)
  nullify(this%srcsink_gas)

  this%nphase = 0
  this%nspecies = 0
  this%intrusion_time_start = UNINITIALIZED_DOUBLE
  this%bh_zero_value = 1.d-40
  this%dt_flow = 0.d0
  this%dt_tran = 0.d0
  this%min_dt_flow = 1.d-15
  this%min_dt_tran = 1.d-15
  this%cumulative_dt_flow = 0.d0
  this%cumulative_dt_tran = 0.d0
  this%transport = PETSC_FALSE
  this%ss_check = PETSC_FALSE
  this%well_on = PETSC_TRUE
  this%well_force_ts_cut = 0
  this%update_for_wippflo_qi_coupling = PETSC_FALSE
  this%flow_coupling = UNINITIALIZED_INTEGER
  this%print_well = PETSC_FALSE
  this%print_output = PETSC_FALSE

  PMWellCreate => this

end function PMWellCreate

! ************************************************************************** !

subroutine PMWellFlowCreate(pm_well)
  !
  ! Creates flow variables.
  !
  ! Author: Michael Nole
  ! Date: 03/04/2023

  implicit none

  class(pm_well_type), pointer :: pm_well

  ! create the well flow solution object:
  allocate(pm_well%flow_soln)

  pm_well%flow_soln%ndof = UNINITIALIZED_INTEGER

  nullify(pm_well%flow_soln%residual)
  nullify(pm_well%flow_soln%Jacobian)
  nullify(pm_well%flow_soln%update)

  pm_well%flow_soln%not_converged = PETSC_TRUE
  pm_well%flow_soln%converged = PETSC_FALSE
  pm_well%flow_soln%cut_timestep = PETSC_FALSE

  pm_well%flow_soln%max_iter = 8
  pm_well%flow_soln%max_ts_cut = 20
  pm_well%flow_soln%ts_cut_factor = 2.0d0
  pm_well%flow_soln%ts_ramp_factor = 1.1d0
  pm_well%flow_soln%itol_abs_res = 1.0d-8
  pm_well%flow_soln%itol_scaled_res = 1.0d-5
  pm_well%flow_soln%n_steps = 0
  pm_well%flow_soln%n_newton = 0

  nullify(pm_well%flow_soln%prev_soln%pl)
  nullify(pm_well%flow_soln%prev_soln%sg)
  nullify(pm_well%flow_soln%soln_save%pl)
  nullify(pm_well%flow_soln%soln_save%sg)
  pm_well%flow_soln%prev_soln%bh_p = UNINITIALIZED_DOUBLE
  pm_well%flow_soln%soln_save%bh_p = UNINITIALIZED_DOUBLE

  pm_well%flow_soln%bh_p = PETSC_FALSE
  pm_well%flow_soln%th_p = PETSC_FALSE
  pm_well%flow_soln%bh_q = PETSC_FALSE
  pm_well%flow_soln%th_q = PETSC_FALSE
  pm_well%flow_soln%itol_abs_update_p = 1.0d0
  pm_well%flow_soln%itol_abs_update_s = 1.0d-5
  pm_well%flow_soln%itol_rel_update_p = 1.0d-4
  pm_well%flow_soln%itol_rel_update_s = 1.0d-4

end subroutine PMWellFlowCreate

! ************************************************************************** !

subroutine PMWellTranCreate(pm_well)
  !
  ! Creates transport variables.
  !
  ! Author: Michael Nole
  ! Date: 03/04/2023

  implicit none

  class(pm_well_type), pointer :: pm_well

  ! create the well transport solution object:
  allocate(pm_well%tran_soln)
  nullify(pm_well%tran_soln%prev_soln%aqueous_conc)
  nullify(pm_well%tran_soln%prev_soln%aqueous_mass)

  pm_well%tran_soln%ndof = UNINITIALIZED_INTEGER
  nullify(pm_well%tran_soln%residual)
  nullify(pm_well%tran_soln%Jacobian)
  nullify(pm_well%tran_soln%update)

  pm_well%tran_soln%not_converged = PETSC_TRUE
  pm_well%tran_soln%converged = PETSC_FALSE
  pm_well%tran_soln%cut_timestep = PETSC_FALSE
  pm_well%tran_soln%max_iter = 8
  pm_well%tran_soln%max_ts_cut = 20
  pm_well%tran_soln%ts_cut_factor = 2.0d0
  pm_well%tran_soln%ts_ramp_factor = 1.1d0
  pm_well%tran_soln%itol_abs_res = 1.0d-8
  pm_well%tran_soln%itol_scaled_res = 1.0d-4
  pm_well%tran_soln%n_steps = 0
  pm_well%tran_soln%n_newton = 0
  pm_well%tran_soln%cut_ts_flag = PETSC_FALSE

  pm_well%tran_soln%itol_abs_update = 1.0d0
  pm_well%tran_soln%itol_rel_update = 1.0d-1

end subroutine PMWellTranCreate

! ************************************************************************** !

subroutine PMWellResCreate(reservoir)
  !
  ! Creates reservoir variables.
  !
  ! Author: Michael Nole
  ! Date: 03/04/2023

  implicit none

  type(well_reservoir_type), pointer :: reservoir

  ! create the reservoir object:
  allocate(reservoir)
  nullify(reservoir%p_l)
  nullify(reservoir%p_g)
  nullify(reservoir%s_l)
  nullify(reservoir%s_g)
  nullify(reservoir%temp)
  nullify(reservoir%mobility_l)
  nullify(reservoir%mobility_g)
  nullify(reservoir%kr_l)
  nullify(reservoir%kr_g)
  nullify(reservoir%den_l)
  nullify(reservoir%den_g)
  nullify(reservoir%visc_l)
  nullify(reservoir%visc_g)
  nullify(reservoir%e_por)
  nullify(reservoir%aqueous_conc)
  nullify(reservoir%aqueous_mass)
  nullify(reservoir%xmass_liq)
  nullify(reservoir%xmass_gas)
  nullify(reservoir%kx)
  nullify(reservoir%ky)
  nullify(reservoir%kz)
  nullify(reservoir%dx)
  nullify(reservoir%dy)
  nullify(reservoir%dz)
  nullify(reservoir%volume)

end subroutine PMWellResCreate

! ************************************************************************** !

subroutine PMWellVarCreate(well)
  !
  ! Creates well variables.
  !
  ! Author: Michael Nole
  ! Date: 03/04/2023

  implicit none

  type(well_type) :: well

  ! create the well_pert object:
  well%well_model_type = UNINITIALIZED_INTEGER
  well%well_constraint_type = UNINITIALIZED_INTEGER
  nullify(well%mass_balance_liq)
  nullify(well%liq)
  nullify(well%gas)
  nullify(well%area)
  nullify(well%diameter)
  nullify(well%volume)
  nullify(well%friction_factor)
  nullify(well%skin)
  nullify(well%WI_base)
  nullify(well%WI)
  nullify(well%r0)
  well%WI_model = PEACEMAN_3D
  nullify(well%pl)
  nullify(well%pg)
  nullify(well%temp)
  well%output_pl = PETSC_FALSE
  well%output_pg = PETSC_FALSE
  well%output_sl = PETSC_FALSE
  well%output_sg = PETSC_FALSE
  nullify(well%ql)
  nullify(well%qg)
  nullify(well%ql_bc)
  nullify(well%qg_bc)
  nullify(well%ql_kmol)
  nullify(well%qg_kmol)
  nullify(well%ql_kmol_bc)
  nullify(well%qg_kmol_bc)
  nullify(well%liq_cum_mass)
  nullify(well%liq_mass)
  nullify(well%species_names)
  nullify(well%species_parent_id)
  nullify(well%species_radioactive)
  nullify(well%species_decay_rate)
  nullify(well%species_parent_decay_rate)
  nullify(well%aqueous_conc)
  nullify(well%aqueous_mass)
  well%output_aqc = PETSC_FALSE
  well%output_aqm = PETSC_FALSE
  nullify(well%aqueous_conc_th)
  well%bh_p_set_by_reservoir = PETSC_FALSE
  well%bh_sg_set_by_reservoir = PETSC_FALSE
  well%bh_p = UNINITIALIZED_DOUBLE
  well%th_p = UNINITIALIZED_DOUBLE
  well%bh_sg = UNINITIALIZED_DOUBLE
  well%th_sg = UNINITIALIZED_DOUBLE
  well%bh_ql = UNINITIALIZED_DOUBLE
  well%bh_qg = UNINITIALIZED_DOUBLE
  well%th_ql = UNINITIALIZED_DOUBLE
  well%th_qg = UNINITIALIZED_DOUBLE
  well%tran_condition_name = ''
  nullify(well%ccid)
  nullify(well%permeability)
  nullify(well%phi)

  ! create the fluid/liquid objects:
  allocate(well%liq)
  well%liq%ifluid = 1
  nullify(well%liq%kr)
  well%liq%den0 = UNINITIALIZED_DOUBLE
  nullify(well%liq%den)
  nullify(well%liq%visc)
  nullify(well%liq%s)
  nullify(well%liq%xmass)
  nullify(well%liq%Q)
  well%liq%output_Q = PETSC_FALSE

  ! create the fluid/gas objects:
  allocate(well%gas)
  well%gas%ifluid = 2
  nullify(well%gas%kr)
  well%gas%den0 = UNINITIALIZED_DOUBLE
  nullify(well%gas%den)
  nullify(well%gas%visc)
  nullify(well%gas%s)
  nullify(well%gas%xmass)
  nullify(well%gas%Q)
  well%gas%output_Q = PETSC_FALSE

  ! create the reservoir object
  allocate(well%reservoir)
  nullify(well%reservoir%p_l)
  nullify(well%reservoir%p_g)
  nullify(well%reservoir%s_l)
  nullify(well%reservoir%s_g)
  nullify(well%reservoir%temp)
  nullify(well%reservoir%mobility_l)
  nullify(well%reservoir%mobility_g)
  nullify(well%reservoir%kr_l)
  nullify(well%reservoir%kr_g)
  nullify(well%reservoir%den_l)
  nullify(well%reservoir%den_g)
  nullify(well%reservoir%visc_l)
  nullify(well%reservoir%visc_g)
  nullify(well%reservoir%e_por)
  nullify(well%reservoir%aqueous_conc)
  nullify(well%reservoir%aqueous_mass)
  nullify(well%reservoir%xmass_liq)
  nullify(well%reservoir%xmass_gas)
  nullify(well%reservoir%kx)
  nullify(well%reservoir%ky)
  nullify(well%reservoir%kz)
  nullify(well%reservoir%dx)
  nullify(well%reservoir%dy)
  nullify(well%reservoir%dz)
  nullify(well%reservoir%volume)

  ! create the reservoir_save object
  allocate(well%reservoir_save)
  nullify(well%reservoir_save%p_l)
  nullify(well%reservoir_save%p_g)
  nullify(well%reservoir_save%s_l)
  nullify(well%reservoir_save%s_g)
  nullify(well%reservoir_save%temp)
  nullify(well%reservoir_save%mobility_l)
  nullify(well%reservoir_save%mobility_g)
  nullify(well%reservoir_save%kr_l)
  nullify(well%reservoir_save%kr_g)
  nullify(well%reservoir_save%den_l)
  nullify(well%reservoir_save%den_g)
  nullify(well%reservoir_save%visc_l)
  nullify(well%reservoir_save%visc_g)
  nullify(well%reservoir_save%e_por)
  nullify(well%reservoir_save%aqueous_conc)
  nullify(well%reservoir_save%aqueous_mass)
  nullify(well%reservoir_save%xmass_liq)
  nullify(well%reservoir_save%xmass_gas)
  nullify(well%reservoir_save%kx)
  nullify(well%reservoir_save%ky)
  nullify(well%reservoir_save%kz)
  nullify(well%reservoir_save%dx)
  nullify(well%reservoir_save%dy)
  nullify(well%reservoir_save%dz)
  nullify(well%reservoir_save%volume)


end subroutine PMWellVarCreate

! ************************************************************************** !

subroutine PMWellSetupGrid(well_grid,res_grid,option)
  !
  ! Sets up a well grid based off of well grid info and reservoir info.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  !

  use Grid_module
  use Utility_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Input_Aux_module
  use Dataset_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Transport_Constraint_NWT_module
  use NW_Transport_Aux_module
  use SCO2_Aux_module

  implicit none

  type(well_grid_type), pointer :: well_grid
  type(grid_type), pointer :: res_grid
  type(option_type), pointer :: option
  
  type(point3d_type) :: dummy_h
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscReal :: diff_x,diff_y,diff_z
  PetscReal :: dh_x,dh_y,dh_z
  PetscReal :: total_length
  PetscReal :: top_of_reservoir, top_of_hole
  PetscReal :: bottom_of_reservoir, bottom_of_hole
  PetscReal :: temp_real, temp_real2
  PetscReal :: cum_z, z_dum
  PetscInt, allocatable :: temp_id_list(:)
  PetscReal, allocatable :: temp_z_list(:)
  PetscInt, allocatable :: temp_repeated_list(:)
  PetscInt :: cur_id, cum_z_int, cur_cum_z_int 
  PetscInt :: repeated
  PetscInt :: num_entries
  PetscReal, allocatable :: dz_list(:), res_dz_list(:)
  PetscReal :: min_dz, dz, z
  PetscInt, allocatable :: cell_id_list(:)
  PetscInt :: local_id
  PetscInt :: res_cell_count
  PetscInt :: nsegments, nsegments_save
  PetscInt :: k, i, j

  num_entries = 10000
  allocate(dz_list(num_entries))
  allocate(res_dz_list(num_entries))
  allocate(cell_id_list(num_entries))

  top_of_reservoir = res_grid%z_max_global
  top_of_hole = well_grid%tophole(3)
  bottom_of_reservoir = res_grid%z_min_global
  bottom_of_hole = well_grid%bottomhole(3)
  if (top_of_reservoir < top_of_hole) then
    option%io_buffer = 'The WELLBORE_MODEL TOP_OF_HOLE coordinates extend &
                       &beyond the top of the reservoir domain. &
                       &You must fix the TOP_OF_HOLE coordinates to align &
                       &with the top face of the reservoir grid cell that &
                       &it occupies.'
    call PrintErrMsg(option)
  endif
  if (top_of_reservoir > top_of_hole) then
    option%io_buffer = 'The WELLBORE_MODEL TOP_OF_HOLE coordinates do not &
                       &reach the top of the reservoir domain. &
                       &You must fix the TOP_OF_HOLE coordinates to align &
                       &with the top face of the reservoir grid cell that &
                       &it occupies.'
    call PrintErrMsg(option)
  endif
  if (bottom_of_reservoir > bottom_of_hole) then
    option%io_buffer = 'The WELLBORE_MODEL BOTTOM_OF_HOLE coordinates extend &
                       &beyond the bottom of the reservoir domain. &
                       &You must fix the BOTTOM_OF_HOLE coordinates so that &
                       &the bottom of the well is aligned with the bottom &
                       &face of the reservoir, or is above the bottom &
                       &face of the reservoir in the vertical column that the &
                       &well occupies.'
    call PrintErrMsg(option)
  endif

  diff_x = well_grid%tophole(1)-well_grid%bottomhole(1)
  diff_y = well_grid%tophole(2)-well_grid%bottomhole(2)
  diff_z = well_grid%tophole(3)-well_grid%bottomhole(3)

  if ((diff_y >= 1.d-10) .or. (diff_x >= 1.d-10)) then
    option%io_buffer = 'WELLBORE_MODEL does not support a tilted &
                        &well geometry. Please ensure that the well is &
                        &perfectly vertical, and the vertical direction is &
                        &set to the z-axis. Tilted well geometry is &
                        &still in development.'
    call PrintErrMsg(option)
  endif

  if (well_grid%match_reservoir) then
  !---------------------------------------------------------------------------
    if (option%driver%comm%size > 1) then
      option%io_buffer = 'WELLBORE_MODEL WELL_GRID,MATCH_RESERVOIR option &
        &is not supported yet for parallel simulations. Use & 
        &WELL_GRID,NUMBER_OF_SEGMENTS option or WELL_GRID,&
        &SEGMENT_CENTER_Z_VALUES with SEGMENT_LENGTH_VALUES option to define &
        &the wellbore model grid.'
      call PrintErrMsg(option)
    endif

    allocate(temp_id_list(10000))
    allocate(temp_repeated_list(10000))
    temp_id_list = -999
    temp_repeated_list = -999
    dummy_h%x = well_grid%bottomhole(1)
    dummy_h%y = well_grid%bottomhole(2) 

    k = 0; j = 0; cur_id = -999; repeated = 0; cum_z = 0; cur_cum_z_int = 0; 
    ! search and peck procedure for finding reservoir z list within well
    do
      j = j + 1
      cum_z = (j)*well_grid%dz_peck
      cum_z_int = int(cum_z)
      dummy_h%z = well_grid%bottomhole(3) + cum_z
      if (dummy_h%z > well_grid%tophole(3)) exit

      call GridGetLocalIDFromCoordinate(res_grid,dummy_h,option,local_id)
      if (cum_z_int > cur_cum_z_int) then
        call PrintProgressBarInt(diff_z,5,cum_z_int)
        cur_cum_z_int = cum_z_int 
      endif
      
      if (j == 1) then
        cur_id = local_id
        k = 1
        temp_id_list(k) = local_id
      endif 

      if (local_id /= cur_id) then
        k = k + 1
        if (k > 10000) then
          option%io_buffer = 'More than 10,000 reservoir cells have been &
                             &counted in the z-direction within the wellbore.'
          call PrintErrMsgToDev(option, &
                           'if reducing to less than 10,000 is not an option.')
        endif
        temp_id_list(k) = local_id 
        temp_repeated_list(k-1) = repeated
        repeated = 0 
        cur_id = local_id
      endif
      repeated = repeated + 1
    enddo
    temp_repeated_list(k) = repeated

    well_grid%nsegments = k 
    nsegments = well_grid%nsegments

    allocate(well_grid%dh(nsegments))
    allocate(well_grid%res_dz(nsegments))
    allocate(well_grid%h(nsegments))
    allocate(well_grid%h_local_id(nsegments))
    allocate(well_grid%h_ghosted_id(nsegments))
    allocate(well_grid%h_global_id(nsegments))
    allocate(well_grid%h_rank_id(nsegments))
    allocate(well_grid%strata_id(nsegments))

    well_grid%h_local_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_ghosted_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_global_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_rank_id(:) = 0

    dh_x = diff_x/nsegments
    dh_y = diff_y/nsegments

    do k = 1,well_grid%nsegments
      well_grid%h(k)%id = k
      well_grid%h(k)%x = well_grid%bottomhole(1)+(dh_x*(k-0.5))
      well_grid%h(k)%y = well_grid%bottomhole(2)+(dh_y*(k-0.5))

      well_grid%dh(k) = temp_repeated_list(k)*well_grid%dz_peck
      if (k == 1) then
        well_grid%h(k)%z = well_grid%bottomhole(3)+(well_grid%dh(k)*(0.5d0))
      else
        well_grid%h(k)%z = well_grid%h(k-1)%z + (well_grid%dh(k-1)*(0.5d0)) + &
                           (well_grid%dh(k)*(0.5d0)) 
      endif
      
      well_grid%h_local_id(k) = temp_id_list(k)
      well_grid%h_ghosted_id(k) = res_grid%nL2G(well_grid%h_local_id(k))
      well_grid%h_global_id(k) = res_grid%nG2A(well_grid%h_ghosted_id(k))
    enddo

    well_grid%res_dz(:) = well_grid%dh(:)

    diff_x = diff_x*diff_x
    diff_y = diff_y*diff_y
    diff_z = diff_z*diff_z

    total_length = sqrt(diff_x+diff_y+diff_z)

    deallocate(temp_repeated_list)
    deallocate(temp_id_list)

  elseif (associated(well_grid%z_list) .or. associated(well_grid%l_list)) then
  !---------------------------------------------------------------------------
  ! Use the provided z list to build the grid

    if (.not.associated(well_grid%l_list)) then
      option%io_buffer = 'When providing SEGMENT_CENTER_Z_VALUES, you must & 
                         &also provide a list of SEGMENT_LENGTH_VALUES that &
                         &contains the same number of values (i.e. one length &
                         &value for every z-center value).'
      call PrintErrMsg(option)
    endif
    if (.not.associated(well_grid%z_list)) then
      option%io_buffer = 'When providing SEGMENT_LENGTH_VALUES, you must & 
                         &also provide a list of SEGMENT_CENTER_Z_VALUES that &
                         &contains the same number of values (i.e. one z-center&
                         & value for every length value).'
      call PrintErrMsg(option)
    endif
    if (size(well_grid%l_list) /= size(well_grid%z_list)) then
      option%io_buffer = 'The length of SEGMENT_LENGTH_VALUES must match the & 
                         &length of SEGMENT_CENTER_Z_VALUES (i.e. one z-center &
                         &value for every length value) provided in the &
                         &WELLBORE_MODEL,WELL_GRID block.'
      call PrintErrMsg(option)
    endif

    nsegments = well_grid%nsegments

    allocate(well_grid%dh(nsegments))
    allocate(well_grid%res_dz(nsegments))
    allocate(well_grid%h(nsegments))
    allocate(well_grid%h_local_id(nsegments))
    allocate(well_grid%h_ghosted_id(nsegments))
    allocate(well_grid%h_global_id(nsegments))
    allocate(well_grid%h_rank_id(nsegments))
    allocate(well_grid%strata_id(nsegments))

    well_grid%h_local_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_ghosted_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_global_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_rank_id(:) = UNINITIALIZED_INTEGER

    ! sort the z-list in ascending order, in case it was not provided that way
    allocate(temp_z_list(nsegments))
    temp_z_list = well_grid%z_list
    do i = 1, nsegments
      do j = i+1, nsegments
        if (well_grid%z_list(i) > well_grid%z_list(j)) then
          z_dum = well_grid%z_list(j)
          well_grid%z_list(j) = well_grid%z_list(i)
          well_grid%z_list(i) = z_dum
        endif
      enddo
    enddo
    ! if re-sorted doesn't match the given list, throw error
    do i = 1, nsegments
      if (well_grid%z_list(i) /= temp_z_list(i)) then
        option%io_buffer = 'The list of SEGMENT_CENTER_Z_VALUES must be & 
          &provided in ascending order. Ensure that the list of corresponding & 
          &SEGMENT_LENGTH_VALUES is also modified accordingly.'
        call PrintErrMsg(option)
      endif
    enddo
    deallocate(temp_z_list)

    dh_x = diff_x/nsegments
    dh_y = diff_y/nsegments

    do k = 1,well_grid%nsegments
      well_grid%h(k)%id = k
      well_grid%h(k)%x = well_grid%bottomhole(1)+(dh_x*(k-0.5))
      well_grid%h(k)%y = well_grid%bottomhole(2)+(dh_y*(k-0.5))
      well_grid%h(k)%z = well_grid%z_list(k)

      call GridGetLocalIDFromCoordinate(res_grid,well_grid%h(k), &
                                        option,local_id)
      if (Initialized(local_id)) then
        well_grid%h_local_id(k) = local_id
        well_grid%h_ghosted_id(k) = res_grid%nL2G(local_id)
        well_grid%h_global_id(k) = res_grid%nG2A(well_grid%h_ghosted_id(k))
        well_grid%h_rank_id(k) = option%myrank
      endif
    enddo

    well_grid%dh = well_grid%l_list

    temp_real = SUM(well_grid%dh)
    total_length = sqrt((diff_x*diff_x)+(diff_y*diff_y)+(diff_z*diff_z))
    write(string,'(1pe12.4)') temp_real
    write(string2,'(1pe12.4)') total_length
    if (temp_real /= total_length) then
      temp_real2 = abs(temp_real - total_length)
      !write(*,*) temp_real2
      if (temp_real2 > 1.d-2) then
        option%io_buffer = 'The sum of the list of SEGMENT_LENGTH_VALUES & 
          &(' // trim(string) // ' m) does not match the total length of the &
          &well according to the coordinates provided by WELLBORE_MODEL,&
          &TOP_OF_HOLE and WELLBORE_MODEL,TOP_OF_HOLE (' // trim(string2) // &
          ' m).'
        call PrintErrMsg(option)
      endif
    endif

    !This could be way off if # of well cells per reservoir cell >> 1
    !This could also lead to inconsisent well indices between the 
    !generated vs. read-in well.
    well_grid%res_dz(:) = well_grid%dh(:)

  elseif (Initialized(well_grid%nsegments)) then 
  !---------------------------------------------------------------------------
  ! Build an equally-spaced grid based on nsegments:

    nsegments = well_grid%nsegments

    allocate(well_grid%dh(nsegments))
    allocate(well_grid%res_dz(nsegments))
    allocate(well_grid%h(nsegments))
    allocate(well_grid%h_local_id(nsegments))
    allocate(well_grid%h_ghosted_id(nsegments))
    allocate(well_grid%h_global_id(nsegments))
    allocate(well_grid%h_rank_id(nsegments))
    allocate(well_grid%strata_id(nsegments))

    well_grid%h_local_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_ghosted_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_global_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_rank_id(:) = UNINITIALIZED_INTEGER

    dh_x = diff_x/nsegments
    dh_y = diff_y/nsegments
    dh_z = diff_z/nsegments

    diff_x = diff_x*diff_x
    diff_y = diff_y*diff_y
    diff_z = diff_z*diff_z

    total_length = sqrt(diff_x+diff_y+diff_z)

    do k = 1,well_grid%nsegments
      well_grid%h(k)%id = k
      well_grid%h(k)%x = well_grid%bottomhole(1)+(dh_x*(k-0.5))
      well_grid%h(k)%y = well_grid%bottomhole(2)+(dh_y*(k-0.5))
      well_grid%h(k)%z = well_grid%bottomhole(3)+(dh_z*(k-0.5))
    enddo

    well_grid%dh(:) = total_length/nsegments
    well_grid%res_dz(:) = well_grid%dh(:)

    do k = 1,well_grid%nsegments
      call GridGetLocalIDFromCoordinate(res_grid,well_grid%h(k), &
                                        option,local_id)
      if (Initialized(local_id)) then
        well_grid%h_local_id(k) = local_id
        well_grid%h_ghosted_id(k) = res_grid%nL2G(local_id)
        well_grid%h_global_id(k) = res_grid%nG2A(well_grid%h_ghosted_id(k))
        well_grid%h_rank_id(k) = option%myrank
      endif
    enddo
  !---------------------------------------------------------------------------
  else
   ! Use reservoir grid info
     dz_list = UNINITIALIZED_DOUBLE
     res_dz_list = UNINITIALIZED_DOUBLE

     if (Initialized(well_grid%min_dz)) then
       min_dz = well_grid%min_dz
     else
       min_dz = 1.d-5
     endif

     dz = min_dz
     z = well_grid%bottomhole(3)
     dummy_h%x = well_grid%bottomhole(1)
     dummy_h%y = well_grid%bottomhole(2)
     dummy_h%z = z
     call GridGetLocalIDFromCoordinate(res_grid,dummy_h,option,local_id)
     cell_id_list(1) = local_id
     cur_id = local_id
     z = z + min_dz
     dummy_h%z = z
     nsegments = 1
     nsegments_save = 0
     res_cell_count = 1
     do
       if (z > well_grid%tophole(3)) exit
       call GridGetLocalIDFromCoordinate(res_grid,dummy_h,option,local_id)

       res_cell_count = res_cell_count + 1
       if (res_cell_count <= well_grid%well_res_ratio .and. &
           cur_id == local_id) then
         nsegments = nsegments + 1
         cell_id_list(nsegments) = local_id
         dz = dz + min_dz
       elseif (cur_id /= local_id) then
         res_dz_list(nsegments_save+1:nsegments) = dz
         dz = dz / (nsegments-nsegments_save)
         dz_list(nsegments_save+1:nsegments) = dz
         res_cell_count = 1
         cur_id = local_id
         nsegments_save = nsegments
         nsegments = nsegments+1
         cell_id_list(nsegments) = local_id
         dz = min_dz
       else
         dz = dz + min_dz
       endif

       z = z + min_dz
       dummy_h%z = z
     enddo
     res_dz_list(nsegments_save+1:nsegments) = dz
     dz = dz / (nsegments-nsegments_save)
     dz_list(nsegments_save+1:nsegments) = dz

     allocate(well_grid%dh(nsegments))
     allocate(well_grid%res_dz(nsegments))
     allocate(well_grid%h(nsegments))
     allocate(well_grid%h_local_id(nsegments))
     allocate(well_grid%h_ghosted_id(nsegments))
     allocate(well_grid%h_global_id(nsegments))
     allocate(well_grid%h_rank_id(nsegments))
     allocate(well_grid%strata_id(nsegments))

     well_grid%h_rank_id(:) = 0
     well_grid%res_dz(1:nsegments) = res_dz_list(1:nsegments)

     well_grid%strata_id(:) = UNINITIALIZED_INTEGER
     well_grid%nsegments = nsegments
     well_grid%nconnections = well_grid%nsegments - 1

     well_grid%h(1)%id = 1
     well_grid%h(1)%x = well_grid%bottomhole(1)
     well_grid%h(1)%y = well_grid%bottomhole(2)
     well_grid%h(1)%z = well_grid%bottomhole(3) + dz_list(1)/2.d0
     well_grid%dh(1) = dz_list(1)

     local_id = cell_id_list(1)
     well_grid%h_local_id(1) = local_id
     well_grid%h_ghosted_id(1) = res_grid%nL2G(local_id)
     well_grid%h_global_id(1) = res_grid%nG2A(well_grid%h_ghosted_id(1))

     do k = 2,nsegments
       well_grid%h(k)%id = k
       well_grid%h(k)%x = well_grid%bottomhole(1)
       well_grid%h(k)%y = well_grid%bottomhole(2)
       well_grid%h(k)%z = well_grid%bottomhole(3) + &
                          sum(dz_list(1:k-1)) + dz_list(k)/2.d0
       well_grid%dh(k) = dz_list(k)

       local_id = cell_id_list(k)
       well_grid%h_local_id(k) = local_id
       well_grid%h_ghosted_id(k) = res_grid%nL2G(local_id)
       well_grid%h_global_id(k) = res_grid%nG2A(well_grid%h_ghosted_id(k))

     enddo
  !---------------------------------------------------------------------------
  endif

  allocate(well_grid%strata_id(nsegments))
  well_grid%strata_id(:) = UNINITIALIZED_INTEGER
  well_grid%nconnections = well_grid%nsegments - 1

end subroutine PMWellSetupGrid

subroutine PMWellSetup(this)
  !
  ! Initializes variables associated with the well process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  !

  use Grid_module
  use Utility_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Input_Aux_module
  use Dataset_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Transport_Constraint_NWT_module
  use NW_Transport_Aux_module
  use SCO2_Aux_module

  implicit none

  class(pm_well_type) :: this

  type(option_type), pointer :: option
  type(grid_type), pointer :: res_grid
  type(well_grid_type), pointer :: well_grid
  type(coupler_type), pointer :: source_sink
  type(input_type) :: input_dummy
  class(realization_subsurface_type), pointer :: realization
  type(point3d_type) :: dummy_h
  class(tran_constraint_coupler_nwt_type), pointer ::tran_constraint_coupler_nwt
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscInt, allocatable :: h_global_id_unique(:)
  PetscInt, allocatable :: h_rank_id_unique(:)
  PetscInt, allocatable :: h_all_rank_id(:)
  PetscInt, allocatable :: h_all_global_id(:)
  PetscReal :: max_diameter, xy_span
  PetscReal :: temp_real
  PetscInt :: local_id
  PetscInt :: local_id_well, local_id_res
  PetscInt :: nsegments
  PetscInt :: max_val, min_val
  PetscInt :: k, j
  PetscInt :: count1, count2_local, count2_global
  PetscBool :: well_grid_res_is_OK = PETSC_FALSE
  PetscBool :: res_grid_cell_within_well_z
  PetscBool :: res_grid_cell_within_well_y
  PetscBool :: res_grid_cell_within_well_x
  PetscErrorCode :: ierr
  PetscInt :: well_bottom_local, well_bottom_ghosted
  PetscInt, allocatable :: temp(:), temp2(:)

  option => this%option
  realization => this%realization
  res_grid => realization%patch%grid
  well_grid => this%well_grid

  well_bottom_local = ZERO_INTEGER
  well_bottom_ghosted = ZERO_INTEGER

  allocate(this%well_pert(option%nflowdof))
  do k = 1,option%nflowdof
    call PMWellVarCreate(this%well_pert(k))
  enddo

  option%io_buffer = ' '
  call PrintMsg(option)
  option%io_buffer = 'WELLBORE_MODEL: Creating well grid discretization.... '
  call PrintMsg(option)

  call PMWellSetupGrid(well_grid,res_grid,option)

  write(string,'(I0.5)') well_grid%nsegments
  option%io_buffer = 'WELLBORE_MODEL: Grid created with ' // trim(string) // &
                     ' segments. '
  call PrintMsg(option)

  nsegments = well_grid%nsegments

  ! Create a well MPI communicator
  this%well_comm%petsc_rank = option%myrank
  allocate(h_all_rank_id(nsegments))
  allocate(h_all_global_id(nsegments))
  allocate(h_rank_id_unique(nsegments))

  call MPI_Allreduce(well_grid%h_rank_id,h_all_rank_id,nsegments, &
                     MPI_INTEGER,MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(well_grid%h_global_id,h_all_global_id,nsegments, &
                     MPI_INTEGER,MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  well_grid%h_rank_id = h_all_rank_id
  well_grid%h_global_id = h_all_global_id

  h_rank_id_unique(:) = UNINITIALIZED_INTEGER  

  min_val = minval(h_all_rank_id)-1
  max_val = maxval(h_all_rank_id)
  !Always include rank 0 in comm
  if (min_val > -1) then
    h_rank_id_unique(1) = 0
    k = 1
  else
    k = 0
  endif
  do while (min_val < max_val)
    k = k + 1
    min_val = minval(h_all_rank_id, mask=h_all_rank_id > min_val)
    h_rank_id_unique(k) = min_val
  enddo

  allocate(this%well_comm%petsc_rank_list(k)) 
  allocate(this%well_comm%well_rank_list(k))
  this%well_comm%petsc_rank_list = UNINITIALIZED_INTEGER
  this%well_comm%well_rank_list = UNINITIALIZED_INTEGER

  this%well_comm%petsc_rank_list = h_rank_id_unique(1:k)
  this%well_comm%commsize = k 
  call MPI_Group_incl(option%driver%comm%group,this%well_comm%commsize, &
                      this%well_comm%petsc_rank_list,this%well_comm%group, &
                      ierr);CHKERRQ(ierr)
  call MPI_Comm_create(option%driver%comm%communicator,this%well_comm%group, &
                       this%well_comm%comm,ierr);CHKERRQ(ierr)

  if (this%well_comm%comm /= MPI_COMM_NULL) then
    call MPI_Comm_rank(this%well_comm%comm,this%well_comm%rank, &
                       ierr);CHKERRQ(ierr)
  endif
  do j = 0,(this%well_comm%commsize-1)
    this%well_comm%well_rank_list(j+1) = j
  enddo
    
  allocate(this%well%r0(nsegments))

  if (size(this%well%diameter) /= nsegments) then
    if (size(this%well%diameter) == 1) then
      temp_real = this%well%diameter(1)
      deallocate(this%well%diameter)
      allocate(this%well%diameter(nsegments))
      this%well%diameter(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,WELL,DIAMETER must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,WELL,DIAMETER, &
        &it will be applied to all segments of the well.'
      call PrintErrMsg(option)
    endif
  endif
  max_diameter = maxval(this%well%diameter)
  xy_span = well_grid%xy_span_multiplier*max_diameter

  if (size(this%well%WI_base) /= nsegments) then
    if (size(this%well%WI_base) == 1) then
      temp_real = this%well%WI_base(1)
      deallocate(this%well%WI_base)
      allocate(this%well%WI_base(nsegments))
      this%well%WI_base(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,WELL,WELL_INDEX must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,WELL,WELL_INDEX, &
        &it will be applied to all segments of the well.'
      call PrintErrMsg(option)
    endif
  endif

  if (size(this%well%friction_factor) /= nsegments) then
    if (size(this%well%friction_factor) == 1) then
      temp_real = this%well%friction_factor(1)
      deallocate(this%well%friction_factor)
      allocate(this%well%friction_factor(nsegments))
      this%well%friction_factor(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,WELL,FRICTION_COEFFICIENT must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,WELL, &
        &FRICTION_COEFFICIENT, it will be applied to all segments of the well.'
      call PrintErrMsg(option)
    endif
  endif

  if (size(this%well%skin) /= nsegments) then
    if (size(this%well%skin) == 1) then
      temp_real = this%well%skin(1)
      deallocate(this%well%skin)
      allocate(this%well%skin(nsegments))
      this%well%skin(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,WELL,SKIN_FACTOR must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,WELL,SKIN_FACTOR, &
        &it will be applied to all segments of the well.'
      call PrintErrMsg(option)
    endif
  endif

  ! Check that no reservoir grid cells were skipped.
  ! Count how many of the h_global_id's are unique.
  ! This sum must be = to the number of reservoir cells that the
  !   well passes through.
  option%io_buffer = 'WELLBORE_MODEL: Checking well grid resolution against &
                     &reservoir grid.... '
  call PrintMsg(option)

  allocate(h_global_id_unique(nsegments))
  h_global_id_unique(:) = 0

  min_val = minval(h_all_global_id)-1
  max_val = maxval(h_all_global_id)
  k = 0
  do while (min_val < max_val)
    k = k + 1
    min_val = minval(h_all_global_id, mask=h_all_global_id > min_val)
    h_global_id_unique(k) = min_val
  enddo

  count1 = 0
  do k = 1,nsegments
    if (h_global_id_unique(k) > 0) then
      count1 = count1 + 1
    endif
  enddo
  write(string,'(I0.5)') count1

  deallocate(h_global_id_unique)
  deallocate(h_rank_id_unique)
  deallocate(h_all_rank_id)
  deallocate(h_all_global_id)

  ! Next, sum up how many grid cells the well passes thru.
  ! Note: This count assumes that the well is vertical and the top and
  !       bottom surfaces do not slope or undulate.
  count2_local = 0
  count2_global = 0
  do k = 1,res_grid%ngmax
    res_grid_cell_within_well_z = PETSC_FALSE
    res_grid_cell_within_well_y = PETSC_FALSE
    res_grid_cell_within_well_x = PETSC_FALSE
    if ( (res_grid%z(k) >= well_grid%bottomhole(3)) .and. &
         (res_grid%z(k) <= well_grid%tophole(3)) ) then
      res_grid_cell_within_well_z = PETSC_TRUE
    endif
    if ( (res_grid%y(k) >= (well_grid%tophole(2)-xy_span)) .and. &
         (res_grid%y(k) <= (well_grid%tophole(2)+xy_span)) ) then
      res_grid_cell_within_well_y = PETSC_TRUE
    endif
    if ( (res_grid%x(k) >= (well_grid%tophole(1)-xy_span)) .and. &
         (res_grid%x(k) <= (well_grid%tophole(1)+xy_span)) ) then
      res_grid_cell_within_well_x = PETSC_TRUE
    endif
    if (res_grid_cell_within_well_z .and. res_grid_cell_within_well_y &
        .and. res_grid_cell_within_well_x) then

      ! What should the local_id of the reservoir cell at this z-level
      ! along the well be?
      dummy_h%z = res_grid%z(k)
      dummy_h%y = well_grid%tophole(2)
      dummy_h%x = well_grid%tophole(1)
      call GridGetLocalIDFromCoordinate(res_grid,dummy_h,option,local_id)
      local_id_well = local_id
      ! What is the ghosted_id of the actual reservoir cell at this z-level?
      dummy_h%z = res_grid%z(k)
      dummy_h%y = res_grid%y(k)
      dummy_h%x = res_grid%x(k)
      call GridGetLocalIDFromCoordinate(res_grid,dummy_h,option,local_id)
      local_id_res = local_id
      ! Does the well occupy this grid cell? If a reservoir cell was skipped
      ! by the well, then the count will never be incremented, and count2
      ! will not equal count1. count2 will be larger than count1.
      if ((local_id_res == local_id_well) .and. Initialized(local_id_res)) then
        count2_local = count2_local + 1
      endif

    endif
  enddo

  ! All of the MPI processes need to sum up their counts and place the
  ! total in count2_global.
  if (this%well_comm%comm /= MPI_COMM_NULL) then
    call MPI_Allreduce(count2_local,count2_global,ONE_INTEGER_MPI,MPI_INTEGER, &
                       MPI_SUM,this%well_comm%comm,ierr);CHKERRQ(ierr)
  endif
  call MPI_Bcast(count2_global,ONE_INTEGER_MPI,MPI_INTEGER, &
                 this%well_comm%petsc_rank_list(1),option%mycomm, &
                 ierr);CHKERRQ(ierr)
  write(string2,'(I0.5)') count2_global

  ! The only way we can ensure that the well discretization did not skip a
  ! reservoir cell, is if the number of unique global_id's that the well
  ! is connected to (count1) matches the number of reservoir grid cells that
  ! the well occupies (count2):
  if (count1 == count2_global) well_grid_res_is_OK = PETSC_TRUE

  do k = 1,nsegments
    if (Uninitialized(well_grid%h_local_id(k))) well_grid%h_local_id(k) = -1
  enddo


  if (.not.well_grid_res_is_OK) then
    option%io_buffer = 'ERROR:  &
      &The number of reservoir grid cells that are occupied by the well &
      &(' // trim(string2) // ') is larger than the number of unique &
      &reservoir grid cells that have a connection to the well (' // &
      trim(string) // '). Therefore, some of the reservoir grid cells &
      &have been skipped and have no connection to the well. You must &
      &increase the resolution of the WELLBORE_MODEL grid. '
    call PrintMsg(option)
    if (well_grid%match_reservoir) then
      option%io_buffer = '(cont.) You should try &
        &decreasing the value set for the dz step parameter &
        &WELLBORE_MODEL,WELL_GRID,MINIMUM_DZ_STEP (default value = 1.0d-2 m).'
      call PrintErrMsg(option)
    else
      option%io_buffer = '(see above) Alternatively, &
        &if you are sure your well grid resolution is fine enough, try &
        &increasing the value set for the x-y search parameter &
        &WELLBORE_MODEL,WELL_GRID,XY_SEARCH_MULTIPLIER (default value = 10).'
      call PrintErrMsg(option)
    endif
    
  else
    option%io_buffer = 'WELLBORE_MODEL: &
      &Well grid resolution is adequate. No reservoir grid cell &
      &connections have been skipped.'
    call PrintMsg(option)
    if (well_grid%match_reservoir) then
      option%io_buffer = 'WELLBORE_MODEL: &
        &For your convenience, the SEGMENT_CENTER_Z_VALUES (meters) are: '
      call PrintMsg(option)
      write(*,*) well_grid%h%z
      option%io_buffer = 'WELLBORE_MODEL: &
        &For your convenience, the SEGMENT_LENGTH_VALUES (meters) are: '
      call PrintMsg(option)
      write(*,*) well_grid%dh
    endif
  endif

  this%flow_soln%ndof = this%nphase

  if (option%itranmode /= NULL_MODE) then
    this%transport = PETSC_TRUE
    if (option%itranmode /= NWT_MODE) then
      option%io_buffer ='The only transport mode allowed with the &
      &WELLBORE_MODEL is NWT_MODE.'
      call PrintErrMsg(option)
    endif
    this%nspecies = realization%reaction_nw%params%nspecies
    this%tran_soln%ndof = this%nspecies
  endif

  ! add a reservoir src/sink coupler for each well segment
  do k = 1,well_grid%nsegments
    if (well_grid%h_rank_id(k) /= option%myrank) cycle

    write(string,'(I0.6)') k
    source_sink => CouplerCreate(SRC_SINK_COUPLER_TYPE)
    source_sink%name = 'well_segment_' // trim(string)

    ! ----- flow ------------------
    source_sink%flow_condition_name = 'well_segment_' // trim(string) // &
                                      '_flow_srcsink'
    source_sink%flow_condition => FlowConditionCreate(option)
    source_sink%flow_condition%name = source_sink%flow_condition_name
    select case (option%iflowmode)
    case(WF_MODE)
      source_sink%flow_condition%general => FlowGeneralConditionCreate(option)
      string = 'RATE'
      source_sink%flow_condition%general%rate => FlowGeneralSubConditionPtr( &
        input_dummy,string,source_sink%flow_condition%general,option)
      source_sink%flow_condition%general%rate%itype = SCALED_MASS_RATE_SS ! [kg/s]
      source_sink%flow_condition%general%liquid_pressure => &
            FlowGeneralSubConditionPtr(input_dummy,string,source_sink% &
                                       flow_condition%general,option)
      source_sink%flow_condition%general%gas_pressure => &
            FlowGeneralSubConditionPtr(input_dummy,string,source_sink% &
                                      flow_condition%general,option)
      allocate(source_sink%flow_condition%general%rate%dataset%rarray(2))
      source_sink%flow_condition%general%rate%dataset%rarray(:) = 0.d0

      source_sink%flow_condition%well => FlowSubConditionCreate(ONE_INTEGER)

    case(SCO2_MODE)
      source_sink%flow_condition%sco2 => FlowSCO2ConditionCreate(option)
      string = 'RATE'
      source_sink%flow_condition%sco2%rate => FlowSCO2SubConditionPtr( &
        input_dummy,string,source_sink%flow_condition%sco2,option)
      source_sink%flow_condition%sco2%rate%itype = SCALED_MASS_RATE_SS ! [kg/s]
      source_sink%flow_condition%sco2%liquid_pressure => &
            FlowSCO2SubConditionPtr(input_dummy,string,source_sink% &
                                       flow_condition%sco2,option)
      source_sink%flow_condition%sco2%gas_pressure => &
            FlowSCO2SubConditionPtr(input_dummy,string,source_sink% &
                                       flow_condition%sco2,option)
      allocate(source_sink%flow_condition%sco2%rate%dataset%rarray(2))
      source_sink%flow_condition%sco2%rate%dataset%rarray(:) = 0.d0

      source_sink%flow_condition%well => FlowSubConditionCreate(ONE_INTEGER)
      
      ! Bottom of hole is special for fully implicit coupling with steady
      ! state well model.
      if (k==1) then
        source_sink%flow_condition%well%ctype = 'well-bottom'
      else
        source_sink%flow_condition%well%ctype = 'well'
      endif

    end select

    ! ----- transport -------------
    if (this%transport) then
      write(string,'(I0.6)') k
      source_sink%tran_condition_name = 'well_segment_' // trim(string) // &
                                        '_tran_srcsink'
      source_sink%tran_condition => TranConditionCreate(option)
      source_sink%tran_condition%name = source_sink%tran_condition_name
      source_sink%tran_condition%itype = WELL_SS
      tran_constraint_coupler_nwt => TranConstraintCouplerNWTCreate(option)
      allocate(tran_constraint_coupler_nwt%nwt_auxvar)
      call NWTAuxVarInit(tran_constraint_coupler_nwt%nwt_auxvar,&
                         this%realization%reaction_nw,option)
      tran_constraint_coupler_nwt%constraint_name = &
                                               'well_segment_' // trim(string)
      source_sink%tran_condition%cur_constraint_coupler => &
                                                   tran_constraint_coupler_nwt
    endif

    source_sink%connection_set => ConnectionCreate(1,SRC_SINK_CONNECTION_TYPE)
    source_sink%connection_set%id_dn = well_grid%h_local_id(k)

    call CouplerAddToList(source_sink,this%realization%patch%source_sink_list)
    nullify(source_sink)
  enddo

  ! For fully-implicit well coupling, resize the matrix zeroing arrays to 
  ! exclude the bottom of the well for hydrostatic well model.

  select case (option%iflowmode)
    case(SCO2_MODE)
      well_bottom_ghosted = well_grid%h_ghosted_id(1)
      well_bottom_local = well_grid%h_local_id(1)
      if (well_bottom_ghosted > 0 .and. &
          realization%patch%aux%sco2%inactive_cells_exist) then
        if (size(realization%patch%aux%sco2%matrix_zeroing%inactive_rows_local) &
            == 1) then
          deallocate(realization%patch%aux%sco2%matrix_zeroing% &
                     inactive_rows_local)
          deallocate(realization%patch%aux%sco2%matrix_zeroing%&
                     inactive_rows_local_ghosted)
          realization%patch%aux%sco2%inactive_cells_exist = PETSC_FALSE
          realization%patch%aux%sco2%matrix_zeroing%n_inactive_rows = 0
        else
          if (allocated(temp)) deallocate(temp)
          if (allocated(temp2)) deallocate(temp2)
          allocate(temp(size(realization%patch%aux%sco2% &
               matrix_zeroing%inactive_rows_local)))
          allocate(temp2(size(realization%patch%aux%sco2% &
               matrix_zeroing%inactive_rows_local)-1))
          realization%patch%aux%sco2%matrix_zeroing%n_inactive_rows = &
               realization%patch%aux%sco2%matrix_zeroing%n_inactive_rows - 1
          temp(:) = realization%patch%aux%sco2%matrix_zeroing% &
                    inactive_rows_local(:)
          temp2 = pack(temp,temp /= well_bottom_local*option%nflowdof)
          deallocate(realization%patch%aux%sco2%matrix_zeroing%inactive_rows_local)
          allocate(realization%patch%aux%sco2%matrix_zeroing% &
                   inactive_rows_local(size(temp2)))
          realization%patch%aux%sco2%matrix_zeroing%inactive_rows_local(:) = temp2(:)
          temp(:) = realization%patch%aux%sco2%matrix_zeroing% &
                    inactive_rows_local_ghosted(:)
          temp2 = pack(temp,temp /= well_bottom_local*option%nflowdof-1)
          deallocate(realization%patch%aux%sco2%matrix_zeroing% &
                     inactive_rows_local_ghosted)
          allocate(realization%patch%aux%sco2%matrix_zeroing% &
                   inactive_rows_local_ghosted(size(temp2)))
          realization%patch%aux%sco2%matrix_zeroing% &
                   inactive_rows_local_ghosted(:) = temp2(:)
        endif
      endif
  end select
end subroutine PMWellSetup

! ************************************************************************** !

subroutine PMWellReadSimOptionsBlock(this,input)
  !
  ! Read simulation options for the well model
  !
  ! Author: Michael Nole
  ! Date: 03/08/24
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword, word
  class(pm_well_type) :: this
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  option => this%option

  error_string = 'Well Model Options'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      case('FLOW_COUPLING')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case (trim(word))
          case('FULLY_IMPLICIT')
            this%flow_coupling = FULLY_IMPLICIT_WELL
          case('QUASI_IMPLICIT')
            this%flow_coupling = QUASI_IMPLICIT_WELL
          case('SEQUENTIAL')
            this%flow_coupling = SEQUENTIAL_WELL
        end select
      case('WELL_MODEL_TYPE','TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !-----------------------------
          case('HYDROSTATIC')
            this%well%well_model_type = WELL_MODEL_HYDROSTATIC
        !-----------------------------
          case('WIPP_DARCY')
            this%well%well_model_type = WELL_MODEL_WIPP_DARCY
        !-----------------------------
          case('FULL_MOMENTUM')
            this%well%well_model_type = WELL_MODEL_FULL_MOMENTUM
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      case default
          call InputKeywordUnrecognized(input,keyword,'Well Model',option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine PMWellReadSimOptionsBlock

! ************************************************************************** !

subroutine PMWellReadPMBlock(this,input)
  !
  ! Reads input file parameters associated with the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  use Input_Aux_module
  use String_module
  use Utility_module

  implicit none

  class(pm_well_type) :: this
  type(input_type), pointer :: input

  type(option_type), pointer :: option
  type(well_grid_type), pointer :: well_grid
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option
  well_grid => this%well_grid
  input%ierr = 0
  error_string = 'WELLBORE_MODEL'

  option%io_buffer = 'pflotran card:: WELLBORE_MODEL'
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

    ! Read keywords within WELLBORE_MODEL block:
    select case(trim(word))
    !-------------------------------------
      case('SKIP_RESTART')
        this%skip_restart = PETSC_TRUE
        cycle
    !-------------------------------------
      case('SINGLE_PHASE')
        this%nphase = 1
        cycle
    !-------------------------------------
      case('TWO_PHASE')
        this%nphase = 2
        cycle
    !-------------------------------------
      case('PRINT_WELL_FILE')
        this%print_well = PETSC_TRUE
        cycle
    !-------------------------------------
      case('CHECK_FOR_SS')
        this%ss_check = PETSC_TRUE
        cycle
    !-------------------------------------
      case('WIPP_INTRUSION_START_TIME')
        call InputReadDouble(input,option,this%intrusion_time_start)
        call InputErrorMsg(input,option,'WIPP_INTRUSION_START_TIME', &
                           error_string)
        call InputReadAndConvertUnits(input,this%intrusion_time_start,'sec', &
                           'WELLBORE_MODEL, WIPP_INTRUSION_START_TIME',option)
        this%well_on = PETSC_FALSE
        cycle
    !-------------------------------------
      case('WIPP_INTRUSION_ZERO_VALUE')  ! [mol/m3-bulk]
        call InputReadDouble(input,option,this%bh_zero_value)
        call InputErrorMsg(input,option,'WIPP_INTRUSION_ZERO_VALUE', &
                           error_string)
        cycle
    !-------------------------------------
    end select

    ! Read sub-blocks within WELLBORE_MODEL block:
    error_string = 'WELLBORE_MODEL'
    call PMWellReadGrid(well_grid,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadWell(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadWellBCs(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadFlowSolver(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadTranSolver(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadWellModelType(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadWellConstraintType(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'WELLBORE_MODEL'
    call PMWellReadWellOutput(this,input,option,word,error_string,found)
    if (found) cycle

    if (.not. found) then
      option%io_buffer = 'Keyword "' // trim(word) // &
                         '" does not exist for WELLBORE_MODEL.'
      call PrintErrMsg(option)
    endif

  enddo
  call InputPopBlock(input,option)

  if (this%nphase == 0) then
    option%io_buffer = 'The number of fluid phases must be indicated &
      &in the WELLBORE_MODEL block using one of these keywords: &
      &SINGLE_PHASE, TWO_PHASE.'
    call PrintErrMsg(option)
  endif

end subroutine PMWellReadPMBlock

! ************************************************************************** !

subroutine PMWellReadGrid(well_grid,input,option,keyword,error_string,found)
  !
  ! Reads input file parameters associated with the well model grid.
  !
  ! Author: Jenn Frederick
  ! Date: 11/03/2021
  !
  use Input_Aux_module
  use String_module

  implicit none

  type(well_grid_type), pointer :: well_grid
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscInt :: num_errors, num_read
  PetscInt :: read_max = 10000
  PetscInt :: k
  character(len=MAXWORDLENGTH) :: word
  PetscReal, pointer :: temp_z_list(:)
  PetscReal, pointer :: temp_l_list(:)

  error_string = trim(error_string) // ',WELL_GRID'
  found = PETSC_TRUE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('WELL_GRID')
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
          case('MATCH_RESERVOIR')
            well_grid%match_reservoir = PETSC_TRUE
        !-------------------------------------
          case('MINIMUM_DZ_STEP')
            call InputReadDouble(input,option,well_grid%dz_peck)
            call InputErrorMsg(input,option,'MINIMUM_DZ_STEP', &
                               error_string)
        !-------------------------------------
          case('MAX_WELL_RESERVOIR_CELL_RATIO')
            call InputReadInt(input,option,well_grid%well_res_ratio)
            call InputErrorMsg(input,option,'MAX_WELL_RESERVOIR_CELL_RATIO', &
                               error_string)
        !-------------------------------------
          case('MIN_DZ')
            call InputReadDouble(input,option,well_grid%min_dz)
            call InputErrorMsg(input,option,'MIN_DZ', &
                               error_string)
        !-----------------------------
          case('NUMBER_OF_SEGMENTS')
            call InputReadInt(input,option,well_grid%nsegments)
            call InputErrorMsg(input,option,'NUMBER_OF_SEGMENTS',error_string)
        !-----------------------------
          case('SEGMENT_CENTER_Z_VALUES')
            num_read = 0
            allocate(temp_z_list(read_max))
            temp_z_list(:) = UNINITIALIZED_DOUBLE
            do k = 1,read_max
              call InputReadDouble(input,option,temp_z_list(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for &
                &SEGMENT_CENTER_Z_VALUES must be &
                &provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(well_grid%z_list(num_read))
            well_grid%z_list(1:num_read) = temp_z_list(1:num_read)
            well_grid%nsegments = num_read
            deallocate(temp_z_list)
        !-----------------------------
          case('SEGMENT_LENGTH_VALUES')
            num_read = 0
            allocate(temp_l_list(read_max))
            temp_l_list(:) = UNINITIALIZED_DOUBLE
            do k = 1,read_max
              call InputReadDouble(input,option,temp_l_list(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for &
                &SEGMENT_LENGTH_VALUES must be &
                &provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(well_grid%l_list(num_read))
            well_grid%l_list(1:num_read) = temp_l_list(1:num_read)
            well_grid%nsegments = num_read
            deallocate(temp_l_list)
        !-----------------------------
          case('TOP_OF_HOLE')
            call InputReadDouble(input,option,well_grid%tophole(1))
            call InputErrorMsg(input,option,'TOP_OF_HOLE x-coordinate', &
                               error_string)
            call InputReadDouble(input,option,well_grid%tophole(2))
            call InputErrorMsg(input,option,'TOP_OF_HOLE y-coordinate', &
                               error_string)
            call InputReadDouble(input,option,well_grid%tophole(3))
            call InputErrorMsg(input,option,'TOP_OF_HOLE z-coordinate', &
                               error_string)
        !-----------------------------
          case('BOTTOM_OF_HOLE')
            call InputReadDouble(input,option,well_grid%bottomhole(1))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE x-coordinate', &
                               error_string)
            call InputReadDouble(input,option,well_grid%bottomhole(2))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE y-coordinate', &
                               error_string)
            call InputReadDouble(input,option,well_grid%bottomhole(3))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE z-coordinate', &
                               error_string)
        !-----------------------------
          case('XY_SEARCH_MULTIPLIER')
            call InputReadInt(input,option,well_grid%xy_span_multiplier)
            call InputErrorMsg(input,option,'XY_SEARCH_MULTIPLIER',error_string)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      ! ----------------- error messaging -------------------------------------
!      if (well_grid%match_reservoir .eqv. PETSC_FALSE) then
!        if (Uninitialized(well_grid%nsegments)) then
!          option%io_buffer = 'ERROR: NUMBER_OF_SEGMENTS must be specified &
!                             &in the ' // trim(error_string) // ' block, or &
!                             &SEGMENT_CENTER_Z_VALUES must be provided.'
!          call PrintMsg(option)
!          num_errors = num_errors + 1
!        endif
!        if (well_grid%nsegments < 3) then
!          option%io_buffer = 'ERROR: The well must consist of >= 3 segments &
!                             &in the ' // trim(error_string) // ' block.'
!          call PrintMsg(option)
!          num_errors = num_errors + 1
!        endif
!      endif
      if (well_grid%match_reservoir) then
        !if (Uninitialized(well_grid%min_dz)) then
        !  option%io_buffer = 'ERROR: The minimum z spacing (MIN_DZ) must be &
        !                     &given in the ' // trim(error_string) // ' block.'
        !  call PrintMsg(option)
        !  num_errors = num_errors + 1
        !endif
      endif
      if (Uninitialized(well_grid%tophole(1)) .or. &
          Uninitialized(well_grid%tophole(2)) .or. &
          Uninitialized(well_grid%tophole(3))) then
        option%io_buffer = 'ERROR: TOP_OF_HOLE must be specified &
                           &in the ' // trim(error_string) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Uninitialized(well_grid%bottomhole(1)) .or. &
          Uninitialized(well_grid%bottomhole(2)) .or. &
          Uninitialized(well_grid%bottomhole(3))) then
        option%io_buffer = 'ERROR: BOTTOM_OF_HOLE must be specified &
                           &in the ' // trim(error_string) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the WELLBORE_MODEL,WELL_GRID block. See above &
                       &error messages.'
    call PrintErrMsg(option)
  endif

  end subroutine PMWellReadGrid

! ************************************************************************** !

subroutine PMWellReadWell(pm_well,input,option,keyword,error_string,found)
  !
  ! Reads input file parameters associated with the well model
  ! well properties.
  !
  ! Author: Jenn Frederick
  ! Date: 11/03/2021
  !
  use Input_Aux_module
  use String_module

  implicit none

  class(pm_well_type) :: pm_well
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscInt :: k
  PetscInt :: num_errors,read_max,num_read
  PetscInt :: at_index,num_segs
  character(len=MAXWORDLENGTH) :: word, index_word
  PetscReal :: index_val
  PetscReal, pointer :: temp_diameter(:)
  PetscReal, pointer :: temp_friction(:)
  PetscReal, pointer :: temp_well_index(:)
  PetscReal, pointer :: temp_well_perm(:)
  PetscReal, pointer :: temp_well_phi(:)
  PetscReal, pointer :: temp(:)

  error_string = trim(error_string) // ',WELL'
  found = PETSC_TRUE
  num_errors = 0

  read_max = 9000
  allocate(temp_diameter(read_max))
  allocate(temp_friction(read_max))
  allocate(temp_well_index(read_max))
  allocate(temp_well_perm(read_max))
  allocate(temp_well_phi(read_max))
  allocate(temp(read_max))
  temp_diameter(:) = UNINITIALIZED_DOUBLE
  temp_friction(:) = UNINITIALIZED_DOUBLE
  temp_well_index(:) = UNINITIALIZED_DOUBLE
  temp_well_perm(:) = UNINITIALIZED_DOUBLE
  temp_well_phi(:) = UNINITIALIZED_DOUBLE
  temp(:) = UNINITIALIZED_DOUBLE

  select case(trim(keyword))
  !-------------------------------------
    case('WELL')
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        num_read = 0
        select case(trim(word))
        !-----------------------------
          case('DIAMETER')
            do k = 1,read_max
              call InputReadDouble(input,option,temp_diameter(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for DIAMETER must be &
                &provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(pm_well%well%diameter(num_read))
            pm_well%well%diameter(1:num_read) = temp_diameter(1:num_read)
        !-----------------------------
          case('FRICTION_COEFFICIENT')
            do k = 1,read_max
              call InputReadDouble(input,option,temp_friction(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for FRICTION_COEFFICIENT &
                &must be provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(pm_well%well%friction_factor(num_read))
            pm_well%well%friction_factor(1:num_read) = temp_friction(1:num_read)
        !-----------------------------
          case('SKIN_FACTOR')
            do k = 1,read_max
              call InputReadDouble(input,option,temp(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for SKIN_FACTOR &
                &must be provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(pm_well%well%skin(num_read))
            pm_well%well%skin(1:num_read) = temp(1:num_read)
        !-----------------------------
          case('WELL_INDEX')
            do k = 1,read_max
              index_word = trim(input%buf)
              call InputReadDouble(input,option,temp_well_index(k))
              if (InputError(input)) then
                at_index = 0; at_index = index(index_word,'@')
                if (at_index > 0) then
                  read(index_word(1:at_index-1), *) num_segs
                  read(index_word(at_index+1:len(index_word)), *) index_val
                  temp_well_index(k:k+num_segs-1) = index_val
                  num_read = num_read + num_segs
                endif
                exit
              endif 
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for WELL_INDEX &
                &must be provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            write(*,*) 'num_read =', num_read
            allocate(pm_well%well%WI_base(num_read))
            pm_well%well%WI_base(1:num_read) = temp_well_index(1:num_read)
        !-----------------------------
          case('WELL_INDEX_MODEL')
            call InputReadWord(input,option,word,PETSC_TRUE)
            select case(word)
              case('PEACEMAN_ISO')
                pm_well%well%WI_model = PEACEMAN_ISO 
              case('PEACEMAN_2D')
                pm_well%well%WI_model = PEACEMAN_2D
              case('PEACEMAN_3D')
                pm_well%well%WI_model = PEACEMAN_3D
              case default
                option%io_buffer = 'Unrecognized option for WELL_INDEX_MODEL &
                &in the ' // trim(error_string) // ' block. Default is 3D &
                &Peaceman (PEACEMAN_3D).'
              call PrintErrMsg(option)
            end select
            call InputErrorMsg(input,option,'WELL_INDEX_MODEL',error_string)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      ! ----------------- error messaging -------------------------------------
      if (.not.associated(pm_well%well%WI_base)) then
        option%io_buffer = 'Keyword WELL_INDEX must be provided in &
                           &the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
      endif
      if (.not.associated(pm_well%well%friction_factor)) then
        option%io_buffer = 'Keyword FRICTION_COEFFICIENT must be provided in &
                           &the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
      endif
      if (.not.associated(pm_well%well%diameter)) then
        option%io_buffer = 'Keyword DIAMETER must be provided in &
                           &the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
      endif

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                  &the WELLBORE_MODEL,WELL block. See above error messages.'
    call PrintErrMsg(option)
  endif

  deallocate(temp_well_index)
  deallocate(temp_friction)
  deallocate(temp_diameter)
  deallocate(temp)

  end subroutine PMWellReadWell

! ************************************************************************** !

subroutine PMWellReadWellBCs(pm_well,input,option,keyword,error_string,found)
  !
  ! Reads input file parameters associated with the well model
  ! boundary conditions.
  !
  ! Author: Jenn Frederick
  ! Date: 12/01/2021
  !
  use Input_Aux_module
  use String_module

  implicit none

  class(pm_well_type) :: pm_well
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  character(len=MAXWORDLENGTH) :: word

  error_string = trim(error_string) // ',WELL_BOUNDARY_CONDITIONS'
  found = PETSC_TRUE

  select case(trim(keyword))
  !-------------------------------------
    case('WELL_BOUNDARY_CONDITIONS')
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
          case('BOTTOM_OF_HOLE')
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
                    case('LIQUID_PRESSURE')
                      call InputReadDouble(input,option,pm_well%well%bh_p)
                      call InputReadAndConvertUnits(input,pm_well%well%bh_p, &
                           'Pa','WELL_BOUNDARY_CONDITIONS,BOTTOM_OF_HOLE,&
                           &LIQUID_PRESSURE',option)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('GAS_SATURATION')
                      call InputReadDouble(input,option,pm_well%well%bh_sg)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('PRESSURE_SET_BY_RESERVOIR')
                      pm_well%well%bh_p_set_by_reservoir = PETSC_TRUE
                  !-----------------------------
                    case('SATURATION_SET_BY_RESERVOIR')
                      pm_well%well%bh_sg_set_by_reservoir = PETSC_TRUE
                  !-----------------------------
                    case('LIQUID_VELOCITY')
                      call InputReadDouble(input,option,pm_well%well%bh_ql)
                      call InputReadAndConvertUnits(input,pm_well%well%bh_ql, &
                           'm/s','WELL_BOUNDARY_CONDITIONS,BOTTOM_OF_HOLE,&
                           &LIQUID_VELOCITY',option)
                  !-----------------------------
                    case('LIQUID_MASS_RATE')
                      call InputReadDouble(input,option,pm_well%well%bh_ql)
                      call InputReadAndConvertUnits(input,pm_well%well%bh_ql, &
                           'kg/s','WELL_BOUNDARY_CONDITIONS,BOTTOM_OF_HOLE,&
                           &LIQUID_MASS_RATE',option)
                  !-----------------------------
                    case('GAS_VELOCITY')
                      call InputReadDouble(input,option,pm_well%well%bh_qg)
                      call InputReadAndConvertUnits(input,pm_well%well%bh_qg, &
                           'm/s','WELL_BOUNDARY_CONDITIONS,BOTTOM_OF_HOLE,&
                           &GAS_VELOCITY',option)
                  !-----------------------------
                    case('GAS_MASS_RATE')
                      call InputReadDouble(input,option,pm_well%well%bh_qg)
                      call InputReadAndConvertUnits(input,pm_well%well%bh_qg, &
                           'kg/s','WELL_BOUNDARY_CONDITIONS,BOTTOM_OF_HOLE,&
                           &GAS_VELOCITY',option)
                  !-----------------------------
                    case default
                      call InputKeywordUnrecognized(input,word, &
                                                    error_string,option)
                  !-----------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case('TOP_OF_HOLE')
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
                    case('LIQUID_PRESSURE')
                      call InputReadDouble(input,option,pm_well%well%th_p)
                      call InputReadAndConvertUnits(input,pm_well%well%th_p, &
                           'Pa','WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE,&
                           &LIQUID_PRESSURE',option)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('GAS_SATURATION')
                      call InputReadDouble(input,option,pm_well%well%th_sg)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('LIQUID_VELOCITY')
                      call InputReadDouble(input,option,pm_well%well%th_ql)
                      call InputReadAndConvertUnits(input,pm_well%well%th_ql, &
                           'm/s','WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE,&
                           &LIQUID_VELOCITY',option)
                  !-----------------------------
                    case('LIQUID_MASS_RATE')
                      call InputReadDouble(input,option,pm_well%well%th_ql)
                      call InputReadAndConvertUnits(input,pm_well%well%th_ql, &
                           'kg/s','WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE,&
                           &LIQUID_MASS_RATE',option)
                  !-----------------------------
                    case('GAS_VELOCITY')
                      call InputReadDouble(input,option,pm_well%well%th_qg)
                      call InputReadAndConvertUnits(input,pm_well%well%th_qg, &
                           'm/s','WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE,&
                           &GAS_VELOCITY',option)
                  !-----------------------------
                    case('GAS_MASS_RATE')
                      call InputReadDouble(input,option,pm_well%well%th_qg)
                      call InputReadAndConvertUnits(input,pm_well%well%th_qg, &
                           'kg/s','WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE,&
                           &GAS_MASS_RATE',option)
                  !-----------------------------
                    case('TRANSPORT_CONDITION')
                      call InputReadWord(input,option, &
                                   pm_well%well%tran_condition_name,PETSC_TRUE)
                      call InputErrorMsg(input,option, &
                                    'TRANSPORT_CONDITION name',error_string)
                  !-----------------------------
                    case default
                      call InputKeywordUnrecognized(input,word, &
                                                    error_string,option)
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

      ! ----------------- error messaging -------------------------------------
      if (Initialized(pm_well%well%bh_p) .and. &
          pm_well%well%bh_p_set_by_reservoir) then
        option%io_buffer = 'Either keyword BOTTOM_OF_HOLE,PRESSURE or keyword &
          &BOTTOM_OF_HOLE,PRESSURE_SET_BY_RESERVOIR must be provided in the ' &
          // trim(error_string) // ' block, but NOT BOTH.'
        call PrintErrMsg(option)
      endif

      if (.not. (Initialized(pm_well%well%bh_ql) .and. &
                 Initialized(pm_well%well%bh_qg))) then
        if ( .not. ((Initialized(pm_well%well%bh_p) .or. &
                     pm_well%well%bh_p_set_by_reservoir) .and. &
                    (Initialized(pm_well%well%bh_sg) .or. &
                     pm_well%well%bh_sg_set_by_reservoir))) then

          option%io_buffer ='WELL_MODEL_WIPP_DARCY well model needs both Dirichlet &
             &LIQUID_PRESSURE and GAS_SATURATION set in the ' &
             // trim(error_string) // ' block.'
          call PrintErrMsg(option)
        endif
      endif

      if ((Initialized(pm_well%well%th_p).or.Initialized(pm_well%well%th_sg)) .and. &
          .not.(Initialized(pm_well%well%th_p).and.Initialized(pm_well%well%th_sg))) &
          then
        option%io_buffer ='WIPP_DARCY well model needs both Dirichlet &
           &LIQUID_PRESSURE and GAS_SATURATION set in the ' &
           // trim(error_string) // ' block.'
        call PrintErrMsg(option)
      endif

      if ((pm_well%option%itranmode /= NULL_MODE) .and. &
          (trim(pm_well%well%tran_condition_name) == '')) then
        option%io_buffer ='A TRANSPORT_CONDITION was not provided in the &
          &WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE block.'
        call PrintErrMsg(option)
      endif

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (Initialized(pm_well%well%bh_p)) pm_well%flow_soln%bh_p = PETSC_TRUE
  if (pm_well%well%bh_p_set_by_reservoir) pm_well%flow_soln%bh_p = PETSC_TRUE
  if (Initialized(pm_well%well%th_p)) pm_well%flow_soln%th_p = PETSC_TRUE

  if (Initialized(pm_well%well%bh_ql)) pm_well%flow_soln%bh_q = PETSC_TRUE
  if (Initialized(pm_well%well%th_ql)) pm_well%flow_soln%th_q = PETSC_TRUE

  end subroutine PMWellReadWellBCs

! ************************************************************************** !

subroutine PMWellReadFlowSolver(pm_well,input,option,keyword,error_string, &
                                found)
  !
  ! Reads input file parameters associated with the well model solver.
  !
  ! Author: Jenn Frederick
  ! Date: 02/03/2022
  !
  use Input_Aux_module
  use String_module

  implicit none

  class(pm_well_type) :: pm_well
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscInt :: num_errors
  character(len=MAXWORDLENGTH) :: word

  error_string = trim(error_string) // ',WELL_FLOW_SOLVER'
  found = PETSC_TRUE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('WELL_FLOW_SOLVER')
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
          case('MAXIMUM_NUMBER_OF_ITERATIONS')
            call InputReadInt(input,option,pm_well%flow_soln%max_iter)
            call InputErrorMsg(input,option,'MAXIMUM_NUMBER_OF_ITERATIONS', &
                               error_string)
        !-----------------------------
          case('MAXIMUM_NUMBER_OF_TS_CUTS')
            call InputReadInt(input,option,pm_well%flow_soln%max_ts_cut)
            call InputErrorMsg(input,option,'MAXIMUM_NUMBER_OF_TS_CUTS', &
                               error_string)
        !-----------------------------
          case('TIMESTEP_CUT_FACTOR')
            call InputReadDouble(input,option,pm_well%flow_soln%ts_cut_factor)
            call InputErrorMsg(input,option,'TIMESTEP_CUT_FACTOR', &
                               error_string)
        !-----------------------------
          case('TIMESTEP_RAMP_FACTOR')
            call InputReadDouble(input,option,pm_well%flow_soln%ts_ramp_factor)
            call InputErrorMsg(input,option,'TIMESTEP_RAMP_FACTOR', &
                               error_string)
        !-----------------------------
          case('ITOL_ABSOLUTE_RESIDUAL')
            call InputReadDouble(input,option,pm_well%flow_soln%itol_abs_res)
            call InputErrorMsg(input,option,'ITOL_ABSOLUTE_RESIDUAL', &
                               error_string)
        !-----------------------------
          case('ITOL_SCALED_RESIDUAL')
            call InputReadDouble(input,option,pm_well%flow_soln%itol_scaled_res)
            call InputErrorMsg(input,option,'ITOL_SCALED_RESIDUAL', &
                               error_string)
        !-----------------------------
          case('ITOL_ABS_UPDATE_PRESSURE')
            call InputReadDouble(input,option,pm_well%flow_soln%itol_abs_update_p)
            call InputErrorMsg(input,option,'ITOL_ABS_UPDATE_PRESSURE', &
                               error_string)
        !-----------------------------
          case('ITOL_ABS_UPDATE_SATURATION')
            call InputReadDouble(input,option,pm_well%flow_soln%itol_abs_update_s)
            call InputErrorMsg(input,option,'ITOL_ABS_UPDATE_SATURATION', &
                               error_string)
        !-----------------------------
          case('ITOL_REL_UPDATE_PRESSURE')
            call InputReadDouble(input,option,pm_well%flow_soln%itol_rel_update_p)
            call InputErrorMsg(input,option,'ITOL_REL_UPDATE_PRESSURE', &
                               error_string)
        !-----------------------------
          case('ITOL_REL_UPDATE_SATURATION')
            call InputReadDouble(input,option,pm_well%flow_soln%itol_rel_update_s)
            call InputErrorMsg(input,option,'ITOL_REL_UPDATE_SATURATION', &
                               error_string)
        !-----------------------------
          case('MIN_FLOW_DT_SCALE')
            call InputReadDouble(input,option,min_flow_dt_scale)
            call InputErrorMsg(input,option,'MIN_FLOW_DT_SCALE', &
                               error_string)
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      ! ----------------- error messaging -------------------------------------


  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the WELLBORE_MODEL,WELL_FLOW_SOLVER block. See above &
                       &error messages.'
    call PrintErrMsg(option)
  endif

  end subroutine PMWellReadFlowSolver

! ************************************************************************** !

subroutine PMWellReadTranSolver(pm_well,input,option,keyword,error_string, &
                                found)
  !
  ! Reads input file parameters associated with the well model solver.
  !
  ! Author: Jenn Frederick
  ! Date: 02/17/2022
  !
  use Input_Aux_module
  use String_module

  implicit none

  class(pm_well_type) :: pm_well
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscInt :: num_errors
  character(len=MAXWORDLENGTH) :: word

  error_string = trim(error_string) // ',WELL_TRANSPORT_SOLVER'
  found = PETSC_TRUE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('WELL_TRANSPORT_SOLVER')
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
          case('MAXIMUM_NUMBER_OF_ITERATIONS')
            call InputReadInt(input,option,pm_well%tran_soln%max_iter)
            call InputErrorMsg(input,option,'MAXIMUM_NUMBER_OF_ITERATIONS', &
                               error_string)
        !-----------------------------
          case('MAXIMUM_NUMBER_OF_TS_CUTS')
            call InputReadInt(input,option,pm_well%tran_soln%max_ts_cut)
            call InputErrorMsg(input,option,'MAXIMUM_NUMBER_OF_TS_CUTS', &
                               error_string)
        !-----------------------------
          case('TIMESTEP_CUT_FACTOR')
            call InputReadDouble(input,option,pm_well%tran_soln%ts_cut_factor)
            call InputErrorMsg(input,option,'TIMESTEP_CUT_FACTOR', &
                               error_string)
        !-----------------------------
          case('TIMESTEP_RAMP_FACTOR')
            call InputReadDouble(input,option,pm_well%tran_soln%ts_ramp_factor)
            call InputErrorMsg(input,option,'TIMESTEP_RAMP_FACTOR', &
                               error_string)
        !-----------------------------
          case('ITOL_ABSOLUTE_RESIDUAL')
            call InputReadDouble(input,option,pm_well%tran_soln%itol_abs_res)
            call InputErrorMsg(input,option,'ITOL_ABSOLUTE_RESIDUAL', &
                               error_string)
        !-----------------------------
          case('ITOL_SCALED_RESIDUAL')
            call InputReadDouble(input,option, &
                                 pm_well%tran_soln%itol_scaled_res)
            call InputErrorMsg(input,option,'ITOL_SCALED_RESIDUAL', &
                               error_string)
        !-----------------------------
          case('ITOL_ABSOLUTE_UPDATE')
            call InputReadDouble(input,option, &
                                 pm_well%tran_soln%itol_abs_update)
            call InputErrorMsg(input,option,'ITOL_ABSOLUTE_UPDATE', &
                               error_string)
        !-----------------------------
          case('ITOL_RELATIVE_UPDATE')
            call InputReadDouble(input,option, &
                                 pm_well%tran_soln%itol_rel_update)
            call InputErrorMsg(input,option,'ITOL_RELATIVE_UPDATE', &
                               error_string)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      ! ----------------- error messaging -------------------------------------


  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the WELLBORE_MODEL,WELL_TRANSPORT_SOLVER block. &
                       &See above error messages.'
    call PrintErrMsg(option)
  endif

  end subroutine PMWellReadTranSolver

! ************************************************************************** !

subroutine PMWellReadWellModelType(pm_well,input,option,keyword,error_string, &
                                   found)
  !
  ! Reads input file parameters associated with the well model type.
  !
  ! Author: Michael Nole
  ! Date: 12/22/2021
  !
  use Input_Aux_module
  use String_module

  implicit none

  class(pm_well_type) :: pm_well
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  character(len=MAXWORDLENGTH) :: word

  error_string = trim(error_string) // ',WELL_MODEL_TYPE'
  found = PETSC_TRUE

  select case(trim(keyword))
  !-------------------------------------
    case('WELL_MODEL_TYPE')
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
          case('HYDROSTATIC')
            pm_well%well%well_model_type = WELL_MODEL_HYDROSTATIC
        !-----------------------------
          case('WIPP_DARCY')
            pm_well%well%well_model_type = WELL_MODEL_WIPP_DARCY
        !-----------------------------
          case('FULL_MOMENTUM')
            pm_well%well%well_model_type = WELL_MODEL_FULL_MOMENTUM
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

end subroutine PMWellReadWellModelType

! ************************************************************************** !

subroutine PMWellReadWellConstraintType(pm_well,input,option,keyword, &
                                        error_string,found)
  !
  ! Reads input file parameters associated with the well constraint type.
  !
  ! Author: Michael Nole
  ! Date: 12/22/2021
  !
  use Input_Aux_module
  use String_module

  implicit none

  class(pm_well_type) :: pm_well
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  character(len=MAXWORDLENGTH) :: word

  error_string = trim(error_string) // ',WELL_CONSTRAINT_TYPE'
  found = PETSC_TRUE

  select case(trim(keyword))
  !-------------------------------------
    case('WELL_CONSTRAINT_TYPE')
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
          case('CONSTANT_PRESSURE')
            pm_well%well%well_constraint_type = CONSTANT_PRESSURE
        !-----------------------------
          case('CONSTANT_PRESSURE_HYDROSTATIC')
            pm_well%well%well_constraint_type = CONSTANT_PRESSURE_HYDROSTATIC
        !-----------------------------
          case('CONSTANT_RATE')
            pm_well%well%well_constraint_type = CONSTANT_RATE
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

end subroutine PMWellReadWellConstraintType

! ************************************************************************** !

subroutine PMWellReadWellOutput(pm_well,input,option,keyword,error_string, &
                                found)
  !
  ! Reads input file parameters associated with the well model type.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 09/29/2022
  !
  use Input_Aux_module
  use String_module

  implicit none

  class(pm_well_type) :: pm_well
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  character(len=MAXWORDLENGTH) :: word

  error_string = trim(error_string) // ',WELL_MODEL_OUTPUT'
  found = PETSC_TRUE

  select case(trim(keyword))
  !-------------------------------------
    case('WELL_MODEL_OUTPUT')
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
          case('WELL_LIQ_PRESSURE')
            pm_well%well%output_pl = PETSC_TRUE
        !-----------------------------
          case('WELL_GAS_PRESSURE')
            pm_well%well%output_pg = PETSC_TRUE
        !-----------------------------
          case('WELL_LIQ_SATURATION')
            pm_well%well%output_sl = PETSC_TRUE
        !-----------------------------
          case('WELL_GAS_SATURATION')
            pm_well%well%output_sg = PETSC_TRUE
        !-----------------------------
          case('WELL_AQ_CONC')
            pm_well%well%output_aqc = PETSC_TRUE
        !-----------------------------
          case('WELL_AQ_MASS')
            pm_well%well%output_aqm = PETSC_TRUE
        !-----------------------------
          case('WELL_LIQ_Q')
            pm_well%well%liq%output_Q = PETSC_TRUE
        !-----------------------------
          case('WELL_GAS_Q')
            pm_well%well%gas%output_Q = PETSC_TRUE
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

end subroutine PMWellReadWellOutput

! ************************************************************************** !

subroutine PMWellReadPass2(input,option)
  !
  ! Reads input file parameters associated with the well process model.
  ! Second pass dummy read.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 11/03/2021

  use Input_Aux_module
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: card

  error_string = 'SUBSURFACE,WELLBORE_MODEL'

  input%ierr = 0
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,card)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      !--------------------
      case('WELL_GRID','WELL','WELL_MODEL_TYPE','WELL_CONSTRAINT_TYPE', &
           'WELL_FLOW_SOLVER','WELL_TRANSPORT_SOLVER','WELL_MODEL_OUTPUT')
        call InputSkipToEND(input,option,card)
      !--------------------
      case('WELL_BOUNDARY_CONDITIONS')
        call InputSkipToEND(input,option,card)
        call InputSkipToEND(input,option,card)
        call InputSkipToEND(input,option,card)
      !--------------------
    end select
  enddo
  call InputPopBlock(input,option)
  !call InputSkipToEND(input,option,card)

end subroutine PMWellReadPass2

! ************************************************************************** !

subroutine PMWellSetRealization(this,realization)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  use Realization_Subsurface_class

  implicit none

  class(pm_well_type) :: this
  class(realization_subsurface_type), pointer :: realization

  this%realization => realization
  this%realization_base => realization

end subroutine PMWellSetRealization

! ************************************************************************** !

recursive subroutine PMWellInitializeRun(this)
  !
  ! Initializes the simulation run for the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  use Condition_module
  use Strata_module
  use Output_Aux_module
  use Option_module
  use Transport_Constraint_Base_module
  use Transport_Constraint_NWT_module
  use Transport_Constraint_module

  implicit none

  class(pm_well_type) :: this

  type(strata_type), pointer :: patch_strata, well_strata
  type(radioactive_decay_rxn_type), pointer :: rad_rxn
  type(species_type), pointer :: species
  type(tran_condition_type), pointer :: tran_condition
  class(tran_constraint_base_type), pointer :: cur_constraint
  type(output_variable_list_type), pointer :: output_var_list
  type(option_type), pointer :: option
  PetscInt :: nsegments, k
  PetscReal :: curr_time

  option => this%option

  curr_time = option%time
  nsegments = this%well_grid%nsegments

  ! srcsink_water/gas is indexed (0,:) unperturbed value
  !                              (1,:) perturbed wrt gas pressure
  !                              (2,:) perturbed wrt gas saturation
  allocate(this%srcsink_water(0:2,nsegments))
  allocate(this%srcsink_gas(0:2,nsegments))

  patch_strata => this%realization%patch%strata_list%first
  do
    if (.not.associated(patch_strata)) exit
    if (patch_strata%well) then
      well_strata => StrataCreateFromStrata(patch_strata)
      well_strata%region => patch_strata%region
      well_strata%material_property => patch_strata%material_property
      call StrataAddToList(well_strata,this%strata_list)
    endif
    patch_strata => patch_strata%next
  enddo

  call PMWellUpdateStrata(this,curr_time)

  call PMWellInitWellVars(this%well,this%well_grid,this%transport, &
                          nsegments,this%nspecies)
  call PMWellInitFluidVars(this%well,nsegments, &
                           this%option%flow%reference_density,option)

  do k = 1,option%nflowdof
    allocate(this%well_pert(k)%diameter(nsegments))
    allocate(this%well_pert(k)%WI_base(nsegments))
    allocate(this%well_pert(k)% &
             friction_factor(size(this%well%friction_factor)))
    allocate(this%well_pert(k)%r0(size(this%well%r0)))
    allocate(this%well_pert(k)%skin(size(this%well%skin)))
    this%well_pert(k)%diameter = this%well%diameter
    this%well_pert(k)%WI_base = this%well%WI_base
    this%well_pert(k)%friction_factor = this%well%friction_factor
    this%well_pert(k)%well_model_type = this%well%well_model_type
    call PMWellInitWellVars(this%well_pert(k),this%well_grid, &
                          this%transport,nsegments,this%nspecies)
    call PMWellInitFluidVars(this%well_pert(k),nsegments, &
                           this%option%flow%reference_density,option)
    call PMWellInitRes(this%well_pert(k)%reservoir,nsegments,this%nspecies,option)
    call PMWellInitRes(this%well_pert(k)%reservoir_save,nsegments,this%nspecies,option)
  enddo

  call PMWellInitRes(this%well%reservoir,nsegments,this%nspecies,option)
  call PMWellInitRes(this%well%reservoir_save,nsegments,this%nspecies,option)
  call PMWellInitFlowSoln(this%flow_soln,this%nphase,nsegments)
  if (this%transport) then
    call PMWellInitTranSoln(this%tran_soln,this%nspecies,nsegments)
  endif


  ! Initialize perturbations
  allocate(this%pert(nsegments,this%nphase))
  this%pert = 0.d0


  if (this%transport) then
    allocate(this%well%reservoir%aqueous_conc(this%nspecies, &
                                         this%well_grid%nsegments))
    allocate(this%well%reservoir%aqueous_mass(this%nspecies, &
                                         this%well_grid%nsegments))

    allocate(this%well%species_names(this%nspecies))
    allocate(this%well%species_parent_id(this%nspecies))
    allocate(this%well%species_radioactive(this%nspecies))
    allocate(this%well%species_decay_rate(this%nspecies))
    allocate(this%well%species_parent_decay_rate(this%nspecies))
    this%well%species_names(:) = ''
    this%well%species_decay_rate(:) = 0.d0
    this%well%species_parent_decay_rate(:) = 0.d0
    this%well%species_parent_id(:) = 0
    k = 0
    species => this%realization%reaction_nw%species_list
    do
      if (.not.associated(species)) exit
      k = k + 1
      this%well%species_names(k) = species%name
      this%well%species_radioactive(k) = species%radioactive
      if (species%radioactive) then
        ! Find the reaction object associated with this species
        rad_rxn => this%realization%reaction_nw%rad_decay_rxn_list
        do
          if (.not.associated(rad_rxn)) exit
          if (rad_rxn%species_id == species%id) exit
          rad_rxn => rad_rxn%next
        enddo
        this%well%species_decay_rate(k) = rad_rxn%rate_constant
        if (rad_rxn%parent_id > 0.d0) then
          this%well%species_parent_id(k) = rad_rxn%parent_id
          ! Find the reaction object associated with the parent species
          rad_rxn => this%realization%reaction_nw%rad_decay_rxn_list
          do
            if (.not.associated(rad_rxn)) exit
            if (rad_rxn%species_id == this%well%species_parent_id(k)) exit
            rad_rxn => rad_rxn%next
          enddo
          this%well%species_parent_decay_rate(k) = rad_rxn%rate_constant
        endif
      endif
      species => species%next
    enddo

    tran_condition => this%realization%transport_conditions%first
    do
      if (.not.associated(tran_condition)) exit
        if (trim(tran_condition%name) == &
            trim(this%well%tran_condition_name)) exit
      tran_condition => tran_condition%next
    enddo
    if (.not.associated(tran_condition)) then
      this%option%io_buffer = 'TRANSPORT_CONDITION ' // &
        trim(this%well%tran_condition_name) // ' for WELLBORE_MODEL,&
        &WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE not found.'
      call PrintErrMsg(this%option)
    endif
    cur_constraint => tran_condition%cur_constraint_coupler%constraint
    select type(constraint=>cur_constraint)
      class is (tran_constraint_nwt_type)
        if (any(constraint%nwt_species%constraint_type /=  &
            CONSTRAINT_AQ_EQUILIBRIUM)) then
          this%option%io_buffer = 'TRANSPORT_CONDITION ' // &
            trim(this%well%tran_condition_name) // ' for WELLBORE_MODEL,&
            &WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE CONSTRAINT must be of &
            &type "AQ".'
          call PrintErrMsg(this%option)
        endif
        this%well%aqueous_conc_th = constraint%nwt_species%constraint_conc
    end select
  endif

  ! Setup the output variables for snapshot and observation files
  output_var_list => this%realization%output_option%output_snap_variable_list
  call PMWellSetPlotVariables(output_var_list,this)
  if (.not.associated( &
                this%realization%output_option%output_snap_variable_list, &
                this%realization%output_option%output_obs_variable_list)) then
    output_var_list => this%realization%output_option%output_obs_variable_list
    call PMWellSetPlotVariables(output_var_list,this)
  endif

  if (this%print_well) then
    call PMWellOutputHeader(this)
  endif

end subroutine PMWellInitializeRun

! ************************************************************************** !

subroutine PMWellInitWellVars(well,well_grid,with_transport,nsegments,nspecies)
  !
  ! Initializes well variables.
  !
  ! Author: Michael Nole
  ! Date: 03/06/2023

  use Option_module

  implicit none

  type(well_type) :: well
  type(well_grid_type), pointer :: well_grid
  PetscBool :: with_transport
  PetscInt :: nsegments
  PetscInt :: nspecies

  allocate(well%WI(nsegments))
  allocate(well%pl(nsegments))
  allocate(well%pg(nsegments))
  allocate(well%temp(nsegments))
  allocate(well%ql(nsegments-1))
  allocate(well%qg(nsegments-1))
  allocate(well%ql_bc(2))
  allocate(well%qg_bc(2))
  allocate(well%ql_kmol(nsegments-1))
  allocate(well%qg_kmol(nsegments-1))
  allocate(well%ql_kmol_bc(2))
  allocate(well%qg_kmol_bc(2))
  allocate(well%liq_cum_mass(nsegments))
  allocate(well%liq_mass(nsegments))
  well%liq_cum_mass = 0.d0
  well%liq_mass = 0.d0

  allocate(well%ccid(nsegments))
  allocate(well%permeability(nsegments))
  allocate(well%phi(nsegments))

  allocate(well%area(nsegments))
  well%area = 3.14159*(well%diameter/2.d0)*(well%diameter/2.d0)

  allocate(well%volume(nsegments))
  well%volume = well%area*well_grid%dh

  allocate(well%mass_balance_liq(nsegments))

  if (with_transport) then
    allocate(well%aqueous_conc(nspecies,nsegments))
    allocate(well%aqueous_mass(nspecies,nsegments))
    allocate(well%aqueous_conc_th(nspecies))
  endif

end subroutine PMWellInitWellVars

! ************************************************************************** !

subroutine PMWellInitFluidVars(well,nsegments,reference_density,option)
  !
  ! Initializes fluid variables.
  !
  ! Author: Michael Nole
  ! Date: 03/06/2023

  use Option_module

  implicit none

  type(well_type) :: well
  PetscInt :: nsegments
  PetscReal :: reference_density(MAX_PHASE)
  type(option_type), pointer :: option

  allocate(well%liq%s(nsegments))
  well%liq%den0 = reference_density(1)
  allocate(well%liq%den(nsegments))
  well%liq%den(:) = well%liq%den0
  allocate(well%liq%visc(nsegments))
  allocate(well%liq%Q(nsegments))
  allocate(well%liq%kr(nsegments))
  allocate(well%liq%xmass(nsegments,option%nflowspec))

  allocate(well%gas%s(nsegments))
  well%gas%den0 = reference_density(2)
  allocate(well%gas%den(nsegments))
  well%gas%den(:) = well%gas%den0
  allocate(well%gas%visc(nsegments))
  allocate(well%gas%Q(nsegments))
  allocate(well%gas%kr(nsegments))
  allocate(well%gas%xmass(nsegments,option%nflowspec))

end subroutine PMWellInitFluidVars

! ************************************************************************** !

subroutine PMWellInitRes(reservoir,nsegments,idof,option)
  !
  ! Initializes well variables.
  !
  ! Author: Michael Nole
  ! Date: 03/06/2023

  use Option_module

  implicit none

  type(well_reservoir_type) :: reservoir
  PetscInt :: nsegments, idof
  type(option_type), pointer :: option

  allocate(reservoir%p_l(nsegments))
  allocate(reservoir%p_g(nsegments))
  allocate(reservoir%s_l(nsegments))
  allocate(reservoir%s_g(nsegments))
  allocate(reservoir%temp(nsegments))
  allocate(reservoir%mobility_l(nsegments))
  allocate(reservoir%mobility_g(nsegments))
  allocate(reservoir%kr_l(nsegments))
  allocate(reservoir%kr_g(nsegments))
  allocate(reservoir%den_l(nsegments))
  allocate(reservoir%den_g(nsegments))
  allocate(reservoir%visc_l(nsegments))
  allocate(reservoir%visc_g(nsegments))
  allocate(reservoir%e_por(nsegments))
  allocate(reservoir%xmass_liq(nsegments,option%nflowspec))
  allocate(reservoir%xmass_gas(nsegments,option%nflowspec))
  allocate(reservoir%kx(nsegments))
  allocate(reservoir%ky(nsegments))
  allocate(reservoir%kz(nsegments))
  allocate(reservoir%dx(nsegments))
  allocate(reservoir%dy(nsegments))
  allocate(reservoir%dz(nsegments))
  allocate(reservoir%volume(nsegments))

  allocate(reservoir%aqueous_conc(idof,nsegments))
  allocate(reservoir%aqueous_mass(idof,nsegments))

end subroutine PMWellInitRes

! ************************************************************************** !

subroutine PMWellInitFlowSoln(flow_soln,nphase,nsegments)
  !
  ! Initializes flow solution.
  !
  ! Author: Michael Nole
  ! Date: 03/06/2023

  type(well_soln_flow_type) :: flow_soln
  PetscInt :: nphase
  PetscInt :: nsegments

  PetscInt :: i,j

  allocate(flow_soln%residual(nsegments*flow_soln%ndof))
  allocate(flow_soln%update(nsegments*flow_soln%ndof))
  flow_soln%residual(:) = UNINITIALIZED_DOUBLE
  flow_soln%update(:) = UNINITIALIZED_DOUBLE

  allocate(flow_soln%Jacobian(nphase*nsegments,nphase*nsegments))
  do i = 1,nphase*nsegments
    do j = 1,nphase*nsegments
      flow_soln%Jacobian(i,j) = UNINITIALIZED_DOUBLE
    enddo
  enddo

  allocate(flow_soln%prev_soln%pl(nsegments))
  allocate(flow_soln%prev_soln%sg(nsegments))
  allocate(flow_soln%soln_save%pl(nsegments))
  allocate(flow_soln%soln_save%sg(nsegments))

end subroutine PMWellInitFlowSoln

! ************************************************************************** !

subroutine PMWellInitTranSoln(tran_soln,nspecies,nsegments)
  !
  ! Initializes transport solution.
  !
  ! Author: Michael Nole
  ! Date: 03/06/2023

  type(well_soln_tran_type) :: tran_soln
  PetscInt :: nspecies
  PetscInt :: nsegments


  allocate(tran_soln%residual(nsegments*tran_soln%ndof))
  allocate(tran_soln%update(nsegments*tran_soln%ndof))
  tran_soln%residual(:) = UNINITIALIZED_DOUBLE
  tran_soln%update(:) = UNINITIALIZED_DOUBLE

  allocate(tran_soln%Jacobian(nspecies*nsegments, &
                              nspecies*nsegments))
  tran_soln%Jacobian = UNINITIALIZED_DOUBLE

  allocate(tran_soln%prev_soln%aqueous_conc(nspecies,nsegments))
  allocate(tran_soln%prev_soln%aqueous_mass(nspecies,nsegments))

end subroutine PMWellInitTranSoln

! ************************************************************************** !

recursive subroutine PMWellFinalizeRun(this)
  !
  ! Finalizes the simulation run for the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this

  ! placeholder

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMWellFinalizeRun

! ************************************************************************** !

subroutine PMWellInitializeTimestep(this)
  !
  ! Initializes and takes the time step for the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this

  PetscReal :: curr_time

  curr_time = this%option%time - this%option%flow_dt
  this%dt_flow = this%realization%option%flow_dt

  if (Initialized(this%intrusion_time_start) .and. &
      (curr_time < this%intrusion_time_start) .and. &
      .not. this%well_on) return
  
  call PMWellInitializeTimestepFlow(this,curr_time)

end subroutine PMWellInitializeTimestep

! ************************************************************************** !

subroutine PMWellInitializeTimestepFlow(pm_well,curr_time)
  !
  ! Initializes and takes the time step for the well process model - flow.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 05/11/2023

  use Option_module

  implicit none

  class(pm_well_type) :: pm_well
  PetscReal :: curr_time
  PetscInt :: i

  ! update the reservoir object with current reservoir properties
  call PMWellCopyReservoir(pm_well%well%reservoir,pm_well%well%reservoir_save, &
                           pm_well%transport)
  select case (pm_well%option%iflowmode)
    case(WF_MODE)
      call PMWellUpdateReservoirWIPP(pm_well,-999)
    case(SCO2_MODE)
      call PMWellUpdateReservoirSCO2(pm_well,-999,-999)
  end select
  
  call PMWellComputeWellIndex(pm_well)

  call PMWellUpdateStrata(pm_well,curr_time)

  if (initialize_well_flow) then
    ! enter here if its the very first timestep
    call PMWellInitializeWellFlow(pm_well)
  endif

  do i = 1,pm_well%option%nflowdof
    call PMWellCopyWell(pm_well%well,pm_well%well_pert(i),pm_well%transport)
  enddo

  if (pm_well%well%bh_p_set_by_reservoir) then
    pm_well%well%bh_p = pm_well%well%reservoir%p_l(1)
  endif
  if (pm_well%well%bh_sg_set_by_reservoir) then
    pm_well%well%bh_sg = pm_well%well%reservoir%s_g(1)
  endif

  select case (pm_well%option%iflowmode)
    case(WF_MODE)
      call PMWellUpdatePropertiesWIPPFlow(pm_well,pm_well%well,&
                        pm_well%realization%patch%characteristic_curves_array,&
                        pm_well%realization%option)
    case(SCO2_MODE)
      call PMWellUpdatePropertiesSCO2Flow(pm_well,pm_well%well,&
                        pm_well%realization%option)
  end select
  

end subroutine PMWellInitializeTimestepFlow

! ************************************************************************** !

subroutine PMWellInitializeWellFlow(pm_well)
  !
  ! Initializes the well for the first time step for flow.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  use SCO2_Aux_module

  implicit none

  class(pm_well_type) :: pm_well

  type(well_reservoir_type), pointer :: reservoir
  type(strata_type), pointer :: strata
  type(sco2_auxvar_type), pointer :: sco2_auxvar
  type(option_type), pointer :: option
  PetscInt :: k, ghosted_id

  option => pm_well%option
  reservoir => pm_well%well%reservoir

  ! set initial flow parameters to the reservoir flow parameters
  pm_well%well%pl = reservoir%p_l
  pm_well%well%pg = reservoir%p_g
  pm_well%well%temp = reservoir%temp
  pm_well%well%liq%s = reservoir%s_l
  pm_well%well%gas%s = reservoir%s_g
  pm_well%well%liq%den = reservoir%den_l
  pm_well%well%gas%den = reservoir%den_g
  pm_well%well%liq%visc = reservoir%visc_l
  pm_well%well%gas%visc = reservoir%visc_g
  pm_well%well%liq%xmass = reservoir%xmass_liq
  pm_well%well%gas%xmass = reservoir%xmass_gas
  ! BHP can be a flow primary variable if fully coupled to flow.
  if (option%iflowmode == SCO2_MODE) then
    if (pm_well%well_grid%h_rank_id(1) == option%myrank) then
        ghosted_id = pm_well%well_grid%h_ghosted_id(1)
        sco2_auxvar => &
          pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)
        pm_well%well%bh_p = sco2_auxvar%pres(option%gas_phase) - &
                            pm_well%well%liq%den(1) * &
                            pm_well%option%gravity(Z_DIRECTION) * &
                            pm_well%well_grid%dh(1)/2.d0
        
    endif
  endif
  ! update the Darcy fluxes within the well
  do k = 1,option%nflowdof

    pm_well%well_pert(k)%pl = reservoir%p_l
    pm_well%well_pert(k)%pg = reservoir%p_g
    pm_well%well_pert(k)%temp = reservoir%temp
    pm_well%well_pert(k)%liq%s = reservoir%s_l
    pm_well%well_pert(k)%gas%s = reservoir%s_g
    pm_well%well_pert(k)%liq%den = reservoir%den_l
    pm_well%well_pert(k)%gas%den = reservoir%den_g
    pm_well%well_pert(k)%liq%visc = reservoir%visc_l
    pm_well%well_pert(k)%gas%visc = reservoir%visc_g
    pm_well%well_pert(k)%liq%xmass = reservoir%xmass_liq
    pm_well%well_pert(k)%gas%xmass = reservoir%xmass_gas
    pm_well%well_pert(k)%bh_p = pm_well%well%bh_p

  enddo

  pm_well%flow_soln%prev_soln%pl = pm_well%well%pl
  pm_well%flow_soln%prev_soln%sg = pm_well%well%gas%s
  pm_well%flow_soln%soln_save%pl = pm_well%well%pl
  pm_well%flow_soln%soln_save%sg = pm_well%well%gas%s
  pm_well%flow_soln%prev_soln%bh_p = pm_well%well%bh_p
  pm_well%flow_soln%soln_save%bh_p = pm_well%well%bh_p

  ! Link well material properties
  if (pm_well%well_comm%comm /= MPI_COMM_NULL) then
    do k = 1,pm_well%well_grid%nsegments
      strata => pm_well%strata_list%first
      do
        if (.not.associated(strata)) exit
        if (strata%id == pm_well%well_grid%strata_id(k)) then
          pm_well%well%ccid(k) = strata%material_property%saturation_function_id
          pm_well%well%permeability(k) = strata%material_property%permeability(3,3)
          pm_well%well%phi(k) = strata%material_property%porosity
          exit
        endif
        strata => strata%next
      enddo
    enddo
  endif

  call PMWellComputeWellIndex(pm_well)

  do k = 1,option%nflowdof
    pm_well%well_pert(k)%permeability(:) = pm_well%well%permeability(:)
    pm_well%well_pert(k)%phi(:) = pm_well%well%phi(:)
    pm_well%well_pert(k)%ccid(:) = pm_well%well%ccid(:)
    pm_well%well_pert(k)%WI(:) = pm_well%well%WI(:)
  enddo

  initialize_well_flow = PETSC_FALSE

end subroutine PMWellInitializeWellFlow

! ************************************************************************** !

subroutine PMWellInitializeWellTran(pm_well)
  !
  ! Initializes the well for the first time step for transport.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: k

  ! set initial transport parameters to the reservoir transport parameters
  if (pm_well%transport) then
    if (Initialized(pm_well%intrusion_time_start)) then
      ! set the borehole concentrations to the borehole zero value now
      do k = 1,pm_well%well_grid%nsegments
        pm_well%well%aqueous_mass(:,k) = pm_well%bh_zero_value*pm_well%well%volume(k)
        pm_well%well%aqueous_conc(:,k) = &
          pm_well%well%aqueous_mass(:,k) / &                           ! [mol]
          (pm_well%well%phi(k)*pm_well%well%volume(k)*pm_well%well%liq%s(k)) ! [m3-liq]
      enddo
    else
      ! set the wellbore concentrations to the reservoir values
      pm_well%well%aqueous_conc = pm_well%well%reservoir%aqueous_conc
      do k = 1,pm_well%well_grid%nsegments
        pm_well%well%aqueous_mass(:,k) = pm_well%well%aqueous_conc(:,k) * &
                pm_well%well%phi(k) * pm_well%well%volume(k) * pm_well%well%liq%s(k)
      enddo
    endif
    pm_well%tran_soln%prev_soln%aqueous_conc = pm_well%well%aqueous_conc
    pm_well%tran_soln%prev_soln%aqueous_mass = pm_well%well%aqueous_mass
  endif

  initialize_well_tran = PETSC_FALSE

end subroutine PMWellInitializeWellTran

! ************************************************************************** !

subroutine PMWellBasePerturb(this)
  !
  ! Base Perturb function, calls flow-mode specific perturbations.
  !
  ! Author: Michael Nole
  ! Date: 04/01/2024
  !

  use Option_module

  implicit none

  class(pm_well_type) :: this

  type(option_type), pointer :: option

  option => this%option

  select case(option%iflowmode)
    case(SCO2_MODE)
      call PMWellSCO2Perturb(this)
  end select

end subroutine PMWellBasePerturb

! ************************************************************************** !

subroutine PMWellSCO2Perturb(pm_well)
  !
  ! Perturb well variables when using SCO2 flow mode.
  !
  ! Author: Michael Nole
  ! Date: 04/01/2024
  !

  use Option_module
  use SCO2_Aux_module
  use Material_Aux_module
  use Grid_module

  implicit none

  class(pm_well_type) :: pm_well

  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: reservoir
  type(sco2_auxvar_type), pointer :: sco2_auxvar
  type(material_auxvar_type), pointer :: material_auxvar
  type(grid_type), pointer :: res_grid
  type(option_type), pointer :: option
  PetscInt :: idof, k
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  option => pm_well%option

  res_grid => pm_well%realization%patch%grid

  ! Go up the well: for each well segment that is on-process, update
  ! perturbed values for reservoir and well variables.

  ! Perturbed fluxes wrt reservoir variables
  do k = 1,pm_well%well_grid%nsegments
    if (pm_well%well_grid%h_rank_id(k) /= option%myrank) cycle

    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    do idof = 1,option%nflowdof
      well => pm_well%well_pert(idof)
      if (idof == SCO2_WELL_DOF) then
        sco2_auxvar => &
          pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)
        well%bh_p = pm_well%well%bh_p + pm_well%realization%patch%aux%sco2% &
                    auxvars(SCO2_WELL_DOF,ghosted_id)%pert
      else
        sco2_auxvar => &
          pm_well%realization%patch%aux%sco2%auxvars(idof,ghosted_id)
        well%bh_p = pm_well%well%bh_p
      endif
      reservoir => well%reservoir
      
      material_auxvar => &
        pm_well%realization%patch%aux%material%auxvars(ghosted_id)

      reservoir%p_l(k) = sco2_auxvar%pres(option%liquid_phase)
      reservoir%p_g(k) = sco2_auxvar%pres(option%gas_phase)
      reservoir%s_l(k) = sco2_auxvar%sat(option%liquid_phase)
      reservoir%s_g(k) = sco2_auxvar%sat(option%gas_phase)
      reservoir%temp(k) = sco2_auxvar%temp
      reservoir%mobility_l(k) = &
        sco2_auxvar%mobility(option%liquid_phase)
      reservoir%mobility_g(k) = sco2_auxvar%mobility(option%gas_phase)
      reservoir%kr_l(k) = sco2_auxvar%kr(option%liquid_phase)
      reservoir%kr_g(k) = sco2_auxvar%kr(option%gas_phase)
      reservoir%den_l(k) = sco2_auxvar%den_kg(option%liquid_phase)
      reservoir%den_g(k) = sco2_auxvar%den_kg(option%gas_phase)
      reservoir%visc_l(k) = sco2_auxvar%visc(option%liquid_phase)
      reservoir%visc_g(k) = sco2_auxvar%visc(option%gas_phase)
      reservoir%xmass_liq(k,:) = sco2_auxvar%xmass(:, &
                                                    option%liquid_phase)
      reservoir%xmass_gas(k,:) = sco2_auxvar%xmass(:, &
                                                    option%gas_phase)
      reservoir%e_por(k) = sco2_auxvar%effective_porosity

      reservoir%kx(k) = material_auxvar%permeability(1)
      reservoir%ky(k) = material_auxvar%permeability(2)
      reservoir%kz(k) = material_auxvar%permeability(3)
      reservoir%volume(k) = material_auxvar%volume

      if (res_grid%itype == STRUCTURED_GRID) then
        reservoir%dx(k) = res_grid%structured_grid%dx(ghosted_id)
        reservoir%dy(k) = res_grid%structured_grid%dy(ghosted_id)
        reservoir%dz(k) = res_grid%structured_grid%dz(ghosted_id)
      else
        reservoir%dz(k) = pm_well%well_grid%res_dz(k)
        reservoir%dx(k) = sqrt(material_auxvar%volume/ &
                                  reservoir%dz(k))
        reservoir%dy(k) = reservoir%dx(k)
      endif
    enddo
  enddo

  ! do k = 1,pm_well%well_grid%nsegments
  !   if (pm_well%well_grid%h_rank_id(k) /= option%myrank) cycle

  !   ghosted_id = pm_well%well_grid%h_ghosted_id(k)

  !   sco2_auxvar => &
  !         pm_well%realization%patch%aux%sco2%auxvars(SCO2_WELL_DOF,ghosted_id)
  !   well => pm_well%well_pert(SCO2_WELL_DOF)
  !   well%bh_p = sco2_auxvar%well%bh_p
  !   exit
  ! enddo

  ! Now update fluxes associated with perturbed values
  do idof = 1,option%nflowdof
    call PMWellSolveFlow(pm_well,idof,ierr)
  enddo

end subroutine PMWellSCO2Perturb

! ************************************************************************** !

subroutine PMWellUpdateStrata(pm_well,curr_time)
  !
  ! Updates the strata at current time for the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 11/23/2022

  use Strata_module

  implicit none

  class(pm_well_type) :: pm_well
  PetscReal :: curr_time

  type(strata_type), pointer :: well_strata
  PetscInt, pointer :: all_strata_id(:)
  PetscInt :: nsegments
  PetscInt :: k
  PetscErrorCode :: ierr

  nsegments = pm_well%well_grid%nsegments

  ! loop thru well stratas and mark them as active or inactive
  well_strata => pm_well%strata_list%first
  do
    if (.not.associated(well_strata)) exit
    if (Initialized(well_strata%start_time) .and. &
        Initialized(well_strata%final_time)) then
      if ((curr_time >= well_strata%start_time) .and. &
          (curr_time <= well_strata%final_time))  then
        well_strata%active = PETSC_TRUE
      else
        well_strata%active = PETSC_FALSE
      endif
    else
      well_strata%active = PETSC_TRUE
    endif
    well_strata => well_strata%next
  enddo

  ! update the well_grid%strata_id assignment for active strata
  do k = 1,nsegments
    if (pm_well%well_grid%h_rank_id(k) /= pm_well%option%myrank) cycle
    well_strata => pm_well%strata_list%first
    do
      if (.not.associated(well_strata)) exit
      ! not all materials are assigned to the regions from the start
      if (associated(well_strata%region%cell_ids)) then
          if ((any(well_strata%region%cell_ids == &
                   pm_well%well_grid%h_local_id(k))) .and. &
              (well_strata%active)) then
            pm_well%well_grid%strata_id(k) = well_strata%id
          endif
      endif    
      well_strata => well_strata%next
    enddo
  enddo

  allocate(all_strata_id(nsegments))
  all_strata_id = UNINITIALIZED_INTEGER
  if (pm_well%well_comm%comm /= MPI_COMM_NULL) then
    call MPI_Allreduce(pm_well%well_grid%strata_id,all_strata_id,nsegments, &
                MPI_INTEGER,MPI_MAX,pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
  endif
  pm_well%well_grid%strata_id = all_strata_id
  if (any(pm_well%well_grid%strata_id == UNINITIALIZED_INTEGER) .and. &
      any(pm_well%well_comm%petsc_rank_list == pm_well%option%myrank)) then
    pm_well%option%io_buffer =  'At least one WELLBORE_MODEL grid segment has &
        &not been assigned with a REGION and MATERIAL_PROPERTY with the use &
        &of the STRATA block.'
    call PrintErrMsg(pm_well%option)
  endif

end subroutine PMWellUpdateStrata

! ************************************************************************** !

subroutine PMWellUpdateReservoirWIPP(pm_well,wippflo_update_index)
  !
  ! Updates the reservoir properties for the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

  use WIPP_Flow_Aux_module
  use Material_Aux_module
  use NW_Transport_Aux_module
  use Grid_module

  implicit none

  class(pm_well_type) :: pm_well
  PetscInt :: wippflo_update_index

  type(well_reservoir_type), pointer :: reservoir
  type(wippflo_auxvar_type), pointer :: wippflo_auxvar
  type(nw_transport_auxvar_type), pointer :: nwt_auxvar
  type(material_auxvar_type), pointer :: material_auxvar
  type(grid_type), pointer :: res_grid
  type(option_type), pointer :: option
  type(well_comm_type), pointer :: well_comm
  PetscInt :: TAG, peer, root_rank
  PetscInt :: k, indx
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  option => pm_well%option
  well_comm => pm_well%well_comm 
  reservoir => pm_well%well%reservoir

  res_grid => pm_well%realization%patch%grid

  if (wippflo_update_index < ZERO_INTEGER) then
    indx = ZERO_INTEGER
  else
    indx = wippflo_update_index
  endif

  do k = 1,pm_well%well_grid%nsegments
    if (pm_well%well_grid%h_rank_id(k) /= option%myrank) cycle

    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    wippflo_auxvar => &
      pm_well%realization%patch%aux%wippflo%auxvars(indx,ghosted_id)
    if (pm_well%transport) then
      nwt_auxvar => &
        pm_well%realization%patch%aux%nwt%auxvars(ghosted_id)
    endif
    material_auxvar => &
      pm_well%realization%patch%aux%material%auxvars(ghosted_id)

    reservoir%p_l(k) = wippflo_auxvar%pres(option%liquid_phase)
    reservoir%p_g(k) = wippflo_auxvar%pres(option%gas_phase)
    reservoir%temp(k) = wippflo_auxvar%temp
    reservoir%s_l(k) = wippflo_auxvar%sat(option%liquid_phase)
    reservoir%s_g(k) = wippflo_auxvar%sat(option%gas_phase)
    reservoir%mobility_l(k) = &
      wippflo_auxvar%mobility(option%liquid_phase)
    reservoir%mobility_g(k) = wippflo_auxvar%mobility(option%gas_phase)
    reservoir%kr_l(k) = wippflo_auxvar%kr(option%liquid_phase)
    reservoir%kr_g(k) = wippflo_auxvar%kr(option%gas_phase)
    reservoir%den_l(k) = wippflo_auxvar%den_kg(option%liquid_phase)
    reservoir%den_g(k) = wippflo_auxvar%den_kg(option%gas_phase)
    reservoir%visc_l(k) = wippflo_auxvar%mu(option%liquid_phase)
    reservoir%visc_g(k) = wippflo_auxvar%mu(option%gas_phase)
    reservoir%e_por(k) = wippflo_auxvar%effective_porosity

    reservoir%kx(k) = material_auxvar%permeability(1)
    reservoir%ky(k) = material_auxvar%permeability(2)
    reservoir%kz(k) = material_auxvar%permeability(3)
    reservoir%volume(k) = material_auxvar%volume

    if (res_grid%itype == STRUCTURED_GRID) then
      reservoir%dx(k) = res_grid%structured_grid%dx(ghosted_id)
      reservoir%dy(k) = res_grid%structured_grid%dy(ghosted_id)
      reservoir%dz(k) = res_grid%structured_grid%dz(ghosted_id)
    else
      reservoir%dz(k) = pm_well%well_grid%res_dz(k)
      reservoir%dx(k) = sqrt(material_auxvar%volume/ &
                                  reservoir%dz(k))
      reservoir%dy(k) = reservoir%dx(k)
    endif

    if (pm_well%transport) then
      reservoir%aqueous_conc(:,k) = nwt_auxvar%aqueous_eq_conc(:)
      reservoir%aqueous_mass(:,k) = &
            reservoir%aqueous_conc(:,k) * reservoir%e_por(k) * &
            reservoir%volume(k) * reservoir%s_l(k)
    endif

    ! unused
    reservoir%xmass_liq(k,:) = 0.d0
    reservoir%xmass_gas(k,:) = 0.d0

  enddo

  if (wippflo_update_index < 0 .or. initialize_well_flow) then
  if (option%myrank == pm_well%well_grid%h_rank_id(1)) then
      root_rank = pm_well%well_comm%rank
  endif 
  call MPI_Bcast(root_rank,1,MPI_INTEGER,pm_well%well_grid%h_rank_id(1), &
                 option%mycomm,ierr);CHKERRQ(ierr)

  if (well_comm%commsize > 1) then
    do k = 1,pm_well%well_grid%nsegments
      TAG = k
      if (pm_well%well_grid%h_rank_id(k) /= pm_well%well_grid%h_rank_id(1)) then
        if (option%myrank == pm_well%well_grid%h_rank_id(k)) then
          peer = pm_well%well_grid%h_rank_id(1)
        endif
        if (option%myrank == pm_well%well_grid%h_rank_id(1)) then
          peer = pm_well%well_grid%h_rank_id(k)
        endif
        if ((option%myrank == pm_well%well_grid%h_rank_id(k)) .or. &
            (option%myrank == pm_well%well_grid%h_rank_id(1))) then
          call MPI_Sendrecv_replace(reservoir%p_l(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%p_g(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%temp(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%s_l(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%s_g(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%mobility_l(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%mobility_g(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%kr_l(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%kr_g(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%den_l(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%den_g(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%visc_l(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%visc_g(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%e_por(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%xmass_liq(k,:),&
                 option%nflowspec,MPI_DOUBLE_PRECISION,peer,TAG,peer, &
                 TAG,option%mycomm,MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%xmass_gas(k,:),&
                 option%nflowspec,MPI_DOUBLE_PRECISION,peer,TAG,peer, &
                 TAG,option%mycomm,MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%kx(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%ky(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%kz(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%volume(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%dx(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%dy(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          call MPI_Sendrecv_replace(reservoir%dz(k),1, &
                 MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                 MPI_STATUS_IGNORE,ierr)
          if (pm_well%transport) then
            call MPI_Sendrecv_replace(reservoir%aqueous_conc(:,k), &
                   pm_well%nspecies,MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG, &
                   option%mycomm,MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%aqueous_mass(:,k), &
                   pm_well%nspecies,MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG, &
                   option%mycomm,MPI_STATUS_IGNORE,ierr)
          endif
        endif
      endif
    enddo
    
    if (pm_well%well_comm%comm /= MPI_COMM_NULL) then
      call MPI_Bcast(reservoir%p_l,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%p_g,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%temp,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%s_l,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%s_g,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%mobility_l,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%mobility_g,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%kr_l,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%kr_g,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%den_l,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%den_g,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%visc_l,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%visc_g,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%e_por,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%kx,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%ky,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%kz,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%volume,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%dx,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%dy,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      call MPI_Bcast(reservoir%dz,pm_well%well_grid%nsegments, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      do k = 1,pm_well%well_grid%nsegments
        if (pm_well%transport) then
            call MPI_Bcast(reservoir%aqueous_conc(:,k),pm_well%nspecies, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
          call MPI_Bcast(reservoir%aqueous_mass(:,k),pm_well%nspecies, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
        endif
        call MPI_Bcast(reservoir%xmass_liq(k,:),option%nflowspec, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%xmass_gas(k,:),option%nflowspec, &
                     MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                     ierr);CHKERRQ(ierr)
      enddo
    endif
  endif
  endif

end subroutine PMWellUpdateReservoirWIPP

! ************************************************************************** !

subroutine PMWellUpdateReservoirSCO2(pm_well,update_index,segment_index)
  !
  ! Updates the SCO2 mode reservoir properties for the well process model.
  !
  ! Author: Michael Nole
  ! Date: 03/06/2024

  use SCO2_Aux_module
  use Material_Aux_module
  use Grid_module

  implicit none

  class(pm_well_type) :: pm_well
  PetscInt :: update_index
  PetscInt :: segment_index
  
  type(well_reservoir_type), pointer :: reservoir
  type(sco2_auxvar_type), pointer :: sco2_auxvar
  type(material_auxvar_type), pointer :: material_auxvar
  type(grid_type), pointer :: res_grid
  type(option_type), pointer :: option
  type(well_comm_type), pointer :: well_comm
  PetscInt :: TAG, peer, root_rank
  PetscInt :: k, indx
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  option => pm_well%option
  well_comm => pm_well%well_comm 
  reservoir => pm_well%well%reservoir

  res_grid => pm_well%realization%patch%grid

  if (update_index < ZERO_INTEGER) then
    indx = ZERO_INTEGER
  else
    indx = update_index
  endif

  do k = 1,pm_well%well_grid%nsegments
    if (pm_well%well_grid%h_rank_id(k) /= option%myrank) cycle

    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    if (Initialized(segment_index)) then
      ! Use reservoir perturbation: only perturb one reservoir variable
      ! at a time.
      if (k == segment_index) then
        sco2_auxvar => &
          pm_well%realization%patch%aux%sco2%auxvars(indx,ghosted_id)
      else
        sco2_auxvar => &
          pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)
      endif
    else
      sco2_auxvar => &
        pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)
    endif

    if (k == 1 .and. indx == 0) pm_well%well%bh_p = sco2_auxvar%well%bh_p
      
    material_auxvar => &
      pm_well%realization%patch%aux%material%auxvars(ghosted_id)

    reservoir%p_l(k) = sco2_auxvar%pres(option%liquid_phase)
    reservoir%p_g(k) = sco2_auxvar%pres(option%gas_phase)
    reservoir%s_l(k) = sco2_auxvar%sat(option%liquid_phase)
    reservoir%s_g(k) = sco2_auxvar%sat(option%gas_phase)
    reservoir%temp(k) = sco2_auxvar%temp
    reservoir%mobility_l(k) = &
      sco2_auxvar%mobility(option%liquid_phase)
    reservoir%mobility_g(k) = sco2_auxvar%mobility(option%gas_phase)
    reservoir%kr_l(k) = sco2_auxvar%kr(option%liquid_phase)
    reservoir%kr_g(k) = sco2_auxvar%kr(option%gas_phase)
    reservoir%den_l(k) = sco2_auxvar%den_kg(option%liquid_phase)
    reservoir%den_g(k) = sco2_auxvar%den_kg(option%gas_phase)
    reservoir%visc_l(k) = sco2_auxvar%visc(option%liquid_phase)
    reservoir%visc_g(k) = sco2_auxvar%visc(option%gas_phase)
    reservoir%xmass_liq(k,:) = sco2_auxvar%xmass(:, &
                                                    option%liquid_phase)
    reservoir%xmass_gas(k,:) = sco2_auxvar%xmass(:, &
                                                    option%gas_phase)
    reservoir%e_por(k) = sco2_auxvar%effective_porosity

    reservoir%kx(k) = material_auxvar%permeability(1)
    reservoir%ky(k) = material_auxvar%permeability(2)
    reservoir%kz(k) = material_auxvar%permeability(3)
    reservoir%volume(k) = material_auxvar%volume

    if (res_grid%itype == STRUCTURED_GRID) then
      reservoir%dx(k) = res_grid%structured_grid%dx(ghosted_id)
      reservoir%dy(k) = res_grid%structured_grid%dy(ghosted_id)
      reservoir%dz(k) = res_grid%structured_grid%dz(ghosted_id)
    else
      reservoir%dz(k) = pm_well%well_grid%res_dz(k)
      reservoir%dx(k) = sqrt(material_auxvar%volume/ &
                                  reservoir%dz(k))
      reservoir%dy(k) = reservoir%dx(k)
    endif
  enddo

  if (update_index < 0 .or. initialize_well_flow) then
    if (option%myrank == pm_well%well_grid%h_rank_id(1)) then
        root_rank = pm_well%well_comm%rank
    endif 
    call MPI_Bcast(root_rank,1,MPI_INTEGER,pm_well%well_grid%h_rank_id(1), &
                   option%mycomm,ierr);CHKERRQ(ierr)

    if (well_comm%commsize > 1) then
      do k = 1,pm_well%well_grid%nsegments
        TAG = k
        if (pm_well%well_grid%h_rank_id(k) /= &
            pm_well%well_grid%h_rank_id(1)) then
          if (option%myrank == pm_well%well_grid%h_rank_id(k)) then
            peer = pm_well%well_grid%h_rank_id(1)
          endif
          if (option%myrank == pm_well%well_grid%h_rank_id(1)) then
            peer = pm_well%well_grid%h_rank_id(k)
          endif
          if ((option%myrank == pm_well%well_grid%h_rank_id(k)) .or. &
              (option%myrank == pm_well%well_grid%h_rank_id(1))) then
            if (k == 1) then
              call MPI_Sendrecv_replace(pm_well%well%bh_p,1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            endif
            call MPI_Sendrecv_replace(reservoir%p_l(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%p_g(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%temp(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%s_l(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%s_g(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%mobility_l(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%mobility_g(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%kr_l(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%kr_g(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%den_l(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%den_g(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%visc_l(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%visc_g(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%xmass_liq(k,:), &
                   option%nflowspec,MPI_DOUBLE_PRECISION,peer,TAG,peer, &
                   TAG,option%mycomm,MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%xmass_gas(k,:), &
                   option%nflowspec,MPI_DOUBLE_PRECISION,peer,TAG,peer, &
                   TAG,option%mycomm,MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%e_por(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%kx(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%ky(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%kz(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%volume(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%dx(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%dy(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            call MPI_Sendrecv_replace(reservoir%dz(k),1, &
                   MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG,option%mycomm, &
                   MPI_STATUS_IGNORE,ierr)
            if (pm_well%transport) then
              call MPI_Sendrecv_replace(reservoir%aqueous_conc(:,k), &
                     pm_well%nspecies,MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG, &
                     option%mycomm,MPI_STATUS_IGNORE,ierr)
              call MPI_Sendrecv_replace(reservoir%aqueous_mass(:,k), &
                     pm_well%nspecies,MPI_DOUBLE_PRECISION,peer,TAG,peer,TAG, &
                     option%mycomm,MPI_STATUS_IGNORE,ierr)
            endif
          endif
        endif
      enddo
    
      if (pm_well%well_comm%comm /= MPI_COMM_NULL) then
        call MPI_Bcast(reservoir%p_l,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%p_g,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%temp,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%s_l,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%s_g,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%mobility_l, &
                       pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%mobility_g, &
                       pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%kr_l,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%kr_g,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%den_l,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%den_g,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%visc_l,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%visc_g,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%e_por,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%kx,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%ky,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%kz,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%volume,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%dx,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%dy,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        call MPI_Bcast(reservoir%dz,pm_well%well_grid%nsegments, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        do k = 1,pm_well%well_grid%nsegments
          if (pm_well%transport) then
            call MPI_Bcast(reservoir%aqueous_conc(:,k), &
                       pm_well%nspecies,MPI_DOUBLE_PRECISION,root_rank, &
                       pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
            call MPI_Bcast(reservoir%aqueous_mass(:,k), &
                       pm_well%nspecies,MPI_DOUBLE_PRECISION,root_rank, &
                       pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
          endif
          call MPI_Bcast(reservoir%xmass_liq(k,:),option%nflowspec, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
          call MPI_Bcast(reservoir%xmass_gas(k,:),option%nflowspec, &
                       MPI_DOUBLE_PRECISION,root_rank,pm_well%well_comm%comm, &
                       ierr);CHKERRQ(ierr)
        enddo
      endif
    endif
  endif

end subroutine PMWellUpdateReservoirSCO2

! ************************************************************************** !

subroutine PMWellUpdateTimestep(this,update_dt, &
                                dt,dt_min,dt_max,iacceleration, &
                                num_newton_iterations,tfac, &
                                time_step_max_growth_factor)
  !
  ! Updates the time step for the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  ! placeholder

end subroutine PMWellUpdateTimestep

! ************************************************************************** !

subroutine PMWellFinalizeTimestep(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_type) :: this

  PetscReal :: curr_time

  curr_time = this%option%time - this%option%flow_dt

  if (Initialized(this%intrusion_time_start) .and. &
      (curr_time < this%intrusion_time_start) .and. &
      .not. this%well_on) return

  call PMWellUpdateReservoirSrcSinkFlow(this)
  if (this%transport) then
    call PMWellUpdateReservoirSrcSinkTran(this)
  endif

  call PMWellUpdateMass(this)

  call PMWellMassBalance(this)

  if (this%print_well) then
    call PMWellOutput(this)
  endif

end subroutine PMWellFinalizeTimestep

! ************************************************************************** !

subroutine PMWellUpdateReservoirSrcSinkFlow(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/21/2021
  !

  use Coupler_module
  use NW_Transport_Aux_module
  use Transport_Constraint_NWT_module
  use Option_module
  use SCO2_Aux_module, only: sco2_well_coupling, &
                             SCO2_FULLY_IMPLICIT_WELL

  implicit none

  class(pm_well_type) :: pm_well

  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: srcsink_name
  type(coupler_type), pointer :: source_sink
  PetscInt :: k, ghosted_id
  PetscReal :: well_delta_liq, well_delta_gas
  PetscReal :: density_avg

  option => pm_well%option

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  do k = 1,pm_well%well_grid%nsegments
    if (pm_well%well_grid%h_rank_id(k) /= option%myrank) cycle

    write(string,'(I0.6)') k
    srcsink_name = 'well_segment_' // trim(string)

    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    ! [kg-liq/m3]
    density_avg = 0.5d0 * (pm_well%well%liq%den(k) + pm_well%well%reservoir%den_l(k))

    source_sink => pm_well%realization%patch%source_sink_list%first
    do
      if (.not.associated(source_sink)) exit

      if (trim(srcsink_name) == trim(source_sink%name)) then

        select case(option%iflowmode)
          case(WF_MODE)
            if (wippflo_well_quasi_imp_coupled) then
              source_sink%flow_condition%general%rate%dataset%rarray(1) = &
                0.d0 ! [kmol/s]
              source_sink%flow_condition%general%rate%dataset%rarray(2) = &
                0.d0 ! [kmol/s]
            else
              source_sink%flow_condition%general%rate%dataset%rarray(1) = &
                -1.d0 * pm_well%well%liq%Q(k) ! [kmol/s]
              source_sink%flow_condition%general%rate%dataset%rarray(2) = &
                -1.d0 * pm_well%well%gas%Q(k) ! [kmol/s]
            endif

            source_sink%flow_condition%general%liquid_pressure%aux_real(1) = &
                                                               pm_well%well%pl(k)
            source_sink%flow_condition%general%gas_pressure%aux_real(1) = &
                                                               pm_well%well%pg(k)
            well_delta_liq = pm_well%well%pl(k) - pm_well%well%reservoir%p_l(k)
            well_delta_gas = pm_well%well%pg(k) - pm_well%well%reservoir%p_g(k)
            source_sink%flow_condition%general%liquid_pressure%aux_real(2) = &
                                                               well_delta_liq
            source_sink%flow_condition%general%gas_pressure%aux_real(2) = &
                                                               well_delta_gas
            source_sink%flow_condition%well%aux_real(1) = density_avg ! kg/m3

            pm_well%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
                 well%pl = pm_well%well%pl(k)
            pm_well%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
                 well%pg = pm_well%well%pg(k)
            pm_well%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
                 well%sl = pm_well%well%liq%s(k)
            pm_well%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
                 well%sg = pm_well%well%gas%s(k)
            pm_well%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
                 well%dpl = well_delta_liq
            pm_well%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
                 well%dpg = well_delta_gas
            pm_well%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
                 well%Ql = pm_well%well%liq%Q(k)
            pm_well%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
                 well%Qg = pm_well%well%gas%Q(k)
          case(SCO2_MODE)
            if (sco2_well_coupling == SCO2_FULLY_IMPLICIT_WELL) then
              source_sink%flow_condition%sco2%rate%dataset%rarray(1) = &
                0.d0 ! [kmol/s]
              source_sink%flow_condition%sco2%rate%dataset%rarray(2) = &
                0.d0 ! [kmol/s]
            else
              source_sink%flow_condition%sco2%rate%dataset%rarray(1) = &
                -1.d0 * pm_well%well%liq%Q(k) ! [kmol/s]
              source_sink%flow_condition%sco2%rate%dataset%rarray(2) = &
                -1.d0 * pm_well%well%gas%Q(k) ! [kmol/s]
            endif

            source_sink%flow_condition%sco2%liquid_pressure%aux_real(1) = &
                                                               pm_well%well%pl(k)
            source_sink%flow_condition%sco2%gas_pressure%aux_real(1) = &
                                                               pm_well%well%pg(k)
            well_delta_liq = pm_well%well%pl(k) - pm_well%well%reservoir%p_l(k)
            well_delta_gas = pm_well%well%pg(k) - pm_well%well%reservoir%p_g(k)
            source_sink%flow_condition%sco2%liquid_pressure%aux_real(2) = &
                                                               well_delta_liq
            source_sink%flow_condition%sco2%gas_pressure%aux_real(2) = &
                                                               well_delta_gas
            source_sink%flow_condition%well%aux_real(1) = density_avg ! kg/m3

            pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)%&
                 well%pl = pm_well%well%pl(k)
            pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)%&
                 well%pg = pm_well%well%pg(k)
            pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)%&
                 well%sl = pm_well%well%liq%s(k)
            pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)%&
                 well%sg = pm_well%well%gas%s(k)
            pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)%&
                 well%dpl = well_delta_liq
            pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)%&
                 well%dpg = well_delta_gas
            pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)%&
                 well%Ql = pm_well%well%liq%Q(k)
            pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)%&
                 well%Qg = pm_well%well%gas%Q(k)
            pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)%&
                 well%bh_p = pm_well%well%bh_p
        end select

        exit
      endif

      source_sink => source_sink%next
    enddo

  enddo

end subroutine PMWellUpdateReservoirSrcSinkFlow

! ************************************************************************** !

subroutine PMWellUpdateReservoirSrcSinkTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/21/2021
  !

  use Coupler_module
  use NW_Transport_Aux_module
  use Transport_Constraint_NWT_module
  use Option_module

  implicit none

  class(pm_well_type) :: pm_well

  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: srcsink_name
  type(nw_transport_auxvar_type), pointer :: nwt_auxvar
  type(coupler_type), pointer :: source_sink
  PetscInt :: k, ghosted_id
  PetscReal :: density_avg
  PetscErrorCode :: ierr

  option => pm_well%option

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return
  ierr = 0

  do k = 1,pm_well%well_grid%nsegments
    if (pm_well%well_grid%h_rank_id(k) /= pm_well%option%myrank) cycle

    write(string,'(I0.6)') k
    srcsink_name = 'well_segment_' // trim(string)

    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    ! [kg-liq/m3]
    density_avg = 0.5d0 * (pm_well%well%liq%den(k) + pm_well%well%reservoir%den_l(k))

    source_sink => pm_well%realization%patch%source_sink_list%first
    do
      if (.not.associated(source_sink)) exit

      if (trim(srcsink_name) == trim(source_sink%name)) then
        select case(option%iflowmode)
          case (WF_MODE)
            source_sink%flow_condition%general%rate%dataset%rarray(1) = &
              -1.d0 * pm_well%well%liq%Q(k) ! [kmol/s]
            source_sink%flow_condition%general%rate%dataset%rarray(2) = &
              -1.d0 * pm_well%well%gas%Q(k) ! [kmol/s]
          case (SCO2_MODE)
            source_sink%flow_condition%sco2%rate%dataset%rarray(1) = &
              -1.d0 * pm_well%well%liq%Q(k) ! [kg/s]
            source_sink%flow_condition%sco2%rate%dataset%rarray(2) = &
              -1.d0 * pm_well%well%gas%Q(k) ! [kg/s]
        end select
        
        source_sink%flow_condition%well%aux_real(1) = density_avg ! kg/m3

        ! access nwt_auxvar from the tran_condition
        nwt_auxvar => &
          TranConstraintNWTGetAuxVar(source_sink%tran_condition% &
                                     cur_constraint_coupler)

        ! modify nwt_auxvar from the tran_condition
        nwt_auxvar%aqueous_eq_conc(:) = pm_well%well%aqueous_conc(:,k)

        pm_well%realization%patch%aux%nwt%auxvars(ghosted_id)%&
          well%AQ_conc(1:pm_well%nspecies) = pm_well%well%aqueous_conc(:,k)
        pm_well%realization%patch%aux%nwt%auxvars(ghosted_id)%&
          well%AQ_mass(1:pm_well%nspecies) = pm_well%well%aqueous_mass(:,k)

        exit
      endif

      source_sink => source_sink%next
    enddo
  enddo
  call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)

end subroutine PMWellUpdateReservoirSrcSinkTran

! ************************************************************************** !

subroutine PMWellUpdateFlowRates(this,well_pert,res_pert,segment_index,ierr)
  !
  ! This subroutine performs the well rate computation when called from the
  ! fully- or quasi-coupled source/sink update.
  !
  ! Author: Michael Nole
  ! Date: 01/16/2023
  !

  Use Option_module

  implicit none

  class(pm_well_type) :: this
  PetscErrorCode :: ierr
  PetscInt :: well_pert
  PetscInt :: res_pert
  PetscInt :: segment_index

  type(option_type), pointer :: option
  PetscReal :: time
  PetscInt :: k


  option => this%option
  time = this%realization%option%time
!  if (this%flow_soln%n_steps < 1) return

  ! Need to limit well model timestepping 
  this%min_dt_flow = option%flow_dt * min_flow_dt_scale

  if (.not. this%well_on .and. Initialized(this%intrusion_time_start) .and. &
      time < this%intrusion_time_start) then
      this%srcsink_water(well_pert,:) = 0.d0
      this%srcsink_gas(well_pert,:) = 0.d0
    return
  elseif (.not. this%well_on) then
    this%well_on = PETSC_TRUE
  endif

  this%print_output = PETSC_FALSE
  select case (option%iflowmode)
    case(WF_MODE)
      call PMWellUpdateReservoirWIPP(this,res_pert)
    case(SCO2_MODE)
      call PMWellUpdateReservoirSCO2(this,res_pert,segment_index)
  end select

  if (initialize_well_flow) then
    call PMWellInitializeWellFlow(this)
    do k = 1,option%nflowdof
      call PMWellCopyWell(this%well,this%well_pert(k),this%transport)
    enddo
  elseif (option%iflowmode == WF_MODE) then
    this%well%pl = this%flow_soln%soln_save%pl
    this%well%gas%s = this%flow_soln%soln_save%sg
    this%well%bh_p = this%flow_soln%soln_save%bh_p
  endif
  if (option%iflowmode == WF_MODE) then
    ! Quasi-implicit
    do k = 1,option%nflowdof
      call PMWellCopyWell(this%well,this%well_pert(k),this%transport)
    enddo
  
    call PMWellUpdatePropertiesWIPPFlow(this,this%well,&
                        this%realization%patch%characteristic_curves_array,&
                        this%realization%option)
  endif
  this%dt_flow = this%realization%option%flow_dt
  call PMWellSolveFlow(this,well_pert,ierr)
  this%print_output = PETSC_TRUE  

  this%srcsink_water(well_pert,:) = -1.d0 * &
                                          this%well%liq%Q(:)! [kmol/s]
  this%srcsink_gas(well_pert,:)   = -1.d0 * &
                                          this%well%gas%Q(:)! [kmol/s]

end subroutine PMWellUpdateFlowRates

! ************************************************************************** !

subroutine PMWellModifyFlowResidual(this,residual,ss_flow_vol_flux)
  !
  ! This subroutine computes the well contribution to the reservoir residual when 
  ! called from the fully- or quasi-coupled source/sink update in WIPP FLOW mode.
  !
  ! Author: Michael Nole
  ! Date: 01/16/2023
  !

  implicit none

  class(pm_well_type) :: this
  PetscReal, pointer :: residual(:)
  PetscReal :: ss_flow_vol_flux(this%option%nphase)

  PetscReal :: Res(this%flow_soln%ndof)
  PetscReal :: Q, sum_q
  PetscInt :: air_comp_id, wat_comp_id
  PetscInt :: local_id, local_start, local_end
  PetscInt :: ghosted_id
  PetscInt :: j,k

  air_comp_id = this%option%air_id
  wat_comp_id = this%option%water_id

  Res(:) = 0.d0

  if (this%well_comm%comm == MPI_COMM_NULL) return

  select case (this%option%iflowmode)
    case(WF_MODE)
      do k = 1,this%well_grid%nsegments
        if (this%well_grid%h_rank_id(k) /= this%option%myrank) cycle

        ghosted_id = this%well_grid%h_ghosted_id(k)
        local_id = this%realization%patch%grid%nG2L(ghosted_id)
        local_end = local_id * this%flow_soln%ndof
        local_start = local_end - this%flow_soln%ndof + 1

        ! kmol/sec
        Res(wat_comp_id) = this%srcsink_water(UNPERT,k)
        Res(air_comp_id) = this%srcsink_gas(UNPERT,k)

        call WIPPFloConvertUnitsToBRAGFlo(Res,this%realization%patch% &
                                          aux%Material% &
                                          auxvars(ghosted_id),this%option)
        residual(local_start:local_end) = residual(local_start:local_end) - &
                                          Res(:)
    
      enddo
    case(SCO2_MODE)
      select case (this%well%well_model_type)

        case(WELL_MODEL_HYDROSTATIC)
          ! Add residual for bottom hole cell
          select case (this%flow_coupling)
            case(FULLY_IMPLICIT_WELL)
              sum_q = 0.d0
              do k = 1,this%well_grid%nsegments
                if (this%well_grid%h_rank_id(k) /= this%option%myrank) &
                   cycle

                ghosted_id = this%well_grid%h_ghosted_id(k)
                local_id = this%realization%patch%grid%nG2L(ghosted_id)
                local_end = local_id * this%option%nflowdof
                local_start = local_end - this%option%nflowdof + 1
                if (k == 1) then
                  ! An extra residual is required for the bottom well cell.
                  ! Residual = Q - sum(q)
                  if (dabs(this%well%th_qg) > 0.d0) then
                    sum_q = sum_q + sum(this%well%gas%q)
                    Q = -1.d0 * (this%well%th_qg)
                  else
                    sum_q = sum(this%well%liq%q)
                    Q = -1.d0 * (this%well%th_ql)
                  endif
                  residual(local_end) = Q - sum_q
                endif
                do j = 0, this%option%nflowdof-2
                  ! Compontent j+1 residual at well segment k
                  residual(local_start + j) = residual(local_start + j) + &
                          (this%well%liq%q(k)* &
                          this%well%liq%xmass(k,j+1) + &
                          this%well%gas%q(k)* &
                          this%well%gas%xmass(k,j+1))
                enddo
              enddo
          end select

        case(WELL_MODEL_FULL_MOMENTUM)
      
      case default

      end select

  end select

end subroutine PMWellModifyFlowResidual

! ************************************************************************** !

subroutine PMWellModifyFlowJacobian(this,Jac,ierr)
  !
  ! This subroutine computes the well contribution to the reservoir Jacobian when 
  ! called from the fully- or quasi-coupled source/sink update.
  !
  ! Author: Michael Nole
  ! Date: 01/16/2023
  !

  use Option_module
  use SCO2_Aux_module, only : SCO2_WELL_DOF

  implicit none

  class(pm_well_type) :: this
  Mat :: Jac
  PetscErrorCode :: ierr

  type(well_type), pointer :: well
  type(well_type), pointer :: well_pert
  type(option_type), pointer :: option
  PetscReal :: pert_pl, pert_sg
  PetscReal :: res, res_pert_pl, res_pert_sg
  PetscReal :: pert, res_pert
  PetscReal :: Q, sum_q
  PetscInt :: ghosted_id
  PetscReal, allocatable :: residual(:,:)
  PetscReal, allocatable :: J_block(:,:)
  PetscInt :: air_comp_id, wat_comp_id
  PetscInt :: j,k,idof,irow
  PetscInt :: local_row_index, local_col_index
  PetscReal :: J_well

  option => this%option

  air_comp_id = this%option%air_id
  wat_comp_id = this%option%water_id

  if (this%well_comm%comm == MPI_COMM_NULL) return

  select case (this%option%iflowmode)
    case(WF_MODE)

      allocate(J_block(this%flow_soln%ndof,this%flow_soln%ndof))
      J_block = 0.d0

      do k = 1,this%well_grid%nsegments
        if (this%well_grid%h_rank_id(k) /= this%option%myrank) cycle

        J_block = 0.d0

        ghosted_id = this%well_grid%h_ghosted_id(k)
        pert_pl = this%realization%patch%aux%WIPPFlo% &
                       auxvars(WIPPFLO_LIQUID_PRESSURE_DOF,ghosted_id)%pert
        pert_sg = this%realization%patch%aux%WIPPFlo% &
                       auxvars(WIPPFLO_GAS_SATURATION_DOF,ghosted_id)%pert

        ! Liquid portion
        res = this%srcsink_water(UNPERT,k)
        res_pert_pl = this%srcsink_water(PERT_WRT_PL,k)
        res_pert_sg = this%srcsink_water(PERT_WRT_SG,k)
        J_block(wat_comp_id,WIPPFLO_LIQUID_PRESSURE_DOF) = &
               (res_pert_pl - res)/pert_pl
        J_block(wat_comp_id,WIPPFLO_GAS_SATURATION_DOF) = &
               (res_pert_sg - res)/pert_sg

        ! Gas portion
        res = this%srcsink_gas(UNPERT,k)
        res_pert_pl = this%srcsink_gas(PERT_WRT_PL,k)
        res_pert_sg = this%srcsink_gas(PERT_WRT_SG,k)
        J_block(air_comp_id,WIPPFLO_LIQUID_PRESSURE_DOF) = &
               (res_pert_pl - res)/pert_pl
        J_block(air_comp_id,WIPPFLO_GAS_SATURATION_DOF) = &
               (res_pert_sg - res)/pert_sg

        J_block = -1.d0*J_block
        call WIPPFloConvertUnitsToBRAGFlo(J_block,this%realization%patch%aux% &
                                          Material%auxvars(ghosted_id),this%option)

        call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1,J_block, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
      enddo
    case(SCO2_MODE)

      allocate(J_block(option%nflowdof,option%nflowdof))
      allocate(residual(this%well_grid%nsegments,option%nflowdof))

      J_block = 0.d0
      residual = 0.d0

      if (this%flow_coupling == FULLY_IMPLICIT_WELL) then
        ! Calculate Jacobian entries wrt reservoir perturbation:
        ! BHP residual and reservoir source/sink residuals.

        ! Unperturbed residual
        well => this%well
        do k = 1,this%well_grid%nsegments
          do irow = 1, option%nflowdof-1
            residual(k,irow) = well%liq%q(k)* &
                               well%liq%xmass(k,irow) + &
                               well%gas%q(k)* &
                               well%gas%xmass(k,irow)
          enddo
          if (k==1) then
            sum_q = 0.d0
            if (dabs(well%th_qg) > 0.d0) then
              sum_q = sum(well%gas%q)
              Q = -1.d0 * (well%th_qg)
            else
              sum_q = sum(well%liq%q)
              Q = -1.d0 * (well%th_ql)
            endif
            residual(k,option%nflowdof) = Q - sum_q
          endif
        enddo

        ! Perturbations
        do k = 1,this%well_grid%nsegments
          if (this%well_grid%h_rank_id(k) /= option%myrank) cycle

          J_block = 0.d0
          ghosted_id = this%well_grid%h_ghosted_id(k)
          do idof = 1,option%nflowdof-1
            well_pert => this%well_pert(idof)
            pert = this%realization%patch%aux%SCO2% &
                         auxvars(idof,ghosted_id)%pert

            ! Compute dRwell / dXres
            if (this%well_grid%h_rank_id(1) == option%myrank) then
              res_pert = 0.d0
              if (dabs(well_pert%th_qg) > 0.d0) then
                res_pert = well_pert%gas%q(k) - well%gas%q(k)
              else
                res_pert = well_pert%liq%q(k) - well%liq%q(k)
              endif

              local_row_index = this%well_grid%h_ghosted_id(1)*option%nflowdof-1
              local_col_index = (ghosted_id-1)*option%nflowdof + idof - 1
              J_well = -(res_pert)/pert
              call MatSetValuesLocal(Jac,1,local_row_index,1,local_col_index,J_well, &
                                     ADD_VALUES,ierr);CHKERRQ(ierr)

            endif

            ! Compute dRres / dXres
            do irow = 1,option%nflowdof-1
              res_pert = well_pert%liq%q(k)* &
                         well_pert%liq%xmass(k,irow) + &
                         well_pert%gas%q(k)* &
                         well_pert%gas%xmass(k,irow)
              J_block(irow,idof) = (res_pert - residual(k,irow))/pert
            enddo
          enddo
          call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1,&
                                       J_block,ADD_VALUES,ierr);CHKERRQ(ierr)

          ! Make sure perturbation generates flow.
          ! pres_bump = 0.d0
          ! i = 0
          ! do
          !   call this%UpdateFlowRates(ONE_INTEGER,ZERO_INTEGER,-999,ierr)
          !   if (dabs(this%well%th_ql) > 0.d0 .or. &
          !       dabs(this%well%th_qg) > 0.d0) then
          !     if (any(dabs(this%well_pert(ONE_INTEGER)%gas%Q) > 0.d0)) exit
          !     pres_bump = this%well_pert(ONE_INTEGER)%bh_p - &
          !            this%well%bh_p
          !     this%well_pert(ONE_INTEGER)%bh_p = &
          !                            this%well_pert(ONE_INTEGER)%bh_p - pres_bump
          !     pres_bump = pres_bump * 1.25d0
          !     this%well_pert(ONE_INTEGER)%bh_p = &
          !                            this%well_pert(ONE_INTEGER)%bh_p + pres_bump
          !     i = i + 1
          !     if (i > 100) then
          !       option%io_buffer = "Maximum number of iterations exceeded to recover &
          !                           &vanishing well Jacobian. "
          !       call PrintMsg(this%option)
          !     endif
          !   else
          !     exit
          !   endif
          ! enddo

          !Perturbed rates wrt well perturbation.
          well_pert => this%well_pert(SCO2_WELL_DOF)

          if (this%well_grid%h_rank_id(k) /= option%myrank) cycle
  
          J_block = 0.d0
  
          ghosted_id = this%well_grid%h_ghosted_id(k)
          pert = well_pert%bh_p - well%bh_p

          ! Compute dRwell / dPwell
          if (k == 1) then
            sum_q = 0.d0
            if (dabs(well%th_qg) > 0.d0) then
              sum_q = sum(well%gas%q)
              Q = -1.d0 * (well%th_qg)
            else
              sum_q = sum(well%liq%q)
              Q = -1.d0 * (well%th_ql)
            endif
            res = Q - sum_q

            sum_q = 0.d0
            if (dabs(well_pert%th_qg) > 0.d0) then
              sum_q = sum(well_pert%gas%q)
              Q = -1.d0 * (well_pert%th_qg)
            else
              sum_q = sum(well_pert%liq%q)
              Q = -1.d0 * (well_pert%th_ql)
            endif
            res_pert = Q - sum_q

            J_block(option%nflowdof,option%nflowdof) = &
                    (res_pert - res)/pert
          endif
  
          ! Compute dRres / dPwell
          do j = 0,option%nflowdof-2
            res = this%well%liq%q(k)*this%well%liq%xmass(k,j+1) + &
                  this%well%gas%q(k)*this%well%gas%xmass(k,j+1)
            res_pert = well_pert%liq%q(k)* &
                       well_pert%liq%xmass(k,j+1) + &
                       well_pert%gas%q(k)* &
                       well_pert%gas%xmass(k,j+1)
            ! Just change 1 column
            local_row_index = (ghosted_id-1)*option%nflowdof + j
            local_col_index = this%well_grid%h_ghosted_id(1)*option%nflowdof-1
            J_well = (res_pert - res)/pert
            call MatSetValuesLocal(Jac,1,local_row_index,1,local_col_index,J_well, &
                                   ADD_VALUES,ierr);CHKERRQ(ierr)
          enddo
  
          call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1, &
                                        J_block,ADD_VALUES,ierr);CHKERRQ(ierr)
        enddo
      endif
    end select

end subroutine PMWellModifyFlowJacobian

! ************************************************************************** !

subroutine PMWellResidualFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: i, iup, idn
  PetscReal :: res_accum(pm_well%nphase)
  PetscReal :: res_src_sink(pm_well%nphase)
  PetscReal :: res_flux(pm_well%nphase)
  PetscReal :: res_flux_bc(2*pm_well%nphase)
  PetscReal :: res_temp(pm_well%flow_soln%ndof*pm_well%well_grid%nsegments)

  res_accum = 0.d0
  res_src_sink = 0.d0
  res_flux = 0.d0
  res_flux_bc = 0.d0

  res_temp(:) = pm_well%flow_soln%residual(:)

  select case(pm_well%well%well_model_type)
    !-------------------------------------------------------------------------
    case(WELL_MODEL_WIPP_DARCY)

      call PMWellBCFlux(pm_well,pm_well%well,res_flux_bc,PETSC_TRUE)

      do i = 1,pm_well%well_grid%nsegments
        iup = i
        idn = i + 1

        ! Accumulation Term
        call PMWellAccumulationFlow(pm_well,pm_well%well,i,res_accum)

        res_temp(pm_well%flow_soln%ndof*(i-1)+1) = &
               res_temp(pm_well%flow_soln%ndof*(i-1)+1) + &
               res_accum(ONE_INTEGER)

        ! Source/Sink Term
        call PMWellSrcSink(pm_well,pm_well%well,i,res_src_sink)

        res_temp(pm_well%flow_soln%ndof*(i-1)+1) = &
             res_temp(pm_well%flow_soln%ndof*(i-1)+1) + &
             res_src_sink(ONE_INTEGER)
  
        ! Flux Term
        if (i < pm_well%well_grid%nsegments) then
          call PMWellFlux(pm_well,pm_well%well,pm_well%well,iup,idn,res_flux,PETSC_TRUE)
        endif

        if (i == 1) then
          ! Water mass residual in cell i+1: Subtract flux to i+1 cell
          res_temp(pm_well%flow_soln%ndof*i+1) = &
               res_temp(pm_well%flow_soln%ndof*i+1) &
               - res_flux(1)

          ! Water mass residual in cell i: Add flux in from BC,
          ! add flux to i+1 cell
          res_temp(pm_well%flow_soln%ndof*(i-1)+1) = &
               res_temp(pm_well%flow_soln%ndof*(i-1)+1) &
               - res_flux_bc(1) + res_flux(1)

        elseif (i < pm_well%well_grid%nsegments) then
          ! Water mass residual in cell i: Subtract flux to i+1 cell
          res_temp(pm_well%flow_soln%ndof*(i-1)+1) = &
               res_temp(pm_well%flow_soln%ndof*(i-1)+1) &
               + res_flux(1)
          ! Water mass residual in cell i+1: Add flux to i+1 cell
          res_temp(pm_well%flow_soln%ndof*i+1) = &
               res_temp(pm_well%flow_soln%ndof*i+1) &
               - res_flux(1)
        else
          ! Water mass residual in cell i: Subtract flux to BC
          res_temp(pm_well%flow_soln%ndof*(i-1)+1) = &
               res_temp(pm_well%flow_soln%ndof*(i-1)+1) &
               - res_flux_bc(3)
        endif

        if (pm_well%nphase == 2) then

          res_temp(pm_well%flow_soln%ndof*(i-1)+2) = &
               res_temp(pm_well%flow_soln%ndof*(i-1)+2) + &
               res_accum(TWO_INTEGER)

          res_temp(pm_well%flow_soln%ndof*(i-1)+2) = &
               res_temp(pm_well%flow_soln%ndof*(i-1)+2) + &
               res_src_sink(TWO_INTEGER)

          if (i == 1) then
            ! Air mass residual in cell i+1: Subtract flux to i+1 cell
            res_temp(pm_well%flow_soln%ndof*i+2) = &
                 res_temp(pm_well%flow_soln%ndof*i+2) &
                 - res_flux(2)
            ! Air mass residual in cell i: Subtract flux in from BC,
            ! add flux to i+1 cell
            res_temp(pm_well%flow_soln%ndof*(i-1)+2) = &
                 res_temp(pm_well%flow_soln%ndof*(i-1)+2) &
                 - res_flux_bc(2) + res_flux(2)
          elseif (i < pm_well%well_grid%nsegments) then
            ! Air mass residual in cell i: Subtract flux to i+1 cell
            res_temp(pm_well%flow_soln%ndof*(i-1)+2) = &
                 res_temp(pm_well%flow_soln%ndof*(i-1)+2) &
                 + res_flux(2)
            ! Air mass residual in cell i+1: Add flux to i+1 cell
            res_temp(pm_well%flow_soln%ndof*i+2) = &
                 res_temp(pm_well%flow_soln%ndof*i+2) &
                 - res_flux(2)
          else
            ! Air mass residual in cell i: Subtract flux to BC
            res_temp(pm_well%flow_soln%ndof*(i-1)+2) = &
                 res_temp(pm_well%flow_soln%ndof*(i-1)+2) &
                 - res_flux_bc(4)
          endif
        endif
      enddo
      pm_well%flow_soln%residual(:) = res_temp(:)
    !-------------------------------------------------------------------------
    case(WELL_MODEL_FULL_MOMENTUM)
    !-------------------------------------------------------------------------
  end select

end subroutine PMWellResidualFlow

! ************************************************************************** !

subroutine PMWellQISolveTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 03/13/2023
  !

  implicit none

  class(pm_well_type) :: pm_well
  
  PetscReal :: curr_time 
  PetscErrorCode :: ierr

  ierr = 0
  pm_well%tran_soln%cut_ts_flag = PETSC_FALSE
  curr_time = pm_well%option%time

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  if (initialize_well_tran) then
    call PMWellInitializeWellTran(pm_well)
  endif

  call PMWellUpdatePropertiesTran(pm_well)
  pm_well%dt_tran = pm_well%option%tran_dt 

  call PMWellSolveTran(pm_well,ierr)
  if (pm_well%tran_soln%cut_ts_flag) return

  call PMWellUpdateReservoirSrcSinkTran(pm_well)

end subroutine PMWellQISolveTran

! ************************************************************************** !

subroutine PMWellResidualTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  ! at this time, the tran_soln%residual vector has been zero'd out and has
  ! been loaded with the fixed accumulation divided by current dt

  ! update the auxiliary variables (runs through the equilibrium
  ! diss/precip/sorb routine) - we only need to update the aqueous_conc from
  ! the aqueous_mass value, or vice versa?

  call PMWellResidualTranAccum(pm_well)

  ! calculate the source/sink terms (in/out of well segments)
  call PMWellResidualTranSrcSink(pm_well)

  ! calculate the rxn terms (decay/ingrowth)
  call PMWellResidualTranRxn(pm_well)

  ! calculate the flux terms
  call PMWellResidualTranFlux(pm_well)

end subroutine PMWellResidualTran

! ************************************************************************** !

subroutine PMWellResidualTranAccum(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend
  PetscReal :: Res(pm_well%nspecies)

  ! porosity in [m^3-void/m^3-bulk]
  ! saturation in [m^3-liq/m^3-void]
  ! volume in [m^3-bulk]
  ! aqueous conc in [mol-species/m^3-liq]
  ! residual in [mol-species/sec]

  ! calculate the accumulation term as:
  ! residual = (Res_accum/tran_dt)
  !            - residual(which is already = fixed_accum/dt)

  do isegment = 1,pm_well%well_grid%nsegments

    offset = (isegment-1)*pm_well%nspecies
    istart = offset + 1
    iend = offset + pm_well%nspecies

    do ispecies = 1,pm_well%nspecies
      k = ispecies
      Res(k) = pm_well%well%volume(isegment) * pm_well%well%phi(isegment) * &
               pm_well%well%liq%s(isegment) * &
               pm_well%well%aqueous_conc(ispecies,isegment)
      Res(k) = Res(k) / pm_well%dt_tran
    enddo

    pm_well%tran_soln%residual(istart:iend) = Res(:) - &
                                         pm_well%tran_soln%residual(istart:iend)
  enddo

end subroutine PMWellResidualTranAccum

! ************************************************************************** !

subroutine PMWellResidualTranSrcSink(pm_well)
  !
  ! Calculates the source sink terms (Q in/out of well) for the transport
  ! residual equation.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/06/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend
  PetscReal :: Res(pm_well%nspecies)

  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: resr
  PetscReal :: coef_Qin, coef_Qout ! into well, out of well
  PetscReal :: Qin, Qout
  PetscReal :: den_avg

  ! Q src/sink is in [kmol-liq/sec]
  ! FMWH2O is in [kg-liq/kmol-liq] where liq = water
  ! density is in [kg-liq/m^3-liq] where liq = water
  ! aqueous conc in [mol-species/m^3-liq]
  ! residual in [mol-species/sec]

  ! From the flow solution:
  ! + Q goes into well from reservoir       
  ! - Q goes out of well into reservoir     

  well => pm_well%well
  resr => pm_well%well%reservoir

  do isegment = 1,pm_well%well_grid%nsegments

    den_avg = 0.5d0*(well%liq%den(isegment)+resr%den_l(isegment))
    ! units of coef = [m^3-liq/sec]
    if (well%liq%Q(isegment) < 0.d0) then ! Q out of well
      coef_Qin = 0.d0
      coef_Qout = well%liq%Q(isegment)*FMWH2O/den_avg
    else ! Q into well
    !            [kmol-liq/sec]*[kg-liq/kmol-liq]/[kg-liq/m^3-liq]  
      coef_Qin = well%liq%Q(isegment)*FMWH2O/den_avg
      coef_Qout = 0.d0
    endif

    offset = (isegment-1)*pm_well%nspecies
    istart = offset + 1
    iend = offset + pm_well%nspecies

    do ispecies = 1,pm_well%nspecies
      k = ispecies
      Qin = coef_Qin*resr%aqueous_conc(ispecies,isegment)
      Qout = coef_Qout*well%aqueous_conc(ispecies,isegment)
      Res(k) = Qin + Qout
    enddo

    pm_well%tran_soln%residual(istart:iend) = &
                          pm_well%tran_soln%residual(istart:iend) + Res(:)
  enddo

end subroutine PMWellResidualTranSrcSink

! ************************************************************************** !

subroutine PMWellResidualTranRxn(pm_well)
  !
  ! Calculates the decay/ingrowth of radioactive species for the transport
  ! residual equation.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/06/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend, parent_id
  PetscReal :: Res(pm_well%nspecies)

  ! decay_rate in [1/sec]
  ! aqueous mass in [mol-species]
  ! residual in [mol-species/sec]

  do isegment = 1,pm_well%well_grid%nsegments

    offset = (isegment-1)*pm_well%nspecies
    istart = offset + 1
    iend = offset + pm_well%nspecies

    do ispecies = 1,pm_well%nspecies
      k = ispecies
      ! Add in species decay
      Res(k) = -(pm_well%well%species_decay_rate(k)* &
                          pm_well%well%aqueous_mass(ispecies,isegment))
      ! Add in contribution from parent (if exists)
      parent_id = pm_well%well%species_parent_id(ispecies)
      if (parent_id > 0) then
        Res(k) = Res(k) + (pm_well%well%species_parent_decay_rate(k)* &
                 pm_well%well%aqueous_mass(parent_id,isegment))
      endif
    enddo

    pm_well%tran_soln%residual(istart:iend) = &
                          pm_well%tran_soln%residual(istart:iend) + Res(:)
  enddo

end subroutine PMWellResidualTranRxn

! ************************************************************************** !

subroutine PMWellResidualTranFlux(pm_well)
  !
  ! Calculates the interior and BC flux terms for the transport residual
  ! equation.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/24/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend
  PetscInt :: n_up, n_dn
  PetscReal :: area_up, area_dn
  PetscReal :: q_up, q_dn
  PetscReal :: conc
  PetscReal :: diffusion
  PetscReal :: Res(pm_well%nspecies)
  PetscReal :: Res_up(pm_well%nspecies), Res_dn(pm_well%nspecies)

  ! residual in [mol-species/sec]
  ! area in [m2-bulk]
  ! q_up, d_dn in [m3-liq/m2-bulk-sec]
  ! conc in [mol-species/m3-liq]

  ! NOTE: The up direction is towards well top, and the dn direction is
  !       towards the well bottom.
  !       +q flows up the well
  !       -q flows down the well

  n_dn = +1
  n_up = -1

  diffusion = 0.d0 ! for now, since WIPP has no diffusion

  ! ----------------------------------------INTERIOR-FLUXES------------------

  do isegment = 2,(pm_well%well_grid%nsegments-1)

    Res(:) = 0.d0; Res_up(:) = 0.d0; Res_dn(:) = 0.d0

    area_up = 0.5d0 * (pm_well%well%area(isegment) + pm_well%well%area(isegment+1))
    area_dn = 0.5d0 * (pm_well%well%area(isegment) + pm_well%well%area(isegment-1))

    q_up = pm_well%well%ql(isegment)
    q_dn = pm_well%well%ql(isegment-1)

    offset = (isegment-1)*pm_well%nspecies
    istart = offset + 1
    iend = offset + pm_well%nspecies

    do ispecies = 1,pm_well%nspecies
      k = ispecies

      ! north surface:
      if (q_up < 0.d0) then ! flow is down well
        conc = pm_well%well%aqueous_conc(k,isegment+1)
        Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
      elseif (q_up > 0.d0) then ! flow is up well
        conc = pm_well%well%aqueous_conc(k,isegment)
        Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
      else ! q_up = 0
        Res_up(k) = (n_up*area_up)*(0.d0 - diffusion)
      endif

      ! south surface:
      if (q_dn < 0.d0) then ! flow is down well
        conc = pm_well%well%aqueous_conc(k,isegment)
        Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
      elseif (q_dn > 0.d0) then ! flow up well
        conc = pm_well%well%aqueous_conc(k,isegment-1)
        Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
      else ! q_dn = 0
        Res_dn(k) = (n_dn*area_dn)*(0.d0 - diffusion)
      endif

      Res(k) = Res_up(k) + Res_dn(k)
    enddo

    pm_well%tran_soln%residual(istart:iend) = &
                      pm_well%tran_soln%residual(istart:iend) + Res(:)
  enddo

  ! ----------------------------------------BOUNDARY-FLUXES------------------

  ! ----- bottom of well -----
  isegment = 1
  Res(:) = 0.d0; Res_up(:) = 0.d0; Res_dn(:) = 0.d0

  area_up = 0.5d0 * (pm_well%well%area(isegment) + pm_well%well%area(isegment+1))
  area_dn = pm_well%well%area(isegment)

  q_up = pm_well%well%ql(isegment)
  q_dn = pm_well%well%ql_bc(1) ! bottom of hole ql

  offset = (isegment-1)*pm_well%nspecies ! = 0
  istart = offset + 1
  iend = offset + pm_well%nspecies

  do ispecies = 1,pm_well%nspecies
    k = ispecies

    ! north surface:
    if (q_up < 0.d0) then ! flow is down the well
      conc = pm_well%well%aqueous_conc(k,isegment+1)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    elseif (q_up > 0.d0) then ! flow is up the well
      conc = pm_well%well%aqueous_conc(k,isegment)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    else ! q_up = 0
      Res_up(k) = (n_up*area_up)*(0.d0 - diffusion)
    endif

    ! south surface:
    if (q_dn < 0.d0) then ! flow is down the well
      conc = pm_well%well%aqueous_conc(k,isegment)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    elseif (q_dn > 0.d0) then ! flow is up the well
      conc = pm_well%well%reservoir%aqueous_conc(k,isegment)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    else ! q_dn = 0
      Res_dn(k) = (n_dn*area_dn)*(0.d0 - diffusion)
    endif

    Res(k) = Res_up(k) + Res_dn(k)
  enddo

  pm_well%tran_soln%residual(istart:iend) = &
                      pm_well%tran_soln%residual(istart:iend) + Res(:)


  ! ----- top of well -----
  isegment = pm_well%well_grid%nsegments
  Res(:) = 0.d0; Res_up(:) = 0.d0; Res_dn(:) = 0.d0

  area_up = pm_well%well%area(isegment)
  area_dn = 0.5d0 * (pm_well%well%area(isegment) + pm_well%well%area(isegment-1))

  q_up = pm_well%well%ql_bc(2) ! top of hole ql
  q_dn = pm_well%well%ql(isegment-1)

  offset = (isegment-1)*pm_well%nspecies
  istart = offset + 1
  iend = offset + pm_well%nspecies

  do ispecies = 1,pm_well%nspecies
    k = ispecies

    ! north surface:
    if (q_up < 0.d0) then ! flow is down the well
      conc = pm_well%well%aqueous_conc_th(ispecies)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    elseif (q_up > 0.d0) then ! flow is up the well
      conc = pm_well%well%aqueous_conc(k,isegment)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    else ! q_up = 0
      Res_up(k) = (n_up*area_up)*(0.d0 - diffusion)
    endif

    ! south surface:
    if (q_dn < 0.d0) then ! flow is down the well
      conc = pm_well%well%aqueous_conc(k,isegment)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    elseif (q_dn > 0.d0) then ! flow is up the well
      conc = pm_well%well%aqueous_conc(k,isegment-1)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    else ! q_dn = 0
      Res_dn(k) = (n_dn*area_dn)*(0.d0 - diffusion)
    endif

    Res(k) = Res_up(k) + Res_dn(k)
  enddo

  pm_well%tran_soln%residual(istart:iend) = &
                      pm_well%tran_soln%residual(istart:iend) + Res(:)

end subroutine PMWellResidualTranFlux

! ************************************************************************** !

subroutine PMWellJacobianFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: local_id
  PetscInt :: iconn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: i,k
  Vec, parameter :: null_vec = tVec(0)
  PetscReal :: Jup(pm_well%nphase,pm_well%nphase), &
               Jdn(pm_well%nphase,pm_well%nphase), &
               Jtop(pm_well%nphase,pm_well%nphase), &
               Jbtm(pm_well%nphase,pm_well%nphase), &
               Jtmp(pm_well%nphase,pm_well%nphase), &
               Jac(pm_well%nphase*pm_well%well_grid%nsegments, &
                   pm_well%nphase*pm_well%well_grid%nsegments)

  pm_well%flow_soln%Jacobian = 0.d0
  Jac = 0.d0
  Jup = 0.d0
  Jtop = 0.d0
  Jbtm = 0.d0
  Jtmp = 0.d0

  select case(pm_well%well%well_model_type)
    case(WELL_MODEL_WIPP_DARCY)
      ! Not expecting ghosting at this time (1D model)
      !if (.not. well_analytical_derivatives) then
      !  call PMWellPerturb(pm_well)
      !endif

      ! Accumulation Term ------------------------------------
      do local_id = 1,pm_well%well_grid%nsegments
        call PMWellAccumDerivative(pm_well,local_id,Jup)
        call PMWellFillJacFlow(pm_well,Jac,Jup,local_id,local_id)
      enddo

      ! Source/Sink Term
      do local_id = 1,pm_well%well_grid%nsegments
        call PMWellSrcSinkDerivative(pm_well,local_id,Jup)
        call PMWellFillJacFlow(pm_well,Jac,Jup,local_id,local_id)
      enddo

      ! Interior Flux Terms -----------------------------------
      do iconn = 1,pm_well%well_grid%nconnections

        local_id_up = iconn
        local_id_dn = iconn+1

        call PMWellFluxDerivative(pm_well,local_id_up,local_id_dn,Jup,Jdn)

        Jtmp = Jup
        call PMWellFillJacFlow(pm_well,Jac,Jtmp,local_id_up,local_id_up)

        Jtmp = Jdn
        call PMWellFillJacFlow(pm_well,Jac,Jtmp,local_id_up,local_id_dn)

        Jup = -Jup
        Jdn = -Jdn
        Jtmp = Jdn
        call PMWellFillJacFlow(pm_well,Jac,Jtmp,local_id_dn,local_id_dn)

        Jtmp = Jup
        call PMWellFillJacFlow(pm_well,Jac,Jtmp,local_id_dn,local_id_up)

      enddo

      ! Boundary Flux Terms -----------------------------------
      local_id = 1
      call PMWellBCFluxDerivative(pm_well,Jtop,Jbtm)
      Jbtm = -Jbtm
      call PMWellFillJacFlow(pm_well,Jac,Jbtm,local_id,local_id)

      local_id = pm_well%well_grid%nsegments
      Jtop = -Jtop
      call PMWellFillJacFlow(pm_well,Jac,Jtop,local_id,local_id)

      !pm_well_ni_count = pm_well_ni_count + 1

    case(WELL_MODEL_FULL_MOMENTUM)

  end select

  do i = 1,pm_well%nphase*pm_well%well_grid%nsegments
    do k = 1,pm_well%nphase*pm_well%well_grid%nsegments
      pm_well%flow_soln%Jacobian(i,k) = Jac(i,k)
    enddo
  enddo

  !pm_well%flow_soln%Jacobian = Jac

end subroutine PMWellJacobianFlow

! ************************************************************************** !
subroutine PMWellFillJacFlow(pm_well,Jac,Jtmp,id1,id2)
  !
  ! Author: Michael Nole
  ! Date: 01/10/2022
  !

  implicit none

  class(pm_well_type) :: pm_well
  PetscReal :: Jtmp(pm_well%nphase,pm_well%nphase), &
               Jac(pm_well%nphase*pm_well%well_grid%nsegments, &
                   pm_well%nphase*pm_well%well_grid%nsegments)
  PetscInt :: id1,id2

  PetscInt :: i,j

  do i = 1,pm_well%nphase
    do j = 1,pm_well%nphase
      Jac((id1-1)*pm_well%nphase+i,(id2-1)*pm_well%nphase+j) = &
      Jac((id1-1)*pm_well%nphase+i,(id2-1)*pm_well%nphase+j) + Jtmp(i,j)
    enddo
  enddo

end subroutine PMWellFillJacFlow

! ************************************************************************** !

subroutine PMWellJacobianTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/14/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscReal :: Jblock(pm_well%nspecies,pm_well%nspecies)
  PetscInt :: k, nspecies
  PetscInt :: jstart, jend

  nspecies = pm_well%nspecies
  pm_well%tran_soln%Jacobian(:,:) = 0.d0

  do k = 1,pm_well%well_grid%nsegments

    Jblock(:,:) = 0.d0

    call PMWellJacTranAccum(pm_well,Jblock,k)

    call PMWellJacTranSrcSink(pm_well,Jblock,k)

    call PMWellJacTranFlux(pm_well,Jblock,k)

    call PMWellJacTranRxn(pm_well,Jblock,k)

    ! place JBlock into full Jac based on isegment
    jstart = (k-1)*nspecies + 1
    jend = jstart + nspecies - 1
    pm_well%tran_soln%Jacobian(jstart:jend,jstart:jend) = Jblock

  enddo

end subroutine PMWellJacobianTran

! ************************************************************************** !

subroutine PMWellJacTranAccum(pm_well,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/14/2022
  !

  implicit none

  class(pm_well_type) :: pm_well
  PetscReal :: Jblock(pm_well%nspecies,pm_well%nspecies)
  PetscInt :: isegment

  PetscReal :: vol_dt
  PetscInt :: istart, iend, ispecies

  ! units of Jac = [m^3-bulk/sec]
  ! units of volume = [m^3-bulk]
  ! units of tran_dt = [sec]

  vol_dt = pm_well%well%volume(isegment)/pm_well%dt_tran

  istart = 1
  iend = pm_well%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + vol_dt

  enddo

end subroutine PMWellJacTranAccum

! ************************************************************************** !

subroutine PMWellJacTranSrcSink(pm_well,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/2022
  !

  implicit none

  class(pm_well_type) :: pm_well
  PetscReal :: Jblock(pm_well%nspecies,pm_well%nspecies)
  PetscInt :: isegment

  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: resr
  PetscInt :: istart, iend, ispecies
  PetscReal :: Qin, Qout
  PetscReal :: SSin, SSout, SS
  PetscReal :: vol, den_avg 

  ! units of Jac = [m^3-bulk/sec]
  ! units of volume = [m^3-bulk]
  ! units of liq%Q = [kmol-liq/sec]
  ! units of FMWH2O = [kg-liq/kmol-liq] 
  ! units of density = [kg-liq/m^3-liq] 
  ! units of Qin = [m^3-liq/sec]
  ! units of SS = [m^3-bulk/sec]

  well => pm_well%well
  resr => pm_well%well%reservoir

  ! From the flow solution:
  ! + Q goes into well from reservoir
  ! - Q goes out of well into reservoir

  vol = pm_well%well%volume(isegment)
  den_avg = 0.5d0*(well%liq%den(isegment)+resr%den_l(isegment))

  ! units of Qin/out = [m^3-liq/sec]
  if (well%liq%Q(isegment) < 0.d0) then ! Q out of well
    Qin = 0.d0
    Qout = well%liq%Q(isegment)*FMWH2O/den_avg
    if (well%liq%s(isegment) < 1.d-40) then
      pm_well%option%io_buffer = 'HINT: The liquid saturation is zero. &
        &Division by zero will occur in PMWellJacTranSrcSink().'
      call PrintMsg(pm_well%option)
    endif
  else ! Q into well
    Qin = 0.d0 !well%liq%Q(isegment)*FMWH2O/den_avg
    Qout = 0.d0
    if (resr%s_l(isegment) < 1.d-40) then
      pm_well%option%io_buffer = 'HINT: The liquid saturation is zero. &
        &Division by zero will occur in PMWellJacTranSrcSink().'
      call PrintMsg(pm_well%option)
    endif
  endif

  SSin = Qin / (resr%e_por(isegment)*resr%s_l(isegment))    ! [m3-bulk/sec]
  SSout = Qout / (well%phi(isegment)*well%liq%s(isegment))  ! [m3-bulk/sec]
  SS = SSin + SSout

  istart = 1
  iend = pm_well%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + vol*(SS/vol)

  enddo

end subroutine PMWellJacTranSrcSink

! ************************************************************************** !

subroutine PMWellJacTranFlux(pm_well,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/2022
  !

  implicit none

  class(pm_well_type) :: pm_well
  PetscReal :: Jblock(pm_well%nspecies,pm_well%nspecies)
  PetscInt :: isegment

  type(well_type), pointer :: well
  PetscInt :: istart, iend, ispecies
  PetscInt :: n_up, n_dn
  PetscReal :: d_diffusion_dM
  PetscReal :: J_up, J_dn
  PetscReal :: area_up, area_dn
  PetscReal :: sat_up, sat_dn, por_up, por_dn
  PetscReal :: u_up, u_dn

  ! units of Jac = [m^3-bulk/sec]
  ! area in [m2-bulk]
  ! q in [m3-liq/m2-bulk-sec]
  ! u in [m-liq/sec]
  ! sat in [m2-liq/m2-void] 

  well => pm_well%well

  ! NOTE: The up direction is towards well top, and the dn direction is
  !       towards the well bottom.
  !       +q flows up the well
  !       -q flows down the well

  n_dn = +1
  n_up = -1

  d_diffusion_dM = 0.d0 ! for now, since WIPP has no diffusion

  if ((isegment > 1) .and. (isegment < pm_well%well_grid%nsegments)) then
  ! ----------------------------------------INTERIOR-FLUXES------------------

    ! define face values with arithmetic averages:
    area_up = 0.5d0 * (well%area(isegment) + well%area(isegment+1))
    area_dn = 0.5d0 * (well%area(isegment) + well%area(isegment-1))
    por_up = 0.5d0 * (well%phi(isegment) + well%phi(isegment+1)) 
    por_dn = 0.5d0 * (well%phi(isegment) + well%phi(isegment-1))    
    sat_up = 0.5d0 * (well%liq%s(isegment) + well%liq%s(isegment+1)) 
    sat_dn = 0.5d0 * (well%liq%s(isegment) + well%liq%s(isegment-1)) 

    u_up = well%ql(isegment)/(sat_up*por_up)
    u_dn = well%ql(isegment-1)/(sat_dn*por_dn)

  ! ----------------------------------------BOUNDARY-FLUXES------------------
  else if (isegment == 1) then
    ! ----- bottom of well -----

    ! define face values with arithmetic averages:
    area_up = 0.5d0 * (well%area(isegment) + well%area(isegment+1))
    area_dn = well%area(isegment)
    por_up = 0.5d0 * (well%phi(isegment) + well%phi(isegment+1)) 
    por_dn = well%phi(isegment)   
    sat_up = 0.5d0 * (well%liq%s(isegment) + well%liq%s(isegment+1)) 
    sat_dn = well%liq%s(isegment)

    u_up = well%ql(isegment)/(sat_up*por_up)
    u_dn = well%ql_bc(1)/(sat_dn*por_dn)        ! bottom of hole ql = ql_bc(1)

  else if (isegment == pm_well%well_grid%nsegments) then
    ! ----- top of well -----

    ! define face values with arithmetic averages:
    area_up = well%area(isegment)
    area_dn = 0.5d0 * (well%area(isegment) + well%area(isegment-1))
    por_up = well%phi(isegment)
    por_dn = 0.5d0 * (well%phi(isegment) + well%phi(isegment-1))    
    sat_up = well%liq%s(isegment)
    sat_dn = 0.5d0 * (well%liq%s(isegment) + well%liq%s(isegment-1)) 

    u_up = well%ql_bc(2)/(sat_up*por_up)      ! top of hole ql = ql_bc(2)
    u_dn = well%ql(isegment-1)/(sat_dn*por_dn)

  endif

  ! north surface:
  J_up = (n_up*area_up)*(u_up - d_diffusion_dM)

  ! south surface:
  J_dn = (n_dn*area_dn)*(u_dn - d_diffusion_dM)

  istart = 1
  iend = pm_well%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + (J_up + J_dn)

  enddo

end subroutine PMWellJacTranFlux

! ************************************************************************** !

subroutine PMWellJacTranRxn(pm_well,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/2022
  !

  implicit none

  class(pm_well_type) :: pm_well
  PetscReal :: Jblock(pm_well%nspecies,pm_well%nspecies)
  PetscInt :: isegment

  PetscInt :: istart, iend, ispecies, parent_id
  PetscReal :: vol

  ! units of Jac = [m^3-bulk/sec]
  ! units of volume = [m^3-bulk]
  ! decay_rate in [1/sec]

  vol = pm_well%well%volume(isegment)

  istart = 1
  iend = pm_well%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + &
                               (vol*pm_well%well%species_decay_rate(ispecies))

    parent_id = pm_well%well%species_parent_id(ispecies)
    if (parent_id > 0) then
      Jblock(ispecies,parent_id) = Jblock(ispecies,parent_id) - &
                        (vol*pm_well%well%species_parent_decay_rate(ispecies))
    endif

  enddo

end subroutine PMWellJacTranRxn

! ************************************************************************** !

subroutine PMWellPreSolve(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_type) :: this

  ! placeholder

end subroutine PMWellPreSolve

! ************************************************************************** !

subroutine PMWellPreSolveFlow(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_type) :: pm_well

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: cur_time, cur_time_converted
  PetscReal :: dt_converted

  pm_well%flow_soln%not_converged = PETSC_TRUE
  pm_well%flow_soln%converged = PETSC_FALSE

  cur_time = pm_well%option%time + pm_well%cumulative_dt_flow
  cur_time_converted = cur_time/pm_well%output_option%tconv
  dt_converted = pm_well%dt_flow/pm_well%output_option%tconv

  if (pm_well%print_output) then
    write(out_string,'(" FLOW Step ",i6,"   Time =",1pe12.5,"   Dt =", &
                       1pe12.5," ",a4)') &
                     (pm_well%flow_soln%n_steps+1),cur_time_converted, &
                     dt_converted,pm_well%output_option%tunit
    call PrintMsg(pm_well%option,out_string)
  endif

end subroutine PMWellPreSolveFlow

! ************************************************************************** !

subroutine PMWellPreSolveTran(pm_well,master_dt)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/22/2022
  !

  implicit none

  class(pm_well_type) :: pm_well
  PetscReal :: master_dt

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: cur_time, cur_time_converted
  PetscReal :: dt_converted

  pm_well%tran_soln%not_converged = PETSC_TRUE
  pm_well%tran_soln%converged = PETSC_FALSE

  cur_time = pm_well%option%time + pm_well%cumulative_dt_tran
  pm_well%tran_soln%tran_time = cur_time

  cur_time_converted = cur_time/pm_well%output_option%tconv
  dt_converted = pm_well%dt_tran/pm_well%output_option%tconv

  write(out_string,'(" WELL TRAN Step ",i6,"   Time =",1pe12.5,"   Dt =", &
                   1pe12.5," ",a4)') &
                   (pm_well%tran_soln%n_steps+1),cur_time_converted, &
                   dt_converted,pm_well%output_option%tunit

  call PrintMsg(pm_well%option,out_string)

end subroutine PMWellPreSolveTran

! ************************************************************************** !

subroutine PMWellSolve(this,time,ierr)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021
  !

  implicit none

  class(pm_well_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: curr_time, curr_time_converted

  curr_time = this%option%time + this%cumulative_dt_tran
  curr_time_converted = curr_time/this%output_option%tconv

  ierr = 0 ! If this is not set to zero, TS_STOP_FAILURE occurs if the solve
           ! routines are not entered, either due to an inactive well or due
           ! to being on a process that doesn't contain a well segment.          


  if (Initialized(this%intrusion_time_start) .and. &
      (curr_time < this%intrusion_time_start)) then
    write(out_string,'(" Inactive.    Time =",1pe12.5," ",a4)') &
          curr_time_converted,this%output_option%tunit
    call PrintMsg(this%option,out_string)
    return
  endif

  if (this%update_for_wippflo_qi_coupling) then
    write(out_string,'(" FLOW Step          Quasi-implicit wellbore flow &
                      &coupling is being used.")')
    call PrintMsg(this%option,out_string)
    this%update_for_wippflo_qi_coupling = PETSC_FALSE
  elseif (this%flow_coupling == FULLY_IMPLICIT_WELL) then
    write(out_string,'(" FLOW Step          Fully-implicit wellbore flow &
                      &coupling is being used.")')
    call PrintMsg(this%option,out_string)
  else
    call PMWellSolveFlow(this,-999,ierr)
  endif

  !Debugging
  !call MPI_Barrier(this%option%comm%communicator,ierr);CHKERRQ(ierr)
  if (this%transport) then
    write(out_string,'(" TRAN Step          Quasi-implicit wellbore &
                     &transport coupling is being used.")')
    call PrintMsg(this%option,out_string)
    this%tran_soln%prev_soln%aqueous_conc = this%well%aqueous_conc
    this%tran_soln%prev_soln%aqueous_mass = this%well%aqueous_mass
  endif

end subroutine PMWellSolve

! ************************************************************************** !

subroutine PMWellSolveFlow(pm_well,perturbation_index,ierr)
  !
  ! Author: Michael Nole
  ! Date: 12/01/2021
  !

  use Option_module
  use EOS_Water_module
  use EOS_Gas_module
  use SCO2_Aux_module, only : fmw_comp

  implicit none

  class(pm_well_type) :: pm_well
  PetscInt :: perturbation_index
  PetscErrorCode :: ierr

  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: reservoir
  type(option_type), pointer :: option
  type(well_soln_flow_type), pointer :: flow_soln
  character(len=MAXSTRINGLENGTH) :: out_string
  PetscLogDouble :: log_start_time, log_end_time
  PetscInt :: n_iter,ts_cut,easy_converge_count
  PetscInt :: istart, iend
  PetscReal :: res(pm_well%flow_soln%ndof)
  PetscReal :: res_fixed(pm_well%flow_soln%ndof*pm_well%well_grid%nsegments)
  PetscReal :: Q_liq(pm_well%well_grid%nsegments,pm_well%well_grid%nsegments),&
               Q_gas(pm_well%well_grid%nsegments,pm_well%well_grid%nsegments)
  PetscReal :: v_darcy
  PetscBool :: at_steady_state, upwind
  PetscReal :: ss_check_p(pm_well%well_grid%nsegments,2), &
               ss_check_s(pm_well%well_grid%nsegments,2)
  PetscReal :: mobility
  PetscReal :: area, mass_conserved_liq, mass_conserved_gas
  PetscInt :: i
  PetscInt :: ss_step_count, steps_to_declare_ss
  PetscReal :: pl, pl0, pg, pg0, temperature, den_mol
  PetscReal :: rho_kg_liq, rho_zero_liq, rho_one_liq, &
               rho_kg_gas, rho_zero_gas, rho_one_gas
  PetscReal :: dist_z, delta_z
  PetscReal :: gravity, den_ave
  PetscInt :: num_iteration
  PetscReal :: dummy, dummy2
  PetscReal :: aux(2)
  PetscReal, parameter :: threshold_p = 0.d0

  option => pm_well%realization%option

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  flow_soln => pm_well%flow_soln

  ierr = 0
  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

  ts_cut = 0
  easy_converge_count = 0

  pm_well%cumulative_dt_flow = 0.d0
  flow_soln%converged = PETSC_FALSE
  flow_soln%not_converged = PETSC_TRUE

  ss_check_p(:,1) = pm_well%well%pl(:)
  ss_check_s(:,1) = pm_well%well%gas%s(:)
  at_steady_state = PETSC_FALSE
  ss_step_count = 0
  steps_to_declare_ss = 10

  gravity = option%gravity(Z_DIRECTION)

  ! update well index
  call PMWellComputeWellIndex(pm_well)

  if (pm_well%well%well_model_type == WELL_MODEL_HYDROSTATIC) then
    ! Compute hydrostatic pressure relative to bottom-hole pressure,
    ! then compute well fluxes from those pressures.

    if (perturbation_index > 0) then
      well => pm_well%well_pert(perturbation_index)
    else
      well => pm_well%well
    endif

    reservoir => well%reservoir

    aux(:) = 0.d0
    area = 0.d0
    v_darcy = 0.d0
    Q_liq = 0.d0
    Q_gas = 0.d0
    well%liq%Q = 0.d0
    well%gas%Q = 0.d0

    temperature = well%temp(1)
    ! Start from bottom hole pressure
    pl = well%bh_p
    ! No capillarity in the well
    pg = pl

    !MAN: pure brine and pure gas columns for now, needs updating.
    aux(1) = well%liq%xmass(1,option%salt_id)
    call EOSWaterDensityExt(temperature,pl, &
                            aux,rho_kg_liq,dummy,ierr)
    pl0 = pl

    call EOSGasDensity(temperature,pg, &
                       den_mol,dummy,dummy2,ierr)
    rho_kg_gas = den_mol * fmw_comp(TWO_INTEGER)
    
    pg0 = pg
  
    ! compute pressures above datum
    dist_z = 0.d0
    rho_zero_liq = rho_kg_liq
    rho_zero_gas = rho_kg_gas

    do i = 1,pm_well%well_grid%nsegments

      ! Compute well cell pressures based off of bottom segment pressure
      temperature = well%temp(i)
      aux(1) = well%liq%xmass(i,option%salt_id)
      call EOSWaterDensityExt(temperature,pl0, &
                              aux,rho_kg_liq,dummy,ierr)
      call EOSGasDensity(temperature,pg0, &
                       den_mol,dummy,dummy2,ierr)
      rho_kg_gas = den_mol * fmw_comp(TWO_INTEGER)

      num_iteration = 0
      delta_z = pm_well%well_grid%h(i)%z - pm_well%well_grid%h(i-1)%z
      do
        pl = pl0 + 0.5d0*(rho_kg_liq+rho_zero_liq) * &
              gravity * delta_z
        pg = pg0 + 0.5d0*(rho_kg_gas+rho_zero_gas) * &
              gravity * delta_z
        call EOSWaterDensityExt(temperature,pl,aux,rho_one_liq,dummy,ierr)
        call EOSGasDensity(temperature,pg,den_mol,dummy,dummy2,ierr)
        rho_one_gas = den_mol * fmw_comp(TWO_INTEGER)

        ! STOMP doesn't iterate
        exit

        if (dabs(rho_kg_liq-rho_one_liq) < 1.d-5 .and. &
            dabs(rho_kg_gas-rho_one_gas) < 1.d-5) exit
        rho_kg_liq = rho_one_liq
        rho_kg_gas = rho_one_gas
        num_iteration = num_iteration + 1
        if (num_iteration > 100) then
          option%io_buffer = 'Hydrostatic iteration failed to converge in well model'
          call PrintErrMsg(option)
        endif
      enddo
      rho_zero_liq = rho_kg_liq
      rho_zero_gas = rho_kg_gas
      well%pl(i) = pl
      well%pg(i) = pg
      well%liq%den(i) = rho_kg_liq
      well%gas%den(i) = rho_kg_gas
      pl0 = pl
      pg0 = pg

      if (pg < 0.d0) then
        option%io_buffer = 'BHP cannot support the gas column: &
                            &hydrostatic well pressure calcuation &
                            &results in negative gas pressure.'
        call PrintErrMsg(option)
      endif
    enddo

    ! Update well properties based off new pressures
    select case(option%iflowmode)
      case(SCO2_MODE)
        call PMWellUpdatePropertiesSCO2Flow(pm_well,well,option)
    end select

    ! Compute fluxes in/out of well
    well%liq%Q = 0.d0
    well%gas%Q = 0.d0

    do i = 1,pm_well%well_grid%nsegments
      
      if (well%WI(i) == 0) cycle

      if (well%th_qg > 0.d0) then
        ! Rate-controlled gas injection well. Can potentially have
        ! under-pressure in some well segments.
        upwind = reservoir%p_g(i) > well%pg(i)
        if (upwind) then
          mobility = reservoir%kr_g(i)/reservoir%visc_g(i)
          den_ave = reservoir%den_g(i)
        else
          mobility = 1.d0 / well%gas%visc(i)
          den_ave = well%gas%den(i)
        endif

        ! Flowrate in kg/s
        well%gas%Q(i) = den_ave*mobility*well%WI(i)* &
                        (reservoir%p_g(i)-well%pg(i))
      elseif (well%th_ql > 0.d0) then
        ! Rate-controlled water injection well
        upwind = reservoir%p_l(i) > well%pl(i)
        if (upwind) then
          mobility = reservoir%kr_l(i)/reservoir%visc_l(i)
          den_ave = reservoir%den_l(i)
        else
          mobility = 1.d0 / well%liq%visc(i)
          den_ave = well%liq%den(i)
        endif
        ! Flowrate in kg/s
        well%liq%Q(i) = den_ave*mobility*well%WI(i)* &
                        (reservoir%p_l(i)-well%pl(i))
      elseif (well%th_qg < 0.d0 .or. well%th_ql < 0.d0) then
        ! Extraction well
        upwind = reservoir%p_l(i) > well%pl(i)
        if (upwind) then
          mobility = reservoir%kr_l(i)/reservoir%visc_l(i)
          den_ave = reservoir%den_l(i)
        else
          mobility = dabs(well%th_ql)/ (dabs(well%th_ql + well%th_qg))/ &
                     well%liq%visc(i)
          den_ave = well%liq%den(i)
        endif
        ! Flowrate in kg/s
        well%liq%Q(i) = den_ave*mobility*well%WI(i)* &
                        (reservoir%p_l(i)-well%pl(i))

        upwind = reservoir%p_g(i) > well%pg(i)
        if (upwind) then
          mobility = reservoir%kr_g(i)/reservoir%visc_g(i)
          den_ave = reservoir%den_g(i)
        else
          mobility = dabs(well%th_qg)/ (dabs(well%th_ql + well%th_qg))/ &
                     well%gas%visc(i)
          den_ave = well%gas%den(i)
        endif

        ! Flowrate in kg/s
        well%gas%Q(i) = den_ave*mobility*well%WI(i)* &
                        (reservoir%p_g(i)-well%pg(i))
      endif
    enddo

    ! Should equal the total flux out the top of the domain
    mass_conserved_liq = sum(well%liq%Q)
    mass_conserved_gas = sum(well%gas%Q)

    pm_well%cumulative_dt_flow = pm_well%realization%option%flow_dt

    ! Update transport
    ! flow_soln%n_steps = flow_soln%n_steps + 1

  else

  do while (pm_well%cumulative_dt_flow < pm_well%realization%option%flow_dt)

    ! update the well src/sink Q vector at start of time step
    call PMWellUpdateWellQ(pm_well%well,pm_well%well%reservoir)

    call PMWellPreSolveFlow(pm_well)

    ! Fixed accumulation term
    res_fixed = 0.d0
    res = 0.d0
    do i = 1,pm_well%well_grid%nsegments
      call PMWellAccumulationFlow(pm_well,pm_well%well,i,res)
      istart = flow_soln%ndof*(i-1)+1
      iend = flow_soln%ndof*i
      res_fixed(istart:iend) = -1.d0 * res * pm_well%dt_flow
    enddo

    n_iter = 0

    do while (flow_soln%not_converged)

      if (n_iter > (flow_soln%max_iter-1)) then
        flow_soln%cut_timestep = PETSC_TRUE
        if (pm_well%print_output) then
          out_string = ' Maximum number of FLOW Newton iterations reached. &
                        &Cutting timestep!'
          call PrintMsg(pm_well%option,out_string)
        endif
        call PMWellCutTimestepFlow(pm_well)
        n_iter = 0
        ts_cut = ts_cut + 1
        easy_converge_count = 0

        if (ss_step_count > 2 .and. pm_well%ss_check) then
          at_steady_state = PETSC_TRUE
          pm_well%cumulative_dt_flow = pm_well%realization%option%flow_dt
          WRITE(out_string,'(" PM Well FLOW convergence declared due to &
            &automatic time step control criterion. ")')
          call PrintMsg(pm_well%option,out_string)
        endif


        exit
      endif
      if (ts_cut > flow_soln%max_ts_cut) then
        pm_well%realization%option%io_buffer = &
          ' Maximum timestep cuts reached in PM Well FLOW. Solution has not &
           &converged. Exiting.'
        if (pm_well%print_well) then
          call PMWellOutput(pm_well)  
        endif
        call PrintErrMsg(pm_well%realization%option)
      endif

      if (pm_well%dt_flow <= pm_well%min_dt_flow) then
        pm_well%well_force_ts_cut = 1
        call PMWellCopyReservoir(pm_well%well%reservoir_save, &
                                 pm_well%well%reservoir,&
                                 pm_well%transport)
        return
      endif

      flow_soln%residual = 0.d0
      flow_soln%residual = res_fixed / pm_well%dt_flow

      easy_converge_count = easy_converge_count + 1

      call PMWellNewtonFlow(pm_well)

      if (pm_well%well_force_ts_cut > 0) return

      call PMWellCheckConvergenceFlow(pm_well,n_iter,res_fixed)

    enddo

    if (easy_converge_count > 4 ) then
      if (pm_well%cumulative_dt_flow + pm_well%dt_flow * &
          flow_soln%ts_cut_factor < pm_well%realization%option%flow_dt) then
        pm_well%dt_flow = pm_well%dt_flow * flow_soln%ts_cut_factor
        flow_soln%cut_timestep = PETSC_FALSE
      endif
    endif

    if (pm_well%cumulative_dt_flow + pm_well%dt_flow > &
        pm_well%realization%option%flow_dt) then
      pm_well%dt_flow = pm_well%realization%option%flow_dt - pm_well%cumulative_dt_flow
    endif

    ts_cut = 0
    flow_soln%n_steps = flow_soln%n_steps + 1

    ! Check if we're at steady-state

    ! TOUGH way:
    if (n_iter == 1 .and. flow_soln%converged) then
      ss_step_count = ss_step_count+1
    else
      ss_step_count = 0
    endif

    if (pm_well%ss_check) then
      if (ss_step_count >= steps_to_declare_ss) then
        at_steady_state = PETSC_TRUE
        pm_well%cumulative_dt_flow = pm_well%realization%option%flow_dt
      endif
    endif

    ! Other way:
    !if (pm_well%ss_check .and. flow_soln%converged) then
    !  ss_check_p(:,2) = pm_well%well%pl(:)
    !  ss_check_s(:,2) = pm_well%well%gas%s(:)

    !  dpdt = (ss_check_p(:,2) - ss_check_p(:,1)) / pm_well%dt_flow
    !  dsdt = (ss_check_s(:,2) - ss_check_s(:,1)) / pm_well%dt_flow

    !  if (maxval(abs(dpdt)) < eps_p) then
    !    if (maxval(abs(dsdt)) < eps_s) then
    !      ss_step_count = ss_step_count + 1
    !      if (ss_step_count > steps_to_declare_ss) at_steady_state = PETSC_TRUE
    !    else
    !      ss_step_count = 0
    !    endif
    !  endif
    !  if (at_steady_state) then
    !    pm_well%cumulative_dt_flow = pm_well%realization%option%flow_dt
    !  endif
    !  ss_check_p(:,1) = pm_well%well%pl(:)
    !  ss_check_s(:,1) = pm_well%well%gas%s(:)
    !endif

  enddo

  endif
  
  if (pm_well%flow_coupling /= FULLY_IMPLICIT_WELL .and. &
     pm_well%flow_coupling /= QUASI_IMPLICIT_WELL) then
    call PMWellPostSolveFlow(pm_well)
    call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
  endif

end subroutine PMWellSolveFlow

! ************************************************************************** !

subroutine PMWellSolveTran(pm_well,ierr)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/22/2022
  !

  implicit none

  class(pm_well_type) :: pm_well
  PetscErrorCode :: ierr

  type(well_soln_tran_type), pointer :: soln
  character(len=MAXSTRINGLENGTH) :: out_string
  PetscLogDouble :: log_start_time, log_end_time
  PetscReal :: res_fixed(pm_well%tran_soln%ndof*pm_well%well_grid%nsegments)
  PetscReal :: master_dt
  PetscInt :: n_iter, ts_cut
  PetscInt :: istart, iend
  PetscInt :: k

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  soln => pm_well%tran_soln

  ierr = 0
  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

  ts_cut = 0

  pm_well%cumulative_dt_tran = 0.d0
  soln%converged = PETSC_FALSE
  soln%not_converged = PETSC_TRUE

  master_dt = pm_well%option%tran_dt

  do while (pm_well%cumulative_dt_tran < master_dt)

    call PMWellPreSolveTran(pm_well,master_dt)

    n_iter = 0

    do while (soln%not_converged)
      if (n_iter > (soln%max_iter-1)) then
        soln%cut_timestep = PETSC_TRUE
        soln%cut_ts_flag = PETSC_TRUE
        out_string = ' Maximum number of TRAN Newton iterations reached. &
                      &Cutting timestep!'
        call PrintMsg(pm_well%option,out_string)
        call PMWellCutTimestepTran(pm_well)
        return 
      endif
      if (ts_cut > soln%max_ts_cut) then
        pm_well%realization%option%io_buffer = &
          ' Maximum timestep cuts reached in PM Well TRAN. Solution has not &
           &converged. Exiting.'
        if (pm_well%print_well) then
          call PMWellOutput(pm_well)
        endif
        call PrintErrMsg(pm_well%realization%option)
      endif

      soln%residual = 0.d0
      if (any(pm_well%option%myrank == pm_well%well_grid%h_rank_id)) then
        ! Get fixed accumulation term (not yet divided by dt)
        do k = 1,pm_well%well_grid%nsegments
          istart = soln%ndof*(k-1)+1
          iend = soln%ndof*k
          call PMWellAccumulationTran(pm_well,k,res_fixed(istart:iend))
        enddo
        soln%residual = res_fixed / pm_well%dt_tran
      endif
      call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call PMWellNewtonTran(pm_well,n_iter)
      call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call PMWellCheckConvergenceTran(pm_well,n_iter,res_fixed)

    enddo

    ! try to increase the time step, if possible
    if (soln%converged) then
      pm_well%dt_tran = soln%ts_ramp_factor * pm_well%dt_tran
    endif 
    if (pm_well%dt_tran > master_dt) then
      pm_well%dt_tran = master_dt
    endif

    ! if pm_well next time step will overstep master_dt, then correct it
    if (pm_well%cumulative_dt_tran + pm_well%dt_tran > master_dt) then
      pm_well%dt_tran = master_dt - pm_well%cumulative_dt_tran
    endif

    soln%n_steps = soln%n_steps + 1

  enddo

  call PMWellPostSolveTran(pm_well)

  call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

end subroutine PMWellSolveTran

! ************************************************************************** !

subroutine PMWellUpdateSolutionFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 01/21/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: i
  PetscInt :: idof
  PetscReal, parameter :: MIN_SAT = 1.d-6
  PetscReal, parameter :: MAX_SAT = 0.99999

  select case(pm_well%well%well_model_type)
    case(WELL_MODEL_WIPP_DARCY)
      do i = 1,pm_well%well_grid%nsegments
        idof = pm_well%flow_soln%ndof*(i-1)
        pm_well%well%pl(i) = pm_well%well%pl(i) +  &
                             pm_well%flow_soln%update(idof+1)
        idof = pm_well%flow_soln%ndof*i
        if (pm_well%well%gas%s(i) + pm_well%flow_soln%update(idof) < &
            MIN_SAT) then
          pm_well%flow_soln%update(idof) = MIN_SAT - pm_well%well%gas%s(i)
          pm_well%well%gas%s(i) = MIN_SAT
        elseif (pm_well%well%gas%s(i) + pm_well%flow_soln%update(idof) > &
                MAX_SAT) then
          pm_well%flow_soln%update(idof) = MAX_SAT - pm_well%well%gas%s(i)
          pm_well%well%gas%s(i) = MAX_SAT
        else
          pm_well%well%gas%s(i) = pm_well%well%gas%s(i) + &
                                  pm_well%flow_soln%update(idof)
        endif

      enddo
      call PMWellUpdatePropertiesWIPPFlow(pm_well,pm_well%well, &
                     pm_well%realization%patch%characteristic_curves_array, &
                     pm_well%realization%option)
  end select
end subroutine PMWellUpdateSolutionFlow

! ************************************************************************** !

subroutine PMWellUpdateSolutionTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/22
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: isegment
  PetscInt :: offset, istart, iend
  PetscInt :: nspecies 
  PetscReal :: vol

  ! update in [mol/m3-bulk]
  ! volume in [m3-bulk]

  nspecies = pm_well%nspecies 

  do isegment = 1,pm_well%well_grid%nsegments

    offset = (isegment-1)*nspecies
    istart = offset + 1
    iend = offset + nspecies

    vol = pm_well%well%volume(isegment)

    pm_well%well%aqueous_mass(1:nspecies,isegment) = &              ! [mol]
                           pm_well%well%aqueous_mass(1:nspecies,isegment) + &
                           (pm_well%tran_soln%update(istart:iend) * vol)
    pm_well%well%aqueous_conc(1:nspecies,isegment) = &
          pm_well%well%aqueous_mass(1:nspecies,isegment) / &        ! [mol]
          (pm_well%well%phi(isegment)*vol*pm_well%well%liq%s(isegment)) ! [m3-liq]

  enddo

end subroutine PMWellUpdateSolutionTran

! ************************************************************************** !

subroutine PMWellCutTimestepFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 01/24/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  select case(pm_well%well%well_model_type)
    case(WELL_MODEL_WIPP_DARCY)
      ! could make this smarter or call smarter timestepping routines
      pm_well%dt_flow = pm_well%dt_flow / pm_well%flow_soln%ts_cut_factor
      pm_well%dt_flow = max(pm_well%dt_flow,pm_well%min_dt_flow)
      pm_well%well%pl = pm_well%flow_soln%prev_soln%pl
      pm_well%well%gas%s = pm_well%flow_soln%prev_soln%sg
      call PMWellUpdatePropertiesWIPPFlow(pm_well,pm_well%well, &
                     pm_well%realization%patch%characteristic_curves_array, &
                     pm_well%realization%option)
    case(WELL_MODEL_HYDROSTATIC)
      pm_well%well%bh_p = pm_well%flow_soln%prev_soln%bh_p
  end select

end subroutine PMWellCutTimestepFlow

! ************************************************************************** !

subroutine PMWellCutTimestepTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  implicit none

  class(pm_well_type) :: pm_well

  pm_well%well%aqueous_mass = pm_well%tran_soln%prev_soln%aqueous_mass
  pm_well%well%aqueous_conc = pm_well%tran_soln%prev_soln%aqueous_conc
  call PMWellUpdatePropertiesTran(pm_well)

end subroutine PMWellCutTimestepTran

! ************************************************************************** !

subroutine PMWellNewtonFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 01/20/2022
  !

  use Utility_module

  implicit none

  class(pm_well_type) :: pm_well

  PetscReal :: identity(pm_well%nphase*pm_well%well_grid%nsegments,&
                        pm_well%nphase*pm_well%well_grid%nsegments)
  PetscReal :: new_dx(pm_well%nphase*pm_well%well_grid%nsegments)
  PetscInt :: indx(pm_well%nphase*pm_well%well_grid%nsegments)
  PetscInt :: i,j
  PetscInt :: d

  call PMWellUpdateWellQ(pm_well%well,pm_well%well%reservoir)

  call PMWellPerturb(pm_well)

  call PMWellResidualFlow(pm_well)

  call PMWellJacobianFlow(pm_well)

  select case (pm_well%well%well_model_type)
  !--------------------------------------
  case(WELL_MODEL_HYDROSTATIC)
    ! No capillarity yet
    ! pm_well%well%pl(:) = pm_well%well%bh_p
    ! pm_well%well%pg(:) = pm_well%well%bh_p
  !--------------------------------------
  case(WELL_MODEL_WIPP_DARCY)

    do i = 1,pm_well%nphase*pm_well%well_grid%nsegments
      do j = 1,pm_well%nphase*pm_well%well_grid%nsegments
        if (i==j) then
          identity(i,j) = 1.d0
        else
          identity(i,j) = 0.d0
        endif
      enddo
    enddo
    call LUDecomposition(pm_well%flow_soln%Jacobian,pm_well%nphase*pm_well% &
                         well_grid%nsegments,indx,d)
    call LUBackSubstitution(pm_well%flow_soln%Jacobian, &
                            pm_well%nphase*pm_well%well_grid%nsegments,&
                            indx,pm_well%flow_soln%residual)
    new_dx = -1.d0 * pm_well%flow_soln%residual


    do i = 1,pm_well%well_grid%nsegments
      if (dabs(new_dx(i)) > 1.d15) then
        pm_well%well_force_ts_cut = 1
        return        
      endif
      if (isnan(new_dx(i))) then
        pm_well%well_force_ts_cut = 1
        return
      endif
    enddo

    pm_well%flow_soln%update = new_dx

    call PMWellUpdateSolutionFlow(pm_well)

  !--------------------------------------
  end select

end subroutine PMWellNewtonFlow

! ************************************************************************** !

subroutine PMWellNewtonTran(pm_well,n_iter)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  use Utility_module

  implicit none

  class(pm_well_type) :: pm_well
  PetscInt :: n_iter

  PetscInt :: nm, dummy
  PetscInt :: indx(pm_well%nspecies*pm_well%well_grid%nsegments)

  if (.not. any(pm_well%option%myrank == pm_well%well_grid%h_rank_id)) return

  nm = pm_well%nspecies * pm_well%well_grid%nsegments

  ! at this time, the tran_soln%residual vector has been zero'd out and has
  ! been loaded with the fixed accumulation divided by current dt

  call PMWellResidualTran(pm_well)

  call PMWellJacobianTran(pm_well)

  ! J dx = -R     => dx = J^(-1)(-R)
  ! [m3-bulk/sec] dx = -[mol/sec]
  ! dx in [mol/m3-bulk]

  call LUDecomposition(pm_well%tran_soln%Jacobian,nm,indx,dummy)

  call LUBackSubstitution(pm_well%tran_soln%Jacobian,nm,indx, &
                          pm_well%tran_soln%residual)

  pm_well%tran_soln%update = +1.0d0 * pm_well%tran_soln%residual ! [mol/m3-bulk]

  call PMWellUpdateSolutionTran(pm_well)

end subroutine PMWellNewtonTran

! ************************************************************************** !

subroutine PMWellPostSolve(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_type) :: this

  ! placeholder

end subroutine PMWellPostSolve

! ************************************************************************** !

subroutine PMWellPostSolveFlow(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_type) :: pm_well

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: cur_time_converted

  cur_time_converted = pm_well%option%time/pm_well%output_option%tconv

  WRITE(out_string,'(" PM Well FLOW Step Complete!    Time=",1pe12.5," &
                    &",a4,"Total Newton Its =",i8)') &
                    cur_time_converted,pm_well%output_option%tunit, &
                    pm_well%flow_soln%n_newton
  call PrintMsg(pm_well%option,out_string)
  call PrintMsg(pm_well%option,'')

end subroutine PMWellPostSolveFlow

! ************************************************************************** !

subroutine PMWellPostSolveTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/22/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscReal :: cur_time, cur_time_converted

  cur_time = pm_well%option%time + pm_well%option%tran_dt
  pm_well%tran_soln%tran_time = cur_time
  cur_time_converted = cur_time/pm_well%output_option%tconv

end subroutine PMWellPostSolveTran

! ************************************************************************** !

subroutine PMWellCheckConvergenceFlow(pm_well,n_iter,fixed_accum)
  !
  ! Checks flow solution convergence against prescribed tolerances.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 01/20/2022
  !

  implicit none

  class(pm_well_type) :: pm_well
  PetscInt :: n_iter
  PetscReal :: fixed_accum(pm_well%flow_soln%ndof*pm_well%well_grid%nsegments)

  PetscReal, parameter :: zero_saturation = 1.d-15
  PetscReal, parameter :: zero_accumulation = 1.d-15

  type(well_soln_flow_type), pointer :: flow_soln
  character(len=MAXSTRINGLENGTH) :: out_string
  character(len=MAXSTRINGLENGTH) :: rsn_string
  PetscBool :: cnvgd_due_to_residual(pm_well%well_grid%nsegments* &
                                     pm_well%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_abs_res(pm_well%well_grid%nsegments* &
                                    pm_well%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_scaled_res(pm_well%well_grid%nsegments* &
                                       pm_well%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_update(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update_p(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update_s(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update_p(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update_s(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_on_pressure(pm_well%well_grid%nsegments)
  PetscBool :: cnvgd_on_saturation(pm_well%well_grid%nsegments)
  PetscReal :: update_p(pm_well%well_grid%nsegments) ! liquid pressure
  PetscReal :: update_s(pm_well%well_grid%nsegments) ! gas saturation
  PetscReal :: temp_real
  PetscReal :: max_scaled_residual,max_absolute_residual
  PetscReal :: max_relative_update_p,max_relative_update_s
  PetscReal :: max_absolute_update_p,max_absolute_update_s
  PetscInt :: loc_max_scaled_residual,loc_max_abs_residual
  PetscInt :: loc_max_rel_update_p,loc_max_rel_update_s
  PetscInt :: loc_max_abs_update_p,loc_max_abs_update_s
  PetscInt :: idof
  PetscInt :: k

  flow_soln => pm_well%flow_soln

  n_iter = n_iter + 1
  flow_soln%n_newton = flow_soln%n_newton + 1

  cnvgd_due_to_residual = PETSC_FALSE
  cnvgd_due_to_abs_res = PETSC_FALSE
  cnvgd_due_to_scaled_res = PETSC_FALSE
  cnvgd_due_to_update = PETSC_FALSE
  cnvgd_due_to_abs_update_p = PETSC_FALSE
  cnvgd_due_to_abs_update_s = PETSC_FALSE
  cnvgd_due_to_abs_update = PETSC_FALSE
  cnvgd_due_to_rel_update_p = PETSC_FALSE
  cnvgd_due_to_rel_update_s = PETSC_FALSE
  cnvgd_due_to_rel_update = PETSC_FALSE
  cnvgd_on_pressure = PETSC_FALSE
  cnvgd_on_saturation = PETSC_FALSE
  update_p = UNINITIALIZED_DOUBLE
  update_s = UNINITIALIZED_DOUBLE
  rsn_string = ''

  ! Update the residual
  flow_soln%residual = fixed_accum / pm_well%dt_flow
  call PMWellResidualFlow(pm_well)

  ! Update mass balance
  call PMWellMassBalance(pm_well)

  do k = 1,pm_well%well_grid%nsegments
    idof = flow_soln%ndof*(k-1)+1
    update_p(k) = flow_soln%update(idof)
    update_s(k) = flow_soln%update(idof+1)

    ! Absolute Solution Updates
    temp_real = dabs(update_p(k))
    if (temp_real < flow_soln%itol_abs_update_p) then
      cnvgd_due_to_abs_update_p(k) = PETSC_TRUE
    endif

    temp_real = dabs(update_s(k))
    if (temp_real > 0.d0) then
      if ((-1.d0*log10(temp_real)) >= &
           (-1.d0*log10(flow_soln%itol_abs_update_s))) then
        cnvgd_due_to_abs_update_s(k) = PETSC_TRUE
      endif
    else
      cnvgd_due_to_abs_update_s(k) = PETSC_TRUE
    endif

    ! Relative Solution Updates
    temp_real = dabs(update_p(k)/pm_well%well%pl(k))
    if (temp_real < flow_soln%itol_rel_update_p) then
      cnvgd_due_to_rel_update_p(k) = PETSC_TRUE
    endif
    temp_real = dabs(update_s(k)/pm_well%well%gas%s(k))
    if (temp_real < flow_soln%itol_rel_update_s) then
      cnvgd_due_to_rel_update_s(k) = PETSC_TRUE
    endif

    ! Liquid (water) Component
    if (dabs(fixed_accum(idof)) > zero_accumulation) then
      ! Absolute Residual
      temp_real = dabs(flow_soln%residual(idof))
      if (temp_real <= flow_soln%itol_abs_res) then
        cnvgd_due_to_abs_res(idof) = PETSC_TRUE
      endif

      ! Scaled Residual
      temp_real = dabs(flow_soln%residual(idof) / &
                       (fixed_accum(idof)/pm_well%dt_flow))
      if (temp_real <= flow_soln%itol_scaled_res) then
        cnvgd_due_to_scaled_res(idof) = PETSC_TRUE
      endif
    else
      cnvgd_due_to_abs_res(idof) = PETSC_TRUE
      cnvgd_due_to_scaled_res(idof) = PETSC_TRUE
    endif

    ! Gas (air) Component
    if (dabs(fixed_accum(idof+1)) > zero_accumulation) then
      ! Absolute Residual
      temp_real = dabs(flow_soln%residual(idof+1))
      if (temp_real <= flow_soln%itol_abs_res) then
        cnvgd_due_to_abs_res(idof+1) = PETSC_TRUE
      endif

      ! Scaled Residual
      temp_real = dabs(flow_soln%residual(idof+1) / &
                       (fixed_accum(idof+1)/pm_well%dt_flow))
      if (temp_real <= flow_soln%itol_scaled_res) then
        cnvgd_due_to_scaled_res(idof+1) = PETSC_TRUE
      endif
    else
      cnvgd_due_to_abs_res(idof+1) = PETSC_TRUE
      cnvgd_due_to_scaled_res(idof+1) = PETSC_TRUE
    endif
  enddo

  max_absolute_residual = maxval(dabs(flow_soln%residual))
  loc_max_abs_residual = maxloc(dabs(flow_soln%residual),1)

  max_scaled_residual = maxval(dabs(flow_soln%residual/ &
                                    (fixed_accum/pm_well%dt_flow)))
  loc_max_scaled_residual = maxloc(dabs(flow_soln%residual/ &
                                        (fixed_accum/pm_well%dt_flow)),1)

  max_absolute_update_p = maxval(dabs(update_p))
  loc_max_abs_update_p = maxloc(dabs(update_p),1)

  max_absolute_update_s = maxval(dabs(update_s))
  loc_max_abs_update_s = maxloc(dabs(update_s),1)

  max_relative_update_p = maxval(dabs(update_p/pm_well%well%pl))
  loc_max_rel_update_p = maxloc(dabs(update_p/pm_well%well%pl),1)

  max_relative_update_s = maxval(dabs(update_s/pm_well%well%gas%s))
  loc_max_rel_update_s = maxloc(dabs(update_s/pm_well%well%gas%s),1)

  do k = 1,pm_well%well_grid%nsegments*flow_soln%ndof
    if (cnvgd_due_to_scaled_res(k) .or. cnvgd_due_to_abs_res(k)) then
      cnvgd_due_to_residual(k) = PETSC_TRUE
    endif
  enddo
  do k = 1,pm_well%well_grid%nsegments
    if (cnvgd_due_to_abs_update_p(k) .or. &
        cnvgd_due_to_rel_update_p(k)) then
      cnvgd_on_pressure(k) = PETSC_TRUE
    endif
    if (cnvgd_due_to_abs_update_s(k) .or. &
        cnvgd_due_to_rel_update_s(k)) then
      cnvgd_on_saturation(k) = PETSC_TRUE
    endif
    if (cnvgd_on_pressure(k) .and. cnvgd_on_saturation(k)) then
      cnvgd_due_to_update(k) = PETSC_TRUE
    endif
  enddo

  if (all(cnvgd_due_to_abs_res)) then
    rsn_string = trim(rsn_string) // ' R '
  endif
  if (all(cnvgd_due_to_scaled_res)) then
    rsn_string = trim(rsn_string) // ' sR '
  endif
  if (all(cnvgd_due_to_abs_update)) then
    rsn_string = trim(rsn_string) // ' uP&uS '
  endif
  if (all(cnvgd_due_to_rel_update)) then
    rsn_string = trim(rsn_string) // ' ruP&ruS '
  endif

  if (all(cnvgd_due_to_residual) .and. all(cnvgd_due_to_update)) then
    flow_soln%converged = PETSC_TRUE
    flow_soln%not_converged = PETSC_FALSE
    pm_well%cumulative_dt_flow = pm_well%cumulative_dt_flow + pm_well%dt_flow
    pm_well%flow_soln%prev_soln%pl = pm_well%well%pl
    pm_well%flow_soln%prev_soln%sg = pm_well%well%gas%s
    !call PMWellUpdateSolutionFlow(pm_well)
  else
    flow_soln%converged = PETSC_FALSE
    flow_soln%not_converged = PETSC_TRUE
  endif

  if (pm_well%print_output) then
    write(out_string,'(i2," aR:",es10.2,"  sR:",es10.2,"  uP:", es10.2," &
          &  uS:",es10.2,"  ruP:",es10.2,"  ruS:",es10.2)') &
          n_iter,max_absolute_residual,max_scaled_residual, &
          max_absolute_update_p,max_absolute_update_s, &
          max_relative_update_p,max_relative_update_s
    call PrintMsg(pm_well%option,out_string)
    if (flow_soln%converged) then
      out_string = ' WELL FLOW Solution converged!  ---> ' // trim(rsn_string)
      call PrintMsg(pm_well%option,out_string)
    endif
  endif

end subroutine PMWellCheckConvergenceFlow

! ************************************************************************** !

subroutine PMWellCheckConvergenceTran(pm_well,n_iter,fixed_accum)
  !
  ! Checks transport solution convergence against prescribed tolerances.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  implicit none

  class(pm_well_type) :: pm_well
  PetscInt :: n_iter
  PetscReal :: fixed_accum(pm_well%tran_soln%ndof*pm_well%well_grid%nsegments)

  type(well_soln_tran_type), pointer :: soln
  character(len=MAXSTRINGLENGTH) :: out_string
  character(len=MAXSTRINGLENGTH) :: rsn_string
  PetscReal :: temp_real
  PetscBool :: cnvgd_due_to_residual(pm_well%well_grid%nsegments* &
                                     pm_well%tran_soln%ndof)
  PetscBool :: cnvgd_due_to_abs_res(pm_well%well_grid%nsegments* &
                                    pm_well%tran_soln%ndof)
  PetscBool :: cnvgd_due_to_scaled_res(pm_well%well_grid%nsegments* &
                                       pm_well%tran_soln%ndof)
  PetscBool :: cnvgd_due_to_update(pm_well%well_grid%nsegments*pm_well%tran_soln%ndof)
  PetscReal :: vol_vec(pm_well%well_grid%nsegments*pm_well%tran_soln%ndof)
  PetscReal :: aq_mass_vec(pm_well%well_grid%nsegments*pm_well%tran_soln%ndof)
  PetscReal :: max_scaled_residual,max_absolute_residual
  PetscReal :: max_update
  PetscInt :: loc_max_scaled_residual,loc_max_abs_residual
  PetscInt :: loc_max_update
  PetscInt :: k,n,j,S,TAG,last_rank
  PetscInt :: isegment, ispecies
  PetscErrorCode :: ierr

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  soln => pm_well%tran_soln

  n_iter = n_iter + 1
  soln%n_newton = soln%n_newton + 1
  S = pm_well%well_grid%nsegments*pm_well%tran_soln%ndof

  cnvgd_due_to_residual = PETSC_FALSE
  cnvgd_due_to_abs_res = PETSC_FALSE
  cnvgd_due_to_scaled_res = PETSC_FALSE
  cnvgd_due_to_update = PETSC_FALSE
  rsn_string = ''

  ! Update the residual
  soln%residual = 0.d0
  if (any(pm_well%option%myrank == pm_well%well_grid%h_rank_id)) then
    soln%residual = fixed_accum/pm_well%dt_tran 
    call PMWellResidualTran(pm_well)

    do k = 1,(pm_well%well_grid%nsegments*soln%ndof)
      ! Absolute Residual
      temp_real = dabs(soln%residual(k))
      if (temp_real < soln%itol_abs_res) then
        cnvgd_due_to_abs_res(k) = PETSC_TRUE
      endif
      ! Scaled Residual
      temp_real = dabs(soln%residual(k)/(fixed_accum(k)/pm_well%dt_tran))
      if (temp_real < soln%itol_scaled_res) then
        cnvgd_due_to_scaled_res(k) = PETSC_TRUE
      endif
    enddo

    ! Relative Update
    do n = 1,pm_well%well_grid%nsegments
      isegment = n
      do k = 1, soln%ndof
        ispecies = k
        j = ((isegment-1)*soln%ndof) + ispecies
        vol_vec(j) = pm_well%well%volume(isegment)
        aq_mass_vec(j) = pm_well%well%aqueous_mass(ispecies,isegment)
        temp_real = dabs(soln%update(j)*vol_vec(j)/aq_mass_vec(j))
        if (temp_real < soln%itol_rel_update) then
          cnvgd_due_to_update(j) = PETSC_TRUE
        endif
      enddo
    enddo
  
    max_absolute_residual = maxval(dabs(soln%residual))
    loc_max_abs_residual = maxloc(dabs(soln%residual),1)

    max_scaled_residual = maxval(dabs(soln%residual/ &
                                      (fixed_accum/pm_well%dt_tran)))
    loc_max_scaled_residual = maxloc(dabs(soln%residual/ &
                                          (fixed_accum/pm_well%dt_tran)),1)

    max_update = maxval(dabs(soln%update*vol_vec/aq_mass_vec))
    loc_max_update = maxloc(dabs(soln%update*vol_vec/aq_mass_vec),1)

    do k = 1,(pm_well%well_grid%nsegments*soln%ndof)
      if (cnvgd_due_to_scaled_res(k) .or. cnvgd_due_to_abs_res(k)) then
        cnvgd_due_to_residual(k) = PETSC_TRUE
      endif
    enddo
  endif

  call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
  last_rank = pm_well%well_comm%well_rank_list(pm_well%well_comm%commsize)
  if (pm_well%well_comm%commsize > 1) then
    TAG = 0
    if (pm_well%well_comm%rank == last_rank) then
      call MPI_Send(cnvgd_due_to_abs_res,S,MPI_LOGICAL,0,TAG, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call MPI_Send(cnvgd_due_to_scaled_res,S,MPI_LOGICAL,0,TAG+1, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call MPI_Send(cnvgd_due_to_update,S,MPI_LOGICAL,0,TAG+2, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
    endif
    if (pm_well%well_comm%rank == 0) then
      call MPI_Recv(cnvgd_due_to_abs_res,S,MPI_LOGICAL, &
                    last_rank,TAG,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
      call MPI_Recv(cnvgd_due_to_scaled_res,S,MPI_LOGICAL, &
                    last_rank,TAG+1,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
      call MPI_Recv(cnvgd_due_to_update,S,MPI_LOGICAL, &
                    last_rank,TAG+2,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
    endif
  endif

  call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
  if (all(cnvgd_due_to_abs_res)) then
    rsn_string = trim(rsn_string) // ' aR '
  endif
  if (all(cnvgd_due_to_scaled_res)) then
    rsn_string = trim(rsn_string) // ' sR '
  endif
  if (all(cnvgd_due_to_update)) then
    rsn_string = trim(rsn_string) // ' rU '
  endif

  if (pm_well%well_comm%commsize > 1) then
    TAG = 0
    if (pm_well%well_comm%rank == last_rank) then
      call MPI_Send(max_update,1,MPI_DOUBLE_PRECISION,0,TAG, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call MPI_Send(max_absolute_residual,1,MPI_DOUBLE_PRECISION,0,TAG+1, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call MPI_Send(max_scaled_residual,1,MPI_DOUBLE_PRECISION,0,TAG+2, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
    endif
    if (pm_well%well_comm%rank == 0) then
      call MPI_Recv(max_update,1,MPI_DOUBLE_PRECISION, &
                    last_rank,TAG,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
      call MPI_Recv(max_absolute_residual,1,MPI_DOUBLE_PRECISION, &
                    last_rank,TAG+1,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
      call MPI_Recv(max_scaled_residual,1,MPI_DOUBLE_PRECISION, &
                    last_rank,TAG+2,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
    endif
  endif

  call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
  write(out_string,'(i4,"    aR:",es10.3,"    sR:",es10.3,"    rU:", es10.3)')&
        n_iter,max_absolute_residual,max_scaled_residual, &
        max_update
  call PrintMsg(pm_well%option,out_string)

  if (pm_well%well_comm%commsize > 1) then
    TAG = 0
    if (pm_well%well_comm%rank == last_rank) then
      call MPI_Send(cnvgd_due_to_residual,S,MPI_LOGICAL,0,TAG, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
      call MPI_Send(cnvgd_due_to_update,S,MPI_LOGICAL,0,TAG+1, &
                    pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
    endif
    if (pm_well%well_comm%rank == 0) then
      call MPI_Recv(cnvgd_due_to_residual,S,MPI_LOGICAL, &
                    last_rank,TAG,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
      call MPI_Recv(cnvgd_due_to_update,S,MPI_LOGICAL, &
                    last_rank,TAG+1,pm_well%well_comm%comm,MPI_STATUS_IGNORE, &
                    ierr);CHKERRQ(ierr)
    endif
  endif

  call MPI_Barrier(pm_well%well_comm%comm,ierr);CHKERRQ(ierr)
  if (all(cnvgd_due_to_residual) .and. all(cnvgd_due_to_update)) then
    soln%converged = PETSC_TRUE
    soln%not_converged = PETSC_FALSE
    out_string = ' WELL TRAN Solution converged!  ---> ' // trim(rsn_string)
    call PrintMsg(pm_well%option,out_string)
    pm_well%cumulative_dt_tran = pm_well%cumulative_dt_tran + pm_well%dt_tran
  else
    soln%converged = PETSC_FALSE
    soln%not_converged = PETSC_TRUE
  endif

  pm_well%tran_soln%cut_ts_flag = PETSC_FALSE

end subroutine PMWellCheckConvergenceTran

! ************************************************************************** !

subroutine PMWellUpdateWellQ(well,reservoir)
  !
  ! Updates the src/sink vector for the fluid object.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021
  !

  use SCO2_Aux_module, only: sco2_fmw => fmw_comp

  implicit none

  type(well_type) :: well
  type(well_reservoir_type), pointer :: reservoir

  type(well_fluid_type), pointer :: liq
  type(well_fluid_type), pointer :: gas

  PetscReal, parameter :: threshold_p = 0.d0 !1.d-2 !1.d-1
  PetscReal :: mobility, den_ave
  PetscBool :: upwind
  PetscInt :: i, nsegments

  liq => well%liq
  gas => well%gas

  nsegments = size(well%liq%Q)

  ! + Q goes out of well to reservoir
  ! - Q goes into well from reservoir

  select case (well%well_model_type)
    !------------------------------------------------------------------------
    case(WELL_MODEL_WIPP_DARCY)
      do i = 1,nsegments
        if (dabs((reservoir%p_l(i)-well%pl(i)))/well%pl(i) > threshold_p) then
          upwind = reservoir%p_l(i) > well%pl(i)
          if (upwind) then
            mobility = reservoir%kr_l(i)/reservoir%visc_l(i)
          else
            mobility = liq%kr(i)/liq%visc(i)
          endif
          den_ave = 0.5d0 * (liq%den(i) + reservoir%den_l(i)) / FMWH2O
          ! Flowrate in kmol/s
          liq%Q(i) = den_ave*mobility*well%WI(i)* &
                     (reservoir%p_l(i)-well%pl(i))

          upwind = reservoir%p_g(i) > well%pg(i)
          if (upwind) then
            mobility = reservoir%kr_g(i)/reservoir%visc_g(i)
          else
            mobility = gas%kr(i)/gas%visc(i)
          endif
          den_ave = 0.5d0 * (gas%den(i) + reservoir%den_g(i)) / &
                             fmw_comp(TWO_INTEGER)

          ! Flowrate in kmol/s
          gas%Q(i) = den_ave*mobility*well%WI(i)* &
                     (reservoir%p_g(i)-well%pg(i))
        else
          liq%Q(i) = 0.d0
          gas%Q(i) = 0.d0
        endif
      enddo
    !------------------------------------------------------------------------
  end select

end subroutine PMWellUpdateWellQ

! ************************************************************************** !

subroutine PMWellComputeWellIndex(pm_well)
  !
  ! Computes the well index.
  !
  ! Author: Michael Nole
  ! Date: 12/22/2021
  !

  implicit none

  class(pm_well_type) :: pm_well
  type(well_reservoir_type), pointer :: reservoir

  PetscReal :: r0
  PetscReal, parameter :: PI=3.141592653589793d0
  PetscReal :: temp_real
  type(option_type), pointer :: option
  PetscInt :: k
  PetscReal :: kxy,kyx,kxz,kzx,kyz,kzy
  PetscReal :: dx,dy,dz,dh_x,dh_y,dh_z
  PetscReal :: r0x,r0y,r0z
  PetscReal :: wix,wiy,wiz
  PetscReal :: dx_tot, dy_tot, dz_tot
  character(len=8) :: diameter_string, dx_string

  if (pm_well%well_comm%comm == MPI_COMM_NULL) then
    pm_well%well%WI = UNINITIALIZED_DOUBLE
    return
  endif

  option => pm_well%option
  reservoir => pm_well%well%reservoir

  ! Peaceman Model: default = anisotropic
  ! This assumes z is vertical (not true for WIPP)
  select case(pm_well%well%WI_model)
    case(PEACEMAN_ISO)
      do k = 1,pm_well%well_grid%nsegments
        write(diameter_string,'(F7.4)') pm_well%well%diameter(k)
        write(dx_string,'(F7.4)') reservoir%dx(k)
        temp_real = log(2.079d-1*reservoir%dx(k)/ &
                        (pm_well%well%diameter(k)/2.d0))

        if (temp_real <= 0.d0) then
          option%io_buffer = 'Wellbore diameter (' // diameter_string // '&
          & m) is too large relative to reservoir dx (' // dx_string  //  '&
          & m). For the PEACEMAN_ISO model, wellbore diameter must be &
          &smaller than 0.4158 * reservoir dx.'
          call PrintErrMsg(option)
        endif

        pm_well%well%WI(k) = 2.d0*PI*reservoir%kx(k)*pm_well%well_grid%dh(k)/ &
                          temp_real
       enddo

    case(PEACEMAN_2D)
      do k = 1,pm_well%well_grid%nsegments
        write(diameter_string,'(F7.4)') pm_well%well%diameter(k)
        r0 = 2.8d-1*(sqrt(sqrt(reservoir%ky(k)/reservoir%kx(k))* &
             reservoir%dx(k)**2 + sqrt(reservoir%kx(k)/ &
             reservoir%ky(k))*reservoir%dy(k)**2) / &
             ((reservoir%ky(k)/reservoir%kx(k))**2.5d-1 + &
             (reservoir%kx(k)/reservoir%ky(k))**2.5d-1))

        temp_real = log(r0/(pm_well%well%diameter(k)/2.d0))

        if (temp_real <= 0.d0) then
          option%io_buffer = 'Wellbore diameter (' // diameter_string // ' m)&
          & is too large relative to reservoir discretization and &
          &permeability for the anisotropic PEACEMAN_2D well model.'
          call PrintErrMsg(option)
        endif

        pm_well%well%WI(k) = 2.d0*PI*sqrt(reservoir%kx(k)* &
                          reservoir%ky(k))*pm_well%well_grid%dh(k)/temp_real
      enddo
    case(PEACEMAN_3D)
      dx_tot = 0.d0
      dy_tot = 0.d0
      dz_tot = 0.d0
      do k = 1,pm_well%well_grid%nsegments
        write(diameter_string,'(F7.4)') pm_well%well%diameter(k)
        kyz = sqrt(reservoir%ky(k)/reservoir%kz(k))
        kzy = sqrt(reservoir%kz(k)/reservoir%ky(k))
        kzx = sqrt(reservoir%kz(k)/reservoir%kx(k))
        kxz = sqrt(reservoir%kx(k)/reservoir%kz(k))
        kyx = sqrt(reservoir%ky(k)/reservoir%kx(k))
        kxy = sqrt(reservoir%kx(k)/reservoir%ky(k))
        dx = reservoir%dx(k) / pm_well%well%friction_factor(k)
        dy = reservoir%dy(k) / pm_well%well%friction_factor(k)
        dz = reservoir%dz(k) / pm_well%well%friction_factor(k)
        r0x = 2.8d-1*sqrt(kyz*(dz**2) + kzy*(dy**2)) / &
              (sqrt(kyz)+sqrt(kzy))
        r0y = 2.8d-1*sqrt(kzx*(dx**2) + kxz*(dz**2)) / &
              (sqrt(kzx)+sqrt(kxz))
        r0z = 2.8d-1*sqrt(kyx*(dx**2) + kxy*(dy**2)) / &
              (sqrt(kyx)+sqrt(kxy))

        if (k == 1) then
          dh_x = pm_well%well_grid%h(k)%x - pm_well%well_grid%bottomhole(1)
          dh_y = pm_well%well_grid%h(k)%y - pm_well%well_grid%bottomhole(2)
          dh_z = pm_well%well_grid%h(k)%z - pm_well%well_grid%bottomhole(3)
        else
          dh_x = pm_well%well_grid%h(k)%x - dx_tot - &
                 pm_well%well_grid%bottomhole(1)
          dh_y = pm_well%well_grid%h(k)%y - dy_tot - &
                 pm_well%well_grid%bottomhole(2)
          dh_z = pm_well%well_grid%h(k)%z - dz_tot - &
                 pm_well%well_grid%bottomhole(3)
        endif
        dh_x = 2.d0 * dh_x
        dh_y = 2.d0 * dh_y
        dh_z = 2.d0 * dh_z
        dx_tot = dx_tot + dh_x
        dy_tot = dy_tot + dh_y
        dz_tot = dz_tot + dh_z

        wix = 2.d0 * PI * sqrt(reservoir%ky(k) * &
              reservoir%kz(k)) * dh_x / &
              (log(r0x/(pm_well%well%diameter(k)/2.d0)) + pm_well%well%skin(k))
        wiy = 2.d0 * PI * sqrt(reservoir%kx(k) * &
              reservoir%kz(k)) * dh_y / &
              (log(r0y/(pm_well%well%diameter(k)/2.d0)) + pm_well%well%skin(k))
        wiz = 2.d0 * PI * sqrt(reservoir%kx(k) * &
              reservoir%ky(k)) * dh_z / &
              (log(r0z/(pm_well%well%diameter(k)/2.d0)) + pm_well%well%skin(k))
        
        if (wix < 0.d0 .or. wiy < 0.d0 .or. wiz < 0.d0) then
          option%io_buffer = 'Wellbore diameter (' // diameter_string // ' m)&
          & is too large relative to reservoir discretization and &
          &permeability for the default anisotropic PEACEMAN_3D well model.'
          call PrintErrMsg(option)
        endif

        pm_well%well%WI(k) = sqrt((wix**2) + (wiy**2) + (wiz**2))
      enddo
      
  end select

  pm_well%well%WI = pm_well%well%WI*pm_well%well%WI_base

end subroutine PMWellComputeWellIndex

! ************************************************************************** !

subroutine PMWellAccumulationFlow(pm_well,well,id,Res)
  !
  ! Computes the accumulation term for the flow residual based on the
  ! chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 12/23/2021
  !

  use SCO2_Aux_module, only: sco2_fmw => fmw_comp

  implicit none

  type(pm_well_type) :: pm_well
  type(well_type) :: well
  PetscInt :: id
  PetscReal :: Res(pm_well%nphase)

  Res = 0.d0

  select case(well%well_model_type)
    !---------------------------------------------
    case(WELL_MODEL_WIPP_DARCY)
      ! liquid accumulation term
      Res(1) = Res(1) + well%liq%s(id) * well%liq%den(id) / FMWH2O * &
                well%phi(id) * well%volume(id) / pm_well%dt_flow
      ! gas accumulation term
      Res(2) = Res(2) + well%gas%s(id) * well%gas%den(id) / &
                fmw_comp(TWO_INTEGER) * well%phi(id) * well%volume(id) / &
                pm_well%dt_flow
    !---------------------------------------------
    case default
    !---------------------------------------------
  end select

end subroutine PMWellAccumulationFlow

! ************************************************************************** !

subroutine PMWellSrcSink(pm_well,well,id,Res)
  !
  ! Computes the source/sink term for the flow residual based on the
  ! chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 02/26/2022
  !

  implicit none

  type(pm_well_type) :: pm_well
  type(well_type) :: well
  PetscInt :: id
  PetscReal :: Res(pm_well%nphase)

  Res = 0.d0

  call PMWellUpdateWellQ(pm_well%well,pm_well%well%reservoir)

  select case(well%well_model_type)
    !---------------------------------------------
    case(WELL_MODEL_WIPP_DARCY)
      ! kmol/s
      Res(1) = Res(1) - well%liq%Q(id)
      ! kmol/s
      Res(2) = Res(2) - well%gas%Q(id)
    !---------------------------------------------
    case default
    !---------------------------------------------
  end select

end subroutine PMWellSrcSink

! ************************************************************************** !

subroutine PMWellAccumulationTran(pm_well,isegment,Res)
  !
  ! Computes the fixed accumulation term for the transport residual.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  implicit none

  type(pm_well_type) :: pm_well
  PetscInt :: isegment
  PetscReal :: Res(pm_well%nspecies)

  PetscInt :: ispecies, k

  Res = 0.d0

  ! porosity in [m^3-void/m^3-bulk]
  ! saturation in [m^3-liq/m^3-void]
  ! volume in [m^3-bulk]
  ! aqueous conc in [mol-species/m^3-liq]
  ! Res(:) in [mol-species]

  ! NOTE: division by dt occurs later, when the fixed accumulation is needed

  ! NOTE: this calculation uses the converged solution because it is using
  !       this%well%liq%s() and this%well%prev_soln%aqueous_conc() values

  do ispecies = 1,pm_well%nspecies
    k = ispecies
    Res(k) = pm_well%well%volume(isegment) * pm_well%well%phi(isegment) * &
             pm_well%well%liq%s(isegment) * &
             pm_well%tran_soln%prev_soln%aqueous_conc(ispecies,isegment)
  enddo

end subroutine PMWellAccumulationTran

! ************************************************************************** !

subroutine PMWellAccumDerivative(pm_well,local_id,Jac)
  !
  ! Computes the derivative of the accumulation term for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 01/05/2022
  !

  implicit none

  type(pm_well_type) :: pm_well
  PetscInt :: local_id
  PetscReal :: Jac(pm_well%nphase,pm_well%nphase)

  PetscInt :: idof, irow
  PetscReal :: res(pm_well%nphase),res_pert(pm_well%nphase)

  call PMWellAccumulationFlow(pm_well,pm_well%well,local_id,res)

  do idof = 1, pm_well%nphase
    call PMWellAccumulationFlow(pm_well,pm_well%well_pert(idof),local_id, &
                                res_pert)
    do irow = 1, pm_well%nphase
      Jac(irow,idof) = (res_pert(irow)-res(irow))/pm_well%pert(local_id,idof)
    enddo !irow
  enddo ! idof

end subroutine PMWellAccumDerivative
! ************************************************************************** !

subroutine PMWellSrcSinkDerivative(pm_well,local_id,Jac)
  !
  ! Computes the derivative of the source/sink term for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 02/26/2022
  !

  implicit none

  type(pm_well_type) :: pm_well
  PetscInt :: local_id
  PetscReal :: Jac(pm_well%nphase,pm_well%nphase)

  PetscInt :: idof, irow
  PetscReal :: res(pm_well%nphase),res_pert(pm_well%nphase)

  call PMWellSrcSink(pm_well,pm_well%well,local_id,res)

  do idof = 1, pm_well%nphase
    call PMWellSrcSink(pm_well,pm_well%well_pert(idof),local_id,res_pert)
    do irow = 1, pm_well%nphase
      Jac(irow,idof) = (res_pert(irow)-res(irow))/pm_well%pert(local_id,idof)
    enddo !irow
  enddo ! idof

end subroutine PMWellSrcSinkDerivative

! ************************************************************************** !

subroutine PMWellFlux(pm_well,well_up,well_dn,iup,idn,Res,save_flux)
  !
  ! Computes the internal flux terms for the residual based on
  ! the chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 12/23/2021
  !

  use SCO2_Aux_module, only: sco2_fmw => fmw_comp

  implicit none

  type(pm_well_type) :: pm_well
  type(well_type) :: well_up, well_dn
  PetscInt :: iup, idn
  PetscReal :: Res(pm_well%nphase)
  PetscBool :: save_flux

  type(well_grid_type), pointer :: well_grid

  PetscReal :: perm_den_mu_area_ave_over_dist(2), perm_den_mu_area_up(2), &
               perm_den_mu_area_dn(2)
  PetscReal :: perm_up, perm_dn, dist_up, dist_dn, density_kg_ave, rel_perm
  PetscReal :: gravity_term, delta_pressure
  PetscReal :: v_darcy
  PetscReal :: density_ave_kmol, tot_mole_flux
  PetscReal :: up_scale, dn_scale
  PetscBool :: upwind


  well_grid => pm_well%well_grid

  Res(:) = 0.d0

  select case(pm_well%well%well_model_type)
    case(WELL_MODEL_WIPP_DARCY)
      ! Vertical well, no Klinkenberg

        perm_up = well_up%permeability(iup)
        perm_dn = well_dn%permeability(idn)
        dist_up = well_grid%dh(iup)/2.d0
        dist_dn = well_grid%dh(idn)/2.d0

        perm_den_mu_area_up(1) = perm_up * well_up%liq%den(iup) / &
                              FMWH2O / well_up%liq%visc(iup) * &
                              PI * (well_up%diameter(iup)/2.d0)**2
        perm_den_mu_area_up(2) = perm_up * well_up%gas%den(iup) / &
                              fmw_comp(TWO_INTEGER) / well_up%gas%visc(iup) * &
                              PI * (well_up%diameter(iup)/2.d0)**2
        perm_den_mu_area_dn(1) = perm_dn * well_dn%liq%den(idn) / &
                              FMWH2O / well_dn%liq%visc(idn) * &
                              PI * (well_dn%diameter(idn)/2.d0)**2
        perm_den_mu_area_dn(2) = perm_dn * well_dn%gas%den(idn) / &
                              fmw_comp(TWO_INTEGER) / well_dn%gas%visc(idn) * &
                              PI * (well_dn%diameter(idn)/2.d0)**2

        perm_den_mu_area_ave_over_dist(1) = &
               (perm_den_mu_area_up(1) * perm_den_mu_area_dn(1)) / &
               (dist_up*perm_den_mu_area_dn(1) + dist_dn*perm_den_mu_area_up(1))

        perm_den_mu_area_ave_over_dist(2) = &
               (perm_den_mu_area_up(2) * perm_den_mu_area_dn(2)) / &
               (dist_up*perm_den_mu_area_dn(2) + dist_dn*perm_den_mu_area_up(2))

        ! Liquid flux
        density_kg_ave = 0.5d0*(well_up%liq%den(iup)+well_dn%liq%den(idn))
        ! Assuming the well is always vertical and gravity is in the
        ! (-) direction
        gravity_term = density_kg_ave * gravity * &
                       0.5d0*(well_grid%dh(iup)+well_grid%dh(idn))
        delta_pressure = well_up%pl(iup) - well_dn%pl(idn) + &
                         gravity_term
        up_scale = 0.d0
        dn_scale = 0.d0

        upwind = delta_pressure > 0.d0

        ! Only upwinding the perm here, not mobility ratio
        ! following the LumpedHarmonic BRAGFLO formulation
        if (upwind) then
          up_scale = 1.d0
          rel_perm = well_up%liq%kr(iup)
        else
          dn_scale = 1.d0
          rel_perm = well_dn%liq%kr(idn)
        endif

        !kmol/sec
        tot_mole_flux = perm_den_mu_area_ave_over_dist(1) * rel_perm * &
                  delta_pressure
        density_ave_kmol = density_kg_ave / fmw_comp(ONE_INTEGER)
        ! v_darcy = kmol/sec / kmol/m^3 / area[m^2]
        v_darcy = tot_mole_flux/density_ave_kmol/(5.d-1*(well_up%area(iup)+ &
                  well_dn%area(idn)))
        ! Store flux calculation for consistency with transport
        if (save_flux) then
          well_up%ql(iup) = v_darcy
          well_up%ql_kmol(iup) = tot_mole_flux
        endif

        Res(1) = Res(1) + tot_mole_flux

        ! Gas flux
        density_kg_ave = 0.5d0*(well_up%gas%den(iup)+well_dn%gas%den(idn))
        ! Assuming the well is always vertical and gravity is in the
        ! (-) direction
        gravity_term = density_kg_ave * gravity * &
                       0.5d0*(well_grid%dh(iup)+well_grid%dh(idn))
        delta_pressure = well_up%pg(iup) - well_dn%pg(idn) + &
                         gravity_term

        up_scale = 0.d0
        dn_scale = 0.d0

        upwind = delta_pressure > 0.d0

        ! Only upwinding the perm here, not mobility ratio
        ! following the LumpedHarmonic BRAGFLO formulation
        if (upwind) then
          up_scale = 1.d0
          rel_perm = well_up%gas%kr(iup)
        else
          dn_scale = 1.d0
          rel_perm = well_dn%gas%kr(idn)
        endif

        !kmol/sec
        tot_mole_flux = perm_den_mu_area_ave_over_dist(2) * rel_perm * &
                        delta_pressure
        density_ave_kmol = density_kg_ave / fmw_comp(TWO_INTEGER)
        ! v_darcy [m/sec] = mole flux [kmol/sec] / den [kmol/m^3] / area[m^2]
        v_darcy = tot_mole_flux/density_ave_kmol/(5.d-1*(well_up%area(iup) + &
                            well_dn%area(idn)))
        ! Store flux calculation for consistency with transport
        if (save_flux) then
          well_up%qg(iup) = v_darcy
          well_up%qg_kmol(iup) = tot_mole_flux
        endif

        Res(2) = Res(2) + tot_mole_flux

    case default

  end select

end subroutine PMWellFlux

! ************************************************************************** !

subroutine PMWellFluxDerivative(pm_well,iup,idn,Jup,Jdn)
  !
  ! Computes the derivative of the internal flux terms for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 01/05/2022
  !

  implicit none

  type(pm_well_type) :: pm_well
  PetscInt :: iup, idn, iphase, irow
  PetscReal :: Jup(pm_well%nphase,pm_well%nphase), &
               Jdn(pm_well%nphase,pm_well%nphase)

  PetscReal :: res_up(pm_well%nphase),res_dn(pm_well%nphase), &
               res_pert(pm_well%nphase)

  call PMWellFlux(pm_well,pm_well%well,pm_well%well,iup,idn,res_up,PETSC_FALSE)

  res_dn = res_up

  ! upgradient derivatives
  do iphase = 1,pm_well%nphase
    call PMWellFlux(pm_well,pm_well%well_pert(iphase),pm_well%well,iup,idn, &
                    res_pert,PETSC_FALSE)
    do irow = 1, pm_well%nphase
      Jup(irow,iphase) = (res_pert(irow)-res_up(irow)) / &
                         pm_well%pert(iup,iphase)
    enddo !irow
  enddo

  ! downgradient derivatives
  do iphase = 1,pm_well%nphase
    call PMWellFlux(pm_well,pm_well%well,pm_well%well_pert(iphase),iup,idn, &
                    res_pert,PETSC_FALSE)
    do irow = 1, pm_well%nphase
      Jdn(irow,iphase) = (res_pert(irow)-res_dn(irow)) / &
                         pm_well%pert(idn,iphase)
    enddo !irow
  enddo

end subroutine PMWellFluxDerivative

! ************************************************************************** !

subroutine PMWellBCFlux(pm_well,well,Res,save_flux)
  !
  ! Computes the boundary flux terms for the residual based on the
  ! chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 12/24/2021
  !

  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_WIPP_module
  use SCO2_Aux_module, only: sco2_fmw => fmw_comp

  implicit none

  type(pm_well_type) :: pm_well
  type(well_type) :: well
  PetscReal :: Res(2*pm_well%nphase)
  PetscBool :: save_flux

  type(option_type), pointer :: option
  type(well_grid_type), pointer :: well_grid
  type(well_reservoir_type), pointer :: reservoir
  class(characteristic_curves_type), pointer :: characteristic_curves
  class(sat_func_base_type), pointer :: saturation_function

  !MAN: clean these up
  PetscReal :: perm_ave_over_dist
  PetscReal :: gravity_term, delta_pressure
  PetscReal :: density_ave, tot_mole_flux
  PetscReal :: boundary_pressure, boundary_den
  PetscReal :: boundary_pg, boundary_krg, dn_scale
  PetscReal :: t,dwmol,dwp,dwt,Psat,visl,visg
  PetscReal :: Pc,dpc_dsatl,dkrl_dsatl,dkrg_dsatl
  PetscReal :: v_darcy,q,rel_perm
  PetscBool :: upwind
  PetscInt :: itop
  PetscErrorCode :: ierr

  option => pm_well%option

  well_grid => pm_well%well_grid
  reservoir => pm_well%well%reservoir

  t = 25.d0 !Constant temperature

  Res(:) = 0.d0

  select case(well%well_model_type)
    case(WELL_MODEL_WIPP_DARCY)
      ! Vertical well, no Klinkenberg.
      itop = pm_well%well_grid%nsegments
      if (pm_well%flow_soln%bh_p) then
        !Dirichlet pressure and saturation at the bottom

        characteristic_curves => pm_well%realization%patch% &
                                characteristic_curves_array(well%ccid(1))%ptr
        saturation_function => characteristic_curves%saturation_function

        ! Water Residual

        perm_ave_over_dist = well%permeability(1) / (well_grid%dh(1)/2.d0)
        boundary_pressure = well%bh_p
        gravity_term = well%liq%den(1) * gravity * &
                       well_grid%dh(1)/2.d0
        delta_pressure = boundary_pressure - well%pl(1) + gravity_term

        call EOSWaterSaturationPressure(t,Psat,ierr)
        call EOSWaterDensityBRAGFLO(t,boundary_pressure,PETSC_FALSE, &
                                boundary_den,dwmol,dwp,dwt,ierr)

        upwind = delta_pressure > 0.d0
        if (upwind) then
          call characteristic_curves%liq_rel_perm_function% &
               RelativePermeability(1.d0-well%bh_sg,rel_perm,dkrl_dsatl,option)
          call EOSWaterViscosity(t,boundary_pressure,Psat,visl,ierr)
        else
          dn_scale = 1.d0
          rel_perm = well%liq%kr(1)
          visl = well%liq%visc(1)
        endif

        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy = perm_ave_over_dist * rel_perm / visl * &
                  delta_pressure
        if (upwind) then
          density_ave = boundary_den
        else
          density_ave = (well%liq%den(1)+boundary_den) / &
                        (2.d0 * fmw_comp(ONE_INTEGER))
        endif
        q = v_darcy * well%area(1)
        tot_mole_flux = q * density_ave
        ! Store boundary flux for consistency with transport
        if (save_flux) then
          well%ql_bc(1) = v_darcy
          well%ql_kmol_bc(1) = tot_mole_flux
        endif
        Res(1) = Res(1) + tot_mole_flux

        ! Gas Residual

        ! Capillary Pressure
        select type(sat_func => saturation_function)
          class is (sat_func_KRP3_type)
            if (.not. option%flow%pct_updated) then
              sat_func%pct = sat_func%pct_a * well% permeability(1) ** &
                             sat_func%pct_exp
              option%flow%pct_updated = PETSC_TRUE
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            else
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            endif
          class is (sat_func_KRP4_type)
            if (.not. option%flow%pct_updated) then
              sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
              option%flow%pct_updated = PETSC_TRUE
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            else
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            endif
          class is (sat_func_KRP5_type)
            if (.not. option%flow%pct_updated) then
              sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
              option%flow%pct_updated = PETSC_TRUE
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            else
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            endif
          class default
            call sat_func% &
               CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
        end select
        boundary_pg = boundary_pressure + Pc

        call characteristic_curves%gas_rel_perm_function% &
               RelativePermeability(1.d0-well%bh_sg,boundary_krg, &
               dkrg_dsatl,option)

        gravity_term = well%gas%den(1) * gravity * &
                       well_grid%dh(1)/2.d0
        delta_pressure = boundary_pg - well%pg(1) + gravity_term

        call EOSGasDensity(t,boundary_pg,boundary_den,ierr)

        upwind = delta_pressure > 0.d0
        if (upwind) then
          call characteristic_curves%gas_rel_perm_function% &
               RelativePermeability(1.d0-well%bh_sg,rel_perm,dkrl_dsatl,option)
          call EOSGasViscosity(t,boundary_pg,boundary_pg,boundary_den,visg,ierr)
        else
          dn_scale = 1.d0
          rel_perm = well%gas%kr(1)
          visg = well%gas%visc(1)
        endif

        v_darcy = perm_ave_over_dist * rel_perm/visg * &
                  delta_pressure
        if (upwind) then
          density_ave = boundary_den
        else
          density_ave = (well%gas%den(1)+boundary_den) / &
                        (2.d0 *fmw_comp(TWO_INTEGER))
        endif
        q = v_darcy * well%area(1)
        tot_mole_flux = q * density_ave
        ! Store boundary flux for consistency with transport
        if (save_flux) then
          well%qg_bc(1) = v_darcy
          well%qg_kmol_bc(1) = tot_mole_flux
        endif
        Res(2) = Res(2) + tot_mole_flux

      else if (pm_well%flow_soln%bh_q) then
        !Neumann flux at the bottom
        v_darcy = well%bh_ql

        if (v_darcy > 0.d0) then
          density_ave = reservoir%den_l(1) / fmw_comp(ONE_INTEGER)
        else
          density_ave = well%liq%den(1) / fmw_comp(ONE_INTEGER)
        endif
        q = v_darcy * well%area(1)
        tot_mole_flux = q * density_ave
        ! Store boundary flux for consistency with transport
        if (save_flux) then
          well%ql_bc(1) = v_darcy
          well%ql_kmol_bc(1) = tot_mole_flux
        endif
        Res(1) = Res(1) + tot_mole_flux

        v_darcy = well%bh_qg

        if (v_darcy > 0.d0) then
          density_ave = reservoir%den_g(1) / fmw_comp(TWO_INTEGER)
        else
          density_ave = well%gas%den(1) / fmw_comp(TWO_INTEGER)
        endif
        q = v_darcy * well%area(1)
        tot_mole_flux = q * density_ave
        ! Store boundary flux for consistency with transport
        if (save_flux) then
          well%qg_bc(1) = v_darcy
          well%qg_kmol_bc(1) = tot_mole_flux
        endif
        Res(2) = Res(2) + tot_mole_flux
      else
        ! this should not happen once error messaging is updated
      endif

      if (pm_well%flow_soln%th_p) then
        !Dirichlet pressure and saturation at the top
        characteristic_curves => pm_well%realization%patch% &
                                characteristic_curves_array(well%ccid(itop))%ptr
        saturation_function => characteristic_curves%saturation_function

        ! Water Residual

        perm_ave_over_dist = well%permeability(itop) / (well_grid%dh(itop)/2.d0)
        boundary_pressure = well%th_p
        gravity_term = well%liq%den(itop) * gravity * &
                       well_grid%dh(itop)/2.d0
        delta_pressure = well%pl(itop) - boundary_pressure + gravity_term

        call EOSWaterSaturationPressure(t,Psat,ierr)
        call EOSWaterDensityBRAGFLO(t,boundary_pressure,PETSC_FALSE, &
                                boundary_den,dwmol,dwp,dwt,ierr)

        upwind = delta_pressure < 0.d0
        if (upwind) then
          call characteristic_curves%liq_rel_perm_function% &
               RelativePermeability(1.d0-well%th_sg,rel_perm,dkrl_dsatl,option)
          call EOSWaterViscosity(t,boundary_pressure,Psat,visl,ierr)
        else
          dn_scale = 1.d0
          rel_perm = well%liq%kr(itop)
          visl = well%liq%visc(itop)
        endif

        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy = perm_ave_over_dist * rel_perm / visl * &
                  delta_pressure

        density_ave = (well%liq%den(itop)+boundary_den) / &
                      (2.d0 * fmw_comp(ONE_INTEGER))
        q = v_darcy * well%area(itop)
        tot_mole_flux = q * density_ave
        ! Store boundary flux for consistency with transport
        if (save_flux) then
          well%ql_bc(2) = v_darcy
          well%ql_kmol_bc(2) = tot_mole_flux
        endif
        Res(3) = Res(3) - tot_mole_flux

        ! Gas Residual

        ! Capillary Pressure
        select type(sat_func => saturation_function)
          class is (sat_func_KRP3_type)
            if (.not. option%flow%pct_updated) then
              sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
              option%flow%pct_updated = PETSC_TRUE
              call sat_func% &
                   CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
            else
              call sat_func% &
                   CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
            endif
          class is (sat_func_KRP4_type)
            if (.not. option%flow%pct_updated) then
              sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
              option%flow%pct_updated = PETSC_TRUE
              call sat_func% &
                   CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
            else
              call sat_func% &
                   CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
            endif
          class is (sat_func_KRP5_type)
            if (.not. option%flow%pct_updated) then
              sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
              option%flow%pct_updated = PETSC_TRUE
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            else
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            endif
          class default
            call sat_func% &
               CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
        end select
        boundary_pg = boundary_pressure + Pc

        call characteristic_curves%gas_rel_perm_function% &
               RelativePermeability(1.d0-well%th_sg,boundary_krg, &
               dkrg_dsatl,option)

        gravity_term = well%gas%den(itop) * gravity * &
                       well_grid%dh(itop)/2.d0
        delta_pressure = well%pg(itop) - boundary_pg + gravity_term

        call EOSGasDensity(t,boundary_pg,boundary_den,ierr)

        upwind = delta_pressure < 0.d0
        if (upwind) then
          call characteristic_curves%gas_rel_perm_function% &
               RelativePermeability(1.d0-well%th_sg,rel_perm,dkrl_dsatl,option)
          call EOSGasViscosity(t,boundary_pg,boundary_pg,boundary_den,visg,ierr)
        else
          dn_scale = 1.d0
          rel_perm = well%gas%kr(itop)
          visg = well%gas%visc(itop)
        endif

        v_darcy = perm_ave_over_dist * rel_perm/visg * &
                  delta_pressure

        density_ave = (well%gas%den(itop)+boundary_den) / &
                      (2.d0 *fmw_comp(TWO_INTEGER))
        q = v_darcy * well%area(itop)
        tot_mole_flux = q * density_ave
        ! Store boundary flux for consistency with transport
        if (save_flux) then
          well%qg_bc(2) = v_darcy
          well%qg_kmol_bc(2) = tot_mole_flux
        endif
        Res(4) = Res(4) - tot_mole_flux
      else
        !Neumann flux at the top
        v_darcy = -well%th_ql

        ! Always take well density with tophole flux bc
        density_ave = well%liq%den(itop) / fmw_comp(ONE_INTEGER)
        q = v_darcy * well%area(itop)
        tot_mole_flux = q * density_ave
        ! Store boundary flux for consistency with transport
        if (save_flux) then
          well%ql_bc(2) = v_darcy
          well%ql_kmol_bc(2) = tot_mole_flux
        endif
        Res(3) = Res(3) - tot_mole_flux

        v_darcy = -well%th_qg

        ! Always take well density with tophole flux bc
        density_ave = well%gas%den(itop) / fmw_comp(TWO_INTEGER)
        q = v_darcy * well%area(itop)
        tot_mole_flux = q * density_ave
        ! Store boundary flux for consistency with transport
        if (save_flux) then
          well%qg_bc(2) = v_darcy
          well%qg_kmol_bc(2) = tot_mole_flux
        endif
        Res(4) = Res(4) - tot_mole_flux
      endif
    case default

  end select

end subroutine PMWellBCFlux

! ************************************************************************** !

subroutine PMWellBCFluxDerivative(pm_well,Jtop,Jbtm)
  !
  ! Computes the derivative of the boundary flux terms for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 01/05/2022
  !

  implicit none

  type(pm_well_type) :: pm_well
  PetscReal :: Jtop(pm_well%flow_soln%ndof,pm_well%flow_soln%ndof), &
               Jbtm(pm_well%flow_soln%ndof,pm_well%flow_soln%ndof)

  PetscInt :: idof, irow
  PetscReal :: res(2*pm_well%flow_soln%ndof),res_pert(2*pm_well%flow_soln%ndof)

  call PMWellBCFlux(pm_well,pm_well%well,res,PETSC_FALSE)

  Jtop = 0.d0
  Jbtm = 0.d0

  ! downgradient derivatives
  do idof = 1,pm_well%nphase
    call PMWellBCFlux(pm_well,pm_well%well_pert(idof),res_pert,PETSC_FALSE)
    do irow = 1, pm_well%nphase
      Jbtm(irow,idof) = (res_pert(irow)-res(irow)) / &
                        pm_well%pert(1,idof)
    enddo
    do irow = 1, pm_well%nphase
      Jtop(irow,idof) = (res_pert(irow + pm_well%flow_soln%ndof)- &
                         res(irow + pm_well%flow_soln%ndof)) / &
                         pm_well%pert(pm_well%well_grid%nsegments,idof)
    enddo
  enddo


end subroutine PMWellBCFluxDerivative

! ************************************************************************** !

subroutine PMWellPerturb(pm_well)
  !
  ! Calculates the state variables for the perturbed well system.
  !
  ! Author: Michael Nole
  ! Date: 01/06/2022
  !

  implicit none

  type(pm_well_type) :: pm_well

  PetscReal :: x(pm_well%well_grid%nsegments,pm_well%nphase), &
               pert(pm_well%well_grid%nsegments,pm_well%nphase)
  PetscInt :: i

  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-10

  x(:,ONE_INTEGER) = pm_well%well%pl
  x(:,TWO_INTEGER) = pm_well%well%gas%s

  ! Re-initialize perturbations
  pm_well%well_pert(ONE_INTEGER)%pl = pm_well%well%pl
  pm_well%well_pert(TWO_INTEGER)%pl = pm_well%well%pl
  pm_well%well_pert(ONE_INTEGER)%gas%s = pm_well%well%gas%s
  pm_well%well_pert(TWO_INTEGER)%gas%s = pm_well%well%gas%s

  pert(:,ONE_INTEGER) = perturbation_tolerance*x(:,ONE_INTEGER) + &
                        min_perturbation
  do i = 1,pm_well%well_grid%nsegments
    if (x(i,TWO_INTEGER) > 0.5d0) then
      pert(i,TWO_INTEGER) = -1.d0 * perturbation_tolerance
    else
      pert(i,TWO_INTEGER) = perturbation_tolerance
    endif
  enddo

  pm_well%well_pert(ONE_INTEGER)%pl = x(:,ONE_INTEGER) + pert(:,ONE_INTEGER)
  pm_well%well_pert(TWO_INTEGER)%gas%s = x(:,TWO_INTEGER) + pert(:,TWO_INTEGER)

  ! Update perturbed well properties
  select case (pm_well%option%iflowmode)
    case(WF_MODE)
      call PMWellUpdatePropertiesWIPPFlow(pm_well, &
                        pm_well%well_pert(ONE_INTEGER), &
                        pm_well%realization%patch%characteristic_curves_array,&
                        pm_well%realization%option)
      call PMWellUpdatePropertiesWIPPFlow(pm_well, &
                        pm_well%well_pert(TWO_INTEGER), &
                        pm_well%realization%patch%characteristic_curves_array,&
                        pm_well%realization%option)
    case(SCO2_MODE)
      call PMWellUpdatePropertiesSCO2Flow(pm_well, &
                        pm_well%well_pert(ONE_INTEGER), &
                        pm_well%realization%option)
      call PMWellUpdatePropertiesSCO2Flow(pm_well, &
                        pm_well%well_pert(TWO_INTEGER), &
                        pm_well%realization%option)
  end select
  ! Update perturbed source/sink term from the reservoir
  call PMWellUpdateWellQ(pm_well%well_pert(ONE_INTEGER), &
                         pm_well%well_pert(ONE_INTEGER)%reservoir)
  call PMWellUpdateWellQ(pm_well%well_pert(TWO_INTEGER), &
                         pm_well%well_pert(TWO_INTEGER)%reservoir)

  pm_well%pert = pert

end subroutine PMWellPerturb

! ************************************************************************** !

subroutine PMWellUpdatePropertiesWIPPFlow(pm_well,well, &
                                          characteristic_curves_array, &
                                          option)
  !
  ! Updates flow well object properties, when WIPP_FLOW is the flow mode.
  !
  ! Author: Michael Nole
  ! Date: 01/06/2022
  !

  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_WIPP_module

  implicit none

  class(pm_well_type) :: pm_well
  type(well_type) :: well
  type(characteristic_curves_ptr_type), pointer::characteristic_curves_array(:)
  type(option_type) :: option

  class(characteristic_curves_type), pointer :: characteristic_curves
  class(sat_func_base_type), pointer :: saturation_function
  type(strata_type), pointer :: strata
  PetscInt :: i,nsegments
  PetscReal :: T,dw,dg,dwmol,dwp,dwt,Psat,visl,visg
  PetscReal :: Pc,dpc_dsatl,krl,dkrl_dsatl,krg,dkrg_dsatl
  PetscErrorCode :: ierr

  T = option%flow%reference_temperature

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  nsegments =pm_well%well_grid%nsegments

  do i = 1,nsegments
    ! Material Properties
    strata => pm_well%strata_list%first
    do
      if (.not.associated(strata)) exit
      if (strata%id == pm_well%well_grid%strata_id(i)) then
        well%ccid(i) = strata%material_property%saturation_function_id
        well%permeability(i) = strata%material_property%permeability(3,3)
        well%phi(i) = strata%material_property%porosity
        exit
      endif
      strata => strata%next
    enddo

    ! Saturations
    well%gas%s(i) = max(well%gas%s(i),0.d0)
    well%gas%s(i) = min(well%gas%s(i),1.d0)
    well%liq%s(i) = 1.d0 - well%gas%s(i)

    ! Capillary Pressure
    characteristic_curves => characteristic_curves_array(well%ccid(i))%ptr
    saturation_function => characteristic_curves%saturation_function
    select type(sat_func => saturation_function)
      class is (sat_func_KRP3_type)
        if (.not. option%flow%pct_updated) then
          sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
                   CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
        else
          call sat_func% &
                   CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
        endif
      class is (sat_func_KRP4_type)
        if (.not. option%flow%pct_updated) then
          sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
               CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
        else
          call sat_func% &
               CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
        endif
      class is (sat_func_KRP5_type)
            if (.not. option%flow%pct_updated) then
              sat_func%pct = sat_func%pct_a * well%permeability(1) ** &
                             sat_func%pct_exp
              option%flow%pct_updated = PETSC_TRUE
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            else
              call sat_func% &
                   CapillaryPressure(1.d0-well%bh_sg,Pc,dpc_dsatl,option)
            endif
      class default
        call sat_func% &
               CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
    end select
    well%pg(i) = well%pl(i) + Pc

    ! Relative Permeabilities
    call characteristic_curves%liq_rel_perm_function% &
               RelativePermeability(well%liq%s(i),krl,dkrl_dsatl,option)
    well%liq%kr(i) = krl

    call characteristic_curves%gas_rel_perm_function% &
               RelativePermeability(well%liq%s(i),krg,dkrg_dsatl,option)
    well%gas%kr(i) = krg

    !Density
    call EOSWaterDensityBRAGFLO(T,well%pl(i),PETSC_FALSE, &
                                dw,dwmol,dwp,dwt,ierr)
    call EOSGasDensity(T,well%pg(i),dg,ierr)

    well%liq%den(i) = dw
    !No water vapor in WIPP_Darcy mode
    well%gas%den(i) = dg

    !Viscosity
    call EOSWaterSaturationPressure(T,Psat,ierr)
    call EOSWaterViscosity(T,well%pl(i),Psat,visl,ierr)
    call EOSGasViscosity(T,well%pg(i),well%pg(i),dg,visg,ierr)

    well%liq%visc(i) = visl
    well%gas%visc(i) = visg

  enddo


end subroutine PMWellUpdatePropertiesWIPPFlow

! ************************************************************************** !
subroutine PMWellUpdatePropertiesSCO2Flow(pm_well,well,option)
  !
  ! Updates flow well object properties, when SCO2 mode is the flow mode.
  !
  ! Author: Michael Nole
  ! Date: 01/06/2022
  !

  use EOS_Water_module
  use EOS_Gas_module
  use SCO2_Aux_module

  implicit none

  class(pm_well_type) :: pm_well
  type(well_type) :: well
  type(option_type) :: option

  PetscInt :: i,nsegments
  PetscReal :: drho_dT,drho_dP
  PetscReal :: xsl, Pco2, Pvap, Pva, Ps, Prvap 
  PetscReal :: den_kg_water, den_kg_steam, &
               den_kg_brine, den_kg_liq, &
               den_kg_gas, den_mol_co2, &
               den_kg_co2
  PetscReal :: visc_co2, visc_water, visc_liq, &
               visc_brine, visc_gas
  PetscReal :: xco2g, xwg, xco2l,xwl, &
               xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl
  PetscInt :: wid, co2_id, sid
  PetscErrorCode :: ierr

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id
  
  nsegments =pm_well%well_grid%nsegments

  if (well%th_qg > 0.d0) then
    ! CO2 Injection well. Need to update to flexibly accommodate humidity.
    well%gas%xmass(:,:) = 0.d0
    well%gas%xmass(:,TWO_INTEGER) = 1.d0
  elseif (well%th_ql > 0.d0) then
    ! Liquid Injection well: Need to update to flexibly 
    ! accommodate dissolved gas.
    well%liq%xmass(:,:) = 0.d0
    well%liq%xmass(:,ONE_INTEGER) = 1.d0
  endif

  do i = 1,nsegments

    !Liquid Density
    xsl = well%liq%xmass(i,option%salt_id)
    call SCO2BrineSaturationPressure(well%temp(i), &
                                     xsl,Ps)
    call SCO2BrineDensity(well%temp(i), well%pg(i), &
                          xsl, den_kg_brine, option)
    call SCO2VaporPressureBrine(well%temp(i), Ps, &
                                0.d0, den_kg_brine, &
                                xsl, Prvap)
    call SCO2WaterDensity(well%temp(i),Prvap, &
                          TWO_INTEGER,den_kg_water, &
                          den_kg_steam,option)
    call SCO2DensityCompositeLiquid(well%temp(i),den_kg_brine, &
                                  well%liq%xmass(i,option%co2_id), &
                                  den_kg_liq)

    well%liq%den(i) = den_kg_liq

    call SCO2Equilibrate(well%temp(i),well%pg(i), &
                         Pco2, Pvap, Ps, Prvap, &
                         xco2g, xwg, xco2l, xsl, xwl, &
                         xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl, option)

    xmolco2g = (well%gas%xmass(i,co2_id)/fmw_comp(2)) / &
               ((well%gas%xmass(i,co2_id)/fmw_comp(2)) + &
                (well%gas%xmass(i,wid)/fmw_comp(1)))
    xmolwg = 1.d0 - xmolco2g
    xmolco2l = (well%liq%xmass(i,co2_id)/fmw_comp(2)) / &
               ((well%gas%xmass(i,co2_id)/fmw_comp(2)) + &
               (well%gas%xmass(i,wid)/fmw_comp(1)) + &
               (well%gas%xmass(i,sid)/fmw_comp(3)))
    xmolwl = (well%liq%xmass(i,wid)/fmw_comp(2)) / &
               ((well%gas%xmass(i,co2_id)/fmw_comp(2)) + &
               (well%gas%xmass(i,wid)/fmw_comp(1)) + &
               (well%gas%xmass(i,sid)/fmw_comp(3)))
    xmolsl = 1.d0 - xmolco2l - xmolwl

    !Gas Density
    Pva = max(well%pg(i),Prvap)
    call EOSGasDensity(well%temp(i),Pva, &
                       den_mol_co2,drho_dT,drho_dP,ierr)
    den_kg_co2 = den_mol_co2 * fmw_comp(2)
    den_kg_gas = well%gas%xmass(i,option%co2_id) * &
                 den_kg_co2 + &
                 well%gas%xmass(i,option%water_id) * &
                 den_kg_steam
    well%gas%den(i) = den_kg_gas

    ! Liquid Viscosity
    call SCO2ViscosityWater(well%temp(i),well%pg(i), &
                           den_kg_water,visc_water,option)
    call SCO2ViscosityCO2(well%temp(i), den_kg_co2, &
                          visc_co2)
    call SCO2ViscosityBrine(well%temp(i), xsl, &
                           visc_water, visc_brine)
    call SCO2ViscosityLiquid(xmolco2l, visc_brine, &
                             visc_co2, visc_liq)

    well%liq%visc(i) = visc_liq

    ! Gas Viscosity
    call SCO2ViscosityGas(visc_water,visc_co2,xmolwg, &
                          xmolco2g,visc_gas)
    
    well%gas%visc(i) = visc_gas

  enddo


end subroutine PMWellUpdatePropertiesSCO2Flow

! ************************************************************************** !

subroutine PMWellUpdatePropertiesTran(pm_well)
  !
  ! Updates transport related well object properties for each time step.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/13/2022
  !

  use Transport_Constraint_Base_module
  use Transport_Constraint_NWT_module
  use Transport_Constraint_module
  use Condition_module

  implicit none

  class(pm_well_type) :: pm_well

  type(tran_condition_type), pointer :: tran_condition
  class(tran_constraint_base_type), pointer :: cur_constraint

  ! update the top of hole boundary condition with current constraint
  tran_condition => pm_well%realization%transport_conditions%first
  do
    if (.not.associated(tran_condition)) exit
      if (trim(tran_condition%name) == &
          trim(pm_well%well%tran_condition_name)) exit
    tran_condition => tran_condition%next
  enddo
  cur_constraint => tran_condition%cur_constraint_coupler%constraint
  select type(constraint=>cur_constraint)
    class is (tran_constraint_nwt_type)
      if (any(constraint%nwt_species%constraint_type /=  &
          CONSTRAINT_AQ_EQUILIBRIUM)) then
        pm_well%option%io_buffer = 'TRANSPORT_CONDITION ' // &
          trim(pm_well%well%tran_condition_name) // ' for WELLBORE_MODEL,&
          &WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE CONSTRAINT must be of &
          &type "AQ".'
        call PrintErrMsg(pm_well%option)
      endif
      pm_well%well%aqueous_conc_th = constraint%nwt_species%constraint_conc
  end select

end subroutine PMWellUpdatePropertiesTran

! ************************************************************************** !
subroutine PMWellCopyWell(well,well_copy,transport)
  !
  ! Copies well object properties from one to another.
  !
  ! Author: Michael Nole
  ! Date: 02/25/2022
  !

  implicit none

  type(well_type) :: well
  type(well_type) :: well_copy
  PetscBool :: transport

  well_copy%liq%den(:) = well%liq%den0
  well_copy%liq%visc(:) = well%liq%visc(:)
  well_copy%liq%kr = well%liq%kr
  well_copy%gas%den(:) = well%gas%den0
  well_copy%gas%visc(:) = well%gas%visc(:)
  well_copy%gas%kr = well%gas%kr
  well_copy%diameter(:) = well%diameter(:)
  well_copy%WI_base(:) = well%WI_base(:)
  well_copy%permeability(:) = well%permeability(:)
  well_copy%phi(:) = well%phi(:)
  well_copy%friction_factor(:) = well%friction_factor(:)
  well_copy%area(:) = well%area(:)
  well_copy%volume(:) = well%volume(:)
  well_copy%bh_p = well%bh_p
  well_copy%th_p = well%th_p
  well_copy%bh_sg = well%bh_sg
  well_copy%th_sg = well%th_sg
  well_copy%bh_ql = well%bh_ql
  well_copy%bh_qg = well%bh_qg
  well_copy%th_ql = well%th_ql
  well_copy%th_qg = well%th_qg
  well_copy%liq%visc(:) = well%liq%visc(:)
  well_copy%gas%visc(:) = well%gas%visc(:)
  call PMWellCopyReservoir(well%reservoir, well_copy%reservoir,transport)


end subroutine PMWellCopyWell

! ************************************************************************** !

subroutine PMWellCopyReservoir(reservoir,reservoir_copy,transport)
  !
  ! Copies reservoir object properties from one to another.
  !
  ! Author: Michael Nole
  ! Date: 05/15/2023
  !

  implicit none

  type(well_reservoir_type) :: reservoir
  type(well_reservoir_type) :: reservoir_copy
  PetscBool :: transport

  reservoir_copy%p_l(:) = reservoir%p_l(:)
  reservoir_copy%p_g(:) = reservoir%p_g(:)
  reservoir_copy%temp(:) = reservoir%temp(:)
  reservoir_copy%s_l(:) = reservoir%s_l(:)
  reservoir_copy%s_g(:) = reservoir%s_g(:)
  reservoir_copy%mobility_l(:) = reservoir%mobility_l(:)
  reservoir_copy%mobility_g(:) = reservoir%mobility_g(:)
  reservoir_copy%kr_l(:) = reservoir%kr_l(:)
  reservoir_copy%kr_g(:) = reservoir%kr_g(:)
  reservoir_copy%den_l(:) = reservoir%den_l(:)
  reservoir_copy%den_g(:) = reservoir%den_g(:)
  reservoir_copy%visc_l(:) = reservoir%visc_l(:)
  reservoir_copy%visc_g(:) = reservoir%visc_g(:)
  reservoir_copy%e_por(:) = reservoir%e_por(:)
  reservoir_copy%kx(:) = reservoir%kx(:)
  reservoir_copy%ky(:) = reservoir%ky(:)
  reservoir_copy%kz(:) = reservoir%kz(:)
  reservoir_copy%dx(:) = reservoir%dx(:)
  reservoir_copy%dy(:) = reservoir%dy(:)
  reservoir_copy%dz(:) = reservoir%dz(:)
  reservoir_copy%volume(:) = reservoir%volume(:)

  if (transport) then
    reservoir_copy%aqueous_conc(:,:) = reservoir%aqueous_conc(:,:)
    reservoir_copy%aqueous_mass(:,:) = reservoir%aqueous_mass(:,:)
  endif

end subroutine PMWellCopyReservoir
! ************************************************************************** !

subroutine PMWellSetPlotVariables(list,pm_well)
  !
  ! Adds variables to be printed for plotting.
  !
  ! Author: Jenn Frederick
  ! Date: 09/29/2022
  !
  use Output_Aux_module
  use Variables_module

  type(output_variable_list_type), pointer :: list
  class(pm_well_type) :: pm_well

  character(len=MAXWORDLENGTH) :: name,  units
  PetscInt :: i

  if (pm_well%well%output_pl) then
    name = 'Well Liq. Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                 WELL_LIQ_PRESSURE) 
  endif
  
  if (pm_well%well%output_pg) then
    name = 'Well Gas Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                 WELL_GAS_PRESSURE)
  endif

  if (pm_well%well%output_sl) then
    name = 'Well Liq. Saturation'
    units = '-'
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                 WELL_LIQ_SATURATION)
  endif

  if (pm_well%well%output_sg) then
    name = 'Well Gas Saturation'
    units = '-'
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                 WELL_GAS_SATURATION)
  endif

  if (pm_well%well%output_aqc .and. pm_well%transport) then
    do i=1,pm_well%nspecies
      name = 'Well AQ Conc. ' // trim(pm_well%well%species_names(i))
      units = 'mol/m^3-liq'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   WELL_AQ_CONC,i)
    enddo
  endif

  if (pm_well%well%output_aqm .and. pm_well%transport) then
    do i=1,pm_well%nspecies
      name = 'Well AQ Mass ' // trim(pm_well%well%species_names(i))
      units = 'mol'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   WELL_AQ_MASS,i)
    enddo  
  endif

  if (pm_well%well%liq%output_Q) then
    name = 'Well Liq. Q'
    units = 'kmol/sec'
    call OutputVariableAddToList(list,name,OUTPUT_RATE,units,WELL_LIQ_Q)
  endif

  if (pm_well%well%gas%output_Q) then
    name = 'Well Gas Q'
    units = 'kmol/sec'
    call OutputVariableAddToList(list,name,OUTPUT_RATE,units,WELL_GAS_Q)
  endif

end subroutine PMWellSetPlotVariables

! ************************************************************************** !

function PMWellOutputFilename(option)
  !
  ! Generates a filename for wellbore model output
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021
  !

  implicit none

  type(option_type), pointer :: option

  character(len=MAXSTRINGLENGTH) :: PMWellOutputFilename

  PMWellOutputFilename = trim(option%global_prefix) // &
                         trim(option%group_prefix) // '.well'

end function PMWellOutputFilename

! ************************************************************************** !

subroutine PMWellOutputHeader(pm_well)
  !
  ! Writes the header for wellbore model output file.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021
  !

  use Output_Aux_module
  use Grid_module
  use Utility_module

  implicit none

  class(pm_well_type) :: pm_well

  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename, word
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: cell_string
  PetscBool :: exist
  PetscInt :: fid
  PetscInt :: icolumn
  PetscInt :: k, j

  output_option => pm_well%realization%output_option

  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif

  fid = 555
  filename = PMWellOutputFilename(pm_well%option)
  exist = FileExists(trim(filename))
  if (pm_well%option%restart_flag .and. exist) return
  open(unit=fid,file=filename,action="write",status="replace")

  ! First write out the well grid information
  write(fid,'(a)',advance="yes") '========= WELLBORE MODEL GRID INFORMATION &
                                  &==================' 
  write(word,'(i5)') pm_well%well_grid%nsegments
  write(fid,'(a)',advance="yes") ' Number of segments: ' // trim(word) 
  write(word,'(i5)') pm_well%well_grid%nconnections 
  write(fid,'(a)',advance="yes") ' Number of connections: ' // trim(word)
  write(word,'(es10.3,es10.3,es10.3)') pm_well%well_grid%tophole(1), &
                          pm_well%well_grid%tophole(2), pm_well%well_grid%tophole(3)
  write(fid,'(a)',advance="yes") ' Top of hole (x,y,z) [m]: ' // trim(word)
  write(word,'(es10.3,es10.3,es10.3)') pm_well%well_grid%bottomhole(1), &
                    pm_well%well_grid%bottomhole(2), pm_well%well_grid%bottomhole(3)
  write(fid,'(a)',advance="yes") ' Bottom of hole (x,y,z) [m]: ' // trim(word)
  write(fid,'(a)',advance="yes") '===========================================&
                                  &================='
  write(fid,'(a)',advance="yes") ' Segment Number: Center coordinate (x,y,z) [m] '
  do k = 1,pm_well%well_grid%nsegments
    write(word,'(i4,a3,es10.3,es10.3,es10.3,a1)') k,': (', &
      pm_well%well_grid%h(k)%x,pm_well%well_grid%h(k)%y,pm_well%well_grid%h(k)%z,')'
    write(fid,'(a)',advance="yes") trim(word)
  enddo
  write(fid,'(a)',advance="yes") '===========================================&
                                  &================='   
  write(fid,'(a)',advance="yes") ' Segment Number: Segment length [m] '
  do k = 1,pm_well%well_grid%nsegments
    write(word,'(i4,a3,es10.3,a1)') k,': (',pm_well%well_grid%dh(k),')'
    write(fid,'(a)',advance="yes") trim(word)
  enddo
  write(fid,'(a)',advance="yes") '===========================================&
                                  &=================' 

  write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'
  cell_string = ''

  do k = 1,pm_well%well_grid%nsegments
    variable_string = 'Seg.#'
    units_string = ''
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'X'
    units_string = 'm'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Y'
    units_string = 'm'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Z'
    units_string = 'm'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well P-liq'
    units_string = 'Pa'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Res P-liq'
    units_string = 'Pa'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well P-gas'
    units_string = 'Pa'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Res P-gas'
    units_string = 'Pa'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well S-liq'
    units_string = '-'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well S-gas'
    units_string = '-'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well Q-liq'
    units_string = 'kg/s'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well Q-gas'
    units_string = 'kg/s'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well q-liq'
    units_string = 'm/s'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well q-gas'
    units_string = 'm/s'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well Liq Mass Bal'
    units_string = 'kmol/s'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well Liq Mass'
    units_string = 'kmol'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Well P. Index'
    units_string = '-'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    if (pm_well%transport) then
      do j = 1,pm_well%nspecies
        variable_string = 'Well Aqueous Conc. ' // &
                          trim(pm_well%well%species_names(j))
        units_string = 'mol/m^3-liq'
        call OutputWriteToHeader(fid,variable_string,units_string, &
                                 cell_string,icolumn)
        variable_string = 'Well Aqueous Mass. ' &
                           // trim(pm_well%well%species_names(j))
        units_string = 'mol'
        call OutputWriteToHeader(fid,variable_string,units_string, &
                                 cell_string,icolumn)
        variable_string = 'Res Aqueous Conc. ' // &
                          trim(pm_well%well%species_names(j))
        units_string = 'mol/m^3-liq'
        call OutputWriteToHeader(fid,variable_string,units_string, &
                                 cell_string,icolumn)
      enddo
    endif
  enddo

  close(fid)

end subroutine PMWellOutputHeader

! ************************************************************************** !

subroutine PMWellOutput(pm_well)
  !
  ! Sets up output for the wellbore process model
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021
  !

  use Output_Aux_module
  use Global_Aux_module
  use Grid_module

  implicit none

  class(pm_well_type) :: pm_well

  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: k, j

100 format(100es18.8)
101 format(1I6.1)

  option => pm_well%realization%option
  output_option => pm_well%realization%output_option

  fid = 555
  filename = PMWellOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")

  ! pm_well time is set at the end of the wellbore step????
  write(fid,100,advance="no") option%time / output_option%tconv

  do k = 1,pm_well%well_grid%nsegments
    write(fid,101,advance="no") k
    write(fid,100,advance="no") pm_well%well_grid%h(k)%x, &
                                pm_well%well_grid%h(k)%y, &
                                pm_well%well_grid%h(k)%z, &
                                pm_well%well%pl(k), &
                                pm_well%well%reservoir%p_l(k), &
                                pm_well%well%pg(k), &
                                pm_well%well%reservoir%p_g(k), &
                                pm_well%well%liq%s(k), &
                                pm_well%well%gas%s(k), &
                                pm_well%well%liq%Q(k), &
                                pm_well%well%gas%Q(k)
    if (k == 1) then
      write(fid,100,advance="no") pm_well%well%ql_bc(1), &
                                  pm_well%well%qg_bc(1)
    endif
    if (k > 1) then
      write(fid,100,advance="no") pm_well%well%ql(k-1), &
                                  pm_well%well%qg(k-1)
    endif
    write(fid,100,advance="no") pm_well%well%mass_balance_liq(k), &
                                pm_well%well%liq_mass(k), &
                                pm_well%well%WI(k)
    if (pm_well%transport) then
      do j = 1,pm_well%nspecies
        write(fid,100,advance="no") pm_well%well%aqueous_conc(j,k), &
                                    pm_well%well%aqueous_mass(j,k), &
                                    pm_well%well%reservoir%aqueous_conc(j,k)
      enddo
    endif
  enddo

  close(fid)

end subroutine PMWellOutput

! ************************************************************************** !

subroutine PMWellMassBalance(pm_well)
  !
  ! Verify the mass balance in the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 06/28/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  type(well_type), pointer :: well
  PetscInt :: isegment, nsegments
  PetscInt :: n_up, n_dn
  PetscReal :: mass_rate_up, mass_rate_dn

  ! q in [m3-liq/m2-bulk-sec]
  ! area in [m2-bulk]

  ! (-) mass balance rate means net mass is being lost
  ! (+) mass balance rate means net mass is being gained

  well => pm_well%well
  nsegments = pm_well%well_grid%nsegments

  if (pm_well%well_comm%comm == MPI_COMM_NULL) then 
    well%mass_balance_liq = UNINITIALIZED_DOUBLE
    return
  endif

  n_dn = +1
  n_up = -1

  do isegment = 1,nsegments

    if ((isegment > 1) .and. (isegment < nsegments)) then
    ! ----------------------------------------INTERIOR-FLUXES----------------

      ! [kmol/sec] 
      mass_rate_up = well%ql_kmol(isegment)
      mass_rate_dn = well%ql_kmol(isegment-1)

      !WRITE(*,*) 'isegment = ', isegment, ' ----------------------'
      !WRITE(*,*) 'q_up =', mass_rate_up, ' kmol/sec'
      !WRITE(*,*) 'q_dn =', mass_rate_dn, ' kmol/sec'
      !WRITE(*,*) 'Q =', well%liq%Q(isegment), ' kmol/sec'

    ! ----------------------------------------BOUNDARY-FLUXES----------------
    else if (isegment == 1) then
    ! ----- bottom of well -----

      ! [kmol/sec] 
      mass_rate_up = well%ql_kmol(isegment)
      mass_rate_dn = well%ql_kmol_bc(1)

      !WRITE(*,*) 'isegment = ', isegment, ' ----------------------'
      !WRITE(*,*) 'q_up =', mass_rate_up, ' kmol/sec'
      !WRITE(*,*) 'q_dn =', mass_rate_dn, ' kmol/sec'
      !WRITE(*,*) 'Q =', well%liq%Q(isegment), ' kmol/sec'
      !WRITE(*,*) 'BOTTOM ratio (Q/up) =', well%liq%Q(isegment)/mass_rate_up

    else if (isegment == pm_well%well_grid%nsegments) then
    ! ----- top of well -----

      ! [kmol/sec] 
      mass_rate_up = well%ql_kmol_bc(2)
      mass_rate_dn = well%ql_kmol(isegment-1)

      !WRITE(*,*) 'isegment = ', isegment, ' ----------------------'
      !WRITE(*,*) 'q_up =', mass_rate_up, ' kmol/sec'
      !WRITE(*,*) 'q_dn =', mass_rate_dn, ' kmol/sec'
      !WRITE(*,*) 'Q =', well%liq%Q(isegment), ' kmol/sec'
      !WRITE(*,*) 'TOP ratio (up/dn) =', mass_rate_up/mass_rate_dn

    endif

    well%mass_balance_liq(isegment) = mass_rate_up - mass_rate_dn &
                                      - well%liq%Q(isegment)

    !WRITE(*,*) 'MASS BALANCE: ', well%mass_balance_liq(isegment), ' kmol/sec'

  enddo

end subroutine PMWellMassBalance

! ************************************************************************** !

subroutine PMWellUpdateMass(pm_well)
  !
  ! Verify the mass balance in the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/17/2022
  !

  implicit none

  class(pm_well_type) :: pm_well

  type(well_type), pointer :: well
  PetscInt :: isegment, nsegments
  PetscReal :: inst_mass

  well => pm_well%well
  nsegments = pm_well%well_grid%nsegments

  if (pm_well%well_comm%comm == MPI_COMM_NULL) then
    well%liq_mass = UNINITIALIZED_DOUBLE
    well%liq_cum_mass = UNINITIALIZED_DOUBLE
    return
  endif

  do isegment = 1,nsegments

    ! [kmol]
    inst_mass = well%volume(isegment)*well%phi(isegment)* &
                well%liq%s(isegment)*well%liq%den(isegment)*FMWH2O

    well%liq_mass(isegment) = inst_mass
    well%liq_cum_mass(isegment) = well%liq_cum_mass(isegment) + inst_mass

  enddo

end subroutine PMWellUpdateMass

! *************************************************************************** !

subroutine PMWellInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 01/03/2023
  !

  implicit none

  class(pm_well_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMWellInputRecord

! ************************************************************************** !

subroutine PMWellDestroy(this)
  !
  ! Destroys the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !
  
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_well_type) :: this

  PetscInt :: s, i

  call PMBaseDestroy(this)

  call WellGridDestroy(this%well_grid)

  call DeallocateArray(this%well%reservoir%p_l)
  call DeallocateArray(this%well%reservoir%p_g)
  call DeallocateArray(this%well%reservoir%s_l)
  call DeallocateArray(this%well%reservoir%s_g)
  call DeallocateArray(this%well%reservoir%temp)
  call DeallocateArray(this%well%reservoir%mobility_l)
  call DeallocateArray(this%well%reservoir%mobility_g)
  call DeallocateArray(this%well%reservoir%kr_l)
  call DeallocateArray(this%well%reservoir%kr_g)
  call DeallocateArray(this%well%reservoir%den_l)
  call DeallocateArray(this%well%reservoir%den_g)
  call DeallocateArray(this%well%reservoir%visc_l)
  call DeallocateArray(this%well%reservoir%visc_g)
  call DeallocateArray(this%well%reservoir%e_por)
  call DeallocateArray(this%well%reservoir%kx)
  call DeallocateArray(this%well%reservoir%ky)
  call DeallocateArray(this%well%reservoir%kz)
  call DeallocateArray(this%well%reservoir%dx)
  call DeallocateArray(this%well%reservoir%dy)
  call DeallocateArray(this%well%reservoir%dz)
  call DeallocateArray(this%well%reservoir%volume)
  if (this%transport) then
    call DeallocateArray(this%well%reservoir%aqueous_conc)
    call DeallocateArray(this%well%reservoir%aqueous_mass)
  endif
  nullify(this%well%reservoir)

  call DeallocateArray(this%well%reservoir_save%p_l)
  call DeallocateArray(this%well%reservoir_save%p_g)
  call DeallocateArray(this%well%reservoir_save%s_l)
  call DeallocateArray(this%well%reservoir_save%s_g)
  call DeallocateArray(this%well%reservoir_save%temp)
  call DeallocateArray(this%well%reservoir_save%mobility_l)
  call DeallocateArray(this%well%reservoir_save%mobility_g)
  call DeallocateArray(this%well%reservoir_save%kr_l)
  call DeallocateArray(this%well%reservoir_save%kr_g)
  call DeallocateArray(this%well%reservoir_save%den_l)
  call DeallocateArray(this%well%reservoir_save%den_g)
  call DeallocateArray(this%well%reservoir_save%visc_l)
  call DeallocateArray(this%well%reservoir_save%visc_g)
  call DeallocateArray(this%well%reservoir_save%e_por)
  call DeallocateArray(this%well%reservoir_save%kx)
  call DeallocateArray(this%well%reservoir_save%ky)
  call DeallocateArray(this%well%reservoir_save%kz)
  call DeallocateArray(this%well%reservoir_save%dx)
  call DeallocateArray(this%well%reservoir_save%dy)
  call DeallocateArray(this%well%reservoir_save%dz)
  call DeallocateArray(this%well%reservoir_save%volume)
  if (this%transport) then
    call DeallocateArray(this%well%reservoir_save%aqueous_conc)
    call DeallocateArray(this%well%reservoir_save%aqueous_mass)
  endif
  nullify(this%well%reservoir_save)


  call DeallocateArray(this%well%area)
  call DeallocateArray(this%well%diameter)
  call DeallocateArray(this%well%volume)
  call DeallocateArray(this%well%friction_factor)
  call DeallocateArray(this%well%WI)
  call DeallocateArray(this%well%pl)
  call DeallocateArray(this%well%pg)
  call DeallocateArray(this%well%ql)
  call DeallocateArray(this%well%qg)
  call DeallocateArray(this%well%ql_bc)
  call DeallocateArray(this%well%qg_bc)
  if (this%transport) then
    call DeallocateArray(this%well%aqueous_conc)
    call DeallocateArray(this%well%aqueous_mass)
  endif
  call DeallocateArray(this%well%liq%den)
  call DeallocateArray(this%well%liq%visc)
  call DeallocateArray(this%well%liq%s)
  call DeallocateArray(this%well%liq%Q)
  call DeallocateArray(this%well%gas%den)
  call DeallocateArray(this%well%gas%visc)
  call DeallocateArray(this%well%gas%s)
  call DeallocateArray(this%well%gas%Q)
  nullify(this%well%liq)
  nullify(this%well%gas)
  nullify(this%well)

  s = size(this%well_pert)

  do i = 1,s
    call DeallocateArray(this%well_pert(i)%area)
    call DeallocateArray(this%well_pert(i)%diameter)
    call DeallocateArray(this%well_pert(i)%volume)
    call DeallocateArray(this%well_pert(i)%friction_factor)
    call DeallocateArray(this%well_pert(i)%WI)
    call DeallocateArray(this%well_pert(i)%pl)
    call DeallocateArray(this%well_pert(i)%pg)
    call DeallocateArray(this%well_pert(i)%liq%den)
    call DeallocateArray(this%well_pert(i)%liq%visc)
    call DeallocateArray(this%well_pert(i)%liq%s)
    call DeallocateArray(this%well_pert(i)%liq%Q)
    call DeallocateArray(this%well_pert(i)%gas%den)
    call DeallocateArray(this%well_pert(i)%gas%visc)
    call DeallocateArray(this%well_pert(i)%gas%s)
    call DeallocateArray(this%well_pert(i)%gas%Q)
    nullify(this%well_pert(i)%liq)
    nullify(this%well_pert(i)%gas)
    call DeallocateArray(this%well_pert(i)%reservoir%p_l)
    call DeallocateArray(this%well_pert(i)%reservoir%p_g)
    call DeallocateArray(this%well_pert(i)%reservoir%s_l)
    call DeallocateArray(this%well_pert(i)%reservoir%s_g)
    call DeallocateArray(this%well_pert(i)%reservoir%temp)
    call DeallocateArray(this%well_pert(i)%reservoir%mobility_l)
    call DeallocateArray(this%well_pert(i)%reservoir%mobility_g)
    call DeallocateArray(this%well_pert(i)%reservoir%kr_l)
    call DeallocateArray(this%well_pert(i)%reservoir%kr_g)
    call DeallocateArray(this%well_pert(i)%reservoir%den_l)
    call DeallocateArray(this%well_pert(i)%reservoir%den_g)
    call DeallocateArray(this%well_pert(i)%reservoir%visc_l)
    call DeallocateArray(this%well_pert(i)%reservoir%visc_g)
    call DeallocateArray(this%well_pert(i)%reservoir%e_por)
    call DeallocateArray(this%well_pert(i)%reservoir%kx)
    call DeallocateArray(this%well_pert(i)%reservoir%ky)
    call DeallocateArray(this%well_pert(i)%reservoir%kz)
    call DeallocateArray(this%well_pert(i)%reservoir%dx)
    call DeallocateArray(this%well_pert(i)%reservoir%dy)
    call DeallocateArray(this%well_pert(i)%reservoir%dz)
    call DeallocateArray(this%well_pert(i)%reservoir%volume)
    if (this%transport) then
      call DeallocateArray(this%well_pert(i)%reservoir%aqueous_conc)
      call DeallocateArray(this%well_pert(i)%reservoir%aqueous_mass)
    endif
    nullify(this%well_pert(i)%reservoir)
    call DeallocateArray(this%well_pert(i)%reservoir_save%p_l)
    call DeallocateArray(this%well_pert(i)%reservoir_save%p_g)
    call DeallocateArray(this%well_pert(i)%reservoir_save%s_l)
    call DeallocateArray(this%well_pert(i)%reservoir_save%s_g)
    call DeallocateArray(this%well_pert(i)%reservoir_save%temp)
    call DeallocateArray(this%well_pert(i)%reservoir_save%mobility_l)
    call DeallocateArray(this%well_pert(i)%reservoir_save%mobility_g)
    call DeallocateArray(this%well_pert(i)%reservoir_save%kr_l)
    call DeallocateArray(this%well_pert(i)%reservoir_save%kr_g)
    call DeallocateArray(this%well_pert(i)%reservoir_save%den_l)
    call DeallocateArray(this%well_pert(i)%reservoir_save%den_g)
    call DeallocateArray(this%well_pert(i)%reservoir_save%visc_l)
    call DeallocateArray(this%well_pert(i)%reservoir_save%visc_g)
    call DeallocateArray(this%well_pert(i)%reservoir_save%e_por)
    call DeallocateArray(this%well_pert(i)%reservoir_save%kx)
    call DeallocateArray(this%well_pert(i)%reservoir_save%ky)
    call DeallocateArray(this%well_pert(i)%reservoir_save%kz)
    call DeallocateArray(this%well_pert(i)%reservoir_save%dx)
    call DeallocateArray(this%well_pert(i)%reservoir_save%dy)
    call DeallocateArray(this%well_pert(i)%reservoir_save%dz)
    call DeallocateArray(this%well_pert(i)%reservoir_save%volume)
    if (this%transport) then
      call DeallocateArray(this%well_pert(i)%reservoir_save%aqueous_conc)
      call DeallocateArray(this%well_pert(i)%reservoir_save%aqueous_mass)
    endif
    nullify(this%well_pert(i)%reservoir_save)
  enddo
  nullify(this%well_pert)

  call DeallocateArray(this%flow_soln%residual)
  call DeallocateArray(this%flow_soln%Jacobian)
  call DeallocateArray(this%flow_soln%update)
  call DeallocateArray(this%flow_soln%prev_soln%pl)
  call DeallocateArray(this%flow_soln%prev_soln%sg)
  call DeallocateArray(this%flow_soln%soln_save%pl)
  call DeallocateArray(this%flow_soln%soln_save%sg)
  nullify(this%flow_soln)

  if (this%transport) then
    call DeallocateArray(this%tran_soln%residual)
    call DeallocateArray(this%tran_soln%Jacobian)
    call DeallocateArray(this%tran_soln%update)
    call DeallocateArray(this%tran_soln%prev_soln%aqueous_conc)
    call DeallocateArray(this%tran_soln%prev_soln%aqueous_mass)
    nullify(this%tran_soln)
  endif

  nullify(this%strata_list)

end subroutine PMWellDestroy

! ************************************************************************** !

end module PM_Well_class
