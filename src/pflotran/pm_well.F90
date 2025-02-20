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
  use Condition_module

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


  PetscBool, public :: initialize_well_flow = PETSC_TRUE
  PetscBool, public :: well_output = PETSC_FALSE
  PetscReal, public :: min_flow_dt_scale = 1.d-3

  PetscInt, parameter, public :: PEACEMAN_NONE = 0
  PetscInt, parameter, public :: PEACEMAN_ISO = 1
  PetscInt, parameter, public :: PEACEMAN_2D = 2
  PetscInt, parameter, public :: PEACEMAN_3D = 3

  ! Well model types
  PetscInt, parameter, public :: WELL_MODEL_WIPP_SEQUENTIAL = 1
  PetscInt, parameter, public :: WELL_MODEL_WIPP_QI = 2
  PetscInt, parameter, public :: WELL_MODEL_FULL_MOMENTUM = 3
  PetscInt, parameter, public :: WELL_MODEL_HYDROSTATIC = 4
  PetscInt, parameter, public :: WELL_MODEL_U_SHAPE = 5
  PetscInt, parameter, public :: WELL_MODEL_COAXIAL = 6

  ! Well constraint types
  PetscInt, parameter :: CONSTANT_PRESSURE = 1
  PetscInt, parameter :: CONSTANT_PRESSURE_HYDROSTATIC = 2
  PetscInt, parameter :: CONSTANT_RATE = 3

  ! Well coupling options
  PetscInt, parameter, public :: FULLY_IMPLICIT_WELL = ONE_INTEGER
  PetscInt, parameter, public :: QUASI_IMPLICIT_WELL = TWO_INTEGER
  PetscInt, parameter, public :: SEQUENTIAL_WELL = THREE_INTEGER


  type, public :: well_reservoir_type
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
    ! reservoir liquid enthalpy
    PetscReal, pointer :: H_l(:)
    ! reservoir gas enthalpy
    PetscReal, pointer :: H_g(:)
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
    ! temp array
    PetscReal, pointer :: tmp_flow(:,:)
    PetscReal, pointer :: tmp_tran(:,:,:)
  end type

  type, public :: well_type
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
    ! flag for .well file segment output
    PetscBool, pointer :: output_in_well_file(:)
    ! list of segments for printing in the .well file
    PetscInt, pointer :: segments_for_output(:)
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
    ! well species q cumulative aqueous mass [mol] (ispecies,mass@segment)
    PetscReal, pointer :: aqueous_mass_q_cumulative(:,:)
    ! well species Q cumulative aqueous mass [mol] (ispecies,mass@segment)
    PetscReal, pointer :: aqueous_mass_Qcumulative(:,:)
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
    PetscReal :: total_rate
    ! well transport constraint name
    character(len=MAXWORDLENGTH) :: tran_condition_name
    ! Link to characteristic curves
    PetscInt, pointer :: ccid(:)
    ! permeability along the well [m2]
    PetscReal, pointer :: permeability(:)
    ! porosity
    PetscReal, pointer :: phi(:)
  end type well_type

  type, public :: well_fluid_type
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
    ! fluid enthalpy
    PetscReal, pointer :: H(:)
    ! fluid source/sink in/out of well [kmol/s]
    PetscReal, pointer :: Q(:)
    ! flag for output
    PetscBool :: output_Q
    ! cumulative fluid source/sink in/out of well [kmol/s]
    PetscReal, pointer :: Q_cumulative(:)
  end type well_fluid_type

  ! primary variables necessary to reset flow solution
  type, public :: well_flow_save_type
    PetscReal, pointer :: pl(:)
    PetscReal, pointer :: sg(:)
    PetscReal :: bh_p
  end type well_flow_save_type

  ! primary variables necessary to reset transport solution
  type, public :: well_tran_save_type
    PetscReal, pointer :: aqueous_conc(:,:)
    PetscReal, pointer :: aqueous_mass(:,:)
    PetscReal, pointer :: resr_aqueous_conc(:,:)
  end type well_tran_save_type

  type, public :: well_soln_base_type
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

  type, public, extends(well_soln_base_type) :: well_soln_flow_type
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

  type, public, extends(well_soln_base_type) :: well_soln_tran_type
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

  type, public :: well_comm_type
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
    type(strata_list_type), pointer :: strata_list
    type(well_comm_type), pointer :: well_comm
    type(flow_condition_type), pointer :: flow_condition
    PetscReal, pointer :: pert(:,:)
    PetscReal, pointer :: srcsink_water(:,:)
    PetscReal, pointer :: srcsink_gas(:,:)
    PetscBool :: split_output_file
    PetscInt :: flow_coupling
    PetscBool :: print_well
    PetscBool :: print_output
    PetscBool :: use_well_coupler
    PetscBool :: pressure_controlled
    PetscReal :: pressure_threshold_min
    PetscReal :: pressure_threshold_max
    class(pm_well_type), pointer :: next_well

  contains
    procedure, public :: SetRealization => PMWellSetRealization
    procedure, public :: FinalizeRun => PMWellFinalizeRun
    procedure, public :: UpdateTimestep => PMWellUpdateTimestep
    procedure, public :: PreSolve => PMWellPreSolve
    procedure, public :: PostSolve => PMWellPostSolve
    procedure, public :: InputRecord => PMWellInputRecord
    procedure, public :: Perturb => PMWellBasePerturb
    ! These get extended based off well model type
    procedure, public :: ReadPMBlock => PMWellReadPMBlockBase
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMWellReadSimOptionsBlockBase
    procedure, public :: InitializeRun => PMWellInitializeRunBase
    procedure, public :: InitializeTimestep => PMWellInitializeTimestepBase
    procedure, public :: Setup => PMWellSetupBase
    procedure, public :: Solve => PMWellSolveBase
    procedure, public :: SolveFlow => PMWellSolveFlowBase
    procedure, public :: ModifyFlowResidual => PMWellModifyFlowResBase
    procedure, public :: ModifyFlowJacobian => PMWellModifyFlowJacBase
    procedure, public :: UpdateFlowRates => PMWellUpdateFlowRatesBase
    procedure, public :: UpdateFlowProperties => PMWellUpdateFlowPropertiesBase
    procedure, public :: FinalizeTimestep => PMWellFinalizeTimestepBase
    procedure, public :: Destroy => PMWellDestroyBase
  end type pm_well_type

  ! Two base well model types: sequential and implicit

  ! Sequential well types have their own flow and transport solvers
  ! and interface with the reservoir through source/sink terms
  type, public, extends(pm_well_type) :: pm_well_sequential_type
    PetscInt :: nphase
    PetscInt :: nspecies
    PetscReal :: min_dt_flow       ! [sec]
    PetscReal :: min_dt_tran       ! [sec]
    PetscReal :: dt_flow           ! [sec]
    PetscReal :: dt_tran           ! [sec]
    PetscBool :: transport
    PetscBool :: ss_check
    PetscReal :: cumulative_dt_flow             ! [sec]
    PetscReal :: cumulative_dt_tran             ! [sec]
    PetscReal :: intrusion_time_start           ! [sec]
    PetscReal :: bh_zero_value                  ! [mol/m3-bulk]
    PetscInt  :: well_force_ts_cut
    PetscBool :: well_on    !Turns the well on, regardless of other checks
    type(well_soln_flow_type), pointer :: flow_soln
    type(well_soln_tran_type), pointer :: tran_soln
  contains
    procedure, public :: ReadPMBlock => PMWellReadPMBlockSeq
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMWellReadSimOptionsBlockSeq
    procedure, public :: InitializeRun => PMWellInitializeRunSequential
    procedure, public :: InitializeTimestep => PMWellInitializeTimestepSeq
    procedure, public :: Setup => PMWellSetupSequential
    procedure, public :: Solve => PMWellSolveSequential
    procedure, public :: SolveFlow => PMWellSolveFlowSequential
    procedure, public :: FinalizeTimestep => PMWellFinalizeTimestepSequential
    procedure, public :: Destroy => PMWellDestroySequential
  end type pm_well_sequential_type

  ! Implicit well types have their own entries in the reservoir
  ! residual and Jacobian
  type, public, extends(pm_well_type) :: pm_well_implicit_type
  contains
    procedure, public :: ReadPMBlock => PMWellReadPMBlockImplicit
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMWellReadSimOptionsBlockImp
    procedure, public :: Setup => PMWellSetupImplicit
    procedure, public :: Solve => PMWellSolveImplicit
    procedure, public :: SolveFlow => PMWellSolveFlowImplicit
  end type pm_well_implicit_type

  ! Extensions of the sequential well model type

  ! Quasi-implicit well model
  type, public, extends(pm_well_sequential_type) :: pm_well_qi_type
    PetscBool :: update_for_flow_qi_coupling
  contains
    procedure, public :: ReadPMBlock => PMWellReadPMBlockQI
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMWellReadSimOptionsBlockQI
    procedure, public :: Setup => PMWellSetupQI
    procedure, public :: Solve => PMWellSolveQI
    procedure, public :: SolveFlow => PMWellSolveFlowQI
  end type pm_well_qi_type

  ! Extensions of the implicit well model type

  ! Hyrostatic well model
  type, public, extends(pm_well_implicit_type) :: pm_well_hydrostatic_type
  contains
    procedure, public :: ReadPMBlock => PMWellReadPMBlockHydrostatic
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMWellReadSimOptionsBlockHyd
    procedure, public :: InitializeRun => PMWellInitializeRunHydrostatic
    procedure, public :: Setup => PMWellSetupHydrostatic
    procedure, public :: Solve => PMWellSolveHydrostatic
    procedure, public :: SolveFlow => PMWellSolveFlowHydrostatic
    procedure, public :: ModifyFlowResidual => PMWellModifyFlowResHydrostatic
    procedure, public :: ModifyFlowJacobian => PMWellModifyFlowJacHydrostatic
    procedure, public :: UpdateFlowRates => PMWellUpdateFlowRatesHydrostatic
    procedure, public :: FinalizeTimestep => PMWellFinalizeTimestepHydrostatic
    procedure, public :: Destroy => PMWellDestroyHydrostatic
  end type pm_well_hydrostatic_type

  ! Closed Loop Well Model
  type, public, extends(pm_well_implicit_type) :: pm_well_closed_loop_type
    ! Placeholder for now
  end type pm_well_closed_loop_type

  ! Extensions of the closed loop well model type

  ! U-shaped well model
  type, public, extends(pm_well_closed_loop_type) :: pm_well_u_shape_type
    ! Placeholder for now
  end type pm_well_u_shape_type

  ! Coaxial well model
  type, public, extends(pm_well_closed_loop_type) :: pm_well_coaxial_type
    ! Placeholder for now
  end type pm_well_coaxial_type

  public :: PMWellHydrostaticCreate, &
            PMWellUShapeCreate, &
            PMWellCoaxialCreate, &
            PMWellSetupGrid, &
            PMWellReadGrid, &
            PMWellReadPass2, &
            PMWellUpdateReservoirSrcSinkFlow, &
            PMWellModifyDummyFlowJacobian, &
            PMWellCopyWell, &
            PMWellReadWellOutput, &
            PMWellComputeWellIndex, &
            PMWellOutputSequential, &
            PMWellUpdateReservoirSrcSinkTran, &
            PMWellUpdateReservoirConcTran, &
            PMWellMassBalance, &
            PMWellCalcCumulativeQFlux, &
            PMWellCalcCumulativeTranFlux, &
            PMWellInitializeWellFlow, &
            PMWellUpdateStrata, &
            PMWellReadWell, &
            PMWellReadWellBCs, &
            PMWellReadFlowSolver, &
            PMWellReadTranSolver, &
            PMWellReadWellConstraintType, &
            PMWellReadSimOptionsBlockBase, &
            PMWellFlowCreate, &
            PMWellTranCreate, &
            PMWellBaseInit, &
            PMWellCopyReservoir, &
            PMWellVarCreate, &
            PMWellSequentialInit, &
            PMWellQIInit
  contains

! ************************************************************************** !

subroutine PMWellBaseInit(pm_well)
  !
  ! Creates the well process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: pm_well

  call PMBaseInit(pm_well)

  pm_well%header = 'WELLBORE MODEL'
  allocate(pm_well%well)
  allocate(pm_well%well_comm)
  nullify(pm_well%realization)
  pm_well%well_grid => WellGridCreate()
  call PMWellVarCreate(pm_well%well)

  nullify(pm_well%flow_condition)
  nullify(pm_well%well_comm%petsc_rank_list)
  nullify(pm_well%well_comm%well_rank_list)
  pm_well%well_comm%petsc_rank = 0
  pm_well%well_comm%comm = 0
  pm_well%well_comm%group = 0
  pm_well%well_comm%rank = 0
  pm_well%well_comm%commsize = 0

  nullify(pm_well%pert)
  nullify(pm_well%srcsink_water)
  nullify(pm_well%srcsink_gas)

  pm_well%split_output_file = PETSC_FALSE

  pm_well%flow_coupling = UNINITIALIZED_INTEGER
  pm_well%print_well = PETSC_FALSE
  pm_well%print_output = PETSC_FALSE
  pm_well%use_well_coupler = PETSC_FALSE
  pm_well%pressure_controlled = PETSC_FALSE
  pm_well%pressure_threshold_min = 0.d0
  pm_well%pressure_threshold_max = MAX_DOUBLE

end subroutine PMWellBaseInit

! ************************************************************************** !

subroutine PMWellSequentialInit(pm_well)
  !
  ! Creates the sequential well process model.
  !
  ! Author: Michael Nole
  ! Date: 02/19/2025

  implicit none

  class(pm_well_sequential_type) :: pm_well

  call PMWellBaseInit(pm_well)
  call PMWellFlowCreate(pm_well)
  call PMWellTranCreate(pm_well)

  ! strata list specific to well
  allocate(pm_well%strata_list)
  nullify(pm_well%strata_list%first)
  nullify(pm_well%strata_list%last)
  nullify(pm_well%strata_list%array)
  pm_well%strata_list%num_strata = 0

  pm_well%nphase = 0
  pm_well%nspecies = 0
  pm_well%dt_flow = 0.d0
  pm_well%dt_tran = 0.d0
  pm_well%min_dt_flow = 1.d-15
  pm_well%min_dt_tran = 1.d-15
  pm_well%cumulative_dt_flow = 0.d0
  pm_well%cumulative_dt_tran = 0.d0
  pm_well%transport = PETSC_FALSE
  pm_well%ss_check = PETSC_FALSE
  pm_well%split_output_file = PETSC_FALSE
  pm_well%well_on = PETSC_TRUE
  pm_well%well_force_ts_cut = 0
  pm_well%flow_coupling = UNINITIALIZED_INTEGER
  pm_well%print_well = PETSC_FALSE
  pm_well%print_output = PETSC_FALSE
  pm_well%use_well_coupler = PETSC_FALSE

  pm_well%split_output_file = PETSC_FALSE
  pm_well%print_well = PETSC_FALSE
  pm_well%print_output = PETSC_FALSE

end subroutine PMWellSequentialInit

! ************************************************************************** !

subroutine PMWellQIInit(pm_well)
  !
  ! Creates the sequential well process model.
  !
  ! Author: Michael Nole
  ! Date: 02/19/2025

  implicit none

  class(pm_well_qi_type) :: pm_well

  call PMWellSequentialInit(pm_well)
  pm_well%update_for_flow_qi_coupling = PETSC_FALSE

end subroutine PMWellQIInit

! ************************************************************************** !

function PMWellHydrostaticCreate()
  !
  ! Creates the hydrostatic well process model.
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025

  implicit none

  class(pm_well_hydrostatic_type), pointer :: PMWellHydrostaticCreate
  class(pm_well_hydrostatic_type), pointer :: pm_well

  allocate(pm_well)
  call PMWellBaseInit(pm_well)

  pm_well%well%well_model_type = WELL_MODEL_HYDROSTATIC
  pm_well%flow_coupling = FULLY_IMPLICIT_WELL

  nullify(pm_well%next_well)

  ! strata list specific to well
  allocate(pm_well%strata_list)
  nullify(pm_well%strata_list%first)
  nullify(pm_well%strata_list%last)
  nullify(pm_well%strata_list%array)
  pm_well%strata_list%num_strata = 0

  pm_well%split_output_file = PETSC_FALSE
  pm_well%print_well = PETSC_FALSE
  pm_well%print_output = PETSC_FALSE

  PMWellHydrostaticCreate => pm_well

end function PMWellHydrostaticCreate

! ************************************************************************** !

function PMWellUShapeCreate()
  !
  ! Creates the U-shape well process model.
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025

  implicit none

  class(pm_well_u_shape_type), pointer :: PMWellUShapeCreate
  class(pm_well_u_shape_type), pointer :: pm_well

  allocate(pm_well)
  call PMBaseInit(pm_well)
  allocate(pm_well%well)
  allocate(pm_well%well_comm)
  call PMWellBaseInit(pm_well)

  pm_well%well%well_model_type = WELL_MODEL_U_SHAPE
  pm_well%flow_coupling = FULLY_IMPLICIT_WELL

  nullify(pm_well%next_well)

  PMWellUShapeCreate => pm_well

end function PMWellUShapeCreate

! ************************************************************************** !

function PMWellCoaxialCreate()
  !
  ! Creates the U-shape well process model.
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025

  implicit none

  class(pm_well_coaxial_type), pointer :: PMWellCoaxialCreate
  class(pm_well_coaxial_type), pointer :: pm_well

  allocate(pm_well)
  call PMBaseInit(pm_well)
  allocate(pm_well%well)
  allocate(pm_well%well_comm)
  call PMWellBaseInit(pm_well)

  pm_well%well%well_model_type = WELL_MODEL_COAXIAL
  pm_well%flow_coupling = FULLY_IMPLICIT_WELL

  nullify(pm_well%next_well)

  PMWellCoaxialCreate => pm_well

end function PMWellCoaxialCreate

! ************************************************************************** !

subroutine PMWellFlowCreate(pm_well)
  !
  ! Creates flow variables.
  !
  ! Author: Michael Nole
  ! Date: 03/04/2023

  implicit none

  class(pm_well_sequential_type) :: pm_well

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

  class(pm_well_sequential_type) :: pm_well

  ! create the well transport solution object:
  allocate(pm_well%tran_soln)
  nullify(pm_well%tran_soln%prev_soln%aqueous_conc)
  nullify(pm_well%tran_soln%prev_soln%aqueous_mass)
  nullify(pm_well%tran_soln%prev_soln%resr_aqueous_conc)

  pm_well%tran_soln%ndof = UNINITIALIZED_INTEGER
  nullify(pm_well%tran_soln%residual)
  nullify(pm_well%tran_soln%Jacobian)
  nullify(pm_well%tran_soln%update)

  pm_well%tran_soln%not_converged = PETSC_TRUE
  pm_well%tran_soln%converged = PETSC_FALSE
  pm_well%tran_soln%cut_timestep = PETSC_FALSE
  pm_well%tran_soln%max_iter = 8
  pm_well%tran_soln%max_ts_cut = 20
  pm_well%tran_soln%ts_cut_factor = 2.0d0  ! these are not used in transport
  pm_well%tran_soln%ts_ramp_factor = 1.1d0 ! these are not used in transport
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
  nullify(reservoir%H_l)
  nullify(reservoir%H_g)
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
  nullify(reservoir%tmp_flow)
  nullify(reservoir%tmp_tran)

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
  nullify(well%output_in_well_file)
  nullify(well%segments_for_output)
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
  nullify(well%aqueous_mass_q_cumulative)
  nullify(well%aqueous_mass_Qcumulative)
  well%output_aqc = PETSC_FALSE
  well%output_aqm = PETSC_FALSE
  nullify(well%aqueous_conc_th)
  well%bh_p_set_by_reservoir = PETSC_FALSE
  well%bh_sg_set_by_reservoir = PETSC_FALSE
  well%bh_p = UNINITIALIZED_DOUBLE
  well%th_p = UNINITIALIZED_DOUBLE
  well%bh_sg = UNINITIALIZED_DOUBLE
  well%th_sg = UNINITIALIZED_DOUBLE
  !MAN: this might break regression tests
  well%bh_ql = UNINITIALIZED_DOUBLE
  well%bh_qg = UNINITIALIZED_DOUBLE
  well%th_ql = 0.d0
  well%th_qg = 0.d0
  well%total_rate = 999.d0
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
  nullify(well%liq%H)
  nullify(well%liq%Q)
  well%liq%output_Q = PETSC_FALSE
  nullify(well%liq%Q_cumulative)

  ! create the fluid/gas objects:
  allocate(well%gas)
  well%gas%ifluid = 2
  nullify(well%gas%kr)
  well%gas%den0 = UNINITIALIZED_DOUBLE
  nullify(well%gas%den)
  nullify(well%gas%visc)
  nullify(well%gas%s)
  nullify(well%gas%xmass)
  nullify(well%gas%H)
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
  nullify(well%reservoir%H_l)
  nullify(well%reservoir%H_g)
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
  nullify(well%reservoir_save%H_l)
  nullify(well%reservoir_save%H_g)
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
  use String_module

  implicit none

  type(well_grid_type), pointer :: well_grid
  type(grid_type), pointer :: res_grid
  type(option_type), pointer :: option

  type(point3d_type) :: dummy_h
  type(point3d_type), allocatable :: well_nodes(:), face_centroids(:)
  character(len=MAXSTRINGLENGTH) :: dim
  PetscReal :: diff_x,diff_y,diff_z
  PetscReal :: dh_x,dh_y,dh_z
  PetscReal :: total_length
  PetscReal :: top_of_reservoir, top_of_hole
  PetscReal :: bottom_of_reservoir, bottom_of_hole
  PetscReal :: temp_real, temp_real2
  PetscReal :: cum_z, z_dum
  PetscInt, allocatable :: temp_id_list(:)
  PetscReal, allocatable :: temp_z_list(:)
  PetscReal, allocatable :: collect_x_temp(:)
  PetscReal, allocatable :: collect_y_temp(:)
  PetscReal, allocatable :: collect_z_temp(:)
  PetscInt, allocatable :: temp_repeated_list(:), collect_rank(:)
  PetscInt :: cur_id, cum_z_int, cur_cum_z_int
  PetscInt :: repeated
  PetscInt :: num_entries
  PetscReal, allocatable :: dz_list(:), res_dz_list(:)
  PetscReal :: min_dz, dz, z
  PetscInt, allocatable :: cell_id_list(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: res_cell_count
  PetscInt :: nsegments, nsegments_save
  PetscInt :: k, i, j
  PetscInt :: local_index
  PetscReal :: ds
  PetscReal, allocatable :: well_trajectory(:,:), temp_trajectory(:,:)
  PetscBool, allocatable :: well_casing(:), temp_casing(:)
  PetscReal :: cur_location(3)
  PetscReal :: delta_segment, zmin, zmax
  PetscInt :: temp_segments, total_segments, dims, index
  PetscReal :: proj, r, angle, angle_to_horiz
  PetscReal, allocatable :: temp_dx(:), temp_dy(:), temp_dz(:)
  PetscReal, allocatable :: temp_x(:), temp_y(:), temp_z(:)
  PetscReal :: circle_center(3), surface_origin(3)
  type(deviated_well_type), pointer :: well_segment
  PetscReal, allocatable :: segment_coordinates(:,:), segment_dxyz(:,:)
  PetscReal :: progress, well_length
  PetscErrorCode :: ierr
  PetscBool :: find_segments
  PetscReal, parameter :: epsilon = 1.d-10

  find_segments = PETSC_TRUE
  num_entries = 10000

  if (associated(well_grid%deviated_well_segment_list)) then
    well_segment => well_grid%deviated_well_segment_list
    ! Collect points on a line along the well trajectory.
    ! This has been tested for downward trajectory wells.

    ! Line increment default
    ! MAN: make this a knob in the future.
    ds = 1.d-2
    dims = 3
    nsegments = 0
    well_length = 0.d0
    ! Angle from vertical (downward)
    angle = 0.d0
    proj = 0.d0
    do
      if (.not. associated(well_segment)) exit

      if (Initialized(well_segment%surface_origin(1))) then
        cur_location = well_segment%surface_origin
        surface_origin = well_segment%surface_origin
      elseif (allocated(well_segment%segment_coordinates)) then
        find_segments = PETSC_FALSE
        well_grid%save_well_segment_list = PETSC_FALSE
        nsegments = size(well_segment%segment_coordinates(:,1))
        allocate(segment_coordinates(nsegments,dims))
        allocate(well_casing(nsegments))
        allocate(segment_dxyz(nsegments,dims))
        segment_coordinates(:,:) = well_segment%segment_coordinates(:,:)
        segment_dxyz(:,:) = well_segment%segment_dxyz(:,:)
        well_casing(:) = well_segment%casing(:)
        allocate(well_trajectory(nsegments,dims))
        well_trajectory(nsegments,1) = segment_coordinates(nsegments,1) + &
                  sign(segment_dxyz(nsegments,1), &
                  segment_coordinates(nsegments,1) - &
                  segment_coordinates(nsegments-1,1)) / 2.d0
        well_trajectory(nsegments,2) = segment_coordinates(nsegments,2) + &
                  sign(segment_dxyz(nsegments,2), &
                  segment_coordinates(nsegments,2) - &
                  segment_coordinates(nsegments-1,2)) / 2.d0
        well_trajectory(nsegments,3) = segment_coordinates(nsegments,3) + &
                  sign(segment_dxyz(nsegments,3), &
                  segment_coordinates(nsegments,3) - &
                  segment_coordinates(nsegments-1,3)) / 2.d0

      elseif (Initialized(well_segment%dxyz(1)) .and. find_segments) then
        delta_segment = sqrt(well_segment%dxyz(1)**2 + &
                             well_segment%dxyz(2)**2 + &
                             well_segment%dxyz(3)**2)
        temp_segments = floor(delta_segment / ds)
        total_segments = nsegments + temp_segments
        allocate(temp_trajectory(total_segments,dims))
        allocate(temp_casing(total_segments))
        if (nsegments > 0) then
          temp_trajectory(1:nsegments,:) = well_trajectory(1:nsegments,:)
          temp_casing(1:nsegments) = well_casing(1:nsegments)
          deallocate(well_trajectory)
          deallocate(well_casing)
        endif
        allocate(well_trajectory(total_segments,dims))
        allocate(well_casing(total_segments))
        do i = nsegments + 1, nsegments + temp_segments
          temp_trajectory(i,1) = cur_location(1) + well_segment%dxyz(1) / &
                                 temp_segments
          temp_trajectory(i,2) = cur_location(2) + well_segment%dxyz(2) / &
                                 temp_segments
          temp_trajectory(i,3) = cur_location(3) + well_segment%dxyz(3) / &
                                 temp_segments
          cur_location(1) = temp_trajectory(i,1)
          cur_location(2) = temp_trajectory(i,2)
          cur_location(3) = temp_trajectory(i,3)
          temp_casing(i) = well_segment%cased
          well_length = well_length + ds
        enddo
        nsegments = nsegments + temp_segments
        well_trajectory(:,:) = temp_trajectory(:,:)
        well_casing(:) = temp_casing(:)
        deallocate(temp_trajectory)
        deallocate(temp_casing)
        proj = sqrt((cur_location(1) - surface_origin(1)) **2 + &
                    (cur_location(2) - surface_origin(2)) **2)
        dz = cur_location(3) - surface_origin(3)
        angle = atan(proj,-dz)

      elseif (Initialized(well_segment%radius_to_horizontal_x) .and. &
              find_segments) then
        r = well_segment%radius_to_horizontal_x
        circle_center(:) = cur_location(:)
        circle_center(1) = circle_center(1) + r
        delta_segment = (sign(PI / 2.d0, r) &
                        - angle) * r
        temp_segments = floor(delta_segment / ds)
        allocate(temp_trajectory(nsegments + temp_segments,3))
        allocate(temp_casing(nsegments + temp_segments))
        if (nsegments > 0) then
          temp_trajectory(1:nsegments,:) = well_trajectory(:,:)
          temp_casing(1:nsegments) = well_casing(:)
          deallocate(well_trajectory)
          deallocate(well_casing)
        endif
        allocate(well_trajectory(nsegments + temp_segments,3))
        allocate(well_casing(nsegments + temp_segments))
        j = 0
        do i = nsegments + 1, nsegments + temp_segments
          j = j+1
          temp_trajectory(i,1) = circle_center(1) - r * &
                                 cos(j * (sign(PI / 2.d0, r)-angle)/ &
                                 temp_segments)
          temp_trajectory(i,2) = cur_location(2)
          temp_trajectory(i,3) = circle_center(3) - dabs(r * &
                                 sin(j * (sign(PI / 2.d0, r)-angle)/ &
                                 temp_segments))
          cur_location(1) = temp_trajectory(i,1)
          cur_location(2) = temp_trajectory(i,2)
          cur_location(3) = temp_trajectory(i,3)
          temp_casing(i) = well_segment%cased
          well_length = well_length + ds
        enddo
        nsegments = nsegments + temp_segments
        well_trajectory(:,:) = temp_trajectory(:,:)
        well_casing(:) = temp_casing(:)
        deallocate(temp_trajectory)
        deallocate(temp_casing)

      elseif (Initialized(well_segment%radius_to_horizontal_y) .and. &
              find_segments) then
        r = well_segment%radius_to_horizontal_y
        circle_center(:) = cur_location(:)
        circle_center(2) = circle_center(2) + r
        delta_segment = (sign(PI / 2.d0, r) &
                        - angle) * r
        temp_segments = floor(delta_segment / ds)
        allocate(temp_trajectory(nsegments + temp_segments,3))
        allocate(temp_casing(nsegments + temp_segments))
        if (nsegments > 0) then
          temp_trajectory(1:nsegments,:) = well_trajectory(:,:)
          temp_casing(1:nsegments) = well_casing(:)
          deallocate(well_trajectory)
          deallocate(well_casing)
        endif
        allocate(well_trajectory(nsegments + temp_segments,3))
        allocate(well_casing(nsegments + temp_segments))
        j = 0
        do i = nsegments + 1, nsegments + temp_segments
          j = j+1
          temp_trajectory(i,1) = circle_center(1)
          temp_trajectory(i,2) = cur_location(2) - r * &
                                 cos(j * (sign(PI / 2.d0, r)-angle)/ &
                                 temp_segments)
          temp_trajectory(i,3) = circle_center(3) - dabs(r * &
                                 sin(j * (sign(PI / 2.d0, r)-angle)/ &
                                 temp_segments))
          cur_location(1) = temp_trajectory(i,1)
          cur_location(2) = temp_trajectory(i,2)
          cur_location(3) = temp_trajectory(i,3)
          temp_casing(i) = well_segment%cased
          well_length = well_length + ds
        enddo
        nsegments = nsegments + temp_segments
        well_trajectory(:,:) = temp_trajectory(:,:)
        well_casing(:) = temp_casing(:)
        deallocate(temp_trajectory)
        deallocate(temp_casing)

      elseif (Initialized(well_segment%radius_to_horizontal_angle(1)) .and. &
              find_segments) then
        r = well_segment%radius_to_horizontal_angle(1)
        angle_to_horiz = well_segment%radius_to_horizontal_angle(2)
        circle_center(:) = cur_location(:)
        circle_center(1) = circle_center(1) + r * cos(angle_to_horiz * PI / &
                           180.d0)
        circle_center(2) = circle_center(1) + r * sin(angle_to_horiz * PI / &
                           180.d0)
        delta_segment = (sign(PI / 2.d0, r) &
                        - angle) * r
        temp_segments = floor(delta_segment / ds)
        allocate(temp_trajectory(nsegments + temp_segments,3))
        allocate(temp_casing(nsegments + temp_segments))
        if (nsegments > 0) then
          temp_trajectory(1:nsegments,:) = well_trajectory(:,:)
          temp_casing(1:nsegments) = well_casing(:)
          deallocate(well_trajectory)
          deallocate(well_casing)
        endif
        allocate(well_trajectory(nsegments + temp_segments,3))
        allocate(well_casing(nsegments + temp_segments))
        j = 0
        do i = nsegments + 1, nsegments + temp_segments
          j = j+1
          temp_trajectory(i,1) = circle_center(1) - r * &
                                 cos(angle_to_horiz * PI / 180.d0) * &
                                 cos(j * (sign(PI / 2.d0, r)-angle)/ &
                                 temp_segments)
          temp_trajectory(i,2) = cur_location(2) - r * &
                                 sin(angle_to_horiz * PI / 180.d0) * &
                                 cos(j * (sign(PI / 2.d0, r)-angle)/ &
                                 temp_segments)
          temp_trajectory(i,3) = circle_center(3) - dabs(r * &
                                 sin(j * (sign(PI / 2.d0, r)-angle)/ &
                                 temp_segments))
          cur_location(1) = temp_trajectory(i,1)
          cur_location(2) = temp_trajectory(i,2)
          cur_location(3) = temp_trajectory(i,3)
          temp_casing(i) = well_segment%cased
          well_length = well_length + ds
        enddo
        nsegments = nsegments + temp_segments
        well_trajectory(:,:) = temp_trajectory(:,:)
        well_casing(:) = temp_casing(:)
        deallocate(temp_trajectory)
        deallocate(temp_casing)

      elseif (Initialized(well_segment%radius_to_vertical)) then
        option%io_buffer = "RADIUS_TO_VERTICAL is not yet supported."
        call PrintErrMsg(option)
      else
        option%io_buffer = "Error in constructing well grid using &
           &WELL_TRAJECTORY: please specify a SURFACE_ORIGIN and &
           &{SEGMENT_LIST or [DXYZ, RADIUS_TO_HORIZONTAL, &
           &RADIUS_TO_VERTICAL]}."
        call PrintErrMsg(option)
      endif
      well_segment => well_segment%next
    enddo

    well_grid%tophole(:) = surface_origin
    well_grid%bottomhole(:) = well_trajectory(nsegments,:)

    ! Now refine based off of the grid.

    allocate(temp_id_list(nsegments))
    allocate(temp_dx(nsegments))
    allocate(temp_dy(nsegments))
    allocate(temp_dz(nsegments))
    allocate(temp_x(nsegments))
    allocate(temp_y(nsegments))
    allocate(temp_z(nsegments))
    allocate(temp_casing(nsegments))
    temp_id_list = UNINITIALIZED_INTEGER
    temp_dx = UNINITIALIZED_DOUBLE
    temp_dy = UNINITIALIZED_DOUBLE
    temp_dz = UNINITIALIZED_DOUBLE
    temp_x = -MAX_DOUBLE
    temp_y =-MAX_DOUBLE
    temp_z = -MAX_DOUBLE

    k = 0
    j = 0
    cur_id = UNINITIALIZED_INTEGER
    cum_z = 0
    cur_cum_z_int = 0
    ! search procedure for finding reservoir cell list within well
    ! MAN: This needs to be tested for a case where well goes from
    ! on-process to off-process and then back on-process.

    if (find_segments) then
      option%io_buffer = 'Well Length: ' // StringWrite(well_length) // 'm'
      call PrintMsg(option)
      option%io_buffer = 'WELL_TRAJECTORY can take a long time on large &
                          &grids. Try specifying SEGMENT_LIST &
                          &to reduce time spent constructing the well. &
                          &Current progress (roughly): '
      call PrintMsg(option)
      do j = 1, nsegments
        dummy_h%x = well_trajectory(j,1)
        dummy_h%y = well_trajectory(j,2)
        dummy_h%z = well_trajectory(j,3)

        call GridGetLocalIDFromCoordinate(res_grid,dummy_h,option,local_id)
        progress = j/nsegments * well_length
        if (option%comm%io_rank == option%myrank) then
          call PrintProgressBarInt(nsegments+1.d-10,5,j)
        endif

        if (k == 0 .and. Initialized(local_id)) then
          cur_id = local_id
          k = 1
          temp_id_list(k) = local_id
          temp_casing(k) = well_casing(j)
          dh_x = dummy_h%x
          dh_y = dummy_h%y
          dh_z = dummy_h%z
        endif

        if (local_id /= cur_id) then
          k = k + 1
          if (k > nsegments) then
            option%io_buffer = 'The well intersects a lot of reservoir cells, &
                              &exceeding the current buffer size.'
            call PrintErrMsgToDev(option, &
                            'if coarsening resolution is not an option.')
          endif
          temp_id_list(k) = local_id
          dh_x = well_trajectory(j,1) - dh_x
          dh_y = well_trajectory(j,2) - dh_y
          dh_z = well_trajectory(j,3) - dh_z
          temp_dx(k-1) = dh_x
          temp_dy(k-1) = dh_y
          temp_dz(k-1) = dh_z
          temp_x(k-1) = well_trajectory(j,1) - dh_x/2.d0
          temp_y(k-1) = well_trajectory(j,2) - dh_y/2.d0
          temp_z(k-1) = well_trajectory(j,3) - dh_z/2.d0
          dh_x = well_trajectory(j,1)
          dh_y = well_trajectory(j,2)
          dh_z = well_trajectory(j,3)
          temp_casing(k) = well_casing(j)
          ! Don't count off-process cell ID's, but still compute
          ! where the change from on- to off- process occurs.
          if (Uninitialized(cur_id)) then
            k = k - 1
            temp_id_list(k) = temp_id_list(k+1)
          endif
          cur_id = local_id
        endif

        if (j == nsegments .and. Initialized(local_id)) then
          dh_x = well_trajectory(j,1) - dh_x
          dh_y = well_trajectory(j,2) - dh_y
          dh_z = well_trajectory(j,3) - dh_z
          temp_x(k) = well_trajectory(j,1) - dh_x/2.d0
          temp_y(k) = well_trajectory(j,2) - dh_y/2.d0
          temp_z(k) = well_trajectory(j,3) - dh_z/2.d0
          temp_dx(k) = dh_x
          temp_dy(k) = dh_y
          temp_dz(k) = dh_z
        endif
      enddo

    else
      k = 0
      do j = 1, nsegments

        dummy_h%x = segment_coordinates(j,1)
        dummy_h%y = segment_coordinates(j,2)
        dummy_h%z = segment_coordinates(j,3)

        call GridGetLocalIDFromCoordinate(res_grid,dummy_h,option,local_id)


        dh_x = segment_dxyz(j,1)
        dh_y = segment_dxyz(j,2)
        dh_z = segment_dxyz(j,3)

        if (Initialized(local_id)) then
          k = k + 1
          temp_casing(k) = well_casing(j)
          temp_dx(k) = dh_x
          temp_dy(k) = dh_y
          temp_dz(k) = dh_z
          temp_x(k) = dummy_h%x
          temp_y(k) = dummy_h%y
          temp_z(k) = dummy_h%z
          temp_id_list(k) = local_id
        endif

      enddo

    endif

    k = 0
    do j = 1, nsegments
      if (temp_z(j) > -1.d19) k = k+1
    enddo

    well_grid%nsegments = k

    ! Num local number of segments across processes,
    ! sorted well nodes by z-location in ascending order.
    call MPI_Allreduce(well_grid%nsegments,nsegments,ONE_INTEGER_MPI, &
                       MPI_INTEGER,MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
    allocate(well_nodes(nsegments))
    allocate(collect_x_temp(nsegments*option%comm%size))
    allocate(collect_y_temp(nsegments*option%comm%size))
    allocate(collect_z_temp(nsegments*option%comm%size))
    allocate(collect_rank(nsegments*option%comm%size))
    collect_x_temp = -MAX_DOUBLE
    collect_x_temp(nsegments*option%myrank+1:nsegments*option%myrank+k) = &
              temp_x(1:well_grid%nsegments)
    collect_y_temp = -MAX_DOUBLE
    collect_y_temp(nsegments*option%myrank+1:nsegments*option%myrank+k) = &
              temp_y(1:well_grid%nsegments)
    collect_z_temp = -MAX_DOUBLE
    collect_z_temp(nsegments*option%myrank+1:nsegments*option%myrank+k) = &
              temp_z(1:well_grid%nsegments)
    collect_rank = UNINITIALIZED_INTEGER
    collect_rank(nsegments*option%myrank+1:nsegments*option%myrank+k) = &
              option%myrank
    call MPI_Allreduce(MPI_IN_PLACE,collect_x_temp, &
                       nsegments*option%comm%size,MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Allreduce(MPI_IN_PLACE,collect_y_temp, &
                       nsegments*option%comm%size,MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Allreduce(MPI_IN_PLACE,collect_z_temp, &
                       nsegments*option%comm%size,MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Allreduce(MPI_IN_PLACE,collect_rank, &
                       nsegments*option%comm%size,MPI_INTEGER, &
                       MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
    ! MAN: not sure if getting the ordering exactly right for all possible
    ! cases is strictly necessary for hydrostatic well model, but should try.
    k = 0
    well_nodes(:)%id = UNINITIALIZED_INTEGER
    well_nodes(:)%x = UNINITIALIZED_DOUBLE
    well_nodes(:)%y = UNINITIALIZED_DOUBLE
    well_nodes(:)%z = UNINITIALIZED_DOUBLE
    do j = 1,nsegments*option%comm%size
      if (collect_z_temp(j) > -1.d19) then
        k = k + 1
        well_nodes(k)%id = collect_rank(j)
        well_nodes(k)%x = collect_x_temp(j)
        well_nodes(k)%y = collect_y_temp(j)
        well_nodes(k)%z = collect_z_temp(j)
      endif
    enddo
    dim = 'z'
    call UtilitySortArray(well_nodes,dim)

    ! Need to write the well in the reservoir coordinates, from bottom to top
    well_grid%nsegments = nsegments
    allocate(well_grid%casing(nsegments))
    allocate(well_grid%dh(nsegments))
    allocate(well_grid%h(nsegments))
    allocate(well_grid%h_local_id(nsegments))
    allocate(well_grid%h_ghosted_id(nsegments))
    allocate(well_grid%h_global_id(nsegments))
    allocate(well_grid%h_rank_id(nsegments))
    allocate(well_grid%strata_id(nsegments))
    allocate(well_grid%dx(nsegments))
    allocate(well_grid%dy(nsegments))
    allocate(well_grid%dz(nsegments))
    allocate(well_grid%res_z(nsegments))
    allocate(well_grid%res_dz(nsegments))

    well_grid%casing(:) = UNINITIALIZED_DOUBLE
    well_grid%dh(:) = UNINITIALIZED_DOUBLE
    well_grid%h(:)%x = -MAX_DOUBLE
    well_grid%h(:)%y = -MAX_DOUBLE
    well_grid%h(:)%z = -MAX_DOUBLE
    well_grid%h_local_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_ghosted_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_global_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_rank_id(:) = UNINITIALIZED_INTEGER
    well_grid%strata_id(:) = UNINITIALIZED_INTEGER
    well_grid%dx(:) = UNINITIALIZED_DOUBLE
    well_grid%dy(:) = UNINITIALIZED_DOUBLE
    well_grid%dz(:) = UNINITIALIZED_DOUBLE
    well_grid%res_z(:) = -MAX_DOUBLE
    well_grid%res_dz(:) = UNINITIALIZED_DOUBLE

    do k = 1,well_grid%nsegments
      if (well_nodes(k)%id == option%myrank) then
        do local_index = 1,size(temp_id_list)
          if (well_nodes(k)%x == temp_x(local_index) .and. &
              well_nodes(k)%y == temp_y(local_index) .and. &
              well_nodes(k)%z == temp_z(local_index)) then
            index = k
            well_grid%h_rank_id(index) = option%myrank
            well_grid%h_local_id(index) = temp_id_list(local_index)
            well_grid%h_ghosted_id(index) = res_grid%nL2G(well_grid% &
                                            h_local_id(index))
            ghosted_id = well_grid%h_ghosted_id(index)
            well_grid%h_global_id(index) = res_grid%nG2A(ghosted_id)
            well_grid%res_z(index) = res_grid%z(ghosted_id)
            if (temp_casing(local_index)) then
              well_grid%casing(index) = 0.d0
            else
              well_grid%casing(index) = 1.d0
            endif

            well_grid%h(index)%id = local_index
            well_grid%h(index)%x = temp_x(local_index)
            well_grid%h(index)%y = temp_y(local_index)
            well_grid%h(index)%z = temp_z(local_index)

            well_grid%dx(index) = dabs(temp_dx(local_index))
            well_grid%dy(index) = dabs(temp_dy(local_index))
            well_grid%dz(index) = dabs(temp_dz(local_index))
            well_grid%dh(index) = sqrt(well_grid%dx(index)**2 + &
                          well_grid%dy(index)**2 + well_grid%dz(index)**2)

            if (StringCompare(res_grid%ctype,'STRUCTURED')) then
              well_grid%res_dz(index) = res_grid%structured_grid%dz(index)
            elseif (StringCompare(res_grid%ctype,'IMPLICIT UNSTRUCTURED')) then
              if (associated(res_grid%unstructured_grid%face_centroid)) then
                allocate(face_centroids(size(res_grid%unstructured_grid% &
                        cell_to_face_ghosted(:,ghosted_id))))
                do j = 1,size(face_centroids)
                  face_centroids(j) = res_grid%unstructured_grid% &
                                      face_centroid(res_grid% &
                                      unstructured_grid% &
                                      cell_to_face_ghosted(j,ghosted_id))
                enddo
                dim = 'z'
                call UtilitySortArray(face_centroids,dim)
                well_grid%res_dz(index) = dabs(face_centroids(size( &
                                                face_centroids))%z - &
                                                face_centroids(1)%z)
                deallocate(face_centroids)
              endif
            else
              zmin = MAX_DOUBLE
              zmax = -MAX_DOUBLE
              do j = 1,size(res_grid%unstructured_grid%explicit_grid% &
                            face_centroids)
                if (any(res_grid%unstructured_grid%explicit_grid% &
                       connections(:,j) ==  ghosted_id)) then
                  zmin = min(zmin,res_grid%unstructured_grid%explicit_grid% &
                             face_centroids(j)%z)
                  zmax = max(zmax,res_grid%unstructured_grid%explicit_grid% &
                             face_centroids(j)%z)
                endif
              enddo
              well_grid%res_dz(index) = dabs(zmax - zmin)
            endif
            exit
          endif
        enddo
      endif
    enddo

    if (allocated(well_trajectory)) deallocate(well_trajectory)
    if (allocated(well_casing)) deallocate(well_casing)
    if (allocated(temp_id_list)) deallocate(temp_id_list)
    if (allocated(collect_x_temp)) deallocate(collect_x_temp)
    if (allocated(collect_y_temp)) deallocate(collect_y_temp)
    if (allocated(collect_z_temp)) deallocate(collect_z_temp)
    if (allocated(temp_dx)) deallocate(temp_dx)
    if (allocated(temp_dy)) deallocate(temp_dy)
    if (allocated(temp_dz)) deallocate(temp_dz)
    if (allocated(temp_x)) deallocate(temp_x)
    if (allocated(temp_y)) deallocate(temp_y)
    if (allocated(temp_z)) deallocate(temp_z)
    if (allocated(temp_casing)) deallocate(temp_casing)
    if (allocated(well_nodes)) deallocate(well_nodes)
    if (allocated(collect_rank)) deallocate(collect_rank)

  else
  top_of_reservoir = res_grid%z_max_global
  top_of_hole = well_grid%tophole(3)
  bottom_of_reservoir = res_grid%z_min_global
  bottom_of_hole = well_grid%bottomhole(3)
  if (top_of_reservoir < top_of_hole) then
    ! MAN: this might be hard to get right, not sure it is necessary.
    ! option%io_buffer = 'The WELLBORE_MODEL TOP_OF_HOLE coordinates extend &
    !                    &beyond the top of the reservoir domain. &
    !                    &You must fix the TOP_OF_HOLE coordinates to align &
    !                    &with the top face of the reservoir grid cell that &
    !                    &it occupies.'
    ! call PrintErrMsg(option)
    well_grid%tophole(3) = top_of_reservoir
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
  ! MAN: this might be hard to get right, not sure it is necessary.
  ! option%io_buffer = 'The WELLBORE_MODEL BOTTOM_OF_HOLE coordinates extend &
  !                    &beyond the bottom of the reservoir domain. &
  !                    &You must fix the BOTTOM_OF_HOLE coordinates so that &
  !                    &the bottom of the well is aligned with the bottom &
  !                    &face of the reservoir, or is above the bottom &
  !                    &face of the reservoir in the vertical column that the &
  !                    &well occupies.'
  ! call PrintErrMsg(option)
    well_grid%bottomhole(3) = bottom_of_reservoir
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
    temp_id_list = UNINITIALIZED_INTEGER
    temp_repeated_list = UNINITIALIZED_INTEGER
    dummy_h%x = well_grid%bottomhole(1)
    dummy_h%y = well_grid%bottomhole(2)

    k = 0
    j = 0
    cur_id = UNINITIALIZED_INTEGER
    repeated = 0
    cum_z = 0
    cur_cum_z_int = 0
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
    allocate(well_grid%res_z(nsegments))

    well_grid%dh(:) = UNINITIALIZED_DOUBLE
    well_grid%res_dz(:) = UNINITIALIZED_DOUBLE
    well_grid%h(:)%x = -MAX_DOUBLE
    well_grid%h(:)%y = -MAX_DOUBLE
    well_grid%h(:)%z = -MAX_DOUBLE
    well_grid%h_local_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_ghosted_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_global_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_rank_id(:) = UNINITIALIZED_INTEGER
    well_grid%strata_id(:) = UNINITIALIZED_INTEGER
    well_grid%res_z(:) = -MAX_DOUBLE

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
      well_grid%res_z(k) = res_grid%z(well_grid%h_ghosted_id(k))
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
                         &contains the same number of values (i.e. one &
                         &z-center value for every length value).'
      call PrintErrMsg(option)
    endif
    if (size(well_grid%l_list) /= size(well_grid%z_list)) then
      option%io_buffer = 'The length of SEGMENT_LENGTH_VALUES must match the &
                         &length of SEGMENT_CENTER_Z_VALUES (i.e. one z-center&
                         & value for every length value) provided in the &
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
    allocate(well_grid%res_z(nsegments))

    well_grid%dh(:) = UNINITIALIZED_DOUBLE
    well_grid%res_dz(:) = UNINITIALIZED_DOUBLE
    well_grid%h(:)%x = -MAX_DOUBLE
    well_grid%h(:)%y = -MAX_DOUBLE
    well_grid%h(:)%z = -MAX_DOUBLE
    well_grid%h_local_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_ghosted_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_global_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_rank_id(:) = UNINITIALIZED_INTEGER
    well_grid%strata_id(:) = UNINITIALIZED_INTEGER
    well_grid%res_z(:) = -MAX_DOUBLE

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
        well_grid%res_z(k) = res_grid%z(well_grid%h_ghosted_id(k))
      endif
    enddo

    well_grid%dh = well_grid%l_list

    temp_real = SUM(well_grid%dh)
    total_length = sqrt((diff_x*diff_x)+(diff_y*diff_y)+(diff_z*diff_z))
    if (temp_real /= total_length) then
      temp_real2 = abs(temp_real - total_length)
      !write(*,*) temp_real2
      if (temp_real2 > 1.d-2) then
        option%io_buffer = 'The sum of the list of SEGMENT_LENGTH_VALUES &
          &(' // StringWrite(temp_real) // ' m) does not match the total &
          &length of the well according to the coordinates provided by &
          &WELLBORE_MODEL, TOP_OF_HOLE and WELLBORE_MODEL,TOP_OF_HOLE &
          &(' // StringWrite(total_length) // ' m).'
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
    allocate(well_grid%res_z(nsegments))

    well_grid%dh(:) = UNINITIALIZED_DOUBLE
    well_grid%res_dz(:) = UNINITIALIZED_DOUBLE
    well_grid%h(:)%x = -MAX_DOUBLE
    well_grid%h(:)%y = -MAX_DOUBLE
    well_grid%h(:)%z = -MAX_DOUBLE
    well_grid%h_local_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_ghosted_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_global_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_rank_id(:) = UNINITIALIZED_INTEGER
    well_grid%strata_id(:) = UNINITIALIZED_INTEGER
    well_grid%res_z(:) = -MAX_DOUBLE



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
        well_grid%res_z(k) = res_grid%z(well_grid%h_ghosted_id(k))
      endif
    enddo
  !---------------------------------------------------------------------------
  else
   ! Use reservoir grid info
    allocate(dz_list(num_entries))
    allocate(res_dz_list(num_entries))
    allocate(cell_id_list(num_entries))
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
     allocate(well_grid%res_z(nsegments))

    well_grid%dh(:) = UNINITIALIZED_DOUBLE
    well_grid%res_dz(:) = UNINITIALIZED_DOUBLE
    well_grid%h(:)%x = -MAX_DOUBLE
    well_grid%h(:)%y = -MAX_DOUBLE
    well_grid%h(:)%z = -MAX_DOUBLE
    well_grid%h_local_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_ghosted_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_global_id(:) = UNINITIALIZED_INTEGER
    well_grid%h_rank_id(:) = UNINITIALIZED_INTEGER
    well_grid%strata_id(:) = UNINITIALIZED_INTEGER
    well_grid%res_z(:) = -MAX_DOUBLE

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
     well_grid%res_z(1) = res_grid%z(well_grid%h_global_id(1))

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
       well_grid%res_z(k) = res_grid%z(well_grid%h_ghosted_id(k))
     enddo
  !---------------------------------------------------------------------------
  endif
  endif

  allocate(well_grid%strata_id(nsegments))
  well_grid%strata_id(:) = UNINITIALIZED_INTEGER
  well_grid%nconnections = well_grid%nsegments - 1

end subroutine PMWellSetupGrid

! ************************************************************************** !

subroutine PMWellSetupBase(this)
  !
  ! Initializes variables associated with the base well process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_type) :: this
  this%option%io_buffer = 'PMWellSetupBase must be extended.'
  call PrintErrMsg(this%option)

end subroutine PMWellSetupBase

! ************************************************************************** !

subroutine PMWellSetupSequential(this)
  !
  ! Initializes variables associated with the base well process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_sequential_type) :: this

  this%option%io_buffer = 'PMWellSetupSequential must be extended.'
  call PrintErrMsg(this%option)

end subroutine PMWellSetupSequential

! ************************************************************************** !

subroutine PMWellSetupImplicit(this)
  !
  ! Initializes variables associated with the base well process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_implicit_type) :: this

  this%option%io_buffer = 'PMWellSetupImplicit must be extended.'
  call PrintErrMsg(this%option)

end subroutine PMWellSetupImplicit

! ************************************************************************** !

subroutine PMWellSetupQI(this)
  !
  ! Initializes variables associated with the base well process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_qi_type) :: this
  this%option%io_buffer = 'PMWellSetupQI must be extended.'
  call PrintErrMsg(this%option)

end subroutine PMWellSetupQI

! ************************************************************************** !

subroutine PMWellSetupHydrostatic(this)
  !
  ! Initializes variables associated with the hydrostatic well process model.
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  use Grid_module
  use Utility_module
  use String_module
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
  use Hydrate_Aux_module
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use Matrix_Zeroing_module

  implicit none

  class(pm_well_hydrostatic_type) :: this

  type(option_type), pointer :: option
  type(grid_type), pointer :: res_grid
  type(well_grid_type), pointer :: well_grid
  type(coupler_type), pointer :: source_sink, coupler
  type(input_type) :: input_dummy
  class(realization_subsurface_type), pointer :: realization
  type(point3d_type) :: dummy_h
  character(len=MAXSTRINGLENGTH) :: string, filename
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
  PetscInt :: k, j, i, index
  PetscInt :: count1, count2_local, count2_global
  PetscBool :: well_grid_res_is_OK = PETSC_FALSE
  PetscBool :: res_grid_cell_within_well_z
  PetscBool :: res_grid_cell_within_well_y
  PetscBool :: res_grid_cell_within_well_x
  PetscErrorCode :: ierr
  PetscInt :: well_bottom_local, well_bottom_ghosted
  PetscInt :: iconn, sum_connection, global_ss_connections
  PetscInt :: num_new_source_sinks, offset
  PetscInt, allocatable :: temp(:), temp2(:)
  PetscInt :: fid = 86
  type(global_auxvar_type), pointer :: auxvars_ss(:), auxvars_ss_temp(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_ss(:), &
                                                   rt_auxvars_ss_temp(:)
  type(matrix_zeroing_type), pointer :: matrix_zeroing

  call this%SetRealization()
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

  option%io_buffer = 'WELLBORE_MODEL: Grid created with ' // &
                      StringWrite(well_grid%nsegments)// &
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
  call MPI_Allreduce(MPI_IN_PLACE,well_grid%casing,size(well_grid%casing), &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                     CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE,well_grid%dh,nsegments, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                     CHKERRQ(ierr)
  if (associated(well_grid%dx)) then
    call MPI_Allreduce(MPI_IN_PLACE,well_grid%dx,nsegments, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                      CHKERRQ(ierr)
  endif
  if (associated(well_grid%dy)) then
    call MPI_Allreduce(MPI_IN_PLACE,well_grid%dy,nsegments, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                      CHKERRQ(ierr)
  endif
  if (associated(well_grid%dz)) then
    call MPI_Allreduce(MPI_IN_PLACE,well_grid%dz,nsegments, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                      CHKERRQ(ierr)
  endif
  call MPI_Allreduce(MPI_IN_PLACE,well_grid%res_dz,nsegments, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                     CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE,well_grid%res_z,nsegments, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                     CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE,well_grid%h(:)%x,nsegments, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                     CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE,well_grid%h(:)%y,nsegments, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                     CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE,well_grid%h(:)%z,nsegments, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                     CHKERRQ(ierr)
  well_grid%h_rank_id = h_all_rank_id
  well_grid%h_global_id = h_all_global_id

  do k = 1,res_grid%ngmax
    do j = 1,well_grid%nsegments
      if (res_grid%nG2A(k) == well_grid%h_global_id(j)) then
        well_grid%h_ghosted_id(j) = k
      endif
    enddo
  enddo

  h_rank_id_unique(:) = UNINITIALIZED_INTEGER

  min_val = minval(h_all_rank_id)-1
  max_val = maxval(h_all_rank_id)

  k = 0
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
  if (this%well%WI_model == PEACEMAN_NONE) then
    if (.not. associated(this%well%diameter)) then
      allocate(this%well%diameter(nsegments))
    endif
    if (.not. associated(this%well%friction_factor)) then
      allocate(this%well%friction_factor(nsegments))
    endif
    if (.not. associated(this%well%skin)) then
      allocate(this%well%skin(nsegments))
    endif
    this%well%diameter = 0.d0
    this%well%friction_factor = 0.d0
    this%well%skin = 0.d0
  endif

  allocate(this%well%output_in_well_file(nsegments))
  ! by default, print all segments in the .well file:
  this%well%output_in_well_file = PETSC_TRUE
  if (associated(this%well%segments_for_output)) then
    this%well%output_in_well_file = PETSC_FALSE
    do k = 1,size(this%well%segments_for_output)
      if (this%well%segments_for_output(k) > nsegments) then
        option%io_buffer = 'Invalid segment number encountered. Please review &
          &the list of segment numbers provided after the PRINT_WELL_FILE &
          &SEGMENTS keyword. At least one segment number exceeds the total &
          &number of well segments, as defined in the WELL_GRID block.'
        call PrintErrMsg(option)
      endif
      this%well%output_in_well_file(this%well%segments_for_output(k)) = &
        PETSC_TRUE
    enddo
  endif

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

  if (size(this%well_grid%casing) /= nsegments) then
    if (size(this%well_grid%casing) == 1) then
      temp_real = this%well_grid%casing(1)
      deallocate(this%well_grid%casing)
      allocate(this%well_grid%casing(nsegments))
      this%well_grid%casing(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,GRID,CASING must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,GRID,CASING, &
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

  if (res_grid%itype == STRUCTURED_GRID) then
    ! MAN: I think this check does not work for explicit unstructured grids.

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
      call MPI_Allreduce(count2_local,count2_global,ONE_INTEGER_MPI, &
                        MPI_INTEGER, MPI_SUM,this%well_comm%comm,ierr); &
                        CHKERRQ(ierr)
    endif
    call MPI_Bcast(count2_global,ONE_INTEGER_MPI,MPI_INTEGER, &
                  this%well_comm%petsc_rank_list(1),option%mycomm, &
                  ierr);CHKERRQ(ierr)

    ! The only way we can ensure that the well discretization did not skip a
    ! reservoir cell, is if the number of unique global_id's that the well
    ! is connected to (count1) matches the number of reservoir grid cells that
    ! the well occupies (count2):
    if (count1 == count2_global) then
      well_grid_res_is_OK = PETSC_TRUE
    elseif (associated(well_grid%deviated_well_segment_list)) then
      !Basically bypass this check for now, since well grid was generated from
      !res grid
      well_grid_res_is_OK = PETSC_TRUE
    endif

    do k = 1,nsegments
      if (Uninitialized(well_grid%h_local_id(k))) well_grid%h_local_id(k) = -1
    enddo


    if (.not.well_grid_res_is_OK) then
      option%io_buffer = 'ERROR:  &
        &The number of reservoir grid cells that are occupied by the well &
        &(' // StringWrite(count2_global) // ') is larger than the number of &
        &unique reservoir grid cells that have a connection to the well (' // &
        StringWrite(count1) // '). Therefore, some of the reservoir grid cells &
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
  endif

  if (this%well_grid%save_well_segment_list .and. &
      associated(well_grid%deviated_well_segment_list)) then
110 format(100es24.16)
    filename = trim(this%name) // '_well-segments.dat'
    if (OptionIsIORank(option)) then
      open(unit=fid,file=filename,action="write",status="replace")
      do i = 1,this%well_grid%nsegments
        index = this%well_grid%nsegments - i + 1
        if (this%well_grid%casing(index) < 5.d-1) then
          string = 'CASED'
        else
          string = 'SCREENED'
        endif
        write(fid,'(a)',advance="no") trim(string)
        write(fid,110,advance="no") this%well_grid%h(index)%x
        write(fid,110,advance="no") this%well_grid%h(index)%y
        write(fid,110,advance="no") this%well_grid%h(index)%z
        write(fid,110,advance="no") this%well_grid%dx(index)
        write(fid,110,advance="no") this%well_grid%dy(index)
        write(fid,110,advance="yes") this%well_grid%dz(index)
      enddo
      close(fid)
    endif
  endif

  ! add a reservoir src/sink coupler for each well segment
  num_new_source_sinks = 0
  if (associated(realization%patch%aux%Global%auxvars_ss)) then
    offset = realization%patch%aux%Global%num_aux_ss
  else
    offset = 0
  endif
  do k = 1,well_grid%nsegments

    source_sink => CouplerCreate()
    source_sink%itype = SRC_SINK_COUPLER_TYPE
    source_sink%name = trim(this%name) // '_well_segment_' // StringWrite(k)

    ! ----- flow ------------------
    source_sink%flow_condition_name = trim(this%name) // '_well_segment_' // &
                                      StringWrite(k) // '_flow_srcsink'
    source_sink%flow_condition => FlowConditionCreate(option)
    source_sink%flow_condition%name = source_sink%flow_condition_name
    select case (option%iflowmode)
      case(SCO2_MODE)
        source_sink%flow_condition%sco2 => FlowSCO2ConditionCreate(option)
        string = 'RATE'
        source_sink%flow_condition%sco2%rate => FlowSCO2SubConditionPtr( &
          input_dummy,string,source_sink%flow_condition%sco2,option)
        source_sink%flow_condition%sco2%rate%itype = SCALED_MASS_RATE_SS
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

        if (this%use_well_coupler) then
          coupler => realization%patch%well_coupler_list%first
          do
            if (.not.associated(coupler)) exit
            if (StringCompare(coupler%well_name,this%name)) then
              source_sink%well_name = coupler%well_name
              if (option%ntrandof > 0) then
                if (option%itranmode == RT_MODE) then
                  source_sink%tran_condition => coupler%tran_condition
                endif
              endif
            endif
            coupler => coupler%next
          enddo
        endif

      case(H_MODE)
        source_sink%flow_condition%hydrate => FlowHydrateConditionCreate(option)
        string = 'RATE'
        source_sink%flow_condition%hydrate%rate => FlowHydrateSubConditionPtr( &
          input_dummy,string,source_sink%flow_condition%hydrate,option)
        source_sink%flow_condition%hydrate%rate%itype = SCALED_MASS_RATE_SS
        source_sink%flow_condition%hydrate%liquid_pressure => &
              FlowHydrateSubConditionPtr(input_dummy,string,source_sink% &
                                        flow_condition%hydrate,option)
        source_sink%flow_condition%hydrate%gas_pressure => &
              FlowHydrateSubConditionPtr(input_dummy,string,source_sink% &
                                        flow_condition%hydrate,option)
        allocate(source_sink%flow_condition%hydrate%rate%dataset%rarray(2))
        source_sink%flow_condition%hydrate%rate%dataset%rarray(:) = 0.d0
        source_sink%flow_condition%well => FlowSubConditionCreate(ONE_INTEGER)
        ! Bottom of hole is special for fully implicit coupling with steady
        ! state well model.
        if (k==1) then
          source_sink%flow_condition%well%ctype = 'well-bottom'
        else
          source_sink%flow_condition%well%ctype = 'well'
        endif

        if (this%use_well_coupler) then
          coupler => realization%patch%well_coupler_list%first
          do
            if (.not.associated(coupler)) exit
            if (StringCompare(coupler%well_name,this%name)) then
              source_sink%well_name = coupler%well_name
              if (option%ntrandof > 0) then
                if (option%itranmode == RT_MODE) then
                  source_sink%tran_condition => coupler%tran_condition
                endif
              endif
            endif
            coupler => coupler%next
          enddo
        endif

    end select

    source_sink%connection_set => &
      ConnectionCreate(1,GENERIC_CONNECTION_TYPE,res_grid%itype)
    source_sink%connection_set%id_dn = well_grid%h_local_id(k)
    if (well_grid%h_local_id(k) < 0) then
      source_sink%connection_set%num_connections = 0
    else
      num_new_source_sinks = num_new_source_sinks + 1
      source_sink%connection_set%offset = num_new_source_sinks + offset - 1
    endif

    call CouplerAddToList(source_sink,this%realization%patch%source_sink_list)
    nullify(source_sink)
  enddo

  if (this%use_well_coupler) then
    coupler => realization%patch%well_coupler_list%first
    do
      if (.not.associated(coupler)) exit
      if (StringCompare(coupler%well_name,this%name)) then
        this%flow_condition => coupler%flow_condition
      endif
      coupler => coupler%next
    enddo
    if (.not. associated(this%flow_condition)) then
      option%io_buffer = 'Well model invokes USE_WELL_COUPLER &
        &for flow conditions, but no WELL_COUPLERs were found &
        &linked to well named '// trim(this%name) // '.'
      call PrintErrMsg(option)
    endif
  endif

  ! Enable allocation of mass balance array
  if (associated(realization%patch%aux%Global%auxvars_ss)) then
    auxvars_ss => realization%patch%aux%Global%auxvars_ss
    global_ss_connections = realization%patch%aux%Global%num_aux_ss
    sum_connection = num_new_source_sinks + global_ss_connections
    allocate(auxvars_ss_temp(sum_connection))
    auxvars_ss_temp(1:global_ss_connections) = auxvars_ss(:)
    realization%patch%aux%Global%auxvars_ss => auxvars_ss_temp
    deallocate(auxvars_ss)
    nullify(auxvars_ss_temp)
    auxvars_ss => realization%patch%aux%Global%auxvars_ss
    if (option%ntrandof > 0) then
      rt_auxvars_ss => realization%patch%aux%RT%auxvars_ss
      allocate(rt_auxvars_ss_temp(sum_connection))
      rt_auxvars_ss_temp(1:global_ss_connections) = rt_auxvars_ss(:)
      realization%patch%aux%RT%auxvars_ss => rt_auxvars_ss_temp
      deallocate(rt_auxvars_ss)
      nullify(rt_auxvars_ss_temp)
      realization%patch%aux%RT%num_aux_ss = &
            realization%patch%aux%RT%num_aux_ss + num_new_source_sinks
    endif
  else
    sum_connection = num_new_source_sinks
    global_ss_connections = 0
    allocate (auxvars_ss(sum_connection))
    realization%patch%aux%Global%auxvars_ss => auxvars_ss
    if (option%ntrandof > 0) then
      allocate(realization%patch%aux%RT%auxvars_ss(sum_connection))
      realization%patch%aux%RT%num_aux_ss = &
            realization%patch%aux%RT%num_aux_ss + num_new_source_sinks
    endif
  endif
  select case (option%iflowmode)
  case(SCO2_MODE)
    realization%patch%aux%SCO2%num_aux_ss = realization%patch%aux%SCO2% &
                                            num_aux_ss + num_new_source_sinks
  case(H_MODE)
    realization%patch%aux%Hydrate%num_aux_ss = realization%patch%aux% &
                                    Hydrate%num_aux_ss + num_new_source_sinks
  end select
  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array
    do iconn = global_ss_connections + 1, sum_connection
      call GlobalAuxVarInit(auxvars_ss(iconn),option)
      if (option%ntrandof > 0) then
        call RTAuxVarInit(realization%patch%aux%RT%auxvars_ss(iconn), &
                          realization%reaction,PETSC_FALSE,option)
      endif
    enddo
  endif
  realization%patch%aux%Global%num_aux_ss = sum_connection

  if (sum_connection > 0) then
    ! flow
    if (option%nflowdof > 0) then
      allocate(realization%patch%ss_flow_fluxes(option%nflowdof, &
                sum_connection))
      realization%patch%ss_flow_fluxes = 0.d0
    endif
    ! transport
    if (option%ntrandof > 0) then
      allocate(realization%patch%ss_tran_fluxes(option%ntrandof, &
                sum_connection))
      realization%patch%ss_tran_fluxes = 0.d0
      ! only needed by transport
      allocate(realization%patch%ss_flow_vol_fluxes(option%nphase, &
                sum_connection))
      realization%patch%ss_flow_vol_fluxes = 0.d0
    endif
  endif

  ! For fully-implicit well coupling, resize the matrix zeroing arrays to
  ! exclude the bottom of the well for hydrostatic well model.

  well_bottom_ghosted = well_grid%h_ghosted_id(1)
  well_bottom_local = well_grid%h_local_id(1)
  nullify(matrix_zeroing)
  select case (option%iflowmode)
    case(SCO2_MODE)
      matrix_zeroing => realization%patch%aux%SCO2%matrix_zeroing
    case(H_MODE)
      matrix_zeroing => realization%patch%aux%Hydrate%matrix_zeroing
  end select
  if (associated(matrix_zeroing)) then
    if (well_bottom_local > 0 .and. matrix_zeroing%zero_rows_exist) then
      if (size(matrix_zeroing%zero_rows_local) == 1) then
        deallocate(matrix_zeroing%zero_rows_local)
        deallocate(matrix_zeroing%zero_rows_local_ghosted)
        matrix_zeroing%zero_rows_exist = PETSC_FALSE
        matrix_zeroing%n_zero_rows = 0
      else
        if (allocated(temp)) deallocate(temp)
        if (allocated(temp2)) deallocate(temp2)
        allocate(temp(size(matrix_zeroing%zero_rows_local)))
        allocate(temp2(size(matrix_zeroing%zero_rows_local)-1))
        matrix_zeroing%n_zero_rows = matrix_zeroing%n_zero_rows - 1
        temp(:) = matrix_zeroing%zero_rows_local(:)
        temp2 = pack(temp,temp /= well_bottom_local*option%nflowdof)
        deallocate(matrix_zeroing%zero_rows_local)
        allocate(matrix_zeroing%zero_rows_local(size(temp2)))
        matrix_zeroing%zero_rows_local(:) = temp2(:)
        temp(:) = matrix_zeroing%zero_rows_local_ghosted(:)
        temp2 = pack(temp,temp /= well_bottom_ghosted*option%nflowdof-1)
        deallocate(matrix_zeroing%zero_rows_local_ghosted)
        allocate(matrix_zeroing%zero_rows_local_ghosted(size(temp2)))
        matrix_zeroing%zero_rows_local_ghosted(:) = temp2(:)
      endif
    endif
  endif

end subroutine PMWellSetupHydrostatic

! ************************************************************************** !

subroutine PMWellReadSimOptionsBlockBase(this,input)
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

  input%ierr = INPUT_ERROR_NONE
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
      case default
          call InputKeywordUnrecognized(input,keyword,'Well Model',option)
    end select

  enddo
  call InputPopBlock(input,option)

  if (Uninitialized(this%well%well_model_type)) then
    option%io_buffer = 'WELL_MODEL_TYPE must be specified in the &
                        &WELL_MODEL OPTIONS block.'
    call PrintErrMsg(option)
  endif

end subroutine PMWellReadSimOptionsBlockBase

! ************************************************************************** !

subroutine PMWellReadSimOptionsBlockSeq(this,input)
  !
  ! Author: Michael Nole
  ! Date: 03/08/24
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  class(pm_well_sequential_type) :: this

  call PMWellReadSimOptionsBlockBase(this,input)

end subroutine PMWellReadSimOptionsBlockSeq

! ************************************************************************** !

subroutine PMWellReadSimOptionsBlockImp(this,input)
  !
  ! Author: Michael Nole
  ! Date: 03/08/24
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  class(pm_well_implicit_type) :: this

  call PMWellReadSimOptionsBlockBase(this,input)

end subroutine PMWellReadSimOptionsBlockImp

! ************************************************************************** !

subroutine PMWellReadSimOptionsBlockQI(this,input)
  !
  ! Author: Michael Nole
  ! Date: 03/08/24
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  class(pm_well_qi_type) :: this

  call PMWellReadSimOptionsBlockBase(this,input)

end subroutine PMWellReadSimOptionsBlockQI

! ************************************************************************** !

subroutine PMWellReadSimOptionsBlockHyd(this,input)
  !
  ! Author: Michael Nole
  ! Date: 03/08/24
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  class(pm_well_hydrostatic_type) :: this

  call PMWellReadSimOptionsBlockBase(this,input)

end subroutine PMWellReadSimOptionsBlockHyd

! ************************************************************************** !

subroutine PMWellReadPMBlockBase(this,input)
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
  PetscInt, pointer :: temp_seg_nums(:)
  PetscBool :: found
  PetscInt :: k, num_read
  PetscInt :: read_max = 10000

  option => this%option
  well_grid => this%well_grid
  input%ierr = INPUT_ERROR_NONE
  error_string = 'WELLBORE_MODEL'

  option%io_buffer = 'pflotran card:: WELLBORE_MODEL'
  call PrintMsg(option)

  allocate(temp_seg_nums(read_max))

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
      case('PRINT_WELL_FILE')
        this%print_well = PETSC_TRUE
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) cycle
        call StringToUpper(word)
        if (StringCompare(word,'SEGMENTS')) then
          ! count the segment numbers
          num_read = 0
          do k = 1,read_max
            call InputReadInt(input,option,temp_seg_nums(k))
            if (InputError(input)) exit
            if (temp_seg_nums(k) <= 0) then
              option%io_buffer = 'A value provided for SEGMENTS &
                &after the ' // trim(error_string) // ', PRINT_WELL_FILE &
                &keyword was 0 or negative. Only positive integers are valid.'
              call PrintErrMsg(option)
            endif
            num_read = num_read + 1
          enddo
          if (num_read == 0) then
            option%io_buffer = 'At least one value for SEGMENTS &
              &must be provided after the ' // trim(error_string) // ', &
              &PRINT_WELL_FILE keyword.'
            call PrintErrMsg(option)
          endif
        else
          option%io_buffer = 'Keyword ' // trim(word) // ' not recognized &
            &after the ' // trim(error_string) // ', PRINT_WELL_FILE keyword&
            &. Did you mean "SEGMENTS"?'
          call PrintErrMsg(option)
        endif
        allocate(this%well%segments_for_output(num_read))
        this%well%segments_for_output(1:num_read) = temp_seg_nums(1:num_read)
        cycle
    !-------------------------------------
      case('SPLIT_WELL_FILE')
        this%split_output_file = PETSC_TRUE
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
    call PMWellReadWellOutput(this,input,option,word,error_string,found)
    if (found) cycle

    if (.not. found) then
      option%io_buffer = 'Keyword "' // trim(word) // &
                         '" does not exist for WELLBORE_MODEL.'
      call PrintErrMsg(option)
    endif

  enddo

  deallocate(temp_seg_nums)

end subroutine PMWellReadPMBlockBase

! ************************************************************************** !

subroutine PMWellReadPMBlockSeq(this,input)
  !
  ! Author: Michael Nole
  ! Date: 03/08/24
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  class(pm_well_sequential_type) :: this

  call PMWellReadPMBlockBase(this,input)

end subroutine PMWellReadPMBlockSeq

! ************************************************************************** !

subroutine PMWellReadPMBlockImplicit(this,input)
  !
  ! Author: Michael Nole
  ! Date: 03/08/24
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  class(pm_well_implicit_type) :: this

  call PMWellReadPMBlockBase(this,input)

end subroutine PMWellReadPMBlockImplicit

! ************************************************************************** !

subroutine PMWellReadPMBlockQI(this,input)
  !
  ! Author: Michael Nole
  ! Date: 03/08/24
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  class(pm_well_qi_type) :: this

  call PMWellReadPMBlockBase(this,input)

end subroutine PMWellReadPMBlockQI

! ************************************************************************** !

subroutine PMWellReadPMBlockHydrostatic(this,input)
  !
  ! Author: Michael Nole
  ! Date: 03/08/24
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  class(pm_well_hydrostatic_type) :: this
  type(option_type), pointer :: option

  option => this%option

  call InputPushBlock(input,option)
  call PMWellReadPMBlockBase(this,input)
  call InputPopBlock(input,option)

end subroutine PMWellReadPMBlockHydrostatic

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
  use Utility_module

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
  PetscInt :: at_index, nsegments, length
  PetscReal :: index_val
  PetscReal, pointer :: temp_well_index(:)
  character(len=MAXWORDLENGTH) :: word, index_word
  PetscReal, pointer :: temp_z_list(:)
  PetscReal, pointer :: temp_l_list(:)
  type(deviated_well_type), pointer :: deviated_well_segment
  type(point3d_type) :: temp_coordinate
  PetscReal, allocatable :: temp_segment_list(:,:)
  PetscBool, allocatable :: temp_casing(:)
  type(input_type), pointer :: input2
  character(len=MAXSTRINGLENGTH) :: filename
  PetscBool :: data_file_read
  PetscReal :: epsilon = 1.d-7

  allocate(temp_well_index(read_max))
  temp_well_index(:) = UNINITIALIZED_DOUBLE

  error_string = trim(error_string) // ',WELL_GRID'
  found = PETSC_TRUE
  data_file_read = PETSC_FALSE
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
            call InputErrorMsg(input,option,'XY_SEARCH_MULTIPLIER', &
                               error_string)
        !-----------------------------
          case('CASING')
            num_read = 0
            do k = 1,read_max
              index_word = trim(input%buf)
              call InputReadDouble(input,option,temp_well_index(k))
              if (InputError(input)) then
                at_index = 0; at_index = index(index_word,'@')
                if (at_index > 0) then
                  read(index_word(1:at_index-1), *) nsegments
                  read(index_word(at_index+1:len(index_word)), *) index_val
                  temp_well_index(k:k+nsegments-1) = index_val
                  num_read = num_read + nsegments
                endif
                exit
              endif
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for CASING &
                &must be provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(well_grid%casing(num_read))
            well_grid%casing(1:num_read) = 1.d0 - temp_well_index(1:num_read)
          !-----------------------------
          case('WELL_TRAJECTORY')
            call InputPushBlock(input,option)
            do
              call InputReadPFLOTRANString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'keyword', &
                                 'WELL_TRAJECTORY')
              call StringToUpper(word)

              select case(trim(word))
                case('SURFACE_ORIGIN')
                  if (associated(well_grid%deviated_well_segment_list)) then
                    option%io_buffer = "Error in constructing &
                      &WELL_TRAJECTORY: please only specify &
                      &one SURFACE_ORIGIN."
                    call PrintErrMsg(option)
                  endif
                  allocate(well_grid%deviated_well_segment_list)
                  deviated_well_segment => well_grid%deviated_well_segment_list
                  call WellSegmentInit(deviated_well_segment)
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       surface_origin(1))
                  call InputErrorMsg(input,option, &
                                     'SURFACE_ORIGIN x-coordinate', &
                                     error_string)
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       surface_origin(2))
                  call InputErrorMsg(input,option, &
                                     'SURFACE_ORIGIN y-coordinate', &
                                     error_string)
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       surface_origin(3))
                  call InputErrorMsg(input,option, &
                                     'SURFACE_ORIGIN z-coordinate', &
                                     error_string)
                case('SEGMENT_LIST')
                  if (.not. associated(deviated_well_segment)) then
                    option%io_buffer = "Error in constructing &
                      &WELL_TRAJECTORY: please first specify &
                      &SURFACE_ORIGIN before SEGMENT_LIST."
                    call PrintErrMsg(option)
                  endif
                  allocate(deviated_well_segment%next)
                  deviated_well_segment => deviated_well_segment%next
                  call WellSegmentInit(deviated_well_segment)
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call StringToUpper(word)
                  length = len_trim(word)
                  if (length > 0) then
                    if (StringCompare(word,'FILE',FOUR_INTEGER)) then
                      call InputReadFilename(input,option,filename)
                      input2 => InputCreate(IUNIT_TEMP,filename,option)
                      data_file_read = PETSC_TRUE
                      call InputPushBlock(input2,option)
                    endif
                  else
                    call InputPushBlock(input,option)
                    input2 => input
                  endif
                  nsegments = 0
                  do
                    call InputReadPFLOTRANString(input2,option)
                    if (InputError(input2)) exit
                    if (InputCheckExit(input2,option)) exit
                    nsegments = nsegments + 1
                    call InputReadWord(input2,option,word,PETSC_TRUE)
                    allocate(temp_casing(nsegments))
                    allocate(temp_segment_list(nsegments,3))
                    if (allocated(deviated_well_segment &
                                  %segment_coordinates)) then
                      temp_segment_list(1:nsegments-1,:) = &
                                  deviated_well_segment% &
                                  segment_coordinates(:,:)
                      deallocate(deviated_well_segment%segment_coordinates)
                    endif
                    if (allocated(deviated_well_segment%casing)) then
                      temp_casing(1:nsegments-1) = &
                                      deviated_well_segment%casing(:)
                      deallocate(deviated_well_segment%casing)
                    endif
                    select case(word)
                      case('CASED')
                        temp_casing(nsegments) = PETSC_TRUE
                      case('UNCASED','SCREENED','PERFORATED')
                        temp_casing(nsegments) = PETSC_FALSE
                      case default
                        option%io_buffer = "Error in constructing &
                            &WELL_TRAJECTORY: please specify &
                            &whether each segment is CASED or SCREENED/UNCASED."
                        call PrintErrMsg(option)
                    end select
                    allocate(deviated_well_segment%casing(nsegments))
                    deviated_well_segment%casing(:) = temp_casing(:)
                    deallocate(temp_casing)
                    call GeometryReadCoordinate(input2,option, &
                                              temp_coordinate, &
                                              'WELL_TRAJECTORY: SEGMENT_LIST')
                    temp_segment_list(nsegments,1) = temp_coordinate%x
                    temp_segment_list(nsegments,2) = temp_coordinate%y
                    temp_segment_list(nsegments,3) = temp_coordinate%z
                    allocate(deviated_well_segment% &
                             segment_coordinates(nsegments,3))
                    deviated_well_segment%segment_coordinates(:,:) = &
                                      temp_segment_list(:,:)
                    deallocate(temp_segment_list)
                    allocate(temp_segment_list(nsegments,3))
                    if (allocated(deviated_well_segment%segment_dxyz)) then
                      temp_segment_list(1:nsegments-1,:) = &
                                      deviated_well_segment%segment_dxyz(:,:)
                      deallocate(deviated_well_segment%segment_dxyz)
                    endif
                    call GeometryReadCoordinate(input2,option, &
                                               temp_coordinate, &
                                               'WELL_TRAJECTORY: SEGMENT_LIST')
                    temp_segment_list(nsegments,1) = temp_coordinate%x
                    temp_segment_list(nsegments,2) = temp_coordinate%y
                    temp_segment_list(nsegments,3) = temp_coordinate%z
                    allocate(deviated_well_segment%segment_dxyz(nsegments,3))
                    deviated_well_segment%segment_dxyz(:,:) = &
                                      temp_segment_list(:,:)
                    deallocate(temp_segment_list)
                  enddo
                  if (data_file_read) then
                    call InputPopBlock(input2,option)
                    call InputDestroy(input2)
                  else
                    call InputPopBlock(input2,option)
                    nullify(input2)
                  endif
                case('SEGMENT_DXYZ')
                  if (.not. associated(deviated_well_segment)) then
                    option%io_buffer = "Error in constructing &
                      &WELL_TRAJECTORY: please first specify &
                      &SURFACE_ORIGIN before SEGMENT_DXYZ."
                    call PrintErrMsg(option)
                  endif
                  allocate(deviated_well_segment%next)
                  deviated_well_segment => deviated_well_segment%next
                  call WellSegmentInit(deviated_well_segment)
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  select case(word)
                    case('CASED')
                      deviated_well_segment%cased = PETSC_TRUE
                    case('UNCASED','SCREENED','PERFORATED')
                      deviated_well_segment%cased = PETSC_FALSE
                    case default
                      option%io_buffer = "Error in constructing &
                          &WELL_TRAJECTORY: please specify &
                          &whether each segment is CASED or SCREENED/UNCASED."
                      call PrintErrMsg(option)
                  end select
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       dxyz(1))
                  call InputErrorMsg(input,option, &
                                     'SEGMENT_DXYZ x-coordinate', &
                                     error_string)
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       dxyz(2))
                  call InputErrorMsg(input,option, &
                                     'SEGMENT_DXYZ y-coordinate', &
                                     error_string)
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       dxyz(3))
                  if (dabs(deviated_well_segment%dxyz(3)) < epsilon) then
                    ! This is to ensure all procs properly order the
                    ! well segments. Epsilon in m/m.
                    deviated_well_segment%dxyz(3) = -1.d0 * epsilon * &
                        (deviated_well_segment%dxyz(1) + &
                        deviated_well_segment%dxyz(2))
                  endif
                  call InputErrorMsg(input,option, &
                                     'SEGMENT_DXYZ z-coordinate', &
                                     error_string)
                case('SEGMENT_RADIUS_TO_HORIZONTAL_X')
                  if (.not. associated(deviated_well_segment)) then
                    option%io_buffer = "Error in constructing &
                      &WELL_TRAJECTORY: please first specify &
                      &SURFACE_ORIGIN before SEGMENT_RADIUS_TO_HORIZONTAL_X."
                    call PrintErrMsg(option)
                  endif
                  allocate(deviated_well_segment%next)
                  deviated_well_segment => deviated_well_segment%next
                  call WellSegmentInit(deviated_well_segment)
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  select case(word)
                    case('CASED')
                      deviated_well_segment%cased = PETSC_TRUE
                    case('UNCASED','SCREENED','PERFORATED')
                      deviated_well_segment%cased = PETSC_FALSE
                    case default
                      option%io_buffer = "Error in constructing &
                          &WELL_TRAJECTORY: please specify &
                          &whether each segment is CASED or SCREENED/UNCASED."
                      call PrintErrMsg(option)
                  end select
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       radius_to_horizontal_x)
                  call InputErrorMsg(input,option, &
                                     'SEGMENT_RADIUS_TO_HORIZONTAL_X', &
                                     error_string)
                case('SEGMENT_RADIUS_TO_HORIZONTAL_Y')
                  if (.not. associated(deviated_well_segment)) then
                    option%io_buffer = "Error in constructing &
                      &WELL_TRAJECTORY: please first specify &
                      &SURFACE_ORIGIN before SEGMENT_RADIUS_TO_HORIZONTAL_Y."
                    call PrintErrMsg(option)
                  endif
                  allocate(deviated_well_segment%next)
                  deviated_well_segment => deviated_well_segment%next
                  call WellSegmentInit(deviated_well_segment)
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  select case(word)
                    case('CASED')
                      deviated_well_segment%cased = PETSC_TRUE
                    case('UNCASED','SCREENED','PERFORATED')
                      deviated_well_segment%cased = PETSC_FALSE
                    case default
                      option%io_buffer = "Error in constructing &
                          &WELL_TRAJECTORY: please specify &
                          &whether each segment is CASED or SCREENED/UNCASED."
                      call PrintErrMsg(option)
                  end select
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       radius_to_horizontal_y)
                  call InputErrorMsg(input,option, &
                                     'SEGMENT_RADIUS_TO_HORIZONTAL_Y', &
                                     error_string)
                case('SEGMENT_RADIUS_TO_HORIZONTAL_ANGLE')
                  if (.not. associated(deviated_well_segment)) then
                    option%io_buffer = "Error in constructing &
                      &WELL_TRAJECTORY: please first specify &
                      &SURFACE_ORIGIN before &
                      &SEGMENT_RADIUS_TO_HORIZONTAL_ANGLE."
                    call PrintErrMsg(option)
                  endif
                  allocate(deviated_well_segment%next)
                  deviated_well_segment => deviated_well_segment%next
                  call WellSegmentInit(deviated_well_segment)
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  select case(word)
                    case('CASED')
                      deviated_well_segment%cased = PETSC_TRUE
                    case('UNCASED','SCREENED','PERFORATED')
                      deviated_well_segment%cased = PETSC_FALSE
                    case default
                      option%io_buffer = "Error in constructing &
                          &WELL_TRAJECTORY: please specify &
                          &whether each segment is CASED or SCREENED/UNCASED."
                      call PrintErrMsg(option)
                  end select
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       radius_to_horizontal_angle(1))
                  call InputErrorMsg(input,option, &
                                     'SEGMENT_RADIUS_TO_HORIZONTAL_ANGLE &
                                     &radius', error_string)
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       radius_to_horizontal_angle(2))
                  call InputErrorMsg(input,option, &
                                     'SEGMENT_RADIUS_TO_HORIZONTAL_ANGLE &
                                     &angle',error_string)
                case('SEGMENT_RADIUS_TO_VERTICAL')
                  if (.not. associated(deviated_well_segment)) then
                    option%io_buffer = "Error in constructing &
                      &WELL_TRAJECTORY: please first specify &
                      &SURFACE_ORIGIN before SEGMENT_RADIUS_TO_VERTICAL."
                    call PrintErrMsg(option)
                  endif
                  allocate(deviated_well_segment%next)
                  deviated_well_segment => deviated_well_segment%next
                  call WellSegmentInit(deviated_well_segment)
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  select case(word)
                    case('CASED')
                      deviated_well_segment%cased = PETSC_TRUE
                    case('UNCASED','SCREENED','PERFORATED')
                      deviated_well_segment%cased = PETSC_FALSE
                    case default
                      option%io_buffer = "Error in constructing &
                          &WELL_TRAJECTORY: please specify &
                          &whether each segment is CASED or SCREENED/UNCASED."
                      call PrintErrMsg(option)
                  end select
                  call InputReadDouble(input,option,deviated_well_segment% &
                                       radius_to_vertical)
                  call InputErrorMsg(input,option, &
                                     'SEGMENT_RADIUS_TO_VERTICAL', &
                                     error_string)
                case default
                  call InputKeywordUnrecognized(input,word,&
                                             'WELL_TRAJECTORY',option)
              end select
            enddo
            if (associated(deviated_well_segment)) &
                nullify(deviated_well_segment)
            call InputPopBlock(input,option)
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
      if (.not.associated(well_grid%casing) .and. .not. &
          associated(well_grid%deviated_well_segment_list)) then
        option%io_buffer = 'Either keyword CASING or WELL_TRAJECTORY &
                           &must be provided in &
                           &the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
      endif
      if ((Uninitialized(well_grid%tophole(1)) .or. &
          Uninitialized(well_grid%tophole(2)) .or. &
          Uninitialized(well_grid%tophole(3))) .and. .not. &
          associated(well_grid%deviated_well_segment_list)) then
        option%io_buffer = 'ERROR: Either keyword TOP_OF_HOLE &
                           &or WELL_TRAJECTORY &
                           &must be provided in &
                           &the ' // trim(error_string) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if ((Uninitialized(well_grid%bottomhole(1)) .or. &
          Uninitialized(well_grid%bottomhole(2)) .or. &
          Uninitialized(well_grid%bottomhole(3))) .and. .not. &
          associated(well_grid%deviated_well_segment_list)) then
        option%io_buffer = 'ERROR: Either keywordk BOTTOM_OF_HOLE &
                           &or WELL_TRAJECTORY &
                           &must be provided in &
                           &the ' // trim(error_string) // ' block.'
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

  deallocate(temp_well_index)

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
  character(len=MAXWORDLENGTH) :: word, keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscInt :: k
  PetscInt :: num_errors,read_max,num_read
  PetscReal, pointer :: temp_diameter(:)
  PetscReal, pointer :: temp_friction(:)
  PetscReal, pointer :: temp_well_perm(:)
  PetscReal, pointer :: temp_well_phi(:)
  PetscReal, pointer :: temp(:)

  error_string = trim(error_string) // ',WELL'
  found = PETSC_TRUE
  num_errors = 0

  read_max = 9000
  allocate(temp_diameter(read_max))
  allocate(temp_friction(read_max))
  allocate(temp_well_perm(read_max))
  allocate(temp_well_phi(read_max))
  allocate(temp(read_max))
  temp_diameter(:) = UNINITIALIZED_DOUBLE
  temp_friction(:) = UNINITIALIZED_DOUBLE
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
          case('WELL_INDEX_MODEL')
            call InputReadWord(input,option,word,PETSC_TRUE)
            select case(word)
              case('PEACEMAN_ISO')
                pm_well%well%WI_model = PEACEMAN_ISO
              case('PEACEMAN_2D')
                pm_well%well%WI_model = PEACEMAN_2D
              case('PEACEMAN_3D')
                pm_well%well%WI_model = PEACEMAN_3D
              case('SCALE_BY_PERM')
                pm_well%well%WI_model = PEACEMAN_NONE
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
      if (.not.associated(pm_well%well%friction_factor) .and. &
          pm_well%well%WI_model /= PEACEMAN_NONE) then
        option%io_buffer = 'Keyword FRICTION_COEFFICIENT must be provided in &
                           &the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
      endif
      if (.not.associated(pm_well%well%diameter) .and. &
          pm_well%well%WI_model /= PEACEMAN_NONE) then
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

  deallocate(temp_friction)
  deallocate(temp_diameter)
  deallocate(temp_well_phi)
  deallocate(temp_well_perm)
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
  PetscReal :: temp

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
                           &GAS_MASS_RATE',option)
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
          case('TEMPERATURE')
            call InputReadDouble(input,option,temp)
            allocate(pm_well%well%temp(pm_well%well_grid%nsegments))
            pm_well%well%temp = temp
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

          option%io_buffer ='WELL_MODEL_WIPP_DARCY well model needs both &
              &Dirichlet LIQUID_PRESSURE and GAS_SATURATION set in the ' &
              // trim(error_string) // ' block.'
          call PrintErrMsg(option)
        endif
      endif

      if ((Initialized(pm_well%well%th_p).or.Initialized(pm_well%well%th_sg)) &
          .and. .not.(Initialized(pm_well%well%th_p) .and. &
          Initialized(pm_well%well%th_sg))) then
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
    case ('USE_WELL_COUPLER')
      pm_well%use_well_coupler = PETSC_TRUE
    case ('FRACTURE_PRESSURE')
      call InputReadDouble(input,option,pm_well%pressure_threshold_max)
    case ('MINIMUM_PRESSURE')
      call InputReadDouble(input,option,pm_well%pressure_threshold_min)
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

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

  class(pm_well_sequential_type) :: pm_well
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
            call InputReadDouble(input,option,pm_well%flow_soln% &
                                 itol_scaled_res)
            call InputErrorMsg(input,option,'ITOL_SCALED_RESIDUAL', &
                               error_string)
        !-----------------------------
          case('ITOL_ABS_UPDATE_PRESSURE')
            call InputReadDouble(input,option,pm_well%flow_soln% &
                                 itol_abs_update_p)
            call InputErrorMsg(input,option,'ITOL_ABS_UPDATE_PRESSURE', &
                               error_string)
        !-----------------------------
          case('ITOL_ABS_UPDATE_SATURATION')
            call InputReadDouble(input,option,pm_well%flow_soln% &
                                 itol_abs_update_s)
            call InputErrorMsg(input,option,'ITOL_ABS_UPDATE_SATURATION', &
                               error_string)
        !-----------------------------
          case('ITOL_REL_UPDATE_PRESSURE')
            call InputReadDouble(input,option,pm_well%flow_soln% &
                                 itol_rel_update_p)
            call InputErrorMsg(input,option,'ITOL_REL_UPDATE_PRESSURE', &
                               error_string)
        !-----------------------------
          case('ITOL_REL_UPDATE_SATURATION')
            call InputReadDouble(input,option,pm_well%flow_soln% &
                                 itol_rel_update_s)
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

  class(pm_well_sequential_type) :: pm_well
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
            option%io_buffer = 'WELL_TRANSPORT_SOLVER, TIMESTEP_CUT_FACTOR &
              &has been depreciated for the wellbore model. All time step &
              &cuts are controlled by the transport mode for the reservoir.'
            call PrintErrMsg(option)
        !-----------------------------
          case('TIMESTEP_RAMP_FACTOR')
            option%io_buffer = 'WELL_TRANSPORT_SOLVER, TIMESTEP_RAMP_FACTOR &
              &has been depreciated for the wellbore model. All time step &
              &cuts are controlled by the transport mode for the reservoir.'
            call PrintErrMsg(option)
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

  class(pm_well_sequential_type) :: pm_well
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

  input%ierr = INPUT_ERROR_NONE
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
      case('WELL_GRID')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadCard(input,option,keyword)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(keyword)
          select case('keyword')
            case('WELL_TRAJECTORY')
              call InputPopBlock(input,option)
          end select
        enddo
        call InputSkipToEND(input,option,card)
      case('WELL','WELL_MODEL_TYPE','WELL_CONSTRAINT_TYPE', &
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

subroutine PMWellSetRealization(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  use Realization_Subsurface_class

  implicit none

  class(pm_well_type) :: this

  this%realization => RealizationCast(this%realization_base)

end subroutine PMWellSetRealization

! ************************************************************************** !

recursive subroutine PMWellInitializeRunBase(this)
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

  type(option_type), pointer :: option

  option => this%option

  option%io_buffer = "PMWellInitializeRunBase must be extended."
  call PrintErrMsg(option)


end subroutine PMWellInitializeRunBase

! ************************************************************************** !

recursive subroutine PMWellInitializeRunSequential(this)
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

  class(pm_well_sequential_type) :: this

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
  this%srcsink_water = 0.d0
  this%srcsink_gas = 0.d0

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
    allocate(this%well_pert(k)% &
             friction_factor(size(this%well%friction_factor)))
    allocate(this%well_pert(k)%r0(size(this%well%r0)))
    allocate(this%well_pert(k)%skin(size(this%well%skin)))
    this%well_pert(k)%diameter = this%well%diameter
    this%well_pert(k)%friction_factor = this%well%friction_factor
    this%well_pert(k)%well_model_type = this%well%well_model_type
    call PMWellInitWellVars(this%well_pert(k),this%well_grid, &
                          this%transport,nsegments,this%nspecies)
    call PMWellInitFluidVars(this%well_pert(k),nsegments, &
                           this%option%flow%reference_density,option)
    call PMWellInitRes(this%well_pert(k)%reservoir,nsegments,this%nspecies, &
                       option)
    call PMWellInitRes(this%well_pert(k)%reservoir_save,nsegments, &
                       this%nspecies,option)
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
  call PMWellSetPlotVariablesSequential(output_var_list,this)
  if (.not.associated( &
                this%realization%output_option%output_snap_variable_list, &
                this%realization%output_option%output_obs_variable_list)) then
    output_var_list => this%realization%output_option%output_obs_variable_list
    call PMWellSetPlotVariablesSequential(output_var_list,this)
  endif

  if (this%print_well) then
    call PMWellOutputHeaderSequential(this)
  endif

end subroutine PMWellInitializeRunSequential

! ************************************************************************** !

recursive subroutine PMWellInitializeRunHydrostatic(this)
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

  class(pm_well_hydrostatic_type) :: this

  type(output_variable_list_type), pointer :: output_var_list
  type(option_type), pointer :: option
  PetscInt :: nsegments, k

  option => this%option

  nsegments = this%well_grid%nsegments

  ! srcsink_water/gas is indexed (0,:) unperturbed value
  !                              (1,:) perturbed wrt gas pressure
  !                              (2,:) perturbed wrt gas saturation
  allocate(this%srcsink_water(0:2,nsegments))
  allocate(this%srcsink_gas(0:2,nsegments))
  this%srcsink_water = 0.d0
  this%srcsink_gas = 0.d0

  call PMWellInitWellVars(this%well,this%well_grid,PETSC_FALSE, &
                          nsegments,ZERO_INTEGER)
  call PMWellInitFluidVars(this%well,nsegments, &
                           this%option%flow%reference_density,option)

  do k = 1,option%nflowdof
    allocate(this%well_pert(k)%diameter(nsegments))
    allocate(this%well_pert(k)% &
             friction_factor(size(this%well%friction_factor)))
    allocate(this%well_pert(k)%r0(size(this%well%r0)))
    allocate(this%well_pert(k)%skin(size(this%well%skin)))
    this%well_pert(k)%diameter = this%well%diameter
    this%well_pert(k)%friction_factor = this%well%friction_factor
    this%well_pert(k)%well_model_type = this%well%well_model_type
    call PMWellInitWellVars(this%well_pert(k),this%well_grid, &
                          PETSC_FALSE,nsegments,ZERO_INTEGER)
    call PMWellInitFluidVars(this%well_pert(k),nsegments, &
                           this%option%flow%reference_density,option)
    call PMWellInitRes(this%well_pert(k)%reservoir,nsegments,ZERO_INTEGER, &
                       option)
    call PMWellInitRes(this%well_pert(k)%reservoir_save,nsegments, &
                       ZERO_INTEGER,option)
  enddo

  call PMWellInitRes(this%well%reservoir,nsegments,ZERO_INTEGER,option)
  call PMWellInitRes(this%well%reservoir_save,nsegments,ZERO_INTEGER,option)



  ! Initialize perturbations
  allocate(this%pert(nsegments,option%nflowdof))
  this%pert = 0.d0

  ! Setup the output variables for snapshot and observation files
  output_var_list => this%realization%output_option%output_snap_variable_list
  call PMWellSetPlotVariablesHydrostatic(output_var_list,this)
  if (.not.associated( &
                this%realization%output_option%output_snap_variable_list, &
                this%realization%output_option%output_obs_variable_list)) then
    output_var_list => this%realization%output_option%output_obs_variable_list
    call PMWellSetPlotVariablesHydrostatic(output_var_list,this)
  endif

  if (this%print_well) then
    call PMWellOutputHeaderHydrostatic(this)
  endif

end subroutine PMWellInitializeRunHydrostatic

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
  if (.not. associated(well%temp)) then
    allocate(well%temp(nsegments))
    well%temp = UNINITIALIZED_DOUBLE
  endif
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

  well%ql = 0.d0
  well%qg = 0.d0
  well%ql_bc = 0.d0
  well%qg_bc = 0.d0
  well%ql_kmol = 0.d0
  well%qg_kmol = 0.d0
  well%ql_kmol_bc = 0.d0
  well%qg_kmol_bc = 0.d0
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
    allocate(well%aqueous_mass_q_cumulative(nspecies,nsegments))
    allocate(well%aqueous_mass_Qcumulative(nspecies,nsegments))
    well%aqueous_mass_q_cumulative(:,:) = 0.d0
    well%aqueous_mass_Qcumulative(:,:) = 0.d0
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
  well%liq%visc = 0.d0
  allocate(well%liq%H(nsegments))
  well%liq%H = 0.d0
  allocate(well%liq%Q(nsegments))
  well%liq%Q = 0.d0
  allocate(well%liq%Q_cumulative(nsegments))
  well%liq%Q_cumulative = 0.d0
  allocate(well%liq%kr(nsegments))
  well%liq%kr = 0.d0
  allocate(well%liq%xmass(nsegments,option%nflowspec))
  well%liq%xmass = 0.d0

  allocate(well%gas%s(nsegments))
  well%gas%den0 = reference_density(2)
  allocate(well%gas%den(nsegments))
  well%gas%den(:) = well%gas%den0
  allocate(well%gas%visc(nsegments))
  well%gas%visc = 0.d0
  allocate(well%gas%H(nsegments))
  well%gas%H = 0.d0
  allocate(well%gas%Q(nsegments))
  well%gas%Q = 0.d0
  allocate(well%gas%kr(nsegments))
  well%gas%kr = 0.d0
  allocate(well%gas%xmass(nsegments,option%nflowspec))
  well%gas%xmass = 0.d0

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
  allocate(reservoir%H_l(nsegments))
  allocate(reservoir%H_g(nsegments))
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
  allocate(reservoir%tmp_flow(nsegments,20))

  allocate(reservoir%aqueous_conc(idof,nsegments))
  allocate(reservoir%aqueous_mass(idof,nsegments))
  allocate(reservoir%tmp_tran(idof,nsegments,2))

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
  allocate(tran_soln%prev_soln%resr_aqueous_conc(nspecies,nsegments))

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

subroutine PMWellInitializeTimestepBase(this)
  !
  ! Initializes and takes the time step for the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this

end subroutine PMWellInitializeTimestepBase

! ************************************************************************** !

subroutine PMWellInitializeTimestepSeq(this)
  !
  ! Initializes and takes the time step for the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_sequential_type) :: this

  this%option%io_buffer = "PMWellInitializeTimestepSeq must be extended."
  call PrintErrMsg(this%option)

end subroutine PMWellInitializeTimestepSeq

! ************************************************************************** !

subroutine PMWellInitializeWellFlow(pm_well)
  !
  ! Initializes the well for the first time step for flow.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  use SCO2_Aux_module
  use Hydrate_Aux_module

  implicit none

  class(pm_well_type) :: pm_well

  type(well_reservoir_type), pointer :: reservoir
  type(strata_type), pointer :: strata
  type(sco2_auxvar_type), pointer :: sco2_auxvar
  type(hydrate_auxvar_type), pointer :: hyd_auxvar
  type(option_type), pointer :: option
  PetscInt :: k, ghosted_id

  option => pm_well%option
  reservoir => pm_well%well%reservoir

  ! set initial flow parameters to the reservoir flow parameters
  pm_well%well%pl = reservoir%p_l
  pm_well%well%pg = reservoir%p_g
  if (Uninitialized(pm_well%well%temp(1))) pm_well%well%temp = reservoir%temp
  pm_well%well%liq%s = reservoir%s_l
  pm_well%well%gas%s = reservoir%s_g
  pm_well%well%liq%den = reservoir%den_l
  pm_well%well%gas%den = reservoir%den_g
  pm_well%well%liq%visc = reservoir%visc_l
  pm_well%well%gas%visc = reservoir%visc_g
  pm_well%well%liq%xmass = reservoir%xmass_liq
  pm_well%well%gas%xmass = reservoir%xmass_gas
  pm_well%well%liq%H = reservoir%H_l
  pm_well%well%gas%H = reservoir%H_g
  ! BHP can be a flow primary variable if fully coupled to flow.
  select case (option%iflowmode)
    case (SCO2_MODE)
      ghosted_id = pm_well%well_grid%h_ghosted_id(1)
      sco2_auxvar => &
        pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)
        pm_well%well%bh_p = sco2_auxvar%pres(option%gas_phase) - &
                            pm_well%well%liq%den(1) * &
                            pm_well%option%gravity(Z_DIRECTION) * &
                            pm_well%well_grid%dh(1)/2.d0
    case (H_MODE)
      ghosted_id = pm_well%well_grid%h_ghosted_id(1)
      hyd_auxvar => &
        pm_well%realization%patch%aux%hydrate%auxvars(ZERO_INTEGER,ghosted_id)
        pm_well%well%bh_p = hyd_auxvar%pres(option%gas_phase) - &
                            pm_well%well%liq%den(1) * &
                            pm_well%option%gravity(Z_DIRECTION) * &
                            pm_well%well_grid%dh(1)/2.d0
  end select
  ! update the Darcy fluxes within the well
  do k = 1,option%nflowdof

    pm_well%well_pert(k)%pl = reservoir%p_l
    pm_well%well_pert(k)%pg = reservoir%p_g
    if (Uninitialized(pm_well%well_pert(k)%temp(1))) then
      pm_well%well_pert(k)%temp = reservoir%temp
    endif
    pm_well%well_pert(k)%liq%s = reservoir%s_l
    pm_well%well_pert(k)%gas%s = reservoir%s_g
    pm_well%well_pert(k)%liq%den = reservoir%den_l
    pm_well%well_pert(k)%gas%den = reservoir%den_g
    pm_well%well_pert(k)%liq%visc = reservoir%visc_l
    pm_well%well_pert(k)%gas%visc = reservoir%visc_g
    pm_well%well_pert(k)%liq%xmass = reservoir%xmass_liq
    pm_well%well_pert(k)%gas%xmass = reservoir%xmass_gas
    pm_well%well_pert(k)%liq%H = reservoir%H_l
    pm_well%well_pert(k)%gas%H = reservoir%H_g
    pm_well%well_pert(k)%bh_p = pm_well%well%bh_p

  enddo

  ! Link well material properties
  if (pm_well%well_comm%comm /= MPI_COMM_NULL) then
    do k = 1,pm_well%well_grid%nsegments
      strata => pm_well%strata_list%first
      do
        if (.not.associated(strata)) exit
        if (strata%id == pm_well%well_grid%strata_id(k)) then
          pm_well%well%ccid(k) = strata%material_property% &
                                 saturation_function_id
          pm_well%well%permeability(k) = strata%material_property% &
                                         permeability(3,3)
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

  if (option%iflowmode /= SCO2_MODE .and. option%iflowmode /= H_MODE) &
     initialize_well_flow = PETSC_FALSE

end subroutine PMWellInitializeWellFlow

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
    case (H_MODE)
      call PMWellHydratePerturb(this)
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
  PetscInt :: idof, k, i
  PetscInt :: ghosted_id
  PetscInt :: lid
  PetscReal :: pres_bump
  PetscErrorCode :: ierr

  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-15
  ! PetscReal :: dpl, dpg

  option => pm_well%option
  res_grid => pm_well%realization%patch%grid

  lid = option%liquid_phase

  ! Make sure well is actually flowing if injection/production rates
  ! are specified.
  pres_bump = 0.d0
  i = 0
  well => pm_well%well
  do
    call pm_well%SolveFlow(ZERO_INTEGER,ierr)
    if (any(dabs(well%gas%Q) > 0.d0) .or. &
          any(dabs(well%liq%Q) > 0.d0)) exit
    if (well%th_ql > 0.d0 .or. &
        well%th_qg > 0.d0) then
      ! Injection well
      pres_bump = 1.25d0 * (pres_bump + 1.d0)
    elseif (well%total_rate < 0.d0) then
      ! Extraction well
      pres_bump = 1.25d0 * (pres_bump - 1.d0)
    else
      ! No need for this
      exit
    endif
    well%bh_p = well%bh_p + pres_bump
    i = i + 1
    if (i > 100) then
      option%io_buffer = "Exceeded maximum number of iterations (100) to &
                          &recover vanishing well Jacobian."
      call PrintErrMsg(option)
    endif
  enddo

  ! Go up the well: for each well segment that is on-process, update
  ! perturbed values for reservoir and well variables.

  ! Perturbed fluxes wrt reservoir variables
  do k = 1,pm_well%well_grid%nsegments
    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    do idof = 1,option%nflowdof
      well => pm_well%well_pert(idof)
      well%bh_p = pm_well%well%bh_p
      if (idof == SCO2_WELL_DOF) then
        sco2_auxvar => &
          pm_well%realization%patch%aux%sco2%auxvars(ZERO_INTEGER,ghosted_id)
        well%bh_p = well%bh_p + pm_well%realization%patch%aux%sco2% &
                    auxvars(SCO2_WELL_DOF,ghosted_id)%pert
        if (.not. Initialized(sco2_auxvar%well%pressure_bump)) then
          if (k == 1) then
            sco2_auxvar%well%pressure_bump = (well%bh_p - &
                      sco2_auxvar%pres(lid)) + pres_bump
          else
            sco2_auxvar%well%pressure_bump = 0.d0
          endif
        else
          sco2_auxvar%well%pressure_bump = pres_bump
        endif
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
      reservoir%H_l(k) = sco2_auxvar%H(option%liquid_phase)
      reservoir%H_g(k) = sco2_auxvar%H(option%gas_phase)
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

  ! Now update fluxes associated with perturbed values
  do idof = 1,option%nflowdof
    call pm_well%SolveFlow(idof,ierr)
  enddo

end subroutine PMWellSCO2Perturb

! ************************************************************************** !

subroutine PMWellHydratePerturb(pm_well)
  !
  ! Perturb well variables when using Hydrate flow mode.
  !
  ! Author: Michael Nole
  ! Date: 06/06/2024
  !
  use Option_module
  use Hydrate_Aux_module
  use Material_Aux_module
  use Grid_module
  implicit none
  class(pm_well_type) :: pm_well
  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: reservoir
  type(hydrate_auxvar_type), pointer :: hyd_auxvar
  type(material_auxvar_type), pointer :: material_auxvar
  type(grid_type), pointer :: res_grid
  type(option_type), pointer :: option
  PetscInt :: idof, k, i
  PetscInt :: ghosted_id
  PetscInt :: lid
  PetscReal :: pres_bump
  PetscErrorCode :: ierr
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-15
  ! PetscReal :: dpl, dpg
  option => pm_well%option
  res_grid => pm_well%realization%patch%grid
  lid = option%liquid_phase
  ! Make sure well is actually flowing if injection/production rates
  ! are specified.
  pres_bump = 0.d0
  i = 0
  well => pm_well%well
  do
    call pm_well%SolveFlow(ZERO_INTEGER,ierr)
    if (any(dabs(well%gas%Q) > 0.d0) .or. &
          any(dabs(well%liq%Q) > 0.d0)) exit
    if (well%th_ql > 0.d0 .or. &
        well%th_qg > 0.d0) then
      ! Injection well
      pres_bump = 1.25d0 * (pres_bump + 1.d0)
    elseif (well%th_ql < 0.d0 .or. &
            well%th_qg < 0.d0) then
      ! Extraction well
      pres_bump = 1.25d0 * (pres_bump - 1.d0)
    else
      ! No need for this
      exit
    endif
    well%bh_p = well%bh_p + pres_bump
    i = i + 1
    if (i > 100) then
      option%io_buffer = "Exceeded maximum number of iterations (100) to &
                          &recover vanishing well Jacobian."
      call PrintErrMsg(option)
    endif
  enddo
  ! Go up the well: for each well segment that is on-process, update
  ! perturbed values for reservoir and well variables.

  ! Perturbed fluxes wrt reservoir variables
  do k = 1,pm_well%well_grid%nsegments
    ghosted_id = pm_well%well_grid%h_ghosted_id(k)
    do idof = 1,option%nflowdof
      well => pm_well%well_pert(idof)
      well%bh_p = pm_well%well%bh_p
      if (idof == HYDRATE_WELL_DOF) then
        hyd_auxvar => &
          pm_well%realization%patch%aux%hydrate%auxvars(ZERO_INTEGER,ghosted_id)
        well%bh_p = well%bh_p + pm_well%realization%patch%aux%hydrate% &
                    auxvars(HYDRATE_WELL_DOF,ghosted_id)%pert
        if (.not. Initialized(hyd_auxvar%well%pressure_bump)) then
          if (k == 1) then
            hyd_auxvar%well%pressure_bump = (well%bh_p - &
                      hyd_auxvar%pres(lid)) + pres_bump
          else
            hyd_auxvar%well%pressure_bump = 0.d0
          endif
        else
          hyd_auxvar%well%pressure_bump = pres_bump
        endif
      else
        hyd_auxvar => &
          pm_well%realization%patch%aux%hydrate%auxvars(idof,ghosted_id)
        well%bh_p = pm_well%well%bh_p
      endif
      reservoir => well%reservoir
      material_auxvar => &
        pm_well%realization%patch%aux%material%auxvars(ghosted_id)
      reservoir%p_l(k) = hyd_auxvar%pres(option%liquid_phase)
      reservoir%p_g(k) = hyd_auxvar%pres(option%gas_phase)
      reservoir%s_l(k) = hyd_auxvar%sat(option%liquid_phase)
      reservoir%s_g(k) = hyd_auxvar%sat(option%gas_phase)
      reservoir%temp(k) = hyd_auxvar%temp
      reservoir%mobility_l(k) = &
        hyd_auxvar%mobility(option%liquid_phase)
      reservoir%mobility_g(k) = hyd_auxvar%mobility(option%gas_phase)
      reservoir%kr_l(k) = hyd_auxvar%kr(option%liquid_phase)
      reservoir%kr_g(k) = hyd_auxvar%kr(option%gas_phase)
      reservoir%den_l(k) = hyd_auxvar%den_kg(option%liquid_phase)
      reservoir%den_g(k) = hyd_auxvar%den_kg(option%gas_phase)
      reservoir%visc_l(k) = hyd_auxvar%visc(option%liquid_phase)
      reservoir%visc_g(k) = hyd_auxvar%visc(option%gas_phase)
      reservoir%H_l(k) = hyd_auxvar%H(option%liquid_phase)
      reservoir%H_g(k) = hyd_auxvar%H(option%gas_phase)
      reservoir%xmass_liq(k,:) = hyd_auxvar%xmass(:, &
                                                  option%liquid_phase)
      reservoir%xmass_gas(k,:) = hyd_auxvar%xmass(:, &
                                                  option%gas_phase)
      reservoir%e_por(k) = hyd_auxvar%effective_porosity
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
  ! Now update fluxes associated with perturbed values
  do idof = 1,option%nflowdof
    call pm_well%SolveFlow(idof,ierr)
  enddo
end subroutine PMWellHydratePerturb

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

  if (pm_well%well%well_model_type == WELL_MODEL_HYDROSTATIC .or. &
      pm_well%well%well_model_type == WELL_MODEL_COAXIAL .or. &
      pm_well%well%well_model_type == WELL_MODEL_U_SHAPE) return

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
      any(pm_well%well_comm%petsc_rank_list == pm_well%option%myrank) .and. &
      pm_well%option%iflowmode == WF_MODE) then
    pm_well%option%io_buffer =  'At least one WELLBORE_MODEL grid segment has &
        &not been assigned with a REGION and MATERIAL_PROPERTY with the use &
        &of the STRATA block.'
    call PrintErrMsg(pm_well%option)
  endif

end subroutine PMWellUpdateStrata

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
  PetscInt :: k, indx
  PetscInt :: ghosted_id
  PetscInt :: tag
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

    if (k == 1 .and. indx == 0 .and. Initialized(sco2_auxvar%well%bh_p)) &
      pm_well%well%bh_p = sco2_auxvar%well%bh_p

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
    reservoir%H_l(k) = sco2_auxvar%H(option%liquid_phase)
    reservoir%H_g(k) = sco2_auxvar%H(option%gas_phase)
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
    tag = 0
    do k = 2,pm_well%well_grid%nsegments
      if (option%myrank == pm_well%well_grid%h_rank_id(1)) then
        if (option%myrank == pm_well%well_grid%h_rank_id(k)) cycle
        call MPI_Send(pm_well%well%bh_p,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                      pm_well%well_grid%h_rank_id(k),tag,option%mycomm,ierr); &
                      CHKERRQ(ierr)
      elseif (option%myrank == pm_well%well_grid%h_rank_id(k)) then
        call MPI_Recv(pm_well%well%bh_p,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                      pm_well%well_grid%h_rank_id(1),tag,option%mycomm, &
                      MPI_STATUS_IGNORE,ierr); CHKERRQ(ierr)
      endif
    enddo
  endif

end subroutine PMWellUpdateReservoirSCO2

! ************************************************************************** !
subroutine PMWellUpdateReservoirHydrate(pm_well,update_index,segment_index)
  !
  ! Updates the Hydrate mode reservoir properties for the well process model.
  !
  ! Author: Michael Nole
  ! Date: 03/06/2024

  use Hydrate_Aux_module
  use Material_Aux_module
  use Grid_module

  implicit none

  class(pm_well_type) :: pm_well
  PetscInt :: update_index
  PetscInt :: segment_index

  type(well_reservoir_type), pointer :: reservoir
  type(hydrate_auxvar_type), pointer :: hyd_auxvar
  type(material_auxvar_type), pointer :: material_auxvar
  type(grid_type), pointer :: res_grid
  type(option_type), pointer :: option
  type(well_comm_type), pointer :: well_comm
  PetscInt :: k, indx
  PetscInt :: ghosted_id
  PetscInt :: tag
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

    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    if (Initialized(segment_index)) then
      ! Use reservoir perturbation: only perturb one reservoir variable
      ! at a time.
      if (k == segment_index) then
        hyd_auxvar => &
          pm_well%realization%patch%aux%hydrate%auxvars(indx,ghosted_id)
      else
        hyd_auxvar => &
          pm_well%realization%patch%aux%hydrate%auxvars(ZERO_INTEGER,ghosted_id)
      endif
    else
      hyd_auxvar => &
        pm_well%realization%patch%aux%hydrate%auxvars(ZERO_INTEGER,ghosted_id)
    endif

    if (k == 1 .and. indx == 0 .and. Initialized(hyd_auxvar%well%bh_p)) &
        pm_well%well%bh_p = hyd_auxvar%well%bh_p

    material_auxvar => &
      pm_well%realization%patch%aux%material%auxvars(ghosted_id)

    reservoir%p_l(k) = hyd_auxvar%pres(option%liquid_phase)
    reservoir%p_g(k) = hyd_auxvar%pres(option%gas_phase)
    reservoir%s_l(k) = hyd_auxvar%sat(option%liquid_phase)
    reservoir%s_g(k) = hyd_auxvar%sat(option%gas_phase)
    reservoir%temp(k) = hyd_auxvar%temp
    reservoir%mobility_l(k) = &
      hyd_auxvar%mobility(option%liquid_phase)
    reservoir%mobility_g(k) = hyd_auxvar%mobility(option%gas_phase)
    reservoir%kr_l(k) = hyd_auxvar%kr(option%liquid_phase)
    reservoir%kr_g(k) = hyd_auxvar%kr(option%gas_phase)
    reservoir%den_l(k) = hyd_auxvar%den_kg(option%liquid_phase)
    reservoir%den_g(k) = hyd_auxvar%den_kg(option%gas_phase)
    reservoir%visc_l(k) = hyd_auxvar%visc(option%liquid_phase)
    reservoir%visc_g(k) = hyd_auxvar%visc(option%gas_phase)
    reservoir%H_l(k) = hyd_auxvar%H(option%liquid_phase)
    reservoir%H_g(k) = hyd_auxvar%H(option%gas_phase)
    reservoir%xmass_liq(k,:) = hyd_auxvar%xmass(:, &
                                                    option%liquid_phase)
    reservoir%xmass_gas(k,:) = hyd_auxvar%xmass(:, &
                                                    option%gas_phase)
    reservoir%e_por(k) = hyd_auxvar%effective_porosity

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
    tag = 0
    do k = 2,pm_well%well_grid%nsegments
      if (option%myrank == pm_well%well_grid%h_rank_id(1)) then
        if (option%myrank == pm_well%well_grid%h_rank_id(k)) cycle
        call MPI_Send(pm_well%well%bh_p,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                      pm_well%well_grid%h_rank_id(k),TAG,option%mycomm,ierr); &
                      CHKERRQ(ierr)
      elseif (option%myrank == pm_well%well_grid%h_rank_id(k)) then
        call MPI_Recv(pm_well%well%bh_p,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                      pm_well%well_grid%h_rank_id(1),TAG,option%mycomm, &
                      MPI_STATUS_IGNORE,ierr); CHKERRQ(ierr)
      endif
    enddo
  endif

end subroutine PMWellUpdateReservoirHydrate
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

subroutine PMWellFinalizeTimestepBase(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_type) :: this

  this%option%io_buffer = "PMWellFinalizeTimestepBase must be extended."
  call PrintErrMsg(this%option)

end subroutine PMWellFinalizeTimestepBase

! ************************************************************************** !

subroutine PMWellFinalizeTimestepSequential(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_sequential_type) :: this

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
    call PMWellOutputSequential(this)
  endif

end subroutine PMWellFinalizeTimestepSequential

! ************************************************************************** !

subroutine PMWellFinalizeTimestepHydrostatic(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  implicit none

  class(pm_well_hydrostatic_type) :: this

  call PMWellUpdateReservoirSrcSinkFlow(this)

  if (this%print_well) then
    call PMWellOutputHydrostatic(this)
  endif

end subroutine PMWellFinalizeTimestepHydrostatic

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
  use String_module
  use Hydrate_Aux_module, only: hydrate_well_coupling, &
                             HYDRATE_FULLY_IMPLICIT_WELL

  implicit none

  class(pm_well_type) :: pm_well

  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: srcsink_name
  type(coupler_type), pointer :: source_sink
  PetscInt :: k, ghosted_id
  PetscReal :: well_delta_liq, well_delta_gas
  PetscReal :: density_avg

  option => pm_well%option

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  do k = 1,pm_well%well_grid%nsegments
    if (pm_well%well_grid%h_rank_id(k) /= option%myrank) cycle
    srcsink_name = trim(pm_well%name) // '_well_segment_' // StringWrite(k)

    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    ! [kg-liq/m3]
    density_avg = 0.5d0 * (pm_well%well%liq%den(k) + pm_well%well%reservoir% &
                  den_l(k))

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
              ! access to Q is needed in NWT Mode, so load Q into 3 & 4:
              source_sink%flow_condition%general%rate%dataset%rarray(3) = &
                -1.d0 * pm_well%well%liq%Q(k) ! [kmol/s]
              source_sink%flow_condition%general%rate%dataset%rarray(4) = &
                -1.d0 * pm_well%well%gas%Q(k) ! [kmol/s]
            else
              source_sink%flow_condition%general%rate%dataset%rarray(1) = &
                -1.d0 * pm_well%well%liq%Q(k) ! [kmol/s]
              source_sink%flow_condition%general%rate%dataset%rarray(2) = &
                -1.d0 * pm_well%well%gas%Q(k) ! [kmol/s]
              ! access to Q is needed in NWT Mode, so load Q into 3 & 4:
              source_sink%flow_condition%general%rate%dataset%rarray(3) = &
                -1.d0 * pm_well%well%liq%Q(k) ! [kmol/s]
              source_sink%flow_condition%general%rate%dataset%rarray(4) = &
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
            source_sink%flow_condition%sco2%rate%dataset%rarray(1) = &
              -1.d0 * pm_well%well%liq%Q(k) ! [kg/s]
            source_sink%flow_condition%sco2%rate%dataset%rarray(2) = &
              -1.d0 * pm_well%well%gas%Q(k) ! [kg/s]

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

            pm_well%realization%patch%aux%sco2% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%pl = pm_well%well%pl(k)
            pm_well%realization%patch%aux%sco2% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%pg = pm_well%well%pg(k)
            pm_well%realization%patch%aux%sco2% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%sl = pm_well%well%liq%s(k)
            pm_well%realization%patch%aux%sco2% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%sg = pm_well%well%gas%s(k)
            pm_well%realization%patch%aux%sco2% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%dpl = well_delta_liq
            pm_well%realization%patch%aux%sco2% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%dpg = well_delta_gas
            pm_well%realization%patch%aux%sco2% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%Ql = pm_well%well%liq%Q(k)
            pm_well%realization%patch%aux%sco2% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%Qg = pm_well%well%gas%Q(k)
            pm_well%realization%patch%aux%sco2% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%bh_p = pm_well%well%bh_p

          case(H_MODE)
            if (hydrate_well_coupling == HYDRATE_FULLY_IMPLICIT_WELL) then
              source_sink%flow_condition%hydrate%rate%dataset%rarray(1) = &
                                      -1.d0 * pm_well%well%liq%Q(k) ! [kg/s]
              source_sink%flow_condition%hydrate%rate%dataset%rarray(2) = &
                                      -1.d0 * pm_well%well%gas%Q(k) ! [kg/s]
            else
              source_sink%flow_condition%hydrate%rate%dataset%rarray(1) = &
                -1.d0 * pm_well%well%liq%Q(k) ! [kmol/s]
              source_sink%flow_condition%hydrate%rate%dataset%rarray(2) = &
                -1.d0 * pm_well%well%gas%Q(k) ! [kmol/s]
            endif
            source_sink%flow_condition%hydrate%liquid_pressure%aux_real(1) = &
                                                            pm_well%well%pl(k)
            source_sink%flow_condition%hydrate%gas_pressure%aux_real(1) = &
                                                            pm_well%well%pg(k)
            well_delta_liq = pm_well%well%pl(k) - pm_well%well%reservoir%p_l(k)
            well_delta_gas = pm_well%well%pg(k) - pm_well%well%reservoir%p_g(k)
            source_sink%flow_condition%hydrate%liquid_pressure%aux_real(2) = &
                                                               well_delta_liq
            source_sink%flow_condition%hydrate%gas_pressure%aux_real(2) = &
                                                               well_delta_gas
            source_sink%flow_condition%well%aux_real(1) = density_avg ! kg/m3
            pm_well%realization%patch%aux%hydrate% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%pl = pm_well%well%pl(k)
            pm_well%realization%patch%aux%hydrate% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%pg = pm_well%well%pg(k)
            pm_well%realization%patch%aux%hydrate% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%sl = pm_well%well%liq%s(k)
            pm_well%realization%patch%aux%hydrate% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%sg = pm_well%well%gas%s(k)
            pm_well%realization%patch%aux%hydrate% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%dpl = well_delta_liq
            pm_well%realization%patch%aux%hydrate% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%dpg = well_delta_gas
            pm_well%realization%patch%aux%hydrate% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%Ql = pm_well%well%liq%Q(k)
            pm_well%realization%patch%aux%hydrate% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%Qg = pm_well%well%gas%Q(k)
            pm_well%realization%patch%aux%hydrate% &
              auxvars(ZERO_INTEGER,ghosted_id)%well%bh_p = pm_well%well%bh_p
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
  use String_module

  implicit none

  class(pm_well_sequential_type) :: pm_well

  type(option_type), pointer :: option
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

    srcsink_name = trim(pm_well%name) // '_well_segment_' // StringWrite(k)

    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    ! [kg-liq/m3]
    density_avg = 0.5d0 * (pm_well%well%liq%den(k) + &
                  pm_well%well%reservoir%den_l(k))

    source_sink => pm_well%realization%patch%source_sink_list%first
    do
      if (.not.associated(source_sink)) exit

      if (trim(srcsink_name) == trim(source_sink%name)) then
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

subroutine PMWellUpdateReservoirConcTran(pm_well)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 01/16/2025
  !
  use NW_Transport_Aux_module
  use Material_Aux_module

  implicit none

  class(pm_well_type) :: pm_well

  type(nw_transport_auxvar_type), pointer :: nwt_auxvar
  PetscInt :: ghosted_id, k

  do k = 1,pm_well%well_grid%nsegments
    if (pm_well%well_grid%h_rank_id(k) /= pm_well%option%myrank) cycle
    ghosted_id = pm_well%well_grid%h_ghosted_id(k)
    nwt_auxvar => pm_well%realization%patch%aux%nwt%auxvars(ghosted_id)

    ! aqueous concentration [mol-species/m3-liq]
    pm_well%well%reservoir%aqueous_conc(:,k) = nwt_auxvar%aqueous_eq_conc(:)
    ! aqueous_mass = aq_conc * e_por * volume * s_l
    pm_well%well%reservoir%aqueous_mass(:,k) = nwt_auxvar%aqueous_eq_conc(:) * &
      pm_well%well%reservoir%volume(k) * pm_well%well%reservoir%e_por(k) * &
      pm_well%well%reservoir%s_l(k)
  enddo

  end subroutine PMWellUpdateReservoirConcTran

! ************************************************************************** !

subroutine PMWellUpdateFlowPropertiesBase(this,pert,index)

  implicit none

  class(pm_well_type) :: this
  PetscBool :: pert
  PetscInt :: index

  this%option%io_buffer = "PMWellUpdateFlowPropertiesBase must be extended."
  call PrintErrMsg(this%option)

end subroutine PMWellUpdateFlowPropertiesBase

! ************************************************************************** !

subroutine PMWellUpdateFlowRatesBase(this,well_pert,res_pert, &
                                     segment_index,ierr)

  implicit none

  class(pm_well_type) :: this
  PetscErrorCode :: ierr
  PetscInt :: well_pert
  PetscInt :: res_pert
  PetscInt :: segment_index

  this%option%io_buffer = "PMWellUpdateFlowRatesBase must be extended."
  call PrintErrMsg(this%option)

end subroutine PMWellUpdateFlowRatesBase

! ************************************************************************** !

subroutine PMWellUpdateFlowRatesHydrostatic(this,well_pert,res_pert, &
                                            segment_index,ierr)
  !
  ! This subroutine performs the well rate computation when called from the
  ! fully- or quasi-coupled source/sink update.
  !
  ! Author: Michael Nole
  ! Date: 01/16/2023
  !

  Use Option_module

  implicit none

  class(pm_well_hydrostatic_type) :: this
  PetscErrorCode :: ierr
  PetscInt :: well_pert
  PetscInt :: res_pert
  PetscInt :: segment_index

  type(option_type), pointer :: option
  PetscInt :: k

  option => this%option

  this%print_output = PETSC_FALSE
  select case (option%iflowmode)
    case(SCO2_MODE)
      call PMWellUpdateReservoirSCO2(this,res_pert,segment_index)
    case(H_MODE)
      call PMWellUpdateReservoirHydrate(this,res_pert,segment_index)
  end select

  if (initialize_well_flow) then
    call PMWellInitializeWellFlow(this)
    do k = 1,option%nflowdof
      call PMWellCopyWell(this%well,this%well_pert(k),PETSC_FALSE)
    enddo
  endif
  call this%SolveFlow(well_pert,ierr)
  this%print_output = PETSC_TRUE

end subroutine PMWellUpdateFlowRatesHydrostatic

! ************************************************************************** !

subroutine PMWellModifyFlowResBase(this,residual)
  !
  ! Author: Michael Nole
  ! Date: 02/04/2025
  !

  implicit none

  class(pm_well_type) :: this
  PetscReal, pointer :: residual(:)

  this%option%io_buffer = "PMWellModifyFlowResBase must be extended."
  call PrintErrMsg(this%option)

end subroutine PMWellModifyFlowResBase

! ************************************************************************** !

subroutine PMWellModifyFlowResHydrostatic(this,residual)
  !
  ! This subroutine computes the well contribution to the reservoir residual
  ! when called from the Hydrostatic well model.
  !
  ! Author: Michael Nole
  ! Date: 01/16/2023
  !

  use SCO2_Aux_module, only : sco2_thermal, SCO2_TEMPERATURE_DOF
  use Hydrate_Aux_module, only : HYDRATE_ENERGY_DOF, hydrate_fmw_comp

  class(pm_well_hydrostatic_type) :: this
  PetscReal, pointer :: residual(:)

  PetscReal :: Res(this%option%nflowdof)
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
    case(SCO2_MODE)
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
          sum_q = sum_q + sum(this%well%gas%q) + sum(this%well%liq%q)
          if (this%well%total_rate < 0.d0) then
            Q = -1.d0 * this%well%total_rate
          else
            Q = -1.d0 * (this%well%th_qg + this%well%th_ql)
          endif
          residual(local_end) = Q - sum_q
        endif
        do j = 0,2
          ! Compontent j+1 residual at well segment k
          residual(local_start + j) = residual(local_start + j) + &
                  (this%well%liq%q(k)* &
                  this%well%liq%xmass(k,j+1) + &
                  this%well%gas%q(k)* &
                  this%well%gas%xmass(k,j+1))
        enddo
        if (sco2_thermal) then
          ! Energy contribution
          j = SCO2_TEMPERATURE_DOF - 1
          residual(local_start + j) = residual(local_start + j) + &
                  (this%well%liq%q(k)* &
                  this%well%liq%H(k) + &
                  this%well%gas%q(k)* &
                  this%well%gas%H(k))
        endif
      enddo
    case(H_MODE)
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
        do j = 0,2
          ! Compontent j+1 residual at well segment k
          if (j == TWO_INTEGER) then
            ! Salt is offset by 1
            residual(local_start + j + 1) = &
                  residual(local_start + j + 1) + &
                  (this%well%liq%q(k)* &
                  this%well%liq%xmass(k,j+1) / &
                  hydrate_fmw_comp(j+1) + &
                  this%well%gas%q(k)* &
                  this%well%gas%xmass(k,j+1)) / &
                  hydrate_fmw_comp(j+1)
          else
            residual(local_start + j) = residual(local_start + j) + &
                  (this%well%liq%q(k)* &
                  this%well%liq%xmass(k,j+1) / &
                  hydrate_fmw_comp(j+1) + &
                  this%well%gas%q(k)* &
                  this%well%gas%xmass(k,j+1)) / &
                  hydrate_fmw_comp(j+1)
          endif
        enddo
        j = HYDRATE_ENERGY_DOF - 1
        residual(local_start + j) = residual(local_start + j) + &
                  (this%well%liq%q(k)* &
                  this%well%liq%H(k) + &
                  this%well%gas%q(k)* &
                  this%well%gas%H(k))
      enddo
  end select

end subroutine PMWellModifyFlowResHydrostatic

! ************************************************************************** !

subroutine PMWellModifyFlowJacBase(this,Jac,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/04/2025
  !

  implicit none

  class(pm_well_type) :: this
  Mat :: Jac
  PetscErrorCode :: ierr

  this%option%io_buffer = "PMWellModifyFlowJacBase must be extended."
  call PrintErrMsg(this%option)

end subroutine PMWellModifyFlowJacBase

! ************************************************************************** !

subroutine PMWellModifyFlowJacHydrostatic(this,Jac,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/04/2025
  !

  use Option_module
  use SCO2_Aux_module, only : SCO2_WELL_DOF, SCO2_TEMPERATURE_DOF, sco2_thermal
  use Hydrate_Aux_module, only : HYDRATE_WELL_DOF, HYDRATE_ENERGY_DOF, &
                                 hydrate_fmw_comp

  implicit none

  class(pm_well_hydrostatic_type) :: this
  Mat :: Jac
  PetscErrorCode :: ierr

  type(well_type), pointer :: well
  type(well_type), pointer :: well_pert
  type(option_type), pointer :: option
  PetscReal :: res
  PetscReal :: pert, res_pert
  PetscReal :: Q, sum_q
  PetscInt :: ghosted_id
  PetscReal, allocatable :: residual(:,:)
  PetscReal, allocatable :: J_block(:,:)
  PetscInt :: air_comp_id, wat_comp_id
  PetscInt :: j,k,idof,irow
  PetscInt :: local_row_index, local_col_index
  ! PetscInt :: global_row_index, global_col_index
  ! PetscReal :: dpl, dpg
  ! PetscReal :: global_id
  PetscReal :: J_well

  option => this%option

  air_comp_id = this%option%air_id
  wat_comp_id = this%option%water_id

  if (this%well_comm%comm == MPI_COMM_NULL) return

  select case (this%option%iflowmode)

    case(SCO2_MODE)

      allocate(J_block(option%nflowdof,option%nflowdof))
      allocate(residual(this%well_grid%nsegments,option%nflowdof))

      J_block = 0.d0
      residual = 0.d0
      ! Calculate Jacobian entries wrt reservoir perturbation:
      ! BHP residual and reservoir source/sink residuals.

      ! Unperturbed residuals
      well => this%well
      do k = 1,this%well_grid%nsegments
        do irow = 1, option%nflowspec
          residual(k,irow) = well%liq%q(k)* &
                              well%liq%xmass(k,irow) + &
                              well%gas%q(k)* &
                              well%gas%xmass(k,irow)
        enddo
        if (k==1) then
          sum_q = 0.d0
          sum_q = sum(well%gas%q) + sum(well%liq%q)
          if (well%total_rate < 0.d0) then
            Q = -1.d0 * well%total_rate
          else
            Q = -1.d0 * (well%th_qg + well%th_ql)
          endif
          residual(k,option%nflowdof) = Q - sum_q
        endif
        ! Energy residual
        if (sco2_thermal) then
          irow = SCO2_TEMPERATURE_DOF
          residual(k,irow) = well%liq%q(k)* &
                              well%liq%H(k) + &
                              well%gas%q(k)* &
                              well%gas%H(k)
        endif
      enddo

      ! Perturbations
      do k = 1,this%well_grid%nsegments
        J_block = 0.d0
        ghosted_id = this%well_grid%h_ghosted_id(k)
        do idof = 1,option%nflowdof-1
          well_pert => this%well_pert(idof)
          pert = this%realization%patch%aux%SCO2% &
                        auxvars(idof,ghosted_id)%pert
          ! Compute dRwell / dXres: only in the cell that owns the
          ! well bottom.
          if (this%well_grid%h_rank_id(1) == option%myrank) then
            res_pert = 0.d0
            res_pert = (well_pert%gas%q(k) + well_pert%liq%q(k)) - &
                        (well%gas%q(k) + well%liq%q(k))

            local_row_index = this%well_grid%h_ghosted_id(1)* &
                              option%nflowdof-1
            local_col_index = (ghosted_id-1)*option%nflowdof + idof - 1
            J_well = -(res_pert)/pert
            if (dabs(J_well) > 0.d0) then
              call MatSetValuesLocal(Jac,1,local_row_index,1, &
                                      local_col_index,J_well, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
            endif

          endif
          if (this%well_grid%h_rank_id(k) == option%myrank) then
            ! Compute dRres / dXres at a given BHP
            do irow = 1,option%nflowspec
              res_pert = well_pert%liq%q(k)* &
                          well_pert%liq%xmass(k,irow) + &
                          well_pert%gas%q(k)* &
                          well_pert%gas%xmass(k,irow)
              J_block(irow,idof) = (res_pert - residual(k,irow))/pert
            enddo
            if (sco2_thermal) then
              irow = SCO2_TEMPERATURE_DOF
              res_pert = well_pert%liq%q(k)* &
                          well_pert%liq%H(k) + &
                          well_pert%gas%q(k)* &
                          well_pert%gas%H(k)
              J_block(irow,idof) = (res_pert - residual(k,irow))/pert
            endif
          endif
        enddo
        if (any(dabs(J_block) > 0.d0)) then
          call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1,&
                                      J_block,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif

        !Perturbed rates wrt well perturbation.
        well_pert => this%well_pert(SCO2_WELL_DOF)

        if (this%well_grid%h_rank_id(k) /= option%myrank) cycle

        J_block = 0.d0

        ghosted_id = this%well_grid%h_ghosted_id(k)
        pert = well_pert%bh_p - well%bh_p

        ! Compute dRwell / dPwell
        if (k == 1) then
          sum_q = 0.d0
          sum_q = sum(well%gas%q) + sum(well%liq%q)
          if (well%total_rate < 0.d0) then
            Q = -1.d0 * well%total_rate
          else
            Q = -1.d0 * (well%th_qg + well%th_ql)
          endif
          res = Q - sum_q

          sum_q = 0.d0
          sum_q = sum(well_pert%gas%q) + sum(well_pert%liq%q)
          if (well_pert%total_rate < 0.d0) then
            Q = -1.d0 * well_pert%total_rate
          else
            Q = -1.d0 * (well_pert%th_qg + well_pert%th_ql)
          endif
          res_pert = Q - sum_q

          J_block(option%nflowdof,option%nflowdof) = &
                  (res_pert - res)/pert
        endif

        ! Compute dRres / dPwell
        do j = 0,option%nflowspec-1
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
          if (dabs(J_well) > 0.d0) then
            call MatSetValuesLocal(Jac,1,local_row_index,1,local_col_index, &
                                    J_well,ADD_VALUES,ierr);CHKERRQ(ierr)
          endif
        enddo
        if (sco2_thermal) then
          ! Energy contribution
          j = SCO2_TEMPERATURE_DOF - 1
          res = this%well%liq%q(k)*this%well%liq%H(k) + &
                this%well%gas%q(k)*this%well%gas%H(k)
          res_pert = well_pert%liq%q(k)* &
                      well_pert%liq%H(k) + &
                      well_pert%gas%q(k)* &
                      well_pert%gas%H(k)
          ! Just change 1 column
          local_row_index = (ghosted_id-1)*option%nflowdof + j
          local_col_index = this%well_grid%h_ghosted_id(1)* &
                            option%nflowdof-1
          J_well = (res_pert - res)/pert
          if (dabs(J_well) > 0.d0) then
            call MatSetValuesLocal(Jac,1,local_row_index,1,local_col_index, &
                                    J_well,ADD_VALUES,ierr);CHKERRQ(ierr)
          endif
        endif

        if (any(dabs(J_block) > 0.d0)) then
          call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1, &
                                      J_block,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif
      enddo
    case(H_MODE)

      allocate(J_block(option%nflowdof,option%nflowdof))
      allocate(residual(this%well_grid%nsegments,option%nflowdof))

      J_block = 0.d0
      residual = 0.d0

      ! Calculate Jacobian entries wrt reservoir perturbation:
      ! BHP residual and reservoir source/sink residuals.

      ! Unperturbed residual
      well => this%well
      do k = 1,this%well_grid%nsegments
        do irow = 1, option%nflowspec
          ! Salt is 3rd component, but 4th row
          if (irow == option%nflowspec) then
            residual(k,irow+1) = well%liq%q(k)* &
                              well%liq%xmass(k,irow) / &
                              hydrate_fmw_comp(irow) + &
                              well%gas%q(k)* &
                              well%gas%xmass(k,irow) / &
                              hydrate_fmw_comp(irow)
          else
            residual(k,irow) = well%liq%q(k)* &
                              well%liq%xmass(k,irow) / &
                              hydrate_fmw_comp(irow) + &
                              well%gas%q(k)* &
                              well%gas%xmass(k,irow) / &
                              hydrate_fmw_comp(irow)
          endif
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
        irow = HYDRATE_ENERGY_DOF
        residual(k,irow) = well%liq%q(k)* &
                              well%liq%H(k)+ &
                              well%gas%q(k)* &
                              well%gas%H(k)
      enddo

      ! Perturbations
      do k = 1,this%well_grid%nsegments
        J_block = 0.d0
        ghosted_id = this%well_grid%h_ghosted_id(k)
        do idof = 1,option%nflowdof-1
          well_pert => this%well_pert(idof)
          pert = this%realization%patch%aux%hydrate% &
                        auxvars(idof,ghosted_id)%pert

          ! Compute dRwell / dXres: only in the cell that owns the
          ! well bottom
          if (this%well_grid%h_rank_id(1) == option%myrank) then
            res_pert = 0.d0
            if (dabs(well_pert%th_qg) > 0.d0) then
              res_pert = well_pert%gas%q(k) - well%gas%q(k)
            else
              res_pert = well_pert%liq%q(k) - well%liq%q(k)
            endif

            local_row_index = this%well_grid%h_ghosted_id(1)* &
                              option%nflowdof-1
            local_col_index = (ghosted_id-1)*option%nflowdof + idof - 1
            J_well = -(res_pert)/pert
            if (dabs(J_well) > 0.d0) then
              call MatSetValuesLocal(Jac,1,local_row_index,1, &
                                      local_col_index,J_well, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
            endif

          endif

          if (this%well_grid%h_rank_id(k) == option%myrank) then
            ! Compute dRres / dXres at a given P_BHP
            do irow = 1,option%nflowspec
              res_pert = well_pert%liq%q(k)* &
                        well_pert%liq%xmass(k,irow) / &
                        hydrate_fmw_comp(irow) + &
                        well_pert%gas%q(k)* &
                        well_pert%gas%xmass(k,irow) / &
                        hydrate_fmw_comp(irow)
              ! Salt is 3rd component, but 4th row
              if (irow == option%nflowspec) then
                J_block(irow+1,idof) = (res_pert - residual(k,irow+1))/pert
              else
                J_block(irow,idof) = (res_pert - residual(k,irow))/pert
              endif
            enddo
            irow = HYDRATE_ENERGY_DOF
            res_pert = well_pert%liq%q(k)* &
                        well_pert%liq%H(k) + &
                        well_pert%gas%q(k)* &
                        well_pert%gas%H(k)
            J_block(irow,idof) = (res_pert - residual(k,irow))/pert
          endif
        enddo
        if (any(dabs(J_block) > 0.d0)) then
          call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1,&
                                      J_block,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif

        !Perturbed rates wrt well perturbation.
        well_pert => this%well_pert(HYDRATE_WELL_DOF)

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
        do j = 0,option%nflowspec-1
          res = this%well%liq%q(k)*this%well%liq%xmass(k,j+1) / &
                hydrate_fmw_comp(j+1) + &
                this%well%gas%q(k)*this%well%gas%xmass(k,j+1) / &
                hydrate_fmw_comp(j+1)
          res_pert = well_pert%liq%q(k)* &
                      well_pert%liq%xmass(k,j+1) / &
                      hydrate_fmw_comp(j+1) + &
                      well_pert%gas%q(k)* &
                      well_pert%gas%xmass(k,j+1) / &
                      hydrate_fmw_comp(j+1)
          ! Just change 1 column
          local_row_index = (ghosted_id-1)*option%nflowdof + j
          ! Salt is 3rd component, but 4th row
          if (j == option%nflowspec - 1) &
            local_row_index = local_row_index + 1
          local_col_index = this%well_grid%h_ghosted_id(1)* &
                            option%nflowdof-1
          J_well = (res_pert - res)/pert
          if (dabs(J_well) > 0.d0) then
            call MatSetValuesLocal(Jac,1,local_row_index,1,local_col_index, &
                                    J_well,ADD_VALUES,ierr);CHKERRQ(ierr)
          endif
        enddo
        ! Energy contribution
        j = HYDRATE_ENERGY_DOF - 1
        res = this%well%liq%q(k)*this%well%liq%H(k) + &
              this%well%gas%q(k)*this%well%gas%H(k)
        res_pert = well_pert%liq%q(k)* &
                    well_pert%liq%H(k) + &
                    well_pert%gas%q(k)* &
                    well_pert%gas%H(k)
        ! Just change 1 column
        local_row_index = (ghosted_id-1)*option%nflowdof + j
        local_col_index = this%well_grid%h_ghosted_id(1)* &
                          option%nflowdof-1
        J_well = (res_pert - res)/pert
        if (dabs(J_well) > 0.d0) then
          call MatSetValuesLocal(Jac,1,local_row_index,1,local_col_index, &
                                  J_well,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif

        if (any(dabs(J_block) > 0.d0)) then
          call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1, &
                                      J_block,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif
      enddo
  end select

  if (allocated(J_block)) deallocate(J_block)

end subroutine PMWellModifyFlowJacHydrostatic

! ************************************************************************** !

subroutine PMWellModifyDummyFlowJacobian(this,Jac,ierr)
  !
  ! This subroutine computes a dummy Jacobian for initializing a
  ! simulation with enough connectivity.
  !
  ! Author: Michael Nole
  ! Date: 04/12/2024
  !

  use Option_module

  implicit none

  class(pm_well_type) :: this
  Mat :: Jac
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  PetscInt :: ghosted_id
  PetscReal, allocatable :: J_block(:,:)
  PetscInt :: j,k,idof,irow
  PetscInt :: local_row_index, local_col_index
  PetscReal :: J_well

  option => this%option

  if (this%flow_coupling /= FULLY_IMPLICIT_WELL) return
  if (this%well_comm%comm == MPI_COMM_NULL) return

  select case (this%option%iflowmode)
    case(SCO2_MODE)

      allocate(J_block(option%nflowdof,option%nflowdof))

      J_block = 0.d0

      do k = 1,this%well_grid%nsegments
        J_block = 0.d0
        ghosted_id = this%well_grid%h_ghosted_id(k)
        do idof = 1,option%nflowdof-1
          ! Compute dRwell / dXres: only in the cell that owns the well
          if (this%well_grid%h_rank_id(1) == option%myrank) then
            local_row_index = this%well_grid%h_ghosted_id(1)* &
                              option%nflowdof-1
            local_col_index = (ghosted_id-1)*option%nflowdof + idof - 1
            J_well = UNINITIALIZED_DOUBLE
            if (dabs(J_well) > 0.d0) then
              call MatSetValuesLocal(Jac,1,local_row_index,1, &
                                      local_col_index,J_well, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
            endif
          endif
          if (this%well_grid%h_rank_id(k) == option%myrank) then
            ! Compute dRres / dXres at a given P_BHP
            do irow = 1,option%nflowspec
              J_block(irow,idof) = UNINITIALIZED_DOUBLE
            enddo
          endif
        enddo
        if (any(dabs(J_block) > 0.d0)) then
          call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1,&
                                      J_block,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif

        !Perturbed rates wrt well perturbation.

        if (this%well_grid%h_rank_id(k) /= option%myrank) cycle

        J_block = 0.d0

        ghosted_id = this%well_grid%h_ghosted_id(k)

        ! Compute dRwell / dPwell
        if (k == 1) then
          J_block(option%nflowdof,option%nflowdof) = UNINITIALIZED_DOUBLE
        endif

        ! Compute dRres / dPwell
        do j = 0,option%nflowspec-1
          ! Just change 1 column
          local_row_index = (ghosted_id-1)*option%nflowdof + j
          local_col_index = this%well_grid%h_ghosted_id(1)*option%nflowdof-1
          J_well = UNINITIALIZED_DOUBLE
          if (dabs(J_well) > 0.d0) then
            call MatSetValuesLocal(Jac,1,local_row_index,1,local_col_index, &
                                    J_well,ADD_VALUES,ierr);CHKERRQ(ierr)
          endif
        enddo
        if (option%nflowdof > 4) then
          j = 3
          local_row_index = (ghosted_id-1)*option%nflowdof + j
          local_col_index = this%well_grid%h_ghosted_id(1)* &
                              option%nflowdof-1
          J_well = UNINITIALIZED_DOUBLE
            if (dabs(J_well) > 0.d0) then
              call MatSetValuesLocal(Jac,1,local_row_index,1,local_col_index, &
                                      J_well,ADD_VALUES,ierr);CHKERRQ(ierr)
            endif
        endif

        if (any(dabs(J_block) > 0.d0)) then
          call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1, &
                                      J_block,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif
      enddo
    case(H_MODE)
      allocate(J_block(option%nflowdof,option%nflowdof))
      J_block = 0.d0
      ! Perturbations
      do k = 1,this%well_grid%nsegments
        J_block = 0.d0
        ghosted_id = this%well_grid%h_ghosted_id(k)
        do idof = 1,option%nflowdof-1
          ! Compute dRwell / dXres
          if (this%well_grid%h_rank_id(1) == option%myrank) then
            local_row_index = this%well_grid%h_ghosted_id(1)* &
                              option%nflowdof-1
            local_col_index = (ghosted_id-1)*option%nflowdof + idof - 1
            J_well = UNINITIALIZED_DOUBLE
            call MatSetValuesLocal(Jac,1,local_row_index,1, &
                                    local_col_index,J_well, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
          endif
          ! Compute dRres / dXres
          do irow = 1,option%nflowdof-1
            J_block(irow,idof) = UNINITIALIZED_DOUBLE
          enddo
        enddo
        if (any(dabs(J_block) > 0.d0)) then
          call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1,&
                                      J_block,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif
        !Perturbed rates wrt well perturbation.
        if (this%well_grid%h_rank_id(k) /= option%myrank) cycle
        J_block = 0.d0
        ghosted_id = this%well_grid%h_ghosted_id(k)
        ! Compute dRwell / dPwell
        if (k == 1) then
          J_block(option%nflowdof,option%nflowdof) = UNINITIALIZED_DOUBLE
        endif
        ! Compute dRres / dPwell
        do j = 0,option%nflowdof-2
          ! Just change 1 column
          local_row_index = (ghosted_id-1)*option%nflowdof + j
          local_col_index = this%well_grid%h_ghosted_id(1)*option%nflowdof-1
          J_well = UNINITIALIZED_DOUBLE
          call MatSetValuesLocal(Jac,1,local_row_index,1,local_col_index, &
                                  J_well,ADD_VALUES,ierr);CHKERRQ(ierr)
        enddo
        if (any(dabs(J_block) > 0.d0)) then
          call MatSetValuesBlockedLocal(Jac,1,ghosted_id-1,1,ghosted_id-1, &
                                        J_block,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif
      enddo
  end select

end subroutine PMWellModifyDummyFlowJacobian

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

subroutine PMWellSolveBase(this,time,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  implicit none

  class(pm_well_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  this%realization%option%io_buffer = 'PMWellSolveBase must be extended.'
  call PrintErrMsg(this%realization%option)

end subroutine PMWellSolveBase

! ************************************************************************** !

subroutine PMWellSolveSequential(this,time,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  implicit none

  class(pm_well_sequential_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  this%realization%option%io_buffer = 'PMWellSolveSequential must be extended.'
  call PrintErrMsg(this%realization%option)

end subroutine PMWellSolveSequential

! ************************************************************************** !

subroutine PMWellSolveQI(this,time,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  implicit none

  class(pm_well_qi_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  this%realization%option%io_buffer = 'PMWellSolveQI must be extended.'
  call PrintErrMsg(this%realization%option)

end subroutine PMWellSolveQI

! ************************************************************************** !

subroutine PMWellSolveImplicit(this,time,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  implicit none

  class(pm_well_implicit_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  this%realization%option%io_buffer = 'PMWellSolveImplicit must be extended.'
  call PrintErrMsg(this%realization%option)

end subroutine PMWellSolveImplicit

! ************************************************************************** !

subroutine PMWellSolveHydrostatic(this,time,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  implicit none

  class(pm_well_hydrostatic_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  ierr = 0
  call this%SolveFlow(UNINITIALIZED_INTEGER,ierr)
  call PMWellCalcCumulativeQFlux(this)

end subroutine PMWellSolveHydrostatic

! ************************************************************************** !

subroutine PMWellSolveFlowSequential(this,perturbation_index,ierr)
  !
  ! Author: Michael Nole
  ! Date: 12/01/2021
  !

  implicit none

  class(pm_well_sequential_type) :: this
  PetscInt :: perturbation_index
  PetscErrorCode :: ierr

  this%realization%option%io_buffer = 'PMWellSolveFlowSequential must be extended.'
  call PrintErrMsg(this%realization%option)

end subroutine PMWellSolveFlowSequential

! ************************************************************************** !

subroutine PMWellSolveFlowQI(this,perturbation_index,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  implicit none

  class(pm_well_qi_type) :: this
  PetscInt :: perturbation_index
  PetscErrorCode :: ierr

  this%realization%option%io_buffer = 'PMWellSolveFlowQI must be extended.'
  call PrintErrMsg(this%realization%option)

end subroutine PMWellSolveFlowQI

! ************************************************************************** !

subroutine PMWellSolveFlowBase(this,perturbation_index,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  implicit none

  class(pm_well_type) :: this
  PetscInt :: perturbation_index
  PetscErrorCode :: ierr

  this%realization%option%io_buffer = 'PMWellSolveFlowBase must be extended.'
  call PrintErrMsg(this%realization%option)

end subroutine PMWellSolveFlowBase

! ************************************************************************** !

subroutine PMWellSolveFlowHydrostatic(this,perturbation_index,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/03/2025
  !

  use Option_module
  use Grid_module
  use EOS_Water_module
  use EOS_Gas_module
  use SCO2_Aux_module, only : fmw_comp, SCO2BrineDensity
  use Hydrate_Aux_module, only : hydrate_fmw_comp

  implicit none

  class(pm_well_hydrostatic_type) :: this
  PetscInt :: perturbation_index
  PetscErrorCode :: ierr

  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: reservoir
  type(grid_type), pointer :: reservoir_grid
  type(option_type), pointer :: option
  PetscLogDouble :: log_start_time
  PetscReal :: Q_liq(this%well_grid%nsegments,this%well_grid%nsegments),&
               Q_gas(this%well_grid%nsegments,this%well_grid%nsegments)
  PetscReal :: v_darcy
  PetscBool :: upwind
  PetscReal :: mobility, enthalpy
  PetscReal :: area, mass_conserved_liq, mass_conserved_gas
  PetscInt :: i
  PetscReal :: pl, pl0, pg, pg0, temperature, den_mol
  PetscReal :: rho_kg_liq, rho_zero_liq, rho_one_liq, &
               rho_kg_gas, rho_zero_gas, rho_one_gas
  PetscReal :: dist_z, delta_z
  PetscReal :: gravity, den_ave
  PetscInt :: num_iteration
  PetscReal :: dummy, dummy2
  PetscReal :: aux(2)
  PetscReal :: res_pg_temp, res_pl_temp, res_z
  PetscReal :: fmw_comp_temp(3)
  PetscReal :: mixture_ratio
  PetscReal, parameter :: threshold_p = 0.d0
  PetscReal, parameter :: epsilon = 1.d-14

  option => this%realization%option
  reservoir_grid => this%realization%patch%grid

  if (this%well_comm%comm == MPI_COMM_NULL) return

  ierr = 0
  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

  gravity = option%gravity(Z_DIRECTION)

  ! update well index
  call PMWellComputeWellIndex(this)

  ! Compute hydrostatic pressure relative to bottom-hole pressure,
  ! then compute well fluxes from those pressures.

  select case (option%iflowmode)
    case (SCO2_MODE)
      fmw_comp_temp(:) = fmw_comp(:)
    case (H_MODE)
      fmw_comp_temp(:) = hydrate_fmw_comp(:)
  end select

  if (perturbation_index > 0) then
    well => this%well_pert(perturbation_index)
  else
    well => this%well
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
  select case (option%iflowmode)
    case (SCO2_MODE)
      call SCO2BrineDensity(temperature, pl, &
                            aux(1), rho_kg_liq, option)
    case default
      call EOSWaterDensityExt(temperature,pl, &
                                aux,rho_kg_liq,dummy,ierr)
  end select

  pl0 = pl

  call EOSGasDensity(temperature,pg, &
                      den_mol,dummy,dummy2,ierr)
  rho_kg_gas = den_mol * fmw_comp_temp(TWO_INTEGER)

  pg0 = pg

  ! compute pressures above datum
  dist_z = 0.d0
  rho_zero_liq = rho_kg_liq
  rho_zero_gas = rho_kg_gas

  do i = 1,this%well_grid%nsegments

    ! Compute well cell pressures based off of bottom segment pressure
    temperature = well%temp(i)
    aux(1) = well%liq%xmass(1,option%salt_id)
    select case (option%iflowmode)
      case (SCO2_MODE)
        call SCO2BrineDensity(temperature, pl, &
                              aux(1), rho_kg_liq, option)
      case default
        call EOSWaterDensityExt(temperature,pl, &
                                  aux,rho_kg_liq,dummy,ierr)
    end select
    call EOSGasDensity(temperature,pg0, &
                      den_mol,dummy,dummy2,ierr)
    rho_kg_gas = den_mol * fmw_comp_temp(TWO_INTEGER)

    num_iteration = 0
    if (i == 1) then
      delta_z = this%well_grid%h(i)%z - this%well_grid%bottomhole(3)
    else
      delta_z = this%well_grid%h(i)%z - this%well_grid%h(i-1)%z
    endif
    do
      pl = pl0 + 0.5d0*(rho_kg_liq+rho_zero_liq) * &
            gravity * delta_z
      pg = pg0 + 0.5d0*(rho_kg_gas+rho_zero_gas) * &
            gravity * delta_z
      call EOSWaterDensityExt(temperature,pl,aux,rho_one_liq,dummy,ierr)
      call EOSGasDensity(temperature,pg,den_mol,dummy,dummy2,ierr)
      rho_one_gas = den_mol * fmw_comp_temp(TWO_INTEGER)

      ! STOMP doesn't iterate
      exit

      if (dabs(rho_kg_liq-rho_one_liq) < 1.d-5 .and. &
          dabs(rho_kg_gas-rho_one_gas) < 1.d-5) exit
      rho_kg_liq = rho_one_liq
      rho_kg_gas = rho_one_gas
      num_iteration = num_iteration + 1
      if (num_iteration > 100) then
        option%io_buffer = 'Hydrostatic iteration failed to &
                            &converge in the well model'
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

    if (pg < 0.d0 .and. well%th_qg > 0.d0) then
      option%io_buffer = 'BHP cannot support the gas column: &
                          &hydrostatic well pressure calcuation &
                          &results in negative gas pressure.'
      call PrintErrMsg(option)
    endif
  enddo

  ! Update well properties based off new pressures
  select case(option%iflowmode)
    case(SCO2_MODE)
      call PMWellUpdatePropertiesSCO2Flow(this,well,option)
    case(H_MODE)
      call PMWellUpdatePropertiesHydrateFlow(this,well,option)
  end select

  ! Compute fluxes in/out of well
  well%liq%Q = 0.d0
  well%gas%Q = 0.d0

  do i = 1,this%well_grid%nsegments

    if (well%WI(i) == 0) cycle

    res_z = this%well_grid%res_z(i)
    delta_z = res_z - this%well_grid%h(i)%z

    if (well%th_qg > 0.d0 .and. well%th_ql > 0.d0) then
      ! mixture injection: assume constant mass ratio
      mixture_ratio = well%th_ql / well%th_qg

      if (reservoir%s_g(i) > epsilon) then
        res_pg_temp = reservoir%p_g(i) + reservoir%den_g(i) * gravity * &
                      delta_z
      else
        res_pg_temp = reservoir%p_g(i) + reservoir%den_l(i) * gravity * &
                      delta_z
      endif
      upwind = res_pg_temp > well%pg(i)
      if (upwind) then
        mobility = reservoir%kr_g(i)/reservoir%visc_g(i)
        den_ave = reservoir%den_g(i)
        enthalpy = reservoir%H_g(i)
      else
        mobility = 1.d0 / well%gas%visc(i)
        den_ave = well%gas%den(i)
        enthalpy = well%gas%H(i)
      endif

      ! Flowrate in kg/s
      well%gas%Q(i) = den_ave*mobility*well%WI(i)* &
                      (res_pg_temp-well%pg(i))

      ! Assuming the mass ratio of water and CO2 remains the same
      ! everywhere.
      well%liq%Q(i) = well%gas%Q(i) * mixture_ratio

    elseif (well%th_qg > 0.d0) then
      ! Rate-controlled gas injection well. Can potentially have
      ! under-pressure in some well segments.
      ! Compute reservoir pressure at well cell center
      if (reservoir%s_g(i) > epsilon) then
        res_pg_temp = reservoir%p_g(i) + reservoir%den_g(i) * gravity * delta_z
      else
        res_pg_temp = reservoir%p_g(i) + reservoir%den_l(i) * gravity * delta_z
      endif
      upwind = res_pg_temp > well%pg(i)
      if (upwind) then
        mobility = reservoir%kr_g(i)/reservoir%visc_g(i)
        den_ave = reservoir%den_g(i)
        enthalpy = reservoir%H_g(i)
      else
        mobility = 1.d0 / well%gas%visc(i)
        den_ave = well%gas%den(i)
        enthalpy = well%gas%H(i)
      endif

      ! Flowrate in kg/s
      well%gas%Q(i) = den_ave*mobility*well%WI(i)* &
                      (res_pg_temp-well%pg(i))
      ! if (well%gas%Q(i) < 0.d0) well%gas%Q(i) = 0.d0
    elseif (well%th_ql > 0.d0) then
      ! Rate-controlled water injection well
      ! Compute reservoir pressure at well cell center
      res_pl_temp = reservoir%p_l(i) + reservoir%den_l(i) * gravity * delta_z
      upwind = res_pl_temp > well%pl(i)
      if (upwind) then
        mobility = reservoir%kr_l(i)/reservoir%visc_l(i)
        den_ave = reservoir%den_l(i)
        enthalpy = reservoir%H_l(i)
      else
        mobility = 1.d0 / well%liq%visc(i)
        den_ave = well%liq%den(i)
        enthalpy = well%liq%H(i)
      endif
      ! Flowrate in kg/s
      well%liq%Q(i) = den_ave*mobility*well%WI(i)* &
                      (res_pl_temp-well%pl(i))
      ! if (well%liq%Q(i) > 0.d0) well%liq%Q(i) = 0.d0
    elseif (well%total_rate < 0.d0) then
      ! Extraction well
      ! Compute reservoir pressure at well cell center
        well%liq%Q(i) = 0.d0
        well%gas%Q(i) = 0.d0
        res_pl_temp = reservoir%p_l(i) + reservoir%den_l(i) * gravity * delta_z
        mobility = reservoir%kr_l(i)/reservoir%visc_l(i)
        den_ave = reservoir%den_l(i)
        ! Flowrate in kg/s
        if (res_pl_temp > well%pl(i)) then
          well%liq%Q(i) = den_ave*mobility*well%WI(i) * &
                          (res_pl_temp-well%pl(i)) * &
                          reservoir%xmass_liq(i,ONE_INTEGER)
          well%gas%Q(i) = den_ave*mobility*well%WI(i)* &
                          (res_pl_temp-well%pl(i)) * &
                          reservoir%xmass_liq(i,TWO_INTEGER)
        endif

        res_pg_temp = reservoir%p_g(i) + reservoir%den_g(i) * gravity * delta_z
        mobility = reservoir%kr_g(i)/reservoir%visc_g(i)
        den_ave = reservoir%den_g(i)
        ! Flowrate in kg/s
        if (res_pg_temp > well%pg(i)) then
          well%liq%Q(i) = well%liq%Q(i) + den_ave*mobility*well%WI(i) * &
                          (res_pg_temp-well%pg(i)) * &
                          reservoir%xmass_gas(i,ONE_INTEGER)
          well%gas%Q(i) = well%gas%Q(i) + den_ave*mobility*well%WI(i) * &
                          (res_pg_temp-well%pg(i)) * &
                          reservoir%xmass_gas(i,TWO_INTEGER)
        endif
    endif
  enddo

  ! Should equal the total flux out the top of the domain
  mass_conserved_liq = sum(well%liq%Q)
  mass_conserved_gas = sum(well%gas%Q)

end subroutine PMWellSolveFlowHydrostatic

! ! ************************************************************************** !

subroutine PMWellSolveFlowImplicit(this,perturbation_index,ierr)
  !
  ! Author: Michael Nole
  ! Date: 12/01/2021
  !

  implicit none

  class(pm_well_implicit_type) :: this
  PetscInt :: perturbation_index
  PetscErrorCode :: ierr

  this%realization%option%io_buffer = 'PMWellSolveFlowImplicit must be &
                                       &extended.'
  call PrintErrMsg(this%realization%option)

end subroutine PMWellSolveFlowImplicit

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
        if (pm_well%well_grid%casing(k) <= 0.d0) then
          pm_well%well%WI(k) = 0.d0
          cycle
        endif
        temp_real = log(2.079d-1*reservoir%dx(k)/ &
                        (pm_well%well%diameter(k)/2.d0)) + pm_well%well%skin(k)

        if (temp_real <= 0.d0) then
          write(diameter_string,'(F7.4)') pm_well%well%diameter(k)
          write(dx_string,'(F7.4)') reservoir%dx(k)
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
        if (pm_well%well_grid%casing(k) <= 0.d0) then
          pm_well%well%WI(k) = 0.d0
          cycle
        endif
        r0 = 2.8d-1*(sqrt(sqrt(reservoir%ky(k)/reservoir%kx(k))* &
             reservoir%dx(k)**2 + sqrt(reservoir%kx(k)/ &
             reservoir%ky(k))*reservoir%dy(k)**2) / &
             ((reservoir%ky(k)/reservoir%kx(k))**2.5d-1 + &
             (reservoir%kx(k)/reservoir%ky(k))**2.5d-1))

        temp_real = log(r0/(pm_well%well%diameter(k)/2.d0)) + &
                    pm_well%well%skin(k)

        if (temp_real <= 0.d0) then
          write(diameter_string,'(F7.4)') pm_well%well%diameter(k)
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
        if (associated(pm_well%well_grid%dx)) then
          dh_x = pm_well%well_grid%dx(k)
          dh_y = pm_well%well_grid%dy(k)
          dh_z = pm_well%well_grid%dz(k)
        else
          ! Infer the segment lengths from vertical well.
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
        endif

        if (pm_well%well_grid%casing(k) <= 0.d0) then
          pm_well%well%WI(k) = 0.d0
          cycle
        endif

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
          write(diameter_string,'(F7.4)') pm_well%well%diameter(k)
          option%io_buffer = 'Wellbore diameter (' // diameter_string // ' m)&
          & is too large relative to reservoir discretization and &
          &permeability for the default anisotropic PEACEMAN_3D well model.'
          call PrintErrMsg(option)
        endif

        pm_well%well%WI(k) = sqrt((wix**2) + (wiy**2) + (wiz**2))
      enddo

    case(PEACEMAN_NONE)
      do k = 1,pm_well%well_grid%nsegments
        ! Assume a vertical well
        ! Assume connection between segment and reservoir has height
        ! this%well_grid%dh(k)
        pm_well%well%WI(k) = sqrt(reservoir%kx(k)*reservoir%ky(k)) * &
                             pm_well%well_grid%dh(k)
      enddo
  end select

  pm_well%well%WI = pm_well%well%WI*pm_well%well_grid%casing

end subroutine PMWellComputeWellIndex

! ************************************************************************** !

subroutine PMWellUpdatePropertiesSCO2Flow(pm_well,well,option)
  !
  ! Updates flow well object properties, when SCO2 mode is the flow mode.
  !
  ! Author: Michael Nole
  ! Date: 01/06/2024
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
               den_kg_co2, den_steam
  PetscReal :: visc_co2, visc_water, visc_liq, &
               visc_brine, visc_gas
  PetscReal :: xco2g, xwg, xco2l,xwl, &
               xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl
  PetscInt :: wid, co2_id, sid
  PetscReal :: H_temp, U_temp, den_co2, H_steam
  PetscErrorCode :: ierr

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id

  nsegments =pm_well%well_grid%nsegments

  if (well%total_rate < 0.d0) then
    ! Extraction well: use reservoir fluid properties
    well%liq%xmass(:,:) = well%reservoir%xmass_liq(:,:)
    well%gas%xmass(:,:) = well%reservoir%xmass_gas(:,:)
    well%liq%den(:) = well%reservoir%den_l(:)
    well%gas%den(:) = well%reservoir%den_g(:)
    well%liq%visc(:) = well%reservoir%visc_l(:)
    well%gas%visc(:) = well%reservoir%visc_g(:)
    well%liq%H(:) = well%reservoir%h_l(:)
    well%gas%H(:) = well%reservoir%h_g(:)
  else
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
                                  xsl, Prvap,PETSC_FALSE)
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

      if (sco2_thermal) then
        ! Energy calculations

        ! Brine enthalpy
        call EOSWaterEnthalpy(well%temp(i),well%pg(i), &
                              well%liq%H(i),ierr)
        ! J/kmol --> J/kg
        well%liq%H(i) = well%liq%H(i) / fmw_comp(wid)
        call SCO2BrineEnthalpy(well%temp(i), well%liq%xmass(i,sid), &
                              well%liq%H(i),H_temp)
        ! CO2 density, internal energy, enthalpy
        call EOSGasDensityEnergy(well%temp(i),well%pg(i),den_co2, &
                                well%gas%H(i),U_temp,ierr)
        ! J/kmol --> J/kg
        well%gas%H(i) = well%gas%H(i) / fmw_comp(co2_id)

        ! Liquid phase enthalpy
        well%liq%H(i) = SCO2EnthalpyCompositeLiquid(well%temp(i), &
                                      well%liq%xmass(i,sid), &
                                      well%liq%xmass(i,co2_id), &
                                      H_temp, well%gas%H(i))

        well%liq%H(i) = well%liq%H(i) * 1.d-6 ! J/kg -> MJ/kg
        well%gas%H(i) = well%gas%H(i)  * 1.d-6 ! MJ/kg
        call EOSWaterSteamDensityEnthalpy(well%temp(i), &
                                          well%pg(i), &
                                          den_kg_steam, &
                                          den_steam, &
                                          H_steam,ierr)
        ! J/kmol -> MJ/kg
        H_steam = H_steam / fmw_comp(wid) * 1.d-6
      else
        den_steam = 0.d0
        H_steam = 0.d0
      endif

      ! Gas phase enthalpy
      well%gas%H(i) = well%gas%xmass(i,wid) * H_steam + &
                      well%gas%xmass(i,co2_id) * well%gas%H(i)

    enddo
  endif

end subroutine PMWellUpdatePropertiesSCO2Flow

! ************************************************************************** !

subroutine PMWellUpdatePropertiesHydrateFlow(pm_well,well,option)
  !
  ! Updates flow well object properties, when Hydrate mode is the flow mode.
  ! Can only inject CO2 right now.
  !
  ! Author: Michael Nole
  ! Date: 06/06/2024
  !

  use EOS_Water_module
  use EOS_Gas_module
  use Hydrate_Aux_module

  implicit none

  class(pm_well_type) :: pm_well
  type(well_type) :: well
  type(option_type) :: option

  PetscInt :: i,nsegments
  PetscReal :: drho_dT,drho_dP
  PetscReal :: xsl, Pa, Pvap, Pva, Ps, Prvap
  PetscReal :: cell_pressure
  PetscReal :: den_kg_water, den_kg_steam, &
               den_kg_brine, den_kg_liq, &
               den_kg_gas, den_mol_a, &
               den_kg_a, den_steam
  PetscReal :: visc_a, visc_water, visc_liq, &
               visc_brine, visc_gas
  PetscReal :: xag, xwg, xal,xwl, &
               xmolag, xmolwg, xmolal, xmolsl, xmolwl
  PetscInt :: wid, acid, sid
  PetscReal :: H_temp, U_temp, den_a, H_steam
  PetscErrorCode :: ierr

  if (pm_well%well_comm%comm == MPI_COMM_NULL) return

  wid = option%water_id
  acid = option%air_id
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
    call HydrateBrineSaturationPressure(well%temp(i), &
                                     xsl,Ps)
    call HydrateBrineDensity(well%temp(i), well%pg(i), &
                          xsl, den_kg_brine, option)
    call HydrateVaporPressureBrine(well%temp(i), Ps, &
                                0.d0, den_kg_brine, &
                                xsl, Prvap)
    call HydrateWaterDensity(well%temp(i),Prvap, &
                          TWO_INTEGER,den_kg_water, &
                          den_kg_steam,option)
    call HydrateDensityCompositeLiquid(well%temp(i),den_kg_brine, &
                                  well%liq%xmass(i,acid), &
                                  den_kg_liq)

    well%liq%den(i) = den_kg_liq

    call HydrateEquilibrate(well%temp(i),well%pg(i), FIVE_INTEGER, &
                         0.d0, Pa, Pvap, Ps, Prvap, &
                         xag, xwg, xal, xsl, xwl, &
                         xmolag, xmolwg, xmolal, xmolsl, xmolwl, option)

    xmolag = (well%gas%xmass(i,acid)/hydrate_fmw_comp(2)) / &
               ((well%gas%xmass(i,acid)/hydrate_fmw_comp(2)) + &
                (well%gas%xmass(i,wid)/hydrate_fmw_comp(1)))
    xmolwg = 1.d0 - xmolag
    xmolal = (well%liq%xmass(i,acid)/hydrate_fmw_comp(2)) / &
               ((well%gas%xmass(i,acid)/hydrate_fmw_comp(2)) + &
               (well%gas%xmass(i,wid)/hydrate_fmw_comp(1)) + &
               (well%gas%xmass(i,sid)/hydrate_fmw_comp(3)))
    xmolwl = (well%liq%xmass(i,wid)/hydrate_fmw_comp(2)) / &
               ((well%gas%xmass(i,acid)/hydrate_fmw_comp(2)) + &
               (well%gas%xmass(i,wid)/hydrate_fmw_comp(1)) + &
               (well%gas%xmass(i,sid)/hydrate_fmw_comp(3)))
    xmolsl = 1.d0 - xmolal - xmolwl

    !Gas Density
    Pva = max(well%pg(i),Prvap)
    call EOSGasDensity(well%temp(i),Pva, &
                       den_mol_a,drho_dT,drho_dP,ierr)
    den_kg_a = den_mol_a * hydrate_fmw_comp(2)
    den_kg_gas = well%gas%xmass(i,option%air_id) * &
                 den_kg_a + &
                 well%gas%xmass(i,option%water_id) * &
                 den_kg_steam
    well%gas%den(i) = den_kg_gas

    ! Liquid Viscosity
    call HydrateViscosityWater(well%temp(i),well%pg(i), &
                           den_kg_water,visc_water,option)
    call HydrateViscosityCO2(well%temp(i), den_kg_a, &
                          visc_a)
    call HydrateViscosityBrine(well%temp(i), xsl, &
                           visc_water, visc_brine)
    call HydrateViscosityLiquid(xmolal, visc_brine, &
                             visc_a, visc_liq)

    well%liq%visc(i) = visc_liq

    ! Gas Viscosity
    call HydrateViscosityGas(visc_water,visc_a,xmolwg, &
                          xmolag,visc_gas)

    well%gas%visc(i) = visc_gas

    ! Energy calculations

    ! Brine enthalpy
    call EOSWaterEnthalpy(well%temp(i),well%pg(i), &
                          well%liq%H(i),ierr)
    ! J/kmol --> J/kg
    well%liq%H(i) = well%liq%H(i) / hydrate_fmw_comp(wid)
    call HydrateBrineEnthalpy(well%temp(i), well%liq%xmass(i,sid), &
                           well%liq%H(i),H_temp)
    ! CO2 density, internal energy, enthalpy
    call EOSGasDensityEnergy(well%temp(i),well%pg(i),den_a, &
                             well%gas%H(i),U_temp,ierr)
    ! J/kmol --> J/kg
    well%gas%H(i) = well%gas%H(i) / hydrate_fmw_comp(acid)

    ! Liquid phase enthalpy
    well%liq%H(i) = HydrateEnthalpyCompositeLiquid(well%temp(i), &
                                   well%liq%xmass(i,sid), &
                                   well%liq%xmass(i,acid), &
                                   H_temp, well%gas%H(i))

    well%liq%H(i) = well%liq%H(i) * 1.d-6 ! J/kg -> MJ/kg
    well%gas%H(i) = well%gas%H(i)  * 1.d-6 ! MJ/kg
    cell_pressure = min(Ps, Prvap)
    call EOSWaterSteamDensityEnthalpy(well%temp(i), &
                                      cell_pressure, &
                                      den_kg_steam, &
                                      den_steam, &
                                      H_steam,ierr)
    ! J/kmol -> MJ/kg
    H_steam = H_steam / hydrate_fmw_comp(wid) * 1.d-6

    ! Gas phase enthalpy
    well%gas%H(i) = well%gas%xmass(i,wid) * H_steam + &
                    well%gas%xmass(i,acid) * well%gas%H(i)
  enddo

end subroutine PMWellUpdatePropertiesHydrateFlow

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
  well_copy%total_rate = well%total_rate
  well_copy%liq%visc(:) = well%liq%visc(:)
  well_copy%gas%visc(:) = well%gas%visc(:)
  well_copy%temp(:) = well%temp(:)
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
  reservoir_copy%H_l(:) = reservoir%H_l(:)
  reservoir_copy%H_g(:) = reservoir%H_g(:)
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

subroutine PMWellSetPlotVariablesBase(list,pm_well)
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

  select case (pm_well%option%iflowmode)
    case (WF_MODE)
      units = 'kmol/sec'
    case default
      units = 'kg/sec'
  end select
  if (pm_well%well%liq%output_Q) then
    name = 'Well Liq. Q'
    call OutputVariableAddToList(list,name,OUTPUT_RATE,units,WELL_LIQ_Q)
  endif
  if (pm_well%well%gas%output_Q) then
    name = 'Well Gas Q'
    call OutputVariableAddToList(list,name,OUTPUT_RATE,units,WELL_GAS_Q)
  endif

end subroutine PMWellSetPlotVariablesBase

! ************************************************************************** !

subroutine PMWellSetPlotVariablesSequential(list,pm_well)
  !
  ! Adds variables to be printed for plotting.
  !
  ! Author: Jenn Frederick
  ! Date: 09/29/2022
  !
  use Output_Aux_module
  use Variables_module

  type(output_variable_list_type), pointer :: list
  class(pm_well_sequential_type) :: pm_well

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

  select case (pm_well%option%iflowmode)
    case (WF_MODE)
      units = 'kmol/sec'
    case default
      units = 'kg/sec'
  end select
  if (pm_well%well%liq%output_Q) then
    name = 'Well Liq. Q'
    call OutputVariableAddToList(list,name,OUTPUT_RATE,units,WELL_LIQ_Q)
  endif
  if (pm_well%well%gas%output_Q) then
    name = 'Well Gas Q'
    call OutputVariableAddToList(list,name,OUTPUT_RATE,units,WELL_GAS_Q)
  endif
  select case (pm_well%option%iflowmode)
    case (SCO2_MODE, H_MODE)
      name = 'Well BHP'
      units = 'Pa'
      call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units,WELL_BHP)
  end select

end subroutine PMWellSetPlotVariablesSequential

! ************************************************************************** !

subroutine PMWellSetPlotVariablesHydrostatic(list,pm_well)
  !
  ! Adds variables to be printed for plotting.
  !
  ! Author: Jenn Frederick
  ! Date: 09/29/2022
  !
  use Output_Aux_module
  use Variables_module

  type(output_variable_list_type), pointer :: list
  class(pm_well_hydrostatic_type) :: pm_well

  character(len=MAXWORDLENGTH) :: name,  units

  call PMWellSetPlotVariablesBase(list,pm_well)

  name = 'Well BHP'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units,WELL_BHP)

end subroutine PMWellSetPlotVariablesHydrostatic

! ************************************************************************** !

function PMWellOutputFilename(name,option)
  !
  ! Generates a filename for wellbore model output
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021
  !

  implicit none

  character(len=MAXWORDLENGTH) :: name
  type(option_type), pointer :: option

  character(len=MAXSTRINGLENGTH) :: PMWellOutputFilename

  PMWellOutputFilename = trim(name) // '.well'

end function PMWellOutputFilename

! ************************************************************************** !

subroutine PMWellOutputHeaderSequential(pm_well)
  !
  ! Writes the header for wellbore model output file.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021, updated 01/09/2025
  !

  use Output_Aux_module
  use Grid_module
  use Utility_module

  implicit none

  class(pm_well_sequential_type) :: pm_well

  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename,word
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: cell_string
  PetscBool :: exist
  PetscInt :: fid
  PetscInt :: icolumn
  PetscInt :: k, j, i

  if (pm_well%option%myrank /= pm_well%well_grid%h_rank_id(1)) return

  output_option => pm_well%realization%output_option

  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif

  do k = 1,pm_well%well_grid%nsegments
    ! only print the segments that are being requested in the .well file
    if (.not. pm_well%well%output_in_well_file(k)) cycle

    fid = k + 100
    filename = PMWellOutputFilename(pm_well%name,pm_well%option)
    write(word,'(i5)') k
    if (k < 10000) word(1:1) = '0'; if (k < 1000) word(1:2) = '00'
    if (k < 100) word(1:3) = '000'; if (k < 10) word(1:4) = '0000'
    if (pm_well%split_output_file) then
      filename = trim(filename) // trim(word)
    endif

    exist = FileExists(trim(filename))
    if (pm_well%option%restart_flag .and. exist) return
    open(unit=fid,file=filename,action="write",status="replace")

    ! First write out the well grid information
    write(fid,'(a)',advance="yes") '========= WELLBORE MODEL GRID INFORMATION &
                                    &======================='
    write(word,'(i5)') pm_well%well_grid%nsegments
    write(fid,'(a)',advance="yes") ' Number of segments: ' // trim(word)
    write(word,'(i5)') pm_well%well_grid%nconnections
    write(fid,'(a)',advance="yes") ' Number of connections: ' // trim(word)
    write(word,'(es10.3,es10.3,es10.3)') pm_well%well_grid%tophole(1), &
                            pm_well%well_grid%tophole(2), pm_well%well_grid% &
                            tophole(3)
    write(fid,'(a)',advance="yes") ' Top of hole (x,y,z) [m]: ' // trim(word)
    write(word,'(es10.3,es10.3,es10.3)') pm_well%well_grid%bottomhole(1), &
                      pm_well%well_grid%bottomhole(2), pm_well%well_grid% &
                      bottomhole(3)
    write(fid,'(a)',advance="yes") ' Bottom of hole (x,y,z) [m]: ' // trim(word)
    write(fid,'(a)',advance="yes") '===========================================&
                                    &======================'
    write(fid,'(a)',advance="yes") ' Segment Number: Center coordinate (x,y,z) [m] '
    do j = 1,pm_well%well_grid%nsegments
      write(word,'(i4,a3,es10.3,es10.3,es10.3,a1)') j,': (', &
        pm_well%well_grid%h(j)%x,pm_well%well_grid%h(j)%y,pm_well% &
        well_grid%h(j)%z,')'
      write(fid,'(a)',advance="yes") trim(word)
    enddo
    write(fid,'(a)',advance="yes") '===========================================&
                                    &======================'
    write(fid,'(a)',advance="yes") ' Segment Number: Segment length [m] '
    do j = 1,pm_well%well_grid%nsegments
      write(word,'(i4,a3,es10.3,a1)') j,': (',pm_well%well_grid%dh(j),')'
      write(fid,'(a)',advance="yes") trim(word)
    enddo
    write(fid,'(a)',advance="yes") '===========================================&
                                    &======================'
    write(fid,'(a)',advance="no") ' Segment Numbers Requested for Output: '
    if (associated(pm_well%well%segments_for_output)) then
      do j = 1,size(pm_well%well%segments_for_output)
        write(word,'(i4)') pm_well%well%segments_for_output(j)
        write(fid,'(a)',advance="no") trim(word)
      enddo
      write(fid,'(a)',advance="yes") ' '
    else
      write(fid,'(a)',advance="yes") 'ALL SEGMENTS'
    endif
    write(fid,'(a)',advance="yes") '===========================================&
                                    &======================'

    close(fid)
    if (.not. pm_well%split_output_file) exit
  enddo

  do k = 1,pm_well%well_grid%nsegments
    ! only print the segments that are being requested in the .well file
    if (.not. pm_well%well%output_in_well_file(k)) cycle

    fid = k + 100
    filename = PMWellOutputFilename(pm_well%name,pm_well%option)
    write(word,'(i5)') k
    if (k < 10000) word(1:1) = '0'; if (k < 1000) word(1:2) = '00'
    if (k < 100) word(1:3) = '000'; if (k < 10) word(1:4) = '0000'
    if (pm_well%split_output_file) then
      filename = trim(filename) // trim(word)
    endif
    open(unit=fid,file=filename,action="write",status="old", &
         position="append")

    write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'
    cell_string = ''

    ! here a j while loop that depends on if the file is split or not
    j = k
    do while (j <= pm_well%well_grid%nsegments)
      if (.not. pm_well%well%output_in_well_file(j)) then
        j = j + 1
        cycle
      endif

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
      units_string = 'kmol/s'
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      variable_string = 'Well Q-gas'
      units_string = 'kmol/s'
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      variable_string = 'Well Q-liq-cumu'
      units_string = 'kmol'
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
        do i = 1,pm_well%nspecies
          variable_string = 'Well Aqueous Conc. ' // &
                            trim(pm_well%well%species_names(i))
          units_string = 'mol/m^3-liq'
          call OutputWriteToHeader(fid,variable_string,units_string, &
                                   cell_string,icolumn)
          variable_string = 'Well Aqueous Mass. ' &
                             // trim(pm_well%well%species_names(i))
          units_string = 'mol'
          call OutputWriteToHeader(fid,variable_string,units_string, &
                                   cell_string,icolumn)
          variable_string = 'Well q-Cumu Aqueous Mass. ' &
                             // trim(pm_well%well%species_names(i))
          units_string = 'mol'
          call OutputWriteToHeader(fid,variable_string,units_string, &
                                   cell_string,icolumn)
          variable_string = 'Well Q-Cumu Aqueous Mass. ' &
                             // trim(pm_well%well%species_names(i))
          units_string = 'mol'
          call OutputWriteToHeader(fid,variable_string,units_string, &
                                   cell_string,icolumn)
          variable_string = 'Res Aqueous Conc. ' // &
                            trim(pm_well%well%species_names(i))
          units_string = 'mol/m^3-liq'
          call OutputWriteToHeader(fid,variable_string,units_string, &
                                   cell_string,icolumn)
        enddo
      endif

      if (pm_well%split_output_file) exit
      j = j + 1
    enddo

    close(fid)
    if (.not. pm_well%split_output_file) exit
  enddo

end subroutine PMWellOutputHeaderSequential

! ************************************************************************** !

subroutine PMWellOutputHeaderHydrostatic(pm_well)
  !
  ! Writes the header for wellbore model output file.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021, updated 01/09/2025
  !

  use Output_Aux_module
  use Grid_module
  use Utility_module

  implicit none

  class(pm_well_hydrostatic_type) :: pm_well

  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename,word
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: cell_string
  PetscBool :: exist
  PetscInt :: fid
  PetscInt :: icolumn
  PetscInt :: k, j

  if (pm_well%option%myrank /= pm_well%well_grid%h_rank_id(1)) return

  output_option => pm_well%realization%output_option

  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif

  do k = 1,pm_well%well_grid%nsegments
    ! only print the segments that are being requested in the .well file
    if (.not. pm_well%well%output_in_well_file(k)) cycle

    fid = k + 100
    filename = PMWellOutputFilename(pm_well%name,pm_well%option)
    write(word,'(i5)') k
    if (k < 10000) word(1:1) = '0'; if (k < 1000) word(1:2) = '00'
    if (k < 100) word(1:3) = '000'; if (k < 10) word(1:4) = '0000'
    if (pm_well%split_output_file) then
      filename = trim(filename) // trim(word)
    endif

    exist = FileExists(trim(filename))
    if (pm_well%option%restart_flag .and. exist) return
    open(unit=fid,file=filename,action="write",status="replace")

    ! First write out the well grid information
    write(fid,'(a)',advance="yes") '========= WELLBORE MODEL GRID INFORMATION &
                                    &======================='
    write(word,'(i5)') pm_well%well_grid%nsegments
    write(fid,'(a)',advance="yes") ' Number of segments: ' // trim(word)
    write(word,'(i5)') pm_well%well_grid%nconnections
    write(fid,'(a)',advance="yes") ' Number of connections: ' // trim(word)
    write(word,'(es10.3,es10.3,es10.3)') pm_well%well_grid%tophole(1), &
                            pm_well%well_grid%tophole(2), pm_well%well_grid% &
                            tophole(3)
    write(fid,'(a)',advance="yes") ' Top of hole (x,y,z) [m]: ' // trim(word)
    write(word,'(es10.3,es10.3,es10.3)') pm_well%well_grid%bottomhole(1), &
                      pm_well%well_grid%bottomhole(2), pm_well%well_grid% &
                      bottomhole(3)
    write(fid,'(a)',advance="yes") ' Bottom of hole (x,y,z) [m]: ' // trim(word)
    write(fid,'(a)',advance="yes") '===========================================&
                                    &======================'
    write(fid,'(a)',advance="yes") ' Segment Number: Center coordinate (x,y,z) [m] '
    do j = 1,pm_well%well_grid%nsegments
      write(word,'(i4,a3,es10.3,es10.3,es10.3,a1)') j,': (', &
        pm_well%well_grid%h(j)%x,pm_well%well_grid%h(j)%y,pm_well% &
        well_grid%h(j)%z,')'
      write(fid,'(a)',advance="yes") trim(word)
    enddo
    write(fid,'(a)',advance="yes") '===========================================&
                                    &======================'
    write(fid,'(a)',advance="yes") ' Segment Number: Segment length [m] '
    do j = 1,pm_well%well_grid%nsegments
      write(word,'(i4,a3,es10.3,a1)') j,': (',pm_well%well_grid%dh(j),')'
      write(fid,'(a)',advance="yes") trim(word)
    enddo
    write(fid,'(a)',advance="yes") '===========================================&
                                    &======================'
    write(fid,'(a)',advance="no") ' Segment Numbers Requested for Output: '
    if (associated(pm_well%well%segments_for_output)) then
      do j = 1,size(pm_well%well%segments_for_output)
        write(word,'(i4)') pm_well%well%segments_for_output(j)
        write(fid,'(a)',advance="no") trim(word)
      enddo
      write(fid,'(a)',advance="yes") ' '
    else
      write(fid,'(a)',advance="yes") 'ALL SEGMENTS'
    endif
    write(fid,'(a)',advance="yes") '===========================================&
                                    &======================'

    close(fid)
    if (.not. pm_well%split_output_file) exit
  enddo

  do k = 1,pm_well%well_grid%nsegments
    ! only print the segments that are being requested in the .well file
    if (.not. pm_well%well%output_in_well_file(k)) cycle

    fid = k + 100
    filename = PMWellOutputFilename(pm_well%name,pm_well%option)
    write(word,'(i5)') k
    if (k < 10000) word(1:1) = '0'; if (k < 1000) word(1:2) = '00'
    if (k < 100) word(1:3) = '000'; if (k < 10) word(1:4) = '0000'
    if (pm_well%split_output_file) then
      filename = trim(filename) // trim(word)
    endif
    open(unit=fid,file=filename,action="write",status="old", &
         position="append")

    write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'
    cell_string = ''

    variable_string = 'Well BHP'
    units_string = 'Pa'
    call OutputWriteToHeader(fid,variable_string,units_string, &
                              cell_string,icolumn)

    ! here a j while loop that depends on if the file is split or not
    j = k
    do while (j <= pm_well%well_grid%nsegments)
      if (.not. pm_well%well%output_in_well_file(j)) then
        j = j + 1
        cycle
      endif

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
      units_string = 'kmol/s'
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      variable_string = 'Well Q-gas'
      units_string = 'kmol/s'
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      variable_string = 'Well Q-liq-cumu'
      units_string = 'kmol'
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

      if (pm_well%split_output_file) exit
      j = j + 1
    enddo

    close(fid)
    if (.not. pm_well%split_output_file) exit
  enddo

end subroutine PMWellOutputHeaderHydrostatic

! ************************************************************************** !

subroutine PMWellOutputSequential(pm_well)
  !
  ! Sets up .well file output for the wellbore process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021, updated 01/09/2025
  !

  use Output_Aux_module
  use Global_Aux_module
  use Grid_module

  implicit none

  class(pm_well_sequential_type) :: pm_well

  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename,word
  PetscInt :: fid
  PetscInt :: k, j, i

100 format(100es18.8)
101 format(1I6.1)

  if (pm_well%option%myrank /= pm_well%well_grid%h_rank_id(1)) return

  option => pm_well%realization%option
  output_option => pm_well%realization%output_option

  do k = 1,pm_well%well_grid%nsegments
    ! only print the segments that are being requested in the .well file
    if (.not. pm_well%well%output_in_well_file(k)) cycle

    fid = k + 100
    filename = PMWellOutputFilename(pm_well%name,pm_well%option)
    write(word,'(i5)') k
    if (k < 10000) word(1:1) = '0'; if (k < 1000) word(1:2) = '00'
    if (k < 100) word(1:3) = '000'; if (k < 10) word(1:4) = '0000'
    if (pm_well%split_output_file) then
      filename = trim(filename) // trim(word)
    endif
    open(unit=fid,file=filename,action="write",status="old", &
         position="append")

    ! pm_well time is set at the end of the wellbore step????
    write(fid,100,advance="no") option%time / output_option%tconv

    ! here a j while loop that depends on if the file is split or not
    j = k
    do while (j <= pm_well%well_grid%nsegments)
      if (.not. pm_well%well%output_in_well_file(j)) then
        j = j + 1
        cycle
      endif

      write(fid,101,advance="no") j
      write(fid,100,advance="no") pm_well%well_grid%h(j)%x, &
                                  pm_well%well_grid%h(j)%y, &
                                  pm_well%well_grid%h(j)%z, &
                                  pm_well%well%pl(j), &
                                  pm_well%well%reservoir%p_l(j), &
                                  pm_well%well%pg(j), &
                                  pm_well%well%reservoir%p_g(j), &
                                  pm_well%well%liq%s(j), &
                                  pm_well%well%gas%s(j), &
                                  pm_well%well%liq%Q(j), &
                                  pm_well%well%gas%Q(j), &
                                  pm_well%well%liq%Q_cumulative(j)
      if (j == 1) then
        write(fid,100,advance="no") pm_well%well%ql_bc(1), &
                                    pm_well%well%qg_bc(1)
      endif
      if (j > 1) then
         write(fid,100,advance="no") pm_well%well%ql(j-1), &
                                     pm_well%well%qg(j-1)
      endif
      write(fid,100,advance="no") pm_well%well%mass_balance_liq(j), &
                                  pm_well%well%liq_mass(j), &
                                  pm_well%well%WI(j)
      if (pm_well%transport) then
        do i = 1,pm_well%nspecies
          write(fid,100,advance="no") pm_well%well%aqueous_conc(i,j), &
                                pm_well%well%aqueous_mass(i,j), &
                                pm_well%well%aqueous_mass_q_cumulative(i,j), &
                                pm_well%well%aqueous_mass_Qcumulative(i,j), &
                                pm_well%well%reservoir%aqueous_conc(i,j)
        enddo
      endif

      if (pm_well%split_output_file) exit
      j = j + 1
    enddo

    close(fid)
    if (.not. pm_well%split_output_file) exit
  enddo

end subroutine PMWellOutputSequential

! ************************************************************************** !

subroutine PMWellOutputHydrostatic(pm_well)
  !
  ! Sets up .well file output for the wellbore process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021, updated 01/09/2025
  !

  use Output_Aux_module
  use Global_Aux_module
  use Grid_module

  implicit none

  class(pm_well_hydrostatic_type) :: pm_well

  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename,word
  PetscInt :: fid
  PetscInt :: k, j

100 format(100es18.8)
101 format(1I6.1)

  if (pm_well%option%myrank /= pm_well%well_grid%h_rank_id(1)) return

  option => pm_well%realization%option
  output_option => pm_well%realization%output_option

  do k = 1,pm_well%well_grid%nsegments
    ! only print the segments that are being requested in the .well file
    if (.not. pm_well%well%output_in_well_file(k)) cycle

    fid = k + 100
    filename = PMWellOutputFilename(pm_well%name,pm_well%option)
    write(word,'(i5)') k
    if (k < 10000) word(1:1) = '0'; if (k < 1000) word(1:2) = '00'
    if (k < 100) word(1:3) = '000'; if (k < 10) word(1:4) = '0000'
    if (pm_well%split_output_file) then
      filename = trim(filename) // trim(word)
    endif
    open(unit=fid,file=filename,action="write",status="old", &
         position="append")

    ! pm_well time is set at the end of the wellbore step????
    write(fid,100,advance="no") option%time / output_option%tconv

    write(fid,100,advance="no") pm_well%well%bh_p

    ! here a j while loop that depends on if the file is split or not
    j = k
    do while (j <= pm_well%well_grid%nsegments)
      if (.not. pm_well%well%output_in_well_file(j)) then
        j = j + 1
        cycle
      endif

      write(fid,101,advance="no") j
      write(fid,100,advance="no") pm_well%well_grid%h(j)%x, &
                                  pm_well%well_grid%h(j)%y, &
                                  pm_well%well_grid%h(j)%z, &
                                  pm_well%well%pl(j), &
                                  pm_well%well%reservoir%p_l(j), &
                                  pm_well%well%pg(j), &
                                  pm_well%well%reservoir%p_g(j), &
                                  pm_well%well%liq%s(j), &
                                  pm_well%well%gas%s(j), &
                                  pm_well%well%liq%Q(j), &
                                  pm_well%well%gas%Q(j), &
                                  pm_well%well%liq%Q_cumulative(j)
      if (j == 1) then
        write(fid,100,advance="no") pm_well%well%ql_bc(1), &
                                    pm_well%well%qg_bc(1)
      endif
      if (j > 1) then
         write(fid,100,advance="no") pm_well%well%ql(j-1), &
                                     pm_well%well%qg(j-1)
      endif
      write(fid,100,advance="no") pm_well%well%mass_balance_liq(j), &
                                  pm_well%well%liq_mass(j), &
                                  pm_well%well%WI(j)

      if (pm_well%split_output_file) exit
      j = j + 1
    enddo

    close(fid)
    if (.not. pm_well%split_output_file) exit
  enddo

end subroutine PMWellOutputHydrostatic

! ************************************************************************** !

subroutine PMWellCalcCumulativeQFlux(pm_well)
  !
  ! Calculates the cumulative Q flux in the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 01/10/2025
  !

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: k,nsegments
  PetscReal :: dt

  nsegments = pm_well%well_grid%nsegments

  dt = pm_well%option%flow_dt

  ! (+) Q is into the well [kmol-liq/sec]
  ! (-) Q is out of the well [kmol-liq/sec]

  do k = 1,nsegments
    ! Q_cumulative [kmol-liq]
    pm_well%well%liq%Q_cumulative(k) = pm_well%well%liq%Q_cumulative(k) + &
                                       (pm_well%well%liq%Q(k)*dt)
  enddo

end subroutine PMWellCalcCumulativeQFlux

! ************************************************************************** !

subroutine PMWellCalcCumulativeTranFlux(pm_well)
  !
  ! Calculates the cumulative flux of species into and out of the well.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 01/16/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: resr
  PetscReal :: coef_Qin,coef_Qout ! into well, out of well
  PetscReal :: Qin,Qout
  PetscInt :: k,nsegments
  PetscInt :: i,nspecies
  PetscReal :: dt,den_avg,area,q_liq,conc
  PetscReal :: diffusion

  nsegments = pm_well%well_grid%nsegments
  nspecies = pm_well%nspecies

  well => pm_well%well
  resr => pm_well%well%reservoir

  dt = pm_well%option%flow_dt

  ! (+) Q is into the well [kmol-liq/sec]
  ! (-) Q is out of the well [kmol-liq/sec]
  ! (+) aqueous_mass_Qcumulative = net [mol-species] is into well
  ! (-) aqueous_mass_Qcumulative = net [mol-species] is out of well

  diffusion = 0.d0 ! for now, since WIPP has no diffusion

  do k = 1,nsegments
    do i = 1,nspecies
      if (k == nsegments) then
        area = pm_well%well%area(k) ! [m2]
        q_liq = pm_well%well%ql_bc(2) ! [m3-liq/m2-bulk-sec]=[m/s]
      else
        area = 0.5d0*(pm_well%well%area(k)+pm_well%well%area(k+1))
        q_liq = pm_well%well%ql(k) ! [m3-liq/m2-bulk-sec]=[m/s]
      endif

      if (q_liq < 0.d0) then ! flow is down the well
        if (k == nsegments) then
          conc = pm_well%well%aqueous_conc_th(i) ! [mol-species/m3-liq]
        else
          conc = pm_well%tran_soln%prev_soln%aqueous_conc(i,k+1) ! [mol-species/m3-liq]
        endif
      elseif (q_liq > 0.d0) then ! flow is up the well
        conc = pm_well%tran_soln%prev_soln%aqueous_conc(i,k) ! [mol-species/m3-liq]
      else ! q=0
        conc = 0.d0 ! [mol-species/m3-liq]
      endif

      well%aqueous_mass_q_cumulative(i,k) = &
        ! [mol-species]                      + [m2*[m/s]*[mol/m3]*s]
        well%aqueous_mass_q_cumulative(i,k) + (dt*area*(q_liq*conc-diffusion))
    enddo

    den_avg = 0.5d0*(well%liq%den(k)+resr%den_l(k))
    ! units of coef = [m^3-liq/sec]
    if (well%liq%Q(k) < 0.d0) then ! Q out of well
      coef_Qin = 0.d0
      coef_Qout = well%liq%Q(k)*FMWH2O/den_avg
    else ! Q into well
    !            [kmol-liq/sec]*[kg-liq/kmol-liq]/[kg-liq/m^3-liq]
      coef_Qin = well%liq%Q(k)*FMWH2O/den_avg
      coef_Qout = 0.d0
    endif
    do i = 1,nspecies
      !     [m^3-liq/sec]*[mol-species/m^3-liq] = [mol-species/sec]
      Qin = coef_Qin*pm_well%tran_soln%prev_soln%resr_aqueous_conc(i,k)
      Qout = coef_Qout*pm_well%tran_soln%prev_soln%aqueous_conc(i,k)
      well%aqueous_mass_Qcumulative(i,k) = &
      ! [mol-species]                      + [mol-species/sec]*[sec]
        well%aqueous_mass_Qcumulative(i,k) + ((Qin+Qout)*dt)
    enddo
  enddo

end subroutine PMWellCalcCumulativeTranFlux

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
  PetscInt :: isegment,nsegments
  PetscInt :: n_up,n_dn
  PetscReal :: mass_rate_up,mass_rate_dn

  ! q in [m3-liq/m2-bulk-sec]
  ! area in [m2-bulk]

  ! (-) mass balance rate means net mass is being lost
  ! (+) mass balance rate means net mass is being gained

  well => pm_well%well
  nsegments = pm_well%well_grid%nsegments

  if (pm_well%well_comm%comm == MPI_COMM_NULL .or. &
      pm_well%option%iflowmode /= WF_MODE) then
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

  if (pm_well%well_comm%comm == MPI_COMM_NULL .or. &
        pm_well%option%iflowmode /= WF_MODE) then
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

subroutine PMWellDestroyBase(this)
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
  call DeallocateArray(this%well%reservoir%H_l)
  call DeallocateArray(this%well%reservoir%H_g)
  call DeallocateArray(this%well%reservoir%e_por)
  call DeallocateArray(this%well%reservoir%kx)
  call DeallocateArray(this%well%reservoir%ky)
  call DeallocateArray(this%well%reservoir%kz)
  call DeallocateArray(this%well%reservoir%dx)
  call DeallocateArray(this%well%reservoir%dy)
  call DeallocateArray(this%well%reservoir%dz)
  call DeallocateArray(this%well%reservoir%volume)
  call DeallocateArray(this%well%reservoir%tmp_flow)

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
  call DeallocateArray(this%well%output_in_well_file)
  call DeallocateArray(this%well%segments_for_output)
  call DeallocateArray(this%well%liq%den)
  call DeallocateArray(this%well%liq%visc)
  call DeallocateArray(this%well%liq%s)
  call DeallocateArray(this%well%liq%Q)
  call DeallocateArray(this%well%liq%Q_cumulative)
  call DeallocateArray(this%well%gas%den)
  call DeallocateArray(this%well%gas%visc)
  call DeallocateArray(this%well%gas%s)
  call DeallocateArray(this%well%gas%Q)
  nullify(this%well%liq)
  nullify(this%well%gas)

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
    call DeallocateArray(this%well_pert(i)%reservoir%H_l)
    call DeallocateArray(this%well_pert(i)%reservoir%H_g)
    call DeallocateArray(this%well_pert(i)%reservoir%e_por)
    call DeallocateArray(this%well_pert(i)%reservoir%kx)
    call DeallocateArray(this%well_pert(i)%reservoir%ky)
    call DeallocateArray(this%well_pert(i)%reservoir%kz)
    call DeallocateArray(this%well_pert(i)%reservoir%dx)
    call DeallocateArray(this%well_pert(i)%reservoir%dy)
    call DeallocateArray(this%well_pert(i)%reservoir%dz)
    call DeallocateArray(this%well_pert(i)%reservoir%volume)
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
  enddo

end subroutine PMWellDestroyBase

! ************************************************************************** !

subroutine PMWellDestroySequential(this)
  !
  ! Destroys the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_well_sequential_type) :: this

  PetscInt :: s, i

  call PMWellDestroyBase(this)

  if (this%transport) then
    call DeallocateArray(this%well%reservoir%aqueous_conc)
    call DeallocateArray(this%well%reservoir%aqueous_mass)
    call DeallocateArray(this%well%reservoir%tmp_tran)
  endif
  nullify(this%well%reservoir)

  if (this%transport) then
    call DeallocateArray(this%well%reservoir_save%aqueous_conc)
    call DeallocateArray(this%well%reservoir_save%aqueous_mass)
  endif
  nullify(this%well%reservoir_save)

  if (this%transport) then
    call DeallocateArray(this%well%aqueous_conc)
    call DeallocateArray(this%well%aqueous_mass)
    call DeallocateArray(this%well%aqueous_mass_q_cumulative)
    call DeallocateArray(this%well%aqueous_mass_Qcumulative)
  endif
  nullify(this%well)

  s = size(this%well_pert)

  do i = 1,s
    if (this%transport) then
      call DeallocateArray(this%well_pert(i)%reservoir%aqueous_conc)
      call DeallocateArray(this%well_pert(i)%reservoir%aqueous_mass)
    endif
    nullify(this%well_pert(i)%reservoir)
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

end subroutine PMWellDestroySequential

! ************************************************************************** !

subroutine PMWellDestroyHydrostatic(this)
  !
  ! Destroys the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !

  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_well_hydrostatic_type) :: this

  PetscInt :: s, i

  call PMWellDestroyBase(this)

  nullify(this%well%reservoir)

  nullify(this%well%reservoir_save)

  nullify(this%well)

  s = size(this%well_pert)

  do i = 1,s
    nullify(this%well_pert(i)%reservoir)
    nullify(this%well_pert(i)%reservoir_save)
  enddo
  nullify(this%well_pert)

  nullify(this%strata_list)

end subroutine PMWellDestroyHydrostatic

! ************************************************************************** !

end module PM_Well_class
