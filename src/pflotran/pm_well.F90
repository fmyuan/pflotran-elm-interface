module PM_Well_class

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscsnes.h"
  use petscsys
  use petscsnes
  use PM_Base_class
  use Option_module
  use Geometry_module
  use Strata_module
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module
  use WIPP_Flow_Aux_module
  use NW_Transport_Aux_module

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


  PetscBool :: initialize_well = PETSC_TRUE
  PetscReal, parameter :: gravity = -9.80665d0 ! [m/s2]

  PetscInt, parameter :: PEACEMAN_ISO = 1
  PetscInt, parameter :: PEACEMAN_ANISOTROPIC = 2

  type :: well_grid_type
    ! number of well segments
    PetscInt :: nsegments
    ! number of well connections
    PetscInt :: nconnections
    ! delta h discretization of each segment center [m]
    PetscReal, pointer :: dh(:)
    ! reservoir dz
    PetscReal, pointer :: res_dz(:)
    ! h coordinate of each segment center [m]
    type(point3d_type), pointer :: h(:)
    ! the local id of the reservoir grid cell within which each segment
    ! center resides
    PetscInt, pointer :: h_local_id(:)
    ! the ghosted id of the reservoir grid cell within which each segment
    ! center resides
    PetscInt, pointer :: h_ghosted_id(:)
    ! the global id of the reservoir grid cell within which each segment
    ! center resides
    PetscInt, pointer :: h_global_id(:)
    ! the strata id number associated with each well segment
    PetscInt, pointer :: strata_id(:)
    ! coordinate of the top/bottom of the well [m]
    PetscReal :: tophole(3)
    PetscReal :: bottomhole(3)
    ! x-y span search distance multipier
    PetscInt :: xy_span_multiplier
  end type well_grid_type

  type :: well_reservoir_type
    ! reservoir liquid pressure [Pa]
    PetscReal, pointer :: p_l(:)
    ! reservoir gas pressure [Pa]
    PetscReal, pointer :: p_g(:)
    ! reservoir liquid saturation [-]
    PetscReal, pointer :: s_l(:)
    ! reservoir gas saturation [-]
    PetscReal, pointer :: s_g(:)
    ! reservoir liquid mobility
    PetscReal, pointer :: mobility_l(:)
    ! reservoir gas mobility
    PetscReal, pointer :: mobility_g(:)
    ! reservoir liquid relative permeability [-]
    PetscReal, pointer :: kr_l(:)
    ! reservoir gas relative permeability [-]
    PetscReal, pointer :: kr_g(:)
    ! reservoir liquid density [kg/m3]
    PetscReal, pointer :: rho_l(:)
    ! reservoir gas density [kg/m3]
    PetscReal, pointer :: rho_g(:)
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
    ! type of well model to use
    character(len=MAXWORDLENGTH) :: well_model_type
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
    PetscReal, pointer :: f(:)
    ! well index of each well segment [0,1]  0 = cased; 1 = open
    PetscReal, pointer :: WI_base(:)
    ! total well index for each well segment (including reservoir effects)
    PetscReal, pointer :: WI(:)
    ! Peaceman equivalent radius
    PetscReal, pointer :: r0(:)
    ! well index model (probably has to get moved out of well_type)
    PetscInt :: WI_model = PEACEMAN_ISO
    ! well liquid pressure [Pa]
    PetscReal, pointer :: pl(:)
    ! well gas pressure [Pa]
    PetscReal, pointer :: pg(:)
    ! flag for output
    PetscBool :: output_pl 
    ! flag for output
    PetscBool :: output_pg 
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
    PetscReal :: rho0
    ! fluid density [kg/m3]
    PetscReal, pointer :: rho(:)
    ! fluid viscosity [Pa-s]
    PetscReal, pointer :: visc(:)
    ! fluid saturation
    PetscReal, pointer :: s(:)
    ! fluid source/sink in/out of well [kmol/s]
    PetscReal, pointer :: Q(:)
    ! flag for output
    PetscBool :: output_Q
  end type well_fluid_type

  ! primary variables necessary to reset flow solution
  type :: well_flow_save_type
    PetscReal, pointer :: pl(:)
    PetscReal, pointer :: sg(:)
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
  end type well_soln_tran_type

  type, public, extends(pm_base_type) :: pm_well_type
    class(realization_subsurface_type), pointer :: realization
    type(well_grid_type), pointer :: well_grid
    type(well_type), pointer :: well
    type(well_type), pointer :: well_pert(:)
    type(well_reservoir_type), pointer :: reservoir
    type(well_soln_flow_type), pointer :: flow_soln
    type(well_soln_tran_type), pointer :: tran_soln
    type(strata_list_type), pointer :: strata_list
    PetscReal, pointer :: pert(:,:)
    PetscInt :: nphase
    PetscInt :: nspecies
    PetscReal :: intrusion_time_start           ! [sec]
    PetscReal :: bh_zero_value                  ! [mol/m3-bulk]
    PetscReal :: dt_flow, dt_tran               ! [sec]
    PetscReal :: min_dt_flow, min_dt_tran       ! [sec]
    PetscReal :: cumulative_dt_flow             ! [sec]
    PetscReal :: cumulative_dt_tran             ! [sec]
    PetscReal :: min_dz
    PetscInt :: well_res_ratio
    PetscBool :: transport
    PetscBool :: ss_check
    PetscBool :: print_well
  contains
    procedure, public :: Setup => PMWellSetup
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
    procedure, public :: Destroy => PMWellDestroy
  end type pm_well_type

  public :: PMWellCreate, PMWellReadPMBlock, PMWellReadPass2

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

  allocate(PMWellCreate)
  call PMBaseInit(PMWellCreate)

  PMWellCreate%header = 'WELLBORE MODEL'

  nullify(PMWellCreate%realization)
  PMWellCreate%min_dt_flow = 1.d-15
  PMWellCreate%min_dt_tran = 1.d-15
  PMWellCreate%bh_zero_value = 1.d-40
  PMWellCreate%intrusion_time_start = UNINITIALIZED_DOUBLE
  PMWellCreate%nphase = 0
  PMWellCreate%nspecies = 0
  PMWellCreate%transport = PETSC_FALSE
  PMWellCreate%ss_check = PETSC_FALSE
  PMWellCreate%print_well = PETSC_FALSE
  PMWellCreate%min_dz = UNINITIALIZED_DOUBLE
  PMWellCreate%well_res_ratio = 1

  nullify(PMWellCreate%pert)

  ! create the well grid object:
  allocate(PMWellCreate%well_grid)
  PMWellCreate%well_grid%nsegments = UNINITIALIZED_INTEGER
  PMWellCreate%well_grid%nconnections = UNINITIALIZED_INTEGER
  nullify(PMWellCreate%well_grid%dh)
  nullify(PMWellCreate%well_grid%h)
  nullify(PMWellCreate%well_grid%h_local_id)
  nullify(PMWellCreate%well_grid%h_ghosted_id)
  nullify(PMWellCreate%well_grid%h_global_id)
  nullify(PMWellCreate%well_grid%strata_id)
  PMWellCreate%well_grid%tophole(:) = UNINITIALIZED_DOUBLE
  PMWellCreate%well_grid%bottomhole(:) = UNINITIALIZED_DOUBLE
  PMWellCreate%well_grid%xy_span_multiplier = 10

  ! create the well object:
  allocate(PMWellCreate%well)
  nullify(PMWellCreate%well%mass_balance_liq)
  nullify(PMWellCreate%well%area)
  nullify(PMWellCreate%well%diameter)
  nullify(PMWellCreate%well%volume)
  nullify(PMWellCreate%well%f)
  nullify(PMWellCreate%well%WI_base)
  nullify(PMWellCreate%well%WI)
  nullify(PMWellCreate%well%r0)
  nullify(PMWellCreate%well%pl)
  nullify(PMWellCreate%well%pg)
  PMWellCreate%well%output_pl = PETSC_FALSE
  PMWellCreate%well%output_pg = PETSC_FALSE
  nullify(PMWellCreate%well%ql)
  nullify(PMWellCreate%well%qg)
  nullify(PMWellCreate%well%ql_bc)
  nullify(PMWellCreate%well%qg_bc)
  nullify(PMWellCreate%well%ql_kmol)
  nullify(PMWellCreate%well%qg_kmol)
  nullify(PMWellCreate%well%ql_kmol_bc)
  nullify(PMWellCreate%well%qg_kmol_bc)
  nullify(PMWellCreate%well%liq_cum_mass)
  nullify(PMWellCreate%well%liq_mass)
  nullify(PMWellCreate%well%species_names)
  nullify(PMWellCreate%well%species_parent_id)
  nullify(PMWellCreate%well%species_radioactive)
  nullify(PMWellCreate%well%species_decay_rate)
  nullify(PMWellCreate%well%species_parent_decay_rate)
  nullify(PMWellCreate%well%aqueous_conc)
  nullify(PMWellCreate%well%aqueous_mass)
  PMWellCreate%well%output_aqc = PETSC_FALSE
  PMWellCreate%well%output_aqm = PETSC_FALSE
  nullify(PMWellCreate%well%aqueous_conc_th)
  nullify(PMWellCreate%well%ccid)
  nullify(PMWellCreate%well%permeability)
  nullify(PMWellCreate%well%phi)
  PMWellCreate%well%bh_p_set_by_reservoir = PETSC_FALSE
  PMWellCreate%well%bh_sg_set_by_reservoir = PETSC_FALSE
  PMWellCreate%well%bh_p = UNINITIALIZED_DOUBLE
  PMWellCreate%well%th_p = UNINITIALIZED_DOUBLE
  PMWellCreate%well%bh_sg = UNINITIALIZED_DOUBLE
  PMWellCreate%well%th_sg = UNINITIALIZED_DOUBLE
  PMWellCreate%well%bh_ql = UNINITIALIZED_DOUBLE
  PMWellCreate%well%bh_qg = UNINITIALIZED_DOUBLE
  PMWellCreate%well%th_ql = UNINITIALIZED_DOUBLE
  PMWellCreate%well%th_qg = UNINITIALIZED_DOUBLE
  PMWellCreate%well%tran_condition_name = ''

  ! create the well_pert object:
  allocate(PMWellCreate%well_pert(TWO_INTEGER))
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%mass_balance_liq)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%area)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%diameter)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%volume)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%f)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%WI_base)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%WI)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%r0)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%pl)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%pg)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%ql)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%qg)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%ql_bc)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%qg_bc)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%ql_kmol)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%qg_kmol)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%ql_kmol_bc)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%qg_kmol_bc)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%liq_cum_mass)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%liq_mass)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%aqueous_conc)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%aqueous_mass)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%ccid)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%permeability)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%phi)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%mass_balance_liq)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%area)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%diameter)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%volume)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%f)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%WI_base)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%WI)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%r0)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%pl)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%pg)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%ql)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%qg)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%ql_bc)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%qg_bc)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%ql_kmol)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%qg_kmol)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%ql_kmol_bc)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%qg_kmol_bc)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%liq_cum_mass)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%liq_mass)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%aqueous_conc)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%aqueous_mass)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%ccid)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%permeability)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%phi)
  PMWellCreate%well_pert(:)%bh_p_set_by_reservoir = PETSC_FALSE
  PMWellCreate%well_pert(:)%bh_sg_set_by_reservoir = PETSC_FALSE
  PMWellCreate%well_pert(:)%bh_p = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(:)%th_p = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(:)%bh_sg = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(:)%th_sg = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(:)%bh_ql = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(:)%bh_qg = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(:)%th_ql = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(:)%th_qg = UNINITIALIZED_DOUBLE

  ! create the reservoir object:
  allocate(PMWellCreate%reservoir)
  nullify(PMWellCreate%reservoir%p_l)
  nullify(PMWellCreate%reservoir%p_g)
  nullify(PMWellCreate%reservoir%s_l)
  nullify(PMWellCreate%reservoir%s_g)
  nullify(PMWellCreate%reservoir%mobility_l)
  nullify(PMWellCreate%reservoir%mobility_g)
  nullify(PMWellCreate%reservoir%kr_l)
  nullify(PMWellCreate%reservoir%kr_g)
  nullify(PMWellCreate%reservoir%rho_l)
  nullify(PMWellCreate%reservoir%rho_g)
  nullify(PMWellCreate%reservoir%visc_l)
  nullify(PMWellCreate%reservoir%visc_g)
  nullify(PMWellCreate%reservoir%e_por)
  nullify(PMWellCreate%reservoir%aqueous_conc)
  nullify(PMWellCreate%reservoir%aqueous_mass)
  nullify(PMWellCreate%reservoir%kx)
  nullify(PMWellCreate%reservoir%ky)
  nullify(PMWellCreate%reservoir%kz)
  nullify(PMWellCreate%reservoir%dx)
  nullify(PMWellCreate%reservoir%dy)
  nullify(PMWellCreate%reservoir%dz)
  nullify(PMWellCreate%reservoir%volume)


  ! create the fluid/liq objects:
  allocate(PMWellCreate%well%liq)
  PMWellCreate%well%liq%ifluid = 1
  nullify(PMWellCreate%well%liq%kr)
  PMWellCreate%well%liq%rho0 = UNINITIALIZED_DOUBLE
  nullify(PMWellCreate%well%liq%rho)
  nullify(PMWellCreate%well%liq%visc)
  nullify(PMWellCreate%well%liq%s)
  nullify(PMWellCreate%well%liq%Q)
  PMWellCreate%well%liq%output_Q = PETSC_FALSE
  allocate(PMWellCreate%well_pert(ONE_INTEGER)%liq)
  allocate(PMWellCreate%well_pert(TWO_INTEGER)%liq)
  PMWellCreate%well_pert(ONE_INTEGER)%liq%ifluid = 1
  PMWellCreate%well_pert(TWO_INTEGER)%liq%ifluid = 1
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%liq%kr)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%liq%kr)
  PMWellCreate%well_pert(ONE_INTEGER)%liq%rho0 = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(TWO_INTEGER)%liq%rho0 = UNINITIALIZED_DOUBLE

  ! create the fluid/gas objects:
  allocate(PMWellCreate%well%gas)
  PMWellCreate%well%gas%ifluid = 2
  nullify(PMWellCreate%well%gas%kr)
  PMWellCreate%well%gas%rho0 = UNINITIALIZED_DOUBLE
  nullify(PMWellCreate%well%gas%rho)
  nullify(PMWellCreate%well%gas%visc)
  nullify(PMWellCreate%well%gas%s)
  nullify(PMWellCreate%well%gas%Q)
  PMWellCreate%well%gas%output_Q = PETSC_FALSE
  allocate(PMWellCreate%well_pert(ONE_INTEGER)%gas)
  allocate(PMWellCreate%well_pert(TWO_INTEGER)%gas)
  PMWellCreate%well_pert(ONE_INTEGER)%gas%ifluid = 2
  PMWellCreate%well_pert(TWO_INTEGER)%gas%ifluid = 2
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%gas%kr)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%gas%kr)
  PMWellCreate%well_pert(ONE_INTEGER)%gas%rho0 = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(TWO_INTEGER)%gas%rho0 = UNINITIALIZED_DOUBLE

  ! create the well flow solution object:
  allocate(PMWellCreate%flow_soln)
  nullify(PMWellCreate%flow_soln%prev_soln%pl)
  nullify(PMWellCreate%flow_soln%prev_soln%sg)
  nullify(PMWellCreate%flow_soln%residual)
  nullify(PMWellCreate%flow_soln%Jacobian)
  nullify(PMWellCreate%flow_soln%update)
  PMWellCreate%flow_soln%ndof = UNINITIALIZED_INTEGER
  PMWellCreate%flow_soln%bh_p = PETSC_FALSE
  PMWellCreate%flow_soln%th_p = PETSC_FALSE
  PMWellCreate%flow_soln%bh_q = PETSC_FALSE
  PMWellCreate%flow_soln%th_q = PETSC_FALSE
  PMWellCreate%flow_soln%not_converged = PETSC_TRUE
  PMWellCreate%flow_soln%converged = PETSC_FALSE
  PMWellCreate%flow_soln%cut_timestep = PETSC_FALSE
  PMWellCreate%flow_soln%max_iter = 8
  PMWellCreate%flow_soln%max_ts_cut = 20
  PMWellCreate%flow_soln%ts_cut_factor = 2.0d0
  PMWellCreate%flow_soln%ts_ramp_factor = 1.1d0
  PMWellCreate%flow_soln%itol_abs_res = 1.0d-8
  PMWellCreate%flow_soln%itol_scaled_res = 1.0d-5
  PMWellCreate%flow_soln%itol_abs_update_p = 1.0d0
  PMWellCreate%flow_soln%itol_abs_update_s = 1.0d-5
  PMWellCreate%flow_soln%itol_rel_update_p = 1.0d-4
  PMWellCreate%flow_soln%itol_rel_update_s = 1.0d-4
  PMWellCreate%flow_soln%n_steps = 0
  PMWellCreate%flow_soln%n_newton = 0

  ! create the well transport solution object:
  allocate(PMWellCreate%tran_soln)
  nullify(PMWellCreate%tran_soln%prev_soln%aqueous_conc)
  nullify(PMWellCreate%tran_soln%prev_soln%aqueous_mass)
  nullify(PMWellCreate%tran_soln%residual)
  nullify(PMWellCreate%tran_soln%Jacobian)
  nullify(PMWellCreate%tran_soln%update)
  PMWellCreate%tran_soln%ndof = UNINITIALIZED_INTEGER
  PMWellCreate%tran_soln%not_converged = PETSC_TRUE
  PMWellCreate%tran_soln%converged = PETSC_FALSE
  PMWellCreate%tran_soln%cut_timestep = PETSC_FALSE
  PMWellCreate%tran_soln%max_iter = 8
  PMWellCreate%tran_soln%max_ts_cut = 20
  PMWellCreate%tran_soln%ts_cut_factor = 2.0d0
  PMWellCreate%tran_soln%ts_ramp_factor = 1.1d0
  PMWellCreate%tran_soln%itol_abs_res = 1.0d-8
  PMWellCreate%tran_soln%itol_scaled_res = 1.0d-4
  PMWellCreate%tran_soln%itol_abs_update = 1.0d0
  PMWellCreate%tran_soln%itol_rel_update = 1.0d-1
  PMWellCreate%tran_soln%n_steps = 0
  PMWellCreate%tran_soln%n_newton = 0

  ! strata list specific to well
  allocate(PMWellCreate%strata_list)
  nullify(PMWellCreate%strata_list%first)
  nullify(PMWellCreate%strata_list%last)
  nullify(PMWellCreate%strata_list%array)

end function PMWellCreate

! ************************************************************************** !

subroutine PMWellSetup(this)
  !
  ! Initializes variables associated with the well process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  !

  use Grid_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Input_Aux_module
  use Dataset_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Transport_Constraint_NWT_module
  use NW_Transport_Aux_module

  implicit none

  class(pm_well_type) :: this

  type(option_type), pointer :: option
  type(grid_type), pointer :: res_grid
  type(well_grid_type), pointer :: well_grid
  type(coupler_type), pointer :: source_sink
  type(input_type) :: input_dummy
  class(realization_subsurface_type), pointer :: realization
  class(dataset_ascii_type), pointer :: dataset_ascii
  type(point3d_type) :: dummy_h
  class(tran_constraint_coupler_nwt_type), pointer ::tran_constraint_coupler_nwt
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscInt, pointer :: h_global_id_unique(:)
  PetscInt, pointer :: h_rank_id(:)
  PetscReal :: diff_x,diff_y,diff_z
  PetscReal :: dh_x,dh_y,dh_z
  PetscReal :: total_length
  PetscReal :: top_of_reservoir, top_of_hole
  PetscReal :: bottom_of_reservoir, bottom_of_hole
  PetscReal :: max_diameter, xy_span
  PetscReal :: temp_real
  PetscReal :: dz_list(10000), res_dz_list(10000)
  PetscReal :: min_dz, dz, z
  PetscInt :: cell_id_list(10000)
  PetscInt :: local_id, ghosted_id, cur_id
  PetscInt :: local_id_well, local_id_res, res_cell_count
  PetscInt :: nsegments, nsegments_save
  PetscInt :: max_val, min_val
  PetscInt :: k, i
  PetscInt :: count1_local, count2_local
  PetscInt :: count1_global, count2_global
  PetscBool :: well_grid_res_is_OK = PETSC_FALSE
  PetscBool :: res_grid_cell_within_well_z
  PetscBool :: res_grid_cell_within_well_y
  PetscBool :: res_grid_cell_within_well_x
  PetscErrorCode :: ierr

  option => this%option
  realization => this%realization
  res_grid => realization%patch%grid
  well_grid => this%well_grid
  nsegments = this%well_grid%nsegments

  option%io_buffer = ' '
  call PrintMsg(option)
  option%io_buffer = 'WELLBORE_MODEL: Creating well grid discretization.... '
  call PrintMsg(option)

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

  if (nsegments == UNINITIALIZED_INTEGER) then

  ! Use reservoir grid info: 1 well cell per reservoir cell
    dz_list = UNINITIALIZED_DOUBLE
    res_dz_list = UNINITIALIZED_DOUBLE

    if (Initialized(this%min_dz)) then
      min_dz = this%min_dz
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
      if (res_cell_count <= this%well_res_ratio .and. cur_id == local_id) then
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
    allocate(well_grid%h(nsegments))
    allocate(well_grid%res_dz(nsegments))
    allocate(well_grid%h_local_id(nsegments))
    allocate(well_grid%h_ghosted_id(nsegments))
    allocate(well_grid%h_global_id(nsegments))
    allocate(well_grid%strata_id(nsegments))

    well_grid%res_dz(1:nsegments) = res_dz_list(1:nsegments)

    well_grid%strata_id(:) = UNINITIALIZED_INTEGER
    this%well_grid%nsegments = nsegments
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

  else
  ! Build an equally-spaced grid
    allocate(well_grid%dh(nsegments))
    allocate(well_grid%h(nsegments))
    allocate(well_grid%h_local_id(nsegments))
    allocate(well_grid%h_ghosted_id(nsegments))
    allocate(well_grid%h_global_id(nsegments))
    allocate(well_grid%strata_id(nsegments))

    well_grid%strata_id(:) = UNINITIALIZED_INTEGER

    well_grid%nconnections = well_grid%nsegments - 1

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

    ! Get the local_id for each well segment center from the reservoir grid
    well_grid%h_local_id(:) = -1

    do k = 1,well_grid%nsegments
      call GridGetLocalIDFromCoordinate(res_grid,well_grid%h(k), &
                                        option,local_id)
      well_grid%h_local_id(k) = local_id
      well_grid%h_ghosted_id(k) = res_grid%nL2G(local_id)
      well_grid%h_global_id(k) = res_grid%nG2A(well_grid%h_ghosted_id(k))
    enddo
    !write(*,*) 'well_grid%h_local_id', well_grid%h_local_id
    !write(*,*) 'well_grid%h_ghosted_id', well_grid%h_ghosted_id
    !write(*,*) 'well_grid%h_global_id', well_grid%h_global_id

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

  ! Check that no reservoir grid cells were skipped.
  ! Count how many of the h_global_id's are unique.
  ! This sum must be = to the number of reservoir cells that the
  !   well passes through.
  option%io_buffer = 'WELLBORE_MODEL: Checking well grid resolution.... '
  call PrintMsg(option)

  allocate(h_global_id_unique(nsegments))
  h_global_id_unique(:) = 0
  allocate(h_rank_id(nsegments))
  h_rank_id(:) = -999
  do k = 1,nsegments
    if (well_grid%h_local_id(k) > -999) then
      h_rank_id(k) = option%myrank
    endif
  enddo
  !write(*,*) 'h_rank_id(k)', h_rank_id

  min_val = minval(well_grid%h_global_id)-1
  max_val = maxval(well_grid%h_global_id)
  k = 0
  do while (min_val < max_val)
    k = k + 1
    min_val = minval(well_grid%h_global_id, &
                     mask=well_grid%h_global_id > min_val)
    h_global_id_unique(k) = min_val
  enddo

  !write(*,*) 'h_global_id_unique', h_global_id_unique

  count1_local = 0
  count1_global = 0
  do k = 1,nsegments
    if (h_global_id_unique(k) > 0) then
      count1_local = count1_local + 1
    endif
  enddo

  !write(*,*) '---> count1_local =', count1_local

  ! count1_local is the number of unique reservoir grid cells that the well has
  ! a connection to per MPI process. Next, all of the MPI processes need to
  ! sum up their counts and place the total in count1_global.
  call MPI_Allreduce(count1_local,count1_global,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  !write(*,*) '---> count1_global =', count1_global
  write(string,'(I0.5)') count1_global

  !write(*,*) '---> Passed 9'

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
      if (local_id_res == local_id_well) then
        count2_local = count2_local + 1
      endif
    endif
  enddo

  ! All of the MPI processes need to sum up their counts and place the
  ! total in count2_global.
  call MPI_Allreduce(count2_local,count2_global,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  write(string2,'(I0.5)') count2_global

  ! The only way we can ensure that the well discretization did not skip a
  ! reservoir cell, is if the number of unique global_id's that the well
  ! is connected to (count1) matches the number of reservoir grid cells that
  ! the well occupies (count2):
  if (count1_global == count2_global) well_grid_res_is_OK = PETSC_TRUE


  if (.not.well_grid_res_is_OK) then
    option%io_buffer = 'ERROR:  &
      &The number of reservoir grid cells that are occupied by the well &
      &(' // trim(string2) // ') is larger than the number of unique &
      &reservoir grid cells that have a connection to the well (' // &
      trim(string) // '). Therefore, some of the reservoir grid cells &
      &have been skipped and have no connection to the well. You must &
      &increase the resolution of the WELLBORE_MODEL grid. '
    call PrintMsg(option)
    option%io_buffer = '(see above)   Alternatively, &
      &if you are sure your well grid resolution is fine enough, try &
      &increasing the value set for the x-y search parameter &
      &WELLBORE_MODEL,WELL_GRID,XY_SEARCH_MULTIPLIER (default value = 10).'
    call PrintErrMsg(option)
  else
    option%io_buffer = 'WELLBORE_MODEL: &
      &Well grid resolution is adequate. No reservoir grid cell &
      &connections have been skipped.'
    call PrintMsg(option)
  endif

  allocate(this%well%ccid(nsegments))
  allocate(this%well%permeability(nsegments))
  allocate(this%well%phi(nsegments))

  allocate(this%well%area(nsegments))
  this%well%area = 3.14159*(this%well%diameter/2.d0)*(this%well%diameter/2.d0)

  allocate(this%well%volume(nsegments))
  this%well%volume = this%well%area*well_grid%dh

  allocate(this%well%liq%visc(nsegments))
  allocate(this%well%gas%visc(nsegments))

  allocate(this%well%mass_balance_liq(nsegments))

  this%flow_soln%ndof = this%nphase

  if (this%option%itranmode /= NULL_MODE) then
    this%transport = PETSC_TRUE
    if (this%option%itranmode /= NWT_MODE) then
      option%io_buffer ='The only transport mode allowed with the &
      &WELLBORE_MODEL is NWT_MODE.'
      call PrintErrMsg(option)
    endif
    this%nspecies = realization%reaction_nw%params%nspecies
    this%tran_soln%ndof = this%nspecies
  endif

  ! add a reservoir src/sink coupler for each well segment
  do k = 1,well_grid%nsegments
    write(string,'(I0.6)') k
    source_sink => CouplerCreate(SRC_SINK_COUPLER_TYPE)
    source_sink%name = 'well_segment_' // trim(string)

    ! ----- flow ------------------
    source_sink%flow_condition_name = 'well_segment_' // trim(string) // &
                                      '_flow_srcsink'
    source_sink%flow_condition => FlowConditionCreate(option)
    source_sink%flow_condition%name = source_sink%flow_condition_name
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
    source_sink%flow_condition%well%ctype = 'liq_average_density_kg/m3'

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

end subroutine PMWellSetup

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
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option
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
        cycle
    !-------------------------------------
      case('WIPP_INTRUSION_ZERO_VALUE')  ! [mol/m3-bulk]
        call InputReadDouble(input,option,this%bh_zero_value)
        call InputErrorMsg(input,option,'WIPP_INTRUSION_ZERO_VALUE', &
                           error_string)
        cycle
    !-------------------------------------
      case('MIN_DZ')
        call InputReadDouble(input,option,this%min_dz)
        call InputErrorMsg(input,option,'MIN_DZ', &
                           error_string)
        cycle
    !-------------------------------------
      case('MAX_WELL_RESERVOIR_CELL_RATIO')
        call InputReadInt(input,option,this%well_res_ratio)
        call InputErrorMsg(input,option,'MAX_WELL_RESERVOIR_CELL_RATIO', &
                           error_string)
        cycle
    !-------------------------------------
    end select

    ! Read sub-blocks within WELLBORE_MODEL block:
    error_string = 'WELLBORE_MODEL'
    call PMWellReadGrid(this,input,option,word,error_string,found)
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

subroutine PMWellReadGrid(pm_well,input,option,keyword,error_string,found)
  !
  ! Reads input file parameters associated with the well model grid.
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

  PetscInt :: num_errors
  character(len=MAXWORDLENGTH) :: word

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
          case('NUMBER_OF_SEGMENTS')
            call InputReadInt(input,option,pm_well%well_grid%nsegments)
            call InputErrorMsg(input,option,'NUMBER_OF_SEGMENTS',error_string)
        !-----------------------------
          case('TOP_OF_HOLE')
            call InputReadDouble(input,option,pm_well%well_grid%tophole(1))
            call InputErrorMsg(input,option,'TOP_OF_HOLE x-coordinate', &
                               error_string)
            call InputReadDouble(input,option,pm_well%well_grid%tophole(2))
            call InputErrorMsg(input,option,'TOP_OF_HOLE y-coordinate', &
                               error_string)
            call InputReadDouble(input,option,pm_well%well_grid%tophole(3))
            call InputErrorMsg(input,option,'TOP_OF_HOLE z-coordinate', &
                               error_string)
        !-----------------------------
          case('BOTTOM_OF_HOLE')
            call InputReadDouble(input,option,pm_well%well_grid%bottomhole(1))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE x-coordinate', &
                               error_string)
            call InputReadDouble(input,option,pm_well%well_grid%bottomhole(2))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE y-coordinate', &
                               error_string)
            call InputReadDouble(input,option,pm_well%well_grid%bottomhole(3))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE z-coordinate', &
                               error_string)
        !-----------------------------
          case('XY_SEARCH_MULTIPLIER')
            call InputReadInt(input,option,pm_well%well_grid%xy_span_multiplier)
            call InputErrorMsg(input,option,'XY_SEARCH_MULTIPLIER',error_string)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      ! ----------------- error messaging -------------------------------------
      !if (Uninitialized(pm_well%well_grid%nsegments)) then
      !  option%io_buffer = 'ERROR: NUMBER_OF_SEGMENTS must be specified &
      !                     &in the ' // trim(error_string) // ' block.'
      !  call PrintMsg(option)
      !  num_errors = num_errors + 1
      !endif
      !if (pm_well%well_grid%nsegments < 3) then
      !  option%io_buffer = 'ERROR: The well must consist of >= 3 segments &
      !                     &in the ' // trim(error_string) // ' block.'
      !  call PrintMsg(option)
      !  num_errors = num_errors + 1
      !endif
      if (Uninitialized(pm_well%well_grid%tophole(1)) .or. &
          Uninitialized(pm_well%well_grid%tophole(2)) .or. &
          Uninitialized(pm_well%well_grid%tophole(3))) then
        option%io_buffer = 'ERROR: TOP_OF_HOLE must be specified &
                           &in the ' // trim(error_string) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Uninitialized(pm_well%well_grid%bottomhole(1)) .or. &
          Uninitialized(pm_well%well_grid%bottomhole(2)) .or. &
          Uninitialized(pm_well%well_grid%bottomhole(3))) then
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
  character(len=MAXWORDLENGTH) :: word
  PetscReal, pointer :: temp_diameter(:)
  PetscReal, pointer :: temp_friction(:)
  PetscReal, pointer :: temp_well_index(:)
  PetscReal, pointer :: temp_well_perm(:)
  PetscReal, pointer :: temp_well_phi(:)

  error_string = trim(error_string) // ',WELL'
  found = PETSC_TRUE
  num_errors = 0

  read_max = 200
  allocate(temp_diameter(read_max))
  allocate(temp_friction(read_max))
  allocate(temp_well_index(read_max))
  allocate(temp_well_perm(read_max))
  allocate(temp_well_phi(read_max))
  temp_diameter(:) = UNINITIALIZED_DOUBLE
  temp_friction(:) = UNINITIALIZED_DOUBLE
  temp_well_index(:) = UNINITIALIZED_DOUBLE
  temp_well_perm(:) = UNINITIALIZED_DOUBLE
  temp_well_phi(:) = UNINITIALIZED_DOUBLE

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
            allocate(pm_well%well%f(num_read))
            pm_well%well%f(1:num_read) = temp_friction(1:num_read)
        !-----------------------------
          case('WELL_INDEX')
            do k = 1,read_max
              call InputReadDouble(input,option,temp_well_index(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for WELL_INDEX &
                &must be provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(pm_well%well%WI_base(num_read))
            pm_well%well%WI_base(1:num_read) = temp_well_index(1:num_read)
        !-----------------------------
          case('WELL_INDEX_MODEL')
            call InputReadWord(input,option,word,PETSC_TRUE)
            select case(word)
              case('PEACEMAN_ISO')
                pm_well%well%WI_model = PEACEMAN_ISO 
              case('PEACEMAN_ANISOTROPIC')
                pm_well%well%WI_model = PEACEMAN_ANISOTROPIC
              case default
                option%io_buffer = 'Unrecognized option for WELL_INDEX_MODEL &
                &in the ' // trim(error_string) // ' block. Default is isotropic &
                &Peaceman (PEACEMAN_ISO).'
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
      if (.not.associated(pm_well%well%f)) then
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

  end subroutine PMWellReadWell

! ************************************************************************** !

subroutine PMWellReadWellBCs(this,input,option,keyword,error_string,found)
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

  class(pm_well_type) :: this
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
                      call InputReadDouble(input,option,this%well%bh_p)
                      call InputReadAndConvertUnits(input,this%well%bh_p, &
                           'Pa','WELL_BOUNDARY_CONDITIONS,BOTTOM_OF_HOLE,&
                           &LIQUID_PRESSURE',option)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('GAS_SATURATION')
                      call InputReadDouble(input,option,this%well%bh_sg)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('PRESSURE_SET_BY_RESERVOIR')
                      this%well%bh_p_set_by_reservoir = PETSC_TRUE
                  !-----------------------------
                    case('SATURATION_SET_BY_RESERVOIR')
                      this%well%bh_sg_set_by_reservoir = PETSC_TRUE
                  !-----------------------------
                    case('LIQUID_RATE')
                      call InputReadDouble(input,option,this%well%bh_ql)
                      call InputReadAndConvertUnits(input,this%well%bh_ql, &
                           'm/s','WELL_BOUNDARY_CONDITIONS,BOTTOM_OF_HOLE,&
                           &LIQUID_RATE',option)
                  !-----------------------------
                    case('GAS_RATE')
                      call InputReadDouble(input,option,this%well%bh_qg)
                      call InputReadAndConvertUnits(input,this%well%bh_qg, &
                           'm/s','WELL_BOUNDARY_CONDITIONS,BOTTOM_OF_HOLE,&
                           &GAS_RATE',option)
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
                      call InputReadDouble(input,option,this%well%th_p)
                      call InputReadAndConvertUnits(input,this%well%th_p, &
                           'Pa','WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE,&
                           &LIQUID_PRESSURE',option)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('GAS_SATURATION')
                      call InputReadDouble(input,option,this%well%th_sg)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('LIQUID_RATE')
                      call InputReadDouble(input,option,this%well%th_ql)
                      call InputReadAndConvertUnits(input,this%well%th_ql, &
                           'm/s','WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE,&
                           &LIQUID_RATE',option)
                  !-----------------------------
                    case('GAS_RATE')
                      call InputReadDouble(input,option,this%well%th_qg)
                      call InputReadAndConvertUnits(input,this%well%th_qg, &
                           'm/s','WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE,&
                           &GAS_RATE',option)
                  !-----------------------------
                    case('TRANSPORT_CONDITION')
                      call InputReadWord(input,option, &
                                   this%well%tran_condition_name,PETSC_TRUE)
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
      if (Initialized(this%well%bh_p) .and. &
          this%well%bh_p_set_by_reservoir) then
        option%io_buffer = 'Either keyword BOTTOM_OF_HOLE,PRESSURE or keyword &
          &BOTTOM_OF_HOLE,PRESSURE_SET_BY_RESERVOIR must be provided in the ' &
          // trim(error_string) // ' block, but NOT BOTH.'
        call PrintErrMsg(option)
      endif

      if (.not. (Initialized(this%well%bh_ql) .and. &
                 Initialized(this%well%bh_qg))) then
        if ( .not. ((Initialized(this%well%bh_p) .or. &
                     this%well%bh_p_set_by_reservoir) .and. &
                    (Initialized(this%well%bh_sg) .or. &
                     this%well%bh_sg_set_by_reservoir))) then

          option%io_buffer ='WIPP_DARCY well model needs both Dirichlet &
             &LIQUID_PRESSURE and GAS_SATURATION set in the ' &
             // trim(error_string) // ' block.'
          call PrintErrMsg(option)
        endif
      endif

      if ((Initialized(this%well%th_p).or.Initialized(this%well%th_sg)) .and. &
          .not.(Initialized(this%well%th_p).and.Initialized(this%well%th_sg))) &
          then
        option%io_buffer ='WIPP_DARCY well model needs both Dirichlet &
           &LIQUID_PRESSURE and GAS_SATURATION set in the ' &
           // trim(error_string) // ' block.'
        call PrintErrMsg(option)
      endif

      if ((this%option%itranmode /= NULL_MODE) .and. &
          (trim(this%well%tran_condition_name) == '')) then
        option%io_buffer ='A TRANSPORT_CONDITION was not provided in the &
          &WELL_BOUNDARY_CONDITIONS,TOP_OF_HOLE block.'
        call PrintErrMsg(option)
      endif

  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (Initialized(this%well%bh_p)) this%flow_soln%bh_p = PETSC_TRUE
  if (this%well%bh_p_set_by_reservoir) this%flow_soln%bh_p = PETSC_TRUE
  if (Initialized(this%well%th_p)) this%flow_soln%th_p = PETSC_TRUE

  if (Initialized(this%well%bh_ql)) this%flow_soln%bh_q = PETSC_TRUE
  if (Initialized(this%well%th_ql)) this%flow_soln%th_q = PETSC_TRUE

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

subroutine PMWellReadWellModelType(this,input,option,keyword,error_string, &
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

  class(pm_well_type) :: this
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
          case('CONSTANT_PRESSURE')
            this%well%well_model_type = 'CONSTANT_PRESSURE'
        !-----------------------------
          case('CONSTANT_PRESSURE_HYDROSTATIC')
            this%well%well_model_type = 'CONSTANT_PRESSURE_HYDROSTATIC'
        !-----------------------------
          case('CONSTANT_RATE')
            this%well%well_model_type = 'CONSTANT_RATE'
        !-----------------------------
          case('WIPP_DARCY')
            this%well%well_model_type = 'WIPP_DARCY'
            !call PMWellReadDarcyInput(this,input,option,keyword,&
            !                           error_string,found)
        !-----------------------------
          case('STEADY_STATE')
            this%well%well_model_type = 'STEADY_STATE'
        !-----------------------------
          case('FULL_MOMENTUM')
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

subroutine PMWellReadWellOutput(this,input,option,keyword,error_string, &
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

  class(pm_well_type) :: this
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
            this%well%output_pl = PETSC_TRUE
        !-----------------------------
          case('WELL_GAS_PRESSURE')
            this%well%output_pg = PETSC_TRUE
        !-----------------------------
          case('WELL_AQ_CONC')
            this%well%output_aqc = PETSC_TRUE
        !-----------------------------
          case('WELL_AQ_MASS')
            this%well%output_aqm = PETSC_TRUE
        !-----------------------------
          case('WELL_LIQ_Q')
            this%well%liq%output_Q = PETSC_TRUE
        !-----------------------------
          case('WELL_GAS_Q')
            this%well%gas%output_Q = PETSC_TRUE
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
      case('WELL_GRID','WELL','WELL_MODEL_TYPE','WELL_FLOW_SOLVER', &
           'WELL_TRANSPORT_SOLVER','WELL_MODEL_OUTPUT')
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

  use Condition_module
  use Strata_module
  use Output_Aux_module
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
  PetscInt :: nsegments, k, i, j
  PetscReal :: curr_time

  curr_time = this%option%time
  nsegments = this%well_grid%nsegments

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

  ! loop thru well stratas and mark them as active or inactive
  well_strata => this%strata_list%first
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

  do k = 1,nsegments
    well_strata => this%strata_list%first
    do
      if (.not.associated(well_strata)) exit
      if ((any(well_strata%region%cell_ids == &
               this%well_grid%h_ghosted_id(k))) .and. &
          (well_strata%active)) then
        this%well_grid%strata_id(k) = well_strata%id
      endif
      well_strata => well_strata%next
    enddo
    if (Uninitialized(this%well_grid%strata_id(k))) then
      this%option%io_buffer =  'At least one WELLBORE_MODEL grid segment has &
        &not been assigned with a REGION and MATERIAL_PROPERTY with the use &
        &of the STRATA block.'
      call PrintErrMsg(this%option)
    endif
  enddo

  allocate(this%flow_soln%residual(nsegments*this%flow_soln%ndof))
  allocate(this%flow_soln%update(nsegments*this%flow_soln%ndof))
  this%flow_soln%residual(:) = UNINITIALIZED_DOUBLE
  this%flow_soln%update(:) = UNINITIALIZED_DOUBLE

  allocate(this%flow_soln%Jacobian(this%nphase*nsegments,this%nphase*nsegments))
  do i = 1,this%nphase*nsegments
    do j = 1,this%nphase*nsegments
      this%flow_soln%Jacobian(i,j) = UNINITIALIZED_DOUBLE
    enddo
  enddo

  allocate(this%flow_soln%prev_soln%pl(nsegments))
  allocate(this%flow_soln%prev_soln%sg(nsegments))

  if (this%transport) then
    allocate(this%tran_soln%residual(nsegments*this%tran_soln%ndof))
    allocate(this%tran_soln%update(nsegments*this%tran_soln%ndof))
    this%tran_soln%residual(:) = UNINITIALIZED_DOUBLE
    this%tran_soln%update(:) = UNINITIALIZED_DOUBLE

    allocate(this%tran_soln%Jacobian(this%nspecies*nsegments, &
                                     this%nspecies*nsegments))
    this%tran_soln%Jacobian = UNINITIALIZED_DOUBLE

    allocate(this%tran_soln%prev_soln%aqueous_conc(this%nspecies, &
                                                   this%well_grid%nsegments))
    allocate(this%tran_soln%prev_soln%aqueous_mass(this%nspecies, &
                                                   this%well_grid%nsegments))
  endif

  allocate(this%well%WI(nsegments))
  allocate(this%well%r0(nsegments))
  allocate(this%well%pl(nsegments))
  allocate(this%well%pg(nsegments))
  allocate(this%well%ql(nsegments-1))
  allocate(this%well%qg(nsegments-1))
  allocate(this%well%ql_bc(2))
  allocate(this%well%qg_bc(2))
  allocate(this%well%ql_kmol(nsegments-1))
  allocate(this%well%qg_kmol(nsegments-1))
  allocate(this%well%ql_kmol_bc(2))
  allocate(this%well%qg_kmol_bc(2))
  allocate(this%well%liq_cum_mass(nsegments))
  allocate(this%well%liq_mass(nsegments))
  this%well%liq_cum_mass = 0.d0 
  this%well%liq_mass = 0.d0
  
  if (this%transport) then
    allocate(this%well%aqueous_conc(this%nspecies,nsegments))
    allocate(this%well%aqueous_mass(this%nspecies,nsegments))
    allocate(this%well%aqueous_conc_th(this%nspecies))
  endif

  allocate(this%pert(nsegments,this%nphase))
  this%pert = 0.d0

  allocate(this%well%liq%s(nsegments))
  this%well%liq%rho0 = this%option%flow%reference_density(1)
  allocate(this%well%liq%rho(nsegments))
  this%well%liq%rho(:) = this%well%liq%rho0
  allocate(this%well%liq%visc(nsegments))
  allocate(this%well%liq%Q(nsegments))
  allocate(this%well%liq%kr(nsegments))
  allocate(this%well%gas%s(nsegments))
  this%well%gas%rho0 = this%option%flow%reference_density(2)
  allocate(this%well%gas%rho(nsegments))
  this%well%gas%rho(:) = this%well%gas%rho0
  allocate(this%well%gas%visc(nsegments))
  allocate(this%well%gas%Q(nsegments))
  allocate(this%well%gas%kr(nsegments))

  ! Initialize perturbations
  allocate(this%pert(nsegments,this%nphase))
  this%pert = 0.d0

  ! Initialize stored aux variables at well perturbations
  do k = 1,this%nphase
    this%well_pert(k)%well_model_type = this%well%well_model_type
    this%well_pert(k)%wi_model = this%well%wi_model

    allocate(this%well_pert(k)%WI(nsegments))
    allocate(this%well_pert(k)%r0(nsegments))
    allocate(this%well_pert(k)%pl(nsegments))
    allocate(this%well_pert(k)%pg(nsegments))
    allocate(this%well_pert(k)%diameter(nsegments))
    allocate(this%well_pert(k)%WI_base(nsegments))
    allocate(this%well_pert(k)%ccid(nsegments))
    allocate(this%well_pert(k)%permeability(nsegments))
    allocate(this%well_pert(k)%phi(nsegments))
    allocate(this%well_pert(k)%f(nsegments))
    allocate(this%well_pert(k)%area(nsegments))
    allocate(this%well_pert(k)%volume(nsegments))
    allocate(this%well_pert(k)%liq%visc(nsegments))
    allocate(this%well_pert(k)%gas%visc(nsegments))
    allocate(this%well_pert(k)%liq%s(nsegments))
    allocate(this%well_pert(k)%gas%s(nsegments))
    allocate(this%well_pert(k)%liq%rho(nsegments))
    allocate(this%well_pert(k)%gas%rho(nsegments))
    allocate(this%well_pert(k)%liq%Q(nsegments))
    allocate(this%well_pert(k)%gas%Q(nsegments))
    allocate(this%well_pert(k)%liq%kr(nsegments))
    allocate(this%well_pert(k)%gas%kr(nsegments))

    call PMWellCopyWell(this%well,this%well_pert(k))

    this%well_pert(k)%liq%rho0 = this%option%flow%reference_density(1)
    this%well_pert(k)%gas%rho0 = this%option%flow%reference_density(2)
  enddo

  allocate(this%reservoir%p_l(nsegments))
  allocate(this%reservoir%p_g(nsegments))
  allocate(this%reservoir%s_l(nsegments))
  allocate(this%reservoir%s_g(nsegments))
  allocate(this%reservoir%mobility_l(nsegments))
  allocate(this%reservoir%mobility_g(nsegments))
  allocate(this%reservoir%kr_l(nsegments))
  allocate(this%reservoir%kr_g(nsegments))
  allocate(this%reservoir%rho_l(nsegments))
  allocate(this%reservoir%rho_g(nsegments))
  allocate(this%reservoir%visc_l(nsegments))
  allocate(this%reservoir%visc_g(nsegments))
  allocate(this%reservoir%e_por(nsegments))
  allocate(this%reservoir%kx(nsegments))
  allocate(this%reservoir%ky(nsegments))
  allocate(this%reservoir%kz(nsegments))
  allocate(this%reservoir%dx(nsegments))
  allocate(this%reservoir%dy(nsegments))
  allocate(this%reservoir%dz(nsegments))
  allocate(this%reservoir%volume(nsegments))
  if (this%transport) then
    allocate(this%reservoir%aqueous_conc(this%nspecies, &
                                         this%well_grid%nsegments))
    allocate(this%reservoir%aqueous_mass(this%nspecies, &
                                         this%well_grid%nsegments))
  endif

  if (this%transport) then
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
  endif

  if (this%transport) then
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

  use Strata_module

  implicit none

  class(pm_well_type) :: this

  type(strata_type), pointer :: strata
  PetscInt :: k
  PetscReal :: curr_time

  curr_time = this%option%time - this%option%flow_dt

  if (Initialized(this%intrusion_time_start) .and. &
      (curr_time < this%intrusion_time_start)) return

  ! update the reservoir object with current reservoir properties
  call PMWellUpdateReservoir(this)
  call PMWellComputeWellIndex(this)

  ! loop thru strata and mark them as active or inactive
  strata => this%strata_list%first
  do
    if (.not.associated(strata)) exit
    if (Initialized(strata%start_time) .and. &
        Initialized(strata%final_time)) then
      if ((curr_time >= strata%start_time) .and. &
          (curr_time <= strata%final_time))  then
        strata%active = PETSC_TRUE
      else
        strata%active = PETSC_FALSE
      endif
    else
      strata%active = PETSC_TRUE
    endif
    strata => strata%next
  enddo

  ! update the well_grid%strata_id assignment for active strata
  this%well_grid%strata_id(:) = 0
  do k = 1,this%well_grid%nsegments
    strata => this%strata_list%first
    do
      if (.not.associated(strata)) exit
      if ((any(strata%region%cell_ids == &
               this%well_grid%h_ghosted_id(k))) .and. &
          (strata%active)) then
        this%well_grid%strata_id(k) = strata%id
      endif
      strata => strata%next
    enddo
  enddo
  if (any(this%well_grid%strata_id == 0)) then
    this%option%io_buffer =  'At least one WELLBORE_MODEL grid segment has not &
        &been assigned with a REGION and MATERIAL_PROPERTY with the use of the &
        &STRATA block. Check the STRATA START_TIME/FINAL_TIME cards associated &
        &with the WELL keyword.'
    call PrintErrMsg(this%option)
  endif

  if (initialize_well) then
    ! enter here if its the very first timestep
    call PMWellInitializeWell(this)
  endif

  if (this%well%bh_p_set_by_reservoir) then
    this%well%bh_p = this%reservoir%p_l(1)
  endif
  if (this%well%bh_sg_set_by_reservoir) then
    this%well%bh_sg = this%reservoir%s_g(1)
  endif

  call PMWellUpdatePropertiesFlow(this,this%well,&
                        this%realization%patch%characteristic_curves_array,&
                        this%realization%option)
  this%dt_flow = this%realization%option%flow_dt

  if (this%transport) then
    call PMWellUpdatePropertiesTran(this)
    this%dt_tran = this%dt_flow
  endif

end subroutine PMWellInitializeTimestep

! ************************************************************************** !

subroutine PMWellInitializeWell(this)
  !
  ! Initializes the well for the first time step.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  implicit none

  class(pm_well_type) :: this
  type(strata_type), pointer :: strata

  PetscInt :: k

  ! set initial flow parameters to the reservoir flow parameters
  this%well%pl = this%reservoir%p_l
  this%well%pg = this%reservoir%p_g
  this%well%liq%s = this%reservoir%s_l
  this%well%gas%s = this%reservoir%s_g
  this%well%liq%rho = this%reservoir%rho_l
  this%well%gas%rho = this%reservoir%rho_g
  this%well%liq%visc = this%reservoir%visc_l
  this%well%gas%visc = this%reservoir%visc_g
  ! update the Darcy fluxes within the well

  this%well_pert(ONE_INTEGER)%pl = this%reservoir%p_l
  this%well_pert(ONE_INTEGER)%pg = this%reservoir%p_g
  this%well_pert(ONE_INTEGER)%liq%s = this%reservoir%s_l
  this%well_pert(ONE_INTEGER)%gas%s = this%reservoir%s_g
  this%well_pert(ONE_INTEGER)%liq%rho = this%reservoir%rho_l
  this%well_pert(ONE_INTEGER)%gas%rho = this%reservoir%rho_g
  this%well_pert(ONE_INTEGER)%liq%visc = this%reservoir%visc_l
  this%well_pert(ONE_INTEGER)%gas%visc = this%reservoir%visc_g

  this%well_pert(TWO_INTEGER)%pl = this%reservoir%p_l
  this%well_pert(TWO_INTEGER)%pg = this%reservoir%p_g
  this%well_pert(TWO_INTEGER)%liq%s = this%reservoir%s_l
  this%well_pert(TWO_INTEGER)%gas%s = this%reservoir%s_g
  this%well_pert(TWO_INTEGER)%liq%rho = this%reservoir%rho_l
  this%well_pert(TWO_INTEGER)%gas%rho = this%reservoir%rho_g
  this%well_pert(TWO_INTEGER)%liq%visc = this%reservoir%visc_l
  this%well_pert(TWO_INTEGER)%gas%visc = this%reservoir%visc_g

  this%flow_soln%prev_soln%pl = this%well%pl
  this%flow_soln%prev_soln%sg = this%well%gas%s

  ! Link well material properties
  do k = 1,this%well_grid%nsegments
    strata => this%strata_list%first
    do
      if (.not.associated(strata)) exit
      if (strata%id == this%well_grid%strata_id(k)) then
        this%well%ccid(k) = strata%material_property%saturation_function_id
        this%well%permeability(k) = strata%material_property%permeability(3,3)
        this%well%phi(k) = strata%material_property%porosity
        exit
      endif
      strata => strata%next
    enddo
  enddo

  call PMWellComputeWellIndex(this)

  do k = 1,this%nphase
    this%well_pert(k)%permeability(:) = this%well%permeability(:)
    this%well_pert(k)%phi(:) = this%well%phi(:)
    this%well_pert(k)%ccid(:) = this%well%ccid(:)
    this%well_pert(k)%WI(:) = this%well%WI(:)
    this%well_pert(k)%r0(:) = this%well%r0(:)
  enddo

  ! set initial transport parameters to the reservoir transport parameters
  if (this%transport) then
    if (Initialized(this%intrusion_time_start)) then
      ! set the borehole concentrations to the borehole zero value now
      do k = 1,this%well_grid%nsegments
        this%well%aqueous_mass(:,k) = this%bh_zero_value*this%well%volume(k)
        this%well%aqueous_conc(:,k) = &
          this%well%aqueous_mass(:,k) / &                           ! [mol]
          (this%well%phi(k)*this%well%volume(k)*this%well%liq%s(k)) ! [m3-liq]
      enddo
    else
      ! set the wellbore concentrations to the reservoir values
      this%well%aqueous_conc = this%reservoir%aqueous_conc
      do k = 1,this%well_grid%nsegments
        this%well%aqueous_mass(:,k) = this%well%aqueous_conc(:,k) * &
                this%well%phi(k) * this%well%volume(k) * this%well%liq%s(k)
      enddo
    endif
    this%tran_soln%prev_soln%aqueous_conc = this%well%aqueous_conc
    this%tran_soln%prev_soln%aqueous_mass = this%well%aqueous_mass
  endif

  initialize_well = PETSC_FALSE

end subroutine PMWellInitializeWell

! ************************************************************************** !

subroutine PMWellUpdateReservoir(this)
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

  class(pm_well_type) :: this

  type(wippflo_auxvar_type), pointer :: wippflo_auxvar
  type(nw_transport_auxvar_type), pointer :: nwt_auxvar
  type(material_auxvar_type), pointer :: material_auxvar
  type(grid_type), pointer :: res_grid
  type(option_type), pointer :: option
  PetscInt :: k
  PetscInt :: ghosted_id

  option => this%option

  res_grid => this%realization%patch%grid

  do k=1,size(this%well_grid%h_ghosted_id)
    ghosted_id = this%well_grid%h_ghosted_id(k)

    wippflo_auxvar => &
      this%realization%patch%aux%wippflo%auxvars(0,ghosted_id)
    if (this%transport) then
      nwt_auxvar => &
        this%realization%patch%aux%nwt%auxvars(ghosted_id)
    endif
    material_auxvar => &
      this%realization%patch%aux%material%auxvars(ghosted_id)

    this%reservoir%p_l(k) = wippflo_auxvar%pres(option%liquid_phase)
    this%reservoir%p_g(k) = wippflo_auxvar%pres(option%gas_phase)
    this%reservoir%s_l(k) = wippflo_auxvar%sat(option%liquid_phase)
    this%reservoir%s_g(k) = wippflo_auxvar%sat(option%gas_phase)
    this%reservoir%mobility_l(k) = &
      wippflo_auxvar%mobility(option%liquid_phase)
    this%reservoir%mobility_g(k) = wippflo_auxvar%mobility(option%gas_phase)
    this%reservoir%kr_l(k) = wippflo_auxvar%kr(option%liquid_phase)
    this%reservoir%kr_g(k) = wippflo_auxvar%kr(option%gas_phase)
    this%reservoir%rho_l(k) = wippflo_auxvar%den_kg(option%liquid_phase)
    this%reservoir%rho_g(k) = wippflo_auxvar%den_kg(option%gas_phase)
    this%reservoir%visc_l(k) = wippflo_auxvar%mu(option%liquid_phase)
    this%reservoir%visc_g(k) = wippflo_auxvar%mu(option%gas_phase)
    this%reservoir%e_por(k) = wippflo_auxvar%effective_porosity

    this%reservoir%kx(k) = material_auxvar%permeability(1)
    this%reservoir%ky(k) = material_auxvar%permeability(2)
    this%reservoir%kz(k) = material_auxvar%permeability(3)
    this%reservoir%volume(k) = material_auxvar%volume

    if (res_grid%itype == STRUCTURED_GRID) then
      this%reservoir%dx(k) = res_grid%structured_grid%dx(ghosted_id)
      this%reservoir%dy(k) = res_grid%structured_grid%dy(ghosted_id)
      this%reservoir%dz(k) = res_grid%structured_grid%dz(ghosted_id)
    else
      this%reservoir%dz(k) = this%well_grid%res_dz(k)
      this%reservoir%dx(k) = sqrt(material_auxvar%volume/ &
                                  this%reservoir%dz(k))
      this%reservoir%dy(k) = this%reservoir%dx(k)
    endif

    if (this%transport) then
      this%reservoir%aqueous_conc(:,k) = nwt_auxvar%aqueous_eq_conc(:)
      this%reservoir%aqueous_mass(:,k) = &
            this%reservoir%aqueous_conc(:,k) * this%reservoir%e_por(k) * &
            this%reservoir%volume(k) * this%reservoir%s_l(k)
    endif

  enddo

end subroutine PMWellUpdateReservoir

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

  implicit none

  class(pm_well_type) :: this

  PetscReal :: curr_time

  curr_time = this%option%time - this%option%flow_dt

  if (Initialized(this%intrusion_time_start) .and. &
      (curr_time < this%intrusion_time_start)) return

  call PMWellUpdateReservoirSrcSink(this)

  call PMWellUpdateMass(this)

  call PMWellMassBalance(this)

  if (this%print_well) then
    call PMWellOutput(this)
  endif

end subroutine PMWellFinalizeTimestep

! ************************************************************************** !

subroutine PMWellUpdateReservoirSrcSink(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/21/2021

  use Coupler_module
  use NW_Transport_Aux_module
  use Transport_Constraint_NWT_module

  implicit none

  class(pm_well_type) :: this

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: srcsink_name
  type(nw_transport_auxvar_type), pointer :: nwt_auxvar
  type(coupler_type), pointer :: source_sink
  PetscInt :: k, ghosted_id
  PetscReal :: well_delta_liq, well_delta_gas
  PetscReal :: density_avg

  do k = 1,this%well_grid%nsegments
    write(string,'(I0.6)') k
    srcsink_name = 'well_segment_' // trim(string)

    ghosted_id = this%well_grid%h_ghosted_id(k)

    ! [kg-liq/m3]
    density_avg = 0.5d0 * (this%well%liq%rho(k) + this%reservoir%rho_l(k))

    source_sink => this%realization%patch%source_sink_list%first
    do
      if (.not.associated(source_sink)) exit

      if (trim(srcsink_name) == trim(source_sink%name)) then
        source_sink%flow_condition%general%rate%dataset%rarray(1) = &
          -1.d0 * this%well%liq%Q(k) * FMWH2O ! [kg/s]
        source_sink%flow_condition%general%rate%dataset%rarray(2) = &
          -1.d0 * this%well%gas%Q(k) * fmw_comp(TWO_INTEGER) ! [kg/s]

        source_sink%flow_condition%general%liquid_pressure%aux_real(1) = &
                                                           this%well%pl(k)
        source_sink%flow_condition%general%gas_pressure%aux_real(1) = &
                                                           this%well%pg(k)
        well_delta_liq = this%well%pl(k) - this%reservoir%p_l(k)
        well_delta_gas = this%well%pg(k) - this%reservoir%p_g(k)
        source_sink%flow_condition%general%liquid_pressure%aux_real(2) = &
                                                           well_delta_liq
        source_sink%flow_condition%general%gas_pressure%aux_real(2) = &
                                                           well_delta_gas
        source_sink%flow_condition%well%aux_real(1) = density_avg ! kg/m3

        this%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
             well%pl = this%well%pl(k)
        this%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
             well%pg = this%well%pg(k)
        this%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
             well%dpl = well_delta_liq
        this%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
             well%dpg = well_delta_gas
        this%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
             well%Ql = this%well%liq%Q(k)
        this%realization%patch%aux%wippflo%auxvars(:,ghosted_id)%&
             well%Qg = this%well%gas%Q(k)

        if (this%transport) then
          ! access nwt_auxvar from the tran_condition
          nwt_auxvar => &
            TranConstraintNWTGetAuxVar(source_sink%tran_condition% &
                                       cur_constraint_coupler)
          ! modify nwt_auxvar from the tran_condition
          nwt_auxvar%aqueous_eq_conc(:) = this%well%aqueous_conc(:,k)

          this%realization%patch%aux%nwt%auxvars(ghosted_id)%&
            well%AQ_conc(1:this%nspecies) = this%well%aqueous_conc(:,k)
          this%realization%patch%aux%nwt%auxvars(ghosted_id)%&
            well%AQ_mass(1:this%nspecies) = this%well%aqueous_mass(:,k)
        endif
        exit
      endif

      source_sink => source_sink%next
    enddo

  enddo

end subroutine PMWellUpdateReservoirSrcSink

! ************************************************************************** !

subroutine PMWellResidualFlow(this)
  !
  ! Author: Michael Nole
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this

  PetscInt :: i, iup, idn
  PetscReal :: res_accum(this%nphase)
  PetscReal :: res_src_sink(this%nphase)
  PetscReal :: res_flux(this%nphase)
  PetscReal :: res_flux_bc(2*this%nphase)

  res_accum = 0.d0
  res_src_sink = 0.d0
  res_flux = 0.d0
  res_flux_bc = 0.d0

  select case(this%well%well_model_type)
    !-------------------------------------------------------------------------
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    !-------------------------------------------------------------------------
    case('WIPP_DARCY')

      call PMWellBCFlux(this,this%well,res_flux_bc,PETSC_TRUE)

      do i = 1,this%well_grid%nsegments
        iup = i
        idn = i + 1

        ! Accumulation Term
        call PMWellAccumulationFlow(this,this%well,i,res_accum)

        this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) + &
               res_accum(ONE_INTEGER)

        ! Source/Sink Term
        call PMWellSrcSink(this,this%well,i,res_src_sink)

        this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) = &
             this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) + &
             res_src_sink(ONE_INTEGER)
  
        ! Flux Term
        if (i < this%well_grid%nsegments) then
          call PMWellFlux(this,this%well,this%well,iup,idn,res_flux,PETSC_TRUE)
        endif

        if (i == 1) then
          ! Water mass residual in cell i+1: Subtract flux to i+1 cell
          this%flow_soln%residual(this%flow_soln%ndof*i+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*i+1) &
               - res_flux(1)

          ! Water mass residual in cell i: Add flux in from BC,
          ! add flux to i+1 cell
          this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) &
               - res_flux_bc(1) + res_flux(1)

        elseif (i < this%well_grid%nsegments) then
          ! Water mass residual in cell i: Subtract flux to i+1 cell
          this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) &
               + res_flux(1)
          ! Water mass residual in cell i+1: Add flux to i+1 cell
          this%flow_soln%residual(this%flow_soln%ndof*i+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*i+1) &
               - res_flux(1)
        else
          ! Water mass residual in cell i: Subtract flux to BC
          this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) &
               - res_flux_bc(3)
        endif

        if (this%nphase == 2) then

          this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) = &
               this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) + &
               res_accum(TWO_INTEGER)

          this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) = &
               this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) + &
               res_src_sink(TWO_INTEGER)

          if (i == 1) then
            ! Air mass residual in cell i+1: Subtract flux to i+1 cell
            this%flow_soln%residual(this%flow_soln%ndof*i+2) = &
                 this%flow_soln%residual(this%flow_soln%ndof*i+2) &
                 - res_flux(2)
            ! Air mass residual in cell i: Subtract flux in from BC,
            ! add flux to i+1 cell
            this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) = &
                 this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) &
                 - res_flux_bc(2) + res_flux(2)
          elseif (i < this%well_grid%nsegments) then
            ! Air mass residual in cell i: Subtract flux to i+1 cell
            this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) = &
                 this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) &
                 + res_flux(2)
            ! Air mass residual in cell i+1: Add flux to i+1 cell
            this%flow_soln%residual(this%flow_soln%ndof*i+2) = &
                 this%flow_soln%residual(this%flow_soln%ndof*i+2) &
                 - res_flux(2)
          else
            ! Air mass residual in cell i: Subtract flux to BC
            this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) = &
                 this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) &
                 - res_flux_bc(4)
          endif
        endif
      enddo
    !-------------------------------------------------------------------------
    case('FULL_MOMENTUM')
    !-------------------------------------------------------------------------
  end select

end subroutine PMWellResidualFlow

! ************************************************************************** !

subroutine PMWellResidualTran(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  implicit none

  class(pm_well_type) :: this

  ! at this time, the tran_soln%residual vector has been zero'd out and has
  ! been loaded with the fixed accumulation divided by current dt

  ! update the auxiliary variables (runs through the equilibrium
  ! diss/precip/sorb routine) - we only need to update the aqueous_conc from
  ! the aqueous_mass value, or vice versa?

  call PMWellResidualTranAccum(this)

  ! calculate the source/sink terms (in/out of well segments)
  call PMWellResidualTranSrcSink(this)

  ! calculate the rxn terms (decay/ingrowth)
  call PMWellResidualTranRxn(this)

  ! calculate the flux terms
  call PMWellResidualTranFlux(this)

end subroutine PMWellResidualTran

! ************************************************************************** !

subroutine PMWellResidualTranAccum(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  implicit none

  class(pm_well_type) :: this

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend
  PetscReal :: Res(this%nspecies)

  ! porosity in [m^3-void/m^3-bulk]
  ! saturation in [m^3-liq/m^3-void]
  ! volume in [m^3-bulk]
  ! aqueous conc in [mol-species/m^3-liq]
  ! residual in [mol-species/sec]

  ! calculate the accumulation term as:
  ! residual = (Res_accum/tran_dt)
  !            - residual(which is already = fixed_accum/dt)

  do isegment = 1,this%well_grid%nsegments

    offset = (isegment-1)*this%nspecies
    istart = offset + 1
    iend = offset + this%nspecies

    do ispecies = 1,this%nspecies
      k = ispecies
      Res(k) = this%well%volume(isegment) * this%well%phi(isegment) * &
               this%well%liq%s(isegment) * &
               this%well%aqueous_conc(ispecies,isegment)
      Res(k) = Res(k) / this%dt_tran
    enddo

    this%tran_soln%residual(istart:iend) = Res(:) - &
                                         this%tran_soln%residual(istart:iend)
  enddo

end subroutine PMWellResidualTranAccum

! ************************************************************************** !

subroutine PMWellResidualTranSrcSink(this)
  !
  ! Calculates the source sink terms (Q in/out of well) for the transport
  ! residual equation.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/06/2022

  implicit none

  class(pm_well_type) :: this

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend
  PetscReal :: Res(this%nspecies)

  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: resr
  PetscReal :: coef_Qin, coef_Qout ! into well, out of well
  PetscReal :: Qin, Qout
  PetscReal :: rho_avg

  ! Q src/sink is in [kmol-liq/sec]
  ! FMWH2O is in [kg-liq/kmol-liq] where liq = water
  ! density is in [kg-liq/m^3-liq] where liq = water
  ! aqueous conc in [mol-species/m^3-liq]
  ! residual in [mol-species/sec]

  ! From the flow solution:
  ! + Q goes into well from reservoir       
  ! - Q goes out of well into reservoir     

  well => this%well
  resr => this%reservoir

  do isegment = 1,this%well_grid%nsegments

    rho_avg = 0.5d0*(well%liq%rho(isegment)+resr%rho_l(isegment))
    ! units of coef = [m^3-liq/sec]
    if (well%liq%Q(isegment) < 0.d0) then ! Q out of well
      coef_Qin = 0.d0
      coef_Qout = well%liq%Q(isegment)*FMWH2O/rho_avg
    else ! Q into well
    !            [kmol-liq/sec]*[kg-liq/kmol-liq]/[kg-liq/m^3-liq]  
      coef_Qin = well%liq%Q(isegment)*FMWH2O/rho_avg
      coef_Qout = 0.d0
    endif

    offset = (isegment-1)*this%nspecies
    istart = offset + 1
    iend = offset + this%nspecies

    do ispecies = 1,this%nspecies
      k = ispecies
      Qin = coef_Qin*resr%aqueous_conc(ispecies,isegment)
      Qout = coef_Qout*well%aqueous_conc(ispecies,isegment)
      Res(k) = Qin + Qout
    enddo

    this%tran_soln%residual(istart:iend) = &
                          this%tran_soln%residual(istart:iend) + Res(:)
  enddo

end subroutine PMWellResidualTranSrcSink

! ************************************************************************** !

subroutine PMWellResidualTranRxn(this)
  !
  ! Calculates the decay/ingrowth of radioactive species for the transport
  ! residual equation.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/06/2022

  implicit none

  class(pm_well_type) :: this

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend, parent_id
  PetscReal :: Res(this%nspecies)

  ! decay_rate in [1/sec]
  ! aqueous mass in [mol-species]
  ! residual in [mol-species/sec]

  do isegment = 1,this%well_grid%nsegments

    offset = (isegment-1)*this%nspecies
    istart = offset + 1
    iend = offset + this%nspecies

    do ispecies = 1,this%nspecies
      k = ispecies
      ! Add in species decay
      Res(k) = -(this%well%species_decay_rate(k)* &
                          this%well%aqueous_mass(ispecies,isegment))
      ! Add in contribution from parent (if exists)
      parent_id = this%well%species_parent_id(ispecies)
      if (parent_id > 0) then
        Res(k) = Res(k) + (this%well%species_parent_decay_rate(k)* &
                 this%well%aqueous_mass(parent_id,isegment))
      endif
    enddo

    this%tran_soln%residual(istart:iend) = &
                          this%tran_soln%residual(istart:iend) - Res(:)
  enddo

end subroutine PMWellResidualTranRxn

! ************************************************************************** !

subroutine PMWellResidualTranFlux(this)
  !
  ! Calculates the interior and BC flux terms for the transport residual
  ! equation.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/24/2022

  implicit none

  class(pm_well_type) :: this

  PetscInt :: ispecies, isegment, k
  PetscInt :: offset, istart, iend
  PetscInt :: n_up, n_dn
  PetscReal :: area_up, area_dn
  PetscReal :: q_up, q_dn, q_test
  PetscReal :: conc
  PetscReal :: diffusion
  PetscReal :: Res(this%nspecies)
  PetscReal :: Res_up(this%nspecies), Res_dn(this%nspecies)

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

  do isegment = 2,(this%well_grid%nsegments-1)

    Res(:) = 0.d0; Res_up(:) = 0.d0; Res_dn(:) = 0.d0

    area_up = 0.5d0 * (this%well%area(isegment) + this%well%area(isegment+1))
    area_dn = 0.5d0 * (this%well%area(isegment) + this%well%area(isegment-1))

    q_up = this%well%ql(isegment)
    q_dn = this%well%ql(isegment-1)

    offset = (isegment-1)*this%nspecies
    istart = offset + 1
    iend = offset + this%nspecies

    do ispecies = 1,this%nspecies
      k = ispecies

      ! north surface:
      if (q_up < 0.d0) then ! flow is down well
        conc = this%well%aqueous_conc(k,isegment+1)
        Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
      elseif (q_up > 0.d0) then ! flow is up well
        conc = this%well%aqueous_conc(k,isegment)
        Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
      else ! q_up = 0
        Res_up(k) = (n_up*area_up)*(0.d0 - diffusion)
      endif

      ! south surface:
      if (q_dn < 0.d0) then ! flow is down well
        conc = this%well%aqueous_conc(k,isegment)
        Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
      elseif (q_dn > 0.d0) then ! flow up well
        conc = this%well%aqueous_conc(k,isegment-1)
        Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
      else ! q_dn = 0
        Res_dn(k) = (n_dn*area_dn)*(0.d0 - diffusion)
      endif

      Res(k) = Res_up(k) + Res_dn(k)
    enddo

    this%tran_soln%residual(istart:iend) = &
                      this%tran_soln%residual(istart:iend) + Res(:)
  enddo

  ! ----------------------------------------BOUNDARY-FLUXES------------------

  ! ----- bottom of well -----
  isegment = 1
  Res(:) = 0.d0; Res_up(:) = 0.d0; Res_dn(:) = 0.d0

  area_up = 0.5d0 * (this%well%area(isegment) + this%well%area(isegment+1))
  area_dn = this%well%area(isegment)

  q_up = this%well%ql(isegment)
  q_dn = this%well%ql_bc(1) ! bottom of hole ql

  offset = (isegment-1)*this%nspecies ! = 0
  istart = offset + 1
  iend = offset + this%nspecies

  do ispecies = 1,this%nspecies
    k = ispecies

    ! north surface:
    if (q_up < 0.d0) then ! flow is down the well
      conc = this%well%aqueous_conc(k,isegment+1)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    elseif (q_up > 0.d0) then ! flow is up the well
      conc = this%well%aqueous_conc(k,isegment)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    else ! q_up = 0
      Res_up(k) = (n_up*area_up)*(0.d0 - diffusion)
    endif

    ! south surface:
    if (q_dn < 0.d0) then ! flow is down the well
      conc = this%well%aqueous_conc(k,isegment)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    elseif (q_dn > 0.d0) then ! flow is up the well
      conc = this%reservoir%aqueous_conc(k,isegment)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    else ! q_dn = 0
      Res_dn(k) = (n_dn*area_dn)*(0.d0 - diffusion)
    endif

    Res(k) = Res_up(k) + Res_dn(k)
  enddo

  this%tran_soln%residual(istart:iend) = &
                      this%tran_soln%residual(istart:iend) + Res(:)


  ! ----- top of well -----
  isegment = this%well_grid%nsegments
  Res(:) = 0.d0; Res_up(:) = 0.d0; Res_dn(:) = 0.d0

  area_up = this%well%area(isegment)
  area_dn = 0.5d0 * (this%well%area(isegment) + this%well%area(isegment-1))

  q_up = this%well%ql_bc(2) ! top of hole ql
  q_dn = this%well%ql(isegment-1)

  offset = (isegment-1)*this%nspecies
  istart = offset + 1
  iend = offset + this%nspecies

  do ispecies = 1,this%nspecies
    k = ispecies

    ! north surface:
    if (q_up < 0.d0) then ! flow is down the well
      conc = this%well%aqueous_conc_th(ispecies)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    elseif (q_up > 0.d0) then ! flow is up the well
      conc = this%well%aqueous_conc(k,isegment)
      Res_up(k) = (n_up*area_up)*(q_up*conc - diffusion)
    else ! q_up = 0
      Res_up(k) = (n_up*area_up)*(0.d0 - diffusion)
    endif

    ! south surface:
    if (q_dn < 0.d0) then ! flow is down the well
      conc = this%well%aqueous_conc(k,isegment)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    elseif (q_dn > 0.d0) then ! flow is up the well
      conc = this%well%aqueous_conc(k,isegment-1)
      Res_dn(k) = (n_dn*area_dn)*(q_dn*conc - diffusion)
    else ! q_dn = 0
      Res_dn(k) = (n_dn*area_dn)*(0.d0 - diffusion)
    endif

    Res(k) = Res_up(k) + Res_dn(k)
  enddo

  this%tran_soln%residual(istart:iend) = &
                      this%tran_soln%residual(istart:iend) + Res(:)

end subroutine PMWellResidualTranFlux

! ************************************************************************** !

subroutine PMWellJacobianFlow(this)
  !
  ! Author: Michael Nole
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this

  PetscInt :: local_id
  PetscInt :: iconn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: i,k
  Vec, parameter :: null_vec = tVec(0)
  PetscReal :: Jup(this%nphase,this%nphase), &
               Jdn(this%nphase,this%nphase), &
               Jtop(this%nphase,this%nphase), &
               Jbtm(this%nphase,this%nphase), &
               Jtmp(this%nphase,this%nphase), &
               Jac(this%nphase*this%well_grid%nsegments, &
                   this%nphase*this%well_grid%nsegments)

  this%flow_soln%Jacobian = 0.d0
  Jac = 0.d0
  Jup = 0.d0
  Jtop = 0.d0
  Jbtm = 0.d0
  Jtmp = 0.d0

  select case(this%well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('WIPP_DARCY')
      ! Not expecting ghosting at this time (1D model)
      !if (.not. well_analytical_derivatives) then
      !  call PMWellPerturb(this)
      !endif

      ! Accumulation Term ------------------------------------
      do local_id = 1,this%well_grid%nsegments
        call PMWellAccumDerivative(this,local_id,Jup)
        call PMWellFillJacFlow(this,Jac,Jup,local_id,local_id)
      enddo

      ! Source/Sink Term
      do local_id = 1,this%well_grid%nsegments
        call PMWellSrcSinkDerivative(this,local_id,Jup)
        call PMWellFillJacFlow(this,Jac,Jup,local_id,local_id)
      enddo

      ! Interior Flux Terms -----------------------------------
      do iconn = 1,this%well_grid%nconnections

        local_id_up = iconn
        local_id_dn = iconn+1

        call PMWellFluxDerivative(this,local_id_up,local_id_dn,Jup,Jdn)

        Jtmp = Jup
        call PMWellFillJacFlow(this,Jac,Jtmp,local_id_up,local_id_up)

        Jtmp = Jdn
        call PMWellFillJacFlow(this,Jac,Jtmp,local_id_up,local_id_dn)

        Jup = -Jup
        Jdn = -Jdn
        Jtmp = Jdn
        call PMWellFillJacFlow(this,Jac,Jtmp,local_id_dn,local_id_dn)

        Jtmp = Jup
        call PMWellFillJacFlow(this,Jac,Jtmp,local_id_dn,local_id_up)

      enddo

      ! Boundary Flux Terms -----------------------------------
      local_id = 1
      call PMWellBCFluxDerivative(this,Jtop,Jbtm)
      Jbtm = -Jbtm
      call PMWellFillJacFlow(this,Jac,Jbtm,local_id,local_id)

      local_id = this%well_grid%nsegments
      Jtop = -Jtop
      call PMWellFillJacFlow(this,Jac,Jtop,local_id,local_id)

      !pm_well_ni_count = pm_well_ni_count + 1

    case('FULL_MOMENTUM')

  end select

  do i = 1,this%nphase*this%well_grid%nsegments
    do k = 1,this%nphase*this%well_grid%nsegments
      this%flow_soln%Jacobian(i,k) = Jac(i,k)
    enddo
  enddo

  !this%flow_soln%Jacobian = Jac

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

subroutine PMWellJacobianTran(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/14/2022
  !

  implicit none

  class(pm_well_type) :: this

  PetscReal :: Jblock(this%nspecies,this%nspecies)
  PetscInt :: k, nspecies
  PetscInt :: jstart, jend

  nspecies = this%nspecies
  this%tran_soln%Jacobian(:,:) = 0.d0

  do k = 1,this%well_grid%nsegments

    Jblock(:,:) = 0.d0

    call PMWellJacTranAccum(this,Jblock,k)

    call PMWellJacTranSrcSink(this,Jblock,k)

    call PMWellJacTranFlux(this,Jblock,k)

    call PMWellJacTranRxn(this,Jblock,k)

    ! place JBlock into full Jac based on isegment
    jstart = (k-1)*nspecies + 1
    jend = jstart + nspecies - 1
    this%tran_soln%Jacobian(jstart:jend,jstart:jend) = Jblock

  enddo

end subroutine PMWellJacobianTran

! ************************************************************************** !

subroutine PMWellJacTranAccum(this,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/14/2022
  !

  implicit none

  class(pm_well_type) :: this
  PetscReal :: Jblock(this%nspecies,this%nspecies)
  PetscInt :: isegment

  PetscReal :: vol_dt
  PetscInt :: istart, iend, ispecies

  ! units of Jac = [m^3-bulk/sec]
  ! units of volume = [m^3-bulk]
  ! units of tran_dt = [sec]

  vol_dt = this%well%volume(isegment)/this%dt_tran

  istart = 1
  iend = this%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + vol_dt

  enddo

end subroutine PMWellJacTranAccum

! ************************************************************************** !

subroutine PMWellJacTranSrcSink(this,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/2022
  !

  implicit none

  class(pm_well_type) :: this
  PetscReal :: Jblock(this%nspecies,this%nspecies)
  PetscInt :: isegment

  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: resr
  PetscInt :: istart, iend, ispecies
  PetscReal :: Qin, Qout
  PetscReal :: SSin, SSout, SS
  PetscReal :: vol, rho_avg 

  ! units of Jac = [m^3-bulk/sec]
  ! units of volume = [m^3-bulk]
  ! units of liq%Q = [kmol-liq/sec]
  ! units of FMWH2O = [kg-liq/kmol-liq] 
  ! units of density = [kg-liq/m^3-liq] 
  ! units of Qin = [m^3-liq/sec]
  ! units of SS = [m^3-bulk/sec]

  well => this%well
  resr => this%reservoir

  ! From the flow solution:
  ! + Q goes into well from reservoir
  ! - Q goes out of well into reservoir

  vol = this%well%volume(isegment)
  rho_avg = 0.5d0*(well%liq%rho(isegment)+resr%rho_l(isegment))

  ! units of Qin/out = [m^3-liq/sec]
  if (well%liq%Q(isegment) < 0.d0) then ! Q out of well
    Qin = 0.d0
    Qout = well%liq%Q(isegment)*FMWH2O/rho_avg
    if (well%liq%s(isegment) < 1.d-40) then
      this%option%io_buffer = 'HINT: The liquid saturation is zero. &
        &Division by zero will occur in PMWellJacTranSrcSink().'
      call PrintMsg(this%option)
    endif
  else ! Q into well
    Qin = well%liq%Q(isegment)*FMWH2O/rho_avg
    Qout = 0.d0
    if (resr%s_l(isegment) < 1.d-40) then
      this%option%io_buffer = 'HINT: The liquid saturation is zero. &
        &Division by zero will occur in PMWellJacTranSrcSink().'
      call PrintMsg(this%option)
    endif
  endif

  SSin = Qin / (resr%e_por(isegment)*resr%s_l(isegment))    ! [m3-bulk/sec]
  SSout = Qout / (well%phi(isegment)*well%liq%s(isegment))  ! [m3-bulk/sec]
  SS = SSin + SSout

  istart = 1
  iend = this%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + vol*(SS/vol)

  enddo

end subroutine PMWellJacTranSrcSink

! ************************************************************************** !

subroutine PMWellJacTranFlux(this,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/2022
  !

  implicit none

  class(pm_well_type) :: this
  PetscReal :: Jblock(this%nspecies,this%nspecies)
  PetscInt :: isegment

  type(well_type), pointer :: well
  PetscInt :: istart, iend, ispecies
  PetscInt :: n_up, n_dn
  PetscReal :: d_diffusion_dM
  PetscReal :: J_up, J_dn
  PetscReal :: area_up, area_dn
  PetscReal :: sat_up, sat_dn, por_up, por_dn
  PetscReal :: u_up, u_dn, q_test

  ! units of Jac = [m^3-bulk/sec]
  ! area in [m2-bulk]
  ! q in [m3-liq/m2-bulk-sec]
  ! u in [m-liq/sec]
  ! sat in [m2-liq/m2-void] 

  well => this%well

  ! NOTE: The up direction is towards well top, and the dn direction is
  !       towards the well bottom.
  !       +q flows up the well
  !       -q flows down the well

  n_dn = +1
  n_up = -1

  d_diffusion_dM = 0.d0 ! for now, since WIPP has no diffusion

  if ((isegment > 1) .and. (isegment < this%well_grid%nsegments)) then
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

  else if (isegment == this%well_grid%nsegments) then
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
  iend = this%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + (J_up + J_dn)

  enddo

end subroutine PMWellJacTranFlux

! ************************************************************************** !

subroutine PMWellJacTranRxn(this,Jblock,isegment)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/2022
  !

  implicit none

  class(pm_well_type) :: this
  PetscReal :: Jblock(this%nspecies,this%nspecies)
  PetscInt :: isegment

  PetscInt :: istart, iend, ispecies, parent_id
  PetscReal :: vol

  ! units of Jac = [m^3-bulk/sec]
  ! units of volume = [m^3-bulk]
  ! decay_rate in [1/sec]

  vol = this%well%volume(isegment)

  istart = 1
  iend = this%nspecies

  do ispecies = istart,iend

    Jblock(ispecies,ispecies) = Jblock(ispecies,ispecies) + &
                               (vol*this%well%species_decay_rate(ispecies))

    parent_id = this%well%species_parent_id(ispecies)
    if (parent_id > 0) then
      Jblock(ispecies,parent_id) = Jblock(ispecies,parent_id) - &
                        (vol*this%well%species_parent_decay_rate(ispecies))
    endif

  enddo

end subroutine PMWellJacTranRxn

! ************************************************************************** !

subroutine PMWellPreSolve(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this

  ! pseudo-code for solve order of operations:
  !
  ! call pm%solve()
  !   call pm_well%pmwellsolve()
  !     call pm_well%pmwellsolveflow()
  !       call pm_well%pmwellpresolveflow()
  !       flow solve occurs
  !       call pm_well%pmwellpostsolveflow()
  !     call pm_well%pmwellsolvetran()       ---> if (transport)
  !       call pm_well%pmwellpresolvetran()  ---> if (transport)
  !       transport solve occurs             ---> if (transport)
  !       call pm_well%pmwellpostsolvetran() ---> if (transport)


end subroutine PMWellPreSolve

! ************************************************************************** !

subroutine PMWellPreSolveFlow(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: cur_time, cur_time_converted
  PetscReal :: dt_converted

  this%flow_soln%not_converged = PETSC_TRUE
  this%flow_soln%converged = PETSC_FALSE

  cur_time = this%option%time - this%option%flow_dt + this%cumulative_dt_flow
  cur_time_converted = cur_time/this%output_option%tconv
  dt_converted = this%dt_flow/this%output_option%tconv

  write(out_string,'(" FLOW Step ",i6,"   Time =",1pe12.5,"   Dt =", &
                     1pe12.5," ",a4)') &
                   (this%flow_soln%n_steps+1),cur_time_converted, &
                   dt_converted,this%output_option%tunit
  call PrintMsg(this%option,out_string)

end subroutine PMWellPreSolveFlow

! ************************************************************************** !

subroutine PMWellPreSolveTran(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/22/2022

  implicit none

  class(pm_well_type) :: this

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: cur_time, cur_time_converted
  PetscReal :: dt_converted

  this%tran_soln%not_converged = PETSC_TRUE
  this%tran_soln%converged = PETSC_FALSE

  cur_time = this%option%time - this%option%flow_dt + this%cumulative_dt_tran
  cur_time_converted = cur_time/this%output_option%tconv
  dt_converted = this%dt_tran/this%output_option%tconv

  write(out_string,'(" TRAN Step ",i6,"   Time =",1pe12.5,"   Dt =", &
                     1pe12.5," ",a4)') &
                   (this%tran_soln%n_steps+1),cur_time_converted, &
                   dt_converted,this%output_option%tunit
  call PrintMsg(this%option,out_string)

end subroutine PMWellPreSolveTran

! ************************************************************************** !

subroutine PMWellSolve(this,time,ierr)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

  implicit none

  class(pm_well_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: curr_time

  curr_time = this%option%time - this%option%flow_dt

  if (Initialized(this%intrusion_time_start) .and. &
      (curr_time < this%intrusion_time_start)) then
    write(out_string,'(" Inactive.    Time =",1pe12.5," sec.")') curr_time
    call PrintMsg(this%option,out_string)
    ierr = 0 ! If this is not set to zero, TS_STOP_FAILURE occurs!
    return
  endif

  call PMWellSolveFlow(this,time,ierr)

  if (this%transport) then
    call PMWellSolveTran(this,time,ierr)
  endif

end subroutine PMWellSolve

! ************************************************************************** !

subroutine PMWellSolveFlow(this,time,ierr)
  !
  ! Author: Michael Nole
  ! Date: 12/01/2021

  implicit none

  class(pm_well_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(well_soln_flow_type), pointer :: flow_soln
  character(len=MAXSTRINGLENGTH) :: out_string
  PetscLogDouble :: log_start_time, log_end_time
  PetscInt :: n_iter,ts_cut,easy_converge_count
  PetscInt :: istart, iend
  PetscReal :: res(this%flow_soln%ndof)
  PetscReal :: res_fixed(this%flow_soln%ndof*this%well_grid%nsegments)
  PetscReal :: Q_liq(this%well_grid%nsegments,this%well_grid%nsegments), &
               Q_gas(this%well_grid%nsegments,this%well_grid%nsegments)
  PetscReal :: v_darcy
  PetscBool :: steady_state, upwind
  PetscReal :: ss_check_p(this%well_grid%nsegments,2), &
               ss_check_s(this%well_grid%nsegments,2)
  PetscReal :: dpdt(this%well_grid%nsegments), dsdt(this%well_grid%nsegments)
  PetscReal :: den_kg_ave, mobility, delta_pressure, perm_factor, dx, well_perm
  PetscReal :: gravity_term, area, mass_conserved_liq, mass_conserved_gas
  PetscInt :: i, j, k
  PetscInt :: ss_step_count, steps_to_declare_ss

  flow_soln => this%flow_soln

  ierr = 0
  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

  ts_cut = 0
  easy_converge_count = 0

  this%cumulative_dt_flow = 0.d0
  flow_soln%converged = PETSC_FALSE
  flow_soln%not_converged = PETSC_TRUE

  ss_check_p(:,1) = this%well%pl(:)
  ss_check_s(:,1) = this%well%gas%s(:)
  steady_state = PETSC_FALSE
  ss_step_count = 0
  steps_to_declare_ss = 10

  ! update well index (should not need to be updated every time if
  ! grid permeability and discretization are static
  call PMWellComputeWellIndex(this)

  if (this%well%well_model_type == 'STEADY_STATE') then

    call PMWellComputeWellIndex(this)
    area = 0.d0
    v_darcy = 0.d0
    Q_liq = 0.d0
    Q_gas = 0.d0
    this%well%liq%Q = 0.d0
    this%well%gas%Q = 0.d0
    do i = 1,this%well_grid%nsegments
      if (this%well%WI(i) == 0) cycle
      do j = i,this%well_grid%nsegments

        if (i==j) cycle
        if (this%well%WI(j) == 0) cycle

        area = pi*0.5d0*(this%well%diameter(i)*this%well_grid%dh(i) + &
                         this%well%diameter(j)*this%well_grid%dh(j))
        well_perm = this%well%permeability(i)
        dx = this%well_grid%dh(i)/2.d0
        if ((this%well%permeability(i) > this%reservoir%kz(i)) .and. &
            (this%well%permeability(j) > this%reservoir%kz(j)) ) then
          do k = i+1,j
            if (k==j) then
              well_perm = (dx+this%well_grid%dh(k)/2.d0)/ &
                          (dx/well_perm + &
                          this%well_grid%dh(k)/(2.d0*this%well%permeability(k)))
              dx = dx+this%well_grid%dh(k)
            else
              well_perm = (dx+this%well_grid%dh(k))/ &
                          (dx/well_perm + &
                           this%well_grid%dh(k)/this%well%permeability(k))
              dx = dx+this%well_grid%dh(k)
            endif
          enddo

          ! Take the harmonic mean of well indicies, factoring out dh
          perm_factor = (dx + this%well%r0(i) + this%well%r0(j)) / &
                      (this%well%r0(i)/(this%well%WI(i)/this%well_grid%dh(i))+ &
                       this%well%r0(j)/(this%well%WI(j)/this%well_grid%dh(j))+ &
                       dx / (well_perm))
          dx = dx + this%well%r0(i) + this%well%r0(j)

          ! Liquid Phase

          !den_kg_ave = sum(this%reservoir%rho_l(i:j))/(i-j+1)
          !gravity_term = 0.d0
          !do k = i+1,j
          !  gravity_term = gravity_term + (this%reservoir%rho_l(k) + &
          !                 this%reservoir%rho_l(k-1))/2.d0 * gravity * &
          !                 dabs(this%well_grid%h(k)%z-this%well_grid%h(k-1)%z)
          !enddo

          ! Gravity term might get weird if you have multiple well cells per
          ! reservoir cell, so with this method I think you have to only allow
          ! 1 well cell per reservoir cell.
          ! Gravity will be negative for z upward
          den_kg_ave = (this%reservoir%rho_l(i)+this%reservoir%rho_l(j))/2.d0
          gravity_term = den_kg_ave * gravity * &
                         dabs(this%well_grid%h(j)%z-this%well_grid%h(i)%z)
          ! dP = Plow - rho * g * z - Phigh
          delta_pressure = this%reservoir%p_l(i) + gravity_term  - &
                           this%reservoir%p_l(j)
          upwind = delta_pressure > 0
          if (upwind) then
            ! Flow is upward
            mobility = this%reservoir%kr_l(i)/this%reservoir%visc_l(i)
          else
            ! Flow is downward
            mobility = this%reservoir%kr_l(j)/this%reservoir%visc_l(j)
          endif

          ! Flowrate in kg/s: Positive is from i to j (upward)
          Q_liq(i,j) = den_kg_ave*mobility*delta_pressure/dx*perm_factor*area
          Q_liq(j,i) = -1.d0 * Q_liq(i,j)


          ! Gas Phase
          den_kg_ave = (this%reservoir%rho_g(i)+this%reservoir%rho_g(j))/2.d0
          ! Gravity term might get weird if you have multiple well cells per
          ! reservoir cell, so with this method I think you have to only allow
          ! 1 well cell per reservoir cell.
          ! Gravity will be negative for z upward
          gravity_term = den_kg_ave * gravity * &
                         dabs(this%well_grid%h(j)%z-this%well_grid%h(i)%z)
          ! dP = Plow - rho * g * z - Phigh
          delta_pressure = this%reservoir%p_g(i) + gravity_term  - &
                           this%reservoir%p_g(j)
          upwind = delta_pressure > 0
          if (upwind) then
            ! Flow is upward
            mobility = this%reservoir%kr_g(i)/this%reservoir%visc_g(i)
          else
            ! Flow is downward
            mobility = this%reservoir%kr_g(j)/this%reservoir%visc_g(j)
          endif

          ! Flowrate in kg/s: Positive is from i to j (upward)
          Q_gas(i,j) = den_kg_ave*mobility*delta_pressure/dx*perm_factor*area
          Q_gas(j,i) = -1.d0 * Q_gas(i,j)

        else
          Q_liq(i,j) = 0.d0
          Q_liq(j,i) = 0.d0
          Q_gas(i,j) = 0.d0
          Q_gas(j,i) = 0.d0
        endif

          this%well%liq%Q(i) = this%well%liq%Q(i) + Q_liq(i,j)
          this%well%liq%Q(j) = this%well%liq%Q(j) + Q_liq(j,i)

          this%well%gas%Q(i) = this%well%gas%Q(i) + Q_gas(i,j)
          this%well%gas%Q(j) = this%well%gas%Q(j) + Q_gas(j,i) 
      enddo

      if (this%flow_soln%th_p) then
        well_perm = this%well%permeability(i)
        dx = this%well_grid%dh(i)/2.d0
        !Compute the effective permeability of the segment
        !from cell i to the top of the well
        if (i /= this%well_grid%nsegments) then
          do k = i+1,this%well_grid%nsegments
            well_perm = (dx+this%well_grid%dh(k))/ &
                        (dx/well_perm + &
                         this%well_grid%dh(k)/this%well%permeability(k))
            dx = dx+this%well_grid%dh(k)
          enddo
        endif
        !Just the well index of cell i over dh(i), since flux is out
        !the top of the well
        perm_factor = this%well%WI(i)/this%well_grid%dh(i)

        !Average the surface area and cross-sectional area of the well
        !(flux is into/out of the well in the radial direction at 
        !cell i, and into/out of the well in the axial direction at top.
        area = pi*0.5*(this%well%diameter(i)*this%well_grid%dh(i) + &
                       (this%well%diameter(i)/2.d0)**2)

        !Liquid Phase

        !Take average density between cell i and top cell
        den_kg_ave = (this%reservoir%rho_l(i)+ &
                      this%reservoir%rho_l(this%well_grid%nsegments))/2.d0

        gravity_term = den_kg_ave * gravity * &
                       dabs(this%well_grid%h(i)%z- &
                            (this%well_grid%h(this%well_grid%nsegments)%z + &
                             this%well_grid%dh(this%well_grid%nsegments)/2.d0))
        ! dP = Plow - rho * g * z - Phigh
        delta_pressure = this%reservoir%p_l(i) + gravity_term  - &
                         this%well%th_p
        ! For flow out
        mobility = this%reservoir%kr_l(i)/this%reservoir%visc_l(i)

        this%well%liq%Q(i) = this%well%liq%Q(i) + &
                             den_kg_ave*mobility*delta_pressure/dx* &
                             well_perm * perm_factor*area

        !Gas Phase

        den_kg_ave = (this%reservoir%rho_g(i)+ &
                      this%reservoir%rho_g(this%well_grid%nsegments))/2.d0

        gravity_term = den_kg_ave * gravity * &
                       dabs(this%well_grid%h(i)%z- &
                            (this%well_grid%h(this%well_grid%nsegments)%z  + &
                             this%well_grid%dh(this%well_grid%nsegments)/2.d0))
        ! dP = Plow - rho * g * z - Phigh
        delta_pressure = this%reservoir%p_g(i) + gravity_term  - &
                         this%well%th_p

        ! For flow out
        mobility = this%reservoir%kr_g(i)/this%reservoir%visc_g(i)

        this%well%gas%Q(i) = this%well%gas%Q(i) + &
                             den_kg_ave*mobility*delta_pressure/dx* &
                             well_perm * perm_factor*area

      endif
    enddo

    ! Should equal the total flux out the top of the domain
    mass_conserved_liq = sum(this%well%liq%Q)
    mass_conserved_gas = sum(this%well%gas%Q)

    this%cumulative_dt_flow = this%realization%option%flow_dt

    ! Update transport

    call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

  endif


  do while (this%cumulative_dt_flow < this%realization%option%flow_dt)

    ! update the well src/sink Q vector at start of time step
    call PMWellUpdateWellQ(this%well,this%reservoir)

    call PMWellPreSolveFlow(this)

    ! Fixed accumulation term
    res_fixed = 0.d0
    res = 0.d0
    do i = 1,this%well_grid%nsegments
      call PMWellAccumulationFlow(this,this%well,i,res)
      istart = flow_soln%ndof*(i-1)+1
      iend = flow_soln%ndof*i
      res_fixed(istart:iend) = -1.d0 * res * this%dt_flow
    enddo

    n_iter = 0

    do while (flow_soln%not_converged)

      if (n_iter > (flow_soln%max_iter-1)) then
        flow_soln%cut_timestep = PETSC_TRUE
        out_string = ' Maximum number of FLOW Newton iterations reached. &
                      &Cutting timestep!'
        call PrintMsg(this%option,out_string); WRITE(*,*) ""
        call PMWellCutTimestepFlow(this)
        n_iter = 0
        ts_cut = ts_cut + 1
        easy_converge_count = 0

        if (ss_step_count > 2 .and. this%ss_check) then
          steady_state = PETSC_TRUE
          this%cumulative_dt_flow = this%realization%option%flow_dt
          WRITE(out_string,'(" PM Well FLOW convergence declared due to &
            &automatic time step control criterion. ")')
          call PrintMsg(this%option,out_string)
        endif


        exit
      endif
      if (ts_cut > flow_soln%max_ts_cut) then
        this%realization%option%io_buffer = &
          ' Maximum timestep cuts reached in PM Well FLOW. Solution has not &
           &converged. Exiting.'
        if (this%print_well) then
          call PMWellOutput(this)  
        endif
        call PrintErrMsg(this%realization%option)
      endif

      flow_soln%residual = 0.d0
      flow_soln%residual = res_fixed / this%dt_flow

      easy_converge_count = easy_converge_count + 1

      call PMWellNewtonFlow(this)

      call PMWellCheckConvergenceFlow(this,n_iter,res_fixed)

    enddo

    if (easy_converge_count > 4 ) then
      if (this%cumulative_dt_flow + this%dt_flow * &
          flow_soln%ts_cut_factor < this%realization%option%flow_dt) then
        this%dt_flow = this%dt_flow * flow_soln%ts_cut_factor
        flow_soln%cut_timestep = PETSC_FALSE
      endif
    endif

    if (this%cumulative_dt_flow + this%dt_flow > &
        this%realization%option%flow_dt) then
      this%dt_flow = this%realization%option%flow_dt - this%cumulative_dt_flow
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

    if (this%ss_check) then
      if (ss_step_count >= steps_to_declare_ss) then
        steady_state = PETSC_TRUE
        this%cumulative_dt_flow = this%realization%option%flow_dt
      endif
    endif

    ! Other way:
    !if (this%ss_check .and. flow_soln%converged) then
    !  ss_check_p(:,2) = this%well%pl(:)
    !  ss_check_s(:,2) = this%well%gas%s(:)

    !  dpdt = (ss_check_p(:,2) - ss_check_p(:,1)) / this%dt_flow
    !  dsdt = (ss_check_s(:,2) - ss_check_s(:,1)) / this%dt_flow

    !  if (maxval(abs(dpdt)) < eps_p) then
    !    if (maxval(abs(dsdt)) < eps_s) then
    !      ss_step_count = ss_step_count + 1
    !      if (ss_step_count > steps_to_declare_ss) steady_state = PETSC_TRUE
    !    else
    !      ss_step_count = 0
    !    endif
    !  endif
    !  if (steady_state) then
    !    this%cumulative_dt_flow = this%realization%option%flow_dt
    !  endif
    !  ss_check_p(:,1) = this%well%pl(:)
    !  ss_check_s(:,1) = this%well%gas%s(:)
    !endif

  enddo

  call PMWellPostSolveFlow(this)

  call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

end subroutine PMWellSolveFlow

! ************************************************************************** !

subroutine PMWellSolveTran(this,time,ierr)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/22/2022

  implicit none

  class(pm_well_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(well_soln_tran_type), pointer :: soln
  character(len=MAXSTRINGLENGTH) :: out_string
  PetscLogDouble :: log_start_time, log_end_time
  PetscReal :: res_fixed(this%tran_soln%ndof*this%well_grid%nsegments)
  PetscInt :: n_iter, ts_cut
  PetscInt :: istart, iend
  PetscInt :: k

  soln => this%tran_soln

  ierr = 0
  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

  ts_cut = 0

  this%cumulative_dt_tran = 0.d0
  soln%converged = PETSC_FALSE
  soln%not_converged = PETSC_TRUE

  do while (this%cumulative_dt_tran < this%realization%option%flow_dt)

    call PMWellPreSolveTran(this)

    n_iter = 0

    do while (soln%not_converged)
      if (n_iter > (soln%max_iter-1)) then
        soln%cut_timestep = PETSC_TRUE
        out_string = ' Maximum number of TRAN Newton iterations reached. &
                      &Cutting timestep!'
        call PrintMsg(this%option,out_string); WRITE(*,*) ""
        call PMWellCutTimestepTran(this)
        n_iter = 0
        ts_cut = ts_cut + 1
        exit
      endif
      if (ts_cut > soln%max_ts_cut) then
        this%realization%option%io_buffer = &
          ' Maximum timestep cuts reached in PM Well TRAN. Solution has not &
           &converged. Exiting.'
        if (this%print_well) then
          call PMWellOutput(this)
        endif
        call PrintErrMsg(this%realization%option)
      endif

      soln%residual = 0.d0
      ! Get fixed accumulation term (not yet divided by dt)
      do k = 1,this%well_grid%nsegments
        istart = soln%ndof*(k-1)+1
        iend = soln%ndof*k
        call PMWellAccumulationTran(this,k,res_fixed(istart:iend))
      enddo
      soln%residual = res_fixed / this%dt_tran

      call PMWellNewtonTran(this,n_iter)

      call PMWellCheckConvergenceTran(this,n_iter,res_fixed)

    enddo

    ! try to increase the time step, if possible
    if (soln%converged) then
      this%dt_tran = soln%ts_ramp_factor * this%dt_tran
    endif 
    if (this%dt_tran > this%realization%option%flow_dt) then
      this%dt_tran = this%realization%option%flow_dt
    endif

    ! if this next time step will overstep flow_dt, then correct it
    if (this%cumulative_dt_tran + this%dt_tran > &
        this%realization%option%flow_dt) then
      this%dt_tran = this%realization%option%flow_dt - this%cumulative_dt_tran
    endif

    soln%n_steps = soln%n_steps + 1

  enddo

  call PMWellPostSolveTran(this)

  call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

end subroutine PMWellSolveTran

! ************************************************************************** !

subroutine PMWellUpdateSolutionFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 01/21/22

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: i
  PetscInt :: idof
  PetscReal, parameter :: MIN_SAT = 1.d-6
  PetscReal, parameter :: MAX_SAT = 0.99999

  select case(pm_well%well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('WIPP_DARCY')
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
      call PMWellUpdatePropertiesFlow(pm_well,pm_well%well, &
                     pm_well%realization%patch%characteristic_curves_array, &
                     pm_well%realization%option)
  end select
end subroutine PMWellUpdateSolutionFlow

! ************************************************************************** !

subroutine PMWellUpdateSolutionTran(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/15/22

  implicit none

  class(pm_well_type) :: this

  PetscInt :: isegment
  PetscInt :: offset, istart, iend
  PetscInt :: nspecies 
  PetscReal :: vol

  ! update in [mol/m3-bulk]
  ! volume in [m3-bulk]

  nspecies = this%nspecies 

  do isegment = 1,this%well_grid%nsegments

    offset = (isegment-1)*nspecies
    istart = offset + 1
    iend = offset + nspecies

    vol = this%well%volume(isegment)

    this%well%aqueous_mass(1:nspecies,isegment) = &              ! [mol]
                           this%well%aqueous_mass(1:nspecies,isegment) + &
                           (this%tran_soln%update(istart:iend) * vol)
    this%well%aqueous_conc(1:nspecies,isegment) = &
          this%well%aqueous_mass(1:nspecies,isegment) / &        ! [mol]
          (this%well%phi(isegment)*vol*this%well%liq%s(isegment)) ! [m3-liq]

  enddo

end subroutine PMWellUpdateSolutionTran

! ************************************************************************** !

subroutine PMWellCutTimestepFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 01/24/2022

  implicit none

  class(pm_well_type) :: pm_well

  select case(pm_well%well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('WIPP_DARCY')
      ! could make this smarter or call smarter timestepping routines
      pm_well%dt_flow = pm_well%dt_flow / pm_well%flow_soln%ts_cut_factor
      pm_well%dt_flow = max(pm_well%dt_flow,pm_well%min_dt_flow)
      pm_well%well%pl = pm_well%flow_soln%prev_soln%pl
      pm_well%well%gas%s = pm_well%flow_soln%prev_soln%sg
      call PMWellUpdatePropertiesFlow(pm_well,pm_well%well, &
                     pm_well%realization%patch%characteristic_curves_array, &
                     pm_well%realization%option)
  end select

end subroutine PMWellCutTimestepFlow

! ************************************************************************** !

subroutine PMWellCutTimestepTran(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  implicit none

  class(pm_well_type) :: this

  this%dt_tran = this%dt_tran / this%tran_soln%ts_cut_factor
  this%dt_tran = max(this%dt_tran,this%min_dt_tran)
  this%well%aqueous_mass = this%tran_soln%prev_soln%aqueous_mass
  this%well%aqueous_conc = this%tran_soln%prev_soln%aqueous_conc
  call PMWellUpdatePropertiesTran(this)

end subroutine PMWellCutTimestepTran

! ************************************************************************** !

subroutine PMWellNewtonFlow(this)
  !
  ! Author: Michael Nole
  ! Date: 01/20/2022

  use Utility_module

  implicit none

  class(pm_well_type) :: this

  PetscReal :: identity(this%nphase*this%well_grid%nsegments,&
                        this%nphase*this%well_grid%nsegments)
  PetscReal :: new_dx(this%nphase*this%well_grid%nsegments)
  PetscInt :: indx(this%nphase*this%well_grid%nsegments)
  PetscInt :: i,j
  PetscInt :: d

  call PMWellUpdateWellQ(this%well,this%reservoir)

  call PMWellPerturb(this)

  call PMWellResidualFlow(this)

  call PMWellJacobianFlow(this)

  select case (this%well%well_model_type)
  !--------------------------------------
  case('CONSTANT_PRESSURE')
    ! No capillarity yet
    this%well%pl(:) = this%well%bh_p
    this%well%pg(:) = this%well%bh_p
  !--------------------------------------
  case('WIPP_DARCY')
    do i = 1,this%nphase*this%well_grid%nsegments
      do j = 1,this%nphase*this%well_grid%nsegments
        if (i==j) then
          identity(i,j) = 1.d0
        else
          identity(i,j) = 0.d0
        endif
      enddo
    enddo
    call LUDecomposition(this%flow_soln%Jacobian,this%nphase*this%well_grid% &
                         nsegments,indx,d)
    call LUBackSubstitution(this%flow_soln%Jacobian, &
                            this%nphase*this%well_grid%nsegments,&
                            indx,this%flow_soln%residual)
    new_dx = -1.d0 * this%flow_soln%residual
    this%flow_soln%update = new_dx

    call PMWellUpdateSolutionFlow(this)

  !--------------------------------------
  end select

end subroutine PMWellNewtonFlow

! ************************************************************************** !

subroutine PMWellNewtonTran(this,n_iter)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  use Utility_module

  implicit none

  class(pm_well_type) :: this
  PetscInt :: n_iter

  PetscInt :: nm, dummy
  PetscInt :: indx(this%nspecies*this%well_grid%nsegments)
  PetscReal :: res_fixed(this%nspecies*this%well_grid%nsegments)

  nm = this%nspecies * this%well_grid%nsegments

  ! at this time, the tran_soln%residual vector has been zero'd out and has
  ! been loaded with the fixed accumulation divided by current dt

  call PMWellResidualTran(this)

  call PMWellJacobianTran(this)

  ! J dx = -R     => dx = J^(-1)(-R)
  ! [m3-bulk/sec] dx = -[mol/sec]
  ! dx in [mol/m3-bulk]

  call LUDecomposition(this%tran_soln%Jacobian,nm,indx,dummy)

  call LUBackSubstitution(this%tran_soln%Jacobian,nm,indx, &
                          this%tran_soln%residual)

  this%tran_soln%update = +1.0d0 * this%tran_soln%residual ! [mol/m3-bulk]

  call PMWellUpdateSolutionTran(this)

end subroutine PMWellNewtonTran

! ************************************************************************** !

subroutine PMWellPostSolve(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this

  ! placeholder

end subroutine PMWellPostSolve

! ************************************************************************** !

subroutine PMWellPostSolveFlow(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: cur_time_converted

  cur_time_converted = this%option%time/this%output_option%tconv

  WRITE(out_string,'(" PM Well FLOW Step Complete!    Time=",1pe12.5," &
                    &",a4,"Total Newton Its =",i8)') &
                    cur_time_converted,this%output_option%tunit, &
                    this%flow_soln%n_newton
  call PrintMsg(this%option,out_string)
  WRITE(*,*) ""

end subroutine PMWellPostSolveFlow

! ************************************************************************** !

subroutine PMWellPostSolveTran(this)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/22/2022

  implicit none

  class(pm_well_type) :: this

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: cur_time_converted

  cur_time_converted = this%option%time/this%output_option%tconv

  WRITE(out_string,'(" PM Well TRAN Step Complete!    Time=",1pe12.5," &
                    &",a4,"Total Newton Its =",i8)') &
                    cur_time_converted,this%output_option%tunit, &
                    this%flow_soln%n_newton
  call PrintMsg(this%option,out_string)
  WRITE(*,*) ""

end subroutine PMWellPostSolveTran

! ************************************************************************** !

subroutine PMWellCheckConvergenceFlow(this,n_iter,fixed_accum)
  !
  ! Checks flow solution convergence against prescribed tolerances.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 01/20/2022

  implicit none

  class(pm_well_type) :: this
  PetscInt :: n_iter
  PetscReal :: fixed_accum(this%flow_soln%ndof*this%well_grid%nsegments)

  PetscReal, parameter :: zero_saturation = 1.d-15
  PetscReal, parameter :: zero_accumulation = 1.d-15

  type(well_soln_flow_type), pointer :: flow_soln
  character(len=MAXSTRINGLENGTH) :: out_string
  character(len=MAXSTRINGLENGTH) :: rsn_string
  PetscBool :: cnvgd_due_to_residual(this%well_grid%nsegments* &
                                     this%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_abs_res(this%well_grid%nsegments* &
                                    this%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_scaled_res(this%well_grid%nsegments* &
                                       this%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_update(this%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update_p(this%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update_s(this%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update(this%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update_p(this%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update_s(this%well_grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update(this%well_grid%nsegments)
  PetscBool :: cnvgd_on_pressure(this%well_grid%nsegments)
  PetscBool :: cnvgd_on_saturation(this%well_grid%nsegments)
  PetscReal :: update_p(this%well_grid%nsegments) ! liquid pressure
  PetscReal :: update_s(this%well_grid%nsegments) ! gas saturation
  PetscReal :: temp_real
  PetscReal :: max_scaled_residual,max_absolute_residual
  PetscReal :: max_relative_update_p,max_relative_update_s
  PetscReal :: max_absolute_update_p,max_absolute_update_s
  PetscInt :: loc_max_scaled_residual,loc_max_abs_residual
  PetscInt :: loc_max_rel_update_p,loc_max_rel_update_s
  PetscInt :: loc_max_abs_update_p,loc_max_abs_update_s
  PetscInt :: idof
  PetscInt :: k

  flow_soln => this%flow_soln

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
  flow_soln%residual = fixed_accum / this%dt_flow
  call PMWellResidualFlow(this)

  ! Update mass balance
  call PMWellMassBalance(this)

  do k = 1,this%well_grid%nsegments
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
    temp_real = dabs(update_p(k)/this%well%pl(k))
    if (temp_real < flow_soln%itol_rel_update_p) then
      cnvgd_due_to_rel_update_p(k) = PETSC_TRUE
    endif
    temp_real = dabs(update_s(k)/this%well%gas%s(k))
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
                       (fixed_accum(idof)/this%dt_flow))
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
                       (fixed_accum(idof+1)/this%dt_flow))
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
                                    (fixed_accum/this%dt_flow)))
  loc_max_scaled_residual = maxloc(dabs(flow_soln%residual/ &
                                        (fixed_accum/this%dt_flow)),1)

  max_absolute_update_p = maxval(dabs(update_p))
  loc_max_abs_update_p = maxloc(dabs(update_p),1)

  max_absolute_update_s = maxval(dabs(update_s))
  loc_max_abs_update_s = maxloc(dabs(update_s),1)

  max_relative_update_p = maxval(dabs(update_p/this%well%pl))
  loc_max_rel_update_p = maxloc(dabs(update_p/this%well%pl),1)

  max_relative_update_s = maxval(dabs(update_s/this%well%gas%s))
  loc_max_rel_update_s = maxloc(dabs(update_s/this%well%gas%s),1)

  do k = 1,this%well_grid%nsegments*flow_soln%ndof
    if (cnvgd_due_to_scaled_res(k) .or. cnvgd_due_to_abs_res(k)) then
      cnvgd_due_to_residual(k) = PETSC_TRUE
    endif
  enddo
  do k = 1,this%well_grid%nsegments
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

  write(out_string,'(i2," aR:",es10.2,"  sR:",es10.2,"  uP:" &
        &,es10.2,"  uS:",es10.2,"  ruP:",es10.2,"  ruS:",es10.2)') &
        n_iter,max_absolute_residual,max_scaled_residual, &
        max_absolute_update_p,max_absolute_update_s, &
        max_relative_update_p,max_relative_update_s
  call PrintMsg(this%option,out_string)

  if (all(cnvgd_due_to_residual) .and. all(cnvgd_due_to_update)) then
    flow_soln%converged = PETSC_TRUE
    flow_soln%not_converged = PETSC_FALSE
    out_string = ' FLOW Solution converged!  ---> ' // trim(rsn_string)
    call PrintMsg(this%option,out_string); WRITE(*,*) ""
    this%cumulative_dt_flow = this%cumulative_dt_flow + this%dt_flow
    this%flow_soln%prev_soln%pl = this%well%pl
    this%flow_soln%prev_soln%sg = this%well%gas%s
    !call PMWellUpdateSolutionFlow(this)
  else
    flow_soln%converged = PETSC_FALSE
    flow_soln%not_converged = PETSC_TRUE
  endif


end subroutine PMWellCheckConvergenceFlow

! ************************************************************************** !

subroutine PMWellCheckConvergenceTran(this,n_iter,fixed_accum)
  !
  ! Checks transport solution convergence against prescribed tolerances.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022

  implicit none

  class(pm_well_type) :: this
  PetscInt :: n_iter
  PetscReal :: fixed_accum(this%tran_soln%ndof*this%well_grid%nsegments)

  type(well_soln_tran_type), pointer :: soln
  character(len=MAXSTRINGLENGTH) :: out_string
  character(len=MAXSTRINGLENGTH) :: rsn_string
  PetscReal :: temp_real
  PetscBool :: cnvgd_due_to_residual(this%well_grid%nsegments* &
                                     this%tran_soln%ndof)
  PetscBool :: cnvgd_due_to_abs_res(this%well_grid%nsegments* &
                                    this%tran_soln%ndof)
  PetscBool :: cnvgd_due_to_scaled_res(this%well_grid%nsegments* &
                                       this%tran_soln%ndof)
  PetscBool :: cnvgd_due_to_update(this%well_grid%nsegments*this%tran_soln%ndof)
  PetscReal :: vol_vec(this%well_grid%nsegments*this%tran_soln%ndof)
  PetscReal :: aq_mass_vec(this%well_grid%nsegments*this%tran_soln%ndof)
  PetscReal :: max_scaled_residual,max_absolute_residual
  PetscReal :: max_update
  PetscInt :: loc_max_scaled_residual,loc_max_abs_residual
  PetscInt :: loc_max_update
  PetscInt :: k, n, j
  PetscInt :: isegment, ispecies

  soln => this%tran_soln

  n_iter = n_iter + 1
  soln%n_newton = soln%n_newton + 1

  cnvgd_due_to_residual = PETSC_FALSE
  cnvgd_due_to_abs_res = PETSC_FALSE
  cnvgd_due_to_scaled_res = PETSC_FALSE
  cnvgd_due_to_update = PETSC_FALSE
  rsn_string = ''

  ! Update the residual
  soln%residual = 0.d0
  soln%residual = fixed_accum/this%dt_tran 
  call PMWellResidualTran(this)

  do k = 1,(this%well_grid%nsegments*soln%ndof)
    ! Absolute Residual
    temp_real = dabs(soln%residual(k))
    if (temp_real < soln%itol_abs_res) then
      cnvgd_due_to_abs_res(k) = PETSC_TRUE
    endif
    ! Scaled Residual
    temp_real = dabs(soln%residual(k)/(fixed_accum(k)/this%dt_tran))
    if (temp_real < soln%itol_scaled_res) then
      cnvgd_due_to_scaled_res(k) = PETSC_TRUE
    endif
  enddo

  ! Relative Update
  do n = 1,this%well_grid%nsegments
    isegment = n
    do k = 1, soln%ndof
      ispecies = k
      j = ((isegment-1)*soln%ndof) + ispecies
      vol_vec(j) = this%well%volume(isegment)
      aq_mass_vec(j) = this%well%aqueous_mass(ispecies,isegment)
      temp_real = dabs(soln%update(j)*vol_vec(j)/aq_mass_vec(j))
      if (temp_real < soln%itol_rel_update) then
        cnvgd_due_to_update(j) = PETSC_TRUE
      endif
    enddo
  enddo

  max_absolute_residual = maxval(dabs(soln%residual))
  loc_max_abs_residual = maxloc(dabs(soln%residual),1)

  max_scaled_residual = maxval(dabs(soln%residual/ &
                                    (fixed_accum/this%dt_tran)))
  loc_max_scaled_residual = maxloc(dabs(soln%residual/ &
                                        (fixed_accum/this%dt_tran)),1)

  max_update = maxval(dabs(soln%update*vol_vec/aq_mass_vec))
  loc_max_update = maxloc(dabs(soln%update*vol_vec/aq_mass_vec),1)

  do k = 1,(this%well_grid%nsegments*soln%ndof)
    if (cnvgd_due_to_scaled_res(k) .or. cnvgd_due_to_abs_res(k)) then
      cnvgd_due_to_residual(k) = PETSC_TRUE
    endif
  enddo
  if (all(cnvgd_due_to_abs_res)) then
    rsn_string = trim(rsn_string) // ' aR '
  endif
  if (all(cnvgd_due_to_scaled_res)) then
    rsn_string = trim(rsn_string) // ' sR '
  endif
  if (all(cnvgd_due_to_update)) then
    rsn_string = trim(rsn_string) // ' rU '
  endif

  write(out_string,'(i2," aR:",es10.2,"  sR:",es10.2,"  rU:",es10.2)') &
        n_iter,max_absolute_residual,max_scaled_residual,max_update
  call PrintMsg(this%option,out_string)

  if (all(cnvgd_due_to_residual) .and. all(cnvgd_due_to_update)) then
    soln%converged = PETSC_TRUE
    soln%not_converged = PETSC_FALSE
    out_string = ' TRAN Solution converged!  ---> ' // trim(rsn_string)
    call PrintMsg(this%option,out_string)
    call PrintMsg(this%option,'')
    this%cumulative_dt_tran = this%cumulative_dt_tran + this%dt_tran
    soln%prev_soln%aqueous_conc = this%well%aqueous_conc
    soln%prev_soln%aqueous_mass = this%well%aqueous_mass
  else
    soln%converged = PETSC_FALSE
    soln%not_converged = PETSC_TRUE
  endif

end subroutine PMWellCheckConvergenceTran

! ************************************************************************** !

subroutine PMWellUpdateWellQ(well,reservoir)
  !
  ! Updates the src/sink vector for the fluid object.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

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
    case('CONSTANT_RATE')
      ! weight the rate by the reservoir permeability via the well index
      ! not fully flexible to accommodate single-phase gas
      if (well%bh_ql /= UNINITIALIZED_DOUBLE) then
        liq%Q = well%bh_ql*well%WI/(abs(sum(well%WI)))
        gas%Q = well%bh_qg*well%WI/(abs(sum(well%WI)))
      elseif (well%th_ql /= UNINITIALIZED_DOUBLE) then
        liq%Q = well%th_ql*well%WI/(abs(sum(well%WI)))
        gas%Q = well%th_qg*well%WI/(abs(sum(well%WI)))
      endif
    !------------------------------------------------------------------------
    case('WIPP_DARCY')
      do i = 1,nsegments
        if (dabs((reservoir%p_l(i)-well%pl(i)))/well%pl(i) > threshold_p) then
          upwind = reservoir%p_l(i) > well%pl(i)
          if (upwind) then
            mobility = reservoir%kr_l(i)/reservoir%visc_l(i)
          else
            mobility = liq%kr(i)/liq%visc(i)
          endif
          den_ave = 0.5d0 * (liq%rho(i) + reservoir%rho_l(i)) / FMWH2O
          ! Flowrate in kmol/s
          liq%Q(i) = den_ave*mobility*well%WI(i)* &
                     (reservoir%p_l(i)-well%pl(i))

          upwind = reservoir%p_g(i) > well%pg(i)
          if (upwind) then
            mobility = reservoir%kr_g(i)/reservoir%visc_g(i)
          else
            mobility = gas%kr(i)/gas%visc(i)
          endif
          den_ave = 0.5d0 * (gas%rho(i) + reservoir%rho_g(i)) / &
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

subroutine PMWellComputeWellIndex(this)
  !
  ! Computes the well index.
  !
  ! Author: Michael Nole
  ! Date: 12/22/21

  implicit none

  class(pm_well_type) :: this

  PetscReal :: r0
  PetscReal, parameter :: PI=3.141592653589793d0
  PetscReal :: temp_real
  type(option_type), pointer :: option
  PetscInt :: k
  character(len=8) :: diameter_string

  option => this%option

  ! Peaceman Model: default = anisotropic
  ! This assumes z is vertical (not true for WIPP)
  select case(this%well%WI_model)
    case(PEACEMAN_ISO)

      do k = 1,this%well_grid%nsegments
        write(diameter_string,'(F7.4)') this%well%diameter(k)
        temp_real = log(2.079d-1*this%reservoir%dx(k)/ &
                        (this%well%diameter(k)/2.d0))

        if (temp_real <= 0.d0) then
          option%io_buffer = 'Wellbore diameter (' // diameter_string // '&
          & m) is too large relative to reservoir dx. For the PEACEMAN_ISO &
          &model, wellbore diameter must be smaller than 0.4158 * reservoir dx.'
          call PrintErrMsg(option)
        endif

        this%well%WI(k) = 2.d0*PI*this%reservoir%kx(k)*this%well_grid%dh(k)/ &
                          temp_real
        this%well%r0(k) = 2.079d-1*this%reservoir%dx(k)
       enddo

    case(PEACEMAN_ANISOTROPIC)
      do k = 1,this%well_grid%nsegments
        write(diameter_string,'(F7.4)') this%well%diameter(k)
        r0 = 2.8d-1*(sqrt(sqrt(this%reservoir%ky(k)/this%reservoir%kx(k))* &
             this%reservoir%dx(k)**2 + sqrt(this%reservoir%kx(k)/ &
             this%reservoir%ky(k))*this%reservoir%dy(k)**2) / &
             ((this%reservoir%ky(k)/this%reservoir%kx(k))**2.5d-1 + &
             (this%reservoir%kx(k)/this%reservoir%ky(k))**2.5d-1))

        temp_real = log(r0/(this%well%diameter(k)/2.d0))

        if (temp_real <= 0.d0) then
          option%io_buffer = 'Wellbore diameter (' // diameter_string // ' m)&
          & is too large relative to reservoir discretization and &
          &permeability for the PEACEMAN_ANISOTROPIC well model.'
          call PrintErrMsg(option)
        endif

        this%well%WI(k) = 2.d0*PI*sqrt(this%reservoir%kx(k)* &
                          this%reservoir%ky(k))*this%well_grid%dh(k)/temp_real
        this%well%r0(k) = r0
      enddo
      
  end select

  this%well%WI = this%well%WI*this%well%WI_base

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

  implicit none

  type(pm_well_type) :: pm_well
  type(well_type) :: well
  PetscInt :: id
  PetscReal :: Res(pm_well%nphase)

  Res = 0.d0

  select case(well%well_model_type)
    !---------------------------------------------
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    !---------------------------------------------
    case('WIPP_DARCY')
      ! liquid accumulation term
      Res(1) = Res(1) + well%liq%s(id) * well%liq%rho(id) / FMWH2O * &
                well%phi(id) * well%volume(id) / pm_well%dt_flow
      ! gas accumulation term
      Res(2) = Res(2) + well%gas%s(id) * well%gas%rho(id) / &
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

  call PMWellUpdateWellQ(pm_well%well,pm_well%reservoir)

  select case(well%well_model_type)
    !---------------------------------------------
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    !---------------------------------------------
    case('WIPP_DARCY')
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

subroutine PMWellAccumulationTran(this,isegment,Res)
  !
  ! Computes the fixed accumulation term for the transport residual.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 02/23/2022
  !

  implicit none

  type(pm_well_type) :: this
  PetscInt :: isegment
  PetscReal :: Res(this%nspecies)

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

  do ispecies = 1,this%nspecies
    k = ispecies
    Res(k) = this%well%volume(isegment) * this%well%phi(isegment) * &
             this%well%liq%s(isegment) * &
             this%tran_soln%prev_soln%aqueous_conc(ispecies,isegment)
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

  implicit none

  type(pm_well_type) :: pm_well
  type(well_type) :: well_up, well_dn
  PetscInt :: iup, idn
  PetscReal :: Res(pm_well%nphase)
  PetscBool :: save_flux

  type(well_grid_type), pointer :: well_grid

  PetscInt :: i, ghosted_id
  PetscReal :: pres_up, pres_dn
  PetscReal :: perm_rho_mu_area_ave_over_dist(2), perm_rho_mu_area_up(2), &
               perm_rho_mu_area_dn(2)
  PetscReal :: perm_up, perm_dn, dist_up, dist_dn, density_kg_ave, rel_perm
  PetscReal :: gravity_term, delta_pressure
  PetscReal :: density_ave_kmol, q, tot_mole_flux
  PetscReal :: up_scale, dn_scale
  PetscBool :: upwind


  well_grid => pm_well%well_grid

  Res(:) = 0.d0

  select case(pm_well%well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('WIPP_DARCY')
      ! This is good for either: single-phase liquid, or two-phase liquid/gas.
      ! Vertical well, no Klinkenberg, no capillary pressure, constant mobility.

        perm_up = well_up%permeability(iup)
        perm_dn = well_dn%permeability(idn)
        dist_up = well_grid%dh(iup)/2.d0
        dist_dn = well_grid%dh(idn)/2.d0

        perm_rho_mu_area_up(1) = perm_up * well_up%liq%rho(iup) / &
                              FMWH2O / well_up%liq%visc(iup) * &
                              PI * (well_up%diameter(iup)/2.d0)**2
        perm_rho_mu_area_up(2) = perm_up * well_up%gas%rho(iup) / &
                              fmw_comp(TWO_INTEGER) / well_up%gas%visc(iup) * &
                              PI * (well_up%diameter(iup)/2.d0)**2
        perm_rho_mu_area_dn(1) = perm_dn * well_dn%liq%rho(idn) / &
                              FMWH2O / well_dn%liq%visc(idn) * &
                              PI * (well_dn%diameter(idn)/2.d0)**2
        perm_rho_mu_area_dn(2) = perm_dn * well_dn%gas%rho(idn) / &
                              fmw_comp(TWO_INTEGER) / well_dn%gas%visc(idn) * &
                              PI * (well_dn%diameter(idn)/2.d0)**2

        perm_rho_mu_area_ave_over_dist(1) = &
               (perm_rho_mu_area_up(1) * perm_rho_mu_area_dn(1)) / &
               (dist_up*perm_rho_mu_area_dn(1) + dist_dn*perm_rho_mu_area_up(1))

        perm_rho_mu_area_ave_over_dist(2) = &
               (perm_rho_mu_area_up(2) * perm_rho_mu_area_dn(2)) / &
               (dist_up*perm_rho_mu_area_dn(2) + dist_dn*perm_rho_mu_area_up(2))

        ! Liquid flux
        density_kg_ave = 0.5d0*(well_up%liq%rho(iup)+well_dn%liq%rho(idn))
        ! Assuming the well is always vertical and gravity is in the
        ! (-) direction
        ! And assuming dh is the connection length
        gravity_term = density_kg_ave * gravity * &
                       0.5d0*(well_grid%dh(iup)+well_grid%dh(idn))
        ! No capillary pressure yet.
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
        tot_mole_flux = perm_rho_mu_area_ave_over_dist(1) * rel_perm * &
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
        ! Mobility may need to be non-constant?

        density_kg_ave = 0.5d0*(well_up%gas%rho(iup)+well_dn%gas%rho(idn))
        ! Assuming the well is always vertical and gravity is in the
        ! (-) direction
        gravity_term = density_kg_ave * gravity * &
                       0.5d0*(well_grid%dh(iup)+well_grid%dh(idn))
        ! No capillary pressure yet.
        delta_pressure = well_up%pg(iup) - well_dn%pg(idn) + &
                         gravity_term

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
        tot_mole_flux = perm_rho_mu_area_ave_over_dist(2) * rel_perm * &
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
  PetscReal :: boundary_pressure, boundary_rho
  PetscReal :: boundary_pg, boundary_krg, dn_scale
  PetscReal :: t,dwmol,dwp,dwt,Psat,visl,visg
  PetscReal :: Pc,dpc_dsatl,dkrl_dsatl,dkrg_dsatl
  PetscReal :: v_darcy,q,rel_perm
  PetscBool :: upwind
  PetscInt :: itop
  PetscErrorCode :: ierr

  option => pm_well%option

  well_grid => pm_well%well_grid
  reservoir => pm_well%reservoir

  t = 25.d0 !Constant temperature

  Res(:) = 0.d0

  select case(well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('WIPP_DARCY')
      ! This is good for either: single-phase liquid, or two-phase liquid/gas.
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
        gravity_term = well%liq%rho(1) * gravity * &
                       well_grid%dh(1)/2.d0
        delta_pressure = boundary_pressure - well%pl(1) + gravity_term

        call EOSWaterSaturationPressure(t,Psat,ierr)
        call EOSWaterDensityBRAGFLO(t,boundary_pressure,PETSC_FALSE, &
                                boundary_rho,dwmol,dwp,dwt,ierr)

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
          density_ave = boundary_rho
        else
          density_ave = (well%liq%rho(1)+boundary_rho) / &
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

        gravity_term = well%gas%rho(1) * gravity * &
                       well_grid%dh(1)/2.d0
        delta_pressure = boundary_pg - well%pg(1) + gravity_term

        call EOSGasDensity(t,boundary_pg,boundary_rho,ierr)

        upwind = delta_pressure > 0.d0
        if (upwind) then
          call characteristic_curves%gas_rel_perm_function% &
               RelativePermeability(1.d0-well%bh_sg,rel_perm,dkrl_dsatl,option)
          call EOSGasViscosity(t,boundary_pg,boundary_pg,boundary_rho,visg,ierr)
        else
          dn_scale = 1.d0
          rel_perm = well%gas%kr(1)
          visg = well%gas%visc(1)
        endif

        v_darcy = perm_ave_over_dist * rel_perm/visg * &
                  delta_pressure
        if (upwind) then
          density_ave = boundary_rho
        else
          density_ave = (well%gas%rho(1)+boundary_rho) / &
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
          density_ave = reservoir%rho_l(1) / fmw_comp(ONE_INTEGER)
        else
          density_ave = well%liq%rho(1) / fmw_comp(ONE_INTEGER)
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
          density_ave = reservoir%rho_g(1) / fmw_comp(TWO_INTEGER)
        else
          density_ave = well%gas%rho(1) / fmw_comp(TWO_INTEGER)
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
        gravity_term = well%liq%rho(itop) * gravity * &
                       well_grid%dh(itop)/2.d0
        delta_pressure = well%pl(itop) - boundary_pressure + gravity_term

        call EOSWaterSaturationPressure(t,Psat,ierr)
        call EOSWaterDensityBRAGFLO(t,boundary_pressure,PETSC_FALSE, &
                                boundary_rho,dwmol,dwp,dwt,ierr)

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

        density_ave = (well%liq%rho(itop)+boundary_rho) / &
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
              option%flow%pct_updated = PETSC_TRUE
              call sat_func% &
                   CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
            else
              call sat_func% &
                   CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
            endif
          class default
            call sat_func% &
               CapillaryPressure(1.d0-well%th_sg,Pc,dpc_dsatl,option)
        end select
        boundary_pg = boundary_pressure + Pc

        call characteristic_curves%gas_rel_perm_function% &
               RelativePermeability(1.d0-well%th_sg,boundary_krg, &
               dkrg_dsatl,option)

        gravity_term = well%gas%rho(itop) * gravity * &
                       well_grid%dh(itop)/2.d0
        delta_pressure = well%pg(itop) - boundary_pressure + gravity_term

        call EOSGasDensity(t,boundary_pg,boundary_rho,ierr)

        upwind = delta_pressure < 0.d0
        if (upwind) then
          call characteristic_curves%gas_rel_perm_function% &
               RelativePermeability(1.d0-well%th_sg,rel_perm,dkrl_dsatl,option)
          call EOSGasViscosity(t,boundary_pg,boundary_pg,boundary_rho,visg,ierr)
        else
          dn_scale = 1.d0
          rel_perm = well%gas%kr(itop)
          visg = well%gas%visc(itop)
        endif

        v_darcy = perm_ave_over_dist * rel_perm/visg * &
                  delta_pressure

        density_ave = (well%gas%rho(itop)+boundary_rho) / (2.d0 *fmw_comp(TWO_INTEGER))
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
        density_ave = well%liq%rho(itop) / fmw_comp(ONE_INTEGER)
        q = v_darcy * well%area(itop)
        tot_mole_flux = q * density_ave
        ! Store boundary flux for consistency with transport
        if (save_flux) then
          well%ql_bc(2) = v_darcy
          well%ql_kmol_bc(2) = tot_mole_flux
        endif
        Res(3) = Res(3) - tot_mole_flux

        v_darcy = well%th_qg

        ! Always take well density with tophole flux bc
        density_ave = well%gas%rho(itop) / fmw_comp(TWO_INTEGER)
        q = v_darcy * well%area(itop)
        tot_mole_flux = q * density_ave
        ! Store boundary flux for consistency with transport
        if (save_flux) then
          well%qg_bc(2) = -v_darcy
          well%qg_kmol_bc(2) = -tot_mole_flux
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

  ! I don't like how pressure and saturation are attributes of
  ! different objects
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
  call PMWellUpdatePropertiesFlow(pm_well,pm_well%well_pert(ONE_INTEGER), &
                        pm_well%realization%patch%characteristic_curves_array, &
                        pm_well%realization%option)
  call PMWellUpdatePropertiesFlow(pm_well,pm_well%well_pert(TWO_INTEGER), &
                        pm_well%realization%patch%characteristic_curves_array, &
                        pm_well%realization%option)

  ! Update perturbed source/sink term from the reservoir
  call PMWellUpdateWellQ(pm_well%well_pert(ONE_INTEGER),pm_well%reservoir)
  call PMWellUpdateWellQ(pm_well%well_pert(TWO_INTEGER),pm_well%reservoir)

  pm_well%pert = pert

end subroutine PMWellPerturb

! ************************************************************************** !

subroutine PMWellUpdatePropertiesFlow(this,well,characteristic_curves_array, &
                                      option)
  !
  ! Updates flow well object properties.
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

  class(pm_well_type) :: this
  type(well_type) :: well
  type(characteristic_curves_ptr_type), pointer::characteristic_curves_array(:)
  type(option_type) :: option

  class(characteristic_curves_type), pointer :: characteristic_curves
  class(sat_func_base_type), pointer :: saturation_function
  type(strata_type), pointer :: strata
  PetscInt :: i,nsegments
  PetscReal :: t,dw,dg,dwmol,dwp,dwt,Psat,visl,visg
  PetscReal :: Pc,dpc_dsatl,krl,dkrl_dsatl,krg,dkrg_dsatl
  PetscErrorCode :: ierr

  nsegments =this%well_grid%nsegments

  do i = 1,nsegments
    ! Material Properties
    strata => this%strata_list%first
    do
      if (.not.associated(strata)) exit
      if (strata%id == this%well_grid%strata_id(i)) then
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
          option%flow%pct_updated = PETSC_TRUE
          call sat_func% &
                   CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
        else
          call sat_func% &
                   CapillaryPressure(well%liq%s(i),Pc,dpc_dsatl,option)
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
    call EOSWaterDensityBRAGFLO(t,well%pl(i),PETSC_FALSE, &
                                dw,dwmol,dwp,dwt,ierr)
    call EOSGasDensity(t,well%pg(i),dg,ierr)

    well%liq%rho(i) = dw
    !No water vapor in WIPP_Darcy mode
    well%gas%rho(i) = dg

    !Viscosity
    call EOSWaterSaturationPressure(t,Psat,ierr)
    call EOSWaterViscosity(t,well%pl(i),Psat,visl,ierr)
    call EOSGasViscosity(t,well%pg(i),well%pg(i),dg,visg,ierr)

    well%liq%visc(i) = visl
    well%gas%visc(i) = visg

  enddo


end subroutine PMWellUpdatePropertiesFlow

! ************************************************************************** !

subroutine PMWellUpdatePropertiesTran(this)
  !
  ! Updates transport related well object properties for each time step.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/13/2022

  use Transport_Constraint_Base_module
  use Transport_Constraint_NWT_module
  use Transport_Constraint_module
  use Condition_module

  implicit none

  class(pm_well_type) :: this

  type(tran_condition_type), pointer :: tran_condition
  class(tran_constraint_base_type), pointer :: cur_constraint

  ! update the top of hole boundary condition with current constraint
  tran_condition => this%realization%transport_conditions%first
  do
    if (.not.associated(tran_condition)) exit
      if (trim(tran_condition%name) == &
          trim(this%well%tran_condition_name)) exit
    tran_condition => tran_condition%next
  enddo
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

end subroutine PMWellUpdatePropertiesTran

! ************************************************************************** !
subroutine PMWellCopyWell(well,well_copy)
  !
  ! Copies well object properties from one to another.
  !
  ! Author: Michael Nole
  ! Date: 02/25/2022
  !

  implicit none

  type(well_type) :: well
  type(well_type) :: well_copy

  well_copy%liq%rho(:) = well%liq%rho0
  well_copy%liq%visc(:) = well%liq%visc(:)
  well_copy%liq%kr = well%liq%kr
  well_copy%gas%rho(:) = well%gas%rho0
  well_copy%gas%visc(:) = well%gas%visc(:)
  well_copy%gas%kr = well%gas%kr
  well_copy%diameter(:) = well%diameter(:)
  well_copy%WI_base(:) = well%WI_base(:)
  well_copy%permeability(:) = well%permeability(:)
  well_copy%phi(:) = well%phi(:)
  well_copy%f(:) = well%f(:)
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


end subroutine PMWellCopyWell

! ************************************************************************** !

subroutine PMWellSetPlotVariables(list,this)
  !
  ! Adds variables to be printed for plotting.
  !
  ! Author: Jenn Frederick
  ! Date: 09/29/2022
  !
  use Output_Aux_module
  use Variables_module

  type(output_variable_list_type), pointer :: list
  class(pm_well_type) :: this

  character(len=MAXWORDLENGTH) :: name,  units
  PetscInt :: i

  if (this%well%output_pl) then
    name = 'Well Liq. Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                 WELL_LIQ_PRESSURE) 
  endif
  
  if (this%well%output_pg) then
    name = 'Well Gas Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                 WELL_GAS_PRESSURE)
  endif

  if (this%well%output_aqc .and. this%transport) then
    do i=1,this%nspecies
      name = 'Well AQ Conc. ' // trim(this%well%species_names(i))
      units = 'mol/m^3-liq'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   WELL_AQ_CONC,i)
    enddo
  endif

  if (this%well%output_aqm .and. this%transport) then
    do i=1,this%nspecies
      name = 'Well AQ Mass ' // trim(this%well%species_names(i))
      units = 'mol'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   WELL_AQ_MASS,i)
    enddo  
  endif

  if (this%well%liq%output_Q) then
    name = 'Well Liq. Q'
    units = 'kmol/sec'
    call OutputVariableAddToList(list,name,OUTPUT_RATE,units,WELL_LIQ_Q)
  endif

  if (this%well%gas%output_Q) then
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

  implicit none

  type(option_type), pointer :: option

  character(len=MAXSTRINGLENGTH) :: PMWellOutputFilename

  PMWellOutputFilename = trim(option%global_prefix) // &
                         trim(option%group_prefix) // '.well'

end function PMWellOutputFilename

! ************************************************************************** !

subroutine PMWellOutputHeader(this)
  !
  ! Writes the header for wellbore model output file.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021

  use Output_Aux_module
  use Grid_module
  use Utility_module

  implicit none

  class(pm_well_type) :: this

  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: cell_string
  PetscBool :: exist
  PetscInt :: fid
  PetscInt :: icolumn
  PetscInt :: k, j

  output_option => this%realization%output_option

  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif

  fid = 555
  filename = PMWellOutputFilename(this%option)
  exist = FileExists(trim(filename))
  if (this%option%restart_flag .and. exist) return
  open(unit=fid,file=filename,action="write",status="replace")

  write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'
  cell_string = ''

  do k = 1,this%well_grid%nsegments
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
    if (this%transport) then
      do j = 1,this%nspecies
        variable_string = 'Well Aqueous Conc. ' // &
                          trim(this%well%species_names(j))
        units_string = 'mol/m^3-liq'
        call OutputWriteToHeader(fid,variable_string,units_string, &
                                 cell_string,icolumn)
        variable_string = 'Well Aqueous Mass. ' &
                           // trim(this%well%species_names(j))
        units_string = 'mol'
        call OutputWriteToHeader(fid,variable_string,units_string, &
                                 cell_string,icolumn)
        variable_string = 'Res Aqueous Conc. ' // &
                          trim(this%well%species_names(j))
        units_string = 'mol/m^3-liq'
        call OutputWriteToHeader(fid,variable_string,units_string, &
                                 cell_string,icolumn)
      enddo
    endif
  enddo

  close(fid)

end subroutine PMWellOutputHeader

! ************************************************************************** !

subroutine PMWellOutput(this)
  !
  ! Sets up output for the wellbore process model
  !
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021

  use Output_Aux_module
  use Global_Aux_module
  use Grid_module

  implicit none

  class(pm_well_type) :: this

  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: k, j

100 format(100es18.8)
101 format(1I6.1)

  option => this%realization%option
  output_option => this%realization%output_option

  fid = 555
  filename = PMWellOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")

  ! this time is set at the end of the wellbore step????
  write(fid,100,advance="no") option%time / output_option%tconv

  do k = 1,this%well_grid%nsegments
    write(fid,101,advance="no") k
    write(fid,100,advance="no") this%well_grid%h(k)%x, &
                                this%well_grid%h(k)%y, &
                                this%well_grid%h(k)%z, &
                                this%well%pl(k), &
                                this%reservoir%p_l(k), &
                                this%well%liq%s(k), &
                                this%well%gas%s(k), &
                                this%well%liq%Q(k), &
                                this%well%gas%Q(k)
    if (k == 1) then
      write(fid,100,advance="no") this%well%ql_bc(1), &
                                  this%well%qg_bc(1)
    endif
    if (k > 1) then
      write(fid,100,advance="no") this%well%ql(k-1), &
                                  this%well%qg(k-1)
    endif
    write(fid,100,advance="no") this%well%mass_balance_liq(k), &
                                this%well%liq_mass(k)
    if (this%transport) then
      do j = 1,this%nspecies
        write(fid,100,advance="no") this%well%aqueous_conc(j,k), &
                                    this%well%aqueous_mass(j,k), &
                                    this%reservoir%aqueous_conc(j,k)
      enddo
    endif
  enddo

  close(fid)

end subroutine PMWellOutput

! ************************************************************************** !

subroutine PMWellMassBalance(this)
  !
  ! Verify the mass balance in the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 06/28/2022

  implicit none

  class(pm_well_type) :: this

  type(well_type), pointer :: well
  PetscInt :: isegment, nsegments
  PetscInt :: n_up, n_dn
  PetscReal :: area_up, area_dn 
  PetscReal :: rho_kg_up, rho_kg_dn
  PetscReal :: mass_rate_up, mass_rate_dn
  PetscReal :: q_up, q_dn

  ! q in [m3-liq/m2-bulk-sec]
  ! area in [m2-bulk]

  ! (-) mass balance rate means net mass is being lost
  ! (+) mass balance rate means net mass is being gained

  well => this%well
  nsegments = this%well_grid%nsegments

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

    else if (isegment == this%well_grid%nsegments) then
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

subroutine PMWellUpdateMass(this)
  !
  ! Verify the mass balance in the well process model.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/17/2022

  implicit none

  class(pm_well_type) :: this

  type(well_type), pointer :: well
  PetscInt :: isegment, nsegments
  PetscReal :: inst_mass

  well => this%well
  nsegments = this%well_grid%nsegments

  do isegment = 1,nsegments

    ! [kmol]
    inst_mass = well%volume(isegment)*well%phi(isegment)* &
                well%liq%s(isegment)*well%liq%rho(isegment)*FMWH2O

    well%liq_mass(isegment) = inst_mass
    well%liq_cum_mass(isegment) = well%liq_cum_mass(isegment) + inst_mass

  enddo

end subroutine PMWellUpdateMass

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

  call PMBaseDestroy(this)

  call DeallocateArray(this%well_grid%h_local_id)
  call DeallocateArray(this%well_grid%h_ghosted_id)
  call DeallocateArray(this%well_grid%dh)
  nullify(this%well_grid)

  call DeallocateArray(this%reservoir%p_l)
  call DeallocateArray(this%reservoir%p_g)
  call DeallocateArray(this%reservoir%s_l)
  call DeallocateArray(this%reservoir%s_g)
  call DeallocateArray(this%reservoir%mobility_l)
  call DeallocateArray(this%reservoir%mobility_g)
  call DeallocateArray(this%reservoir%kr_l)
  call DeallocateArray(this%reservoir%kr_g)
  call DeallocateArray(this%reservoir%rho_l)
  call DeallocateArray(this%reservoir%rho_g)
  call DeallocateArray(this%reservoir%visc_l)
  call DeallocateArray(this%reservoir%visc_g)
  call DeallocateArray(this%reservoir%e_por)
  if (this%transport) then
    call DeallocateArray(this%reservoir%aqueous_conc)
    call DeallocateArray(this%reservoir%aqueous_mass)
  endif
  nullify(this%reservoir)


  call DeallocateArray(this%well%area)
  call DeallocateArray(this%well%diameter)
  call DeallocateArray(this%well%volume)
  call DeallocateArray(this%well%f)
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
  call DeallocateArray(this%well%liq%rho)
  call DeallocateArray(this%well%liq%visc)
  call DeallocateArray(this%well%liq%s)
  call DeallocateArray(this%well%liq%Q)
  call DeallocateArray(this%well%gas%rho)
  call DeallocateArray(this%well%gas%visc)
  call DeallocateArray(this%well%gas%s)
  call DeallocateArray(this%well%gas%Q)
  nullify(this%well%liq)
  nullify(this%well%gas)
  nullify(this%well)

  call DeallocateArray(this%well_pert(ONE_INTEGER)%area)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%diameter)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%volume)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%f)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%WI)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%pl)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%pg)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%liq%rho)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%liq%visc)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%liq%s)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%liq%Q)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%gas%rho)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%gas%visc)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%gas%s)
  call DeallocateArray(this%well_pert(ONE_INTEGER)%gas%Q)
  nullify(this%well_pert(ONE_INTEGER)%liq)
  nullify(this%well_pert(ONE_INTEGER)%gas)

  call DeallocateArray(this%well_pert(TWO_INTEGER)%area)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%diameter)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%volume)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%f)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%WI)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%pl)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%pg)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%liq%rho)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%liq%visc)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%liq%s)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%liq%Q)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%gas%rho)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%gas%visc)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%gas%s)
  call DeallocateArray(this%well_pert(TWO_INTEGER)%gas%Q)
  nullify(this%well_pert(TWO_INTEGER)%liq)
  nullify(this%well_pert(TWO_INTEGER)%gas)
  nullify(this%well_pert)

  call DeallocateArray(this%flow_soln%residual)
  call DeallocateArray(this%flow_soln%Jacobian)
  call DeallocateArray(this%flow_soln%update)
  call DeallocateArray(this%flow_soln%prev_soln%pl)
  call DeallocateArray(this%flow_soln%prev_soln%sg)
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
