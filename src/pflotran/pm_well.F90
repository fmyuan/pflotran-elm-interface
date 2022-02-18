module PM_Well_class

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscsnes.h"
  use petscsys
  use petscsnes
  use PM_Base_class
  use Option_module
  use Geometry_module
  use Realization_Subsurface_class
  
  use PFLOTRAN_Constants_module

  use WIPP_Flow_Aux_module

  implicit none

  private

  PetscBool :: initialize_well = PETSC_TRUE
  PetscReal, parameter :: gravity = -9.80665d0 ! [m/s2]

  type :: well_grid_type
    ! number of well segments
    PetscInt :: nsegments      
    ! number of well connections
    PetscInt :: nconnections
    ! delta h discretization of each segment center [m]       
    PetscReal, pointer :: dh(:)
    ! h coordinate of each segment center [m]
    type(point3d_type), pointer :: h(:)  
    ! the local id of the reservoir grid cell within which each segment 
    ! center resides 
    PetscInt, pointer :: h_local_id(:) 
    ! the ghosted id of the reservoir grid cell within which each segment 
    ! center resides 
    PetscInt, pointer :: h_ghosted_id(:)
    ! coordinate of the top/bottom of the well [m]
    PetscReal :: tophole(3)
    PetscReal :: bottomhole(3)     
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
  end type

  type :: well_type
    ! type of well model to use
    character(len=MAXWORDLENGTH) :: well_model_type
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
    ! well index model (probably has to get moved out of well_type)
    character(len=MAXWORDLENGTH) :: WI_model
    ! well liquid pressure [Pa]     
    PetscReal, pointer :: pl(:) 
    ! well gas pressure [Pa]
    PetscReal, pointer :: pg(:)
    ! well liquid Darcy flux [m3/m2-s] of interior interfaces  
    PetscReal, pointer :: ql(:)
    ! well gas Darcy flux [m3/m2-s] of interior interfaces  
    PetscReal, pointer :: qg(:)
    ! well species aqueous concentration [mol/m3-liq] (idof,conc@segment)
    PetscReal, pointer :: aqueous_conc(:,:)
    ! well species aqueous mass [mol] (idof,mass@segment)
    PetscReal, pointer :: aqueous_mass(:,:)
    ! well bottom of hole pressure BC flag
    PetscBool :: bh_p_set_by_reservoir
    ! well bottom of hole pressure BC [Pa]
    PetscReal :: bh_p 
    ! well top of hole pressure BC [Pa]
    PetscReal :: th_p
    ! well bottom of hole rate BC [kg/s]
    PetscReal :: bh_ql
    PetscReal :: bh_qg
    ! well top of hole rate BC [kg/s]
    PetscReal :: th_ql
    PetscReal :: th_qg
    ! permeability along the well [m2]
    PetscReal, pointer :: permeability(:)
    ! porosity
    PetscReal, pointer :: phi(:)
  end type well_type

  type :: well_fluid_type
    ! fluid phase ID (liq/gas)
    PetscInt :: ifluid
    ! fluid mobility (phase_permeability/viscosity)
    PetscReal :: mobility
    ! fluid reference density [kg/m3]
    PetscReal :: rho0
    ! fluid density [kg/m3]
    PetscReal, pointer :: rho(:)
    ! fluid viscosity [Pa-s]
    PetscReal, pointer :: visc(:)
    ! fluid saturation
    PetscReal, pointer :: s(:)
    ! fluid source/sink in/out of well [kg/s]   
    PetscReal, pointer :: Q(:)
    ! equation of state for density
    !procedure(rho_interface), pointer, nopass :: update_rho_ptr => null()
  end type well_fluid_type

  !interface
  !  subroutine rho_interface(fluid)
  !    type(well_fluid_type) :: fluid
  !  end subroutine rho_interface
  !end interface

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
    ! time step cuts factor (divides the time step)
    PetscInt :: ts_cut_factor
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
    type(well_grid_type), pointer :: grid
    type(well_type), pointer :: well
    type(well_type), pointer :: well_pert(:)
    type(well_reservoir_type), pointer :: reservoir 
    type(well_soln_flow_type), pointer :: flow_soln
    type(well_soln_tran_type), pointer :: tran_soln
    PetscReal, pointer :: pert(:,:)
    PetscInt :: nphase 
    PetscInt :: nspecies 
    PetscReal :: dt
    PetscReal :: min_dt
    PetscReal :: cumulative_dt
    PetscBool :: transport
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
  PMWellCreate%min_dt = 1.d-15
  PMWellCreate%nphase = 0
  PMWellCreate%nspecies = 0
  PMWellCreate%transport = PETSC_FALSE

  nullify(PMWellCreate%pert)

  ! create the well grid object:
  allocate(PMWellCreate%grid)
  PMWellCreate%grid%nsegments = UNINITIALIZED_INTEGER
  PMWellCreate%grid%nconnections = UNINITIALIZED_INTEGER
  nullify(PMWellCreate%grid%dh)
  nullify(PMWellCreate%grid%h)
  nullify(PMWellCreate%grid%h_local_id)
  nullify(PMWellCreate%grid%h_ghosted_id)
  PMWellCreate%grid%tophole(:) = UNINITIALIZED_DOUBLE
  PMWellCreate%grid%bottomhole(:) = UNINITIALIZED_DOUBLE

  ! create the well object:
  allocate(PMWellCreate%well)
  nullify(PMWellCreate%well%area)
  nullify(PMWellCreate%well%diameter)
  nullify(PMWellCreate%well%volume)
  nullify(PMWellCreate%well%f)
  nullify(PMWellCreate%well%WI_base)
  nullify(PMWellCreate%well%WI)
  nullify(PMWellCreate%well%pl)
  nullify(PMWellCreate%well%pg)
  nullify(PMWellCreate%well%ql)
  nullify(PMWellCreate%well%qg)
  nullify(PMWellCreate%well%aqueous_conc)
  nullify(PMWellCreate%well%aqueous_mass)
  nullify(PMWellCreate%well%permeability)
  nullify(PMWellCreate%well%phi)
  PMWellCreate%well%bh_p_set_by_reservoir = PETSC_FALSE
  PMWellCreate%well%bh_p = UNINITIALIZED_DOUBLE
  PMWellCreate%well%th_p = UNINITIALIZED_DOUBLE
  PMWellCreate%well%bh_ql = UNINITIALIZED_DOUBLE
  PMWellCreate%well%bh_qg = UNINITIALIZED_DOUBLE
  PMWellCreate%well%th_ql = UNINITIALIZED_DOUBLE
  PMWellCreate%well%th_qg = UNINITIALIZED_DOUBLE

  ! create the well_pert object:
  allocate(PMWellCreate%well_pert(TWO_INTEGER))
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%area)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%diameter)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%volume)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%f)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%WI_base)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%WI)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%pl)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%pg)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%ql)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%qg)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%aqueous_conc)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%aqueous_mass)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%permeability)
  nullify(PMWellCreate%well_pert(ONE_INTEGER)%phi)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%area)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%diameter)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%volume)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%f)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%WI_base)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%WI)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%pl)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%pg)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%ql)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%qg)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%aqueous_conc)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%aqueous_mass)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%permeability)
  nullify(PMWellCreate%well_pert(TWO_INTEGER)%phi)
  PMWellCreate%well_pert(:)%bh_p_set_by_reservoir = PETSC_FALSE
  PMWellCreate%well_pert(:)%bh_p = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(:)%th_p = UNINITIALIZED_DOUBLE
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


  ! create the fluid/liq objects:
  allocate(PMWellCreate%well%liq)
  PMWellCreate%well%liq%ifluid = 1 
  PMWellCreate%well%liq%mobility = UNINITIALIZED_DOUBLE
  PMWellCreate%well%liq%rho0 = UNINITIALIZED_DOUBLE
  nullify(PMWellCreate%well%liq%rho)
  nullify(PMWellCreate%well%liq%visc)
  nullify(PMWellCreate%well%liq%s)
  nullify(PMWellCreate%well%liq%Q)
  allocate(PMWellCreate%well_pert(ONE_INTEGER)%liq)
  allocate(PMWellCreate%well_pert(TWO_INTEGER)%liq)
  PMWellCreate%well_pert(ONE_INTEGER)%liq%ifluid = 1
  PMWellCreate%well_pert(TWO_INTEGER)%liq%ifluid = 1
  PMWellCreate%well_pert(ONE_INTEGER)%liq%mobility = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(TWO_INTEGER)%liq%mobility = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(ONE_INTEGER)%liq%rho0 = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(TWO_INTEGER)%liq%rho0 = UNINITIALIZED_DOUBLE

  ! create the fluid/gas objects:
  allocate(PMWellCreate%well%gas)
  PMWellCreate%well%gas%ifluid = 2
  PMWellCreate%well%gas%mobility = UNINITIALIZED_DOUBLE
  PMWellCreate%well%gas%rho0 = UNINITIALIZED_DOUBLE
  nullify(PMWellCreate%well%gas%rho)
  nullify(PMWellCreate%well%gas%visc)
  nullify(PMWellCreate%well%gas%s)
  nullify(PMWellCreate%well%gas%Q)
  allocate(PMWellCreate%well_pert(ONE_INTEGER)%gas)
  allocate(PMWellCreate%well_pert(TWO_INTEGER)%gas)
  PMWellCreate%well_pert(ONE_INTEGER)%gas%ifluid = 2
  PMWellCreate%well_pert(TWO_INTEGER)%gas%ifluid = 2
  PMWellCreate%well_pert(ONE_INTEGER)%gas%mobility = UNINITIALIZED_DOUBLE
  PMWellCreate%well_pert(TWO_INTEGER)%gas%mobility = UNINITIALIZED_DOUBLE
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
  PMWellCreate%flow_soln%ts_cut_factor = 2
  PMWellCreate%flow_soln%itol_abs_res = 1.0d-8
  PMWellCreate%flow_soln%itol_scaled_res = 1.0d-4
  PMWellCreate%flow_soln%itol_abs_update_p = 1.0d0
  PMWellCreate%flow_soln%itol_abs_update_s = 1.0d-3
  PMWellCreate%flow_soln%itol_rel_update_p = 1.0d-1
  PMWellCreate%flow_soln%itol_rel_update_s = 1.0d-1
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
  PMWellCreate%tran_soln%ts_cut_factor = 2
  PMWellCreate%tran_soln%itol_abs_res = 1.0d-8
  PMWellCreate%tran_soln%itol_scaled_res = 1.0d-4
  PMWellCreate%tran_soln%itol_abs_update = 1.0d0
  PMWellCreate%tran_soln%itol_rel_update = 1.0d-1
  PMWellCreate%tran_soln%n_steps = 0
  PMWellCreate%tran_soln%n_newton = 0


end function PMWellCreate

! ************************************************************************** !
  
subroutine PMWellSetup(this)
  ! 
  ! Initializes variables associated with the well process model.
  ! 
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  ! 

  use Option_module
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
  type(coupler_type), pointer :: source_sink
  type(input_type) :: input_dummy
  class(dataset_ascii_type), pointer :: dataset_ascii
  character(len=MAXSTRINGLENGTH) :: string
  type(tran_constraint_coupler_nwt_type), pointer :: tran_constraint_coupler_nwt
  PetscReal :: diff_x,diff_y,diff_z
  PetscReal :: dh_x,dh_y,dh_z
  PetscReal :: total_length
  PetscReal :: temp_real
  PetscInt :: local_id
  PetscInt :: nsegments
  PetscInt :: k

  option => this%option
  res_grid => this%realization%patch%grid
  nsegments = this%grid%nsegments

  allocate(this%grid%dh(nsegments))
  allocate(this%grid%h(nsegments))
  allocate(this%grid%h_local_id(nsegments))
  allocate(this%grid%h_ghosted_id(nsegments))

  this%grid%nconnections = this%grid%nsegments - 1

  diff_x = this%grid%tophole(1)-this%grid%bottomhole(1)
  diff_y = this%grid%tophole(2)-this%grid%bottomhole(2)
  diff_z = this%grid%tophole(3)-this%grid%bottomhole(3)

  dh_x = diff_x/nsegments
  dh_y = diff_y/nsegments
  dh_z = diff_z/nsegments

  diff_x = diff_x*diff_x
  diff_y = diff_y*diff_y
  diff_z = diff_z*diff_z

  total_length = sqrt(diff_x+diff_y+diff_z)

  do k = 1,this%grid%nsegments
    this%grid%h(k)%id = k
    this%grid%h(k)%x = this%grid%bottomhole(1)+(dh_x*(k-0.5))
    this%grid%h(k)%y = this%grid%bottomhole(2)+(dh_y*(k-0.5))
    this%grid%h(k)%z = this%grid%bottomhole(3)+(dh_z*(k-0.5))
    !if (k < this%grid%nsegments) then
    !  this%grid%h(k)%id = k
    !  this%grid%h(k)%x = this%grid%h(k)%x + dh_x/2.d0
    !  this%grid%h(k)%y = this%grid%h(k)%y + dh_y/2.d0
    !  this%grid%h(k)%z = this%grid%h(k)%z + dh_z/2.d0
    !endif
  enddo

  this%grid%dh(:) = total_length/nsegments

  ! Get the local_id for each well segment center from the reservoir grid
  this%grid%h_local_id(:) = -1  
  do k = 1,this%grid%nsegments
    call GridGetLocalIDFromCoordinate(res_grid,this%grid%h(k), &
                                      option,local_id)
    this%grid%h_local_id(k) = local_id
    this%grid%h_ghosted_id(k) = res_grid%nL2G(local_id)
  enddo

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
  if (associated(this%well%permeability) .and. size(this%well%permeability) &
      /= nsegments) then
    if (size(this%well%permeability) == 1) then
      temp_real = this%well%permeability(1)
      deallocate(this%well%permeability)
      allocate(this%well%permeability(nsegments))
      this%well%permeability(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,WELL,PERMEABILITY must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,WELL,PERMEABILITY, &
        &it will be applied to all segments of the well.'
      call PrintErrMsg(option)
    endif
  else
    allocate(this%well%permeability(nsegments))
    this%well%permeability(:) = UNINITIALIZED_DOUBLE
  endif
  if (associated(this%well%phi) .and. size(this%well%phi) &
      /= nsegments) then
    if (size(this%well%phi) == 1) then
      temp_real = this%well%phi(1)
      deallocate(this%well%phi)
      allocate(this%well%phi(nsegments))
      this%well%phi(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,WELL,POROSITY must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,WELL,PERMEABILITY, &
        &it will be applied to all segments of the well.'
      call PrintErrMsg(option)
    endif
  else
    allocate(this%well%phi(nsegments))
    this%well%phi(:) = UNINITIALIZED_DOUBLE
  endif
  if (size(this%well%f) /= nsegments) then
    if (size(this%well%f) == 1) then
      temp_real = this%well%f(1)
      deallocate(this%well%f)
      allocate(this%well%f(nsegments))
      this%well%f(:) = temp_real
    else
      option%io_buffer = 'The number of values provided in &
        &WELLBORE_MODEL,WELL,FRICTION_COEFFICIENT must match the number &
        &of well segments provided in &
        &WELLBORE_MODEL,GRID,NUMBER_OF_SEGMENTS. Alternatively, if &
        &a single value is provided in WELLBORE_MODEL,WELL,FRICTION_COEFFICIENT, &
        &it will be applied to all segments of the well.'
      call PrintErrMsg(option)
    endif
  endif

  allocate(this%well%area(nsegments))
  this%well%area = 3.14159*(this%well%diameter/2.d0)*(this%well%diameter/2.d0)

  allocate(this%well%volume(nsegments))
  this%well%volume = this%well%area*this%grid%dh

  allocate(this%well%liq%visc(nsegments))
  allocate(this%well%gas%visc(nsegments))

  this%flow_soln%ndof = this%nphase

  if (this%option%itranmode /= NULL_MODE) then
    this%transport = PETSC_TRUE
    this%nspecies = this%realization%reaction_nw%params%nspecies
    this%tran_soln%ndof = this%nspecies 
  endif

  ! add a reservoir src/sink coupler for each well segment
  do k = 1,this%grid%nsegments
    write(string,'(I0.4)') k
    source_sink => CouplerCreate(SRC_SINK_COUPLER_TYPE)
    source_sink%name = 'well_segment_' // trim(string)

    ! ----- flow ----------------
    source_sink%flow_condition_name = 'well_segment_' // trim(string) // &
                                      '_flow_srcsink'
    source_sink%flow_condition => FlowConditionCreate(option)
    source_sink%flow_condition%name = source_sink%flow_condition_name
    source_sink%flow_condition%general => FlowGeneralConditionCreate(option)
    string = 'RATE'
    source_sink%flow_condition%general%rate => FlowGeneralSubConditionPtr( &
      input_dummy,string,source_sink%flow_condition%general,option)
    source_sink%flow_condition%general%rate%itype = MASS_RATE_SS ! [kg/s]
    allocate(source_sink%flow_condition%general%rate%dataset%rarray(2))
    source_sink%flow_condition%general%rate%dataset%rarray(:) = 0.d0 

    ! ----- transport -----------
    if (this%transport) then
      write(string,'(I0.4)') k
      source_sink%tran_condition_name = 'well_segment_' // trim(string) // &
                                        '_tran_srcsink'
      source_sink%tran_condition => TranConditionCreate(option) 
      source_sink%tran_condition%name = source_sink%tran_condition_name
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
    source_sink%connection_set%id_dn = this%grid%h_local_id(k)

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
    !-------------------------------------
      case('SINGLE_PHASE')
        this%nphase = 1
        cycle
    !-------------------------------------
      case('TWO_PHASE')
        this%nphase = 2
        cycle
    !-------------------------------------
      case('LIQUID_MOBILITY')
        call InputReadDouble(input,option,this%well%liq%mobility)
        call InputErrorMsg(input,option,'LIQUID_MOBILITY value', &
                           error_string)
        cycle
    !-------------------------------------
      case('GAS_MOBILITY')
        call InputReadDouble(input,option,this%well%gas%mobility)
        call InputErrorMsg(input,option,'GAS_MOBILITY value', &
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

  if (Uninitialized(this%well%liq%mobility)) then
    option%io_buffer = 'LIQUID_MOBILITY must be provided in the &
                       &WELLBORE_MODEL block.'
    call PrintErrMsg(option)
  endif

  if ((this%nphase == 2) .and. (Uninitialized(this%well%gas%mobility))) then
    option%io_buffer = 'GAS_MOBILITY must be provided in the &
                       &WELLBORE_MODEL block.'
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
  use Option_module
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

  error_string = trim(error_string) // ',GRID'
  found = PETSC_TRUE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('GRID')
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
            call InputReadInt(input,option,pm_well%grid%nsegments)
            call InputErrorMsg(input,option,'NUMBER_OF_SEGMENTS',error_string)
        !-----------------------------
          case('TOP_OF_HOLE')
            call InputReadDouble(input,option,pm_well%grid%tophole(1))
            call InputErrorMsg(input,option,'TOP_OF_HOLE x-coordinate', &
                               error_string)
            call InputReadDouble(input,option,pm_well%grid%tophole(2))
            call InputErrorMsg(input,option,'TOP_OF_HOLE y-coordinate', &
                               error_string)
            call InputReadDouble(input,option,pm_well%grid%tophole(3))
            call InputErrorMsg(input,option,'TOP_OF_HOLE z-coordinate', &
                               error_string)
        !-----------------------------
          case('BOTTOM_OF_HOLE')
            call InputReadDouble(input,option,pm_well%grid%bottomhole(1))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE x-coordinate', &
                               error_string)
            call InputReadDouble(input,option,pm_well%grid%bottomhole(2))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE y-coordinate', &
                               error_string)
            call InputReadDouble(input,option,pm_well%grid%bottomhole(3))
            call InputErrorMsg(input,option,'BOTTOM_OF_HOLE z-coordinate', &
                               error_string)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      ! ----------------- error messaging -------------------------------------
      if (Uninitialized(pm_well%grid%nsegments)) then
        option%io_buffer = 'ERROR: NUMBER_OF_SEGMENTS must be specified &
                           &in the ' // trim(error_string) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Uninitialized(pm_well%grid%tophole(1)) .or. &
          Uninitialized(pm_well%grid%tophole(2)) .or. &
          Uninitialized(pm_well%grid%tophole(3))) then
        option%io_buffer = 'ERROR: TOP_OF_HOLE must be specified &
                           &in the ' // trim(error_string) // ' block.'
        call PrintMsg(option)
        num_errors = num_errors + 1
      endif
      if (Uninitialized(pm_well%grid%bottomhole(1)) .or. &
          Uninitialized(pm_well%grid%bottomhole(2)) .or. &
          Uninitialized(pm_well%grid%bottomhole(3))) then
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
                       &the WELLBORE_MODEL,GRID block. See above error messages.'
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
  use Option_module
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
          !Ideally (at least for Darcy model) we will link phys props to 
          !the material properties block
          case('PERMEABILITY')
            do k = 1,read_max
              call InputReadDouble(input,option,temp_well_perm(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for PERMEABILITY &
                &must be provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(pm_well%well%permeability(num_read))
            pm_well%well%permeability(1:num_read) = temp_well_perm(1:num_read)
        !-----------------------------
          case('POROSITY')
            do k = 1,read_max
              call InputReadDouble(input,option,temp_well_phi(k))
              if (InputError(input)) exit
              num_read = num_read + 1
            enddo
            if (num_read == 0) then
              option%io_buffer = 'At least one value for POROSITY &
                &must be provided in the ' // trim(error_string) // ' block.'
              call PrintErrMsg(option)
            endif
            allocate(pm_well%well%phi(num_read))
            pm_well%well%phi(1:num_read) = temp_well_phi(1:num_read)
        !-----------------------------
          case('WELL_INDEX_MODEL')
            call InputReadWord(input,option,pm_well%well%WI_model,PETSC_TRUE)
            call InputErrorMsg(input,option,'WELL_INDEX_MODEL',error_string)
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
  use Option_module
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
                    case('PRESSURE')
                      call InputReadDouble(input,option,this%well%bh_p)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('PRESSURE_SET_BY_RESERVOIR')
                      this%well%bh_p_set_by_reservoir = PETSC_TRUE
                  !-----------------------------
                    case('LIQUID_RATE')
                      call InputReadDouble(input,option,this%well%bh_ql)
                  !-----------------------------
                    case('GAS_RATE')
                      call InputReadDouble(input,option,this%well%bh_qg)
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
                    case('PRESSURE')
                      call InputReadDouble(input,option,this%well%th_p)
                      if (InputError(input)) exit
                  !-----------------------------
                    case('LIQUID_RATE')
                      call InputReadDouble(input,option,this%well%th_ql)
                  !-----------------------------
                    case('GAS_RATE')
                      call InputReadDouble(input,option,this%well%th_qg)
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
      !
      ! The following error checks are too restrictive:
      !
      !if (Uninitialized(this%well%th_p)) then
      !  if (Uninitialized(this%well%bh_p) .and. &
      !      .not.this%well%bh_p_set_by_reservoir)  then
      !    option%io_buffer = 'Keyword BOTTOM_OF_HOLE,PRESSURE/PRESSURE_SET_&
      !    &BY_RESERVOIR or keyword TOP_OF_HOLE,PRESSURE must be provided in &
      !    &the ' // trim(error_string) // ' block.'
      !    call PrintErrMsg(option) 
      !  endif
      !endif
      !if (Initialized(this%well%bh_p) .and. Initialized(this%well%th_p)) then
      !  option%io_buffer = 'Either keyword BOTTOM_OF_HOLE,PRESSURE or keyword &
      !    &TOP_OF_HOLE,PRESSURE must be provided in the ' &
      !    // trim(error_string) // ' block, but NOT BOTH.'
      !  call PrintErrMsg(option)
      !endif
      !if (this%well%bh_p_set_by_reservoir .and. &
      !    Initialized(this%well%th_p)) then
      !  option%io_buffer = 'Either keyword BOTTOM_OF_HOLE,PRESSURE_SET_BY_&
      !    &RESERVOIR or keyword TOP_OF_HOLE,PRESSURE must be provided in &
      !    &the ' // trim(error_string) // ' block, but NOT BOTH.'
      !  call PrintErrMsg(option)
      !endif

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
  use Option_module
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
            call InputReadInt(input,option,pm_well%flow_soln%ts_cut_factor)
            call InputErrorMsg(input,option,'TIMESTEP_CUT_FACTOR', &
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
  use Option_module
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
            call InputReadInt(input,option,pm_well%tran_soln%ts_cut_factor)
            call InputErrorMsg(input,option,'TIMESTEP_CUT_FACTOR', &
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

subroutine PMWellReadWellModelType(this,input,option,keyword,error_string,found)
  ! 
  ! Reads input file parameters associated with the well model 
  ! type.
  ! 
  ! Author: Michael Nole
  ! Date: 12/22/2021
  !
  use Input_Aux_module
  use Option_module
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
          case('CONSTANT_PRESSURE_HYDROSTATIC')
            this%well%well_model_type = 'CONSTANT_PRESSURE_HYDROSTATIC'
          case('CONSTANT_RATE')
            this%well%well_model_type = 'CONSTANT_RATE'
          case('DARCY')
            this%well%well_model_type = 'DARCY'
            !call PMWellReadDarcyInput(this,input,option,keyword,&
            !                           error_string,found)
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
    case('GRID','WELL','WELL_MODEL_TYPE','WELL_FLOW_SOLVER', &
         'WELL_TRANSPORT_SOLVER')
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
  call InputSkipToEND(input,option,card)

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
  
  implicit none

  class(pm_well_type) :: this
  
  PetscInt :: nsegments, k

  nsegments = this%grid%nsegments 

  allocate(this%flow_soln%residual(nsegments*this%flow_soln%ndof))
  allocate(this%flow_soln%update(nsegments*this%flow_soln%ndof))
  this%flow_soln%residual(:) = UNINITIALIZED_DOUBLE
  this%flow_soln%update(:) = UNINITIALIZED_DOUBLE

  allocate(this%flow_soln%Jacobian(this%nphase*nsegments,this%nphase*nsegments))
  this%flow_soln%Jacobian = UNINITIALIZED_DOUBLE

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
                                                   this%grid%nsegments))
    allocate(this%tran_soln%prev_soln%aqueous_mass(this%nspecies, &
                                                   this%grid%nsegments))
  endif

  allocate(this%well%WI(nsegments))
  allocate(this%well%pl(nsegments))
  allocate(this%well%pg(nsegments))
  allocate(this%well%ql(nsegments-1))
  allocate(this%well%qg(nsegments-1))
  if (this%transport) then
    allocate(this%well%aqueous_conc(this%nspecies,nsegments))
    allocate(this%well%aqueous_mass(this%nspecies,nsegments))
  endif

  allocate(this%pert(nsegments,this%nphase))
  this%pert = 0.d0

  allocate(this%well%liq%s(nsegments))
  this%well%liq%rho0 = this%option%flow%reference_density(1)
  allocate(this%well%liq%rho(nsegments))
  this%well%liq%rho(:) = this%well%liq%rho0
  allocate(this%well%liq%visc(nsegments))
  this%well%liq%visc(:) = 2.1d-3
  allocate(this%well%liq%Q(nsegments))
  if (this%nphase == 2) then
    allocate(this%well%gas%s(nsegments))
    this%well%gas%rho0 = this%option%flow%reference_density(2)
    allocate(this%well%gas%rho(nsegments))
    this%well%gas%rho(:) = this%well%gas%rho0
    allocate(this%well%gas%visc(nsegments))
    this%well%gas%visc(:) = 8.9339d-6
    allocate(this%well%gas%Q(nsegments))
  endif

  ! Initialize perturbations
  allocate(this%pert(nsegments,this%nphase))
  this%pert = 0.d0

  ! Initialize stored aux variables at well perturbations
  do k = 1,this%nphase
    this%well_pert(k)%well_model_type = this%well%well_model_type
    this%well_pert(k)%wi_model = this%well%wi_model

    allocate(this%well_pert(k)%WI(nsegments))
    allocate(this%well_pert(k)%pl(nsegments))
    allocate(this%well_pert(k)%pg(nsegments))
    allocate(this%well_pert(k)%diameter(nsegments))
    allocate(this%well_pert(k)%WI_base(nsegments))
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

    this%well_pert(k)%liq%rho0 = this%option%flow%reference_density(1)
    this%well_pert(k)%liq%rho(:) = this%well%liq%rho0
    this%well_pert(k)%liq%visc(:) = this%well%liq%visc(:)
    this%well_pert(k)%liq%mobility = this%well%liq%mobility
    this%well_pert(k)%gas%rho0 = this%option%flow%reference_density(2)
    this%well_pert(k)%gas%rho(:) = this%well%gas%rho0
    this%well_pert(k)%gas%visc(:) = this%well%gas%visc(:)
    this%well_pert(k)%gas%mobility = this%well%gas%mobility
    this%well_pert(k)%diameter(:) = this%well%diameter(:)
    this%well_pert(k)%WI_base(:) = this%well%WI_base(:)
    this%well_pert(k)%permeability(:) = this%well%permeability(:)
    this%well_pert(k)%phi(:) = this%well%phi(:)
    this%well_pert(k)%f(:) = this%well%f(:)
    this%well_pert(k)%area(:) = this%well%area(:)
    this%well_pert(k)%volume(:) = this%well%volume(:)
    this%well_pert(k)%bh_p = this%well%bh_p
    this%well_pert(k)%th_p = this%well%th_p
    this%well_pert(k)%bh_ql = this%well%bh_ql
    this%well_pert(k)%bh_qg = this%well%bh_qg
    this%well_pert(k)%th_ql = this%well%th_ql
    this%well_pert(k)%th_qg = this%well%th_qg
    this%well_pert(k)%liq%visc(:) = this%well%liq%visc(:)
    this%well_pert(k)%gas%visc(:) = this%well%gas%visc(:)
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
  if (this%transport) then
    allocate(this%reservoir%aqueous_conc(this%nspecies,this%grid%nsegments))
    allocate(this%reservoir%aqueous_mass(this%nspecies,this%grid%nsegments))
  endif

  call PMWellOutputHeader(this)
  
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
  
  implicit none
  
  class(pm_well_type) :: this

  ! update the reservoir object's pressure and saturations
  call PMWellUpdateReservoir(this)

  if (initialize_well) then
    this%well%pl = this%reservoir%p_l
    this%well%pg = this%reservoir%p_g
    this%well%liq%s = this%reservoir%s_l
    this%well%gas%s = this%reservoir%s_g
    this%well%liq%rho = this%reservoir%rho_l
    this%well%gas%rho = this%reservoir%rho_g
    this%well%liq%visc = this%reservoir%visc_l
    this%well%gas%visc = this%reservoir%visc_g
    this%well_pert(ONE_INTEGER)%pl = this%reservoir%p_l
    this%well_pert(TWO_INTEGER)%pl = this%reservoir%p_l
    this%well_pert(ONE_INTEGER)%pg = this%reservoir%p_g
    this%well_pert(TWO_INTEGER)%pg = this%reservoir%p_g
    this%well_pert(ONE_INTEGER)%liq%s = this%reservoir%s_l
    this%well_pert(TWO_INTEGER)%liq%s = this%reservoir%s_l
    this%well_pert(ONE_INTEGER)%gas%s = this%reservoir%s_g
    this%well_pert(TWO_INTEGER)%gas%s = this%reservoir%s_g
    this%well_pert(ONE_INTEGER)%liq%rho = this%reservoir%rho_l
    this%well_pert(ONE_INTEGER)%gas%rho = this%reservoir%rho_g
    this%well_pert(TWO_INTEGER)%liq%rho = this%reservoir%rho_l
    this%well_pert(TWO_INTEGER)%gas%rho = this%reservoir%rho_g
    this%well_pert(ONE_INTEGER)%liq%visc = this%reservoir%visc_l
    this%well_pert(ONE_INTEGER)%gas%visc = this%reservoir%visc_g
    this%well_pert(TWO_INTEGER)%liq%visc = this%reservoir%visc_l
    this%well_pert(TWO_INTEGER)%gas%visc = this%reservoir%visc_g

    initialize_well = PETSC_FALSE
    this%flow_soln%prev_soln%pl = this%well%pl
    this%flow_soln%prev_soln%sg = this%well%gas%s
    
  !else
  endif

  call PMWellUpdateProperties(this%well)
  !MAN: need to add adaptive timestepping
  this%dt = this%realization%option%flow_dt

end subroutine PMWellInitializeTimestep

! ************************************************************************** !

subroutine PMWellUpdateReservoir(this)
  !
  ! Updates the reservoir properties for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

  use WIPP_Flow_Aux_module
  use Material_Aux_class
  use Grid_module

  implicit none
  
  class(pm_well_type) :: this

  type(wippflo_auxvar_type), pointer :: wippflo_auxvar
  type(material_auxvar_type), pointer :: material_auxvar
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscInt :: k
  PetscInt :: ghosted_id

  option => this%option 

  grid => this%realization%patch%grid

  do k=1,size(this%grid%h_ghosted_id)
    ghosted_id = this%grid%h_ghosted_id(k)

    wippflo_auxvar => &
      this%realization%patch%aux%wippflo%auxvars(0,ghosted_id)

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
    ! note: the following other things available in wippflo_auxvar object:
    ! xmol, alpha, elevation, fracture_perm_scaling_factor, klinkenberg 

    this%reservoir%kx(k) = material_auxvar%permeability(1)
    this%reservoir%ky(k) = material_auxvar%permeability(2)
    this%reservoir%kz(k) = material_auxvar%permeability(3)

    if (associated(grid%structured_grid)) then
      this%reservoir%dx(k) = grid%structured_grid%dx(ghosted_id) 
      this%reservoir%dy(k) = grid%structured_grid%dy(ghosted_id)
      this%reservoir%dz(k) = grid%structured_grid%dz(ghosted_id)
    else
      option%io_buffer = 'Well model is not yet compatible with an unstructured &
          &grid reservoir. '
      call PrintErrMsg(option)
    endif

    if ((k == 1) .and. this%well%bh_p_set_by_reservoir) then
      this%well%bh_p = this%reservoir%p_l(k)
    endif
  enddo

end subroutine PMWellUpdateReservoir

! ************************************************************************** !

subroutine PMWellUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                               num_newton_iterations,tfac, &
                               time_step_max_growth_factor)
  !
  ! Updates the time step for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  
  implicit none
  
  class(pm_well_type) :: this
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

  call PMWellUpdateReservoirSrcSink(this)

  call PMWellOutput(this)
  
end subroutine PMWellFinalizeTimestep

! ************************************************************************** !

subroutine PMWellUpdateReservoirSrcSink(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 12/21/2021

  use Coupler_module

  implicit none
  
  class(pm_well_type) :: this

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: srcsink_name
  type(coupler_type), pointer :: source_sink
  PetscInt :: k

  do k = 1,this%grid%nsegments
    write(string,'(I0.4)') k
    srcsink_name = 'well_segment_' // trim(string)

    source_sink => this%realization%patch%source_sink_list%first 
    do 
      if (.not.associated(source_sink)) exit

      if (trim(srcsink_name) == trim(source_sink%name)) then
        source_sink%flow_condition%general%rate%dataset%rarray(1) = &
          this%well%liq%Q(k) ! [kg/s]
        source_sink%flow_condition%general%rate%dataset%rarray(2) = &
          this%well%gas%Q(k) ! [kg/s]
        exit
      endif

      source_sink => source_sink%next
    enddo

  enddo

end subroutine PMWellUpdateReservoirSrcSink

! ************************************************************************** !

subroutine PMWellRhoIncompress(fluid)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 11/09/2021

  implicit none
  
  type(well_fluid_type) :: fluid

  fluid%rho(:) = fluid%rho0

end subroutine PMWellRhoIncompress

! ************************************************************************** !

subroutine PMWellResidual(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  !SNES :: snes
  !Vec :: xx
  !Vec :: r
  !PetscErrorCode :: ierr
  PetscInt :: i, k, iup, idn
  PetscReal :: res_accum(this%nphase)
  PetscReal :: res_flux(this%nphase)
  PetscReal :: res_flux_bc(2*this%nphase)

  res_accum = 0.d0
  res_flux = 0.d0
  res_flux_bc = 0.d0

  select case(this%well%well_model_type)
    !-------------------------------------------------------------------------
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    !-------------------------------------------------------------------------
    case('DARCY')

      call PMWellBCFlux(this,this%well,res_flux_bc)

      do i = 1,this%grid%nsegments
        iup = i
        idn = i + 1
        call PMWellFlux(this,this%well,this%well,iup,idn,res_flux)
        call PMWellAccumulation(this,this%well,i,res_accum)

        this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) = &
             this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) + &
             res_accum(ONE_INTEGER)
        if (i == 1) then
          ! Water mass residual in cell i: Add flux in from BC, 
          ! subtract flux to i+1 cell
          this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) &
               + res_flux_bc(1) - res_flux(1)
          ! Water mass residual in cell i+1: Add flux to i+1 cell
          this%flow_soln%residual(this%flow_soln%ndof*i+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*i+1) &
               + res_flux(1)

        elseif (i < this%grid%nsegments) then
          ! Water mass residual in cell i: Subtract flux to i+1 cell
          this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) &
               - res_flux(1) 
          ! Water mass residual in cell i+1: Add flux to i+1 cell
          this%flow_soln%residual(this%flow_soln%ndof*i+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*i+1) &
               + res_flux(1)
        else
          ! Water mass residual in cell i: Subtract flux to BC
          this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) = &
               this%flow_soln%residual(this%flow_soln%ndof*(i-1)+1) &
               - res_flux_bc(1)
        endif

        if (this%nphase == 2) then
          this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) = &
               this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) + &
               res_accum(TWO_INTEGER)
          if (i == 1) then
            ! Air mass residual in cell i: Add flux in from BC, 
            ! subtract flux to i+1 cell
            this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) = &
                 this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) &
                 + res_flux_bc(2) - res_flux(2)
            ! Air mass residual in cell i+1: Add flux to i+1 cell
            this%flow_soln%residual(this%flow_soln%ndof*i+2) = &
                 this%flow_soln%residual(this%flow_soln%ndof*i+2) &
                 + res_flux(2)
          elseif (i < this%grid%nsegments) then
            ! Air mass residual in cell i: Subtract flux to i+1 cell
            this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) = &
                 this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) & 
                 - res_flux(2) 
            ! Air mass residual in cell i+1: Add flux to i+1 cell
            this%flow_soln%residual(this%flow_soln%ndof*i+2) = &
                 this%flow_soln%residual(this%flow_soln%ndof*i+2) &
                 + res_flux(2)
          else
            ! Air mass residual in cell i: Subtract flux to BC
            this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) = &
                 this%flow_soln%residual(this%flow_soln%ndof*(i-1)+2) &
                 - res_flux_bc(2)
          endif
        endif
      enddo
    !-------------------------------------------------------------------------
    case('FULL_MOMENTUM')
    !-------------------------------------------------------------------------
  end select

end subroutine PMWellResidual

! ************************************************************************** !

subroutine PMWellJacobian(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none

  class(pm_well_type) :: this

  !SNES :: snes
  !Vec :: xx
  Mat :: A, B
  !type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscReal :: norm
  PetscViewer :: viewer

  PetscReal :: qsrc, scale
  PetscInt :: imat, imat_up, imat_dn
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: irow, iconn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  Vec, parameter :: null_vec = tVec(0)

  PetscReal :: Jup(this%nphase,this%nphase), &
               Jdn(this%nphase,this%nphase), &
               Jtop(this%nphase,this%nphase), &
               Jbtm(this%nphase,this%nphase), &
               Jtmp(this%nphase,this%nphase), &
               Jac(this%nphase*this%grid%nsegments, &
                   this%nphase*this%grid%nsegments)
               
  this%flow_soln%Jacobian = 0.d0
  Jac = 0.d0

  select case(this%well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('DARCY')
      ! Not expecting ghosting at this time (1D model)
      !if (.not. well_analytical_derivatives) then
      !  call PMWellPerturb(this)
      !endif
      ! Accumulation terms ------------------------------------
      do local_id = 1,this%grid%nsegments
        call PMWellAccumDerivative(this,local_id,Jup)
        call PMWellFillJac(this,Jac,Jup,local_id,local_id)
      enddo

      ! Interior Flux Terms ----------------------------------- 
      do iconn = 1,this%grid%nconnections

        local_id_up = iconn
        local_id_dn = iconn+1

        call PMWellFluxDerivative(this,local_id_up,local_id_dn,Jup,Jdn)

        Jtmp = Jup
        call PMWellFillJac(this,Jac,Jtmp,local_id_up,local_id_up)

        Jtmp = Jdn
        call PMWellFillJac(this,Jac,Jtmp,local_id_up,local_id_dn)

        Jup = -Jup
        Jdn = -Jdn
        Jtmp = Jdn
        call PMWellFillJac(this,Jac,Jtmp,local_id_dn,local_id_dn)

        Jtmp = Jup
        call PMWellFillJac(this,Jac,Jtmp,local_id_dn,local_id_up)

      enddo

      ! Boundary Flux Terms -----------------------------------
      local_id = 1
      call PMWellBCFluxDerivative(this,Jtop,Jbtm)
      Jbtm = -Jbtm
      call PMWellFillJac(this,Jac,Jbtm,local_id,local_id)

      local_id = this%grid%nsegments
      Jtop = -Jtop
      call PMWellFillJac(this,Jac,Jtop,local_id,local_id)

      !pm_well_ni_count = pm_well_ni_count + 1

    case('FULL_MOMENTUM')

  end select

  this%flow_soln%Jacobian = Jac

end subroutine PMWellJacobian

! ************************************************************************** !
subroutine PMWellFillJac(pm_well,Jac,Jtmp,id1,id2)
  ! 
  ! Author: Michael Nole
  ! Date: 01/10/2022
  !

  implicit none

  class(pm_well_type) :: pm_well
  PetscReal :: Jtmp(pm_well%nphase,pm_well%nphase), &
               Jac(pm_well%nphase*pm_well%grid%nsegments, &
                   pm_well%nphase*pm_well%grid%nsegments)
  PetscInt :: id1,id2

  PetscInt :: i,j

  do i = 1,pm_well%nphase
    do j = 1,pm_well%nphase
      Jac((id1-1)*pm_well%nphase+i,(id2-1)*pm_well%nphase+j) = &
      Jac((id1-1)*pm_well%nphase+i,(id2-1)*pm_well%nphase+j) + Jtmp(i,j)
    enddo
  enddo

end subroutine PMWellFillJac

! ************************************************************************** !

subroutine PMWellPreSolve(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: cur_time
  
  this%flow_soln%not_converged = PETSC_TRUE
  this%flow_soln%converged = PETSC_FALSE

  cur_time = this%option%time - this%option%flow_dt + this%cumulative_dt

  write(out_string,'(" Step ",i6,"   Time =",1pe12.5,"   Dt =", &
                     1pe12.5," sec.")') &
                   (this%flow_soln%n_steps+1),cur_time,this%dt  
  call OptionPrint(out_string,this%option)
  
end subroutine PMWellPreSolve

! ************************************************************************** !

subroutine PMWellSolve(this,time,ierr)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

  implicit none
  
  class(pm_well_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(well_soln_flow_type), pointer :: flow_soln
  character(len=MAXSTRINGLENGTH) :: out_string
  PetscLogDouble :: log_start_time, log_end_time
  PetscInt :: n_iter,ts_cut,i,easy_converge_count
  PetscReal :: res(this%flow_soln%ndof)
  PetscReal :: res_fixed(this%flow_soln%ndof*this%grid%nsegments)

  flow_soln => this%flow_soln

  ierr = 0
  call PetscTime(log_start_time, ierr);CHKERRQ(ierr)

  n_iter = 0
  ts_cut = 0
  easy_converge_count = 0

  this%cumulative_dt = 0.d0
  flow_soln%converged = PETSC_FALSE
  flow_soln%not_converged = PETSC_TRUE

  do while (this%cumulative_dt < this%realization%option%flow_dt) 

    call PMWellPreSolve(this)

    ! Fixed accumulation term
    res_fixed = 0.d0
    do i = 1,this%grid%nsegments
      call PMWellAccumulation(this,this%well,i,res)
      res_fixed(flow_soln%ndof*(i-1)+1:flow_soln%ndof*i) = -1.d0 * res 
    enddo

    n_iter = 0

    do while (flow_soln%not_converged)

      if (n_iter > (flow_soln%max_iter-1)) then
        flow_soln%cut_timestep = PETSC_TRUE
        out_string = ' Maximum number of Newton iterations reached. Cutting &
                      &timestep!'
        call OptionPrint(out_string,this%option); WRITE(*,*) ""
        call PMWellCutTimestep(this)
        n_iter = 0
        ts_cut = ts_cut + 1
        easy_converge_count = 0
        exit
      endif
      if (ts_cut > flow_soln%max_ts_cut) then
        this%realization%option%io_buffer = &
          'Maximum timestep cuts reached in PM Well. Solution has not &
           &converged. Exiting.'
        call PrintErrMsg(this%realization%option)
      endif

      flow_soln%residual = 0.d0
      flow_soln%residual = res_fixed

      easy_converge_count = easy_converge_count + 1

      ! use Newton's method to solve for the well pressure and saturation
      call PMWellNewton(this)
      ! Declare convergence
      call PMWellCheckConvergence(this,n_iter,res_fixed)
    enddo

    if (.not. flow_soln%not_converged .and. flow_soln%cut_timestep .and. &
         easy_converge_count > 10 ) then
      if (this%cumulative_dt + this%dt * this%flow_soln%ts_cut_factor < &
          this%realization%option%flow_dt) then
        this%dt = this%dt * this%flow_soln%ts_cut_factor 
        flow_soln%cut_timestep = PETSC_FALSE
      endif
    endif

    if (this%cumulative_dt + this%dt > this%realization%option%flow_dt) then
      this%dt = this%realization%option%flow_dt - this%cumulative_dt
    endif

    ts_cut = 0
    flow_soln%n_steps = flow_soln%n_steps + 1
 
  enddo

  call PMWellPostSolve(this)

  ! update well index (should not need to be updated every time if
  ! grid permeability and discretization are static
  call PMWellComputeWellIndex(this)

  ! update the well src/sink Q vector
  call PMWellUpdateWellQ(this%well,this%reservoir)

  ! update the Darcy fluxes within the well
  call PMWellCalcVelocity(this)

  call PetscTime(log_end_time, ierr);CHKERRQ(ierr)
  
end subroutine PMWellSolve

! ************************************************************************** !

subroutine PMWellUpdateSolution(pm_well)
  ! 
  ! Author: Michael Nole
  ! Date: 01/21/22

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: i
  PetscInt :: idof

  select case(pm_well%well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('DARCY')
      do i = 1,pm_well%grid%nsegments
        idof = pm_well%flow_soln%ndof*(i-1)
        pm_well%well%pl(i) = pm_well%well%pl(i) +  &
                             pm_well%flow_soln%update(idof+1)
        idof = pm_well%flow_soln%ndof*i
        pm_well%well%gas%s(i) = pm_well%well%gas%s(i) + &
                                pm_well%flow_soln%update(idof)
      enddo
      call PMWellUpdateProperties(pm_well%well)
  end select
end subroutine PMWellUpdateSolution

! ************************************************************************** !

subroutine PMWellCutTimestep(pm_well)
  ! 
  ! Author: Michael Nole
  ! Date: 01/24/22

  implicit none

  class(pm_well_type) :: pm_well

  PetscInt :: i

  select case(pm_well%well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('DARCY')
      ! could make this smarter or call smarter timestepping routines
      pm_well%dt = pm_well%dt / pm_well%flow_soln%ts_cut_factor
      pm_well%dt = max(pm_well%dt,pm_well%min_dt)
      pm_well%well%pl = pm_well%flow_soln%prev_soln%pl
      pm_well%well%gas%s = pm_well%flow_soln%prev_soln%sg
      call PMWellUpdateProperties(pm_well%well)
  end select

end subroutine PMWellCutTimestep

! ************************************************************************** !

subroutine PMWellNewton(this)
  ! 
  ! Author: Michael Nole
  ! Date: 01/20/2022

  use Utility_module

  implicit none
  
  class(pm_well_type) :: this

  PetscReal :: identity(this%nphase*this%grid%nsegments,&
                        this%nphase*this%grid%nsegments)
  PetscReal :: inv_Jac(this%nphase*this%grid%nsegments, &
                       this%nphase*this%grid%nsegments), &
               Jac_save(this%nphase*this%grid%nsegments, &
                       this%nphase*this%grid%nsegments)
  PetscReal :: new_dx(this%nphase*this%grid%nsegments), &
               res_save(this%nphase*this%grid%nsegments)
  PetscInt :: indx(this%nphase*this%grid%nsegments)
  PetscInt :: i,j
  PetscInt :: d

  res_save = 0.d0
  Jac_save = 0.d0 

  call PMWellPerturb(this) 
 
  call PMWellResidual(this)

  call PMWellJacobian(this)

  select case (this%well%well_model_type)
  !--------------------------------------
  case('CONSTANT_PRESSURE')
    ! No capillarity yet
    this%well%pl(:) = this%well%bh_p
    this%well%pg(:) = this%well%bh_p
  !--------------------------------------
  case('DARCY')
    do i = 1,this%nphase*this%grid%nsegments
      do j = 1,this%nphase*this%grid%nsegments
        if (i==j) then
          identity(i,j) = 1.d0
        else
          identity(i,j) = 0.d0
        endif
      enddo
    enddo
    Jac_save = this%flow_soln%Jacobian
    call LUDecomposition(this%flow_soln%Jacobian,this%nphase*this%grid%nsegments, &
                         indx,d)
    res_save = this%flow_soln%residual
    call LUBackSubstitution(this%flow_soln%Jacobian, &
                            this%nphase*this%grid%nsegments,&
                            indx,this%flow_soln%residual)
    new_dx = this%flow_soln%residual
    this%flow_soln%update(:) = new_dx(:)

    this%flow_soln%residual = res_save
    this%flow_soln%Jacobian = Jac_save

    call PMWellUpdateSolution(this)

  !--------------------------------------
  end select

end subroutine PMWellNewton

! ************************************************************************** !

subroutine PMWellPostSolve(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this

  character(len=MAXSTRINGLENGTH) :: out_string

  WRITE(out_string,'(" PM Well Step Complete!    Time=",1pe12.5," sec &
                     & Total Newton Iterations =",i8)') &
                    this%option%time,this%flow_soln%n_newton 
  call OptionPrint(out_string,this%option)
  WRITE(*,*) ""
  
end subroutine PMWellPostSolve

! ************************************************************************** !

subroutine PMWellCheckConvergence(this,n_iter,fixed_accum)
  ! 
  ! Checks solution convergence against prescribed tolerances.
  !
  ! Author: Jennifer M. Frederick
  ! Date: 01/20/2022

  implicit none
  
  class(pm_well_type) :: this
  PetscInt :: n_iter
  PetscReal :: fixed_accum(this%flow_soln%ndof*this%grid%nsegments)
  
  type(well_soln_flow_type), pointer :: flow_soln
  character(len=MAXSTRINGLENGTH) :: out_string
  character(len=MAXSTRINGLENGTH) :: rsn_string
  PetscBool :: cnvgd_due_to_residual(this%grid%nsegments*this%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_abs_res(this%grid%nsegments*this%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_scaled_res(this%grid%nsegments*this%flow_soln%ndof)
  PetscBool :: cnvgd_due_to_update(this%grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update_p(this%grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update_s(this%grid%nsegments)
  PetscBool :: cnvgd_due_to_abs_update(this%grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update_p(this%grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update_s(this%grid%nsegments)
  PetscBool :: cnvgd_due_to_rel_update(this%grid%nsegments)
  PetscReal :: update_p(this%grid%nsegments) ! liquid pressure
  PetscReal :: update_s(this%grid%nsegments) ! gas saturation
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
  rsn_string = ''

  do k = 1,this%grid%nsegments
    idof = this%flow_soln%ndof*(k-1)
    update_s(k) = flow_soln%update(idof)
    update_p(k) = flow_soln%update(idof+1)

    ! Absolute Solution Updates
    temp_real = dabs(update_p(k))
    if (temp_real < flow_soln%itol_abs_update_p) then
      cnvgd_due_to_abs_update_p(k) = PETSC_TRUE
    endif
    temp_real = dabs(update_s(k))
    if (temp_real < flow_soln%itol_abs_update_s) then
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
  enddo

  do k = 1,(this%grid%nsegments*flow_soln%ndof)
    ! Absolute Residual
    temp_real = dabs(flow_soln%residual(k))
    if (temp_real < flow_soln%itol_abs_res) then
      cnvgd_due_to_abs_res(k) = PETSC_TRUE
    endif
 
    ! Scaled Residual
    temp_real = dabs(flow_soln%residual(k)/(fixed_accum(k)/this%dt))
    if (temp_real < flow_soln%itol_scaled_res) then
      cnvgd_due_to_scaled_res(k) = PETSC_TRUE
    endif
  enddo

  max_absolute_residual = maxval(dabs(flow_soln%residual))
  loc_max_abs_residual = maxloc(dabs(flow_soln%residual),1)

  max_scaled_residual = maxval(dabs(flow_soln%residual/ &
                                    (fixed_accum/this%dt)))
  loc_max_scaled_residual = maxloc(dabs(flow_soln%residual/ &
                                        (fixed_accum/this%dt)),1)

  max_absolute_update_p = maxval(dabs(update_p))
  loc_max_abs_update_p = maxloc(dabs(update_p),1)

  max_absolute_update_s = maxval(dabs(update_s))
  loc_max_abs_update_s = maxloc(dabs(update_s),1)

  max_relative_update_p = maxval(dabs(update_p/this%well%pl))
  loc_max_rel_update_p = maxloc(dabs(update_p/this%well%pl),1)

  max_relative_update_s = maxval(dabs(update_s/this%well%gas%s))
  loc_max_rel_update_s = maxloc(dabs(update_s/this%well%gas%s),1)

  do k = 1,(this%grid%nsegments*flow_soln%ndof)
    if (cnvgd_due_to_scaled_res(k) .or. cnvgd_due_to_abs_res(k)) then
      cnvgd_due_to_residual(k) = PETSC_TRUE
    endif
  enddo
  do k = 1,(this%grid%nsegments)
    if (cnvgd_due_to_abs_update_p(k) .and. &
        cnvgd_due_to_abs_update_s(k)) then
      cnvgd_due_to_abs_update(k) = PETSC_TRUE
    endif
    if (cnvgd_due_to_rel_update_p(k) .and. &
        cnvgd_due_to_rel_update_s(k)) then
      cnvgd_due_to_rel_update(k) = PETSC_TRUE
    endif
    if (cnvgd_due_to_abs_update(k) .or. cnvgd_due_to_rel_update(k)) then
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
        ,es10.2,"  uS:",es10.2,"  ruP:",es10.2,"  ruS:",es10.2)') &
        n_iter,max_absolute_residual,max_scaled_residual, &
        max_absolute_update_p,max_absolute_update_s, &
        max_relative_update_p,max_relative_update_s  
  call OptionPrint(out_string,this%option)

  if (all(cnvgd_due_to_residual) .and. all(cnvgd_due_to_update)) then
    flow_soln%converged = PETSC_TRUE
    flow_soln%not_converged = PETSC_FALSE
    out_string = ' Solution converged!  ---> ' // trim(rsn_string)
    call OptionPrint(out_string,this%option); WRITE(*,*) ""
    this%cumulative_dt = this%cumulative_dt + this%dt
  else
    flow_soln%converged = PETSC_FALSE
    flow_soln%not_converged = PETSC_TRUE
  endif

  
end subroutine PMWellCheckConvergence

! ************************************************************************** !

subroutine PMWellUpdateWellQ(well,reservoir)
  !
  ! Updates the src/sink vector for the fluid object.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 12/01/2021

  implicit none
  
  type(well_type), pointer :: well
  type(well_reservoir_type), pointer :: reservoir

  type(well_fluid_type), pointer :: liq
  type(well_fluid_type), pointer :: gas

  liq => well%liq
  gas => well%gas

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
    case default 
      ! Flowrate in kg/s
      liq%Q = liq%rho*liq%mobility*well%WI*(reservoir%p_l-well%pl)
      gas%Q = gas%rho*gas%mobility*well%WI*(reservoir%p_g-well%pg)
    !------------------------------------------------------------------------
  end select
  
end subroutine PMWellUpdateWellQ

! ************************************************************************** !

subroutine PMWellCalcVelocity(this)
  !
  ! Calculates the Darcy flux in the well given the well pressure.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 02/16/2022

  implicit none
  
  class(pm_well_type) :: this

  type(well_grid_type), pointer :: grid
  type(well_type), pointer :: well
  PetscReal :: perm_up, perm_dn
  PetscReal :: dist_up, dist_dn
  PetscReal :: perm_ave_over_dist
  PetscReal :: density_kg_ave
  PetscReal :: gravity_term
  PetscReal :: delta_pressure
  PetscInt :: iup, idn
  PetscInt :: k

  PetscReal, parameter :: eps = 1.d-8

  grid => this%grid
  well => this%well 

  ! iup/idn convention:
  ! up direction goes towards well bottom (towards lower k value)
  ! dn direction goes towards well top (towards higher k value)

  ! the routines for calculating velocity must match PMWellFlux() exactly
  
  select case(well%well_model_type)
    !-------------------------------------------------------------------------
    case('CONSTANT_PRESSURE_HYDROSTATIC')
    !-------------------------------------------------------------------------
    case('CONSTANT_PRESSURE')
    !-------------------------------------------------------------------------
    case('CONSTANT_RATE')
    !-------------------------------------------------------------------------
    case('DARCY')
      do k = 1,(grid%nsegments-1)
        iup = k
        idn = k+1

        perm_up = well%permeability(iup)
        perm_dn = well%permeability(idn)
        dist_up = grid%dh(iup)/2.d0
        dist_dn = grid%dh(idn)/2.d0

        perm_ave_over_dist = (perm_up * perm_dn) / &
                             (dist_up*perm_dn + dist_dn*perm_up)

        ! ------- liquid Darcy flux -------
        density_kg_ave = 0.5d0*(well%liq%rho(iup)+well%liq%rho(idn))
        gravity_term = density_kg_ave*gravity*grid%dh(iup)
        delta_pressure = well%pl(iup) - well%pl(idn) + gravity_term
        well%ql(k) = perm_ave_over_dist*well%liq%mobility*delta_pressure


        ! ------- gas Darcy flux ----------
        density_kg_ave = 0.5d0*(well%gas%rho(iup)+well%gas%rho(idn))
        gravity_term = density_kg_ave*gravity*grid%dh(iup)
        delta_pressure = well%pg(iup) - well%pg(idn) + gravity_term
        well%qg(k) = perm_ave_over_dist*well%gas%mobility*delta_pressure

      enddo
    !-------------------------------------------------------------------------
    case('FULL_MOMENTUM')
    !-------------------------------------------------------------------------
  end select
  
end subroutine PMWellCalcVelocity

! ************************************************************************** !

subroutine PMWellComputeWellIndex(this)
  !
  ! Computes the well index.
  ! 
  ! Author: Michael Nole
  ! Date: 12/22/21

  implicit none

  class(pm_well_type) :: this
  PetscReal :: r0(this%grid%nsegments)
  PetscReal, parameter :: PI=3.141592653589793d0

  ! Peaceman Model: default = anisotropic
  ! This assumes z is vertical (not true for WIPP)
  select case(this%well%WI_model)
    case('PEACEMAN_ISO')
      this%well%WI = 2.d0*PI*this%reservoir%kx*this%grid%dh/ &
                     (log(2.079d-1*this%reservoir%dx/ &
                      (this%well%diameter/2.d0)))
    case('PEACEMAN_ANISOTROPIC')
      r0 = 2.8d-1*(sqrt(sqrt(this%reservoir%ky/this%reservoir%kx)* &
           this%reservoir%dx**2 + sqrt(this%reservoir%kx/this%reservoir%ky)* &
           this%reservoir%dy**2) / ((this%reservoir%ky/this%reservoir%kx)** &
           2.5d-1 + (this%reservoir%kx/this%reservoir%ky)**2.5d-1))
      this%well%WI = 2.d0*PI*sqrt(this%reservoir%kx*this%reservoir%ky)* &
                     this%grid%dh/log((this%well%diameter/2.d0)/r0)
  end select

  this%well%WI = this%well%WI*this%well%WI_base


end subroutine PMWellComputeWellIndex

! ************************************************************************** !

subroutine PMWellAccumulation(pm_well,well,id,Res)
  !
  ! Computes the accumulation term for the residual based on the 
  ! chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 12/23/2021
  !

  implicit none

  type(pm_well_type) :: pm_well
  type(well_type) :: well

  PetscReal :: Res(pm_well%nphase)
  PetscInt :: id,iphase

  Res = 0.d0

  select case(well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('DARCY')
      ! liquid accumulation term
      Res(1) = (Res(1) + well%liq%s(id) * well%liq%rho(id)) / FMWH2O * &
                well%phi(id) * well%volume(id) / pm_well%dt 
      ! gas accumulation term
      Res(2) = (Res(2) + well%gas%s(id) * well%gas%rho(id)) / &
                fmw_comp(TWO_INTEGER) * well%phi(id) * well%volume(id) / &
                pm_well%dt
    case default
  end select

end subroutine PMWellAccumulation

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

  call PMWellAccumulation(pm_well,pm_well%well,local_id,res)

  do idof = 1, pm_well%nphase
    call PMWellAccumulation(pm_well,pm_well%well_pert(idof),local_id,res_pert)
    do irow = 1, pm_well%nphase
      Jac(irow,idof) = (res_pert(irow)-res(irow))/pm_well%pert(local_id,idof)
    enddo !irow
  enddo ! idof

end subroutine PMWellAccumDerivative

! ************************************************************************** !

subroutine PMWellFlux(pm_well,well_up,well_dn,iup,idn,Res)
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

  type(well_grid_type), pointer :: grid

  PetscInt :: i, ghosted_id
  PetscReal :: pres_up, pres_dn

  PetscReal :: perm_ave_over_dist(2), perm_rho_mu_area_up(2), &
               perm_rho_mu_area_dn(2)
  PetscReal :: perm_up, perm_dn, dist_up, dist_dn, density_kg_ave
  PetscReal :: gravity_term, delta_pressure, v_darcy
  PetscReal :: density_ave_kmol, q, tot_mole_flux
  PetscReal :: up_scale, dn_scale
  PetscBool :: upwind

  PetscReal, parameter :: eps = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

  grid => pm_well%grid

  Res(:) = 0.d0

  select case(pm_well%well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('DARCY')
      ! This is good for either: single-phase liquid, or two-phase liquid/gas.
      ! Vertical well, no Klinkenberg, no capillary pressure, constant mobility.
        !MAN: need to check that this is the correct up/dn orientation

        perm_up = well_up%permeability(iup)
        perm_dn = well_dn%permeability(idn)
        dist_up = grid%dh(iup)/2.d0
        dist_dn = grid%dh(idn)/2.d0

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

        perm_ave_over_dist = (perm_rho_mu_area_up * perm_rho_mu_area_dn) / &
                     (dist_up*perm_rho_mu_area_dn + dist_dn*perm_rho_mu_area_up)

        ! Liquid flux
        ! Mobility may need to be non-constant?
        if (well_up%liq%mobility > eps) then

          density_kg_ave = 0.5d0*(well_up%liq%rho(iup)+well_dn%liq%rho(idn))
          ! Assuming the well is always vertical and gravity is in the
          ! (-) directiup
          ! And assuming dh is the connection length
          gravity_term = density_kg_ave * gravity * grid%dh(iup)
          ! No capillary pressure yet.
          delta_pressure = well_up%pl(iup) - well_dn%pl(idn) + &
                           gravity_term
          up_scale = 0.d0
          dn_scale = 0.d0

          ! This only applies if mobility isn't constant
          !upwind = delta_pressure > 0.d0

          !if (upwind) then
          !  up_scale = 1.d0
          !  mobility = well%liq%mobility(iup)
          !else
          !  dn_scale = 1.d0
          !  mobility = well%liq%mobility(idn)
          !endif

          ! For constant mobility
          if (well_up%liq%mobility > floweps ) then
          ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
          !                    dP[Pa]]
            v_darcy = perm_ave_over_dist(1) * well_up%liq%mobility * &
                      delta_pressure
            density_ave_kmol = density_kg_ave * fmw_comp(ONE_INTEGER)
          ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
            q = v_darcy * 5.d-1*(well_up%area(iup) + &
                                 well_dn%area(idn))
          ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
          !                             density_ave[kmol phase/m^3 phase]        
            tot_mole_flux = q*density_ave_kmol
          ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
          !                                 xmol[kmol comp/kmol phase]
            Res(1) = Res(1) + tot_mole_flux

          endif
        endif

        ! Gas flux
        ! Mobility may need to be non-constant?
        if (well_up%gas%mobility > eps) then

          density_kg_ave = 0.5d0*(well_up%gas%rho(iup)+well_dn%gas%rho(idn))
          ! Assuming the well is always vertical and gravity is in the
          ! (-) direction
          gravity_term = density_kg_ave * gravity * grid%dh(iup)
          ! No capillary pressure yet.
          delta_pressure = well_up%pg(iup) - well_dn%pg(idn) + &
                           gravity_term

          ! This only applies if mobility isn't constant
          !up_scale = 0.d0
          !dn_scale = 0.d0

          !upwind = delta_pressure > 0.d0

          !if (upwind) then
          !  up_scale = 1.d0
          !  mobility = well%gas%mobility(iup)
          !else
          !  dn_scale = 1.d0
          !  mobility = well%gas%mobility(idn)
          !endif

          if (well_up%gas%mobility > floweps ) then
          ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
          !                    dP[Pa]]
            v_darcy = perm_ave_over_dist(1) * well_up%gas%mobility * &
                      delta_pressure
            density_ave_kmol = density_kg_ave * fmw_comp(TWO_INTEGER)
          ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
            q = v_darcy * 5.d-1*(well_up%area(iup) + &
                                 well_dn%area(idn))
          ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
          !                             density_ave[kmol phase/m^3 phase]        
            tot_mole_flux = q*density_ave_kmol
          ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
          !                                 xmol[kmol comp/kmol phase]
            Res(2) = Res(2) + tot_mole_flux
          endif
        endif
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
  PetscInt :: iup, idn, idof, irow
  PetscReal :: Jup(pm_well%nphase,pm_well%nphase), &
               Jdn(pm_well%nphase,pm_well%nphase)

  PetscReal :: res_up(pm_well%nphase),res_dn(pm_well%nphase), &
               res_pert(pm_well%nphase)

  call PMWellFlux(pm_well,pm_well%well,pm_well%well,iup,idn,res_up)

  res_dn = res_up

  ! upgradient derivatives
  do idof = 1,pm_well%nphase
    call PMWellFlux(pm_well,pm_well%well_pert(idof),pm_well%well,iup,idn, &
                    res_pert)
    do irow = 1, pm_well%nphase
      Jup(irow,idof) = (res_pert(irow)-res_up(irow)) / &
                       pm_well%pert(iup,idof)
    enddo !irow
  enddo

  ! downgradient derivatives
  do idof = 1,pm_well%nphase
    call PMWellFlux(pm_well,pm_well%well,pm_well%well_pert(idof),iup,idn, &
                    res_pert)
    do irow = 1, pm_well%nphase
      Jdn(irow,idof) = (res_pert(irow)-res_dn(irow)) / &
                       pm_well%pert(idn,idof)
    enddo !irow
  enddo

end subroutine PMWellFluxDerivative

! ************************************************************************** !

subroutine PMWellBCFlux(pm_well,well,Res)
  !
  ! Computes the boundary flux terms for the residual based on the 
  ! chosen well model.
  !
  ! Author: Michael Nole
  ! Date: 12/24/2021
  !

  implicit none

  type(pm_well_type) :: pm_well
  type(well_type) :: well
  PetscReal :: Res(2*pm_well%nphase)

  type(well_grid_type), pointer :: grid
  type(well_reservoir_type), pointer :: reservoir

  PetscInt :: i, iup, idn, ghosted_id
  PetscReal :: pres_up, pres_dn

  !MAN: clean these up
  PetscReal :: perm_ave_over_dist, density_kg_ave
  PetscReal :: perm_up, perm_dn, dist_up, dist_dn
  PetscReal :: gravity_term, delta_pressure
  PetscReal :: density_ave, density_ave_kmol, tot_mole_flux
  PetscReal :: boundary_pressure, viscosity
  PetscReal :: v_darcy,q
  PetscBool :: upwind
  PetscInt :: idof, irow, itop

  PetscReal, parameter :: eps = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

  grid => pm_well%grid
  reservoir => pm_well%reservoir

  Res(:) = 0.d0

  select case(well%well_model_type)
    case('CONSTANT_RATE', 'CONSTANT_PRESSURE','CONSTANT_PRESSURE_HYDROSTATIC')
      ! No nonlinear solve needed.
    case('DARCY')
      ! This is good for either: single-phase liquid, or two-phase liquid/gas.
      ! Vertical well, no Klinkenberg, no capillary pressure, constant mobility.
      itop = pm_well%grid%nsegments
      if (pm_well%flow_soln%bh_p) then
        !Dirichlet pressure and saturation at the bottom
        perm_ave_over_dist = well%permeability(1) / (grid%dh(1)/2.d0)
        boundary_pressure = well%bh_p
        gravity_term = well%liq%rho(1) * gravity * &
                       grid%dh(1)/2.d0
        delta_pressure = boundary_pressure - well%pl(1) + gravity_term

        ! This only applies if mobility isn't constant
        !dn_scale = 0.d0
        !upwind = delta_pressure > 0.d0
        !if (upwind) then
        !  rel_perm = wippflo_auxvar_up%kr(iphase)
        !else
        !  dn_scale = 1.d0
        !  rel_perm = wippflo_auxvar_dn%kr(iphase)
        !endif

        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy = perm_ave_over_dist * well%liq%mobility * delta_pressure
        density_ave = well%liq%rho(1) / fmw_comp(ONE_INTEGER)
        q = v_darcy * well%area(1)
        tot_mole_flux = q * density_ave
        Res(1) = Res(1) + tot_mole_flux

        v_darcy = perm_ave_over_dist * well%gas%mobility * delta_pressure
        density_ave = well%gas%rho(1) / fmw_comp(TWO_INTEGER)
        q = v_darcy * well%area(1)
        tot_mole_flux = q * density_ave
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
        Res(1) = Res(1) + tot_mole_flux
 
        v_darcy = well%bh_qg
        if (v_darcy > 0.d0) then
          density_ave = reservoir%rho_g(1) / fmw_comp(TWO_INTEGER)
        else
          density_ave = well%gas%rho(1) / fmw_comp(TWO_INTEGER)
        endif
        q = v_darcy * well%area(1)
        tot_mole_flux = q * density_ave
        Res(2) = Res(2) + tot_mole_flux
      else
        ! this should not happen once error messaging is updated
      endif

      if (pm_well%flow_soln%th_p) then
        !Dirichlet pressure and saturation at the top
        perm_ave_over_dist = well%permeability(itop) / &
                             (grid%dh(itop)/2.d0)
        boundary_pressure = well%th_p
        gravity_term = well%liq%rho(itop) * gravity * &
                       grid%dh(itop)/2.d0
        delta_pressure = boundary_pressure - well%pl(itop) + gravity_term

        ! This only applies if mobility isn't constant
        !dn_scale = 0.d0
        !upwind = delta_pressure > 0.d0
        !if (upwind) then
        !  rel_perm = wippflo_auxvar_up%kr(iphase)
        !else
        !  dn_scale = 1.d0
        !  rel_perm = wippflo_auxvar_dn%kr(iphase)
        !endif

        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy = perm_ave_over_dist * well%liq%mobility * delta_pressure
        density_ave = well%liq%rho(itop) / fmw_comp(ONE_INTEGER)
        q = v_darcy * well%area(itop)
        tot_mole_flux = q * density_ave
        Res(3) = Res(3) + tot_mole_flux

        v_darcy = perm_ave_over_dist * well%gas%mobility * delta_pressure
        density_ave = well%gas%rho(itop) / fmw_comp(TWO_INTEGER)
        q = v_darcy * well%area(itop)
        tot_mole_flux = q * density_ave
        Res(4) = Res(4) + tot_mole_flux
      else
        !Neumann flux at the top
        v_darcy = well%th_ql
        if (v_darcy > 0.d0) then
          density_ave = reservoir%rho_l(itop) / fmw_comp(ONE_INTEGER)
        else
          density_ave = well%liq%rho(itop) / fmw_comp(ONE_INTEGER)
        endif
        q = v_darcy * well%area(itop)
        tot_mole_flux = q * density_ave
        Res(3) = Res(3) + tot_mole_flux

        v_darcy = well%th_qg
        if (v_darcy > 0.d0) then
          density_ave = reservoir%rho_g(itop) / fmw_comp(TWO_INTEGER)
        else
          density_ave = well%gas%rho(itop) / fmw_comp(TWO_INTEGER)
        endif
        q = v_darcy * well%area(itop)
        tot_mole_flux = q * density_ave
        Res(4) = Res(4) + tot_mole_flux
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
  PetscInt :: idn
  PetscReal :: Jtop(pm_well%flow_soln%ndof,pm_well%flow_soln%ndof), &
               Jbtm(pm_well%flow_soln%ndof,pm_well%flow_soln%ndof)

  PetscInt :: idof, irow
  PetscReal :: res(2*pm_well%flow_soln%ndof),res_pert(2*pm_well%flow_soln%ndof)

  call PMWellBCFlux(pm_well,pm_well%well,res)
  
  Jtop = 0.d0
  Jbtm = 0.d0

  ! downgradient derivatives
  do idof = 1,pm_well%nphase
    call PMWellBCFlux(pm_well,pm_well%well_pert(idof),res_pert)
    do irow = 1, pm_well%nphase
      Jbtm(irow,idof) = (res_pert(irow)-res(irow)) / &
                        pm_well%pert(1,idof)
    enddo
    do irow = 1, pm_well%nphase
      Jtop(irow,idof) = (res_pert(irow + pm_well%flow_soln%ndof)- &
                         res(irow + pm_well%flow_soln%ndof)) / &
                         pm_well%pert(pm_well%grid%nsegments,idof)
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

  use Option_module

  implicit none

  type(pm_well_type) :: pm_well

  PetscReal :: x(pm_well%grid%nsegments,pm_well%nphase), &
               x_pert(pm_well%grid%nsegments,pm_well%nphase), &
               pert(pm_well%grid%nsegments,pm_well%nphase)
  PetscBool :: pert_indices(pm_well%grid%nsegments)
  PetscInt :: idof,i

  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-10

  ! I don't like how pressure and saturation are attributes of 
  ! different objects
  x(:,ONE_INTEGER) = pm_well%well%pl
  x(:,TWO_INTEGER) = pm_well%well%gas%s

  pert(:,ONE_INTEGER) = perturbation_tolerance*x(:,ONE_INTEGER) + &
                        min_perturbation
  do i = 1,pm_well%grid%nsegments
    if (x(i,TWO_INTEGER) > 0.5d0) then
      pert(i,TWO_INTEGER) = -1.d0 * perturbation_tolerance
    else    
      pert(i,TWO_INTEGER) = perturbation_tolerance
    endif
  enddo

  pm_well%well_pert(ONE_INTEGER)%pl = x(:,ONE_INTEGER) + pert(:,ONE_INTEGER)
  pm_well%well_pert(TWO_INTEGER)%gas%s = x(:,TWO_INTEGER) + pert(:,TWO_INTEGER)

  call PMWellUpdateProperties(pm_well%well_pert(ONE_INTEGER))
  call PMWellUpdateProperties(pm_well%well_pert(TWO_INTEGER))

  pm_well%pert = pert

end subroutine PMWellPerturb

! ************************************************************************** !

subroutine PMWellUpdateProperties(well)
  !
  ! Updates well object properties.
  !
  ! Author: Michael Nole
  ! Date: 01/06/2022
  !

  use Option_module
  use EOS_Water_module

  implicit none

  type(well_type) :: well

  PetscInt :: i
  PetscReal :: t,dw,dwmol,dwp,dwt
  PetscErrorCode :: ierr

  ! Saturations
  do i = 1,size(well%gas%s)
    well%gas%s(i) = max(well%gas%s(i),0.d0)
    well%gas%s(i) = min(well%gas%s(i),1.d0)
  enddo
  well%liq%s = 1.d0 - well%gas%s

  !no capillarity right now
  well%pg = well%pl

  !Densities
  do i = 1,size(well%pl)
    call EOSWaterDensityBRAGFLO(t,well%pl(i),PETSC_FALSE, &
                                dw,dwmol,dwp,dwt,ierr)
    well%liq%rho(i) = dw
  enddo

end subroutine PMWellUpdateProperties

! ************************************************************************** !

function PMWellOutputFilename(option)
  ! 
  ! Generates a filename for wellbore model output
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 12/16/2021

  use Option_module

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
  type(grid_type), pointer :: grid
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: cell_string
  PetscBool :: exist
  PetscInt :: fid
  PetscInt :: icolumn
  PetscInt :: k
  
  output_option => this%realization%output_option
  grid => this%realization%patch%grid

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

  do k = 1,this%grid%nsegments
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
    variable_string = 'P'
    units_string = 'Pa' 
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'S-liq'
    units_string = '-' 
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'S-gas'
    units_string = '-' 
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
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

  use Option_module
  use Output_Aux_module
  use Global_Aux_module
  use Grid_module

  implicit none
  
  class(pm_well_type) :: this
  
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(grid_type), pointer :: grid
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: k
  
100 format(100es18.8)
101 format(1I6.1)

  option => this%realization%option
  output_option => this%realization%output_option
  grid => this%realization%patch%grid
  
  fid = 555
  filename = PMWellOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")

  ! this time is set at the end of the wellbore step????
  write(fid,100,advance="no") option%time / output_option%tconv

  do k = 1,this%grid%nsegments
    write(fid,101,advance="no") k
    write(fid,100,advance="no") this%grid%h(k)%x, &
                                this%grid%h(k)%y, &
                                this%grid%h(k)%z, &
                                this%well%pl(k), &
                                this%well%liq%s(k), &
                                this%well%gas%s(k) 
  enddo
  
  close(fid)
  
end subroutine PMWellOutput

! ************************************************************************** !

subroutine PMWellDestroy(this)
  ! 
  ! Destroys the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !
  use Utility_module, only : DeallocateArray
  use Option_module

  implicit none
  
  class(pm_well_type) :: this
    
  call PMBaseDestroy(this)

  call DeallocateArray(this%grid%h_local_id)
  call DeallocateArray(this%grid%h_ghosted_id)
  call DeallocateArray(this%grid%dh)
  nullify(this%grid)

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
  
end subroutine PMWellDestroy

! ************************************************************************** !

end module PM_Well_class
