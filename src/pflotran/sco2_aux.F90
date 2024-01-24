module SCO2_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none

  private

  ! Physical Properties
  ! MAN: change FMWNACL and density to generic salt?
  PetscReal, parameter, public :: fmw_comp(3) = [FMWH2O,FMWCO2,FMWNACL]
  PetscReal, parameter, public :: SALT_DENSITY_KG = 2.170D3 !kg/m^3
  ! PetscReal, parameter, public :: SALT_THERMAL_CONDUCTIVITY = 6.d0 !W/m-K
  PetscReal, parameter :: CO2_REFERENCE_SURFACE_TENSION = 0.072d0 ! N/m
  PetscReal, parameter :: SALT_REFERENCE_TEMPERATURE = 293.15d0
  PetscReal, parameter :: LIQUID_REFERENCE_VISCOSITY = 1.01764892595942d-3
  PetscReal, parameter, public :: LIQUID_REFERENCE_DENSITY = 998.32142721500441
  PetscReal, parameter, public :: SCO2_REFERENCE_PRESSURE = 101325.d0

  ! Solution Control
  PetscReal, public :: sco2_phase_chng_epsilon = 1.d-6
  PetscBool, public :: sco2_chk_max_dpl_liq_state_only = PETSC_FALSE
  PetscBool, public :: sco2_restrict_state_chng = PETSC_FALSE
  PetscReal, public :: sco2_window_epsilon = 0.d0
  PetscInt, public :: sco2_diffusion_model = ONE_INTEGER
  PetscBool, public :: sco2_spycher_simple = PETSC_FALSE !PETSC_TRUE
  PetscBool, public :: sco2_central_diff_jacobian = PETSC_FALSE
  PetscBool, public :: sco2_harmonic_diff_density = PETSC_TRUE
  ! PetscBool, public :: sco2_high_temp_ts_cut = PETSC_FALSE
  PetscBool, public :: sco2_allow_state_change = PETSC_TRUE
  PetscBool, public :: sco2_state_changed = PETSC_FALSE
  PetscBool, public :: sco2_force_iteration = PETSC_FALSE
  PetscBool, public :: sco2_newtontrdc_hold_inner = PETSC_FALSE
  PetscInt, public :: sco2_newton_iteration_number = 0
  PetscInt, public :: sco2_sub_newton_iter_num = 0
  PetscInt, public :: sco2_newtontrdc_prev_iter_num = 0
  PetscInt, public :: sco2_max_it_before_damping = UNINITIALIZED_INTEGER
  PetscReal, public :: sco2_max_pressure_change = 5.d4
  PetscReal, public :: sco2_isothermal_temperature = 25.d0


  ! Output Control
  PetscBool, public :: sco2_print_state_transition = PETSC_TRUE

  ! Debugging
  PetscInt, public :: sco2_debug_cell_id = UNINITIALIZED_INTEGER
  PetscInt, public :: sco2_ni_count
  PetscInt, public :: sco2_ts_cut_count
  PetscInt, public :: sco2_ts_count

  ! Thermodynamic States
  PetscInt, parameter, public :: SCO2_NULL_STATE = 0
  PetscInt, parameter, public :: SCO2_LIQUID_STATE = 1
  PetscInt, parameter, public :: SCO2_GAS_STATE = 2
  PetscInt, parameter, public :: SCO2_TRAPPED_GAS_STATE = 3
  PetscInt, parameter, public :: SCO2_LIQUID_GAS_STATE = 4
  PetscInt, parameter, public :: SCO2_MAX_STATE = 4

  ! For over/under-determined flow conditions
  PetscInt, parameter, public :: SCO2_ANY_STATE = 8
  PetscInt, parameter, public :: SCO2_MULTI_STATE = 9

  PetscInt, parameter, public :: PREV_TS = 1
  PetscInt, parameter, public :: PREV_IT = 2
  
  ! Indexing the primary variables
  ! Gas Pressure is always DOF 1
  PetscInt, parameter, public :: SCO2_GAS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: SCO2_LIQUID_PRESSURE_DOF = 1

  ! DOF 2 in Gas State
  PetscInt, parameter, public :: SCO2_CO2_PRESSURE_DOF = 2
  ! DOF 2 in Liquid State:
  PetscInt, parameter, public :: SCO2_CO2_MASS_FRAC_DOF = 2
  ! DOF 2 in Gas State:
  ! PetscInt, parameter, public :: SCO2_LIQUID_SATURATION_DOF = 2
  ! DOF 2 in Trapped Gas State:
  PetscInt, parameter, public :: SCO2_GAS_SATURATION_DOF = 2
  PetscInt, parameter, public :: SCO2_TWO_PHASE_GAS_PRES_DOF = 2
  ! Without energy, salt mass is DOF 3
  PetscInt, parameter, public :: SCO2_SALT_MASS_FRAC_DOF = 3

  ! Auxvar mapping from flow condition couplers
  PetscInt, parameter, public :: SCO2_LIQUID_PRESSURE_INDEX = 2
  PetscInt, parameter, public :: SCO2_GAS_PRESSURE_INDEX = 3
  PetscInt, parameter, public :: SCO2_AIR_PRESSURE_INDEX = 4
  PetscInt, parameter, public :: SCO2_MOLE_FRACTION_INDEX = 5
  ! PetscInt, parameter, public :: SCO2_TEMPERATURE_INDEX = 6
  PetscInt, parameter, public :: SCO2_GAS_SATURATION_INDEX = 7
  PetscInt, parameter, public :: SCO2_LIQUID_FLUX_INDEX = 8
  PetscInt, parameter, public :: SCO2_GAS_FLUX_INDEX = 9
  ! PetscInt, parameter, public :: SCO2_ENERGY_FLUX_INDEX = 10
  PetscInt, parameter, public :: SCO2_LIQUID_CONDUCTANCE_INDEX = 11
  PetscInt, parameter, public :: SCO2_GAS_CONDUCTANCE_INDEX = 12
  PetscInt, parameter, public :: SCO2_GAS_WATER_MOL_FRAC_INDEX = 13
  PetscInt, parameter, public :: SCO2_SALT_INDEX = 14
  PetscInt, parameter, public :: SCO2_MAX_INDEX = 16

  
  ! Temperature is always DOF 4
  ! PetscInt, parameter, public :: SCO2_TEMPERATURE_DOF = 4

  ! Indexing the equations/residuals
  PetscInt, parameter, public :: SCO2_WATER_EQUATION_INDEX = 1
  PetscInt, parameter, public :: SCO2_CO2_EQUATION_INDEX = 2
  ! PetscInt, parameter, public :: SCO2_ENERGY_EQUATION_INDEX = 3
  ! PetscInt, parameter, public :: SCO2_SALT_EQUATION_INDEX = 4
  PetscInt, parameter, public :: SCO2_SALT_EQUATION_INDEX = 3

  ! Update flags
  PetscInt, parameter, public :: SCO2_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: SCO2_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: SCO2_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: SCO2_UPDATE_FOR_BOUNDARY = 2
  PetscInt, parameter, public :: SCO2_UPDATE_FOR_SS = 3

  ! Physics Options
  ! PetscBool, public :: sco2_isothermal = PETSC_FALSE
  ! MAN: might want to move these elsewhere
  PetscBool, public :: sco2_update_permeability = PETSC_FALSE
  PetscReal, public :: permeability_func_porosity_exp = 1.d0
  PetscInt, public :: permeability_reduction_model = TWO_INTEGER

  ! DOF map
  PetscInt, public, pointer :: dof_to_primary_variable(:,:)

  type, public :: sco2_auxvar_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    PetscBool :: istatechng
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal :: max_gas_sat ! Store max gas sat separately
    PetscReal, pointer :: den(:) ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: den_kg(:) ! (iphase) kg/m^3 phase
    PetscReal :: temp
    PetscReal :: m_salt(2)     ! (kg NaCl / kg (brine + precipitate)), kg NaCl
    PetscReal :: sl_min  ! min liquid saturation for hysteresis
    PetscReal :: sg_trapped !keep track of trapped gas at beginning of timestep
    PetscReal, pointer :: xmass(:,:) ! (icomp,iphase)
    PetscReal, pointer :: xmol(:,:)
    PetscReal, pointer :: effective_diffusion_coeff(:,:) ! (icomp,iphase)
    PetscReal, pointer :: dispersivity(:,:) ! (icomp,iphase)
    PetscReal, pointer :: H(:) ! MJ/kg
    PetscReal, pointer :: U(:) ! MJ/kg
    PetscReal, pointer :: kr(:) ! (iphase)
    PetscReal, pointer :: visc(:) ! Pa-s
    PetscReal, pointer :: mobility(:) ! relative perm / kinematic viscosity
    PetscReal :: effective_porosity ! factors in compressibility
    PetscReal :: effective_permeability ! scaling factor on permeability
    PetscReal, pointer :: tortuosity(:) ! (iphase)
    PetscReal :: perm_base
    PetscReal :: pert
  end type sco2_auxvar_type

  type, public :: sco2_parameter_type
    PetscReal, pointer :: diffusion_coefficient(:,:) ! (icomp,iphase)
    PetscReal :: newton_inf_scaled_res_tol
    PetscBool :: check_post_converged
  end type sco2_parameter_type

  type, public :: sco2_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(sco2_parameter_type), pointer :: sco2_parameter
    type(sco2_auxvar_type), pointer :: auxvars(:,:)
    type(sco2_auxvar_type), pointer :: auxvars_bc(:)
    type(sco2_auxvar_type), pointer :: auxvars_ss(:,:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type sco2_type

  interface SCO2AuxVarDestroy
    module procedure SCO2AuxVarSingleDestroy
    module procedure SCO2AuxVarArray1Destroy
    module procedure SCO2AuxVarArray2Destroy
  end interface SCO2AuxVarDestroy

  interface SCO2OutputAuxVars
    module procedure SCO2OutputAuxVars1
    module procedure SCO2OutputAuxVars2
  end interface SCO2OutputAuxVars

  public :: SCO2AuxCreate, &
            SCO2AuxVarInit, &
            SCO2AuxVarCopy, &
            SCO2AuxVarPerturb, &
            SCO2AuxVarUpdateState, &
            SCO2AuxDestroy, &
            SCO2AuxVarCompute, &
            SCO2OutputAuxVars, &
            SCO2AuxVarDestroy, &
            SCO2AuxVarStrip, &
            SCO2ComputeSaltSolubility, &
            SCO2BrineSaturationPressure, &
            SCO2VaporPressureBrine, &
            SCO2BrineDensity, &
            SCO2Henry, &
            SCO2ComputeSaltDensity, &
            SCO2Equilibrate
            

contains

! ************************************************************************** !

function SCO2AuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 11/22/23
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option

  type(sco2_type), pointer :: SCO2AuxCreate

  type(sco2_type), pointer :: aux


  allocate(dof_to_primary_variable(option%nflowdof,SCO2_MAX_STATE))
  
  dof_to_primary_variable(1:option%nflowdof,1:SCO2_MAX_STATE) = &
               ! Liquid State
      reshape([SCO2_LIQUID_PRESSURE_DOF, SCO2_CO2_MASS_FRAC_DOF, & !L
               SCO2_SALT_MASS_FRAC_DOF, &
               ! Gas State
               SCO2_GAS_PRESSURE_DOF, SCO2_CO2_PRESSURE_DOF, &   !G
               SCO2_SALT_MASS_FRAC_DOF, &
               ! Trapped Gas State
               SCO2_LIQUID_PRESSURE_DOF, SCO2_GAS_SATURATION_DOF, &   !TG
               SCO2_SALT_MASS_FRAC_DOF, &
               ! Liquid-Gas State
               SCO2_LIQUID_PRESSURE_DOF, SCO2_TWO_PHASE_GAS_PRES_DOF, & !LG
               SCO2_SALT_MASS_FRAC_DOF],&
               shape(dof_to_primary_variable))
  !With energy             
  ! dof_to_primary_variable(1:option%nflowdof,1:SCO2_MAX_STATE) = &
  !              ! Liquid State
  !     reshape([SCO2_LIQUID_PRESSURE_DOF, SCO2_CO2_MASS_FRAC_DOF, & !L
  !              SCO2_TEMPERATURE_DOF, SCO2_SALT_MASS_FRAC_DOF, &
  !              ! Gas State
  !              SCO2_GAS_PRESSURE_DOF, SCO2_CO2_PRESSURE_DOF, &   !G
  !              SCO2_TEMPERATURE_DOF, SCO2_SALT_MASS_FRAC_DOF, &
  !              ! Trapped Gas State
  !              SCO2_LIQUID_PRESSURE_DOF, SCO2_GAS_SATURATION_DOF, &   !TG
  !              SCO2_TEMPERATURE_DOF, SCO2_SALT_MASS_FRAC_DOF, &
  !              ! Liquid-Gas State
  !              SCO2_LIQUID_PRESSURE_DOF, SCO2_TWO_PHASE_GAS_PRES_DOF, & !LG
  !              SCO2_TEMPERATURE_DOF, SCO2_SALT_MASS_FRAC_DOF],&
  !              shape(dof_to_primary_variable))

  allocate(aux)
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%matrix_zeroing)

  allocate(aux%sco2_parameter)
  allocate(aux%sco2_parameter%diffusion_coefficient(option%nflowspec, &
                                                    option%nflowdof))

  aux%sco2_parameter%diffusion_coefficient(:,LIQUID_PHASE) = &
                                                   UNINITIALIZED_DOUBLE
  aux%sco2_parameter%diffusion_coefficient(:,GAS_PHASE) = 2.13d-05
  
  aux%sco2_parameter%newton_inf_scaled_res_tol = 1.d-50
  aux%sco2_parameter%check_post_converged = PETSC_FALSE

  SCO2AuxCreate => aux

end function SCO2AuxCreate

! ************************************************************************** !

subroutine SCO2AuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 11/22/23
  !

  use Option_module

  implicit none

  type(sco2_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%istate_store = SCO2_NULL_STATE
  auxvar%istatechng = PETSC_FALSE
  auxvar%temp = sco2_isothermal_temperature
  auxvar%effective_porosity = 0.d0
  auxvar%effective_permeability = 0.d0
  auxvar%pert = 0.d0
  auxvar%perm_base = -999.9d0
  auxvar%m_salt = 0.d0
  auxvar%sl_min = 1.d0
  auxvar%sg_trapped = 0.d0

  allocate(auxvar%pres(option%nphase+FOUR_INTEGER))
  auxvar%pres = 0.d0
  allocate(auxvar%sat(option%trapped_gas_phase))
  auxvar%sat = 0.d0
  auxvar%max_gas_sat = 0.d0
  allocate(auxvar%den_kg(3+option%nphase)) ! Pure component and mixture
  auxvar%den_kg = 0.d0
  allocate(auxvar%den(3+option%nphase)) ! Pure component and mixture
  auxvar%den = 0.d0
  allocate(auxvar%xmass(option%nflowspec,option%nphase))
  auxvar%xmass = 0.d0
  allocate(auxvar%xmol(option%nflowspec,option%nphase))
  auxvar%xmol = 0.d0  
  allocate(auxvar%effective_diffusion_coeff(option%nflowspec,option%nphase))
  auxvar%effective_diffusion_coeff = 0.d0
  allocate(auxvar%dispersivity(option%nflowspec,option%nphase))
  auxvar%dispersivity = 0.d0
  allocate(auxvar%H(3+option%nphase)) ! Pure component and mixture
  auxvar%H = 0.d0
  allocate(auxvar%U(3+option%nphase)) ! Pure component and mixture
  auxvar%U = 0.d0
  allocate(auxvar%mobility(option%nphase))
  auxvar%mobility = 0.d0
  allocate(auxvar%kr(option%nphase))
  auxvar%kr = 0.d0
  allocate(auxvar%visc(option%nphase))
  auxvar%visc = 0.d0
  allocate(auxvar%tortuosity(option%nphase))
  auxvar%tortuosity = 1.d0
  

end subroutine SCO2AuxVarInit

! ************************************************************************** !

subroutine SCO2AuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Michael Nole
  ! Date: 11/22/23
  !

  use Option_module

  implicit none

  type(sco2_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%istate_store = auxvar%istate_store
  auxvar2%istatechng = auxvar%istatechng
  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%xmass = auxvar%xmass
  auxvar2%xmol = auxvar%xmol
  auxvar2%m_salt = auxvar%m_salt
  auxvar2%sl_min = auxvar%sl_min
  auxvar2%sg_trapped = auxvar%sg_trapped
  auxvar2%effective_diffusion_coeff = auxvar%effective_diffusion_coeff
  auxvar2%dispersivity = auxvar%dispersivity
  auxvar2%H = auxvar%H
  auxvar2%U = auxvar%U
  auxvar2%mobility = auxvar%mobility
  auxvar2%kr = auxvar%kr
  auxvar2%visc = auxvar%visc
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%effective_permeability = auxvar%effective_permeability
  auxvar2%tortuosity = auxvar%tortuosity
  auxvar2%pert = auxvar%pert

end subroutine SCO2AuxVarCopy

! ************************************************************************** !

subroutine SCO2AuxVarPerturb(sco2_auxvar, global_auxvar, material_auxvar, &
                             characteristic_curves, sco2_parameter, &
                             natural_id, option)
  !
  ! Computes perturbations in the primary variables for SCO2 mode
  !
  ! Author: Michael Nole
  ! Date: 12/05/23
  !

  use Option_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  type(sco2_auxvar_type) :: sco2_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type), pointer :: characteristic_curves
  type(sco2_parameter_type), pointer :: sco2_parameter
  PetscInt :: natural_id
  type(option_type) :: option

  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-15

  PetscReal :: x(option%nflowdof), x_pert_plus(option%nflowdof), &
               pert(option%nflowdof), x_pert_minus(option%nflowdof)
  PetscReal :: xco2g, xwg, xco2l, xsl, xwl, xmolco2g, xmolwg, xmolco2l, &
               xmolsl, xmolwl, salt_mass
  PetscReal :: sigma, beta_gl
  PetscReal :: Pv, Psat, Prvap, Pco2
  PetscReal :: dpl, dpg, dpco2, dxco2, dxs, dsg ! ,dt
  PetscReal :: cell_pressure, sgt_max
  PetscInt :: idof

  ! Phase ID's
  PetscInt :: lid, gid, pid, pwid, pbid, spid, tgid
  ! Component ID's
  PetscInt :: wid, co2_id, co2_pressure_id, sid, pgid
  ! Other ID's
  PetscInt :: cpid, vpid, rvpid
  
  lid = option%liquid_phase
  gid = option%gas_phase
  pid = option%precipitate_phase
  pwid = option%pure_water_phase
  pbid = option%pure_brine_phase
  tgid = option%trapped_gas_phase
  pgid = option%trapped_gas_phase
  spid = option%saturation_pressure_id

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id

  cpid = option%capillary_pressure_id
  co2_pressure_id = option%co2_pressure_id
  vpid = option%vapor_pressure_id
  rvpid = option%reduced_vapor_pressure_id
        
  call SCO2ComputeSaltSolubility(sco2_auxvar(ZERO_INTEGER)%temp,xsl)
  dxs = 1.0d-5 * xsl
  !MAN: need to make sure total salt mass is updated in AuxVarCompute
  salt_mass = sco2_auxvar(ZERO_INTEGER)%m_salt(ONE_INTEGER)
  xsl = min(salt_mass,xsl)

  ! dt = -1.d0 * perturbation_tolerance * (sco2_auxvar(ZERO_INTEGER)%temp + &
  !      min_perturbation)

  call SCO2ComputeSurfaceTension(sco2_auxvar(ZERO_INTEGER)%temp, xsl, sigma)
  beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma

  sgt_max = characteristic_curves%saturation_function%Sgt_max

  select case(global_auxvar%istate)
    case(SCO2_LIQUID_STATE)
      dpl = max(1.d-1, 1.d-7 * sco2_auxvar(ZERO_INTEGER)%pres(lid))
      cell_pressure = sco2_auxvar(ZERO_INTEGER)%pres(lid)
      ! cell_pressure = min(cell_pressure, 1.d8)          
      call SCO2BrineSaturationPressure(sco2_auxvar(ZERO_INTEGER)%temp, xsl, &
                                       Psat)
      Prvap = Psat
      call SCO2Equilibrate(sco2_auxvar(ZERO_INTEGER)%temp,cell_pressure, &
                           Pco2, Pv, Psat, Prvap, &
                           xco2g, xwg, xco2l, xsl, xwl, &
                           xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl, option)
      if (sco2_auxvar(ZERO_INTEGER)%xmass(co2_id,lid) > (1.d-2 * xco2l)) then
        dxco2 = sign(1.d-4 * xco2l, &
                     5.d-1 * xco2l - &
                     sco2_auxvar(ZERO_INTEGER)%xmass(co2_id,lid))
      else
        dxco2 = sign(1.d-3 * xco2l, &
                     5.d-1 * xco2l - &
                     sco2_auxvar(ZERO_INTEGER)%xmass(co2_id,lid))
      endif
    case(SCO2_GAS_STATE)

      dpg = 1.d-3

      dpco2 = max(1.d-2, &
                  1.d-7 * dabs(sco2_auxvar(ZERO_INTEGER)%pres(gid) - &
                               sco2_auxvar(ZERO_INTEGER)%pres(lid)))

    case(SCO2_TRAPPED_GAS_STATE)

      dpl = max(1.d-1, 1.d-6 * sco2_auxvar(ZERO_INTEGER)%pres(lid))

      dsg = sign(1.d-7, 5.d-1 * sgt_max - sco2_auxvar(ZERO_INTEGER)%sat(gid))

    case(SCO2_LIQUID_GAS_STATE)
      ! Perturb into the state
      dpl = max(1.d-2, &
                1.d-6 * dabs(sco2_auxvar(ZERO_INTEGER)%pres(gid) - &
                             sco2_auxvar(ZERO_INTEGER)%pres(lid)))
      dpl = sign(dpl, (5.d-1 - sco2_auxvar(ZERO_INTEGER)%sat(lid)))
      dpg = -1.d0 * dpl
  end select

  dpg = dpg 
  dpl = dpl
  dpco2 = dpco2
  dsg = dsg 
  dxs = dxs 
  dxco2 = dxco2

  select case(global_auxvar%istate)
    case(SCO2_LIQUID_STATE)

      x(SCO2_LIQUID_PRESSURE_DOF) = sco2_auxvar(ZERO_INTEGER)%pres(lid)
      x(SCO2_CO2_MASS_FRAC_DOF) = sco2_auxvar(ZERO_INTEGER)%xmass(co2_id,lid)
      ! x(SCO2_TEMPERATURE_DOF) = sco2_auxvar(ZERO_INTEGER)%temp
      x(SCO2_SALT_MASS_FRAC_DOF) = sco2_auxvar(ZERO_INTEGER)%m_salt(1)

      pert(SCO2_LIQUID_PRESSURE_DOF) = dpl
      pert(SCO2_CO2_MASS_FRAC_DOF) = dxco2
      ! pert(SCO2_TEMPERATURE_DOF) = dt
      pert(SCO2_SALT_MASS_FRAC_DOF) = dxs

    case(SCO2_GAS_STATE)

      x(SCO2_GAS_PRESSURE_DOF) = sco2_auxvar(ZERO_INTEGER)%pres(gid)
      x(SCO2_CO2_PRESSURE_DOF) = sco2_auxvar(ZERO_INTEGER)%pres(co2_pressure_id) 
      ! x(SCO2_TEMPERATURE_DOF) = sco2_auxvar(ZERO_INTEGER)%temp
      x(SCO2_SALT_MASS_FRAC_DOF) = sco2_auxvar(ZERO_INTEGER)%m_salt(2)

      pert(SCO2_GAS_PRESSURE_DOF) = dpg
      pert(SCO2_CO2_PRESSURE_DOF) = -dpco2
      ! pert(SCO2_TEMPERATURE_DOF) = -1.d0 * dt
      pert(SCO2_SALT_MASS_FRAC_DOF) = dxs

    case(SCO2_TRAPPED_GAS_STATE)

      x(SCO2_LIQUID_PRESSURE_DOF) = sco2_auxvar(ZERO_INTEGER)%pres(lid)
      x(SCO2_GAS_SATURATION_DOF) = sco2_auxvar(ZERO_INTEGER)%sat(gid)
      ! x(SCO2_TEMPERATURE_DOF) = sco2_auxvar(ZERO_INTEGER)%temp
      x(SCO2_SALT_MASS_FRAC_DOF) = sco2_auxvar(ZERO_INTEGER)%m_salt(1)

      pert(SCO2_LIQUID_PRESSURE_DOF) = dpl
      pert(SCO2_GAS_SATURATION_DOF) = dsg
      ! pert(SCO2_TEMPERATURE_DOF) = sign(dt, dsg)
      pert(SCO2_SALT_MASS_FRAC_DOF) = dxs

    case(SCO2_LIQUID_GAS_STATE)

      x(SCO2_LIQUID_PRESSURE_DOF) = sco2_auxvar(ZERO_INTEGER)%pres(lid)
      x(SCO2_TWO_PHASE_GAS_PRES_DOF) = sco2_auxvar(ZERO_INTEGER)%pres(gid)
      ! x(SCO2_TEMPERATURE_DOF) = sco2_auxvar(ZERO_INTEGER)%temp
      x(SCO2_SALT_MASS_FRAC_DOF) = sco2_auxvar(ZERO_INTEGER)%m_salt(1)

      pert(SCO2_LIQUID_PRESSURE_DOF) = dpl
      pert(SCO2_TWO_PHASE_GAS_PRES_DOF) = dpg
      ! pert(SCO2_TEMPERATURE_DOF) = sign(dt, dpg)
      pert(SCO2_SALT_MASS_FRAC_DOF) = dxs

    case default

      write(option%io_buffer,*) global_auxvar%istate
      option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
                         ') not recognized in SCO2AuxVarPerturb.'
      call PrintErrMsgByRank(option)

  end select

  ! SCO2_UPDATE_FOR_DERIVATIVE indicates call from perturbation

  option%iflag = SCO2_UPDATE_FOR_DERIVATIVE

  do idof = 1, option%nflowdof

    if (sco2_central_diff_jacobian) then
      !pert(idof) = max(1.d-7 * x(idof),1.d-7)

      x_pert_minus = x
      x_pert_minus(idof) = x(idof) - pert(idof)
      call SCO2AuxVarCompute(x_pert_minus, &
             sco2_auxvar(idof+option%nflowdof),global_auxvar,material_auxvar, &
             characteristic_curves,sco2_parameter,natural_id,option)

      sco2_auxvar(idof+option%nflowdof)%pert = pert(idof)

    endif

    sco2_auxvar(idof)%pert = pert(idof)
    x_pert_plus = x
    x_pert_plus(idof) = x(idof) + pert(idof)

    call SCO2AuxVarCompute(x_pert_plus,sco2_auxvar(idof),global_auxvar, &
                           material_auxvar, &
                           characteristic_curves,sco2_parameter, &
                           natural_id,option)
  enddo

end subroutine SCO2AuxVarPerturb

! ************************************************************************** !

subroutine SCO2AuxVarUpdateState(x, sco2_auxvar, global_auxvar, &
                                 material_auxvar, characteristic_curves, &
                                 sco2_parameter, natural_id, option)
  !
  ! Computes state transitions and swaps primary variables for SCO2 mode
  !
  ! Author: Michael Nole
  ! Date: 12/05/23
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  
  type(sco2_auxvar_type) :: sco2_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  type(sco2_parameter_type), pointer :: sco2_parameter
  PetscInt :: natural_id
  type(option_type) :: option
  PetscReal :: x(option%nflowdof)

  PetscReal, parameter :: epsilon = 1.d-14
  PetscReal, parameter :: eps_sl = 1.d-4
  PetscReal, parameter :: peta = 1.d-1

  PetscReal :: Pc_entry, cell_pressure, sg_min
  PetscReal :: sl_temp, sgt_temp, sg_est, Slr, Sgt
  PetscReal :: Pc, Pv, Prvap, Pco2, Psat, Pg
  PetscReal :: beta_gl, sgt_max
  PetscReal :: xco2g, xwg, xco2l, xsl, xwl, xmolco2g, xmolwg, xmolco2l, &
               xmolsl, xmolwl
  PetscReal :: salt_solubility, sigma 
  PetscReal :: den_mol, den_co2, den_brine, den_liq
  PetscReal :: drho_dP, drho_dT
  PetscErrorCode :: ierr
  ! Phase ID's
  PetscInt :: lid, gid, pid, pwid, pbid, spid, tgid
  ! Component ID's
  PetscInt :: wid, co2_id, co2_pressure_id, sid
  ! Other ID's
  PetscInt :: cpid, vpid, rvpid
  PetscReal :: state_change_threshold

  PetscInt :: old_state, new_state
  character(len=MAXSTRINGLENGTH) :: state_change_string, append
  PetscBool :: istatechng

  state_change_threshold = 0.d0

  istatechng = PETSC_FALSE

  lid = option%liquid_phase
  gid = option%gas_phase
  pid = option%precipitate_phase
  pwid = option%pure_water_phase
  pbid = option%pure_brine_phase
  tgid = option%trapped_gas_phase
  spid = option%saturation_pressure_id

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id

  cpid = option%capillary_pressure_id
  co2_pressure_id = option%co2_pressure_id
  vpid = option%vapor_pressure_id
  rvpid = option%reduced_vapor_pressure_id

  if (sco2_auxvar%istatechng) return

  call SCO2ComputeSaltSolubility(sco2_auxvar%temp, salt_solubility)
  if (global_auxvar%istate == SCO2_GAS_STATE) then
    if (sco2_auxvar%m_salt(2) > epsilon) then
      xsl = salt_solubility
    else
      xsl = 0.d0
    endif
  else
    xsl = min(salt_solubility,sco2_auxvar%m_salt(1))
  endif
  call SCO2ComputeSurfaceTension(sco2_auxvar%temp, &
                                 xsl, sigma)
  beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma

  ! Max effective trapped gas saturation
  ! MAN: need to check how this works with Pc function Webb extensions
  sgt_max = characteristic_curves%saturation_function%Sgt_max

  !Compute capillary entry pressure. MAN: need to expand this list
  !This also probably does not need to be computed over and over
  sg_min = 1.d-3
  Pc_entry = 0.d0
  select type(sf => characteristic_curves%saturation_function)
    class is (sat_func_vg_type)
      sg_min = 1.0d1**(-3.d0+log10(1.d0/ &
               characteristic_curves%saturation_function%GetAlpha_()))
      sg_min = min(max(sg_min,1.d-4),1.d-3)
    class default
  end select

  old_state = global_auxvar%istate

  ! Decide if state changes
  select case(global_auxvar%istate)
    case(SCO2_LIQUID_STATE)
      cell_pressure = max(sco2_auxvar%pres(lid),sco2_auxvar%pres(vpid))
      ! cell_pressure = min(cell_pressure, 1.d8)          
      call SCO2Equilibrate(sco2_auxvar%temp,sco2_auxvar%pres(lid), &
                           Pco2,Pv, &
                           sco2_auxvar%pres(spid), &
                           sco2_auxvar%pres(rvpid), &
                           xco2g, xwg, xco2l, xsl, xwl, &
                           xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl, option)
      ! Check if dissolved CO2 exceeds solubility
      if (sco2_auxvar%xmass(co2_id,lid) > xco2l * (1.d0 + &
          state_change_threshold)) then 
        ! Compute what gas saturation would be
        call EOSGasDensity(sco2_auxvar%temp,Pco2, &
                     den_mol,drho_dT,drho_dP,ierr)
        den_co2 = den_mol * fmw_comp(2)
        call SCO2BrineDensity(sco2_auxvar%temp,cell_pressure, &
                             xsl, den_brine, option)
        call SCO2DensityCompositeLiquid(sco2_auxvar%temp, &
                                  den_brine,xco2l,den_liq)
        sg_est =  (sco2_auxvar%xmass(co2_id,lid) - xco2l * (1.d0 + &
                   state_change_threshold)) * &
                   den_liq / den_co2
        ! Check to see if gas bubbles out
        if (sg_est < sg_min) then
          ! No state change
          istatechng = PETSC_FALSE
          ! MAN: STOMP updates current gas pressure here. This would 
          !      just get overwritten in AuxVarCompute
          sco2_auxvar%pres(gid) = sco2_auxvar%pres(lid) + &
                  Pc_entry / beta_gl - eps_sl
        else
          sl_temp = 1.d0 - min(sg_est, 1.d-1)
          sgt_temp = 0.d0
          !MAN: sco2_auxvar%sat(tgid) should be 0
          call SCO2ComputePcHysteresis(characteristic_curves, &
                                       sl_temp, &
                                       sco2_auxvar%sat(tgid), &
                                       beta_gl, Pc, option)
          Pc = min(Pc,Pc_entry / beta_gl + 1.d5)
          
          ! State has changed, so update state and one primary variable
          sco2_auxvar%pres(gid) = sco2_auxvar%pres(lid) + Pc
          istatechng = PETSC_TRUE
          global_auxvar%istate = SCO2_LIQUID_GAS_STATE
        endif

      else

        Pg = sco2_auxvar%pres(lid) + &
             (Pc_entry / beta_gl) - eps_sl
        ! No state change
        istatechng = PETSC_FALSE
        
      endif
    case(SCO2_GAS_STATE)
      Pv = sco2_auxvar%pres(gid) - &
                sco2_auxvar%pres(co2_pressure_id)
      prvap = sco2_auxvar%pres(rvpid)

      if (pv > prvap * (1.d0 + state_change_threshold)) then
        ! Aqueous phase appears, transition state. Update primary variables
        ! for salt mass and liquid pressure
        istatechng = PETSC_TRUE
        global_auxvar%istate = SCO2_LIQUID_GAS_STATE

        sl_temp = (pv - prvap * (1.d0 + &
                   state_change_threshold)) / (prvap * (1.d0 + &
                   state_change_threshold))
        sco2_auxvar%m_salt(1) = sco2_auxvar%m_salt(2) / &
                                (sco2_auxvar%den_kg(lid) * &
                                sl_temp * &
                                material_auxvar%porosity)
        call SCO2ComputeSurfaceTension(sco2_auxvar%temp, &
                                       xsl, sigma)
        beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma
        sgt_temp = 0.d0
        call SCO2ComputePcHysteresis(characteristic_curves, &
                                     sco2_auxvar%sat(lid), &
                                     sco2_auxvar%sat(tgid), &
                                     beta_gl, Pc, option)
        sco2_auxvar%pres(lid) = sco2_auxvar%pres(gid) - &
                                Pc
      endif

    case(SCO2_TRAPPED_GAS_STATE)
      call SCO2ComputeSurfaceTension(sco2_auxvar%temp, &
                                     xsl, sigma)
      beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma
      Slr = characteristic_curves%saturation_function%Sr
      Sgt = sco2_auxvar%sat(tgid)

      if ((sgt_max - Sgt) > epsilon) then
        sl_temp = ((sgt_max - Sgt) / (sgt_max + sgt_max * Sgt - Sgt))
      else
        sl_temp = 0.d0
      endif
      cell_pressure = sco2_auxvar%pres(lid)
      ! cell_pressure = min(cell_pressure, 1.d8)          

      call SCO2BrineSaturationPressure(sco2_auxvar%temp, xsl, &
                                       Psat)
      call SCO2Equilibrate(sco2_auxvar%temp,cell_pressure, &
                           sco2_auxvar%pres(co2_pressure_id), &
                           Pv, Psat, &
                           sco2_auxvar%pres(rvpid), &
                           xco2g, xwg, xco2l, xsl, xwl, &
                           xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl, option)

      if (sco2_auxvar%sat(gid) < epsilon) then
        ! Transition to fully saturated without trapped gas
        istatechng = PETSC_TRUE
        global_auxvar%istate = SCO2_LIQUID_STATE

      elseif (((Sgt - sco2_auxvar%sg_trapped) > epsilon) .or. &
              ((Sgt - sgt_max) > epsilon)) then
        ! Trapped gas saturation is increasing or exceeds max trapped gas sat
        sgt_temp = sco2_auxvar%sg_trapped
        call SCO2ComputePcHysteresis(characteristic_curves, &
                                sco2_auxvar%sat(lid), Sgt_temp, &
                                beta_gl, Pc, option)
        Pc = min(Pc,(Pc_entry / beta_gl) + 1.d5)
        sco2_auxvar%pres(gid) = sco2_auxvar%pres(lid) + Pc
        istatechng = PETSC_TRUE
        global_auxvar%istate = SCO2_LIQUID_GAS_STATE

      else
        ! Remain in liquid saturated, trapped gas state
        istatechng = PETSC_FALSE
      endif

    case(SCO2_LIQUID_GAS_STATE)

      ! Compute Saturation including Hysteresis
      call SCO2ComputeSatHysteresis(characteristic_curves, &
                                    sco2_auxvar%pres(cpid), &
                                    sco2_auxvar%sl_min, &
                                    beta_gl, sco2_auxvar%den_kg(lid), &
                                    sl_temp, &
                                    sgt_temp, &
                                    option)
      sl_temp = sl_temp + sgt_temp
      if (dabs(1.d0 - sl_temp) < epsilon) then

        if (sgt_max > epsilon .and. sgt_temp > epsilon) then
          ! Fully liquid saturated, with trapped gas, update 1 
          ! primary variable

          istatechng = PETSC_TRUE
          global_auxvar%istate = SCO2_TRAPPED_GAS_STATE
          
          sl_temp = sl_temp - sgt_temp
          sco2_auxvar%sat(gid) = 1.d0 - sco2_auxvar%sat(lid)

        else
          ! No trapped gas, just liquid. Update 1 primary variable.
          istatechng = PETSC_TRUE
          global_auxvar%istate = SCO2_LIQUID_STATE

          cell_pressure = sco2_auxvar%pres(lid)
          ! cell_pressure = min(cell_pressure, 1.d8)          
          call SCO2BrineSaturationPressure(sco2_auxvar%temp,xsl, &
                                           Psat)
          call SCO2Equilibrate(sco2_auxvar%temp,cell_pressure, &
                               sco2_auxvar%pres(co2_pressure_id), &
                               Pv, Psat, &
                               sco2_auxvar%pres(rvpid), &
                               xco2g, xwg, xco2l, xsl, xwl, &
                               xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl, &
                               option)
          sco2_auxvar%xmass(co2_id,lid) = xco2l
        endif

      elseif (sl_temp < epsilon .and. &
              (1.d0 - sco2_auxvar%sat(lid)) > epsilon) then
        ! Transition to fully gas-saturated. Update 2 primary variables.
        istatechng = PETSC_TRUE
        global_auxvar%istate = SCO2_GAS_STATE

      else
        ! No state transition
        istatechng = PETSC_FALSE
      endif
  end select

  new_state = global_auxvar%istate

  ! Update primary variables
  if (istatechng) then

    select case(old_state)
      case(SCO2_LIQUID_STATE)
        state_change_string = 'Liquid --> '
      case(SCO2_GAS_STATE)
        state_change_string = 'Gas --> '
      case(SCO2_TRAPPED_GAS_STATE) 
        state_change_string = 'Trapped Gas --> '
      case(SCO2_LIQUID_GAS_STATE)
        state_change_string = 'Liquid & Mobile Gas -->'
    end select

    select case(new_state)
      case(SCO2_LIQUID_STATE)
        state_change_string = trim(state_change_string) // ' Liquid'
      case(SCO2_GAS_STATE)
        state_change_string = trim(state_change_string) // ' Gas'
      case(SCO2_TRAPPED_GAS_STATE) 
        state_change_string = trim(state_change_string) // ' Trapped Gas'
      case(SCO2_LIQUID_GAS_STATE)
        state_change_string = trim(state_change_string) // &
                              ' Liquid & Mobile Gas'
    end select

    if (option%iflag == SCO2_UPDATE_FOR_ACCUM) then
      write(append,'('' at Cell '',i8)') natural_id
    else if (option%iflag == SCO2_UPDATE_FOR_DERIVATIVE) then
      write(append, &
            '(''(due to perturbation) '',i8)') natural_id
    else
      write(append, &
             '('' at Boundary Face '', i8)') natural_id
    endif

    state_change_string = trim(state_change_string) // trim(append)

    if (sco2_restrict_state_chng) sco2_auxvar%istatechng = PETSC_TRUE

    select case(global_auxvar%istate)
      case(SCO2_LIQUID_STATE)

        x(SCO2_LIQUID_PRESSURE_DOF) = sco2_auxvar%pres(lid)
        x(SCO2_CO2_MASS_FRAC_DOF) = sco2_auxvar%xmass(co2_id,lid)
        ! x(SCO2_TEMPERATURE_DOF) = sco2_auxvar%temp
        x(SCO2_SALT_MASS_FRAC_DOF) = sco2_auxvar%m_salt(1)

      case(SCO2_GAS_STATE)

        x(SCO2_GAS_PRESSURE_DOF) = sco2_auxvar%pres(gid)
        x(SCO2_CO2_PRESSURE_DOF) = sco2_auxvar%pres(co2_pressure_id) 
        ! x(SCO2_TEMPERATURE_DOF) = sco2_auxvar%temp
        x(SCO2_SALT_MASS_FRAC_DOF) = sco2_auxvar%m_salt(2)

      case(SCO2_TRAPPED_GAS_STATE)

        x(SCO2_LIQUID_PRESSURE_DOF) = sco2_auxvar%pres(lid)
        x(SCO2_GAS_SATURATION_DOF) = sco2_auxvar%sat(gid)
        ! x(SCO2_TEMPERATURE_DOF) = sco2_auxvar%temp
        x(SCO2_SALT_MASS_FRAC_DOF) = sco2_auxvar%m_salt(1)

      case(SCO2_LIQUID_GAS_STATE)

        x(SCO2_LIQUID_PRESSURE_DOF) = sco2_auxvar%pres(lid)
        x(SCO2_TWO_PHASE_GAS_PRES_DOF) = sco2_auxvar%pres(gid)
        ! x(SCO2_TEMPERATURE_DOF) = sco2_auxvar%temp
        x(SCO2_SALT_MASS_FRAC_DOF) = sco2_auxvar%m_salt(1)

      case default

        write(option%io_buffer,*) global_auxvar%istate
        option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
                           ') not recognized in SCO2AuxVarUpdateState.'
        call PrintErrMsgByRank(option)

    end select

    ! Update secondary variables
    call SCO2AuxVarCompute(x, sco2_auxvar, global_auxvar, material_auxvar, &
                           characteristic_curves, sco2_parameter, natural_id, &
                           option)
    
    state_change_string = 'State Transition: ' // trim(state_change_string)
    if (sco2_print_state_transition) then
      call PrintMsgByRank(option,state_change_string)
    endif

  endif


end subroutine SCO2AuxVarUpdateState

! ************************************************************************** !

subroutine SCO2AuxVarCompute(x,sco2_auxvar,global_auxvar,material_auxvar, &
                             characteristic_curves,sco2_parameter, &
                             natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell
  !
  ! Author: Michael Nole
  ! Date: 11/22/23
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_module
  use Fracture_module

  implicit none

  type(sco2_auxvar_type) :: sco2_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  type(sco2_parameter_type), pointer :: sco2_parameter
  PetscInt :: natural_id
  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  
  ! Phase ID's
  PetscInt :: lid, gid, pid, pwid, pgid, pbid
  ! Component ID's
  PetscInt :: wid, co2_id, co2_pressure_id, sid
  ! Other ID's
  PetscInt :: cpid, vpid, rvpid, spid, tgid

  PetscReal :: cell_pressure
  PetscReal :: xco2g, xwg, xco2l, xsl, xwl, xmolco2g, xmolwg, xmolco2l, &
               xmolsl, xmolwl
  PetscReal :: mw_mix
  PetscReal :: den_mol, den_steam_kg
  PetscReal :: den_steam
  PetscReal :: salt_solubility, x_salt_dissolved
  PetscReal :: beta_gl
  PetscReal :: sigma
  PetscReal :: dkrl_dsatl, dkrg_dsatl
  PetscReal :: visc_water, visc_brine, visc_co2
  PetscReal :: sl_temp, pva
  PetscErrorCode :: ierr

  ! Unused
  PetscReal :: dpor_dp, drho_dp, drho_dT

  PetscReal, parameter :: epsilon = 1.d-14

  lid = option%liquid_phase
  gid = option%gas_phase
  pid = option%precipitate_phase
  pwid = option%pure_water_phase
  pbid = option%pure_brine_phase
  tgid = option%trapped_gas_phase
  pgid = tgid ! pure gas

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id

  cpid = option%capillary_pressure_id
  co2_pressure_id = option%co2_pressure_id
  vpid = option%vapor_pressure_id
  rvpid = option%reduced_vapor_pressure_id
  spid = option%saturation_pressure_id

  sco2_auxvar%H = 0.d0
  sco2_auxvar%U = 0.d0
  sco2_auxvar%pres = 0.d0
  sco2_auxvar%sat = 0.d0
  sco2_auxvar%den_kg = 0.d0
  sco2_auxvar%xmass = 0.d0
  sco2_auxvar%xmol = 0.d0
  sco2_auxvar%effective_diffusion_coeff = 0.d0
  sco2_auxvar%dispersivity = 0.d0
  sco2_auxvar%effective_porosity = material_auxvar%porosity_base
  sco2_auxvar%effective_permeability = 0.d0
  sco2_auxvar%tortuosity = 0.d0

  sco2_auxvar%mobility = 0.d0
  sco2_auxvar%kr = 0.d0

  ! Precipitate Phase has constant pure salt properties
  sco2_auxvar%den_kg(pid) = SALT_DENSITY_KG
  sco2_auxvar%den(pid) = sco2_auxvar%den_kg(pid) / fmw_comp(3)
  sco2_auxvar%xmass(wid,pid) = 0.d0
  sco2_auxvar%xmass(co2_id,pid) = 0.d0
  sco2_auxvar%xmass(sid,pid) = 1.d0
  !sco2_auxvar%m_salt = 0.d0

  sco2_auxvar%temp = sco2_isothermal_temperature

  select case(global_auxvar%istate)
  ! Potential primary variables need to be given values here:
  ! liquid pressure, gas pressure, CO2 partial pressure, 
  ! dissolved CO2 mass fraction, temperature, total NaCl brine mass fraction.
    case(SCO2_LIQUID_STATE)

      ! State: Saturated system without trapped gas
      ! Primary Variables: 
      !               Liquid Pressure, Aqueous CO2 mass fraction,
      !               Temperature, total NaCl brine fraction (kg NaCl/kg brine)

      ! MAN: I wonder if this way runs into issues with lagging dissolved
      !      salt mass fraction

      ! Primary Variables
      sco2_auxvar%pres(lid) = x(SCO2_LIQUID_PRESSURE_DOF)
      sco2_auxvar%xmass(co2_id,lid) = x(SCO2_CO2_MASS_FRAC_DOF)
      ! sco2_auxvar%temp = x(SCO2_TEMPERATURE_DOF)

      ! This is the total salt mass fraction including precipitate phase, 
      ! neglecting the mass of dissolved CO2.
      sco2_auxvar%m_salt(1) = x(SCO2_SALT_MASS_FRAC_DOF)

      sco2_auxvar%pres(cpid) = 0.d0
      ! Starting guess for Equilibrate
      sco2_auxvar%xmass(sid,lid) = sco2_auxvar%m_salt(1)
      ! Secondary Variables
      
      ! kg NaCl/kg liquid
      call SCO2ComputeSaltSolubility(sco2_auxvar%temp, salt_solubility)
      ! Dissolved salt mass fraction
      x_salt_dissolved = min(sco2_auxvar%m_salt(1),salt_solubility)
      call SCO2BrineSaturationPressure(sco2_auxvar%temp, &
                                         x_salt_dissolved, &
                                         sco2_auxvar%pres(spid))
      ! Brine density
      call SCO2BrineDensity(sco2_auxvar%temp, sco2_auxvar%pres(lid), &
                            x_salt_dissolved, sco2_auxvar%den_kg(pbid), option)
      !aux(1) = x_salt_dissolved
      !call EOSWaterDensityExt(sco2_auxvar%temp,sco2_auxvar%pres(lid), &
      !                         aux, sco2_auxvar%den_kg(pbid), &
      !                         den_mol,ierr)
      ! Brine vapor pressure
      call SCO2VaporPressureBrine(sco2_auxvar%temp, sco2_auxvar%pres(spid), &
                                   sco2_auxvar%pres(cpid), &
                                   sco2_auxvar%den_kg(pbid), &
                                   x_salt_dissolved, sco2_auxvar%pres(rvpid))

      ! Pure water density
      call SCO2WaterDensity(sco2_auxvar%temp,sco2_auxvar%pres(rvpid), &
                            TWO_INTEGER,sco2_auxvar%den_kg(pwid), &
                            den_steam_kg,option)
      pva = max(sco2_auxvar%pres(lid) - sco2_auxvar%pres(rvpid), 0.d0)             
      xsl = x_salt_dissolved
      call SCO2Equilibrate(sco2_auxvar%temp,sco2_auxvar%pres(lid), &
                           sco2_auxvar%pres(co2_pressure_id), &
                           sco2_auxvar%pres(vpid), &
                           sco2_auxvar%pres(spid), &
                           sco2_auxvar%pres(rvpid), &
                           xco2g, xwg, xco2l, xsl, xwl, &
                           xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl, option)
      
      sco2_auxvar%xmass(sid,lid) = x_salt_dissolved + &
                              (xsl-x_salt_dissolved) * &
                              (sco2_auxvar%xmass(co2_id,lid) / xco2l)

      sco2_auxvar%xmass(wid,lid) = 1.d0 - sco2_auxvar%xmass(sid,lid) - &
                                   sco2_auxvar%xmass(co2_id,lid)
      sco2_auxvar%xmass(wid,lid) = max(sco2_auxvar%xmass(wid,lid),0.d0)
      sco2_auxvar%sat(lid) = 1.d0 !- sco2_auxvar%sat(pid)
      sco2_auxvar%sat(gid) = 0.d0

      ! Populate all pressures, even though gas phase is not present.
      sco2_auxvar%pres(gid) = sco2_auxvar%pres(lid)
      sco2_auxvar%xmass(co2_id,gid) = 1.d0

      ! Update the liquid mole fractions
      mw_mix = 1.d0 / (sco2_auxvar%xmass(wid,lid)/fmw_comp(1) + &
               sco2_auxvar%xmass(co2_id,lid)/fmw_comp(2) + &
               sco2_auxvar%xmass(sid,lid)/fmw_comp(3))
      sco2_auxvar%xmol(wid,lid) = sco2_auxvar%xmass(wid,lid)* &
                                  mw_mix/fmw_comp(1)
      sco2_auxvar%xmol(co2_id,lid) = sco2_auxvar%xmass(co2_id,lid)* &
                                     mw_mix/fmw_comp(2)
      sco2_auxvar%xmol(sid,lid) = sco2_auxvar%xmass(sid,lid)* &
                                  mw_mix/fmw_comp(3)

    case (SCO2_GAS_STATE)
      ! Fully unsaturated system with or without trapped gas.
      ! Primary Variables: 
      !               Gas Pressure, CO2 Partial Pressure,
      !               Temperature, total NaCl mass

      ! Primary Variables
      sco2_auxvar%pres(gid) = x(SCO2_GAS_PRESSURE_DOF)
      sco2_auxvar%pres(co2_pressure_id) = x(SCO2_CO2_PRESSURE_DOF)
      ! sco2_auxvar%temp = x(SCO2_TEMPERATURE_DOF)
      ! This PV is now total salt mass, not salt mass fraction in brine
      sco2_auxvar%m_salt(2) = x(SCO2_SALT_MASS_FRAC_DOF)

      ! Secondary Variables
      ! kg NaCl/kg liquid
      call SCO2ComputeSaltSolubility(sco2_auxvar%temp, salt_solubility)
      if (sco2_auxvar%m_salt(2) > epsilon) then
        x_salt_dissolved = salt_solubility
      else
        x_salt_dissolved = 0.d0
      endif
      sco2_auxvar%xmass(sid,lid) = x_salt_dissolved
      call SCO2ComputeSurfaceTension(sco2_auxvar%temp, &
                                     x_salt_dissolved, sigma)
      ! MAN: check the reference surface tension
      !MAN: sco2_auxvar%sat(tgid) should be small?
      beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma 
      call SCO2ComputePcHysteresis(characteristic_curves, &
                                   sco2_auxvar%sat(lid), &
                                   sco2_auxvar%sat(tgid), &
                                   beta_gl,sco2_auxvar%pres(cpid), option)
      sco2_auxvar%pres(cpid) = sco2_auxvar%pres(cpid) / beta_gl

      
      cell_pressure = max(sco2_auxvar%pres(gid),sco2_auxvar%pres(spid))
      ! cell_pressure = min(cell_pressure, 1.d8)          

      sco2_auxvar%pres(rvpid) = max(sco2_auxvar%pres(gid) - &
                               sco2_auxvar%pres(co2_pressure_id), 0.d0)
      sco2_auxvar%pres(vpid) = sco2_auxvar%pres(rvpid)
      pva = max(sco2_auxvar%pres(gid) - sco2_auxvar%pres(rvpid),0.d0)

      call SCO2WaterDensity(sco2_auxvar%temp,sco2_auxvar%pres(rvpid), &
                            TWO_INTEGER,sco2_auxvar%den_kg(pwid), &
                            den_steam_kg,option)
      ! Compute equilibrium mass and mole fractions for all components
      xsl = x_salt_dissolved
      call SCO2BrineSaturationPressure(sco2_auxvar%temp, &
                                         x_salt_dissolved, &
                                         sco2_auxvar%pres(spid))
      call SCO2Equilibrate(sco2_auxvar%temp,cell_pressure, &
                           pva, sco2_auxvar%pres(vpid), &
                           sco2_auxvar%pres(spid), &
                           sco2_auxvar%pres(rvpid), &
                           xco2g, xwg, xco2l, xsl, xwl, &
                           xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl, option)
      
      call SCO2ComputeSaltDensity(sco2_auxvar%temp, cell_pressure, &
                                  sco2_auxvar%den_kg(pid))
      
      sco2_auxvar%sat(pid) = sco2_auxvar%m_salt(2) / &
                             (sco2_auxvar%den_kg(pid) * &
                             material_auxvar%volume)
      sco2_auxvar%sat(pid) = max(min(sco2_auxvar%sat(pid),1.d0),0.d0)
      sco2_auxvar%sat(lid) = 0.d0
      sco2_auxvar%sat(gid) = 1.d0

      !sco2_auxvar%xmass(sid,lid) = sco2_auxvar%m_salt(2) * &
      !                             sco2_auxvar%den_kg(pid) * &
      !                             sco2_auxvar%sat(lid) * &
      !                             material_auxvar%porosity

      !sco2_auxvar%xmass(co2_id,lid) = xco2l
      !sco2_auxvar%xmass(wid,lid) = 1.d0 - sco2_auxvar%xmass(gid,lid) - &
      !                             sco2_auxvar%xmass(sid,lid)

      !sco2_auxvar%xmass(wid,gid) = xwg
      !sco2_auxvar%xmass(co2_id,gid) = xco2g

      ! Update mass fractions
      sco2_auxvar%xmass(co2_id,lid) = xco2l
      sco2_auxvar%xmass(wid,lid) = xwl
      sco2_auxvar%xmass(sid,lid) = xsl
      sco2_auxvar%xmass(co2_id,gid) = xco2g
      sco2_auxvar%xmass(wid,gid) = xwg     
      
      ! Update mole fractions
      sco2_auxvar%xmol(co2_id,lid) = xmolco2l
      sco2_auxvar%xmol(wid,lid) = xmolwl
      sco2_auxvar%xmol(sid,lid) = xmolsl
      sco2_auxvar%xmol(co2_id,gid) = xmolco2g
      sco2_auxvar%xmol(wid,gid) = xmolwg

      ! Brine density
      call SCO2BrineDensity(sco2_auxvar%temp, cell_pressure, &
                            x_salt_dissolved, sco2_auxvar%den_kg(pbid), option)
      !aux(1) = x_salt_dissolved
      !call EOSWaterDensityExt(sco2_auxvar%temp,cell_pressure, &
      !                         aux, sco2_auxvar%den_kg(pbid), &
      !                         den_mol,ierr)
      ! Liquid phase density (including CO2)
      call SCO2DensityCompositeLiquid(sco2_auxvar%temp,sco2_auxvar%den_kg(pbid), &
                                  sco2_auxvar%xmass(co2_id,lid), &
                                  sco2_auxvar%den_kg(lid))
      sco2_auxvar%pres(lid) = sco2_auxvar%pres(gid) - sco2_auxvar%pres(cpid) 
      call SCO2BrineSaturationPressure(sco2_auxvar%temp, &
                                         x_salt_dissolved, &
                                         sco2_auxvar%pres(spid))
    case (SCO2_TRAPPED_GAS_STATE)
      ! Fully liquid saturated system with trapped gas
      ! Primary Variables: 
      !               Liquid Pressure, Trapped Gas Saturation,
      !               Temperature, total NaCl brine fraction (kg NaCl/kg brine)
      sco2_auxvar%pres(lid) = x(SCO2_LIQUID_PRESSURE_DOF)
      sco2_auxvar%sat(gid) = x(SCO2_GAS_SATURATION_DOF)
      ! sco2_auxvar%temp = x(SCO2_TEMPERATURE_DOF)
      sco2_auxvar%m_salt(1) = x(SCO2_SALT_MASS_FRAC_DOF)

      ! Starting guess for Equilibrate
      sco2_auxvar%xmass(sid,lid) = sco2_auxvar%m_salt(1)

      ! Secondary Variables
      ! Trapped gas is disconnected
      sco2_auxvar%pres(cpid) = 0.d0 
      ! kg NaCl/kg liquid
      call SCO2ComputeSaltSolubility(sco2_auxvar%temp, salt_solubility)
      ! Dissolved salt mass fraction
      x_salt_dissolved = min(sco2_auxvar%m_salt(1),salt_solubility)
      call SCO2BrineSaturationPressure(sco2_auxvar%temp, &
                                         x_salt_dissolved, &
                                         sco2_auxvar%pres(spid))
      cell_pressure = max(sco2_auxvar%pres(gid),sco2_auxvar%pres(spid))
      ! cell_pressure = min(cell_pressure, 1.d8)

      ! Brine density
      call SCO2BrineDensity(sco2_auxvar%temp, cell_pressure, &
                            x_salt_dissolved, sco2_auxvar%den_kg(pbid), option)
      !aux(1) = x_salt_dissolved
      !call EOSWaterDensityExt(sco2_auxvar%temp,cell_pressure, &
      !                        aux, sco2_auxvar%den_kg(pbid), &
      !                        den_mol,ierr)
      ! Brine vapor pressure
      call SCO2VaporPressureBrine(sco2_auxvar%temp, sco2_auxvar%pres(spid), &
                                   sco2_auxvar%pres(cpid), &
                                   sco2_auxvar%den_kg(pbid), &
                                   x_salt_dissolved, sco2_auxvar%pres(rvpid))

      ! Pure water density
      call SCO2WaterDensity(sco2_auxvar%temp,sco2_auxvar%pres(rvpid), &
                            TWO_INTEGER,sco2_auxvar%den_kg(pwid), &
                            den_steam_kg,option)
      pva = max(sco2_auxvar%pres(gid) - sco2_auxvar%pres(rvpid), 0.d0)
      sco2_auxvar%pres(vpid) = sco2_auxvar%pres(rvpid)
      ! Compute equilibrium mass and mole fractions for all components
      xsl = x_salt_dissolved
      call SCO2Equilibrate(sco2_auxvar%temp,cell_pressure, &
                           sco2_auxvar%pres(co2_pressure_id), &
                           sco2_auxvar%pres(vpid), &
                           sco2_auxvar%pres(spid), &
                           sco2_auxvar%pres(rvpid), &
                           xco2g, xwg, xco2l, xsl, xwl, &
                           xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl, option)

      ! Update trapped gas
      sco2_auxvar%sat(tgid) = sco2_auxvar%sat(gid)

      ! Update mass fractions
      sco2_auxvar%xmass(co2_id,lid) = xco2l
      sco2_auxvar%xmass(wid,lid) = xwl
      sco2_auxvar%xmass(sid,lid) = xsl
      sco2_auxvar%xmass(co2_id,gid) = xco2g
      sco2_auxvar%xmass(wid,gid) = xwg

    case (SCO2_LIQUID_GAS_STATE)
      ! State: Unsaturated, with or without trapped gas
      ! Primary Variables: 
      !               Liquid Pressure, Gas Pressure,
      !               Temperature, total NaCl brine fraction (kg NaCl/kg brine)
      sco2_auxvar%pres(lid) = x(SCO2_LIQUID_PRESSURE_DOF)
      sco2_auxvar%pres(gid) = x(SCO2_TWO_PHASE_GAS_PRES_DOF)
      ! sco2_auxvar%temp = x(SCO2_TEMPERATURE_DOF)
      sco2_auxvar%m_salt(1) = x(SCO2_SALT_MASS_FRAC_DOF)

      ! Starting guess for Equilibrate
      sco2_auxvar%xmass(sid,lid) = sco2_auxvar%m_salt(1)
      
      ! Secondary Variables

      sco2_auxvar%pres(cpid) = max(sco2_auxvar%pres(gid) - &
                                   sco2_auxvar%pres(lid), 0.d0)

      ! kg NaCl/kg liquid
      call SCO2ComputeSaltSolubility(sco2_auxvar%temp, salt_solubility)
      ! Dissolved salt mass fraction
      x_salt_dissolved = min(sco2_auxvar%m_salt(1),salt_solubility)
      call SCO2BrineSaturationPressure(sco2_auxvar%temp, &
                                         x_salt_dissolved, &
                                         sco2_auxvar%pres(spid))
      cell_pressure = max(sco2_auxvar%pres(gid),sco2_auxvar%pres(spid))
      ! cell_pressure = min(cell_pressure, 1.d8)          
      ! Brine density
      call SCO2BrineDensity(sco2_auxvar%temp, cell_pressure, &
                           x_salt_dissolved, sco2_auxvar%den_kg(pbid), option)
      ! aux(1) = x_salt_dissolved
      ! call EOSWaterDensityExt(sco2_auxvar%temp,cell_pressure, &
      !                         aux, sco2_auxvar%den_kg(pbid), &
      !                         den_mol,ierr)
      ! Brine vapor pressure
      call SCO2VaporPressureBrine(sco2_auxvar%temp, sco2_auxvar%pres(spid), &
                                   sco2_auxvar%pres(cpid), &
                                   sco2_auxvar%den_kg(pbid), &
                                   x_salt_dissolved, sco2_auxvar%pres(rvpid))

      ! Pure water density
      xsl = x_salt_dissolved
      call SCO2WaterDensity(sco2_auxvar%temp,sco2_auxvar%pres(rvpid), &
                            TWO_INTEGER,sco2_auxvar%den_kg(pwid), &
                            den_steam_kg,option)                             
      pva = max(sco2_auxvar%pres(gid) - sco2_auxvar%pres(rvpid), 0.d0)
      sco2_auxvar%pres(vpid) = sco2_auxvar%pres(rvpid)
      call SCO2Equilibrate(sco2_auxvar%temp,cell_pressure, &
                           sco2_auxvar%pres(co2_pressure_id), &
                           sco2_auxvar%pres(vpid), &
                           sco2_auxvar%pres(spid), &
                           sco2_auxvar%pres(rvpid), &
                           xco2g, xwg, xco2l, xsl, xwl, &
                           xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl, option)
      
      ! Update mass fractions
      sco2_auxvar%xmass(co2_id,lid) = xco2l
      sco2_auxvar%xmass(wid,lid) = xwl
      sco2_auxvar%xmass(sid,lid) = xsl
      sco2_auxvar%xmass(co2_id,gid) = xco2g
      sco2_auxvar%xmass(wid,gid) = xwg

      ! Update mole fractions
      sco2_auxvar%xmol(co2_id,lid) = xmolco2l
      sco2_auxvar%xmol(wid,lid) = xmolwl
      sco2_auxvar%xmol(sid,lid) = xmolsl
      sco2_auxvar%xmol(co2_id,gid) = xmolco2g
      sco2_auxvar%xmol(wid,gid) = xmolwg

    case default

      write(option%io_buffer,*) global_auxvar%istate
      option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
        ') not recognized in SCO2AuxVarCompute.'
      call PrintErrMsgByRank(option)

  end select

  ! ! Update the liquid mole fractions
  ! mw_mix = 1.d0 / (sco2_auxvar%xmass(wid,lid)/fmw_comp(1) + &
  !                  sco2_auxvar%xmass(co2_id,lid)/fmw_comp(2) + &
  !                  sco2_auxvar%xmass(sid,lid)/fmw_comp(3))
  ! sco2_auxvar%xmol(wid,lid) = sco2_auxvar%xmass(wid,lid)* &
  !                             mw_mix/fmw_comp(1)
  ! sco2_auxvar%xmol(co2_id,lid) = sco2_auxvar%xmass(co2_id,lid)* &
  !                                mw_mix/fmw_comp(2)
  ! sco2_auxvar%xmol(sid,lid) = sco2_auxvar%xmass(sid,lid)* &
  !                             mw_mix/fmw_comp(3)

  ! ! Update the gas mole fractions
  ! mw_mix = 1.d0 / (sco2_auxvar%xmass(wid,gid)/fmw_comp(1) + &
  !                  sco2_auxvar%xmass(co2_id,gid)/fmw_comp(2))
  ! sco2_auxvar%xmol(wid,gid) = sco2_auxvar%xmass(wid,gid)* &
  !                             mw_mix/fmw_comp(1)
  ! sco2_auxvar%xmol(co2_id,gid) = sco2_auxvar%xmass(co2_id,gid)* &
  !                                mw_mix/fmw_comp(2)

  cell_pressure = max(sco2_auxvar%pres(lid),sco2_auxvar%pres(gid), &
                      sco2_auxvar%pres(spid))
  ! cell_pressure = min(cell_pressure, 1.d8)          
  sco2_auxvar%xmass(co2_id,gid) = xco2g
  sco2_auxvar%xmass(wid,gid) = 1.d0 - xco2g
  sco2_auxvar%xmol(co2_id,gid) = xmolco2g
  sco2_auxvar%xmol(wid,gid) = 1.d0 - xmolco2g 

  ! Update Porosity
  ! MAN: need to update the porosity compressibility model.
  !if (option%iflag /= SCO2_UPDATE_FOR_BOUNDARY) then
    dpor_dp = 0.d0
    sco2_auxvar%effective_porosity = material_auxvar%porosity_base
    ! creep_closure, fracture, and soil_compressibility are mutually exclusive
    if (associated(material_auxvar%fracture)) then
      call FracturePoroEvaluate(material_auxvar,cell_pressure, &
                                sco2_auxvar%effective_porosity,dpor_dp)
    else if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                sco2_auxvar%effective_porosity,dpor_dp)
    endif
    if (option%iflag /= SCO2_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = sco2_auxvar%effective_porosity
    endif
  !endif

  ! Gas phase density
  call EOSGasDensity(sco2_auxvar%temp,pva, &
                     den_mol,drho_dT,drho_dP,ierr)
  sco2_auxvar%den_kg(co2_pressure_id) = den_mol * fmw_comp(2)
  sco2_auxvar%den_kg(gid) = sco2_auxvar%xmass(co2_id,gid) * &
                            sco2_auxvar%den_kg(co2_pressure_id) + &
                            sco2_auxvar%xmass(wid,gid) * &
                            den_steam_kg
  ! Gas phase viscosity
  call SCO2ViscosityWater(sco2_auxvar%temp,sco2_auxvar%pres(vpid), &
                          den_steam_kg,visc_water,option)
  call SCO2ViscosityCO2(sco2_auxvar%temp,sco2_auxvar%den_kg(co2_pressure_id), &
                        visc_co2)
  call SCO2ViscosityGas(visc_water,visc_co2,sco2_auxvar%xmol(wid,gid), &
                        sco2_auxvar%xmol(co2_id,gid),sco2_auxvar%visc(gid))
  
  ! Liquid phase density (including CO2)
  call SCO2DensityCompositeLiquid(sco2_auxvar%temp,sco2_auxvar%den_kg(pbid), &
                                  sco2_auxvar%xmass(co2_id,lid), &
                                  sco2_auxvar%den_kg(lid))
  ! Liquid phase viscosity
  call SCO2WaterDensity(sco2_auxvar%temp, cell_pressure, ONE_INTEGER, &
                        sco2_auxvar%den_kg(pwid), den_steam, &
                        option)
  call SCO2ViscosityWater(sco2_auxvar%temp,cell_pressure, &
                          sco2_auxvar%den_kg(pwid),visc_water,option)
  call SCO2ViscosityBrine(sco2_auxvar%temp, sco2_auxvar%xmass(sid,lid), &
                          visc_water, visc_brine)
  call SCO2ViscosityLiquid(sco2_auxvar%xmass(co2_id,lid), visc_brine, &
                           visc_co2, sco2_auxvar%visc(lid))

  ! CO2-water surface tension 
  ! MAN: doesn't do anything here right now
  call SCO2ComputeSurfaceTension(sco2_auxvar%temp,sco2_auxvar%xmass(sid,lid), &
                                 sigma)
  beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma

  if (global_auxvar%istate /= SCO2_GAS_STATE) then
    call SCO2ComputeSatHysteresis(characteristic_curves, &
                                    sco2_auxvar%pres(cpid), &
                                    sco2_auxvar%sl_min, &
                                    beta_gl, sco2_auxvar%den_kg(lid), &
                                    sl_temp, &
                                    sco2_auxvar%sat(tgid), &
                                    option)
    sco2_auxvar%sat(lid) = sl_temp
    sco2_auxvar%sat(gid) = 1.d0 - sco2_auxvar%sat(lid)
  endif
  ! sco2_auxvar%sat(gid) is always mobile gas + trapped gas
  !sco2_auxvar%sat(gid) = (1.d0 - sl_temp) + sco2_auxvar%sat(tgid)
  ! Compute relative permeabilities
  ! MAN: Check if surface tension needs to be incorported into rel perm.
  ! MAN: Need to check that trapped gas presence is treated properly in 
  !      rel perm calcs
  call characteristic_curves%liq_rel_perm_function% &
           RelativePermeability(sco2_auxvar%sat(lid),sco2_auxvar%kr(lid), &
                                dkrl_dsatl,option)
  !sco2_auxvar%kr(lid) = min(max(sco2_auxvar%kr(lid),1.d-24),1.d0)
  call characteristic_curves%gas_rel_perm_function% &
           RelativePermeability(sco2_auxvar%sat(lid),sco2_auxvar%kr(gid), &
                                dkrg_dsatl,option)

  ! Convert to molar density: liquid
  mw_mix = sco2_auxvar%xmol(wid,lid) * fmw_comp(1) + &
          sco2_auxvar%xmol(co2_id,lid) * fmw_comp(2) + &
          sco2_auxvar%xmol(sid,lid) * fmw_comp(3)
  sco2_auxvar%den(lid) = sco2_auxvar%den_kg(lid) / mw_mix
  
  ! Convert to molar density: gas
  mw_mix = sco2_auxvar%xmol(wid,gid) * fmw_comp(1) + &
          sco2_auxvar%xmol(co2_id,gid) * fmw_comp(2) + &
          sco2_auxvar%xmol(sid,gid) * fmw_comp(3)
  sco2_auxvar%den(gid) = sco2_auxvar%den_kg(gid) / mw_mix

  ! Tortuosity
  call SCO2Tortuosity(sco2_auxvar%sat(lid), sco2_auxvar%sat(gid), &
                      sco2_auxvar%effective_porosity, &
                      sco2_auxvar%tortuosity(lid), &
                      sco2_auxvar%tortuosity(gid))
  ! Update Diffusivities: water vapor, dissolved CO2, dissolved salt
  call SCO2DiffusionCoeff(sco2_auxvar%temp, cell_pressure, &
                          sco2_auxvar%xmass(sid,lid), &
                          sco2_auxvar%visc(lid), &
                          sco2_parameter, option)
  call SCO2ComputeEffectiveDiffusion(sco2_parameter, sco2_auxvar, option)

  ! Precipitate salt
  sco2_auxvar%xmass(sid,pid) = 1.d0
  sco2_auxvar%xmol(sid,pid) = 1.d0

  ! Salt precipitate density and saturation
  call SCO2ComputeSaltDensity(sco2_auxvar%temp, cell_pressure, &
                              sco2_auxvar%den_kg(pid))
  if (global_auxvar%istate == SCO2_GAS_STATE .or. &
      sco2_auxvar%sat(lid) == 0.d0) then
    sco2_auxvar%sat(pid) = sco2_auxvar%m_salt(2) / (sco2_auxvar%den_kg(pid) * &
                           sco2_auxvar%effective_porosity)
    !MAN: not sure why we need this:
    sco2_auxvar%m_salt(1) = sco2_auxvar%m_salt(2) * sco2_auxvar%den_kg(pbid) * &
                            epsilon * sco2_auxvar%effective_porosity
  else
    sco2_auxvar%sat(pid) = max(sco2_auxvar%m_salt(1) - salt_solubility, &
                           0.d0) * sco2_auxvar%den_kg(pbid) * &
                           sco2_auxvar%sat(lid) / &
                           sco2_auxvar%den_kg(pid)
    sco2_auxvar%m_salt(2) = sco2_auxvar%m_salt(1) * sco2_auxvar%den_kg(pbid) * &
                            sco2_auxvar%sat(lid) * &
                            sco2_auxvar%effective_porosity 
  endif

  ! Permeability and porosity reduction with salt precipitate effects
  call SCO2ScalePermPhi(sco2_auxvar, material_auxvar, global_auxvar, option)

  sco2_auxvar%mobility(lid) = sco2_auxvar%kr(lid) / sco2_auxvar%visc(lid)
  sco2_auxvar%mobility(gid) = sco2_auxvar%kr(gid) / sco2_auxvar%visc(gid)

  ! Energy calculations

  ! Brine enthalpy
  ! aux(1) = sco2_auxvar%xmass(sid,lid)
  ! call EOSWaterEnthalpyExt(sco2_auxvar%temp,cell_pressure, &
  !                          aux,sco2_auxvar%H(pwid),ierr)
  ! ! CO2 density, internal energy, enthalpy
  ! call EOSGasDensityEnergy(sco2_auxvar%temp,sco2_auxvar% &
  !                          pres(co2_pressure_id),den_co2, &
  !                          sco2_auxvar%H(pgid),sco2_auxvar%U(pgid),ierr)
  ! ! Liquid phase enthalpy
  ! sco2_auxvar%H(lid) = SCO2EnthalpyCompositeLiquid(sco2_auxvar%temp, &
  !                                  sco2_auxvar%xmass(sid,lid), &
  !                                  sco2_auxvar%xmass(co2_id,lid), &
  !                                  sco2_auxvar%H(pwid), sco2_auxvar%H(pgid))

  ! sco2_auxvar%H(pwid) = sco2_auxvar%H(pwid) * 1.d-6 ! J/kg -> MJ/kg
  ! sco2_auxvar%H(lid) = sco2_auxvar%H(lid) * 1.d-6 ! J/kg -> MJ/kg
  ! ! MJ/kg comp
  ! sco2_auxvar%U(lid) = (sco2_auxvar%H(lid) - &
  !                       ! Pa / kg/m^3 * 1.e-6 = MJ/kg
  !                       (cell_pressure / sco2_auxvar%den_kg(lid) * &
  !                       1.d-6))
  ! sco2_auxvar%U(pwid) = (sco2_auxvar%H(pwid) - &
  !                       ! Pa / kg/m^3 * 1.e-6 = MJ/kg
  !                       (cell_pressure / sco2_auxvar%den_kg(pwid) * &
  !                       1.d-6))

  ! sco2_auxvar%H(pgid) = sco2_auxvar%H(pgid) / fmw_comp(co2_id) * 1.d-6 ! MJ/kg
  ! sco2_auxvar%U(pgid) = sco2_auxvar%U(pgid) / fmw_comp(co2_id) * 1.d-6 ! MJ/kg
  ! if (sco2_auxvar%pres(rvpid) > 0.d0) then
  !   call EOSWaterSteamDensityEnthalpy(sco2_auxvar%temp, &
  !                                   sco2_auxvar%pres(rvpid), &
  !                                   den_steam_kg, &
  !                                   sco2_auxvar%den(stid), &
  !                                   sco2_auxvar%H(stid),ierr)
  ! else
  !   sco2_auxvar%den(stid) = 0.d0
  !   sco2_auxvar%H(stid) = 0.d0
  ! endif
  ! ! J/kmol -> MJ/kg 
  ! sco2_auxvar%H(stid) = sco2_auxvar%H(stid) / fmw_comp(wid) * 1.d-6                                  
  ! sco2_auxvar%U(stid) = sco2_auxvar%H(stid) - &
  !                       sco2_auxvar%pres(vpid) / den_steam_kg
  
  ! ! Gas phase enthalpy
  ! sco2_auxvar%H(gid) = sco2_auxvar%xmass(wid,gid) * sco2_auxvar%H(stid) + &
  !                      sco2_auxvar%xmass(co2_id,gid) * sco2_auxvar%H(pgid)
  ! sco2_auxvar%U(gid) = sco2_auxvar%U(gid) - &
  !                      ! Pa / kg/m^3 * 1.e-6 = MJ/kg
  !                      sco2_auxvar%pres(gid) / sco2_auxvar%den_kg(gid) * 1.d-6

  ! ! Precipitate phase enthalpy
  ! call SCO2SaltEnthalpy(sco2_auxvar%temp,sco2_auxvar%H(pid))
  ! ! MJ/kg
  ! sco2_auxvar%H(pid) = sco2_auxvar%H(pid) * 1.d-6
  ! sco2_auxvar%U(pid) = sco2_auxvar%H(pid)

  ! sco2_auxvar%mobility(lid) = sco2_auxvar%kr(lid) / sco2_auxvar%visc(lid)
  ! sco2_auxvar%mobility(gid) = sco2_auxvar%kr(gid) / sco2_auxvar%visc(gid)

end subroutine SCO2AuxVarCompute

! ************************************************************************** !

subroutine SCO2VaporPressureBrine(T,P_sat,Pc,rho_kg,x_salt,P_vap)
  !
  ! Computes the reduced vapor pressure of water following the Kelvin equation.
  !
  ! Author: Michael Nole
  ! Date: 11/30/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P_sat
  PetscReal, intent(in) :: Pc
  PetscReal, intent(in) :: rho_kg
  PetscReal, intent(in) :: x_salt
  PetscReal, intent(out) :: P_vap

  PetscReal, parameter :: epsilon = 1.d-14
  PetscReal :: T_k, mw_mix

  T_k = T + 273.15d0 ! K
  mw_mix = x_salt*fmw_comp(3) + (1.d0-x_salt)*fmw_comp(1)

  if (Pc > epsilon) then
    ! P_vap = P_sat * exp(- fmw_comp(1) * (Pc**1.25d0) / &
    !          (rho_kg * IDEAL_GAS_CONSTANT * 1.d3 * T_k))
    P_vap = P_sat * exp(-1.d0 *fmw_comp(1) * (Pc ** 1.25d0) / &
            (rho_kg * IDEAL_GAS_CONSTANT * 1.d3 * T_k))
  else
    P_vap = P_sat
  endif

end subroutine SCO2VaporPressureBrine

! ************************************************************************** !


subroutine SCO2Equilibrate(T,P,p_co2,p_vap,p_sat,p_vap_brine, &
                           xco2g, xwg, xco2l, xsl, xwl, &
                           xmolco2g, xmolwg, xmolco2l, xmolsl, xmolwl, option)
  !
  ! Computes equilibrium partitioning between CO2 and water following 
  ! Spycher and Pruess, 2010
  !
  ! Author: Michael Nole
  ! Date: 11/29/23
  !

  use Option_module

  implicit none

  PetscReal, intent(in) :: T ! temperature (C)
  PetscReal, intent(in) :: P ! liquid or gas pressure (Pa)
  PetscReal, intent(out) :: p_co2 ! partial pressure of CO2 (Pa)
  PetscReal, intent(out) :: p_vap ! partial pressure of water (Pa)
  PetscReal, intent(in) :: p_sat ! saturated brine vapor pressure (Pa)
  PetscReal, intent(in) :: p_vap_brine ! reduced vapor pressure (Pa)
  PetscReal, intent(out) :: xco2g ! mass fraction of CO2 in gas phase
  PetscReal, intent(out) :: xwg ! mass fraction of water in gas phase
  PetscReal, intent(out) :: xco2l ! mass fraction of CO2 in liquid phase
  PetscReal, intent(inout) :: xsl ! mass fraction of salt in liquid phase
  PetscReal, intent(out) :: xwl ! mass fraction of water in liquid phase
  PetscReal, intent(out) :: xmolco2g ! mole fraction of CO2 in gas phase
  PetscReal, intent(out) :: xmolwg ! mole fraction of water in gas phase
  PetscReal, intent(out) :: xmolco2l ! mole fraction of CO2 in liquid phase
  PetscReal, intent(out) :: xmolsl ! mole fraction of salt in liquid phase
  PetscReal, intent(out) :: xmolwl ! mole fraction of water in liquid phase
  type(option_type) :: option

  PetscReal, parameter :: cac(2) = [7.54d7,-4.13d4]
  PetscReal, parameter :: caw(2) = [0.d0, 0.d0]
  PetscReal, parameter :: cacw(2) = [7.89d7,0.d0]
  PetscReal, parameter :: cbc = 27.8d0
  PetscReal, parameter :: CBW = 18.18d0
  PetscReal, parameter :: clkw(4) = [-2.209d0,3.097d-02,-1.098d-04,2.048d-07]
  PetscReal, parameter :: clkcn(4) = [1.169d0,1.368d-02,-5.380d-05,0.d0]
  PetscReal, parameter :: clkcg(4) = [1.189d0,1.304d-02,-5.446d-05,0.d0]
  PetscReal, parameter :: cvc = 32.6d0
  PetscReal, parameter :: cvw = 18.1d0
  PetscReal, parameter :: cam = 0.d0
  PetscReal, parameter :: cpr = 1.d0
  PetscReal, parameter :: clmb(3) = [2.217d-04,1.074d0,2.648d3]
  PetscReal, parameter :: cxi(3) = [1.3d-05,-2.012d1,5.259d3]
  PetscReal, parameter :: wtmna = 22.9898d0
  PetscReal, parameter :: wtmcl = 35.453d0
  PetscReal, parameter :: dac(2) = [8.008d7,-4.984d4]
  PetscReal, parameter :: daw(2) = [1.337d8,-1.4d4]
  PetscReal, parameter :: dbc = 28.25d0
  PetscReal, parameter :: dbw = 15.70d0
  PetscReal, parameter :: dkwc(2) = [1.427d-02,-4.037d-04]
  PetscReal, parameter :: dkcw(2) = [0.4228d0,-7.422d-04]
  PetscReal, parameter :: dlkw(5) = [-2.1077d0,2.8127d-02,-8.4298d-05, &
                                     1.4969d-07, -1.1812d-10]
  PetscReal, parameter :: dlkc(5) = [1.668d0,3.992d-03,-1.156d-05, &
                                     1.593d-09,0.d0]
  PetscReal, parameter :: dvmc(2) = [32.6d0,3.413d-02]
  PetscReal, parameter :: dvmw(2) = [18.1d0,3.137d-02]
  PetscReal, parameter :: dam(2) = [-3.084d-02,1.927d-05]
  PetscReal, parameter :: dpr(5) = [-1.9906d-01,2.0471d-03,1.0152d-04, &
                                    -1.4234d-6,1.4168d-08]
  PetscReal, parameter :: cmgw(17) = [3.2217d-03,1.225d-08,3.067d0, &
                                       -9.7198d-03,5.1621d0,3.485d2, &
                                       7.7053d1,1.0928d-02,3.663d2, &
                                       -1.9472d0,1.3937d0,2.4992d1, &
                                       2.5343d2,1.4677d1,3.7952d-02, &
                                       2.2122d3,-1.8936d0]
  PetscReal, parameter :: cmla(8) = [4.1265d-02,1.0715d1,2.6542d1, &
                                      2.8708d2,2.5478d-02,-3.0218d-04, &
                                      1.3776d-6,-2.2457d-09]


  PetscReal :: T_k, P_bar
  PetscReal :: T_bound(2)
  PetscReal :: nacl_param, cl_param
  PetscReal :: y0,a1,a2,tau1,tau2,Hc
  PetscReal :: a, b, coeff_a, coeff_b, coeff_c, coeff_d
  PetscReal :: xmol_na, xmol_cl
  PetscReal :: r1, r2, r3, w1, w2
  PetscReal :: vg, vn, v
  PetscReal :: a_mat(2,2)
  PetscReal :: b_vec(2), y_vec(2)
  PetscReal :: sum_fug
  PetscReal :: fugacity(2)
  PetscReal :: eqkw, eqkco2, pref
  PetscReal :: apc
  PetscReal :: fmw_gas, fmw_liq
  PetscReal :: pva

  PetscInt :: wid, co2_id, sid, lid, gid
  PetscInt :: i, k

  PetscReal, parameter :: epsilon = 1.d-14

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id
  lid = option%liquid_phase
  gid = option%gas_phase

  p_co2 = 0.d0
  xco2g = 0.d0
  xwg = 0.d0
  xco2l = 0.d0
  xwl = 0.d0
  xmolco2g = 0.d0
  xmolwg = 0.d0
  xmolco2l = 0.0d0
  xmolsl = 0.d0
  xmolwl = 0.d0


  T_k = T + 273.15d0
  P_bar = max(P,1.01325d5)*1.d-5

  T_bound(1) = 99.d0
  T_bound(2) = 101.d0 !109.d0

  T_bound = T_bound + 273.15d0

  ! Salinity offset
  nacl_param = cxi(1)*T_k + cxi(2)/T_k + cxi(3)/(T_k**2)
  cl_param = clmb(1)*T_k + clmb(2)/T_k + clmb(3)/(T_k**2)

  ! Nacl mass fraction to molality
  xmolsl = 1.d3*(xsl/(1.d0-xsl))/fmw_comp(3)
  xmol_na = xmolsl
  xmol_cl = xmolsl
  xmolwl = 1.d3 / fmw_comp(1)
  xmolsl = xmolsl / (xmolsl + xmolwl)
  apc = (1.d0 + (xmol_na + xmol_cl)/xmolwl)*exp(2.d0*cl_param*xmol_na + &
         nacl_param*xmol_cl*xmol_na)
  
  ! Simplified solution, up to 275 C
  if (sco2_spycher_simple) then
    ! Temperature in C
    y0 = cmgw(1) + cmgw(2)*(T**cmgw(3))
    a1 = cmgw(4) + cmgw(5)/(1.d0 + exp(-(T-cmgw(6))/cmgw(7)))
    tau1 = cmgw(8) + cmgw(9)*(T**cmgw(10))
    a2 = cmgw(11) + cmgw(12)/(1.d0 + exp(-(T-cmgw(13))/cmgw(14)))
    tau2 = cmgw(15) + cmgw(16)*(T**cmgw(17))

    ! mole fraction of water in the gas phase
    xmolwg = y0 + a1*exp(-tau1*P_bar) + a2*exp(-tau2*P_bar)
    xmolwg = max(min(xmolwg,1.d0),0.d0)

    Hc = 6.305d-4 * exp(2.4d3*((1.d0/T_k)-(1.d0/298.15d0)))
    pva = max(P-p_vap_brine,0.d0)/1.d5

    ! mole fraction of CO2 in the liquid phase
    xmolco2l = pva * Hc
    xmolco2l = max(min(xmolco2l,1.d0),0.d0)
    xmolsl = (1.d0 - xmolco2l)*xmolsl
    xmolwl = 1.d0 - xmolco2l - xmolsl

    xmolwg = xmolwg * p_vap / p_sat
    xmolco2g = 1.d0 - xmolwg
    

  elseif (T_k < T_bound(1)) then !Low temperature regime
    a = cac(1) + cac(2)*T_k
    b = cbc

    ! RKS EOS coefficients
    ! MAN: should compute gas density elsewhere and using Gas EOS routines?
    coeff_a = 1.d0 
    coeff_b = -1.d0 * IDEAL_GAS_CONSTANT * 1.d1 * T_k / P_bar
    coeff_c = -1.d0 * ((IDEAL_GAS_CONSTANT * 1.d1 * T_k * b / P_bar) - &
                        (a / (P_bar * sqrt(T_k))) + b**2)
    coeff_d = -1.d0 * (a*b/(P_bar*sqrt(T_k)))

    call CubicRootsNickalls(coeff_a, coeff_b, coeff_c, coeff_d, r1, r2, r3)

    vg = max(r1,r2,r3)
    vn = min(r1,r2,r3)

    if (T_k > CO2_CRITICAL_TEMPERATURE) then
      v = vg
    else
      w1 = P_bar*(vg-vn)
      w2 = IDEAL_GAS_CONSTANT * 1.d1 * T_k * log((vg-b)/(vn-b)) + &
           (a/(sqrt(T_k)*b)) * log((vg+b)*vn/(vn+b)*vg)
      if ((w2-w1)/epsilon >= epsilon) then
        !Gas Phase
        v = vg
      elseif ((w1-w2)/epsilon >= epsilon) then
        ! Liquid phase
        v = vn
      else
        ! Gas and liquid, use gas molar volume
        v = vg
      endif
    endif

    ! Fugacity coefficients
    a_mat(1,1) = caw(1)
    a_mat(2,1) = cacw(1)
    a_mat(1,2) = cacw(1)
    a_mat(2,2) = cac(1) + cac(2)*T_k
    b_vec(1) = cbw
    b_vec(2) = cbc
    y_vec(1) = 0.d0
    y_vec(2) = 1.d0
    do k = 1,2
      sum_fug = 0.d0
      do i = 1,2
        sum_fug = sum_fug + y_vec(i)*a_mat(i,k)
      enddo
      fugacity(k) = log(v/(v-b)) + (b_vec(k)/(v-b)) - &
          (2.d0*sum_fug/(IDEAL_GAS_CONSTANT*1.d1*(T_k**1.5d0)*b))* &
          log((v+b)/v) + &
          (a*b_vec(k)/(IDEAL_GAS_CONSTANT*1.d1*(T_k**1.5d0)*(b**2)))* &
          (log((v+b)/v) - b/(v+b)) - log(P_bar*v/(IDEAL_GAS_CONSTANT*1.d1*T_k))
      fugacity(k) = exp(fugacity(k))
    enddo

    ! A coefficient
    eqkw = 0.d0
    pref = cpr

    do i = 1,5
      eqkw = eqkw + dlkw(i)*(T**(i-1))
    enddo
    eqkw = (1.d1**eqkw)*exp((P_bar-pref)*cvw/(IDEAL_GAS_CONSTANT * 1.d1 * T_k))
    coeff_a = (eqkw/(fugacity(1)*P_bar))

    ! B coefficient
    eqkco2 = 0.d0
    pref = cpr
    if (T < 31.d0 .and. v < 94.d0) then
      do i = 1,4
        eqkco2 = eqkco2 + clkcn(i)*(T**(i-1))
      enddo
    else
      do i = 1,4
        eqkco2 = eqkco2 + clkcg(i)*(T**(i-1))
      enddo
    endif
    eqkco2 = (1.d1**eqkco2)* &
             exp((P_bar-pref)*cvc/(IDEAL_GAS_CONSTANT * 1.d1 * T_k))
    coeff_b = (fugacity(2)*P_bar) / (xmolwl * apc * eqkco2)

    ! Mole fractions
    xmolwg = (1.d0 - coeff_b) * xmolwl / (((1.d0/coeff_a)-coeff_b)* &
          (xmol_na + xmol_cl + xmolwl) + (xmol_na + xmol_cl)*coeff_b)
    xmolco2l = coeff_b * (1.d0 - xmolwg)
    xmolsl = (1.d0 - xmolco2l)*xmolsl
    xmolwl = 1.d0 - xmolco2l - xmolsl

    ! Vapor pressure lowering
    xmolwg = xmolwg * p_vap / p_sat
    xmolco2g = 1.d0 - xmolwg

  elseif (T_k > T_bound(2)) then ! High temperature regime
    ! Iterative solution 
    ! MAN: not yet implemented

    option%io_buffer = 'System temperature entered the high temperature &
                        &regime, T > 101 C. High temperature regime for &
                        &CO2 has not yet been implemented.'
    call PrintErrMsg(option)

  else ! Intermediate temperature regime
    ! Interpolate betweeen high and low temp regimes for smoothness
    ! MAN: not yet implemented

    option%io_buffer = 'System temperature entered the intermediate &
                        &temperature regime, 99-101 C. Intermediate &
                        &temperature regime for CO2 has not &
                        &yet been implemented.'
    call PrintErrMsg(option)

  endif

  ! Truncate mole fractions
  if (xmolwl < 1.d-16) xmolwl = 0.d0
  if (xmolco2l < 1.d-16) xmolco2l = 0.d0
  if (xmolsl < 1.d-16) xmolsl = 0.d0
  if (xmolwg < 1.d-16) xmolwg = 0.d0
  if (xmolco2g < 1.d-16) xmolco2g = 0.d0

  ! Component partial pressures
  p_co2 = xmolco2g * P
  p_vap = xmolwg * P

  ! Mass Fractions

  fmw_gas = xmolwg * fmw_comp(1) + xmolco2g * fmw_comp(2)
  fmw_liq = xmolwl * fmw_comp(1) + xmolco2l * fmw_comp(2) + xmolsl * fmw_comp(3)

  ! Gas Phase
  xwg = xmolwg * fmw_comp(1)/fmw_gas
  xco2g = xmolco2g * fmw_comp(2)/fmw_gas
  ! Liquid Phase
  xwl = xmolwl * fmw_comp(1) / fmw_liq
  xco2l = xmolco2l * fmw_comp(2) / fmw_liq
  ! MAN: This might cause problems since mass frac was used previously to
  !      compute intermediate variables.
  xsl = xmolsl * fmw_comp(3) / fmw_liq

  
end subroutine SCO2Equilibrate

! ************************************************************************** !

subroutine CubicRootsNickalls(a,b,c,d,r1,r2,r3)
  !
  ! Computes roots of a cubic polynomial following Nickalls, 1993, 
  ! A new approach to solving the cubic Cardans solution
  !
  ! Author: Michael Nole
  ! Date: 11/30/23
  !

  implicit none

  PetscReal, intent(in) :: a
  PetscReal, intent(in) :: b
  PetscReal, intent(in) :: c
  PetscReal, intent(in) :: d
  PetscReal, intent(out) :: r1
  PetscReal, intent(out) :: r2
  PetscReal, intent(out) :: r3

  PetscReal, parameter :: epsilon = 1.d-14
  PetscReal :: xn, yn, yn2, del, del2, h, h2, theta

  xn = -b/(3.d0*a)
  yn = a*(xn**3) + b*(xn**2) + c*xn + d
  yn2 = yn**2
  del2 = ((b**2) - (3.d0*a*c))/((3.d0*a)**2)

  if (del2 <= 0.d0) then
    h = 0.d0
  else
    del = sqrt(del2)
    h = -2.d0 * (del**3)
  endif

  h2 = 4.d0*(a**2)*(del2**3)

  if ((yn2 - h2) > epsilon) then

    r1 = (5.d-1/a)*(-yn + sqrt(yn2-h2))
    r2 = (5.d-1/a)*(-yn - sqrt(yn2-h2))
    r3 = xn + sign(dabs(r1)**(1.d0/3.d0),r1) + sign(dabs(r2)**(1.d0/3.d0),r2)
    r1 = r3
    r2 = r3
    
  elseif ((yn2-h2) < -epsilon) then

    theta = acos(yn/h)/3.d0
    r1 = xn + 2.d0 * del * cos(theta)
    r2 = xn + 2.d0 * del * cos(2.d0*PI/3.d0 + theta)
    r3 = xn + 2.d0 * del * cos(4.d0*PI/3.d0 + theta)

  else

    if (dabs(h)/epsilon > epsilon) then

      del = yn / (2.d0 *a)
      del = sign(dabs(del)**(1.d0/3.d0), del)
      r1 = xn + del
      r2 = r1
      r3 = xn - 2.d0*del

    else
      
      r1 = xn
      r2 = r1
      r3 = r2

    endif

  endif

end subroutine CubicRootsNickalls

! ************************************************************************** !

subroutine SCO2ComputeSaltDensity(T,P,rho_s)
  !
  ! Computes NaCl density following Battistelli et al., 1997
  !
  ! Author: Michael Nole
  ! Date: 11/30/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscReal, intent(out) :: rho_s

  rho_s = 2.165d3 * exp(-1.2d-4 * T + 4.d-11 * P)

end subroutine SCO2ComputeSaltDensity

! ************************************************************************** !

subroutine SCO2ComputeSurfaceTension(T,x_nacl,surface_tension)
  !
  ! Computes CO2-Water surface tension as a function of temperature
  ! and salt concentration following Abramzon and Gaukhberg, 1993.
  !
  ! Author: Michael Nole
  ! Date: 12/01/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: x_nacl
  PetscReal, intent(out) :: surface_tension

  PetscReal :: molality

  molality = 1.d3*x_nacl/(fmw_comp(3)*(1.d0-x_nacl))

  ! Pure water
  surface_tension = 1.d-3*(75.6592d0 - 1.40959d-1*T - 2.66317d-4*(T**2))
  ! With salt
  surface_tension = surface_tension + 1.57d-3*molality

  ! MAN: turn off surface tension effects
  surface_tension = CO2_REFERENCE_SURFACE_TENSION

end subroutine SCO2ComputeSurfaceTension

! ************************************************************************** !

subroutine SCO2WaterSaturationPressure(T,P_sat)
  !
  ! Computes pure water saturation pressure following Meyer et al., 1993
  !
  ! Author: Michael Nole
  ! Date: 12/03/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: P_sat

  PetscReal, parameter :: k(9) = [-7.691234564d0, -2.608023696d1, &
          -1.681706546d2, 6.423285504d1, -1.189646225d2, 4.167117320d0, &
          2.097506760d1, 1.d9, 6.d0]
  PetscReal, parameter :: T_c = 647.096d0
  PetscReal, parameter :: P_c = 22.064d6

  PetscReal :: T_r, T_rx
  PetscInt :: i

  T_r = (T + 273.15) / T_c
  T_rx = 1.d0 - T_r
  
  P_sat = 0.d0
  do i = 1,5
    P_sat = P_sat + k(i) * (T_rx ** i)
  enddo
  P_sat = P_sat / ((1.d0 + k(6) * T_rx + k(7) * (T_rx **2)) * T_r)
  P_sat = P_sat - T_rx / (k(8) * (T_rx ** 2) + k(9))
  P_sat = exp(P_sat) * P_c

end subroutine SCO2WaterSaturationPressure

! ************************************************************************** !

subroutine SCO2BrineSaturationPressure(T, x_salt, P_sat)
  !
  ! Computes brine saturation pressure following Haas, 1976
  !
  ! Author: Michael Nole
  ! Date: 12/05/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: x_salt
  PetscReal, intent(out) :: P_sat

  PetscReal, parameter :: s_a(3) = [5.93582d-6, -5.19386d-5, 1.23156d-5]
  PetscReal, parameter :: s_b(5) = [1.15420d-6, 1.41254d-7, -1.92476d-8, &
                                    -1.70717d-9, 1.05390d-10]

  PetscReal :: T_k, T_eq
  PetscReal :: x_salt_molal
  PetscReal :: a, b, c
  PetscInt :: i

  T_k = T + 273.15d0
  x_salt_molal = 1.d3 * x_salt / (fmw_comp(3) * (1.d0 - x_salt))

  a = 1.d0
  do i = 1,3
    a = a + s_a(i) * (x_salt_molal ** i)
  enddo

  b = 0.d0
  do i = 1,5
    b = b + s_b(i) * (x_salt_molal ** i)
  enddo

  c = 1.d0 / (a + b * T_k)

  T_eq = exp(c * log(T_k)) - 273.15d0

  call SCO2WaterSaturationPressure(T_eq, P_sat)

  end subroutine SCO2BrineSaturationPressure
! ************************************************************************** !

subroutine SCO2WaterSubregion(T,P,isubr)
  !
  ! Computes subregion of water EOS following Meyer et al., 1993
  !
  ! Author: Michael Nole
  ! Date: 12/03/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscInt, intent(out) :: isubr

  PetscReal, parameter :: L_coeff(3) = [1.574373327d1, -3.417061978d1, &
                                        1.931380707d1]
  PetscReal, parameter :: tol = 1.d-2
  PetscReal, parameter :: T_c = 647.096d0
  PetscReal, parameter :: P_c = 22.064d6 

  PetscReal :: T_k, P_sat, T_r

  T_k = T + 273.15

  if (T_k <= T_c) then
    call SCO2WaterSaturationPressure(T,P_sat)
    ! Subregion 5-6 and 1-4 boundary
    if (T <= 350.d0) then
      if ( (P-P_sat) >= tol) then
        isubr = 1
      elseif ( (P-P_sat) <= -tol) then
        isubr = 2
      else
        isubr = 6
      endif
    else
      if ((P-P_sat) >= tol) then
        isubr = 4
      elseif ((P-P_sat) <= -tol) then
        isubr = 2
      else
        isubr = 5
      endif
    endif
  else
    T_r = T_k / T_c
    P_sat = P_c * (L_coeff(1) + L_coeff(2)*T_r + L_coeff(3)*(T_r**2))
    if (P > P_sat) then
      isubr = 3
    else
      isubr = 2
    endif
  endif


end subroutine SCO2WaterSubregion

! ************************************************************************** !

subroutine SCO2WaterDensity(T,P,isubr,rho_l,rho_v,option)
  !
  ! Computes pure water density following Meyer et al., 1993
  !
  ! Author: Michael Nole
  ! Date: 12/03/23
  !

  use Option_module

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscInt, intent(in) :: isubr
  PetscReal, intent(out) :: rho_l
  PetscReal, intent(out) :: rho_v
  type(option_type) :: option

  PetscReal, parameter :: c_a(23) = &
                          [6.824687741d3, -5.422063673d2, -2.096666205d4, &
                           3.941286787d4, -6.733277739d4, 9.902381028d4, &
                           -1.093911774d5, 8.590841667d4, -4.511168742d4, &
                           1.418138926d4, -2.017271113d3, 7.982692717d0, &
                           -2.616571843d-2, 1.522411790d-3, 2.284279054d-2, &
                           2.421647003d2, 1.269716088d-10, 2.074838328d-7, &
                           2.174020350d-8, 1.105710498d-9, 1.293441934d1, &
                           1.308119072d-5, 6.047626338d-14]
  PetscReal, parameter :: s_a(12) = &
                          [8.438375405d-1, 5.362162162d-4, 1.720000000d0, &
                           7.342278489d-2, 4.975858870d-2, 6.537154300d-1, &
                           1.150000000d-6, 1.150800000d-5, 1.418800000d-1, &
                           7.002753165d0, 2.995284926d-4, 2.040000000d-1]
  PetscReal, parameter :: c_b(31) = &
                          [1.683599274d1, 2.856067796d1, -5.438923329d1, &
                           4.330662834d-1, -6.547711697d-1, 8.565182058d-2, &
                           6.670375918d-2, 1.388983801d0, 8.390104328d-2, &
                           2.614670893d-2, -3.373439453d-2, 4.520918904d-1, &
                           1.069036614d-1, -5.975336707d-1, -8.847535804d-2, &
                           5.958051609d-1, -5.159303373d-1, 2.075021122d-1, &
                           1.190610271d-1, -9.867174132d-2, 1.683998803d-1, &
                           -5.809438001d-2, 6.552390126d-3, 5.710218649d-4, &
                           1.936587558d2, -1.388522425d3, 4.126607219d3, &
                           -6.508211677d3, 5.745984054d3, -2.693088365d3, &
                           5.235718623d2]
  PetscReal, parameter :: s_b(5) = &
                          [7.633333333d-1, 4.006073948d-1, 8.636081627d-2, &
                          -8.532322921d-1, 3.460208861d-1]
  PetscReal, parameter :: L_coeff(3) = & 
                          [1.574373327d1, -3.417061978d1, 1.931380707d1]
  PetscReal, parameter :: s_l = 4.260321148d0
  PetscInt, parameter :: i_n(8) = [2, 3, 2, 2, 3, 2, 2, 2]
  PetscInt, parameter :: i_z(8,3) = reshape([13, 18, 18, 25, 32, 12, 24, 24, &
                                             3, 2, 10,14, 28, 11, 18, 14, 0, &
                                             1, 0, 0, 24, 0, 0, 0],shape(i_z))
  PetscInt, parameter :: i_t(8) = [0, 0, 0, 0, 0, 1, 1, 2]
  PetscInt, parameter :: i_x(8,2) = reshape([0, 0, 0, 0, 0, 14, 19, 54, 0, 0, &
                                             0, 0, 0, 0, 0, 27],shape(i_x))
  PetscReal, parameter :: epsilon = 1.d-14
  PetscReal, parameter :: T_c = 647.096d0
  PetscReal, parameter :: P_c = 22.064d6
  PetscReal, parameter :: v_c = 57.1075d0

  PetscReal :: T_r, P_r, beta_l
  PetscReal :: c_x, c_y, c_z, r_v, r_va, r_vb
  PetscInt :: i, j, indx, indx2

  ! MAN: if truncation is required, either do it elsewhere or throw error.
  ! if (T < 1.d-2 .or. T > 8.d2) then
  !   T = max(T,1.d-2)
  !   T = min(T,8.d2)
  ! endif
  ! if (P < 0.d0 .or. P > 1.d8) then
  !   P = max(P,0.d0)
  !   P = min(P,1.d8)
  ! endif

  if (P < epsilon) then
    rho_l = 0.d0
    rho_v = 0.d0
    return
  endif

  rho_l = 0.d0
  rho_v = 0.d0
  T_r = (T + 273.15)/T_c
  P_r = P/P_c

  beta_l = L_coeff(1) + L_coeff(2) * T_r + L_coeff(3) * (T_r ** 2)

  if (isubr == 1 .or. isubr == 6) then
    c_y = 1.d0 - s_a(1) * (T_r ** 2) - s_a(2) / (T_r ** 6)
    c_z = c_y + sqrt((s_a(3) * (c_y**2)) - (2.d0 * s_a(4) * T_r) + &
          (2.d0 * s_a(5) * P_r))
    r_v = c_a(12) * s_a(5) * (c_z **(-5.d0/17.d0)) + (c_a(13) + c_a(14) * &
          T_r + c_a(15) * (T_r **2) + c_a(16) * ((s_a(6)-T_r)**10) + c_a(17)/ &
          (s_a(7) + (T_r**19)))
    r_v = r_v - (c_a(18) + 2.d0*c_a(19)*P_r + 3.d0 * c_a(20) * &
          (P_r**2)) / (s_a(8) + (T_r **11)) - c_a(21) * (T_r **18) * &
          (s_a(9) + (T_r **2)) * (-3.d0 / ((s_a(10) + P_r) ** 4))
    r_v = r_v + 3.d0 * c_a(22) * (s_a(12) - T_r) * (P_r **2) + &
          4.d0 * c_a(23) * (P_r**3) / (T_r **20)
    rho_l = 1.d3 * fmw_comp(1) / (r_v * v_c)
  endif

  if (isubr == 2 .or. isubr == 6) then
    c_x = exp(s_b(1)*(1.d0-T_r))
    r_v = s_l * T_r / P_r
    indx = 6
    do i = 1,5
      r_va = 0.d0
      do j = 1,i_n(i)
        indx = indx + 1
        r_va = r_va + c_b(indx) * (c_x ** i_z(i,j))
      enddo
      r_v = r_v - i * (P_r**(i-1)) * r_va
    enddo
    indx = 18
    indx2 = 1
    do i = 6,8

      r_va = 0.d0
      do j = 1,i_n(i)
        indx = indx + 1
        r_va = r_va + c_b(indx) * (c_x ** i_z(i,j))
      enddo

      r_vb = 0.d0
      do j = 1,i_t(i)
        indx2 = indx2 + 1
        r_vb = r_vb + s_b(indx2) * (c_x ** i_x(i,j))
      enddo

      r_v = r_v - ((i-2) * (P_r ** (1-i)) * r_va) / &
            (((P_r ** (2-i)) + r_vb) ** 2)

    enddo

    indx = 24
    r_va = 0.d0
    do i = 0,6
      indx = indx + 1
      r_va = r_va + c_b(indx) * (c_x ** i)
    enddo
    r_v = r_v + 1.1d1 * ((P_R / beta_l) ** 10) * r_va
    rho_v = 1.d3 * fmw_comp(1) / (r_v * v_c)
  endif

  if (isubr == 3 .or. isubr == 5) then
    option%io_buffer = "Temperature range exceeded in Water EOS steam table: &
                        &entered subregion 3/5."
    call PrintErrMsg(option)
  endif

  if (isubr == 4 .or. isubr == 5) then
    option%io_buffer = "Temperature range exceeded in Water EOS steam table &
                        &entered subregion 4/5."
    call PrintErrMsg(option)
  endif


end subroutine SCO2WaterDensity

! ************************************************************************** !

subroutine SCO2BrineDensity(T, P, x_s, rho_b, option)
  !
  ! Computes brine density following Phillips et al., 1983
  !
  ! Author: Michael Nole
  ! Date: 12/04/23
  !

  use Option_module
  use EOS_Water_module

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscReal, intent(in) :: x_s ! Salt mass fraction
  PetscReal, intent(out) :: rho_b
  type(option_type) :: option

  PetscReal, parameter :: c_c(4) = [-3.033405D+0, 10.128163D+0, -8.750567D+0, &
                                    2.663107D+0]
  PetscReal, parameter :: c_h(10) = [-167.219D+0, 448.55D+0, -261.07D+0, &
                                     -13.644D+0, 13.97D+0, -0.315154D+0, &
                                     -1.203374D-3, 7.48908D-13, 0.1342489D+0, &
                                     -3.946963D-3]
  PetscReal, parameter :: s_c(3) = [-9.9559D+0, 7.0845D+0, 3.9093D+0]
  PetscReal, parameter :: s_a(3) = [-4.539D-3, -1.638D-4, 2.551D-5]
  PetscReal, parameter :: v_c = 3.1975D+0

  PetscReal :: P_bar, P_sat, P_w
  PetscReal :: rho_l, rho_v, spec_vol, x_s_molal
  PetscReal :: phi0, phi

  P_bar = P * 1.d-5
  x_s_molal = 1.d3 * x_s / (fmw_comp(3) * (1.d0 - x_s))

  call SCO2WaterSaturationPressure(T, P_sat)
  P_w = max(P,P_sat)
  call SCO2WaterDensity(T, P_w, ONE_INTEGER, rho_l, rho_v, option)

  rho_l = 1.d-3 * rho_l
  spec_vol = 1.d0 / rho_l

  phi0 = c_h(1) + c_h(2) * spec_vol + c_h(3) * (spec_vol ** 2)
  phi = phi0 + (c_h(4) + c_h(5) * spec_vol) * &
               ((spec_vol/(v_c - spec_vol)) ** 2) * sqrt(x_s_molal)

  rho_b = (1.d3 + x_s_molal * fmw_comp(3)) / &
          (1.d3 * spec_vol + x_s_molal * phi)

  ! kg/m^3
  rho_b = rho_b * 1.d3

end subroutine SCO2BrineDensity

! ************************************************************************** !

subroutine SCO2DensityCompositeLiquid(T,rho_b,x_co2, rho_l)
  !
  ! Computes density of the liquid phase as a funtion of brine density and
  ! CO2 concentration, Alendal and Drange, 2001
  !
  ! Author: Michael Nole
  ! Date: 12/04/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: rho_b
  PetscReal, intent(in) :: x_co2
  PetscReal, intent(out) :: rho_l

  PetscReal, parameter :: pv_coeff(5) = [37.36d-3, -7.109d-5, -3.812d-8, &
                                         3.296d-9, -3.702d-12]

  PetscReal :: pv_co2, c_co2
  PetscInt :: i

  pv_co2 = 0.d0
  do i = 1,5
    pv_co2 = pv_co2 + pv_coeff(i) * (T ** (i-1))
  enddo

  c_co2 = pv_co2 * rho_b * x_co2 / fmw_comp(2)

  rho_l = rho_b / (1.d0 + c_co2 - x_co2)


end subroutine SCO2DensityCompositeLiquid

! ************************************************************************** !

subroutine SCO2ViscosityWater(T, P, rho_w, visc, option)
  !
  ! Computes viscosity of pure water following Meyer et al., 1993
  !
  ! Author: Michael Nole
  ! Date: 12/03/23
  !

  use Option_module

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscReal, intent(in) :: rho_w
  PetscReal, intent(out) :: visc
  type(option_type) :: option

  PetscReal, parameter :: T_ref = 647.27d0
  PetscReal, parameter :: rho_ref = 317.763d0
  PetscReal, parameter :: P_ref = 2.2115d7
  PetscReal, parameter :: visc_ref = 5.5071d1
  PetscReal, parameter :: coeff(46) = [ 1.d0, 9.78197d-1, 5.79829d-1, &
          -2.02354d-1, 5.132047d-1, 3.205656d-1, 0.d0, 0.d0, -7.782567d-1, &
          1.885447d-1, 2.151778d-1, 7.317883d-1, 1.241044d0, 1.476783d0, &
          0.d0, 0.d0, -2.818107d-1, -1.070786d0, -1.263184d0, 0.d0, &
          0.d0, 0.d0, 1.778064d-1, 4.605040d-1, 2.340379d-1, -4.924179d-1, &
          0.d0, 0.d0, -4.176610d-2, 0.d0, 0.d0, 1.600435d-1, 0.d0, &
          0.d0, 0.d0, -1.578386d-2, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
          0.d0, 0.d0, -3.629481d-3,  0.d0, 0.d0 ]

  PetscReal :: T_r
  PetscReal :: rho_r, rho_r2, chi
  PetscReal :: P_r, dP, dP_r, P_inc
  PetscReal :: visc_a, rho_l, rho_vap
  PetscInt :: i,j,ix,isubr

  T_r = (T + 273.15) / T_ref
  rho_r = rho_w / rho_ref
  P_r = P / P_ref

  visc = 0.d0
  do i = 0,3
    visc = visc + coeff(i+1) / (T_r**i)
  enddo
  visc = sqrt(T_r)/visc

  visc_a = 0.d0
  ix = 0
  do i = 0,5
    do j = 0,6
      ix = 6 * j + i + 5
      visc_a = visc_a + coeff(ix) * (((1.d0/T_r) - 1.d0) ** i) * &
               ((rho_r - 1.d0) ** j)
    enddo
  enddo
  visc = visc * exp(rho_r*visc_a)

  if (T_r >= 0.997 .and. T_r <= 1.0082 .and. &
      rho_r >= 0.755 .and. rho_r <= 1.290) then
    dP = 1.d-1
    dP_r = dP / P_ref
    P_inc = P + dP
    call SCO2WaterSubregion(T,P_inc,isubr)
    call SCO2WaterDensity(T,P_inc,isubr,rho_l,rho_vap,option)
    if ((1.d0 - dabs(rho_l / rho_w)) < (1.d0 - dabs(rho_vap/rho_w))) then
      rho_r2 = rho_l / rho_ref
    else
      rho_r2 = rho_vap / rho_ref
    endif
    chi = rho_r * (rho_r2 - rho_r) / dP_r
    if (chi >= 21.93d0) visc = visc * 0.922d0 * (chi ** 0.0263d0)
  endif

  visc = 1.d-6 * visc * visc_ref

end subroutine SCO2ViscosityWater

! ************************************************************************** !

subroutine SCO2ViscosityCO2(T, rho_co2, visc)
  !
  ! Computes viscosity of CO2 following Fenghour et al., 1998
  !
  ! Author: Michael Nole
  ! Date: 12/03/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: rho_co2
  PetscReal, intent(out) :: visc

  PetscReal, parameter :: s_a(5) = [0.235156d0, -0.491266d0, 5.211155d-2, &
                                    5.347906d-2, -1.537102d-2]
  PetscReal, parameter :: s_b(5) = [0.4071119d-2, 0.7198037d-4, &
                                    0.2411697d-16, 0.2971072d-22, &
                                    -0.1627888d-22]
  PetscReal, parameter :: T_ref = 251.196d0

  PetscReal :: T_k, T_r
  PetscReal :: ecs, visc_0, visc_ex
  PetscInt :: i

  T_k = T + 273.15
  T_r = T_k / T_ref

  ecs = 0.d0
  do i = 0,4
    ecs = ecs + s_a(i+1)*(log(T_r) ** i)
  enddo
  ecs = exp(ecs)
  visc_0 = 1.00697d0 * sqrt(T_k)/ecs

  visc_ex = s_b(1) * rho_co2 + s_b(2) * (rho_co2 ** 2) + &
           s_b(3) * (rho_co2 ** 6) / (T_r ** 3) + &
           s_b(4) * (rho_co2 ** 8) + s_b(5) * (rho_co2 **8) / T_r
           
  visc = (visc_0 + visc_ex) * 1.d-6

end subroutine SCO2ViscosityCO2

! ************************************************************************** !

subroutine SCO2ViscosityGas(visc_w, visc_co2, xwg, xco2g, visc)
  !
  ! Computes gas phase viscosity following Reid et al., 1987
  !
  ! Author: Michael Nole
  ! Date: 12/03/23
  !

  implicit none

  PetscReal, intent(in) :: visc_w
  PetscReal, intent(in) :: visc_co2
  PetscReal, intent(in) :: xwg ! mole fraction
  PetscReal, intent(in) :: xco2g ! mole fraction
  PetscReal, intent(out) :: visc

  PetscReal :: phi_w, phi_co2, chi_w, chi_co2

  phi_w = ((1.d0 + sqrt(visc_co2/visc_w) * &
        ((fmw_comp(1)/fmw_comp(2)) ** 2.5d-1)) **2) / &
        sqrt(8.d0 * (1.d0 + fmw_comp(2)/fmw_comp(1)))
  phi_co2 = ((1.d0 + sqrt(visc_w/visc_co2) * &
        ((fmw_comp(2)/fmw_comp(1)) ** 2.5d-1)) **2) / &
        sqrt(8.d0 * (1.d0 + fmw_comp(1)/fmw_comp(2)))
  chi_w = xwg + xco2g * phi_co2
  chi_co2 = xwg * phi_w + xco2g
  visc = xwg * visc_w / chi_w + xco2g * visc_co2 / chi_co2

end subroutine SCO2ViscosityGas

! ************************************************************************** !

subroutine SCO2ViscosityBrine(T, x_salt, visc_w, visc_b)
  !
  ! Computes brine viscosity following Phillips et al., 1981
  !
  ! Author: Michael Nole
  ! Date: 12/04/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: x_salt
  PetscReal, intent(in) :: visc_w
  PetscReal, intent(out) :: visc_b

  PetscReal, parameter :: s_a(5) = [0.0816d0, 0.0122d0, 0.000128d0, &
                                    0.000629d0,-0.7d0]
  PetscReal :: x_salt_molal

  x_salt_molal = 1.d3 * x_salt / (fmw_comp(3) * (1.d0 - x_salt))

  visc_b = visc_w * (1.d0 + s_a(1) * x_salt_molal + s_a(2) * &
           (x_salt_molal ** 2) + s_a(3) * (x_salt_molal ** 3) + &
           s_a(4) * T * (1.d0 - exp(s_a(5) * x_salt_molal)))

end subroutine SCO2ViscosityBrine

! ************************************************************************** !

subroutine SCO2ViscosityLiquid(x_co2, visc_b, visc_co2, visc_l)
  !
  ! Computes composite liquid phase viscosity, Kimagai and Yokoyama, 1999
  !
  ! Author: Michael Nole
  ! Date: 12/04/23
  !

  implicit none

  PetscReal, intent(in) :: x_co2
  PetscReal, intent(in) :: visc_b
  PetscReal, intent(in) :: visc_co2
  PetscReal, intent(out) :: visc_l

  PetscReal, parameter :: epsilon = 1.d-14

  visc_l = (1.d0 - x_co2) * log(visc_b)
  if (visc_co2 > 1.d-14) visc_l = visc_l + x_co2 * log(visc_co2)
  visc_l = exp(visc_l)

end subroutine SCO2ViscosityLiquid

! ************************************************************************** !

function SCO2EnthalpyCompositeLiquid(T, x_salt, x_co2, h_brine, h_co2)
  !
  ! Computes composite liquid phase enthalpy, Battistelli et al., 1997
  !
  ! Author: Michael Nole
  ! Date: 12/08/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: x_salt ! mass fraction
  PetscReal, intent(in) :: x_co2 ! mass fraction
  PetscReal, intent(in) :: h_brine ! J/kg
  PetscReal, intent(in) :: h_co2 ! J/kg
  PetscReal :: SCO2EnthalpyCompositeLiquid

  PetscReal :: dT, Hc, T_pert, Hc_pert, dHc, T_k, h_sol

  dT = 1.d-6
  Hc = SCO2Henry(T, x_salt)
  T_pert = T + dT
  Hc_pert = SCO2Henry(T_pert, x_salt)
  dHc = log(Hc_pert / Hc) / dT

  T_k = T + 273.15d0
  h_sol = -IDEAL_GAS_CONSTANT * 1.d3 * (T_k **2) * dHc / fmw_comp(2)

  ! J/kg
  SCO2EnthalpyCompositeLiquid = max(1.d0-x_co2,0.d0) * h_brine + &
                                x_co2 * (h_co2 + h_sol)

end function SCO2EnthalpyCompositeLiquid

! ************************************************************************** !

function SCO2Henry(T, x_salt)
  !
  ! Computes Henry's coefficient for CO2 in brine, Battistelli et al., 1997
  ! Eq. 29
  !
  ! Author: Michael Nole
  ! Date: 12/08/23
  !

  implicit none

  PetscReal, intent(in) :: T ! C
  PetscReal, intent(in) :: x_salt ! mass fraction
  PetscReal :: SCO2Henry ! Henrys constant

  PetscReal, parameter :: c_b(6) = [7.83666d7, 1.96025d6, 8.20574d4, &
                                    -7.40674d2, 2.18380d0, -2.20999d-3]
  PetscReal, parameter :: c_C(5) = [1.19784d-1, -7.17823d-4, 4.93854d-6, &
                                    -1.03826d-8,1.08233d-11]
  PetscInt :: i
  PetscReal :: skb
  PetscReal :: Hc
                                    
  Hc = 0.d0
  do i = 0,5
    Hc = Hc + c_b(i+1) * (T ** i)
  enddo

  skb = 0.d0
  do i = 0,4
    skb = skb + c_c(i+1) * (T ** i)
  enddo

  SCO2Henry = Hc * &
              (1.d1 ** (1.d3 * x_salt / (fmw_comp(3) * (1.d0 - x_salt)) * skb))
  

end function SCO2Henry

! ************************************************************************** !

subroutine SCO2DiffusionCoeff(T,P,xsl,viscl,sco2_parameter,option)
  !
  ! Computes CO2-Water diffusion coefficient, Cadogan et al., 2014
  !
  ! Author: Michael Nole
  ! Date: 12/04/23
  !

  use Option_module

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscReal, intent(in) :: xsl
  PetscReal, intent(in) :: viscl
  type(sco2_parameter_type) :: sco2_parameter
  type(option_type) :: option

  PetscReal, parameter :: c_a(8) = [1.06036d0,1.5610d-1,1.9300d-1,4.7635d-1, &
                                   1.03587d0,1.52996d0,1.76474d0,3.89411d0]
  PetscReal, parameter :: c_w(2) = [3.190008977d0, 429.18d0]
  PetscReal, parameter :: c_co2(2) = [3.795165630d0, 95.85245d0]
  PetscReal, parameter :: epsilon = 1.d-14
  PetscReal, parameter :: viscw_ref = 0.8904339807d-3
  PetscReal, parameter :: T_ref = 25.d0

  PetscReal :: T_k , P_bar, s_molal
  PetscReal :: eps, sig, T_r, omega, w_mix
  PetscReal :: Dco2l, Dwg, Dnacl
  PetscInt :: lid, gid, co2_id, sid, wid
  PetscReal :: dlng, viscb, viscbr

  lid = option%liquid_phase
  gid = option%gas_phase

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id

  ! CO2 diffusion through the gas phase
  T_k = T + 273.15
  P_bar = P * 1.d-5

  eps = sqrt(c_w(2)*c_co2(2))
  sig = 5.d-1 * (c_w(1) + c_co2(1))
  T_r = T_k / eps
  omega = (c_a(1) / (T_r ** c_a(2))) + (c_a(3) / (exp(c_a(4) * T_r))) + &
          (c_a(5) / (exp(c_a(6) * T_r))) + (c_a(7) / (exp(c_a(8) * T_r)))
  w_mix = 2.d0 / ((1.d0 / fmw_comp(2)) + (1.d0 / fmw_comp(1)))
  Dwg = (3.03d0 - (9.8d-1 / sqrt(w_mix))) * 1.d-3 * &
                    (T_k ** 1.5d0) / (P_bar * sqrt(w_mix) * (sig ** 2) * &
                     omega) * 1.d-4

  ! CO2 diffusion through the liquid phase
  Dco2l = 3.5984d0 - 6.5113d-2*T_k + 2.0282D-4*(T_k**2)
  Dco2l = Dco2l*1.D-9

  ! Correct for NaCl
  Dco2l = Dco2l*(1.6678d0 - 1.2531d-1*(1.d3*(xsl/fmw_comp(THREE_INTEGER)) / &
         (1.d0-xsl))) / 1.6678d0

  ! Salt diffusion through the liquid phase
  s_molal = 1.d3 * xsl / (fmw_comp(THREE_INTEGER) * (1.d0 - xsl))

  if (s_molal > epsilon) then
    dlng = (-0.2555d0/(sqrt(s_molal)*(1.d0+sqrt(s_molal))) + &
         0.2555d0/((1.d0+sqrt(s_molal))**2) + &
         (6.d-2 + 6.d-1*0.0547d0)/((1.d0+1.5d0*s_molal)**2) - &
         3.d0*(6.d-2 + 6.d-1*0.0547d0)*s_molal/((1.d0+1.5d0*s_molal)**3) + &
         0.0547d0)*2.302585d0
  else
    dlng = 0.d0
  endif

  Dnacl = 2.254d-9
  call SCO2ViscosityBrine(T_ref,xsl,viscw_ref,viscbr)
  Dnacl = Dnacl * (viscw_ref/viscbr) * (1.d0 + s_molal*(dlng))
  call SCO2ViscosityBrine(T_k,xsl,viscl,viscb)
  Dnacl = Dnacl * (T_k/2.9815d2)*(viscbr/viscb)

  sco2_parameter%diffusion_coefficient(co2_id,lid) = Dco2l
  sco2_parameter%diffusion_coefficient(wid,gid) = Dwg
  sco2_parameter%diffusion_coefficient(sid,lid) = Dnacl

end subroutine SCO2DiffusionCoeff

! ************************************************************************** !

subroutine SCO2ScalePermPhi(sco2_auxvar, material_auxvar, global_auxvar, &
                            option)
  !
  ! Computes effective permeability and porosity as a function of precipitate
  ! saturations (Verma and Pruess, 1998)
  !
  ! Author: Michael Nole
  ! Date: 12/04/23
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(sco2_auxvar_type) :: sco2_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(global_auxvar_type) :: global_auxvar ! Stores mineral info
  type(option_type) :: option

  PetscReal :: phi_r ! zero-permeability limit fraction of porosity
  PetscReal :: f ! geometric factor
  PetscReal :: phi_0 ! initial fraction of porosity
  PetscReal :: theta
  PetscReal :: tao, omega
  PetscInt :: pid

  pid = option%precipitate_phase

  !MAN: hard-code for now
  tao = 1.5d0
  phi_0 = material_auxvar%porosity_base
  phi_r = 8.d-1
  f = phi_r

  sco2_auxvar%effective_porosity = max(sco2_auxvar%effective_porosity * &
                                      (1.d0 - sco2_auxvar%sat(pid)), &
                                      sco2_auxvar%effective_porosity * phi_r, &
                                      1.d-12)

  select case(permeability_reduction_model)

    case(ONE_INTEGER)
      ! Simplified Verma & Pruess
      sco2_auxvar%effective_permeability = ((sco2_auxvar%effective_porosity / &
                                             phi_0 - phi_r ) / &
                                             (1.d0 - phi_r )) ** tao
    case(TWO_INTEGER)
      ! Verma & Pruess model
      theta = max((1.d0 - sco2_auxvar%sat(pid) - phi_r) / (1.d0 - phi_r) , &
                  0.d0)
      omega = 1.d0 + (1.d0/f)/(1.d0/phi_r - 1.d0)
      sco2_auxvar%effective_permeability = &
                          (theta ** 2) * (1.d0 - f +  f/(omega**2)) / &
                          (1.d0 - f + f * (theta / (theta + omega - 1.d0)) ** 2)

  end select

end subroutine SCO2ScalePermPhi

! ************************************************************************** !

subroutine SCO2Tortuosity(s_l, s_g, phi, tao_l, tao_g)
  !
  ! Computes tortuosity of the rock to gas and liquid phases
  !
  ! Author: Michael Nole
  ! Date: 12/04/23
  !

  implicit none

  PetscReal, intent(in) :: s_l
  PetscReal, intent(in) :: s_g
  PetscReal, intent(in) :: phi
  PetscReal, intent(out) :: tao_l
  PetscReal, intent(out) :: tao_g

  PetscReal, parameter :: epsilon = 1.d-14

  ! Right now, just Millington & Quirk

  if (s_l * phi < epsilon) then
    tao_l = 0.d0
  else
    tao_l = (phi * (s_l ** 7)) ** (1.d0 / 3.d0)
  endif
  if (s_g * phi < epsilon) then
    tao_g = 0.d0
  else
    tao_g = (phi * (s_g ** 7)) ** (1.d0 / 3.d0)
  endif

end subroutine SCO2Tortuosity

! ************************************************************************** !

subroutine SCO2ComputeSatHysteresis(characteristic_curves, Pc, Sl_min, &
                                    beta_gl,rho_l, Sl, Sgt, option)
  !
  ! Compute saturation as a function of Pc including hysteretic effects.
  ! I believe this only works for monontonic Pc functions.
  !
  ! Author: Michael Nole
  ! Date: 12/07/23
  !

  use Option_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module

  implicit none

  class(characteristic_curves_type), intent(in) :: characteristic_curves
  PetscReal, intent(in) :: Pc
  PetscReal, intent(in) :: Sl_min
  PetscReal, intent(in) :: beta_gl
  PetscReal, intent(in) :: rho_l
  PetscReal, intent(out) :: Sl
  PetscReal, intent(inout) :: Sgt
  type(option_type) :: option

  PetscReal, parameter :: gravity = EARTH_GRAVITY

  PetscReal :: dsat_dpres
  PetscReal :: R
  PetscReal :: Sgt_max_bar, Sl_min_bar, Sgr_bar, Sgt_bar, Se
  PetscReal :: Sgt_max, Srl
  PetscReal :: capillary_head

  Sgt_max = characteristic_curves%saturation_function%sgt_max
  Srl = characteristic_curves%saturation_function%Sr

  select type(sf => characteristic_curves%saturation_function)
    class is (sat_func_VG_STOMP_type)
      capillary_head = max(beta_gl * Pc / &
                       (LIQUID_REFERENCE_DENSITY * gravity),1.d-14)
      call characteristic_curves%saturation_function% &
                    Saturation(capillary_head,Sl,dsat_dpres,option)
    class default
      call characteristic_curves%saturation_function% &
                    Saturation(Pc,Sl,dsat_dpres,option)
  end select
  
  if (Uninitialized(Sgt_max)) return

  ! Check if we're on a scanning path
  if (Sl > Sl_min .and. Sgt_max > 0.d0) then
    ! Hysteresis adjustment: Update Sl and Sgt
    Sgt_max_bar = Sgt_max / (1.d0 - Srl)
    Sl_min_bar = (Sl_min - Srl) / (1.d0 - Srl)
    R = 1.d0 / (Sgt_max_bar) - 1.d0
    Sgr_bar = (1.d0 - Sl_min_bar) / (1.d0 + R * (1.d0 - Sl_min_bar))

    Se = (Sl - Srl) / (1.d0 - Srl)
    Sgt_bar = Sgr_bar - (1.d0 - Se) / &
              (1.d0 + R * (1.d0 - Se))

    Sgt = Sgt_bar * (1.d0 - Srl)
    Sl = Sl - Sgt
  else
    Sgt = 0.d0
  endif

  if (Sgt > Sgt_max) Sgt = Sgt_max

end subroutine SCO2ComputeSatHysteresis

! ************************************************************************** !

subroutine SCO2ComputePcHysteresis(characteristic_curves, Sl, Sgt, beta_gl, &
                                   Pc, option)
  !
  ! Compute Pc as a function of saturation including trapped gas.
  ! I believe this only works for monontonic Pc functions.
  !
  ! Author: Michael Nole
  ! Date: 12/07/23
  !

use Option_module
use Characteristic_Curves_module
use Characteristic_Curves_Common_module

implicit none

class(characteristic_curves_type), intent(in) :: characteristic_curves
PetscReal, intent(in) :: Sl
PetscReal, intent(in) :: Sgt
PetscReal, intent(in) :: beta_gl
PetscReal, intent(out) :: Pc
type(option_type) :: option

PetscReal, parameter :: gravity = EARTH_GRAVITY
PetscReal :: dpc_dsatl
PetscReal :: Sl_eff

! Add trapped gas to the liquid saturation before sending into the Pc function
Sl_eff = Sl + Sgt

call characteristic_curves%saturation_function%CapillaryPressure(Sl_eff, Pc, &
                                                            dpc_dsatl,option)
select type(sf => characteristic_curves%saturation_function)
    class is (sat_func_VG_STOMP_type)
      ! Pc is the capillary head
      Pc = Pc * LIQUID_REFERENCE_DENSITY * gravity / beta_gl
    class default
      
end select


end subroutine SCO2ComputePcHysteresis

! ************************************************************************** !

subroutine SCO2ComputeEffectiveDiffusion(sco2_parameter, sco2_auxvar, option)
  !
  ! Compute effective diffusion coefficients
  !
  ! Author: Michael Nole
  ! Date: 12/08/23
  !

  use Option_module

  implicit none

  type(sco2_parameter_type) :: sco2_parameter
  type(sco2_auxvar_type) :: sco2_auxvar
  type(option_type) :: option

  PetscInt :: lid, gid, wid, co2_id, sid, tgid
  PetscReal :: T_scaled

  lid = option%liquid_phase
  gid = option%gas_phase
  tgid = option%trapped_gas_phase

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id

  select case(sco2_diffusion_model)

    case(ONE_INTEGER)

      ! Salt effective_diffusion_coeff in liquid
      T_scaled = (sco2_auxvar%temp + 273.15d0) / SALT_REFERENCE_TEMPERATURE 
      sco2_auxvar%effective_diffusion_coeff(sid,lid) = &
                 sco2_parameter%diffusion_coefficient(sid,lid) * T_scaled * &
                (LIQUID_REFERENCE_VISCOSITY / sco2_auxvar%visc(lid)) * &
                sco2_auxvar%tortuosity(lid) * sco2_auxvar%sat(lid) * &
                sco2_auxvar%effective_porosity

    case(TWO_INTEGER)
      ! MAN: not complete

    case(THREE_INTEGER)
      ! Salt effective_diffusion_coeff in liquid
      sco2_auxvar%effective_diffusion_coeff(sid,lid) = &
                 sco2_auxvar%tortuosity(lid) * &
                 sco2_auxvar%sat(lid) * sco2_auxvar%effective_porosity * &
                 sco2_parameter%diffusion_coefficient(sid,lid)
                 

  end select

  sco2_auxvar%effective_diffusion_coeff(co2_id,lid) = &
                 sco2_auxvar%tortuosity(lid) * &
                 sco2_auxvar%sat(lid) * sco2_auxvar%effective_porosity * &
                 sco2_parameter%diffusion_coefficient(co2_id,lid)

  sco2_auxvar%effective_diffusion_coeff(wid,gid) = &
                 sco2_auxvar%tortuosity(gid) * &
                 (sco2_auxvar%sat(gid) - sco2_auxvar%sat(tgid)) * &
                 sco2_auxvar%effective_porosity * &
                 sco2_parameter%diffusion_coefficient(wid,gid)

  sco2_auxvar%effective_diffusion_coeff(co2_id,gid) = &
                               sco2_auxvar%effective_diffusion_coeff(wid,gid)

end subroutine SCO2ComputeEffectiveDiffusion

! ************************************************************************** !

subroutine SCO2SaltEnthalpy(T,H)

  ! Calculate internal energy of solid salt (Lide et al., 1994)
  !
  ! Author: Michael Nole
  ! Date 12/13/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: H

  PetscReal, parameter :: c_a(5) = [25.19d0, 0.1973d0, -6.0114d-4, &
                                    8.81505d-7,-4.765d-10]
  PetscReal :: T_k
  PetscInt :: i

  T_k = T + 273.15d0

  H = -1.24858d-4
  
  do i = 1,5
    H = H + c_a(i) * (T_k ** i) / i
  enddo

  ! J/kg
  H = 1.d3 * H / fmw_comp(THREE_INTEGER)

end subroutine SCO2SaltEnthalpy

! ************************************************************************** !

subroutine SCO2ComputeSaltSolubility(T, x_salt)
  !
  ! Computes solubility of NaCl in water. McKibbin and McNabb, 1993.
  !
  ! Author: Michael Nole
  ! Date: 12/17/23
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: x_salt

  PetscReal, parameter :: coeff(3) = [2.6218d-1, 7.2d-5, 1.06d-6]

  x_salt = coeff(1) + coeff(2) * T + coeff(3) * T ** 2

end subroutine SCO2ComputeSaltSolubility

! ************************************************************************** !

subroutine SCO2OutputAuxVars1(sco2_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,string,append,option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Michael Nole
  ! Date: 12/15/23
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(sco2_auxvar_type) :: sco2_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  PetscBool :: append
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string2
  PetscInt :: apid, cpid, vpid, spid, tgid
  PetscInt :: gid, lid, acid, wid, sid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density, salt_density
  PetscReal :: liquid_energy, gas_energy, salt_energy
  PetscReal :: liquid_saturation, gas_saturation, &
               trapped_gas_saturation, precipitate_saturation

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id
  tgid = option%trapped_gas_phase

  acid = option%air_id ! air component id
  wid = option%water_id
  sid = option%salt_id

  liquid_density = 0.d0
  gas_density = 0.d0
  salt_density = 0.d0
  liquid_energy = 0.d0
  gas_energy = 0.d0
  salt_energy = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0
  trapped_gas_saturation = 0.d0
  precipitate_saturation = 0.d0

  write(string2,*) natural_id
  string2 = trim(adjustl(string)) // '_' // trim(adjustl(string2)) // '.txt'
  if (append) then
    open(unit=86,file=string2,position='append')
  else
    open(unit=86,file=string2)
  endif

  write(86,*) '--------------------------------------------------------'
  write(86,*) trim(string)
  write(86,*) '             cell id: ', natural_id
  select case(global_auxvar%istate)
    case(SCO2_LIQUID_STATE)
      write(86,*) ' Thermodynamic state: Liquid phase'
      liquid_density = sco2_auxvar%den(lid)
      liquid_energy = sco2_auxvar%U(lid)
      liquid_saturation = sco2_auxvar%sat(lid)
    case(SCO2_GAS_STATE)
      write(86,*) ' Thermodynamic state: Gas phase'
      gas_density = sco2_auxvar%den(gid)
      gas_energy = sco2_auxvar%U(gid)
      gas_saturation = sco2_auxvar%sat(gid)
    case(SCO2_LIQUID_GAS_STATE)
      write(86,*) ' Thermodynamic state: Liquid-Gas phase'
      liquid_density = sco2_auxvar%den(lid)
      gas_density = sco2_auxvar%den(gid)
      liquid_energy = sco2_auxvar%U(lid)
      gas_energy = sco2_auxvar%U(gid)
      liquid_saturation = sco2_auxvar%sat(lid)
      gas_saturation = sco2_auxvar%sat(gid)
    case(SCO2_TRAPPED_GAS_STATE)
      write(86,*) ' Thermodynamic state: Trapped Gas phase'

  end select
  liquid_mass = (liquid_density*sco2_auxvar%xmol(lid,lid)* &
                 liquid_saturation+ &
                 gas_density*sco2_auxvar%xmol(lid,gid)* &
                 gas_saturation)* &
                 sco2_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*sco2_auxvar%xmol(gid,lid)* &
              liquid_saturation+ &
              gas_density*sco2_auxvar%xmol(gid,gid)* &
              gas_saturation)* &
              sco2_auxvar%effective_porosity*material_auxvar%volume
  write(86,*) 'tot liq comp mass [kmol]: ', liquid_mass
  write(86,*) 'tot gas comp mass [kmol]: ', gas_mass
  write(86,*) '             energy [MJ]: ', liquid_mass*liquid_energy + &
                                            gas_mass*gas_energy
  write(86,*) '         liquid pressure: ', sco2_auxvar%pres(lid)
  write(86,*) '            gas pressure: ', sco2_auxvar%pres(gid)
  write(86,*) '            air pressure: ', sco2_auxvar%pres(apid)
  write(86,*) '      capillary pressure: ', sco2_auxvar%pres(cpid)
  write(86,*) '          vapor pressure: ', sco2_auxvar%pres(vpid)
  write(86,*) '     saturation pressure: ', sco2_auxvar%pres(spid)
  write(86,*) '         temperature [C]: ', sco2_auxvar%temp
  write(86,*) '       liquid saturation: ', sco2_auxvar%sat(lid)
  write(86,*) '          gas saturation: ', sco2_auxvar%sat(gid)
  write(86,*) '   liquid density [kmol]: ', sco2_auxvar%den(lid)
  write(86,*) '     liquid density [kg]: ', sco2_auxvar%den_kg(lid)
  write(86,*) '      gas density [kmol]: ', sco2_auxvar%den(gid)
  write(86,*) '        gas density [kg]: ', sco2_auxvar%den_kg(gid)
  write(86,*) '     X (water in liquid): ', sco2_auxvar%xmol(lid,lid)
  write(86,*) '       X (air in liquid): ', sco2_auxvar%xmol(gid,lid)
  write(86,*) '        X (water in gas): ', sco2_auxvar%xmol(lid,gid)
  write(86,*) '          X (air in gas): ', sco2_auxvar%xmol(gid,gid)
  write(86,*) '      liquid H [MJ/kmol]: ', sco2_auxvar%H(lid)
  write(86,*) '         gas H [MJ/kmol]: ', sco2_auxvar%H(gid)
  write(86,*) '      liquid U [MJ/kmol]: ', sco2_auxvar%U(lid)
  write(86,*) '         gas U [MJ/kmol]: ', sco2_auxvar%U(gid)
  write(86,*) '         liquid mobility: ', sco2_auxvar%mobility(lid)
  write(86,*) '            gas mobility: ', sco2_auxvar%mobility(gid)
  write(86,*) '      effective porosity: ', sco2_auxvar%effective_porosity
  write(86,*) '...'
  write(86,*) liquid_mass
  write(86,*) gas_mass
  write(86,*) liquid_mass*sco2_auxvar%U(lid) + &
              gas_mass*sco2_auxvar%U(gid)
  write(86,*) sco2_auxvar%pres(lid)
  write(86,*) sco2_auxvar%pres(gid)
  write(86,*) sco2_auxvar%pres(apid)
  write(86,*) sco2_auxvar%pres(cpid)
  write(86,*) sco2_auxvar%pres(vpid)
  write(86,*) sco2_auxvar%pres(spid)
  write(86,*) sco2_auxvar%temp
  write(86,*) sco2_auxvar%sat(lid)
  write(86,*) sco2_auxvar%sat(gid)
  write(86,*) sco2_auxvar%den(lid)
  write(86,*) sco2_auxvar%den_kg(lid)
  write(86,*) sco2_auxvar%den(gid)
  write(86,*) sco2_auxvar%den_kg(gid)
  write(86,*) sco2_auxvar%xmol(lid,lid)
  write(86,*) sco2_auxvar%xmol(gid,lid)
  write(86,*) sco2_auxvar%xmol(lid,gid)
  write(86,*) sco2_auxvar%xmol(gid,gid)
  write(86,*) sco2_auxvar%H(lid)
  write(86,*) sco2_auxvar%H(gid)
  write(86,*) sco2_auxvar%U(lid)
  write(86,*) sco2_auxvar%U(gid)
  write(86,*) ''
  write(86,*) sco2_auxvar%mobility(lid)
  write(86,*) sco2_auxvar%mobility(gid)
  write(86,*) sco2_auxvar%effective_porosity
  write(86,*) '--------------------------------------------------------'

  close(86)

end subroutine SCO2OutputAuxVars1

! ************************************************************************** !

subroutine SCO2OutputAuxVars2(sco2_auxvars,global_auxvars,option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Glenn Hammond
  ! Date: 02/18/13
  !

  use Global_Aux_module
  use Option_module

  implicit none

  type(sco2_auxvar_type) :: sco2_auxvars(0:,:)
  type(global_auxvar_type) :: global_auxvars(:)
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: apid, cpid, vpid
  PetscInt :: gid, lid, acid, wid
  PetscInt :: i, n, idof

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id

  string = 'sco2_auxvar.txt'
  open(unit=86,file=string)

  n = size(global_auxvars)

100 format(a,100('','',i9))

  write(86,'(a,100('','',i9))') '             cell id: ', &
    ((i,i=1,n),idof=0,3)
  write(86,'(a,100('','',i2))') '                idof: ', &
    ((idof,i=1,n),idof=0,3)
  write(86,'(a,100('','',i2))') '               state: ', &
    (global_auxvars(i)%istate,i=1,n)
  write(86,100) '      liquid pressure: ', &
    ((sco2_auxvars(idof,i)%pres(lid),i=1,n),idof=0,3)
  write(86,100) '         gas pressure: ', &
    ((sco2_auxvars(idof,i)%pres(gid),i=1,n),idof=0,3)
  write(86,100) '         air pressure: ', &
    ((sco2_auxvars(idof,i)%pres(apid),i=1,n),idof=0,3)
  write(86,100) '   capillary pressure: ', &
    ((sco2_auxvars(idof,i)%pres(cpid),i=1,n),idof=0,3)
  write(86,100) '       vapor pressure: ', &
    ((sco2_auxvars(idof,i)%pres(vpid),i=1,n),idof=0,3)
  write(86,100) '      temperature [C]: ', &
    ((sco2_auxvars(idof,i)%temp,i=1,n),idof=0,3)
  write(86,100) '    liquid saturation: ', &
    ((sco2_auxvars(idof,i)%sat(lid),i=1,n),idof=0,3)
  write(86,100) '       gas saturation: ', &
    ((sco2_auxvars(idof,i)%sat(gid),i=1,n),idof=0,3)
  write(86,100) 'liquid density [kmol]: ', &
    ((sco2_auxvars(idof,i)%den(lid),i=1,n),idof=0,3)
  write(86,100) '  liquid density [kg]: ', &
    ((sco2_auxvars(idof,i)%den_kg(lid),i=1,n),idof=0,3)
  write(86,100) '   gas density [kmol]: ', &
    ((sco2_auxvars(idof,i)%den(gid),i=1,n),idof=0,3)
  write(86,100) '     gas density [kg]: ', &
    ((sco2_auxvars(idof,i)%den_kg(gid),i=1,n),idof=0,3)
  write(86,100) '  X (water in liquid): ', &
    ((sco2_auxvars(idof,i)%xmol(lid,lid),i=1,n),idof=0,3)
  write(86,100) '    X (air in liquid): ', &
    ((sco2_auxvars(idof,i)%xmol(gid,lid),i=1,n),idof=0,3)
  write(86,100) '     X (water in gas): ', &
    ((sco2_auxvars(idof,i)%xmol(lid,gid),i=1,n),idof=0,3)
  write(86,100) '       X (air in gas): ', &
    ((sco2_auxvars(idof,i)%xmol(gid,gid),i=1,n),idof=0,3)
  write(86,100) '   liquid H [MJ/kmol]: ', &
    ((sco2_auxvars(idof,i)%H(lid),i=1,n),idof=0,3)
  write(86,100) '      gas H [MJ/kmol]: ', &
    ((sco2_auxvars(idof,i)%H(gid),i=1,n),idof=0,3)
  write(86,100) '   liquid U [MJ/kmol]: ', &
    ((sco2_auxvars(idof,i)%U(lid),i=1,n),idof=0,3)
  write(86,100) '      gas U [MJ/kmol]: ', &
    ((sco2_auxvars(idof,i)%U(gid),i=1,n),idof=0,3)
  write(86,*)
  write(86,100) '      liquid mobility: ', &
    ((sco2_auxvars(idof,i)%mobility(lid),i=1,n),idof=0,3)
  write(86,100) '         gas mobility: ', &
    ((sco2_auxvars(idof,i)%mobility(gid),i=1,n),idof=0,3)
  write(86,100) '   effective porosity: ', &
    ((sco2_auxvars(idof,i)%effective_porosity,i=1,n),idof=0,3)

  close(86)

end subroutine SCO2OutputAuxVars2

! ************************************************************************** !

subroutine SCO2AuxDestroy(aux)
  !
  ! Deallocates an SCO2 auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 12/05/23
  !

  use Utility_module, only : DeallocateArray
  
  implicit none

  type(sco2_type), pointer :: aux

  if (.not. associated(aux)) return

  call SCO2AuxVarDestroy(aux%auxvars)
  call SCO2AuxVarDestroy(aux%auxvars_bc)
  call SCO2AuxVarDestroy(aux%auxvars_ss)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%sco2_parameter)) then
    call DeallocateArray(aux%sco2_parameter%diffusion_coefficient)
    deallocate(aux%sco2_parameter)
  endif
  nullify(aux%sco2_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine SCO2AuxDestroy

! ************************************************************************** !
subroutine SCO2AuxVarSingleDestroy(auxvar)
  !
  ! Deallocates an SCO2 mode auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 12/05/23
  !

  implicit none

  type(sco2_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call SCO2AuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine SCO2AuxVarSingleDestroy

! ************************************************************************** !

subroutine SCO2AuxVarArray1Destroy(auxvars)
  !
  ! Deallocates an SCO2 mode auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 12/05/23
  !

  implicit none

  type(sco2_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call SCO2AuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine SCO2AuxVarArray1Destroy

! ************************************************************************** !

subroutine SCO2AuxVarArray2Destroy(auxvars)
  !
  ! Deallocates an SCO2 mode auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 12/05/23
  !

  implicit none

  type(sco2_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call SCO2AuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine SCO2AuxVarArray2Destroy

! ************************************************************************** !

subroutine SCO2AuxVarStrip(auxvar)
  !
  ! Deallocates an SCO2 auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 12/05/23
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(sco2_auxvar_type) :: auxvar

  call DeallocateArray(auxvar%pres)
  call DeallocateArray(auxvar%sat)
  call DeallocateArray(auxvar%den)
  call DeallocateArray(auxvar%den_kg)
  call DeallocateArray(auxvar%xmass)
  call DeallocateArray(auxvar%xmol)
  call DeallocateArray(auxvar%H)
  call DeallocateArray(auxvar%U)
  call DeallocateArray(auxvar%kr)
  call DeallocateArray(auxvar%mobility)
  call DeallocateArray(auxvar%effective_diffusion_coeff)
  call DeallocateArray(auxvar%visc)
  call DeallocateArray(auxvar%tortuosity)
  call DeallocateArray(auxvar%dispersivity)

end subroutine SCO2AuxVarStrip

! ************************************************************************** !
end module SCO2_Aux_module
