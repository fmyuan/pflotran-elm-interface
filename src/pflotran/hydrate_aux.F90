module Hydrate_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none

  private

  PetscBool, public :: hydrate_full_convergence = PETSC_FALSE
  PetscBool, public :: hydrate_use_governors = PETSC_FALSE
  PetscBool, public :: hydrate_check_updates = PETSC_FALSE
  PetscBool, public :: hydrate_truncate_updates = PETSC_TRUE
  PetscBool, public :: hydrate_use_henry_co2 = PETSC_FALSE
  PetscBool, public :: hydrate_print_state_transition = PETSC_TRUE
  PetscBool, public :: hydrate_analytical_derivatives = PETSC_FALSE
  PetscBool, public :: hydrate_immiscible = PETSC_FALSE
  PetscBool, public :: hydrate_central_diff_jacobian = PETSC_FALSE
  PetscBool, public :: hydrate_restrict_state_chng = PETSC_FALSE
  PetscReal, public :: window_epsilon = 0.d0 !1.d-4 !0.d0
  PetscReal, public :: hydrate_phase_chng_epsilon = 0.d0 !1.d-6
  PetscReal, public :: hydrate_max_pressure_change = 5.d4
  PetscInt, public :: hydrate_max_it_before_damping = UNINITIALIZED_INTEGER
  PetscReal, public :: hydrate_damping_factor = 0.6d0
  PetscInt, public :: hydrate_debug_cell_id = UNINITIALIZED_INTEGER
  PetscBool, public :: hydrate_diffuse_xmol = PETSC_TRUE
  PetscBool, public :: hydrate_temp_dep_gas_air_diff = PETSC_TRUE
  PetscInt, public :: hydrate_diffusion_model = ONE_INTEGER
  PetscBool, public :: hydrate_spycher_simple = PETSC_FALSE !PETSC_TRUE
  PetscBool, public :: hydrate_harmonic_diff_density = PETSC_TRUE
  PetscInt, public :: hydrate_newton_iteration_number = 0
  PetscInt, public :: hydrate_sub_newton_iter_num = 0
  PetscInt, public :: hydrate_newtontrdc_prev_iter_num = 0
  PetscBool, public :: hydrate_newtontrdc_hold_inner = PETSC_FALSE
  PetscBool, public :: hydrate_allow_state_change = PETSC_TRUE
  PetscBool, public :: hydrate_force_iteration = PETSC_FALSE
  PetscBool, public :: hydrate_state_changed = PETSC_FALSE
  PetscReal, public :: hydrate_bc_reference_pressure = 101325
  PetscReal, public :: hydrate_min_xmol = 1.d-10

  !Salinity
  PetscReal, public :: hydrate_xmass_nacl = 0.d0
  PetscReal, public :: hydrate_xmol_nacl = 0.d0 
  PetscInt, parameter, public :: HYDRATE_FORMER_NULL = ZERO_INTEGER
  PetscInt, parameter, public :: HYDRATE_FORMER_CH4 = ONE_INTEGER
  PetscInt, parameter, public :: HYDRATE_FORMER_CO2 = TWO_INTEGER

  !Boolean for the gas used
  PetscInt, public :: hydrate_former = HYDRATE_FORMER_CH4

  PetscBool, public :: hyd_chk_max_dpl_liq_state_only = PETSC_FALSE
  PetscBool, public :: hydrate_high_temp_ts_cut = PETSC_FALSE

  ! debugging
  PetscInt, public :: hydrate_ni_count
  PetscInt, public :: hydrate_ts_cut_count
  PetscInt, public :: hydrate_ts_count

  ! thermodynamic state of fluid ids
  !PetscInt, parameter, public :: NULL_STATE = 0
  !PetscInt, parameter, public :: L_STATE = 1
  !PetscInt, parameter, public :: G_STATE = 2
  !PetscInt, parameter, public :: GA_STATE = 3
  !PetscInt, parameter, public :: ANY_STATE = 4
  !PetscInt, parameter, public :: MULTI_STATE = 5

  PetscInt, parameter, public :: PREV_TS = 1
  PetscInt, parameter, public :: PREV_IT = 2

  PetscInt, parameter, public :: HYDRATE_LIQUID_PRESSURE_DOF = 1
  PetscInt, parameter, public :: HYDRATE_GAS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: HYDRATE_2PH_STATE_AIR_PRESSURE_DOF = 3
  PetscInt, parameter, public :: HYDRATE_G_STATE_AIR_PRESSURE_DOF = 2
  PetscInt, parameter, public :: HYDRATE_GAS_SATURATION_DOF = 2

  PetscInt, parameter, public :: HYDRATE_ENERGY_DOF = 3
  PetscInt, parameter, public :: HYDRATE_L_STATE_X_MOLE_DOF = 2

  PetscInt, parameter, public :: HYDRATE_STATE_INDEX = 1
  PetscInt, parameter, public :: HYDRATE_LIQUID_EQUATION_INDEX = 1
  PetscInt, parameter, public :: HYDRATE_GAS_EQUATION_INDEX = 2
  PetscInt, parameter, public :: HYDRATE_ENERGY_EQUATION_INDEX = 3

  PetscInt, parameter, public :: HYDRATE_LIQUID_PRESSURE_INDEX = 2
  PetscInt, parameter, public :: HYDRATE_GAS_PRESSURE_INDEX = 3
  PetscInt, parameter, public :: HYDRATE_AIR_PRESSURE_INDEX = 4
  PetscInt, parameter, public :: HYDRATE_LIQ_MOLE_FRACTION_INDEX = 5
  PetscInt, parameter, public :: HYDRATE_HYD_MOLE_FRACTION_INDEX = 6
  PetscInt, parameter, public :: HYDRATE_ICE_MOLE_FRACTION_INDEX = 7
  PetscInt, parameter, public :: HYDRATE_TEMPERATURE_INDEX = 8
  PetscInt, parameter, public :: HYDRATE_GAS_SATURATION_INDEX = 9
  PetscInt, parameter, public :: HYDRATE_LIQ_SATURATION_INDEX = 10
  PetscInt, parameter, public :: HYDRATE_HYD_SATURATION_INDEX = 11
  PetscInt, parameter, public :: HYDRATE_ICE_SATURATION_INDEX = 12
  PetscInt, parameter, public :: HYDRATE_LIQUID_FLUX_INDEX = 13
  PetscInt, parameter, public :: HYDRATE_GAS_FLUX_INDEX = 14
  PetscInt, parameter, public :: HYDRATE_ENERGY_FLUX_INDEX = 15
  PetscInt, parameter, public :: HYDRATE_LIQUID_CONDUCTANCE_INDEX = 16
  PetscInt, parameter, public :: HYDRATE_GAS_CONDUCTANCE_INDEX = 17
  PetscInt, parameter, public :: HYDRATE_ONE_INDEX = 18
  PetscInt, parameter, public :: HYDRATE_TWO_INDEX = 19
  PetscInt, parameter, public :: HYDRATE_THREE_INDEX = 20
  PetscInt, parameter, public :: HYDRATE_MAX_INDEX = 21

  PetscInt, parameter, public :: HYDRATE_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: HYDRATE_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: HYDRATE_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: HYDRATE_UPDATE_FOR_BOUNDARY = 2
  PetscInt, parameter, public :: HYDRATE_UPDATE_FOR_SS = 3

  PetscReal, parameter, public :: HYDRATE_IMMISCIBLE_VALUE = 1.d-10
  PetscReal, parameter, public :: HYDRATE_PRESSURE_SCALE = 1.d0
  PetscReal, parameter, public :: HYDRATE_REFERENCE_PRESSURE = 101325.d0

  ! these variables, which are global to hydrate, can be modified
  PetscInt, public :: dof_to_primary_variable(3,15)
  PetscInt, public :: hydrate_2ph_energy_dof = HYDRATE_TEMPERATURE_INDEX

  PetscInt, parameter, public :: HYD_MULTI_STATE = -2
  PetscInt, parameter, public :: HYD_ANY_STATE = -1
  PetscInt, parameter, public :: NULL_STATE = 0
  PetscInt, parameter, public :: L_STATE = 1
  PetscInt, parameter, public :: G_STATE = 2
  PetscInt, parameter, public :: H_STATE = 3 !5 (4 and 5 conflict with
  PetscInt, parameter, public :: I_STATE = 4 !4 ANY_STATE and MULTI_STATE)
  PetscInt, parameter, public :: GA_STATE = 5
  PetscInt, parameter, public :: HG_STATE = 6
  PetscInt, parameter, public :: HA_STATE = 7
  PetscInt, parameter, public :: HI_STATE = 8
  PetscInt, parameter, public :: GI_STATE = 9
  PetscInt, parameter, public :: AI_STATE = 10
  PetscInt, parameter, public :: HGA_STATE = 11
  PetscInt, parameter, public :: HAI_STATE = 12
  PetscInt, parameter, public :: HGI_STATE = 13
  PetscInt, parameter, public :: GAI_STATE = 14
  PetscInt, parameter, public :: HGAI_STATE = 15

  PetscInt, parameter :: lid = 1
  PetscInt, parameter :: gid = 2
  PetscInt, parameter :: hid = 3
  PetscInt, parameter :: iid = 4

  !Structure 1 methane hydrate:
  PetscReal, parameter :: CH4_HYDRATION_NUMBER = 5.75d0
  PetscReal, parameter :: CO2_HYDRATION_NUMBER = 5.75d0
  PetscReal, parameter :: CH4_HYDRATE_DENSITY_KG = 900.d0 !kg/m^3
  PetscReal, parameter :: CH4_HYDRATE_DENSITY = 52.15551276d0 !kmol/m^3
  PetscReal, parameter :: CO2_HYDRATE_DENSITY = 63.7456267067 !7.45d0 !kmol/m^3
  PetscReal, parameter :: CO2_HYDRATE_DENSITY_KG = 1100.d0 !kg/m^3
  PetscReal, parameter :: MW_CH4 = 16.04d0
  PetscReal, parameter :: MW_H2O = 18.01d0
  PetscReal, parameter :: MW_NACL = 58.44277d0
  PetscReal, parameter :: L_GH = 7161.d0 ! enthalpy of fusion, J/mol
  PetscReal, public :: hydrate_fmw_comp(3) = [MW_H2O,MW_CH4,MW_NACL]
  PetscReal, parameter, public :: MOL_RATIO_METH = 0.14285714285d0
  PetscReal, parameter :: MOL_RATIO_H2O = 1.d0 - MOL_RATIO_METH

  PetscReal, public :: TQD = 0.d0 !1.d-2 !Quad point temperature (C)

  !Ice:
  PetscReal, parameter :: ICE_DENSITY_KG = 920.d0 !kg/m^3
  PetscReal, parameter :: ICE_DENSITY = 50.86d0 !mol/L
  PetscReal, parameter :: L_ICE = 6033.54 !J/mol

  PetscReal, parameter :: CO2_REFERENCE_SURFACE_TENSION = 0.072d0 ! N/m
  PetscReal, parameter :: SALT_REFERENCE_TEMPERATURE = 293.15d0
  PetscReal, parameter :: LIQUID_REFERENCE_VISCOSITY = 1.01764892595942d-3
  PetscReal, parameter, public :: LIQUID_REFERENCE_DENSITY = 998.32142721500441
  PetscReal, parameter, public :: HYD_REFERENCE_PRESSURE = 101325.d0


  PetscReal, parameter :: lambda_hyd = 0.49d0 !W/m-K

  PetscInt, public :: permeability_reduction_model = TWO_INTEGER
  PetscInt, public :: hydrate_perm_scaling_function = ONE_INTEGER
  PetscInt, public :: hydrate_phase_boundary = ONE_INTEGER
  PetscInt, public :: hydrate_henrys_constant = ONE_INTEGER
  PetscInt, public :: hydrate_tcond = TWO_INTEGER
  PetscBool, public :: hydrate_perm_scaling = PETSC_TRUE
  PetscBool, public :: hydrate_eff_sat_scaling = PETSC_TRUE
  PetscBool, public :: hydrate_no_ice_density_change = PETSC_FALSE
  PetscBool, public :: hydrate_with_gibbs_thomson = PETSC_FALSE
  PetscBool, public :: hydrate_gt_3phase = PETSC_FALSE
  PetscBool, public :: hydrate_adjust_ghsz_solubility = PETSC_FALSE
  PetscBool, public :: hydrate_with_sedimentation = PETSC_FALSE
  PetscBool, public :: hydrate_no_pc = PETSC_FALSE
  PetscBool, public :: hydrate_with_methanogenesis = PETSC_FALSE
  PetscBool, public :: hydrate_compute_surface_tension = PETSC_FALSE

  type, public :: hydrate_auxvar_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal :: sl_min  ! min liquid saturation for hysteresis
    PetscReal, pointer :: den(:)    ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: den_kg(:) ! (iphase) kg/m^3 phase
    PetscReal :: temp
    PetscReal :: m_salt(2)      ! (kg NaCl / kg (brine + precipitate)), kg NaCl
    PetscReal, pointer :: xmass(:,:) ! (icomp,iphase)
    PetscReal, pointer :: xmol(:,:) ! (icomp,iphase)
    PetscReal, pointer :: effective_diffusion_coeff(:,:) ! (icomp,iphase)
    PetscReal, pointer :: dispersivity(:,:) ! (icomp,iphase)
    PetscReal, pointer :: H(:) ! MJ/kmol
    PetscReal, pointer :: U(:) ! MJ/kmol
    PetscReal, pointer :: kr(:)
    PetscReal, pointer :: visc(:) ! Pa-s
    PetscReal, pointer :: mobility(:) ! relative perm / kinematic viscosity
    PetscReal :: effective_porosity ! factors in compressibility
    PetscReal :: effective_permeability
    PetscReal, pointer :: tortuosity(:) ! (iphase)
    PetscReal :: perm_base
    PetscReal :: v_sed
    PetscReal :: srl
    PetscReal :: srg
    PetscReal :: pert
    PetscBool :: istatechng
  end type hydrate_auxvar_type

  type, public :: hydrate_parameter_type
    PetscReal, pointer :: diffusion_coefficient(:,:) ! (icomp,iphase)
    PetscReal :: newton_inf_scaled_res_tol
    PetscBool :: check_post_converged
    type(methanogenesis_type), pointer :: methanogenesis
  end type hydrate_parameter_type

  type, public :: methanogenesis_type
    character(len=MAXWORDLENGTH) :: source_name
    PetscReal :: alpha
    PetscReal :: k_alpha
    PetscReal :: lambda
    PetscReal :: omega
    PetscReal :: z_smt
    type(methanogenesis_type), pointer :: next
  end type methanogenesis_type

  type, public :: hydrate_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(hydrate_parameter_type), pointer :: hydrate_parameter
    type(hydrate_auxvar_type), pointer :: auxvars(:,:)
    type(hydrate_auxvar_type), pointer :: auxvars_bc(:)
    type(hydrate_auxvar_type), pointer :: auxvars_ss(:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type hydrate_type

  interface HydrateAuxVarDestroy
    module procedure HydrateAuxVarSingleDestroy
    module procedure HydrateAuxVarArray1Destroy
    module procedure HydrateAuxVarArray2Destroy
  end interface HydrateAuxVarDestroy

  interface HydrateOutputAuxVars
    module procedure HydrateOutputAuxVars1
    module procedure HydrateOutputAuxVars2
  end interface HydrateOutputAuxVars

  public :: HydrateAuxCreate, &
            HydrateMethanogenesisCreate, &
            HydrateAuxDestroy, &
            HydrateAuxSetEnergyDOF, &
            HydrateAuxVarCompute, &
            HydrateAuxVarInit, &
            HydrateAuxVarCopy, &
            HydrateAuxVarDestroy, &
            HydrateAuxVarStrip, &
            HydrateAuxVarUpdateState, &
            HydrateAuxVarPerturb, &
            HydratePrintAuxVars, &
            HydrateOutputAuxVars, &
            HydrateCompositeThermalCond,&
            HydratePE, &
            HydrateMethanogenesis, &
            HydrateGHSZSolubilityCorrection, &
            CalcFreezingTempDepression, &
            EOSHydrateEnthalpy, &
            HydrateComputeSaltSolubility, &
            HydrateBrineSaturationPressure, &
            HydrateVaporPressureBrine, &
            HydrateBrineDensity, &
            HydrateHenryCO2, &
            HydrateComputeSaltDensity, &
            HydrateEquilibrate


contains

! ************************************************************************** !

function HydrateAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Option_module

  implicit none

  type(option_type) :: option

  type(hydrate_type), pointer :: HydrateAuxCreate

  type(hydrate_type), pointer :: aux

  ! L_STATE,G_STATE,H_STATE,I_STATE,GA_STATE,HG_STATE,HA_STATE,HI_STATE,
  ! GI_STATE,AI_STATE,HGA_STATE,HAI_STATE,HGI_STATE,GAI_STATE,HGAI_STATE
  dof_to_primary_variable(1:3,1:15) = &
             !L_STATE
    reshape([HYDRATE_LIQUID_PRESSURE_INDEX, HYDRATE_LIQ_MOLE_FRACTION_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !G_STATE
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_AIR_PRESSURE_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !H_STATE
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_HYD_MOLE_FRACTION_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !I_STATE
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_ICE_MOLE_FRACTION_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !GA_STATE
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_GAS_SATURATION_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !HG_STATE
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_GAS_SATURATION_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !HA_STATE
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_HYD_SATURATION_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !HI_STATE
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_HYD_SATURATION_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !GI_STATE 2INDEX = Si
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_TWO_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !AI_STATE 3INDEX = Sl
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_LIQ_MOLE_FRACTION_INDEX, &
             HYDRATE_THREE_INDEX, &
             !HGA_STATE
             HYDRATE_LIQ_SATURATION_INDEX, HYDRATE_HYD_SATURATION_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !HAI_STATE 2INDEX = Sl
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_TWO_INDEX, &
             HYDRATE_ICE_SATURATION_INDEX, &
             !HGI_STATE 1INDEX = Si
             HYDRATE_ONE_INDEX, HYDRATE_HYD_SATURATION_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !GAI_STATE 2INDEX = Sg
             HYDRATE_GAS_PRESSURE_INDEX, HYDRATE_TWO_INDEX, &
             HYDRATE_TEMPERATURE_INDEX, &
             !HGAI_STATE
             HYDRATE_LIQ_SATURATION_INDEX, HYDRATE_GAS_SATURATION_INDEX, &
             HYDRATE_ICE_SATURATION_INDEX &
             ],shape(dof_to_primary_variable))

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

  allocate(aux%hydrate_parameter)
  allocate(aux%hydrate_parameter%diffusion_coefficient(option%nflowspec, &
                                                       option%nphase))
  aux%hydrate_parameter%diffusion_coefficient(:,LIQUID_PHASE) = &
                                                           UNINITIALIZED_DOUBLE
  aux%hydrate_parameter%diffusion_coefficient(:,GAS_PHASE) = 2.13d-5
  aux%hydrate_parameter%newton_inf_scaled_res_tol = 1.d-50
  aux%hydrate_parameter%check_post_converged = PETSC_FALSE

  HydrateAuxCreate => aux

end function HydrateAuxCreate

! ************************************************************************** !

function HydrateMethanogenesisCreate()

  !
  ! Allocate and initialize methanogenesis object
  !
  ! Author: Michael Nole
  ! Date: 11/21/19
  !

  type(methanogenesis_type), pointer :: HydrateMethanogenesisCreate
  type(methanogenesis_type), pointer :: methanogenesis

  allocate(methanogenesis)

  methanogenesis%source_name = ''
  methanogenesis%alpha = UNINITIALIZED_DOUBLE
  methanogenesis%k_alpha = UNINITIALIZED_DOUBLE
  methanogenesis%lambda = UNINITIALIZED_DOUBLE
  methanogenesis%omega = UNINITIALIZED_DOUBLE
  methanogenesis%z_smt = UNINITIALIZED_DOUBLE

  HydrateMethanogenesisCreate => methanogenesis

end function HydrateMethanogenesisCreate

! ************************************************************************** !
subroutine HydrateAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%istate_store = NULL_STATE
  auxvar%temp = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%effective_permeability = -999.d0
  auxvar%perm_base = -999.9d0
  auxvar%v_sed = 0.d0
  auxvar%srl = 0.d0
  auxvar%srg = 0.d0
  auxvar%pert = 0.d0
  auxvar%m_salt = 0.d0
  auxvar%sl_min = 1.d0
  auxvar%istatechng = PETSC_FALSE

  allocate(auxvar%pres(option%nphase+FOUR_INTEGER))
  auxvar%pres = 0.d0
  allocate(auxvar%sat(option%nphase+1))
  auxvar%sat = 0.d0
  allocate(auxvar%den(3 + option%nphase)) ! Pure component and mixture
  auxvar%den = 0.d0
  allocate(auxvar%den_kg(3 + option%nphase)) ! Pure component and mixture
  auxvar%den_kg = 0.d0
  allocate(auxvar%xmass(option%nflowspec,option%nphase))
  auxvar%xmass = 0.d0
  allocate(auxvar%xmol(option%nflowspec,option%nphase))
  auxvar%xmol = 0.d0
  allocate(auxvar%effective_diffusion_coeff(option%nflowspec,option%nphase))
  auxvar%effective_diffusion_coeff = 0.d0
  allocate(auxvar%dispersivity(option%nflowspec,option%nphase))
  auxvar%dispersivity = 0.d0
  allocate(auxvar%H(3 + option%nphase)) ! Pure component and mixture
  auxvar%H = 0.d0
  allocate(auxvar%U(3 + option%nphase)) ! Pure component and mixture
  auxvar%U = 0.d0
  allocate(auxvar%mobility(option%nphase))
  auxvar%mobility = 0.d0
  allocate(auxvar%kr(option%nphase))
  auxvar%kr = 0.d0
  allocate(auxvar%visc(option%nphase))
  auxvar%visc = 0.d0
  allocate(auxvar%tortuosity(option%nphase))
  auxvar%tortuosity = 1.d0
  

end subroutine HydrateAuxVarInit

! ************************************************************************** !

subroutine HydrateAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%istate_store = auxvar%istate_store
  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%xmass = auxvar%xmass
  auxvar2%xmol = auxvar%xmol
  auxvar2%effective_diffusion_coeff = auxvar%effective_diffusion_coeff
  auxvar2%dispersivity = auxvar%dispersivity
  auxvar2%H = auxvar%H
  auxvar2%U = auxvar%U
  auxvar2%mobility = auxvar%mobility
  auxvar2%kr = auxvar%kr
  auxvar2%visc = auxvar%visc
  auxvar2%tortuosity = auxvar%tortuosity
  auxvar2%perm_base = auxvar%perm_base
  auxvar2%v_sed = auxvar%v_sed
  auxvar2%srl = auxvar%srl
  auxvar2%srg = auxvar%srg
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%effective_permeability = auxvar%effective_permeability
  auxvar2%pert = auxvar%pert

end subroutine HydrateAuxVarCopy

! ************************************************************************** !

subroutine HydrateAuxSetEnergyDOF(energy_keyword,option)
  !
  ! Sets the two phase primary dependent variable for energy based on user
  ! input.
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !
  use Option_module
  use String_module

  implicit none

  character(len=MAXWORDLENGTH) :: energy_keyword
  type(option_type) :: option

  call StringToUpper(energy_keyword)
  select case(energy_keyword)
    case('TEMPERATURE')
      hydrate_2ph_energy_dof = HYDRATE_TEMPERATURE_INDEX
    case('AIR_PRESSURE')
      hydrate_2ph_energy_dof = HYDRATE_AIR_PRESSURE_INDEX
    case default
      option%io_buffer = 'Energy Keyword: ' // trim(energy_keyword) // &
                          ' not recognized in Hydrate Mode'
      call PrintErrMsg(option)
  end select

end subroutine HydrateAuxSetEnergyDOF

! ************************************************************************** !
subroutine HydrateAuxVarCompute(x,hyd_auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,hydrate_parameter, &
                                natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell, with gas hydrate physics
  ! Author: Michael Nole
  ! Date: 04/04/19
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  
  PetscReal :: x(option%nflowdof)
  type(hydrate_auxvar_type) :: hyd_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  type(hydrate_parameter_type), pointer :: hydrate_parameter
  PetscInt :: natural_id

  ! Phase ID's
  PetscInt :: lid, gid, hid, pid, pwid, pgid, pbid
  ! Component ID's
  PetscInt :: wid, acid, sid
  ! Other ID's
  PetscInt :: cpid, vpid, rvpid, spid, tgid, apid
  PetscReal :: xag, xwg, xal, xsl, xwl, xmolag, xmolwg, xmolal, &
               xmolsl, xmolwl
  PetscReal :: mw_mix
  PetscReal :: cell_pressure, water_vapor_pressure
  PetscReal :: den_water_vapor, den_kg_water_vapor
  PetscReal :: u_water_vapor, h_water_vapor
  PetscReal :: den_air, h_air, u_air
  PetscReal :: den_mol, den_steam_kg
  PetscReal :: den_steam
  PetscReal :: salt_solubility, x_salt_dissolved
  PetscReal :: drho_dp, drho_dT
  PetscReal :: visc_water, visc_brine, visc_a
  PetscReal :: sl_temp, pva
  PetscReal :: xmol_air_in_gas, xmol_water_in_gas
  PetscReal :: dkrl_dsatl
  PetscReal :: dkrg_dsatl
  PetscReal :: beta_gl
  PetscReal :: K_H_tilde
  PetscReal :: Hg_mixture_fractioned
  PetscReal :: H_hyd, U_ice, PE_hyd, du_ice_dT, du_ice_dP
  PetscReal :: aux(1)
  PetscReal :: hw
  PetscReal :: dpor_dp
  PetscReal :: dTf, Pc, dTfs, h_sat_eff, i_sat_eff, l_sat_eff, g_sat_eff
  PetscReal :: Tf_ice, T_temp
  PetscReal :: sigma, dP
  PetscReal :: sat_temp
  PetscErrorCode :: ierr
  PetscReal, parameter :: epsilon = 1.d-14

  ierr = 0

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  pid = option%precipitate_phase
  pwid = option%pure_water_phase
  pbid = option%pure_brine_phase
  tgid = option%trapped_gas_phase
  pgid = tgid ! pure gas

  wid = option%water_id
  acid = option%air_id
  sid = option%salt_id

  cpid = option%capillary_pressure_id
  apid = option%air_pressure_id
  vpid = option%vapor_pressure_id
  rvpid = option%reduced_vapor_pressure_id
  spid = option%saturation_pressure_id


  hyd_auxvar%H = 0.d0
  hyd_auxvar%U = 0.d0
  hyd_auxvar%pres = 0.d0
  hyd_auxvar%sat = 0.d0
  hyd_auxvar%den = 0.d0
  hyd_auxvar%den_kg = 0.d0
  hyd_auxvar%xmass = 0.d0
  hyd_auxvar%xmol = 0.d0
  hyd_auxvar%effective_diffusion_coeff = 0.d0
  hyd_auxvar%dispersivity = 0.d0
  hyd_auxvar%effective_porosity = material_auxvar%porosity_base
  hyd_auxvar%effective_permeability = 1.d0
  hyd_auxvar%tortuosity = 0.d0
  hyd_auxvar%mobility = 0.d0
  hyd_auxvar%kr = 0.d0

#if 0
  if (option%iflag >= HYDRATE_UPDATE_FOR_ACCUM) then
    if (option%iflag == HYDRATE_UPDATE_FOR_ACCUM) then
      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
        natural_id, x(1:3)
    else
    endif
  endif
#endif

  if (hydrate_former == HYDRATE_FORMER_CO2) then
    hyd_auxvar%den(hid) = CO2_HYDRATE_DENSITY
    hyd_auxvar%den_kg(hid) = CO2_HYDRATE_DENSITY_KG
    hyd_auxvar%xmol(acid,hid) = 1.d0 / (1.d0 + CO2_HYDRATION_NUMBER)
    hyd_auxvar%xmol(wid,hid) = 1.d0 - hyd_auxvar%xmol(acid,hid)
  elseif (hydrate_former == HYDRATE_FORMER_CH4) then
    hyd_auxvar%den(hid) = CH4_HYDRATE_DENSITY
    hyd_auxvar%den_kg(hid) = CH4_HYDRATE_DENSITY_KG
    hyd_auxvar%xmol(acid,hid) = 1.d0 / (1.d0 + CH4_HYDRATION_NUMBER)
    hyd_auxvar%xmol(wid,hid) = 1.d0 - hyd_auxvar%xmol(acid,hid)
  endif

  if (option%flow%density_depends_on_salinity) then
    hydrate_xmass_nacl = global_auxvar%m_nacl(1)
  endif

  if (hydrate_xmass_nacl > 0.d0) then
    call IceSalinityOffset(hydrate_xmass_nacl,dTfs)
  else
    dTfs = 0.d0
  endif

  Tf_ice = dTfs !Bulk freezing temperature

  select case(global_auxvar%istate)
    case(L_STATE)
!     ********* Aqueous State (A) ********************************
!     Primary variables: Pl, Xma, T
!
      hyd_auxvar%pres(lid) = x(HYDRATE_LIQUID_PRESSURE_DOF)
      hyd_auxvar%xmass(acid,lid) = x(HYDRATE_L_STATE_X_MOLE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      T_temp = hyd_auxvar%temp - Tf_ice

      if (T_temp <= 0.d0) then
        ! Clausius-Clayperon equation
        Pc = -(T_temp) * (L_ICE * ICE_DENSITY * 1.d6) / (Tf_ice + T273K)
        ! Get the corresponding liquid saturation
        call HydrateComputeSatHysteresis(characteristic_curves, &
                                    Pc, &
                                    hyd_auxvar%sl_min, &
                                    1.d0, hyd_auxvar%den_kg(lid), &
                                    hyd_auxvar%sat(lid), &
                                    hyd_auxvar%sat(tgid), &
                                    option)
        hyd_auxvar%sat(iid) = 1.d0 - hyd_auxvar%sat(lid)
      else
        hyd_auxvar%sat(lid) = 1.d0
      endif
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%pres(cpid) = 0.d0

      ! MAN: Replace with dissolved salt concentration:
      hyd_auxvar%m_salt(1) = hydrate_xmass_nacl
      ! kg NaCl/kg liquid
      call HydrateComputeSaltSolubility(hyd_auxvar%temp, salt_solubility)
      ! Dissolved salt mass fraction
      x_salt_dissolved = min(hyd_auxvar%m_salt(1),salt_solubility)
      call HydrateBrineSaturationPressure(hyd_auxvar%temp, &
                                         x_salt_dissolved, &
                                         hyd_auxvar%pres(spid))
      ! Brine density
      call HydrateBrineDensity(hyd_auxvar%temp, hyd_auxvar%pres(lid), &
                            x_salt_dissolved, hyd_auxvar%den_kg(pbid), option)
      ! Brine vapor pressure
      call HydrateVaporPressureBrine(hyd_auxvar%temp, hyd_auxvar%pres(spid), &
                                   hyd_auxvar%pres(cpid), &
                                   hyd_auxvar%den_kg(pbid), &
                                   x_salt_dissolved, hyd_auxvar%pres(rvpid))

      ! Pure water density
      call HydrateWaterDensity(hyd_auxvar%temp,hyd_auxvar%pres(rvpid), &
                            TWO_INTEGER,hyd_auxvar%den_kg(pwid), &
                            den_steam_kg,option)
      pva = max(hyd_auxvar%pres(lid) - hyd_auxvar%pres(rvpid), 0.d0)
      xsl = x_salt_dissolved
      call HydrateEquilibrate(hyd_auxvar%temp,hyd_auxvar%pres(lid), &
                           global_auxvar%istate, hyd_auxvar%sat(hid), &
                           hyd_auxvar%pres(apid), &
                           hyd_auxvar%pres(vpid), &
                           hyd_auxvar%pres(spid), &
                           hyd_auxvar%pres(rvpid), &
                           xag, xwg, xal, xsl, xwl, &
                           xmolag, xmolwg, xmolal, xmolsl, xmolwl, &
                           characteristic_curves, material_auxvar, option)

      hyd_auxvar%xmass(sid,lid) = x_salt_dissolved + &
                              (xsl-x_salt_dissolved) * &
                              (hyd_auxvar%xmass(acid,lid) / xal)

      hyd_auxvar%xmass(wid,lid) = 1.d0 - hyd_auxvar%xmass(sid,lid) - &
                                   hyd_auxvar%xmass(acid,lid)
      hyd_auxvar%xmass(wid,lid) = max(hyd_auxvar%xmass(wid,lid),0.d0)

      ! Populate all pressures, even though gas phase is not present.
      hyd_auxvar%pres(gid) = hyd_auxvar%pres(lid)
      hyd_auxvar%xmass(acid,gid) = 1.d0

      ! Update the liquid mole fractions
      mw_mix = 1.d0 / (hyd_auxvar%xmass(wid,lid)/hydrate_fmw_comp(1) + &
               hyd_auxvar%xmass(acid,lid)/hydrate_fmw_comp(2) + &
               hyd_auxvar%xmass(sid,lid)/hydrate_fmw_comp(3))
      hyd_auxvar%xmol(wid,lid) = hyd_auxvar%xmass(wid,lid)* &
                                  mw_mix/hydrate_fmw_comp(1)
      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%xmass(acid,lid)* &
                                     mw_mix/hydrate_fmw_comp(2)
      hyd_auxvar%xmol(sid,lid) = hyd_auxvar%xmass(sid,lid)* &
                                  mw_mix/hydrate_fmw_comp(3)

    case (G_STATE)
!     ********* Gas State (G) ********************************
!     Primary variables: Pg, Pa, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)

      hyd_auxvar%pres(apid) = x(HYDRATE_G_STATE_AIR_PRESSURE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      T_temp = hyd_auxvar%temp - Tf_ice
      
      ! hyd_auxvar%m_salt(2) = x(HYD_SALT_MASS_FRAC_DOF)

      ! Secondary Variables
      ! kg NaCl/kg liquid
      call HydrateComputeSaltSolubility(hyd_auxvar%temp, salt_solubility)
      if (hyd_auxvar%m_salt(2) > epsilon) then
        x_salt_dissolved = salt_solubility
      else
        x_salt_dissolved = 0.d0
      endif
      hyd_auxvar%xmass(sid,lid) = x_salt_dissolved
      call HydrateComputeSurfaceTension(hyd_auxvar%temp, &
                                     x_salt_dissolved, sigma)
      beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma
      call HydrateComputePcHysteresis(characteristic_curves, &
                                   hyd_auxvar%sat(lid), &
                                   hyd_auxvar%sat(tgid), &
                                   beta_gl,hyd_auxvar%pres(cpid), option)
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid) / beta_gl


      cell_pressure = max(hyd_auxvar%pres(gid),hyd_auxvar%pres(spid))
      ! cell_pressure = min(cell_pressure, 1.d8)

      hyd_auxvar%pres(rvpid) = max(hyd_auxvar%pres(gid) - &
                               hyd_auxvar%pres(apid), 0.d0)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(rvpid)
      pva = max(hyd_auxvar%pres(gid) - hyd_auxvar%pres(rvpid),0.d0)

      call HydrateWaterDensity(hyd_auxvar%temp,hyd_auxvar%pres(rvpid), &
                            TWO_INTEGER,hyd_auxvar%den_kg(pwid), &
                            den_steam_kg,option)
      ! Compute equilibrium mass and mole fractions for all components
      xsl = x_salt_dissolved
      call HydrateBrineSaturationPressure(hyd_auxvar%temp, &
                                         x_salt_dissolved, &
                                         hyd_auxvar%pres(spid))
      call HydrateEquilibrate(hyd_auxvar%temp,cell_pressure, &
                           global_auxvar%istate, hyd_auxvar%sat(hid), &
                           pva, hyd_auxvar%pres(vpid), &
                           hyd_auxvar%pres(spid), &
                           hyd_auxvar%pres(rvpid), &
                           xag, xwg, xal, xsl, xwl, &
                           xmolag, xmolwg, xmolal, xmolsl, xmolwl, &
                           characteristic_curves, material_auxvar, option)

      call HydrateComputeSaltDensity(hyd_auxvar%temp, cell_pressure, &
                                  hyd_auxvar%den_kg(pid))

      hyd_auxvar%sat(pid) = hyd_auxvar%m_salt(2) / &
                             (hyd_auxvar%den_kg(pid) * &
                             material_auxvar%volume)
      hyd_auxvar%sat(pid) = max(min(hyd_auxvar%sat(pid),1.d0),0.d0)
      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 1.d0

      ! Update mass fractions
      hyd_auxvar%xmass(acid,lid) = xal
      hyd_auxvar%xmass(wid,lid) = xwl
      hyd_auxvar%xmass(sid,lid) = xsl
      hyd_auxvar%xmass(acid,gid) = xag
      hyd_auxvar%xmass(wid,gid) = xwg

      ! Update mole fractions
      hyd_auxvar%xmol(acid,lid) = xmolal
      hyd_auxvar%xmol(wid,lid) = xmolwl
      hyd_auxvar%xmol(sid,lid) = xmolsl
      hyd_auxvar%xmol(acid,gid) = xmolag
      hyd_auxvar%xmol(wid,gid) = xmolwg

      ! Brine density
      call HydrateBrineDensity(hyd_auxvar%temp, cell_pressure, &
                            x_salt_dissolved, hyd_auxvar%den_kg(pbid), option)
      call HydrateDensityCompositeLiquid(hyd_auxvar%temp,hyd_auxvar%den_kg(pbid), &
                                  hyd_auxvar%xmass(acid,lid), &
                                  hyd_auxvar%den_kg(lid))
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      call HydrateBrineSaturationPressure(hyd_auxvar%temp, &
                                         x_salt_dissolved, &
                                         hyd_auxvar%pres(spid))

    case (H_STATE)
!     ********* Hydrate State (H) ********************************
!     Primary variables: Pg, Xmh, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      !x(HYDRATE_GAS_SATURATION_DOF) = MOL_RATIO_METH
      hyd_auxvar%xmol(acid,hid) = MOL_RATIO_METH
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      T_temp = hyd_auxvar%temp - Tf_ice

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 1.d0
      hyd_auxvar%sat(iid) = 0.d0

      call HydratePE(T_temp,hyd_auxvar%sat(hid),PE_hyd,dP, &
              characteristic_curves, material_auxvar, option)
      hyd_auxvar%pres(apid) = PE_hyd
      call EOSGasHenry(T_temp,hyd_auxvar%pres(spid),K_H_tilde, &
                       ierr)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

    case(I_STATE)
!     ********* Ice State (I) ********************************
!     Primary variables: Pg, Xmi, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      !x(HYDRATE_GAS_SATURATION_DOF) = 0.d0
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      T_temp = hyd_auxvar%temp - Tf_ice

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 1.d0

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid)
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      call EOSWaterSaturationPressure(T_temp, &
                                          hyd_auxvar%pres(spid),ierr)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)

      hyd_auxvar%xmol(acid,lid) = 0.d0
      hyd_auxvar%xmol(wid,lid) = 1.d0
      hyd_auxvar%xmol(wid,gid) = 0.d0
      hyd_auxvar%xmol(acid,gid) = 1.d0

    case(GA_STATE)
!     ********* Gas & Aqueous State (GA) ********************************
!     Primary variables: Pg, Sg, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      T_temp = hyd_auxvar%temp - Tf_ice

      if (T_temp <= 0.d0) then
        ! Clausius-Clayperon equation
        Pc = -(T_temp) * (L_ICE * ICE_DENSITY * 1.d6) / (Tf_ice + 273.15d0)
        ! Get the corresponding liquid saturation
        call HydrateComputeSatHysteresis(characteristic_curves, &
                                    Pc, &
                                    hyd_auxvar%sl_min, &
                                    1.d0, hyd_auxvar%den_kg(lid), &
                                    hyd_auxvar%sat(lid), &
                                    hyd_auxvar%sat(tgid), &
                                    option)
        hyd_auxvar%sat(iid) = max(0.d0, 1.d0 - hyd_auxvar%sat(lid) - &
                              hyd_auxvar%sat(gid))
      else
        hyd_auxvar%sat(lid) = 1.d0 - hyd_auxvar%sat(gid)
        hyd_auxvar%sat(iid) = 0.d0
      endif
      hyd_auxvar%sat(hid) = 0.d0

      ! Secondary Variables

      call HydrateComputeSaltSolubility(hyd_auxvar%temp, salt_solubility)
      ! Dissolved salt mass fraction
      x_salt_dissolved = min(hyd_auxvar%m_salt(1),salt_solubility)

      call HydrateComputeSurfaceTension(hyd_auxvar%temp, &
                                     x_salt_dissolved, sigma)
      beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma
      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        sl_temp = 1.d0 - hyd_auxvar%sat(gid)
        call HydrateComputePcHysteresis(characteristic_curves, &
                                        sl_temp, &
                                        hyd_auxvar%sat(tgid), &
                                        beta_gl, hyd_auxvar%pres(cpid), &
                                        option)
      endif

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      call HydrateBrineSaturationPressure(hyd_auxvar%temp, &
                                         x_salt_dissolved, &
                                         hyd_auxvar%pres(spid))
      cell_pressure = max(hyd_auxvar%pres(gid),hyd_auxvar%pres(spid))
      ! Brine density
      call HydrateBrineDensity(hyd_auxvar%temp, cell_pressure, &
                           x_salt_dissolved, hyd_auxvar%den_kg(pbid), option)
      ! Brine vapor pressure
      call HydrateVaporPressureBrine(hyd_auxvar%temp, hyd_auxvar%pres(spid), &
                                   hyd_auxvar%pres(cpid), &
                                   hyd_auxvar%den_kg(pbid), &
                                   x_salt_dissolved, hyd_auxvar%pres(rvpid))

      ! Pure water density
      xsl = x_salt_dissolved
      call HydrateWaterDensity(hyd_auxvar%temp,hyd_auxvar%pres(rvpid), &
                            TWO_INTEGER,hyd_auxvar%den_kg(pwid), &
                            den_steam_kg,option)
      pva = max(hyd_auxvar%pres(gid) - hyd_auxvar%pres(rvpid), 0.d0)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(rvpid)
      call HydrateEquilibrate(hyd_auxvar%temp,cell_pressure, &
                           global_auxvar%istate, hyd_auxvar%sat(hid), &
                           hyd_auxvar%pres(apid), &
                           hyd_auxvar%pres(vpid), &
                           hyd_auxvar%pres(spid), &
                           hyd_auxvar%pres(rvpid), &
                           xag, xwg, xal, xsl, xwl, &
                           xmolag, xmolwg, xmolal, xmolsl, xmolwl, &
                           characteristic_curves, material_auxvar, option)

      ! Update mass fractions
      hyd_auxvar%xmass(acid,lid) = xal
      hyd_auxvar%xmass(wid,lid) = xwl
      hyd_auxvar%xmass(sid,lid) = xsl
      hyd_auxvar%xmass(acid,gid) = xag
      hyd_auxvar%xmass(wid,gid) = xwg

      ! Update mole fractions
      hyd_auxvar%xmol(acid,lid) = xmolal
      hyd_auxvar%xmol(wid,lid) = xmolwl
      hyd_auxvar%xmol(sid,lid) = xmolsl
      hyd_auxvar%xmol(acid,gid) = xmolag
      hyd_auxvar%xmol(wid,gid) = xmolwg

    case(HG_STATE)
!     ********* Hydrate & Gas State (HG) ********************************
!     Primary variables: Pg, Sg, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      T_temp = hyd_auxvar%temp - Tf_ice

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(hid) = 1.d0 - hyd_auxvar%sat(gid)
      hyd_auxvar%sat(iid) = 0.d0

      call HydrateComputeSaltSolubility(hyd_auxvar%temp, salt_solubility)
      ! Dissolved salt mass fraction
      x_salt_dissolved = min(hyd_auxvar%m_salt(1),salt_solubility)

      call HydrateComputeSurfaceTension(hyd_auxvar%temp, &
                                     x_salt_dissolved, sigma)
      beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma
      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        sl_temp = 1.d0 - hyd_auxvar%sat(gid)
        call HydrateComputePcHysteresis(characteristic_curves, &
                                        sl_temp, &
                                        hyd_auxvar%sat(tgid), &
                                        beta_gl, hyd_auxvar%pres(cpid), &
                                        option)
      endif

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      call HydrateBrineSaturationPressure(hyd_auxvar%temp, &
                                         x_salt_dissolved, &
                                         hyd_auxvar%pres(spid))
      cell_pressure = max(hyd_auxvar%pres(gid),hyd_auxvar%pres(spid))
      ! Brine density
      call HydrateBrineDensity(hyd_auxvar%temp, cell_pressure, &
                           x_salt_dissolved, hyd_auxvar%den_kg(pbid), option)
      ! Brine vapor pressure
      call HydrateVaporPressureBrine(hyd_auxvar%temp, hyd_auxvar%pres(spid), &
                                   hyd_auxvar%pres(cpid), &
                                   hyd_auxvar%den_kg(pbid), &
                                   x_salt_dissolved, hyd_auxvar%pres(rvpid))

      ! Pure water density
      xsl = x_salt_dissolved
      call HydrateWaterDensity(hyd_auxvar%temp,hyd_auxvar%pres(rvpid), &
                            TWO_INTEGER,hyd_auxvar%den_kg(pwid), &
                            den_steam_kg,option)
      pva = max(hyd_auxvar%pres(gid) - hyd_auxvar%pres(rvpid), 0.d0)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(rvpid)
      call HydrateEquilibrate(hyd_auxvar%temp,cell_pressure, &
                           global_auxvar%istate, hyd_auxvar%sat(hid), &
                           hyd_auxvar%pres(apid), &
                           hyd_auxvar%pres(vpid), &
                           hyd_auxvar%pres(spid), &
                           hyd_auxvar%pres(rvpid), &
                           xag, xwg, xal, xsl, xwl, &
                           xmolag, xmolwg, xmolal, xmolsl, xmolwl, &
                           characteristic_curves, material_auxvar, option)

      ! Update mass fractions
      hyd_auxvar%xmass(acid,lid) = xal
      hyd_auxvar%xmass(wid,lid) = xwl
      hyd_auxvar%xmass(sid,lid) = xsl
      hyd_auxvar%xmass(acid,gid) = xag
      hyd_auxvar%xmass(wid,gid) = xwg

      ! Update mole fractions
      hyd_auxvar%xmol(acid,lid) = xmolal
      hyd_auxvar%xmol(wid,lid) = xmolwl
      hyd_auxvar%xmol(sid,lid) = xmolsl
      hyd_auxvar%xmol(acid,gid) = xmolag
      hyd_auxvar%xmol(wid,gid) = xmolwg
      
    case(HA_STATE)
!     ********* Hydrate & Aqueous State (HA) ********************************
!     Primary variables: Pg, Sh, T
!

      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid)

      T_temp = hyd_auxvar%temp - Tf_ice

      if (T_temp <= 0.d0) then
        ! Clausius-Clayperon equation
        Pc = -(T_temp) * (L_ICE * ICE_DENSITY * 1.d6) / (Tf_ice + 273.15d0)
        ! Get the corresponding liquid saturation
        call HydrateComputeSatHysteresis(characteristic_curves, &
                                    Pc, &
                                    hyd_auxvar%sl_min, &
                                    1.d0, hyd_auxvar%den_kg(lid), &
                                    hyd_auxvar%sat(lid), &
                                    hyd_auxvar%sat(tgid), &
                                    option)
        hyd_auxvar%sat(iid) = max(0.d0, 1.d0 - hyd_auxvar%sat(lid) - &
                              hyd_auxvar%sat(hid))
      else
        hyd_auxvar%sat(lid) = 1.d0 - hyd_auxvar%sat(hid)
        hyd_auxvar%sat(iid) = 0.d0
      endif

      hyd_auxvar%sat(gid) = 0.d0

      hyd_auxvar%pres(cpid) = 0.d0

      Pc = 0.d0

      ! MAN: Replace with dissolved salt concentration:
      hyd_auxvar%m_salt(1) = hydrate_xmass_nacl
      ! kg NaCl/kg liquid
      call HydrateComputeSaltSolubility(hyd_auxvar%temp, salt_solubility)
      ! Dissolved salt mass fraction
      x_salt_dissolved = min(hyd_auxvar%m_salt(1),salt_solubility)
      call HydrateBrineSaturationPressure(hyd_auxvar%temp, &
                                         x_salt_dissolved, &
                                         hyd_auxvar%pres(spid))
      ! Brine density
      call HydrateBrineDensity(hyd_auxvar%temp, hyd_auxvar%pres(lid), &
                            x_salt_dissolved, hyd_auxvar%den_kg(pbid), option)
      ! Brine vapor pressure
      call HydrateVaporPressureBrine(hyd_auxvar%temp, hyd_auxvar%pres(spid), &
                                   Pc, hyd_auxvar%den_kg(pbid), &
                                   x_salt_dissolved, hyd_auxvar%pres(rvpid))

      ! Pure water density
      call HydrateWaterDensity(hyd_auxvar%temp,hyd_auxvar%pres(rvpid), &
                            TWO_INTEGER,hyd_auxvar%den_kg(pwid), &
                            den_steam_kg,option)
      pva = max(hyd_auxvar%pres(lid) - hyd_auxvar%pres(rvpid), 0.d0)
      xsl = x_salt_dissolved
      call HydrateEquilibrate(hyd_auxvar%temp,hyd_auxvar%pres(lid), &
                           global_auxvar%istate, hyd_auxvar%sat(hid), &
                           hyd_auxvar%pres(apid), &
                           hyd_auxvar%pres(vpid), &
                           hyd_auxvar%pres(spid), &
                           hyd_auxvar%pres(rvpid), &
                           xag, xwg, xal, xsl, xwl, &
                           xmolag, xmolwg, xmolal, xmolsl, xmolwl, &
                           characteristic_curves, material_auxvar, option)

      ! Update mass fractions
      hyd_auxvar%xmass(acid,lid) = xal
      hyd_auxvar%xmass(wid,lid) = xwl
      hyd_auxvar%xmass(sid,lid) = xsl
      hyd_auxvar%xmass(acid,gid) = xag
      hyd_auxvar%xmass(wid,gid) = xwg

      ! Update mole fractions
      hyd_auxvar%xmol(acid,lid) = xmolal
      hyd_auxvar%xmol(wid,lid) = xmolwl
      hyd_auxvar%xmol(sid,lid) = xmolsl
      hyd_auxvar%xmol(acid,gid) = xmolag
      hyd_auxvar%xmol(wid,gid) = xmolwg

    case(HI_STATE)
!     ********* Hydrate & Ice State (HI) ********************************
!     Primary variables: Pg, Sh, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      T_temp = hyd_auxvar%temp - Tf_ice

      hyd_auxvar%sat(hid) = min(max(hyd_auxvar%sat(hid),0.d0),1.d0)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(iid) = 1.d0 - hyd_auxvar%sat(hid)

      call HydrateComputeEffectiveSat(hyd_auxvar,g_sat_eff,&
                                      h_sat_eff,i_sat_eff)

      call HydratePE(T_temp, h_sat_eff, PE_hyd, &
                dP, characteristic_curves, material_auxvar,option)

      hyd_auxvar%pres(apid) = PE_hyd
      call EOSGasHenry(T_temp,hyd_auxvar%pres(spid),K_H_tilde, &
                       ierr)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

    case(GI_STATE)
!     ********* Gas & Ice State (GI) ********************************
!     Primary variables: Pg, Si, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      T_temp = hyd_auxvar%temp - Tf_ice

      hyd_auxvar%sat(iid) = min(max(hyd_auxvar%sat(iid),0.d0),1.d0)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 1.d0 - hyd_auxvar%sat(iid)
      hyd_auxvar%sat(hid) = 0.d0

      call EOSWaterSaturationPressure(T_temp, &
                                          hyd_auxvar%pres(spid),ierr)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(vpid)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = 0.d0
      hyd_auxvar%xmol(wid,lid) = 1.d0
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

    case(AI_STATE)
!     ********* Aqueous & Ice State (AI) ********************************
!     Primary variables: Pl, Xma, Sl
!
      hyd_auxvar%pres(lid) = x(HYDRATE_LIQUID_PRESSURE_DOF)
      hyd_auxvar%xmol(acid,lid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%sat(lid) = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%xmol(acid,lid) = max(0.d0,hyd_auxvar%xmol(acid,lid))

      hyd_auxvar%sat(lid) = max(0.d0,min(1.d0,hyd_auxvar%sat(lid)))

      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 1.d0 - hyd_auxvar%sat(lid)

      call CalcFreezingTempDepression(hyd_auxvar%sat(lid),Tf_ice, &
                                      characteristic_curves, dTf,option)

      hyd_auxvar%temp = Tf_ice-dTf

      T_temp = hyd_auxvar%temp - Tf_ice

      call EOSWaterSaturationPressure(T_temp, &
                                          hyd_auxvar%pres(spid),ierr)
      call EOSGasHenry(T_temp,hyd_auxvar%pres(spid),K_H_tilde, &
                       ierr)

      hyd_auxvar%pres(gid) = max(hyd_auxvar%pres(lid),hyd_auxvar%pres(spid))
      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(apid) = K_H_tilde*hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(lid) - hyd_auxvar%pres(apid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(acid,gid) = 0.d0
      hyd_auxvar%xmol(wid,gid) = 0.d0

    case(HGA_STATE)
!     ********* Hydrate, Gas, & Aqueous State (HGA) **************************
!     Primary variables: Sl, Sh, T
!
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      !if (hyd_auxvar%sat(lid) + hyd_auxvar%sat(hid) > 1.d0) then
      !  hyd_auxvar%sat(lid) = 1.d0 - hyd_auxvar%sat(hid)
      !endif

      T_temp = hyd_auxvar%temp - Tf_ice

      hyd_auxvar%sat(gid) = 1.d0 - hyd_auxvar%sat(lid) - hyd_auxvar%sat(hid)
      hyd_auxvar%sat(gid) = max(hyd_auxvar%sat(gid),0.d0)
      hyd_auxvar%sat(iid) = 0.d0

      call HydrateComputeEffectiveSat(hyd_auxvar,g_sat_eff,&
                                    h_sat_eff,i_sat_eff)

      call HydratePE(T_temp, h_sat_eff, PE_hyd, dP,&
                      characteristic_curves, material_auxvar,option)
      call HydrateComputeSaltSolubility(hyd_auxvar%temp, salt_solubility)
      x_salt_dissolved = min(hyd_auxvar%m_salt(1),salt_solubility)
      call HydrateComputeSurfaceTension(hyd_auxvar%temp, &
                                     x_salt_dissolved, sigma)

      beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma
      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        sl_temp = 1.d0 - g_sat_eff
        call HydrateComputePcHysteresis(characteristic_curves, &
                                   sl_temp, &
                                   hyd_auxvar%sat(tgid), &
                                   beta_gl,hyd_auxvar%pres(cpid), option)
      endif

      hyd_auxvar%pres(apid) = PE_hyd
      
      call HydrateBrineSaturationPressure(hyd_auxvar%temp, &
                                         x_salt_dissolved, &
                                         hyd_auxvar%pres(spid))
      cell_pressure = max(hyd_auxvar%pres(apid),hyd_auxvar%pres(spid))
      ! Brine density
      call HydrateBrineDensity(hyd_auxvar%temp, cell_pressure, &
                           x_salt_dissolved, hyd_auxvar%den_kg(pbid), option)
      ! Brine vapor pressure
      call HydrateVaporPressureBrine(hyd_auxvar%temp, hyd_auxvar%pres(spid), &
                                   hyd_auxvar%pres(cpid), &
                                   hyd_auxvar%den_kg(pbid), &
                                   x_salt_dissolved, hyd_auxvar%pres(rvpid))

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(rvpid)
      hyd_auxvar%pres(gid) = hyd_auxvar%pres(apid) + hyd_auxvar%pres(vpid)
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      ! Pure water density
      xsl = x_salt_dissolved
      call HydrateWaterDensity(hyd_auxvar%temp,hyd_auxvar%pres(rvpid), &
                            TWO_INTEGER,hyd_auxvar%den_kg(pwid), &
                            den_steam_kg,option)
      pva = max(hyd_auxvar%pres(gid) - hyd_auxvar%pres(rvpid), 0.d0)
      call HydrateEquilibrate(hyd_auxvar%temp,cell_pressure, &
                           global_auxvar%istate, h_sat_eff, &
                           hyd_auxvar%pres(apid), &
                           hyd_auxvar%pres(vpid), &
                           hyd_auxvar%pres(spid), &
                           hyd_auxvar%pres(rvpid), &
                           xag, xwg, xal, xsl, xwl, &
                           xmolag, xmolwg, xmolal, xmolsl, xmolwl, &
                           characteristic_curves, material_auxvar, option)

      ! Update mass fractions
      hyd_auxvar%xmass(acid,lid) = xal
      hyd_auxvar%xmass(wid,lid) = xwl
      hyd_auxvar%xmass(sid,lid) = xsl
      hyd_auxvar%xmass(acid,gid) = xag
      hyd_auxvar%xmass(wid,gid) = xwg

      ! Update mole fractions
      hyd_auxvar%xmol(acid,lid) = xmolal
      hyd_auxvar%xmol(wid,lid) = xmolwl
      hyd_auxvar%xmol(sid,lid) = xmolsl
      hyd_auxvar%xmol(acid,gid) = xmolag
      hyd_auxvar%xmol(wid,gid) = xmolwg
    case(HAI_STATE)
!     ********* Hydrate, Aqueous, & Ice State (HAI) **************************
!     Primary variables: Pg, Sl, Si
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = max(0.d0,min(1.d0,hyd_auxvar%sat(lid)))
      hyd_auxvar%sat(iid) = min(max(0.d0,hyd_auxvar%sat(iid)),1.d0)

      if (hyd_auxvar%sat(lid) + hyd_auxvar%sat(iid) > 1.d0) then
        sat_temp = hyd_auxvar%sat(lid) + hyd_auxvar%sat(iid)
        hyd_auxvar%sat(lid) = hyd_auxvar%sat(lid)/sat_temp
        hyd_auxvar%sat(iid) = hyd_auxvar%sat(iid)/sat_temp
      endif

      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 1.d0 - hyd_auxvar%sat(lid) - hyd_auxvar%sat(iid)

      hyd_auxvar%sat(hid) = min(max(0.d0,hyd_auxvar%sat(hid)),1.d0)


      call HydrateComputeEffectiveSat(hyd_auxvar,g_sat_eff,&
                                    h_sat_eff,i_sat_eff)

      call CalcFreezingTempDepression(1.d0-i_sat_eff, Tf_ice, &
                                     characteristic_curves,dTf,option)

      hyd_auxvar%temp = Tf_ice-dTf

      T_temp = hyd_auxvar%temp - Tf_ice

      call EOSWaterSaturationPressure(T_temp, &
                                          hyd_auxvar%pres(spid),ierr)
      call HydratePE(T_temp,h_sat_eff, PE_hyd, dP,&
           characteristic_curves, material_auxvar,option)
      call EOSGasHenry(T_temp,hyd_auxvar%pres(spid),K_H_tilde, &
                       ierr)
      hyd_auxvar%pres(cpid) = 0.d0

      hyd_auxvar%pres(apid) = PE_hyd !hyd_auxvar%pres(gid)
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = 0.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0
    case(HGI_STATE)
!     ********* Hydrate, Gas, & Ice State (HGI) ******************************
!     Primary variables: Si, Sh, T
!
      hyd_auxvar%sat(iid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      T_temp = hyd_auxvar%temp - Tf_ice

      hyd_auxvar%sat(iid) = min(max(hyd_auxvar%sat(iid),0.d0),1.d0)
      hyd_auxvar%sat(hid) = min(max(hyd_auxvar%sat(hid),0.d0),1.d0)

      if (hyd_auxvar%sat(iid) + hyd_auxvar%sat(hid) > 1.d0) then
        sat_temp = hyd_auxvar%sat(iid) + hyd_auxvar%sat(hid)
        hyd_auxvar%sat(iid) = hyd_auxvar%sat(iid)/sat_temp
        hyd_auxvar%sat(hid) = hyd_auxvar%sat(hid)/sat_temp
      endif

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 1.d0 - hyd_auxvar%sat(hid) - hyd_auxvar%sat(iid)

      call HydrateComputeEffectiveSat(hyd_auxvar,g_sat_eff,&
                                    h_sat_eff,i_sat_eff)

      call HydratePE(T_temp,h_sat_eff, PE_hyd, dP, &
          characteristic_curves, material_auxvar, option)

      hyd_auxvar%pres(apid) = PE_hyd

      call EOSWaterSaturationPressure(T_temp, &
                                    hyd_auxvar%pres(spid),ierr)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(gid) = hyd_auxvar%pres(apid) + hyd_auxvar%pres(vpid)

      hyd_auxvar%xmol(acid,lid) = 0.d0
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

    case(GAI_STATE)
!     ********* Gas, Aqueous, & Ice State (GAI) ******************************
!     Primary variables: Pg, Sg, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      T_temp = hyd_auxvar%temp - Tf_ice

      ! Clausius-Clayperon equation
      Pc = -(T_temp) * (L_ICE * ICE_DENSITY * 1.d6) / (Tf_ice + T273K)
      ! Get the corresponding liquid saturation
      call HydrateComputeSatHysteresis(characteristic_curves, &
                                    Pc, &
                                    hyd_auxvar%sl_min, &
                                    1.d0, hyd_auxvar%den_kg(lid), &
                                    hyd_auxvar%sat(lid), &
                                    hyd_auxvar%sat(tgid), &
                                    option)

      hyd_auxvar%sat(gid) = min(max(hyd_auxvar%sat(gid),0.d0),1.d0)
      hyd_auxvar%sat(lid) = min(max(hyd_auxvar%sat(lid),0.d0),1.d0)

      if (hyd_auxvar%sat(lid) + hyd_auxvar%sat(gid) > 1.d0) then
        sat_temp = hyd_auxvar%sat(lid) + hyd_auxvar%sat(gid)
        hyd_auxvar%sat(lid) = hyd_auxvar%sat(lid)/sat_temp
        hyd_auxvar%sat(gid) = hyd_auxvar%sat(gid)/sat_temp
      endif

      hyd_auxvar%sat(iid) = 1.d0 - hyd_auxvar%sat(lid) - hyd_auxvar%sat(gid)
      hyd_auxvar%sat(hid) = 0.d0

      call HydrateComputeEffectiveSat(hyd_auxvar,g_sat_eff,&
                                    h_sat_eff,i_sat_eff)

      call EOSWaterSaturationPressure(T_temp, &
                                    hyd_auxvar%pres(spid),ierr)

      call EOSGasHenry(T_temp,hyd_auxvar%pres(spid),K_H_tilde, &
                       ierr)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(vpid)

      call HydrateComputeSaltSolubility(hyd_auxvar%temp, salt_solubility)
      if (hyd_auxvar%m_salt(2) > epsilon) then
        x_salt_dissolved = salt_solubility
      else
        x_salt_dissolved = 0.d0
      endif
      hyd_auxvar%xmass(sid,lid) = x_salt_dissolved
      call HydrateComputeSurfaceTension(hyd_auxvar%temp, &
                                     x_salt_dissolved, sigma)
      ! MAN: check the reference surface tension
      !MAN: hyd_auxvar%sat(tgid) should be small?
      beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma
      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call HydrateComputePcHysteresis(characteristic_curves, &
                                   1.d0-g_sat_eff, &
                                   hyd_auxvar%sat(tgid), &
                                   beta_gl,hyd_auxvar%pres(cpid), option)
      endif

      !IFT calculation
      sigma=1.d0
      if (hydrate_compute_surface_tension) then
       call EOSWaterSurfaceTension(T_temp,sigma)
      endif
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid)*sigma

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

    case(HGAI_STATE)
!     ********* 4-Phase (HGAI) ********************************
!     Primary variables: Sl, Sg, Si
!
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = min(max(hyd_auxvar%sat(lid),0.d0),1.d0)
      hyd_auxvar%sat(gid) = min(max(hyd_auxvar%sat(gid),0.d0),1.d0)
      hyd_auxvar%sat(iid) = min(max(hyd_auxvar%sat(iid),0.d0),1.d0)

      if (hyd_auxvar%sat(lid) + hyd_auxvar%sat(gid) + hyd_auxvar%sat(iid) > &
          1.d0) then
        sat_temp = hyd_auxvar%sat(lid) + hyd_auxvar%sat(gid) + &
                   hyd_auxvar%sat(iid)
        hyd_auxvar%sat(lid) = hyd_auxvar%sat(lid)/sat_temp
        hyd_auxvar%sat(gid) = hyd_auxvar%sat(gid)/sat_temp
        hyd_auxvar%sat(iid) = hyd_auxvar%sat(iid)/sat_temp
      endif
      hyd_auxvar%sat(hid) = 1.d0 - hyd_auxvar%sat(lid) - hyd_auxvar%sat(gid) &
                            - hyd_auxvar%sat(iid)

      call HydrateComputeEffectiveSat(hyd_auxvar,g_sat_eff,&
                                    h_sat_eff,i_sat_eff)

      call CalcFreezingTempDepression(1.d0-i_sat_eff, Tf_ice, &
                                      characteristic_curves,dTf,option)

      hyd_auxvar%temp = Tf_ice - dTf

      T_temp = hyd_auxvar%temp - Tf_ice

      call HydratePE(T_temp,h_sat_eff, PE_hyd, dP, &
          characteristic_curves, material_auxvar, option)
      hyd_auxvar%pres(apid) = PE_hyd

      call EOSWaterSaturationPressure(T_temp, &
                                    hyd_auxvar%pres(spid),ierr)

      call EOSGasHenry(T_temp,hyd_auxvar%pres(spid),K_H_tilde, &
                       ierr)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(gid) = hyd_auxvar%pres(apid) + hyd_auxvar%pres(vpid)

      call HydrateComputeSaltSolubility(hyd_auxvar%temp, salt_solubility)
      if (hyd_auxvar%m_salt(2) > epsilon) then
        x_salt_dissolved = salt_solubility
      else
        x_salt_dissolved = 0.d0
      endif
      hyd_auxvar%xmass(sid,lid) = x_salt_dissolved
      call HydrateComputeSurfaceTension(hyd_auxvar%temp, &
                                     x_salt_dissolved, sigma)
      ! MAN: check the reference surface tension
      !MAN: hyd_auxvar%sat(tgid) should be small?
      beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma
      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call HydrateComputePcHysteresis(characteristic_curves, &
                                   1.d0 - g_sat_eff, &
                                   hyd_auxvar%sat(tgid), &
                                   beta_gl,hyd_auxvar%pres(cpid), option)
      endif
      !IFT calculation
      sigma=1.d0
      if (hydrate_compute_surface_tension) then
       call EOSWaterSurfaceTension(T_temp,sigma)
      endif
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid)*sigma

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)


    case default
      write(option%io_buffer,*) global_auxvar%istate
      option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
        ') not recognized in HydrateAuxVarCompute.'
      call PrintErrMsgByRank(option)

  end select

  cell_pressure = max(hyd_auxvar%pres(lid),hyd_auxvar%pres(gid), &
                      hyd_auxvar%pres(spid))
  hyd_auxvar%xmass(acid,gid) = xag
  hyd_auxvar%xmass(wid,gid) = 1.d0 - xag
  hyd_auxvar%xmol(acid,gid) = xmolag
  hyd_auxvar%xmol(wid,gid) = 1.d0 - xmolag

  ! calculate effective porosity as a function of pressure
  if (option%iflag /= HYDRATE_UPDATE_FOR_BOUNDARY) then
    dpor_dp = 0.d0
    hyd_auxvar%effective_porosity = material_auxvar%porosity_base
    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                hyd_auxvar%effective_porosity,dpor_dp)
    endif
    if (option%iflag /= HYDRATE_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = hyd_auxvar%effective_porosity
    endif

  endif

  ! Gas phase density
  call EOSGasDensity(hyd_auxvar%temp,pva, &
                     den_mol,drho_dT,drho_dP,ierr)
  hyd_auxvar%den_kg(apid) = den_mol * hydrate_fmw_comp(2)
  hyd_auxvar%den_kg(gid) = hyd_auxvar%xmass(acid,gid) * &
                            hyd_auxvar%den_kg(apid) + &
                            hyd_auxvar%xmass(wid,gid) * &
                            den_steam_kg
  select case (hydrate_former)
    case(HYDRATE_FORMER_CO2)
      ! Gas phase viscosity
      call HydrateViscosityWater(hyd_auxvar%temp,hyd_auxvar%pres(rvpid), &
                                 den_steam_kg,visc_water,option)
      call HydrateViscosityCO2(hyd_auxvar%temp,hyd_auxvar%den_kg(apid), &
                               visc_a)
      call HydrateViscosityGas(visc_water,visc_a,hyd_auxvar%xmol(wid,gid), &
                               hyd_auxvar%xmol(acid,gid),hyd_auxvar%visc(gid))
      ! Liquid phase density (including air)
      call HydrateDensityCompositeLiquid(hyd_auxvar%temp,hyd_auxvar%den_kg(pbid), &
                                  hyd_auxvar%xmass(acid,lid), &
                                  hyd_auxvar%den_kg(lid))
      ! Liquid phase viscosity
      call HydrateWaterDensity(hyd_auxvar%temp, cell_pressure, ONE_INTEGER, &
                               hyd_auxvar%den_kg(pwid), den_steam, &
                               option)
      call HydrateViscosityWater(hyd_auxvar%temp,cell_pressure, &
                                 hyd_auxvar%den_kg(pwid),visc_water,option)
      call HydrateViscosityBrine(hyd_auxvar%temp, hyd_auxvar%xmass(sid,lid), &
                                 visc_water, visc_brine)
      call HydrateViscosityLiquid(hyd_auxvar%xmol(acid,lid), visc_brine, &
                                  visc_a, hyd_auxvar%visc(lid))
    case default
      if (.not.option%flow%density_depends_on_salinity) then
        call EOSWaterDensity(T_temp,cell_pressure, &
                             hyd_auxvar%den_kg(lid),hyd_auxvar%den(lid),ierr)
      else
        aux(1) = hyd_auxvar%xmass(sid,lid)
        call EOSWaterDensityExt(T_temp,cell_pressure,aux, &
                                  hyd_auxvar%den_kg(lid),hyd_auxvar%den(lid),ierr)
      endif
      call EOSWaterViscosity(T_temp,cell_pressure, &
                            hyd_auxvar%pres(spid),hyd_auxvar%visc(lid),ierr)
      call EOSGasViscosity(hyd_auxvar%temp,hyd_auxvar%pres(apid), &
                           hyd_auxvar%pres(gid),hyd_auxvar%den_kg(apid), & 
                           hyd_auxvar%visc(gid),ierr)
      visc_a = hyd_auxvar%visc(gid)
  end select

  ! CO2-water surface tension
  ! MAN: doesn't do anything here right now
  call HydrateComputeSurfaceTension(hyd_auxvar%temp,hyd_auxvar%xmass(sid,lid), &
                                 sigma)
  beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma

  if (hydrate_eff_sat_scaling) then
    l_sat_eff = hyd_auxvar%sat(lid)/(hyd_auxvar%sat(lid)+ &
                hyd_auxvar%sat(gid))
    g_sat_eff = 1.d0 - l_sat_eff
  else
    l_sat_eff = hyd_auxvar%sat(lid)
    g_sat_eff = hyd_auxvar%sat(gid)
  endif
  
  ! Relative Permeability
  if (hyd_auxvar%sat(lid) > 0.d0) then
    if (hyd_auxvar%sat(lid) >= 1.d0) then
      hyd_auxvar%kr(lid) = 1.d0
    else
      call characteristic_curves%liq_rel_perm_function% &
               RelativePermeability(l_sat_eff,hyd_auxvar%kr(lid), &
                                    dkrl_dsatl,option)
    endif
  else
    hyd_auxvar%kr(lid) = 0.d0
  endif
  hyd_auxvar%kr(lid) = min(max(hyd_auxvar%kr(lid),1.d-24),1.d0)
  hyd_auxvar%mobility(lid) = hyd_auxvar%kr(lid) / hyd_auxvar%visc(lid)

  if (hyd_auxvar%sat(gid) > 0.d0) then
    if (hyd_auxvar%sat(gid) >=1.d0) then
      hyd_auxvar%kr(gid) = 1.d0
    else
      call characteristic_curves%gas_rel_perm_function% &
           RelativePermeability(1.d0 - g_sat_eff,hyd_auxvar%kr(gid),dkrg_dsatl,option)
      hyd_auxvar%kr(gid) = max(0.d0,hyd_auxvar%kr(gid))
    endif
  else
    hyd_auxvar%kr(gid) = 0.d0
  endif
  hyd_auxvar%mobility(gid) = hyd_auxvar%kr(gid) / hyd_auxvar%visc(gid)
  
  ! Convert to molar density: liquid
  mw_mix = hyd_auxvar%xmol(wid,lid) * hydrate_fmw_comp(1) + &
          hyd_auxvar%xmol(acid,lid) * hydrate_fmw_comp(2) + &
          hyd_auxvar%xmol(sid,lid) * hydrate_fmw_comp(3)
  hyd_auxvar%den(lid) = hyd_auxvar%den_kg(lid) / mw_mix

  ! Convert to molar density: gas
  mw_mix = hyd_auxvar%xmol(wid,gid) * hydrate_fmw_comp(1) + &
          hyd_auxvar%xmol(acid,gid) * hydrate_fmw_comp(2) + &
          hyd_auxvar%xmol(sid,gid) * hydrate_fmw_comp(3)
  hyd_auxvar%den(gid) = hyd_auxvar%den_kg(gid) / mw_mix

  ! Tortuosity
  call HydrateTortuosity(hyd_auxvar%sat(lid), hyd_auxvar%sat(gid), &
                      hyd_auxvar%effective_porosity, &
                      hyd_auxvar%tortuosity(lid), &
                      hyd_auxvar%tortuosity(gid))
  ! Update Diffusivities: water vapor, dissolved CO2, dissolved salt
  call HydrateDiffusionCoeff(hyd_auxvar%temp, cell_pressure, &
                          hyd_auxvar%xmass(sid,lid), &
                          hyd_auxvar%visc(lid), &
                          hydrate_parameter, option)
  call HydrateComputeEffectiveDiffusion(hydrate_parameter, hyd_auxvar, option)

  ! Precipitate salt
  hyd_auxvar%xmass(sid,pid) = 1.d0
  hyd_auxvar%xmol(sid,pid) = 1.d0

  ! Salt precipitate density and saturation
  call HydrateComputeSaltDensity(hyd_auxvar%temp, cell_pressure, &
                              hyd_auxvar%den_kg(pid))
  if (global_auxvar%istate == G_STATE .or. &
      hyd_auxvar%sat(lid) == 0.d0) then
    hyd_auxvar%sat(pid) = hyd_auxvar%m_salt(2) / (hyd_auxvar%den_kg(pid) * &
                           hyd_auxvar%effective_porosity)
    !MAN: not sure why we need this:
    hyd_auxvar%m_salt(1) = hyd_auxvar%m_salt(2) * hyd_auxvar%den_kg(pbid) * &
                            epsilon * hyd_auxvar%effective_porosity
  else
    hyd_auxvar%sat(pid) = max(hyd_auxvar%m_salt(1) - salt_solubility, &
                           0.d0) * hyd_auxvar%den_kg(pbid) * &
                           hyd_auxvar%sat(lid) / &
                           hyd_auxvar%den_kg(pid)
    hyd_auxvar%m_salt(2) = hyd_auxvar%m_salt(1) * hyd_auxvar%den_kg(pbid) * &
                            hyd_auxvar%sat(lid) * &
                            hyd_auxvar%effective_porosity
  endif

  ! Permeability and porosity reduction with salt precipitate effects
  call HydrateScalePermPhi(hyd_auxvar, material_auxvar, global_auxvar, option)

  hyd_auxvar%mobility(lid) = hyd_auxvar%kr(lid) / hyd_auxvar%visc(lid)
  hyd_auxvar%mobility(gid) = hyd_auxvar%kr(gid) / hyd_auxvar%visc(gid)

  ! Thermal properties
  call EOSWaterEnthalpy(T_temp,cell_pressure,hw,ierr)

  hyd_auxvar%H(lid) = hw * 1.d-6 ! J/kmol -> MJ/kmol
  ! MJ/kmol comp
  hyd_auxvar%U(lid) = hyd_auxvar%H(lid) - &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        (cell_pressure / hyd_auxvar%den(lid) * &
                        1.d-6)
  if (global_auxvar%istate /= L_STATE) then
    water_vapor_pressure = hyd_auxvar%pres(vpid)
    call EOSGasDensityEnergy(T_temp,hyd_auxvar%pres(apid),den_air, &
                               h_air,u_air,ierr)
    h_air = h_air * 1.d-6 ! J/kmol -> MJ/kmol
    u_air = u_air * 1.d-6 ! J/kmol -> MJ/kmol
    if (water_vapor_pressure > 0.d0) then
      call EOSWaterSteamDensityEnthalpy(T_temp,water_vapor_pressure, &
                                        den_kg_water_vapor,den_water_vapor, &
                                        h_water_vapor,ierr)
      u_water_vapor = h_water_vapor - &
                    ! Pa / kmol/m^3 = J/kmol
                    water_vapor_pressure / den_water_vapor
    else
      h_water_vapor = 0.d0
      u_water_vapor = 0.d0
    endif
    
    
    h_water_vapor = h_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    u_water_vapor = u_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    xmol_air_in_gas = hyd_auxvar%xmol(acid,gid)
    xmol_water_in_gas = hyd_auxvar%xmol(wid,gid)

    ! MJ/kmol
    hyd_auxvar%U(gid) = xmol_water_in_gas * u_water_vapor + &
                        xmol_air_in_gas * u_air
    Hg_mixture_fractioned = xmol_water_in_gas*h_water_vapor + &
                            xmol_air_in_gas*h_air
    hyd_auxvar%H(gid) = hyd_auxvar%U(gid) + &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        hyd_auxvar%pres(gid)/hyd_auxvar%den(gid) * 1.d-6

  endif

  ! Hydrate phase properties
  call EOSHydrateEnthalpy(T_temp, H_hyd)
  hyd_auxvar%U(hid) = H_hyd !- cell_pressure/hyd_auxvar%den(hid)*1.d-6
  hyd_auxvar%H(hid) = H_hyd
  hyd_auxvar%mobility(hid) = 0.d0

  ! Ice phase properties
  call EOSWaterInternalEnergyIce(T_temp, U_ice, du_ice_dT, du_ice_dP,ierr)
  U_ice = U_ice * 1.d-3
  hyd_auxvar%xmol(wid,iid) = 1.d0
  hyd_auxvar%xmol(gid,iid) = 0.d0
  hyd_auxvar%den(iid) = ICE_DENSITY
  hyd_auxvar%den_kg(iid) = ICE_DENSITY_KG
  if (hydrate_no_ice_density_change) then
    hyd_auxvar%den(iid) = hyd_auxvar%den(lid)
    hyd_auxvar%den_kg(iid) = hyd_auxvar%den_kg(lid)
  else
    hyd_auxvar%den(iid) = ICE_DENSITY
    hyd_auxvar%den_kg(iid) = ICE_DENSITY_KG
  endif
  hyd_auxvar%U(iid) = U_ice
  hyd_auxvar%H(iid) = U_ice
  hyd_auxvar%mobility(iid) = 0.d0

  hyd_auxvar%srl = characteristic_curves%gas_rel_perm_function%sr
  hyd_auxvar%srg = characteristic_curves%gas_rel_perm_function%srg

end subroutine HydrateAuxVarCompute

! ************************************************************************** !

subroutine HydrateEOSGasError(natural_id,ierr,hyd_auxvar,option)

  !
  !  HydrateEOSGasError: Elaborates when variables exceeds the bounds of
  !                      the equation of state.
  !
  !  Author: Heeho Park
  !  Date: 5/29/19
  !

  use Option_module
  use String_module
  use EOS_Gas_module

  implicit none

  PetscInt :: natural_id
  PetscInt :: ierr
  type(hydrate_auxvar_type) :: hyd_auxvar
  type(option_type) :: option


  call PrintMsgByCell(option,natural_id, &
                      'Error in HydrateAuxVarCompute->EOSGasHenry')
  if (ierr == EOS_GAS_TEMP_BOUND_EXCEEDED) then
    option%io_buffer = 'Temperature at cell ID ' // trim(StringWrite(natural_id)) // &
                               ' exceeds the equation of state temperature bound with ' // &
                               trim(StringWrite(hyd_auxvar%temp)) // ' [C].'
    call PrintErrMsgByRank(option)
    hydrate_high_temp_ts_cut = PETSC_TRUE
  endif


end subroutine HydrateEOSGasError

! ************************************************************************** !

subroutine HydrateAuxVarUpdateState(x,hyd_auxvar,global_auxvar, &
                                    material_auxvar, &
                                    characteristic_curves,hydrate_parameter, &
                                    natural_id,option)

  !
  ! Decides on state changes and adds epsilons to new primary variables
  ! accordingly.
  !
  ! Author: Michael Nole
  ! Date: 01/28/18
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  class(characteristic_curves_type) :: characteristic_curves
  type(hydrate_auxvar_type) :: hyd_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(hydrate_parameter_type), pointer :: hydrate_parameter

  PetscReal, parameter :: epsilon = 1.d-14
  PetscReal, parameter :: eps_sl = 1.d-4
  PetscReal, parameter :: peta = 1.d-1

  PetscReal :: liq_epsilon, gas_epsilon, hyd_epsilon, two_phase_epsilon
  PetscReal :: ha_epsilon
  PetscReal :: x(option%nflowdof)
  PetscReal :: PE_hyd, dP, Tf_ice, dTfs, T_temp
  PetscReal :: h_sat_eff,g_sat_eff,i_sat_eff
  PetscReal :: K_H_tilde
  PetscReal :: Pc_entry, cell_pressure, sg_min, sh_min
  PetscReal :: sl_temp, sgt_temp, sg_est
  PetscReal :: Pc, Pv, Prvap, Pa, Psat
  PetscReal :: beta_gl, sgt_max
  PetscReal :: xag, xwg, xal, xsl, xwl, xmolag, xmolwg, xmolal, &
               xmolsl, xmolwl
  PetscReal :: salt_solubility, sigma
  PetscReal :: den_mol, den_a, den_brine, den_liq
  PetscReal :: drho_dP, drho_dT
  PetscInt :: apid, cpid, vpid, spid, rvpid
  PetscInt :: gid, lid, hid, iid, acid, wid, tgid, pwid, pbid, pid, sid
  PetscReal :: state_change_threshold
  PetscInt :: old_state,new_state
  PetscBool :: istatechng
  PetscReal :: dpc_dsatl
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: state_change_string, append

  if (hydrate_immiscible .or. hyd_auxvar%istatechng) return

  ierr = 0
  state_change_threshold = 0.d0

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = 3
  iid = 4
  pid = option%precipitate_phase
  pwid = option%pure_water_phase
  pbid = option%pure_brine_phase
  tgid = option%trapped_gas_phase
  spid = option%saturation_pressure_id

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  rvpid = option%reduced_vapor_pressure_id

  acid = option%air_id
  wid = option%water_id
  sid = option%salt_id

  hyd_auxvar%istate_store(PREV_IT) = global_auxvar%istate
  istatechng = PETSC_FALSE

  gas_epsilon = 0.d0
  liq_epsilon = 0.d0
  hyd_epsilon = 0.d0
  two_phase_epsilon = 0.d0

  !man: right now comparing hydrate equilib pressure to gas
  !pressure (assuming low water solubility in methane).
  !Ideally would compare to partial pressure of methane.

  if (option%flow%density_depends_on_salinity) then
    hydrate_xmass_nacl = global_auxvar%m_nacl(1)
  endif

  call HydrateComputeSaltSolubility(hyd_auxvar%temp, salt_solubility)
  if (global_auxvar%istate == G_STATE) then
    if (hyd_auxvar%m_salt(2) > epsilon) then
      xsl = salt_solubility
    else
      xsl = 0.d0
    endif
  else
    xsl = min(salt_solubility,hyd_auxvar%m_salt(1))
  endif
  call HydrateComputeSurfaceTension(hyd_auxvar%temp, &
                                 xsl, sigma)
  beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma

  ! Max effective trapped gas saturation
  ! MAN: need to check how this works with Pc function Webb extensions
  sgt_max = characteristic_curves%saturation_function%Sgt_max

  sh_min = 1.d-5
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

  !MAN: why is this here:
  if (global_auxvar%istate == ZERO_INTEGER .and. hyd_auxvar%sat(gid) &
       < 0.d0) then
    global_auxvar%istate = HA_STATE
    hyd_auxvar%sat(hid) = -1.d0 * hyd_auxvar%sat(gid)
    hyd_auxvar%sat(gid) = 0.d0
  endif

  if (global_auxvar%istate == ZERO_INTEGER) global_auxvar% &
                              istate = global_auxvar%istate
  !MAN
  hyd_auxvar%istate_store(PREV_IT) = global_auxvar%istate

  call HydrateComputeEffectiveSat(hyd_auxvar,g_sat_eff,&
                                    h_sat_eff,i_sat_eff)

  if (hydrate_xmass_nacl > 0.d0) then
    call IceSalinityOffset(hydrate_xmass_nacl,dTfs)
  else
    dTfs = 0.d0
  endif

  Tf_ice = dTfs !Bulk freezing temperature

  T_temp = hyd_auxvar%temp - Tf_ice

  call HydratePE(T_temp,h_sat_eff, PE_hyd, dP,&
          characteristic_curves, material_auxvar, option)
  call HydrateEquilibrate(hyd_auxvar%temp,hyd_auxvar%pres(lid), &
                          global_auxvar%istate, &
                          hyd_auxvar%sat(hid), &
                          Pa,Pv, &
                          hyd_auxvar%pres(spid), &
                          hyd_auxvar%pres(rvpid), &
                          xag, xwg, xal, xsl, xwl, &
                          xmolag, xmolwg, xmolal, xmolsl, xmolwl, &
                          characteristic_curves, material_auxvar, option)

  !Update State

  old_state = global_auxvar%istate

  select case(global_auxvar%istate)
    case(L_STATE)
      if (hyd_auxvar%sat(iid) == 0.d0) then !Not frozen
        cell_pressure = max(hyd_auxvar%pres(lid),hyd_auxvar%pres(vpid))
        ! Check if dissolved air exceeds solubility
        if (hyd_auxvar%xmass(acid,lid) > xal * (1.d0 + &
            state_change_threshold)) then
          call HydrateBrineDensity(hyd_auxvar%temp,cell_pressure, &
                                   xsl, den_brine, option)
          call HydrateDensityCompositeLiquid(hyd_auxvar%temp, &
                                             den_brine,xal,den_liq)
          if (hyd_auxvar%pres(lid) >= PE_hyd) then
            istatechng = PETSC_TRUE
            global_auxvar%istate = HA_STATE
          else
            ! Compute what gas saturation would be
            call EOSGasDensity(hyd_auxvar%temp,Pa, &
                    den_mol,drho_dT,drho_dP,ierr)
            den_a= den_mol * hydrate_fmw_comp(2)
            call HydrateBrineDensity(hyd_auxvar%temp,cell_pressure, &
                             xsl, den_brine, option)
            call HydrateDensityCompositeLiquid(hyd_auxvar%temp, &
                                  den_brine,xal,den_liq)
            sg_est =  (hyd_auxvar%xmass(acid,lid) - xal * (1.d0 + &
                   state_change_threshold)) * &
                   den_liq / den_a
            if (sg_est < sg_min) then
                    ! No state change
              istatechng = PETSC_FALSE
            else
              sl_temp = 1.d0 - min(sg_est, 1.d-1)
              sgt_temp = 0.d0
              !MAN: hyd_auxvar%sat(tgid) should be 0
              call characteristic_curves%saturation_function% &
                              CapillaryPressure(sl_temp,Pc,dpc_dsatl,option)
              Pc = min(Pc,Pc_entry / beta_gl + 1.d5)

              ! State has changed, so update state and one primary variable
              hyd_auxvar%pres(gid) = hyd_auxvar%pres(lid) + Pc
              hyd_auxvar%sat(gid) = 1.d0 - sl_temp

              istatechng = PETSC_TRUE
              global_auxvar%istate = GA_STATE
            endif
            istatechng = PETSC_TRUE
            global_auxvar%istate = GA_STATE
          endif
        else
          ! No state change
          istatechng = PETSC_FALSE
        endif
      else !Frozen
        if (hyd_auxvar%xmass(acid,lid) > xal * (1.d0 + &
          state_change_threshold)) then
          if (hyd_auxvar%pres(lid) >= PE_hyd) then
            istatechng = PETSC_TRUE
            global_auxvar%istate = HA_STATE
          else
            ! Compute what gas saturation would be
            call EOSGasDensity(hyd_auxvar%temp,Pa, &
                  den_mol,drho_dT,drho_dP,ierr)
            den_a= den_mol * hydrate_fmw_comp(2)
            call HydrateBrineDensity(hyd_auxvar%temp,cell_pressure, &
                                     xsl, den_brine, option)
            call HydrateDensityCompositeLiquid(hyd_auxvar%temp, &
                                  den_brine,xal,den_liq)
            sg_est =  (hyd_auxvar%xmass(acid,lid) - xal * (1.d0 + &
                       state_change_threshold)) * &
                       den_liq / den_a
            if (sg_est < sg_min) then
                    ! No state change
              istatechng = PETSC_FALSE
            else
              sl_temp = 1.d0 - min(sg_est, 1.d-1)
              sgt_temp = 0.d0
              !MAN: hyd_auxvar%sat(tgid) should be 0
              call characteristic_curves%saturation_function% &
                            CapillaryPressure(sl_temp,Pc,dpc_dsatl,option)
              Pc = min(Pc,Pc_entry / beta_gl + 1.d5)

              ! State has changed, so update state and one primary variable
              hyd_auxvar%pres(gid) = hyd_auxvar%pres(lid) + Pc
              hyd_auxvar%sat(gid) = 1.d0 - sl_temp

              istatechng = PETSC_TRUE
              global_auxvar%istate = GA_STATE
            endif
          endif
        else
          istatechng = PETSC_FALSE
        endif
      endif
    case(G_STATE)
      Pv = hyd_auxvar%pres(gid) - &
                hyd_auxvar%pres(apid)
      prvap = hyd_auxvar%pres(rvpid)

      if (pv > prvap * (1.d0 + state_change_threshold)) then
        if (hyd_auxvar%temp >= Tf_ice .and. &
            hyd_auxvar%pres(apid) < PE_hyd) then
          ! Aqueous phase appears, transition state. Update primary variables
          ! for salt mass and liquid pressure
          istatechng = PETSC_TRUE
          global_auxvar%istate = GA_STATE

          sl_temp = (pv - prvap * (1.d0 + &
                     state_change_threshold)) / (prvap * (1.d0 + &
                     state_change_threshold))
          hyd_auxvar%m_salt(1) = hyd_auxvar%m_salt(2) / &
                                (hyd_auxvar%den_kg(lid) * &
                                sl_temp * &
                                material_auxvar%porosity)
          call HydrateComputeSurfaceTension(hyd_auxvar%temp, &
                                       xsl, sigma)
          beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma
          sgt_temp = 0.d0
          call HydrateComputePcHysteresis(characteristic_curves, &
                                          hyd_auxvar%sat(lid), &
                                          hyd_auxvar%sat(tgid), &
                                          beta_gl, Pc, option)
          hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - &
                                 Pc
        elseif (hyd_auxvar%pres(apid) < PE_hyd) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = GI_STATE
          gas_epsilon = hydrate_phase_chng_epsilon
        else
          istatechng = PETSC_TRUE
          global_auxvar%istate = HG_STATE
          gas_epsilon = hydrate_phase_chng_epsilon
        endif
      else
        istatechng = PETSC_FALSE
      endif

    case(H_STATE)
      if (hyd_auxvar%pres(apid) < PE_hyd) then
      !if (hyd_auxvar%pres(gid) < PE_hyd) then

        if (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HGA_STATE
        elseif (hyd_auxvar%temp == Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HGAI_STATE
        else
          istatechng = PETSC_TRUE
          global_auxvar%istate = HGI_STATE
        endif

      else

          istatechng = PETSC_FALSE

      endif

    case(I_STATE)
      if (hyd_auxvar%temp > Tf_ice) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = L_STATE
      else
        istatechng = PETSC_FALSE
      endif

    case(GA_STATE)
      ! Compute Saturation including Hysteresis
      if (.not. hydrate_no_pc) then
        call HydrateComputeSatHysteresis(characteristic_curves, &
                                    hyd_auxvar%pres(cpid), &
                                    hyd_auxvar%sl_min, &
                                    beta_gl, hyd_auxvar%den_kg(lid), &
                                    sl_temp, &
                                    sgt_temp, &
                                    option)
        sl_temp = sl_temp + sgt_temp
      else
        sl_temp = hyd_auxvar%sat(lid)
      endif
      if (dabs(1.d0 - sl_temp) < epsilon) then
        ! Gas goes away, just liquid. Update 1 primary variable.
        istatechng = PETSC_TRUE
        global_auxvar%istate = L_STATE

        cell_pressure = hyd_auxvar%pres(lid)
        ! cell_pressure = min(cell_pressure, 1.d8)
        call HydrateBrineSaturationPressure(hyd_auxvar%temp,xsl, &
                                              Psat)
        call HydrateEquilibrate(hyd_auxvar%temp,cell_pressure, &
                               global_auxvar%istate, &
                               hyd_auxvar%sat(hid), &
                               hyd_auxvar%pres(apid), &
                               Pv, Psat, &
                               hyd_auxvar%pres(rvpid), &
                               xag, xwg, xal, xsl, xwl, &
                               xmolag, xmolwg, xmolal, xmolsl, xmolwl, &
                               characteristic_curves, material_auxvar, option)
        hyd_auxvar%xmass(acid,lid) = xal
      elseif (sl_temp < epsilon .and. &
              (1.d0 - hyd_auxvar%sat(lid)) > epsilon) then
      ! Transition to fully gas-saturated. Update 2 primary variables.
        istatechng = PETSC_TRUE
        global_auxvar%istate = G_STATE
      else
        if (hyd_auxvar%pres(apid) < PE_hyd) then
          if (hyd_auxvar%temp > Tf_ice) then
          ! No state transition
            istatechng = PETSC_FALSE
          else
            !istatechng = PETSC_TRUE
            !global_auxvar%istate = GAI_STATE
          endif
        else
          if (hyd_auxvar%temp > Tf_ice) then
            !sh_est = 1.d0 - dabs(PE_hyd / hyd_auxvar%pres(apid))
            !sh_est = min(hyd_auxvar%sat(gid)-1.d-3,sh_est)
            !hyd_auxvar%sat(hid) = sh_est
            !hyd_auxvar%sat(gid) = hyd_auxvar%sat(gid) - sh_est
            istatechng = PETSC_TRUE
            global_auxvar%istate = HGA_STATE
            two_phase_epsilon = hydrate_phase_chng_epsilon
          else
            istatechng = PETSC_TRUE
            global_auxvar%istate = HGAI_STATE
            two_phase_epsilon = hydrate_phase_chng_epsilon
          endif
        endif
      endif
    case(HG_STATE)
      if (hyd_auxvar%pres(apid) > PE_hyd) then
      !if (hyd_auxvar%pres(gid) > PE_hyd) then

        if (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(hid) > 0.d0) then
          ! istatechng = PETSC_TRUE
          ! global_auxvar%istate = H_STATE
          ! two_phase_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = G_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        endif

      else

        if (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HGA_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%temp == Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HGAI_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        else
          istatechng = PETSC_TRUE
          global_auxvar%istate = HGI_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        endif

      endif

    case(HA_STATE)
      !if (hyd_auxvar%pres(gid) > PE_hyd .and. hyd_auxvar%temp > Tf_ice) then
      if (hyd_auxvar%pres(apid) > PE_hyd .and. hyd_auxvar%temp > Tf_ice) then

        if (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(hid) > 0.d0) then
          ! istatechng = PETSC_TRUE
          ! global_auxvar%istate = H_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = L_STATE
        endif

      elseif (hyd_auxvar%temp > Tf_ice) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = HGA_STATE
        ha_epsilon = hydrate_phase_chng_epsilon

      elseif (hyd_auxvar%pres(apid) > PE_hyd) then
      !elseif (hyd_auxvar%pres(gid) > PE_hyd) then
       istatechng = PETSC_TRUE
       global_auxvar%istate = HAI_STATE

      else
       istatechng = PETSC_TRUE
       global_auxvar%istate = HGAI_STATE

      endif

    case(HI_STATE)

      if (hyd_auxvar%pres(apid) >= PE_hyd) then
      !if (hyd_auxvar%pres(gid) > PE_hyd) then

        if (hyd_auxvar%temp < Tf_ice) then
          istatechng = PETSC_FALSE
        else
          istatechng = PETSC_TRUE
          global_auxvar%istate = HAI_STATE
        endif

      else

        if (hyd_auxvar%temp < Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HGI_STATE
        else
          istatechng = PETSC_TRUE
          global_auxvar%istate = HGAI_STATE
        endif

      endif
    case(GI_STATE)
      if (hyd_auxvar%temp < Tf_ice .and. hyd_auxvar%pres(apid) < PE_hyd) then
      !if (hyd_auxvar%temp < Tf_ice .and. hyd_auxvar%pres(gid) < PE_hyd &
      !    .and. hydrate_former == HYDRATE_FORMER_CH4) then
        if (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = G_STATE
        else
          istatechng = PETSC_FALSE
          !istatechng = PETSC_TRUE
          !global_auxvar%istate = I_STATE
        endif

      elseif (hyd_auxvar%temp < Tf_ice .and. hydrate_former /= HYDRATE_FORMER_NULL) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = HGI_STATE

      elseif (hyd_auxvar%pres(apid) < PE_hyd) then
      !elseif ((hyd_auxvar%pres(gid) < PE_hyd .and. hydrate_former == HYDRATE_FORMER_CH4) &
      !        .or. (hydrate_former == HYDRATE_FORMER_NULL)) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = GA_STATE

      elseif (hydrate_former /= HYDRATE_FORMER_NULL) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = HGAI_STATE
      endif

     case(AI_STATE)
       !MAN: state transitions do not naturally enter this state anymore.
       !     May still be useful for initialization.
       !if (hyd_auxvar%pres(apid) >= hyd_auxvar% &
       !       pres(lid)*(1.d0-window_epsilon)) then
       !  if (hyd_auxvar%pres(apid) < PE_hyd) then
       if (hyd_auxvar%pres(lid) >= PE_hyd .and. &
           hyd_auxvar%xmass(acid,lid) >= xal * (1.d0 + &
           state_change_threshold)) then
         istatechng = PETSC_TRUE
         global_auxvar%istate = HAI_STATE
       elseif (hyd_auxvar%pres(lid) <= PE_hyd .and. &
               hyd_auxvar%xmass(acid,lid) >= xal * (1.d0 + &
               state_change_threshold)) then
         istatechng = PETSC_TRUE
         global_auxvar%istate = GA_STATE
       else
         if (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
           istatechng = PETSC_FALSE
         elseif (hyd_auxvar%sat(lid) > 0.d0) then
           istatechng = PETSC_TRUE
           global_auxvar%istate = L_STATE
         else
           istatechng = PETSC_FALSE
         endif
       endif
    case(HGA_STATE)
      if (hyd_auxvar%temp > Tf_ice) then

        if (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) > 0.d0 &
            .and. hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HA_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HG_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = GA_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = L_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = G_STATE
        else
          ! istatechng = PETSC_TRUE
          ! global_auxvar%istate = H_STATE
        endif

      else
        istatechng = PETSC_TRUE
        global_auxvar%istate = HGAI_STATE
      endif

    case(HAI_STATE)

      if (hyd_auxvar%pres(apid) >= PE_hyd*(1.d0-window_epsilon)) then
      !if (hyd_auxvar%pres(gid) > PE_hyd*(1.d0-window_epsilon)) then
        if (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(hid) > 0.d0 &
            .and. hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = L_STATE
          !global_auxvar%istate = AI_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(hid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HA_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HI_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = L_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0) then
          ! istatechng = PETSC_TRUE
          ! global_auxvar%istate = H_STATE
        elseif (hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
          !istatechng = PETSC_TRUE
          !global_auxvar%istate = I_STATE
        endif

      else
        istatechng = PETSC_TRUE
        global_auxvar%istate = HGAI_STATE
      endif

    case(HGI_STATE)

      if (hyd_auxvar%temp < Tf_ice) then
        if (hyd_auxvar%sat(iid) > 0.d0 .and. hyd_auxvar%sat(hid) > 0.d0 &
            .and. hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(hid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HG_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = HI_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = GI_STATE
        elseif (hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
          !istatechng = PETSC_TRUE
          !global_auxvar%istate = I_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%istate = G_STATE
        else
          ! istatechng = PETSC_TRUE
          ! global_auxvar%istate = H_STATE
        endif

      else
        istatechng = PETSC_TRUE
        global_auxvar%istate = HGAI_STATE
      endif

    case(GAI_STATE)
      !MAN: state transition doesn't naturally enter this state anymore.
      !     May be useful for initialization.
      if (hyd_auxvar%pres(apid) < PE_hyd) then !Gas phase is stable
        !if (hyd_auxvar%pres(gid) < PE_hyd) then !Gas phase is stable
        if (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) > 0.d0 &
            .and. hyd_auxvar%sat(iid) > 0.d0) then

          istatechng = PETSC_FALSE

        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) &
                  > 0.d0) then

          istatechng = PETSC_TRUE
          global_auxvar%istate = GA_STATE

        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) &
                  > 0.d0) then

          istatechng = PETSC_TRUE
          global_auxvar%istate = GI_STATE

        elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) &
                > 0.d0) then

          istatechng = PETSC_TRUE
          global_auxvar%istate = L_STATE

        elseif (hyd_auxvar%sat(gid) > 0.d0) then

          istatechng = PETSC_TRUE
          global_auxvar%istate = G_STATE

        elseif (hyd_auxvar%sat(lid) > 0.d0) then

          istatechng = PETSC_TRUE
          global_auxvar%istate = L_STATE

        endif
      else
        istatechng = PETSC_TRUE
        global_auxvar%istate = HGAI_STATE
      endif
    case(HGAI_STATE)

      if (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(gid) > 0.d0 &
          .and. hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) &
          > 0.d0) then
        istatechng = PETSC_FALSE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) &
               > 0.d0 .and. hyd_auxvar%sat(lid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = HGA_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = HAI_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = GA_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) &
              > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = HGI_STATE
      elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = L_STATE
        !global_auxvar%istate = AI_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = GA_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = HA_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = HI_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = GI_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(hid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = HG_STATE
      elseif (hyd_auxvar%sat(lid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = L_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0) then
        ! istatechng = PETSC_TRUE
        ! global_auxvar%istate = H_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%istate = G_STATE
      else
        istatechng = PETSC_FALSE
        !istatechng = PETSC_TRUE
        !global_auxvar%istate = I_STATE
      endif
  end select

  hydrate_state_changed = istatechng

  new_state = global_auxvar%istate

  !Update primary variables

  if (istatechng) then

    select case (old_state)
      case(L_STATE)
        state_change_string = 'Liquid --> '
      case(G_STATE)
        state_change_string = 'Gas --> '
      case(H_STATE)
        state_change_string = 'Hydrate --> '
      case(I_STATE)
        state_change_string = 'Ice --> '
      case(GA_STATE)
        state_change_string = 'Gas-Aqueous --> '
      case(HG_STATE)
        state_change_string = 'Hydrate-Gas --> '
      case(HA_STATE)
        state_change_string = 'Hydrate-Aqueous --> '
      case(HI_STATE)
        state_change_string = 'Hydrate-Ice --> '
      case(GI_STATE)
        state_change_string = 'Gas-Ice --> '
      case(AI_STATE)
        state_change_string = 'Aqueous-Ice --> '
      case(HGA_STATE)
        state_change_string = 'Hydrate-Gas-Aqueous --> '
      case(HAI_STATE)
        state_change_string = 'Hydrate-Aqueous-Ice --> '
      case(HGI_STATE)
        state_change_string = 'Hydrate-Gas-Ice --> '
      case(GAI_STATE)
        state_change_string = 'Gas-Aqueous-Ice --> '
      case(HGAI_STATE)
        state_change_string = '4-Phase --> '

    end select

    select case (new_state)
      case(L_STATE)
        state_change_string = trim(state_change_string) // ' Liquid'
      case(G_STATE)
        state_change_string = trim(state_change_string) // ' Gas'
      case(H_STATE)
        state_change_string = trim(state_change_string) // ' Hydrate'
      case(I_STATE)
        state_change_string = trim(state_change_string) // ' Ice'
      case(GA_STATE)
        state_change_string = trim(state_change_string) // ' Gas-Aqueous'
      case(HG_STATE)
        state_change_string = trim(state_change_string) // ' Hydrate-Gas'
      case(HA_STATE)
        state_change_string = trim(state_change_string) // ' Hydrate-Aqueous'
      case(HI_STATE)
        state_change_string = trim(state_change_string) // ' Hydrate-Ice'
      case(GI_STATE)
        state_change_string = trim(state_change_string) // ' Gas-Ice'
      case(AI_STATE)
        state_change_string = trim(state_change_string) // ' Aqueous-Ice'
      case(HGA_STATE)
        state_change_string = trim(state_change_string) // &
                              ' Hydrate-Gas-Aqueous'
      case(HAI_STATE)
        state_change_string = trim(state_change_string) // &
                              ' Hydrate-Aqueous-Ice'
      case(HGI_STATE)
        state_change_string = trim(state_change_string) // ' Hydrate-Gas-Ice'
      case(GAI_STATE)
        state_change_string = trim(state_change_string) // ' Gas-Aqueous-Ice'
      case(HGAI_STATE)
        state_change_string = trim(state_change_string) // ' 4-Phase'
    end select

    if (option%iflag == HYDRATE_UPDATE_FOR_ACCUM) then
      write(append,'('' at Cell '',i8)') natural_id
    else if (option%iflag == HYDRATE_UPDATE_FOR_DERIVATIVE) then
      write(append, &
            '(''(due to perturbation) '',i8)') natural_id
    else
      write(append, &
             '('' at Boundary Face '', i8)') natural_id
    endif

    state_change_string = trim(state_change_string) // trim(append)

    if (hydrate_restrict_state_chng) hyd_auxvar%istatechng = PETSC_TRUE

    select case(global_auxvar%istate)

      case(L_STATE)
!     ********* Aqueous State (A) ********************************
!     Primary variables: Pl, Xma, T
!

        x(HYDRATE_LIQUID_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_L_STATE_X_MOLE_DOF) = hyd_auxvar%xmass(acid,lid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(G_STATE)
!     ********* Gas State (G) ********************************
!     Primary variables: Pg, Pa, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_G_STATE_AIR_PRESSURE_DOF) = hyd_auxvar%pres(apid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(H_STATE)
!     ********* Hydrate State (H) ********************************
!     Primary variables: Pg, Xmh, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%xmol(acid,hid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(I_STATE)
!     ********* Ice State (I) ********************************
!     Primary variables: Pg, Xmi, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = 0.d0
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(GA_STATE)
!     ********* Gas & Aqueous State (GA) ********************************
!     Primary variables: Pg, Sg, T
!

        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid) 
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(HG_STATE)
!     ********* Hydrate & Gas State (HG) ********************************
!     Primary variables: Pg, Sg, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(HA_STATE)
!     ********* Hydrate & Aqueous State (HA) ********************************
!     Primary variables: Pg, Sh, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(HI_STATE)
!     ********* Hydrate & Ice State (HI) ********************************
!     Primary variables: Pg, Sh, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(GI_STATE)
!     ********* Gas & Ice State (GI) ********************************
!     Primary variables: Pg, Si, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(iid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(AI_STATE)
!     ********* Aqueous & Ice State (AI) ********************************
!     Primary variables: Pg, Xma, Sl
!
        x(HYDRATE_LIQUID_PRESSURE_DOF) = hyd_auxvar%pres(lid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%xmass(acid,lid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(lid)

      case(HGA_STATE)
!     ********* Hydrate, Gas, & Aqueous State (HGA) **************************
!     Primary variables: Sl, Sh, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%sat(lid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(HAI_STATE)
!     ********* Hydrate, Aqueous, & Ice State (HAI) **************************
!     Primary variables: Pg, Sl, Si
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(lid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(iid)

      case(HGI_STATE)
!     ********* Hydrate, Gas, & Ice State (HGI) ******************************
!     Primary variables: Sh, Si, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%sat(iid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(GAI_STATE)
!     ********* Gas, Aqueous, & Ice State (GAI) ******************************
!     Primary variables: Pg, Sg, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(HGAI_STATE)
!     ********* 4-Phase (HGAI) ********************************
!     Primary variables: Sl, Sg, Si
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%sat(lid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(iid)

      case default
        write(option%io_buffer,*) global_auxvar%istate
        option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
          ') not recognized in HydrateAuxVarUpdateState.'
        call PrintErrMsgByRank(option)

    end select

    call HydrateAuxVarCompute(x,hyd_auxvar, global_auxvar,material_auxvar, &
          characteristic_curves,hydrate_parameter,natural_id,option)

    state_change_string = 'State Transition: ' // trim(state_change_string)
    if (hydrate_print_state_transition) then
      call PrintMsgByRank(option,state_change_string)
    endif

  endif

end subroutine HydrateAuxVarUpdateState


! ************************************************************************** !

subroutine HydrateAuxVarPerturb(hyd_auxvar,global_auxvar, &
                                material_auxvar, &
                                characteristic_curves,hydrate_parameter, &
                                natural_id, option)
  !
  ! Calculates auxiliary variables for perturbed system
  !
  ! Author: Michael Nole
  ! Date: 01/30/19
  !

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(hydrate_auxvar_type) :: hyd_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  type(hydrate_parameter_type), pointer :: hydrate_parameter

  PetscReal :: x(option%nflowdof), x_pert_plus(option%nflowdof), &
               pert(option%nflowdof), x_pert_minus(option%nflowdof)

  PetscReal, parameter :: perturbation_tolerance = 1.d-8
!  PetscReal, parameter :: perturbation_tolerance = 1.d-11
  PetscReal, parameter :: min_perturbation = 1.d-10

  PetscReal, parameter :: min_mole_fraction_pert = 1.d-12
  PetscReal, parameter :: min_pres_pert = 1.d-3
  PetscReal, parameter :: min_temp_pert = 8.66d-9
  PetscReal, parameter :: min_xmass_pert = 1.d-14
  PetscReal, parameter :: min_sat_pert = 3.16d-11

  PetscReal :: xag, xwg, xal, xsl, xwl, xmolag, xmolwg, xmolal, &
               xmolsl, xmolwl, salt_mass
  PetscReal :: sigma, beta_gl
  PetscReal :: Pv, Psat, Prvap, Pa
  PetscReal :: dpl, dpg, dpa, dxa, dxs, dsg, dt
  PetscReal :: cell_pressure, sgt_max
  PetscInt :: idof

  ! Phase ID's
  PetscInt :: lid, gid, hid, pid, pwid, pbid, spid, tgid
  ! Component ID's
  PetscInt :: wid, acid, air_pressure_id, sid, pgid
  ! Other ID's
  PetscInt :: cpid, vpid, rvpid

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  pid = option%precipitate_phase
  pwid = option%pure_water_phase
  pbid = option%pure_brine_phase
  tgid = option%trapped_gas_phase
  pgid = option%trapped_gas_phase
  spid = option%saturation_pressure_id

  wid = option%water_id
  acid = option%air_id
  sid = option%salt_id

  cpid = option%capillary_pressure_id
  air_pressure_id = option%air_pressure_id
  vpid = option%vapor_pressure_id
  rvpid = option%reduced_vapor_pressure_id

  call HydrateComputeSaltSolubility(hyd_auxvar(ZERO_INTEGER)%temp,xsl)
  dxs = 1.0d-5 * xsl
  salt_mass = hyd_auxvar(ZERO_INTEGER)%m_salt(ONE_INTEGER)
  xsl = min(salt_mass,xsl)

  dt = -1.d0 * perturbation_tolerance * (hyd_auxvar(ZERO_INTEGER)%temp + &
       min_perturbation)

  call HydrateComputeSurfaceTension(hyd_auxvar(ZERO_INTEGER)%temp, xsl, sigma)
  beta_gl = CO2_REFERENCE_SURFACE_TENSION / sigma

  sgt_max = characteristic_curves%saturation_function%Sgt_max
  if(.not. Initialized(sgt_max)) then
    sgt_max = 1.d0
  endif


  select case(global_auxvar%istate)
    case(L_STATE)
      dpl = max(1.d-1, 1.d-7 * hyd_auxvar(ZERO_INTEGER)%pres(lid))
      cell_pressure = hyd_auxvar(ZERO_INTEGER)%pres(lid)
      ! cell_pressure = min(cell_pressure, 1.d8)
      call HydrateBrineSaturationPressure(hyd_auxvar(ZERO_INTEGER)%temp, xsl, &
                                       Psat)
      Prvap = Psat
      call HydrateEquilibrate(hyd_auxvar(ZERO_INTEGER)%temp,cell_pressure, &
                           global_auxvar%istate, &
                           hyd_auxvar(ZERO_INTEGER)%sat(hid), &
                           Pa, Pv, Psat, Prvap, &
                           xag, xwg, xal, xsl, xwl, &
                           xmolag, xmolwg, xmolal, xmolsl, xmolwl, &
                           characteristic_curves, material_auxvar, option)
      if (hyd_auxvar(ZERO_INTEGER)%xmass(acid,lid) > (1.d-2 * xal)) then
        dxa = sign(1.d-4 * xal, &
                     5.d-1 * xal - &
                     hyd_auxvar(ZERO_INTEGER)%xmass(acid,lid))
      else
        dxa = sign(1.d-3 * xal, &
                     5.d-1 * xal - &
                     hyd_auxvar(ZERO_INTEGER)%xmass(acid,lid))
      endif
      x(HYDRATE_LIQUID_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%liquid_phase)
      x(HYDRATE_L_STATE_X_MOLE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%xmass(option%air_id,option%liquid_phase)
      x(HYDRATE_ENERGY_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_LIQUID_PRESSURE_DOF) = dpl
      pert(HYDRATE_L_STATE_X_MOLE_DOF) = dxa
      pert(HYDRATE_ENERGY_DOF) = dt

    case(G_STATE)
      dpg = 1.d-3

      dpa = max(1.d-2, &
                  1.d-7 * dabs(hyd_auxvar(ZERO_INTEGER)%pres(gid) - &
                               hyd_auxvar(ZERO_INTEGER)%pres(lid)))

      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_G_STATE_AIR_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = dpg
      pert(HYDRATE_G_STATE_AIR_PRESSURE_DOF) = dpa
      pert(HYDRATE_ENERGY_DOF) = dt

    case(H_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar(ZERO_INTEGER)%xmol(acid,hid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = -1.d0 * &
         (perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF) + min_perturbation)
      pert(HYDRATE_GAS_SATURATION_DOF) = 999.d0 !dR should = 0
      pert(HYDRATE_ENERGY_DOF) = &
         perturbation_tolerance*x(HYDRATE_ENERGY_DOF) + min_perturbation

    case(I_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = 0.d0
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
        perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)
      pert(HYDRATE_GAS_SATURATION_DOF) = 999.d0 !dR should = 0
      pert(HYDRATE_ENERGY_DOF) = &
          perturbation_tolerance*x(HYDRATE_ENERGY_DOF)

    case(GA_STATE)
      dpl = max(1.d-2, &
                1.d-6 * dabs(hyd_auxvar(ZERO_INTEGER)%pres(gid) - &
                             hyd_auxvar(ZERO_INTEGER)%pres(lid)))
      dpl = sign(dpl, (5.d-1 - hyd_auxvar(ZERO_INTEGER)%sat(lid)))
      dpg = -1.d0 * dpl
      dsg = sign(1.d-7, 5.d-1 * sgt_max - hyd_auxvar(ZERO_INTEGER)%sat(gid))

      x(HYDRATE_GAS_PRESSURE_DOF) = &
       hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%sat(option%gas_phase)
      x(HYDRATE_ENERGY_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = dpg
      pert(HYDRATE_GAS_SATURATION_DOF) = dsg
      pert(HYDRATE_ENERGY_DOF) = dt

    case(HG_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(option%gas_phase)
      x(HYDRATE_ENERGY_DOF) = &
           hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)+min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
           perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(HA_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(hid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
        perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)+min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
         perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(HI_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(hid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)+min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
           perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(GI_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(iid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)+min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
           perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(AI_STATE)
      x(HYDRATE_LIQUID_PRESSURE_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%pres(option%liquid_phase)
      x(HYDRATE_L_STATE_X_MOLE_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%xmass(option%air_id,option%liquid_phase)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(lid)
      pert(HYDRATE_LIQUID_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_LIQUID_PRESSURE_DOF) + &
         min_perturbation

      if (x(HYDRATE_L_STATE_X_MOLE_DOF) > &
           1.d3 * perturbation_tolerance) then
        pert(HYDRATE_L_STATE_X_MOLE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_L_STATE_X_MOLE_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_ENERGY_DOF) > 0.5d0) then
        pert(HYDRATE_ENERGY_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_ENERGY_DOF) = perturbation_tolerance
      endif

    case(HGA_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%sat(lid)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(hid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_ENERGY_DOF) = dt
      if (x(HYDRATE_GAS_SATURATION_DOF) + x(HYDRATE_GAS_PRESSURE_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
        pert(HYDRATE_GAS_PRESSURE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
        pert(HYDRATE_GAS_PRESSURE_DOF) = perturbation_tolerance
      endif

    case(HAI_STATE)

      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(lid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(iid)

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF) + &
         min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_ENERGY_DOF) > 0.5d0) then
        pert(HYDRATE_ENERGY_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_ENERGY_DOF) = perturbation_tolerance
      endif

    case(HGI_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%sat(hid)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(iid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      if (x(HYDRATE_GAS_PRESSURE_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_PRESSURE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_PRESSURE_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
           perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(GAI_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(gid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF) + &
         min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
           perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(HGAI_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%sat(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(lid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(iid)

      if (x(HYDRATE_GAS_PRESSURE_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_PRESSURE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_PRESSURE_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_ENERGY_DOF) > 0.5d0) then
        pert(HYDRATE_ENERGY_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_ENERGY_DOF) = perturbation_tolerance
      endif

  end select
  ! HYDRATE_UPDATE_FOR_DERIVATIVE indicates call from perturbation

  option%iflag = HYDRATE_UPDATE_FOR_DERIVATIVE

  do idof = 1, option%nflowdof

    if (hydrate_central_diff_jacobian) then
      ! pert(idof) = max(1.d-7 * x(idof),1.d-7)

      x_pert_minus = x
      x_pert_minus(idof) = x(idof) - pert(idof)
      call HydrateAuxVarCompute(x_pert_minus, &
             hyd_auxvar(idof+option%nflowdof),global_auxvar,material_auxvar, &
             characteristic_curves,hydrate_parameter,natural_id,option)

    endif

    hyd_auxvar(idof)%pert = pert(idof)
    hyd_auxvar(idof+option%nflowdof)%pert = pert(idof)
    x_pert_plus = x
    x_pert_plus(idof) = x(idof) + pert(idof)
    call HydrateAuxVarCompute(x_pert_plus,hyd_auxvar(idof),global_auxvar, &
                              material_auxvar, characteristic_curves, &
                              hydrate_parameter,natural_id,option)
  enddo

end subroutine HydrateAuxVarPerturb


! ************************************************************************** !

subroutine HydratePrintAuxVars(hydrate_auxvar,global_auxvar,material_auxvar, &
                               natural_id,string,option)
  !
  ! Prints out the contents of an auxvar
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: hydrate_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_energy, gas_energy
  PetscReal :: liquid_saturation, gas_saturation

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_energy = 0.d0
  gas_energy = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  select case(global_auxvar%istate)
    case(L_STATE)
      print *, '     Thermodynamic state: Liquid phase'
      liquid_density = hydrate_auxvar%den(lid)
      liquid_energy = hydrate_auxvar%U(lid)
      liquid_saturation = hydrate_auxvar%sat(lid)
    case(G_STATE)
      print *, '     Thermodynamic state: Gas phase'
      gas_density = hydrate_auxvar%den(gid)
      gas_energy = hydrate_auxvar%U(gid)
      gas_saturation = hydrate_auxvar%sat(gid)
    case(GA_STATE)
      print *, '     Thermodynamic state: Two phase'
      liquid_density = hydrate_auxvar%den(lid)
      gas_density = hydrate_auxvar%den(gid)
      liquid_energy = hydrate_auxvar%U(lid)
      gas_energy = hydrate_auxvar%U(gid)
      liquid_saturation = hydrate_auxvar%sat(lid)
      gas_saturation = hydrate_auxvar%sat(gid)
  end select
  liquid_mass = (liquid_density*hydrate_auxvar%xmol(lid,lid)* &
                 liquid_saturation+ &
                 gas_density*hydrate_auxvar%xmol(lid,gid)* &
                 gas_saturation)* &
                 hydrate_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*hydrate_auxvar%xmol(gid,lid)* &
              liquid_saturation+ &
              gas_density*hydrate_auxvar%xmol(gid,gid)* &
              gas_saturation)* &
              hydrate_auxvar%effective_porosity*material_auxvar%volume
  print *, 'tot liq comp mass [kmol]: ', liquid_mass
  print *, 'tot gas comp mass [kmol]: ', gas_mass
  print *, '             energy [MJ]: ', liquid_mass*liquid_energy + &
                                         gas_mass*gas_energy
  print *, '         liquid pressure: ', hydrate_auxvar%pres(lid)
  print *, '            gas pressure: ', hydrate_auxvar%pres(gid)
  print *, '            air pressure: ', hydrate_auxvar%pres(apid)
  print *, '      capillary pressure: ', hydrate_auxvar%pres(cpid)
  print *, '          vapor pressure: ', hydrate_auxvar%pres(vpid)
  print *, '     saturation pressure: ', hydrate_auxvar%pres(spid)
  print *, '       liquid saturation: ', hydrate_auxvar%sat(lid)
  print *, '          gas saturation: ', hydrate_auxvar%sat(gid)
  print *, '   liquid density [kmol]: ', hydrate_auxvar%den(lid)
  print *, '      gas density [kmol]: ', hydrate_auxvar%den(gid)
  print *, '     liquid density [kg]: ', hydrate_auxvar%den_kg(lid)
  print *, '        gas density [kg]: ', hydrate_auxvar%den_kg(gid)
  print *, '         temperature [C]: ', hydrate_auxvar%temp
  print *, '      liquid H [MJ/kmol]: ', hydrate_auxvar%H(lid)
  print *, '         gas H [MJ/kmol]: ', hydrate_auxvar%H(gid)
  print *, '      liquid U [MJ/kmol]: ', hydrate_auxvar%U(lid)
  print *, '         gas U [MJ/kmol]: ', hydrate_auxvar%U(gid)
  print *, '     X (water in liquid): ', hydrate_auxvar%xmol(lid,lid)
  print *, '       X (air in liquid): ', hydrate_auxvar%xmol(gid,lid)
  print *, '        X (water in gas): ', hydrate_auxvar%xmol(lid,gid)
  print *, '          X (air in gas): ', hydrate_auxvar%xmol(gid,gid)
  print *, '         liquid mobility: ', hydrate_auxvar%mobility(lid)
  print *, '            gas mobility: ', hydrate_auxvar%mobility(gid)
  print *, '      effective porosity: ', hydrate_auxvar%effective_porosity
  print *, '--------------------------------------------------------'

end subroutine HydratePrintAuxVars

! ************************************************************************** !

subroutine HydrateOutputAuxVars1(hydrate_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,string,append,option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: hydrate_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  PetscBool :: append
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string2
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_energy, gas_energy
  PetscReal :: liquid_saturation, gas_saturation

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_energy = 0.d0
  gas_energy = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0

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
    case(L_STATE)
      write(86,*) ' Thermodynamic state: Liquid phase'
      liquid_density = hydrate_auxvar%den(lid)
      liquid_energy = hydrate_auxvar%U(lid)
      liquid_saturation = hydrate_auxvar%sat(lid)
    case(G_STATE)
      write(86,*) ' Thermodynamic state: Gas phase'
      gas_density = hydrate_auxvar%den(gid)
      gas_energy = hydrate_auxvar%U(gid)
      gas_saturation = hydrate_auxvar%sat(gid)
    case(GA_STATE)
      write(86,*) ' Thermodynamic state: Two phase'
      liquid_density = hydrate_auxvar%den(lid)
      gas_density = hydrate_auxvar%den(gid)
      liquid_energy = hydrate_auxvar%U(lid)
      gas_energy = hydrate_auxvar%U(gid)
      liquid_saturation = hydrate_auxvar%sat(lid)
      gas_saturation = hydrate_auxvar%sat(gid)
  end select
  liquid_mass = (liquid_density*hydrate_auxvar%xmol(lid,lid)* &
                 liquid_saturation+ &
                 gas_density*hydrate_auxvar%xmol(lid,gid)* &
                 gas_saturation)* &
                 hydrate_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*hydrate_auxvar%xmol(gid,lid)* &
              liquid_saturation+ &
              gas_density*hydrate_auxvar%xmol(gid,gid)* &
              gas_saturation)* &
              hydrate_auxvar%effective_porosity*material_auxvar%volume
  write(86,*) 'tot liq comp mass [kmol]: ', liquid_mass
  write(86,*) 'tot gas comp mass [kmol]: ', gas_mass
  write(86,*) '             energy [MJ]: ', liquid_mass*liquid_energy + &
                                            gas_mass*gas_energy
  write(86,*) '         liquid pressure: ', hydrate_auxvar%pres(lid)
  write(86,*) '            gas pressure: ', hydrate_auxvar%pres(gid)
  write(86,*) '            air pressure: ', hydrate_auxvar%pres(apid)
  write(86,*) '      capillary pressure: ', hydrate_auxvar%pres(cpid)
  write(86,*) '          vapor pressure: ', hydrate_auxvar%pres(vpid)
  write(86,*) '     saturation pressure: ', hydrate_auxvar%pres(spid)
  write(86,*) '         temperature [C]: ', hydrate_auxvar%temp
  write(86,*) '       liquid saturation: ', hydrate_auxvar%sat(lid)
  write(86,*) '          gas saturation: ', hydrate_auxvar%sat(gid)
  write(86,*) '   liquid density [kmol]: ', hydrate_auxvar%den(lid)
  write(86,*) '     liquid density [kg]: ', hydrate_auxvar%den_kg(lid)
  write(86,*) '      gas density [kmol]: ', hydrate_auxvar%den(gid)
  write(86,*) '        gas density [kg]: ', hydrate_auxvar%den_kg(gid)
  write(86,*) '     X (water in liquid): ', hydrate_auxvar%xmol(lid,lid)
  write(86,*) '       X (air in liquid): ', hydrate_auxvar%xmol(gid,lid)
  write(86,*) '        X (water in gas): ', hydrate_auxvar%xmol(lid,gid)
  write(86,*) '          X (air in gas): ', hydrate_auxvar%xmol(gid,gid)
  write(86,*) '      liquid H [MJ/kmol]: ', hydrate_auxvar%H(lid)
  write(86,*) '         gas H [MJ/kmol]: ', hydrate_auxvar%H(gid)
  write(86,*) '      liquid U [MJ/kmol]: ', hydrate_auxvar%U(lid)
  write(86,*) '         gas U [MJ/kmol]: ', hydrate_auxvar%U(gid)
  write(86,*) '         liquid mobility: ', hydrate_auxvar%mobility(lid)
  write(86,*) '            gas mobility: ', hydrate_auxvar%mobility(gid)
  write(86,*) '      effective porosity: ', hydrate_auxvar%effective_porosity
  write(86,*) '...'
  write(86,*) liquid_mass
  write(86,*) gas_mass
  write(86,*) liquid_mass*hydrate_auxvar%U(lid) + &
              gas_mass*hydrate_auxvar%U(gid)
  write(86,*) hydrate_auxvar%pres(lid)
  write(86,*) hydrate_auxvar%pres(gid)
  write(86,*) hydrate_auxvar%pres(apid)
  write(86,*) hydrate_auxvar%pres(cpid)
  write(86,*) hydrate_auxvar%pres(vpid)
  write(86,*) hydrate_auxvar%pres(spid)
  write(86,*) hydrate_auxvar%temp
  write(86,*) hydrate_auxvar%sat(lid)
  write(86,*) hydrate_auxvar%sat(gid)
  write(86,*) hydrate_auxvar%den(lid)
  write(86,*) hydrate_auxvar%den_kg(lid)
  write(86,*) hydrate_auxvar%den(gid)
  write(86,*) hydrate_auxvar%den_kg(gid)
  write(86,*) hydrate_auxvar%xmol(lid,lid)
  write(86,*) hydrate_auxvar%xmol(gid,lid)
  write(86,*) hydrate_auxvar%xmol(lid,gid)
  write(86,*) hydrate_auxvar%xmol(gid,gid)
  write(86,*) hydrate_auxvar%H(lid)
  write(86,*) hydrate_auxvar%H(gid)
  write(86,*) hydrate_auxvar%U(lid)
  write(86,*) hydrate_auxvar%U(gid)
  write(86,*) ''
  write(86,*) hydrate_auxvar%mobility(lid)
  write(86,*) hydrate_auxvar%mobility(gid)
  write(86,*) hydrate_auxvar%effective_porosity
  write(86,*) '--------------------------------------------------------'

  close(86)

end subroutine HydrateOutputAuxVars1

! ************************************************************************** !

subroutine HydrateOutputAuxVars2(hydrate_auxvars,global_auxvars,option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Global_Aux_module
  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: hydrate_auxvars(0:,:)
  type(global_auxvar_type) :: global_auxvars(:)
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: apid, cpid, vpid
  PetscInt :: gid, lid, acid, wid, eid
  PetscInt :: i, n, idof

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  string = 'hydrate_auxvar.txt'
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
    ((hydrate_auxvars(idof,i)%pres(lid),i=1,n),idof=0,3)
  write(86,100) '         gas pressure: ', &
    ((hydrate_auxvars(idof,i)%pres(gid),i=1,n),idof=0,3)
  write(86,100) '         air pressure: ', &
    ((hydrate_auxvars(idof,i)%pres(apid),i=1,n),idof=0,3)
  write(86,100) '   capillary pressure: ', &
    ((hydrate_auxvars(idof,i)%pres(cpid),i=1,n),idof=0,3)
  write(86,100) '       vapor pressure: ', &
    ((hydrate_auxvars(idof,i)%pres(vpid),i=1,n),idof=0,3)
  write(86,100) '      temperature [C]: ', &
    ((hydrate_auxvars(idof,i)%temp,i=1,n),idof=0,3)
  write(86,100) '    liquid saturation: ', &
    ((hydrate_auxvars(idof,i)%sat(lid),i=1,n),idof=0,3)
  write(86,100) '       gas saturation: ', &
    ((hydrate_auxvars(idof,i)%sat(gid),i=1,n),idof=0,3)
  write(86,100) 'liquid density [kmol]: ', &
    ((hydrate_auxvars(idof,i)%den(lid),i=1,n),idof=0,3)
  write(86,100) '  liquid density [kg]: ', &
    ((hydrate_auxvars(idof,i)%den_kg(lid),i=1,n),idof=0,3)
  write(86,100) '   gas density [kmol]: ', &
    ((hydrate_auxvars(idof,i)%den(gid),i=1,n),idof=0,3)
  write(86,100) '     gas density [kg]: ', &
    ((hydrate_auxvars(idof,i)%den_kg(gid),i=1,n),idof=0,3)
  write(86,100) '  X (water in liquid): ', &
    ((hydrate_auxvars(idof,i)%xmol(lid,lid),i=1,n),idof=0,3)
  write(86,100) '    X (air in liquid): ', &
    ((hydrate_auxvars(idof,i)%xmol(gid,lid),i=1,n),idof=0,3)
  write(86,100) '     X (water in gas): ', &
    ((hydrate_auxvars(idof,i)%xmol(lid,gid),i=1,n),idof=0,3)
  write(86,100) '       X (air in gas): ', &
    ((hydrate_auxvars(idof,i)%xmol(gid,gid),i=1,n),idof=0,3)
  write(86,100) '   liquid H [MJ/kmol]: ', &
    ((hydrate_auxvars(idof,i)%H(lid),i=1,n),idof=0,3)
  write(86,100) '      gas H [MJ/kmol]: ', &
    ((hydrate_auxvars(idof,i)%H(gid),i=1,n),idof=0,3)
  write(86,100) '   liquid U [MJ/kmol]: ', &
    ((hydrate_auxvars(idof,i)%U(lid),i=1,n),idof=0,3)
  write(86,100) '      gas U [MJ/kmol]: ', &
    ((hydrate_auxvars(idof,i)%U(gid),i=1,n),idof=0,3)
  write(86,*)
  write(86,100) '      liquid mobility: ', &
    ((hydrate_auxvars(idof,i)%mobility(lid),i=1,n),idof=0,3)
  write(86,100) '         gas mobility: ', &
    ((hydrate_auxvars(idof,i)%mobility(gid),i=1,n),idof=0,3)
  write(86,100) '   effective porosity: ', &
    ((hydrate_auxvars(idof,i)%effective_porosity,i=1,n),idof=0,3)

  close(86)

end subroutine HydrateOutputAuxVars2

! ************************************************************************** !

subroutine HydrateAuxVarSingleDestroy(auxvar)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  implicit none

  type(hydrate_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call HydrateAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine HydrateAuxVarSingleDestroy

! ************************************************************************** !

subroutine HydrateAuxVarArray1Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  implicit none

  type(hydrate_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call HydrateAuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine HydrateAuxVarArray1Destroy

! ************************************************************************** !

subroutine HydrateAuxVarArray2Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  implicit none

  type(hydrate_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call HydrateAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine HydrateAuxVarArray2Destroy

! ************************************************************************** !

subroutine HydrateAuxVarStrip(auxvar)
  !
  ! HydrateAuxVarDestroy: Deallocates a hydrate auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(hydrate_auxvar_type) :: auxvar

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

end subroutine HydrateAuxVarStrip

! ************************************************************************** !

subroutine HydrateAuxDestroy(aux)
  !
  ! Deallocates a hydrate auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(hydrate_type), pointer :: aux

  if (.not.associated(aux)) return

  call HydrateAuxVarDestroy(aux%auxvars)
  call HydrateAuxVarDestroy(aux%auxvars_bc)
  call HydrateAuxVarDestroy(aux%auxvars_ss)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%hydrate_parameter)) then
    call DeallocateArray(aux%hydrate_parameter%diffusion_coefficient)
    deallocate(aux%hydrate_parameter)
  endif
  nullify(aux%hydrate_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine HydrateAuxDestroy

! ************************************************************************** !

subroutine HydrateCompositeThermalCond(phi,sat,kdry,kwet,keff)

  implicit none

  PetscReal :: phi, kdry, kwet
  PetscReal, pointer :: sat(:)

  PetscReal :: keff
  PetscReal :: k_h2o,k_ch4,k_hyd,k_ice
  PetscInt :: lid, gid, hid, iid

  lid = 1
  gid = 2
  hid = 3
  iid = 4

  k_h2o = 0.59d0 !W/m-K
  k_ch4 = 30.d-3 !W/m-K
  k_hyd = 0.58d0 !W/m-K
  k_ice = 2.2d0   !W/m-K


  select case(hydrate_tcond)
    case(0)
      keff = sqrt(sat(lid)) * (kwet - kdry)
    case(1)
      ! IGHCC2 function (seems odd if phi = 1 and sat(lid) = 1)
      keff = kdry + phi * (sat(lid)*kwet + sat(hid)*k_hyd + sat(iid) * k_ice &
             + sat(gid)*k_ch4)
    case(2) ! Default function
      keff = kdry + phi * (sat(lid)*k_h2o + sat(hid)*k_hyd + &
             sat(iid)*k_ice + sat(gid)*k_ch4)
    case(3) ! Grenier et al. benchmark
      keff = phi * sat(lid) * 0.6 + phi * sat(iid) * 2.14 + (1-phi) * 9.d0
  end select


end subroutine HydrateCompositeThermalCond

! ************************************************************************** !

subroutine HydratePE(T, sat, PE, dP, characteristic_curves, material_auxvar, &
                     option)
  !
  ! This subroutine calculates the 3-phase equilibrium pressure of methane
  ! hydrate in pure water, from polynomial fit (Moridis, 2003)
  !
  ! Author: Michael Nole
  ! Date: 01/22/19
  !

  use Characteristic_Curves_module
  use Option_module
  use Material_Aux_module

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: sat
  PetscReal, intent(out) :: PE
  PetscReal, intent(out) :: dP

  class(characteristic_curves_type) :: characteristic_curves
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option

  PetscReal :: T_temp, dTf, Tf_ice
  PetscReal :: a, b, c, T_k

  T_temp = T + T273K
  T_k = T_temp
  dP = 0.d0

  select case(hydrate_former)
    case(HYDRATE_FORMER_NULL)
      ! Hydrate does not form
      PE = 1.d20 
    case(HYDRATE_FORMER_CH4)
      if (hydrate_with_gibbs_thomson) then
        call GibbsThomsonHydrate(1.d0-sat,L_GH, CH4_HYDRATE_DENSITY, T, dTf, &
                                  characteristic_curves, material_auxvar, option)
      else
        dTf = 0.d0
      endif
      Tf_ice = dTf
      if (T < Tf_ice) then
        select case(hydrate_phase_boundary)
          case(1)
            !Kamath, 1984
            PE = exp(1.4717d1-1.88679d3/T_temp)*1.d-3
            if (hydrate_adjust_ghsz_solubility) then
              dP = PE - exp(1.4717d1-1.88679d3/T_temp)*1.d-3
              dP = dP * 1.d6
            endif
          case(2)
            !Moridis, 2003
            PE = exp(-43.8921173434628 + 0.776302133739303 * T_temp &
                 - 7.27291427030502d-3 * T_temp**2 + 3.85413985900724d-5 * &
                 T_temp**3 - 1.03669656828834d-7 * T_temp**4 + &
                 1.09882180475307d-10 * T_temp**5)
            if (hydrate_adjust_ghsz_solubility) then
              dP = PE - exp(-43.8921173434628 + 0.776302133739303 * &
                   (T_temp-dTf) - 7.27291427030502d-3 * (T_temp-dTf)**2 + &
                   3.85413985900724d-5 * T_temp**3 - 1.03669656828834d-7 * &
                   (T_temp-dTf)**4 + 1.09882180475307d-10 * (T_temp-dTf)**5)
              dP = dP * 1.d6
            endif
          case(3)
            !Moridis, 2003 simple
            PE = exp(0.0334940999*T_temp - 8.1938174346)
            if (hydrate_adjust_ghsz_solubility) then
              dP = PE - exp(0.0334940999*(T_temp-dTf) - 8.1938174346)
              dP = dP * 1.d6
            endif
        end select
      else
        select case(hydrate_phase_boundary)
          case(1)
            !Kamath, 1984
            PE = exp(3.898d1-8.533d3/T_temp)*1.d-3
            if (hydrate_adjust_ghsz_solubility) then
              dP = PE - exp(3.898d1 - 8.533d3/(T_temp - dTf))* 1.d-3
              dP = dP * 1.d6
            endif
          case(2)
            !Moridis, 2003
            PE = exp(-1.9413850446456d5 + 3.31018213397926d3 * T_temp &
                 - 22.5540264493806* T_temp**2 + 0.0767559117787059 * &
                 T_temp**3 - 1.30465829788791d-4 * T_temp**4 + &
                 8.86065316687571d-8 * T_temp**5)
            if (hydrate_adjust_ghsz_solubility) then
              dP = PE - exp(-1.9413850446456d5 + 3.31018213397926d3 * &
                   (T_temp-dTf) - 22.5540264493806*(T_temp-dTf)**2 + &
                   0.0767559117787059 * (T_temp-dTf)**3 - &
                   1.30465829788791d-4 * (T_temp-dTf)**4 + &
                   8.86065316687571d-8 * (T_temp-dTf)**5)
              dP = dP * 1.d6
            endif
          case(3)
            !Moridis, 2003 simple
            PE = exp(0.1100383278*T_temp - 29.1133440975)
            if (hydrate_adjust_ghsz_solubility) then
              dP = PE - exp(0.1100383278*(T_temp-dTf) - 29.1133440975)
              dP = dP * 1.d6
            endif
        end select
      endif
    case(HYDRATE_FORMER_CO2)
      if (hydrate_with_gibbs_thomson) then
        call GibbsThomsonHydrate(1.d0-sat,L_GH, CO2_HYDRATE_DENSITY, T, dTf, &
                                  characteristic_curves, material_auxvar, option)
      else
        dTf = 0.d0
      endif
      Tf_ice = dTf
      !Sloan compilation fit (Clathrate Hydrates of Natural Gases)
    !   if (T < TQD) then
    !     PE = 1.1046 + 0.04449 * (T_temp - 273.15) + 0.000629 * &
    !          (T_temp - 273.15) ** 2
    !   else
    !     PE = 1.2241 + 0.13700 * (T_temp - 273.15) ** 2 - 0.0015018 * (T_temp - &
    !           273.15) ** 3 + 0.0001733 * (T_temp - 273.15) ** 4
    !  endif
      if (T_k < 275.4d0) then
        a = 1.29916148d0
        b = 2.539998564d1
        c = 2.7033007553d2
        PE = (sqrt(T_k+(b**2/(4.d0*a))-c) - b/(2*sqrt(a))) / sqrt(a)
        dP = exp(PE) - &
             exp(((sqrt((T_k-dTf)+(b**2/(4.d0*a))-c) - b/(2*sqrt(a))) / sqrt(a)))
      else
        a = 9.13936d-3
        b = -4.86611852
        c = 6.4721366487d2
        PE = a * T_k**2 + b * T_k + c
        dP = exp(PE) - exp((a * (T_k-dTf)**2 + b * (T_k-dTf) + c))
     endif
     PE = exp(PE)
     dP = dP * 1.d6
  end select

  PE = PE * 1.d6

end subroutine HydratePE

! ************************************************************************** !
subroutine HydrateMethanogenesis(z,offset,hydrate_parameter,q_meth)

  ! A simple methanogenesis source parameterized as a function of depth
  ! assuming top of domain is the seafloor.
  ! Author: Michael Nole
  ! Date: 03/05/19
  !

  implicit none

  PetscReal :: z, offset
  type(hydrate_parameter_type), pointer :: hydrate_parameter
  PetscReal :: q_meth

  type(methanogenesis_type), pointer :: methanogenesis
  PetscReal :: alpha, k_alpha, lambda, omega, z_smt

  methanogenesis => hydrate_parameter%methanogenesis

  alpha = methanogenesis%alpha
  k_alpha = methanogenesis%k_alpha
  lambda = methanogenesis%lambda
  omega = methanogenesis%omega
  z_smt = methanogenesis%z_smt

  if (offset - z > z_smt) then
    ! Malinverno, 2011
    q_meth = k_alpha * lambda * alpha * exp(-lambda/omega * (offset - &
                                    z - z_smt))
  else
    q_meth = 0.d0
  endif

  !kg/m^3/s to kmol/s
  q_meth = q_meth / MW_CH4

end subroutine HydrateMethanogenesis
! ************************************************************************** !

subroutine HydrateComputeEffectiveSat(hyd_auxvar,g_sat_eff,&
                                        h_sat_eff,i_sat_eff)
  !
  ! Computes effective saturation assuming equal nonwetting phase
  ! distribution in large pores.
  !
  ! Author: Michael Nole
  ! Date: 09/01/2022
  !

  implicit none

  type(hydrate_auxvar_type) :: hyd_auxvar
  PetscReal :: g_sat_eff,h_sat_eff,i_sat_eff

  if (hydrate_gt_3phase) then
    if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
        h_sat_eff = 1.d0 - hyd_auxvar%sat(lid)
        if (hyd_auxvar%sat(iid) > hyd_auxvar%sat(gid)) then
          ! Sh > Si > Sg
          i_sat_eff = 2.d0 * hyd_auxvar%sat(iid) + hyd_auxvar%sat(gid)
          g_sat_eff = 3.d0 * hyd_auxvar%sat(gid)
        else
          ! Sh > Sg > Si
          g_sat_eff = 2.d0 * hyd_auxvar%sat(gid) + hyd_auxvar%sat(iid)
          i_sat_eff = 3.d0 * hyd_auxvar%sat(iid)
        endif
      else
        ! Sg > Sh > Si
        g_sat_eff = 1.d0 - hyd_auxvar%sat(lid)
        h_sat_eff = 2.d0 * hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid)
        i_sat_eff = 3.d0 * hyd_auxvar%sat(iid)
      endif
    elseif (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
      ! Si > Sh > Sg
      i_sat_eff = 1.d0 - hyd_auxvar%sat(lid)
      h_sat_eff = 2.d0 * hyd_auxvar%sat(hid) + hyd_auxvar%sat(gid)
      g_sat_eff = 3.d0 * hyd_auxvar%sat(gid)
    elseif (hyd_auxvar%sat(iid) > hyd_auxvar%sat(gid)) then
      ! Si > Sg > Sh
      i_sat_eff = 1.d0 - hyd_auxvar%sat(lid)
      g_sat_eff = 2.d0 * hyd_auxvar%sat(gid) + hyd_auxvar%sat(hid)
      h_sat_eff = 3.d0 * hyd_auxvar%sat(hid)
    else
      ! Sg > Si > Sh
      g_sat_eff = 1.d0 - hyd_auxvar%sat(lid)
      i_sat_eff = 2.d0 * hyd_auxvar%sat(iid) + hyd_auxvar%sat(hid)
      h_sat_eff = 3.d0 * hyd_auxvar%sat(hid)
    endif
  else
    g_sat_eff = hyd_auxvar%sat(gid)
    h_sat_eff = hyd_auxvar%sat(hid)
    i_sat_eff = hyd_auxvar%sat(iid)
  endif

end subroutine HydrateComputeEffectiveSat

! ************************************************************************** !

subroutine HydrateGHSZSolubilityCorrection(T,P,dP,K_H)

  !Adjusts methane solubility within the hydrate stabilty zone, following
  !Davie et al., 2004
  !
  !Author: Michael Nole
  !

  implicit none

  PetscReal, intent(in) :: T, P, dP
  PetscReal, intent(inout) :: K_H

  PetscReal, parameter :: alpha = 14.4d0 !C
  PetscReal :: logP, P_MPa
  PetscReal :: delta_pressure
  PetscReal :: T3, T_k
  PetscReal :: a, b, c, P_ln  

  T_k = T + T273K
  P_MPa = P * 1.d6
  delta_pressure = P-dP

  if (dP < 0.d0 .or. delta_pressure < 0.d0) return

  P_ln = log(delta_pressure)

  ! Inverting the phase boundary
  select case(hydrate_former)
    case(HYDRATE_FORMER_NULL)
      ! Hydrate does not form
      T3 = -999.d0
    case(HYDRATE_FORMER_CH4) 
      if (P_MPa > 2.4638d0) then
        select case (hydrate_phase_boundary)
          case(1)
            !Kamath
            T3 = -8.533d3/(log((delta_pressure)*1.d-6*1.d3)-3.898d1)
          case(2)
            !Moridis
            !Lower-order
            T3 = 9.0622d0 * log((delta_pressure)*1.d-6) + 264.66d0

            !Higher-order
            !logP = log((P-dP)*1.d-6)
            !T3 = -0.0109874018d0*logP**6 + 0.17330155005d0*logP**5 &
            !     - 0.9678974011d0*logP**4 + 2.3491936188d0*logP**3 &
            !     - 2.7714662486d0*logP**2 + 11.3889445128d0*logP + 263.4959590135d0

          case(3)
            !Moridis, 2003 simple
            logP = log((delta_pressure)*1.d-6)
            T3 = (logP + 29.1133440975)/0.1100383278
          end select
      else
        select case (hydrate_phase_boundary)
          case(1)
            T3 = -1.88679d3/(log((delta_pressure)*1.d-6*1.d3) - 1.4717d1)
          case(3)
            T3 = (log((delta_pressure)*1.d-6) + 8.1938174346)/0.0334940999
        end select
        T3 = T + T273K
      endif
    case(HYDRATE_FORMER_CO2)
      if (P_ln < 0.d0) then
        a = 1.29916148d0
        b = 2.539998564d1
        c = 2.7033007553d2
        T3 = a * P_ln**2 + b * P_ln + c
      else
        a = 9.13936d-3
        b = -4.86611852d0
        c = 6.4721366487d2
        T3 = (sqrt(P_ln+(b**2/(4.d0*a))-c) - b/(2*sqrt(a))) / sqrt(a)
      endif
  end select

  if (T_k < T3) K_H = K_H / exp((T_k-T3)/alpha)

end subroutine HydrateGHSZSolubilityCorrection

! ************************************************************************** !

subroutine CalcFreezingTempDepression(sat,Tf_ice,characteristic_curves,dTf,option)

  !This subroutine ties the capillary pressure function to a
  !subcooling required to form ice in pores.
  !
  !Author: Michael Nole
  !Date: 12/05/22
  !

  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Option_module

  implicit none

  PetscReal, intent(in) :: sat !liquid saturation
  PetscReal, intent(in) :: Tf_ice
  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal, intent(out) :: dTf

  PetscReal, parameter :: gravity = EARTH_GRAVITY
  PetscReal :: Pc,dw,dpc_dsatl
  PetscReal :: sigma, theta, beta

  sigma = 0.073d0 !interfacial tension
  theta = 0.d0 !wetting angle
  dw = ICE_DENSITY !density of water
  beta = 1.d0

  !Clausius-Clapeyron derivation
  call characteristic_curves%saturation_function% &
         CapillaryPressure(sat,Pc,dpc_dsatl,option)
  select type(sf => characteristic_curves%saturation_function)
    class is (sat_func_VG_STOMP_type)
      ! Pc is the capillary head
      Pc = Pc * LIQUID_REFERENCE_DENSITY * gravity / beta
    class default
   end select
  dTf = Pc/(L_ICE * dw * 1.d6) * (Tf_ice + T273K)

end subroutine CalcFreezingTempDepression

! ************************************************************************** !

subroutine GibbsThomsonHydrate(sat,Hf,rho,Tb,dTf,characteristic_curves,&
                                material_auxvar,option)

  !This subroutine ties the capillary pressure function to a Gibbs-Thomson
  !subcooling required to precipitate a solid in pores.
  !
  !Author: Michael Nole
  !Date: 04/04/19
  !

  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Option_module
  use Material_Aux_module

  implicit none

  PetscReal, intent(in) :: sat
  PetscReal, intent(in) :: Hf
  PetscReal, intent(in) :: rho
  PetscReal, intent(in) :: Tb
  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  type(material_auxvar_type) :: material_auxvar
  PetscReal, intent(out) :: dTf

  PetscReal, parameter :: gravity = EARTH_GRAVITY
  PetscReal :: Pc,sat_temp,dpc_dsatl,sigma,theta,beta

  sigma = 0.073d0
  theta = 0.d0
  beta = 1.d0

  sat_temp = sat !- hydrate_phase_chng_epsilon !accounting for buffer

  call characteristic_curves%saturation_function% &
            CapillaryPressure(sat_temp,Pc,dpc_dsatl,option)
  select type(sf => characteristic_curves%saturation_function)
    class is (sat_func_VG_STOMP_type)
      ! Pc is the capillary head
      Pc = Pc * LIQUID_REFERENCE_DENSITY * gravity / beta
    class default
  end select
  dTf = (Tb+T273K)*Pc/(Hf * rho * 1000.d0)

end subroutine GibbsThomsonHydrate

! ************************************************************************** !

subroutine EOSHydrateEnthalpy(T,H)

  !Enthalpy of gas hydrate as f(Temperature) (Handa, 1998)
  !
  !Author: Michael Nole, David Fukuyama
  !Date: 01/22/19 
  !
  implicit none

  PetscReal, intent(in):: T
  PetscReal, intent(out) :: H

  PetscReal, parameter :: Hh0 = -54734.d0 ! J/mol
  PetscReal, parameter :: Hh0_c = -66800.d0  !J/mol
  PetscReal :: Cph, T_temp

  T_temp = T !TQD in C

  ! Integral of Cph * dT ; Cph from Handa, 1998

  ! Units: J/mol
  !H = Hh0 + 6.6d0 * (T_temp-T273K) + 7.269d-1 * (T_temp**2 - T273K**2) - 1.21333d-3 * &
  !      (T_temp**3 - T273K**3)  + 1.578d-6 * (T_temp**4 - T273K**4)
  ! Units: MJ/kmol
  !H = H / 1.d3

  !H = H / (Nhyd+1.d0)

  select case(hydrate_former)
    case(HYDRATE_FORMER_CH4)
      !Constant Cp
      Cph = 1620.d0*(FMWH2O*CH4_HYDRATION_NUMBER + FMWCH4)/1.d3
      H = Cph * (T-TQD) + Hh0 / (CH4_HYDRATION_NUMBER + 1.d0)
      H = H / 1.d3
    case(HYDRATE_FORMER_CO2)
      Cph = 1620.d0*(FMWH2O*CO2_HYDRATION_NUMBER + FMWCO2)/1.d3
      H = Cph * (T-TQD) + Hh0_c / (CO2_HYDRATION_NUMBER + 1.d0)
      H = H / 1.d3
    case default
      H = 0.d0
  end select

end subroutine EOSHydrateEnthalpy

! ************************************************************************** !

subroutine HydrateSalinityOffset(xmol,dTd)

  !
  ! This ties salinity to the subcooling required to precipitate
  ! hydrate in pores, similar to the GibbsThomsonHydrate subroutine
  ! Based on the TOUGH method of PE-salinity behavior, using reference
  ! parameters from Dickens and Quinby-Hunt, 1997
  !
  ! Author: David Fukuyama
  ! Date: 5/1/2020
  !

  implicit none

  PetscReal, intent(in) :: xmol
  PetscReal, intent(out) :: dTd

  PetscReal :: dTr = -0.37d0
  PetscReal :: xmol_ref = 0.03311d0

  dTd = dTr * log(1-xmol)/log(1-xmol_ref)

end subroutine HydrateSalinityOffset

! ************************************************************************** !

subroutine IceSalinityOffset(xmass,dTd)
  !
  ! Author: Michael Nole
  ! Date: 03/21/22
  !
  ! From Fujino et al., 1974

  implicit none

  PetscReal, intent(in) :: xmass
  PetscReal, intent(out) :: dTd

  dTd = -0.0575d0*(xmass*1000) + 0.000112d0*(xmass*1000)**2

end subroutine IceSalinityOffset

! ************************************************************************** !

subroutine HydrateEquilibrate(T,P,state,s_h,p_a,p_vap,p_sat,p_vap_brine, &
                           xag, xwg, xal, xsl, xwl, &
                           xmolag, xmolwg, xmolal, xmolsl, xmolwl, &
                           characteristic_curves, material_auxvar,option)
  !
  ! Computes equilibrium partitioning between CO2 and water following
  ! Spycher and Pruess, 2010, and between CH4 and water using a Henry's 
  ! constant.
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  use Option_module
  use EOS_Gas_module
  use Material_Aux_module
  use Characteristic_Curves_module

  implicit none

  PetscReal, intent(in) :: T ! temperature (C)
  PetscReal, intent(in) :: P ! liquid or gas pressure (Pa)
  PetscInt, intent(in) :: state ! State of the system
  PetscReal, intent(in) :: s_h ! Hydrate saturation
  PetscReal, intent(out) :: p_a ! partial pressure of Air (Pa)
  PetscReal, intent(out) :: p_vap ! partial pressure of water (Pa)
  PetscReal, intent(in) :: p_sat ! saturated brine vapor pressure (Pa)
  PetscReal, intent(in) :: p_vap_brine ! reduced vapor pressure (Pa)
  PetscReal, intent(out) :: xag ! mass fraction of air in gas phase
  PetscReal, intent(out) :: xwg ! mass fraction of water in gas phase
  PetscReal, intent(out) :: xal ! mass fraction of air in liquid phase
  PetscReal, intent(inout) :: xsl ! mass fraction of salt in liquid phase
  PetscReal, intent(out) :: xwl ! mass fraction of water in liquid phase
  PetscReal, intent(out) :: xmolag ! mole fraction of air in gas phase
  PetscReal, intent(out) :: xmolwg ! mole fraction of water in gas phase
  PetscReal, intent(out) :: xmolal ! mole fraction of air in liquid phase
  PetscReal, intent(out) :: xmolsl ! mole fraction of salt in liquid phase
  PetscReal, intent(out) :: xmolwl ! mole fraction of water in liquid phase
  class(characteristic_curves_type) :: characteristic_curves
  type(material_auxvar_type) :: material_auxvar
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
  PetscReal :: K_H

  PetscInt :: wid, acid, sid, lid, gid
  PetscInt :: i, k

  PetscReal, parameter :: epsilon = 1.d-14
  
  PetscErrorCode :: ierr

  wid = option%water_id
  acid = option%air_id
  sid = option%salt_id
  lid = option%liquid_phase
  gid = option%gas_phase

  p_a = 0.d0
  xag = 0.d0
  xwg = 0.d0
  xal = 0.d0
  xwl = 0.d0
  xmolag = 0.d0
  xmolwg = 0.d0
  xmolal = 0.0d0
  xmolsl = 0.d0
  xmolwl = 0.d0

  select case(hydrate_former)

  case (HYDRATE_FORMER_CO2)
    T_k = T + 273.15d0
    P_bar = max(P,1.01325d5)*1.d-5

    T_bound(1) = 99.d0
    T_bound(2) = 101.d0 !109.d0

    T_bound = T_bound + 273.15d0

    ! Salinity offset
    nacl_param = cxi(1)*T_k + cxi(2)/T_k + cxi(3)/(T_k**2)
    cl_param = clmb(1)*T_k + clmb(2)/T_k + clmb(3)/(T_k**2)

    ! Nacl mass fraction to molality
    xmolsl = 1.d3*(xsl/(1.d0-xsl))/hydrate_fmw_comp(3)
    xmol_na = xmolsl
    xmol_cl = xmolsl
    xmolwl = 1.d3 / hydrate_fmw_comp(1)
    xmolsl = xmolsl / (xmolsl + xmolwl)
    apc = (1.d0 + (xmol_na + xmol_cl)/xmolwl)*exp(2.d0*cl_param*xmol_na + &
           nacl_param*xmol_cl*xmol_na)

    ! Simplified solution, up to 275 C
    if (hydrate_spycher_simple) then
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
      ! if (state == L_STATE .or. state == HA_STATE) then
      !   call HydratePE(T,s_h, PE_hyd, dP,&
      !     characteristic_curves, material_auxvar, option)
      !   call HydrateGHSZSolubilityCorrection(T,P,dP,Hc)
      ! endif
      xmolal = pva * Hc
      xmolal = max(min(xmolal,1.d0),0.d0)
      xmolsl = (1.d0 - xmolal)*xmolsl
      xmolwl = 1.d0 - xmolal - xmolsl

      xmolwg = xmolwg * p_vap_brine / p_sat
      xmolag = 1.d0 - xmolwg
    elseif (hydrate_use_henry_co2) then

      call EOSGasHenry(T,p_sat,K_H,ierr)
      K_H = HydrateHenryCO2(T, xsl)
    
      ! if (state == L_STATE .or. state == HA_STATE) then
      !   call HydratePE(T,s_h, PE_hyd, dP,&
      !       characteristic_curves, material_auxvar, option)
      !   call HydrateGHSZSolubilityCorrection(T,P,dP,K_H)
      ! endif

      p_vap = p_vap_brine
      p_a = P - p_vap
    
      xmolal = p_a / K_H
      xmolwl = 1.d0 - xmolal
      xmolag =  p_a / P
      xmolwg = 1.d0 - xmolag

      xmolsl = 1.d3*(xsl/(1.d0-xsl))/hydrate_fmw_comp(3)
      xmol_na = xmolsl
      xmol_cl = xmolsl
      xmolsl = xmolsl / (xmolsl + 1.d3 / hydrate_fmw_comp(1))
      xmolsl = (1.d0 - xmolal)*xmolsl
      xmolwl = 1.d0 - xmolal - xmolsl
    elseif (T_k < T_bound(1)) then !Low temperature regime
      a = cac(1) + cac(2)*T_k
      b = cbc
  
      ! RKS EOS coefficients
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
      xmolal = coeff_b * (1.d0 - xmolwg)
  
      ! Vapor pressure lowering
      xmolwg = xmolwg * p_vap_brine / p_sat
      xmolag = 1.d0 - xmolwg

      K_H = 1.d0
      ! if (state == L_STATE .or. state == HA_STATE) then
      !   call HydratePE(T,s_h, PE_hyd, dP,&
      !     characteristic_curves, material_auxvar, option)
      !   call HydrateGHSZSolubilityCorrection(T,P,dP,K_H)
      ! endif
      xmolal = xmolal * K_H
      xmolsl = (1.d0 - xmolal)*xmolsl
      xmolwl = 1.d0 - xmolal - xmolsl
  
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
    if (xmolal < 1.d-16) xmolal = 0.d0
    if (xmolsl < 1.d-16) xmolsl = 0.d0
    if (xmolwg < 1.d-16) xmolwg = 0.d0
    if (xmolag < 1.d-16) xmolag = 0.d0

    ! Component partial pressures. 
    ! MAN: This part messes up hmode
    p_a = xmolag * P
    p_vap = xmolwg * P

  case default
    call EOSGasHenry(T,p_sat,K_H,ierr)

    p_vap = p_vap_brine
    p_a = P - p_vap

    xmolal = p_a / K_H
    xmolwl = 1.d0 - xmolal
    xmolag =  p_a / P
    xmolwg = 1.d0 - xmolag

    xmolsl = 1.d3*(xsl/(1.d0-xsl))/hydrate_fmw_comp(3)
    xmol_na = xmolsl
    xmol_cl = xmolsl
    xmolsl = xmolsl / (xmolsl + 1.d3 / hydrate_fmw_comp(1))
    xmolsl = (1.d0 - xmolal)*xmolsl
    xmolwl = 1.d0 - xmolal - xmolsl

  end select

  ! Mass Fractions

  fmw_gas = xmolwg * hydrate_fmw_comp(1) + xmolag * hydrate_fmw_comp(2)
  fmw_liq = xmolwl * hydrate_fmw_comp(1) + xmolal * hydrate_fmw_comp(2) + &
            xmolsl * hydrate_fmw_comp(3)

  ! Gas Phase
  xwg = xmolwg * hydrate_fmw_comp(1) / fmw_gas
  xag = xmolag * hydrate_fmw_comp(2) / fmw_gas
  ! Liquid Phase
  xwl = xmolwl * hydrate_fmw_comp(1) / fmw_liq
  xal = xmolal * hydrate_fmw_comp(2) / fmw_liq
  xsl = xmolsl * hydrate_fmw_comp(3) / fmw_liq

end subroutine HydrateEquilibrate

! ************************************************************************** !

subroutine CubicRootsNickalls(a,b,c,d,r1,r2,r3)
  !
  ! Computes roots of a cubic polynomial following Nickalls, 1993,
  ! A new approach to solving the cubic Cardans solution
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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

subroutine HydrateComputeSaltDensity(T,P,rho_s)
  !
  ! Computes NaCl density following Battistelli et al., 1997
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscReal, intent(out) :: rho_s

  rho_s = 2.165d3 * exp(-1.2d-4 * T + 4.d-11 * P)

end subroutine HydrateComputeSaltDensity

! ************************************************************************** !

subroutine HydrateComputeSurfaceTension(T,x_nacl,surface_tension)
  !
  ! Computes CO2-Water surface tension as a function of temperature
  ! and salt concentration following Abramzon and Gaukhberg, 1993.
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: x_nacl
  PetscReal, intent(out) :: surface_tension

  PetscReal :: molality

  molality = 1.d3*x_nacl/(hydrate_fmw_comp(3)*(1.d0-x_nacl))

  ! Pure water
  surface_tension = 1.d-3*(75.6592d0 - 1.40959d-1*T - 2.66317d-4*(T**2))
  ! With salt
  surface_tension = surface_tension + 1.57d-3*molality

  ! MAN: turn off surface tension effects
  surface_tension = CO2_REFERENCE_SURFACE_TENSION

end subroutine HydrateComputeSurfaceTension

! ************************************************************************** !

subroutine HydrateWaterSaturationPressure(T,P_sat)
  !
  ! Computes pure water saturation pressure following Meyer et al., 1993 and 
  ! Huang, 2018
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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

  ! Meyer et al.
  T_r = (T + 273.15) / T_c
  T_rx = 1.d0 - T_r

  P_sat = 0.d0
  do i = 1,5
    P_sat = P_sat + k(i) * (T_rx ** i)
  enddo
  P_sat = P_sat / ((1.d0 + k(6) * T_rx + k(7) * (T_rx **2)) * T_r)
  P_sat = P_sat - T_rx / (k(8) * (T_rx ** 2) + k(9))
  P_sat = exp(P_sat) * P_c

  ! Huang
  ! if (T < 0.d0) then
  !   P_sat = exp(43.494d0 - 6545.8d0 / (T + 278.d0)) / ((T + 868)**2)
  ! else
  !   P_sat = exp(34.494d0 - 4924.99d0 / (T + 273.1.d0)) / ((T + 105)**1.57d0)
  ! endif

end subroutine HydrateWaterSaturationPressure

! ************************************************************************** !

subroutine HydrateBrineSaturationPressure(T, x_salt, P_sat)
  !
  ! Computes brine saturation pressure following Haas, 1976
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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
  x_salt_molal = 1.d3 * x_salt / (hydrate_fmw_comp(3) * (1.d0 - x_salt))

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

  call HydrateWaterSaturationPressure(T_eq, P_sat)

  end subroutine HydrateBrineSaturationPressure
! ************************************************************************** !

subroutine HydrateWaterSubregion(T,P,isubr)
  !
  ! Computes subregion of water EOS following Meyer et al., 1993
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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
    call HydrateWaterSaturationPressure(T,P_sat)
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


end subroutine HydrateWaterSubregion

! ************************************************************************** !

subroutine HydrateWaterDensity(T,P,isubr,rho_l,rho_v,option)
  !
  ! Computes pure water density following Meyer et al., 1993
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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
    rho_l = 1.d3 * hydrate_fmw_comp(1) / (r_v * v_c)
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
    rho_v = 1.d3 * hydrate_fmw_comp(1) / (r_v * v_c)
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


end subroutine HydrateWaterDensity

! ************************************************************************** !

subroutine HydrateBrineDensity(T, P, x_s, rho_b, option)
  !
  ! Computes brine density following Phillips et al., 1983
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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
  x_s_molal = 1.d3 * x_s / (hydrate_fmw_comp(3) * (1.d0 - x_s))

  call HydrateWaterSaturationPressure(T, P_sat)
  P_w = max(P,P_sat)
  call HydrateWaterDensity(T, P_w, ONE_INTEGER, rho_l, rho_v, option)

  rho_l = 1.d-3 * rho_l
  spec_vol = 1.d0 / rho_l

  phi0 = c_h(1) + c_h(2) * spec_vol + c_h(3) * (spec_vol ** 2)
  phi = phi0 + (c_h(4) + c_h(5) * spec_vol) * &
               ((spec_vol/(v_c - spec_vol)) ** 2) * sqrt(x_s_molal)

  rho_b = (1.d3 + x_s_molal * hydrate_fmw_comp(3)) / &
          (1.d3 * spec_vol + x_s_molal * phi)

  ! kg/m^3
  rho_b = rho_b * 1.d3

end subroutine HydrateBrineDensity

! ************************************************************************** !

subroutine HydrateDensityCompositeLiquid(T,rho_b,x_co2, rho_l)
  !
  ! Computes density of the liquid phase as a funtion of brine density and
  ! CO2 concentration, Alendal and Drange, 2001
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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

  c_co2 = pv_co2 * rho_b * x_co2 / hydrate_fmw_comp(2)

  rho_l = rho_b / (1.d0 + c_co2 - x_co2)


end subroutine HydrateDensityCompositeLiquid

! ************************************************************************** !

subroutine HydrateViscosityWater(T, P, rho_w, visc, option)
  !
  ! Computes viscosity of pure water following Meyer et al., 1993
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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
    call HydrateWaterSubregion(T,P_inc,isubr)
    call HydrateWaterDensity(T,P_inc,isubr,rho_l,rho_vap,option)
    if ((1.d0 - dabs(rho_l / rho_w)) < (1.d0 - dabs(rho_vap/rho_w))) then
      rho_r2 = rho_l / rho_ref
    else
      rho_r2 = rho_vap / rho_ref
    endif
    chi = rho_r * (rho_r2 - rho_r) / dP_r
    if (chi >= 21.93d0) visc = visc * 0.922d0 * (chi ** 0.0263d0)
  endif

  visc = 1.d-6 * visc * visc_ref

end subroutine HydrateViscosityWater

! ************************************************************************** !

subroutine HydrateViscosityCO2(T, rho_co2, visc)
  !
  ! Computes viscosity of CO2 following Fenghour et al., 1998
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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

end subroutine HydrateViscosityCO2

! ************************************************************************** !

subroutine HydrateViscosityGas(visc_w, visc_a, xwg, xag, visc)
  !
  ! Computes gas phase viscosity following Reid et al., 1987
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  implicit none

  PetscReal, intent(in) :: visc_w
  PetscReal, intent(in) :: visc_a
  PetscReal, intent(in) :: xwg ! mole fraction
  PetscReal, intent(in) :: xag ! mole fraction
  PetscReal, intent(out) :: visc

  PetscReal :: phi_w, phi_a, chi_w, chi_a

  phi_w = ((1.d0 + sqrt(visc_a/visc_w) * &
        ((hydrate_fmw_comp(1)/hydrate_fmw_comp(2)) ** 2.5d-1)) **2) / &
        sqrt(8.d0 * (1.d0 + hydrate_fmw_comp(2)/hydrate_fmw_comp(1)))
  phi_a = ((1.d0 + sqrt(visc_w/visc_a) * &
        ((hydrate_fmw_comp(2)/hydrate_fmw_comp(1)) ** 2.5d-1)) **2) / &
        sqrt(8.d0 * (1.d0 + hydrate_fmw_comp(1)/hydrate_fmw_comp(2)))
  chi_w = xwg + xag * phi_a
  chi_a = xwg * phi_w + xag
  visc = xwg * visc_w / chi_w + xag * visc_a / chi_a

end subroutine HydrateViscosityGas

! ************************************************************************** !

subroutine HydrateViscosityBrine(T, x_salt, visc_w, visc_b)
  !
  ! Computes brine viscosity following Phillips et al., 1981
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: x_salt
  PetscReal, intent(in) :: visc_w
  PetscReal, intent(out) :: visc_b

  PetscReal, parameter :: s_a(5) = [0.0816d0, 0.0122d0, 0.000128d0, &
                                    0.000629d0,-0.7d0]
  PetscReal :: x_salt_molal

  x_salt_molal = 1.d3 * x_salt / (hydrate_fmw_comp(3) * (1.d0 - x_salt))

  visc_b = visc_w * (1.d0 + s_a(1) * x_salt_molal + s_a(2) * &
           (x_salt_molal ** 2) + s_a(3) * (x_salt_molal ** 3) + &
           s_a(4) * T * (1.d0 - exp(s_a(5) * x_salt_molal)))

end subroutine HydrateViscosityBrine

! ************************************************************************** !

subroutine HydrateViscosityLiquid(x_a, visc_b, visc_a, visc_l)
  !
  ! Computes composite liquid phase viscosity, Kimagai and Yokoyama, 1999
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  implicit none

  PetscReal, intent(in) :: x_a
  PetscReal, intent(in) :: visc_b
  PetscReal, intent(in) :: visc_a
  PetscReal, intent(out) :: visc_l

  PetscReal, parameter :: epsilon = 1.d-14

  visc_l = (1.d0 - x_a) * log(visc_b)
  if (visc_a > 1.d-14) visc_l = visc_l + x_a * log(visc_a)
  visc_l = exp(visc_l)

end subroutine HydrateViscosityLiquid

! ************************************************************************** !

function HydrateEnthalpyCompositeLiquid(T, x_salt, x_a, h_brine, h_a)
  !
  ! Computes composite liquid phase enthalpy, Battistelli et al., 1997
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: x_salt ! mass fraction
  PetscReal, intent(in) :: x_a ! mass fraction
  PetscReal, intent(in) :: h_brine ! J/kg
  PetscReal, intent(in) :: h_a ! J/kg
  PetscReal :: HydrateEnthalpyCompositeLiquid

  PetscReal :: dT, Hc, T_pert, Hc_pert, dHc, T_k, h_sol

  dT = 1.d-6
  Hc = HydrateHenryCO2(T, x_salt)
  T_pert = T + dT
  Hc_pert = HydrateHenryCO2(T_pert, x_salt)
  dHc = log(Hc_pert / Hc) / dT

  T_k = T + 273.15d0
  h_sol = -IDEAL_GAS_CONSTANT * 1.d3 * (T_k **2) * dHc / hydrate_fmw_comp(2)

  ! J/kg
  HydrateEnthalpyCompositeLiquid = max(1.d0-x_a,0.d0) * h_brine + &
                                x_a * (h_a + h_sol)

end function HydrateEnthalpyCompositeLiquid

! ************************************************************************** !

function HydrateHenryCO2(T, x_salt)
  !
  ! Computes Henry's coefficient for CO2 in brine, Battistelli et al., 1997
  ! Eq. 29
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  implicit none

  PetscReal, intent(in) :: T ! C
  PetscReal, intent(in) :: x_salt ! mass fraction
  PetscReal :: HydrateHenryCO2 ! Henrys constant

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

  HydrateHenryCO2 = Hc * &
              (1.d1 ** (1.d3 * x_salt / (hydrate_fmw_comp(3) * (1.d0 - x_salt)) * skb))


end function HydrateHenryCO2

! ************************************************************************** !

subroutine HydrateDiffusionCoeff(T,P,xsl,viscl,hydrate_parameter,option)
  !
  ! Computes CO2-Water diffusion coefficient, Cadogan et al., 2014
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  use Option_module

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscReal, intent(in) :: xsl
  PetscReal, intent(in) :: viscl
  type(hydrate_parameter_type) :: hydrate_parameter
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
  PetscInt :: lid, gid, acid, sid, wid
  PetscReal :: dlng, viscb, viscbr

  lid = option%liquid_phase
  gid = option%gas_phase

  wid = option%water_id
  acid = option%air_id
  sid = option%salt_id

  ! CO2 diffusion through the gas phase
  T_k = T + 273.15
  P_bar = P * 1.d-5

  eps = sqrt(c_w(2)*c_co2(2))
  sig = 5.d-1 * (c_w(1) + c_co2(1))
  T_r = T_k / eps
  omega = (c_a(1) / (T_r ** c_a(2))) + (c_a(3) / (exp(c_a(4) * T_r))) + &
          (c_a(5) / (exp(c_a(6) * T_r))) + (c_a(7) / (exp(c_a(8) * T_r)))
  w_mix = 2.d0 / ((1.d0 / hydrate_fmw_comp(2)) + (1.d0 / hydrate_fmw_comp(1)))
  Dwg = (3.03d0 - (9.8d-1 / sqrt(w_mix))) * 1.d-3 * &
                    (T_k ** 1.5d0) / (P_bar * sqrt(w_mix) * (sig ** 2) * &
                     omega) * 1.d-4

  ! CO2 diffusion through the liquid phase
  Dco2l = 3.5984d0 - 6.5113d-2*T_k + 2.0282D-4*(T_k**2)
  Dco2l = Dco2l*1.D-9

  ! Correct for NaCl
  Dco2l = Dco2l*(1.6678d0 - 1.2531d-1*(1.d3*(xsl/hydrate_fmw_comp(THREE_INTEGER)) / &
         (1.d0-xsl))) / 1.6678d0

  ! Salt diffusion through the liquid phase
  s_molal = 1.d3 * xsl / (hydrate_fmw_comp(THREE_INTEGER) * (1.d0 - xsl))

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
  call HydrateViscosityBrine(T_ref,xsl,viscw_ref,viscbr)
  Dnacl = Dnacl * (viscw_ref/viscbr) * (1.d0 + s_molal*(dlng))
  call HydrateViscosityBrine(T_k,xsl,viscl,viscb)
  Dnacl = Dnacl * (T_k/2.9815d2)*(viscbr/viscb)

  hydrate_parameter%diffusion_coefficient(acid,lid) = Dco2l
  hydrate_parameter%diffusion_coefficient(wid,gid) = Dwg
  hydrate_parameter%diffusion_coefficient(sid,lid) = Dnacl

end subroutine HydrateDiffusionCoeff

! ************************************************************************** !

subroutine HydrateScalePermPhi(hyd_auxvar, material_auxvar, global_auxvar, &
                               option)
  !
  ! Computes effective permeability and porosity as a function of precipitate
  ! saturations (Verma and Pruess, 1998)
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: hyd_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(global_auxvar_type) :: global_auxvar ! Stores mineral info
  type(option_type) :: option

  PetscReal :: phi_r ! zero-permeability limit fraction of porosity
  PetscReal :: f ! geometric factor
  PetscReal :: phi_0 ! initial fraction of porosity
  ! PetscReal :: theta
  ! PetscReal :: tao, omega
  PetscReal :: solid_sat_eff
  PetscInt :: pid

  pid = option%precipitate_phase

  !MAN: hard-code for now
  ! tao = 1.5d0
  phi_0 = material_auxvar%porosity_base
  phi_r = 8.d-1
  f = phi_r

  solid_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid) + hyd_auxvar%sat(pid)

   if (hydrate_perm_scaling) then
     select case (hydrate_perm_scaling_function)
       case(1) ! Dai and Seol, 2014
         hyd_auxvar%effective_permeability = max(1.d-5, &
                     (1.d0-solid_sat_eff)**3/(1.d0+2.d0*solid_sat_eff)**2)
       case default
     end select
   else
     hyd_auxvar%effective_permeability = 1.d0
   endif

  !hyd_auxvar%effective_porosity = max(hyd_auxvar%effective_porosity * &
  !                                    (1.d0 - solid_sat_eff), &
  !                                    hyd_auxvar%effective_porosity * phi_r, &
  !                                    1.d-12)

  !select case(permeability_reduction_model)

  !  case(ONE_INTEGER)
      ! Simplified Verma & Pruess
  !    hyd_auxvar%effective_permeability = ((hyd_auxvar%effective_porosity / &
  !                                           phi_0 - phi_r ) / &
  !                                           (1.d0 - phi_r )) ** tao
  !  case(TWO_INTEGER)
  !    ! Verma & Pruess model
  !    theta = max((1.d0 - solid_sat_eff - phi_r) / (1.d0 - phi_r) , &
  !                0.d0)
  !    omega = 1.d0 + (1.d0/f)/(1.d0/phi_r - 1.d0)
  !    hyd_auxvar%effective_permeability = &
  !                        (theta ** 2) * (1.d0 - f +  f/(omega**2)) / &
  !                        (1.d0 - f + f * (theta / (theta + omega - 1.d0)) ** 2)

  !end select

end subroutine HydrateScalePermPhi

! ************************************************************************** !

subroutine HydrateTortuosity(s_l, s_g, phi, tao_l, tao_g)
  !
  ! Computes tortuosity of the rock to gas and liquid phases
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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

end subroutine HydrateTortuosity

! ************************************************************************** !

subroutine HydrateComputeSatHysteresis(characteristic_curves, Pc, Sl_min, &
                                    beta_gl,rho_l, Sl, Sgt, option)
  !
  ! Compute saturation as a function of Pc including hysteretic effects.
  ! I believe this only works for monontonic Pc functions.
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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
  Sgt = 0.d0

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

end subroutine HydrateComputeSatHysteresis

! ************************************************************************** !

subroutine HydrateComputePcHysteresis(characteristic_curves, Sl, Sgt, beta_gl,&
                                      Pc, option)
  !
  ! Compute Pc as a function of saturation including trapped gas.
  ! I believe this only works for monontonic Pc functions.
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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

  if (hydrate_no_pc) then
    Pc = 0.d0 
    return
  endif
  
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


end subroutine HydrateComputePcHysteresis

! ************************************************************************** !

subroutine HydrateComputeEffectiveDiffusion(hydrate_parameter, hyd_auxvar, option)
  !
  ! Compute effective diffusion coefficients
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  use Option_module

  implicit none

  type(hydrate_parameter_type) :: hydrate_parameter
  type(hydrate_auxvar_type) :: hyd_auxvar
  type(option_type) :: option

  PetscInt :: lid, gid, wid, acid, sid, tgid
  PetscReal :: T_scaled

  lid = option%liquid_phase
  gid = option%gas_phase
  tgid = option%trapped_gas_phase

  wid = option%water_id
  acid = option%air_id
  sid = option%salt_id

  select case(hydrate_diffusion_model)

    case(ONE_INTEGER)

      ! Salt effective_diffusion_coeff in liquid
      T_scaled = (hyd_auxvar%temp + 273.15d0) / SALT_REFERENCE_TEMPERATURE
      hyd_auxvar%effective_diffusion_coeff(sid,lid) = &
                 hydrate_parameter%diffusion_coefficient(sid,lid) * T_scaled * &
                (LIQUID_REFERENCE_VISCOSITY / hyd_auxvar%visc(lid)) * &
                hyd_auxvar%tortuosity(lid) * hyd_auxvar%sat(lid) * &
                hyd_auxvar%effective_porosity

    case(TWO_INTEGER)
      ! MAN: not complete

    case(THREE_INTEGER)
      ! Salt effective_diffusion_coeff in liquid
      hyd_auxvar%effective_diffusion_coeff(sid,lid) = &
                 hyd_auxvar%tortuosity(lid) * &
                 hyd_auxvar%sat(lid) * hyd_auxvar%effective_porosity * &
                 hydrate_parameter%diffusion_coefficient(sid,lid)


  end select

  hyd_auxvar%effective_diffusion_coeff(acid,lid) = &
                 hyd_auxvar%tortuosity(lid) * &
                 hyd_auxvar%sat(lid) * hyd_auxvar%effective_porosity * &
                 hydrate_parameter%diffusion_coefficient(acid,lid)

  hyd_auxvar%effective_diffusion_coeff(wid,gid) = &
                 hyd_auxvar%tortuosity(gid) * &
                 (hyd_auxvar%sat(gid) - hyd_auxvar%sat(tgid)) * &
                 hyd_auxvar%effective_porosity * &
                 hydrate_parameter%diffusion_coefficient(wid,gid)

  hyd_auxvar%effective_diffusion_coeff(acid,gid) = &
                               hyd_auxvar%effective_diffusion_coeff(wid,gid)

end subroutine HydrateComputeEffectiveDiffusion

! ************************************************************************** !

subroutine HydrateSaltEnthalpy(T,H)

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
  H = 1.d3 * H / hydrate_fmw_comp(THREE_INTEGER)

end subroutine HydrateSaltEnthalpy

! ************************************************************************** !

subroutine HydrateComputeSaltSolubility(T, x_salt)
  !
  ! Computes solubility of NaCl in water. McKibbin and McNabb, 1993.
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: x_salt

  PetscReal, parameter :: coeff(3) = [2.6218d-1, 7.2d-5, 1.06d-6]

  x_salt = coeff(1) + coeff(2) * T + coeff(3) * T ** 2

end subroutine HydrateComputeSaltSolubility

! ************************************************************************** !

subroutine HydrateVaporPressureBrine(T,P_sat,Pc,rho_kg,x_salt,P_vap)
  !
  ! Computes the reduced vapor pressure of water following the Kelvin equation.
  !
  ! Author: Michael Nole
  ! Date: 01/29/24
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
  mw_mix = x_salt*hydrate_fmw_comp(3) + (1.d0-x_salt)*hydrate_fmw_comp(1)

  if (Pc > epsilon) then
    ! P_vap = P_sat * exp(- hydrate_fmw_comp(1) * (Pc**1.25d0) / &
    !          (rho_kg * IDEAL_GAS_CONSTANT * 1.d3 * T_k))
    P_vap = P_sat * exp(-1.d0 *hydrate_fmw_comp(1) * (Pc ** 1.25d0) / &
            (rho_kg * IDEAL_GAS_CONSTANT * 1.d3 * T_k))
  else
    P_vap = P_sat
  endif

end subroutine HydrateVaporPressureBrine


! ************************************************************************** !
end module Hydrate_Aux_module
