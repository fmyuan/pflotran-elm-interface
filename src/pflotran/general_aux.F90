module General_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none
  
  private 


  PetscBool, public :: general_print_state_transition = PETSC_TRUE
  PetscBool, public :: general_analytical_derivatives = PETSC_FALSE
  PetscBool, public :: general_immiscible = PETSC_FALSE
  PetscBool, public :: general_non_darcy_flow = PETSC_FALSE
  PetscReal, public :: general_phase_chng_epsilon = 1.d-6
  PetscBool, public :: general_restrict_state_chng = PETSC_FALSE
  PetscReal, public :: window_epsilon = 1.d-4 !0.d0
  PetscReal, public :: fmw_comp(2) = [FMWH2O,FMWAIR]
  PetscReal, public :: general_max_pressure_change = 5.d4
  PetscInt, public :: general_max_it_before_damping = UNINITIALIZED_INTEGER
  PetscReal, public :: general_damping_factor = 0.6d0
  PetscInt, public :: general_debug_cell_id = UNINITIALIZED_INTEGER
  PetscReal, public :: non_darcy_B = -0.78d0
#if defined(MATCH_TOUGH2)
  PetscBool, public :: general_temp_dep_gas_air_diff = PETSC_FALSE
  PetscBool, public :: general_diffuse_xmol = PETSC_FALSE
  PetscBool, public :: general_harmonic_diff_density = PETSC_TRUE
#else
  PetscBool, public :: general_diffuse_xmol = PETSC_TRUE
  PetscBool, public :: general_temp_dep_gas_air_diff = PETSC_TRUE
  PetscBool, public :: general_harmonic_diff_density = PETSC_TRUE
#endif
  PetscInt, public :: general_newton_iteration_number = 0
  PetscInt, public :: general_sub_newton_iter_num = 0

  PetscBool, public :: general_high_temp_ts_cut = PETSC_FALSE
  PetscBool, public :: general_using_newtontr = PETSC_FALSE
  PetscBool, public :: general_allow_state_change = PETSC_TRUE
  PetscBool, public :: general_state_changed = PETSC_FALSE
  PetscBool, public :: general_force_iteration = PETSC_FALSE
  PetscBool, public :: gen_chk_max_dpl_liq_state_only = PETSC_FALSE

  ! debugging
  PetscInt, public :: general_ni_count
  PetscInt, public :: general_ts_cut_count
  PetscInt, public :: general_ts_count

  ! thermodynamic state of fluid ids
  PetscInt, parameter, public :: NULL_STATE = 0
  PetscInt, parameter, public :: LIQUID_STATE = 1
  PetscInt, parameter, public :: GAS_STATE = 2
  PetscInt, parameter, public :: TWO_PHASE_STATE = 3
  PetscInt, parameter, public :: ANY_STATE = 4
  PetscInt, parameter, public :: MULTI_STATE = 5
  
  PetscInt, parameter, public :: PREV_TS = 1
  PetscInt, parameter, public :: PREV_IT = 2

  PetscInt, parameter, public :: GENERAL_LIQUID_PRESSURE_DOF = 1
  PetscInt, parameter, public :: GENERAL_GAS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: GENERAL_2PH_STATE_AIR_PRESSURE_DOF = 3
  PetscInt, parameter, public :: GENERAL_GAS_STATE_AIR_PRESSURE_DOF = 2
  PetscInt, parameter, public :: GENERAL_GAS_SATURATION_DOF = 2
  
  PetscInt, parameter, public :: GENERAL_ENERGY_DOF = 3
  PetscInt, parameter, public :: GENERAL_LIQUID_STATE_X_MOLE_DOF = 2
  
  PetscInt, parameter, public :: GENERAL_STATE_INDEX = 1
  PetscInt, parameter, public :: GENERAL_LIQUID_EQUATION_INDEX = 1
  PetscInt, parameter, public :: GENERAL_GAS_EQUATION_INDEX = 2
  PetscInt, parameter, public :: GENERAL_ENERGY_EQUATION_INDEX = 3
  
  PetscInt, parameter, public :: GENERAL_LIQUID_PRESSURE_INDEX = 2
  PetscInt, parameter, public :: GENERAL_GAS_PRESSURE_INDEX = 3
  PetscInt, parameter, public :: GENERAL_AIR_PRESSURE_INDEX = 4
  PetscInt, parameter, public :: GENERAL_MOLE_FRACTION_INDEX = 5
  PetscInt, parameter, public :: GENERAL_TEMPERATURE_INDEX = 6
  PetscInt, parameter, public :: GENERAL_GAS_SATURATION_INDEX = 7
  PetscInt, parameter, public :: GENERAL_LIQUID_FLUX_INDEX = 8
  PetscInt, parameter, public :: GENERAL_GAS_FLUX_INDEX = 9
  PetscInt, parameter, public :: GENERAL_ENERGY_FLUX_INDEX = 10
  PetscInt, parameter, public :: GENERAL_LIQUID_CONDUCTANCE_INDEX = 11
  PetscInt, parameter, public :: GENERAL_GAS_CONDUCTANCE_INDEX = 12
  PetscInt, parameter, public :: GENERAL_GAS_WATER_MOL_FRAC_INDEX = 13
  PetscInt, parameter, public :: GENERAL_MAX_INDEX = 13
  
  PetscInt, parameter, public :: GENERAL_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: GENERAL_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: GENERAL_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: GENERAL_UPDATE_FOR_BOUNDARY = 2
  PetscInt, parameter, public :: GENERAL_UPDATE_FOR_SS = 3

  PetscReal, parameter, public :: GENERAL_IMMISCIBLE_VALUE = 1.d-10
  PetscReal, parameter, public :: GENERAL_PRESSURE_SCALE = 1.d0
  
  ! these variables, which are global to general, can be modified
  PetscInt, public :: dof_to_primary_variable(3,3)
  PetscInt, public :: general_2ph_energy_dof = GENERAL_TEMPERATURE_INDEX
  PetscInt, public :: general_gas_air_mass_dof = GENERAL_AIR_PRESSURE_INDEX
  PetscBool, public :: general_isothermal = PETSC_FALSE
  PetscBool, public :: general_no_air = PETSC_FALSE
  
  type, public :: general_auxvar_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    PetscBool :: istatechng
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal, pointer :: den(:)    ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: den_kg(:) ! (iphase) kg/m^3 phase
    PetscReal :: temp
    PetscReal, pointer :: xmol(:,:) ! (icomp,iphase)
    PetscReal, pointer :: H(:) ! MJ/kmol
    PetscReal, pointer :: U(:) ! MJ/kmol
!    PetscReal, pointer :: dsat_dp(:,:)
!    PetscReal, pointer :: dden_dp(:,:)
!    PetscReal, pointer :: dsat_dT(:)
!    PetscReal, pointer :: dden_dT(:)
    PetscReal, pointer :: kr(:)
    PetscReal, pointer :: mobility(:) ! relative perm / kinematic viscosity
    PetscReal :: effective_porosity ! factors in compressibility
    PetscReal :: pert
!    PetscReal, pointer :: dmobility_dp(:)
    type(general_derivative_auxvar_type), pointer :: d
  end type general_auxvar_type

  type, public :: general_derivative_auxvar_type
    PetscReal :: pc_satg
    PetscReal :: por_p
    PetscReal :: denl_pl
    PetscReal :: denl_T
    PetscReal :: deng_pg
    PetscReal :: deng_pa
    PetscReal :: deng_T    
    PetscReal :: dengkg_pg
    PetscReal :: dengkg_T    
    PetscReal :: Ul_pl
    PetscReal :: Ul_T
    PetscReal :: Hl_pl
    PetscReal :: Hl_T
    PetscReal :: Ug_pg
    PetscReal :: Ug_pa
    PetscReal :: Ug_T
    PetscReal :: Hg_pg
    PetscReal :: Hg_pa
    PetscReal :: Hg_T 

    PetscReal :: Hv
    PetscReal :: Ha
    PetscReal :: Uv
    PetscReal :: Ua
    PetscReal :: Hv_pg
    PetscReal :: Hv_pa
    PetscReal :: Hv_T 
    PetscReal :: Ha_pg
    PetscReal :: Ha_pa
    PetscReal :: Ha_T 
    PetscReal :: Uv_pg
    PetscReal :: Uv_pa
    PetscReal :: Uv_T 
    PetscReal :: Ua_pg
    PetscReal :: Ua_pa
    PetscReal :: Ua_T 
    PetscReal :: denv
    PetscReal :: dena 
    PetscReal :: denv_T 
    PetscReal :: dena_T 
    PetscReal :: denv_pg
    PetscReal :: dena_pg
    PetscReal :: Hc 
    
    PetscReal :: psat_p
    PetscReal :: psat_T
    PetscReal :: pv_p
    PetscReal :: pv_pa
    PetscReal :: pv_T
    PetscReal :: Hc_p
    PetscReal :: Hc_T
    PetscReal :: mobilityl_pl
    PetscReal :: mobilityl_T
    PetscReal :: mobilityl_satg
    PetscReal :: mobilityg_pg
    PetscReal :: mobilityg_T
    PetscReal :: mobilityg_satg
    PetscReal :: mobilityg_pa
    PetscReal :: mug
    PetscReal :: mug_T
    PetscReal :: mug_pg
    PetscReal :: xmol_p(2,2)
    PetscReal :: xmol_T(2,2)    
  end type general_derivative_auxvar_type
  
  type, public :: general_parameter_type
    PetscReal, pointer :: diffusion_coefficient(:) ! (iphase)
    PetscReal :: newton_inf_scaled_res_tol
    PetscBool :: check_post_converged
  end type general_parameter_type

  type, public :: general_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(general_parameter_type), pointer :: general_parameter
    type(general_auxvar_type), pointer :: auxvars(:,:)
    type(general_auxvar_type), pointer :: auxvars_bc(:)
    type(general_auxvar_type), pointer :: auxvars_ss(:,:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type general_type

  interface GeneralAuxVarDestroy
    module procedure GeneralAuxVarSingleDestroy
    module procedure GeneralAuxVarArray1Destroy
    module procedure GeneralAuxVarArray2Destroy
  end interface GeneralAuxVarDestroy
  
  interface GeneralOutputAuxVars
    module procedure GeneralOutputAuxVars1
    module procedure GeneralOutputAuxVars2
  end interface GeneralOutputAuxVars
  
  public :: GeneralAuxCreate, &
            GeneralAuxDestroy, &
            GeneralAuxSetEnergyDOF, &
            GeneralAuxSetAirMassDOF, &
            GeneralAuxVarCompute, &
            GeneralAuxVarInit, &
            GeneralAuxVarCopy, &
            GeneralAuxVarDestroy, &
            GeneralAuxVarStrip, &
            GeneralAuxVarUpdateState, &
            GeneralAuxVarPerturb, &
            GeneralPrintAuxVars, &
            GeneralOutputAuxVars

contains

! ************************************************************************** !

function GeneralAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  use Option_module

  implicit none

  type(option_type) :: option
    
  type(general_type), pointer :: GeneralAuxCreate
  
  type(general_type), pointer :: aux

  dof_to_primary_variable(1:3,1:3) = &
    reshape([GENERAL_LIQUID_PRESSURE_INDEX, GENERAL_MOLE_FRACTION_INDEX, &
             GENERAL_TEMPERATURE_INDEX, &
             GENERAL_GAS_PRESSURE_INDEX, general_gas_air_mass_dof, &
             GENERAL_TEMPERATURE_INDEX, &
             GENERAL_GAS_PRESSURE_INDEX, GENERAL_GAS_SATURATION_INDEX, &
             general_2ph_energy_dof],shape(dof_to_primary_variable))
  
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

  allocate(aux%general_parameter)
  allocate(aux%general_parameter%diffusion_coefficient(option%nphase))
  !geh: there is no point in setting default lquid diffusion coeffcient values 
  !     here as they will be overwritten by the fluid property defaults.
  aux%general_parameter%diffusion_coefficient(LIQUID_PHASE) = &
                                                           UNINITIALIZED_DOUBLE
  aux%general_parameter%diffusion_coefficient(GAS_PHASE) = 2.13d-5
  aux%general_parameter%newton_inf_scaled_res_tol = 1.d-50
  aux%general_parameter%check_post_converged = PETSC_FALSE
  
  GeneralAuxCreate => aux
  
end function GeneralAuxCreate

! ************************************************************************** !

subroutine GeneralAuxVarInit(auxvar,allocate_derivative,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  use Option_module

  implicit none
  
  type(general_auxvar_type) :: auxvar
  PetscBool :: allocate_derivative
  type(option_type) :: option

  auxvar%istate_store = NULL_STATE
  auxvar%istatechng = PETSC_FALSE
  auxvar%temp = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%pert = 0.d0
  
  allocate(auxvar%pres(option%nphase+FOUR_INTEGER))
  auxvar%pres = 0.d0
  allocate(auxvar%sat(option%nphase))
  auxvar%sat = 0.d0
  allocate(auxvar%den(option%nphase))
  auxvar%den = 0.d0
  allocate(auxvar%den_kg(option%nphase))
  auxvar%den_kg = 0.d0
  allocate(auxvar%xmol(option%nflowspec,option%nphase))
  auxvar%xmol = 0.d0
  allocate(auxvar%H(option%nphase))
  auxvar%H = 0.d0
  allocate(auxvar%U(option%nphase))
  auxvar%U = 0.d0
  allocate(auxvar%mobility(option%nphase))
  auxvar%mobility = 0.d0
  allocate(auxvar%kr(option%nphase))
  auxvar%kr = 0.d0
  if (allocate_derivative) then
    allocate(auxvar%d)
    auxvar%d%pc_satg = 0.d0
    auxvar%d%por_p = 0.d0
    auxvar%d%denl_pl = 0.d0
    auxvar%d%denl_T = 0.d0
    auxvar%d%deng_pg = 0.d0
    auxvar%d%deng_pa = 0.d0
    auxvar%d%deng_T = 0.d0
    auxvar%d%dengkg_pg = 0.d0
    auxvar%d%dengkg_T = 0.d0
    auxvar%d%Ul_pl = 0.d0
    auxvar%d%Ul_T = 0.d0
    auxvar%d%Hl_pl = 0.d0
    auxvar%d%Hl_T = 0.d0
    auxvar%d%Ug_pg = 0.d0
    auxvar%d%Ug_pa = 0.d0
    auxvar%d%Ug_T = 0.d0
    auxvar%d%Hg_pg = 0.d0
    auxvar%d%Hg_pa = 0.d0
    auxvar%d%Hg_T = 0.d0

    auxvar%d%Hv = 0.d0
    auxvar%d%Ha = 0.d0
    auxvar%d%Uv = 0.d0
    auxvar%d%Ua = 0.d0
    auxvar%d%Hv_pg = 0.d0
    auxvar%d%Hv_pa = 0.d0
    auxvar%d%Hv_T = 0.d0
    auxvar%d%Ha_pg = 0.d0
    auxvar%d%Ha_pa = 0.d0
    auxvar%d%Ha_T = 0.d0
    auxvar%d%Uv_pg = 0.d0
    auxvar%d%Uv_pa = 0.d0
    auxvar%d%Uv_T = 0.d0
    auxvar%d%Ua_pg = 0.d0
    auxvar%d%Ua_pa = 0.d0
    auxvar%d%Ua_T = 0.d0
    auxvar%d%denv = 0.d0
    auxvar%d%dena = 0.d0
    auxvar%d%denv_T = 0.d0
    auxvar%d%dena_T = 0.d0
    auxvar%d%denv_pg = 0.d0
    auxvar%d%dena_pg = 0.d0
    auxvar%d%Hc = 0.d0
        
    auxvar%d%psat_p = 0.d0
    auxvar%d%psat_T = 0.d0
    auxvar%d%pv_p = 0.d0
    auxvar%d%pv_pa = 0.d0
    auxvar%d%pv_T = 0.d0
    auxvar%d%Hc_p = 0.d0
    auxvar%d%Hc_T = 0.d0
    auxvar%d%mug = 0.d0
    auxvar%d%mug_T = 0.d0
    auxvar%d%mug_pg = 0.d0
    auxvar%d%mobilityl_pl = 0.d0
    auxvar%d%mobilityl_T = 0.d0
    auxvar%d%mobilityl_satg = 0.d0
    auxvar%d%mobilityg_pg = 0.d0
    auxvar%d%mobilityg_T = 0.d0
    auxvar%d%mobilityg_satg = 0.d0
    auxvar%d%mobilityg_pa = 0.d0
    auxvar%d%xmol_p = 0.d0
    auxvar%d%xmol_T = 0.d0
  else
    nullify(auxvar%d)
  endif
  
end subroutine GeneralAuxVarInit

! ************************************************************************** !

subroutine GeneralAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  use Option_module

  implicit none
  
  type(general_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%istate_store = auxvar%istate_store
  auxvar2%istatechng = auxvar%istatechng
  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%xmol = auxvar%xmol
  auxvar2%H = auxvar%H
  auxvar2%U = auxvar%U
  auxvar2%mobility = auxvar%mobility
  auxvar2%kr = auxvar%kr
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%pert = auxvar%pert

end subroutine GeneralAuxVarCopy

! ************************************************************************** !

subroutine GeneralAuxSetEnergyDOF(energy_keyword,option)
  ! 
  ! Sets the two phase primary dependent variable for energy based on user 
  ! input.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/25/14
  ! 
  use Option_module
  use String_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: energy_keyword
  type(option_type) :: option

  call StringToUpper(energy_keyword)
  select case(energy_keyword)
    case('TEMPERATURE')
      general_2ph_energy_dof = GENERAL_TEMPERATURE_INDEX
    case('AIR_PRESSURE')
      general_2ph_energy_dof = GENERAL_AIR_PRESSURE_INDEX
    case default
      option%io_buffer = 'Energy Keyword: ' // trim(energy_keyword) // &
                          ' not recognized in General Mode'    
      call PrintErrMsg(option)
  end select

end subroutine GeneralAuxSetEnergyDOF

! ************************************************************************** !
subroutine GeneralAuxSetAirMassDOF(air_mass_keyword,option)
  !
  !Sets the gas phase primary variable for the air mass equation based on
  !user input.
  !
  !Author: Michael Nole
  !Date: 06/13/19
  !

  use Option_module
  use String_module

  implicit none

  character(len=MAXWORDLENGTH) :: air_mass_keyword
  type(option_type) :: option

  call StringToUpper(air_mass_keyword)
  select case(air_mass_keyword)
    case('WATER_MOL_FRAC')
      general_gas_air_mass_dof = GENERAL_GAS_WATER_MOL_FRAC_INDEX
    case('AIR_PRESSURE')
      general_gas_air_mass_dof = GENERAL_AIR_PRESSURE_INDEX
  end select

end subroutine GeneralAuxSetAirMassDOF
! ************************************************************************** !

subroutine GeneralAuxVarCompute(x,gen_auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use Creep_Closure_module
  use Fracture_module
  use WIPP_module
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(creep_closure_type), pointer :: creep_closure
  PetscInt :: natural_id

  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: cell_pressure, water_vapor_pressure
  PetscReal :: den_water_vapor, den_kg_water_vapor
  PetscReal :: u_water_vapor, h_water_vapor
  PetscReal :: den_air, h_air, u_air
  PetscReal :: xmol_air_in_gas, xmol_water_in_gas
  PetscReal :: krl, visl, dvis_dp, dvis_dT, dvis_dpa
  PetscReal :: dkrl_dsatl, dkrl_dsatg
  PetscReal :: dkrg_dsatl, dkrg_dsatg
  PetscReal :: krg, visg
  PetscReal :: K_H_tilde
  PetscReal :: guess, dummy
  PetscInt :: apid, cpid, vpid, spid
  PetscReal :: NaN
  PetscReal :: creep_closure_time
  PetscReal :: xmass_air_in_gas
  PetscReal :: Ugas_J_kg, Hgas_J_kg
  PetscReal :: Uair_J_kg, Hair_J_kg
  PetscReal :: Uvapor_J_kg, Hvapor_J_kg
  PetscReal :: Hg_mixture_fractioned  
  PetscReal :: aux(1)
  PetscReal :: hw, hw_dp, hw_dT
  PetscReal :: dpor_dp
  PetscReal :: one_over_dw
  PetscReal :: tempreal, tempreal2, tempreal3
  PetscReal :: dpair_dT, dpair_dpgas
  PetscReal :: dden_air_dT, dden_air_dpa, dden_air_dpg
  PetscReal :: du_air_dT, dh_air_dT
  PetscReal :: du_air_dpa, dh_air_dpa
  PetscReal :: dden_water_vapor_dpv, dden_water_vapor_dT
  PetscReal :: dh_water_vapor_dpv, dh_water_vapor_dT
  PetscReal :: du_water_vapor_dpv, du_water_vapor_dT
  PetscReal :: dpc_dsatl
  character(len=8) :: state_char
  PetscErrorCode :: ierr
  PetscErrorCode :: eos_henry_ierr
  
  PetscReal :: sigma

  ! from init.F90
!  option%nphase = 2
!  option%liquid_phase = 1  ! liquid_pressure
!  option%gas_phase = 2     ! gas_pressure

!  option%air_pressure_id = 3
!  option%capillary_pressure_id = 4
!  option%vapor_pressure_id = 5

!  option%water_id = 1
!  option%air_id = 2
!  option%energy_id = 3

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id
  
  eos_henry_ierr = 0
  
  
#ifdef DEBUG_GENERAL  
  ! create a NaN
  NaN = InitToNan()

  gen_auxvar%H = NaN
  gen_auxvar%U = NaN
  gen_auxvar%pres = NaN
  gen_auxvar%sat = NaN
  gen_auxvar%den = NaN
  gen_auxvar%den_kg = NaN
  gen_auxvar%xmol = NaN
  gen_auxvar%effective_porosity = NaN
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      state_char = 'L'
    case(GAS_STATE)
      state_char = 'G'
    case(TWO_PHASE_STATE)
      state_char = '2P'
  end select
#else
  !geh: do not initialize gen_auxvar%temp as the previous value is used as the
  !     initial guess for two phase.
  gen_auxvar%H = 0.d0
  gen_auxvar%U = 0.d0
  gen_auxvar%pres = 0.d0
  gen_auxvar%sat = 0.d0
  gen_auxvar%den = 0.d0
  gen_auxvar%den_kg = 0.d0
  gen_auxvar%xmol = 0.d0
  gen_auxvar%effective_porosity = 0.d0
#endif  
  if (associated(gen_auxvar%d)) then
    call GeneralZeroAnalyticalDeriv(gen_auxvar)
  endif
  gen_auxvar%mobility = 0.d0
  gen_auxvar%kr = 0.d0

#if 0
  if (option%iflag >= GENERAL_UPDATE_FOR_ACCUM) then
    if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
        natural_id, x(1:3), trim(state_char)
    else
!      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
!        -1*natural_id, x(1:3), trim(state_char)
    endif
  endif
#endif
  
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      gen_auxvar%pres(lid) = x(GENERAL_LIQUID_PRESSURE_DOF)
      gen_auxvar%xmol(acid,lid) = x(GENERAL_LIQUID_STATE_X_MOLE_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)

      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      ! with the gas state, we must calculate the mole fraction of air in 
      ! in the liquid phase, even though the liquid phase does not exist
      ! due to air diffusion between neighboring GAS and LIQUID cells (this
      ! is more of an issue on a boundary condition).  this is not 
      ! necessary for water since we do not calculate water diffusion 
      ! explicitly.  set mole fractions to zero in gas phase.
      gen_auxvar%xmol(:,gid) = 0.d0
      gen_auxvar%sat(lid) = 1.d0
      gen_auxvar%sat(gid) = 0.d0

      if (associated(gen_auxvar%d)) then
        call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                        gen_auxvar%pres(spid), &
                                        gen_auxvar%d%psat_T,ierr)
        gen_auxvar%d%psat_p = 0.d0
        call EOSGasHenry(gen_auxvar%temp,gen_auxvar%pres(spid), &
                          gen_auxvar%d%psat_p,gen_auxvar%d%psat_T, &
                          K_H_tilde,gen_auxvar%d%Hc_p,gen_auxvar%d%Hc_T, &
                          eos_henry_ierr)
        gen_auxvar%d%Hc = K_H_tilde
      else
        call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                        gen_auxvar%pres(spid),ierr)
      !geh: Henry_air_xxx returns K_H in units of Pa, but I am not confident
      !     that K_H is truly K_H_tilde (i.e. p_g * K_H).
        call EOSGasHenry(gen_auxvar%temp,gen_auxvar%pres(spid),K_H_tilde, &
                         eos_henry_ierr)
    endif
      gen_auxvar%pres(gid) = max(gen_auxvar%pres(lid),gen_auxvar%pres(spid))
      gen_auxvar%pres(apid) = K_H_tilde*gen_auxvar%xmol(acid,lid)
      ! need vpres for liq -> 2ph check
      
      ! at this point, if the liquid pressure is negative, we have to go to 
      ! two phase or everything blows up:
      if (gen_auxvar%pres(gid) <= 0.d0) then
        write(option%io_buffer,'(''Negative gas pressure at cell '', &
          & i8,'' in GeneralAuxVarCompute(LIQUID_STATE).  Attempting bailout.'')') &
          natural_id
!        call PrintErrMsgByRank(option)
        call PrintMsgByRank(option)
        ! set vapor pressure to just under saturation pressure
        gen_auxvar%pres(vpid) = 0.5d0*gen_auxvar%pres(spid)
        ! set gas pressure to vapor pressure + air pressure
        gen_auxvar%pres(gid) = gen_auxvar%pres(vpid) + gen_auxvar%pres(apid)
        ! capillary pressure won't matter here.        
      else
        gen_auxvar%pres(vpid) = gen_auxvar%pres(lid) - gen_auxvar%pres(apid)
      endif
      gen_auxvar%pres(cpid) = 0.d0
      
      if (associated(gen_auxvar%d)) then
        gen_auxvar%d%pv_p = 1.d0 - gen_auxvar%d%Hc_p*gen_auxvar%xmol(acid,lid)
        gen_auxvar%d%pv_T = -gen_auxvar%d%Hc_T*gen_auxvar%xmol(acid,lid)
      endif
      
    case(GAS_STATE)
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      
      !man
      !gen_auxvar%pres(gid) = max(0.d0,gen_auxvar%pres(gid))
      
      if (general_gas_air_mass_dof == GENERAL_AIR_PRESSURE_INDEX) then
        gen_auxvar%pres(apid) = x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
        gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / &
                                   gen_auxvar%pres(gid)
        gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)
      else
        gen_auxvar%xmol(wid,gid) = x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
        gen_auxvar%xmol(acid,gid) = 1.d0 - gen_auxvar%xmol(wid,gid)
        gen_auxvar%pres(apid) = gen_auxvar%xmol(acid,gid) * &
                                gen_auxvar%pres(gid)
      endif

      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)

      gen_auxvar%sat(lid) = 0.d0
      gen_auxvar%sat(gid) = 1.d0
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)
      ! need to set mole fractions in liquid phase in equilibrium with
      ! water saturated with air in order to accommodate air diffusion between
      ! GAS_STATE cell and TWO_PHASE/LIQUID_STATE cells as air should still
      ! diffuse through the liquid phase.
      if (associated(gen_auxvar%d)) then
        call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                        gen_auxvar%pres(spid), &
                                        gen_auxvar%d%psat_T,ierr)
        gen_auxvar%d%psat_p = 0.d0
        call EOSGasHenry(gen_auxvar%temp,gen_auxvar%pres(spid), &
                          gen_auxvar%d%psat_p,gen_auxvar%d%psat_T, &
                          K_H_tilde,gen_auxvar%d%Hc_p,gen_auxvar%d%Hc_T, &
                          eos_henry_ierr)
        gen_auxvar%d%Hc = K_H_tilde
      else
        call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                        gen_auxvar%pres(spid),ierr)
        call EOSGasHenry(gen_auxvar%temp,gen_auxvar%pres(spid),K_H_tilde, &
                         eos_henry_ierr)
      endif
      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      ! set water mole fraction to zero as there is no water in liquid phase
      gen_auxvar%xmol(wid,lid) = 0.d0
      
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)

       
      ! we have to have a liquid pressure to counter a neighboring 
      ! liquid pressure.  Set to gas pressure.
!      gen_auxvar%pres(lid) = gen_auxvar%pres(gid)
!      gen_auxvar%pres(cpid) = 0.d0

      call characteristic_curves%saturation_function% &
             CapillaryPressure(gen_auxvar%sat(lid), &
                               gen_auxvar%pres(cpid),dpc_dsatl,option)                             
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - &
                             gen_auxvar%pres(cpid)
                             
      if (associated(gen_auxvar%d)) then
        option%io_buffer = 'Analytical derivatives for gas state cells in &
          &GENERAL mode have a bug. Please use numerical derivatives until &
          &the bug is resolved.'
        call PrintErrMsg(option)
        dpair_dT = 0.d0
        dpair_dpgas = 0.d0
        gen_auxvar%d%pv_p = 1.d0
        gen_auxvar%d%pv_pa = -1.d0
        gen_auxvar%d%pv_T = 0.d0
      
        gen_auxvar%d%xmol_p(acid,gid) = -gen_auxvar%xmol(acid,gid)/gen_auxvar%pres(gid)
        gen_auxvar%d%xmol_p(wid,gid) = -gen_auxvar%d%xmol_p(acid,gid)
        ! we hijack the liquid phase for air pressure
        ! this could be pushed to where it is used
        gen_auxvar%d%xmol_p(acid,lid) = 1.d0/gen_auxvar%pres(gid)
        gen_auxvar%d%xmol_p(wid,lid) = -gen_auxvar%d%xmol_p(acid,lid)
      endif                             

    case(TWO_PHASE_STATE)

      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      
      !man
      !gen_auxvar%pres(gid) = max(0.d0,gen_auxvar%pres(gid))
      
      gen_auxvar%sat(gid) = x(GENERAL_GAS_SATURATION_DOF)
      
      if (gen_auxvar%istatechng) then
        gen_auxvar%sat(gid) = max(0.d0,gen_auxvar%sat(gid))
        gen_auxvar%sat(gid) = min(1.d0,gen_auxvar%sat(gid))
      endif
      
      if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
        gen_auxvar%temp = x(GENERAL_ENERGY_DOF)
        if (associated(gen_auxvar%d)) then
          call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                          gen_auxvar%pres(spid), &
                                          gen_auxvar%d%psat_T,ierr)
          gen_auxvar%d%psat_p = 0.d0
          call EOSGasHenry(gen_auxvar%temp,gen_auxvar%pres(spid), &
                           gen_auxvar%d%psat_p,gen_auxvar%d%psat_T, &
                           K_H_tilde,gen_auxvar%d%Hc_p,gen_auxvar%d%Hc_T, &
                           eos_henry_ierr)
          gen_auxvar%d%Hc = K_H_tilde
        else
          call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                          gen_auxvar%pres(spid),ierr)
          call EOSGasHenry(gen_auxvar%temp,gen_auxvar%pres(spid),K_H_tilde, &
                           eos_henry_ierr)
        endif
        if (general_immiscible) then
          gen_auxvar%pres(spid) = GENERAL_IMMISCIBLE_VALUE
        endif
        gen_auxvar%pres(vpid) = gen_auxvar%pres(spid)
        gen_auxvar%pres(apid) = gen_auxvar%pres(gid) - gen_auxvar%pres(vpid)
        if (associated(gen_auxvar%d)) then
          dpair_dT = -1.d0*gen_auxvar%d%psat_T
          dpair_dpgas = 1.d0
        endif
      else
        gen_auxvar%pres(apid) = x(GENERAL_ENERGY_DOF)
        gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)
      
        gen_auxvar%pres(spid) = gen_auxvar%pres(vpid)
        guess = gen_auxvar%temp
        call EOSWaterSaturationTemperature(gen_auxvar%temp, &
                                           gen_auxvar%pres(spid),dummy, &
                                           guess,ierr)
      endif

      gen_auxvar%sat(lid) = 1.d0 - gen_auxvar%sat(gid)
      
      call characteristic_curves%saturation_function% &
             CapillaryPressure(gen_auxvar%sat(lid), &
                               gen_auxvar%pres(cpid),dpc_dsatl,option) 
      
      !man: IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(gen_auxvar%temp,sigma)
      endif
      gen_auxvar%pres(cpid) = gen_auxvar%pres(cpid)*sigma
      
      
      if (associated(gen_auxvar%d)) then
        ! for now, calculate derivative through finite differencing
#if 0
      !TODO(geh): make an analytical derivative
        tempreal = 1.d-6 * gen_auxvar%sat(lid)
        tempreal2 = gen_auxvar%sat(lid) + tempreal
        call characteristic_curves%saturation_function% &
            CapillaryPressure(material_auxvar,tempreal2,tempreal3,dpc_dsatl, &
                               option)
        gen_auxvar%d%pc_satg = -1.d0*(tempreal3-gen_auxvar%pres(cpid))/tempreal
#else
        gen_auxvar%d%pc_satg = -1.d0*dpc_dsatl
#endif
      endif
!      gen_auxvar%pres(cpid) = 0.d0
      
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      if (general_immiscible) then
        gen_auxvar%xmol(acid,lid) = GENERAL_IMMISCIBLE_VALUE
      endif
      
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / gen_auxvar%pres(gid)
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)
      
      if (associated(gen_auxvar%d)) then
        gen_auxvar%d%pv_T = gen_auxvar%d%psat_T
        gen_auxvar%d%pv_p = gen_auxvar%d%psat_p
        gen_auxvar%d%xmol_p(acid,lid) = dpair_dpgas/K_H_tilde - &
          gen_auxvar%xmol(acid,lid) / K_H_tilde * gen_auxvar%d%Hc_p
        gen_auxvar%d%xmol_p(wid,lid) = -1.d0*gen_auxvar%d%xmol_p(acid,lid)
        gen_auxvar%d%xmol_p(acid,gid) = dpair_dpgas/gen_auxvar%pres(gid) - &
          gen_auxvar%xmol(acid,gid)/gen_auxvar%pres(gid)
        gen_auxvar%d%xmol_p(wid,gid) = -1.d0*gen_auxvar%d%xmol_p(acid,gid)
        
        gen_auxvar%d%xmol_T(acid,lid) = dpair_dT/K_H_tilde - &
          gen_auxvar%xmol(acid,lid) / K_H_tilde * gen_auxvar%d%Hc_T
        gen_auxvar%d%xmol_T(wid,lid) = -1.d0*gen_auxvar%d%xmol_T(acid,lid)
        gen_auxvar%d%xmol_T(acid,gid) = dpair_dT/gen_auxvar%pres(gid)
        gen_auxvar%d%xmol_T(wid,gid) = -1.d0*gen_auxvar%d%xmol_T(acid,gid)
      endif
    case default
      write(option%io_buffer,*) global_auxvar%istate
      option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
        ') not recognized in GeneralAuxVarCompute.'
      call PrintErrMsgByRank(option)

  end select

  if (eos_henry_ierr /= 0) then
     call GeneralEOSGasError(natural_id,eos_henry_ierr,gen_auxvar,option)
  endif

  
  cell_pressure = max(gen_auxvar%pres(lid),gen_auxvar%pres(gid), &
                      gen_auxvar%pres(spid))
        
  ! calculate effective porosity as a function of pressure
  if (option%iflag /= GENERAL_UPDATE_FOR_BOUNDARY) then
    dpor_dp = 0.d0
    gen_auxvar%effective_porosity = material_auxvar%porosity_base
#if 0
!geh this code is no longer valid
    if (associated(material_auxvar%fracture) .and. & 
      material_auxvar%fracture%setup) then
      ! The initiating pressure and maximum pressure must be calculated
      ! before fracture function applies - Heeho
      call FractureInitialSetup(material_auxvar,cell_pressure)
    endif
    if (soil_compressibility_index > 0 .and. &
      material_auxvar%setup_reference_pressure) then
      call MaterialReferencePressureSetup(material_auxvar,cell_pressure)
    endif
#endif
    ! creep_closure, fracture, and soil_compressibility are mutually exclusive
    if (option%flow%creep_closure_on) then
      creep_closure => wipp%creep_closure_tables_array( &
                         material_auxvar%creep_closure_id )%ptr
      if ( associated(creep_closure) ) then
        ! option%time here is the t time, not t + dt time.
        creep_closure_time = option%time
        if (option%iflag /= GENERAL_UPDATE_FOR_FIXED_ACCUM) then
          creep_closure_time = creep_closure_time + option%flow_dt
        endif
        
        gen_auxvar%effective_porosity = &
          creep_closure%Evaluate(creep_closure_time,cell_pressure)
      else if (associated(material_auxvar%fracture)) then
          call FracturePoroEvaluate(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dpor_dp)
      else if (soil_compressibility_index > 0) then
          call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dpor_dp)
      endif
    else if (associated(material_auxvar%fracture)) then
      call FracturePoroEvaluate(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dpor_dp)
    else if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dpor_dp)
    endif
    if (option%iflag /= GENERAL_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = gen_auxvar%effective_porosity
    endif
  endif
  if (associated(gen_auxvar%d)) then
    gen_auxvar%d%por_p = dpor_dp
  endif

  ! ALWAYS UPDATE THERMODYNAMIC PROPERTIES FOR BOTH PHASES!!!

  ! Liquid phase thermodynamic properties
  ! must use cell_pressure as the pressure, not %pres(lid)
  if (.not.option%flow%density_depends_on_salinity) then
    if (associated(gen_auxvar%d)) then
      call EOSWaterDensity(gen_auxvar%temp,cell_pressure, &
                           gen_auxvar%den_kg(lid),gen_auxvar%den(lid), &
                           gen_auxvar%d%denl_pl,gen_auxvar%d%denl_T,ierr)
    else
      call EOSWaterDensity(gen_auxvar%temp,cell_pressure, &
                           gen_auxvar%den_kg(lid),gen_auxvar%den(lid),ierr)
    endif
  else
    aux(1) = global_auxvar%m_nacl(1)
    if (associated(gen_auxvar%d)) then
      call EOSWaterDensityExt(gen_auxvar%temp,cell_pressure,aux, &
                              gen_auxvar%den_kg(lid),gen_auxvar%den(lid), &
                              gen_auxvar%d%denl_pl,gen_auxvar%d%denl_T,ierr)
    else
      call EOSWaterDensityExt(gen_auxvar%temp,cell_pressure,aux, &
                              gen_auxvar%den_kg(lid),gen_auxvar%den(lid),ierr)
    endif
  endif
  if (associated(gen_auxvar%d)) then
    call EOSWaterEnthalpy(gen_auxvar%temp,cell_pressure,hw,hw_dp,hw_dT,ierr)
    one_over_dw = 1.d0/gen_auxvar%den(lid)
    !TODO(geh): merge the common terms in dUl_pl and dUl_T equations
    gen_auxvar%d%Ul_pl = hw_dp - &
                         (one_over_dw - &
                          cell_pressure * one_over_dw * one_over_dw * &
                          gen_auxvar%d%denl_pl)
    gen_auxvar%d%Ul_T = hw_dT - &
                        (one_over_dw - &
                         cell_pressure * one_over_dw * one_over_dw * &
                         gen_auxvar%d%denl_T)
    gen_auxvar%d%Hl_pl = hw_dp * 1.d-6
    gen_auxvar%d%Hl_T = hw_dT * 1.d-6
    gen_auxvar%d%Ul_T = gen_auxvar%d%Ul_T * 1.d-6 ! J/kmol-C -> MJ/kmol-C
    gen_auxvar%d%Ul_pl = gen_auxvar%d%Ul_pl * 1.d-6 ! J/kmol-Pa -> MJ/kmol-Pa
  else
    call EOSWaterEnthalpy(gen_auxvar%temp,cell_pressure,hw,ierr)
  endif
  gen_auxvar%H(lid) = hw * 1.d-6 ! J/kmol -> MJ/kmol
  ! MJ/kmol comp
  gen_auxvar%U(lid) = gen_auxvar%H(lid) - &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        (cell_pressure / gen_auxvar%den(lid) * &
                        1.d-6)

  ! Gas phase thermodynamic properties
  ! we cannot use %pres(vpid) as vapor pressure in the liquid phase, since
  ! it can go negative
  if (global_auxvar%istate /= LIQUID_STATE) then
    if (global_auxvar%istate == GAS_STATE) then
      water_vapor_pressure = gen_auxvar%pres(vpid)
    else
      water_vapor_pressure = gen_auxvar%pres(spid)
    endif
    if (associated(gen_auxvar%d)) then
      call EOSGasDensityEnergy(gen_auxvar%temp,gen_auxvar%pres(apid),den_air, &
                               dden_air_dT,dden_air_dpa, &
                               h_air,dh_air_dT,dh_air_dpa, &
                               u_air,du_air_dT,du_air_dpa,ierr)
      ! add in partial w/respec to pa_T
      dden_air_dT = dden_air_dT + dden_air_dpa * dpair_dT
      dden_air_dpg = dden_air_dpa*dpair_dpgas
      ! J/kmol -> MJ/kmol                         
      dh_air_dT = dh_air_dT * 1.d-6
      dh_air_dpa = dh_air_dpa * 1.d-6
      du_air_dT = du_air_dT * 1.d-6
      du_air_dpa = du_air_dpa * 1.d-6
    else
      call EOSGasDensityEnergy(gen_auxvar%temp,gen_auxvar%pres(apid),den_air, &
                               h_air,u_air,ierr)
    endif
    h_air = h_air * 1.d-6 ! J/kmol -> MJ/kmol
    u_air = u_air * 1.d-6 ! J/kmol -> MJ/kmol
    if (associated(gen_auxvar%d)) then
      call EOSWaterSteamDensityEnthalpy(gen_auxvar%temp,water_vapor_pressure, &
                                        den_kg_water_vapor,den_water_vapor, &
                                        h_water_vapor, &
                                        dden_water_vapor_dpv,dden_water_vapor_dT, &
                                        dh_water_vapor_dpv,dh_water_vapor_dT,ierr)
      ! add in partial w/respec to pv_T
      dden_water_vapor_dT = dden_water_vapor_dT + dden_water_vapor_dpv * gen_auxvar%d%pv_T  
      dh_water_vapor_dT = dh_water_vapor_dT + dh_water_vapor_dpv * gen_auxvar%d%pv_T
      !geh: the numerical derivatives with respect to water vapor calculated 
      !     through the chain rule can be very sensitive to the perturbation.
      !     Try decreasing the perturbation to see the effect.
      du_water_vapor_dpv = dh_water_vapor_dpv - &
        (1.d0/den_water_vapor- &
         water_vapor_pressure/(den_water_vapor*den_water_vapor)*dden_water_vapor_dpv)
      du_water_vapor_dT = dh_water_vapor_dT - &
        (gen_auxvar%d%pv_T/den_water_vapor - &
         water_vapor_pressure/(den_water_vapor*den_water_vapor)*dden_water_vapor_dT)
      ! J/kmol -> MJ/kmol                                     
      dh_water_vapor_dpv = dh_water_vapor_dpv * 1.d-6
      dh_water_vapor_dT = dh_water_vapor_dT * 1.d-6
      du_water_vapor_dpv = du_water_vapor_dpv * 1.d-6
      du_water_vapor_dT = du_water_vapor_dT * 1.d-6
    else
      call EOSWaterSteamDensityEnthalpy(gen_auxvar%temp,water_vapor_pressure, &
                                        den_kg_water_vapor,den_water_vapor, &
                                        h_water_vapor,ierr)
    endif
    
    u_water_vapor = h_water_vapor - &
                    ! Pa / kmol/m^3 = J/kmol
                    water_vapor_pressure / den_water_vapor
    h_water_vapor = h_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    u_water_vapor = u_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    gen_auxvar%den(gid) = den_water_vapor + den_air
    gen_auxvar%den_kg(gid) = den_kg_water_vapor + den_air*fmw_comp(gid)
    ! if xmol not set for gas phase, as is the case for LIQUID_STATE, 
    ! set based on densities
!    if (gen_auxvar%xmol(acid,gid) < 1.d-40) then
!      xmol_air_in_gas = den_air / gen_auxvar%den(gid)
!      xmol_water_in_gas = 1.d0 - gen_auxvar%xmol(acid,gid)
!    else    
     xmol_air_in_gas = gen_auxvar%xmol(acid,gid)
     xmol_water_in_gas = gen_auxvar%xmol(wid,gid)
!    endif
      
#ifdef DEBUG_GENERAL  
    xmass_air_in_gas = xmol_air_in_gas*fmw_comp(gid) / &
                       (xmol_water_in_gas*FMWH2O + &
                        xmol_air_in_gas*fmw_comp(gid))
    Hair_J_kg = h_air*1.d6/fmw_comp(gid)
    Uair_J_kg = u_air*1.d6/fmw_comp(gid)
    Hvapor_J_kg = h_water_vapor*1.d6/FMWH2O
    Uvapor_J_kg = u_water_vapor*1.d6/FMWH2O
    Ugas_J_kg = xmass_air_in_gas*Uair_J_kg + &
                (1.d0-xmass_air_in_gas)*Uvapor_J_kg
    Hgas_J_kg = Ugas_J_kg + &
                gen_auxvar%pres(gid)/gen_auxvar%den_kg(gid)
#endif
    
    ! MJ/kmol
    gen_auxvar%U(gid) = xmol_water_in_gas * u_water_vapor + &
                        xmol_air_in_gas * u_air
    Hg_mixture_fractioned = xmol_water_in_gas*h_water_vapor + &
                            xmol_air_in_gas*h_air
    gen_auxvar%H(gid) = gen_auxvar%U(gid) + &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        gen_auxvar%pres(gid)/gen_auxvar%den(gid) * 1.d-6
    
    if (associated(gen_auxvar%d)) then
      gen_auxvar%d%Uv = u_water_vapor
      gen_auxvar%d%Hv = h_water_vapor
      gen_auxvar%d%Uv_pg = du_water_vapor_dpv*gen_auxvar%d%pv_p
      gen_auxvar%d%Uv_pa = du_water_vapor_dpv*gen_auxvar%d%pv_pa
      gen_auxvar%d%Hv_pg = dh_water_vapor_dpv*gen_auxvar%d%pv_p
      gen_auxvar%d%Hv_pa = dh_water_vapor_dpv*gen_auxvar%d%pv_pa
      gen_auxvar%d%Uv_T = du_water_vapor_dT
      gen_auxvar%d%Hv_T = dh_water_vapor_dT
      gen_auxvar%d%Ua = u_air
      gen_auxvar%d%Ha = h_air
      gen_auxvar%d%Ua_pg = du_air_dpa*dpair_dpgas
      gen_auxvar%d%Ua_pa = du_air_dpa
      gen_auxvar%d%Ha_pg = dh_air_dpa*dpair_dpgas
      gen_auxvar%d%Ha_pa = dh_air_dpa
      gen_auxvar%d%Ua_T = du_air_dT
      gen_auxvar%d%Ha_T = dh_air_dT
      gen_auxvar%d%denv = den_water_vapor
      gen_auxvar%d%dena = den_air
      gen_auxvar%d%denv_T = dden_water_vapor_dT
      gen_auxvar%d%dena_T = dden_air_dT
      gen_auxvar%d%denv_pg = dden_water_vapor_dpv * gen_auxvar%d%pv_p 
      gen_auxvar%d%dena_pg = dden_air_dpa * dpair_dpgas
      
      gen_auxvar%d%deng_pg = dden_water_vapor_dpv * gen_auxvar%d%pv_p + &
                             dden_air_dpa * dpair_dpgas
      gen_auxvar%d%deng_pa = dden_water_vapor_dpv * gen_auxvar%d%pv_pa + &
                             dden_air_dpa
      gen_auxvar%d%dengkg_pg = dden_water_vapor_dpv * gen_auxvar%d%pv_p * FMWH2O + &
                               dden_air_dpa * dpair_dpgas * fmw_comp(2)
      gen_auxvar%d%deng_T = dden_water_vapor_dT + dden_air_dT
      gen_auxvar%d%dengkg_T = dden_water_vapor_dT*FMWH2O + dden_air_dT*fmw_comp(2)
      gen_auxvar%d%Ug_pg = xmol_water_in_gas * du_water_vapor_dpv * gen_auxvar%d%pv_p + &
                           gen_auxvar%d%xmol_p(wid,gid) * u_water_vapor + &
                           xmol_air_in_gas * du_air_dpa * dpair_dpgas + &
                           gen_auxvar%d%xmol_p(acid,gid) * u_air
      ! when Ug_pa matters, xmol_pa is in xmol_p(:,lid)
      gen_auxvar%d%Ug_pa = xmol_water_in_gas * du_water_vapor_dpv * gen_auxvar%d%pv_pa + &
                           gen_auxvar%d%xmol_p(wid,lid) * u_water_vapor + &
                           xmol_air_in_gas * du_air_dpa  + &
                           gen_auxvar%d%xmol_p(acid,lid) * u_air
      gen_auxvar%d%Ug_T = xmol_water_in_gas * du_water_vapor_dT + &
                          gen_auxvar%d%xmol_T(wid,gid) * u_water_vapor + &
                          xmol_air_in_gas * du_air_dT + &
                          gen_auxvar%d%xmol_T(acid,gid) * u_air
!geh H = U + pressure / density
      tempreal = 1.d0/gen_auxvar%den(gid)
      gen_auxvar%d%Hg_pg = gen_auxvar%d%Ug_pg + &
        (tempreal-gen_auxvar%pres(gid)*tempreal*tempreal*gen_auxvar%d%deng_pg)*1.d-6
      gen_auxvar%d%Hg_pa = gen_auxvar%d%Ug_pa + &
        (0.d0-gen_auxvar%pres(gid)*tempreal*tempreal*gen_auxvar%d%deng_pa)*1.d-6
      gen_auxvar%d%Hg_T = gen_auxvar%d%Ug_T - &
        gen_auxvar%pres(gid)*tempreal*tempreal*gen_auxvar%d%deng_T*1.d-6
#if 0
!geh U = H - pressure / density
      tempreal = 1.d0/gen_auxvar%den(gid)
      gen_auxvar%d%Ug_pg = gen_auxvar%d%Hg_pg - &
        (tempreal-gen_auxvar%pres(gid)*tempreal*tempreal*gen_auxvar%d%deng_pg)
      gen_auxvar%d%Ug_T = gen_auxvar%d%Hg_T - &
        gen_auxvar%pres(gid)*tempreal*tempreal*gen_auxvar%d%deng_T
#endif
    endif
  endif ! istate /= LIQUID_STATE
  
  if (global_auxvar%istate == LIQUID_STATE .or. &
      global_auxvar%istate == TWO_PHASE_STATE) then
    ! this does not need to be calculated for LIQUID_STATE (=1)
    call characteristic_curves%liq_rel_perm_function% &
           RelativePermeability(gen_auxvar%sat(lid),krl,dkrl_dsatl,option)
    ! dkrl_sat is with respect to liquid pressure, but the primary dependent
    ! variable is gas pressure.  therefore, negate
    dkrl_dsatg = -1.d0 * dkrl_dsatl
    ! use cell_pressure; cell_pressure - psat calculated internally
    if (.not.option%flow%density_depends_on_salinity) then
      if (associated(gen_auxvar%d)) then
        call EOSWaterViscosity(gen_auxvar%temp,cell_pressure, &
                               gen_auxvar%pres(spid), &
                               gen_auxvar%d%psat_T,visl, &
                               dvis_dT,dvis_dp,ierr)
      else
        call EOSWaterViscosity(gen_auxvar%temp,cell_pressure, &
                               gen_auxvar%pres(spid),visl,ierr)
      endif
    else
      aux(1) = global_auxvar%m_nacl(1)
      if (associated(gen_auxvar%d)) then
        call EOSWaterViscosityExt(gen_auxvar%temp,cell_pressure, &
                                  gen_auxvar%pres(spid), &
                                  gen_auxvar%d%psat_T,aux,visl, &
                                  dvis_dT,dvis_dp,ierr)
      else
        call EOSWaterViscosityExt(gen_auxvar%temp,cell_pressure, &
                                  gen_auxvar%pres(spid),aux,visl,ierr)
      endif
    endif
    gen_auxvar%mobility(lid) = krl/visl
    gen_auxvar%kr(lid) = krl
    if (associated(gen_auxvar%d)) then
      ! use chainrule for derivative
      tempreal = -1.d0*krl/(visl*visl)
      gen_auxvar%d%mobilityl_pl = tempreal*dvis_dp
      gen_auxvar%d%mobilityl_T = tempreal*dvis_dT
      gen_auxvar%d%mobilityl_satg = dkrl_dsatg/visl
    endif
  endif

  if (global_auxvar%istate == GAS_STATE .or. &
      global_auxvar%istate == TWO_PHASE_STATE) then
    ! this does not need to be calculated for GAS_STATE (=1)
    call characteristic_curves%gas_rel_perm_function% &
           RelativePermeability(gen_auxvar%sat(lid),krg,dkrg_dsatl,option)                            
    ! dkrl_sat is with respect to liquid pressure, but the primary dependent
    ! variable is gas pressure.  therefore, negate
    dkrg_dsatg = -1.d0 * dkrg_dsatl
    ! STOMP uses separate functions for calculating viscosity of vapor and
    ! and air (WATGSV,AIRGSV) and then uses GASVIS to calculate mixture 
    ! viscosity.
    if (associated(gen_auxvar%d)) then
      call EOSGasViscosity(gen_auxvar%temp, gen_auxvar%pres(apid), &
                           gen_auxvar%pres(gid), den_air, &
    ! viscosity routine says it should be derivative of rho_comp wrt pg
!geh                           dden_air_dT, gen_auxvar%d%deng_pg, &
                           dden_air_dT, dden_air_dpa, dden_air_dpg, &
                           dpair_dT, dpair_dpgas, &
                           visg,dvis_dT,dvis_dpa,dvis_dp,ierr)      
    else    
      call EOSGasViscosity(gen_auxvar%temp,gen_auxvar%pres(apid), &
                           gen_auxvar%pres(gid),den_air,visg,ierr)
    endif
    gen_auxvar%mobility(gid) = krg/visg
    gen_auxvar%kr(gid) = krg
    if (associated(gen_auxvar%d)) then
      ! use chainrule for derivative
      tempreal = -1.d0*krg/(visg*visg)
      gen_auxvar%d%mobilityg_pg = tempreal*dvis_dp
      gen_auxvar%d%mobilityg_T = tempreal*dvis_dT
      gen_auxvar%d%mobilityg_satg = dkrg_dsatg/visg
      gen_auxvar%d%mobilityg_pa = tempreal*dvis_dpa
      gen_auxvar%d%mug = visg
      gen_auxvar%d%mug_T = dvis_dT
      gen_auxvar%d%mug_pg = dvis_dp
    endif    
  endif

#if 0
  if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
    write(*,'(a,i3,2f13.4,es13.6,f13.4,a3)') 'i/l/g/x[a/l]/t: ', &
      natural_id, gen_auxvar%pres(1:2), gen_auxvar%xmol(acid,lid), &
      gen_auxvar%temp, trim(state_char)
    if (natural_id == 5) then
#if 0
    write(*,'(a,i3,8f13.4,a3)') 'i/l/g/a/c/v/s/t: ', &
      natural_id, gen_auxvar%pres(1:6), gen_auxvar%sat(1), gen_auxvar%temp, &
      trim(state_char)
#endif
#if 0
    if (gen_auxvar%sat(2) > 0.d0) then
      write(*,'(a,7es13.6)') 'kmol/kmol/kmol/MJ/MJ/MJ: ', &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(1,1) + &
        gen_auxvar%den(2)*gen_auxvar%sat(2)*gen_auxvar%xmol(1,2),  &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(2,1) + &
        gen_auxvar%den(2)*gen_auxvar%sat(2)*gen_auxvar%xmol(2,2),  &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(1,1) + &
        gen_auxvar%den(2)*gen_auxvar%sat(2)*gen_auxvar%xmol(1,2) + &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(2,1) + &
        gen_auxvar%den(2)*gen_auxvar%sat(2)*gen_auxvar%xmol(2,2),  &
        gen_auxvar%sat(1)*gen_auxvar%den(1)*gen_auxvar%U(1),  &
        gen_auxvar%sat(2)*gen_auxvar%den(2)*gen_auxvar%U(2),  &
        gen_auxvar%sat(1)*gen_auxvar%den(1)*gen_auxvar%U(1) + &
        gen_auxvar%sat(2)*gen_auxvar%den(2)*gen_auxvar%U(2)
    else
      write(*,'(a,7es13.6)') 'kmol/kmol/kmol/MJ/MJ/MJ: ', &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(1,1), &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(2,1), &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(1,1) + &
        gen_auxvar%den(1)*gen_auxvar%sat(1)*gen_auxvar%xmol(2,1), &
        gen_auxvar%sat(1)*gen_auxvar%den(1)*gen_auxvar%U(1), 0.d0, &
        gen_auxvar%sat(1)*gen_auxvar%den(1)*gen_auxvar%U(1)
    endif
#endif
    endif
  endif
#endif

end subroutine GeneralAuxVarCompute

! ************************************************************************** !

subroutine GeneralEOSGasError(natural_id,ierr,gen_auxvar,option)

  !
  !  GeneralEOSGasError: Elaborates when variables exceeds the bounds of
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
  type(general_auxvar_type) :: gen_auxvar
  type(option_type) :: option
  
  
  call PrintMsgByCell(option,natural_id, &
                      'Error in GeneralAuxVarCompute->EOSGasHenry')
  if (ierr == EOS_GAS_TEMP_BOUND_EXCEEDED) then
    option%io_buffer = 'Temperature at cell ID ' // trim(StringWrite(natural_id)) // &
                               ' exceeds the equation of state temperature bound with ' // &
                               trim(StringWrite(gen_auxvar%temp)) // ' [C].'
    call PrintMsgByRank(option)
    general_high_temp_ts_cut = PETSC_TRUE
  endif

  
end subroutine GeneralEOSGasError

! ************************************************************************** !

subroutine GeneralAuxVarUpdateState(x,gen_auxvar,global_auxvar, &
                                    material_auxvar, &
                                    characteristic_curves,natural_id, &
                                    option)
  !
  ! GeneralUpdateState: Updates the state and swaps primary variables
  !
  ! Author: Glenn Hammond
  ! Date: 05/25/11
  ! Modified by Michael Nole
  ! Date: 12/03/18
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use Characteristic_Curves_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  class(characteristic_curves_type) :: characteristic_curves
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscReal, parameter :: epsilon = 0.d0
  PetscReal :: liq_epsilon, gas_epsilon, two_phase_epsilon
  PetscReal :: x(option%nflowdof)
  PetscReal :: Sg_new
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscBool :: istatechng, gas_flag
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: state_change_string


  if (general_immiscible .or. gen_auxvar%istatechng) return

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  gas_flag = PETSC_FALSE

  gen_auxvar%istate_store(PREV_IT) = global_auxvar%istate
  istatechng = PETSC_FALSE

  gas_epsilon = 0.d0
  liq_epsilon = 0.d0
  two_phase_epsilon = 0.d0

  ! Change state
  select case(global_auxvar%istate)

   case(LIQUID_STATE)

      !if (gen_auxvar%pres(vpid) <= gen_auxvar%pres(spid)*(1.d0+ &
      !    window_epsilon) .or. gen_auxvar%pres(apid) >= gen_auxvar% &
      !    pres(lid)*(1.d0-window_epsilon)) then
      if (gen_auxvar%pres(vpid) <= gen_auxvar%pres(spid)*(1.d0- &
          window_epsilon)) then

          global_auxvar%istate = TWO_PHASE_STATE
          liq_epsilon = general_phase_chng_epsilon
          istatechng = PETSC_TRUE

        if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
          write(state_change_string,'(''Liquid -> 2 Phase at Cell '',i8)') &
                natural_id
        else if (option%iflag == GENERAL_UPDATE_FOR_DERIVATIVE) then
          write(state_change_string, &
              '(''Liquid -> 2 Phase at Cell (due to perturbation) '',i8)') &
              natural_id
        else
          write(state_change_string,'(''Liquid -> 2 Phase at Boundary Face '', &
                                    & i8)') natural_id
        endif
      endif

    case(GAS_STATE)
      if (gen_auxvar%pres(vpid) >= gen_auxvar%pres(spid)* &
         (1.d0+window_epsilon)) then

        global_auxvar%istate = TWO_PHASE_STATE
        gas_epsilon = general_phase_chng_epsilon
        gas_flag = PETSC_TRUE
        istatechng = PETSC_TRUE

!#ifdef DEBUG_GENERAL
#ifdef DEBUG_GENERAL_INFO
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,'Before Update',option)
#endif
        if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
          write(state_change_string,'(''Gas -> 2 Phase at Cell '',i8)') &
            natural_id
        else if (option%iflag == GENERAL_UPDATE_FOR_DERIVATIVE) then
          write(state_change_string, &
            '(''Gas -> 2 Phase at Cell (due to perturbation) '',i8)') &
            natural_id
        else
          write(state_change_string,'(''Gas -> 2 Phase at Boundary Face '', &
                                    & i8)') natural_id
        endif
      endif

    case(TWO_PHASE_STATE)
      Sg_new = x(GENERAL_GAS_SATURATION_DOF)
      if (Sg_new < 0.d0) then

        global_auxvar%istate = LIQUID_STATE
        two_phase_epsilon = general_phase_chng_epsilon
        istatechng = PETSC_TRUE

#ifdef DEBUG_GENERAL_INFO
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,'Before Update',option)
#endif
        if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
          write(state_change_string,'(''2 Phase -> Liquid at Cell '',i8)') &
            natural_id
        else if (option%iflag == GENERAL_UPDATE_FOR_DERIVATIVE) then
          write(state_change_string, &
            '(''2 Phase -> Liquid at Cell (due to perturbation) '',i8)') &
            natural_id
        else
          write(state_change_string,'(''2 Phase -> Liquid at Boundary Face '', &
                                    & i8)') natural_id
        endif

      elseif (Sg_new > 1.d0 ) then

        global_auxvar%istate = GAS_STATE
        two_phase_epsilon = general_phase_chng_epsilon
        istatechng = PETSC_TRUE

#ifdef DEBUG_GENERAL_INFO
        call GeneralPrintAuxVars(gen_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,'Before Update',option)
#endif
        if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
          write(state_change_string,'(''2 Phase -> Gas at Cell '',i8)') &
            natural_id
        else if (option%iflag == GENERAL_UPDATE_FOR_DERIVATIVE) then
          write(state_change_string, &
            '(''2 Phase -> Gas at Cell (due to perturbation) '',i8)') &
            natural_id
        else
          write(state_change_string,'(''2 Phase -> Gas at Boundary Face '', &
                                    & i8)') natural_id
        endif
      endif

  end select

  !Update the primary variables
  if (istatechng) then
    general_state_changed = istatechng
    if (general_restrict_state_chng) gen_auxvar%istatechng = PETSC_TRUE

    select case(global_auxvar%istate)
      case(LIQUID_STATE)
        x(GENERAL_LIQUID_PRESSURE_DOF) = gen_auxvar%pres(lid) * (1.d0 - epsilon)
        x(GENERAL_LIQUID_STATE_X_MOLE_DOF) = max(0.d0,gen_auxvar% &
          xmol(acid,lid))*(1.d0 + epsilon)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp*(1.d0-epsilon)
      case(GAS_STATE)
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid) * (1.d0 - epsilon)
        if (general_gas_air_mass_dof == GENERAL_AIR_PRESSURE_INDEX) then
          x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = gen_auxvar%pres(apid)* &
                                                   (1.d0+epsilon)
        else
          x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = max(0.d0,gen_auxvar% &
                                          xmol(wid,gid)*(1.d0 + epsilon))
        endif 
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp*(1.d0-epsilon)
      case(TWO_PHASE_STATE)
        if (gas_flag) then
          x(GENERAL_GAS_SATURATION_DOF) = 1.d0-gas_epsilon
        else
          x(GENERAL_GAS_PRESSURE_DOF) = max(gen_auxvar%pres(gid),gen_auxvar% &
                                 pres(spid))*(1.d0+epsilon)
          x(GENERAL_GAS_SATURATION_DOF) = liq_epsilon
        endif

        if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
          if (.not.general_isothermal) then
            X(GENERAL_ENERGY_DOF) = X(GENERAL_ENERGY_DOF) * &
                                   (1.d0 + epsilon - epsilon)
          endif
        else
          if (gas_flag) then
            X(GENERAL_ENERGY_DOF) = gen_auxvar%pres(apid)*(1.d0-epsilon)
          else
            x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = max(0.01d0* &
                 x(GENERAL_GAS_PRESSURE_DOF), x(GENERAL_GAS_PRESSURE_DOF) - &
                   gen_auxvar%pres(spid))
          endif

          if (x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) <= 0.d0) then
            write(state_change_string,*) natural_id
            option%io_buffer = 'Negative air pressure during state change ' // &
              'at ' // trim(adjustl(state_change_string))
            call PrintMsgByRank(option)
            x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
              0.01d0*x(GENERAL_GAS_PRESSURE_DOF)
          endif
        endif

    end select

    call GeneralAuxVarCompute(x,gen_auxvar, global_auxvar,material_auxvar, &
          characteristic_curves,natural_id,option)
    state_change_string = 'State Transition: ' // trim(state_change_string)
    if (general_print_state_transition) then
      call PrintMsgByRank(option,state_change_string)
    endif
#ifdef DEBUG_GENERAL_INFO
    call GeneralPrintAuxVars(gen_auxvar,global_auxvar,material_auxvar, &
                             natural_id,'After Update',option)
#endif

  endif


end subroutine GeneralAuxVarUpdateState


! ************************************************************************** !

subroutine GeneralAuxVarPerturb(gen_auxvar,global_auxvar, &
                                material_auxvar, &
                                characteristic_curves,natural_id, &
                                option)
  ! 
  ! Calculates auxiliary variables for perturbed system
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(general_auxvar_type) :: gen_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
     
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), &
               pert(option%nflowdof), x_pert_save(option%nflowdof)

  PetscReal :: tempreal
!#define LEGACY_PERTURBATION
!#define HEEHO_PERTURBATION
!#define HP_HARMONIC

#ifdef HEEHO_PERTURBATION
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
#ifdef HP_HARMONIC
  PetscReal, parameter :: min_pres_tol = 0.002d0
  PetscReal, parameter :: min_temp_tol = 4.61d-7
!  PetscReal, parameter :: min_xmol_tol = 2.d-15
!  PetscReal, parameter :: min_airl_tol = 2.d-13
  PetscReal, parameter :: min_xmol_tol = 2.d-11
  PetscReal, parameter :: min_airl_tol = 2.d-11
  PetscReal, parameter :: min_airu_tol = 2.d-8
  PetscReal, parameter :: min_sat_tol = 2.d-10
#else
  PetscReal, parameter :: min_pres_tol = 0.01d0
  PetscReal, parameter :: min_temp_tol = 8.66d-7
  PetscReal, parameter :: min_xmol_tol = 1.d-13
  PetscReal, parameter :: min_airl_tol = 3.16d-11
  PetscReal, parameter :: min_airu_tol = 3.16d-6
  PetscReal, parameter :: min_sat_tol = 3.16d-10
#endif

#else
#ifdef LEGACY_PERTURBATION
  PetscReal, parameter :: perturbation_tolerance = 1.d-5
#else
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
! 1.d-11 works well for Emily's 1D nacl2
!  PetscReal, parameter :: perturbation_tolerance = 1.d-11
#endif
#endif

  PetscReal, parameter :: min_mole_fraction_pert = 1.d-13
  PetscReal, parameter :: min_perturbation = 1.d-10
  PetscInt :: idof

#ifdef DEBUG_GENERAL
  type(global_auxvar_type) :: global_auxvar_debug
  type(general_auxvar_type) :: general_auxvar_debug
  call GlobalAuxVarInit(global_auxvar_debug,option)
  call GeneralAuxVarInit(general_auxvar_debug,PETSC_FALSE,option)
#endif

  select case(global_auxvar%istate)
    case(LIQUID_STATE)
       x(GENERAL_LIQUID_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%liquid_phase)
       x(GENERAL_LIQUID_STATE_X_MOLE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%xmol(option%air_id,option%liquid_phase)
       x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%temp
#ifdef HEEHO_PERTURBATION
       pert(GENERAL_LIQUID_PRESSURE_DOF) = &
            max(perturbation_tolerance*x(GENERAL_LIQUID_PRESSURE_DOF), &
            min_pres_tol)
       tempreal = max(perturbation_tolerance*x(GENERAL_LIQUID_STATE_X_MOLE_DOF), &
            min_xmol_tol)
       if (x(GENERAL_LIQUID_STATE_X_MOLE_DOF) > 1.d-5) then
          pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = -1.d0 * tempreal
       else
          pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = tempreal
       endif
       pert(GENERAL_ENERGY_DOF) = -1.d0 * &
            max(perturbation_tolerance*x(GENERAL_ENERGY_DOF), min_temp_tol)
#else
#ifdef LEGACY_PERTURBATION
       ! if the liquid state, the liquid pressure will always be greater
       ! than zero.
       pert(GENERAL_LIQUID_PRESSURE_DOF) = &
         max(perturbation_tolerance*x(GENERAL_LIQUID_PRESSURE_DOF), &
             perturbation_tolerance)
       ! if the air mole fraction perturbation is too small, the derivatives
       ! can be poor.
       pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = &
         -1.d0*max(perturbation_tolerance*x(GENERAL_LIQUID_STATE_X_MOLE_DOF), &
                   min_mole_fraction_pert)
       pert(GENERAL_ENERGY_DOF) = &
         -1.d0*perturbation_tolerance*x(GENERAL_ENERGY_DOF)
#else
       pert(GENERAL_LIQUID_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_LIQUID_PRESSURE_DOF) + &
         min_perturbation
       if (x(GENERAL_LIQUID_STATE_X_MOLE_DOF) > &
           1.d3 * perturbation_tolerance) then
         pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = -1.d0 * perturbation_tolerance
       else
         pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = perturbation_tolerance
       endif
       pert(GENERAL_ENERGY_DOF) = -1.d0 * &
         (perturbation_tolerance*x(GENERAL_ENERGY_DOF) + min_perturbation)
#endif
#endif
    case(GAS_STATE)
       x(GENERAL_GAS_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
       if (general_gas_air_mass_dof == GENERAL_AIR_PRESSURE_INDEX) then
         x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
           gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
       else
         x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
           gen_auxvar(ZERO_INTEGER)%xmol(option%water_id,option%air_id)
       endif
       x(GENERAL_ENERGY_DOF) = gen_auxvar(ZERO_INTEGER)%temp
#ifdef HEEHO_PERTURBATION
       pert(GENERAL_GAS_PRESSURE_DOF) = &
            max(perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF), &
            min_pres_tol)
       if (x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) > 1.d0) then
          tempreal = max(perturbation_tolerance*x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF), &
               min_airu_tol)
       else
          tempreal = max(perturbation_tolerance*x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF), &
               min_airl_tol)
       endif
       if (x(GENERAL_GAS_PRESSURE_DOF) - &
            x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) > tempreal) then
          pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = tempreal
       else
          pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
       endif
       pert(GENERAL_ENERGY_DOF) = max(perturbation_tolerance*x(GENERAL_ENERGY_DOF), &
            min_temp_tol)
#else
#ifdef LEGACY_PERTURBATION
       ! gas pressure [p(g)] must always be perturbed down as p(v) = p(g) - p(a)
       ! and p(v) >= Psat (i.e. an increase in p(v)) results in two phase.
       !pert(GENERAL_GAS_PRESSURE_DOF) = &
       !  -1.d0*perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)
       
       ! perturb air pressure towards gas pressure unless the perturbed
       ! air pressure exceeds the gas pressure
       if (general_gas_air_mass_dof == GENERAL_AIR_PRESSURE_INDEX) then
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
             x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) > &
             perturbation_tolerance* &
             x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)) then
           pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
             perturbation_tolerance*x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
         else
           pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
             -1.d0*perturbation_tolerance*x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
         endif
       else
         if (x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) > &
            1.d3 * perturbation_tolerance) then
            pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = -1.d0 * & 
                 perturbation_tolerance * x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
         else
            pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = perturbation_tolerance*&
                 x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
         endif
       endif
       pert(GENERAL_ENERGY_DOF) = &
         perturbation_tolerance*x(GENERAL_ENERGY_DOF)
#else
       ! gas pressure [p(g)] must always be perturbed down as p(v) = p(g) - p(a)
       ! and p(v) >= Psat (i.e. an increase in p(v)) results in two phase.
       !pert(GENERAL_GAS_PRESSURE_DOF) = -1.d0 * &
       !  (perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF) + min_perturbation)

       !MAN: Try perturbing upward, because lower gas pressure is associated with
       !     lower gas saturation (i.e. two-phase)
       pert(GENERAL_GAS_PRESSURE_DOF) = perturbation_tolerance* &
                                x(GENERAL_GAS_PRESSURE_DOF) + min_perturbation

       ! perturb air pressure towards gas pressure unless the perturbed
       ! air pressure exceeds the gas pressure
       tempreal = perturbation_tolerance* &
                  x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) + min_perturbation
       if (general_gas_air_mass_dof == GENERAL_AIR_PRESSURE_INDEX) then
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
             x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) > tempreal) then
           pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = tempreal
         else
           pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
         endif
       else
         if (x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) > &
            1.d3 * perturbation_tolerance) then
            pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = -1.d0 * &
                 perturbation_tolerance * x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
         else
            pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = perturbation_tolerance*&
                 x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
         endif
       endif
       
       pert(GENERAL_ENERGY_DOF) = &
         perturbation_tolerance*x(GENERAL_ENERGY_DOF) + min_perturbation
#endif
#endif
    case(TWO_PHASE_STATE)
       x(GENERAL_GAS_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
!       x(GENERAL_AIR_PRESSURE_DOF) = &
!         gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
       x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(option%gas_phase)
#ifdef HEEHO_PERTURBATION
       if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
          x(GENERAL_ENERGY_DOF) = gen_auxvar(ZERO_INTEGER)%temp
          pert(GENERAL_ENERGY_DOF) = max(perturbation_tolerance*x(GENERAL_ENERGY_DOF), &
               min_temp_tol)
       else
         x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
              gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
         if (x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) > 1.d0) then
            tempreal = max(perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF), &
                 min_airu_tol)
         else
            tempreal = max(perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF), &
                 min_airl_tol)
         endif
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
              x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) > tempreal) then
            pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = tempreal
         else
            pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
         endif
      endif
      
      pert(GENERAL_GAS_PRESSURE_DOF) = &
           max(perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF), &
           min_pres_tol)
      tempreal = max(perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF), &
                 min_sat_tol)
      if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then 
         pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * tempreal
      else
         pert(GENERAL_GAS_SATURATION_DOF) = tempreal
      endif
#else
#ifdef LEGACY_PERTURBATION
       if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%temp
         pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)
       else
         ! here GENERAL_2PH_STATE_AIR_PRESSURE_DOF = GENERAL_ENERGY_DOF
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
         ! perturb air pressure towards gas pressure unless the perturbed
         ! air pressure exceeds the gas pressure
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
             x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) > &
             perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)) then
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
             perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)
         else
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
             -1.d0*perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)
         endif
       endif
       pert(GENERAL_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)
       ! always perturb toward 0.5
       if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then 
         pert(GENERAL_GAS_SATURATION_DOF) = &
           -1.d0*(perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF)      
!         I think it should be what's below, has to confirm with M. Nole. -hdp
!           -1.d0*(perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF) + &
!           min_perturbation)
       else
         pert(GENERAL_GAS_SATURATION_DOF) = &
           perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF)  
!         I think it should be what's below, has to confirm with M. Nole. -hdp
!           perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF) + min_perturbation
       endif
#else
       if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%temp
         if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
           pert(GENERAL_ENERGY_DOF) = -1.d0 * &
             (perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation)
         else 
           pert(GENERAL_ENERGY_DOF) = &
             perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation
         endif
       else
         ! here GENERAL_2PH_STATE_AIR_PRESSURE_DOF = GENERAL_ENERGY_DOF
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
         ! perturb air pressure towards gas pressure unless the perturbed
         ! air pressure exceeds the gas pressure
         tempreal = perturbation_tolerance* &
                    x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) + min_perturbation
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
             x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) > tempreal) then
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = tempreal
         else
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
         endif
       endif
       if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then 
!         I think it should be what's below, has to confirm with M. Nole. -hdp
         pert(GENERAL_GAS_PRESSURE_DOF) = -1.d0 *(perturbation_tolerance* &
                      x(GENERAL_GAS_PRESSURE_DOF)+min_perturbation)
         pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * (perturbation_tolerance * &
                                x(GENERAL_GAS_SATURATION_DOF)+min_perturbation)
       else
!         I think it should be what's below, has to confirm with M. Nole. -hdp
         pert(GENERAL_GAS_PRESSURE_DOF) = perturbation_tolerance* &
                      x(GENERAL_GAS_PRESSURE_DOF)+min_perturbation
         pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance * &
                        x(GENERAL_GAS_SATURATION_DOF)+min_perturbation
       endif
#endif
#endif
  end select
  
  ! GENERAL_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = GENERAL_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    gen_auxvar(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    x_pert_save = x_pert
    call GeneralAuxVarCompute(x_pert,gen_auxvar(idof),global_auxvar, &
                              material_auxvar, &
                              characteristic_curves,natural_id,option)
#ifdef DEBUG_GENERAL
    call GlobalAuxVarCopy(global_auxvar,global_auxvar_debug,option)
    call GeneralAuxVarCopy(gen_auxvar(idof),general_auxvar_debug,option)
    call GeneralAuxVarUpdateState(x_pert,general_auxvar_debug, &
                                  global_auxvar_debug, &
                                  material_auxvar, &
                                  characteristic_curves, &
                                  natural_id,option)
    if (global_auxvar%istate /= global_auxvar_debug%istate) then
      write(option%io_buffer, &
            &'(''Change in state due to perturbation: '',i3,'' -> '',i3, &
            &'' at cell '',i3,'' for dof '',i3)') &
        global_auxvar%istate, global_auxvar_debug%istate, natural_id, idof
      call PrintMsg(option)
      write(option%io_buffer,'(''orig: '',6es17.8)') x(1:3)
      call PrintMsg(option)
      write(option%io_buffer,'(''pert: '',6es17.8)') x_pert_save(1:3)
      call PrintMsg(option)
    endif
#endif

  enddo

  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      gen_auxvar(GENERAL_LIQUID_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_LIQUID_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(TWO_PHASE_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
      if (general_2ph_energy_dof == GENERAL_AIR_PRESSURE_INDEX) then
        gen_auxvar(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)%pert = &
          gen_auxvar(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)%pert / &
          GENERAL_PRESSURE_SCALE
      endif
    case(GAS_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
      gen_auxvar(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)%pert / &
        GENERAL_PRESSURE_SCALE
  end select
  
#ifdef DEBUG_GENERAL
  call GlobalAuxVarStrip(global_auxvar_debug)
  call GeneralAuxVarStrip(general_auxvar_debug)
#endif
 
end subroutine GeneralAuxVarPerturb

! ************************************************************************** !

subroutine GeneralPrintAuxVars(general_auxvar,global_auxvar,material_auxvar, &
                               natural_id,string,option)
  ! 
  ! Prints out the contents of an auxvar
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/13
  ! 

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(general_auxvar_type) :: general_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
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
    case(LIQUID_STATE)
      print *, '     Thermodynamic state: Liquid phase'
      liquid_density = general_auxvar%den(lid)
      liquid_energy = general_auxvar%U(lid)
      liquid_saturation = general_auxvar%sat(lid)
    case(GAS_STATE)
      print *, '     Thermodynamic state: Gas phase'
      gas_density = general_auxvar%den(gid)
      gas_energy = general_auxvar%U(gid)
      gas_saturation = general_auxvar%sat(gid)
    case(TWO_PHASE_STATE)
      print *, '     Thermodynamic state: Two phase'
      liquid_density = general_auxvar%den(lid)
      gas_density = general_auxvar%den(gid)
      liquid_energy = general_auxvar%U(lid)
      gas_energy = general_auxvar%U(gid)
      liquid_saturation = general_auxvar%sat(lid)
      gas_saturation = general_auxvar%sat(gid)
  end select
  liquid_mass = (liquid_density*general_auxvar%xmol(lid,lid)* & 
                 liquid_saturation+ &
                 gas_density*general_auxvar%xmol(lid,gid)* & 
                 gas_saturation)* & 
                 general_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*general_auxvar%xmol(gid,lid)* & 
              liquid_saturation+ &
              gas_density*general_auxvar%xmol(gid,gid)* & 
              gas_saturation)* & 
              general_auxvar%effective_porosity*material_auxvar%volume
  print *, 'tot liq comp mass [kmol]: ', liquid_mass
  print *, 'tot gas comp mass [kmol]: ', gas_mass
  print *, '             energy [MJ]: ', liquid_mass*liquid_energy + &
                                         gas_mass*gas_energy
  print *, '         liquid pressure: ', general_auxvar%pres(lid)
  print *, '            gas pressure: ', general_auxvar%pres(gid)
  print *, '            air pressure: ', general_auxvar%pres(apid)
  print *, '      capillary pressure: ', general_auxvar%pres(cpid)
  print *, '          vapor pressure: ', general_auxvar%pres(vpid)
  print *, '     saturation pressure: ', general_auxvar%pres(spid)
  print *, '       liquid saturation: ', general_auxvar%sat(lid)
  print *, '          gas saturation: ', general_auxvar%sat(gid)
  print *, '   liquid density [kmol]: ', general_auxvar%den(lid)
  print *, '      gas density [kmol]: ', general_auxvar%den(gid)
  print *, '     liquid density [kg]: ', general_auxvar%den_kg(lid)
  print *, '        gas density [kg]: ', general_auxvar%den_kg(gid)
  print *, '         temperature [C]: ', general_auxvar%temp
  print *, '      liquid H [MJ/kmol]: ', general_auxvar%H(lid)
  print *, '         gas H [MJ/kmol]: ', general_auxvar%H(gid)
  print *, '      liquid U [MJ/kmol]: ', general_auxvar%U(lid)
  print *, '         gas U [MJ/kmol]: ', general_auxvar%U(gid)
  print *, '     X (water in liquid): ', general_auxvar%xmol(lid,lid)
  print *, '       X (air in liquid): ', general_auxvar%xmol(gid,lid)
  print *, '        X (water in gas): ', general_auxvar%xmol(lid,gid)
  print *, '          X (air in gas): ', general_auxvar%xmol(gid,gid)
  print *, '         liquid mobility: ', general_auxvar%mobility(lid)
  print *, '            gas mobility: ', general_auxvar%mobility(gid)
  print *, '      effective porosity: ', general_auxvar%effective_porosity
  print *, '--------------------------------------------------------'

end subroutine GeneralPrintAuxVars

! ************************************************************************** !

subroutine GeneralOutputAuxVars1(general_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,string,append,option)
  ! 
  ! Prints out the contents of an auxvar to a file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/13
  ! 

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(general_auxvar_type) :: general_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
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
    case(LIQUID_STATE)
      write(86,*) ' Thermodynamic state: Liquid phase'
      liquid_density = general_auxvar%den(lid)
      liquid_energy = general_auxvar%U(lid)
      liquid_saturation = general_auxvar%sat(lid)
    case(GAS_STATE)
      write(86,*) ' Thermodynamic state: Gas phase'
      gas_density = general_auxvar%den(gid)
      gas_energy = general_auxvar%U(gid)
      gas_saturation = general_auxvar%sat(gid)
    case(TWO_PHASE_STATE)
      write(86,*) ' Thermodynamic state: Two phase'
      liquid_density = general_auxvar%den(lid)
      gas_density = general_auxvar%den(gid)
      liquid_energy = general_auxvar%U(lid)
      gas_energy = general_auxvar%U(gid)
      liquid_saturation = general_auxvar%sat(lid)
      gas_saturation = general_auxvar%sat(gid)
  end select
  liquid_mass = (liquid_density*general_auxvar%xmol(lid,lid)* & 
                 liquid_saturation+ &
                 gas_density*general_auxvar%xmol(lid,gid)* & 
                 gas_saturation)* & 
                 general_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*general_auxvar%xmol(gid,lid)* & 
              liquid_saturation+ &
              gas_density*general_auxvar%xmol(gid,gid)* & 
              gas_saturation)* & 
              general_auxvar%effective_porosity*material_auxvar%volume
  write(86,*) 'tot liq comp mass [kmol]: ', liquid_mass
  write(86,*) 'tot gas comp mass [kmol]: ', gas_mass
  write(86,*) '             energy [MJ]: ', liquid_mass*liquid_energy + &
                                            gas_mass*gas_energy
  write(86,*) '         liquid pressure: ', general_auxvar%pres(lid)
  write(86,*) '            gas pressure: ', general_auxvar%pres(gid)
  write(86,*) '            air pressure: ', general_auxvar%pres(apid)
  write(86,*) '      capillary pressure: ', general_auxvar%pres(cpid)
  write(86,*) '          vapor pressure: ', general_auxvar%pres(vpid)
  write(86,*) '     saturation pressure: ', general_auxvar%pres(spid)
  write(86,*) '         temperature [C]: ', general_auxvar%temp
  write(86,*) '       liquid saturation: ', general_auxvar%sat(lid)
  write(86,*) '          gas saturation: ', general_auxvar%sat(gid)
  write(86,*) '   liquid density [kmol]: ', general_auxvar%den(lid)
  write(86,*) '     liquid density [kg]: ', general_auxvar%den_kg(lid)
  write(86,*) '      gas density [kmol]: ', general_auxvar%den(gid)
  write(86,*) '        gas density [kg]: ', general_auxvar%den_kg(gid)
  write(86,*) '     X (water in liquid): ', general_auxvar%xmol(lid,lid)
  write(86,*) '       X (air in liquid): ', general_auxvar%xmol(gid,lid)
  write(86,*) '        X (water in gas): ', general_auxvar%xmol(lid,gid)
  write(86,*) '          X (air in gas): ', general_auxvar%xmol(gid,gid)
  write(86,*) '      liquid H [MJ/kmol]: ', general_auxvar%H(lid)
  write(86,*) '         gas H [MJ/kmol]: ', general_auxvar%H(gid)
  write(86,*) '      liquid U [MJ/kmol]: ', general_auxvar%U(lid)
  write(86,*) '         gas U [MJ/kmol]: ', general_auxvar%U(gid)
  write(86,*) '         liquid mobility: ', general_auxvar%mobility(lid)
  write(86,*) '            gas mobility: ', general_auxvar%mobility(gid)
  write(86,*) '      effective porosity: ', general_auxvar%effective_porosity
  write(86,*) '...'
  write(86,*) liquid_mass
  write(86,*) gas_mass
  write(86,*) liquid_mass*general_auxvar%U(lid) + &
              gas_mass*general_auxvar%U(gid)
  write(86,*) general_auxvar%pres(lid)
  write(86,*) general_auxvar%pres(gid)
  write(86,*) general_auxvar%pres(apid)
  write(86,*) general_auxvar%pres(cpid)
  write(86,*) general_auxvar%pres(vpid)
  write(86,*) general_auxvar%pres(spid)
  write(86,*) general_auxvar%temp
  write(86,*) general_auxvar%sat(lid)
  write(86,*) general_auxvar%sat(gid)
  write(86,*) general_auxvar%den(lid)
  write(86,*) general_auxvar%den_kg(lid)
  write(86,*) general_auxvar%den(gid)
  write(86,*) general_auxvar%den_kg(gid)
  write(86,*) general_auxvar%xmol(lid,lid)
  write(86,*) general_auxvar%xmol(gid,lid)
  write(86,*) general_auxvar%xmol(lid,gid)
  write(86,*) general_auxvar%xmol(gid,gid)
  write(86,*) general_auxvar%H(lid)
  write(86,*) general_auxvar%H(gid)
  write(86,*) general_auxvar%U(lid)
  write(86,*) general_auxvar%U(gid)
  write(86,*) ''
  write(86,*) general_auxvar%mobility(lid)
  write(86,*) general_auxvar%mobility(gid)
  write(86,*) general_auxvar%effective_porosity
  write(86,*) '--------------------------------------------------------'
  
  close(86)

end subroutine GeneralOutputAuxVars1

! ************************************************************************** !

subroutine GeneralOutputAuxVars2(general_auxvars,global_auxvars,option)
  ! 
  ! Prints out the contents of an auxvar to a file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/13
  ! 

  use Global_Aux_module
  use Option_module

  implicit none

  type(general_auxvar_type) :: general_auxvars(0:,:)
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
  
  string = 'general_auxvar.txt'
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
    ((general_auxvars(idof,i)%pres(lid),i=1,n),idof=0,3)
  write(86,100) '         gas pressure: ', &
    ((general_auxvars(idof,i)%pres(gid),i=1,n),idof=0,3)
  write(86,100) '         air pressure: ', &
    ((general_auxvars(idof,i)%pres(apid),i=1,n),idof=0,3)
  write(86,100) '   capillary pressure: ', &
    ((general_auxvars(idof,i)%pres(cpid),i=1,n),idof=0,3)
  write(86,100) '       vapor pressure: ', &
    ((general_auxvars(idof,i)%pres(vpid),i=1,n),idof=0,3)
  write(86,100) '      temperature [C]: ', &
    ((general_auxvars(idof,i)%temp,i=1,n),idof=0,3)
  write(86,100) '    liquid saturation: ', &
    ((general_auxvars(idof,i)%sat(lid),i=1,n),idof=0,3)
  write(86,100) '       gas saturation: ', &
    ((general_auxvars(idof,i)%sat(gid),i=1,n),idof=0,3)
  write(86,100) 'liquid density [kmol]: ', &
    ((general_auxvars(idof,i)%den(lid),i=1,n),idof=0,3)
  write(86,100) '  liquid density [kg]: ', &
    ((general_auxvars(idof,i)%den_kg(lid),i=1,n),idof=0,3)
  write(86,100) '   gas density [kmol]: ', &
    ((general_auxvars(idof,i)%den(gid),i=1,n),idof=0,3)
  write(86,100) '     gas density [kg]: ', &
    ((general_auxvars(idof,i)%den_kg(gid),i=1,n),idof=0,3)
  write(86,100) '  X (water in liquid): ', &
    ((general_auxvars(idof,i)%xmol(lid,lid),i=1,n),idof=0,3)
  write(86,100) '    X (air in liquid): ', &
    ((general_auxvars(idof,i)%xmol(gid,lid),i=1,n),idof=0,3)
  write(86,100) '     X (water in gas): ', &
    ((general_auxvars(idof,i)%xmol(lid,gid),i=1,n),idof=0,3)
  write(86,100) '       X (air in gas): ', &
    ((general_auxvars(idof,i)%xmol(gid,gid),i=1,n),idof=0,3)
  write(86,100) '   liquid H [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%H(lid),i=1,n),idof=0,3)
  write(86,100) '      gas H [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%H(gid),i=1,n),idof=0,3)
  write(86,100) '   liquid U [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%U(lid),i=1,n),idof=0,3)
  write(86,100) '      gas U [MJ/kmol]: ', &
    ((general_auxvars(idof,i)%U(gid),i=1,n),idof=0,3)
  write(86,*)
  write(86,100) '      liquid mobility: ', &
    ((general_auxvars(idof,i)%mobility(lid),i=1,n),idof=0,3)
  write(86,100) '         gas mobility: ', &
    ((general_auxvars(idof,i)%mobility(gid),i=1,n),idof=0,3)
  write(86,100) '   effective porosity: ', &
    ((general_auxvars(idof,i)%effective_porosity,i=1,n),idof=0,3)
  
  close(86)

end subroutine GeneralOutputAuxVars2

! ************************************************************************** !

subroutine GeneralSetBogusDeriv(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/29/19
  ! 

  implicit none

  type(general_auxvar_type) :: auxvar

  auxvar%d%pc_satg = 1234.56789d0
  auxvar%d%por_p = 1234.56789d0
  auxvar%d%denl_pl = 1234.56789d0
  auxvar%d%denl_T = 1234.56789d0
  auxvar%d%deng_pg = 1234.56789d0
  auxvar%d%deng_pa = 1234.56789d0
  auxvar%d%deng_T = 1234.56789d0
  auxvar%d%dengkg_pg = 1234.56789d0
  auxvar%d%dengkg_T = 1234.56789d0
  auxvar%d%Ul_pl = 1234.56789d0
  auxvar%d%Ul_T = 1234.56789d0
  auxvar%d%Hl_pl = 1234.56789d0
  auxvar%d%Hl_T = 1234.56789d0
  auxvar%d%Ug_pg = 1234.56789d0
  auxvar%d%Ug_pa = 1234.56789d0
  auxvar%d%Ug_T = 1234.56789d0
  auxvar%d%Hg_pg = 1234.56789d0
  auxvar%d%Hg_pa = 1234.56789d0
  auxvar%d%Hg_T = 1234.56789d0

  auxvar%d%Hv = 1234.56789d0
  auxvar%d%Ha = 1234.56789d0
  auxvar%d%Uv = 1234.56789d0
  auxvar%d%Ua = 1234.56789d0
  auxvar%d%Hv_pg = 1234.56789d0
  auxvar%d%Hv_pa = 1234.56789d0
  auxvar%d%Hv_T = 1234.56789d0
  auxvar%d%Ha_pg = 1234.56789d0
  auxvar%d%Ha_pa = 1234.56789d0
  auxvar%d%Ha_T = 1234.56789d0
  auxvar%d%Uv_pg = 1234.56789d0
  auxvar%d%Uv_pa = 1234.56789d0
  auxvar%d%Uv_T = 1234.56789d0
  auxvar%d%Ua_pg = 1234.56789d0
  auxvar%d%Ua_pa = 1234.56789d0
  auxvar%d%Ua_T = 1234.56789d0
  auxvar%d%denv = 1234.56789d0
  auxvar%d%dena = 1234.56789d0
  auxvar%d%denv_T = 1234.56789d0
  auxvar%d%dena_T = 1234.56789d0
  auxvar%d%denv_pg = 1234.56789d0
  auxvar%d%dena_pg = 1234.56789d0
  auxvar%d%Hc = 1234.56789d0
    
  auxvar%d%psat_p = 1234.56789d0
  auxvar%d%psat_T = 1234.56789d0
  auxvar%d%pv_p = 1234.56789d0
  auxvar%d%pv_pa = 1234.56789d0
  auxvar%d%pv_T = 1234.56789d0
  auxvar%d%Hc_p = 1234.56789d0
  auxvar%d%Hc_T = 1234.56789d0
  auxvar%d%mobilityl_pl = 1234.56789d0
  auxvar%d%mobilityl_T = 1234.56789d0
  auxvar%d%mobilityl_satg = 1234.56789d0
  auxvar%d%mobilityg_pg = 1234.56789d0
  auxvar%d%mobilityg_T = 1234.56789d0
  auxvar%d%mobilityg_satg = 1234.56789d0
  auxvar%d%mobilityg_pa = 1234.56789d0
  auxvar%d%mug = 1234.56789d0
  auxvar%d%mug_T = 1234.56789d0
  auxvar%d%mug_pg = 1234.56789d0
  auxvar%d%xmol_p(:,:) = 1234.56789d0
  auxvar%d%xmol_T(:,:) = 1234.56789d0

end subroutine GeneralSetBogusDeriv

! ************************************************************************** !

subroutine GeneralZeroAnalyticalDeriv(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/29/19
  ! 

  implicit none

  type(general_auxvar_type) :: auxvar

  auxvar%d%pc_satg = 0.d0
  auxvar%d%por_p = 0.d0
  auxvar%d%denl_pl = 0.d0
  auxvar%d%denl_T = 0.d0
  auxvar%d%deng_pg = 0.d0
  auxvar%d%deng_pa = 0.d0
  auxvar%d%deng_T = 0.d0
  auxvar%d%dengkg_pg = 0.d0
  auxvar%d%dengkg_T = 0.d0
  auxvar%d%Ul_pl = 0.d0
  auxvar%d%Ul_T = 0.d0
  auxvar%d%Hl_pl = 0.d0
  auxvar%d%Hl_T = 0.d0
  auxvar%d%Ug_pg = 0.d0
  auxvar%d%Ug_pa = 0.d0
  auxvar%d%Ug_T = 0.d0
  auxvar%d%Hg_pg = 0.d0
  auxvar%d%Hg_pa = 0.d0
  auxvar%d%Hg_T = 0.d0

  auxvar%d%Hv = 0.d0
  auxvar%d%Ha = 0.d0
  auxvar%d%Uv = 0.d0
  auxvar%d%Ua = 0.d0
  auxvar%d%Hv_pg = 0.d0
  auxvar%d%Hv_pa = 0.d0
  auxvar%d%Hv_T = 0.d0
  auxvar%d%Ha_pg = 0.d0
  auxvar%d%Ha_pa = 0.d0
  auxvar%d%Ha_T = 0.d0
  auxvar%d%Uv_pg = 0.d0
  auxvar%d%Uv_pa = 0.d0
  auxvar%d%Uv_T = 0.d0
  auxvar%d%Ua_pg = 0.d0
  auxvar%d%Ua_pa = 0.d0
  auxvar%d%Ua_T = 0.d0
  auxvar%d%denv = 0.d0
  auxvar%d%dena = 0.d0
  auxvar%d%denv_T = 0.d0
  auxvar%d%dena_T = 0.d0
  auxvar%d%denv_pg = 0.d0
  auxvar%d%dena_pg = 0.d0
  auxvar%d%Hc = 0.d0
    
  auxvar%d%psat_p = 0.d0
  auxvar%d%psat_T = 0.d0
  auxvar%d%pv_p = 0.d0
  auxvar%d%pv_pa = 0.d0
  auxvar%d%pv_T = 0.d0
  auxvar%d%Hc_p = 0.d0
  auxvar%d%Hc_T = 0.d0
  auxvar%d%mobilityl_pl = 0.d0
  auxvar%d%mobilityl_T = 0.d0
  auxvar%d%mobilityl_satg = 0.d0
  auxvar%d%mobilityg_pg = 0.d0
  auxvar%d%mobilityg_T = 0.d0
  auxvar%d%mobilityg_satg = 0.d0
  auxvar%d%mobilityg_pa = 0.d0
  auxvar%d%mug = 0.d0
  auxvar%d%mug_T = 0.d0
  auxvar%d%mug_pg = 0.d0
  auxvar%d%xmol_p(:,:) = 0.d0
  auxvar%d%xmol_T(:,:) = 0.d0

end subroutine GeneralZeroAnalyticalDeriv

! ************************************************************************** !

subroutine GeneralAuxVarSingleDestroy(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(general_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call GeneralAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)  

end subroutine GeneralAuxVarSingleDestroy

! ************************************************************************** !

subroutine GeneralAuxVarArray1Destroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(general_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call GeneralAuxVarStrip(auxvars(iaux))
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)  

end subroutine GeneralAuxVarArray1Destroy

! ************************************************************************** !

subroutine GeneralAuxVarArray2Destroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(general_auxvar_type), pointer :: auxvars(:,:)
  
  PetscInt :: iaux, idof
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call GeneralAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)  

end subroutine GeneralAuxVarArray2Destroy

! ************************************************************************** !

subroutine GeneralAuxVarStrip(auxvar)
  ! 
  ! GeneralAuxVarDestroy: Deallocates a general auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  type(general_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%pres)  
  call DeallocateArray(auxvar%sat)  
  call DeallocateArray(auxvar%den)  
  call DeallocateArray(auxvar%den_kg)  
  call DeallocateArray(auxvar%xmol)  
  call DeallocateArray(auxvar%H)  
  call DeallocateArray(auxvar%U)  
  call DeallocateArray(auxvar%mobility)  
  if (associated(auxvar%d)) then
    deallocate(auxvar%d)
    nullify(auxvar%d)
  endif
  
end subroutine GeneralAuxVarStrip

! ************************************************************************** !

subroutine GeneralAuxDestroy(aux)
  ! 
  ! Deallocates a general auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/11
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  type(general_type), pointer :: aux
  
  if (.not.associated(aux)) return
  
  call GeneralAuxVarDestroy(aux%auxvars)
  call GeneralAuxVarDestroy(aux%auxvars_bc)
  call GeneralAuxVarDestroy(aux%auxvars_ss)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%general_parameter)) then
    call DeallocateArray(aux%general_parameter%diffusion_coefficient)
    deallocate(aux%general_parameter)
  endif
  nullify(aux%general_parameter)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine GeneralAuxDestroy

end module General_Aux_module
