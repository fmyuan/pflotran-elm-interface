module MpFlow_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use AuxVars_Flow_Energy_module

  implicit none

  private 

  PetscReal, public :: flow_itol_scaled_res = 1.d-15
  PetscReal, public :: flow_itol_rel_update = UNINITIALIZED_DOUBLE
  !PetscBool, public :: flow_using_newtontr = PETSC_FALSE

  ! currently only liq. water and air gas as fluids (option%nfluids = 2)
  PetscReal,public :: FMW_FLUIDS(2) = [FMWH2O,FMWAIR]
  ! currently only liq. water as flow species
  PetscInt, public :: IFLOW1 = ONE_INTEGER       ! liq. water
  PetscInt, public :: IFLOW2 = TWO_INTEGER       ! air (not-yet, AND TODO)

  ! 2 dofs currently
  PetscInt, public :: IFDOF1 = MPFLOW_PRESSURE_DOF
  PetscInt, public :: IFDOF2 = MPFLOW_TEMPERATURE_DOF
  PetscInt, public :: IFDOF3 = MPFLOW_CONDUCTANCE_DOF ! not-yet
  PetscInt, public :: IFDOF4 = MPFLOW_ENTHALPY_DOF    ! not-yet

  type, public, extends(auxvar_flow_energy_type) :: flow_auxvar_type

    PetscReal :: transient_por
    PetscReal :: air_pressure    ! unit: Pa ( air pressure of water-air interface. It's reference_pressure if an open system)

    PetscReal :: molv_air                 ! vapor mole fraction in air mixture
    PetscReal, pointer :: dmolv_air(:)    ! (nflowdof)

    PetscReal :: pres_fh2o                ! liq. water pressure on ice-included soil matrix
    PetscReal, pointer :: dpres_fh2o(:)   ! (nflowdof)

    PetscReal :: vis
    PetscReal :: kvr
    PetscReal, pointer :: dkvr(:)         ! (nflowdof)

    PetscReal :: Dk_eff
    PetscReal, pointer :: dDk_eff(:)      ! (nflowdof)

  end type flow_auxvar_type


  type, public :: flow_parameter_type
    PetscReal :: dencpr   ! Rock (soil particle) densityXspecific_heat_capacity: MJ/m^3-K
    PetscReal :: ckdry    ! Thermal conductivity (dry)
    PetscReal :: ckwet    ! Thermal conductivity (saturated - wet)
    PetscReal :: alpha
    PetscReal :: ckfrozen ! Thermal conductivity (frozen/saturated soil)
    PetscReal :: alpha_fr ! exponent frozen
    PetscReal, pointer :: sir(:)   ! (nphase)
    PetscReal, pointer :: diffusion_coefficient(:)         !(nphase)
    PetscReal, pointer :: diffusion_activation_energy(:)   !(nphase)
  end type flow_parameter_type
  
  type, public :: MpFlow_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(flow_parameter_type), pointer :: flow_parameters(:)    ! material_id
    type(flow_auxvar_type), pointer :: auxvars(:)               ! grid_id
  end type MpFlow_type

  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

  public :: MpFlowAuxCreate, MpFlowAuxDestroy, &
            MpFlowAuxVarCompute, MpFlowAuxVarInit, &
            MpFlowAuxVarCopy, MpFlowAuxVarDestroy

  public :: DarcyFlowDerivative, &
            AdvectionDerivative, &
            ConductionDerivative

contains

! ************************************************************************** !

function MpFlowAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(MpFlow_type), pointer :: MpFlowAuxCreate
  
  type(MpFlow_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  aux%n_zero_rows = 0

  nullify(aux%flow_parameters)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  MpFlowAuxCreate => aux
  
end function MpFlowAuxCreate

! ************************************************************************** !

subroutine MpFlowAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module
  use PFLOTRAN_Constants_module, only : UNINITIALIZED_DOUBLE

  implicit none
  
  type(flow_auxvar_type) :: auxvar
  type(option_type) :: option
  
  auxvar%transient_por = UNINITIALIZED_DOUBLE
  auxvar%air_pressure  = UNINITIALIZED_DOUBLE

  auxvar%molv_air      = UNINITIALIZED_DOUBLE
  allocate(auxvar%dmolv_air(option%nflowspec))
  auxvar%dmolv_air(:)  = UNINITIALIZED_DOUBLE

  auxvar%pres_fh2o     = UNINITIALIZED_DOUBLE
  allocate(auxvar%dpres_fh2o(option%nflowspec))
  auxvar%dpres_fh2o(:) = UNINITIALIZED_DOUBLE

  auxvar%vis       = UNINITIALIZED_DOUBLE
  auxvar%kvr       = UNINITIALIZED_DOUBLE
  allocate(auxvar%dkvr(option%nflowspec))
  auxvar%dkvr(:)   = UNINITIALIZED_DOUBLE

  auxvar%Dk_eff    = UNINITIALIZED_DOUBLE
  allocate(auxvar%dDk_eff(option%nflowspec))
  auxvar%dDK_eff(:)= UNINITIALIZED_DOUBLE

  option%flow%numerical_derivatives = PETSC_FALSE

  call AuxVarFlowEnergyInit(auxvar,option)

end subroutine MpFlowAuxVarInit

! ************************************************************************** !

subroutine MpFlowAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 

  use Option_module

  implicit none
  
  type(flow_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%air_pressure = auxvar%air_pressure
  auxvar2%transient_por= auxvar%transient_por

  auxvar2%molv_air     = auxvar%molv_air
  auxvar2%dmolv_air    = auxvar%dmolv_air

  auxvar2%pres_fh2o    = auxvar%pres_fh2o
  auxvar2%dpres_fh2o   = auxvar%dpres_fh2o

  auxvar2%vis  = auxvar%vis
  auxvar2%kvr  = auxvar%kvr
  auxvar2%dkvr = auxvar%dkvr

  auxvar2%Dk_eff  = auxvar%Dk_eff
  auxvar2%dDk_eff = auxvar%dDk_eff

  !
  auxvar2%temp = auxvar%temp
  auxvar2%H = auxvar%H
  auxvar2%U = auxvar%U

  auxvar2%has_derivs = auxvar%has_derivs
  if(auxvar%has_derivs) then
    auxvar2%D_H = auxvar%D_H
    auxvar2%D_U = auxvar%D_U
  endif

  auxvar2%pres = auxvar%pres
  auxvar2%pc = auxvar%pc
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%mobility = auxvar%mobility
  auxvar2%viscosity = auxvar%viscosity

  if (auxvar%has_derivs) then
    auxvar2%D_pres = auxvar%D_pres
    auxvar2%D_sat  = auxvar%D_sat
    auxvar2%D_pc   = auxvar%D_pc
    auxvar2%D_den  = auxvar%D_den
    auxvar2%D_den_kg   = auxvar%D_den_kg
    auxvar2%D_mobility = auxvar%D_mobility
    auxvar2%D_por      = auxvar%D_por
  endif

end subroutine MpFlowAuxVarCopy

! ************************************************************************** !

subroutine MpFlowAuxVarCompute(x, auxvar,              &
                           auxvar_update, global_auxvar, &
                           material_auxvar,              &
                           characteristic_curves,        &
                           th_parameter,                 &
                           option)
  ! 
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 11/16/11
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type), pointer :: characteristic_curves
  PetscReal, intent(in) :: x(option%nflowdof)
  type(flow_auxvar_type) :: auxvar
  PetscBool, intent(in) :: auxvar_update
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(flow_parameter_type) :: th_parameter

  ! local variables
  PetscErrorCode :: ierr

  PetscReal :: pres_l, dpresl_dp
  PetscReal :: temperature

  PetscReal :: dw_kg, dw_mol, dw_dp, dw_dt            ! liq. water density
  PetscReal :: den_ice, dden_ice_dp, dden_ice_dt      ! ice density
  PetscReal :: den_air, dden_air_dp, dden_air_dt      ! air density
  PetscReal :: psat, dpsat_dt                         ! vapor pressure
  PetscReal :: molv_air, dmolv_air_dt, dmolv_air_dp   ! vapor mole fraction in air

  PetscReal :: sl, dsl_dp, dsl_dt
  PetscReal :: sg, dsg_dp, dsg_dt
  PetscReal :: si, dsi_dp, dsi_dt

  PetscReal :: kr, dkr_dp, dkr_dt                ! liq. water permissivity
  PetscReal :: visl, dvis_dp, dvis_dt            ! liq. water viscosity

  PetscReal :: h, h_dp, h_dt                     ! liq. water enthalpy
  PetscReal :: u, u_dp, u_dt                     ! liq. water internal energy
  PetscReal :: u_ice, du_ice_dp, du_ice_dt       ! ice water internal energy
  PetscReal :: u_air, du_air_dp, du_air_dt       ! air (inc. vapor) internal energy
  PetscReal :: p_g, Tk_g
  PetscReal :: C_g
  PetscReal, parameter :: C_a  = 1.860d-3        ! in MJ/kg/K at 300K, air
  PetscReal, parameter :: C_wv = 1.005d-3        ! in MJ/kg/K at 300K, vapor

  PetscReal :: Ke, dKe_dp, dKe_dt
  PetscReal :: Ke_fr, dKe_fr_dp, dKe_fr_dt
  PetscReal :: alpha
  PetscReal :: alpha_fr
  PetscReal :: Dk_wet
  PetscReal :: Dk_dry
  PetscReal :: Dk_ice

  PetscBool :: DTRUNC_FLAG = PETSC_TRUE  ! option for truncating deriatives to zero at bounds (default: TRUE)
  PetscReal :: tcmin, tcmax, pcmax, pcmin
  PetscReal :: dt_trunc, dp_trunc

  !----------------------------------------------------------------------------------------------------------------

  ! air_pressure at water-air interface
  auxvar%air_pressure = option%reference_pressure ! (TODO: in a closed-system, this is NOT the case)
  !
  pres_l      = x(1)
  ! Check if the capillary pressure is less than -100MPa
  pcmax = abs(characteristic_curves%saturation_function%pcmax)  ! non-negative
  if (pres_l - auxvar%air_pressure <= -pcmax) then
    pres_l = auxvar%air_pressure - pcmax
  endif
  auxvar%pres(LIQUID_PHASE) = pres_l
  auxvar%pc = max(0.d0, auxvar%air_pressure - pres_l)  ! always non-negative

  pcmin  = 0.d0  ! near-saturation PC zone (hint: a small positive value may be helpful ?)
  if (auxvar%pc <= pcmin) then
    auxvar%pc = 0.d0
    dpresl_dp = 1.d0
  else
    dpresl_dp = 0.d0
  endif
  !
  auxvar%temp = x(2)
  temperature = x(2)

  !----------------------------------------------------------------
  auxvar%kvr = 0.d0

!***************  P/T Bounds ***********************************************

  ! do the trunction of derivatives when assigning temporary values (local) to global variables (auxvar%, global_auxvar%)
  ! for P/T of out of bounds
  dt_trunc = 1.d0
  dp_trunc = 1.d0

  ! isothermal, i.e. hydrology only (for simplifying water density calculation)
  if(option%flow%isothermal) temperature = option%reference_temperature

  !--------------------------------------------------------------------

  ! general bounds of P/T. We may need to further limit bounds for 3-phase water properties.
  if (temperature>=100.d0 .or. temperature<=-273.d0) dt_trunc = 0.d0
  temperature    = min(100.d0, max(-273.d0, temperature))

  if (pres_l>=16.54d6) dp_trunc = 0.d0                         ! 16.54 MPa is upper limit for using IFC67 EOS.
  if (pres_l<=auxvar%air_pressure-pcmax) dp_trunc = 0.d0
  

!***************  Characteristic Curves ***********************************
  
  ! using modules in 'characteristic_curves.F90' to calculate needed variables:
  ! saturations and derivatives for 3 phases
  call MpFlowAuxVarComputeCharacteristicCurves(pres_l, temperature,    &
                                           global_auxvar,  auxvar,       &
                                           characteristic_curves,        &
                                           sl,  dsl_dp, dsl_dt,          &
                                           si,  dsi_dp, dsi_dt,          &
                                           sg,  dsg_dp, dsg_dt,          &
                                           auxvar%pres_fh2o,             &
                                           auxvar%dpres_fh2o(IFDOF1),    &
                                           auxvar%dpres_fh2o(IFDOF2),    &
                                           kr,  dkr_dp, dkr_dt,          &
                                           option)
  if(DTRUNC_FLAG) then
    dsl_dp = dsl_dp * dp_trunc
    dsl_dt = dsl_dt * dt_trunc
    dsi_dp = dsi_dp * dp_trunc
    dsi_dt = dsi_dt * dt_trunc
    dsg_dp = dsg_dp * dp_trunc
    dsg_dt = dsg_dt * dt_trunc

    auxvar%dpres_fh2o(IFDOF1) = auxvar%dpres_fh2o(IFDOF1) * dp_trunc
    auxvar%dpres_fh2o(IFDOF2) = auxvar%dpres_fh2o(IFDOF2) * dt_trunc
    dkr_dp = dkr_dp * dp_trunc
    dkr_dt = dkr_dt * dt_trunc

  endif


!***************  3-phase water properties **********************************************************

  ! ----- Liq. water ---------------------------------------------------------
  auxvar%sat(LIQUID_PHASE)           = sl
  auxvar%D_sat(LIQUID_PHASE, IFDOF1) = dsl_dp
  auxvar%D_sat(LIQUID_PHASE, IFDOF2) = dsl_dt


  ! Liq. water density/energy from available EOS (EOS_water.F90) could be crazy when below 0oC,
  ! So a trunction at/below 0.1oC, it may be helpful to avoid those
  tcmin = 0.1d0
  !--------- LIQ. Water Density
  call EOSWaterDensity(max(tcmin, temperature),               &
                       max(auxvar%air_pressure,pres_l),       &
                        dw_kg, dw_mol, dw_dp, dw_dt,ierr)
  if (DTRUNC_FLAG) dw_dp = dw_dp * dpresl_dp
  if (DTRUNC_FLAG .and. temperature<tcmin) dw_dt = 0.d0

  auxvar%den(LIQUID_PHASE)           = dw_mol
  auxvar%D_den(LIQUID_PHASE, IFDOF1) = dw_dp
  auxvar%D_den(LIQUID_PHASE, IFDOF2) = dw_dt

  auxvar%den_kg(LIQUID_PHASE)           = dw_mol*FMWH2O  ! Not-use 'dw_kg' in case inconsistent 'FMWH2O'
  auxvar%D_den_kg(LIQUID_PHASE, IFDOF1) = dw_dp*FMWH2O
  auxvar%D_den_kg(LIQUID_PHASE, IFDOF2) = dw_dt*FMWH2O

  !----------- LIQ. Water Energy
  call EOSWaterEnthalpy(max(tcmin, temperature),               &
                        max(auxvar%air_pressure,pres_l),       &
                        h, h_dp, h_dt, ierr)
  if (DTRUNC_FLAG) h_dp = h_dp * dpresl_dp
  if (DTRUNC_FLAG .and. temperature<tcmin) h_dt = 0.d0

  auxvar%H(LIQUID_PHASE)           = h * 1.d-6         ! J/kmol -> MJ/kmol
  auxvar%D_H(LIQUID_PHASE, IFDOF1) = h_dp * 1.d-6
  auxvar%D_H(LIQUID_PHASE, IFDOF2) = h_dt * 1.d-6

  u = h - pres_l/dw_mol
  u_dp = h_dp - (dpresl_dp/dw_mol - pres_l/(dw_mol*dw_mol)*dw_dp)
  u_dt = h_dt + pres_l/(dw_mol*dw_mol)*dw_dt

  auxvar%U(LIQUID_PHASE)           = u * 1.d-6         ! J/kmol -> MJ/kmol
  auxvar%D_U(LIQUID_PHASE, IFDOF1) = u_dp * 1.d-6
  auxvar%D_U(LIQUID_PHASE, IFDOF2) = u_dt * 1.d-6

  !----------- LIQ. Water Viscosity
  ! A note here (F.-M. Yuan: 2017-01-17)
  ! The Viscosity Eq. shows that: temp< ~ -63oC (1atm), 'visl' sharply increases starting from ~ 1.e-2 order.
  ! 'visl' ~ 0. around -133oC, which produces 'inf' for 'kr'
  tcmin = -63.d0

  call EOSWaterSaturationPressure(max(tcmin,temperature), psat, dpsat_dt, ierr)
  ! the lowest Tk of 200 for vapor exists in EOS-h2o phase-diagram, but here make it consistent with 'Viscosity'
  if(DTRUNC_FLAG .and. temperature<=tcmin) dpsat_dt = 0.d0

  call EOSWaterViscosity(max(tcmin,temperature),       &
                         pres_l,                       &
                         psat, dpsat_dt,               &
                         visl, dvis_dt,dvis_dp, ierr)
  if(DTRUNC_FLAG .and. temperature<=tcmin) dvis_dt = 0.d0
  if(DTRUNC_FLAG) dvis_dp = dvis_dp*dpresl_dp

  auxvar%vis = visl
  auxvar%kvr = kr/visl
  auxvar%dkvr(IFDOF1) = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  auxvar%dkvr(IFDOF2) = -kr/(visl*visl)*dvis_dt + dkr_dt/visl

  ! ----- ice water ---------------------------------------------------------
  auxvar%sat(SOLID_PHASE)           = si
  auxvar%D_sat(SOLID_PHASE, IFDOF1) = dsi_dp
  auxvar%D_sat(SOLID_PHASE, IFDOF2) = dsi_dt

  ! for ice Ih, Tk limit is ~127K (-146oC) in EOS-h2o phase-diagram
  tcmin = -146.d0
  tcmax = 0.d0
  !---------------
  call EOSWaterIceDensity(max(tcmin, min(tcmax,temperature)),    &
                          max(auxvar%air_pressure,pres_l), &
                          den_ice, dden_ice_dt, dden_ice_dp, ierr)
  if(DTRUNC_FLAG) dden_ice_dp = dden_ice_dp * dpresl_dp
  if(DTRUNC_FLAG .and. (temperature<=tcmin .or. temperature>=tcmax)) dden_ice_dT = 0.d0

  auxvar%den(SOLID_PHASE)           = den_ice
  auxvar%D_den(SOLID_PHASE, IFDOF1) = dden_ice_dp
  auxvar%D_den(SOLID_PHASE, IFDOF2) = dden_ice_dt

  auxvar%den_kg(SOLID_PHASE)           = den_ice*FMWH2O  ! Not-use 'dw_kg' in case inconsistent 'FMWH2O'
  auxvar%D_den_kg(SOLID_PHASE, IFDOF1) = dden_ice_dp*FMWH2O
  auxvar%D_den_kg(SOLID_PHASE, IFDOF2) = dden_ice_dt*FMWH2O

  !--------------
  call EOSWaterIceInternalEnergy(max(tcmin, min(tcmax,temperature)), &
                                 max(auxvar%air_pressure,pres_l), &
                                 PETSC_TRUE, &
                                 u_ice, du_ice_dp, du_ice_dt)
  if(DTRUNC_FLAG .and. (temperature<=tcmin .or. temperature>=tcmax)) du_ice_dt = 0.d0

  auxvar%U(SOLID_PHASE)           = u_ice * 1.d-6         ! J/kmol -> MJ/kmol
  auxvar%D_U(SOLID_PHASE, IFDOF1) = du_ice_dp * 1.d-6
  auxvar%D_U(SOLID_PHASE, IFDOF2) = du_ice_dt * 1.d-6

  ! ----- Air (incl. vapor) ------------------------------------------------------------------------------------------------

  auxvar%sat(GAS_PHASE)           = sg
  auxvar%D_sat(GAS_PHASE, IFDOF1) = dsg_dp
  auxvar%D_sat(GAS_PHASE, IFDOF2) = dsg_dt

  ! Calculate the values and derivatives for air density and its internal energy

  p_g      = max(pres_l, auxvar%air_pressure)
  ! the lowest Tk of ~200K for vapor exists in EOS-h2o phase-diagram
  tcmin    = -73.d0
  Tk_g     = max(tcmin,temperature)+TC2TK

  !---------
  den_air     = p_g/(IDEAL_GAS_CONSTANT*tk_g)*1.d-3        ! in kmol/m3 for all air-mixture
  dden_air_dp = -p_g/(IDEAL_GAS_CONSTANT*tk_g**2)*1.d-3
  dden_air_dt = 1.d0/(IDEAL_GAS_CONSTANT*tk_g)*1.d-3
  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    ! the following is a MUST for reducing tiny-time step (NOT sure why).
    dden_air_dp = dden_air_dp*dpresl_dp    ! w.r.t from 'pw' to 'pres(1)' upon soil total saturation
    if (tk_g<=tcmin) dden_air_dt = 0.d0
  endif

  !----------
  ! NOTE: vapor 'molv_air' is included in 'den_air', 'u_air'
  ! (because 'molv_air', fraction of vapor in air-mixture, going to be as multiplier in all 'air' calculations in 'MpFlow.F90')
#ifdef NO_VAPOR_DIFFUSION
  molv_air     = 0.d0   ! no vapor assumed in open-air
  dmolv_air_dt = 0.d0
  dmolv_air_dp = 0.d0

#else
  call EOSWaterSaturationPressure(Tk_g-TC2TK, psat, dpsat_dt, ierr)
  molv_air     = psat/p_g
  dmolv_air_dt = dpsat_dt/p_g
  dmolv_air_dp = -psat/p_g/p_g
  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    dmolv_air_dp = dmolv_air_dp * dpresl_dp
    if (tk_g<=tcmin) then
      dpsat_dt     = 0.d0
      dmolv_air_dt = 0.d0
    endif
  endif
#endif

  if(DTRUNC_FLAG) then
    dmolv_air_dp = dmolv_air_dp * dp_trunc
    dmolv_air_dt = dmolv_air_dt * dt_trunc
  endif
  auxvar%molv_air          = molv_air
  auxvar%dmolv_air(IFDOF1) = dmolv_air_dp
  auxvar%dmolv_air(IFDOF2) = dmolv_air_dt

  auxvar%den(GAS_PHASE)           = den_air*molv_air            ! Vapor only
  auxvar%D_den(GAS_PHASE, IFDOF1) = den_air*dmolv_air_dp + dden_air_dp*molv_air
  auxvar%D_den(GAS_PHASE, IFDOF2) = den_air*dmolv_air_dt + dden_air_dt*molv_air

  !-----------------------
  C_g       = C_wv*molv_air*FMWH2O + C_a*(1.d0 - molv_air)*FMWAIR          ! in MJ/kmol/K
  u_air     = C_g*tk_g                                                     ! in MJ/kmol
  du_air_dp = (C_wv*FMWH2O-C_a*FMWAIR)*dmolv_air_dp*tk_g
  du_air_dt = C_g + (C_wv*dmolv_air_dt*FMWH2O - C_a*dmolv_air_dt*FMWAIR)*tk_g
  if(DTRUNC_FLAG) then
    du_air_dp = du_air_dp * dp_trunc
    du_air_dt = du_air_dt * dt_trunc
  endif
  auxvar%U(GAS_PHASE)           = u_air
  auxvar%D_U(GAS_PHASE, IFDOF1) = du_air_dp
  auxvar%D_U(GAS_PHASE, IFDOF2) = du_air_dt




  ! ----- Thermal conductivity (effective) ---------------------------------------------------------

  ! Parameters for computation of effective thermal conductivity
  alpha    = th_parameter%alpha
  alpha_fr = th_parameter%alpha_fr
  Dk_wet = th_parameter%ckwet
  Dk_dry = th_parameter%ckdry
  Dk_ice = th_parameter%ckfrozen

  !Soil Kersten number
  Ke = (sl + eps)**(alpha)
  Ke_fr = (si + eps)**(alpha_fr)

  ! Derivative of Kersten number
  dKe_dp = alpha*(sl+eps)**(alpha-1.d0)*dsl_dp
  dKe_dt = alpha*(sl+eps)**(alpha-1.d0)*dsl_dt
  dKe_fr_dt = alpha_fr*(si + eps)**(alpha_fr - 1.d0)*dsi_dt
  dKe_fr_dp = alpha_fr*(si + eps)**(alpha_fr - 1.d0)*dsi_dp

  ! Effective thermal conductivity
  auxvar%Dk_eff = Dk_wet*Ke + Dk_ice*Ke_fr + (1.d0 - Ke - Ke_fr)*Dk_dry
  ! derivative of 'Dk_eff': Dk_eff = Dk_dry + (Dk_wet-Dk_dry)*Ke + (Dk_ice-Dk_dry)*Ke_fr
  auxvar%dDk_eff(IFDOF1) = (Dk_wet-Dk_dry)*dKe_dp + (Dk_ice-Dk_dry)*dKe_fr_dp
  auxvar%dDk_eff(IFDOF2) = (Dk_wet-Dk_dry)*dKe_dt + (Dk_ice-Dk_dry)*dKe_fr_dt


  ! isothermal, i.e. hydrology only
  if(option%flow%isothermal) then
    ! NO thermal conductivity and zeroing thermal states
    auxvar%Dk_eff = 0.d0
    auxvar%H      = 0.d0
    auxvar%U      = 0.d0

    ! zeroing thermal states derivatives
    auxvar%dDk_eff = 0.d0
    auxvar%D_H     = 0.d0
    auxvar%D_U     = 0.d0

  endif

 ! thermal process only
  if(option%flow%onlythermal) then
    ! no flow (BUT do allowing water states changing with thermal process, if any)
    auxvar%vis = 0.d0
    auxvar%kvr = 0.d0
    auxvar%dkvr= 0.d0
  endif

  ! Assign values to GLOBAL_auxvar, which shared with other MODEs in PFLOTRAN
  ! when instruted so
  if (auxvar_update) then

    global_auxvar%pres(LIQUID_PHASE) = x(1)
    global_auxvar%pres(SOLID_PHASE)  = auxvar%pres_fh2o
    global_auxvar%pres(GAS_PHASE)    = p_g
    global_auxvar%temp               = x(2)

    global_auxvar%sat    = auxvar%sat
    global_auxvar%den    = auxvar%den
    global_auxvar%den_kg = auxvar%den_kg
  endif

end subroutine MpFlowAuxVarCompute

! ************************************************************************** !
subroutine MpFlowAuxVarComputeCharacteristicCurves( presl,  tc,           &
                                    global_auxvar, auxvar,                  &
                                    characteristic_curves,                  &
                                    sl,  dsl_dpl, dsl_dt,                   &
                                    si,  dsi_dpl, dsi_dt,                   &
                                    sg,  dsg_dpl, dsg_dt,                   &
                                    ice_presl, ice_presl_dpl, ice_presl_dt, &
                                    kr,  dkr_dpl, dkr_dt,                   &
                                    option)
  !
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present
  !
  ! Revised by fengming Yuan @03-08-2016/CCSI-ONRL
  !       (1) ANY saturation_function, from 'Characteristic_Curves_module'
  !       (2) ANY permissivity function, from 'Characteristci_Curves_module' as well

  use Option_module
  use Global_Aux_module
  use Characteristic_Curves_module
  use EOS_Water_module

  implicit none

  type(option_type) :: option
  PetscReal, intent(in) :: presl    ! unit: Pa (liq water-air interface pressure: -pc+air_pressure)
  PetscReal, intent(in) :: tc       ! unit: oC
  type(flow_auxvar_type), intent(in) :: auxvar
  type(global_auxvar_type), intent(in) :: global_auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal, intent(out) :: sl,  dsl_dpl, dsl_dt
  PetscReal, intent(out) :: si,  dsi_dpl, dsi_dt
  PetscReal, intent(out) :: ice_presl, ice_presl_dpl, ice_presl_dt
  PetscReal, intent(out) :: sg,  dsg_dpl, dsg_dt
  PetscReal, intent(out) :: kr,  dkr_dpl, dkr_dt

  ! local variables

  PetscReal :: pc
  PetscReal :: sli, dsli_dpl, dsli_dt, pcli, presli
  PetscReal :: pcice, dpcice_dpl, dpcice_dt, slx, dslx_dx
  PetscReal :: dkr_dsl

  PetscReal :: sr, se, dse_dpl, dse_dt
  PetscReal :: funcB, dfuncB_dpl, funcA, dfuncA_dpl, dfuncA_dt

  PetscReal :: icedpore, icedpore_sl, icedpore_sl_dpl, icedpore_sl_dt
  PetscReal :: icedpore_pc, icedpore_pc_dsl

  ! ----------------
  ! (0) inputs
  pc = max(0.d0, auxvar%air_pressure - presl)   ! always non-negative (0 = saturated)
  if (pc > abs(characteristic_curves%saturation_function%pcmax)) then
    pc = characteristic_curves%saturation_function%pcmax
  endif

  !
  ! (1) saturations
  sl      = 0.d0  ! init to zero
  dsl_dpl = 0.d0
  dsl_dt  = 0.d0
  si      = 0.d0
  dsi_dpl = 0.d0
  dsi_dt  = 0.d0
  sg      = 0.d0
  dsg_dpl = 0.d0
  dsg_dt  = 0.d0

  !--------------------------------------------------------------------------
  ! liq.+ice saturation and its derivatives, under current 'presl' and 'tc'
  ! NOTE: 'presl' IS that from DOF1, currently is total water pressure (inc. ice).
  call characteristic_curves%saturation_function%Saturation(pc, sl, dsl_dpl, option)
  dsl_dt = 0.d0

  ! liq. only water PC at soil-ice/pore interface (positive, if no ice it's pc), under 'presl' and 'tc'
  ! NOTE: this is the liq. PC on iced pore/soil interface, and it must be not over total water PC (i.e. <=pc)
  call characteristic_curves%saturation_function%IceCapillaryPressure(presl, tc, &
                                   pcice, dpcice_dpl, dpcice_dt, option) ! w.r.t 'pressure' already in %IceCapillaryPressure()
  ice_presl    = (presl+pc) - pcice     ! 'pres_l+pc' is the liq. water column head (i.e. saturated column+atm.P)
  ice_presl_dpl= -dpcice_dpl            ! dpresl_dpl = 1, dpc_dpl = -1
  ice_presl_dt = -dpcice_dt             ! dpresl_dt = 0, dpc_dt = 0

  !--------------------------------------------------------------------------
  ! if ice module turns on, 3-phase saturation recalculated (liq. and ice)
  if (option%flow%ice_model /= UNINITIALIZED_INTEGER) then

    !----------------
    ! assuming all 'sli' in liq. at first, and separate them later
    sli      = sl
    pcli    = pc
    presli  = presl
    dsli_dt = 0.d0
    dsli_dpl= dsl_dpl

    !----------------
    ! update liq. saturation and its derivatives, under ice-adjusted capillary pressure, 'pcice'
    call characteristic_curves%saturation_function%Saturation(pcice, slx, dslx_dx, option)   ! pc on ice ---> sl on ice, but '_dx' is w.r.t. '_dpres'
    sl     = slx
    dsl_dpl= dslx_dx * ice_presl_dpl       ! or, dslx_dpl * dpl_dpc * dpcice_dpl, i.e. here '_dx' is '_dpl', 'pcice' is 'pc', and 'dpl_dpc=-1'
    dsl_dt = dslx_dx * ice_presl_dt

    icedpore_sl    = sl
    icedpore_sl_dpl= dsl_dpl
    icedpore_sl_dt = dsl_dt

    !----------------
    ! ice satuation and its derivatives
    select case (option%flow%ice_model)
      case (PAINTER_EXPLICIT)

        funcB = 1.d0
        dfuncB_dpl= 0.d0
        sr = characteristic_curves%saturation_function%Sr
        if (pc>0.d0) then
          se = (sli - sr)/(1.0d0 - sr)
          dse_dpl = dsli_dpl/(1.0d0 - sr)
          ! note: dsli_dt = 0, so ignored there for funcB
          funcB     = 1.d0/(se+1.d-20)
            ! must be inclusive of se=0, then in which case must be added a tiny value
            ! (1.d-20 is good, but any tiny <1.d-15 could cause error-checking below)
          dfuncB_dpl= -1.d0/((se+1.d-20)*(se+1.d-20))* dse_dpl
        endif
        !
        funcA = 1.d0
        dfuncA_dpl= 0.d0
        dfuncA_dt = 0.d0
        if(tc<=0.d0) then
          se = (sl - sr)/(1.0d0 - sr)
          dse_dpl = dsl_dpl/(1.0d0 - sr)
          dse_dt  = dsl_dt/(1.0d0 - sr)

          if (se>=0.d0) then
            funcA = 1.d0/(se+1.d-20)
            ! must be inclusive of se=0, then in which case must be added a tiny value
            ! (1.d-20 is good, but any tiny <1.d-15 could cause error-checking below)
            dfuncA_dpl = -1.d0/((se+1.d-20)*(se+1.d-20))* dse_dpl
            dfuncA_dt  = -1.d0/((se+1.d-20)*(se+1.d-20))* dse_dt
          endif
        endif

        se = 1.d0/(funcA + funcB - 1.d0)
        sl = se*(1.0d0 - sr) + sr
        dsl_dpl = - 1.d0/((funcA + funcB - 1.d0)*(funcA + funcB - 1.d0)) * &
                  (dfuncB_dpl + dfuncA_dpl)
        dsl_dpl = dsl_dpl*(1.0d0 - sr)
        dsl_dt  = - 1.d0/((funcA + funcB - 1.d0)*(funcA + funcB - 1.d0)) * &
                  dfuncA_dt
        dsl_dt  = dsl_dt*(1.0d0 - sr)

#if 0
        ! the following is from the original, but has mass-balance issue when freezing
        ! so now is off
        si = sl*(funcA - 1.d0)
        dsi_dpl = dsl_dpl*(funcA - 1.d0) + sl*dfuncA_dpl
        dsi_dt  = dsl_dt*(funcA - 1.d0)  + sl*dfuncA_dt
#else
        si = sli - sl
        dsi_dpl = dsli_dpl - dsl_dpl
        dsi_dt  = dsli_dt  - dsl_dt
#endif

      case (PAINTER_KARRA_EXPLICIT, PAINTER_KARRA_EXPLICIT_SMOOTH)

        ! ice satuation and its derivatives
#if 0
        ! the following is not going to work mathmatically
        ! e.g. assuming frozen soil, but with a very small ice-water mixture, i.e. sli small
        ! it's likely calculated into very high 'si' value over total 'sli'.
        if (sli>0.d0) then
          si     = 1.d0 - sl/sli             ! P.-K. Eq.(19)
          dsi_dt = -1.d0/sli*dsl_dt
          dsi_dpl= (sl*dsli_dpl-sli*dsl_dpl)/(sli**2)
        else
          si     = 0.d0
          dsi_dt = 0.d0
          dsi_dpl= 0.d0
        endif
#else
        si = sli - sl
        dsi_dpl = dsli_dpl - dsl_dpl
        dsi_dt  = dsli_dt  - dsl_dt
#endif

      case (DALL_AMICO)
        ! Model from Dall'Amico (2010) and Dall' Amico et al. (2011)
        ! ice satuation and its derivatives
        si     = sli - sl
        dsi_dpl= dsli_dpl - dsl_dpl
        dsi_dt = dsli_dt  - dsl_dt

      case default
        option%io_buffer = 'Ice module NOT recognized'
        call printErrMsg(option)

    end select

    ! re-calculating ice_presl
    !  there is a controversial of 'ice_presl' above and 'sl'
    !  (1) 'ice_presl' is that from 'sl', but with porosity of total
    !  (2) assuming a fully saturated (both ice/liq) situation, boundary with a water table
    !         'sl' is small, and 'ice_presl' is small too.
    !          So the water will be moving from boundaried water table into soil, which definitely
    !           is a stiff problem to solve
    !  (3) here we re-calculate 'ice_presl' assuming porosity reduced by ice

    if(si>0.d0 .and. si<1.d0) then
      ! sl maybe enlarge due to shrinking pore by ice
      icedpore = (1.0d0 - si)
      icedpore_sl = sl/icedpore

      icedpore_sl_dpl = (dsl_dpl*icedpore - sl*(-dsi_dpl)) &
                           /icedpore/icedpore
      icedpore_sl_dt  = (dsl_dt *icedpore - sl*(-dsi_dt))  &
                           /icedpore/icedpore
      call characteristic_curves%saturation_function%CapillaryPressure(icedpore_sl, &
                                icedpore_pc, icedpore_pc_dsl, option)

      ice_presl     = (presl+pc) - icedpore_pc
      ice_presl_dpl = -icedpore_pc_dsl * icedpore_sl_dpl
      ice_presl_dt  = -icedpore_pc_dsl * icedpore_sl_dt
    endif

  ! endif ice module turns on, 3-phase saturation recalculated (liq. and ice)
  endif
  !--------------------------------------------------------------------------

  !----------------
  ! air saturation as difference
  sg      = 1.d0 - sl - si
  dsg_dpl = -dsl_dpl - dsi_dpl
  dsg_dt  = -dsl_dt - dsi_dt

  ! some checking
  if(sl/=sl .or. abs(sl)>huge(sl) .or. &
     si/=si .or. abs(si)>huge(si) .or. &
     sg/=sg .or. abs(sg)>huge(sg)) then
    print *, 'checking Saturation cal. in MpFlow: ', sl, sg, si, presl, pc, ice_presl, pcice
    option%io_buffer = 'MpFlow with characteristic curve: NaN or INF properties'
    call printErrMsg(option)
  endif

  ! Check for bounds on saturations
  if ((sl-1.d0)>1.d-15 .or. sl<-1.d-15) then
    print *, tc, pc, pcice, sli, sl, si, sg
    option%io_buffer = 'MpFlow with ice mode: LIQUID Saturation error: >1 or <0'
    call printErrMsg(option)
  endif

  if ((si-1.d0)>1.d-15 .or. si<-1.d-15) then
    print *, tc, pc, pcice, sli, sl, si, sg
    option%io_buffer = 'MpFlow with ice mode: ICE Saturation error:  >1 or <0'
    call printErrMsg(option)
  endif

  if ((sg-1.d0)>1.d-15 .or. sg<-1.d-15) then
    print *, tc, pc, pcice, sli, sl, si, sg
    option%io_buffer = 'MpFlow with ice mode: AIR Saturation error:  >1 or <0'
    call printErrMsg(option)
  endif
  if (abs((sl + si + sg)-1.d0)>1.d-15) then
    print *, tc, pc, pcice, sli, sl, si, sg
    option%io_buffer = 'MpFlow with ice mode: Saturation not summed to 1 '
    call printErrMsg(option)
  endif

  !--------------------------------------------------------------------------------------------------

  ! (2) relative permissivity of liq. only water in multiple-phase mixture
  kr      = 0.d0  !all initialized to zero
  dkr_dsl = 0.d0
  dkr_dt  = 0.d0
  dkr_dpl = 0.d0

#if 0
  ! the following will greatly reduce any liq. water flow in iced pore even if unfilled or saturated
  ! causing problem when initially freezing (which would allow expanding and pushing rest liq. water flow, but unless 'kr' large enough)
  call characteristic_curves%liq_rel_perm_function%RelativePermeability(sl, kr, dkr_dsl, option)
  dkr_dpl= dkr_dsl*dsl_dpl
  dkr_dt = dkr_dsl*dsl_dt

#else
  call characteristic_curves%liq_rel_perm_function%RelativePermeability(icedpore_sl, kr, dkr_dsl, option)
  dkr_dpl= dkr_dsl*icedpore_sl_dpl
  dkr_dt = dkr_dsl*icedpore_sl_dt

#endif

end subroutine MpFlowAuxVarComputeCharacteristicCurves


! ************************************************************************** !
subroutine DarcyFlowDerivative(flow_auxvar_up,      flow_auxvar_dn,       &
                               material_auxvar_up,  material_auxvar_dn,   &
                               flow_parameter_up,   flow_parameter_dn,    &
                               area, dist,                                &
                               option,                                    &
                               v_darcy,                                   &
                               q,                                         &
                               dq_up, dq_dn, ifderivative)

  !
  ! Computes the derivatives of Liquid Water Darcy Flow
  !
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
  !

  use Option_module
  use Material_Aux_class, only : material_auxvar_type
  use Connection_module

  implicit none

  type(Flow_auxvar_type)      :: flow_auxvar_up, flow_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(Flow_parameter_type)   :: flow_parameter_up, flow_parameter_dn
  PetscReal, intent(in) :: area, dist(-1:3)
  type(option_type) :: option
  PetscBool, intent(in)  :: ifderivative
  PetscReal, intent(out) :: v_darcy                ! unit: m/sec (default)
  PetscReal, intent(out) :: q                      ! unit: kmol/sec, + from up -> dn; - from dn -> up
  PetscReal, intent(out) :: dq_up(option%nflowdof)
  PetscReal, intent(out) :: dq_dn(option%nflowdof)

  ! local variables
  PetscReal :: dd_up, dd_dn, upweight
  PetscReal :: dist_gravity                                       ! distance along gravity vector
  PetscReal :: Dq, perm_up, perm_dn

  PetscReal :: presl_up
  PetscReal :: presl_dn

  PetscInt :: idof
  PetscReal :: gravity, dgravity_up(option%nflowdof), dgravity_dn(option%nflowdof)
  PetscReal :: dphi, ddphi_up(option%nflowdof), ddphi_dn(option%nflowdof)
  PetscReal :: ukvr, dukvr_up(option%nflowdof), dukvr_dn(option%nflowdof)
  PetscReal :: den, dden_up(option%nflowdof), dden_dn(option%nflowdof)  ! density in kmol/m3

  !--------------------------------------------------------------------------------------------
  q        = 0.d0
  dq_up(:) = 0.d0
  dq_dn(:) = 0.d0

  presl_up = flow_auxvar_up%pres(LIQUID_PHASE)
  presl_dn = flow_auxvar_dn%pres(LIQUID_PHASE)

  !-------------------------------------------------------------------------
  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)
  if (flow_auxvar_up%sat(LIQUID_PHASE) < flow_parameter_up%sir(LIQUID_PHASE)) then
    upweight = 0.d0
  elseif (flow_auxvar_dn%sat(LIQUID_PHASE) < flow_parameter_dn%sir(LIQUID_PHASE)) then
    upweight = 1.d0
  end if

  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  v_darcy = 0.d0
  ! Darcy Flow of liquid water
  if (flow_auxvar_up%sat(LIQUID_PHASE) > flow_parameter_up%sir(LIQUID_PHASE) .or.  &
      flow_auxvar_dn%sat(LIQUID_PHASE) > flow_parameter_dn%sir(LIQUID_PHASE)) then

    gravity = (upweight        *flow_auxvar_up%den_kg(LIQUID_PHASE)+  &
               (1.d0-upweight) *flow_auxvar_dn%den_kg(LIQUID_PHASE))  &
              * dist_gravity

#if 0
    dphi = presl_up - presl_dn + gravity
#else
    dphi = flow_auxvar_up%pres_fh2o - flow_auxvar_dn%pres_fh2o + gravity
#endif

    if (dphi>=0.d0) then
      den  = flow_auxvar_up%den(LIQUID_PHASE)
      ukvr = flow_auxvar_up%kvr
    else
      den  = flow_auxvar_dn%den(LIQUID_PHASE)
      ukvr = flow_auxvar_dn%kvr
    endif

    if (dabs(ukvr)>floweps) then
      v_darcy = Dq * ukvr * dphi  ! needed for output (m/s)
    endif
    q = v_darcy*area*den          ! m/s * m2 * kmol/m3 = kmol/s

    !-----------------------------------
    if(ifderivative) then
      do idof = 1, option%nflowdof
        dgravity_up(idof) = upweight*dist_gravity* &
          flow_auxvar_up%D_den_kg(LIQUID_PHASE,idof)

        dgravity_dn(idof) = (1.d0-upweight)*dist_gravity* &
          flow_auxvar_dn%D_den_kg(LIQUID_PHASE,idof)

#if 0
        ddphi_up(idof) = dgravity_up(idof) + flow_auxvar_up%D_pres(LIQUID_PHASE, idof)
        ddphi_dn(idof) = dgravity_dn(idof) - flow_auxvar_dn%D_pres(LIQUID_PHASE, idof)
#else
        ddphi_up(idof) = dgravity_up(idof) + flow_auxvar_up%dpres_fh2o(idof)
        ddphi_dn(idof) = dgravity_dn(idof) - flow_auxvar_dn%dpres_fh2o(idof)
#endif

        if (dphi>0.D0) then
          dukvr_up(idof) = flow_auxvar_up%dkvr(idof)
          dukvr_dn(idof) = 0.d0

          dden_up(idof)  = flow_auxvar_up%D_den(LIQUID_PHASE, idof)
          dden_dn(idof)  = 0.d0

        else
          dukvr_up(idof) = 0.d0
          dukvr_dn(idof) = flow_auxvar_dn%dkvr(idof)

          dden_up(idof)  = 0.d0
          dden_dn(idof)  = flow_auxvar_up%D_den(LIQUID_PHASE, idof)

        end if

        if (dabs(q)>floweps) then
          !q = Dq*area*ukvr*dphi*den
          dq_up(idof) = Dq*area* &
                       ( dukvr_up(idof)*dphi           *den          + &
                         ukvr          *ddphi_up(idof) *den          + &
                         ukvr          *dphi           *dden_up(idof) )


          dq_dn(idof) = Dq*area* &
                       ( dukvr_dn(idof)*dphi           *den          + &
                         ukvr          *ddphi_dn(idof) *den          + &
                         ukvr          *dphi           *dden_dn(idof) )

        end if

      end do

    end if

  end if

end subroutine DarcyFlowDerivative

! ************************************************************************** !
subroutine AdvectionDerivative(flow_auxvar_up,   flow_auxvar_dn,   &
                              q, dq_up, dq_dn,                     &
                              option,                              &
                              qe,                                  &
                              dqe_up, dqe_dn, ifderivative)

  !
  ! Computes the derivatives of energy advections, given flux and its derivatives
  !
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(Flow_auxvar_type):: flow_auxvar_up, flow_auxvar_dn
  PetscReal, intent(in) :: q    ! unit: kmol/sec, + from up -> dn; - from dn -> up
  PetscReal, intent(in) :: dq_up(option%nflowdof), dq_dn(option%nflowdof)
  PetscBool, intent(in) :: ifderivative
  PetscReal, intent(out) :: qe  ! unit: MJ/sec, + from up -> dn; - from dn -> up
  PetscReal, intent(out) :: dqe_up(option%nflowdof), dqe_dn(option%nflowdof)

  ! local variables

  PetscInt :: idof
  PetscReal :: uh,  duh_up(option%nflowdof), duh_dn(option%nflowdof)

  !--------------------------------------------------------------------------------------------
  qe        = 0.d0
  dqe_up(:) = 0.d0
  dqe_dn(:) = 0.d0

  !-------------------------------------------------------------------------
  if (q>=0.d0) then
    uh   = flow_auxvar_up%H(LIQUID_PHASE)    ! U or H (needs more thinking)
  else
    uh   = flow_auxvar_dn%H(LIQUID_PHASE)
  endif
  qe     = q*uh

  !-----------------------------------
  if(ifderivative) then
    do idof = 1, option%nflowdof
      if (q>=0.D0) then
        duh_up(idof)   = flow_auxvar_up%D_H(LIQUID_PHASE, idof)
        duh_dn(idof)   = 0.d0
      else
        duh_up(idof)   = 0.d0
        duh_dn(idof)   = flow_auxvar_dn%D_H(LIQUID_PHASE, idof)
      end if

      !qe = q*uh
      dqe_up(idof) = q           *duh_up(idof) + &
                     dq_up(idof) *uh

      dqe_dn(idof) = q           *duh_dn(idof) + &
                     dq_dn(idof) *uh
    end do
  end if

end subroutine AdvectionDerivative

! ************************************************************************** !

subroutine ConductionDerivative(flow_auxvar_up,   flow_auxvar_dn,      &
                              material_auxvar_up, material_auxvar_dn,  &
                              area, dist,                              &
                              option,                                  &
                              qe,                                      &
                              dqe_up, dqe_dn, ifderivative)

  use Option_module
  use Material_Aux_class, only : material_auxvar_type
  use Connection_module

  implicit none

  type(option_type) :: option
  type(Flow_auxvar_type):: flow_auxvar_up, flow_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  PetscBool, intent(in) :: ifderivative
  PetscReal, intent(in) :: area, dist(-1:3)
  PetscReal, intent(out) :: qe  ! unit: MJ/sec, + from up -> dn; - from dn -> up
  PetscReal, intent(out) :: dqe_up(option%nflowdof), dqe_dn(option%nflowdof)

  ! local variables
  PetscReal :: tc_up, tc_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: dist_gravity, upweight

  PetscInt :: idof
  PetscReal :: Dke_up, Dke_dn, Dke
  PetscReal :: dDke_up(option%nflowdof), dDke_dn(option%nflowdof)
  PetscReal :: dDk_up(option%nflowdof), dDk_dn(option%nflowdof)

  !--------------------------------------------------------------------------------------------
  tc_up = flow_auxvar_up%temp
  tc_dn = flow_auxvar_dn%temp

  qe        = 0.d0
  dqe_up(:) = 0.d0
  dqe_dn(:) = 0.d0

  !-------------------------------------------------------------------------
  call ConnectionCalculateDistances(dist, option%gravity, dd_up, dd_dn, &
                                    dist_gravity, upweight)

  Dke_up = flow_auxvar_up%Dk_eff
  Dke_dn = flow_auxvar_dn%Dk_eff

  Dke = 0.d0
  if(Dke_up /= 0.d0 .or. Dke_dn /= 0.d0) then
    ! effective thermal conductivity, Dk
    ! 1/Dke = dd_dn/Dke_dn + dd_up/Dke_up
    Dke = (Dke_up * Dke_dn) / (dd_dn*Dke_up + dd_up*Dke_dn)

    qe = Dke*area*(tc_up-tc_dn)
  end if

  if (ifderivative) then
    if(Dke /= 0.d0) then

      dDke_up = flow_auxvar_up%dDk_eff
      dDke_dn = flow_auxvar_dn%dDk_eff

      do idof = 1, option%nflowdof
       ! 1/Dke = dd_dn/Dke_dn + dd_up/Dke_up

        dDk_up(idof) = Dke*Dke * &
          dd_up*(dDke_up(idof)/Dke_up/Dke_up)  ! d(1/Dk)

        dDk_dn(idof) = Dke*Dke * &
          dd_dn*(dDke_dn(idof)/Dke_dn/Dke_dn)  ! d(1/Dk)

      end do
      ! cond = Dke*area*(flow_auxvar_up%temp-flow_auxvar_dn%temp)
      dqe_up(IFDOF1) = area * dDk_up(IFDOF1) * &
                      (tc_up-tc_dn)
      dqe_dn(IFDOF1) = area * dDk_dn(IFDOF1) * &
                      (tc_up-tc_dn)

      dqe_up(IFDOF2) = area * &
                      (dDk_up(IFDOF2) * (tc_up - tc_dn) + Dke)
      dqe_dn(IFDOF2) = area * &
                      (dDk_dn(IFDOF2) * (tc_up - tc_dn) - Dke)

    endif
  endif

end subroutine ConductionDerivative

! ************************************************************************** !
#if 0
subroutine VaporDiffusionDerivative(auxvar_up,global_auxvar_up, &
                            material_auxvar_up, &
                            auxvar_dn,global_auxvar_dn, &
                            material_auxvar_dn, &
                            area, &
                            dist, &
                            option, &
                            Jup,Jdn)

  !
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  !
  ! Author: ???
  ! Date: 12/13/07
  !

  use Option_module
  use Characteristic_Curves_module
  use Connection_module
  use EOS_Water_module
  use Utility_module

  implicit none

  type(Flow_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  PetscReal :: sir_up, sir_dn
  PetscReal :: area, dist(-1:3)
  type(option_type) :: option
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)

  ! local variables
  PetscReal :: dd_up, dd_dn, upweight
  PetscReal :: dist_gravity           ! distance along gravity vector

  ! ice/air variables
  PetscReal :: Ddiffgas_avg, Ddiffgas_up, Ddiffgas_dn
  PetscReal :: p_g
  PetscReal :: deng_up, deng_dn
  PetscReal :: molg_up, molg_dn
  PetscReal :: satg_up, satg_dn
  PetscReal :: Diffg_up, Diffg_dn
  PetscReal :: ddeng_dt_up, ddeng_dt_dn, ddeng_dp_up, ddeng_dp_dn
  PetscReal :: dmolg_dt_up, dmolg_dt_dn
  PetscReal :: dDiffg_dt_up, dDiffg_dt_dn
  PetscReal :: dDiffg_dp_up, dDiffg_dp_dn
  PetscReal :: dsatg_dp_up, dsatg_dp_dn, dsatg_dt_up, dsatg_dt_dn
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscReal :: dmolg_dp_up, dmolg_dp_dn
  PetscReal :: ugas_ave, dugas_ave_dt, dugas_ave_dp, fdiffgas, fdiffgas_dx

  !-------------------------------------------------------------------------
  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)

  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0

  Jup = 0.d0
  Jdn = 0.d0

  dmolg_dp_up = 0.d0
  dmolg_dp_dn = 0.d0
  dmolg_dt_up = 0.d0
  dmolg_dt_dn = 0.d0


    satg_up = auxvar_up%sat(GAS_PHASE)
    satg_dn = auxvar_dn%sat(GAS_PHASE)
    if ((satg_up > eps) .and. (satg_dn > eps)) then

      p_g = auxvar_up%pres(GAS_PHASE)
      deng_up = auxvar_up%ice%den_air
      deng_dn = auxvar_dn%ice%den_air

      Diffg_ref = 2.13D-5                 ! Reference diffusivity, need to read from input file
      p_ref = option%reference_pressure   ! in Pa
      T_ref = 25.d0       ! in deg C

      Diffg_up = Diffg_ref*(p_ref/p_g)*((max(-50d0, global_auxvar_up%temp) + TC2TK) &
                 /(T_ref + TC2TK))**(1.8d0)
      Diffg_dn = Diffg_ref*(p_ref/p_g)*((max(-50d0, global_auxvar_dn%temp) + TC2TK) &
                 /(T_ref + TC2TK))**(1.8d0)

      Ddiffgas_up = por_up*tor_up*satg_up*Diffg_up
      Ddiffgas_dn = por_dn*tor_dn*satg_dn*Diffg_dn

      molg_up = auxvar_up%ice%molv_air
      molg_dn = auxvar_dn%ice%molv_air
      dmolg_dt_up = auxvar_up%ice%dmolv_air_dt
      dmolg_dt_dn = auxvar_dn%ice%dmolv_air_dt
      dmolg_dp_up = auxvar_up%ice%dmolv_air_dp
      dmolg_dp_dn = auxvar_dn%ice%dmolv_air_dp

      ddeng_dt_up = auxvar_up%ice%dden_air_dt
      ddeng_dt_dn = auxvar_dn%ice%dden_air_dt
      ddeng_dp_up = auxvar_up%ice%dden_air_dp
      ddeng_dp_dn = auxvar_dn%ice%dden_air_dp

      dDiffg_dt_up = 1.8d0*Diffg_up/(max(-50.0d0,global_auxvar_up%temp) + TC2TK)
      dDiffg_dt_dn = 1.8d0*Diffg_dn/(max(-50.0d0,global_auxvar_dn%temp) + TC2TK)
      dDiffg_dp_up = 0.d0
      dDiffg_dp_dn = 0.d0

      dsatg_dp_up = auxvar_up%ice%dsat_air_dp
      dsatg_dp_dn = auxvar_dn%ice%dsat_air_dp
      dsatg_dt_up = auxvar_up%ice%dsat_air_dt
      dsatg_dt_dn = auxvar_dn%ice%dsat_air_dt

      if (deng_up*molg_up > deng_dn*molg_dn) then
      ! fmyuan: 'molg' is mole fraction of vapor in air, NOT vapor density
        upweight = 0.d0
        ugas_ave = auxvar_up%ice%u_air
        dugas_ave_dt = auxvar_up%ice%du_air_dt
        dugas_ave_dp = auxvar_up%ice%du_air_dp
      else
        upweight = 1.d0
        ugas_ave = auxvar_dn%ice%u_air
        dugas_ave_dt = auxvar_dn%ice%du_air_dt
        dugas_ave_dp = auxvar_dn%ice%du_air_dp
      endif

      Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn

      ! derivaives: Ddiffgas_avg*area*(deng_up*molg_up - deng_dn*molg_dn)/(dd_up + dd_dn)
      fdiffgas = area/(dd_up + dd_dn)*Ddiffgas_avg*(deng_up*molg_up - deng_dn*molg_dn)

      fdiffgas_dx = area/(dd_up+dd_dn) * &                                             !_up, _dp
                 ( Ddiffgas_avg * (ddeng_dp_up*molg_up + deng_up*dmolg_dp_up) + &
                   (deng_up*molg_up - deng_dn*molg_dn) * upweight * por_up*tor_up* &
                    (satg_up*dDiffg_dp_up + dsatg_dp_up*Diffg_up) )
      Jup(1,1) = Jup(1,1) + fdiffgas_dx
      Jup(2,1) = Jup(2,1) + (fdiffgas*dugas_ave_dp + ugas_ave*fdiffgas_dx)


      fdiffgas_dx = area/(dd_up+dd_dn) * &                                             !_up, _dt
                 ( Ddiffgas_avg * (ddeng_dt_up*molg_up + deng_up*dmolg_dt_up) + &
                   (deng_up*molg_up - deng_dn*molg_dn) * upweight * por_up*tor_up* &
                    (satg_up*dDiffg_dt_up + dsatg_dt_up*Diffg_up) )
      Jup(1,2) = Jup(1,2) + fdiffgas_dx
      Jup(2,2) = Jup(2,2) + (fdiffgas*dugas_ave_dt + ugas_ave*fdiffgas_dx)


      fdiffgas_dx = area/(dd_up+dd_dn) * &                                             !_dn, _dp
                 ( Ddiffgas_avg * (-ddeng_dp_dn*molg_dn - deng_dn*dmolg_dp_dn) + &
                   (deng_up*molg_up - deng_dn*molg_dn) * (1.d0-upweight) * por_dn*tor_dn* &
                    (satg_dn*dDiffg_dp_dn + dsatg_dp_dn*Diffg_dn) )
      Jdn(1,1) = Jdn(1,1) + fdiffgas_dx
      Jdn(2,1) = Jdn(2,1) + (fdiffgas*dugas_ave_dp + ugas_ave*fdiffgas_dx)


      fdiffgas_dx = area/(dd_up+dd_dn) * &                                             !_dn, _dt
                 ( Ddiffgas_avg * (-ddeng_dt_dn*molg_dn - deng_dn*dmolg_dt_dn) + &
                   (deng_up*molg_up - deng_dn*molg_dn) * (1.d0-upweight) * por_dn*tor_dn* &
                    (satg_dn*dDiffg_dt_dn + dsatg_dt_dn*Diffg_dn) )
      Jdn(1,2) = Jdn(1,2) + fdiffgas_dx
      Jdn(2,2) = Jdn(2,2) + (fdiffgas*dugas_ave_dt + ugas_ave*fdiffgas_dx)

    endif
  !

end subroutine VarporDiffusionDerivative
#endif

! ************************************************************************** !

subroutine MpFlowAuxVarDestroy(auxvar)
  ! 
  ! Deallocates a MpFlow auxiliary object's auxvar
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  use Utility_module, only: DeallocateArray

  implicit none

  type(flow_auxvar_type), pointer :: auxvar
  
  if (.not.associated(auxvar)) return

  call DeallocateArray(auxvar%pres)
  call DeallocateArray(auxvar%dmolv_air)
  call DeallocateArray(auxvar%dpres_fh2o)
  call DeallocateArray(auxvar%dkvr)
  call DeallocateArray(auxvar%dDk_eff)

  call AuxVarFlowEnergyStrip(auxvar)

  !if(associated(auxvar)) deallocate(auxvar)
  nullify(auxvar)
  
end subroutine MpFlowAuxVarDestroy

! ************************************************************************** !

subroutine MpFlowAuxDestroy(aux)
  ! 
  ! Deallocates all MpFlow auxiliary objects
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  implicit none

  type(MpFlow_type), pointer :: aux
  type(flow_auxvar_type), pointer :: auxvar
  type(flow_parameter_type), pointer :: flow_parameter
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) then
    do iaux = 1, size(aux%auxvars)
      auxvar => aux%auxvars(iaux)
      call MpFlowAuxVarDestroy(auxvar)
    end do
  endif
  nullify(aux%auxvars)

  if (associated(aux%flow_parameters)) then
    do iaux = 1, size(aux%flow_parameters)
      flow_parameter => aux%flow_parameters(iaux)

      if (associated(flow_parameter%sir)) deallocate(flow_parameter%sir)
      nullify(flow_parameter%sir)

      if (associated(flow_parameter%diffusion_coefficient)) &
        deallocate(flow_parameter%diffusion_coefficient)
      nullify(flow_parameter%diffusion_coefficient)

      if (associated(flow_parameter%diffusion_activation_energy)) &
        deallocate(flow_parameter%diffusion_activation_energy)
      nullify(flow_parameter%diffusion_activation_energy)

    enddo

    deallocate(aux%flow_parameters)
    nullify(aux%flow_parameters)

  endif
  
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)

  deallocate(aux)
  nullify(aux)  

  end subroutine MpFlowAuxDestroy

end module MpFlow_Aux_module
