module MpFlow_Aux_module
  ! data types and modules for MPFLOW MODE

  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use AuxVars_Flow_Energy_module
  use Option_module

  implicit none

  private 

  PetscReal, public :: flow_itol_scaled_res = 1.d-15
  PetscReal, public :: flow_itol_rel_update = UNINITIALIZED_DOUBLE

  ! indexing

  ! fluid(s)
  PetscInt, public :: IFLUID1 = 1      ! fluid 1: LIQ_FLUID (1)
  PetscInt, public :: IFLUID2 = 2      ! fluid 2: AIR_FLUID (2)

  ! currently only liq. water as flow species
  PetscInt, public :: IFLOWSPEC1 = 1       ! fluid spec 1: liq. water species
  PetscInt, public :: IFLOWSPEC2 = 2       ! fluid spec 2: air species
  !PetscInt, public :: IFLOWSPEC3 = 3       ! fluid spec 3: another AIR_FLUID, if any
 PetscReal, public :: FMW_FLUIDS(2) = [FMWH2O,FMWAIR]  ! should be appended, if more fluid identified

  ! 2 dofs currently: pressure, temperature
  PetscInt, public :: IFDOF1 = PRESSURE_DOF
  PetscInt, public :: IFDOF2 = TEMPERATURE_DOF
  PetscInt, public :: IFDOF3 = CONDUCTANCE_DOF ! thermal
  PetscInt, public :: IFDOF4 = ENTHALPY_DOF    ! not-yet

  type, public, extends(auxvar_flow_energy_type) :: mpflow_auxvar_type

    PetscReal :: air_pressure             ! unit: Pa ( air pressure of liq.water/ice-air interface. It's reference_pressure if an open system)

    PetscReal :: pres_fh2o                ! liq. water pressure on ice-included soil matrix
    PetscReal, pointer :: dpres_fh2o(:)   ! (nflowdof)

    !relative permissivity for all fluid(s)
    PetscReal, pointer :: kvr(:)          ! (nfluid)
    PetscReal, pointer :: dkvr(:,:)       ! (nfluid, nflowdof)

    ! effective thermal conductivity integrated for all fluid(s), including fluid changed to immobile or solid, and porous media for fluid(s)
    PetscReal :: Dk_eff
    PetscReal, pointer :: dDk_eff(:)      ! (nflowdof)

  end type mpflow_auxvar_type

  type, public :: mpflow_parameter_type
    PetscReal :: dencpr   ! porous media (e.g. soil/rock) particle density X specific_heat_capacity: MJ/m^3-K
    PetscReal :: ckdry    ! Thermal conductivity (dry)
    PetscReal :: ckwet    ! Thermal conductivity (saturated - wet)
    PetscReal :: alpha
    PetscReal :: ckfrozen ! Thermal conductivity (frozen/saturated soil)
    PetscReal :: alpha_fr ! exponent frozen

    PetscReal, pointer :: sr(:)                            !(nfluid) residue saturation
    PetscReal, pointer :: diffusion_coefficient(:)         !(nfluid)
    PetscReal, pointer :: diffusion_activation_energy(:)   !(nfluid)
  end type mpflow_parameter_type
  
  ! the following data-type is for domain-wide use
  type, public :: MpFlow_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(mpflow_parameter_type), pointer :: mpflow_parameters(:)    ! (material number)
    type(mpflow_auxvar_type), pointer :: auxvars(:)                 ! (gridcell number)
    type(mpflow_auxvar_type), pointer :: auxvars_bc(:)              ! (gridcell number)
    type(mpflow_auxvar_type), pointer :: auxvars_ss(:)              ! (gridcell number)
  end type MpFlow_type

  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

  public :: MpFlowAuxCreate, MpFlowAuxDestroy, &
            MpFlowAuxVarCompute, MpFlowAuxVarInit, &
            MpFlowAuxVarCopy, MpFlowAuxVarStrip

  public :: DarcyFlowDerivative, &
            AdvectionDerivative, &
            ConductionDerivative

contains

! ************************************************************************** !

function MpFlowAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019

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
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  aux%n_zero_rows = 0

  nullify(aux%mpflow_parameters)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  MpFlowAuxCreate => aux
  
end function MpFlowAuxCreate

! ************************************************************************** !

subroutine MpFlowAuxVarInit(auxvar,option)
  ! 
  ! Initialize a single auxiliary object
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  ! 

  use Option_module
  use PFLOTRAN_Constants_module, only : UNINITIALIZED_DOUBLE

  use AuxVars_Flow_Energy_module

  implicit none
  
  type(mpflow_auxvar_type) :: auxvar
  type(option_type) :: option

  !
  call AuxVarFlowEnergyInit(auxvar,option)

  !
  auxvar%air_pressure  = UNINITIALIZED_DOUBLE
  auxvar%pres_fh2o     = UNINITIALIZED_DOUBLE
  nullify(auxvar%kvr)
  allocate(auxvar%kvr(option%flow%nfluid))
  auxvar%kvr           = UNINITIALIZED_DOUBLE
  auxvar%Dk_eff        = UNINITIALIZED_DOUBLE

  nullify(auxvar%dpres_fh2o)
  nullify(auxvar%dkvr)
  nullify(auxvar%dDk_eff)
  if (.not.option%flow%numerical_derivatives) then
    auxvar%has_derivs = PETSC_TRUE    ! not needed, but just in case

    allocate(auxvar%dpres_fh2o(option%nflowdof))
    auxvar%dpres_fh2o    = UNINITIALIZED_DOUBLE

    allocate(auxvar%dkvr(option%flow%nfluid, option%nflowdof))
    auxvar%dkvr          = UNINITIALIZED_DOUBLE
    allocate(auxvar%dDk_eff(option%nflowdof))
    auxvar%dDK_eff       = UNINITIALIZED_DOUBLE
  endif

end subroutine MpFlowAuxVarInit

! ************************************************************************** !

subroutine MpFlowAuxVarCopy(auxvar,auxvar2)
  ! 
  ! Copies an auxiliary object
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  ! 

  implicit none
  
  class(mpflow_auxvar_type) :: auxvar, auxvar2

  call AuxVarFlowEnergyCopy(auxvar, auxvar2)

  auxvar2%air_pressure = auxvar%air_pressure
  auxvar2%pres_fh2o    = auxvar%pres_fh2o
  auxvar2%kvr  = auxvar%kvr
  auxvar2%Dk_eff  = auxvar%Dk_eff

  if (associated(auxvar%dpres_fh2o)) &
    auxvar2%dpres_fh2o   = auxvar%dpres_fh2o

  if (associated(auxvar%dkvr)) &
    auxvar2%dkvr = auxvar%dkvr

  if (associated(auxvar%dDk_eff)) &
    auxvar2%dDk_eff = auxvar%dDk_eff


end subroutine MpFlowAuxVarCopy

! ************************************************************************** !

subroutine MpFlowAuxVarCompute(x,                 &
                           characteristic_curves, &
                           mpflow_parameter,      &
                           auxvar,                &
                           option)
  ! 
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present for water fluid
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use Characteristic_Curves_module
  
  implicit none

  type(option_type) :: option
  PetscReal, intent(in) :: x(option%nflowdof)
  class(characteristic_curves_type), pointer :: characteristic_curves
  type(mpflow_parameter_type) :: mpflow_parameter
  type(mpflow_auxvar_type) :: auxvar

  ! local variables
  PetscErrorCode :: ierr

  PetscReal :: presl, presg, dpresl_dp
  PetscReal :: Tk, tc

  PetscReal :: dw_kg, dw_mol, dw_dp, dw_dt            ! (liq.) water fluid density
  PetscReal :: den_air, dden_air_dp, dden_air_dt      ! air fluid density
  PetscReal :: psat, dpsat_dt                         ! vapor pressure
  PetscReal :: den_ice, dden_ice_dp, dden_ice_dt      ! ice density

  PetscReal :: sl, dsl_dp, dsl_dt
  PetscReal :: sg, dsg_dp, dsg_dt
  PetscReal :: si, dsi_dp, dsi_dt
  PetscReal :: ice_presl, ice_presl_dp, ice_presl_dt

  PetscReal :: kr, dkr_dp, dkr_dt                 ! fluid rel. permissivity
  PetscReal :: vis, dvis_dp, dvis_dt              ! fluid viscosity

  PetscReal :: h, h_dp, h_dt                      ! enthalpy
  PetscReal :: u, u_dp, u_dt                      ! internal energy
  PetscReal :: C_g

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

  ! dofs
  presl  = x(IFDOF1)   ! presl: liq. pressure for liq_fluid movement
  presg  = option%reference_pressure   ! presg: air pressure for air_fluid movement (TODO: the diff is PC)
  Tk     = x(IFDOF2)+TC2TK

  ! isothermal, i.e. hydrology only
  if(option%flow%isothermal) TK = option%reference_temperature+TC2TK
  ! option%flow%isobaric, i.e. thermal process only
  if(option%flow%isobaric) then
    presl = option%reference_pressure
    presg = option%reference_pressure
  endif

  !----------------------------------------------------------------
  ! air pressure at water-air interface
  auxvar%air_pressure = presg
  !
  !----------------------------------------------------------------
  pcmax = abs(characteristic_curves%saturation_function%pcmax)  ! non-negative
  if (presl - auxvar%air_pressure <= -pcmax) then
    presl = auxvar%air_pressure - pcmax
  endif
  auxvar%pc = max(0.d0, auxvar%air_pressure - presl)  ! always non-negative
  auxvar%pres(LIQ_FLUID) = presl
  auxvar%pres(AIR_FLUID) = presg

  pcmin  = 0.d0  ! near-saturation PC zone (hint: a small positive value may be helpful ?)
  if (auxvar%pc <= pcmin) then
    auxvar%pc = 0.d0
    dpresl_dp = 1.d0
  else
    dpresl_dp = 0.d0
  endif
  !
  !----------------------------------------------------------------
  auxvar%tc   = Tk-TC2TK
  auxvar%TK   = Tk-TC2TK


  !----------------------------------------------------------------
  auxvar%kvr = 0.d0

!***************  P/T Bounds ***********************************************

  ! do the trunction of derivatives when assigning temporary values (local) to global variables (auxvar%, global_auxvar%)
  ! for P/T of out of bounds
  dt_trunc = 1.d0
  dp_trunc = 1.d0

  ! general bounds of P/T. We may need to further limit bounds for 3-phase water properties.
  if (Tk<=0.d0) dt_trunc = 0.d0
  Tk = max(0.d0, TK)

  if (presl>=16.54d6) dp_trunc = 0.d0  ! 16.54 MPa is upper limit for using IFC67 EOS.

!***************  Characteristic Curves ***********************************
  
  ! using modules in 'characteristic_curves.F90' to calculate needed variables:
  ! saturations and derivatives for 3 phases
  tc = Tk - TC2TK
  call MpFlowAuxVarComputeCharacteristicCurves(presl, presg, tc,         &
                                           characteristic_curves,        &
                                           sl,  dsl_dp, dsl_dt,          &
                                           si,  dsi_dp, dsi_dt,          &
                                           sg,  dsg_dp, dsg_dt,          &
                                  ice_presl, ice_presl_dp, ice_presl_dt, &
                                           kr,  dkr_dp, dkr_dt,          &
                                           option)
  if(DTRUNC_FLAG) then
    dsl_dp = dsl_dp * dp_trunc
    dsl_dt = dsl_dt * dt_trunc
    dsi_dp = dsi_dp * dp_trunc
    dsi_dt = dsi_dt * dt_trunc
    dsg_dp = dsg_dp * dp_trunc
    dsg_dt = dsg_dt * dt_trunc

    ice_presl_dp = ice_presl_dp * dp_trunc
    ice_presl_dt = ice_presl_dt * dt_trunc
    dkr_dp = dkr_dp * dp_trunc
    dkr_dt = dkr_dt * dt_trunc

  endif

  auxvar%pres_fh2o = ice_presl
  if (auxvar%has_derivs) then
    auxvar%dpres_fh2o(IFDOF1) = ice_presl_dp
    auxvar%dpres_fh2o(IFDOF2) = ice_presl_dt
  endif

!***************  3-phase water properties **********************************************************

  ! ----- Liq. water ---------------------------------------------------------
  auxvar%sat(LIQ_FLUID)           = sl
  if (auxvar%has_derivs) then
    auxvar%D_sat(LIQ_FLUID, IFDOF1) = dsl_dp
    auxvar%D_sat(LIQ_FLUID, IFDOF2) = dsl_dt
  endif

  ! Liq. water density/energy from available EOS (EOS_water.F90) could be crazy when below 0oC,
  ! So a trunction at/below 0.1oC, it may be helpful to avoid those
  tcmin = 0.1d0
  !--------- LIQ. Water Density
  call EOSWaterDensity(max(tcmin, tc),                    &
                       max(auxvar%air_pressure,presl),    &
                        dw_kg, dw_mol, dw_dp, dw_dt,ierr)
  if (DTRUNC_FLAG) dw_dp = dw_dp * dpresl_dp
  if (DTRUNC_FLAG .and. tc<tcmin) dw_dt = 0.d0

  auxvar%den(LIQ_FLUID)           = dw_mol
  if (auxvar%has_derivs) then
    auxvar%D_den(LIQ_FLUID, IFDOF1) = dw_dp
    auxvar%D_den(LIQ_FLUID, IFDOF2) = dw_dt
  endif

  auxvar%den_kg(LIQ_FLUID)           = dw_mol*FMWH2O  ! Not-use 'dw_kg' in case inconsistent 'FMWH2O'
  if (auxvar%has_derivs) then
    auxvar%D_den_kg(LIQ_FLUID, IFDOF1) = dw_dp*FMWH2O
    auxvar%D_den_kg(LIQ_FLUID, IFDOF2) = dw_dt*FMWH2O
  endif

  !----------- LIQ. Water Energy
  call EOSWaterEnthalpy(max(tcmin, tc),                   &
                        max(auxvar%air_pressure,presl),   &
                        h, h_dp, h_dt, ierr)
  if (DTRUNC_FLAG) h_dp = h_dp * dpresl_dp
  if (DTRUNC_FLAG .and. tc<tcmin) h_dt = 0.d0

  auxvar%H(LIQ_FLUID)           = h * 1.d-6         ! J/kmol -> MJ/kmol
  if (auxvar%has_derivs) then
    auxvar%D_H(LIQ_FLUID, IFDOF1) = h_dp * 1.d-6
    auxvar%D_H(LIQ_FLUID, IFDOF2) = h_dt * 1.d-6
  endif

  u = h - presl/dw_mol
  u_dp = h_dp - (dpresl_dp/dw_mol - presl/(dw_mol*dw_mol)*dw_dp)
  u_dt = h_dt + presl/(dw_mol*dw_mol)*dw_dt

  auxvar%U(LIQ_FLUID)           = u * 1.d-6         ! J/kmol -> MJ/kmol
  if (auxvar%has_derivs) then
    auxvar%D_U(LIQ_FLUID, IFDOF1) = u_dp * 1.d-6
    auxvar%D_U(LIQ_FLUID, IFDOF2) = u_dt * 1.d-6
  endif

  !----------- LIQ. Water Viscosity
  ! A note here (F.-M. Yuan: 2017-01-17)
  ! The Viscosity Eq. shows that: temp< ~ -63oC (1atm), 'visl' sharply increases starting from ~ 1.e-2 order.
  ! 'visl' ~ 0. around -133oC, which produces 'inf' for 'kr'
  tcmin = -63.d0

  call EOSWaterSaturationPressure(max(tcmin,tc), psat, dpsat_dt, ierr)
  ! the lowest Tk of 200 for vapor exists in EOS-h2o phase-diagram, but here make it consistent with 'Viscosity'
  if(DTRUNC_FLAG .and. tc<=tcmin) dpsat_dt = 0.d0

  call EOSWaterViscosity(max(tcmin,tc),                &
                         presl,                        &
                         psat, dpsat_dt,               &
                         vis,  dvis_dt, dvis_dp, ierr)
  if(DTRUNC_FLAG .and. tc<=tcmin) dvis_dt = 0.d0
  if(DTRUNC_FLAG) dvis_dp = dvis_dp*dpresl_dp

  auxvar%viscosity = vis
  auxvar%kvr(LIQ_FLUID) = kr/vis
  if (auxvar%has_derivs) then
    auxvar%dkvr(LIQ_FLUID,IFDOF1) = dkr_dp/vis - kr/(vis*vis)*dvis_dp
    auxvar%dkvr(LIQ_FLUID,IFDOF2) = -kr/(vis*vis)*dvis_dt + dkr_dt/vis
  endif

  ! ----- ice water ---------------------------------------------------------

  ! for ice Ih, Tk limit is ~127K (-146oC) in EOS-h2o phase-diagram
  tcmin = -146.d0
  tcmax = 0.d0
  !---------------
  call EOSWaterIceDensity(max(tcmin, min(tcmax,tc)),      &
                          max(auxvar%air_pressure,presl), &
                          den_ice, dden_ice_dt, dden_ice_dp, ierr)
  if(DTRUNC_FLAG) dden_ice_dp = dden_ice_dp * dpresl_dp
  if(DTRUNC_FLAG .and. (tc<=tcmin .or. tc>=tcmax)) dden_ice_dT = 0.d0

  auxvar%den_kg(option%flow%nfluid+1)          = den_ice*FMWH2O  ! '%den' in kmol/m3 ONLY defined for fluids
  if (auxvar%has_derivs) then
    auxvar%D_den_kg(option%flow%nfluid+1,IFDOF1) = dden_ice_dp*FMWH2O
    auxvar%D_den_kg(option%flow%nfluid+1,IFDOF2) = dden_ice_dt*FMWH2O
  endif

  !--------------
  call EOSWaterIceInternalEnergy(max(tcmin, min(tcmax,tc)),      &
                                 max(auxvar%air_pressure,presl), &
                                 PETSC_TRUE,                     &
                                 u, u_dp, u_dt)
  if(DTRUNC_FLAG .and. (tc<=tcmin .or. tc>=tcmax)) u_dt = 0.d0

  auxvar%U(option%flow%nfluid+1)           = u * 1.d-6         ! J/kmol -> MJ/kmol
  if (auxvar%has_derivs) then
    auxvar%D_U(option%flow%nfluid+1,IFDOF1) = u_dp * 1.d-6
    auxvar%D_U(option%flow%nfluid+1,IFDOF2) = u_dt * 1.d-6
  endif

  call EOSWaterIceEnthalpy(max(tcmin, min(tcmax,tc)),      &
                           max(auxvar%air_pressure,presl), &
                           h, h_dp, h_dt, ierr)
  if(DTRUNC_FLAG .and. (tc<=tcmin .or. tc>=tcmax)) h_dt = 0.d0

  auxvar%H(option%flow%nfluid+1)           = h * 1.d-6         ! J/kmol -> MJ/kmol
  if (auxvar%has_derivs) then
    auxvar%D_H(option%flow%nfluid+1, IFDOF1) = h_dp * 1.d-6
    auxvar%D_H(option%flow%nfluid+1,IFDOF2)  = h_dt * 1.d-6
  endif

  ! ----- Air (more checking NEEDED?) ------------------------------------------------------------------------------------------------

  auxvar%sat(AIR_FLUID)           = sg
  if (auxvar%has_derivs) then
    auxvar%D_sat(AIR_FLUID, IFDOF1) = dsg_dp
    auxvar%D_sat(AIR_FLUID, IFDOF2) = dsg_dt
  endif

  ! Calculate the values and derivatives for air density and its internal energy

  ! the lowest Tk of ~200K for vapor exists in EOS-h2o phase-diagram
  tcmin    = -73.d0

  !---------
  den_air     = presg/(IDEAL_GAS_CONSTANT*Tk)*1.d-3        ! in kmol/m3 for all air-mixture
  dden_air_dp = -presg/(IDEAL_GAS_CONSTANT*Tk*Tk)*1.d-3
  dden_air_dt = 1.d0/(IDEAL_GAS_CONSTANT*Tk)*1.d-3
  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    ! the following is a MUST for reducing tiny-time step (NOT sure why).
    dden_air_dp = dden_air_dp*dpresl_dp    ! w.r.t from 'pw' to 'pres(1)' upon soil total saturation
    if (Tk<=tcmin) dden_air_dt = 0.d0
  endif

  auxvar%den(AIR_FLUID)           = den_air
  if (auxvar%has_derivs) then
    auxvar%D_den(AIR_FLUID, IFDOF1) = dden_air_dp
    auxvar%D_den(AIR_FLUID, IFDOF2) = dden_air_dt
  endif

  !-----------------------
  C_g   = C_AIR*FMWAIR          ! in MJ/kmol/K, vapor included in air
  u     = C_g*Tk
  u_dp = 0.d0 ! right ?
  u_dt = C_g
  if(DTRUNC_FLAG) then
    u_dp = u_dp * dp_trunc
    u_dt = u_dt * dt_trunc
  endif
  auxvar%U(AIR_FLUID)           = u
  if (auxvar%has_derivs) then
    auxvar%D_U(AIR_FLUID, IFDOF1) = u_dp
    auxvar%D_U(AIR_FLUID, IFDOF2) = u_dt
  endif

  ! ----- Thermal conductivity (effective, i.e. integrated over fluids and its porous media ) ------

  ! Parameters for computation of effective thermal conductivity
  alpha    = mpflow_parameter%alpha
  alpha_fr = mpflow_parameter%alpha_fr
  Dk_wet   = mpflow_parameter%ckwet
  Dk_dry   = mpflow_parameter%ckdry
  Dk_ice   = mpflow_parameter%ckfrozen

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
  if (auxvar%has_derivs) then
    ! derivative of 'Dk_eff': Dk_eff = Dk_dry + (Dk_wet-Dk_dry)*Ke + (Dk_ice-Dk_dry)*Ke_fr
    auxvar%dDk_eff(IFDOF1) = (Dk_wet-Dk_dry)*dKe_dp + (Dk_ice-Dk_dry)*dKe_fr_dp
    auxvar%dDk_eff(IFDOF2) = (Dk_wet-Dk_dry)*dKe_dt + (Dk_ice-Dk_dry)*dKe_fr_dt
  endif

  ! isothermal, i.e. hydrology only
  if(option%flow%isothermal) then
    ! NO thermal conductivity and zeroing thermal states
    auxvar%Dk_eff = 0.d0
    auxvar%H      = 0.d0
    auxvar%U      = 0.d0

    ! zeroing thermal states derivatives
    if (auxvar%has_derivs) then
      auxvar%dDk_eff = 0.d0
      auxvar%D_H     = 0.d0
      auxvar%D_U     = 0.d0
    endif
  endif

 ! thermal process only
  if(option%flow%isobaric) then
    ! no flow (BUT do allowing water states changing with thermal process, if any)
    auxvar%viscosity = 0.d0
    auxvar%kvr       = 0.d0
    if (auxvar%has_derivs) then
      auxvar%dkvr    = 0.d0
    endif
  endif

end subroutine MpFlowAuxVarCompute

! ************************************************************************** !
subroutine MpFlowAuxVarComputeCharacteristicCurves( presl, presg, tc,       &
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
  use Characteristic_Curves_module
  use EOS_Water_module

  implicit none

  type(option_type) :: option
  PetscReal, intent(in) :: presl    ! unit: Pa (liq water-air interface pressure: -pc+air_pressure)
  PetscReal, intent(in) :: presg    ! unit: Pa (air pressure: air_pressure)
  PetscReal, intent(in) :: tc       ! unit: oC
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
  pc = max(0.d0, presg - presl)   ! always non-negative (0 = saturated)
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
! ************************************************************************** !
! ************************************************************************** !

subroutine DarcyFlowDerivative(mpflow_auxvar_up,    mpflow_auxvar_dn,     &
                               material_auxvar_up,  material_auxvar_dn,   &
                               mpflow_parameter_up, mpflow_parameter_dn,  &
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

  type(mpflow_auxvar_type)    :: mpflow_auxvar_up, mpflow_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(mpflow_parameter_type) :: mpflow_parameter_up, mpflow_parameter_dn
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

  presl_up = mpflow_auxvar_up%pres(LIQ_FLUID)
  presl_dn = mpflow_auxvar_dn%pres(LIQ_FLUID)

  !-------------------------------------------------------------------------
  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)
  if (mpflow_auxvar_up%sat(LIQ_FLUID) < mpflow_parameter_up%sr(LIQ_FLUID)) then
    upweight = 0.d0
  elseif (mpflow_auxvar_dn%sat(LIQ_FLUID) < mpflow_parameter_dn%sr(LIQ_FLUID)) then
    upweight = 1.d0
  end if

  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  v_darcy = 0.d0
  ! Darcy Flow of liquid water
  if (mpflow_auxvar_up%sat(LIQ_FLUID) > mpflow_parameter_up%sr(LIQ_FLUID) .or.  &
      mpflow_auxvar_dn%sat(LIQ_FLUID) > mpflow_parameter_dn%sr(LIQ_FLUID)) then

    gravity = (upweight        *mpflow_auxvar_up%den_kg(LIQ_FLUID)+  &
               (1.d0-upweight) *mpflow_auxvar_dn%den_kg(LIQ_FLUID))  &
              * dist_gravity

#if 0
    dphi = presl_up - presl_dn + gravity
#else
    dphi = mpflow_auxvar_up%pres_fh2o - mpflow_auxvar_dn%pres_fh2o + gravity
#endif

    if (dphi>=0.d0) then
      den  = mpflow_auxvar_up%den(LIQ_FLUID)
      ukvr = mpflow_auxvar_up%kvr(LIQ_FLUID)
    else
      den  = mpflow_auxvar_dn%den(LIQ_FLUID)
      ukvr = mpflow_auxvar_dn%kvr(LIQ_FLUID)
    endif

    if (dabs(ukvr)>floweps) then
      v_darcy = Dq * ukvr * dphi  ! needed for output (m/s)
    endif
    q = v_darcy*area*den          ! m/s * m2 * kmol/m3 = kmol/s

    !-----------------------------------
    if(ifderivative) then
      do idof = 1, option%nflowdof
        dgravity_up(idof) = upweight*dist_gravity* &
          mpflow_auxvar_up%D_den_kg(LIQ_FLUID,idof)

        dgravity_dn(idof) = (1.d0-upweight)*dist_gravity* &
          mpflow_auxvar_dn%D_den_kg(LIQ_FLUID,idof)

#if 0
        ddphi_up(idof) = dgravity_up(idof) + mpflow_auxvar_up%D_pres(LIQ_FLUID, idof)
        ddphi_dn(idof) = dgravity_dn(idof) - mpflow_auxvar_dn%D_pres(LIQ_FLUID, idof)
#else
        ddphi_up(idof) = dgravity_up(idof) + mpflow_auxvar_up%dpres_fh2o(idof)
        ddphi_dn(idof) = dgravity_dn(idof) - mpflow_auxvar_dn%dpres_fh2o(idof)
#endif

        if (dphi>0.D0) then
          dukvr_up(idof) = mpflow_auxvar_up%dkvr(LIQ_FLUID, idof)
          dukvr_dn(idof) = 0.d0

          dden_up(idof)  = mpflow_auxvar_up%D_den(LIQ_FLUID, idof)
          dden_dn(idof)  = 0.d0

        else
          dukvr_up(idof) = 0.d0
          dukvr_dn(idof) = mpflow_auxvar_dn%dkvr(LIQ_FLUID, idof)

          dden_up(idof)  = 0.d0
          dden_dn(idof)  = mpflow_auxvar_up%D_den(LIQ_FLUID, idof)

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

subroutine AdvectionDerivative(mpflow_auxvar_up,   mpflow_auxvar_dn,   &
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
  type(mpflow_auxvar_type):: mpflow_auxvar_up, mpflow_auxvar_dn
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
    uh   = mpflow_auxvar_up%H(LIQ_FLUID)    ! U or H (needs more thinking)
  else
    uh   = mpflow_auxvar_dn%H(LIQ_FLUID)
  endif
  qe     = q*uh

  !-----------------------------------
  if(ifderivative) then
    do idof = 1, option%nflowdof
      if (q>=0.D0) then
        duh_up(idof)   = mpflow_auxvar_up%D_H(LIQ_FLUID, idof)
        duh_dn(idof)   = 0.d0
      else
        duh_up(idof)   = 0.d0
        duh_dn(idof)   = mpflow_auxvar_dn%D_H(LIQ_FLUID, idof)
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

subroutine ConductionDerivative(mpflow_auxvar_up,   mpflow_auxvar_dn,  &
                              material_auxvar_up, material_auxvar_dn,  &
                              area, dist,                              &
                              option,                                  &
                              qe,                                      &
                              dqe_up, dqe_dn, ifderivative)
  !
  ! Computes the derivatives of energy conduction, given flux and its derivatives
  !
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
  !

  use Option_module
  use Material_Aux_class, only : material_auxvar_type
  use Connection_module

  implicit none

  type(option_type) :: option
  type(mpflow_auxvar_type):: mpflow_auxvar_up, mpflow_auxvar_dn
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
  PetscReal, pointer :: dDke_up(:), dDke_dn(:)
  PetscReal :: dDk_up(option%nflowdof), dDk_dn(option%nflowdof)

  !--------------------------------------------------------------------------------------------
  tc_up = mpflow_auxvar_up%tc
  tc_dn = mpflow_auxvar_dn%tc

  qe        = 0.d0
  dqe_up(:) = 0.d0
  dqe_dn(:) = 0.d0

  !-------------------------------------------------------------------------
  call ConnectionCalculateDistances(dist, option%gravity, dd_up, dd_dn, &
                                    dist_gravity, upweight)

  Dke_up = mpflow_auxvar_up%Dk_eff
  Dke_dn = mpflow_auxvar_dn%Dk_eff

  Dke = 0.d0
  if(Dke_up /= 0.d0 .or. Dke_dn /= 0.d0) then
    ! effective thermal conductivity, Dk
    ! 1/Dke = dd_dn/Dke_dn + dd_up/Dke_up
    Dke = (Dke_up * Dke_dn) / (dd_dn*Dke_up + dd_up*Dke_dn)

    qe = Dke*area*(tc_up-tc_dn)
  end if

  if (ifderivative) then
    if(Dke /= 0.d0) then

      dDke_up => mpflow_auxvar_up%dDk_eff
      dDke_dn => mpflow_auxvar_dn%dDk_eff

      do idof = 1, option%nflowdof
       ! 1/Dke = dd_dn/Dke_dn + dd_up/Dke_up

        dDk_up(idof) = Dke*Dke * &
          dd_up*(dDke_up(idof)/Dke_up/Dke_up)  ! d(1/Dk)

        dDk_dn(idof) = Dke*Dke * &
          dd_dn*(dDke_dn(idof)/Dke_dn/Dke_dn)  ! d(1/Dk)

      end do
      ! cond = Dke*area*(mpflow_auxvar_up%temp-mpflow_auxvar_dn%temp)
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
  ! Computes the derivatives of vapor diffusion, given flux and its derivatives
  !
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
  !
  !

  use Option_module
  use Characteristic_Curves_module
  use Connection_module
  use EOS_Water_module
  use Utility_module

  implicit none

  type(mpflow_auxvar_type) :: auxvar_up, auxvar_dn
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


    satg_up = auxvar_up%sat(AIR_FLUID)
    satg_dn = auxvar_dn%sat(AIR_FLUID)
    if ((satg_up > eps) .and. (satg_dn > eps)) then

      p_g = auxvar_up%pres(AIR_FLUID)
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

subroutine MpFlowAuxVarStrip(auxvar)
  ! 
  ! Deallocates a MpFlow auxiliary object's auxvar
  !
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
  !
  ! 

  use Utility_module, only: DeallocateArray

  implicit none

  type(mpflow_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%pres)

  if (auxvar%has_derivs) then
    call DeallocateArray(auxvar%dpres_fh2o)
    call DeallocateArray(auxvar%dkvr)
    call DeallocateArray(auxvar%dDk_eff)
  endif

  call AuxVarFlowEnergyStrip(auxvar)
  
end subroutine MpFlowAuxVarStrip

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
  type(mpflow_auxvar_type), pointer :: auxvar, auxvar_bc, auxvar_ss
  type(mpflow_parameter_type), pointer :: mpflow_parameter
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) then
    do iaux = 1, size(aux%auxvars)
      auxvar => aux%auxvars(iaux)
      call MpFlowAuxVarStrip(auxvar)
    end do
  endif
  nullify(aux%auxvars)

  if (associated(aux%auxvars_bc)) then
    do iaux = 1, size(aux%auxvars_bc)
      auxvar_bc => aux%auxvars_bc(iaux)
      call MpFlowAuxVarStrip(auxvar_bc)
    end do
  endif
  nullify(aux%auxvars_bc)

  if (associated(aux%auxvars_ss)) then
    do iaux = 1, size(aux%auxvars_ss)
      auxvar_ss => aux%auxvars_ss(iaux)
      call MpFlowAuxVarStrip(auxvar_ss)
    end do
  endif
  nullify(aux%auxvars_ss)

  if (associated(aux%mpflow_parameters)) then
    do iaux = 1, size(aux%mpflow_parameters)
      mpflow_parameter => aux%mpflow_parameters(iaux)

      if (associated(mpflow_parameter%sr)) deallocate(mpflow_parameter%sr)
      nullify(mpflow_parameter%sr)

      if (associated(mpflow_parameter%diffusion_coefficient)) &
        deallocate(mpflow_parameter%diffusion_coefficient)
      nullify(mpflow_parameter%diffusion_coefficient)

      if (associated(mpflow_parameter%diffusion_activation_energy)) &
        deallocate(mpflow_parameter%diffusion_activation_energy)
      nullify(mpflow_parameter%diffusion_activation_energy)

    enddo

    deallocate(aux%mpflow_parameters)
    nullify(aux%mpflow_parameters)

  endif
  
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)

  deallocate(aux)
  nullify(aux)  

  end subroutine MpFlowAuxDestroy

end module MpFlow_Aux_module
