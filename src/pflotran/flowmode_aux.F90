module Flowmode_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private 

#include "petsc/finclude/petscsys.h"

  PetscReal, public :: th_itol_scaled_res = 1.d-15
  PetscReal, public :: th_itol_rel_update = UNINITIALIZED_DOUBLE

  type, public :: TH_auxvar_type
    PetscReal :: avgmw
    PetscReal :: transient_por

    PetscReal :: air_pressure    ! unit: Pa ( air pressure of water-air interface. It's reference_pressure if an open system)

    PetscReal :: h   ! enthalpy
    PetscReal :: u   ! internal energy
    PetscReal :: pc
    PetscReal :: vis
    PetscReal :: kvr
    PetscReal :: dsat_dp
    PetscReal :: dsat_dt
    PetscReal :: dden_dp
    PetscReal :: dden_dt
    PetscReal :: dkvr_dp
    PetscReal :: dkvr_dt

    PetscReal :: dh_dp
    PetscReal :: dh_dt
    PetscReal :: du_dp
    PetscReal :: du_dt
    PetscReal :: Dk_eff
    PetscReal :: dDk_eff_dp
    PetscReal :: dDK_eff_dt

    type(th_ice_type), pointer :: ice
  end type TH_auxvar_type

  type, public :: th_ice_type
    PetscReal :: sat_ice
    PetscReal :: sat_air
    PetscReal :: dsat_ice_dp
    PetscReal :: dsat_air_dp
    PetscReal :: dsat_ice_dt
    PetscReal :: dsat_air_dt
    PetscReal :: den_ice
    PetscReal :: dden_ice_dp
    PetscReal :: dden_ice_dt
    PetscReal :: u_ice
    PetscReal :: du_ice_dt
    PetscReal :: du_ice_dp
    PetscReal :: den_air
    PetscReal :: dden_air_dt
    PetscReal :: dden_air_dp
    PetscReal :: u_air
    PetscReal :: du_air_dt
    PetscReal :: du_air_dp
    PetscReal :: molv_air        ! mole fraction of vapor in air
    PetscReal :: dmolv_air_dt
    PetscReal :: dmolv_air_dp
    PetscReal :: pres_fh2o
    PetscReal :: dpres_fh2o_dp
    PetscReal :: dpres_fh2o_dt
  end type th_ice_type
  
  type, public :: TH_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckdry(:) ! Thermal conductivity (dry)
    PetscReal, pointer :: ckwet(:) ! Thermal conductivity (wet)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: ckfrozen(:) ! Thermal conductivity (frozen soil)
    PetscReal, pointer :: alpha_fr(:) ! exponent frozen
    PetscReal, pointer :: sir(:,:)
    PetscReal, pointer :: diffusion_coefficient(:)
    PetscReal, pointer :: diffusion_activation_energy(:)
  end type TH_parameter_type
  
  type, public :: TH_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(TH_parameter_type), pointer :: TH_parameter
    type(TH_auxvar_type), pointer :: auxvars(:)
  end type TH_type

  PetscReal, parameter :: epsilon = 1.d-6

  public :: THAuxCreate, THAuxDestroy, &
            THAuxVarCompute, THAuxVarInit, &
            THAuxVarCopy, THAuxVarDestroy

contains

! ************************************************************************** !

function THAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(TH_type), pointer :: THAuxCreate
  
  type(TH_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  aux%n_zero_rows = 0

  allocate(aux%TH_parameter)
  nullify(aux%TH_parameter%dencpr)
  nullify(aux%TH_parameter%ckdry)
  nullify(aux%TH_parameter%ckwet)
  nullify(aux%TH_parameter%alpha)
  nullify(aux%TH_parameter%ckfrozen)
  nullify(aux%TH_parameter%alpha_fr)
  nullify(aux%TH_parameter%sir)
  nullify(aux%TH_parameter%diffusion_coefficient)
  nullify(aux%TH_parameter%diffusion_activation_energy)
  
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  allocate(aux%TH_parameter%diffusion_coefficient(option%nphase))
  allocate(aux%TH_parameter%diffusion_activation_energy(option%nphase))
  aux%TH_parameter%diffusion_coefficient = 1.d-9
  aux%TH_parameter%diffusion_activation_energy = 0.d0
 
  THAuxCreate => aux
  
end function THAuxCreate

! ************************************************************************** !

subroutine THAuxVarInit(auxvar,option)
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
  
  type(TH_auxvar_type) :: auxvar
  type(option_type) :: option
  
  PetscReal :: uninit_value
  uninit_value     = UNINITIALIZED_DOUBLE

  auxvar%air_pressure = uninit_value
  auxvar%avgmw     = uninit_value
  auxvar%h         = uninit_value
  auxvar%u         = uninit_value
  auxvar%pc        = uninit_value
  auxvar%vis       = uninit_value
  auxvar%kvr       = uninit_value
  auxvar%dsat_dp   = uninit_value
  auxvar%dsat_dt   = uninit_value
  auxvar%dden_dp   = uninit_value
  auxvar%dden_dt   = uninit_value
  auxvar%dkvr_dp   = uninit_value
  auxvar%dkvr_dt   = uninit_value
  auxvar%dh_dp     = uninit_value
  auxvar%dh_dt     = uninit_value
  auxvar%du_dp     = uninit_value
  auxvar%du_dt     = uninit_value    
  auxvar%transient_por = uninit_value
  auxvar%Dk_eff    = uninit_value
  auxvar%dDK_eff_dp= uninit_value
  auxvar%dDK_eff_dt= uninit_value

  ! (TODO - fully 3-phase in all situations)
    allocate(auxvar%ice)
    auxvar%ice%sat_ice       = uninit_value
    auxvar%ice%sat_air       = uninit_value
    auxvar%ice%dsat_ice_dp   = uninit_value
    auxvar%ice%dsat_air_dp   = uninit_value
    auxvar%ice%dsat_ice_dt   = uninit_value
    auxvar%ice%dsat_air_dt   = uninit_value
    auxvar%ice%den_ice       = uninit_value
    auxvar%ice%dden_ice_dp   = uninit_value
    auxvar%ice%dden_ice_dt   = uninit_value
    auxvar%ice%u_ice         = uninit_value
    auxvar%ice%du_ice_dt     = uninit_value
    auxvar%ice%du_ice_dp     = uninit_value
    auxvar%ice%den_air       = uninit_value
    auxvar%ice%dden_air_dt   = uninit_value
    auxvar%ice%dden_air_dp   = uninit_value
    auxvar%ice%u_air         = uninit_value
    auxvar%ice%du_air_dt     = uninit_value
    auxvar%ice%du_air_dp     = uninit_value
    auxvar%ice%molv_air       = uninit_value
    auxvar%ice%dmolv_air_dt   = uninit_value
    auxvar%ice%dmolv_air_dp   = uninit_value
    auxvar%ice%pres_fh2o     = uninit_value
    auxvar%ice%dpres_fh2o_dp = uninit_value
    auxvar%ice%dpres_fh2o_dt = uninit_value
  !

end subroutine THAuxVarInit

! ************************************************************************** !

subroutine THAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 

  use Option_module

  implicit none
  
  type(TH_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%air_pressure = auxvar%air_pressure
  auxvar2%avgmw = auxvar%avgmw
  auxvar2%h = auxvar%h
  auxvar2%u = auxvar%u
  auxvar2%pc = auxvar%pc
  auxvar2%vis = auxvar%vis
  auxvar2%kvr = auxvar%kvr
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dsat_dt = auxvar%dsat_dt
  auxvar2%dden_dp = auxvar%dden_dp
  auxvar2%dden_dt = auxvar%dden_dt
  auxvar2%dkvr_dp = auxvar%dkvr_dp
  auxvar2%dkvr_dt = auxvar%dkvr_dt
  auxvar2%dh_dp = auxvar%dh_dp
  auxvar2%dh_dt = auxvar%dh_dt
  auxvar2%du_dp = auxvar%du_dp
  auxvar2%du_dt = auxvar%du_dt  
  auxvar2%transient_por = auxvar%transient_por
  auxvar2%Dk_eff = auxvar%Dk_eff
  if (associated(auxvar%ice)) then
    auxvar2%ice%sat_ice = auxvar%ice%sat_ice 
    auxvar2%ice%sat_air = auxvar%ice%sat_air
    auxvar2%ice%dsat_ice_dp = auxvar%ice%dsat_ice_dp
    auxvar2%ice%dsat_air_dp = auxvar%ice%dsat_air_dp
    auxvar2%ice%dsat_ice_dt = auxvar%ice%dsat_ice_dt
    auxvar2%ice%dsat_air_dt = auxvar%ice%dsat_air_dt
    auxvar2%ice%den_ice = auxvar%ice%den_ice
    auxvar2%ice%dden_ice_dp = auxvar%ice%dden_ice_dp
    auxvar2%ice%dden_ice_dt = auxvar%ice%dden_ice_dt
    auxvar2%ice%u_ice = auxvar%ice%u_ice
    auxvar2%ice%du_ice_dt = auxvar%ice%du_ice_dt
    auxvar2%ice%du_ice_dp = auxvar%ice%du_ice_dp
    auxvar2%ice%pres_fh2o = auxvar%ice%pres_fh2o
    auxvar2%ice%dpres_fh2o_dp = auxvar%ice%dpres_fh2o_dp
    auxvar2%ice%dpres_fh2o_dt = auxvar%ice%dpres_fh2o_dt
    auxvar2%ice%den_air = auxvar%ice%den_air
    auxvar2%ice%dden_air_dt = auxvar%ice%dden_air_dt
    auxvar2%ice%dden_air_dp = auxvar%ice%dden_air_dp
    auxvar2%ice%u_air = auxvar%ice%u_air
    auxvar2%ice%du_air_dt = auxvar%ice%du_air_dt
    auxvar2%ice%du_air_dp = auxvar%ice%du_air_dp
    auxvar2%ice%molv_air = auxvar%ice%molv_air
    auxvar2%ice%dmolv_air_dt = auxvar%ice%dmolv_air_dt
    auxvar2%ice%dmolv_air_dp = auxvar%ice%dmolv_air_dp
  endif

end subroutine THAuxVarCopy

! ************************************************************************** !

! ************************************************************************** !

subroutine THAuxVarCompute(x, auxvar, global_auxvar, &
                           material_auxvar,          &
                           iphase,                   &
                           characteristic_curves,    &
                           th_parameter, ithrm,      &
                           option)
  ! 
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 11/16/11
  ! 

!sk: Not sure if we need por, perm
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
  PetscReal :: x(option%nflowdof)
  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(TH_parameter_type) :: th_parameter
  PetscInt :: ithrm
  PetscInt :: iphase

  ! local variables
  PetscErrorCode :: ierr

  PetscReal :: pres_l, dpresl_dp
  PetscReal :: temperature

  PetscReal :: dw_kg, dw_mol, dw_dp, dw_dt            ! liq. water density
  PetscReal :: den_ice, dden_ice_dp, dden_ice_dt      ! ice density
  PetscReal :: psat, dpsat_dt                         ! vapor pressure
  PetscReal :: molv_air, dmolv_air_dt, dmolv_air_dp   ! vapor mole fraction in air

  PetscReal :: sl, dsl_dp, dsl_dt
  PetscReal :: sg, dsg_dp, dsg_dt
  PetscReal :: si, dsi_dp, dsi_dt

  PetscReal :: kr, dkr_dp, dkr_dt                ! liq. water permissivity
  PetscReal :: visl, dvis_dp, dvis_dt            ! liq. water viscosity

  PetscReal :: hw, hw_dp, hw_dt                  ! liq. water enthalpy
  PetscReal :: u_ice, du_ice_dt                  ! liq. water internal energy
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

  ! air_pressure at water-air interface
  auxvar%air_pressure = option%reference_pressure ! (TODO: in a closed-system, this is NOT the case)
  !
  global_auxvar%pres(1) = x(1)
  ! Check if the capillary pressure is less than -100MPa, which also limit the lower limit of pres(1).
  pcmax = abs(characteristic_curves%saturation_function%pcmax)  ! non-negative
  if (x(1) - auxvar%air_pressure <= -pcmax) then
    global_auxvar%pres(1) = auxvar%air_pressure - pcmax
  endif
  auxvar%pc = max(0.d0, auxvar%air_pressure - global_auxvar%pres(1))  ! always non-negative
  pres_l = global_auxvar%pres(1)

  !
  global_auxvar%temp    = x(2)
  temperature = global_auxvar%temp

  !----------------------------------------------------------------
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%kvr = 0.d0

!***************  P/T Bounds ***********************************************

  ! do the trunction of derivatives when assigning temporary values (local) to global variables (auxvar%, global_auxvar%)
  ! for P/T of out of bounds
  dt_trunc = 1.d0
  dp_trunc = 1.d0

  ! isothermal, i.e. hydrology only (for simplifying water density calculation)
  if(option%flow%isothermal_eq) temperature = option%reference_temperature

  !--------------------------------------------------------------------

  ! general bounds of P/T. We may need to further limit bounds for 3-phase water properties.
  if (temperature>=100.d0 .or. temperature<=-273.d0) dt_trunc = 0.d0
  temperature    = min(100.d0, max(-273.d0, temperature))

  if (pres_l>=16.54d6) dp_trunc = 0.d0                         ! 16.54 MPa is upper limit for using IFC67 EOS.
  if (pres_l<=auxvar%air_pressure-pcmax) dp_trunc = 0.d0
  

!***************  Characteristic Curves ***********************************
  
  ! using modules in 'characteristic_curves.F90' to calculate needed variables:
  ! saturations and derivatives for 3 phases
  call THAuxVarComputeCharacteristicCurves(pres_l, temperature, auxvar,  &
                                           characteristic_curves,        &
                                           sl,  dsl_dp, dsl_dt,          &
                                           si,  dsi_dp, dsi_dt,          &
                                           sg,  dsg_dp, dsg_dt,          &
                                           auxvar%ice%pres_fh2o,         &
                                           auxvar%ice%dpres_fh2o_dp,     &
                                           auxvar%ice%dpres_fh2o_dt,     &
                                           kr,  dkr_dp, dkr_dt,          &
                                           option)


!***************  3-phase water properties **********************************************************
  auxvar%avgmw = FMWH2O

  ! ----- Liq. water ---------------------------------------------------------

  pcmin  = 0.d0  ! near-saturation PC zone (hint: a small positive value may be helpful ?)
  if (auxvar%pc > pcmin) then
    iphase    = 3
    dpresl_dp = 0.d0    ! assuming no pc depression effect on density
  else
    iphase    = 1       ! saturated
    auxvar%pc = 0.d0
    dpresl_dp = 1.d0
  endif

  global_auxvar%sat(1) = sl
  auxvar%dsat_dp = dsl_dp
  auxvar%dsat_dt = dsl_dt
  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    auxvar%dsat_dp = auxvar%dsat_dp * dp_trunc
    auxvar%dsat_dt = auxvar%dsat_dt * dt_trunc
  endif

  ! F.-M. Yuan (2016-07-10)
  ! Liq. water density/energy from available EOS (EOS_water.F90) could be crazy when below 0oC,
  ! So a trunction at/below 0.1oC, it may be helpful to avoid those
  tcmin = 0.1d0
  !---------
  call EOSWaterDensity(max(tcmin, temperature),               &
                       max(auxvar%air_pressure,pres_l),       &
                        dw_kg, dw_mol, dw_dp, dw_dt,ierr)
  if (DTRUNC_FLAG) dw_dp = dw_dp * dpresl_dp
  if (DTRUNC_FLAG .and. temperature<tcmin) dw_dt = 0.d0

  global_auxvar%den(1)    = dw_mol
  global_auxvar%den_kg(1) = dw_mol * FMWH2O                 ! in case water mole weight not always in a consistent way.
  auxvar%dden_dt = dw_dt * FMWH2O
  auxvar%dden_dp = dw_dp * FMWH2O

  !-----------
  call EOSWaterEnthalpy(max(tcmin, temperature),               &
                        max(auxvar%air_pressure,pres_l),       &
                        hw, hw_dp, hw_dt, ierr)
  if (DTRUNC_FLAG .and. temperature<tcmin) hw_dt = 0.d0
  hw = hw * option%scale         ! J/kmol -> MJ/kmol
  hw_dp = hw_dp * option%scale
  hw_dt = hw_dt * option%scale
  if(DTRUNC_FLAG) hw_dp = hw_dp * dpresl_dp

  auxvar%h = hw
  auxvar%dh_dp = hw_dp
  auxvar%dh_dt = hw_dt

  auxvar%u = hw - pres_l/dw_mol * option%scale
  auxvar%du_dp = hw_dp - (dpresl_dp/dw_mol - pres_l/(dw_mol*dw_mol)*dw_dp)* option%scale
  auxvar%du_dt = hw_dt + pres_l/(dw_mol*dw_mol)*option%scale*dw_dt

  !
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
  auxvar%dkvr_dt = -kr/(visl*visl)*dvis_dt + dkr_dt/visl
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp

  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    auxvar%dsat_dp = auxvar%dsat_dp * dp_trunc
    auxvar%dsat_dt = auxvar%dsat_dt * dt_trunc
    auxvar%dden_dt = auxvar%dden_dt * dt_trunc
    auxvar%dden_dp = auxvar%dden_dp * dp_trunc
    auxvar%dkvr_dt = auxvar%dkvr_dt * dt_trunc
    auxvar%dkvr_dp = auxvar%dkvr_dp * dp_trunc
    auxvar%dh_dt = auxvar%dh_dt * dt_trunc
    auxvar%dh_dp = auxvar%dh_dp * dp_trunc
    auxvar%du_dt = auxvar%du_dt * dt_trunc
    auxvar%du_dp = auxvar%du_dp * dp_trunc
  endif

  ! ----- ice water ---------------------------------------------------------

  auxvar%ice%sat_ice     = si
  auxvar%ice%dsat_ice_dp = dsi_dp
  auxvar%ice%dsat_ice_dt = dsi_dt

  ! for ice Ih, Tk limit is ~127K (-146oC) in EOS-h2o phase-diagram
  tcmin = -146.d0
  tcmax = 0.d0


  !---------------
  call EOSWaterDensityIce(max(tcmin, min(tcmax,temperature)),    &
                          max(auxvar%air_pressure,pres_l), &
                          den_ice, dden_ice_dT, dden_ice_dp, ierr)
  if(DTRUNC_FLAG) dden_ice_dp = dden_ice_dp * dpresl_dp
  if(DTRUNC_FLAG .and. (temperature<=tcmin .or. temperature>=tcmax)) dden_ice_dT = 0.d0

  auxvar%ice%den_ice     = den_ice                ! in kmol/m3
  auxvar%ice%dden_ice_dt = dden_ice_dT
  auxvar%ice%dden_ice_dp = dden_ice_dP

  !--------------
  call EOSWaterInternalEnergyIce(max(tcmin, min(tcmax,temperature)), &
                                 u_ice, du_ice_dT)
  if(DTRUNC_FLAG .and. (temperature<=tcmin .or. temperature>=tcmax)) du_ice_dT = 0.d0

  auxvar%ice%u_ice     = u_ice*1.d-3              !kJ/kmol --> MJ/kmol
  auxvar%ice%du_ice_dt = du_ice_dT*1.d-3          !kJ/kmol/K --> MJ/kmol/K 
  auxvar%ice%du_ice_dp = 0.d0

  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    auxvar%ice%dsat_ice_dp = auxvar%ice%dsat_ice_dp * dp_trunc
    auxvar%ice%dsat_ice_dt = auxvar%ice%dsat_ice_dt * dt_trunc
    auxvar%ice%dden_ice_dp = auxvar%ice%dden_ice_dp * dp_trunc
    auxvar%ice%dden_ice_dt = auxvar%ice%dden_ice_dt * dt_trunc
    auxvar%ice%du_ice_dp   = auxvar%ice%du_ice_dp * dp_trunc
    auxvar%ice%du_ice_dt   = auxvar%ice%du_ice_dt * dt_trunc
  endif

  ! ----- Air (incl. vapor) ---------------------------------------------------------

  auxvar%ice%sat_air = sg
  auxvar%ice%dsat_air_dp = dsg_dp
  auxvar%ice%dsat_air_dt = dsg_dt

  ! Calculate the values and derivatives for air density and its internal energy

  p_g      = max(pres_l, auxvar%air_pressure)
  ! the lowest Tk of ~200K for vapor exists in EOS-h2o phase-diagram
  tcmin    = -73.d0
  Tk_g     = max(tcmin,temperature)+TC2TK

  !---------
  auxvar%ice%den_air     = p_g/(IDEAL_GAS_CONSTANT*tk_g)*1.d-3                ! in kmol/m3 for all air-mixture
  auxvar%ice%dden_air_dt = -p_g/(IDEAL_GAS_CONSTANT*tk_g**2)*1.d-3
  auxvar%ice%dden_air_dp = 1.d0/(IDEAL_GAS_CONSTANT*tk_g)*1.d-3
  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    ! the following is a MUST for reducing tiny-time step (NOT sure why).
    auxvar%ice%dden_air_dp = auxvar%ice%dden_air_dp*dpresl_dp    ! w.r.t from 'pw' to 'pres(1)' upon soil total saturation
    if (tk_g<=tcmin) auxvar%ice%dden_air_dt = 0.d0
  endif

  !----------
  ! NOTE: vapor 'molv_air' is included in 'den_air', 'u_air'
  ! (because 'molv_air', fraction of vapor in air-mixture, going to be as multiplier in all 'air' calculations in 'th.F90')
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

  C_g                    = C_wv*molv_air*FMWH2O + C_a*(1.d0 - molv_air)*FMWAIR          ! in MJ/kmol/K
  auxvar%ice%u_air       = C_g*tk_g                                                     ! in MJ/kmol

  auxvar%ice%molv_air     = molv_air
  auxvar%ice%dmolv_air_dt = dmolv_air_dt
  auxvar%ice%dmolv_air_dp = dmolv_air_dp

  auxvar%ice%du_air_dt   = C_g + (C_wv*dmolv_air_dt*FMWH2O - C_a*dmolv_air_dt*FMWAIR)*tk_g
  auxvar%ice%du_air_dp   = (C_wv*FMWH2O-C_a*FMWAIR)*dmolv_air_dp*tk_g


  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    ! something happened in the following, causing difficulties to re-freeze soils (TODO - checking: 2017-01-18)
    ! (maybe: except for DALL_AMICO model)
    ! 02-07-2017: It may be relevant to soil compressibility.
    !             This also causes LARGE temperature oscillation during F/T - with one layer supper HOT while other supper COLD.

    auxvar%ice%dsat_air_dp = auxvar%ice%dsat_air_dp * dp_trunc
    auxvar%ice%dsat_air_dt = auxvar%ice%dsat_air_dt * dt_trunc

    auxvar%ice%dden_air_dp = auxvar%ice%dden_air_dp * dp_trunc
    auxvar%ice%dden_air_dt = auxvar%ice%dden_air_dt * dt_trunc
    auxvar%ice%dmolv_air_dp = auxvar%ice%dmolv_air_dp * dp_trunc
    auxvar%ice%dmolv_air_dt = auxvar%ice%dmolv_air_dt * dt_trunc
    auxvar%ice%du_air_dp = auxvar%ice%du_air_dp * dp_trunc
    auxvar%ice%du_air_dt = auxvar%ice%du_air_dt * dt_trunc
  endif


  ! ----- Thermal conductivity (effective) ---------------------------------------------------------

  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(ithrm)
  alpha_fr = th_parameter%alpha_fr(ithrm)
  Dk_wet = th_parameter%ckwet(ithrm)
  Dk_dry = th_parameter%ckdry(ithrm)
  Dk_ice = th_parameter%ckfrozen(ithrm)

  !Soil Kersten number
  Ke = (global_auxvar%sat(1) + epsilon)**(alpha)
  Ke_fr = (auxvar%ice%sat_ice + epsilon)**(alpha_fr)

  ! Derivative of Kersten number
  dKe_dp = alpha*(global_auxvar%sat(1)+epsilon)**(alpha-1.d0)*auxvar%dsat_dp
  dKe_dt = alpha*(global_auxvar%sat(1)+epsilon)**(alpha-1.d0)*auxvar%dsat_dt
  dKe_fr_dt = alpha_fr*(auxvar%ice%sat_ice + epsilon)**(alpha_fr - 1.d0)*auxvar%ice%dsat_ice_dt
  dKe_fr_dp = alpha_fr*(auxvar%ice%sat_ice + epsilon)**(alpha_fr-1.d0)*auxvar%ice%dsat_ice_dp

  ! Effective thermal conductivity
  auxvar%Dk_eff = Dk_wet*Ke + Dk_ice*Ke_fr + (1.d0 - Ke - Ke_fr)*Dk_dry
  ! derivative of 'Dk_eff': Dk_eff = Dk_dry + (Dk_wet-Dk_dry)*Ke + (Dk_ice-Dk_dry)*Ke_fr
  auxvar%dDk_eff_dp = (Dk_wet-Dk_dry)*dKe_dp + (Dk_ice-Dk_dry)*dKe_fr_dp
  auxvar%dDk_eff_dt = (Dk_wet-Dk_dry)*dKe_dt + (Dk_ice-Dk_dry)*dKe_fr_dt


  ! isothermal, i.e. hydrology only
  if(option%flow%isothermal_eq) then
    ! NO thermal conductivity and zeroing thermal states
    auxvar%Dk_eff = 0.d0
    auxvar%h      = 0.d0
    auxvar%u      = 0.d0
    auxvar%ice%u_ice       = 0.d0
    auxvar%ice%u_air       = 0.d0

  endif

 ! thermal process only
  if(option%flow%only_thermal_eq) then
    ! no flow
    auxvar%vis = 0.d0
    auxvar%kvr = 0.d0
  endif




end subroutine THAuxVarCompute

! ************************************************************************** !
subroutine THAuxVarComputeCharacteristicCurves( presl,  tc,  auxvar,        &
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
  PetscReal, intent(in) :: tc       ! unit: oC
  type(TH_auxvar_type), intent(in) :: auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal, intent(out) :: sl,  dsl_dpl, dsl_dt
  PetscReal, intent(out) :: si,  dsi_dpl, dsi_dt
  PetscReal, intent(out) :: ice_presl, ice_presl_dpl, ice_presl_dt
  PetscReal, intent(out) :: sg,  dsg_dpl, dsg_dt
  PetscReal, intent(out) :: kr,  dkr_dpl, dkr_dt

  ! local variables

  PetscReal :: pc
  PetscReal :: sli, dsli_dpl, dsli_dt, pcli, dpcli_dsli, presli
  PetscReal :: xplice, dxplice_dpl, dxplice_dt, slx, dslx_dx
  PetscReal :: dkr_dsl

  PetscReal :: rhol, rhol_mol, rhoi, rhoi_mol, rholXsl
  PetscReal :: drhol_dp, drhol_dt, drhoi_dp, drhoi_dt, drholXsl_dp, drholXsl_dt

  PetscErrorCode :: ierr

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
  ! liq. saturation and derivatives without considering ice
  call characteristic_curves%saturation_function%Saturation(pc, sl, dsl_dpl, option)
  ! a note: in CC modules, already w.r.t 'pres_l' by dpc_dpres = -1.d0
  dsl_dt  = 0.d0
  sg      = 1.d0 - sl
  dsg_dpl = -dsl_dpl
  dsg_dt  = -dsl_dt

  ! liq. water PC at ice-liq interface (positive, if no ice it's pc), under 'presl' and 'tc'
  call characteristic_curves%saturation_function%IceCapillaryPressure(presl, tc, &
                                   xplice, dxplice_dpl, dxplice_dt, option) ! w.r.t 'pressure' already in %IceCapillaryPressure()
  ice_presl    = (presl+pc) - xplice     ! 'pres_l+pc' is the liq. water column head (i.e. saturated column+atm.P)
  ice_presl_dpl= -dxplice_dpl            ! dpresl_dpl = 1, dpc_dpl = -1
  ice_presl_dt = -dxplice_dt             ! dpresl_dt = 0, dpc_dt = 0


  !--------------------------------------------------------------------------
  ! if ice module turns on, 3-phase saturation recalculated (liq. and ice)

  ! assuming all 'sli' in liq. at first, and separate them later
  if(auxvar%ice%sat_ice==UNINITIALIZED_DOUBLE) auxvar%ice%sat_ice = 0.d0

  !sli      = auxvar%ice%sat_ice + sl
  sli      = sl
  ! NOTES: sl + 'sat_ice' is correct one, but has issue to cal. ice sat in coming time-step;
  !        only 'sl' has mass-conservation when freezing/thawing likely due to liq/ice water density changing.

  call characteristic_curves%saturation_function%CapillaryPressure(sli, &
                                pcli,dpcli_dsli,option)

  call characteristic_curves%saturation_function%Saturation(pcli, sli, dsli_dpl, option)
  dsli_dt = 0.d0
  presli = (presl+pc) - pcli

  ! update liq. saturation and its derivatives, under ice-adjusted capillary pressure, 'xplice'
  call characteristic_curves%saturation_function%Saturation(xplice, slx, dslx_dx, option)   ! pc on ice ---> sl on ice, but '_dx' is w.r.t. '_dpres'
  sl     = slx
  dsl_dpl= dslx_dx * ice_presl_dpl       ! or, dslx_dpl * dpl_dpc * dxplice_dpl, i.e. here '_dx' is '_dpl', 'xplice' is 'pc', and 'dpl_dpc=-1'
  dsl_dt = dslx_dx * ice_presl_dt


  ! ice satuation and its derivatives
  ! by mass-conservation approach, i.e. simply we have
  ! rhoi*si=delta(rhol*sl), assuming constant volume and porosity
  ! delta(rhol*sl)
  call EOSWaterDensity(tc, max(ice_presl, auxvar%air_pressure),  &    ! under 'ice_presl'
                       rhol, rhol_mol, drhol_dp, drhol_dt,ierr)
  rholXsl     = rhol_mol*sl
  drholXsl_dp = (drhol_dp*ice_presl_dpl)*sl + rholXsl*dsl_dpl         ! w.r.t. (original) 'pres_l' from 'ice_presl'
  drholXsl_dt = (drhol_dt*ice_presl_dt) *sl + rholXsl*dsl_dt

  call EOSWaterDensity(tc, max(presli, auxvar%air_pressure),    &     ! under assumed 'presli'
                       rhol, rhol_mol, drhol_dp, drhol_dt,ierr)

  rholXsl     = (-rholXsl+rhol_mol*sli)                               ! Delta(rhol*sl)
  drholXsl_dp = (-drholXsl_dp)+(drhol_dp*sli +rhol_mol*dsli_dpl)      ! sli ~ presli
  drholXsl_dt = (-drholXsl_dt)+(drhol_dt*sli +rhol_mol*dsli_dt)

  ! updated ice (rhoi*si)
  call EOSWaterDensityIce(tc, max(ice_presl, auxvar%air_pressure), &
                          rhoi_mol, drhoi_dt, drhoi_dp, ierr)
  rhoi = rhoi_mol
  drhoi_dp = drhoi_dp*ice_presl_dpl
  drhoi_dt = drhoi_dt*ice_presl_dt

  si      = rholXsl/rhoi
  dsi_dpl = (rholXsl*drhoi_dp - drholXsl_dp*rhoi)/rhoi/rhoi
  dsi_dt  = (rholXsl*drhoi_dt - drholXsl_dt*rhoi)/rhoi/rhoi

  ! air saturation as difference
  sg      = 1.d0 - sl - si
  dsg_dpl = -dsl_dpl - dsi_dpl
  dsg_dt  = -dsl_dt - dsi_dt


  ! some checking
  if(sl/=sl .or. abs(sl)>huge(sl) .or. &
     si/=si .or. abs(si)>huge(si) .or. &
     sg/=sg .or. abs(sg)>huge(sg)) then
    print *, 'checking Saturation cal. in TH: ', sl, sg, si, presl, pc, ice_presl, xplice
    option%io_buffer = 'TH with characteristic curve: NaN or INF properties'
    call printErrMsg(option)
  endif

  ! Check for bounds on saturations
  if ((sl-1.d0)>1.d-15 .or. sl<-1.d-15) then
    print *, tc, pc, sl, si, sg, sli, xplice
    option%io_buffer = 'TH with ice mode: Liquid Saturation error: >1 or <0'
    call printErrMsg(option)
  endif
  if ((si-1.d0)>1.d-15 .or. si<-1.d-15) then
    print *, tc, pc, sl, si, sg, sli, xplice
    option%io_buffer = 'TH with ice mode: ICE Saturation error:  >1 or <0'
    call printErrMsg(option)
  endif
  if ((sg-1.d0)>1.d-15 .or. sg<-1.d-15) then
    print *, tc, pc, sl, si, sg, sli, xplice
    option%io_buffer = 'TH with ice mode: Air Saturation error:  >1 or <0'
    call printErrMsg(option)
  endif
  if (abs((sl + si + sg)-1.d0)>1.d-15) then
    option%io_buffer = 'TH with ice mode: Saturation not summed to 1 '
    call printErrMsg(option)
  endif

  !--------------------------------------------------------------------------------------------------

  ! (2) relative permissivity of liq. water in multiple-phase mixture
  kr      = 0.d0  !all initialized to zero
  dkr_dsl = 0.d0
  dkr_dt  = 0.d0
  dkr_dpl = 0.d0

  call characteristic_curves%liq_rel_perm_function%RelativePermeability(sl, kr, dkr_dsl, option)
  dkr_dpl= dkr_dsl*dsl_dpl
  dkr_dt = dkr_dsl*dsl_dt

end subroutine THAuxVarComputeCharacteristicCurves

! ************************************************************************** !

subroutine THAuxVarDestroy(auxvar)
  ! 
  ! Deallocates a TH auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  implicit none

  type(TH_auxvar_type) :: auxvar
  
  if (associated(auxvar%ice)) deallocate(auxvar%ice)
  nullify(auxvar%ice)
  
end subroutine THAuxVarDestroy

! ************************************************************************** !

subroutine THAuxDestroy(aux)
  ! 
  ! Deallocates a TH auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  implicit none

  type(TH_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    call THAuxVarDestroy(aux%auxvars(iaux))
  enddo  
  
  if (associated(aux%auxvars)) deallocate(aux%auxvars)
  nullify(aux%auxvars)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
  if (associated(aux%TH_parameter)) then
    if (associated(aux%TH_parameter%diffusion_coefficient)) &
      deallocate(aux%TH_parameter%diffusion_coefficient)
    nullify(aux%TH_parameter%diffusion_coefficient)
    if (associated(aux%TH_parameter%diffusion_activation_energy)) &
      deallocate(aux%TH_parameter%diffusion_activation_energy)
    nullify(aux%TH_parameter%diffusion_activation_energy)
    if (associated(aux%TH_parameter%dencpr)) deallocate(aux%TH_parameter%dencpr)
    nullify(aux%TH_parameter%dencpr)
    if (associated(aux%TH_parameter%ckwet)) deallocate(aux%TH_parameter%ckwet)
    nullify(aux%TH_parameter%ckwet)
    if (associated(aux%TH_parameter%ckdry)) deallocate(aux%TH_parameter%ckdry)
    nullify(aux%TH_parameter%ckdry)
    if (associated(aux%TH_parameter%alpha)) deallocate(aux%TH_parameter%alpha)
    nullify(aux%TH_parameter%alpha)

    if (associated(aux%TH_parameter%ckfrozen)) deallocate(aux%TH_parameter%ckfrozen)
    nullify(aux%TH_parameter%ckfrozen)
    if (associated(aux%TH_parameter%alpha_fr)) deallocate(aux%TH_parameter%alpha_fr)
    nullify(aux%TH_parameter%alpha_fr)

    if (associated(aux%TH_parameter%sir)) deallocate(aux%TH_parameter%sir)
    nullify(aux%TH_parameter%sir)
  endif
  nullify(aux%TH_parameter)
  
  deallocate(aux)
  nullify(aux)  

  end subroutine THAuxDestroy

end module Flowmode_Aux_module
