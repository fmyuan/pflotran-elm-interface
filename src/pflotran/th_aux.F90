module TH_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscReal, public :: th_itol_scaled_res = 1.d-5
  PetscReal, public :: th_itol_rel_update = UNINITIALIZED_DOUBLE

  type, public :: TH_auxvar_type
    PetscReal :: avgmw
    PetscReal :: h
    PetscReal :: u
    PetscReal :: pc
    PetscReal :: vis
!    PetscReal :: dvis_dp
!    PetscReal :: kr
!    PetscReal :: dkr_dp
    PetscReal :: kvr
    PetscReal :: dsat_dp
    PetscReal :: dden_dp
    PetscReal :: dden_dt
    PetscReal :: dkvr_dp
    PetscReal :: dkvr_dt
    PetscReal :: dh_dp
    PetscReal :: dh_dt
    PetscReal :: du_dp
    PetscReal :: du_dt
    PetscReal :: transient_por
    PetscReal :: Dk_eff
    PetscReal :: Ke
    PetscReal :: Ke_fr
    PetscReal :: dKe_dp
    PetscReal :: dKe_dt
    PetscReal :: dKe_fr_dp
    PetscReal :: dKe_fr_dt
    ! ice
    PetscReal :: sat_ice
    PetscReal :: sat_gas
    PetscReal :: dsat_dt
    PetscReal :: dsat_ice_dp
    PetscReal :: dsat_gas_dp
    PetscReal :: dsat_ice_dt
    PetscReal :: dsat_gas_dt
    PetscReal :: den_ice
    PetscReal :: dden_ice_dp
    PetscReal :: dden_ice_dt
    PetscReal :: u_ice
    PetscReal :: du_ice_dt
    PetscReal :: den_gas
    PetscReal :: dden_gas_dt
    PetscReal :: u_gas
    PetscReal :: du_gas_dt
    PetscReal :: mol_gas
    PetscReal :: dmol_gas_dt
    ! For DallAmico model
    PetscReal :: pres_fh2o
    PetscReal :: dpres_fh2o_dp
    PetscReal :: dpres_fh2o_dt
    ! For surface-flow
    PetscBool :: surf_wat
    PetscReal :: P_min
    PetscReal :: P_max
    PetscReal :: coeff_for_cubic_approx(4)
    PetscReal :: coeff_for_deriv_cubic_approx(4)
    PetscReal :: range_for_linear_approx(4)
    PetscReal :: dlinear_slope_dT
    PetscBool :: bcflux_default_scheme

  end type TH_auxvar_type

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
    type(TH_auxvar_type), pointer :: auxvars_bc(:)
    type(TH_auxvar_type), pointer :: auxvars_ss(:)
  end type TH_type

  PetscReal, parameter :: epsilon = 1.d-6

  public :: THAuxCreate, THAuxDestroy, &
            THAuxVarComputeNoFreezing, THAuxVarInit, &
            THAuxVarCopy, THAuxVarDestroy

  public :: THAuxVarComputeFreezing
  public :: THAuxVarComputeFreezing2

contains

! ************************************************************************** !

function THAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

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
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
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

  use Option_module
  use PFLOTRAN_Constants_module, only : UNINITIALIZED_DOUBLE

  implicit none
  
  type(TH_auxvar_type) :: auxvar
  type(option_type) :: option
  
  PetscReal :: uninit_value
  uninit_value     = UNINITIALIZED_DOUBLE

  auxvar%avgmw     = uninit_value
  auxvar%h         = uninit_value
  auxvar%u         = uninit_value
  auxvar%pc        = uninit_value
  !auxvar%kr       = uninit_value
  !auxvar%dkr_dp   = uninit_value
  auxvar%vis       = uninit_value
  !auxvar%dvis_dp  = uninit_value
  auxvar%kvr       = uninit_value
  auxvar%dsat_dp   = uninit_value
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
  auxvar%Ke        = uninit_value
  auxvar%Ke_fr     = uninit_value
  auxvar%dKe_dp    = uninit_value
  auxvar%dKe_dt    = uninit_value
  auxvar%dKe_fr_dp = uninit_value
  auxvar%dKe_fr_dt = uninit_value
  ! NOTE(bja, 2013-12) always initialize ice variables to zero, even if not used!
  auxvar%sat_ice       = uninit_value
  auxvar%sat_gas       = uninit_value
  auxvar%dsat_dt       = uninit_value
  auxvar%dsat_ice_dp   = uninit_value
  auxvar%dsat_gas_dp   = uninit_value
  auxvar%dsat_ice_dt   = uninit_value
  auxvar%dsat_gas_dt   = uninit_value
  auxvar%den_ice       = uninit_value
  auxvar%dden_ice_dp   = uninit_value
  auxvar%dden_ice_dt   = uninit_value
  auxvar%u_ice         = uninit_value
  auxvar%du_ice_dt     = uninit_value
  auxvar%den_gas       = uninit_value
  auxvar%dden_gas_dt   = uninit_value
  auxvar%u_gas         = uninit_value
  auxvar%du_gas_dt     = uninit_value
  auxvar%mol_gas       = uninit_value
  auxvar%dmol_gas_dt   = uninit_value
  auxvar%pres_fh2o     = uninit_value
  auxvar%dpres_fh2o_dp = uninit_value
  auxvar%dpres_fh2o_dt = uninit_value
  auxvar%surf_wat      = PETSC_FALSE
  auxvar%P_min         = uninit_value
  auxvar%P_max         = uninit_value
  auxvar%coeff_for_cubic_approx(:)       = uninit_value
  auxvar%coeff_for_deriv_cubic_approx(:) = uninit_value
  auxvar%range_for_linear_approx(:)      = uninit_value
  auxvar%dlinear_slope_dT                = uninit_value
  auxvar%bcflux_default_scheme           = PETSC_FALSE

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

! auxvar2%pres = auxvar%pres
! auxvar2%temp = auxvar%temp
! auxvar2%den = auxvar%den
! auxvar2%den_kg = auxvar%den_kg
    
  auxvar2%avgmw = auxvar%avgmw
  auxvar2%h = auxvar%h
  auxvar2%u = auxvar%u
  auxvar2%pc = auxvar%pc
!  auxvar2%kr = auxvar%kr
!  auxvar2%dkr_dp = auxvar%dkr_dp
  auxvar2%vis = auxvar%vis
!  auxvar2%dvis_dp = auxvar%dvis_dp
  auxvar2%kvr = auxvar%kvr
  auxvar2%dsat_dp = auxvar%dsat_dp
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
  auxvar2%Ke = auxvar%Ke
  auxvar2%Ke_fr = auxvar%Ke_fr
  auxvar2%dKe_dp = auxvar%dKe_dp
  auxvar2%dKe_fr_dp = auxvar%dKe_fr_dp
  auxvar2%dKe_dt = auxvar%dKe_dt
  auxvar2%dKe_fr_dt = auxvar%dKe_fr_dt
  if (option%use_th_freezing) then
     auxvar2%sat_ice = auxvar%sat_ice 
     auxvar2%sat_gas = auxvar%sat_gas
     auxvar2%dsat_dt = auxvar%dsat_dt
     auxvar2%dsat_ice_dp = auxvar%dsat_ice_dp
     auxvar2%dsat_gas_dp = auxvar%dsat_gas_dp
     auxvar2%dsat_ice_dt = auxvar%dsat_ice_dt
     auxvar2%dsat_gas_dt = auxvar%dsat_gas_dt
     auxvar2%den_ice = auxvar%den_ice
     auxvar2%dden_ice_dp = auxvar%dden_ice_dp
     auxvar2%dden_ice_dt = auxvar%dden_ice_dt
     auxvar2%u_ice = auxvar%u_ice
     auxvar2%du_ice_dt = auxvar%du_ice_dt
     auxvar2%pres_fh2o = auxvar%pres_fh2o
     auxvar2%dpres_fh2o_dp = auxvar%dpres_fh2o_dp
     auxvar2%dpres_fh2o_dt = auxvar%dpres_fh2o_dt
     auxvar2%den_gas = auxvar%den_gas
     auxvar2%dden_gas_dt = auxvar%dden_gas_dt
     auxvar2%u_gas = auxvar%u_gas
     auxvar2%du_gas_dt = auxvar%du_gas_dt
     auxvar2%mol_gas = auxvar%mol_gas
     auxvar2%dmol_gas_dt = auxvar%dmol_gas_dt
  endif

  auxvar2%surf_wat = auxvar%surf_wat
  auxvar2%P_min = auxvar%P_min
  auxvar2%P_max = auxvar%P_max
  auxvar2%coeff_for_cubic_approx(:) = auxvar%coeff_for_cubic_approx(:)
  auxvar2%coeff_for_deriv_cubic_approx(:) = auxvar%coeff_for_deriv_cubic_approx(:)
  auxvar2%range_for_linear_approx(:) = auxvar%range_for_linear_approx(:)
  auxvar2%dlinear_slope_dT = auxvar%dlinear_slope_dT
  auxvar2%bcflux_default_scheme = auxvar%bcflux_default_scheme

end subroutine THAuxVarCopy

! ************************************************************************** !

subroutine THAuxVarComputeNoFreezing(x,auxvar,global_auxvar, &
                                     material_auxvar, &
                                     iphase,saturation_function, &
                                     th_parameter, ithrm, &
                                     option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: ???
  ! Date: 02/22/08
  ! 

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Saturation_Function_module  
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: iphase
  type(TH_parameter_type) :: th_parameter
  PetscInt :: ithrm
  class(material_auxvar_type) :: material_auxvar

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  PetscReal :: Ke
  PetscReal :: alpha
  PetscReal :: Dk
  PetscReal :: Dk_dry

! auxvar%den = 0.d0
! auxvar%den_kg = 0.d0
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%kvr = 0.d0
  kr = 0.d0
 
! auxvar%pres = x(1)  
! auxvar%temp = x(2)
  global_auxvar%pres = x(1)  
  global_auxvar%temp = x(2)
 
! auxvar%pc = option%reference_pressure - auxvar%pres
  auxvar%pc = option%reference_pressure - global_auxvar%pres(1)

!***************  Liquid phase properties **************************
  auxvar%avgmw = FMWH2O

  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
!  if (auxvar%pc > 0.d0) then
  if (auxvar%pc > 1.d0) then
    iphase = 3
    call SaturationFunctionCompute(auxvar%pc,global_auxvar%sat(1), &
                                   kr,ds_dp,dkr_dp, &
                                   saturation_function, &
                                   material_auxvar%porosity, &
                                   material_auxvar%permeability(perm_xx_index), &
                                   option)
    dpw_dp = 0.d0
  else
    iphase = 1
    auxvar%pc = 0.d0
    global_auxvar%sat(1) = 1.d0  
    kr = 1.d0    
!   pw = auxvar%pres
    pw = global_auxvar%pres(1)
    dpw_dp = 1.d0
  endif  

!  call wateos_noderiv(option%temp,pw,dw_kg,dw_mol,hw,option%scale,ierr)
  call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol,dw_dp,dw_dt,ierr)
  call EOSWaterEnthalpy(global_auxvar%temp,pw,hw,hw_dp,hw_dt,ierr)
  ! J/kmol -> whatever units
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dt = hw_dt * option%scale
  
! may need to compute dpsat_dt to pass to VISW
  call EOSWaterSaturationPressure(global_auxvar%temp,sat_pressure,dpsat_dt,ierr)
  
!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
  call EOSWaterViscosity(global_auxvar%temp,pw,sat_pressure,dpsat_dt,visl, &
                         dvis_dt,dvis_dp,dvis_dpsat,ierr)
  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif

! auxvar%den = dw_mol
! auxvar%den_kg = dw_kg
  global_auxvar%den = dw_mol
  global_auxvar%den_kg = dw_kg
  
  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_mol * option%scale
  auxvar%kvr = kr/visl
  
  auxvar%vis = visl
!  auxvar%dvis_dp = dvis_dp
!  auxvar%kr = kr
!  auxvar%dkr_dp = dkr_dp
  auxvar%dsat_dp = ds_dp
  auxvar%dden_dt = dw_dt

  auxvar%dden_dp = dw_dp
  
!geh: contribution of dvis_dpsat is now added in EOSWaterViscosity
!  auxvar%dkvr_dt = -kr/(visl*visl)*(dvis_dt+dvis_dpsat*dpsat_dt)
  auxvar%dkvr_dt = -kr/(visl*visl)*dvis_dt
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  if (iphase < 3) then !kludge since pw is constant in the unsat zone
    auxvar%dh_dp = hw_dp
    auxvar%du_dp = hw_dp - (dpw_dp/dw_mol-pw/(dw_mol*dw_mol)*dw_dp)*option%scale
  else
    auxvar%dh_dp = 0.d0
    auxvar%du_dp = 0.d0
  endif

  auxvar%dh_dt = hw_dt
  auxvar%du_dt = hw_dt + pw/(dw_mol*dw_mol)*option%scale*dw_dt
  
  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(ithrm)
  Dk = th_parameter%ckwet(ithrm)
  Dk_dry = th_parameter%ckdry(ithrm)

  !unfrozen soil Kersten number
  Ke = (global_auxvar%sat(1) + epsilon)**(alpha)
  auxvar%Ke = Ke

  ! Effective thermal conductivity
  auxvar%Dk_eff = Dk_dry + (Dk - Dk_dry)*Ke

  ! Derivative of soil Kersten number
  auxvar%dKe_dp = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)*auxvar%dsat_dp
  auxvar%dKe_dt = 0.d0

end subroutine THAuxVarComputeNoFreezing

! ************************************************************************** !

subroutine THAuxVarComputeFreezing(x, auxvar, global_auxvar, &
                                   material_auxvar, &
                                   iphase, &
                                   saturation_function, &
                                   th_parameter, ithrm, &
                                   option)
  ! 
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 11/16/11
  ! 

!sk: Not sure if we need por, perm

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Saturation_Function_module  
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(TH_parameter_type) :: th_parameter
  PetscInt :: ithrm
  PetscInt :: iphase

  PetscErrorCode :: ierr
  PetscReal :: pw, dw_kg, dw_mol, hw, sat_pressure, visl
  PetscReal :: kr, ds_dp, dkr_dp, dkr_dt
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  PetscReal :: ice_saturation, gas_saturation
  PetscReal :: dsl_temp
  PetscReal :: dsg_pl, dsg_temp
  PetscReal :: dsi_pl, dsi_temp
  PetscReal :: den_ice, dden_ice_dT, dden_ice_dP
  PetscReal :: u_ice, du_ice_dT
  PetscBool :: out_of_table_flag
  PetscReal :: p_th

  PetscReal :: p_g
  PetscReal :: p_sat
  PetscReal :: mol_g
  PetscReal :: C_g
  PetscReal :: dmolg_dt
  PetscReal, parameter :: C_a = 1.86d-3 ! in MJ/kg/K at 300K
  PetscReal, parameter :: C_wv = 1.005d-3 ! in MJ/kg/K

  PetscReal :: Ke
  PetscReal :: Ke_fr
  PetscReal :: alpha
  PetscReal :: alpha_fr
  PetscReal :: Dk
  PetscReal :: Dk_dry
  PetscReal :: Dk_ice

  out_of_table_flag = PETSC_FALSE
 
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%kvr = 0.d0
   
  global_auxvar%pres = x(1)  
  global_auxvar%temp = x(2)
  
  ! Check if the capillary pressure is less than -100MPa
  
  if (global_auxvar%pres(1) - option%reference_pressure < -1.d8 + 1.d0) then
    global_auxvar%pres(1) = -1.d8 + option%reference_pressure + 1.d0
  endif

 
  auxvar%pc = option%reference_pressure - global_auxvar%pres(1)

!***************  Liquid phase properties **************************
  auxvar%avgmw = FMWH2O

  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
  if (auxvar%pc > 1.d0) then
    iphase = 3
    dpw_dp = 0.d0
  else
    iphase = 1
    auxvar%pc = 0.d0
    pw = global_auxvar%pres(1)
    dpw_dp = 1.d0
  endif  
  
  call CapillaryPressureThreshold(saturation_function,p_th,option)

  select case (option%ice_model)
    case (PAINTER_EXPLICIT)
      ! Model from Painter, Comp. Geosci. (2011)
      call SatFuncComputeIcePExplicit(global_auxvar%pres(1), & 
                                      global_auxvar%temp, ice_saturation, &
                                      global_auxvar%sat(1), gas_saturation, &
                                      kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                      dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                      saturation_function, p_th, option)    
    case (PAINTER_KARRA_IMPLICIT)
      ! Implicit model from Painter & Karra, VJZ (2013)
      call SatFuncComputeIcePKImplicit(global_auxvar%pres(1), & 
                                       global_auxvar%temp, ice_saturation, &
                                       global_auxvar%sat(1), gas_saturation, &
                                       kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                       dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                       saturation_function, p_th, option)    
    case (PAINTER_KARRA_EXPLICIT)
      ! Explicit model from Painter & Karra, VJZ (2013)
      call SatFuncComputeIcePKExplicit(global_auxvar%pres(1), & 
                                       global_auxvar%temp, ice_saturation, &
                                       global_auxvar%sat(1), gas_saturation, &
                                       kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                       dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                       saturation_function, p_th, option) 
    case (DALL_AMICO)
      ! Model from Dall'Amico (2010) and Dall' Amico et al. (2011)
      call SatFuncComputeIceDallAmico(global_auxvar%pres(1), &
                                      global_auxvar%temp, &
                                      auxvar%pres_fh2o, &
                                      auxvar%dpres_fh2o_dp, &
                                      auxvar%dpres_fh2o_dt, &
                                      ice_saturation, &
                                      global_auxvar%sat(1), gas_saturation, &
                                      kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                      dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                      saturation_function, option)
    case (PAINTER_KARRA_EXPLICIT_NOCRYO)
      ! Explicit model from Painter & Karra, VJZ (2013) and removed cryosuction
      call SatFuncComputeIcePKExplicitNoCryo(global_auxvar%pres(1), & 
                                       global_auxvar%temp, ice_saturation, &
                                       global_auxvar%sat(1), gas_saturation, &
                                       kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                       dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                       saturation_function, p_th, option) 
    case default
      option%io_buffer = 'THCAuxVarComputeIce: Ice model not recognized.'
      call printErrMsg(option)
  end select

  call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol,dw_dp,dw_dt,ierr)
  call EOSWaterEnthalpy(global_auxvar%temp,pw,hw,hw_dp,hw_dt,ierr)
  ! J/kmol -> MJ/kmol
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dt = hw_dt * option%scale
                         
  call EOSWaterSaturationPressure(global_auxvar%temp, sat_pressure, &
                                  dpsat_dt, ierr)
  call EOSWaterViscosity(global_auxvar%temp, pw, sat_pressure, dpsat_dt, &
                         visl, dvis_dt,dvis_dp, dvis_dpsat, ierr)

  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif

  global_auxvar%den = dw_mol
  global_auxvar%den_kg = dw_kg
  
  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_mol * option%scale
  auxvar%kvr = kr/visl
  auxvar%vis = visl
  auxvar%dsat_dp = ds_dp
  auxvar%dden_dt = dw_dt
  auxvar%dden_dp = dw_dp
!geh: contribution of dvis_dpsat is now added in EOSWaterViscosity  
!  auxvar%dkvr_dt = -kr/(visl*visl)*(dvis_dt + dvis_dpsat*dpsat_dt) + dkr_dt/visl
  auxvar%dkvr_dt = -kr/(visl*visl)*dvis_dt + dkr_dt/visl
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  auxvar%dh_dp = hw_dp
  auxvar%du_dp = hw_dp - (dpw_dp/dw_mol - pw/(dw_mol*dw_mol)*dw_dp)* &
                  option%scale
  auxvar%dh_dt = hw_dt
  auxvar%du_dt = hw_dt + pw/(dw_mol*dw_mol)*option%scale*dw_dt

  auxvar%sat_ice = ice_saturation
  auxvar%sat_gas = gas_saturation
  auxvar%dsat_dt = dsl_temp
  auxvar%dsat_ice_dp = dsi_pl
  auxvar%dsat_gas_dp = dsg_pl
  auxvar%dsat_ice_dt = dsi_temp
  auxvar%dsat_gas_dt = dsg_temp
  
  ! Calculate the density, internal energy and derivatives for ice
  call EOSWaterDensityIce(global_auxvar%temp, global_auxvar%pres(1), &
                          den_ice, dden_ice_dT, dden_ice_dP, ierr)

  call EOSWaterInternalEnergyIce(global_auxvar%temp, u_ice, du_ice_dT)

  auxvar%den_ice = den_ice
  auxvar%dden_ice_dt = dden_ice_dT
  auxvar%dden_ice_dp = dden_ice_dP
  auxvar%u_ice = u_ice*1.d-3                  !kJ/kmol --> MJ/kmol
  auxvar%du_ice_dt = du_ice_dT*1.d-3          !kJ/kmol/K --> MJ/kmol/K 

  ! Calculate the values and derivatives for density and internal energy
  call EOSWaterSaturationPressure(global_auxvar%temp, p_sat, ierr)

  p_g            = option%reference_pressure
  auxvar%den_gas = p_g/(IDEAL_GAS_CONSTANT*(global_auxvar%temp + 273.15d0))*1.d-3 !in kmol/m3
  mol_g          = p_sat/p_g
  C_g            = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR ! in MJ/kmol/K
  auxvar%u_gas   = C_g*(global_auxvar%temp + 273.15d0)           ! in MJ/kmol
  auxvar%mol_gas = mol_g

  auxvar%dden_gas_dt = - p_g/(IDEAL_GAS_CONSTANT*(global_auxvar%temp + 273.15d0)**2)*1.d-3
  dmolg_dt           = dpsat_dt/p_g
  auxvar%du_gas_dt   = C_g + (C_wv*dmolg_dt*FMWH2O - C_a*dmolg_dt*FMWAIR)* &
                       (global_auxvar%temp + 273.15d0)
  auxvar%dmol_gas_dt = dmolg_dt

  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(ithrm)
  alpha_fr = th_parameter%alpha_fr(ithrm)
  Dk = th_parameter%ckwet(ithrm)
  Dk_dry = th_parameter%ckdry(ithrm)
  Dk_ice = th_parameter%ckfrozen(ithrm)

  !Soil Kersten number
  Ke = (global_auxvar%sat(1) + epsilon)**(alpha)
  Ke_fr = (auxvar%sat_ice + epsilon)**(alpha_fr)
  auxvar%Ke = Ke
  auxvar%Ke_fr = Ke_fr

  ! Effective thermal conductivity
  auxvar%Dk_eff = Dk*Ke + Dk_ice*Ke_fr + (1.d0 - Ke - Ke_fr)*Dk_dry

  ! Derivative of Kersten number
  auxvar%dKe_dp = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)*auxvar%dsat_dp
  auxvar%dKe_dt = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)*auxvar%dsat_dt
  auxvar%dKe_fr_dt = alpha_fr*(auxvar%sat_ice + epsilon)**(alpha_fr - 1.d0)*auxvar%dsat_ice_dt
  auxvar%dKe_fr_dp = alpha_fr*(auxvar%sat_ice + epsilon)**(alpha_fr - 1.d0)*auxvar%dsat_ice_dp

  if (option%ice_model == DALL_AMICO) then
    auxvar%den_ice = dw_mol
    auxvar%dden_ice_dt = auxvar%dden_dt
    auxvar%dden_ice_dp = auxvar%dden_dp
!    auxvar%u_ice = auxvar%u  ! commented out by S.Karra 06/02/14. setting
!    internal energy of ice and water might not be correct.
!    auxvar%du_ice_dt = auxvar%du_dt

    auxvar%sat_gas       = 0.d0
    auxvar%dsat_gas_dp   = 0.d0
    auxvar%dsat_gas_dt   = 0.d0
    auxvar%den_gas       = 0.d0
    auxvar%dden_gas_dt   = 0.d0
    auxvar%u_gas         = 0.d0
    auxvar%du_gas_dt     = 0.d0
    auxvar%mol_gas       = 0.d0
    auxvar%dmol_gas_dt   = 0.d0
  endif

end subroutine THAuxVarComputeFreezing

! ************************************************************************** !

subroutine THAuxVarComputeFreezing2(x, auxvar, global_auxvar, &
                                   material_auxvar, &
                                   isat_state, &
                                   characteristic_curves, &
                                   th_parameter, ithrm, &
                                   option)
  !
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present
  !
  ! Revised by fengming Yuan @03-08-2016/CCSI-ONRL (renamed as version 2)
  ! NOTE: (1) ice_model 'PAINTER_KARRA_EXPLICIT'
  !       (2) VG saturation_function, from 'Characteristic_Curves_module'
  !       (3) smoothed MULEM permissivity function, from 'Characteristci_Curves_module' as well

  use Option_module
  use Global_Aux_module

  use EOS_Water_module
  use Characteristic_Curves_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscReal,intent(in)        :: x(option%nflowdof)
  type(TH_auxvar_type)        :: auxvar
  type(global_auxvar_type)    :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  type(TH_parameter_type)           :: th_parameter
  PetscInt, intent(in)  :: ithrm
  PetscInt, intent(inout) :: isat_state

  ! locals
  PetscErrorCode :: ierr

  ! liq. water
  PetscReal :: pw, dpw_dp                            ! pressure
  PetscReal :: denw_kg, denw_mol, ddenw_dp, ddenw_dt ! density
  PetscReal :: hw, dhw_dp, dhw_dt                    ! enthalpy (h)
  PetscReal :: u, du_dp, du_dt                       ! internal energy (u)
  PetscReal :: kr, dkr_dse, dkr_dp, dkr_dt           ! permissivity (rel.)
  PetscReal :: visl, dvisl_dp, dvisl_dt              ! viscosity

  ! ice
  PetscReal :: den_ice, dden_ice_dt, dden_ice_dp
  PetscReal :: u_ice, du_ice_dt

  ! air(gas)
  PetscReal :: p_g
  PetscReal :: mol_g
  PetscReal :: C_g
  PetscReal :: dmolg_dt
  PetscReal :: krg, dkrg_dp, dkrg_dt                 ! permissivity (rel.)

  ! saturations (multiple-phase mixture in por)
  PetscReal :: sl, dsl_dp, dsl_dt
  PetscReal :: sg, dsg_dp, dsg_dt
  PetscReal :: si, dsi_dp, dsi_dt

  ! thermal properties.
  PetscReal :: Ke, Ke_fr         ! Kersten number: liq, ice
  PetscReal :: alpha, alpha_fr   !
  PetscReal :: dKe_dp, dKe_dt, dKe_fr_dp, dKe_fr_dt

  ! misc.
  PetscReal :: psat
  PetscReal :: dpsat_dt, dvisl_dpsat

  ! some thermal properties
  PetscReal :: Dk
  PetscReal :: Dk_dry
  PetscReal :: Dk_ice
  PetscReal, parameter :: C_a = 1.86d-3   ! air heat capacity: in MJ/kg/K at 300K
  PetscReal, parameter :: C_wv = 1.005d-3 ! liq. water heat capacity: in MJ/kg/K

  PetscReal, parameter :: beta = 2.20d0          ! dimensionless -- ratio of soil ice surf. tension
  PetscReal, parameter :: Tf_0 = 273.15d0        ! freezing-point at standard pressure: in K
  PetscReal, parameter :: L_f = 3.34d5           ! in J/kg
  PetscReal :: Tf, tftheta, sli, dsli_dp, xplice, dxplice_dp, dxplice_dt, slx, dslx_dx, rho_l

  PetscReal :: tc, Tk, pres_l, pc
  PetscInt  :: i,j
  PetscBool :: saturated
  !----------------------------------
  pres_l = x(TH_PRESSURE_DOF)
  tc = x(TH_TEMPERATURE_DOF)

#if 0
  ! The following 2 do-loops are for testing 3-phase saturations and derivatives regarding to P and T
   write(option%myrank+100, *) "pc, pw, xplice, tc, Tf, denw_mol, ddenw_dp, ddenw_dt, den_ice, " // &
    "dden_ice_dp, dden_ice_dt, sl, dsl_dp, dsl_dt, si, dsi_dp, dsi_dt, sli, dsli_dp, slx, dslx_dx"
  do i=-1,9
  do j=-300,100
  pres_l = -(10**i-0.1d0)+option%reference_pressure
  tc = j*0.1d0
#endif

  global_auxvar%pres   = pres_l
  global_auxvar%temp   = tc
  global_auxvar%sat    = 0.d0    ! will be updated below
  global_auxvar%den    = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%h   = 0.d0
  auxvar%u   = 0.d0
  auxvar%kvr = 0.d0
  auxvar%avgmw = FMWH2O

  !
  pc = max(0.d0, option%reference_pressure - global_auxvar%pres(1))   ! always non-negative (0 = saturated)
  if (pc > abs(characteristic_curves%saturation_function%pcmax)) then
    pc = characteristic_curves%saturation_function%pcmax
    global_auxvar%pres(1) = -pc + option%reference_pressure
  endif
  auxvar%pc =  pc

  !***************  phase-relevant properties and derivatives, regarding to P/T *********************
  if (auxvar%pc > 0.d0) then
    isat_state = 3      ! unsaturated
    saturated  = PETSC_FALSE
  else
    isat_state = 1      ! saturated
    saturated  = PETSC_TRUE
  endif

  ! assuming 3-phase H2O is under ONE common water head pressure ('pw') and temperature ('tc')
  pw = option%reference_pressure
  dpw_dp = 0.d0
  if (isat_state == 1) then
    pw = global_auxvar%pres(1)
    dpw_dp = 1.d0
  endif

  ! (1) densities (calculated at first, so that may be used later)
  ! liq. water
  call EOSWaterDensity   (tc, pw, &
                          denw_kg, denw_mol, ddenw_dp, ddenw_dt, ierr)
  global_auxvar%den    = denw_mol
  global_auxvar%den_kg = denw_kg

  auxvar%dden_dt       = ddenw_dt
  if (.not.saturated) ddenw_dp = 0.d0 !kludge since pw is constant in the unsat zone
  auxvar%dden_dp       = ddenw_dp

  ! ice water, if any
  call EOSWaterDensityIce(tc, pw, &
                          den_ice, dden_ice_dt, dden_ice_dp, ierr)
  auxvar%den_ice = den_ice
  auxvar%dden_ice_dt = dden_ice_dt
  auxvar%dden_ice_dp = dden_ice_dp

  ! (2) saturations of multiple-phase water
  sl     = 0.d0  !all init to zero
  dsl_dp = 0.d0
  dsl_dt = 0.d0
  si     = 0.d0
  dsi_dp = 0.d0
  dsi_dt = 0.d0
  sg     = 0.d0
  dsg_dp = 0.d0
  dsg_dt = 0.d0

  select type (sf => characteristic_curves%saturation_function)
    ! class is(sat_func_BC_type)
       ! not-yet
    class is(sat_func_VG_type)
      call sf%Saturation(pc, sli, dsli_dp, option)

      ! default liq. saturation and derivatives
      sl = sli
      dsl_dp = dsli_dp
      dsl_dt = 0.d0

      ! if ice module turns on, 2-phase saturation recalculated (liq. and ice) under common 'pw' and 't'
      if (option%use_th_freezing) then
        select case (option%ice_model)
          case (PAINTER_KARRA_EXPLICIT)
           ! explicit model from Painter & Karra, VJZ (2014)
            rho_l = denw_mol*FMWH2O    ! kg/m3: kmol/m3*kg/kmol
            Tk = tc + Tf_0                                       ! convert to K
            Tf = Tf_0 - 1.d0/beta*Tf_0/L_f/rho_l*pc     ! not used in 'xplice' ?
            tftheta = (Tk-Tf_0)/Tf_0

            xplice = -beta*tftheta*L_f*rho_l
            dxplice_dt = -beta*L_f*                  &
                          ( rho_l*(1.d0/Tf_0)          &
                           +tftheta*(FMWH2O*ddenw_dt))
            dxplice_dp = -beta*L_f*tftheta*          &
                          (FMWH2O*ddenw_dp)

            call sf%Saturation(xplice, slx, dslx_dx, option)

            if (slx<sli .and. Tk<=Tf) then
              sl = slx
              si = 1.d0 - sl/sli

              dsl_dt = dslx_dx*dxplice_dt
              dsl_dp = dslx_dx*dxplice_dp

              dsi_dt = -1.d0/sli*dsl_dt     ! dsli_dt = 0 (see above)
              dsi_dp = -1.d0*(-sl*dsli_dp+sli*dsl_dp)/(sli**2)   ! slx = sl

            endif

          case default
            option%io_buffer = 'THAuxVarComputeFreezing2: Ice model now only support: PAINTER_KARRA_EXPLICIT.'
            call printErrMsg(option)
        end select ! select case (option%ice_model)

      endif ! 'option%use_th_freezing'

    class default
      option%io_buffer = 'Only van Genuchten supported with ice model option'
      call printErrMsg(option)
  end select ! select type (sf=>characteristic_curves%saturation_function)
  !
  sg = 1.d0 - sl - si
  dsg_dp = -dsl_dp - dsi_dp
  dsg_dt = -dsl_dt - dsi_dt

#if 0
  ! The following print will produce a txt file, which can be imported to other tools for checking
  write(option%myrank+100, *) pc, pw, xplice, tc, Tf,                &
    denw_mol, ddenw_dp, ddenw_dt, den_ice, dden_ice_dp, dden_ice_dt, &
    sl, dsl_dp, dsl_dt, si, dsi_dp, dsi_dt,                          &
    sli, dsli_dp, slx, dslx_dx
  ! end of 2-do-loops for testing 3-phase water saturation and derivatives
  enddo
  enddo
  option%io_buffer = 'testing freezing-thawing curve'
  call printErrMsg(option)
#endif

  ! Check for bounds on saturations
  if (sl > 1.d0 .or. sl<0.d0) then
    option%io_buffer = 'TH with ice mode: Liquid Saturation error: >1 or <0'
    call printErrMsg(option)
  endif
  if (si > 1.d0 .or. si<0.d0) then
    option%io_buffer = 'TH with ice mode: ICE Saturation error:  >1 or <0'
    call printErrMsg(option)
  endif
  if (sg > 1.d0 .or. sg<0.d0) then
    option%io_buffer = 'TH with ice mode: Gas Saturation error:  >1 or <0'
    call printErrMsg(option)
  endif
  if (sl + si + sg /= 1.d0) then
    option%io_buffer = 'TH with ice mode: Saturation not summed to 1 '
    call printErrMsg(option)
  endif

  ! (3) relative permissivity of multiple-phase water
  kr     = 0.d0  !all init to zero
  dkr_dp = 0.d0
  dkr_dt = 0.d0
  krg    = 0.d0  !not yet done (TODO)
  dkrg_dp= 0.d0
  dkrg_dt= 0.d0

  select type (rpf_liq => characteristic_curves%liq_rel_perm_function)
    class is(rpf_Mualem_VG_liq_type)

      call rpf_liq%RelativePermeability(sl, kr, dkr_dse, option)
      dkr_dp = rpf_liq%DRelPerm_DPressure(dsl_dp, dkr_dse)

      dkr_dt = dkr_dse*dsl_dt

    class default
      option%io_buffer = 'THAuxVarComputeFreezing2: ' // &
        'TH with Ice model now only support rpf_liq: Mualem_VG_liq_type.'
      call printErrMsg(option)

  end select

  !***************  liq. phase properties **************************
  ! saturation
  global_auxvar%sat(1) = sl
  auxvar%dsat_dp = dsl_dp
  auxvar%dsat_dt = dsl_dt

  ! enthalpy (h)
  call EOSWaterEnthalpy(tc, pw, hw, dhw_dp, dhw_dt, ierr)
  hw    = hw * option%scale     ! J/kmol -> MJ/kmol
  dhw_dt = dhw_dt * option%scale
  dhw_dp = dhw_dp * option%scale
  if (.not.saturated) dhw_dp = 0.d0     !kludge since pw is constant in the unsat zone

  auxvar%h     = hw
  auxvar%dh_dp = dhw_dp
  auxvar%dh_dt = dhw_dt

  ! internal energy (u)
  u = hw - pw / denw_mol * option%scale
  du_dp = dhw_dp - &
          (dpw_dp/denw_mol - pw/(denw_mol*denw_mol)*ddenw_dp)* option%scale
  du_dt = dhw_dt + &
          pw/(denw_mol*denw_mol)*option%scale*ddenw_dt

  auxvar%u = u
  auxvar%du_dp = du_dp
  auxvar%du_dt = du_dt

  ! viscosity
  call EOSWaterSaturationPressure(tc,                     &
                                  psat, dpsat_dt, ierr)
  call EOSWaterViscosity(tc, pw, psat, dpsat_dt,          &
              visl, dvisl_dt, dvisl_dp, dvisl_dpsat, ierr)
  if (.not.saturated) dvisl_dp = 0.d0

  auxvar%vis = visl
  auxvar%kvr = kr/visl

  auxvar%dkvr_dt = dkr_dt/visl - kr/(visl*visl)*dvisl_dt
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvisl_dp

  !***************  ice phase properties **************************

  auxvar%sat_ice = si
  auxvar%dsat_ice_dp = dsi_dp
  auxvar%dsat_ice_dt = dsi_dt

  ! Calculate internal energy and derivatives for ice
  call EOSWaterInternalEnergyIce(tc, u_ice, du_ice_dt)

  auxvar%u_ice = u_ice*1.d-3                  !kJ/kmol --> MJ/kmol
  auxvar%du_ice_dt = du_ice_dt*1.d-3          !kJ/kmol/K --> MJ/kmol/K

  !***************  (air) gas phase properties **************************
  auxvar%sat_gas = sg
  auxvar%dsat_gas_dp = dsg_dp
  auxvar%dsat_gas_dt = dsg_dt

  Tk = tc + Tf_0
  p_g            = pw
  auxvar%den_gas = p_g/(IDEAL_GAS_CONSTANT*Tk)*1.d-3 !in kmol/m3
  mol_g          = psat/p_g
  C_g            = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR ! in MJ/kmol/K
  auxvar%u_gas   = C_g*Tk           ! in MJ/kmol
  auxvar%mol_gas = mol_g

  auxvar%dden_gas_dt = - p_g/(IDEAL_GAS_CONSTANT*Tk**2)*1.d-3
  dmolg_dt           = dpsat_dt/p_g
  auxvar%du_gas_dt   = C_g + (C_wv*dmolg_dt*FMWH2O - C_a*dmolg_dt*FMWAIR)*Tk
  auxvar%dmol_gas_dt = dmolg_dt

  !***************  mixture properties **************************
  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(ithrm)
  alpha_fr = th_parameter%alpha_fr(ithrm)
  Dk = th_parameter%ckwet(ithrm)
  Dk_dry = th_parameter%ckdry(ithrm)
  Dk_ice = th_parameter%ckfrozen(ithrm)

  !Soil Kersten number
  Ke = (sl+epsilon)**(alpha)
  Ke_fr = (si+epsilon)**(alpha_fr)
  auxvar%Ke = Ke
  auxvar%Ke_fr = Ke_fr

  ! Effective thermal conductivity
  auxvar%Dk_eff = Dk*Ke + Dk_ice*Ke_fr + (1.d0 - Ke - Ke_fr)*Dk_dry

  ! Derivative of Kersten number
  dKe_dp = alpha*((sl+epsilon)**(alpha - 1.d0))*dsl_dp
  dKe_dt = alpha*((sl+epsilon)**(alpha - 1.d0))*dsl_dt
  dKe_fr_dt = alpha_fr*((si+epsilon)**(alpha_fr - 1.d0))*dsi_dt
  dKe_fr_dp = alpha_fr*((si+epsilon)**(alpha_fr - 1.d0))*dsi_dp

  auxvar%dKe_dp = 0.d0
  auxvar%dKe_dt = 0.d0
  auxvar%dKe_fr_dt = 0.d0
  auxvar%dKe_fr_dp = 0.d0
  if (sl>0.d0) then
    auxvar%dKe_dp = dKe_dp
    auxvar%dKe_dt = dKe_dt
  endif
  if (si>0.d0) then
    auxvar%dKe_fr_dp = dKe_fr_dp
    auxvar%dKe_fr_dt = dKe_fr_dt
  endif

#if 0
  write(option%myrank+100, *) 'Sats:', sl, dsl_dp, dsl_dt, si, dsi_dp, dsi_dt, sg, dsg_dp, dsg_dt
  write(option%myrank+100, *) 'U:', auxvar%u, auxvar%du_dp, auxvar%du_dt, &
    auxvar%u_ice, auxvar%du_ice_dt, auxvar%u_gas, auxvar%du_gas_dt
  write(option%myrank+100, *) 'kr/kvr:', kr, dkr_dp, dkr_dt, &
    auxvar%kvr, auxvar%dkvr_dp, auxvar%dkvr_dt
  write(option%myrank+100, *) 'Ke/Ke_fr:', auxvar%Ke, auxvar%dKe_dp, auxvar%dKe_dt, &
    auxvar%Ke_fr, auxvar%dKe_fr_dp, auxvar%dKe_fr_dt
#endif

end subroutine THAuxVarComputeFreezing2

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
  do iaux = 1, aux%num_aux_bc
    call THAuxVarDestroy(aux%auxvars_bc(iaux))
  enddo  
  do iaux = 1, aux%num_aux_ss
    call THAuxVarDestroy(aux%auxvars_ss(iaux))
  enddo  
  
  if (associated(aux%auxvars)) deallocate(aux%auxvars)
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) deallocate(aux%auxvars_bc)
  nullify(aux%auxvars_bc)
  if (associated(aux%auxvars_ss)) deallocate(aux%auxvars_ss)
  nullify(aux%auxvars_ss)
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
    ! ice
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

end module TH_Aux_module
