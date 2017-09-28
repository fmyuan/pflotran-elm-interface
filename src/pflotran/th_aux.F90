module TH_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none
  
  private 

  PetscInt, public :: TH_ni_count
  PetscInt, public :: TH_ts_cut_count
  PetscInt, public :: TH_ts_count
  PetscInt, public :: th_ice_model

  type, public :: TH_auxvar_type
    PetscReal :: avgmw
    PetscReal :: h
    PetscReal :: u
    PetscReal :: pc
    PetscReal :: vis
    PetscReal :: kvr
    PetscReal :: dsat_dp
    PetscReal :: dsat_dT
    PetscReal :: dden_dp
    PetscReal :: dden_dT
    PetscReal :: dkvr_dp
    PetscReal :: dkvr_dT
    PetscReal :: dh_dp
    PetscReal :: dh_dT
    PetscReal :: du_dp
    PetscReal :: du_dT
    PetscReal :: Dk_eff
    PetscReal :: Ke
    PetscReal :: dKe_dp
    PetscReal :: dKe_dT
    PetscReal :: dpres_dtime ! for TS
    PetscReal :: dtemp_dtime  ! for TS
    PetscReal :: d2sat_dp2
    PetscReal :: d2den_dp2
    PetscReal :: d2u_dp2
    PetscReal :: d2sat_dT2
    PetscReal :: d2den_dT2
    PetscReal :: d2u_dT2
    PetscReal :: d2sat_dTdp
    PetscReal :: d2den_dTdp
    PetscReal :: d2u_dTdp
    
    PetscReal :: dDk_eff_dp
    PetscReal :: dDK_eff_dt
    ! for ice
    type(th_ice_type), pointer :: ice
  end type TH_auxvar_type

  type, public :: th_ice_type
    PetscReal :: Ke_fr
    PetscReal :: dKe_fr_dp
    PetscReal :: dKe_fr_dT
    ! ice
    PetscReal :: sat_ice
    PetscReal :: sat_gas
    PetscReal :: dsat_ice_dp
    PetscReal :: dsat_gas_dp
    PetscReal :: dsat_ice_dT
    PetscReal :: dsat_gas_dT
    PetscReal :: den_ice
    PetscReal :: dden_ice_dp
    PetscReal :: dden_ice_dT
    PetscReal :: u_ice
    PetscReal :: du_ice_dp
    PetscReal :: du_ice_dT
    PetscReal :: den_gas
    PetscReal :: dden_gas_dp
    PetscReal :: dden_gas_dT
    PetscReal :: u_gas
    PetscReal :: du_gas_dp
    PetscReal :: du_gas_dT
    PetscReal :: mol_gas
    PetscReal :: dmol_gas_dp
    PetscReal :: dmol_gas_dT
    ! For DallAmico model
    PetscReal :: pres_fh2o
    PetscReal :: dpres_fh2o_dp
    PetscReal :: dpres_fh2o_dT
  end type th_ice_type
  
  type, public :: th_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckdry(:) ! Thermal conductivity (dry)
    PetscReal, pointer :: ckwet(:) ! Thermal conductivity (wet)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: ckfrozen(:) ! Thermal conductivity (frozen soil)
    PetscReal, pointer :: alpha_fr(:) ! exponent frozen
    PetscReal, pointer :: sir(:,:)
    PetscReal, pointer :: diffusion_coefficient(:)
    PetscReal, pointer :: diffusion_activation_energy(:)
  end type th_parameter_type
  
  type, public :: TH_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(th_parameter_type), pointer :: th_parameter
    type(TH_auxvar_type), pointer :: auxvars(:)
    type(TH_auxvar_type), pointer :: auxvars_bc(:)
    type(TH_auxvar_type), pointer :: auxvars_ss(:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type TH_type

  PetscReal, parameter :: epsilon = 1.d-6
  PetscReal, parameter :: perturbation_tolerance = 1.d-8

  public :: THAuxCreate, THAuxDestroy, &
            THAuxVarComputeNoFreezing, THAuxVarInit, &
            THAuxVarCopy, THAuxVarDestroy, &
            THAuxVarCompute2ndOrderDeriv

  public :: THAuxVarComputeFreezing

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
  nullify(aux%matrix_zeroing)

  allocate(aux%th_parameter)
  nullify(aux%th_parameter%dencpr)
  nullify(aux%th_parameter%ckdry)
  nullify(aux%th_parameter%ckwet)
  nullify(aux%th_parameter%alpha)
  nullify(aux%th_parameter%ckfrozen)
  nullify(aux%th_parameter%alpha_fr)
  nullify(aux%th_parameter%sir)
  nullify(aux%th_parameter%diffusion_coefficient)
  nullify(aux%th_parameter%diffusion_activation_energy)
  
  allocate(aux%th_parameter%diffusion_coefficient(option%nphase))
  allocate(aux%th_parameter%diffusion_activation_energy(option%nphase))
  aux%th_parameter%diffusion_coefficient = 1.d-9
  aux%th_parameter%diffusion_activation_energy = 0.d0
 
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
  auxvar%dsat_dT   = uninit_value
  auxvar%dden_dp   = uninit_value
  auxvar%dden_dT   = uninit_value
  auxvar%dkvr_dp   = uninit_value
  auxvar%dkvr_dT   = uninit_value
  auxvar%dh_dp     = uninit_value
  auxvar%dh_dT     = uninit_value
  auxvar%du_dp     = uninit_value
  auxvar%du_dT     = uninit_value    
  auxvar%Dk_eff    = uninit_value
  auxvar%dDk_eff_dp = uninit_value
  auxvar%dDk_eff_dt = uninit_value
  auxvar%Ke        = uninit_value
  auxvar%dKe_dp    = uninit_value
  auxvar%dKe_dT    = uninit_value
 if (option%flow%th_freezing) then
    allocate(auxvar%ice)
    auxvar%ice%Ke_fr     = uninit_value
    auxvar%ice%dKe_fr_dp = uninit_value
    auxvar%ice%dKe_fr_dT = uninit_value
    ! NOTE(bja, 2013-12) always initialize ice variables to zero, even if 
    !                    not used!
    auxvar%ice%sat_ice       = uninit_value
    auxvar%ice%sat_gas       = uninit_value
    auxvar%ice%dsat_ice_dp   = uninit_value
    auxvar%ice%dsat_gas_dp   = uninit_value
    auxvar%ice%dsat_ice_dT   = uninit_value
    auxvar%ice%dsat_gas_dT   = uninit_value
    auxvar%ice%den_ice       = uninit_value
    auxvar%ice%dden_ice_dp   = uninit_value
    auxvar%ice%dden_ice_dT   = uninit_value
    auxvar%ice%u_ice         = uninit_value
    auxvar%ice%du_ice_dp     = uninit_value
    auxvar%ice%du_ice_dT     = uninit_value
    auxvar%ice%den_gas       = uninit_value
    auxvar%ice%dden_gas_dp   = uninit_value
    auxvar%ice%dden_gas_dT   = uninit_value
    auxvar%ice%u_gas         = uninit_value
    auxvar%ice%du_gas_dp     = uninit_value
    auxvar%ice%du_gas_dT     = uninit_value
    auxvar%ice%mol_gas       = uninit_value
    auxvar%ice%dmol_gas_dp   = uninit_value
    auxvar%ice%dmol_gas_dT   = uninit_value
    auxvar%ice%pres_fh2o     = uninit_value
    auxvar%ice%dpres_fh2o_dp = uninit_value
    auxvar%ice%dpres_fh2o_dT = uninit_value
  else
    nullify(auxvar%ice)
  endif
  
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
  auxvar2%dsat_dT = auxvar%dsat_dT
  auxvar2%dden_dp = auxvar%dden_dp
  auxvar2%dden_dT = auxvar%dden_dT
  auxvar2%dkvr_dp = auxvar%dkvr_dp
  auxvar2%dkvr_dT = auxvar%dkvr_dT
  auxvar2%dh_dp = auxvar%dh_dp
  auxvar2%dh_dT = auxvar%dh_dT
  auxvar2%du_dp = auxvar%du_dp
  auxvar2%du_dT = auxvar%du_dT  
  auxvar2%Dk_eff = auxvar%Dk_eff
  auxvar2%dDk_eff_dp = auxvar%dDk_eff_dp
  auxvar2%dDk_eff_dt = auxvar%dDk_eff_dt
  auxvar2%Ke = auxvar%Ke
  auxvar2%dKe_dp = auxvar%dKe_dp
  auxvar2%dKe_dT = auxvar%dKe_dT
  if (associated(auxvar%ice)) then
    auxvar2%ice%Ke_fr = auxvar%ice%Ke_fr
    auxvar2%ice%dKe_fr_dp = auxvar%ice%dKe_fr_dp
    auxvar2%ice%dKe_fr_dT = auxvar%ice%dKe_fr_dT
    auxvar2%ice%sat_ice = auxvar%ice%sat_ice 
    auxvar2%ice%sat_gas = auxvar%ice%sat_gas
    auxvar2%ice%dsat_ice_dp = auxvar%ice%dsat_ice_dp
    auxvar2%ice%dsat_gas_dp = auxvar%ice%dsat_gas_dp
    auxvar2%ice%dsat_ice_dT = auxvar%ice%dsat_ice_dT
    auxvar2%ice%dsat_gas_dT = auxvar%ice%dsat_gas_dT
    auxvar2%ice%den_ice = auxvar%ice%den_ice
    auxvar2%ice%dden_ice_dp = auxvar%ice%dden_ice_dp
    auxvar2%ice%dden_ice_dT = auxvar%ice%dden_ice_dT
    auxvar2%ice%u_ice = auxvar%ice%u_ice
    auxvar2%ice%du_ice_dp = auxvar%ice%du_ice_dp
    auxvar2%ice%du_ice_dT = auxvar%ice%du_ice_dT
    auxvar2%ice%pres_fh2o = auxvar%ice%pres_fh2o
    auxvar2%ice%dpres_fh2o_dp = auxvar%ice%dpres_fh2o_dp
    auxvar2%ice%dpres_fh2o_dT = auxvar%ice%dpres_fh2o_dT
    auxvar2%ice%den_gas = auxvar%ice%den_gas
    auxvar2%ice%dden_gas_dp = auxvar%ice%dden_gas_dp
    auxvar2%ice%dden_gas_dT = auxvar%ice%dden_gas_dT
    auxvar2%ice%u_gas = auxvar%ice%u_gas
    auxvar2%ice%du_gas_dp = auxvar%ice%du_gas_dp
    auxvar2%ice%du_gas_dT = auxvar%ice%du_gas_dT
    auxvar2%ice%mol_gas = auxvar%ice%mol_gas
    auxvar2%ice%dmol_gas_dp = auxvar%ice%dmol_gas_dp
    auxvar2%ice%dmol_gas_dT = auxvar%ice%dmol_gas_dT
  endif

end subroutine THAuxVarCopy

! ************************************************************************** !

subroutine THAuxVarComputeNoFreezing(x,auxvar,global_auxvar, &
                                     material_auxvar, &
                                     iphase,characteristic_curves, &
                                     thermal_cc, &
                                     th_parameter, icct, natural_id, &
                                     update_porosity,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: ???
  ! Date: 02/22/08
  ! 

  use Option_module
  use Global_Aux_module
 
  use EOS_Water_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module  
  use Characteristic_Curves_Thermal_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  class(cc_thermal_type) :: thermal_cc
  PetscReal :: x(option%nflowdof)
  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: iphase
  type(th_parameter_type) :: th_parameter
  PetscInt :: icct
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  PetscBool :: update_porosity

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dT, dvis_dp
  PetscReal :: dw_dp, dw_dT, hw_dp, hw_dT
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dT
  PetscReal :: Ke
  PetscReal :: alpha
  PetscReal :: Dk
  PetscReal :: Dk_dry
  PetscReal :: aux(1)
  PetscReal :: dkr_dsat1
  PetscReal :: dk_ds, dk_dT
  
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
 
  if (update_porosity) then
    call MaterialAuxVarCompute(material_auxvar,global_auxvar%pres(1))
  endif

  auxvar%pc = min(option%flow%reference_pressure - global_auxvar%pres(1), &
                  characteristic_curves%saturation_function%pcmax)

!***************  Liquid phase properties **************************
  auxvar%avgmw = FMWH2O

  pw = option%flow%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0

!  if (auxvar%pc > 0.d0) then
  if (auxvar%pc > 1.d0) then
     iphase = 3

    call characteristic_curves%saturation_function% &
        Saturation(auxvar%pc,global_auxvar%sat(1), &
                   ds_dp, option)

    if (ds_dp < 1.d-40) then
      iphase = 1
      auxvar%pc = 0.d0
      global_auxvar%sat(1) = 1.d0
      kr = 1.d0
      pw = global_auxvar%pres(1)
      dpw_dp = 1.d0
    else  
      call characteristic_curves%liq_rel_perm_function% &
             RelativePermeability(global_auxvar%sat(1),kr, &
                                  dkr_dsat1,option) 

      dkr_dp = ds_dp * dkr_dsat1
      dpw_dp = 0.d0
    endif

  else
    iphase = 1
    auxvar%pc = 0.d0
    global_auxvar%sat(1) = 1.d0  
    kr = 1.d0    
!   pw = auxvar%pres
    pw = global_auxvar%pres(1)
    dpw_dp = 1.d0
  endif  

  ! may need to compute dpsat_dT to pass to VISW
  call EOSWaterSaturationPressure(global_auxvar%temp,sat_pressure,dpsat_dT,ierr)
  call EOSWaterEnthalpy(global_auxvar%temp,pw,hw,hw_dp,hw_dT,ierr)
  if (.not.option%flow%density_depends_on_salinity) then
    call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol,dw_dp,dw_dT,ierr)
    if (ierr /= 0) then
      call PrintMsgByCell(option,natural_id, &
                       'Error in THAuxVarComputeNoFreezing->EOSWaterDensity')
    endif
    call EOSWaterViscosity(global_auxvar%temp,pw,sat_pressure,dpsat_dT,visl, &
                           dvis_dT,dvis_dp,ierr)
  else
    aux(1) = global_auxvar%m_nacl(1)
    call EOSWaterDensityExt(global_auxvar%temp,pw,aux, &
                            dw_kg,dw_mol,dw_dp,dw_dT,ierr)
    if (ierr /= 0) then
      call PrintMsgByCell(option,natural_id, &
                     'Error in THAuxVarComputeNoFreezing->EOSWaterDensityExt')
    endif
    call EOSWaterViscosityExt(global_auxvar%temp,pw,sat_pressure,dpsat_dT,aux, &
                              visl,dvis_dT,dvis_dp,ierr)
  endif
  ! J/kmol -> whatever units
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dT = hw_dT * option%scale
  
!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
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
  auxvar%dden_dT = dw_dT
  auxvar%dsat_dT = 0.d0
  auxvar%dden_dp = dw_dp
  
!geh: contribution of dvis_dpsat is now added in EOSWaterViscosity
!  auxvar%dkvr_dT = -kr/(visl*visl)*(dvis_dT+dvis_dpsat*dpsat_dT)
  auxvar%dkvr_dT = -kr/(visl*visl)*dvis_dT
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  if (iphase < 3) then !kludge since pw is constant in the unsat zone
    auxvar%dh_dp = hw_dp
    auxvar%du_dp = hw_dp - (dpw_dp/dw_mol-pw/(dw_mol*dw_mol)*dw_dp)*option%scale
  else
    auxvar%dh_dp = 0.d0
    auxvar%du_dp = 0.d0
  endif

  auxvar%dh_dT = hw_dT
  auxvar%du_dT = hw_dT + pw/(dw_mol*dw_mol)*option%scale*dw_dT
  
  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(icct)
  Dk = th_parameter%ckwet(icct)
  Dk_dry = th_parameter%ckdry(icct)

  !unfrozen soil Kersten number
  Ke = (global_auxvar%sat(1) + epsilon)**(alpha)
  auxvar%Ke = Ke

  ! Effective thermal conductivity
  call thermal_cc%thermal_conductivity_function%CalculateTCond( &
       global_auxvar%sat(1),global_auxvar%temp,auxvar%Dk_eff,dk_ds,dk_dT,option)

  ! Derivative of soil Kersten number
  auxvar%dKe_dp = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dp
  auxvar%dKe_dT = 0.d0

end subroutine THAuxVarComputeNoFreezing

! ************************************************************************** !

subroutine THAuxVarComputeFreezing(x, auxvar, global_auxvar, &
                                   material_auxvar, &
                                   iphase, &
                                   characteristic_curves,    &
                                   thermal_cc, &
                                   th_parameter, icct, natural_id, &
                                   update_porosity,option)
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
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  class(cc_thermal_type) :: thermal_cc
  PetscReal :: x(option%nflowdof)
  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(th_parameter_type) :: th_parameter
  PetscInt :: icct
  PetscInt :: iphase
  PetscInt :: natural_id
  PetscBool :: update_porosity

  PetscErrorCode :: ierr
  PetscReal :: pw, dw_kg, dw_mol, hw, sat_pressure, visl
  PetscReal :: kr, dkr_dp, dkr_dt
  PetscReal :: dvis_dt, dvis_dp
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  PetscReal :: liq_saturation, ice_saturation, gas_saturation
  PetscReal :: dsl_dp, dsl_dt
  PetscReal :: dsg_dp, dsg_dt
  PetscReal :: dsi_dp, dsi_dt
  PetscReal :: den_ice, dden_ice_dt, dden_ice_dp
  PetscReal :: u_ice, du_ice_dt
  PetscReal :: p_th

  PetscReal :: p_g, tk_g
  PetscReal :: p_sat
  PetscReal :: mol_g
  PetscReal :: C_g
  PetscReal :: dmolg_dt, dmolg_dp
  PetscReal, parameter :: C_a  = 1.86d-3   ! in MJ/kg/K at 300K
  PetscReal, parameter :: C_wv = 1.005d-3  ! in MJ/kg/K

  PetscReal :: Ke
  PetscReal :: Ke_fr
  PetscReal :: alpha
  PetscReal :: alpha_fr
  PetscReal :: Dk
  PetscReal :: Dk_dry
  PetscReal :: Dk_ice
  PetscReal :: pcmax, pcmin
  PetscReal :: dk_ds, dK_di, dk_dT

  !--------------------------------------------------------------------------------------------------
 
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%kvr = 0.d0
   
  global_auxvar%pres = x(1)  
  global_auxvar%temp = x(2)

  ! Check if the capillary pressure is less than -100MPa, which also limit the lower limit of pres(1).
  pcmax = abs(characteristic_curves%saturation_function%pcmax)  ! non-negative
  if (global_auxvar%pres(1) - option%flow%reference_pressure <= -pcmax) then
    global_auxvar%pres(1) = -pcmax + option%flow%reference_pressure
  endif

  if (update_porosity) then
    call MaterialAuxVarCompute(material_auxvar,global_auxvar%pres(1))
  endif

  auxvar%pc = option%flow%reference_pressure - global_auxvar%pres(1)  ! always non-negative


!***************  Characteristic Curves ***********************************

  ! using modules in 'characteristic_curves.F90' to calculate needed variables:
  ! saturations and derivatives for 3 phases
  call THAuxVarComputeCharacteristicCurves(global_auxvar%pres(1), global_auxvar%temp,   &
                                           characteristic_curves,                       &
                                           liq_saturation,  dsl_dp, dsl_dt,             &
                                           ice_saturation,  dsi_dp, dsi_dt,             &
                                           auxvar%ice%pres_fh2o,                        &
                                           auxvar%ice%dpres_fh2o_dp,                    &
                                           auxvar%ice%dpres_fh2o_dt,                    &
                                           gas_saturation,  dsg_dp, dsg_dt,             &
                                           kr,              dkr_dp, dkr_dt,             &
                                           option)


!***************  3-phase water properties **********************************************************

  ! ----- Liq. water ---------------------------------------------------------
  auxvar%avgmw = FMWH2O

  pw     = option%flow%reference_pressure
  pcmin  = 0.d0  ! near-saturation PC zone (hint: a small positive value may be helpful ?)
  dkr_dp = 0.d0
  if (auxvar%pc > 1.d0) then
    iphase = 3
    dpw_dp = 0.d0
  else
    iphase    = 1
    auxvar%pc = 0.d0
    pw        = global_auxvar%pres(1)
    dpw_dp    = 1.d0
  endif

  global_auxvar%sat(1) = liq_saturation
  auxvar%dsat_dp = dsl_dp
  auxvar%dsat_dt = dsl_dt

  call EOSWaterDensity(global_auxvar%temp,      &
                        pw, dw_kg, dw_mol, dw_dp, dw_dt,ierr)

  call EOSWaterEnthalpy(global_auxvar%temp,     &
                        pw, hw,hw_dp,hw_dt,ierr)

  ! J/kmol -> MJ/kmol
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dt = hw_dt * option%scale


  call EOSWaterSaturationPressure(global_auxvar%temp, sat_pressure, dpsat_dt, ierr)

  call EOSWaterViscosity(global_auxvar%temp, pw,   &
                         sat_pressure, dpsat_dt,   &
                         visl, dvis_dt,dvis_dp, ierr)

  global_auxvar%den = dw_mol
  global_auxvar%den_kg = dw_kg
  
  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_mol * option%scale
  auxvar%kvr = kr/visl
  auxvar%vis = visl
  auxvar%dden_dt = dw_dt
  auxvar%dden_dp = dw_dp
!geh: contribution of dvis_dpsat is now added in EOSWaterViscosity  
!  auxvar%dkvr_dt = -kr/(visl*visl)*(dvis_dt + dvis_dpsat*dpsat_dt) + dkr_dt/visl
  auxvar%dkvr_dt = -kr/(visl*visl)*dvis_dt + dkr_dt/visl
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  auxvar%dh_dp = hw_dp
  auxvar%du_dp = hw_dp - (dpw_dp/dw_mol - pw/(dw_mol*dw_mol)*dw_dp)* option%scale
  auxvar%dh_dt = hw_dt
  auxvar%du_dt = hw_dt + pw/(dw_mol*dw_mol)*option%scale*dw_dt

  ! ----- ice water ---------------------------------------------------------

  auxvar%ice%sat_ice     = ice_saturation
  auxvar%ice%dsat_ice_dp = dsi_dp
  auxvar%ice%dsat_ice_dt = dsi_dt

  call EOSWaterDensityIce(global_auxvar%temp, pw,                    &
                          den_ice, dden_ice_dT, dden_ice_dp, ierr)

  call EOSWaterInternalEnergyIce(global_auxvar%temp, u_ice, du_ice_dT)

  auxvar%ice%den_ice     = den_ice
  auxvar%ice%dden_ice_dt = dden_ice_dT
  auxvar%ice%dden_ice_dp = dden_ice_dP
  auxvar%ice%u_ice     = u_ice*1.d-3              !kJ/kmol --> MJ/kmol
  auxvar%ice%du_ice_dt = du_ice_dT*1.d-3          !kJ/kmol/K --> MJ/kmol/K
  auxvar%ice%du_ice_dp = 0.d0


  ! ----- Air (incl. vapor) ---------------------------------------------------------

  auxvar%ice%sat_gas = gas_saturation
  auxvar%ice%dsat_gas_dp = dsg_dp
  auxvar%ice%dsat_gas_dt = dsg_dt

  ! Calculate the values and derivatives for vapor density and internal energy

  p_g      = pw
  tk_g     = global_auxvar%temp+273.15d0

  auxvar%ice%den_gas     = p_g/(IDEAL_GAS_CONSTANT*tk_g)*1.d-3                ! in kmol/m3 for all air-mixture
  auxvar%ice%dden_gas_dt = -p_g/(IDEAL_GAS_CONSTANT*tk_g**2)*1.d-3
  auxvar%ice%dden_gas_dp = 1.d0/(IDEAL_GAS_CONSTANT*tk_g)*1.d-3

  call EOSWaterSaturationPressure(tk_g-273.15d0, p_sat, dpsat_dt, ierr)
  mol_g    = p_sat/p_g
  dmolg_dt = dpsat_dt/p_g
  dmolg_dp = -p_sat/p_g/p_g

  C_g                    = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR          ! in MJ/kmol/K
  auxvar%ice%u_gas       = C_g*tk_g                                               ! in MJ/kmol

#if 0
! 2017-01-30: tests show that the following may cause tiny-time recovering issue when re-freezing (so OFF now).
  ! NOTE: vapor 'mol_gas' is included in 'den_gas', 'u_gas'
  ! (because 'mol_g', fraction of vapor in air-mixture, going to be as multiplier in all 'gas' calculations in 'th.F90')
  ! so, if for air-mixture, may assign this to 1.0 ( appears very helpful to reduce time-step)
  !     if totally ignore gas (all), may assign this to 0.0
#ifdef NO_VAPOR_DIFFUSION
  auxvar%ice%mol_gas     = 0.d0   ! no gas (inc. vapor)
#else
  auxvar%ice%mol_gas     = 1.d0   ! air-mixture
#endif
  auxvar%ice%dmol_gas_dt = 0.d0
  auxvar%ice%dmol_gas_dp = 0.d0
#else
  ! vapor-only (because 'mol_g', fraction of vapor in air-mixture, going to be as multiplier in all 'gas' calculations in 'th.F90'
  auxvar%ice%mol_gas     = mol_g
  auxvar%ice%dmol_gas_dt = dmolg_dt
  auxvar%ice%dmol_gas_dp = dmolg_dp
#endif
  auxvar%ice%du_gas_dt   = C_g + (C_wv*dmolg_dt*FMWH2O - C_a*dmolg_dt*FMWAIR)*tk_g
  auxvar%ice%du_gas_dp   = (C_wv*FMWH2O-C_a*FMWAIR)*dmolg_dp*tk_g


  ! ----- Thermal conductivity ---------------------------------------------------------

  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(icct)
  alpha_fr = th_parameter%alpha_fr(icct)
  Dk = th_parameter%ckwet(icct)
  Dk_dry = th_parameter%ckdry(icct)
  Dk_ice = th_parameter%ckfrozen(icct)

  !Soil Kersten number
  Ke = (global_auxvar%sat(1) + epsilon)**(alpha)
  Ke_fr = (auxvar%ice%sat_ice + epsilon)**(alpha_fr)
  auxvar%Ke = Ke
  auxvar%ice%Ke_fr = Ke_fr

  ! Effective thermal conductivity
  call thermal_cc%thermal_conductivity_function%CalculateFTCond( &
       global_auxvar%sat(1),auxvar%ice%sat_ice,global_auxvar%temp, &
       auxvar%Dk_eff,dk_ds,dK_di,dk_dT,option)

  ! Derivative of Kersten number
  auxvar%dKe_dp = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dp
  auxvar%dKe_dt = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dt
  auxvar%ice%dKe_fr_dt = alpha_fr* &
                         (auxvar%ice%sat_ice + epsilon)**(alpha_fr - 1.d0)* &
                         auxvar%ice%dsat_ice_dt
  auxvar%ice%dKe_fr_dp = alpha_fr* &
                         (auxvar%ice%sat_ice + epsilon)**(alpha_fr - 1.d0)* &
                         auxvar%ice%dsat_ice_dp

  ! derivative of 'Dk_eff': Dk_eff = Dk_dry + (Dk-Dk_dry)*Ke + (Dk_ice-Dk_dry)*Ke_fr
  auxvar%dDk_eff_dp = (Dk-Dk_dry)*auxvar%dKe_dp + (Dk_ice-Dk_dry)*auxvar%ice%dKe_fr_dp
  auxvar%dDk_eff_dt = (Dk-Dk_dry)*auxvar%dKe_dt + (Dk_ice-Dk_dry)*auxvar%ice%dKe_fr_dt

  if (option%th_ice_model == DALL_AMICO) then
    auxvar%ice%den_ice = dw_mol
    auxvar%ice%dden_ice_dt = auxvar%dden_dt
    auxvar%ice%dden_ice_dp = auxvar%dden_dp
!    auxvar%ice%u_ice = auxvar%u  ! commented out by S.Karra 06/02/14. setting
!    internal energy of ice and water might not be correct.
!    auxvar%ice%du_ice_dt = auxvar%du_dt

    auxvar%ice%sat_gas       = 0.d0
    auxvar%ice%dsat_gas_dp   = 0.d0
    auxvar%ice%dsat_gas_dt   = 0.d0
    auxvar%ice%den_gas       = 0.d0
    auxvar%ice%dden_gas_dt   = 0.d0
    auxvar%ice%u_gas         = 0.d0
    auxvar%ice%du_gas_dt     = 0.d0
    auxvar%ice%mol_gas       = 0.d0
    auxvar%ice%dmol_gas_dt   = 0.d0
  endif

end subroutine THAuxVarComputeFreezing

! ************************************************************************** !
subroutine THAuxVarCompute2ndOrderDeriv(TH_auxvar,global_auxvar, &
                                        material_auxvar,th_parameter, &
                                        icct,characteristic_curves,&
                                        thermal_cc,&
                                        option)

  ! Computes 2nd order derivatives auxiliary variables for each grid cell
  ! 
  ! Author: Satish Karra
  ! Date: 06/06/2019
  ! 
  
  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  class(cc_thermal_type) :: thermal_cc
  type(TH_auxvar_type) :: TH_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar  
  PetscInt :: icct
  PetscErrorCode :: ierr
  
  type(th_parameter_type) :: th_parameter
  PetscInt :: iphase, ideriv
  type(TH_auxvar_type) :: TH_auxvar_pert
  type(global_auxvar_type) :: global_auxvar_pert
  ! leave as type
  type(material_auxvar_type) :: material_auxvar_pert
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), pert

  call GlobalAuxVarInit(global_auxvar_pert,option)
  call MaterialAuxVarInit(material_auxvar_pert,option)

  call THAuxVarCopy(TH_auxvar,TH_auxvar_pert,option)
  call GlobalAuxVarCopy(global_auxvar,global_auxvar_pert,option)
  call MaterialAuxVarCopy(material_auxvar,material_auxvar_pert,option)

  x(1) = global_auxvar%pres(1)
  x(2) = global_auxvar%temp

!  if (option%flow%th_freezing) then
!    option%io_buffer = 'ERROR: TH_TS MODE not implemented with freezing'
!    call PrintErrMsg(option)
!  else
!    call THAuxVarComputeNoFreezing(x,TH_auxvar,&
!                      global_auxvar,material_auxvar,&
!                      iphase,sat_func, &
!                      th_parameter,icct, &
!                      -999,option)
!  endif
  
  TH_auxvar%d2sat_dp2 = 0.d0 
  TH_auxvar%d2den_dp2 = 0.d0 
  TH_auxvar%d2u_dp2 = 0.d0 
  TH_auxvar%d2sat_dT2 = 0.d0 
  TH_auxvar%d2den_dT2 = 0.d0 
  TH_auxvar%d2u_dT2 = 0.d0 


  do ideriv = 1,option%nflowdof
    pert = x(ideriv)*perturbation_tolerance
    x_pert = x
    if (option%flow%th_freezing) then
       if (ideriv == 1) then
          if (x_pert(ideriv) < option%flow%reference_pressure) then
             pert = - pert
          endif
          x_pert(ideriv) = x_pert(ideriv) + pert
       endif

       if (ideriv == 2) then
          if (x_pert(ideriv) < 0.d0) then
             pert = - 1.d-8
          else
             pert =  1.d-8
          endif
          x_pert(ideriv) = x_pert(ideriv) + pert
       endif
    else
       x_pert(ideriv) = x_pert(ideriv) + pert
    endif

    if (option%flow%th_freezing) then
      option%io_buffer = 'ERROR: TH_TS MODE not implemented with freezing'
      call PrintErrMsg(option)
    else
      call THAuxVarComputeNoFreezing(x_pert,TH_auxvar_pert,&
                            global_auxvar_pert,material_auxvar_pert,&
                            iphase,characteristic_curves, &
                            thermal_cc, &
                            th_parameter,icct, &
                            -999,PETSC_TRUE,option)
    endif
    
    
    if (ideriv == 1) then
      TH_auxvar%d2sat_dp2 = (TH_auxvar_pert%dsat_dp - TH_auxvar%dsat_dp)/pert
      TH_auxvar%d2sat_dTdp = (TH_auxvar_pert%dsat_dT - TH_auxvar%dsat_dT)/pert
      TH_auxvar%d2den_dp2 = (TH_auxvar_pert%dden_dp - TH_auxvar%dden_dp)/pert
      TH_auxvar%d2den_dTdp = (TH_auxvar_pert%dden_dT - TH_auxvar%dden_dT)/pert
      TH_auxvar%d2u_dp2 = (TH_auxvar_pert%du_dp - TH_auxvar%du_dp)/pert
      TH_auxvar%d2u_dTdp = (TH_auxvar_pert%du_dT - TH_auxvar%du_dT)/pert
    endif  
    
    if (ideriv == 2) then
      TH_auxvar%d2sat_dT2 = (TH_auxvar_pert%dsat_dT - TH_auxvar%dsat_dT)/pert
      TH_auxvar%d2den_dT2 = (TH_auxvar_pert%dden_dT - TH_auxvar%dden_dT)/pert
      TH_auxvar%d2u_dT2 = (TH_auxvar_pert%du_dT - TH_auxvar%du_dT)/pert
    endif   
    
  enddo

end subroutine THAuxVarCompute2ndOrderDeriv  

! ************************************************************************** !

subroutine THAuxVarComputeCharacteristicCurves( pres_l,  tc,                &
                                    characteristic_curves,                  &
                                    sl,  dsl_dpl, dsl_dt,                   &
                                    si,  dsi_dpl, dsi_dt,                   &
                                    ice_presl, ice_presl_dpl, ice_presl_dt, &
                                    sg,  dsg_dpl, dsg_dt,                   &
                                    kr,  dkr_dpl, dkr_dt,                   &
                                    option)
  !
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present
  !
  ! Revised by fengming Yuan @03-08-2016/CCSI-ONRL
  ! NOTE: (1) ice_model 'PAINTER_EXPLICIT', or, 'PAINTER_KARRA_EXPLICIT',
  !             or, 'PAINTER_KARRA_EXPLICIT_SMOOTH', or, 'DALL_AMICO'
  !       (2) ANY saturation_function, from 'Characteristic_Curves_module'
  !       (3) ANY permissivity function, from 'Characteristci_Curves_module' as well

  use Option_module
  use Characteristic_Curves_module
  implicit none

  type(option_type) :: option
  PetscReal, intent(in) :: pres_l   ! unit: Pa (liq water-air interface pressure: -pc+reference_pressure)
  PetscReal, intent(in) :: tc       ! unit: oC
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal, intent(out) :: sl,  dsl_dpl, dsl_dt
  PetscReal, intent(out) :: si,  dsi_dpl, dsi_dt
  PetscReal, intent(out) :: ice_presl, ice_presl_dpl, ice_presl_dt
  PetscReal, intent(out) :: sg,  dsg_dpl, dsg_dt
  PetscReal, intent(out) :: kr,  dkr_dpl, dkr_dt

  ! local variables

  PetscReal :: pc
  PetscReal :: sli, dsli_dpl, xplice, dxplice_dpl, dxplice_dt, slx, dslx_dx
  PetscReal :: dkr_dsl

  PetscReal :: se, dse_dpc, function_A, dfunc_A_dt, function_B, dfunc_B_dpl

  ! ----------------
  ! (0) inputs
  pc = max(0.d0, option%flow%reference_pressure - pres_l)   ! always non-negative (0 = saturated)
  if (pc > abs(characteristic_curves%saturation_function%pcmax)) then
    pc = characteristic_curves%saturation_function%pcmax
  endif

  !
  ! (1) saturations
  sl      = 0.d0  !all init to zero
  dsl_dpl = 0.d0
  dsl_dt  = 0.d0
  si      = 0.d0
  dsi_dpl = 0.d0
  dsi_dt  = 0.d0
  sg      = 1.d0
  dsg_dpl = 0.d0
  dsg_dt  = 0.d0

  ice_presl    = pres_l
  ice_presl_dpl= 1.d0
  ice_presl_dt = 0.d0

  ! initial liq. saturation and derivatives
  call characteristic_curves%saturation_function%Saturation(pc, sli, dsli_dpl, option)  ! a note: in CC modules, already w.r.t 'pres_l' by dpc_dpres = -1.d0
  sl      = sli
  dsl_dpl = dsli_dpl
  dsl_dt  = 0.d0

  ! if ice module turns on, 2-phase saturation recalculated (liq. and ice) under common 'pres_l' and 'tc'
  xplice = pc   ! liq. water PC at ice-liq interface (positive, if no ice it's pc by default)
  if (option%flow%th_freezing) then

    call characteristic_curves%saturation_function%IceCapillaryPressure(pres_l, tc, &
                                   xplice, dxplice_dpl, dxplice_dt, option) ! w.r.t 'pressure' already in %IceCapillaryPressure()

     select case (option%th_ice_model)
       case (PAINTER_EXPLICIT)

        ! NOTE: in 'saturation_function.F90': SatFuncComputeIcePExplicit(), 'se' and 'dse_dpc' are used
        ! while in 'characteristic_curves.F90': SF_VG_Saturaton(), the following are output:
        !   sl = this%Sr + (1.d0-this%Sr)*Se
        !   dsl_dpl = -(1.d0-this%Sr)*dSe_dpc

        function_B = 1.d0
        dfunc_B_dpl= 0.d0
        if (pc>0.d0) then
          se = (sl - characteristic_curves%saturation_function%Sr) &
               /(1.0d0 - characteristic_curves%saturation_function%Sr)
          dse_dpc = (-dsl_dpl)/(1.0d0 - characteristic_curves%saturation_function%Sr)

          function_B = 1.d0/se
          dfunc_B_dpl= 1.d0/(se*se)*dse_dpc
        endif
        !
        function_A = 1.d0
        dfunc_A_dt = 0.d0
        if(tc<0.d0) then
          call characteristic_curves%saturation_function%Saturation(xplice, slx, dslx_dx, option)
          se = (slx - characteristic_curves%saturation_function%Sr) &
              /(1.0d0 - characteristic_curves%saturation_function%Sr)
          dse_dpc = (-dslx_dx)/(1.0d0 - characteristic_curves%saturation_function%Sr)

          function_A = 1.d0/se
          ! dfunc_A_dt = dfuncA_dslx * dslx_dx * dxplice_dt   (note: dslx_dx = dslx_dxplice)
          dfunc_A_dt = 1.d0/(se*se)* dse_dpc
          dfunc_A_dt = dfunc_A_dt * (-dxplice_dt)

        endif

         sl = 1.d0/(function_A + function_B - 1.d0)
        sg = sl*(function_B - 1.d0)
        si = sl*(function_A - 1.d0)

        dsl_dpl = - 1.d0/(function_A + function_B - 1.d0)**(2.d0)*(dfunc_B_dpl)
        dsl_dt  = - 1.d0/(function_A + function_B - 1.d0)**(2.d0)*(dfunc_A_dt)

        dsg_dpl = dsl_dpl*(function_B - 1.d0) + sl*dfunc_B_dpl
        dsg_dt  = dsl_dt*(function_B - 1.d0)

        dsi_dpl = dsl_dpl*(function_A - 1.d0)
        dsi_dt  = dsl_dt*(function_A - 1.d0) + sl*dfunc_A_dt

      case (PAINTER_KARRA_EXPLICIT, PAINTER_KARRA_EXPLICIT_SMOOTH)

#if 0
        ! (TODO-checking) when coupled with CLM, the following not works well
        ! although it's theoretically better (because it's used for calculating 'dphi' for water Darcy flux)
        ice_presl    = (pres_l - pc) - xplice
        ice_presl_dpl= dxplice_dpl
        ice_presl_dt = dxplice_dt
#else
        ice_presl    = pres_l   ! this implies that actually not use ice-adjusted pressure
        ice_presl_dpl= 1.d0
        ice_presl_dt = 0.d0
#endif

        call characteristic_curves%saturation_function%Saturation(xplice, slx, dslx_dx, option)

        ! liq. saturation and its derivatives, with ice-adjusted capillary pressure
        sl      = slx
        dsl_dt  = -dslx_dx*dxplice_dt
        dsl_dpl = dslx_dx*dxplice_dpl

        ! ice satuation and its derivatives
        si     = 1.d0 - sl/sli             ! P.-K. Eq.(19)
        dsi_dt = -1.d0/sli*dsl_dt          ! dsli_dt = 0 (see above)
        dsi_dpl= (sl*dsli_dpl-sli*dsl_dpl)/(sli**2)

        ! gas component (as difference)
        sg      = 1.d0 - sl - si
        dsg_dpl = -dsl_dpl - dsi_dpl
        dsg_dt  = -dsl_dt - dsi_dt

      case (DALL_AMICO)

        ! Model from Dall'Amico (2010) and Dall' Amico et al. (2011)
        ! rewritten following 'saturation_function.F90:SatFuncComputeIceDallAmico()'
        ! NOTE: here calculate 'saturations and its derivatives'
#if 0
        ice_presl    = pres_l   ! this implies that actually not use ice-adjusted pressure
        ice_presl_dpl= 1.d0
        ice_presl_dt = 0.d0
#endif
        call characteristic_curves%saturation_function%Saturation(xplice, slx, dslx_dx, option)   ! Pc1 ---> S1, but '_dx' is w.r.t. '_dpres'

        !
        ! liq. saturation and its derivatives, with ice-adjusted capillary pressure
        sl     = slx
        dsl_dpl= dslx_dx * dxplice_dpl
        dsl_dt = -dslx_dx * dxplice_dt

        ! ice satuation and its derivatives
        si     = sli - sl
        dsi_dpl= dsli_dpl - dsl_dpl
        dsi_dt = -dsl_dt                 ! dsli_dt = 0 (see above)

        ! gas phase
        sg      = 1.d0 - sl - si
        dsg_dpl = -dsl_dpl - dsi_dpl
        dsg_dt  = -dsl_dt - dsi_dt

      case default
        option%io_buffer = 'Ice module NOT recognized'
        call printErrMsg(option)

    end select
  endif ! 'option%flow%th_freezing'

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
    option%io_buffer = 'TH with ice mode: Gas Saturation error:  >1 or <0'
    call printErrMsg(option)
  endif
  if (abs((sl + si + sg)-1.d0)>1.d-10) then
    option%io_buffer = 'TH with ice mode: Saturation not summed to 1 '
    call printErrMsg(option)
  endif

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

subroutine THPrintAuxVars(file_unit,th_auxvar,global_auxvar, &
                          material_auxvar,natural_id,string,option)
  
  ! Prints the content of an TH auxvar to the designated file unit

  ! Author: Glenn Hammond
  ! Date: 07/16/20

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  PetscInt :: file_unit
  type(th_auxvar_type) :: th_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  PetscReal ::  liquid_mass

  liquid_mass = material_auxvar%volume*material_auxvar%porosity* &
                global_auxvar%sat(1)*global_auxvar%den(1)

  write(file_unit,*) '--------------------------------------------------------'
  if (len_trim(string) > 0) write(file_unit,*) trim(string)
  write(file_unit,*) '                  cell id: ', natural_id
  write(file_unit,*) '       liquid mass [kmol]: ', liquid_mass
  write(file_unit,*) '          liquid pressure: ', global_auxvar%pres(1)
  write(file_unit,*) '       capillary pressure: ', th_auxvar%pc
  write(file_unit,*) '          temperature [C]: ', global_auxvar%pres(1)
  write(file_unit,*) '          liquid enthalpy: ', th_auxvar%h
  write(file_unit,*) '   liquid internal energy: ', th_auxvar%u
  write(file_unit,*) '         liquid viscosity: ', th_auxvar%vis
  write(file_unit,*) '          liquid mobility: ', th_auxvar%kvr
  write(file_unit,*) '     liquid relative perm: ', th_auxvar%kvr * &
                                                    th_auxvar%vis
  write(file_unit,*) '    liquid_density [kmol]: ', global_auxvar%den(1)
  write(file_unit,*) '      liquid_density [kg]: ', global_auxvar%den_kg(1)
  write(file_unit,*) 'eff. thermal conductivity: ', th_auxvar%Dk_eff
  write(file_unit,*) '                porosity : ', material_auxvar%porosity
  write(file_unit,*) '             volume [m^3]: ', material_auxvar%volume
  write(file_unit,*) '--------------------------------------------------------'

end subroutine THPrintAuxVars
  
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

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%th_parameter)) then
    if (associated(aux%th_parameter%diffusion_coefficient)) &
      deallocate(aux%th_parameter%diffusion_coefficient)
    nullify(aux%th_parameter%diffusion_coefficient)
    if (associated(aux%th_parameter%diffusion_activation_energy)) &
      deallocate(aux%th_parameter%diffusion_activation_energy)
    nullify(aux%th_parameter%diffusion_activation_energy)
    if (associated(aux%th_parameter%dencpr)) deallocate(aux%th_parameter%dencpr)
    nullify(aux%th_parameter%dencpr)
    if (associated(aux%th_parameter%ckwet)) deallocate(aux%th_parameter%ckwet)
    nullify(aux%th_parameter%ckwet)
    if (associated(aux%th_parameter%ckdry)) deallocate(aux%th_parameter%ckdry)
    nullify(aux%th_parameter%ckdry)
    if (associated(aux%th_parameter%alpha)) deallocate(aux%th_parameter%alpha)
    nullify(aux%th_parameter%alpha)
    ! ice
    if (associated(aux%th_parameter%ckfrozen)) &
      deallocate(aux%th_parameter%ckfrozen)
    nullify(aux%th_parameter%ckfrozen)
    if (associated(aux%th_parameter%alpha_fr)) &
      deallocate(aux%th_parameter%alpha_fr)
    nullify(aux%th_parameter%alpha_fr)

    if (associated(aux%th_parameter%sir)) deallocate(aux%th_parameter%sir)
    nullify(aux%th_parameter%sir)
  endif
  nullify(aux%th_parameter)
  
  deallocate(aux)
  nullify(aux)  

  end subroutine THAuxDestroy

end module TH_Aux_module
