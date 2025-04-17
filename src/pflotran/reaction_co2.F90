module Reaction_CO2_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Global_Aux_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module

  implicit none

  private

  public :: RCO2AqActCoeff, &
            RCO2MoleFraction, &
            RCO2CalculateSCO2Solubility, &
            RCO2TotalCO2

contains

! ************************************************************************** !

subroutine RCO2AqActCoeff(rt_auxvar,global_auxvar,reaction,option)
  !
  ! Computes activity coefficients of aqueous CO2
  !
  ! Author: Chuan Lu
  ! Date: 07/13/09
  !

  use Option_module
  use Option_Transport_module
  use co2eos_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscReal :: m_na, m_cl, tc, co2aqact, lngamco2, henry, pco2
  PetscReal :: sat_pressure

  tc = global_auxvar%temp
  pco2 = global_auxvar%pres(2)
  sat_pressure =0D0

  m_na = option%m_nacl; m_cl = m_na
  if (option%transport%couple_co2_salinity) then
    m_na = rt_auxvar%pri_molal(reaction%species_idx%na_ion_id)
    m_cl = rt_auxvar%pri_molal(reaction%species_idx%cl_ion_id)
  endif

  call Henry_duan_sun(tc,pco2*1D-5,henry,lngamco2, &
                      m_na,m_cl,co2aqact)

  if (reaction%species_idx%co2_aq_id > 0) then
    rt_auxvar%pri_act_coef(reaction%species_idx%co2_aq_id) = co2aqact
  elseif (reaction%species_idx%co2_aq_id < 0) then
    rt_auxvar%sec_act_coef(-reaction%species_idx%co2_aq_id) = co2aqact
  else
    co2aqact = 1.d0
  endif

end subroutine RCO2AqActCoeff

! ************************************************************************** !

function RCO2MoleFraction(rt_auxvar,global_auxvar,reaction,option)
  !
  ! Sums the total moles of primary and secondary aqueous species
  !
  ! Author: Glenn Hammond
  ! Date: 12/01/14
  !

  use Option_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  PetscReal :: RCO2MoleFraction

  PetscInt :: ico2
  PetscReal :: sum_co2, sum_mol

  if (associated(reaction%species_idx)) then
    ico2 = reaction%species_idx%co2_aq_id
  else
    option%io_buffer = 'reaction%species_idx not set in RCO2MoleFraction(). &
      &Have you defined CO2(aq) as a species in the input deck?'
    call PrintErrMsg(option)
  endif
  if (ico2 == 0) then
    option%io_buffer = 'CO2 is not set in RCO2MoleFraction(). Have you &
      &defined CO2(aq) as a species in the input deck?'
    call PrintErrMsg(option)
  endif

  if (ico2 > 0) then
    sum_co2 = rt_auxvar%pri_molal(ico2)
  else
    sum_co2 = rt_auxvar%sec_molal(-ico2)
  endif
  sum_mol = RTSumMoles(rt_auxvar,reaction,option)
  ! sum_co2 and sum_mol are both in units mol/kg water
  ! FMWH2O is in units g/mol
  ! therefore, scale by 1.d-3 to convert from mol/kg water - g water/mol water
  ! to mol/mol water -- kg water / 1000g water
  RCO2MoleFraction = sum_co2 * FMWH2O * 1.d-3 / &
                     (1.d0 + FMWH2O * sum_mol * 1.d-3)

end function RCO2MoleFraction

! ************************************************************************** !

subroutine RCO2CalculateSCO2Solubility(rt_auxvar,global_auxvar,reaction, &
                                       gas_solubility,option)

  ! RCalculateSCO2Solubility: Calculates the solubility of supercritical CO2
  ! in water
  !
  ! Author: Glenn Hammond
  ! Date: 11/12/24

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  use co2eos_module, only: Henry_duan_sun
  use co2_span_wagner_module, only: co2_span_wagner, co2_sw_itable

  implicit none

  class(reaction_rt_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: gas_solubility
  type(option_type) :: option

  PetscReal :: temperature_C
  PetscReal :: temperature_K
  PetscReal :: gas_pressure_Pa
  PetscReal :: gas_pressure_MPa
  PetscReal :: gas_pressure_bar
  PetscReal :: Henrys_constant
  PetscReal :: fugacity_MPa
  PetscReal :: na_molality, cl_molality
  PetscReal :: lngamco2
#if 0
  PetscReal :: fugacity_coefficient
  PetscReal :: saturation_pressure_Pa
  PetscReal :: saturation_pressure_MPa
  PetscReal :: co2_molality
  PetscReal :: co2_partial_pressure_bar
  PetscReal :: effective_Henrys_constant
  PetscReal :: mole_fraction_co2
  PetscErrorCode :: ierr
#endif
  PetscReal :: dummy
  PetscInt :: iflag

  temperature_C = global_auxvar%temp
  temperature_K = temperature_C + T273K
  gas_pressure_Pa = global_auxvar%pres(2)
  gas_pressure_MPa = gas_pressure_Pa*1.d-6
  gas_pressure_bar = gas_pressure_Pa*1.d-5
  if (option%transport%couple_co2_salinity) then
    na_molality = rt_auxvar%pri_molal(reaction%species_idx%na_ion_id)
    cl_molality = rt_auxvar%pri_molal(reaction%species_idx%cl_ion_id)
  else
    na_molality = option%m_nacl
    cl_molality = option%m_nacl
  endif

  call Henry_duan_sun(temperature_C,gas_pressure_bar,Henrys_constant, &
                      lngamco2,na_molality,cl_molality)
  iflag = 1
  call co2_span_wagner(gas_pressure_MPa,temperature_K,dummy,dummy,dummy, &
                       fugacity_MPa,dummy,dummy,dummy,dummy,dummy,dummy, &
                       dummy,dummy,dummy,iflag,co2_sw_itable)

#if 0
  call EOSWaterSaturationPressure(temperature_C,saturation_pressure_Pa,ierr)

  ! yco2 in Duan and Sun, 2003
  mole_fraction_co2 = 1.d0-saturation_pressure_Pa/gas_pressure_Pa
  !
  fugacity_coefficient = fugacity_MPa/gas_pressure_MPa/mole_fraction_co2
  effective_Henrys_constant = Henrys_constant*fugacity_coefficient

!     sat_pressure = sat_pressure * 1.D5
  co2_partial_pressure_bar = (gas_pressure_Pa - saturation_pressure_Pa)*1.d-5
  co2_molality = co2_partial_pressure_bar * effective_Henrys_constant


  co2_molality = co2_partial_pressure_bar * effective_Henrys_constant
  co2_molality = (gas_pressure_Pa - saturation_pressure_Pa)*1.d-5 * &
                 Henrys_constant*fugacity_coefficient
  co2_molality = (gas_pressure_Pa - saturation_pressure_Pa)*1.d-5 * &
                 Henrys_constant*fugacity_MPa / &
                 (gas_pressure_MPa*(1.d0-saturation_pressure_Pa/ &
                                         gas_pressure_Pa ))
  co2_molality = (gas_pressure_Pa - saturation_pressure_Pa)*1.d-5 * &
                 Henrys_constant*fugacity_MPa / &
                 (gas_pressure_MPa-saturation_pressure_MPa)
  co2_molality = (gas_pressure_Pa - saturation_pressure_Pa)*1.d-5 * &
                 Henrys_constant*fugacity_MPa / 7
                 ((gas_pressure_Pa-saturation_pressure_Pa)*1.d-6)
  co2_molality = 1.d-5 * Henrys_constant*fugacity_MPa/(1.d-6)
  co2_molality = Henrys_constant*fugacity_MPa*10
#endif

  ! based on comments in co2eos.F90, Henry's constant has lngamco2 built in.
  gas_solubility = fugacity_MPa*10.d0*Henrys_constant

end subroutine RCO2CalculateSCO2Solubility

! ************************************************************************** !

subroutine RCO2TotalCO2(rt_auxvar,global_auxvar,reaction,option)
  !
  ! Computes the total component concentrations and derivative with
  ! respect to free-ion for CO2 modes; this is legacy cod3
  !
  ! Author: Glenn Hammond, but originally by Chuan Lu
  ! Date: 08/01/16
  !

  use Option_module
  use EOS_Water_module
  use co2eos_module, only: Henry_duan_sun

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscErrorCode :: ierr
  PetscInt :: iphase
  PetscInt :: icomp
  PetscReal :: tempreal
  PetscReal :: pco2,sat_pressure,lngamco2
  PetscReal :: den_kg_per_L
  PetscReal :: den
  PetscReal :: lnQK
  PetscReal :: m_cl, m_na, muco2, xmass
  PetscReal :: pressure, temperature, xphico2
  PetscInt :: iactgas

! *********** Add SC phase and gas contributions ***********************
  ! CO2-specific
  iphase = 2

  rt_auxvar%total(:,iphase) = 0.D0
  rt_auxvar%aqueous%dtotal(:,:,iphase) = 0.D0

  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  den_kg_per_L = global_auxvar%den_kg(iphase)*xmass*1.d-3

  if (global_auxvar%sat(iphase) > 1.D-20) then
    do iactgas = 1, reaction%gas%nactive_gas ! all gas phase species are secondary

      pressure = global_auxvar%pres(2)
      temperature = global_auxvar%temp
      xphico2 = global_auxvar%fugacoeff(1)
      den = global_auxvar%den(2)

      call EOSWaterSaturationPressure(temperature, sat_pressure, ierr)
      pco2 = pressure - sat_pressure
!     call co2_span_wagner(pressure*1.D-6,temperature+T273K,dg,dddt,dddp,fg, &
!              dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
!
!            fg = fg*1D6
!            xphico2 = fg / pco2
!            global_auxvar%fugacoeff(1) = xphico2


      if (abs(reaction%species_idx%co2_gas_id) == iactgas ) then

        if (option%transport%couple_co2_salinity) then
          m_na = rt_auxvar%pri_molal(reaction%species_idx%na_ion_id)
          m_cl = rt_auxvar%pri_molal(reaction%species_idx%cl_ion_id)
          call Henry_duan_sun(temperature,pressure*1D-5,muco2, &
                lngamco2,m_na,m_cl)
        else
          call Henry_duan_sun(temperature,pressure*1D-5,muco2, &
                lngamco2,option%m_nacl,option%m_nacl)
        endif
        !lnQk = - log(muco2)
        lnQk = - log(muco2)-lngamco2

      else
        lngamco2 = 0.d0
        lnQK = -reaction%gas%acteqlogK(iactgas)*LOG_TO_LN
      endif

      if (reaction%gas%acteqh2oid(iactgas) > 0) then
        lnQK = lnQK + reaction%gas%acteqh2ostoich(iactgas)*rt_auxvar%ln_act_h2o
      endif

   ! contribute to %total
   !     do i = 1, ncomp
   ! removed loop over species, suppose only one primary species is related
      icomp = reaction%gas%acteqspecid(1,iactgas)
      pressure = pressure * 1.D-5

!     rt_auxvar%gas_pp(iactgas) = &
!         exp(lnQK+lngamco2)*rt_auxvar%pri_molal(icomp) &
!         /(IDEAL_GAS_CONSTANT*1.d-2*(temperature+T273K)*xphico2)

!     This form includes factor Z in pV = ZRT for nonideal gas
      rt_auxvar%gas_pp(iactgas) = &
          exp(lnQK)*rt_auxvar%pri_act_coef(icomp)*rt_auxvar%pri_molal(icomp)* &
          den/pressure/xphico2

      rt_auxvar%total(icomp,iphase) = rt_auxvar%total(icomp,iphase) + &
          reaction%gas%acteqstoich(1,iactgas)* &
          rt_auxvar%gas_pp(iactgas)

!       print *,'RTotal: ',icomp,iactgas,pressure, temperature, xphico2, &
!         global_auxvar%sat(iphase),rt_auxvar%gas_pp(iactgas), &
!         rt_auxvar%pri_act_coef(icomp)*exp(lnQK)*rt_auxvar%pri_molal(icomp) &
!         /pressure/xphico2*den


   ! contribute to %dtotal
   !      tempreal = exp(lnQK+lngamco2)/pressure/xphico2*den
!     tempreal = rt_auxvar%pri_act_coef(icomp)*exp(lnQK) &
!         /pressure/xphico2*den
      tempreal = rt_auxvar%gas_pp(iactgas)/rt_auxvar%pri_molal(icomp)
      rt_auxvar%aqueous%dtotal(icomp,icomp,iphase) = &
          rt_auxvar%aqueous%dtotal(icomp,icomp,iphase) + &
          reaction%gas%acteqstoich(1,iactgas)*tempreal
    enddo
  ! rt_auxvar%total(:,iphase) = rt_auxvar%total(:,iphase)!*den_kg_per_L
  ! units of dtotal = kg water/L water
  ! rt_auxvar%dtotal(:, :,iphase) = rt_auxvar%dtotal(:,:,iphase)!*den_kg_per_L
  endif

end subroutine RCO2TotalCO2

end module Reaction_CO2_module
