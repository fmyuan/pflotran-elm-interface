module ELM_RspFuncs_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

! constants  
  PetscReal, parameter, public :: rpi = 3.14159265358979323846d0

! temperature response function
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_OFF  = 0
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_CLMCN= 1
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_Q10  = 2
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_DLEM = 3
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_ARRHENIUS = 4

! moisture response function
  PetscInt, parameter, public :: MOISTURE_RESPONSE_FUNCTION_OFF   = 0
  PetscInt, parameter, public :: MOISTURE_RESPONSE_FUNCTION_CLMCN = 1
  PetscInt, parameter, public :: MOISTURE_RESPONSE_FUNCTION_DLEM  = 2
  PetscInt, parameter, public :: MOISTURE_RESPONSE_FUNCTION_LOGTHETA = 3

! pH response function
  PetscInt, parameter, public :: PH_RESPONSE_FUNCTION_OFF     = 0
  PetscInt, parameter, public :: PH_RESPONSE_FUNCTION_CENTURY = 1 
  PetscInt, parameter, public :: PH_RESPONSE_FUNCTION_DLEM    = 2

! molecular weight
  PetscReal, parameter, public :: N_molecular_weight = 14.0067d0
  PetscReal, parameter, public :: C_molecular_weight = 12.0110d0
  PetscReal, parameter, public :: CN_ratio_mass_to_mol = 1.166156023644992d0
  ! A NOTE here: when coupled with CLM-CN, make sure that the above constants ARE consistent with CLM
  ! otherwise may cause some tiny but detectable mass-balance errors due to unit conversion.

  ! aerobic condition function TYPE
  PetscInt, parameter, public :: Ox_RESPONSE_FUNCTION_OFF   = 0   !
  PetscInt, parameter, public :: Ox_RESPONSE_FUNCTION_MONOD = 1   ! O2 concentration based Monod function
  PetscInt, parameter, public :: Ox_RESPONSE_FUNCTION_WFPS  = 2   ! Water-Filled-Pore-Space (Parton et al. 1998, CENTURY)

  PetscInt, parameter, public :: INHIBITION_THRESHOLD     = 1
  PetscInt, parameter, public :: INHIBITION_THERMODYNAMIC = 2   ! (TODO)
  PetscInt, parameter, public :: INHIBITION_MONOD         = 3   ! this is actually an 'inhibition' equation: k/(k+s)
  PetscInt, parameter, public :: INHIBITION_INVERSE_MONOD = 4   ! this is actually a 'MONOD' equation: s/(k+s)

  public :: GetTemperatureResponse, &
            GetMoistureResponse, &
            GetpHResponse, &
            GetAerobicCondition, &
            FuncMonod, &
            FuncInhibition

contains
! ************************************************************************** !
! temperature response function 

function GetTemperatureResponse(tc, itype, Q10orEA)

  implicit none
  
  PetscInt  :: itype
  PetscReal :: tc
  PetscReal,optional :: Q10orEA   ! it's upon 'itype'
  PetscReal :: Ft, tk, Q10, EA

  PetscReal, parameter :: one_over_71_02 = 1.408054069d-2
  PetscReal :: GetTemperatureResponse
  PetscReal, parameter :: Frz_Q10 = 2.0d0

  select case(itype)
!     CLM4.5 temperature response function
      case(TEMPERATURE_RESPONSE_FUNCTION_Q10)
        Q10 = Q10orEA
        if(tc>0.d0) then
          Ft = Q10 ** ((tc - 25.0d0) / 10.0d0)
        else
          Ft = (Q10**(-25.0d0/10.0d0))*(Frz_Q10**((tc/10.0d0)))
        endif

!     CLM-CN
!     Equation: F_t = exp(308.56*(1/17.02 - 1/(T - 227.13)))
      case(TEMPERATURE_RESPONSE_FUNCTION_CLMCN)
  
        tk = tc + 273.15d0

        if(tk > 227.15d0) then
          Ft = exp(308.56d0*(one_over_71_02 - 1.d0/(tk - 227.13d0)))
        else
          Ft = 0.d0
        endif
!     DLEM temperature response function for methane oxidation
!     Tian et al. 2010 Biogeosciences, 7, 2673-2694 Eq. 12
      case(TEMPERATURE_RESPONSE_FUNCTION_DLEM)
        Q10 = Q10orEA
        if(tc < -5.0d0) then
          Ft = 0.0d0
        elseif(tc >= 30.0d0) then
          Ft = 1.0d0
        else
          Ft = Q10 ** ((tc - 30.0d0) / 10.0d0)
        endif

      ! Arrhenius equation
      case(TEMPERATURE_RESPONSE_FUNCTION_ARRHENIUS)
        EA = Q10orEA
        Ft = exp(EA/IDEAL_GAS_CONSTANT* &
                (1.d0/298.15d0-1.d0/(tc+273.15d0)))

      ! no temperature dependence
      case default
        Ft = 1.0d0
  end select

  GetTemperatureResponse = Ft 

end function GetTemperatureResponse
  
! ************************************************************************** !

Function GetMoistureResponse(theta, ghosted_id, itype)


#ifdef ELM_PFLOTRAN
  use petscvec
  use elmpf_interface_data
#endif
  
  implicit none
  
  PetscReal :: F_theta
  ! theta IS soil VWC: 0 - porosity (not adjusted)
  PetscReal :: theta
  PetscInt  :: ghosted_id, itype
  PetscReal :: GetMoistureResponse

  PetscReal, parameter :: theta_min = 0.01d0
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1

  PetscReal :: maxpsi, psi, lsat
  PetscReal, parameter :: minpsi = -10.0d6    ! Pa

#ifdef ELM_PFLOTRAN
  PetscScalar, pointer :: sucsat_pf_loc(:)    !
  PetscScalar, pointer :: watfc_pf_loc(:)     !
  PetscScalar, pointer :: porosity_pf_loc(:)  !
  PetscScalar, pointer :: bd_dry_pf_loc(:)    !
  PetscScalar, pointer :: bsw_pf_loc(:)    !
  PetscReal :: thetar, thetas, se
#endif

  PetscErrorCode :: ierr

#ifdef ELM_PFLOTRAN

  select case(itype)
!   CLM-CN
    case(MOISTURE_RESPONSE_FUNCTION_CLMCN)
      call VecGetArrayReadF90(elm_pf_idata%sucsat_pfs, sucsat_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecGetArrayReadF90(elm_pf_idata%bulkdensity_dry_pfs, bd_dry_pf_loc, ierr)   ! 'bd' (kg/m3)
      CHKERRQ(ierr)
      call VecGetArrayReadF90(elm_pf_idata%bsw_pfs, bsw_pf_loc, ierr)
      CHKERRQ(ierr)
      ! sucsat [mm of H20] from CLM is the suction (positive) at water saturated (called air-entry pressure)
      ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
      maxpsi = sucsat_pf_loc(ghosted_id) * (-EARTH_GRAVITY)                         ! mmH2O --> -Pa
      lsat = theta/min(1.d0, 1.d0-min(0.9999d0,bd_dry_pf_loc(ghosted_id)/2.70d3))     ! bd = (1._r8-dry_porosity)*2.7d3

      ! soil matric potential by Clapp-Hornburger method (this is the default used by CLM)
      psi = sucsat_pf_loc(ghosted_id) * (-EARTH_GRAVITY) * (lsat**(-bsw_pf_loc(ghosted_id)))  ! mmH2O --> -Pa
      psi = min(psi, maxpsi)
      if(psi > minpsi) then
        F_theta = log(minpsi/psi)/log(minpsi/maxpsi)
        ! very wet soil (close to saturated)
        if(psi>(maxpsi-1.d02)) then
          F_theta = F_theta*0.10d0   ! 0.10 is an arbitrary value here (but NOT totaly shut off decomp)
        endif
      else
        F_theta = 0.0d0
      endif

      call VecRestoreArrayReadF90(elm_pf_idata%sucsat_pfs, sucsat_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayReadF90(elm_pf_idata%bulkdensity_dry_pfs, bd_dry_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayReadF90(elm_pf_idata%bsw_pfs, bsw_pf_loc, ierr)
      CHKERRQ(ierr)

!     DLEM
!     Tian et al. 2010 Biogeosciences, 7, 2673-2694 Eq. 13
    case(MOISTURE_RESPONSE_FUNCTION_DLEM) 
      call VecGetArrayReadF90(elm_pf_idata%effporosity_pfs, porosity_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecGetArrayReadF90(elm_pf_idata%watfc_pfs, watfc_pf_loc, ierr)
      CHKERRQ(ierr)
      thetas = porosity_pf_loc(ghosted_id)
      thetar = watfc_pf_loc(ghosted_id)
      if(theta >= thetas) then
        F_theta = 1.0d0
      elseif (theta <= thetar) then
        F_theta = 0.0d0
      else
        se = (theta - thetar)/(thetas - thetar)
        F_theta = 1.0 - se * se * 0.368 * exp(se)

        if(F_theta < 0.0d0) then
           F_theta = 0.0d0
        endif

        if(F_theta > 1.0d0) then
           F_theta = 1.0d0
        endif
      endif
      call VecRestoreArrayReadF90(elm_pf_idata%effporosity_pfs, porosity_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayReadF90(elm_pf_idata%watfc_pfs, watfc_pf_loc, ierr)
      CHKERRQ(ierr)
    case default
        F_theta = 1.0d0
  end select
#else

  ! inhibition due to moisture content
  ! Equation: F_theta = log(theta_min/theta) / log(theta_min/theta_max)
  ! Assumptions: theta is saturation
  !              theta_min = 0.01, theta_max = 1.
  F_theta = 1.0d0 !log(theta_min/max(theta_min,theta)) * one_over_log_theta_min 
  
#endif

  GetMoistureResponse = F_theta

end function GetMoistureResponse

! ************************************************************************** !

Function GetpHResponse(pH, itype)

  PetscReal :: f_ph, pH, GetpHResponse
  PetscInt :: itype

! ph function from Parton et al., (2001, 1996)
!  k_nitr_ph_vr(c,j) = 0.56 + atan(rpi * 0.45 * (-5.+ pH(c)))/rpi
#ifdef ELM_PFLOTRAN
  select case(itype)
      case(PH_RESPONSE_FUNCTION_CENTURY)
          f_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + pH))/rpi

!     DLEM temperature response function for methane oxidation
! Tian et al. 2010 Biogeosciences, 7, 2673-2694 Eq. 12
      case(PH_RESPONSE_FUNCTION_DLEM)
          if(pH <= 4.0d0 .or. pH >= 10.0d0) then
            f_ph = 0.0d0
          elseif(pH < 7.0d0) then
            f_ph = 1.02d0 /(1.0d0 + 1.0d6 * exp(-2.5d0 * pH))
          else
            f_ph = 1.02d0 /(1.0d0 + 1.0d6 * exp(-2.5d0 * (14.0d0 - pH)))
          endif
      case default
          f_ph = 1.0d0
  end select
#else
  f_ph = 1.0
#endif

  if(f_ph < 0.0d0) then
     f_ph = 0.0d0
  endif

  if(f_ph > 1.0d0) then
     f_ph = 1.0d0
  endif
  GetpHResponse = f_ph

end function GetpHResponse

! ************************************************************************** !
Function GetAerobicCondition(OXorWFPS, K_Ox, itype, compute_derivative)

  implicit none

  PetscReal :: F_Ox    ! 0 - 1
  PetscReal :: OXorWFPS  ! upon 'itype', oxidizer concentration or equavalent (as long as its unit consistent with K_Ox)
                         ! OR, water filled pore space (WFPS), i.e. liquid water saturation
  PetscInt  :: itype
  PetscReal :: K_Ox    ! O2 or other alternative oxidier effect expressed as Ox/(K_Ox+Ox), i.e. Monod function
  PetscBool :: compute_derivative

  PetscReal :: GetAerobicCondition

  !PetscReal :: Sg_crit = 0.10d0    ! critical level of air saturation below which anaerobic getting strong
  !
  select case(itype)
    case(Ox_RESPONSE_FUNCTION_MONOD)
      F_Ox = FuncMonod(OXorWFPS, K_Ox, compute_derivative)

    case(Ox_RESPONSE_FUNCTION_WFPS)
      F_Ox = ((1.27d0 - OXorWFPS)/0.67d0)**(3.1777d0) * &
        ((OXorWFPS - 0.0012d0)/0.5988d0)**2.84d0

      if(compute_derivative) F_Ox = 0.d0  ! independent of Ox species on saturation

    case default
      !
      F_Ox = 1.d0
      if(compute_derivative) F_Ox = 0.d0

  end select

  GetAerobicCondition = F_Ox

end function GetAerobicCondition

! ************************************************************************** !
! Monod function
Function FuncMonod(conc, monod_k, compute_derivative)

  implicit none

  PetscBool :: compute_derivative
  PetscReal :: conc, monod_k
  PetscReal :: FuncMonod

  !----------------------------------------------------------
  if (.not.compute_derivative) then
    FuncMonod = conc/(conc+monod_k)
  else
    FuncMonod = monod_k/(conc+monod_k)/(conc+monod_k)
  endif

end function FuncMonod
! ************************************************************************** !
  ! Microbial-type inhibition function
  ! Followed what in 'reaction_microbial.F90', excluding 'INVERSE_MONOD' which actually is 'MONOD' as above
  function FuncInhibition(compute_derivative, conc, inhibition_C, inhibition_C2)
    implicit none

    PetscReal :: FuncInhibition
    PetscBool :: compute_derivative
    PetscReal :: conc
    PetscReal :: inhibition_C
    PetscReal, optional:: inhibition_C2
    PetscReal :: tempreal

    ! no inhibition
    FuncInhibition = 1.d0
    if(compute_derivative) FuncInhibition = 0.d0

    ! typical inhibition function
    if (.not.present(inhibition_C2)) then
      FuncInhibition = inhibition_C / (conc + inhibition_C)

      if(compute_derivative) then
        FuncInhibition = -inhibition_C                                    &
                        / (conc + inhibition_C) / (conc + inhibition_C)
      endif

     ! THRESHOLD inhibition function
     else
       FuncInhibition = 0.5d0 + &
                   atan((conc - inhibition_C) * inhibition_C2) / PI

       if(compute_derivative) then
         tempreal = (conc - inhibition_C) * inhibition_C2
         FuncInhibition = (inhibition_C2 / (1.d0 + tempreal*tempreal)) / PI
       endif

    end if

  end function FuncInhibition

! ************************************************************************** !

end module ELM_RspFuncs_module
