module CLM_BGC_module

  implicit none

  private

#include "finclude/petscsys.h"

! constants  
  PetscReal, parameter, public :: rpi = 3.14159265358979323846

! unit conversion
                          ! 14.00674d0 / 12.011d0
  PetscReal, parameter, public :: CN_ratio_mass_to_mol = 1.16616d0 

! temperature response function
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_CLM4 = 1 
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_Q10 = 2 
  PetscInt, parameter, public :: TEMPERATURE_RESPONSE_FUNCTION_DLEM = 3

! moisture response function
  PetscInt, parameter, public :: MOISTURE_RESPONSE_FUNCTION_CLM4 = 1 
  PetscInt, parameter, public :: MOISTURE_RESPONSE_FUNCTION_DLEM = 2

! pH response function
  PetscInt, parameter, public :: PH_RESPONSE_FUNCTION_CENTURY = 1 
  PetscInt, parameter, public :: PH_RESPONSE_FUNCTION_DLEM = 2

! molecular weight
  PetscReal, parameter, public :: N_molecular_weight = 14.0067d0

  public :: GetTemperatureResponse, &
            GetMoistureResponse, &
            GetpHResponse

contains
! ************************************************************************** !
! temperature response function 

function GetTemperatureResponse(tc, itype, Q10)

  implicit none
  
  PetscInt  :: itype
  PetscReal :: Ft, tc, Q10, tk

  PetscReal, parameter :: one_over_71_02 = 1.408054069d-2
  PetscReal :: GetTemperatureResponse

  select case(itype)
!     CLM4.5 temperature response function
      case(TEMPERATURE_RESPONSE_FUNCTION_Q10)
          Ft = Q10 ** ((tc - 25.0d0) / 10.0d0)

!     CLM-CN
!     Equation: F_t = exp(308.56*(1/17.02 - 1/(T - 227.13)))
      case(TEMPERATURE_RESPONSE_FUNCTION_CLM4) 
  
          tk = tc + 273.15d0

          if(tk > 227.15d0) then
              Ft = exp(308.56d0*(one_over_71_02 - 1.d0/(tk - 227.13d0)))
          else
              Ft = 0.d0
          endif
!     DLEM temperature response function for methane oxidation
! Tian et al. 2010 Biogeosciences, 7, 2673-2694 Eq. 12
      case(TEMPERATURE_RESPONSE_FUNCTION_DLEM)
          if(tc < -5.0d0) then
            Ft = 0.0d0
          elseif(tc >= 30.0d0) then
            Ft = 1.0d0
          else
            Ft = Q10 ** ((tc - 30.0d0) / 10.0d0)
          endif
       case default
            Ft = 1.0d0
  end select
  GetTemperatureResponse = Ft 

end function GetTemperatureResponse
  
! ************************************************************************** !

Function GetMoistureResponse(thetapsi, local_id, itype)

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif
  
  implicit none
  
#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif

  PetscReal :: F_theta
! thetapsi IS either 'theta' (soil VWC: 0 - porosity) or 'psi' (matric potential, -Pa), upon 'moisture_response_function type'
  PetscReal :: thetapsi
  PetscReal :: GetMoistureResponse
  PetscReal, parameter :: theta_min = 0.01d0     ! 1/nat log(0.01d0)
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1
  PetscReal, parameter :: twelve_over_14 = 0.857142857143d0

  PetscInt :: local_id, itype
  PetscReal :: maxpsi, psi, theta, tc
  PetscReal, parameter :: minpsi = -10.0d6    ! Pa

#ifdef CLM_PFLOTRAN
  PetscScalar, pointer :: sucsat_pf_loc(:)    !
  PetscScalar, pointer :: watfc_pf_loc(:)     !
  PetscScalar, pointer :: watsat_pf_loc(:)    !
  PetscReal :: thetar, thetas, se
#endif

  PetscErrorCode :: ierr

#ifdef CLM_PFLOTRAN

  select case(itype)
!   CLM-CN
    case(MOISTURE_RESPONSE_FUNCTION_CLM4) 
      call VecGetArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)
      ! sucsat [mm of H20] from CLM is the suction (positive) at water saturated (called air-entry pressure)
      ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
      maxpsi = sucsat_pf_loc(local_id) * (-9.81d0)
      psi = min(thetapsi, maxpsi)                     ! thetapsi IS psi (-Pa)
      if(psi > minpsi) then
        F_theta = log(minpsi/psi)/log(minpsi/maxpsi)
      else
        F_theta = 0.0d0
      endif
      call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)

! DLEM 
! Tian et al. 2010 Biogeosciences, 7, 2673-2694 Eq. 13
    case(MOISTURE_RESPONSE_FUNCTION_DLEM) 
      call VecGetArrayReadF90(clm_pf_idata%watsat_pf, watsat_pf_loc, ierr)
      call VecGetArrayReadF90(clm_pf_idata%watfc_pf, watfc_pf_loc, ierr)
      thetas = watsat_pf_loc(local_id)
      thetar = watfc_pf_loc(local_id)
      theta = thetapsi                           ! thetapsi IS 'theta'
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
      call VecGetArrayReadF90(clm_pf_idata%watsat_pf, watsat_pf_loc, ierr)
      call VecGetArrayReadF90(clm_pf_idata%watfc_pf, watfc_pf_loc, ierr)
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

Function GetpHResponse(pH, itype)

  PetscReal :: f_ph, pH, GetpHResponse
  PetscInt :: itype

! ph function from Parton et al., (2001, 1996)
!  k_nitr_ph_vr(c,j) = 0.56 + atan(rpi * 0.45 * (-5.+ pH(c)))/rpi
#ifdef CLM_PFLOTRAN
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

end module CLM_BGC_module
