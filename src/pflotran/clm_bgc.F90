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
  end select
  GetTemperatureResponse = Ft 

end function GetTemperatureResponse
  
! ************************************************************************** !

Function GetMoistureResponse(theta, local_id)

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif
  
  implicit none
  
#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif

  PetscReal :: F_theta, theta
  PetscReal :: GetMoistureResponse
  PetscReal, parameter :: theta_min = 0.01d0     ! 1/nat log(0.01d0)
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1
  PetscReal, parameter :: twelve_over_14 = 0.857142857143d0

  PetscInt :: local_id
  PetscReal :: maxpsi, psi, tc
  PetscReal, parameter :: minpsi = -10.0d0  ! MPa

#ifdef CLM_PFLOTRAN
  PetscScalar, pointer :: sucsat_pf_loc(:)   !
  PetscScalar, pointer :: soilpsi_pf_loc(:)   !
#endif

  PetscErrorCode :: ierr


#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)
  call VecGetArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi_pf_loc, ierr)

  maxpsi = sucsat_pf_loc(local_id) * (-9.8d-6)
  psi = min(soilpsi_pf_loc(local_id), maxpsi)

  if(psi > minpsi) then
     F_theta = log(minpsi/psi)/log(minpsi/maxpsi)
  else
     F_theta = 0.0d0
     call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi_pf_loc, ierr)
  endif

  call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi_pf_loc, ierr)

#else

  ! inhibition due to moisture content
  ! Equation: F_theta = log(theta_min/theta) / log(theta_min/theta_max)
  ! Assumptions: theta is saturation
  !              theta_min = 0.01, theta_max = 1.
  F_theta = 1.0d0 !log(theta_min/max(theta_min,theta)) * one_over_log_theta_min 
  
#endif

  GetMoistureResponse = F_theta

end function GetMoistureResponse

Function GetpHResponse(pH)

  PetscReal :: f_ph, pH, GetpHResponse

! ph function from Parton et al., (2001, 1996)
!  k_nitr_ph_vr(c,j) = 0.56 + atan(rpi * 0.45 * (-5.+ pH(c)))/rpi
#ifdef CLM_PFLOTRAN
  f_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + pH))/rpi
#else
  f_ph = 1.0
#endif

  GetpHResponse = f_ph

end function GetpHResponse

end module CLM_BGC_module
