module Reaction_Sandbox_Nitrification_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
#include "finclude/petscsys.h"

  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_CLM4 = 1 
  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_Q10 = 2 

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_nitrification_type
    PetscInt :: ispec_nh3
    PetscInt :: ispec_no3
    PetscReal :: rate_constant
    PetscReal :: half_saturation
    PetscInt :: temperature_response_function
    PetscReal :: Q10

  contains
    procedure, public :: ReadInput => NitrificationRead
    procedure, public :: Setup => NitrificationSetup
    procedure, public :: Evaluate => NitrificationReact
    procedure, public :: Destroy => NitrificationDestroy
  end type reaction_sandbox_nitrification_type

  public :: NitrificationCreate

contains

! ************************************************************************** !
!
! NitrificationCreate: Allocates nitrification reaction object.
! author: Guoping Tang (replace in all subroutine headers with name of developer) 
! date: 09/09/2013 (replace in all subroutine headers with current date)
!
! ************************************************************************** !
function NitrificationCreate()

  implicit none
  
  class(reaction_sandbox_nitrification_type), pointer :: NitrificationCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(NitrificationCreate)
  NitrificationCreate%ispec_nh3 = 0
  NitrificationCreate%ispec_no3 = 0
  NitrificationCreate%rate_constant = 0.d0
  NitrificationCreate%half_saturation = 1.0d-10
  NitrificationCreate%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4
  NitrificationCreate%Q10 = 1.5d0
  nullify(NitrificationCreate%next)  
      
end function NitrificationCreate

! ************************************************************************** !
!
! NitrificationRead: Reads input deck for nitrification reaction parameters (if any)
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine NitrificationRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_nitrification_type) :: this
  type(input_type) :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION')
    call StringToUpper(word)   

    select case(trim(word))
      case('TEMPERATURE_RESPONSE_FUNCTION')
        do
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit

         call InputReadWord(input,option,word,PETSC_TRUE)
         call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION,TEMPERATURE RESPONSE FUNCTION')
         call StringToUpper(word)   

         select case(trim(word))
              case('CLM4')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4    
              case('Q10')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_Q10    
                  call InputReadDouble(input,option,this%Q10)  
                  call InputErrorMsg(input,option,'Q10', &
                        'CHEMISTRY,REACTION_SANDBOX_NITRIFICATION,TEMPERATURE RESPONSE FUNCTION')
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION,TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
        enddo 

      case('RATE_CONSTANT')
          call InputReadDouble(input,option,this%rate_constant)
          call InputErrorMsg(input,option,'rate constant', &
                 'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION,REACTION')
      case('N_INHIBITION')
          call InputReadDouble(input,option,this%half_saturation)
          call InputErrorMsg(input,option,'inhibition coefficient', &
                 'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION,REACTION')
      case default
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION,' // &
            'REACTION keyword: ' // trim(word) // ' not recognized.'
          call printErrMsg(option)
    end select
  enddo
  
end subroutine NitrificationRead

! ************************************************************************** !
!
! NitrificationSetup: Sets up the nitrification reaction either with parameters either
!                read from the input deck or hardwired.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine NitrificationSetup(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module
  use Immobile_Aux_module, only : GetImmobileSpeciesIDFromName 

  implicit none
  
  class(reaction_sandbox_nitrification_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
 
  word = 'N'
  this%ispec_nh3 = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'Nitrate'
  this%ispec_no3 = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)
      
end subroutine NitrificationSetup

! ************************************************************************** !
!
! NitrificationReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine NitrificationReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,porosity,volume,reaction, &
                         option,local_id)

  use Option_module
  use Reaction_Aux_module

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data !, only : rate_ndeni_decomp_pf 
#endif
  
  implicit none

#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif
  
  class(reaction_sandbox_nitrification_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscBool :: compute_derivative

  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: rate, drate, concN, rate0
  PetscReal :: volume, porosity, saturation
  PetscInt :: local_id, ires_nh3, ires_no3
  PetscErrorCode     :: ierr

  ! inhibition variables
  PetscReal :: tc, tk
  PetscReal :: F_t
  PetscReal :: F_theta
  PetscReal :: F_ph
  PetscReal :: tmp_real 
  PetscReal :: anaerobic_frac
  PetscReal :: k_nitr_max

  PetscReal, parameter :: one_over_71_02 = 1.408054069d-2
  PetscReal, parameter :: theta_min = 0.01d0     ! 1/nat log(0.01d0)
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1
  PetscReal, parameter :: twelve_over_14 = 0.857142857143d0
  PetscReal, parameter :: minpsi = -10.0d0  ! MPa

#ifdef CLM_PFLOTRAN
  PetscReal, parameter :: rpi = 3.14159265358979323846

  PetscInt, parameter :: iphase = 1
  PetscReal, parameter :: WT_Saturation = 0.95 !volumetric soil water defining top of water table clm4.5
  PetscReal, parameter :: spval = 1.0d36
 
  PetscReal :: psi, maxpsi

  PetscReal :: surface_tension_water
  PetscReal :: rij_kro_a                  !  Arah and Vinten 1995
  PetscReal :: rij_kro_alpha              !  Arah and Vinten 1995
  PetscReal :: rij_kro_beta               !  Arah and Vinten 1995
  PetscReal :: rij_kro_gamma              !  Arah and Vinten 1995
  PetscReal :: rij_kro_delta              !  Arah and Vinten 1995
  PetscReal :: organic_max                ! organic matter content (kg/m3) where soil is assumed to act like peat
  PetscReal :: pH                         ! placeholder
  PetscReal :: co2diff_con1, co2diff_con2 ! diffusion constants for CO2
 
  PetscReal :: d_con_w_1_o2, d_con_w_2_o2, d_con_w_3_o2           ! water diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-4)
  PetscReal :: d_con_g_1_o2, d_con_g_2_o2           ! gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)

  PetscReal :: f_a   !
  PetscReal :: e_a   ! air filled fraction of total soil volume
  PetscReal :: om_frac
  PetscReal :: diffus
  PetscReal :: r_min, r_max, r_psi
  PetscReal :: rho_w, grav
  PetscReal :: ratio_diffusivity_water_gas
  PetscReal :: h2osoi_vol

  PetscScalar, pointer :: watsat(:)
  PetscScalar, pointer :: sucsat(:)
  PetscScalar, pointer :: watfc(:)
  PetscScalar, pointer :: cellorg(:)
  PetscScalar, pointer :: bsw(:)
  PetscScalar, pointer :: soilpsi(:)
  PetscScalar, pointer :: o2_decomp_depth_unsat(:)
  PetscScalar, pointer :: conc_o2_unsat(:)
  PetscScalar, pointer :: o2_decomp_depth_sat(:)
  PetscScalar, pointer :: conc_o2_sat(:)
#endif

  ! temperature response function 
  tc = global_auxvar%temp(1)
  tk = tc + 273.15d0
  if(tk <= 273.15d0) then
  ! limit to non-frozen soil layers
    return
!  else
!    print *, 'tc = ', tc 
  endif

  select case(this%temperature_response_function)
      case(TEMPERATURE_RESPONSE_FUNCTION_Q10)
! CLM4.5 temperature response function
          F_t = this%Q10 ** ((tc - 25.0d0) / 10.0d0)
      case(TEMPERATURE_RESPONSE_FUNCTION_CLM4) 
  
  ! inhibition due to temperature
  ! Equation: F_t = exp(308.56*(1/17.02 - 1/(T - 227.13)))

          if(tk > 227.15d0) then
             F_t = exp(308.56d0*(one_over_71_02 - 1.d0/(tk - 227.13d0)))
          else
             F_t = 0.0
             return
          endif
  end select

!moisture response function
#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%sucsat_pf, sucsat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi, ierr)

  maxpsi = sucsat(local_id) * (-9.8d-6)
  psi = min(soilpsi(local_id), maxpsi)

  if(psi > minpsi) then
     F_theta = log(minpsi/psi)/log(minpsi/maxpsi)
  else
     F_theta = 0.0d0
     call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat, ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi, ierr)
     return
  endif

  call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi, ierr)
#else

  ! inhibition due to moisture content
  ! Equation: F_theta = log(theta_min/theta) / log(theta_min/theta_max)
  ! Assumptions: theta is saturation
  !              theta_min = 0.01, theta_max = 1.
  F_theta = log(theta_min/max(theta_min,global_auxvar%sat(1))) * one_over_log_theta_min 
  
#endif

! anaerobic fraction
#ifdef CLM_PFLOTRAN
  saturation = global_auxvar%sat(iphase) 
  h2osoi_vol = porosity * saturation  !? 


  surface_tension_water = 73.d-3   ! (J/m^2), Arah and Vinten 1995

  ! Set parameters from simple-structure model to calculate anoxic fratction (Arah and Vinten 1995)
  rij_kro_a = 1.5d-10              !  Arah and Vinten 1995
  rij_kro_alpha = 1.26             !  Arah and Vinten 1995
  rij_kro_beta = 0.6               !  Arah and Vinten 1995
  rij_kro_gamma = 0.6              !  Arah and Vinten 1995
  rij_kro_delta = 0.85             !  Arah and Vinten 1995

  organic_max  = 130.0             ! organic matter content (kg/m3) where soil is assumed to act like peat
  ! for diffusion. Very large values will lead to all soil being treated as mineral. Negative values will lead   ! to all soil being treated as peat.


  pH = 6.5  !!! set all soils with the same pH as placeholder here
  co2diff_con1 =   0.1325
  co2diff_con2 =   0.0009

!  data (d_con_w(1,i),i=1,3) /0.9798_r8, 0.02986_r8, 0.0004381_r8/ ! CH4
!  data (d_con_w(2,i),i=1,3) /1.172_r8, 0.03443_r8, 0.0005048_r8/ ! O2
!  data (d_con_w(3,i),i=1,3) /0.939_r8, 0.02671_r8, 0.0004095_r8/ ! CO2

  
  d_con_g_1_o2 = 0.1759
  d_con_g_2_o2 = 0.00117

  d_con_w_1_o2 = 1.172
  d_con_w_2_o2 = 0.03443
  d_con_w_3_o2 = 0.0005048

  rho_w  = 1.d3                   ! (kg/m3)
  grav   = 9.80616                ! acceleration of gravity ~ m/s^2
 

  call VecGetArrayReadF90(clm_pf_idata%watsat_pf, watsat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%sucsat_pf, sucsat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%watfc_pf, watfc, ierr)
  call VecGetArrayReadF90(clm_pf_idata%cellorg_pf, cellorg, ierr)
  call VecGetArrayReadF90(clm_pf_idata%bsw_pf, bsw, ierr)
  call VecGetArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi, ierr)
  call VecGetArrayReadF90(clm_pf_idata%o2_decomp_depth_unsat_pf, o2_decomp_depth_unsat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%o2_decomp_depth_sat_pf, o2_decomp_depth_sat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%conc_o2_unsat_pf, conc_o2_unsat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%conc_o2_sat_pf, conc_o2_sat, ierr)

  f_a = 1.0 - watfc(local_id) / watsat(local_id)
  e_a = watsat(local_id) - watfc(local_id)

  if (clm_pf_idata%use_lch4) then
     if (organic_max > 0.0) then
        om_frac = min(cellorg(local_id)/organic_max, 1.0)
     else
        om_frac = 1.0
     end if

     diffus = (d_con_g_1_o2 + d_con_g_2_o2*tk) * 1.d-4 * &
           (om_frac * f_a**(10.0/3.0) / watsat(local_id)**2 + &
           (1.0 - om_frac) * e_a**2 * f_a**(3.0 / bsw(local_id)))

  ! calculate anoxic fraction of soils
  ! use rijtema and kroess model after Riley et al., 2000
  ! caclulated r_psi as a function of psi

!     if(saturation < WT_saturation) then
        tmp_real = soilpsi(local_id)
        r_min = 2.0 * surface_tension_water / (rho_w * grav * abs(soilpsi(local_id)))
        r_max = 2.0 * surface_tension_water / (rho_w * grav * 0.1)
        r_psi = sqrt(r_min * r_max)

        ratio_diffusivity_water_gas = (d_con_g_1_o2 + d_con_g_2_o2*tk) * 1.d-4 / &
             ((d_con_w_1_o2 + d_con_w_2_o2*tk + d_con_w_3_o2*tk**2) * 1.d-9)

        if (o2_decomp_depth_unsat(local_id) .ne. spval .and. &
           conc_o2_unsat(local_id) .ne. spval .and. &
           o2_decomp_depth_unsat(local_id) > 0.0) then
           anaerobic_frac = exp(-rij_kro_a * r_psi**(-rij_kro_alpha) * &
                       o2_decomp_depth_unsat(local_id)**(-rij_kro_beta) * &
                       conc_o2_unsat(local_id)**rij_kro_gamma * (h2osoi_vol + &
                       ratio_diffusivity_water_gas * &
                       watsat(local_id))**rij_kro_delta)

        else
           anaerobic_frac = 0.0
        endif
!     else
! anoxia_wtsat = .false by default, NaN in o2_decomp_depth_sat
!     if (anoxia_wtsat) then ! Average saturated fraction values into anaerobic_frac(c,j).
!         r_min = 2.0 * surface_tension_water / (rho_w * grav * abs(grav * 1.e-6 * sucsat(local_id)))
!         r_max = 2.0 * surface_tension_water / (rho_w * grav * 0.1)
!         r_psi = sqrt(r_min * r_max)
!         ratio_diffusivity_water_gas = (d_con_g_1_o2 + d_con_g_2_o2*tk) * 1.d-4 / &
!             ((d_con_w_1_o2 + d_con_w_2_o2*tk + d_con_w_3_o2*tk**2) * 1.d-9)

 !        if (o2_decomp_depth_sat(local_id) .ne. spval .and. &
 !            conc_o2_sat(local_id) .ne. spval .and. &
 !            o2_decomp_depth_sat(local_id) > 0.0) then
 !            anaerobic_frac = exp(-rij_kro_a * r_psi**(-rij_kro_alpha) * &
 !                      o2_decomp_depth_sat(local_id)**(-rij_kro_beta) * &
 !                      conc_o2_sat(local_id)**rij_kro_gamma * (watsat(local_id) +  &
 !                      ratio_diffusivity_water_gas * watsat(local_id))**rij_kro_delta)
!             anaerobic_frac = 0.0
!         else
!             anaerobic_frac = 0.0
!         endif
!               anaerobic_frac(c,j) = (1._r8 - finundated(c))*anaerobic_frac(c,j) + finundated(c)*anaerobic_frac_sat
!     end if
  else
      anaerobic_frac = 0.0
  endif

  call VecRestoreArrayReadF90(clm_pf_idata%watsat_pf, watsat, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%watfc_pf, watfc, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%cellorg_pf, cellorg, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%bsw_pf, bsw, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%o2_decomp_depth_unsat_pf, o2_decomp_depth_unsat, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%o2_decomp_depth_sat_pf, o2_decomp_depth_sat, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%conc_o2_unsat_pf, conc_o2_unsat, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%conc_o2_sat_pf, conc_o2_sat, ierr)
#else
  anaerobic_frac = 1.0
#endif

!  anaerobic_frac = 1.0

!---------------- nitrification
! follows CENTURY nitrification scheme (Parton et al., (2001, 1996))

! assume nitrification temp function equal to the HR scalar
!  k_nitr_t_vr(c,j) = min(t_scalar(c,j), 1._r8)
  F_t = min(F_t, 1.0)
 
! ph function from Parton et al., (2001, 1996)
!  k_nitr_ph_vr(c,j) = 0.56 + atan(rpi * 0.45 * (-5.+ pH(c)))/rpi
#ifdef CLM_PFLOTRAN
  F_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + pH))/rpi
#else
  F_ph = 1.0
#endif
         ! moisture function-- assume the same moisture function as limits heterotrophic respiration
         ! Parton et al. base their nitrification- soil moisture rate constants based on heterotrophic rates-- can we do the same?
!         k_nitr_h2o_vr(c,j) = w_scalar(c,j)

  ! Set maximum nitrification rate constant 
  k_nitr_max =  0.1 / 86400.0   ! [1/sec] 10%/day  Parton et al., 2001 
  ! Todo:  SPM - the explicit divide gives different results than when that
  ! value is placed in the parameters netcdf file.  To get bfb, keep the 
  ! divide in source.
  !k_nitr_max = CNNitrifDenitrifParamsInst%k_nitr_max

!  k_nitr_max =  0.1 / 86400.0   ! [1/sec] 10%/day  Parton et al., 2001 
! nitrification constant is a set scalar * temp, moisture, and ph scalars
!         k_nitr_vr(c,j) = k_nitr_max * k_nitr_t_vr(c,j) * k_nitr_h2o_vr(c,j) * k_nitr_ph_vr(c,j)

  concN = rt_auxvar%immobile(this%ispec_nh3)
  rate = k_nitr_max * F_t * F_ph * F_theta * anaerobic_frac * concN

  if(this%half_saturation > 1.0d-20) then
!    f = k N N/(s + N)
!    df/dN = k (2N(s + N) - N^2)/(s + N)^2 = k (2Ns + N^2)/(s + N)^2
!                                          = k N (2s + N)/(S + N)^2
     tmp_real = this%half_saturation + concN
     drate =  rate  * (2.0 * this%half_saturation + concN) / tmp_real / tmp_real
     rate = rate * concN / (concN + this%half_saturation) 
  else
     drate = k_nitr_max * F_t * F_ph * F_theta * anaerobic_frac
  endif

!         ! limit to non-frozen soil layers
!         if ( t_soisno(c,j) <= SHR_CONST_TKFRZ .and. no_frozen_nitrif_denitrif) then
!            pot_f_nit_vr(c,j) = 0._r8
!         endif

  ires_nh3 = this%ispec_nh3 + reaction%offset_immobile      
  ires_no3 = this%ispec_no3 + reaction%offset_immobile      
  Residual(ires_nh3) = Residual(ires_nh3) - (-1.0) * rate
  Residual(ires_no3) = Residual(ires_no3) - rate

  if (compute_derivative) return

    ! always add contribution to Jacobian
  Jacobian(ires_nh3,ires_nh3) = &
      Jacobian(ires_nh3,ires_nh3) - drate

  Jacobian(ires_no3,ires_nh3) = Jacobian(ires_no3,ires_nh3) + drate

end subroutine NitrificationReact

! ************************************************************************** !
!
! NitrificationDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine NitrificationDestroy(this)

  implicit none
  
  class(reaction_sandbox_nitrification_type) :: this  

end subroutine NitrificationDestroy

end module Reaction_Sandbox_Nitrification_class
