module Reaction_Sandbox_Denitrification_class

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
    extends(reaction_sandbox_base_type) :: reaction_sandbox_denitrification_type
    PetscInt :: ispec_no3
    PetscInt :: ispec_n2
    PetscInt :: ispec_n2o

!   for co2 respiration calculation
    PetscInt :: ispec_n
    PetscInt :: ispec_lit1c
    PetscInt :: ispec_lit2c
    PetscInt :: ispec_lit3c
    PetscInt :: ispec_lit1n
    PetscInt :: ispec_lit2n
    PetscInt :: ispec_lit3n
    PetscInt :: ispec_som1
    PetscInt :: ispec_som2
    PetscInt :: ispec_som3
    PetscInt :: ispec_som4

    PetscReal :: rate_constant
    PetscReal :: half_saturation
    PetscInt :: temperature_response_function
    PetscReal :: Q10

  contains
    procedure, public :: ReadInput => DenitrificationRead
    procedure, public :: Setup => DenitrificationSetup
    procedure, public :: Evaluate => DenitrificationReact
    procedure, public :: Destroy => DenitrificationDestroy
  end type reaction_sandbox_denitrification_type

  public :: DenitrificationCreate

contains

! ************************************************************************** !
!
! DenitrificationCreate: Allocates denitrification reaction object.
! author: Guoping Tang (replace in all subroutine headers with name of developer) 
! date: 09/09/2013 (replace in all subroutine headers with current date)
!
! ************************************************************************** !
function DenitrificationCreate()

  implicit none
  
  class(reaction_sandbox_denitrification_type), pointer :: DenitrificationCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(DenitrificationCreate)
  DenitrificationCreate%ispec_no3 = 0
  DenitrificationCreate%ispec_n2o = 0
  DenitrificationCreate%ispec_n2 = 0

  DenitrificationCreate%ispec_n = 0
  DenitrificationCreate%ispec_lit1c = 0
  DenitrificationCreate%ispec_lit2c = 0
  DenitrificationCreate%ispec_lit3c = 0
  DenitrificationCreate%ispec_lit1n = 0
  DenitrificationCreate%ispec_lit2n = 0
  DenitrificationCreate%ispec_lit3n = 0
  DenitrificationCreate%ispec_som1 = 0
  DenitrificationCreate%ispec_som2 = 0
  DenitrificationCreate%ispec_som3 = 0
  DenitrificationCreate%ispec_som4 = 0
  DenitrificationCreate%rate_constant = 0.d0
  DenitrificationCreate%half_saturation = 1.0d-10
  DenitrificationCreate%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4
  DenitrificationCreate%Q10 = 1.5d0
  nullify(DenitrificationCreate%next)  
      
end function DenitrificationCreate

! ************************************************************************** !
!
! DenitrificationRead: Reads input deck for denitrification reaction parameters (if any)
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine DenitrificationRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_denitrification_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION')
    call StringToUpper(word)   

    select case(trim(word))
      case('TEMPERATURE_RESPONSE_FUNCTION')
        do
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit

         call InputReadWord(input,option,word,PETSC_TRUE)
         call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,TEMPERATURE RESPONSE FUNCTION')
         call StringToUpper(word)   

         select case(trim(word))
              case('CLM4')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4    
              case('Q10')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_Q10    
                  call InputReadDouble(input,option,this%Q10)  
                  call InputErrorMsg(input,option,'Q10', &
                        'CHEMISTRY,REACTION_SANDBOX_DENITRIFICATION,TEMPERATURE RESPONSE FUNCTION')
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
        enddo 

      case('RATE_CONSTANT')
          call InputReadDouble(input,option,this%rate_constant)
          call InputErrorMsg(input,option,'rate constant', &
                 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,REACTION')
      case('N_INHIBITION')
          call InputReadDouble(input,option,this%half_saturation)
          call InputErrorMsg(input,option,'inhibition coefficient', &
                 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,REACTION')
      case default
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,' // &
            'REACTION keyword: ' // trim(word) // ' not recognized.'
          call printErrMsg(option)
    end select
  enddo
  
end subroutine DenitrificationRead

! ************************************************************************** !
!
! DenitrificationSetup: Sets up the denitrification reaction either with parameters either
!                read from the input deck or hardwired.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine DenitrificationSetup(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module
  use Immobile_Aux_module, only : GetImmobileSpeciesIDFromName 

  implicit none
  
  class(reaction_sandbox_denitrification_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
 
  word = 'Nitrate'
  this%ispec_no3 = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)
      
  word = 'N2O'
  this%ispec_n2o = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'N2'
  this%ispec_n2 = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'N'
  this%ispec_n = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'LabileC'
  this%ispec_lit1c = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'LabileN'
  this%ispec_lit1n = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'CelluloseC'
  this%ispec_lit2c = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'CelluloseN'
  this%ispec_lit2n = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'LigninC'
  this%ispec_lit3c = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'LigninN'
  this%ispec_lit3n = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'SOM1'
  this%ispec_som1 = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'SOM2'
  this%ispec_som2 = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'SOM3'
  this%ispec_som3 = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'SOM4'
  this%ispec_som4 = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

end subroutine DenitrificationSetup

! ************************************************************************** !
!
! DenitrificationReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 01/09/2014
!
! ************************************************************************** !
subroutine DenitrificationReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,porosity,volume,reaction, &
                         option,local_id)

  use Option_module
  use Reaction_Aux_module

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif
  
  implicit none

#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif
  
  class(reaction_sandbox_denitrification_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscBool :: compute_derivative

  ! the following arrays must be declared after reaction
  PetscInt :: local_id 
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: rate, rate_no3, drate_no3, rate_co2, drate_co2, rate_tmp
  PetscReal :: rate_no3_co2
  PetscReal :: volume, porosity, saturation
  PetscErrorCode     :: ierr

  ! inhibition variables
  PetscReal :: tc, tk
  PetscReal :: F_t
  PetscReal :: F_theta
  PetscReal :: F_ph
  PetscReal :: tmp_real 
  PetscReal :: anaerobic_frac
  PetscReal :: k_nitr_max
  PetscReal :: soil_bulkdensity  ! clm4.5.35 denitrification, dry bulk density + water

  PetscReal, parameter :: one_over_71_02 = 1.408054069d-2
  PetscReal, parameter :: theta_min = 0.01d0     ! 1/nat log(0.01d0)
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1
  PetscReal, parameter :: twelve_over_14 = 0.857142857143d0
  PetscReal, parameter :: minpsi = -10.0d0  ! MPa
  PetscReal, parameter :: C_molecular_weight = 12.0107d0
  PetscReal, parameter :: N_molecular_weight = 14.0067d0

  PetscReal :: diffus
  PetscReal :: h2osoi_vol
  PetscInt, parameter :: iphase = 1

#ifdef CLM_PFLOTRAN
  PetscReal, parameter :: rpi = 3.14159265358979323846

!  PetscReal, parameter :: WT_Saturation = 0.95 !volumetric soil water defining top of water table clm4.5
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
  PetscReal :: r_min, r_max, r_psi
  PetscReal :: rho_w, grav
  PetscReal :: ratio_diffusivity_water_gas

  PetscScalar, pointer :: sucsat(:)
  PetscScalar, pointer :: watfc(:)
  PetscScalar, pointer :: cellorg(:)
  PetscScalar, pointer :: bsw(:)
  PetscScalar, pointer :: soilpsi(:)
  PetscScalar, pointer :: bulkdensity_dry(:)
  PetscScalar, pointer :: o2_decomp_depth_unsat(:)
  PetscScalar, pointer :: conc_o2_unsat(:)
  PetscScalar, pointer :: o2_decomp_depth_sat(:)
  PetscScalar, pointer :: conc_o2_sat(:)
#endif

  ! soil respiration
  PetscInt :: ires_n2, ires_no3, ires_n2o
  PetscReal :: k_n
  PetscReal :: kg_water_m3
  PetscReal :: ratio_k1, ratio_no3_co2
  PetscReal :: klit1, klit2, klit3, ksom1, ksom2, ksom3, ksom4
  PetscReal :: conc_lit1c, conc_lit2c, conc_lit3c
  PetscReal :: conc_lit1n, conc_lit2n, conc_lit3n
  PetscReal :: conc_som1, conc_som2, conc_som3, conc_som4
  PetscReal :: conc_no3, conc_n
  PetscReal :: respiration_fraction_lit1, respiration_fraction_lit2
  PetscReal :: respiration_fraction_lit3, respiration_fraction_som1
  PetscReal :: respiration_fraction_som2, respiration_fraction_som3
  PetscReal :: respiration_fraction_som4
  PetscReal :: cn_downstream_lit1, cn_downstream_lit2, cn_downstream_lit3
  PetscReal :: mol_per_m3__to__ug_per_gsoil
  PetscReal :: fr_WFPS, wfps_vr
  PetscReal :: n_inhibition
  PetscReal :: n2_n2o_ratio, stoich_n2o, stoich_n2, stoich_n
  PetscReal ::  d_stoich_n2o_d_no3, d_stoich_n2_d_no3

  ! temperature response function 
  tc = global_auxvar%temp(1)
  tk = tc + 273.15d0
  if(tk <= 273.15d0) then
  ! limit to non-frozen soil layers
    return
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
  saturation = global_auxvar%sat(iphase) 
  h2osoi_vol = porosity * saturation  !? 
#ifdef CLM_PFLOTRAN


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
 

  call VecGetArrayReadF90(clm_pf_idata%sucsat_pf, sucsat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%watfc_pf, watfc, ierr)
  call VecGetArrayReadF90(clm_pf_idata%cellorg_pf, cellorg, ierr)
  call VecGetArrayReadF90(clm_pf_idata%bsw_pf, bsw, ierr)
  call VecGetArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi, ierr)
  call VecGetArrayReadF90(clm_pf_idata%o2_decomp_depth_unsat_pf, o2_decomp_depth_unsat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%o2_decomp_depth_sat_pf, o2_decomp_depth_sat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%conc_o2_unsat_pf, conc_o2_unsat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%conc_o2_sat_pf, conc_o2_sat, ierr)

  f_a = 1.0 - watfc(local_id) / porosity
  e_a = porosity - watfc(local_id)

  if (clm_pf_idata%use_lch4) then
     if (organic_max > 0.0) then
        om_frac = min(cellorg(local_id)/organic_max, 1.0)
     else
        om_frac = 1.0
     end if

     diffus = (d_con_g_1_o2 + d_con_g_2_o2*tk) * 1.d-4 * &
           (om_frac * f_a**(10.0/3.0) / porosity**2 + &
           (1.0 - om_frac) * e_a**2 * f_a**(3.0 / bsw(local_id)))

  ! calculate anoxic fraction of soils
  ! use rijtema and kroess model after Riley et al., 2000
  ! caclulated r_psi as a function of psi

!     if(saturation < WT_saturation) then
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
                       porosity**rij_kro_delta))

        else
           anaerobic_frac = 0.0
           diffus = 1.d-6
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
 !                      conc_o2_sat(local_id)**rij_kro_gamma * (porosity +  &
 !                      ratio_diffusivity_water_gas * porosity)**rij_kro_delta)
!             anaerobic_frac = 0.0
!         else
!             anaerobic_frac = 0.0
!         endif
!               anaerobic_frac(c,j) = (1._r8 - finundated(c))*anaerobic_frac(c,j) + finundated(c)*anaerobic_frac_sat
!     end if
  else
      anaerobic_frac = 0.1
      diffus = 1.d-6
  endif

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
  anaerobic_frac = 0.1
  diffus = 1.d-6
#endif

!  anaerobic_frac = 1.0

!---------------- denitrification
! follows CENTURY denitrification scheme (Parton et al., (2001, 1996))

! assume denitrification temp function equal to the HR scalar
!  k_nitr_t_vr(c,j) = min(t_scalar(c,j), 1._r8)
  F_t = min(F_t, 1.0)
 
#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%bulkdensity_dry_pf, bulkdensity_dry, ierr)
  soil_bulkdensity = bulkdensity_dry(local_id)    
  call VecRestoreArrayReadF90(clm_pf_idata%bulkdensity_dry_pf, bulkdensity_dry, ierr)
#else
  soil_bulkdensity = 1.0d0
#endif

  kg_water_m3 = porosity*global_auxvar%sat(iphase)*global_auxvar%den_kg(iphase)
  soil_bulkdensity = soil_bulkdensity + kg_water_m3

  mol_per_m3__to__ug_per_gsoil = 1.e3 / soil_bulkdensity * N_molecular_weight
  
!! maximum potential denitrification rates based on heterotrophic respiration rates or nitrate concentrations, 
!! from (del Grosso et al., 2000)
! rate due to no3
  conc_no3 = rt_auxvar%immobile(this%ispec_no3)
  rate_no3 = 1.15 * (conc_no3 * mol_per_m3__to__ug_per_gsoil)**0.57  ! ug/(g soil day)
  rate_no3 = rate_no3/mol_per_m3__to__ug_per_gsoil/86400.0 

  if(this%half_saturation > 1.0d-20) then
     tmp_real  = this%half_saturation + conc_no3
     rate_no3  = rate_no3 * conc_no3 / tmp_real 
!     rate_no3 = 1.15 * mol_per_m3__to__ug_per_gsoil)**(-0.43)/86400.0 * conc_no3**1.57/(k + c)
! d(x^n/(k+x))/dx = [nx^(n-1)(k+x) - x^n)]/(k + x)^2

     drate_no3 = 1.15 * mol_per_m3__to__ug_per_gsoil**(-0.43)/86400.0 * &
                (1.57 * conc_no3**0.57 * tmp_real - conc_no3**1.57)/tmp_real/tmp_real
  else
!     rate_no3 = 1.15 * mol_per_m3__to__ug_per_gsoil)**(-0.43)/86400.0 * conc_no3**0.57
     drate_no3 = 1.15 * mol_per_m3__to__ug_per_gsoil**(-0.43)/86400.0 * 0.57 / conc_no3**0.43
  endif

! respiration rate /co2 production rate from c decomposition
! hardwired for CLM-CN
!  
  conc_n = rt_auxvar%immobile(this%ispec_n)
  k_n = 1.0d-15
  n_inhibition = conc_n/(conc_n + k_n)

  rate_co2 = 0.0
  rate_tmp = 0.0
!  drate_co2_dn = 0.0

! lit1
  conc_lit1c = rt_auxvar%immobile(this%ispec_lit1c)
  klit1 = 1.204/86400.0
  respiration_fraction_lit1 = 0.39
  cn_downstream_lit1 = 12.0
  if(conc_lit1c > 1.0d-20) then
     rate_tmp = klit1 * conc_lit1c * respiration_fraction_lit1
     conc_lit1n = rt_auxvar%immobile(this%ispec_lit1n)
     stoich_n   = conc_lit1n / conc_lit1c - (1.0 - respiration_fraction_lit1) &
                * C_molecular_weight / cn_downstream_lit1 / N_molecular_weight  

     if(stoich_n < 0.0) then
        rate_co2 = rate_co2 + rate_tmp * n_inhibition
     else
        rate_co2 = rate_co2 + rate_tmp     
     endif     
  endif

  conc_lit2c = rt_auxvar%immobile(this%ispec_lit2c)
  klit2 = 0.0726/86400.0
  respiration_fraction_lit2 = 0.55
  cn_downstream_lit2 = 12.0      
  if(conc_lit2c > 1.0d-20) then
     rate_tmp = klit2 * conc_lit2c * respiration_fraction_lit2
     conc_lit2n = rt_auxvar%immobile(this%ispec_lit2n)
     stoich_n   = conc_lit2n / conc_lit2c - (1.0 - respiration_fraction_lit2) &
                * C_molecular_weight / cn_downstream_lit2 / N_molecular_weight  
     if(stoich_n < 0.0) then
        rate_co2 = rate_co2 + rate_tmp * n_inhibition     
     else
        rate_co2 = rate_co2 + rate_tmp     
     endif     
  endif

  conc_lit3c = rt_auxvar%immobile(this%ispec_lit3c)
  klit3 = 0.0141/86400.0
  respiration_fraction_lit3 = 0.29
  cn_downstream_lit3 = 10.0 
   
  if(conc_lit3c > 1.0d-20) then
     rate_tmp = klit3 * conc_lit3c * respiration_fraction_lit3
     conc_lit3n = rt_auxvar%immobile(this%ispec_lit3n)
     stoich_n   = conc_lit3n / conc_lit3c - (1.0 - respiration_fraction_lit3) &
                * C_molecular_weight / cn_downstream_lit3 / N_molecular_weight  

     if(stoich_n < 0.0) then
        rate_co2 = rate_co2 + rate_tmp * n_inhibition    
     else
        rate_co2 = rate_co2 + rate_tmp    
     endif     
  endif

  conc_som1 = rt_auxvar%immobile(this%ispec_som1)
  ksom1 = 0.0726/86400.0
  respiration_fraction_som1 = 0.28  
  rate_co2 = rate_co2 + ksom1 * conc_som1 * respiration_fraction_som1     

  conc_som2 = rt_auxvar%immobile(this%ispec_som2)
  ksom2 = 0.0141/86400.0
  respiration_fraction_som2 = 0.46  
  rate_co2 = rate_co2 + ksom2 * conc_som2 * respiration_fraction_som2     

  conc_som3 = rt_auxvar%immobile(this%ispec_som3)
  ksom3 = 0.007/86400.0
  respiration_fraction_som3 = 0.55  
  rate_co2 = rate_co2 + ksom3 * conc_som3 * respiration_fraction_som3     

  conc_som4 = rt_auxvar%immobile(this%ispec_som4)
  ksom4 = 0.007/86400.0
  rate_co2 = rate_co2 + ksom4 * conc_som4     

  rate_co2 = rate_co2 * F_t * F_theta

  rate_no3_co2 = 0.1 * (rate_co2 * mol_per_m3__to__ug_per_gsoil*86400.0)**1.3/mol_per_m3__to__ug_per_gsoil/86400.0

! now calculate the ratio of N2O to N2 from denitrifictaion, following Del Grosso et al., 2000
! diffusivity constant (figure 6b)
  ratio_k1 = max(1.7, 38.4 - 350.0 * diffus)

! ratio function (figure 7c)
  if ( rate_co2 > 0.0 ) then
       ratio_no3_co2 = conc_no3 / rate_co2 / 86400.0    !  not unitless?
  else
! fucntion saturates at large no3/co2 ratios, so set as some nominally large number
       ratio_no3_co2 = 100.0
  endif

! total water limitation function (Del Grosso et al., 2000, figure 7a)
  wfps_vr = max(min(h2osoi_vol/porosity, 1.0), 0.0) * 100.0
  fr_WFPS = max(0.1, 0.015 * wfps_vr - 0.32)
!         if (use_lch4) then
!            if (anoxia_wtsat) then
!               fr_WFPS(c,j) = fr_WFPS(c,j)*(1._r8 - finundated(c)) + finundated(c)*1.18_r8
!            end if
!         end if

         ! final ratio expression 
  n2_n2o_ratio = ratio_k1 * max(0.16, exp(-0.8 * ratio_no3_co2)) * fr_WFPS
  stoich_n2 = n2_n2o_ratio/(1.0 + n2_n2o_ratio)
  stoich_n2o = 1.0 - stoich_n2

  if(rate_co2 > 0.0) then
     d_stoich_n2_d_no3 = 1.0 / (1.0 + n2_n2o_ratio) / (1.0 + n2_n2o_ratio) * &
            ratio_k1 * fr_WFPS * (-0.8 / rate_co2 ) * exp (-0.8*ratio_no3_co2)
     d_stoich_n2o_d_no3 = 1.0 - d_stoich_n2_d_no3
  endif

  ires_no3 = this%ispec_no3 + reaction%offset_immobile      
  ires_n2o = this%ispec_n2o + reaction%offset_immobile      
  ires_n2  = this%ispec_n2  + reaction%offset_immobile      

  if(rate_no3 > rate_no3_co2) then
     rate = rate_no3_co2 * anaerobic_frac
  else
     rate = rate_no3 * anaerobic_frac
  endif 

  Residual(ires_no3) = Residual(ires_no3) - (-1.0) * rate
  Residual(ires_n2o) = Residual(ires_n2o) - rate * stoich_n2o
  Residual(ires_n2) = Residual(ires_n2) - rate * stoich_n2

  if (compute_derivative) return

    ! always add contribution to Jacobian
  tmp_real =  exp(-0.8 * ratio_no3_co2)
  if(rate_no3 < rate_no3_co2) then
    drate_no3 = drate_no3 * anaerobic_frac
    Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) - drate_no3
    if(tmp_real > 0.16) then
       Jacobian(ires_n2o,ires_no3) = Jacobian(ires_n2o,ires_no3) &
                                   + drate_no3 * stoich_n2o &
                                   + rate * d_stoich_n2_d_no3 
       Jacobian(ires_n2,ires_no3)  = Jacobian(ires_n2,ires_no3) &
                                   + drate_no3 * stoich_n2 &
                                   + rate * d_stoich_n2o_d_no3
    else
       Jacobian(ires_n2o,ires_no3) = Jacobian(ires_n2o,ires_no3) + drate_no3 * stoich_n2o
       Jacobian(ires_n2,ires_no3) = Jacobian(ires_n2,ires_no3) + drate_no3 * stoich_n2
    endif
  else
    if(tmp_real > 0.16) then
       Jacobian(ires_n2o,ires_no3) = Jacobian(ires_n2o,ires_no3) + rate * d_stoich_n2_d_no3 
       Jacobian(ires_n2,ires_no3) = Jacobian(ires_n2,ires_no3) + rate * d_stoich_n2o_d_no3
    endif
  endif

end subroutine DenitrificationReact

! ************************************************************************** !
!
! DenitrificationDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine DenitrificationDestroy(this)

  implicit none
  
  class(reaction_sandbox_denitrification_type) :: this  

end subroutine DenitrificationDestroy

end module Reaction_Sandbox_Denitrification_class
