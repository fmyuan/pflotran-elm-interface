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
    PetscInt :: ispec_proton
    PetscInt :: ispec_nh3
    PetscInt :: ispec_nh4sorb
    PetscInt :: ispec_no3
    PetscInt :: ispec_n2o
    PetscInt :: ispec_ngasnit
    PetscReal :: k_nitr_max
    PetscReal :: k_nitr_n2o
    PetscReal :: half_saturation
    PetscInt :: temperature_response_function
    PetscReal :: Q10
    PetscReal :: x0eps

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
  NitrificationCreate%ispec_proton = 0
  NitrificationCreate%ispec_nh3 = 0
  NitrificationCreate%ispec_nh4sorb = 0
  NitrificationCreate%ispec_no3 = 0
  NitrificationCreate%ispec_ngasnit = 0
  NitrificationCreate%k_nitr_max = 1.d-6
  NitrificationCreate%k_nitr_n2o = 3.5d-8
  NitrificationCreate%half_saturation = 1.0d-10
  NitrificationCreate%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4
  NitrificationCreate%Q10 = 1.5d0
  NitrificationCreate%x0eps = 1.0d-20
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
     case('X0EPS')
        call InputReadDouble(input,option,this%x0eps)
        call InputErrorMsg(input,option,'x0eps', &
                  'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION,REACTION')
     case('NITRIFICATION_RATE_COEF')
         call InputReadDouble(input,option,this%k_nitr_max)
         call InputErrorMsg(input,option,'nitrification rate coefficient', &
                     'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION,REACTION')
     case('N2O_RATE_COEF_NITRIFICATION')
         call InputReadDouble(input,option,this%k_nitr_n2o)
         call InputErrorMsg(input,option,'N2O rate coefficient from nirification', &
                     'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION,REACTION')
     case('AMMONIA_HALF_SATURATION')
          call InputReadDouble(input,option,this%half_saturation)
          call InputErrorMsg(input,option,'ammonia half-saturation', &
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

  word = 'H+'
  this%ispec_proton = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'NH4+'
  this%ispec_nh3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if(this%ispec_nh3 > 0) then
     word = 'NH4sorb'   ! this is the immobile species from 'reaction_sandbox_langmuir'
     this%ispec_nh4sorb = GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                        PETSC_FALSE,option)
  endif

  word = 'NO3-'
  this%ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'N2O(aq)'
  this%ispec_n2o = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if(this%ispec_nh3 < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION: ' // &
                        ' NH4+ is not specified in the input file.'
     call printErrMsg(option)
  endif

  if(this%ispec_no3 < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION: ' // &
                        ' NO3- is not specified in the input file.'
     call printErrMsg(option)
  endif

  if(this%ispec_n2o < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION: ' // &
                        ' N2O(aq) is not specified in the input file.'
     call printErrMsg(option)
  endif

  word = 'NGASnitr'
  this%ispec_ngasnit = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
end subroutine NitrificationSetup

! ************************************************************************** !
!
! NitrificationReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 09/09/2013
!
! Rewritten by Fengming Yuan @Aug-14-2014. The orginal was messed-up with 'patches',
! which caused a lot of issues.
!
! ************************************************************************** !
subroutine NitrificationReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only : material_auxvar_type

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
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
  class(material_auxvar_type) :: material_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscReal :: saturation
  PetscReal :: tc
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  PetscInt, parameter :: iphase = 1
  PetscReal, parameter :: rpi = 3.14159265358979323846
  PetscReal, parameter :: N_molecular_weight = 14.0067d0
  PetscReal :: M_2_ug_per_g

  PetscInt :: ires_nh3s, ires_nh3, ires_no3, ires_n2o, ires
  PetscInt :: ires_ngasnit

  PetscScalar, pointer :: bulkdensity(:)
  PetscReal :: c_nh3      ! mole/L
  PetscReal :: s_nh3      ! mole/m3
  PetscReal :: c_nh3_ugg  ! ug ammonium N / g soil
  PetscReal :: ph
  PetscReal :: rate_n2o, drate_n2o_dnh3
  PetscReal :: rate_nitri, drate_nitri_dnh3
  PetscReal :: f_t, f_w, f_ph
  PetscReal :: rho_b
  PetscReal :: theta
  PetscReal :: L_water
  PetscReal :: temp_real, feps0, dfeps0_dx

!---------------------------------------------------------------------------------

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  saturation = global_auxvar%sat(1)
  theta = saturation * porosity
  L_water = theta * volume * 1.0d3      ! Litres of H2O

  tc = global_auxvar%temp

  ! indices for C and N species
  ires_nh3 = this%ispec_nh3
  ires_no3 = this%ispec_no3
  ires_n2o = this%ispec_n2o
  ires_nh3s = this%ispec_nh4sorb + reaction%offset_immobile
  ires_ngasnit = this%ispec_ngasnit + reaction%offset_immobile

  c_nh3 = rt_auxvar%total(this%ispec_nh3, iphase) * L_water     ! moles/L --> moles
! (TODO) not sure if 'absorbed NH4' should be included in the reaction.
!  if (associated(rt_auxvar%total_sorb_eq)) then           ! original absorption-reactions in PF used
!     s_nh3 = rt_auxvar%total_sorb_eq(this%ispec_nh3)
!  elseif (this%ispec_nh4sorb>0) then                      ! 'reaction_sandbox_langmuir' used
!     s_nh3 = rt_auxvar%immobile(this%ispec_nh4sorb)
!  else
!     s_nh3 = 1.d-20
!  endif
! c_nh3 = c_nh3 + s_nh3*volume


  !if (c_nh3 < this%x0eps) return     ! this may bring in 'oscillation' around 'this%x0eps'
  if(this%x0eps>0.d0) then
    feps0 = c_nh3/(c_nh3+this%x0eps)      ! using these two for trailer smoothing, alternatively
    dfeps0_dx = this%x0eps/(c_nh3+this%x0eps)/(c_nh3+this%x0eps)
  else
    feps0 = 1.0d0
    dfeps0_dx = 0.d0
  endif

  ! nitrification (Dickinson et al. 2002)
  if(this%ispec_no3 > 0) then
    ! Eq. 28 for nitrification in Dickinson et al. 2002, can be separated into three parts:
    ! rate = k * f(trz) * (s*(1-s))/(0.25+1/nh4)
    ! (1) f_t = f(trz)
    ! (2) f_w = s*(1-s)/0.25  [ranging from 0 - 1]
    ! (3) f_nh4 = nh4/(nh4 + 4.0)

    f_t = exp(0.08d0 * (tc - 25.0d0))

    saturation = max(0.d0,min(saturation,1.d0))
    f_w = saturation * (1.0d0 - saturation)/0.25d0

    temp_real = this%k_nitr_max * f_t * f_w * c_nh3 / (c_nh3 + 4.0d0)                        ! 1/sec.
    rate_nitri = temp_real * (c_nh3 * feps0)                                                 ! moles/sec.

    Residual(ires_nh3) = Residual(ires_nh3) + rate_nitri
    Residual(ires_no3) = Residual(ires_no3) - rate_nitri

    if (compute_derivative) then

     ! d(rate_nh3))/d(nh3), in which, rate(c_nh3) = c_nh3*c_nh3/(c_nh3+4.0)*feps0
     temp_real = c_nh3*c_nh3/(c_nh3+4.0d0) * dfeps0_dx + &
                 c_nh3*(c_nh3+8.0d0)/(c_nh3+4.0d0)/(c_nh3+4.0d0) * feps0
     drate_nitri_dnh3 = f_t*f_w*this%k_nitr_max* temp_real
 
     Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) + drate_nitri_dnh3 * &
        rt_auxvar%aqueous%dtotal(this%ispec_nh3,this%ispec_nh3,iphase)

     Jacobian(ires_no3,ires_nh3) = Jacobian(ires_no3,ires_nh3) - drate_nitri_dnh3 * &
        rt_auxvar%aqueous%dtotal(this%ispec_no3,this%ispec_nh3,iphase)
    endif
  endif

! N2O production during nitrification (Parton et al. 1996)
#ifdef CLM_PFLOTRAN
  ghosted_id = option%iflag

  call VecGetArrayReadF90(clm_pf_idata%bulkdensity_dry_pf, bulkdensity, ierr)
  CHKERRQ(ierr)
  rho_b = bulkdensity(ghosted_id) ! kg/m3
  call VecRestoreArrayReadF90(clm_pf_idata%bulkdensity_dry_pf, bulkdensity, ierr)
  CHKERRQ(ierr)
#else
  rho_b = 1.25d3
#endif

! moleN_Lwater*g_mol*Lwater *1.d6 = ugN
  temp_real = N_molecular_weight*1.0d6      ! unit: ugN/(mol)
! g soil = m3_soil*kg_m3*1.d3                      ! unit: gSoil
  M_2_ug_per_g = temp_real/(volume*rho_b*1.d3)      ! unit: (ugN/gSoil)/(mol)
  !c_nh3_ugg = (c_nh3 + s_nh3 / theta / 1000.0d0)* M_2_ug_per_g
  c_nh3_ugg = c_nh3 * M_2_ug_per_g                 ! c_nh3 already in moles

  if(this%ispec_n2o > 0.0d0 .and. c_nh3_ugg > 3.0d0) then
    ! temperature response function (Parton et al. 1996)
    f_t = -0.06d0 + 0.13d0 * exp( 0.07d0 * tc )

    f_w = ((1.27d0 - saturation)/0.67d0)**(3.1777d0) * &
        ((saturation - 0.0012d0)/0.5988d0)**2.84d0

    ph = 6.5d0       ! default
    if (this%ispec_proton > 0) then
      if (reaction%species_idx%h_ion_id > 0) then
        ph = &
          -log10(rt_auxvar%pri_molal(reaction%species_idx%h_ion_id)* &
                 rt_auxvar%pri_act_coef(reaction%species_idx%h_ion_id))
      else if (reaction%species_idx%h_ion_id < 0) then
        ph = &
          -log10(rt_auxvar%sec_molal(abs(reaction%species_idx%h_ion_id))* &
                 rt_auxvar%sec_act_coef(abs(reaction%species_idx%h_ion_id)))
      endif
    endif
    f_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + ph))/rpi

    if(f_t > 0.0d0 .and. f_w > 0.0d0 .and. f_ph > 0.0d0) then
       f_t  = min(f_t, 1.0d0)
       f_w  = min(f_w, 1.0d0)
       f_ph = min(f_ph, 1.0d0)

       temp_real = (1.0 - exp(-0.0105d0 * c_nh3_ugg))  &
                  * f_t * f_w * f_ph * this%k_nitr_n2o                     ! 1/s
       rate_n2o = temp_real * (c_nh3*feps0)                                ! moles/s

       Residual(ires_nh3) = Residual(ires_nh3) + rate_n2o
       Residual(ires_n2o) = Residual(ires_n2o) - 0.5d0 * rate_n2o
       
       if(this%ispec_ngasnit > 0) then
          Residual(ires_ngasnit) = Residual(ires_ngasnit) - rate_n2o
       endif

       ! a note here: in case that want to track total nitrification rate,
       !              it should be = rate_nitri + rate_n2o

       if (compute_derivative) then

           ! d(rate_nh3))/d(nh3), in which,
           ! rate(nh3) = (c_nh3*feps0) * (1.0-exp(-0.0105d0*c_nh3*m_2_ug_per_g))
           temp_real = (c_nh3*dfeps0_dx + feps0) * &                   !d(c_nh3*feps0)/d(nh3)
                       (1.0 - exp(-0.0105d0*c_nh3_ugg))

           temp_real = temp_real + (c_nh3*feps0) *  &
                       0.0105d0*M_2_ug_per_g*exp(-0.0105d0*c_nh3_ugg)  !d(1.0-exp(-0.0105d0*c_nh3*m_2_ug_per_g))/d(nh3)

           drate_n2o_dnh3 = temp_real *this%k_nitr_n2o* f_t * f_w * f_ph
 
           Jacobian(ires_nh3,ires_nh3)=Jacobian(ires_nh3,ires_nh3) + &
             drate_n2o_dnh3 * &
             rt_auxvar%aqueous%dtotal(this%ispec_nh3,this%ispec_nh3,iphase)

           Jacobian(ires_n2o,ires_nh3)=Jacobian(ires_n2o,ires_nh3) - &
             0.5d0 * drate_n2o_dnh3 * &
             rt_auxvar%aqueous%dtotal(this%ispec_n2o,this%ispec_nh3,iphase)
      
       !    if(this%ispec_ngasnit > 0) then
       !      Jacobian(ires_ngasnit,ires_nh3)=Jacobian(ires_ngasnit,ires_nh3)- &
       !        drate_n2o_dnh3
       !    endif

       endif  !if (compute_derivative)

     endif    !if(f_t > 0.0d0 .and. f_w > 0.0d0 .and. f_ph > 0.0d0)

  endif       !if(this%ispec_n2o > 0.0d0 .and. c_nh3_ugg > 3.0d0)

#ifdef DEBUG
  do ires=1, reaction%ncomp
    temp_real = Residual(ires)

    if (abs(temp_real) > huge(temp_real)) then
      write(option%fid_out, *) 'infinity of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: NITRIFICATION'
      option%io_buffer = ' checking infinity of Residuals matrix @ NitrificationReact '
      call printErrMsg(option)
    endif
  enddo
#endif

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
