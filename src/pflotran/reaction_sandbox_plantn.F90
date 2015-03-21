module Reaction_Sandbox_PlantN_class

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_plantn_type
    PetscInt  :: ispec_nh4
    PetscInt  :: ispec_no3
    PetscInt  :: ispec_plantn
    PetscInt  :: ispec_plantndemand
    PetscInt  :: ispec_plantnuptake
    PetscReal :: rate_plantndemand
    PetscReal :: half_saturation_nh4
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh4_no3
    PetscReal :: x0eps_nh4
    PetscReal :: x0eps_no3

  contains
    procedure, public :: ReadInput => PlantNRead
    procedure, public :: Setup => PlantNSetup
    procedure, public :: Evaluate => PlantNReact
    procedure, public :: Destroy => PlantNDestroy
  end type reaction_sandbox_plantn_type

  public :: PlantNCreate

contains

! ************************************************************************** !
!
! PlantNCreate: Allocates plantn reaction object.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
function PlantNCreate()

  implicit none
  
  class(reaction_sandbox_plantn_type), pointer :: PlantNCreate

  allocate(PlantNCreate)
  PlantNCreate%ispec_nh4 = 0
  PlantNCreate%ispec_no3 = 0
  PlantNCreate%ispec_plantn = 0
  PlantNCreate%ispec_plantndemand = 0
  PlantNCreate%ispec_plantnuptake = 0
  PlantNCreate%rate_plantndemand = 0.d0
  PlantNCreate%half_saturation_nh4 = 1.d-15
  PlantNCreate%half_saturation_no3 = 1.d-15
  PlantNCreate%inhibition_nh4_no3  = 1.d0
  PlantNCreate%x0eps_nh4  = 1.d-20
  PlantNCreate%x0eps_no3  = 1.d-20
  nullify(PlantNCreate%next)
      
end function PlantNCreate

! ************************************************************************** !
!
! PlantNRead: Reads input deck for plantn reaction parameters
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_plantn_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,PLANTN')
    call StringToUpper(word)   

    select case(trim(word))
      case('RATE_PLANTNDEMAND')
          call InputReadDouble(input,option,this%rate_plantndemand)
          call InputErrorMsg(input,option,'rate_plantndemand', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('AMMONIUM_HALF_SATURATION')
          call InputReadDouble(input,option,this%half_saturation_nh4)
          call InputErrorMsg(input,option,'half saturation for ammonium', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('NITRATE_HALF_SATURATION')
          call InputReadDouble(input,option,this%half_saturation_no3)
          call InputErrorMsg(input,option,'half saturation for nitrate', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('AMMONIUM_INHIBITION_NITRATE')
          call InputReadDouble(input,option,this%inhibition_nh4_no3)
          call InputErrorMsg(input,option,'ammonium inhibition on nitrate', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
          if (this%inhibition_nh4_no3<0.d-20) then
            option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN,' // &
              'AMMONIUM_INHIBITION_NITRATE cannot be too small to close to 0'
            call printErrMsg(option)
          endif
      case('X0EPS_NH4')
          call InputReadDouble(input,option,this%x0eps_nh4)
          call InputErrorMsg(input,option,'x0eps_nh4', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('X0EPS_NO3')
          call InputReadDouble(input,option,this%x0eps_no3)
          call InputErrorMsg(input,option,'x0eps_no3', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case default
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
          call printErrMsg(option)
    end select
  enddo
  
end subroutine PlantNRead

! ************************************************************************** !
!
! PlantNSetup: Sets up the plantn reaction with parameters either
!                read from the input deck or hardwired.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNSetup(this,reaction,option)

  use Reaction_Aux_module
  use Option_module
  use Reaction_Immobile_Aux_module

  implicit none
  
  class(reaction_sandbox_plantn_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

!------------------------------------------------------------------------------------
  word = 'NH4+'
  this%ispec_nh4 = GetPrimarySpeciesIDFromName(word, reaction, PETSC_FALSE, option)

  word = 'NO3-'
  this%ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction, PETSC_FALSE,option)

  if(this%ispec_nh4 < 0 .and. this%ispec_no3 < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN: ' // &
       ' at least one of NH4+ and NO3- must be specified as primary species in the input file.'
     call printErrMsg(option)
  endif

  word = 'PlantN'
  this%ispec_plantn = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  if(this%ispec_plantn < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN: ' // &
       ' PlantN is not specified as immobile species in the input file.'
     call printErrMsg(option)
  endif

  word = 'Plantndemand'
  this%ispec_plantndemand = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  word = 'Plantnuptake'
  this%ispec_plantnuptake = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)
#ifdef CLM_PFLOTRAN
  if(this%ispec_plantnuptake < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN: ' // &
       'Plantnuptake is not specified as immobile species in the ' // &
       'input file, but It is required when coupled with CLM.'
     call printErrMsg(option)
  endif
#endif

end subroutine PlantNSetup

! ************************************************************************** !
!
! PlantNReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 09/09/2013
!
! Rewritten by Fengming Yuan @Aug-14-2014. The orginal was messed-up with 'patches',
! which caused a lot of issues.
!
! ************************************************************************** !
subroutine PlantNReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)

  use Option_module
  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module
  use Material_Aux_class, only : material_auxvar_type

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif
  
  implicit none

#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif
  
  class(reaction_sandbox_plantn_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscBool :: compute_derivative

  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: rate_plantndemand, drate_dn, nconc
  PetscReal :: volume, porosity, saturation, tc
  PetscReal :: theta, L_water
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word

  PetscInt, parameter :: iphase = 1
  PetscInt :: ires_nh4, ires_no3, ires_plantn
  PetscInt :: ires_plantndemand, ires_plantnuptake
  PetscInt :: ires

  PetscReal :: c_nh4         ! concentration (mole/L)
  PetscReal :: fnh4          ! nh4 / (half_saturation + nh4): rate dependence on substrate
  PetscReal :: dfnh4_dnh4    ! d(fnh4)/d(nh4)

  PetscReal :: c_no3         ! concentration (mole/L)
  PetscReal :: fno3          ! no3 / (half_saturation + no3): rate dependence on substrate
  PetscReal :: dfno3_dno3    ! d(fno3)/d(no3)

  ! nh4 inhibition on no3 uptake, or plant N uptake preference btw nh4 and no3
  ! (Currently it's similar function as microbial N immobilization)
  ! crate_nh4 = fnh4*fnh4_inhibit_no3, while crate_no3 = 1.-fnh4*fnh4_inhibition_no3
  ! by the following eq., if 'inhibition=1', uptake will be equal for both (if same conc.);
  ! higher coef, strong NH4 inhibition on NO3 (i.e., more NH4 uptake over NO3)
  PetscReal :: fnh4_inhibit_no3 ! inhibition_coef/(inhibition_coef + no3/nh4):
  PetscReal :: dfnh4_inhibit_no3_dnh4 ! d(fnh4_inhibit_no3)/dnh4
  PetscReal :: dfnh4_inhibit_no3_dno3 ! d(fnh4_inhibit_no3)/dno3

  PetscReal :: temp_real, feps0, dfeps0_dx

  PetscReal :: nrate_nh4
  PetscReal :: nrate_no3
  PetscReal :: dnrate_nh4_dnh4       !d(nrate_nh4)/d(nh4)
  PetscReal :: dnrate_nh4_dno3       !d(nrate_nh4)/d(no3)
  PetscReal :: dnrate_no3_dnh4       !d(nrate_no3)/d(nh4)
  PetscReal :: dnrate_no3_dno3       !d(nrate_no3)/d(no3)

#ifdef CLM_PFLOTRAN
  PetscScalar, pointer :: rate_plantndemand_pf_loc(:)   !
#endif 

  !--------------------------------------------------------------------------------------------
  volume = material_auxvar%volume
  porosity = material_auxvar%porosity
  saturation = global_auxvar%sat(1)
  if(saturation < 0.01d0) return

  theta = saturation * porosity
  L_water = theta * 1.0d3      ! Litres of H2O /m3 bulk soil

  tc = global_auxvar%temp
  if(tc < 0.01d0) return

  ires_plantn = this%ispec_plantn + reaction%offset_immobile

  !--------------------------------------------------------------------------------------------
  if (this%ispec_nh4 > 0) then
    ires_nh4 = this%ispec_nh4 + reaction%offset_aqueous

    c_nh4     = rt_auxvar%total(this%ispec_nh4, iphase) * L_water
    temp_real = c_nh4 + this%half_saturation_nh4
    fnh4      = c_nh4 / temp_real
    dfnh4_dnh4= this%half_saturation_nh4 / temp_real / temp_real

    ! the following may not be needed, but just in case
    if(this%x0eps_nh4>0.d0) then
      !feps0  = c_nh4/(c_nh4 + this%x0eps_nh4)         ! for trailer smoothing
      !dfeps0_dx = this%x0eps_nh4/(c_nh4 + this%x0eps_nh4)/(c_nh4 + this%x0eps_nh4)

      ! GP's cut-off approach (from 'x0eps*10' to 'x0eps')
      if (c_nh4 <= this%x0eps_nh4) then
        feps0     = 0.0d0
        dfeps0_dx = 0.0d0
      elseif (c_nh4 >= this%x0eps_nh4*1.d1) then
        feps0     = 1.0d0
        dfeps0_dx = 0.0d0
      else
        feps0 = 1.0d0 - ( 1.0d0-(c_nh4-this%x0eps_nh4)*(c_nh4-this%x0eps_nh4)       &
                                /(81.0d0*this%x0eps_nh4*this%x0eps_nh4) ) ** 2
        dfeps0_dx = 4.0d0 * (1.0d0 - (c_nh4-this%x0eps_nh4)*(c_nh4-this%x0eps_nh4)  &
                                     /(81.0d0*this%x0eps_nh4*this%x0eps_nh4) )      &
                   * (c_nh4-this%x0eps_nh4)/(81.0d0*this%x0eps_nh4*this%x0eps_nh4)
      endif

      dfnh4_dnh4 = dfnh4_dnh4*feps0 + dfeps0_dx*fnh4  ! do the derivative first
      fnh4 = fnh4 * feps0
    endif

    !
    ! nh4 inhibition on no3 uptake, if any ('this%inhibition_nh4_no3')
    ! this is for quantifying plant N uptake preference btw NH4 and NO3
    ! (very similar as N immobilization now).
    if (this%ispec_no3 > 0) then
      c_no3     = rt_auxvar%total(this%ispec_no3, iphase) * L_water
       ! (DON'T change the 'rate' and 'derivatives' after this)
      if(c_nh4>this%x0eps_nh4 .and. c_no3>this%x0eps_no3 &
        .and. this%inhibition_nh4_no3>0.d0) then
        temp_real = this%inhibition_nh4_no3 + c_no3/c_nh4
        fnh4_inhibit_no3 = this%inhibition_nh4_no3/temp_real
        dfnh4_inhibit_no3_dnh4 = this%inhibition_nh4_no3*(c_no3/c_nh4/c_nh4)  &
                               /temp_real/temp_real     ! over 'd_nh4'
        dfnh4_inhibit_no3_dno3 = -this%inhibition_nh4_no3/c_nh4  &
                               /temp_real/temp_real     ! over 'd_no3'
      else
        if (c_nh4>this%x0eps_nh4 .and. c_no3<=this%x0eps_no3) then
          fnh4_inhibit_no3 = 0.999d0
        elseif (c_nh4<=this%x0eps_nh4 .and. c_no3>this%x0eps_no3) then
          fnh4_inhibit_no3 = 0.001d0
        else
          fnh4_inhibit_no3 = 0.50d0
        endif
        dfnh4_inhibit_no3_dnh4 = 0.d0
        dfnh4_inhibit_no3_dno3 = 0.d0
      endif

    endif

  endif

  if (this%ispec_no3 > 0) then
    ires_no3 = this%ispec_no3

    c_no3     = rt_auxvar%total(this%ispec_no3, iphase)
    temp_real = c_no3 + this%half_saturation_no3
    fno3      = c_no3 / temp_real
    dfno3_dno3= this%half_saturation_no3 / temp_real / temp_real

    ! the following may not be needed, but just in case
    if(this%x0eps_no3>0.d0) then
      !feps0  = c_no3/(c_no3 + this%x0eps_no3)         ! for trailer smoothing
      !dfeps0_dx = this%x0eps_no3 &
      !           /(c_no3 + this%x0eps_no3)/(c_no3 + this%x0eps_no3)

      ! GP's cut-off approach (from 'x0eps*10' to 'x0eps')
      if (c_no3 <= this%x0eps_no3) then
        feps0     = 0.0d0
        dfeps0_dx = 0.0d0
      elseif (c_no3 >= this%x0eps_no3*1.d1) then
        feps0     = 1.0d0
        dfeps0_dx = 0.0d0
      else
        feps0 = 1.0d0 - ( 1.0d0-(c_no3-this%x0eps_no3)*(c_no3-this%x0eps_no3)       &
                                /(81.0d0*this%x0eps_no3*this%x0eps_no3) ) ** 2
        dfeps0_dx = 4.0d0 * (1.0d0 - (c_no3-this%x0eps_no3)*(c_no3-this%x0eps_no3)  &
                                     /(81.0d0*this%x0eps_no3*this%x0eps_no3) )      &
                   * (c_no3-this%x0eps_no3)/(81.0d0*this%x0eps_no3*this%x0eps_no3)
      endif

      dfno3_dno3 = dfno3_dno3*feps0 + dfeps0_dx*fno3  ! do the derivative first
      fno3 = fno3 * feps0
    endif


  endif

#ifdef CLM_PFLOTRAN
  ghosted_id = option%iflag

  call VecGetArrayReadF90(clm_pf_idata%rate_plantndemand_pfs, &
       rate_plantndemand_pf_loc, ierr)

  this%rate_plantndemand = rate_plantndemand_pf_loc(ghosted_id) * volume          ! moles/m3/s * m3

  call VecRestoreArrayReadF90(clm_pf_idata%rate_plantndemand_pfs, &
       rate_plantndemand_pf_loc, ierr)
#endif
  if (this%ispec_plantndemand > 0) then  ! for tracking
    ires_plantndemand = this%ispec_plantndemand + reaction%offset_immobile
    Residual(ires_plantndemand) = Residual(ires_plantndemand) - this%rate_plantndemand
  endif

  if(this%ispec_nh4 > 0) then

    ! rates
    nrate_nh4 = this%rate_plantndemand * fnh4
    if(this%ispec_no3 > 0) then
    ! splitting (fractioning) potential uptake rate by the 'fnh4_inhibition_no3' for NH4 uptake
      nrate_nh4 = this%rate_plantndemand * fnh4 * fnh4_inhibit_no3
    endif

    ! residuals
    Residual(ires_nh4) = Residual(ires_nh4) + nrate_nh4
    Residual(ires_plantn) = Residual(ires_plantn) - nrate_nh4

    if (this%ispec_plantnuptake>0) then   ! for tracking
      ires_plantnuptake = this%ispec_plantnuptake + reaction%offset_immobile
      Residual(ires_plantnuptake) = Residual(ires_plantnuptake) - nrate_nh4
    endif

    ! jacobians
    if(compute_derivative) then

      dnrate_nh4_dnh4 = this%rate_plantndemand * dfnh4_dnh4
      if(this%ispec_no3 > 0) then
        temp_real = fnh4 * dfnh4_inhibit_no3_dnh4 + &
                    fnh4_inhibit_no3 * dfnh4_dnh4
        dnrate_nh4_dnh4 = this%rate_plantndemand * temp_real

        temp_real = fnh4 * dfnh4_inhibit_no3_dno3     ! dfnh4_dno3 = 0
        dnrate_nh4_dno3 = this%rate_plantndemand * temp_real
      endif

      Jacobian(ires_nh4,ires_nh4) = Jacobian(ires_nh4,ires_nh4) + &
        dnrate_nh4_dnh4 * &
        rt_auxvar%aqueous%dtotal(this%ispec_nh4,this%ispec_nh4,iphase)

      Jacobian(ires_plantn,ires_nh4) = Jacobian(ires_plantn,ires_nh4) - &
        dnrate_nh4_dnh4

      if(this%ispec_no3 > 0) then
        Jacobian(ires_nh4,ires_no3)=Jacobian(ires_nh4,ires_no3) + &       ! may need a checking of the sign (+/-) here
          dnrate_nh4_dno3 * &
          rt_auxvar%aqueous%dtotal(this%ispec_nh4,this%ispec_no3,iphase)
      endif

    endif ! if(compute_derivative)

  endif !if(this%ispec_nh4 > 0)

  if(this%ispec_no3 > 0) then

    ! rates
    nrate_no3 = this%rate_plantndemand * fno3
    if(this%ispec_nh4 > 0) then
    ! splitting (fractioning) potential uptake rate by the rest of nrate_nh4,
    ! which adjusted by 'fnh4_inhibition_no3'
    ! i.e., 1.0-fnh4*fnh4_inhibit_no3
      nrate_no3 = this%rate_plantndemand * fno3 * &
                  (1.0d0-fnh4*fnh4_inhibit_no3)
    endif

    ! residuals
    Residual(ires_no3) = Residual(ires_no3) + nrate_no3
    Residual(ires_plantn) = Residual(ires_plantn) - nrate_no3

    if (compute_derivative) then

      dnrate_no3_dno3 = this%rate_plantndemand * dfno3_dno3
      if(this%ispec_nh4 > 0) then
        temp_real = dfno3_dno3 * (1.d0-fnh4*fnh4_inhibit_no3) + &
                    fno3 * (-1.0d0*fnh4*dfnh4_inhibit_no3_dno3)              ! 'dfnh4_dno3=0'
        dnrate_no3_dno3 = this%rate_plantndemand * temp_real

        temp_real = fno3 * (-1.0d0) * &                                      ! 'dfno3_dnh4=0'
                    ( fnh4*dfnh4_inhibit_no3_dnh4 + &
                      dfnh4_dnh4*fnh4_inhibit_no3 )
        dnrate_no3_dnh4 = this%rate_plantndemand * temp_real
      endif

      Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) + &
        dnrate_no3_dno3 * &
        rt_auxvar%aqueous%dtotal(this%ispec_no3,this%ispec_no3,iphase)

      Jacobian(ires_plantn,ires_no3) = Jacobian(ires_plantn,ires_no3) - &
        dnrate_no3_dno3

      if(this%ispec_nh4 > 0) then
        Jacobian(ires_no3,ires_nh4)=Jacobian(ires_no3,ires_nh4) + &      ! may need a checking of sign (+/-) here
          dnrate_no3_dnh4 * &
          rt_auxvar%aqueous%dtotal(this%ispec_no3,this%ispec_nh4,iphase)
      endif

    endif
  endif

#ifdef DEBUG
  do ires=1, reaction%ncomp
    temp_real = Residual(ires)

    if (abs(temp_real) > huge(temp_real)) then
      write(option%fid_out, *) 'infinity of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: PLANT N UPTAKE'
      option%io_buffer = ' checking infinity of Residuals matrix @ PlantNReact '
      call printErrMsg(option)
    endif
  enddo
#endif

end subroutine PlantNReact

! ************************************************************************** !
!
! PlantNDestroy: Destroys allocatable or pointer objects created in this
!                  module
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNDestroy(this)

  implicit none
  
  class(reaction_sandbox_plantn_type) :: this

end subroutine PlantNDestroy

end module Reaction_Sandbox_PlantN_class
