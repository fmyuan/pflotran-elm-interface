module Reaction_Sandbox_PlantN_class

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_plantntake_type
    PetscReal :: rate
    PetscReal :: half_saturation_nh3
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh3_no3
    PetscReal :: x0eps_nh4
    PetscReal :: x0eps_no3

  contains
    procedure, public :: ReadInput => PlantNRead
    procedure, public :: Setup => PlantNSetup
    procedure, public :: Evaluate => PlantNReact
    procedure, public :: Destroy => PlantNDestroy
  end type reaction_sandbox_plantntake_type

  public :: PlantNCreate

contains

! ************************************************************************** !
!
! PlantNCreate: Allocates plantntake reaction object.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
function PlantNCreate()

  implicit none
  
  class(reaction_sandbox_plantntake_type), pointer :: PlantNCreate

  allocate(PlantNCreate)
  PlantNCreate%rate = 0.d0
  PlantNCreate%half_saturation_nh3 = 1.d-15
  PlantNCreate%half_saturation_no3 = 1.d-15
  PlantNCreate%inhibition_nh3_no3  = 1.d-15
  PlantNCreate%x0eps_nh4  = 1.d-20
  PlantNCreate%x0eps_no3  = 1.d-20
  nullify(PlantNCreate%next)
      
end function PlantNCreate

! ************************************************************************** !
!
! PlantNRead: Reads input deck for plantntake reaction parameters
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
  
  class(reaction_sandbox_plantntake_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE')
    call StringToUpper(word)   

    select case(trim(word))
      case('RATE')
          call InputReadDouble(input,option,this%rate)
          call InputErrorMsg(input,option,'rate', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case('AMMONIA_HALF_SATURATION')
          call InputReadDouble(input,option,this%half_saturation_nh3)
          call InputErrorMsg(input,option,'half saturation for ammonia', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case('NITRATE_HALF_SATURATION')
          call InputReadDouble(input,option,this%half_saturation_no3)
          call InputErrorMsg(input,option,'half saturation for nitrate', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case('AMMONIA_INHIBITION_NITRATE')
          call InputReadDouble(input,option,this%inhibition_nh3_no3)
          call InputErrorMsg(input,option,'ammonia inhibition on nitrate', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case('X0EPS_NH4')
          call InputReadDouble(input,option,this%x0eps_nh4)
          call InputErrorMsg(input,option,'x0eps_nh4', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case('X0EPS_NO3')
          call InputReadDouble(input,option,this%x0eps_no3)
          call InputErrorMsg(input,option,'x0eps_no3', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case default
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
          call printErrMsg(option)
    end select
  enddo
  
end subroutine PlantNRead

! ************************************************************************** !
!
! PlantNSetup: Sets up the plantntake reaction with parameters either
!                read from the input deck or hardwired.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNSetup(this,reaction,option)

  use Reaction_Aux_module
  use Option_module
  use Immobile_Aux_module

  implicit none
  
  class(reaction_sandbox_plantntake_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
 
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
  use Immobile_Aux_module
  use Material_Aux_class, only : material_auxvar_type

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif
  
  implicit none

#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif
  
  class(reaction_sandbox_plantntake_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscBool :: compute_derivative

  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: rate, drate_dn, nconc
  PetscReal :: volume, porosity
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word

  PetscInt, parameter :: iphase = 1
  PetscInt :: ispec_nh3, ispec_no3, ispec_plantn
  PetscInt :: ires_nh3, ires_no3, ires_plantn
  PetscInt :: ispec_plantndemand, ires_plantndemand
  PetscInt :: ispec_plantnuptake, ires_plantnuptake
  PetscInt :: ires

  PetscReal :: c_nh3         ! concentration (mole/L)
  PetscReal :: fnh3          ! nh3 / (half_saturation + nh3): rate dependence on substrate
  PetscReal :: dfnh3_dnh3    ! d(fnh3)/d(nh3)

  PetscReal :: c_no3         ! concentration (mole/L)
  PetscReal :: fno3          ! no3 / (half_saturation + no3): rate dependence on substrate
  PetscReal :: dfno3_dno3    ! d(fno3)/d(no3)
  PetscReal :: fnh3_inhibit_no3        ! inhibition of nh3 on no3 uptake (plant preference over nh3/no3): inhibition_coef/(inhibition_coef + nh3)
  PetscReal :: dfnh3_inhibit_no3_dnh3  ! derivative of d(fnh3_inhibit_no3)/d(nh3)
  PetscReal :: dfnh3_inhibit_no3_dno3  ! derivative of d(fnh3_inhibit_no3)/d(no3)

  PetscReal :: temp_real, feps0, dfeps0_dx

  PetscReal :: nrate_nh3
  PetscReal :: nrate_no3
  PetscReal :: dnrate_nh3_dnh3       !d(nrate_nh3)/d(nh3)
  PetscReal :: dnrate_nh3_dno3       !d(nrate_nh3)/d(no3)
  PetscReal :: dnrate_no3_dnh3       !d(nrate_no3)/d(nh3)
  PetscReal :: dnrate_no3_dno3       !d(nrate_no3)/d(no3)

#ifdef CLM_PFLOTRAN
  PetscScalar, pointer :: rate_plantndemand_pf_loc(:)   !
#endif 

!------------------------------------------------------------------------------------
  word = 'NH4+'
  ispec_nh3 = GetPrimarySpeciesIDFromName(word, reaction, PETSC_FALSE, option)

  word = 'NO3-'
  ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction, PETSC_FALSE,option)

  if(ispec_nh3 < 0 .and. ispec_no3 < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE: ' // &
       ' at least one of NH4+ and NO3- must be specified as primary species in the input file.'
     call printErrMsg(option)
  endif

  word = 'PlantN'
  ispec_plantn = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  if(ispec_plantn < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE: ' // &
       ' PlantN is not specified as immobile species in the input file.'
     call printErrMsg(option)
  else
     ires_plantn = ispec_plantn + reaction%offset_immobile
  endif

  word = 'Plantndemand'
  ispec_plantndemand = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)
  ires_plantndemand = ispec_plantndemand + reaction%offset_immobile

  word = 'Plantnuptake'
  ispec_plantnuptake = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)
  ires_plantnuptake = ispec_plantnuptake + reaction%offset_immobile
#ifdef CLM_PFLOTRAN
  if(ispec_plantnuptake < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE: ' // &
       'Plantnuptake is not specified as immobile species in the ' // &
       'input file, but It is required when coupled with CLM.'
     call printErrMsg(option)
  endif
#endif

  !--------------------------------------------------------------------------------------------
  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  if (ispec_nh3 > 0) then
    ires_nh3 = ispec_nh3

    c_nh3     = rt_auxvar%total(ispec_nh3, iphase)
    temp_real = c_nh3 + this%half_saturation_nh3
    fnh3      = c_nh3 / temp_real
    dfnh3_dnh3= this%half_saturation_nh3 / temp_real / temp_real

    if (ispec_no3 > 0) then
      temp_real = this%inhibition_nh3_no3 + c_nh3
      fnh3_inhibit_no3 = this%inhibition_nh3_no3/temp_real
      dfnh3_inhibit_no3_dnh3 = -this%inhibition_nh3_no3/temp_real/temp_real
      dfnh3_inhibit_no3_dno3 = 0.d0

    endif

  endif

  if (ispec_no3 > 0) then
    ires_no3 = ispec_no3

    c_no3     = rt_auxvar%total(ispec_no3, iphase)
    temp_real = c_no3 + this%half_saturation_no3
    fno3      = c_no3 / temp_real
    dfno3_dno3= this%half_saturation_no3 / temp_real / temp_real

  endif

#ifdef CLM_PFLOTRAN
  ghosted_id = option%iflag

  call VecGetArrayReadF90(clm_pf_idata%rate_plantndemand_pfs, &
       rate_plantndemand_pf_loc, ierr)

  this%rate = rate_plantndemand_pf_loc(ghosted_id) * volume          ! moles/m3/s * m3

  call VecRestoreArrayReadF90(clm_pf_idata%rate_plantndemand_pfs, &
       rate_plantndemand_pf_loc, ierr)
#endif
  if (ispec_plantndemand > 0) then  ! for tracking
    Residual(ires_plantndemand) = Residual(ires_plantndemand) - this%rate
  endif

  if(ispec_nh3 > 0) then

    if(this%x0eps_nh4>0.d0) then
      feps0  = c_nh3/(c_nh3 + this%x0eps_nh4)         ! for trailer smoothing
      dfeps0_dx = this%x0eps_nh4/(c_nh3 + this%x0eps_nh4)/(c_nh3 + this%x0eps_nh4)
    else
      feps0 = 1.d0
      dfeps0_dx = 0.d0
    endif

    ! rates
    nrate_nh3 = this%rate * fnh3 * feps0
    if(ispec_no3 > 0) then
    ! splitting (fractioning) potential uptake rate by the '1-fnh3_inhibition_no3'
      nrate_nh3 = this%rate * fnh3 * feps0 * &
                  (1.0-fnh3_inhibit_no3)
    endif

    ! residuals
    Residual(ires_nh3) = Residual(ires_nh3) + nrate_nh3
    Residual(ires_plantn) = Residual(ires_plantn) - nrate_nh3

    if (ispec_plantnuptake>0) then   ! for tracking
      Residual(ires_plantnuptake) = Residual(ires_plantnuptake) - nrate_nh3
    endif

    ! jacobians
    if(compute_derivative) then

      dnrate_nh3_dnh3 = this%rate * (dfnh3_dnh3*feps0 + fnh3*dfeps0_dx)
      if(ispec_no3 > 0) then
        temp_real = fnh3 * feps0 * (-1.0d0*dfnh3_inhibit_no3_dnh3) + &
                    (1.0-fnh3_inhibit_no3) * &
                    (fnh3*dfeps0_dx + dfnh3_dnh3*feps0)
        dnrate_nh3_dnh3 = this%rate * temp_real

        dnrate_nh3_dno3 = 0.d0    ! no rate dependence on 'no3' currently
      endif

      Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) + &
        dnrate_nh3_dnh3 * &
        rt_auxvar%aqueous%dtotal(ispec_nh3,ispec_nh3,iphase)

      Jacobian(ires_plantn,ires_nh3) = Jacobian(ires_plantn,ires_nh3) - &
        dnrate_nh3_dnh3

      if(ispec_no3 > 0) then
        Jacobian(ires_nh3,ires_no3)=Jacobian(ires_nh3,ires_no3) + &       ! may need a checking of the sign (+/-) here
          dnrate_nh3_dno3 * &
          rt_auxvar%aqueous%dtotal(ispec_nh3,ispec_no3,iphase)
      endif

    endif ! if(compute_derivative)

  endif !if(ispec_nh3 > 0)

  if(ispec_no3 > 0) then

    if(this%x0eps_no3>0.d0) then
      feps0  = c_no3/(c_no3 + this%x0eps_no3)         ! for trailer smoothing
      dfeps0_dx = this%x0eps_no3 &
                 /(c_no3 + this%x0eps_no3)/(c_no3 + this%x0eps_no3)
    else
      feps0 = 1.d0
      dfeps0_dx = 0.d0
    endif

    ! rates
    nrate_no3 = this%rate * fno3 * feps0
    if(ispec_nh3 > 0) then
    ! splitting (fractioning) potential uptake rate by the 'fnh3_inhibition_no3'
      nrate_no3 = this%rate * fno3 * feps0 * &
                  fnh3_inhibit_no3
    endif

    ! residuals
    Residual(ires_no3) = Residual(ires_no3) + nrate_no3
    Residual(ires_plantn) = Residual(ires_plantn) - nrate_no3

    if (compute_derivative) then

      dnrate_no3_dno3 = this%rate * (dfno3_dno3*feps0 + fno3*dfeps0_dx)
      if(ispec_nh3 > 0) then
        temp_real = fno3 * feps0 *dfnh3_inhibit_no3_dno3 + &             ! 'dfnh3_inhibit_no3_dno3'=0
                    fnh3_inhibit_no3 * &
                    (fno3*dfeps0_dx + dfno3_dno3*feps0)
        dnrate_no3_dno3 = this%rate * temp_real

        dnrate_no3_dnh3 = this%rate*fno3*feps0*dfnh3_inhibit_no3_dnh3    ! both 'fno3' and 'feps0' independent of 'nh3'
      endif

      Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) + &
        dnrate_no3_dno3 * &
        rt_auxvar%aqueous%dtotal(ispec_no3,ispec_no3,iphase)

      Jacobian(ires_plantn,ires_no3) = Jacobian(ires_plantn,ires_no3) - &
        dnrate_no3_dno3

      if(ispec_nh3 > 0) then
        Jacobian(ires_no3,ires_nh3)=Jacobian(ires_no3,ires_nh3) + &      ! may need a checking of sign (+/-) here
          dnrate_no3_dnh3 * &
          rt_auxvar%aqueous%dtotal(ispec_no3,ispec_nh3,iphase)
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
  
  class(reaction_sandbox_plantntake_type) :: this  

end subroutine PlantNDestroy

end module Reaction_Sandbox_PlantN_class
