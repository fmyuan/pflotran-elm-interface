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
  PetscReal :: rate, drate, concN, rate0
  PetscReal :: volume, porosity
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word

  PetscInt, parameter :: iphase = 1
  PetscInt :: ispec_nh3, ispec_no3, ispec_plantn
  PetscInt :: ires_nh3, ires_no3, ires_plantn
  PetscInt :: ispec_plantndemand, ires_plantndemand, ires

  PetscReal :: c_nh3         ! concentration (mole/L)
  PetscReal :: f_nh3         ! nh3 / (half_saturation + nh3)
  PetscReal :: d_nh3         ! half_saturation / (half_saturation + nh3)^2
  PetscReal :: c_no3         ! concentration (mole/L)
  PetscReal :: f_no3         ! no3 / (half_saturation + no3)
  PetscReal :: d_no3         ! half_saturation/(no3 + half_saturation)^2 
  PetscReal :: f_nh3_inhibit_no3 ! inhibition of nh3 on no3 uptake: inhibition_coef/(inhibition_coef + nh3)
  PetscReal :: d_nh3_inhibit_no3 ! derivative of inhibition of nh3 on no3 uptake (/nh3): -inhibition_coef/(inhibition_coef + nh3)^2
  PetscReal :: d_no3inhibit_nh3  ! derivative of f_no3*f_nh3_inhibit_no3 over dnh3
  PetscReal :: temp_real

  PetscReal :: rate_nplant_nh3
  PetscReal :: rate_nplant_no3
  PetscReal :: drate_nplant_nh3       !drate_nh3/dnh3
  PetscReal :: drate_nplant_no3       !drate_no3/dno3
  PetscReal :: drate_nplant_no3_nh3   !drate_no3/dnh3
  PetscReal :: rate_nh3, rate_no3
  PetscReal :: c_plantn, c_plantndemand

#ifdef CLM_PFLOTRAN
  PetscScalar, pointer :: rate_plantnuptake_pf_loc(:)   !
#endif 

!------------------------------------------------------------------------------------
  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  word = 'NH4+'
  ispec_nh3 = GetPrimarySpeciesIDFromName(word, reaction, PETSC_FALSE, option)

  if (ispec_nh3 < 0) then
      word = 'NH3(aq)'
      ispec_nh3 = GetPrimarySpeciesIDFromName(word, reaction, PETSC_FALSE, option)
  endif

  word = 'NO3-'
  ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction, PETSC_FALSE,option)

  word = 'PlantN'
  ispec_plantn = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  word = 'Plantndemand'
  ispec_plantndemand = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  ires_plantndemand = ispec_plantndemand + reaction%offset_immobile

  if(ispec_plantn < 0) then
    write(*, *) 'Warning: PlantN is not specified in the chemical species!'
    return
  else
     ires_plantn = ispec_plantn + reaction%offset_immobile
  endif

  if (ispec_nh3 > 0) then
     c_nh3     = rt_auxvar%total(ispec_nh3, iphase)
     temp_real = c_nh3 + this%half_saturation_nh3
     f_nh3     = c_nh3 / temp_real
     d_nh3     = this%half_saturation_nh3 / temp_real / temp_real

     ires_nh3 = ispec_nh3

    if (ispec_no3 > 0) then
      temp_real = this%inhibition_nh3_no3 + c_nh3
      f_nh3_inhibit_no3 = this%inhibition_nh3_no3/temp_real
      d_nh3_inhibit_no3 = -this%inhibition_nh3_no3/temp_real/temp_real
      !f_nh3_inhibit_no3 = 1.0d0 - f_nh3
      !d_nh3_inhibit_no3 = -d_nh3

      ! f_nh3 should be adjusted by (1-f_nh3_inhibit_no3) so that it not overlapping of rates
      d_nh3 = d_nh3 * (1.0d0-f_nh3_inhibit_no3) + f_nh3 * (-d_nh3_inhibit_no3)
      f_nh3 = f_nh3 * (1.0d0-f_nh3_inhibit_no3)

    endif

  endif

  if (ispec_no3 > 0) then
    c_no3 = rt_auxvar%total(ispec_no3, iphase)
    temp_real = c_no3 + this%half_saturation_no3
    f_no3 = c_no3 / temp_real
    d_no3 = this%half_saturation_no3 / temp_real / temp_real

    ires_no3 = ispec_no3

    if(ispec_nh3 > 0) then
      ! f_no3 should be adjusted by f_nh3_inhibit_no3 so that it is inhibited by nh3
      d_no3 = d_no3 * f_nh3_inhibit_no3     ! d(f_no3*f_nh3_inhibit_no3)/d(no3) (so do the derivative first, then 'fno3' is from above)
      d_no3inhibit_nh3 = f_no3*d_nh3_inhibit_no3   ! d(f_no3*f_nh3_inhibit_no3)/d(nh3)
      f_no3 = f_no3 * f_nh3_inhibit_no3

    endif

  endif

#ifdef CLM_PFLOTRAN
  ghosted_id = option%iflag

  call VecGetArrayReadF90(clm_pf_idata%rate_plantnuptake_pfs, &
       rate_plantnuptake_pf_loc, ierr)

  this%rate = rate_plantnuptake_pf_loc(ghosted_id) * volume ! mol/m3/s * m3

  call VecRestoreArrayReadF90(clm_pf_idata%rate_plantnuptake_pfs, &
       rate_plantnuptake_pf_loc, ierr)
#endif
  if (ispec_plantndemand > 0) then
    Residual(ires_plantndemand) = Residual(ires_plantndemand) - this%rate
  endif

  if(ispec_nh3 > 0 .and. c_nh3 > this%x0eps_nh4) then

    drate_nplant_nh3 = this%rate * d_nh3
    rate_nplant_nh3 = this%rate * f_nh3

    Residual(ires_nh3) = Residual(ires_nh3) + rate_nplant_nh3
    Residual(ires_plantn) = Residual(ires_plantn) - rate_nplant_nh3

    if (compute_derivative) then
       Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3)+drate_nplant_nh3 * &
         rt_auxvar%aqueous%dtotal(ispec_nh3,ispec_nh3,iphase)

       Jacobian(ires_plantn,ires_nh3)=Jacobian(ires_plantn,ires_nh3)-drate_nplant_nh3

    endif
  endif

  if(ispec_no3 > 0 .and. c_no3 > this%x0eps_no3) then
    rate_nplant_no3 = this%rate * f_no3
    drate_nplant_no3 = this%rate * d_no3
    drate_nplant_no3_nh3 = this%rate * d_no3inhibit_nh3

    Residual(ires_no3) = Residual(ires_no3) + rate_nplant_no3
    Residual(ires_plantn) = Residual(ires_plantn) - rate_nplant_no3

    if (compute_derivative) then
      Jacobian(ires_no3,ires_no3)=Jacobian(ires_no3,ires_no3)+drate_nplant_no3* &
        rt_auxvar%aqueous%dtotal(ispec_no3,ispec_no3,iphase)
      Jacobian(ires_plantn,ires_no3)=Jacobian(ires_plantn,ires_no3)-drate_nplant_no3

      !(TODO) need some thoughts/consulations to guoping/glenn for the following?
      Jacobian(ires_no3,ires_nh3)=Jacobian(ires_no3,ires_nh3)+drate_nplant_no3_nh3* &
        rt_auxvar%aqueous%dtotal(ispec_no3,ispec_nh3,iphase)
      Jacobian(ires_nh3,ires_no3)=Jacobian(ires_nh3,ires_no3)-drate_nplant_no3_nh3* &
        rt_auxvar%aqueous%dtotal(ispec_nh3,ispec_no3,iphase)

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
