module Reaction_Sandbox_PlantNTake_class

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
    PetscReal :: x0eps
  contains
    procedure, public :: ReadInput => PlantNTakeRead
    procedure, public :: Setup => PlantNTakeSetup
    procedure, public :: Evaluate => PlantNTakeReact
    procedure, public :: Destroy => PlantNTakeDestroy
  end type reaction_sandbox_plantntake_type

  public :: PlantNTakeCreate

contains

! ************************************************************************** !
!
! PlantNTakeCreate: Allocates plantntake reaction object.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
function PlantNTakeCreate()

  implicit none
  
  class(reaction_sandbox_plantntake_type), pointer :: PlantNTakeCreate

  allocate(PlantNTakeCreate)
  PlantNTakeCreate%rate = 0.d0
  PlantNTakeCreate%half_saturation_nh3 = 1.d-15
  PlantNTakeCreate%half_saturation_no3 = 1.d-15
  PlantNTakeCreate%inhibition_nh3_no3  = 1.d-15
  PlantNTakeCreate%x0eps  = 1.d-20
  nullify(PlantNTakeCreate%next)  
      
end function PlantNTakeCreate

! ************************************************************************** !
!
! PlantNTakeRead: Reads input deck for plantntake reaction parameters
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNTakeRead(this,input,option)

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
      case('HALF_SATURATION_AMMONIA')
          call InputReadDouble(input,option,this%half_saturation_nh3)
          call InputErrorMsg(input,option,'half saturation for ammonia', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case('HALF_SATURATION_NITRATE')
          call InputReadDouble(input,option,this%half_saturation_no3)
          call InputErrorMsg(input,option,'half saturation for nitrate', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case('AMMONIA_INHIBITION_NITRATE')
          call InputReadDouble(input,option,this%inhibition_nh3_no3)
          call InputErrorMsg(input,option,'ammonia inhibition on nitrate', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case default
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
          call printErrMsg(option)
    end select
  enddo
  
end subroutine PlantNTakeRead

! ************************************************************************** !
!
! PlantNTakeSetup: Sets up the plantntake reaction with parameters either
!                read from the input deck or hardwired.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNTakeSetup(this,reaction,option)

  use Reaction_Aux_module
  use Option_module
  use Immobile_Aux_module

  implicit none
  
  class(reaction_sandbox_plantntake_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
 
end subroutine PlantNTakeSetup

! ************************************************************************** !
!
! PlantNTakeReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNTakeReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,porosity,volume,reaction, &
                         option,local_id)

  use Option_module
  use Reaction_Aux_module
  use Immobile_Aux_module

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
  PetscBool :: compute_derivative

  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: rate, drate, concN, rate0
  PetscReal :: volume, porosity
  PetscInt :: local_id
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word

  PetscInt, parameter :: iphase = 1
  PetscInt :: ispec_nh3, ispec_no3, ispec_plantn
  PetscInt :: ires_nh3, ires_no3, ires_plantn
  PetscInt :: ispec_nh4in, ispec_no3in, ires_nh4in, ires_no3in

  PetscReal :: c_nh3         ! concentration (mole/L)
  PetscReal :: f_nh3         ! nh3 / (half_saturation + nh3)
  PetscReal :: d_nh3         ! half_saturation / (half_saturation + nh3)^2
  PetscReal :: f_nh3_inhibit ! inhibition_coef/(inhibition_coef + nh3)
  PetscReal :: c_no3         ! concentration (mole/L)
  PetscReal :: f_no3         ! no3 / (half_saturation + no3)
  PetscReal :: d_no3         ! half_saturation/(no3 + half_saturation)^2 
  PetscReal :: temp_real

  PetscReal :: rate_nplant
  PetscReal :: rate_nplant_no3
  PetscReal :: drate_nplant       !drate/dnh4+ 
  PetscReal :: drate_nplant_no3   !drate/dno3-
  PetscReal :: rate_nh4, rate_no3

#ifdef CLM_PFLOTRAN
  PetscScalar, pointer :: rate_plantnuptake_pf_loc(:)   !
  PetscScalar, pointer :: rate_smin_nh4_pf_loc(:)   !
  PetscScalar, pointer :: rate_smin_no3_pf_loc(:)   !
#endif 

  word = 'AmmoniaH4+'
  ispec_nh3 = GetPrimarySpeciesIDFromName(word, reaction, PETSC_FALSE, option)

  if (ispec_nh3 < 0) then
      word = 'NH4+'
      ispec_nh3 = GetPrimarySpeciesIDFromName(word, reaction, PETSC_FALSE, option)
  endif

  if (ispec_nh3 < 0) then
      word = 'NH3(aq)'
      ispec_nh3 = GetPrimarySpeciesIDFromName(word, reaction, PETSC_FALSE, option)
  endif

  word = 'NO3-'
  ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction, PETSC_FALSE,option)

  word = 'PlantN'
  ispec_plantn = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  word = 'Ain'
  ispec_nh4in = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  word = 'Tin'
  ispec_no3in = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  ires_nh4in = ispec_nh4in + reaction%offset_immobile
  ires_no3in = ispec_no3in + reaction%offset_immobile

  if(ispec_plantn < 0) then
    write(*, *) 'Warning: PlantN is not specified in the chemical species!'
    return
  else
     ires_plantn = ispec_plantn + reaction%offset_immobile
  endif

  if (ispec_nh3 > 0) then
     c_nh3     = rt_auxvar%total(ispec_nh3, iphase)
     temp_real = c_nh3 - this%x0eps + this%half_saturation_nh3
     f_nh3     = (c_nh3 - this%x0eps) / temp_real
     d_nh3     = this%half_saturation_nh3 / temp_real / temp_real
     ires_nh3 = ispec_nh3
  endif

  if (ispec_no3 > 0) then
      c_no3 = rt_auxvar%total(ispec_no3, iphase)
      temp_real = c_no3 -this%x0eps + this%half_saturation_no3
      f_no3 = (c_no3 - this%x0eps) / temp_real
      d_no3 = this%half_saturation_no3 / temp_real / temp_real

      ires_no3 = ispec_no3

      if(ispec_nh3 > 0) then
        temp_real = this%inhibition_nh3_no3 + c_nh3
        f_nh3_inhibit = this%inhibition_nh3_no3/temp_real
      else
        f_nh3_inhibit = 1.0d0
      endif
  endif


#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%rate_plantnuptake_pf, &
       rate_plantnuptake_pf_loc, ierr)

  rate_nplant = rate_plantnuptake_pf_loc(local_id) * volume ! mol/m3/s * m3

  call VecRestoreArrayReadF90(clm_pf_idata%rate_plantnuptake_pf, &
       rate_plantnuptake_pf_loc, ierr)
#else
  rate_nplant = this%rate
#endif

  rate_nplant_no3 = rate_nplant * f_nh3_inhibit   ! NO3 uptake rate

  if(ispec_nh3 > 0) then
    drate_nplant = rate_nplant * d_nh3
    rate_nplant = rate_nplant * f_nh3

    Residual(ires_nh3) = Residual(ires_nh3) + rate_nplant
    Residual(ires_plantn) = Residual(ires_plantn) - rate_nplant

    if (compute_derivative) then
       Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3)+drate_nplant * &
         rt_auxvar%aqueous%dtotal(ispec_nh3,ispec_nh3,iphase)

       Jacobian(ires_plantn,ires_nh3)=Jacobian(ires_plantn,ires_nh3)-drate_nplant
    endif
  endif

  if(ispec_no3 > 0) then
    drate_nplant_no3 = rate_nplant_no3 * d_no3
    rate_nplant_no3  = rate_nplant_no3 * f_no3

    Residual(ires_no3) = Residual(ires_no3) + rate_nplant_no3
    Residual(ires_plantn) = Residual(ires_plantn) - rate_nplant_no3

    if (compute_derivative) then
     Jacobian(ires_no3,ires_no3)=Jacobian(ires_no3,ires_no3)+drate_nplant_no3* &
       rt_auxvar%aqueous%dtotal(ispec_no3,ispec_no3,iphase)

     Jacobian(ires_plantn,ires_no3)=Jacobian(ires_plantn,ires_no3)-drate_nplant_no3
    endif
  endif

#ifdef CLM_PFLOTRAN
  if(ispec_nh4in > 0) then
     call VecGetArrayReadF90(clm_pf_idata%rate_smin_nh4_pfs, &
       rate_smin_nh4_pf_loc, ierr)

     rate_nh4 = rate_smin_nh4_pf_loc(local_id) * volume ! mol/m3/s * m3
     Residual(ires_nh4in) = Residual(ires_nh4in) - rate_nh4

     call VecRestoreArrayReadF90(clm_pf_idata%rate_smin_nh4_pfs, &
       rate_smin_nh4_pf_loc, ierr)
  endif

  if(ispec_no3in > 0) then
     call VecGetArrayReadF90(clm_pf_idata%rate_smin_no3_pfs, &
       rate_smin_no3_pf_loc, ierr)

     rate_no3 = rate_smin_no3_pf_loc(local_id) * volume ! mol/m3/s * m3
     Residual(ires_no3in) = Residual(ires_no3in) - rate_no3

     call VecRestoreArrayReadF90(clm_pf_idata%rate_smin_no3_pfs, &
       rate_smin_no3_pf_loc, ierr)
  endif

#endif

end subroutine PlantNTakeReact

! ************************************************************************** !
!
! PlantNTakeDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNTakeDestroy(this)

  implicit none
  
  class(reaction_sandbox_plantntake_type) :: this  

end subroutine PlantNTakeDestroy

end module Reaction_Sandbox_PlantNTake_class
