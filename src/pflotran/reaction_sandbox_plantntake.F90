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
    PetscInt :: ispec_mineralN
    PetscInt :: ispec_plantN
    PetscReal :: rate
    PetscReal :: half_saturation
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
! author: Guoping Tang (replace in all subroutine headers with name of developer) 
! date: 09/09/2013 (replace in all subroutine headers with current date)
!
! ************************************************************************** !
function PlantNTakeCreate()

  implicit none
  
  class(reaction_sandbox_plantntake_type), pointer :: PlantNTakeCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(PlantNTakeCreate)
  PlantNTakeCreate%ispec_mineralN = 0
  PlantNTakeCreate%ispec_plantN = 0
  PlantNTakeCreate%rate = 0.d0
  PlantNTakeCreate%half_saturation = -10.d0
  nullify(PlantNTakeCreate%next)  
      
end function PlantNTakeCreate

! ************************************************************************** !
!
! PlantNTakeRead: Reads input deck for plantntake reaction parameters (if any)
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
           case('N_INHIBITION')
              call InputReadDouble(input,option,this%half_saturation)
              call InputErrorMsg(input,option,'inhibition coefficient', &
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
! PlantNTakeSetup: Sets up the plantntake reaction either with parameters either
!                read from the input deck or hardwired.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNTakeSetup(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module
  use Immobile_Aux_module, only : GetImmobileSpeciesIDFromName 

  implicit none
  
  class(reaction_sandbox_plantntake_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
 
  word = 'N'
  this%ispec_mineralN = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  word = 'PlantN'
  this%ispec_plantN = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)
      
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

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data !, only : rate_plantnuptake_pf 
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

  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: rate, drate, concN, rate0
  PetscReal :: volume, porosity
  PetscInt :: local_id
  PetscInt :: ires_mineralN, ires_plantN
  PetscErrorCode :: ierr

#ifdef CLM_PFLOTRAN
  PetscScalar, pointer :: rate_plantnuptake_pf_loc(:)   !
  call VecGetArrayReadF90(clm_pf_idata%rate_plantnuptake_pf, rate_plantnuptake_pf_loc, ierr)

  rate0 = rate_plantnuptake_pf_loc(local_id) * volume ! mol/m3/s * m3

  call VecRestoreArrayReadF90(clm_pf_idata%rate_plantnuptake_pf, rate_plantnuptake_pf_loc, ierr)
#else
  rate0 = this%rate
#endif

  if(this%half_saturation .GT. 1.0d-20) then
     concN = rt_auxvar%immobile(this%ispec_mineralN)
     rate = rate0 * concN/(concN + this%half_saturation) 
  endif

  ires_mineralN = this%ispec_mineralN + reaction%offset_immobile      
  ires_plantN = this%ispec_plantN + reaction%offset_immobile      

  Residual(ires_mineralN) = Residual(ires_mineralN) - (-1.0) * rate
  Residual(ires_plantN) = Residual(ires_plantN) - rate

  if (compute_derivative) return

  if(this%half_saturation .LT. 1.0d-20) return

  drate = rate0 * this%half_saturation / (concN + this%half_saturation) & 
                                      / (concN + this%half_saturation) 

    ! always add contribution to Jacobian
  Jacobian(ires_mineralN,ires_mineralN) = Jacobian(ires_mineralN,ires_mineralN) - drate

  Jacobian(ires_plantN,ires_mineralN) = Jacobian(ires_plantN,ires_mineralN) + drate

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
