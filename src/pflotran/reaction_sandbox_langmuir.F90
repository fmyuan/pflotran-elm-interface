module Reaction_Sandbox_Langmuir_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_langmuir_type
    PetscInt  :: ispec_nh4a
    PetscInt  :: ispec_nh4s
    PetscReal :: k_kinetic
    PetscReal :: k_equilibrium
    PetscReal :: s_max
    PetscReal :: x0eps

  contains
    procedure, public :: ReadInput => LangmuirRead
    procedure, public :: Setup => LangmuirSetup
    procedure, public :: Evaluate => LangmuirReact
    procedure, public :: Destroy => LangmuirDestroy
  end type reaction_sandbox_langmuir_type

  public :: LangmuirCreate

contains

! ************************************************************************** !
!
! LangmuirCreate: Allocates langmuir reaction object.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
function LangmuirCreate()

  implicit none
  
  class(reaction_sandbox_langmuir_type), pointer :: LangmuirCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(LangmuirCreate)
  LangmuirCreate%ispec_nh4a = 0
  LangmuirCreate%ispec_nh4s = 0
  LangmuirCreate%k_kinetic = 1.d-5
  LangmuirCreate%k_equilibrium = 2.5d+3
  LangmuirCreate%s_max = 1.0d-3
  LangmuirCreate%x0eps = 1.0d-20
  nullify(LangmuirCreate%next)  
      
end function LangmuirCreate

! ************************************************************************** !
!
! LangmuirRead: Reads input deck for langmuir reaction parameters (if any)
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine LangmuirRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_langmuir_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,LANGMUIR')
    call StringToUpper(word)   

    select case(trim(word))
     case('EQUILIBRIUM_CONSTANT')
         call InputReadDouble(input,option,this%k_equilibrium)
         call InputErrorMsg(input,option,'langmuir equilibrium constant', &
                     'CHEMISTRY,REACTION_SANDBOX,LANGMUIR,REACTION')
     case('KINETIC_CONSTANT')
         call InputReadDouble(input,option,this%k_kinetic)
         call InputErrorMsg(input,option,'Langmuir kinetic rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,LANGMUIR,REACTION')
     case('S_MAX')
         call InputReadDouble(input,option,this%s_max)
         call InputErrorMsg(input,option,'Langmuir sorption capacity', &
                     'CHEMISTRY,REACTION_SANDBOX,LANGMUIR,REACTION')
      case default
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,LANGMUIR,' // &
            'REACTION keyword: ' // trim(word) // ' not recognized.'
          call printErrMsg(option)
    end select
  enddo
  
end subroutine LangmuirRead

! ************************************************************************** !
!
! LangmuirSetup: Sets up the langmuir reaction either with parameters either
!                read from the input deck or hardwired.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine LangmuirSetup(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module
  use Immobile_Aux_module, only : GetImmobileSpeciesIDFromName 

  implicit none
  
  class(reaction_sandbox_langmuir_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  word = 'NH4+'
  this%ispec_nh4a = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'Nsorb'
  this%ispec_nh4s = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
end subroutine LangmuirSetup

! ************************************************************************** !
!
! LangmuirReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine LangmuirReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option,local_id)

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
  
  class(reaction_sandbox_langmuir_type) :: this  
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
  PetscInt :: local_id
  PetscErrorCode :: ierr

  PetscInt, parameter :: iphase = 1

  PetscInt :: ires_nh4a, ires_nh4s

  PetscReal :: c_nh4      ! mole/L
  PetscReal :: s_nh4      ! mole/m3
  PetscReal :: s_nh4_eq      ! mole/m3
  PetscReal :: kc, rate, drate_a, drate_s

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  ires_nh4a = this%ispec_nh4a
  ires_nh4s = this%ispec_nh4s + reaction%offset_immobile

  c_nh4 = rt_auxvar%total(this%ispec_nh4a, iphase)
  s_nh4 = rt_auxvar%immobile(this%ispec_nh4s)
  c_nh4 = c_nh4 - this%x0eps

  kc = this%k_equilibrium * c_nh4
  s_nh4_eq = this%s_max * kc / (1.0d0 + kc)

  if(s_nh4 >= s_nh4_eq) then
    return
  endif

  rate = this%k_kinetic * (s_nh4_eq - s_nh4) * volume 

  Residual(ires_nh4a) = Residual(ires_nh4a) + rate
  Residual(ires_nh4s) = Residual(ires_nh4s) - rate

  if (compute_derivative) then
     drate_a = this%k_kinetic * this%k_equilibrium * this%s_max / kc / kc
     drate_s = -1.0d0 * this%k_kinetic

     Jacobian(ires_nh4a,ires_nh4a) = Jacobian(ires_nh4a,ires_nh4a) + drate_a * &
        rt_auxvar%aqueous%dtotal(this%ispec_nh4a,this%ispec_nh4a,iphase)

     Jacobian(ires_nh4s,ires_nh4a) = Jacobian(ires_nh4s,ires_nh4a) - drate_a

     Jacobian(ires_nh4a,ires_nh4s) = Jacobian(ires_nh4a,ires_nh4s) + drate_s

     Jacobian(ires_nh4s,ires_nh4s) = Jacobian(ires_nh4s,ires_nh4s) - drate_s

  endif

end subroutine LangmuirReact

! ************************************************************************** !
!
! LangmuirDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine LangmuirDestroy(this)

  implicit none
  
  class(reaction_sandbox_langmuir_type) :: this  

end subroutine LangmuirDestroy

end module Reaction_Sandbox_Langmuir_class
