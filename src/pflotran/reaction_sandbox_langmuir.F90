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
    character(len=MAXWORDLENGTH) :: name_aq
    character(len=MAXWORDLENGTH) :: name_sorb
    PetscInt  :: ispec_aq
    PetscInt  :: ispec_sorb
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
! Revised by Fengming YUAN
! ************************************************************************** !
function LangmuirCreate()

  implicit none
  
  class(reaction_sandbox_langmuir_type), pointer :: LangmuirCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(LangmuirCreate)
  LangmuirCreate%name_aq = ''
  LangmuirCreate%name_sorb = ''
  LangmuirCreate%ispec_aq = 0
  LangmuirCreate%ispec_sorb = 0
  LangmuirCreate%k_kinetic = 1.d-5
  LangmuirCreate%k_equilibrium = 2.5d+3
  LangmuirCreate%s_max = 1.0d-3
  LangmuirCreate%x0eps = 1.0d-21
  nullify(LangmuirCreate%next)  
      
end function LangmuirCreate

! ************************************************************************** !
!
! LangmuirRead: Reads input deck for langmuir reaction parameters (if any)
! author: Guoping Tang
! date: 09/09/2013
! revised by Fengming Yuan @07-10-2015
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
  character(len=MAXWORDLENGTH) :: word, name
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,LANGMUIR')
    call StringToUpper(word)   

    select case(trim(word))
     case('NAME_AQ')
         call InputReadWord(input,option,name, PETSC_TRUE)
         call InputErrorMsg(input,option,'desorbed (aqueous) specie name', &
                     'CHEMISTRY,REACTION_SANDBOX,LANGMUIR,REACTION')
         this%name_aq = trim(name)
     case('NAME_SORB')
         call InputReadWord(input,option,name, PETSC_TRUE)
         call InputErrorMsg(input,option,'sorbed (immobile) specie name', &
                     'CHEMISTRY,REACTION_SANDBOX,LANGMUIR,REACTION')
         this%name_sorb = trim(name)
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
! revised by Fengming Yuan @07-10-2015
!
! ************************************************************************** !
subroutine LangmuirSetup(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName

  implicit none
  
  class(reaction_sandbox_langmuir_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  word = this%name_aq
  this%ispec_aq = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = this%name_sorb
  this%ispec_sorb = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)

  if(this%ispec_aq < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,LANGMUIR,' // &
            'desorbed specie name is not defined as primary aqueous species!.'
     call printErrMsg(option)
  endif

  if(this%ispec_sorb < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,LANGMUIR,' // &
            'sorbed specie name is not defined as immobile species!.'
     call printErrMsg(option)
  endif

end subroutine LangmuirSetup

! ************************************************************************** !
!
! LangmuirReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 09/09/2013
! revised by Fengming Yuan @07-10-2015
!
! ************************************************************************** !
subroutine LangmuirReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only : material_auxvar_type
  use CLM_RspFuncs_module

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
  PetscErrorCode :: ierr

  PetscInt, parameter :: iphase = 1

  PetscInt :: ires_aq, ires_sorb, ires

  PetscReal :: c_aq           ! mole/L
  PetscReal :: c_aq_eq        ! mole/L
  PetscReal :: c_sorb         ! mole/m3
  PetscReal :: c_sorb_eq      ! mole/m3
  PetscReal :: rate, drate_daq, drate_dsorb
  PetscReal :: Lwater         ! litres h2o
  PetscReal :: temp_real

  !-------------------------------------------------------------------------------

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume
  Lwater = volume * 1000.0d0 * porosity * global_auxvar%sat(1)     ! Litres H2O

  ires_aq = this%ispec_aq
  ires_sorb = this%ispec_sorb + reaction%offset_immobile

  c_aq    = rt_auxvar%total(this%ispec_aq, iphase)      ! moles/L
  c_sorb  = rt_auxvar%immobile(this%ispec_sorb)         ! moles/m3

  ! c_sorb_eq = s_max * (Kl*c_eq/(1+Kl*c_eq))
  !c_aq_eq = c_sorb / (this%s_max - c_sorb) / this%k_equilibrium        ! here will cause infinity error
  c_sorb_eq = this%s_max * (this%k_equilibrium*c_aq/(1.d0+this%k_equilibrium*c_aq))

  ! the following IS wrong - there is no way that the aq. concentration over eq. conc WILL throw away to absorbed state
  ! esp. when there is a max. capacity of sorption
  !rate = this%k_kinetic * (c_aq - c_aq_eq) * Lwater      ! moles/L/s --> moles/s
  rate = this%k_kinetic * (c_sorb_eq - c_sorb) * volume   ! positive, when sorption; negative when desorption. ANd moles/m3/s --> moles/s

  if (abs(rate)<this%x0eps/option%tran_dt) return

  ! --------------------------------------------------
  Residual(ires_aq)   = Residual(ires_aq) - rate
  Residual(ires_sorb) = Residual(ires_sorb) + rate

  if (compute_derivative) then
     !
     ! drate/d(Csorb*volume)
     drate_dsorb = -1.d0 * this%k_kinetic
     !
     ! drate/d(Caq*Lwater): this%k_kinetic*this%s_max*this%k_equilibrium*volume*
     !               d(Caq*Lwater/(Lwater+k*Caq*Lwater)) / d(Caq*Lwater)
     temp_real = this%k_kinetic*this%s_max*this%k_equilibrium*volume
     drate_daq = temp_real * Lwater/(this%k_equilibrium*c_aq+Lwater)/(this%k_equilibrium*c_aq+Lwater)

     Jacobian(ires_aq,ires_aq) = Jacobian(ires_aq,ires_aq) - drate_daq * &
        rt_auxvar%aqueous%dtotal(this%ispec_aq,this%ispec_aq,iphase)

     Jacobian(ires_sorb,ires_aq) = Jacobian(ires_sorb,ires_aq) + drate_daq

     Jacobian(ires_aq,ires_sorb) = Jacobian(ires_aq,ires_sorb) - drate_dsorb
     Jacobian(ires_sorb,ires_sorb) = Jacobian(ires_sorb,ires_sorb) + drate_dsorb

  endif

!#ifdef DEBUG
  do ires=1, reaction%ncomp
    temp_real = Residual(ires)

    if (abs(temp_real) > huge(temp_real)) then
      write(option%fid_out, *) 'infinity of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: LANGMUIR'
      option%io_buffer = ' checking infinity of Residuals matrix  @ LangmuiReact '
      call printErrMsg(option)
    endif

    if (temp_real /= temp_real) then
      write(option%fid_out, *) 'NaN of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: LANGMUIR'
      option%io_buffer = ' checking NaN of Residuals matrix  @ LangmuiReact '
      call printErrMsg(option)
    endif

  enddo
!#endif

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
