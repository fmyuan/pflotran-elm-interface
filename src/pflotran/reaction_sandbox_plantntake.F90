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
    PetscReal :: x0eps_no3
    ! additional down regulation for plant NO3- uptake with NO3 uptake f(NO3-) 
    ! f(NO3-) = 0 for NO3- <= downreg_no3_0 
    ! f(NO3-) = 1 for NO3- >= downreg_no3_1
    ! f(NO3-) = 1 - [1 - (x/d)^2]^2 
    ! with x = c - downreg_no3_0, d = downreg_no3_1 - downreg_0
    PetscReal :: downreg_no3_0  ! shut off
    PetscReal :: downreg_no3_1  ! start to decrease from 1
    PetscReal :: downreg_nh3_0  ! shut off
    PetscReal :: downreg_nh3_1  ! start to decrease from 1

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
  PlantNTakeCreate%x0eps_no3  = 1.d-20
  PlantNTakeCreate%downreg_no3_0 = -1.0d-9
  PlantNTakeCreate%downreg_no3_1 = 1.0d-7
  PlantNTakeCreate%downreg_nh3_0 = -1.0d-9
  PlantNTakeCreate%downreg_nh3_1 = 1.0d-7
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
          call InputReadDouble(input,option,this%x0eps)
          call InputErrorMsg(input,option,'x0eps_nh4', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case('X0EPS_NO3')
          call InputReadDouble(input,option,this%x0eps_no3)
          call InputErrorMsg(input,option,'x0eps_no3', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
      case('DOWNREGULATE_NH4')
        call InputReadDouble(input,option,this%downreg_nh3_0)
        call InputErrorMsg(input,option,'downreg_nh3_0', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
        call InputReadDouble(input,option,this%downreg_nh3_1)
        call InputErrorMsg(input,option,'downreg_nh3_1', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
        if (this%downreg_nh3_0 > this%downreg_nh3_1) then
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,' // &
            'NH4+ down regulation cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
      case('DOWNREGULATE_NO3')
        call InputReadDouble(input,option,this%downreg_no3_0)
        call InputErrorMsg(input,option,'downreg_no3_0', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
        call InputReadDouble(input,option,this%downreg_no3_1)
        call InputErrorMsg(input,option,'downreg_no3_1', &
          'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,REACTION')
        if (this%downreg_no3_0 > this%downreg_no3_1) then
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTNTAKE,' // &
            'NO3- down regulation cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
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
  PetscInt :: ispec_nh4in, ispec_no3in, ires_nh4in, ires_no3in
  PetscInt :: ispec_plantndemand, ires_plantndemand, ires

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
  PetscReal :: drate_nplant       !drate_nh4/dnh4+
  PetscReal :: drate_nplant_no3   !drate_no3/dno3-
  PetscReal :: rate_nh4, rate_no3
  PetscReal :: c_plantn, c_plantno3, c_plantnh4, c_plantndemand
  PetscReal :: xxx, delta, regulator, dregulator

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

  word = 'Ain'
  ispec_nh4in = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  word = 'Tin'
  ispec_no3in = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  word = 'Plantndemand'
  ispec_plantndemand = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  ires_nh4in = ispec_nh4in + reaction%offset_immobile
  ires_no3in = ispec_no3in + reaction%offset_immobile
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

    if (this%downreg_nh3_0 > 0.0d0) then
      ! additional down regulation for plant NH4+ uptake
      if (c_nh3 <= this%downreg_nh3_0) then
        regulator = 0.0d0
        dregulator = 0.0d0
      elseif (c_nh3 >= this%downreg_nh3_1) then
        regulator = 1.0d0
        dregulator = 0.0d0
      else
        xxx = c_nh3 - this%downreg_nh3_0
        delta = this%downreg_nh3_1 - this%downreg_nh3_0
        regulator = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
        dregulator = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx / delta
      endif
    
      ! rate = rate_orginal * regulator
      ! drate = drate_original * regulator + rate_orginal * dregulator
      d_nh3 = d_nh3 * regulator + f_nh3 * dregulator

      f_nh3 = f_nh3 * regulator

    endif

    ires_nh3 = ispec_nh3
  endif

  if (ispec_no3 > 0) then
    c_no3 = rt_auxvar%total(ispec_no3, iphase)
    temp_real = c_no3 + this%half_saturation_no3
    f_no3 = c_no3 / temp_real
    d_no3 = this%half_saturation_no3 / temp_real / temp_real

    if (this%downreg_no3_0 > 0.0d0) then
      ! additional down regulation for plant NO3- uptake
      if (c_no3 <= this%downreg_no3_0) then
        regulator = 0.0d0
        dregulator = 0.0d0
      elseif (c_no3 >= this%downreg_no3_1) then
        regulator = 1.0d0
        dregulator = 0.0d0
      else
        xxx = c_no3 - this%downreg_no3_0
        delta = this%downreg_no3_1 - this%downreg_no3_0
        regulator = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
        dregulator = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx / delta
      endif

      ! rate = rate_orginal * regulator
      ! drate = drate_original * regulator + rate_orginal * dregulator
      d_no3 = d_no3 * regulator + f_no3 * dregulator
      f_no3 = f_no3 * regulator

    endif

    ires_no3 = ispec_no3

    if (ispec_nh3 > 0 .and. c_nh3 > this%x0eps) then
      !temp_real = this%inhibition_nh3_no3 + c_nh3
      !f_nh3_inhibit = this%inhibition_nh3_no3/temp_real
      f_nh3_inhibit = 1.0d0 - f_nh3
    else
      f_nh3_inhibit = 1.0d0
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

  rate_nplant = this%rate

  if (ispec_plantndemand > 0) then
    Residual(ires_plantndemand) = Residual(ires_plantndemand) - rate_nplant
  endif

  if(ispec_nh3 > 0 .and. c_nh3 > this%x0eps) then

    drate_nplant = rate_nplant * d_nh3
    rate_nplant = rate_nplant * f_nh3

    Residual(ires_nh3) = Residual(ires_nh3) + rate_nplant
    Residual(ires_plantn) = Residual(ires_plantn) - rate_nplant

    if (ispec_nh4in > 0) then
      Residual(ires_nh4in) = Residual(ires_nh4in) - rate_nplant
    endif

    if (compute_derivative) then
       Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3)+drate_nplant * &
         rt_auxvar%aqueous%dtotal(ispec_nh3,ispec_nh3,iphase)

       Jacobian(ires_plantn,ires_nh3)=Jacobian(ires_plantn,ires_nh3)-drate_nplant

      if (ispec_nh4in > 0) then
        Jacobian(ires_nh4in,ires_nh3)=Jacobian(ires_nh4in,ires_nh3)-drate_nplant
      endif
    endif
  endif

  if(ispec_no3 > 0 .and. c_no3 > this%x0eps_no3) then
    rate_nplant_no3 = this%rate * f_nh3_inhibit * f_no3
    drate_nplant_no3 = this%rate * f_nh3_inhibit * d_no3

    Residual(ires_no3) = Residual(ires_no3) + rate_nplant_no3
    Residual(ires_plantn) = Residual(ires_plantn) - rate_nplant_no3

    if (ispec_no3in > 0) then
      Residual(ires_no3in) = Residual(ires_no3in) - rate_nplant_no3
    endif

    if (compute_derivative) then
     Jacobian(ires_no3,ires_no3)=Jacobian(ires_no3,ires_no3)+drate_nplant_no3* &
       rt_auxvar%aqueous%dtotal(ispec_no3,ispec_no3,iphase)

     Jacobian(ires_plantn,ires_no3)=Jacobian(ires_plantn,ires_no3)-drate_nplant_no3

      if (ispec_no3in > 0) then
        Jacobian(ires_no3in,ires_no3)=Jacobian(ires_no3in,ires_no3)-drate_nplant_no3
      endif

    endif
  endif

  do ires=1, reaction%ncomp
    temp_real = Residual(ires)

    if (abs(temp_real) > huge(temp_real)) then
      write(option%fid_out, *) 'infinity of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: PLANT N UPTAKE'
      option%io_buffer = ' checking infinity of Residuals matrix @ PlantNTakeReact '
      call printErrMsg(option)
    endif
  enddo

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
