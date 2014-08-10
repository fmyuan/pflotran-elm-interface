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
    PetscInt :: ispec_ngasdeni

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

    PetscReal :: half_saturation
    PetscInt :: temperature_response_function
    PetscReal :: Q10
    PetscReal :: k_deni_max                 ! denitrification rate
    PetscReal :: x0eps
    PetscReal :: downreg_no3_0  ! shut off
    PetscReal :: downreg_no3_1  ! start to decrease from 1

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
  DenitrificationCreate%ispec_ngasdeni = 0

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
  DenitrificationCreate%half_saturation = -1.0d-15
  DenitrificationCreate%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4
  DenitrificationCreate%Q10 = 1.5d0
  DenitrificationCreate%k_deni_max = 2.5d-6  ! denitrification rate
  DenitrificationCreate%x0eps = 1.0d-20
  DenitrificationCreate%downreg_no3_0 = -1.0d-9
  DenitrificationCreate%downreg_no3_1 = 1.0d-7

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
        call InputReadDouble(input,option,this%k_deni_max)
        call InputErrorMsg(input,option,'k_deni_max', &
                 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,REACTION')
      case('NITRATE_HALF_SATURATION')
        call InputReadDouble(input,option,this%half_saturation)
        call InputErrorMsg(input,option,'nitrate half-saturation', &
                 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,REACTION')
      case('X0EPS')
        call InputReadDouble(input,option,this%x0eps)
        call InputErrorMsg(input,option,'x0eps', &
                  'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION,REACTION')
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
 
  word = 'NO3-'
  this%ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  if(this%ispec_no3 < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION: ' // &
                        ' NO3- is not specified in the input file.'
     call printErrMsg(option)
  endif

  word = 'N2O(aq)'
  this%ispec_n2o = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'N2(aq)'
  this%ispec_n2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if(this%ispec_n2 < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION: ' // &
                        ' N2(aq) is not specified in the input file.'
     call printErrMsg(option)
  endif

  word = 'NGASdeni'
  this%ispec_ngasdeni = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
end subroutine DenitrificationSetup

!********************************************************************************************!
subroutine DenitrificationReact(this,Residual,Jacobian,compute_derivative, &
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

  class(reaction_sandbox_denitrification_type) :: this
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
  PetscReal :: L_water
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  PetscReal :: temp_real

  PetscInt :: ires_no3, ires_n2o, ires_n2
  PetscInt :: ires_ngasdeni, ires

  PetscScalar, pointer :: bsw(:)
  PetscScalar, pointer :: bulkdensity(:)

  PetscReal :: s_min
  PetscReal :: tc
  PetscReal :: f_t, f_w

  PetscReal :: c_no3      ! mole/L
  PetscReal :: f_no3         ! no3 / (half_saturation + no3)
  PetscReal :: d_no3         ! half_saturation/(no3 + half_saturation)^2
  PetscReal :: rate_deni, drate_deni
  PetscReal :: saturation
  PetscInt, parameter :: iphase = 1
  PetscReal :: xxx, delta, regulator, dregulator

!---------------------------------------------------------------------------------

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3

  ! indices for C and N species
  ires_no3 = this%ispec_no3
  ires_n2o = this%ispec_n2o
  ires_n2 = this%ispec_n2
  ires_ngasdeni = this%ispec_ngasdeni + reaction%offset_immobile

! denitrification (Dickinson et al. 2002)
  if(this%ispec_n2 < 0) return

#ifdef CLM_PFLOTRAN
  ghosted_id = option%iflag
  call VecGetArrayReadF90(clm_pf_idata%bsw_pf, bsw, ierr)
  CHKERRQ(ierr)
  temp_real = bsw(ghosted_id)
  call VecRestoreArrayReadF90(clm_pf_idata%bsw_pf, bsw, ierr)
  CHKERRQ(ierr)

#if defined(CHECK_DATAPASSING) && defined(CLM_PFLOTRAN)
  write(option%myrank+200,*) 'checking pflotran-bgc-denitr:', &
    'rank=',option%myrank, 'ghosted_id=',ghosted_id, 'porosity=', porosity, &
    'lsat=',global_auxvar%sat(iphase), &
    'soilt=',global_auxvar%temp, &
    'bsw(ghosted_id)=',temp_real
#endif

#else
  temp_real = 1.0d0
#endif

  tc = global_auxvar%temp
  f_t = exp(0.08d0 * (tc - 25.d0))

  saturation = global_auxvar%sat(1)
  s_min = 0.6d0
  f_w = 0.d0
  if(saturation > s_min) then
     f_w = (saturation - s_min)/(1.0d0 - s_min)
     f_w = f_w ** temp_real
  endif

  c_no3 = rt_auxvar%total(this%ispec_no3, iphase)

  if (this%half_saturation > 0.0d0) then
    temp_real = c_no3 + this%half_saturation
    f_no3 = c_no3 * c_no3 / temp_real
    d_no3 = c_no3 * (c_no3 + &
             2.d0 * this%half_saturation) / temp_real /temp_real
  else
    f_no3 = c_no3 - this%x0eps
    d_no3 = 1.0d0
  endif

  if (this%downreg_no3_0 > 0.0d0) then
    ! additional down regulation for denitrification
    if (c_no3 <= this%downreg_no3_0) then
      regulator = 0.0d0
      dregulator = 0.0d0
    elseif (c_no3 >= this%downreg_no3_1) then
      if (this%downreg_no3_1 - this%downreg_no3_0 > 1.0d-20) then
        regulator = 1.0d0
        dregulator = 0.0d0
      else
        regulator = 0.0d0
        dregulator = 0.0d0
      endif
    else
      xxx = c_no3 - this%downreg_no3_0
      delta = this%downreg_no3_1 - this%downreg_no3_0
      if (delta >= 1.0d-20) then
        regulator = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
        dregulator = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx / delta / delta
      else
        regulator = 0.0d0
        dregulator = 0.0d0
      endif
    endif

    ! rate = rate_orginal * regulator
    ! drate = drate_original * regulator + rate_orginal * dregulator
    d_no3 = d_no3 * regulator + f_no3 * dregulator
    f_no3 = f_no3 * regulator

  endif

  if(f_t > 0.d0 .and. f_w > 0.d0 .and. c_no3>this%x0eps) then
     rate_deni = this%k_deni_max * f_t * f_w * L_water * f_no3

     Residual(ires_no3) = Residual(ires_no3) + rate_deni
     Residual(ires_n2) = Residual(ires_n2) - 0.5d0*rate_deni
    
     if(this%ispec_ngasdeni > 0) then
        Residual(ires_ngasdeni) = Residual(ires_ngasdeni) - 0.5d0*rate_deni
     endif

    if (compute_derivative) then

       drate_deni = this%k_deni_max * f_t * f_w * L_water * d_no3 

       Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) + drate_deni * &
        rt_auxvar%aqueous%dtotal(this%ispec_no3,this%ispec_no3,iphase)

       Jacobian(ires_n2,ires_no3)=Jacobian(ires_n2,ires_no3)-0.5d0*drate_deni * &
        rt_auxvar%aqueous%dtotal(this%ispec_n2,this%ispec_no3,iphase)
    
       if(this%ispec_ngasdeni > 0) then
         Jacobian(ires_ngasdeni,ires_no3)=Jacobian(ires_ngasdeni,ires_no3)-0.5d0*drate_deni
       endif
    endif
  endif

  do ires=1, reaction%ncomp
    temp_real = Residual(ires)

    if (abs(temp_real) > huge(temp_real)) then
      write(option%fid_out, *) 'infinity of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: DENITRIFICATION'
      option%io_buffer = ' checking infinity of Residuals matrix @ DenitrificationReact '
      call printErrMsg(option)
    endif
  enddo

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
