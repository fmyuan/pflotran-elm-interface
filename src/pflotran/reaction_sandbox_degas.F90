module Reaction_Sandbox_degas_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_degas_type
    PetscInt  :: ispec_co2a
    PetscInt  :: ispec_co2g    
    PetscInt  :: ispec_proton
    PetscInt  :: ispec_himm
    PetscReal :: k_kinetic_co2
    PetscReal :: k_kinetic_h
    PetscBool :: b_fixph
    PetscReal :: fixph

  contains
    procedure, public :: ReadInput => degasRead
    procedure, public :: Setup => degasSetup
    procedure, public :: Evaluate => degasReact
    procedure, public :: Destroy => degasDestroy
  end type reaction_sandbox_degas_type

  public :: degasCreate

contains

! ************************************************************************** !
!
! degasCreate: Allocates degas reaction object.
! author: Guoping Tang
! date: 04/03/2014
!
! ************************************************************************** !
function degasCreate()

  implicit none
  
  class(reaction_sandbox_degas_type), pointer :: degasCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(degasCreate)
  degasCreate%ispec_co2a = 0
  degasCreate%ispec_co2g = 0
  degasCreate%ispec_proton = 0
  degasCreate%k_kinetic_co2 = 1.d-5
  degasCreate%k_kinetic_h = 1.d-5
  degasCreate%fixph = 6.5d0
  degasCreate%b_fixph = PETSC_FALSE
  nullify(degasCreate%next)  
      
end function degasCreate

! ************************************************************************** !
!
! degasRead: Reads input deck for degas reaction parameters (if any)
! author: Guoping Tang
! date: 04/03/2014
!
! ************************************************************************** !
subroutine degasRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_degas_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,DEGAS')
    call StringToUpper(word)   

    select case(trim(word))
      case('KINETIC_CONSTANT_CO2')
         call InputReadDouble(input,option,this%k_kinetic_co2)
         call InputErrorMsg(input,option,'CO2 degas kinetic rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,DEGAS,REACTION')
      case('KINETIC_CONSTANT_H+')
         call InputReadDouble(input,option,this%k_kinetic_h)
         call InputErrorMsg(input,option,'H+ fix pH  kinetic rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,DEGAS,REACTION')
      case('FIXPH')
         call InputReadDouble(input,option,this%fixph)
         call InputErrorMsg(input,option,'fix ph', &
                     'CHEMISTRY,REACTION_SANDBOX,DEGAS,REACTION')
         this%b_fixph = PETSC_TRUE
      case default
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DEGAS,' // &
            'REACTION keyword: ' // trim(word) // ' not recognized.'
          call printErrMsg(option)
    end select
  enddo
  
end subroutine degasRead

! ************************************************************************** !
!
! degasSetup: Sets up the degas reaction either with parameters either
!                read from the input deck or hardwired.
! author: Guoping Tang
! date: 04/03/2014
!
! ************************************************************************** !
subroutine degasSetup(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module
  use Immobile_Aux_module, only : GetImmobileSpeciesIDFromName 

  implicit none
  
  class(reaction_sandbox_degas_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  word = 'HCO3-'
  this%ispec_co2a = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'H+'
  this%ispec_proton = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'Cimm'
  this%ispec_co2g = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  word = 'Himm'
  this%ispec_himm = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
end subroutine degasSetup

! ************************************************************************** !
!
! degasReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 04/03/2014
!
! ************************************************************************** !
subroutine degasReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option,local_id)

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only : material_auxvar_type
  use co2eos_module, only: duanco2, HENRY_co2_noderiv     ! co2eos.F90


#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif
  
  implicit none

#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif
  
  class(reaction_sandbox_degas_type) :: this  
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
  PetscReal, parameter :: H2O_kg_mol = 18.01534d-3 ! kg mol-1 for pure water


  PetscInt :: ires_co2a, ires_co2g
  PetscInt :: ires_proton, ires_himm

  PetscReal :: c_hco3     ! mole/L
  PetscReal :: c_hco3_eq      ! mole/m3
  PetscReal :: tc, pres, rate, drate
  PetscReal :: co2_p, co2_rho, co2_fg, co2_xphi
  PetscReal :: xmole_co2, xmass_co2, co2_henry, co2_poyn
  PetscReal :: temp_real, s_hco3
  PetscReal :: c_h, c_h_fix
  PetscReal :: convert_molal_to_molar
  PetscReal :: xmass
 
  xmass = 1.d0

  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  
  if (reaction%initialize_with_molality) then
    convert_molal_to_molar = global_auxvar%den_kg(iphase)*xmass/1000.d0
!    convert_molar_to_molal = 1.d0
  else
    convert_molal_to_molar = 1.d0
!    convert_molar_to_molal = 1000.d0/global_auxvar%den_kg(iphase)/xmass
  endif

#ifdef CLM_PFLOTRAN 
  tc = global_auxvar%temp(1) 
  pres = global_auxvar%pres(1) 
#else
  tc = option%reference_temperature
  pres = option%reference_pressure
#endif

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  ires_co2a = this%ispec_co2a
  ires_co2g = this%ispec_co2g + reaction%offset_immobile

  c_hco3 = rt_auxvar%total(this%ispec_co2a, iphase)

! co2_p should be calculated based on pres, Fengming uses atmospheric CO2 forc_pco2
! here co2_p = 350.0d-6 for demonstration 
  co2_p = 350.0d-6 * 101325.d0

  call duanco2(tc,co2_p, co2_rho, co2_fg, co2_xphi)     ! only need 'co2_xphi' (fugacity coefficient) for the follwong call

  call Henry_CO2_noderiv(xmole_co2,xmass_co2,tc,co2_p,co2_xphi,co2_henry,co2_poyn)   ! 'xmolco2': mol fraction; 'xmco2': mass fraction (CO2:CO2+H2O)

  c_hco3_eq =  xmole_co2/H2O_kg_mol

  temp_real = volume * 1000.0d0 * porosity * global_auxvar%sat(1)
!  rate = this%k_kinetic_co2 * (c_hco3 - c_hco3_eq) * temp_real 
  rate = this%k_kinetic_co2 * (c_hco3/c_hco3_eq - 1.0d0) * temp_real 

!  if(rate > 0.0) then

    Residual(ires_co2a) = Residual(ires_co2a) + rate
    Residual(ires_co2g) = Residual(ires_co2g) - rate

    if (compute_derivative) then
!     drate = -1.0d0 * this%k_kinetic_co2 * temp_real 
     drate = this%k_kinetic_co2 /c_hco3_eq * temp_real 

!     if(rate < 0.0) then
!        drate = drate * s_hco3/(1.0d-15 + s_hco3)
!     endif

     Jacobian(ires_co2a,ires_co2a) = Jacobian(ires_co2a,ires_co2a) + drate * &
        rt_auxvar%aqueous%dtotal(this%ispec_co2a,this%ispec_co2a,iphase)

     Jacobian(ires_co2g,ires_co2a) = Jacobian(ires_co2g,ires_co2a) - drate

!     if(rate < 0.0) then
!        drate = this%k_kinetic_co2 * (c_hco3 - c_hco3_eq) * temp_real * &
!                1.0d-15 /(1.0d-15 + s_hco3)/(1.0d-15+s_hco3) 
!
!        Jacobian(ires_co2a,ires_co2g) = Jacobian(ires_co2a,ires_co2g) + drate

!        Jacobian(ires_co2g,ires_co2g) = Jacobian(ires_co2g,ires_co2g) - drate
!     endif

    endif
! endif

  if(this%b_fixph) then
    ires_proton = this%ispec_proton
    ires_himm = this%ispec_himm + reaction%offset_immobile

    c_h = rt_auxvar%pri_molal(this%ispec_proton) * convert_molal_to_molar

    c_h_fix = 10.0d0 ** (-1.0d0 * this%fixph) / rt_auxvar%pri_act_coef(this%ispec_proton)

!    rate = this%k_kinetic_h * (c_h - c_h_fix) * temp_real
    rate = this%k_kinetic_h * (c_h/c_h_fix - 1.0d0) * temp_real
     
    Residual(ires_proton) = Residual(ires_proton) + rate
    Residual(ires_himm) = Residual(ires_himm) - rate

    if (compute_derivative) then

!       drate = -1.0d0 * this%k_kinetic_h * temp_real 
       drate = this%k_kinetic_h * convert_molal_to_molar /c_h_fix * temp_real 

       Jacobian(ires_proton,ires_proton) = Jacobian(ires_proton,ires_proton) + drate

       Jacobian(ires_himm,ires_proton) = Jacobian(ires_proton,ires_himm) - drate

    endif

  endif 

end subroutine degasReact

! ************************************************************************** !
!
! degasDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Guoping Tang
! date: 04/03/2014
!
! ************************************************************************** !
subroutine degasDestroy(this)

  implicit none
  
  class(reaction_sandbox_degas_type) :: this  

end subroutine degasDestroy

end module Reaction_Sandbox_degas_class
