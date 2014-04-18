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

  word = 'CO2(aq)'
  this%ispec_co2a = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)


  if(this%ispec_co2a < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DEGAS,' // &
            'aqueous species, either CO2(aq) is not defined!'
     call printErrMsg(option)
  endif

  ! (TODO) currently, gas CO2 related process not ready, so it's assumed as immobile species.
  word = 'CO2imm'
  this%ispec_co2g = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  if(this%ispec_co2g < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DEGAS,' // &
            'gas species CO2imm is not defined!'
     call printErrMsg(option)
  endif

  if(this%b_fixph) then
     word = 'H+'
     this%ispec_proton = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

     word = 'Himm'
     this%ispec_himm = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
   
     if(this%ispec_proton < 0) then
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DEGAS,' // &
            'H+ is not defined even though pH needs to be fixed!'
        call printErrMsg(option)
     endif
 
     if(this%ispec_himm < 0) then
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DEGAS,' // &
            'Himm is not defined even though pH needs to be fixed!'
        call printErrMsg(option)
     endif
  endif

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
  use TH_Aux_module
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
  type(TH_auxvar_type) :: th_auxvar
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
  PetscReal, parameter :: rgas = 8.3144621d0       ! m3 Pa K-1 mol-1

  PetscInt :: ires_co2a, ires_co2g
  PetscInt :: ires_proton, ires_himm

  PetscReal :: c_hco3         ! mole/L
  PetscReal :: c_hco3_eq      ! mole/m3
  PetscReal :: tc, pres, lsat, isat, rate, drate
  PetscReal :: co2_p, co2_rho, co2_fg, co2_xphi
  PetscReal :: air_vol, air_molar, co2_molar
  PetscReal :: xmole_co2, xmass_co2, co2_henry, co2_poyn
  PetscReal :: temp_real, s_hco3
  PetscReal :: c_h, c_h_fix
  PetscReal :: convert_molal_to_molar
  PetscReal :: xmass
 
  xmass = 1.d0

  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  
  if (reaction%initialize_with_molality) then
    convert_molal_to_molar = global_auxvar%den_kg(iphase)*xmass/1000.d0
  else
    convert_molal_to_molar = 1.d0
  endif

  !default values for calculating gas solubility
  tc = option%reference_temperature
  pres = option%reference_pressure
  lsat = 0.50d0  ! 50% saturation assumed as default
  isat = 0.d0
  co2_p = 350.0d-6 * option%reference_pressure
  if (option%iflowmode == RICHARDS_MODE .or. &
      option%iflowmode == TH_MODE) then

      pres = global_auxvar%pres(1)      ! water pressure or total (air)gas pressure
      lsat = global_auxvar%sat(1)
      if (option%iflowmode == TH_MODE) then
         tc = global_auxvar%temp(1)
         !isat = th_auxvar%sat_ice       ! (TODO) not yet figure out how to point to 'th_auxvar'
      endif
#ifdef CLM_PFLOTRAN
  elseif (option%ntrandof.gt.0 ) then
      pres = global_auxvar%pres(1)      ! water pressure or total (air)gas pressure
      lsat = global_auxvar%sat(1)
      tc = global_auxvar%temp(1)
  else
      option%io_buffer='reaction_sandbox_degas ' // &
                 'not supported for the modes applied for.'
      call printErrMsg(option)
#endif
  endif

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  ires_co2a = this%ispec_co2a
  ires_co2g = this%ispec_co2g + reaction%offset_immobile

  c_hco3 = rt_auxvar%total(this%ispec_co2a, iphase)

#ifdef CLM_PFLOTRAN
  ! resetting 'co2g' from CLM after adjusting
  if (this%ispec_co2g > 0) then
     !air_vol = max(0.01d0, porosity * (1.d0-lsat-isat))               ! min. 0.01 to avoid math. issue
     air_vol = 1.d0                                                   ! atm. air volume fraction (no adjusting of soil air volume)
     co2_molar = rt_auxvar%immobile(this%ispec_co2g)/air_vol          ! molCO2/m3 bulk soil --> mol/m3 air space
     air_molar = pres/rgas/(tc+273.15d0)                              ! molAir/m3
     co2_p = co2_molar/air_molar*pres                                 ! mole fraction --> Pa
  endif
#endif

  temp_real = max(min(tc,65.d0), 1.d-20)                                   ! 'duanco2' only functions from 0 - 65oC

  call duanco2(temp_real,co2_p, co2_rho, co2_fg, co2_xphi)     ! only need 'co2_xphi' (fugacity coefficient) for the follwong call

  call Henry_CO2_noderiv(xmole_co2,xmass_co2,tc,co2_p,co2_xphi,co2_henry,co2_poyn)   ! 'xmolco2': mol fraction; 'xmco2': mass fraction (CO2:CO2+H2O)

  c_hco3_eq =  xmole_co2/H2O_kg_mol

  if (PETSC_FALSE) then
    rate = this%k_kinetic_co2 * (c_hco3/c_hco3_eq - 1.0d0) * temp_real
  else
    rate = this%k_kinetic_co2 * (c_hco3 - c_hco3_eq) * temp_real
  endif

! degas occurs only when oversaturated, not considering uptake of CO2 from atmosphere
  if(abs(rate) > 1.0d-20) then

    Residual(ires_co2a) = Residual(ires_co2a) + rate
    Residual(ires_co2g) = Residual(ires_co2g) - rate

    if (compute_derivative) then
        if (PETSC_FALSE) then
            drate = this%k_kinetic_co2 /c_hco3_eq * temp_real 
        else
            drate = this%k_kinetic_co2 * temp_real 
        endif

        Jacobian(ires_co2a,ires_co2a) = Jacobian(ires_co2a,ires_co2a) + drate * &
        rt_auxvar%aqueous%dtotal(this%ispec_co2a,this%ispec_co2a,iphase)

        Jacobian(ires_co2g,ires_co2a) = Jacobian(ires_co2g,ires_co2a) - drate

    endif
  endif

  ! the following is for setting pH at a fixed value from input
  ! (TODO) if PFLOTRAN H+/OH- related processes are in place, this should be removed
  if(this%b_fixph) then
    ires_proton = this%ispec_proton
    ires_himm = this%ispec_himm + reaction%offset_immobile

    c_h = rt_auxvar%pri_molal(this%ispec_proton) * convert_molal_to_molar

    c_h_fix = 10.0d0 ** (-1.0d0 * this%fixph) / rt_auxvar%pri_act_coef(this%ispec_proton)

    if(PETSC_FALSE) then
      rate = this%k_kinetic_h * (c_h/c_h_fix - 1.0d0) * temp_real
    else
      rate = this%k_kinetic_h * (c_h - c_h_fix) * temp_real
    endif  
     
    if(abs(rate) > 1.0d-20) then
       Residual(ires_proton) = Residual(ires_proton) + rate
       Residual(ires_himm) = Residual(ires_himm) - rate

       if (compute_derivative) then

       if (PETSC_FALSE) then
           drate = this%k_kinetic_h * convert_molal_to_molar /c_h_fix * temp_real
       else
           drate = this%k_kinetic_h * convert_molal_to_molar * temp_real
       endif


       Jacobian(ires_proton,ires_proton) = Jacobian(ires_proton,ires_proton) + drate

         Jacobian(ires_himm,ires_proton) = Jacobian(ires_proton,ires_himm) - drate

      endif
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

! ************************************************************************** !

  subroutine weiss_price_n2o (tt, tp, ts, pn2o, xmole, xmass, kh, fg, phi)
  !
  ! Weiss and Price, 1980. Nitrous oxide solubility in water and seawater. Marine
  ! chemistry, 8(1980)347-359
  !
  ! Input: tt   [C]        temperature
  !        tp   [Pa]       total air pressure
  !        ts   [%o]       total salinity (parts per thousands)
  !        pn2o [Pa]       N2O partial pressure
  ! Output: fg  [Pa]       N2O fugacity
  !         phi [-]        N2O fugacity coefficient
  !         kh  [Pa]       N2O Henry's Law Constant
  !         xmole  [-]     mole fraction of N2O solubility
  !         xmass  [-]     mass fraction of N2O solubility

      implicit none

#include "finclude/petscsys.h"

      PetscReal, intent(in) :: tt, tp, ts, pn2o
      PetscReal, intent(out):: xmole, xmass, kh, fg, phi

      PetscReal, parameter :: atm  = 1.0314d5         ! Pa of 1 atm
      PetscReal, parameter :: rgas = 0.08205601       ! L atm mol-1 K-1
      PetscReal, parameter :: xmwn2o  = 44.01287d-3   ! kg mol-1
      PetscReal, parameter :: xmwh2o  = 18.01534d-3   ! kg mol-1

      PetscReal :: tk, tk2, tk_100k
      PetscReal :: p_rt
      PetscReal :: x1,x2
      PetscReal :: k0, epsilon, bt
      PetscReal :: vbar, cn2o

      PetscReal :: a1,a2,a3,b1,b2,b3

      ! empirical parameters for temperature effect
      data a1    /-62.7076d0/
      data a2    / 97.3066d0/
      data a3    / 24.1406d0/
      ! empirical parameters for salt-water effect
      data b1    /-0.058420d0/
      data b2    / 0.033193d0/
      data b3    /-0.0051313d0/

      ! variable values transformation
      tk = tt + 273.15d0
      tk2= tk*tk
      tk_100k = tk/100.d0
      p_rt = (tp/atm)/rgas/tk                       ! mol L-1
      p_rt = p_rt/1000.d0                           ! mol cm-3  (this conversion needed? - the original paper did say)
      x1 = pn2o/tp                                  ! mole fractions of binary mixture of n2o-air
      x2 = 1.d0 - x1

      ! fugacity of n2o gas
      epsilon = 65.0d0 - 0.1338d0*tk                ! eq.(6): cm3/mol
      bt = -905.95d0+4.1685d0*tk-0.0052734d0*tk2   ! eq.(4): cm3/mol
      fg = pn2o*dexp((bt+2.d0*x2*x2*epsilon)*p_rt)   ! eq.(5): Pa (so, pn2o is moist-air based)

      ! fugacity coefficient
      phi = fg / pn2o

      ! K0 in Weiss and Price (1980)'s paper
      k0 = a1+a2/tk_100k+a3*dlog(tk_100k) + &
           ts*(b1+b2*tk_100k+b3*tk_100k*tk_100k)    ! eq. (12): ln(mol/L/atm)
      k0 = dexp(k0)                                 ! mol/L/atm
      k0 = k0/atm                                   ! mol/L/pa

      ! solubility
      vbar = dexp((1.0d0-tp/atm)*32.3d-3/rgas/tk)    ! 32.3 is in cm3/mol (unit inconsistency??)
      cn2o = k0*fg*vbar                              ! eq. (1): mol/L
      xmole = cn2o/(1.d0/xmwh2o)                     ! mole fraction: assuming solution volume is all water (1L=1kgH2O)
      xmass = xmole*xmwn2o/(xmole*xmwn2o+(1.d0-xmole)*xmwh2o)    ! mass fraction

      ! Henry's Law constant: kh = pn2o/[xmole]
      kh = 1.d0/k0                      ! Pa L mol-1: pn2o in pa, solubility in mol/L
      kh = kh*(1.0d0/xmwh2o)            ! Pa (pa mol mol-1): pn2o in pa, solubility in mole fraction (1Lsolution = 1kgH2O)

      return

    end subroutine weiss_price_n2o

end module Reaction_Sandbox_degas_class
