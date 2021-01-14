module Characteristic_Curves_Base_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private
  
  PetscReal, parameter, public :: DEFAULT_PCMAX = 1.d9
  
  type, public :: polynomial_type
    PetscReal :: low
    PetscReal :: high
    PetscReal :: coefficients(4)
  end type polynomial_type  
  
!-----------------------------------------------------------------------------
!-- Saturation Functions -----------------------------------------------------
!-----------------------------------------------------------------------------
  type, public :: sat_func_base_type
    type(polynomial_type), pointer :: sat_poly
    type(polynomial_type), pointer :: pres_poly
    PetscBool :: analytical_derivative_available
    PetscBool :: calc_int_tension
    PetscReal :: pcmax
    PetscReal :: Sr
  contains

    procedure, public :: Init => SFBaseInit
    procedure, public :: Verify => SFBaseVerify
    procedure, public :: Test => SFBaseTest
    procedure, public :: SetupPolynomials => SFBaseSetupPolynomials
    procedure, public :: CapillaryPressure => SFBaseCapillaryPressure
    procedure, public :: Saturation => SFBaseSaturation
    procedure, public :: D2SatDP2 => SFBaseD2SatDP2
    procedure, public :: CalcInterfacialTension => SFBaseSurfaceTension
  end type sat_func_base_type

!-----------------------------------------------------------------------------
!-- Relative Permeability Functions ------------------------------------------
!-----------------------------------------------------------------------------  
  type, public :: rel_perm_func_base_type
    type(polynomial_type), pointer :: poly
    PetscReal :: Sr
    PetscReal :: Srg
    PetscBool :: analytical_derivative_available
  contains
    procedure, public :: Init => RPFBaseInit
    procedure, public :: Verify => RPFBaseVerify
    procedure, public :: Test => RPFBaseTest
    procedure, public :: SetupPolynomials => RPFBaseSetupPolynomials
    procedure, public :: RelativePermeability => RPFBaseRelPerm
  end type rel_perm_func_base_type
  
  public :: PolynomialCreate, &
            PolynomialDestroy, &
            SFBaseInit, &
            SFBaseVerify, &
            SFBaseCapillaryPressure, &
            SFBaseSaturation, &
            RPFBaseInit, &
            RPFBaseVerify, &
            RPFBaseRelPerm, &
            SaturationFunctionDestroy, &
            PermeabilityFunctionDestroy

contains

! ************************************************************************** !

function PolynomialCreate()

  implicit none
  
  type(polynomial_type), pointer :: PolynomialCreate  

  allocate(PolynomialCreate)
  PolynomialCreate%low = 0.d0
  PolynomialCreate%high = 0.d0
  PolynomialCreate%coefficients(:) = 0.d0
  
end function PolynomialCreate

! ************************************************************************** !
! ************************************************************************** !

subroutine SFBaseInit(this)

  implicit none
  
  class(sat_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%sat_poly)
  nullify(this%pres_poly)
  this%Sr = UNINITIALIZED_DOUBLE
  this%pcmax = DEFAULT_PCMAX
  this%analytical_derivative_available = PETSC_FALSE
  this%calc_int_tension = PETSC_FALSE
  
end subroutine SFBaseInit

! ************************************************************************** !

subroutine SFBaseVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  
  
  if (Uninitialized(this%Sr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                            name)
    call PrintErrMsg(option)
  endif
  
  if ((.not.this%analytical_derivative_available) .and. &
      (.not.option%flow%numerical_derivatives)) then
    option%io_buffer = 'Analytical derivatives are not available for the &
      &capillary pressure - saturation function chosen: ' // &
      trim(name)
    call PrintErrMsg(option)
  endif
  
end subroutine SFBaseVerify

! ************************************************************************** !

subroutine RPFBaseInit(this)

  implicit none
  
  class(rel_perm_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%poly)
  this%Sr = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine RPFBaseInit

! ************************************************************************** !

subroutine RPFBaseVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rel_perm_func_base_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  if (Uninitialized(this%Sr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                            name)
    call PrintErrMsg(option)
  endif
  
  if ((.not.this%analytical_derivative_available) .and. &
      (.not.option%flow%numerical_derivatives)) then
    option%io_buffer = 'Analytical derivatives are not available for the &
      &relative permeability function chosen: ' // trim(name)
    call PrintErrMsg(option)
  endif
  
end subroutine RPFBaseVerify

! ************************************************************************** !

subroutine SFBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing saturation functions

  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  option%io_buffer = 'SF Smoothing not supported for ' // trim(error_string)
  call PrintErrMsg(option)
  
end subroutine SFBaseSetupPolynomials

! ************************************************************************** !

subroutine RPFBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing relative permeability functions

  use Option_module
  
  implicit none
  
  class(rel_perm_func_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  option%io_buffer = 'RPF Smoothing not supported for ' // trim(error_string)
  call PrintErrMsg(option)
  
end subroutine RPFBaseSetupPolynomials

! ************************************************************************** !

subroutine SFBaseCapillaryPressure(this,liquid_saturation, & 
                                   capillary_pressure,dpc_dsatl,option)
  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseCapillaryPressure must be extended.'
  call PrintErrMsg(option)
  
end subroutine SFBaseCapillaryPressure

! ************************************************************************** !

subroutine SFBaseSaturation(this,capillary_pressure, &
                            liquid_saturation,dsat_dpres,option)
  use Option_module

  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseSaturation must be extended.'
  call PrintErrMsg(option)
  
end subroutine SFBaseSaturation

! ************************************************************************** !

subroutine SFBaseD2SatDP2(this,pc,d2s_dp2,option)

  use Option_module

  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseD2SatDP2 must be extended.'
  call PrintErrMsg(option)
  
end subroutine SFBaseD2SatDP2

! ************************************************************************** !

subroutine SFBaseTest(this,cc_name,option)
  use Option_module

  implicit none
  
  class(sat_func_base_type), intent(in) :: this
  character(len=MAXWORDLENGTH), intent(in) :: cc_name
  type(option_type), intent(inout) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: num_values = 1000
  PetscReal, parameter :: dSl = 10d0*epsilon(dSl)
  PetscReal, parameter :: dPc = 1d-6  
  PetscInt :: i ! Loop iterator
  PetscReal ::  Pc,  Sl, dSl_dPc, dPc_dSl ! Analytical values
  PetscReal :: nPc, nSl, nSl_nPc, nPc_nSl ! Numerical derivative values
  PetscReal :: slogPc, sSl ! Spacing

  ! Open sat_from_pc file and write header
  write(string,*) cc_name
  string = trim(cc_name) // '_sat_from_pc.dat'
  open(unit=86,file=string)
  write(86,*) '"capillary pressure", "saturation", "dsat/dpres", &
              &"dsat/dpres_numerical"'

  ! Evaluate saturation function at zero capillary pressure
  ! Evaluate first at 0
  Pc = 0d0
  call this%Saturation(Pc,Sl,dSl_dPc,option)

  ! Use dPc as minimum perturbation at saturation
  nPc = dPc
  call this%Saturation(nPc,nSl,nSl_nPc,option)
  nSl_nPc = (nSl - Sl) / (nPc - Pc)

  write(86,'(4es14.6)') Pc, Sl, dSl_dPc, nSl_nPc

  ! Continue at a minimum capillary pressure of 1 Pa
  Pc = 1d0
  slogPc = exp(log(this%Pcmax/Pc)/num_values) 
  do i = 1, num_values
    call this%Saturation(Pc,Sl,dSl_dPc,option)

    ! Always forward difference with fractional perturbation
    nPc = Pc*dPc
    call this%Saturation(nPc,nSl,nSl_nPc,option)
    nSl_nPc = (nSl - Sl) / (nPc - Pc)

    write(86,'(4es14.6)') Pc, Sl, dSl_dPc, nSl_nPc

    ! Update using logarithmic spacing
    Pc = Pc * slogPc
  end do
  close(86)

  ! Open pc_from_sat file and write header
  write(string,*) cc_name
  string = trim(cc_name) // '_pc_from_sat.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "capillary pressure", "dpc/dsat", &
             &dpc_dsat_numerical"'

  ! Evaluate capillary pressure function
  ! Using linear spacing over 0 to 1
  Sl = 0d0
  sSl = 1d0/num_values
  do i = 0, num_values
    call this%CapillaryPressure(Sl,Pc,dPc_dSl,option)

    ! Use backward difference below 1/2, forward above 1/2
    if (Sl > 0.5d0) then
      nSl = Sl-dSl
    else
      nSl = Sl+dSl
    end if
    call this%CapillaryPressure(nSl,nPc,nPc_nSl,option)
    nPc_nSl = (nPc - Pc) / (nSl - Sl)

    write(86,'(4es14.6)') Sl, Pc, dPc_dSl, nPc_nSl

    ! Update using linear spacing
    Sl = Sl + sSl
  end do

  close(86)

end subroutine

! ************************************************************************** !
! ************************************************************************** !

subroutine RPFBaseRelPerm(this,liquid_saturation,relative_permeability, &
                            dkr_sat,option)
  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'RPFBaseRelPerm must be extended.'
  call PrintErrMsg(option)
  
end subroutine RPFBaseRelPerm

! ************************************************************************** !

subroutine RPFBaseTest(this,cc_name,phase,option)
  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  character(len=MAXWORDLENGTH) :: phase
  type(option_type), intent(inout) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: num_values = 1000
  PetscReal, parameter :: dSl = 1d-6

  PetscInt :: i ! Loop iterator
  PetscReal ::  Sl,  Kr, dKr_dSl ! Analytical values
  PetscReal :: nSl, nKr, nKr_nSl ! Numerical Values
  PetscReal :: sSl ! Spacing
  
  ! Open file and write header 
  write(string,*) cc_name
  string = trim(cc_name) // '_' // trim(phase) // '_rel_perm.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "' // trim(phase) // ' relative permeability", "' &
              // trim(phase) // ' dkr/dsat", "' // trim(phase) // &
              ' dkr/dsat_numerical"'

  ! Evaluate relative permeability function
  ! Using linear spacing over 0 to 1
  sSl = 1d0/num_values
  Sl = 0d0
  nSl = dSl
  do i = 0, num_values
    call this%RelativePermeability(Sl,Kr,dKr_dSl,option)
    call this%RelativePermeability(nSl,nKr,nKr_nSl,option)
    nKr_nSl = (nKr-Kr)/(nSl-Sl)
    write(86,'(4es14.6)') Sl, Kr, dKr_dSl, nKr_nSl

    Sl = Sl + sSl
    if (Sl < 0.5d0) then
      nSl = Sl + dSl
    else
      nSl = Sl - dSl
    end if
  end do

  ! Close file
  close(86)

end subroutine RPFBaseTest

! ************************************************************************** !

subroutine PolynomialDestroy(poly)
  !
  ! Destroys a polynomial smoother
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  !

  implicit none

  type(polynomial_type), pointer :: poly

  if (.not.associated(poly)) return

  deallocate(poly)
  nullify(poly)

end subroutine PolynomialDestroy
  
! ************************************************************************** !

subroutine SaturationFunctionDestroy(sf)
  !
  ! Destroys a saturuation function
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  !

  implicit none

  class(sat_func_base_type), pointer :: sf

  if (.not.associated(sf)) return
  
  call PolynomialDestroy(sf%sat_poly)
  call PolynomialDestroy(sf%sat_poly)
  deallocate(sf)
  nullify(sf)

end subroutine SaturationFunctionDestroy

! ************************************************************************** !

subroutine PermeabilityFunctionDestroy(rpf)
  ! 
  ! Destroys a saturuation function
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(rel_perm_func_base_type), pointer :: rpf
  
  if (.not.associated(rpf)) return
  
  call PolynomialDestroy(rpf%poly)
  deallocate(rpf)
  nullify(rpf)

end subroutine PermeabilityFunctionDestroy

subroutine SFBaseSurfaceTension(this,T,sigma)
  
  !Surface tension of water equation from Revised Release on Surface
  !Tension of Ordinary Water Substance, June 2014. Valid from -25C to
  !373 C
  
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: sigma
  
  PetscReal, parameter :: Tc = 647.096d0
  PetscReal, parameter :: B = 235.8d0
  PetscReal, parameter :: b_2 = -0.625d0
  PetscReal, parameter :: mu = 1.256d0
  PetscReal, parameter :: sigma_base = 0.073d0
  PetscReal :: Temp
  PetscReal :: tao
  
  Temp=T+273.15d0
  
  if (T <= 373.d0) then
    tao = 1.d0-Temp/Tc
    sigma = B*(tao**mu)*(1+b_2*tao)
    sigma = sigma * 1.d-3
  else
    sigma = 0.d0
  endif
  sigma= sigma/sigma_base

  !TOUGH3 way (not pressure-dependent)
  !if (Temp >= 101) sigma = 0
  
end subroutine SFBaseSurfaceTension
  
end module Characteristic_Curves_Base_module
