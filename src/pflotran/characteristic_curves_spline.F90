module Characteristic_Curves_spline_module
#include "petsc/finclude/petscsys.h"

use petscsys ! Necessary for PETSC_TRUE / PETSC_FALSE
use Characteristic_Curves_Base_module ! Needed to define base type
use Option_module ! Needed for unused arguments in Pc and Sl
use slatec_pchip_module
use PFLOTRAN_constants_module ! Needed for UNINITIALIZED_DOUBLE

implicit none
private

! **************************************************************************** !
!
! Author: Matthew Paul
! Date:   04/25/2023
!
! This module contains the methods to construct, evaluate, and deconstruct 
! cubic splines for any monotonic capillary pressure or relative permeability
! function or data set. This can be done to use available data or to accelerate
! computationally expensive but analytical functions.

! Cubic splines are the simpliest arithmetic expression that ensure continuity
! and smoothness over an arbitrary data set. When hysteresis in imbibition and
! drainage are not modeled, the model seeks the minimum free energy for a given
! degree of saturation, and conseqently, capillary pressure is necessarially a
! monotonic function of saturation.
!
! Monotonic splines are calculated using the SLATEC Piecewise Cubic Hermite
! Interpolation Polynomial (PCHIP) package, a public domain library.
! PCHIP was optimized for graphing, thus the evaluation function assumes the
! points being evaluated are in an ordered list, and a linear search was used
! to evaluate the entire list. Here, the PCHIP splines are cached as standard
! form polynomials and a binary search is used to locate the the correct spline.
! Additionally, binary search can be performed branchlessly, enabling future
! vectorization, should vector calls to these routines be implemented. For
! clarity, they are left in conditional branched form for scalar calls.

! Caveat, monotonic splines are not recommended for use with Richard's mode.
! Different degrees of saturation can yield identical matrix potentials. Where
! this occurs, capillary pressure and saturation are not uniquely invertible.
! Consequently, Richard's mode has degenerate roots and is ill-defined.
! Furthermore, while cubic polynomials can be efficiently evaluated, the cubic
! is computationally more complex.
!
! References:
!
! Ross, P (1992) "Cubic approximation of hydraulic properties for simulations of
! unsaturated flow", Water Resour. Res. 28(10):2617-2620.
! doi: 10.1029/92WR01310a
!
! Fritsch, FN and Carlson, RE (1980) "Monotone piecewise cubic interpolation",
! SIAM J. Numer. Anal., 17(2):238-246.
! doi: 10.1137/0717021
!
! Fritsch, FN and Butland, J (1984) "A method for constucting local monotone
! piecewise cubic interpolants", SIAM J. Sci. Stat. Comput., 5(2):300-304.
! doi: 10.1137/0905021
!
! **************************************************************************** !

type, private :: knot_type
  PetscReal :: x, y
  type(knot_type), pointer :: next
end type knot_type

! **************************************************************************** !

type, public :: knot_queue_type
! This type exists to support the parser
  private 
    PetscInt :: N = 0
    type(knot_type), pointer :: front
    type(knot_type), pointer :: back
  contains
    procedure, public :: enqueue
    procedure, public :: dequeue
    procedure, public :: depth
    final :: knot_queue_dtor
end type knot_queue_type

! **************************************************************************** !

type, private :: cubic_type
! Polynomial coefficients are stored in a structure for data locality
! x-reference points are stored separately to for efficient binary search
  PetscReal :: y, dy, c2, c3
end type cubic_type
 
! **************************************************************************** !

type, public, extends(sat_func_base_type) :: sf_pchip_type
  private
    PetscInt  :: N ! Number of knots, 1 more than splines
    PetscReal, dimension(:), allocatable :: x ! Reference points
    type(cubic_type), dimension(:), allocatable :: coef ! Coefficients
  contains
    procedure, public :: Init              => SFPCHIPInit
    procedure, public :: CapillaryPressure => SFPCHIPCapillaryPressure
    procedure, public :: Saturation        => SFPCHIPSaturation
    procedure, public :: D2SatDP2          => SFPCHIPD2SatDP2
    procedure, public :: Test              => SFPCHIPTest
    procedure, public :: Label             => SFPCHIPLabel
    final :: SFPCHIPDtor
end type sf_pchip_type

! **************************************************************************** !

type, public, extends(rel_perm_func_base_type) :: rpf_pchip_type
  private
    PetscInt :: N ! Number of knots, 1 more than splines
    PetscReal, dimension(:), allocatable :: x ! Reference points
    type(cubic_type), dimension(:), allocatable :: coef ! Coefficients
  contains
    procedure, public :: Init                 => RPFPCHIPInit
    procedure, public :: RelativePermeability => RPFPCHIPRelativePermeability
    procedure, public :: Test                 => RPFPCHIPTest
    procedure, public :: Label                => RPFPCHIPLabel
    final :: RPFPCHIPDtor
end type rpf_pchip_type

! **************************************************************************** !

private PCHIPCoefficients ! Calculates SF and RPF coefficients via SLATEC/PCHIP
 
private SFPCHIPAllocate
public  SFPCHIPCreate
public  SFPCHIPCtorFunction
public  SFPCHIPCtorQueue

private RPFPCHIPAllocate
public  RPFPCHIPCreate
public  RPFPCHIPCtorFunction
public  RPFPCHIPCtorQueue

contains

! **************************************************************************** !
! Knot queue methods to support the parser
! **************************************************************************** !

subroutine enqueue(this, x, y)
  class(knot_queue_type) :: this
  PetscReal, intent(in) :: x, y

  type(knot_type), pointer :: new

  allocate(new)
  new%x = x
  new%y = y
  nullify(new%next) ! Null pointer at back end

  if (this%N == 0) then ! First node
    this%front => new
    this%back  => new
  else ! Enqueue at back
    this%back%next => new
    this%back => new
  end if

  this%N = this%N + 1

end subroutine enqueue

! **************************************************************************** !

subroutine dequeue(this, x, y)
  class(knot_queue_type) :: this
  PetscReal, intent(out) :: x, y

  type(knot_type), pointer :: old

  if (this%N == 0) return

  this%N = this%N - 1

  old => this%front
  this%front => old%next

  x = old%x
  y = old%y
  deallocate(old)

end subroutine dequeue

! **************************************************************************** !

function depth(this)
  class(knot_queue_type) :: this
  PetscInt :: depth

  depth = this%N
end function

! **************************************************************************** !

subroutine knot_queue_dtor(this)
  implicit none
  type(knot_queue_type) :: this

  type(knot_type), pointer :: old

  do while (this%N > 0)
    this%N = this%N - 1

    old => this%front
    this%front => old%next

    deallocate(old)
  end do

end subroutine

! **************************************************************************** !
! SLATEC/PCHIP/PCHIM Wrapper
! **************************************************************************** !

subroutine PCHIPCoefficients(N, x, y, dy, c2, c3)
  PetscInt, intent(in)   :: N
  PetscReal, intent(in)  :: x(N), y(N)
  PetscReal, intent(out) :: dy(N-1), c2(N-1), c3(N-1)
  PetscInt :: I
  PetscReal :: h, delta, del1, del2

  ! Use SLATEC/PCHIP/PCHIM to calculate first derivatives for monotonic splines
  ! I is part of the SLATEC error handling, but it is not used here.
  I = 0
  call PCHIM(N, x, y, dy, 1, I)

  ! 2nd and 3rd derivatives are defined by adjacent Hermite polynomials.
  ! Here, they are cached in standard polynomial form.
  ! Note, N knots generates N-1 splines
  do I = 1, N-1
    h     =  x(i+1)  - x(i)
    delta = (y(i+1)  - y(i))/h
    del1  = (dy(i)   - delta)/h
    del2  = (dy(i+1) - delta)/h

    c2(i) = -(del1+del1+del2)
    c3(i) =  (del1     +del2)/h
  end do
end subroutine

! **************************************************************************** !
! Saturation Function Spline Methods
! **************************************************************************** !

function SFPCHIPCreate() result (new)
  class(sf_pchip_type), pointer :: new
! Create a placeholder object of the PCHIP type for the parser
  allocate(new)
  call new%Init()
end function

! ************************************************************************** !

subroutine SFPCHIPInit(this)
  class(sf_pchip_type) :: this
! Unused objects in base class that must be nullified
  nullify(this%sat_poly)
  nullify(this%pres_poly)
! Unitialized values, which ought to be set in ctor
  this%Sr = UNINITIALIZED_DOUBLE
  this%pcmax = UNINITIALIZED_DOUBLE
! Common values regardless of ctor
  this%analytical_derivative_available = PETSC_TRUE
  this%calc_int_tension = PETSC_FALSE
  this%calc_vapor_pressure = PETSC_FALSE
end subroutine SFPCHIPInit

! **************************************************************************** !

function SFPCHIPAllocate(N) result (new)
  implicit none
  PetscInt, intent(in) :: N
  class(sf_pchip_type), pointer :: new

  nullify(new)

  ! Return null as no interpolation is possible with less than 2 knots
  if (N < 2) return

  ! Allocate spline object
  allocate(new)
  if (.not. associated(new)) return ! Return null as allocation failed

  ! Default initialization
  call new%Init()
  new%N = N

  ! Allocate subordinate dynamic objects
  allocate(new%x(N))
  allocate(new%coef(N)) ! Allocating for N knots/N-1 splines
  if (.not. allocated(new%x) .or. .not. allocated(new%coef)) then
    ! Note, final method will deallocate if only one is allocated
    deallocate(new)
    nullify(new)
  end if ! Return null if any subordinate allocation failed

end function

! **************************************************************************** !

subroutine SFPCHIPDtor(this)
  implicit none
  type(sf_pchip_type) :: this
  ! Deallocate subordinate dynamic objects
  ! These will not have been allocated if "create" is used instead of a ctor
  if (allocated(this%x)) deallocate(this%x)
  if (allocated(this%coef)) deallocate(this%coef)
end subroutine

! **************************************************************************** !

function SFPCHIPCtorFunction(N, sf_analytic) result (new)
  implicit none
  PetscInt, intent(in) :: N
  class(sat_func_base_type), intent(in) :: sf_analytic
  class(sf_pchip_type), pointer :: new

  PetscInt  :: I
  type(option_type) :: option ! Placeholder

  new => SFPCHIPAllocate(N+1) ! Allocate space for N+1 knots
  if (.not. associated(new)) return

  ! Copy attributes of saturation function to be approximated
  new%Sr    = sf_analytic%Sr
  new%Pcmax = sf_analytic%Pcmax

  ! Generate N+1 evenly spaced knots starting at 0 and ending at 1
  do I = 1, N+1
    new%x(i) = dble(I-1)/dble(N)
    call sf_analytic%CapillaryPressure(new%x(i), new%coef(i)%y, new%coef(i)%dy, option)
  end do
 
! Calculate the PCHIP splines
  call PCHIPCoefficients(N+1, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function SFPCHIPCtorFunction

! **************************************************************************** !

function SFPCHIPCtorQueue(queue) result (new)
  implicit none
  type(knot_queue_type), intent(inout) :: queue
  class(sf_pchip_type), pointer :: new

  PetscInt :: i
  PetscReal :: x, y

  print *, "Entering queue ctor with depth ", queue%depth()

  new => SFPCHIPAllocate(queue%depth())
  if (.not. associated(new)) return

! Pack queue into the arrays, assuming queue was in order of increasing saturation
! A sort method could go here if desired
  do i = 1, new%N
    call queue%dequeue(new%x(i), new%coef(i)%y)
  end do

! Set base class attributes with limiting knots
  new%Sr    = new%x(1)
  new%Pcmax = new%coef(1)%y

! Calculate the PCHIP splines
  call PCHIPCoefficients(new%N, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function SFPCHIPCtorQueue

! **************************************************************************** !

function SFPCHIPLabel(this) result (Label)
  class(sf_pchip_type) :: this
  character(len=MAXSTRINGLENGTH) :: Label
! Class label to simplify parser
  Label = 'PCHIP Splines'
end function

! **************************************************************************** !

subroutine SFPCHIPCapillaryPressure(this, liquid_saturation, &
                                   capillary_pressure, dPc_dSatl, option)
  implicit none
  class(sf_pchip_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option

  PetscInt :: i, j, k
  PetscReal :: x

! Truncate saturation to be within bounds
  x = min(max(liquid_saturation,this%x(1)),this%x(this%N))

! Binary search for correct polynomial
  i = 1
  j = this%n
  do while (j - i > 1)
    k = (i + j)/2
    if (this%x(k) > x) then
      j = k
    else
      i = k
    end if
  end do

! Polynomials are defined with respect to nearest knot for precision
  x  = x - this%x(i)

! Horner's method leverages fused multiply-add
  capillary_pressure    =    this%coef(i)%y + x* &
                        (   this%coef(i)%dy + x* &
                        (   this%coef(i)%c2 + x*this%coef(i)%c3))
  dPc_dSatl             =   this%coef(i)%dy + x* &
                        ( 2*this%coef(i)%c2 + 3*x*this%coef(i)%c3)

end subroutine SFPCHIPCapillaryPressure

! **************************************************************************** !

subroutine SFPCHIPSaturation(this, capillary_pressure, &
                            liquid_saturation, dsat_dpres, option)
  implicit none
  class(sf_pchip_type) :: this
  PetscReal, intent(in)  :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  PetscReal, parameter :: pi = 4*atan(1d0)
  PetscInt :: i, j, k
  PetscReal :: x, dx, y
  PetscReal :: a, b, c
  PetscReal :: q, r, s, t

  ! Truncate capillary pressure to be within bounds
  y = min(max(capillary_pressure, this%coef(this%N)%y), this%Pcmax)

  ! Binary search to find correct polynomial
  i = 1
  j = this%n
  do while (j - i > 1)
    k = (i + j) / 2
    if (this%coef(k)%y < y) then
     j = k
    else
     i = k
    end if
  end do

  ! Find root of polynomial of up to 3rd order
  if (this%coef(i)%c3 /= 0d0) then ! Cubic
    a = this%coef(i)%c2       / this%coef(i)%c3
    b = this%coef(i)%dy       / this%coef(i)%c3
    c = (this%coef(i)%y - y)  / this%coef(i)%c3

    q = (    a**2 - 3d0*  b         ) /  9d0
    r = (2d0*a**3 - 9d0*a*b + 27d0*c) / 54d0
    s = r**2 - q**3

    if (s < 0d0) then
      ! Three real roots, but only one in the domain
      t = acos(r/sqrt(q**3)) ! Angle in radians
      ! 1st root
      x = -2d0*sqrt(q)*cos(t/3d0) - a/3d0
      if (x < 0d0 .or. x > (this%x(i+1) - this%x(i))) then
        ! 2nd root
        x = -2d0*sqrt(q)*cos((t + 2d0*pi)/3d0) - a/3d0
        if (x < 0d0 .or. x > (this%x(i+1) - this%x(i))) then
          ! 3rd root
          x = -2d0*sqrt(q)*cos((t - 2d0*pi)/3d0) - a/3d0
        end if
      end if
    else ! One real root, complex roots are ignored
      t = sign( (abs(r)+sqrt(s))**(1d0/3d0), -r)
      if (t /= 0d0) then
        x = t + q/t - a/3d0
      else
        x = -a/3d0
      end if
    end if

    dx =   this%coef(i)%dy + x* ( 2*this%coef(i)%c2 + 3*x*this%coef(i)%c3)
  else if (this%coef(i)%c2 /= 0d0) then ! Quadratic
    a = this%coef(i)%c2
    b = this%coef(i)%dy
    c = this%coef(i)%y - y

    q = (b + sign(sqrt(b**2 - 4d0*a*c), b)) / (-2d0)

    x = q/a ! Check if this root is within spline domain
    if (x < 0d0 .or. x > (this%x(i+1) - this%x(i))) then
      x = c/q ! If not, the other must be
    end if

    dx = this%coef(i)%dy + 2*x*this%coef(i)%c2
  else if (this%coef(i)%dy /= 0d0) then ! Linear
    x = (y - this%coef(i)%y) / this%coef(i)%dy
    dx = this%coef(i)%dy
  else ! Degenerate
    x = 0d0
    dx = tiny(this%coef(i)%dy)
  end if

  liquid_saturation = x + this%x(i) ! Linear translation to saturation space
  dsat_dpres = 1d0 / dx ! Invert 1st derivative

end subroutine SFPCHIPSaturation

! **************************************************************************** !a

subroutine SFPCHIPD2SatDP2(this,Pc, d2s_dp2, option)
  implicit none
  class(sf_pchip_type) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option

  PetscInt :: i, j, k
  PetscReal :: Sw, dSw_dPc

  ! Using inverse/first derivative rather than repeating cubic formula etc.
  call this%Saturation(Pc, Sw, dSw_dPc, option)

  ! Repeat binary search to get index for 2nd derivative
  i = 1
  j = this%n
  do while (j - 1 > 1)
    k = (i + j) /2
    if (this%coef(i)%y > Pc) then
     j = k
    else
     i = k
    end if
  end do

  d2s_dp2 = -(2*this%coef(i)%c2 + 6*Sw*this%coef(i)%c3) / dSw_dPc**3

end subroutine SFPCHIPD2SatDP2

! **************************************************************************** !

subroutine SFPCHIPTest(this,cc_name,option)
  use Option_module
  use Material_Aux_module
  use PFLOTRAN_Constants_module
  implicit none
  class(sf_pchip_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: i

! Write knots and spline derivatives to file
  write(string,*) cc_name
  string = trim(cc_name) // '_pc_knots.dat'
  open(unit=87,file=string)
  write(87,*) '#Index, Sw, Pc, dPc/dSw'

  do i = 1, this%N
    write(87,*) i, this%x(i), this%coef(i)%y, this%coef(i)%dy
  end do

! Also call base test
  call SFBaseTest(this, cc_name, option)

end subroutine SFPCHIPTest

! **************************************************************************** !
! Relative Permeability PCHIP Methods
! **************************************************************************** !

function RPFPCHIPCreate() result (new)
  class(rpf_pchip_type), pointer :: new
  allocate(new)
  call new%Init()
end function

! **************************************************************************** !
subroutine RPFPCHIPInit(this)
  implicit none
  class(rpf_pchip_type) :: this
! Unused objects in base class that must be nullified
  nullify(this%poly)
! Unitialized values, which ought to be set in ctor
  this%Sr = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
! Common values regardless of ctor
  this%analytical_derivative_available = PETSC_TRUE
end subroutine RPFPCHIPInit

! **************************************************************************** !

function RPFPCHIPAllocate(N) result (new)
  ! Perform allocation for both data or function defined splines
  implicit none
  PetscInt, intent(in) :: N
  class(rpf_pchip_type), pointer :: new

  nullify(new)

  ! Return null as no interpolation is possible with less than 2 knots
  if (N < 2) return

  ! Allocate spline object
  allocate(new)
  if (.not. associated(new)) return ! Return null as allocation failed
  call new%Init()
  new%N = N

  ! Allocate subordinate dynamic objects
  allocate(new%x(N))
  allocate(new%coef(N)) ! Allocating for N knots/N-1 splines)
  if (.not. allocated(new%x) .or. .not. allocated(new%coef)) then
    ! Note, final method will deallocate if only one is allocated
    deallocate(new)
    nullify(new)
  end if ! Return null if any subordinate allocation failed
end function

! **************************************************************************** !

subroutine RPFPCHIPDtor(this)
  type(rpf_pchip_type) :: this
  ! Deallocate subordinate dynamic objects
  ! These will not have been allocated if "create" is used instead of a ctor
  if (allocated(this%x)) deallocate(this%x)
  if (allocated(this%coef)) deallocate(this%coef)
end subroutine RPFPCHIPDtor

! **************************************************************************** !

function RPFPCHIPCtorFunction(N, rpf_analytic) result (new)
  implicit none

  class(rpf_pchip_type), pointer :: new
  PetscInt, intent(in) :: N
  class(rel_perm_func_base_type), intent(in) :: rpf_analytic

  PetscInt  :: I
  type(option_type) :: option

  new => RPFPCHIPAllocate(N+1) ! Allocate space for N+1 knots
  if (.not. associated(new)) return

! Copy attributes of relative permeability function to be approximated
  new%Sr  = rpf_analytic%Sr
  new%Srg = rpf_analytic%Srg
  ! Caveat, liquid relative perm functions may have an "uninitialized" gas residual
  if (rpf_analytic%Srg == UNINITIALIZED_DOUBLE) new%Srg = 0d0

! Generate N+1 evenly spaced knots starting at Sr and ending at 1-Srg
  do I = 1, N
    new%x(i) = (1d0 - new%Srg - new%Sr) * dble(I-1)/dble(N) + new%Sr
    call rpf_analytic%RelativePermeability(new%x(I), new%coef(I)%y, new%coef(I)%dy, option)
  end do
! Ensure the gas residual endpoint is precise
  new%x(N+1) = 1d0 - new%Srg
  call rpf_analytic%RelativePermeability(new%x(N+1), new%coef(N+1)%y, new%coef(I+1)%dy, option)

! Calculate the PCHIP splines
  call PCHIPCoefficients(N+1, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function RPFPCHIPCtorFunction

! **************************************************************************** !

function RPFPCHIPCtorQueue(queue) result (new)
  implicit none
  class(knot_queue_type), intent(inout) :: queue
  class(rpf_pchip_type), pointer :: new

  PetscInt :: i, N

  N = queue%depth()

  new => RPFPCHIPAllocate(N)
  if (.not. associated(new)) return

! Pack queue into the arrays, assuming queue is in order of increasing saturation
! A sort method could go here if desired
  do i = 1, N
    call queue%dequeue(new%x(i), new%coef(i)%y)
  end do

! Set base class attributes with limiting knots
  new%Sr  = new%x(1)
  new%Srg = 1d0 - new%x(N)

! Calculate the PCHIP splines
  call PCHIPCoefficients(N, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function RPFPCHIPCtorQueue

! **************************************************************************** !

function RPFPCHIPLabel(this) result (Label)
  class(rpf_pchip_type) :: this
  character(len=MAXSTRINGLENGTH) :: Label
! Class label to simplify parser
  Label = 'PCHIP Splines'
end function RPFPCHIPLabel

! **************************************************************************** !

subroutine RPFPCHIPRelativePermeability(this, liquid_saturation, &
                                   relative_permeability, dkr_sat, option)
  implicit none
  class(rpf_pchip_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: relative_permeability, dkr_sat
  type(option_type), intent(inout) :: option

  PetscInt :: i, j, k
  PetscReal :: x

! Truncate saturation to be within bounds
  x = min(max(liquid_saturation,this%x(1)), this%x(this%N))

! Binary search for the correct polynomial
  i = 1
  j = this%n
  do while (j - i > 1)
    k = (i + j)/2
    if (this%x(k) > x) then
      j = k
    else
      i = k
    end if
  end do

! Polynomials are defined with respoect to the nearest knto for precision
  x  = x - this%x(i)

! Horner's method leverages fused multiply-add
  relative_permeability =    this%coef(i)%y + x* &
                        (   this%coef(i)%dy + x* &
                        (   this%coef(i)%c2 + x*this%coef(i)%c3))
  dkr_sat               =   this%coef(i)%dy + x* &
                        ( 2*this%coef(i)%c2 + 3*x*this%coef(i)%c3)

end subroutine RPFPCHIPRelativePermeability

! **************************************************************************** !

subroutine RPFPCHIPTest(this,cc_name,phase,option)
  use Option_module
  use Material_Aux_module
  use PFLOTRAN_Constants_module
  implicit none
  class(rpf_pchip_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option
  character(len=MAXWORDLENGTH) :: phase
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: i

! Write knots and spline derivatives to file
  write(string,*) cc_name
  string = trim(cc_name) // '_' //  trim(phase) // '_Kr_knots.dat'
  open(unit=87,file=string)
  write(87,*) '#Index, Sa, Kra, dKra/dSa'
  do i = 1, this%N
    write(87,*) i, this%x(i), this%coef(i)%y, this%coef(i)%dy
  end do
  close(87)

! Also call base test
  call RPFBaseTest(this, cc_name, phase, option)
    
end subroutine RPFPCHIPTest

! **************************************************************************** !

end module
