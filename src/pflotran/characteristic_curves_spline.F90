module Characteristic_Curves_spline_module
#include "petsc/finclude/petscsys.h"

use petscsys ! Necessary for PETSC_TRUE / PETSC_FALSE
use Characteristic_Curves_Base_module   ! Needed to define base type
use Option_module ! Needed for Verify and unused arguments in Pc and Sl
use spline_module

implicit none
private

! **************************************************************************** !
!
! References:
!
! Ross, P (1992) "Cubic approximation of hydraulic properties for simulations of
! unsaturated flow", Water Resource Research 28(10):2617-2620.
! doi: 10.1029/92WR01310a
!
! TODO Inline the function in spline_module to use data structures, avoiding
! overhead. Also, a hash function would be better than bisection.
!
! Numerical Recipes in Fortran 77: The Art of Scientific Computing, 2nd ed. 
!
! **************************************************************************** !
!
! sat_func_base_type        External base type, quasi-abstract
! |
! |-->sf_spline_type        Internal spline type
!
! **************************************************************************** !

public :: SFSplineCreate & ! "create" needed for CharacteristicCurvesRead
        , SFSplineCtor     ! Constructor

type, private :: sf_spline
  ! TODO replace primitive arrays with an array of structures
  ! Idea here is to keep values used together, together
  ! If the arrays are large, x, y, and dy2 could be on different pages of memory
  ! resulting in page faults
  PetscReal :: x, y, dy2
end type

! **************************************************************************** !
type, public, extends(sat_func_base_type) :: sf_spline_type
  private
    PetscInt :: N  ! Number of knots
    PetscReal :: h ! Saturation interval width
    Real*8, dimension(:), allocatable :: x, y, dy2
  contains
! Definition of base type methods
    procedure, public  :: Init                  => SFSplineInit
    procedure, public  :: CapillaryPressure     => SFSplineCapillaryPressure
    procedure, public  :: Saturation            => SFSplineSaturation
    procedure, public  :: D2SatDP2              => SFSplineD2SatDP2
    procedure, public  :: Test                  => SFSplineTest
    procedure, public  :: Verify                => SFSplineVerify
end type

contains

! **************************************************************************** !
! Spline Methods
! **************************************************************************** !

function SFSplineCreate() result (new)
  implicit none
  class(sf_spline_type), pointer :: new

  allocate(new)
end function SFSplineCreate

subroutine SFSplineInit(this)
  implicit none
  class(sf_spline_type) :: this
end subroutine SFSplineInit

subroutine SFSplineVerify(this,name,option)
  use PFLOTRAN_Constants_module
  use Option_module
  implicit none
  class(sf_spline_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

end subroutine SFSplineVerify

! **************************************************************************** !

function SFSplineCtor(sf_analytic, N) result (new)
  implicit none
  class(sf_spline_type), pointer :: new
  class(sat_func_base_type), intent(in) :: sf_analytic
  PetscInt, intent(in) :: N
  PetscInt :: I
  PetscReal :: buffer
  type(option_type) :: option

  ! Passing an array of knots, we seek the cubic polynomial splines to match
  ! the knots exactly. 
  allocate(new)
  if (.not. associated(new)) return
  allocate(new%x(N))
  allocate(new%y(N))
  allocate(new%dy2(N))

  new%N = N
  new%h = 1d0/dble(N-1)
  do I = 1, new%N
    new%x(I) = dble(I-1)/dble(N-1)
    call sf_analytic%CapillaryPressure(new%x(I), new%y(I), buffer, option)
  end do

! Calculate 2nd derivitives, as per:
!     press, w.h., b.p. flannery, s.a. teukolsky, and w.t. vetterling.
!     1986.  numerical recipes, the art of scientific computing,
!     cambridge university press, cambridge.  pp. 86-89.

  call spline(new%x, new%y, new%N, new%dy2)

  new%analytical_derivative_available = PETSC_TRUE

end function SFSplineCtor

! **************************************************************************** !

subroutine SFSplineCapillaryPressure(this, liquid_saturation, &
                                   capillary_pressure, dPc_dSatl, option)
  implicit none

  class(sf_spline_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option
  PetscInt :: klo,khi
  PetscReal :: h,a,b,x,y, dy

! Evalute the cubic spline at index I
  
!     cubic spline interpolation.

!     press, w.h., b.p. flannery, s.a. teukolsky, and w.t. vetterling.
!     1986.  numerical recipes, the art of scientific computing,
!     cambridge university press, cambridge.  pp. 86-89.
  x = liquid_saturation

! Table index is found using a hash function instead of bisection
  klo = ceiling(x/this%h)
  if (klo <= 0) klo = 1
  khi = klo+1

  a = (this%x(khi)-x)/this%h
  b = (x-this%x(klo))/this%h
  c = (a**3-a)*this%h**2/6d0
  d = (b**3-b)*this%h**2/6d0

  y = a*this%y(klo) + b*this%y(khi) + c*this%dy2(klo) + d*this%dy2(khi)

  dy = (this%y(khi)-this%y(klo))/this%h + this%h/6d0 * &
     (-(3d0*a**2-1)*this%dy2(klo) + (3d0*b**2-1)*this%dy2(khi))

  capillary_pressure = y
  dPc_dSatl = dy

end subroutine SFSplineCapillaryPressure

! **************************************************************************** !

subroutine SFSplineSaturation(this, capillary_pressure, &
                            liquid_saturation, dsat_dpres, option)
  implicit none

  class(sf_spline_type) :: this
  PetscReal, intent(in)  :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  ! TODO
  liquid_saturation = 0d0
  dsat_dpres = 0d0
end subroutine SFSplineSaturation

! **************************************************************************** !

subroutine SFSplineD2SatDP2(this,Pc, d2s_dp2, option)

  implicit none

  class(sf_spline_type) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option

  ! TODO
  d2s_dp2 = 0d0
end subroutine SFSplineD2SatDP2

subroutine SFSplineTest(this,cc_name,option)

  use Option_module
  use Material_Aux_module
  use PFLOTRAN_Constants_module
  implicit none
  class(sf_spline_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt, parameter :: num_values = 1000
  PetscReal :: Pc, Sw, dPc_dSw
  PetscInt :: i

  write(string,*) cc_name
  string = trim(cc_name) // '_pc_from_sat.dat'
  open(unit=86,file=string)
  write(86,*) '#Sw,       Pc,       dPc/dSw'

 ! calculate capillary pressure as a function of saturation
  do i = 0, num_values
    Sw = dble(i)/dble(num_values)
    call this%CapillaryPressure(Sw, Pc, dPc_dSw, option)
    write(86,'(4es14.6)') Sw, Pc, dPc_dSw
  enddo
  close(86)

  write(string,*) cc_name
  string = trim(cc_name) // '_splines.dat'
  open(unit=87,file=string)
  write(87,*) '# Index, x, y, dy2'
  do i = 1, 101
    write(87, *) i, this%x(i), this%y(i), this%dy2(i)
  end do
  close(87)
  


end subroutine SFSplineTest

end module

