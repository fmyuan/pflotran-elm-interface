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

public SFSplineCtor     ! Constructor

! **************************************************************************** !
type, public, extends(sat_func_base_type) :: sf_spline_type
  private
    PetscInt  :: N ! Number of knots
    PetscReal :: h ! Saturation interval width
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

subroutine SFSplineInit(this)
  implicit none
  class(sf_spline_type) :: this

  nullify(this%sat_poly)
  nullify(this%pres_poly)
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
  PetscReal :: u(N)
  PetscReal :: sig, p, qn, un
  type(option_type) :: option

  ! Passing an array of knots, we seek the cubic polynomial splines to match
  ! the knots exactly. 
  allocate(new)
  if (.not. associated(new)) return
  allocate(new%spline(N))

  new%N = N
  new%h = 1d0/dble(N-1)
  do I = 1, new%N
    new%spline(I)%x = dble(I-1)/dble(N-1)
    call sf_analytic%CapillaryPressure(new%spline(I)%x, new%spline(I)%y, buffer, option)
  end do

! Calculate 2nd derivitives, as per:
!     press, w.h., b.p. flannery, s.a. teukolsky, and w.t. vetterling.
!     1986.  numerical recipes, the art of scientific computing,
!     cambridge university press, cambridge.  pp. 86-89.

  new%spline(1)%dy2 = 0d0
  u(1) = 0d0
  do i = 2,n-1
   sig = (new%spline(I)%x-new%spline(I-1)%x)/(new%spline(I+1)%x-new%spline(I-1)%x)
   p = sig*new%spline(I-1)%dy2+2d0
   new%spline(I)%dy2 = (sig-1d0)/p

   u(i) = (6d0*((new%spline(I+1)%y-new%spline(I)%y)/(new%spline(I+1)%x-new%spline(I)%x) - &
   (new%spline(I)%y-new%spline(I-1)%y)/(new%spline(I)%x-new%spline(I-1)%x))/ &
   (new%spline(I+1)%x-new%spline(I-1)%x) - sig*u(i-1))/p
  enddo
  qn = 0d0
  un = 0d0
  new%spline(N)%dy2 = (un-qn*u(n-1))/(qn*new%spline(N-1)%dy2+1d0)

  do i = n-1,1,-1
    new%spline(I)%dy2 = new%spline(I)%dy2*new%spline(I+1)%dy2+u(i)
  enddo

  new%analytical_derivative_available = PETSC_TRUE
  nullify(new%sat_poly)
  nullify(new%pres_poly)

end function SFSplineCtor

! **************************************************************************** !

subroutine SFSPlinePcV(this, N, Sw, Pc)
  implicit none
  class(sf_spline_type)  :: this
  PetscInt,  intent(in)  :: N
  PetscReal, intent(in)  :: Sw(N)
  PetscReal, intent(out) :: Pc(N)
  PetscInt :: klo(N), khi(N)
  PetscReal :: a(N), b(N), c(N), d(N)

  ! Possibly, like BLAS-LAPACK, we may need to cut the local variables to some register width.
  ! E.g. set to be 4x wide, do the modulo 4 as scalar, then do the rest 4 at a time
  ! The problem is Sw and Pc may have no logical packing.

  klo = max(ceiling(Sw/this%h),1)
  khi = klo + 1

  a = (this%spline(khi)%x-Sw)/this%h
  b = (Sw-this%spline(klo)%x)/this%h
  c = (a**3-a)*this%h*this%h/6d0
  d = (b**3-b)*this%h*this%h/6d0

  Pc = a*this%spline(klo)%y + b*this%spline(khi)%y + c*this%spline(klo)%dy2 + d*this%spline(khi)%dy2

end subroutine

subroutine SFSplineCapillaryPressure(this, liquid_saturation, &
                                   capillary_pressure, dPc_dSatl, option)
  implicit none

  class(sf_spline_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option
  PetscInt :: klo,khi
  PetscReal :: h,a,b,c,d,x,y,dy

! Evalute the cubic spline at index I
  
!     cubic spline interpolation.

!     press, w.h., b.p. flannery, s.a. teukolsky, and w.t. vetterling.
!     1986.  numerical recipes, the art of scientific computing,
!     cambridge university press, cambridge.  pp. 86-89.
  x = liquid_saturation

! Table index is found using a hash function instead of bisection
  klo = max(ceiling(x/this%h),1)
  khi = klo+1

  a = (this%spline(khi)%x-x)/this%h
  b = (x-this%spline(klo)%x)/this%h
  c = (a**3-a)*this%h*this%h/6d0
  d = (b**3-b)*this%h*this%h/6d0

  y = a*this%spline(klo)%y + b*this%spline(khi)%y + c*this%spline(klo)%dy2 + d*this%spline(khi)%dy2

  dy = (this%spline(khi)%y-this%spline(klo)%y)/this%h + this%h/6d0 * &
     (-(3d0*a*a-1)*this%spline(klo)%dy2 + (3d0*b*b-1)*this%spline(khi)%dy2)

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
    write(87, *) i, this%spline(i)
  end do
  close(87)
  
end subroutine SFSplineTest

end module

