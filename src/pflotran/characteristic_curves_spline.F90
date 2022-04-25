module Characteristic_Curves_spline_module
#include "petsc/finclude/petscsys.h"

use petscsys ! Necessary for PETSC_TRUE / PETSC_FALSE
use Characteristic_Curves_Base_module   ! Needed to define base type
use Option_module ! Needed for Verify and unused arguments in Pc and Sl
use spline_module
use PFLOTRAN_constants_module ! Needed for UNINITIALIZED_DOUBLE

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
! **************************************************************************** !
!
! sat_func_base_type        External base type, quasi-abstract
! |
! |-->sf_spline_type        Hashed cubic spline type
!
! **************************************************************************** !

public SFSplineCtor

! **************************************************************************** !
type, public, extends(sat_func_base_type) :: sf_spline_type
  private
    PetscInt  :: N ! Number of splines
    PetscReal :: h ! Saturation interval width
  contains
    procedure, public :: Init              => SFSplineInit
    procedure, public :: CapillaryPressure => SFSplineCapillaryPressure
!   Inverse of cubic splines is far more complex. Better to make separate
!   set iff they are needed.
!   procedure, public :: Saturation        => SFSplineSaturation
!   procedure, public :: D2SatDP2          => SFSplineD2SatDP2
    procedure, public :: Test              => SFSplineTest
end type

! **************************************************************************** !

public RPFSplineCtor

! **************************************************************************** !

type, public, extends(rel_perm_func_base_type) :: rpf_spline_type
  private
    PetscInt  :: N ! Number of splines
    PetscReal :: h ! Saturation interval width
  contains
    procedure, public :: Init                 => RPFSplineInit
    procedure, public :: RelativePermeability => RPFSplineRelativePermeability
end type

contains

! **************************************************************************** !
! Saturation Function Spline Methods
! **************************************************************************** !

subroutine SFSplineInit(this)
  implicit none
  class(sf_spline_type) :: this
  ! This method is intentionally left blank.
end subroutine SFSplineInit

! **************************************************************************** !

function SFSplineCtor(sf_analytic, N) result (new)
  implicit none
  class(sf_spline_type), pointer :: new
  class(sat_func_base_type), intent(in) :: sf_analytic
  PetscInt, intent(in) :: N
  PetscInt  :: I
  PetscReal :: buffer
  PetscReal :: x(N+1), y(N+1), dy2(N+1)
  PetscReal :: A, B, C, D
  type(option_type) :: option

  allocate(new)
  if (.not. associated(new)) return
  allocate(new%spline(N))

  if (.not. allocated(new%spline)) then
    deallocate(new)
    nullify(new)
    return
  end if

  ! Generate N knots at equal linear spacing h TODO could be truncation error
  new%N = N
  new%h = 1d0/dble(N-1)

  ! Find N+1 knots for N splines
  do I = 1, new%N+1
    x(I) = dble(I)/dble(N)
    call sf_analytic%CapillaryPressure(x(I), y(I), buffer, option)
  end do
  ! Calculate 2nd derivatives
  call spline(x, y, N+1, dy2)

! Store splines in standard polynomial form - more memory but less floating point ops
  do i = 1,n-1
    A =   - dy2(i)  /(6d0*new%h)
    A = A + dy2(i+1)/(6d0*new%h)

    B =     x(i+1)*dy2(i)  /(2d0*new%h)
    B = B - x(i)  *dy2(i+1)/(2d0*new%h)

    C =   - x(i+1)**2*dy2(i)  /(2d0*new%h) + new%h*dy2(i)  /6d0 - y(i)  /new%h
    C = C + x(i)  **2*dy2(i+1)/(2d0*new%h) - new%h*dy2(i+1)/6d0 + y(i+1)/new%h

    D =     x(i+1)**3*dy2(i  )/(6d0*new%h) - x(i+1)*new%h*dy2(i)  /6d0 + x(i+1)*y(i)/new%h
    D = D - x(i)  **3*dy2(i+1)/(6d0*new%h) + x(i)  *new%h*dy2(i+1)/6d0 - x(i)*y(i+1)/new%h

    new%spline(i)%A = A
    new%spline(i)%B = B
    new%spline(i)%C = C
    new%spline(i)%D = D
  end do

  new%analytical_derivative_available = PETSC_TRUE

  ! Unused pointers and attributes in base class
  nullify(new%sat_poly)
  nullify(new%pres_poly)
  new%Sr = 0d0
  call new%CapillaryPressure(0d0, new%Pcmax, buffer, option) ! Popuplate Pcmax just in case
  new%calc_int_tension = PETSC_FALSE ! Default to false unless overriden elsewhere

end function SFSplineCtor

! **************************************************************************** !

subroutine SFSplineCapillaryPressure(this, liquid_saturation, &
                                   capillary_pressure, dPc_dSatl, option)
  implicit none

  class(sf_spline_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option
  PetscInt :: i ! This really should not be a huge integer
  PetscReal :: x, y, dy

! Define alias because this is problematically long
  x = liquid_saturation

! Natural table index is found using a hash function
  i = max(ceiling(x/this%h),1) ! TODO determine if min/max bounds checking are necessary
  
  y  = ((this%spline(i)%A*x + this%spline(i)%B)*x + this%spline(i)%C)*x + this%spline(i)%D
  dy = (3d0*this%spline(i)%A*x + 2d0*this%spline(i)%B) * x + this%spline(i)%C

! Again, aliases because dummy variable names are absurdly long. Compiler should inline
  capillary_pressure = y
  dPc_dSatl = dy

end subroutine SFSplineCapillaryPressure

! **************************************************************************** !

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

! **************************************************************************** !
! Relative Permeability Spline Methods
! **************************************************************************** !

subroutine RPFSplineInit(this)
  implicit none
  class(rpf_spline_type) :: this
  ! This method is intentionally left blank.
end subroutine RPFSplineInit

! **************************************************************************** !

function RPFSplineCtor(rpf_analytic, N) result (new)
  implicit none
  class(rpf_spline_type), pointer :: new
  class(rel_perm_func_base_type), intent(in) :: rpf_analytic
  PetscInt, intent(in) :: N
  PetscInt  :: I
  PetscReal :: buffer
  PetscReal :: x(N+1), y(N+1), dy(2), d2y(N+1)
  PetscReal :: A, B, C, D
  type(option_type) :: option

  allocate(new)
  if (.not. associated(new)) return

  allocate(new%spline(N+2))
  if (.not. allocated(new%spline)) then
    deallocate(new)
    nullify(new)
    return
  end if

  ! Initialize public attributes
  nullify(new%poly)
  new%Sr = rpf_analytic%Sr
  new%Srg = rpf_analytic%Srg
  new%analytical_derivative_available = PETSC_TRUE
  ! Liquid relative perm functions may have an "uninitialized" gas residual
  if (rpf_analytic%Srg == UNINITIALIZED_DOUBLE) new%Srg = 0d0

  ! Calculate private attributes
  new%N = N
  new%h = (1d0 - new%Sr - new%Srg)/dble(N) ! Width of splines between the residuals

  ! N+1 knots (x, y) for N internal splines
  x(1) = new%Sr
  call rpf_analytic%RelativePermeability(x(1), y(1), dy(1), option)
  do I = 2, new%N
    x(I) = new%Sr + new%h*dble(I-1)
    call rpf_analytic%RelativePermeability(x(I), y(I), buffer, option)
  end do
  x(N+1) = 1d0 - new%Srg
  call rpf_analytic%RelativePermeability(x(N+1), y(N+1), dy(2), option)

  ! Calculate 2nd derivatives
  call spline(x, y, N+1, d2y)
! 1st derivatives are incorrect at the end-points in some analytical functions
! For now, assume "natural" splines with 0 2nd derivatives at end
! call RPFspline(x, y, N+1, dy(1), dy(2), d2y)

  ! Store splines in standard polynomial form for speed
  ! External "spline" below liquid residual
  new%spline(1)%A = 0d0
  new%spline(1)%B = 0d0
  new%spline(1)%C = 0d0
  new%spline(1)%D = y(1)
  ! Internal splines
  do i = 2,n+1
    A =   - d2y(i)  /(6d0*new%h)
    A = A + d2y(i+1)/(6d0*new%h)

    B =     x(i+1)*d2y(i)  /(2d0*new%h)
    B = B - x(i)  *d2y(i+1)/(2d0*new%h)

    C =   - x(i+1)**2*d2y(i)  /(2d0*new%h) + new%h*d2y(i)  /6d0 - y(i)  /new%h
    C = C + x(i)  **2*d2y(i+1)/(2d0*new%h) - new%h*d2y(i+1)/6d0 + y(i+1)/new%h

    D =     x(i+1)**3*d2y(i  )/(6d0*new%h) - x(i+1)*new%h*d2y(i)  /6d0 + x(i+1)*y(i)/new%h
    D = D - x(i)  **3*d2y(i+1)/(6d0*new%h) + x(i)  *new%h*d2y(i+1)/6d0 - x(i)*y(i+1)/new%h

    new%spline(i)%A = A
    new%spline(i)%B = B
    new%spline(i)%C = C
    new%spline(i)%D = D
  end do
! External "spline" above gas residual 
  new%spline(N+2)%A = 0d0
  new%spline(N+2)%B = 0d0
  new%spline(N+2)%C = 0d0
  new%spline(N+2)%D = y(N+1)

end function RPFSplineCtor

! **************************************************************************** !

subroutine RPFSplineRelativePermeability(this, liquid_saturation, &
                                   relative_permeability, dkr_sat, option)
  implicit none

  class(rpf_spline_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: relative_permeability, dkr_sat
  type(option_type), intent(inout) :: option
  PetscInt :: i
  PetscReal :: x, y, dy

! If provided with an array, these functions can be done in parallel
! Ceiling, min, and max are hardware instructions

! Aliases because dummy variable names are absurdly long. Compiler should inline
  x = liquid_saturation

! Hash function for constant spacing between the residuals 
! 1   corresponds to "spline" below liquid residual
! N+2 corresponds to "spline" above gas residual
  i = min(max(ceiling((x - this%Sr)/this%h),1),this%N+2)

! Cubic polynomial by Horner's method
  y  = ((this%spline(i)%A*x + this%spline(i)%B)*x + this%spline(i)%C)*x + this%spline(i)%D
  dy = (3d0*this%spline(i)%A*x + 2d0*this%spline(i)%B) * x + this%spline(i)%C

! Aliases because dummy variable names are absurdly long. Compiler should inline
  relative_permeability = y
  dkr_sat = dy

end subroutine RPFSplineRelativePermeability

! **************************************************************************** !

subroutine RPFSpline(x, y, n, yp1, ypn, y2)
  PetscInt :: n
  PetscReal :: yp1, ypn, x(n), y(n), y2(n)
  PetscInt :: i

  PetscReal :: p, qn, sig, un,u(n)

! This subroutine implies the 1st derivatives at the endpoints are set
  y2(1) = -0.5d0
  u(1)  = (3d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  do i =2, n-1
    sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
    p   = sig*y2(i-1)+2
    y2(i) = (sig-1d0)/p
    u(i)  = (6d0*((Y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
          / (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  end do
  qn = 0.5d0
  un = (3d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))

! Back substitution to find 2nd derivatives
  y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1d0)
  do i = n-1, 1, -1
    y2(i) = y2(i)*y2(i+1)+u(i)
  end do
end subroutine

end module

