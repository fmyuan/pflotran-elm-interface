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
! |-->sf_spline_type        Generated cubic spline type
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
    procedure, public :: Test              => SFSplineTest
end type

! **************************************************************************** !
public RPFSplineCtor
! **************************************************************************** !

type, public, extends(rel_perm_func_base_type) :: rpf_spline_type
  private
    PetscInt  :: N ! Number of splines
    PetscReal :: h ! Saturation interval width
    PetscReal :: span
  contains
    procedure, public :: Init                 => RPFSplineInit
    procedure, public :: RelativePermeability => RPFSplineRelativePermeability
    procedure, public :: Test                 => RPFSplineTest
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
  PetscReal :: dy1, dyn
  type(option_type) :: option

  ! Validate the number of knots - could set a default here
  if (N < 2) then
    nullify(new)
    return
  end if

  ! Allocate space for base structure
  allocate(new)
  if (.not. associated(new)) return

  ! Allocate knot table
  allocate(new%spline(N))
  if (.not. allocated(new%spline)) then
    deallocate(new)
    nullify(new)
    return
  end if

  ! Standard attributes of base class
  nullify(new%sat_poly)
  nullify(new%pres_poly)
  new%Sr = 0d0
  new%calc_int_tension = PETSC_FALSE  ! Default, can be mutated
  new%analytical_derivative_available = PETSC_TRUE
  new%Pcmax = sf_analytic%Pcmax ! Copy by not really used

  ! N knots
  new%N = N
  new%spline(1)%x = 0d0
  call sf_analytic%CapillaryPressure(new%spline(1)%x, new%spline(1)%y, dy1, option)
  do I = 2, N
    new%spline(i)%x = dble(I-1)/dble(N-1)
    call sf_analytic%CapillaryPressure(new%spline(i)%x, new%spline(i)%y, dyn, option)
  end do

  ! Calculate 2nd derivatives assuming "natural" splines
  ! Considered using derivatives at end points, but some "analytic" functions
  ! have singularities and/or errors at end points
  call spline(new%spline%x, new%spline%y, new%N, new%spline%d2y)

end function SFSplineCtor

! **************************************************************************** !

subroutine SFSplineCapillaryPressure(this, liquid_saturation, &
                                   capillary_pressure, dPc_dSatl, option)
  implicit none

  class(sf_spline_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option
  PetscInt :: i, j, k
  PetscReal :: h, a, b
  PetscReal :: x, y, dy

  x = min(max(liquid_saturation,0d0),1d0)

! Consider replacing with hash function
! i = 1
! j = this%n
! do while (j - i > 1)
!   k = (i + j)/2
!   if (this%spline(k)%x > x) then
!     j = k
!   else
!     i = k
!   end if
! end do
  i = min(max(ceiling(x * dble(this%N-1)),1),this%N-1)
  j = i + 1

  h = this%spline(j)%x - this%spline(i)%x
  a = (this%spline(j)%x-x)/h
  b = (x-this%spline(i)%x)/h
  y = a*this%spline(i)%y + b*this%spline(j)%y &
    + ((a**3-a)*this%spline(i)%d2y + (b**3-b)*this%spline(j)%d2y)*h**2/6d0

  dy = (this%spline(j)%y - this%spline(i)%y)/h &
     + (-(3d0*a**2-1d0)*this%spline(i)%d2y + (3d0*b**2-1d0)*this%spline(j)%d2y)*h/6d0

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
  PetscReal :: Sw, Pc, dPc_dSw
  PetscInt :: i

  write(string,*) cc_name
  string = trim(cc_name) // '_pc_from_sat.dat'
  open(unit=86,file=string)
  write(86,*) '#Sw,       Pc,       dPc/dSw'

 ! calculate capillary pressure as a function of saturation
  do i = 0, num_values
    Sw = dble(i)/dble(num_values)
    call this%CapillaryPressure(Sw, Pc, dPc_dSw, option)
    write(86,*) Sw, Pc, dPc_dSw
  enddo
  close(86)
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
  PetscReal :: dy1, dyn
  PetscReal :: span
  type(option_type) :: option

  ! Validate the number of knots - could set a default here
  if (N < 2) then
    nullify(new)
    return
  end if

  allocate(new)
  if (.not. associated(new)) return

  allocate(new%spline(N))
  if (.not. allocated(new%spline)) then
    deallocate(new)
    nullify(new)
    return
  end if

  ! Standard attributes of base class
  nullify(new%poly)
  new%Sr = rpf_analytic%Sr
  new%Srg = rpf_analytic%Srg
  ! Liquid relative perm functions may have an "uninitialized" gas residual
  new%analytical_derivative_available = PETSC_TRUE
  ! Initialize public attributes
  if (rpf_analytic%Srg == UNINITIALIZED_DOUBLE) new%Srg = 0d0

  ! N knots
  new%N = N
  new%span = 1d0 - new%Sr - new%Srg

  ! Generate N knots between Sr and 1-Srg, saving derivatives at endpoints
  new%spline(1)%x = new%Sr
  call rpf_analytic%RelativePermeability(new%spline(1)%x, new%spline(1)%y, dy1, option)
  do I = 2, N - 1
    new%spline(i)%x = new%Sr + new%span*dble(I-1)/dble(N-1)
    call rpf_analytic%RelativePermeability(new%spline(I)%x, new%spline(I)%y, dyn, option)
  end do
  ! Need to be certain this point is exact
  new%spline(N)%x = 1d0 - new%Srg
  call rpf_analytic%RelativePermeability(new%spline(N)%x, new%spline(N)%y, dyn, option)

  call spline(new%spline%x, new%spline%y, N, new%spline%d2y)

end function RPFSplineCtor

! **************************************************************************** !

subroutine RPFSplineRelativePermeability(this, liquid_saturation, &
                                   relative_permeability, dkr_sat, option)
  implicit none

  class(rpf_spline_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: relative_permeability, dkr_sat
  type(option_type), intent(inout) :: option
  PetscInt :: i, j, k
  PetscReal :: h, a, b
  PetscReal :: x, y, dy

! If provided with an array, these functions can be done in parallel
! Ceiling, min, and max are hardware instructions

! Aliases because dummy variable names are absurdly long. Compiler should inline
! Cut liquid saturation off at residuals

  x = min(max(liquid_saturation,this%Sr),1d0-this%Srg)

!  i = 1
!  j = this%n
!  do while (j - i > 1)
!    k = (i + j)/2
!   if (this%spline(k)%x > x) then
!     j = k
!   else
!     i = k
!   end if
! end do

  i = min(max(ceiling((x-this%Sr)*dble(this%N-1)/this%span), 1), this%N-1)
  j = i + 1

  h = this%spline(j)%x - this%spline(i)%x
  a = (this%spline(j)%x-x)/h
  b = (x-this%spline(i)%x)/h
  y = a*this%spline(i)%y + b*this%spline(j)%y &
    + ((a**3-a)*this%spline(i)%d2y + (b**3-b)*this%spline(j)%d2y)*h**2/6d0
  dy = (this%spline(j)%y - this%spline(i)%y)/h &
     + (-(3d0*a**2-1d0)*this%spline(i)%d2y + (3d0*b**2-1d0)*this%spline(j)%d2y)*h/6d0

  relative_permeability = y
  dkr_sat = dy

end subroutine RPFSplineRelativePermeability

! **************************************************************************** !

subroutine RPFSplineTest(this,cc_name,phase,option)
  use Option_module
  use Material_Aux_module
  use PFLOTRAN_Constants_module
  implicit none
  class(rpf_spline_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option
  character(len=MAXWORDLENGTH) :: phase
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt, parameter :: num_values = 10000
  PetscReal :: Sw, Kr, dKr_dSw
  PetscInt :: i

  write(string,*) cc_name
  string = trim(cc_name) // '_' //  trim(phase) // '_Kr.dat'
  open(unit=86,file=string)
  write(86,*) '#Sw,       Kr,       dKr/dSw'

 ! calculate capillary pressure as a function of saturation
  do i = 0, num_values
    Sw = dble(i)/dble(num_values)
    call this%RelativePermeability(Sw, Kr, dKr_dSw, option)
    write(86,*) Sw, Kr, dKr_dSw
  enddo
  close(86)

  write(string,*) cc_name
  string = trim(cc_name) // '_' //  trim(phase) // '_Kr_spline.dat'
  open(unit=87,file=string)
  write(86,*) '#Index, x_max, A, B, C, D'
  do i = 1, this%N
    write(87,*) this%spline(i)%x, this%spline(i)%y, this%spline(i)%d2y
  end do
  close(87)
    
end subroutine RPFSplineTest

end module

