module Characteristic_Curves_spline_module
#include "petsc/finclude/petscsys.h"

use petscsys ! Necessary for PETSC_TRUE / PETSC_FALSE
use Characteristic_Curves_Base_module   ! Needed to define base type
use Option_module ! Needed for Verify and unused arguments in Pc and Sl
use slatec_pchip_module
use PFLOTRAN_constants_module ! Needed for UNINITIALIZED_DOUBLE

implicit none
private

! **************************************************************************** !
!
! Author: Matthew Paul
! Date:   05/05/2022
!
! This module contains the methods to construct, evaluate, and deconstruct 
! cubic spline interpolation for any monotonic capillary pressure or relative
! permeability function or data set.  Cubic splines are the simpliest
! arithmeticl expression that provides continuity and smoothness over any data.
! In the abscense of hysterisis , the free energy and hence capillary pressure
! must be a monotonic function of saturation. In addition, monotonicity is
! necessary to avoid introducing degenerate roots for the Newton solver.
!
! Monotonic splines are calculated using the Piecewise Cubic Hermite
! Interpolation Polynomial (PCHIP) package from SLATEC, a public domain library.
! PCHIP was optimized for graphing, thus the evaluation function assumes the
! points being evaluated are in an ordered list, and a linear search was used
! to evaluate the entire list.

! Here, the PCHIP calculated polynomials are cached and a binary search
! is used to search the key values, which is more efficient when there is no
! correlation between consequentive function evaluations. Binary search can be
! make branchless and vectorized using conditional move / merge() intrinsic,
! but is left as the branched form for scalar function calls.
!
! References:
!
! Ross, P (1992) "Cubic approximation of hydraulic properties for simulations of
! unsaturated flow", Water Resour. Res. 28(10):2617-2620.
! doi: 10.1029/92WR01310a
!
! Fritsch, FN and Carlson, RE (1980) "Monotone piecewise cubic interpolation",
! SIAM J. Numer. Anal., 17(2):238-246.
! doi: 10.11377/0717021
!
! **************************************************************************** !

type, private :: cubic_type
! Polynomial coefficients are stored in a structure for data locality
! x-reference values are stored separately to reduce cache size during search
  PetscReal :: y, dy, c2, c3
end type
 
! **************************************************************************** !
type, public, extends(sat_func_base_type) :: sf_pchip_type
  private
    PetscInt  :: N ! Number of knots
    PetscReal, dimension(:), allocatable :: x ! Reference points
    type(cubic_type), dimension(:), allocatable :: coef ! Coefficients
  contains
    procedure, public :: Init              => SFPCHIPInit
    procedure, public :: CapillaryPressure => SFPCHIPCapillaryPressure
    procedure, public :: Test              => SFPCHIPTest
    procedure, public :: Dtor              => SFPCHIPDtor
end type

! **************************************************************************** !

type, public, extends(rel_perm_func_base_type) :: rpf_pchip_type
  private
    PetscInt :: N ! Number of knots
    PetscReal, dimension(:), allocatable :: x ! Reference points
    type(cubic_type), dimension(:), allocatable :: coef ! Coefficients
  contains
    procedure, public :: Init                 => RPFPCHIPInit
    procedure, public :: RelativePermeability => RPFPCHIPRelativePermeability
    procedure, public :: Test                 => RPFPCHIPTest
    procedure, public :: Dtor                 => RPFPCHIPDtor
end type

! **************************************************************************** !

private PCHIPCoefficients
 
private SFPCHIPAllocate
public  SFPCHIPCtorFunction
public  SFPCHIPCtorData

private RPFPCHIPAllocate
public  RPFPCHIPCtorFunction
public  RPFPCHIPCtorData

contains

! **************************************************************************** !
! PCHIP wrapper
! **************************************************************************** !

subroutine PCHIPCoefficients(N, x, y, dy, c2, c3)
  PetscInt, intent(in)   :: N
  PetscReal, intent(in)  :: x(N), y(N)
  PetscReal, intent(out) :: dy(N-1), c2(N-1), c3(N-1)
  PetscInt :: I
  PetscReal :: h, delta, del1, del2

  ! First derivatives are calculating using PCHIP library
  I = 0
  call PCHIM(N, x, y, dy, 1, I)

  ! Polynomial coefficients are cached
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

subroutine SFPCHIPInit(this)
  implicit none
  class(sf_pchip_type) :: this
  ! This method is intentionally left blank.
end subroutine SFPCHIPInit

! **************************************************************************** !

function SFPCHIPAllocate(N) result (new)
  ! Perform allocation for both data or function defined splines
  implicit none
  PetscInt, intent(in) :: N
  class(sf_pchip_type), pointer :: new

  ! No interpolation is possible with less than 2 knots
  if (N < 2) then
    nullify(new)
    return
  end if

  ! Allocate spline object
  allocate(new)
  if (.not. associated(new)) return

  new%N = N

  ! Allocate subordinate arrays
  allocate(new%x(N))
  if (.not. allocated(new%x)) then
    deallocate(new)
    nullify(new)
    return
  end if

  allocate(new%coef(N))
  if (.not. allocated(new%coef)) then
    deallocate(new)
    nullify(new)
    return
  end if

! Unused pointers in base class that must be nullified
! Could be eliminated if base class is refactored
  nullify(new%sat_poly)
  nullify(new%pres_poly)

end function

! **************************************************************************** !

function SFPCHIPCtorFunction(N, sf_analytic) result (new)
  implicit none
  class(sf_pchip_type), pointer :: new
  class(sat_func_base_type), intent(in) :: sf_analytic
  PetscInt, intent(in) :: N
  PetscInt  :: I, ierr
  PetscReal :: dyn
  type(option_type) :: option

  new => SFPCHIPAllocate(N)
  if (.not. associated(new)) return

  ! Standard attributes of base class
  new%Sr = 0d0
  new%calc_int_tension = PETSC_FALSE  ! Default, can be mutated
  new%analytical_derivative_available = PETSC_TRUE
  new%Pcmax = sf_analytic%Pcmax

  ! Generate N knots starting at 0 and ending at 1
  do I = 1, N
    new%x(i) = dble(I-1)/dble(N-1)
    call sf_analytic%CapillaryPressure(new%x(i), new%coef(i)%y, dyn, option)
  end do
 
! Calculate the polynomial coefficients
  call PCHIPCoefficients(N, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function SFPCHIPCtorFunction

! **************************************************************************** !

function SFPCHIPCtorData(N, x, y) result (new)
  implicit none
  class(sf_pchip_type), pointer :: new
  PetscInt, intent(in) :: N
  PetscReal, intent(in) :: x(N), y(N)

  new => SFPCHIPAllocate(N)
  if (.not. associated(new)) return

  ! Standard attributes of base class
  new%Sr = x(1)
  new%calc_int_tension = PETSC_FALSE  ! Default, can be mutated
  new%analytical_derivative_available = PETSC_TRUE
  new%Pcmax = y(1)

  new%x = x
  new%coef%y = y

! Calculate the polynomial coefficients
  call PCHIPCoefficients(N, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function SFPCHIPCtorData

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

  x = min(max(liquid_saturation,this%x(1)),this%x(this%N))

! Binary search, compiler can optimize this to be branchless
! Replace branch with merge() when vectorizing
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

  x  = x - this%x(i)

  capillary_pressure    =    this%coef(i)%y + x* &
                        (   this%coef(i)%dy + x* &
                        (   this%coef(i)%c2 + x*this%coef(i)%c3))
  dPc_dSatl             =   this%coef(i)%dy + x* &
                        ( 2*this%coef(i)%c2 + 3*x*this%coef(i)%c3)

end subroutine SFPCHIPCapillaryPressure

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

  PetscInt, parameter :: num_values = 1000
  PetscReal :: Sw, Pc, dPc_dSw
  PetscInt :: i

  write(string,*) cc_name
  string = trim(cc_name) // '_pc_from_sat.dat'
  open(unit=86,file=string)
  write(86,*) '#Sw, Pc, dPc/dSw'
  do i = 0, num_values
    Sw = dble(i)/dble(num_values)
    call this%CapillaryPressure(Sw, Pc, dPc_dSw, option)
    write(86,*) Sw, Pc, dPc_dSw
  enddo
  close(86)

! Write knots and derivatives to file
  write(string,*) cc_name
  string = trim(cc_name) // '_pc_knots.dat'
  open(unit=87,file=string)
  write(87,*) '#Index, Sw, Pc, dPc/dSw'
  do i = 1, this%N
    write(87,*) i, this%x(i), this%coef(i)%y, this%coef(i)%dy
  end do

end subroutine SFPCHIPTest

! **************************************************************************** !

subroutine SFPCHIPDtor(this)
  implicit none
  class(sf_pchip_type) :: this
  deallocate(this%x)
  deallocate(this%coef)
end subroutine SFPCHIPDtor

! **************************************************************************** !
! Relative Permeability PCHIP Methods
! **************************************************************************** !

subroutine RPFPCHIPInit(this)
  implicit none
  class(rpf_pchip_type) :: this
  ! This method is intentionally left blank.
end subroutine RPFPCHIPInit

! **************************************************************************** !

function RPFPCHIPAllocate(N) result (new)
  ! Perform allocation for both data or function defined splines
  implicit none
  PetscInt, intent(in) :: N
  class(rpf_pchip_type), pointer :: new

  ! No interpolation is possible with less than 2 knots
  if (N < 2) then
    nullify(new)
    return
  end if

  ! Allocate spline base
  allocate(new)
  if (.not. associated(new)) return

  new%N = N

  ! Allocate reference point array
  allocate(new%x(N))
  if (.not. allocated(new%x)) then
    deallocate(new)
    nullify(new)
    return
  end if

  ! Allocate coefficient array
  allocate(new%coef(N))
  if (.not. allocated(new%coef)) then
    deallocate(new)
    nullify(new)
    return
  end if

! Unused pointers in base class that must be nullified
  nullify(new%poly)

end function

! **************************************************************************** !

function RPFPCHIPCtorFunction(N, rpf_analytic) result (new)
  implicit none
  class(rpf_pchip_type), pointer :: new
  class(rel_perm_func_base_type), intent(in) :: rpf_analytic
  PetscInt, intent(in) :: N
  PetscInt  :: I
  PetscReal :: dyn, span
  type(option_type) :: option

  new => RPFPCHIPAllocate(N)
  if (.not. associated(new)) return

  ! Standard attributes of base class
  new%Sr  = rpf_analytic%Sr
  new%Srg = rpf_analytic%Srg
  ! Liquid relative perm functions may have an "uninitialized" gas residual
  if (rpf_analytic%Srg == UNINITIALIZED_DOUBLE) new%Srg = 0d0
  new%analytical_derivative_available = PETSC_TRUE

  ! N knots starting at Sr and ending at 1-Srg
  span = 1d0 - new%Srg - new%Sr
  do I = 1, N-1
    new%x(i) = new%Sr + span*dble(I-1)/dble(N-1)
    call rpf_analytic%RelativePermeability(new%x(I), new%coef(I)%y, dyn, option)
  end do
  ! Ensure the endpoint is exact
  new%x(N) = 1d0 - new%Srg
  call rpf_analytic%RelativePermeability(new%x(N), new%coef(N)%y, dyn, option)

  ! Calculate the quadratic and cubic coefficients
  call PCHIPCoefficients(N, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function RPFPCHIPCtorFunction

! **************************************************************************** !

function RPFPCHIPCtorData(N, x, y) result (new)
  implicit none
  class(rpf_pchip_type), pointer :: new
  PetscInt, intent(in) :: N
  PetscReal, intent(in) :: x(N), y(N)

  new => RPFPCHIPAllocate(N)
  if (.not. associated(new)) return

  ! Base class assignments
  new%Sr  = x(1)
  new%Srg = 1d0 - x(N)
  new%analytical_derivative_available = PETSC_TRUE

  new%x = x
  new%coef%y = y

  ! Calculate polynomial coefficients
  call PCHIPCoefficients(N, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function RPFPCHIPCtorData
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

  x = min(max(liquid_saturation,this%x(1)), this%x(this%N))

! Binary search, compiler can optimize this to be branchless
! Replace branch with merge() when vectorizing
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

  x  = x - this%x(i)

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
  string = trim(cc_name) // '_' //  trim(phase) // '_Kr_knots.dat'
  open(unit=87,file=string)
  write(87,*) '#Index, x, y, dy'
  do i = 1, this%N
    write(87,*) this%x(i), this%coef(i)%y, this%coef(i)%dy
  end do
  close(87)
    
end subroutine RPFPCHIPTest

! **************************************************************************** !

subroutine RPFPCHIPDtor(this)
  class(rpf_pchip_type) :: this
  deallocate(this%x)
  deallocate(this%coef) 
end subroutine RPFPCHIPDtor

end module

