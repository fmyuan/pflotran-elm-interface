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
! Date:   05/05/2022
!
! This module contains the methods to construct, evaluate, and deconstruct 
! cubic spline interpolation for any monotonic capillary pressure or relative
! permeability function or data set.  Cubic splines are the simpliest
! arithmeticl expression that provides continuity and smoothness over any data.
! In the absence of hysteresis , the free energy and hence capillary pressure
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

type, private :: knot_type
  PetscReal :: x, y
  type(knot_type), pointer :: next
end type knot_type

! **************************************************************************** !

type, public :: knot_queue_type
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
! x-reference points are stored separately to reduce cache size during search
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
    procedure, public :: Test              => SFPCHIPTest
    procedure, public :: Name              => SFPCHIPName
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
    procedure, public :: Name                 => RPFPCHIPName
    final :: RPFPCHIPDtor
end type rpf_pchip_type

! **************************************************************************** !

private PCHIPCoefficients ! Calculates SF and RPF coefficients via SLATEC/PCHIP
 
private SFPCHIPAllocate
public  SFPCHIPCreate
public  SFPCHIPCtorData
public  SFPCHIPCtorFunction
public  SFPCHIPCtorQueue

private RPFPCHIPAllocate
public  RPFPCHIPCreate
public  RPFPCHIPCtorFunction
public  RPFPCHIPCtorData
public  RPFPCHIPCtorQueue

contains

! **************************************************************************** !
! Knot queue to support the parser
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
! Common PCHIP Coefficient Routine
! **************************************************************************** !

subroutine PCHIPCoefficients(N, x, y, dy, c2, c3)
  PetscInt, intent(in)   :: N
  PetscReal, intent(in)  :: x(N), y(N)
  PetscReal, intent(out) :: dy(N-1), c2(N-1), c3(N-1)
  PetscInt :: I
  PetscReal :: h, delta, del1, del2

  ! Use SLATEC/PCHIP/PCHIM to calculate first derivatives for monotonic splines
  I = 0
  call PCHIM(N, x, y, dy, 1, I)

  ! 2nd and 3rd derivatives are defined by the Hermite polynomials, however
  ! unlike SLATEC/PCHIP/CHFDV, polynomial coefficients will be cached
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
  allocate(new)
  call new%Init()
end function

! ************************************************************************** !

subroutine SFPCHIPInit(this)
  class(sf_pchip_type) :: this
! Common values regardless of ctor
  this%analytical_derivative_available = PETSC_TRUE
  this%calc_int_tension = PETSC_FALSE
! Unused objects in base class that must be nullified
  nullify(this%sat_poly)
  nullify(this%pres_poly)
end subroutine SFPCHIPInit

! **************************************************************************** !

function SFPCHIPAllocate(N) result (new)
  implicit none
  PetscInt, intent(in) :: N
  class(sf_pchip_type), pointer :: new

  nullify(new)

  ! No interpolation is possible with less than 2 knots
  if (N < 2) return

  ! Allocate spline object
  allocate(new)
  if (.not. associated(new)) return
  call new%Init()
  new%N = N

  ! Allocate subordinate dynamic objects
  allocate(new%x(N))
  allocate(new%coef(N))
  if (.not. allocated(new%x) .or. .not. allocated(new%coef)) then
    ! Final method will deallocate if only one is allocated
    deallocate(new)
    nullify(new)
  end if

end function

! **************************************************************************** !

subroutine SFPCHIPDtor(this)
  implicit none
  type(sf_pchip_type) :: this
  ! Deallocate subordinate dynamic objects
  ! These will not be allocated if "create" is used instead of a ctor
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
  PetscReal :: dyi
  type(option_type) :: option

  new => SFPCHIPAllocate(N+1)
  if (.not. associated(new)) return

  ! Copy attributes of base class
  new%Sr    = sf_analytic%Sr
  new%Pcmax = sf_analytic%Pcmax

  ! Generate N knots starting at 0 and ending at 1
  do I = 1, N+1
    new%x(i) = dble(I-1)/dble(N)
    call sf_analytic%CapillaryPressure(new%x(i), new%coef(i)%y, dyi, option)
  end do
 
! Calculate the polynomial coefficients
  call PCHIPCoefficients(N+1, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function SFPCHIPCtorFunction

! **************************************************************************** !

function SFPCHIPCtorQueue(queue) result (new)
  implicit none
  type(knot_queue_type), intent(inout) :: queue
  class(sf_pchip_type), pointer :: new

  PetscInt :: i, N
  PetscReal :: x, y

  N = queue%depth()
  print *, "Queue depth", N

  new => SFPCHIPAllocate(N)
  if (.not. associated(new)) return

! Pack queue into the arrays
  do i = 1, N
    call queue%dequeue(new%x(i), new%coef(i)%y)
    print *, "Knot:", new%x(i), new%coef(i)%y
  end do

! Standard attributes of base class
  new%Sr    = new%x(1)
  new%Pcmax = new%coef(1)%y

! Calculate the polynomial coefficients
  call PCHIPCoefficients(N, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function SFPCHIPCtorQueue

! **************************************************************************** !

function SFPCHIPCtorData(N, x, y) result (new)
  implicit none
  PetscInt, intent(in) :: N
  PetscReal, intent(in) :: x(N), y(N)
  class(sf_pchip_type), pointer :: new

  new => SFPCHIPAllocate(N)
  if (.not. associated(new)) return

  ! Standard attributes of base class
  new%Sr = x(1)
  new%Pcmax = y(1)

  ! Vector copy into object space
  new%x = x
  new%coef%y = y

! Calculate the polynomial coefficients
  call PCHIPCoefficients(N, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function SFPCHIPCtorData

! **************************************************************************** !

function SFPCHIPName(this) result (name)
  class(sf_pchip_type) :: this
  character(len=MAXSTRINGLENGTH) :: name

  name = 'PCHIP Splines'
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

  x = min(max(liquid_saturation,this%x(1)),this%x(this%N))

! Binary search, compiler can optimize this to be branchless
! To vectorize, replace branch with merge()
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
! Common values regardless of ctor
  this%analytical_derivative_available = PETSC_TRUE
! Unused objects in base class that must be nullified
  nullify(this%poly)
end subroutine RPFPCHIPInit

! **************************************************************************** !

function RPFPCHIPAllocate(N) result (new)
  ! Perform allocation for both data or function defined splines
  implicit none
  PetscInt, intent(in) :: N
  class(rpf_pchip_type), pointer :: new

  nullify(new)

  ! No interpolation is possible with less than 2 knots
  if (N < 2) return

  ! Allocate spline base
  allocate(new)
  if (.not. associated(new)) return
  call new%Init()
  new%N = N

  ! Allocate subordinate dynamic objects
  allocate(new%x(N))
  allocate(new%coef(N))
  if (.not. allocated(new%x) .or. .not. allocated(new%coef)) then
    ! Final method will deallocate if only one is allocated
    deallocate(new)
    nullify(new)
  end if
end function

! **************************************************************************** !

subroutine RPFPCHIPDtor(this)
  type(rpf_pchip_type) :: this
  ! Deallocate coefficient array
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
  PetscReal :: dyi, span
  type(option_type) :: option

  new => RPFPCHIPAllocate(N+1)
  if (.not. associated(new)) return

  ! Standard attributes of base class
  new%Sr  = rpf_analytic%Sr
  new%Srg = rpf_analytic%Srg
  ! Liquid relative perm functions may have an "uninitialized" gas residual
  if (rpf_analytic%Srg == UNINITIALIZED_DOUBLE) new%Srg = 0d0

  ! N+1 knots starting at Sr and ending at 1-Srg
  span = 1d0 - new%Srg - new%Sr
  do I = 1, N
    new%x(i) = new%Sr + span*dble(I-1)/dble(N)
    call rpf_analytic%RelativePermeability(new%x(I), new%coef(I)%y, dyi, option)
  end do
  ! Ensure the endpoint is exact
  new%x(N+1) = 1d0 - new%Srg
  call rpf_analytic%RelativePermeability(new%x(N+1), new%coef(N+1)%y, dyi, option)

  ! Calculate the quadratic and cubic coefficients
  call PCHIPCoefficients(N+1, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function RPFPCHIPCtorFunction

! **************************************************************************** !

function RPFPCHIPCtorData(N, x, y) result (new)
  implicit none
  PetscInt, intent(in) :: N
  PetscReal, intent(in) :: x(N), y(N)
  class(rpf_pchip_type), pointer :: new

  new => RPFPCHIPAllocate(N)
  if (.not. associated(new)) return

  ! Base class assignments
  new%Sr  = x(1)
  new%Srg = 1d0 - x(N)

  ! Vector copy into object space
  new%x = x
  new%coef%y = y

  ! Calculate polynomial coefficients
  call PCHIPCoefficients(N, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function RPFPCHIPCtorData

! **************************************************************************** !

function RPFPCHIPCtorQueue(queue) result (new)
  implicit none
  class(knot_queue_type), intent(inout) :: queue
  class(rpf_pchip_type), pointer :: new

  PetscInt :: i, N

  N = queue%depth()

! Allocate object and arrays to fit queue
  new => RPFPCHIPAllocate(N)
  if (.not. associated(new)) return

! Pack queue into the arrays
  do i = 1, N
    call queue%dequeue(new%x(i), new%coef(i)%y)
  end do

! Base class assignments
  new%Sr  = new%x(1)
  new%Srg = 1d0 - new%x(N)

! Calculate polynomial coefficients
  call PCHIPCoefficients(N, new%x, new%coef%y, new%coef%dy, new%coef%c2, new%coef%c3)

end function RPFPCHIPCtorQueue

! **************************************************************************** !

function RPFPCHIPName(this) result (name)
  class(rpf_pchip_type) :: this
  character(len=MAXSTRINGLENGTH) :: name

  name = 'PCHIP Splines'
end function

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
! To vectorize, replace branch with merge()
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
    write(87,*) i, this%x(i), this%coef(i)%y, this%coef(i)%dy
  end do
  close(87)
    
end subroutine RPFPCHIPTest

! **************************************************************************** !

end module
