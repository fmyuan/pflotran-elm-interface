module Derivatives_utilities_module 

! Collection of simple routines for applying common calculus
! rules.

#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none

  private

  Public ::  ProdRule, &
             ProdRule3, &
             DivRule, &
             DivRule1, &
             PowerRule!, &
             !ProdDivrule

! ************************************************************************** !
!                 General explaination for schema here
!
! Consider two (or more) variables X, Y, could be for example density or
! enthalpy or any auxvar type variable in general.
!
! Assume we have corresponding arrays D_X, D_Y, where D_X(i) is the 
! partial derivative of variable X with respect to some other variable
! indexed by i. Typically these would be the solution varibles like
! pressure or saturation etc., but we do not care here what these 
! variables are here, only that they are indexed.

! In this module we define simple routines for obtaining corresponding
! arrays of partial derivatives D_A, where A is related to X and Y by:

! A = X*Y          ( see ProdRule() )
! A = X/Y          ( see DivRule()  )
! A = X*Y*Z        ( see ProdRule3())
! A = 1/X          ( see DivRule1() )

! ************************************************************************** !

  contains

! ************************************************************************** !
function  ProdRule(x, dx, y, dy,n)

! Compute derivatives of x*y.
! Assume x and y are functions and there are some n variables
! t_1, t_2, ... 2_n, so that 
!
! x = x(t_1, .... t_n)
! y = y(t_1, .... t_n)
!
! Assume we are given the partial derivatives in the arrays
! dx and dy, i.e.,
!
! dx = [dx/dt_1,  .... , dx_dt_n]
! dy = [dy/dt_1,  .... , dy_dt_n]
!
! Then we compute and return an array:
!
! [d(xy)/dt_1, ... , d(xy)/dt_n]

implicit none

  PetscInt :: n
  PetscReal :: x, y
  PetscReal, dimension(1:n) :: dx, dy
  PetscReal, dimension(1:n) :: ProdRule

  PetscInt :: i

  Prodrule = 0.d0

  do i = 1,n
    ProdRule(i) = x*dy(i) + y*dx(i)  
  enddo

end function ProdRule

! ************************************************************************** !

function  ProdRule3(x, dx, y, dy, z, dz, n)

implicit none

  PetscInt :: n
  PetscReal :: x, y, z
  PetscReal, dimension(1:n) :: dx, dy, dz
  PetscReal, dimension(1:n) :: d_yz
  PetscReal, dimension(1:n) :: ProdRule3
  PetscReal :: yz

  d_yz = ProdRule(y, dy, z, dz, n)
  yz = y*z
  ProdRule3 = ProdRule(x, dx, yz, d_yz, n)

end function  ProdRule3

! ************************************************************************** !

function scal_d_prodrule(x,dx,y,dy)

  implicit none
  PetscReal :: x,dx,y,dy
  PetscReal :: scal_d_prodrule

  scal_d_prodrule = dx*y + x*dy

end function scal_d_prodrule

! ************************************************************************** !

function DivRule(x, dx, y, dy, n)

!! derivatives of x/y, see comments in d_prodrule_2 for general schema
!! of partial derivatives and etc.

implicit none
  PetscInt :: n
  PetscReal :: x, y
  PetscReal, dimension(1:n) :: dx, dy
  PetscReal, dimension(1:n) :: DivRule 

  PetscReal :: denompart
  PetscInt :: i

  DivRule = 0.d0
  if( abs(y)<EPSILON(x) ) then
    return
  endif

  denompart = 1.d0/y/y
  do i = 1,n
    DivRule(i) =  dx(i)/y - x*dy(i)*denompart
  end do


end function DivRule 

! ************************************************************************** !

function DivRule1(x, dx, n)

! derivatives of 1/x, see comments in above for general schema
! of partial derivatives and so on.

implicit none
  PetscInt :: n
  PetscReal :: x
  PetscReal, dimension(1:n) :: dx
  PetscReal, dimension(1:n) :: DivRule1

  PetscReal :: denompart
  PetscInt :: i

  DivRule1 = 0.d0
  if( abs(x)<EPSILON(x) ) then
    return
  endif

  denompart = 1.d0/x/x
  do i = 1,n
    DivRule1(i) = -1.d0*dx(i)*denompart
  end do


end function DivRule1

! ************************************************************************** !

! not currently used but may come in handy at some point:

function PowerRule(x, dx, alpha, n)

! derivatives of: x^alpha

implicit none
  PetscInt :: n
  PetscReal :: x, alpha
  PetscReal, dimension(1:n) :: dx
  PetscReal, dimension(1:n) :: PowerRule

  PetscReal :: constpart
  PetscInt :: i

  constpart = alpha*(x**(alpha-1.d0))
  do i=1,n
    PowerRule(i) = dx(i)*constpart 
  end do

end function PowerRule

! ************************************************************************** !

!!! to be completed
#if 0
function  ProdDivRule(x, dx, y, dy, z, dz, n)


implicit none

  PetscInt :: n
  PetscReal :: x, y, z
  PetscReal, dimension(1:n) :: dx, dy, dz
  PetscReal, dimension(1:n) :: ProdDivRule

  PetscInt :: i

  ProdDivRule = 0.d0

  do i = 1,n
    ProdDivRule(i) = x*dy(i) + y*dx(i)  
  enddo

end function ProdDivRule
#endif

! ************************************************************************** !


end module Derivatives_utilities_module 
