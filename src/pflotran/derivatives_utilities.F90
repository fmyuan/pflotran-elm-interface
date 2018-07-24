module Derivatives_utilities_module 

! Collection of simple routines for applying common calculus
! rules.

#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none

  private

  public  :: d_prodrule_2, &
             d_prodrule_3, &
             !d_divrule, &
             d_power, &
             scal_d_prodrule, &
             ProdRule, &
             ProdRule3, &
             DivRule, &
             DivRule1

  contains

! ************************************************************************** !
function  ProdRule(x, dx, y, dy,n)

!! Compute derivatives of x*y.
!! Assume x and y are functions and there are some n variables
!! t_1, t_2, ... 2_n, so that 
!!
!! x = x(t_1, .... t_n)
!! y = y(t_1, .... t_n)
!!
!! Assume we are given the partial derivatives in the arrays
!! dx and dy, i.e.,
!!
!! dx = [dx/dt_1,  .... , dx_dt_n]
!! dy = [dy/dt_1,  .... , dy_dt_n]
!!
!! Then we compute and return an array:
!!
!! [d(xy)/dt_1, ... , d(xy)/dt_n]

implicit none

  PetscInt :: n
  PetscReal :: x, y
  PetscReal, dimension(1:n) :: dx, dy
  PetscReal, dimension(1:n) :: ProdRule

  PetscInt :: i

  !!! TODO: error checking on lengths?

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
function  d_prodrule_2(x, dx, y, dy, n)

!! Compute derivatives of x*y.
!! Assume x and y are functions and there are some n variables
!! t_1, t_2, ... 2_n, so that 
!!
!! x = x(t_1, .... t_n)
!! y = y(t_1, .... t_n)
!!
!! Assume we are given the partial derivatives in the arrays
!! dx and dy, i.e.,
!!
!! dx = [dx/dt_1,  .... , dx_dt_n]
!! dy = [dy/dt_1,  .... , dy_dt_n]
!!
!! Then we compute and return an array:
!!
!! [d(xy)/dt_1, ... , d(xy)/dt_n]

implicit none

  PetscInt :: n
  PetscReal :: x, y
  PetscReal, dimension(1:n) :: dx, dy
  PetscReal, dimension(1:n) :: d_prodrule_2

  PetscInt :: i

  !!! TODO: error checking on lengths?

  d_prodrule_2 = 0.d0

  do i = 1,n
    d_prodrule_2(i) = x*dy(i) + y*dx(i)  
  enddo

end function d_prodrule_2

! ************************************************************************** !

function  d_prodrule_3(x, dx, y, dy, z, dz, n)

implicit none

  PetscInt :: n
  PetscReal :: x, y, z
  PetscReal, dimension(1:n) :: dx, dy, dz
  PetscReal, dimension(1:n) :: d_yz
  PetscReal, dimension(1:n) :: d_prodrule_3
  PetscReal :: yz

  d_yz = d_prodrule_2(y, dy, z, dz, n)
  yz = y*z
  d_prodrule_3 = d_prodrule_2(x, dx, yz, d_yz, n)

end function  d_prodrule_3

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

!! derivatives of 1/x, see comments in above for general schema
!! of partial derivatives and etc.

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

function d_power(x, dx, alpha, n)

!! x^alpha

implicit none
  PetscInt :: n
  PetscReal :: x, alpha
  PetscReal, dimension(1:n) :: dx
  PetscReal, dimension(1:n) :: d_power

  PetscReal :: constpart
  PetscInt :: i

  constpart = alpha*(x**(alpha-1.d0))
  do i=1,n
    d_power(i) = dx(i)*constpart 
  end do

end function d_power

! ************************************************************************** !

end module Derivatives_utilities_module 
