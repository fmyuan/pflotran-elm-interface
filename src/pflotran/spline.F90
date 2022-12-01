module Spline_module

#include "petsc/finclude/petscsys.h"

  public

contains

! ************************************************************************** !

subroutine SplineSecondDeriv(t,s,n,s_dp)

  !
  !     s_dp(1:n), s double prime, is the second derivative of the
  !     interpolating function at the tabluated point t(n),
  !     where s(1:n) = function(t(1:n))
  !
  !     Reference: press, w.h., b.p. flannery, s.a. teukolsky, and
  !     w.t. vetterling.
  !     1986.  numerical recipes, the art of scientific computing,
  !     cambridge university press, cambridge.  pp. 86-89.
  !
  !     Modified by Heeho Park
  !     Date: 11/10/2023

  use PFLOTRAN_Constants_module

  implicit none

  PetscInt :: i,n,k
  PetscReal :: t(n),s(n),s_dp(n),temp_array(n)
  PetscReal :: alpha,temp,bound1,bound2

  s_dp(1) = 0.d0
  temp_array(1) = 0.d0
  
  do i = 2,n-1
    alpha = (t(i)-t(i-1))/(t(i+1)-t(i-1))
    temp = alpha*s_dp(i-1)+2.d0
    s_dp(i) = (alpha-1.d0)/temp
    temp_array(i) = (6.d0*((s(i+1)-s(i))/(t(i+1)-t(i)) - &
                    (s(i)-s(i-1))/(t(i)-t(i-1)))/ &
                    (t(i+1)-t(i-1)) - alpha*temp_array(i-1))/temp
  enddo
  
  bound1 = 0.d0
  bound2 = 0.d0
  s_dp(n) = (bound2-bound1*temp_array(n-1))/(bound1*s_dp(n-1)+1.d0)

  do k = n-1,1,-1
    s_dp(k) = s_dp(k)*s_dp(k+1)+temp_array(k)
  enddo

  return
  
end subroutine SplineSecondDeriv

! ************************************************************************** !

subroutine SplineInterp(arr1,arr2,s_dp,n,t,s_int)

  !     two arrays arr1(1:n) and arr2(1:n) of length n and the second
  !     derivative from SplineSecondDeriv s_dp(1:n), give a value of t,
  !     this returns cubic-spline interpolated value of s_int. 
  !
  !     Reference: press, w.h., b.p. flannery, s.a. teukolsky, and
  !     w.t. vetterling.
  !     1986.  numerical recipes, the art of scientific computing,
  !     cambridge university press, cambridge.  pp. 86-89.
  !
  !     Modified by Heeho Park
  !     Date: 11/10/2023

  implicit none

  PetscInt :: n,k,lower,upper
  PetscReal :: arr1(n),arr2(n),s_dp(n)
  PetscReal :: denom,c1,c2,t,s_int

  lower = 1
  upper = n
  
  do while (upper-lower > 1)
    k = (upper+lower)/2
    if (arr1(k) > t) then
      upper = k
    else
      lower = k
    endif
  end do
  
  denom = arr1(upper)-arr1(lower)
  c1 = (arr1(upper)-t)/denom
  c2 = (t-arr1(lower))/denom
  s_int = c1*arr2(lower)+c2*arr2(upper)+ &
          ((c1**3-c1)*s_dp(lower)+(c2**3-c2)*s_dp(upper))*(denom**2)/6.d0

  return

end subroutine SplineInterp

! ************************************************************************** !

subroutine BisectionSearch(arr,n,val,ind)

  !     given an array arr of length n, and given a value val, returns a
  !     value ind such that val is between arr(ind) and arr(ind+1).
  !     arr must be monotonic, either increasing or decreasing.
  !     ind=0 or ind=n is returned to indicate that val is out of range.
  !
  !     Reference press, w.h., b.p. flannery, s.a. teukolsky, and
  !     w.t. vetterling.
  !     1986.  numerical recipes, the art of scientific computing.
  !     cambridge university press, cambridge.
  !
  !     Modified by Heeho Park
  !     Date: 11/10/2023

  implicit none

  PetscInt :: lower,upper,mid,ind,n
  PetscReal :: arr(n)
  PetscReal :: val

  lower = 0
  upper = n+1
  do while (upper-lower > 1)
    mid = (upper+lower)/2
    if ((arr(n) > arr(1)).eqv.(val > arr(mid))) then
      lower = mid
    else
      upper = mid
    endif
  end do
  ind = lower

  return
  
end subroutine BisectionSearch

end module Spline_module
