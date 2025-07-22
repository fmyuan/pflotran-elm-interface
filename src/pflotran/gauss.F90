module Gauss_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Grid_Unstructured_Cell_module
  use PFLOTRAN_Constants_module

  implicit none

  private
  PetscInt, parameter, public :: LINE_TYPE          = 7

  type, public :: gauss_type
    PetscInt :: dim                 ! dimension
    PetscInt :: element_type        ! Element type
    PetscInt :: num_gauss_pts       ! Number of gauss points
    PetscReal, pointer :: r(:,:)    ! location of points
    PetscReal, pointer :: w(:)      ! weights
  end type gauss_type

  public :: GaussCalculatePoints, GaussDestroy, GaussInitialize

  contains

! ************************************************************************** !

subroutine GaussInitialize(gauss)
  !
  ! Initializes Gauss type
  !
  ! Author: Satish Karra, LANL
  ! Date: 6/19/2013
  !

  type(gauss_type) :: gauss

  gauss%dim = 0
  gauss%element_type = 0
  gauss%num_gauss_pts = 0
  nullify(gauss%r)
  nullify(gauss%w)

end subroutine GaussInitialize

! ************************************************************************** !

subroutine GaussCalculatePoints(gauss)
  !
  ! Calculates Gauss points
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !

  use Utility_module, only: DeallocateArray

  type(gauss_type) :: gauss
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)

  select case(gauss%dim)
    case(ONE_DIM_GRID)
      call Gauss1D(gauss%element_type,gauss%num_gauss_pts,r,w)
    case(TWO_DIM_GRID)
      call Gauss2D(gauss%element_type,gauss%num_gauss_pts,r,w)
    case(THREE_DIM_GRID)
      call Gauss3D(gauss%element_type,gauss%num_gauss_pts,r,w)
    case default
      print *, 'Error: Invalid dimension for Gauss point calculation'
      stop
  end select

  allocate(gauss%r(size(r,1),size(r,2)))
  allocate(gauss%w(size(w)))

  gauss%r = r
  gauss%w = w

  deallocate(r)
  deallocate(w)

end subroutine GaussCalculatePoints

! ************************************************************************** !

subroutine Gauss1D(element_type,num_gauss_pts,r,w)
  !
  ! Calculates Gauss points for 1D elements
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !

  PetscInt :: element_type
  PetscInt :: num_gauss_pts
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)

  allocate(r(num_gauss_pts,1))
  allocate(w(num_gauss_pts))

  if (element_type /= LINE_TYPE) then
    print *, 'Error: in Element type. Only L2 ' // &
             '(line type) can be used for 1D Gauss quadrature.'
  endif

  select case(num_gauss_pts)
    case(1)

    !------------------------------
    ! No of Gauss Points = 1
    !------------------------------

      r(1,1) = 0.d0
      !
      w(1) = 2.d0

    !------------------------------
    ! No of Gauss Points = 2
    !------------------------------

    case(2)

      r(1,1) = -1.d0/sqrt(3.d0)
      r(2,1) = -r(1,1)
      !
      w(1) = 1.d0
      w(2) = 1.d0

    !-------------------------------
    ! No of Gauss Points = 3
    !-------------------------------

    case(3)

      r(1,1) = -sqrt(0.6d0)
      r(2,1) = 0.d0
      r(3,1) = -r(1,1)
      !
      w(1) = 5.d0/9.d0
      w(2) = 8.d0/9.d0
      w(3) = w(1)

    !--------------------------------
    ! No of Gauss Points = 4
    !--------------------------------

    case(4)

      r(1,1) = -0.861136311594053d0
      r(2,1) = -0.339981043584856d0
      r(3,1) =  0.339981043584856d0
      r(4,1) =  0.861136311594053d0
      !
      w(1) = 0.347854845137454d0
      w(2) = 0.652145154862546d0
      w(3) = 0.652145154862546d0
      w(4) = 0.347854845137454d0

    !----------------------------------
    ! No of Gauss Points = 5
    !----------------------------------

    case(5)

      r(1,1) = -0.906179845938664d0
      r(2,1) = -0.538469310105683d0
      r(3,1) =  0.000000000000000d0
      r(4,1) =  0.538469310105683d0
      r(5,1) =  0.906179845938664d0
      !
      w(1) =  0.236926885056189d0
      w(2) =  0.478628670499366d0
      w(3) =  0.568888888888889d0
      w(4) =  0.478628670499366d0
      w(5) =  0.236926885056189d0

    !----------------------------------
    ! No of Gauss Points = 6
    !----------------------------------

    case(6)

      r(1,1) = -0.932469514203152d0
      r(2,1) = -0.661209386466265d0
      r(3,1) = -0.238619186083197d0
      r(4,1) =  0.238619186083197d0
      r(5,1) =  0.661209386466265d0
      r(6,1) =  0.932469514203152d0
      !
      w(1) =  0.171324492379170d0
      w(2) =  0.360761573048139d0
      w(3) =  0.467913934572691d0
      w(4) =  0.467913934572691d0
      w(5) =  0.360761573048139d0
      w(6) =  0.171324492379170d0

    !------------------------------------
    ! No of Gauss Points = 7
    !------------------------------------

    case(7)

      r(1,1) = -0.949107912342759d0
      r(2,1) = -0.741531185599394d0
      r(3,1) = -0.405845151377397d0
      r(4,1) =  0.000000000000000d0
      r(5,1) =  0.405845151377397d0
      r(6,1) =  0.741531185599394d0
      r(7,1) =  0.949107912342759d0
      !
      w(1) =  0.129484966168870d0
      w(2) =  0.279705391489277d0
      w(3) =  0.381830050505119d0
      w(4) =  0.417959183673469d0
      w(5) =  0.381830050505119d0
      w(6) =  0.279705391489277d0
      w(7) =  0.129484966168870d0

    !------------------------------------
    ! No of Gauss Points = 8
    !------------------------------------

    case(8)

      r(1,1) = -0.960289856497536d0
      r(2,1) = -0.796666477413627d0
      r(3,1) = -0.525532409916329d0
      r(4,1) = -0.183434642495650d0
      r(5,1) =  0.183434642495650d0
      r(6,1) =  0.525532409916329d0
      r(7,1) =  0.796666477413627d0
      r(8,1) =  0.960289856497536d0
      !
      w(1) =  0.101228536290376d0
      w(2) =  0.222381034453374d0
      w(3) =  0.313706645877887d0
      w(4) =  0.362683783378362d0
      w(5) =  0.362683783378362d0
      w(6) =  0.313706645877887d0
      w(7) =  0.222381034453374d0
      w(8) =  0.101228536290376d0

    !------------------------------------
    ! No of Gauss Points = 9
    !------------------------------------

    case(9)

      r(1,1) = -0.968160239507626d0
      r(2,1) = -0.836031170326636d0
      r(3,1) = -0.613371432700590d0
      r(4,1) = -0.324253423403809d0
      r(5,1) =  0.000000000000000d0
      r(6,1) =  0.324253423403809d0
      r(7,1) =  0.613371432700590d0
      r(8,1) =  0.836031107326636d0
      r(9,1) =  0.968160239507626d0

      w(1) =  0.081274388361574d0
      w(2) =  0.180648160694857d0
      w(3) =  0.260610696402935d0
      w(4) =  0.312347077040003d0
      w(5) =  0.330239355001260d0
      w(6) =  0.312347077040003d0
      w(7) =  0.260610696402935d0
      w(8) =  0.180648160694857d0
      w(9) =  0.081274388361574d0

   case default
     print *, 'Error in num_gauss_pts for 1D Gauss quadrature'
     stop
   end select

end subroutine Gauss1D

! ************************************************************************** !

subroutine Gauss2D(element_type,num_gauss_pts,r,w)
  !
  ! Calculates Gauss points for 2D elements
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !

  PetscInt :: element_type
  PetscInt :: num_gauss_pts
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)

  select case(element_type)
    case(QUAD_TYPE)
      call GaussSquare(num_gauss_pts,r,w)
    case(TRI_TYPE)
      call GaussTriangle(num_gauss_pts,r,w)
    case default
      print *, 'Error: Only T3 and Q4 elements available for 2D.'
      stop
  end select

end subroutine Gauss2D

! ************************************************************************** !

subroutine GaussSquare(num_gauss_pts,r,w)
  !
  ! Calculates Gauss points for Q4 element
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !

  PetscInt :: num_gauss_pts
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  PetscReal, pointer :: l(:,:)
  PetscReal, pointer :: m(:)
  PetscInt :: counter,i,j

  allocate(r(num_gauss_pts*num_gauss_pts,2))
  allocate(w(num_gauss_pts*num_gauss_pts))
  allocate(l(num_gauss_pts,1))
  allocate(m(num_gauss_pts))

  ! 1D Gauss points are stored in l vector and weights are stored in m vector

  call Gauss1D(LINE_TYPE,num_gauss_pts,l,m)

  ! Generate the Q4 Gauss points and weights using for loops
  counter = 1
  do i = 1, num_gauss_pts
    do j = 1, num_gauss_pts
      r(counter,1) = l(i,1)
      r(counter,2) = l(j,1)
      w(counter) = m(i)*m(j)
      counter = counter + 1
    enddo
  enddo

  deallocate(l)
  deallocate(m)

end subroutine GaussSquare

! ************************************************************************** !

subroutine GaussTriangle(num_gauss_pts,r,w)
  !
  ! Calculates Gauss points for T3 elements
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !

  PetscInt :: num_gauss_pts
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)

  allocate(r(num_gauss_pts,2))
  allocate(w(num_gauss_pts))

  select case(num_gauss_pts)
    case(1)

    !------------------------------
    ! No of Gauss Points = 1
    !------------------------------
      r(1,:) = 1.d0/3.d0*(/1.d0,1.d0/)
      !
      w(1) = 0.5d0

    !-------------------------------
    ! No of Gauss Points = 3
    !-------------------------------

    case(3)

      r(1,:) = 0.5d0*(/1.d0,1.d0/)
      r(2,:) = 0.5d0*(/1.d0,0.d0/)
      r(3,:) = 0.5d0*(/0.d0,1.d0/)
      !
      w = 1.d0/6.d0*(/1.d0,1.d0,1.d0/)

    !--------------------------------
    ! No of Gauss Points = 4
    !--------------------------------

    case(4)

      r(1,:) = (/0.333333333333333d0,0.333333333333333d0/)
      r(2,:) = (/0.600000000000000d0,0.200000000000000d0/)
      r(3,:) = (/0.200000000000000d0,0.600000000000000d0/)
      r(4,:) = (/0.200000000000000d0,0.200000000000000d0/)
      !
      w = 0.5d0*(/-0.562500000000000d0, &
                 0.520833333333333d0, &
                 0.520833333333333d0, &
                 0.520833333333333d0/)

    !------------------------------------
    ! No of Gauss Points = 7
    !------------------------------------

    case(7)

      r(1,:) = (/0.333333333333333d0,0.333333333333333d0/)
      r(2,:) = (/0.797426985353087d0,0.101286507323456d0/)
      r(3,:) = (/0.101286507323456d0,0.797426985353087d0/)
      r(4,:) = (/0.101286507323456d0,0.101286507323456d0/)
      r(5,:) = (/0.470142064105115d0,0.059715871789770d0/)
      r(6,:) = (/0.059715871789770d0,0.470142064105115d0/)
      r(7,:) = (/0.470142064105115d0,0.470142064105115d0/)
      !
      w = 0.5d0*(/0.225000000000000d0, &
                0.125939180544827d0, &
                0.125939180544827d0, &
                0.125939180544827d0, &
                0.132394152788506d0, &
                0.132394152788506d0, &
                0.132394152788506d0/)

    !------------------------------------
    ! No of Gauss Points = 9
    !------------------------------------

    case(9)

      r(1,:) = (/0.124949503233232d0,0.437525248383384d0/)
      r(2,:) = (/0.437525248383384d0,0.124949503233232d0/)
      r(3,:) = (/0.437525248383384d0,0.437525248383384d0/)
      r(4,:) = (/0.797112651860071d0,0.165409927389841d0/)
      r(5,:) = (/0.797112651860071d0,0.037477420750088d0/)
      r(6,:) = (/0.165409927389841d0,0.797112651860071d0/)
      r(7,:) = (/0.165409927389841d0,0.037477420750088d0/)
      r(8,:) = (/0.037477420750088d0,0.797112651860071d0/)
      r(9,:) = (/0.037477420750088d0,0.165409927389841d0/)

      w = 0.5d0*(/0.205950504760887d0, &
                0.205950504760887d0, &
                0.205950504760887d0, &
                0.063691414286223d0, &
                0.063691414286223d0, &
                0.063691414286223d0, &
                0.063691414286223d0, &
                0.063691414286223d0, &
                0.063691414286223d0/)

    !------------------------------------
    ! No of Gauss Points = 12
    !------------------------------------

    case(12)

      r(1,:) = (/0.873821971016996d0,0.063089014491502d0/)
      r(2,:) = (/0.063089014491502d0,0.873821971016996d0/)
      r(3,:) = (/0.063089014491502d0,0.063089014491502d0/)
      r(4,:) = (/0.501426509658179d0,0.249286745170910d0/)
      r(5,:) = (/0.249286745170910d0,0.501426509658179d0/)
      r(6,:) = (/0.249286745170910d0,0.249286745170910d0/)
      r(7,:) = (/0.636502499121399d0,0.310352451033785d0/)
      r(8,:) = (/0.636502499121399d0,0.053145049844816d0/)
      r(9,:) = (/0.310352451033785d0,0.636502499121399d0/)
      r(10,:) = (/0.310352451033785d0,0.053145049844816d0/)
      r(11,:) = (/0.053145049844816d0,0.636502499121399d0/)
      r(12,:) = (/0.053145049844816d0,0.310352451033785d0/)

      w = 0.5d0*(/0.050844906370207d0, &
                0.050844906370207d0, &
                0.050844906370207d0, &
                0.116786275726379d0, &
                0.116786275726379d0, &
                0.116786275726379d0, &
                0.082851075618374d0, &
                0.082851075618374d0, &
                0.082851075618374d0, &
                0.082851075618374d0, &
                0.082851075618374d0, &
                0.082851075618374d0/)

    !------------------------------------
    ! No of Gauss Points = 13
    !------------------------------------

    case(13)

      r(1,:) = (/0.333333333333333d0,0.333333333333333d0/)
      r(2,:) = (/0.479308067841923d0,0.260345966079038d0/)
      r(3,:) = (/0.260345966079038d0,0.479308067841923d0/)
      r(4,:) = (/0.260345966079038d0,0.260345966079038d0/)
      r(5,:) = (/0.869739794195568d0,0.065130102902216d0/)
      r(6,:) = (/0.065130102902216d0,0.869739794195568d0/)
      r(7,:) = (/0.065130102902216d0,0.065130102902216d0/)
      r(8,:) = (/0.638444188569809d0,0.312865496004875d0/)
      r(9,:) = (/0.638444188569809d0,0.086903154253160d0/)
      r(10,:) = (/0.312865496004875d0,0.638444188569809d0/)
      r(11,:) = (/0.312865496004875d0,0.086903154253160d0/)
      r(12,:) = (/0.086903154253160d0,0.638444188569809d0/)
      r(13,:) = (/0.086903154253160d0,0.312865496004875d0/)

      w = 0.5d0*(/-0.149570044467670d0, &
                0.175615257433204d0, &
                0.175615257433204d0, &
                0.175615257433204d0, &
                0.053347235608839d0, &
                0.053347235608839d0, &
                0.053347235608839d0, &
                0.077113760890257d0, &
                0.077113760890257d0, &
                0.077113760890257d0, &
                0.077113760890257d0, &
                0.077113760890257d0, &
                0.077113760890257d0/)

   case default
     print *, 'Invalid num_gauss_pts for T3 Gauss quadrature'
     stop
   end select


end subroutine GaussTriangle

! ************************************************************************** !

subroutine GaussTetrahedra(num_gauss_pts,r,w)
  !
  ! Calculates Gauss points for tetrahedra elements
  !
  ! Author: Satish Karra, LANL
  ! Date: 7/11/2013
  !

  PetscInt :: num_gauss_pts
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  PetscReal :: p1,p2,p3,p4,p5,p6,p7
  PetscReal :: q1,q2,q3,q4

  allocate(r(num_gauss_pts,3))
  allocate(w(num_gauss_pts))

  select case(num_gauss_pts)
    case(1)

    !------------------------------
    ! No of Gauss Points = 1
    !------------------------------
      r(1,:) = 1.d0/4.d0*(/1.d0,1.d0,1.d0/)
      !
      w(1) = 1.d0/6.d0

    !--------------------------------
    ! No of Gauss Points = 4
    !--------------------------------

    case(4)

      p1 = 0.5854101966249638d0
      p2 = 0.1381966011250105d0
      r(1,:) = (/p1,p2,p2/)
      r(2,:) = (/p2,p1,p2/)
      r(3,:) = (/p2,p2,p1/)
      r(4,:) = (/p2,p2,p2/)
      !
      w = 1.0/(6.d0*4.d0)*(/1.d0,1.d0,1.d0,1.d0/)

    !------------------------------------
    ! No of Gauss Points = 5
    !------------------------------------

    case(5)

      r(1,:) = 1.d0/4.d0*(/1.d0,1.d0,1.d0/)
      r(2,:) = (/1.d0/2.d0,1.d0/6.d0,1.d0/6.d0/)
      r(3,:) = (/1.d0/6.d0,1.d0/2.d0,1.d0/6.d0/)
      r(4,:) = (/1.d0/6.d0,1.d0/6.d0,1.d0/2.d0/)
      r(5,:) = (/1.d0/6.d0,1.d0/6.d0,1.d0/6.d0/)
      !
      w = 1.d0/6.d0*(/-4.d0/5.d0,9.d0/20.d0,9.d0/20.d0,9.d0/20.d0,9.d0/20.d0/)

    !------------------------------------
    ! No of Gauss Points = 11
    !------------------------------------

    case(11)

      p1 = 0.250000000000000d0
      p2 = 0.785714285714286d0
      p3 = 0.071428571428571d0
      p4 = 0.399403576166799d0
      p5 = 0.100596423833201d0

      r(1,:) = (/p1,p1,p1/)
      r(2,:) = (/p2,p3,p3/)
      r(3,:) = (/p3,p2,p3/)
      r(4,:) = (/p3,p3,p2/)
      r(5,:) = (/p3,p3,p3/)
      r(6,:)  = (/p4,p5,p5/)
      r(7,:)  = (/p5,p4,p5/)
      r(8,:)  = (/p5,p5,p4/)
      r(9,:)  = (/p5,p4,p4/)
      r(10,:) = (/p4,p5,p4/)
      r(11,:) = (/p4,p4,p5/)
      ! Gauss weights
      q1 = -0.013155555555556d0
      q2 =  0.007622222222222d0
      q3 =  0.024888888888889d0

      w = (/q1,q2,q2,q2,q2,q3,q3,q3,q3,q3,q3/)


    !------------------------------------
    ! No of Gauss Points = 15
    !------------------------------------

    case(15)

      p1 = 0.250000000000000d0
      p2 = 0.000000000000000d0
      p3 = 0.333333333333333d0
      p4 = 0.727272727272727d0
      p5 = 0.090909090909091d0
      p6 = 0.066550153573664d0
      p7 = 0.433449846426336d0

      r(1,:) = (/p1,p1,p1/)
      r(2,:) = (/p2,p3,p3/)
      r(3,:) = (/p3,p2,p3/)
      r(4,:) = (/p3,p3,p2/)
      r(5,:) = (/p3,p3,p3/)
      r(6,:) = (/p4,p5,p5/)
      r(7,:) = (/p5,p4,p5/)
      r(8,:) = (/p5,p5,p4/)
      r(9,:) = (/p5,p5,p5/)
      r(10,:) = (/p6,p7,p7/)
      r(11,:) = (/p7,p6,p7/)
      r(12,:) = (/p7,p7,p6/)
      r(13,:) = (/p7,p6,p6/)
      r(14,:) = (/p6,p7,p6/)
      r(15,:) = (/p6,p6,p7/)
      !
      q1 = 0.030283678097089d0
      q2 = 0.006026785714286d0
      q3 = 0.011645249086029d0
      q4 = 0.010949141561386d0

      w = (/q1,q2,q2,q2, q2,q3,q3,q3,q3,q4,q4,q4,q4,q4,q4/)

    case default
     print *, 'Invalid num_gauss_pts for Tetrahedra Gauss quadrature'
     stop
   end select


end subroutine GaussTetrahedra

! ************************************************************************** !

subroutine GaussPyramid(num_gauss_pts,r,w)
  !
  ! Calculates Gauss points for tetrahedra elements
  ! Reference:
  ! http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_pyramid/quadrature_rules_pyramid.html
  !
  ! Author: Satish Karra, LANL
  ! Date: 7/11/2013
  !
  PetscInt :: num_gauss_pts
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)

  allocate(r(num_gauss_pts,3))
  allocate(w(num_gauss_pts))

  select case(num_gauss_pts)
    case(1)

    !------------------------------
    ! No of Gauss Points = 1
    !------------------------------
      r(1,:) = (/0.d0,0.d0,0.25d0/)
      !
      w(1) = 1.d0

    !--------------------------------
    ! No of Gauss Points = 5
    !--------------------------------

    case(5)


      r(1,:) = (/-0.48686449556014765641d0,-0.48686449556014765641d0,0.16666666666666666667d0/)
      r(2,:) = (/ 0.48686449556014765641d0,-0.48686449556014765641d0,0.16666666666666666667d0/)
      r(3,:) = (/ 0.48686449556014765641d0, 0.48686449556014765641d0,0.16666666666666666667d0/)
      r(4,:) = (/-0.48686449556014765641d0, 0.48686449556014765641d0,0.16666666666666666667d0/)
      r(5,:) = (/ 0.00000000000000000000d0, 0.00000000000000000000d0,0.70000000000000000000d0/)
      !
      w = (/0.2109375d0, &
            0.2109375d0, &
            0.2109375d0, &
            0.2109375d0, &
            0.15625d0/)

    !------------------------------------
    ! No of Gauss Points = 6
    !------------------------------------

    case(6)

      r(1,:) = (/-0.48795003647426658968d0,-0.48795003647426658968d0,0.16666666666666666667d0/)
      r(2,:) = (/ 0.48795003647426658968d0,-0.48795003647426658968d0,0.16666666666666666667d0/)
      r(3,:) = (/ 0.48795003647426658968d0, 0.48795003647426658968d0,0.16666666666666666667d0/)
      r(4,:) = (/-0.48795003647426658968d0, 0.48795003647426658968d0,0.16666666666666666667d0/)
      r(5,:) = (/ 0.00000000000000000000d0, 0.00000000000000000000d0,0.58333333333333333333d0/)
      r(6,:) = (/ 0.00000000000000000000d0, 0.00000000000000000000d0,0.75000000000000000000d0/)
      !
      w = (/0.21000000000000000000d0, &
            0.21000000000000000000d0, &
            0.21000000000000000000d0, &
            0.21000000000000000000d0, &
            0.21000000000000000000d0, &
            0.06000000000000000000d0, &
            0.10000000000000000000d0/)

     case default
     print *, 'Invalid num_gauss_pts for Tetrahedra Gauss quadrature'
     stop
   end select


end subroutine GaussPyramid

! ************************************************************************** !

subroutine Gauss3D(element_type,num_gauss_pts,r,w)
  !
  ! Calculates Gauss points for 3D element
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !

  PetscInt :: element_type
  PetscInt :: num_gauss_pts
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)

  select case(element_type)
    case(HEX_TYPE)
      call GaussBrick(num_gauss_pts,r,w)
    case(WEDGE_TYPE)
      call GaussWedge(num_gauss_pts,r,w)
    case(TET_TYPE)
      call GaussTetrahedra(num_gauss_pts,r,w)
    case(PYR_TYPE)
      call GaussPyramid(num_gauss_pts,r,w)
    case default
      print *, 'Error: Only B8, W6, P5 and TET4 elements available for 3D.'
      stop
  end select

end subroutine Gauss3D

! ************************************************************************** !

subroutine GaussBrick(num_gauss_pts,r,w)
  !
  ! Calculates Gauss points for B8 element
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !

  PetscInt :: num_gauss_pts
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  PetscReal, pointer :: l(:,:)
  PetscReal, pointer :: m(:)
  PetscInt :: counter, i, j, k

  allocate(r(num_gauss_pts*num_gauss_pts*num_gauss_pts,3))
  allocate(w(num_gauss_pts*num_gauss_pts*num_gauss_pts))

  call Gauss1D(LINE_TYPE,num_gauss_pts,l,m)

  ! Generate the B8 Gauss points and weights using for loops

  counter = 1
  do i = 1, num_gauss_pts
    do j = 1, num_gauss_pts
      do k = 1, num_gauss_pts
        r(counter,1) = l(i,1)
        r(counter,2) = l(j,1)
        r(counter,3) = l(k,1)
        w(counter) = m(i)*m(j)*m(k)
        counter = counter + 1
      enddo
    enddo
  enddo

  deallocate(l)
  deallocate(m)

end subroutine GaussBrick

! ************************************************************************** !

subroutine GaussWedge(num_gauss_pts,r,w)
  !
  ! Calculates Gauss points for wedge element
  !
  ! Author: Satish Karra, LANL
  ! Date: 7/10//2013
  !

  PetscInt :: num_gauss_pts
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  PetscReal, pointer :: rT3(:,:),rL2(:,:)
  PetscReal, pointer :: wT3(:),wL2(:)
  PetscInt :: i, j

  allocate(r(num_gauss_pts*num_gauss_pts,3))
  allocate(w(num_gauss_pts*num_gauss_pts))

  call Gauss1D(LINE_TYPE,num_gauss_pts,rL2,wL2)
  call Gauss2D(TRI_TYPE,num_gauss_pts,rT3,wT3)

  ! Generate the wedge Gauss points and weights using for loops
  do i = 1, num_gauss_pts
    do j = 1, num_gauss_pts
      r((i-1)*num_gauss_pts+j,1) = rT3(i,1)
      r((i-1)*num_gauss_pts+j,2) = rT3(i,2)
      r((i-1)*num_gauss_pts+j,3) = rL2(j,1)
      w((i-1)*num_gauss_pts+j) = wT3(i)*wL2(j)
    enddo
  enddo

  deallocate(rL2)
  deallocate(rT3)
  deallocate(wL2)
  deallocate(wT3)

end subroutine GaussWedge

! ************************************************************************** !

subroutine GaussDestroy(gauss)
  !
  ! Deallocate gauss type
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !

  type(gauss_type) :: gauss

  deallocate(gauss%r)
  nullify(gauss%r)
  deallocate(gauss%w)
  nullify(gauss%w)

end subroutine GaussDestroy

end module Gauss_module
