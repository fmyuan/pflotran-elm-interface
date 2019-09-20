module Grid_Grdecl_Util_module

! A set of small utilities used by the grdecl and Output_Eclipse_module modules

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  ! Unit numbers for reading and writing reservoir engineering format files
  ! 50-59 are reserved for reservoir files
  PetscInt, parameter, public :: UNIT_GRDECL_READ = 50
  PetscInt, parameter, public :: UNIT_SPEC_WRITE  = 51
  PetscInt, parameter, public :: UNIT_SUMM_WRITE  = 52
  PetscInt, parameter, public :: UNIT_GRID_WRITE  = 53
  PetscInt, parameter, public :: UNIT_INIT_WRITE  = 54
  PetscInt, parameter, public :: UNIT_REST_WRITE  = 55

  public :: GetCorners, GetOtherDirections, GetIntersection
  public :: GetMDtoM2Conv, GetTriangleArea
  public :: GetM2toMDConv, qarea

  PetscReal, parameter :: e_atm = 1.01325

  contains

! *************************************************************************** !

subroutine GetCorners( ix, iy, iz, &
                       x000, x100, x010, x110, &
                       x001, x101, x011, x111, &
                       coord, zcorn, nx, ny )
  !
  ! Return the vertices of a hexahederal grid cell with location (ix,iy,iz)
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscInt , intent(in)  :: ix, iy, iz
  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3)
  PetscReal, intent(in) :: coord(:)
  PetscReal, intent(in) :: zcorn(:)
  PetscInt , intent(in)  :: nx, ny

  PetscInt  :: ibcx, ibcy, ibcz, ixo, iyo, ixp, iyp, &
               ibase, inx, iny, icl, icu, nxp
  PetscReal :: xl, yl, zl, xu, yu, zu, d0, d1

  nxp = nx+1

  ibcx = 2*(ix-1)
  ibcy = 2*(iy-1)
  ibcz = 2*(iz-1)

  ! Now loop over the four pillars of this cell

  do ixo = 1, 2
    do iyo = 1, 2

      ! Find pillar coordinates in the (nx+1).(ny+1)
      ! coord array (ix,ix+1),(iy,iy+1)

      ixp = ix+ixo-1
      iyp = iy+iyo-1

      ! Find base point in the coord data and extract the six pillar points

      ibase = 6*(nxp*(iyp-1)+(ixp-1))

      xl = coord(ibase+1)
      yl = coord(ibase+2)
      zl = coord(ibase+3)
      xu = coord(ibase+4)
      yu = coord(ibase+5)
      zu = coord(ibase+6)

      ! Find coordinates of the depths of the two corners on this pillar

      inx = ibcx + ixo
      iny = ibcy + iyo

      ! Find the coordinates of these corners in the depth array

      icl = GetLocationInDepthBlock(inx, iny, ibcz+1, nx, ny)
      icu = GetLocationInDepthBlock(inx, iny, ibcz+2, nx, ny)

      d0 = zcorn(icl)
      d1 = zcorn(icu)

      !  Fill in the cell corner locations

       if (ixo == 1) then
         if (iyo == 1) call fillGeoCorner(x000, x001, d0, d1, xl, yl, zl, &
                                                              xu, yu, zu)
         if (iyo == 2) call fillGeoCorner(x010, x011, d0, d1, xl, yl, zl, &
                                                              xu, yu, zu)
       endif
       if (ixo == 2) then
         if (iyo == 1) call fillGeoCorner(x100, x101, d0, d1, xl, yl, zl, &
                                                              xu, yu, zu)
         if (ixo == 2) call fillGeoCorner(x110, x111, d0, d1, xl, yl, zl, &
                                                              xu, yu, zu)
       endif

     enddo
  enddo

end subroutine GetCorners

! *************************************************************************** !

subroutine fillGeoCorner(x0, x1, d0, d1, xl, yl, zl, xu , yu, zu)
  !
  ! For coordinate line from xl,yl,zl to xu,yu,zu,
  ! find xyz locations at depths d0 and d1
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscReal, intent(out) :: x0(3), x1(3)
  PetscReal, intent(in ) :: d0, d1, xl, yl, zl, xu, yu, zu

  PetscReal :: f0, f1

  ! Fractions through the zl -> zu interval

  if (abs(zu-zl)>0.0) then
    f0 = (d0-zl)/(zu-zl)
    f1 = (d1-zl)/(zu-zl)
  else
    f0 = 0.5
    f1 = 0.5
  endif

  ! Fill in values (if d0 = zu, then f0 = 1 and x0->(xu,yu,zu)

  x0(1) = xl+f0*(xu-xl)
  x1(1) = xl+f1*(xu-xl)

  x0(2) = yl+f0*(yu-yl)
  x1(2) = yl+f1*(yu-yl)

  x0(3) = d0
  x1(3) = d1

end subroutine fillGeoCorner

! *************************************************************************** !

function GetLocationInDepthBlock(icx, icy, icz, nx, ny)
  !
  ! Get location of a node (icx,icy,icz) in the ZCORN 8.nx.ny.nzblock of depths
  !
  ! Author: Dave Ponting
  ! Date: 01/21/19

  implicit none

  PetscInt :: GetLocationInDepthBlock

  PetscInt, intent(in) :: icx, icy, icz, nx, ny
  PetscInt :: ncx, ncy, ncxy

  ncx  = 2*nx
  ncy  = 2*ny
  ncxy = ncx*ncy

  GetLocationInDepthBlock = ncxy*(icz-1)+ncx*(icy-1)+icx

end function GetLocationInDepthBlock

! *************************************************************************** !

function GetMDtoM2Conv()
  !
  ! Get conversion from mD to m2 (for handing Eclipse permeabilities)
  !
  ! Author: Dave Ponting
  ! Date: 01/21/19

  implicit none

  PetscReal :: GetMDtoM2Conv

  GetMDtoM2Conv = (1.0E-15)/e_atm

end function GetMDtoM2Conv

! *************************************************************************** !

function GetM2toMDConv()
  !
  ! Get conversion from m2 to mD (for handing Eclipse permeabilities)
  !
  ! Author: Dave Ponting
  ! Date: 01/21/19

  implicit none

  PetscReal :: GetM2toMDConv

  GetM2toMDConv = (1.0E+15)*e_atm

end function GetM2toMDConv

! *************************************************************************** !

subroutine GetIntersection(alx, aly, aux, auy, &
                           blx, bly, bux, buy, px, py, qinter)

  !
  ! Given two lines in two dimensions, one from (alx,aly) to (aux,auy) and
  ! one from (blx,bly) to (bux,buy), determine if they intersect and return
  ! intersection point (px,py) and flag indicating intersection
  !
  ! Author: Dave Ponting
  ! Date: 02/05/18

  implicit none

  PetscReal, intent(in)  :: alx, aly, aux, auy, blx, bly, bux, buy
  PetscReal, intent(out) :: px, py
  PetscBool, intent(out) :: qinter

  PetscReal :: apx, apy, bpx, bpy, den, t, p,r0, r1, &
               mai00, mai10, mai01, mai11, dt, dp
  PetscReal, parameter :: eps    = 1.0E-6
  PetscReal, parameter :: epseps = 1.0E-12

  ! Initialise returned values

  px     = 0.0
  py     = 0.0
  qinter = PETSC_FALSE

  ! Find diffences in x and y

  apx = aux - alx
  apy = auy - aly

  bpx = bux - blx
  bpy = buy - bly

  ! Check determinant - will be zero if parallel

  den = apx*bpy - apy*bpx

  if (abs(den) > 0.0) then

    ! Parameter values

    t = ( (aly-bly)*bpx - (alx-blx)*bpy )/den
    p = ( (blx-alx)*apy - (bly-aly)*apx )/den

    ! Residual errors

    r0 = blx - alx - t*apx + p*bpx
    r1 = bly - aly - t*apy + p*bpy

    ! If any residual left, do a cycle of iterative refinement

    if ( (abs(r0)>epseps) .or. (abs(r1)>epseps) ) then

      mai00 =  bpy/den
      mai01 = -bpx/den
      mai10 =  apy/den
      mai11 = -apx/den

      dt = mai00*r0 + mai01*r1
      dp = mai10*r0 + mai11*r1

      t=t+dt
      p=p+dp

    endif

    ! Check for intersection within (0,1)

    if (      (t >= -eps) .and. (t <= 1.0+eps) &
        .and. (p >= -eps) .and. (p <= 1.0+eps)) qinter = PETSC_TRUE

    ! Find intersection

    px  = alx + t*apx
    py  = aly + t*apy

  endif

end subroutine GetIntersection

! *************************************************************************** !

function GetTriangleArea(xa, ya, xb, yb, xc, yc)

  !
  ! Utility routine to get area of 2D triangle with corners a,b,c
  ! Uses the cross product expression
  !
  !         | 1  1  1  |
  ! A = 1/2*| xa xb xc | = 1/2*|xb.yc-xc.yb - xa.yc + xc.ya + xa.yb - xb.ya|
  !         | ya yb yc !
  !
  ! But (xa-xc).(yb-ya)-(xa-xb).(yc-ya) = xa.yb-xa.ya-xc.yb+xc.ya
  !                                      -xa.yc+xa.ya+xb.yc-xb.ya
  !                                     = xa.yb-xc.yb+xc.ya-xa.yc+xb.yc-xb.ya
  !                                     = xb.yc-xc.yb-xa.yc+xc.ya+xa.yb-xb.ya
  !
  ! Author: Dave Ponting
  ! Date: 02/05/18

  implicit none

  PetscReal GetTriangleArea
  PetscReal, intent(in) :: xa, ya, xb, yb, xc, yc

  GetTriangleArea = 0.5 * abs((xa-xc)*(yb-ya)-(xa-xb)*(yc-ya))

end function GetTriangleArea

! *************************************************************************** !

subroutine GetOtherDirections(idir, jdir, kdir)
  !
  ! Given three directions in cyclic order, retun the two following idir
  !
  ! Author: Dave Ponting
  ! Date: 02/05/18

  implicit none

  PetscInt, intent(in)  :: idir
  PetscInt, intent(out) :: jdir, kdir

  jdir = idir + 1
  if (jdir>3) jdir = 1
  kdir = jdir +1
  if (kdir>3) kdir = 1

end subroutine GetOtherDirections

! *************************************************************************** !

function qarea(f, idir)

  !
  ! Routine to get component id of the vector area of the quadrilateral f
  !
  ! Author: Dave Ponting
  ! Date: 02/05/18

  implicit none

  PetscReal qarea
  PetscReal, intent(in) :: f(0:1, 0:1, 3)
  PetscInt , intent(in) :: idir
  PetscReal :: x00, x01, x10, x11
  PetscReal :: y00, y01, y10, y11
  PetscReal :: at1, at2

  PetscInt :: jdir, kdir

  !  Find other directions

  call GetOtherDirections(idir, jdir, kdir)

  !  Get projection into the 2-plane orthogonal to e(idir)

  x00 = f(0, 0, jdir);y00 = f(0, 0, kdir)
  x10 = f(1, 0, jdir);y10 = f(1, 0, kdir)
  x01 = f(0, 1, jdir);y01 = f(0, 1, kdir)
  x11 = f(1, 1, jdir);y11 = f(1, 1, kdir)

  !  Obtain the area as two triangles

  at1 = GetTriangleArea(x00, y00, x01, y01, x11, y11)
  at2 = GetTriangleArea(x00, y00, x11, y11, x10, y10)

  !  Sum to get result

  qarea = at1+at2

end function qarea

! *************************************************************************** !

subroutine GetQuadSegment(iseg, idir, f, alx, aly, aux, auy)

  !
  ! Given a quadrilateral stored in terms of offsets io,jo
  ! return idir th projection of the iseg th segment from al to au
  !
  ! Author: Dave Ponting
  ! Date: 02/05/18

  implicit none

  PetscInt , intent(in) :: iseg, idir
  PetscReal, intent(in) :: f(0:1, 0:1, 3)
  PetscReal, intent(out) :: alx, aly, aux, auy
  PetscInt :: jdir, kdir

  ! Want idir direction projection, so find the two other transverse directions

  call GetOtherDirections(idir, jdir, kdir)

  ! Get segments going cyclically around the quadrilateral

  if (iseg == 1) then
    ! Segment 1 is (0,0) to (1,0)
    alx = f(0, 0, jdir);aly = f(0, 0 , kdir)
    aux = f(1, 0, jdir);auy = f(1, 0 , kdir)
  endif

  if (iseg == 2) then
    ! Segment 2 is (1,0) to (1,1)
    alx = f(1, 0, jdir);aly = f(1, 0, kdir)
    aux = f(1, 1, jdir);auy = f(1, 1, kdir)
  endif

  if (iseg == 3) then
    ! Segment 3 is (1,1) to (0,1)
    alx = f(1, 1, jdir);aly = f(1, 1, kdir)
    aux = f(0, 1, jdir);auy = f(0, 1, kdir)
  endif

  if (iseg == 4) then
    ! Segment 4 is (0,1) to (0,0)
    alx = f(0, 1, jdir);aly = f(0, 1, kdir)
    aux = f(0, 0, jdir);auy = f(0, 0, kdir)
  endif

end subroutine GetQuadSegment

! *************************************************************************** !

subroutine getCentre( x000, x100, x010, x110, &
                      x001, x101, x011, x111, centre )
  !
  ! Obtain the centre of a hexahederal grid cell
  !  with corners x000,..,x111
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscReal, intent(in) :: x000(3), x100(3), x010(3), x110(3), &
                           x001(3), x101(3), x011(3), x111(3)
  PetscReal,  intent(out) :: centre(3)

  PetscInt :: id

  do id = 1, 3
    centre(id) = 0.125*( x000(id) + x100(id) + x010(id) + x110(id) + &
                         x001(id) + x101(id) + x011(id) + x111(id) )
  enddo

end subroutine getCentre

end module Grid_Grdecl_Util_module
