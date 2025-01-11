module Grid_Eclipse_Util_module

! A set of small utilities used by the grdecl and Output_Eclipse_module modules

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  public :: GetCorners, GetOtherDirections, GetIntersection
  public :: GetMDtoM2Conv,GetETtoTrConv, GetTriangleArea
  public :: GetM2toMDConv,GetTrtoETConv, QArea, GetCentre, GetQuadSegment
  public :: SetMapAxes, GetMapAxes, SetMapAxesUnits, GetMapAxesUnits, &
            SetGDORIENT, GetGDORIENT
  public :: BroadcastInt, BroadcastIntN ,BroadcastRealN, BroadcastCharN
  public :: ParallelSendToIO, ParallelRecvFromProc
  public :: SetProbDims, GetProbDims, &
            SetDefaultLocalsPerGlobal, MakeLocalBasePointers, &
            SetNLGR, SetLGRName, SetLGRDimensions, GetLGRLocation, GetIsLGR, GetLGRName8, &
            GetNLGR, GetLGRDimensions, GetNRepReg, SplitRptreg, &
            GetFipaName, GetBaseFip, GetNinFip
  public :: GetDRPSubindex, GetDRPSubindexName, AllocArgs, DeallocArgs
  public :: StoreCGL, GetNCGL ,FindCGL, ReleaseCGL
  public :: CopyString, PrintTidy8, PrintInt8

  ! Index of the global grid and max number of lgrs
  PetscInt, public, parameter :: g_ifld = 1
  PetscInt, public, parameter :: g_mlgr = 10

  ! Number of directional rel perm num arrays and pointers to contents

  PetscInt, public , parameter :: idrpa_krx = 1
  PetscInt, public , parameter :: idrpa_kry = 2
  PetscInt, public , parameter :: idrpa_krz = 3
  PetscInt, public , parameter :: idrpa_imx = 4
  PetscInt, public , parameter :: idrpa_imy = 5
  PetscInt, public , parameter :: idrpa_imz = 6
  PetscInt, public , parameter :: mdrpa     = 6
  PetscInt, public , parameter :: mfipa     = 6

  ! Up and down orientations

  PetscInt, public , parameter :: g_npon   = 2
  PetscInt, public , parameter :: g_ipon_l = 1
  PetscInt, public , parameter :: g_ipon_u = 2

  ! Cartesian frame pointers

  PetscInt, public , parameter :: g_xdir = 1
  PetscInt, public , parameter :: g_ydir = 2
  PetscInt, public , parameter :: g_zdir = 3

  PetscInt, public , parameter :: g_ndir = 3

  PetscInt, parameter, public :: UNIT_GRDECL_READ = 50
  PetscInt, parameter, public :: UNIT_EGRID_READ  = 60

  PetscReal, parameter :: e_atm = 1.01325

  private

  PetscReal :: g_mapaxes(6) = 0.0
  PetscBool :: g_mapaxes_set = PETSC_FALSE
  character(len=MAXWORDLENGTH) :: g_mapaxesunits = 'METRES'

  PetscBool :: g_setgdorient    = PETSC_FALSE
  PetscBool :: g_islefthandgrid = PETSC_FALSE

  PetscInt :: gu_nx = 1
  PetscInt :: gu_ny = 1
  PetscInt :: gu_nz = 1

  PetscInt :: gu_nlgr = 1
  character(len=MAXWORDLENGTH) :: gu_zlgr(g_mlgr) = ""
  PetscInt :: gu_lgrd(3,g_mlgr) = 1

  PetscInt :: g_total_fip = 1

  PetscInt :: g_b_fipreg(mfipa)
  PetscInt :: g_n_fipreg(mfipa)

  PetscInt :: g_ncgl  = 0
  PetscInt :: g_ntclg = 0
  PetscInt, allocatable :: g_cg(:)
  PetscInt, allocatable :: g_cl(:)
  PetscInt, allocatable :: g_io(:)
  PetscInt :: g_ibclg(g_mlgr) = 0
  PetscInt :: g_ndclg(g_mlgr) = 0

  character(len=MAXWORDLENGTH) :: g_fipname(mfipa) = ''

  contains

! *************************************************************************** !

subroutine GetCorners(ix, iy, iz, &
                       x000, x100, x010, x110, &
                       x001, x101, x011, x111, &
                       coord, zcorn, nx, ny)
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
         if (iyo == 1) call FillGeoCorner(x000, x001, d0, d1, xl, yl, zl, &
                                                              xu, yu, zu)
         if (iyo == 2) call FillGeoCorner(x010, x011, d0, d1, xl, yl, zl, &
                                                              xu, yu, zu)
       endif
       if (ixo == 2) then
         if (iyo == 1) call FillGeoCorner(x100, x101, d0, d1, xl, yl, zl, &
                                                              xu, yu, zu)
         if (iyo == 2) call FillGeoCorner(x110, x111, d0, d1, xl, yl, zl, &
                                                              xu, yu, zu)
       endif

     enddo
  enddo

end subroutine GetCorners

! *************************************************************************** !

subroutine FillGeoCorner(x0, x1, d0, d1, xl, yl, zl, xu , yu, zu)
  !
  ! For coordinate line from xl,yl,zl to xu,yu,zu,
  ! find xyz locations at depths d0 and d1
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscReal, intent(out) :: x0(3), x1(3)
  PetscReal, intent(in) :: d0, d1, xl, yl, zl, xu, yu, zu

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

end subroutine FillGeoCorner

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

function GetETtoTrConv()
  !
  ! Get conversion from cp.rm3/day/bar to PaS.rm3/sec/Pa
  !
  ! Author: Dave Ponting
  ! Date: 11/24/20

  implicit none

  PetscReal :: GetETtoTrConv

  GetETtoTrConv = 0.001/(3600.0*24.0*1.0D5)

end function GetETtoTrConv

! *************************************************************************** !

function GetTrtoETConv()
  !
  ! Get conversion from PaS.rm3/sec/Pa to cp.rm3/day/bar
  !
  ! Author: Dave Ponting
  ! Date: 11/24/20

  implicit none

  PetscReal :: GetTrtoETConv

  GetTrtoETConv = 1.0/GetETtoTrConv()

end function GetTrtoETConv

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

  PetscReal :: apx, apy, bpx, bpy, den, t, p, pxc, pyc
  PetscReal :: eps = 0.001

  !  Initialise returned values

  px     = 0.0
  py     = 0.0
  qinter = PETSC_FALSE

  !  Find gradients wrt x and y

  apx = aux - alx
  apy = auy - aly

  bpx = bux - blx
  bpy = buy - bly

  !  Check determinant - will be zero if parallel

  den = apx*bpy - apy*bpx

  if (abs(den) > 0.0) then

    !  Parameter values

    t = ((aly-bly)*bpx - (alx-blx)*bpy)/den
    p = ((blx-alx)*apy - (bly-aly)*apx)/den

    !  Check for intersection within (0,1)

    if ((t >= -eps) .and. (t <= 1.0+eps) &
        .and. (p >= -eps) .and. (p <= 1.0+eps)) qinter = PETSC_TRUE

    !  Find intersection (two methods)

    px  = alx + t*apx
    py  = aly + t*apy
    pxc = blx + p*bpx
    pyc = bly + p*bpy

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

  GetTriangleArea = 0.5 * ((xa-xc)*(yb-ya)-(xa-xb)*(yc-ya))

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

function QArea(f, idir)

  !
  ! Routine to get component id of the vector area of the quadrilateral f
  !
  ! Author: Dave Ponting
  ! Date: 02/05/18

  implicit none

  PetscReal QArea
  PetscReal, intent(in) :: f(0:1, 0:1, 3)
  PetscInt , intent(in) :: idir
  PetscReal :: x00, x01, x10, x11
  PetscReal :: y00, y01, y10, y11
  PetscReal :: at1, at2

  PetscInt :: jdir, kdir

  !  Find other directions

  call GetOtherDirections(idir, jdir, kdir)

  !  Get projection into the 2-plane orthogonal to e(idir)

  x00 = f(0, 0, jdir)
  y00 = f(0, 0, kdir)
  x10 = f(1, 0, jdir)
  y10 = f(1, 0, kdir)
  x01 = f(0, 1, jdir)
  y01 = f(0, 1, kdir)
  x11 = f(1, 1, jdir)
  y11 = f(1, 1, kdir)

  !  Obtain the area as two triangles

  at1 = GetTriangleArea(x00, y00, x01, y01, x11, y11)
  at2 = GetTriangleArea(x00, y00, x11, y11, x10, y10)

  !  Sum to get result

  QArea = at1+at2

end function QArea

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
    alx = f(0, 0, jdir)
    aly = f(0, 0 , kdir)
    aux = f(1, 0, jdir)
    auy = f(1, 0 , kdir)
  endif

  if (iseg == 2) then
    ! Segment 2 is (1,0) to (1,1)
    alx = f(1, 0, jdir)
    aly = f(1, 0, kdir)
    aux = f(1, 1, jdir)
    auy = f(1, 1, kdir)
  endif

  if (iseg == 3) then
    ! Segment 3 is (1,1) to (0,1)
    alx = f(1, 1, jdir)
    aly = f(1, 1, kdir)
    aux = f(0, 1, jdir)
    auy = f(0, 1, kdir)
  endif

  if (iseg == 4) then
    ! Segment 1 is (0,1) to (0,0)
    alx = f(0, 1, jdir)
    aly = f(0, 1, kdir)
    aux = f(0, 0, jdir)
    auy = f(0, 0, kdir)
  endif

end subroutine GetQuadSegment

! *************************************************************************** !

subroutine GetCentre(x000, x100, x010, x110, &
                      x001, x101, x011, x111, centre)
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
    centre(id) = 0.125*(x000(id) + x100(id) + x010(id) + x110(id) + &
                         x001(id) + x101(id) + x011(id) + x111(id))
  enddo

end subroutine GetCentre

! *************************************************************************** !

subroutine SetMapAxes(mapaxes)

  implicit none

  PetscReal, intent(in) :: mapaxes(6)

  g_mapaxes = mapaxes
  g_mapaxes_set = PETSC_TRUE

end subroutine SetMapAxes

! *************************************************************************** !

function GetMapAxes(mapaxes)

  PetscBool :: GetMapAxes
  PetscReal, intent(out) :: mapaxes(6)

  mapaxes    = g_mapaxes
  GetMapAxes = g_mapaxes_set

end function GetMapAxes

! *************************************************************************** !

subroutine SetMapAxesUnits(mapaxesunits)

  implicit none

  character(len=*), intent(in) :: mapaxesunits

  g_mapaxesunits = mapaxesunits

end subroutine SetMapAxesUnits

! *************************************************************************** !

subroutine SetGDORIENT(zbuf,qerr,zmess)
  !
  ! Check and store GDORIENT keyword info
  !
  ! Author: Dave Ponting
  ! Date: 08/12/21

  use String_module, only : StringCompareIgnoreCase

  implicit none

  character(len = MAXSTRINGLENGTH), intent(in)  :: zbuf(:)
  PetscBool,                        intent(out) :: qerr
  character(len = MAXSTRINGLENGTH), intent(out) :: zmess

  PetscInt  :: nzbuf
  PetscBool :: is_inc_i, is_inc_j, is_inc_k, is_down, is_left_hand

! Initialise

  qerr  = PETSC_FALSE
  zmess = 'OK'
  nzbuf = size(zbuf)

!  Check that buffer is large enough

  if(nzbuf>=5) then

    !  Extract and check values

    is_inc_i = StringCompareIgnoreCase(zbuf(1), 'inc')
    is_inc_j = StringCompareIgnoreCase(zbuf(2), 'inc')
    is_inc_k = StringCompareIgnoreCase(zbuf(3), 'inc')

    if (.not. (is_inc_i .and. is_inc_j .and. is_inc_k)) then
      zmess = 'All three GDORIENT orders must be INC'
      qerr  = PETSC_TRUE
    endif

    is_down = StringCompareIgnoreCase(zbuf(4), 'down')
    if (.not.is_down) then
      zmess = 'GDORIENT depth convention must be DOWN'
      qerr  = PETSC_TRUE
    endif

    is_left_hand = StringCompareIgnoreCase(zbuf(5), 'left')

    g_SetGDORIENT    = PETSC_TRUE
    g_islefthandgrid = is_left_hand

  else

    zmess = 'GDORIENT expected with at least 5 items of data'
    qerr  = PETSC_TRUE

  endif

end subroutine SetGDORIENT

! *************************************************************************** !

subroutine GetMapAxesUnits(mapaxesunits)

  character(len=*), intent(out) :: mapaxesunits

  mapaxesunits = g_mapaxesunits

end subroutine GetMapAxesUnits

! *************************************************************************** !

function GetGDORIENT(islefthandgrid)

  PetscBool :: GetGDORIENT
  PetscBool , intent(out) :: islefthandgrid

  GetGDORIENT     = g_SetGDORIENT
  islefthandgrid  = g_islefthandgrid

end function GetGDORIENT

! *************************************************************************** !

subroutine BroadcastInt(ival, option)
  !
  ! Send a PetscInt value from the I/O rank to all the others
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Option_module

  implicit none

  PetscInt, intent(inout) :: ival
  type(option_type) :: option
  PetscMPIInt :: int_mpi(1)
  PetscInt    :: ierr

  ierr    = 0
  int_mpi(1) = ival

  call MPI_Bcast(int_mpi, ONE_INTEGER_MPI, MPI_INTEGER, option%comm%io_rank, &
                 option%mycomm, ierr)

  ival = int_mpi(1)

end subroutine BroadcastInt

! *************************************************************************** !

subroutine BroadcastIntN(ival, option)
  !
  ! Send a PetscInt array from the I/O rank to all the others
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Option_module

  implicit none

  PetscMPIInt, intent(inout) :: ival(:)
  type(option_type) :: option

  PetscInt    :: n, ierr
  PetscMPIInt :: nmpi

  ierr = 0
  n    = size(ival)
  nmpi = n

  call MPI_Bcast(ival, nmpi, MPI_INTEGER, option%comm%io_rank, &
                 option%mycomm, ierr)

end subroutine BroadcastIntN

! *************************************************************************** !

subroutine BroadcastRealN(rval, option)
  !
  ! Send a PetscReal array from the I/O rank to all the others
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Option_module
  implicit none

  PetscReal, intent(inout) :: rval(:)
  type(option_type) :: option

  PetscInt    :: n, ierr
  PetscMPIInt :: nmpi

  ierr = 0
  n    = size(rval)
  nmpi = n

  call MPI_Bcast(rval, nmpi, MPI_DOUBLE_PRECISION, &
                  option%comm%io_rank, option%mycomm, ierr)

end subroutine BroadcastRealN

! *************************************************************************** !

subroutine BroadcastCharN(zval, option)
  !
  ! Send a character array from the I/O rank to all the others
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Option_module
  implicit none

  character(len=*), intent(inout) :: zval(:)
  type(option_type) :: option

  PetscInt    :: n, nlen, ierr
  PetscMPIInt :: nmpi

  ierr = 0
  n    = size(zval)
  nlen = len (zval)

  nmpi = n*nlen

  call MPI_Bcast(zval, nmpi , MPI_CHARACTER, &
                 option%comm%io_rank, option%mycomm, ierr)

end subroutine BroadcastCharN

! *************************************************************************** !

subroutine ParallelSendToIO(rval, option)
  !
  ! Send a PetscReal value from this rank to the I/O proc
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Option_module
  implicit none

  PetscReal, intent(in) :: rval
  type(option_type) :: option
  PetscReal :: temp_real_array(1)

  PetscInt    :: ierr
  PetscMPIInt :: nmpi, irank
  PetscMPIInt, parameter :: tag_mpi = 0

  ierr = 0
  nmpi = 1
  irank = option%comm%io_rank
  temp_real_array(1) = rval

  call MPI_Send(temp_real_array, nmpi, MPI_DOUBLE_PRECISION, &
                 irank, tag_mpi, option%mycomm, ierr)

end subroutine ParallelSendToIO

! *************************************************************************** !

function ParallelRecvFromProc(iproc, option)
  !
  ! Receive a PetscReal value from rank iproc (will be to I/O proc)
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Option_module

  implicit none

  PetscReal :: ParallelRecvFromProc
  PetscInt, intent(in) :: iproc
  type(option_type) :: option
  PetscReal :: rval1(1)
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)

  PetscInt    :: ierr
  PetscMPIInt :: nmpi
  PetscReal   :: ParallelRecvFromReal

  ParallelRecvFromReal = 0.0

  ierr = 0
  nmpi = 1

  call MPI_Recv(rval1, nmpi, MPI_DOUBLE_PRECISION, iproc, &
                 MPI_ANY_TAG, option%mycomm, status_mpi, ierr)

  ParallelRecvFromProc = rval1(1)

end function ParallelRecvFromProc

! *************************************************************************** !

subroutine SetProbDims(nx,ny,nz)

  !
  ! Store the dimensions so can be obtained from other modules
  ! without dependency issues
  !
  ! Author: Dave Ponting
  ! Date: 12/08/20

  implicit none

  PetscInt, intent(in) :: nx
  PetscInt, intent(in) :: ny
  PetscInt, intent(in) :: nz

  gu_nx=nx
  gu_ny=ny
  gu_nz=nz

end subroutine SetProbDims

! *************************************************************************** !

subroutine GetProbDims(nx,ny,nz)

  !
  ! Provide the dimensions so can be obtained from other modules
  ! without dependency issues
  !
  ! Author: Dave Ponting
  ! Date: 12/08/20

  implicit none

  PetscInt, intent(out) :: nx
  PetscInt, intent(out) :: ny
  PetscInt, intent(out) :: nz

  nx = gu_nx
  ny = gu_ny
  nz = gu_nz

end subroutine GetProbDims

! *************************************************************************** !

subroutine SetDefaultLocalsPerGlobal(nl,ng,nlpg,qisz)

  ! Set up default local per global split
  !
  ! Author: Dave Ponting
  ! Date: 10/12/21
  !

  use Characteristic_Curves_module

  implicit none

  PetscInt ,intent(in)  :: nl, ng
  PetscInt ,intent(out) :: nlpg(:)
  PetscBool, intent(in) :: qisz
  PetscInt :: irem, nrem, j

  ! Basic cells on an even split and any remainder

  nlpg=    nl/ng
  nrem=mod(nl,ng)

  !  Distribute any remaining cells

  do irem=1,nrem
    j = irem
    if(qisz) j = nrem - irem + 1
    nlpg(j) = nlpg(j) + 1
  enddo

end subroutine SetDefaultLocalsPerGlobal

! **************************************************************************** !

subroutine MakeLocalBasePointers(nlpg,blpg,ndl,ndg,qerr,ndlis)

  implicit none

  PetscInt,intent (in) :: nlpg(:)
  PetscInt,intent (out) :: blpg(:)
  PetscInt,intent (in)  :: ndl,ndg
  PetscBool,intent(out) :: qerr
  PetscInt,intent (out) :: ndlis

  PetscInt :: idg

  ! Initialise

  qerr  = PETSC_FALSE

  !  Find ndl implied by data

  ndlis = 0
  do idg=1,ndg
    blpg(idg)=ndlis
    ndlis = ndlis + nlpg(idg)
  enddo

  if(ndlis /= ndl) qerr = PETSC_TRUE

end subroutine MakeLocalBasePointers

! **************************************************************************** !

subroutine SetNLGR(nlgr)

  implicit none

  PetscInt, intent(in) :: nlgr

  gu_nlgr = nlgr

end subroutine SetNLGR

! **************************************************************************** !

subroutine SetLGRName(ilgr,zlgr)

  implicit none

  PetscInt,        intent(in) :: ilgr
  character(len=*),intent(in) :: zlgr

  if(ilgr>0 .and. ilgr<=gu_nlgr) then
    gu_zlgr(ilgr) = zlgr
  endif

end subroutine SetLGRName

! **************************************************************************** !

subroutine SetLGRDimensions(ilgr,nx,ny,nz)

  implicit none

  PetscInt,intent(in) :: ilgr
  PetscInt,intent(in) :: nx,ny,nz

  if(ilgr>0 .and. ilgr<=gu_nlgr) then
    gu_lgrd(1,ilgr) = nx
    gu_lgrd(2,ilgr) = ny
    gu_lgrd(3,ilgr) = nz
  endif

end subroutine SetLGRDimensions

! **************************************************************************** !

function GetLGRLocation(zlgr)

  use String_module,only : StringCompareIgnoreCase

  implicit none

  PetscInt :: GetLGRLocation
  character(len=*),intent(in) :: zlgr
  PetscInt :: ilgr,jlgr

  ilgr = g_ifld

  do jlgr=1,gu_nlgr
    if (StringCompareIgnoreCase(zlgr,gu_zlgr(jlgr))) then
      ilgr = jlgr
      exit
    endif
  enddo

  GetLGRLocation = ilgr

end function GetLGRLocation

! **************************************************************************** !

function GetIsLGR()

  implicit none

  PetscBool :: GetIsLGR

  GetisLGR = PETSC_FALSE
  if(gu_nlgr > 1) GetIsLGR=PETSC_TRUE

end function GetIsLGR

! **************************************************************************** !

function GetNLGR()

  implicit none

  PetscInt :: GetNLGR

  GetNLGR=gu_nlgr

end function GetNLGR

! **************************************************************************** !

function GetLGRName8(ilgr)

  implicit none

  character(len=8) :: GetLGRName8
  PetscInt, intent(in) :: ilgr
  character(len=MAXWORDLENGTH) :: zlgr

  GetLGRName8 = ':+:+:+:+'
  if(ilgr>0 .and. ilgr<=gu_nlgr) then
    zlgr = gu_zlgr(ilgr)
    call CopyString(GetLGRName8,zlgr)
  endif

end function GetLGRName8

! **************************************************************************** !

subroutine GetLGRDimensions(ilgr,nx,ny,nz)

  implicit none

  PetscInt,intent(in) :: ilgr
  PetscInt,intent(out) :: nx,ny,nz

  nx = 1
  ny = 1
  nz = 1

  if(ilgr>0 .and. ilgr<=gu_nlgr) then
    nx = gu_lgrd(1,ilgr)
    ny = gu_lgrd(2,ilgr)
    nz = gu_lgrd(3,ilgr)
  endif

end subroutine GetLGRDimensions

! **************************************************************************** !

function GetDRPSubindex(zkey)
  !
  ! Get sub-array for a given directional rel. perm. keyword
  !
  ! Author: Dave Ponting
  ! Date: 05/16/18

  use String_module,only : StringCompareIgnoreCase

  implicit none

  PetscInt :: GetDRPSubindex

  character(len=*),intent(in) :: zkey

  GetDRPSubindex = 0

  if(StringCompareIgnoreCase(zkey,'KRNUMX')) GetDRPSubindex = idrpa_krx
  if(StringCompareIgnoreCase(zkey,'KRNUMY')) GetDRPSubindex = idrpa_kry
  if(StringCompareIgnoreCase(zkey,'KRNUMZ')) GetDRPSubindex = idrpa_krz

  if(StringCompareIgnoreCase(zkey,'IMBNUMX')) GetDRPSubindex = idrpa_imx
  if(StringCompareIgnoreCase(zkey,'IMBNUMY')) GetDRPSubindex = idrpa_imy
  if(StringCompareIgnoreCase(zkey,'IMBNUMZ')) GetDRPSubindex = idrpa_imz

end function GetDRPSubindex

! **************************************************************************** !

function GetDRPSubindexName(idrpa)
  !
  ! Get sub-array for a given directional rel. perm. keyword
  !
  ! Author: Dave Ponting
  ! Date: 05/16/18

  use String_module,only : StringCompareIgnoreCase

  implicit none

  character(len=8) :: GetDRPSubindexName
  PetscInt, intent(in) :: idrpa

  GetDRPSubindexName = ' '

  if(idrpa == idrpa_krx) GetDRPSubindexName = 'KRNUMX '
  if(idrpa == idrpa_kry) GetDRPSubindexName = 'KRNUMY '
  if(idrpa == idrpa_krz) GetDRPSubindexName = 'KRNUMZ '

  if(idrpa == idrpa_imx) GetDRPSubindexName = 'IMBNUMX'
  if(idrpa == idrpa_imy) GetDRPSubindexName = 'IMBNUMY'
  if(idrpa == idrpa_imz) GetDRPSubindexName = 'IMBNUMZ'

end function GetDRPSubindexName

! **************************************************************************** !

function GetNRepReg()

  implicit none

  PetscInt :: GetNRepReg

  GetNRepReg = g_total_fip + 1

end function GetNRepReg

! **************************************************************************** !

subroutine SplitRptreg(irepreg,ifipa,ifipv)

  implicit none

  PetscInt, intent(in) :: irepreg
  PetscInt, intent(out) :: ifipa,ifipv
  PetscInt :: ifipreg

  PetscInt :: jfipa, ib, nv

  ifipa = 0 ! For whole field - deemed fip region 0
  ifipv = 1

  if(irepreg > 1) then
    ifipreg=irepreg-1
    do jfipa=1,mfipa
      ib = g_b_fipreg(jfipa) ! Note this is a zero-based pointer
      nv = g_n_fipreg(jfipa)
      if(ifipreg>=ib+1 .and. ifipreg<=(ib+nv)) then
        ifipa = jfipa
        ifipv = ifipreg - ib
      endif
    enddo
  endif

end subroutine SplitRptreg

! **************************************************************************** !

function GetBaseFip(ifipa)

  implicit none

  PetscInt,intent(in) :: ifipa
  PetscInt :: ib,GetBaseFip

  ib = 0
  if(ifipa>0 .and. ifipa<mfipa) then
    ib = g_b_fipreg(ifipa)
  endif

  GetBaseFip = ib

end function GetBaseFip

! **************************************************************************** !

function GetNinFip(ifipa)

  implicit none

  PetscInt,intent(in) :: ifipa
  PetscInt :: n,GetNinFip

  n = 0
  if(ifipa>0 .and. ifipa<mfipa) then
    n = g_n_fipreg(ifipa)
  endif

  GetNinFip = n

end function GetNinFip

! **************************************************************************** !

subroutine GetFipaName(ifipa,fipaname)

  implicit none

  PetscInt, intent(in) :: ifipa
  character(len=MAXWORDLENGTH),intent(out) :: fipaname

  fipaname = 'Field'
  if(ifipa>0 .and. ifipa<=mfipa) then
    fipaname = g_fipname(ifipa)
  endif

end subroutine GetFipaName

! **************************************************************************** !

subroutine AllocArgs(argi,argr,argz,nargi,nargr,nargz)

  implicit none

  PetscInt                      ,intent(inout) ,allocatable :: argi(:)
  PetscReal                     ,intent(inout) ,allocatable :: argr(:)
  character(len=MAXSTRINGLENGTH),intent(inout) ,allocatable :: argz(:)
  PetscInt, intent(in) :: nargi, nargr, nargz

  if(nargi>0) then
    allocate(argi(nargi))
  else
    if(allocated(argi)) deallocate(argi)
  endif

  if(nargr>0) then
    allocate(argr(nargr))
  else
    if(allocated(argr)) deallocate(argr)
  endif

  if(nargz>0) then
    allocate(argz(nargz))
  else
    if(allocated(argz)) deallocate(argz)
  endif

end subroutine AllocArgs

! **************************************************************************** !

subroutine DeallocArgs(argi,argr,argz)

  implicit none

  PetscInt                      ,intent(inout) ,allocatable :: argi(:)
  PetscReal                     ,intent(inout) ,allocatable :: argr(:)
  character(len=MAXSTRINGLENGTH),intent(inout) ,allocatable :: argz(:)

  if(allocated(argi)) deallocate(argi)
  if(allocated(argr)) deallocate(argr)
  if(allocated(argz)) deallocate(argz)

end subroutine DeallocArgs

! **************************************************************************** !

subroutine StoreCGL(cg,cl,io,ncgl,ibclg,ndclg,ntclg)

  !
  ! Store mappings used to set up local grids
  ! if missing from base run of restart
  !
  ! Author: Dave Ponting
  ! Date: 08/25/22
  !
  ! cg    - global compressed natural location of global cell
  ! cl    - local  compressed natural location of local  cell
  ! io    - offset into buffer for each local  cell
  ! ncgl  - number of values used in cg,cl,io arrays
  !
  ! ibclg - base pointers  into work array in compressed natural order by LGR
  ! ndclg - number of values in work array in compressed natural order by LGR
  ! ntclg - size of required work array

  implicit none

  PetscInt,intent(in) :: cg(:),cl(:),io(:),ncgl
  PetscInt,intent(in) :: ibclg(g_mlgr),ndclg(g_mlgr),ntclg

  PetscInt :: icgl

!---------------------------------------------------------------------------------
!  First, store the g->l mappings.
!---------------------------------------------------------------------------------

  g_ncgl  = ncgl

  allocate(g_cg(ncgl))
  allocate(g_cl(ncgl))
  allocate(g_io(ncgl))

  do icgl=1,ncgl
    g_cg(icgl) = cg(icgl)
    g_cl(icgl) = cl(icgl)
    g_io(icgl) = io(icgl)
  enddo

!--------------------------------------------------------------------------------
!  Now store compressed Eclipse natural order LG buffer base pointers, LG block
!  sizes and total buffer size. These hold the LG solutions when main grid read.
!--------------------------------------------------------------------------------

  g_ibclg = ibclg
  g_ndclg = ndclg
  g_ntclg = ntclg

end subroutine StoreCGL

! **************************************************************************** !

function GetNCGL(ibclg,ndclg,ntclg)

  !
  ! Supply mappings used to set up local grids
  ! if missing from base run of restart
  !
  ! Author: Dave Ponting
  ! Date: 08/25/22

  implicit none

  PetscInt,intent(out) :: ibclg(g_mlgr)
  PetscInt,intent(out) :: ndclg(g_mlgr)
  PetscInt,intent(out) :: ntclg

  PetscInt :: GetNCGL

  GetNCGL = g_ncgl

  ibclg   = g_ibclg
  ndclg   = g_ndclg
  ntclg   = g_ntclg

end function GetNCGL

! **************************************************************************** !

subroutine FindCGL(icgl,icg,icl,iob)

  implicit none

  PetscInt, intent(in) :: icgl
  PetscInt, intent(out) :: icg, icl, iob

  icg = g_cg(icgl)
  icl = g_cl(icgl)
  iob = g_io(icgl)

end subroutine FindCGL

! **************************************************************************** !

subroutine ReleaseCGL

  implicit none

  if(allocated(g_cg)) deallocate(g_cg)
  if(allocated(g_cl)) deallocate(g_cl)
  if(allocated(g_io)) deallocate(g_io)

end subroutine ReleaseCGL

! ************************************************************************** !

subroutine CopyString(a,b)

  !
  ! Copy string and truncate result if required
  ! Avoids compiler warnings with default Fortran version
  !
  ! Author: Dave Ponting
  ! Date  : 10/05/21

  character(len=*),intent(out) :: a
  character(len=*),intent(in ) :: b

  PetscInt :: la, lb, ia, ib
  character(len=1) :: v

  ! Initialise output value

  a    = ''

  la = len(a)
  lb = len_trim(b)

  ! Initialise output count

  ia = 0

  ! Add in characters of b

  do ib = 1, lb
    v = b(ib:ib)
    if ( ia<la ) then
      ia = ia + 1
      a(ia:ia)=v
    else
      exit
    endif
  enddo

end subroutine CopyString

! ************************************************************************** !

subroutine PrintTidy8(va, s)
  !
  ! Write out real number in neat 8-character format
  !
  ! Author: Dave Ponting
  ! Date  : 03/08/19

  implicit none
  PetscReal, intent(in) :: va
  character(len=8), intent(out) :: s
  character(len=8) :: zn
  PetscBool :: isNeg
  PetscInt  :: exponent, ilog10
  PetscReal :: v, vl, vu

  v  = va
  vu = 9999.99D0
  vl = 1.0D-3
  isNeg    = PETSC_FALSE
  exponent = 0

  s = '     0.0' ! Default

  if (v<0.0) then
    isNeg = PETSC_TRUE;v=-v
  endif

  if (v>1.0D-10) then

    iLog10 = int(log10(1.0000001D0*v))
    if ((v >= vu) .or. (v < vl)) then
      exponent = iLog10
      v = v/10.0**iLog10
      if (isNeg) then
        if( v>9.995D0 ) then
          write(zn, '(F4.1)') v ! Will round up to 10.0
        else
          write(zn, '(F4.2)') v ! Write value as X.XX
        endif
      else
        if( v>9.9995D0 ) then
          write(zn, '(F5.2)') v ! Will round to 10.00
        else
          write(zn, '(F5.3)') v ! Write value as X.XXX
        endif
      endif
      s = trim(adjustl(zn))//'E' ! Set to 4 or 5 chars+exp => 5 or 6 chars
      write(zn, '(I2)') exponent
      s = trim(adjustl(s))//trim(adjustl(zn)) ! Add two chars exponent
    else
      if (isNeg) then
        if( v>9999.995D0 ) then
          write(s, '(F7.1)') v ! Will round to 10000.0
        else
          write(s, '(F7.2)') v ! Allow space for -9999.99
        endif
      else
        if( v>9999.9995D0 ) then
          write(s, '(F8.2)') v ! Will round to 10000.00
        else
          write(s, '(F8.3)') v ! Use all the chars up to 9999.999
        endif
      endif
    endif

    if (isNeg) then
     s = '-' // trim(adjustl(s))
    endif

    s = trim(adjustr(s))

  else
    s = '     0.0'
  endif

end subroutine PrintTidy8

! ************************************************************************** !

subroutine PrintInt8(ia, s)
  !
  ! Write out integer number in neat 8-character format
  !
  ! Author: Dave Ponting
  ! Date  : 03/08/19

  implicit none
  PetscInt, intent(in) :: ia
  character(len=8), intent(out) :: s
  character(len=8) :: zn

  zn = '       0' ! Default
  write(zn,'(I8)') ia
  s = adjustl(zn)

end subroutine PrintInt8

end module Grid_Eclipse_Util_module
