module Grid_Grdecl_Util_module

! A set of small utilities used by the grdecl and Output_Eclipse_module modules

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  public:: GetCorners
  public:: GetMDtoM2Conv
  public:: GetM2toMDConv

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
  PetscReal:: x000(3), x100(3), x010(3), x110(3), &
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

subroutine fillGeoCorner(x0, x1, d0, d1, xl, yl, zl, xu ,yu ,zu)
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

  implicit none

  PetscInt :: GetLocationInDepthBlock

  PetscInt, intent(in) :: icx, icy, icz, nx, ny
  PetscInt :: ncx, ncy, ncxy

  ncx  = 2*nx
  ncy  = 2*ny
  ncxy = ncx*ncy

  GetLocationInDepthBlock = ncxy*(icz-1)+ncx*(icy-1)+icx

end function GetLocationInDepthBlock

! ************************************************************************** !

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

! ************************************************************************** !

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

end module Grid_Grdecl_Util_module
