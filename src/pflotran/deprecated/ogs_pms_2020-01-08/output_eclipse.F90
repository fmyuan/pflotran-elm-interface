module Output_Eclipse_module

!  Module of routines which write out industry-standard Eclipse files
!  Eclipse is a trademark of the Schlumberger Corporation

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Grid_Grdecl_Util_module
  use, intrinsic :: iso_fortran_env, only: int32, real32, real64

  implicit none

  private

  ! Block types in Eclipse files: integer, single, double, bool and char

  PetscInt, parameter :: e_typeI = 1
  PetscInt, parameter :: e_typeS = 2
  PetscInt, parameter :: e_typeD = 3
  PetscInt, parameter :: e_typeB = 4
  PetscInt, parameter :: e_typeC = 5

  ! Size of headers of various types

  PetscInt, parameter :: nInteHead = 500
  PetscInt, parameter :: nLogihead = 200
  PetscInt, parameter :: nDoubhead = 200

  ! Status flags and counts

  PetscBool :: e_opened = PETSC_FALSE

  PetscInt :: e_istep_summ = 0
  PetscInt :: e_sequn_rest = 0

  PetscBool :: e_firstSummaryWrite = PETSC_TRUE

  PetscInt :: e_fileunit = 0

  ! Flags and problem dimensions

  PetscBool :: e_formatted = PETSC_FALSE

  PetscInt :: e_nx   = 1
  PetscInt :: e_ny   = 1
  PetscInt :: e_nz   = 1
  PetscInt :: e_na   = 1
  PetscInt :: e_nxy  = 1
  PetscInt :: e_nxyz = 1

  ! Number of summary items to be written

  PetscInt :: e_nwell  = 0
  PetscInt :: e_nwelmx = 1
  PetscInt :: e_ncwmax = 1
  PetscInt :: e_ngroup = 0
  PetscInt :: e_ngrpmx = 1
  PetscInt :: e_nwgmax = 1

  ! Dimensions of group/well/completion arrays

  PetscInt :: e_nigrp = 50
  PetscInt :: e_nsgrp = 10
  PetscInt :: e_nxgrp = 10
  PetscInt :: e_nzgrp = 3

  PetscInt :: e_niwel = 120
  PetscInt :: e_nswel = 10
  PetscInt :: e_nxwel = 10
  PetscInt :: e_nzwel = 3

  PetscInt :: e_nicon = 20 
  PetscInt :: e_nscon = 10
  PetscInt :: e_nxcon = 10

  PetscInt :: e_nlmax = 1
  PetscInt :: e_mlmax = 1

  ! Mapping arrays (if allocated, released from ReleaseEwriterBuffers)

  PetscInt, allocatable :: e_atoc(:)
  PetscBool :: e_atoc_allocated = PETSC_FALSE
  PetscInt, allocatable :: e_ltocp(:,:)
  PetscInt, allocatable :: e_nlmaxp(:)
  PetscBool :: e_ltoap_allocated = PETSC_FALSE

  ! The public interface to this module

  public :: selectFormattedFiles
  public :: WriteEclipseFilesGrid
  public :: WriteEclipseFilesInit
  public :: WriteEclipseFilesSpec
  public :: WriteEclipseFilesSumm
  public :: WriteEclipseFilesRest
  public :: ReleaseEwriterBuffers
  public :: SetupRestMaps, GetMlmax

  private ::  SetProblemSize, &
              WriteSpecFile, &
              WriteSummFile, &
              WriteGridFile, &
              WriteInitFile, &
              WriteRestFile

contains

! *************************************************************************** !

subroutine SelectFormattedFiles()
  !
  ! Eclipse files can be formatted or unformatted.
  ! The default is unformatted, but this routine selects formatted
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  e_formatted = PETSC_TRUE

end subroutine SelectFormattedFiles

! *************************************************************************** !

subroutine WriteEclipseFilesGrid(efilename, nx, ny, nz, &
                                 coord, zcorn, gtoa, nw, mcpw)
  !
  ! The first Eclipse file to be written is the grid file
  ! This call opens all five types of file and sets e_opened
  ! The problem dimensions (nx etc) are stored, and the grid file is written
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  character(len = *), intent(in) :: efilename
  PetscInt, intent(in) :: nx, ny, nz, nw, mcpw
  PetscReal, intent(in) :: coord(:)
  PetscReal, intent(in) :: zcorn(:)
  PetscInt , intent(in) :: gtoa (:)

  character(len = MAXSTRINGLENGTH) :: specfn, summfn, gridfn, initfn, restfn

  if (e_formatted) then

    specfn = trim(efilename) // '.fsmspec'
    summfn = trim(efilename) // '.funsmry'
    gridfn = trim(efilename) // '.fgrid'
    initfn = trim(efilename) // '.finit'
    restfn = trim(efilename) // '.funrst'

    open(unit = UNIT_SPEC_WRITE, status = 'replace', form = 'formatted', &
         file = specfn)
    open(unit = UNIT_SUMM_WRITE, status = 'replace', form = 'formatted', &
         file = summfn)
    open(unit = UNIT_GRID_WRITE, status = 'replace', form = 'formatted', &
         file = gridfn)
    open(unit = UNIT_INIT_WRITE, status = 'replace', form = 'formatted', &
         file = initfn)
    open(unit = UNIT_REST_WRITE, status = 'replace', form = 'formatted', &
         file = restfn)

  else

    specfn = trim(efilename) // '.smspec'
    summfn = trim(efilename) // '.unsmry'
    gridfn = trim(efilename) // '.grid'
    initfn = trim(efilename) // '.init'
    restfn = trim(efilename) // '.unrst'

    open(unit = UNIT_SPEC_WRITE, status = 'replace', form = 'unformatted', &
         file = specfn, convert = 'big_endian')
    open(unit = UNIT_SUMM_WRITE, status = 'replace', form = 'unformatted', &
         file = summfn, convert = 'big_endian')
    open(unit = UNIT_GRID_WRITE, status = 'replace', form = 'unformatted', &
         file = gridfn, convert = 'big_endian')
    open(unit = UNIT_INIT_WRITE, status = 'replace', form = 'unformatted', &
         file = initfn, convert = 'big_endian')
    open(unit = UNIT_REST_WRITE, status = 'replace', form = 'unformatted', &
         file = restfn, convert = 'big_endian')

  endif

  ! Set flag indicating all the files have been opened

  e_opened = PETSC_TRUE

  ! Store problem size

  call SetProblemSize(nx, ny, nz, nw, mcpw)

  ! Write grid file

  call WriteGridFile(coord, zcorn, gtoa)

  ! Grid file output complete, so can close

  close(unit = UNIT_GRID_WRITE)

end subroutine WriteEclipseFilesGrid

! *************************************************************************** !

subroutine WriteEclipseFilesInit(kx, ky, kz, mx, my, mz, &
                                 depth, poro, ntg, bvol, gtoa, atoc)
  !
  ! Writes init file: this contains static values like pore volumes and perms
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  PetscReal, intent(in) :: kx(:), ky(:), kz(:), mx(:), my(:), mz(:), &
                           depth(:), poro(:), ntg(:), bvol(:)
  PetscInt , intent(in) :: gtoa(:)
  PetscInt , intent(in) :: atoc(:)

  if (e_opened) then

  ! Write init file

    call WriteInitFile(kx, ky, kz, &
                       mx, my, mz, depth, poro, ntg, bvol, gtoa, atoc)

  ! Init file complete, so can close

    close(unit = UNIT_INIT_WRITE)

  endif

end subroutine WriteEclipseFilesInit

! *************************************************************************** !

subroutine WriteEclipseFilesSpec(zm, zn, zu, ni, is_restart, restart_filename)
  !
  ! Writes spec file: this is an index to the values in the summary files
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  character(len = 8), intent(in) :: zm(:), zn(:), zu(:)
  PetscInt, intent(in) :: ni
  PetscBool, intent(in) :: is_restart
  character(len=8), intent(in) :: restart_filename(9)

  if (e_opened) then

  ! Write spec file

    call WriteSpecFile(zm, zn, zu, ni, is_restart, restart_filename)

  ! Spec file complete, so can close

    close(unit = UNIT_SPEC_WRITE)

  endif

end subroutine WriteEclipseFilesSpec

! *************************************************************************** !

subroutine WriteEclipseFilesSumm(vd, nd)
  !
  ! Writes summary file for this step
  ! Contains well, group and field information, mainly for line graph plotting
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  PetscReal, intent(in) :: vd(:)
  PetscInt , intent(in) :: nd

  if (e_opened) then

  ! Write summary file

    call WriteSummFile(vd, nd)

  ! Output is for this step, so flush so that files can be read at run-time

    flush(UNIT_SUMM_WRITE)

  endif

end subroutine WriteEclipseFilesSumm

! *************************************************************************** !

subroutine WriteEclipseFilesRest( vsoll, nsol, zsol, time, is_ioproc, &
                                  wname, wtype, wncmpl, &
                                  ixcmpl, iycmpl, izcmpl, idcmpl, option)
  !
  ! Writes restart file for this step
  ! Contains arrays like saturations, plus well and completion status
  ! Not really restart information as such: is to drive visualisation
  ! Note that (unlike other WriteEclipse routines) this is called by all ranks
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  use Option_module

  implicit none

  PetscReal, pointer :: vsoll(:,:)
  PetscInt           :: nsol
  character(len = 8), pointer :: zsol(:)
  PetscReal, intent(in) :: time
  PetscBool, intent(in) :: is_ioproc
  character(len = 8), intent(in) :: wname(:)
  PetscInt, intent(in) :: wtype(:), wncmpl(:), &
                         ixcmpl(:), iycmpl(:), izcmpl(:), idcmpl(:)
  type(option_type), intent(in), pointer :: option

  ! Write the file

  call WriteRestFile( vsoll, nsol, zsol, time, is_ioproc, &
                      wname, wtype, wncmpl, &
                      ixcmpl, iycmpl, izcmpl, idcmpl, option )

  ! Flush on the I/O processor

  if (is_ioproc) flush(UNIT_REST_WRITE)

end subroutine WriteEclipseFilesRest

! *************************************************************************** !

subroutine SetProblemSize(nx, ny, nz, nw, ncpw)
  !
  ! Store problem size variables:
  ! grid dimensions, number of wells, number of completions per well
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  PetscInt, intent(in) :: nx, ny, nz, nw, ncpw

  e_nx = nx
  e_ny = ny
  e_nz = nz

  e_nxy  = e_nx*e_ny
  e_nxyz = e_nx*e_ny*e_nz

  e_nwell  = nw
  e_nwelmx = max(e_nwelmx, e_nwell)
  e_ncwmax = max(3*nz, ncpw)

  e_ngroup = 1
  e_ngrpmx = max(e_ngrpmx, e_ngroup)
  e_nwgmax = nw

  e_nigrp  = e_nwgmax+50

end subroutine SetProblemSize

! *************************************************************************** !

subroutine WriteSpecFile(vmnem, vwgname, vunits, ni, &
                         is_restart, restart_filename)
  !
  ! Write the spec file, containing mnemonics, well/group names and units
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  character(len = 8), intent(in) :: vmnem(:), vwgname(:), vunits(:)
  PetscInt, intent(in) :: ni
  PetscBool, intent(in) :: is_restart
  character(len=8), intent(in) :: restart_filename(9)

  PetscInt, parameter :: nDimens    = 6
  PetscInt, parameter :: nStartData = 6

  PetscInt :: vdimens  (nDimens   )
  PetscInt :: vstartdat(nStartData)

  PetscInt, allocatable :: vnums(:)

  ! Store the file pointer

  e_fileunit = UNIT_SPEC_WRITE

  ! Allocations

  allocate(vnums(ni))

  ! Set up DIMENS for _ni summary items

  vdimens = 0

  vdimens(1) = ni
  vdimens(2) = e_nx
  vdimens(3) = e_ny
  vdimens(4) = e_nz

  ! Set up NUMS

   vnums = 1

  ! Set up STARTDAT (not really specified in Pflotran, so set to 1/1/2000

  vstartdat    = 1
  vstartdat(3) = 2000

  ! Write out values

  if (is_restart) then
    call WriteBlockC(restart_filename, 'RESTART', 9)
  endif
  call WriteBlockI(vdimens  , 'DIMENS'  , nDimens   )
  call WriteBlockC(vmnem    , 'KEYWORDS', ni)
  call WriteBlockC(vwgname  , 'WGNAMES' , ni)
  call WriteBlockI(vnums    , 'NUMS'    , ni)
  call WriteBlockC(vunits   , 'UNITS'   , ni)
  call WriteBlockI(vstartdat, 'STARTDAT', nStartData)

  ! Deallocations

  deallocate(vnums)

end subroutine WriteSpecFile

! *************************************************************************** !

subroutine WriteSummFile(vd, nd)
  !
  ! Write the summ file, containing values to match the mnemonics in spec file
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  PetscReal, intent(in) :: vd(:)
  PetscInt , intent(in) :: nd
  PetscInt              :: vi(1)

  e_fileunit = UNIT_SUMM_WRITE

  e_istep_summ = e_istep_summ + 1

  ! Set up integer and real buffers for the summary data at each step

  vi = 0

  ! Write out SEQHDR first step and then MINISTEP and values at each step

  if (e_firstSummaryWrite) then
    call WriteBlockI(vi, 'SEQHDR', 1)
    e_firstSummaryWrite = PETSC_FALSE
  endif

  ! Step counter (starting from 1)

  vi(1) = e_istep_summ
  call WriteBlockI(vi, 'MINISTEP', 1)

  ! Summary values for this step

  call WriteBlockS(vd, 'PARAMS', nd)

end subroutine WriteSummFile

! *************************************************************************** !

subroutine WriteGridFile(coord, zcorn, gtoa)
  !
  ! Write the grid file, containing cell vertex information
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  PetscReal, intent(in) :: coord(:)
  PetscReal, intent(in) :: zcorn(:)
  PetscInt , intent(in) :: gtoa(:)

  PetscInt :: vdimens(3)
  character(len = 8) :: vgridunit(2)

  PetscInt, parameter :: nCoords  =  7
  PetscInt, parameter :: nCorners = 24

  PetscInt  :: vcoords (nCoords )
  PetscReal :: vcorners(nCorners)

  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3)

  PetscInt  :: ix, iy, izp, ize, igp, ige, ia

  ! Store the file pointer

  e_fileunit = UNIT_GRID_WRITE

  ! Set up and write out DIMENS and GRIDUNIT

  vdimens(1) = e_nx
  vdimens(2) = e_ny
  vdimens(3) = e_nz

  call WriteBlockI(vdimens, 'DIMENS', 3)

  vgridunit(1) = 'METRIC'
  vgridunit(2) = ' '

  call WriteBlockC(vgridunit, 'GRIDUNIT', 2)

  ! Set up storage for coord and corners data

  vcoords  = 0 
  vcorners = 0.0

  ! Loop over all the cells in the grid in Eclipse order (natural, k-down)

  igp = 0
  ige = 0

  do ize = 1, e_nz
   izp = e_nz-ize+1 ! Get the Pflotran k-up index
    do iy = 1, e_ny
      do ix = 1, e_nx

        ige = ige +1

        igp = e_nxy*(izp-1) + e_nx*(iy-1) + ix
        ia  = gtoa(igp)

  ! Get the corners of this cell

        call GetCorners( ix, iy, izp, &
                         x000, x100, x010, x110, &
                         x001, x101, x011, x111, &
                         coord, zcorn, e_nx, e_ny )

  ! Fill in COORDS data for this cell (location in ijk grid)

        vcoords(1) = ix  ! x-location (in Fortran convention)
        vcoords(2) = iy  ! y-location
        vcoords(3) = ize ! z-location
        vcoords(4) = ige  ! Cell count
        if (ia >= 1) then
          vcoords(5) = 1  ! Active indicator
        else
          vcoords(5) = 0  ! Inactive indicator
        endif

  ! Fill in Eclipse CORNERS data for this cell (location in xyz space)

  ! Eclipse tops are k-high values in Pflotran (sign-flipped)
        vcorners( 1) = x001(1);vcorners( 2) = x001(2);vcorners( 3) = -x001(3)
        vcorners( 4) = x101(1);vcorners( 5) = x101(2);vcorners( 6) = -x101(3)
        vcorners( 7) = x011(1);vcorners( 8) = x011(2);vcorners( 9) = -x011(3)
        vcorners(10) = x111(1);vcorners(11) = x111(2);vcorners(12) = -x111(3)

  ! Eclipse bottoms are k-low values in Pflotran (sign-flipped)
        vcorners(13) = x000(1);vcorners(14) = x000(2);vcorners(15) = -x000(3)
        vcorners(16) = x100(1);vcorners(17) = x100(2);vcorners(18) = -x100(3)
        vcorners(19) = x010(1);vcorners(20) = x010(2);vcorners(21) = -x010(3)
        vcorners(22) = x110(1);vcorners(23) = x110(2);vcorners(24) = -x110(3)

  ! Write out values and increment cell counter

        call WriteBlockI(vcoords  , 'COORDS' , nCoords )
        call WriteBlockS(vcorners , 'CORNERS', nCorners)

      enddo
    enddo
  enddo

end subroutine WriteGridFile

! *************************************************************************** !

subroutine WriteInitFile(kx, ky, kz, &
                         mx, my, mz, depth, poro, ntg, bvol, gtoa, atoc)
  !
  ! Write the init file, containing static cell information, like perms
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  PetscReal, intent(in) :: kx(:), ky(:), kz(:), mx(:), my(:), mz(:), &
                           depth(:), poro(:), ntg(:), bvol(:)
  PetscInt, intent(in) :: gtoa(:)
  PetscInt, intent(in) :: atoc(:)

  PetscReal, parameter :: flip = -1.0
  PetscReal, parameter :: nofl =  1.0

  PetscInt  :: intehead(nInteHead)
  PetscBool :: logihead(nLogiHead)
  PetscReal :: doubhead(nDoubHead)

  PetscInt :: vi(1)

  PetscReal, allocatable :: porv(:)
  PetscReal, allocatable :: buf(:)

  PetscInt  :: ix, iy, ize, izp, ig, ia, na
  PetscReal :: cprm, cdef

  ! Store the file pointer

  e_fileunit = UNIT_INIT_WRITE

  cdef = 1.0
  cprm = GetM2toMDConv()

  ! Set up headers

  intehead = 0
  logihead = PETSC_FALSE
  doubhead = 0.0

  ! Set up properties

  e_nxyz = e_nx*e_ny*e_nz

  allocate(porv(e_nxyz))
  allocate(buf (e_nxyz))

  ! Prepare full (all cells) porv array

  e_na = 0
  do ize = 1, e_nz
    izp = e_nz-ize+1
    do iy = 1, e_ny
      do ix = 1, e_nx
        ig = e_nx*e_ny*(izp-1) + e_nx*(iy-1) + ix
        ia = gtoa(ig)
        if (ia>0) then
          porv(ig) = poro(ig)*bvol(ia)
          e_na = e_na+1
        else
          porv(ig) = 0.0
        endif
      enddo
    enddo
  enddo

  ! Set up integer header now that active count known

  call SetInteHead(intehead)

  ! Write operations

  vi(1) = 1
  call WriteBlockI(vi, 'SEQNUM', 1)

  call WriteBlockI(intehead, 'INTEHEAD', nIntehead)
  call WriteBlockB(logihead, 'LOGIHEAD', nLogihead)
  call WriteBlockD(doubhead, 'DOUBHEAD', nDoubhead)

  call WriteBlockS(porv, 'PORV', e_nxyz)

  call CmpToCNOBuf(buf, kx   , porv, nofl, cprm)
  call WriteBlockS(buf, 'PERMX', e_na)
  call CmpToCNOBuf(buf, ky   , porv, nofl, cprm)
  call WriteBlockS(buf, 'PERMY', e_na)
  call CmpToCNOBuf(buf, kz   , porv, nofl, cprm)
  call WriteBlockS(buf, 'PERMZ', e_na)

  call CmpToCNOBuf(buf, mx   , porv, nofl, cdef)
  call WriteBlockS(buf, 'MULTX', e_na)
  call CmpToCNOBuf(buf, my   , porv, nofl, cdef)
  call WriteBlockS(buf, 'MULTY', e_na)
  call CmpToCNOBuf(buf, mz   , porv, nofl, cdef)
  call WriteBlockS(buf, 'MULTZ', e_na)

  call CmpToCNOBuf(buf, poro , porv, nofl, cdef)
  call WriteBlockS(buf, 'PORO' , e_na)
  call CmpToCNOBuf(buf, depth, porv, flip, cdef)
  call WriteBlockS(buf, 'DEPTH', e_na)
  call CmpToCNOBuf(buf, ntg  , porv, nofl, cdef)
  call WriteBlockS(buf, 'NTOG' , e_na)

  ! Deallocate

  deallocate(porv)
  deallocate(buf)

  ! Allocate the active to compressed natural mapping,
  ! check active count agrees

  allocate(e_atoc(e_na))
  e_atoc_allocated = PETSC_TRUE
  na = size(atoc)
  if (na == e_na) then
    e_atoc = atoc
  else
   call ThrowEwriterException('Active cell mismatch in WriteInitFile')
  endif

end subroutine WriteInitFile

! *************************************************************************** !

subroutine WriteRestFile(vsoll, nsol, zsol, time, is_ioproc, &
                         wname, wtype, wncmpl, &
                         ixcmpl, iycmpl, izcmpl, idcmpl, option)
  !
  ! Write the rest file, containing dynamic cell information, like saturations
  ! This routine is called by all ranks, but only the I/O rank writes the file
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  use Option_module
  use String_module, only : StringCompareIgnoreCase

  implicit none

  PetscReal, pointer        :: vsoll(:,:)
  PetscInt                  :: nsol
  character(len = 8), pointer :: zsol(:)
  PetscReal :: time, tdays
  PetscBool, intent(in)       :: is_ioproc
  character(len = 8), intent(in) :: wname(:)
  PetscInt, intent(in) :: wtype(:), wncmpl(:), &
                          ixcmpl(:), iycmpl(:), izcmpl(:), idcmpl(:)

  type(option_type), intent(in), pointer :: option

  PetscInt  :: intehead(nInteHead)
  PetscBool :: logihead(nLogiHead)
  PetscReal :: doubhead(nDoubHead)

  PetscInt :: vi(1)

  PetscReal, allocatable :: varr(:)
  PetscReal, allocatable :: vbuf(:)

  PetscReal :: conv, tconv
  PetscBool :: is_pressure,is_psat

  PetscInt :: isol, hours, mins, microsecs, years, months, days, &
              ic, il, nproc, &
              iproco, nlmaxo, iproct, ioproc, lioproc, liproct, liproco

  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE), ierr, itag

  ! Set up useful scalars

  itag = 0
  ioproc  = option%io_rank
  iproct  = option%myrank
  lioproc = ioproc+1
  liproct = iproct+1
  nproc   = option%mycommsize
  tconv   = 3600.0*24.0

  ! Write headers and wells on the I/O proc only

  if (is_ioproc) then

    e_fileunit = UNIT_REST_WRITE

  ! Set up headers

    intehead = 0
    logihead = PETSC_FALSE
    doubhead = 0.0

    call SetInteHead(intehead)
    tdays = time/tconv
    doubhead(  1) = tdays
    doubhead(162) = tdays

  ! Get date and time

    call GetYMDHMMS(tdays, years, months, days, hours, mins, microsecs)

    intehead(65)  = days
    intehead(66)  = months
    intehead(67)  = years
    intehead(207) = hours
    intehead(208) = mins
    intehead(411) = microsecs

  ! Set up array to hold values

    allocate(varr(e_na))

  ! Headers

    e_sequn_rest = e_sequn_rest+1
    vi(1) = e_sequn_rest
    call WriteBlockI(vi, 'SEQNUM', 1)

    call WriteBlockI(intehead, 'INTEHEAD' , nIntehead)
    call WriteBlockB(logihead, 'LOGIHEAD' , nLogihead)
    call WriteBlockD(doubhead, 'DOUBHEAD' , nDoubhead)

  ! Write wells

    call WriteWells(wname, wtype, wncmpl, ixcmpl, iycmpl, izcmpl, idcmpl)

  endif

  ! Write arrays : these need to be collected from all the procs

  allocate(vbuf(e_mlmax))

  ! Loop over the solution arrays to be written

  do isol = 1, nsol

    is_pressure = StringCompareIgnoreCase(zsol(isol), 'PRESSURE')
    is_psat     = StringCompareIgnoreCase(zsol(isol), 'PSAT'    )

! Conv. from Pflotran Pa to Eclipse Bars
    conv = 1.0
    if (is_pressure .or. is_psat) conv = 1.0E-5

    if (is_ioproc) then
      ! Add the local bit on this proc
      do il = 1, e_nlmax
        ic = e_ltocp(il, liproct)
        varr(ic) = conv*vsoll(il, isol)
      enddo
      ! Receive the values from other procs
      do iproco = 0, nproc-1

        if (iproco /= option%io_rank) then
          liproco = iproco+1
          nlmaxo  = e_nlmaxp(liproco)
          call MPI_Recv(vbuf, nlmaxo, MPI_DOUBLE_PRECISION, iproco, &
                        MPI_ANY_TAG, option%mycomm, status_mpi, ierr)
          do il = 1, nlmaxo
            ic = e_ltocp(il, liproco)
            varr(ic) = conv*vbuf(il)
          enddo
        endif
      enddo
      ! Write out the whole thing
      call WriteBlockS(varr, zsol(isol), e_na)
    else
    ! Send values to the IO proc
     do il = 1, e_nlmax
       vbuf(il) = vsoll(il, isol)
     enddo
     call MPI_Send(vbuf, e_nlmax, MPI_DOUBLE_PRECISION, ioproc, &
                   itag, option%mycomm, ierr)
    endif

  enddo

  ! Deallocate

  if (is_ioproc) then
    deallocate(varr)
  endif
  deallocate(vbuf)

end subroutine WriteRestFile

! *************************************************************************** !

subroutine WriteWells(wname, wtype, wncmpl, ixcmpl, iycmpl, izcmpl, idcmpl)
  !
  ! Write well information onto the rest file
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  use String_module
  use Well_Type_class

  implicit none

  character(len = 8), intent(in) :: wname(:)
  PetscInt, intent(in) :: wtype(:), wncmpl(:), &
                          ixcmpl(:), iycmpl(:), izcmpl(:), idcmpl(:)

  PetscInt , allocatable :: igrp(:)
  PetscReal, allocatable :: sgrp(:)
  PetscReal, allocatable :: xgrp(:)
  character(len = 8), allocatable :: zgrp(:)

  PetscInt , allocatable :: iwel(:)
  PetscReal, allocatable :: swel(:)
  PetscReal, allocatable :: xwel(:)
  character(len = 8), allocatable :: zwel(:)

  PetscInt , allocatable :: icon(:)
  PetscReal, allocatable :: scon(:)
  PetscReal, allocatable :: xcon(:)

  ! Set up default values

  PetscInt     , parameter :: idef = 0
  PetscReal    , parameter :: sdef = 0.0
  PetscReal    , parameter :: ddef = 0.0
  character(len = 8), parameter :: zdef = ' '

  PetscInt :: ntigrp, ntsgrp, ntxgrp, ntzgrp, &
              ntiwel, ntswel, ntxwel, ntzwel, &
              nticon, ntscon, ntxcon, ig, iwg, iw, ik

  PetscInt :: iwpx, iwpy, iwpz, iewtype, icpx, icpy, icpz, &
              ibigrp, ibzgrp, ibiwel, ibzwel, ibicon, ibscon, welltype, &
              cdd, ibcmpl, ncmpl

  PetscReal :: sccf, skh

  character(len = 8) :: name8
  character(len = MAXSTRINGLENGTH) :: name

  ! Initialise local scalars

  iwpz = 1
  sccf = 1.0
  skh  = 1000.0
  cdd  = 3

  ! Set up total size of the well and completion arrays

  ntigrp = e_ngrpmx*e_nigrp
  ntsgrp = e_ngrpmx*e_nsgrp
  ntxgrp = e_ngrpmx*e_nxgrp
  ntzgrp = e_ngrpmx*e_nzgrp

  ntiwel = e_nwelmx*e_niwel
  ntswel = e_nwelmx*e_nswel
  ntxwel = e_nwelmx*e_nxwel
  ntzwel = e_nwelmx*e_nzwel

  nticon = e_nwelmx*e_ncwmax*e_nicon
  ntscon = e_nwelmx*e_ncwmax*e_nscon
  ntxcon = e_nwelmx*e_ncwmax*e_nxcon

  ! Allocate the well arrays

  allocate(igrp(ntigrp))
  allocate(sgrp(ntsgrp))
  allocate(xgrp(ntxgrp))
  allocate(zgrp(ntzgrp))

  allocate(iwel(ntiwel))
  allocate(swel(ntswel))
  allocate(xwel(ntxwel))
  allocate(zwel(ntzwel))

  allocate(icon(nticon))
  allocate(scon(ntscon))
  allocate(xcon(ntxcon))

  ! Set up the group, well and completion data structures

  igrp = idef
  sgrp = sdef
  xgrp = ddef
  zgrp = zdef

  iwel = idef
  swel = sdef
  xwel = ddef
  zwel = zdef

  icon = idef
  scon = sdef
  xcon = ddef

  ! Fill in required values

  do ig = 1, e_ngroup

  ! Find required base pointers into well structures

     ibigrp = (ig-1)*e_nigrp
     ibzgrp = (ig-1)*e_nzgrp

  ! Set the required values

     do iwg = 1, e_nwgmax
       igrp(ibigrp+iwg) = iwg
     enddo
     igrp(ibigrp+e_nwgmax+1) = e_nwgmax
     zgrp(ibzgrp+1) = 'FIELD'

  enddo

  ! Loop over the well and completions filling up the data structure

  ibcmpl = 0
  do iw = 1, e_nwell

    ! Get name and type of this well

    name = wname(iw)
    call StringToUpper(name)
    welltype = wtype(iw)
    name8 = name(1:8)

    iewtype = 4
    if (welltype ==  PROD_WELL_TYPE   ) iewtype = 1
    if (welltype ==  OIL_INJ_WELL_TYPE) iewtype = 2
    if (welltype ==  WAT_INJ_WELL_TYPE) iewtype = 3

    ncmpl = wncmpl(iw)

    ! Set up well locations

   if (ncmpl > 0) then
      iwpx = ixcmpl(ibcmpl+1)
      iwpy = iycmpl(ibcmpl+1)
    else
      iwpx = 1
      iwpy = 1
    endif

    ! Find required base pointers into well structures

    ibiwel = (iw-1)*e_niwel
    ibzwel = (iw-1)*e_nzwel

    ! Set the required values

    iwel(ibiwel+ 1) = iwpx
    iwel(ibiwel+ 2) = iwpy
    iwel(ibiwel+ 3) = iwpz
    iwel(ibiwel+ 5) = ncmpl
    iwel(ibiwel+ 6) = 1       ! Group index (1->Field)
    iwel(ibiwel+ 7) = iewtype ! Well type (water injector or producer)
    iwel(ibiwel+11) = 1       ! Well is open

    zwel(ibzwel+ 1) = name(1:8)

    ! Loop over completions

    do ik = 1, ncmpl

      ! Set up completion location

      icpx = ixcmpl(ibcmpl+ik)
      icpy = iycmpl(ibcmpl+ik)
      icpz = izcmpl(ibcmpl+ik)
      icpz = e_nz-icpz+1 ! Convert back to Eclipse k-order

      cdd  = idcmpl(ibcmpl+ik)

      ! Find base pointers into completion structures

      ibicon = ((iw-1)*e_ncwmax+(ik-1))*e_nicon
      ibscon = ((iw-1)*e_ncwmax+(ik-1))*e_nscon

      ! Set the required values

      icon(ibicon+ 1) = ik
      icon(ibicon+ 2) =  icpx
      icon(ibicon+ 3) = icpy
      icon(ibicon+ 4) = icpz
      icon(ibicon+ 6) = 1    ! Completion is open
      icon(ibicon+14) = cdd  ! Drilling direction

      scon(ibscon+ 1) = sccf
      scon(ibscon+ 4) = skh

    enddo

    ibcmpl = ibcmpl + ncmpl

  enddo

  ! Write out data

  call WriteBlockI(iwel, 'IWEL' , ntiwel)
  call WriteBlockS(swel, 'SWEL' , ntswel)
  call WriteBlockD(xwel, 'XWEL' , ntxwel)
  call WriteBlockC(zwel, 'ZWEL' , ntzwel)

  call WriteBlockI(igrp, 'IGRP' , ntigrp)
  call WriteBlockS(sgrp, 'SGRP' , ntsgrp)
  call WriteBlockD(xgrp, 'XGRP' , ntxgrp)
  call WriteBlockC(zgrp, 'ZGRP' , ntzgrp)

  call WriteBlockI(icon, 'ICON' , nticon)
  call WriteBlockS(scon, 'SCON' , ntscon)
  call WriteBlockD(xcon, 'XCON' , ntxcon)

  ! Release the well arrays

  deallocate(iwel)
  deallocate(swel)
  deallocate(xwel)
  deallocate(zwel)

  deallocate(igrp)
  deallocate(sgrp)
  deallocate(xgrp)
  deallocate(zgrp)

  deallocate(icon)
  deallocate(scon)
  deallocate(xcon)

end subroutine WriteWells

! *************************************************************************** !

subroutine WriteBlockI(a, mnem, n)
  !
  ! Write out a block of int*32 data
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscInt, intent(in) :: a(:)
  character(len = *), intent(in) :: mnem
  PetscInt, intent(in) :: n

  PetscInt, parameter :: blksize  = 6    ! Values/line (formatted)
  PetscInt, parameter :: mrecsize = 1000 ! Values/rec (unformatted)

  PetscInt :: il, iu, j, irec, nrec, ninrec
  integer(kind = int32), allocatable :: ibuf(:)
 
  ! Write out an integer header line or record

  call WriteHeader(e_typeI, mnem, n)

  ! Write out n values

  if (e_formatted) then

  ! Formatted case: write out in lines of 6 values per line, I11 format

    do il = 1, n, blksize
      iu = min(il+blksize-1, n)
      write(e_fileunit, '(6(1x,i11))') (a(j), j = il, iu)
    enddo

  else

  ! Allocate an int32 buffer

    allocate(ibuf(mrecsize))

    ! Loop over the number of records required

    nrec = GetNumberOfRecords(n, mrecsize)
    do irec = 1, nrec
      ! Get record details
      call GetRecordDetails(irec, n, mrecsize, il, iu, ninrec)
      ! Copy to buffer
      call CopyToBufferI(a, ibuf, il, iu)
      ! Write out the record
      write(e_fileunit) (ibuf(j), j = 1, ninrec)
    enddo

    ! Delete the buffer

    deallocate(ibuf)
  endif

end subroutine WriteBlockI

! *************************************************************************** !

subroutine WriteBlockS(a, mnem, n)
  !
  ! Write out a block of real*32 data
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscReal, intent(in) :: a(:)
  character(len = *), intent(in) :: mnem
  PetscInt, intent(in) :: n

  PetscInt, parameter :: blksize = 4 ! Values/line (formatted)
  PetscInt, parameter :: mrecsize = 1000 ! Values/rec (unformatted)

  PetscInt :: il, iu, j, irec, nrec, ninrec
  real(kind = real32), allocatable :: fbuf(:)
 
  ! Write out an single precision header line or record

  call WriteHeader(e_typeS, mnem, n)

  ! Write out n values

  if (e_formatted) then

    ! Formatted case: write out in lines of 6 values per line, I11 format

    do il = 1, n, blksize
      iu = min(il+blksize-1, n)
      write(e_fileunit, '(4(1x,e16.8))') (a(j), j = il, iu)
    enddo

  else

    ! Allocate a real32 buffer

    allocate(fbuf(mrecsize))

    ! Loop over the number of records required

    nrec = GetNumberOfRecords(n, mrecsize)
    do irec = 1, nrec
      ! Get record details
      call GetRecordDetails(irec, n, mrecsize, il, iu, ninrec)
      ! Copy to buffer
      call CopyToBufferS(a, fbuf, il, iu)
      ! Write out the record
      write(e_fileunit) (fbuf(j), j = 1, ninrec)
    enddo

    ! Delete the buffer

    deallocate(fbuf)
  endif

end subroutine WriteBlockS

! *************************************************************************** !

subroutine WriteBlockD(a, mnem, n)
  !
  ! Write out a block of real*64 data
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscReal, intent(in) :: a(:)
  character(len = *), intent(in) :: mnem
  PetscInt, intent(in) :: n

  PetscInt, parameter :: blksize = 3 ! Values/line (formatted)
  PetscInt, parameter :: mrecsize = 1000 ! Values/rec (unformatted)

  PetscInt :: il, iu, j, irec, nrec, ninrec
  real(kind = real64), allocatable :: dbuf(:)

  ! Write out an integer header line or record

  call WriteHeader(e_typeD, mnem, n)

  ! Write out n values

  if (e_formatted) then

    ! Formatted case: write out in lines of 6 values per line, I11 format

    do il = 1, n, blksize
      iu = min(il+blksize-1, n)
      write(e_fileunit, '(3(1x,e22.14))') (a(j), j = il,iu)
    enddo

  else

    ! Allocate an double precision buffer

    allocate(dbuf(mrecsize))

    ! Loop over the number of records required

    nrec = GetNumberOfRecords(n, mrecsize)
    do irec = 1, nrec
      ! Get record details
      call GetRecordDetails(irec, n, mrecsize, il, iu, ninrec)
      ! Copy to buffer
      call CopyToBufferD(a, dbuf, il, iu)
      ! Write out the record
      write(e_fileunit) (dbuf(j), j = 1, ninrec)
    enddo

    ! Delete the buffer

    deallocate(dbuf)

  endif

end subroutine WriteBlockD

! *************************************************************************** !

subroutine WriteBlockB(a, mnem, n)
  !
  ! Write out a block of bool*32 data
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscBool, intent(in) :: a(:)
  character(len = *), intent(in) :: mnem
  PetscInt, intent(in) :: n

  PetscInt, parameter :: blksize = 25 ! Values/line (formatted)
  PetscInt, parameter :: mrecsize = 1000 ! Values/rec (unformatted)

  PetscInt :: il, iu, j, irec, nrec, ninrec
  integer(kind = int32), allocatable :: ibuf(:)

  ! Write out an boolean header line or record

  call WriteHeader(e_typeB, mnem, n)

  ! Write out n bool values

  if (e_formatted) then

    ! Formatted case: write out in lines of 6 values per line, I11 format

    do il = 1, n, blksize
      iu = min(il+blksize-1, n)
      write(e_fileunit, '(25(1X,L2))') (a(j), j = il,iu)
    enddo

  else

    ! Allocate an int32 buffer

    allocate(ibuf(mrecsize))

! Loop over the number of records required

    nrec = GetNumberOfRecords(n, mrecsize)
    do irec = 1, nrec
      ! Get record details
      call GetRecordDetails(irec, n, mrecsize, il, iu, ninrec)
      ! Copy to buffer
      call CopyToBufferB(a, ibuf, il, iu)
      ! Write out the record
      write(e_fileunit) (ibuf(j), j = 1, ninrec)
    enddo

    ! Delete the buffer

    deallocate(ibuf)
  endif

end subroutine WriteBlockB

! *************************************************************************** !

subroutine WriteBlockC(a, mnem, n)
  !
  ! Write out a block of char*8 data
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  character(len = 8), intent(in) :: a(:)
  character(len = *), intent(in) :: mnem
  PetscInt, intent(in) :: n

  PetscInt, parameter :: blksize  = 7   ! Values/line (formatted)
  PetscInt, parameter :: mrecsize = 105 ! Values/rec  (unformatted)

  PetscInt :: il, iu, j, irec, nrec, ninrec

10 format(7(1X, "'", A8, "'"))

  ! Write out an integer header line or record

  call WriteHeader(e_typeC, mnem, n)

  ! Write out n integer values

  if (e_formatted) then

  ! Formatted case: write out in lines of 7 values per line, A8 format

    do il = 1, n, blksize
      iu = min(il+blksize-1, n)
      write(e_fileunit, 10) (a(j), j = il,iu)
    enddo

  else

  ! Loop over the number of records required

    nrec = GetNumberOfRecords(n, mrecsize)
    do irec = 1, nrec
  ! Get record details
      call GetRecordDetails(irec, n, mrecsize, il, iu, ninrec)
  ! Write out the record
      write(e_fileunit) (a(j), j = il, iu)
    enddo

  endif

end subroutine WriteBlockC

! *************************************************************************** !

subroutine WriteHeader(itype, mnem, n)
  !
  ! Write out a block header of type itype
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscInt, intent(in) :: itype
  character(len = *), intent(in) :: mnem
  character(len = 8) :: zmnem8
  character(len = 4 ) :: ztype4

  PetscInt, intent(in) :: n
  integer(kind = int32) :: n4

10 format(1X, "'", A8, "'", 1X, I11, 1X, "'", A4, "'")

  zmnem8 = mnem

  n4 = n

  ztype4 = '    '
  if (itype == e_typeI) ztype4 = 'INTE'
  if (itype == e_typeS) ztype4 = 'REAL'
  if (itype == e_typeD) ztype4 = 'DOUB'
  if (itype == e_typeB) ztype4 = 'LOGI'
  if (itype == e_typeC) ztype4 = 'CHAR'

  if (e_formatted) then
    write(e_fileunit, 10) zmnem8, n, ztype4
  else
    write(e_fileunit) zmnem8, n4, ztype4
  endif

end subroutine WriteHeader

! *************************************************************************** !

subroutine CmpToCNOBuf(buff, arr, porv, flip, conv)
  !
  ! Compress an array to compressed natural order
  ! This is x-fastest, z-slowest, with inactive cells removed
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscReal, intent(out) :: buff(:)
  PetscReal, intent(in)  :: arr (:)
  PetscReal, intent(in)  :: porv(:)
  PetscReal, intent(in)  :: flip
  PetscReal, intent(in)  :: conv

  PetscInt :: ig, ia, ix, iy, ize, izp

  ! Loop in natural order, skipping inactive cells and holding a count

  ia = 0
  do ize = 1, e_nz
    izp = e_nz-ize+1
    do iy = 1, e_ny
      do ix = 1, e_nx
        ig = e_nx*e_ny*(izp-1) + e_nx*(iy-1) + ix
        if (porv(ig)>0.0) then
          ia = ia+1
          buff(ia) = flip*conv*arr(ig)
        endif
      enddo
    enddo
  enddo

end subroutine CmpToCNOBuf

! *************************************************************************** !

subroutine CopyToBufferI(a, ibuf, il, iu)
  !
  ! Copy values to an int*32 buffer
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscInt, intent(in) :: a(:)
  integer(kind = int32), intent(out) :: ibuf(:)
  PetscInt, intent(in) :: il, iu

  PetscInt :: i

  do i = il, iu
    ibuf(i-il+1) = a(i)
  enddo

end subroutine CopyToBufferI

! *************************************************************************** !

subroutine CopyToBufferS(a, sbuf, il, iu)
  !
  ! Copy values to an real*32 buffer
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscReal, intent(in) :: a(:)
  real(kind = real32), intent(out) :: sbuf(:)
  PetscInt, intent(in) :: il, iu

  PetscInt :: i

  do i = il, iu
    sbuf(i-il+1) = real(a(i), real32)
  enddo

end subroutine CopyToBufferS

! *************************************************************************** !

subroutine CopyToBufferD(a, dbuf, il, iu)
  !
  ! Copy values to an real*64 buffer
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscReal, intent(in) :: a(:)
  real(kind = real64), intent(out) :: dbuf(:)
  PetscInt, intent(in) :: il, iu

  PetscInt :: i

  do i = il, iu
    dbuf(i-il+1) = real(a(i), real64)
  enddo

end subroutine CopyToBufferD

! *************************************************************************** !

subroutine CopyToBufferB(a, ibuf, il, iu)
  !
  ! Copy bool to an int*32 buffer (false->0;true->1)
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscBool, intent(in) :: a(:)
  integer(kind = int32), intent(out) :: ibuf(:)
  PetscInt, intent(in) :: il, iu

  PetscInt :: i

  do i = il, iu
    if (a(i)) then
      ibuf(i-il+1) = 1
    else
      ibuf(i-il+1) = 0
    endif
  enddo

end subroutine CopyToBufferB

! *************************************************************************** !

subroutine GetRecordDetails(iRec, n, mInRec, il, iu, nInRec)
  !
  ! Given the number of values to be written, get details of (iRec)th record
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscInt, intent(in ) :: iRec, n, mInRec
  PetscInt, intent(out) :: il, iu, nInRec

  PetscInt :: iumax

  iumax = n
  il = mInRec*(iRec-1)+1
  iu = il+mInRec-1
  if (iu>iumax) iu = iumax
  nInRec = iu-il+1

end subroutine GetRecordDetails

! *************************************************************************** !

function GetNumberOfRecords(n, mrecsize)
  !
  ! Given the number of values, and the record size, return number of records
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscInt :: GetNumberOfRecords

  PetscInt, intent(in) :: n, mrecsize

  GetNumberOfRecords = (n+mrecsize-1)/mrecsize

end function GetNumberOfRecords

! *************************************************************************** !

subroutine SetInteHead(intehead)
  !
  ! Set up the 'intehead' header record
  ! This contains a range of problem dimensions and values
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscInt :: intehead(:)

  intehead( 3) = 1           ! Metric units
  intehead( 9) = e_nx         ! x-dimension
  intehead(10) = e_ny         ! y-dimension
  intehead(11) = e_nz         ! z-dimension
  intehead(12) = e_na         ! Cell count

  intehead(17) = e_nwell      ! Number of wells
  intehead(18) = e_ncwmax     ! Max completions/well
  intehead(20) = e_nwelmx     ! Max wells/group
  intehead(21) = e_ngroup     ! Number of groups

  intehead(25) = e_niwel      ! Number ints/well
  intehead(26) = e_nswel      ! Number reals/well
  intehead(27) = e_nxwel      ! Number solution values/well
  intehead(28) = e_nzwel      ! Number char*8/well

  intehead(33) = e_nicon      ! Number ints/completion
  intehead(34) = e_nscon      ! Number reals/completion
  intehead(35) = e_nxcon      ! Number solution values/completion

  intehead(37) = e_nigrp      ! Number ints/group
  intehead(38) = e_nsgrp      ! Number reals/group
  intehead(39) = e_nxgrp      ! Number solution values/group
  intehead(40) = e_nzgrp      ! Number char*8/group

  intehead(65) = 1            ! Date day
  intehead(66) = 1            ! Date month
  intehead(67) = 2000         ! Date year
  intehead(95) = -1           ! Not Eclipse writing this

end subroutine SetInteHead

! *************************************************************************** !

subroutine ReleaseEwriterBuffers()
  !
  ! Release buffers held by Output_Eclipse_module and
  ! close summary and restart files
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  call DeleteAtoC()

  if (e_ltoap_allocated) then
    deallocate(e_ltocp )
    deallocate(e_nlmaxp)
    e_ltoap_allocated = PETSC_FALSE
  endif

  if (e_opened) then
    close(unit = UNIT_SUMM_WRITE)
    close(unit = UNIT_REST_WRITE)
  endif

end subroutine ReleaseEwriterBuffers

! *************************************************************************** !

subroutine SetupRestMaps(ltoa, option, nlmax, mlmax)
  !
  ! Set up the local to compressed natural mapping by proc
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  use Option_module

  implicit none

  PetscInt :: ltoa(:)
  type(option_type) :: option
  PetscInt, intent(in) ::nlmax, mlmax

  PetscInt :: iproct, nproc, ioproc, iproco, lproco, &
              ibuf1(1), nbuf1, nlmaxo, &
              liproct, lioproc, il, ia, ic
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE), ierr, itag

  ! Set up useful scalars

  itag  = 1
  nbuf1 = 1
  ibuf1(1) = nlmax

  e_nlmax = nlmax
  e_mlmax = mlmax

  ierr = 0

  nproc   = option%mycommsize
  iproct  = option%myrank
  ioproc  = option%io_rank
  liproct = iproct+1
  lioproc = ioproc+1

  ! For I/O proc, store own values and receive from others
  ! For other procs, set to I/O proc

  if (option%myrank == option%io_rank) then

    ! Allocate the all-proc arrays

    allocate(e_ltocp (mlmax, nproc))
    allocate(e_nlmaxp(       nproc))

    ! Store this-proc values

    e_ltoap_allocated = PETSC_TRUE
    do il = 1, nlmax
      ia =   ltoa(il)
      ic = e_atoc(ia)
      e_ltocp (il, liproct) = ic
    enddo
    e_nlmaxp(liproct) = nlmax

    ! Receive from other procs

    do iproco = 0, nproc-1
      if (iproco /= option%io_rank) then

        ! Receive nlmax value from other proc

        call MPI_Recv(ibuf1, nbuf1, MPI_INTEGER, iproco, MPI_ANY_TAG, &
                      option%mycomm, status_mpi, ierr)

        ! Store other-proc nlmax value

        lproco = iproco+1
        nlmaxo = ibuf1(1)
        e_nlmaxp(lproco) = nlmaxo

        ! Receive ltoa map from other proc (temporary store in ltoa)

        call MPI_Recv(ltoa, nlmaxo, MPI_INTEGER, iproco, MPI_ANY_TAG, &
                      option%mycomm, status_mpi, ierr)

        ! Copy ltoa into all-proc array

        do il = 1, nlmaxo
          ia =   ltoa(il)
          ic = e_atoc(ia)
          e_ltocp (il, lproco) = ic
        enddo

      endif
    enddo
  else

    ! Send to the IO proc

    ibuf1(1) = nlmax
    call MPI_Send(ibuf1, nbuf1, MPI_INTEGER, ioproc, itag, option%mycomm, ierr)
    call MPI_Send(ltoa , nlmax, MPI_INTEGER, ioproc, itag, option%mycomm, ierr)

  endif

  ! No need for e_atoc now, so delete it

  call DeleteAtoC()

end subroutine SetupRestMaps

! *************************************************************************** !

function GetMlmax()
  !
  ! Get the maximum value of nlmax over all ranks
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscInt :: GetMlmax

  GetMlMax = e_mlmax

end function GetMlmax

! *************************************************************************** !

subroutine GetYMDHMMS(tdays, years, months, days, hours, mins, microsecs)
  !
  ! For a time since 2000, find the date.
  ! A somewhat approximate calendar, as it
  ! ignores the special (100 and 400 year) leap years
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscReal, intent(in)  :: tdays
  PetscInt , intent(out) :: years, months, days, hours, mins, microsecs

  PetscReal :: remhours, remmins, remmsecs, tdib, tdiy, tend
  PetscInt  :: dim ! Days in month

  ! Days/month normal and leap years

  !                        Ja  Fe  Ma  Ap  Ma  Ju  Jl  Au  Sp  Oc  No  De
  PetscInt, parameter :: dimn(12) = &
                         (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  PetscInt, parameter :: diml(12) = &
                         (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

  ! Days in each 4-year leap year block

  PetscReal, parameter :: dp4y = 1461.0
  PetscReal, parameter :: dp3y = 1095.0
  PetscReal, parameter :: dp2y =  730.0
  PetscReal, parameter :: dp1y =  365.0

  PetscReal, parameter :: epst =  0.000001

  PetscInt  :: n4y, yinb, imon
  PetscBool :: is_leap

  ! Find number of 4-year blocks, and subtract this to get days in block

  n4y  = int(tdays/dp4y)
  tdib = tdays-n4y*dp4y

  ! If in the last year of a block, is a leap year

  is_leap = PETSC_FALSE
  if (tdib > dp3y) is_leap = PETSC_TRUE

  ! Set up days in year, subtract off 1, 2 or 3
  ! normal years to find time in year

  tdiy = tdib
  yinb = 0
  if (tdib > dp3y) then
    tdiy = tdib - dp3y
    yinb = 3
  else if (tdib > dp2y) then
    tdiy = tdib - dp2y
    yinb = 2
  else if (tdib > dp1y) then
    tdiy = tdib - dp1y
    yinb = 1
  endif

  ! Assume 1st month, go through year and subtract days in each elapsed month

  months = 1
  days   = int(tdiy)+1

  tend = 0.0
  do imon = 1, 11
    if (is_leap) then
      dim = diml(imon)
    else
      dim = dimn(imon)
    endif
    tend = tend+dim
    if (tdiy >= tend) then
      months = imon+1
      days   = days-dim
    endif
  enddo
  if (days > 31) days = 31

  ! Set up years

  years = 4*n4y+yinb+2000

  remhours = 24.0      *(tdays   -int(tdays   ))
  remmins  = 60.0      *(remhours-int(remhours))
  remmsecs = 60000000.0*(remmins -int(remmins ))

  hours     = int(remhours)
  mins      = int(remmins )
  microsecs = int(remmsecs)

end subroutine GetYMDHMMS

! *************************************************************************** !

subroutine DeleteAtoC()
  !
  ! Delete the active to compressed active map
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  if (e_atoc_allocated) then
    deallocate(e_atoc)
    e_atoc_allocated = PETSC_FALSE
  endif

end subroutine DeleteAtoC

! *************************************************************************** !

subroutine ThrowEwriterException(message)
  !
  ! Throw a serious error (should never be called)
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  character(len = *) :: message
  print *, message
  stop
end subroutine ThrowEwriterException

end module Output_Eclipse_module
