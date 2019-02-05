module Grid_Grdecl_module

!  Module to read a grid using industry-standard Eclipse keyword syntax
!  and convert it to a Pflotran explicit unstructured grid

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Grid_Unstructured_Aux_module
  use Option_module
  use Input_Aux_module

  implicit none

  private

  ! ASCII codes for tab, line feed and carriage return

  PetscInt, parameter :: g_ictab = 9
  PetscInt, parameter :: g_iclf  = 10
  PetscInt, parameter :: g_iccr  = 13

  ! Cartesian frame pointers

  PetscInt, parameter :: g_xdir = 1
  PetscInt, parameter :: g_ydir = 2
  PetscInt, parameter :: g_zdir = 3

  PetscInt, parameter :: g_ndir = 3

  ! Master flag indicating grdecl system is being used

  PetscBool :: g_is_grdecl = PETSC_FALSE

  ! Problem dimensions (set defaults for single cell)

  PetscInt :: g_nx     = 1
  PetscInt :: g_ny     = 1
  PetscInt :: g_nxp    = 2
  PetscInt :: g_nyp    = 2
  PetscInt :: g_nz     = 1
  PetscInt :: g_nxy    = 1
  PetscInt :: g_nxpnyp = 4
  PetscInt :: g_nxyz   = 1

  PetscBool :: g_dimens_read = PETSC_FALSE

  ! z-flip from Eclipse to Pflotran convention

  PetscBool :: g_iscpg = PETSC_FALSE
  PetscReal :: z_flip  = -1.0

  PetscBool :: g_cpgallocated = PETSC_FALSE
  PetscBool :: g_isnewtran    = PETSC_FALSE

  ! Counters for active cells, connections and max connections.

  PetscInt :: g_na = 1
  PetscInt :: g_nc = 1
  PetscInt :: g_mc = 1

  ! Grid-sized arrays

  PetscReal, pointer :: g_coord(:) => null()
  PetscReal, pointer :: g_zcorn(:) => null()

  PetscReal, pointer :: g_dx  (:) => null()
  PetscReal, pointer :: g_dy  (:) => null()
  PetscReal, pointer :: g_dz  (:) => null()

  PetscReal, pointer :: g_kx  (:) => null()
  PetscReal, pointer :: g_ky  (:) => null()
  PetscReal, pointer :: g_kz  (:) => null()

  PetscReal, pointer :: g_mx  (:) => null()
  PetscReal, pointer :: g_my  (:) => null()
  PetscReal, pointer :: g_mz  (:) => null()

  PetscReal, pointer :: g_tops(:) => null()
  PetscReal, pointer :: g_poro(:) => null()
  PetscReal, pointer :: g_ntg (:) => null()

  PetscReal, pointer :: g_xloc(:) => null()
  PetscReal, pointer :: g_yloc(:) => null()
  PetscReal, pointer :: g_zloc(:) => null()

  PetscInt , pointer :: g_actn(:) => null()

  PetscInt , pointer :: g_gtoa(:) => null()
  PetscInt , pointer :: g_atog(:) => null()
  PetscInt , pointer :: g_atoc(:) => null()

  ! Active sized arrays

  PetscReal, pointer :: g_bvol(:) => null()

  PetscReal, pointer :: g_x   (:) => null()
  PetscReal, pointer :: g_y   (:) => null()
  PetscReal, pointer :: g_z   (:) => null()

  PetscInt , pointer :: g_cia (:) => null()
  PetscInt , pointer :: g_cja (:) => null()

  PetscReal, pointer :: g_ccx  (:) => null()
  PetscReal, pointer :: g_ccy  (:) => null()
  PetscReal, pointer :: g_ccz  (:) => null()
  PetscReal, pointer :: g_carea(:) => null()

  ! Column pointer for reader

  PetscInt :: g_column = 1

  character(len = MAXSTRINGLENGTH) :: g_error_string = 'OK'
  PetscInt :: g_error_flag = 0

  ! Public access to this module

  public  :: SetIsGrdecl
  public  :: GetIsGrdecl
  public  :: UGrdEclExplicitRead
  public  :: SetUGrdEclCmplLocation
  public  :: UGrdEclWellCmplCleanup
  public  :: FindWellIndex
  public  :: GetGrdNCmpl
  public  :: GetCmplData
  public  :: DeallocatePoroPermArrays
  public  :: GetPoroPermValues
  public  :: WriteStaticDataAndCleanup
  public  :: PermPoroExchangeAndSet

  private :: GrdeclReader

  ! Type containing well location data

  type, public :: well_locn_type

    private

    character(len = MAXSTRINGLENGTH), public :: w_name ! well name

    PetscInt :: ikl
    PetscInt :: iku

  end type well_locn_type

  ! Counters for well location instances

  PetscInt :: g_nwell_data = 0
  PetscInt :: g_mwell_data = 0

  ! Array of well location instances

  type(well_locn_type), allocatable :: g_well_data(:)

  ! Type containing completion location data

  type, public :: cmpl_data_type

    private

    PetscInt  :: ci
    PetscInt  :: cj
    PetscInt  :: ck
    PetscInt  :: ik
    PetscInt  :: iw
    PetscInt  :: ig
    PetscInt  :: ia

    PetscReal :: dx
    PetscReal :: dy
    PetscReal :: dz

    PetscReal :: z

  end type cmpl_data_type

  ! Counters for completion location instances

  PetscInt :: g_ncmpl_data = 0
  PetscInt :: g_mcmpl_data = 0

  ! Array of completion location instances

  type(cmpl_data_type), allocatable :: g_cmpl_data(:)

contains

! ************************************************************************** !

subroutine SetIsGrdecl()
  !
  ! Set flag indicating that a combined grid/well Jacobian is required
  !
  ! Author: Dave Ponting
  ! Date: 01/30/18
  !
  implicit none

  g_is_grdecl = PETSC_TRUE

end subroutine SetIsGrdecl

! ************************************************************************** !

function GetIsGrdecl()
  !
  ! Get flag indicating that a combined grid/well Jacobian is being used
  !
  ! Author: Dave Ponting
  ! Date: 01/30/18
  !
  implicit none

  PetscBool :: GetIsGrdecl

  GetIsGrdecl = g_is_grdecl

end function GetIsGrdecl

! ************************************************************************** !

subroutine UGrdEclExplicitRead(unstructured_grid, filename, option)
  !
  ! Reads an Eclgrid file, stores the data, finds connections
  ! and divides them across the ranks
  !
  ! Author: Dave Ponting
  ! Date: 11/05/18
  !
  use String_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid 
  type(unstructured_explicit_type), pointer :: explicit_grid
  character(len = MAXSTRINGLENGTH) :: filename
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscInt :: fileid, ierr

  ! Set the grid pointer

  explicit_grid => unstructured_grid%explicit_grid

  ! Read the grdecl file on i/o rank, generate cells and connections

  ierr = 0

  ! Create the input system

  fileid = UNIT_GRDECL_READ
  input => InputCreate(fileid, filename, option)

  ! Read on io_rank only

  if (option%myrank == option%io_rank) then
    call GrdeclReader(input, option)
  endif
  call MPI_Bcast(g_error_flag, ONE_INTEGER_MPI, MPI_INTEGER, &
                 option%io_rank, option%mycomm, ierr)
  if (g_error_flag>0) then
    input%ierr = 1
    call MPI_Bcast(g_error_string, MAXSTRINGLENGTH, MPI_CHARACTER, &
                   option%io_rank, option%mycomm, ierr)
    call InputErrorMsg(input, option, 'GRDECL read error', g_error_string)
  endif

  ! Destroy the input system

  call InputDestroy(input)

  ! Sycn after read

  call MPI_Barrier(option%mycomm, ierr)

  ! Distribute the well data

  call DistributeWells(option)

  ! Distribute and store the cell data

  call DistributeCells(explicit_grid, option)

  ! Distribute and store the connection data

  call DistributeConnections(explicit_grid, option)

  ! Distribute and store the poro and perm data

  call DistributePoroPerm(option)

  ! Set up the cell locations

  if (option%myrank == option%io_rank) then
    call CreateElements(unstructured_grid, explicit_grid)
  endif

end subroutine UGrdEclExplicitRead

!*****************************************************************************!

subroutine WriteStaticDataAndCleanup(write_ecl,eclipse_options,option)
  !
  ! If Eclipse files are required, output grid and init files
  ! Then clean up all the static data arrays which are no longer needed
  ! Author: Dave Ponting
  ! Date: 11/05/18
  !
  use Output_Aux_module, only : output_option_eclipse_type
  use Output_Eclipse_module, only: WriteEclipseFilesGrid, &
                                   WriteEclipseFilesInit, &
                                   SelectFormattedFiles

  implicit none

  PetscBool, intent(in) :: write_ecl
  type(output_option_eclipse_type), pointer :: eclipse_options
  type(option_type) :: option

  PetscInt :: iw, nw, nctw, mcpw

  character(len=MAXSTRINGLENGTH) :: efilename

  ! Output the Eclipse grid and init files

  if (option%myrank == option%io_rank) then
    if (write_ecl) then
      nw = g_nwell_data
      mcpw = 1
      do iw = 1, nw
        nctw = GetGrdNCmpl(iw)
        if (nctw > mcpw) mcpw = nctw
      enddo
      efilename = trim(option%output_file_name_prefix)
      if (eclipse_options%write_ecl_form) call SelectFormattedFiles()
      call WriteEclipseFilesGrid(efilename, g_nx, g_ny, g_nz, &
                                 g_coord, g_zcorn, g_gtoa, nw, mcpw)
      call WriteEclipseFilesInit(g_kx, g_ky, g_kz, g_mx, g_my, g_mz, &
                                 g_zloc, g_poro, g_ntg, g_bvol, g_gtoa, g_atoc)
    endif
  endif

  ! Free allocated memory

  call DeallocateGridArrays  ()
  call DeallocateActiveArrays()

end subroutine WriteStaticDataAndCleanup

!*****************************************************************************!

subroutine SetUGrdEclCmplLocation(wname, ci, cj, ckuser, cijk_d, qerr)
  !
  ! Set up a completion for well wname at location (ci, cj, ck)
  !
  ! Author: Dave Ponting
  ! Date: 11/21/18
  !

  implicit none

  character(len = *), intent(in)    :: wname
  PetscInt        , intent(in)    :: ci, cj, ckuser
  PetscBool       , intent(in)    :: cijk_d
  PetscBool       , intent(inout) :: qerr

  PetscInt :: ck, ik, iw
  PetscBool :: found

  ! Set up Pflotran ck value

  ck = ckuser
  if (cijk_d) ck = -ckuser

  ! Find next completion location

  ik = g_ncmpl_data+1

  ! Have we heard of this well? Allocate if not

  found = findWellIndex(wname, iw)
  if (.not. found) then
    ! Not found, so assume is new
    iw = g_nwell_data+1
    if (iw > g_mwell_data) then
      call ExtendWellNames()
    endif
    g_well_data(iw)%w_name = wname
    g_well_data(iw)%ikl    = ik   ! First completion
    g_nwell_data           = iw
  else
    ! Found, so should be current WELL_DATA block
    if (iw /= g_nwell_data) then
      call setError('Same WELL_DATA well has appeared more than once')
      qerr = PETSC_TRUE
    endif
  endif

  g_well_data(iw)%iku = ik !  Latest completion

  if (ik > g_mcmpl_data) then
    call ExtendCmplData()
  endif

  g_cmpl_data(ik)%ci = ci
  g_cmpl_data(ik)%cj = cj
  g_cmpl_data(ik)%ck = ck
  g_cmpl_data(ik)%iw = iw
  g_cmpl_data(ik)%ik = ik
  g_cmpl_data(ik)%ig = -1
  g_cmpl_data(ik)%ia = -1

  g_ncmpl_data = ik

end subroutine SetUGrdEclCmplLocation

!*****************************************************************************!

subroutine ExtendWellNames()
  !
  ! Extend the g_well_data structure
  !
  ! Author: Dave Ponting
  ! Date: 11/21/18
  !

  implicit none

  type(well_locn_type), allocatable :: temp_well_data(:)

  PetscInt :: mwell_data_new, iw

  ! If existing well locations, copy them

  if (g_nwell_data > 0) then

    allocate(temp_well_data(g_nwell_data))

    do iw = 1, g_nwell_data
      call CopyWellData(temp_well_data(iw), g_well_data(iw))
    enddo

  endif

  ! Deallocate, increase and re-allocate the well data store

  if (allocated(g_well_data)) deallocate(g_well_data)

  mwell_data_new = g_mwell_data + 10
  allocate(g_well_data(mwell_data_new))

  ! Copy over the  existing well locations

  do iw = 1, g_nwell_data
    call CopyWellData(g_well_data(iw), temp_well_data(iw))
  enddo

  ! Deallocate the temporary well locations

  if (allocated(temp_well_data)) deallocate(temp_well_data)

  ! Store the new dimension

  g_mwell_data = mwell_data_new

end subroutine ExtendWellNames

! ************************************************************************** !

subroutine ExtendCmplData()
  !
  ! Extend the g_cmpl_data structure
  !
  ! Author: Dave Ponting
  ! Date: 11/21/18
  !

  implicit none

  type(cmpl_data_type), allocatable :: temp_cmpl_data(:)

  PetscInt :: mcmpl_data_new, ik

  ! If existing well locations, copy them

  if (g_ncmpl_data > 0) then

    allocate(temp_cmpl_data(g_ncmpl_data))

    do ik = 1, g_ncmpl_data
      call CopyCmplData(temp_cmpl_data(ik), g_cmpl_data(ik))
    enddo

  endif

  ! Deallocate, increase and re-allocate the well locations store

  if (allocated(g_cmpl_data)) deallocate(g_cmpl_data)
  mcmpl_data_new = g_mcmpl_data + 10
  allocate(g_cmpl_data(mcmpl_data_new))

  ! Copy over the  existing well locations

  do ik = 1, g_ncmpl_data
    call CopyCmplData(g_cmpl_data(ik), temp_cmpl_data(ik))
  enddo

  ! Deallocate the temporary well locations

  if (allocated(temp_cmpl_data)) deallocate(temp_cmpl_data)

  ! Store the new dimension

  g_mcmpl_data = mcmpl_data_new

end subroutine ExtendCmplData

!*****************************************************************************!

subroutine UGrdEclWellCmplCleanup
  !
  ! Clean up the g_well_data and g_cmpl_data structures
  !
  ! Author: Dave Ponting
  ! Date: 11/21/18
  !

  implicit none

  ! g_well_data

  if (g_mwell_data > 0) then
    if (allocated(g_well_data)) then
      deallocate(g_well_data)
      g_mwell_data = 0
    endif
  endif

  ! g_cmpl_data

  if (g_mcmpl_data > 0) then
    if (allocated(g_cmpl_data)) then
      deallocate(g_cmpl_data)
      g_mcmpl_data = 0
    endif
  endif

end subroutine UGrdEclWellCmplCleanup

!*****************************************************************************!

subroutine GrdeclReader(input, option)
  !
  ! Reads an Eclgrid file on I/O proc
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  PetscBool, parameter :: isdep_no   = PETSC_FALSE
  PetscBool, parameter :: isdep_yes  = PETSC_TRUE
  PetscBool, parameter :: isperm_no  = PETSC_FALSE
  PetscBool, parameter :: isperm_yes = PETSC_TRUE

  PetscInt  :: nx, ny, nz, ierr
  PetscBool :: qerr
  character(len = MAXWORDLENGTH) :: word
  character(len = MAXSTRINGLENGTH) :: zmess

  ierr = 0
  qerr = PETSC_FALSE

  ! Read through items in the grdecl file

  do

    call InputReadPflotranString(input, option)
    if (input%ierr /= 0) exit                   ! Detect eof and leave

  ! Keyword found

    call InputReadWord(input, option, word, PETSC_TRUE)
    call CheckError(input, 'ECLGRD subkeyword', qerr);if ( qerr ) exit
    call StringToUpper(word)

    select case(trim(word))
      case('DIMENS')
        call InputReadPflotranString(input, option)
        call CheckError(input, 'Data following DIMENS', qerr);if (qerr) exit
        call InputReadInt(input, option, nx)
        call CheckError(input, 'NX in DIMENS', qerr);if (qerr) exit
        call InputReadInt(input, option, ny)
        call CheckError(input, 'NY in DIMENS', qerr);if (qerr) exit
        call InputReadInt(input, option, nz)
        call CheckError(input, 'NZ in DIMENS', qerr);if (qerr) exit
        call SetDimens(nx, ny, nz)
      case('COORD')
        call IsCPG()
        call checkDimensRead(qerr);if (qerr) exit
        call ReadECoordArray(g_coord, ierr, input, option, qerr)
      case('ZCORN')
        call IsCPG()
        call checkDimensRead(qerr);if (qerr) exit
        call ReadEZcornArray(g_zcorn, ierr, input, option, qerr)
      case('DX')
        call ReadEGridArrayR(g_dx, 'DX', ierr, input, option, &
                             isdep_no, isperm_no, qerr)
        if (qerr) then
          call SetError('DX in GRDECL');exit
        endif
      case('DY')
        call ReadEGridArrayR(g_dy, 'DY', ierr, input, option, &
                             isdep_no, isperm_no, qerr)
        if (qerr) then
          call SetError('DY in GRDECL');exit
        endif
      case('DZ')
        call ReadEGridArrayR(g_dz, 'DZ', ierr, input, option, &
                             isdep_no, isperm_no, qerr)
        if (qerr) then
          call SetError('DZ in GRDECL');exit
        endif
      case('PERMX')
        call ReadEGridArrayR(g_kx, 'PERMX', ierr, input, option, &
                             isdep_no, isperm_yes, qerr)
        if (qerr) then
          call SetError('PERMX in GRDECL');exit
        endif
      case('PERMY')
        call ReadEGridArrayR(g_ky, 'PERMY', ierr, input, option, &
                             isdep_no, isperm_yes, qerr)
        if (qerr) then
          call SetError('PERMY in GRDECL');exit
        endif
      case('PERMZ')
        call ReadEGridArrayR(g_kz, 'PERMZ', ierr, input, option, &
                             isdep_no, isperm_yes, qerr)
        if (qerr) then
          call SetError('PERMZ in GRDECL');exit
        endif
      case('MULTX')
        call ReadEGridArrayR(g_kx, 'MULTX', ierr, input, option, &
                             isdep_no, isperm_yes, qerr)
        if (qerr) then
          call SetError('MULTX in GRDECL');exit
        endif
      case('MULTY')
        call ReadEGridArrayR(g_ky, 'MULTY', ierr, input, option, &
                             isdep_no, isperm_yes, qerr)
        if (qerr) then
          call SetError('MULTY in GRDECL');exit
        endif
      case('MULTZ')
        call ReadEGridArrayR(g_kz, 'MULTZ', ierr, input, option, &
                             isdep_no, isperm_yes, qerr)
        if (qerr) then
          call SetError('MULTZ in GRDECL');exit
        endif
      case('PORO')
        call ReadEGridArrayR(g_poro, 'PORO', ierr, input, option, &
                             isdep_no, isperm_no, qerr)
        if (qerr) then
          call SetError('PORO in GRDECL');exit
        endif
      case('TOPS')
        call ReadEGridArrayR(g_tops, 'TOPS', ierr, input, option, &
                             isdep_yes, isperm_no, qerr)
        if (qerr) then
          call SetError('TOPS in GRECL');exit
        endif
      case('NTG')
        call ReadEGridArrayR(g_ntg, 'NTG', ierr, input, option, &
                             isdep_no, isperm_no, qerr)
        if (qerr) then
          call SetError('NTG in GRECL');exit
        endif
      case('ACTNUM')
        call ReadEGridArrayI(g_actn, 'ACTNUM', ierr, input, option, qerr)
        if (qerr) then
          call SetError('ACTNUM in GRECL');exit
        endif
      case('NEWTRAN')
        g_isnewtran   = PETSC_TRUE
      case('OLDTRAN')
        g_isnewtran   = PETSC_FALSE
      case('/') ! Isolated un-used terminator on its own line is harmless
      case default
        qerr  = PETSC_TRUE
        zmess = 'GRDECL sub-keyword ' // trim(word) // ' not recognised'
        call SetError(zmess)
        if (qerr) exit
    end select
  enddo

  ! Process the grid data read

  call ProcessGridData()

  ! Process the well data read

  call ProcessWellData(qerr)

end subroutine GrdeclReader

! ************************************************************************** !

subroutine CheckError(input, zerr, qerr)
  !
  ! Check the error code on the input stream
  ! Set an error message and return qerr as true if an error is found
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  type(input_type), pointer :: input
  character(len = *):: zerr
  PetscBool, intent(inout) :: qerr

  PetscInt :: ierr

  ierr = input%ierr

  if (ierr == 0) then
    qerr = PETSC_FALSE
  else
    qerr = PETSC_TRUE
    call SetError(zerr)
  endif

end subroutine CheckError

! ************************************************************************** !

subroutine SetError(zerr)
  !
  ! Set an error message and flag
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  character(len = *) :: zerr

  g_error_flag = 1
  g_error_string = zerr

end subroutine SetError

!*****************************************************************************!

subroutine ProcessGridData()
  !
  ! Process the data found in the grdecl file
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18
  !

  implicit none

  PetscInt  :: ix, iy, iz, izp, ize, ig, ia, iact, ng, ic
  PetscReal :: porv, bvol, eps

  ! Pore volume for rock-filled cell

  eps = 1.0e-6

  ! Fill in missing data values

  if (g_iscpg) then

    ! If corner point data supplied, extract cell dimensions like dx

    call extractCellDimensionsAndLocationsFromCPG()
  else

    ! Fill in x and y locations for each z-layer

    do iz = 1, g_nz
     call FillXYPositionsForLayer(iz)
    enddo

    ! Fill in tops and z locations for each ix/iy column

    do ix = 1, g_nx
      do iy = 1, g_ny
        call FillZPositionsForColumn(ix, iy)
      enddo
    enddo

    ! Allocate and fill in coord and zcorn

    if (.not.g_cpgallocated) then
      allocate(g_coord(6*g_nxpnyp));g_coord =  0.0
      allocate(g_zcorn(8*g_nxyz  ));g_zcorn =  0.0
      g_cpgallocated = PETSC_TRUE
    endif

    call ExtractCPGFromCellDimensionsAndLocations()

  endif

  ! Find pore volumes and set up active order

  g_na = 0
  do ix = 1, g_nx
    do iy = 1, g_ny
      do iz = 1, g_nz

        ig = GetNaturalIndex(ix, iy, iz)

        ! Find pore volume

        bvol = 0.0
        porv = 0.0

        ! Uses non-dip dx.dy.dz volume calculation

        iact = g_actn(ig)

        bvol = g_dx(ig)*g_dy(ig)*g_dz(ig)
        if (iact == 1) then
          porv = g_poro(ig)*g_ntg(ig)*bvol
        endif
        if (iact == 2) then
          porv = bvol
        endif
        if (iact == 3) then
          porv = eps
        endif

        if (porv > 0.0) then
          g_na = g_na+1
          g_gtoa(ig) = g_na
        endif

      enddo
    enddo
  enddo

  allocate(g_atog(g_na));g_atog = -1
  allocate(g_atoc(g_na));g_atoc = -1

  ! Set up a to g

  ng = g_nx*g_ny*g_nz
  do ig = 1, ng
    ia = g_gtoa(ig)
    if (ia>-1) then
      g_atog(ia) = ig
    endif
  enddo

  ! Set up atoc (loop over cells in Eclipse compresed natural order)

  ic = 0
  do ize = 1, g_nz
    izp = g_nz-ize+1
    do iy = 1, g_ny
      do ix = 1, g_nx
        ig = GetNaturalIndex(ix, iy, izp)
        ia = g_gtoa(ig)
        if (ia>-1) then
          ic = ic+1
          g_atoc(ia) = ic
        endif
      enddo
    enddo
  enddo

  ! Initial estmate for number of connections
  ! Sufficient for normal connections + some space for wells

  g_mc = 3*g_na
  g_nc = 0

  ! Allocate active arrays

  call AllocateActiveArrays()

  ! Find the cell locations and connections

  call GenerateGridConnections()

end subroutine ProcessGridData

!*****************************************************************************!

subroutine ProcessWellData(qerr)
  !
  ! Process the well data found in the WELL_DATA sections
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18
  !

  ! Generate connections for well completions

  implicit none

  PetscInt :: iw, ix, iy, iz, ia, ja, ig
  PetscInt :: ik, jk, jw, ikl, iku, ncmpl_active
  PetscBool, intent(inout) :: qerr

  character(len = MAXSTRINGLENGTH) :: wname, zmess

  ! Loop over wells found, loading completion data into single big array

  outer:do iw = 1, g_nwell_data

    wname = g_well_data(iw)%w_name

    ikl = g_well_data(iw)%ikl
    iku = g_well_data(iw)%iku

    do ik = ikl, iku

      ix = g_cmpl_data(ik)%ci
      iy = g_cmpl_data(ik)%cj
      iz = g_cmpl_data(ik)%ck

      ! Deal with negative iz case (input using CIJK_D)

      if (iz < 0) then
        iz = g_nz+iz+1
        g_cmpl_data(ik)%ck = iz
      endif

      ! Check for out-of-range well locations

      if ((ix < 1) .or. (ix > g_nx)) then
         zmess = 'Completion I-location out of range, well ' // trim(wname)
         call SetError(zmess)
         qerr = PETSC_TRUE
         exit outer
      endif

      if ((iy < 1) .or. (iy > g_ny)) then
         zmess = 'Completion J-location out of range, well ' // trim(wname)
         call SetError(zmess)
         qerr = PETSC_TRUE
         exit outer
      endif

      if ((iz < 1) .or. (iz > g_nz)) then
         zmess = 'Completion K-location out of range, well ' // trim(wname)
         call SetError(zmess)
         qerr = PETSC_TRUE
         exit outer
      endif

      ig = GetNaturalIndex(ix, iy, iz)
      ia = g_gtoa(ig)

      g_cmpl_data(ik)%ia = ia
      g_cmpl_data(ik)%ig = ig
      g_cmpl_data(ik)%iw = iw

      g_cmpl_data(ik)%dx = g_dx(ig)
      g_cmpl_data(ik)%dy = g_dy(ig)
      g_cmpl_data(ik)%dz = g_dz(ig)

      if (ia>-1) then
        g_cmpl_data(ik)%z = g_z (ia)
      else
        g_cmpl_data(ik)%z = 0.0
      endif

    enddo 

  enddo outer

  if (.not.qerr) then

    ! Scan the array and compress out the inactive cells

    ncmpl_active = 0
    do ik = 1, g_ncmpl_data

      ia = g_cmpl_data(ik)%ia
      if (ia > -1) then
        ncmpl_active = ncmpl_active + 1
        if (ncmpl_active < ik) then
          call CopyCmplData( g_cmpl_data(ncmpl_active), g_cmpl_data(ik) )
        endif
      endif

    enddo

    g_ncmpl_data = ncmpl_active

    ! Reset the lower and upper pointers into the
    ! compressed completion array for each well

    do iw = 1, g_nwell_data

      ! Default range will zero-trip

      ikl =  0
      iku = -1

      do ik = 1, g_ncmpl_data
        jw = g_cmpl_data(ik)%iw
        if (iw .eq. jw) then
          if (ikl == 0) ikl = ik
                        iku = ik
        endif
      enddo

      g_well_data(iw)%ikl = ikl
      g_well_data(iw)%iku = iku

    enddo

    ! Now add the well connections to the list of connections

    do iw = 1, g_nwell_data

      ikl = g_well_data(iw)%ikl
      iku = g_well_data(iw)%iku

      do ik = ikl, iku

        ia = g_cmpl_data(ik)%ia

        do jk = ik+1, iku

          ja = g_cmpl_data(jk)%ia

          g_nc = g_nc + 1
          if (g_nc == g_mc) call ReallocateConnectionArrays()
          g_cia  (g_nc) = ia
          g_cja  (g_nc) = ja
          g_ccx  (g_nc) = 0.5*(g_x(ia)+g_x(ja))
          g_ccy  (g_nc) = 0.5*(g_y(ia)+g_y(ja))
          g_ccz  (g_nc) = 0.5*(g_z(ia)+g_z(ja))
          g_carea(g_nc) = 0.0

        enddo
      enddo

    enddo

  endif

end subroutine ProcessWellData

!*****************************************************************************!

subroutine DistributeWells(option)
  !
  ! Distribute the well information from the I/O ranks to other ranks
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  type(option_type) :: option

  PetscInt :: iw, ik
  PetscMPIInt :: ibuf(12)
  PetscReal   :: rbuf(12)

  ! First, clean up any existing well and cmpl structures on non IO proc

  if (option%myrank /= option%io_rank) then
    deallocate(g_cmpl_data)
  endif

  call BroadcastInt(g_ncmpl_data, option)

  if (option%myrank /= option%io_rank) then
    if (g_ncmpl_data>0) then
      allocate(g_cmpl_data(g_ncmpl_data))
    endif
  endif

  do iw = 1, g_nwell_data
    call BroadcastInt(g_well_data(iw)%ikl, option)
    call BroadcastInt(g_well_data(iw)%iku, option)
  enddo

  ibuf = 0
  rbuf = 0.0
  do ik = 1, g_ncmpl_data

    if (option%myrank == option%io_rank) then

      ibuf(1) = g_cmpl_data(ik)%ci
      ibuf(2) = g_cmpl_data(ik)%cj
      ibuf(3) = g_cmpl_data(ik)%ck
      ibuf(4) = g_cmpl_data(ik)%ik
      ibuf(5) = g_cmpl_data(ik)%iw
      ibuf(6) = g_cmpl_data(ik)%ig
      ibuf(7) = g_cmpl_data(ik)%ia

      rbuf(1) = g_cmpl_data(ik)%dx
      rbuf(2) = g_cmpl_data(ik)%dy
      rbuf(3) = g_cmpl_data(ik)%dz
      rbuf(4) = g_cmpl_data(ik)%z

    endif

    call BroadcastIntN (ibuf, option)
    call BroadcastRealN(rbuf, option)

    if (option%myrank /= option%io_rank) then

      g_cmpl_data(ik)%ci = ibuf(1)
      g_cmpl_data(ik)%cj = ibuf(2)
      g_cmpl_data(ik)%ck = ibuf(3)
      g_cmpl_data(ik)%ik = ibuf(4)
      g_cmpl_data(ik)%iw = ibuf(5)
      g_cmpl_data(ik)%ig = ibuf(6)
      g_cmpl_data(ik)%ia = ibuf(7)

      g_cmpl_data(ik)%dx = rbuf(1)
      g_cmpl_data(ik)%dy = rbuf(2)
      g_cmpl_data(ik)%dz = rbuf(3)
      g_cmpl_data(ik)%z  = rbuf(4)

    endif

  enddo

end subroutine DistributeWells

!*****************************************************************************!

subroutine DistributeCells(explicit_grid, option)
  !
  ! Split the generated connections between ranks
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !

  implicit none

  type(unstructured_explicit_type), pointer :: explicit_grid
  type(option_type) :: option

  PetscInt :: ierr, ial, nal, nals, rem, ia, irank, narank, iabase
  PetscReal, pointer :: temp_real_array(:,:)
  PetscMPIInt :: int_mpi
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)

  ! Distribute the global active count

  call BroadcastInt(g_na, option)

  explicit_grid%num_cells_global = g_na

  ! Divide cells across ranks

  call GetLocalCount(g_na, nal, nals, rem, option)

  ! Allocate cells

  call AllocateCells(nal, explicit_grid)

  ! Now do the send/recv and load

  if (option%myrank == option%io_rank) then

    allocate(temp_real_array(5, nals+1))

    iabase = 0
    do irank = 0, option%mycommsize-1

      narank = nals
      if (irank < rem ) narank = narank + 1

      if (irank == option%io_rank) then

        ! if the cells reside on io_rank

        do ial = 1, narank
          ia = ial + iabase
          explicit_grid%cell_ids      (ial)   =        ia
          explicit_grid%cell_centroids(ial)%x = g_x   (ia)
          explicit_grid%cell_centroids(ial)%y = g_y   (ia)
          explicit_grid%cell_centroids(ial)%z = g_z   (ia)
          explicit_grid%cell_volumes  (ial)   = g_bvol(ia)
        enddo
      else

        ! if the cells reside on another rank

        int_mpi = narank*5
        do ial = 1, narank
          ia = ial + iabase
          temp_real_array(1, ial) = real  (ia)
          temp_real_array(2, ial) = g_x   (ia)
          temp_real_array(3, ial) = g_y   (ia)
          temp_real_array(4, ial) = g_z   (ia)
          temp_real_array(5, ial) = g_bvol(ia)
        enddo
        call MPI_Send(temp_real_array, int_mpi, MPI_DOUBLE_PRECISION, irank, &
                      nal, option%mycomm, ierr)
      endif
      iabase = iabase + narank
    enddo
  else
    ! other ranks post the recv
    allocate(temp_real_array(5, nal))
    int_mpi = nal*5
    call MPI_Recv(temp_real_array, int_mpi, &
                  MPI_DOUBLE_PRECISION, option%io_rank, &
                  MPI_ANY_TAG, option%mycomm, status_mpi, ierr)
    do ia = 1, nal
      explicit_grid%cell_ids      (ia)   = int(temp_real_array(1, ia))
      explicit_grid%cell_centroids(ia)%x =     temp_real_array(2, ia)
      explicit_grid%cell_centroids(ia)%y =     temp_real_array(3, ia)
      explicit_grid%cell_centroids(ia)%z =     temp_real_array(4, ia)
      explicit_grid%cell_volumes  (ia)   =     temp_real_array(5, ia)
    enddo

  endif

  deallocate(temp_real_array)

end subroutine DistributeCells

!*****************************************************************************!

subroutine DistributeConnections(explicit_grid, option)
  !
  ! Split the generated connections between ranks
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !

  implicit none

  type(unstructured_explicit_type), pointer :: explicit_grid
  type(option_type) :: option

  PetscInt :: ierr, rem, irank, icbase, ic, icl, ncl, ncls, ncrank
  PetscReal, pointer :: temp_real_array(:,:)
  PetscMPIInt :: int_mpi
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)

  ! Distribute the global connection count

  call BroadcastInt(g_nc, option)

  ! Divide connections across ranks

  call GetLocalCount(g_nc, ncl, ncls, rem, option)

  ! Allocate connections

  call AllocateConnections(ncl, explicit_grid)

  ! for now, read all cells from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_real_array(6, ncls+1))
    ! read for other processors
    icbase = 0
    do irank = 0, option%mycommsize-1
      temp_real_array = UNINITIALIZED_DOUBLE
      ncrank = ncls
      if (irank < rem) ncrank = ncrank + 1

      ! if the cells reside on io_rank
      if (irank == option%io_rank) then

        do icl = 1, ncl
          ic = icl + icbase
          explicit_grid%connections   (1, icl)   = g_cia  (ic)
          explicit_grid%connections   (2, icl)   = g_cja  (ic)
          explicit_grid%face_centroids(   icl)%x = g_ccx  (ic)
          explicit_grid%face_centroids(   icl)%y = g_ccy  (ic)
          explicit_grid%face_centroids(   icl)%z = g_ccz  (ic)
          explicit_grid%face_areas   (    icl)   = g_carea(ic)
        enddo
      else
        ! otherwise communicate to other ranks
        int_mpi = ncrank*6
        do icl = 1, ncrank
          ic = icl + icbase
          temp_real_array(1, icl) = real(g_cia  (ic))
          temp_real_array(2, icl) = real(g_cja  (ic))
          temp_real_array(3, icl) =      g_ccx  (ic)
          temp_real_array(4, icl) =      g_ccy  (ic)
          temp_real_array(5, icl) =      g_ccz  (ic)
          temp_real_array(6, icl) =      g_carea(ic)
        enddo
        call MPI_Send(temp_real_array, int_mpi, MPI_DOUBLE_PRECISION, irank, &
                      ncrank, option%mycomm, ierr)
      endif
      icbase = icbase + ncrank
    enddo
  else
    ! other ranks post the recv
    allocate(temp_real_array(6, ncl))
    int_mpi = ncl*6
    call MPI_Recv(temp_real_array, int_mpi, &
                  MPI_DOUBLE_PRECISION, option%io_rank, &
                  MPI_ANY_TAG, option%mycomm, status_mpi, ierr)
    do icl = 1, ncl
      explicit_grid%connections   (1, icl)   = int(temp_real_array(1, icl))
      explicit_grid%connections   (2, icl)   = int(temp_real_array(2, icl))
      explicit_grid%face_centroids(   icl)%x =     temp_real_array(3, icl)
      explicit_grid%face_centroids(   icl)%y =     temp_real_array(4, icl)
      explicit_grid%face_centroids(   icl)%z =     temp_real_array(5, icl)
      explicit_grid%face_areas    (   icl)   =     temp_real_array(6, icl)
    enddo

  endif

  deallocate(temp_real_array)

end subroutine DistributeConnections

! *************************************************************************** !

subroutine DistributePoroPerm(option)
  !
  ! Split the porosity and permeability between ranks
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !

  implicit none

  type(option_type) :: option

  call BroadcastInt(g_nxyz, option)
  call BroadcastInt(g_na  , option)

end subroutine DistributePoroPerm

!*****************************************************************************!

subroutine CreateElements(unstructured_grid, explicit_grid)

  use Grid_Grdecl_Util_module, only : getCorners

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(unstructured_explicit_type), pointer :: explicit_grid

  PetscInt, parameter :: nVertPerCell = 8

  PetscInt :: num_vertices, icorn, ivert, ia, ig, &
              ix, iy, iz, iox, ioy, ioz

  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3), xw(3)

  num_vertices = nVertPerCell*g_na

  explicit_grid%num_elems = g_na
  explicit_grid%num_vertices = num_vertices
  unstructured_grid%max_nvert_per_cell = nVertPerCell

  allocate(explicit_grid%cell_vertices(0:unstructured_grid% &
                                       max_nvert_per_cell, g_na))
  allocate(explicit_grid%vertex_coordinates(explicit_grid%num_vertices))

  do ivert = 1, explicit_grid%num_vertices
    explicit_grid%vertex_coordinates(ivert)%x = 0.0
    explicit_grid%vertex_coordinates(ivert)%y = 0.0
    explicit_grid%vertex_coordinates(ivert)%z = 0.0
  enddo

  ivert = 0
  do ia = 1, g_na

    ig = g_atog(ia)
    call getCellCoordinates(ig, ix, iy, iz)

    ! Get the corners

    call GetCorners( ix, iy, iz, &
                     x000, x100, x010, x110, &
                     x001, x101, x011, x111, &
                     g_coord, g_zcorn, g_nx, g_ny )

    explicit_grid%cell_vertices(0, ia) = nVertPerCell

    do iox = 1, 2
      do ioy = 1, 2
        do ioz = 1, 2

          if (iox == 1 .and. ioy == 1 .and. ioz == 1) xw = x000
          if (iox == 2 .and. ioy == 1 .and. ioz == 1) xw = x100
          if (iox == 1 .and. ioy == 2 .and. ioz == 1) xw = x010
          if (iox == 2 .and. ioy == 2 .and. ioz == 1) xw = x110
          if (iox == 1 .and. ioy == 1 .and. ioz == 2) xw = x001
          if (iox == 2 .and. ioy == 1 .and. ioz == 2) xw = x101
          if (iox == 1 .and. ioy == 2 .and. ioz == 2) xw = x011
          if (iox == 2 .and. ioy == 2 .and. ioz == 2) xw = x111

          ivert = ivert + 1

          icorn = 0
          if (ioz == 2) icorn = 4
          if ((iox == 1) .and. (ioy == 1)) icorn = icorn + 1
          if ((iox == 2) .and. (ioy == 1)) icorn = icorn + 2
          if ((iox == 2) .and. (ioy == 2)) icorn = icorn + 3
          if ((iox == 1) .and. (ioy == 2)) icorn = icorn + 4

          explicit_grid%cell_vertices(icorn, ia) = ivert

          explicit_grid%vertex_coordinates(ivert)%x = xw(g_xdir)
          explicit_grid%vertex_coordinates(ivert)%y = xw(g_ydir)
          explicit_grid%vertex_coordinates(ivert)%z = xw(g_zdir)

        enddo
      enddo
    enddo
  enddo

  explicit_grid%output_mesh_type = CELL_CENTERED_OUTPUT_MESH

end subroutine CreateElements

!*****************************************************************************!

subroutine SetDimens(nx, ny, nz)
  !
  ! Dimens read, so allocate the grid-sized arrays
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !

  implicit none
  PetscInt, intent(in) :: nx, ny, nz

  g_nx     = nx
  g_ny     = ny
  g_nz     = nz
  g_nxp    = nx + 1
  g_nyp    = ny + 1
  g_nxy    = nx*ny
  g_nxpnyp = g_nxp*g_nyp
  g_nxyz   = nx*ny*nz

  g_dimens_read = PETSC_TRUE

  call AllocateGridArrays()

end subroutine SetDimens

!*****************************************************************************!

subroutine CheckDimensRead(qerr)
  !
  ! Check that Dimens has been read, throw error if not
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !

  implicit none

  PetscBool, intent(inout) :: qerr

  qerr = PETSC_FALSE
  if (.not. g_dimens_read) then
    qerr = PETSC_TRUE
  endif

end subroutine CheckDimensRead

!*****************************************************************************!

subroutine AllocateGridArrays()
  !
  ! Allocate basic grid arrays:
  ! cell sizes, perms, mults, poro, ntg, actn, cell locations
  ! These arrays have the full grid dimension (nx.ny.nz)
  ! Also allocate map top active order
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !

  implicit none

  allocate(g_dx (g_nxyz  ));g_dx  =  0.0
  allocate(g_dy (g_nxyz  ));g_dy  =  0.0
  allocate(g_dz (g_nxyz  ));g_dz  =  0.0

  allocate(g_kx (g_nxyz  ));g_kx  =  0.0
  allocate(g_ky (g_nxyz  ));g_ky  =  0.0
  allocate(g_kz (g_nxyz  ));g_kz  =  0.0

  allocate(g_mx (g_nxyz  ));g_mx  =  1.0
  allocate(g_my (g_nxyz  ));g_my  =  1.0
  allocate(g_mz (g_nxyz  ));g_mz  =  1.0

  allocate(g_tops(g_nxyz));g_tops = UNINITIALIZED_DOUBLE ! To trap unset layers
  allocate(g_poro(g_nxyz));g_poro =  0.0
  allocate(g_ntg (g_nxyz));g_ntg  =  1.0 ! Net to gross ration (default to 1)

  allocate(g_xloc(g_nxyz ));g_xloc = 0.0
  allocate(g_yloc(g_nxyz ));g_yloc = 0.0
  allocate(g_zloc(g_nxyz ));g_zloc = 0.0

  allocate(g_actn(g_nxyz));g_actn =  1 ! Actnum (used to set cells inactive)
  allocate(g_gtoa(g_nxyz));g_gtoa = -1

end subroutine AllocateGridArrays

!*****************************************************************************!

subroutine AllocateActiveArrays()
  !
  ! Allocate active grid arrays (cell x/y/z position and bulk volume)
  ! Also allocate connection data (active cells connected, centroids, area)
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  allocate(g_bvol(g_na));g_bvol   = 0.0

  allocate(g_x   (g_na));g_x      = 0.0
  allocate(g_y   (g_na));g_y      = 0.0
  allocate(g_z   (g_na));g_z      = 0.0

  allocate(g_cia (g_mc));g_cia    = -1
  allocate(g_cja (g_mc));g_cja    = -1

  allocate(g_ccx  (g_mc));g_ccx   = 0.0
  allocate(g_ccy  (g_mc));g_ccy   = 0.0
  allocate(g_ccz  (g_mc));g_ccz   = 0.0
  allocate(g_carea(g_mc));g_carea = 0.0

end subroutine AllocateActiveArrays

!*****************************************************************************!

subroutine ReallocateConnectionArrays()
  !
  ! Reallocate the connection arrays (used to extend thee arrays)
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Utility_module, only: ReallocateArray

  implicit none

  PetscInt :: mc

  ! Note keep resetting mc as Reallocate keeps increasing it

  mc = g_mc;call ReallocateArray(g_cia  , mc)
  mc = g_mc;call ReallocateArray(g_cja  , mc)

  mc = g_mc;call ReallocateArray(g_ccx  , mc)
  mc = g_mc;call ReallocateArray(g_ccy  , mc)
  mc = g_mc;call ReallocateArray(g_ccz  , mc)
  mc = g_mc;call ReallocateArray(g_carea, mc)

  g_mc = mc

end subroutine ReAllocateConnectionArrays

!*****************************************************************************!

subroutine GenerateGridConnections
  !
  ! Given the cell locations, generate the connections
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt  :: ix, iy, iz, ig, ia, jg, ja
  PetscReal :: x, y, z, dx, dy, dz, mx, my, mz

  ! Set up grid connections

  do ix = 1, g_nx
    do iy = 1, g_ny
      do iz = 1, g_nz

        ig = GetNaturalIndex(ix, iy, iz)
        ia = g_gtoa(ig)
        if (ia>0) then

          dx = g_dx(ig)
          dy = g_dy(ig)
          dz = g_dz(ig)

          mx = g_mx(ig)
          my = g_my(ig)
          mz = g_mz(ig)

          g_bvol(ia) = dx*dy*dz

          g_x(ia) = g_xloc(ig)
          g_y(ia) = g_yloc(ig)
          g_z(ia) = g_zloc(ig)

          x = g_xloc(ig)
          y = g_yloc(ig)
          z = g_zloc(ig)

          ! x-direction

          if (ix < g_nx) then
            jg = GetNaturalIndex(ix+1, iy, iz)
            ja = g_gtoa(jg)
            if (ja >0) then
              g_nc = g_nc + 1
              if (g_nc == g_mc) call ReallocateConnectionArrays()
              g_cia  (g_nc) = ia
              g_cja  (g_nc) = ja
              g_ccx  (g_nc) = x + 0.5*dx
              g_ccy  (g_nc) = y
              g_ccz  (g_nc) = z
              g_carea(g_nc) = mx*dy*dz
            endif
          endif

          ! y-direction

          if (iy< g_ny) then
            jg = GetNaturalIndex(ix, iy+1, iz)
            ja = g_gtoa(jg)
            if (ja > 0) then
              g_nc = g_nc + 1
              if (g_nc == g_mc) call ReallocateConnectionArrays()
              g_cia  (g_nc) = ia
              g_cja  (g_nc) = ja
              g_ccx  (g_nc) = x
              g_ccy  (g_nc) = y + 0.5*dy
              g_ccz  (g_nc) = z
              g_carea(g_nc) = my*dz*dx
            endif
          endif

          ! z-direction

          if (iz < g_nz) then
            jg = GetNaturalIndex(ix, iy, iz+1)
            mz = g_mz(jg)
            ja = g_gtoa(jg)
            if (ja > 0) then
              g_nc = g_nc + 1
              if (g_nc == g_mc) call ReallocateConnectionArrays()
              g_cia  (g_nc) = ia
              g_cja  (g_nc) = ja
              g_ccx  (g_nc) = x
              g_ccy  (g_nc) = y
              g_ccz  (g_nc) = z + 0.5*dz
              g_carea(g_nc) = mz*dx*dy
            endif
          endif

        endif

      enddo ! enddo iz
    enddo ! enddo iy
  enddo ! enddo ix

end subroutine GenerateGridConnections

!*****************************************************************************!

subroutine DeallocateGridArrays()
  !
  ! Deallocate the grid arrays
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Utility_module, only : DeallocateArray

  implicit none

  if (g_cpgallocated) then
    call DeallocateArray(g_coord)
    call DeallocateArray(g_zcorn)
    g_cpgallocated = PETSC_FALSE
  endif

  call DeallocateArray(g_dx  )
  call DeallocateArray(g_dy  )
  call DeallocateArray(g_dz  )

  call DeallocateArray(g_tops)
  call DeallocateArray(g_ntg )

  call DeallocateArray(g_xloc)
  call DeallocateArray(g_yloc)
  call DeallocateArray(g_zloc)

  call DeallocateArray(g_actn)
  call DeallocateArray(g_gtoa)

end subroutine DeallocateGridArrays

!*****************************************************************************!

subroutine DeallocateActiveArrays
  !
  ! Deallocate the active arrays
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Utility_module, only : DeallocateArray

  implicit none

  call DeallocateArray(g_bvol)

  call DeallocateArray(g_x   )
  call DeallocateArray(g_y   )
  call DeallocateArray(g_z   )

  call DeallocateArray(g_mx  )
  call DeallocateArray(g_my  )
  call DeallocateArray(g_mz  )

  call DeallocateArray(g_cia  )
  call DeallocateArray(g_cja  )
  call DeallocateArray(g_ccx  )
  call DeallocateArray(g_ccy  )
  call DeallocateArray(g_ccz  )
  call DeallocateArray(g_carea)

  call DeallocateArray(g_atoc )

end subroutine DeallocateActiveArrays

!*****************************************************************************!

subroutine DeallocatePoroPermArrays(option)
  !
  ! Deallocate the perma and porosity arrays
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Utility_module, only : DeallocateArray

  implicit none

  type(option_type) :: option

  if (option%myrank == option%io_rank) then
    call DeallocateArray(g_atog)
    call DeallocateArray(g_kx  )
    call DeallocateArray(g_ky  )
    call DeallocateArray(g_kz  )
    call DeallocateArray(g_poro)
  endif

end subroutine DeallocatePoroPermArrays

!*****************************************************************************!

subroutine GetPoroPermValues(ia, poro, permx, permy, permz)
  !
  ! Get the perm and porosity values of a given active cell
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt , intent(in)  :: ia
  PetscReal, intent(out) :: poro, permx, permy, permz

  PetscInt :: ig

  ig = g_atog(ia)

  poro  = g_poro(ig)
  permx = g_kx  (ig)
  permy = g_ky  (ig)
  permz = g_kz  (ig)

end subroutine GetPoroPermValues

!*****************************************************************************!

function GetNaturalIndex(ix, iy, iz)
  !
  ! Get the natural cell index (in ix/iy/iz order with ix varying fastest)
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt :: GetNaturalIndex
  PetscInt, intent(in) :: ix, iy, iz

  GetNaturalIndex = g_nxy*(iz-1) + g_nx*(iy-1) + ix

end function GetNaturalIndex

!*****************************************************************************!

subroutine GetCellCoordinates(ig, ix, iy, iz)
  !
  ! For a given natural cell index, extract the ix, iy, iz location
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(in ) :: ig
  PetscInt, intent(out) :: ix, iy, iz

  PetscInt :: ir

  iz = (ig-1)/g_nxy+1
  ir = ig-(iz-1)*g_nxy
  iy = (ir-1)/g_nx+1
  ix = ir-(iy-1)*g_nx

end subroutine GetCellCoordinates

!*****************************************************************************!

subroutine BroadcastInt(ival, option)
  !
  ! Send a PetscInt value from the I/O rank to all the others
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(inout) :: ival
  type(option_type) :: option
  PetscMPIInt :: int_mpi
  PetscInt    :: ierr

  ierr    = 0
  int_mpi = ival

  call MPI_Bcast(int_mpi, ONE_INTEGER_MPI, MPI_INTEGER, option%io_rank, &
                 option%mycomm, ierr)

  ival = int_mpi

end subroutine BroadcastInt

!*****************************************************************************!

subroutine BroadcastIntN(ival, option)
  !
  ! Send a PetscInt array from the I/O rank to all the others
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscMPIInt, intent(inout) :: ival(:)
  type(option_type) :: option

  PetscInt    :: n, ierr
  PetscMPIInt :: nmpi

  ierr = 0
  n    = size(ival)
  nmpi = n

  call MPI_Bcast(ival, nmpi, MPI_INTEGER, option%io_rank, &
                 option%mycomm, ierr)

end subroutine BroadcastIntN

!*****************************************************************************!

subroutine BroadcastRealN(rval, option)
  !
  ! Send a PetscReal array from the I/O rank to all the others
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscReal, intent(inout) :: rval(:)
  type(option_type) :: option

  PetscInt    :: n, ierr
  PetscMPIInt :: nmpi

  ierr = 0
  n    = size(rval)
  nmpi = n

  call MPI_Bcast( rval, nmpi, MPI_DOUBLE_PRECISION, &
                  option%io_rank, option%mycomm, ierr )

end subroutine BroadcastRealN

!*****************************************************************************!

subroutine GetLocalCount(ng, nl, nls, rem, option)
  !
  ! Find the local count size when splitting an array over all the ranks
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(in ) :: ng
  PetscInt, intent(out) :: nl
  PetscInt, intent(out) :: nls
  PetscInt, intent(out) :: rem
  type(option_type)     :: option

  nl  = ng/option%mycommsize
  nls = nl

  rem = ng - nl*option%mycommsize
  if (option%myrank < rem) nl = nl + 1

end subroutine GetLocalCount

!*****************************************************************************!

subroutine AllocateCells(nal, explicit_grid)
  !
  ! Allocate a set of arrays to hold cell properties
  ! Used when distributing cells over ranks
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(in) :: nal
  type(unstructured_explicit_type), pointer :: explicit_grid
  PetscInt :: ial

  allocate(explicit_grid%cell_ids(nal))
  explicit_grid%cell_ids = 0
  allocate(explicit_grid%cell_volumes(nal))
  explicit_grid%cell_volumes = 0
  allocate(explicit_grid%cell_centroids(nal))
  do ial = 1, nal
    explicit_grid%cell_centroids(ial)%x = 0.d0
    explicit_grid%cell_centroids(ial)%y = 0.d0
    explicit_grid%cell_centroids(ial)%z = 0.d0
  enddo

end subroutine AllocateCells

!*****************************************************************************!

subroutine AllocateConnections(ncl, explicit_grid)
  !
  ! Allocate a set of arrays to hold connection properties
  ! Used when distributing connections over ranks
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(in) :: ncl
  type(unstructured_explicit_type), pointer :: explicit_grid

  PetscInt :: icl

  allocate(explicit_grid%connections(2, ncl))
  explicit_grid%connections = 0
  allocate(explicit_grid%face_areas(ncl))
  explicit_grid%face_areas = 0
  allocate(explicit_grid%face_centroids(ncl))
  do icl = 1, ncl
    explicit_grid%face_centroids(icl)%x = 0.d0
    explicit_grid%face_centroids(icl)%y = 0.d0
    explicit_grid%face_centroids(icl)%z = 0.d0
  enddo

end subroutine AllocateConnections

!*****************************************************************************!

subroutine FillXYPositionsForLayer(iz)
  !
  ! When building a grid from cell size data (dx,dy,dz) set up an iz-layer
  ! by filling in the x and y locations
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(in) :: iz

  PetscInt :: ix, iy, ig
  PetscReal :: sum, dx, dy

  ! Fill in x-locations for each y-row

  do iy = 1, g_ny
    sum = 0.0
    do ix = 1, g_nx
      ig  = GetNaturalIndex(ix, iy, iz)
      dx  = g_dx(ig)
      g_xloc(ig) = sum + 0.5*dx
      sum        = sum +     dx
    enddo
  enddo

  ! Fill in y-locations for each x-row

  do ix = 1, g_nx
    sum = 0.0
    do iy = 1, g_ny
      ig  = GetNaturalIndex(ix, iy, iz)
      dy  = g_dy(ig)
      g_yloc(ig) = sum + 0.5*dy
      sum        = sum +     dy
    enddo
  enddo

end subroutine FillXYPositionsForLayer

! *************************************************************************** !

subroutine FillZPositionsForColumn(ix, iy)
  !
  ! When building a grid from cell-size data (dx,dy,dz)
  ! set up an (ix,iy)-column by filling in the z locations
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(in) :: ix, iy

  PetscInt :: iz, izabove, izbelow, ig, igabove, igbelow
  PetscReal :: dz

  ! First pass for tops (starting at top, going down)

  do izabove = g_nz, 2, -1
     izbelow = izabove-1

     igabove = GetNaturalIndex(ix, iy, izabove)
     igbelow = GetNaturalIndex(ix, iy, izbelow)

     if (Uninitialized(g_tops(igbelow))) then
       g_tops(igbelow) = g_tops(igabove) - g_dz(igabove)
     endif

   enddo

  ! Second pass for z-locations

   do iz = 1, g_nz
     ig = GetNaturalIndex(ix, iy, iz)
     dz = g_dz(ig)
     g_zloc(ig) = g_tops(ig) - 0.5*dz
   enddo

end subroutine FillZPositionsForColumn

! *************************************************************************** !

subroutine ExtractCellDimensionsAndLocationsFromCPG()
  !
  ! Extract cell dimensions and locations from COORD/ZCORN data
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  use Grid_Grdecl_Util_module, only : GetCorners

  implicit none

  PetscInt  :: ix, iy, iz, ig
  PetscReal :: vdx, vdy, vdz, vpx, vpy, vpz

  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3)

  ! For each cell, find base offsets in the
  ! (2.nx).(2.ny).(2.nz) corner z-val array

  do ix = 1, g_nx
    do iy = 1, g_ny
      do iz = 1, g_nz

        ! Get the corners

        call GetCorners( ix, iy, iz, x000, x100, x010, x110, &
                         x001, x101, x011, x111, &
                         g_coord, g_zcorn, g_nx, g_ny )

        ! Given the cell corner locations,
        ! find the cell dimensions and locations

        call GetHexDims( x000, x100, x010, x110, &
                         x001, x101, x011, x111, &
                         vdx, vdy, vdz, vpx, vpy, vpz )

        ! Find cell index and store locations

        ig = GetNaturalIndex(ix, iy, iz)

        g_dx  (ig) = vdx
        g_dy  (ig) = vdy
        g_dz  (ig) = vdz

        g_xloc(ig) = vpx
        g_yloc(ig) = vpy
        g_zloc(ig) = vpz

      enddo
    enddo
  enddo

end subroutine ExtractCellDimensionsAndLocationsFromCPG

! *************************************************************************** !

subroutine ExtractCPGFromCellDimensionsAndLocations()
  !
  ! Extract corner point locations from cell size data
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscReal :: xsum, xl, xu, ysum, yl, yu, dx, dy, dz, xc, yc, zc, xn, yn, zn
  PetscInt  :: ix, iy, iz, ig, ixnb, iynb, iznb, iox, ioy, ioz, jn

  ! Do the coords, looping over columns (use dx,dy in layer 1)

  iz = 1

  xsum = 0.0
  do ix = 1, g_nx

    iy = 1 ! (Use first line of y-values in finding x-spacing)
    ig = GetNaturalIndex(ix, iy, iz)
    xl = xsum
    xu = xsum+g_dx(ig)
    xsum = xu

    ysum = 0.0

    do iy = 1, g_ny

      ig = GetNaturalIndex(ix, iy, iz)
      yl = ysum
      yu = ysum+g_dy(ig)
      ysum = yu

                                       call setCoordLine(ix  , iy  , xl, yl)
      if (ix == g_nx)                  call setCoordLine(ix+1, iy  , xu, yl)
      if (iy == g_ny)                  call setCoordLine(ix  , iy+1, xl, yu)
      if (ix == g_nx .and. iy == g_ny) call setCoordLine(ix+1, iy+1, xu, yu)

    enddo
  enddo

  ! Now the zcorns

  do ix = 1, g_nx
    do iy = 1, g_ny
      do iz = 1, g_nz

        ig = GetNaturalIndex(ix, iy, iz)

        xc = g_xloc(ig)
        yc = g_yloc(ig)
        zc = g_zloc(ig)

        dx = g_dx  (ig)
        dy = g_dy  (ig)
        dz = g_dz  (ig)

        ixnb = 2*(ix-1)
        iynb = 2*(iy-1)
        iznb = 2*(iz-1)

        do iox = 0, 1
          xn = xc + 0.5*(2*iox-1)*dx
          do ioy = 0, 1
            yn = yc + 0.5*(2*ioy-1)*dy
            do ioz = 0, 1
              zn = zc + 0.5*(2*ioz-1)*dz
              jn = 4*g_nxy*(iznb+ioz)+2*g_nx*(iynb+ioy)+ixnb+iox+1
              g_zcorn(jn) = zn
            enddo
          enddo
        enddo

      enddo
    enddo
  enddo

end subroutine ExtractCPGFromCellDimensionsAndLocations

! *************************************************************************** !

subroutine SetCoordLine(ixp, iyp, x, y)
  !
  ! Set up a vertical coordinate line (corner of a pillar of cells)
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscInt , intent(in) :: ixp, iyp
  PetscReal, intent(in) :: x, y

  PetscReal, parameter :: zl = -10000.0
  PetscReal, parameter :: zu =  10000.0

  PetscInt :: ibase

  ibase = 6*((iyp-1)*g_nxp + (ixp-1))

  g_coord(ibase+1) = x
  g_coord(ibase+2) = y
  g_coord(ibase+3) = zl
  g_coord(ibase+4) = x
  g_coord(ibase+5) = y
  g_coord(ibase+6) = zu

end subroutine SetCoordLine

! *************************************************************************** !

subroutine GetHexDims( x000, x100, x010, x110, &
                       x001, x101, x011, x111, &
                       vdx, vdy, vdz, vpx, vpy, vpz )
  !
  ! Given cell corner locations, extract cell dimensions and locations
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscReal, intent(in ) :: x000(3), x100(3), x010(3), x110(3), &
                            x001(3), x101(3), x011(3), x111(3)
  PetscReal, intent(out) :: vdx, vdy, vdz, vpx, vpy, vpz

  PetscInt  :: idir
  PetscReal :: dx(3), dy(3), dz(3)

  ! Cell centre is just the average

  vpx = 0.125*( x000(g_xdir)+x100(g_xdir)+x010(g_xdir)+x110(g_xdir) &
               +x001(g_xdir)+x101(g_xdir)+x011(g_xdir)+x111(g_xdir) )

  vpy = 0.125*( x000(g_ydir)+x100(g_ydir)+x010(g_ydir)+x110(g_ydir) &
               +x001(g_ydir)+x101(g_ydir)+x011(g_ydir)+x111(g_ydir) )

  vpz = 0.125*( x000(g_zdir)+x100(g_zdir)+x010(g_zdir)+x110(g_zdir) &
               +x001(g_zdir)+x101(g_zdir)+x011(g_zdir)+x111(g_zdir) )

  ! Cell dimensions are distances between faces: form x,y and z intervals

  do idir = 1, g_ndir

                     ! 1**        0**
    dx(idir) = 0.25*( x100(idir)-x000(idir) &
                     +x110(idir)-x010(idir) &
                     +x101(idir)-x001(idir) &
                     +x111(idir)-x011(idir) )

                     ! *1*        *0*
    dy(idir) = 0.25*( x010(idir)-x000(idir) &
                     +x110(idir)-x100(idir) &
                     +x011(idir)-x001(idir) &
                     +x111(idir)-x101(idir) )

                     ! **1        **0
    dz(idir) = 0.25*( x001(idir)-x000(idir) &
                     +x101(idir)-x100(idir) &
                     +x011(idir)-x010(idir) &
                     +x111(idir)-x110(idir) )

  enddo

  ! Find scalar distances

  vdx = GetScalarLength(dx)
  vdy = GetScalarLength(dy)
  vdz = GetScalarLength(dz)

end subroutine GetHexDims

! *************************************************************************** !

function GetScalarLength(v)
  !
  ! Get the scalar length of a three-vector
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscReal :: GetScalarLength
  PetscReal, intent(in) :: v(g_ndir)
  PetscReal :: sum

  sum = v(g_xdir)*v(g_xdir) &
       +v(g_ydir)*v(g_ydir) &
       +v(g_zdir)*v(g_zdir)

  GetScalarLength = sqrt(sum)

end function GetScalarLength

! *************************************************************************** !

subroutine GetNextWord(word, exitTime, input, option)
  !
  ! Get the next white-space-delimited item on the input file
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  character(len = *), intent(out  ) :: word
  PetscBool, intent(out) :: exitTime
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len = MAXSTRINGLENGTH) :: line
  character(len = MAXWORDLENGTH  ) :: test

  character c
  PetscInt :: i, l, j, startcol, ic
  PetscBool :: lastCharRead, somethingRead, iws

  PetscBool :: started

  started       = PETSC_FALSE
  somethingRead = PETSC_FALSE

  line    = input%buf
  l       = len(line)
  word = ' '
  j = 0
  exitTime = PETSC_FALSE

  ! Loop over attempts to find some characters

  do

  ! Rest of current line

    lastCharRead = PETSC_FALSE

    startcol = g_column
    do i = startcol, l
      c = line(i:i)
      iws = IsWhiteSpace(c)
      ic = ichar(c)
      if ((.not. started) .and. (.not.iws)) started = PETSC_TRUE
      if (       started  .and. ( iws .or. (i == l) )) exit
      if (started) then
        j = j + 1
        word(j:j) = c
        somethingRead = PETSC_TRUE
      endif
      g_column = i + 1
      if ((i == l) .or. (ic == g_iclf)) lastCharRead = PETSC_TRUE
    enddo

  ! Action on end of line: either exit or read a new line

    if (lastCharRead) then
      if (somethingRead) then
        ! If something read, is terminator
        exit
      else
        ! If nothing read, read another line
        call InputReadPflotranString(input, option)
        if (InputCheckExit(input, option)) exitTime = PETSC_FALSE
        g_column = 1
        line     = input%buf
        l        = len(line)
      endif
    endif

    if (exitTime .or. somethingRead) exit

  enddo

  ! Check anything found is right-aligned and look for / teminator

  if (somethingRead  .and. (.not. exitTime)) then

    c = word(1:1)
    if (IsWhiteSpace(c)) then
      test = adjustl(word)
      word = test
    endif

  ! Check for / terminator

    c = word(1:1)
    if (c == '/') then
      exitTime = PETSC_TRUE
    endif

  endif

end subroutine GetNextWord

! *************************************************************************** !

function IsWhiteSpace(c)
  !
  ! Indicates that a character is white space
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscBool :: isWhiteSpace
  character, intent(in) :: c
  PetscInt :: ic

  IsWhiteSpace = PETSC_FALSE

  ic = ichar(c)
  if (( c == ' '    ) .or. &
      (ic == g_ictab) .or. &
      (ic == g_iccr ) .or. &
      (ic == g_iclf )      ) isWhiteSpace = PETSC_TRUE

end function IsWhiteSpace

! *************************************************************************** !

subroutine ReadEGridArrayI(a, keyword, ierr, input, option, qerr)
  !
  ! Reads an Eclgrid grid array
  ! If it is a depth array, then flip sign to Pflotran elevation convention
  ! If is is a permeability array, convert from mD to m2
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18
  !
  implicit none

  PetscInt, intent(inout) :: a(:)
  character(len = *) :: keyword
  PetscInt, intent(out) :: ierr
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr

  PetscInt , allocatable :: column_buffer(:)
  PetscReal, allocatable :: ar(:)
  PetscInt  :: ix, iy, iz, izpft, ig, igpft, i

  ! Do the actual read operation

  allocate(ar(g_nxyz))
  ar = 1.0

  call ReadEvalues(ar, g_nxyz, keyword, 'GRID' , ierr, input, option, qerr)

  do i = 1, g_nxyz
   a(i) = nint(ar(i))
  enddo

  ! Allocate a column buffer to reorder to bottom-up convention

  allocate(column_buffer(g_nz))

  ! Re-order the input values (Eclipse goes down the columns, PFT goes up)

  do ix = 1, g_nx
    do iy = 1, g_ny

      ! Load a column of Eclipse values

      do iz = 1, g_nz
        ig = GetNaturalIndex(ix, iy, iz)
        column_buffer(iz) = a(ig)
      enddo

      ! Place into Pflotran-converion array with appropriate data modifications

      do iz = 1, g_nz
        izpft = g_nz - iz + 1
        igpft = GetNaturalIndex(ix, iy, izpft)
        a(igpft) = column_buffer(iz)
      enddo
    enddo

  enddo

  ! Free the column buffer

  deallocate(column_buffer)
  deallocate(ar)

end subroutine ReadEGridArrayI

! *************************************************************************** !

subroutine ReadEGridArrayR(a, keyword, ierr, input, option, &
                           is_dep, is_perm, qerr)
  !
  ! Reads an Eclgrid grid array
  ! If it is a depth array, then flip sign to Pflotran elevation convention
  ! If is is a permeability array, convert from mD to m2
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  use Grid_Grdecl_Util_module, only : GetMDtoM2Conv

  implicit none

  PetscReal, intent(inout) :: a(:)
  character(len = *) :: keyword
  PetscInt, intent(out) :: ierr
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(in) :: is_dep
  PetscBool, intent(in) :: is_perm
  PetscBool, intent(inout) :: qerr

  PetscReal, allocatable :: column_buffer(:)
  PetscInt  :: ix, iy, iz, izpft, ig, igpft, nread
  PetscReal :: flip, conv, v

  qerr = PETSC_FALSE

  ! Check DIMENS read

  call checkDimensRead(qerr)

  if (.not.qerr) then

    ! Do the actual read operation

    nread = size(a)
    if (nread /= g_nxyz) then
      call SetError(keyword)
      qerr = PETSC_TRUE
    else
      call ReadEvalues(a, nread, keyword, 'GRID', ierr, input, option, qerr)
    endif
  endif

  if (.not.qerr) then

    ! If is a depth, change signs

    flip = 1.0
    if (is_dep) flip = z_flip

    ! If is perm, change units

    conv = 1.0
    if (is_perm) conv = GetMDtoM2Conv()

    ! Allocate a column buffer to reorder to bottom-up convention

    allocate(column_buffer(g_nz))

    ! Re-order the input values (Eclipse goes down the columns, PFT goes up)

    do ix = 1, g_nx
      do iy = 1, g_ny

        ! Load a column of Eclipse values

        do iz = 1, g_nz
          ig = GetNaturalIndex(ix, iy, iz)
          column_buffer(iz) = a(ig)
        enddo

        ! Place into Pflotran-conversion array with
        ! appropriate data modifications

        do iz = 1, g_nz
          izpft = g_nz - iz + 1
          igpft = GetNaturalIndex(ix, iy, izpft)
          v     = column_buffer(iz)
          if (Uninitialized(v)) then
            a(igpft) = v
          else
            a(igpft) = flip*conv*v
          endif
        enddo
      enddo

    enddo

  ! Free the column buffer

    deallocate(column_buffer)

  endif

end subroutine ReadEGridArrayR

! ************************************************************************** !

subroutine ReadECoordArray(a, ierr, input, option, qerr)
  !
  ! Reads an Eclgrid coord array
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscReal, intent(inout) :: a(:)
  PetscInt,  intent(out) :: ierr
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr

  PetscInt  :: i, nread

  nread = size(a)
  if (nread /= 6*g_nxpnyp) then
    call SetError('Coord read')
    qerr = PETSC_TRUE
  else
    call ReadEvalues(a, nread, 'COORD', 'GRID', ierr, input, option, qerr)
  endif

  do i = 3, nread, 3
    a(i) = z_flip*a(i)
  enddo

end subroutine ReadECoordArray

! *************************************************************************** !

subroutine ReadEZcornArray(a, ierr, input, option, qerr)
  !
  ! Reads an Eclgrid zcorn array
  ! If it is a depth array, then flip sign to Pflotran elevation convention
  ! If is is a permeability array, convert from mD to m2
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscReal, intent(inout) :: a(:)
  PetscInt , intent(out) :: ierr
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr

  PetscReal, allocatable :: column_buffer(:)
  PetscInt  :: inx , iny, inz, nnx, nny, nnz, nnxy, inode, nread, inzpft

  qerr = PETSC_TRUE

  ! Do the actual read operation

  nread = size(a)
  if (nread /= 8*g_nxyz) then
    call SetError('Zcorn read')
    qerr = PETSC_TRUE
  else
    call ReadEvalues(a, nread, 'ZCORN', 'GRID', ierr, input, option, qerr)
  endif

  if (.not.qerr) then

    nnx = 2*g_nx
    nny = 2*g_ny
    nnz = 2*g_nz
    nnxy = nnx*nny

    ! Allocate a column buffer to reorder to bottom-up convention

    allocate(column_buffer(nnz))

    ! Re-order the input values (Eclipse goes down the columns, PFT goes up)

    do inx = 1, nnx
      do iny = 1, nny

        ! Load a column of Eclipse values

        do inz = 1, nnz
          inode = nnxy*(inz-1) + nnx*(iny-1) + inx
          column_buffer(inz) = a(inode)
        enddo

        ! Unload the values in reverse z-order
        ! and change sign to Pflotran convention

        do inz = 1, nnz
          inzpft = nnz - inz + 1
          inode = nnxy*(inzpft-1) + nnx*(iny-1) + inx
          a(inode) = z_flip*column_buffer(inz)
        enddo

      enddo
    enddo

  ! Free the column buffer

    deallocate(column_buffer)

  endif

end subroutine ReadEZcornArray

! *************************************************************************** !

subroutine ReadEvalues(a, n, keyword, section, ierr, input, option, qerr)
  !
  ! Read a series of values from an Eclipse syntax file
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscReal, intent(inout) :: a(:)
  PetscInt , intent(in)    :: n
  character(len = *) :: keyword
  character(len = *) :: section
  PetscInt, intent(out) :: ierr
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr
  
  PetscInt  :: i, iostat, istar, nstack
  PetscReal :: dval
  PetscBool :: exittime
  character(len = MAXWORDLENGTH) :: word
  character(len = MAXWORDLENGTH) :: repc
  character(len = MAXWORDLENGTH) :: hold

  ierr   = 0
  qerr   = PETSC_FALSE

  i      = 0
  iostat = 0
  dval   = 0.0
  exittime = PETSC_FALSE
  nstack = 0

  call InputReadPflotranString(input, option)
  if (InputCheckExit(input, option)) exittime = PETSC_TRUE

  if (.not. exittime) then

    g_column = 1
    do i = 1, n

      if (nstack > 0) then

        ! Case of value in stack

        a(i)   = dval
        nstack = nstack - 1

      else

        ! Case of no value if stack

        call getNextWord(word, exitTime, input, option)
        if (exitTime) exit

        ! Look at word, test for repeat count

        istar = scan(word, '*')
        if (istar > 0) then
          repc = word(:istar-1)
          hold = word(istar+1:)
          read(repc, *, iostat = ierr) nstack
          word = hold
          nstack = nstack - 1
        endif

        ! Read the actual value

        read(word, *, iostat = ierr) dval
        if (ierr /= 0) then
          call InputErrorMsg(input, option, keyword, section)
          ierr = 1
          exit
        else
          a(i) = dval
        endif

      endif

    enddo

  endif

  if (ierr == 1) then
    qerr = PETSC_TRUE
  endif

end subroutine ReadEvalues

!*****************************************************************************!

subroutine CopyWellData(wto, wfrom)
  !
  ! Copy well data from one structure to another
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  type(well_locn_type), intent(out) :: wto
  type(well_locn_type), intent(in ) :: wfrom

  wto%ikl    = wfrom%ikl
  wto%iku    = wfrom%iku
  wto%w_name = wfrom%w_name

end subroutine CopyWellData

!*****************************************************************************!

subroutine CopyCmplData(clto, clfrom)
  !
  ! Copy completion data from one structure to another
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  type(cmpl_data_type), intent(out) :: clto
  type(cmpl_data_type), intent(in ) :: clfrom

  clto%ci = clfrom%ci
  clto%cj = clfrom%cj
  clto%ck = clfrom%ck
  clto%ik = clfrom%ik

  clto%iw = clfrom%iw
  clto%ia = clfrom%ia

  clto%dx = clfrom%dx
  clto%dy = clfrom%dy
  clto%dz = clfrom%dz

  clto%z  = clfrom%z
  
end subroutine CopyCmplData

!*****************************************************************************!

function FindWellIndex(name, iw)
  !
  ! Given a well name, find the well index
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscBool :: FindWellIndex
  character(len = *), intent(in ) :: name
  PetscInt          , intent(out) :: iw

  PetscInt :: jw

  FindWellIndex = PETSC_FALSE
  iw = -1

  do jw = 1 , g_nwell_data

    if (g_well_data(jw)%w_name == name) then
      FindWellIndex = PETSC_TRUE
      iw = jw
      exit
    endif
  enddo

end function findWellIndex

!*****************************************************************************!

function GetGrdNCmpl(iw)
  !
  ! Given a well name, find the number of completions
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscInt             :: GetGrdNCmpl
  PetscInt, intent(in) :: iw

 if (iw <= g_nwell_data) then
   GetGrdNCmpl = g_well_data(iw)%iku - g_well_data(iw)%ikl + 1
 else
   GetGrdNCmpl = 0
 endif

end function GetGrdNCmpl

!*****************************************************************************!

subroutine GetCmplData(iw, ik, ci, cj, ck, ia, dx, dy, dz, z)
  !
  ! Given a well index and completion index, return completion data
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscInt, intent(in)   :: iw, ik
  PetscInt, intent(out)  :: ci, cj, ck, ia
  PetscReal, intent(out) :: dx, dy, dz, z

  PetscInt :: ikl, k, ig

 if (iw <= g_nwell_data) then
   ikl = g_well_data(iw)%ikl
 else
   ikl = 0
 endif

 k = ikl + ik - 1

 ci = g_cmpl_data(k)%ci
 cj = g_cmpl_data(k)%cj
 ck = g_cmpl_data(k)%ck

 ig = g_cmpl_data(k)%ig
 ia = g_cmpl_data(k)%ia

 dx = g_cmpl_data(k)%dx
 dy = g_cmpl_data(k)%dy
 dz = g_cmpl_data(k)%dz

 z  = g_cmpl_data(k)%z

end subroutine GetCmplData

subroutine IsCPG()

  implicit none

  g_iscpg = PETSC_TRUE

  if (.not.g_cpgallocated) then
    allocate(g_coord(6*g_nxpnyp));g_coord =  0.0
    allocate(g_zcorn(8*g_nxyz  ));g_zcorn =  0.0
    g_cpgallocated = PETSC_TRUE
  endif

end subroutine isCPG

! ************************************************************************** !

subroutine PermPoroExchangeAndSet(poro_p, permx_p, permy_p, permz_p, &
                                  inatsend, nlmax, option)
  !
  ! Perm and porosity values from the grdecl file are known on the IO rank
  ! These are needed on the other ranks, but only for the cells on those ranks
  ! To avoid over-sending, the other ranks send lists to the IO ranks,
  ! and the IO ranks send back the required values.
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscReal, pointer :: poro_p  (:)
  PetscReal, pointer :: permx_p (:)
  PetscReal, pointer :: permy_p (:)
  PetscReal, pointer :: permz_p (:)
  PetscInt , pointer :: inatsend(:)
  PetscInt, intent(in) :: nlmax
  type(option_type) :: option

  PetscInt :: irank, iorank, t_rank, ibp
  PetscInt :: ilt, ilo, ino
  PetscInt :: nlo, ierr
  PetscMPIInt :: int_mpi, temp_int_array(1), status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt, parameter :: tag_mpi = 0
  PetscReal :: poro, permx, permy, permz

  ! Work arrays

  PetscInt , allocatable :: winat(:)
  PetscReal, allocatable :: wporo(:)
  PetscReal, allocatable :: wperm(:)

  ! Set scalars

  ierr   = 0
  iorank = option%io_rank
  t_rank = option%myrank

  ! Loop over exchange operations between IO ranks and non-IO ranks (irank)

  do irank = 0, option%mycommsize-1

    if (irank /= iorank) then

      ! Consider exchange between irank and iorank

      if (t_rank  == irank) then

        ! This is irank: send nlmax value and then irequest array to ioproc

        temp_int_array(1) = nlmax
        call MPI_Send(temp_int_array, ONE_INTEGER_MPI, MPI_INTEGER, &
                      iorank, tag_mpi, option%mycomm, ierr)
        int_mpi = nlmax
        call MPI_Send(inatsend      , int_mpi        , MPI_INTEGER, &
                      iorank, tag_mpi, option%mycomm, ierr)

        ! Allocate buffers to hold response poro and perm values from ioproc

        allocate(wporo(  nlmax))
        allocate(wperm(3*nlmax))

        ! Receive poro and perm values from ioproc

        int_mpi =   nlmax
        call MPI_Recv(wporo, int_mpi, MPI_DOUBLE_PRECISION, iorank, &
                      MPI_ANY_TAG, option%mycomm, status_mpi, ierr)
        int_mpi = 3*nlmax
        call MPI_Recv(wperm, int_mpi, MPI_DOUBLE_PRECISION, iorank, &
                      MPI_ANY_TAG, option%mycomm, status_mpi, ierr)

        ! Copy buffers into correct storage locations

        do ilt = 1, nlmax
          ibp = 3*(ilt-1)
          poro_p (ilt) = wporo(ilt  )
          permx_p(ilt) = wperm(ibp+1)
          permy_p(ilt) = wperm(ibp+2)
          permz_p(ilt) = wperm(ibp+3)
        enddo

        ! Free buffers

        deallocate(wporo)
        deallocate(wperm)

      endif

      if (t_rank  == iorank) then

       !This is IO rank - receive the other rank nlmax value (nlo) from irank

        temp_int_array(1) = 0
        call MPI_Recv(temp_int_array, ONE_INTEGER_MPI, MPI_INTEGER, &
                      irank, MPI_ANY_TAG, option%mycomm, status_mpi, ierr)
        nlo    = temp_int_array(1)

        ! Allocate work array to hold the received natural addresses

        allocate(winat(  nlo))
        allocate(wporo(  nlo))
        allocate(wperm(3*nlo))

        ! Receive the natural address array fom irank

        int_mpi = nlo
        call MPI_Recv(winat, int_mpi, MPI_INTEGER, irank, MPI_ANY_TAG, &
                      option%mycomm, status_mpi, ierr)

        ! On this (io) proc, set up the perm and poro values to be returned

        do ilo = 1, nlo
          ino = winat(ilo)
          call GetPoroPermValues(ino, poro, permx, permy, permz)
          ibp = 3*(ilo-1)
          wporo(ilo  ) = poro
          wperm(ibp+1) = permx
          wperm(ibp+2) = permy
          wperm(ibp+3) = permz
        enddo

        ! Send back the poro and perm values to irank

        int_mpi = nlo
        call MPI_Send(wporo, int_mpi, MPI_DOUBLE_PRECISION, &
                      irank, tag_mpi, option%mycomm, ierr)
        int_mpi = 3*nlo
        call MPI_Send(wperm, int_mpi, MPI_DOUBLE_PRECISION, &
                      irank, tag_mpi, option%mycomm, ierr)

        ! Free the work arrays

        deallocate(winat)
        deallocate(wporo)
        deallocate(wperm)

      endif
    endif

    ! Barrier call to stop other proces running ahead

    call MPI_Barrier(option%mycomm, ierr)
  enddo

end subroutine PermPoroExchangeAndSet

end module Grid_Grdecl_module
