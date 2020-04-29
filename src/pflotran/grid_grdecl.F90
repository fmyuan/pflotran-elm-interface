module Grid_Grdecl_module

!  Module to read a grid using industry-standard Eclipse keyword syntax
!  and convert it to a Pflotran explicit unstructured grid

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Grid_Unstructured_Aux_module
  use Option_module
  use Input_Aux_module
  use Grid_Grdecl_Util_module

  implicit none

  private

  ! ASCII codes for tab, line feed, carriage return and quotes

  PetscInt, parameter :: g_ictab  = 9
  PetscInt, parameter :: g_iclf   = 10
  PetscInt, parameter :: g_iccr   = 13
  PetscInt, parameter :: g_squote = 39
  PetscInt, parameter :: g_dquote = 34

  ! Cartesian frame pointers

  PetscInt, parameter :: g_xdir = 1
  PetscInt, parameter :: g_ydir = 2
  PetscInt, parameter :: g_zdir = 3

  PetscInt, parameter :: g_ndir = 3

  ! Face types for (FAULTS/MULTFLT)

  PetscInt, parameter :: g_facei  = 1
  PetscInt, parameter :: g_faceim = 2
  PetscInt, parameter :: g_facej  = 3
  PetscInt, parameter :: g_facejm = 4
  PetscInt, parameter :: g_facek  = 5
  PetscInt, parameter :: g_facekm = 6

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

  ! Flag indicating SATNUM satn. table indices read and max. SATNUM value

  PetscInt :: g_rsatn  = 0
  PetscInt :: g_maxsatn= 0

  ! Flag indicating DIMENS read (contains problem dimensions)

  PetscBool :: g_dimens_read = PETSC_FALSE

  ! z-flip from Eclipse to Pflotran convention and corner-point flags

  PetscBool :: g_iscpg = PETSC_FALSE
  PetscReal :: z_flip  = -1.0

  PetscBool :: g_cpgallocated = PETSC_FALSE
  PetscBool :: g_isnewtran    = PETSC_FALSE

  ! Counters for active cells, connections, max connections, faults, pinchouts

  PetscInt :: g_na  = 1
  PetscInt :: g_nc  = 0
  PetscInt :: g_mc  = 1
  PetscInt :: g_nf  = 0
  PetscInt :: g_npo = 0

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
  PetscInt , pointer :: g_satn(:) => null()

  PetscInt , pointer :: g_gtoa(:) => null()

  ! Active sized arrays

  PetscInt , pointer :: g_atog(:) => null()
  PetscInt , pointer :: g_atoc(:) => null()

  PetscReal, pointer :: g_bvol(:) => null()

  PetscReal, pointer :: g_x   (:) => null()
  PetscReal, pointer :: g_y   (:) => null()
  PetscReal, pointer :: g_z   (:) => null()

  ! Connection sized arrays

  PetscInt , pointer :: g_cia (:) => null()
  PetscInt , pointer :: g_cja (:) => null()

  PetscReal, pointer :: g_ccx  (:) => null()
  PetscReal, pointer :: g_ccy  (:) => null()
  PetscReal, pointer :: g_ccz  (:) => null()
  PetscReal, pointer :: g_carea(:) => null()

  PetscReal :: g_mapaxes(6) = 0.0
  PetscReal :: g_minpv(1)   = 1.0D-3
  PetscReal :: g_pinch(1)   = 1.0D-3

  ! Column pointer for reader

  PetscInt :: g_column = 1

  character(len = MAXSTRINGLENGTH) :: g_error_string = 'OK'
  PetscInt :: g_error_flag = 0

!  Null test flag

  PetscBool :: g_geometry_test = PETSC_FALSE

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
  public  :: GetSatnumSet, GetSatnumValue, SatnumExchangeAndSet

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

  ! Type containing fault data

  type, public :: fault_data_type

    private

    character(len=8) :: zname

    PetscInt  :: il    = 1
    PetscInt  :: iu    = 0
    PetscInt  :: jl    = 1
    PetscInt  :: ju    = 0
    PetscInt  :: kl    = 1
    PetscInt  :: ku    = 1
    PetscInt  :: iface = 0

  end type fault_data_type

  ! Counters for fault instances

  PetscInt :: g_nfault_data = 0
  PetscInt :: g_mfault_data = 0

  ! Array of fault instances

  type(fault_data_type), allocatable :: g_fault_data(:)

  ! Type containing multflt data

  type, public :: multflt_data_type

    private

    character(len=8) :: zname

    PetscReal  :: vm

  end type multflt_data_type

  ! Counters for multflt instances

  PetscInt :: g_nmultflt_data = 0
  PetscInt :: g_mmultflt_data = 0

  ! Array of multflt_data instances

  type(multflt_data_type), allocatable :: g_multflt_data(:)

contains

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

function GetSatnumSet(maxsatn)
  !
  ! Get flag indicating that SATNUM array is being used
  !
  ! Author: Dave Ponting
  ! Date: 02/14/18
  !

  implicit none

  PetscBool :: GetSatnumSet
  PetscInt, intent(out) :: maxsatn

  GetSatnumSet = PETSC_FALSE
  if (g_rsatn == 1) GetSatnumSet = PETSC_TRUE
  maxsatn = g_maxsatn

end function GetSatnumSet

! *************************************************************************** !

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

  !  Set non-blocking error detection only (reader will not exit)

  call OptionSetBlocking(option, PETSC_FALSE)
  option%error_while_nonblocking = PETSC_FALSE

  !  Now read grdecl file on the io-rank only

  if (option%myrank == option%io_rank) then
    call GrdeclReader(input, option)
  endif
  call MPI_Bcast(g_error_flag, ONE_INTEGER_MPI, MPI_INTEGER, &
                 option%io_rank, option%mycomm, ierr)
  if (g_error_flag>0) then
    input%ierr = 1
    call MPI_Bcast(g_error_string, MAXSTRINGLENGTH, MPI_CHARACTER, &
                   option%io_rank, option%mycomm, ierr)
    call InputErrorMsg(input, option, 'GRDECL file', g_error_string)
  endif

  ! Reset default blocking error state

  call OptionSetBlocking(option, PETSC_TRUE)
  call OptionCheckNonBlockingError(option)

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

  ! Distribute and store cell counts and satnum flag

  call DistributeCounts(option)

  ! Set up the cell locations

  if (option%myrank == option%io_rank) then
    call CreateElements(unstructured_grid, explicit_grid)
  endif

end subroutine UGrdEclExplicitRead

! *************************************************************************** !

subroutine WriteStaticDataAndCleanup(write_ecl, eclipse_options, option)
  !
  ! If Eclipse files are required, output grid and init files
  ! Then clean up all the static data arrays which are no longer needed
  !
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

  PetscInt  :: iw, nw, nctw, mcpw

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
      efilename = trim(option%global_prefix)//trim(option%group_prefix)
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

! *************************************************************************** !

subroutine SetUGrdEclCmplLocation(wname, ci, cj, ckuser, cijk_d, qerr)
  !
  ! Set up a completion for well wname at location (ci, cj, ckuser)
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

  ! Set up Pflotran ck value.
  ! If entered using cijk_d, is count from top, set negative to flag this later

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
      call setError('Same WELL_DATA well has appeared more than once', qerr)
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

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

subroutine ExtendFaultData()
  !
  ! Extend the g_fault_data structure
  !
  ! Author: Dave Ponting
  ! Date: 02/11/18
  !

  implicit none

  type(fault_data_type), allocatable :: temp_fault_data(:)

  PetscInt :: mfault_data_new, ifault

  ! If existing fault locations, copy them

  if (g_nfault_data > 0) then

    allocate(temp_fault_data(g_nfault_data))

    do ifault = 1, g_nfault_data
      call CopyFaultData(temp_fault_data(ifault), g_fault_data(ifault))
    enddo

  endif

  ! Deallocate, increase and re-allocate the fault store

  if (allocated(g_fault_data)) deallocate(g_fault_data)

  if (g_mfault_data == 0 ) then
    mfault_data_new = 10
  else
    mfault_data_new = 2 * g_mfault_data
  endif

  allocate(g_fault_data(mfault_data_new))

  ! Copy over the existing fault locations

  do ifault = 1, g_nfault_data
    call CopyFaultData(g_fault_data(ifault), temp_fault_data(ifault))
  enddo

  ! Deallocate the temporary fault locations

  if (allocated(temp_fault_data)) deallocate(temp_fault_data)

  ! Store the new dimension

  g_mfault_data = mfault_data_new

end subroutine ExtendFaultData

! *************************************************************************** !

subroutine ExtendMultfltData()
  !
  ! Extend the g_multflt_data structure
  !
  ! Author: Dave Ponting
  ! Date: 02/11/18
  !

  implicit none

  type(multflt_data_type), allocatable :: temp_multflt_data(:)

  PetscInt :: mmultflt_data_new, imultflt

  ! If existing multflt values, copy them

  if (g_nmultflt_data > 0) then

    allocate(temp_multflt_data(g_nmultflt_data))

    do imultflt = 1, g_nmultflt_data
      call CopyMultfltData(temp_multflt_data(imultflt), &
                              g_multflt_data(imultflt))
    enddo

  endif

  ! Deallocate, increase and re-allocate the multflt store

  if (allocated(g_multflt_data)) deallocate(g_multflt_data)

  if (g_mmultflt_data == 0) then
    mmultflt_data_new = 10
  else
    mmultflt_data_new = 2*g_nmultflt_data
  endif

  allocate(g_multflt_data(mmultflt_data_new))

  ! Copy over the existing multflt store

  do imultflt = 1, g_nmultflt_data
    call CopyMultfltData(g_multflt_data(imultflt), temp_multflt_data(imultflt))
  enddo

  ! Deallocate the temporary multflt store

  if (allocated(temp_multflt_data)) deallocate(temp_multflt_data)

  ! Store the new dimension

  g_mmultflt_data = mmultflt_data_new

end subroutine ExtendMultfltData

! *************************************************************************** !

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

! *************************************************************************** !

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

  PetscBool, parameter :: dep_no  = PETSC_FALSE
  PetscBool, parameter :: dep_yes = PETSC_TRUE
  PetscBool, parameter :: prm_no  = PETSC_FALSE
  PetscBool, parameter :: prm_yes = PETSC_TRUE
  PetscBool, parameter :: nn_no   = PETSC_FALSE
  PetscBool, parameter :: nn_yes  = PETSC_TRUE

  PetscInt, parameter :: mzbuf = 10

  PetscInt  :: ierr, nread, nxsg, nysg, nzsg
  PetscBool :: qerr, is_metric, newtran_force, oldtran_force
  PetscReal :: dimens(3)
  character(len = MAXWORDLENGTH) :: word
  character(len = MAXWORDLENGTH) :: zbuf(mzbuf)
  character(len = MAXSTRINGLENGTH) :: zmess

  ierr = 0
  qerr = PETSC_FALSE
  zbuf = ' '
  dimens = 1.0

  newtran_force = PETSC_FALSE
  oldtran_force = PETSC_FALSE

  ! Read through items in the grdecl file
  call InputPushBlock(input,option)
  do

    call InputReadPflotranStringNotComment(input, option)
    if (option%error_while_nonblocking) then
      ! I/O error detected so set qerr and exit
      call SetError('Unable to read from GRDECL file', qerr)
      exit
    endif
    if (input%ierr /= 0) exit                   ! Detect eof and leave

  ! Keyword found

    call InputReadCard(input, option, word)
    call CheckError(input, 'ECLGRD subkeyword', qerr);if ( qerr ) exit
    call StringToUpper(word)

    if (word(1:1) /= '/') print *, 'pflotran card:: ', word

    select case(trim(word))
      case('DIMENS')
        call ReadEvalues(dimens, 3, 'DIMENS', 'GRID', &
                         input, option, nn_yes, qerr)
        if (.not. qerr) call SetDimens(dimens)
      case('COORD')
        call IsCPG()
        g_isnewtran = PETSC_TRUE
        call checkDimensRead(qerr)
        if (.not. qerr) call ReadECoordArray(g_coord, input, option, qerr)
      case('ZCORN')
        call IsCPG()
        g_isnewtran = PETSC_TRUE
        call checkDimensRead(qerr)
        if (.not. qerr) call ReadEZcornArray(g_zcorn, input, option, qerr)
      case('DX')
        call ReadEGridArrR(g_dx, 'DX'   , &
                           input, option, dep_no , prm_no , nn_yes, qerr)
      case('DY')
        call ReadEGridArrR(g_dy, 'DY'   , &
                           input, option, dep_no , prm_no , nn_yes, qerr)
      case('DZ')
        call ReadEGridArrR(g_dz, 'DZ'   , &
                           input, option, dep_no , prm_no , nn_yes, qerr)
      case('PERMX')
        call ReadEGridArrR(g_kx, 'PERMX', &
                           input, option, dep_no, prm_yes, nn_yes, qerr)
      case('PERMY')
        call ReadEGridArrR(g_ky, 'PERMY', &
                           input, option, dep_no, prm_yes, nn_yes, qerr)
      case('PERMZ')
        call ReadEGridArrR(g_kz, 'PERMZ', &
                           input, option, dep_no, prm_yes, nn_yes, qerr)
      case('MULTX')
        call ReadEGridArrR(g_mx, 'MULTX', &
                           input, option, dep_no, prm_no , nn_yes, qerr)
      case('MULTY')
        call ReadEGridArrR(g_my, 'MULTY', &
                           input, option, dep_no, prm_no , nn_yes, qerr)
      case('MULTZ')
        call ReadEGridArrR(g_mz, 'MULTZ', &
                           input, option, dep_no, prm_no , nn_yes, qerr)
      case('PORO')
        call ReadEGridArrR(g_poro, 'PORO', &
                           input, option, dep_no, prm_no , nn_yes, qerr)
      case('TOPS')
        call ReadEGridArrR(g_tops, 'TOPS', &
                           input, option, dep_yes, prm_no , nn_no , qerr)
      case('NTG')
        call ReadEGridArrR(g_ntg, 'NTG', &
                           input, option, dep_no , prm_no , nn_yes, qerr)
      case('ACTNUM')
        call ReadEGridArrI(g_actn, 'ACTNUM', &
                           input, option, nn_yes, qerr)
      case('SATNUM')
        g_rsatn = 1
        call ReadEGridArrI(g_satn, 'SATNUM', &
                           input, option, nn_yes, qerr)
        g_maxsatn= maxval(g_satn)
      case('MINPV')
        call ReadEvalues(g_minpv, 1, 'MINPV', 'GRID' , &
                         input, option, nn_yes, qerr)
      case('PINCH')
        call ReadEvalues(g_pinch, 1, 'PINCH', 'GRID' , &
                         input, option, nn_yes, qerr)
      case('NEWTRAN')
        g_isnewtran   = PETSC_TRUE
        newtran_force = PETSC_TRUE
      case('OLDTRAN')
        g_isnewtran   = PETSC_FALSE
        oldtran_force = PETSC_TRUE
      case('ADD', 'COPY', 'EQUALS', 'MULTIPLY')
       call HandleOpKeyword(word, input, option, qerr)
      case('FAULTS')
        call ReadFaults(word, input, option, qerr)
      case('MULTFLT')
        call ReadMultflt(word, input, option, qerr)
      case('ECHO', 'NOECHO', 'INIT', 'EDIT', 'REGIONS')  ! No need to action
      case('MAPUNITS')
        call ReadEstrings(word, zbuf, 1, input, option, qerr)
        is_metric = StringCompareIgnoreCase(zbuf(1), 'metres')
        if (.not.is_metric) &
          call SetError('MAXUNITS units must be metres', qerr)
      case('GRIDUNIT')
        call ReadEstrings(word, zbuf, 2, input, option, qerr)
        is_metric = StringCompareIgnoreCase(zbuf(1), 'metres')
        if (.not.is_metric) &
          call SetError('GRIDUNIT units must be metres', qerr)
      case('GDORIENT')
        call ReadEstrings(word, zbuf, 5, input, option, qerr)
      case('GRIDFILE')
        call ReadEstrings(word, zbuf, 2, input, option, qerr)
      case('SPECGRID')
        call ReadEstrings(word, zbuf, 5, input, option, qerr)
         word = zbuf(1)
         read(word, *, iostat=ierr) nxsg
         if (nxsg /= g_nx .or. ierr /= 0) then
           call SetError('SPECGRID NX ' &
                         // trim(word) // ' does not match DIMENS NX', qerr)
         endif
         word = zbuf(2)
         read(word, *, iostat=ierr) nysg
         if (nysg /= g_ny .or. ierr /= 0) then
           call SetError('SPECGRID NY ' &
                         // trim(word) // ' does not match DIMENS NY', qerr)
         endif
         word = zbuf(3)
         read(word, *, iostat=ierr) nzsg
         if (nzsg /= g_nz .or. ierr /= 0) then
           call SetError('SPECGRID NZ ' &
                         // trim(word) // ' does not match DIMENS NZ', qerr)
         endif
      case('MAPAXES')
        nread = size(g_mapaxes)
        call ReadEvalues(g_mapaxes, nread, 'MAPAXES', 'GRID', &
                         input, option, nn_no, qerr)
      case('/')  ! Isolated un-used terminator on its own line is harmless
      case default
       zmess = 'GRDECL sub-keyword ' // trim(word) // ' not recognised'
       call SetError(zmess, qerr)
    end select

   !  Leave if error flag set

    if (qerr) exit
  enddo
  call InputPopBlock(input,option)

  !  Case in which NEWTRAN or OLDTRAN specified explicitly

  if (newtran_force .and. oldtran_force) then
    zmess = 'Both NEWTRAN and OLDTRAN specified'
    call SetError(zmess, qerr)
  elseif (newtran_force) then
    g_isnewtran = PETSC_TRUE
  elseif (oldtran_force) then
    g_isnewtran = PETSC_FALSE
  endif

  !  Check for errors and process data

  if (.not.qerr) then

    ! Process the grid data read

    call ProcessGridData()
    if (g_na == 0) then

!  Case of no active cells: error

      zmess = 'Problem has no active cells'
      call SetError(zmess, qerr)

    else

    ! Process the well data read

      call ProcessWellData(qerr)

    endif

  endif

end subroutine GrdeclReader

! *************************************************************************** !

subroutine ReadFaults(zkey, input, option, qerr)
  !
  ! Handle a FAULTS keyword
  !
  ! Author: Dave Ponting
  ! Date: 02/10/19

  use String_module

  implicit none

  character(len = *), intent(in) :: zkey

  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr

  PetscInt, parameter :: ma = 8
  character(len = MAXWORDLENGTH) :: za(ma),zface

  PetscInt  :: ixl, ixu, iyl, iyu, izl, izu , izle, izue

!  Initialise values (note read in Eclipse iz convention)

  qerr = PETSC_FALSE

  ixl = 1
  ixu = g_nx

  iyl = 1
  iyu = g_ny

  izle = 1
  izue = g_nz

!  Loop through data until a null record found

  do

    za    = ' '
    za(1) = '/'
    call ReadEstrings(zkey, za, ma, input, option, qerr)
    if (qerr) exit

    if (StringCompareIgnoreCase(za(1), '/')) exit

   ! Read and error-check data

    call ProcessArgToInt(ixl , za(2), zkey, 1   , g_nx, qerr); if (qerr) exit
    call ProcessArgToInt(ixu , za(3), zkey, ixl , g_nx, qerr); if (qerr) exit

    call ProcessArgToInt(iyl , za(4), zkey, 1   , g_ny, qerr); if (qerr) exit
    call ProcessArgToInt(iyu , za(5), zkey, iyl , g_ny, qerr); if (qerr) exit

    call ProcessArgToInt(izle, za(6), zkey, 1   , g_nz, qerr); if (qerr) exit
    call ProcessArgToInt(izue, za(7), zkey, izle, g_nz, qerr); if (qerr) exit

    zface = za(8)

!   Convert to Pflotran layer convention
!   Note lower E-conv. layer becomes the upper P-conv. layer and v.v.

    izl = g_nz-izue+1
    izu = g_nz-izle+1

    call StoreFault(za(1), ixl, ixu, iyl, iyu, izl, izu, zface, qerr)
    if (qerr) exit

  enddo

end subroutine ReadFaults

! *************************************************************************** !

subroutine ReadMultflt(zkey, input, option, qerr)
  !
  ! Handle a MULTFLT keyword
  !
  ! Author: Dave Ponting
  ! Date: 02/10/19

  use String_module

  implicit none

  character(len = *), intent(in) :: zkey
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr

  PetscInt, parameter :: ma = 2
  character(len = MAXWORDLENGTH) :: za(ma)
  PetscReal :: vm

!  Initialise values

  qerr = PETSC_FALSE
  vm   = 1.0

!  Loop through data until a null record found

  do

    za    = ' '
    za(1) = '/'
    call ReadEstrings(zkey, za, ma, input, option, qerr);if (qerr) exit

    if (StringCompareIgnoreCase(za(1), '/')) exit

    call ProcessArgToReal(vm, za(2), zkey, qerr); if (qerr) exit

    call StoreMultflt(za(1), vm)

  enddo

end subroutine ReadMultflt

! *************************************************************************** !

subroutine StoreFault( zname, il, iu, jl, ju, kl, ku, zface, qerr )

  !
  ! Store a fault structure, expanding the array as required
  !
  ! Author: Dave Ponting
  ! Date: 02/10/19

  use String_module

  implicit none

  character(len=*), intent(in) :: zname, zface
  PetscInt, intent(in) :: il, iu, jl, ju, kl, ku
  PetscBool, intent(out) :: qerr

  PetscInt :: ifault, iface
  character(len = MAXSTRINGLENGTH) :: zmess

  qerr = PETSC_FALSE
  iface = 0

  ifault = g_nfault_data + 1

  if (ifault > g_mfault_data) then
    call ExtendFaultData()
  endif

  ! Positive or negative X-face (I accepted as alias)

  if (     StringCompareIgnoreCase(zface, 'I' ) &
      .or. StringCompareIgnoreCase(zface, 'X' )) iface = g_facei

  if (     StringCompareIgnoreCase(zface, 'I-') &
      .or. StringCompareIgnoreCase(zface, 'X-')) iface = g_faceim

  ! Positive or negative Y-face (J accepted as alias)

  if (     StringCompareIgnoreCase(zface, 'J' ) &
      .or. StringCompareIgnoreCase(zface, 'Y' )) iface = g_facej

  if (     StringCompareIgnoreCase(zface, 'J-') &
      .or. StringCompareIgnoreCase(zface, 'Y-')) iface = g_facejm

  ! Positive or negative Z-face (K accepted as alias)

  if (     StringCompareIgnoreCase(zface, 'K' ) &
      .or. StringCompareIgnoreCase(zface, 'Z' )) iface = g_facek

  if (     StringCompareIgnoreCase(zface, 'K-') &
      .or. StringCompareIgnoreCase(zface, 'Z-')) iface = g_facekm

  ! Error trap on face string (a harmless null fault will be stored)

  if (iface == 0 ) then
    zmess = 'Face '//trim(zface)//' not recognised'
    call SetError(zmess, qerr)
  endif

  ! Store values

  g_fault_data(ifault)%zname = zname
  g_fault_data(ifault)%il    = il
  g_fault_data(ifault)%iu    = iu
  g_fault_data(ifault)%jl    = jl
  g_fault_data(ifault)%ju    = ju
  g_fault_data(ifault)%kl    = kl
  g_fault_data(ifault)%ku    = ku
  g_fault_data(ifault)%iface = iface

  g_nfault_data = ifault

end subroutine StoreFault

! *************************************************************************** !

subroutine StoreMultflt( zf, vm )

  !
  ! Store a multflt structure, expanding the array as required
  !
  ! Author: Dave Ponting
  ! Date: 02/10/19

  implicit none

  character(len=*), intent(in) :: zf
  PetscReal, intent(in) :: vm

  PetscInt :: imultflt

  imultflt = g_nmultflt_data + 1

  if (imultflt > g_mmultflt_data) then
    call ExtendMultfltData()
  endif

  g_multflt_data(imultflt)%zname = zf
  g_multflt_data(imultflt)%vm    = vm

  g_nmultflt_data = imultflt

end subroutine StoreMultflt

! *************************************************************************** !

subroutine HandleOpKeyword(zkey, input, option, qerr)
  !
  ! Handle ADD, COPY, EQUALS or MULTIPLY
  !
  ! Author: Dave Ponting
  ! Date: 02/10/19

  use String_module

  implicit none

  character(len=*), intent(in) :: zkey
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr

  PetscInt, parameter :: ma = 8
  character(len = MAXWORDLENGTH) :: za(ma)

  PetscBool :: isadd, iscpy, iseql, ismlt
  PetscInt  :: ixl, ixu, iyl, iyu, izl, izu, izle, izue

!  Initialise values

  qerr = PETSC_FALSE

  isadd = StringCompareIgnoreCase(zkey, 'add'     )
  iscpy = StringCompareIgnoreCase(zkey, 'copy'    )
  iseql = StringCompareIgnoreCase(zkey, 'equals'  )
  ismlt = StringCompareIgnoreCase(zkey, 'multiply')

!  Initialise values (note read in Eclipse iz convention)

  ixl  = 1
  ixu  = g_nx

  iyl  = 1
  iyu  = g_ny

  izle = 1
  izue = g_nz

  do

    za    = ' '
    za(1) = '/'
    call ReadEstrings(zkey, za, ma, input, option, qerr)
    if (qerr) exit

    if (StringCompareIgnoreCase(za(1), '/')) exit

    call ProcessArgToInt(ixl , za(3), zkey, 1   , g_nx, qerr); if (qerr) exit
    call ProcessArgToInt(ixu , za(4), zkey, ixl , g_nx, qerr); if (qerr) exit

    call ProcessArgToInt(iyl , za(5), zkey, 1   , g_ny, qerr); if (qerr) exit
    call ProcessArgToInt(iyu , za(6), zkey, iyl , g_ny, qerr); if (qerr) exit

    call ProcessArgToInt(izle, za(7), zkey, 1   , g_nz, qerr); if (qerr) exit
    call ProcessArgToInt(izue, za(8), zkey, izle, g_nz, qerr); if (qerr) exit

!   Convert to Pflotran layer convention
!   Note lower E-conv. layer becomes the upper P-conv. layer and v.v.

    izl = g_nz-izue+1
    izu = g_nz-izle+1

    call BoxOp( za(1), za(2), ixl, ixu, iyl, iyu, izl, izu, &
                isadd, iscpy, iseql, ismlt, zkey, qerr )
    if (qerr) exit

  enddo

end subroutine HandleOpKeyword

! *************************************************************************** !

subroutine BoxOp( zto, zfr, ixl, ixu, &
                            iyl, iyu, &
                            izl, izu, isadd, iscpy, iseql, ismlt, zk, qerr )

  !
  ! Do a add/copy/equals/multiply operation
  !
  ! Author: Dave Ponting
  ! Date: 02/11/19
  !

  implicit none

  character(len = *), intent(in) :: zto
  character(len = *), intent(in) :: zfr
  PetscInt , intent(in) :: ixl, ixu, iyl, iyu, izl, izu
  PetscBool, intent(in) :: isadd, iscpy, iseql, ismlt
  character(len = *), intent(in) :: zk
  PetscBool, intent(out) :: qerr

  PetscBool :: isint, isreal, isintb, isrealb, &
               is_mult_only , is_mult_dummy, &
               is_depth     , is_depth_dummy, &
               is_perm      , is_perm_dummy

  PetscInt , pointer :: ai(:), bi(:)
  PetscReal, pointer :: ar(:), br(:)

  PetscReal :: v,conv
  PetscInt  :: ix, iy, iz, ig, iv

!  Initialise values

  isint   = PETSC_FALSE
  isreal  = PETSC_FALSE
  isintb  = PETSC_FALSE
  isrealb = PETSC_FALSE

  is_mult_only  = PETSC_FALSE
  is_mult_dummy = PETSC_FALSE

  is_depth       = PETSC_FALSE
  is_depth_dummy = PETSC_FALSE

  is_perm        = PETSC_FALSE
  is_perm_dummy  = PETSC_FALSE

  nullify(ai, ar, bi, br)

  v  = 0.0
  iv = 0

  qerr = PETSC_FALSE

!  Get the required arrays

  call GetGridArrayPointer( zto, isint, isreal, ai, ar, qerr, &
                            is_mult_only, is_depth, is_perm )
  if (iscpy) then
    call GetGridArrayPointer( zfr, isintb, isrealb, bi, br, qerr, &
                              is_mult_dummy, is_depth_dummy, is_perm_dummy)
  else
    call ProcessArgToReal(v, zfr, zk, qerr)
  endif

!  Check if operation is OK

  if (is_mult_only .and. (.not. ismlt)) then
    call SetError('Only MULT supported for ' // trim(zto), qerr)
  else

!  Carry out operations

    if (.not.qerr) then

!     Case of depth array like TOPS: convert set-or-add value to an elevation

      if (is_depth) then
        if (isadd .or. iseql) v = -v
      endif

!     Case of perm array like PERMX: convert set-or-add value from mD to m^2

      if (is_perm) then
        conv = GetMDtoM2Conv()
        if (isadd .or. iseql) v = v*conv
      endif

      do iz = izl, izu
        do iy = iyl, iyu
          do ix = ixl, ixu

            ig = GetNaturalIndex(ix, iy, iz)

!  Case of integer array
            if (isint) then
              iv = nint(v)
              if (isadd) ai(ig) = ai(ig) + iv
              if (iscpy) then
                if (isintb ) bi(ig) =      ai(ig)
                if (isrealb) bi(ig) = nint(ar(ig))
              endif
              if (iseql) ai(ig) = iv
              if (ismlt) ai(ig) = ai(ig) * iv
            endif

! Case of real array
            if (isreal) then
              if (isadd) ar(ig) = ar(ig) + v
              if (iscpy) then
                if (isintb ) br(ig) = real(ai(ig))
                if (isrealb) br(ig) =      ar(ig)
              endif
              if (iseql) ar(ig) = v
              if (ismlt) ar(ig) = ar(ig) * v
            endif

          enddo
        enddo
      enddo
    endif
  endif

end subroutine BoxOp

! *************************************************************************** !

subroutine GetGridArrayPointer(za, isint, isreal, ai, ar, qerr, &
                               is_mult_only, is_depth, is_perm)
  !
  ! Get a pointer to the grid array with mnemonic za
  !
  ! Author: Dave Ponting
  ! Date: 02/11/19
  !

  use String_module

  implicit none

  character(len = *), intent(in) :: za

  PetscBool, intent(out) :: isint
  PetscBool, intent(out) :: isreal

  PetscInt , pointer, intent(out) :: ai(:)
  PetscReal, pointer, intent(out) :: ar(:)

  PetscBool, intent(out) :: is_mult_only
  PetscBool, intent(out) :: is_depth
  PetscBool, intent(out) :: is_perm

  PetscBool, intent(out) :: qerr

  isint  = PETSC_FALSE
  isreal = PETSC_FALSE

  qerr = PETSC_FALSE

  is_mult_only = PETSC_FALSE
  is_depth     = PETSC_FALSE

  nullify(ai, ar)

! Cell dimensions

  if (StringCompareIgnoreCase(za, 'dx'    )) ar => g_dx
  if (StringCompareIgnoreCase(za, 'dy'    )) ar => g_dy
  if (StringCompareIgnoreCase(za, 'dz'    )) ar => g_dz

!  Permeabilities

  if (StringCompareIgnoreCase(za, 'permx' )) then
    ar => g_kx
    is_perm = PETSC_TRUE
  endif

  if (StringCompareIgnoreCase(za, 'permy' )) then
    ar => g_ky
    is_perm = PETSC_TRUE
  endif

  if (StringCompareIgnoreCase(za, 'permz' )) then
    ar => g_kz
    is_perm = PETSC_TRUE
  endif

! Multipliers

  if (StringCompareIgnoreCase(za, 'multx' )) ar => g_mx
  if (StringCompareIgnoreCase(za, 'multy' )) ar => g_my
  if (StringCompareIgnoreCase(za, 'multz' )) ar => g_mz

! Porosity

  if (StringCompareIgnoreCase(za, 'poro'  )) ar => g_poro

! Tops (note this is a depth, so mark as such)

  if (StringCompareIgnoreCase(za, 'tops'  )) then
    ar => g_tops
    is_depth = PETSC_TRUE
  endif

! Net to gross

  if (StringCompareIgnoreCase(za, 'ntg'   )) ar => g_ntg

!  Integer actnum and satnum

  if (StringCompareIgnoreCase(za, 'actnum')) ai => g_actn
  if (StringCompareIgnoreCase(za, 'satnum')) ai => g_satn

!  Transmissibilities (can only be multiplied)

  if (StringCompareIgnoreCase(za, 'tranx' )) then
    ar => g_mx
    is_mult_only = PETSC_TRUE
  endif

  if (StringCompareIgnoreCase(za, 'trany' )) then
    ar => g_my
    is_mult_only = PETSC_TRUE
  endif

  if (StringCompareIgnoreCase(za, 'tranz' )) then
    ar => g_mz
    is_mult_only = PETSC_TRUE
  endif

!  Pore volumes (can only be multiplied)

  if (StringCompareIgnoreCase(za, 'porv'  )) then
    ar => g_poro
    is_mult_only = PETSC_TRUE
  endif

! Set up the type

  if (associated(ai)) then
    isint  = PETSC_TRUE
  else if (associated(ar)) then
    isreal = PETSC_TRUE
  else
    call SetError('Unable to find array' // trim(za), qerr)
  endif

end subroutine GetGridArrayPointer

! *************************************************************************** !

subroutine ProcessArgToInt(iv, za, zk, il, iu, qerr)

  !
  ! Convert argument za of keyword zk to an integer
  ! Note iv is left unchanged if the argument is blank
  !
  ! Author: Dave Ponting
  ! Date: 02/11/19
  !

  implicit none

  PetscInt, intent(inout) :: iv
  character(len=*), intent(in) :: za
  character(len=*), intent(in) :: zk
  PetscInt, intent(in) :: il, iu
  PetscBool, intent(out) :: qerr

  PetscInt :: ierr

!  Initialise values

  ierr = 0
  qerr = PETSC_FALSE

!  Do the read operation

  if (len_trim(za) > 0) then
    read(za, *, iostat = ierr) iv
  endif

!  Check for errors

  if (ierr /= 0) then
    call SetError(trim(zk) // ' (error in reading integer value)', qerr)
  else if (iv < il) then
    call SetError(trim(zk) // ' (value below lower bound)' // trim(zk), qerr)
  else if (iv > iu) then
    call SetError(trim(zk) // ' (value above upper bound)' // trim(zk), qerr)
  endif

end subroutine ProcessArgToInt

! *************************************************************************** !

subroutine ProcessArgToReal(rv, za, zk, qerr)

  !
  ! Convert argument za of keyword zk to a real
  ! Note rv is left unchanged if the argument is blank
  !
  ! Author: Dave Ponting
  ! Date: 02/11/19
  !

  implicit none

  PetscReal, intent(inout) :: rv
  character(len=*), intent(in) :: za
  character(len=*), intent(in) :: zk
  PetscBool, intent(out) :: qerr

  PetscInt  :: ierr

!  Initialise values

  ierr = 0
  qerr = PETSC_FALSE

!  Read value

  if (len_trim(za) > 0) then
    read(za, *, iostat = ierr) rv
  endif

!  Check for errors

  if (ierr /= 0) then
    call SetError( trim(zk) // ' (error in reading real value)', qerr)
  endif

end subroutine ProcessArgToReal

! *************************************************************************** !

subroutine CheckError(input, zerr, qerr)
  !
  ! Check the error code on the input stream
  ! Set an error message and return qerr as true if an error is found
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  type(input_type), pointer :: input
  character(len = *) :: zerr
  PetscBool, intent(inout) :: qerr

  PetscInt :: ierr

  ierr = input%ierr

  if (ierr == 0) then
    qerr = PETSC_FALSE
  else
    call SetError(zerr, qerr)
  endif

end subroutine CheckError

! *************************************************************************** !

subroutine SetError(zerr, qerr)
  !
  ! Set an error message and flag
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  character(len = *) :: zerr
  PetscBool, intent(out) :: qerr

  g_error_flag = 1
  g_error_string = zerr
  qerr = PETSC_TRUE

end subroutine SetError

! *************************************************************************** !

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
  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3)
  PetscReal, allocatable :: bvg(:)

  ! Set up bulk volume by grid index work array

  allocate(bvg(g_nxyz));bvg =  0.0

  !  Process faults

  call ProcessFaults()

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

  print *, 'Calculating pore volumes'

  g_na = 0
  do ix = 1, g_nx
    do iy = 1, g_ny
      do iz = 1, g_nz

        ig = GetNaturalIndex(ix, iy, iz)

        ! Find pore volume

        bvol = 0.0
        porv = 0.0

        if (g_iscpg) then
          call GetCorners( ix, iy, iz, &
                           x000, x100, x010, x110, &
                           x001, x101, x011, x111, &
                           g_coord, g_zcorn, g_nx, g_ny )
          bvg(ig) = findVolume( x000, x100, x010, x110, &
                                x001, x101, x011, x111 )
        else
          bvg(ig) = g_dx(ig)*g_dy(ig)*g_dz(ig)
        endif

        bvol = bvg(ig)

        ! iact is the ACTNUM array value

        iact = g_actn(ig)

        if (iact == 1) then
          porv = g_poro(ig)*g_ntg(ig)*bvol
        endif
        if (iact == 2) then ! Rock volume only, tiny pore vol, no fluid flow
          porv     = eps
          g_kx(ig) = 0.0
          g_ky(ig) = 0.0
          g_kz(ig) = 0.0
        endif
        if (iact == 3) then ! Pore volume only, but ntg honoured
          porv = bvol*g_ntg(ig)
        endif

        if (                        (porv > g_minpv(1))  &
            .or. ((iact == 2) .and. (bvol > 0.0       )) ) then
          g_na = g_na+1
          g_gtoa(ig) = g_na
        endif

      enddo
    enddo
  enddo

  ! Only continue if some active cells

  if (g_na>0) then

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

    ! Set up atoc (loop over cells in Eclipse compressed natural order)

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

    ! Initial estimate for number of connections
    ! Sufficient for normal connections + some space for wells

    g_mc = 3*g_na
    g_nc = 0

    ! Allocate active arrays

    call AllocateActiveArrays()

    ! Find the cell locations and connections

    call GenerateGridConnections(bvg)

  endif

  !  Deallocate the bulk volume by grid indexc work array

  if (allocated(bvg)) deallocate(bvg)

end subroutine ProcessGridData

! *************************************************************************** !

subroutine ProcessFaults()

  !
  ! Process the FAULTS/MULTFLT data
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18
  !

  use String_module

  implicit none

  PetscInt  :: imultflt, ifault, il, iu, jl, ju, kl, ku, i, j, k, ig, jg, iface
  PetscReal :: vm
  character(len=8) :: zf, zname

  do imultflt = 1, g_nmultflt_data

    zf = g_multflt_data(imultflt)%zname
    vm = g_multflt_data(imultflt)%vm

   do ifault = 1, g_nfault_data

    zname = g_fault_data(ifault)%zname

    if (StringCompareIgnoreCase(zf, zname)) then

      il = g_fault_data(ifault)%il
      iu = g_fault_data(ifault)%iu

      jl = g_fault_data(ifault)%jl
      ju = g_fault_data(ifault)%ju

      kl = g_fault_data(ifault)%kl
      ku = g_fault_data(ifault)%ku

      iface = g_fault_data(ifault)%iface

      do i = il, iu
        do j = jl, ju
          do k = kl, ku

            ig =  GetNaturalIndex(i, j, k)

            ! Positive values set MULTX/Y/Z for the fault cell

            if (iface == g_facei) g_mx(ig) = g_mx(ig)*vm
            if (iface == g_facej) g_my(ig) = g_my(ig)*vm
            if (iface == g_facek) g_mz(ig) = g_mz(ig)*vm

            ! Negative values set MULTX/Y/Z for the neighbouring cell

            if (iface == g_faceim .and. i>1) then
              jg =  GetNaturalIndex(i-1, j, k)
              g_mx(jg) = g_mx(jg)*vm
            endif

            if (iface == g_facejm .and. j>1) then
              jg =  GetNaturalIndex(i, j-1, k)
              g_my(jg) = g_my(jg)*vm
            endif

            !  Care needed to modify MULTZ for correct cell:
            !  Lower in Eclipse terms is upper in Pflotran terms

            if (iface == g_facekm .and. k<g_nz) then
              jg =  GetNaturalIndex(i, j, k+1)
              g_mz(jg) = g_mz(jg)*vm
            endif

          enddo
        enddo
      enddo

    endif

   enddo

  enddo

end subroutine processFaults

! *************************************************************************** !

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
         call SetError(zmess, qerr)
         exit outer
      endif

      if ((iy < 1) .or. (iy > g_ny)) then
         zmess = 'Completion J-location out of range, well ' // trim(wname)
         call SetError(zmess, qerr)
         exit outer
      endif

      if ((iz < 1) .or. (iz > g_nz)) then
         zmess = 'Completion K-location out of range, well ' // trim(wname)
         call SetError(zmess, qerr)
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

! *************************************************************************** !

subroutine DistributeWells(option)
  !
  ! Distribute the well information from the I/O rank to other ranks
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
    if (g_mcmpl_data > 0) then
      if (allocated(g_cmpl_data)) then
        deallocate(g_cmpl_data)
      endif
    endif
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

! *************************************************************************** !

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

        ! if the cells reside on io_rank, load them up

        do ial = 1, narank
          ia = ial + iabase
          explicit_grid%cell_ids      (ial)   =        ia
          explicit_grid%cell_centroids(ial)%x = g_x   (ia)
          explicit_grid%cell_centroids(ial)%y = g_y   (ia)
          explicit_grid%cell_centroids(ial)%z = g_z   (ia)
          explicit_grid%cell_volumes  (ial)   = g_bvol(ia)
        enddo

      else

        ! if the cells reside on another rank, send them

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
      explicit_grid%cell_ids      (ia)   = nint(temp_real_array(1, ia))
      explicit_grid%cell_centroids(ia)%x =      temp_real_array(2, ia)
      explicit_grid%cell_centroids(ia)%y =      temp_real_array(3, ia)
      explicit_grid%cell_centroids(ia)%z =      temp_real_array(4, ia)
      explicit_grid%cell_volumes  (ia)   =      temp_real_array(5, ia)
    enddo

  endif

  deallocate(temp_real_array)

end subroutine DistributeCells

! *************************************************************************** !

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
      explicit_grid%connections   (1, icl)   = nint(temp_real_array(1, icl))
      explicit_grid%connections   (2, icl)   = nint(temp_real_array(2, icl))
      explicit_grid%face_centroids(   icl)%x =      temp_real_array(3, icl)
      explicit_grid%face_centroids(   icl)%y =      temp_real_array(4, icl)
      explicit_grid%face_centroids(   icl)%z =      temp_real_array(5, icl)
      explicit_grid%face_areas    (   icl)   =      temp_real_array(6, icl)
    enddo

  endif

  deallocate(temp_real_array)

end subroutine DistributeConnections

! *************************************************************************** !

subroutine DistributeCounts(option)
  !
  ! Broadcast various counts from the io proc to all the other procs
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !

  implicit none

  type(option_type) :: option

  call BroadcastInt(g_nxyz   , option)
  call BroadcastInt(g_na     , option)
  call BroadcastInt(g_rsatn  , option)
  call BroadcastInt(g_maxsatn, option)

end subroutine DistributeCounts

! *************************************************************************** !

subroutine CreateElements(unstructured_grid, explicit_grid)
  !
  ! Set up the set locations as hexahedra in xyz space
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !

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

! *************************************************************************** !

subroutine SetDimens(dimens)
  !
  ! Dimens read, so allocate the grid-sized arrays
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !

  implicit none
  PetscReal, intent(in) :: dimens(3)

  g_nx     = nint(dimens(1))
  g_ny     = nint(dimens(2))
  g_nz     = nint(dimens(3))

  g_nxp    = g_nx + 1
  g_nyp    = g_ny + 1
  g_nxy    = g_nx*g_ny
  g_nxpnyp = g_nxp*g_nyp
  g_nxyz   = g_nx*g_ny*g_nz

  g_dimens_read = PETSC_TRUE

  call AllocateGridArrays()

end subroutine SetDimens

! *************************************************************************** !

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

! *************************************************************************** !

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
  allocate(g_satn(g_nxyz));g_satn =  1 ! Satnum (used to set sat. tab. numbers)
  allocate(g_gtoa(g_nxyz));g_gtoa = -1

end subroutine AllocateGridArrays

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

subroutine GenerateGridConnections(bvg)
  !
  ! Given the cell locations, generate the connections
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscReal, allocatable, intent(in) :: bvg(:)

  PetscInt  :: ix, iy, iz, ig, ia
  PetscReal :: x, y, z, dx, dy, dz, mx, my, ntgi

  ! Set up grid connections

  print *, 'Calculating connections'

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

          ntgi = g_ntg(ig)

          g_bvol(ia) = bvg(ig)

          g_x(ia) = g_xloc(ig)
          g_y(ia) = g_yloc(ig)
          g_z(ia) = g_zloc(ig)

          x = g_xloc(ig)
          y = g_yloc(ig)
          z = g_zloc(ig)

          ! x-direction

          if (ix < g_nx) then
            if (g_isnewtran) then
              call getArea1( ia, ix, iy, iz, g_xdir, mx, my, ntgi )
            else
              call getArea0( ia, ig, ix, iy, iz, g_xdir, &
                             dx, dy, dz, mx, my, ntgi )
            endif
          endif

          ! y-direction

          if (iy< g_ny) then
            if (g_isnewtran) then
              call getArea1( ia, ix, iy, iz, g_ydir, mx, my, ntgi )
            else
              call getArea0( ia, ig, ix, iy, iz, g_ydir, &
                             dx, dy, dz, mx, my, ntgi )
            endif
          endif

          ! z-direction

          if (iz < g_nz) then
            if (g_isnewtran) then
              call getArea1( ia, ix, iy, iz, g_zdir, mx, my, ntgi )
            else
              call getArea0( ia, ig, ix, iy, iz, g_zdir, &
                             dx, dy, dz, mx, my, ntgi )
            endif
          endif

        endif

      enddo ! enddo iz
    enddo ! enddo iy
  enddo ! enddo ix

  print *, 'Ncell=', g_nxyz, ' Nact=', g_na, ' Nconn=', g_nc, &
           ' Nflt=', g_nf  , ' Npo=' , g_npo

end subroutine GenerateGridConnections

! *************************************************************************** !

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
  call DeallocateArray(g_satn)

  call DeallocateArray(g_gtoa)

end subroutine DeallocateGridArrays

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

function GetSatnumValue(ia)
  !
  ! Get the satnum values of a given active cell
  !
  ! Author: Dave Ponting
  ! Date: 02/14/19

  implicit none

  PetscInt :: GetSatnumValue

  PetscInt, intent(in)  :: ia
  PetscInt :: ig

  ! Get the grid order (ix-fastest)

  ig = g_atog(ia)

  ! Extract and return satnum

  GetSatnumValue  = g_satn(ig)

end function GetSatnumValue

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

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
  PetscReal :: dx(3), dy(3), dz(3), prx, pry

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

  ! Find dx and dy scalar distances from x- and y-direction projections

  prx = dx(g_xdir)
  pry = dx(g_ydir)
  vdx = sqrt(prx**2+pry**2)

  prx = dy(g_xdir)
  pry = dy(g_ydir)
  vdy = sqrt(prx**2+pry**2)

  ! Find dz from the block vertical thickness

  vdz = abs(dz(g_zdir))

end subroutine GetHexDims

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

  PetscBool :: started,inquotes,quote

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
    started      = PETSC_FALSE
    inquotes     = PETSC_FALSE

    startcol = g_column
    do i = startcol, l
      c = line(i:i)
      iws = IsWhiteSpace(c)
      ic = ichar(c)
      quote=IsQuote(c)
      if ( quote ) inquotes = .not.inquotes
      if ((.not. started) .and. (.not.iws)) started = PETSC_TRUE
      if (       started  .and. ((iws.and. (.not.inquotes)) .or. (i == l) )) exit
      if (started) then
        if (.not.quote) then
          j = j + 1
          word(j:j) = c
          somethingRead = PETSC_TRUE
        endif
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
        call InputReadPflotranStringNotComment(input, option)
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

function IsQuote(c)
  !
  ! Indicates that a character is a quote (' or ")
  ! Author: Dave Ponting
  ! Date: 02/11/19

  implicit none

  PetscBool :: isQuote
  character, intent(in) :: c
  PetscInt :: ic

  IsQuote = PETSC_FALSE

  ic = ichar(c)
  if ((ic == g_squote) .or. &
      (ic == g_dquote) ) isQuote = PETSC_TRUE

end function IsQuote

! *************************************************************************** !

subroutine ReadEGridArrI(a, keyword, input, option, is_nn, qerr)
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
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(in) :: is_nn
  PetscBool, intent(inout) :: qerr

  PetscInt , allocatable :: column_buffer(:)
  PetscReal, allocatable :: ar(:)
  PetscInt  :: ix, iy, iz, izpft, ig, igpft, i

  ! Do the actual read operation

  allocate(ar(g_nxyz))
  ar = 1.0

  call ReadEvalues(ar, g_nxyz, keyword, 'GRID' , input, option, is_nn, qerr)

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

      ! Place into Pflotran-conversion array with appropriate data modifications

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

end subroutine ReadEGridArrI

! *************************************************************************** !

subroutine ReadEGridArrR(a, keyword, &
                         input, option, is_dep, is_perm, is_nn, qerr)
  !
  ! Reads an Eclgrid grid array
  ! If it is a depth array, then flip sign to Pflotran elevation convention
  ! If is is a permeability array, convert from mD to m2
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscReal, intent(inout) :: a(:)
  character(len = *) :: keyword
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(in) :: is_dep
  PetscBool, intent(in) :: is_perm
  PetscBool, intent(in) :: is_nn
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
      call SetError(keyword, qerr)
    else
      call ReadEvalues(a, nread, keyword, 'GRID', input, option, is_nn, qerr)
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

end subroutine ReadEGridArrR

! *************************************************************************** !

subroutine ReadECoordArray(a, input, option, qerr)
  !
  ! Reads an Eclgrid coord array
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscReal, intent(inout) :: a(:)
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr
  PetscBool, parameter :: nn_no = PETSC_FALSE

  PetscInt  :: i, nread

  nread = size(a)
  if (nread /= 6*g_nxpnyp) then
    call SetError('Coord read', qerr)
  else
    call ReadEvalues(a, nread, 'COORD', 'GRID', input, option, nn_no, qerr)
  endif

  do i = 3, nread, 3
    a(i) = z_flip*a(i)
  enddo

end subroutine ReadECoordArray

! *************************************************************************** !

subroutine ReadEZcornArray(a, input, option, qerr)
  !
  ! Reads an Eclgrid zcorn array
  ! If it is a depth array, then flip sign to Pflotran elevation convention
  ! If is is a permeability array, convert from mD to m2
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscReal, intent(inout) :: a(:)
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr
  PetscBool, parameter :: nn_no = PETSC_FALSE


  PetscReal, allocatable :: column_buffer(:)
  PetscInt  :: inx , iny, inz, nnx, nny, nnz, nnxy, inode, nread, inzpft

  qerr = PETSC_TRUE

  ! Do the actual read operation

  nread = size(a)
  if (nread /= 8*g_nxyz) then
    call SetError('Zcorn read', qerr)
  else
    call ReadEvalues(a, nread, 'ZCORN', 'GRID', input, option, nn_no, qerr)
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

subroutine ReadEvalues(a, n, keyword, section, input, option, is_nn, qerr)
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
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(in) :: is_nn
  PetscBool, intent(inout) :: qerr

  PetscInt  :: i, iostat, istar, nstack, ierr
  PetscReal :: dval
  PetscBool :: exittime
  character(len = MAXWORDLENGTH) :: word, repc, hold
  character(len = MAXWORDLENGTH) :: zmess

  PetscReal, parameter :: vmargin = -1.0E-10

  qerr   = PETSC_FALSE

  i      = 0
  iostat = 0
  dval   = 0.0
  exittime = PETSC_FALSE
  nstack = 0
  ierr   = 0
  zmess = ' '

  call InputReadPflotranStringNotComment(input, option)
  if (InputCheckExit(input, option)) exittime = PETSC_TRUE

  if (.not. exittime) then

    g_column = 1
    do i = 1, n

      if (nstack > 0) then

        ! Case of value in stack

        a(i)   = dval
        nstack = nstack - 1

      else

        ! Case of no value in stack

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

      if (is_nn) then
        if (a(i) < vmargin) then
          zmess = trim(keyword) //', value ' // trim(word) &
                                // ' cannot be negative'
          call SetError(zmess, qerr)
        endif
      endif

    enddo

  endif

  if (ierr == 1 ) then
    zmess = 'Unable to read ' // trim(keyword)
    call SetError(zmess, qerr)
  endif

end subroutine ReadEvalues

! *************************************************************************** !

subroutine ReadEstrings(zkey, a, n, input, option, qerr)
  !
  ! Read a series of strings from an Eclipse syntax file
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  character(len = MAXWORDLENGTH), intent(in) :: zkey
  character(len = MAXWORDLENGTH), intent(inout) :: a(:)
  PetscInt , intent(in)    :: n
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr

  PetscInt  :: i, iostat, istar, nstack, m, nr, ierr
  PetscBool :: exittime
  character(len = MAXWORDLENGTH) :: word
  character(len = MAXWORDLENGTH) :: repc
  character(len = MAXWORDLENGTH) :: hold
  character(len = MAXWORDLENGTH) :: zval
  character(len = MAXWORDLENGTH) :: zmess

  ierr   = 0
  qerr   = PETSC_FALSE

  i      = 0
  iostat = 0
  zval   = ' '
  exittime = PETSC_FALSE
  nstack = 0

  call InputReadPflotranStringNotComment(input, option)
  if (input%ierr /= 0) then
    exittime = PETSC_TRUE
  else
    if (InputCheckExit(input, option)) exittime = PETSC_TRUE
  endif

  if (.not. exittime) then

    g_column = 1
    m  = size(a)
    nr = min(n, m)

    do i = 1, nr

      if (nstack > 0) then

        ! Case of value in stack

        a(i)   = zval
        nstack = nstack - 1

      else

        ! Case of no value in stack

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

        ! Store the word

        zval = word
        a(i) = zval

      endif

    enddo

  endif

  if (ierr == 1 ) then
    zmess = 'Unable to read ' // trim(zkey)
    call SetError(zmess, qerr)
  endif

end subroutine ReadEstrings

! *************************************************************************** !

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

! *************************************************************************** !

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

! *************************************************************************** !

subroutine CopyFaultData(fto, ffrom)
  !
  ! Copy completion data from one structure to another
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  type(fault_data_type), intent(out) :: fto
  type(fault_data_type), intent(in ) :: ffrom

  fto%zname = ffrom%zname

  fto%il    = ffrom%il
  fto%iu    = ffrom%iu
  fto%jl    = ffrom%jl
  fto%ju    = ffrom%ju
  fto%kl    = ffrom%kl
  fto%ku    = ffrom%ku
  fto%iface = ffrom%iface

end subroutine CopyFaultData

! *************************************************************************** !

subroutine CopyMultfltData(fto, ffrom)
  !
  ! Copy multflt data from one structure to another
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  type(multflt_data_type), intent(out) :: fto
  type(multflt_data_type), intent(in ) :: ffrom

  fto%zname = ffrom%zname
  fto%vm    = ffrom%vm

end subroutine CopyMultfltData

! *************************************************************************** !

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

! *************************************************************************** !

function GetGrdNCmpl(iw)
  !
  ! Given a well name, find the number of completions
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscInt             :: GetGrdNCmpl
  PetscInt, intent(in) :: iw

 if ((iw>0) .and. (iw <= g_nwell_data)) then
   GetGrdNCmpl = g_well_data(iw)%iku - g_well_data(iw)%ikl + 1
 else
   GetGrdNCmpl = 0
 endif

end function GetGrdNCmpl

! *************************************************************************** !

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

! *************************************************************************** !

subroutine IsCPG()
  !
  ! If a keyword encountered indicating corner point input, (ZORN or COORD)
  ! set flags and allocate the ccord and zcorn arrays
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  g_iscpg = PETSC_TRUE

  if (.not.g_cpgallocated) then
    allocate(g_coord(6*g_nxpnyp));g_coord =  0.0
    allocate(g_zcorn(8*g_nxyz  ));g_zcorn =  0.0
    g_cpgallocated = PETSC_TRUE
  endif

end subroutine isCPG

! *************************************************************************** !

subroutine PermPoroExchangeAndSet(poro_p, permx_p, permy_p, permz_p, &
                                  permxy_p, permxz_p, permyz_p, &
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
  PetscReal, pointer :: permxy_p (:)
  PetscReal, pointer :: permxz_p (:)
  PetscReal, pointer :: permyz_p (:)
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
          permxy_p(ilt) = 0.d0
          permxz_p(ilt) = 0.d0
          permyz_p(ilt) = 0.d0
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

        ! Receive the natural address array from irank

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

    ! Barrier call to stop other procs running ahead

    call MPI_Barrier(option%mycomm, ierr)
  enddo

end subroutine PermPoroExchangeAndSet

! *************************************************************************** !

subroutine SatnumExchangeAndSet(satnum, inatsend, nlmax, nl2g, option)
  !
  ! Satnum (and potentially other arrays like PVTNUM or IMBNUM)
  ! values from the grdecl file are known on the IO rank
  ! These are needed on the other ranks, but only for the cells on those ranks
  ! To avoid over-sending, the other ranks send lists to the IO ranks,
  ! and the IO ranks send back the required values.
  !
  ! Author: Dave Ponting
  ! Date: 02/21/19

  implicit none

  PetscInt , pointer :: satnum  (:)
  PetscInt , pointer :: inatsend(:)
  PetscInt, intent(in) :: nlmax
  PetscInt, intent(in) :: nl2g(:)
  type(option_type) :: option

  PetscInt :: irank, iorank, t_rank
  PetscInt :: ilt, ilo, ino, igt
  PetscInt :: nlo, ierr
  PetscMPIInt :: int_mpi, temp_int_array(1), status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt, parameter :: tag_mpi = 0

  ! Work arrays

  PetscInt, allocatable :: winat(:)
  PetscInt, allocatable :: wsatn(:)

  ! Set scalars

  ierr   = 0
  iorank = option%io_rank
  t_rank = option%myrank

  ! Loop over exchange operations between IO rank and non-IO ranks (irank)

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

        ! Allocate buffers to hold response satnum values from ioproc

        allocate(wsatn(nlmax))

        ! Receive satnum values from ioproc

        int_mpi = nlmax
        call MPI_Recv(wsatn, int_mpi, MPI_INTEGER, iorank, &
                      MPI_ANY_TAG, option%mycomm, status_mpi, ierr)

        ! Copy buffers into correct storage locations

        do ilt = 1, nlmax
          igt = nl2g(ilt)
          satnum(igt) = wsatn(ilt)
        enddo

        ! Free buffers

        deallocate(wsatn)

      endif

      if (t_rank  == iorank) then

       !This is IO rank - receive the other rank nlmax value (nlo) from irank

        temp_int_array(1) = 0
        call MPI_Recv(temp_int_array, ONE_INTEGER_MPI, MPI_INTEGER, &
                      irank, MPI_ANY_TAG, option%mycomm, status_mpi, ierr)
        nlo = temp_int_array(1)

        ! Allocate work array to hold the received natural addresses

        allocate(winat(nlo))
        allocate(wsatn(nlo))

        ! Receive the natural address array from irank

        int_mpi = nlo
        call MPI_Recv(winat, int_mpi, MPI_INTEGER, irank, MPI_ANY_TAG, &
                      option%mycomm, status_mpi, ierr)

        ! On this (io) proc, set up the satnum values to be returned

        do ilo = 1, nlo
          ino = winat(ilo)
          wsatn(ilo) = GetSatnumValue(ino)
        enddo

        ! Send back the satn to irank

        int_mpi = nlo
        call MPI_Send(wsatn, int_mpi, MPI_INTEGER, &
                      irank, tag_mpi, option%mycomm, ierr)

        ! Free the work arrays

        deallocate(winat)
        deallocate(wsatn)

      endif
    endif

    ! Barrier call to stop other procs running ahead

    call MPI_Barrier(option%mycomm, ierr)
  enddo

end subroutine SatnumExchangeAndSet

! *************************************************************************** !

function findVolume( x000, x100, x010, x110, &
                     x001, x101, x011, x111 )

  !
  ! Obtain the volume of a hexahederal grid cell
  ! with corners x000,..,x111
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscReal findVolume
  PetscReal, intent(in) :: x000(3), x100(3), x010(3), x110(3), &
                           x001(3), x101(3), x011(3), x111(3)
  PetscInt  :: id
  PetscReal :: r000, r100, &
               r010, r110, &
               r001, r101, &
               r011, r111
  PetscReal :: c(0:1, 0:1, 0:1, 3), vp0, vp1, vp2, vp3, vp4, vp5, vp

  ! Loop over xyz directions loading the coefficients such that:
  !
  ! r = r000 + a  *r100 + b  *r010 + c  *r001
  !          + a*b*r110 + b*g*r110 + g*a*1001 + a*b*g*r111
  !
  ! where r is the location in cell if alpha, beta, gamma each lie in (0,1)

  do id = 1, 3

    r000 = x000(id)
    r100 = x100(id)
    r010 = x010(id)
    r110 = x110(id)
    r001 = x001(id)
    r101 = x101(id)
    r011 = x011(id)
    r111 = x111(id)

    c(0, 0, 0, id) = r000
    c(1, 0, 0, id) = r100 - r000
    c(0, 1, 0, id) = r010 - r000
    c(0, 0, 1, id) = r001 - r000

    c(1, 1, 0, id) = r110 - (r100 + r010) + r000
    c(0, 1, 1, id) = r011 - (r010 + r001) + r000
    c(1, 0, 1, id) = r101 - (r100 + r001) + r000

    c(1, 1, 1, id) = r111 - (r110 + r101 + r011) + (r100 + r010 + r001) - r000

  enddo

  !  Form the 6 components of the volume corresponding
  !  to the 6 components of the Jacobian d(xyz)/d(abg)

  vp0 = findVolume1(1, 2, 3, c)
  vp1 = findVolume1(1, 3, 2, c)

  vp2 = findVolume1(2, 1, 3, c)
  vp3 = findVolume1(2, 3, 1, c)

  vp4 = findVolume1(3, 1, 2, c)
  vp5 = findVolume1(3, 2, 1, c)

  !  Sum with appropriate signs into the final pore volume

  vp = vp0 - vp1 - ( vp2 - vp3 ) + vp4 - vp5

  ! Take absolute value to get pore volume

  findvolume = abs(vp)

end function findVolume

! *************************************************************************** !

function findVolume1(id0, id1, id2, c)

  !
  ! Obtain a volume contribution integral(J) where J
  ! one of the six permutation elements of d(px,py,pz)/d(a,b,g)
  ! where px,py,pz are a permutation of x,y,z indicated by id0,id1,id2
  ! So id0,id1,id2=1,2,3 indicates px,py,pz=x,y,z
  !    id0,id1,id2=3,2,1 indicates px,py,pz=z,y,x etc.
  ! The permutation sign is added externally.
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscReal :: findVolume1
  PetscInt , intent(in) :: id0, id1, id2
  PetscReal, intent(in) :: c(0:1, 0:1, 0:1, 3)

  PetscInt :: ja, jb, jc, jd, je, jf
  PetscReal :: den

  findVolume1 = 0.0

  do ja = 0, 1
    do jb = 0, 1
      do jc = 0, 1
        do jd = 0, 1
          do je = 0, 1
            do jf = 0, 1

              den = real( (jc + je + 1) * (ja + jf + 1) * (jb + jd + 1) )
              findVolume1 = findVolume1             &
                           + c(1 , ja, jb, id0)     &
                           * c(jc, 1 , jd, id1)     &
                           * c(je, jf, 1 , id2) / den

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end function  findVolume1

! *************************************************************************** !

subroutine getArea0(ia, ig, ix, iy, iz, idir, &
                    dxi, dyi, dzi, mx, my, ntgi)
  !
  ! Do a block-based area calculation between cell (ix,iy,iz) and its
  ! positive neighbour in the idir-direction (g_xdir,..,g_z_dir)
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt , intent(in) :: ia, ig, ix, iy, iz, idir
  PetscReal, intent(in) :: dxi, dyi, dzi, mx, my, ntgi

  PetscInt  :: jx, jy, jz, jg, ja
  PetscReal :: xf, yf, zf, a, ai, aj, mz, ntgj, dxj, dyj, dzj, wi, wj, dd, dh
  PetscReal :: dxij, dyij, dzij, dip, dipc

  !  Default mz

  mz = 1.0

  ! Find neighbour

  jx = ix
  jy = iy
  jz = iz

  if (idir == g_xdir) jx = jx + 1
  if (idir == g_ydir) jy = jy + 1
  if (idir == g_zdir) jz = jz + 1

  jg = GetNaturalIndex(jx, jy, jz)
  ja = g_gtoa(jg)

  if (ja > 0) then

!   If z-dir, multz of lower cell in Eclipse is from cell above in Pflotran

    if (idir == g_zdir) mz = g_mz(jg)

!   Get ntg for cell j

    ntgj = g_ntg(jg)

!   Get size values in each direction

    dxj  = g_dx(jg)
    dyj  = g_dy(jg)
    dzj  = g_dz(jg)

    dxij = dxi + dxj
    dyij = dyi + dyj
    dzij = dzi + dzj

    ! Weights (if dxi>>dxj, wi->1, wj->0)

    wi = 0.5
    wj = 0.5

    ! Depth difference for dip corrections

    dd   = abs(g_zloc(ig)-g_zloc(jg))
    dipc = 1.0

    ! x-direction

    if (idir == g_xdir) then
      if (dxij>0.0) then
        dh   = 0.5*dxij
        dip  = atan(dd/dh)
        dipc = cos(dip)
        wi   = dxi/dxij
        wj   = dxj/dxij
      endif
      ai = dyi * dzi * ntgi
      aj = dyj * dzj * ntgj
      a  = mx * (ai*wj + aj*wi)*dipc
    endif

    if (idir == g_ydir) then
      if (dyij>0.0) then
        dh   = 0.5*dyij
        dip  = atan(dd/dh)
        dipc = cos(dip)
        wi   = dyi/dyij
        wj   = dyj/dyij
      endif
      ai  = dzi * dxi * ntgi
      aj  = dzj * dxj * ntgj
      a  = my * (ai*wj + aj*wi)*dipc
    endif

    if (idir == g_zdir) then
      if (dzij>0.0) then
        wi = dzi/dzij
        wj = dzj/dzij
      endif
      ai = dxi * dyi
      aj = dxj * dyj
      a  = mz *  (ai*wj + aj*wi)
    endif

    ! Set up the interface along the line from i to j, weighted by cell size

    xf = wj*g_xloc(ig)+wi*g_xloc(jg)
    yf = wj*g_yloc(ig)+wi*g_yloc(jg)
    zf = wj*g_zloc(ig)+wi*g_zloc(jg)

    ! Add in the connection

    if (a>0.0) then
      call AddConnection(ia, ja, xf, yf, zf, a)
    endif

  endif

end subroutine getArea0

! *************************************************************************** !

subroutine getArea1(ia, ix, iy, iz, idir, mx, my, ntgi)
  !
  ! Do a corner-point-based area calculation between cell (ix,iy,iz) and its
  ! positive neighbour in the id-direction (g_xdir,..,g_z_dir)
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt , intent(in) :: ia, ix, iy, iz, idir
  PetscReal, intent(in) :: mx, my, ntgi

  PetscInt  :: jx, jy, jz, jg, ja, kz, kg, ka
  PetscReal :: a, xf, yf, zf, fm, fntg , af, pod
  PetscBool :: match, matchingNeighbourFound, found

  ! Set up natural neighbour location

  jx = ix
  jy = iy
  jz = iz

  if (idir == g_xdir) jx = jx + 1
  if (idir == g_ydir) jy = jy + 1
  if (idir == g_zdir) jz = jz + 1

  ! Initialise flags

  match                  = PETSC_FALSE
  matchingNeighbourFound = PETSC_FALSE
  found                  = PETSC_FALSE

  !  Set mult factor for first cell

  fm  = 1.0
  if (idir == g_xdir) fm = mx
  if (idir == g_ydir) fm = my

  ! Look at the natural neighbour

  jg = GetNaturalIndex(jx, jy, jz)

  !  Find average ntg and set ntg factor

  fntg = 1.0
  if (     (idir == g_xdir) &
      .or. (idir == g_ydir) ) then
    fntg  = 0.5*(ntgi+g_ntg(jg))
  endif

  ! Take z-mult from upper cell (downwards dirn, lower cell in Eclipse terms)

  if (idir == g_zdir) fm = g_mz(jg)

  ja = g_gtoa(jg)
  if (ja > 0) then
    call GetMif(ix, iy, iz, jx, jy, jz, idir, a, xf, yf, zf, match, found)
    if (match) matchingNeighbourFound = PETSC_TRUE
    if (found) then
      af = a*fm*fntg
      if (af>0.0) then
        call AddConnection(ia, ja, xf, yf, zf, af)
      endif
    endif
  endif

  ! Not a perfect match - look for other connections

  if (.not.matchingNeighbourFound ) then

    ! In x- and y-dir, look for displacement faults to the jx,jy column

    if (idir == g_xdir .or. idir == g_ydir) then
      do kz = 1, g_nz
        if (kz /= jz) then !  Have looked already
          kg = GetNaturalIndex(jx, jy, kz)

  !       Find ntg factor

          fntg  = 0.5*(ntgi+g_ntg(kg))

          ka = g_gtoa(kg)
          if (ka > 0) then
            call GetMif(ix, iy, iz, jx, jy, kz, idir, &
                        a, xf, yf, zf, match, found)
            if (found) then
              g_nf = g_nf + 1
              af = a*fm*fntg
              if (af > 0.0) then
                call AddConnection(ia, ka, xf, yf, zf, af)
              endif
            endif
          endif
        endif
      enddo
    endif

    !  In z-dir, look for pinch-outs from iz to cells above jz in ix,iy column

    if (idir == g_zdir .and. iz<=(g_nz-2) ) then

      ! Search up column from cell above z-neighbour

      column:do kz = iz+2, g_nz

        ! Look for active cell

        kg = GetNaturalIndex(ix, iy, kz)
        ka = g_gtoa(kg)
        if (ka > 0) then

          ! z-mult from upper cell (downwrds dirn, lower cell in Eclipse terms)

          fm = g_mz(kg)

          ! No ntg in vertical

          fntg = 1.0

          ! Find pinch-out distance

          call GetPod(ix, iy, iz, kz, pod)
          if (pod<g_pinch(1)) then
            call GetMif(ix, iy, iz, jx, jy, kz, idir, &
                        a, xf, yf, zf, match, found)
            if (found) then
              af = a*fm*fntg
              if (af > 0.0) then
                 g_npo = g_npo + 1
                call AddConnection(ia, ka, xf, yf, zf, af)
              endif
            endif
          endif

          !  Active cell found, finish search

          exit column
        endif
      enddo column
    endif

  endif

end subroutine getArea1

! *************************************************************************** !

subroutine GetMif(ix, iy, iz, jx, jy, jz, idir, a, xf, yf, zf, match, found)

  !
  ! Obtain the interface area between cell (ix,iy,iz) and cell (jx,jy,jz)
  ! in the direction idir
  ! Return the face centre and the scalar interface area
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt, intent(in) :: ix, iy, iz, jx, jy, jz, idir
  PetscReal, intent(out) :: a, xf, yf, zf
  PetscBool, intent(out) :: match, found

  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3), &
               fl(0:1, 0:1, 3), fu(0:1, 0:1, 3), &
               fcl(3), fcu(3), cl(3), cu(3), dlu(3), &
               d, ax, ay, az, ax1, ay1, az1, &
               srayl, srayu, srayl2, srayu2, &
               sdlu, sdlu2, elu(3), sdlui, f, sraylu, auns

  PetscReal, parameter :: missdistance = 0.01   ! 1 cm
  PetscReal, parameter :: minfacearea  = 0.0001 ! 1 cm2
  PetscBool, parameter :: positive_yes = PETSC_TRUE
  PetscBool, parameter :: positive_no  = PETSC_FALSE

  PetscBool :: missed

  PetscInt  :: i, j, k
  PetscReal :: diff,dl,du,topl,botl,topu,botu
  PetscReal :: axl,ayl,azl,axu,ayu,azu,areal,areau

  ! Default return values

     a  = 0.0
    xf  = 0.0
    yf  = 0.0
    zf  = 0.0
  match = PETSC_TRUE
  found = PETSC_FALSE

  ! Get the lower and upper cell centres and faces

  call GetCorners( ix, iy, iz, &
                   x000, x100, x010, x110, &
                   x001, x101, x011, x111, &
                   g_coord, g_zcorn, g_nx, g_ny )

  call GetCentre( x000, x100, x010, x110, &
                  x001, x101, x011, x111, cl )

  call GetFace( idir, positive_yes, &
                x000, x100, x010, x110, &
                x001, x101, x011, x111, fl )

  call GetCorners( jx, jy, jz, &
                   x000, x100, x010, x110, &
                   x001, x101, x011, x111, &
                   g_coord, g_zcorn, g_nx, g_ny )

  call GetCentre( x000, x100, x010, x110, &
                  x001, x101, x011, x111, cu )

  call GetFace( idir, positive_no, &
                x000, x100, x010, x110, &
                x001, x101, x011, x111, fu )

  ! Check for depth miss if x or y direction area

  missed = PETSC_FALSE

  if (idir /= g_zdir) then

    dl=fl(0, 0, g_zdir)
    du=fu(0, 0, g_zdir)

    topl = dl
    botl = dl
    topu = du
    botu = du

    do i = 0, 1
      do j = 0, 1
        dl=fl(i, j, g_zdir)
        du=fu(i, j, g_zdir)
        topl = max(topl,dl)
        botl = min(botl,dl)
        topu = max(topu,du)
        botu = min(botu,du)
      enddo
    enddo

    if ( (botl>topu) .or. (topl<botu) ) missed=PETSC_TRUE

  endif

  ! If not a depth miss, calculate area

  if (missed) then
    match=PETSC_FALSE
  else

    ! Store values

    do i = 1, 3
      dlu(i) = cu(i)-cl(i)
      fcl(i) = 0.25*(fl(0, 0, i)+fl(0, 1, i)+fl(1, 0, i)+fl(1, 1, i))
      fcu(i) = 0.25*(fu(0, 0, i)+fu(0, 1, i)+fu(1, 0, i)+fu(1, 1, i))
    enddo

    srayl2 = 0.0
    srayu2 = 0.0
    sdlu2  = 0.0

    do i = 1, 3
      d = fcl(i)-cl(i)
      srayl2 = srayl2 + d*d
      d = fcu(i)-cu(i)
      srayu2 = srayu2 + d*d
      sdlu2 = sdlu2 + dlu(i)*dlu(i)
    enddo

    srayl = sqrt(srayl2)
    srayu = sqrt(srayu2)
    sdlu  = sqrt(sdlu2)

    sdlui = 0.0
    if (sdlu > 0.0) sdlui = 1.0/sdlu

    do i = 1, 3
      elu(i) = dlu(i)*sdlui
    enddo

  !  Check for a match

    do i = 0, 1
      do j = 0, 1
        do k = 1, 3
          diff = fu(i, j, k)-fl(i, j, k)
          if (abs(diff)>missdistance) match = PETSC_FALSE
        enddo
      enddo
    enddo

  !  If a match, can use simple quad area calculation
  !  otherwise do mutual interface area calculation

    if (match) then

      ax = qarea(fl, g_xdir)
      ay = qarea(fl, g_ydir)
      az = qarea(fl, g_zdir)

    else

      ! Check that the individual quads are not of zero area

      axl = qarea(fl, g_xdir)
      ayl = qarea(fl, g_ydir)
      azl = qarea(fl, g_zdir)

      axu = qarea(fu, g_xdir)
      ayu = qarea(fu, g_ydir)
      azu = qarea(fu, g_zdir)

      areal=sqrt(axl*axl+ayl*ayl+azl*azl)
      areau=sqrt(axu*axu+ayu*ayu+azu*azu)

      ! If both have significant area, find overlap

      if ( (areal>minfacearea) .and. (areau>minfacearea) ) then
        ax = marea(fl, fu, g_xdir)
        ay = marea(fl, fu, g_ydir)
        az = marea(fl, fu, g_zdir)
      else
        ax = 0.0
        ay = 0.0
        az = 0.0
      endif

    endif

    !  Set up projection of area onto cell centre connection
    !  shifted pro rata the centre to face distances

    auns = ax*elu(g_xdir) + ay*elu(g_ydir) + az*elu(g_zdir)
    a  = abs(auns)
    sraylu = srayl+srayu
    f = 0.5
    if (sraylu > 0.0) then
      f = srayl/(srayl+srayu)
    endif

    xf = cl(g_xdir) + f*dlu(g_xdir)
    yf = cl(g_ydir) + f*dlu(g_ydir)
    zf = cl(g_zdir) + f*dlu(g_zdir)

    if (a>0.0) then
      found = PETSC_TRUE
    endif

  endif

end subroutine GetMif

! *************************************************************************** !

subroutine GetPod(ix, iy, iz, jz, pod)

  !
  ! Obtain gap between top of cell (ix,iy,iz) and bottom of cell (ix,iy,jz)
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt , intent(in)  :: ix, iy, iz, jz
  PetscReal, intent(out) :: pod

  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3), d, &
               fl(0:1, 0:1, 3), fu(0:1, 0:1, 3)

  PetscBool, parameter :: positive_yes = PETSC_TRUE
  PetscBool, parameter :: positive_no  = PETSC_FALSE

  PetscInt  :: i, j

  ! Get the lower and upper cell centres and faces

  call GetCorners( ix, iy, iz, &
                   x000, x100, x010, x110, &
                   x001, x101, x011, x111, &
                   g_coord, g_zcorn, g_nx, g_ny )

  call GetFace( g_zdir, positive_yes, &
                x000, x100, x010, x110, &
                x001, x101, x011, x111, fl )

  call GetCorners( ix, iy, jz, &
                   x000, x100, x010, x110, &
                   x001, x101, x011, x111, &
                   g_coord, g_zcorn, g_nx, g_ny )

  call GetFace( g_zdir, positive_no, &
                x000, x100, x010, x110, &
                x001, x101, x011, x111, fu )

  pod = 0.0
  do i = 0, 1
    do j = 0, 1
      d = abs(fu(i, j, g_zdir)-fl(i, j, g_zdir))
      pod = max(d, pod)
    enddo
  enddo

end subroutine GetPod

! *************************************************************************** !

subroutine AddConnection(ia, ja, xf, yf, zf, a)

  !
  !  Store connection between cells ia and ja, face location (xf,yf,zf), area a
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt , intent(in) :: ia, ja
  PetscReal, intent(in) :: xf, yf, zf, a

  PetscReal, parameter :: rm6 = 1.0E-6
  PetscInt :: ig, jg
  PetscReal :: xa, ya, za, xb, yb, zb, f, xt, yt, zt

  g_nc = g_nc + 1

  if (g_nc == g_mc) call ReallocateConnectionArrays()

  if (ja>ia) then
    g_cia(g_nc) = ia
    g_cja(g_nc) = ja
  elseif (ia>ja) then
    g_cia(g_nc) = ja
    g_cja(g_nc) = ia
  endif

  g_ccx(g_nc) = xf
  g_ccy(g_nc) = yf
  g_ccz(g_nc) = zf

  g_carea(g_nc) = a

  if (g_geometry_test) then

    ig = g_atog(ia)
    jg = g_atog(ja)

    xa = g_xloc(ig)
    ya = g_yloc(ig)
    za = g_zloc(ig)

    xb = g_xloc(jg)
    yb = g_yloc(jg)
    zb = g_zloc(jg)

    f= -1.0
    if (abs(xb-xa)>rm6) then
      f = (xf-xa)/(xb-xa)
    else if (abs(yb-ya)>rm6) then
      f = (yf-ya)/(yb-ya)
    else if (abs(zb-za)>rm6) then
      f = (zf-za)/(zb-za)
    endif

    xt = xa+f*(xb-xa)
    yt = ya+f*(yb-ya)
    zt = za+f*(zb-za)

    if (abs(xt-xf)>0.001) then
      print *,'x:xt,yt,zt,xf,yf,zf ', xt, yt, zt, xf, yf, zf
      stop
    endif

    if (abs(yt-yf)>0.001) then
      print *,'y:xt,yt,zt,xf,yf,zf ', xt, yt, zt, xf, yf, zf
      stop
    endif

    if (abs(zt-zf)>0.001) then
      print *,'z:xt,yt,zt,xf,yf,zf ', xt, yt, zt, xf, yf, zf
      stop
    endif

  endif

end subroutine AddConnection

! *************************************************************************** !

subroutine getFace( idir, positive, &
                    x000, x100, x010, x110, &
                    x001, x101, x011, x111, face )

  !
  ! Obtain the volume of a hexahederal grid cell
  !  with corners x000,..,x111
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt , intent(in) :: idir
  PetscBool, intent(in) :: positive
  PetscReal, intent(in) :: x000(3), x100(3), x010(3), x110(3), &
                           x001(3), x101(3), x011(3), x111(3)
  PetscReal,  intent(out) :: face(0:1, 0:1, 3)

  PetscInt :: id

  do id = 1, 3

    if (idir == g_xdir) then
      if (positive) then
        face(0, 0, id) = x100(id)
        face(1, 0, id) = x110(id)
        face(0, 1, id) = x101(id)
        face(1, 1, id) = x111(id)
      else
        face(0, 0, id) = x000(id)
        face(1, 0, id) = x010(id)
        face(0, 1, id) = x001(id)
        face(1, 1, id) = x011(id)
      endif
    endif

    if (idir == g_ydir) then
      if (positive) then
        face(0, 0, id) = x010(id)
        face(1, 0, id) = x110(id)
        face(0, 1, id) = x011(id)
        face(1, 1, id) = x111(id)
      else
        face(0, 0, id) = x000(id)
        face(1, 0, id) = x100(id)
        face(0, 1, id) = x001(id)
        face(1, 1, id) = x101(id)
      endif
    endif

    if (idir == g_zdir) then
      if (positive) then
        face(0, 0, id) = x001(id)
        face(1, 0, id) = x101(id)
        face(0, 1, id) = x011(id)
        face(1, 1, id) = x111(id)
      else
        face(0, 0, id) = x000(id)
        face(1, 0, id) = x100(id)
        face(0, 1, id) = x010(id)
        face(1, 1, id) = x110(id)
      endif
    endif
  enddo

end subroutine getFace

! *************************************************************************** !

function marea(fl, fu, idir)

  !
  ! Routine to get component id of the vector area of the quadrilateral f
  !
  ! Author: Dave Ponting
  ! Date: 02/05/18

  implicit none

  PetscReal marea
  PetscReal, intent(in) :: fl(0:1, 0:1, 3)
  PetscReal, intent(in) :: fu(0:1, 0:1, 3)
  PetscInt , intent(in) :: idir
  PetscReal :: xv(28), xvdd(28), x , y, diff, dy1, dy2

  PetscInt  :: jdir, kdir, nxv, nxvd, i, j, ierr, ifacepair, icell
  PetscReal :: alx, aly, aux, auy, &
               blx, bly, bux, buy, px, py, yll, yuu, x1, x2, dx, dy
  PetscBool :: qinter, qoverlap1, qoverlap2

  alx = 0.0
  aly = 0.0
  aux = 1.0
  auy = 1.0
  blx = 0.0
  bly = 1.0
  bux = 1.0
  buy = 0.0

  px = -1.0
  py = -1.0

  diff  = 1.0E-9
  marea = 0.0

  !  Find other directions

  call GetOtherDirections(idir, jdir, kdir)

  !  Set up y-range limits

  yll = fl(0, 0, kdir)
  yuu = fl(0, 0, kdir)

  !  Find all the distinct x-values and the y-limits

  nxv = 0
  do i = 0, 1
    do j = 0, 1

      x = fl(i, j, jdir)
      y = fl(i, j, kdir)
      nxv = nxv + 1
      xv(nxv) = x

      if (y>yuu) yuu = y
      if (y<yll) yll = y

      x = fu(i, j, jdir)
      y = fu(i, j, kdir)
      nxv = nxv + 1
      xv(nxv) = x

      if (y>yuu) yuu = y
      if (y<yll) yll = y

    enddo
  enddo

  !  Look for intersections between quads

  do i = 1, 4
    call GetQuadSegment(i, idir, fl, alx, aly, aux, auy)
    do j = 1, 4
      call GetQuadSegment(j, idir, fu, blx, bly, bux, buy)
      call GetIntersection(alx, aly, aux, auy, &
                           blx, bly, bux, buy, px, py, qinter)
      if (qinter) then
        nxv = nxv +1
        xv(nxv) = px
      endif
    enddo
  enddo

  !  Look for self-crossed quads

  do icell=1,2
    do ifacepair=1,2  ! Segment pairs (1 and 3) and (2 and 4)
      if (icell == 1) then
        call GetQuadSegment(ifacepair  , idir, fl, alx, aly, aux, auy)
        call GetQuadSegment(ifacepair+2, idir, fl, blx, bly, bux, buy)
      else
        call GetQuadSegment(ifacepair  , idir, fu, alx, aly, aux, auy)
        call GetQuadSegment(ifacepair+2, idir, fu, blx, bly, bux, buy)
      endif
      call GetIntersection(alx, aly, aux, auy, &
                           blx, bly, bux, buy, px, py, qinter)
      if (qinter) then
        nxv = nxv +1
        xv(nxv) = px
      endif
    enddo
  enddo

  !  Sort x-values into order

  ierr = 0
  call PetscSortReal(nxv, xv, ierr)

  !  Set up and count distinct depths

  nxvd = 1
  xvdd(1) = xv(1)

  do i = 2, nxv
    if (abs(xv(i)-xv(i-1))>diff) then
      nxvd = nxvd +1
      xvdd(nxvd) = xv(i)
    endif
  enddo

  !  There is no area unless at least two different x-values

  if (nxvd > 1) then

    marea = 0.0

    do i = 1, nxvd-1

      x1 = xvdd(i  )
      x2 = xvdd(i+1)
      dx = x2 - x1

      call FindYLimitsAtThisX(fl, fu, idir, x1, yll, yuu, dy1, qoverlap1)
      call FindYLimitsAtThisX(fl, fu, idir, x2, yll, yuu, dy2, qoverlap2)

      if (qoverlap1 .and. qoverlap2) then
        dy = 0.5*(dy1+dy2)
        if ( (dx > 0.0) .and. (dy > 0.0) ) then
          marea = marea + dx*dy
        endif
      endif

    enddo

  endif

end function marea

! *************************************************************************** !

subroutine FindYLimitsAtThisX(fl, fu, idir, x, yll, yuu, dy, qoverlap)

  !
  ! Given two quads fl and fu, find width of overlap region in y at x
  !
  ! Author: Dave Ponting
  ! Date: 02/05/18

  implicit none

  PetscReal, intent(in)  :: fl(0:1, 0:1, 3), fu(0:1, 0:1, 3)
  PetscInt , intent(in)  :: idir
  PetscReal, intent(in)  :: x, yll, yuu
  PetscReal, intent(out) :: dy
  PetscBool, intent(out) :: qoverlap

  PetscInt  :: j, ilu
  PetscReal :: alx, aly, aux, auy, blx, bly, bux, buy, px, py, &
               ylumin(2), ylumax(2), ymin, ymax
  PetscBool :: qfound(2), qinter

!  Set up the test line

  alx = x
  aux = x
  aly = yll - 10.0
  auy = yuu + 10.0

  !  For each quad, find the min and max intersections

  qfound = PETSC_FALSE

  do ilu = 1, 2
    do j = 1, 4
      if (ilu == 1) call GetQuadSegment(j, idir, fl, blx, bly, bux, buy)
      if (ilu == 2) call GetQuadSegment(j, idir, fu, blx, bly, bux, buy)
      call GetIntersection(alx, aly, aux, auy, &
                           blx, bly, bux, buy, px, py, qinter)
      if (qinter) then
        if (.not.qfound(ilu)) then
          ylumin(ilu) = py
          ylumax(ilu) = py
          qfound(ilu) = PETSC_TRUE
        else
          ylumin(ilu) = min(ylumin(ilu), py)
          ylumax(ilu) = max(ylumax(ilu), py)
        endif
      endif
    enddo
  enddo

  !  For an overlap, the min of the two maxes must be above the max of two mins

  ymin = max(ylumin(1), ylumin(2))
  ymax = min(ylumax(1), ylumax(2))

  dy = 0.0
  if (ymax>ymin) dy = ymax-ymin

  qoverlap = qfound(1) .and. qfound(2)

end subroutine FindYLimitsAtThisX

! *************************************************************************** !

subroutine InputReadPflotranStringNotComment(input, option)

  !
  ! Read a new input line, throwing away Eclipse comment lines (start with --)
  !
  ! Author: Dave Ponting
  ! Date: 02/11/18

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  character(len = MAXSTRINGLENGTH) :: word

  word = ' '

  do
    call InputReadPflotranString(input, option)
    if (input%ierr /= 0) exit                   ! Detect eof and leave
    word = adjustl(input%buf)
    if (word(1:2) /= '--') exit
  enddo

end subroutine InputReadPflotranStringNotComment

end module Grid_Grdecl_module
