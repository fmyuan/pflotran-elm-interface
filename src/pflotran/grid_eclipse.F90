module Grid_Eclipse_module

  !  Module to read a grid using industry-standard Eclipse keyword syntax
  !  and convert it to a Pflotran explicit unstructured grid

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Grid_Unstructured_Aux_module
  use Option_module
  use Input_Aux_module
  use Grid_Eclipse_Util_module

  implicit none

  character(len = MAXSTRINGLENGTH) :: g_word = ' '
  PetscBool, parameter :: g_callwaste = PETSC_FALSE

  ! Dimension of coordinate line arrays

  PetscInt, parameter :: g_mcl = 6

  ! Values to be used with kind() function

  PetscInt , parameter :: g_petsc_int  = 0
  PetscReal, parameter :: g_petsc_real = 0.0

  ! Identifiers for the num arrays

  PetscInt, public , parameter :: g_id_satnum_arr = 1
  PetscInt, public , parameter :: g_id_imbnum_arr = 2
  PetscInt, public , parameter :: g_id_eqlnum_arr = 3
  PetscInt, public , parameter :: g_id_mltnum_arr = 4
  PetscInt, public , parameter :: g_id_drpnum_arr(mdrpa) = [5 , 6, 7, 8, 9,10]
  PetscInt, public , parameter :: g_id_fipnum_arr(mfipa) = [11,12,13,14,15,16]

  !  Defaults for the num arrays

  PetscInt, public , parameter :: g_satnum_dflt = 1
  PetscInt, public , parameter :: g_imbnum_dflt = 1
  PetscInt, public , parameter :: g_eqlnum_dflt = 1
  PetscInt, public , parameter :: g_mltnum_dflt = 1
  PetscInt, public , parameter :: g_drpnum_dflt = 1
  PetscInt, public , parameter :: g_fipnum_dflt = 1

  private

  PetscInt :: e_fileunit = 0
  PetscBool :: e_formatted_r = PETSC_FALSE
  PetscInt, parameter :: e_typeB = 4
  PetscInt, parameter :: e_typeC = 5
  PetscInt, parameter :: e_typeD = 3
  PetscInt, parameter :: e_typeI = 1
  PetscInt, parameter :: e_typeS = 2
  PetscMPIInt, parameter, public :: tag_mpi_0 = 0

  ! ASCII codes for tab, line feed, carriage return and quotes

  PetscInt, parameter :: g_ictab  = 9
  PetscInt, parameter :: g_iclf   = 10
  PetscInt, parameter :: g_iccr   = 13
  PetscInt, parameter :: g_squote = 39
  PetscInt, parameter :: g_dquote = 34

  ! Face types for (FAULTS/MULTFLT)

  PetscInt, parameter :: g_facei  = 1
  PetscInt, parameter :: g_faceim = 2
  PetscInt, parameter :: g_facej  = 3
  PetscInt, parameter :: g_facejm = 4
  PetscInt, parameter :: g_facek  = 5
  PetscInt, parameter :: g_facekm = 6

  ! Direction and operation types for MULTREGT

  PetscInt, parameter  :: e_multregt_dir_xyz = 1
  PetscInt, parameter  :: e_multregt_op_all  = 1
  PetscInt             :: m_multregt = 0
  PetscInt             :: n_multregt = 0
  PetscInt ,allocatable:: g_multregt_ir1 (:)
  PetscInt ,allocatable:: g_multregt_ir2 (:)
  PetscReal,allocatable:: g_multregt_mult(:)

  PetscInt :: g_nlgr = 1

  !--Integer pointer-container-------------------------------------------------

  type, public :: int_ptr_type
    PetscInt, pointer :: p(:) => null()
    PetscBool         :: pset =  PETSC_FALSE
  end type int_ptr_type

  type, public :: lgr_type

    character(len = MAXSTRINGLENGTH), public :: ref_name = 'Field'

    ! Problem dimensions (set defaults for single cell)

    PetscInt :: nx     = 1
    PetscInt :: ny     = 1
    PetscInt :: nz     = 1
    PetscInt :: nxp    = 2
    PetscInt :: nyp    = 2
    PetscInt :: nxy    = 1
    PetscInt :: nxpnyp = 4
    PetscInt :: nxyz   = 1

    PetscInt :: ixl = 1
    PetscInt :: ixu = 1

    PetscInt :: iyl = 1
    PetscInt :: iyu = 1

    PetscInt :: izl = 1
    PetscInt :: izu = 1

    PetscInt :: nxg = 1
    PetscInt :: nyg = 1
    PetscInt :: nzg = 1

  ! Grid-sized arrays

    PetscReal, pointer :: coord(:) => null()
    PetscReal, pointer :: zcorn(:) => null()
    PetscInt , pointer :: iglob(:) => null()

    PetscInt , pointer :: nlpgx(:) => null()
    PetscInt , pointer :: nlpgy(:) => null()
    PetscInt , pointer :: nlpgz(:) => null()

    PetscInt , pointer :: ibpgx(:) => null()
    PetscInt , pointer :: ibpgy(:) => null()
    PetscInt , pointer :: ibpgz(:) => null()

    PetscReal, pointer :: hrefx(:) => null()
    PetscReal, pointer :: hrefy(:) => null()
    PetscReal, pointer :: hrefz(:) => null()

    PetscReal, pointer :: dx  (:) => null()
    PetscReal, pointer :: dy  (:) => null()
    PetscReal, pointer :: dz  (:) => null()

    PetscReal, pointer :: kx  (:) => null()
    PetscReal, pointer :: ky  (:) => null()
    PetscReal, pointer :: kz  (:) => null()

    PetscReal, pointer :: mv  (:) => null()

    PetscReal, pointer :: mx  (:) => null()
    PetscReal, pointer :: my  (:) => null()
    PetscReal, pointer :: mz  (:) => null()

    PetscReal, pointer :: mxn (:) => null()
    PetscReal, pointer :: myn (:) => null()
    PetscReal, pointer :: mzn (:) => null()

    PetscReal, pointer :: tx  (:) => null()
    PetscReal, pointer :: ty  (:) => null()
    PetscReal, pointer :: tz  (:) => null()

    PetscReal, pointer :: tops(:) => null()
    PetscReal, pointer :: poro(:) => null()
    PetscReal, pointer :: ntg (:) => null()

    PetscReal, pointer :: xloc(:) => null()
    PetscReal, pointer :: yloc(:) => null()
    PetscReal, pointer :: zloc(:) => null()

    PetscInt , pointer :: actn(:) => null()

    PetscInt , pointer :: satn(:) => null()
    PetscInt , pointer :: imbn(:) => null()
    PetscInt , pointer :: eqln(:) => null()
    PetscInt , pointer :: mltn(:) => null()
    PetscInt , pointer :: tbcn(:) => null()
    PetscInt , pointer :: prcn(:) => null()
    type(int_ptr_type) :: drpn(mdrpa)
    type(int_ptr_type) :: fipn(mfipa)

    PetscInt , pointer :: gtoa(:) => null()
    PetscInt , pointer :: ihost(:) => null()

    PetscBool :: cpgallocated  = PETSC_FALSE
    PetscBool :: hdset(g_ndir) = PETSC_FALSE

    PetscBool :: drpallocated  = PETSC_FALSE
    PetscBool :: fipallocated(mfipa) = PETSC_FALSE

  end type lgr_type

  type(lgr_type) ::  g_g(g_mlgr)

  PetscInt :: g_grid_iof(g_mlgr) = 0
  PetscInt :: g_cnat_iof(g_mlgr) = 0

  ! Total count

  PetscInt :: g_nxyzt  = 0

  ! Flag indicating SATNUM table indices read and max. SATNUM value

  PetscInt :: g_rsatn  = 0
  PetscInt :: g_maxsatn= 0

  ! Flag indicating IMBNUM table indices read and max. IMBNUM value

  PetscInt :: g_rimbn  = 0
  PetscInt :: g_maximbn= 0

  ! Flag indicating EQLNUM table indices read and max. EQLNUM value

  PetscInt :: g_reqln  = 0
  PetscInt :: g_maxeqln= 0

  ! Flag indicating MULTNUM table indices read and max. MULTNUM value

  PetscInt :: g_rmltn  = 0
  PetscInt :: g_maxmltn= 0

  ! Markers for procnum

  PetscInt :: g_rprcn  = 0

  ! Flag indicating DRP table indices read and max. DRP value

  PetscInt :: g_rdrpn  (mdrpa) = 0
  PetscInt :: g_maxdrpn(mdrpa) = 0

  ! Flag indicating FIP indices read and max. FIP value

  PetscInt :: g_rfipn  (mfipa) = 0
  PetscInt :: g_maxfipn(mfipa) = 0

  ! Flag indicating DIMENS read (contains problem dimensions)

  PetscBool :: g_dimens_read = PETSC_FALSE

  ! z-flip from Eclipse to Pflotran convention and corner-point flags

  PetscBool :: g_iscpg = PETSC_FALSE
  PetscReal :: z_flip  = -1.0

  PetscBool :: g_isnewtran    = PETSC_FALSE
  PetscBool :: g_isdpcf       = PETSC_FALSE
  PetscBool :: g_create_tc    = PETSC_FALSE

  ! Counters for active cells, connections, max connections, faults, pinchouts

  PetscInt :: g_na  = 1
  PetscInt :: g_nc  = 0
  PetscInt :: g_mc  = 1
  PetscInt :: g_nf  = 0
  PetscInt :: g_npo = 0
  PetscInt :: g_nac = 0

  ! Active sized arrays

  PetscInt , pointer :: g_atog(:) => null()
  PetscInt , pointer :: g_atol(:) => null()
  PetscInt , pointer :: g_atoc(:) => null()

  PetscInt , pointer :: g_atox(:) => null()
  PetscInt , pointer :: g_atoy(:) => null()
  PetscInt , pointer :: g_atoz(:) => null()

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

  PetscReal, pointer :: g_tran   (:) => null()   ! Transmisibility
  PetscReal, pointer :: g_tran_th(:) => null()   ! Thermal transmisibility

  ! Direction (g_xdir, g_ydir, g_zdir for upper connection (lower to upper)
  !            -g_xdir,-g_ydir,-g_zdir for lower connection (upper to lower)
  !            0 if not known
  PetscInt , pointer :: g_idir   (:) => null()

  PetscReal, allocatable :: g_ttx(:)
  PetscReal, allocatable :: g_tty(:)
  PetscReal, allocatable :: g_ttz(:)
  PetscBool              :: g_tt_allocated = PETSC_FALSE

  PetscReal :: g_mapaxes(6) = 0.0
  PetscReal :: g_minpv(1)   = 1.0D-3
  PetscReal :: g_pinch(1)   = 1.0D-3
  PetscReal :: g_dpcf (2)   = 0.0
  PetscReal :: g_pvfloor(1) = 0.0

  PetscBool :: g_atog_held = PETSC_FALSE

  ! Column pointer for reader

  PetscInt :: g_column = 1

  character(len = MAXSTRINGLENGTH) :: g_error_string = 'OK'
  PetscInt :: g_error_flag = 0

!  Null test flag

  PetscBool :: g_geometry_test = PETSC_FALSE

! GRIDFILE option (defaulted to 1)

  PetscInt :: g_gridfileoption=1

! Arrays to hold init array limits

  PetscInt, parameter :: g_mina = 21

  public  :: GridEclipseRead, &
             GridEclipseGetNumArrSet, &
             GridEclipseGetPorPerm, &
             GridEcilpseGetProcnumFlag, &
             GridEclipsePorPermExchangeAndSet, &
             GridEclipseProcnumExchangeAndSet, &
             GridEcipseDeallocatePorPermArrays

  private :: GrdeclReader

contains

! *************************************************************************** !

subroutine GridEclipseRead(unstructured_grid, filename, option)
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
  PetscInt :: fileid
  PetscErrorCode :: ierr

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

  if (option%myrank == option%comm%io_rank) then
    call GrdeclReader(input, option)
  endif

  call BroadcastInt(g_error_flag, option)
  if (g_error_flag>0) then
    input%ierr = INPUT_ERROR_DEFAULT
    call MPI_Bcast(g_error_string, MAXSTRINGLENGTH, MPI_CHARACTER, &
                    option%comm%io_rank, option%mycomm, ierr)
    call InputErrorMsg(input, option, 'GRDECL file', g_error_string)
  endif

  ! Reset default blocking error state

  call OptionSetBlocking(option, PETSC_TRUE)
  call OptionCheckNonBlockingError(option)

  ! Destroy the input system

  call InputDestroy(input)

  ! Sycn after read

  call MPI_Barrier(option%mycomm, ierr)

  ! Distribute and store the cell data

  call DistributeCells(explicit_grid, option)

  ! Distribute and store the connection data

  call DistributeConnections(explicit_grid, option)

  ! Distribute and store cell counts and num flags

  call DistributeCounts(option)

  ! Set up the cell locations

  if (option%myrank == option%comm%io_rank) then
    call CreateElements(unstructured_grid, explicit_grid)
  endif

end subroutine GridEclipseRead

! *************************************************************************** !

subroutine ACEMOp(zto, zfr, ixl, ixu, &
                   iyl, iyu, &
                   izl, izu, isadd, iscpy, iseql, ismlt, &
                   zk, qerr, convdist, convtran, ilgr, option)

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
  PetscReal, intent(in)  :: convdist, convtran
  PetscInt , intent(in)  :: ilgr
  type(option_type) :: option

  PetscBool :: isint, isreal, isintb, isrealb, &
               is_mult_only , is_mult_dummy, &
               is_depth     , is_depth_dummy, &
               is_perm      , is_perm_dummy, &
               is_dist      , is_dist_dummy, &
               is_tran      , is_tran_dummy

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

  is_dist        = PETSC_FALSE
  is_dist_dummy  = PETSC_FALSE

  is_tran        = PETSC_FALSE
  is_tran_dummy  = PETSC_FALSE

  nullify(ai, ar, bi, br)

  v  = 0.0
  iv = 0

  qerr = PETSC_FALSE

!  Get the required arrays

  call GetGridArrayPointer(zto, isint, isreal, ai, ar, qerr, iseql, &
                            is_mult_only, is_depth, is_perm, is_dist, &
                            is_tran, ilgr, option)
  if (iscpy) then
    call GetGridArrayPointer(zfr, isintb, isrealb, bi, br, qerr, iseql, &
                              is_mult_dummy, is_depth_dummy, &
                              is_perm_dummy, is_dist_dummy, is_tran_dummy, &
                              ilgr, option)
  else
    call ProcessArgToReal(v, zfr, zk, qerr)
  endif

!  Check if operation is OK

  if (is_mult_only .and. (.not. ismlt)) then
    call SetError('Only MULT supported for ' // trim(zto), qerr)
  else

!  Carry out operations

    if (.not.qerr) then

!     Case of distance array (DX,DY,DZ,TOPS): unit conversion (ft to m if req)

      if (is_dist) then
        if (isadd .or. iseql) v = v*convdist
      endif

!     Case of depth array like TOPS: convert set-or-add value to an elevation

      if (is_depth) then
        if (isadd .or. iseql) v = -v
      endif

!     Case of perm array like PERMX: convert set-or-add value from mD to m^2

      if (is_perm) then
        conv = GetMDtoM2Conv()
        if (isadd .or. iseql) v = v*conv
      endif

!     Case of tran array like TRANX: convert set-or-add value
!     from cp.rm3/day/bar to PaS.rm3/sec/Pa

      if (is_tran) then
        conv = GetETtoTrConv()*convtran
        if (isadd .or. iseql) v = v*conv
      endif

      do iz = izl, izu
        do iy = iyl, iyu
          do ix = ixl, ixu

            ig = GetNaturalIndex(ix, iy, iz, ilgr)

!  Case of integer array
            if (isint) then
              iv = nint(v,kind(g_petsc_int))
              if (isadd) ai(ig) = ai(ig) + iv
              if (iscpy) then
                if (isintb) bi(ig) =      ai(ig)
                if (isrealb) bi(ig) = nint(ar(ig),kind(g_petsc_int))
              endif
              if (iseql) ai(ig) = iv
              if (ismlt) ai(ig) = ai(ig) * iv
            endif

! Case of real array
            if (isreal) then
              if (isadd) ar(ig) = ar(ig) + v
              if (iscpy) then
                if (isintb) br(ig) = real(ai(ig),kind(g_petsc_real))
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

end subroutine ACEMOp

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

subroutine AllocateGridArrays(nx,ny,nz,ilgr)
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

  PetscInt, intent(in) ::nx, ny, nz, ilgr

  PetscInt  :: nxyz, def_an, nxg, nyg, nzg, ndlis
  PetscReal :: def_pp, def_mn
  PetscBool :: qerr

  PetscBool, parameter :: qisz_true  = PETSC_TRUE
  PetscBool, parameter :: qisz_false = PETSC_FALSE

  ! Basic dimensions

  g_g(ilgr)%nx = nx
  g_g(ilgr)%ny = ny
  g_g(ilgr)%nz = nz

  g_g(ilgr)%nxp    = nx+1
  g_g(ilgr)%nyp    = ny+1
  g_g(ilgr)%nxy    = nx*ny
  g_g(ilgr)%nxpnyp = (nx+1)*(ny+1)
  g_g(ilgr)%nxyz   = nx*ny*nz

  nxyz = g_g(ilgr)%nxyz

  ! Set standard defaults

  def_pp = 0.0D0 ! Default perm and poro
  def_mn = 1.0D0 ! Default mults and ntg
  def_an = 1     ! Default actnum

  ! LGR defaults (to be filled in later)

  if(ilgr>1) then
    def_pp = UNINITIALIZED_DOUBLE
    def_mn = UNINITIALIZED_DOUBLE
    def_an = -1

    nxg = g_g(ilgr)%nxg
    nyg = g_g(ilgr)%nyg
    nzg = g_g(ilgr)%nzg

    allocate(g_g(ilgr)%nlpgx(nxg))
    allocate(g_g(ilgr)%nlpgy(nyg))
    allocate(g_g(ilgr)%nlpgz(nzg))

    allocate(g_g(ilgr)%ibpgx(nxg))
    allocate(g_g(ilgr)%ibpgy(nyg))
    allocate(g_g(ilgr)%ibpgz(nzg))

    allocate(g_g(ilgr)%hrefx(nx))
    allocate(g_g(ilgr)%hrefy(ny))
    allocate(g_g(ilgr)%hrefz(nz))

    g_g(ilgr)%hrefx=1.0d0
    g_g(ilgr)%hrefy=1.0d0
    g_g(ilgr)%hrefz=1.0d0

    ! Set up default local cells/global cell

    call SetDefaultLocalsPerGlobal(nx,nxg,g_g(ilgr)%nlpgx,qisz_false)
    call SetDefaultLocalsPerGlobal(ny,nyg,g_g(ilgr)%nlpgy,qisz_false)
    call SetDefaultLocalsPerGlobal(nz,nzg,g_g(ilgr)%nlpgz,qisz_true)

   ! Make base pointers from local cells/global cell

    call MakeLocalBasePointers(g_g(ilgr)%nlpgx,g_g(ilgr)%ibpgx,nx,nxg,qerr, &
                               ndlis)
    call MakeLocalBasePointers(g_g(ilgr)%nlpgy,g_g(ilgr)%ibpgy,ny,nyg,qerr, &
                               ndlis)
    call MakeLocalBasePointers(g_g(ilgr)%nlpgz,g_g(ilgr)%ibpgz,nz,nzg,qerr, &
                               ndlis)

  endif

  ! Allocate and set

  allocate(g_g(ilgr)%dx(nxyz))
  allocate(g_g(ilgr)%dy(nxyz))
  allocate(g_g(ilgr)%dz(nxyz))

  allocate(g_g(ilgr)%kx(nxyz))
  allocate(g_g(ilgr)%ky(nxyz))
  allocate(g_g(ilgr)%kz(nxyz))

  allocate(g_g(ilgr)%mv(nxyz))

  allocate(g_g(ilgr)%mx(nxyz))
  allocate(g_g(ilgr)%my(nxyz))
  allocate(g_g(ilgr)%mz(nxyz))

  allocate(g_g(ilgr)%mxn(nxyz))
  allocate(g_g(ilgr)%myn(nxyz))
  allocate(g_g(ilgr)%mzn(nxyz))

  allocate(g_g(ilgr)%tx(nxyz))
  allocate(g_g(ilgr)%ty(nxyz))
  allocate(g_g(ilgr)%tz(nxyz))

  allocate(g_g(ilgr)%tops(nxyz))
  allocate(g_g(ilgr)%poro(nxyz))
  allocate(g_g(ilgr)%ntg (nxyz))

  allocate(g_g(ilgr)%xloc(nxyz))
  allocate(g_g(ilgr)%yloc(nxyz))
  allocate(g_g(ilgr)%zloc(nxyz))

  ! Int arrays set to del_an=1 for global, -1 for locals,
  ! so that later can copy over global values to unset locals

  allocate(g_g(ilgr)%actn(nxyz))
  allocate(g_g(ilgr)%satn(nxyz))
  allocate(g_g(ilgr)%imbn(nxyz))
  allocate(g_g(ilgr)%eqln(nxyz))
  allocate(g_g(ilgr)%mltn(nxyz))

  if(g_rprcn==1) then
    allocate(g_g(ilgr)%prcn(nxyz))
    g_g(ilgr)%prcn = 0    ! Processor numbers
  endif

  if(ilgr>1) then
    allocate(g_g(ilgr)%ihost(nxyz))
    g_g(ilgr)%ihost =  0
  endif

  allocate(g_g(ilgr)%gtoa(nxyz))


  g_g(ilgr)%dx  =  0.d0
  g_g(ilgr)%dy  =  0.d0
  g_g(ilgr)%dz  =  0.d0
  g_g(ilgr)%kx  =  def_pp
  g_g(ilgr)%ky  =  def_pp
  g_g(ilgr)%kz  =  def_pp
  g_g(ilgr)%mv  =  def_mn
  g_g(ilgr)%mx  =  def_mn
  g_g(ilgr)%my  =  def_mn
  g_g(ilgr)%mz  =  def_mn
  g_g(ilgr)%mxn =  def_mn
  g_g(ilgr)%myn =  def_mn
  g_g(ilgr)%mzn =  def_mn
  g_g(ilgr)%tx  = -1.d0
  g_g(ilgr)%ty  = -1.d0
  g_g(ilgr)%tz  = -1.d0
  g_g(ilgr)%tops =  UNINITIALIZED_DOUBLE ! To trap unset layers
  g_g(ilgr)%poro =  def_pp
  g_g(ilgr)%ntg  =  def_mn ! Net to gross ratio (default to 1)
  g_g(ilgr)%xloc = 0.d0
  g_g(ilgr)%yloc = 0.d0
  g_g(ilgr)%zloc = 0.d0
  g_g(ilgr)%actn = def_an ! Actnum (used to set cells inactive)
  g_g(ilgr)%satn = def_an ! Saturation table numbers
  g_g(ilgr)%imbn = def_an ! Imbibition table numbers
  g_g(ilgr)%eqln = def_an ! Equilibration table numbers
  g_g(ilgr)%mltn = def_an ! Multnum numbers
  g_g(ilgr)%gtoa = -1

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

  allocate(g_bvol(g_na))

  allocate(g_x(g_na))
  allocate(g_y(g_na))
  allocate(g_z(g_na))

  allocate(g_cia(g_mc))
  allocate(g_cja(g_mc))

  allocate(g_ccx(g_mc))
  allocate(g_ccy(g_mc))
  allocate(g_ccz(g_mc))
  allocate(g_carea(g_mc))

  allocate(g_tran(g_mc))
  allocate(g_tran_th(g_mc))
  allocate(g_idir(g_mc))

  g_bvol = 0.d0
  g_x = 0.d0
  g_y = 0.d0
  g_z = 0.d0
  g_cia = -1
  g_cja = -1
  g_ccx = 0.d0
  g_ccy = 0.d0
  g_ccz = 0.d0
  g_carea = 0.d0
  g_tran = 0.d0
  g_tran_th = 0.d0
  g_idir = 0.d0

end subroutine AllocateActiveArrays

! *************************************************************************** !

subroutine CheckAllocCPGA(ilgr)

  !
  ! Check if the ccord and zcorn arrays allocated and allocate if not
  !
  ! Author: Dave Ponting
  ! Date: 08/04/18

  implicit none

  PetscInt, intent(in) :: ilgr

  if (.not.g_g(ilgr)%cpgallocated) then
    allocate(g_g(ilgr)%coord(6*g_g(ilgr)%nxpnyp))
    allocate(g_g(ilgr)%zcorn(8*g_g(ilgr)%nxyz))
    allocate(g_g(ilgr)%iglob(g_g(ilgr)%nxyz))
    g_g(ilgr)%coord = 0.d0
    g_g(ilgr)%zcorn = 0.d0
    g_g(ilgr)%iglob = 0
    g_g(ilgr)%cpgallocated = PETSC_TRUE
  endif

end subroutine CheckAllocCPGA

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

  if (.not.InputError(ierr)) then
    qerr = PETSC_FALSE
  else
    call SetError(zerr, qerr)
  endif

end subroutine CheckError

! *************************************************************************** !

subroutine CheckLGRLocation(ixl, ixu, iyl, iyu, izle, izue, zilgr, qerr)
  !
  ! Check that a local grid does not overlap or adjoin another
  !
  ! Author: Dave Ponting
  ! Date: 09/21/21
  !

  use String_module,only : StringCompareIgnoreCase

  implicit none

  PetscInt , intent(in)  :: ixl, ixu, iyl, iyu, izle, izue
  character(len=*), intent(in) :: zilgr
  PetscBool, intent(out) :: qerr

  PetscInt :: nx, ny, nz
  PetscInt :: jlgr, jxl, jxu, jyl, jyu, jzle, jzue
  PetscInt :: ix, iy, ize
  character(len=MAXSTRINGLENGTH) :: zjlgr

!--Initialise-----------------------------------------------------------------

  qerr = PETSC_FALSE
  nx   = g_g(g_ifld)%nx
  ny   = g_g(g_ifld)%ny
  nz   = g_g(g_ifld)%nz

  do jlgr = 2, g_nlgr

    zjlgr = g_g(jlgr)%ref_name

    jxl   = g_g(jlgr)%ixl
    jxu   = g_g(jlgr)%ixu

    jyl   = g_g(jlgr)%iyl
    jyu   = g_g(jlgr)%iyu

    jzle  = nz + 1 - g_g(jlgr)%izu
    jzue  = nz + 1 - g_g(jlgr)%izl

    ! Check if name already used

    if(StringCompareIgnoreCase(zilgr,zjlgr)) then
      call SetError('Duplicate LGR names: '//trim(zilgr)// &
                    ' and '//trim(zjlgr),qerr)
      return
    endif

    ! Check if any of the cells in the new LGR are in or adjoint existing lgr

    do ix = ixl, ixu
      do iy = iyl, iyu
        do ize = izle, izue
          qerr = GetInOrAdjBox(ix,iy,ize, &
                               jxl,jxu, jyl,jyu, jzle,jzue, zilgr,zjlgr)
          if(qerr) return
        enddo
      enddo
    enddo

  enddo

end subroutine CheckLGRLocation

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

subroutine CheckForDimensionErr(zdirn,nd,ndfile,zmess,qerr)

  !
  ! Set up a dimension error message
  !
  ! Author: Dave Ponting
  ! Date: 01/15/19

  implicit none

  character(len=*),intent(in) :: zdirn
  PetscInt,intent(in) :: nd, ndfile
  character(len = MAXSTRINGLENGTH),intent(out) :: zmess
  PetscBool,intent(out) :: qerr

  character(len = MAXWORDLENGTH) :: zdfile, zd

  if(ndfile/=nd) then
    write(zd    ,'(I6)') nd
    write(zdfile,'(I6)') ndfile
    qerr  = PETSC_TRUE
    zmess = trim(zdirn)//'-dimension file mismatch, expect '//trim(zd)//', file has '//trim(zdfile)
  endif

end subroutine CheckForDimensionErr

! *************************************************************************** !

subroutine CheckValidOperateArray(zname,qerr)

  !
  ! Check if an array is valid with OPERATE
  !
  ! Author: Dave Ponting
  ! Date: 09/21/21
  !

  use String_module

  implicit none

  character(len = *), intent(in) :: zname
  PetscBool         , intent(out) :: qerr

  PetscBool :: qok

!--Initialise------------------------------------------------------------------

  qok  = PETSC_FALSE
  qerr = PETSC_FALSE

  if (StringCompareIgnoreCase(zname,'permx')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'permy')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'permz')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'poro')) qok=PETSC_TRUE

  if (StringCompareIgnoreCase(zname,'multx')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'multy')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'multz')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'multv')) qok=PETSC_TRUE

  if (StringCompareIgnoreCase(zname,'multx-')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'multy-')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'multz-')) qok=PETSC_TRUE

  if (StringCompareIgnoreCase(zname,'actnum')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'satnum')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'imbnum')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'eqlnum')) qok=PETSC_TRUE
  if (StringCompareIgnoreCase(zname,'multnum')) qok=PETSC_TRUE

  if (.not.qok) then
    call SetError('Array ' // trim(zname) //' not available with OPERATE', qerr)
    qerr = PETSC_TRUE
  endif

end subroutine CheckValidOperateArray

! *************************************************************************** !

subroutine CopyFromBufferI(a, ibuf, il, iu)
  !
  ! Copy values from a int*32 buffer
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscInt, intent(out) :: a(:)
  PetscInt, intent(in) :: ibuf(:)
  PetscInt, intent(in) :: il, iu

  PetscInt :: i

  do i = il, iu
    a(i) = ibuf(i-il+1)
  enddo

end subroutine CopyFromBufferI

! *************************************************************************** !

subroutine CopyFromBufferS(a, sbuf, il, iu)
  !
  ! Copy values from a real*32 buffer
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscReal, intent(out) :: a(:)
  PetscReal, intent(in) :: sbuf(:)
  PetscInt, intent(in) :: il, iu

  PetscInt :: i

  do i = il, iu
    a(i) = sbuf(i-il+1)
  enddo

end subroutine CopyFromBufferS

! *************************************************************************** !

subroutine CopyFromBufferD(a, dbuf, il, iu)
  !
  ! Copy values from a real*64 buffer
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscReal, intent(out) :: a(:)
  PetscReal, intent(in) :: dbuf(:)
  PetscInt, intent(in) :: il, iu

  PetscInt :: i

  do i = il, iu
    a(i) = dbuf(i-il+1)
  enddo

end subroutine CopyFromBufferD

! *************************************************************************** !

subroutine CopyFromBufferB(a, ibuf, il, iu)
  !
  ! Copy bool from a int*32 buffer (false->0;true->1)
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscBool, intent(out) :: a(:)
  PetscInt, intent(in) :: ibuf(:)
  PetscInt, intent(in) :: il, iu

  PetscInt :: i

  do i = il, iu
   if (ibuf(i-il+1) == 1) then
     a(i)=PETSC_TRUE
   else
     a(i)=PETSC_FALSE
   endif
  enddo

end subroutine CopyFromBufferB

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
              ix, iy, iz, iox, ioy, ioz, ilgr

  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3), xw(3)

  xw = 0.0

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

    ig   = g_atog(ia)
    ilgr = g_atol(ia)
    call GetCellCoordinates(ig, ix, iy, iz, ilgr)

    ! Get the corners

    call GetCorners(ix, iy, iz, &
                     x000, x100, x010, x110, &
                     x001, x101, x011, x111, &
                     g_g(ilgr)%coord, g_g(ilgr)%zcorn, g_g(ilgr)%nx, g_g(ilgr)%ny)

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

subroutine CreateTransT(ix,iy,iz,ig,ia,idir,treq,kdir,ilgr,option)
  !
  ! Create a specified transmissibility
  !
  ! Author: Dave Ponting
  ! Date: 11/24/20
  !

  implicit none

  PetscInt,intent(in)   :: ix,iy,iz,ig,ia,idir,ilgr
  PetscReal, intent(in) :: treq
  PetscReal, pointer    :: kdir(:)
  type(option_type)     :: option

  PetscInt  :: jx,jy,jz,jg,ja,jdir,jlgr
  PetscReal :: tcur, area, di, dj, &
               xi(g_ndir), xj(g_ndir), xa(g_ndir), ti, tj, &
               dxi, dyi, dzi, dxj, dyj, dzj, &
               xf, yf, zf, area_th
  PetscBool :: ingrid

  ! Find the neighbour

  ingrid = GetPosNeighbour(ix,iy,iz,idir,jx,jy,jz,ilgr)

  ! If in grid, check if active

  if (ingrid) then

    jg   = GetNaturalIndex(jx, jy, jz, ilgr)
    jlgr = ilgr

    ! Check for active neighbour

    ja = g_g(ilgr)%gtoa(jg)
    if (ja > 0) then

      ! Find cell dimensions

      dxi  = g_g(ilgr)%dx(ig)
      dyi  = g_g(ilgr)%dy(ig)
      dzi  = g_g(ilgr)%dz(ig)

      dxj  = g_g(ilgr)%dx(jg)
      dyj  = g_g(ilgr)%dy(jg)
      dzj  = g_g(ilgr)%dz(jg)

      ! Find reasonable estimate of area

      area = 0.0
      if (idir == g_xdir) area = 0.5*(dyi*dzi + dyj*dzj)
      if (idir == g_ydir) area = 0.5*(dzi*dxi + dzj*dxj)
      if (idir == g_zdir) area = 0.5*(dxi*dyi + dxj*dyj)

      ! Cell locations

      call SetVec3(xi,g_g(ilgr)%xloc(ig),g_g(ilgr)%yloc(ig),g_g(ilgr)%zloc(ig))
      call SetVec3(xj,g_g(ilgr)%xloc(jg),g_g(ilgr)%yloc(jg),g_g(ilgr)%zloc(jg))

      ! Place the interface at half-way

      do jdir = 1, g_ndir
        xa(jdir) = 0.5*(xi(jdir)+xj(jdir))
      enddo

      xf = xa(g_xdir)
      yf = xa(g_ydir)
      zf = xa(g_zdir)

      ! Estimate current transmissibility

      di = DistVec3(xi,xa)
      dj = DistVec3(xj,xa)

      if ((di>0.0) .and. (dj>0.0)) then

        ti = kdir(ig)*area/di
        tj = kdir(jg)*area/dj

        tcur = 0.0
        if ((ti>0.0) .and. (tj>0.0)) then
          tcur = 1.0/(1.0/ti+1.0/tj)
        endif

        ! If tcur exists (really if perms>0) multiply area to match required and store
        if (tcur>0.0) then
          area_th = area
          area    = area*treq/tcur
          call AddConnection(ia,ja,idir,g_ipon_u,xf,yf,zf,area,ilgr,jlgr,treq,area_th)
        endif

      endif

    endif

  endif

end subroutine CreateTransT

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

  if (option%myrank == option%comm%io_rank) then

    allocate(temp_real_array(5, nals+1))

    iabase = 0
    do irank = 0, option%comm%size-1

      narank = nals
      if (irank < rem) narank = narank + 1

      if (irank == option%comm%io_rank) then

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
          temp_real_array(1, ial) = real  (ia,kind(g_petsc_real))
          temp_real_array(2, ial) = g_x   (ia)
          temp_real_array(3, ial) = g_y   (ia)
          temp_real_array(4, ial) = g_z   (ia)
          temp_real_array(5, ial) = g_bvol(ia)
        enddo
        call MPI_Send(temp_real_array, int_mpi, MPI_DOUBLE_PRECISION, irank, &
                      tag_mpi_0, option%mycomm, ierr)
      endif
      iabase = iabase + narank
    enddo
  else
    ! other ranks post the recv
    allocate(temp_real_array(5, nal))
    int_mpi = nal*5
    call MPI_Recv(temp_real_array, int_mpi, &
                  MPI_DOUBLE_PRECISION, option%comm%io_rank, &
                  MPI_ANY_TAG, option%mycomm, status_mpi, ierr)
    do ia = 1, nal
      explicit_grid%cell_ids      (ia)   = nint(temp_real_array(1, ia), &
                                                kind(g_petsc_int))
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
  if (option%myrank == option%comm%io_rank) then
    allocate(temp_real_array(9, ncls+1))
    ! read for other processors
    icbase = 0
    do irank = 0, option%comm%size-1
      temp_real_array = UNINITIALIZED_DOUBLE
      ncrank = ncls
      if (irank < rem) ncrank = ncrank + 1
      ! if the cells reside on io_rank
      if (irank == option%comm%io_rank) then

        do icl = 1, ncl
          ic = icl + icbase
          explicit_grid%connections   (1, icl)   = g_cia  (ic)
          explicit_grid%connections   (2, icl)   = g_cja  (ic)
          explicit_grid%face_centroids(icl)%x = g_ccx  (ic)
          explicit_grid%face_centroids(icl)%y = g_ccy  (ic)
          explicit_grid%face_centroids(icl)%z = g_ccz  (ic)
          explicit_grid%face_areas    (icl)   = g_carea(ic)
        enddo
      else
        ! otherwise communicate to other ranks
        int_mpi = ncrank*9
        do icl = 1, ncrank
          ic = icl + icbase
          temp_real_array(1, icl) = real(g_cia  (ic),kind(g_petsc_real))
          temp_real_array(2, icl) = real(g_cja  (ic),kind(g_petsc_real))
          temp_real_array(3, icl) =      g_ccx  (ic)
          temp_real_array(4, icl) =      g_ccy  (ic)
          temp_real_array(5, icl) =      g_ccz  (ic)
          temp_real_array(6, icl) =      g_carea(ic)
          temp_real_array(7, icl) =      g_tran   (ic)
          temp_real_array(8, icl) =      g_tran_th(ic)
          temp_real_array(9, icl) =      g_idir   (ic)
        enddo
        call MPI_Send(temp_real_array, int_mpi, MPI_DOUBLE_PRECISION, irank, &
                      tag_mpi_0, option%mycomm, ierr)
      endif
      icbase = icbase + ncrank
    enddo
  else
    ! other ranks post the recv
    allocate(temp_real_array(9, ncl))
    int_mpi = ncl*9
    call MPI_Recv(temp_real_array, int_mpi, &
                  MPI_DOUBLE_PRECISION, option%comm%io_rank, &
                  MPI_ANY_TAG, option%mycomm, status_mpi, ierr)
    do icl = 1, ncl
      explicit_grid%connections   (1, icl)   = nint(temp_real_array(1, icl), &
                                                    kind(g_petsc_int))
      explicit_grid%connections   (2, icl)   = nint(temp_real_array(2, icl), &
                                                    kind(g_petsc_int))
      explicit_grid%face_centroids(icl)%x =      temp_real_array(3, icl)
      explicit_grid%face_centroids(icl)%y =      temp_real_array(4, icl)
      explicit_grid%face_centroids(icl)%z =      temp_real_array(5, icl)
      explicit_grid%face_areas    (icl)   =      temp_real_array(6, icl)
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

  use Grid_Eclipse_Util_module, only : SetProbDims

  implicit none

  type(option_type) :: option

  PetscInt :: ilgr

  character(len = MAXSTRINGLENGTH) :: zarr(1)

  zarr(1) = 'Field'

  call BroadcastInt (g_na      , option)
  call BroadcastInt (g_rsatn   , option)
  call BroadcastInt (g_maxsatn , option)
  call BroadcastInt (g_rimbn   , option)
  call BroadcastInt (g_maximbn , option)
  call BroadcastIntN(g_rdrpn   , option)
  call BroadcastIntN(g_maxdrpn , option)
  call BroadcastIntN(g_rfipn   , option)
  call BroadcastIntN(g_maxfipn , option)
  call BroadcastInt (g_reqln   , option)
  call BroadcastInt (g_maxeqln , option)
  call BroadcastInt (g_rprcn   , option)

  call BroadcastInt(g_nlgr    , option)

  do ilgr = 1, g_nlgr
    call BroadcastInt(g_g(ilgr)%nx   , option)
    call BroadcastInt(g_g(ilgr)%ny   , option)
    call BroadcastInt(g_g(ilgr)%nz   , option)
    call BroadcastInt(g_g(ilgr)%nxy  , option)
    call BroadcastInt(g_g(ilgr)%nxyz , option)
    zarr(1) = g_g(ilgr)%ref_name
    call BroadcastCharN(zarr, option)
    g_g(ilgr)%ref_name = zarr(1)
  enddo

  call SetProbDims(g_g(g_ifld)%nx,g_g(g_ifld)%ny,g_g(g_ifld)%nz)

end subroutine DistributeCounts

! *************************************************************************** !

function DistVec3(va,vb)
  !
  ! Distance between two 3-vector points
  !
  ! Author: Dave Ponting
  ! Date: 05/05/20

  implicit none

  PetscReal :: DistVec3
  PetscReal, intent(in) :: va(g_ndir),vb(g_ndir)
  PetscReal :: dx, dy, dz

  dx = va(g_xdir) - vb(g_xdir)
  dy = va(g_ydir) - vb(g_ydir)
  dz = va(g_zdir) - vb(g_zdir)

  DistVec3 = sqrt(dx*dx + dy*dy + dz*dz)

end function DistVec3

! *************************************************************************** !

! *************************************************************************** !

subroutine DoOperate(ztarg, ixl, ixu, iyl, iyu, izl, izu, &
                      rarg1, qarg1, rarg2, qarg2, &
                      z_op_type, z_op_arg, qerr, ilgr,option)

  !
  ! Do an OPERATE operation
  !
  ! Author: Dave Ponting
  ! Date: 09/21/21
  !

  use String_module ,only : StringCompareIgnoreCase

  implicit none

  character(len = *), intent(in) :: ztarg
  PetscInt , intent(in) :: ixl, ixu, iyl, iyu, izl, izu
  PetscReal, intent(in) :: rarg1, rarg2
  PetscBool, intent(in) :: qarg1, qarg2
  character(len = *), intent(in) :: z_op_type, z_op_arg
  PetscBool, intent(out) :: qerr
  PetscInt , intent(in)  :: ilgr
  type(option_type) :: option

  PetscInt , pointer :: targi(:), operi(:)
  PetscReal, pointer :: targr(:), operr(:)

  PetscBool :: isintt, isrealt, isinta, isreala, &
               iseql, is_mult_only, is_depth, is_permt, is_perma ,is_dist, is_tran
  PetscBool :: ismulta, ismultiply, isminlim, ismaxlim , &
               iscopy , islog10   , isloge  , isabs    , &
               isinv  , ismultx   , isaddx  , ismultp  , &
               ispoly , isslog    , arg_not_perm, arg1_required, arg2_required
  PetscInt  :: ix, iy, iz, ig
  PetscReal :: vt, va, vr, mdtom2, m2tomd
  character(len = MAXSTRINGLENGTH) :: zmess
  character(len = 8) :: zv

!--Initialise------------------------------------------------------------------

  qerr = PETSC_FALSE

  isintt   = PETSC_FALSE
  isrealt  = PETSC_FALSE

  isinta   = PETSC_FALSE
  isreala  = PETSC_FALSE

  iseql        = PETSC_FALSE
  is_mult_only = PETSC_FALSE
  is_depth     = PETSC_FALSE
  is_permt     = PETSC_FALSE
  is_perma     = PETSC_FALSE
  is_dist      = PETSC_FALSE
  is_tran      = PETSC_FALSE

  mdtom2 = GetMDtoM2Conv()
  m2tomd = 1.0D0/mdtom2

!--Check that these are valid arrays for OPERATE-------------------------------

  call CheckValidOperateArray(ztarg   ,qerr)
  if(g_error_flag>0) return

  call CheckValidOperateArray(z_op_arg,qerr)
  if(g_error_flag>0) return

!--Get the target array--------------------------------------------------------

  call GetGridArrayPointer(ztarg, isintt, isrealt, targi, targr, qerr, iseql, &
                            is_mult_only, is_depth, is_permt, is_dist, is_tran, ilgr, option)
  if(g_error_flag>0) then
    qerr = PETSC_TRUE
    return
  endif

!--Get the argument array------------------------------------------------------

  call GetGridArrayPointer(z_op_arg, isinta, isreala, operi, operr, qerr, iseql, &
                            is_mult_only, is_depth, is_perma, is_dist, is_tran, ilgr, option)
  if(g_error_flag>0) then
    qerr = PETSC_TRUE
    return
  endif

!--Sort out the required operation---------------------------------------------

  ismulta    = StringCompareIgnoreCase(z_op_type, 'multa')
  iscopy     = StringCompareIgnoreCase(z_op_type, 'copy')
  ismultiply = StringCompareIgnoreCase(z_op_type, 'multiply')
  isminlim   = StringCompareIgnoreCase(z_op_type, 'minlim')
  ismaxlim   = StringCompareIgnoreCase(z_op_type, 'maxlim')
  islog10    = StringCompareIgnoreCase(z_op_type, 'log10')
  isloge     = StringCompareIgnoreCase(z_op_type, 'loge')
  isabs      = StringCompareIgnoreCase(z_op_type, 'abs')
  isinv      = StringCompareIgnoreCase(z_op_type, 'inv')
  ismultx    = StringCompareIgnoreCase(z_op_type, 'multx')
  isaddx     = StringCompareIgnoreCase(z_op_type, 'addx')
  ismultp    = StringCompareIgnoreCase(z_op_type, 'multp')
  ispoly     = StringCompareIgnoreCase(z_op_type, 'poly')
  isslog     = StringCompareIgnoreCase(z_op_type, 'slog')

!--Check for cases in which argument cannot be dimensionful-------------------

  arg_not_perm = isabs .or. islog10 .or. isloge .or. &
                 isinv .or. ismultp .or. ispoly .or. isslog
  if(arg_not_perm .and. is_perma) then
    zmess = 'OPERATE arg. '      //trim(z_op_arg) // &
            ' not possible with '//trim(z_op_arg) // &
            ' as permeability has dimension'
    call SetError(zmess, qerr)
  endif

!--Check for cases in which arguments expected---------------------------------

  arg1_required = ismulta .or. ismultx .or. isaddx .or. ismultp .or. ispoly .or. isslog
  if(arg1_required .and. (.not. qarg1)) then
    zmess = 'OPERATE arg. '// trim(z_op_arg) // ' requires first scalar value set'
    call SetError(zmess, qerr)
  endif

  arg2_required = ismulta .or. ismultp .or. ispoly .or. isslog
  if(arg2_required .and. (.not. qarg2)) then
    zmess = 'OPERATE arg. '// trim(z_op_arg) // ' requires second scalar value set'
    call SetError(zmess, qerr)
  endif

!--Carry out operation if recognised------------------------------------------

  if(ismulta  .or. iscopy  .or. ismultiply .or. isminlim .or. &
      ismaxlim .or. islog10 .or. isloge     .or. isabs    .or. &
      isinv    .or. ismultx .or. isaddx     .or. ismultp) then

!--Loop over the specified box of cells----------------------------------------

    do iz = izl, izu
      do iy = iyl, iyu
        do ix = ixl, ixu

          ig = GetNaturalIndex(ix, iy, iz, ilgr)

          ! Get current targ value
          ! For perms, convert to md

          if(isrealt) then
            if(is_permt) then
              vt = m2tomd*targr(ig)
            else
              vt =        targr(ig)
            endif
          else
            vt = real(targi(ig))
          endif

          ! Get current arg value
          ! For perms, convert to md

          if(isreala) then
            if(is_perma) then
              va = m2tomd*operr(ig)
            else
              va =        operr(ig)
            endif
          else
            va = real(operi(ig))
          endif

          ! Default is result to target

          vr = vt

          ! Perform operation if recognised

          if(ismulta) vr = rarg1*va + rarg2
          if(ismultiply) vr = vt*va
          if(isminlim) vr = min(rarg1,va)
          if(ismaxlim) vr = max(rarg1,va)
          if(iscopy) vr =       va
          if(isabs) vr = abs  (va)
          if(islog10) vr = log10(va)
          if(isloge) vr = log  (va)
          if(isinv) then
            if(abs(va)>0.0D0) vr = 1.0/va
          endif
          if(ismultx) vr = rarg1*va
          if(isaddx) vr = va + rarg1
          if(ismultp) vr =      rarg1*va**rarg2
          if(ispoly) vr = vr + rarg1*va**rarg2
          if(isslog) vr = 10.0**(rarg1+va*rarg2)

        ! Store result back in target array
        ! For perms, convert back to m2

          if(isrealt) then
            if(vr<0.0) then
              call PrintTidy8(vr,zv)
              zmess = 'OPERATE yields neg. value ' // zv
              call SetError(zmess, qerr)
              return
            endif
            if(is_permt) then
              targr(ig) = mdtom2*vr
            else
              targr(ig) =        vr
            endif
          else
            targi(ig) = nint(vr)
          endif

        enddo
      enddo
    enddo
  else
    zmess = 'OPERATE action ' // trim(z_op_type) // ' not recognised'
    call SetError(zmess, qerr)
  endif

end subroutine DoOperate

! *************************************************************************** !

subroutine ExtractCellDimensionsAndLocationsFromCPG(ilgr)
  !
  ! Extract cell dimensions and locations from COORD/ZCORN data
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscInt,intent(in) :: ilgr

  PetscInt  :: ix, iy, iz, ig
  PetscReal :: vdx, vdy, vdz, vpx, vpy, vpz

  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3)

  ! For each cell, find base offsets in the
  ! (2.nx).(2.ny).(2.nz) corner z-val array

  do ix = 1, g_g(ilgr)%nx
    do iy = 1, g_g(ilgr)%ny
      do iz = 1, g_g(ilgr)%nz

        ! Get the corners

        call GetCorners(ix, iy, iz, x000, x100, x010, x110, &
                         x001, x101, x011, x111, &
                         g_g(ilgr)%coord, g_g(ilgr)%zcorn, g_g(ilgr)%nx, g_g(ilgr)%ny)

        ! Given the cell corner locations,
        ! find the cell dimensions and locations

        call GetHexDims(x000, x100, x010, x110, &
                         x001, x101, x011, x111, &
                         vdx, vdy, vdz, vpx, vpy, vpz)

        ! Find cell index and store locations

        ig = GetNaturalIndex(ix, iy, iz, ilgr)

        g_g(ilgr)%dx  (ig) = vdx
        g_g(ilgr)%dy  (ig) = vdy
        g_g(ilgr)%dz  (ig) = vdz

        g_g(ilgr)%xloc(ig) = vpx
        g_g(ilgr)%yloc(ig) = vpy
        g_g(ilgr)%zloc(ig) = vpz

      enddo
    enddo
  enddo

end subroutine ExtractCellDimensionsAndLocationsFromCPG

! *************************************************************************** !

subroutine ExtractCPGFromCellDimensionsAndLocations(ilgr)
  !
  ! Extract corner point locations from cell size data
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscInt, intent(in) :: ilgr

  PetscReal :: xsum, xl, xu, ysum, yl, yu, dx, dy, dz, xc, yc, zc, xn, yn, zn
  PetscInt  :: ix, iy, iz, ig, ixnb, iynb, iznb, iox, ioy, ioz, jn

  ! Do the coords, looping over columns (use dx,dy in layer 1)

  iz = 1

  xsum = 0.0
  do ix = 1, g_g(ilgr)%nx

    iy = 1 ! (Use first line of y-values in finding x-spacing)
    ig = GetNaturalIndex(ix, iy, iz, ilgr)
    xl = xsum
    xu = xsum+g_g(ilgr)%dx(ig)
    xsum = xu

    ysum = 0.0

    do iy = 1, g_g(ilgr)%ny

      ig = GetNaturalIndex(ix, iy, iz, ilgr)
      yl = ysum
      yu = ysum+g_g(ilgr)%dy(ig)
      ysum = yu

                                       call setCoordLine(ix  , iy  , xl, yl,ilgr)
      if (ix == g_g(ilgr)%nx)          call setCoordLine(ix+1, iy  , xu, yl,ilgr)
      if (iy == g_g(ilgr)%ny)          call setCoordLine(ix  , iy+1, xl, yu,ilgr)
      if (ix == g_g(ilgr)%nx .and. &
          iy == g_g(ilgr)%ny)          call setCoordLine(ix+1, iy+1, xu, yu,ilgr)

    enddo
  enddo

  ! Now the zcorns

  do ix = 1, g_g(ilgr)%nx
    do iy = 1, g_g(ilgr)%ny
      do iz = 1, g_g(ilgr)%nz

        ig = GetNaturalIndex(ix, iy, iz, ilgr)

        xc = g_g(ilgr)%xloc(ig)
        yc = g_g(ilgr)%yloc(ig)
        zc = g_g(ilgr)%zloc(ig)

        dx = g_g(ilgr)%dx  (ig)
        dy = g_g(ilgr)%dy  (ig)
        dz = g_g(ilgr)%dz  (ig)

        ixnb = 2*(ix-1)
        iynb = 2*(iy-1)
        iznb = 2*(iz-1)

        do iox = 0, 1
          xn = xc + 0.5*(2*iox-1)*dx
          do ioy = 0, 1
            yn = yc + 0.5*(2*ioy-1)*dy
            do ioz = 0, 1
              zn = zc + 0.5*(2*ioz-1)*dz
              jn = 4*g_g(ilgr)%nxy*(iznb+ioz)+2*g_g(ilgr)%nx*(iynb+ioy)+ixnb+iox+1
              g_g(ilgr)%zcorn(jn) = zn
            enddo
          enddo
        enddo

      enddo
    enddo
  enddo

end subroutine ExtractCPGFromCellDimensionsAndLocations

! *************************************************************************** !

subroutine FillXYPositionsForLayer(iz, ilgr)
  !
  ! When building a grid from cell size data (dx,dy,dz) set up an iz-layer
  ! by filling in the x and y locations
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(in) :: iz, ilgr

  PetscInt :: ix, iy, ig
  PetscReal :: sum, dx, dy

  ! Fill in x-locations for each y-row

  do iy = 1, g_g(ilgr)%ny
    sum = 0.0
    do ix = 1, g_g(ilgr)%nx
      ig  = GetNaturalIndex(ix, iy, iz, ilgr)
      dx  = g_g(ilgr)%dx(ig)
      g_g(ilgr)%xloc(ig) = sum + 0.5*dx
      sum                = sum +     dx
    enddo
  enddo

  ! Fill in y-locations for each x-row

  do ix = 1, g_g(ilgr)%nx
    sum = 0.0
    do iy = 1, g_g(ilgr)%ny
      ig  = GetNaturalIndex(ix, iy, iz, ilgr)
      dy  = g_g(ilgr)%dy(ig)
      g_g(ilgr)%yloc(ig) = sum + 0.5*dy
      sum                = sum +     dy
    enddo
  enddo

end subroutine FillXYPositionsForLayer

! *************************************************************************** !

subroutine FillZPositionsForColumn(ix, iy, ilgr)
  !
  ! When building a grid from cell-size data (dx,dy,dz)
  ! set up an (ix,iy)-column by filling in the z locations
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(in) :: ix, iy, ilgr

  PetscInt :: iz, izabove, izbelow, ig, igabove, igbelow
  PetscReal :: dz

  ! First pass for tops (starting at top, going down)

  do izabove = g_g(ilgr)%nz, 2, -1
     izbelow = izabove-1

     igabove = GetNaturalIndex(ix, iy, izabove, ilgr)
     igbelow = GetNaturalIndex(ix, iy, izbelow, ilgr)

     if (Uninitialized(g_g(ilgr)%tops(igbelow))) then
       g_g(ilgr)%tops(igbelow) = g_g(ilgr)%tops(igabove) - g_g(ilgr)%dz(igabove)
     endif

   enddo

  ! Second pass for z-locations

   do iz = 1, g_g(ilgr)%nz
     ig = GetNaturalIndex(ix, iy, iz, ilgr)
     dz = g_g(ilgr)%dz(ig)
     g_g(ilgr)%zloc(ig) = g_g(ilgr)%tops(ig) - 0.5*dz
   enddo

end subroutine FillZPositionsForColumn

! *************************************************************************** !

subroutine GenerateGridConnections(bvg,option)
  !
  ! Given the cell locations, generate the connections
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscReal, allocatable, intent(in) :: bvg(:)
  type(option_type) :: option

  PetscInt  :: ix, iy, iz, jx, jy, jz, ig, igt, ia, icbase, nxyzt, ilgr
  PetscReal :: x, y, z, dx, dy, dz, mxi, myi, mzi, ntgi

  ! Set up grid connections

  print *, 'Calculating connections'

  nxyzt = 0
  do ilgr = 1,g_nlgr

  ! Allocate and set the temporary store for generated transmissilities

    allocate(g_ttx(g_g(ilgr)%nxyz))
    allocate(g_tty(g_g(ilgr)%nxyz))
    allocate(g_ttz(g_g(ilgr)%nxyz))

    g_ttx = 0.d0
    g_tty = 0.d0
    g_ttz = 0.d0

    g_tt_allocated = PETSC_TRUE

    do ix = 1, g_g(ilgr)%nx
      do iy = 1, g_g(ilgr)%ny
        do iz = 1, g_g(ilgr)%nz

          ig = GetNaturalIndex(ix, iy, iz, ilgr)
          igt = ig + nxyzt
          ia = g_g(ilgr)%gtoa(ig)
          if (ia>0) then

            dx  = g_g(ilgr)%dx(ig)
            dy  = g_g(ilgr)%dy(ig)
            dz  = g_g(ilgr)%dz(ig)

            mxi = g_g(ilgr)%mx (ig)
            myi = g_g(ilgr)%my (ig)
            mzi = g_g(ilgr)%mzn(ig)

            ntgi = g_g(ilgr)%ntg(ig)

            g_bvol(ia) = bvg(igt)

            g_x(ia) = g_g(ilgr)%xloc(ig)
            g_y(ia) = g_g(ilgr)%yloc(ig)
            g_z(ia) = g_g(ilgr)%zloc(ig)

            x = g_g(ilgr)%xloc(ig)
            y = g_g(ilgr)%yloc(ig)
            z = g_g(ilgr)%zloc(ig)

            ! Base count for this cell (held as g_nc will be incremented)

            icbase = g_nc

            jx = ix
            jy = iy
            jz = iz

            ! x-direction

            if (ix < g_g(ilgr)%nx) then
              jx = ix + 1
              jy = iy
              jz = iz
              ! Find tranx
              if (g_isnewtran) then
                call GetArea1(ia, ig, ix, iy, iz, ilgr, &
                                       jx, jy, jz, ilgr, g_xdir, &
                                       g_ipon_u, mxi, myi, mzi, ntgi,option)
              else
                call GetArea0(ia, ig, ix, iy, iz, g_xdir, &
                               dx, dy, dz, mxi, myi, mzi, ntgi, ilgr,option)
              endif
            endif

            ! y-direction

            if (iy< g_g(ilgr)%ny) then
              jx = ix
              jy = iy + 1
              jz = iz
              ! Find trany
              if (g_isnewtran) then
                call GetArea1(ia, ig, ix, iy, iz, ilgr, &
                                       jx, jy, jz, ilgr, g_ydir, &
                                       g_ipon_u, mxi, myi, mzi, ntgi,option)
              else
                call GetArea0(ia, ig, ix, iy, iz, g_ydir, &
                               dx, dy, dz, mxi, myi, mzi, ntgi, ilgr,option)
              endif
            endif

          ! z-direction

            if (iz < g_g(ilgr)%nz) then
              jx = ix
              jy = iy
              jz = iz + 1
              ! Find tranz
              if (g_isnewtran) then
                call GetArea1(ia, ig, ix, iy, iz, ilgr, &
                                       jx, jy, jz, ilgr, g_zdir, &
                                       g_ipon_u, mxi, myi, mzi, ntgi,option)
              else
                call GetArea0(ia, ig, ix, iy, iz, g_zdir, &
                               dx, dy, dz, mxi, myi, mzi, ntgi, ilgr,option)
              endif
            endif

          endif

        enddo ! enddo iz
      enddo ! enddo iy
    enddo ! enddo ix

  ! Overwrite any unset trans with the temporary generated ones

    do ig = 1, g_g(ilgr)%nxyz
      if (g_g(ilgr)%tx(ig)<0.0) g_g(ilgr)%tx(ig)=g_ttx(ig)
      if (g_g(ilgr)%ty(ig)<0.0) g_g(ilgr)%ty(ig)=g_tty(ig)
      if (g_g(ilgr)%tz(ig)<0.0) g_g(ilgr)%tz(ig)=g_ttz(ig)
    enddo

  ! Deallocate the temporary transmissibilities

    deallocate(g_ttx)
    deallocate(g_tty)
    deallocate(g_ttz)

    g_tt_allocated = PETSC_FALSE

    nxyzt = nxyzt + g_g(ilgr)%nx*g_g(ilgr)%ny*g_g(ilgr)%nz

  enddo ! endo g_ilgr

end subroutine GenerateGridConnections

! *************************************************************************** !

subroutine GenerateLGRConnections(option)
  !
  ! Generate connections between local grids
  !
  ! Author: Dave Ponting
  ! Date: 08/13/21

  implicit none

  type(option_type) :: option

  PetscInt  :: idir, ipon, ilgr
  PetscBool :: exists
  PetscInt  :: igxl, igxu ,igyl, igyu ,igzl, igzu, nrx, nry, nrz
  PetscInt  :: ibxl, ibxu ,ibyl, ibyu, ibzl, ibzu
  PetscInt  :: igx, igy, igz

  ! Loop over the refinements (not the global grid)

  print *, 'Calculating local grid connections'

  do ilgr = 2,g_nlgr

    allocate(g_ttx(g_g(ilgr)%nxyz))
    allocate(g_tty(g_g(ilgr)%nxyz))
    allocate(g_ttz(g_g(ilgr)%nxyz))

    g_ttx = 0.d0
    g_tty = 0.d0
    g_ttz = 0.d0

    g_tt_allocated = PETSC_TRUE

    ! Get limits of global grid cells in this LGR

    igxl = g_g(ilgr)%ixl
    igxu = g_g(ilgr)%ixu

    igyl = g_g(ilgr)%iyl
    igyu = g_g(ilgr)%iyu

    igzl = g_g(ilgr)%izl
    igzu = g_g(ilgr)%izu

    nrx  = g_g(ilgr)%nx
    nry  = g_g(ilgr)%ny
    nrz  = g_g(ilgr)%nz

!--Loop over faces of the refinement-------------------------------------------

    do idir=1,g_ndir
      do ipon = 1,2

        ! Check if this box can generate connections (may be at edge)

        exists = PETSC_TRUE

        ! Box of global cells to be considered

        ibxl = igxl
        ibxu = igxu
        ibyl = igyl
        ibyu = igyu
        ibzl = igzl
        ibzu = igzu

        if(idir == g_xdir) then
          if(ipon == g_ipon_l) then
            if (igxl == 1) exists = PETSC_FALSE
            ibxl = igxl
            ibxu = ibxl
          endif
          if(ipon == g_ipon_u) then
            if (igxu == g_g(g_ifld)%nx) exists = PETSC_FALSE
            ibxu = igxu
            ibxl = ibxu
          endif
        endif

        if(idir == g_ydir) then
          if(ipon == g_ipon_l) then
            if (igyl == 1) exists = PETSC_FALSE
            ibyl = igyl
            ibyu = ibyl
          endif
          if(ipon == g_ipon_u) then
            if (igyu == g_g(g_ifld)%ny) exists = PETSC_FALSE
            ibyu = igyu
            ibyl = ibyu
          endif
        endif

        if(idir == g_zdir) then
          if(ipon == g_ipon_l) then
            if (igzl == 1) exists = PETSC_FALSE
            ibzl = igzl
            ibzu = ibzl
          endif
          if(ipon == g_ipon_u) then
            if (igzu == g_g(g_ifld)%nz) exists = PETSC_FALSE
            ibzu = igzu
            ibzl = ibzu
          endif
        endif

!--If it exists, loop over the global cells in the refinement edge-------------

        if(exists) then

          do igx=ibxl, ibxu
            do igy=ibyl, ibyu
              do igz=ibzl, ibzu
                call GenerateLGRConnections1(idir , ipon , ilgr , &
                                             igx  , igy  , igz  , &
                                             igxl , igyl , igzl , option)
              enddo
            enddo
          enddo

        endif

      enddo
    enddo

    deallocate(g_ttx)
    deallocate(g_tty)
    deallocate(g_ttz)

    g_tt_allocated = PETSC_FALSE

  enddo

end subroutine GenerateLGRConnections

! *************************************************************************** !

subroutine GenerateLGRConnections1(idir, ipon  , ilgr, &
                                   igx  , igy  , igz , &
                                   igxl , igyl , igzl, option)
  !
  ! Given an edge cell in the current refinement (g_ilgr), find connections to
  ! global grid cells in the direction idir and the orientation ipon = g_ipon_l/u
  !
  ! Author: Dave Ponting
  ! Date: 08/13/21

  implicit none

  PetscInt,intent(in) :: idir, ipon  , ilgr , &
                         igx  , igy  , igz  , &
                         igxl , igyl , igzl
  type(option_type) :: option

  PetscInt :: irx , iry , irz , irg
  PetscInt :: irxl, iryl, irzl
  PetscInt :: irxu, iryu, irzu
  PetscInt :: nrx , nry , nrz , nrxy, ia
  PetscInt :: jgx , jgy , jgz
  PetscInt :: ibpgx, nlpgx
  PetscInt :: ibpgy, nlpgy
  PetscInt :: ibpgz, nlpgz

  PetscReal :: mxi , myi , mzi, ntgi

!--Data for this refinement-----------------------------------------------------

  nrx = g_g(ilgr)%nx
  nry = g_g(ilgr)%ny
  nrz = g_g(ilgr)%nz

  nrxy = nrx*nry

!--Get refined limits for this cell--------------------------------------------

  call GetLRForGC(igx,igxl,g_g(ilgr)%ibpgx,g_g(ilgr)%nlpgx,ibpgx,nlpgx)
  call GetLRForGC(igy,igyl,g_g(ilgr)%ibpgy,g_g(ilgr)%nlpgy,ibpgy,nlpgy)
  call GetLRForGC(igz,igzl,g_g(ilgr)%ibpgz,g_g(ilgr)%nlpgz,ibpgz,nlpgz)

  irxl = ibpgx+1
  irxu = ibpgx+nlpgx

  iryl = ibpgy+1
  iryu = ibpgy+nlpgy

  irzl = ibpgz+1
  irzu = ibpgz+nlpgz

!--Set neighbour in global grid--------------------------------------------------

   jgx = igx
   jgy = igy
   jgz = igz

!--Set to cells in the required face of the refinement---------------------------

  if((idir == g_xdir) .and. (ipon==g_ipon_l)) then
    irxu = 1            ! Move upper to the lower to get lower only
    jgx  = igx - 1
  endif

  if((idir == g_xdir) .and. (ipon==g_ipon_u)) then
    irxl = g_g(ilgr)%nx ! Move lower to the upper to get upper only
    jgx  = igx + 1
  endif

  if((idir == g_ydir) .and. (ipon==g_ipon_l)) then
    iryu = 1            ! Move upper to the lower to get lower only
    jgy  = igy - 1
  endif

  if((idir == g_ydir) .and. (ipon==g_ipon_u)) then
    iryl = g_g(ilgr)%ny ! Move lower to the upper to get upper only
    jgy  = igy + 1
  endif

  if((idir == g_zdir) .and. (ipon==g_ipon_l)) then
    irzu = 1            ! Move upper to the lower to get lower only
    jgz  = igz - 1
  endif

  if((idir == g_zdir) .and. (ipon==g_ipon_u)) then
    irzl = g_g(ilgr)%nz ! Move lower to the upper to get upper only
    jgz  = igz + 1
  endif

!--------------------------------------------------------------------------------------
!  Loop over the face of refined cells
!  In each case look for the connection from (irx,iry,irz) to (jgx,jgy,jgz)
!  and any other connections in the column (jgx,jgy,*)
!--------------------------------------------------------------------------------------

  do irx=irxl, irxu
    do iry=iryl, iryu
      do irz=irzl, irzu

        irg = nrxy*(irz-1)+nrx*(iry-1)+irx
        ia  = g_g(ilgr)%gtoa(irg)
        if (ia>0) then

          mxi  = g_g(ilgr)%mx (irg)
          myi  = g_g(ilgr)%my (irg)
          mzi  = g_g(ilgr)%mzn(irg)
          ntgi = g_g(ilgr)%ntg(irg)

          call GetArea1(ia, irg, irx, iry, irz, ilgr, jgx, jgy, jgz, g_ifld, idir, ipon, mxi, myi, mzi, ntgi, option)

        endif
      enddo
    enddo
  enddo

end subroutine GenerateLGRConnections1

! *************************************************************************** !

subroutine GetArea0(ia, ig, ix, iy, iz, idir, &
                    dxi, dyi, dzi, mxi, myi, mzi, ntgi, ilgr, option)
  !
  ! Do a block-based area calculation between cell (ix,iy,iz) and its
  ! positive neighbour in the idir-direction (g_xdir,..,g_z_dir)
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt , intent(in) :: ia, ig, ix, iy, iz, idir, ilgr
  PetscReal, intent(in) :: dxi, dyi, dzi, mxi, myi, mzi, ntgi
  type(option_type)     :: option

  PetscInt  :: jx, jy, jz, jg, ja, jlgr
  PetscReal :: xf, yf, zf, a, ai, aj, mxj, myj, mzj, ntgj, dxj, dyj, dzj, wi, wj, dd, dh
  PetscReal :: dxij, dyij, dzij, dip, dipc, tran

  tran = 0.d0

  ! Find neighbour

  jx   = ix
  jy   = iy
  jz   = iz
  jlgr = ilgr

  if (idir == g_xdir) jx = jx + 1
  if (idir == g_ydir) jy = jy + 1
  if (idir == g_zdir) jz = jz + 1

  jg = GetNaturalIndex(jx, jy, jz, ilgr)

  mxj = g_g(ilgr)%mxn(jg)
  myj = g_g(ilgr)%myn(jg)
  mzj = g_g(ilgr)%mz (jg)

  ja = g_g(ilgr)%gtoa(jg)

  if (ja > 0) then

!   Get ntg for cell j

    ntgj = g_g(ilgr)%ntg(jg)

!   Get size values in each direction

    dxj  = g_g(ilgr)%dx(jg)
    dyj  = g_g(ilgr)%dy(jg)
    dzj  = g_g(ilgr)%dz(jg)

    dxij = dxi + dxj
    dyij = dyi + dyj
    dzij = dzi + dzj

    ! Weights (if dxi>>dxj, wi->1, wj->0)

    wi = 0.5
    wj = 0.5

    ! Depth difference for dip corrections

    dd   = abs(g_g(ilgr)%zloc(ig)-g_g(ilgr)%zloc(jg))
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
      a  = mxi* mxj * (ai*wj + aj*wi)*dipc
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
      a  = myi *myj * (ai*wj + aj*wi)*dipc
    endif

    if (idir == g_zdir) then
      if (dzij>0.0) then
        wi = dzi/dzij
        wj = dzj/dzij
      endif
      ai = dxi * dyi
      aj = dxj * dyj
      a  = mzi* mzj *  (ai*wj + aj*wi)
    endif

    ! Set up the interface along the line from i to j, weighted by cell size

    xf = wj*g_g(ilgr)%xloc(ig)+wi*g_g(ilgr)%xloc(jg)
    yf = wj*g_g(ilgr)%yloc(ig)+wi*g_g(ilgr)%yloc(jg)
    zf = wj*g_g(ilgr)%zloc(ig)+wi*g_g(ilgr)%zloc(jg)

    ! Add in the connection

    if (a>0.0) then
      call   AddConnection(ia,ja,idir,g_ipon_u,xf,yf,zf,a,ilgr,jlgr,tran,a)
    endif

  endif

end subroutine GetArea0

! *************************************************************************** !

subroutine GetArea1(ia, ig, ix, iy, iz, ilgr, jx, jy, jz, jlgr, idir, ipon, mxi, myi, mzi, ntgi, option)
  !
  ! Do a corner-point-based area calculation between cell (ix,iy,iz) and its
  ! positive neighbour in the id-direction (g_xdir,..,g_z_dir)
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt , intent(in) :: ia, ig, ix, iy, iz, ilgr, jx, jy, jz, jlgr, idir, ipon
  PetscReal, intent(in) :: mxi, myi, mzi, ntgi
  type(option_type)     :: option

  PetscInt  :: jg, ja, kz, kg, ka
  PetscReal :: a, xf, yf, zf, fm, fntg , af, pod, mxj, myj, mzj, tran
  PetscBool :: match, matchingNeighbourFound, found

  tran = 0.d0

  ! Initialise flags

  match                  = PETSC_FALSE
  matchingNeighbourFound = PETSC_FALSE
  found                  = PETSC_FALSE

  ! Look at the natural neighbour

  jg = GetNaturalIndex(jx, jy, jz, jlgr)

  mxj = g_g(jlgr)%mxn(jg)
  myj = g_g(jlgr)%myn(jg)
  mzj = g_g(jlgr)%mz (jg)

  !  Find average ntg and set ntg factor

  fntg = 1.0
  if ((idir == g_xdir) &
      .or. (idir == g_ydir)) then
    fntg  = 0.5*(ntgi+g_g(jlgr)%ntg(jg))
  endif

  !  Set mult factor for first cell

  fm  = 1.0

  if (idir == g_xdir) then
   fm = mxi*mxj
  endif

  if (idir == g_ydir) then
   fm = myi*myj
  endif

  if (idir == g_zdir) then
    fm = mzi*mzj
  endif

  ja = g_g(jlgr)%gtoa(jg)
  if (ja > 0) then
    call GetMif(ix, iy, iz, ilgr, jx, jy, jz, jlgr, idir, ipon, a, xf, yf, zf, match, found)
    if (match) matchingNeighbourFound = PETSC_TRUE
    if (found) then
      af = a*fm*fntg
      if (af>0.0) then
        call   AddConnection(ia,ja,idir,g_ipon_u,xf,yf,zf,af,ilgr,jlgr,tran,a)
      endif
    endif
  endif

  ! Not a perfect match - look for other connections

  if (.not.matchingNeighbourFound) then

    ! In x- and y-dir, look for displacement faults to the jx,jy column

    if (idir == g_xdir .or. idir == g_ydir) then
      do kz = 1, g_g(jlgr)%nz
        if (kz /= jz) then !  Have looked already
          kg = GetNaturalIndex(jx, jy, kz, jlgr)

  !       Find ntg factor

          fntg  = 0.5*(ntgi+g_g(jlgr)%ntg(kg))

          ka = g_g(jlgr)%gtoa(kg)
          if (ka > 0) then
            call GetMif(ix, iy, iz, ilgr, jx, jy, kz, jlgr, idir, ipon, &
                        a, xf, yf, zf, match, found)
            if (found) then
              g_nf = g_nf + 1
              af = a*fm*fntg
              if (af > 0.0) then
                call   AddConnection(ia,ka,idir,g_ipon_u,xf,yf,zf,af,ilgr,jlgr,tran,a)
              endif
            endif
          endif
        endif
      enddo
    endif

    !  In z-dir, look for pinch-outs from iz to cells above jz in ix,iy column

    if (idir == g_zdir .and. iz<=(g_g(ilgr)%nz-2) .and. (ilgr==jlgr)) then

      ! Search up column from cell above z-neighbour

      column:do kz = iz+2, g_g(ilgr)%nz

        ! Look for active cell

        kg = GetNaturalIndex(ix, iy, kz, jlgr)
        ka = g_g(jlgr)%gtoa(kg)
        if (ka > 0) then

          ! Lower i-cell in Pft is mzi~g_mzn (backward/upwards   in Ecl terms)
          ! Upper k-cell in Pft is mzk~g_mz  (forwards/downwards in Ecl terms)

          fm = mzi*g_g(ilgr)%mz(kg)

          ! No ntg in vertical

          fntg = 1.0

          ! Find pinch-out distance

          call GetPod(ix, iy, iz, kz, jlgr, pod)
          if (pod<g_pinch(1)) then
            call GetMif(ix, iy, iz, ilgr, jx, jy, kz, jlgr, idir, ipon, &
                        a, xf, yf, zf, match, found)
            if (found) then
              af = a*fm*fntg
              if (af > 0.0) then
                 g_npo = g_npo + 1
                call   AddConnection(ia,ka,idir,g_ipon_u,xf,yf,zf,af,ilgr,jlgr,tran,a)
              endif
            endif
          endif

          !  Active cell found, finish search

          exit column
        endif
      enddo column
    endif

  endif

end subroutine GetArea1

! *************************************************************************** !

subroutine GetCellCoordinates(ig, ix, iy, iz, ilgr)
  !
  ! For a given natural cell index, extract the ix, iy, iz location
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(in) :: ig, ilgr
  PetscInt, intent(out) :: ix, iy, iz

  PetscInt :: ir

  iz = (ig-1)/g_g(ilgr)%nxy+1
  ir = ig-(iz-1)*g_g(ilgr)%nxy
  iy = (ir-1)/g_g(ilgr)%nx+1
  ix = ir-(iy-1)*g_g(ilgr)%nx

end subroutine GetCellCoordinates

! *************************************************************************** !

subroutine GetGridArrayPointer(za, isint, isreal, ai, ar, qerr, iseql, &
                               is_mult_only, is_depth, is_perm, &
                               is_dist, is_tran, ilgr, option)
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

  PetscBool, intent(in) :: iseql

  PetscBool, intent(out) :: is_mult_only
  PetscBool, intent(out) :: is_depth
  PetscBool, intent(out) :: is_perm
  PetscBool, intent(out) :: is_dist
  PetscBool, intent(out) :: is_tran

  PetscInt , intent(in)  :: ilgr
  type(option_type) :: option

  PetscBool, intent(out) :: qerr

  character(len=3) :: z3
  character(len=MAXWORDLENGTH) :: zwork

  zwork=adjustl(za)
  z3  = zwork(1:3)

  isint  = PETSC_FALSE
  isreal = PETSC_FALSE

  qerr = PETSC_FALSE

  is_mult_only = PETSC_FALSE
  is_depth     = PETSC_FALSE
  is_dist      = PETSC_FALSE
  is_tran      = PETSC_FALSE
  is_perm      = PETSC_FALSE

  nullify(ai, ar)

! Cell dimensions

  if (StringCompareIgnoreCase(za, 'dx')) then
    ar => g_g(ilgr)%dx
    is_dist = PETSC_TRUE
  endif
  if (StringCompareIgnoreCase(za, 'dy')) then
    ar => g_g(ilgr)%dy
    is_dist = PETSC_TRUE
  endif
  if (StringCompareIgnoreCase(za, 'dz')) then
    ar => g_g(ilgr)%dz
    is_dist = PETSC_TRUE
  endif

!  Permeabilities

  if (StringCompareIgnoreCase(za, 'permx')) then
    ar => g_g(ilgr)%kx
    is_perm = PETSC_TRUE
  endif

  if (StringCompareIgnoreCase(za, 'permy')) then
    ar => g_g(ilgr)%ky
    is_perm = PETSC_TRUE
  endif

  if (StringCompareIgnoreCase(za, 'permz')) then
    ar => g_g(ilgr)%kz
    is_perm = PETSC_TRUE
  endif

! Multipliers

  if (StringCompareIgnoreCase(za, 'multv')) ar => g_g(ilgr)%mv

  if (StringCompareIgnoreCase(za, 'multx')) ar => g_g(ilgr)%mx
  if (StringCompareIgnoreCase(za, 'multy')) ar => g_g(ilgr)%my
  if (StringCompareIgnoreCase(za, 'multz')) ar => g_g(ilgr)%mz

  if (StringCompareIgnoreCase(za, 'multx-')) ar => g_g(ilgr)%mxn
  if (StringCompareIgnoreCase(za, 'multy-')) ar => g_g(ilgr)%myn
  if (StringCompareIgnoreCase(za, 'multz-')) ar => g_g(ilgr)%mzn

! Porosity

  if (StringCompareIgnoreCase(za, 'poro')) ar => g_g(ilgr)%poro

! Tops (note this is a depth, so mark as such)

  if (StringCompareIgnoreCase(za, 'tops')) then
    ar => g_g(ilgr)%tops
    is_depth = PETSC_TRUE
    is_dist  = PETSC_TRUE
  endif

! Set up the type

  if (associated(ai)) then
    isint  = PETSC_TRUE
  else if (associated(ar)) then
    isreal = PETSC_TRUE
  else
    call SetError('Unable to find array ' // trim(za), qerr)
  endif

end subroutine GetGridArrayPointer

! *************************************************************************** !

function GetInOrAdjBox(ix,iy,ize,ixl,ixu,iyl,iyu,izle,izue,zilgr,zjlgr)
  !
  ! Check if cell (ix,iy,iz) is in or adjoins box (ixl..jzu)
  !
  ! Author: Dave Ponting
  ! Date: 09/21/21
  !

  implicit none

  PetscBool :: GetInOrAdjBox
  PetscInt, intent(in) :: ix,iy,ize,ixl,ixu,iyl,iyu,izle,izue
  character(len=*), intent(in) :: zilgr,zjlgr
  PetscBool :: qerr, in_x_range, in_y_range, in_z_range, &
                     ad_x_range, ad_y_range, ad_z_range, &
                     ad_x_error, ad_y_error, ad_z_error, over_error

!--Initialise------------------------------------------------------------------

  qerr = PETSC_FALSE

!--Set flags in each interval-------------------------------------------------

  ! In each range

  in_x_range = (ix >=ixl) .and. (ix <=ixu)
  in_y_range = (iy >=iyl) .and. (iy <=iyu)
  in_z_range = (ize>=izle) .and. (ize<=izue)

  ! Adjacent to each range

  ad_x_range = (ix ==ixl -1) .or. (ix ==ixu +1)
  ad_y_range = (iy ==iyl -1) .or. (iy ==iyu +1)
  ad_z_range = (ize==izle-1) .or. (ize==izue+1)

  ! Flags for overlap and adjacency conditions

  over_error = in_x_range .and. in_y_range .and. in_z_range

  ad_x_error = ad_x_range .and. in_y_range .and. in_z_range
  ad_y_error = in_x_range .and. ad_y_range .and. in_z_range
  ad_z_error = in_x_range .and. in_y_range .and. ad_z_range

  ! In all three ranges - cell is in LGR

  if(over_error) then
    call SetError('LGR overlap: '//trim(zilgr)//' and '//trim(zjlgr),qerr)
  endif

  ! Check if adjacent in x

  if((.not.qerr) .and. ad_x_error) then
    call SetError('LGRs adjoin in x: '//trim(zilgr)//' and '//trim(zjlgr),qerr)
  endif

  ! Check if adjacent in y

  if((.not.qerr) .and. ad_y_error) then
    call SetError('LGRs adjoin in y: '//trim(zilgr)//' and '//trim(zjlgr),qerr)
  endif

  ! Check if adjacent in z

  if((.not.qerr) .and. ad_z_error) then
    call SetError('LGRs adjoin in z: '//trim(zilgr)//' and '//trim(zjlgr),qerr)
  endif

  GetInOrAdjBox = qerr

end function GetInOrAdjBox

! *************************************************************************** !

subroutine GetLocalCount(ng, nl, nls, rem, option)
  !
  ! Find the local count size when splitting an array over all the ranks
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt, intent(in) :: ng
  PetscInt, intent(out) :: nl
  PetscInt, intent(out) :: nls
  PetscInt, intent(out) :: rem
  type(option_type)     :: option

  nl  = ng/option%comm%size
  nls = nl

  rem = ng - nl*option%comm%size
  if (option%myrank < rem) nl = nl + 1

end subroutine GetLocalCount

! *************************************************************************** !

subroutine GetLRForGC(ig, igl, ibpgv, nlpgv, ibpgc, nlpgc)

  ! For a global cell (igx,igy,igz) in refinement based at (igxl,igyl,igzl),
  ! find the 0-base pointer and the number of local cells
  !
  ! Author: Dave Ponting
  ! Date: 10/12/21
  !

  implicit none

  PetscInt,intent(in) :: ig
  PetscInt,intent(in) :: igl
  PetscInt,intent(in)  :: ibpgv(:)
  PetscInt,intent(in)  :: nlpgv(:)
  PetscInt,intent(out) :: ibpgc
  PetscInt,intent(out) :: nlpgc

  PetscInt :: iginlgr

  ! Set global grid indices within the refinement

  iginlgr = ig-igl+1

  ! Get base pointers for cell

  ibpgc = ibpgv(iginlgr)
  nlpgc = nlpgv(iginlgr)

end subroutine GetLRForGC

! *************************************************************************** !

function GetMultregtMultiplier(ia,ja,ilgr,jlgr)

  !
  ! Apply a MULTREGT muliplier between two cells
  !
  ! Author: Dave Ponting
  ! Date: 11/18/22
  !

  implicit none

  PetscReal :: GetMultregtMultiplier
  PetscInt, intent(in) :: ia,ja,ilgr,jlgr

  PetscInt  :: ig,jg,ira,jra,ir1,ir2, i_multregt
  PetscReal :: mult
  PetscBool :: matchi1, matchi2, &
               matchj1, matchj2, sameregion

  ! Initialise and set up cell natural grid indices

  mult = 1.0

  ig  = g_atog(ia)
  jg  = g_atog(ja)

  ! Cell region values

  ira = g_g(ilgr)%mltn(ig)
  jra = g_g(jlgr)%mltn(jg)

  ! Loop over MULTREGT records

  do i_multregt = 1,n_multregt

    ir1 = g_multregt_ir1(i_multregt)
    ir2 = g_multregt_ir2(i_multregt)

    ! If both regions are positive and the same, special same-region treatment

    sameregion = PETSC_FALSE
    if((ir1==ir2) .and. (ir1>0)) sameregion=PETSC_TRUE

    if(sameregion) then

      ! If same-region, set the multiplier if either cell is in this region
      ! So treats connection in this region, and from this region to others

      if((ira == ir1) .or. (jra == ir1)) then
        mult = g_multregt_mult(i_multregt)
      endif

    else

      ! If not same-region, only continue for cells in different regions

      if(ira /= jra) then

        ! Set up flags for cells matching regions, as positive values or wildcards

        matchi1 = (ira==ir1) .or. (ir1<=0)
        matchi2 = (ira==ir2) .or. (ir2<=0)
        matchj1 = (jra==ir1) .or. (ir1<=0)
        matchj2 = (jra==ir2) .or. (ir2<=0)

        ! If cell i matches region 1 and cell j matches region 2
        ! or cell i matches region 2 and cell j matches region 1, apply muliplier

        if ((matchi1 .and. matchj2) .or. &
             (matchi2 .and. matchj1)) then
          mult = g_multregt_mult(i_multregt)
        endif

      endif
    endif
  enddo

  GetMultRegtMultiplier = mult

end function GetMultregtMultiplier

! *************************************************************************** !

subroutine FindYLimitsAtThisX(fl, fu, idir, x, yll, yuu, dy, qoverlap)

  !
  ! Given two quadrilaterals fl and fu, find width of overlap region at x
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

!  Initialise dy to zero

  dy = 0.0

!  Set up the test line

  alx = x
  aux = x
  aly = yll - 1.0
  auy = yuu + 1.0

  !  For each quad, find the min and max intersections

  qfound = PETSC_FALSE
  ylumin = 0.0D0
  ylumax = 0.0D0

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
  !  or dy remains at zero (although not normally used if no overlap)

  qoverlap = qfound(1) .and. qfound(2)

  if(qoverlap) then
    ymin = max(ylumin(1), ylumin(2))
    ymax = min(ylumax(1), ylumax(2))

    if (ymax>ymin) dy = ymax-ymin
  endif

end subroutine FindYLimitsAtThisX

! *************************************************************************** !

function FindVolume(x000, x100, x010, x110, &
                     x001, x101, x011, x111)

  !
  ! Obtain the volume of a hexahederal grid cell
  ! with corners x000,..,x111
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscReal FindVolume
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

  vp0 = FindVolume1(1, 2, 3, c)
  vp1 = FindVolume1(1, 3, 2, c)

  vp2 = FindVolume1(2, 1, 3, c)
  vp3 = FindVolume1(2, 3, 1, c)

  vp4 = FindVolume1(3, 1, 2, c)
  vp5 = FindVolume1(3, 2, 1, c)

  !  Sum with appropriate signs into the final pore volume

  vp = vp0 - vp1 - (vp2 - vp3) + vp4 - vp5

  ! Take absolute value to get pore volume

  FindVolume = abs(vp)

end function FindVolume

! *************************************************************************** !

function FindVolume1(id0, id1, id2, c)

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

  PetscReal :: FindVolume1
  PetscInt , intent(in) :: id0, id1, id2
  PetscReal, intent(in) :: c(0:1, 0:1, 0:1, 3)

  PetscInt :: ja, jb, jc, jd, je, jf
  PetscReal :: den

  FindVolume1 = 0.0

  do ja = 0, 1
    do jb = 0, 1
      do jc = 0, 1
        do jd = 0, 1
          do je = 0, 1
            do jf = 0, 1

              den = real((jc + je + 1) * (ja + jf + 1) * (jb + jd + 1), &
                          kind(g_petsc_real))
              FindVolume1 = FindVolume1             &
                           + c(1 , ja, jb, id0)     &
                           * c(jc, 1 , jd, id1)     &
                           * c(je, jf, 1 , id2) / den

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end function  FindVolume1

! *************************************************************************** !

subroutine GetHexDims(x000, x100, x010, x110, &
                       x001, x101, x011, x111, &
                       vdx, vdy, vdz, vpx, vpy, vpz)
  !
  ! Given cell corner locations, extract cell dimensions and locations
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscReal, intent(in) :: x000(3), x100(3), x010(3), x110(3), &
                            x001(3), x101(3), x011(3), x111(3)
  PetscReal, intent(out) :: vdx, vdy, vdz, vpx, vpy, vpz

  PetscInt  :: idir
  PetscReal :: dx(3), dy(3), dz(3), prx, pry

  ! Cell centre is just the average

  vpx = 0.125*(x000(g_xdir)+x100(g_xdir)+x010(g_xdir)+x110(g_xdir) &
               +x001(g_xdir)+x101(g_xdir)+x011(g_xdir)+x111(g_xdir))

  vpy = 0.125*(x000(g_ydir)+x100(g_ydir)+x010(g_ydir)+x110(g_ydir) &
               +x001(g_ydir)+x101(g_ydir)+x011(g_ydir)+x111(g_ydir))

  vpz = 0.125*(x000(g_zdir)+x100(g_zdir)+x010(g_zdir)+x110(g_zdir) &
               +x001(g_zdir)+x101(g_zdir)+x011(g_zdir)+x111(g_zdir))

  ! Cell dimensions are distances between faces: form x,y and z intervals

  do idir = 1, g_ndir

                     ! 1**        0**
    dx(idir) = 0.25*(x100(idir)-x000(idir) &
                     +x110(idir)-x010(idir) &
                     +x101(idir)-x001(idir) &
                     +x111(idir)-x011(idir))

                     ! *1*        *0*
    dy(idir) = 0.25*(x010(idir)-x000(idir) &
                     +x110(idir)-x100(idir) &
                     +x011(idir)-x001(idir) &
                     +x111(idir)-x101(idir))

                     ! **1        **0
    dz(idir) = 0.25*(x001(idir)-x000(idir) &
                     +x101(idir)-x100(idir) &
                     +x011(idir)-x010(idir) &
                     +x111(idir)-x110(idir))

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

subroutine GetMif(ix, iy, iz, ilgr, jx, jy, jz, jlgr, idir, ipon, a, xf, yf, zf, match, found)

  !
  ! Obtain the interface area between cell (ix,iy,iz) and cell (jx,jy,jz)
  ! in the direction idir
  ! Return the face centre and the scalar interface area
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt, intent(in) :: ix, iy, iz, ilgr, jx, jy, jz, jlgr, idir, ipon
  PetscReal, intent(out) :: a, xf, yf, zf
  PetscBool, intent(out) :: match, found

  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3), &
               fl(0:1, 0:1, 3), fu(0:1, 0:1, 3), &
               fcl(3), fcu(3), cl(3), cu(3), dlu(3), &
               d, ax, ay, az, srayl, srayu, srayl2, srayu2, &
               sdlu, sdlu2, elu(3), sdlui, f, sraylu, auns, &
               dl,du,topl,botl,topu,botu, &
               signx, signy, signz, au

  PetscReal, parameter :: missdistance = 0.01 ! 1cm
  PetscBool, parameter :: positive_yes = PETSC_TRUE
  PetscBool, parameter :: positive_no  = PETSC_FALSE
  PetscReal, parameter :: micro        = 1.0D-6

  PetscBool :: missed, positivei, positivej

  PetscInt  :: i, j, k, nrxi, nryi, nrzi, irgi, nrxj, nryj, nrzj, irgj

  PetscReal :: diff

  if(ipon ==  g_ipon_u) then
    positivei = positive_yes
    positivej = positive_no
  else
    positivei = positive_no
    positivej = positive_yes
  endif

  ! Default return values

     a  = 0.0
    xf  = 0.0
    yf  = 0.0
    zf  = 0.0
  match = PETSC_TRUE
  found = PETSC_FALSE

  ! Get the lower and upper cell centres and faces

  call GetCorners(ix, iy, iz, &
                   x000, x100, x010, x110, &
                   x001, x101, x011, x111, &
                   g_g(ilgr)%coord, g_g(ilgr)%zcorn, g_g(ilgr)%nx, g_g(ilgr)%ny)

  call GetFace(idir, positivei, &
                x000, x100, x010, x110, &
                x001, x101, x011, x111, fl)

  call GetCorners(jx, jy, jz, &
                   x000, x100, x010, x110, &
                   x001, x101, x011, x111, &
                   g_g(jlgr)%coord, g_g(jlgr)%zcorn, g_g(jlgr)%nx, g_g(jlgr)%ny)

  call GetFace(idir, positivej, &
                x000, x100, x010, x110, &
                x001, x101, x011, x111, fu)

  nrxi = g_g(ilgr)%nx
  nryi = g_g(ilgr)%ny
  nrzi = g_g(ilgr)%nz

  irgi = nrxi*nryi*(iz-1)+nrxi*(iy-1)+ix

  nrxj = g_g(jlgr)%nx
  nryj = g_g(jlgr)%ny
  nrzj = g_g(jlgr)%nz

  irgj = nrxj*nryj*(jz-1)+nrxj*(jy-1)+jx

  cl(g_xdir)=g_g(ilgr)%xloc(irgi)
  cl(g_ydir)=g_g(ilgr)%yloc(irgi)
  cl(g_zdir)=g_g(ilgr)%zloc(irgi)

  cu(g_xdir)=g_g(jlgr)%xloc(irgj)
  cu(g_ydir)=g_g(jlgr)%yloc(irgj)
  cu(g_zdir)=g_g(jlgr)%zloc(irgj)

  ! Check for depth miss if x or y direction face

  missed = PETSC_FALSE

  if (idir /= g_zdir) then

    ! Get first value from elevation of the (0,0) quad element

    dl=fl(0, 0, g_zdir)
    topl = dl
    botl = dl

    du=fu(0, 0, g_zdir)
    topu = du
    botu = du

   ! Loop over the face vertices to find the top and bottom vertex depths

    do i = 0, 1
      do j = 0, 1

        dl=fl(i, j, g_zdir)
        topl = max(topl,dl)
        botl = min(botl,dl)

        du=fu(i, j, g_zdir)
        topu = max(topu,du)
        botu = min(botu,du)

      enddo
    enddo

    ! If the bottom of fl is above the top    of fu, missed
    ! If the top    of fl is below top bottom of fu, missed

    if ((botl >= (topu-micro)) &
        .or. (topl <= (botu+micro))) missed=PETSC_TRUE

  endif

  ! If not a depth miss, calculate area

  if (missed) then
    match=PETSC_FALSE
  else

    do i = 1, 3
      dlu(i) = cu(i)-cl(i)
      fcl(i) = 0.25D0*(fl(0, 0, i)+fl(0, 1, i)+fl(1, 0, i)+fl(1, 1, i))
      fcu(i) = 0.25D0*(fu(0, 0, i)+fu(0, 1, i)+fu(1, 0, i)+fu(1, 1, i))
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

    !  Do simple quad area calculation of fl

    ax = qarea(fl, g_xdir)
    ay = qarea(fl, g_ydir)
    az = qarea(fl, g_zdir)

    !  If no match, use mutual interfce area calculation

    if (.not.match) then

      !  As the mutual interface is a sub-area of each fl projection,
      !  it has the same sign, so use the signs of that projection

      signx = 1.0
      signy = 1.0
      signz = 1.0

      if (ax<0.0) signx = -1.0
      if (ay<0.0) signy = -1.0
      if (az<0.0) signz = -1.0

      !  Only do mutual calc if fl area non-zero (mutual cannot be larger)

      if (abs(ax)>0.0) then
        au = MArea(fl, fu, g_xdir)
        ax = signx*abs(au)
      endif

      if (abs(ay)>0.0) then
        au = MArea(fl, fu, g_ydir)
        ay = signy*abs(au)
      endif

      if (abs(az)>0.0) then
        au = MArea(fl, fu, g_zdir)
        az = signz*abs(au)
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

function GetNaturalIndex(ix, iy, iz, ilgr)
  !
  ! Get the natural cell index (in ix/iy/iz order with ix varying fastest)
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt :: GetNaturalIndex
  PetscInt, intent(in) :: ix, iy, iz, ilgr

  GetNaturalIndex = g_g(ilgr)%nxy*(iz-1) + g_g(ilgr)%nx*(iy-1) + ix

end function GetNaturalIndex

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

subroutine GetPod(ix, iy, iz, jz, ilgr, pod)

  !
  ! Obtain gap between top of cell (ix,iy,iz) and bottom of cell (ix,iy,jz)
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt , intent(in)  :: ix, iy, iz, jz, ilgr
  PetscReal, intent(out) :: pod

  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3), d, &
               fl(0:1, 0:1, 3), fu(0:1, 0:1, 3)

  PetscBool, parameter :: positive_yes = PETSC_TRUE
  PetscBool, parameter :: positive_no  = PETSC_FALSE

  PetscInt  :: i, j

  ! Get the lower and upper cell centres and faces

  call GetCorners(ix, iy, iz, &
                   x000, x100, x010, x110, &
                   x001, x101, x011, x111, &
                   g_g(ilgr)%coord, g_g(ilgr)%zcorn, g_g(ilgr)%nx, g_g(ilgr)%ny)

  call GetFace(g_zdir, positive_yes, &
                x000, x100, x010, x110, &
                x001, x101, x011, x111, fl)

  call GetCorners(ix, iy, jz, &
                   x000, x100, x010, x110, &
                   x001, x101, x011, x111, &
                   g_g(ilgr)%coord, g_g(ilgr)%zcorn, g_g(ilgr)%nx, g_g(ilgr)%ny)

  call GetFace(g_zdir, positive_no, &
                x000, x100, x010, x110, &
                x001, x101, x011, x111, fu)

  pod = 0.0
  do i = 0, 1
    do j = 0, 1
      d = abs(fu(i, j, g_zdir)-fl(i, j, g_zdir))
      pod = max(d, pod)
    enddo
  enddo

end subroutine GetPod

! *************************************************************************** !

function GetPosNeighbour(ix, iy, iz, idir, jx, jy, jz, ilgr)
  !
  ! Find the neighbour of (ix,iy,iz) in the direction idir
  !
  ! Author: Dave Ponting
  ! Date: 05/20/20

  implicit none

  PetscBool :: GetPosNeighbour
  PetscInt, intent(in)  :: ix, iy, iz, idir, ilgr
  PetscInt, intent(out) :: jx, jy, jz
  PetscBool :: ingrid

  ingrid = PETSC_TRUE

  jx = ix
  jy = iy
  jz = iz

  if (idir == g_xdir) jx = ix + 1
  if (idir == g_ydir) jy = iy + 1
  if (idir == g_zdir) jz = iz + 1

  if ((jx>g_g(ilgr)%nx) .or. (jy>g_g(ilgr)%ny) .or. (jz>g_g(ilgr)%nz)) ingrid = PETSC_FALSE

  GetPosNeighbour = ingrid

end function GetPosNeighbour

! *************************************************************************** !

subroutine AddConnection(ia, ja, idir, ipon, xf, yf, zf, a, ilgr, jlgr, &
                         tran, tran_th)

  !
  !  Store connection between cells ia and ja, face location (xf,yf,zf), area a
  !
  ! Author: Dave Ponting
  ! Date: 02/03/18

  implicit none

  PetscInt , intent(in) :: ia, ja, idir, ipon, ilgr, jlgr
  PetscReal, intent(in) :: xf, yf, zf, a, tran, tran_th

  PetscReal, parameter :: rm6 = 1.0E-6
  PetscInt :: ig, jg, iloc
  PetscReal :: xa, ya, za, xb, yb, zb, f, xt, yt, zt, err, mult
  PetscReal :: distx, disty, distz
  PetscReal, parameter :: err_eps=1.0D-3

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

  g_tran   (g_nc) = tran
  g_tran_th(g_nc) = tran_th
  if(ipon == g_ipon_u) then
    g_idir(g_nc) =  idir
  else
    g_idir(g_nc) = -idir
  endif

  if (g_rmltn == 1) then
    mult=GetMultregtMultiplier(ia,ja,ilgr,jlgr)
    g_tran (g_nc) = mult*g_tran (g_nc)
    g_carea(g_nc) = mult*g_carea(g_nc)
  endif

  if (g_geometry_test) then

    ig = g_atog(ia)
    jg = g_atog(ja)

    xa = g_g(ilgr)%xloc(ig)
    ya = g_g(ilgr)%yloc(ig)
    za = g_g(ilgr)%zloc(ig)

    xb = g_g(jlgr)%xloc(jg)
    yb = g_g(jlgr)%yloc(jg)
    zb = g_g(jlgr)%zloc(jg)

    distx = abs(xb-xa)
    disty = abs(yb-ya)
    distz = abs(zb-za)

    iloc = g_xdir
    if(disty>distx) iloc=g_ydir
    if(distz>distx .and. distz>disty) iloc=g_zdir

    f= 0.5
    if (iloc==g_xdir) then
      f = (xf-xa)/(xb-xa)
    else if (iloc==g_ydir) then
      f = (yf-ya)/(yb-ya)
    else if (iloc==g_zdir) then
      f = (zf-za)/(zb-za)
    endif

    xt = xa+f*(xb-xa)
    yt = ya+f*(yb-ya)
    zt = za+f*(zb-za)

    err = sqrt((xt-xf)**2+(yt-yf)**2+(zt-zf)**2)

    if (err>err_eps) then
      print *,'x:xt,yt,zt ', xt, yt, zt
      print *,'x:xf,yf,zf ', xf, yf, zf
      print *,'x:xa,ya,za ', xa, ya, za
      print *,'x:xb,yb,zb ', xb, yb, zb
      print *,'f,iloc,err ', f,iloc,err
    endif

  endif

end subroutine AddConnection

! *************************************************************************** !

subroutine getFace(idir, positive, &
                    x000, x100, x010, x110, &
                    x001, x101, x011, x111, face)

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

subroutine GrdeclReader(input, option)
  !
  ! Reads an Eclgrid file on I/O proc
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18
  !
  use String_module
  use PFLOTRAN_Constants_module

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

  PetscInt  :: ierr, nread, nxsg, nysg, nzsg, gridfileoption, &
                ilgr, jlgr
  PetscBool :: qerr, gridunit_is_metric, gridunit_is_feet, ismetric, &
                newtran_force, oldtran_force
  PetscReal :: dimens(3),convdist, convrvol, convtran, convfull

  PetscReal, parameter :: convnull = 1.0

  character(len = MAXWORDLENGTH)   :: word
  character(len = MAXSTRINGLENGTH) :: zbuf(mzbuf)
  character(len = MAXSTRINGLENGTH) :: zmess, efile_name, wordlong

  ierr = INPUT_ERROR_NONE
  qerr = PETSC_FALSE
  zbuf = ' '
  dimens   = 1.0
  convdist = 1.0
  convrvol = 1.0
  convtran = 1.0 ! Conversion of transmissibilities from FIELD to METRIC (not to SI)
  convfull = 1.0 ! Full conversion: FIELD to METRIC if required, then METRIC to Pflotran

  newtran_force = PETSC_FALSE
  oldtran_force = PETSC_FALSE

  ismetric      = PETSC_TRUE

  ilgr = g_ifld

  ! Read through items in the grdecl file

  g_column = 1

  call InputPushBlock(input,option)

  do

    call InputReadPflotranStringNotComment(input, option)
    if (option%error_while_nonblocking) then
      ! I/O error detected so set qerr and exit
      call SetError('Unable to read from GRDECL file', qerr)
      exit
    endif
    if (InputError(input)) exit

  ! Keyword found

    call InputReadCard(input, option, word)
    call CheckError(input, 'ECLGRD subkeyword', qerr)
    if (qerr) exit
    call StringToUpper(word)

    if (word(1:1) /= '/') print *, 'Pflotran card:: ', word

    select case(trim(word))
      case('DIMENS')
        call ReadEvalues(dimens, 3, 'DIMENS', 'GRID', &
                          input, option, convnull, nn_yes, qerr)
        if (.not. qerr) call SetDimens(dimens)
      case('CARFIN')
        ilgr = HandleCarfin(word, input, option, qerr)
        g_gridfileoption = 2
      case('NXFIN')
        call HandleNdfin(word, ilgr, g_xdir, input, option, qerr)
      case('NYFIN')
        call HandleNdfin(word, ilgr, g_ydir, input, option, qerr)
      case('NZFIN')
        call HandleNdfin(word, ilgr, g_zdir, input, option, qerr)
      case('HXFIN')
        call HandleHdfin(word, ilgr, g_xdir, input, option, qerr)
      case('HYFIN')
        call HandleHdfin(word, ilgr, g_ydir, input, option, qerr)
      case('HZFIN')
        call HandleHdfin(word, ilgr, g_zdir, input, option, qerr)
      case('ENDFIN')
        ilgr = g_ifld
      case('COORD')
        call IsCPG(g_ifld)
        g_isnewtran = PETSC_TRUE
        call CheckDimensRead(qerr)
        if (.not.qerr) call ReadECoordArray(g_g(ilgr)%coord,input,option,convdist,qerr,g_ifld)
      case('ZCORN')
        call IsCPG(g_ifld)
        g_isnewtran = PETSC_TRUE
        call CheckDimensRead(qerr)
        if (.not.qerr) call ReadEZcornArray(g_g(ilgr)%zcorn,input,option,convdist,qerr,g_ifld)
      case('GDFILE')
        call IsCPG(g_ifld)
        g_isnewtran = PETSC_TRUE
        call CheckDimensRead(qerr)
        if(.not.qerr) then
          call ReadEstrings(word, zbuf, 1, input, option, qerr)
          if(.not.qerr) then
            efile_name=zbuf(1)
            call ReadEclipseEgridFile(efile_name,g_g(g_ifld)%coord,g_g(g_ifld)%zcorn,g_g(g_ifld)%actn, &
                                      g_g(g_ifld)%nx,g_g(g_ifld)%ny,g_g(g_ifld)%nz,qerr,zmess)
            if(qerr) then
              call SetError(zmess, qerr)
            endif
          endif
        endif
      case('DX')
        call ReadEGridArrR(g_g(ilgr)%dx, 'DX'   , &
                            input, option, convdist,dep_no,prm_no,nn_yes,qerr,ilgr)
      case('DY')
        call ReadEGridArrR(g_g(ilgr)%dy, 'DY'   , &
                            input, option, convdist,dep_no,prm_no, nn_yes,qerr,ilgr)
      case('DZ')
        call ReadEGridArrR(g_g(ilgr)%dz, 'DZ'   , &
                            input, option, convdist,dep_no,prm_no, nn_yes,qerr,ilgr)

      case('PERMX')
        call ReadEGridArrR(g_g(ilgr)%kx, 'PERMX', &
                            input, option, convnull,dep_no,prm_yes,nn_yes,qerr,ilgr)
      case('PERMY')
        call ReadEGridArrR(g_g(ilgr)%ky, 'PERMY', &
                            input, option, convnull,dep_no, prm_yes,nn_yes,qerr,ilgr)
      case('PERMZ')
        call ReadEGridArrR(g_g(ilgr)%kz, 'PERMZ', &
                            input, option, convnull,dep_no, prm_yes,nn_yes,qerr,ilgr)

      case('MULTV')
        call ReadEGridArrR(g_g(ilgr)%mv, 'MULTV', &
                            input, option, convnull,dep_no, prm_no ,nn_yes,qerr,ilgr)

      case('MULTX')
        call ReadEGridArrR(g_g(ilgr)%mx, 'MULTX', &
                            input, option, convnull,dep_no, prm_no ,nn_yes,qerr,ilgr)
      case('MULTY')
        call ReadEGridArrR(g_g(ilgr)%my, 'MULTY', &
                            input, option, convnull,dep_no, prm_no ,nn_yes,qerr,ilgr)
      case('MULTZ')
        call ReadEGridArrR(g_g(ilgr)%mz, 'MULTZ', &
                            input, option, convnull,dep_no, prm_no ,nn_yes,qerr,ilgr)

      case('MULTX-')
        call ReadEGridArrR(g_g(ilgr)%mxn, 'MULTX-', &
                            input, option, convnull,dep_no, prm_no ,nn_yes,qerr,ilgr)
      case('MULTY-')
        call ReadEGridArrR(g_g(ilgr)%myn, 'MULTY-', &
                            input, option, convnull,dep_no, prm_no ,nn_yes,qerr,ilgr)
      case('MULTZ-')
        call ReadEGridArrR(g_g(ilgr)%mzn, 'MULTZ-', &
                            input, option, convnull,dep_no, prm_no ,nn_yes,qerr,ilgr)

      case('PORO')
        call ReadEGridArrR(g_g(ilgr)%poro, 'PORO', &
                            input, option, convnull, dep_no, prm_no,nn_yes,qerr,ilgr)
      case('TOPS')
        call ReadEGridArrR(g_g(ilgr)%tops, 'TOPS', &
                            input, option, convdist, dep_yes,prm_no,nn_no ,qerr,ilgr)
      case('NTG','NTOG')
        call ReadEGridArrR(g_g(ilgr)%ntg, trim(word), &
                            input, option, convnull, dep_no ,prm_no,nn_yes,qerr,ilgr)

      case('TRANX')
        convfull = convtran*GetETtoTrConv()
        call ReadEGridArrR(g_g(ilgr)%tx, 'TRANX', &
                            input, option, convfull, dep_no ,prm_no,nn_yes,qerr,ilgr)
      case('TRANY')
        convfull = convtran*GetETtoTrConv()
        call ReadEGridArrR(g_g(ilgr)%ty, 'TRANY', &
                            input, option, convfull, dep_no ,prm_no,nn_yes,qerr,ilgr)
      case('TRANZ')
        convfull = convtran*GetETtoTrConv()
        call ReadEGridArrR(g_g(ilgr)%tz, 'TRANZ', &
                            input, option, convfull, dep_no ,prm_no,nn_yes,qerr,ilgr)

      case('ACTNUM')
        call ReadEGridArrI(g_g(ilgr)%actn, 'ACTNUM', input, option, nn_yes, qerr,ilgr)
      case('SATNUM')
        g_rsatn = 1
        call ReadEGridArrI(g_g(ilgr)%satn, 'SATNUM', input, option, nn_yes, qerr,ilgr)
      case('IMBNUM')
        g_rimbn = 1
        call ReadEGridArrI(g_g(ilgr)%imbn, 'IMBNUM', input, option, nn_yes, qerr, ilgr)
      case('EQLNUM')
        g_reqln = 1
        call ReadEGridArrI(g_g(ilgr)%eqln, 'EQLNUM', input, option, nn_yes, qerr, ilgr)
      case('MULTNUM')
        g_rmltn = 1
        call ReadEGridArrI(g_g(ilgr)%mltn, 'MULTNUM', input, option, nn_yes, qerr, ilgr)
      case('MINPV')
        call ReadEvalues(g_minpv, 1, 'MINPV', 'GRID' , &
                          input, option, convrvol, nn_yes, qerr)
      case('PVFLOOR')
        call ReadEvalues(g_pvfloor, 1, 'PVFLOOR', 'GRID' , &
                          input, option, convrvol, nn_yes, qerr)
      case('PINCH')
        call ReadEvalues(g_pinch, 1, 'PINCH', 'GRID' , &
                          input, option, convdist, nn_yes, qerr)
      case('DPCF')
        call ReadEvalues(g_dpcf, 2, 'DPCF', 'GRID' , &
                          input, option, convnull, nn_yes, qerr)
        g_isdpcf = PETSC_TRUE
      case('NEWTRAN')
        g_isnewtran   = PETSC_TRUE
        newtran_force = PETSC_TRUE
      case('OLDTRAN')
        g_isnewtran   = PETSC_FALSE
        oldtran_force = PETSC_TRUE
      case('CREATE_TRAN_CONNECTIONS')
        g_create_tc   = PETSC_TRUE
      case('ADD', 'COPY', 'EQUALS', 'MULTIPLY')
        call HandleACEMKeyword(word, input, option, convdist, convtran, qerr,ilgr)
      case('OPERATE')
        call HandleOPERATEKeyword(word, input, option, qerr,ilgr)
      ! case('FAULTS')
      !   call ReadFaults(word, input, option, qerr,g_ifld)
      ! case('MULTFLT')
      !   call ReadMultflt(word, input, option, qerr)
      case('ECHO', 'NOECHO', 'INIT', 'EDIT', 'REGIONS')  ! No need to action
      case('FILEUNIT')
        call ReadEstrings(word, zbuf, 1, input, option, qerr)
      case('MAPUNITS')
        call ReadEstrings(word, zbuf, 1, input, option, qerr)
        call SetMapAxesUnits(zbuf(1))
      case('GRIDUNIT')
        call ReadEstrings(word, zbuf, 2, input, option, qerr)
        gridunit_is_metric = StringCompareIgnoreCase(zbuf(1), 'metres')
        gridunit_is_feet   = StringCompareIgnoreCase(zbuf(1), 'feet')
        if (.not. (gridunit_is_metric .or. gridunit_is_feet)) then
          call SetError('GRIDUNIT units must be metres or feet', qerr)
        else if (gridunit_is_feet) then
          convdist = 0.3048D0                  ! Convert from ft to m
        else if (gridunit_is_metric) then
          convdist = 1.0
        endif
      case('GDORIENT')
        zbuf(1) = 'INC'
        zbuf(2) = 'INC'
        zbuf(3) = 'INC'
        zbuf(4) = 'DOWN'
        zbuf(5) = 'RIGHT'
        call ReadEstrings(word, zbuf, 5, input, option, qerr)
        call SetGDORIENT(zbuf,qerr,zmess)
        if(qerr) then
          call SetError(zmess,qerr)
        endif
      case('GRIDFILE')
        zbuf(1) = '0' ! Keyword default is zero
        gridfileoption = 0
        call ReadEstrings(word, zbuf, 2, input, option, qerr)
        wordlong = zbuf(1)
        read(wordlong,*, iostat=ierr) gridfileoption
        if ((gridfileoption < 0) .or. &
              (gridfileoption > 2) .or. InputError(ierr)) then
          call SetError('GRIDFILE arg 1 must be 0,1 or 2 ', qerr)
        else
          g_gridfileoption = gridfileoption
        endif
      case('SPECGRID')
        call ReadEstrings(word, zbuf, 5, input, option, qerr)
          wordlong = zbuf(1)
          read(wordlong, *, iostat=ierr) nxsg
          if (nxsg /= g_g(g_ifld)%nx .or. InputError(ierr)) then
            call SetError('SPECGRID NX ' &
                          // trim(wordlong) // ' does not match DIMENS NX', qerr)
          endif
          wordlong = zbuf(2)
          read(wordlong, *, iostat=ierr) nysg
          if (nysg /= g_g(g_ifld)%ny .or. InputError(ierr)) then
            call SetError('SPECGRID NY ' &
                          // trim(wordlong) // ' does not match DIMENS NY', qerr)
          endif
          wordlong = zbuf(3)
          read(wordlong, *, iostat=ierr) nzsg
          if (nzsg /= g_g(g_ifld)%nz .or. InputError(ierr)) then
            call SetError('SPECGRID NZ ' &
                          // trim(wordlong) // ' does not match DIMENS NZ', qerr)
          endif
      case('MAPAXES')
        nread = size(g_mapaxes)
        call ReadEvalues(g_mapaxes, nread, 'MAPAXES', 'GRID', &
                          input, option, convnull, nn_no, qerr)
        call SetMapAxes(g_mapaxes)
      ! case('NNC')
      !   call HandleNNC(word, input, option, convtran, qerr)
      ! case('MULTREGT')
      !   call HandleMULTREGT(word, input, option, qerr)
      case('FIELD')
        ismetric = PETSC_FALSE
        convdist = 0.3048D0                    ! Convert from ft to m
        convrvol = 0.158987294928D0            ! Convert from rb to rm3
        convtran = 14.503774D0/6.28981077043D0 ! Convert rb.cp/psi/d to rm3.cp/bar/d
      case('METRIC')
        ismetric = PETSC_TRUE
        convdist = 1.0
        convrvol = 1.0
        convtran = 1.0
      case('/')  ! Isolated un-used terminator on its own line is harmless
      case default
        zmess = 'GRDECL sub-keyword ' // trim(word) // ' not recognised'
        call SetError(zmess, qerr)
    end select

    !  Leave if error flag set

    if (qerr) exit
  enddo
  call InputPopBlock(input,option)

  ! Set up maximum satnum, imbnum and eqlnum values

  if (g_rsatn == 1) then
    g_maxsatn = 1
    do jlgr=1,g_nlgr
      g_maxsatn = max(g_maxsatn,maxval(g_g(jlgr)%satn))
    enddo
  endif

  if (g_rimbn == 1) then
    g_maximbn = 1
    do jlgr=1,g_nlgr
      g_maximbn = max(g_maximbn,maxval(g_g(jlgr)%imbn))
    enddo
  endif

  if (g_reqln == 1) then
    g_maxeqln = 1
    do jlgr=1,g_nlgr
      g_maxeqln = max(g_maxeqln,maxval(g_g(jlgr)%eqln))
    enddo
  endif

  if (g_rmltn == 1) then
    g_maxmltn = 1
    do jlgr=1,g_nlgr
      g_maxmltn = max(g_maxmltn,maxval(g_g(jlgr)%mltn))
    enddo
  endif

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

    call ProcessGridData(qerr,option)
    if (g_na == 0) then

!  Case of no active cells: error

      zmess = 'Problem has no active cells'
      call SetError(zmess, qerr)

    else

    ! Process the well data read

      ! call ProcessWellData(qerr)

    endif

  endif

end subroutine GrdeclReader

! *************************************************************************** !

subroutine GetNextWord(word, exitTime, input, option)
  !
  ! Get the next white-space-delimited item on the input file
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  character(len = *), intent(out) :: word
  PetscBool, intent(out) :: exitTime
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len = MAXSTRINGLENGTH) :: line
  character(len = MAXSTRINGLENGTH) :: test

  character c
  PetscInt :: i, l, j, startcol, ic
  PetscBool :: lastCharRead, somethingRead, iws

  PetscBool :: started,inquotes,quote,started_in_quotes

  started           = PETSC_FALSE
  somethingRead     = PETSC_FALSE
  started_in_quotes = PETSC_FALSE

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
      if (quote) inquotes = .not.inquotes
      ! If not started, any non-whitespace character causes start
      if ((.not.started) .and. (.not.iws)) then
        started = PETSC_TRUE
        if(inquotes) started_in_quotes=PETSC_TRUE
      endif
      ! If started, stop if any whitespace not in quotes or end of line
      if (started .and. ((iws.and. (.not.inquotes)) .or. (i == l))) exit
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
      word = trim   (test)
    endif

  ! Check for / terminator

    c = word(1:1)
    if (c == '/' .and. (.not.started_in_quotes)) then
      exitTime = PETSC_TRUE
    endif

  endif

end subroutine GetNextWord

! *************************************************************************** !

function GridEclipseGetNumArrSet(maxnum,id_num_arr)
  !
  ! Get flag indicating that num array (eg satnum, imbnum, eqlnum, multnum) has been set
  !
  ! Author: Dave Ponting
  ! Date: 02/14/18
  !

  implicit none

  PetscBool :: GridEclipseGetNumArrSet
  PetscInt, intent(in)  :: id_num_arr
  PetscInt, intent(out) :: maxnum

  GridEclipseGetNumArrSet = PETSC_FALSE

  if ((id_num_arr == g_id_satnum_arr) .and. &
       (g_rsatn    == 1)) then
        GridEclipseGetNumArrSet = PETSC_TRUE
    maxnum       = g_maxsatn
  endif

  if ((id_num_arr == g_id_imbnum_arr) .and. &
       (g_rimbn    == 1)) then
    GridEclipseGetNumArrSet = PETSC_TRUE
    maxnum       = g_maximbn
  endif

  if ((id_num_arr == g_id_eqlnum_arr) .and. &
       (g_reqln    == 1)) then
    GridEclipseGetNumArrSet = PETSC_TRUE
    maxnum       = g_maxeqln
  endif

  if ((id_num_arr == g_id_mltnum_arr) .and. &
       (g_rmltn    == 1)) then
    GridEclipseGetNumArrSet = PETSC_TRUE
    maxnum       = g_maxmltn
  endif

end function GridEclipseGetNumArrSet

! *************************************************************************** !

subroutine GridEclipseGetPorPerm(ia, poro, permx, permy, permz)
  !
  ! Get the perm and porosity values of a given active cell
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  implicit none

  PetscInt , intent(in)  :: ia
  PetscReal, intent(out) :: poro, permx, permy, permz

  PetscInt :: ig, ilgr

  ig   = g_atog(ia)
  ilgr = g_atol(ia)

  poro  = g_g(ilgr)%poro(ig)
  permx = g_g(ilgr)%kx  (ig)
  permy = g_g(ilgr)%ky  (ig)
  permz = g_g(ilgr)%kz  (ig)

end subroutine GridEclipseGetPorPerm

! *************************************************************************** !

function GridEcilpseGetProcnumFlag()
  !
  ! Get the procnum required flag
  !
  ! Author: Dave Ponting
  ! Date: 03/10/22

  implicit none

  PetscInt :: GridEcilpseGetProcnumFlag

  GridEcilpseGetProcnumFlag = g_rprcn

end function GridEcilpseGetProcnumFlag

! *************************************************************************** !

subroutine GetRecordDetails(iRec, n, mInRec, il, iu, nInRec)
  !
  ! Given the number of values to be written, get details of (iRec)th record
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscInt, intent(in) :: iRec, n, mInRec
  PetscInt, intent(out) :: il, iu, nInRec

  PetscInt :: iumax

  iumax = n
  il = mInRec*(iRec-1)+1
  iu = il+mInRec-1
  if (iu>iumax) iu = iumax
  nInRec = iu-il+1

end subroutine GetRecordDetails

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
  if ((c == ' ') .or. &
      (ic == g_ictab) .or. &
      (ic == g_iccr) .or. &
      (ic == g_iclf)) isWhiteSpace = PETSC_TRUE

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
      (ic == g_dquote)) isQuote = PETSC_TRUE

end function IsQuote

! *************************************************************************** !

function HandleCarfin(zkey, input, option, qerr)

  !
  ! Handle a CARFIN (local grid refinement) keyword
  !
  ! Author: Dave Ponting
  ! Date: 08/04/21
  !

  use String_module

  implicit none

  PetscInt :: HandleCarfin
  character(len=*), intent(in) :: zkey
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr

  PetscInt, parameter :: ma = 10

  character(len = MAXSTRINGLENGTH) :: za(ma)

  PetscInt :: ixl, ixu, iyl, iyu, izle, izue, izl, izu, nxref, nyref, nzref, &
              nxf, nyf, nzf, ilgr
  PetscInt, parameter :: ibig = 10000000

!  Initialise values

  qerr = PETSC_FALSE

!  Initialise values (note read in Eclipse iz convention)

  ixl   = 1
  ixu   = 1

  iyl   = 1
  iyu   = 1

  izle  = 1
  izue  = 1

  nxref = 1
  nyref = 1
  nzref = 1

  nxf = g_g(g_ifld)%nx
  nyf = g_g(g_ifld)%ny
  nzf = g_g(g_ifld)%nz

  za   = ' '

  ilgr = g_ifld

  call ReadEstrings(zkey, za, ma, input, option, qerr)

  if(.not.qerr) then

    ! Note check location in field grid - no nesting

                    call ProcessArgToInt(ixl , za(2), zkey, 1   , nxf, qerr)
    if(.not.qerr) call ProcessArgToInt(ixu , za(3), zkey, ixl , nxf, qerr)

    if(.not.qerr) call ProcessArgToInt(iyl , za(4), zkey, 1   , nyf, qerr)
    if(.not.qerr) call ProcessArgToInt(iyu , za(5), zkey, iyl , nyf, qerr)

    if(.not.qerr) call ProcessArgToInt(izle, za(6), zkey, 1   , nzf, qerr)
    if(.not.qerr) call ProcessArgToInt(izue, za(7), zkey, izle, nzf, qerr)

    if(.not.qerr) call ProcessArgToInt(nxref, za(8), zkey, 1   , ibig, qerr)
    if(.not.qerr) call ProcessArgToInt(nyref, za(9), zkey, 1   , ibig, qerr)
    if(.not.qerr) call ProcessArgToInt(nzref, za(10), zkey, 1   , ibig, qerr)

    if(.not.qerr) call CheckLGRLocation(ixl, ixu, iyl, iyu, izle, izue, &
                                          za(1), qerr)

    if(.not.qerr) then

!   Convert to Pflotran layer convention
!   Note lower E-conv. layer becomes the upper P-conv. layer and v.v.
!   Note these are in field

      izl = g_g(g_ifld)%nz-izue+1
      izu = g_g(g_ifld)%nz-izle+1

      g_nlgr = g_nlgr + 1
      ilgr   = g_nlgr

      g_g(ilgr)%ref_name = za(1)

      g_g(ilgr)%ixl = ixl
      g_g(ilgr)%ixu = ixu

      g_g(ilgr)%iyl = iyl
      g_g(ilgr)%iyu = iyu

      g_g(ilgr)%izl = izl
      g_g(ilgr)%izu = izu

      g_g(ilgr)%nxg = ixu - ixl + 1
      g_g(ilgr)%nyg = iyu - iyl + 1
      g_g(ilgr)%nzg = izu - izl + 1

      call AllocateGridArrays(nxref,nyref,nzref,ilgr)

    endif

  endif

  HandleCarfin = ilgr

end function HandleCarfin

! *************************************************************************** !

subroutine HandleHdfin(zkey, ilgr, idir, input, option, qerr)

  !
  ! Handle a HdFIN (HXFIN, HYFIN, HZFIN) local grid refinement keyword
  !
  ! Author: Dave Ponting
  ! Date: 12/13/21
  !

  use String_module

  implicit none

  character(len=*), intent(in) :: zkey
  PetscInt, intent(in) :: ilgr
  PetscInt, intent(in) :: idir
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr
  PetscInt :: idl, ndl

  PetscReal, allocatable :: ra(:)
  character(len = MAXSTRINGLENGTH), allocatable :: za(:)

!  Initialise values

  qerr = PETSC_FALSE

  ndl = 1
  if (idir == g_xdir) ndl = g_g(ilgr)%nx
  if (idir == g_ydir) ndl = g_g(ilgr)%ny
  if (idir == g_zdir) ndl = g_g(ilgr)%nz

  allocate(za(ndl))
  allocate(ra(ndl))
  ra = 1.0d0

  call ReadEstrings(zkey, za, ndl, input, option, qerr)

  if(.not.qerr) then

    do idl=1, ndl
      call ProcessArgToReal(ra(idl), za(idl), zkey, qerr)
      if(qerr) exit
    enddo

  endif

  if (idir == g_xdir) g_g(ilgr)%hrefx = ra
  if (idir == g_ydir) g_g(ilgr)%hrefy = ra
  if (idir == g_zdir) then
    do idl=1,ndl
      g_g(ilgr)%hrefz(idl) = ra(ndl-idl+1)
    enddo
  endif

  g_g(ilgr)%hdset(idir) = PETSC_TRUE

  deallocate(za)
  deallocate(ra)

end subroutine HandleHdfin

! *************************************************************************** !

subroutine HandleNdfin(zkey, ilgr, idir, input, option, qerr)

  !
  ! Handle a NdFIN (NXFIN, NYFIN, NZFIN) local grid refinement keyword
  !
  ! Author: Dave Ponting
  ! Date: 10/11/21
  !

  use String_module

  implicit none

  character(len=*), intent(in) :: zkey
  PetscInt, intent(in) :: ilgr
  PetscInt, intent(in) :: idir
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr
  PetscInt :: idg, ndl, ndg, ndlis

  PetscInt, allocatable :: ia(:)
  character(len = MAXSTRINGLENGTH), allocatable :: za(:)
  character(len = 8) :: zndlis,zndl

!  Initialise values

  qerr = PETSC_FALSE

  if (idir == g_xdir) then
    ndl = g_g(ilgr)%nx
    ndg = g_g(ilgr)%nxg
  endif

  if (idir == g_ydir) then
    ndl = g_g(ilgr)%ny
    ndg = g_g(ilgr)%nyg
  endif

  if (idir == g_zdir) then
    ndl = g_g(ilgr)%nz
    ndg = g_g(ilgr)%nzg
  endif

  allocate(za(ndg))
  allocate(ia(ndg))
  ia=0

  call ReadEstrings(zkey, za, ndg, input, option, qerr)

  if(.not.qerr) then

    do idg=1, ndg
      call ProcessArgToInt(ia(idg), za(idg), zkey, 1, ndl, qerr)
    enddo

  endif

  if (idir == g_xdir) then
    g_g(ilgr)%nlpgx = ia
    call MakeLocalBasePointers(g_g(ilgr)%nlpgx,g_g(ilgr)%ibpgx,ndl,ndg,qerr,ndlis)
  endif

  if (idir == g_ydir) then
    g_g(ilgr)%nlpgy = ia
    call MakeLocalBasePointers(g_g(ilgr)%nlpgy,g_g(ilgr)%ibpgy,ndl,ndg,qerr,ndlis)
  endif

  if (idir == g_zdir) then
    do idg=1,ndg
      g_g(ilgr)%nlpgz(idg) = ia(ndg-idg+1)
    enddo
    call MakeLocalBasePointers(g_g(ilgr)%nlpgz,g_g(ilgr)%ibpgz,ndl,ndg,qerr,ndlis)
  endif

  if(qerr) then
    write(zndlis,'(I8)') ndlis
    write(zndl  ,'(I8)') ndl
    call SetError('Local cell total in '     // trim(adjustl(zkey)) // &
                                  ' is '     // trim(adjustl(zndlis)) // &
                                  ', expect '// trim(adjustl(zndl)) , qerr)
   endif

  deallocate(za)
  deallocate(ia)

end subroutine HandleNdfin

! *************************************************************************** !

subroutine HandleACEMKeyword(zkey, input, option, convdist, convtran, qerr, ilgr)
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
  PetscReal, intent(in) :: convdist, convtran
  PetscBool, intent(inout) :: qerr
  PetscInt , intent(in) :: ilgr

  PetscInt, parameter :: ma = 8
  character(len = MAXSTRINGLENGTH) :: za(ma)

  PetscBool :: isadd, iscpy, iseql, ismlt
  PetscInt  :: ixl, ixu, iyl, iyu, izl, izu, izle, izue

!  Initialise values

  qerr = PETSC_FALSE

  isadd = StringCompareIgnoreCase(zkey, 'add')
  iscpy = StringCompareIgnoreCase(zkey, 'copy')
  iseql = StringCompareIgnoreCase(zkey, 'equals')
  ismlt = StringCompareIgnoreCase(zkey, 'multiply')

!  Initialise values (note read in Eclipse iz convention)

  ixl  = 1
  ixu  = g_g(ilgr)%nx

  iyl  = 1
  iyu  = g_g(ilgr)%ny

  izle = 1
  izue = g_g(ilgr)%nz

  do

    za    = ' '
    za(1) = '/'
    call ReadEstrings(zkey, za, ma, input, option, qerr)
    if (qerr) exit

    if (StringCompareIgnoreCase(za(1), '/')) exit

    call ProcessArgToInt(ixl , za(3), zkey, 1   , g_g(ilgr)%nx, qerr)
    if (qerr) exit
    call ProcessArgToInt(ixu , za(4), zkey, ixl , g_g(ilgr)%nx, qerr)
    if (qerr) exit

    call ProcessArgToInt(iyl , za(5), zkey, 1   , g_g(ilgr)%ny, qerr)
    if (qerr) exit
    call ProcessArgToInt(iyu , za(6), zkey, iyl , g_g(ilgr)%ny, qerr)
    if (qerr) exit

    call ProcessArgToInt(izle, za(7), zkey, 1   , g_g(ilgr)%nz, qerr)
    if (qerr) exit
    call ProcessArgToInt(izue, za(8), zkey, izle, g_g(ilgr)%nz, qerr)
    if (qerr) exit

!   Convert to Pflotran layer convention
!   Note lower E-conv. layer becomes the upper P-conv. layer and v.v.

    izl = g_g(ilgr)%nz-izue+1
    izu = g_g(ilgr)%nz-izle+1

    call ACEMOp(za(1), za(2), ixl, ixu, iyl, iyu, izl, izu, &
                 isadd, iscpy, iseql, ismlt, zkey, qerr, &
                 convdist, convtran ,ilgr,option)
    if (qerr) exit

  enddo

end subroutine HandleACEMKeyword

! *************************************************************************** !

subroutine HandleOperateKeyword(zkey, input, option, qerr, ilgr)
  !
  ! Handle OPERATE
  !
  ! Author: Dave Ponting
  ! Date: 09/20/21

  use String_module

  implicit none

  character(len=*), intent(in) :: zkey
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr
  PetscInt , intent(in) :: ilgr

  PetscInt, parameter :: ma = 11
  character(len = MAXSTRINGLENGTH) :: za(ma)

  PetscInt  :: ixl, ixu, iyl, iyu, izl, izu, izle, izue
  PetscReal :: rarg1, rarg2
  character(len = MAXSTRINGLENGTH) :: z_op_type, z_op_arg
  PetscBool :: qarg1, qarg2

!  Initialise values

  qerr = PETSC_FALSE

!  Initialise values (note read in Eclipse iz convention)

  ixl  = 1
  ixu  = g_g(ilgr)%nx

  iyl  = 1
  iyu  = g_g(ilgr)%ny

  izle = 1
  izue = g_g(ilgr)%nz

  do

    za    = ' '
    za(1) = '/'
    call ReadEstrings(zkey, za, ma, input, option, qerr)
    if (qerr) exit

    if (StringCompareIgnoreCase(za(1), '/')) exit

    call ProcessArgToInt(ixl , za(2), zkey, 1   , g_g(ilgr)%nx, qerr)
    if (qerr) exit
    call ProcessArgToInt(ixu , za(3), zkey, ixl , g_g(ilgr)%nx, qerr)
    if (qerr) exit

    call ProcessArgToInt(iyl , za(4), zkey, 1   , g_g(ilgr)%ny, qerr)
    if (qerr) exit
    call ProcessArgToInt(iyu , za(5), zkey, iyl , g_g(ilgr)%ny, qerr)
    if (qerr) exit

    call ProcessArgToInt(izle, za(6), zkey, 1   , g_g(ilgr)%nz, qerr)
    if (qerr) exit
    call ProcessArgToInt(izue, za(7), zkey, izle, g_g(ilgr)%nz, qerr)
    if (qerr) exit

    z_op_type = za(8)
    z_op_arg  = za(9)

    qarg1 = PETSC_FALSE
    rarg1 = 0.0
    if(za(10) /= ' ') then
      call ProcessArgToReal(rarg1, za(10), zkey, qerr)
      if (qerr) exit
      qarg1 = PETSC_TRUE
    endif

    rarg2 = 0.0
    qarg2 = PETSC_FALSE
    if(za(11) /= ' ') then
      call ProcessArgToReal(rarg2, za(11), zkey, qerr)
      if (qerr) exit
      qarg2 = PETSC_TRUE
    endif

!   Convert to Pflotran layer convention
!   Note lower E-conv. layer becomes the upper P-conv. layer and v.v.

    izl = g_g(ilgr)%nz-izue+1
    izu = g_g(ilgr)%nz-izle+1

    call DoOperate(za(1), ixl, ixu, iyl, iyu, izl, izu, &
                    rarg1, qarg1, rarg2, qarg2, &
                    z_op_type, z_op_arg, qerr, ilgr, option)
    if (qerr) exit

  enddo

end subroutine HandleOperateKeyword

subroutine InputReadPflotranStringNotComment(input, option)

  !
  ! Read a new input line ,throwing away Eclipse comment lines (start with --)
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
    if (InputError(input)) exit                   ! Detect eof and leave
    word = adjustl(input%buf)
    if (word(1:2) /= '--') exit
  enddo

  g_word=word

end subroutine InputReadPflotranStringNotComment

! *************************************************************************** !

! *************************************************************************** !

subroutine IsCPG(ilgr)
  !
  ! If a keyword encountered indicating corner point input, (ZCORN or COORD)
  ! set flags and allocate the ccord and zcorn arrays
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscInt, intent(in) :: ilgr

  g_iscpg = PETSC_TRUE

  call CheckAllocCPGA(ilgr)

end subroutine IsCPG

! *************************************************************************** !

function MArea(fl, fu, idir)

  !
  ! Routine to get component id of the vector area of the quadrilateral f
  !
  ! Author: Dave Ponting
  ! Date: 02/05/18

  implicit none

  PetscReal MArea
  PetscReal, intent(in) :: fl(0:1, 0:1, 3)
  PetscReal, intent(in) :: fu(0:1, 0:1, 3)
  PetscInt , intent(in) :: idir
  PetscReal :: xv(24), xvdd(24), x , y, diff, dy1, dy2

  PetscInt  :: jdir, kdir, nxv, nxvd, i, j, ierr
  PetscReal :: alx, aly, aux, auy, &
               blx, bly, bux, buy, px, py, yll, yuu, x1, x2
  PetscBool :: qinter, qoverlap1, qoverlap2
  PetscCount :: sortcount

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

  diff  = 0.01
  MArea = 0.0

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

  !  Look for intersections

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

  !  Sort x-values into order

  ierr = 0
  sortcount = nxv
  call PetscSortReal(sortcount, xv, ierr)

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

    MArea = 0.0

    do i = 1, nxvd-1

      x1 = xvdd(i)
      x2 = xvdd(i+1)

      call FindYLimitsAtThisX(fl, fu, idir, x1, yll, yuu, dy1, qoverlap1)
      call FindYLimitsAtThisX(fl, fu, idir, x2, yll, yuu, dy2, qoverlap2)

      if (qoverlap1 .and. qoverlap2) then
        MArea = MArea + 0.5*(x2-x1)*(dy1+dy2)
      endif

    enddo

  endif

end function MArea

! *************************************************************************** !

subroutine OpenEgridFileToRead(root, ios)
  !
  ! Open an Egrid file to read
  ! This routine is called by only the I/O rank
  !
  ! Author: Dave Ponting
  ! Date: 11/08/21

  use String_module,only : StringToUpper

  implicit none

  character(len = *), intent(in) :: root
  PetscInt, intent(out) :: ios

  PetscInt  :: isuf, isuf_f, isuf_u
  PetscBool :: nosuffix
  character(len = MAXSTRINGLENGTH) :: filename,filenameu

  ! Initialise

  ios           = 0
  e_formatted_r = PETSC_FALSE

  ! Start by assuming root is the full filename

  filename  = root

  ! Check for suffix in case-independent manner

  filenameu = root
  call StringToUpper(filenameu)
  isuf_f = index(filenameu,'.FEGRID')
  isuf_u = index(filenameu,'.EGRID')
  isuf  = isuf_f + isuf_u
  nosuffix = (isuf==0)

  ! Attempt to open specified files

  if (isuf_f>0) then
    ! As .FGRID suffix found, open formatted with full name
    open(unit =  UNIT_EGRID_READ, status = 'old', form = 'formatted', &
         file = filename, iostat = ios)
    e_formatted_r = PETSC_TRUE
  else
    ! Is unformatted - add suffix if not found
    if (nosuffix) filename = trim(root) // '.EGRID'
    open(unit =  UNIT_EGRID_READ, status = 'old', form = 'unformatted', &
         file = filename, convert = 'big_endian', iostat = ios)
    ! Might be suffix case issue, try lower case
    if (nosuffix .and. (ios /= 0)) then
      filename = trim(root) // '.egrid'
      open(unit =  UNIT_EGRID_READ, status = 'old', form = 'unformatted', &
           file = filename, convert = 'big_endian', iostat = ios)
    endif
  endif

end subroutine OpenEgridFileToRead

! *************************************************************************** !

subroutine GridEclipsePorPermExchangeAndSet(poro_p, permx_p, permy_p, permz_p,&
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
  PetscReal :: poro, permx, permy, permz

  ! Work arrays

  PetscInt , allocatable :: winat(:)
  PetscReal, allocatable :: wporo(:)
  PetscReal, allocatable :: wperm(:)

  ! Set scalars

  ierr   = 0
  iorank = option%comm%io_rank
  t_rank = option%myrank

  ! Loop over exchange operations between IO ranks and non-IO ranks (irank)

  do irank = 0, option%comm%size-1

    if (irank /= iorank) then

      ! Consider exchange between irank and iorank

      if (t_rank  == irank) then

        ! This is irank: send nlmax value and then irequest array to ioproc

        temp_int_array(1) = nlmax
        call MPI_Send(temp_int_array, ONE_INTEGER_MPI, MPI_INTEGER, &
                      iorank, tag_mpi_0, option%mycomm, ierr)
        int_mpi = nlmax
        call MPI_Send(inatsend      , int_mpi        , MPI_INTEGER, &
                      iorank, tag_mpi_0, option%mycomm, ierr)

        ! Allocate buffers to hold response poro and perm values from ioproc

        allocate(wporo(nlmax))
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
          poro_p (ilt) = wporo(ilt)
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

        allocate(winat(nlo))
        allocate(wporo(nlo))
        allocate(wperm(3*nlo))

        ! Receive the natural address array from irank

        int_mpi = nlo
        call MPI_Recv(winat, int_mpi, MPI_INTEGER, irank, MPI_ANY_TAG, &
                      option%mycomm, status_mpi, ierr)

        ! On this (io) proc, set up the perm and poro values to be returned

        do ilo = 1, nlo
          ino = winat(ilo)
          call GridEclipseGetPorPerm(ino, poro, permx, permy, permz)
          ibp = 3*(ilo-1)
          wporo(ilo) = poro
          wperm(ibp+1) = permx
          wperm(ibp+2) = permy
          wperm(ibp+3) = permz
        enddo

        ! Send back the poro and perm values to irank

        int_mpi = nlo
        call MPI_Send(wporo, int_mpi, MPI_DOUBLE_PRECISION, &
                      irank,  tag_mpi_0, option%mycomm, ierr)
        int_mpi = 3*nlo
        call MPI_Send(wperm, int_mpi, MPI_DOUBLE_PRECISION, &
                      irank,  tag_mpi_0, option%mycomm, ierr)

        ! Free the work arrays

        deallocate(winat)
        deallocate(wporo)
        deallocate(wperm)

      endif
    endif

    ! Barrier call to stop other procs running ahead

    call MPI_Barrier(option%mycomm, ierr)
  enddo

end subroutine GridEclipsePorPermExchangeAndSet

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

  ierr = INPUT_ERROR_NONE
  qerr = PETSC_FALSE

!  Do the read operation

  if (len_trim(za) > 0) then
    read(za, *, iostat = ierr) iv
  endif

!  Check for errors

  if (InputError(ierr)) then
    call SetError(trim(zk)// &
      ' (error in reading integer value: '//trim(za)//')', qerr)
  else if (iv < il) then
    call SetError(trim(zk)//' (value '//trim(za)//' below lower bound)', qerr)
  else if (iv > iu) then
    call SetError(trim(zk)//' (value '//trim(za)//' above upper bound)', qerr)
  endif

end subroutine ProcessArgToInt

! *************************************************************************** !

subroutine GridEclipseProcnumExchangeAndSet(inatsend, nlmax, option)
  !
  ! Set up procnum array by getting the inat cell indices from each proc
  !
  ! Author: Dave Ponting
  ! Date: 02/21/19

  implicit none

  PetscInt , pointer :: inatsend(:)
  PetscInt, intent(in) :: nlmax
  type(option_type) :: option

  PetscInt :: irank, iorank, t_rank
  PetscInt :: ilo, ino
  PetscInt :: nlo, ierr
  PetscMPIInt :: int_mpi, temp_int_array(1), status_mpi(MPI_STATUS_SIZE)

  ! Work array

  PetscInt, allocatable :: winat(:)

  ! Set scalars

  ierr   = 0
  iorank = option%comm%io_rank
  t_rank = option%myrank

  ! Loop over exchange operations between IO rank and non-IO ranks (irank)

  do irank = 0, option%comm%size-1

    if (irank /= iorank) then

      ! Consider exchange between irank and iorank

      if (t_rank  == irank) then

        ! This is irank: send nlmax value and natural addresses

        temp_int_array(1) = nlmax
        call MPI_Send(temp_int_array, ONE_INTEGER_MPI, MPI_INTEGER, &
                      iorank,  tag_mpi_0, option%mycomm, ierr)

        int_mpi = nlmax
        call MPI_Send(inatsend      , int_mpi        , MPI_INTEGER, &
                      iorank,  tag_mpi_0, option%mycomm, ierr)

      endif

      if (t_rank  == iorank) then

       !This is IO rank - receive the other rank nlmax value (nlo) from irank

        temp_int_array(1) = 0
        call MPI_Recv(temp_int_array, ONE_INTEGER_MPI, MPI_INTEGER, &
                      irank, MPI_ANY_TAG, option%mycomm, status_mpi, ierr)
        nlo = temp_int_array(1)

        ! Allocate work array to hold the received natural addresses

        allocate(winat(nlo))

        ! Receive the natural address array from irank

        int_mpi = nlo
        call MPI_Recv(winat, int_mpi, MPI_INTEGER, irank, MPI_ANY_TAG, &
                      option%mycomm, status_mpi, ierr)

        ! On this (io) proc, set up the procnum values

        do ilo = 1, nlo
          ino = winat(ilo)
          call SetProcnumValue(ino,irank)
        enddo

        ! Free the inat arrays

        deallocate(winat)

      endif
    endif

    ! Barrier call to stop other procs running ahead

    call MPI_Barrier(option%mycomm, ierr)
  enddo

end subroutine GridEclipseProcnumExchangeAndSet

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

  ierr = INPUT_ERROR_NONE
  qerr = PETSC_FALSE

!  Read value

  if (len_trim(za) > 0) then
    read(za, *, iostat = ierr) rv
  endif

!  Check for errors

  if (InputError(ierr)) then
    call SetError(trim(zk)//' (error in reading real value: '//trim(za)//')',qerr)
  endif

end subroutine ProcessArgToReal

! *************************************************************************** !

subroutine ProcessGridData(qerr,option)
  !
  ! Process the data found in the grdecl file
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18
  !

  use Grid_Eclipse_Util_module, only : SetLGRName, SetLGRDimensions

  implicit none

  PetscBool, intent(inout) :: qerr
  type(option_type) :: option

  PetscInt  :: ix, iy, iz, izp, ize, ig, ia, iact, ic, igt, &
               nx, ny, nz, nxyzt, iof, ilgr
  PetscReal :: porv, bvol, eps
  PetscReal :: x000(3), x100(3), x010(3), x110(3), &
               x001(3), x101(3), x011(3), x111(3)
  PetscReal, allocatable :: bvg(:)

! Set up any unset optional grid arrays

  if(g_nlgr>1) then
    call SetupOptionalGridArrays()
  endif

!  Set the local grid data in grdecl_utilities

  call SetNLGR(g_nlgr)
  do ilgr=1,g_nlgr
    call SetLGRName      (ilgr,g_g(ilgr)%ref_name)
    call SetLGRDimensions(ilgr,g_g(ilgr)%nx, &
                               g_g(ilgr)%ny, &
                               g_g(ilgr)%nz)
  enddo

!  Process faults in base grid

  ! call ProcessFaults(g_ifld)

!--Set up geometry for the base grid-------------------------------------------

  ! Pore volume for rock-filled cell

  eps = 1.0e-6

  ! Fill in missing data values

  if (g_iscpg) then

    ! If corner point data supplied, extract cell dimensions like dx

    call extractCellDimensionsAndLocationsFromCPG(g_ifld)
  else

    ! Fill in x and y locations for each z-layer in field grid

    do iz = 1, g_g(g_ifld)%nz
     call FillXYPositionsForLayer(iz,g_ifld)
    enddo

    ! Fill in tops and z locations for each ix/iy column in field grid

    do ix = 1, g_g(g_ifld)%nx
      do iy = 1, g_g(g_ifld)%ny
        call FillZPositionsForColumn(ix, iy, g_ifld)
      enddo
    enddo

    ! Allocate and fill in coord and zcorn

    call CheckAllocCPGA(g_ifld)

    call ExtractCPGFromCellDimensionsAndLocations(g_ifld)

  endif

!--Set up geometry for the local grids-------------------------------------------

  ! if (g_nlgr>1) then
  !   call DoLGRSetup()
  ! endif

  ! Do Dykstra-Parsons if required

  ! if (g_isdpcf) then
  !   call applyDPCF(g_ifld)
  ! endif

!--Find pore volumes and set up active order---------------------------------------

  print *, 'Calculating pore volumes'

  ! Set up bulk volume by grid index work array

  g_na = 0
  g_nxyzt = 0
  do ilgr = 1,g_nlgr

    nx = g_g(ilgr)%nx
    ny = g_g(ilgr)%ny
    nz = g_g(ilgr)%nz

    g_grid_iof(ilgr) = g_nxyzt

    g_nxyzt = g_nxyzt + nx*ny*nz
  enddo

  allocate(bvg(g_nxyzt))
  bvg =  0.d0

!--Do calculations over active grids------------------------------------------------

  nxyzt = 0
  do ilgr = 1,g_nlgr
    do ix = 1, g_g(ilgr)%nx
      do iy = 1, g_g(ilgr)%ny
        do iz = 1, g_g(ilgr)%nz

          ig  = GetNaturalIndex(ix, iy, iz, ilgr)
          igt = ig+nxyzt

          ! Find pore volume

          bvol = 0.0
          porv = 0.0

          if (g_iscpg) then
            call GetCorners(ix, iy, iz, &
                             x000, x100, x010, x110, &
                             x001, x101, x011, x111, &
                             g_g(ilgr)%coord, g_g(ilgr)%zcorn, g_g(ilgr)%nx, g_g(ilgr)%ny)
            bvg(igt) = FindVolume(x000, x100, x010, x110, &
                                   x001, x101, x011, x111)
          else
            bvg(igt) = g_g(ilgr)%dx(ig)* &
                       g_g(ilgr)%dy(ig)* &
                       g_g(ilgr)%dz(ig)
          endif

          bvg(igt) = bvg(igt)*g_g(ilgr)%mv (ig)* &
                              g_g(ilgr)%ntg(ig)

          bvol = bvg(igt)

          ! iact is the ACTNUM array value

          iact = g_g(ilgr)%actn(ig)

          ! Find pore volume as poro*bvol

          if (iact == 1) then
            porv = g_g(ilgr)%poro(ig)*bvol
          endif
          if (iact == 2) then ! Rock volume only, tiny pore vol, no fluid flow
             g_g(ilgr)%poro(ig) = eps
             porv       = eps*bvol
             g_g(ilgr)%kx(ig)   = 0.0
             g_g(ilgr)%ky(ig)   = 0.0
             g_g(ilgr)%kz(ig)   = 0.0
          endif
          if (iact == 3) then ! Pore volume only, but ntg honoured
             g_g(ilgr)%poro(ig) = 1.0
             porv              = bvol
          endif
          if(iact == 4) then
            bvol = 1.0D0
            porv = g_g(ilgr)%poro(ig)*bvol
          endif

          if ((porv > g_minpv(1))  &
              .or. ((iact == 2) .and. (bvol > eps))) then
            g_na = g_na+1
            g_g(ilgr)%gtoa(ig) = g_na
            if(porv<g_pvfloor(1)) bvg(igt)=bvg(igt)*g_pvfloor(1)/porv
          endif

        enddo
      enddo
    enddo

!  Increment the count of cells in previous local grids

    nxyzt = nxyzt + g_g(ilgr)%nx*g_g(ilgr)%ny*g_g(ilgr)%nz
  enddo

!--Only continue if some active cells------------------------------------------

  if (g_na>0) then

    allocate(g_atog(g_na))
    allocate(g_atol(g_na))
    allocate(g_atoc(g_na))

    allocate(g_atox(g_na))
    allocate(g_atoy(g_na))
    allocate(g_atoz(g_na))

    g_atog = -1
    g_atol = -1
    g_atoc = -1
    g_atox = -1
    g_atoy = -1
    g_atoz = -1

    ! Set up a to g

    do ilgr = 1,g_nlgr

      iof= g_grid_iof(ilgr)

      do ix=1,g_g(ilgr)%nx
        do iy=1,g_g(ilgr)%ny
          do iz=1,g_g(ilgr)%nz

            ig = GetNaturalIndex(ix, iy, iz, ilgr)
            ia = g_g(ilgr)%gtoa(ig)
            if (ia>-1) then

              g_atog(ia) = ig
              g_atol(ia) = ilgr

              g_atox(ia) = ix
              g_atoy(ia) = iy
              g_atoz(ia) = iz

            endif

          enddo
        enddo
      enddo

    enddo

    ! Set up atoc (loop over cells in Eclipse compressed natural order)

    ic = 0
    do ilgr = 1,g_nlgr
      g_cnat_iof(ilgr) = ic
      do ize = 1, g_g(ilgr)%nz
        izp = g_g(ilgr)%nz-ize+1
        do iy = 1, g_g(ilgr)%ny
          do ix = 1, g_g(ilgr)%nx
            ig = GetNaturalIndex(ix, iy, izp, ilgr)
            ia = g_g(ilgr)%gtoa(ig)
            if (ia>-1) then
              ic = ic+1
              g_atoc(ia) = ic
            else
              if(g_nlgr>1) then
                if(ilgr == g_ifld) then
                  if(g_g(ilgr)%iglob(ig)>1) then
                    ic = ic + 1
                  endif
                endif
              endif
            endif
          enddo
        enddo
      enddo
      g_nac = ic
    enddo

    ! Initial estimate for number of connections
    ! Sufficient for normal connections + some space for wells

    g_mc = 3*g_na
    g_nc = 0

    ! Allocate active arrays

    call AllocateActiveArrays()

    ! Find the cell locations and connections

    call GenerateGridConnections(bvg,option)

    ! Connections between lgrs

    if (g_nlgr>1) then
      call GenerateLGRConnections(option)
    endif

    call ReportGridConnections()

  endif

  !  Deallocate the bulk volume by grid indexc work array

  if (allocated(bvg)) deallocate(bvg)

end subroutine ProcessGridData

! *************************************************************************** !

subroutine ReadEclipseEgridFile(efilename, coord, zcorn, actn, nx, ny, nz, &
                                qerr, zmess)
  !
  ! Read an Eclipse egrid file
  !
  ! Author: Dave Ponting
  ! Date: 08/11/21

  use String_module,only : StringCompareIgnoreCase, StringToUpper
  use Grid_Eclipse_Util_module

  implicit none

  character(len = *), intent(in) :: efilename

  PetscReal,intent(out) :: coord(:)
  PetscReal,intent(out) :: zcorn(:)
  PetscInt ,intent(out) :: actn (:)
  PetscInt ,intent(in)  :: nx, ny, nz
  PetscBool,intent(out) :: qerr
  character(len = MAXSTRINGLENGTH),intent(out) :: zmess

  PetscInt :: nxp, nyp, nnx, nny, nnz, nxyz, ncoord, nzcorn

  PetscInt :: nxy, ix, iy, izp, ize, ige, igp, i
  PetscInt :: inx, iny, inze, inzp, ine, inp, ixp, iyp, ibase
  PetscInt :: nxfile, nyfile, nzfile

  PetscInt, allocatable           :: vfilehead   (:)
  PetscInt, allocatable           :: vgridhead   (:)
  character(len = 8), allocatable :: vgridunit   (:)
  PetscReal, allocatable          :: vmapaxes    (:)
  character(len = 8), allocatable :: vmapaxesunit(:)
  character(len = 8), allocatable :: vgdorient   (:)
  PetscReal,allocatable           :: vcoord      (:)
  PetscReal,allocatable           :: vzcorn      (:)
  PetscInt ,allocatable           :: vactn       (:)

  character(len = MAXSTRINGLENGTH), allocatable :: vgdorientword(:)

  PetscInt :: ios, isuf_g

  PetscReal :: conv
  PetscBool :: eof,block_read
  PetscInt  :: itype, n, efile_type
  character(len = MAXSTRINGLENGTH) :: efilenameu

  character(len = 8) :: zmnem8, units

  ! Initialise

  ios        = 0
  qerr       = PETSC_FALSE
  zmess      = 'OK'
  conv       = 1.0D0
  efile_type = 0

  ! Basic dimensions

  nxp = nx+1
  nyp = ny+1

  nnx = 2*nx
  nny = 2*ny
  nnz = 2*nz

  nxy    = nx*ny
  nxyz   = nx*ny*nz
  ncoord = nxp*nyp*6
  nzcorn = nnx*nny*nnz

  !  Check file name

  efilenameu = efilename
  call StringToUpper(efilenameu)
  isuf_g = index(efilenameu,'.GRID')

  if(isuf_g > 0) then

    zmess = 'GDFILE only supports egrid files, not grid files '//trim(efilename)
    qerr  = PETSC_TRUE

  else

    !  Open files on the io proc

    call OpenEgridFileToRead(efilename, ios)
    if(ios/=0) then
      zmess = 'Unable to read egrid file '//trim(efilename)
      qerr  = PETSC_TRUE
    endif

  endif

  if(.not. qerr) then

    ! Initialise and set up useful scalars

    itype   = 0
    zmnem8  = ' '

    ! Read on egrid input unit (and only on I/O proc)

    e_fileunit = UNIT_EGRID_READ

    eof = PETSC_FALSE

    do while(.not.eof)

      call ReadHeader(itype,zmnem8, n, eof)
      if (eof) zmnem8 = 'END'
      call StringToUpper(zmnem8)

      ! Exit if end of file

      if (trim(zmnem8) == 'END' .or. qerr) exit

      block_read=PETSC_FALSE

      ! File header

      if (trim(zmnem8) == 'FILEHEAD') then
        allocate(vfilehead(n))
        call ReadBlockI(vfilehead,n)
        if(n>=5) then
          efile_type = vfilehead(5)
          if(efile_type>0) then
            qerr  = PETSC_TRUE
            zmess = 'Cannot read unstructured type egrid files'
          endif
        endif
        deallocate(vfilehead)
        block_read=PETSC_TRUE
      endif

      ! Map units

      if (trim(zmnem8) == 'MAPUNITS') then
        allocate(vmapaxesunit(n))
        call ReadBlockC(vmapaxesunit,n)
        deallocate(vmapaxesunit)
        block_read=PETSC_TRUE
      endif

      ! Map axes

      if (trim(zmnem8) == 'MAPAXES') then
        allocate(vmapaxes(n))
        call ReadBlockS(vmapaxes,n)
        if(n>=6) call SetMapAxes(vmapaxes)
        deallocate(vmapaxes)
        block_read=PETSC_TRUE
      endif

      ! Grid header

      if (trim(zmnem8) == 'GRIDHEAD') then
        allocate(vgridhead(n))
        call ReadBlockI(vgridhead,n)
        if(n>=4) then

          efile_type = vgridhead(1)
          if(efile_type/=1) then
            qerr  = PETSC_TRUE
            zmess = 'Cannot read unstructured type egrid files'
          endif

          nxfile = vgridhead(2)
          nyfile = vgridhead(3)
          nzfile = vgridhead(4)

          call CheckForDimensionErr('X',nx,nxfile,zmess,qerr)
          call CheckForDimensionErr('Y',ny,nyfile,zmess,qerr)
          call CheckForDimensionErr('Z',nz,nzfile,zmess,qerr)

        endif
        deallocate(vgridhead)
        block_read=PETSC_TRUE
      endif

      ! Grid units

      if (trim(zmnem8) == 'GRIDUNIT') then
        allocate(vgridunit(n))
        call ReadBlockC(vgridunit,n)
        units = vgridunit(1)
        if(StringCompareIgnoreCase(units,'feet')) then
          conv=0.3048D0
        else
          if(.not. StringCompareIgnoreCase(units,'metres')) then
            zmess = 'GRIDUNIT units must be metres or feet'
            qerr  = PETSC_TRUE
          endif
        endif
        deallocate(vgridunit)
        block_read=PETSC_TRUE
      endif

    ! Grid orientation

      if (trim(zmnem8) == 'GDORIENT') then
        allocate(vgdorient    (n))
        allocate(vgdorientword(n))
        call ReadBlockC(vgdorient,n)
        do i=1,n
          vgdorientword(i) = vgdorient(i)
        enddo
        call SetGDORIENT(vgdorientword,qerr,zmess)
        deallocate(vgdorient)
        deallocate(vgdorientword)
        block_read=PETSC_TRUE
      endif

     ! COORD values

      if (trim(zmnem8) == 'COORD') then

        if(n/=ncoord) then
          zmess = 'COORD array not of expected size'
          qerr  = PETSC_TRUE
        endif

        if(.not.qerr) then

          allocate(vcoord(ncoord))

          call ReadBlockS(vcoord,ncoord)

          do iyp = 1, nyp
            do ixp = 1, nxp

              ibase = 6*(nxp*(iyp-1) + ixp-1)

              coord(ibase+4) =  conv*vcoord(ibase+1)
              coord(ibase+5) =  conv*vcoord(ibase+2)
              coord(ibase+6) = -conv*vcoord(ibase+3)

              coord(ibase+1) =  conv*vcoord(ibase+4)
              coord(ibase+2) =  conv*vcoord(ibase+5)
              coord(ibase+3) = -conv*vcoord(ibase+6)

            enddo
          enddo

          deallocate(vcoord)

          block_read=PETSC_TRUE

       endif

     endif

    ! ZCORN values

     if (trim(zmnem8) == 'ZCORN') then

       if(n/=nzcorn) then
         zmess = 'ZCORN array not of expected size'
         qerr  = PETSC_TRUE
       endif

       if(.not.qerr) then

         allocate(vzcorn(nzcorn))

         call ReadBlockS(vzcorn,nzcorn)

         do inze = 1, nnz
           inzp = nnz-inze+1 ! Get the Pflotran k-up index
           do iny = 1, nny
             do inx = 1, nnx

               ine = (inze-1)*nnx*nny + (iny-1)*nnx + inx
               inp = (inzp-1)*nnx*nny + (iny-1)*nnx + inx

               zcorn(inp) = -conv*vzcorn(ine)

             enddo
           enddo
         enddo

         deallocate(vzcorn)

         block_read=PETSC_TRUE

       endif

     endif

     ! Active cell indicator

     if (trim(zmnem8) == 'ACTNUM') then

       if(n/=nxyz) then
         zmess = 'ACTNUM array not of expected size'
         qerr  = PETSC_TRUE
       endif

       if(.not.qerr) then

         allocate(vactn(nxyz))

         vactn = 0

         call ReadBlockI(vactn,nxyz)

         do ize = 1, nz
           izp = nz-ize+1 ! Get the Pflotran k-up index
           do iy = 1, ny
             do ix = 1, nx

               igp = nxy*(izp-1) + nx*(iy-1) + ix
               ige = nxy*(ize-1) + nx*(iy-1) + ix
               actn(igp) = vactn(ige)

              enddo
            enddo
          enddo

          deallocate(vactn)

          block_read=PETSC_TRUE

        endif

      endif

      !  Read dummy values if block not required

      if(.not.block_read) then
        if (itype == e_typeI) call readDummyI(n)
        if (itype == e_typeS) call readDummyS(n)
        if (itype == e_typeD) call readDummyD(n)
        if (itype == e_typeB) call readDummyB(n)
        if (itype == e_typeC) call readDummyC(n)
      endif

    enddo

  endif ! End of able to open file check

  ! Close file and reset formatted read flag

  close(unit = UNIT_EGRID_READ, iostat=ios)

  e_formatted_r = PETSC_FALSE

end subroutine ReadEclipseEgridFile

! *************************************************************************** !

subroutine ReadEZcornArray(a, input, option, convfact, qerr, ilgr)
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
  PetscReal, intent(in) :: convfact
  PetscBool, intent(inout) :: qerr
  PetscInt , intent(in) :: ilgr

  PetscBool, parameter :: nn_no = PETSC_FALSE

  PetscReal, allocatable :: column_buffer(:)
  PetscInt  :: inx , iny, inz, nnx, nny, nnz, nnxy, inode, nread, inzpft

  qerr = PETSC_TRUE

  ! Do the actual read operation

  nread = size(a)
  if (nread /= 8*g_g(ilgr)%nxyz) then
    call SetError('Zcorn read', qerr)
  else
    call ReadEvalues(a, nread,'ZCORN','GRID',input,option,convfact,nn_no,qerr)
  endif

  if (.not.qerr) then

    nnx = 2*g_g(ilgr)%nx
    nny = 2*g_g(ilgr)%ny
    nnz = 2*g_g(ilgr)%nz
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

subroutine readBlockI(a,n)
  !
  ! Read an integer record
  !
  ! Author: Dave Ponting
  ! Date: 10/08/19

  implicit none

  PetscInt,intent(out)::a(:)
  PetscInt,intent(in) ::n

  PetscInt, parameter :: blksize  = 6    ! Values/line (formatted)
  PetscInt, parameter :: mrecsize = 1000 ! Values/rec (unformatted)

  PetscInt :: il, iu, j, irec, nrec, ninrec, ios

  PetscInt, allocatable :: ibuf(:)

  !  Initialise

  a = 0

  ! Read n values

  if (e_formatted_r) then

  ! Formatted case: write out in lines of 6 values per line, I11 format

    do il = 1, n, blksize
      iu = min(il+blksize-1, n)
      read(e_fileunit, '(6(1x,i11))') (a(j), j = il, iu)
    enddo

  else

  ! Allocate an integer buffer

    allocate(ibuf(mrecsize))
    ibuf = 0

    ! Loop over the number of records required

    nrec = GetNumberOfRecords(n, mrecsize)
    do irec = 1, nrec
      ! Get record details
      call GetRecordDetails(irec, n, mrecsize, il, iu, ninrec)
      ! Read the record
      read(e_fileunit,iostat=ios) (ibuf(j), j = 1, ninrec)
      ! Copy from buffer
      call CopyFromBufferI(a, ibuf, il, iu)
    enddo

    ! Delete the buffer

    deallocate(ibuf)
  endif

end subroutine readBlockI

! *************************************************************************** !

subroutine readBlockS(a,n)
  !
  ! Read a single precision record
  !
  ! Author: Dave Ponting
  ! Date: 10/08/19

  implicit none

  PetscReal,intent(out)::a(:)
  PetscInt,intent(in) ::n

  PetscInt, parameter :: blksize = 4 ! Values/line (formatted)
  PetscInt, parameter :: mrecsize = 1000 ! Values/rec (unformatted)

  PetscInt :: il, iu, j, irec, nrec, ninrec, ios
  PetscReal, allocatable :: fbuf(:)

  a = 0.0

  ! Read n values

  if (e_formatted_r) then

    ! Formatted case: read in lines of 6 values per line, I11 format

    do il = 1, n, blksize
      iu = min(il+blksize-1, n)
      read(e_fileunit, '(4(1x,e16.8))') (a(j), j = il, iu)
    enddo

  else

    ! Allocate a real32 buffer

    allocate(fbuf(mrecsize))
    fbuf = 0.0

    ! Loop over the number of records required

    nrec = GetNumberOfRecords(n, mrecsize)
    do irec = 1, nrec
      ! Get record details
      call GetRecordDetails(irec, n, mrecsize, il, iu, ninrec)
      ! Read the record
      read(e_fileunit,iostat=ios) (fbuf(j), j = 1, ninrec)
      ! Copy from buffer
      call CopyFromBufferS(a, fbuf, il, iu)
    enddo

    ! Delete the buffer

    deallocate(fbuf)
  endif

end subroutine readBlockS

! *************************************************************************** !

subroutine readBlockD(a,n)
  !
  ! Read a double precision record
  !
  ! Author: Dave Ponting
  ! Date: 10/08/19

  implicit none

  PetscReal,intent(out)::a(:)
  PetscInt,intent(in) ::n

  PetscInt, parameter :: blksize = 3 ! Values/line (formatted)
  PetscInt, parameter :: mrecsize = 1000 ! Values/rec (unformatted)

  PetscInt :: il, iu, j, irec, nrec, ninrec, ios
  PetscReal, allocatable :: dbuf(:)

  a = 0.0

  ! Read n values

  if (e_formatted_r) then

    ! Formatted case: read in lines of 6 values per line, I11 format

    do il = 1, n, blksize
      iu = min(il+blksize-1, n)
      read(e_fileunit, '(3(1x,e22.14))') (a(j), j = il,iu)
    enddo

  else

    ! Allocate a double precision buffer

    allocate(dbuf(mrecsize))
    dbuf = 0.0

    ! Loop over the number of records required

    nrec = GetNumberOfRecords(n, mrecsize)
    do irec = 1, nrec
      ! Get record details
      call GetRecordDetails(irec, n, mrecsize, il, iu, ninrec)
      ! Read the record
      read(e_fileunit,iostat=ios) (dbuf(j), j = 1, ninrec)
      ! Copy from buffer
      call CopyFromBufferD(a, dbuf, il, iu)
    enddo

    ! Delete the buffer

    deallocate(dbuf)

  endif

end subroutine readBlockD

! *************************************************************************** !

subroutine readBlockB(a,n)
  !
  ! Read a boolean record
  !
  ! Author: Dave Ponting
  ! Date: 10/08/19

  implicit none

  PetscBool,intent(out)::a(:)
  PetscInt,intent(in) ::n

  PetscInt, parameter :: blksize = 25 ! Values/line (formatted)
  PetscInt, parameter :: mrecsize = 1000 ! Values/rec (unformatted)

  PetscInt :: il, iu, j, irec, nrec, ninrec, ios
  PetscInt, allocatable :: ibuf(:)

  a = PETSC_FALSE

  ! Read n bool values

  if (e_formatted_r) then

    ! Formatted case: read in lines of 6 values per line, I11 format

    do il = 1, n, blksize
      iu = min(il+blksize-1, n)
      read(e_fileunit, '(25(1X,L2))') (a(j), j = il,iu)
    enddo

  else

    ! Allocate an int32 buffer

    allocate(ibuf(mrecsize))
    ibuf = 0

    ! Loop over the number of records required

    nrec = GetNumberOfRecords(n, mrecsize)
    do irec = 1, nrec
      ! Get record details
      call GetRecordDetails(irec, n, mrecsize, il, iu, ninrec)
      ! Read the record
      read(e_fileunit,iostat=ios) (ibuf(j), j = 1, ninrec)
      ! Copy from buffer
      call CopyFromBufferB(a, ibuf, il, iu)
    enddo

    ! Delete the buffer

    deallocate(ibuf)
  endif

end subroutine readBlockB

! *************************************************************************** !

subroutine readBlockC(a,n)
  !
  ! Read a char*8 record
  !
  ! Author: Dave Ponting
  ! Date: 10/08/19

  implicit none

  character(len=8),intent(out)::a(:)
  PetscInt,intent(in) ::n

  PetscInt, parameter :: blksize  = 7   ! Values/line (formatted)
  PetscInt, parameter :: mrecsize = 105 ! Values/rec  (unformatted)

  PetscInt :: il, iu, j, irec, nrec, ninrec, ios

  a=' '

10 format(7(1X, "'", A8, "'"))

  ! Read n integer values

  if (e_formatted_r) then

  ! Formatted case: read in lines of 7 values per line, A8 format

    do il = 1, n, blksize
      iu = min(il+blksize-1, n)
      read(e_fileunit, 10) (a(j), j = il,iu)
    enddo

  else

  ! Loop over the number of records required

    nrec = GetNumberOfRecords(n, mrecsize)
    do irec = 1, nrec
  ! Get record details
      call GetRecordDetails(irec, n, mrecsize, il, iu, ninrec)
  ! Read the record
      read(e_fileunit,iostat=ios) (a(j), j = il, iu)
    enddo

  endif

end subroutine readBlockC

! *************************************************************************** !

subroutine readDummyI(n)
  !
  ! Read a integer record but do not hold the values
  !
  ! Author: Dave Ponting
  ! Date: 10/08/19

  implicit none

  PetscInt, intent(in) :: n
  PetscInt, allocatable :: a(:)

  allocate(a(n))

  call readBlockI(a,n)

  deallocate(a)

end subroutine readDummyI

! *************************************************************************** !

subroutine readDummyS(n)
  !
  ! Read a single precision record but do not hold the values
  !
  ! Author: Dave Ponting
  ! Date: 10/08/19

  implicit none

  PetscInt, intent(in) :: n
  PetscReal, allocatable :: a(:)

  allocate(a(n))

  call readBlockS(a,n)

  deallocate(a)

end subroutine readDummyS

! *************************************************************************** !

subroutine readDummyD(n)
  !
  ! Read a double precision record but do not hold the values
  !
  ! Author: Dave Ponting
  ! Date: 10/08/19

  implicit none

  PetscInt, intent(in) :: n
  PetscReal, allocatable :: a(:)

  allocate(a(n))

  call readBlockD(a,n)

  deallocate(a)

end subroutine readDummyD

! *************************************************************************** !

subroutine readDummyB(n)
  !
  ! Read a boolean record but do not hold the values
  !
  ! Author: Dave Ponting
  ! Date: 10/08/19

  implicit none

  PetscInt, intent(in) :: n
  PetscBool, allocatable :: a(:)

  allocate(a(n))

  call readBlockB(a,n)

  deallocate(a)

end subroutine readDummyB

! *************************************************************************** !

subroutine readDummyC(n)
  !
  ! Read a single precision record but do not hold the values
  !
  ! Author: Dave Ponting
  ! Date: 10/08/19

  implicit none

  PetscInt, intent(in) :: n
  character(len = 8),allocatable :: a(:)

  allocate(a(n))

  call readBlockC(a,n)

  deallocate(a)

end subroutine readDummyC

! *************************************************************************** !

subroutine ReadHeader(itype, zmnem8, n, eof)
  !
  ! Read a block header
  !
  ! Author: Dave Ponting
  ! Date: 12/15/18

  implicit none

  PetscInt, intent(out) :: itype
  character(len = 8), intent(out) :: zmnem8
  PetscInt, intent(out) :: n
  PetscBool, intent(out) :: eof
  character(len = 4) :: ztype4

  PetscInt :: n4
  PetscInt :: ios

10 format(1X, "'", A8, "'", 1X, I11, 1X, "'", A4, "'")

  zmnem8 = '    '
  n  = 0
  n4 = 0
  ztype4 = '    '
  eof    = PETSC_FALSE

  !  Do the read

  if (e_formatted_r) then
    read(e_fileunit, 10,iostat=ios) zmnem8, n , ztype4
  else
    read(e_fileunit    ,iostat=ios) zmnem8, n4, ztype4
    n = n4
  endif

  if (ios == 0) then

  !  Set up output values

    if (ztype4 == 'INTE') itype = e_typeI
    if (ztype4 == 'REAL') itype = e_typeS
    if (ztype4 == 'DOUB') itype = e_typeD
    if (ztype4 == 'LOGI') itype = e_typeB
    if (ztype4 == 'CHAR') itype = e_typeC

  else
    eof = PETSC_TRUE
  endif

end subroutine ReadHeader

! *************************************************************************** !

subroutine ReadECoordArray(a, input, option, convfact, qerr, ilgr)
  !
  ! Reads an Eclgrid coord array
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  PetscReal, intent(inout) :: a(:)
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscReal, intent(in) :: convfact
  PetscBool, intent(inout) :: qerr
  PetscInt  ,intent(in) :: ilgr
  PetscBool, parameter :: nn_no = PETSC_FALSE

  PetscInt  :: i, nread

  nread = size(a)
  if (nread /= 6*g_g(ilgr)%nxpnyp) then
    call SetError('Coord read', qerr)
  else
    call ReadEvalues(a,nread,'COORD','GRID',input,option,convfact,nn_no,qerr)
  endif

  do i = 3, nread, 3
    a(i) = z_flip*a(i)
  enddo

end subroutine ReadECoordArray

! *************************************************************************** !

subroutine ReadEGridArrR(a, keyword ,&
                         input, option, convfact, is_dep, is_perm, is_nn, qerr, ilgr)
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
  PetscReal, intent(in) :: convfact
  PetscBool, intent(in) :: is_dep
  PetscBool, intent(in) :: is_perm
  PetscBool, intent(in) :: is_nn
  PetscBool, intent(inout) :: qerr
  PetscInt , intent(in) :: ilgr

  PetscReal, allocatable :: column_buffer(:)
  PetscInt  :: ix, iy, iz, izpft, ig, igpft, nread
  PetscReal :: flip, conv, v

  qerr = PETSC_FALSE

  ! Check DIMENS read

  call CheckDimensRead(qerr)

  if (.not.qerr) then

    ! Do the actual read operation

    nread = size(a)
    if (nread /= g_g(ilgr)%nxyz) then
      call SetError(keyword, qerr)
    else
      call ReadEvalues(a,nread,keyword,'GRID',input,option,convfact,is_nn,qerr)
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

    allocate(column_buffer(g_g(ilgr)%nz))

    ! Re-order the input values (Eclipse goes down the columns, PFT goes up)

    do ix = 1, g_g(ilgr)%nx
      do iy = 1, g_g(ilgr)%ny

        ! Load a column of Eclipse values

        do iz = 1, g_g(ilgr)%nz
          ig = GetNaturalIndex(ix, iy, iz, ilgr)
          column_buffer(iz) = a(ig)
        enddo

        ! Place into Pflotran-conversion array with
        ! appropriate data modifications

        do iz = 1, g_g(ilgr)%nz
          izpft = g_g(ilgr)%nz - iz + 1
          igpft = GetNaturalIndex(ix, iy, izpft, ilgr)
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

subroutine ReadEGridArrI(a, keyword, input, option, is_nn, qerr, ilgr)
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
  PetscInt, intent(in) :: ilgr

  PetscInt , allocatable :: column_buffer(:)
  PetscReal, allocatable :: ar(:)
  PetscInt  :: ix, iy, iz, izpft, ig, igpft, i
  PetscReal, parameter :: convnull=1.0

  ! Do the actual read operation

  allocate(ar(g_g(ilgr)%nxyz))
  ar = 1.0

  call ReadEvalues(ar,g_g(ilgr)%nxyz,keyword,'GRID',input,option,convnull,is_nn,qerr)

  do i = 1, g_g(ilgr)%nxyz
    a(i) = nint(ar(i),kind(g_petsc_int))
  enddo

  ! Allocate a column buffer to reorder to bottom-up convention

  allocate(column_buffer(g_g(ilgr)%nz))

  ! Re-order the input values (Eclipse goes down the columns, PFT goes up)

  do ix = 1, g_g(ilgr)%nx
    do iy = 1, g_g(ilgr)%ny

      ! Load a column of Eclipse values

      do iz = 1, g_g(ilgr)%nz
        ig = GetNaturalIndex(ix, iy, iz, ilgr)
        column_buffer(iz) = a(ig)
      enddo

      ! Place into Pflotran-converion array with appropriate data modifications

      do iz = 1, g_g(ilgr)%nz
        izpft = g_g(ilgr)%nz - iz + 1
        igpft = GetNaturalIndex(ix, iy, izpft, ilgr)
        a(igpft) = column_buffer(iz)
      enddo
    enddo

  enddo

  ! Free the column buffer

  deallocate(column_buffer)
  deallocate(ar)

end subroutine ReadEGridArrI

! *************************************************************************** !

subroutine ReadEstrings(zkey, a, n, input, option, qerr)
  !
  ! Read a series of strings from an Eclipse syntax file
  !
  ! Author: Dave Ponting
  ! Date: 11/23/18

  implicit none

  character(len = MAXWORDLENGTH), intent(in)    :: zkey
  character(len = MAXSTRINGLENGTH), intent(inout) :: a(:)
  PetscInt , intent(in)    :: n
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscBool, intent(inout) :: qerr

  PetscInt  :: i, iostat, istar, nstack, m, nr, ierr
  PetscBool :: exittime
  character(len = MAXSTRINGLENGTH) :: word
  character(len = MAXSTRINGLENGTH) :: repc
  character(len = MAXSTRINGLENGTH) :: hold
  character(len = MAXSTRINGLENGTH) :: zval
  character(len = MAXSTRINGLENGTH) :: zmess

  ierr   = INPUT_ERROR_NONE
  qerr   = PETSC_FALSE

  i      = 0
  iostat = 0
  zval   = ' '
  exittime = PETSC_FALSE
  nstack = 0

  call InputReadPflotranStringNotComment(input, option)
  if (InputError(input)) then
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

  if (InputError(ierr)) then
    zmess = 'Unable to read ' // trim(zkey)
    call SetError(zmess, qerr)
  endif

end subroutine ReadEstrings

! *************************************************************************** !

subroutine ReadEvalues(a, n, keyword, section, input, option, &
                       convfact, is_nn, qerr)
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
  PetscReal, intent(in) :: convfact
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
  ierr   = INPUT_ERROR_NONE
  zmess = ' '

  call InputReadPflotranStringNotComment(input, option)
  if (InputCheckExit(input, option)) exittime = PETSC_TRUE

  if (.not. exittime) then

    g_column = 1
    do i = 1, n

      if (nstack > 0) then

        ! Case of value in stack

        a(i)   = dval*convfact
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
        if (InputError(ierr)) then
          call InputErrorMsg(input, option, keyword, section)
          ierr = INPUT_ERROR_DEFAULT
          exit
        else
          a(i) = dval*convfact
        endif

      endif

      if (is_nn) then
        if (a(i)*convfact < vmargin) then
          zmess = trim(keyword) //', value ' // trim(word) &
                                // ' cannot be negative'
          call SetError(zmess, qerr)
        endif
      endif

    enddo

  endif

  if (InputError(ierr)) then
    zmess = 'Unable to read ' // trim(keyword)
    call SetError(zmess, qerr)
  endif

end subroutine ReadEvalues

! *************************************************************************** !

subroutine ReallocateConnectionArrays()
  !
  ! Reallocate the connection arrays (used to extend these arrays)
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Utility_module, only: ReallocateArray

  implicit none

  PetscInt :: mc

  ! Note keep resetting mc as Reallocate keeps increasing it

  mc = g_mc
  call ReallocateArray(g_cia  , mc)
  mc = g_mc
  call ReallocateArray(g_cja  , mc)

  mc = g_mc
  call ReallocateArray(g_ccx  , mc)
  mc = g_mc
  call ReallocateArray(g_ccy  , mc)
  mc = g_mc
  call ReallocateArray(g_ccz  , mc)
  mc = g_mc
  call ReallocateArray(g_carea, mc)

  mc = g_mc
  call ReallocateArray(g_tran   , mc)
  mc = g_mc
  call ReallocateArray(g_tran_th, mc)
  mc = g_mc
  call ReallocateArray(g_idir   , mc)

  g_mc = mc

end subroutine ReAllocateConnectionArrays

! *************************************************************************** !

subroutine SetCoordLine(ixp, iyp, x, y, ilgr)
  !
  ! Set up a vertical coordinate line (corner of a pillar of cells)
  !
  ! Author: Dave Ponting
  ! Date: 11/13/18

  implicit none

  PetscInt , intent(in) :: ixp, iyp, ilgr
  PetscReal, intent(in) :: x, y

  PetscReal, parameter :: zl = -10000.0
  PetscReal, parameter :: zu =  10000.0

  PetscInt :: ibase

  ibase = 6*((iyp-1)*g_g(ilgr)%nxp + (ixp-1))

  g_g(ilgr)%coord(ibase+1) = x
  g_g(ilgr)%coord(ibase+2) = y
  g_g(ilgr)%coord(ibase+3) = zl
  g_g(ilgr)%coord(ibase+4) = x
  g_g(ilgr)%coord(ibase+5) = y
  g_g(ilgr)%coord(ibase+6) = zu

end subroutine SetCoordLine

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
  PetscInt :: nx, ny, nz

  nx = nint(dimens(1),kind(g_petsc_int))
  ny = nint(dimens(2),kind(g_petsc_int))
  nz = nint(dimens(3),kind(g_petsc_int))

  g_dimens_read = PETSC_TRUE

  call AllocateGridArrays(nx, ny, nz, g_ifld)

end subroutine SetDimens

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

subroutine SetProcnumValue(ia,irank)
  !
  ! Set the prcn value of a given active cell
  !
  ! Author: Dave Ponting
  ! Date: 03/10/22

  implicit none

  PetscInt, intent(in) :: ia
  PetscInt, intent(in) :: irank
  PetscInt :: ig, ilgr

  ! Get the grid order (ix-fastest)

  ig   = g_atog(ia)
  ilgr = g_atol(ia)

  ! Set value

  if (g_rprcn == 1) then
    g_g(ilgr)%prcn(ig) = irank
  endif

end subroutine SetProcnumValue

! *************************************************************************** !

subroutine SetupOptionalGridArrays()
  !
  ! Check that optional grids (drps, lgrs) set up over the LGR stack
  !
  ! Author: Dave Ponting
  ! Date: 05/10/22
  !
  implicit none

  PetscInt :: ilgr, nxyz

!--Loop over LGRs--------------------------------------------------------------

  do ilgr=1,g_nlgr

    nxyz = g_g(ilgr)%nxyz

  enddo

end subroutine SetupOptionalGridArrays

! *************************************************************************** !

subroutine SetVec3(vec3,x,y,z)
  !
  ! Load a 3-vector
  !
  ! Author: Dave Ponting
  ! Date: 05/05/20

  implicit none

  PetscReal,intent(out) :: vec3(g_ndir)
  PetscReal,intent(in) :: x, y, z

  vec3(g_xdir) = x
  vec3(g_ydir) = y
  vec3(g_zdir) = z

end subroutine SetVec3

! *************************************************************************** !

subroutine ReportGridConnections
  !
  ! Report all grid connections
  !
  ! Author: Dave Ponting
  ! Date: 08/13/21

  use Utility_module, only : DeallocateArray

  implicit none

  PetscInt :: ilgr, nxyz

  nxyz=0

  do ilgr=1,g_nlgr
    nxyz = nxyz + g_g(ilgr)%nx* &
                  g_g(ilgr)%ny* &
                  g_g(ilgr)%nz
  enddo

  print *, 'Ncell=', nxyz, ' Nact=', g_na, ' Nconn=', g_nc, &
           ' Nflt=', g_nf  , ' Npo=' , g_npo

end subroutine ReportGridConnections

! *************************************************************************** !

subroutine DeallocateGridArrays()
  !
  ! Deallocate the grid arrays
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Utility_module, only : DeallocateArray

  implicit none

  PetscInt :: ilgr

  do ilgr=1,g_nlgr

    if (g_g(ilgr)%cpgallocated) then
      call DeallocateArray(g_g(ilgr)%coord)
      call DeallocateArray(g_g(ilgr)%zcorn)
      g_g(ilgr)%cpgallocated = PETSC_FALSE
    endif

    if(ilgr>1) then

      call DeallocateArray(g_g(ilgr)%nlpgx)
      call DeallocateArray(g_g(ilgr)%nlpgy)
      call DeallocateArray(g_g(ilgr)%nlpgz)

      call DeallocateArray(g_g(ilgr)%ibpgx)
      call DeallocateArray(g_g(ilgr)%ibpgy)
      call DeallocateArray(g_g(ilgr)%ibpgz)

      call DeallocateArray(g_g(ilgr)%hrefx)
      call DeallocateArray(g_g(ilgr)%hrefy)
      call DeallocateArray(g_g(ilgr)%hrefz)

    endif

    call DeallocateArray(g_g(ilgr)%dx)
    call DeallocateArray(g_g(ilgr)%dy)
    call DeallocateArray(g_g(ilgr)%dz)

    call DeallocateArray(g_g(ilgr)%mv)

    call DeallocateArray(g_g(ilgr)%mx)
    call DeallocateArray(g_g(ilgr)%my)
    call DeallocateArray(g_g(ilgr)%mz)

    call DeallocateArray(g_g(ilgr)%mxn)
    call DeallocateArray(g_g(ilgr)%myn)
    call DeallocateArray(g_g(ilgr)%mzn)

    call DeallocateArray(g_g(ilgr)%tx)
    call DeallocateArray(g_g(ilgr)%ty)
    call DeallocateArray(g_g(ilgr)%tz)

    call DeallocateArray(g_g(ilgr)%tops)
    call DeallocateArray(g_g(ilgr)%ntg)

    call DeallocateArray(g_g(ilgr)%xloc)
    call DeallocateArray(g_g(ilgr)%yloc)
    call DeallocateArray(g_g(ilgr)%zloc)

    call DeallocateArray(g_g(ilgr)%actn)
    call DeallocateArray(g_g(ilgr)%satn)
    call DeallocateArray(g_g(ilgr)%imbn)
    call DeallocateArray(g_g(ilgr)%eqln)
    call DeallocateArray(g_g(ilgr)%mltn)
    call DeallocateArray(g_g(ilgr)%tbcn)
    call DeallocateArray(g_g(ilgr)%prcn)

    call DeallocateArray(g_g(ilgr)%gtoa)

    g_g(ilgr)%dx =>null()
    g_g(ilgr)%dy =>null()
    g_g(ilgr)%dz =>null()
    g_g(ilgr)%mv =>null()
    g_g(ilgr)%mx =>null()
    g_g(ilgr)%my =>null()
    g_g(ilgr)%mz =>null()
    g_g(ilgr)%mxn =>null()
    g_g(ilgr)%myn =>null()
    g_g(ilgr)%mzn =>null()
    g_g(ilgr)%tx  =>null()
    g_g(ilgr)%ty  =>null()
    g_g(ilgr)%tz  =>null()
    g_g(ilgr)%tops=>null()
    g_g(ilgr)%ntg =>null()
    g_g(ilgr)%xloc=>null()
    g_g(ilgr)%yloc=>null()
    g_g(ilgr)%zloc=>null()
    g_g(ilgr)%actn=>null()
    g_g(ilgr)%satn=>null()
    g_g(ilgr)%imbn=>null()
    g_g(ilgr)%eqln=>null()
    g_g(ilgr)%mltn=>null()
    g_g(ilgr)%tbcn=>null()
    g_g(ilgr)%prcn=>null()
    g_g(ilgr)%gtoa=>null()

    if(ilgr>1) then
      call DeallocateArray(g_g(ilgr)%ihost)
      g_g(ilgr)%ihost => null()
    endif

  enddo

  if(m_multregt > 0) then
    if(allocated(g_multregt_ir1)) deallocate(g_multregt_ir1)
    if(allocated(g_multregt_ir2)) deallocate(g_multregt_ir2)
    if(allocated(g_multregt_mult)) deallocate(g_multregt_mult)
    m_multregt = 0
  endif

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

  call DeallocateArray(g_x)
  call DeallocateArray(g_y)
  call DeallocateArray(g_z)

  call DeallocateArray(g_cia)
  call DeallocateArray(g_cja)
  call DeallocateArray(g_ccx)
  call DeallocateArray(g_ccy)
  call DeallocateArray(g_ccz)
  call DeallocateArray(g_carea)

  call DeallocateArray(g_tran)
  call DeallocateArray(g_tran_th)
  call DeallocateArray(g_idir)

  call DeallocateArray(g_atoc)

end subroutine DeallocateActiveArrays

! *************************************************************************** !

subroutine GridEcipseDeallocatePorPermArrays(solution_monitor,option)
  !
  ! Deallocate the perma and porosity arrays
  !
  ! Author: Dave Ponting
  ! Date: 12/11/18

  use Utility_module, only : DeallocateArray

  implicit none

  PetscBool, intent(in) :: solution_monitor
  type(option_type) :: option

  PetscInt :: ilgr

  if (option%myrank == option%comm%io_rank) then

    if (.not.solution_monitor) then
      call DeallocateArray(g_atog)
      call DeallocateArray(g_atol)
      g_atog_held = PETSC_FALSE
    else
      g_atog_held = PETSC_TRUE
    endif

    call DeallocateArray(g_atox)
    call DeallocateArray(g_atoy)
    call DeallocateArray(g_atoz)

    do ilgr=1,g_nlgr
      call DeallocateArray(g_g(ilgr)%kx)
      call DeallocateArray(g_g(ilgr)%ky)
      call DeallocateArray(g_g(ilgr)%kz)
      call DeallocateArray(g_g(ilgr)%poro)
    enddo

  endif

end subroutine GridEcipseDeallocatePorPermArrays

! *************************************************************************** !

end module Grid_Eclipse_module
