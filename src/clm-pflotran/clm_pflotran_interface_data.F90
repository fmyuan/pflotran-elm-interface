module clm_pflotran_interface_data
!
! NOTES for convenience:
!        (1) '*_pfp': mpi vecs for PF variables; '_clmp': mpi vecs for CLM variables;
!            '*_pfs': seq vecs for PF variables; '_clms': seq. vecs for CLM variables;
!        (2) '*_': 3D (XYZ) subsurface domain's variables;
!                  with '_x/y/z_' used for different directions of 3D domain.
!        (3) '*_subsurf_': 2D (XY) surface of 3D domain's varialbes;
!                         for bottom, uses '_subbase_', but essentially supposing same shape/area as surface.
!            '*_srf_': for variables with 2D-grids at ground surface, which may or may not same as '_subsurf_'. NOTE that it's not supported now.
! Revised by Fengming Yuan, CCSI-ORNL @May-2015

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
  use petscsys
  use petscvec

  implicit none

  private

  type, public :: clm_pflotran_idata_type

  ! Time invariant data:

  ! num of CLM soil layers that are mapped to/from PFLOTRAN (global constants, not local copy)
  PetscInt :: nzclm_mapped
  PetscInt :: nxclm_mapped
  PetscInt :: nyclm_mapped
  PetscReal :: x0clm_global
  PetscReal :: y0clm_global
  PetscReal :: z0clm_global
  PetscReal, pointer :: dxclm_global(:)              ! this is NOT the 3-D vec 'dxsoil' defined below, rather it's the universal x-direction interval (OR, longitudal degree interval from CLM land surf grids) for all gridcells
  PetscReal, pointer :: dyclm_global(:)              ! this is NOT the 3-D vec 'dysoil' defined below, rather it's the universal y-direction interval (OR, longitudal degree interval from CLM land surf grids)
  PetscReal, pointer :: dzclm_global(:)              ! this is NOT the 3-D vec 'dzsoil' defined below, rather it's the universal soil layer thickness (unit: m) for all gridcells

  ! decompose domain in 3-D (only work with structured PF grid currently)
  ! processors no.
  PetscInt :: npx, npy, npz
  ! domain nodes no. for each processors
  PetscInt, pointer :: clm_lx(:)   ! array size is 'npx'
  PetscInt, pointer :: clm_ly(:)   ! array size is 'npy'
  PetscInt, pointer :: clm_lz(:)   ! array size is 'npz'

  !-------------------------------------------------------------
  ! Number of cells for the 3D subsurface domain
  PetscInt :: nlclm_sub ! num of local clm cells
  PetscInt :: ngclm_sub ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_sub  ! num of local pflotran cells
  PetscInt :: ngpf_sub  ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the top/bottom cells of the 3D subsurface domain
  PetscInt :: nlclm_2dtop  ! num of local clm cells
  PetscInt :: ngclm_2dtop  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2dtop   ! num of local pflotran cells
  PetscInt :: ngpf_2dtop   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  PetscInt :: nlclm_2dbot  ! num of local clm cells
  PetscInt :: ngclm_2dbot  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2dbot   ! num of local pflotran cells
  PetscInt :: ngpf_2dbot   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the 2D surface domain
  PetscInt :: nlclm_srf  ! num of local clm cells
  PetscInt :: ngclm_srf  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_srf   ! num of local pflotran cells
  PetscInt :: ngpf_srf   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Mesh property

  ! sub-surf/sub-base area (in 2D): the toppest/lowest cells of subsurface domain for BCs
  ! At this moment, assumes both surf/base cells are exactly SAME in area.
  Vec :: area_subsurf_clmp ! mpi vec
  Vec :: area_subsurf_pfs  ! seq vec
  Vec :: area_subsurf_pfp  ! mpi vec
  Vec :: area_subsurf_clms ! seq vec

  ! Area of top face of PF cells (in 3D) (note: a PF cell has faces (facets) of top/bottom/east/west/south/north)
  Vec :: area_top_face_clmp ! mpi vec
  Vec :: area_top_face_pfs  ! seq vec
  Vec :: area_top_face_pfp  ! mpi vec
  Vec :: area_top_face_clms ! seq vec

  ! z-axis (in 3D), soil depth of center of a soil cell in unit of meters
  Vec :: zsoil_clmp          ! mpi vec
  Vec :: zsoil_pfs           ! seq vec
  Vec :: zsoil_pfp           ! mpi vec
  Vec :: zsoil_clms          ! seq vec

  ! x/y-axis (in 3D), grid center of a soil cell in unit of meters
  Vec :: xsoil_clmp          ! mpi vec
  Vec :: xsoil_pfs           ! seq vec
  Vec :: ysoil_clmp          ! mpi vec
  Vec :: ysoil_pfs           ! seq vec

  ! soil cell inter-nodes coordinates ('vertex' called in PF mesh; 'interface level' called in CLM soil layers)
  Vec :: zisoil_clmp          ! mpi vec
  Vec :: zisoil_pfs           ! seq vec

  ! length/width/thickness of soil cells (in 3D) in unit of meters
  Vec :: dxsoil_clmp          ! mpi vec
  Vec :: dxsoil_pfs           ! seq vec
  Vec :: dysoil_clmp          ! mpi vec
  Vec :: dysoil_pfs           ! seq vec
  Vec :: dzsoil_clmp          ! mpi vec
  Vec :: dzsoil_pfs           ! seq vec
  ! a NOTE here: Given a 3D-cell's 'area_gtop_face' and 'zsoi' known, it's possible to calculate its volume (may be useful ?)

   ! cell IDs (in 3D) (for tesing meshes)
  Vec :: cellid_clmp          ! mpi vec
  Vec :: cellid_pfs           ! seq vec
  Vec :: cellid_pfp           ! mpi vec
  Vec :: cellid_clms          ! seq vec
   ! top layer cell IDs (in 2D) (for tesing meshes)
  Vec :: cellid_2dtop_clmp ! mpi vec
  Vec :: cellid_2dtop_pfs  ! seq vec
  Vec :: cellid_2dtop_pfp  ! mpi vec
  Vec :: cellid_2dtop_clms ! seq vec

  ! -----TH vecs from CLM to PF --------------------
  ! TH properties
  PetscReal :: pressure_reference

  ! CLM's hydraulic properties
  Vec :: hksat_x_clmp
  Vec :: hksat_y_clmp
  Vec :: hksat_z_clmp
  Vec :: watsat_clmp
  Vec :: watfc_clmp
  Vec :: bulkdensity_dry_clmp
  Vec :: effporosity_clmp

  Vec :: hksat_x_pfs
  Vec :: hksat_y_pfs
  Vec :: hksat_z_pfs
  Vec :: watsat_pfs
  Vec :: watfc_pfs
  Vec :: bulkdensity_dry_pfs
  Vec :: effporosity_pfs

  Vec :: sucsat_clmp   ! clapp-Horburger's function parameters - needed in BGC somehow
  Vec :: bsw_clmp
  Vec :: sucsat_pfs
  Vec :: bsw_pfs

  ! CLM's thermal properties
  Vec :: tkwet_clmp     ! unit: W/m/K
  Vec :: tkdry_clmp
  Vec :: tkfrz_clmp
  Vec :: hcvsol_clmp    ! unit: J/m^3-K

  Vec :: tkwet_pfs
  Vec :: tkdry_pfs
  Vec :: tkfrz_pfs
  Vec :: hcvsol_pfs

  ! -----TH vecs from PF to CLM --------------------
  ! PF's TH properties (useful to do some calculation in the interface)
  Vec :: sr_pcwmax_pfp
  Vec :: pcwmax_pfp
  Vec :: effporosity_pfp
  Vec :: sr_pcwmax_clms
  Vec :: pcwmax_clms
  Vec :: effporosity_clms

  !---------------------------------------------------------------

  end type clm_pflotran_idata_type

  type(clm_pflotran_idata_type) , public, target , save :: clm_pf_idata
  
  public :: CLMPFLOTRANIDataInit, &
            CLMPFLOTRANIDataCreateVec, &
            CLMPFLOTRANIDataDestroy
  
contains

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataInit()
  ! 
  ! This routine initialized the data transfer type.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  ! Revised by Fengming Yuan, CCSI-ORNL @May-2015
  
    implicit none

    nullify(clm_pf_idata%dxclm_global)
    nullify(clm_pf_idata%dyclm_global)
    nullify(clm_pf_idata%dzclm_global)

    clm_pf_idata%npx = 1    !default 'np' for PF mesh decompose is 1x1x1
    clm_pf_idata%npy = 1
    clm_pf_idata%npz = 1
    nullify(clm_pf_idata%clm_lx)
    nullify(clm_pf_idata%clm_ly)
    nullify(clm_pf_idata%clm_lz)

    clm_pf_idata%nzclm_mapped = 0
    clm_pf_idata%nxclm_mapped = 0
    clm_pf_idata%nyclm_mapped = 0

    clm_pf_idata%x0clm_global = 0
    clm_pf_idata%y0clm_global = 0
    clm_pf_idata%z0clm_global = 0

    clm_pf_idata%nlclm_sub = 0
    clm_pf_idata%ngclm_sub = 0
    clm_pf_idata%nlpf_sub  = 0
    clm_pf_idata%ngpf_sub  = 0

    clm_pf_idata%nlclm_2dtop = 0
    clm_pf_idata%ngclm_2dtop = 0
    clm_pf_idata%nlpf_2dtop  = 0
    clm_pf_idata%ngpf_2dtop  = 0

    clm_pf_idata%nlclm_2dbot = 0
    clm_pf_idata%ngclm_2dbot = 0
    clm_pf_idata%nlpf_2dbot  = 0
    clm_pf_idata%ngpf_2dbot  = 0

    clm_pf_idata%nlclm_srf = 0
    clm_pf_idata%ngclm_srf = 0
    clm_pf_idata%nlpf_srf  = 0
    clm_pf_idata%ngpf_srf  = 0

    !
    clm_pf_idata%zsoil_clmp      = PETSC_NULL_VEC
    clm_pf_idata%zsoil_pfs       = PETSC_NULL_VEC
    clm_pf_idata%zsoil_pfp       = PETSC_NULL_VEC
    clm_pf_idata%zsoil_clms      = PETSC_NULL_VEC
    clm_pf_idata%dxsoil_clmp     = PETSC_NULL_VEC
    clm_pf_idata%dxsoil_pfs      = PETSC_NULL_VEC
    clm_pf_idata%dysoil_clmp     = PETSC_NULL_VEC
    clm_pf_idata%dysoil_pfs      = PETSC_NULL_VEC
    clm_pf_idata%dzsoil_clmp     = PETSC_NULL_VEC
    clm_pf_idata%dzsoil_pfs      = PETSC_NULL_VEC
    clm_pf_idata%xsoil_clmp      = PETSC_NULL_VEC
    clm_pf_idata%xsoil_pfs       = PETSC_NULL_VEC
    clm_pf_idata%ysoil_clmp      = PETSC_NULL_VEC
    clm_pf_idata%ysoil_pfs       = PETSC_NULL_VEC
    clm_pf_idata%zisoil_clmp     = PETSC_NULL_VEC
    clm_pf_idata%zisoil_pfs      = PETSC_NULL_VEC

    clm_pf_idata%area_subsurf_clmp     = PETSC_NULL_VEC
    clm_pf_idata%area_subsurf_pfs      = PETSC_NULL_VEC
    clm_pf_idata%area_subsurf_pfp      = PETSC_NULL_VEC
    clm_pf_idata%area_subsurf_clms     = PETSC_NULL_VEC

    clm_pf_idata%area_top_face_clmp = PETSC_NULL_VEC
    clm_pf_idata%area_top_face_pfs  = PETSC_NULL_VEC
    clm_pf_idata%area_top_face_pfp  = PETSC_NULL_VEC
    clm_pf_idata%area_top_face_clms = PETSC_NULL_VEC

    clm_pf_idata%cellid_clmp     = PETSC_NULL_VEC
    clm_pf_idata%cellid_pfs      = PETSC_NULL_VEC
    clm_pf_idata%cellid_pfp      = PETSC_NULL_VEC
    clm_pf_idata%cellid_clms     = PETSC_NULL_VEC

    clm_pf_idata%cellid_2dtop_clmp     = PETSC_NULL_VEC
    clm_pf_idata%cellid_2dtop_pfs      = PETSC_NULL_VEC
    clm_pf_idata%cellid_2dtop_pfp      = PETSC_NULL_VEC
    clm_pf_idata%cellid_2dtop_clms     = PETSC_NULL_VEC

    !-------------
    clm_pf_idata%pressure_reference = 1.01325d5

    !--------------------------------------------------------------------
    clm_pf_idata%hksat_x_clmp = PETSC_NULL_VEC
    clm_pf_idata%hksat_y_clmp = PETSC_NULL_VEC
    clm_pf_idata%hksat_z_clmp = PETSC_NULL_VEC
    clm_pf_idata%watsat_clmp  = PETSC_NULL_VEC
    clm_pf_idata%watfc_clmp   = PETSC_NULL_VEC
    clm_pf_idata%bulkdensity_dry_clmp = PETSC_NULL_VEC
    clm_pf_idata%effporosity_clmp     = PETSC_NULL_VEC

    clm_pf_idata%tkwet_clmp  = PETSC_NULL_VEC
    clm_pf_idata%tkdry_clmp  = PETSC_NULL_VEC
    clm_pf_idata%tkfrz_clmp  = PETSC_NULL_VEC
    clm_pf_idata%hcvsol_clmp = PETSC_NULL_VEC

    clm_pf_idata%hksat_x_pfs = PETSC_NULL_VEC
    clm_pf_idata%hksat_y_pfs = PETSC_NULL_VEC
    clm_pf_idata%hksat_z_pfs = PETSC_NULL_VEC
    clm_pf_idata%watsat_pfs  = PETSC_NULL_VEC
    clm_pf_idata%watfc_pfs   = PETSC_NULL_VEC
    clm_pf_idata%bulkdensity_dry_pfs = PETSC_NULL_VEC
    clm_pf_idata%effporosity_pfs     = PETSC_NULL_VEC

    clm_pf_idata%tkwet_pfs  = PETSC_NULL_VEC
    clm_pf_idata%tkdry_pfs  = PETSC_NULL_VEC
    clm_pf_idata%tkfrz_pfs  = PETSC_NULL_VEC
    clm_pf_idata%hcvsol_pfs = PETSC_NULL_VEC

    clm_pf_idata%sucsat_clmp = PETSC_NULL_VEC
    clm_pf_idata%bsw_clmp    = PETSC_NULL_VEC
    clm_pf_idata%sucsat_pfs  = PETSC_NULL_VEC
    clm_pf_idata%bsw_pfs     = PETSC_NULL_VEC
   
   !--------------------------------------------------------------------
    clm_pf_idata%sr_pcwmax_pfp   = PETSC_NULL_VEC
    clm_pf_idata%pcwmax_pfp      = PETSC_NULL_VEC
    clm_pf_idata%effporosity_pfp = PETSC_NULL_VEC
    clm_pf_idata%sr_pcwmax_clms  = PETSC_NULL_VEC
    clm_pf_idata%pcwmax_clms     = PETSC_NULL_VEC
    clm_pf_idata%effporosity_clms= PETSC_NULL_VEC

  end subroutine CLMPFLOTRANIDataInit

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataCreateVec(mycomm)
  ! 
  ! This routine creates PETSc vectors required for data transfer between
  ! CLM and PFLOTRAN.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! Revised by Fengming Yuan, CCSI-ORNL @May-2015
  
    implicit none
    
    PetscErrorCode :: ierr
    PetscMPIInt    :: mycomm, rank

    call MPI_Comm_rank(mycomm,rank, ierr)

    ! The following block of data definition is for THC coupled clm-pflotran (Currently ONLY subsurface or soil)
    !
    !NOTES (fmy): From mpi vecs To seq. vecs for passing data IS in one-way only at this momment
    !             (1) First, here will create 4 sets of 3D/2D vecs.
    !             (2) then, below will copy these vecs to create vecs for other variables.

    ! -------- FOR CLM (mpi) ==> PFLOTRAN (seq)
    ! CLM(mpi)
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%zsoil_clmp,ierr)             ! 3D Subsurface PFLOTRAN
    call VecSet(clm_pf_idata%zsoil_clmp,0.d0,ierr)
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_2dtop,PETSC_DECIDE,clm_pf_idata%area_subsurf_clmp,ierr)     ! 2D top-cells of 3D Subsurface PFLOTRAN
    call VecSet(clm_pf_idata%area_subsurf_clmp,0.d0,ierr)
    ! PFLOTRAN(seq)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%zsoil_pfs,ierr)                   ! 3D Subsurface CLM
    call VecSet(clm_pf_idata%zsoil_pfs,0.d0,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_2dtop,clm_pf_idata%area_subsurf_pfs,ierr)           ! 2D top-cells of 3D Subsurface CLM
    call VecSet(clm_pf_idata%area_subsurf_pfs,0.d0,ierr)

    ! -------- FOR PFLOTRAN (mpi) ==> CLM (seq)
    ! PFLOTRAN(mpi)
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_sub,PETSC_DECIDE,clm_pf_idata%zsoil_pfp,ierr)               ! 3D Subsurface PFLOTRAN
    call VecSet(clm_pf_idata%zsoil_pfp,0.d0,ierr)
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_2dtop,PETSC_DECIDE,clm_pf_idata%area_subsurf_pfp,ierr)     ! 2D top-cells of 3D Subsurface PFLOTRAN
    call VecSet(clm_pf_idata%area_subsurf_pfp,0.d0,ierr)
    ! CLM(seq)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_sub,clm_pf_idata%zsoil_clms,ierr)                 ! 3D Subsurface CLM
    call VecSet(clm_pf_idata%zsoil_clms,0.d0,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_2dtop,clm_pf_idata%area_subsurf_clms,ierr)       ! 2D top-cells of 3D Subsurface CLM
    call VecSet(clm_pf_idata%area_subsurf_clms,0.d0,ierr)



    !
    ! ---------- For data transfer from CLM to PFLOTRAN -----------------------------------
    !

    ! (by copying) Create MPI Vectors for CLM ---------------------------------

    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%xsoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%ysoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%zisoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%dxsoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%dysoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%dzsoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%area_top_face_clmp,ierr)

    ! soil cell ids (3D) / surface cell ids (2D)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%cellid_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%cellid_2dtop_clmp,ierr)

    ! soil physical properties (3D)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%hksat_x_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%hksat_y_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%hksat_z_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%watsat_clmp,ierr)       ! total vwc at saturation (total 'porosity')
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%watfc_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%bulkdensity_dry_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%effporosity_clmp,ierr)     ! this may/may not same as 'bd'/'watsat' above

    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%tkwet_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%tkdry_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%tkfrz_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%hcvsol_clmp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%sucsat_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%bsw_clmp,ierr)

    ! (by copying) Create Seq. Vectors for PFLOTRAN  ----------------------

    ! 3-D
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%xsoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%ysoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%zisoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%dxsoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%dysoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%dzsoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%area_top_face_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%cellid_pfs,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%cellid_2dtop_pfs,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%hksat_x_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%hksat_y_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%hksat_z_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%watsat_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%watfc_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%bulkdensity_dry_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%effporosity_pfs,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%tkwet_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%tkdry_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%tkfrz_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%hcvsol_pfs,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%sucsat_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%bsw_pfs,ierr)

    !
    ! --------- For data transfer from PFLOTRAN to CLM  -------------------------------
    !

    ! (by copying) Create MPI Vectors for PFLOTRAN ------------

    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%area_top_face_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%cellid_pfp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%cellid_2dtop_pfp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%sr_pcwmax_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%pcwmax_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%effporosity_pfp,ierr)

    ! (by copying) create Seq. Vectors for CLM ---------

    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%area_top_face_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%cellid_clms,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%cellid_2dtop_clms,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%sr_pcwmax_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%pcwmax_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%effporosity_clms,ierr)

    !--------------------------------------------------------------------------------------------------------------------------

  end subroutine CLMPFLOTRANIDataCreateVec

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataDestroy()
  ! 
  ! This routine destroys PETSc vectors that were created for data transfer.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  ! Revised by Fengming Yuan, CCSI-ORNL @May-2015
  
    implicit none
    
    PetscErrorCode :: ierr

    if (associated(clm_pf_idata%dxclm_global)) &
    deallocate(clm_pf_idata%dxclm_global)
    if (associated(clm_pf_idata%dyclm_global)) &
    deallocate(clm_pf_idata%dyclm_global)
    if (associated(clm_pf_idata%dzclm_global)) &
    deallocate(clm_pf_idata%dzclm_global)

    !----------------------------------------------------------------------------------

    if(clm_pf_idata%zsoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zsoil_clmp,ierr)
    if(clm_pf_idata%zsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zsoil_pfs,ierr)
    if(clm_pf_idata%zsoil_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zsoil_pfp,ierr)
    if(clm_pf_idata%zsoil_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zsoil_clms,ierr)

    if(clm_pf_idata%xsoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%xsoil_clmp,ierr)
    if(clm_pf_idata%xsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%xsoil_pfs,ierr)
    if(clm_pf_idata%ysoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%ysoil_clmp,ierr)
    if(clm_pf_idata%ysoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%ysoil_pfs,ierr)
    if(clm_pf_idata%zisoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zisoil_clmp,ierr)
    if(clm_pf_idata%zisoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zisoil_pfs,ierr)

    if(clm_pf_idata%dxsoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dxsoil_clmp,ierr)
    if(clm_pf_idata%dxsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dxsoil_pfs,ierr)
    if(clm_pf_idata%dysoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dysoil_clmp,ierr)
    if(clm_pf_idata%dysoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dysoil_pfs,ierr)
    if(clm_pf_idata%dzsoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dzsoil_clmp,ierr)
    if(clm_pf_idata%dzsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dzsoil_pfs,ierr)

    if(clm_pf_idata%area_subsurf_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_subsurf_clmp,ierr)
    if(clm_pf_idata%area_subsurf_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_subsurf_pfs,ierr)
    if(clm_pf_idata%area_subsurf_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_subsurf_pfp,ierr)
    if(clm_pf_idata%area_subsurf_clms  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_subsurf_clms,ierr)

    if(clm_pf_idata%area_top_face_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_top_face_clmp,ierr)
    if(clm_pf_idata%area_top_face_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_top_face_pfs,ierr)
    if(clm_pf_idata%area_top_face_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_top_face_pfp,ierr)
    if(clm_pf_idata%area_top_face_clms  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_top_face_clms,ierr)

    !----
    if(clm_pf_idata%cellid_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_clmp,ierr)
    if(clm_pf_idata%cellid_2dtop_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_2dtop_clmp,ierr)
    if(clm_pf_idata%cellid_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_pfs,ierr)
    if(clm_pf_idata%cellid_2dtop_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_2dtop_pfs,ierr)
    if(clm_pf_idata%hksat_x_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_x_clmp,ierr)
    if(clm_pf_idata%hksat_y_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_y_clmp,ierr)
    if(clm_pf_idata%hksat_z_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_z_clmp,ierr)
    if(clm_pf_idata%sucsat_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%sucsat_clmp,ierr)
    if(clm_pf_idata%watsat_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%watsat_clmp,ierr)
    if(clm_pf_idata%bsw_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%bsw_clmp,ierr)
    if(clm_pf_idata%watfc_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%watfc_clmp,ierr)
    if(clm_pf_idata%bulkdensity_dry_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%bulkdensity_dry_clmp,ierr)

    if(clm_pf_idata%cellid_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_pfp,ierr)
    if(clm_pf_idata%cellid_2dtop_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_2dtop_pfp,ierr)
    if(clm_pf_idata%cellid_clms  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_clms,ierr)
    if(clm_pf_idata%cellid_2dtop_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_2dtop_clms,ierr)
    if(clm_pf_idata%tkwet_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkwet_clmp,ierr)
    if(clm_pf_idata%tkdry_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkdry_clmp,ierr)
    if(clm_pf_idata%tkfrz_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkfrz_clmp,ierr)
    if(clm_pf_idata%hcvsol_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hcvsol_clmp,ierr)

    if(clm_pf_idata%hksat_x_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_x_pfs,ierr)
    if(clm_pf_idata%hksat_y_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_y_pfs,ierr)
    if(clm_pf_idata%hksat_z_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_z_pfs,ierr)
    if(clm_pf_idata%sucsat_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%sucsat_pfs,ierr)
    if(clm_pf_idata%watsat_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%watsat_pfs,ierr)
    if(clm_pf_idata%bsw_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%bsw_pfs,ierr)
    if(clm_pf_idata%watfc_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%watfc_pfs,ierr)
    if(clm_pf_idata%bulkdensity_dry_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%bulkdensity_dry_pfs,ierr)
    if(clm_pf_idata%effporosity_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%effporosity_clmp,ierr)
    if(clm_pf_idata%effporosity_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%effporosity_pfs,ierr)

    !----
    if(clm_pf_idata%tkwet_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkwet_pfs,ierr)
    if(clm_pf_idata%tkdry_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkdry_pfs,ierr)
    if(clm_pf_idata%tkfrz_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkfrz_pfs,ierr)
    if(clm_pf_idata%hcvsol_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hcvsol_pfs,ierr)

    ! -----
    if(clm_pf_idata%sr_pcwmax_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%sr_pcwmax_pfp,ierr)
    if(clm_pf_idata%pcwmax_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%pcwmax_pfp,ierr)
    if(clm_pf_idata%sr_pcwmax_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%sr_pcwmax_clms,ierr)
    if(clm_pf_idata%pcwmax_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%pcwmax_clms,ierr)
    if(clm_pf_idata%effporosity_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%effporosity_pfp,ierr)
    if(clm_pf_idata%effporosity_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%effporosity_clms,ierr)

    !
    ! -----------------------------------------------------------------------------------------------------------

  end subroutine CLMPFLOTRANIDataDestroy

end module clm_pflotran_interface_data
