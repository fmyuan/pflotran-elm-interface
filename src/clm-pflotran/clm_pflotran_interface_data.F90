module clm_pflotran_interface_data

  implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  private

  type, public :: clm_pflotran_idata_type

  ! Time invariant data:
  
  ! (i) Soil properties -
  ! Local for CLM  - mpi vectors
  Vec :: hksat_x_clm
  Vec :: hksat_y_clm
  Vec :: hksat_z_clm
  Vec :: sucsat_clm
  Vec :: watsat_clm
  Vec :: bsw_clm
  Vec :: press_clm

  ! Local for PFLOTRAN - seq. vec
  Vec :: hksat_x_pf
  Vec :: hksat_y_pf
  Vec :: hksat_z_pf
  Vec :: sucsat_pf
  Vec :: watsat_pf
  Vec :: bsw_pf
  Vec :: press_pf

  ! (ii) Mesh property

  ! Area of top face
  Vec :: area_top_face_clm ! seq vec
  Vec :: area_top_face_pf  ! mpi vec

  ! Time variant data
  
  ! (i) Sink/Source of water for PFLOTRAN's 3D subsurface domain
  Vec :: qflx_clm   ! mpi vec
  Vec :: qflx_pf    ! seq vec
  
  ! (ii) Source of water for PFLOTRAN's 2D surface domain
  Vec :: rain_clm   ! mpi vec
  Vec :: rain_pf    ! seq vec
  
  ! (iii) Ground heat flux BC for PFLOTRAN's subsurface domain
  !       This BC is applied on top surface of the subsurface domain
  Vec :: gflux_clm  ! mpi vec
  Vec :: gflux_pf   ! seq vec

  ! (iv) Saturation
  Vec :: sat_clm    ! seq vec
  Vec :: sat_pf     ! mpi vec

  ! (v) Subsurface temperature
  Vec :: temp_clm   ! seq vec
  Vec :: temp_pf    ! mpi vec

  ! (vi) Ice saturation
  Vec :: sat_ice_clm ! seq vec
  Vec :: sat_ice_pf  ! mpi vec

  ! (vii) ground/soil C/N pools (G.-P. Tang)
  Vec :: decomp_cpools_vr_lit1_clm     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_clm     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_clm     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  !Vec :: decomp_cpools_vr_cwd_clm      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_clm     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_clm     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_clm     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_clm     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_clm     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_clm     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_clm     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  !Vec :: decomp_npools_vr_cwd_clm     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: sminn_vr_clm                  ! (gN/m3) vertically-resolved soil mineral N
  !Vec :: smin_no3_vr_clm               ! (gN/m3) vertically-resolved soil mineral NO3
  !Vec :: smin_nh4_vr_clm               ! (gN/m3) vertically-resolved soil mineral NH4

  Vec :: decomp_cpools_vr_lit1_pf      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_pf      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_pf      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  !Vec :: decomp_cpools_vr_cwd_pf       ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_pf      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_pf      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_pf      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_pf      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_pf      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_pf      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_pf      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  !Vec :: decomp_npools_vr_cwd_pf       ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: sminn_vr_pf                   ! (gN/m3) vertically-resolved soil mineral N
  !Vec :: smin_no3_vr_pf                ! (gN/m3) vertically-resolved soil mineral NO3
  !Vec :: smin_nh4_vr_pf                ! (gN/m3) vertically-resolved soil mineral NH4

  ! (viii) ground/soil C/N rates (G.-P. Tang)
    Vec :: rate_lit1c_clm
    Vec :: rate_lit2c_clm
    Vec :: rate_lit3c_clm
    Vec :: rate_lit1n_clm
    Vec :: rate_lit2n_clm
    Vec :: rate_lit3n_clm
    !Vec :: rate_cwdc_clm
    !Vec :: rate_cwdn_clm
    Vec :: rate_minn_clm
    Vec :: rate_plantnuptake_clm

    Vec :: rate_lit1c_pf
    Vec :: rate_lit2c_pf
    Vec :: rate_lit3c_pf
    Vec :: rate_lit1n_pf
    Vec :: rate_lit2n_pf
    Vec :: rate_lit3n_pf
    !Vec :: rate_cwdc_pf
    !Vec :: rate_cwdn_pf
    Vec :: rate_minn_pf
    Vec :: rate_plantnuptake_pf

    ! 'hrc' is accumulative in 'PFLOTRAN', so needs previous time-step to calculate 'hr' fluxes for CLM-CN
    Vec :: hrc_vr_clm_prv                ! (gN/m3) vertically-resolved soil heterotrophic respiration C at previous time-step
    Vec :: hrc_vr_clm                    ! (gN/m3) vertically-resolved soil heterotrophic respiration C
    Vec :: hrc_vr_pf                     ! (gN/m3) vertically-resolved soil heterotrophic respiration C

  ! Number of cells for the 3D subsurface domain
  PetscInt :: nlclm_3d  ! num of local clm cells
  PetscInt :: ngclm_3d  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_3d   ! num of local pflotran cells
  PetscInt :: ngpf_3d   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the surface of the 3D subsurface domain
  PetscInt :: nlclm_surf_3d  ! num of local clm cells
  PetscInt :: ngclm_surf_3d  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_surf_3d   ! num of local pflotran cells
  PetscInt :: ngpf_surf_3d   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the 2D surface domain
  PetscInt :: nlclm_2d  ! num of local clm cells
  PetscInt :: ngclm_2d  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2d   ! num of local pflotran cells
  PetscInt :: ngpf_2d   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  end type clm_pflotran_idata_type

  type(clm_pflotran_idata_type) , public, target , save :: clm_pf_idata
  
  public :: CLMPFLOTRANIDataInit, &
            CLMPFLOTRANIDataCreateVec, &
            CLMPFLOTRANIDataDestroy
  
contains

  ! ************************************************************************** !
  !> This routine initialized the data transfer type.
  !!
  !> @author
  !! Gautam Bisht, LBNL
  !!
  !! date: 4/10/2013
  ! ************************************************************************** !
  subroutine CLMPFLOTRANIDataInit()
  
    implicit none

    clm_pf_idata%nlclm_3d = 0
    clm_pf_idata%ngclm_3d = 0
    clm_pf_idata%nlpf_3d = 0
    clm_pf_idata%ngpf_3d = 0

    clm_pf_idata%nlclm_surf_3d = 0
    clm_pf_idata%ngclm_surf_3d = 0
    clm_pf_idata%nlpf_surf_3d = 0
    clm_pf_idata%ngpf_surf_3d = 0

    clm_pf_idata%nlclm_2d = 0
    clm_pf_idata%ngclm_2d = 0
    clm_pf_idata%nlpf_2d = 0
    clm_pf_idata%ngpf_2d = 0

    clm_pf_idata%hksat_x_clm = 0
    clm_pf_idata%hksat_y_clm = 0
    clm_pf_idata%hksat_z_clm = 0
    clm_pf_idata%sucsat_clm = 0
    clm_pf_idata%watsat_clm = 0
    clm_pf_idata%bsw_clm = 0
    clm_pf_idata%press_clm = 0

    clm_pf_idata%hksat_x_pf = 0
    clm_pf_idata%hksat_y_pf = 0
    clm_pf_idata%hksat_z_pf = 0
    clm_pf_idata%sucsat_pf = 0
    clm_pf_idata%watsat_pf = 0
    clm_pf_idata%bsw_pf = 0
    clm_pf_idata%press_pf = 0

    clm_pf_idata%qflx_clm = 0
    clm_pf_idata%qflx_pf = 0
    
    clm_pf_idata%rain_clm = 0
    clm_pf_idata%rain_pf = 0
    
    clm_pf_idata%gflux_clm = 0
    clm_pf_idata%gflux_pf = 0

    clm_pf_idata%sat_clm = 0
    clm_pf_idata%sat_pf = 0

    clm_pf_idata%temp_clm = 0
    clm_pf_idata%temp_pf = 0

    clm_pf_idata%sat_ice_clm = 0
    clm_pf_idata%sat_ice_pf = 0

   ! (vii) soil C/N pools
    clm_pf_idata%decomp_cpools_vr_lit1_clm = 0
    clm_pf_idata%decomp_cpools_vr_lit2_clm = 0
    clm_pf_idata%decomp_cpools_vr_lit3_clm = 0
    !clm_pf_idata%decomp_cpools_vr_cwd_clm  = 0
    clm_pf_idata%decomp_cpools_vr_som1_clm = 0
    clm_pf_idata%decomp_cpools_vr_som2_clm = 0
    clm_pf_idata%decomp_cpools_vr_som3_clm = 0
    clm_pf_idata%decomp_cpools_vr_som4_clm = 0
    clm_pf_idata%decomp_npools_vr_lit1_clm = 0
    clm_pf_idata%decomp_npools_vr_lit2_clm = 0
    clm_pf_idata%decomp_npools_vr_lit3_clm = 0
    !clm_pf_idata%decomp_npools_vr_cwd_clm  = 0
    clm_pf_idata%sminn_vr_clm         = 0
    !clm_pf_idata%smin_no3_vr_clm      = 0
    !clm_pf_idata%smin_nh4_vr_clm      = 0

    clm_pf_idata%hrc_vr_clm_prv       = 0
    clm_pf_idata%hrc_vr_clm           = 0

    clm_pf_idata%decomp_cpools_vr_lit1_pf = 0
    clm_pf_idata%decomp_cpools_vr_lit2_pf = 0
    clm_pf_idata%decomp_cpools_vr_lit3_pf = 0
    !clm_pf_idata%decomp_cpools_vr_cwd_pf  = 0
    clm_pf_idata%decomp_cpools_vr_som1_pf = 0
    clm_pf_idata%decomp_cpools_vr_som2_pf = 0
    clm_pf_idata%decomp_cpools_vr_som3_pf = 0
    clm_pf_idata%decomp_cpools_vr_som4_pf = 0
    clm_pf_idata%decomp_npools_vr_lit1_pf = 0
    clm_pf_idata%decomp_npools_vr_lit2_pf = 0
    clm_pf_idata%decomp_npools_vr_lit3_pf = 0
    !clm_pf_idata%decomp_npools_vr_cwd_pf  = 0
    clm_pf_idata%sminn_vr_pf          = 0
    !clm_pf_idata%smin_no3_vr_pf       = 0
    !clm_pf_idata%smin_nh4_vr_pf       = 0

    clm_pf_idata%hrc_vr_pf            = 0

   ! (viii) ground/soil C/N rates (G.-P. Tang)
    clm_pf_idata%rate_lit1c_clm            = 0
    clm_pf_idata%rate_lit2c_clm            = 0
    clm_pf_idata%rate_lit3c_clm            = 0
    clm_pf_idata%rate_lit1n_clm            = 0
    clm_pf_idata%rate_lit2n_clm            = 0
    clm_pf_idata%rate_lit3n_clm            = 0
    !clm_pf_idata%rate_cwdc_clm            = 0
    !clm_pf_idata%rate_cwdn_clm            = 0
    clm_pf_idata%rate_minn_clm             = 0
    clm_pf_idata%rate_plantnuptake_clm     = 0

    clm_pf_idata%rate_lit1c_pf            = 0
    clm_pf_idata%rate_lit2c_pf            = 0
    clm_pf_idata%rate_lit3c_pf            = 0
    clm_pf_idata%rate_lit1n_pf            = 0
    clm_pf_idata%rate_lit2n_pf            = 0
    clm_pf_idata%rate_lit3n_pf            = 0
    !clm_pf_idata%rate_cwdc_pf            = 0
    !clm_pf_idata%rate_cwdn_pf            = 0
    clm_pf_idata%rate_minn_pf             = 0
    clm_pf_idata%rate_plantnuptake_pf     = 0

  end subroutine CLMPFLOTRANIDataInit

  ! ************************************************************************** !
  !> This routine creates PETSc vectors required for data transfer between
  !! CLM and PFLOTRAN.
  !!
  !> @author
  !! Gautam Bisht, ORNL
  !!
  !! date: 2011
  ! ************************************************************************** !
  subroutine CLMPFLOTRANIDataCreateVec(mycomm)
  
    implicit none
    
    PetscErrorCode :: ierr
    PetscMPIInt    :: mycomm, rank
    PetscReal      :: zero = 0.0d0
    Vec :: vec_test

 !-------------------------------------------------------------------------------------------------

    call MPI_Comm_rank(mycomm,rank, ierr)

 ! F.-M. YUAN: all vectors need initialized (reorganized the code)
    !
    ! Create MPI Vectors for CLM ---------------------------------------------------------------
    ! (i) soil TH variables - parameters
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_3d,PETSC_DECIDE,clm_pf_idata%hksat_x_clm,ierr)
    call VecSet(clm_pf_idata%hksat_x_clm,0.d0,ierr)

    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%hksat_y_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%hksat_z_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%sucsat_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%watsat_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%bsw_clm,ierr)

    ! (ii) soil TH variables - states
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%press_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%sat_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%temp_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%sat_ice_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%area_top_face_clm,ierr)

    ! (iii) soil TH variables - fluxes
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%qflx_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%gflux_clm,ierr)

    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_2d,PETSC_DECIDE,clm_pf_idata%rain_clm,ierr)
    call VecSet(clm_pf_idata%rain_clm,0.d0,ierr)

    ! (iv) bgc variables - states
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_cpools_vr_lit1_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_cpools_vr_lit2_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_cpools_vr_lit3_clm,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_cpools_vr_cwd_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_cpools_vr_som1_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_cpools_vr_som2_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_cpools_vr_som3_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_cpools_vr_som4_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_npools_vr_lit1_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_npools_vr_lit2_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_npools_vr_lit3_clm,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%decomp_npools_vr_cwd_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%sminn_vr_clm,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%smin_no3_vr_clm,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%smin_nh4_vr_clm,ierr)

    !(iv) bgc variable - fluxes
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%rate_lit1c_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%rate_lit2c_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%rate_lit3c_clm,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%rate_cwdc_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%rate_lit1n_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%rate_lit2n_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%rate_lit3n_clm,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%rate_cwdn_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%rate_minn_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%rate_plantnuptake_clm,ierr)

    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%hrc_vr_clm_prv,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%hrc_vr_clm,ierr)

! Create Seq. Vectors for PFLOTRAN -----------------------------------------------------
    ! (i) soil TH variables - parameters
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_3d,clm_pf_idata%hksat_x_pf,ierr)
    call VecSet(clm_pf_idata%hksat_x_pf,0.d0,ierr)

    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%hksat_y_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%hksat_z_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%sucsat_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%watsat_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%bsw_pf,ierr)

    ! (ii) soil TH variables - states
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%press_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%sat_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%sat_ice_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%temp_pf,ierr)

    ! (iii) soil TH variables - fluxes
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%qflx_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%gflux_pf,ierr)

    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_2d,clm_pf_idata%rain_pf,ierr)
    call VecSet(clm_pf_idata%rain_pf,0.d0,ierr)

    ! (iv) soil bgc variables - states
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_cpools_vr_lit1_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_cpools_vr_lit2_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_cpools_vr_lit3_pf,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_cpools_vr_cwd_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_cpools_vr_som1_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_cpools_vr_som2_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_cpools_vr_som3_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_cpools_vr_som4_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_npools_vr_lit1_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_npools_vr_lit2_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_npools_vr_lit3_pf,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%decomp_npools_vr_cwd_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%sminn_vr_pf,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%smin_no3_vr_pf,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%smin_nh4_vr_pf,ierr)

    !(v) bgc variable - fluxes
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%rate_lit1c_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%rate_lit2c_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%rate_lit3c_pf,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%rate_cwdc_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%rate_lit1n_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%rate_lit2n_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%rate_lit3n_pf,ierr)
    !call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%rate_cwdn_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%rate_minn_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%rate_plantnuptake_pf,ierr)

    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%hrc_vr_pf,ierr)

  end subroutine CLMPFLOTRANIDataCreateVec

  ! ************************************************************************** !
  !> This routine destroys PETSc vectors that were created for data transfer.
  !!
  !> @author
  !! Gautam Bisht, LBNL
  !!
  !! date: 4/10/2013
  ! ************************************************************************** !
  subroutine CLMPFLOTRANIDataDestroy()
  
    implicit none
    
    PetscErrorCode :: ierr

    if(clm_pf_idata%hksat_x_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_x_clm,ierr)
    if(clm_pf_idata%hksat_y_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_y_clm,ierr)
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)
    if(clm_pf_idata%sucsat_clm  /= 0) call VecDestroy(clm_pf_idata%sucsat_clm,ierr)
    if(clm_pf_idata%watsat_clm  /= 0) call VecDestroy(clm_pf_idata%watsat_clm,ierr)
    if(clm_pf_idata%bsw_clm  /= 0) call VecDestroy(clm_pf_idata%bsw_clm,ierr)
    if(clm_pf_idata%press_clm  /= 0) call VecDestroy(clm_pf_idata%press_clm,ierr)

    if(clm_pf_idata%hksat_x_pf  /= 0) call VecDestroy(clm_pf_idata%hksat_x_pf,ierr)
    if(clm_pf_idata%hksat_y_pf  /= 0) call VecDestroy(clm_pf_idata%hksat_y_pf,ierr)
    if(clm_pf_idata%hksat_z_pf  /= 0) call VecDestroy(clm_pf_idata%hksat_z_pf,ierr)
    if(clm_pf_idata%sucsat_pf  /= 0) call VecDestroy(clm_pf_idata%sucsat_pf,ierr)
    if(clm_pf_idata%watsat_pf  /= 0) call VecDestroy(clm_pf_idata%watsat_pf,ierr)
    if(clm_pf_idata%bsw_pf  /= 0) call VecDestroy(clm_pf_idata%bsw_pf,ierr)
    if(clm_pf_idata%press_pf  /= 0) call VecDestroy(clm_pf_idata%press_pf,ierr)

    if(clm_pf_idata%qflx_clm  /= 0) call VecDestroy(clm_pf_idata%qflx_clm,ierr)
    if(clm_pf_idata%qflx_pf  /= 0) call VecDestroy(clm_pf_idata%qflx_pf,ierr)
    
    if(clm_pf_idata%rain_clm  /= 0) call VecDestroy(clm_pf_idata%rain_clm,ierr)
    if(clm_pf_idata%rain_pf  /= 0) call VecDestroy(clm_pf_idata%rain_pf,ierr)
    
    if(clm_pf_idata%gflux_clm  /= 0) call VecDestroy(clm_pf_idata%gflux_clm,ierr)
    if(clm_pf_idata%gflux_pf  /= 0) call VecDestroy(clm_pf_idata%gflux_pf,ierr)

    if(clm_pf_idata%sat_clm  /= 0) call VecDestroy(clm_pf_idata%sat_clm,ierr)
    if(clm_pf_idata%sat_pf  /= 0) call VecDestroy(clm_pf_idata%sat_pf,ierr)

    if(clm_pf_idata%temp_clm  /= 0) call VecDestroy(clm_pf_idata%temp_clm,ierr)
    if(clm_pf_idata%temp_pf  /= 0) call VecDestroy(clm_pf_idata%temp_pf,ierr)

    if(clm_pf_idata%sat_ice_clm  /= 0) call VecDestroy(clm_pf_idata%sat_ice_clm,ierr)
    if(clm_pf_idata%sat_ice_pf  /= 0) call VecDestroy(clm_pf_idata%sat_ice_pf,ierr)

    if(clm_pf_idata%area_top_face_clm  /= 0) &
      call VecDestroy(clm_pf_idata%area_top_face_clm,ierr)
    if(clm_pf_idata%area_top_face_pf  /= 0) &
      call VecDestroy(clm_pf_idata%area_top_face_pf,ierr)

    ! soil C/N pools
    if(clm_pf_idata%decomp_cpools_vr_lit1_clm /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_clm,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_clm /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_clm,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_clm /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_clm,ierr)
    !if(clm_pf_idata%decomp_cpools_vr_cwd_clm  /= 0) &
    !   call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_clm,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_clm /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_clm,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_clm /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_clm,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_clm /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_clm,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_clm /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_clm,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_clm /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_clm,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_clm /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_clm,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_clm /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_clm,ierr)
    !if(clm_pf_idata%decomp_npools_vr_cwd_clm  /= 0) &
    !   call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_clm,ierr)
    if(clm_pf_idata%sminn_vr_clm /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_clm,ierr)
    !if(clm_pf_idata%smin_no3_vr_clm /= 0) &
    !   call VecDestroy(clm_pf_idata%smin_no3_vr_clm,ierr)
    !if(clm_pf_idata%smin_nh4_vr_clm /= 0) &
    !   call VecDestroy(clm_pf_idata%smin_nh4_vr_clm,ierr)
!
    if(clm_pf_idata%decomp_cpools_vr_lit1_pf /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_pf,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_pf /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_pf,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_pf /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_pf,ierr)
    !if(clm_pf_idata%decomp_cpools_vr_cwd_pf  /= 0) &
    !   call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_pf,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_pf /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_pf,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_pf /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_pf,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_pf /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_pf,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_pf /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_pf,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_pf /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_pf,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_pf /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_pf,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_pf /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_pf,ierr)
    !if(clm_pf_idata%decomp_npools_vr_cwd_pf  /= 0) &
    !   call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_pf,ierr)
    if(clm_pf_idata%sminn_vr_pf /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_pf,ierr)
    !if(clm_pf_idata%smin_no3_vr_pf /= 0) &
    !   call VecDestroy(clm_pf_idata%smin_no3_vr_pf,ierr)
    !if(clm_pf_idata%smin_nh4_vr_pf /= 0) &
    !  call VecDestroy(clm_pf_idata%smin_nh4_vr_pf,ierr)

    ! soil C/N fluxes
    if(clm_pf_idata%rate_lit1c_clm /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1c_clm,ierr)
    if(clm_pf_idata%rate_lit2c_clm /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2c_clm,ierr)
    if(clm_pf_idata%rate_lit3c_clm /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3c_clm,ierr)
    !if(clm_pf_idata%rate_cwdc_clm /= 0) &
    !   call VecDestroy(clm_pf_idata%rate_cwdc_clm,ierr)
    if(clm_pf_idata%rate_lit1n_clm /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1n_clm,ierr)
    if(clm_pf_idata%rate_lit2n_clm /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2n_clm,ierr)
    if(clm_pf_idata%rate_lit3n_clm /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3n_clm,ierr)
    !if(clm_pf_idata%rate_cwdn_clm /= 0) &
    !   call VecDestroy(clm_pf_idata%rate_cwdn_clm,ierr)
    if(clm_pf_idata%rate_minn_clm /= 0) &
       call VecDestroy(clm_pf_idata%rate_minn_clm,ierr)
    if(clm_pf_idata%rate_plantnuptake_clm /= 0) &
       call VecDestroy(clm_pf_idata%rate_plantnuptake_clm,ierr)

    if(clm_pf_idata%hrc_vr_clm_prv /= 0) &
       call VecDestroy(clm_pf_idata%hrc_vr_clm_prv,ierr)
    if(clm_pf_idata%hrc_vr_clm /= 0) &
       call VecDestroy(clm_pf_idata%hrc_vr_clm,ierr)

    if(clm_pf_idata%rate_lit1c_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1c_pf,ierr)
    if(clm_pf_idata%rate_lit2c_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2c_pf,ierr)
    if(clm_pf_idata%rate_lit3c_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3c_pf,ierr)
    !if(clm_pf_idata%rate_cwdc_pf /= 0) &
    !   call VecDestroy(clm_pf_idata%rate_cwdc_pf,ierr)
    if(clm_pf_idata%rate_lit1n_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1n_pf,ierr)
    if(clm_pf_idata%rate_lit2n_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2n_pf,ierr)
    if(clm_pf_idata%rate_lit3n_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3n_pf,ierr)
    !if(clm_pf_idata%rate_cwdn_pf /= 0) &
    !   call VecDestroy(clm_pf_idata%rate_cwdn_pf,ierr)
    if(clm_pf_idata%rate_minn_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_minn_pf,ierr)
    if(clm_pf_idata%rate_plantnuptake_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_plantnuptake_pf,ierr)

    if(clm_pf_idata%hrc_vr_pf /= 0) &
       call VecDestroy(clm_pf_idata%hrc_vr_pf,ierr)

  end subroutine CLMPFLOTRANIDataDestroy

end module clm_pflotran_interface_data
