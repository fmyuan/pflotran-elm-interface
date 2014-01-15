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
  
  ! (ii) Source of water and temperature of rain for PFLOTRAN's 2D surface domain
  Vec :: rain_clm   ! mpi vec
  Vec :: rain_pf    ! seq vec
  Vec :: rain_temp_clm ! seq vec
  Vec :: rain_temp_pf  ! mpi vec
  
  ! (iii) Ground heat flux BC for PFLOTRAN's subsurface domain
  !       This BC is applied on top surface of the subsurface domain
  Vec :: gflux_subsurf_clm  ! mpi vec
  Vec :: gflux_subsurf_pf   ! seq vec
  !       When running PFLOTRAN surface-subsurface simulation, ground heat flux
  !       is a SS for PFLOTRAN's surface domain.
  !
  !       Note: CLM decomposes the domain across processors in a horizontal.
  !       Thus, nlclm_2dsub = nlclm_srf across all processors. Thus, there is
  !       no need for 'gflux_surf_clm'
  !
  Vec :: gflux_surf_pf   ! seq vec

  ! (iv) Saturation
  Vec :: sat_clm    ! seq vec
  Vec :: sat_pf     ! mpi vec

  ! (v) Subsurface temperature
  Vec :: temp_clm   ! seq vec
  Vec :: temp_pf    ! mpi vec

  ! (vi) Ice saturation
  Vec :: sat_ice_clm ! seq vec
  Vec :: sat_ice_pf  ! mpi vec

  ! (vii) Stand water head
  Vec :: h2osfc_clm ! seq vec
  Vec :: h2osfc_pf  ! mpi vec

  ! Number of cells for the 3D subsurface domain
  PetscInt :: nlclm_sub ! num of local clm cells
  PetscInt :: ngclm_sub ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_sub  ! num of local pflotran cells
  PetscInt :: ngpf_sub   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the surface of the 3D subsurface domain
  PetscInt :: nlclm_2dsub  ! num of local clm cells
  PetscInt :: ngclm_2dsub  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2dsub   ! num of local pflotran cells
  PetscInt :: ngpf_2dsub   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the 2D surface domain
  PetscInt :: nlclm_srf  ! num of local clm cells
  PetscInt :: ngclm_srf  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_srf   ! num of local pflotran cells
  PetscInt :: ngpf_srf   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  !-------------------------------------------------------------------------------------
  ! note: mapping ONLY can do from MPI vecs to Seq. vecs now
  ! -----BGC vecs from CLM to PF --------------------
  ! initial TH state vecs from CLM (mpi) to PF (seq) for bgc
  Vec :: soilpsi_clmp                   ! soil matric potential
  Vec :: soillsat_clmp                  ! soil liq. water saturation
  Vec :: soilisat_clmp                  ! soil ice water saturation
  Vec :: soilt_clmp                     ! soil temperature
  Vec :: soilpsi_pf
  Vec :: soillsat_pfs
  Vec :: soilisat_pfs
  Vec :: soilt_pfs
  ! initial ground/soil C/N pools from CLM (mpi) to PF (seq)
  Vec :: decomp_cpools_vr_lit1_clmp     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_clmp     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_clmp     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_cwd_clmp      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_clmp     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_clmp     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_clmp     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_clmp     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_clmp     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_clmp     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_clmp     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_cwd_clmp      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: sminn_vr_clmp                  ! (gN/m3) vertically-resolved soil mineral N
  Vec :: smin_no3_vr_clmp               ! (gN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_clmp               ! (gN/m3) vertically-resolved soil mineral NH4
  Vec :: decomp_cpools_vr_lit1_pfs      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_pfs      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_pfs      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_cwd_pfs       ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_pfs      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_pfs      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_pfs      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_pfs      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_pfs      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_pfs      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_pfs      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_cwd_pfs       ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: sminn_vr_pfs                   ! (gN/m3) vertically-resolved soil mineral N
  Vec :: smin_no3_vr_pfs                ! (gN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_pfs                ! (gN/m3) vertically-resolved soil mineral NH4
  ! time-varying ground/soil C/N rates from CLM (mpi) to PF (seq) (Units: moleC(N)/m3/s)
  Vec :: rate_lit1c_clmp
  Vec :: rate_lit2c_clmp
  Vec :: rate_lit3c_clmp
  Vec :: rate_lit1n_clmp
  Vec :: rate_lit2n_clmp
  Vec :: rate_lit3n_clmp
  Vec :: rate_cwdc_clmp
  Vec :: rate_cwdn_clmp
  Vec :: rate_minn_clmp
  Vec :: rate_plantnuptake_clmp
  Vec :: rate_nleached_clmp
  Vec :: rate_ndenitri_clmp
  Vec :: rate_lit1c_pfs
  Vec :: rate_lit2c_pfs
  Vec :: rate_lit3c_pfs
  Vec :: rate_lit1n_pfs
  Vec :: rate_lit2n_pfs
  Vec :: rate_lit3n_pfs
  Vec :: rate_cwdc_pfs
  Vec :: rate_cwdn_pfs
  Vec :: rate_minn_pfs
  Vec :: rate_plantnuptake_pf
  Vec :: rate_nleached_pf
  Vec :: rate_ndenitri_pf

  ! -----BGC vecs from PF (mpi, ghosted) to CLM (seq, local) --------------------
  Vec :: decomp_cpools_vr_lit1_pfp
  Vec :: decomp_cpools_vr_lit2_pfp
  Vec :: decomp_cpools_vr_lit3_pfp
  Vec :: decomp_cpools_vr_cwd_pfp
  Vec :: decomp_cpools_vr_som1_pfp
  Vec :: decomp_cpools_vr_som2_pfp
  Vec :: decomp_cpools_vr_som3_pfp
  Vec :: decomp_cpools_vr_som4_pfp
  Vec :: decomp_npools_vr_lit1_pfp
  Vec :: decomp_npools_vr_lit2_pfp
  Vec :: decomp_npools_vr_lit3_pfp
  Vec :: decomp_npools_vr_cwd_pfp
  Vec :: sminn_vr_pfp
  Vec :: smin_no3_vr_pfp
  Vec :: smin_nh4_vr_pfp
  !
  Vec :: decomp_cpools_vr_lit1_clms     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_clms     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_clms     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_cwd_clms      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_clms     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_clms     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_clms     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_clms     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_clms     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_clms     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_clms     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_cwd_clms      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: sminn_vr_clms                  ! (gN/m3) vertically-resolved soil mineral N
  Vec :: smin_no3_vr_clms               ! (gN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_clms               ! (gN/m3) vertically-resolved soil mineral NH4
  !
  Vec :: decomp_cpools_vr_lit1_clms_prv ! C/N states at previous time-step for pool-change calculation
  Vec :: decomp_cpools_vr_lit2_clms_prv
  Vec :: decomp_cpools_vr_lit3_clms_prv
  Vec :: decomp_cpools_vr_cwd_clms_prv
  Vec :: decomp_cpools_vr_som1_clms_prv
  Vec :: decomp_cpools_vr_som2_clms_prv
  Vec :: decomp_cpools_vr_som3_clms_prv
  Vec :: decomp_cpools_vr_som4_clms_prv
  Vec :: decomp_npools_vr_lit1_clms_prv
  Vec :: decomp_npools_vr_lit2_clms_prv
  Vec :: decomp_npools_vr_lit3_clms_prv
  Vec :: decomp_npools_vr_cwd_clms_prv
  Vec :: sminn_vr_clms_prv
  Vec :: smin_no3_vr_clms_prv
  Vec :: smin_nh4_vr_clms_prv

  ! 'hrc' is accumulative in 'PFLOTRAN', so needs previous time-step to calculate 'hr' fluxes for CLM-CN
  Vec :: hrc_vr_pfp                     ! (gC/m3) vertically-resolved soil heterotrophic respiration C
  Vec :: hrc_vr_clms_prv                ! (gC/m3) vertically-resolved soil heterotrophic respiration C at previous time-step
  Vec :: hrc_vr_clms                    ! (gC/m3) vertically-resolved soil heterotrophic respiration C
  ! 'accextrn' is accumulative N extract in 'PFLOTRAN', so needs previous time-step to calculate 'sminn_to_plant' fluxes for CLM-CN
  Vec :: accextrn_vr_pfp                ! (gN/m3) vertically-resolved root extraction N at previous time-step
  Vec :: accextrn_vr_clms_prv           ! (gN/m3) vertically-resolved root extraction N at previous time-step
  Vec :: accextrn_vr_clms               ! (gN/m3) vertically-resolved root extraction N at previous time-step

  !---------------------------------------------------------------

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

    clm_pf_idata%nlclm_sub = 0
    clm_pf_idata%ngclm_sub = 0
    clm_pf_idata%nlpf_sub = 0
    clm_pf_idata%ngpf_sub = 0

    clm_pf_idata%nlclm_2dsub = 0
    clm_pf_idata%ngclm_2dsub = 0
    clm_pf_idata%nlpf_2dsub = 0
    clm_pf_idata%ngpf_2dsub = 0

    clm_pf_idata%nlclm_srf = 0
    clm_pf_idata%ngclm_srf = 0
    clm_pf_idata%nlpf_srf = 0
    clm_pf_idata%ngpf_srf = 0

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
    clm_pf_idata%rain_temp_clm = 0
    clm_pf_idata%rain_temp_pf = 0
    
    clm_pf_idata%gflux_subsurf_clm = 0
    clm_pf_idata%gflux_subsurf_pf = 0
    clm_pf_idata%gflux_surf_pf = 0

    clm_pf_idata%sat_clm = 0
    clm_pf_idata%sat_pf = 0

    clm_pf_idata%temp_clm = 0
    clm_pf_idata%temp_pf = 0

    clm_pf_idata%sat_ice_clm = 0
    clm_pf_idata%sat_ice_pf = 0

    clm_pf_idata%h2osfc_clm = 0
    clm_pf_idata%h2osfc_pf = 0
   
   !--------------------------------------------------------------------
   ! (initial) soil TH and C/N pools
    clm_pf_idata%soilpsi_clmp  = 0
    clm_pf_idata%soillsat_clmp = 0
    clm_pf_idata%soilisat_clmp = 0
    clm_pf_idata%soilt_clmp    = 0
    clm_pf_idata%soilpsi_pf   = 0
    clm_pf_idata%soillsat_pfs = 0
    clm_pf_idata%soilisat_pfs = 0
    clm_pf_idata%soilt_pfs    = 0

    clm_pf_idata%decomp_cpools_vr_lit1_clmp = 0
    clm_pf_idata%decomp_cpools_vr_lit2_clmp = 0
    clm_pf_idata%decomp_cpools_vr_lit3_clmp = 0
    clm_pf_idata%decomp_cpools_vr_cwd_clmp  = 0
    clm_pf_idata%decomp_cpools_vr_som1_clmp = 0
    clm_pf_idata%decomp_cpools_vr_som2_clmp = 0
    clm_pf_idata%decomp_cpools_vr_som3_clmp = 0
    clm_pf_idata%decomp_cpools_vr_som4_clmp = 0
    clm_pf_idata%decomp_npools_vr_lit1_clmp = 0
    clm_pf_idata%decomp_npools_vr_lit2_clmp = 0
    clm_pf_idata%decomp_npools_vr_lit3_clmp = 0
    clm_pf_idata%decomp_npools_vr_cwd_clmp  = 0
    clm_pf_idata%sminn_vr_clmp         = 0
    clm_pf_idata%smin_no3_vr_clmp      = 0
    clm_pf_idata%smin_nh4_vr_clmp      = 0

    clm_pf_idata%decomp_cpools_vr_lit1_pfs = 0
    clm_pf_idata%decomp_cpools_vr_lit2_pfs = 0
    clm_pf_idata%decomp_cpools_vr_lit3_pfs = 0
    clm_pf_idata%decomp_cpools_vr_cwd_pfs  = 0
    clm_pf_idata%decomp_cpools_vr_som1_pfs = 0
    clm_pf_idata%decomp_cpools_vr_som2_pfs = 0
    clm_pf_idata%decomp_cpools_vr_som3_pfs = 0
    clm_pf_idata%decomp_cpools_vr_som4_pfs = 0
    clm_pf_idata%decomp_npools_vr_lit1_pfs = 0
    clm_pf_idata%decomp_npools_vr_lit2_pfs = 0
    clm_pf_idata%decomp_npools_vr_lit3_pfs = 0
    clm_pf_idata%decomp_npools_vr_cwd_pfs  = 0
    clm_pf_idata%sminn_vr_pfs          = 0
    clm_pf_idata%smin_no3_vr_pfs       = 0
    clm_pf_idata%smin_nh4_vr_pfs       = 0

    !ground/soil C/N rates as source/sink
    clm_pf_idata%rate_lit1c_clmp            = 0
    clm_pf_idata%rate_lit2c_clmp            = 0
    clm_pf_idata%rate_lit3c_clmp            = 0
    clm_pf_idata%rate_lit1n_clmp            = 0
    clm_pf_idata%rate_lit2n_clmp            = 0
    clm_pf_idata%rate_lit3n_clmp            = 0
    clm_pf_idata%rate_cwdc_clmp             = 0
    clm_pf_idata%rate_cwdn_clmp             = 0
    clm_pf_idata%rate_minn_clmp             = 0
    clm_pf_idata%rate_plantnuptake_clmp     = 0
    clm_pf_idata%rate_nleached_clmp         = 0
    clm_pf_idata%rate_ndenitri_clmp         = 0

    clm_pf_idata%rate_lit1c_pfs            = 0
    clm_pf_idata%rate_lit2c_pfs            = 0
    clm_pf_idata%rate_lit3c_pfs            = 0
    clm_pf_idata%rate_lit1n_pfs            = 0
    clm_pf_idata%rate_lit2n_pfs            = 0
    clm_pf_idata%rate_lit3n_pfs            = 0
    clm_pf_idata%rate_cwdc_pfs             = 0
    clm_pf_idata%rate_cwdn_pfs             = 0
    clm_pf_idata%rate_minn_pfs             = 0
    clm_pf_idata%rate_plantnuptake_pf      = 0
    clm_pf_idata%rate_nleached_pf          = 0
    clm_pf_idata%rate_ndenitri_pf          = 0

    ! for updating bgc states
    clm_pf_idata%decomp_cpools_vr_lit1_pfp = 0
    clm_pf_idata%decomp_cpools_vr_lit2_pfp = 0
    clm_pf_idata%decomp_cpools_vr_lit3_pfp = 0
    clm_pf_idata%decomp_cpools_vr_cwd_pfp  = 0
    clm_pf_idata%decomp_cpools_vr_som1_pfp = 0
    clm_pf_idata%decomp_cpools_vr_som2_pfp = 0
    clm_pf_idata%decomp_cpools_vr_som3_pfp = 0
    clm_pf_idata%decomp_cpools_vr_som4_pfp = 0
    clm_pf_idata%decomp_npools_vr_lit1_pfp = 0
    clm_pf_idata%decomp_npools_vr_lit2_pfp = 0
    clm_pf_idata%decomp_npools_vr_lit3_pfp = 0
    clm_pf_idata%decomp_npools_vr_cwd_pfp  = 0
    clm_pf_idata%sminn_vr_pfp          = 0
    clm_pf_idata%smin_no3_vr_pfp       = 0
    clm_pf_idata%smin_nh4_vr_pfp       = 0
    clm_pf_idata%decomp_cpools_vr_lit1_clms = 0
    clm_pf_idata%decomp_cpools_vr_lit2_clms = 0
    clm_pf_idata%decomp_cpools_vr_lit3_clms = 0
    clm_pf_idata%decomp_cpools_vr_cwd_clms  = 0
    clm_pf_idata%decomp_cpools_vr_som1_clms = 0
    clm_pf_idata%decomp_cpools_vr_som2_clms = 0
    clm_pf_idata%decomp_cpools_vr_som3_clms = 0
    clm_pf_idata%decomp_cpools_vr_som4_clms = 0
    clm_pf_idata%decomp_npools_vr_lit1_clms = 0
    clm_pf_idata%decomp_npools_vr_lit2_clms = 0
    clm_pf_idata%decomp_npools_vr_lit3_clms = 0
    clm_pf_idata%decomp_npools_vr_cwd_clms  = 0
    clm_pf_idata%sminn_vr_clms         = 0
    clm_pf_idata%smin_no3_vr_clms      = 0
    clm_pf_idata%smin_nh4_vr_clms      = 0
    clm_pf_idata%decomp_cpools_vr_lit1_clms_prv = 0
    clm_pf_idata%decomp_cpools_vr_lit2_clms_prv = 0
    clm_pf_idata%decomp_cpools_vr_lit3_clms_prv = 0
    clm_pf_idata%decomp_cpools_vr_cwd_clms_prv  = 0
    clm_pf_idata%decomp_cpools_vr_som1_clms_prv = 0
    clm_pf_idata%decomp_cpools_vr_som2_clms_prv = 0
    clm_pf_idata%decomp_cpools_vr_som3_clms_prv = 0
    clm_pf_idata%decomp_cpools_vr_som4_clms_prv = 0
    clm_pf_idata%decomp_npools_vr_lit1_clms_prv = 0
    clm_pf_idata%decomp_npools_vr_lit2_clms_prv = 0
    clm_pf_idata%decomp_npools_vr_lit3_clms_prv = 0
    clm_pf_idata%decomp_npools_vr_cwd_clms_prv  = 0
    clm_pf_idata%sminn_vr_clms_prv         = 0
    clm_pf_idata%smin_no3_vr_clms_prv      = 0
    clm_pf_idata%smin_nh4_vr_clms_prv      = 0

    ! for soil hr calculation
    clm_pf_idata%hrc_vr_pfp            = 0
    clm_pf_idata%hrc_vr_clms_prv       = 0
    clm_pf_idata%hrc_vr_clms           = 0

    ! for root N extraction calculation
    clm_pf_idata%accextrn_vr_pfp       = 0
    clm_pf_idata%accextrn_vr_clms_prv  = 0
    clm_pf_idata%accextrn_vr_clms      = 0

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

    call MPI_Comm_rank(mycomm,rank, ierr)

    !
    ! For data transfer from CLM to PFLOTRAN
    !

    ! Create MPI Vectors for CLM
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%hksat_x_clm,ierr)
    call VecSet(clm_pf_idata%hksat_x_clm,0.d0,ierr)

    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%hksat_y_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%hksat_z_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%sucsat_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%watsat_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%bsw_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%press_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%qflx_clm,ierr)

    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_2dsub,PETSC_DECIDE,clm_pf_idata%gflux_subsurf_clm,ierr)
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_srf,PETSC_DECIDE,clm_pf_idata%rain_clm,ierr)
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_srf,PETSC_DECIDE,clm_pf_idata%rain_temp_clm,ierr)
    call VecSet(clm_pf_idata%gflux_subsurf_clm,0.d0,ierr)
    call VecSet(clm_pf_idata%rain_clm,0.d0,ierr)
    call VecSet(clm_pf_idata%rain_temp_clm,0.d0,ierr)

    ! Create Seq. Vectors for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%hksat_x_pf,ierr)
    call VecSet(clm_pf_idata%hksat_x_pf,0.d0,ierr)

    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%hksat_y_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%hksat_z_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%sucsat_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%watsat_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%bsw_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%press_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%qflx_pf,ierr)

    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_2dsub,clm_pf_idata%gflux_subsurf_pf,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_srf,clm_pf_idata%gflux_surf_pf,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_srf,clm_pf_idata%rain_pf,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_srf,clm_pf_idata%rain_temp_pf,ierr)
    call VecSet(clm_pf_idata%gflux_subsurf_pf,0.d0,ierr)
    call VecSet(clm_pf_idata%gflux_surf_pf,0.d0,ierr)
    call VecSet(clm_pf_idata%rain_pf,0.d0,ierr)
    call VecSet(clm_pf_idata%rain_temp_pf,0.d0,ierr)

    !
    ! For data transfer from PFLOTRAN to CLM
    !

    ! Create MPI Vectors for PFLOTRAN
    ! 3D Subsurface PFLOTRAN ---to--- 3D Subsurface CLM
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_sub,PETSC_DECIDE,clm_pf_idata%sat_pf,ierr)
    call VecSet(clm_pf_idata%sat_pf,0.d0,ierr)

    call VecDuplicate(clm_pf_idata%sat_pf,clm_pf_idata%temp_pf,ierr)
    call VecDuplicate(clm_pf_idata%sat_pf,clm_pf_idata%sat_ice_pf,ierr)
    call VecDuplicate(clm_pf_idata%sat_pf,clm_pf_idata%area_top_face_pf,ierr)

    ! 2D Surface PFLOTRAN ---to--- 2D Surface CLM
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_srf,PETSC_DECIDE,clm_pf_idata%h2osfc_pf,ierr)
    call VecSet(clm_pf_idata%h2osfc_pf,0.d0,ierr)

    ! Create Seq. Vectors for CLM
    ! 3D Subsurface PFLOTRAN ---to--- 3D Subsurface CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_sub,clm_pf_idata%sat_clm,ierr)
    call VecSet(clm_pf_idata%sat_clm,0.d0,ierr)

    call VecDuplicate(clm_pf_idata%sat_clm,clm_pf_idata%temp_clm,ierr)
    call VecDuplicate(clm_pf_idata%sat_clm,clm_pf_idata%sat_ice_clm,ierr)
    call VecDuplicate(clm_pf_idata%sat_clm,clm_pf_idata%area_top_face_clm,ierr)

    ! 2D Surface PFLOTRAN ---to--- 2D Surface CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%nlclm_2dsub,clm_pf_idata%h2osfc_clm,ierr)
    call VecSet(clm_pf_idata%h2osfc_clm,0.d0,ierr)

    !------------------------------------------------------
    !NOTES (fmy): From mpi vecs To seq. vecs for passing data IS in one-way only
    ! (i) BGC state variables: 3D subsurface CLM ---to--- 3D subsurface PFLOTRAN (e.g., initialization or abiotic factors)
    ! MPI Vecs for CLM
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%decomp_cpools_vr_lit1_clmp,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_lit1_clmp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_cpools_vr_lit2_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_cpools_vr_lit3_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_cpools_vr_cwd_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_cpools_vr_som1_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_cpools_vr_som2_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_cpools_vr_som3_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_cpools_vr_som4_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_npools_vr_lit1_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_npools_vr_lit2_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_npools_vr_lit3_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_npools_vr_cwd_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%sminn_vr_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%smin_no3_vr_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%smin_nh4_vr_clmp,ierr)

    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%soilpsi_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%soillsat_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%soilisat_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%soilt_clmp,ierr)

    ! Seq. Vecs for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%decomp_cpools_vr_lit1_pfs,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_lit1_pfs,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_cpools_vr_lit2_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_cpools_vr_lit3_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_cpools_vr_cwd_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_cpools_vr_som1_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_cpools_vr_som2_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_cpools_vr_som3_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_cpools_vr_som4_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_npools_vr_lit1_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_npools_vr_lit2_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_npools_vr_lit3_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_npools_vr_cwd_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%sminn_vr_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%smin_no3_vr_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%smin_nh4_vr_pfs,ierr)

    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%soilpsi_pf,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%soillsat_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%soilisat_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%soilt_pfs,ierr)

    ! (ii) BGC interface source/sink (rates): 3D subsurface CLM ---to--- 3D subsurface PFLOTRAN
    ! MPI Vecs for CLM
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%rate_lit1c_clmp,ierr)
    call VecSet(clm_pf_idata%rate_lit1c_clmp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_lit2c_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_lit3c_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_lit1n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_lit2n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_lit3n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_cwdc_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_cwdn_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_minn_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_plantnuptake_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_nleached_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_ndenitri_clmp,ierr)

    ! Seq. Vecs for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%rate_lit1c_pfs,ierr)
    call VecSet(clm_pf_idata%rate_lit1c_pfs,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_lit2c_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_lit3c_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_lit1n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_lit2n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_lit3n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_cwdc_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_cwdn_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_minn_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_plantnuptake_pf,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_nleached_pf,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_ndenitri_pf,ierr)

    ! (i) BGC state variables: 3D subsurface PFLOTRAN ---to--- 3D subsurface CLM
    ! MPI Vecs for PFLOTRAN
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_sub,PETSC_DECIDE,clm_pf_idata%decomp_cpools_vr_lit1_pfp,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_lit1_pfp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_cpools_vr_lit2_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_cpools_vr_lit3_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_cpools_vr_cwd_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_cpools_vr_som1_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_cpools_vr_som2_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_cpools_vr_som3_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_cpools_vr_som4_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_npools_vr_lit1_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_npools_vr_lit2_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_npools_vr_lit3_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_npools_vr_cwd_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%sminn_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%smin_no3_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%smin_nh4_vr_pfp,ierr)
    ! Seq. Vecs for CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_sub,clm_pf_idata%decomp_cpools_vr_lit1_clms,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_lit1_clms,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_lit2_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_lit3_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_cwd_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_som1_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_som2_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_som3_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_som4_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_lit1_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_lit2_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_lit3_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_cwd_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%sminn_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%smin_no3_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%smin_nh4_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_lit1_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_lit2_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_lit3_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_cwd_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_som1_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_som2_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_som3_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_cpools_vr_som4_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_lit1_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_lit2_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_lit3_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_cwd_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%sminn_vr_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%smin_no3_vr_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%smin_nh4_vr_clms_prv,ierr)

    ! (iii) BGC flux variables: 3D subsurface PFLOTRAN ---to--- 3D subsurface CLM
    ! MPI Vecs for PFLOTRAN
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_sub,PETSC_DECIDE,clm_pf_idata%hrc_vr_pfp,ierr)
    call VecSet(clm_pf_idata%hrc_vr_pfp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%hrc_vr_pfp,clm_pf_idata%accextrn_vr_pfp,ierr)

    ! Seq. Vecs for CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_sub,clm_pf_idata%hrc_vr_clms,ierr)
    call VecSet(clm_pf_idata%hrc_vr_clms,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%hrc_vr_clms,clm_pf_idata%hrc_vr_clms_prv,ierr)
    call VecDuplicate(clm_pf_idata%hrc_vr_clms,clm_pf_idata%accextrn_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%hrc_vr_clms,clm_pf_idata%accextrn_vr_clms_prv,ierr)

    !---------------------------------------------

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
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)

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
    if(clm_pf_idata%rain_temp_clm  /= 0) call VecDestroy(clm_pf_idata%rain_temp_clm,ierr)
    if(clm_pf_idata%rain_temp_pf  /= 0) call VecDestroy(clm_pf_idata%rain_temp_pf,ierr)
    
    if(clm_pf_idata%gflux_subsurf_clm  /= 0) call VecDestroy(clm_pf_idata%gflux_subsurf_clm,ierr)
    if(clm_pf_idata%gflux_subsurf_pf  /= 0) call VecDestroy(clm_pf_idata%gflux_subsurf_pf,ierr)
    if(clm_pf_idata%gflux_surf_pf  /= 0) call VecDestroy(clm_pf_idata%gflux_surf_pf,ierr)

    if(clm_pf_idata%sat_clm  /= 0) call VecDestroy(clm_pf_idata%sat_clm,ierr)
    if(clm_pf_idata%sat_pf  /= 0) call VecDestroy(clm_pf_idata%sat_pf,ierr)

    if(clm_pf_idata%temp_clm  /= 0) call VecDestroy(clm_pf_idata%temp_clm,ierr)
    if(clm_pf_idata%temp_pf  /= 0) call VecDestroy(clm_pf_idata%temp_pf,ierr)

    if(clm_pf_idata%sat_ice_clm  /= 0) call VecDestroy(clm_pf_idata%sat_ice_clm,ierr)
    if(clm_pf_idata%sat_ice_pf  /= 0) call VecDestroy(clm_pf_idata%sat_ice_pf,ierr)

    if(clm_pf_idata%h2osfc_clm  /= 0) call VecDestroy(clm_pf_idata%h2osfc_clm,ierr)
    if(clm_pf_idata%h2osfc_pf  /= 0) call VecDestroy(clm_pf_idata%h2osfc_pf,ierr)

    if(clm_pf_idata%area_top_face_clm  /= 0) &
      call VecDestroy(clm_pf_idata%area_top_face_clm,ierr)
    if(clm_pf_idata%area_top_face_pf  /= 0) &
      call VecDestroy(clm_pf_idata%area_top_face_pf,ierr)

    ! -----------------------------------------------------------------------------------------------------------
    ! soil C/N pools (initial)
    if(clm_pf_idata%decomp_cpools_vr_lit1_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_clmp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_clmp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_clmp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_cwd_clmp  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_clmp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_clmp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_clmp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_clmp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_clmp,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_clmp,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_clmp,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_clmp,ierr)
    if(clm_pf_idata%decomp_npools_vr_cwd_clmp  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_clmp,ierr)
    if(clm_pf_idata%sminn_vr_clmp /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_clmp,ierr)
    if(clm_pf_idata%smin_no3_vr_clmp /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_clmp,ierr)
    if(clm_pf_idata%smin_nh4_vr_clmp /= 0) &
       call VecDestroy(clm_pf_idata%smin_nh4_vr_clmp,ierr)
    if(clm_pf_idata%soilpsi_clmp /= 0) &
       call VecDestroy(clm_pf_idata%soilpsi_clmp,ierr)
    if(clm_pf_idata%soillsat_clmp /= 0) &
       call VecDestroy(clm_pf_idata%soillsat_clmp,ierr)
    if(clm_pf_idata%soilisat_clmp /= 0) &
       call VecDestroy(clm_pf_idata%soilisat_clmp,ierr)
    if(clm_pf_idata%soilt_clmp /= 0) &
       call VecDestroy(clm_pf_idata%soilt_clmp,ierr)
!
    if(clm_pf_idata%decomp_cpools_vr_lit1_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_pfs,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_pfs,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_pfs,ierr)
    if(clm_pf_idata%decomp_cpools_vr_cwd_pfs  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_pfs,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_pfs,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_pfs,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_pfs,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_pfs,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_pfs,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_pfs,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_pfs,ierr)
    if(clm_pf_idata%decomp_npools_vr_cwd_pfs  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_pfs,ierr)
    if(clm_pf_idata%sminn_vr_pfs /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_pfs,ierr)
    if(clm_pf_idata%smin_no3_vr_pfs /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_pfs,ierr)
    if(clm_pf_idata%smin_nh4_vr_pfs /= 0) &
      call VecDestroy(clm_pf_idata%smin_nh4_vr_pfs,ierr)
    if(clm_pf_idata%soilpsi_pf /= 0) &
      call VecDestroy(clm_pf_idata%soilpsi_pf,ierr)
    if(clm_pf_idata%soillsat_pfs /= 0) &
      call VecDestroy(clm_pf_idata%soillsat_pfs,ierr)
    if(clm_pf_idata%soilisat_pfs /= 0) &
      call VecDestroy(clm_pf_idata%soilisat_pfs,ierr)
    if(clm_pf_idata%soilt_pfs /= 0) &
      call VecDestroy(clm_pf_idata%soilt_pfs,ierr)

    ! soil C/N fluxes at interface (source/sink)
    if(clm_pf_idata%rate_lit1c_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1c_clmp,ierr)
    if(clm_pf_idata%rate_lit2c_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2c_clmp,ierr)
    if(clm_pf_idata%rate_lit3c_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3c_clmp,ierr)
    if(clm_pf_idata%rate_cwdc_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdc_clmp,ierr)
    if(clm_pf_idata%rate_lit1n_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1n_clmp,ierr)
    if(clm_pf_idata%rate_lit2n_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2n_clmp,ierr)
    if(clm_pf_idata%rate_lit3n_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3n_clmp,ierr)
    if(clm_pf_idata%rate_cwdn_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdn_clmp,ierr)
    if(clm_pf_idata%rate_minn_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_minn_clmp,ierr)
    if(clm_pf_idata%rate_plantnuptake_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_plantnuptake_clmp,ierr)
    if(clm_pf_idata%rate_nleached_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_nleached_clmp,ierr)
    if(clm_pf_idata%rate_ndenitri_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_ndenitri_clmp,ierr)
    if(clm_pf_idata%rate_lit1c_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1c_pfs,ierr)
    if(clm_pf_idata%rate_lit2c_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2c_pfs,ierr)
    if(clm_pf_idata%rate_lit3c_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3c_pfs,ierr)
    if(clm_pf_idata%rate_cwdc_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdc_pfs,ierr)
    if(clm_pf_idata%rate_lit1n_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1n_pfs,ierr)
    if(clm_pf_idata%rate_lit2n_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2n_pfs,ierr)
    if(clm_pf_idata%rate_lit3n_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3n_pfs,ierr)
    if(clm_pf_idata%rate_cwdn_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdn_pfs,ierr)
    if(clm_pf_idata%rate_minn_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_minn_pfs,ierr)
    if(clm_pf_idata%rate_plantnuptake_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_plantnuptake_pf,ierr)
    if(clm_pf_idata%rate_nleached_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_nleached_pf,ierr)
    if(clm_pf_idata%rate_ndenitri_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_ndenitri_pf,ierr)

    ! update BGC states
    if(clm_pf_idata%decomp_cpools_vr_lit1_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_pfp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_pfp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_pfp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_cwd_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_pfp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_pfp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_pfp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_pfp,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_pfp,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_pfp,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_pfp,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_pfp,ierr)
    if(clm_pf_idata%decomp_npools_vr_cwd_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_pfp,ierr)
    if(clm_pf_idata%sminn_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_pfp,ierr)
    if(clm_pf_idata%smin_no3_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_pfp,ierr)
    if(clm_pf_idata%smin_nh4_vr_pfp /= 0) &
      call VecDestroy(clm_pf_idata%smin_nh4_vr_pfp,ierr)

    if(clm_pf_idata%decomp_cpools_vr_lit1_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_clms,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_clms,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_clms,ierr)
    if(clm_pf_idata%decomp_cpools_vr_cwd_clms  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_clms,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_clms,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_clms,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_clms,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_clms,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_clms,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_clms,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_clms,ierr)
    if(clm_pf_idata%decomp_npools_vr_cwd_clms  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_clms,ierr)
    if(clm_pf_idata%sminn_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_clms,ierr)
    if(clm_pf_idata%smin_no3_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_clms,ierr)
    if(clm_pf_idata%smin_nh4_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%smin_nh4_vr_clms,ierr)

    if(clm_pf_idata%decomp_cpools_vr_lit1_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_clms_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_clms_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_clms_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_cwd_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_clms_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_clms_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_clms_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_clms_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_clms_prv,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_clms_prv,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_clms_prv,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_clms_prv,ierr)
    if(clm_pf_idata%decomp_npools_vr_cwd_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_clms_prv,ierr)
    if(clm_pf_idata%sminn_vr_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_clms_prv,ierr)
    if(clm_pf_idata%smin_no3_vr_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_clms_prv,ierr)
    if(clm_pf_idata%smin_nh4_vr_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%smin_nh4_vr_clms_prv,ierr)

    ! update BGC fluxes
    if(clm_pf_idata%hrc_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%hrc_vr_pfp,ierr)
    if(clm_pf_idata%hrc_vr_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%hrc_vr_clms_prv,ierr)
    if(clm_pf_idata%hrc_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%hrc_vr_clms,ierr)

    if(clm_pf_idata%accextrn_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%accextrn_vr_pfp,ierr)
    if(clm_pf_idata%accextrn_vr_clms_prv /= 0) &
       call VecDestroy(clm_pf_idata%accextrn_vr_clms_prv,ierr)
    if(clm_pf_idata%accextrn_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%accextrn_vr_clms,ierr)

    !----------------------------------------------------------------------------------

  end subroutine CLMPFLOTRANIDataDestroy

end module clm_pflotran_interface_data
