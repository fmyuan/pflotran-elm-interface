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
  Vec :: soilpsi_clmg                   ! soil matric potential
  Vec :: soillsat_clmg                  ! soil liq. water saturation
  Vec :: soilisat_clmg                  ! soil ice water saturation
  Vec :: soilt_clmg                     ! soil temperature
  Vec :: soilpsi_pf
  Vec :: soillsat_pfl
  Vec :: soilisat_pfl
  Vec :: soilt_pfl
  ! initial ground/soil C/N pools from CLM (mpi) to PF (seq)
  Vec :: decomp_cpools_vr_lit1_clmg     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_clmg     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_clmg     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_cwd_clmg      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_clmg     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_clmg     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_clmg     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_clmg     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_clmg     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_clmg     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_clmg     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_cwd_clmg      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: sminn_vr_clmg                  ! (gN/m3) vertically-resolved soil mineral N
  Vec :: smin_no3_vr_clmg               ! (gN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_clmg               ! (gN/m3) vertically-resolved soil mineral NH4
  Vec :: decomp_cpools_vr_lit1_pfl      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_pfl      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_pfl      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_cwd_pfl       ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_pfl      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_pfl      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_pfl      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_pfl      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_pfl      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_pfl      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_pfl      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_cwd_pfl       ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: sminn_vr_pfl                   ! (gN/m3) vertically-resolved soil mineral N
  Vec :: smin_no3_vr_pfl                ! (gN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_pfl                ! (gN/m3) vertically-resolved soil mineral NH4
  ! time-varying ground/soil C/N rates from CLM (mpi) to PF (seq) (Units: moleC(N)/m3/s)
  Vec :: rate_lit1c_clmg
  Vec :: rate_lit2c_clmg
  Vec :: rate_lit3c_clmg
  Vec :: rate_lit1n_clmg
  Vec :: rate_lit2n_clmg
  Vec :: rate_lit3n_clmg
  Vec :: rate_cwdc_clmg
  Vec :: rate_cwdn_clmg
  Vec :: rate_minn_clmg
  Vec :: rate_plantnuptake_clmg
  Vec :: rate_nleached_clmg
  Vec :: rate_ndenitri_clmg
  Vec :: rate_lit1c_pfl
  Vec :: rate_lit2c_pfl
  Vec :: rate_lit3c_pfl
  Vec :: rate_lit1n_pfl
  Vec :: rate_lit2n_pfl
  Vec :: rate_lit3n_pfl
  Vec :: rate_cwdc_pfl
  Vec :: rate_cwdn_pfl
  Vec :: rate_minn_pfl
  Vec :: rate_plantnuptake_pf
  Vec :: rate_nleached_pf
  Vec :: rate_ndenitri_pf

  ! -----BGC vecs from PF (mpi, ghosted) to CLM (seq, local) --------------------
  Vec :: decomp_cpools_vr_lit1_pfg
  Vec :: decomp_cpools_vr_lit2_pfg
  Vec :: decomp_cpools_vr_lit3_pfg
  Vec :: decomp_cpools_vr_cwd_pfg
  Vec :: decomp_cpools_vr_som1_pfg
  Vec :: decomp_cpools_vr_som2_pfg
  Vec :: decomp_cpools_vr_som3_pfg
  Vec :: decomp_cpools_vr_som4_pfg
  Vec :: decomp_npools_vr_lit1_pfg
  Vec :: decomp_npools_vr_lit2_pfg
  Vec :: decomp_npools_vr_lit3_pfg
  Vec :: decomp_npools_vr_cwd_pfg
  Vec :: sminn_vr_pfg
  Vec :: smin_no3_vr_pfg
  Vec :: smin_nh4_vr_pfg
  !
  Vec :: decomp_cpools_vr_lit1_clml     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_clml     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_clml     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_cwd_clml      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_clml     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_clml     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_clml     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_clml     ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_clml     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_clml     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_clml     ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_cwd_clml      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: sminn_vr_clml                  ! (gN/m3) vertically-resolved soil mineral N
  Vec :: smin_no3_vr_clml               ! (gN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_clml               ! (gN/m3) vertically-resolved soil mineral NH4
  !
  Vec :: decomp_cpools_vr_lit1_clml_prv ! C/N states at previous time-step for pool-change calculation
  Vec :: decomp_cpools_vr_lit2_clml_prv
  Vec :: decomp_cpools_vr_lit3_clml_prv
  Vec :: decomp_cpools_vr_cwd_clml_prv
  Vec :: decomp_cpools_vr_som1_clml_prv
  Vec :: decomp_cpools_vr_som2_clml_prv
  Vec :: decomp_cpools_vr_som3_clml_prv
  Vec :: decomp_cpools_vr_som4_clml_prv
  Vec :: decomp_npools_vr_lit1_clml_prv
  Vec :: decomp_npools_vr_lit2_clml_prv
  Vec :: decomp_npools_vr_lit3_clml_prv
  Vec :: decomp_npools_vr_cwd_clml_prv
  Vec :: sminn_vr_clml_prv
  Vec :: smin_no3_vr_clml_prv
  Vec :: smin_nh4_vr_clml_prv

  ! 'hrc' is accumulative in 'PFLOTRAN', so needs previous time-step to calculate 'hr' fluxes for CLM-CN
  Vec :: hrc_vr_pfg                     ! (gC/m3) vertically-resolved soil heterotrophic respiration C
  Vec :: hrc_vr_clml_prv                ! (gC/m3) vertically-resolved soil heterotrophic respiration C at previous time-step
  Vec :: hrc_vr_clml                    ! (gC/m3) vertically-resolved soil heterotrophic respiration C
  ! 'accextrn' is accumulative N extract in 'PFLOTRAN', so needs previous time-step to calculate 'sminn_to_plant' fluxes for CLM-CN
  Vec :: accextrn_vr_pfg                ! (gN/m3) vertically-resolved root extraction N at previous time-step
  Vec :: accextrn_vr_clml_prv           ! (gN/m3) vertically-resolved root extraction N at previous time-step
  Vec :: accextrn_vr_clml               ! (gN/m3) vertically-resolved root extraction N at previous time-step

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
    clm_pf_idata%soilpsi_clmg  = 0
    clm_pf_idata%soillsat_clmg = 0
    clm_pf_idata%soilisat_clmg = 0
    clm_pf_idata%soilt_clmg    = 0
    clm_pf_idata%soilpsi_pf   = 0
    clm_pf_idata%soillsat_pfl = 0
    clm_pf_idata%soilisat_pfl = 0
    clm_pf_idata%soilt_pfl    = 0

    clm_pf_idata%decomp_cpools_vr_lit1_clmg = 0
    clm_pf_idata%decomp_cpools_vr_lit2_clmg = 0
    clm_pf_idata%decomp_cpools_vr_lit3_clmg = 0
    clm_pf_idata%decomp_cpools_vr_cwd_clmg  = 0
    clm_pf_idata%decomp_cpools_vr_som1_clmg = 0
    clm_pf_idata%decomp_cpools_vr_som2_clmg = 0
    clm_pf_idata%decomp_cpools_vr_som3_clmg = 0
    clm_pf_idata%decomp_cpools_vr_som4_clmg = 0
    clm_pf_idata%decomp_npools_vr_lit1_clmg = 0
    clm_pf_idata%decomp_npools_vr_lit2_clmg = 0
    clm_pf_idata%decomp_npools_vr_lit3_clmg = 0
    clm_pf_idata%decomp_npools_vr_cwd_clmg  = 0
    clm_pf_idata%sminn_vr_clmg         = 0
    clm_pf_idata%smin_no3_vr_clmg      = 0
    clm_pf_idata%smin_nh4_vr_clmg      = 0

    clm_pf_idata%decomp_cpools_vr_lit1_pfl = 0
    clm_pf_idata%decomp_cpools_vr_lit2_pfl = 0
    clm_pf_idata%decomp_cpools_vr_lit3_pfl = 0
    clm_pf_idata%decomp_cpools_vr_cwd_pfl  = 0
    clm_pf_idata%decomp_cpools_vr_som1_pfl = 0
    clm_pf_idata%decomp_cpools_vr_som2_pfl = 0
    clm_pf_idata%decomp_cpools_vr_som3_pfl = 0
    clm_pf_idata%decomp_cpools_vr_som4_pfl = 0
    clm_pf_idata%decomp_npools_vr_lit1_pfl = 0
    clm_pf_idata%decomp_npools_vr_lit2_pfl = 0
    clm_pf_idata%decomp_npools_vr_lit3_pfl = 0
    clm_pf_idata%decomp_npools_vr_cwd_pfl  = 0
    clm_pf_idata%sminn_vr_pfl          = 0
    clm_pf_idata%smin_no3_vr_pfl       = 0
    clm_pf_idata%smin_nh4_vr_pfl       = 0

    !ground/soil C/N rates as source/sink
    clm_pf_idata%rate_lit1c_clmg            = 0
    clm_pf_idata%rate_lit2c_clmg            = 0
    clm_pf_idata%rate_lit3c_clmg            = 0
    clm_pf_idata%rate_lit1n_clmg            = 0
    clm_pf_idata%rate_lit2n_clmg            = 0
    clm_pf_idata%rate_lit3n_clmg            = 0
    clm_pf_idata%rate_cwdc_clmg             = 0
    clm_pf_idata%rate_cwdn_clmg             = 0
    clm_pf_idata%rate_minn_clmg             = 0
    clm_pf_idata%rate_plantnuptake_clmg     = 0
    clm_pf_idata%rate_nleached_clmg         = 0
    clm_pf_idata%rate_ndenitri_clmg         = 0

    clm_pf_idata%rate_lit1c_pfl            = 0
    clm_pf_idata%rate_lit2c_pfl            = 0
    clm_pf_idata%rate_lit3c_pfl            = 0
    clm_pf_idata%rate_lit1n_pfl            = 0
    clm_pf_idata%rate_lit2n_pfl            = 0
    clm_pf_idata%rate_lit3n_pfl            = 0
    clm_pf_idata%rate_cwdc_pfl             = 0
    clm_pf_idata%rate_cwdn_pfl             = 0
    clm_pf_idata%rate_minn_pfl             = 0
    clm_pf_idata%rate_plantnuptake_pf      = 0
    clm_pf_idata%rate_nleached_pf          = 0
    clm_pf_idata%rate_ndenitri_pf          = 0

    ! for updating bgc states
    clm_pf_idata%decomp_cpools_vr_lit1_pfg = 0
    clm_pf_idata%decomp_cpools_vr_lit2_pfg = 0
    clm_pf_idata%decomp_cpools_vr_lit3_pfg = 0
    clm_pf_idata%decomp_cpools_vr_cwd_pfg  = 0
    clm_pf_idata%decomp_cpools_vr_som1_pfg = 0
    clm_pf_idata%decomp_cpools_vr_som2_pfg = 0
    clm_pf_idata%decomp_cpools_vr_som3_pfg = 0
    clm_pf_idata%decomp_cpools_vr_som4_pfg = 0
    clm_pf_idata%decomp_npools_vr_lit1_pfg = 0
    clm_pf_idata%decomp_npools_vr_lit2_pfg = 0
    clm_pf_idata%decomp_npools_vr_lit3_pfg = 0
    clm_pf_idata%decomp_npools_vr_cwd_pfg  = 0
    clm_pf_idata%sminn_vr_pfg          = 0
    clm_pf_idata%smin_no3_vr_pfg       = 0
    clm_pf_idata%smin_nh4_vr_pfg       = 0
    clm_pf_idata%decomp_cpools_vr_lit1_clml = 0
    clm_pf_idata%decomp_cpools_vr_lit2_clml = 0
    clm_pf_idata%decomp_cpools_vr_lit3_clml = 0
    clm_pf_idata%decomp_cpools_vr_cwd_clml  = 0
    clm_pf_idata%decomp_cpools_vr_som1_clml = 0
    clm_pf_idata%decomp_cpools_vr_som2_clml = 0
    clm_pf_idata%decomp_cpools_vr_som3_clml = 0
    clm_pf_idata%decomp_cpools_vr_som4_clml = 0
    clm_pf_idata%decomp_npools_vr_lit1_clml = 0
    clm_pf_idata%decomp_npools_vr_lit2_clml = 0
    clm_pf_idata%decomp_npools_vr_lit3_clml = 0
    clm_pf_idata%decomp_npools_vr_cwd_clml  = 0
    clm_pf_idata%sminn_vr_clml         = 0
    clm_pf_idata%smin_no3_vr_clml      = 0
    clm_pf_idata%smin_nh4_vr_clml      = 0
    clm_pf_idata%decomp_cpools_vr_lit1_clml_prv = 0
    clm_pf_idata%decomp_cpools_vr_lit2_clml_prv = 0
    clm_pf_idata%decomp_cpools_vr_lit3_clml_prv = 0
    clm_pf_idata%decomp_cpools_vr_cwd_clml_prv  = 0
    clm_pf_idata%decomp_cpools_vr_som1_clml_prv = 0
    clm_pf_idata%decomp_cpools_vr_som2_clml_prv = 0
    clm_pf_idata%decomp_cpools_vr_som3_clml_prv = 0
    clm_pf_idata%decomp_cpools_vr_som4_clml_prv = 0
    clm_pf_idata%decomp_npools_vr_lit1_clml_prv = 0
    clm_pf_idata%decomp_npools_vr_lit2_clml_prv = 0
    clm_pf_idata%decomp_npools_vr_lit3_clml_prv = 0
    clm_pf_idata%decomp_npools_vr_cwd_clml_prv  = 0
    clm_pf_idata%sminn_vr_clml_prv         = 0
    clm_pf_idata%smin_no3_vr_clml_prv      = 0
    clm_pf_idata%smin_nh4_vr_clml_prv      = 0

    ! for soil hr calculation
    clm_pf_idata%hrc_vr_pfg            = 0
    clm_pf_idata%hrc_vr_clml_prv       = 0
    clm_pf_idata%hrc_vr_clml           = 0

    ! for root N extraction calculation
    clm_pf_idata%accextrn_vr_pfg       = 0
    clm_pf_idata%accextrn_vr_clml_prv  = 0
    clm_pf_idata%accextrn_vr_clml      = 0

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
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%decomp_cpools_vr_lit1_clmg,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_lit1_clmg,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_cpools_vr_lit2_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_cpools_vr_lit3_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_cpools_vr_cwd_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_cpools_vr_som1_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_cpools_vr_som2_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_cpools_vr_som3_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_cpools_vr_som4_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_npools_vr_lit1_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_npools_vr_lit2_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_npools_vr_lit3_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%decomp_npools_vr_cwd_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%sminn_vr_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%smin_no3_vr_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%smin_nh4_vr_clmg,ierr)

    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%soilpsi_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%soillsat_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%soilisat_clmg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmg,clm_pf_idata%soilt_clmg,ierr)

    ! Seq. Vecs for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%decomp_cpools_vr_lit1_pfl,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_lit1_pfl,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_cpools_vr_lit2_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_cpools_vr_lit3_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_cpools_vr_cwd_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_cpools_vr_som1_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_cpools_vr_som2_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_cpools_vr_som3_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_cpools_vr_som4_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_npools_vr_lit1_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_npools_vr_lit2_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_npools_vr_lit3_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%decomp_npools_vr_cwd_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%sminn_vr_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%smin_no3_vr_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%smin_nh4_vr_pfl,ierr)

    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%soilpsi_pf,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%soillsat_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%soilisat_pfl,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfl,clm_pf_idata%soilt_pfl,ierr)

    ! (ii) BGC interface source/sink (rates): 3D subsurface CLM ---to--- 3D subsurface PFLOTRAN
    ! MPI Vecs for CLM
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%rate_lit1c_clmg,ierr)
    call VecSet(clm_pf_idata%rate_lit1c_clmg,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_lit2c_clmg,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_lit3c_clmg,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_lit1n_clmg,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_lit2n_clmg,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_lit3n_clmg,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_cwdc_clmg,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_cwdn_clmg,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_minn_clmg,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_plantnuptake_clmg,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_nleached_clmg,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmg,clm_pf_idata%rate_ndenitri_clmg,ierr)

    ! Seq. Vecs for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%rate_lit1c_pfl,ierr)
    call VecSet(clm_pf_idata%rate_lit1c_pfl,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_lit2c_pfl,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_lit3c_pfl,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_lit1n_pfl,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_lit2n_pfl,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_lit3n_pfl,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_cwdc_pfl,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_cwdn_pfl,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_minn_pfl,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_plantnuptake_pf,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_nleached_pf,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfl,clm_pf_idata%rate_ndenitri_pf,ierr)

    ! (i) BGC state variables: 3D subsurface PFLOTRAN ---to--- 3D subsurface CLM
    ! MPI Vecs for PFLOTRAN
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_sub,PETSC_DECIDE,clm_pf_idata%decomp_cpools_vr_lit1_pfg,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_lit1_pfg,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_cpools_vr_lit2_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_cpools_vr_lit3_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_cpools_vr_cwd_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_cpools_vr_som1_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_cpools_vr_som2_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_cpools_vr_som3_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_cpools_vr_som4_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_npools_vr_lit1_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_npools_vr_lit2_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_npools_vr_lit3_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%decomp_npools_vr_cwd_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%sminn_vr_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%smin_no3_vr_pfg,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfg,clm_pf_idata%smin_nh4_vr_pfg,ierr)
    ! Seq. Vecs for CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_sub,clm_pf_idata%decomp_cpools_vr_lit1_clml,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_lit1_clml,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_lit2_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_lit3_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_cwd_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_som1_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_som2_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_som3_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_som4_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_npools_vr_lit1_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_npools_vr_lit2_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_npools_vr_lit3_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_npools_vr_cwd_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%sminn_vr_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%smin_no3_vr_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%smin_nh4_vr_clml,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_lit1_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_lit2_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_lit3_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_cwd_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_som1_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_som2_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_som3_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_cpools_vr_som4_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_npools_vr_lit1_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_npools_vr_lit2_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_npools_vr_lit3_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%decomp_npools_vr_cwd_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%sminn_vr_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%smin_no3_vr_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clml,clm_pf_idata%smin_nh4_vr_clml_prv,ierr)

    ! (iii) BGC flux variables: 3D subsurface PFLOTRAN ---to--- 3D subsurface CLM
    ! MPI Vecs for PFLOTRAN
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_sub,PETSC_DECIDE,clm_pf_idata%hrc_vr_pfg,ierr)
    call VecSet(clm_pf_idata%hrc_vr_pfg,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%hrc_vr_pfg,clm_pf_idata%accextrn_vr_pfg,ierr)

    ! Seq. Vecs for CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_sub,clm_pf_idata%hrc_vr_clml,ierr)
    call VecSet(clm_pf_idata%hrc_vr_clml,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%hrc_vr_clml,clm_pf_idata%hrc_vr_clml_prv,ierr)
    call VecDuplicate(clm_pf_idata%hrc_vr_clml,clm_pf_idata%accextrn_vr_clml,ierr)
    call VecDuplicate(clm_pf_idata%hrc_vr_clml,clm_pf_idata%accextrn_vr_clml_prv,ierr)

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
    if(clm_pf_idata%decomp_cpools_vr_lit1_clmg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_clmg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_clmg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_clmg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_clmg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_clmg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_cwd_clmg  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_clmg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_clmg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_clmg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_clmg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_clmg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_clmg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_clmg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_clmg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_clmg,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_clmg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_clmg,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_clmg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_clmg,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_clmg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_clmg,ierr)
    if(clm_pf_idata%decomp_npools_vr_cwd_clmg  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_clmg,ierr)
    if(clm_pf_idata%sminn_vr_clmg /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_clmg,ierr)
    if(clm_pf_idata%smin_no3_vr_clmg /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_clmg,ierr)
    if(clm_pf_idata%smin_nh4_vr_clmg /= 0) &
       call VecDestroy(clm_pf_idata%smin_nh4_vr_clmg,ierr)
    if(clm_pf_idata%soilpsi_clmg /= 0) &
       call VecDestroy(clm_pf_idata%soilpsi_clmg,ierr)
    if(clm_pf_idata%soillsat_clmg /= 0) &
       call VecDestroy(clm_pf_idata%soillsat_clmg,ierr)
    if(clm_pf_idata%soilisat_clmg /= 0) &
       call VecDestroy(clm_pf_idata%soilisat_clmg,ierr)
    if(clm_pf_idata%soilt_clmg /= 0) &
       call VecDestroy(clm_pf_idata%soilt_clmg,ierr)
!
    if(clm_pf_idata%decomp_cpools_vr_lit1_pfl /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_pfl,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_pfl /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_pfl,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_pfl /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_pfl,ierr)
    if(clm_pf_idata%decomp_cpools_vr_cwd_pfl  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_pfl,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_pfl /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_pfl,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_pfl /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_pfl,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_pfl /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_pfl,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_pfl /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_pfl,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_pfl /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_pfl,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_pfl /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_pfl,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_pfl /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_pfl,ierr)
    if(clm_pf_idata%decomp_npools_vr_cwd_pfl  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_pfl,ierr)
    if(clm_pf_idata%sminn_vr_pfl /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_pfl,ierr)
    if(clm_pf_idata%smin_no3_vr_pfl /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_pfl,ierr)
    if(clm_pf_idata%smin_nh4_vr_pfl /= 0) &
      call VecDestroy(clm_pf_idata%smin_nh4_vr_pfl,ierr)
    if(clm_pf_idata%soilpsi_pf /= 0) &
      call VecDestroy(clm_pf_idata%soilpsi_pf,ierr)
    if(clm_pf_idata%soillsat_pfl /= 0) &
      call VecDestroy(clm_pf_idata%soillsat_pfl,ierr)
    if(clm_pf_idata%soilisat_pfl /= 0) &
      call VecDestroy(clm_pf_idata%soilisat_pfl,ierr)
    if(clm_pf_idata%soilt_pfl /= 0) &
      call VecDestroy(clm_pf_idata%soilt_pfl,ierr)

    ! soil C/N fluxes at interface (source/sink)
    if(clm_pf_idata%rate_lit1c_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1c_clmg,ierr)
    if(clm_pf_idata%rate_lit2c_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2c_clmg,ierr)
    if(clm_pf_idata%rate_lit3c_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3c_clmg,ierr)
    if(clm_pf_idata%rate_cwdc_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdc_clmg,ierr)
    if(clm_pf_idata%rate_lit1n_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1n_clmg,ierr)
    if(clm_pf_idata%rate_lit2n_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2n_clmg,ierr)
    if(clm_pf_idata%rate_lit3n_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3n_clmg,ierr)
    if(clm_pf_idata%rate_cwdn_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdn_clmg,ierr)
    if(clm_pf_idata%rate_minn_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_minn_clmg,ierr)
    if(clm_pf_idata%rate_plantnuptake_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_plantnuptake_clmg,ierr)
    if(clm_pf_idata%rate_nleached_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_nleached_clmg,ierr)
    if(clm_pf_idata%rate_ndenitri_clmg /= 0) &
       call VecDestroy(clm_pf_idata%rate_ndenitri_clmg,ierr)
    if(clm_pf_idata%rate_lit1c_pfl /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1c_pfl,ierr)
    if(clm_pf_idata%rate_lit2c_pfl /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2c_pfl,ierr)
    if(clm_pf_idata%rate_lit3c_pfl /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3c_pfl,ierr)
    if(clm_pf_idata%rate_cwdc_pfl /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdc_pfl,ierr)
    if(clm_pf_idata%rate_lit1n_pfl /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1n_pfl,ierr)
    if(clm_pf_idata%rate_lit2n_pfl /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2n_pfl,ierr)
    if(clm_pf_idata%rate_lit3n_pfl /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3n_pfl,ierr)
    if(clm_pf_idata%rate_cwdn_pfl /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdn_pfl,ierr)
    if(clm_pf_idata%rate_minn_pfl /= 0) &
       call VecDestroy(clm_pf_idata%rate_minn_pfl,ierr)
    if(clm_pf_idata%rate_plantnuptake_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_plantnuptake_pf,ierr)
    if(clm_pf_idata%rate_nleached_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_nleached_pf,ierr)
    if(clm_pf_idata%rate_ndenitri_pf /= 0) &
       call VecDestroy(clm_pf_idata%rate_ndenitri_pf,ierr)

    ! update BGC states
    if(clm_pf_idata%decomp_cpools_vr_lit1_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_pfg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_pfg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_pfg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_cwd_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_pfg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_pfg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_pfg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_pfg,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_pfg,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_pfg,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_pfg,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_pfg,ierr)
    if(clm_pf_idata%decomp_npools_vr_cwd_pfg /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_pfg,ierr)
    if(clm_pf_idata%sminn_vr_pfg /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_pfg,ierr)
    if(clm_pf_idata%smin_no3_vr_pfg /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_pfg,ierr)
    if(clm_pf_idata%smin_nh4_vr_pfg /= 0) &
      call VecDestroy(clm_pf_idata%smin_nh4_vr_pfg,ierr)

    if(clm_pf_idata%decomp_cpools_vr_lit1_clml /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_clml,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_clml /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_clml,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_clml /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_clml,ierr)
    if(clm_pf_idata%decomp_cpools_vr_cwd_clml  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_clml,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_clml /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_clml,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_clml /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_clml,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_clml /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_clml,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_clml /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_clml,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_clml /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_clml,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_clml /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_clml,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_clml /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_clml,ierr)
    if(clm_pf_idata%decomp_npools_vr_cwd_clml  /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_clml,ierr)
    if(clm_pf_idata%sminn_vr_clml /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_clml,ierr)
    if(clm_pf_idata%smin_no3_vr_clml /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_clml,ierr)
    if(clm_pf_idata%smin_nh4_vr_clml /= 0) &
       call VecDestroy(clm_pf_idata%smin_nh4_vr_clml,ierr)

    if(clm_pf_idata%decomp_cpools_vr_lit1_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit1_clml_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit2_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit2_clml_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_lit3_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_lit3_clml_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_cwd_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_cwd_clml_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som1_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som1_clml_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som2_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som2_clml_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som3_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som3_clml_prv,ierr)
    if(clm_pf_idata%decomp_cpools_vr_som4_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_clml_prv,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit1_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit1_clml_prv,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit2_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit2_clml_prv,ierr)
    if(clm_pf_idata%decomp_npools_vr_lit3_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_lit3_clml_prv,ierr)
    if(clm_pf_idata%decomp_npools_vr_cwd_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_cwd_clml_prv,ierr)
    if(clm_pf_idata%sminn_vr_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%sminn_vr_clml_prv,ierr)
    if(clm_pf_idata%smin_no3_vr_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_clml_prv,ierr)
    if(clm_pf_idata%smin_nh4_vr_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%smin_nh4_vr_clml_prv,ierr)

    ! update BGC fluxes
    if(clm_pf_idata%hrc_vr_pfg /= 0) &
       call VecDestroy(clm_pf_idata%hrc_vr_pfg,ierr)
    if(clm_pf_idata%hrc_vr_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%hrc_vr_clml_prv,ierr)
    if(clm_pf_idata%hrc_vr_clml /= 0) &
       call VecDestroy(clm_pf_idata%hrc_vr_clml,ierr)

    if(clm_pf_idata%accextrn_vr_pfg /= 0) &
       call VecDestroy(clm_pf_idata%accextrn_vr_pfg,ierr)
    if(clm_pf_idata%accextrn_vr_clml_prv /= 0) &
       call VecDestroy(clm_pf_idata%accextrn_vr_clml_prv,ierr)
    if(clm_pf_idata%accextrn_vr_clml /= 0) &
       call VecDestroy(clm_pf_idata%accextrn_vr_clml,ierr)

    !----------------------------------------------------------------------------------

  end subroutine CLMPFLOTRANIDataDestroy

end module clm_pflotran_interface_data
