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
  Vec :: watfc_clm
  Vec :: bsw_clm
  Vec :: press_clm

  ! pflotran has rock (solid) density, not bulk density,
  ! there could be issues with influence of soil organic content, ice on bulk density
  ! this will be temporaily used for calculation needed in denitrification calculation
  ! following clm4.5.35
  Vec :: bulkdensity_dry_clm

  ! Local for PFLOTRAN - seq. vec
  Vec :: hksat_x_pf
  Vec :: hksat_y_pf
  Vec :: hksat_z_pf
  Vec :: sucsat_pf
  Vec :: watsat_pf
  Vec :: watfc_pf
  Vec :: bsw_pf
  !Vec :: press_pf
  Vec :: bulkdensity_dry_pf   

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

  Vec :: eff_therm_cond_clm ! seq vec
  Vec :: eff_therm_cond_pf  ! mpi vec

  ! Number of cells for the 3D subsurface domain
  PetscInt :: nlclm_sub ! num of local clm cells
  PetscInt :: ngclm_sub ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_sub  ! num of local pflotran cells
  PetscInt :: ngpf_sub  ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the surface of the 3D subsurface domain
  PetscInt :: nlclm_2dsub  ! num of local clm cells
  PetscInt :: ngclm_2dsub  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2dsub   ! num of local pflotran cells
  PetscInt :: ngpf_2dsub   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  PetscInt :: nlclm_bottom  ! num of local clm cells
  PetscInt :: ngclm_bottom  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_bottom   ! num of local pflotran cells
  PetscInt :: ngpf_bottom   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the 2D surface domain
  PetscInt :: nlclm_srf  ! num of local clm cells
  PetscInt :: ngclm_srf  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_srf   ! num of local pflotran cells
  PetscInt :: ngpf_srf   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  PetscInt :: nzclm_mapped ! num of CLM soil layers that are mapped

  !-------------------------------------------------------------------------------------
  ! note: mapping ONLY can do from MPI vecs to Seq. vecs now
  ! -----BGC vecs from CLM to PF --------------------
  ! TH properties
  Vec :: sr_clmp    ! MVM soil hydraulic properties
  Vec :: lamda_clmp
  Vec :: alpha_clmp
  Vec :: pcwmax_clmp
  Vec :: porosity_clmp
  Vec :: press_ref_clmp
  Vec :: zsoi_clmp
  Vec :: sr_pfs    ! MVM soil hydraulic properties
  Vec :: lamda_pfs
  Vec :: alpha_pfs
  Vec :: pcwmax_pfs
  Vec :: porosity_pfs
  Vec :: press_ref_pfs
  Vec :: zsoi_pfs

  ! TH state vecs from CLM (mpi) to PF (seq) for bgc
  Vec :: press_clmp                     ! water pressure head (Pa)
  Vec :: soilpsi_clmp                   ! soil matric potential (Pa)
  Vec :: soillsat_clmp                  ! soil liq. water saturation (0 - 1)
  Vec :: soilisat_clmp                  ! soil ice water saturation (0 - 1)
  Vec :: soilt_clmp                     ! soil temperature (degC)
  Vec :: press_pfs
  Vec :: soilpsi_pfs
  Vec :: soillsat_pfs
  Vec :: soilisat_pfs
  Vec :: soilt_pfs
  ! initial ground/soil C/N pools from CLM (mpi) to PF (seq)
  Vec :: decomp_cpools_vr_lit1_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_cwd_clmp      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_clmp     ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_clmp     ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_clmp     ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_cwd_clmp      ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_som1_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_som2_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_som3_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_som4_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: smin_no3_vr_clmp               ! (moleN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_clmp               ! (moleN/m3) vertically-resolved soil mineral NH4
  Vec :: decomp_cpools_vr_lit1_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_cwd_pfs       ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_pfs      ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_pfs      ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_pfs      ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_cwd_pfs       ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_som1_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_som2_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_som3_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_som4_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: smin_no3_vr_pfs                ! (moleN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_pfs                ! (moleN/m3) vertically-resolved soil mineral NH4
  ! time-varying ground/soil C/N rates from CLM (mpi) to PF (seq) (Units: moleC(N)/m3/s)
  Vec :: rate_lit1c_clmp
  Vec :: rate_lit2c_clmp
  Vec :: rate_lit3c_clmp
  Vec :: rate_cwdc_clmp
  Vec :: rate_som1c_clmp
  Vec :: rate_som2c_clmp
  Vec :: rate_som3c_clmp
  Vec :: rate_som4c_clmp
  Vec :: rate_lit1n_clmp
  Vec :: rate_lit2n_clmp
  Vec :: rate_lit3n_clmp
  Vec :: rate_cwdn_clmp
  Vec :: rate_som1n_clmp
  Vec :: rate_som2n_clmp
  Vec :: rate_som3n_clmp
  Vec :: rate_som4n_clmp
  Vec :: rate_smin_no3_clmp
  Vec :: rate_smin_nh4_clmp
  Vec :: rate_plantndemand_clmp
  Vec :: rate_lit1c_pfs
  Vec :: rate_lit2c_pfs
  Vec :: rate_lit3c_pfs
  Vec :: rate_cwdc_pfs
  Vec :: rate_som1c_pfs
  Vec :: rate_som2c_pfs
  Vec :: rate_som3c_pfs
  Vec :: rate_som4c_pfs
  Vec :: rate_lit1n_pfs
  Vec :: rate_lit2n_pfs
  Vec :: rate_lit3n_pfs
  Vec :: rate_cwdn_pfs
  Vec :: rate_som1n_pfs
  Vec :: rate_som2n_pfs
  Vec :: rate_som3n_pfs
  Vec :: rate_som4n_pfs
  Vec :: rate_smin_no3_pfs
  Vec :: rate_smin_nh4_pfs
  Vec :: rate_plantndemand_pfs

  ! BC: water pressure (Pa) on the top/bottom of 3-D subsurface domain as boundary conditions from CLM to PF
  Vec :: press_subsurf_clmp    ! mpi vec
  Vec :: press_subbase_clmp    ! mpi vec
  Vec :: press_maxponding_clmp ! mpi vec
  Vec :: press_subsurf_pfs     ! seq vec
  Vec :: press_subbase_pfs     ! seq vec
  Vec :: press_maxponding_pfs  ! seq vec
  ! OR, BC: water infiltration/recharge(drainage) (mH2O/sec) on the top/bottom of 3-D subsurface domain as boundary conditions from CLM to PF
  PetscBool :: topbc_seepage   ! use 'seepage' top-BC to calculate surface-overflow by PFLOTRAN
  Vec :: qflux_subsurf_clmp    ! mpi vec
  Vec :: qflux_subbase_clmp    ! mpi vec
  Vec :: qflux_subsurf_pfs     ! seq vec
  Vec :: qflux_subbase_pfs     ! seq vec

  ! if ground temperature at the subsurface interface is known
  Vec :: gtemp_subsurf_clmp  ! mpi vec
  Vec :: gtemp_subsurf_pfs   ! seq vec
  ! if bottom heat flux at the subsurface interface is known
  Vec :: gflux_subbase_clmp  ! mpi vec
  Vec :: gflux_subbase_pfs   ! seq vec
  ! bottom temperature at the subsurface domain interface is known
  Vec :: gtemp_subbase_clmp  ! mpi vec
  Vec :: gtemp_subbase_pfs   ! seq vec

  ! -----BGC vecs from PF (mpi, ghosted) to CLM (seq, local) --------------------
  ! TH properties (useful to do some calculation in the interface)
  Vec :: sr_pcwmax_pfp     ! MVM soil hydraulic properties
  Vec :: pcwmax_pfp
  Vec :: porosity_pfp
  Vec :: sr_pcwmax_clms    ! MVM soil hydraulic properties
  Vec :: pcwmax_clms
  Vec :: porosity_clms

  PetscReal :: pressure_reference

  ! TH state vecs from PF (mpi) to CLM (seq)
  Vec :: press_pfp                     ! water pressure head (Pa)
  Vec :: soilpsi_pfp                   ! soil matric potential (Pa) (negative)
  Vec :: soillsat_pfp                  ! soil liq. water saturation (0-1)
  Vec :: soilisat_pfp                  ! soil ice water saturation (0-1)
  Vec :: soilt_pfp                     ! soil temperature (degC)
  Vec :: press_clms
  Vec :: soilpsi_clms
  Vec :: soillsat_clms
  Vec :: soilisat_clms
  Vec :: soilt_clms

  ! BGC state variables
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
  Vec :: decomp_npools_vr_som1_pfp
  Vec :: decomp_npools_vr_som2_pfp
  Vec :: decomp_npools_vr_som3_pfp
  Vec :: decomp_npools_vr_som4_pfp
  Vec :: smin_no3_vr_pfp
  Vec :: smin_nh4_vr_pfp                ! (moleN/m3) vertically-resolved total soil mineral NH4 (incl. absorbed)
  Vec :: smin_nh4sorb_vr_pfp            ! (moleN/m3) vertically-resolved absorbed NH4-N
  !
  Vec :: decomp_cpools_vr_lit1_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit2_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_lit3_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_cwd_clms      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som1_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som2_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som3_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_cpools_vr_som4_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_lit1_clms     ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit2_clms     ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_lit3_clms     ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_cwd_clms      ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: decomp_npools_vr_som1_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_som2_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_som3_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_som4_clms     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: smin_no3_vr_clms               ! (moleN/m3) vertically-resolved total soil mineral NO3
  Vec :: smin_nh4_vr_clms               ! (moleN/m3) vertically-resolved total soil mineral NH4 (incl. absorbed)
  Vec :: smin_nh4sorb_vr_clms           ! (moleN/m3) vertically-resolved absorbed NH4-N
  !
  ! 'accextrn' is accumulative N extract by plant roots in 'PFLOTRAN' within a CLM timestep
  Vec :: accextrn_vr_pfp                ! (moleN/m3) vertically-resolved root extraction N
  Vec :: accextrn_vr_clms               ! (moleN/m3) vertically-resolved root extraction N

  ! gases in water (aqueous solution of gases)
  ! gases species is accumulative in 'PFLOTRAN', so needs to calculate their fluxes in the CLM-PF interface and reset back to PFLOTRAN
  Vec :: gco2_vr_pfp                   ! (moleC/m3) vertically-resolved soil CO2 C
  Vec :: gco2_vr_clms                  ! (moleC/m3) vertically-resolved soil CO2 C
  Vec :: gco2_vr_clmp                  ! (moleC/m3) vertically-resolved soil CO2 C, after gas emission
  Vec :: gco2_vr_pfs                   ! (moleC/m3) vertically-resolved soil CO2 C, after gas emission

  Vec :: gn2_vr_pfp                    ! (moleN/m3) vertically-resolved N2-N
  Vec :: gn2_vr_clms                   ! (moleN/m3) vertically-resolved N2-N
  Vec :: gn2_vr_clmp                   ! (moleN/m3) vertically-resolved N2-N, after gas emission
  Vec :: gn2_vr_pfs                    ! (moleN/m3) vertically-resolved N2-N, after gas emission

  Vec :: gn2o_vr_pfp                   ! (moleN/m3) vertically-resolved N2O-N
  Vec :: gn2o_vr_clms                  ! (moleN/m3) vertically-resolved N2O-N
  Vec :: gn2o_vr_clmp                  ! (moleN/m3) vertically-resolved N2O-N, after gas emission
  Vec :: gn2o_vr_pfs                   ! (moleN/m3) vertically-resolved N2O-N, after gas emission

  ! some tracking variables from PFLOTRAN bgc to obtain reaction flux rates which needed by CLM
  Vec :: acchr_vr_pfp                 ! (moleC/m3/timestep) vertically-resolved heterotrophic resp. C from decomposition
  Vec :: acchr_vr_clms                ! (moleC/m3/timestep) vertically-resolved heterotrophic resp. C from decomposition

  Vec :: accnmin_vr_pfp                ! (moleN/m3/timestep) vertically-resolved N mineralization
  Vec :: accnmin_vr_clms               ! (moleN/m3/timestep) vertically-resolved N mineralization

  Vec :: accnimm_vr_pfp                ! (moleN/m3/timestep) vertically-resolved N immoblization
  Vec :: accnimm_vr_clms               ! (moleN/m3/timestep) vertically-resolved N immoblization

  Vec :: accngasmin_vr_pfp              ! (moleN/m3/timestep) vertically-resolved N2O-N from mineralization
  Vec :: accngasmin_vr_clms             ! (moleN/m3/timestep) vertically-resolved N2O-N from mineralization

  Vec :: accngasnitr_vr_pfp             ! (moleN/m3/timestep) vertically-resolved N2O-N from nitrification
  Vec :: accngasnitr_vr_clms            ! (moleN/m3/timestep) vertically-resolved N2O-N from nitrification

  Vec :: accngasdeni_vr_pfp             ! (moleN/m3/timestep) vertically-resolved N2O-N from denitrification
  Vec :: accngasdeni_vr_clms            ! (moleN/m3/timestep) vertically-resolved N2O-N from denitrification

  ! actual mass water flow rate (kgH2O/sec) through the top/bottom BC of 3-D subsurface domain
  ! (+ in, - out)
  Vec :: qinfl_subsurf_pfp    ! mpi vec: actual infiltration (+)
  Vec :: qinfl_subsurf_clms   ! seq vec
  Vec :: qsurf_subsurf_pfp    ! mpi vec: actual overland flow - potential-actual infiltration or water upwarding (-)
  Vec :: qsurf_subsurf_clms   ! seq vec
  Vec :: qflux_subbase_pfp    ! mpi vec
  Vec :: qflux_subbase_clms   ! seq vec
  ! actual aqeuous N mass flow rate(moleN/m2/sec) at the top (runoff)/bottom (leaching) of 3-D subsurface domain
  ! (+ in, - out)
  Vec :: f_nh4_subsurf_pfp    ! mpi vec
  Vec :: f_nh4_subsurf_clms   ! seq vec
  Vec :: f_nh4_subbase_pfp    ! mpi vec
  Vec :: f_nh4_subbase_clms   ! seq vec
  Vec :: f_no3_subsurf_pfp    ! mpi vec
  Vec :: f_no3_subsurf_clms   ! seq vec
  Vec :: f_no3_subbase_pfp    ! mpi vec
  Vec :: f_no3_subbase_clms   ! seq vec

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
  ! 
  
    implicit none

    clm_pf_idata%nlclm_sub = 0
    clm_pf_idata%ngclm_sub = 0
    clm_pf_idata%nlpf_sub = 0
    clm_pf_idata%ngpf_sub = 0

    clm_pf_idata%nlclm_2dsub = 0
    clm_pf_idata%ngclm_2dsub = 0
    clm_pf_idata%nlpf_2dsub = 0
    clm_pf_idata%ngpf_2dsub = 0

    clm_pf_idata%nlclm_bottom = 0
    clm_pf_idata%ngclm_bottom = 0
    clm_pf_idata%nlpf_bottom = 0
    clm_pf_idata%ngpf_bottom = 0

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
    clm_pf_idata%watfc_clm = 0
    clm_pf_idata%bulkdensity_dry_clm = 0

    clm_pf_idata%hksat_x_pf = 0
    clm_pf_idata%hksat_y_pf = 0
    clm_pf_idata%hksat_z_pf = 0
    clm_pf_idata%sucsat_pf = 0
    clm_pf_idata%watsat_pf = 0
    clm_pf_idata%bsw_pf = 0
    clm_pf_idata%watfc_pf = 0
    clm_pf_idata%bulkdensity_dry_pf = 0

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
    clm_pf_idata%sr_clmp       = 0
    clm_pf_idata%lamda_clmp    = 0
    clm_pf_idata%alpha_clmp    = 0
    clm_pf_idata%pcwmax_clmp   = 0
    clm_pf_idata%porosity_clmp = 0
    clm_pf_idata%press_ref_clmp= 0
    clm_pf_idata%zsoi_clmp     = 0
    clm_pf_idata%sr_pfs         = 0
    clm_pf_idata%lamda_pfs      = 0
    clm_pf_idata%alpha_pfs      = 0
    clm_pf_idata%pcwmax_pfs     = 0
    clm_pf_idata%porosity_pfs   = 0
    clm_pf_idata%press_ref_pfs  = 0
    clm_pf_idata%zsoi_pfs       = 0

    clm_pf_idata%eff_therm_cond_clm = 0
    clm_pf_idata%eff_therm_cond_pf = 0

    clm_pf_idata%nzclm_mapped = 0

   ! soil TH and C/N pools
    clm_pf_idata%press_clmp    = 0
    clm_pf_idata%soilpsi_clmp  = 0
    clm_pf_idata%soillsat_clmp = 0
    clm_pf_idata%soilisat_clmp = 0
    clm_pf_idata%soilt_clmp    = 0
    clm_pf_idata%press_pfs      = 0
    clm_pf_idata%soilpsi_pfs    = 0
    clm_pf_idata%soillsat_pfs   = 0
    clm_pf_idata%soilisat_pfs   = 0
    clm_pf_idata%soilt_pfs      = 0

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
    clm_pf_idata%decomp_npools_vr_som1_clmp = 0
    clm_pf_idata%decomp_npools_vr_som2_clmp = 0
    clm_pf_idata%decomp_npools_vr_som3_clmp = 0
    clm_pf_idata%decomp_npools_vr_som4_clmp = 0
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
    clm_pf_idata%decomp_npools_vr_som1_pfs = 0
    clm_pf_idata%decomp_npools_vr_som2_pfs = 0
    clm_pf_idata%decomp_npools_vr_som3_pfs = 0
    clm_pf_idata%decomp_npools_vr_som4_pfs = 0
    clm_pf_idata%smin_no3_vr_pfs       = 0
    clm_pf_idata%smin_nh4_vr_pfs       = 0

    !ground/soil C/N rates as source/sink
    clm_pf_idata%rate_lit1c_clmp            = 0
    clm_pf_idata%rate_lit2c_clmp            = 0
    clm_pf_idata%rate_lit3c_clmp            = 0
    clm_pf_idata%rate_cwdc_clmp             = 0
    clm_pf_idata%rate_som1c_clmp            = 0
    clm_pf_idata%rate_som2c_clmp            = 0
    clm_pf_idata%rate_som3c_clmp            = 0
    clm_pf_idata%rate_som4c_clmp            = 0
    clm_pf_idata%rate_lit1n_clmp            = 0
    clm_pf_idata%rate_lit2n_clmp            = 0
    clm_pf_idata%rate_lit3n_clmp            = 0
    clm_pf_idata%rate_cwdn_clmp             = 0
    clm_pf_idata%rate_som1n_clmp            = 0
    clm_pf_idata%rate_som2n_clmp            = 0
    clm_pf_idata%rate_som3n_clmp            = 0
    clm_pf_idata%rate_som4n_clmp            = 0
    clm_pf_idata%rate_smin_no3_clmp         = 0
    clm_pf_idata%rate_smin_nh4_clmp         = 0
    clm_pf_idata%rate_plantndemand_clmp     = 0

    clm_pf_idata%rate_lit1c_pfs            = 0
    clm_pf_idata%rate_lit2c_pfs            = 0
    clm_pf_idata%rate_lit3c_pfs            = 0
    clm_pf_idata%rate_cwdc_pfs             = 0
    clm_pf_idata%rate_som1c_pfs            = 0
    clm_pf_idata%rate_som2c_pfs            = 0
    clm_pf_idata%rate_som3c_pfs            = 0
    clm_pf_idata%rate_som4c_pfs            = 0
    clm_pf_idata%rate_lit1n_pfs            = 0
    clm_pf_idata%rate_lit2n_pfs            = 0
    clm_pf_idata%rate_lit3n_pfs            = 0
    clm_pf_idata%rate_cwdn_pfs             = 0
    clm_pf_idata%rate_som1n_pfs            = 0
    clm_pf_idata%rate_som2n_pfs            = 0
    clm_pf_idata%rate_som3n_pfs            = 0
    clm_pf_idata%rate_som4n_pfs            = 0
    clm_pf_idata%rate_smin_no3_pfs         = 0
    clm_pf_idata%rate_smin_nh4_pfs         = 0
    clm_pf_idata%rate_plantndemand_pfs     = 0

    clm_pf_idata%press_maxponding_clmp = 0
    clm_pf_idata%press_maxponding_pfs  = 0
    clm_pf_idata%press_subsurf_clmp = 0
    clm_pf_idata%press_subsurf_pfs  = 0
    clm_pf_idata%press_subbase_clmp = 0
    clm_pf_idata%press_subbase_pfs  = 0
    clm_pf_idata%topbc_seepage  = PETSC_FALSE
    clm_pf_idata%qflux_subsurf_clmp = 0
    clm_pf_idata%qflux_subsurf_pfs  = 0
    clm_pf_idata%qflux_subbase_clmp = 0
    clm_pf_idata%qflux_subbase_pfs  = 0

    clm_pf_idata%gtemp_subsurf_clmp = 0
    clm_pf_idata%gtemp_subsurf_pfs  = 0
    clm_pf_idata%gflux_subbase_clmp = 0
    clm_pf_idata%gflux_subbase_pfs  = 0
    clm_pf_idata%gtemp_subbase_clmp = 0
    clm_pf_idata%gtemp_subbase_pfs  = 0

    ! for updating bgc/TH states from PF to CLM
    clm_pf_idata%sr_pcwmax_pfp= 0
    clm_pf_idata%pcwmax_pfp   = 0
    clm_pf_idata%porosity_pfp = 0
    clm_pf_idata%sr_pcwmax_clms= 0
    clm_pf_idata%pcwmax_clms   = 0
    clm_pf_idata%porosity_clms = 0

    clm_pf_idata%pressure_reference = 1.01325d5

    clm_pf_idata%press_pfp    = 0
    clm_pf_idata%soilpsi_pfp  = 0
    clm_pf_idata%soillsat_pfp = 0
    clm_pf_idata%soilisat_pfp = 0
    clm_pf_idata%soilt_pfp    = 0
    clm_pf_idata%press_clms     = 0
    clm_pf_idata%soilpsi_clms   = 0
    clm_pf_idata%soillsat_clms  = 0
    clm_pf_idata%soilisat_clms  = 0
    clm_pf_idata%soilt_clms     = 0

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
    clm_pf_idata%decomp_npools_vr_som1_pfp = 0
    clm_pf_idata%decomp_npools_vr_som2_pfp = 0
    clm_pf_idata%decomp_npools_vr_som3_pfp = 0
    clm_pf_idata%decomp_npools_vr_som4_pfp = 0
    clm_pf_idata%smin_no3_vr_pfp       = 0
    clm_pf_idata%smin_nh4_vr_pfp       = 0
    clm_pf_idata%smin_nh4sorb_vr_pfp   = 0
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
    clm_pf_idata%decomp_npools_vr_som1_clms = 0
    clm_pf_idata%decomp_npools_vr_som2_clms = 0
    clm_pf_idata%decomp_npools_vr_som3_clms = 0
    clm_pf_idata%decomp_npools_vr_som4_clms = 0
    clm_pf_idata%smin_no3_vr_clms      = 0
    clm_pf_idata%smin_nh4_vr_clms      = 0
    clm_pf_idata%smin_nh4sorb_vr_clms  = 0

    ! for root N extraction calculation
    clm_pf_idata%accextrn_vr_pfp       = 0
    clm_pf_idata%accextrn_vr_clms      = 0

    ! for soil hr calculation
    clm_pf_idata%gco2_vr_pfp            = 0
    clm_pf_idata%gco2_vr_clms           = 0
    clm_pf_idata%gco2_vr_clmp           = 0
    clm_pf_idata%gco2_vr_pfs            = 0

    ! for N2 gas emission calculation
    clm_pf_idata%gn2_vr_pfp       = 0
    clm_pf_idata%gn2_vr_clms      = 0
    clm_pf_idata%gn2_vr_clmp      = 0
    clm_pf_idata%gn2_vr_pfs       = 0

    ! for N2O gas emission calculation
    clm_pf_idata%gn2o_vr_pfp       = 0
    clm_pf_idata%gn2o_vr_clms      = 0
    clm_pf_idata%gn2o_vr_clmp      = 0
    clm_pf_idata%gn2o_vr_pfs       = 0

    ! tracking variables in C-N cycle
    clm_pf_idata%acchr_vr_pfp       = 0
    clm_pf_idata%acchr_vr_clms      = 0

    clm_pf_idata%accnmin_vr_pfp       = 0
    clm_pf_idata%accnmin_vr_clms      = 0

    clm_pf_idata%accnimm_vr_pfp       = 0
    clm_pf_idata%accnimm_vr_clms      = 0

    clm_pf_idata%accngasmin_vr_pfp       = 0
    clm_pf_idata%accngasmin_vr_clms      = 0

    clm_pf_idata%accngasnitr_vr_pfp       = 0
    clm_pf_idata%accngasnitr_vr_clms      = 0

    clm_pf_idata%accngasdeni_vr_pfp       = 0
    clm_pf_idata%accngasdeni_vr_clms      = 0

    ! water & aq. chemical species boundary flow
    clm_pf_idata%qinfl_subsurf_pfp   = 0
    clm_pf_idata%qinfl_subsurf_clms  = 0
    clm_pf_idata%qsurf_subsurf_pfp   = 0
    clm_pf_idata%qsurf_subsurf_clms  = 0
    clm_pf_idata%qflux_subbase_pfp   = 0
    clm_pf_idata%qflux_subbase_clms  = 0
    clm_pf_idata%f_nh4_subsurf_pfp   = 0
    clm_pf_idata%f_nh4_subsurf_clms  = 0
    clm_pf_idata%f_nh4_subbase_pfp   = 0
    clm_pf_idata%f_nh4_subbase_clms  = 0
    clm_pf_idata%f_no3_subsurf_pfp   = 0
    clm_pf_idata%f_no3_subsurf_clms  = 0
    clm_pf_idata%f_no3_subbase_pfp   = 0
    clm_pf_idata%f_no3_subbase_clms  = 0

  end subroutine CLMPFLOTRANIDataInit

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataCreateVec(mycomm)
  ! 
  ! This routine creates PETSc vectors required for data transfer between
  ! CLM and PFLOTRAN.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
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

    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%watfc_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%bulkdensity_dry_clm,ierr)

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
    !call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%press_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%qflx_pf,ierr)

    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%watfc_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%bulkdensity_dry_pf,ierr)

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
    call VecDuplicate(clm_pf_idata%sat_pf,clm_pf_idata%eff_therm_cond_pf,ierr)

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
    call VecDuplicate(clm_pf_idata%sat_clm,clm_pf_idata%eff_therm_cond_clm,ierr)

    ! 2D Surface PFLOTRAN ---to--- 2D Surface CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%nlclm_2dsub,clm_pf_idata%h2osfc_clm,ierr)
    call VecSet(clm_pf_idata%h2osfc_clm,0.d0,ierr)

    !--------------------------------------------------------------------------------------------------------------------------
    !
    ! The following block of data definition is for THC coupled clm-pflotran
    !
    !NOTES (fmy): From mpi vecs To seq. vecs for passing data IS in one-way only at this momment
    !
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
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_npools_vr_som1_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_npools_vr_som2_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_npools_vr_som3_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%decomp_npools_vr_som4_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%smin_no3_vr_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clmp,clm_pf_idata%smin_nh4_vr_clmp,ierr)

    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%porosity_clmp,ierr)     ! soil physical properties (3D)
    call VecSet(clm_pf_idata%porosity_clmp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clmp,clm_pf_idata%sr_clmp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clmp,clm_pf_idata%lamda_clmp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clmp,clm_pf_idata%alpha_clmp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clmp,clm_pf_idata%pcwmax_clmp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clmp,clm_pf_idata%zsoi_clmp,ierr)

    call VecDuplicate(clm_pf_idata%porosity_clmp,clm_pf_idata%press_clmp,ierr)                        ! TH states (3D)
    call VecDuplicate(clm_pf_idata%porosity_clmp,clm_pf_idata%soilpsi_clmp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clmp,clm_pf_idata%soillsat_clmp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clmp,clm_pf_idata%soilisat_clmp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clmp,clm_pf_idata%soilt_clmp,ierr)

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
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_npools_vr_som1_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_npools_vr_som2_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_npools_vr_som3_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%decomp_npools_vr_som4_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%smin_no3_vr_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfs,clm_pf_idata%smin_nh4_vr_pfs,ierr)

    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%porosity_pfs,ierr)
    call VecSet(clm_pf_idata%porosity_pfs,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfs,clm_pf_idata%sr_pfs,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfs,clm_pf_idata%lamda_pfs,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfs,clm_pf_idata%alpha_pfs,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfs,clm_pf_idata%pcwmax_pfs,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfs,clm_pf_idata%zsoi_pfs,ierr)

    call VecDuplicate(clm_pf_idata%porosity_pfs,clm_pf_idata%press_pfs,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfs,clm_pf_idata%soilpsi_pfs,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfs,clm_pf_idata%soillsat_pfs,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfs,clm_pf_idata%soilisat_pfs,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfs,clm_pf_idata%soilt_pfs,ierr)

    ! (ii) BGC/TH interface source/sink (rate) or BC (pressure/flux): 3D subsurface CLM ---to--- 3D subsurface PFLOTRAN
    ! MPI Vecs for CLM
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%rate_lit1c_clmp,ierr)
    call VecSet(clm_pf_idata%rate_lit1c_clmp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_lit2c_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_lit3c_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_cwdc_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_som1c_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_som2c_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_som3c_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_som4c_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_lit1n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_lit2n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_lit3n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_cwdn_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_som1n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_som2n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_som3n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_som4n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_smin_no3_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_smin_nh4_clmp,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_clmp,clm_pf_idata%rate_plantndemand_clmp,ierr)

    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_2dsub,PETSC_DECIDE,clm_pf_idata%press_subsurf_clmp,ierr)    ! TH top BC (2D)
    call VecSet(clm_pf_idata%press_subsurf_clmp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%press_subsurf_clmp,clm_pf_idata%qflux_subsurf_clmp,ierr)
    call VecDuplicate(clm_pf_idata%press_subsurf_clmp,clm_pf_idata%press_maxponding_clmp,ierr)
    call VecDuplicate(clm_pf_idata%press_subsurf_clmp,clm_pf_idata%gtemp_subsurf_clmp,ierr)
    call VecDuplicate(clm_pf_idata%press_subsurf_clmp,clm_pf_idata%press_ref_clmp,ierr)

    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_bottom,PETSC_DECIDE,clm_pf_idata%press_subbase_clmp,ierr)    ! TH bottom BC (2D)
    call VecSet(clm_pf_idata%press_subbase_clmp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%press_subbase_clmp,clm_pf_idata%qflux_subbase_clmp,ierr)
    call VecDuplicate(clm_pf_idata%press_subbase_clmp,clm_pf_idata%gflux_subbase_clmp,ierr)
    call VecDuplicate(clm_pf_idata%press_subbase_clmp,clm_pf_idata%gtemp_subbase_clmp,ierr)

    ! Seq. Vecs for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%rate_lit1c_pfs,ierr)
    call VecSet(clm_pf_idata%rate_lit1c_pfs,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_lit2c_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_lit3c_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_cwdc_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_som1c_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_som2c_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_som3c_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_som4c_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_lit1n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_lit2n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_lit3n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_cwdn_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_som1n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_som2n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_som3n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_som4n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_smin_no3_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_smin_nh4_pfs,ierr)
    call VecDuplicate(clm_pf_idata%rate_lit1c_pfs,clm_pf_idata%rate_plantndemand_pfs,ierr)

    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_2dsub,clm_pf_idata%press_subsurf_pfs,ierr)   ! H
    call VecSet(clm_pf_idata%press_subsurf_pfs,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%press_subsurf_pfs,clm_pf_idata%qflux_subsurf_pfs,ierr)
    call VecDuplicate(clm_pf_idata%press_subsurf_pfs,clm_pf_idata%press_maxponding_pfs,ierr)
    call VecDuplicate(clm_pf_idata%press_subsurf_pfs,clm_pf_idata%gtemp_subsurf_pfs,ierr)            ! T
    call VecDuplicate(clm_pf_idata%press_subsurf_pfs,clm_pf_idata%press_ref_pfs,ierr)

    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_bottom,clm_pf_idata%press_subbase_pfs,ierr)  ! H
    call VecSet(clm_pf_idata%press_subbase_pfs,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%press_subbase_pfs,clm_pf_idata%qflux_subbase_pfs,ierr)
    call VecDuplicate(clm_pf_idata%press_subbase_pfs,clm_pf_idata%gflux_subbase_pfs,ierr)            ! T
    call VecDuplicate(clm_pf_idata%press_subbase_pfs,clm_pf_idata%gtemp_subbase_pfs,ierr)

    ! (iii) BGC state variables: 3D subsurface PFLOTRAN ---to--- 3D subsurface CLM
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
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_npools_vr_som1_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_npools_vr_som2_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_npools_vr_som3_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%decomp_npools_vr_som4_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%smin_no3_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%smin_nh4_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_pfp,clm_pf_idata%smin_nh4sorb_vr_pfp,ierr)
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
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_som1_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_som2_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_som3_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%decomp_npools_vr_som4_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%smin_no3_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%smin_nh4_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_lit1_clms,clm_pf_idata%smin_nh4sorb_vr_clms,ierr)

    ! (iv) TH parameters: 3D subsurface PFLOTRAN ---to--- 3D subsurface CLM
    ! MPI Vecs for PFLOTRAN
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_sub,PETSC_DECIDE,clm_pf_idata%porosity_pfp,ierr)
    call VecSet(clm_pf_idata%porosity_pfp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfp,clm_pf_idata%sr_pcwmax_pfp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfp,clm_pf_idata%pcwmax_pfp,ierr)

    call VecDuplicate(clm_pf_idata%porosity_pfp,clm_pf_idata%press_pfp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfp,clm_pf_idata%soilpsi_pfp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfp,clm_pf_idata%soillsat_pfp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfp,clm_pf_idata%soilisat_pfp,ierr)
    call VecDuplicate(clm_pf_idata%porosity_pfp,clm_pf_idata%soilt_pfp,ierr)

    ! Seq. Vecs for CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_sub,clm_pf_idata%porosity_clms,ierr)
    call VecSet(clm_pf_idata%porosity_clms,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clms,clm_pf_idata%sr_pcwmax_clms,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clms,clm_pf_idata%pcwmax_clms,ierr)

    call VecDuplicate(clm_pf_idata%porosity_clms,clm_pf_idata%press_clms,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clms,clm_pf_idata%soilpsi_clms,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clms,clm_pf_idata%soillsat_clms,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clms,clm_pf_idata%soilisat_clms,ierr)
    call VecDuplicate(clm_pf_idata%porosity_clms,clm_pf_idata%soilt_clms,ierr)

    ! (v) BGC flux variables: 3D subsurface PFLOTRAN ---to--- 3D subsurface CLM
    ! MPI Vecs for PFLOTRAN
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_sub,PETSC_DECIDE,clm_pf_idata%gco2_vr_pfp,ierr)
    call VecSet(clm_pf_idata%gco2_vr_pfp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_pfp,clm_pf_idata%gn2_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_pfp,clm_pf_idata%gn2o_vr_pfp,ierr)
    !
    call VecDuplicate(clm_pf_idata%gco2_vr_pfp,clm_pf_idata%accextrn_vr_pfp,ierr)
    !
    call VecDuplicate(clm_pf_idata%gco2_vr_pfp,clm_pf_idata%acchr_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_pfp,clm_pf_idata%accnmin_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_pfp,clm_pf_idata%accnimm_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_pfp,clm_pf_idata%accngasmin_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_pfp,clm_pf_idata%accngasnitr_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_pfp,clm_pf_idata%accngasdeni_vr_pfp,ierr)

    ! Seq. Vecs for CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_sub,clm_pf_idata%gco2_vr_clms,ierr)
    call VecSet(clm_pf_idata%gco2_vr_clms,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_clms,clm_pf_idata%gn2_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_clms,clm_pf_idata%gn2o_vr_clms,ierr)
    !
    call VecDuplicate(clm_pf_idata%gco2_vr_clms,clm_pf_idata%accextrn_vr_clms,ierr)
    !
    call VecDuplicate(clm_pf_idata%gco2_vr_clms,clm_pf_idata%acchr_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_clms,clm_pf_idata%accnmin_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_clms,clm_pf_idata%accnimm_vr_clms,ierr)

    call VecDuplicate(clm_pf_idata%gco2_vr_clms,clm_pf_idata%accngasmin_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_clms,clm_pf_idata%accngasnitr_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_clms,clm_pf_idata%accngasdeni_vr_clms,ierr)

    ! MPI Vecs for CLM to pass reset aq. conc back to PF
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%gco2_vr_clmp,ierr)
    call VecSet(clm_pf_idata%gco2_vr_clmp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_clmp,clm_pf_idata%gn2_vr_clmp,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_clmp,clm_pf_idata%gn2o_vr_clmp,ierr)
    ! Seq. Vecs for PFLOTRAN to get reset aq. conc back from CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%gco2_vr_pfs,ierr)
    call VecSet(clm_pf_idata%gco2_vr_pfs,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_pfs,clm_pf_idata%gn2_vr_pfs,ierr)
    call VecDuplicate(clm_pf_idata%gco2_vr_pfs,clm_pf_idata%gn2o_vr_pfs,ierr)

    ! (v) BC flow variables: 2D faces of subsurface PFLOTRAN ---to--- 2D faces of subsurface CLM
    ! MPI Vecs for PFLOTRAN
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_2dsub,PETSC_DECIDE,clm_pf_idata%qinfl_subsurf_pfp,ierr)
    call VecSet(clm_pf_idata%qinfl_subsurf_pfp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%qinfl_subsurf_pfp,clm_pf_idata%qsurf_subsurf_pfp,ierr)
    call VecDuplicate(clm_pf_idata%qinfl_subsurf_pfp,clm_pf_idata%f_nh4_subsurf_pfp,ierr)
    call VecDuplicate(clm_pf_idata%qinfl_subsurf_pfp,clm_pf_idata%f_no3_subsurf_pfp,ierr)

    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_bottom,PETSC_DECIDE,clm_pf_idata%qflux_subbase_pfp,ierr)
    call VecSet(clm_pf_idata%qflux_subbase_pfp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%qflux_subbase_pfp,clm_pf_idata%f_nh4_subbase_pfp,ierr)
    call VecDuplicate(clm_pf_idata%qflux_subbase_pfp,clm_pf_idata%f_no3_subbase_pfp,ierr)

    ! Seq. Vecs for CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_2dsub,clm_pf_idata%qinfl_subsurf_clms,ierr)
    call VecSet(clm_pf_idata%qinfl_subsurf_clms,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%qinfl_subsurf_clms,clm_pf_idata%qsurf_subsurf_clms,ierr)
    call VecDuplicate(clm_pf_idata%qinfl_subsurf_clms,clm_pf_idata%f_nh4_subsurf_clms,ierr)
    call VecDuplicate(clm_pf_idata%qinfl_subsurf_clms,clm_pf_idata%f_no3_subsurf_clms,ierr)

    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_bottom,clm_pf_idata%qflux_subbase_clms,ierr)
    call VecSet(clm_pf_idata%qflux_subbase_clms,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%qflux_subbase_clms,clm_pf_idata%f_nh4_subbase_clms,ierr)
    call VecDuplicate(clm_pf_idata%qflux_subbase_clms,clm_pf_idata%f_no3_subbase_clms,ierr)

    !---------------------------------------------

  end subroutine CLMPFLOTRANIDataCreateVec

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataDestroy()
  ! 
  ! This routine destroys PETSc vectors that were created for data transfer.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  ! 
  
    implicit none
    
    PetscErrorCode :: ierr

    if(clm_pf_idata%hksat_x_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_x_clm,ierr)
    if(clm_pf_idata%hksat_y_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_y_clm,ierr)
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)
    if(clm_pf_idata%sucsat_clm  /= 0) call VecDestroy(clm_pf_idata%sucsat_clm,ierr)
    if(clm_pf_idata%watsat_clm  /= 0) call VecDestroy(clm_pf_idata%watsat_clm,ierr)
    if(clm_pf_idata%bsw_clm  /= 0) call VecDestroy(clm_pf_idata%bsw_clm,ierr)
    if(clm_pf_idata%press_clm  /= 0) call VecDestroy(clm_pf_idata%press_clm,ierr)
    if(clm_pf_idata%qflx_clm  /= 0) call VecDestroy(clm_pf_idata%qflx_clm,ierr)
    if(clm_pf_idata%watfc_clm  /= 0) call VecDestroy(clm_pf_idata%watfc_clm,ierr)
    if(clm_pf_idata%bulkdensity_dry_clm  /= 0) call VecDestroy(clm_pf_idata%bulkdensity_dry_clm,ierr)

    if(clm_pf_idata%hksat_x_pf  /= 0) call VecDestroy(clm_pf_idata%hksat_x_pf,ierr)
    if(clm_pf_idata%hksat_y_pf  /= 0) call VecDestroy(clm_pf_idata%hksat_y_pf,ierr)
    if(clm_pf_idata%hksat_z_pf  /= 0) call VecDestroy(clm_pf_idata%hksat_z_pf,ierr)
    if(clm_pf_idata%sucsat_pf  /= 0) call VecDestroy(clm_pf_idata%sucsat_pf,ierr)
    if(clm_pf_idata%watsat_pf  /= 0) call VecDestroy(clm_pf_idata%watsat_pf,ierr)
    if(clm_pf_idata%bsw_pf  /= 0) call VecDestroy(clm_pf_idata%bsw_pf,ierr)
    !if(clm_pf_idata%press_pf  /= 0) call VecDestroy(clm_pf_idata%press_pf,ierr)
    if(clm_pf_idata%qflx_pf  /= 0) call VecDestroy(clm_pf_idata%qflx_pf,ierr)
    if(clm_pf_idata%watfc_pf  /= 0) call VecDestroy(clm_pf_idata%watfc_pf,ierr)
    if(clm_pf_idata%bulkdensity_dry_pf  /= 0) call VecDestroy(clm_pf_idata%bulkdensity_dry_pf,ierr)
    
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

    if(clm_pf_idata%eff_therm_cond_clm  /= 0) call VecDestroy(clm_pf_idata%eff_therm_cond_clm,ierr)
    if(clm_pf_idata%eff_therm_cond_pf  /= 0) call VecDestroy(clm_pf_idata%eff_therm_cond_pf,ierr)

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
    if(clm_pf_idata%decomp_npools_vr_som1_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som1_clmp,ierr)
    if(clm_pf_idata%decomp_npools_vr_som2_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som2_clmp,ierr)
    if(clm_pf_idata%decomp_npools_vr_som3_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som3_clmp,ierr)
    if(clm_pf_idata%decomp_npools_vr_som4_clmp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som4_clmp,ierr)
    if(clm_pf_idata%smin_no3_vr_clmp /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_clmp,ierr)
    if(clm_pf_idata%smin_nh4_vr_clmp /= 0) &
       call VecDestroy(clm_pf_idata%smin_nh4_vr_clmp,ierr)

    if(clm_pf_idata%sr_clmp /= 0) &
       call VecDestroy(clm_pf_idata%sr_clmp,ierr)
    if(clm_pf_idata%lamda_clmp /= 0) &
       call VecDestroy(clm_pf_idata%lamda_clmp,ierr)
    if(clm_pf_idata%alpha_clmp /= 0) &
       call VecDestroy(clm_pf_idata%alpha_clmp,ierr)
    if(clm_pf_idata%pcwmax_clmp /= 0) &
       call VecDestroy(clm_pf_idata%pcwmax_clmp,ierr)
    if(clm_pf_idata%porosity_clmp /= 0) &
       call VecDestroy(clm_pf_idata%porosity_clmp,ierr)
    if(clm_pf_idata%zsoi_clmp /= 0) &
       call VecDestroy(clm_pf_idata%zsoi_clmp,ierr)

    if(clm_pf_idata%press_clmp /= 0) &
       call VecDestroy(clm_pf_idata%press_clmp,ierr)
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
    if(clm_pf_idata%decomp_npools_vr_som1_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som1_pfs,ierr)
    if(clm_pf_idata%decomp_npools_vr_som2_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som2_pfs,ierr)
    if(clm_pf_idata%decomp_npools_vr_som3_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som3_pfs,ierr)
    if(clm_pf_idata%decomp_npools_vr_som4_pfs /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som4_pfs,ierr)
    if(clm_pf_idata%smin_no3_vr_pfs /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_pfs,ierr)
    if(clm_pf_idata%smin_nh4_vr_pfs /= 0) &
      call VecDestroy(clm_pf_idata%smin_nh4_vr_pfs,ierr)

    if(clm_pf_idata%sr_pfs /= 0) &
       call VecDestroy(clm_pf_idata%sr_pfs,ierr)
    if(clm_pf_idata%lamda_pfs /= 0) &
       call VecDestroy(clm_pf_idata%lamda_pfs,ierr)
    if(clm_pf_idata%alpha_pfs /= 0) &
       call VecDestroy(clm_pf_idata%alpha_pfs,ierr)
    if(clm_pf_idata%pcwmax_pfs /= 0) &
       call VecDestroy(clm_pf_idata%pcwmax_pfs,ierr)
    if(clm_pf_idata%porosity_pfs /= 0) &
       call VecDestroy(clm_pf_idata%porosity_pfs,ierr)
    if(clm_pf_idata%zsoi_pfs /= 0) &
       call VecDestroy(clm_pf_idata%zsoi_pfs,ierr)

    if(clm_pf_idata%press_pfs /= 0) &
      call VecDestroy(clm_pf_idata%press_pfs,ierr)
    if(clm_pf_idata%soilpsi_pfs /= 0) &
      call VecDestroy(clm_pf_idata%soilpsi_pfs,ierr)
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
    if(clm_pf_idata%rate_som1c_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_som1c_clmp,ierr)
    if(clm_pf_idata%rate_som2c_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_som2c_clmp,ierr)
    if(clm_pf_idata%rate_som3c_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_som3c_clmp,ierr)
    if(clm_pf_idata%rate_som4c_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_som4c_clmp,ierr)
    if(clm_pf_idata%rate_lit1n_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1n_clmp,ierr)
    if(clm_pf_idata%rate_lit2n_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2n_clmp,ierr)
    if(clm_pf_idata%rate_lit3n_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3n_clmp,ierr)
    if(clm_pf_idata%rate_cwdn_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdn_clmp,ierr)
    if(clm_pf_idata%rate_som1n_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_som1n_clmp,ierr)
    if(clm_pf_idata%rate_som2n_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_som2n_clmp,ierr)
    if(clm_pf_idata%rate_som3n_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_som3n_clmp,ierr)
    if(clm_pf_idata%rate_som4n_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_som4n_clmp,ierr)
    if(clm_pf_idata%rate_plantndemand_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_plantndemand_clmp,ierr)
    if(clm_pf_idata%rate_smin_no3_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_smin_no3_clmp,ierr)
    if(clm_pf_idata%rate_smin_nh4_clmp /= 0) &
       call VecDestroy(clm_pf_idata%rate_smin_nh4_clmp,ierr)
    if(clm_pf_idata%rate_lit1c_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1c_pfs,ierr)
    if(clm_pf_idata%rate_lit2c_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2c_pfs,ierr)
    if(clm_pf_idata%rate_lit3c_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3c_pfs,ierr)
    if(clm_pf_idata%rate_cwdc_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdc_pfs,ierr)
    if(clm_pf_idata%rate_som1c_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_som1c_pfs,ierr)
    if(clm_pf_idata%rate_som2c_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_som2c_pfs,ierr)
    if(clm_pf_idata%rate_som3c_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_som3c_pfs,ierr)
    if(clm_pf_idata%rate_som4c_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_som4c_pfs,ierr)
    if(clm_pf_idata%rate_lit1n_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit1n_pfs,ierr)
    if(clm_pf_idata%rate_lit2n_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit2n_pfs,ierr)
    if(clm_pf_idata%rate_lit3n_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_lit3n_pfs,ierr)
    if(clm_pf_idata%rate_cwdn_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_cwdn_pfs,ierr)
    if(clm_pf_idata%rate_som1n_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_som1n_pfs,ierr)
    if(clm_pf_idata%rate_som2n_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_som2n_pfs,ierr)
    if(clm_pf_idata%rate_som3n_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_som3n_pfs,ierr)
    if(clm_pf_idata%rate_som4n_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_som4n_pfs,ierr)
    if(clm_pf_idata%rate_plantndemand_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_plantndemand_pfs,ierr)
    if(clm_pf_idata%rate_smin_no3_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_smin_no3_pfs,ierr)
    if(clm_pf_idata%rate_smin_nh4_pfs /= 0) &
       call VecDestroy(clm_pf_idata%rate_smin_nh4_pfs,ierr)

    if(clm_pf_idata%press_maxponding_clmp  /= 0) &
       call VecDestroy(clm_pf_idata%press_maxponding_clmp,ierr)
    if(clm_pf_idata%press_maxponding_pfs  /= 0) &
       call VecDestroy(clm_pf_idata%press_maxponding_pfs,ierr)
    if(clm_pf_idata%press_subsurf_clmp  /= 0) &
       call VecDestroy(clm_pf_idata%press_subsurf_clmp,ierr)
    if(clm_pf_idata%press_subsurf_pfs  /= 0) &
       call VecDestroy(clm_pf_idata%press_subsurf_pfs,ierr)
    if(clm_pf_idata%press_subbase_clmp  /= 0) &
       call VecDestroy(clm_pf_idata%press_subbase_clmp,ierr)
    if(clm_pf_idata%press_subbase_pfs  /= 0) &
       call VecDestroy(clm_pf_idata%press_subbase_pfs,ierr)
    if(clm_pf_idata%qflux_subsurf_clmp  /= 0) &
       call VecDestroy(clm_pf_idata%qflux_subsurf_clmp,ierr)
    if(clm_pf_idata%qflux_subsurf_pfs  /= 0) &
       call VecDestroy(clm_pf_idata%qflux_subsurf_pfs,ierr)
    if(clm_pf_idata%qflux_subbase_clmp  /= 0) &
       call VecDestroy(clm_pf_idata%qflux_subbase_clmp,ierr)
    if(clm_pf_idata%qflux_subbase_pfs  /= 0) &
       call VecDestroy(clm_pf_idata%qflux_subbase_pfs,ierr)

    if(clm_pf_idata%gtemp_subsurf_clmp  /= 0) &
       call VecDestroy(clm_pf_idata%gtemp_subsurf_clmp,ierr)
    if(clm_pf_idata%gtemp_subsurf_pfs  /= 0) &
       call VecDestroy(clm_pf_idata%gtemp_subsurf_pfs,ierr)
    if(clm_pf_idata%gflux_subbase_clmp  /= 0) &
       call VecDestroy(clm_pf_idata%gflux_subbase_clmp,ierr)
    if(clm_pf_idata%gflux_subbase_pfs  /= 0) &
       call VecDestroy(clm_pf_idata%gflux_subbase_pfs,ierr)
    if(clm_pf_idata%gtemp_subbase_clmp  /= 0) &
       call VecDestroy(clm_pf_idata%gtemp_subbase_clmp,ierr)
    if(clm_pf_idata%gtemp_subbase_pfs  /= 0) &
       call VecDestroy(clm_pf_idata%gtemp_subbase_pfs,ierr)

    !
    if(clm_pf_idata%sr_pcwmax_pfp /= 0) &
       call VecDestroy(clm_pf_idata%sr_pcwmax_pfp,ierr)
    if(clm_pf_idata%pcwmax_pfp /= 0) &
       call VecDestroy(clm_pf_idata%pcwmax_pfp,ierr)
    if(clm_pf_idata%porosity_pfp /= 0) &
       call VecDestroy(clm_pf_idata%porosity_pfp,ierr)

    if(clm_pf_idata%sr_pcwmax_clms /= 0) &
       call VecDestroy(clm_pf_idata%sr_pcwmax_clms,ierr)
    if(clm_pf_idata%pcwmax_clms /= 0) &
       call VecDestroy(clm_pf_idata%pcwmax_clms,ierr)
    if(clm_pf_idata%porosity_clms /= 0) &
       call VecDestroy(clm_pf_idata%porosity_clms,ierr)

    ! update TH-BGC states
    if(clm_pf_idata%press_pfp /= 0) &
       call VecDestroy(clm_pf_idata%press_pfp,ierr)
    if(clm_pf_idata%soilpsi_pfp /= 0) &
       call VecDestroy(clm_pf_idata%soilpsi_pfp,ierr)
    if(clm_pf_idata%soillsat_pfp /= 0) &
       call VecDestroy(clm_pf_idata%soillsat_pfp,ierr)
    if(clm_pf_idata%soilisat_pfp /= 0) &
       call VecDestroy(clm_pf_idata%soilisat_pfp,ierr)
    if(clm_pf_idata%soilt_pfp /= 0) &
       call VecDestroy(clm_pf_idata%soilt_pfp,ierr)

    if(clm_pf_idata%press_clms /= 0) &
       call VecDestroy(clm_pf_idata%press_clms,ierr)
    if(clm_pf_idata%soilpsi_clms /= 0) &
       call VecDestroy(clm_pf_idata%soilpsi_clms,ierr)
    if(clm_pf_idata%soillsat_clms /= 0) &
       call VecDestroy(clm_pf_idata%soillsat_clms,ierr)
    if(clm_pf_idata%soilisat_clms /= 0) &
       call VecDestroy(clm_pf_idata%soilisat_clms,ierr)
    if(clm_pf_idata%soilt_clms /= 0) &
       call VecDestroy(clm_pf_idata%soilt_clms,ierr)

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
    if(clm_pf_idata%decomp_npools_vr_som1_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som1_pfp,ierr)
    if(clm_pf_idata%decomp_npools_vr_som2_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som2_pfp,ierr)
    if(clm_pf_idata%decomp_npools_vr_som3_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som3_pfp,ierr)
    if(clm_pf_idata%decomp_npools_vr_som4_pfp /= 0) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_som4_pfp,ierr)
    if(clm_pf_idata%smin_no3_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_pfp,ierr)
    if(clm_pf_idata%smin_nh4_vr_pfp /= 0) &
      call VecDestroy(clm_pf_idata%smin_nh4_vr_pfp,ierr)
    if(clm_pf_idata%smin_nh4sorb_vr_pfp /= 0) &
      call VecDestroy(clm_pf_idata%smin_nh4sorb_vr_pfp,ierr)

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
    if(clm_pf_idata%decomp_npools_vr_som1_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som1_clms,ierr)
    if(clm_pf_idata%decomp_npools_vr_som2_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som2_clms,ierr)
    if(clm_pf_idata%decomp_npools_vr_som3_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som3_clms,ierr)
    if(clm_pf_idata%decomp_npools_vr_som4_clms /= 0) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_som4_clms,ierr)
    if(clm_pf_idata%smin_no3_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_clms,ierr)
    if(clm_pf_idata%smin_nh4_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%smin_nh4_vr_clms,ierr)
    if(clm_pf_idata%smin_nh4sorb_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%smin_nh4sorb_vr_clms,ierr)

    ! update BGC fluxes
    if(clm_pf_idata%accextrn_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%accextrn_vr_pfp,ierr)
    if(clm_pf_idata%accextrn_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%accextrn_vr_clms,ierr)

    if(clm_pf_idata%gco2_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%gco2_vr_pfp,ierr)
    if(clm_pf_idata%gco2_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%gco2_vr_clms,ierr)
    if(clm_pf_idata%gco2_vr_clmp /= 0) &
       call VecDestroy(clm_pf_idata%gco2_vr_clmp,ierr)
    if(clm_pf_idata%gco2_vr_pfs /= 0) &
       call VecDestroy(clm_pf_idata%gco2_vr_pfs,ierr)

    if(clm_pf_idata%gn2_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%gn2_vr_pfp,ierr)
    if(clm_pf_idata%gn2_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%gn2_vr_clms,ierr)
    if(clm_pf_idata%gn2_vr_clmp /= 0) &
       call VecDestroy(clm_pf_idata%gn2_vr_clmp,ierr)
    if(clm_pf_idata%gn2_vr_pfs /= 0) &
       call VecDestroy(clm_pf_idata%gn2_vr_pfs,ierr)

    if(clm_pf_idata%gn2o_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%gn2o_vr_pfp,ierr)
    if(clm_pf_idata%gn2o_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%gn2o_vr_clms,ierr)
    if(clm_pf_idata%gn2o_vr_clmp /= 0) &
       call VecDestroy(clm_pf_idata%gn2o_vr_clmp,ierr)
    if(clm_pf_idata%gn2o_vr_pfs /= 0) &
       call VecDestroy(clm_pf_idata%gn2o_vr_pfs,ierr)

    if(clm_pf_idata%acchr_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%acchr_vr_pfp,ierr)
    if(clm_pf_idata%acchr_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%acchr_vr_clms,ierr)

    if(clm_pf_idata%accnmin_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%accnmin_vr_pfp,ierr)
    if(clm_pf_idata%accnmin_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%accnmin_vr_clms,ierr)

    if(clm_pf_idata%accnimm_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%accnimm_vr_pfp,ierr)
    if(clm_pf_idata%accnimm_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%accnimm_vr_clms,ierr)

    if(clm_pf_idata%accngasmin_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%accngasmin_vr_pfp,ierr)
    if(clm_pf_idata%accngasmin_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%accngasmin_vr_clms,ierr)

    if(clm_pf_idata%accngasnitr_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%accngasnitr_vr_pfp,ierr)
    if(clm_pf_idata%accngasnitr_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%accngasnitr_vr_clms,ierr)

    if(clm_pf_idata%accngasdeni_vr_pfp /= 0) &
       call VecDestroy(clm_pf_idata%accngasdeni_vr_pfp,ierr)
    if(clm_pf_idata%accngasdeni_vr_clms /= 0) &
       call VecDestroy(clm_pf_idata%accngasdeni_vr_clms,ierr)

    !
    if(clm_pf_idata%qinfl_subsurf_pfp /= 0) &
       call VecDestroy(clm_pf_idata%qinfl_subsurf_pfp,ierr)
    if(clm_pf_idata%qinfl_subsurf_clms /= 0) &
       call VecDestroy(clm_pf_idata%qinfl_subsurf_clms,ierr)
    if(clm_pf_idata%qsurf_subsurf_pfp /= 0) &
       call VecDestroy(clm_pf_idata%qsurf_subsurf_pfp,ierr)
    if(clm_pf_idata%qsurf_subsurf_clms /= 0) &
       call VecDestroy(clm_pf_idata%qsurf_subsurf_clms,ierr)
    if(clm_pf_idata%qflux_subbase_pfp /= 0) &
       call VecDestroy(clm_pf_idata%qflux_subbase_pfp,ierr)
    if(clm_pf_idata%qflux_subbase_clms /= 0) &
       call VecDestroy(clm_pf_idata%qflux_subbase_clms,ierr)

    if(clm_pf_idata%f_nh4_subsurf_pfp /= 0) &
       call VecDestroy(clm_pf_idata%f_nh4_subsurf_pfp,ierr)
    if(clm_pf_idata%f_nh4_subsurf_clms /= 0) &
       call VecDestroy(clm_pf_idata%f_nh4_subsurf_clms,ierr)
    if(clm_pf_idata%f_nh4_subbase_pfp /= 0) &
       call VecDestroy(clm_pf_idata%f_nh4_subbase_pfp,ierr)
    if(clm_pf_idata%f_nh4_subbase_clms /= 0) &
       call VecDestroy(clm_pf_idata%f_nh4_subbase_clms,ierr)

    if(clm_pf_idata%f_no3_subsurf_pfp /= 0) &
       call VecDestroy(clm_pf_idata%f_no3_subsurf_pfp,ierr)
    if(clm_pf_idata%f_no3_subsurf_clms /= 0) &
       call VecDestroy(clm_pf_idata%f_no3_subsurf_clms,ierr)
    if(clm_pf_idata%f_no3_subbase_pfp /= 0) &
       call VecDestroy(clm_pf_idata%f_no3_subbase_pfp,ierr)
    if(clm_pf_idata%f_no3_subbase_clms /= 0) &
       call VecDestroy(clm_pf_idata%f_no3_subbase_clms,ierr)

    !----------------------------------------------------------------------------------

  end subroutine CLMPFLOTRANIDataDestroy

end module clm_pflotran_interface_data
