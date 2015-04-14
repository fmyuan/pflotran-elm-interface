module clm_pflotran_interface_data
!
! NOTES for convenience:
!        (1) '_pfp': mpi vecs for PF variables; '_clmp': mpi vecs for CLM variables;
!            '_pfs': seq vecs for PF variables; '_clms': seq. vecs for CLM variables;
!        (2) '_sub_': 3D (XYZ) subsurface domain's variables;
!                     if '_x/y/z_' used for different directions of 3D domain, '_sub_' will NOT use.
!        (3) '_subsurf_': 2D (XY) face of 3D domain's varialbes;
!                         for bottom face, uses '_subbase_', but essentially supposing same shape/area as top face
!            '_srf_': for variables with 2D-grids at ground surface, which may or may not same as '_subsurf_'.

  implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  private

  type, public :: clm_pflotran_idata_type

  ! Time invariant data:

  ! (i) domains (grids/cells) and Mesh property

  ! Number of grids for the 2D surface domain
  !(maybe same as as 'subsurf' below OR not, So that having some flexibility useful for future)
  PetscInt :: nlclm_srf  ! num of local clm grids
  PetscInt :: ngclm_srf  ! num of ghosted clm grids (ghosted = local+ghosts)
  PetscInt :: nlpf_srf   ! num of local pflotran grids
  PetscInt :: ngpf_srf   ! num of ghosted pflotran grids (ghosted = local+ghosts)

  ! Area of grids of 2-D surface domain
  Vec :: area_srf_clmp   ! mpi vec
  Vec :: area_srf_pfs    ! seq vec
  Vec :: area_srf_pfp    ! mpi vec
  Vec :: area_srf_clms   ! seq vec

  ! Numbers of soil layers mapped btw CLM and PFLOTRAN for 3-D sub cells
  PetscInt :: nlayer_mapped

  ! Number of cells for the 3D subsurface domain
  PetscInt :: nlclm_subcell ! num of local clm cells
  PetscInt :: ngclm_subcell ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_subcell  ! num of local pflotran cells
  PetscInt :: ngpf_subcell  ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the 2D top/bot face of the 3D subsurface domain
  PetscInt :: nlclm_subsurf  ! num of local clm cells
  PetscInt :: ngclm_subsurf  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_subsurf   ! num of local pflotran cells
  PetscInt :: ngpf_subsurf   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Area of top/bottom (two sets not-yet considered) faces for the 3D cells
  Vec :: area_subsurf_clmp   ! mpi vec
  Vec :: area_subsurf_pfs    ! seq vec
  Vec :: area_subsurf_pfp    ! mpi vec
  Vec :: area_subsurf_clms   ! seq vec

  ! 3-D cell centroid Z depth (m) from ground-zero (+ downward)
  Vec :: zsoi_sub_clmp
  Vec :: zsoi_sub_pfs
  Vec :: zsoi_sub_pfp
  Vec :: zsoi_sub_clms

  ! in case the near-surface air pressure (grided) from CLM to PF
  Vec :: press_ref_clmp
  Vec :: press_ref_pfs
  Vec :: press_ref_pfp
  Vec :: press_ref_clms
  PetscReal :: pressure_reference     ! currently a constant
  
  ! (ii) Soil properties -

  ! for CLM --> PF
  Vec :: hksat_x_clmp                 ! common properties
  Vec :: hksat_y_clmp
  Vec :: hksat_z_clmp
  Vec :: bulkdensity_dry_sub_clmp
  Vec :: effporosity_sub_clmp         ! adjustable (effective poro), so NOT exactly from bulkdensity_dry
  Vec :: watsat_sub_clmp
  Vec :: watfc_sub_clmp
  Vec :: watpcwmax_sub_clmp
  Vec :: sucsat_sub_clmp              ! clapp-Horberger hydraulic properties
  Vec :: bsw_sub_clmp
  Vec :: sr_sub_clmp                  ! MVM soil hydraulic properties
  Vec :: lamda_sub_clmp
  Vec :: alpha_sub_clmp
  Vec :: pcwmax_sub_clmp
  Vec :: tcond_dry_sub_clmp           ! dry thermal conductivity (W/m/K)
  Vec :: tcond_wet_sub_clmp           ! wet thermal conductivity (W/m/K)
  Vec :: hc_dry_sub_clmp              ! bulk dry (soil solid) heat capacity (J/m3/K)

  Vec :: hksat_x_pfs                  ! common properties
  Vec :: hksat_y_pfs
  Vec :: hksat_z_pfs
  Vec :: bulkdensity_dry_sub_pfs
  Vec :: effporosity_sub_pfs
  Vec :: watsat_sub_pfs
  Vec :: watfc_sub_pfs
  Vec :: watpcwmax_sub_pfs
  Vec :: sucsat_sub_pfs              ! clapp-Horberger hydraulic properties
  Vec :: bsw_sub_pfs
  Vec :: sr_sub_pfs                  ! MVM soil hydraulic properties
  Vec :: lamda_sub_pfs
  Vec :: alpha_sub_pfs
  Vec :: pcwmax_sub_pfs
  Vec :: tcond_dry_sub_pfs           ! dry (bulk) thermal conductivity
  Vec :: tcond_wet_sub_pfs           ! wet (bulk) thermal conductivity
  Vec :: hc_dry_sub_pfs              ! dry (bulk) heat capacity

  ! for PF --> CLM
  Vec :: hksat_x_pfp                 ! common properties
  Vec :: hksat_y_pfp
  Vec :: hksat_z_pfp
  Vec :: bulkdensity_dry_sub_pfp
  Vec :: effporosity_sub_pfp
  Vec :: watsat_sub_pfp
  Vec :: watfc_sub_pfp
  Vec :: watpcwmax_sub_pfp
  Vec :: sucsat_sub_pfp              ! clapp-Horberger hydraulic properties
  Vec :: bsw_sub_pfp
  Vec :: sr_sub_pfp                  ! MVM soil hydraulic properties
  Vec :: lamda_sub_pfp
  Vec :: alpha_sub_pfp
  Vec :: pcwmax_sub_pfp
  Vec :: tcond_dry_sub_pfp           ! dry (bulk) thermal conductivity
  Vec :: tcond_wet_sub_pfp           ! dry (bulk) thermal conductivity
  Vec :: hc_dry_sub_pfp              ! dry (bulk) heat capacity

  Vec :: hksat_x_clms             ! common properties
  Vec :: hksat_y_clms
  Vec :: hksat_z_clms
  Vec :: bulkdensity_dry_sub_clms
  Vec :: effporosity_sub_clms
  Vec :: watsat_sub_clms
  Vec :: watfc_sub_clms
  Vec :: watpcwmax_sub_clms
  Vec :: sucsat_sub_clms              ! clapp-Horberger hydraulic properties
  Vec :: bsw_sub_clms
  Vec :: sr_sub_clms                  ! MVM soil hydraulic properties
  Vec :: lamda_sub_clms
  Vec :: alpha_sub_clms
  Vec :: pcwmax_sub_clms
  Vec :: tcond_dry_sub_clms           ! dry (bulk) thermal conductivity
  Vec :: tcond_wet_sub_clms           ! dry (bulk) thermal conductivity
  Vec :: hc_dry_sub_clms              ! dry (bulk) heat capacity

  ! Time variant data
  
  ! (i) Sink/Source of water/heat for PFLOTRAN's 3D subsurface domain
  Vec :: qflux_sub_clmp   ! mpi vec
  Vec :: qflux_sub_pfs    ! seq vec
  Vec :: qflux_sub_pfp    ! mpi vec
  Vec :: qflux_sub_clms   ! seq vec

  Vec :: gflux_sub_clmp   ! mpi vec
  Vec :: gflux_sub_pfs    ! seq vec
  Vec :: gflux_sub_pfp    ! mpi vec
  Vec :: gflux_sub_clms   ! seq vec
  
  ! (ii) surface water w/ heat I/O for PFLOTRAN's 2D surface domain
  Vec :: h2o_srf_clmp        ! mpi vec (water head: mH2O)
  Vec :: h2o_srf_pfs         ! seq vec
  Vec :: qh2o_srf_clmp       ! mpi vec (water flux: m/sec)
  Vec :: qh2o_srf_pfs        ! seq vec
  Vec :: h2otemp_srf_clmp    ! seq vec (temperature of i/o surface water flux: degC)
  Vec :: h2otemp_srf_pfs     ! mpi vec
  Vec :: gh2o_srf_clmp       ! mpi vec (heat flux with surface i/o water flux: W/m???)
  Vec :: gh2o_srf_pfs        ! seq vec

  Vec :: h2o_srf_pfp         ! mpi vec (water head)
  Vec :: h2o_srf_clms        ! seq vec
  Vec :: qh2o_srf_pfp        ! mpi vec (water flux)
  Vec :: qh2o_srf_clms       ! seq vec
  Vec :: h2otemp_srf_pfp     ! mpi vec
  Vec :: h2otemp_srf_clms    ! seq vec
  Vec :: gh2o_srf_pfp        ! mpi vec
  Vec :: gh2o_srf_clms       ! seq vec

  ! BC: water pressure (Pa) on the top/bottom of 3-D subsurface domain
  !     as boundary conditions from CLM to PF
  Vec :: press_subsurf_clmp    ! mpi vec
  Vec :: press_subbase_clmp    ! mpi vec
  Vec :: press_maxponding_clmp ! mpi vec
  Vec :: press_subsurf_pfs     ! seq vec
  Vec :: press_subbase_pfs     ! seq vec
  Vec :: press_maxponding_pfs  ! seq vec
  ! OR, BC: water infiltration/recharge(drainage) (mH2O/sec) on the top/bottom
  !         of 3-D subsurface domain as boundary conditions from CLM to PF
  Vec :: qflux_subsurf_clmp    ! mpi vec
  Vec :: qflux_subbase_clmp    ! mpi vec
  Vec :: qflux_subsurf_pfs     ! seq vec
  Vec :: qflux_subbase_pfs     ! seq vec

  ! actual mass water flow rate (kgH2O/sec) through the top/bottom BC of 3-D subsurface domain
  ! from PF to CLM
  ! (+ in, - out)
  Vec :: qinfl_subsurf_pfp    ! mpi vec: actual infiltration (+)
  Vec :: qinfl_subsurf_clms   ! seq vec
  Vec :: qsurf_subsurf_pfp    ! mpi vec: actual overland flow - potential-actual infiltration or water upwarding (-)
  Vec :: qsurf_subsurf_clms   ! seq vec
  Vec :: qflux_subbase_pfp    ! mpi vec
  Vec :: qflux_subbase_clms   ! seq vec

  ! (iii) Ground/base heat flux BC for PFLOTRAN's subsurface domain
  !       This BC is applied on top/bottom of the 3-D subsurface domain
  Vec :: gflux_subsurf_clmp  ! mpi vec
  Vec :: gflux_subsurf_pfs   ! seq vec
  Vec :: gflux_subbase_clmp  ! mpi vec
  Vec :: gflux_subbase_pfs   ! seq vec

  Vec :: gflux_subsurf_pfp   ! mpi vec
  Vec :: gflux_subsurf_clms  ! seq vec
  Vec :: gflux_subbase_pfp  ! mpi vec
  Vec :: gflux_subbase_clms   ! seq vec

  ! OR, if ground/bottom temperature at the subsurface interface is known
  Vec :: gtemp_subsurf_clmp  ! mpi vec
  Vec :: gtemp_subsurf_pfs   ! seq vec
  Vec :: gtemp_subbase_clmp  ! mpi vec
  Vec :: gtemp_subbase_pfs   ! seq vec

  Vec :: gtemp_subsurf_pfp   ! mpi vec
  Vec :: gtemp_subsurf_clms  ! seq vec
  Vec :: gtemp_subbase_pfp   ! mpi vec
  Vec :: gtemp_subbase_clms  ! seq vec

  ! (iv) TH state vecs from CLM (mpi) to PF (seq) - 3D cells
  Vec :: press_sub_clmp                     ! water pressure head (Pa)
  Vec :: soilpsi_sub_clmp                   ! soil matric potential (Pa)
  Vec :: soillsat_sub_clmp                  ! soil liq. water saturation (0 - 1)
  Vec :: soilisat_sub_clmp                  ! soil ice water saturation (0 - 1)
  Vec :: soilt_sub_clmp                     ! soil temperature (degC)
  Vec :: press_sub_pfs
  Vec :: soilpsi_sub_pfs
  Vec :: soillsat_sub_pfs
  Vec :: soilisat_sub_pfs
  Vec :: soilt_sub_pfs

  ! TH state vecs from PF (mpi) to CLM (seq) - 3D cells
  Vec :: press_sub_pfp                     ! water pressure head (Pa)
  Vec :: soilpsi_sub_pfp                   ! soil matric potential (Pa) (negative)
  Vec :: soillsat_sub_pfp                  ! soil liq. water saturation (0-1)
  Vec :: soilisat_sub_pfp                  ! soil ice water saturation (0-1)
  Vec :: soilt_sub_pfp                     ! soil temperature (degC)
  Vec :: press_sub_clms
  Vec :: soilpsi_sub_clms
  Vec :: soillsat_sub_clms
  Vec :: soilisat_sub_clms
  Vec :: soilt_sub_clms

  !-------------------------------------------------------------------------------------
  ! note: mapping ONLY can do from MPI vecs to Seq. vecs now
  ! -----BGC vecs from CLM to PF --------------------
  PetscInt :: ndecomp_pools                 ! num of decomposition pools

  ! ground/soil C/N pools from CLM (mpi) to PF (seq)
  Vec :: decomp_cpools_sub_clmp          ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_sub_clmp          ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) n pools
  Vec :: smin_no3_sub_clmp               ! (gN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_sub_clmp               ! (gN/m3) vertically-resolved soil mineral NH4
  Vec :: smin_nh4sorb_sub_clmp           ! (gN/m3) vertically-resolved soil mineral NH4 absorbed, if any
  Vec :: decomp_cpools_sub_pfs           ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_sub_pfs           ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) n pools
  Vec :: smin_no3_sub_pfs                ! (gN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_sub_pfs                ! (gN/m3) vertically-resolved soil mineral NH4
  Vec :: smin_nh4sorb_sub_pfs            ! (gN/m3) vertically-resolved soil mineral NH4 absorbed, if any

  ! time-varying ground/soil C/N rates from CLM (mpi) to PF (seq) (Units: moleC(N)/m3/s)
  ! from CLM to PF
  Vec :: rate_decomp_cpools_sub_clmp
  Vec :: rate_decomp_npools_sub_clmp
  Vec :: rate_smin_no3_sub_clmp
  Vec :: rate_smin_nh4_sub_clmp
  Vec :: rate_plantndemand_sub_clmp
  Vec :: rate_decomp_cpools_sub_pfs
  Vec :: rate_decomp_npools_sub_pfs
  Vec :: rate_smin_no3_sub_pfs
  Vec :: rate_smin_nh4_sub_pfs
  Vec :: rate_plantndemand_sub_pfs

  ! actual aqeuous N mass flow rate(gN/m2/s) at the top (runoff)/bottom (leaching) of 3-D subsurface domain
  ! from PF to CLM
  ! (+ in, - out)
  Vec :: f_nh4_subsurf_pfp    ! mpi vec
  Vec :: f_nh4_subsurf_clms   ! seq vec
  Vec :: f_nh4_subbase_pfp    ! mpi vec
  Vec :: f_nh4_subbase_clms   ! seq vec
  Vec :: f_no3_subsurf_pfp    ! mpi vec
  Vec :: f_no3_subsurf_clms   ! seq vec
  Vec :: f_no3_subbase_pfp    ! mpi vec
  Vec :: f_no3_subbase_clms   ! seq vec

  ! -----BGC vecs from PF (mpi, ghosted) to CLM (seq, local) --------------------

  ! BGC state variables
  Vec :: decomp_cpools_sub_pfp          ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_sub_pfp          ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) n pools
  Vec :: smin_no3_sub_pfp
  Vec :: smin_nh4_sub_pfp                   ! (gN/m3) vertically-resolved total soil mineral NH4 (incl. absorbed)
  Vec :: smin_nh4sorb_sub_pfp               ! (gN/m3) vertically-resolved absorbed NH4-N
  !
  Vec :: decomp_cpools_sub_clms         ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_sub_clms         ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) n pools
  Vec :: smin_no3_sub_clms                  ! (gN/m3) vertically-resolved total soil mineral NO3
  Vec :: smin_nh4_sub_clms                  ! (gN/m3) vertically-resolved total soil mineral NH4 (incl. absorbed)
  Vec :: smin_nh4sorb_sub_clms              ! (gN/m3) vertically-resolved absorbed NH4-N
  !
  ! gases in water (aqueous solution of gases)
  ! gases species is accumulative in 'PFLOTRAN', so needs to calculate their fluxes in the CLM-PF interface and reset back to PFLOTRAN
  Vec :: gco2_sub_pfp                   ! (gC/m3) vertically-resolved soil CO2 C
  Vec :: gco2_sub_clms                  ! (gC/m3) vertically-resolved soil CO2 C
  Vec :: gco2_sub_clmp                  ! (gC/m3) vertically-resolved soil CO2 C, after gas emission
  Vec :: gco2_sub_pfs                   ! (gC/m3) vertically-resolved soil CO2 C, after gas emission

  Vec :: gn2_sub_pfp                    ! (gN/m3) vertically-resolved N2-N
  Vec :: gn2_sub_clms                   ! (gN/m3) vertically-resolved N2-N
  Vec :: gn2_sub_clmp                   ! (gN/m3) vertically-resolved N2-N, after gas emission
  Vec :: gn2_sub_pfs                    ! (gN/m3) vertically-resolved N2-N, after gas emission

  Vec :: gn2o_sub_pfp                   ! (gN/m3) vertically-resolved N2O-N
  Vec :: gn2o_sub_clms                  ! (gN/m3) vertically-resolved N2O-N
  Vec :: gn2o_sub_clmp                  ! (gN/m3) vertically-resolved N2O-N, after gas emission
  Vec :: gn2o_sub_pfs                   ! (gN/m3) vertically-resolved N2O-N, after gas emission

  ! 'accextrn' is accumulative N extract by plant roots in 'PFLOTRAN' within a CLM timestep
  Vec :: accextrn_sub_pfp                ! (gN/m3) vertically-resolved root extraction N
  Vec :: accextrn_sub_clms               ! (gN/m3) vertically-resolved root extraction N

  ! some tracking variables from PFLOTRAN bgc to obtain reaction flux rates which needed by CLM
  Vec :: acchr_sub_pfp                 ! (gC/m3/timestep) vertically-resolved heterotrophic resp. C from decomposition
  Vec :: acchr_sub_clms                ! (gC/m3/timestep) vertically-resolved heterotrophic resp. C from decomposition

  Vec :: accnmin_sub_pfp                ! (gN/m3/timestep) vertically-resolved N mineralization
  Vec :: accnmin_sub_clms               ! (gN/m3/timestep) vertically-resolved N mineralization

  Vec :: accnimm_sub_pfp                ! (gN/m3/timestep) vertically-resolved N immoblization
  Vec :: accnimm_sub_clms               ! (gN/m3/timestep) vertically-resolved N immoblization

  Vec :: accngasmin_sub_pfp              ! (gN/m3/timestep) vertically-resolved N2O-N from mineralization
  Vec :: accngasmin_sub_clms             ! (gN/m3/timestep) vertically-resolved N2O-N from mineralization

  Vec :: accngasnitr_sub_pfp             ! (gN/m3/timestep) vertically-resolved N2O-N from nitrification
  Vec :: accngasnitr_sub_clms            ! (gN/m3/timestep) vertically-resolved N2O-N from nitrification

  Vec :: accngasdeni_sub_pfp             ! (gN/m3/timestep) vertically-resolved N2O-N from denitrification
  Vec :: accngasdeni_sub_clms            ! (gN/m3/timestep) vertically-resolved N2O-N from denitrification

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
    implicit none

  !--------------------------------------------------------------------
  ! (i) domains (grids/cells) and Mesh property
  ! Number of grids for the 2D surface domain
    clm_pf_idata%nlclm_srf = 0
    clm_pf_idata%ngclm_srf = 0
    clm_pf_idata%nlpf_srf  = 0
    clm_pf_idata%ngpf_srf  = 0

  ! Area of surface grids
    clm_pf_idata%area_srf_clmp = 0
    clm_pf_idata%area_srf_pfs  = 0
    clm_pf_idata%area_srf_pfp  = 0
    clm_pf_idata%area_srf_clms = 0

  ! Number of soil layers
    clm_pf_idata%nlayer_mapped = 0

  ! Number of cells for the 3D subsurface domain
    clm_pf_idata%nlclm_subcell = 0
    clm_pf_idata%ngclm_subcell = 0
    clm_pf_idata%nlpf_subcell  = 0
    clm_pf_idata%ngpf_subcell  = 0

  ! Number of cells for the 2D top/bot face of the 3D subsurface domain
    clm_pf_idata%nlclm_subsurf = 0
    clm_pf_idata%ngclm_subsurf = 0
    clm_pf_idata%nlpf_subsurf  = 0
    clm_pf_idata%ngpf_subsurf  = 0

    clm_pf_idata%zsoi_sub_clmp    = 0
    clm_pf_idata%zsoi_sub_pfs     = 0
    clm_pf_idata%zsoi_sub_pfp     = 0
    clm_pf_idata%zsoi_sub_clms    = 0

  ! Area of top/bottom face
    clm_pf_idata%area_subsurf_clmp = 0
    clm_pf_idata%area_subsurf_pfs  = 0
    clm_pf_idata%area_subsurf_pfp  = 0
    clm_pf_idata%area_subsurf_clms = 0

   !--------------------------------------------------------------------
   ! Soil properties -
    clm_pf_idata%hksat_x_clmp = 0
    clm_pf_idata%hksat_y_clmp = 0
    clm_pf_idata%hksat_z_clmp = 0
    clm_pf_idata%bulkdensity_dry_sub_clmp = 0
    clm_pf_idata%effporosity_sub_clmp  = 0
    clm_pf_idata%press_ref_clmp = 0
    clm_pf_idata%watsat_sub_clmp    = 0
    clm_pf_idata%watfc_sub_clmp     = 0
    clm_pf_idata%watpcwmax_sub_clmp = 0
    clm_pf_idata%sucsat_sub_clmp = 0
    clm_pf_idata%bsw_sub_clmp    = 0
    clm_pf_idata%sr_sub_clmp     = 0
    clm_pf_idata%lamda_sub_clmp  = 0
    clm_pf_idata%alpha_sub_clmp  = 0
    clm_pf_idata%pcwmax_sub_clmp = 0
    clm_pf_idata%tcond_dry_sub_clmp  = 0
    clm_pf_idata%tcond_wet_sub_clmp  = 0
    clm_pf_idata%hc_dry_sub_clmp     = 0

    clm_pf_idata%hksat_x_pfs  = 0
    clm_pf_idata%hksat_y_pfs  = 0
    clm_pf_idata%hksat_z_pfs  = 0
    clm_pf_idata%bulkdensity_dry_sub_pfs = 0
    clm_pf_idata%effporosity_sub_pfs = 0
    clm_pf_idata%press_ref_pfs= 0
    clm_pf_idata%watsat_sub_pfs   = 0
    clm_pf_idata%watfc_sub_pfs    = 0
    clm_pf_idata%watpcwmax_sub_pfs= 0
    clm_pf_idata%sucsat_sub_pfs = 0
    clm_pf_idata%bsw_sub_pfs    = 0
    clm_pf_idata%sr_sub_pfs    = 0
    clm_pf_idata%lamda_sub_pfs = 0
    clm_pf_idata%alpha_sub_pfs = 0
    clm_pf_idata%pcwmax_sub_pfs= 0
    clm_pf_idata%tcond_dry_sub_pfs = 0
    clm_pf_idata%tcond_wet_sub_pfs  = 0
    clm_pf_idata%hc_dry_sub_pfs    = 0

    clm_pf_idata%hksat_x_pfp = 0
    clm_pf_idata%hksat_y_pfp = 0
    clm_pf_idata%hksat_z_pfp = 0
    clm_pf_idata%bulkdensity_dry_sub_pfp = 0
    clm_pf_idata%effporosity_sub_pfp = 0
    clm_pf_idata%press_ref_pfp = 0
    clm_pf_idata%watsat_sub_pfp = 0
    clm_pf_idata%watfc_sub_pfp = 0
    clm_pf_idata%watpcwmax_sub_pfp = 0
    clm_pf_idata%sucsat_sub_pfp  = 0
    clm_pf_idata%bsw_sub_pfp = 0
    clm_pf_idata%sr_sub_pfp = 0
    clm_pf_idata%lamda_sub_pfp = 0
    clm_pf_idata%alpha_sub_pfp = 0
    clm_pf_idata%pcwmax_sub_pfp = 0
    clm_pf_idata%tcond_dry_sub_pfp  = 0
    clm_pf_idata%tcond_wet_sub_pfp  = 0
    clm_pf_idata%hc_dry_sub_pfp = 0

    clm_pf_idata%hksat_x_clms = 0
    clm_pf_idata%hksat_y_clms = 0
    clm_pf_idata%hksat_z_clms = 0
    clm_pf_idata%bulkdensity_dry_sub_clms = 0
    clm_pf_idata%effporosity_sub_clms = 0
    clm_pf_idata%press_ref_clms = 0
    clm_pf_idata%watsat_sub_clms = 0
    clm_pf_idata%watfc_sub_clms = 0
    clm_pf_idata%watpcwmax_sub_clms = 0
    clm_pf_idata%sucsat_sub_clms  = 0
    clm_pf_idata%bsw_sub_clms = 0
    clm_pf_idata%sr_sub_clms = 0
    clm_pf_idata%lamda_sub_clms = 0
    clm_pf_idata%alpha_sub_clms = 0
    clm_pf_idata%pcwmax_sub_clms = 0
    clm_pf_idata%tcond_dry_sub_clms = 0
    clm_pf_idata%tcond_dry_sub_clms  = 0
    clm_pf_idata%hc_dry_sub_clms = 0

  ! Time variant data

  ! (i) Sink/Source of water/heat for PFLOTRAN's 3D subsurface domain
    clm_pf_idata%qflux_sub_clmp = 0
    clm_pf_idata%qflux_sub_pfs  = 0
    clm_pf_idata%qflux_sub_pfp  = 0
    clm_pf_idata%qflux_sub_clms = 0

    clm_pf_idata%gflux_sub_clmp = 0
    clm_pf_idata%gflux_sub_pfs  = 0
    clm_pf_idata%gflux_sub_pfp  = 0
    clm_pf_idata%gflux_sub_clms = 0

  ! (ii) surface water I/O for PFLOTRAN's 2D surface domain
    clm_pf_idata%h2o_srf_clmp = 0
    clm_pf_idata%h2o_srf_pfs = 0
    clm_pf_idata%qh2o_srf_clmp = 0
    clm_pf_idata%qh2o_srf_pfs = 0
    clm_pf_idata%h2otemp_srf_clmp = 0
    clm_pf_idata%h2otemp_srf_pfs = 0
    clm_pf_idata%gh2o_srf_clmp = 0
    clm_pf_idata%gh2o_srf_pfs = 0

    clm_pf_idata%h2o_srf_pfp = 0
    clm_pf_idata%h2o_srf_clms = 0
    clm_pf_idata%qh2o_srf_pfp = 0
    clm_pf_idata%qh2o_srf_clms = 0
    clm_pf_idata%h2otemp_srf_pfp = 0
    clm_pf_idata%h2otemp_srf_clms = 0
    clm_pf_idata%gh2o_srf_pfp = 0
    clm_pf_idata%gh2o_srf_clms = 0

  ! BC: water pressure (Pa) on the top/bottom of 3-D subsurface domain
  !     as boundary conditions from CLM to PF
    clm_pf_idata%press_subsurf_clmp = 0
    clm_pf_idata%press_subbase_clmp  = 0
    clm_pf_idata%press_maxponding_clmp = 0
    clm_pf_idata%press_subsurf_pfs  = 0
    clm_pf_idata%press_subbase_pfs = 0
    clm_pf_idata%press_maxponding_pfs = 0
  ! OR, BC: water infiltration/recharge(drainage) (mH2O/sec) on the top/bottom
  !         of 3-D subsurface domain as boundary conditions from CLM to PF
    clm_pf_idata%qflux_subsurf_clmp = 0
    clm_pf_idata%qflux_subbase_clmp = 0
    clm_pf_idata%qflux_subsurf_pfs  = 0
    clm_pf_idata%qflux_subbase_pfs  = 0

  ! actual mass water flow rate (kgH2O/sec) through the top/bottom BC of 3-D subsurface domain
    clm_pf_idata%qinfl_subsurf_pfp  = 0
    clm_pf_idata%qinfl_subsurf_clms = 0
    clm_pf_idata%qsurf_subsurf_pfp  = 0
    clm_pf_idata%qsurf_subsurf_clms = 0
    clm_pf_idata%qflux_subbase_pfp  = 0
    clm_pf_idata%qflux_subbase_clms = 0

  ! (iii) Ground/base heat flux BC for PFLOTRAN's subsurface domain
  !       This BC is applied on top/bottom of the 3-D subsurface domain
    clm_pf_idata%gflux_subsurf_clmp = 0
    clm_pf_idata%gflux_subsurf_pfs  = 0
    clm_pf_idata%gflux_subbase_clmp = 0
    clm_pf_idata%gflux_subbase_pfs  = 0

    clm_pf_idata%gflux_subsurf_pfp  = 0
    clm_pf_idata%gflux_subsurf_clms = 0
    clm_pf_idata%gflux_subbase_pfp  = 0
    clm_pf_idata%gflux_subbase_clms = 0

  ! OR, if ground/bottom temperature at the subsurface interface is known
    clm_pf_idata%gtemp_subsurf_clmp = 0
    clm_pf_idata%gtemp_subsurf_pfs  = 0
    clm_pf_idata%gtemp_subbase_clmp = 0
    clm_pf_idata%gtemp_subbase_pfs  = 0

    clm_pf_idata%gtemp_subsurf_pfp  = 0
    clm_pf_idata%gtemp_subsurf_clms = 0
    clm_pf_idata%gtemp_subbase_pfp  = 0
    clm_pf_idata%gtemp_subbase_clms = 0

  ! (iv) TH state vecs from CLM (mpi) to PF (seq) - 3D cells
    clm_pf_idata%press_sub_clmp  = 0
    clm_pf_idata%soilpsi_sub_clmp = 0
    clm_pf_idata%soillsat_sub_clmp  = 0
    clm_pf_idata%soilisat_sub_clmp = 0
    clm_pf_idata%soilt_sub_clmp  = 0
    clm_pf_idata%press_sub_pfs = 0
    clm_pf_idata%soilpsi_sub_pfs = 0
    clm_pf_idata%soillsat_sub_pfs = 0
    clm_pf_idata%soilisat_sub_pfs = 0
    clm_pf_idata%soilt_sub_pfs = 0

  ! TH state vecs from PF (mpi) to CLM (seq) - 3D cells
    clm_pf_idata%press_sub_pfp    = 0
    clm_pf_idata%soilpsi_sub_pfp  = 0
    clm_pf_idata%soillsat_sub_pfp = 0
    clm_pf_idata%soilisat_sub_pfp = 0
    clm_pf_idata%soilt_sub_pfp    = 0
    clm_pf_idata%press_sub_clms   = 0
    clm_pf_idata%soilpsi_sub_clms = 0
    clm_pf_idata%soillsat_sub_clms= 0
    clm_pf_idata%soilisat_sub_clms= 0
    clm_pf_idata%soilt_sub_clms   = 0

  !-------------------------------------------------------------------------------------
  ! ground/soil C/N pools from CLM (mpi) to PF (seq)
    clm_pf_idata%ndecomp_pools = 0
    clm_pf_idata%decomp_cpools_sub_clmp  = 0
    clm_pf_idata%decomp_npools_sub_clmp  = 0
    clm_pf_idata%smin_no3_sub_clmp = 0
    clm_pf_idata%smin_nh4_sub_clmp  = 0
    clm_pf_idata%smin_nh4sorb_sub_clmp = 0
    clm_pf_idata%decomp_cpools_sub_pfs = 0
    clm_pf_idata%decomp_npools_sub_pfs  = 0
    clm_pf_idata%smin_no3_sub_pfs  = 0
    clm_pf_idata%smin_nh4_sub_pfs  = 0
    clm_pf_idata%smin_nh4sorb_sub_pfs = 0

  ! time-varying ground/soil C/N rates from CLM (mpi) to PF (seq) (Units: moleC(N)/m3/s)
    clm_pf_idata%rate_decomp_cpools_sub_clmp = 0
    clm_pf_idata%rate_decomp_npools_sub_clmp = 0
    clm_pf_idata%rate_smin_no3_sub_clmp = 0
    clm_pf_idata%rate_smin_nh4_sub_clmp = 0
    clm_pf_idata%rate_plantndemand_sub_clmp = 0
    clm_pf_idata%rate_decomp_cpools_sub_pfs = 0
    clm_pf_idata%rate_decomp_npools_sub_pfs = 0
    clm_pf_idata%rate_smin_no3_sub_pfs = 0
    clm_pf_idata%rate_smin_nh4_sub_pfs = 0
    clm_pf_idata%rate_plantndemand_sub_pfs = 0

  ! actual aqeuous N mass flow rate(gN/m2/s) at the top (runoff)/bottom (leaching) of 3-D subsurface domain (PF to CLM)
    clm_pf_idata%f_nh4_subsurf_pfp  = 0
    clm_pf_idata%f_nh4_subsurf_clms = 0
    clm_pf_idata%f_nh4_subbase_pfp = 0
    clm_pf_idata%f_nh4_subbase_clms = 0
    clm_pf_idata%f_no3_subsurf_pfp = 0
    clm_pf_idata%f_no3_subsurf_clms = 0
    clm_pf_idata%f_no3_subbase_pfp = 0
    clm_pf_idata%f_no3_subbase_clms = 0

  ! -----BGC vecs from PF (mpi, ghosted) to CLM (seq, local) --------------------
  ! BGC state variables
    clm_pf_idata%decomp_cpools_sub_pfp = 0
    clm_pf_idata%decomp_npools_sub_pfp = 0
    clm_pf_idata%smin_no3_sub_pfp      = 0
    clm_pf_idata%smin_nh4_sub_pfp      = 0
    clm_pf_idata%smin_nh4sorb_sub_pfp  = 0
  !
    clm_pf_idata%decomp_cpools_sub_clms = 0
    clm_pf_idata%decomp_npools_sub_clms = 0
    clm_pf_idata%smin_no3_sub_clms  = 0
    clm_pf_idata%smin_nh4_sub_clms  = 0
    clm_pf_idata%smin_nh4sorb_sub_clms  = 0
  !
  ! gases in water (aqueous solution of gases)
    clm_pf_idata%gco2_sub_pfp  = 0
    clm_pf_idata%gco2_sub_clms   = 0
    clm_pf_idata%gco2_sub_clmp  = 0
    clm_pf_idata%gco2_sub_pfs  = 0

    clm_pf_idata%gn2_sub_pfp  = 0
    clm_pf_idata%gn2_sub_clms  = 0
    clm_pf_idata%gn2_sub_clmp  = 0
    clm_pf_idata%gn2_sub_pfs = 0

    clm_pf_idata%gn2o_sub_pfp = 0
    clm_pf_idata%gn2o_sub_clms = 0
    clm_pf_idata%gn2o_sub_clmp  = 0
    clm_pf_idata%gn2o_sub_pfs = 0

  ! 'accextrn' is accumulative N extract by plant roots in 'PFLOTRAN' within a CLM timestep
    clm_pf_idata%accextrn_sub_pfp = 0
    clm_pf_idata%accextrn_sub_clms = 0

  ! some tracking variables from PFLOTRAN bgc to obtain reaction flux rates which needed by CLM
    clm_pf_idata%acchr_sub_pfp  = 0
    clm_pf_idata%acchr_sub_clms = 0

    clm_pf_idata%accnmin_sub_pfp  = 0
    clm_pf_idata%accnmin_sub_clms = 0

    clm_pf_idata%accnimm_sub_pfp  = 0
    clm_pf_idata%accnimm_sub_clms = 0

    clm_pf_idata%accngasmin_sub_pfp  = 0
    clm_pf_idata%accngasmin_sub_clms = 0

    clm_pf_idata%accngasnitr_sub_pfp  = 0
    clm_pf_idata%accngasnitr_sub_clms = 0

    clm_pf_idata%accngasdeni_sub_pfp  = 0
    clm_pf_idata%accngasdeni_sub_clms = 0

  end subroutine CLMPFLOTRANIDataInit

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataCreateVec(mycomm)
  ! 
  ! This routine creates PETSc vectors required for data transfer between
  ! CLM and PFLOTRAN.
  ! 

    implicit none
    
    PetscErrorCode :: ierr
    PetscMPIInt    :: mycomm, rank
    Vec :: vec_test

    call MPI_Comm_rank(mycomm,rank, ierr)

    !
    ! (i) For data transfer from CLM to PFLOTRAN
    !
    ! Create MPI Vectors for CLM
    ! 2-D surface grids
    call VecCreateMPI(mycomm, clm_pf_idata%nlclm_srf, PETSC_DECIDE, clm_pf_idata%area_srf_clmp, ierr)
    call VecSet(clm_pf_idata%area_srf_clmp,0.d0,ierr)
    ! 3-D sub cells
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_subcell, PETSC_DECIDE, clm_pf_idata%zsoi_sub_clmp, ierr)
    call VecSet(clm_pf_idata%zsoi_sub_clmp,0.d0,ierr)
    ! 2-D top/bottom face of 3-D cells
    call VecCreateMPI(mycomm, clm_pf_idata%nlclm_subsurf, PETSC_DECIDE, clm_pf_idata%area_subsurf_clmp, ierr)
    call VecSet(clm_pf_idata%area_subsurf_clmp,0.d0,ierr)
    ! vectors for multiple decomposition pools
    call VecCreateMPI(mycomm, clm_pf_idata%ndecomp_pools*clm_pf_idata%nlclm_subcell, PETSC_DECIDE, &
                      clm_pf_idata%decomp_cpools_sub_clmp, ierr)
    call VecSet(clm_pf_idata%decomp_cpools_sub_clmp, 0.d0, ierr)

    ! Create Seq. Vectors for PFLOTRAN
    ! 2-D surface grids
    call VecCreateSeq(PETSC_COMM_SELF, clm_pf_idata%ngpf_srf,clm_pf_idata%area_srf_pfs,ierr)
    call VecSet(clm_pf_idata%area_srf_pfs,0.d0,ierr)
    ! 3-D sub cells
    call VecCreateSeq(PETSC_COMM_SELF, clm_pf_idata%ngpf_subcell,clm_pf_idata%zsoi_sub_pfs,ierr)
    call VecSet(clm_pf_idata%zsoi_sub_pfs,0.d0,ierr)
    ! 2-D top/bottom face of 3-D cells
    call VecCreateSeq(PETSC_COMM_SELF, clm_pf_idata%ngpf_subsurf,clm_pf_idata%area_subsurf_pfs,ierr)
    call VecSet(clm_pf_idata%area_subsurf_pfs,0.d0,ierr)
    ! vectors for multiple decomposition pools
    call VecCreateMPI(PETSC_COMM_SELF, clm_pf_idata%ndecomp_pools*clm_pf_idata%ngpf_subcell, &
                      clm_pf_idata%decomp_cpools_sub_pfs, ierr)
    call VecSet(clm_pf_idata%decomp_cpools_sub_pfs, 0.d0, ierr)

    ! (ii) For data transfer from PFLOTRAN to CLM
    !
    ! Create MPI Vectors for PFLOTRAN
    ! 2-D surface grids
    call VecCreateMPI(mycomm, clm_pf_idata%nlpf_srf, PETSC_DECIDE, clm_pf_idata%area_srf_pfp, ierr)
    call VecSet(clm_pf_idata%area_srf_pfp,0.d0,ierr)
    ! 3-D sub cells
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_subcell, PETSC_DECIDE, clm_pf_idata%zsoi_sub_pfp, ierr)
    call VecSet(clm_pf_idata%zsoi_sub_pfp,0.d0,ierr)
    ! 2-D top/bottom face of 3-D cells
    call VecCreateMPI(mycomm, clm_pf_idata%nlpf_subsurf, PETSC_DECIDE, clm_pf_idata%area_subsurf_pfp, ierr)
    call VecSet(clm_pf_idata%area_subsurf_pfp,0.d0,ierr)
    ! vectors for multiple decomposition pools
    call VecCreateMPI(mycomm, clm_pf_idata%ndecomp_pools*clm_pf_idata%nlpf_subcell, PETSC_DECIDE, &
                      clm_pf_idata%decomp_cpools_sub_pfp, ierr)
    call VecSet(clm_pf_idata%decomp_cpools_sub_pfp, 0.d0, ierr)

    ! Create Seq. Vectors for CLM
    ! 2-D surface grids
    call VecCreateSeq(PETSC_COMM_SELF, clm_pf_idata%ngclm_srf,clm_pf_idata%area_srf_clms,ierr)
    call VecSet(clm_pf_idata%area_srf_clms,0.d0,ierr)
    ! 3-D sub cells
    call VecCreateSeq(PETSC_COMM_SELF, clm_pf_idata%ngclm_subcell,clm_pf_idata%zsoi_sub_clms,ierr)
    call VecSet(clm_pf_idata%zsoi_sub_clms,0.d0,ierr)
    ! 2-D top/bottom face of 3-D cells
    call VecCreateSeq(PETSC_COMM_SELF, clm_pf_idata%ngclm_subsurf,clm_pf_idata%area_subsurf_clms,ierr)
    call VecSet(clm_pf_idata%area_subsurf_clms,0.d0,ierr)
    ! vectors for multiple decomposition pools
    call VecCreateMPI(PETSC_COMM_SELF, clm_pf_idata%ndecomp_pools*clm_pf_idata%nlclm_subcell, &
                      clm_pf_idata%decomp_cpools_sub_clms, ierr)
    call VecSet(clm_pf_idata%decomp_cpools_sub_clms, 0.d0, ierr)

    ! (iii) then, duplicate vectors for the rests of ALL

   !--------------------------------------------------------------------
   ! Soil properties -
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%hksat_x_clmp ,ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%hksat_y_clmp ,ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%hksat_z_clmp ,ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%bulkdensity_dry_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%effporosity_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%press_ref_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%watsat_sub_clmp   , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%watfc_sub_clmp    , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%watpcwmax_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%sucsat_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%bsw_sub_clmp   , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%sr_sub_clmp    , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%lamda_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%alpha_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%pcwmax_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%tcond_dry_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%tcond_wet_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%hc_dry_sub_clmp    , ierr)

    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%hksat_x_pfs , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%hksat_y_pfs , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%hksat_z_pfs , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%bulkdensity_dry_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%effporosity_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%press_ref_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%watsat_sub_pfs   , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%watfc_sub_pfs    , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%watpcwmax_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%sucsat_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%bsw_sub_pfs   , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%sr_sub_pfs   , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%lamda_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%alpha_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%pcwmax_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%tcond_dry_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%tcond_wet_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%hc_dry_sub_pfs   , ierr)

    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%hksat_x_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%hksat_y_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%hksat_z_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%bulkdensity_dry_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%effporosity_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%press_ref_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%watsat_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%watfc_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%watpcwmax_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%sucsat_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%bsw_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%sr_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%lamda_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%alpha_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%pcwmax_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%tcond_dry_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%tcond_wet_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%hc_dry_sub_pfp, ierr)

    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%hksat_x_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%hksat_y_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%hksat_z_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%bulkdensity_dry_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%effporosity_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%press_ref_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%watsat_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%watfc_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%watpcwmax_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%sucsat_sub_clms , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%bsw_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%sr_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%lamda_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%alpha_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%pcwmax_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%tcond_dry_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%tcond_wet_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%hc_dry_sub_clms, ierr)

  ! Time variant data

  ! Sink/Source of water/heat for PFLOTRAN's 3D subsurface domain
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%qflux_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%qflux_sub_pfs , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%qflux_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%qflux_sub_clms, ierr)

    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%gflux_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%gflux_sub_pfs , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%gflux_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%gflux_sub_clms, ierr)

  ! (BC) surface water I/O for PFLOTRAN's 2D surface domain
    call VecDuplicate(clm_pf_idata%area_srf_clmp,clm_pf_idata%h2o_srf_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_srf_pfs,clm_pf_idata%h2o_srf_pfs, ierr)
    call VecDuplicate(clm_pf_idata%area_srf_clmp,clm_pf_idata%qh2o_srf_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_srf_pfs,clm_pf_idata%qh2o_srf_pfs, ierr)
    call VecDuplicate(clm_pf_idata%area_srf_clmp,clm_pf_idata%h2otemp_srf_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_srf_pfs,clm_pf_idata%h2otemp_srf_pfs, ierr)

    call VecDuplicate(clm_pf_idata%area_srf_pfp,clm_pf_idata%h2o_srf_pfp, ierr)
    call VecDuplicate(clm_pf_idata%area_srf_clms,clm_pf_idata%h2o_srf_clms, ierr)
    call VecDuplicate(clm_pf_idata%area_srf_pfp,clm_pf_idata%qh2o_srf_pfp, ierr)
    call VecDuplicate(clm_pf_idata%area_srf_clms,clm_pf_idata%qh2o_srf_clms, ierr)
    call VecDuplicate(clm_pf_idata%area_srf_pfp,clm_pf_idata%h2otemp_srf_pfp, ierr)
    call VecDuplicate(clm_pf_idata%area_srf_clms,clm_pf_idata%h2otemp_srf_clms, ierr)

  ! BC: water pressure (Pa) on the top/bottom of 3-D subsurface domain
  !     as boundary conditions from CLM to PF
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%press_subsurf_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%press_subbase_clmp , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%press_maxponding_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%press_subsurf_pfs , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%press_subbase_pfs, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%press_maxponding_pfs, ierr)
  ! OR, BC: water infiltration/recharge(drainage) (mH2O/sec) on the top/bottom
  !         of 3-D subsurface domain as boundary conditions from CLM to PF
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%qflux_subsurf_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%qflux_subbase_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%qflux_subsurf_pfs , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%qflux_subbase_pfs , ierr)

  ! actual mass water flow rate (kgH2O/sec) through the top/bottom BC of 3-D subsurface domain
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp, clm_pf_idata%qinfl_subsurf_pfp , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms, clm_pf_idata%qinfl_subsurf_clms, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp, clm_pf_idata%qsurf_subsurf_pfp , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms, clm_pf_idata%qsurf_subsurf_clms, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp, clm_pf_idata%qflux_subbase_pfp , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms, clm_pf_idata%qflux_subbase_clms, ierr)

  ! Ground/base heat flux BC for PFLOTRAN's subsurface domain
  !       This BC is applied on top/bottom of the 3-D subsurface domain
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%gflux_subsurf_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%gflux_subsurf_pfs , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%gflux_subbase_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%gflux_subbase_pfs , ierr)

    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%gflux_subsurf_pfp , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%gflux_subsurf_clms, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%gflux_subbase_pfp , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%gflux_subbase_clms, ierr)
  ! OR, if ground/bottom temperature at the subsurface interface is known
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%gtemp_subsurf_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%gtemp_subsurf_pfs , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%gtemp_subbase_clmp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%gtemp_subbase_pfs , ierr)

    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%gtemp_subsurf_pfp , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%gtemp_subsurf_clms, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%gtemp_subbase_pfp , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%gtemp_subbase_clms, ierr)

  ! TH state vecs from CLM (mpi) to PF (seq) - 3D cells
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%press_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%soilpsi_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%soillsat_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%soilisat_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%soilt_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%press_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%soilpsi_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%soillsat_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%soilisat_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%soilt_sub_pfs, ierr)

  ! TH state vecs from PF (mpi) to CLM (seq) - 3D cells
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%press_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%soilpsi_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%soillsat_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%soilisat_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%soilt_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%press_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%soilpsi_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%soillsat_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%soilisat_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%soilt_sub_clms, ierr)

  !----------------------------------------------------------------------------------------------------------------
  ! ground/soil C/N pools from CLM (mpi) to PF (seq)
    call VecDuplicate(clm_pf_idata%decomp_cpools_sub_clmp,clm_pf_idata%decomp_npools_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%smin_no3_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%smin_nh4_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%smin_nh4sorb_sub_clmp, ierr)

    call VecDuplicate(clm_pf_idata%decomp_cpools_sub_pfs,clm_pf_idata%decomp_npools_sub_pfs , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%smin_no3_sub_pfs , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%smin_nh4_sub_pfs , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%smin_nh4sorb_sub_pfs, ierr)

  ! time-varying ground/soil C/N rates from CLM (mpi) to PF (seq) (Units: moleC(N)/m3/s)
    call VecDuplicate(clm_pf_idata%decomp_cpools_sub_clmp,clm_pf_idata%rate_decomp_cpools_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_sub_clmp,clm_pf_idata%rate_decomp_npools_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%rate_smin_no3_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%rate_smin_nh4_sub_clmp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%rate_plantndemand_sub_clmp, ierr)

    call VecDuplicate(clm_pf_idata%decomp_cpools_sub_pfs,clm_pf_idata%rate_decomp_cpools_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_sub_pfs,clm_pf_idata%rate_decomp_npools_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%rate_smin_no3_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%rate_smin_nh4_sub_pfs, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%rate_plantndemand_sub_pfs, ierr)

  ! actual aqeuous N mass flow rate(gN/m2/s) at the top (runoff)/bottom (leaching) of 3-D subsurface domain (PF to CLM)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%f_nh4_subsurf_pfp , ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%f_nh4_subsurf_clms, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%f_nh4_subbase_pfp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%f_nh4_subbase_clms, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%f_no3_subsurf_pfp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%f_no3_subsurf_clms, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%f_no3_subbase_pfp, ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%f_no3_subbase_clms, ierr)

  ! -----BGC vecs from PF (mpi, ghosted) to CLM (seq, local) --------------------
  ! BGC state variables
    call VecDuplicate(clm_pf_idata%decomp_cpools_sub_pfp,clm_pf_idata%decomp_npools_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%smin_no3_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%smin_nh4_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%smin_nh4sorb_sub_pfp , ierr)
  !
    call VecDuplicate(clm_pf_idata%decomp_cpools_sub_clms,clm_pf_idata%decomp_npools_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%smin_no3_sub_clms , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%smin_nh4_sub_clms , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%smin_nh4sorb_sub_clms , ierr)
  !
  ! gases in water (aqueous solution of gases)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%gco2_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%gco2_sub_clms  , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%gco2_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%gco2_sub_pfs , ierr)

    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%gn2_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%gn2_sub_clms , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%gn2_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%gn2_sub_pfs, ierr)

    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%gn2o_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%gn2o_sub_clms, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clmp,clm_pf_idata%gn2o_sub_clmp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfs,clm_pf_idata%gn2o_sub_pfs, ierr)

  ! 'accextrn' is accumulative N extract by plant roots in 'PFLOTRAN' within a CLM timestep
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%accextrn_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%accextrn_sub_clms, ierr)

  ! some tracking variables from PFLOTRAN bgc to obtain reaction flux rates which needed by CLM
    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%acchr_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%acchr_sub_clms, ierr)

    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%accnmin_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%accnmin_sub_clms, ierr)

    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%accnimm_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%accnimm_sub_clms, ierr)

    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%accngasmin_sub_pfp , ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%accngasmin_sub_clms  , ierr)

    call VecDuplicate(clm_pf_idata%zsoi_sub_pfp,clm_pf_idata%accngasnitr_sub_pfp, ierr)
    call VecDuplicate(clm_pf_idata%zsoi_sub_clms,clm_pf_idata%accngasnitr_sub_clms, ierr)

    !---------------------------------------------

  end subroutine CLMPFLOTRANIDataCreateVec

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataDestroy()
  ! 
  ! This routine destroys PETSc vectors that were created for data transfer.
  ! 

    implicit none
    
    PetscErrorCode :: ierr

   !--------------------------------------------------------------------

   if(clm_pf_idata%area_srf_clmp /= 0) &
     call VecDestroy(clm_pf_idata%area_srf_clmp, ierr)
   if(clm_pf_idata%area_srf_pfs /= 0) &
     call VecDestroy(clm_pf_idata%area_srf_pfs , ierr)
   if(clm_pf_idata%area_srf_pfp  /= 0) &
     call VecDestroy(clm_pf_idata%area_srf_pfp , ierr)
   if(clm_pf_idata%area_srf_clms  /= 0) &
     call VecDestroy(clm_pf_idata%area_srf_clms, ierr)

   if(clm_pf_idata%zsoi_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%zsoi_sub_clmp, ierr)
   if(clm_pf_idata%zsoi_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%zsoi_sub_pfs, ierr)
   if(clm_pf_idata%zsoi_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%zsoi_sub_pfp, ierr)
   if(clm_pf_idata%zsoi_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%zsoi_sub_clms, ierr)

   if(clm_pf_idata%area_subsurf_clmp /= 0) &
     call VecDestroy(clm_pf_idata%area_subsurf_clmp, ierr)
   if(clm_pf_idata%area_subsurf_pfs /= 0) &
     call VecDestroy(clm_pf_idata%area_subsurf_pfs, ierr)
   if(clm_pf_idata%area_subsurf_pfp /= 0) &
     call VecDestroy(clm_pf_idata%area_subsurf_pfp, ierr)
   if(clm_pf_idata%area_subsurf_clms /= 0) &
     call VecDestroy(clm_pf_idata%area_subsurf_clms, ierr)

   !--------------------------------------------------------------------
   ! Soil properties -
    if(clm_pf_idata%hksat_x_clmp /= 0) &
     call VecDestroy(clm_pf_idata%hksat_x_clmp, ierr)
    if(clm_pf_idata%hksat_y_clmp /= 0) &
     call VecDestroy(clm_pf_idata%hksat_y_clmp, ierr)
    if(clm_pf_idata%hksat_z_clmp  /= 0) &
     call VecDestroy(clm_pf_idata%hksat_z_clmp, ierr)
    if(clm_pf_idata%bulkdensity_dry_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%bulkdensity_dry_sub_clmp, ierr)
    if(clm_pf_idata%effporosity_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%effporosity_sub_clmp , ierr)
    if(clm_pf_idata%press_ref_clmp /= 0) &
     call VecDestroy(clm_pf_idata%press_ref_clmp, ierr)
    if(clm_pf_idata%watsat_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%watsat_sub_clmp, ierr)
    if(clm_pf_idata%watfc_sub_clmp  /= 0) &
     call VecDestroy(clm_pf_idata%watfc_sub_clmp, ierr)
    if(clm_pf_idata%watpcwmax_sub_clmp  /= 0) &
     call VecDestroy(clm_pf_idata%watpcwmax_sub_clmp, ierr)
    if(clm_pf_idata%sucsat_sub_clmp  /= 0) &
     call VecDestroy(clm_pf_idata%sucsat_sub_clmp, ierr)
    if(clm_pf_idata%bsw_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%bsw_sub_clmp   , ierr)
    if(clm_pf_idata%sr_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%sr_sub_clmp    , ierr)
    if(clm_pf_idata%lamda_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%lamda_sub_clmp , ierr)
    if(clm_pf_idata%alpha_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%alpha_sub_clmp , ierr)
    if(clm_pf_idata%pcwmax_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%pcwmax_sub_clmp, ierr)
    if(clm_pf_idata%tcond_dry_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%tcond_dry_sub_clmp , ierr)
    if(clm_pf_idata%tcond_wet_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%tcond_wet_sub_clmp , ierr)
    if(clm_pf_idata%hc_dry_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%hc_dry_sub_clmp    , ierr)

    if(clm_pf_idata%hksat_x_pfs /= 0) &
     call VecDestroy(clm_pf_idata%hksat_x_pfs , ierr)
    if(clm_pf_idata%hksat_y_pfs /= 0) &
     call VecDestroy(clm_pf_idata%hksat_y_pfs , ierr)
    if(clm_pf_idata%hksat_z_pfs /= 0) &
     call VecDestroy(clm_pf_idata%hksat_z_pfs , ierr)
    if(clm_pf_idata%bulkdensity_dry_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%bulkdensity_dry_sub_pfs, ierr)
    if(clm_pf_idata%effporosity_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%effporosity_sub_pfs, ierr)
    if(clm_pf_idata%press_ref_pfs /= 0) &
     call VecDestroy(clm_pf_idata%press_ref_pfs, ierr)
    if(clm_pf_idata%watsat_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%watsat_sub_pfs, ierr)
    if(clm_pf_idata%watfc_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%watfc_sub_pfs, ierr)
    if(clm_pf_idata%watfc_sub_pfs  /= 0) &
     call VecDestroy(clm_pf_idata%watfc_sub_pfs, ierr)
    if(clm_pf_idata%sucsat_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%sucsat_sub_pfs, ierr)
    if(clm_pf_idata%bsw_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%bsw_sub_pfs, ierr)
    if(clm_pf_idata%sr_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%sr_sub_pfs, ierr)
    if(clm_pf_idata%lamda_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%lamda_sub_pfs, ierr)
    if(clm_pf_idata%alpha_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%alpha_sub_pfs, ierr)
    if(clm_pf_idata%pcwmax_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%pcwmax_sub_pfs, ierr)
    if(clm_pf_idata%tcond_dry_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%tcond_dry_sub_pfs, ierr)
    if(clm_pf_idata%tcond_wet_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%tcond_wet_sub_pfs, ierr)
    if(clm_pf_idata%hc_dry_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%hc_dry_sub_pfs, ierr)

    if(clm_pf_idata%hksat_x_pfp /= 0) &
     call VecDestroy(clm_pf_idata%hksat_x_pfp , ierr)
    if(clm_pf_idata%hksat_y_pfp /= 0) &
     call VecDestroy(clm_pf_idata%hksat_y_pfp, ierr)
    if(clm_pf_idata%hksat_z_pfp /= 0) &
     call VecDestroy(clm_pf_idata%hksat_z_pfp, ierr)
    if(clm_pf_idata%bulkdensity_dry_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%bulkdensity_dry_sub_pfp, ierr)
    if(clm_pf_idata%effporosity_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%effporosity_sub_pfp, ierr)
    if(clm_pf_idata%press_ref_pfp /= 0) &
     call VecDestroy(clm_pf_idata%press_ref_pfp, ierr)
    if(clm_pf_idata%watsat_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%watsat_sub_pfp, ierr)
    if(clm_pf_idata%watfc_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%watfc_sub_pfp, ierr)
    if(clm_pf_idata%watpcwmax_sub_pfp  /= 0) &
     call VecDestroy(clm_pf_idata%watpcwmax_sub_pfp, ierr)
    if(clm_pf_idata%sucsat_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%sucsat_sub_pfp , ierr)
    if(clm_pf_idata%bsw_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%bsw_sub_pfp, ierr)
    if(clm_pf_idata%sr_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%sr_sub_pfp, ierr)
    if(clm_pf_idata%lamda_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%lamda_sub_pfp, ierr)
    if(clm_pf_idata%alpha_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%alpha_sub_pfp, ierr)
    if(clm_pf_idata%pcwmax_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%pcwmax_sub_pfp, ierr)
    if(clm_pf_idata%tcond_dry_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%tcond_dry_sub_pfp , ierr)
    if(clm_pf_idata%tcond_wet_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%tcond_wet_sub_pfp , ierr)
    if(clm_pf_idata%hc_dry_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%hc_dry_sub_pfp, ierr)

    if(clm_pf_idata%hksat_x_clms /= 0) &
     call VecDestroy(clm_pf_idata%hksat_x_clms, ierr)
    if(clm_pf_idata%hksat_y_clms /= 0) &
     call VecDestroy(clm_pf_idata%hksat_y_clms, ierr)
    if(clm_pf_idata%hksat_z_clms /= 0) &
     call VecDestroy(clm_pf_idata%hksat_z_clms, ierr)
    if(clm_pf_idata%bulkdensity_dry_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%bulkdensity_dry_sub_clms, ierr)
    if(clm_pf_idata%effporosity_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%effporosity_sub_clms, ierr)
    if(clm_pf_idata%press_ref_clms /= 0) &
     call VecDestroy(clm_pf_idata%press_ref_clms, ierr)
    if(clm_pf_idata%watsat_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%watsat_sub_clms, ierr)
    if(clm_pf_idata%watfc_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%watfc_sub_clms, ierr)
    if(clm_pf_idata%watpcwmax_sub_clms  /= 0) &
     call VecDestroy(clm_pf_idata%watpcwmax_sub_clmp, ierr)
    if(clm_pf_idata%sucsat_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%sucsat_sub_clms , ierr)
    if(clm_pf_idata%bsw_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%bsw_sub_clms, ierr)
    if(clm_pf_idata%sr_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%sr_sub_clms, ierr)
    if(clm_pf_idata%lamda_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%lamda_sub_clms, ierr)
    if(clm_pf_idata%alpha_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%alpha_sub_clms, ierr)
    if(clm_pf_idata%pcwmax_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%pcwmax_sub_clms, ierr)
    if(clm_pf_idata%tcond_dry_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%tcond_dry_sub_clms, ierr)
    if(clm_pf_idata%tcond_wet_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%tcond_wet_sub_clms, ierr)
    if(clm_pf_idata%hc_dry_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%hc_dry_sub_clms, ierr)

  ! Time variant data

  ! (i) Sink/Source of water/heat for PFLOTRAN's 3D subsurface domain
    if(clm_pf_idata%qflux_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%qflux_sub_clmp, ierr)
    if(clm_pf_idata%qflux_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%qflux_sub_pfs , ierr)
    if(clm_pf_idata%qflux_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%qflux_sub_pfp , ierr)
    if(clm_pf_idata%qflux_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%qflux_sub_clms, ierr)

    if(clm_pf_idata%gflux_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%gflux_sub_clmp, ierr)
    if(clm_pf_idata%gflux_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%gflux_sub_pfs , ierr)
    if(clm_pf_idata%gflux_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%gflux_sub_pfp , ierr)
    if(clm_pf_idata%gflux_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%gflux_sub_clms, ierr)

  ! (ii) surface water I/O for PFLOTRAN's 2D surface domain
    if(clm_pf_idata%h2o_srf_clmp /= 0) &
     call VecDestroy(clm_pf_idata%h2o_srf_clmp, ierr)
    if(clm_pf_idata%h2o_srf_pfs /= 0) &
     call VecDestroy(clm_pf_idata%h2o_srf_pfs, ierr)
    if(clm_pf_idata%qh2o_srf_clmp /= 0) &
     call VecDestroy(clm_pf_idata%qh2o_srf_clmp, ierr)
    if(clm_pf_idata%qh2o_srf_pfs /= 0) &
     call VecDestroy(clm_pf_idata%qh2o_srf_pfs, ierr)
    if(clm_pf_idata%h2otemp_srf_clmp /= 0) &
     call VecDestroy(clm_pf_idata%h2otemp_srf_clmp, ierr)
    if(clm_pf_idata%h2otemp_srf_pfs /= 0) &
     call VecDestroy(clm_pf_idata%h2otemp_srf_pfs, ierr)

    if(clm_pf_idata%h2o_srf_pfp /= 0) &
     call VecDestroy(clm_pf_idata%h2o_srf_pfp, ierr)
    if(clm_pf_idata%h2o_srf_clms /= 0) &
     call VecDestroy(clm_pf_idata%h2o_srf_clms, ierr)
    if(clm_pf_idata%qh2o_srf_pfp /= 0) &
     call VecDestroy(clm_pf_idata%qh2o_srf_pfp, ierr)
    if(clm_pf_idata%qh2o_srf_clms /= 0) &
     call VecDestroy(clm_pf_idata%qh2o_srf_clms, ierr)
    if(clm_pf_idata%h2otemp_srf_pfp /= 0) &
     call VecDestroy(clm_pf_idata%h2otemp_srf_pfp, ierr)
    if(clm_pf_idata%h2otemp_srf_clms /= 0) &
     call VecDestroy(clm_pf_idata%h2otemp_srf_clms, ierr)

  ! BC: water pressure (Pa) on the top/bottom of 3-D subsurface domain
  !     as boundary conditions from CLM to PF
    if(clm_pf_idata%press_subsurf_clmp /= 0) &
     call VecDestroy(clm_pf_idata%press_subsurf_clmp, ierr)
    if(clm_pf_idata%press_subbase_clmp /= 0) &
     call VecDestroy(clm_pf_idata%press_subbase_clmp , ierr)
    if(clm_pf_idata%press_maxponding_clmp /= 0) &
     call VecDestroy(clm_pf_idata%press_maxponding_clmp, ierr)
    if(clm_pf_idata%press_subbase_pfs /= 0) &
     call VecDestroy(clm_pf_idata%press_subsurf_pfs , ierr)
    if(clm_pf_idata%press_subbase_pfs /= 0) &
     call VecDestroy(clm_pf_idata%press_subbase_pfs, ierr)
    if(clm_pf_idata%press_maxponding_pfs /= 0) &
     call VecDestroy(clm_pf_idata%press_maxponding_pfs, ierr)
  ! OR, BC: water infiltration/recharge(drainage) (mH2O/sec) on the top/bottom
  !         of 3-D subsurface domain as boundary conditions from CLM to PF
    if(clm_pf_idata%qflux_subsurf_clmp /= 0) &
     call VecDestroy(clm_pf_idata%qflux_subsurf_clmp, ierr)
    if(clm_pf_idata%qflux_subbase_clmp /= 0) &
     call VecDestroy(clm_pf_idata%qflux_subbase_clmp, ierr)
    if(clm_pf_idata%qflux_subsurf_pfs /= 0) &
     call VecDestroy(clm_pf_idata%qflux_subsurf_pfs , ierr)
    if(clm_pf_idata%qflux_subbase_pfs /= 0) &
     call VecDestroy(clm_pf_idata%qflux_subbase_pfs , ierr)

  ! actual mass water flow rate (kgH2O/sec) through the top/bottom BC of 3-D subsurface domain
    if(clm_pf_idata%qinfl_subsurf_pfp /= 0) &
     call VecDestroy(clm_pf_idata%qinfl_subsurf_pfp , ierr)
    if(clm_pf_idata%qinfl_subsurf_clms /= 0) &
     call VecDestroy(clm_pf_idata%qinfl_subsurf_clms, ierr)
    if(clm_pf_idata%qsurf_subsurf_pfp /= 0) &
     call VecDestroy(clm_pf_idata%qsurf_subsurf_pfp , ierr)
    if(clm_pf_idata%qsurf_subsurf_clms /= 0) &
     call VecDestroy(clm_pf_idata%qsurf_subsurf_clms, ierr)
    if(clm_pf_idata%qflux_subbase_pfp /= 0) &
     call VecDestroy(clm_pf_idata%qflux_subbase_pfp , ierr)
    if(clm_pf_idata%qflux_subbase_clms /= 0) &
     call VecDestroy(clm_pf_idata%qflux_subbase_clms, ierr)

  ! (iii) Ground/base heat flux BC for PFLOTRAN's subsurface domain
  !       This BC is applied on top/bottom of the 3-D subsurface domain
    if(clm_pf_idata%gflux_subsurf_clmp /= 0) &
     call VecDestroy(clm_pf_idata%gflux_subsurf_clmp, ierr)
    if(clm_pf_idata%gflux_subsurf_pfs /= 0) &
     call VecDestroy(clm_pf_idata%gflux_subsurf_pfs, ierr)
    if(clm_pf_idata%gflux_subbase_clmp /= 0) &
     call VecDestroy(clm_pf_idata%gflux_subbase_clmp, ierr)
    if(clm_pf_idata%gflux_subbase_pfs /= 0) &
     call VecDestroy(clm_pf_idata%gflux_subbase_pfs, ierr)

    if(clm_pf_idata%gflux_subsurf_pfp /= 0) &
     call VecDestroy(clm_pf_idata%gflux_subsurf_pfp , ierr)
    if(clm_pf_idata%gflux_subsurf_clms /= 0) &
     call VecDestroy(clm_pf_idata%gflux_subsurf_clms, ierr)
    if(clm_pf_idata%gflux_subbase_pfp /= 0) &
     call VecDestroy(clm_pf_idata%gflux_subbase_pfp , ierr)
    if(clm_pf_idata%gflux_subbase_clms /= 0) &
     call VecDestroy(clm_pf_idata%gflux_subbase_clms, ierr)

  ! OR, if ground/bottom temperature at the subsurface interface is known
    if(clm_pf_idata%gtemp_subsurf_clmp /= 0) &
     call VecDestroy(clm_pf_idata%gtemp_subsurf_clmp, ierr)
    if(clm_pf_idata%gtemp_subsurf_pfs /= 0) &
     call VecDestroy(clm_pf_idata%gtemp_subsurf_pfs , ierr)
    if(clm_pf_idata%gtemp_subbase_clmp /= 0) &
     call VecDestroy(clm_pf_idata%gtemp_subbase_clmp, ierr)
    if(clm_pf_idata%gtemp_subbase_pfs /= 0) &
     call VecDestroy(clm_pf_idata%gtemp_subbase_pfs , ierr)

    if(clm_pf_idata%gtemp_subsurf_pfp /= 0) &
     call VecDestroy(clm_pf_idata%gtemp_subsurf_pfp , ierr)
    if(clm_pf_idata%gtemp_subsurf_clms /= 0) &
     call VecDestroy(clm_pf_idata%gtemp_subsurf_clms, ierr)
    if(clm_pf_idata%gtemp_subbase_pfp /= 0) &
     call VecDestroy(clm_pf_idata%gtemp_subbase_pfp , ierr)
    if(clm_pf_idata%gtemp_subbase_clms /= 0) &
     call VecDestroy(clm_pf_idata%gtemp_subbase_clms, ierr)

  ! TH state vecs from CLM (mpi) to PF (seq) - 3D cells
    if(clm_pf_idata%press_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%press_sub_clmp , ierr)
    if(clm_pf_idata%soilpsi_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%soilpsi_sub_clmp, ierr)
    if(clm_pf_idata%soillsat_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%soillsat_sub_clmp , ierr)
    if(clm_pf_idata%soilisat_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%soilisat_sub_clmp, ierr)
    if(clm_pf_idata%soilt_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%soilt_sub_clmp , ierr)
    if(clm_pf_idata%press_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%press_sub_pfs, ierr)
    if(clm_pf_idata%soilpsi_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%soilpsi_sub_pfs, ierr)
    if(clm_pf_idata%soillsat_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%soillsat_sub_pfs, ierr)
    if(clm_pf_idata%soilisat_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%soilisat_sub_pfs, ierr)
    if(clm_pf_idata%soilt_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%soilt_sub_pfs, ierr)

  ! TH state vecs from PF (mpi) to CLM (seq) - 3D cells
    if(clm_pf_idata%press_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%press_sub_pfp , ierr)
    if(clm_pf_idata%soilpsi_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%soilpsi_sub_pfp , ierr)
    if(clm_pf_idata%soillsat_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%soillsat_sub_pfp, ierr)
    if(clm_pf_idata%soilisat_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%soilisat_sub_pfp, ierr)
    if(clm_pf_idata%soilt_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%soilt_sub_pfp , ierr)
    if(clm_pf_idata%press_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%press_sub_clms, ierr)
    if(clm_pf_idata%soilpsi_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%soilpsi_sub_clms, ierr)
    if(clm_pf_idata%soillsat_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%soillsat_sub_clms, ierr)
    if(clm_pf_idata%soilisat_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%soilisat_sub_clms, ierr)
    if(clm_pf_idata%soilt_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%soilt_sub_clms, ierr)

  !------------------------------------------------------------------------------------------------------
  ! ground/soil C/N pools from CLM (mpi) to PF (seq)
    if(clm_pf_idata%decomp_cpools_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%decomp_cpools_sub_clmp , ierr)
    if(clm_pf_idata%decomp_npools_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%decomp_npools_sub_clmp , ierr)
    if(clm_pf_idata%smin_no3_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%smin_no3_sub_clmp, ierr)
    if(clm_pf_idata%smin_nh4_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%smin_nh4_sub_clmp , ierr)
    if(clm_pf_idata%smin_nh4sorb_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%smin_nh4sorb_sub_clmp, ierr)
    if(clm_pf_idata%decomp_cpools_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%decomp_cpools_sub_pfs, ierr)
    if(clm_pf_idata%decomp_npools_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%decomp_npools_sub_pfs , ierr)
    if(clm_pf_idata%smin_no3_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%smin_no3_sub_pfs , ierr)
    if(clm_pf_idata%smin_nh4_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%smin_nh4_sub_pfs , ierr)
    if(clm_pf_idata%smin_nh4sorb_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%smin_nh4sorb_sub_pfs, ierr)

  ! time-varying ground/soil C/N rates from CLM (mpi) to PF (seq) (Units: moleC(N)/m3/s)
    if(clm_pf_idata%rate_decomp_cpools_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%rate_decomp_cpools_sub_clmp, ierr)
    if(clm_pf_idata%rate_decomp_npools_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%rate_decomp_npools_sub_clmp, ierr)
    if(clm_pf_idata%rate_smin_no3_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%rate_smin_no3_sub_clmp, ierr)
    if(clm_pf_idata%rate_smin_nh4_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%rate_smin_nh4_sub_clmp, ierr)
    if(clm_pf_idata%rate_plantndemand_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%rate_plantndemand_sub_clmp, ierr)
    if(clm_pf_idata%rate_decomp_cpools_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%rate_decomp_cpools_sub_pfs, ierr)
    if(clm_pf_idata%rate_decomp_npools_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%rate_decomp_npools_sub_pfs, ierr)
    if(clm_pf_idata%rate_smin_no3_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%rate_smin_no3_sub_pfs, ierr)
    if(clm_pf_idata%rate_smin_nh4_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%rate_smin_nh4_sub_pfs, ierr)
    if(clm_pf_idata%rate_plantndemand_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%rate_plantndemand_sub_pfs, ierr)

  ! actual aqeuous N mass flow rate(gN/m2/s) at the top (runoff)/bottom (leaching) of 3-D subsurface domain (PF to CLM)
    if(clm_pf_idata%f_nh4_subsurf_pfp /= 0) &
     call VecDestroy(clm_pf_idata%f_nh4_subsurf_pfp , ierr)
    if(clm_pf_idata%f_nh4_subsurf_clms /= 0) &
     call VecDestroy(clm_pf_idata%f_nh4_subsurf_clms, ierr)
    if(clm_pf_idata%f_nh4_subbase_pfp /= 0) &
     call VecDestroy(clm_pf_idata%f_nh4_subbase_pfp, ierr)
    if(clm_pf_idata%f_no3_subbase_clms /= 0) &
     call VecDestroy(clm_pf_idata%f_nh4_subbase_clms, ierr)
    if(clm_pf_idata%f_no3_subsurf_pfp /= 0) &
     call VecDestroy(clm_pf_idata%f_no3_subsurf_pfp, ierr)
    if(clm_pf_idata%f_no3_subsurf_clms /= 0) &
     call VecDestroy(clm_pf_idata%f_no3_subsurf_clms, ierr)
    if(clm_pf_idata%f_no3_subbase_pfp /= 0) &
     call VecDestroy(clm_pf_idata%f_no3_subbase_pfp, ierr)
    if(clm_pf_idata%f_no3_subbase_clms /= 0) &
     call VecDestroy(clm_pf_idata%f_no3_subbase_clms, ierr)

  ! -----BGC vecs from PF (mpi, ghosted) to CLM (seq, local) --------------------
  ! BGC state variables
    if(clm_pf_idata%decomp_cpools_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%decomp_cpools_sub_pfp, ierr)
    if(clm_pf_idata%decomp_npools_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%decomp_npools_sub_pfp, ierr)
    if(clm_pf_idata%smin_no3_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%smin_no3_sub_pfp, ierr)
    if(clm_pf_idata%smin_nh4_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%smin_nh4_sub_pfp, ierr)
    if(clm_pf_idata%smin_nh4sorb_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%smin_nh4sorb_sub_pfp, ierr)
  !
    if(clm_pf_idata%decomp_cpools_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%decomp_cpools_sub_clms, ierr)
    if(clm_pf_idata%decomp_npools_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%decomp_npools_sub_clms, ierr)
    if(clm_pf_idata%smin_no3_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%smin_no3_sub_clms, ierr)
    if(clm_pf_idata%smin_nh4_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%smin_nh4_sub_clms, ierr)
    if(clm_pf_idata%smin_nh4sorb_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%smin_nh4sorb_sub_clms, ierr)
  !
  ! gases in water (aqueous solution of gases)
    if(clm_pf_idata%gco2_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%gco2_sub_pfp, ierr)
    if(clm_pf_idata%gco2_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%gco2_sub_clms  , ierr)
    if(clm_pf_idata%gco2_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%gco2_sub_clmp , ierr)
    if(clm_pf_idata%gco2_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%gco2_sub_pfs , ierr)

    if(clm_pf_idata%gn2_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%gn2_sub_pfp , ierr)
    if(clm_pf_idata%gn2_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%gn2_sub_clms , ierr)
    if(clm_pf_idata%gn2_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%gn2_sub_clmp , ierr)
    if(clm_pf_idata%gn2_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%gn2_sub_pfs, ierr)

    if(clm_pf_idata%gn2o_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%gn2o_sub_pfp, ierr)
    if(clm_pf_idata%gn2o_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%gn2o_sub_clms, ierr)
    if(clm_pf_idata%gn2o_sub_clmp /= 0) &
     call VecDestroy(clm_pf_idata%gn2o_sub_clmp , ierr)
    if(clm_pf_idata%gn2o_sub_pfs /= 0) &
     call VecDestroy(clm_pf_idata%gn2o_sub_pfs, ierr)

  ! 'accextrn' is accumulative N extract by plant roots in 'PFLOTRAN' within a CLM timestep
    if(clm_pf_idata%accextrn_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%accextrn_sub_pfp, ierr)
    if(clm_pf_idata%accextrn_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%accextrn_sub_clms, ierr)

  ! some tracking variables from PFLOTRAN bgc to obtain reaction flux rates which needed by CLM
    if(clm_pf_idata%acchr_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%acchr_sub_pfp, ierr)
    if(clm_pf_idata%acchr_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%acchr_sub_clms, ierr)

    if(clm_pf_idata%accnmin_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%accnmin_sub_pfp, ierr)
    if(clm_pf_idata%accnmin_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%accnmin_sub_clms, ierr)

    if(clm_pf_idata%accnimm_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%accnimm_sub_pfp, ierr)
    if(clm_pf_idata%accnimm_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%accnimm_sub_clms, ierr)

    if(clm_pf_idata%accngasmin_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%accngasmin_sub_pfp , ierr)
    if(clm_pf_idata%accngasmin_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%accngasmin_sub_clms  , ierr)

    if(clm_pf_idata%accngasnitr_sub_pfp /= 0) &
     call VecDestroy(clm_pf_idata%accngasnitr_sub_pfp, ierr)
    if(clm_pf_idata%accngasnitr_sub_clms /= 0) &
     call VecDestroy(clm_pf_idata%accngasnitr_sub_clms, ierr)

    !----------------------------------------------------------------------------------

  end subroutine CLMPFLOTRANIDataDestroy

end module clm_pflotran_interface_data
