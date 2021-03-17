module elmpf_interface_data
!
! NOTES for convenience:
!        (1) '*_pfp': mpi vecs for PF variables; '_elmp': mpi vecs for ELM variables;
!            '*_pfs': seq vecs for PF variables; '_elms': seq. vecs for ELM variables;
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

  type, public :: elm_pflotran_idata_type


  !------------------- A few global constants --------------------------------------------
  
  ! final time of pflotran run (i.e. stop) in seconds, assuming (re-)starting from zero
  PetscReal :: final_time

  ! numbers of ELM soil layers, grids that are mapped to/from PFLOTRAN (global constants, not local copy)
  PetscInt :: nzelm_mapped
  PetscInt :: nxelm_mapped
  PetscInt :: nyelm_mapped
  PetscReal :: x0elm_global
  PetscReal :: y0elm_global
  PetscReal :: z0elm_global
  PetscReal, pointer :: dxelm_global(:)              ! this is NOT the 3-D vec 'dxsoil' defined below, rather it's the universal x-direction interval (OR, longitudal degree interval from ELM land surf grids) for all gridcells
  PetscReal, pointer :: dyelm_global(:)              ! this is NOT the 3-D vec 'dysoil' defined below, rather it's the universal y-direction interval (OR, longitudal degree interval from ELM land surf grids)
  PetscReal, pointer :: dzelm_global(:)              ! this is NOT the 3-D vec 'dzsoil' defined below, rather it's the universal soil layer thickness (unit: m) for all gridcells

  ! --- Decompose domain in 3-D (only work with structured PF grid currently) ------------
  
  ! processors no.
  PetscInt :: npx, npy, npz
  ! domain nodes no. for each processors
  PetscInt, pointer :: elm_lx(:)   ! array size is 'npx'
  PetscInt, pointer :: elm_ly(:)   ! array size is 'npy'
  PetscInt, pointer :: elm_lz(:)   ! array size is 'npz'

  !--------------------  Mesh property ---------------------------------------------------
  
  ! Number of cells for the 3D subsurface domain
  PetscInt :: nlelm_sub    ! num of local elm cells
  PetscInt :: ngelm_sub    ! num of ghosted elm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_sub     ! num of local pflotran cells
  PetscInt :: ngpf_sub     ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the top/bottom cells of the 3D subsurface domain
  PetscInt :: nlelm_2dtop  ! num of local elm cells
  PetscInt :: ngelm_2dtop  ! num of ghosted elm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2dtop   ! num of local pflotran cells
  PetscInt :: ngpf_2dtop   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  PetscInt :: nlelm_2dbot  ! num of local elm cells
  PetscInt :: ngelm_2dbot  ! num of ghosted elm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2dbot   ! num of local pflotran cells
  PetscInt :: ngpf_2dbot   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the 2D surface domain
  PetscInt :: nlelm_srf    ! num of local elm cells
  PetscInt :: ngelm_srf    ! num of ghosted elm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_srf     ! num of local pflotran cells
  PetscInt :: ngpf_srf     ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! sub-surf/sub-base area (in 2D): the toppest/lowest cells of subsurface domain for BCs
  ! At this moment, assumes both surf/base cells are exactly SAME in area.
  Vec :: area_subsurf_elmp   ! mpi vec
  Vec :: area_subsurf_pfs    ! seq vec
  Vec :: area_subsurf_pfp    ! mpi vec
  Vec :: area_subsurf_elms   ! seq vec

  ! Area of top face of PF cells (in 3D) (note: a PF cell has faces (facets) of top/bottom/east/west/south/north)
  Vec :: area_top_face_elmp  ! mpi vec
  Vec :: area_top_face_pfs   ! seq vec
  Vec :: area_top_face_pfp   ! mpi vec
  Vec :: area_top_face_elms  ! seq vec

  ! z-axis (in 3D), soil depth of center of a soil cell in unit of meters
  Vec :: zsoil_elmp          ! mpi vec
  Vec :: zsoil_pfs           ! seq vec
  Vec :: zsoil_pfp           ! mpi vec
  Vec :: zsoil_elms          ! seq vec

  ! x/y-axis (in 3D), grid center of a soil cell in unit of meters
  Vec :: xsoil_elmp          ! mpi vec
  Vec :: xsoil_pfs           ! seq vec
  Vec :: ysoil_elmp          ! mpi vec
  Vec :: ysoil_pfs           ! seq vec

  ! soil cell inter-nodes coordinates ('vertex' called in PF mesh; 'interface level' called in ELM soil layers)
  Vec :: zisoil_elmp          ! mpi vec
  Vec :: zisoil_pfs           ! seq vec

  ! length/width/thickness of soil cells (in 3D) in unit of meters
  Vec :: dxsoil_elmp          ! mpi vec
  Vec :: dxsoil_pfs           ! seq vec
  Vec :: dysoil_elmp          ! mpi vec
  Vec :: dysoil_pfs           ! seq vec
  Vec :: dzsoil_elmp          ! mpi vec
  Vec :: dzsoil_pfs           ! seq vec
  ! a NOTE here: Given a 3D-cell's 'area_top_face' and 'zsoi' known, 
  !              it's possible to calculate its volume (may be useful ?)

   ! cell IDs (in 3D) (for tesing meshes)
  Vec :: cellid_elmp          ! mpi vec
  Vec :: cellid_pfs           ! seq vec
  Vec :: cellid_pfp           ! mpi vec
  Vec :: cellid_elms          ! seq vec
   ! top layer cell IDs (in 2D) (for tesing meshes)
  Vec :: cellid_2dtop_elmp    ! mpi vec
  Vec :: cellid_2dtop_pfs     ! seq vec
  Vec :: cellid_2dtop_pfp     ! mpi vec
  Vec :: cellid_2dtop_elms    ! seq vec

  !-------------- TH properties ----------------------------------------------------------
  
  PetscBool :: head_based
  PetscReal :: pressure_reference

  ! ELM's hydraulic properties
  Vec :: hksat_x_elmp
  Vec :: hksat_y_elmp
  Vec :: hksat_z_elmp
  Vec :: watsat_elmp
  Vec :: watfc_elmp
  Vec :: bulkdensity_dry_elmp
  Vec :: effporosity_elmp

  Vec :: hksat_x_pfs
  Vec :: hksat_y_pfs
  Vec :: hksat_z_pfs
  Vec :: watsat_pfs
  Vec :: watfc_pfs
  Vec :: bulkdensity_dry_pfs
  Vec :: effporosity_pfs

  ! clapp-Horburger's function parameters - needed in BGC somehow
  Vec :: sucsat_elmp
  Vec :: bsw_elmp
  Vec :: sucsat_pfs
  Vec :: bsw_pfs

  ! PF's hydraulic properties
  Vec :: sr_pcwmax_pfp
  Vec :: pcwmax_pfp
  Vec :: effporosity_pfp
  Vec :: sr_pcwmax_elms
  Vec :: pcwmax_elms
  Vec :: effporosity_elms

  ! ELM's thermal properties
  Vec :: tkwet_elmp     ! unit: W/m/K
  Vec :: tkdry_elmp
  Vec :: tkfrz_elmp
  Vec :: hcvsol_elmp    ! unit: J/m^3-K

  Vec :: tkwet_pfs
  Vec :: tkdry_pfs
  Vec :: tkfrz_pfs
  Vec :: hcvsol_pfs

  ! TH state vecs from ELM (mpi) to PF (seq)
  Vec :: press_ref_elmp                 ! reference pressure head (Pa)
  Vec :: press_ref_pfs

  Vec :: press_elmp                     ! water pressure head (Pa)
  Vec :: soilpsi_elmp                   ! soil matric potential (Pa)
  Vec :: soillsat_elmp                  ! soil liq. water saturation (0 - 1)
  Vec :: soilisat_elmp                  ! soil ice water saturation (0 - 1)
  Vec :: soilliq_elmp                   ! soil liq. water mass (kg/m3 bulk soil)
  Vec :: soilice_elmp                   ! soil ice water mass (kg/m3 bulk soil)
  Vec :: soilt_elmp                     ! soil temperature (degC)
  Vec :: press_pfs
  Vec :: soilpsi_pfs
  Vec :: soillsat_pfs
  Vec :: soilisat_pfs
  Vec :: soilliq_pfs
  Vec :: soilice_pfs
  Vec :: soilt_pfs

  ! TH state vecs from PF (mpi) to ELM (seq)

  Vec :: press_pfp                     ! water pressure head (Pa)
  Vec :: soilpsi_pfp                   ! soil matric potential (Pa)
  Vec :: soillsat_pfp                  ! soil liq. water saturation (0 - 1)
  Vec :: soilisat_pfp                  ! soil ice water saturation (0 - 1)
  Vec :: soilliq_pfp                   ! soil liq. water mass (kg/m3 bulk soil)
  Vec :: soilice_pfp                   ! soil ice water mass (kg/m3 bulk soil)
  Vec :: soilt_pfp                     ! soil temperature (degC)
  Vec :: press_elms
  Vec :: soilpsi_elms
  Vec :: soillsat_elms
  Vec :: soilisat_elms
  Vec :: soilliq_elms
  Vec :: soilice_elms
  Vec :: soilt_elms
 
  !------------------------------- Soil BGC ----------------------------------------------

  !
  ! ----- BGC constants
  !

  ! the following constants are for consistently coverting mass and moles btw ELM-CN and PF bgc
  ! so that mass balance error would not be caused by these
  PetscReal :: N_molecular_weight
  PetscReal :: C_molecular_weight

  ! Soil BGC decomposing pools
  PetscInt :: ndecomp_pools
  logical, pointer :: floating_cn_ratio(:)           ! TRUE => pool has variable C:N ratio
  PetscInt :: ndecomp_elements                       ! no. of elements considered: C, N
  PetscReal, pointer:: decomp_element_ratios(:,:)    ! ratios of elements in decomposing pools (unit: moles)

  !NOTES: The following is what PF bgc right now using for ELM-PFLOTRAN coupling
  ! if need adding or modifying, it's possible and update BOTH here and subroutine 'pflotranModelGetRTspecies'
  ! (Of course, it must be modifying the PF input card and get those variables and relevant reactions in RT).

  ! RT bgc species 'idof' and 'name'
  PetscInt, pointer:: ispec_decomp_c(:)              ! name: pool_name, OR, pool_name // "C"
  PetscInt, pointer:: ispec_decomp_n(:)              ! name: "", OR, pool_name // "N"
  PetscInt, pointer:: ispec_decomp_hr(:)             ! name: pool_name // "CHR"
  PetscInt, pointer:: ispec_decomp_nmin(:)           ! name: pool_name // "NMIN"
  PetscInt, pointer:: ispec_decomp_nimp(:)           ! name: pool_name // "NIMM"
  PetscInt, pointer:: ispec_decomp_nimm(:)           ! name: pool_name // "NIMP"
  character(len=32), pointer :: decomp_pool_name(:)  ! name of pools
  PetscReal, pointer:: ck_decomp_c(:)                ! K: first-order decomposition rate constant in 1/sec
  PetscReal, pointer:: adfactor_ck_c(:)              ! scalar to adjust K based on decomp pools - currently used for speeding-up decomposition
  PetscReal, pointer:: fr_decomp_c(:,:)              ! fractions of downstream pools (receivers - pools excluding donor but adding CO2, note (k,k) is that fraction of CO2 respired)

  PetscInt:: ispec_hrimm
  character(len=32):: name_hrim   = "HRimm"          ! this is for total HR

  PetscInt :: ispec_nmin, ispec_nimp, ispec_nimm
  character(len=32):: name_nmin  = "Nmin"            ! this is for total Nmin
  character(len=32):: name_nimp  = "Nimp"            ! this is for total Nimmp
  character(len=32):: name_nimm  = "Nimm"            ! this is for total Nimm

  PetscInt:: ispec_nh4, ispec_no3, ispec_nh4s, ispec_no3s, ispec_nh4sorb
  character(len=32):: name_nh4     = "NH4+"
  character(len=32):: name_no3     = "NO3-"
  character(len=32):: name_nh4s    = "Ammonium"
  character(len=32):: name_no3s    = "Nitrate"
  character(len=32):: name_nh4sorb = "NH4sorb"

  PetscInt :: ispec_plantndemand, ispec_plantnh4uptake, ispec_plantno3uptake
  character(len=32):: name_plantndemand   = "Plantndemand"
  character(len=32):: name_plantnh4uptake = "Plantnh4uptake"
  character(len=32):: name_plantno3uptake = "Plantno3uptake"

  PetscInt :: ispec_ngasmin, ispec_ngasnitr, ispec_ngasdeni
  character(len=32):: name_ngasmin = "NGASmin"
  character(len=32):: name_ngasnitr= "NGASnitr"
  character(len=32):: name_ngasdeni= "NGASdeni"

  PetscInt :: ispec_co2aq, ispec_n2aq, ispec_n2oaq
  character(len=32):: name_co2aq   = "CO2(aq)"
  character(len=32):: name_n2aq    = "N2(aq)"
  character(len=32):: name_n2oaq   = "N2O(aq)"

  PetscInt :: ispec_co2, ispec_n2, ispec_n2o   ! this is for gases as a holder (immobile)
  character(len=32):: name_co2 = "CO2imm"
  character(len=32):: name_n2o = "N2Oimm"
  character(len=32):: name_n2  = "N2imm"

  !
  ! -----BGC vecs from ELM (mpi, ghosted) to PF (seq, local)
  !
  ! the following can be directly used to drive BGC
  Vec :: t_scalar_elmp                  ! temperature response function value from ELM for decomposition reations
  Vec :: w_scalar_elmp                  ! soil moisture response function value from ELM for decomposition reations
  Vec :: o_scalar_elmp                  ! soil anoxic response function value from ELM for decomposition and/or CH4/NOx processes
  Vec :: depth_scalar_elmp              ! soil depth response function value from ELM for decomposition reations
  Vec :: t_scalar_pfs
  Vec :: w_scalar_pfs
  Vec :: o_scalar_pfs
  Vec :: depth_scalar_pfs
  
  ! A NOTE here: the folllowing decomposing pool vec (1D) is ordered by 'cell' first, then 'species'

  ! initial ground/soil C/N pools from ELM (mpi) to PF (seq)
  Vec :: decomp_cpools_vr_elmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_elmp     ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: smin_no3_vr_elmp          ! (moleN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_elmp          ! (moleN/m3) vertically-resolved soil mineral NH4
  Vec :: smin_nh4sorb_vr_elmp      ! (moleN/m3) vertically-resolved soil absorbed mineral NH4
  Vec :: decomp_cpools_vr_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_pfs      ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: smin_no3_vr_pfs           ! (moleN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_pfs           ! (moleN/m3) vertically-resolved soil mineral NH4
  Vec :: smin_nh4sorb_vr_pfs       ! (moleN/m3) vertically-resolved soil absorbed mineral NH4

  ! time-varying ground/soil C/N rates from ELM (mpi) to PF (seq) (Units: moleC(N)/m3/s)
  Vec :: kscalar_decomp_c_elmp     ! (unitless) a site scalar (c,j) to adjust SOM decomposition rate constants (default: 1.0)
  Vec :: rate_decomp_c_elmp
  Vec :: rate_decomp_n_elmp
  Vec :: rate_smin_no3_elmp
  Vec :: rate_smin_nh4_elmp
  Vec :: rate_plantndemand_elmp
  Vec :: kscalar_decomp_c_pfs
  Vec :: rate_decomp_c_pfs
  Vec :: rate_decomp_n_pfs
  Vec :: rate_smin_no3_pfs
  Vec :: rate_smin_nh4_pfs
  Vec :: rate_plantndemand_pfs

  !
  ! -----BGC vecs from PF (mpi, ghosted) to ELM (seq, local)
  !
  ! BGC state variables
  Vec :: decomp_cpools_vr_pfp
  Vec :: decomp_npools_vr_pfp
  Vec :: smin_no3_vr_pfp
  Vec :: smin_nh4_vr_pfp                ! (moleN/m3) vertically-resolved total soil mineral NH4 (incl. absorbed)
  Vec :: smin_nh4sorb_vr_pfp            ! (moleN/m3) vertically-resolved absorbed NH4-N
  !
  Vec :: decomp_cpools_vr_elms          ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_elms          ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: smin_no3_vr_elms               ! (moleN/m3) vertically-resolved total soil mineral NO3
  Vec :: smin_nh4_vr_elms               ! (moleN/m3) vertically-resolved total soil mineral NH4 (incl. absorbed)
  Vec :: smin_nh4sorb_vr_elms           ! (moleN/m3) vertically-resolved absorbed NH4-N
  !
  ! 'accextrn' is accumulative N extract by plant roots in 'PFLOTRAN' within a ELM timestep
  Vec :: accextrnh4_vr_pfp                ! (moleN/m3) vertically-resolved root extraction N
  Vec :: accextrnh4_vr_elms               ! (moleN/m3) vertically-resolved root extraction N
  Vec :: accextrno3_vr_pfp                ! (moleN/m3) vertically-resolved root extraction N
  Vec :: accextrno3_vr_elms               ! (moleN/m3) vertically-resolved root extraction N

  ! gases in water (aqueous solution of gases)
  ! gases species is accumulative in 'PFLOTRAN', so needs to calculate their fluxes in the ELM-PF interface and reset back to PFLOTRAN
  Vec :: gco2_vr_pfp                   ! (moleC/m3) vertically-resolved soil CO2 C
  Vec :: gco2_vr_elms                  ! (moleC/m3) vertically-resolved soil CO2 C
  Vec :: gco2_vr_elmp                  ! (moleC/m3) vertically-resolved soil CO2 C, after gas emission
  Vec :: gco2_vr_pfs                   ! (moleC/m3) vertically-resolved soil CO2 C, after gas emission

  Vec :: gn2_vr_pfp                    ! (moleN/m3) vertically-resolved N2-N
  Vec :: gn2_vr_elms                   ! (moleN/m3) vertically-resolved N2-N
  Vec :: gn2_vr_elmp                   ! (moleN/m3) vertically-resolved N2-N, after gas emission
  Vec :: gn2_vr_pfs                    ! (moleN/m3) vertically-resolved N2-N, after gas emission

  Vec :: gn2o_vr_pfp                   ! (moleN/m3) vertically-resolved N2O-N
  Vec :: gn2o_vr_elms                  ! (moleN/m3) vertically-resolved N2O-N
  Vec :: gn2o_vr_elmp                  ! (moleN/m3) vertically-resolved N2O-N, after gas emission
  Vec :: gn2o_vr_pfs                   ! (moleN/m3) vertically-resolved N2O-N, after gas emission

  ! some tracking variables from PFLOTRAN bgc to obtain reaction flux rates which needed by ELM
  Vec :: acchr_vr_pfp                  ! (moleC/m3/timestep) vertically-resolved heterotrophic resp. C from individual decomposition
  Vec :: acchr_vr_elms                 ! (moleC/m3/timestep) vertically-resolved heterotrophic resp. C from individual decomposition
  Vec :: acctothr_vr_pfp               ! (moleC/m3/timestep) vertically-resolved heterotrophic resp. C from all decomposition
  Vec :: acctothr_vr_elms              ! (moleC/m3/timestep) vertically-resolved heterotrophic resp. C from all decomposition

  Vec :: accnmin_vr_pfp                ! (moleN/m3/timestep) vertically-resolved N mineralization
  Vec :: accnmin_vr_elms               ! (moleN/m3/timestep) vertically-resolved N mineralization
  Vec :: acctotnmin_vr_pfp             ! (moleN/m3/timestep) vertically-resolved total N mineralization
  Vec :: acctotnmin_vr_elms            ! (moleN/m3/timestep) vertically-resolved total N mineralization

  Vec :: accnimmp_vr_pfp               ! (moleN/m3/timestep) vertically-resolved potential N immoblization
  Vec :: accnimmp_vr_elms              ! (moleN/m3/timestep) vertically-resolved potential N immoblization
  Vec :: acctotnimmp_vr_pfp            ! (moleN/m3/timestep) vertically-resolved potential total N immoblization
  Vec :: acctotnimmp_vr_elms           ! (moleN/m3/timestep) vertically-resolved potential total N immoblization

  Vec :: accnimm_vr_pfp                ! (moleN/m3/timestep) vertically-resolved N immoblization
  Vec :: accnimm_vr_elms               ! (moleN/m3/timestep) vertically-resolved N immoblization
  Vec :: acctotnimm_vr_pfp             ! (moleN/m3/timestep) vertically-resolved total N immoblization
  Vec :: acctotnimm_vr_elms            ! (moleN/m3/timestep) vertically-resolved total N immoblization

  Vec :: accngasmin_vr_pfp              ! (moleN/m3/timestep) vertically-resolved N2O-N from mineralization
  Vec :: accngasmin_vr_elms             ! (moleN/m3/timestep) vertically-resolved N2O-N from mineralization

  Vec :: accngasnitr_vr_pfp             ! (moleN/m3/timestep) vertically-resolved N2O-N from nitrification
  Vec :: accngasnitr_vr_elms            ! (moleN/m3/timestep) vertically-resolved N2O-N from nitrification

  Vec :: accngasdeni_vr_pfp             ! (moleN/m3/timestep) vertically-resolved N2-N from denitrification
  Vec :: accngasdeni_vr_elms            ! (moleN/m3/timestep) vertically-resolved N2-N from denitrification

  ! actual aqeuous N mass flow rate(moleN/m2/sec) at the top (runoff)/bottom (leaching) of 3-D subsurface domain
  ! (+ in, - out)
  Vec :: f_nh4_subsurf_pfp    ! mpi vec
  Vec :: f_nh4_subsurf_elms   ! seq vec
  Vec :: f_nh4_subbase_pfp    ! mpi vec
  Vec :: f_nh4_subbase_elms   ! seq vec
  Vec :: f_no3_subsurf_pfp    ! mpi vec
  Vec :: f_no3_subsurf_elms   ! seq vec
  Vec :: f_no3_subbase_pfp    ! mpi vec
  Vec :: f_no3_subbase_elms   ! seq vec


  !------------------------------- Soil Thermal-Hydrology ----------------------------------------------

  !
  ! -----TH vecs from ELM (mpi, ghosted) to PF (seq, local)
  !

  ! Sink/Source of water (with thermal) for PFLOTRAN's 3D subsurface domain
  Vec :: qflow_elmp   ! mpi vec (H2O): kgH2O/m3/sec
  Vec :: qflow_pfs    ! seq vec
  Vec :: qflowt_elmp  ! mpi vec (H2O)
  Vec :: qflowt_pfs   ! seq vec

  ! Sink/Source of (non-mass) heat flow rate for 3-D subsurface domain
  Vec :: eflow_elmp   ! mpi vec: all non-mass forms of energy src/sink rate (MJ/m3/sec)
  Vec :: eflow_pfs    ! seq vec

  ! BC: water pressure (Pa) on the 2D top/bottom interface of 3-D subsurface domain as boundary conditions from ELM to PF
  Vec :: press_subsurf_elmp    ! mpi vec
  Vec :: press_subbase_elmp    ! mpi vec
  Vec :: press_subsurf_pfs     ! seq vec
  Vec :: press_subbase_pfs     ! seq vec

  Vec :: press_maxponding_elmp ! mpi vec
  Vec :: press_maxponding_pfs  ! seq vec

  ! BC-h: water infiltration/recharge(drainage) (kgH2O/m2/sec) on the 2D top/bottom interface of 3-D subsurface domain as boundary conditions from ELM to PF
  Vec :: qfluxw_subsurf_elmp    ! mpi vec
  Vec :: qfluxw_subbase_elmp    ! mpi vec
  Vec :: qfluxw_subsurf_pfs     ! seq vec
  Vec :: qfluxw_subbase_pfs     ! seq vec

  ! BC-h: soil evaporation (kgH2O/m2/sec) on the 2D top interface of 3-D subsurface domain as boundary conditions from ELM to PF
  Vec :: qfluxev_subsurf_elmp   ! mpi vec
  Vec :: qfluxev_subsurf_pfs    ! seq vec

  ! BC-t: temperature/eflux at the subsurface interface (2D TOP)
  Vec :: gtemp_subsurf_elmp   ! mpi vec: (1) for specifying temperature (thermal state, like enthalpy), or (2) for heat conductance at BC
  Vec :: gtemp_subsurf_pfs    ! seq vec
  Vec :: eflux_subsurf_elmp   ! mpi vec: all forms of energy (MJ/m^2-sec)
  Vec :: eflux_subsurf_pfs    ! seq vec

  Vec :: efluxr_subsurf_elmp  ! mpi vec: radiation form of energy (MJ/m^2-sec)
  Vec :: efluxr_subsurf_pfs   ! seq vec

  Vec :: efluxl_subsurf_elmp  ! mpi vec: latent heat form of energy (MJ/m^2-sec)
  Vec :: efluxl_subsurf_pfs   ! seq vec
  ! (a note here: if specifying 'gtemp_subsurf' above, sensible heat flux at interface may not be needed)

  ! BC-t: temperature/eflux at the subsurface interface (2D BOTTOM)
  Vec :: eflux_subbase_elmp  ! mpi vec
  Vec :: eflux_subbase_pfs   ! seq vec
  Vec :: gtemp_subbase_elmp  ! mpi vec
  Vec :: gtemp_subbase_pfs   ! seq vec

  !
  ! -----TH vecs from PF (mpi, ghosted) to ELM (seq, local)
  !

  ! actual mass water flow rate (kgH2O/m2/sec) through the top/bottom BC (2-D) of 3-D subsurface domain
  ! (+ in, - out)
  Vec :: qevap_subsurf_pfp    ! mpi vec: actual soil evaporation (-)
  Vec :: qevap_subsurf_elms   ! seq vec
  Vec :: qinfl_subsurf_pfp    ! mpi vec: actual infiltration (+)
  Vec :: qinfl_subsurf_elms   ! seq vec
  Vec :: qsurf_subsurf_pfp    ! mpi vec: actual overland flow - potential-actual infiltration or water upwarding (-)
  Vec :: qsurf_subsurf_elms   ! seq vec
  Vec :: qflux_subbase_pfp    ! mpi vec: actual bottom drainage
  Vec :: qflux_subbase_elms   ! seq vec

  Vec :: eflux_subsurf_pfp    ! mpi vec: actual top interface energy flux
  Vec :: eflux_subsurf_elms   ! seq vec
  Vec :: eflux_subbase_pfp    ! mpi vec: actual bottom interface energy flux
  Vec :: eflux_subbase_elms   ! seq vec

  ! net water (with thermal) flow for PFLOTRAN's 3D subsurface domain (NOTE: this is for root water extraction due to transpiration)
  Vec :: qflow_pfp     ! mpi vec (H2O): kgH2O/m3/sec
  Vec :: qflow_elms    ! seq vec
  Vec :: qflowt_pfp    ! mpi vec (H2O)
  Vec :: qflowt_elms   ! seq vec

  ! net heat exchange (non mass) for 3-D subsurface domain
  Vec :: eflow_pfp     ! mpi vec: all non-mass forms of energy exchange rate (MJ/m3/sec)
  Vec :: eflow_elms    ! seq vec



  !---------------------------------------------------------------

  end type elm_pflotran_idata_type

  type(elm_pflotran_idata_type) , public, target , save :: elm_pf_idata
  
  public :: ELMPFLOTRANIDataInit, &
            ELMPFLOTRANIDataCreateVec, &
            ELMPFLOTRANIDataDestroy
  
contains

! ************************************************************************** !

  subroutine ELMPFLOTRANIDataInit()
  ! 
  ! This routine initialized the data transfer type.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  ! Revised by Fengming Yuan, CCSI-ORNL @May-2015
  
    implicit none

    nullify(elm_pf_idata%dxelm_global)
    nullify(elm_pf_idata%dyelm_global)
    nullify(elm_pf_idata%dzelm_global)

    elm_pf_idata%npx = 1    !default 'np' for PF mesh decompose is 1x1x1
    elm_pf_idata%npy = 1
    elm_pf_idata%npz = 1
    nullify(elm_pf_idata%elm_lx)
    nullify(elm_pf_idata%elm_ly)
    nullify(elm_pf_idata%elm_lz)

    elm_pf_idata%nzelm_mapped = 0
    elm_pf_idata%nxelm_mapped = 0
    elm_pf_idata%nyelm_mapped = 0

    elm_pf_idata%x0elm_global = 0
    elm_pf_idata%y0elm_global = 0
    elm_pf_idata%z0elm_global = 0

    elm_pf_idata%nlelm_sub = 0
    elm_pf_idata%ngelm_sub = 0
    elm_pf_idata%nlpf_sub  = 0
    elm_pf_idata%ngpf_sub  = 0

    elm_pf_idata%nlelm_2dtop = 0
    elm_pf_idata%ngelm_2dtop = 0
    elm_pf_idata%nlpf_2dtop  = 0
    elm_pf_idata%ngpf_2dtop  = 0

    elm_pf_idata%nlelm_2dbot = 0
    elm_pf_idata%ngelm_2dbot = 0
    elm_pf_idata%nlpf_2dbot  = 0
    elm_pf_idata%ngpf_2dbot  = 0

    elm_pf_idata%nlelm_srf = 0
    elm_pf_idata%ngelm_srf = 0
    elm_pf_idata%nlpf_srf  = 0
    elm_pf_idata%ngpf_srf  = 0

    !
    elm_pf_idata%zsoil_elmp      = PETSC_NULL_VEC
    elm_pf_idata%zsoil_pfs       = PETSC_NULL_VEC
    elm_pf_idata%zsoil_pfp       = PETSC_NULL_VEC
    elm_pf_idata%zsoil_elms      = PETSC_NULL_VEC
    elm_pf_idata%dxsoil_elmp     = PETSC_NULL_VEC
    elm_pf_idata%dxsoil_pfs      = PETSC_NULL_VEC
    elm_pf_idata%dysoil_elmp     = PETSC_NULL_VEC
    elm_pf_idata%dysoil_pfs      = PETSC_NULL_VEC
    elm_pf_idata%dzsoil_elmp     = PETSC_NULL_VEC
    elm_pf_idata%dzsoil_pfs      = PETSC_NULL_VEC
    elm_pf_idata%xsoil_elmp      = PETSC_NULL_VEC
    elm_pf_idata%xsoil_pfs       = PETSC_NULL_VEC
    elm_pf_idata%ysoil_elmp      = PETSC_NULL_VEC
    elm_pf_idata%ysoil_pfs       = PETSC_NULL_VEC
    elm_pf_idata%zisoil_elmp     = PETSC_NULL_VEC
    elm_pf_idata%zisoil_pfs      = PETSC_NULL_VEC

    elm_pf_idata%area_subsurf_elmp     = PETSC_NULL_VEC
    elm_pf_idata%area_subsurf_pfs      = PETSC_NULL_VEC
    elm_pf_idata%area_subsurf_pfp      = PETSC_NULL_VEC
    elm_pf_idata%area_subsurf_elms     = PETSC_NULL_VEC

    elm_pf_idata%area_top_face_elmp = PETSC_NULL_VEC
    elm_pf_idata%area_top_face_pfs  = PETSC_NULL_VEC
    elm_pf_idata%area_top_face_pfp  = PETSC_NULL_VEC
    elm_pf_idata%area_top_face_elms = PETSC_NULL_VEC

    elm_pf_idata%cellid_elmp     = PETSC_NULL_VEC
    elm_pf_idata%cellid_pfs      = PETSC_NULL_VEC
    elm_pf_idata%cellid_pfp      = PETSC_NULL_VEC
    elm_pf_idata%cellid_elms     = PETSC_NULL_VEC

    elm_pf_idata%cellid_2dtop_elmp     = PETSC_NULL_VEC
    elm_pf_idata%cellid_2dtop_pfs      = PETSC_NULL_VEC
    elm_pf_idata%cellid_2dtop_pfp      = PETSC_NULL_VEC
    elm_pf_idata%cellid_2dtop_elms     = PETSC_NULL_VEC

    !-------------
    elm_pf_idata%head_based = PETSC_TRUE
    elm_pf_idata%pressure_reference = 1.01325d5

    !--------------------------------------------------------------------
    elm_pf_idata%hksat_x_elmp = PETSC_NULL_VEC
    elm_pf_idata%hksat_y_elmp = PETSC_NULL_VEC
    elm_pf_idata%hksat_z_elmp = PETSC_NULL_VEC
    elm_pf_idata%watsat_elmp  = PETSC_NULL_VEC
    elm_pf_idata%watfc_elmp   = PETSC_NULL_VEC
    elm_pf_idata%bulkdensity_dry_elmp = PETSC_NULL_VEC
    elm_pf_idata%effporosity_elmp     = PETSC_NULL_VEC

    elm_pf_idata%tkwet_elmp  = PETSC_NULL_VEC
    elm_pf_idata%tkdry_elmp  = PETSC_NULL_VEC
    elm_pf_idata%tkfrz_elmp  = PETSC_NULL_VEC
    elm_pf_idata%hcvsol_elmp = PETSC_NULL_VEC

    elm_pf_idata%hksat_x_pfs = PETSC_NULL_VEC
    elm_pf_idata%hksat_y_pfs = PETSC_NULL_VEC
    elm_pf_idata%hksat_z_pfs = PETSC_NULL_VEC
    elm_pf_idata%watsat_pfs  = PETSC_NULL_VEC
    elm_pf_idata%watfc_pfs   = PETSC_NULL_VEC
    elm_pf_idata%bulkdensity_dry_pfs = PETSC_NULL_VEC
    elm_pf_idata%effporosity_pfs     = PETSC_NULL_VEC

    elm_pf_idata%tkwet_pfs  = PETSC_NULL_VEC
    elm_pf_idata%tkdry_pfs  = PETSC_NULL_VEC
    elm_pf_idata%tkfrz_pfs  = PETSC_NULL_VEC
    elm_pf_idata%hcvsol_pfs = PETSC_NULL_VEC

    elm_pf_idata%sucsat_elmp = PETSC_NULL_VEC
    elm_pf_idata%bsw_elmp    = PETSC_NULL_VEC
    elm_pf_idata%sucsat_pfs  = PETSC_NULL_VEC
    elm_pf_idata%bsw_pfs     = PETSC_NULL_VEC
   
   !--------------------------------------------------------------------
    elm_pf_idata%sr_pcwmax_pfp   = PETSC_NULL_VEC
    elm_pf_idata%pcwmax_pfp      = PETSC_NULL_VEC
    elm_pf_idata%effporosity_pfp = PETSC_NULL_VEC
    elm_pf_idata%sr_pcwmax_elms  = PETSC_NULL_VEC
    elm_pf_idata%pcwmax_elms     = PETSC_NULL_VEC
    elm_pf_idata%effporosity_elms= PETSC_NULL_VEC

   !--------------------------------------------------------------------
    elm_pf_idata%press_ref_elmp = PETSC_NULL_VEC
    elm_pf_idata%press_ref_pfs  = PETSC_NULL_VEC
    !
    elm_pf_idata%press_elmp    = PETSC_NULL_VEC
    elm_pf_idata%soilpsi_elmp  = PETSC_NULL_VEC
    elm_pf_idata%soillsat_elmp = PETSC_NULL_VEC
    elm_pf_idata%soilisat_elmp = PETSC_NULL_VEC
    elm_pf_idata%soilliq_elmp  = PETSC_NULL_VEC
    elm_pf_idata%soilice_elmp  = PETSC_NULL_VEC
    elm_pf_idata%soilt_elmp    = PETSC_NULL_VEC
    elm_pf_idata%press_pfs      = PETSC_NULL_VEC
    elm_pf_idata%soilpsi_pfs    = PETSC_NULL_VEC
    elm_pf_idata%soillsat_pfs   = PETSC_NULL_VEC
    elm_pf_idata%soilisat_pfs   = PETSC_NULL_VEC
    elm_pf_idata%soilliq_pfs    = PETSC_NULL_VEC
    elm_pf_idata%soilice_pfs    = PETSC_NULL_VEC
    elm_pf_idata%soilt_pfs      = PETSC_NULL_VEC
    !
    elm_pf_idata%press_pfp    = PETSC_NULL_VEC
    elm_pf_idata%soilpsi_pfp  = PETSC_NULL_VEC
    elm_pf_idata%soillsat_pfp = PETSC_NULL_VEC
    elm_pf_idata%soilisat_pfp = PETSC_NULL_VEC
    elm_pf_idata%soilliq_pfp  = PETSC_NULL_VEC
    elm_pf_idata%soilice_pfp  = PETSC_NULL_VEC
    elm_pf_idata%soilt_pfp    = PETSC_NULL_VEC
    elm_pf_idata%press_elms      = PETSC_NULL_VEC
    elm_pf_idata%soilpsi_elms    = PETSC_NULL_VEC
    elm_pf_idata%soillsat_elms   = PETSC_NULL_VEC
    elm_pf_idata%soilisat_elms   = PETSC_NULL_VEC
    elm_pf_idata%soilliq_elms    = PETSC_NULL_VEC
    elm_pf_idata%soilice_elms    = PETSC_NULL_VEC
    elm_pf_idata%soilt_elms      = PETSC_NULL_VEC

    !------------------------------------------------------------------
    elm_pf_idata%N_molecular_weight = 14.0067d0
    elm_pf_idata%C_molecular_weight = 12.0110d0
    nullify(elm_pf_idata%floating_cn_ratio)
    nullify(elm_pf_idata%decomp_element_ratios)
    nullify(elm_pf_idata%ispec_decomp_c)
    nullify(elm_pf_idata%ispec_decomp_n)
    nullify(elm_pf_idata%ispec_decomp_hr)
    nullify(elm_pf_idata%ispec_decomp_nmin)
    nullify(elm_pf_idata%ispec_decomp_nimp)
    nullify(elm_pf_idata%ispec_decomp_nimm)
    nullify(elm_pf_idata%decomp_pool_name)
    nullify(elm_pf_idata%ck_decomp_c)
    nullify(elm_pf_idata%adfactor_ck_c)
    nullify(elm_pf_idata%fr_decomp_c)

    elm_pf_idata%ndecomp_pools    = 0
    elm_pf_idata%ndecomp_elements = 0

    elm_pf_idata%ispec_hrimm   = 0
    elm_pf_idata%ispec_nmin    = 0
    elm_pf_idata%ispec_nimm    = 0
    elm_pf_idata%ispec_nimp    = 0

    elm_pf_idata%ispec_nh4     = 0
    elm_pf_idata%ispec_no3     = 0
    elm_pf_idata%ispec_nh4s    = 0
    elm_pf_idata%ispec_no3s    = 0
    elm_pf_idata%ispec_nh4sorb = 0
    elm_pf_idata%ispec_plantndemand   = 0
    elm_pf_idata%ispec_plantnh4uptake = 0
    elm_pf_idata%ispec_plantno3uptake = 0
    elm_pf_idata%ispec_ngasmin  = 0
    elm_pf_idata%ispec_ngasnitr = 0
    elm_pf_idata%ispec_ngasdeni = 0
    elm_pf_idata%ispec_co2 = 0
    elm_pf_idata%ispec_n2  = 0
    elm_pf_idata%ispec_n2o = 0

   !--------------------------------------------------------------------
    elm_pf_idata%t_scalar_elmp     = PETSC_NULL_VEC
    elm_pf_idata%w_scalar_elmp     = PETSC_NULL_VEC
    elm_pf_idata%o_scalar_elmp     = PETSC_NULL_VEC
    elm_pf_idata%depth_scalar_elmp = PETSC_NULL_VEC
    elm_pf_idata%t_scalar_pfs      = PETSC_NULL_VEC
    elm_pf_idata%w_scalar_pfs      = PETSC_NULL_VEC
    elm_pf_idata%o_scalar_pfs      = PETSC_NULL_VEC
    elm_pf_idata%depth_scalar_pfs  = PETSC_NULL_VEC

    elm_pf_idata%decomp_cpools_vr_elmp = PETSC_NULL_VEC
    elm_pf_idata%decomp_npools_vr_elmp = PETSC_NULL_VEC
    elm_pf_idata%kscalar_decomp_c_elmp = PETSC_NULL_VEC

    elm_pf_idata%smin_no3_vr_elmp      = PETSC_NULL_VEC
    elm_pf_idata%smin_nh4_vr_elmp      = PETSC_NULL_VEC
    elm_pf_idata%smin_nh4sorb_vr_elmp  = PETSC_NULL_VEC

    elm_pf_idata%decomp_cpools_vr_pfs = PETSC_NULL_VEC
    elm_pf_idata%decomp_npools_vr_pfs = PETSC_NULL_VEC
    elm_pf_idata%kscalar_decomp_c_pfs = PETSC_NULL_VEC

    elm_pf_idata%smin_no3_vr_pfs      = PETSC_NULL_VEC
    elm_pf_idata%smin_nh4_vr_pfs      = PETSC_NULL_VEC
    elm_pf_idata%smin_nh4sorb_vr_pfs  = PETSC_NULL_VEC

    elm_pf_idata%rate_decomp_c_elmp         = PETSC_NULL_VEC
    elm_pf_idata%rate_decomp_n_elmp         = PETSC_NULL_VEC
    elm_pf_idata%rate_smin_no3_elmp         = PETSC_NULL_VEC
    elm_pf_idata%rate_smin_nh4_elmp         = PETSC_NULL_VEC
    elm_pf_idata%rate_plantndemand_elmp     = PETSC_NULL_VEC

    elm_pf_idata%rate_decomp_c_pfs         = PETSC_NULL_VEC
    elm_pf_idata%rate_decomp_n_pfs         = PETSC_NULL_VEC
    elm_pf_idata%rate_smin_no3_pfs         = PETSC_NULL_VEC
    elm_pf_idata%rate_smin_nh4_pfs         = PETSC_NULL_VEC
    elm_pf_idata%rate_plantndemand_pfs     = PETSC_NULL_VEC

   !--------------------------------------------------------------------
    ! for C-N states
    elm_pf_idata%decomp_cpools_vr_pfp  = PETSC_NULL_VEC
    elm_pf_idata%decomp_npools_vr_pfp  = PETSC_NULL_VEC
    elm_pf_idata%smin_no3_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%smin_nh4_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%smin_nh4sorb_vr_pfp   = PETSC_NULL_VEC
    elm_pf_idata%decomp_cpools_vr_elms = PETSC_NULL_VEC
    elm_pf_idata%decomp_npools_vr_elms = PETSC_NULL_VEC
    elm_pf_idata%smin_no3_vr_elms      = PETSC_NULL_VEC
    elm_pf_idata%smin_nh4_vr_elms      = PETSC_NULL_VEC
    elm_pf_idata%smin_nh4sorb_vr_elms  = PETSC_NULL_VEC

    ! for root N extraction calculation
    elm_pf_idata%accextrnh4_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%accextrnh4_vr_elms      = PETSC_NULL_VEC
    elm_pf_idata%accextrno3_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%accextrno3_vr_elms      = PETSC_NULL_VEC

    ! for soil hr calculation
    elm_pf_idata%gco2_vr_pfp            = PETSC_NULL_VEC
    elm_pf_idata%gco2_vr_elms           = PETSC_NULL_VEC
    elm_pf_idata%gco2_vr_elmp           = PETSC_NULL_VEC
    elm_pf_idata%gco2_vr_pfs            = PETSC_NULL_VEC

    ! for N2 gas emission calculation
    elm_pf_idata%gn2_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%gn2_vr_elms      = PETSC_NULL_VEC
    elm_pf_idata%gn2_vr_elmp      = PETSC_NULL_VEC
    elm_pf_idata%gn2_vr_pfs       = PETSC_NULL_VEC

    ! for N2O gas emission calculation
    elm_pf_idata%gn2o_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%gn2o_vr_elms      = PETSC_NULL_VEC
    elm_pf_idata%gn2o_vr_elmp      = PETSC_NULL_VEC
    elm_pf_idata%gn2o_vr_pfs       = PETSC_NULL_VEC

    ! for tracking variables in C-N cycle
    elm_pf_idata%acchr_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%acchr_vr_elms      = PETSC_NULL_VEC
    elm_pf_idata%acctothr_vr_pfp    = PETSC_NULL_VEC
    elm_pf_idata%acctothr_vr_elms   = PETSC_NULL_VEC

    elm_pf_idata%accnmin_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%accnmin_vr_elms      = PETSC_NULL_VEC
    elm_pf_idata%acctotnmin_vr_pfp    = PETSC_NULL_VEC
    elm_pf_idata%acctotnmin_vr_elms   = PETSC_NULL_VEC

    elm_pf_idata%accnimmp_vr_pfp      = PETSC_NULL_VEC
    elm_pf_idata%accnimmp_vr_elms     = PETSC_NULL_VEC
    elm_pf_idata%acctotnimmp_vr_pfp   = PETSC_NULL_VEC
    elm_pf_idata%acctotnimmp_vr_elms  = PETSC_NULL_VEC

    elm_pf_idata%accnimm_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%accnimm_vr_elms      = PETSC_NULL_VEC
    elm_pf_idata%acctotnimm_vr_pfp    = PETSC_NULL_VEC
    elm_pf_idata%acctotnimm_vr_elms   = PETSC_NULL_VEC

    elm_pf_idata%accngasmin_vr_pfp        = PETSC_NULL_VEC
    elm_pf_idata%accngasmin_vr_elms       = PETSC_NULL_VEC

    elm_pf_idata%accngasnitr_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%accngasnitr_vr_elms      = PETSC_NULL_VEC

    elm_pf_idata%accngasdeni_vr_pfp       = PETSC_NULL_VEC
    elm_pf_idata%accngasdeni_vr_elms      = PETSC_NULL_VEC

    ! aq. chemical species boundary flux

    elm_pf_idata%f_nh4_subsurf_pfp   = PETSC_NULL_VEC
    elm_pf_idata%f_nh4_subsurf_elms  = PETSC_NULL_VEC
    elm_pf_idata%f_nh4_subbase_pfp   = PETSC_NULL_VEC
    elm_pf_idata%f_nh4_subbase_elms  = PETSC_NULL_VEC
    elm_pf_idata%f_no3_subsurf_pfp   = PETSC_NULL_VEC
    elm_pf_idata%f_no3_subsurf_elms  = PETSC_NULL_VEC
    elm_pf_idata%f_no3_subbase_pfp   = PETSC_NULL_VEC
    elm_pf_idata%f_no3_subbase_elms  = PETSC_NULL_VEC

   !--------------------------------------------------------------------

    elm_pf_idata%qflow_elmp = PETSC_NULL_VEC
    elm_pf_idata%qflow_pfs  = PETSC_NULL_VEC
    elm_pf_idata%qflowt_elmp= PETSC_NULL_VEC
    elm_pf_idata%qflowt_pfs = PETSC_NULL_VEC
    elm_pf_idata%eflow_elmp = PETSC_NULL_VEC
    elm_pf_idata%eflow_pfs  = PETSC_NULL_VEC

    elm_pf_idata%press_maxponding_elmp = PETSC_NULL_VEC
    elm_pf_idata%press_maxponding_pfs  = PETSC_NULL_VEC
    elm_pf_idata%press_subsurf_elmp = PETSC_NULL_VEC
    elm_pf_idata%press_subsurf_pfs  = PETSC_NULL_VEC
    elm_pf_idata%press_subbase_elmp = PETSC_NULL_VEC
    elm_pf_idata%press_subbase_pfs  = PETSC_NULL_VEC

    elm_pf_idata%qfluxw_subsurf_elmp  = PETSC_NULL_VEC
    elm_pf_idata%qfluxw_subsurf_pfs   = PETSC_NULL_VEC
    elm_pf_idata%qfluxev_subsurf_elmp = PETSC_NULL_VEC
    elm_pf_idata%qfluxev_subsurf_pfs  = PETSC_NULL_VEC
    elm_pf_idata%qfluxw_subbase_elmp  = PETSC_NULL_VEC
    elm_pf_idata%qfluxw_subbase_pfs   = PETSC_NULL_VEC

    elm_pf_idata%gtemp_subsurf_elmp = PETSC_NULL_VEC
    elm_pf_idata%gtemp_subsurf_pfs  = PETSC_NULL_VEC
    elm_pf_idata%eflux_subsurf_elmp = PETSC_NULL_VEC
    elm_pf_idata%eflux_subsurf_pfs  = PETSC_NULL_VEC
    elm_pf_idata%efluxr_subsurf_elmp= PETSC_NULL_VEC
    elm_pf_idata%efluxr_subsurf_pfs = PETSC_NULL_VEC
    elm_pf_idata%efluxl_subsurf_elmp= PETSC_NULL_VEC
    elm_pf_idata%efluxl_subsurf_pfs = PETSC_NULL_VEC
    elm_pf_idata%gtemp_subbase_elmp = PETSC_NULL_VEC
    elm_pf_idata%gtemp_subbase_pfs  = PETSC_NULL_VEC
    elm_pf_idata%eflux_subbase_elmp = PETSC_NULL_VEC
    elm_pf_idata%eflux_subbase_pfs  = PETSC_NULL_VEC

    !---------------
    ! water/energy boundary flux
    elm_pf_idata%qevap_subsurf_pfp   = PETSC_NULL_VEC
    elm_pf_idata%qevap_subsurf_elms  = PETSC_NULL_VEC
    elm_pf_idata%qinfl_subsurf_pfp   = PETSC_NULL_VEC
    elm_pf_idata%qinfl_subsurf_elms  = PETSC_NULL_VEC
    elm_pf_idata%qsurf_subsurf_pfp   = PETSC_NULL_VEC
    elm_pf_idata%qsurf_subsurf_elms  = PETSC_NULL_VEC
    elm_pf_idata%qflux_subbase_pfp   = PETSC_NULL_VEC
    elm_pf_idata%qflux_subbase_elms  = PETSC_NULL_VEC
    elm_pf_idata%eflux_subsurf_pfp   = PETSC_NULL_VEC
    elm_pf_idata%eflux_subsurf_elms  = PETSC_NULL_VEC
    elm_pf_idata%eflux_subbase_pfp   = PETSC_NULL_VEC
    elm_pf_idata%eflux_subbase_elms  = PETSC_NULL_VEC

    elm_pf_idata%qflow_pfp   = PETSC_NULL_VEC
    elm_pf_idata%qflow_elms  = PETSC_NULL_VEC
    elm_pf_idata%qflowt_pfp  = PETSC_NULL_VEC
    elm_pf_idata%qflowt_elms = PETSC_NULL_VEC
    elm_pf_idata%eflow_pfp   = PETSC_NULL_VEC
    elm_pf_idata%eflow_elms  = PETSC_NULL_VEC

  end subroutine ELMPFLOTRANIDataInit

! ************************************************************************** !

  subroutine ELMPFLOTRANIDataCreateVec(mycomm)
  ! 
  ! This routine creates PETSc vectors required for data transfer between
  ! ELM and PFLOTRAN.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! Revised by Fengming Yuan, CCSI-ORNL @May-2015
  
    implicit none
    
    PetscErrorCode :: ierr
    PetscMPIInt    :: mycomm, rank

    call MPI_Comm_rank(mycomm,rank, ierr)

    ! The following block of data definition is for THC coupled elm-pflotran (Currently ONLY subsurface or soil)
    !
    !NOTES (fmy): From mpi vecs To seq. vecs for passing data IS in one-way only at this momment
    !             (1) First, here will create 4 sets of 3D/2D vecs.
    !             (2) then, below will copy these vecs to create vecs for other variables.

    ! -------- FOR ELM (mpi) ==> PFLOTRAN (seq)
    ! ELM(mpi)
    call VecCreateMPI(mycomm,elm_pf_idata%nlelm_sub,PETSC_DECIDE,elm_pf_idata%zsoil_elmp,ierr)             ! 3D Subsurface PFLOTRAN
    call VecSet(elm_pf_idata%zsoil_elmp,0.d0,ierr)
    call VecCreateMPI(mycomm,elm_pf_idata%nlelm_2dtop,PETSC_DECIDE,elm_pf_idata%area_subsurf_elmp,ierr)     ! 2D top-cells of 3D Subsurface PFLOTRAN
    call VecSet(elm_pf_idata%area_subsurf_elmp,0.d0,ierr)
    ! PFLOTRAN(seq)
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ngpf_sub,elm_pf_idata%zsoil_pfs,ierr)                   ! 3D Subsurface ELM
    call VecSet(elm_pf_idata%zsoil_pfs,0.d0,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ngpf_2dtop,elm_pf_idata%area_subsurf_pfs,ierr)           ! 2D top-cells of 3D Subsurface ELM
    call VecSet(elm_pf_idata%area_subsurf_pfs,0.d0,ierr)

    ! -------- FOR PFLOTRAN (mpi) ==> ELM (seq)
    ! PFLOTRAN(mpi)
    call VecCreateMPI(mycomm,elm_pf_idata%nlpf_sub,PETSC_DECIDE,elm_pf_idata%zsoil_pfp,ierr)               ! 3D Subsurface PFLOTRAN
    call VecSet(elm_pf_idata%zsoil_pfp,0.d0,ierr)
    call VecCreateMPI(mycomm,elm_pf_idata%nlpf_2dtop,PETSC_DECIDE,elm_pf_idata%area_subsurf_pfp,ierr)     ! 2D top-cells of 3D Subsurface PFLOTRAN
    call VecSet(elm_pf_idata%area_subsurf_pfp,0.d0,ierr)
    ! ELM(seq)
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ngelm_sub,elm_pf_idata%zsoil_elms,ierr)                 ! 3D Subsurface ELM
    call VecSet(elm_pf_idata%zsoil_elms,0.d0,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ngelm_2dtop,elm_pf_idata%area_subsurf_elms,ierr)       ! 2D top-cells of 3D Subsurface ELM
    call VecSet(elm_pf_idata%area_subsurf_elms,0.d0,ierr)



    !
    ! ----------------------------------------------------------------------------
    !
    ! I. CONSTANTS, INITIALS
    !
    !-----------------------------ELM ==> PFLOTRAN
    ! (by copying) Create MPI Vectors for ELM ---------------------------------
    ! soil dimensions
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%xsoil_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%ysoil_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%zisoil_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%dxsoil_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%dysoil_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%dzsoil_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%area_top_face_elmp,ierr)

    ! soil cell ids (3D) / surface cell ids (2D)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%cellid_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%cellid_2dtop_elmp,ierr)

    ! soil physical properties (3D)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%hksat_x_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%hksat_y_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%hksat_z_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%watsat_elmp,ierr)       ! total vwc at saturation (total 'porosity')
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%watfc_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%bulkdensity_dry_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%effporosity_elmp,ierr)     ! this may/may not same as 'bd'/'watsat' above

    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%tkwet_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%tkdry_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%tkfrz_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%hcvsol_elmp,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%sucsat_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%bsw_elmp,ierr)

    ! TH states (3D)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%press_ref_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%press_elmp,ierr)        ! this depends upon 'reference pressure'
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%soilpsi_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%soillsat_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%soilisat_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%soilliq_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%soilice_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%soilt_elmp,ierr)

    ! (by copying) Create Seq. Vectors for PFLOTRAN  ----------------------
    ! soil dimensions
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%xsoil_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%ysoil_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%zisoil_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%dxsoil_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%dysoil_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%dzsoil_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%area_top_face_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%cellid_pfs,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%cellid_2dtop_pfs,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%hksat_x_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%hksat_y_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%hksat_z_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%watsat_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%watfc_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%bulkdensity_dry_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%effporosity_pfs,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%tkwet_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%tkdry_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%tkfrz_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%hcvsol_pfs,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%sucsat_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%bsw_pfs,ierr)

     ! TH states (3D)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%press_ref_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%press_pfs,ierr)        ! this depends upon 'reference pressure'
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%soilpsi_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%soillsat_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%soilisat_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%soilliq_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%soilice_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%soilt_pfs,ierr)

    !-----------------------------PFLOTRAN ==> ELM
    ! (by copying) Create MPI Vectors for PFLOTRAN ------------
    ! soil dimensions
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%area_top_face_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%cellid_pfp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%cellid_2dtop_pfp,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%sr_pcwmax_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%pcwmax_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%effporosity_pfp,ierr)

     ! TH states (3D)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%press_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%soilpsi_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%soillsat_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%soilisat_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%soilliq_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%soilice_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%soilt_pfp,ierr)

    ! (by copying) create Seq. Vectors for ELM ---------
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%area_top_face_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%cellid_elms,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%cellid_2dtop_elms,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%sr_pcwmax_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%pcwmax_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%effporosity_elms,ierr)

     ! TH states (3D)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%press_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%soilpsi_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%soillsat_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%soilisat_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%soilliq_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%soilice_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%soilt_elms,ierr)

    !--------------------------------------------------------------------------------------------------------------------------
    !
    ! II. BGC VARIABLES
    !
    ! BGC state variables: 3D subsurface ELM ---to--- 3D subsurface PFLOTRAN (e.g., initialization or restarting)
    !-----------------------------ELM ==> PFLOTRAN
    ! MPI Vecs for ELM
    call VecCreateMPI(mycomm, elm_pf_idata%ndecomp_pools*elm_pf_idata%nlelm_sub,   &    ! no. of decomp_pools X 3D Subsurface cells
         PETSC_DECIDE,elm_pf_idata%decomp_cpools_vr_elmp,ierr)
    call VecSet(elm_pf_idata%decomp_cpools_vr_elmp,0.d0,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_elmp, elm_pf_idata%decomp_npools_vr_elmp,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%kscalar_decomp_c_elmp,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%t_scalar_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%w_scalar_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%o_scalar_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%depth_scalar_elmp,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%smin_no3_vr_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%smin_nh4_vr_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%smin_nh4sorb_vr_elmp,ierr)

    ! Seq. Vecs for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ndecomp_pools*elm_pf_idata%ngpf_sub, &   ! no. of decomp_pools X 3D Subsurface cells
          elm_pf_idata%decomp_cpools_vr_pfs,ierr)
    call VecSet(elm_pf_idata%decomp_cpools_vr_pfs,0.d0,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_pfs, elm_pf_idata%decomp_npools_vr_pfs,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfs, elm_pf_idata%kscalar_decomp_c_pfs,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%t_scalar_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%w_scalar_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%o_scalar_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%depth_scalar_pfs,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%smin_no3_vr_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%smin_nh4_vr_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%smin_nh4sorb_vr_pfs,ierr)

    ! BGC/TH interface source/sink (rate): 3D subsurface ELM ---to--- 3D subsurface PFLOTRAN
    ! MPI Vecs for ELM
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_elmp,elm_pf_idata%rate_decomp_c_elmp,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_elmp,elm_pf_idata%rate_decomp_n_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%rate_smin_no3_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%rate_smin_nh4_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%rate_plantndemand_elmp,ierr)
    ! Seq. Vecs for PFLOTRAN
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_pfs,elm_pf_idata%rate_decomp_c_pfs,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_pfs,elm_pf_idata%rate_decomp_n_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%rate_smin_no3_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%rate_smin_nh4_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%rate_plantndemand_pfs,ierr)

    ! MPI Vecs for ELM to pass reset aq. conc back to PF
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%gco2_vr_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%gn2_vr_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%gn2o_vr_elmp,ierr)
    ! Seq. Vecs for PFLOTRAN to get reset aq. conc back from ELM
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%gco2_vr_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%gn2_vr_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%gn2o_vr_pfs,ierr)

    !-----------------------------PFLOTRAN ==> ELM
    ! BGC state variables: 3D subsurface PFLOTRAN ---to--- 3D subsurface ELM
    ! MPI Vecs for PFLOTRAN
    call VecCreateMPI(mycomm, elm_pf_idata%ndecomp_pools*elm_pf_idata%nlpf_sub,   &    ! no. of decomp_pools X 3D Subsurface cells
         PETSC_DECIDE,elm_pf_idata%decomp_cpools_vr_pfp,ierr)
    call VecSet(elm_pf_idata%decomp_cpools_vr_pfp,0.d0,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_pfp, elm_pf_idata%decomp_npools_vr_pfp,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%smin_no3_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%smin_nh4_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%smin_nh4sorb_vr_pfp,ierr)

    ! Seq. Vecs for ELM
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ndecomp_pools*elm_pf_idata%ngelm_sub, &   ! no. of decomp_pools X 3D Subsurface cells
          elm_pf_idata%decomp_cpools_vr_elms,ierr)
    call VecSet(elm_pf_idata%decomp_cpools_vr_elms,0.d0,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_elms,elm_pf_idata%decomp_npools_vr_elms,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%smin_no3_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%smin_nh4_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%smin_nh4sorb_vr_elms,ierr)


    ! BGC flux variables: 3D subsurface PFLOTRAN ---to--- 3D subsurface ELM
    ! MPI Vecs for PFLOTRAN
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%gco2_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%gn2_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%gn2o_vr_pfp,ierr)
    !
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%accextrnh4_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%accextrno3_vr_pfp,ierr)
    !
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_pfp,elm_pf_idata%acchr_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_pfp,elm_pf_idata%accnmin_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_pfp,elm_pf_idata%accnimmp_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_pfp,elm_pf_idata%accnimm_vr_pfp,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%acctothr_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%acctotnmin_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%acctotnimmp_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%acctotnimm_vr_pfp,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%accngasmin_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%accngasnitr_vr_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%accngasdeni_vr_pfp,ierr)

    ! Seq. Vecs for ELM
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%gco2_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%gn2_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%gn2o_vr_elms,ierr)
    !
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%accextrnh4_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%accextrno3_vr_elms,ierr)
    !
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_elms,elm_pf_idata%acchr_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_elms,elm_pf_idata%accnmin_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_elms,elm_pf_idata%accnimmp_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%decomp_cpools_vr_elms,elm_pf_idata%accnimm_vr_elms,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%acctothr_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%acctotnmin_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%acctotnimmp_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%acctotnimm_vr_elms,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%accngasmin_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%accngasnitr_vr_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%accngasdeni_vr_elms,ierr)

    ! BC flow variables: 2D faces of subsurface PFLOTRAN ---to--- 2D faces of subsurface ELM
    ! MPI Vecs for PFLOTRAN
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%f_nh4_subsurf_pfp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%f_no3_subsurf_pfp,ierr)
    ! the following assumes SAME bot-cell number as top-cells
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%f_nh4_subbase_pfp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%f_no3_subbase_pfp,ierr)
    ! Seq. Vecs for ELM
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%f_nh4_subsurf_elms,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%f_no3_subsurf_elms,ierr)
    ! the following assumes SAME bot-cell number as top-cells
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%f_nh4_subbase_elms,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%f_no3_subbase_elms,ierr)

    !
    !--------------------------------------------------------------------------------------------------------------------------
    ! THERMAL-HYDROLOGY VARIABLES
    !
    ! --------- For TH data transfer from ELM to PFLOTRAN
    ! (by copying) Create mpi Vectors for ELM  ----------------------
    ! TH Src/Sink (3D)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%qflow_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%qflowt_elmp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elmp,elm_pf_idata%eflow_elmp,ierr)

    ! TH top BC (2D)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%press_subsurf_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%gtemp_subsurf_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%qfluxw_subsurf_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%qfluxev_subsurf_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%eflux_subsurf_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%efluxr_subsurf_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%efluxl_subsurf_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%press_maxponding_elmp,ierr)

    ! TH bottom BC (2D): supposing SAME bottom-cell numbers as top-cells!
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%press_subbase_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%gtemp_subbase_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%qfluxw_subbase_elmp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elmp,elm_pf_idata%eflux_subbase_elmp,ierr)

    ! (by copying) Create Seq. Vectors for PFLOTRAN  ----------------------
    ! TH Src/Sink (3D)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%qflow_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%qflowt_pfs,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfs,elm_pf_idata%eflow_pfs,ierr)

    ! TH top BC (2D)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%press_subsurf_pfs,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%gtemp_subsurf_pfs,ierr)            ! T
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%qfluxw_subsurf_pfs,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%qfluxev_subsurf_pfs,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%eflux_subsurf_pfs,ierr)            ! T
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%efluxr_subsurf_pfs,ierr)            ! T
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%efluxl_subsurf_pfs,ierr)            ! T
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%press_maxponding_pfs,ierr)

    ! TH bottom BC (2D): supposing SAME bot-cell numbers as top-cells!
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%press_subbase_pfs,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%gtemp_subbase_pfs,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%qfluxw_subbase_pfs,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfs,elm_pf_idata%eflux_subbase_pfs,ierr)            ! T

    ! --------- For TH data transfer from PFLOTRAN to ELM
    ! (by copying) Create MPI Vectors for PFLOTRAN  ----------------------
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%qevap_subsurf_pfp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%qinfl_subsurf_pfp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%qsurf_subsurf_pfp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%eflux_subsurf_pfp,ierr)
    ! the following assumes SAME bot-cell number as top-cells
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%qflux_subbase_pfp,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_pfp,elm_pf_idata%eflux_subbase_pfp,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%qflow_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%qflowt_pfp,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_pfp,elm_pf_idata%eflow_pfp,ierr)

    ! (by copying) Create Seq. Vectors for ELM  ----------------------
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%qevap_subsurf_elms,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%qinfl_subsurf_elms,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%qsurf_subsurf_elms,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%eflux_subsurf_elms,ierr)
    ! the following assumes SAME bot-cell number as top-cells
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%qflux_subbase_elms,ierr)
    call VecDuplicate(elm_pf_idata%area_subsurf_elms,elm_pf_idata%eflux_subbase_elms,ierr)

    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%qflow_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%qflowt_elms,ierr)
    call VecDuplicate(elm_pf_idata%zsoil_elms,elm_pf_idata%eflow_elms,ierr)

    !--------------------------------------------------------------------------------------------------------------------------

  end subroutine ELMPFLOTRANIDataCreateVec

! ************************************************************************** !

  subroutine ELMPFLOTRANIDataDestroy()
  ! 
  ! This routine destroys PETSc vectors that were created for data transfer.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  ! Revised by Fengming Yuan, CCSI-ORNL @May-2015
  
    implicit none
    
    PetscErrorCode :: ierr

    if (associated(elm_pf_idata%dxelm_global)) &
    deallocate(elm_pf_idata%dxelm_global)
    if (associated(elm_pf_idata%dyelm_global)) &
    deallocate(elm_pf_idata%dyelm_global)
    if (associated(elm_pf_idata%dzelm_global)) &
    deallocate(elm_pf_idata%dzelm_global)

    !----------------------------------------------------------------------------------

    if(elm_pf_idata%zsoil_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%zsoil_elmp,ierr)
    if(elm_pf_idata%zsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%zsoil_pfs,ierr)
    if(elm_pf_idata%zsoil_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%zsoil_pfp,ierr)
    if(elm_pf_idata%zsoil_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%zsoil_elms,ierr)

    if(elm_pf_idata%xsoil_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%xsoil_elmp,ierr)
    if(elm_pf_idata%xsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%xsoil_pfs,ierr)
    if(elm_pf_idata%ysoil_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%ysoil_elmp,ierr)
    if(elm_pf_idata%ysoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%ysoil_pfs,ierr)
    if(elm_pf_idata%zisoil_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%zisoil_elmp,ierr)
    if(elm_pf_idata%zisoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%zisoil_pfs,ierr)

    if(elm_pf_idata%dxsoil_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%dxsoil_elmp,ierr)
    if(elm_pf_idata%dxsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%dxsoil_pfs,ierr)
    if(elm_pf_idata%dysoil_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%dysoil_elmp,ierr)
    if(elm_pf_idata%dysoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%dysoil_pfs,ierr)
    if(elm_pf_idata%dzsoil_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%dzsoil_elmp,ierr)
    if(elm_pf_idata%dzsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%dzsoil_pfs,ierr)

    if(elm_pf_idata%area_subsurf_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_subsurf_elmp,ierr)
    if(elm_pf_idata%area_subsurf_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_subsurf_pfs,ierr)
    if(elm_pf_idata%area_subsurf_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_subsurf_pfp,ierr)
    if(elm_pf_idata%area_subsurf_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_subsurf_elms,ierr)

    if(elm_pf_idata%area_top_face_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_top_face_elmp,ierr)
    if(elm_pf_idata%area_top_face_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_top_face_pfs,ierr)
    if(elm_pf_idata%area_top_face_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_top_face_pfp,ierr)
    if(elm_pf_idata%area_top_face_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_top_face_elms,ierr)

    !----
    if(elm_pf_idata%cellid_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%cellid_elmp,ierr)
    if(elm_pf_idata%cellid_2dtop_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%cellid_2dtop_elmp,ierr)
    if(elm_pf_idata%cellid_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%cellid_pfs,ierr)
    if(elm_pf_idata%cellid_2dtop_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%cellid_2dtop_pfs,ierr)
    if(elm_pf_idata%hksat_x_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%hksat_x_elmp,ierr)
    if(elm_pf_idata%hksat_y_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%hksat_y_elmp,ierr)
    if(elm_pf_idata%hksat_z_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%hksat_z_elmp,ierr)
    if(elm_pf_idata%sucsat_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%sucsat_elmp,ierr)
    if(elm_pf_idata%watsat_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%watsat_elmp,ierr)
    if(elm_pf_idata%bsw_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%bsw_elmp,ierr)
    if(elm_pf_idata%watfc_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%watfc_elmp,ierr)
    if(elm_pf_idata%bulkdensity_dry_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%bulkdensity_dry_elmp,ierr)

    if(elm_pf_idata%cellid_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%cellid_pfp,ierr)
    if(elm_pf_idata%cellid_2dtop_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%cellid_2dtop_pfp,ierr)
    if(elm_pf_idata%cellid_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%cellid_elms,ierr)
    if(elm_pf_idata%cellid_2dtop_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%cellid_2dtop_elms,ierr)
    if(elm_pf_idata%tkwet_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%tkwet_elmp,ierr)
    if(elm_pf_idata%tkdry_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%tkdry_elmp,ierr)
    if(elm_pf_idata%tkfrz_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%tkfrz_elmp,ierr)
    if(elm_pf_idata%hcvsol_elmp  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%hcvsol_elmp,ierr)

    if(elm_pf_idata%hksat_x_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%hksat_x_pfs,ierr)
    if(elm_pf_idata%hksat_y_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%hksat_y_pfs,ierr)
    if(elm_pf_idata%hksat_z_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%hksat_z_pfs,ierr)
    if(elm_pf_idata%sucsat_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%sucsat_pfs,ierr)
    if(elm_pf_idata%watsat_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%watsat_pfs,ierr)
    if(elm_pf_idata%bsw_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%bsw_pfs,ierr)
    if(elm_pf_idata%watfc_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%watfc_pfs,ierr)
    if(elm_pf_idata%bulkdensity_dry_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%bulkdensity_dry_pfs,ierr)
    if(elm_pf_idata%effporosity_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%effporosity_elmp,ierr)
    if(elm_pf_idata%effporosity_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%effporosity_pfs,ierr)

    !----
    if(elm_pf_idata%tkwet_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%tkwet_pfs,ierr)
    if(elm_pf_idata%tkdry_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%tkdry_pfs,ierr)
    if(elm_pf_idata%tkfrz_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%tkfrz_pfs,ierr)
    if(elm_pf_idata%hcvsol_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%hcvsol_pfs,ierr)

    ! -----
    if(elm_pf_idata%sr_pcwmax_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%sr_pcwmax_pfp,ierr)
    if(elm_pf_idata%pcwmax_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%pcwmax_pfp,ierr)
    if(elm_pf_idata%sr_pcwmax_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%sr_pcwmax_elms,ierr)
    if(elm_pf_idata%pcwmax_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%pcwmax_elms,ierr)
    if(elm_pf_idata%effporosity_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%effporosity_pfp,ierr)
    if(elm_pf_idata%effporosity_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%effporosity_elms,ierr)

    !----
    if(elm_pf_idata%press_ref_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%press_ref_elmp,ierr)
    if(elm_pf_idata%press_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%press_elmp,ierr)
    if(elm_pf_idata%soilpsi_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soilpsi_elmp,ierr)
    if(elm_pf_idata%soillsat_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soillsat_elmp,ierr)
    if(elm_pf_idata%soilisat_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soilisat_elmp,ierr)
    if(elm_pf_idata%soilliq_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soilliq_elmp,ierr)
    if(elm_pf_idata%soilice_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soilice_elmp,ierr)
    if(elm_pf_idata%soilt_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soilt_elmp,ierr)
    
    if(elm_pf_idata%press_ref_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%press_ref_pfs,ierr)
    if(elm_pf_idata%press_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%press_pfs,ierr)
    if(elm_pf_idata%soilpsi_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soilpsi_pfs,ierr)
    if(elm_pf_idata%soillsat_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soillsat_pfs,ierr)
    if(elm_pf_idata%soilisat_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soilisat_pfs,ierr)
    if(elm_pf_idata%soilliq_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soilliq_pfs,ierr)
    if(elm_pf_idata%soilice_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soilice_pfs,ierr)
    if(elm_pf_idata%soilt_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soilt_pfs,ierr)

    !----
    if(elm_pf_idata%press_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%press_pfp,ierr)
    if(elm_pf_idata%soilpsi_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soilpsi_pfp,ierr)
    if(elm_pf_idata%soillsat_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soillsat_pfp,ierr)
    if(elm_pf_idata%soilisat_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soilisat_pfp,ierr)
    if(elm_pf_idata%soilliq_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soilliq_pfp,ierr)
    if(elm_pf_idata%soilice_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soilice_pfp,ierr)
    if(elm_pf_idata%soilt_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%soilt_pfp,ierr)

    if(elm_pf_idata%press_elms /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%press_elms,ierr)
    if(elm_pf_idata%soilpsi_elms /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soilpsi_elms,ierr)
    if(elm_pf_idata%soillsat_elms /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soillsat_elms,ierr)
    if(elm_pf_idata%soilisat_elms /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soilisat_elms,ierr)
    if(elm_pf_idata%soilliq_elms /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soilliq_elms,ierr)
    if(elm_pf_idata%soilice_elms /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soilice_elms,ierr)
    if(elm_pf_idata%soilt_elms /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%soilt_elms,ierr)

    ! -----------------------------------------------------------------------------------------------------------

    if (associated(elm_pf_idata%decomp_element_ratios)) &
    deallocate(elm_pf_idata%decomp_element_ratios)
    if (associated(elm_pf_idata%floating_cn_ratio)) &
    deallocate(elm_pf_idata%floating_cn_ratio)
    if (associated(elm_pf_idata%ck_decomp_c)) &
    deallocate(elm_pf_idata%ck_decomp_c)
    if (associated(elm_pf_idata%adfactor_ck_c)) &
    deallocate(elm_pf_idata%adfactor_ck_c)
    if (associated(elm_pf_idata%fr_decomp_c)) &
    deallocate(elm_pf_idata%fr_decomp_c)

    if (associated(elm_pf_idata%ispec_decomp_c)) &
    deallocate(elm_pf_idata%ispec_decomp_c)
    if (associated(elm_pf_idata%ispec_decomp_n)) &
    deallocate(elm_pf_idata%ispec_decomp_n)
    if (associated(elm_pf_idata%ispec_decomp_hr)) &
    deallocate(elm_pf_idata%ispec_decomp_hr)
    if (associated(elm_pf_idata%ispec_decomp_nmin)) &
    deallocate(elm_pf_idata%ispec_decomp_nmin)
    if (associated(elm_pf_idata%ispec_decomp_nimp)) &
    deallocate(elm_pf_idata%ispec_decomp_nimp)
    if (associated(elm_pf_idata%ispec_decomp_nimm)) &
    deallocate(elm_pf_idata%ispec_decomp_nimm)
    if (associated(elm_pf_idata%decomp_pool_name)) &
    deallocate(elm_pf_idata%decomp_pool_name)

    ! soil C/N pools (initial)
    if(elm_pf_idata%decomp_cpools_vr_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%decomp_cpools_vr_elmp,ierr)
    if(elm_pf_idata%decomp_npools_vr_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%decomp_npools_vr_elmp,ierr)

    if(elm_pf_idata%kscalar_decomp_c_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%kscalar_decomp_c_elmp,ierr)

    if(elm_pf_idata%t_scalar_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%t_scalar_elmp,ierr)
    if(elm_pf_idata%w_scalar_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%w_scalar_elmp,ierr)
    if(elm_pf_idata%o_scalar_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%o_scalar_elmp,ierr)
    if(elm_pf_idata%depth_scalar_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%depth_scalar_elmp,ierr)

    if(elm_pf_idata%smin_no3_vr_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%smin_no3_vr_elmp,ierr)
    if(elm_pf_idata%smin_nh4_vr_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%smin_nh4_vr_elmp,ierr)
    if(elm_pf_idata%smin_nh4sorb_vr_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%smin_nh4sorb_vr_elmp,ierr)

    !
    if(elm_pf_idata%decomp_cpools_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%decomp_cpools_vr_pfs,ierr)
    if(elm_pf_idata%decomp_npools_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%decomp_npools_vr_pfs,ierr)

    if(elm_pf_idata%kscalar_decomp_c_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%kscalar_decomp_c_pfs,ierr)

    if(elm_pf_idata%t_scalar_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%t_scalar_pfs,ierr)
    if(elm_pf_idata%w_scalar_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%w_scalar_pfs,ierr)
    if(elm_pf_idata%o_scalar_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%o_scalar_pfs,ierr)
    if(elm_pf_idata%depth_scalar_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%depth_scalar_pfs,ierr)

    if(elm_pf_idata%smin_no3_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%smin_no3_vr_pfs,ierr)
    if(elm_pf_idata%smin_nh4_vr_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%smin_nh4_vr_pfs,ierr)
    if(elm_pf_idata%smin_nh4sorb_vr_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%smin_nh4sorb_vr_pfs,ierr)


    ! soil C/N fluxes at interface (source/sink)
    if(elm_pf_idata%rate_decomp_c_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%rate_decomp_c_elmp,ierr)
    if(elm_pf_idata%rate_decomp_n_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%rate_decomp_n_elmp,ierr)

    if(elm_pf_idata%rate_plantndemand_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%rate_plantndemand_elmp,ierr)
    if(elm_pf_idata%rate_smin_no3_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%rate_smin_no3_elmp,ierr)
    if(elm_pf_idata%rate_smin_nh4_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%rate_smin_nh4_elmp,ierr)

    if(elm_pf_idata%rate_decomp_c_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%rate_decomp_c_pfs,ierr)
    if(elm_pf_idata%rate_decomp_n_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%rate_decomp_n_pfs,ierr)

    if(elm_pf_idata%rate_plantndemand_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%rate_plantndemand_pfs,ierr)
    if(elm_pf_idata%rate_smin_no3_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%rate_smin_no3_pfs,ierr)
    if(elm_pf_idata%rate_smin_nh4_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%rate_smin_nh4_pfs,ierr)

    !------
    if(elm_pf_idata%decomp_cpools_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%decomp_cpools_vr_pfp,ierr)
    if(elm_pf_idata%decomp_npools_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%decomp_npools_vr_pfp,ierr)
    if(elm_pf_idata%smin_no3_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%smin_no3_vr_pfp,ierr)
    if(elm_pf_idata%smin_nh4_vr_pfp /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%smin_nh4_vr_pfp,ierr)
    if(elm_pf_idata%smin_nh4sorb_vr_pfp /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%smin_nh4sorb_vr_pfp,ierr)

    if(elm_pf_idata%decomp_cpools_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%decomp_cpools_vr_elms,ierr)
    if(elm_pf_idata%decomp_npools_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%decomp_npools_vr_elms,ierr)
    if(elm_pf_idata%smin_no3_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%smin_no3_vr_elms,ierr)
    if(elm_pf_idata%smin_nh4_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%smin_nh4_vr_elms,ierr)
    if(elm_pf_idata%smin_nh4sorb_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%smin_nh4sorb_vr_elms,ierr)

    ! -----------
    if(elm_pf_idata%accextrnh4_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accextrnh4_vr_pfp,ierr)
    if(elm_pf_idata%accextrno3_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accextrno3_vr_pfp,ierr)
    if(elm_pf_idata%accextrnh4_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accextrnh4_vr_elms,ierr)
    if(elm_pf_idata%accextrno3_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accextrno3_vr_elms,ierr)

    if(elm_pf_idata%gco2_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gco2_vr_pfp,ierr)
    if(elm_pf_idata%gco2_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gco2_vr_elms,ierr)
    if(elm_pf_idata%gco2_vr_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gco2_vr_elmp,ierr)
    if(elm_pf_idata%gco2_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gco2_vr_pfs,ierr)

    if(elm_pf_idata%gn2_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gn2_vr_pfp,ierr)
    if(elm_pf_idata%gn2_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gn2_vr_elms,ierr)
    if(elm_pf_idata%gn2_vr_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gn2_vr_elmp,ierr)
    if(elm_pf_idata%gn2_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gn2_vr_pfs,ierr)

    if(elm_pf_idata%gn2o_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gn2o_vr_pfp,ierr)
    if(elm_pf_idata%gn2o_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gn2o_vr_elms,ierr)
    if(elm_pf_idata%gn2o_vr_elmp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gn2o_vr_elmp,ierr)
    if(elm_pf_idata%gn2o_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gn2o_vr_pfs,ierr)

    if(elm_pf_idata%acchr_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%acchr_vr_pfp,ierr)
    if(elm_pf_idata%acchr_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%acchr_vr_elms,ierr)

    if(elm_pf_idata%acctothr_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%acctothr_vr_pfp,ierr)
    if(elm_pf_idata%acctothr_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%acctothr_vr_elms,ierr)

    if(elm_pf_idata%accnmin_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accnmin_vr_pfp,ierr)
    if(elm_pf_idata%accnmin_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accnmin_vr_elms,ierr)

    if(elm_pf_idata%acctotnmin_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%acctotnmin_vr_pfp,ierr)
    if(elm_pf_idata%acctotnmin_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%acctotnmin_vr_elms,ierr)

    if(elm_pf_idata%accnimmp_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accnimmp_vr_pfp,ierr)
    if(elm_pf_idata%accnimmp_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accnimmp_vr_elms,ierr)

    if(elm_pf_idata%acctotnimmp_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%acctotnimmp_vr_pfp,ierr)
    if(elm_pf_idata%acctotnimmp_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%acctotnimmp_vr_elms,ierr)

    if(elm_pf_idata%accnimm_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accnimm_vr_pfp,ierr)
    if(elm_pf_idata%accnimm_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accnimm_vr_elms,ierr)

    if(elm_pf_idata%acctotnimm_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%acctotnimm_vr_pfp,ierr)
    if(elm_pf_idata%acctotnimm_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%acctotnimm_vr_elms,ierr)

    if(elm_pf_idata%accngasmin_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accngasmin_vr_pfp,ierr)
    if(elm_pf_idata%accngasmin_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accngasmin_vr_elms,ierr)

    if(elm_pf_idata%accngasnitr_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accngasnitr_vr_pfp,ierr)
    if(elm_pf_idata%accngasnitr_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accngasnitr_vr_elms,ierr)

    if(elm_pf_idata%accngasdeni_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accngasdeni_vr_pfp,ierr)
    if(elm_pf_idata%accngasdeni_vr_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%accngasdeni_vr_elms,ierr)

    !-------
    if(elm_pf_idata%f_nh4_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%f_nh4_subsurf_pfp,ierr)
    if(elm_pf_idata%f_nh4_subsurf_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%f_nh4_subsurf_elms,ierr)
    if(elm_pf_idata%f_nh4_subbase_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%f_nh4_subbase_pfp,ierr)
    if(elm_pf_idata%f_nh4_subbase_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%f_nh4_subbase_elms,ierr)

    if(elm_pf_idata%f_no3_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%f_no3_subsurf_pfp,ierr)
    if(elm_pf_idata%f_no3_subsurf_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%f_no3_subsurf_elms,ierr)
    if(elm_pf_idata%f_no3_subbase_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%f_no3_subbase_pfp,ierr)
    if(elm_pf_idata%f_no3_subbase_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%f_no3_subbase_elms,ierr)
    !
    ! -----------------------------------------------------------------------------------------------------------
    !-----
    if(elm_pf_idata%qflow_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qflow_elmp,ierr)
    if(elm_pf_idata%qflow_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qflow_pfs,ierr)
    if(elm_pf_idata%qflowt_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qflowt_elmp,ierr)
    if(elm_pf_idata%qflowt_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qflowt_pfs,ierr)
    if(elm_pf_idata%eflow_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflow_elmp,ierr)
    if(elm_pf_idata%eflow_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflow_pfs,ierr)

    !-----
    if(elm_pf_idata%press_maxponding_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%press_maxponding_elmp,ierr)
    if(elm_pf_idata%press_maxponding_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%press_maxponding_pfs,ierr)
    if(elm_pf_idata%press_subsurf_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%press_subsurf_elmp,ierr)
    if(elm_pf_idata%press_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%press_subsurf_pfs,ierr)
    if(elm_pf_idata%press_subbase_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%press_subbase_elmp,ierr)
    if(elm_pf_idata%press_subbase_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%press_subbase_pfs,ierr)
    if(elm_pf_idata%qfluxw_subsurf_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qfluxw_subsurf_elmp,ierr)
    if(elm_pf_idata%qfluxw_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qfluxw_subsurf_pfs,ierr)
    if(elm_pf_idata%qfluxw_subbase_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qfluxw_subbase_elmp,ierr)
    if(elm_pf_idata%qfluxw_subbase_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qfluxw_subbase_pfs,ierr)
    if(elm_pf_idata%qfluxev_subsurf_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qfluxev_subsurf_elmp,ierr)
    if(elm_pf_idata%qfluxev_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qfluxev_subsurf_pfs,ierr)

    if(elm_pf_idata%efluxr_subsurf_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%efluxr_subsurf_elmp,ierr)
    if(elm_pf_idata%efluxr_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%efluxr_subsurf_pfs,ierr)
    if(elm_pf_idata%efluxl_subsurf_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%efluxl_subsurf_elmp,ierr)
    if(elm_pf_idata%efluxl_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%efluxl_subsurf_pfs,ierr)
    if(elm_pf_idata%eflux_subsurf_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflux_subsurf_elmp,ierr)
    if(elm_pf_idata%eflux_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflux_subsurf_pfs,ierr)
    if(elm_pf_idata%gtemp_subsurf_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gtemp_subsurf_elmp,ierr)
    if(elm_pf_idata%gtemp_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gtemp_subsurf_pfs,ierr)
    if(elm_pf_idata%eflux_subbase_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflux_subbase_elmp,ierr)
    if(elm_pf_idata%eflux_subbase_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflux_subbase_pfs,ierr)
    if(elm_pf_idata%gtemp_subbase_elmp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gtemp_subbase_elmp,ierr)
    if(elm_pf_idata%gtemp_subbase_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%gtemp_subbase_pfs,ierr)

    !------------------
    if(elm_pf_idata%qevap_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qevap_subsurf_pfp,ierr)
    if(elm_pf_idata%qevap_subsurf_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qevap_subsurf_elms,ierr)
    if(elm_pf_idata%qinfl_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qinfl_subsurf_pfp,ierr)
    if(elm_pf_idata%qinfl_subsurf_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qinfl_subsurf_elms,ierr)
    if(elm_pf_idata%qsurf_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qsurf_subsurf_pfp,ierr)
    if(elm_pf_idata%qsurf_subsurf_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qsurf_subsurf_elms,ierr)
    if(elm_pf_idata%qflux_subbase_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qflux_subbase_pfp,ierr)
    if(elm_pf_idata%qflux_subbase_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qflux_subbase_elms,ierr)

    if(elm_pf_idata%eflux_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflux_subsurf_pfp,ierr)
    if(elm_pf_idata%eflux_subsurf_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflux_subsurf_elms,ierr)
    if(elm_pf_idata%eflux_subbase_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflux_subbase_pfp,ierr)
    if(elm_pf_idata%eflux_subbase_elms /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflux_subbase_elms,ierr)

    !-----
    if(elm_pf_idata%qflow_pfp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qflow_pfp,ierr)
    if(elm_pf_idata%qflow_elms  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qflow_elms,ierr)
    if(elm_pf_idata%qflowt_pfp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qflowt_pfp,ierr)
    if(elm_pf_idata%qflowt_elms  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%qflowt_elms,ierr)
    if(elm_pf_idata%eflow_pfp  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflow_pfp,ierr)
    if(elm_pf_idata%eflow_elms  /= PETSC_NULL_VEC) &
       call VecDestroy(elm_pf_idata%eflow_elms,ierr)

    ! -----------------------------------------------------------------------------------------------------------

  end subroutine ELMPFLOTRANIDataDestroy

end module elmpf_interface_data
