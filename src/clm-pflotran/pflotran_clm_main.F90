module pflotran_clm_main_module

  use Option_module, only : option_type
  use Simulation_Base_class, only : simulation_base_type
  use Multi_Simulation_module, only : multi_simulation_type
  use Realization_Base_class, only : realization_base_type
  use Mapping_module, only : mapping_type

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petsclog.h"
#include "finclude/petscsysdef.h"
#include "finclude/petscviewer.h"
#include "finclude/petscvec.h"

  type, public :: pflotran_model_type
    class(simulation_base_type),  pointer :: simulation
    type(multi_simulation_type), pointer :: multisimulation
    type(option_type),      pointer :: option
    PetscReal :: pause_time_1
    PetscReal :: pause_time_2
    type(mapping_type),                pointer :: map_clm_sub_to_pf_sub
    type(mapping_type),                pointer :: map_clm_2dtop_to_pf_2dtop
    type(mapping_type),                pointer :: map_pf_sub_to_clm_sub
    type(mapping_type),                pointer :: map_pf_2dtop_to_clm_2dtop

    type(mapping_type),                pointer :: map_clm_2dbot_to_pf_2dbot
    type(mapping_type),                pointer :: map_pf_2dbot_to_clm_2dbot
     
    PetscInt :: nlclm
    PetscInt :: ngclm

  end type pflotran_model_type
  !
  public::pflotranModelCreate,               &
       ! PF running
       pflotranModelStepperRunInit,          &
       pflotranModelUpdateFinalWaypoint,     &
       pflotranModelStepperRunTillPauseTime, &
       pflotranModelSetupRestart,            &
       pflotranModelStepperRunFinalize,      &
       pflotranModelStepperCheckpoint,       &
       pflotranModelDestroy,                 &
       ! Soil properties
       pflotranModelSetSoilProp,              &
       pflotranModelResetSoilPorosityFromCLM, &
       pflotranModelGetSoilPropFromPF,        &
       ! T/H
       pflotranModelUpdateHSourceSink,          &    ! water src/sink (e.g., ET)
       pflotranModelUpdateSubsurfTCond,         &    ! thermal BC
       pflotranModelSetSoilHbcsFromCLM,         &    ! water BC
       pflotranModelSetInternalTHStatesfromCLM, &    ! T/H states from CLM to PFLOTRAN flow mode's field%**
       pflotranModelUpdateTHfromCLM,            &    ! dynamically update TH states from CLM to PF's global vars to drive PFLOTRAN BGC
       pflotranModelGetTemperatureFromPF,       &
       pflotranModelGetSaturationFromPF,        &
       ! BGC
       pflotranModelSetBGCRatesFromCLM,         &
       pflotranModelUpdateAqConcFromCLM,        &
       pflotranModelUpdateAqGasesFromCLM,       &
       pflotranModelSetBgcConcFromCLM,          &
       pflotranModelGetBgcVariablesFromPF,      &
       ! misc.
       pflotranModelGetBCMassBalanceDeltaFromPF

  private :: &
       pflotranModelInsertWaypoint,          &
       pflotranModelDeleteWaypoint,          &
       pflotranModelGetRTspecies

!------------------------------------------------------------
  !NOTES: The following is what PF bgc right now using for CLM-PFLOTRAN coupling
  ! if need adding or modifying, it's possible and update BOTH here and subroutine 'pflotranModelGetRTspecies'
  ! (Of course, it must be modifying the PF input card and get those variables and relevant reactions in RT).

  ! RT bgc species 'idof' and 'name'
  PetscInt, pointer:: ispec_decomps_c(:)
  PetscInt, pointer:: ispec_decomps_n(:)
  character(len=MAXWORDLENGTH), allocatable :: name_decomps(:)          ! appending 'C' or 'N' for real PF species name

  PetscInt:: ispec_nh4, ispec_no3, ispec_nh4sorb
  character(len=MAXWORDLENGTH):: name_nh4     = "NH4+"
  character(len=MAXWORDLENGTH):: name_no3     = "NO3-"
  character(len=MAXWORDLENGTH):: name_nh4sorb = "NH4sorb"

  PetscInt :: ispec_plantndemand, ispec_plantnh4uptake, ispec_plantno3uptake
  character(len=MAXWORDLENGTH):: name_plantndemand   = "Plantndemand"
  character(len=MAXWORDLENGTH):: name_plantnh4uptake = "Plantnh4uptake"
  character(len=MAXWORDLENGTH):: name_plantno3uptake = "Plantno3uptake"

  PetscInt :: ispec_hr, ispec_nmin, ispec_nimmp, ispec_nimm
  character(len=MAXWORDLENGTH):: name_hr   = "HR"
  character(len=MAXWORDLENGTH):: name_nmin = "nmin"
  character(len=MAXWORDLENGTH):: name_nimmp= "nimmp"
  character(len=MAXWORDLENGTH):: name_nimm = "nimm"

  PetscInt :: ispec_ngasmin, ispec_ngasnitr, ispec_ngasdeni
  character(len=MAXWORDLENGTH):: name_ngasmin = "NGASmin"
  character(len=MAXWORDLENGTH):: name_ngasnitr= "NGASnitr"
  character(len=MAXWORDLENGTH):: name_ngasdeni= "NGASdeni"

  PetscInt :: ispec_co2, ispec_n2, ispec_n2o
  character(len=MAXWORDLENGTH):: name_co2 = "CO2imm"
  character(len=MAXWORDLENGTH):: name_n2o = "N2Oimm"
  character(len=MAXWORDLENGTH):: name_n2  = "N2imm"

  PetscReal, parameter :: xeps0_c = 1.0d-20
  PetscReal, parameter :: xeps0_n = 1.0d-21
!------------------------------------------------------------

contains

! ************************************************************************** !

  subroutine pflotranModelCreate(mpicomm, pflotran_prefix, model)
  ! 
  ! Allocates and initializes the pflotranModel object.
  ! It performs the same sequence of commands as done in pflotran.F90
  ! before model integration is performed by the call to StepperRun()
  ! routine
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! 

    use Option_module
    use Simulation_Base_class
    use Multi_Simulation_module
    use Factory_PFLOTRAN_module
    use Factory_Subsurface_module
    use Factory_Hydrogeophysics_module
    use Factory_Surface_module
    use Factory_Surf_Subsurf_module
    use PFLOTRAN_Constants_module
    use Output_Aux_module, only : INSTANTANEOUS_VARS
    use PFLOTRAN_Provenance_module, only : PrintProvenanceToScreen
  
    implicit none

#include "finclude/petscsys.h"

    PetscInt, intent(in) :: mpicomm
    character(len=256), intent(in) :: pflotran_prefix

    type(pflotran_model_type),      pointer :: model

    allocate(model)

    nullify(model%simulation)
    nullify(model%multisimulation)
    nullify(model%option)

    model%option => OptionCreate()
    call OptionInitMPI(model%option, mpicomm)
    call PFLOTRANInitializePrePetsc(model%multisimulation, model%option)

    ! NOTE(bja) 2013-06-25 : external driver must provide an input
    ! prefix string. If the driver wants to use pflotran.in, then it
    ! should explicitly request that with 'pflotran'.
    if (len(trim(pflotran_prefix)) > 1) then
      model%option%input_prefix = trim(pflotran_prefix)
      model%option%input_filename = trim(model%option%input_prefix) // '.in'
      model%option%global_prefix = model%option%input_prefix
    else
      model%option%io_buffer = 'The external driver must provide the ' // &
           'pflotran input file prefix.'
      call printErrMsg(model%option)
    end if

    call OptionInitPetsc(model%option)
    if (model%option%myrank == model%option%io_rank .and. &
        model%option%print_to_screen) then
      call PrintProvenanceToScreen()
    end if
    call PFLOTRANInitializePostPetsc(model%multisimulation, model%option)

    ! NOTE(bja, 2013-07-19) GB's Hack to get communicator correctly
    ! setup on mpich/mac. should be generally ok, but may need an
    ! apple/mpich ifdef if it cause problems elsewhere.
    PETSC_COMM_SELF = MPI_COMM_SELF

    model%pause_time_1 = -1.0d0
    model%pause_time_2 = -1.0d0

    ! FIXME(bja, 2013-07-17) hard code subsurface for this version of CLM-PFLOTRAN
    model%option%simulation_mode = 'SUBSURFACE'

    select case(model%option%simulation_mode)
      case('SUBSURFACE')
         call SubsurfaceInitialize(model%simulation, model%option)
      case default
         model%option%io_buffer = 'Simulation Mode not supported: ' // model%option%simulation_mode
         call printErrMsg(model%option)
    end select

    ! if BGC is on
    if(model%option%ntrandof > 0) then
      call pflotranModelGetRTspecies(model)
    endif

  end subroutine pflotranModelCreate

! ************************************************************************** !

  subroutine pflotranModelStepperRunInit(model)
  ! 
  ! It performs the same execution of commands
  ! that are carried out in StepperRun() before the model integration
  ! begins over the entire simulation time interval
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! 

    type(pflotran_model_type), pointer, intent(inout) :: model

    call model%simulation%InitializeRun()

  end subroutine pflotranModelStepperRunInit

! ************************************************************************** !

  subroutine pflotranModelStepperCheckpoint(model, id_stamp)
  ! 
  ! wrapper around StepperCheckpoint
  ! NOTE(bja, 2013-06-27) : the date stamp from clm is 32 characters
  ! 

    use Option_module

    implicit none

    type(pflotran_model_type), pointer :: model
    character(len=MAXWORDLENGTH), intent(in) :: id_stamp
    PetscViewer :: viewer

    call model%simulation%process_model_coupler_list%Checkpoint(viewer, -1, id_stamp)

  end subroutine pflotranModelStepperCheckpoint

! ************************************************************************** !

  subroutine pflotranModelSetSoilProp(pflotran_model)
  ! 
  ! Converts hydraulic properties from CLM units
  ! into PFLOTRAN units.
  ! #ifdef CLM_PFLOTRAN
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/22/2010
  ! 
  ! (fmy) 4/28/2014: modifications after updating to pflotran-dev

    use Realization_class
    use Discretization_module
    use Patch_module
    use Grid_module
    use Richards_Aux_module
    use TH_Aux_module
    use Field_module
    use Option_module
    use Material_module
    use Material_Aux_class
    use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, PERMEABILITY_XY, &
                               PERMEABILITY_YZ, PERMEABILITY_XZ, &
                               TORTUOSITY, POROSITY

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use clm_pflotran_interface_data
    use Mapping_module
    use Saturation_Function_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(discretization_type), pointer        :: discretization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(richards_auxvar_type), pointer       :: rich_auxvars(:)
    type(richards_auxvar_type), pointer       :: rich_auxvar
    type(th_auxvar_type), pointer             :: th_auxvars(:)
    type(th_auxvar_type), pointer             :: th_auxvar
    type(simulation_base_type), pointer :: simulation
    type(saturation_function_type), pointer :: saturation_function

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal          :: den, vis, grav
    PetscReal, pointer :: porosity_loc_p(:), vol_ovlap_arr(:)
    PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
    PetscReal          :: bc_lambda, bc_alpha, bc_sr

    PetscScalar, pointer :: hksat_x_pf_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_y_pf_loc(:) ! hydraulic conductivity in y-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_z_pf_loc(:) ! hydraulic conductivity in z-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: watsat_pf_loc(:)  ! volumetric soil water at saturation (porosity)
    PetscScalar, pointer :: sucsat_pf_loc(:)  ! minimum soil suction (mm)
    PetscScalar, pointer :: bsw_pf_loc(:)     ! Clapp and Hornberger "b"
    PetscScalar, pointer :: lamda_pf_loc(:)   ! Clapp and Hornberger "1/b"
    PetscScalar, pointer :: alpha_pf_loc(:)   ! Clapp and Hornberger "alpha"
    PetscScalar, pointer :: sr_pf_loc(:)      ! Clapp and Hornberger "sr"
    PetscScalar, pointer :: pcwmax_pf_loc(:)  ! Clapp and Hornberger "pcwmax"

    PetscScalar, pointer :: zsoi_pf_loc(:)    ! soil depth (m)

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav =GRAVITY_CONSTANT       ! [m/S^2]

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: pflotranModelSetSoilProp only works on subsurface simulations."
         call printErrMsg(pflotran_model%option)
    end select

    patch           => realization%patch
    discretization  => realization%discretization
    grid            => patch%grid
    field           => realization%field

    select case(pflotran_model%option%iflowmode)
      case(RICHARDS_MODE)
        rich_auxvars   => patch%aux%Richards%auxvars
      case(TH_MODE)
        th_auxvars   => patch%aux%TH%auxvars
      case default
        !F.-M. Yuan: if no flowmode AND no reactive-transport, then let crash the run
        ! because 'uniform_velocity' with [0, 0, 0] velocity transport doesn't need specific flowmode,
        ! which is used for BGC coupling only (i.e., no TH coupling)
        if(pflotran_model%option%ntrandof.le.0) then
            call printErrMsg(pflotran_model%option, &
               'Current PFLOTRAN mode not supported by pflotranModelSetSoilProp')
        endif
    end select

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%hksat_x_clmp, &
                                    clm_pf_idata%hksat_x_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%hksat_y_clmp, &
                                    clm_pf_idata%hksat_y_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%hksat_z_clmp, &
                                    clm_pf_idata%hksat_z_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sucsat_clmp, &
                                    clm_pf_idata%sucsat_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%bsw_clmp, &
                                    clm_pf_idata%bsw_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%watsat_clmp, &
                                    clm_pf_idata%watsat_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%alpha_clmp, &
                                    clm_pf_idata%alpha_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%lamda_clmp, &
                                    clm_pf_idata%lamda_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sr_clmp, &
                                    clm_pf_idata%sr_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%pcwmax_clmp, &
                                    clm_pf_idata%pcwmax_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%watfc_clmp, &
                                    clm_pf_idata%watfc_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%bulkdensity_dry_clmp, &
                                    clm_pf_idata%bulkdensity_dry_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%zsoi_clmp, &
                                    clm_pf_idata%zsoi_pfs)

    call VecGetArrayF90(clm_pf_idata%hksat_x_pfs, hksat_x_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_y_pfs, hksat_y_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_z_pfs, hksat_z_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%sucsat_pfs,  sucsat_pf_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_pfs,  watsat_pf_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%bsw_pfs,     bsw_pf_loc,     ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%alpha_pfs,  alpha_pf_loc,   ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%lamda_pfs,  lamda_pf_loc,   ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%sr_pfs,     sr_pf_loc,      ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%pcwmax_pfs, pcwmax_pf_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%zsoi_pfs, zsoi_pf_loc, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(field%porosity0, porosity_loc_p, ierr)
    CHKERRQ(ierr)

    if(pflotran_model%option%iflowmode==RICHARDS_MODE .or. &
      pflotran_model%option%iflowmode==TH_MODE) then
      ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
      call VecGetArrayF90(field%perm0_xx,  perm_xx_loc_p,  ierr)
      CHKERRQ(ierr)
      call VecGetArrayF90(field%perm0_yy,  perm_yy_loc_p,  ierr)
      CHKERRQ(ierr)
      call VecGetArrayF90(field%perm0_zz,  perm_zz_loc_p,  ierr)
      CHKERRQ(ierr)
    endif

    do ghosted_id = 1, grid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (ghosted_id < 0 .or. local_id <= 0) cycle
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      !(TODO) need a better way to generate MVG parameters from CLM inputs

      !F.-M. Yuan: (1) the following IS to pass CLM soil hydraulic data into 'saturation_function';
      !            (2) data-passing IS by from 'ghosted_id' to PF's 'local_id'.
      if(pflotran_model%option%nflowdof > 0) then

       saturation_function => patch%  &
         saturation_function_array(patch%sat_func_id(local_id))%ptr

       ! currently VG-Mulaem saturation function type - NOT YET done in CLM (TODO)
       ! if ( saturation_function%saturation_function_itype == VAN_GENUCHTEN .and. &
       !      saturation_function%permeability_function_itype == MUALEM )         then
       !
       !   bc_alpha  = alpha_pf_loc(ghosted_id)
       !   bc_lambda = lamda_pf_loc(ghosted_id)   ! 'm' in VG function ( or, n=1/(1-m))
       !   bc_Sr(1)  = sr_pf_loc(ghosted_id)      ! currently only for liq. water, NOTE that if at 'Sr', pc = inf
       !   bc_pcwmax = pcwmax_pf_loc(ghosted_id)  ! this parameter IS not corresponding with 'Sr'

       ! currently BC-Burdine saturation function type, with specified values to match with Clapp-Hornberger Eq.
        if ( saturation_function%saturation_function_itype == BROOKS_COREY .and. &
             saturation_function%permeability_function_itype == BURDINE )         then
          ! Clapp-Hornberger: soilpsi = sucsat * (-9.81) * (fsattmp)**(-bsw)  ! mm H2O Head --> -pa
          !                   K = Ks*fsattmp**(3+2/bsw)
          !         vs.
          ! BC-Burdine: pc =  (Se**(-1.d0/lambda))/alpha, with Se=(lsat-Sr)/(1-Sr)
          !             relative_perm = Se**power, with power = 3+2/lamda

          bc_alpha  = 1.d0/(9.81d0*sucsat_pf_loc(ghosted_id))
          bc_lambda = 1.d0/bsw_pf_loc(ghosted_id)   !
          bc_sr     = 0.0d0
        else
            call printErrMsg(pflotran_model%option, &
               'Currently ONLY support Brooks_COREY-Burdine saturation function type when coupled with CLM')

        end if

      select case(pflotran_model%option%iflowmode)
        case(RICHARDS_MODE)
          rich_auxvar => rich_auxvars(ghosted_id)
          rich_auxvar%bc_alpha  = bc_alpha
          rich_auxvar%bc_lambda = bc_lambda
          rich_auxvar%bc_sr1    = bc_sr
        case(TH_MODE)
          th_auxvar => th_auxvars(ghosted_id)
          th_auxvar%bc_alpha  = bc_alpha
          th_auxvar%bc_lambda = bc_lambda
          th_auxvar%bc_sr1    = bc_sr
      end select

      endif

      ! hydraulic conductivity => permissivity IS going to 'field%'
      ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
      ! [m^2]          [mm/sec]
      if(pflotran_model%option%iflowmode==RICHARDS_MODE .or. &
         pflotran_model%option%iflowmode==TH_MODE) then
           ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
           perm_xx_loc_p(local_id) = hksat_x_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0
           perm_yy_loc_p(local_id) = hksat_y_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0
           perm_zz_loc_p(local_id) = hksat_z_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0
      endif

      porosity_loc_p(local_id) = watsat_pf_loc(ghosted_id)

#ifdef CLM_PF_DEBUG
      !F.-M. Yuan: the following IS a checking, comparing CLM passed data (watsat):
      !  (turn it on with similar output in clm_pflotran_interfaceMod.F90 and reaction_sandbox_denitrification.F90)
      ! Conclusions: (1) local_id runs from 1 ~ grid%nlmax; and ghosted_id is obtained by 'nL2G' as corrected above;
      !              OR, ghosted_id runs from 1 ~ grid%ngmax; and local_id is obtained by 'nG2L'.
      !              (2) data-passing IS by from 'ghosted_id' to PF-internal (field%)-vec 'local_id';
      write(pflotran_model%option%myrank+200,*) 'checking pflotran-model prior to set soil properties: ', &
        'rank=',pflotran_model%option%myrank, 'ngmax=',grid%ngmax, 'nlmax=',grid%nlmax, &
        'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'pfp_porosity(local_id)=',porosity_loc_p(local_id), &
        'clms_watsat(ghosted_id)=',watsat_pf_loc(ghosted_id), &
        'saturation_function_alpha=', saturation_function%alpha, &
        'saturation_function_lambda=', saturation_function%lambda, &
        'saturation_function_sr=', saturation_function%sr(1), &
        'saturation_function_pcwmax=', saturation_function%pcwmax
#endif

    enddo

    call VecRestoreArrayF90(clm_pf_idata%hksat_x_pfs, hksat_x_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y_pfs, hksat_y_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z_pfs, hksat_z_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat_pfs,  sucsat_pf_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_pfs,  watsat_pf_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw_pfs,     bsw_pf_loc,     ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%alpha_pfs,  alpha_pf_loc,   ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%lamda_pfs,  lamda_pf_loc,   ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%sr_pfs,     sr_pf_loc,      ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%pcwmax_pfs, pcwmax_pf_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%zsoi_pfs, zsoi_pf_loc, ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayF90(field%porosity0, porosity_loc_p, ierr)
    CHKERRQ(ierr)
    if(pflotran_model%option%iflowmode==RICHARDS_MODE .or. &
      pflotran_model%option%iflowmode==TH_MODE) then
      ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
      call VecRestoreArrayF90(field%perm0_xx,  perm_xx_loc_p,  ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayF90(field%perm0_yy,  perm_yy_loc_p,  ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayF90(field%perm0_zz,  perm_zz_loc_p,  ierr)
      CHKERRQ(ierr)
    endif

    ! update ghosted values after resetting soil physical properties from CLM
    call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                               field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,0)
    if (pflotran_model%option%nflowdof > 0) then
        call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_X,0)
        call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Y,0)
        call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Z,0)
    endif

  end subroutine pflotranModelSetSoilProp

! ************************************************************************** !

  subroutine pflotranModelStepperRunTillPauseTime(model, pause_time, dtime, isprintout)
  ! 
  ! It performs the model integration
  ! till the specified pause_time.
  ! NOTE: It is assumed 'pause_time' is in seconds.
  ! NOTE(bja, 2013-07) the strange waypoint insertion of t+30min /
  ! deletion of t is to ensure that we always have a valid waypoint in
  ! front of us, but pflotran does not delete them, so we don't want
  ! to accumulate too large of a list.
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! Revised by Fengming YUAN

    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer :: model
    PetscReal, intent(in) :: pause_time
    PetscReal, intent(in) :: dtime
    PetscBool, intent(in) :: isprintout

    PetscReal :: pause_time1

    if(isprintout) then
      if (model%option%io_rank == model%option%myrank) then
        write(model%option%fid_out, *) '>>>> Inserting waypoint at pause_time (s) = ', pause_time
        write(model%option%fid_out, *) '>>>> for CLM timestep: ', pause_time/dtime
      endif
    endif

    pause_time1 = pause_time + dtime!1800.0d0
    call pflotranModelInsertWaypoint(model, pause_time1)

    call model%simulation%RunToTime(pause_time)

    call pflotranModelDeleteWaypoint(model, pause_time)

  end subroutine pflotranModelStepperRunTillPauseTime

! ************************************************************************** !

  subroutine pflotranModelInsertWaypoint(model, waypoint_time)
  ! 
  ! Inserts a waypoint within the waypoint list
  ! so that the model integration can be paused when that waypoint is
  ! reached
  ! NOTE: It is assumed the 'waypoint_time' is in seconds
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! Revised by Fengming YUAN

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surface_class, only : surface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use Realization_class, only : realization_type
    use Surface_Realization_class, only : surface_realization_type

    use Waypoint_module, only : waypoint_type, WaypointCreate, WaypointInsertInList
    use Units_module, only : UnitsConvertToInternal
    use Option_module, only : printErrMsg

    implicit none

    type(pflotran_model_type), pointer :: model
    type(waypoint_type), pointer       :: waypoint
    type(waypoint_type), pointer       :: waypoint2
    PetscReal                          :: waypoint_time
    character(len=MAXWORDLENGTH)       :: word

    class(realization_type), pointer    :: realization
    class(surface_realization_type), pointer :: surf_realization

    select type (simulation => model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
         nullify(surf_realization)
      class is (surface_simulation_type)
         nullify(realization)
         surf_realization => simulation%surf_realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
         surf_realization => simulation%surf_realization
      class default
         nullify(realization)
         nullify(surf_realization)
         model%option%io_buffer = "pflotranModelInsertWaypoint only " // &
              "works on combinations of surface and subsurface simulations."
         call printErrMsg(model%option)
    end select

    word = 's'
    waypoint => WaypointCreate()
    waypoint%time              = waypoint_time * UnitsConvertToInternal(word, model%option)
    waypoint%update_conditions = PETSC_TRUE
    waypoint%print_output      = PETSC_FALSE
    waypoint%final             = PETSC_FALSE
    waypoint%dt_max            = waypoint_time * UnitsConvertToInternal(word, model%option)!3153600.d0

    if (associated(realization)) then
      call WaypointInsertInList(waypoint, realization%waypoint_list)
    end if

    if (associated(surf_realization)) then
      waypoint2 => WaypointCreate(waypoint)
      ! just in case that subsurface realization NOT there (e.g. only surf_simulation),
      ! the above 'WaypointCreate()' will allocate memory for 'waypoint',
      ! which never have a chance to deallocate and potentially leak memory
      if(associated(waypoint)) deallocate(waypoint)

      call WaypointInsertInList(waypoint2, surf_realization%waypoint_list)
    end if

  end subroutine pflotranModelInsertWaypoint

! ************************************************************************** !

  subroutine pflotranModelDeleteWaypoint(model, waypoint_time)

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surface_class, only : surface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use Realization_class, only : realization_type
    use Surface_Realization_class, only : surface_realization_type
    use Surface_Realization_class, only : surface_realization_type

    use Waypoint_module, only : waypoint_type, WaypointCreate, WaypointDeleteFromList
    use Units_module, only : UnitsConvertToInternal
    use Option_module, only : printErrMsg

    implicit none

    type(pflotran_model_type), pointer :: model
    type(waypoint_type), pointer       :: waypoint
    PetscReal                          :: waypoint_time
    character(len=MAXWORDLENGTH)       :: word

    class(realization_type), pointer    :: realization
    class(surface_realization_type), pointer :: surf_realization

    select type (simulation => model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
         nullify(surf_realization)
      class is (surface_simulation_type)
         nullify(realization)
         surf_realization => simulation%surf_realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
         surf_realization => simulation%surf_realization
      class default
         nullify(realization)
         nullify(surf_realization)
         model%option%io_buffer = "pflotranModelDeleteWaypoint only " // &
              "works on combinations of surface and subsurface simulations."
         call printErrMsg(model%option)
    end select

    word = 's'
    waypoint => WaypointCreate()
    waypoint%time              = waypoint_time * UnitsConvertToInternal(word, model%option)
    waypoint%update_conditions = PETSC_TRUE
    waypoint%dt_max            = 3153600

    if (associated(realization)) then
       call WaypointDeleteFromList(waypoint, realization%waypoint_list)
    end if

    if (associated(surf_realization)) then
       call WaypointDeleteFromList(waypoint, surf_realization%waypoint_list)
    end if

    ! when call 'WaypointCreate()', 'waypoint' will be allocated memory,
    ! and the 'WaypointDeleteFromList''s destroy appears not work( TODO checking)
    ! which causes memory leak continuously - it's very dangerous to system if runs long
    if(associated(waypoint)) deallocate(waypoint)
  end subroutine pflotranModelDeleteWaypoint

! ************************************************************************** !

  subroutine pflotranModelSetupRestart(model, restart_stamp)
  !
  ! pflotranModelSetupRestart()
  ! This checks to see if a restart file stamp was provided by the
  ! driver. If so, we set the restart flag and reconstruct the
  ! restart file name. The actual restart is handled by the standard
  ! pflotran mechanism in TimeStepperInitializeRun()
  ! NOTE: this must be called between pflotranModelCreate() and
  ! pflotranModelStepperRunInit()
  !

    use Option_module
    use String_module

    implicit none

    type(pflotran_model_type), pointer :: model
    character(len=MAXWORDLENGTH) :: restart_stamp

    model%option%io_buffer = 'restart is not implemented in clm-pflotran.' // &
       'AND, pflotran will be initialized from CLM'

    if (.not. StringNull(restart_stamp)) then
       model%option%restart_flag = PETSC_TRUE
       model%option%restart_filename = &
            trim(model%option%global_prefix) // &
            trim(model%option%group_prefix) // &
            '-' // trim(restart_stamp) // '.chk'

       model%option%io_buffer = 'restart file is: ' // &
            trim(model%option%restart_filename)

    end if

    call printWrnMsg(model%option)

  end subroutine pflotranModelSetupRestart

! ************************************************************************** !

  subroutine pflotranModelUpdateHSourceSink(pflotran_model)
  !
  ! Update the source/sink term of hydrology
  !
  ! Author: Gautam Bisht
  ! Date: 11/22/2011
  ! Revised by Fengming YUAN

    use clm_pflotran_interface_data
    use Connection_module
    use Coupler_module
    use Grid_module
    use Mapping_module
    use Option_module
    use Realization_class, only : realization_type
    use Simulation_Base_class, only : simulation_base_type
    use String_module
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type

    implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model

    class(realization_type), pointer          :: subsurf_realization
    type(coupler_type), pointer               :: source_sink
    type(grid_type), pointer                  :: grid
    type(connection_set_type), pointer        :: cur_connection_set
    PetscScalar, pointer                      :: qflx_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn, local_id, ghosted_id
    PetscErrorCode                            :: ierr
    PetscInt                                  :: press_dof

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%qflux_clmp, &
                                    clm_pf_idata%qflux_pfs)

    ! Get pointer to subsurface-realization
    select type (simulation => pflotran_model%simulation)
      class is (surfsubsurface_simulation_type)
         subsurf_realization => simulation%realization
      class is (subsurface_simulation_type)
         subsurf_realization => simulation%realization
      class default
         pflotran_model%option%io_buffer = " Unsupported simulation_type " // &
            " in pflotranModelUpdateSourceSink."
         call printErrMsg(pflotran_model%option)
    end select

    ! Find value of pressure-dof depending on flow mode
    select case (pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        press_dof = RICHARDS_PRESSURE_DOF
      case (TH_MODE)
        press_dof = TH_PRESSURE_DOF
      case default
        pflotran_model%option%io_buffer = 'Unsupported Flow mode'
        call printErrMsg(pflotran_model%option)
    end select

    ! Update the 'clm_et_ss' source/sink term
    call VecGetArrayF90(clm_pf_idata%qflux_pfs,qflx_pf_loc,ierr)
    CHKERRQ(ierr)
    found = PETSC_FALSE

    source_sink => subsurf_realization%patch%source_sink_list%first
    grid        => subsurf_realization%patch%grid

    do
      if (.not.associated(source_sink)) exit

      cur_connection_set => source_sink%connection_set

      ! Find appropriate Source/Sink from the list of Source/Sinks
      if(StringCompare(source_sink%name,'clm_et_ss')) then

        found = PETSC_TRUE
        if (source_sink%flow_condition%rate%itype /= HET_MASS_RATE_SS) then
          call printErrMsg(pflotran_model%option,'clm_et_ss is not of ' // &
                           'HET_MASS_RATE_SS')
        endif

        do iconn = 1, cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)

          source_sink%flow_aux_real_var(press_dof,iconn) = qflx_pf_loc(ghosted_id)

#ifdef CLM_PF_DEBUG
      ! the following checking shows data passing IS from 'ghosted_id' to 'iconn (local_id)' (multiple processors)
      write(pflotran_model%option%myrank+200,*) 'checking H-et ss. -pf_model-UpdateSrcSink:', &
        'rank=',pflotran_model%option%myrank, 'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'iconn=',iconn, 'qflx_pfs_loc(iconn)=',qflx_pf_loc(iconn), &
        'qflx_pfs_loc(ghosted_id)=',qflx_pf_loc(ghosted_id)
#endif

        enddo
      endif

      source_sink => source_sink%next
    enddo
    call VecRestoreArrayF90(clm_pf_idata%qflux_pfs,qflx_pf_loc,ierr)
    CHKERRQ(ierr)

    if(.not.found) &
      call printErrMsg(pflotran_model%option,'clm_et_ss not found in ' // &
                       'source-sink list of subsurface model.')

  end subroutine pflotranModelUpdateHSourceSink


! ************************************************************************** !

  subroutine pflotranModelUpdateSubsurfTCond(pflotran_model)
  ! 
  ! This routine updates subsurface boundary condtions of PFLOTRAN related to
  ! energy equation.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 11/08/2013
  ! 

    use clm_pflotran_interface_data
    use Connection_module
    use Coupler_module
    use Mapping_module
    use Option_module
    use Realization_class, only : realization_type
    use Simulation_Base_class, only : simulation_base_type
    use String_module
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type

    implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model

    class(realization_type), pointer          :: subsurf_realization
    type(coupler_type), pointer               :: boundary_condition
    type(connection_set_type), pointer        :: cur_connection_set
    PetscScalar, pointer                      :: gflux_subsurf_pf_loc(:)
    PetscScalar, pointer                      :: gtemp_subsurf_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn
    PetscErrorCode                            :: ierr

    if (clm_pf_idata%nlpf_2dtop <= 0 .and. clm_pf_idata%ngpf_2dtop <= 0 ) return

    ! Map ground-heat flux from CLM--to--PF grid
    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    pflotran_model%option, &
                                    clm_pf_idata%gflux_subsurf_clmp, &
                                    clm_pf_idata%gflux_subsurf_pfs)

    ! Map ground temperature from CLM--to--PF grid
    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    pflotran_model%option, &
                                    clm_pf_idata%gtemp_subsurf_clmp, &
                                    clm_pf_idata%gtemp_subsurf_pfs)

    ! Get pointer to subsurface-realization
    select type (simulation => pflotran_model%simulation)
      class is (surfsubsurface_simulation_type)
         subsurf_realization => simulation%realization
      class is (subsurface_simulation_type)
         subsurf_realization => simulation%realization
      class default
         pflotran_model%option%io_buffer = " Unsupported simulation_type " // &
            " in pflotranModelUpdateSourceSink."
         call printErrMsg(pflotran_model%option)
    end select

    ! Update the 'clm_gflux_bc' ground heat flux BC term
    call VecGetArrayF90(clm_pf_idata%gflux_subsurf_pfs,gflux_subsurf_pf_loc,ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%gtemp_subsurf_pfs,gtemp_subsurf_pf_loc,ierr)
    CHKERRQ(ierr)
    found = PETSC_FALSE
    boundary_condition => subsurf_realization%patch%boundary_condition_list%first
    do
      if (.not.associated(boundary_condition)) exit

      cur_connection_set => boundary_condition%connection_set

      ! Find appropriate BC from the list of boundary conditions
      if(StringCompare(boundary_condition%name,'clm_gflux_bc')) then

        !if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
        !    /= NEUMANN_BC) then
          !call printErrMsg(pflotran_model%option,'clm_gflux_bc is not of ' // &
          !                 'NEUMANN_BC')
        !endif
        found = PETSC_TRUE

        do iconn = 1, cur_connection_set%num_connections
            if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == NEUMANN_BC) then
                 boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          gflux_subsurf_pf_loc(iconn)

            elseif (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == DIRICHLET_BC) then
                 boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          gtemp_subsurf_pf_loc(iconn)

            end if

        enddo
      endif

      boundary_condition => boundary_condition%next
    enddo
    call VecRestoreArrayF90(clm_pf_idata%gflux_subsurf_pfs,gflux_subsurf_pf_loc,ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%gtemp_subsurf_pfs,gtemp_subsurf_pf_loc,ierr)
    CHKERRQ(ierr)

    if(.not.found) &
      call printErrMsg(pflotran_model%option,'clm_gflux_bc not found in ' // &
                       'boundary-condition list of subsurface model.')

  end subroutine pflotranModelUpdateSubsurfTCond

! ************************************************************************** !

  subroutine pflotranModelGetSaturationFromPF(pflotran_model)
  ! 
  ! Extract soil saturation values simulated by
  ! PFLOTRAN in a PETSc vector.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/22/2011
  ! 
  ! 4/28/2014: fmy - updates

    use Option_module
    use Realization_class
    use Patch_module
    use Grid_module
    use Field_module
    use Global_Aux_module
    use Richards_module
    use Richards_Aux_module
    use TH_module
    use TH_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surface_class, only : surface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
!    use Surface_Realization_class, only : surface_realization_type
    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(global_auxvar_type), pointer         :: global_auxvars(:)
    type(richards_auxvar_type), pointer       :: rich_auxvars(:)
    type(TH_auxvar_type), pointer             :: th_auxvars(:)
    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: soillsat_pf_p(:)
    PetscReal, pointer :: soilisat_pf_p(:)
    PetscReal, pointer :: press_pf_p(:)
    PetscReal, pointer :: soilpsi_pf_p(:)
    PetscReal, pointer :: porosity0_loc_p(:)     ! soil porosity in field%porosity0
    PetscScalar, pointer :: porosity_loc_pfp(:)  ! soil porosity saved in clm-pf-idata

    PetscInt :: i
    PetscReal, pointer :: vec_loc(:)
    PetscViewer :: viewer

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call printErrMsg(pflotran_model%option)
    end select

    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field
    global_auxvars  => patch%aux%Global%auxvars
    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
         call RichardsUpdateAuxVars(realization)
         rich_auxvars => patch%aux%Richards%auxvars
      case (TH_MODE)
         call THUpdateAuxVars(realization)
         th_auxvars => patch%aux%TH%auxvars
      case default
        !F.-M. Yuan: if no flowmode AND no reactive-transport, then let crash the run
        ! because 'uniform_velocity' with [0, 0, 0] velocity transport doesn't need specific flowmode,
        ! which is used for BGC coupling only (i.e., no TH coupling)
        if (pflotran_model%option%ntrandof .le. 0) then
            pflotran_model%option%io_buffer='pflotranModelGetUpdatedStates ' // &
             'implmentation in this mode is not supported!'

            call printErrMsg(pflotran_model%option)
        endif

    end select

    ! Save the saturation/pc/pressure values
    call VecGetArrayF90(clm_pf_idata%soillsat_pfp, soillsat_pf_p, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%press_pfp, press_pf_p, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soilpsi_pfp, soilpsi_pf_p, ierr)
    CHKERRQ(ierr)

    ! save porosity for estimating actual water content from saturation, when needed
    call VecGetArrayF90(field%porosity0,porosity0_loc_p,ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%effporosity_pfp, porosity_loc_pfp, ierr)
    CHKERRQ(ierr)

    do local_id=1, grid%nlmax
      ghosted_id=grid%nL2G(local_id)
      if (ghosted_id <=0 ) cycle

      soillsat_pf_p(local_id)=global_auxvars(ghosted_id)%sat(1)
      press_pf_p(local_id)   =global_auxvars(ghosted_id)%pres(1)

      ! PF's field porosity pass to clm-pf-idata and saved
       porosity_loc_pfp(local_id) = porosity0_loc_p(local_id)

#ifdef CLM_PF_DEBUG
! F.-M. Yuan: the following check proves DATA-passing from PF to CLM MUST BE done by ghosted_id --> local_id
! if passing from 'global_auxvars'
write(pflotran_model%option%myrank+200,*) 'checking pflotran-model 2 (PF->CLM lsat):  ', &
        'local_id=',local_id, 'ghosted_id=',ghosted_id,  &
        'sat_globalvar(ghosted_id)=',global_auxvars(ghosted_id)%sat(1), &
        'idata%sat_pfp(local_id)=',soillsat_pf_p(local_id)
#endif


      if (pflotran_model%option%iflowmode == RICHARDS_MODE) then
        soilpsi_pf_p(local_id) = rich_auxvars(ghosted_id)%pc

      else if (pflotran_model%option%iflowmode == TH_MODE) then
        soilpsi_pf_p(local_id) = th_auxvars(ghosted_id)%pc

      endif
    enddo
    call VecRestoreArrayF90(clm_pf_idata%soillsat_pfp, soillsat_pf_p, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_pfp, press_pf_p, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soilpsi_pfp, soilpsi_pf_p, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(field%porosity0,porosity0_loc_p,ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%effporosity_pfp, porosity_loc_pfp, ierr)
    CHKERRQ(ierr)

    ! mapping to CLM vecs (seq)
    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%soillsat_pfp, &
                                    clm_pf_idata%soillsat_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_pfp, &
                                    clm_pf_idata%press_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%soilpsi_pfp, &
                                    clm_pf_idata%soilpsi_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%effporosity_pfp, &
                                    clm_pf_idata%effporosity_clms)

    if (pflotran_model%option%iflowmode == TH_MODE .and. &
        pflotran_model%option%use_th_freezing) then

      TH_auxvars => patch%aux%TH%auxvars

      call VecGetArrayF90(clm_pf_idata%soilisat_pfp, soilisat_pf_p, ierr)
      CHKERRQ(ierr)
      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (ghosted_id <=0 ) cycle
        soilisat_pf_p(local_id) = TH_auxvars(ghosted_id)%sat_ice
      enddo
      call VecRestoreArrayF90(clm_pf_idata%soilisat_pfp, soilisat_pf_p, ierr)
      CHKERRQ(ierr)

      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                      pflotran_model%option, &
                                      clm_pf_idata%soilisat_pfp, &
                                      clm_pf_idata%soilisat_clms)
    endif

  end subroutine pflotranModelGetSaturationFromPF

! ************************************************************************** !

  subroutine pflotranModelGetTemperatureFromPF(pflotran_model)
  ! 
  ! This routine get updated states evoloved by PFLOTRAN.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 5/14/2013
  ! 

    use Option_module
    use Realization_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use TH_module
    use TH_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surface_class, only : surface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
    use Surface_Realization_class, only : surface_realization_type
    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(global_auxvar_type), pointer         :: global_auxvars(:)
    type(th_auxvar_type), pointer             :: th_auxvars(:)

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: soilt_pf_p(:)

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call printErrMsg(pflotran_model%option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    global_auxvars  => patch%aux%Global%auxvars

    select case(pflotran_model%option%iflowmode)
      case (TH_MODE)
         call THUpdateAuxVars(realization)
         th_auxvars => patch%aux%TH%auxvars
      case default
        !F.-M. Yuan: if no flowmode AND no reactive-transport, then let crash the run
        ! because 'uniform_velocity' with [0, 0, 0] velocity transport doesn't need specific flowmode,
        ! which is used for BGC coupling only (i.e., no TH coupling)
        if (pflotran_model%option%ntrandof .le. 0) then
            pflotran_model%option%io_buffer='pflotranModelGetUpdatedStates ' // &
             'implmentation in this mode is not supported!'
            call printErrMsg(pflotran_model%option)
        endif
    end select

    call VecGetArrayF90(clm_pf_idata%soilt_pfp, soilt_pf_p, ierr)
    CHKERRQ(ierr)
    do local_id=1,grid%nlmax
      ghosted_id = grid%nL2G(ghosted_id)
      if (ghosted_id>0) then
        soilt_pf_p(local_id) = global_auxvars(ghosted_id)%temp
      endif
    enddo
    call VecRestoreArrayF90(clm_pf_idata%soilt_pfp, soilt_pf_p, ierr)
    CHKERRQ(ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%soilt_pfp, &
                                    clm_pf_idata%soilt_clms)

  end subroutine pflotranModelGetTemperatureFromPF

! ************************************************************************** !

  subroutine pflotranModelStepperRunFinalize(model)
  ! 
  ! It performs the same execution of commands
  ! that are carried out in StepperRun() once the model integration is
  ! finished
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! 

    implicit none

    type(pflotran_model_type), pointer :: model

    call model%simulation%FinalizeRun()

  end subroutine pflotranModelStepperRunFinalize

! ************************************************************************** !

  subroutine pflotranModelDestroy(model)
  ! 
  ! Deallocates the pflotranModel object
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! 

    use Factory_PFLOTRAN_module, only : PFLOTRANFinalize
    use Option_module, only : OptionFinalize
    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(pflotran_model_type), pointer :: model

    ! FIXME(bja, 2013-07) none of the mapping information appears to
    ! be cleaned up, so we are leaking memory....

    call model%simulation%FinalizeRun()
    !call model%simulation%Strip()    ! this causes petsc error of seq. fault issue, although doesn't matter.
    if(associated(model%simulation)) deallocate(model%simulation)
    nullify(model%simulation)

    call PFLOTRANFinalize(model%option)
    call OptionFinalize(model%option)

    call CLMPFLOTRANIDataDestroy()

    if (associated(model%map_clm_sub_to_pf_sub)) &
      call MappingDestroy(model%map_clm_sub_to_pf_sub)
    if (associated(model%map_clm_2dtop_to_pf_2dtop)) &
      call MappingDestroy(model%map_clm_2dtop_to_pf_2dtop)
    if (associated(model%map_clm_2dbot_to_pf_2dbot)) &
      call MappingDestroy(model%map_clm_2dbot_to_pf_2dbot)

    if (associated(model%map_pf_sub_to_clm_sub)) &
      call MappingDestroy(model%map_pf_sub_to_clm_sub)
    if (associated(model%map_pf_2dtop_to_clm_2dtop)) &
      call MappingDestroy(model%map_pf_2dtop_to_clm_2dtop)
    if (associated(model%map_pf_2dbot_to_clm_2dbot)) &
      call MappingDestroy(model%map_pf_2dbot_to_clm_2dbot)

    if (associated(model)) deallocate(model)
    nullify(model)

    if (associated(ispec_decomps_c)) deallocate(ispec_decomps_c)
    nullify(ispec_decomps_c)
    if (associated(ispec_decomps_n)) deallocate(ispec_decomps_n)
    nullify(ispec_decomps_n)
    if (associated(name_decomps)) deallocate(name_decomps)
    nullify(name_decomps)

  end subroutine pflotranModelDestroy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
! THE FOLLOWING BLOCKS OF CODES ARE FOR CLM-PFLOTRAN BGC COUPLING
!
!
! ************************************************************************** !
  !
  ! pflotranModelUpdateFinalWaypoint:
  !  Get CLM final timestep and converts it to PFLOTRAN final way point.
  !  And also set an option for turning on/off PF printing
! ************************************************************************** !

  subroutine pflotranModelUpdateFinalWaypoint(model, waypoint_time, isprintout)

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surface_class, only : surface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use Realization_class, only : realization_type, &
      RealizationAddWaypointsToList
    use Surface_Realization_class, only : surface_realization_type, &
      SurfRealizAddWaypointsToList

    use Waypoint_module, only : waypoint_type, WaypointCreate, &
      WaypointDeleteFromList, WaypointInsertInList, &
      WaypointListFillIn, WaypointListRemoveExtraWaypnts
    use Units_module, only : UnitsConvertToInternal
    use Option_module, only : printErrMsg

    implicit none

    type(pflotran_model_type), pointer :: model
    type(waypoint_type), pointer       :: waypoint, waypoint1, waypoint2
    PetscReal, intent(in)              :: waypoint_time
    PetscBool, intent(in)              :: isprintout
    character(len=MAXWORDLENGTH)       :: word

    class(realization_type), pointer    :: realization
    class(surface_realization_type), pointer :: surf_realization

!-------------------------------------------------------------------------
    select type (simulation => model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
         nullify(surf_realization)
      class is (surface_simulation_type)
         nullify(realization)
         surf_realization => simulation%surf_realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
         surf_realization => simulation%surf_realization
      class default
         nullify(realization)
         nullify(surf_realization)
         model%option%io_buffer = "pflotranModelUPdateFinalWaypoint is " // &
              "Not support in this mode."
         call printErrMsg(model%option)
    end select

    ! new final waypoint
    word = 's'
    waypoint1 => WaypointCreate()
    waypoint1%time          = waypoint_time * UnitsConvertToInternal(word, model%option)
    waypoint1%print_output  = PETSC_TRUE
    waypoint1%final         = PETSC_TRUE
    waypoint1%dt_max        = waypoint_time * UnitsConvertToInternal(word, model%option)

    ! update subsurface-realization final waypoint
    if (associated(realization)) then
      ! remove original final waypoint
      waypoint => realization%waypoint_list%first
      do
        if (.not.associated(waypoint)) exit
        if (waypoint%final) then
          waypoint%final = PETSC_FALSE
          exit

        else
           waypoint => waypoint%next
        endif
      enddo

      ! insert new final waypoint
      call WaypointInsertInList(waypoint1, realization%waypoint_list)
      call RealizationAddWaypointsToList(realization)
      call WaypointListFillIn(model%option,realization%waypoint_list)
      call WaypointListRemoveExtraWaypnts(model%option,realization%waypoint_list)

      ! turn off the 'print out' if required from CLM
      if(.not.isprintout) then
        if (model%option%io_rank == model%option%myrank) then
          write(model%option%fid_out, *) 'NOTE: h5 output at input-defined interval ' // &
            'for subsurface flow from PFLOTRAN IS OFF! '
        endif
        waypoint => realization%waypoint_list%first
        do
          if (.not.associated(waypoint)) exit
          waypoint%print_output = PETSC_FALSE
          waypoint => waypoint%next
        enddo
      endif

    endif  !if (associated(realization))

    ! update surface-realization final waypoint
    if (associated(surf_realization)) then
      ! remove original final waypoint
      waypoint => surf_realization%waypoint_list%first
      do
        if (.not.associated(waypoint)) exit
        if (waypoint%final) then
          waypoint%final = PETSC_FALSE
          exit

        else
           waypoint => waypoint%next
        endif
      enddo

      ! insert new final waypoint2 (copied from waypoint1)
      waypoint2 => WaypointCreate(waypoint1)
      ! just in case that subsurface realization NOT there (e.g. only surf_simulation),
      ! the above 'WaypointCreate()' will allocate memory for 'waypoint1',
      ! which never have a chance to deallocate and potentially leak memory
      if(associated(waypoint1)) deallocate(waypoint1)

      call WaypointInsertInList(waypoint2, surf_realization%waypoint_list)
      call SurfRealizAddWaypointsToList(surf_realization)
      call WaypointListFillIn(model%option,surf_realization%waypoint_list)
      call WaypointListRemoveExtraWaypnts(model%option,surf_realization%waypoint_list)

      ! turn off the 'print out' if required from CLM
      if(.not.isprintout) then
        if (model%option%io_rank == model%option%myrank) then
          write(model%option%fid_out, *) 'NOTE: h5 output at input-defined interval ' // &
            'for surface flow from PFLOTRAN IS OFF! '
        endif
        waypoint => surf_realization%waypoint_list%first
        do
          if (.not.associated(waypoint)) exit
          waypoint%print_output = PETSC_FALSE
          waypoint => waypoint%next
        enddo
      endif

    end if !if (associated(surf_realization))

  end subroutine pflotranModelUpdateFinalWaypoint

  ! ************************************************************************** !

  subroutine pflotranModelResetSoilPorosityFromCLM(pflotran_model)
  !
  ! Resetting soil porosity in pflotran's internal vecs due to changes from CLM
  ! Note: this is used to adjust porosity of ice from total, when Thermal mode is NOT used in PFLOTRAN
  ! F.-M. YUAN:  4/28/2014

    use Realization_class
    use Discretization_module
    use Patch_module
    use Grid_module
    use Field_module
    use Option_module
    use Material_module
    use Material_Aux_class
    use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, POROSITY

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(discretization_type), pointer        :: discretization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(simulation_base_type), pointer :: simulation

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: porosity0_loc_p(:)    ! this is from 'field%porosity0'
    PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
    PetscReal          :: unitconv, perm_adj, tempreal

    PetscScalar, pointer :: porosity_pfs_loc(:), porosity_pfp_loc(:)  ! these are from 'clm-pf-idata%'
    PetscScalar, pointer :: hksat_x_pf_loc(:), hksat_y_pf_loc(:), hksat_z_pf_loc(:)
    PetscScalar, pointer :: watsat_pf_loc(:), bsw_pf_loc(:)

    !---------------------------------------------------------------------------------

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: pflotranModelResetSoilPorosity only works on subsurface simulations."
         call printErrMsg(pflotran_model%option)
    end select

    patch           => realization%patch
    discretization  => realization%discretization
    grid            => patch%grid
    field           => realization%field

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%effporosity_clmp, &
                                    clm_pf_idata%effporosity_pfs)
    ! for adjusting porosity
    call VecGetArrayF90(clm_pf_idata%effporosity_pfs,  porosity_pfs_loc,  ierr)
    CHKERRQ(ierr)   !seq. vec (to receive '_clmp' vec)
    call VecGetArrayF90(field%porosity0, porosity0_loc_p, ierr)
    CHKERRQ(ierr)

    ! for adjusting permissivity
    if (pflotran_model%option%nflowdof > 0) then

        unitconv  = 0.001002d0/(998.2d0*GRAVITY_CONSTANT)/1000.d0    ! from hydraulic conductivity (mmH2O/sec) to permissivity (kg/sec)
        perm_adj  = 1.0d0

        call VecGetArrayF90(clm_pf_idata%hksat_x_pfs, hksat_x_pf_loc, ierr)
        CHKERRQ(ierr)
        call VecGetArrayF90(clm_pf_idata%hksat_y_pfs, hksat_y_pf_loc, ierr)
        CHKERRQ(ierr)
        call VecGetArrayF90(clm_pf_idata%hksat_z_pfs, hksat_z_pf_loc, ierr)
        CHKERRQ(ierr)
        call VecGetArrayF90(clm_pf_idata%watsat_pfs,  watsat_pf_loc,  ierr)
        CHKERRQ(ierr)
        call VecGetArrayF90(clm_pf_idata%bsw_pfs,  bsw_pf_loc,  ierr)
        CHKERRQ(ierr)

        call VecGetArrayF90(field%perm0_xx,  perm_xx_loc_p,  ierr)
        CHKERRQ(ierr)
        call VecGetArrayF90(field%perm0_yy,  perm_yy_loc_p,  ierr)
        CHKERRQ(ierr)
        call VecGetArrayF90(field%perm0_zz,  perm_zz_loc_p,  ierr)
        CHKERRQ(ierr)
    endif

    do ghosted_id = 1, grid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (ghosted_id < 0 .or. local_id < 0) cycle
      if (patch%imat(ghosted_id) <= 0) cycle

#ifdef CLM_PF_DEBUG
      !F.-M. Yuan: the following IS a checking, comparing CLM passed data (ice-adjusted porosity):
      ! Conclusions: (1) local_id runs from 1 ~ grid%nlmax; and ghosted_id is obtained by 'nL2G' as corrected above;
      !              OR, ghosted_id runs from 1 ~ grid%ngmax; and local_id is obtained by 'nG2L'.
      !              (2) data-passing IS by from 'ghosted_id' to 'local_id'
      write(pflotran_model%option%myrank+200,*) 'checking pflotran-model prior to resetting porosity:', &
        'rank=',pflotran_model%option%myrank, &
        'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'porosity0(local_id)=',porosity0_loc_p(local_id),'adjporo(ghosted_id)=',porosity_pfs_loc(ghosted_id)
#endif

      porosity0_loc_p(local_id) = porosity_pfs_loc(ghosted_id)

      if (pflotran_model%option%nflowdof > 0) then
           ! Ksat is based on actaul porosity, so when porosity is using the effective one, Ksat should be effective as well
           ! This will prevent large hydraulic conductivity in PFLOTRAN when shrinking pore size
           ! because PFLOTRAN uses pressure (saturation) in its rel. perm calculation.
           tempreal = porosity_pfs_loc(ghosted_id)/watsat_pf_loc(ghosted_id)
           perm_adj = tempreal**(2.0d0*bsw_pf_loc(ghosted_id)+3.0d0)        ! assuming shrunk pore as VWC to estimate K, by Clapp-Hornberger Eq.
           perm_adj = max(0.d0, min(perm_adj*perm_adj, 1.0d0))
           perm_xx_loc_p(local_id) = perm_adj*hksat_x_pf_loc(ghosted_id)*unitconv
           perm_yy_loc_p(local_id) = perm_adj*hksat_y_pf_loc(ghosted_id)*unitconv
           perm_zz_loc_p(local_id) = perm_adj*hksat_z_pf_loc(ghosted_id)*unitconv

      endif

    enddo

    call VecRestoreArrayF90(clm_pf_idata%effporosity_pfs,  porosity_pfs_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(field%porosity0, porosity0_loc_p, ierr)
    CHKERRQ(ierr)
    !
    if (pflotran_model%option%nflowdof > 0) then
        call VecRestoreArrayF90(clm_pf_idata%hksat_x_pfs, hksat_x_pf_loc, ierr)
        CHKERRQ(ierr)
        call VecRestoreArrayF90(clm_pf_idata%hksat_y_pfs, hksat_y_pf_loc, ierr)
        CHKERRQ(ierr)
        call VecRestoreArrayF90(clm_pf_idata%hksat_z_pfs, hksat_z_pf_loc, ierr)
        CHKERRQ(ierr)
        call VecRestoreArrayF90(clm_pf_idata%watsat_pfs,  watsat_pf_loc,  ierr)
        CHKERRQ(ierr)
        call VecRestoreArrayF90(clm_pf_idata%bsw_pfs,  bsw_pf_loc,  ierr)
        CHKERRQ(ierr)

        call VecRestoreArrayF90(field%perm0_xx,  perm_xx_loc_p,  ierr)
        CHKERRQ(ierr)
        call VecRestoreArrayF90(field%perm0_yy,  perm_yy_loc_p,  ierr)
        CHKERRQ(ierr)
        call VecRestoreArrayF90(field%perm0_zz,  perm_zz_loc_p,  ierr)
        CHKERRQ(ierr)
    endif

    !
    call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                               field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,0)

    if (pflotran_model%option%nflowdof > 0) then
        call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_X,0)
        call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Y,0)
        call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Z,0)
    endif

  end subroutine pflotranModelResetSoilPorosityFromCLM


! ************************************************************************** !

  subroutine pflotranModelGetSoilPropFromPF(pflotran_model)
  !
  ! Pass soil physical properties from PFLOTRAN to CLM
  !
  ! Author: Fengming Yuan
  ! Date: 1/30/2014
  !

    use Realization_class
    use Patch_module
    use Grid_module
    use Richards_Aux_module
    use TH_Aux_module
    use Field_module
    use Option_module
    use Saturation_Function_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(simulation_base_type), pointer :: simulation
    type(saturation_function_type), pointer :: saturation_function

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal          :: tempreal

    ! pf internal variables
    PetscReal, pointer :: porosity_loc_p(:)

    ! clm-pf-interface Vecs for PF thermal-hydroloical parameters used in CLM-PFLOTRAN interface
    PetscScalar, pointer :: porosity_loc_pfp(:)  ! soil porosity
    PetscScalar, pointer :: sr_pcwmax_loc_pfp(:) ! soil vwc at max. capillary pressure (note: not 'Sr')
    PetscScalar, pointer :: pcwmax_loc_pfp(:)    ! max. capillary pressure

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: pflotranModelGetSoilProp doesn't support."
         call printErrMsg(pflotran_model%option)
    end select

    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field

    call VecGetArrayF90(field%porosity0,porosity_loc_p,ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(clm_pf_idata%effporosity_pfp, porosity_loc_pfp, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%sr_pcwmax_pfp, sr_pcwmax_loc_pfp, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%pcwmax_pfp, pcwmax_loc_pfp, ierr)
    CHKERRQ(ierr)

    do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (associated(patch%imat)) then
           if (patch%imat(ghosted_id) <= 0) cycle
        endif

        ! PF's porosity
        porosity_loc_pfp(local_id) = porosity_loc_p(local_id)

        saturation_function => patch%    &
            saturation_function_array(patch%sat_func_id(ghosted_id))%ptr

#ifdef CLM_PF_DEBUG
      !F.-M. Yuan: the following IS a checking, comparing CLM passing data (bsw ~ 1/lambda --> PFsat_func --> CLM):
      write(pflotran_model%option%myrank+200,*) &
        'checking pflotran-model prior to get soil properties from PF: ', &
        'rank=',pflotran_model%option%myrank, 'ngmax=',grid%ngmax, 'nlmax=',grid%nlmax, &
        'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'sat_funcid(ghosted_id)=',patch%sat_func_id(ghosted_id), &
        'pfsatfun_alpha=',saturation_function%alpha, &
        'richauxvars_alpha= ',patch%aux%Richards%auxvars(ghosted_id)%bc_alpha, &
        'pfsatfun_lambda=',saturation_function%lambda, &
        'richauxvars_lambda= ',patch%aux%Richards%auxvars(ghosted_id)%bc_lambda, &
        'pfsatfun_sr1=',saturation_function%sr(1), &
        'richauxvars_sr1= ',patch%aux%Richards%auxvars(ghosted_id)%bc_sr1
#endif

        if (pflotran_model%option%iflowmode==RICHARDS_MODE) then
          saturation_function%alpha  = patch%aux%Richards%auxvars(ghosted_id)%bc_alpha
          saturation_function%lambda = patch%aux%Richards%auxvars(ghosted_id)%bc_lambda
          saturation_function%sr(1)  = patch%aux%Richards%auxvars(ghosted_id)%bc_sr1
        elseif (pflotran_model%option%iflowmode==TH_MODE) then
          saturation_function%alpha  = patch%aux%TH%auxvars(ghosted_id)%bc_alpha
          saturation_function%lambda = patch%aux%TH%auxvars(ghosted_id)%bc_lambda
          saturation_function%sr(1)  = patch%aux%TH%auxvars(ghosted_id)%bc_sr1
        endif
       ! needs to re-calculate some extra variables for 'saturation_function', if changed above
       call SatFunctionComputePolynomial(pflotran_model%option,saturation_function)
       call PermFunctionComputePolynomial(pflotran_model%option,saturation_function)

        ! PF's limits on soil matrix potential (Capillary pressure)
        pcwmax_loc_pfp(local_id) = saturation_function%pcwmax

        ! PF's limits on soil water at pcwmax (NOT: not 'Sr', at which PC is nearly 'inf')
        call SaturationFunctionCompute(saturation_function%pcwmax, tempreal, &
                                      saturation_function, &
                                      pflotran_model%option)
        sr_pcwmax_loc_pfp(local_id) = tempreal

    enddo

    call VecRestoreArrayF90(field%porosity0,porosity_loc_p,ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%effporosity_pfp, porosity_loc_pfp, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%sr_pcwmax_pfp, sr_pcwmax_loc_pfp, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%pcwmax_pfp, pcwmax_loc_pfp, ierr)
    CHKERRQ(ierr)

    !
    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%effporosity_pfp, &
                                    clm_pf_idata%effporosity_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sr_pcwmax_pfp, &
                                    clm_pf_idata%sr_pcwmax_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%pcwmax_pfp, &
                                    clm_pf_idata%pcwmax_clms)

    ! reference pressure
    clm_pf_idata%pressure_reference = pflotran_model%option%reference_pressure

  end subroutine pflotranModelGetSoilPropFromPF


  ! ************************************************************************** !
  !> This routine Updates TH drivers (PF global vars) for PFLOTRAN bgc that are from CLM
  !! for testing PFLOTRAN-BGC mode
  !!
  ! ************************************************************************** !
  subroutine pflotranModelUpdateTHfromCLM(pflotran_model, pf_hmode, pf_tmode)

    use Option_module
    use Realization_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surface_class, only : surface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
    use Surface_Realization_class, only : surface_realization_type
    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(global_auxvar_type), pointer         :: global_auxvars(:)
    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: soillsat_pf_loc(:), soilisat_pf_loc(:)
    PetscReal, pointer :: soilt_pf_loc(:)
    PetscReal, pointer :: soilpress_pf_loc(:)

    logical, intent(in):: pf_hmode, pf_tmode
    !---------------------------------------------------------------------
    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface/surfsubsurface simulations."
         call printErrMsg(pflotran_model%option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    global_auxvars  => patch%aux%Global%auxvars

    ! Save the liq saturation values from CLM to PFLOTRAN, if needed
    if (.not.pf_hmode) then
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%soillsat_clmp, &
                                    clm_pf_idata%soillsat_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%soilpsi_clmp, &
                                    clm_pf_idata%soilpsi_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_clmp, &
                                    clm_pf_idata%press_pfs)

        call VecGetArrayF90(clm_pf_idata%soillsat_pfs, soillsat_pf_loc, ierr)
        CHKERRQ(ierr)
        call VecGetArrayF90(clm_pf_idata%press_pfs, soilpress_pf_loc, ierr)
        CHKERRQ(ierr)

        do ghosted_id=1, grid%ngmax
          local_id=grid%nG2L(ghosted_id)
          if (ghosted_id<=0 .or. local_id <=0) cycle
          !if (patch%imat(local_id) <= 0) cycle  !(TODO) imat IS 0 for some cells when decomposing domain in X and Y directions.

#ifdef CLM_PF_DEBUG
       ! F.-M. Yuan: the following check proves DATA-passing from CLM to PF MUST BE done by ghosted_id --> ghosted_id
       ! if passing to 'global_auxvars'
          write(pflotran_model%option%myrank+200,*) 'checking pflotran-model 1 (CLM->PF lsat): ', &
              'local_id=',local_id, 'ghosted_id=',ghosted_id, &
              'sat_globalvars(ghosted_id)=',global_auxvars(ghosted_id)%sat(1), &
              'sat_pfs(ghosted_id)=',soillsat_pf_loc(ghosted_id)
#endif
          global_auxvars(ghosted_id)%sat(1)=soillsat_pf_loc(ghosted_id)
          global_auxvars(ghosted_id)%pres(1)=soilpress_pf_loc(ghosted_id)
        enddo

        call VecRestoreArrayF90(clm_pf_idata%soillsat_pfs, soillsat_pf_loc, ierr)
        CHKERRQ(ierr)
        call VecRestoreArrayF90(clm_pf_idata%press_pfs, soilpress_pf_loc, ierr)
        CHKERRQ(ierr)
    endif

    ! Save soil temperature values from CLM to PFLOTRAN, if needed
    if (.not.pf_tmode) then
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%soilt_clmp, &
                                    clm_pf_idata%soilt_pfs)
        call VecGetArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
        CHKERRQ(ierr)

        do ghosted_id=1, grid%ngmax
            local_id=grid%nG2L(ghosted_id)
            if (ghosted_id<=0 .or. local_id<=0) cycle
            !if (patch%imat(local_id) <= 0) cycle  !(TODO) imat IS 0 for some cells when decomposing domain in X and Y directions.

            global_auxvars(ghosted_id)%temp=soilt_pf_loc(ghosted_id)
        enddo
        call VecRestoreArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
        CHKERRQ(ierr)
    endif

  end subroutine pflotranModelUpdateTHfromCLM

  ! ************************************************************************** !
  !
  ! pflotranModelSetInternalTHStatesfromCLM: Set initial TH States from CLM
  !
  ! Note: This subroutine directly set initial soil temperature and saturation from CLM
  !       It's needed because of uniform initialization of TH states in PFLOTRAN, which
  !       are from the input card.
  ! (This is different from the 'pflotranModelUpdateTHfromCLM', which pass TH from CLM to
  !   pflotran's global variables and will not affect the internal vec of TH mode).

  ! author: Fengming YUAN
  ! date: 9/23/2013
  ! ************************************************************************** !
subroutine pflotranModelSetInternalTHStatesfromCLM(pflotran_model)

    use Realization_class
    use Patch_module
    use Grid_module
    use Field_module
    use Discretization_module
    use TH_Aux_module
    use TH_module
    use Richards_module
    use Richards_Aux_module
    use Option_module

    use clm_pflotran_interface_data

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(simulation_base_type), pointer :: simulation

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id, istart, iend, vecsize
    PetscReal, pointer :: xx_loc_p(:)

    PetscScalar, pointer :: soilt_pf_loc(:)      ! temperature [oC]
    PetscScalar, pointer :: soilpress_pf_loc(:)  ! water pressure (Pa)

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = 'ERROR: pflotranModelSetInitialTStatesfromCLM ' // &
               'only works on subsurface/surfsubsurface simulations.'
         call printErrMsg(pflotran_model%option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field

    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_clmp, &
                                    clm_pf_idata%press_pfs)

      case (TH_MODE)
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%soilt_clmp, &
                                    clm_pf_idata%soilt_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_clmp, &
                                    clm_pf_idata%press_pfs)
      case default
        if(pflotran_model%option%ntrandof.le.0) then
            pflotran_model%option%io_buffer='pflotranModelSetInitialTHStatesfromCLM ' // &
              'not implmented for this mode.'
            call printErrMsg(pflotran_model%option)
        endif
    end select

    call VecGetArrayF90(field%flow_xx, xx_loc_p, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%press_pfs, soilpress_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
    CHKERRQ(ierr)

    do local_id = 1, grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif

       iend = local_id*pflotran_model%option%nflowdof
       istart = iend-pflotran_model%option%nflowdof+1

       xx_loc_p(istart)  = soilpress_pf_loc(ghosted_id)
       if (pflotran_model%option%iflowmode .eq. TH_MODE)  then
            xx_loc_p(istart+1)= soilt_pf_loc(ghosted_id)
       end if
    enddo

    call VecRestoreArrayF90(field%flow_xx, xx_loc_p, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_pfs, soilpress_pf_loc, ierr)
    CHKERRQ(ierr)

    call DiscretizationGlobalToLocal(realization%discretization, field%flow_xx, &
         field%flow_xx_loc, NFLOWDOF)
    call VecCopy(field%flow_xx, field%flow_yy, ierr)
    CHKERRQ(ierr)

    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        call RichardsUpdateAuxVars(realization)
      case (TH_MODE)
        call THUpdateAuxVars(realization)
      case default
        if(pflotran_model%option%ntrandof.le.0) then
           pflotran_model%option%io_buffer='pflotranModelSetInitialTHStatesfromCLM ' // &
                 'not implmented for this mode.'
           call printErrMsg(pflotran_model%option)
        endif
    end select

end subroutine pflotranModelSetInternalTHStatesfromCLM

  ! ************************************************************************** !
  ! pflotranModelSetSoilHbcs()
  ! refresh Hydrological BC variables from CLM to PF
  !
  ! by 1-18-2013: only water pressure-head type (dirichlet) available
  ! by 4-11-2013: dirichlet/neumman both available
  ! ************************************************************************** !
  subroutine pflotranModelSetSoilHbcsFromCLM(pflotran_model)

    use Realization_class
    use Option_module
    use Patch_module
    use Grid_module
    use Coupler_module
    use Connection_module

    use TH_Aux_module
    use TH_module
    use Richards_module
    use Richards_Aux_module

    use String_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(simulation_base_type), pointer :: simulation

    type(coupler_type), pointer :: boundary_condition
    type(connection_set_type), pointer :: cur_connection_set
    PetscInt :: ghosted_id, local_id, press_dof, iconn

    PetscErrorCode     :: ierr

    PetscScalar, pointer :: press_maxponding_pf_loc(:)  ! subsurface top boundary max. ponding pressure (Pa) (seepage BC)
    PetscScalar, pointer :: press_subsurf_pf_loc(:)     ! subsurface top boundary pressure-head (Pa) (dirichlet BC)
    PetscScalar, pointer :: qflux_subsurf_pf_loc(:)     ! subsurface top boundary infiltration rate (m/s) (neumann BC)
    PetscScalar, pointer :: press_subbase_pf_loc(:)     ! bottom boundary pressure-head (Pa) (dirichlet BC)
    PetscScalar, pointer :: qflux_subbase_pf_loc(:)     ! botoom boundary drainage flow rate (m/s) (neumann BC)

    PetscScalar, pointer :: toparea_p(:)                ! subsurface top area saved

    !------------------------------------------------------------------------------------

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: pflotranModelSetSoilTHbcs only works on subsurface simulations."
         call printErrMsg(pflotran_model%option)
    end select

    patch           => realization%patch
    grid            => patch%grid

    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_subsurf_clmp, &
                                    clm_pf_idata%press_subsurf_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    pflotran_model%option, &
                                    clm_pf_idata%qflux_subsurf_clmp, &
                                    clm_pf_idata%qflux_subsurf_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_maxponding_clmp, &
                                    clm_pf_idata%press_maxponding_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dbot_to_pf_2dbot, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_subbase_clmp, &
                                    clm_pf_idata%press_subbase_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dbot_to_pf_2dbot, &
                                    pflotran_model%option, &
                                    clm_pf_idata%qflux_subbase_clmp, &
                                    clm_pf_idata%qflux_subbase_pfs)

    ! interface vecs of PF
    call VecGetArrayF90(clm_pf_idata%press_subsurf_pfs,  press_subsurf_pf_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%qflux_subsurf_pfs,  qflux_subsurf_pf_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%press_subbase_pfs,  press_subbase_pf_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%qflux_subbase_pfs, qflux_subbase_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%press_maxponding_pfs, press_maxponding_pf_loc, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(clm_pf_idata%area_top_face_pfp, toparea_p, ierr)
    CHKERRQ(ierr)

    ! passing from interface to internal
    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        press_dof = RICHARDS_PRESSURE_DOF
      case (TH_MODE)
        press_dof = TH_PRESSURE_DOF
      case default
        pflotran_model%option%io_buffer='pflotranModelSetTHbcs ' // &
          'not implmented for this mode.'
        call printErrMsg(pflotran_model%option)
    end select

    boundary_condition => patch%boundary_condition_list%first
    do
       if (.not.associated(boundary_condition)) exit

       cur_connection_set => boundary_condition%connection_set

       do iconn = 1, cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          if (patch%imat(ghosted_id) <= 0) cycle

          if(StringCompare(boundary_condition%name,'clm_gflux_bc') .and. &
             boundary_condition%flow_condition%itype(press_dof) == NEUMANN_BC) then
                   boundary_condition%flow_aux_real_var(press_dof,iconn)= &
                       qflux_subsurf_pf_loc(iconn)

             cur_connection_set%area(iconn) = toparea_p(local_id)     ! normally it's ON (MPI vec, it's from 'local_id')
             if(press_subsurf_pf_loc(iconn) > clm_pf_idata%pressure_reference) then         ! shut-off the BC by resetting the BC 'area' to a tiny value
                cur_connection_set%area(iconn) = 0.d0
             endif

          endif

          if(StringCompare(boundary_condition%name,'clm_gpress_bc') .and. &
             boundary_condition%flow_condition%itype(press_dof) == DIRICHLET_BC) then
                   boundary_condition%flow_aux_real_var(press_dof,iconn)= &
                       press_subsurf_pf_loc(iconn)

             cur_connection_set%area(iconn) = 0.d0               ! normally shut-off this BC
             if(press_subsurf_pf_loc(iconn) > clm_pf_idata%pressure_reference) then         ! turn on the BC by resetting the BC 'area' to real value
                cur_connection_set%area(iconn) = toparea_p(local_id)

#ifdef CLM_PF_DEBUG
     ! the following shows BC connection IS matching up exactly with surface control volume id from CLM
     ! probably because it's in 2D. but for toparea_p, it's in 3D (all cells, not only surface)
      write(pflotran_model%option%myrank+200,*) 'checking H-PRESS. -pf_model-setSoilHbc:', &
        'rank=',pflotran_model%option%myrank, 'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'iconn=',iconn, 'press_top(iconn)=',press_subsurf_pf_loc(iconn), &
        'toparea_p(iconn)=', toparea_p(local_id),&
        'press_dof=',press_dof, &
        'bc_itype=',boundary_condition%flow_condition%itype(press_dof)
#endif
             endif

          endif

          if(StringCompare(boundary_condition%name,'clm_bflux_bc')) then
              if (boundary_condition%flow_condition%itype(press_dof) == DIRICHLET_BC) then
                   boundary_condition%flow_aux_real_var(press_dof,iconn)= &
                       press_subbase_pf_loc(iconn)
              else if (boundary_condition%flow_condition%itype(press_dof) == NEUMANN_BC) then
                   boundary_condition%flow_aux_real_var(press_dof,iconn)= &
                       qflux_subbase_pf_loc(iconn)
              end if
          endif

       enddo

       boundary_condition => boundary_condition%next

    enddo

    call VecRestoreArrayF90(clm_pf_idata%press_subsurf_pfs, press_subsurf_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%qflux_subsurf_pfs, qflux_subsurf_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_subbase_pfs, press_subbase_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%qflux_subbase_pfs, qflux_subbase_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_maxponding_pfs, press_maxponding_pf_loc, ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayF90(clm_pf_idata%area_top_face_pfp, toparea_p, ierr)
    CHKERRQ(ierr)

    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        call RichardsUpdateAuxVars(realization)
      case (TH_MODE)
        call THUpdateAuxVars(realization)
      case default
        pflotran_model%option%io_buffer='pflotranModelSetTHbcs ' // &
          'not implmented for this mode.'
        call printErrMsg(pflotran_model%option)
    end select

  end subroutine pflotranModelSetSoilHbcsFromCLM


  ! ************************************************************************** !
  ! pflotranModelGetRTspecies:
  !  PF RT bgc species Name and index (idof)
  !  Then, all indices is saved for using in this module
  !  so if modification needed, only this subroutine
  !     and the 'ispec_*' and 'name_*' put in the header of this module are modified.
  !
  ! ************************************************************************** !
  subroutine pflotranModelGetRTspecies(pflotran_model)

    use Option_module
    use Realization_class
    use Reaction_Aux_module
    use Reaction_Immobile_Aux_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use clm_pflotran_interface_data

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer  :: pflotran_model
    class(realization_type), pointer    :: realization
    type(simulation_base_type), pointer :: simulation

    PetscInt           :: k
    PetscErrorCode     :: ierr

    character(len=MAXWORDLENGTH) :: word

    !----------------------------------------

    if (pflotran_model%option%ntrandof <= 0) return

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer &
          = "ERROR: SetInitialConcentrations not supported under this mode.."
         call printErrMsg(pflotran_model%option)
    end select

    !
    !immobile species for liter and SOM (decomposing pools)
    if (.not.associated(ispec_decomps_c)) &
      allocate(ispec_decomps_c(clm_pf_idata%ndecomp_pools))
    ispec_decomps_c(:) = 0
    if (.not.associated(ispec_decomps_n)) &
      allocate(ispec_decomps_n(clm_pf_idata%ndecomp_pools))
    ispec_decomps_n(:) = 0
    if (.not.associated(name_decomps)) &
    allocate(name_decomps(clm_pf_idata%ndecomp_pools))
    name_decomps(:) = ""

    do k=1, clm_pf_idata%ndecomp_pools

      ! NOTE: the PF soil bgc sandbox 'SomDec' has a naming protocol as following
      ! (1) fixed-CN ratio decomposing pool: only C pool name is defined, while N is not needed;
      ! (2) varying-CN ratio decomposing pool: 2 pool names are defined with ending letter 'C' or 'N'

      if (clm_pf_idata%floating_cn_ratio(k)) then
        word = trim(name_decomps(k)) // "C"
        ispec_decomps_c(k) = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
        !
        if (ispec_decomps_c(k) <= 0) then
          pflotran_model%option%io_buffer = 'CLM decomposing pool ' // &
            trim(word) // &
            'in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
          call printErrMsg(option)
        endif

        word = trim(name_decomps(k)) // "N"
        ispec_decomps_n(k) = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

        if (ispec_decomps_n(k) <= 0) then
          pflotran_model%option%io_buffer = 'CLM decomposing pool ' // &
            trim(word) // &
            'in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
          call printErrMsg(option)
        endif
      !

      else

        word = trim(name_decomps(k))
        ispec_decomps_c(k) = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

        if (ispec_decomps_c(k) <= 0) then
          pflotran_model%option%io_buffer = 'CLM decomposing pool ' // &
            trim(word) // &
            'in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
          call printErrMsg(option)
        endif

        !
        ispec_decomps_n(k) = UNINITIALIZED_INTEGER

      endif

    end do

    ! aq. species in soil solution/absorbed
    word = name_nh4
    ispec_nh4  = GetPrimarySpeciesIDFromName(word, &
                  realization%reaction,PETSC_FALSE,realization%option)
    if (ispec_nh4 > 0) then
      word = name_nh4sorb
      ispec_nh4sorb  = GetPrimarySpeciesIDFromName(word, &
                  realization%reaction,PETSC_FALSE,realization%option)
    else
      ispec_nh4sorb  = UNINITIALIZED_INTEGER
    endif

    word = name_no3
    ispec_no3  = GetPrimarySpeciesIDFromName(word, &
                  realization%reaction,PETSC_FALSE,realization%option)

    !
    ! species for gases
    word = name_co2
    ispec_co2  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    word = name_n2
    ispec_n2  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = name_n2o
    ispec_n2o = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    ! CO2 productin from C decomposition reaction network for tracking (immoble species)
    word = name_hr
    ispec_hr = GetImmobileSpeciesIDFromName(word, &
            realization%reaction%immobile,PETSC_FALSE,realization%option)

    ! N bgc reaction fluxes for tracking (immoble species)
    word = name_plantndemand
    ispec_plantndemand  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    word = name_plantnh4uptake
    ispec_plantnh4uptake  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    word = name_plantno3uptake
    ispec_plantno3uptake  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = name_nmin
    ispec_nmin = GetImmobileSpeciesIDFromName(word, &
            realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = name_nimmp
    ispec_nimmp = GetImmobileSpeciesIDFromName(word, &
            realization%reaction%immobile,PETSC_FALSE,realization%option)
    word = name_nimm
    ispec_nimm = GetImmobileSpeciesIDFromName(word, &
            realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = name_ngasmin
    ispec_ngasmin = GetImmobileSpeciesIDFromName(word, &
           realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = name_ngasnitr
    ispec_ngasnitr = GetImmobileSpeciesIDFromName( word, &
           realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = name_ngasdeni
    ispec_ngasdeni = GetImmobileSpeciesIDFromName( word, &
           realization%reaction%immobile,PETSC_FALSE,realization%option)

  end subroutine pflotranModelGetRTspecies

  ! ************************************************************************** !
  ! TEMPORARILY OFF!
  !> This routine Pass CLM SOM decomposition rate constants for PFLOTRAN bgc
  !! So that both are consistent
  !!
  !> @author
  !! F.-M. Yuan
  !!
  !! Date: 12/5/2014
  ! ************************************************************************** !
!  subroutine pflotranModelSetSOMKfromCLM(pflotran_model, pf_cmode)
!
!    use Option_module
!    use Reaction_module
!    use Reaction_Aux_module
!    use Reactive_Transport_Aux_module
!    use Global_Aux_module
!    use Material_Aux_class, only: material_auxvar_type
!    use Reaction_Sandbox_module
!    use Reaction_Sandbox_SomDec_class
!
!    use Simulation_Base_class, only : simulation_base_type
!    use Simulation_Subsurface_class, only : subsurface_simulation_type
!    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
!
!    implicit none
!
!    type(option_type) :: option
!    type(reaction_type) :: reaction
!    type(reactive_transport_auxvar_type) :: rt_auxvar
!    type(global_auxvar_type) :: global_auxvar
!    class(material_auxvar_type) :: material_auxvar
!
!    class(reaction_sandbox_base_type), pointer :: cur_reaction
!
!    cur_reaction => rxn_sandbox_list
!    do
!        if (.not.associated(cur_reaction)) exit
!
!
!
!        cur_reaction => cur_reaction%next
!    enddo
!
!  end subroutine pflotranModelSetSOMKfromCLM

  ! ************************************************************************** !
  !
  ! pflotranModelSetBgcConc:
  !  Get CLM concentrations (C, N), converts from CLM units into PFLOTRAN units.
  !  and set concentrations in PFLOTRAN
  ! ************************************************************************** !
  subroutine pflotranModelSetBgcConcFromCLM(pflotran_model)

    use Global_Aux_module
    use Realization_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_Aux_module
    use Reaction_Immobile_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module
    use Discretization_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(field_type), pointer                 :: field
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(simulation_base_type), pointer       :: simulation
    type(global_auxvar_type), pointer :: global_auxvars(:)

    PetscErrorCode     :: ierr
    PetscInt           :: local_id
    PetscInt           :: ghosted_id
    PetscReal, pointer :: xx_p(:)

    PetscScalar, pointer :: decomp_cpools_vr_pf_loc(:) ! (molesC/m3)
    PetscScalar, pointer :: decomp_npools_vr_pf_loc(:) ! (molesN/m3)
    PetscScalar, pointer :: smin_no3_vr_pf_loc(:)      ! (molesN/m3)
    PetscScalar, pointer :: smin_nh4_vr_pf_loc(:)      ! (molesN/m3)
    PetscScalar, pointer :: smin_nh4sorb_vr_pf_loc(:)  ! (molesN/m3)

    PetscReal, pointer :: porosity_loc_p(:)

    PetscInt  :: offset, offsetim
    PetscReal :: porosity, saturation, theta ! for concentration conversion from mol/m3 to mol/L

    PetscInt, pointer :: idecomp_clmp_index(:), idecomp_pfs_index(:)
    Vec               :: vec_clmp, vec_pfs
    PetscInt          :: j, k, vec_ioffset

!------------------------------------------------------------------------------------------------

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer &
          = "ERROR: SetInitialConcentrations not supported under this mode.."
         call printErrMsg(pflotran_model%option)
    end select

    patch => realization%patch
    grid  => patch%grid
    field => realization%field
    global_auxvars  => patch%aux%Global%auxvars

    !-----------------------------------------------------------------

    ! create temporary vecs/arrays for each 'decomp_pool' data-mapping
    call VecDuplicate(clm_pf_idata%zsoi_clmp, vec_clmp,ierr)
    CHKERRQ(ierr)
    call VecDuplicate(clm_pf_idata%zsoi_pfs, vec_pfs,ierr)
    CHKERRQ(ierr)
    allocate(idecomp_clmp_index(clm_pf_idata%ngclm_sub))
    do j=1, clm_pf_idata%ngclm_sub
      idecomp_clmp_index(j) = j-1
    enddo
    allocate(idecomp_pfs_index(clm_pf_idata%nlpf_sub))
    do j=1, clm_pf_idata%nlpf_sub
      idecomp_pfs_index(j) = j-1
    enddo

    ! decomps'C'
    if (associated(ispec_decompc)) then
      do k=1,clm_pf_idata%ndecomp_pools
        ! get a seg. of data from the whole '_clmp' vec for the 'k'th pool
        call VecGetValues(clm_pf_idata%decomp_cpools_vr_clmp,               &
                          idecomp_clmp_index+(k-1)*clm_pf_idata%ngclm_sub,  &
                          clm_pf_idata%ngclm_sub,                           &
                          vec_clmp, INSERT_VALUES, ierr)
        CHKERRQ(ierr)
        ! assembly the 'vec_clmp'
        call VecAssemblyBegin(vec_clmp, ierr)
        CHKERRQ(ierr)
        call VecAssemblyEnd(vec_clmp, ierr)
        CHKERRQ(ierr)

        ! mapping
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option,                    &
                                    vec_clmp,                                 &
                                    vec_pfs)

        ! insert 'vec_pfs' into the whole '_pfs' vec
        call VecSetValues(clm_pf_idata%decomp_cpools_vr_pfs,              &
                          idecomp_pfs_index+(k-1)*clm_pf_idata%nlpf_sub,  &
                          clm_pf_idata%nlpf_sub,                          &
                          vec_pfs, INSERT_VALUES, ierr)
        CHKERRQ(ierr)

      enddo

      ! assembly the whole '_pfs' vec
      call VecAssemblyBegin(clm_pf_idata%decomp_cpools_vr_pfs, ierr)
      CHKERRQ(ierr)
      call VecAssemblyEnd(clm_pf_idata%decomp_cpools_vr_pfs, ierr)
      CHKERRQ(ierr)

    endif

    ! decomps_'N'
    if (associated(ispec_decompn)) then
      do k=1,clm_pf_idata%ndecomp_pools
        ! get a seg. of data from the whole '_clmp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%ngclm_sub       ! decomps_clmp vec: 'cell' first, then 'species'
        call VecGetValues(clm_pf_idata%decomp_npools_vr_clmp,   &
                          idecomp_clmp_index+vec_offset,        &       ! decomps_clmp vec: 'cell' first, then 'species'
                          clm_pf_idata%ngclm_sub,               &
                          vec_clmp, INSERT_VALUES, ierr)
        CHKERRQ(ierr)
        ! assembly the 'vec_clmp'
        call VecAssemblyBegin(vec_clmp, ierr)
        CHKERRQ(ierr)
        call VecAssemblyEnd(vec_clmp, ierr)
        CHKERRQ(ierr)

        ! mapping
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option,                    &
                                    vec_clmp,                                 &
                                    vec_pfs)

        ! insert 'vec_pfs' into the whole '_pfs' vec
        vec_offset = (k-1)*clm_pf_idata%nlpf_sub       ! decomps_pfs vec: 'cell' first, then 'species'
        call VecSetValues(clm_pf_idata%decomp_npools_vr_pfs,    &
                          idecomp_pfs_index+vec_offset,         &       ! decomps_pfs vec: 'cell' first, then 'species'
                          clm_pf_idata%nlpf_sub,                &
                          vec_pfs, INSERT_VALUES, ierr)
        CHKERRQ(ierr)

      enddo

      ! assembly the whole '_pfs' vec
      call VecAssemblyBegin(clm_pf_idata%decomp_npools_vr_pfs, ierr)
      CHKERRQ(ierr)
      call VecAssemblyEnd(clm_pf_idata%decomp_npools_vr_pfs, ierr)
      CHKERRQ(ierr)
      enddo
    endif

    ! clear-up of temporary vecs/arrarys
    call VecDestroy(vec_clmp,ierr)
    CHKERRQ(ierr)
    call VecDestroy(vec_pfs,ierr)
    CHKERRQ(ierr)
    deallocate(idecomp_clmp_index)
    deallocate(idecomp_pfs_index)



    ! for individual vecs

    if(ispec_no3 > 0) then
       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%smin_no3_vr_clmp, &
                                    clm_pf_idata%smin_no3_vr_pfs)
    endif

    if(ispec_nh4 > 0) then
       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%smin_nh4_vr_clmp, &
                                    clm_pf_idata%smin_nh4_vr_pfs)
    endif

    !----------------------------------------------------------------------------

    if(associated(ispec_decomps_c)) then
      call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_pfs, &
                        decomp_cpools_vr_pf_loc, ierr)
      CHKERRQ(ierr)
    endif

    if(associated(ispec_decomps_n)) then
      call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_pfs, &
                        decomp_npools_vr_pf_loc, ierr)
      CHKERRQ(ierr)
    endif

    if(ispec_no3 > 0) then
      call VecGetArrayF90(clm_pf_idata%smin_no3_vr_pfs, smin_no3_vr_pf_loc, ierr)
      CHKERRQ(ierr)
    endif

    if(ispec_nh4 > 0) then
      call VecGetArrayF90(clm_pf_idata%smin_nh4_vr_pfs, smin_nh4_vr_pf_loc, ierr)
      CHKERRQ(ierr)
    endif

    call VecGetArrayF90(field%tran_xx,xx_p,ierr)
    CHKERRQ(ierr)

    call VecGetArrayReadF90(field%porosity0, porosity_loc_p, ierr)
    CHKERRQ(ierr)

    !----------------------------------------------------------------------------------------------------
    do local_id = 1, grid%nlmax
      ghosted_id=grid%nL2G(local_id)
      if (ghosted_id<=0 .or. local_id <= 0) cycle ! bypass ghosted corner cells
      !if (patch%imat(local_id) <= 0) cycle  !(TODO) imat IS 0 for some cells when decomposing domain in X and Y directions.

      offset = (local_id - 1)*realization%reaction%ncomp

      saturation = global_auxvars(ghosted_id)%sat(1)          ! using 'ghosted_id' if from 'global_auxvars'
      porosity = porosity_loc_p(local_id)                     ! using 'local_id' if from 'field%???'
      theta = saturation * porosity

      if(ispec_no3 > 0) then
         xx_p(offset + ispec_no3) = max(xeps0_n, smin_no3_vr_pf_loc(ghosted_id)  &      ! from 'ghosted_id' to field%xx_p's local
                                                 / theta / 1000.0d0)
      endif

      if(ispec_nh4 > 0) then
         xx_p(offset + ispec_nh4) = max(xeps0_n, smin_nh4_vr_pf_loc(ghosted_id)  &
                                                 / theta / 1000.d0)
      endif

      !
      offsetim = offset + realization%reaction%offset_immobile

      do k=1, clm_pf_idata%ndecomp_pools
        vec_offset = (k-1)*clm_pf_idata%nlpf_sub       ! decomps_pfs vec: 'cell' first, then 'species'
        if(ispec_decomps_c(k) > 0) then
          xx_p(offsetim + ispec_decomps_c(k)) = max( xeps0_c, &
                       decomp_cpools_vr_pf_loc(ghosted_id+vec_offset) )
        endif

        if(ispec_decomps_n(k) > 0) then
          xx_p(offsetim + ispec_decomps_n(k)) = max( xeps0_n, &
                       decomp_npools_vr_pf_loc(ghosted_i+vec_offsetd) )
        endif
      enddo

#ifdef CLM_PF_DEBUG
      !F.-M. Yuan: the following IS a checking, comparing CLM passed data (NH4):
      ! Conclusions: (1) local_id runs from 1 ~ grid%nlmax; and ghosted_id is obtained by 'nL2G' as corrected above;
      !              OR, ghosted_id runs from 1 ~ grid%ngmax; and local_id is obtained by 'nG2L'.
      !              (2) data-passing IS by from 'ghosted_id' to 'local_id'
      write(pflotran_model%option%myrank+200,*) 'checking bgc - pflotran-model setting init. conc.: ', &
        'rank=',pflotran_model%option%myrank, &
        'local_id=',local_id, 'ghosted_id=',ghosted_id, 'xxp_nh4_id', offset+ispec_nh4, &
        'smin_nh4_pfs(ghosted_id)=',smin_nh4_vr_pf_loc(ghosted_id), &
        'xx_p(xxp_nh4_id)=',xx_p(offset + ispec_nh4), &
        'sat_glob(ghosted_id)=',global_auxvars(ghosted_id)%sat(1), &
        'poro(local_id)=',porosity_loc_p(local_id)

#endif

    enddo

    !----------------------------------------------------------------------------------------------------
    if(associated(ispec_decomps_c)) then
      call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_pfs, &
                        decomp_cpools_vr_pf_loc, ierr)
      CHKERRQ(ierr)
    endif

    if(associated(ispec_decomps_n)) then
      call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_pfs, &
                        decomp_npools_vr_pf_loc, ierr)
      CHKERRQ(ierr)
    endif

    if(ispec_no3 > 0) then
      call VecRestoreArrayF90(clm_pf_idata%smin_no3_vr_pfs, smin_no3_vr_pf_loc, ierr)
      CHKERRQ(ierr)
    endif

    if(ispec_nh4 > 0) then
      call VecRestoreArrayF90(clm_pf_idata%smin_nh4_vr_pfs, smin_nh4_vr_pf_loc, ierr)
      CHKERRQ(ierr)
    endif

    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(field%porosity0, porosity_loc_p, ierr)
    CHKERRQ(ierr)

    call DiscretizationGlobalToLocal(realization%discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

    call VecCopy(field%tran_xx,field%tran_yy,ierr)
    CHKERRQ(ierr)
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

  end subroutine pflotranModelSetBgcConcFromCLM

  ! ************************************************************************** !
  ! pflotranModelSetBGCRatesFromCLM:
  !  Get CLM litter, som, mineral N production and plant demand rates
  !  Convert from CLM units into PFLOTRAN units.
  !  Set values in PFLOTRAN
  ! ************************************************************************** !
  subroutine pflotranModelSetBGCRatesFromCLM(pflotran_model)

    use Realization_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_Aux_module
    use Reaction_Immobile_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module
    use Discretization_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
    use Mass_Transfer_module, only : mass_transfer_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer  :: pflotran_model
    class(realization_type), pointer    :: realization
    type(field_type), pointer           :: field
    type(patch_type), pointer           :: patch
    type(grid_type), pointer            :: grid
    type(simulation_base_type), pointer :: simulation
    type(mass_transfer_type), pointer   :: cur_mass_transfer

    PetscErrorCode     :: ierr
    PetscInt           :: ghosted_id,local_id
    PetscInt           :: offset, offsetim

    PetscScalar, pointer :: rate_pf_loc(:)   !

    PetscReal, pointer :: volume_p(:)

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer &
          = "ERROR: SetInitialConcentrations not supported under this mode.."
         call printErrMsg(pflotran_model%option)
    end select

    patch => realization%patch
    grid  => patch%grid
    field => realization%field

    ! mapping data from CLM to PFLOTRAN
    if(ispec_lit1c >0) then
      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_decompsc_clmp, &
                                    clm_pf_idata%rate_decompsc_pfs)
    endif


    if(ispec_lit1n >0) then
      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_decompsn_clmp, &
                                    clm_pf_idata%rate_decompsn_pfs)
    endif


    ! NOTE: direct data passing from interface to PF for N demand
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_plantndemand_clmp, &
                                    clm_pf_idata%rate_plantndemand_pfs)

    if(ispec_no3 >0) then
      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_smin_no3_clmp, &
                                    clm_pf_idata%rate_smin_no3_pfs)
    endif

    if(ispec_nh4 >0) then
      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_smin_nh4_clmp, &
                                    clm_pf_idata%rate_smin_nh4_pfs)
    endif

    !----------------------------------------------------------------------

    ! get cell volume to convert mass transfer rate unit from moles/m3/s to moles/s
    call VecGetArrayReadF90(field%volume0,volume_p,ierr)
    CHKERRQ(ierr)

    if (associated(realization%rt_mass_transfer_list)) then

       offsetim = realization%reaction%offset_immobile

       cur_mass_transfer => realization%rt_mass_transfer_list

       do
         if (.not.associated(cur_mass_transfer)) exit

         if(cur_mass_transfer%idof == ispec_nh4) then
           call VecGetArrayReadF90(clm_pf_idata%rate_smin_nh4_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_no3) then
           call VecGetArrayReadF90(clm_pf_idata%rate_smin_no3_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)

         elseif(cur_mass_transfer%idof == ispec_lit1c+offsetim) then
           call VecGetArrayReadF90(clm_pf_idata%rate_lit1c_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_lit1n+offsetim) then
           call VecGetArrayReadF90(clm_pf_idata%rate_lit1n_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         endif

         do local_id = 1, grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (local_id<=0 .or. ghosted_id<=0) cycle
            !if (patch%imat(local_id) <= 0) cycle  !(TODO) imat IS 0 for some cells when decomposing domain in X and Y directions.

            if(cur_mass_transfer%idof == ispec_nh4                  &
               .or. cur_mass_transfer%idof == ispec_no3             &
               .or. cur_mass_transfer%idof == ispec_lit1c+offsetim  &
               .or. cur_mass_transfer%idof == ispec_lit2c+offsetim  &
               .or. cur_mass_transfer%idof == ispec_lit3c+offsetim  &
               .or. cur_mass_transfer%idof == ispec_lit1n+offsetim  &
               .or. cur_mass_transfer%idof == ispec_lit2n+offsetim  &
               .or. cur_mass_transfer%idof == ispec_lit3n+offsetim  &
               .or. cur_mass_transfer%idof == ispec_cwdc+offsetim  &
               .or. cur_mass_transfer%idof == ispec_cwdn+offsetim  &
               .or. cur_mass_transfer%idof == ispec_som1c+offsetim  &
               .or. cur_mass_transfer%idof == ispec_som2c+offsetim  &
               .or. cur_mass_transfer%idof == ispec_som3c+offsetim  &
               .or. cur_mass_transfer%idof == ispec_som4c+offsetim  &
               .or. cur_mass_transfer%idof == ispec_som1n+offsetim  &
               .or. cur_mass_transfer%idof == ispec_som2n+offsetim  &
               .or. cur_mass_transfer%idof == ispec_som3n+offsetim  &
               .or. cur_mass_transfer%idof == ispec_som4n+offsetim  &
              ) then

#ifdef CLM_PF_DEBUG
      !F.-M. Yuan: the following IS a checking, comparing CLM passed data (mass transfer rate):
      ! Conclusions: (1) local_id runs from 1 ~ grid%nlmax; and ghosted_id is obtained by 'nL2G' as corrected above;
      !              OR, ghosted_id runs from 1 ~ grid%ngmax; and local_id is obtained by 'nG2L'.
      !              (2) data-passing IS by from 'ghosted_id' to 'local_id'
      if (cur_mass_transfer%idof == ispec_nh4) &
      write(pflotran_model%option%myrank+200,*) 'checking bgc-mass-rate - pflotran_model: ', &
        'rank=',pflotran_model%option%myrank, 'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'rate_nh4_pfs(ghosted_id)=',rate_pf_loc(ghosted_id), &
        'masstransfer_nh4_predataset(local_id)=',cur_mass_transfer%dataset%rarray(local_id)
#endif

               cur_mass_transfer%dataset%rarray(local_id) = &
                        rate_pf_loc(ghosted_id)*volume_p(local_id)  ! mol/m3s * m3

            endif
         enddo

         if(cur_mass_transfer%idof == ispec_nh4) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_smin_nh4_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_no3) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_smin_no3_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)

         elseif(cur_mass_transfer%idof == ispec_lit1c+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_lit1c_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_lit2c+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_lit2c_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_lit3c+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_lit3c_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)

         elseif(cur_mass_transfer%idof == ispec_lit1n+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_lit1n_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_lit2n+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_lit2n_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_lit3n+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_lit3n_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)

         elseif(cur_mass_transfer%idof == ispec_cwdc+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_cwdc_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_cwdn+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_cwdn_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)

         elseif(cur_mass_transfer%idof == ispec_som1c+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_som1c_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_som2c+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_som2c_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_som3c+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_som3c_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_som4c+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_som4c_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)

        elseif(cur_mass_transfer%idof == ispec_som1n+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_som1n_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_som2n+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_som2n_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_som3n+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_som3n_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)
         elseif(cur_mass_transfer%idof == ispec_som4n+offsetim) then
           call VecRestoreArrayReadF90(clm_pf_idata%rate_som4n_pfs, rate_pf_loc, ierr)
           CHKERRQ(ierr)

         endif

         cur_mass_transfer => cur_mass_transfer%next
       enddo
    endif

    call VecRestoreArrayReadF90(field%volume0,volume_p,ierr)
    CHKERRQ(ierr)

  end subroutine pflotranModelSetBGCRatesFromCLM

  ! ************************************************************************** !
  !
  ! pflotranModelUpdateAqConcFromCLM:
  !  Get CLM aqueous nutrient concentrations (NH4, NO3 at this momment),
  !  converts from CLM units into PFLOTRAN units, and reset concentrations in PFLOTRAN
  !
  !   notes: when NOT coupled with PF Hydrology, forcing CLM water saturation to reset
  !          PF's global saturation status WOULD cause aq. phase element mass balance issue
  !
  ! ************************************************************************** !
  subroutine pflotranModelUpdateAqConcFromCLM(pflotran_model)

    use Global_Aux_module
    use Realization_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module
    use Discretization_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(field_type), pointer                 :: field
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(simulation_base_type), pointer       :: simulation
    type(global_auxvar_type), pointer         :: global_auxvars(:)
    type(aq_species_type), pointer            :: cur_aq_spec

    PetscErrorCode     :: ierr
    PetscInt           :: local_id
    PetscInt           :: ghosted_id
    PetscReal, pointer :: xx_p(:)

    PetscReal, pointer   :: porosity0_loc_p(:)              ! current CLM-updated porosity

    PetscReal, pointer :: porosity_pre_pf_loc(:)            !previous time-step porosity (m3/m3 bulk soil)
    PetscReal, pointer :: soillsat_pre_pf_loc(:)            !previous time-step soil liq. water saturation (0 - 1)

    PetscInt :: offset, offset_aq

    character(len=MAXWORDLENGTH) :: word
    PetscReal :: porosity, saturation, theta ! for concentration conversion from mol/m3 to mol/L
    PetscReal :: theta_pre

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer &
          = "ERROR: SetInitialConcentrations not supported under this mode.."
         call printErrMsg(pflotran_model%option)
    end select

    patch => realization%patch
    grid  => patch%grid
    field => realization%field
    global_auxvars  => patch%aux%Global%auxvars
    call VecGetArrayF90(field%tran_xx,xx_p,ierr)
    CHKERRQ(ierr)

    call VecGetArrayReadF90(field%porosity0, porosity0_loc_p, ierr)
    CHKERRQ(ierr)

    ! the previous time-step porosity and lsat from PFLOTRAN, saved in clm_pf_idata
    ! NOTE: make sure NOT modified by CLM
    call VecGetArrayReadF90(clm_pf_idata%effporosity_pfp, porosity_pre_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayReadF90(clm_pf_idata%soillsat_pfp, soillsat_pre_pf_loc, ierr)
    CHKERRQ(ierr)

    do local_id = 1, grid%nlmax
      ghosted_id=grid%nL2G(local_id)
      if (ghosted_id<=0 .or. local_id <= 0) cycle ! bypass ghosted corner cells
      !if (patch%imat(local_id) <= 0) cycle  !(TODO) imat IS 0 for some cells when decomposing domain in X and Y directions.

      offset = (local_id - 1)*realization%reaction%ncomp
      offset_aq = offset + realization%reaction%offset_aqueous

      ! at this moment, both 'saturtion' and 'porosity' have been updated from CLM and pass to PF
      ! otherwise, it's no use.
      saturation = global_auxvars(ghosted_id)%sat(1)          ! using 'ghosted_id' if from 'global_auxvars'
      porosity = porosity0_loc_p(local_id)                    ! using 'local_id' if from 'field%???'
      theta = saturation * porosity

      ! adjusting aq. species conc. due to CLM pass-in theta (porosity X saturation)

      ! previous timestep saved PF's 'porosity' and 'saturation' in clm-pf-idata%???_pfp
      theta_pre = porosity_pre_pf_loc(local_id)*  &
                  soillsat_pre_pf_loc(local_id)

      cur_aq_spec => realization%reaction%primary_species_list
      do
        if (.not.associated(cur_aq_spec)) exit

        if (theta_pre> 1.d-20 .and. theta > 1.0d-20) then
          xx_p(offset_aq + cur_aq_spec%id) = xx_p(offset_aq + cur_aq_spec%id) &
                                            * theta_pre / theta
        end if
        cur_aq_spec => cur_aq_spec%next
      enddo

    enddo

    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(field%porosity0, porosity0_loc_p, ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayReadF90(clm_pf_idata%effporosity_pfp, porosity_pre_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayReadF90(clm_pf_idata%soillsat_pfp, soillsat_pre_pf_loc, ierr)
    CHKERRQ(ierr)

    call DiscretizationGlobalToLocal(realization%discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

    call VecCopy(field%tran_xx,field%tran_yy,ierr)
    CHKERRQ(ierr)
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

  end subroutine pflotranModelUpdateAqConcFromCLM

  ! ************************************************************************** !
  !> This routine reset gas concenctration (as immobile species) after emission adjusting
  !> in CLM-PFLOTRAN interface
  !! Note: this is a temporary work-around, because PF doesn't have BGC gas
  !!       processes at this moment
  !> @author
  !! F.-M. Yuan
  !!
  !! date: 3/19/2014
  ! ************************************************************************** !
  subroutine pflotranModelUpdateAqGasesFromCLM(pflotran_model)

    use Global_Aux_module
    use Realization_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Discretization_module
    use Reaction_Aux_module
    use Reaction_Immobile_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(field_type), pointer                 :: field
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(simulation_base_type), pointer       :: simulation
    type(global_auxvar_type), pointer :: global_auxvars(:)

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: xx_p(:)

    PetscScalar, pointer :: gco2_vr_pf_loc(:)              ! (molC/m3 bulk soil)
    PetscScalar, pointer :: gn2_vr_pf_loc(:)               ! (molN2-N/m3 bulk soil)
    PetscScalar, pointer :: gn2o_vr_pf_loc(:)              ! (molN2O-N/m3 bulk soil)

    PetscInt :: offset

    character(len=MAXWORDLENGTH) :: word

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer &
          = "ERROR: SetInitialConcentrations not supported under this mode.."
         call printErrMsg(pflotran_model%option)
    end select

    patch => realization%patch
    grid  => patch%grid
    field => realization%field
    global_auxvars  => patch%aux%Global%auxvars

    ! mapping CLM vecs to PF vecs
    if(ispec_co2 > 0) then
       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%gco2_vr_clmp, &
                                    clm_pf_idata%gco2_vr_pfs)
    endif

    if(ispec_n2o > 0) then
       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%gn2o_vr_clmp, &
                                    clm_pf_idata%gn2o_vr_pfs)
    endif

    if(ispec_n2 > 0) then
       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%gn2_vr_clmp, &
                                    clm_pf_idata%gn2_vr_pfs)
    endif

    ! (iii) get the 'PF' vecs for resetting data
    call VecGetArrayF90(clm_pf_idata%gco2_vr_pfs, gco2_vr_pf_loc,ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%gn2_vr_pfs, gn2_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%gn2o_vr_pfs, gn2o_vr_pf_loc, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(field%tran_xx,xx_p,ierr)
    CHKERRQ(ierr)  ! extract data from pflotran internal portion

    do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (associated(patch%imat)) then
           if (patch%imat(ghosted_id) <= 0) cycle
        endif

        offset = (local_id-1) * realization%reaction%ncomp &
                + realization%reaction%offset_immobile

        if(ispec_co2 > 0) then
            xx_p(offset + ispec_co2) = max(gco2_vr_pf_loc(ghosted_id), xeps0_c)
        endif

        if(ispec_n2 > 0) then
            xx_p(offset + ispec_n2) = max(gn2_vr_pf_loc(ghosted_id), xeps0_n)
        endif

        if(ispec_n2o > 0) then
            xx_p(offset + ispec_n2o) = max(gn2o_vr_pf_loc(ghosted_id), xeps0_n)
        endif

    enddo

    call VecRestoreArrayF90(clm_pf_idata%gco2_vr_pfs, gco2_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%gn2_vr_pfs, gn2_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%gn2o_vr_pfs, gn2o_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    !
    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr)
    CHKERRQ(ierr)

    !
    call DiscretizationGlobalToLocal(realization%discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

    call VecCopy(field%tran_xx,field%tran_yy,ierr)
    CHKERRQ(ierr)
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

  end subroutine pflotranModelUpdateAqGasesFromCLM


  ! ************************************************************************** !
  !> This routine get updated bgc states/fluxes evoloved by PFLOTRAN.
  !!
  ! ************************************************************************** !
  subroutine pflotranModelGetBgcVariablesFromPF(pflotran_model)

    use Global_Aux_module
    use Realization_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_module
    use Reaction_Aux_module
    use Reactive_Transport_module
    use Reactive_Transport_Aux_module
    use Reaction_Immobile_Aux_module
    use Discretization_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer          :: realization
    type(field_type), pointer                 :: field
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(simulation_base_type), pointer       :: simulation
    type(reaction_type), pointer              :: reaction
    type(global_auxvar_type), pointer :: global_auxvar
    type(reactive_transport_auxvar_type), pointer :: rt_auxvar

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: xx_p(:)

    PetscScalar, pointer :: decomp_cpools_vr_lit1_pf_loc(:) ! (molesC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_lit2_pf_loc(:) ! (molesC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_lit3_pf_loc(:) ! (molesC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_cwd_pf_loc(:)  ! (molesC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som1_pf_loc(:) ! (molesC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som2_pf_loc(:) ! (molesC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som3_pf_loc(:) ! (molesC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som4_pf_loc(:) ! (molesC/m3)
    PetscScalar, pointer :: decomp_npools_vr_lit1_pf_loc(:) ! (molesN/m3)
    PetscScalar, pointer :: decomp_npools_vr_lit2_pf_loc(:) ! (molesN/m3)
    PetscScalar, pointer :: decomp_npools_vr_lit3_pf_loc(:) ! (molesN/m3)
    PetscScalar, pointer :: decomp_npools_vr_cwd_pf_loc(:)  ! (molesN/m3)
    PetscScalar, pointer :: decomp_npools_vr_som1_pf_loc(:) ! (molesN/m3)
    PetscScalar, pointer :: decomp_npools_vr_som2_pf_loc(:) ! (molesN/m3)
    PetscScalar, pointer :: decomp_npools_vr_som3_pf_loc(:) ! (molesN/m3)
    PetscScalar, pointer :: decomp_npools_vr_som4_pf_loc(:) ! (molesN/m3)
    PetscScalar, pointer :: accextrnh4_vr_pf_loc(:)         ! (molesN/m3)
    PetscScalar, pointer :: accextrno3_vr_pf_loc(:)         ! (molesN/m3)
    PetscScalar, pointer :: smin_no3_vr_pf_loc(:)           ! (molesN/m3)
    PetscScalar, pointer :: smin_nh4_vr_pf_loc(:)           ! (molesN/m3)
    PetscScalar, pointer :: smin_nh4sorb_vr_pf_loc(:)       ! (molesN/m3)
    PetscScalar, pointer :: gco2_vr_pf_loc(:)               ! (molC/m3)
    PetscScalar, pointer :: gn2_vr_pf_loc(:)                ! (molN2/m3)
    PetscScalar, pointer :: gn2o_vr_pf_loc(:)               ! (molN2O/m3)
    PetscScalar, pointer :: acchr_vr_pf_loc(:)              ! (molesC/m3)
    PetscScalar, pointer :: accnmin_vr_pf_loc(:)            ! (molesN/m3)
    PetscScalar, pointer :: accnimmp_vr_pf_loc(:)           ! (molesN/m3)
    PetscScalar, pointer :: accnimm_vr_pf_loc(:)            ! (molesN/m3)
    PetscScalar, pointer :: accngasmin_vr_pf_loc(:)         ! (molesN/m3)
    PetscScalar, pointer :: accngasnitr_vr_pf_loc(:)        ! (molesN/m3)
    PetscScalar, pointer :: accngasdeni_vr_pf_loc(:)        ! (molesN/m3)

    PetscReal, pointer :: porosity_loc_p(:)

    PetscReal :: porosity, saturation, theta ! for concentration conversion from mol/m3 to mol/L
    PetscReal :: conc
    PetscInt  :: offset, offsetim

    PetscReal, parameter :: zeroing_conc = 1.0d-10

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer &
          = "ERROR: SetInitialConcentrations not supported under this mode.."
         call printErrMsg(pflotran_model%option)
    end select

    patch => realization%patch
    grid  => patch%grid
    field => realization%field
    reaction => realization%reaction

    ! (ii) get the original 'pf' vecs
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit1_pfp, &
                        decomp_cpools_vr_lit1_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit2_pfp, &
                        decomp_cpools_vr_lit2_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit3_pfp, &
                        decomp_cpools_vr_lit3_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_cwd_pfp,  &
                        decomp_cpools_vr_cwd_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som1_pfp, &
                        decomp_cpools_vr_som1_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som2_pfp, &
                        decomp_cpools_vr_som2_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som3_pfp, &
                        decomp_cpools_vr_som3_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som4_pfp, &
                        decomp_cpools_vr_som4_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit1_pfp, &
                        decomp_npools_vr_lit1_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit2_pfp, &
                        decomp_npools_vr_lit2_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit3_pfp, &
                        decomp_npools_vr_lit3_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_cwd_pfp, &
                        decomp_npools_vr_cwd_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_som1_pfp, &
                        decomp_npools_vr_som1_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_som2_pfp, &
                        decomp_npools_vr_som2_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_som3_pfp, &
                        decomp_npools_vr_som3_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_som4_pfp, &
                        decomp_npools_vr_som4_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%smin_no3_vr_pfp, smin_no3_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%smin_nh4_vr_pfp, smin_nh4_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%smin_nh4sorb_vr_pfp, smin_nh4sorb_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    !
    call VecGetArrayF90(clm_pf_idata%gco2_vr_pfp, gco2_vr_pf_loc,ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%gn2_vr_pfp, gn2_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%gn2o_vr_pfp, gn2o_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    !
    call VecGetArrayF90(clm_pf_idata%accextrnh4_vr_pfp, accextrnh4_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%accextrno3_vr_pfp, accextrno3_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%acchr_vr_pfp, acchr_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%accnmin_vr_pfp, accnmin_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%accnimmp_vr_pfp, accnimmp_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%accnimm_vr_pfp, accnimm_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%accngasmin_vr_pfp, accngasmin_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%accngasnitr_vr_pfp, accngasnitr_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%accngasdeni_vr_pfp, accngasdeni_vr_pf_loc, ierr)
    CHKERRQ(ierr)

    ! (iii) pass the data from internal to PFLOTRAN vecs

    call VecGetArrayF90(field%tran_xx,xx_p,ierr)
    CHKERRQ(ierr)  ! extract data from pflotran internal portion
    call VecGetArrayReadF90(field%porosity0, porosity_loc_p, ierr)
    CHKERRQ(ierr)

    do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (ghosted_id <= 0) cycle
        if (associated(patch%imat)) then
           if (patch%imat(ghosted_id) <= 0) cycle
        endif

        global_auxvar  => patch%aux%Global%auxvars(ghosted_id)
        rt_auxvar => patch%aux%RT%auxvars(ghosted_id)

        saturation = global_auxvar%sat(1)
        porosity = porosity_loc_p(local_id)
        theta = saturation * porosity

        offset = (local_id - 1)*realization%reaction%ncomp

        offsetim = offset + realization%reaction%offset_immobile

        if (ispec_lit1c > 0) then
          decomp_cpools_vr_lit1_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit1c), 0.d0)
        endif

        if (ispec_lit2c > 0) then
          decomp_cpools_vr_lit2_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit2c), 0.d0)
        endif

        if (ispec_lit3c > 0) then
          decomp_cpools_vr_lit3_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit3c), 0.d0)
        endif

        if (ispec_lit1n > 0) then
          decomp_npools_vr_lit1_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit1n), 0.d0)
        endif

        if (ispec_lit2n > 0) then
          decomp_npools_vr_lit2_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit2n), 0.d0)
        endif

        if (ispec_lit3n > 0) then
          decomp_npools_vr_lit3_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit3n), 0.d0)
        endif

        if (ispec_som1c > 0) then
          decomp_cpools_vr_som1_pf_loc(local_id) = max(xx_p(offsetim + ispec_som1c), 0.d0)
        endif

        if (ispec_som2c > 0) then
          decomp_cpools_vr_som2_pf_loc(local_id) = max(xx_p(offsetim + ispec_som2c), 0.d0)
        endif

        if (ispec_som3c > 0) then
          decomp_cpools_vr_som3_pf_loc(local_id) = max(xx_p(offsetim + ispec_som3c), 0.d0)
        endif

        if (ispec_som4c > 0) then
          decomp_cpools_vr_som4_pf_loc(local_id) = max(xx_p(offsetim + ispec_som4c), 0.d0)
        endif

        if (ispec_som1n > 0) then
          decomp_npools_vr_som1_pf_loc(local_id) = max(xx_p(offsetim + ispec_som1n), 0.d0)
        endif

        if (ispec_som2n > 0) then
          decomp_npools_vr_som2_pf_loc(local_id) = max(xx_p(offsetim + ispec_som2n), 0.d0)
        endif

        if (ispec_som3n > 0) then
          decomp_npools_vr_som3_pf_loc(local_id) = max(xx_p(offsetim + ispec_som3n), 0.d0)
        endif

        if (ispec_som4n > 0) then
          decomp_npools_vr_som4_pf_loc(local_id) = max(xx_p(offsetim + ispec_som4n), 0.d0)
        endif

        if(ispec_nh4 > 0) then
!           conc = xx_p(offset + ispec_nh4) * theta * 1000.0d0

            ! the following approach appears more like what output module does in pflotran
            ! but needs further checking if it efficient as directly read from 'xx_p' as above
            conc = rt_auxvar%pri_molal(ispec_nh4) * theta * 1000.0d0
            smin_nh4_vr_pf_loc(local_id) = max(conc, 0.d0)

            if (associated(rt_auxvar%total_sorb_eq)) then    ! equilibrium-sorption reactions used
                conc = rt_auxvar%total_sorb_eq(ispec_nh4)
                smin_nh4sorb_vr_pf_loc(local_id) = max(conc, 0.d0)
            else if (ispec_nh4sorb>0) then    ! kinetic-languir adsorption reaction used for soil NH4+ absorption
                conc = xx_p(offsetim + ispec_nh4sorb)                 ! unit: M (molC/m3)
                smin_nh4sorb_vr_pf_loc(local_id) = max(conc, 0.d0)
            endif

        endif

        if(ispec_no3 > 0) then
           conc = xx_p(offset + ispec_no3) * theta * 1000.0d0
           smin_no3_vr_pf_loc(local_id)   = max(conc, 0.d0)
        endif

        ! immobile gas conc in mol/m3 bulk soil to aovid 'theta' inconsistence (due to porosity) during unit conversion
        if(ispec_co2 > 0) then
           conc = xx_p(offsetim + ispec_co2)                    ! unit: M (molC/m3)
           gco2_vr_pf_loc(local_id)   = max(conc, 0.d0)
        endif

        if(ispec_n2 > 0) then
           conc = xx_p(offsetim + ispec_n2)                     ! unit: M (molN2/m3)
           gn2_vr_pf_loc(local_id)   = max(conc, 0.d0)
        endif

        if(ispec_n2o > 0) then
           conc = xx_p(offsetim + ispec_n2o)                    ! unit: M (molN2O/m3)
           gn2o_vr_pf_loc(local_id)   = max(conc, 0.d0)
        endif

        ! tracking N bgc reaction fluxes

        if (ispec_plantnh4uptake > 0) then
           conc = xx_p(offsetim + ispec_plantnh4uptake)
           accextrnh4_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step only
           xx_p(offsetim + ispec_plantnh4uptake) = zeroing_conc
        endif

        if (ispec_plantno3uptake > 0) then
           conc = xx_p(offsetim + ispec_plantno3uptake)
           accextrno3_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step only
           xx_p(offsetim + ispec_plantno3uptake) = zeroing_conc
        endif

        if(ispec_hr > 0) then
           conc = xx_p(offsetim + ispec_hr)
           acchr_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step
           xx_p(offsetim + ispec_hr) = zeroing_conc
        endif

        if(ispec_nmin > 0) then
           conc = xx_p(offsetim + ispec_nmin)
           accnmin_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step
           xx_p(offsetim + ispec_nmin) = zeroing_conc
        endif

        if(ispec_nimmp > 0) then
           conc = xx_p(offsetim + ispec_nimmp)
           accnimmp_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step
           xx_p(offsetim + ispec_nimmp) = zeroing_conc

        endif

        if(ispec_nimm > 0) then
           conc = xx_p(offsetim + ispec_nimm)
           accnimm_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step
           xx_p(offsetim + ispec_nimm) = zeroing_conc

        endif

        if(ispec_ngasmin > 0) then
           conc = xx_p(offsetim + ispec_ngasmin)
           accngasmin_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step
           xx_p(offsetim + ispec_ngasmin) = zeroing_conc
        endif

        if(ispec_ngasnitr > 0) then
           conc = xx_p(offsetim + ispec_ngasnitr)
           accngasnitr_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step
           xx_p(offsetim + ispec_ngasnitr) = zeroing_conc
        endif

        if(ispec_ngasdeni > 0) then
           conc = xx_p(offsetim + ispec_ngasdeni)
           accngasdeni_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step
           xx_p(offsetim + ispec_ngasdeni) = zeroing_conc
        endif

    enddo

    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit1_pfp, decomp_cpools_vr_lit1_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit2_pfp, decomp_cpools_vr_lit2_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit3_pfp, decomp_cpools_vr_lit3_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_cwd_pfp,  decomp_cpools_vr_cwd_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som1_pfp, decomp_cpools_vr_som1_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som2_pfp, decomp_cpools_vr_som2_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som3_pfp, decomp_cpools_vr_som3_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som4_pfp, decomp_cpools_vr_som4_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit1_pfp, decomp_npools_vr_lit1_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit2_pfp, decomp_npools_vr_lit2_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit3_pfp, decomp_npools_vr_lit3_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_cwd_pfp,  decomp_npools_vr_cwd_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_som1_pfp, decomp_npools_vr_som1_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_som2_pfp, decomp_npools_vr_som2_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_som3_pfp, decomp_npools_vr_som3_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_som4_pfp, decomp_npools_vr_som4_pf_loc, ierr)
    CHKERRQ(ierr)
    !
    call VecRestoreArrayF90(clm_pf_idata%smin_no3_vr_pfp, smin_no3_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%smin_nh4_vr_pfp, smin_nh4_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%smin_nh4sorb_vr_pfp, smin_nh4sorb_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    !
    call VecRestoreArrayF90(clm_pf_idata%gco2_vr_pfp, gco2_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%gn2_vr_pfp, gn2_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%gn2o_vr_pfp, gn2o_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    !
    call VecRestoreArrayF90(clm_pf_idata%accextrnh4_vr_pfp, accextrnh4_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%accextrno3_vr_pfp, accextrno3_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%acchr_vr_pfp, acchr_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%accnmin_vr_pfp, accnmin_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%accnimmp_vr_pfp, accnimmp_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%accnimm_vr_pfp, accnimm_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%accngasmin_vr_pfp, accngasmin_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%accngasnitr_vr_pfp, accngasnitr_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%accngasdeni_vr_pfp, accngasdeni_vr_pf_loc, ierr)
    CHKERRQ(ierr)
    !
    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%porosity0, porosity_loc_p, ierr)
    CHKERRQ(ierr)

    ! resetting the tracked variable states
    call DiscretizationGlobalToLocal(realization%discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)
    call VecCopy(field%tran_xx,field%tran_yy,ierr)
    CHKERRQ(ierr)
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)


   ! (iv) pass the 'pf' vecs to 'clm' vecs, which then can be passed to CLMCN (implemented in 'clm_pflotran_interfaceMod'

    if (ispec_lit1c > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_lit1_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_lit1_clms)
    endif

    if (ispec_lit2c > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_lit2_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_lit2_clms)
    endif

    if (ispec_lit3c > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_lit3_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_lit3_clms)
    endif

    if (ispec_cwdc > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_cwd_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_cwd_clms)
    endif

    if (ispec_som1c > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som1_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_som1_clms)
    endif

    if (ispec_som2c > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som2_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_som2_clms)
    endif

    if (ispec_som3c > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som3_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_som3_clms)
    endif

    if (ispec_som4c > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som4_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_som4_clms)
    endif

    if (ispec_lit1n > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_lit1_pfp, &
                                    clm_pf_idata%decomp_npools_vr_lit1_clms)
    endif

    if (ispec_lit2n > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_lit2_pfp, &
                                    clm_pf_idata%decomp_npools_vr_lit2_clms)
    endif

    if (ispec_lit3n > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_lit3_pfp, &
                                    clm_pf_idata%decomp_npools_vr_lit3_clms)
    endif

    if (ispec_cwdn > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_cwd_pfp, &
                                    clm_pf_idata%decomp_npools_vr_cwd_clms)
    endif

    if (ispec_som1n > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_som1_pfp, &
                                    clm_pf_idata%decomp_npools_vr_som1_clms)
    endif

    if (ispec_som2n > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_som2_pfp, &
                                    clm_pf_idata%decomp_npools_vr_som2_clms)
    endif

    if (ispec_som3n > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_som3_pfp, &
                                    clm_pf_idata%decomp_npools_vr_som3_clms)
    endif

    if (ispec_som4n > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_som4_pfp, &
                                    clm_pf_idata%decomp_npools_vr_som4_clms)
    endif

    if (ispec_co2 > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%gco2_vr_pfp, &
                                    clm_pf_idata%gco2_vr_clms)
    endif

    if(ispec_no3 > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%smin_no3_vr_pfp, &
                                    clm_pf_idata%smin_no3_vr_clms)
    endif

    if(ispec_nh4 > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%smin_nh4_vr_pfp, &
                                    clm_pf_idata%smin_nh4_vr_clms)

      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%smin_nh4sorb_vr_pfp, &
                                    clm_pf_idata%smin_nh4sorb_vr_clms)
    endif

    if(ispec_n2 > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%gn2_vr_pfp, &
                                    clm_pf_idata%gn2_vr_clms)
    endif

    if(ispec_n2o > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%gn2o_vr_pfp, &
                                    clm_pf_idata%gn2o_vr_clms)
    endif

    if(ispec_nmin > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%accnmin_vr_pfp, &
                                    clm_pf_idata%accnmin_vr_clms)
    endif

    if(ispec_plantnh4uptake > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%accextrnh4_vr_pfp, &
                                    clm_pf_idata%accextrnh4_vr_clms)
    endif

    if(ispec_plantno3uptake > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%accextrno3_vr_pfp, &
                                    clm_pf_idata%accextrno3_vr_clms)
    endif

    if(ispec_hr > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%acchr_vr_pfp, &
                                    clm_pf_idata%acchr_vr_clms)
    endif

    if(ispec_nimmp > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%accnimmp_vr_pfp, &
                                    clm_pf_idata%accnimmp_vr_clms)
    endif

    if(ispec_nimm > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%accnimm_vr_pfp, &
                                    clm_pf_idata%accnimm_vr_clms)
    endif

    if(ispec_ngasmin > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%accngasmin_vr_pfp, &
                                    clm_pf_idata%accngasmin_vr_clms)
    endif

    if(ispec_ngasnitr > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%accngasnitr_vr_pfp, &
                                    clm_pf_idata%accngasnitr_vr_clms)
    endif

    if(ispec_ngasdeni > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%accngasdeni_vr_pfp, &
                                    clm_pf_idata%accngasdeni_vr_clms)
    endif

  end subroutine pflotranModelGetBgcVariablesFromPF


! ************************************************************************** !

  subroutine pflotranModelGetBCMassBalanceDeltaFromPF(pflotran_model)
  !
  ! Calculate mass balance at BC for passing flow rates to CLM
  !
  ! Author: Fengming Yuan
  ! Date: 03/14/2014
  !
    use Realization_class
    use Patch_module
    use Grid_module
    use Option_module
    use Connection_module
    use Coupler_module
    use Utility_module
    use String_module

    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type

    use Global_Aux_module
    use Reactive_Transport_Aux_module
    use Reaction_Aux_module

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer :: pflotran_model
    class(realization_type), pointer :: realization

    type(option_type), pointer :: option
    type(patch_type), pointer :: patch
    type(grid_type), pointer :: grid

    type(coupler_type), pointer :: boundary_condition
    type(connection_set_type), pointer :: cur_connection_set
    type(global_auxvar_type), pointer :: global_auxvars_bc(:)
    type(reactive_transport_auxvar_type), pointer :: rt_auxvars_bc(:)

    PetscReal, pointer :: qinfl_subsurf_pf_loc(:)
    PetscReal, pointer :: qsurf_subsurf_pf_loc(:)
    PetscReal, pointer :: qflux_subbase_pf_loc(:)
    PetscReal, pointer :: f_nh4_subsurf_pf_loc(:)
    PetscReal, pointer :: f_no3_subsurf_pf_loc(:)
    PetscReal, pointer :: f_nh4_subbase_pf_loc(:)
    PetscReal, pointer :: f_no3_subbase_pf_loc(:)

    PetscInt :: local_id, ghosted_id, iconn
    PetscInt :: offset
    PetscErrorCode :: ierr

    !------------------------------------------------------------------------------------

    if (clm_pf_idata%nlpf_2dtop <= 0 .and. clm_pf_idata%ngpf_2dtop <= 0    &
        .and. clm_pf_idata%nlpf_2dbot <= 0 .and. clm_pf_idata%ngpf_2dbot <= 0) then
        return
    endif

    !
    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer &
          = "ERROR: GetMassBalance not supported under this mode.."
         call printErrMsg(pflotran_model%option)
    end select

    patch => realization%patch
    grid => patch%grid
    option => realization%option
    !
    call VecGetArrayF90(clm_pf_idata%qinfl_subsurf_pfp, qinfl_subsurf_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%qsurf_subsurf_pfp, qsurf_subsurf_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%qflux_subbase_pfp, qflux_subbase_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%f_nh4_subsurf_pfp, f_nh4_subsurf_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%f_nh4_subbase_pfp, f_nh4_subbase_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%f_no3_subsurf_pfp, f_no3_subsurf_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%f_no3_subbase_pfp, f_no3_subbase_pf_loc, ierr)
    CHKERRQ(ierr)

    qinfl_subsurf_pf_loc(:) = 0.d0
    qsurf_subsurf_pf_loc(:) = 0.d0
    qflux_subbase_pf_loc(:) = 0.d0

    !
    boundary_condition => patch%boundary_condition_list%first
    global_auxvars_bc => patch%aux%Global%auxvars_bc
    if (option%ntrandof > 0) then
       rt_auxvars_bc => patch%aux%RT%auxvars_bc
    endif

    do
      if (.not.associated(boundary_condition)) exit

      cur_connection_set => boundary_condition%connection_set

      offset = cur_connection_set%offset

      if (option%nflowdof > 0) then

          ! retrieving H2O flux at top BC
          do iconn = 1, cur_connection_set%num_connections
             local_id = cur_connection_set%id_dn(iconn)
             ghosted_id = grid%nL2G(local_id)
             if (patch%imat(ghosted_id) <= 0) cycle

             if(StringCompare(boundary_condition%name,'clm_gpress_bc')) then          ! infilitration (+)
                qinfl_subsurf_pf_loc(iconn) = &
                               -global_auxvars_bc(offset+iconn)%mass_balance(1,1)

                ! 'mass_balance' IS accumulative, so need to reset to Zero for next desired time-step
                global_auxvars_bc(offset+iconn)%mass_balance(1,1) = 0.d0

             endif

             if(StringCompare(boundary_condition%name,'clm_gflux_overflow')) then    ! surface overflow (-)
                qsurf_subsurf_pf_loc(iconn) = &
                               -global_auxvars_bc(offset+iconn)%mass_balance(1,1)

                 ! 'mass_balance' IS accumulative, so need to reset to Zero for next desired time-step
                global_auxvars_bc(offset+iconn)%mass_balance(1,1) = 0.d0

             endif

             ! retrieving H2O flux at bottom BC
             if(StringCompare(boundary_condition%name,'clm_bflux_bc')) then          ! bottom water flux
                qflux_subbase_pf_loc(iconn) = &
                               -global_auxvars_bc(offset+iconn)%mass_balance(1,1)

                ! 'mass_balance' IS accumulative, so need to reset to Zero for next desired time-step
                global_auxvars_bc(offset+iconn)%mass_balance(1,1) = 0.d0

             endif

          enddo


      endif

! (TODO) the following needs checking
!    if (option%ntrandof > 0) then
!
!      ! boundary chemical flux
!      do iconn = 1, boundary_condition%connection_set%num_connections
!        sum_mol = rt_aux_vars_bc(offset+iconn)%mass_balance_delta
!      enddo
!
!    endif

       boundary_condition => boundary_condition%next

    enddo

    call VecRestoreArrayF90(clm_pf_idata%qinfl_subsurf_pfp, qinfl_subsurf_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%qsurf_subsurf_pfp, qsurf_subsurf_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%qflux_subbase_pfp, qflux_subbase_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%f_nh4_subsurf_pfp, f_nh4_subsurf_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%f_nh4_subbase_pfp, f_nh4_subbase_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%f_no3_subsurf_pfp, f_no3_subsurf_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%f_no3_subbase_pfp, f_no3_subbase_pf_loc, ierr)
    CHKERRQ(ierr)

    ! pass vecs to CLM
    if (clm_pf_idata%nlpf_2dtop > 0 .and. clm_pf_idata%ngpf_2dtop > 0 ) then
      call MappingSourceToDestination(pflotran_model%map_pf_2dtop_to_clm_2dtop, &
                                    pflotran_model%option, &
                                    clm_pf_idata%qinfl_subsurf_pfp, &
                                    clm_pf_idata%qinfl_subsurf_clms)

      call MappingSourceToDestination(pflotran_model%map_pf_2dtop_to_clm_2dtop, &
                                    pflotran_model%option, &
                                    clm_pf_idata%qsurf_subsurf_pfp, &
                                    clm_pf_idata%qsurf_subsurf_clms)
    endif

    if (clm_pf_idata%nlpf_2dbot > 0 .and. clm_pf_idata%ngpf_2dbot > 0 ) then
      call MappingSourceToDestination(pflotran_model%map_pf_2dbot_to_clm_2dbot, &
                                    pflotran_model%option, &
                                    clm_pf_idata%qflux_subbase_pfp, &
                                    clm_pf_idata%qflux_subbase_clms)
    endif
  end subroutine pflotranModelGetBCMassBalanceDeltaFromPF

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  
end module pflotran_clm_main_module

