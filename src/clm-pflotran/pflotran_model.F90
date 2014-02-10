module pflotran_model_module

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

  ! Note:
  !
  ! CLM has the following:
  !   (i) 3D subsurface grid (CLM_SUB);
  !   (ii) 2D surface grid (CLM_SRF).
  ! CLM decomposes the 3D subsurface grid across processors in a 2D (i.e.
  ! cells in Z are not split across processors). Thus, the surface cells of
  ! 3D subsurface grid are on the same processors as the 2D surface grid.
  !
  ! PFLOTRAN has the following:
  !   (i) 3D subsurface grid (PF_SUB);
  !   (ii) surface control volumes of 3D subsurface grid (PF_2DSUB);
  !   (iii) 2D surface grid (PF_SRF).
  ! In PFLOTRAN, control volumes in PF_2DSUB and PF_SRF may reside on different
  ! processors. PF_SUB and PF_2DSUB are derived from simulation%realization;
  ! while PF_SRF refers to simulation%surf_realization.

  ! map level constants
  PetscInt, parameter, public :: CLM_SUB_TO_PF_SUB           = 1 ! 3D --> 3D
  PetscInt, parameter, public :: CLM_SUB_TO_PF_EXTENDED_SUB  = 2 ! 3D --> extended 3D
  PetscInt, parameter, public :: CLM_SRF_TO_PF_2DSUB         = 3 ! 2D --> SURF of 3D grid
  PetscInt, parameter, public :: CLM_SRF_TO_PF_SRF           = 4 ! 2D --> 2D SURF grid
  PetscInt, parameter, public :: PF_SUB_TO_CLM_SUB           = 5 ! 3D --> 3D
  PetscInt, parameter, public :: PF_SRF_TO_CLM_SRF           = 6 ! 2D SURF grid --> 2D

  ! mesh ids
  PetscInt, parameter, public :: CLM_SUB_MESH   = 1
  PetscInt, parameter, public :: CLM_SRF_MESH   = 2
  PetscInt, parameter, public :: PF_SUB_MESH    = 3
  PetscInt, parameter, public :: PF_2DSUB_MESH  = 4
  PetscInt, parameter, public :: PF_SRF_MESH    = 5

  type, public :: inside_each_overlapped_cell
     PetscInt           :: id
     PetscInt           :: ocell_count
     PetscInt,  pointer :: ocell_id(:)
     PetscReal, pointer :: perc_vol_overlap(:)
     PetscReal          :: total_vol_overlap
  end type inside_each_overlapped_cell

  type, public :: pflotran_model_type
    class(simulation_base_type),  pointer :: simulation
    type(multi_simulation_type), pointer :: multisimulation
    type(option_type),      pointer :: option
    PetscReal :: pause_time_1
    PetscReal :: pause_time_2
    type(inside_each_overlapped_cell), pointer :: pf_cells(:)
    type(inside_each_overlapped_cell), pointer :: clm_cells(:)
    type(mapping_type),                pointer :: map_clm_sub_to_pf_sub
    type(mapping_type),                pointer :: map_clm_sub_to_pf_extended_sub
    type(mapping_type),                pointer :: map_clm_srf_to_pf_2dsub
    type(mapping_type),                pointer :: map_clm_srf_to_pf_srf
    type(mapping_type),                pointer :: map_pf_sub_to_clm_sub
    type(mapping_type),                pointer :: map_pf_srf_to_clm_srf
     
    Vec :: hksat_x_clm
    Vec :: hksat_y_clm
    Vec :: hksat_z_clm
    Vec :: sucsat_clm
    Vec :: watsat_clm
    Vec :: bsw_clm
    Vec :: qflx_clm
    Vec :: sat_clm

    Vec :: hksat_x_pf
    Vec :: hksat_y_pf
    Vec :: hksat_z_pf
    Vec :: sucsat_pf
    Vec :: watsat_pf
    Vec :: bsw_pf
    Vec :: qflx_pf
    Vec :: sat_pf

    PetscInt :: nlclm
    PetscInt :: ngclm

    PetscLogDouble :: timex(4), timex_wall(4)

    PetscBool :: b_out_bgc_rate

  end type pflotran_model_type

  public::pflotranModelCreate,               &
       pflotranModelInitMapping,             &
       pflotranModelSetSoilProp,             &
       pflotranModelSetICs,                  &
       pflotranModelUpdateFlowConds,         &
       pflotranModelGetUpdatedStates,        &
       pflotranModelStepperRunInit,          &
       pflotranModelStepperRunTillPauseTime, &
       pflotranModelInsertWaypoint,          &
       pflotranModelDeleteWaypoint,          &
       pflotranModelSetupRestart,            &
       pflotranModelStepperRunFinalize,      &
       pflotranModelStepperCheckpoint,       &
       pflotranModelNSurfCells3DDomain,      &
       pflotranModelGetTopFaceArea,          &
       pflotranModelDestroy

  private :: &
       pflotranModelSetupMappingFiles

  public ::  &
       pflotranModelGetSoilProp,               &
       pflotranModelSetInitialTHStatesfromCLM, &    ! initializing TH from CLM to PFLOTRAN flow mode
       pflotranModelSetInitialConcentrations,  &
       pflotranModelUpdateTHfromCLM,           &     ! dynamically update TH from CLM to drive PFLOTRAN BGC
       pflotranModelUpdateO2fromCLM,           &
       pflotranModelSetBGCRates,               &
       pflotranModelGetBgcVariables,           &
       pflotranModelSetSoilHbcs


contains

! ************************************************************************** !

  function pflotranModelCreate(mpicomm, pflotran_prefix)
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
    use PFLOTRAN_Factory_module
    use Subsurface_Factory_module, only : SubsurfaceInitialize
    use Hydrogeophysics_Factory_module
    use Surface_Factory_module
    use Surf_Subsurf_Factory_module
  
    implicit none

#include "finclude/petscsys.h"

    PetscInt, intent(in) :: mpicomm
    character(len=256), intent(in) :: pflotran_prefix

    type(pflotran_model_type), pointer :: pflotranModelCreate

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
    call PFLOTRANInitializePostPetsc(model%multisimulation, model%option)

    ! NOTE(bja, 2013-07-19) GB's Hack to get communicator correctly
    ! setup on mpich/mac. should be generally ok, but may need an
    ! apple/mpich ifdef if it cause problems elsewhere.
    PETSC_COMM_SELF = MPI_COMM_SELF

    ! TODO(bja, 2013-07-15) this needs to be left alone for pflotran
    ! to deal with, or we need a valid unit number from the driver as
    ! a function parameter.
!!$    model%option%fid_out = 16

    model%pause_time_1 = -1.0d0
    model%pause_time_2 = -1.0d0

    ! FIXME(bja, 2013-07-17) hard code subsurface for now....
    !model%option%simulation_mode = 'SURFACE_SUBSURFACE'
    model%option%simulation_mode = 'SUBSURFACE'
    select case(model%option%simulation_mode)
      case('SUBSURFACE')
         call SubsurfaceInitialize(model%simulation, model%option)
      case('HYDROGEOPHYSICS')
         call HydrogeophysicsInitialize(model%simulation, model%option)
      case('SURFACE')
         call SurfaceInitialize(model%simulation, model%option)
      case('SURFACE_SUBSURFACE')
         call SurfSubsurfaceInitialize(model%simulation, model%option)
      case default
         model%option%io_buffer = 'Simulation Mode not recognized : ' // model%option%simulation_mode
         call printErrMsg(model%option)
    end select
    ! NOTE(bja, 2013-07-15) needs to go before InitializeRun()...?
    call pflotranModelSetupMappingFiles(model)

!  for debug only, to be removed or replaced with h5 file
    model%b_out_bgc_rate = PETSC_FALSE
    if (model%b_out_bgc_rate) then
        open(unit=100,file="bgc_rate_clm_to_pflotran.txt",form="formatted",status="replace")
        write(100, *) '# Lit1C Lit2C Lit3C Lit1N Lit2N Lit3N mineralN plantN leachedN DenitrifiedN'
    endif

    pflotranModelCreate => model

  end function pflotranModelCreate

! ************************************************************************** !

  subroutine pflotranModelSetupMappingFiles(model)
  ! 
  ! pflotranModelSetupMappingFiles
  ! create the mapping objects, reopen the input file and read the file names
  ! before model integration is performed by the call to StepperRun()
  ! routine
  ! NOTE(bja, 2013-07) this really needs to be moved out of pflotran
  ! CLM should be responsible for passing data in the correct
  ! format. That may require pflotran to provide a call back function
  ! for grid info.
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! 

    use String_module
    use Option_module
    use Input_Aux_module
    use Mapping_module

    implicit none

#include "finclude/petscsys.h"

    type(pflotran_model_type), pointer, intent(inout) :: model
    type(input_type), pointer :: input

    PetscBool :: clm2pf_flux_file
    PetscBool :: clm2pf_soil_file
    PetscBool :: clm2pf_gflux_file
    PetscBool :: clm2pf_rflux_file
    PetscBool :: pf2clm_flux_file
    PetscBool :: pf2clm_surf_file
    character(len=MAXSTRINGLENGTH) :: string
    character(len=MAXWORDLENGTH) :: word

    nullify(model%pf_cells)
    nullify(model%clm_cells)
    nullify(model%map_clm_sub_to_pf_sub)
    nullify(model%map_clm_sub_to_pf_extended_sub)
    nullify(model%map_clm_srf_to_pf_2dsub)
    nullify(model%map_clm_srf_to_pf_srf)
    nullify(model%map_pf_sub_to_clm_sub)
    nullify(model%map_pf_srf_to_clm_srf)

    model%map_clm_sub_to_pf_sub          => MappingCreate()
    model%map_clm_sub_to_pf_extended_sub => MappingCreate()
    model%map_clm_srf_to_pf_2dsub        => MappingCreate()
    model%map_clm_srf_to_pf_srf          => MappingCreate()
    model%map_pf_sub_to_clm_sub          => MappingCreate()
    model%map_pf_srf_to_clm_srf          => MappingCreate()

    model%nlclm = -1
    model%ngclm = -1

    input => InputCreate(15, &
                    model%option%input_filename, model%option)

    ! Read names of mapping file
    clm2pf_flux_file=PETSC_FALSE
    clm2pf_soil_file=PETSC_FALSE
    clm2pf_gflux_file=PETSC_FALSE
    clm2pf_rflux_file=PETSC_FALSE
    pf2clm_flux_file=PETSC_FALSE
    pf2clm_surf_file=PETSC_FALSE
    
    string = "MAPPING_FILES"
    call InputFindStringInFile(input,model%option,string)

    do
      call InputReadPflotranString(input, model%option)
      if (InputCheckExit(input, model%option)) exit
      if (input%ierr /= 0) exit

      call InputReadWord(input, model%option, word, PETSC_TRUE)
      call InputErrorMsg(input, model%option, 'keyword', 'MAPPING_FILES')
      call StringToUpper(word)

      select case(trim(word))
        case('CLM2PF_FLUX_FILE')
          call InputReadNChars(input, model%option, model%map_clm_sub_to_pf_sub%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_clm_sub_to_pf_sub%filename = trim(model%map_clm_sub_to_pf_sub%filename)//CHAR(0)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')   
          clm2pf_flux_file=PETSC_TRUE
        case('CLM2PF_SOIL_FILE')
          call InputReadNChars(input, model%option, model%map_clm_sub_to_pf_extended_sub%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')   
          clm2pf_soil_file=PETSC_TRUE
        case('CLM2PF_GFLUX_FILE')
          call InputReadNChars(input, model%option, model%map_clm_srf_to_pf_2dsub%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          clm2pf_gflux_file=PETSC_TRUE
        case('CLM2PF_RFLUX_FILE')
          call InputReadNChars(input, model%option, model%map_clm_srf_to_pf_srf%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          clm2pf_rflux_file=PETSC_TRUE
        case('PF2CLM_SURF_FILE')
          call InputReadNChars(input, model%option, model%map_pf_srf_to_clm_srf%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          pf2clm_surf_file=PETSC_TRUE
        case('PF2CLM_FLUX_FILE')
          call InputReadNChars(input, model%option, model%map_pf_sub_to_clm_sub%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')   
          pf2clm_flux_file=PETSC_TRUE
        case default
          model%option%io_buffer='Keyword ' // trim(word) // &
            ' in input file not recognized'
          call printErrMsg(model%option)
      end select

    enddo
    call InputDestroy(input)

    if ((.not. clm2pf_soil_file) .or. (.not. clm2pf_flux_file) .or. &
        (.not. pf2clm_flux_file) ) then
      model%option%io_buffer='One of the mapping files not found'
      call printErrMsg(model%option)
    endif
    
    if(model%option%iflowmode==TH_MODE.and.(.not.clm2pf_gflux_file)) then
      model%option%io_buffer='Running in TH_MODE without a CLM2PF_GFLUX_FILE'
      call printErrMsg(model%option)
    endif

#ifdef SURFACE_FLOW
    if( (model%option%nsurfflowdof>0)) then
       if ((.not. clm2pf_rflux_file)) then
        model%option%io_buffer='Running in surface flow without a ' // &
          'CLM2PF_RFLUX_FILE'
        call printErrMsg(model%option)
       endif
       if ((.not. pf2clm_surf_file)) then
        model%option%io_buffer='Running in surface flow without a ' // &
          'PF2CLM_SURF_FILE'
        call printErrMsg(model%option)
       endif
    endif
#endif
  end subroutine pflotranModelSetupMappingFiles

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

subroutine pflotranModelSetICs(pflotran_model)
  ! 
  ! Set initial conditions
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/22/2010
  ! 

    use Realization_class
    use Patch_module
    use Grid_module
    use Richards_Aux_module
    use Field_module
    use clm_pflotran_interface_data
    use Global_Aux_module
    use Discretization_module
    use Richards_module
    use TH_module
    use Option_module

    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer           :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(richards_auxvar_type), pointer       :: aux_var
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    type(simulation_base_type), pointer :: simulation

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal          :: den, vis, grav
    PetscReal, pointer :: xx_loc_p(:)

    PetscScalar, pointer :: press_pf_loc(:) ! Pressure [Pa]

    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: pflotranModelSetICs only works on subsurface simulations."
         call printErrMsg(pflotran_model%option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field
    global_aux_vars  => patch%aux%Global%aux_vars

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_clm, &
                                    clm_pf_idata%press_pf)

    if (pflotran_model%option%iflowmode .ne. RICHARDS_MODE .and. &
        pflotran_model%option%iflowmode .ne. TH_MODE) then
        pflotran_model%option%io_buffer='pflotranModelSetICs ' // &
          'not implmented for this mode.'
        call printErrMsg(pflotran_model%option)
    endif

    call VecGetArrayF90(field%flow_xx, xx_loc_p, ierr)
    call VecGetArrayF90(clm_pf_idata%press_pf, press_pf_loc, ierr)

    do local_id = 1, grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
       xx_loc_p(ghosted_id)=press_pf_loc(local_id)
    enddo

    call VecRestoreArrayF90(field%flow_xx, xx_loc_p, ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_pf, press_pf_loc, ierr)

    ! update dependent vectors: Saturation
    call DiscretizationGlobalToLocal(realization%discretization, field%flow_xx, &
         field%flow_xx_loc, NFLOWDOF)
    call VecCopy(field%flow_xx, field%flow_yy, ierr)

    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        call RichardsUpdateAuxVars(realization)
      case (TH_MODE)
        call THUpdateAuxVars(realization)
      case default
        pflotran_model%option%io_buffer='pflotranModelSetICs ' // &
          'not implmented for this mode.'
        call printErrMsg(pflotran_model%option)
    end select

end subroutine pflotranModelSetICs

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

    use Realization_class
    use Patch_module
    use Grid_module
    use Richards_Aux_module
    use TH_Aux_module
    use Field_module
    use Option_module

    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer           :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(richards_auxvar_type), pointer       :: rich_aux_vars(:)
    type(richards_auxvar_type), pointer       :: rich_aux_var
    type(th_auxvar_type), pointer             :: th_aux_vars(:)
    type(th_auxvar_type), pointer             :: th_aux_var
    type(simulation_base_type), pointer :: simulation

    PetscErrorCode     :: ierr
    PetscInt           :: local_id
    PetscReal          :: den, vis, grav
    PetscReal, pointer :: porosity_loc_p(:), vol_ovlap_arr(:)
    PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
    PetscReal          :: bc_lambda, bc_alpha

    PetscScalar, pointer :: hksat_x_pf_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_y_pf_loc(:) ! hydraulic conductivity in y-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_z_pf_loc(:) ! hydraulic conductivity in z-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: watsat_pf_loc(:)  ! minimum soil suction (mm)
    PetscScalar, pointer :: sucsat_pf_loc(:)  ! volumetric soil water at saturation (porosity)
    PetscScalar, pointer :: bsw_pf_loc(:)     ! Clapp and Hornberger "b"
!    PetscScalar, pointer :: bsw_clm_loc(:)    ! Clapp and Hornberger "b"

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav = 9.81d0       ! [m/S^2]

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
    grid            => patch%grid
    field           => realization%field

    select case(pflotran_model%option%iflowmode)
      case(RICHARDS_MODE)
        rich_aux_vars   => patch%aux%Richards%aux_vars
      case(TH_MODE)
        th_aux_vars   => patch%aux%TH%aux_vars
      case default
        !F.-M. Yuan: if no flowmode AND no reactive-transport, then let crash the run
        ! because 'uniform_velocity' with [0, 0, 0] velocity transport doesn't need specific flowmode,
        ! which is used for BGC coupling only (i.e., no TH coupling)
        if(pflotran_model%option%ntrandof.le.0) then
            call printErrMsg(pflotran_model%option, &
               'Current PFLOTRAN mode not supported by pflotranModelSetSoilProp')
        endif
    end select

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_extended_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%hksat_x_clm, &
                                    clm_pf_idata%hksat_x_pf)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_extended_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%hksat_y_clm, &
                                    clm_pf_idata%hksat_y_pf)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_extended_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%hksat_z_clm, &
                                    clm_pf_idata%hksat_z_pf)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_extended_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sucsat_clm, &
                                    clm_pf_idata%sucsat_pf)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_extended_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%bsw_clm, &
                                    clm_pf_idata%bsw_pf)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_extended_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%watsat_clm, &
                                    clm_pf_idata%watsat_pf)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_extended_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%watfc_clm, &
                                    clm_pf_idata%watfc_pf)
! for denitrification
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_extended_sub, &    ! 'extended' for 'physical' (?)
                                    pflotran_model%option, &
                                    clm_pf_idata%bulkdensity_dry_clm, &
                                    clm_pf_idata%bulkdensity_dry_pf)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &    ! 'non-extended' for 'bgc' (?)
                                    pflotran_model%option, &
                                    clm_pf_idata%cellorg_clm, &
                                    clm_pf_idata%cellorg_pf)
!    if(clm_pf_idata%use_lch4) then
!       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
!                                    pflotran_model%option, &
!                                    clm_pf_idata%o2_decomp_depth_unsat_clm, &
!                                    clm_pf_idata%o2_decomp_depth_unsat_pf)

!       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
!                                    pflotran_model%option, &
!                                    clm_pf_idata%conc_o2_unsat_clm, &
!                                    clm_pf_idata%conc_o2_unsat_pf)

!for some reason, o2_decomp_depth_sat is nan
!       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
!                                    pflotran_model%option, &
!                                    clm_pf_idata%o2_decomp_depth_sat_clm, &
!                                    clm_pf_idata%o2_decomp_depth_sat_pf)

!       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
!                                    pflotran_model%option, &
!                                    clm_pf_idata%conc_o2_sat_clm, &
!                                    clm_pf_idata%conc_o2_sat_pf)
!
!    endif

    call VecGetArrayF90(clm_pf_idata%hksat_x_pf, hksat_x_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_y_pf, hksat_y_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_z_pf, hksat_z_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%sucsat_pf,  sucsat_pf_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_pf,  watsat_pf_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%bsw_pf,     bsw_pf_loc,     ierr)
!    call VecGetArrayF90(clm_pf_idata%bsw_clm,    bsw_clm_loc,    ierr)

    call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)

    if(pflotran_model%option%iflowmode==RICHARDS_MODE .or. &
       pflotran_model%option%iflowmode==TH_MODE) then
       ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
         call VecGetArrayF90(field%perm_xx_loc,  perm_xx_loc_p,  ierr)
         call VecGetArrayF90(field%perm_yy_loc,  perm_yy_loc_p,  ierr)
         call VecGetArrayF90(field%perm_zz_loc,  perm_zz_loc_p,  ierr)
    endif

    do local_id = 1, grid%ngmax

      ! bc_alpha [1/Pa]; while sucsat [mm of H20]
      ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
      bc_alpha = 1.d0/(sucsat_pf_loc(local_id)*grav)

      ! bc_lambda = 1/bsw
      bc_lambda = 1.d0/bsw_pf_loc(local_id)
      
      select case(pflotran_model%option%iflowmode)
        case(RICHARDS_MODE)
          rich_aux_var => rich_aux_vars(local_id)
          rich_aux_var%bc_alpha = bc_alpha
          rich_aux_var%bc_lambda = bc_lambda
        case(TH_MODE)
          th_aux_var => th_aux_vars(local_id)
          th_aux_var%bc_alpha = min(bc_alpha,10.d-4)
          th_aux_var%bc_lambda = bc_lambda
      end select

      ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
      ! [m^2]          [mm/sec]
      if(pflotran_model%option%iflowmode==RICHARDS_MODE .or. &
         pflotran_model%option%iflowmode==TH_MODE) then
           ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
           perm_xx_loc_p(local_id) = hksat_x_pf_loc(local_id)*vis/(den*grav)/1000.d0
           perm_yy_loc_p(local_id) = hksat_y_pf_loc(local_id)*vis/(den*grav)/1000.d0
           perm_zz_loc_p(local_id) = hksat_z_pf_loc(local_id)*vis/(den*grav)/1000.d0
      endif

      porosity_loc_p(local_id) = watsat_pf_loc(local_id)

    enddo

    call VecRestoreArrayF90(clm_pf_idata%hksat_x_pf, hksat_x_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y_pf, hksat_y_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z_pf, hksat_z_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat_pf,  sucsat_pf_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_pf,  watsat_pf_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw_pf,     bsw_pf_loc,     ierr)
!    call VecRestoreArrayF90(clm_pf_idata%bsw_clm,    bsw_clm_loc,    ierr)

    call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
    if(pflotran_model%option%iflowmode==RICHARDS_MODE .or. &
       pflotran_model%option%iflowmode==TH_MODE) then
           ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
        call VecRestoreArrayF90(field%perm_xx_loc,  perm_xx_loc_p,  ierr)
        call VecRestoreArrayF90(field%perm_yy_loc,  perm_yy_loc_p,  ierr)
        call VecRestoreArrayF90(field%perm_zz_loc,  perm_zz_loc_p,  ierr)
    endif

  end subroutine pflotranModelSetSoilProp

! ************************************************************************** !

  subroutine pflotranModelInitMapping(pflotran_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id)
  ! 
  ! #endif
  ! Initialize mapping between the two model grid
  ! (CLM and PFLTORAN)
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/24/2011
  ! 

    use Input_Aux_module
    use Option_module
    use Realization_class
    use Grid_module
    use Patch_module
    use Coupler_module
    use Connection_module
    use String_module
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Mapping_module

    implicit none

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename
    
    select case (map_id)
      case (CLM_SUB_TO_PF_SUB, CLM_SUB_TO_PF_EXTENDED_SUB, PF_SUB_TO_CLM_SUB)
        call pflotranModelInitMappingSub2Sub(pflotran_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id)
      case (CLM_SRF_TO_PF_2DSUB)
        call pflotranModelInitMapSrfTo2DSub(pflotran_model,  &
                                            grid_clm_cell_ids_nindex, &
                                            grid_clm_npts_local, &
                                            map_id)
      case (CLM_SRF_TO_PF_SRF, PF_SRF_TO_CLM_SRF)
        call pflotranModelInitMapSrfToSrf(pflotran_model,  &
                                            grid_clm_cell_ids_nindex, &
                                            grid_clm_npts_local, &
                                            map_id)
      case default
        pflotran_model%option%io_buffer = 'Invalid map_id argument to ' // &
          'pflotranModelInitMapping'
        call printErrMsg(pflotran_model%option)
    end select

  end subroutine pflotranModelInitMapping

! ************************************************************************** !

  subroutine pflotranModelInitMappingSub2Sub(pflotran_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id)
  ! 
  ! Initialize mapping between 3D subsurface
  ! CLM grid and 3D subsurface PFLOTRAN grid.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/09/2013
  ! 

    use Input_Aux_module
    use Option_module
    use Realization_class
    use Grid_module
    use Patch_module
    use Coupler_module
    use Connection_module
    use String_module
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Mapping_module

    implicit none

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename

    ! local
    PetscInt                           :: local_id, grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_clm_npts_ghost, source_mesh_id
    PetscInt                           :: dest_mesh_id
    PetscInt, pointer                  :: grid_pf_cell_ids_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_local_nindex(:)

    type(mapping_type), pointer        :: map
    type(option_type), pointer         :: option
    class(realization_type), pointer    :: realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch
    type(coupler_type), pointer        :: boundary_condition
    type(connection_set_type), pointer :: cur_connection_set

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
    option          => pflotran_model%option
    patch           => realization%patch
    grid            => patch%grid

    !
    ! Mapping to/from entire PFLOTRAN 3D subsurface domain
    !

    ! Choose the appriopriate map
    select case(map_id)
      case(CLM_SUB_TO_PF_SUB)
        map => pflotran_model%map_clm_sub_to_pf_sub
        source_mesh_id = CLM_SUB_MESH
        dest_mesh_id = PF_SUB_MESH
      case(CLM_SUB_TO_PF_EXTENDED_SUB)
        map => pflotran_model%map_clm_sub_to_pf_extended_sub
        source_mesh_id = CLM_SUB_MESH
        dest_mesh_id = PF_SUB_MESH
      case(PF_SUB_TO_CLM_SUB)
        map => pflotran_model%map_pf_sub_to_clm_sub
        source_mesh_id = PF_SUB_MESH
      case default
        option%io_buffer = 'Invalid map_id argument to pflotranModelInitMapping'
        call printErrMsg(option)
    end select

    ! Read mapping file
    if (index(map%filename, '.h5') > 0) then
      call MappingReadHDF5(map, map%filename, option)
    else
      call MappingReadTxtFile(map, map%filename, option)
    endif

    grid_clm_npts_ghost=0

    ! Allocate memory to identify if CLM cells are local or ghosted.
    ! Note: Presently all CLM cells are local
    allocate(grid_clm_local_nindex(grid_clm_npts_local))
    do local_id = 1, grid_clm_npts_local
      grid_clm_local_nindex(local_id) = 1 ! LOCAL
    enddo

    ! Find cell IDs for PFLOTRAN grid
    grid_pf_npts_local = grid%nlmax
    grid_pf_npts_ghost = grid%ngmax - grid%nlmax

    allocate(grid_pf_cell_ids_nindex(grid%ngmax))
    do local_id = 1, grid%ngmax
      grid_pf_cell_ids_nindex(local_id) = grid%nG2A(local_id)-1
    enddo

    allocate(grid_pf_local_nindex(grid%ngmax))
    do local_id = 1, grid%ngmax
      if (grid%nG2L(local_id) == 0) then
        grid_pf_local_nindex(local_id) = 0 ! GHOST
      else
        grid_pf_local_nindex(local_id) = 1 ! LOCAL
      endif
    enddo

    select case(source_mesh_id)
      case(CLM_SUB_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_cell_ids_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_SUB_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_pf_npts_local, &
                                        grid_pf_cell_ids_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_clm_npts_local, &
                                              grid_clm_npts_ghost, &
                                              grid_clm_cell_ids_nindex, &
                                              grid_clm_local_nindex)
      case default
        option%io_buffer = 'Invalid argument source_mesh_id passed to pflotranModelInitMapping'
        call printErrMsg(option)
    end select

    deallocate(grid_pf_cell_ids_nindex)
    deallocate(grid_pf_local_nindex)

    call MappingDecompose(map, option)
    call MappingFindDistinctSourceMeshCellIds(map, option)
    call MappingCreateWeightMatrix(map, option)
    call MappingCreateScatterOfSourceMesh(map, option)

  end subroutine pflotranModelInitMappingSub2Sub

! ************************************************************************** !

  subroutine pflotranModelInitMapSrfTo2DSub(pflotran_model,  &
                                            grid_clm_cell_ids_nindex, &
                                            grid_clm_npts_local, &
                                            map_id)
  ! 
  ! This routine maps CLM surface grid onto surface of PFLOTRAN 3D subsurface
  ! grid.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/09/13
  ! 

    use Input_Aux_module
    use Option_module
    use Realization_class
    use Grid_module
    use Patch_module
    use Coupler_module
    use Connection_module
    use String_module
    use clm_pflotran_interface_data
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename
    
    ! local
    PetscInt                           :: local_id, grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_clm_npts_ghost, source_mesh_id
    PetscInt                           :: dest_mesh_id
    PetscInt, pointer                  :: grid_pf_cell_ids_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_cell_ids_nindex_copy(:)
    PetscInt                           :: count
    PetscInt                           :: sum_connection
    PetscInt                           :: ghosted_id
    PetscInt                           :: iconn
    PetscInt                           :: istart
    PetscInt, pointer                  :: int_array(:)
    PetscBool                          :: found
    PetscScalar,pointer                :: v_loc(:)
    PetscErrorCode                     :: ierr

    Vec                                :: surf_ids
    Vec                                :: surf_ids_loc
    IS                                 :: is_from
    IS                                 :: is_to
    VecScatter                         :: vec_scat

    type(mapping_type), pointer        :: map
    type(option_type), pointer         :: option
    class(realization_type), pointer    :: realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch
    type(coupler_type), pointer        :: boundary_condition
    type(coupler_type), pointer        :: source_sink
    type(connection_set_type), pointer :: cur_connection_set

    option          => pflotran_model%option

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

    allocate(grid_clm_cell_ids_nindex_copy(grid_clm_npts_local))
    grid_clm_cell_ids_nindex_copy = grid_clm_cell_ids_nindex

    ! Choose the appriopriate map
    select case(map_id)
      case(CLM_SRF_TO_PF_2DSUB)
        map => pflotran_model%map_clm_srf_to_pf_2dsub
        source_mesh_id = CLM_SUB_MESH
        dest_mesh_id = PF_2DSUB_MESH
      case default
        option%io_buffer = 'Invalid map_id argument to ' // &
          'pflotranModelInitMappingSurf3D'
        call printErrMsg(option)
    end select

    ! Read mapping file
    if (index(map%filename, '.h5') > 0) then
      call MappingReadHDF5(map, map%filename, option)
    else
      call MappingReadTxtFile(map, map%filename, option)
    endif

    grid_clm_npts_ghost=0

    ! Allocate memory to identify if CLM cells are local or ghosted.
    ! Note: Presently all CLM cells are local
    allocate(grid_clm_local_nindex(grid_clm_npts_local))
    do local_id = 1, grid_clm_npts_local
      grid_clm_local_nindex(local_id) = 1 ! LOCAL
    enddo

    ! Mapping to/from surface of PFLOTRAN domain
    found=PETSC_FALSE
    grid_pf_npts_local = 0
    grid_pf_npts_ghost = 0

    select case (dest_mesh_id)

      case(PF_2DSUB_MESH)

        patch => realization%patch
        grid => patch%grid

        ! Destination mesh is PF_2DSUB_MESH
        boundary_condition => patch%boundary_conditions%first
        sum_connection = 0
        do 
          if (.not.associated(boundary_condition)) exit
          cur_connection_set => boundary_condition%connection_set

          if(StringCompare(boundary_condition%name,'clm_gflux_bc')) then

            found=PETSC_TRUE

            ! Allocate memory to save cell ids and flag for local cells
            allocate(grid_pf_cell_ids_nindex(cur_connection_set%num_connections))
            allocate(grid_pf_local_nindex(cur_connection_set%num_connections))
            grid_pf_npts_local = cur_connection_set%num_connections

            ! Save cell ids in application order 0-based
            do iconn=1,cur_connection_set%num_connections
              sum_connection = sum_connection + 1
              local_id = cur_connection_set%id_dn(iconn)
              ghosted_id = grid%nL2G(local_id)
              if (patch%imat(ghosted_id) <= 0) cycle
              grid_pf_cell_ids_nindex(iconn) = grid%nG2A(ghosted_id) - 1
              grid_pf_local_nindex(iconn) = 1
            enddo
          else
            sum_connection = sum_connection + cur_connection_set%num_connections
          endif
          boundary_condition => boundary_condition%next
        enddo

        ! Setting the number of cells constituting the surface of the 3D
        ! subsurface domain for each model.
        clm_pf_idata%nlclm_2dsub = grid_clm_npts_local
        clm_pf_idata%ngclm_2dsub = grid_clm_npts_local
        clm_pf_idata%nlpf_2dsub  = grid_pf_npts_local
        clm_pf_idata%ngpf_2dsub  = grid_pf_npts_local

      case default
        option%io_buffer='Unknown source mesh'
        call printErrMsg(option)
      
    end select

    if(.not.found) then
      pflotran_model%option%io_buffer = 'clm_gflux_bc not found in boundary conditions'
      call printErrMsg(pflotran_model%option)
    endif
    
    !
    ! Step-1: Find surface cells-ids of PFLOTRAN subsurface domain
    !
    allocate(v_loc(grid_pf_npts_local))
    v_loc = 1.d0
    call VecCreateSeq(PETSC_COMM_SELF, grid_pf_npts_local, surf_ids_loc, ierr)
    call VecCreateMPI(option%mycomm, grid%nlmax, PETSC_DECIDE, surf_ids, ierr)
    call VecSet(surf_ids, -1.d0, ierr)
    
    ! Set 1.0 to all cells that make up surface of PFLOTRAN subsurface domain
    call VecSetValues(surf_ids, grid_pf_npts_local, grid_pf_cell_ids_nindex, &
                      v_loc, INSERT_VALUES, ierr)
    deallocate(v_loc)
    call VecAssemblyBegin(surf_ids, ierr)
    call VecAssemblyEnd(surf_ids, ierr)

    call VecGetArrayF90(surf_ids, v_loc, ierr)
    count = 0
    do local_id=1,grid%nlmax
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr)

    count = 0
    do local_id=1,grid%nlmax
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif
    enddo
    call VecRestoreArrayF90(surf_ids, v_loc, ierr)

    !
    allocate(int_array(grid_pf_npts_local))
    do iconn = 1, grid_pf_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, grid_pf_npts_local, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm, grid_pf_npts_local, grid_pf_cell_ids_nindex, &
                         PETSC_COPY_VALUES, is_from, ierr)


    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    call ISDestroy(is_from, ierr)
    call ISDestroy(is_to, ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    count = 0
    do iconn = 1, grid_pf_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_pf_cell_ids_nindex(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    call VecDestroy(surf_ids_loc, ierr)

    !
    ! Step-2: Recompute 'map%s2d_iscr'
    !
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_nwts, surf_ids_loc, ierr)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)


    do iconn = 1, map%s2d_nwts
      int_array(iconn) = map%s2d_icsr(iconn)
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_from, ierr)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    call ISDestroy(is_from, ierr)
    call ISDestroy(is_to, ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        map%s2d_icsr(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    call VecDestroy(surf_ids_loc, ierr)
    
    if(count /= map%s2d_nwts) then
      option%io_buffer='No. of surface cells in mapping dataset does not ' // &
        'match surface cells on which BC is applied.'
      call printErrMsg(option)
    endif
    call VecDestroy(surf_ids, ierr)

    !
    ! Step-3: Find surface cells-ids of CLM subsurface domain
    !
    allocate(v_loc(grid_clm_npts_local))
    v_loc = 1.d0
    call VecCreateSeq(PETSC_COMM_SELF, grid_clm_npts_local, surf_ids_loc, ierr)
    call VecCreateMPI(option%mycomm, clm_pf_idata%nlclm_sub, PETSC_DECIDE, surf_ids, ierr)
    call VecSet(surf_ids, -1.d0, ierr)

    ! Set 1.0 to all cells that make up surface of CLM subsurface domain
    call VecSetValues(surf_ids, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                      v_loc, INSERT_VALUES, ierr)

    deallocate(v_loc)
    call VecAssemblyBegin(surf_ids, ierr)
    call VecAssemblyEnd(surf_ids, ierr)

    call VecGetArrayF90(surf_ids, v_loc, ierr)
    count = 0
    do local_id=1,clm_pf_idata%nlclm_sub
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr)

    count = 0
    do local_id=1,clm_pf_idata%nlclm_sub
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif
    enddo
    call VecRestoreArrayF90(surf_ids, v_loc, ierr)

    !
    allocate(int_array(grid_clm_npts_local))
    do iconn = 1, grid_clm_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, grid_clm_npts_local, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                         PETSC_COPY_VALUES, is_from, ierr)


    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    call ISDestroy(is_from, ierr)
    call ISDestroy(is_to, ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    count = 0
    do iconn = 1, grid_clm_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_clm_cell_ids_nindex_copy(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    call VecDestroy(surf_ids_loc, ierr)


    !
    ! Step-4: Recompute 'map%s2d_jscr'
    !
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_nwts, surf_ids_loc, ierr)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)


    do iconn = 1, map%s2d_nwts
      int_array(iconn) = map%s2d_jcsr(iconn)
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_from, ierr)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    call ISDestroy(is_from, ierr)
    call ISDestroy(is_to, ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        map%s2d_jcsr(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    call VecDestroy(surf_ids_loc, ierr)
    
    if(count /= map%s2d_nwts) then
      option%io_buffer='No. of surface cells in mapping dataset does not ' // &
        'match surface cells on which BC is applied.'
      call printErrMsg(option)
    endif
    call VecDestroy(surf_ids, ierr)

    select case(source_mesh_id)
      case(CLM_SUB_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_cell_ids_nindex_copy)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_SUB_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_pf_npts_local, &
                                        grid_pf_cell_ids_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_clm_npts_local, &
                                              grid_clm_npts_ghost, &
                                              grid_clm_cell_ids_nindex_copy, &
                                              grid_clm_local_nindex)
      case default
        option%io_buffer = 'Invalid argument source_mesh_id passed to ' // &
          'pflotranModelInitMappingSurf3D'
        call printErrMsg(option)
    end select

    deallocate(grid_pf_cell_ids_nindex)
    deallocate(grid_pf_local_nindex)

    call MappingDecompose(map, option)
    call MappingFindDistinctSourceMeshCellIds(map, option)
    call MappingCreateWeightMatrix(map, option)
    call MappingCreateScatterOfSourceMesh(map, option)

  end subroutine pflotranModelInitMapSrfTo2DSub

! ************************************************************************** !

  subroutine pflotranModelInitMapSrfToSrf(pflotran_model,  &
                                          grid_clm_cell_ids_nindex, &
                                          grid_clm_npts_local, &
                                          map_id)
  ! 
  ! This routine maps CLM surface grid onto PFLOTRAN 2D surface grid or
  ! vice-versa.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/09/13
  ! 

    use Input_Aux_module
    use Option_module
    use Realization_class
    use Grid_module
    use Patch_module
    use Coupler_module
    use Connection_module
    use String_module
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Surface_Realization_class, only : surface_realization_type
    use clm_pflotran_interface_data
    use Realization_class, only : realization_type
    use Surface_Realization_class, only : surface_realization_type
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename
    
    ! local
    PetscInt                           :: local_id, grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_clm_npts_ghost, source_mesh_id
    PetscInt                           :: dest_mesh_id
    PetscInt, pointer                  :: grid_pf_cell_ids_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_cell_ids_nindex_copy(:)
    PetscInt                           :: count
    PetscInt                           :: sum_connection
    PetscInt                           :: ghosted_id
    PetscInt                           :: iconn
    PetscInt                           :: istart
    PetscInt, pointer                  :: int_array(:)
    PetscBool                          :: found
    PetscScalar,pointer                :: v_loc(:)
    PetscErrorCode                     :: ierr

    Vec                                :: surf_ids
    Vec                                :: surf_ids_loc
    IS                                 :: is_from
    IS                                 :: is_to
    VecScatter                         :: vec_scat

    type(mapping_type), pointer        :: map
    type(option_type), pointer         :: option
    class(realization_type), pointer    :: realization
    class(surface_realization_type), pointer :: surf_realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch
    type(coupler_type), pointer        :: boundary_condition
    type(coupler_type), pointer        :: source_sink
    type(connection_set_type), pointer :: cur_connection_set

    option          => pflotran_model%option

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

#ifndef SURFACE_FLOW
    option%io_buffer='To support dest_mesh == PF_SRF_MESH, need to '// &
         'compiled with -DSURFACE_FLOW.'
    call printErrMsg(option)
#else
    allocate(grid_clm_cell_ids_nindex_copy(grid_clm_npts_local))
    grid_clm_cell_ids_nindex_copy = grid_clm_cell_ids_nindex

    ! Choose the appriopriate map
    select case(map_id)
      case(CLM_SRF_TO_PF_SRF)
        map => pflotran_model%map_clm_srf_to_pf_srf
        source_mesh_id = CLM_SRF_MESH
        dest_mesh_id = PF_SRF_MESH
      case(PF_SRF_TO_CLM_SRF)
        map => pflotran_model%map_pf_srf_to_clm_srf
        source_mesh_id = PF_SRF_MESH
        dest_mesh_id = CLM_SRF_MESH
      case default
        option%io_buffer = 'Invalid map_id argument to ' // &
          'pflotranModelInitMappingSurf2D'
        call printErrMsg(option)
    end select

    ! Read mapping file
    if (index(map%filename, '.h5') > 0) then
      call MappingReadHDF5(map, map%filename, option)
    else
      call MappingReadTxtFile(map, map%filename, option)
    endif

    grid_clm_npts_ghost=0

    ! Allocate memory to identify if CLM cells are local or ghosted.
    ! Note: Presently all CLM cells are local
    allocate(grid_clm_local_nindex(grid_clm_npts_local))
    do local_id = 1, grid_clm_npts_local
      grid_clm_local_nindex(local_id) = 1 ! LOCAL
    enddo

    ! Mapping to/from surface of PFLOTRAN domain
    ! Destination mesh is surface-mesh
    select type (simulation => pflotran_model%simulation)
      class is (subsurface_simulation_type)
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on surface simulations."
         call printErrMsg(pflotran_model%option)
      class is (surfsubsurface_simulation_type)
         surf_realization => simulation%surf_realization
      class is (surface_simulation_type)
         surf_realization => simulation%surf_realization
      class default
    end select
    patch => surf_realization%patch
    grid => patch%grid

    !
    ! Step-1: Find surface cells-ids of PFLOTRAN surface domain
    !
    grid_pf_npts_local = grid%nlmax
    grid_pf_npts_ghost = 0

    allocate(v_loc(grid%nlmax))
    allocate(grid_pf_cell_ids_nindex(grid%nlmax))
    allocate(grid_pf_local_nindex(grid%nlmax))
    
    grid_pf_local_nindex = 1
    call VecCreateMPI(option%mycomm, &
                      realization%patch%grid%nlmax, &
                      PETSC_DECIDE, &
                      surf_ids, &
                      ierr)
    call VecSet(surf_ids, -1.d0, ierr)

    do local_id = 1,grid%nlmax
      v_loc(local_id) = grid%unstructured_grid%cell_ids_natural(local_id)-1
      grid_pf_cell_ids_nindex(local_id) = &
        grid%unstructured_grid%nat_ids_of_other_grid(local_id)-1
    enddo

    !
    call VecSetValues(surf_ids, grid_pf_npts_local, grid_pf_cell_ids_nindex, &
                      v_loc, INSERT_VALUES, ierr)
    call VecAssemblyBegin(surf_ids, ierr)
    call VecAssemblyEnd(surf_ids, ierr)

    do local_id = 1,grid%nlmax
      grid_pf_cell_ids_nindex(local_id) = &
        grid%unstructured_grid%cell_ids_natural(local_id)-1
    enddo   
    
    !
    ! Step-2: Recompute 'map%s2d_icsr'
    !
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_nwts, surf_ids_loc, ierr)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)

    do iconn = 1, map%s2d_nwts
      if (source_mesh_id == PF_SRF_MESH) then
        int_array(iconn) = map%s2d_jcsr(iconn)
      else
        int_array(iconn) = map%s2d_icsr(iconn)
      endif
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_from, ierr)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    call ISDestroy(is_from, ierr)
    call ISDestroy(is_to, ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        if (source_mesh_id == PF_SRF_MESH) then
          map%s2d_jcsr(count) = INT(v_loc(iconn))
        else
          map%s2d_icsr(count) = INT(v_loc(iconn))
        endif
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    call VecDestroy(surf_ids_loc, ierr)
    
    if(count /= map%s2d_nwts) then
      write(*,*),'count = ',option%myrank,count,map%s2d_nwts
      option%io_buffer='No. of surface cells in mapping dataset does not ' // &
        'match surface cells on which BC is applied. [pflotranModelInitMappingSurf2D]'
      call printErrMsg(option)
    endif
    call VecDestroy(surf_ids, ierr)

    !
    ! Step-3: Find surface cells-ids of CLM subsurface domain
    !
    allocate(v_loc(grid_clm_npts_local))
    v_loc = 1.d0
    call VecCreateSeq(PETSC_COMM_SELF, grid_clm_npts_local, surf_ids_loc, ierr)
    call VecCreateMPI(option%mycomm, clm_pf_idata%nlclm_sub, PETSC_DECIDE, surf_ids, ierr)
    call VecSet(surf_ids, -1.d0, ierr)

    ! Set 1.0 to all cells that make up surface of CLM subsurface domain
    call VecSetValues(surf_ids, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                      v_loc, INSERT_VALUES, ierr)

    deallocate(v_loc)
    call VecAssemblyBegin(surf_ids, ierr)
    call VecAssemblyEnd(surf_ids, ierr)

    call VecGetArrayF90(surf_ids, v_loc, ierr)
    count = 0
    do local_id=1,clm_pf_idata%nlclm_sub
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr)

    count = 0
    do local_id=1,clm_pf_idata%nlclm_sub
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif
    enddo
    call VecRestoreArrayF90(surf_ids, v_loc, ierr)

    !
    allocate(int_array(grid_clm_npts_local))
    do iconn = 1, grid_clm_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, grid_clm_npts_local, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                         PETSC_COPY_VALUES, is_from, ierr)


    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    call ISDestroy(is_from, ierr)
    call ISDestroy(is_to, ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    count = 0
    do iconn = 1, grid_clm_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_clm_cell_ids_nindex_copy(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    call VecDestroy(surf_ids_loc, ierr)

    !
    ! Step-4: Recompute 'map%s2d_jscr'
    !
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_nwts, surf_ids_loc, ierr)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)


    do iconn = 1, map%s2d_nwts
      if (source_mesh_id == PF_SRF_MESH) then
        int_array(iconn) = map%s2d_icsr(iconn)
      else
        int_array(iconn) = map%s2d_jcsr(iconn)
      endif
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_from, ierr)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    call ISDestroy(is_from, ierr)
    call ISDestroy(is_to, ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vec_scat, ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        if (source_mesh_id == PF_SRF_MESH) then
          map%s2d_icsr(count) = INT(v_loc(iconn))
        else
          map%s2d_jcsr(count) = INT(v_loc(iconn))
        endif
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    
    if(count /= map%s2d_nwts) then
      write(*,*),'count = ',option%myrank,count,map%s2d_nwts
      option%io_buffer='No. of surface cells in mapping dataset does not ' // &
        'match surface cells on which BC is applied. [pflotranModelInitMappingSurf2D]'
      call printErrMsgByRank(option)
    endif
    call VecDestroy(surf_ids, ierr)

    select case(source_mesh_id)
      case(CLM_SRF_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_cell_ids_nindex_copy)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_SRF_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_pf_npts_local, &
                                        grid_pf_cell_ids_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_clm_npts_local, &
                                              grid_clm_npts_ghost, &
                                              grid_clm_cell_ids_nindex_copy, &
                                              grid_clm_local_nindex)
      case default
        option%io_buffer = 'Invalid argument source_mesh_id passed to ' // &
          'pflotranModelInitMappingSurf2D'
        call printErrMsg(option)
    end select

    deallocate(grid_pf_local_nindex)
    deallocate(grid_clm_cell_ids_nindex_copy)

    call MappingDecompose(map, option)
    call MappingFindDistinctSourceMeshCellIds(map, option)
    call MappingCreateWeightMatrix(map, option)
    call MappingCreateScatterOfSourceMesh(map, option)
#endif

  end subroutine pflotranModelInitMapSrfToSrf

! ************************************************************************** !

  subroutine pflotranModelStepperRunTillPauseTime(model, pause_time, dtime)
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
  ! 


    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer :: model
    PetscReal, intent(in) :: pause_time
    PetscReal, intent(in) :: dtime

    PetscReal :: pause_time1

    if (model%option%io_rank == model%option%myrank) then
       write(model%option%fid_out, *), '>>>> Inserting waypoint at pause_time = ', pause_time
    endif

    pause_time1 = pause_time + dtime!1800.0d0
    call pflotranModelInsertWaypoint(model, pause_time1)

    call model%simulation%RunToTime(pause_time)

    call pflotranModelDeleteWaypoint(model, pause_time)

    ! TODO(GB): Use XXXUpdateMassBalancePatch() to ensure mass balance
    ! betweent CLM calls

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
  ! 

    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type

    use Realization_class, only : realization_type
    use Surface_Realization_class, only : surface_realization_type
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
    waypoint%dt_max            = waypoint_time * UnitsConvertToInternal(word, model%option)!3153600.d0
    waypoint2 => WaypointCreate(waypoint)

    if (associated(realization)) then
       call WaypointInsertInList(waypoint, realization%waypoints)
    end if

    if (associated(surf_realization)) then
       call WaypointInsertInList(waypoint2, surf_realization%waypoints)
    end if

  end subroutine pflotranModelInsertWaypoint

! ************************************************************************** !

  subroutine pflotranModelDeleteWaypoint(model, waypoint_time)

    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type

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
       call WaypointDeleteFromList(waypoint, realization%waypoints)
    end if

    if (associated(surf_realization)) then
       call WaypointDeleteFromList(waypoint, surf_realization%waypoints)
    end if

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

    model%option%io_buffer = 'restart is not implemented in clm-pflotran.'

    if (.not. StringNull(restart_stamp)) then
       model%option%restart_flag = PETSC_TRUE
       model%option%restart_filename = &
            trim(model%option%global_prefix) // &
            trim(model%option%group_prefix) // &
            '-' // trim(restart_stamp) // '.chk'

    else
       call printWrnMsg(model%option)

    end if


  end subroutine pflotranModelSetupRestart

! ************************************************************************** !

  subroutine pflotranModelUpdateSourceSink(pflotran_model)
  ! 
  ! Update the source/sink term
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/22/2011
  ! 

    use clm_pflotran_interface_data
    use Connection_module
    use Coupler_module
    use Mapping_module
    use Option_module
    use Realization_class, only : realization_type
    use Simulation_Base_class, only : simulation_base_type
    use String_module
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type

    implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model

    class(realization_type), pointer          :: subsurf_realization
    type(coupler_type), pointer               :: source_sink
    type(connection_set_type), pointer        :: cur_connection_set
    PetscScalar, pointer                      :: qflx_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn
    PetscErrorCode                            :: ierr
    PetscInt                                  :: press_dof

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%qflx_clm, &
                                    clm_pf_idata%qflx_pf)

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
    call VecGetArrayF90(clm_pf_idata%qflx_pf,qflx_pf_loc,ierr)
    found = PETSC_FALSE
    source_sink => subsurf_realization%patch%source_sinks%first
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
          source_sink%flow_aux_real_var(press_dof,iconn) = qflx_pf_loc(iconn)
        enddo
      endif

      source_sink => source_sink%next
    enddo
    call VecRestoreArrayF90(clm_pf_idata%qflx_pf,qflx_pf_loc,ierr)

    if(.not.found) &
      call printErrMsg(pflotran_model%option,'clm_et_ss not found in ' // &
                       'source-sink list of subsurface model.')

  end subroutine pflotranModelUpdateSourceSink

! ************************************************************************** !

  subroutine pflotranModelUpdateFlowConds(pflotran_model)
  ! 
  ! This routine Updates boundary and source/sink condtions for PFLOTRAN that
  ! are prescribed by CLM
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  ! 

    use clm_pflotran_interface_data
    use Mapping_module
    use Option_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model

    call pflotranModelUpdateSourceSink(pflotran_model)
    if (pflotran_model%option%nsurfflowdof > 0) then
      call pflotranModelUpdateSurfSource(pflotran_model)
    endif

    if (pflotran_model%option%iflowmode == TH_MODE) then
      if (pflotran_model%option%nsurfflowdof == 0) then
        call pflotranModelUpdateSubsurfTCond(pflotran_model)
      else
        call pflotranModelUpdateSurfTCond(pflotran_model)
      endif
    endif

  end subroutine pflotranModelUpdateFlowConds

! ************************************************************************** !

  subroutine pflotranModelUpdateSurfSource(pflotran_model)
  ! 
  ! This routine updates surface source condition related to mass equation.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 11/11/2013
  ! 

    use clm_pflotran_interface_data
    use Connection_module
    use Coupler_module
    use Mapping_module
    use Option_module
    use Surface_Realization_class, only : surface_realization_type
    use Simulation_Base_class, only : simulation_base_type
    use String_module
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type

    implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model

    class(surface_realization_type), pointer  :: surf_realization
    type(coupler_type), pointer               :: source_sink
    type(connection_set_type), pointer        :: cur_connection_set
    PetscScalar, pointer                      :: rain_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn
    PetscErrorCode                            :: ierr
    PetscInt                                  :: press_dof

    call MappingSourceToDestination(pflotran_model%map_clm_srf_to_pf_srf, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rain_clm, &
                                    clm_pf_idata%rain_pf)

    ! Get pointer to surface-realization
    select type (simulation => pflotran_model%simulation)
      class is (surfsubsurface_simulation_type)
         surf_realization => simulation%surf_realization
      class default
         pflotran_model%option%io_buffer = " Unsupported simulation_type " // &
            " in pflotranModelUpdateSurfSource."
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
    call VecGetArrayF90(clm_pf_idata%rain_pf,rain_pf_loc,ierr)
    found = PETSC_FALSE
    source_sink => surf_realization%patch%source_sinks%first
    do
      if (.not.associated(source_sink)) exit

      cur_connection_set => source_sink%connection_set

      ! Find appropriate Source/Sink from the list of Source/Sinks
      if(StringCompare(source_sink%name,'clm_rain_srf_ss')) then

        found = PETSC_TRUE
        if (source_sink%flow_condition%rate%itype /= HET_VOL_RATE_SS) then
          call printErrMsg(pflotran_model%option,'clm_et_ss is not of ' // &
                           'HET_VOL_RATE_SS')
        endif

        do iconn = 1, cur_connection_set%num_connections
          source_sink%flow_aux_real_var(press_dof,iconn) = rain_pf_loc(iconn)
        enddo
      endif

      source_sink => source_sink%next
    enddo
    call VecRestoreArrayF90(clm_pf_idata%rain_pf,rain_pf_loc,ierr)

    if(.not.found) &
      call printErrMsg(pflotran_model%option,'clm_rain_srf_ss not found in ' // &
                       'source-sink list of surface model.')

  end subroutine pflotranModelUpdateSurfSource

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
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type

    implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model

    class(realization_type), pointer          :: subsurf_realization
    type(coupler_type), pointer               :: boundary_condition
    type(connection_set_type), pointer        :: cur_connection_set
    PetscScalar, pointer                      :: gflux_subsurf_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn
    PetscErrorCode                            :: ierr
    PetscInt                                  :: press_dof

    ! Map ground-heat flux from CLM--to--PF grid
    call MappingSourceToDestination(pflotran_model%map_clm_srf_to_pf_2dsub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%gflux_subsurf_clm, &
                                    clm_pf_idata%gflux_subsurf_pf)

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
    call VecGetArrayF90(clm_pf_idata%gflux_subsurf_pf,gflux_subsurf_pf_loc,ierr)
    found = PETSC_FALSE
    boundary_condition => subsurf_realization%patch%boundary_conditions%first
    do
      if (.not.associated(boundary_condition)) exit

      cur_connection_set => boundary_condition%connection_set

      ! Find appropriate BC from the list of boundary conditions
      if(StringCompare(boundary_condition%name,'clm_gflux_bc')) then

        if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
            /= NEUMANN_BC) then
          call printErrMsg(pflotran_model%option,'clm_gflux_bc is not of ' // &
                           'NEUMANN_BC')
        endif
        found = PETSC_TRUE

        do iconn = 1, cur_connection_set%num_connections
          boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
            gflux_subsurf_pf_loc(iconn)
        enddo
      endif

      boundary_condition => boundary_condition%next
    enddo
    call VecRestoreArrayF90(clm_pf_idata%gflux_subsurf_pf,gflux_subsurf_pf_loc,ierr)

    if(.not.found) &
      call printErrMsg(pflotran_model%option,'clm_gflux_bc not found in ' // &
                       'boundary-condition list of subsurface model.')

  end subroutine pflotranModelUpdateSubsurfTCond

! ************************************************************************** !

  subroutine pflotranModelUpdateSurfTCond(pflotran_model)
  ! 
  ! This routine updates surface source condition related to mass equation.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 11/11/2013
  ! 

    use clm_pflotran_interface_data
    use Connection_module
    use Coupler_module
    use Mapping_module
    use Option_module
    use Surface_Realization_class, only : surface_realization_type
    use Simulation_Base_class, only : simulation_base_type
    use String_module
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type

    implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model

    class(surface_realization_type), pointer  :: surf_realization
    type(coupler_type), pointer               :: source_sink
    type(connection_set_type), pointer        :: cur_connection_set
    PetscScalar, pointer                      :: gflux_surf_pf_loc(:)
    PetscScalar, pointer                      :: rain_temp_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn
    PetscErrorCode                            :: ierr
    PetscInt                                  :: press_dof
    PetscInt                                  :: temp_dof

    ! 1) Mapping energy flux
    call MappingSourceToDestination(pflotran_model%map_clm_srf_to_pf_srf, &
                                    pflotran_model%option, &
                                    clm_pf_idata%gflux_subsurf_clm, &
                                    clm_pf_idata%gflux_surf_pf)

    ! Get pointer to surface-realization
    select type (simulation => pflotran_model%simulation)
      class is (surfsubsurface_simulation_type)
         surf_realization => simulation%surf_realization
      class default
         pflotran_model%option%io_buffer = " Unsupported simulation_type " // &
            " in pflotranModelUpdateSurfSource."
         call printErrMsg(pflotran_model%option)
    end select

    ! Find value of pressure-dof depending on flow mode
    select case (pflotran_model%option%iflowmode)
      case (TH_MODE)
        press_dof = TH_PRESSURE_DOF
        temp_dof = TH_TEMPERATURE_DOF
      case default
        pflotran_model%option%io_buffer = 'Unsupported Flow mode'
        call printErrMsg(pflotran_model%option)
    end select

    ! Update the 'clm_et_ss' source/sink term
    call VecGetArrayF90(clm_pf_idata%gflux_surf_pf,gflux_surf_pf_loc,ierr)
    found = PETSC_FALSE
    source_sink => surf_realization%patch%source_sinks%first
    do
      if (.not.associated(source_sink)) exit

      cur_connection_set => source_sink%connection_set

      ! Find appropriate Source/Sink from the list of Source/Sinks
      if(StringCompare(source_sink%name,'clm_energy_srf_ss')) then

        found = PETSC_TRUE
        if (source_sink%flow_condition%energy_rate%itype /= HET_ENERGY_RATE_SS) then
          call printErrMsg(pflotran_model%option,'clm_et_ss is not of ' // &
                           'HET_ENERGY_RATE_SS')
        endif

        do iconn = 1, cur_connection_set%num_connections
          source_sink%flow_aux_real_var(press_dof,iconn) = 0.d0
          source_sink%flow_aux_real_var(temp_dof,iconn) = gflux_surf_pf_loc(iconn)
        enddo
      endif

      source_sink => source_sink%next
    enddo
    call VecRestoreArrayF90(clm_pf_idata%gflux_surf_pf,gflux_surf_pf_loc,ierr)

    if (.not.found) &
      call printErrMsg(pflotran_model%option,'clm_energy_srf_ss not found in ' // &
                       'source-sink list of surface model.')

    ! 2) Map temperature of rain water
    write(*,*),'call MappingSourceToDestination()'
    call MappingSourceToDestination(pflotran_model%map_clm_srf_to_pf_srf, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rain_temp_clm, &
                                    clm_pf_idata%rain_temp_pf)

    ! Update the 'clm_rain_srf_ss' source/sink term
    call VecGetArrayF90(clm_pf_idata%rain_temp_pf,rain_temp_pf_loc,ierr)
    found = PETSC_FALSE
    source_sink => surf_realization%patch%source_sinks%first
    do
      if (.not.associated(source_sink)) exit

      cur_connection_set => source_sink%connection_set

      ! Find appropriate Source/Sink from the list of Source/Sinks
      if(StringCompare(source_sink%name,'clm_rain_srf_ss')) then

        found = PETSC_TRUE
        if (source_sink%flow_condition%temperature%itype /= HET_DIRICHLET) then
          call printErrMsg(pflotran_model%option,'clm_rain_srf_ss is not of ' // &
                           'HET_DIRICHLET')
        endif

        do iconn = 1, cur_connection_set%num_connections
          source_sink%flow_aux_real_var(temp_dof,iconn) = rain_temp_pf_loc(iconn)
        enddo
      endif

      source_sink => source_sink%next
    enddo
    call VecRestoreArrayF90(clm_pf_idata%rain_temp_pf,rain_temp_pf_loc,ierr)

    if (.not.found) &
      call printErrMsg(pflotran_model%option,'clm_rain_srf_ss not found in ' // &
                       'source-sink list of surface model.')


  end subroutine pflotranModelUpdateSurfTCond

! ************************************************************************** !

  subroutine pflotranModelSurfaceSource(pflotran_model)
  ! 
  ! This routine updates source condtion for 'mass' equation of PFLOTRAN
  ! surface-flow model from CLM.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 9/18/2013
  ! 

    use clm_pflotran_interface_data
    use Coupler_module
    use Connection_module
    use Mapping_module
    use Option_module
    use Realization_class, only : realization_type
    use String_module
    use Simulation_Base_class, only : simulation_base_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Surface_Realization_class, only : surface_realization_type

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model

    type(coupler_type), pointer               :: source_sink
    type(connection_set_list_type), pointer   :: connection_set_list
    type(connection_set_type), pointer        :: cur_connection_set
    class(surface_realization_type), pointer  :: surf_realization
    PetscScalar, pointer                      :: rain_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn
    PetscErrorCode                            :: ierr

    call MappingSourceToDestination(pflotran_model%map_clm_srf_to_pf_srf, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rain_clm, &
                                    clm_pf_idata%rain_pf)

    select type (simulation => pflotran_model%simulation)
      class is (surfsubsurface_simulation_type)
         surf_realization => simulation%surf_realization
      class default
         nullify(surf_realization)
         pflotran_model%option%io_buffer = "pflotranModelSurfaceSource only " // &
              "works on combinations of surfsubsurface_simulation_type."
         call printErrMsg(pflotran_model%option)
    end select

    ! Source/sink terms -------------------------------------
    call VecGetArrayF90(clm_pf_idata%rain_pf,rain_pf_loc,ierr)
    found = PETSC_FALSE
    source_sink => surf_realization%patch%source_sinks%first
    do
      if (.not.associated(source_sink)) exit

      cur_connection_set => source_sink%connection_set

      ! Find appropriate Source/Sink from the list of Source/Sinks
      if(StringCompare(source_sink%name,'rain_from_clm_ss')) then

        found = PETSC_TRUE
        if (source_sink%flow_condition%rate%itype /= HET_VOL_RATE_SS) then
          call printErrMsg(pflotran_model%option,'rain_from_clm_ss is not of ' // &
                           'HET_VOL_RATE_SS')
        endif

        do iconn = 1, cur_connection_set%num_connections
          source_sink%flow_aux_real_var(ONE_INTEGER,iconn) = rain_pf_loc(iconn)
        enddo
      endif

      source_sink => source_sink%next
    enddo
    call VecRestoreArrayF90(clm_pf_idata%rain_pf,rain_pf_loc,ierr)

    if(.not.found) &
      call printErrMsg(pflotran_model%option,'rain_from_clm_ss not found in ' // &
                       'source-sink list of surface-flow model.')

  end subroutine pflotranModelSurfaceSource

! ************************************************************************** !

  subroutine pflotranModelGetUpdatedStates(pflotran_model)
  ! 
  ! This routine get updated states evoloved by PFLOTRAN.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 5/14/2013
  ! 

    use Option_module
    use Richards_module
    use Richards_Aux_module
    use TH_module
    use TH_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Surface_Realization_class, only : surface_realization_type
    use clm_pflotran_interface_data
    use Realization_class, only : realization_type

    type(pflotran_model_type), pointer  :: pflotran_model
    class(realization_type), pointer     :: realization

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

    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        call RichardsUpdateAuxVars(realization)
        call pflotranModelGetSaturation(pflotran_model)
      case (TH_MODE)
        call THUpdateAuxVars(realization)
        call pflotranModelGetSaturation(pflotran_model)
        call pflotranModelGetTemperature(pflotran_model)
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

  end subroutine pflotranModelGetUpdatedStates

! ************************************************************************** !

  subroutine pflotranModelGetSaturation(pflotran_model)
  ! 
  ! Extract soil saturation values simulated by
  ! PFLOTRAN in a PETSc vector.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/22/2011
  ! 

    use Option_module
    use Realization_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
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
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: sat_pf_p(:)
    PetscReal, pointer :: press_pf_p(:)

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
    global_aux_vars => patch%aux%Global%aux_vars

    ! Save the saturation values
    call VecGetArrayF90(clm_pf_idata%sat_pf, sat_pf_p, ierr)
    do local_id=1, grid%nlmax
      ghosted_id=grid%nL2G(local_id)
      sat_pf_p(local_id)=global_aux_vars(ghosted_id)%sat(1)
    enddo
    call VecRestoreArrayF90(clm_pf_idata%sat_pf, sat_pf_p, ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sat_pf, &
                                    clm_pf_idata%sat_clm)

    ! Save the pressure values
    call VecGetArrayF90(clm_pf_idata%press_pfp, press_pf_p, ierr)
    do local_id=1, grid%nlmax
      ghosted_id=grid%nL2G(local_id)
      press_pf_p(local_id)=global_aux_vars(ghosted_id)%pres(1)
    enddo
    call VecRestoreArrayF90(clm_pf_idata%press_pfp, press_pf_p, ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_pfp, &
                                    clm_pf_idata%press_clms)

  end subroutine pflotranModelGetSaturation

! ************************************************************************** !

  subroutine pflotranModelGetSurfaceFlowHead(pflotran_model)
  ! 
  ! This routine returns updated surface-flow standing head of water evoloved
  ! by PFLOTRAN.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 9/18/2013
  ! 

    use Option_module
    use Realization_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Surface_Realization_class, only : surface_realization_type
    use clm_pflotran_interface_data
    use Surface_Global_Aux_module
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model

    class(surface_realization_type), pointer  :: surf_realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(surface_global_auxvar_type), pointer   :: surf_global_aux_vars(:)
    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: h2osfc_pf_p(:)

    select type (simulation => pflotran_model%simulation)
      class is (surfsubsurface_simulation_type)
         surf_realization => simulation%surf_realization
      class default
         nullify(surf_realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on " // &
            "surfsubsurface_simulation_type simulations."
         call printErrMsg(pflotran_model%option)
    end select
    patch           => surf_realization%patch
    grid            => patch%grid
    surf_global_aux_vars => patch%surf_aux%SurfaceGlobal%aux_vars

    ! Save the standing head of water values
    call VecGetArrayF90(clm_pf_idata%h2osfc_pf, h2osfc_pf_p, ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ! Convert 'm' to 'mm'
      h2osfc_pf_p(local_id) = surf_global_aux_vars(ghosted_id)%head(1)*1000.d0
    enddo
    call VecRestoreArrayF90(clm_pf_idata%h2osfc_pf, h2osfc_pf_p, ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_srf_to_clm_srf, &
                                    pflotran_model%option, &
                                    clm_pf_idata%h2osfc_pf, &
                                    clm_pf_idata%h2osfc_clm)

  end subroutine pflotranModelGetSurfaceFlowHead

! ************************************************************************** !

  subroutine pflotranModelGetTemperature(pflotran_model)
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
    use TH_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Surface_Realization_class, only : surface_realization_type
    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_type), pointer           :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    type(th_auxvar_type), pointer             :: th_aux_vars(:)

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscScalar, pointer :: temp_pf_p(:)
    PetscReal, pointer :: sat_ice_pf_p(:)

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
    global_aux_vars => patch%aux%Global%aux_vars
    th_aux_vars     => patch%aux%TH%aux_vars

    call VecGetArrayF90(clm_pf_idata%temp_pf, temp_pf_p, ierr)
    do ghosted_id=1,grid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (local_id>0) then
        !call VecSetValues(clm_pf_idata%sat_ice_pf,1,local_id-1, &
        !                 global_aux_vars(ghosted_id)%sat_ice(1),INSERT_VALUES,ierr)
        temp_pf_p(local_id) = global_aux_vars(ghosted_id)%temp(1)
      endif
    enddo
    call VecRestoreArrayF90(clm_pf_idata%temp_pf, temp_pf_p, ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%temp_pf, &
                                    clm_pf_idata%temp_clm)

!    call VecAssemblyBegin(clm_pf_idata%sat_ice_pf,ierr)
!    call VecAssemblyEnd(clm_pf_idata%sat_ice_pf,ierr)
!    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
!                                    pflotran_model%option, &
!                                    clm_pf_idata%sat_ice_pf, &
!                                    clm_pf_idata%sat_ice_clm)

  end subroutine pflotranModelGetTemperature

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

  function pflotranModelNSurfCells3DDomain(pflotran_model)
  ! 
  ! This function returns the number of control volumes forming surface of
  ! the sub-surface domain. The subroutines assumes the following boundary
  ! condition is specified in inputdeck:
  ! - 'clm_gflux_bc': when running subsurface only simulation.
  ! - 'from_surface_bc': when running surface-subsurface simulation.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 6/03/2013
  ! 

    use Option_module
    use Coupler_module
    use String_module
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Realization_class

    implicit none

    type(pflotran_model_type), pointer :: pflotran_model

    class(realization_type), pointer :: realization
    type(coupler_list_type), pointer :: coupler_list
    type(coupler_type), pointer :: coupler
    type(simulation_base_type), pointer :: simulation
    character(len=MAXWORDLENGTH) :: condition_name
    PetscInt :: pflotranModelNSurfCells3DDomain
    PetscBool :: found

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

    ! Determine the BC coupler name to search from list of BCs depending on
    ! subsurface or surface-subsurface simulation.
    if (pflotran_model%option%nsurfflowdof == 0) then
      condition_name = 'clm_gflux_bc'
    else
      condition_name = 'from_surface_bc'
    endif

    coupler_list => realization%patch%boundary_conditions
    coupler => coupler_list%first
    found = PETSC_FALSE

    do
      if (.not.associated(coupler)) exit
      if (StringCompare(coupler%name,trim(condition_name))) then
        pflotranModelNSurfCells3DDomain=coupler%connection_set%num_connections
        found = PETSC_TRUE
      endif
      coupler => coupler%next
    enddo

    if(.not.found)  &
      call printErrMsg(pflotran_model%option, &
            'Missing from the input deck a BC named clm_gflux_bc')

  end function pflotranModelNSurfCells3DDomain

! ************************************************************************** !

  subroutine pflotranModelGetTopFaceArea(pflotran_model)
  ! 
  ! This subroutine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 6/10/2013
  ! 

    use Option_module
    use Patch_module
    use Discretization_module
    use Unstructured_Grid_Aux_module
    use Unstructured_Cell_module
    use Unstructured_Grid_module
    use Grid_module
    use clm_pflotran_interface_data
    use Utility_module, only : DotProduct, CrossProduct
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Surface_Realization_class, only : surface_realization_type
    use Realization_class, only : realization_type
    use Mapping_module

    implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer :: pflotran_model

    type(option_type), pointer :: option
    class(realization_type), pointer :: realization
    type(discretization_type), pointer :: discretization
    type(patch_type), pointer :: patch
    type(grid_type), pointer :: grid
    type(point_type) :: point1, point2, point3, point4

    PetscInt :: local_id
    PetscInt :: ghosted_id
    PetscInt :: iface
    PetscInt :: face_id
    PetscInt :: cell_type
    PetscInt :: vertex_ids(4)

    PetscReal :: area1

    PetscScalar, pointer :: area_p(:)
    PetscErrorCode :: ierr

    option => pflotran_model%option
    select type (simulation => pflotran_model%simulation)
      type is (subsurface_simulation_type)
         realization => simulation%realization
      type is (surfsubsurface_simulation_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call printErrMsg(pflotran_model%option)
    end select
    discretization => realization%discretization
    patch => realization%patch
    grid => patch%grid

    call VecGetArrayF90(clm_pf_idata%area_top_face_pf, area_p, ierr)
    if(grid%discretization_itype == STRUCTURED_GRID) then
      ! Structured grid
      do ghosted_id=1,grid%ngmax
        local_id = grid%nG2L(ghosted_id)
        if(local_id>0) then
          area1 = grid%structured_grid%dx(ghosted_id)* &
                  grid%structured_grid%dy(ghosted_id)
          area_p(local_id) = area1
        endif
      enddo
    else if (grid%discretization_itype == UNSTRUCTURED_GRID) then
      ! Unstructured grid
      do local_id = 1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        cell_type = grid%unstructured_grid%cell_type(local_id)

        ! Find iface
        if (cell_type == HEX_TYPE) then
          iface = 6
        else if (cell_type == WEDGE_TYPE) then
          iface = 5
        else
          call printErrMsg(pflotran_model%option, &
            'Only hex and wedge cell_type supported in CLM-PFLOTRAN')
        endif

        ! Get face-id
        face_id = grid%unstructured_grid%cell_to_face_ghosted(iface, ghosted_id)

        ! Save face area
        area_p(local_id) = grid%unstructured_grid%face_area(face_id)
      enddo
    endif
    call VecRestoreArrayF90(clm_pf_idata%area_top_face_pf, area_p, ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%area_top_face_pf, &
                                    clm_pf_idata%area_top_face_clm)

  end subroutine pflotranModelGetTopFaceArea

! ************************************************************************** !

  subroutine pflotranModelDestroy(model)
  ! 
  ! Deallocates the pflotranModel object
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! 

    use PFLOTRAN_Factory_module, only : PFLOTRANFinalize
    use Option_module, only : OptionFinalize

    implicit none

    type(pflotran_model_type), pointer :: model

    ! FIXME(bja, 2013-07) none of the mapping information appears to
    ! be cleaned up, so we are leaking memory....

    call model%simulation%FinalizeRun()
    call model%simulation%Strip()
    deallocate(model%simulation)
    nullify(model%simulation)

    if (model%b_out_bgc_rate) then
    close(100)
    endif

    call PFLOTRANFinalize(model%option)
    call OptionFinalize(model%option)

    deallocate(model)

  end subroutine pflotranModelDestroy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
! THE FOLLOWING BLOCKS OF CODES ARE FOR CLM-PFLOTRAN BGC COUPLING
!
!
  ! ************************************************************************** !
  !
  ! pflotranModelSetInitialConcentrations:
  !  Get CLM concentrations (C, N), converts from CLM units into PFLOTRAN units.
  !  and set concentrations in PFLOTRAN
  ! author: Guoping Tang
  ! date: 08/16/2013
  ! ************************************************************************** !
  subroutine pflotranModelSetInitialConcentrations(pflotran_model)

    use Realization_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_Aux_module
    use Immobile_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module
    use Discretization_module

    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type

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
!    type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
!    type(reactive_transport_auxvar_type), pointer :: rt_aux_var

    PetscErrorCode     :: ierr
    PetscInt           :: local_id
    PetscReal, pointer :: xx_p(:)

    PetscScalar, pointer :: decomp_cpools_vr_lit1_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_lit2_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_lit3_pf_loc(:) ! (gC/m3)
    !PetscScalar, pointer :: decomp_cpools_vr_cwd_pf_loc(:)  ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som1_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som2_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som3_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som4_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_npools_vr_lit1_pf_loc(:) ! (gN/m3)
    PetscScalar, pointer :: decomp_npools_vr_lit2_pf_loc(:) ! (gN/m3)
    PetscScalar, pointer :: decomp_npools_vr_lit3_pf_loc(:) ! (gN/m3)
    !PetscScalar, pointer :: decomp_npools_vr_cwd_pf_loc(:)  ! (gN/m3)
    PetscScalar, pointer :: sminn_vr_pf_loc(:)              ! (gN/m3)
    !PetscScalar, pointer :: smin_no3_vr_pf_loc(:)           ! (gN/m3)
    !PetscScalar, pointer :: smin_nh4_vr_pf_loc(:)           ! (gN/m3)

    PetscInt :: ispec_c, ispec_n, ispec_no3, ispec_nh4, offset, offsetim
    PetscInt :: ispec_lit1c, ispec_lit2c, ispec_lit3c
    PetscInt :: ispec_lit1n, ispec_lit2n, ispec_lit3n
    PetscInt :: ispec_som1, ispec_som2, ispec_som3, ispec_som4
    PetscInt :: ispec_plantn

    character(len=MAXWORDLENGTH) :: word
    PetscReal, parameter :: C_molecular_weight = 12.0107d0
    PetscReal, parameter :: N_molecular_weight = 14.0067d0

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

    word = "LabileC"
    ispec_lit1c = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "CelluloseC"
    ispec_lit2c = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "LigninC"
    ispec_lit3c = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "LabileN"
    ispec_lit1n = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "CelluloseN"
    ispec_lit2n = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "LigninN"
    ispec_lit3n = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "SOM1"
    ispec_som1  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "SOM2"
    ispec_som2  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "SOM3"
    ispec_som3  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)


    word = "SOM4"
    ispec_som4  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "N"
    ispec_n  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "C"
    ispec_c  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

!    word = "NO3-"
!    ispec_no3  = GetPrimarySpeciesIDFromName(word, &
!                  realization%reaction,realization%option)

!    word = "NH4+"
!    ispec_nh4  = GetPrimarySpeciesIDFromName(word, &
!                  realization%reaction,realization%option)

    word = "PlantN"
    ispec_plantn  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_lit1_clmp, &
                                    clm_pf_idata%decomp_cpools_vr_lit1_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_lit2_clmp, &
                                    clm_pf_idata%decomp_cpools_vr_lit2_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_lit3_clmp, &
                                    clm_pf_idata%decomp_cpools_vr_lit3_pfs)

    !call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%decomp_cpools_vr_cwd_clmp, &
    !                                clm_pf_idata%decomp_cpools_vr_cwd_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som1_clmp, &
                                    clm_pf_idata%decomp_cpools_vr_som1_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som2_clmp, &
                                    clm_pf_idata%decomp_cpools_vr_som2_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som3_clmp, &
                                    clm_pf_idata%decomp_cpools_vr_som3_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som4_clmp, &
                                    clm_pf_idata%decomp_cpools_vr_som4_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_lit1_clmp, &
                                    clm_pf_idata%decomp_npools_vr_lit1_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_lit2_clmp, &
                                    clm_pf_idata%decomp_npools_vr_lit2_pfs)
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_lit3_clmp, &
                                    clm_pf_idata%decomp_npools_vr_lit3_pfs)

    !call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%decomp_npools_vr_cwd_clmp, &
    !                                clm_pf_idata%decomp_npools_vr_cwd_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sminn_vr_clmp, &
                                    clm_pf_idata%sminn_vr_pfs)

    !call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%smin_no3_vr_clmp, &
    !                                clm_pf_idata%smin_no3_vr_pfs)

    !call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%smin_nh4_vr_clmp, &
    !                                clm_pf_idata%smin_nh4_vr_pfs)

    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit1_pfs, &
                        decomp_cpools_vr_lit1_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit2_pfs, &
                        decomp_cpools_vr_lit2_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit3_pfs, &
                        decomp_cpools_vr_lit3_pf_loc, ierr)
    !call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_cwd_pfs,  &
    !                    decomp_cpools_vr_cwd_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som1_pfs, &
                        decomp_cpools_vr_som1_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som2_pfs, &
                        decomp_cpools_vr_som2_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som3_pfs, &
                        decomp_cpools_vr_som3_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som4_pfs, &
                        decomp_cpools_vr_som4_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit1_pfs, &
                        decomp_npools_vr_lit1_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit2_pfs, &
                        decomp_npools_vr_lit2_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit3_pfs, &
                        decomp_npools_vr_lit3_pf_loc, ierr)
    !call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_cwd_pfs, &
    !                    decomp_npools_vr_cwd_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%sminn_vr_pfs, sminn_vr_pf_loc, ierr)
    !call VecGetArrayF90(clm_pf_idata%smin_no3_vr_pfs, smin_no3_vr_pf_loc, ierr)
    !call VecGetArrayF90(clm_pf_idata%smin_nh4_vr_pfs, smin_nh4_vr_pf_loc, ierr)

    call VecGetArrayF90(field%tran_xx,xx_p,ierr)

    do local_id = 1, grid%nlmax
      if (grid%nG2L(local_id) < 0) cycle ! bypass ghosted corner cells
      !geh - Ignore inactive cells with inactive materials

      if (patch%imat(local_id) <= 0) cycle

      offset = (local_id - 1)*realization%reaction%ncomp

!      xx_p(offset + ispec_no3) = smin_no3_vr_pf_loc(local_id)
!      xx_p(offset + ispec_nh4) = smin_nh4_vr_pf_loc(local_id)

      offsetim = offset + realization%reaction%offset_immobile

      xx_p(offsetim + ispec_lit1c) = max(decomp_cpools_vr_lit1_pf_loc(local_id) &
                                        / C_molecular_weight, 1.0d-20)
      xx_p(offsetim + ispec_lit2c) = max(decomp_cpools_vr_lit2_pf_loc(local_id) &
                                        / C_molecular_weight, 1.0d-20)
      xx_p(offsetim + ispec_lit3c) = max(decomp_cpools_vr_lit3_pf_loc(local_id) &
                                        / C_molecular_weight, 1.0d-20)
      xx_p(offsetim + ispec_lit1n) = max(decomp_npools_vr_lit1_pf_loc(local_id) &
                                        / N_molecular_weight, 1.0d-20)
      xx_p(offsetim + ispec_lit2n) = max(decomp_npools_vr_lit2_pf_loc(local_id) &
                                        / N_molecular_weight, 1.0d-20)
      xx_p(offsetim + ispec_lit3n) = max(decomp_npools_vr_lit3_pf_loc(local_id) &
                                        / N_molecular_weight, 1.0d-20)
      xx_p(offsetim + ispec_som1) = max(decomp_cpools_vr_som1_pf_loc(local_id) &
                                        / C_molecular_weight, 1.0d-20)
      xx_p(offsetim + ispec_som2) = max(decomp_cpools_vr_som2_pf_loc(local_id) &
                                        / C_molecular_weight, 1.0d-20)
      xx_p(offsetim + ispec_som3) = max(decomp_cpools_vr_som3_pf_loc(local_id) &
                                        / C_molecular_weight, 1.0d-20)
      xx_p(offsetim + ispec_som4) = max(decomp_cpools_vr_som4_pf_loc(local_id) &
                                        / C_molecular_weight, 1.0d-20)
      xx_p(offsetim + ispec_n) = max(sminn_vr_pf_loc(local_id) &
                                        / N_molecular_weight, 1.0d-20)

    enddo

    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit1_pfs, decomp_cpools_vr_lit1_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit2_pfs, decomp_cpools_vr_lit2_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit3_pfs, decomp_cpools_vr_lit3_pf_loc, ierr)
    !call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_cwd_pfs,  decomp_cpools_vr_cwd_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som1_pfs, decomp_cpools_vr_som1_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som2_pfs, decomp_cpools_vr_som2_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som3_pfs, decomp_cpools_vr_som3_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som4_pfs, decomp_cpools_vr_som4_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit1_pfs, decomp_npools_vr_lit1_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit2_pfs, decomp_npools_vr_lit2_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit3_pfs, decomp_npools_vr_lit3_pf_loc, ierr)
    !call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_cwd_pfs,  decomp_npools_vr_cwd_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%sminn_vr_pfs, sminn_vr_pf_loc, ierr)
    !call VecRestoreArrayF90(clm_pf_idata%smin_no3_vr_pfs, smin_no3_vr_pf_loc, ierr)
    !call VecRestoreArrayF90(clm_pf_idata%smin_nh4_vr_pfs, smin_nh4_vr_pf_loc, ierr)

    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr)

    call DiscretizationGlobalToLocal(realization%discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

    call VecCopy(field%tran_xx,field%tran_yy,ierr)
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

  end subroutine pflotranModelSetInitialConcentrations

  ! ************************************************************************** !
  !> This routine Updates TH drivers for PFLOTRAN bgc that are from CLM
  !! for testing PFLOTRAN-BGC mode
  !!
  !> @author
  !! F.-M. Yuan
  !!
  !! date: 9/3/2013
  ! ************************************************************************** !
  subroutine pflotranModelUpdateTHfromCLM(pflotran_model, pf_hmode, pf_tmode)

    use Option_module
    use Realization_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
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
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: soillsat_pf_loc(:), soilisat_pf_loc(:)
    PetscReal, pointer :: soilt_pf_loc(:)

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
    global_aux_vars => patch%aux%Global%aux_vars

    ! Save the liq saturation values from CLM to PFLOTRAN, if needed
    if (.not.pf_hmode) then
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%soillsat_clmp, &
                                    clm_pf_idata%soillsat_pfs)
        call VecGetArrayF90(clm_pf_idata%soillsat_pfs, soillsat_pf_loc, ierr)
        do local_id=1, grid%nlmax
            ghosted_id=grid%nL2G(local_id)
            global_aux_vars(ghosted_id)%sat(1)=soillsat_pf_loc(local_id)
        enddo
        call VecRestoreArrayF90(clm_pf_idata%soillsat_pfs, soillsat_pf_loc, ierr)
    endif

    ! Save soil temperature values from CLM to PFLOTRAN, if needed
    if (.not.pf_tmode) then
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%soilt_clmp, &
                                    clm_pf_idata%soilt_pfs)
        call VecGetArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
        do local_id=1, grid%nlmax
            ghosted_id=grid%nL2G(local_id)
            global_aux_vars(ghosted_id)%temp(1)=soilt_pf_loc(local_id)
        enddo
        call VecRestoreArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
    endif

  end subroutine pflotranModelUpdateTHfromCLM

  ! ************************************************************************** !
  !> This routine Updates O2 for PFLOTRAN bgc that are from CLM
  !! for testing PFLOTRAN-BGC mode
  !!
  !> @author
  !! Guoping Tang
  !!
  !! date: 1/6/2014
  ! ************************************************************************** !
  subroutine pflotranModelUpdateO2fromCLM(pflotran_model)

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%o2_decomp_depth_unsat_clm, &
                                    clm_pf_idata%o2_decomp_depth_unsat_pf)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%conc_o2_unsat_clm, &
                                    clm_pf_idata%conc_o2_unsat_pf)

!    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
!                                    pflotran_model%option, &
!                                    clm_pf_idata%o2_decomp_depth_sat_clmp, &
!                                    clm_pf_idata%o2_decomp_depth_sat_pfs)
!   because of NaN t6g
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%conc_o2_sat_clm, &
                                    clm_pf_idata%conc_o2_sat_pf)

  end subroutine pflotranModelUpdateO2fromCLM

  ! ************************************************************************** !
  !
  ! pflotranModelSetInitialTHStatesfromCLM: Set initial TH States from CLM
  !
  ! Note: This subroutine directly set initial soil temperature and saturation from CLM
  !       It's needed because of uniform initialization of TH states in PFLOTRAN, which
  !       are from the input card.
  ! (This is different from the 'pflotranModelUpdateTHfromCLM', which pass TH from CLM to
  !   pflotran's global variables and will not affect the internal vec of TH mode).

  ! author: Fengming YUAN
  ! date: 9/23/2013
  ! ************************************************************************** !
subroutine pflotranModelSetInitialTHStatesfromCLM(pflotran_model)

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
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
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
    PetscInt           :: local_id, ghosted_id, istart, iend
    PetscReal, pointer :: xx_loc_p(:)

    PetscScalar, pointer :: soilt_pf_loc(:)    ! temperature [oC]
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
                                    clm_pf_idata%press_clm, &
                                    clm_pf_idata%press_pf)
      case (TH_MODE)
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%soilt_clmp, &
                                    clm_pf_idata%soilt_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_clm, &
                                    clm_pf_idata%press_pf)
      case default
        pflotran_model%option%io_buffer='pflotranModelSetInitialTHStatesfromCLM ' // &
          'not implmented for this mode.'
        call printErrMsg(pflotran_model%option)
    end select

    call VecGetArrayF90(field%flow_xx, xx_loc_p, ierr)
    call VecGetArrayF90(clm_pf_idata%press_pf, soilpress_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)

    do local_id = 1, grid%ngmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif

       iend = ghosted_id*pflotran_model%option%nflowdof
       istart = iend-pflotran_model%option%nflowdof+1

       xx_loc_p(istart)  = soilpress_pf_loc(local_id)
       if (pflotran_model%option%iflowmode .eq. TH_MODE)  then
            xx_loc_p(istart+1)= soilt_pf_loc(local_id)
       end if
    enddo

    call VecRestoreArrayF90(field%flow_xx, xx_loc_p, ierr)
    call VecRestoreArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_pf, soilpress_pf_loc, ierr)

    call DiscretizationGlobalToLocal(realization%discretization, field%flow_xx, &
         field%flow_xx_loc, NFLOWDOF)
    call VecCopy(field%flow_xx, field%flow_yy, ierr)

    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        call RichardsUpdateAuxVars(realization)
      case (TH_MODE)
        call THUpdateAuxVars(realization)
      case default
        pflotran_model%option%io_buffer='pflotranModelSetInitialTHStatesfromCLM ' // &
          'not implmented for this mode.'
        call printErrMsg(pflotran_model%option)
    end select

end subroutine pflotranModelSetInitialTHStatesfromCLM

  ! ************************************************************************** !
  ! pflotranModelSetSoilHbcs()
  ! refresh Hydrological BC variables from CLM to PF
  !
  ! by 1-18-2013: only water pressure-head type (dirichlet) available
  ! ************************************************************************** !
subroutine pflotranModelSetSoilHbcs(pflotran_model)

    use Realization_class
    use Option_module
    use Patch_module
    use Grid_module
    use Field_module
    use Coupler_module
    use Connection_module

    use TH_Aux_module
    use TH_module
    use Richards_module
    use Richards_Aux_module

    use String_module

    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type

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

    type(coupler_type), pointer :: boundary_condition
    type(connection_set_type), pointer :: cur_connection_set
    PetscInt :: ghosted_id, local_id, sum_connection, press_dof, iconn

    PetscErrorCode     :: ierr

    PetscScalar, pointer :: press_subsurf_pf_loc(:)     ! subsurface top boundary pressure-head (Pa) (dirichlet BC)
    PetscScalar, pointer :: qflux_subsurf_pf_loc(:)     ! subsurface top boundary infiltration rate (m/s) (neumann BC)
    !PetscScalar, pointer :: press_subbase_pf_loc(:)    ! bottom boundary pressure-head (Pa)

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
    field           => realization%field

    call MappingSourceToDestination(pflotran_model%map_clm_srf_to_pf_2dsub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_subsurf_clmp, &
                                    clm_pf_idata%press_subsurf_pfs)
    call MappingSourceToDestination(pflotran_model%map_clm_srf_to_pf_2dsub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%qflux_subsurf_clmp, &
                                    clm_pf_idata%qflux_subsurf_pfs)

    !call MappingSourceToDestination(pflotran_model%map_clm_srf_to_pf_2dsub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%press_subbase_clmp, &
    !                                clm_pf_idata%press_subbase_pfs)

    ! interface vecs of PF
    call VecGetArrayF90(clm_pf_idata%press_subsurf_pfs,  press_subsurf_pf_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%qflux_subsurf_pfs,  qflux_subsurf_pf_loc,  ierr)
    !call VecGetArrayF90(clm_pf_idata%press_subbase_pfs,  press_subbase_pf_loc,  ierr)

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

    boundary_condition => patch%boundary_conditions%first
    do
       if (.not.associated(boundary_condition)) exit

       cur_connection_set => boundary_condition%connection_set
       if(StringCompare(boundary_condition%name,'clm_gflux_bc')) then

           do iconn = 1, cur_connection_set%num_connections
              local_id = cur_connection_set%id_dn(iconn)
              ghosted_id = grid%nL2G(local_id)
              if (patch%imat(ghosted_id) <= 0) cycle

              if (boundary_condition%flow_condition%itype(press_dof) == DIRICHLET_BC) then
                   boundary_condition%flow_aux_real_var(press_dof,iconn)= &
                       press_subsurf_pf_loc(iconn)
              else if (boundary_condition%flow_condition%itype(press_dof) == NEUMANN_BC) then
                   boundary_condition%flow_aux_real_var(press_dof,iconn)= &
                       qflux_subsurf_pf_loc(iconn)
              end if

           enddo

       endif

       boundary_condition => boundary_condition%next

    enddo

    call VecRestoreArrayF90(clm_pf_idata%press_subsurf_pfs, press_subsurf_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%qflux_subsurf_pfs, qflux_subsurf_pf_loc, ierr)
    !call VecRestoreArrayF90(clm_pf_idata%press_subbase_pfs, press_subbase_pf_loc, ierr)

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

end subroutine pflotranModelSetSoilHbcs


! ************************************************************************** !

subroutine pflotranModelGetSoilProp(pflotran_model)
  !
  ! Pass physical properties from PFLOTRAN to CLM
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
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type

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

    ! pf interval variables
    PetscReal, pointer :: porosity_loc_p(:)

    ! clm-pf-interface Vecs for MVM parameters
    PetscScalar, pointer :: porosity_loc_pfp(:)  ! soil porosity
    PetscScalar, pointer :: sr_loc_pfp(:)        ! residual soil vwc
    PetscScalar, pointer :: alpha_loc_pfp(:)     ! alfa
    PetscScalar, pointer :: lamda_loc_pfp(:)     ! lamda
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

    call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)

    call VecGetArrayF90(clm_pf_idata%porosity_pfp, porosity_loc_pfp, ierr)
    call VecGetArrayF90(clm_pf_idata%sr_pfp, sr_loc_pfp, ierr)
    call VecGetArrayF90(clm_pf_idata%alpha_pfp, alpha_loc_pfp, ierr)
    call VecGetArrayF90(clm_pf_idata%lamda_pfp, lamda_loc_pfp, ierr)
    call VecGetArrayF90(clm_pf_idata%pcwmax_pfp, pcwmax_loc_pfp, ierr)

    do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (associated(patch%imat)) then
           if (patch%imat(ghosted_id) <= 0) cycle
        endif

        porosity_loc_pfp(local_id) = porosity_loc_p(ghosted_id)

        saturation_function => patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr

        sr_loc_pfp(local_id) = saturation_function%Sr(pflotran_model%option%nphase)
        alpha_loc_pfp(local_id) = saturation_function%alpha
        lamda_loc_pfp(local_id) = saturation_function%lambda
        pcwmax_loc_pfp(local_id) = saturation_function%pcwmax

    enddo

    call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
    call VecRestoreArrayF90(clm_pf_idata%porosity_pfp, porosity_loc_pfp, ierr)
    call VecRestoreArrayF90(clm_pf_idata%sr_pfp, sr_loc_pfp, ierr)
    call VecRestoreArrayF90(clm_pf_idata%alpha_pfp, alpha_loc_pfp, ierr)
    call VecRestoreArrayF90(clm_pf_idata%lamda_pfp, lamda_loc_pfp, ierr)
    call VecRestoreArrayF90(clm_pf_idata%pcwmax_pfp, pcwmax_loc_pfp, ierr)

    !
    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%porosity_pfp, &
                                    clm_pf_idata%porosity_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sr_pfp, &
                                    clm_pf_idata%sr_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%alpha_pfp, &
                                    clm_pf_idata%alpha_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%lamda_pfp, &
                                    clm_pf_idata%lamda_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%pcwmax_pfp, &
                                    clm_pf_idata%pcwmax_clms)

  end subroutine pflotranModelGetSoilProp


  ! ************************************************************************** !
  ! pflotranModelSetBGCRates:
  !  Get CLM litter, som, mineral N production and plant demand rates
  !  Convert from CLM units into PFLOTRAN units.
  !  Set values in PFLOTRAN
  ! author: Guoping Tang
  ! date: 09/03/2013
  ! ************************************************************************** !
  subroutine pflotranModelSetBGCRates(pflotran_model)

    use Realization_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_Aux_module
    use Immobile_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module
    use Discretization_module

    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
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
    PetscInt           :: local_id

    PetscInt:: ghosted_id, offset, offsetim

    PetscScalar, pointer :: rate_pf_loc(:)   !
    PetscScalar, pointer :: rate_lit1c_pf_loc(:)   !
    PetscScalar, pointer :: rate_lit2c_pf_loc(:)   !
    PetscScalar, pointer :: rate_lit3c_pf_loc(:)   !
    PetscScalar, pointer :: rate_lit1n_pf_loc(:)   !
    PetscScalar, pointer :: rate_lit2n_pf_loc(:)   !
    PetscScalar, pointer :: rate_lit3n_pf_loc(:)   !
    PetscScalar, pointer :: rate_minn_pf_loc(:)   !
    PetscScalar, pointer :: rate_plantnuptake_pf_loc(:)   !
    PetscScalar, pointer :: rate_nleached_pf_loc(:)   !
    PetscScalar, pointer :: rate_ndenitri_pf_loc(:)   !

    PetscReal, pointer :: volume_p(:)

    character(len=MAXWORDLENGTH) :: word
    PetscReal, parameter :: C_molecular_weight = 12.0107d0
    PetscReal, parameter :: N_molecular_weight = 14.0067d0

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

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_lit1c_clmp, &
                                    clm_pf_idata%rate_lit1c_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_lit2c_clmp, &
                                    clm_pf_idata%rate_lit2c_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_lit3c_clmp, &
                                    clm_pf_idata%rate_lit3c_pfs)

    !call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%rate_cwdc_clmp, &
    !                                clm_pf_idata%rate_cwdc_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_lit1n_clmp, &
                                    clm_pf_idata%rate_lit1n_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_lit2n_clmp, &
                                    clm_pf_idata%rate_lit2n_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_lit3n_clmp, &
                                    clm_pf_idata%rate_lit3n_pfs)

    !call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%rate_cwdn_clmp, &
    !                                clm_pf_idata%rate_cwdn_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_minn_clmp, &
                                    clm_pf_idata%rate_minn_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_plantnuptake_clmp, &
                                    clm_pf_idata%rate_plantnuptake_pf)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_nleached_clmp, &
                                    clm_pf_idata%rate_nleached_pf)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%rate_ndenitri_clmp, &
                                    clm_pf_idata%rate_ndenitri_pf)

!   get cell volume as mass transfer rate unit is mol/s
    call VecGetArrayReadF90(field%volume,volume_p,ierr)

!   output rate for debug
    if (pflotran_model%b_out_bgc_rate) then
        call VecGetArrayReadF90(clm_pf_idata%rate_lit1c_pfs, rate_lit1c_pf_loc, ierr)
        call VecGetArrayReadF90(clm_pf_idata%rate_lit2c_pfs, rate_lit2c_pf_loc, ierr)
        call VecGetArrayReadF90(clm_pf_idata%rate_lit3c_pfs, rate_lit3c_pf_loc, ierr)
        call VecGetArrayReadF90(clm_pf_idata%rate_lit1n_pfs, rate_lit1n_pf_loc, ierr)
        call VecGetArrayReadF90(clm_pf_idata%rate_lit2n_pfs, rate_lit2n_pf_loc, ierr)
        call VecGetArrayReadF90(clm_pf_idata%rate_lit3n_pfs, rate_lit3n_pf_loc, ierr)
        call VecGetArrayReadF90(clm_pf_idata%rate_minn_pfs,  rate_minn_pf_loc, ierr)
        call VecGetArrayReadF90(clm_pf_idata%rate_plantnuptake_pf, rate_plantnuptake_pf_loc, ierr)
        call VecGetArrayReadF90(clm_pf_idata%rate_nleached_pf,  rate_nleached_pf_loc, ierr)
        call VecGetArrayReadF90(clm_pf_idata%rate_ndenitri_pf,  rate_ndenitri_pf_loc, ierr)

        do local_id = 1, grid%nlmax
            write(100, *) rate_lit1c_pf_loc(local_id), &
                     rate_lit2c_pf_loc(local_id), &
                     rate_lit3c_pf_loc(local_id), &
                     rate_lit1n_pf_loc(local_id), &
                     rate_lit2n_pf_loc(local_id), &
                     rate_lit3n_pf_loc(local_id), &
                     rate_minn_pf_loc(local_id), &
                     rate_plantnuptake_pf_loc(local_id), &
                     rate_nleached_pf_loc(local_id), &
                     rate_ndenitri_pf_loc(local_id)
        enddo

        call VecRestoreArrayReadF90(clm_pf_idata%rate_lit1c_pfs, rate_lit1c_pf_loc, ierr)
        call VecRestoreArrayReadF90(clm_pf_idata%rate_lit2c_pfs, rate_lit2c_pf_loc, ierr)
        call VecRestoreArrayReadF90(clm_pf_idata%rate_lit3c_pfs, rate_lit3c_pf_loc, ierr)
        call VecRestoreArrayReadF90(clm_pf_idata%rate_lit1n_pfs, rate_lit1n_pf_loc, ierr)
        call VecRestoreArrayReadF90(clm_pf_idata%rate_lit2n_pfs, rate_lit2n_pf_loc, ierr)
        call VecRestoreArrayReadF90(clm_pf_idata%rate_lit3n_pfs, rate_lit3n_pf_loc, ierr)
        call VecRestoreArrayReadF90(clm_pf_idata%rate_minn_pfs,  rate_minn_pf_loc, ierr)
        call VecRestoreArrayReadF90(clm_pf_idata%rate_plantnuptake_pf, rate_plantnuptake_pf_loc, ierr)
        call VecRestoreArrayReadF90(clm_pf_idata%rate_nleached_pf, rate_nleached_pf_loc, ierr)
        call VecRestoreArrayReadF90(clm_pf_idata%rate_ndenitri_pf, rate_ndenitri_pf_loc, ierr)

    endif

    if (associated(realization%rt_mass_transfer_list)) then
       cur_mass_transfer => realization%rt_mass_transfer_list
       do
         if (.not.associated(cur_mass_transfer)) exit

         select case (cur_mass_transfer%idof)
            case(2)
             call VecGetArrayReadF90(clm_pf_idata%rate_minn_pfs, rate_pf_loc, ierr)
            case(8)
             call VecGetArrayReadF90(clm_pf_idata%rate_lit1c_pfs, rate_pf_loc, ierr)
            case(9)
             call VecGetArrayReadF90(clm_pf_idata%rate_lit2c_pfs, rate_pf_loc, ierr)
            case(10)
             call VecGetArrayReadF90(clm_pf_idata%rate_lit3c_pfs, rate_pf_loc, ierr)
            case(11)
             call VecGetArrayReadF90(clm_pf_idata%rate_lit1n_pfs, rate_pf_loc, ierr)
            case(12)
             call VecGetArrayReadF90(clm_pf_idata%rate_lit2n_pfs, rate_pf_loc, ierr)
            case(13)
             call VecGetArrayReadF90(clm_pf_idata%rate_lit3n_pfs, rate_pf_loc, ierr)
            case default
                    pflotran_model%option%io_buffer = 'Error: set PFLOTRAN BGC rates using CLM'
                    call printErrMsg(pflotran_model%option)
         end select

         do local_id = 1, grid%nlmax
            if (grid%nG2L(local_id) < 0) cycle ! bypass ghosted corner cells
            if (patch%imat(local_id) <= 0) cycle
            if(cur_mass_transfer%idof .eq. 2 .or. &
               cur_mass_transfer%idof .eq. 8 .or. &
               cur_mass_transfer%idof .eq. 9 .or. &
               cur_mass_transfer%idof .eq. 10 .or. &
               cur_mass_transfer%idof .eq. 11 .or. &
               cur_mass_transfer%idof .eq. 12 .or. &
               cur_mass_transfer%idof .eq. 13 &
               ) then
!             if(rate_pf_loc(local_id) > 0.0d0) then
               cur_mass_transfer%dataset%rarray(local_id) = &
                        rate_pf_loc(local_id)*volume_p(local_id)  ! mol/m3s * m3
!             endif
            endif

         enddo

         select case (cur_mass_transfer%idof)
            case(2)
             call VecRestoreArrayReadF90(clm_pf_idata%rate_minn_pfs, rate_pf_loc, ierr)
            case(8)
             call VecRestoreArrayReadF90(clm_pf_idata%rate_lit1c_pfs, rate_pf_loc, ierr)
            case(9)
             call VecRestoreArrayReadF90(clm_pf_idata%rate_lit2c_pfs, rate_pf_loc, ierr)
            case(10)
             call VecRestoreArrayReadF90(clm_pf_idata%rate_lit3c_pfs, rate_pf_loc, ierr)
            case(11)
             call VecRestoreArrayReadF90(clm_pf_idata%rate_lit1n_pfs, rate_pf_loc, ierr)
            case(12)
             call VecRestoreArrayReadF90(clm_pf_idata%rate_lit2n_pfs, rate_pf_loc, ierr)
            case(13)
             call VecRestoreArrayReadF90(clm_pf_idata%rate_lit3n_pfs, rate_pf_loc, ierr)
            case default
               pflotran_model%option%io_buffer = 'Error: set PFLOTRAN BGC rates using CLM'
               call printErrMsg(pflotran_model%option)
         end select

         cur_mass_transfer => cur_mass_transfer%next
       enddo
    endif

    call VecRestoreArrayReadF90(field%volume,volume_p,ierr)

  end subroutine pflotranModelSetBGCRates

    ! ************************************************************************** !
  !> This routine get updated bgc states/fluxes evoloved by PFLOTRAN.
  !!
  !> @author
  !! F.-M. Yuan
  !!
  !! date: 9/9/2013
  ! ************************************************************************** !
  subroutine pflotranModelGetBgcVariables(pflotran_model)
    use Realization_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_Aux_module
    use Immobile_Aux_module
    use Discretization_module

    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type

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

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: xx_p(:)

    PetscScalar, pointer :: decomp_cpools_vr_lit1_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_lit2_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_lit3_pf_loc(:) ! (gC/m3)
    !PetscScalar, pointer :: decomp_cpools_vr_cwd_pf_loc(:)  ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som1_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som2_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som3_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_cpools_vr_som4_pf_loc(:) ! (gC/m3)
    PetscScalar, pointer :: decomp_npools_vr_lit1_pf_loc(:) ! (gN/m3)
    PetscScalar, pointer :: decomp_npools_vr_lit2_pf_loc(:) ! (gN/m3)
    PetscScalar, pointer :: decomp_npools_vr_lit3_pf_loc(:) ! (gN/m3)
    !PetscScalar, pointer :: decomp_npools_vr_cwd_pf_loc(:)  ! (gN/m3)
    PetscScalar, pointer :: sminn_vr_pf_loc(:)              ! (gN/m3)
    PetscScalar, pointer :: hrc_vr_pf_loc(:)                ! (gN/m3)
    !PetscScalar, pointer :: smin_no3_vr_pf_loc(:)           ! (gN/m3)
    !PetscScalar, pointer :: smin_nh4_vr_pf_loc(:)           ! (gN/m3)
    PetscScalar, pointer :: accextrn_vr_pf_loc(:)            ! (gN/m3)

    PetscInt :: ispec_c, ispec_n, ispec_no3, ispec_nh4, offset, offsetim
    PetscInt :: ispec_lit1c, ispec_lit2c, ispec_lit3c
    PetscInt :: ispec_lit1n, ispec_lit2n, ispec_lit3n
    PetscInt :: ispec_som1, ispec_som2, ispec_som3, ispec_som4
    PetscInt :: ispec_plantn

    character(len=MAXWORDLENGTH) :: word
    PetscReal, parameter :: C_molecular_weight = 12.0107d0
    PetscReal, parameter :: N_molecular_weight = 14.0067d0

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

    ! (i) indices of bgc variables in PFLOTRAN immobiles
    word = "LabileC"
    ispec_lit1c = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "CelluloseC"
    ispec_lit2c = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "LigninC"
    ispec_lit3c = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "LabileN"
    ispec_lit1n = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "CelluloseN"
    ispec_lit2n = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "LigninN"
    ispec_lit3n = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "SOM1"
    ispec_som1  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "SOM2"
    ispec_som2  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "SOM3"
    ispec_som3  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)


    word = "SOM4"
    ispec_som4  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "N"
    ispec_n  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = "C"
    ispec_c  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

!    word = "NO3-"
!    ispec_no3  = GetPrimarySpeciesIDFromName(word, &
!                  realization%reaction,realization%option)

!    word = "NH4+"
!    ispec_no3  = GetPrimarySpeciesIDFromName(word, &
!                  realization%reaction,realization%option)

    word = "PlantN"
    ispec_plantn  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    ! (ii) get the original 'pf' vecs
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit1_pfp, &
                        decomp_cpools_vr_lit1_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit2_pfp, &
                        decomp_cpools_vr_lit2_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit3_pfp, &
                        decomp_cpools_vr_lit3_pf_loc, ierr)
    !call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_cwd_pfp,  &
    !                    decomp_cpools_vr_cwd_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som1_pfp, &
                        decomp_cpools_vr_som1_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som2_pfp, &
                        decomp_cpools_vr_som2_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som3_pfp, &
                        decomp_cpools_vr_som3_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som4_pfp, &
                        decomp_cpools_vr_som4_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit1_pfp, &
                        decomp_npools_vr_lit1_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit2_pfp, &
                        decomp_npools_vr_lit2_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit3_pfp, &
                        decomp_npools_vr_lit3_pf_loc, ierr)
    !call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_cwd_pfp, &
    !                    decomp_npools_vr_cwd_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%sminn_vr_pfp, sminn_vr_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hrc_vr_pfp, hrc_vr_pf_loc,ierr)
    !call VecGetArrayF90(clm_pf_idata%smin_no3_vr_pfp, smin_no3_vr_pf_loc, ierr)
    !call VecGetArrayF90(clm_pf_idata%smin_nh4_vr_pfp, smin_nh4_vr_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%accextrn_vr_pfp, accextrn_vr_pf_loc, ierr)

    ! (iii) pass the data from internal to PFLOTRAN vecs
    call VecGetArrayF90(field%tran_xx,xx_p,ierr)  ! extract data from pflotran internal portion

    do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (associated(patch%imat)) then
           if (patch%imat(ghosted_id) <= 0) cycle
        endif

        offset = (local_id - 1)*realization%reaction%ncomp

        offsetim = offset + realization%reaction%offset_immobile

        decomp_cpools_vr_lit1_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit1c), 1.0d-20) &
                                        * C_molecular_weight
        decomp_cpools_vr_lit2_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit2c), 1.0d-20) &
                                        * C_molecular_weight
        decomp_cpools_vr_lit3_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit3c), 1.0d-20) &
                                        * C_molecular_weight
        decomp_npools_vr_lit1_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit1n), 1.0d-20) &
                                        * N_molecular_weight
        decomp_npools_vr_lit2_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit2n), 1.0d-20) &
                                        * N_molecular_weight
        decomp_npools_vr_lit3_pf_loc(local_id) = max(xx_p(offsetim + ispec_lit3n), 1.0d-20) &
                                        * N_molecular_weight
        decomp_cpools_vr_som1_pf_loc(local_id) = max(xx_p(offsetim + ispec_som1), 1.0d-20) &
                                        * C_molecular_weight
        decomp_cpools_vr_som2_pf_loc(local_id) = max(xx_p(offsetim + ispec_som2), 1.0d-20) &
                                        * C_molecular_weight
        decomp_cpools_vr_som3_pf_loc(local_id) = max(xx_p(offsetim + ispec_som3), 1.0d-20) &
                                        * C_molecular_weight
        decomp_cpools_vr_som4_pf_loc(local_id) = max(xx_p(offsetim + ispec_som4), 1.0d-20) &
                                        * C_molecular_weight
        sminn_vr_pf_loc(local_id) = max(xx_p(offsetim + ispec_n), 1.0d-20) &
                                        * N_molecular_weight
        hrc_vr_pf_loc(local_id)   = max(xx_p(offsetim + ispec_c), 1.0d-20) &
                                        * C_molecular_weight
        accextrn_vr_pf_loc(local_id) = max(xx_p(offsetim + ispec_plantn), 1.0d-20) &
                                        * N_molecular_weight
    enddo

    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit1_pfp, decomp_cpools_vr_lit1_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit2_pfp, decomp_cpools_vr_lit2_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit3_pfp, decomp_cpools_vr_lit3_pf_loc, ierr)
    !call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_cwd_pfp,  decomp_cpools_vr_cwd_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som1_pfp, decomp_cpools_vr_som1_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som2_pfp, decomp_cpools_vr_som2_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som3_pfp, decomp_cpools_vr_som3_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som4_pfp, decomp_cpools_vr_som4_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit1_pfp, decomp_npools_vr_lit1_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit2_pfp, decomp_npools_vr_lit2_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit3_pfp, decomp_npools_vr_lit3_pf_loc, ierr)
    !call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_cwd_pfp,  decomp_npools_vr_cwd_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%sminn_vr_pfp, sminn_vr_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hrc_vr_pfp, hrc_vr_pf_loc, ierr)
    !call VecRestoreArrayF90(clm_pf_idata%smin_no3_vr_pfp, smin_no3_vr_pf_loc, ierr)
    !call VecRestoreArrayF90(clm_pf_idata%smin_nh4_vr_pfp, smin_nh4_vr_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%accextrn_vr_pfp, accextrn_vr_pf_loc, ierr)

   ! (iv) pass the 'pf' vecs to 'clm' vecs, which then can be passed to CLMCN (implemented in 'clm_pflotran_interfaceMod'
    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_lit1_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_lit1_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_lit2_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_lit2_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_lit3_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_lit3_clms)

    !call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%decomp_cpools_vr_cwd_pfp, &
    !                                clm_pf_idata%decomp_cpools_vr_cwd_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som1_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_som1_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som2_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_som2_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som3_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_som3_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_cpools_vr_som4_pfp, &
                                    clm_pf_idata%decomp_cpools_vr_som4_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_lit1_pfp, &
                                    clm_pf_idata%decomp_npools_vr_lit1_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_lit2_pfp, &
                                    clm_pf_idata%decomp_npools_vr_lit2_clms)
    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%decomp_npools_vr_lit3_pfp, &
                                    clm_pf_idata%decomp_npools_vr_lit3_clms)

    !call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%decomp_npools_vr_cwd_pfp, &
    !                                clm_pf_idata%decomp_npools_vr_cwd_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sminn_vr_pfp, &
                                    clm_pf_idata%sminn_vr_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%hrc_vr_pfp, &
                                    clm_pf_idata%hrc_vr_clms)

    !call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%smin_no3_vr_pfp, &
    !                                clm_pf_idata%smin_no3_vr_clms)

    !call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
    !                                pflotran_model%option, &
    !                                clm_pf_idata%smin_nh4_vr_pfp, &
    !                                clm_pf_idata%smin_nh4_vr_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%accextrn_vr_pfp, &
                                    clm_pf_idata%accextrn_vr_clms)

  end subroutine pflotranModelGetBgcVariables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  
end module pflotran_model_module

