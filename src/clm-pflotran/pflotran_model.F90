module pflotran_model_module

  use Option_module, only : option_type
  use Simulation_Base_class, only : simulation_base_type
  use Realization_Base_class, only : realization_base_type
  use Mapping_module, only : mapping_type

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"
#include "finclude/petscsysdef.h"
#include "finclude/petscviewer.h"
#include "finclude/petscvec.h"
!#include "finclude/petscvec.h90"

  ! module level constants
  PetscInt, parameter :: CLM2PF_FLUX_MAP_ID  = 1 ! 3D --> 3D
  PetscInt, parameter :: CLM2PF_SOIL_MAP_ID  = 2 ! 3D --> extended 3D
  PetscInt, parameter :: PF2CLM_FLUX_MAP_ID  = 3 ! 3D --> 3D
  PetscInt, parameter :: CLM2PF_GFLUX_MAP_ID = 4 ! 3D --> SURF-3D
  PetscInt, parameter :: CLM2PF_RFLUX_MAP_ID = 5 ! 3D --> SURF-2D
  PetscInt, parameter :: PF2CLM_SURF_MAP_ID  = 6 ! SURF-2D --> 3D

  PetscInt, parameter :: CLM_3D_MESH      = 1
  PetscInt, parameter :: PF_3D_MESH       = 2
  PetscInt, parameter :: PF_SURF_3D_MESH  = 3
  PetscInt, parameter :: PF_SURF_2D_MESH  = 4

  type, public :: inside_each_overlapped_cell
     PetscInt           :: id
     PetscInt           :: ocell_count
     PetscInt,  pointer :: ocell_id(:)
     PetscReal, pointer :: perc_vol_overlap(:)
     PetscReal          :: total_vol_overlap
  end type inside_each_overlapped_cell

  type, public :: pflotran_model_type
    class(simulation_base_type),  pointer :: simulation
    type(option_type),      pointer :: option
    PetscReal :: pause_time_1
    PetscReal :: pause_time_2
    type(inside_each_overlapped_cell), pointer :: pf_cells(:)
    type(inside_each_overlapped_cell), pointer :: clm_cells(:)
    type(mapping_type),                pointer :: map_clm2pf
    type(mapping_type),                pointer :: map_clm2pf_soils
    type(mapping_type),                pointer :: map_clm2pf_gflux
    type(mapping_type),                pointer :: map_clm2pf_rflux
    type(mapping_type),                pointer :: map_pf2clm
    type(mapping_type),                pointer :: map_pf2clm_surf
     
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

  end type pflotran_model_type

  public::pflotranModelCreate,               &
       pflotranModelInitMapping,             &
       pflotranModelSetSoilProp,             &
       pflotranModelSetICs,                  &
       pflotranModelUpdateFlowConds,         &
       pflotranModelGetUpdatedStates,        &
       pflotranModelStepperRunInit,          &
       pflotranModelStepperRunTillPauseTime, &
       pflotranModelSetupRestart,            &
       pflotranModelStepperRunFinalize,      &
       pflotranModelStepperCheckpoint,       &
       pflotranModelNSurfCells3DDomain,      &
       pflotranModelGetTopFaceArea,          &
       pflotranModelDestroy

  private :: &
       pflotranModelSetupMappingFiles


contains

  ! ************************************************************************** !
  !
  ! pflotranModelCreate: Allocates and initializes the pflotranModel object.
  !   It performs the same sequence of commands as done in pflotran.F90
  !   before model integration is performed by the call to StepperRun()
  !   routine
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  function pflotranModelCreate(mpicomm, pflotran_prefix)

    use Option_module
    use Simulation_Base_class
    use PFLOTRAN_Factory_module
    use Subsurface_Factory_module, only : SubsurfaceInitialize
    use Hydrogeophysics_Factory_module
    use Surface_Factory_module
    use Surf_Subsurf_Factory_module
  
    implicit none

#include "definitions.h"

    PetscInt, intent(in) :: mpicomm
    character(len=256), intent(in) :: pflotran_prefix

    type(pflotran_model_type), pointer :: pflotranModelCreate

    type(pflotran_model_type),      pointer :: model

    allocate(model)

    nullify(model%simulation)
    nullify(model%option)

    model%option => OptionCreate()
    call OptionInitMPI(model%option, mpicomm)
    call PFLOTRANInitialize(model%option)

    ! NOTE(bja) 2013-06-25 : external driver must provide an input
    ! prefix string. If the driver wants to use pflotran.in, then it
    ! should explicitly request that with 'pflotran'.
    if (len(trim(pflotran_prefix)) > 1) then
      model%option%input_prefix = trim(pflotran_prefix)
      model%option%input_filename = trim(model%option%input_prefix) // '.in'
    else
      model%option%io_buffer = 'The external driver must provide the ' // &
           'pflotran input file prefix.'
      call printErrMsg(model%option)
    end if

    ! TODO(bja, 2013-07-15) this needs to be left alone for pflotran
    ! to deal with, or we need a valid unit number from the driver as
    ! a function parameter.
!!$    model%option%fid_out = 16

    model%pause_time_1 = -1.0d0
    model%pause_time_2 = -1.0d0

    ! FIXME(bja, 2013-07-17) hard code subsurface for now....
    model%option%simulation_mode = 'SURFACE_SUBSURFACE'
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

    pflotranModelCreate => model

    ! NOTE(bja, 2013-07-19) GB's Hack to get communicator correctly
    ! setup on mpich/mac. should be generally ok, but may need an
    ! apple/mpich ifdef if it cause problems elsewhere.
    PETSC_COMM_SELF = MPI_COMM_SELF

  end function pflotranModelCreate


  ! ************************************************************************** !
  !
  ! pflotranModelSetupMappingFiles
  !   create the mapping objects, reopen the input file and read the file names
  !   before model integration is performed by the call to StepperRun()
  !   routine
  !
  !  NOTE(bja, 2013-07) this really needs to be moved out of pflotran
  !  CLM should be responsible for passing data in the correct
  !  format. That may require pflotran to provide a call back function
  !  for grid info.
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelSetupMappingFiles(model)

    use String_module
    use Option_module
    use Input_module
    use Mapping_module

    implicit none

#include "definitions.h"

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
    nullify(model%map_clm2pf)
    nullify(model%map_clm2pf_soils)
    nullify(model%map_clm2pf_gflux)
    nullify(model%map_clm2pf_rflux)
    nullify(model%map_pf2clm)
    nullify(model%map_pf2clm_surf)

    model%map_clm2pf       => MappingCreate()
    model%map_clm2pf_soils => MappingCreate()
    model%map_clm2pf_gflux => MappingCreate()
    model%map_clm2pf_rflux => MappingCreate()
    model%map_pf2clm       => MappingCreate()
    model%map_pf2clm_surf  => MappingCreate()

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
      call InputReadFlotranString(input, model%option)
      if (InputCheckExit(input, model%option)) exit
      if (input%ierr /= 0) exit

      call InputReadWord(input, model%option, word, PETSC_TRUE)
      call InputErrorMsg(input, model%option, 'keyword', 'MAPPING_FILES')
      call StringToUpper(word)

      select case(trim(word))
        case('CLM2PF_FLUX_FILE')
          call InputReadNChars(input, model%option, model%map_clm2pf%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_clm2pf%filename = trim(model%map_clm2pf%filename)//CHAR(0)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')   
          clm2pf_flux_file=PETSC_TRUE
        case('CLM2PF_SOIL_FILE')
          call InputReadNChars(input, model%option, model%map_clm2pf_soils%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')   
          clm2pf_soil_file=PETSC_TRUE
        case('CLM2PF_GFLUX_FILE')
          call InputReadNChars(input, model%option, model%map_clm2pf_gflux%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          clm2pf_gflux_file=PETSC_TRUE
        case('CLM2PF_RFLUX_FILE')
          call InputReadNChars(input, model%option, model%map_clm2pf_rflux%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          clm2pf_rflux_file=PETSC_TRUE
        case('PF2CLM_SURF_FILE')
          call InputReadNChars(input, model%option, model%map_pf2clm_surf%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          pf2clm_surf_file=PETSC_TRUE
        case('PF2CLM_FLUX_FILE')
          call InputReadNChars(input, model%option, model%map_pf2clm%filename, &
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
  !
  ! pflotranModelStepperRunInit: It performs the same execution of commands
  !   that are carried out in StepperRun() before the model integration
  !   begins over the entire simulation time interval
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelStepperRunInit(model)

    type(pflotran_model_type), pointer, intent(inout) :: model

    call model%simulation%InitializeRun()

  end subroutine pflotranModelStepperRunInit


  ! ************************************************************************** !
  !
  ! pflotranModelStepperCheckpoint: wrapper around StepperCheckpoint
  !
  ! NOTE(bja, 2013-06-27) : the date stamp from clm is 32 characters
  !
  ! **************************************************************************
  subroutine pflotranModelStepperCheckpoint(model, date_stamp)

    use Option_module
    use Timestepper_module, only : StepperCheckpoint

    implicit none

    type(pflotran_model_type), pointer :: model
    character(len=MAXWORDLENGTH), intent(in) :: date_stamp

    model%option%io_buffer = 'checkpoint is not implemented for clm-pflotran.'
    call printWrnMsg(model%option)
!!$    call StepperCheckpoint(model%realization, &
!!$         model%simulation%flow_stepper, &
!!$         model%simulation%tran_stepper, &
!!$         NEG_ONE_INTEGER, date_stamp)

  end subroutine pflotranModelStepperCheckpoint


  ! ************************************************************************** !
  !
  ! pflotranModelSetICs: Set initial conditions
  !
  ! author: Gautam Bisht
  ! date: 10/22/2010
  ! ************************************************************************** !
subroutine pflotranModelSetICs(pflotran_model)

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

    call MappingSourceToDestination(pflotran_model%map_clm2pf, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_clm, &
                                    clm_pf_idata%press_pf)

    if (pflotran_model%option%iflowmode .ne. RICHARDS_MODE) then
        pflotran_model%option%io_buffer='pflotranModelSetICs ' // &
          'not implmented for this mode.'
        call printErrMsg(pflotran_model%option)
    endif

    call GridVecGetArrayF90(grid, field%flow_xx, xx_loc_p, ierr)
    call VecGetArrayF90(clm_pf_idata%press_pf, press_pf_loc, ierr)

    do local_id = 1, grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
       xx_loc_p(ghosted_id)=press_pf_loc(local_id)
    enddo

    call GridVecRestoreArrayF90(grid, field%flow_xx, xx_loc_p, ierr)
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
  !
  ! pflotranModelSetSoilProp: Converts hydraulic properties from CLM units
  !  into PFLOTRAN units.
  !
  !
  ! author: Gautam Bisht
  ! date: 10/22/2010
  ! ************************************************************************** !
!#ifdef CLM_PFLOTRAN
  subroutine pflotranModelSetSoilProp(pflotran_model)

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
    PetscScalar, pointer :: bsw_clm_loc(:)    ! Clapp and Hornberger "b"

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
        call printErrMsg(pflotran_model%option, &
          'Current PFLOTRAN mode not supported by pflotranModelSetSoilProp')
    end select

    call MappingSourceToDestination(pflotran_model%map_clm2pf_soils, &
                                    pflotran_model%option, &
                                    clm_pf_idata%hksat_x_clm, &
                                    clm_pf_idata%hksat_x_pf)

    call MappingSourceToDestination(pflotran_model%map_clm2pf_soils, &
                                    pflotran_model%option, &
                                    clm_pf_idata%hksat_y_clm, &
                                    clm_pf_idata%hksat_y_pf)

    call MappingSourceToDestination(pflotran_model%map_clm2pf_soils, &
                                    pflotran_model%option, &
                                    clm_pf_idata%hksat_z_clm, &
                                    clm_pf_idata%hksat_z_pf)

    call MappingSourceToDestination(pflotran_model%map_clm2pf_soils, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sucsat_clm, &
                                    clm_pf_idata%sucsat_pf)

    call MappingSourceToDestination(pflotran_model%map_clm2pf_soils, &
                                    pflotran_model%option, &
                                    clm_pf_idata%bsw_clm, &
                                    clm_pf_idata%bsw_pf)

    call MappingSourceToDestination(pflotran_model%map_clm2pf_soils, &
                                    pflotran_model%option, &
                                    clm_pf_idata%watsat_clm, &
                                    clm_pf_idata%watsat_pf)

    call VecGetArrayF90(clm_pf_idata%hksat_x_pf, hksat_x_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_y_pf, hksat_y_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_z_pf, hksat_z_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%sucsat_pf,  sucsat_pf_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_pf,  watsat_pf_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%bsw_pf,     bsw_pf_loc,     ierr)
    call VecGetArrayF90(clm_pf_idata%bsw_clm,    bsw_clm_loc,    ierr)

    call GridVecGetArrayF90(grid, field%porosity_loc, porosity_loc_p, ierr)
    call GridVecGetArrayF90(grid, field%perm_xx_loc,  perm_xx_loc_p,  ierr)
    call GridVecGetArrayF90(grid, field%perm_yy_loc,  perm_yy_loc_p,  ierr)
    call GridVecGetArrayF90(grid, field%perm_zz_loc,  perm_zz_loc_p,  ierr)

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
      perm_xx_loc_p(local_id) = hksat_x_pf_loc(local_id)*vis/(den*grav)/1000.d0
      perm_yy_loc_p(local_id) = hksat_y_pf_loc(local_id)*vis/(den*grav)/1000.d0
      perm_zz_loc_p(local_id) = hksat_z_pf_loc(local_id)*vis/(den*grav)/1000.d0

      porosity_loc_p(local_id) = watsat_pf_loc(local_id)

    enddo

    call VecRestoreArrayF90(clm_pf_idata%hksat_x_pf, hksat_x_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y_pf, hksat_y_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z_pf, hksat_z_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat_pf,  sucsat_pf_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_pf,  watsat_pf_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw_pf,     bsw_pf_loc,     ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw_clm,    bsw_clm_loc,    ierr)

    call GridVecRestoreArrayF90(grid, field%porosity_loc, porosity_loc_p, ierr)
    call GridVecRestoreArrayF90(grid, field%perm_xx_loc,  perm_xx_loc_p,  ierr)
    call GridVecRestoreArrayF90(grid, field%perm_yy_loc,  perm_yy_loc_p,  ierr)
    call GridVecRestoreArrayF90(grid, field%perm_zz_loc,  perm_zz_loc_p,  ierr)

  end subroutine pflotranModelSetSoilProp
!#endif

  ! ************************************************************************** !
  !
  ! pflotranModelInitMapping: Initialize mapping between the two model grid
  !  (CLM and PFLTORAN)
  !
  ! author: Gautam Bisht
  ! date: 03/24/2011
  ! ************************************************************************** !
  subroutine pflotranModelInitMapping(pflotran_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id)

    use Input_module
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
!    use clm_pflotran_interface_data

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

    if(map_id==CLM2PF_GFLUX_MAP_ID) then
      call pflotranModelInitMappingSurf3D(pflotran_model,  &
                                          grid_clm_cell_ids_nindex, &
                                          grid_clm_npts_local, &
                                          map_id)
      return
    endif

    if(map_id==CLM2PF_RFLUX_MAP_ID .or. map_id==PF2CLM_SURF_MAP_ID) then
      call pflotranModelInitMappingSurf2D(pflotran_model,  &
                                          grid_clm_cell_ids_nindex, &
                                          grid_clm_npts_local, &
                                          map_id)
      return
    endif

    !
    ! Mapping to/from entire PFLOTRAN 3D subsurface domain
    !

    ! Choose the appriopriate map
    select case(map_id)
      case(CLM2PF_FLUX_MAP_ID)
        map => pflotran_model%map_clm2pf
        source_mesh_id = CLM_3D_MESH
        dest_mesh_id = PF_3D_MESH
      case(CLM2PF_SOIL_MAP_ID)
        map => pflotran_model%map_clm2pf_soils
        source_mesh_id = CLM_3D_MESH
        dest_mesh_id = PF_3D_MESH
      case(PF2CLM_FLUX_MAP_ID)
        map => pflotran_model%map_pf2clm
        source_mesh_id = PF_3D_MESH
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
      case(CLM_3D_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_cell_ids_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_3D_MESH)
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

  end subroutine pflotranModelInitMapping

  ! ************************************************************************** !
  !> This routine maps CLM surface grid onto surface of PFLOTRAN 3D subsurface 
  !! grid.
  !!
  !> @author
  !! Gautam Bisht, LBNL
  !!
  !! date: 04/09/13
  ! ************************************************************************** !
  subroutine pflotranModelInitMappingSurf3D(pflotran_model,  &
                                            grid_clm_cell_ids_nindex, &
                                            grid_clm_npts_local, &
                                            map_id)

    use Input_module
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
      case(CLM2PF_GFLUX_MAP_ID)
        map => pflotran_model%map_clm2pf_gflux
        source_mesh_id = CLM_3D_MESH
        dest_mesh_id = PF_SURF_3D_MESH
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

      case(PF_SURF_3D_MESH)

        patch => realization%patch
        grid => patch%grid

        ! Destination mesh is PF_SURF_3D_MESH
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
        clm_pf_idata%nlclm_surf_3d = grid_clm_npts_local
        clm_pf_idata%ngclm_surf_3d = grid_clm_npts_local
        clm_pf_idata%nlpf_surf_3d  = grid_pf_npts_local
        clm_pf_idata%ngpf_surf_3d  = grid_pf_npts_local

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
    call VecCreateMPI(option%mycomm, clm_pf_idata%nlclm_3d, PETSC_DECIDE, surf_ids, ierr)
    call VecSet(surf_ids, -1.d0, ierr)

    ! Set 1.0 to all cells that make up surface of CLM subsurface domain
    call VecSetValues(surf_ids, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                      v_loc, INSERT_VALUES, ierr)

    deallocate(v_loc)
    call VecAssemblyBegin(surf_ids, ierr)
    call VecAssemblyEnd(surf_ids, ierr)

    call VecGetArrayF90(surf_ids, v_loc, ierr)
    count = 0
    do local_id=1,clm_pf_idata%nlclm_3d
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr)

    count = 0
    do local_id=1,clm_pf_idata%nlclm_3d
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
      case(CLM_3D_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_cell_ids_nindex_copy)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_3D_MESH)
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

  end subroutine pflotranModelInitMappingSurf3D

  ! ************************************************************************** !
  !> This routine maps CLM surface grid onto PFLOTRAN 2D surface grid.
  !!
  !> @author
  !! Gautam Bisht, LBNL
  !!
  !! date: 04/09/13
  ! ************************************************************************** !
  subroutine pflotranModelInitMappingSurf2D(pflotran_model,  &
                                            grid_clm_cell_ids_nindex, &
                                            grid_clm_npts_local, &
                                            map_id)

    use Input_module
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
    option%io_buffer='To support dest_mesh == PF_SURF_2D_MESH, need to '// &
         'compiled with -DSURFACE_FLOW.'
    call printErrMsg(option)
#else
    allocate(grid_clm_cell_ids_nindex_copy(grid_clm_npts_local))
    grid_clm_cell_ids_nindex_copy = grid_clm_cell_ids_nindex

    ! Choose the appriopriate map
    select case(map_id)
      case(CLM2PF_RFLUX_MAP_ID)
        map => pflotran_model%map_clm2pf_rflux
        source_mesh_id = CLM_3D_MESH
        dest_mesh_id = PF_SURF_2D_MESH
      case(PF2CLM_SURF_MAP_ID)
        map => pflotran_model%map_pf2clm_surf
        source_mesh_id = PF_SURF_2D_MESH
        dest_mesh_id = CLM_3D_MESH
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
                      patch%grid%nlmax, &
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
    call VecCreateMPI(option%mycomm, clm_pf_idata%nlclm_3d, PETSC_DECIDE, surf_ids, ierr)
    call VecSet(surf_ids, -1.d0, ierr)

    ! Set 1.0 to all cells that make up surface of CLM subsurface domain
    call VecSetValues(surf_ids, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                      v_loc, INSERT_VALUES, ierr)

    deallocate(v_loc)
    call VecAssemblyBegin(surf_ids, ierr)
    call VecAssemblyEnd(surf_ids, ierr)

    call VecGetArrayF90(surf_ids, v_loc, ierr)
    count = 0
    do local_id=1,clm_pf_idata%nlclm_3d
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr)

    count = 0
    do local_id=1,clm_pf_idata%nlclm_3d
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
    
    if(count /= map%s2d_nwts) then
      write(*,*),'count = ',option%myrank,count,map%s2d_nwts
      option%io_buffer='No. of surface cells in mapping dataset does not ' // &
        'match surface cells on which BC is applied. [pflotranModelInitMappingSurf2D]'
      call printErrMsgByRank(option)
    endif
    call VecDestroy(surf_ids, ierr)

    select case(source_mesh_id)
      case(CLM_3D_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_cell_ids_nindex_copy)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_SURF_2D_MESH)
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

  end subroutine pflotranModelInitMappingSurf2D

  ! ************************************************************************** !
  !
  ! pflotranModelStepperRunTillPauseTime: It performs the model integration
  !             till the specified pause_time.
  !
  ! NOTE: It is assumed 'pause_time' is in seconds
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelStepperRunTillPauseTime(model, pause_time)


    implicit none

#include "definitions.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer :: model
    PetscReal, intent(in) :: pause_time

    call model%simulation%RunToTime(pause_time)

    ! TODO(GB): Use XXXUpdateMassBalancePatch() to ensure mass balance
    ! betweent CLM calls

  end subroutine pflotranModelStepperRunTillPauseTime

  ! ************************************************************************** !
  !
  !  pflotranModelSetupRestart()
  !
  !  This checks to see if a restart file stamp was provided by the
  !  driver. If so, we set the restart flag and reconstruct the
  !  restart file name. The actual restart is handled by the standard
  !  pflotran mechanism in TimeStepperInitializeRun()
  !
  !  NOTE: this must be called between pflotranModelCreate() and
  !  pflotranModelStepperRunInit()
  !
  ! ************************************************************************** !
  subroutine pflotranModelSetupRestart(model, restart_stamp)

    use Option_module
    use String_module

    implicit none

    type(pflotran_model_type), pointer :: model
    character(len=MAXWORDLENGTH) :: restart_stamp

    model%option%io_buffer = 'restart is not implemented in clm-pflotran.'
    call printWrnMsg(model%option)

!!$    if (.not. StringNull(restart_stamp)) then
!!$       model%option%restart_flag = PETSC_TRUE
!!$       model%option%restart_filename = &
!!$            trim(model%option%global_prefix) // &
!!$            trim(model%option%group_prefix) // &
!!$            '.' // trim(restart_stamp) // '.chk'
!!$    end if

  end subroutine pflotranModelSetupRestart

  ! ************************************************************************** !
  ! pflotranModelUpdateSourceSink: Update the source/sink term
  !
  ! author: Gautam Bisht
  ! date: 11/22/2011
  ! ************************************************************************** !
  subroutine pflotranModelUpdateSourceSink(pflotran_model)

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model

    call MappingSourceToDestination(pflotran_model%map_clm2pf, &
                                    pflotran_model%option, &
                                    clm_pf_idata%qflx_clm, &
                                    clm_pf_idata%qflx_pf)

  end subroutine pflotranModelUpdateSourceSink

  ! ************************************************************************** !
  !> This routine Updates boundary and source/sink condtions for PFLOTRAN that
  !! are prescribed by CLM
  !!
  !> @author
  !! Gautam Bisht, LBNL
  !!
  !! date: 4/10/2013
  ! ************************************************************************** !
  subroutine pflotranModelUpdateFlowConds(pflotran_model)

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model

    call pflotranModelUpdateSourceSink(pflotran_model)

    if(pflotran_model%option%iflowmode==TH_MODE) then
      call MappingSourceToDestination(pflotran_model%map_clm2pf_gflux, &
                                      pflotran_model%option, &
                                      clm_pf_idata%gflux_clm, &
                                      clm_pf_idata%gflux_pf)
    endif

  end subroutine pflotranModelUpdateFlowConds

  ! ************************************************************************** !
  !> This routine get updated states evoloved by PFLOTRAN.
  !!
  !> @author
  !! Gautam Bisht, LBNL
  !!
  !! date: 5/14/2013
  ! ************************************************************************** !
  subroutine pflotranModelGetUpdatedStates(pflotran_model)

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
        pflotran_model%option%io_buffer='pflotranModelGetSaturation ' // &
          'not implmented for this mode.'
        call printErrMsg(pflotran_model%option)
    end select

  end subroutine pflotranModelGetUpdatedStates


  ! ************************************************************************** !
  ! pflotranModelGetSaturation: Extract soil saturation values simulated by 
  !   PFLOTRAN in a PETSc vector.
  !
  ! author: Gautam Bisht
  ! date: 11/22/2011
  ! ************************************************************************** !
  subroutine pflotranModelGetSaturation(pflotran_model)

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
    PetscReal, pointer :: sat_clm_p(:)

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

    call MappingSourceToDestination(pflotran_model%map_pf2clm, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sat_pf, &
                                    clm_pf_idata%sat_clm)

  end subroutine pflotranModelGetSaturation

  ! ************************************************************************** !
  !> This routine get updated states evoloved by PFLOTRAN.
  !!
  !> @author
  !! Gautam Bisht, LBNL
  !!
  !! date: 5/14/2013
  ! ************************************************************************** !
  subroutine pflotranModelGetTemperature(pflotran_model)

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

    do ghosted_id=1,grid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (local_id>0) then
        call VecSetValues(clm_pf_idata%temp_pf,1,local_id-1, &
                        global_aux_vars(ghosted_id)%temp(1),INSERT_VALUES,ierr)
        !call VecSetValues(clm_pf_idata%sat_ice_pf,1,local_id-1, &
        !                 global_aux_vars(ghosted_id)%sat_ice(1),INSERT_VALUES,ierr)
      endif
    enddo

    call VecAssemblyBegin(clm_pf_idata%temp_pf,ierr)
    call VecAssemblyEnd(clm_pf_idata%temp_pf,ierr)
    call MappingSourceToDestination(pflotran_model%map_pf2clm, &
                                    pflotran_model%option, &
                                    clm_pf_idata%temp_pf, &
                                    clm_pf_idata%temp_clm)

!    call VecAssemblyBegin(clm_pf_idata%sat_ice_pf,ierr)
!    call VecAssemblyEnd(clm_pf_idata%sat_ice_pf,ierr)
!    call MappingSourceToDestination(pflotran_model%map_pf2clm, &
!                                    pflotran_model%option, &
!                                    clm_pf_idata%sat_ice_pf, &
!                                    clm_pf_idata%sat_ice_clm)

  end subroutine pflotranModelGetTemperature

  ! ************************************************************************** !
  !
  ! pflotranModelStepperRunFinalize: It performs the same execution of commands
  !             that are carried out in StepperRun() once the model integration is
  !             finished
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelStepperRunFinalize(model)

    implicit none

    type(pflotran_model_type), pointer :: model

    call model%simulation%FinalizeRun()

  end subroutine pflotranModelStepperRunFinalize



  ! ************************************************************************** !
  !> This function returns the number of control volumes forming surface of
  !! the sub-surface domain. It assumes a boundary condition named 'clm_gflux_bc'
  !! is defined in the inputdeck.
  !!
  !> @author
  !! Gautam Bisht, LBNL
  !!
  !! date: 6/03/2013
  ! ************************************************************************** !
  function pflotranModelNSurfCells3DDomain(pflotran_model)

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

    coupler_list => realization%patch%boundary_conditions
    coupler => coupler_list%first
    found = PETSC_FALSE

    do
      if (.not.associated(coupler)) exit
      if(StringCompare(coupler%name,'clm_gflux_bc')) then
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
  !> This subroutine
  !!
  !> @author
  !! Gautam Bisht, LBNL
  !!
  !! date: 6/10/2013
  ! ************************************************************************** !
  subroutine pflotranModelGetTopFaceArea(pflotran_model)

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
    PetscInt :: cell_type
    PetscInt :: vertex_ids(4)

    PetscReal :: v1(3), v2(3), v3(3), n1(3), n2(3), n_up_dn(3)
    PetscReal :: vcross(3), magnitude
    PetscReal :: area1, area2

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

    if(grid%discretization_itype == STRUCTURED_GRID) then
      ! Structured grid
      do ghosted_id=1,grid%ngmax
        local_id = grid%nG2L(ghosted_id)
        if(local_id>0) then
          area1 = grid%structured_grid%dx(ghosted_id)* &
                  grid%structured_grid%dy(ghosted_id)
          call VecSetValues(clm_pf_idata%area_top_face_pf,1,local_id-1, &
                        area1,INSERT_VALUES,ierr)
        endif
      enddo
    else if (grid%discretization_itype == UNSTRUCTURED_GRID) then
      ! Unstructured grid
      do local_id=1,grid%nlmax
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

        ! Get vertex_id
        call UCellGetFaceVertices(option,cell_type,iface,vertex_ids)

        point1 = grid%unstructured_grid%vertices(vertex_ids(1))
        point2 = grid%unstructured_grid%vertices(vertex_ids(2))
        point3 = grid%unstructured_grid%vertices(vertex_ids(3))

        v1(1) = point3%x-point2%x
        v1(2) = point3%y-point2%y
        v1(3) = point3%z-point2%z
        v2(1) = point1%x-point2%x
        v2(2) = point1%y-point2%y
        v2(3) = point1%z-point2%z
        !geh: area = 0.5 * |v1 x v2|
        vcross = CrossProduct(v1,v2)
        !geh: but then we have to project the area onto the vector between
        !     the cell centers (n_up_dn)
        magnitude = sqrt(DotProduct(vcross,vcross))
        n1 = vcross/magnitude
        area1 = 0.5d0*magnitude
        area1 = dabs(area1*DotProduct(n1,n_up_dn))

        if(cell_type == HEX_TYPE) then
          point4 = grid%unstructured_grid%vertices(vertex_ids(4))
          v1(1) = point1%x-point4%x
          v1(2) = point1%y-point4%y
          v1(3) = point1%z-point4%z
          v2(1) = point3%x-point4%x
          v2(2) = point3%y-point4%y
          v2(3) = point3%z-point4%z
          magnitude = sqrt(DotProduct(vcross,vcross))
          n2 = vcross/magnitude
          area2 = 0.5d0*magnitude
          area2 = dabs(area2*DotProduct(n2,n_up_dn))
        else
          area2 = 0.0d0
        endif

        call VecSetValues(clm_pf_idata%area_top_face_pf,1,local_id-1, &
                       area1+area2,INSERT_VALUES,ierr)
      enddo
    endif

    call VecAssemblyBegin(clm_pf_idata%area_top_face_pf,ierr)
    call VecAssemblyEnd(clm_pf_idata%area_top_face_pf,ierr)
    call MappingSourceToDestination(pflotran_model%map_pf2clm, &
                                    pflotran_model%option, &
                                    clm_pf_idata%area_top_face_pf, &
                                    clm_pf_idata%area_top_face_clm)

  end subroutine pflotranModelGetTopFaceArea

  ! ************************************************************************** !
  !
  ! pflotranModelDestroy: Deallocates the pflotranModel object
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelDestroy(model)

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
  
    call PFLOTRANFinalize(model%option)
    call OptionFinalize(model%option)

    deallocate(model)

  end subroutine pflotranModelDestroy
  
end module pflotran_model_module

