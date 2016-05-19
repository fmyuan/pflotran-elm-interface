module pflotran_clm_main_module

  use Option_module, only : option_type
  use Simulation_Base_class, only : simulation_base_type
  use Multi_Simulation_module, only : multi_simulation_type
  use Realization_Base_class, only : realization_base_type
  use Mapping_module, only : mapping_type

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petsclog.h"
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscvec.h"

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
       ! soil domain
       pflotranModelSetSoilDimension,        &
       ! Soil properties
       pflotranModelSetSoilProp,              &
       pflotranModelGetSoilPropFromPF
       ! T/H
       ! BGC

  private :: &
       pflotranModelInsertWaypoint,          &
       pflotranModelDeleteWaypoint

!------------------------------------------------------------

  character(len=MAXWORDLENGTH) :: subname
!------------------------------------------------------------

contains

! ************************************************************************************ !

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
    use PFLOTRAN_Constants_module
    use Output_Aux_module, only : INSTANTANEOUS_VARS
    use PFLOTRAN_Provenance_module, only : PrintProvenanceToScreen
  
    implicit none

#include "petsc/finclude/petscsys.h"

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

    ! NOTE(bja, 2013-07-19) GB's Hack to get communicator correctly
    ! setup on mpich/mac. should be generally ok, but may need an
    ! apple/mpich ifdef if it cause problems elsewhere.
    PETSC_COMM_SELF = MPI_COMM_SELF
    PETSC_COMM_WORLD = MPI_COMM_WORLD

    call PFLOTRANInitializePostPetsc(model%simulation, model%multisimulation, model%option)

    model%pause_time_1 = -1.0d0
    model%pause_time_2 = -1.0d0

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
  !
  ! pflotranModelUpdateFinalWaypoint:
  !  Get CLM final timestep and converts it to PFLOTRAN final way point.
  !  And also set an option for turning on/off PF printing

  subroutine pflotranModelUpdateFinalWaypoint(model, waypoint_time, waypoint_dtmax, isprintout)

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use Realization_Subsurface_class, only : realization_subsurface_type, &
      RealizationAddWaypointsToList

    use Waypoint_module, only : waypoint_type, waypoint_list_type, &
      WaypointCreate, &
      WaypointDeleteFromList, WaypointInsertInList, &
      WaypointListFillIn, WaypointListRemoveExtraWaypnts
    use Units_module, only : UnitsConvertToInternal
    use Option_module, only : printErrMsg

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer :: model
    PetscReal, intent(in)              :: waypoint_time, waypoint_dtmax
    PetscBool, intent(in)              :: isprintout

    type(waypoint_list_type), pointer  :: waypoint_list
    type(waypoint_type), pointer       :: waypoint, waypoint1
    character(len=MAXWORDLENGTH)       :: word, internal_units

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

!-------------------------------------------------------------------------
    option => model%option
    select type (modsim => model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modsim
        realization => simulation%realization
        waypoint_list => simulation%waypoint_list_subsurface

      class default
        option%io_buffer = "pflotranModelUPdateFinalWaypoint is " // &
              "Not support in this mode."
        call printErrMsg(option)
    end select

    ! new final waypoint
    word = 's'
    internal_units = 'sec'
    waypoint1 => WaypointCreate()
    waypoint1%time          = waypoint_time * UnitsConvertToInternal(word, internal_units, option)
    waypoint1%print_snap_output  = PETSC_TRUE
    waypoint1%final         = PETSC_TRUE
    waypoint1%dt_max        = waypoint_dtmax * UnitsConvertToInternal(word, internal_units, option)

    ! update subsurface-realization final waypoint
    if (associated(realization) .and. associated(simulation)) then
      ! remove original final waypoint
      waypoint => waypoint_list%first
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
      call WaypointInsertInList(waypoint1, waypoint_list)
      call RealizationAddWaypointsToList(realization, waypoint_list)
      call WaypointListFillIn(waypoint_list, option)
      call WaypointListRemoveExtraWaypnts(waypoint_list, option)

      ! turn off the 'print out' if required from CLM
      if(.not.isprintout) then
        if (model%option%io_rank == model%option%myrank) then
          write(model%option%fid_out, *) 'NOTE: h5 output at input-defined interval ' // &
            'for subsurface flow from PFLOTRAN IS OFF! '
        endif
        waypoint => waypoint_list%first
        do
          if (.not.associated(waypoint)) exit
          waypoint%print_snap_output = PETSC_FALSE
          waypoint => waypoint%next
        enddo
      endif

    endif

  end subroutine pflotranModelUpdateFinalWaypoint

  !-------------------------------------------------------------------!

  subroutine pflotranModelInsertWaypoint(model, waypoint_time, waypoint_dtmax, waypoint_final)
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
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use Realization_Base_class, only : realization_base_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use Waypoint_module, only : waypoint_type, WaypointCreate, WaypointInsertInList
    use Units_module, only : UnitsConvertToInternal
    use Option_module, only : printErrMsg

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer :: model
    PetscReal, intent(in)              :: waypoint_time
    PetscReal, intent(in)              :: waypoint_dtmax
    PetscBool, intent(in)              :: waypoint_final

    type(waypoint_type), pointer       :: waypoint
    character(len=MAXWORDLENGTH) :: word, internal_units

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    !------------------------
    option => model%option

    word = 's'
    internal_units = 'sec'
    waypoint => WaypointCreate()
    waypoint%time              = waypoint_time * UnitsConvertToInternal(word, internal_units, option)
    waypoint%update_conditions = PETSC_TRUE
    waypoint%print_snap_output = PETSC_FALSE
    waypoint%print_obs_output  = PETSC_FALSE
    waypoint%print_checkpoint  = PETSC_FALSE
    waypoint%final             = waypoint_final
    waypoint%dt_max            = waypoint_dtmax * UnitsConvertToInternal(word, internal_units, option)

    select type (modsim => model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modsim
        realization => simulation%realization

        call WaypointInsertInList(waypoint, simulation%waypoint_list_subsurface)

      class default
         model%option%io_buffer = "pflotranModelInsertWaypoint only " // &
              "works on subsurface simulations."
         call printErrMsg(model%option)

    end select

  end subroutine pflotranModelInsertWaypoint

  !-------------------------------------------------------------------!

  subroutine pflotranModelDeleteWaypoint(model, waypoint_time)

    use Option_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use Realization_Base_class, only : realization_base_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use Waypoint_module, only : waypoint_type, WaypointCreate, WaypointDeleteFromList
    use Units_module, only : UnitsConvertToInternal

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer :: model
    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    type(waypoint_type), pointer       :: waypoint
    PetscReal                          :: waypoint_time
    character(len=MAXWORDLENGTH) :: word, internal_units

    ! ------------------------------------------ !
    option => model%option
    select type (modsim => model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modsim
        realization => simulation%realization

      class default
         option%io_buffer = "pflotranModelDeleteWaypoint only " // &
              "works on subsurface simulations."
         call printErrMsg(option)
    end select
    !
    word = 's'
    internal_units = 'sec'
    waypoint => WaypointCreate()
    waypoint%time              = waypoint_time * UnitsConvertToInternal(word, internal_units, option)
    waypoint%update_conditions = PETSC_TRUE
    waypoint%dt_max            = 86400.0d0

    if (associated(realization) .and. associated(simulation)) then
       call WaypointDeleteFromList(waypoint, simulation%waypoint_list_subsurface)
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

    type(Option_type), pointer :: option
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

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h90"

    type(Option_type), pointer :: option
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

    pause_time1 = pause_time + dtime
    !call pflotranModelUpdateFinalWaypoint(model, pause_time1, dtime, PETSC_FALSE)
    call pflotranModelInsertWaypoint(model, pause_time1, dtime, PETSC_FALSE)

    call model%simulation%RunToTime(pause_time)

    call pflotranModelDeleteWaypoint(model, pause_time)

  end subroutine pflotranModelStepperRunTillPauseTime

! ************************************************************************** !

  subroutine pflotranModelStepperCheckpoint(model, id_stamp)
  !
  ! wrapper around StepperCheckpoint
  ! NOTE(bja, 2013-06-27) : the date stamp from clm is 32 characters
  !

    use Option_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer :: model
    character(len=MAXWORDLENGTH), intent(in) :: id_stamp
    PetscViewer :: viewer

    if (associated(model%simulation%process_model_coupler_list%checkpoint_option)) then
      call model%simulation%process_model_coupler_list%Checkpoint(id_stamp)
    endif

  end subroutine pflotranModelStepperCheckpoint

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

    type(Option_type), pointer :: option
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

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer :: model

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

  end subroutine pflotranModelDestroy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!  Soil dimension/physical properties
!  (normally NOT variable, unless specified)
!
! ************************************************************************** !

  subroutine pflotranModelSetSoilDimension(pflotran_model)
  !
  ! Force soil dimension from CLM to PFLOTRAN subsurface domain.

    use Discretization_module
    use Patch_module
    use Field_module
    use Grid_module
    use Grid_Structured_module
    use Option_module
    use Material_module
    use Material_Aux_class

    !use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    !use Realization_Base_class, only : realization_base_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use clm_pflotran_interface_data
    use Mapping_module
    use Saturation_Function_module

    use Variables_module, only : VOLUME

    implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "geodesic.inc"

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(discretization_type), pointer        :: discretization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization


    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id, i, j, k

    PetscScalar, pointer :: dlon_pf_loc(:),dlat_pf_loc(:), dzsoil_pf_loc(:)           ! soil cell length/width/thickness (deg/deg/m) for 3-D PF cells
    PetscScalar, pointer :: lonc_pf_loc(:),latc_pf_loc(:), zisoil_pf_loc(:)           ! soil cell coordinates (deg/deg/m) for 3-D PF cells
    PetscScalar, pointer :: topface_pf_loc(:)
    PetscReal            :: lon_c, lat_c, lon_e, lon_w, lat_s, lat_n
    PetscReal            :: x_global, y_global
    PetscReal            :: tempreal

    ! for calling functions in 'geodesic.for'
    double precision a, f, lat1, lon1, azi1, lat2, lon2, azi2, s12,   &
        dummy1, dummy2, dummy3, dummy4, dummy5
    double precision lats(4), lons(4)
    integer omask
    a = 6378137.0d0           ! major-axis length of Earth Ellipsoid in metres in WGS-84
    f = 1.d0/298.257223563d0  ! flatening of Earth Ellipsoid in WGS-84
    omask = 0

    subname = 'pflotranModelSetSoilDimension'

!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select

    patch           => realization%patch
      grid          => patch%grid
    field           => realization%field
    discretization  => realization%discretization

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%dzsoil_clmp, &
                                    clm_pf_idata%dzsoil_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%zisoil_clmp, &
                                    clm_pf_idata%zisoil_pfs)
    call VecGetArrayReadF90(clm_pf_idata%dzsoil_pfs, dzsoil_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayReadF90(clm_pf_idata%zisoil_pfs, zisoil_pf_loc, ierr)
    CHKERRQ(ierr)

    !
    if ( (clm_pf_idata%nxclm_mapped >= 1 .and. clm_pf_idata%nyclm_mapped >= 1) .and. &
         (.not.option%mapping_files) ) then

      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%dxsoil_clmp, &
                                    clm_pf_idata%dxsoil_pfs)

      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%dysoil_clmp, &
                                    clm_pf_idata%dysoil_pfs)

      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%xsoil_clmp, &
                                    clm_pf_idata%xsoil_pfs)

      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%ysoil_clmp, &
                                    clm_pf_idata%ysoil_pfs)

    !
      call VecGetArrayReadF90(clm_pf_idata%dxsoil_pfs, dlon_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecGetArrayReadF90(clm_pf_idata%dysoil_pfs, dlat_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecGetArrayReadF90(clm_pf_idata%xsoil_pfs, lonc_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecGetArrayReadF90(clm_pf_idata%ysoil_pfs, latc_pf_loc, ierr)
      CHKERRQ(ierr)
    end if

    !
    select case(grid%itype)
      case(STRUCTURED_GRID)
         option%io_buffer = "INFO: CLM column dimension will over-ride PF structured grid mesh."
         call printMsg(option)

         select case(grid%structured_grid%itype)
           case(CARTESIAN_GRID)
             option%io_buffer = "INFO: CLM column dimension will over-ride PF structured CARTESIAN_GRID."
             call printMsg(option)

             if (option%io_rank == option%myrank) then
               write(option%fid_out, &
              '(/," Requested processors and decomposition = ", i5,", npx,y,z= ",3i4)') &
               option%mycommsize, &
               grid%structured_grid%npx, &
               grid%structured_grid%npy, &
               grid%structured_grid%npz
               write(option%fid_out,'(" Actual decomposition: npx,y,z= ",3i4,/)') &
               grid%structured_grid%npx_final, &
               grid%structured_grid%npy_final, &
               grid%structured_grid%npz_final
             endif

           case default
             option%io_buffer = "ERROR: Currently only works on structured  CARTESIAN_GRID mesh."
             call printErrMsg(option)
         end select

      case default
         option%io_buffer = "ERROR: Currently only works on structured grids."
         call printErrMsg(option)
    end select

    do ghosted_id = 1, grid%ngmax
      local_id = grid%nG2L(ghosted_id)

      ! hack cell's 3-D dimensions

      ! adjusting (x,y) if runs with 2D CLM grid (usually in lat/lon paired grid)
      if ( (clm_pf_idata%nxclm_mapped >= 1 .and. clm_pf_idata%nyclm_mapped >= 1) .and. &
           (.not.option%mapping_files) ) then

        call StructGridGetIJKFromGhostedID(grid%structured_grid,ghosted_id,i,j,k)

        !        4--(lonc, latc+1/2*dy)-3              ^
        !       /            |           \             |
        !      /             |            \            |
        !     +--------(lonc, latc)--------+     dyg_local(j)
        !    <--------dxg_local(i)---------->          |
        !   /                |               \         |
        !  /                 |                \        |
        ! 1---------(lonc, latc-1/2*dy)--------2       V

        lon_c = lonc_pf_loc(ghosted_id)
        lat_c = latc_pf_loc(ghosted_id)
        if(dlon_pf_loc(ghosted_id) >0.d0 .and. dlat_pf_loc(ghosted_id) >0.d0) then
          ! the following assumes an isoceles trapezoid grid from CLM unstructured or 1-D gridcells
          ! may have distortion for area estimation, but should be good for distance calculation
          lon_e = lon_c + dlon_pf_loc(ghosted_id)/2.0d0  ! East
          lon_w = lon_c - dlon_pf_loc(ghosted_id)/2.0d0  ! West
          lat_s = lat_c - dlat_pf_loc(ghosted_id)/2.0d0  ! South
          lat_n = lat_c + dlat_pf_loc(ghosted_id)/2.0d0  ! North

        else
          lon_e = lon_c + grid%structured_grid%dxg_local(i)/2.0d0  ! East
          lon_w = lon_c - grid%structured_grid%dxg_local(i)/2.0d0  ! West
          lat_s = lat_c - grid%structured_grid%dyg_local(j)/2.0d0  ! South
          lat_n = lat_c + grid%structured_grid%dyg_local(j)/2.0d0  ! North
        endif

        lats(1) = lat_s
        lons(1) = lon_w
        lats(2) = lat_s
        lons(2) = lon_e
        lats(3) = lat_n
        lons(3) = lon_e
        lats(4) = lat_n
        lons(4) = lon_w

        ! mid-longitudal length of trapezoid (x-axis) -
        lat1 = lat_c
        lon1 = lon_w
        lat2 = lat_c
        lon2 = lon_e
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        grid%structured_grid%dx(ghosted_id) = s12

        ! mid-latitudal height of trapezoid (y-axis) -
        lat1 = lat_s
        lon1 = lon_c
        lat2 = lat_n
        lon2 = lon_c
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        grid%structured_grid%dy(ghosted_id) = s12

        ! some checking
        ! areas of grid (x,y)
        call area(a, f, lats, lons, 4, dummy1, dummy2)
        tempreal = grid%structured_grid%dx(ghosted_id)*grid%structured_grid%dy(ghosted_id)/dummy1
        if (abs(tempreal-1.d0)>1.e-5) then
          option%io_buffer = "Warning: remarkably large gaps in grid areas btw two approaches FOR cell: "
          call printMsg(option)
        end if

        ! bottom/top segment line length
        lat1 = lats(1)
        lon1 = lons(1)
        lat2 = lats(2)
        lon2 = lons(2)
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        tempreal = s12

        lat1 = lats(3)
        lon1 = lons(3)
        lat2 = lats(4)
        lon2 = lons(4)
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        tempreal = tempreal+s12
        tempreal = 0.5d0*tempreal/grid%structured_grid%dx(ghosted_id)
        if (abs(tempreal-1.d0)>1.e-5) then   ! mathematically, dx = 0.5*(a+b)
          option%io_buffer = "Warning: remarkably large gaps in longitudal-length FOR a cell: "
          call printMsg(option)
        end if

        ! isoscele side line length
        lat1 = lats(2)
        lon1 = lons(2)
        lat2 = lats(3)
        lon2 = lons(3)
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        tempreal = s12

        lat1 = lats(1)
        lon1 = lons(1)
        lat2 = lats(4)
        lon2 = lons(4)
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        tempreal = tempreal/s12
        if (abs(tempreal-1.d0)>1.e-5) then   ! mathematically, c=d
          option%io_buffer = "Warning: remarkably large gaps in isoscele latitudal-length FOR a cell: "
          call printMsg(option)
        end if

      end if ! if (clm_pf_idata%nxclm_mapped >= 1 .and. clm_pf_idata%nyclm_mapped >= 1 .and. .not.mapping_files)

      ! vertical (z-axis) (always from CLM to PF)
      grid%structured_grid%dz(ghosted_id) = dzsoil_pf_loc(ghosted_id)
      grid%z(ghosted_id)   = zisoil_pf_loc(ghosted_id)    ! directly over-ride PF 'z' coordinate from CLM soil layer 'zi'
    enddo

    if ( (clm_pf_idata%nxclm_mapped >= 1 .and. clm_pf_idata%nyclm_mapped >= 1) .and. &
         (.not.option%mapping_files) ) then
      call VecRestoreArrayReadF90(clm_pf_idata%dxsoil_pfs, dlon_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayReadF90(clm_pf_idata%dysoil_pfs, dlat_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayReadF90(clm_pf_idata%xsoil_pfs, lonc_pf_loc, ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayReadF90(clm_pf_idata%ysoil_pfs, latc_pf_loc, ierr)
      CHKERRQ(ierr)
    end if

    call VecRestoreArrayReadF90(clm_pf_idata%dzsoil_pfs, dzsoil_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayReadF90(clm_pf_idata%zisoil_pfs, zisoil_pf_loc, ierr)
    CHKERRQ(ierr)

    ! re-do some dimension calculations after changes above
    call MPI_Barrier(option%mycomm,ierr)

    call GridComputeVolumes(grid, field%volume0, option)      ! cell volumes
    call GridComputeInternalConnect(grid, option)             ! cell internal connection distances
    call PatchProcessCouplers(patch,realization%flow_conditions,             &  ! BC/IC/SrcSink connection (face) areas
                              realization%transport_conditions,              &
                              realization%option)

    ! re-assign updated field%volume0 to material_auxvar%volume
    call DiscretizationGlobalToLocal(discretization,field%volume0, &
                                   field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material, &
                               field%work_loc,VOLUME,ZERO_INTEGER)

    ! the following variable is directly used in 'sandbox_somdec'
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%zsoil_clmp, &
                                    clm_pf_idata%zsoil_pfs)

    ! the following are 'wtgcell' and 'landfrac' adjusted 'TOP' face area (not yet figured out how to use it for multiple columns in ONE grid)
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%area_top_face_clmp, &
                                    clm_pf_idata%area_top_face_pfs)

    if ( (clm_pf_idata%nxclm_mapped >= 1 .and. clm_pf_idata%nyclm_mapped >= 1) .and. &
         (.not.option%mapping_files) ) then
      ! inactive cells with weighted top-surface area of 0 (i.e. CLM grid of inactive or zero coverage of land)
      ! by setting their 'material' id to -999
      call VecGetArrayReadF90(clm_pf_idata%area_top_face_pfs, topface_pf_loc, ierr)
      CHKERRQ(ierr)

      do ghosted_id = 1, grid%ngmax
        local_id = grid%nG2L(ghosted_id)
        if (ghosted_id <= 0 .or. local_id <= 0) cycle

        if (topface_pf_loc(ghosted_id)<=1.d-20 .and. associated(patch%imat)) then
          patch%imat(ghosted_id) = -999
        endif
      end do

      call VecRestoreArrayReadF90(clm_pf_idata%area_top_face_pfs, topface_pf_loc, ierr)
      CHKERRQ(ierr)
    end if

  end subroutine pflotranModelSetSoilDimension

! ************************************************************************** !

  subroutine pflotranModelSetSoilProp(pflotran_model)
  !
  ! Converts hydraulic properties from CLM units
  ! into PFLOTRAN units.
  !
  ! Author: Gautam Bisht
  ! Date: 10/22/2010
  !
  ! (fmy) 4/28/2014: modifications after updating to pflotran-dev

    use Option_module
    use Realization_Base_class
    use Simulation_Base_class
    use Discretization_module
    use Patch_module
    use Grid_module
    use Field_module
    use Material_module
    use Material_Aux_class
    use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, PERMEABILITY_XY, &
                               PERMEABILITY_YZ, PERMEABILITY_XZ, &
                               TORTUOSITY, POROSITY

    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use TH_Aux_module
    use Saturation_Function_module    ! this now still with 'TH_module'
    use Richards_Aux_module
    use Characteristic_Curves_module   ! this is used by Richards_module

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(discretization_type), pointer        :: discretization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field

    type(richards_auxvar_type), pointer       :: rich_auxvars(:)
    type(richards_auxvar_type), pointer       :: rich_auxvar
    class(characteristic_curves_type), pointer :: characteristic_curves

    type(th_auxvar_type), pointer             :: th_auxvars(:)
    type(th_auxvar_type), pointer             :: th_auxvar
    type(saturation_function_type), pointer :: saturation_function

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal          :: den, vis, grav
    PetscReal, pointer :: porosity_loc_p(:)
    PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)

    PetscScalar, pointer :: hksat_x_pf_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_y_pf_loc(:) ! hydraulic conductivity in y-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_z_pf_loc(:) ! hydraulic conductivity in z-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: watsat_pf_loc(:)  ! volumetric soil water at saturation (porosity)

    ! -------------------------------------------

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav =GRAVITY_CONSTANT       ! [m/S^2]


    subname = 'pflotranModelSetSoilProp'
!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select

    patch           => realization%patch
    discretization  => realization%discretization
    grid            => patch%grid
    field           => realization%field

    select case(option%iflowmode)
      case(RICHARDS_MODE)
        rich_auxvars   => patch%aux%Richards%auxvars
      case(TH_MODE)
        th_auxvars   => patch%aux%TH%auxvars
      case default
        !F.-M. Yuan: if no flowmode AND no reactive-transport, then let crash the run
        ! because 'uniform_velocity' with [0, 0, 0] velocity transport doesn't need specific flowmode,
        ! which is used for BGC coupling only (i.e., no TH coupling)
        if(option%ntrandof.le.0) then
            option%io_buffer =  &
               'Current PFLOTRAN mode not supported by pflotranModelSetSoilProp'
            call printErrMsg(option)
        endif
    end select

    ! ---------------------

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%hksat_x_clmp, &
                                    clm_pf_idata%hksat_x_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%hksat_y_clmp, &
                                    clm_pf_idata%hksat_y_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%hksat_z_clmp, &
                                    clm_pf_idata%hksat_z_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%watsat_clmp, &
                                    clm_pf_idata%watsat_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%bulkdensity_dry_clmp, &
                                    clm_pf_idata%bulkdensity_dry_pfs)

    !
    !---------------------
    !
    call VecGetArrayF90(clm_pf_idata%hksat_x_pfs, hksat_x_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_y_pfs, hksat_y_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_z_pfs, hksat_z_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_pfs,  watsat_pf_loc,  ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(field%porosity0, porosity_loc_p, ierr)
    CHKERRQ(ierr)

    if(option%iflowmode==RICHARDS_MODE .or. &
      option%iflowmode==TH_MODE) then
      ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
      call VecGetArrayF90(field%perm0_xx,  perm_xx_loc_p,  ierr)
      CHKERRQ(ierr)
      call VecGetArrayF90(field%perm0_yy,  perm_yy_loc_p,  ierr)
      CHKERRQ(ierr)
      call VecGetArrayF90(field%perm0_zz,  perm_zz_loc_p,  ierr)
      CHKERRQ(ierr)
    endif

    !
    do ghosted_id = 1, grid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (ghosted_id <= 0 .or. local_id <= 0) cycle
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) < 0) cycle    ! imat maybe 0, which causes issue
      endif

      select case(option%iflowmode)
        case(RICHARDS_MODE)
          rich_auxvar => rich_auxvars(ghosted_id)
        case(TH_MODE)
          th_auxvar => th_auxvars(ghosted_id)
      end select

      ! hydraulic conductivity => permissivity IS going to 'field%'
      ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
      ! [m^2]          [mm/sec]
      if(option%iflowmode==RICHARDS_MODE .or. &
         option%iflowmode==TH_MODE) then
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
      write(option%myrank+200,*) 'checking pflotran-model prior to set soil properties: ', &
        'rank=',option%myrank, 'ngmax=',grid%ngmax, 'nlmax=',grid%nlmax, &
        'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'pfp_porosity(local_id)=',porosity_loc_p(local_id), &
        'clms_watsat(ghosted_id)=',watsat_pf_loc(ghosted_id)
#endif

    enddo

    call VecRestoreArrayF90(clm_pf_idata%hksat_x_pfs, hksat_x_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y_pfs, hksat_y_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z_pfs, hksat_z_pf_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_pfs,  watsat_pf_loc,  ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayF90(field%porosity0, porosity_loc_p, ierr)
    CHKERRQ(ierr)
    if(option%iflowmode==RICHARDS_MODE .or. &
      option%iflowmode==TH_MODE) then
      ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
      call VecRestoreArrayF90(field%perm0_xx,  perm_xx_loc_p,  ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayF90(field%perm0_yy,  perm_yy_loc_p,  ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayF90(field%perm0_zz,  perm_zz_loc_p,  ierr)
      CHKERRQ(ierr)
    endif

    call MPI_Barrier(option%mycomm,ierr)

    ! update ghosted values after resetting soil physical properties from CLM
    call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                               field%work_loc,ONEDOF)

    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_MINERAL)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_CURRENT)

    if (option%nflowdof > 0) then
        call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_X,ZERO_INTEGER)
        call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Y,ZERO_INTEGER)
        call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Z,ZERO_INTEGER)
    endif

  end subroutine pflotranModelSetSoilProp

! ************************************************************************** !

  subroutine pflotranModelGetSoilPropFromPF(pflotran_model)
  !
  ! Pass soil physical properties from PFLOTRAN to CLM, if needed
  !
  ! Author: Fengming Yuan
  ! Date: 1/30/2014
  !

    use Realization_Base_class
    use Patch_module
    use Grid_module
    use Field_module
    use Option_module

    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use Richards_Aux_module
    use Characteristic_Curves_module
    use TH_Aux_module
    use Saturation_Function_module

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    type(saturation_function_type), pointer :: saturation_function
    type(th_auxvar_type), pointer :: th_auxvar

    class(characteristic_curves_type), pointer :: characteristic_curves
    type(richards_auxvar_type), pointer :: rich_auxvar

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal          :: tempreal, tempreal2

    ! pf internal variables
    PetscReal, pointer :: porosity_loc_p(:)

    ! clm-pf-interface Vecs for PF thermal-hydroloical parameters used in CLM-PFLOTRAN interface
    PetscScalar, pointer :: porosity_loc_pfp(:)  ! soil porosity
    PetscScalar, pointer :: sr_pcwmax_loc_pfp(:) ! soil vwc at max. capillary pressure (note: not 'Sr')
    PetscScalar, pointer :: pcwmax_loc_pfp(:)    ! max. capillary pressure

    character(len=MAXSTRINGLENGTH) :: error_string
    PetscInt :: cur_sat_func_id

    subname = 'pflotranModelGetSoilPropFromPF'

!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
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
      if (ghosted_id <= 0 .or. local_id <= 0) cycle
      if (associated(patch%imat)) then
         if (patch%imat(ghosted_id) < 0) cycle    ! imat maybe 0, which causes issue
      endif

      ! PF's porosity
      porosity_loc_pfp(local_id) = porosity_loc_p(local_id)

      ! soil hydraulic properties ID for current cell
      cur_sat_func_id = patch%sat_func_id(ghosted_id)

      !
      if (option%iflowmode==RICHARDS_MODE) then
        ! Richards_MODE now is using 'charateristic_curves' module

        characteristic_curves => patch% &
          characteristic_curves_array(cur_sat_func_id)%ptr
        rich_auxvar => patch%aux%Richards%auxvars(ghosted_id)

        select type(sf => characteristic_curves%saturation_function)
          !class is(sat_func_VG_type)
             ! not-yet (TODO)
          class is(sat_func_BC_type)
            ! PF's limits on soil matrix potential (Capillary pressure)
            pcwmax_loc_pfp(local_id) = sf%pcmax

            ! PF's limits on soil water at pcwmax (NOT: not 'Sr', at which PC is nearly 'inf')
            call sf%Saturation(sf%pcmax, tempreal, tempreal2, option)
            sr_pcwmax_loc_pfp(local_id) = tempreal

          class default
            option%io_buffer = 'Currently ONLY support Brooks_COREY saturation function type' // &
              ' when coupled with CLM.'
            call printErrMsg(option)
        end select

      !
      elseif (option%iflowmode==TH_MODE) then
        ! 'saturation_function' module NOW only with TH_mode
        saturation_function => patch%saturation_function_array(cur_sat_func_id)%ptr
        th_auxvar => patch%aux%TH%auxvars(ghosted_id)

        ! PF's limits on soil matrix potential (Capillary pressure)
        pcwmax_loc_pfp(local_id) = saturation_function%pcwmax

        ! PF's limits on soil water at pcwmax (NOT: not 'Sr', at which PC is nearly 'inf')
        call SaturationFunctionCompute(saturation_function%pcwmax, tempreal, &
                                      saturation_function, &
                                      option)
        sr_pcwmax_loc_pfp(local_id) = tempreal
      endif

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
                                    option, &
                                    clm_pf_idata%effporosity_pfp, &
                                    clm_pf_idata%effporosity_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%sr_pcwmax_pfp, &
                                    clm_pf_idata%sr_pcwmax_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%pcwmax_pfp, &
                                    clm_pf_idata%pcwmax_clms)

    ! reference pressure
    clm_pf_idata%pressure_reference = option%reference_pressure

  end subroutine pflotranModelGetSoilPropFromPF
  ! ************************************************************************** !


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  
end module pflotran_clm_main_module

