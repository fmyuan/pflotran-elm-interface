module pflotran_model_module

  use Simulation_module
  use Realization_class
  use Timestepper_module
  use Option_module
  use Input_module
  use Init_module
  use Logging_module
  use Stochastic_module
  use Stochastic_Aux_module
  use Waypoint_module
  use Units_module

#if defined (CLM_PFLOTRAN)
  ! some CLM internals that we want to use (e.g. i/o)
  use clm_varctl, only            : iulog
  use spmdMod, only : masterproc
#endif  

#if defined (CLM_PFLOTRAN)
  use Mapping_module
#endif  

#if defined (CLM_PFLOTRAN) || defined(CLM_OFFLINE)  
  use Richards_Aux_module
  !use pflotran_clm_interface_type
  !use clm_pflotran_interface_data
#endif

  use Mapping_module

  implicit none

#include "definitions.h"
!  include "piof.h"
#include "finclude/petsclog.h"
#include "finclude/petscsysdef.h"
#include "finclude/petscviewer.h"
#include "finclude/petscvec.h"
!#include "finclude/petscvec.h90"

#ifndef CLM_PFLOTRAN
  ! don't have CLM, so we need to mock a few things
  integer, public :: iulog = 6       ! "stdout" log file unit number, default is 6
  logical, public :: masterproc = .false. ! master io processor
#endif

  PetscLogDouble :: timex(4), timex_wall(4)

  PetscBool :: truth
  PetscBool :: option_found
  PetscBool :: single_inputfile
  PetscInt  :: i
  PetscInt  :: temp_int

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)

  !#ifdef CLM_PFLOTRAN
  type, public :: inside_each_overlapped_cell
     PetscInt           :: id
     PetscInt           :: ocell_count
     PetscInt,  pointer :: ocell_id(:)
     PetscReal, pointer :: perc_vol_overlap(:)
     PetscReal          :: total_vol_overlap
  end type inside_each_overlapped_cell

  !#endif

  type, public :: pflotran_model_type
    type(stochastic_type),  pointer :: stochastic
    type(simulation_type),  pointer :: simulation
    type(realization_type), pointer :: realization
    type(option_type),      pointer :: option
#ifdef CLM_PFLOTRAN
    !type(mapping_type),     pointer :: mapping
#endif
    PetscReal :: pause_time_1
    PetscReal :: pause_time_2
    type(inside_each_overlapped_cell), pointer :: pf_cells(:)
    type(inside_each_overlapped_cell), pointer :: clm_cells(:)
    type(mapping_type),                pointer :: map_clm2pf
    type(mapping_type),                pointer :: map_clm2pf_soils
    type(mapping_type),                pointer :: map_pf2clm
    !type(clm_pflotran_interface_data_type),pointer            :: clm_pf_idata
     
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

  end type pflotran_model_type

  public::pflotranModelCreate,               &
       pflotranModelInitMapping3,            &
       pflotranModelSetSoilProp3,            & !
       pflotranModelSetICs3,                 &
       pflotranModelUpdateSourceSink3,       & !
       pflotranModelGetSaturation3,          & !
       pflotranModelStepperRunInit,          &
       pflotranModelStepperRunTillPauseTime, &
       pflotranModelUpdateTopBCHomogeneous,  &
       pflotranModelStepperRunFinalize,      &
       pflotranModelInsertWaypoint,          &
       pflotranModelDeleteWaypoint,          &
       pflotranModelDestroy

contains

  ! ************************************************************************** !
  !
  ! pflotranModelCreate: Allocates and initializes the pflotranModel object.
  !             It performs the same sequence of commands as done in pflotran.F90
  !     before model integration is performed by the call to StepperRun()
  !             routine
  !
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  function pflotranModelCreate(mpicomm)

    use Simulation_module
    use Realization_class
    use Timestepper_module
    use Option_module
    use Input_module
    use Init_module
    use Logging_module
    use Stochastic_module
    use Stochastic_Aux_module
    use String_module
#ifdef CLM_PFLOTRAN
    !use pflotran_clm_interface_type
#endif
    implicit none

    PetscInt, intent(in) :: mpicomm

    type(pflotran_model_type), pointer :: pflotranModelCreate

    PetscLogDouble :: timex(4), timex_wall(4)

    PetscBool :: truth
    PetscBool :: option_found
    PetscBool :: input_prefix_option_found
    PetscBool :: pflotranin_option_found
    PetscBool :: single_inputfile
    PetscBool :: clm2pf_soil_file
    PetscBool :: clm2pf_flux_file
    PetscBool :: pf2clm_flux_file
    
    PetscInt   :: temp_int

    PetscErrorCode :: ierr
    character(len=MAXSTRINGLENGTH)          :: string
    character(len=MAXSTRINGLENGTH), pointer :: filenames(:)
    character(len=MAXWORDLENGTH) :: word

    type(pflotran_model_type),      pointer :: pflotran_model


    allocate(pflotran_model)
    allocate(pflotran_model%stochastic)
    allocate(pflotran_model%simulation)
    allocate(pflotran_model%realization)
    allocate(pflotran_model%option)
#ifdef CLM_PFLOTRAN
    !allocate(pflotran_model%mapping)
#endif
    allocate(pflotran_model%map_clm2pf)
    allocate(pflotran_model%map_clm2pf_soils)
    allocate(pflotran_model%map_pf2clm)

    ! TODO(bja): creating memory leaks?
    nullify(pflotran_model%stochastic)
    nullify(pflotran_model%simulation)
    nullify(pflotran_model%realization)
    nullify(pflotran_model%option)
#ifdef CLM_PFLOTRAN
    !nullify(pflotran_model%mapping)
#endif
    nullify(pflotran_model%pf_cells)
    nullify(pflotran_model%clm_cells)
    nullify(pflotran_model%map_clm2pf)
    nullify(pflotran_model%map_clm2pf_soils)
    nullify(pflotran_model%map_pf2clm)

    pflotran_model%option => OptionCreate()
    pflotran_model%option%fid_out = 16
    single_inputfile = PETSC_TRUE

    pflotran_model%pause_time_1 = -1.0d0
    pflotran_model%pause_time_2 = -1.0d0


    pflotran_model%option%global_comm = mpicomm

    call MPI_Comm_rank(MPI_COMM_WORLD,pflotran_model%option%global_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,pflotran_model%option%global_commsize,ierr)
    call MPI_Comm_group(MPI_COMM_WORLD,pflotran_model%option%global_group,ierr)
    pflotran_model%option%mycomm = pflotran_model%option%global_comm
    pflotran_model%option%myrank = pflotran_model%option%global_rank
    pflotran_model%option%mycommsize = pflotran_model%option%global_commsize
    pflotran_model%option%mygroup = pflotran_model%option%global_group

#ifndef CLM_PFLOTRAN
    ! mock the master i/o processor flag
    if (pflotran_model%option%myrank == pflotran_model%option%io_rank) then
       masterproc = .true.
    endif
#endif


    ! check for non-default input filename
    pflotran_model%option%input_filename = "pflotran.in"
    string = '-pflotranin'
    call InputGetCommandLineString(string,pflotran_model%option%input_filename, &
                                  pflotranin_option_found,pflotran_model%option)

    string = '-input_prefix'
    call InputGetCommandLineString(string,pflotran_model%option%input_prefix, &
                                  input_prefix_option_found,pflotran_model%option)
  
    if (pflotranin_option_found .and. input_prefix_option_found) then
      pflotran_model%option%io_buffer = 'Cannot specify both "-pflotranin" and ' // &
        '"-input_prefix" on the command lines.'
      call printErrMsg(pflotran_model%option)
    else if (pflotranin_option_found) then
      !TODO(geh): replace this with StringSplit()
      i = index(pflotran_model%option%input_filename,'.',PETSC_TRUE)
      if (i > 1) then
        i = i-1
      else
        ! for some reason len_trim doesn't work on MS Visual Studio in
        ! this location
        i = len(trim(pflotran_model%option%input_filename))
      endif
      pflotran_model%option%input_prefix = pflotran_model%option%input_filename(1:i)
    else if (input_prefix_option_found) then
      pflotran_model%option%input_filename = trim(pflotran_model%option%input_prefix) // '.in'
    endif

    string = '-output_prefix'
    call InputGetCommandLineString(string,pflotran_model%option%global_prefix,option_found,pflotran_model%option)
    if (.not.option_found) pflotran_model%option%global_prefix = pflotran_model%option%input_prefix

    string = '-screen_output'
    call InputGetCommandLineTruth(string,pflotran_model%option%print_to_screen,option_found,pflotran_model%option)

    string = '-file_output'
    call InputGetCommandLineTruth(string,pflotran_model%option%print_to_file,option_found,pflotran_model%option)

    string = '-v'
    call InputGetCommandLineTruth(string,truth,option_found,pflotran_model%option)
    if (option_found) pflotran_model%option%verbosity = 1

    string = '-multisimulation'
    call InputGetCommandLineTruth(string,truth,option_found,pflotran_model%option)
    if (option_found) then
       single_inputfile = PETSC_FALSE
    endif

    string = '-stochastic'
    call InputGetCommandLineTruth(string,truth,option_found,pflotran_model%option)
    if (option_found) pflotran_model%stochastic => StochasticCreate()

    call InitReadStochasticCardFromInput(pflotran_model%stochastic,pflotran_model%option)

    if (associated(pflotran_model%stochastic)) then
       !  call StochasticInit(stochastic,option)
       !  call StochasticRun(stochastic,option)
    endif

    if (single_inputfile) then
       if (masterproc) then
          write(iulog, *),'single_inputfile'
       endif

       PETSC_COMM_WORLD = MPI_COMM_WORLD
       call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    else
       call InitReadInputFilenames(pflotran_model%option,filenames)
       temp_int = size(filenames)
       call SimulationCreateProcessorGroups(pflotran_model%option,temp_int)
       pflotran_model%option%input_filename = filenames(pflotran_model%option%mygroup_id)
       i = index(pflotran_model%option%input_filename,'.',PETSC_TRUE)
       if (i > 1) then
          i = i-1
       else
          ! for some reason len_trim doesn't work on MS Visual Studio in
          ! this location
          i = len(trim(pflotran_model%option%input_filename))
       endif
       pflotran_model%option%global_prefix = pflotran_model%option%input_filename(1:i)
       write(string,*) pflotran_model%option%mygroup_id
       pflotran_model%option%group_prefix = 'G' // trim(adjustl(string))
    endif

    ! TODO(bja): debuging output, or do we want this formally written to pflotran/clm output?
    write(*, *),'option%global_prefix = ',pflotran_model%option%global_prefix
    write(*, *),'option%group_prefix = ',pflotran_model%option%group_prefix

    if (pflotran_model%option%verbosity > 0) then
       call PetscLogBegin(ierr)
       string = '-log_summary'
       call PetscOptionsInsertString(string, ierr)
    endif
    call LoggingCreate()

    pflotran_model%simulation => SimulationCreate(pflotran_model%option)
    pflotran_model%realization => pflotran_model%simulation%realization

    call OptionCheckCommandLine(pflotran_model%option)

    call PetscGetCPUTime(timex(1), ierr)
    call PetscTime(timex_wall(1), ierr)
    pflotran_model%option%start_time = timex_wall(1)

    call Init(pflotran_model%simulation)

    pflotran_model%map_clm2pf       => MappingCreate()
    pflotran_model%map_clm2pf_soils => MappingCreate()
    pflotran_model%map_pf2clm       => MappingCreate()

    pflotran_model%nlclm = -1
    pflotran_model%ngclm = -1

    pflotran_model%realization%input => InputCreate(15,pflotran_model%option%input_filename,pflotran_model%option)

    ! Read names of mapping file
    clm2pf_soil_file=PETSC_FALSE
    clm2pf_flux_file=PETSC_FALSE
    pf2clm_flux_file=PETSC_FALSE
    do
      call InputReadFlotranString(pflotran_model%realization%input, pflotran_model%option)
      if (pflotran_model%realization%input%ierr /= 0) exit

      call InputReadWord(pflotran_model%realization%input, pflotran_model%option, word, PETSC_TRUE)
      call InputErrorMsg(pflotran_model%realization%input, pflotran_model%option, 'keyword', 'MAPPING_FILES')
      call StringToUpper(word)

      select case(trim(word))
        case('CLM2PF_FLUX_FILE')
          call InputReadWord(pflotran_model%realization%input, &
                             pflotran_model%option, &
                             pflotran_model%map_clm2pf%filename, PETSC_TRUE)
          pflotran_model%map_clm2pf%filename = trim(pflotran_model%map_clm2pf%filename)//CHAR(0)
          call InputErrorMsg(pflotran_model%realization%input, &
                             pflotran_model%option,'type', 'MAPPING_FILES')   
          call StringToLower(pflotran_model%map_clm2pf%filename)
          clm2pf_flux_file=PETSC_TRUE
        case('CLM2PF_SOIL_FILE')
          call InputReadWord(pflotran_model%realization%input, &
                             pflotran_model%option, &
                             pflotran_model%map_clm2pf_soils%filename, PETSC_TRUE)
          call InputErrorMsg(pflotran_model%realization%input, &
                             pflotran_model%option,'type', 'MAPPING_FILES')   
          call StringToLower(pflotran_model%map_clm2pf_soils%filename)
          clm2pf_soil_file=PETSC_TRUE
        case('PF2CLM_FLUX_FILE')
          call InputReadWord(pflotran_model%realization%input, &
                             pflotran_model%option, &
                             pflotran_model%map_pf2clm%filename, PETSC_TRUE)
          call InputErrorMsg(pflotran_model%realization%input, &
                             pflotran_model%option,'type', 'MAPPING_FILES')   
          call StringToLower(pflotran_model%map_pf2clm%filename)
          pf2clm_flux_file=PETSC_TRUE
      end select

    enddo
    call InputDestroy(pflotran_model%realization%input)

    if((.not.clm2pf_soil_file).or.(.not.clm2pf_flux_file).or.(.not.pf2clm_flux_file)) then
      pflotran_model%option%io_buffer='One of the mapping files not found'
      call printErrMsg(pflotran_model%option)
    endif

    pflotranModelCreate => pflotran_model

  end function pflotranModelCreate


  ! ************************************************************************** !
  !
  ! pflotranModelStepperRunInit: It performs the same execution of commands
  !             that are carried out in StepperRun() before the model integration
  !             begins over the entire simulation time interval
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelStepperRunInit(pflotran_model)

    type(pflotran_model_type), pointer :: pflotran_model

#ifndef SURFACE_FLOW
    call StepperRunInit(pflotran_model%simulation%realization, &
         pflotran_model%simulation%flow_stepper, &
         pflotran_model%simulation%tran_stepper)
#else
    call StepperRunInit(pflotran_model%simulation%realization, &
         pflotran_model%simulation%surf_realization, &
         pflotran_model%simulation%flow_stepper, &
         pflotran_model%simulation%tran_stepper)
#endif

  end subroutine pflotranModelStepperRunInit


  ! ************************************************************************** !
  !
  ! pflotranModelSetICs3:
  !
  !
  ! author: Gautam Bisht
  ! date: 10/22/2010
  ! ************************************************************************** !
subroutine pflotranModelSetICs3(pflotran_model)

    use Realization_class
    use Patch_module
    use Grid_module
    use Richards_Aux_module
    use Field_module
    use clm_pflotran_interface_data
    use Global_Aux_module
    use Discretization_module
    use Richards_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    type(realization_type),pointer            :: realization
    type(patch_type),pointer                  :: patch
    type(grid_type),pointer                   :: grid
    type(field_type),pointer                  :: field
    type(richards_auxvar_type), pointer       :: aux_var
    type(global_auxvar_type), pointer         :: global_aux_vars(:)


    PetscErrorCode     :: ierr
    PetscInt           :: local_id,ghosted_id
    PetscReal          :: den, vis, grav
    PetscReal,pointer  :: xx_loc_p(:)

    PetscScalar, pointer :: press_pf_loc(:) ! Pressure [Pa]

    realization     => pflotran_model%simulation%realization
    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field
    global_aux_vars  => patch%aux%Global%aux_vars

    call MappingSourceToDestination(pflotran_model%map_clm2pf, &
                                    pflotran_model%option, &
                                    clm_pf_idata%press_clm, &
                                    clm_pf_idata%press_pf)

    call GridVecGetArrayF90(grid,field%flow_xx,xx_loc_p,ierr)
    call VecGetArrayF90(clm_pf_idata%press_pf,press_pf_loc,ierr)

    do local_id = 1, grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
       xx_loc_p(ghosted_id) = press_pf_loc(local_id)
    enddo

    call GridVecRestoreArrayF90(grid,field%flow_xx,xx_loc_p,ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_pf,press_pf_loc,ierr)

    ! update dependent vectors: Saturation
    call DiscretizationGlobalToLocal(realization%discretization,field%flow_xx, &
         field%flow_xx_loc,NFLOWDOF)
    call VecCopy(field%flow_xx, field%flow_yy, ierr)

    call RichardsUpdateAuxVars(realization)

end subroutine pflotranModelSetICs3

  ! ************************************************************************** !
  !
  ! pflotranModelSetSoilProp3:
  !
  !
  ! author: Gautam Bisht
  ! date: 10/22/2010
  ! ************************************************************************** !
!#ifdef CLM_PFLOTRAN
  subroutine pflotranModelSetSoilProp3(pflotran_model)

    use Realization_class
    use Patch_module
    use Grid_module
    use Richards_Aux_module
    use Field_module
    use clm_pflotran_interface_data

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    type(realization_type),pointer            :: realization
    type(patch_type),pointer                  :: patch
    type(grid_type),pointer                   :: grid
    type(field_type),pointer                  :: field
    type(richards_auxvar_type), pointer       :: rich_aux_vars(:)
    type(richards_auxvar_type), pointer       :: aux_var

    PetscErrorCode     :: ierr
    PetscInt           :: local_id
    PetscReal          :: den, vis, grav
    PetscReal,pointer  :: porosity_loc_p(:), vol_ovlap_arr(:)
    PetscReal,pointer  :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)

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

    realization     => pflotran_model%simulation%realization
    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field
    rich_aux_vars   => patch%aux%Richards%aux_vars

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

    call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
    call GridVecGetArrayF90(grid,field%perm_xx_loc,  perm_xx_loc_p,  ierr)
    call GridVecGetArrayF90(grid,field%perm_yy_loc,  perm_yy_loc_p,  ierr)
    call GridVecGetArrayF90(grid,field%perm_zz_loc,  perm_zz_loc_p,  ierr)

    do local_id = 1,grid%ngmax

      aux_var => rich_aux_vars(local_id)

      ! bc_alpha [1/Pa]; while sucsat [mm of H20]
      ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
      aux_var%bc_alpha = 1.d0/(sucsat_pf_loc(local_id)*grav) 

      ! bc_lambda = 1/bsw
      aux_var%bc_lambda = 1.d0/bsw_pf_loc(local_id)
      
      ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
      ! [m^2]          [mm/sec]
      perm_xx_loc_p(local_id) = hksat_x_pf_loc(local_id) * vis/ (den * grav) / 1000.d0
      perm_yy_loc_p(local_id) = hksat_y_pf_loc(local_id) * vis/ (den * grav) / 1000.d0
      perm_zz_loc_p(local_id) = hksat_z_pf_loc(local_id) * vis/ (den * grav) / 1000.d0

      porosity_loc_p(local_id) = watsat_pf_loc(local_id)

    enddo

    call VecRestoreArrayF90(clm_pf_idata%hksat_x_pf, hksat_x_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y_pf, hksat_y_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z_pf, hksat_z_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat_pf,  sucsat_pf_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_pf,  watsat_pf_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw_pf,     bsw_pf_loc,     ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw_clm,    bsw_clm_loc,    ierr)

    call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
    call GridVecRestoreArrayF90(grid,field%perm_xx_loc,  perm_xx_loc_p,  ierr)
    call GridVecRestoreArrayF90(grid,field%perm_yy_loc,  perm_yy_loc_p,  ierr)
    call GridVecRestoreArrayF90(grid,field%perm_zz_loc,  perm_zz_loc_p,  ierr)

  end subroutine pflotranModelSetSoilProp3
!#endif

  ! ************************************************************************** !
  !
  ! pflotranModelInitMapping3:
  !
  !
  ! author: Gautam Bisht
  ! date: 03/24/2011
  ! ************************************************************************** !
  subroutine pflotranModelInitMapping3(pflotran_model,  &
                                       grid_clm_cell_ids_ghosted_nindex, &
                                       grid_clm_npts_local, &
                                       map_id)

    use Input_module
    use Option_module
    use Realization_class
    use Grid_module
    use Patch_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in),pointer                      :: grid_clm_cell_ids_ghosted_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename
    
    ! local
    PetscInt                           :: local_id, grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_clm_npts_ghost, source_id
    PetscInt, pointer                  :: grid_pf_cell_ids_ghosted_nindex(:)
    PetscInt, pointer                  :: grid_pf_cell_ids_local_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_or_ghost_nindex(:)
    PetscInt, pointer                  :: grid_clm_local_or_ghost_nindex(:)
    PetscInt                           :: count
    PetscErrorCode                     :: ierr
    type(mapping_type),pointer         :: map

    type(option_type), pointer         :: option
    type(realization_type), pointer    :: realization
    type(grid_type),pointer            :: grid
    type(patch_type), pointer          :: patch

    option          => pflotran_model%option
    realization     => pflotran_model%simulation%realization
    patch           => realization%patch
    patch           => realization%patch
    grid            => patch%grid

    select case(map_id)
      case(CLM2PF_FLUX_MAP_ID)
        map => pflotran_model%map_clm2pf
        source_id = CLM_MESH
      case(CLM2PF_SOIL_MAP_ID)
        map => pflotran_model%map_clm2pf_soils
        source_id = CLM_MESH
      case(PF2CLM_FLUX_MAP_ID)
        map => pflotran_model%map_pf2clm
        source_id = PF_MESH
      case default
        option%io_buffer = 'Invalid map_id argument to pflotranModelInitMapping3'
        call printErrMsg(option)
    end select

    call MPI_Barrier(option%mycomm, ierr)
    grid_pf_npts_local = grid%nlmax
    grid_pf_npts_ghost = grid%ngmax - grid%nlmax
    grid_clm_npts_ghost= 0

    allocate(grid_clm_local_or_ghost_nindex (grid_clm_npts_local))
    do local_id = 1, grid_clm_npts_local
      grid_clm_local_or_ghost_nindex(local_id) = 1 ! LOCAL
    enddo

    allocate(grid_pf_cell_ids_ghosted_nindex(grid%ngmax))
    do local_id = 1,grid%ngmax
      grid_pf_cell_ids_ghosted_nindex(local_id) = grid%nG2A(local_id)-1
    enddo

    select case(source_id)
      case(CLM_MESH)

        allocate(grid_pf_local_or_ghost_nindex  (grid%ngmax))

        do local_id = 1,grid%ngmax
          if (grid%nG2L(local_id) == 0) then
            grid_pf_local_or_ghost_nindex(local_id) = 0 ! GHOST
          else
            grid_pf_local_or_ghost_nindex(local_id) = 1 ! LOCAL
          endif
        enddo

        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_cell_ids_ghosted_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_ghosted_nindex, &
                                              grid_pf_local_or_ghost_nindex)
        deallocate(grid_pf_local_or_ghost_nindex)
      case(PF_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_pf_npts_local, &
                                        grid_pf_cell_ids_ghosted_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_clm_npts_local, &
                                              grid_clm_npts_ghost, &
                                              grid_clm_cell_ids_ghosted_nindex, &
                                              grid_clm_local_or_ghost_nindex)
      case default
        option%io_buffer = 'Invalid argument source_id passed to pflotranModelInitMapping3'
        call printErrMsg(option)
    end select
    deallocate(grid_pf_cell_ids_ghosted_nindex)
    
    call MappingReadTxtFile(map, map%filename, option)
    call MappingDecompose(map,option)
    call MappingFindDistinctSourceMeshCellIds(map,option)
    call MappingCreateWeightMatrix(map,option)
    call MappingCreateScatterOfSourceMesh(map, option) 

end subroutine pflotranModelInitMapping3

  ! ************************************************************************** !
  !
  ! pflotranModelFindOverlapCells:
  !
  !
  ! author: Gautam Bisht
  ! date: 01/28/2011
  ! ************************************************************************** !
  subroutine pflotranModelFindOverlapCells(grid_cells, grid_npts_ghosted, &
       grid_cell_ids_ghosted_nindex, ocell_ids, ocell_vol, ocell_count, grid_ocell_count)

    implicit none

    PetscInt                           :: ii,jj,tmp_idx
    PetscInt                           :: cell_id
    PetscInt                           :: grid_npts_ghosted
    PetscInt,pointer                   :: grid_cell_ids_ghosted_nindex(:)
    PetscInt,intent(out),pointer       :: ocell_ids(:)
    PetscReal,intent(out),pointer      :: ocell_vol(:)
    PetscInt,intent(out)               :: ocell_count
    PetscInt,intent(out),pointer       :: grid_ocell_count(:)

    type(inside_each_overlapped_cell),pointer :: grid_cells(:)

    ocell_count = 0
    do ii = 1,grid_npts_ghosted
      cell_id              = grid_cell_ids_ghosted_nindex(ii) + 1
      ocell_count          = ocell_count + grid_cells(cell_id)%ocell_count
      grid_ocell_count(ii) = grid_cells(cell_id)%ocell_count
    enddo
    
    allocate(ocell_ids(ocell_count))
    allocate(ocell_vol(ocell_count))

    tmp_idx = 1
    do ii = 1,grid_npts_ghosted
      cell_id = grid_cell_ids_ghosted_nindex(ii) + 1
      do jj = 1,grid_cells(cell_id)%ocell_count
        ocell_ids(tmp_idx) = grid_cells(cell_id)%ocell_id(jj)
        ocell_vol(tmp_idx) = grid_cells(cell_id)%perc_vol_overlap(jj)
        tmp_idx = tmp_idx + 1
      enddo
    enddo

  end subroutine pflotranModelFindOverlapCells

  ! ************************************************************************** !
  !
  ! pflotranModelReadMappingFile:
  !
  !
  ! author: Gautam Bisht
  ! date: 01/28/2011
  ! ************************************************************************** !
  subroutine pflotranModelReadMappingFile(option, filename, ocells, num_cells)

    use Input_module
    use Option_module

    implicit none

    PetscInt  :: cell_id, ocell_id
    PetscInt  :: num_cells, num_ocells
    PetscReal :: ocell_vol
    PetscInt  :: fileid
    PetscInt  :: ii, jj

    character(len=MAXSTRINGLENGTH)     :: filename
    character(len=MAXWORDLENGTH)       :: card
    type(input_type), pointer          :: input
    type(option_type), pointer         :: option
    type(inside_each_overlapped_cell),pointer :: ocells(:)

    fileid   = 20
    card     = 'pflotran_interface_main'

    input => InputCreate(fileid,filename,option)
    call InputReadFlotranString(input,option)
    call InputErrorMsg(input,option,'number of cells',card)

    num_cells = -1
    call InputReadInt(input,option,num_cells)

    allocate(ocells(num_cells))

    do ii = 1,num_cells

       ! Read cell-id and #cells overlapped with
       call InputReadFlotranString(input,option)
       call InputReadInt(input,option,cell_id)
       call InputReadInt(input,option,num_ocells)
       ocells(ii)%id          = cell_id
       ocells(ii)%ocell_count = num_ocells

       ! Nullify and allocate memory
       nullify(ocells(ii)%ocell_id)
       nullify(ocells(ii)%perc_vol_overlap)

       if ( num_ocells.gt.0) then
          allocate(ocells(ii)%ocell_id(        num_ocells) )
          allocate(ocells(ii)%perc_vol_overlap(num_ocells) )
       endif

       ! Initialize
       ocells(ii)%total_vol_overlap = 0.0d0

       ! Read overlapped cell ids and volume overlapped
       do jj = 1,num_ocells
          call InputReadInt(input,option,ocell_id)
          call InputReadDouble(input,option,ocell_vol)

          ocells(ii)%ocell_id(jj)         = ocell_id
          ocells(ii)%perc_vol_overlap(jj) = ocell_vol
          ocells(ii)%total_vol_overlap  = &
               ocells(ii)%total_vol_overlap + ocell_vol
       enddo
    enddo

  end subroutine pflotranModelReadMappingFile

  !#endif

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
  subroutine pflotranModelStepperRunTillPauseTime(pflotran_model, pause_time)

    use Richards_module
    use TH_module
    use Grid_module
    use Patch_module
    use Field_module
    use Global_Aux_module

    implicit none
#include "definitions.h"

    type(pflotran_model_type), pointer :: pflotran_model
    type(realization_type), pointer    :: realization
    type(grid_type),pointer            :: grid
    type(field_type),pointer           :: field
    type(global_auxvar_type), pointer  :: global_aux_vars(:)
    type(patch_type), pointer          :: patch

    PetscErrorCode :: ierr
    PetscReal  :: pause_time
    PetscReal  :: dtime
    PetscReal  :: liq_vol_start ! [m^3/m^3]
    PetscReal  :: liq_vol_end   ! [m^3/m^3]
    PetscReal  :: del_liq_vol, dz, sat, source_sink
    PetscInt   :: local_id, ghosted_id
    PetscReal,pointer  :: porosity_loc_p(:)

    realization     => pflotran_model%simulation%realization
    patch           => realization%patch
    grid            => patch%grid
    global_aux_vars => patch%aux%Global%aux_vars
    field           => realization%field

    if (masterproc) then
       write(iulog,*), '>>>> Inserting waypoint at pause_time = ',pause_time
    endif

    call pflotranModelInsertWaypoint(pflotran_model, pause_time)
    call pflotranModelInsertWaypoint(pflotran_model, pause_time + 100.0d0)

    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        call RichardsUpdateAuxVars(pflotran_model%simulation%realization)
      case (TH_MODE)
        call THUpdateAuxVars(pflotran_model%simulation%realization)
      case default
        pflotran_model%option%io_buffer='pflotranModelStepperRunTillPauseTime ' // &
          'not implmented for this mode.'
        call printErrMsg(pflotran_model%option)
    end select

    call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
#if 1
    liq_vol_start = 0.d0
    liq_vol_end   = 0.d0
    source_sink   = 0.d0

    ! Volume of Water before PFLOTRAN is called
    do local_id= 1,grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
       !sat           = global_aux_vars(ghosted_id)%sat(1)
       !dz            = grid%structured_grid%dz(ghosted_id)
       !del_liq_vol   = sat * porosity_loc_p(ghosted_id) * dz
       !liq_vol_start = liq_vol_start + del_liq_vol
    enddo
#endif

#ifndef SURFACE_FLOW
    call StepperRunOneDT(pflotran_model%simulation%realization, &
         pflotran_model%simulation%flow_stepper, &
         pflotran_model%simulation%tran_stepper, pause_time)
#else
    call StepperRunOneDT(pflotran_model%simulation%realization, &
         pflotran_model%simulation%surf_realization, &
         pflotran_model%simulation%flow_stepper, &
         pflotran_model%simulation%surf_flow_stepper, &
         pflotran_model%simulation%tran_stepper, pause_time)
#endif

#ifdef CLM_PFLOTRAN
    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        call RichardsUpdateAuxVars(pflotran_model%simulation%realization)
      case (TH_MODE)
        call THUpdateAuxVars(pflotran_model%simulation%realization)
      case default
        pflotran_model%option%io_buffer='pflotranModelStepperRunTillPauseTime ' // &
          'not implmented for this mode.'
        call printErrMsg(pflotran_model%option)
    end select

    ! Volume of Water after PFLOTRAN finished
    do local_id= 1,grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
       sat           = global_aux_vars(ghosted_id)%sat(1)
      !dz            = grid%structured_grid%dz(ghosted_id)
       del_liq_vol   = sat * porosity_loc_p(ghosted_id) * dz
       liq_vol_end   = liq_vol_end + del_liq_vol

       !source_sink   = source_sink + pf_clm_data%qsrc_flx(local_id)
    enddo

    !write(iulog, *), '===================================================='
    !write(iulog, *), 'Volume of LIQUID water:'
    !write(iulog, *), 'Before PFLOTRAN call: ',liq_vol_start,                        '[m^3/area]'
    !write(iulog, *), 'After  PFLOTRAN call: ',liq_vol_end  ,                        '[m^3/area]'
    !write(iulog, *), 'Change              : ',liq_vol_start-liq_vol_end,            '[m^3/area]'
    !write(iulog, *), 'Rate of change      : ',(liq_vol_end-liq_vol_start)/1800.0d0, '[m^3/area/sec]'
    !write(iulog, *), 'Source_sink         : ',source_sink/998.2d0,                  '[m/sec]'
    !write(iulog, *), '====================================================?'

    call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
#endif


    if(pflotran_model%pause_time_1.gt.0.0d0) then
       call pflotranModelDeleteWaypoint(pflotran_model, pflotran_model%pause_time_1)
    endif

    if(pflotran_model%pause_time_2.gt.0.0d0) then
       call pflotranModelDeleteWaypoint(pflotran_model, pflotran_model%pause_time_2)
    endif

    pflotran_model%pause_time_1 = pause_time
    pflotran_model%pause_time_2 = pause_time + 100.0d0

  end subroutine pflotranModelStepperRunTillPauseTime

  ! ************************************************************************** !
  !
  !
  ! author: Gautam Bisht
  ! date: 11/22/2011
  ! ************************************************************************** !
  subroutine pflotranModelUpdateSourceSink3 (pflotran_model)

    use clm_pflotran_interface_data

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    !type(clm_pflotran_interface_data_type),pointer :: clm_pf_idata
#if 0
    call MappingSourceToDestination(pflotran_model%map_clm2pf, &
                                    pflotran_model%option, &
                                    clm_pf_idata%qflx_clm, &
                                    clm_pf_idata%qflx_pf)
#endif
  end subroutine pflotranModelUpdateSourceSink3

  ! ************************************************************************** !
  !
  !
  ! author: Gautam Bisht
  ! date: 11/22/2011
  ! ************************************************************************** !
  subroutine pflotranModelGetSaturation3 (pflotran_model)

    use clm_pflotran_interface_data
    use Realization_class
    use Patch_module
    use Grid_module
    !use clm_pflotran_interface_type
    use Richards_module
    use Richards_Aux_module
    use Global_Aux_module
    use Field_module

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(pflotran_model_type), pointer        :: pflotran_model
    type(realization_type),pointer            :: realization
    type(patch_type),pointer                  :: patch
    type(grid_type),pointer                   :: grid
    type(field_type),pointer                  :: field
    type(mapping_type),pointer                :: map
    type(richards_auxvar_type), pointer       :: aux_var
    !type(inside_each_pflotran_cell), pointer  :: pf_cell
    !type(inside_each_clm_cell),pointer        :: clm_cell
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    !PetscInt           :: i, g, lev
    !PetscReal          :: sat, vol_ovp, tmp
    !PetscReal,pointer  :: porosity_loc_p(:)
    PetscReal,pointer  :: sat_pf_p(:)
    PetscReal,pointer  :: sat_clm_p(:)

    realization     => pflotran_model%simulation%realization
    patch           => realization%patch
    grid            => patch%grid
    global_aux_vars => patch%aux%Global%aux_vars
    field           => realization%field

    call RichardsUpdateAuxVars(realization)
    !call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
    call VecGetArrayF90(clm_pf_idata%sat_pf, sat_pf_p, ierr)  

    !write(*,*),'sat_pf:'
    do local_id= 1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      sat_pf_p(local_id) = global_aux_vars(ghosted_id)%sat(1)
      !write(*,*),local_id,ghosted_id,sat_pf_p(local_id)
    enddo
    
    call VecRestoreArrayF90(clm_pf_idata%sat_pf, sat_pf_p, ierr)
    call MappingSourceToDestination(pflotran_model%map_pf2clm, &
                                    pflotran_model%option, &
                                    clm_pf_idata%sat_pf, &
                                    clm_pf_idata%sat_clm)
  end subroutine pflotranModelGetSaturation3

  ! ************************************************************************** !
  !
  !
  ! author: Gautam Bisht
  ! date: 9/22/2010
  ! ************************************************************************** !
  subroutine pflotranModelUpdateTopBCHomogeneous(pflotran_model, flux_value)

    use Condition_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    type(realization_type),pointer            :: realization
    type(condition_list_type),pointer         :: condition_list
    type(flow_condition_type), pointer        :: condition
    type(flow_sub_condition_type), pointer    :: sub_condition
    type(flow_condition_dataset_type),pointer :: dataset
    type(option_type),pointer                 :: option


    PetscReal                                 :: flux_value
    PetscInt                                  :: isub_condition

    realization => pflotran_model%realization
    condition_list => realization%flow_conditions
    option => realization%option

    condition => condition_list%first
    do
       if (.not.associated(condition)) exit

       if (trim(condition%name) == 'top') then

          do isub_condition = 1, condition%num_sub_conditions

             sub_condition => condition%sub_condition_ptr(isub_condition)%ptr

             if (associated(sub_condition)) then
                !dataset => sub_condition%dataset
                !dataset%values(1,1) = flux_value
                !dataset%rate = flux_value
                print *,'In UpdateTopBCHomogeneous: ', flux_value
                !call FlowSubConditionUpdateDataset(option,1.0d0,sub_condition%dataset)
             endif

          enddo
       endif

       condition => condition%next

    enddo

    !call StepperUpdateSolution(realization)

  end subroutine pflotranModelUpdateTopBCHomogeneous

  ! ************************************************************************** !
  !
  ! pflotranModelStepperRunFinalize: It performs the same execution of commands
  !             that are carried out in StepperRun() once the model integration is
  !             finished
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelStepperRunFinalize(pflotran_model)

    implicit none

    type(pflotran_model_type), pointer :: pflotran_model

    call StepperRunFinalize(pflotran_model%simulation%realization, &
         pflotran_model%simulation%flow_stepper, &
         pflotran_model%simulation%tran_stepper)

  end subroutine pflotranModelStepperRunFinalize



  ! ************************************************************************** !
  !
  ! pflotranModelInsertWaypoint: Inserts a waypoint within the waypoint list
  !             so that the model integration can be paused when that waypoint is
  !             reached
  !
  ! NOTE: It is assumed the 'waypoint_time' is in seconds
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelInsertWaypoint(pflotran_model, waypoint_time)

    implicit none

    type(pflotran_model_type),pointer :: pflotran_model
    type(waypoint_type)      ,pointer :: waypoint
    type(option_type)        ,pointer :: option
    PetscReal                         :: waypoint_time
    character(len=MAXWORDLENGTH)      :: word

    option => pflotran_model%realization%option
    word = 's'
    waypoint => WaypointCreate()
    waypoint%time              = waypoint_time * UnitsConvertToInternal(word,option)
    waypoint%update_conditions = PETSC_TRUE
    waypoint%dt_max            = 3153600

    call WaypointInsertInList(waypoint,pflotran_model%realization%waypoints)


  end subroutine pflotranModelInsertWaypoint

  subroutine pflotranModelDeleteWaypoint(pflotran_model, waypoint_time)

    type(pflotran_model_type),pointer :: pflotran_model
    type(waypoint_type)      ,pointer :: waypoint
    type(option_type)        ,pointer :: option
    PetscReal                         :: waypoint_time
    character(len=MAXWORDLENGTH)      :: word

    option => pflotran_model%realization%option
    word = 's'
    waypoint => WaypointCreate()
    waypoint%time              = waypoint_time * UnitsConvertToInternal(word,option)
    waypoint%update_conditions = PETSC_TRUE
    waypoint%dt_max            = 3153600

    call WaypointDeleteFromList(waypoint,pflotran_model%realization%waypoints)


  end subroutine pflotranModelDeleteWaypoint

  ! ************************************************************************** !
  !
  ! pflotranModelDestroy: Deallocates the pflotranModel object
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelDestroy(pflotran_model)

    implicit none

    type(pflotran_model_type), pointer :: pflotran_model
    PetscInt :: ii, jj
    PetscErrorCode :: ierr
#ifdef CLM_PFLOTRAN
    type(mapping_type),pointer         :: map

    !map         => pflotran_model%mapping

    !deallocate ( pf_clm_data%zwt      )
    !deallocate ( pf_clm_data%alpha    )
    !deallocate ( pf_clm_data%lambda   )
    !deallocate ( pf_clm_data%qsrc_flx )
    !deallocate ( pf_clm_data%sat_new  )
    !deallocate ( map%pf2clm )
    !deallocate ( map%clm2pf )
#endif

    ! Clean things up.
    call SimulationDestroy(pflotran_model%simulation)

    ! Final Time
    call PetscGetCPUTime(timex(2), ierr)
    call PetscTime(timex_wall(2), ierr)

    if (pflotran_model%option%myrank == pflotran_model%option%io_rank) then

       if (pflotran_model%option%print_to_screen) then
          if (masterproc) then
             write(iulog,'(/," CPU Time:", 1pe12.4, " [sec] ", &
                  & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
                  timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
                  (timex(2)-timex(1))/3600.d0

             write(iulog,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
                  & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
                  timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
                  (timex_wall(2)-timex_wall(1))/3600.d0
          endif
       endif
       if (pflotran_model%option%print_to_file) then
          write(pflotran_model%option%fid_out,'(/," CPU Time:", 1pe12.4, " [sec] ", &
               & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
               timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
               (timex(2)-timex(1))/3600.d0

          write(pflotran_model%option%fid_out,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
               & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
               timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
               (timex_wall(2)-timex_wall(1))/3600.d0
       endif
    endif

    if (pflotran_model%option%myrank == pflotran_model%option%io_rank .and. pflotran_model%option%print_to_file) &
         close(pflotran_model%option%fid_out)

    call LoggingDestroy()

    call PetscOptionsSetValue('-options_left','no',ierr);

    call OptionDestroy(pflotran_model%option)
    call PetscFinalize(ierr)

  end subroutine pflotranModelDestroy
  
end module pflotran_model_module

