module pflotran_model_module

  use Simulation_module
  use Realization_module
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
  use clm_varctl, only            : iulog
  use Mapping_module
#endif  
#if defined (CLM_PFLOTRAN) || defined(CLM_OFFLINE)  
  use Richards_Aux_module
  use pflotran_clm_interface_type
  use clm_pflotran_interface_data
#endif

#if defined CLM_PFLOTRAN
  use clm_pflotran_interface_type
#endif
  use Mapping_module

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  PetscLogDouble :: timex(4), timex_wall(4)

  PetscBool :: truth
  PetscBool :: option_found
  PetscBool :: single_inputfile
  PetscInt  :: i
  PetscInt  :: temp_int
  PetscErrorCode :: ierr
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
#ifdef WITH_CLM
     type(mapping_type),     pointer :: mapping
#endif
     PetscReal :: pause_time_1
     PetscReal :: pause_time_2
     type(inside_each_overlapped_cell), pointer :: pf_cells(:)
     type(inside_each_overlapped_cell), pointer :: clm_cells(:)
     type(mapping_type),                pointer :: map_clm2pf
     type(mapping_type),                pointer :: map_clm2pf_soils
     type(mapping_type),                pointer :: map_pf2clm

     PetscInt :: num_pf_cells
     PetscInt :: num_clm_cells

  end type pflotran_model_type

  public::pflotranModelCreate,               &
#if defined(CLM_PFLOTRAN)
       pflotranModelSetICs,                  & !
       pflotranModelSetSoilProp,             & !
       pflotranModelSetSoilProp2,            & !
       pflotranModelUpdateSourceSink,        & !
       pflotranModelUpdateSaturation,        & !
       pflotranModelGetSaturation,           & !
       pflotranModelInitMapping,             &
#endif
       pflotranModelInitMapping2,            &
       pflotranModelInitMapping3,            &
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
  function pflotranModelCreate()

    use Simulation_module
    use Realization_module
    use Timestepper_module
    use Option_module
    use Input_module
    use Init_module
    use Logging_module
    use Stochastic_module
    use Stochastic_Aux_module
#ifdef CLM_PFLOTRAN
    use pflotran_clm_interface_type
#endif
    implicit none


    type(pflotran_model_type), pointer :: pflotranModelCreate


    PetscLogDouble :: timex(4), timex_wall(4)

    PetscBool  :: truth
    PetscBool  :: option_found
    PetscBool  :: single_inputfile
    PetscInt   :: i
    PetscInt   :: temp_int
    PetscErrorCode :: ierr
    character(len=MAXSTRINGLENGTH)          :: string
    character(len=MAXSTRINGLENGTH), pointer :: filenames(:)

    type(pflotran_model_type),      pointer :: pflotran_model


    allocate(pflotran_model)
    allocate(pflotran_model%stochastic)
    allocate(pflotran_model%simulation)
    allocate(pflotran_model%realization)
    allocate(pflotran_model%option)
#ifdef CLM_PFLOTRAN
    allocate(pflotran_model%mapping)
#endif
    allocate(pflotran_model%map_clm2pf)
    allocate(pflotran_model%map_clm2pf_soils)
    allocate(pflotran_model%map_pf2clm)

    nullify(pflotran_model%stochastic)
    nullify(pflotran_model%simulation)
    nullify(pflotran_model%realization)
    nullify(pflotran_model%option)
#ifdef CLM_PFLOTRAN
    nullify(pflotran_model%mapping)
#endif
    nullify(pflotran_model%pf_cells)
    nullify(pflotran_model%clm_cells)
    nullify(pflotran_model%map_clm2pf)
    nullify(pflotran_model%map_clm2pf_soils)
    nullify(pflotran_model%map_pf2clm)

    pflotran_model%option => OptionCreate()
    pflotran_model%option%fid_out = IUNIT2
    single_inputfile = PETSC_TRUE

    pflotran_model%pause_time_1 = -1.0d0
    pflotran_model%pause_time_2 = -1.0d0


#ifndef CLM_PFLOTRAN
    call MPI_Init(ierr)
#endif
    pflotran_model%option%global_comm = MPI_COMM_WORLD
    call MPI_Comm_rank(MPI_COMM_WORLD,pflotran_model%option%global_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,pflotran_model%option%global_commsize,ierr)
    call MPI_Comm_group(MPI_COMM_WORLD,pflotran_model%option%global_group,ierr)
    pflotran_model%option%mycomm = pflotran_model%option%global_comm
    pflotran_model%option%myrank = pflotran_model%option%global_rank
    pflotran_model%option%mycommsize = pflotran_model%option%global_commsize
    pflotran_model%option%mygroup = pflotran_model%option%global_group


    ! check for non-default input filename
    pflotran_model%option%input_filename = "pflotran.in"
    string = '-pflotranin'
    call InputGetCommandLineString(string,pflotran_model%option%input_filename,option_found,pflotran_model%option)

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
#ifdef CLM_PFLOTRAN
       write(iulog,*),'single_inputfile'
#else
       if(pflotran_model%option%myrank == pflotran_model%option%io_rank) write(*,*),'single_inputfile'
#endif
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
    call PetscGetTime(timex_wall(1), ierr)
    pflotran_model%option%start_time = timex_wall(1)

    call Init(pflotran_model%simulation)
#if defined (CLM_PFLOTRAN) || defined (CLM_OFFLINE)
    select case( pflotran_model%realization%patch%grid%itype)
      case(STRUCTURED_GRID)
        ! dimension = nlxy
        allocate ( pf_clm_data%zwt     ( pflotran_model%realization%patch%grid%structured_grid%nlxy ) )

        ! dimension = nlmax
        allocate ( pf_clm_data%qsrc_flx( pflotran_model%realization%patch%grid%structured_grid%nlmax) )
        allocate ( pf_clm_data%sat_new ( pflotran_model%realization%patch%grid%structured_grid%nlmax) )

        ! dimension = ngmax
        allocate ( pf_clm_data%alpha   ( pflotran_model%realization%patch%grid%ngmax ) )
        allocate ( pf_clm_data%lambda  ( pflotran_model%realization%patch%grid%ngmax ) )
      case(UNSTRUCTURED_GRID)
        ! dimension = nlxy
        !allocate ( pf_clm_data%zwt     ( pflotran_model%realization%patch%grid%unstructured_grid%nlxy ) )
        nullify(pf_clm_data%zwt)

        ! dimension = nlmax
        allocate ( pf_clm_data%qsrc_flx( pflotran_model%realization%patch%grid%unstructured_grid%nlmax) )
        allocate ( pf_clm_data%sat_new ( pflotran_model%realization%patch%grid%unstructured_grid%nlmax) )        

        ! dimension = ngmax
        allocate ( pf_clm_data%alpha   ( pflotran_model%realization%patch%grid%ngmax ) )
        allocate ( pf_clm_data%lambda  ( pflotran_model%realization%patch%grid%ngmax ) )
        
        pf_clm_data%qsrc_flx = 0.0d0
        pf_clm_data%sat_new  = 0.0d0
        pf_clm_data%alpha    = 0.0d0
        pf_clm_data%lambda   = 0.0d0
        
        !write(*,*), 'size(pf_clm_data%alpha   ) = ',pflotran_model%realization%patch%grid%ngmax
        !write(*,*), 'size(pf_clm_data%qsrc_flx) = ',pflotran_model%realization%patch%grid%unstructured_grid%nlmax
      
    end select
#endif

#ifdef CLM_PFLOTRAN
    pflotran_model%mapping => MappingCreate()
#endif
    pflotran_model%map_clm2pf       => MappingCreate()
    pflotran_model%map_clm2pf_soils => MappingCreate()
    pflotran_model%map_pf2clm       => MappingCreate()

    pflotran_model%num_pf_cells = -1
    pflotran_model%num_clm_cells= -1

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

    call StepperRunInit(pflotran_model%simulation%realization, &
         pflotran_model%simulation%flow_stepper, &
         pflotran_model%simulation%tran_stepper)

  end subroutine pflotranModelStepperRunInit

  ! ************************************************************************** !
  !
  ! pflotranModelSetICs:
  !
  !
  ! author: Gautam Bisht
  ! date: 10/12/2010
  ! ************************************************************************** !
#ifdef CLM_PFLOTRAN
  subroutine pflotranModelSetICs( pflotran_model )

    use Realization_module
    use Patch_module
    use Grid_module
    use Field_module
    use Richards_Aux_module
    use Global_Aux_module
    use Option_module
    use Richards_module
    use Discretization_module

    implicit none
    type(pflotran_model_type), pointer        :: pflotran_model
    type(option_type), pointer                :: option
    type(realization_type),pointer            :: realization
    type(patch_type),pointer                  :: patch
    type(grid_type),pointer                   :: grid
    type(field_type),pointer                  :: field
    type(richards_auxvar_type), pointer       :: rich_aux_vars(:)
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    type(discretization_type), pointer        :: discretization


    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal          :: head_constant, grav, den,tmp, sat
    PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:)
    PetscReal, pointer :: porosity_loc_p(:), perm_xx_loc_p(:)

    realization      => pflotran_model%simulation%realization
    discretization   => realization%discretization
    patch            => realization%patch
    grid             => patch%grid
    field            => realization%field
    rich_aux_vars    => patch%aux%Richards%aux_vars
    global_aux_vars  => patch%aux%Global%aux_vars
    option           => realization%option

    write (iulog,*), 'In pflotranModelSetICs'

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    grav = 9.81d0       ! [m/S^2]

    head_constant = den*grav*(grid%z_max_global - clm_pf_data%zwt(1)) + 101325.0d0

    write(iulog,*), 'z_max_global = ',grid%z_max_global
    write(iulog,*), 'zwt          = ',clm_pf_data%zwt(1)
    write(iulog,*), 'head_constant= ',head_constant

    call GridVecGetArrayF90(grid,field%flow_xx      ,xx_loc_p       ,ierr)

    do local_id = 1, grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
       xx_loc_p(ghosted_id) = head_constant - grid%z(ghosted_id)*den*grav
       !xx_loc_p(ghosted_id) = 0.50d0*xx_loc_p(ghosted_id)
    enddo

    call GridVecRestoreArrayF90(grid,field%flow_xx       ,xx_loc_p       ,ierr)
    ! update dependent vectors
    call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
         field%flow_xx_loc,NFLOWDOF)
    call VecCopy(field%flow_xx, field%flow_yy, ierr)

    call RichardsUpdateAuxVars(realization)

    call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
    tmp = 0.d0
    write(iulog,*), 'Saturation:'
    do local_id= 1,grid%nlmax

       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
       sat      = global_aux_vars(ghosted_id)%sat(1)
       tmp = tmp + sat*porosity_loc_p(ghosted_id)
       write(iulog,*), local_id,sat
    enddo
    write(iulog, *),'total_volume_liq = ', tmp

  end subroutine pflotranModelSetICs


  subroutine pflotranModelSetICs2( pflotran_model )

    use Realization_module
    use Patch_module
    use Grid_module
    use Field_module
    use Richards_Aux_module
    use Global_Aux_module
    use Option_module
    use Richards_module
    use Discretization_module

    implicit none
    type(pflotran_model_type), pointer        :: pflotran_model
    type(option_type), pointer                :: option
    type(realization_type),pointer            :: realization
    type(patch_type),pointer                  :: patch
    type(grid_type),pointer                   :: grid
    type(field_type),pointer                  :: field
    type(richards_auxvar_type), pointer       :: rich_aux_vars(:)
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    type(discretization_type), pointer        :: discretization


    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal          :: head_constant, grav, den,tmp, sat
    PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:)
    PetscReal, pointer :: porosity_loc_p(:), perm_xx_loc_p(:)

    realization      => pflotran_model%simulation%realization
    discretization   => realization%discretization
    patch            => realization%patch
    grid             => patch%grid
    field            => realization%field
    rich_aux_vars    => patch%aux%Richards%aux_vars
    global_aux_vars  => patch%aux%Global%aux_vars
    option           => realization%option

    write (iulog,*), 'In pflotranModelSetICs'

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    grav = 9.81d0       ! [m/S^2]

    head_constant = den*grav*(grid%z_max_global - clm_pf_data%zwt(1)) + 101325.0d0

    call GridVecGetArrayF90(grid,field%flow_xx      ,xx_loc_p       ,ierr)

    do local_id = 1, grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
       xx_loc_p(ghosted_id) = head_constant - grid%z(ghosted_id)*den*grav
    enddo

    call GridVecRestoreArrayF90(grid,field%flow_xx       ,xx_loc_p       ,ierr)
    ! update dependent vectors
    call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
         field%flow_xx_loc,NFLOWDOF)
    call VecCopy(field%flow_xx, field%flow_yy, ierr)

    call RichardsUpdateAuxVars(realization)

    call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
    tmp = 0.d0
    do local_id= 1,grid%nlmax

       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
       sat      = global_aux_vars(ghosted_id)%sat(1)
       tmp = tmp + sat*porosity_loc_p(ghosted_id)
    enddo

  end subroutine pflotranModelSetICs2

#endif

  ! ************************************************************************** !
  !
  ! pflotranModelSetSoilProp:
  !
  !
  ! author: Gautam Bisht
  ! date: 10/22/2010
  ! ************************************************************************** !
#ifdef CLM_PFLOTRAN
  subroutine pflotranModelSetSoilProp( pflotran_model )

    use Realization_module
    use Patch_module
    use Grid_module
    use clm_pflotran_interface_type
    use Richards_Aux_module
    use Field_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    type(realization_type),pointer            :: realization
    type(patch_type),pointer                  :: patch
    type(grid_type),pointer                   :: grid
    type(field_type),pointer                  :: field
    type(mapping_type),pointer                :: map
    type(richards_auxvar_type), pointer       :: rich_aux_vars(:)
    type(richards_auxvar_type), pointer       :: aux_var
    type(inside_each_pflotran_cell), pointer  :: pf_cell


    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id, clm_id
    PetscInt           :: i
    PetscReal          :: den, vis, grav, vol_ovlap
    PetscReal,pointer  :: porosity_loc_p(:), vol_ovlap_arr(:)
    PetscReal,pointer  :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)


    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav = 9.81d0       ! [m/S^2]

    realization     => pflotran_model%simulation%realization
    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field
    rich_aux_vars   => patch%aux%Richards%aux_vars
    map             => pflotran_model%mapping


    call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
    call GridVecGetArrayF90(grid,field%perm_xx_loc,  perm_xx_loc_p,  ierr)
    call GridVecGetArrayF90(grid,field%perm_yy_loc,  perm_yy_loc_p,  ierr)
    call GridVecGetArrayF90(grid,field%perm_zz_loc,  perm_zz_loc_p,  ierr)

    allocate( vol_ovlap_arr(1:grid%ngmax ))

    do local_id = 1,grid%ngmax
       vol_ovlap_arr(local_id) = 0.0d0
    enddo


    write (iulog,*), 'In pflotranModelSetSoilProp: ',grid%ngmax


    ! Initialize permeability and porosity to ZERO
    do local_id = 1,grid%nlmax

       ghosted_id = grid%nL2G(local_id)

       if (associated(patch%imat)) then
          if(patch%imat(ghosted_id) <= 0) cycle
       endif

       perm_xx_loc_p(ghosted_id) = 0.0d0
       perm_yy_loc_p(ghosted_id) = 0.0d0
       perm_zz_loc_p(ghosted_id) = 0.0d0

       porosity_loc_p(ghosted_id) = 0.0d0

    enddo


    ! Map
    do local_id = 1,grid%nlmax

       ghosted_id = grid%nL2G(local_id)

       if (associated(patch%imat)) then
          if(patch%imat(ghosted_id) <= 0) cycle
       endif


       aux_var => rich_aux_vars(ghosted_id)
       pf_cell => map%pf2clm(local_id)

       !write(iulog,*), 'ghosted_id',ghosted_id, pf_cell%num_clm_cells

       do i = 1,pf_cell%num_clm_cells

          clm_id    = pf_cell%id_clm_cells(i)
          vol_ovlap = pf_cell%perc_vol_overlap(i)

          !
          ! bc_alpha [1/Pa]; while sucsat [mm of H20]
          !
          ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
          !
          aux_var%bc_alpha  = aux_var%bc_alpha  + vol_ovlap * 1.d0/(clm_pf_data%sucsat(1,clm_id)*grav)

          !
          ! bc_lambda = 1/bsw
          !
          aux_var%bc_lambda = aux_var%bc_lambda + vol_ovlap * (1.d0/clm_pf_data%bsw(   1,clm_id))

          !
          ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
          ! [m^2]          [mm/sec]
          perm_xx_loc_p(ghosted_id) = perm_xx_loc_p(ghosted_id) + vol_ovlap * clm_pf_data%hksat_x(1,clm_id) * vis/ (den * grav) / 1000.d0
          perm_yy_loc_p(ghosted_id) = perm_yy_loc_p(ghosted_id) + vol_ovlap * clm_pf_data%hksat_y(1,clm_id) * vis/ (den * grav) / 1000.d0
          perm_zz_loc_p(ghosted_id) = perm_zz_loc_p(ghosted_id) + vol_ovlap * clm_pf_data%hksat_z(1,clm_id) * vis/ (den * grav) / 1000.d0

          !
          ! porosity = vol. soil moisture @ saturation
          !
          porosity_loc_p(ghosted_id) = porosity_loc_p(ghosted_id) + vol_ovlap * clm_pf_data%watsat(1,clm_id)

          vol_ovlap_arr(ghosted_id)  = vol_ovlap_arr(ghosted_id) + vol_ovlap

       enddo

       if (pf_cell%num_clm_cells.eq.0) then

          !
          ! Soil layer in PFLOTRAN is lower than CLM's lowest soil layer,
          ! use the values of lowest soil layer
          !

          aux_var%bc_alpha  = 1.d0/(clm_pf_data%sucsat(1,clm_pf_data%nlevsoi)*grav)
          aux_var%bc_lambda = 1.d0/(clm_pf_data%bsw(   1,clm_pf_data%nlevsoi))

          perm_xx_loc_p(ghosted_id) = clm_pf_data%hksat_x(1,clm_pf_data%nlevsoi) * vis/ (den * grav) / 1000.d0
          perm_yy_loc_p(ghosted_id) = clm_pf_data%hksat_y(1,clm_pf_data%nlevsoi) * vis/ (den * grav) / 1000.d0
          perm_zz_loc_p(ghosted_id) = clm_pf_data%hksat_z(1,clm_pf_data%nlevsoi) * vis/ (den * grav) / 1000.d0

          porosity_loc_p(ghosted_id) = clm_pf_data%watsat(1,clm_pf_data%nlevsoi)

          vol_ovlap_arr(ghosted_id) = 1.0d0

       endif

    enddo

    ! Scale the soil properties by total_vol_overlap
    do local_id = 1,grid%nlmax

       ghosted_id = grid%nL2G(local_id)

       if (associated(patch%imat)) then
          if(patch%imat(ghosted_id) <= 0) cycle
       endif

       aux_var => rich_aux_vars(ghosted_id)

       if (map%pf2clm(local_id)%num_clm_cells.gt.0) then

          !aux_var%bc_alpha  = aux_var%bc_alpha /vol_ovlap_arr(ghosted_id)
          !aux_var%bc_lambda = aux_var%bc_lambda/vol_ovlap_arr(ghosted_id)
          !
          !perm_xx_loc_p(ghosted_id) = perm_xx_loc_p(ghosted_id) / vol_ovlap_arr(ghosted_id)
          !perm_yy_loc_p(ghosted_id) = perm_yy_loc_p(ghosted_id) / vol_ovlap_arr(ghosted_id)
          !perm_zz_loc_p(ghosted_id) = perm_zz_loc_p(ghosted_id) / vol_ovlap_arr(ghosted_id)
          !
          !porosity_loc_p(ghosted_id) = porosity_loc_p(ghosted_id) / vol_ovlap_arr(ghosted_id)

          aux_var%bc_alpha  = aux_var%bc_alpha /map%pf2clm(local_id)%total_vol_overlap
          aux_var%bc_lambda = aux_var%bc_lambda/map%pf2clm(local_id)%total_vol_overlap

          perm_xx_loc_p(ghosted_id) = perm_xx_loc_p(ghosted_id) / map%pf2clm(local_id)%total_vol_overlap
          perm_yy_loc_p(ghosted_id) = perm_yy_loc_p(ghosted_id) / map%pf2clm(local_id)%total_vol_overlap
          perm_zz_loc_p(ghosted_id) = perm_zz_loc_p(ghosted_id) / map%pf2clm(local_id)%total_vol_overlap

          porosity_loc_p(ghosted_id) = porosity_loc_p(ghosted_id) / map%pf2clm(local_id)%total_vol_overlap

       endif

    enddo


    call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
    call GridVecRestoreArrayF90(grid,field%perm_xx_loc,  perm_xx_loc_p,  ierr)
    call GridVecRestoreArrayF90(grid,field%perm_yy_loc,  perm_yy_loc_p,  ierr)
    call GridVecRestoreArrayF90(grid,field%perm_zz_loc,  perm_zz_loc_p,  ierr)


    deallocate ( vol_ovlap_arr)

  end subroutine pflotranModelSetSoilProp


  subroutine pflotranModelSetSoilProp2( pflotran_model )

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model

  end subroutine pflotranModelSetSoilProp2


#endif

  ! ************************************************************************** !
  !
  ! pflotranModelGetSaturation:
  !
  !
  ! author: Gautam Bisht
  ! date: 11/01/2010
  ! ************************************************************************** !
#if defined(CLM_PFLOTRAN)
  subroutine pflotranModelGetSaturation (pflotran_model)

    use Realization_module
    use Patch_module
    use Grid_module
    use clm_pflotran_interface_type
    use Richards_module
    use Richards_Aux_module
    use Global_Aux_module
    use Field_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    type(realization_type),pointer            :: realization
    type(patch_type),pointer                  :: patch
    type(grid_type),pointer                   :: grid
    type(field_type),pointer                  :: field
    type(mapping_type),pointer                :: map
    type(richards_auxvar_type), pointer       :: aux_var
    type(inside_each_pflotran_cell), pointer  :: pf_cell
    type(inside_each_clm_cell),pointer        :: clm_cell
    type(global_auxvar_type), pointer         :: global_aux_vars(:)

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id, clm_id
    PetscInt           :: i, g, lev
    PetscReal          :: sat, vol_ovp, tmp
    PetscReal,pointer  :: porosity_loc_p(:)

    realization     => pflotran_model%simulation%realization
    patch           => realization%patch
    grid            => patch%grid
    map             => pflotran_model%mapping
    global_aux_vars => patch%aux%Global%aux_vars
    field           => realization%field

    write(iulog, *), 'pflotranModelGetSaturation'

    call RichardsUpdateAuxVars(realization)
    call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)

    do g = 1, clm_pf_data%ngrids
       do lev = 1, clm_pf_data%nlevsoi
          clm_pf_data%sat(g,lev) = 0.0d0
       enddo
    enddo

    do clm_id = 1,clm_pf_data%nlevsoi

       clm_cell => map%clm2pf(clm_id)

       do i = 1,clm_cell%num_pflotran_cells

          local_id = clm_cell%id_pflotran_cells(i)
          vol_ovp  = clm_cell%perc_vol_overlap( i)

          ghosted_id = grid%nL2G(local_id)
          sat      = global_aux_vars(ghosted_id)%sat(1)

          clm_pf_data%sat(1,clm_id) = clm_pf_data%sat(1,clm_id) + sat*vol_ovp
       enddo

       clm_pf_data%sat(1,clm_id) = min( clm_pf_data%sat(1,clm_id), 1.0d0)

    enddo


    call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)

  end subroutine pflotranModelGetSaturation
#endif

  ! ************************************************************************** !
  !
  ! pflotranModelInitMapping:
  !
  !
  ! author: Gautam Bisht
  ! date: 10/18/2010
  ! ************************************************************************** !
#ifdef CLM_PFLOTRAN
  subroutine pflotranModelInitMapping( pflotran_model)

    use Realization_module
    use Patch_module
    use Grid_module
    use clm_pflotran_interface_type

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    type(realization_type),pointer            :: realization
    type(patch_type),pointer                  :: patch
    type(grid_type),pointer                   :: grid
    type(mapping_type),pointer                :: map

    PetscInt  :: local_id, ghosted_id, lev, i,j, id, tmp_int
    PetscInt  :: z_up_lev, z_dw_lev
    PetscReal :: z_frm_top, zsurf_frm_top,z_up_frm_top,z_dw_frm_top;
    PetscReal :: tmp_a, tmp_b

    realization => pflotran_model%simulation%realization
    patch       => realization%patch
    grid        => patch%grid
    map         => pflotran_model%mapping


    allocate(map%pf2clm(1:grid%nlmax                               ))
    allocate(map%clm2pf(1:clm_pf_data%nlevsoi * clm_pf_data%ngrids ))

    write (iulog,*), 'In pflotranModelInitMapping nlmax = ',grid%nlmax
    do i = 1, clm_pf_data%nlevsoi
       map%clm2pf(i)%num_pflotran_cells = 0
       map%clm2pf(i)%total_vol_overlap  = 0.0d0
    enddo



    do local_id = 1, grid%nlmax

       ghosted_id = grid%nL2G(local_id)

       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif


       !                                     --------
       !    ===================    z_up_frm_top  |
       !    ||               ||                  |
       !    ||               ||                 \ /  +ve
       !    ||               ||
       !    ||       o       ||    z_frm_top
       !    ||               ||
       !    ||               ||
       !    ||               ||
       !    ===================    z_dw_frm_top
       !
       !

       z_frm_top       = grid%z_max_global - grid%z(ghosted_id)
       z_up_frm_top    = z_frm_top - grid%structured_grid%dz(ghosted_id) * 0.5d0
       z_dw_frm_top    = z_frm_top + grid%structured_grid%dz(ghosted_id) * 0.5d0

       if ( z_up_frm_top.lt.0.d0) then
          z_up_frm_top = 0.0d0
       endif


       if( z_up_frm_top.ge.clm_pf_data%zisoi(clm_pf_data%nlevsoi)) then
          ! Top of the CV cell in the PFLOTRAN is below the lowest soil layer witin the CLM
          map%pf2clm(local_id)%num_clm_cells = 0;
       else


          do lev = 0,clm_pf_data%nlevsoi-1
             z_up_lev = lev+1
             if( (clm_pf_data%zisoi(lev).le.z_up_frm_top).and.(clm_pf_data%zisoi(lev+1).ge.z_up_frm_top )) then
                goto 10
             endif
          enddo

10        continue

          do lev = 0,clm_pf_data%nlevsoi-1
             z_dw_lev = lev+1
             if( (clm_pf_data%zisoi(lev).le.z_dw_frm_top).and.(clm_pf_data%zisoi(lev+1).ge.z_dw_frm_top )) then
                goto 20
             endif
          enddo
20        continue

          map%pf2clm(local_id)%num_clm_cells     = z_dw_lev - z_up_lev + 1d0
          map%pf2clm(local_id)%total_vol_overlap = 0.0d0
          allocate( map%pf2clm(local_id)%id_clm_cells(    1:z_dw_lev - z_up_lev + 1 ) )
          allocate( map%pf2clm(local_id)%perc_vol_overlap(1:z_dw_lev - z_up_lev + 1 ) )

          do lev = 1, map%pf2clm(local_id)%num_clm_cells

             map%pf2clm(local_id)%id_clm_cells(lev) = z_up_lev + (lev-1)

             tmp_int = map%clm2pf( z_up_lev + lev-1 )%num_pflotran_cells + 1
             map%clm2pf( z_up_lev + lev-1 )%num_pflotran_cells =  &
                  map%clm2pf( z_up_lev + lev-1 )%num_pflotran_cells + 1

             tmp_a = max( clm_pf_data%zisoi(z_up_lev + lev - 2), z_up_frm_top)
             tmp_b = min( clm_pf_data%zisoi(z_up_lev + lev - 1), z_dw_frm_top)

             map%pf2clm(local_id)%perc_vol_overlap(lev) = (tmp_b-tmp_a)/(z_dw_frm_top - z_up_frm_top)

             map%pf2clm(local_id)%total_vol_overlap = map%pf2clm(local_id)%total_vol_overlap + &
                  map%pf2clm(local_id)%perc_vol_overlap(lev)
          enddo

       endif

    enddo


    do i = 1,clm_pf_data%nlevsoi

       if (map%clm2pf(i)%num_pflotran_cells.gt.0) then
          allocate(map%clm2pf(i)%id_pflotran_cells(1:map%clm2pf(i)%num_pflotran_cells))
          allocate(map%clm2pf(i)%perc_vol_overlap( 1:map%clm2pf(i)%num_pflotran_cells))

          do j=1,map%clm2pf(i)%num_pflotran_cells
             map%clm2pf(i)%id_pflotran_cells(j) = -999
          enddo
       endif
    enddo



    do local_id = 1, grid%nlmax

       ghosted_id = grid%nL2G(local_id)

       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif

       if ( map%pf2clm(local_id)%num_clm_cells.gt.0) then

          do i=1, map%pf2clm(local_id)%num_clm_cells

             id = map%pf2clm(local_id)%id_clm_cells(i)

             do j=1, map%clm2pf( id )%num_pflotran_cells
                lev = j
                if ( map%clm2pf( id )%id_pflotran_cells(j).eq.-999) then
                   goto 30
                endif
             enddo
30           continue


             z_frm_top       = grid%z_max_global - grid%z(ghosted_id)
             z_up_frm_top    = z_frm_top - grid%structured_grid%dz(ghosted_id) * 0.5d0
             z_dw_frm_top    = z_frm_top + grid%structured_grid%dz(ghosted_id) * 0.5d0

             tmp_a = max( clm_pf_data%zisoi(id-1), z_up_frm_top)
             tmp_b = min( clm_pf_data%zisoi(id  ), z_dw_frm_top)

             map%clm2pf(id)%id_pflotran_cells(lev) = local_id
             map%clm2pf(id)%perc_vol_overlap(lev)  = (tmp_b - tmp_a)/&
                  ( clm_pf_data%zisoi(id) - clm_pf_data%zisoi(id-1))
             map%clm2pf(id)%total_vol_overlap = map%clm2pf(id)%total_vol_overlap + &
                  map%clm2pf(id)%perc_vol_overlap(lev)

          enddo
       endif
    enddo


    !write(iulog,*), 'CLM2PF: '
    do i = 1,clm_pf_data%nlevsoi
       !write(iulog,*), ' clm_id =',i,' num_cells = ', map%clm2pf(i)%num_pflotran_cells
       do j=1,map%clm2pf(i)%num_pflotran_cells
          !write(iulog,*), ' pf_id=',map%clm2pf(i)%id_pflotran_cells(j), &
          !     map%clm2pf(i)%perc_vol_overlap(j)
       enddo
    enddo

    write(iulog,*), '=============================================== '
    write(iulog,*), '=============================================== '
    write(iulog,*), 'PF2CLM: '
    do i = 1,grid%nlmax
       !write(iulog,*), ' pf_id =',i,' num_cells = ', map%pf2clm(i)%num_clm_cells
       do j=1,map%pf2clm(i)%num_clm_cells
          !write(iulog,*), ' clm_id=',map%pf2clm(i)%id_clm_cells(j), &
          !     map%pf2clm(i)%perc_vol_overlap(j)
       enddo
    enddo
    write(iulog,*), '=============================================== '


  end subroutine pflotranModelInitMapping
#endif

  ! ************************************************************************** !
  !
  ! pflotranModelInitMapping2:
  !
  !
  ! author: Gautam Bisht
  ! date: 01/28/2011
  ! ************************************************************************** !
  !#ifdef CLM_PFLOTRAN
  subroutine pflotranModelInitMapping2( pflotran_model, clm_pf_data, &
       grid_clm_cell_ids_ghosted_nindex, grid_clm_npts_local)

    use Input_module
    use Option_module
    use Realization_module
    use Grid_module
    use Patch_module
    use clm_pflotran_interface_data

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

    PetscInt, intent(in),pointer       :: grid_clm_cell_ids_ghosted_nindex(:)
    PetscInt, intent(in)               :: grid_clm_npts_local

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    type(clm_pflotran_data), intent(inout),pointer    :: clm_pf_data

    ! Local variables:
    PetscInt                           :: ii,jj, irank
    PetscInt                           :: local_id, clm_ocell_count, pf_ocell_count
    PetscInt, pointer                  :: grid_pf_cell_ids_ghosted_nindex(:)
    PetscInt, pointer                  :: grid_clm_cell_ids_ghosted_nindex_tmp(:)
    PetscInt                           :: grid_pf_npts_local, grid_pf_npts_ghost
  PetscInt                           :: grid_clm_npts_ghost
    PetscInt                           :: grid_pf_npts_ghosted, grid_clm_npts_ghosted
    PetscInt                           :: grid_clm_npts_ghosted_tmp
    PetscInt, pointer                  :: pf_ocell_ids(:),clm_ocell_ids(:)
    PetscReal,pointer                  :: pf_ocell_vol(:),clm_ocell_vol(:)
    PetscInt, pointer                  :: grid_clm_active(:), grid_pf_active(:)
  PetscInt, pointer                  :: pf_ocell_count_array(:), clm_ocell_count_array(:)

    character(len=MAXSTRINGLENGTH)     :: filename
    type(option_type), pointer         :: option
    type(realization_type), pointer    :: realization
    type(grid_type),pointer            :: grid
    type(patch_type), pointer          :: patch

    PetscErrorCode :: ierr
    PetscMPIInt    :: status_mpi(MPI_STATUS_SIZE)


    PetscInt, allocatable :: tmp_int_array(:)

    PetscViewer :: viewer
    Vec :: hksat_x, hksat_y, hksat_z
    Vec :: hksat_x_clmloc, hksat_y_clmloc, hksat_z_clmloc
    Vec :: hksat_x_pfloc, hksat_y_pfloc, hksat_z_pfloc

    type(mapping_type),     pointer :: m_clm2pf
    type(mapping_type),     pointer :: m_pf2clm

#if 0
    m_clm2pf        => pflotran_model%map_clm2pf
    m_pf2clm        => pflotran_model%map_pf2clm
    option          => pflotran_model%option
    realization     => pflotran_model%simulation%realization
    patch           => realization%patch
    patch           => realization%patch
    grid            => patch%grid

   if ( option%myrank == option%io_rank) then
   
       ! Read the mapping files
     
       filename = 'mapping_clm2pf.txt'
       call pflotranModelReadMappingFile( option, filename, &
            pflotran_model%clm_cells, pflotran_model%num_clm_cells)

       filename = 'mapping_pf2clm.txt'
       call pflotranModelReadMappingFile( option, filename, &
            pflotran_model%pf_cells, pflotran_model%num_pf_cells)

       ! TODO: Check the two mapping files are consistent
       !       (i.e. ensure that the cell-ids in one mapping file
       !        are present in the other mapping file and vice-versa)
       allocate(grid_clm_active(pflotran_model%num_clm_cells))
       allocate(grid_pf_active( pflotran_model%num_pf_cells ))
       grid_clm_active = 0
       grid_pf_active  = 0
   
       ! Receiving/Sending data from/to other processors
       do irank = 0,option%mycommsize-1
          if (irank == option%io_rank) cycle

          ! Get the "grid_clm_npts_ghosted"
          call MPI_Recv(grid_clm_npts_ghosted, 1, MPI_INTEGER, &
               irank,MPI_ANY_TAG,option%mycomm,status_mpi,ierr)

          ! Allocate memory         
          allocate(grid_clm_cell_ids_ghosted_nindex_tmp(grid_clm_npts_ghosted))
      allocate(clm_ocell_count_array            (grid_clm_npts_ghosted))
      
          call MPI_Recv(grid_clm_cell_ids_ghosted_nindex_tmp, grid_clm_npts_ghosted, &
               MPI_INTEGER,irank,MPI_ANY_TAG,option%mycomm,status_mpi,ierr)

          ! 
          do local_id = 1,grid_clm_npts_ghosted
        grid_clm_active( grid_clm_cell_ids_ghosted_nindex_tmp(local_id) + 1 ) = 1
      enddo

          call pflotranModelFindOverlapCells(pflotran_model%clm_cells, grid_clm_npts_ghosted,&
               grid_clm_cell_ids_ghosted_nindex_tmp, clm_ocell_ids,clm_ocell_vol,clm_ocell_count,&
         clm_ocell_count_array)

          call MPI_Send(clm_ocell_count,1,MPI_INTEGER,irank,0,option%mycomm,ierr)
          call MPI_Send(clm_ocell_ids,clm_ocell_count,MPI_INTEGER,irank,0,option%mycomm,ierr)
          call MPI_Send(clm_ocell_vol,clm_ocell_count,MPI_DOUBLE_PRECISION,irank,0,option%mycomm,ierr)

          ! Free memory
          deallocate(grid_clm_cell_ids_ghosted_nindex_tmp)
          deallocate(clm_ocell_count_array)
          deallocate(clm_ocell_ids)
          deallocate(clm_ocell_vol)


          call MPI_Recv(grid_pf_npts_ghosted,1,MPI_INTEGER, &
               irank,MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
         
      ! Allocate data
          allocate(grid_pf_cell_ids_ghosted_nindex(grid_pf_npts_ghosted))
      allocate(pf_ocell_count_array            (grid_pf_npts_ghosted))

          call MPI_Recv(grid_pf_cell_ids_ghosted_nindex, grid_pf_npts_ghosted, &
               MPI_INTEGER,irank,MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
          
      !
          do local_id = 1,grid_pf_npts_ghosted
        grid_pf_active( grid_pf_cell_ids_ghosted_nindex(local_id) + 1 ) = 1
      enddo

          call pflotranModelFindOverlapCells(pflotran_model%pf_cells, grid_pf_npts_ghosted,&
               grid_pf_cell_ids_ghosted_nindex,pf_ocell_ids,pf_ocell_vol,pf_ocell_count,&
         pf_ocell_count_array)

          call MPI_Send(pf_ocell_count,1,MPI_INTEGER,irank,0,option%mycomm,ierr)
          call MPI_Send(pf_ocell_ids,pf_ocell_count,MPI_INTEGER,irank,0,option%mycomm,ierr)
          call MPI_Send(pf_ocell_vol,pf_ocell_count,MPI_DOUBLE_PRECISION,irank,0,option%mycomm,ierr)

          ! Free memory
          deallocate(pf_ocell_count_array)
          deallocate(pf_ocell_ids)
          deallocate(pf_ocell_vol)

       enddo
     
       ! Doing local work
       grid_clm_npts_ghost   = 0
       grid_clm_npts_ghosted = grid_clm_npts_local + grid_clm_npts_ghost
       do local_id = 1,grid_clm_npts_ghosted
      grid_clm_active( grid_clm_cell_ids_ghosted_nindex(local_id) + 1 ) = 1
       enddo
     
       allocate(clm_ocell_count_array(grid_clm_npts_ghosted))
       call pflotranModelFindOverlapCells(pflotran_model%clm_cells, grid_clm_npts_ghosted,&
            grid_clm_cell_ids_ghosted_nindex, clm_ocell_ids,clm_ocell_vol,clm_ocell_count,&
      clm_ocell_count_array)

       allocate(grid_pf_cell_ids_ghosted_nindex(grid%ngmax))
       allocate(pf_ocell_count_array            (grid%ngmax))

       do local_id = 1,grid%ngmax
          grid_pf_cell_ids_ghosted_nindex(local_id) = grid%nG2A(local_id)
      grid_pf_active( grid%nG2A(local_id) ) = 1
       enddo
       grid_pf_npts_local   = grid%nlmax
       grid_pf_npts_ghost   = grid%ngmax - grid%nlmax
       grid_pf_npts_ghosted = grid_pf_npts_local + grid_pf_npts_ghost

       call pflotranModelFindOverlapCells(pflotran_model%pf_cells, grid_pf_npts_ghosted,&
            grid_pf_cell_ids_ghosted_nindex,pf_ocell_ids,pf_ocell_vol,pf_ocell_count,&
            pf_ocell_count_array)

       call MappingSetNumCells(m_clm2pf, grid_clm_npts_local, grid_clm_npts_ghost, &
          grid_clm_npts_local+grid_clm_npts_ghost, pf_ocell_count)
       call MappingAllocateMemory(m_clm2pf)
       call MappingSetOcells(m_clm2pf,pf_ocell_ids,pf_ocell_vol)
       call MappingSetTmeshCellIds( m_clm2pf, grid_clm_cell_ids_ghosted_nindex)


     call MappingSetNumCells(m_pf2clm, grid_pf_npts_local, grid_pf_npts_ghost, &
          grid%ngmax,clm_ocell_count)     
       call MappingAllocateMemory(m_pf2clm)
       call MappingSetOcells(m_pf2clm,clm_ocell_ids,clm_ocell_vol)
       call MappingSetTmeshCellIds( m_pf2clm,grid_pf_cell_ids_ghosted_nindex)

       ! Error checking to ensure all CLM cells that are in the 
     ! mapping file are active on some processor
       do local_id = 1,pflotran_model%num_clm_cells
       if ( grid_clm_active(local_id).eq.0 ) then
       write(*,*), 'All CLM grid '
           pflotran_model%option%io_buffer = 'All CLM cells present within mapping file ' //&
       ' are not active!'
          call printErrMsg(option)
     endif
     enddo

       ! Error checking to ensure all PF cells that are in the 
     ! mapping file are active on some processor
       do local_id = 1,pflotran_model%num_pf_cells
          if ( grid_pf_active(local_id).eq.0 ) then
             write(*,*), 'All CLM grid '
             pflotran_model%option%io_buffer = 'All PF cells present within mapping file ' //&
                  ' are not active!'
             call printErrMsg(option)
          endif
       enddo
     
       ! Free memory
     deallocate(clm_ocell_count_array)
       deallocate(clm_ocell_ids)
       deallocate(clm_ocell_vol)
     deallocate(pf_ocell_count_array)
       deallocate(pf_ocell_ids)
       deallocate(pf_ocell_vol)

   else
  
       grid_clm_npts_ghost = 0
       grid_clm_npts_ghosted = grid_clm_npts_local + grid_clm_npts_ghost

       allocate(clm_ocell_count_array(grid_clm_npts_local+grid_clm_npts_ghost))
       call MPI_Send(grid_clm_npts_local+grid_clm_npts_ghost,1,MPI_INTEGER,&
            option%io_rank,option%myrank,option%mycomm,ierr)
       call MPI_Send(grid_clm_cell_ids_ghosted_nindex,&
            grid_clm_npts_local+grid_clm_npts_ghost, &
            MPI_INTEGER,option%io_rank,option%myrank,option%mycomm,ierr)

       call MPI_Recv(clm_ocell_count, 1,MPI_INTEGER,0,MPI_ANY_TAG,option%mycomm,&
            status_mpi,ierr)
       allocate(clm_ocell_ids(clm_ocell_count))
       allocate(clm_ocell_vol(clm_ocell_count))
       call MPI_Recv(clm_ocell_ids,clm_ocell_count,MPI_INTEGER,0,MPI_ANY_TAG,option%mycomm, &
            status_mpi,ierr)
       call MPI_Recv(clm_ocell_vol,clm_ocell_count,MPI_DOUBLE_PRECISION,0,MPI_ANY_TAG,option%mycomm, &
            status_mpi,ierr)

       ! Allocate memory
       allocate(grid_pf_cell_ids_ghosted_nindex(grid%ngmax))
       allocate(pf_ocell_count_array(            grid%ngmax))

       do local_id = 1,grid%ngmax
          grid_pf_cell_ids_ghosted_nindex(local_id) = grid%nG2A(local_id)
       enddo
       grid_pf_npts_local = grid%nlmax
       grid_pf_npts_ghost = grid%ngmax - grid%nlmax

       call MPI_Send(grid%ngmax,1,MPI_INTEGER,option%io_rank,option%myrank,option%mycomm,&
            ierr)
       call MPI_Send(grid_pf_cell_ids_ghosted_nindex,grid%ngmax,MPI_INTEGER, &
            option%io_rank,option%myrank,option%mycomm,ierr)

       call MPI_Recv(pf_ocell_count,1,MPI_INTEGER,0,MPI_ANY_TAG,option%mycomm,&
            status_mpi,ierr)

       allocate(pf_ocell_ids(pf_ocell_count))
       allocate(pf_ocell_vol(pf_ocell_count))

       call MPI_Recv(pf_ocell_ids,pf_ocell_count,MPI_INTEGER,0,MPI_ANY_TAG,option%mycomm, &
            status_mpi,ierr)
       call MPI_Recv(pf_ocell_vol,pf_ocell_count,MPI_DOUBLE_PRECISION,0,MPI_ANY_TAG,option%mycomm, &
            status_mpi,ierr)

       call MappingSetNumCells(m_clm2pf, grid_clm_npts_local, grid_clm_npts_ghost, &
          grid_clm_npts_local+grid_clm_npts_ghost,pf_ocell_count)
       call MappingAllocateMemory(m_clm2pf)
       call MappingSetOcells(m_clm2pf,pf_ocell_ids,pf_ocell_vol)
       call MappingSetTmeshCellIds( m_clm2pf, grid_clm_cell_ids_ghosted_nindex)

     call MappingSetNumCells(m_pf2clm, grid_pf_npts_local, grid_pf_npts_ghost, &
          grid%ngmax,clm_ocell_count)     
       call MappingAllocateMemory(m_pf2clm)
       call MappingSetOcells(m_pf2clm,clm_ocell_ids,clm_ocell_vol)
       call MappingSetTmeshCellIds( m_pf2clm,grid_pf_cell_ids_ghosted_nindex)

     write(*,*), 'm_clm2pf: '
     write(*,*), '          grid_clm_npts_local   = ',grid_clm_npts_local
     write(*,*), '          grid_clm_npts_ghost   = ',grid_clm_npts_ghost
     write(*,*), '          grid_clm_npts_ghosted = ',grid_clm_npts_local+grid_clm_npts_ghost
     write(*,*), '          num_ocells_with_fmesh = ',pf_ocell_count
       write(*,*), '          num_docells_with_fmesh= ',m_clm2pf%num_docells_with_fmesh     
     write(*,*), 'm_pf2clm: '
     write(*,*), '          grid_clm_npts_local   = ',grid_pf_npts_local
     write(*,*), '          grid_clm_npts_ghost   = ',grid_pf_npts_ghost
     write(*,*), '          grid_clm_npts_ghosted = ',grid_pf_npts_ghosted
     write(*,*), '          num_ocells_with_fmesh = ',clm_ocell_count
       write(*,*), '          num_docells_with_fmesh= ',m_pf2clm%num_docells_with_fmesh     

       ! Free up memory
       deallocate(clm_ocell_ids)
       deallocate(clm_ocell_vol)
       deallocate(pf_ocell_ids)
       deallocate(pf_ocell_vol)

   endif

    call MPI_Barrier(pflotran_model%option%global_comm,ierr)

    ! Creating Vectors
    clm_pf_data%clm_num_local    = grid_clm_npts_local
    clm_pf_data%pf_num_local     = grid_pf_npts_local
    clm_pf_data%clm_num_docells  = m_pf2clm%num_docells_with_fmesh
    clm_pf_data%pf_num_docells   = m_clm2pf%num_docells_with_fmesh

    !write(*,*), 'clm_pf_data: ',option%myrank, clm_pf_data%clm_num_local, &
  !     clm_pf_data%pf_num_local, clm_pf_data%clm_num_docells, &
  !   clm_pf_data%pf_num_docells

    call clm_pf_data_allocate_memory(clm_pf_data, option%mycomm)


    ! Creating Index Sets
    call MappingCreateIS(m_clm2pf, option%mycomm, grid_clm_npts_local, &
         grid_clm_cell_ids_ghosted_nindex)
    call MappingCreateIS(m_pf2clm, option%mycomm, grid_pf_npts_local, &
         grid_pf_cell_ids_ghosted_nindex)

  if(option%myrank.eq.0) then
       write(*,*), 'm_pf2clm: '
     write(*,*), '         grid_pf_npts_local = ',grid_pf_npts_local
  endif   
     
    ! Creating VecScatter
    call MappingCreateVecScatter(m_clm2pf, option%mycomm, grid_clm_npts_local, &
         clm_pf_data%hksat_x_clmloc, clm_pf_data%hksat_x, clm_pf_data%hksat_x_pfloc)

    !call MappingCreateVecScatter(m_pf2clm, option%mycomm, grid_pf_npts_local, &
    !     clm_pf_data%sat_pfloc, clm_pf_data%sat, clm_pf_data%sat_clmloc)
#endif
  end subroutine pflotranModelInitMapping2



  ! ************************************************************************** !
  !
  ! pflotranModelInitMapping3:
  !
  !
  ! author: Gautam Bisht
  ! date: 03/24/2011
  ! ************************************************************************** !
  subroutine pflotranModelInitMapping3(pflotran_model, filename, &
       grid_clm_cell_ids_ghosted_nindex, grid_clm_npts_local, map_id,source_id)

    use Input_module
    use Option_module
    use Realization_module
    use Grid_module
    use Patch_module
    use clm_pflotran_interface_data

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    !type(clm_pflotran_data), intent(inout),pointer    :: clm_pf_data
    PetscInt, intent(in),pointer                      :: grid_clm_cell_ids_ghosted_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id, source_id
    character(len=MAXSTRINGLENGTH)                    :: filename
    
    ! local
    PetscInt                           :: local_id, grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_clm_npts_ghost
    PetscInt, pointer                  :: grid_pf_cell_ids_ghosted_nindex(:)
    PetscInt, pointer                  :: grid_pf_cell_ids_local_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_or_ghost_nindex(:)
    PetscInt                           :: count
    PetscErrorCode                     :: ierr
    type(mapping_type),pointer         :: m_clm2pf
    type(mapping_type),pointer         :: m_pf2clm
    type(mapping_type),pointer         :: map

    type(option_type), pointer         :: option
    type(realization_type), pointer    :: realization
    type(grid_type),pointer            :: grid
    type(patch_type), pointer          :: patch

    !write(*,*), 'In pflotranModelInitMapping3'
    option          => pflotran_model%option
    realization     => pflotran_model%simulation%realization
    patch           => realization%patch
    patch           => realization%patch
    grid            => patch%grid

    m_clm2pf        => pflotran_model%map_clm2pf
    m_pf2clm        => pflotran_model%map_pf2clm
    select case(map_id)
      case(1)
        map => pflotran_model%map_clm2pf
      case(2)
        map => pflotran_model%map_clm2pf_soils
      case(3)
        map => pflotran_model%map_pf2clm
      case default
        option%io_buffer = 'Invalid map_id argument to pflotranModelInitMapping3'
        call printErrMsg(option)
    end select

    allocate(grid_pf_cell_ids_ghosted_nindex(grid%ngmax))
    allocate(grid_pf_cell_ids_local_nindex  (grid%nlmax))
    allocate(grid_pf_local_or_ghost_nindex  (grid%ngmax))

    do local_id = 1,grid%nlmax
       grid_pf_cell_ids_local_nindex(local_id)   = grid%nL2A(local_id)
    enddo
    do local_id = 1,grid%ngmax
      grid_pf_cell_ids_ghosted_nindex(local_id) = grid%nG2A(local_id)-1
      if (grid%nG2L(local_id).eq.0) then
        grid_pf_local_or_ghost_nindex(local_id) = 0 ! GHOST
      else
        grid_pf_local_or_ghost_nindex(local_id) = 1 ! LOCAL
      endif
    enddo

    grid_pf_npts_local = grid%nlmax
    grid_pf_npts_ghost = grid%ngmax - grid%nlmax
    grid_clm_npts_ghost= 0
    
    select case(source_id)
      case(1)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_cell_ids_ghosted_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_ghosted_nindex, &
                                              grid_pf_local_or_ghost_nindex)
      case(2)
        call MappingSetSourceMeshCellIds(map, option, grid_pf_npts_local, &
          grid_pf_cell_ids_ghosted_nindex)
        option%io_buffer = 'Source_id = 2 in pflotranModelInitMapping3 not completed'
        call printErrMsg(option)
        !call MappingSetDestinationMeshCellIds(map, option, grid_clm_npts_local, &
        !  grid_clm_npts_ghost, grid_clm_cell_ids_ghosted_nindex, grid_clm_local_or_ghost_nindex)
      case default
        option%io_buffer = 'Invalid argument source_id passed to pflotranModelInitMapping3'
        call printErrMsg(option)
    end select
    
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MappingReadTxtFileMPI(map, trim(filename), option)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MappingFindDistinctSourceMeshCellIds(map,option)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MappingCreateWeightMatrix(map,option)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MappingCreateScatterOfSourceMesh(map, option)        
    if(option%myrank == PRINT_RANK) write(*,*),' ... done'

    deallocate(grid_pf_cell_ids_ghosted_nindex)
    deallocate(grid_pf_cell_ids_local_nindex)
    deallocate(grid_pf_local_or_ghost_nindex)

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

    input => InputCreate(fileid,filename)
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
    use Grid_module
    use Patch_module
    use Field_module
    use Global_Aux_module

    implicit none

    type(pflotran_model_type), pointer :: pflotran_model
    type(realization_type), pointer    :: realization
    type(grid_type),pointer            :: grid
    type(field_type),pointer           :: field
    type(global_auxvar_type), pointer  :: global_aux_vars(:)
    type(patch_type), pointer          :: patch

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

#ifdef CLM_PFLOTRAN
    write(iulog,*), '>>>> Inserting waypoint at pause_time = ',pause_time
#else
    !print *, '>>>> Inserting waypoint at pause_time = ',pause_time
#endif

    call pflotranModelInsertWaypoint(pflotran_model, pause_time)
    call pflotranModelInsertWaypoint(pflotran_model, pause_time + 100.0d0)

    call RichardsUpdateAuxVars(pflotran_model%simulation%realization)

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

    call StepperRunOneDT(pflotran_model%simulation%realization, &
         pflotran_model%simulation%flow_stepper, &
         pflotran_model%simulation%tran_stepper, pause_time)

#ifdef CLM_PFLOTRAN
    call RichardsUpdateAuxVars(pflotran_model%simulation%realization)
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

       source_sink   = source_sink + pf_clm_data%qsrc_flx(local_id)
    enddo

    write(iulog, *), '===================================================='
    write(iulog, *), 'Volume of LIQUID water:'
    write(iulog, *), 'Before PFLOTRAN call: ',liq_vol_start,                        '[m^3/area]'
    write(iulog, *), 'After  PFLOTRAN call: ',liq_vol_end  ,                        '[m^3/area]'
    write(iulog, *), 'Change              : ',liq_vol_start-liq_vol_end,            '[m^3/area]'
    write(iulog, *), 'Rate of change      : ',(liq_vol_end-liq_vol_start)/1800.0d0, '[m^3/area/sec]'
    write(iulog, *), 'Source_sink         : ',source_sink/998.2d0,                  '[m/sec]'
    write(iulog, *), '====================================================?'

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
  ! date: 10/28/2010
  ! ************************************************************************** !
#ifdef CLM_PFLOTRAN
  subroutine pflotranModelUpdateSourceSink (pflotran_model)

    use Realization_module
    use Patch_module
    use Grid_module
    use clm_pflotran_interface_type
    use Richards_Aux_module
    use Field_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    type(realization_type),pointer            :: realization
    type(patch_type),pointer                  :: patch
    type(grid_type),pointer                   :: grid
    type(field_type),pointer                  :: field
    type(mapping_type),pointer                :: map
    type(inside_each_clm_cell), pointer       :: clm_cell


    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id, pf_id
    PetscInt           :: i
    PetscReal          :: vol_ovlap, tmp

    realization     => pflotran_model%simulation%realization
    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field
    map             => pflotran_model%mapping


    ! Initialize source/sink to ZERO
    do local_id = 1,grid%nlmax

       ghosted_id = grid%nL2G(local_id)

       if (associated(patch%imat)) then
          if(patch%imat(ghosted_id) <= 0) cycle
       endif

       pf_clm_data%qsrc_flx(local_id) = 0.0d0

    enddo

    tmp = 0.0d0
    do i = 1,clm_pf_data%nlevsoi

       tmp = tmp + clm_pf_data%qflx_sink(1,i)
       do local_id = 1,map%clm2pf(i)%num_pflotran_cells

          pf_id     = map%clm2pf(i)%id_pflotran_cells(local_id)
          vol_ovlap = map%clm2pf(i)%perc_vol_overlap( local_id)

          pf_clm_data%qsrc_flx(pf_id) = pf_clm_data%qsrc_flx(pf_id) &
               - vol_ovlap * clm_pf_data%qflx_sink(1,i)
       enddo
    enddo

    write(iulog,*), 'pflotranModelUpdateSourceSink: tmp       [mm/sec] = ',tmp*1000.0d0 /998.2d0
    write(iulog,*), 'pflotranModelUpdateSourceSink: tmp*30*60 [mm    ] = ',tmp*30*60*1000.0d0 /998.2d0

  end subroutine pflotranModelUpdateSourceSink
#endif

  ! ************************************************************************** !
  !
  ! pflotranModelUpdateSaturation:
  !
  !
  ! author: Gautam Bisht
  ! date: 10/12/2010
  ! ************************************************************************** !
#ifdef CLM_PFLOTRAN
  subroutine pflotranModelUpdateSaturation( pflotran_model )

    use Realization_module
    use Patch_module
    use Grid_module
    use Field_module
    use Richards_Aux_module
    use Global_Aux_module
    use Option_module
    use Richards_module
    use Saturation_Function_module
    use Discretization_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    type(option_type), pointer                :: option
    type(realization_type),pointer            :: realization
    type(patch_type),pointer                  :: patch
    type(grid_type),pointer                   :: grid
    type(field_type),pointer                  :: field
    type(richards_auxvar_type), pointer       :: rich_aux_vars(:)
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    type(inside_each_clm_cell), pointer       :: clm_cell
    type(mapping_type),pointer                :: map
    type(inside_each_pflotran_cell), pointer  :: pf_cell
    type(saturation_function_type)            :: sat_func
    type(richards_auxvar_type)                :: aux_var
    type(global_auxvar_type)                  :: global_aux_var
    type(discretization_type), pointer        :: discretization

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id, clm_id, i
    PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:)
    PetscReal          :: sat_old, sat_new, sat
    PetscReal          :: pc_old, pc_new, pc
    PetscReal          :: lambda, alpha, Sr, vol_ovlap

    realization      => pflotran_model%simulation%realization
    discretization   => realization%discretization
    patch            => realization%patch
    grid             => patch%grid
    field            => realization%field
    rich_aux_vars    => patch%aux%Richards%aux_vars
    global_aux_vars  => patch%aux%Global%aux_vars
    option           => realization%option
    map              => pflotran_model%mapping

    write (iulog,*), 'In pflotranModelUpdateSaturation'

    ! Initialize sat_new to ZERO
    do local_id = 1,grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if(patch%imat(ghosted_id) <= 0) cycle
       endif
       pf_clm_data%sat_new(local_id) = 0.0d0
    enddo

    do local_id = 1,grid%nlmax

       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if(patch%imat(ghosted_id) <= 0) cycle
       endif

       pf_cell => map%pf2clm(local_id)

       do i = 1,pf_cell%num_clm_cells

          clm_id    = pf_cell%id_clm_cells(i)
          vol_ovlap = pf_cell%perc_vol_overlap(i)

          pf_clm_data%sat_new(local_id) = pf_clm_data%sat_new(local_id) +&
               vol_ovlap * clm_pf_data%sat(1,clm_id)

       enddo
    enddo

    call GridVecGetArrayF90(grid,field%flow_xx,  xx_loc_p,   ierr)
    call GridVecGetArrayF90(grid,field%icap_loc, icap_loc_p, ierr)

    do local_id = 1, grid%nlmax

       ghosted_id = grid%nL2G(local_id)

       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif

       if (map%pf2clm(local_id)%num_clm_cells.gt.0) then
          sat_new = pf_clm_data%sat_new(ghosted_id) / map%pf2clm(local_id)%total_vol_overlap

          sat_func = realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr
          aux_var        = rich_aux_vars(ghosted_id)
          global_aux_var = global_aux_vars(ghosted_id)
          lambda         = aux_var%bc_lambda
          alpha          = aux_var%bc_alpha
          Sr             = sat_func%Sr(1)

          pc = option%reference_pressure - xx_loc_p(ghosted_id)

          sat_old = Sr + (1.0d0 - Sr)* ((pc*alpha)**(-lambda))
          pc_old  = (1.0d0/alpha) * (((1.0d0 - Sr)/(sat_old - Sr))**(1.0/lambda))

          !sat_new = 0.99d0*sat_old
          sat_new = min(sat_new, 1.0d0)
          if ( sat_new.lt.Sr ) then
             write(iulog,*), 'sat_new < Sr'
             pc_new = 1.0d0/alpha
          else
             pc_new  = (1.0d0/alpha) * (((1.0d0 - Sr)/(sat_new - Sr))**(1.0/lambda))
          endif

          !write(iulog, *), 'sat diff: ',sat_new, sat_old
          !write(iulog, *), 'pc  diff: ',pc_new ,  pc_old

          xx_loc_p(ghosted_id) = option%reference_pressure - pc_new

       endif

    enddo

    call GridVecRestoreArrayF90(grid,field%flow_xx,  xx_loc_p,   ierr)
    call GridVecRestoreArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
    ! update dependent vectors
    call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
         field%flow_xx_loc,NFLOWDOF)
    call VecCopy(field%flow_xx, field%flow_yy, ierr)


    local_id = grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    write (iulog, *), 'SAT:', global_aux_vars(ghosted_id)%sat(1)

    call RichardsUpdateAuxVars(realization)

    local_id = grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    write (iulog, *), 'SAT:', global_aux_vars(ghosted_id)%sat(1)

  end subroutine pflotranModelUpdateSaturation
#endif

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
                dataset => sub_condition%dataset
                dataset%values(1,1) = flux_value
                print *,'In UpdateTopBCHomogeneous: ', flux_value
                !call FlowSubConditionUpdateDataset(option,1.0d0,sub_condition%dataset)
             endif

          enddo
       endif

       condition => condition%next

    enddo

    call StepperUpdateSolution(realization)

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

#ifdef CLM_PFLOTRAN
    type(mapping_type),pointer         :: map

    map         => pflotran_model%mapping

    deallocate ( pf_clm_data%zwt      )
    deallocate ( pf_clm_data%alpha    )
    deallocate ( pf_clm_data%lambda   )
    deallocate ( pf_clm_data%qsrc_flx )
    deallocate ( pf_clm_data%sat_new  )
    deallocate ( map%pf2clm )
    deallocate ( map%clm2pf )
#endif

    ! Clean things up.
    call SimulationDestroy(pflotran_model%simulation)

    ! Final Time
    call PetscGetCPUTime(timex(2), ierr)
    call PetscGetTime(timex_wall(2), ierr)

    if (pflotran_model%option%myrank == pflotran_model%option%io_rank) then

       if (pflotran_model%option%print_to_screen) then
#ifdef CLM_PFLOTRAN
          write(iulog,'(/," CPU Time:", 1pe12.4, " [sec] ", &
               & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
               timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
               (timex(2)-timex(1))/3600.d0

          write(iulog,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
               & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
               timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
               (timex_wall(2)-timex_wall(1))/3600.d0
#else
          write(*,'(/," CPU Time:", 1pe12.4, " [sec] ", &
               & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
               timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
               (timex(2)-timex(1))/3600.d0

          write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
               & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
               timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
               (timex_wall(2)-timex_wall(1))/3600.d0
#endif
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
#ifndef CLM_PFLOTRAN
    call MPI_Finalize(ierr)
#endif

  end subroutine pflotranModelDestroy

end module pflotran_model_module
