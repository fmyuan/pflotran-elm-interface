module pflotran_clm_setmapping_module

  use Option_module, only : option_type
  use Simulation_Base_class, only : simulation_base_type
  use Multi_Simulation_module, only : multi_simulation_type
  use Realization_Base_class, only : realization_base_type
  use Mapping_module, only : mapping_type
  use pflotran_clm_main_module, only : pflotran_model_type

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petsclog.h"
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscvec.h"

  ! Note:
  !
  ! CLM has the following:
  !   (i) 3D subsurface grid (CLM_SUB);
  !   (ii) 2D top-cell grid (CLM_2DTOP).
  !   (iii) 2D bottom-cell grid (CLM_2DBOT)
  ! CLM decomposes the 3D subsurface grid across processors in a 2D (i.e.
  ! cells in Z are not split across processors). Thus, the surface cells of
  ! 3D subsurface grid are on the same processors as the 2D surface grid.
  !
  ! PFLOTRAN has the following:
  !   (i) 3D subsurface grid (PF_SUB);
  !   (ii) top-face control volumes of 3D subsurface grid (PF_2DTOP);
  !   (iii) bottom face control volumes of 3D subsurface grid (PF_2DBOT);
  ! In PFLOTRAN, control volumes in PF_2DSUB and PF_SRF may reside on different
  ! processors. PF_SUB and PF_2DSUB are derived from simulation%realization;
  ! while PF_SRF refers to simulation%surf_realization.

  ! map level constants
  PetscInt, parameter, public :: CLM_SUB_TO_PF_SUB           = 1 ! 3D --> 3D
  PetscInt, parameter, public :: PF_SUB_TO_CLM_SUB           = 2 ! 3D --> 3D

  PetscInt, parameter, public :: CLM_2DTOP_TO_PF_2DTOP       = 3 ! TOP face of 3D cell
  PetscInt, parameter, public :: PF_2DTOP_TO_CLM_2DTOP       = 4 ! TOP face of 3D cell
  PetscInt, parameter, public :: CLM_2DBOT_TO_PF_2DBOT       = 5 ! BOTTOM face of 3D cell
  PetscInt, parameter, public :: PF_2DBOT_TO_CLM_2DBOT       = 6 ! BOTTOM face of 3D cell

  ! mesh ids
  PetscInt, parameter, public :: CLM_SUB_MESH   = 1
  PetscInt, parameter, public :: PF_SUB_MESH    = 2
  PetscInt, parameter, public :: CLM_FACE_MESH  = 3
  PetscInt, parameter, public :: PF_FACE_MESH   = 4

  PetscInt, parameter, public :: CLM_2DTOP_MESH = 5
  PetscInt, parameter, public :: PF_2DTOP_MESH  = 6

!  type, public :: inside_each_overlapped_cell
!     PetscInt           :: id
!     PetscInt           :: ocell_count
!     PetscInt,  pointer :: ocell_id(:)
!     PetscReal, pointer :: perc_vol_overlap(:)
!     PetscReal          :: total_vol_overlap
!  end type inside_each_overlapped_cell

  public ::                                  &
       ! mesh-mapping
       pflotranModelSetupMappingFiles,       &
       pflotranModelInitMapping,             &
       pflotranModelGetTopFaceArea

  private ::                                 &
       pflotranModelInitMappingSub2Sub,      &
       pflotranModelInitMapTopTo2DSub,       &
       pflotranModelNSurfCells3DDomain,      &
       pflotranModelInitMapFaceToFace

contains

! ************************************************************************** !

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
  ! Revised by Fengming YUAN

    use String_module
    use Option_module
    use Input_Aux_module
    use Mapping_module

    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Realization_class, only : realization_type

    implicit none

#include "petsc/finclude/petscsys.h"

    type(pflotran_model_type), pointer, intent(inout) :: model
    type(input_type), pointer :: input
    type(option_type), pointer :: option
    class(realization_type), pointer    :: realization

    PetscBool :: clm2pf_3dsub_file
    PetscBool :: pf2clm_3dsub_file

    PetscBool :: clm2pf_bctop_file
    PetscBool :: pf2clm_bctop_file

    PetscBool :: clm2pf_bcbot_file
    PetscBool :: pf2clm_bcbot_file

    character(len=MAXSTRINGLENGTH) :: string
    character(len=MAXWORDLENGTH) :: word

    !nullify(model%pf_cells)
    !nullify(model%clm_cells)
    nullify(model%map_clm_sub_to_pf_sub)
    nullify(model%map_clm_2dtop_to_pf_2dtop)
    nullify(model%map_pf_sub_to_clm_sub)
    !
    nullify(model%map_clm_2dbot_to_pf_2dbot)
    nullify(model%map_pf_2dbot_to_clm_2dbot)
    nullify(model%map_pf_2dtop_to_clm_2dtop)

    model%map_clm_sub_to_pf_sub            => MappingCreate()
    model%map_clm_2dtop_to_pf_2dtop        => MappingCreate()
    model%map_pf_sub_to_clm_sub            => MappingCreate()
    model%map_pf_2dtop_to_clm_2dtop        => MappingCreate()
    model%map_clm_2dbot_to_pf_2dbot        => MappingCreate()
    model%map_pf_2dbot_to_clm_2dbot        => MappingCreate()

    model%nlclm = -1
    model%ngclm = -1

    input => InputCreate(IUNIT_TEMP, &
                    model%option%input_filename, model%option)

    ! Read names of mapping file
    clm2pf_3dsub_file=PETSC_FALSE
    pf2clm_3dsub_file=PETSC_FALSE
    clm2pf_bctop_file=PETSC_FALSE
    pf2clm_bctop_file=PETSC_FALSE
    clm2pf_bcbot_file=PETSC_FALSE
    pf2clm_bcbot_file=PETSC_FALSE
    
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
        case('CLM2PF_SUB_FILE')
          call InputReadNChars(input, model%option, model%map_clm_sub_to_pf_sub%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_clm_sub_to_pf_sub%filename = &
            trim(model%map_clm_sub_to_pf_sub%filename)//CHAR(0)
          model%map_clm_sub_to_pf_sub%id = ONE_INTEGER
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')   
          clm2pf_3dsub_file=PETSC_TRUE
        case('PF2CLM_SUB_FILE')
          call InputReadNChars(input, model%option, model%map_pf_sub_to_clm_sub%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_pf_sub_to_clm_sub%filename = &
            trim(model%map_pf_sub_to_clm_sub%filename)//CHAR(0)
          model%map_pf_sub_to_clm_sub%id = TWO_INTEGER
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          pf2clm_3dsub_file=PETSC_TRUE
        case('CLM2PF_BCTOP_FILE')
          call InputReadNChars(input, model%option, model%map_clm_2dtop_to_pf_2dtop%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_clm_2dtop_to_pf_2dtop%filename = &
            trim(model%map_clm_2dtop_to_pf_2dtop%filename)//CHAR(0)
          model%map_clm_2dtop_to_pf_2dtop%id = THREE_INTEGER
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          clm2pf_bctop_file=PETSC_TRUE
        case('PF2CLM_BCTOP_FILE')
          call InputReadNChars(input, model%option, model%map_pf_2dtop_to_clm_2dtop%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_pf_2dtop_to_clm_2dtop%filename = &
              trim(model%map_pf_2dtop_to_clm_2dtop%filename)//CHAR(0)
          model%map_pf_2dtop_to_clm_2dtop%id = FOUR_INTEGER
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          pf2clm_bctop_file=PETSC_TRUE
        case('CLM2PF_BCBOT_FILE')
          call InputReadNChars(input, model%option, model%map_clm_2dbot_to_pf_2dbot%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_clm_2dbot_to_pf_2dbot%filename = &
            trim(model%map_clm_2dbot_to_pf_2dbot%filename)//CHAR(0)
          model%map_clm_2dbot_to_pf_2dbot%id = FIVE_INTEGER
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          clm2pf_bcbot_file=PETSC_TRUE
        case('PF2CLM_BCBOT_FILE')
          call InputReadNChars(input, model%option, model%map_pf_2dbot_to_clm_2dbot%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_pf_2dbot_to_clm_2dbot%filename = &
            trim(model%map_pf_2dbot_to_clm_2dbot%filename)//CHAR(0)
          model%map_pf_2dbot_to_clm_2dbot%id = SIX_INTEGER
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          pf2clm_bcbot_file=PETSC_TRUE
        case default
          model%option%io_buffer='Keyword ' // trim(word) // &
            ' in input file not recognized and ignored'
          call printMsg(model%option)
      end select

    enddo
    call InputDestroy(input)

    if ((.not. clm2pf_3dsub_file) .or. &
        (.not. pf2clm_3dsub_file) ) then
      model%option%io_buffer='One of the 3D soil-mesh mapping files not found'
      call printErrMsg(model%option)
    endif
    
    if(   (model%option%iflowmode==TH_MODE .or. model%option%iflowmode==RICHARDS_MODE)  &
     .and.(.not.clm2pf_bctop_file .or. .not.pf2clm_bctop_file .or. &
           .not.clm2pf_bcbot_file .or. .not.pf2clm_bcbot_file) ) then
      model%option%io_buffer='Running in TH_MODE/Richards_MODE without one of 4 top/bottom-cell mesh files'
      call printErrMsg(model%option)
    endif

  end subroutine pflotranModelSetupMappingFiles

! ************************************************************************** !

  subroutine pflotranModelInitMapping(pflotran_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id)
  ! 
  ! Initialize mapping between the two model grid
  ! (CLM and PFLTORAN)
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/24/2011
  ! Revised by Fengming YUAN

    !use Input_Aux_module
    use Option_module
    !use Realization_class
    !use Grid_module
    !use Patch_module
    !use Coupler_module
    !use Connection_module
    use String_module
    !use Simulation_Base_class, only : simulation_base_type
    !use Simulation_Subsurface_class, only : subsurface_simulation_type
    !use Simulation_Surface_class, only : surface_simulation_type
    !use Simulation_Surf_subsurf_class, only : surfsubsurface_simulation_type
    use Mapping_module

    implicit none

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename
    
    select case (map_id)
      case (CLM_SUB_TO_PF_SUB, PF_SUB_TO_CLM_SUB)
        call pflotranModelInitMappingSub2Sub(pflotran_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id)
      ! A more generalized Mapping for Faces (sidesets) is now implemented (F.-M. Yuan)
      !case (CLM_2DTOP_TO_PF_2DTOP, PF_2DTOP_TO_CLM_2DTOP)
      !  call pflotranModelInitMapTopTo2DSub(pflotran_model,  &
      !                                      grid_clm_cell_ids_nindex, &
      !                                      grid_clm_npts_local, &
      !                                      map_id)
      case (CLM_2DTOP_TO_PF_2DTOP, PF_2DTOP_TO_CLM_2DTOP, &
            CLM_2DBOT_TO_PF_2DBOT, PF_2DBOT_TO_CLM_2DBOT)
        call pflotranModelInitMapFaceToFace(pflotran_model,  &
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
  ! Revised by Fengming YUAN

    !use Input_Aux_module
    use Option_module
    use Realization_class
    use Grid_module
    use Patch_module
    use String_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
    use Mapping_module
    use clm_pflotran_interface_data

    implicit none

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename

    ! local
    PetscInt                           :: local_id, ghosted_id, grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_clm_npts_ghost, source_mesh_id
    PetscInt                           :: dest_mesh_id
    PetscInt, pointer                  :: grid_pf_cell_ids_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_local_shuffle_nindex(:)
    PetscInt                           :: numg, numk, g, k
    PetscErrorCode:: ierr

    type(mapping_type), pointer        :: map
    type(option_type), pointer         :: option
    class(realization_type), pointer   :: realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch

!-----------------------------------------------------------------------------------

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
      case(PF_SUB_TO_CLM_SUB)
        map => pflotran_model%map_pf_sub_to_clm_sub
        source_mesh_id = PF_SUB_MESH
        dest_mesh_id = CLM_SUB_MESH
      case default
        option%io_buffer = 'Invalid map_id argument to pflotranModelInitMapping'
        call printErrMsg(option)
    end select

    ! Read mapping file
    if (index(map%filename, '.h5') > 0) then
      call MappingReadHDF5(map, map%filename, option)
    else
      call MappingReadTxtFile(map, map%filename, option)

      ! checking if the CLM-PF has same number of soil layers for mapping
      if (map%pflotran_nlev /= clm_pf_idata%nzclm_mapped) then
         option%io_buffer = 'Invalid mapping soil layers between CLM and PFLOTRAN!'
        call printErrMsg(option)
      endif

    endif

    grid_clm_npts_ghost=0

    ! Allocate memory to identify if CLM cells are local or ghosted.
    ! Note: Presently all CLM cells are local
    allocate(grid_clm_local_nindex(grid_clm_npts_local))
    allocate(grid_clm_local_shuffle_nindex(grid_clm_npts_local))
    do local_id = 1, grid_clm_npts_local
      grid_clm_local_nindex(local_id) = local_id     ! LOCAL ID

      ! the following is a hack for passing Data from CLMp vec -> PF array --> CLMs vec directly for 2-D domain
      ! it will twist horizontal and vertical axis of CLM domain.
      numk = clm_pf_idata%nzclm_mapped
      numg = grid_clm_npts_local/numk
      g = mod(local_id-1, numg) + 1
      k = (local_id - g )/numg + 1
      grid_clm_local_shuffle_nindex(local_id) = (g-1)*numk + k

     enddo

    ! Find cell IDs for PFLOTRAN grid
    grid_pf_npts_local = grid%nlmax
    grid_pf_npts_ghost = grid%ngmax - grid%nlmax

    allocate(grid_pf_cell_ids_nindex(grid%ngmax))
    allocate(grid_pf_local_nindex(grid%ngmax))
    do ghosted_id = 1, grid%ngmax
      local_id = grid%nG2L(ghosted_id)
      grid_pf_cell_ids_nindex(ghosted_id) = grid%nG2A(ghosted_id)-1
      if (local_id <= 0) then
        grid_pf_local_nindex(ghosted_id) = 0        ! GHOST
      else
        grid_pf_local_nindex(ghosted_id) = local_id        ! LOCAL ID

      endif
    enddo

    select case(source_mesh_id)
      case(CLM_SUB_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_npts_ghost, &
                                         grid_clm_cell_ids_nindex, &
                                         grid_clm_local_shuffle_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_SUB_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_pf_npts_local, &
                                        grid_pf_npts_ghost, &
                                        grid_pf_cell_ids_nindex, &
                                        grid_pf_local_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_clm_npts_local, &
                                              grid_clm_npts_ghost, &
                                              grid_clm_cell_ids_nindex, &
                                              grid_clm_local_shuffle_nindex)
      case default
        option%io_buffer = 'Invalid argument source_mesh_id passed to pflotranModelInitMapping'
        call printErrMsg(option)
    end select

    call MappingDecompose(map, option)
    call MappingFindDistinctSourceMeshCellIds(map, option)
    call MappingCreateWeightMatrix(map, option)
    call MappingCreateScatterOfSourceMesh(map, option)
    call MappingFreeNotNeeded(map)

    deallocate(grid_pf_cell_ids_nindex)
    deallocate(grid_pf_local_nindex)
    deallocate(grid_clm_local_nindex)
    deallocate(grid_clm_local_shuffle_nindex)

  end subroutine pflotranModelInitMappingSub2Sub

! ************************************************************************** !

  subroutine pflotranModelInitMapTopTo2DSub(pflotran_model,           &
                                            grid_clm_cell_ids_nindex, &
                                            grid_clm_npts_local,      &
                                            map_id)
  ! 
  ! This routine maps CLM top face grid onto surface of PFLOTRAN 3D subsurface
  ! grid.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/09/13
  ! Revised by Fengming Yuan, CCSI-ORNL

    !use Input_Aux_module
    use Option_module
    use Realization_class
    use Grid_module
    use Patch_module
    use Coupler_module
    use Connection_module
    use String_module
    use clm_pflotran_interface_data
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    !use Simulation_Surface_class, only : surface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
    use Mapping_module

    implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscviewer.h"

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
    class(realization_type), pointer   :: realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch
    type(coupler_type), pointer        :: boundary_condition
    !type(coupler_type), pointer        :: source_sink
    type(connection_set_type), pointer :: cur_connection_set

!-----------------------------------------------------------------------------

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
      case(CLM_2DTOP_TO_PF_2DTOP)
        map => pflotran_model%map_clm_2dtop_to_pf_2dtop
        source_mesh_id = CLM_2DTOP_MESH
        dest_mesh_id = PF_2DTOP_MESH
      case(PF_2DTOP_TO_CLM_2DTOP)
        map => pflotran_model%map_pf_2dtop_to_clm_2dtop
        source_mesh_id = PF_2DTOP_MESH
        dest_mesh_id = CLM_2DTOP_MESH
      case default
        option%io_buffer = 'Invalid map_id argument to ' // &
          'pflotranModelInitMappingSurfTo2DSub'
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
    clm_pf_idata%nlclm_2dtop = grid_clm_npts_local
    clm_pf_idata%ngclm_2dtop = grid_clm_npts_local

    ! Find cell IDs for PFLOTRAN 3D-grid's surface
    ! Mapping to/from surface of PFLOTRAN domain
    found=PETSC_FALSE
    grid_pf_npts_local = 0
    grid_pf_npts_ghost = 0

    !select case (dest_mesh_id)
    if (source_mesh_id == PF_2DTOP_MESH .or.  &    !  case(PF_2DSUB_MESH)
        dest_mesh_id == PF_2DTOP_MESH) then        ! mesh is PF_2DSUB_MESH

        patch => realization%patch
        grid => patch%grid

        boundary_condition => patch%boundary_condition_list%first
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
        clm_pf_idata%nlpf_2dtop  = grid_pf_npts_local
        clm_pf_idata%ngpf_2dtop  = grid_pf_npts_local

      else !case default
        option%io_buffer='Unknown mesh'
        call printErrMsg(option)
      
    endif !end select

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
    CHKERRQ(ierr)
    call VecCreateMPI(option%mycomm, grid%nlmax, PETSC_DECIDE, surf_ids, ierr)
    CHKERRQ(ierr)
    call VecSet(surf_ids, -1.d0, ierr)
    CHKERRQ(ierr)
    
    ! Set 1.0 to all cells that make up surface of PFLOTRAN subsurface domain
    call VecSetValues(surf_ids, grid_pf_npts_local, grid_pf_cell_ids_nindex, &
                      v_loc, INSERT_VALUES, ierr)
    CHKERRQ(ierr)
    deallocate(v_loc)
    call VecAssemblyBegin(surf_ids, ierr)
    CHKERRQ(ierr)
    call VecAssemblyEnd(surf_ids, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids, v_loc, ierr)
    CHKERRQ(ierr)
    count = 0
    do local_id=1,grid%nlmax
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr)
    CHKERRQ(ierr)

    count = 0
    do local_id=1,grid%nlmax
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif
    enddo
    call VecRestoreArrayF90(surf_ids, v_loc, ierr)
    CHKERRQ(ierr)

    !
    allocate(int_array(grid_pf_npts_local))
    do iconn = 1, grid_pf_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, grid_pf_npts_local, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    CHKERRQ(ierr)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm, grid_pf_npts_local, grid_pf_cell_ids_nindex, &
                         PETSC_COPY_VALUES, is_from, ierr)
    CHKERRQ(ierr)


    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_from, ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_to, ierr)
    CHKERRQ(ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    count = 0
    do iconn = 1, grid_pf_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_pf_cell_ids_nindex(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    call VecDestroy(surf_ids_loc, ierr)
    CHKERRQ(ierr)

    !
    ! Step-2: Recompute 'map%s2d_i/jscr' for pf mesh
    !
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_nwts, surf_ids_loc, ierr)
    CHKERRQ(ierr)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    CHKERRQ(ierr)


    do iconn = 1, map%s2d_nwts
      if (dest_mesh_id == PF_2DTOP_MESH) then
         int_array(iconn) = map%s2d_icsr(iconn)
      elseif (source_mesh_id == PF_2DTOP_MESH) then
         int_array(iconn) = map%s2d_jcsr(iconn)
      endif
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_from, ierr)
    CHKERRQ(ierr)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_from, ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_to, ierr)
    CHKERRQ(ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        if (dest_mesh_id == PF_2DTOP_MESH) then
            map%s2d_icsr(count) = INT(v_loc(iconn))
        elseif (source_mesh_id == PF_2DTOP_MESH) then
            map%s2d_jcsr(count) = INT(v_loc(iconn))
        endif
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    call VecDestroy(surf_ids_loc, ierr)
    CHKERRQ(ierr)
    
    if(count /= map%s2d_nwts) then
      option%io_buffer='No. of surface cells in mapping dataset does not ' // &
        'match surface cells on which BC is applied.'
      call printErrMsg(option)
    endif
    call VecDestroy(surf_ids, ierr)
    CHKERRQ(ierr)

    !
    ! Step-3: Find surface cells-ids of CLM subsurface domain
    !
    allocate(v_loc(grid_clm_npts_local))
    v_loc = 1.d0
    call VecCreateSeq(PETSC_COMM_SELF, grid_clm_npts_local, surf_ids_loc, ierr)
    CHKERRQ(ierr)
    call VecCreateMPI(option%mycomm, clm_pf_idata%nlclm_sub, PETSC_DECIDE, surf_ids, ierr)
    CHKERRQ(ierr)
    call VecSet(surf_ids, -1.d0, ierr)
    CHKERRQ(ierr)

    ! Set 1.0 to all cells that make up surface of CLM subsurface domain
    call VecSetValues(surf_ids, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                      v_loc, INSERT_VALUES, ierr)
    CHKERRQ(ierr)

    deallocate(v_loc)
    call VecAssemblyBegin(surf_ids, ierr)
    CHKERRQ(ierr)
    call VecAssemblyEnd(surf_ids, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids, v_loc, ierr)
    CHKERRQ(ierr)
    count = 0
    do local_id=1,clm_pf_idata%nlclm_sub
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr)
    CHKERRQ(ierr)

    count = 0
    do local_id=1,clm_pf_idata%nlclm_sub
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif
    enddo
    call VecRestoreArrayF90(surf_ids, v_loc, ierr)
    CHKERRQ(ierr)

    !
    allocate(int_array(grid_clm_npts_local))
    do iconn = 1, grid_clm_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, grid_clm_npts_local, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    CHKERRQ(ierr)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                         PETSC_COPY_VALUES, is_from, ierr)
    CHKERRQ(ierr)


    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_from, ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_to, ierr)
    CHKERRQ(ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    count = 0
    do iconn = 1, grid_clm_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_clm_cell_ids_nindex_copy(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    call VecDestroy(surf_ids_loc, ierr)
    CHKERRQ(ierr)

    !
    ! Step-4: Recompute 'map%s2d_i/jscr' for clm mesh
    !
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_nwts, surf_ids_loc, ierr)
    CHKERRQ(ierr)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    CHKERRQ(ierr)

    do iconn = 1, map%s2d_nwts
      if (source_mesh_id == CLM_SUB_MESH) then
         int_array(iconn) = map%s2d_jcsr(iconn)
      elseif (dest_mesh_id == CLM_SUB_MESH) then
         int_array(iconn) = map%s2d_icsr(iconn)
      endif
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_from, ierr)
    CHKERRQ(ierr)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(surf_ids, is_from, surf_ids_loc, is_to, vec_scat, &
                          ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_from, ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_to, ierr)
    CHKERRQ(ierr)

    call VecScatterBegin(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterEnd(vec_scat, surf_ids, surf_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        if (source_mesh_id == CLM_SUB_MESH) then
           map%s2d_jcsr(count) = INT(v_loc(iconn))
        elseif (dest_mesh_id == CLM_SUB_MESH) then
           map%s2d_icsr(count) = INT(v_loc(iconn))
        endif
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    call VecDestroy(surf_ids_loc, ierr)
    CHKERRQ(ierr)
    
    if(count /= map%s2d_nwts) then
      option%io_buffer='No. of surface cells in mapping dataset does not ' // &
        'match surface cells on which BC is applied.'
      call printErrMsg(option)
    endif
    call VecDestroy(surf_ids, ierr)
    CHKERRQ(ierr)

    select case(source_mesh_id)
      case(CLM_2DTOP_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_npts_ghost, &
                                         grid_clm_cell_ids_nindex_copy, &
                                         grid_clm_local_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_2DTOP_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_pf_npts_local, &
                                        grid_pf_npts_ghost, &
                                        grid_pf_cell_ids_nindex, &
                                        grid_pf_local_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_clm_npts_local, &
                                              grid_clm_npts_ghost, &
                                              grid_clm_cell_ids_nindex_copy, &
                                              grid_clm_local_nindex)
      case default
        option%io_buffer = 'Invalid argument source_mesh_id passed to ' // &
          'pflotranModelInitMappingSurfTo2DSub'
        call printErrMsg(option)
    end select

    deallocate(grid_pf_cell_ids_nindex)
    deallocate(grid_pf_local_nindex)
    deallocate(grid_clm_cell_ids_nindex_copy)
    deallocate(grid_clm_local_nindex)

    call MappingDecompose(map, option)
    call MappingFindDistinctSourceMeshCellIds(map, option)
    call MappingCreateWeightMatrix(map, option)
    call MappingCreateScatterOfSourceMesh(map, option)
    call MappingFreeNotNeeded(map)

  end subroutine pflotranModelInitMapTopTo2DSub

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
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
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

    coupler_list => realization%patch%boundary_condition_list
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
      call printMsg(pflotran_model%option, &
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
    use Grid_Unstructured_Aux_module
    use Grid_Unstructured_Cell_module
    use Grid_Unstructured_module
    use Grid_module
    use clm_pflotran_interface_data
    use Utility_module, only : DotProduct, CrossProduct
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    !use Simulation_Surface_class, only : surface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
    use Surface_Realization_class, only : surface_realization_type
    use Realization_class, only : realization_type
    use Mapping_module

    implicit none
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

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

    call VecGetArrayF90(clm_pf_idata%area_top_face_pfp, area_p, ierr)
    CHKERRQ(ierr)
    if(grid%itype == STRUCTURED_GRID) then
      ! Structured grid
      do ghosted_id=1,grid%ngmax
        local_id = grid%nG2L(ghosted_id)
        if(local_id>0) then
          area1 = grid%structured_grid%dx(ghosted_id)* &
                  grid%structured_grid%dy(ghosted_id)
          area_p(local_id) = area1
        endif
      enddo
    else if (grid%itype == UNSTRUCTURED_GRID) then
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
    call VecRestoreArrayF90(clm_pf_idata%area_top_face_pfp, area_p, ierr)
    CHKERRQ(ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%area_top_face_pfp, &
                                    clm_pf_idata%area_top_face_clms)

  end subroutine pflotranModelGetTopFaceArea

! ************************************************************************** !

  subroutine pflotranModelInitMapFaceToFace(pflotran_model,  &
                                            grid_clm_cell_ids_nindex, &
                                            grid_clm_npts_local, &
                                            map_id)
  !
  ! This routine maps CLM grids/columns structure onto BC faces
  ! (TOP, BOTTOM, EAST, WEST, NORTH, or, SOUTH, which type depends on BC condition-name in PF input cards)
  ! of PFLOTRAN 3D Domain grid,
  ! by extending GB's code - Fengming Yuan, ORNL
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/09/13
  !
  ! 02/14/2014 - TOP/BOTTOM faces, from CLM => PF, finished
  ! 05/12/2014 - TOP/BOTTOM faces, from PF => CLM, finished
  !
  ! NOTE: for TOP face, BC condition name: 'clm_gflux_bc') ('g' for ground);
  !       for BOTTOM face, BC condition name: 'clm_bflux_bc') ('b' for bottom);

    !use Input_Aux_module
    use Option_module
    use Realization_class
    use Grid_module
    use Patch_module
    use Coupler_module
    use Connection_module
    use String_module
    use clm_pflotran_interface_data
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    !use Simulation_Surface_class, only : surface_simulation_type
    use Simulation_Surf_Subsurf_class, only : surfsubsurface_simulation_type
    use Mapping_module

    implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscviewer.h"

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id

    ! local
    PetscInt                           :: local_id, grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_clm_npts_ghost, source_mesh_id
    PetscInt                           :: dest_mesh_id
    PetscInt, pointer                  :: grid_pf_cell_ids_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_cell_ids_nindex_copy(:)
    PetscInt                           :: count
    PetscInt                           :: ghosted_id
    PetscInt                           :: iconn
    PetscInt                           :: istart
    PetscInt, pointer                  :: int_array(:)
    PetscBool                          :: found
    PetscScalar,pointer                :: v_loc(:)
    PetscErrorCode                     :: ierr

    Vec                                :: face_ids
    Vec                                :: face_ids_loc
    IS                                 :: is_from
    IS                                 :: is_to
    VecScatter                         :: vec_scat

    type(mapping_type), pointer        :: map
    type(option_type), pointer         :: option
    class(realization_type), pointer   :: realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch
    type(coupler_type), pointer        :: boundary_condition
    type(coupler_type), pointer        :: source_sink
    type(connection_set_type), pointer :: cur_connection_set
    character(len=MAXSTRINGLENGTH)     :: condition_name

!-------------------------------------------------------------------
!
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
      case(CLM_2DTOP_TO_PF_2DTOP)
        map => pflotran_model%map_clm_2dtop_to_pf_2dtop
        source_mesh_id = CLM_FACE_MESH
        dest_mesh_id = PF_FACE_MESH
        condition_name = 'clm_gflux_bc'

      case(PF_2DTOP_TO_CLM_2DTOP)
        map => pflotran_model%map_pf_2dtop_to_clm_2dtop
        source_mesh_id = PF_FACE_MESH
        dest_mesh_id = CLM_FACE_MESH
        condition_name = 'clm_gflux_bc'

      case(CLM_2DBOT_TO_PF_2DBOT)
        map => pflotran_model%map_clm_2dbot_to_pf_2dbot
        source_mesh_id = CLM_FACE_MESH
        dest_mesh_id = PF_FACE_MESH
        condition_name = 'clm_bflux_bc'

      case(PF_2DBOT_TO_CLM_2DBOT)
        map => pflotran_model%map_pf_2dbot_to_clm_2dbot
        source_mesh_id = PF_FACE_MESH
        dest_mesh_id = CLM_FACE_MESH
        condition_name = 'clm_bflux_bc'

      case default
        option%io_buffer = 'map_id argument NOT yet support ' // &
          'pflotranModelInitMappingFaceToFace'
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

    found=PETSC_FALSE
    grid_pf_npts_local = 0
    grid_pf_npts_ghost = 0

    ! find the specified face of PFLOTRAN domain, by BC condition name
    if (source_mesh_id == PF_FACE_MESH .or. &
        dest_mesh_id == PF_FACE_MESH) then

        patch => realization%patch
        grid => patch%grid

        ! mesh is PF_FACE_MESH, 'FACE' type is defined by BC 'condition_name'
        boundary_condition => patch%boundary_condition_list%first

        do
          if (.not.associated(boundary_condition)) exit
          cur_connection_set => boundary_condition%connection_set

          if(StringCompare(trim(boundary_condition%name),trim(condition_name))) then

            found=PETSC_TRUE

            ! Allocate memory to save cell ids and flag for local cells
            allocate(grid_pf_cell_ids_nindex(cur_connection_set%num_connections))
            allocate(grid_pf_local_nindex(cur_connection_set%num_connections))
            grid_pf_npts_local = cur_connection_set%num_connections

            ! Save cell ids in application order 0-based
            do iconn=1,cur_connection_set%num_connections

              local_id = cur_connection_set%id_dn(iconn)
              ghosted_id = grid%nL2G(local_id)
              if (patch%imat(ghosted_id) <= 0) cycle
              grid_pf_cell_ids_nindex(iconn) = grid%nG2A(ghosted_id) - 1
              grid_pf_local_nindex(iconn) = 1
            enddo

          endif
          boundary_condition => boundary_condition%next
        enddo

        if(.not.found) then
          pflotran_model%option%io_buffer = 'condition name not found in boundary conditions'
          call printErrMsg(pflotran_model%option)
        endif

      else
        option%io_buffer='Unknown PFLOTRAN face mesh id'
        call printErrMsg(option)

    endif

    !
    ! Step-1: Find face cells-ids of PFLOTRAN subsurface domain
    !
    call VecCreateMPI(option%mycomm, grid%nlmax, PETSC_DETERMINE, face_ids, ierr)
    CHKERRQ(ierr)
    call VecSet(face_ids, -1.d0, ierr); CHKERRQ(ierr)

    ! Set 1.0 to all cells that make up a face of PFLOTRAN subsurface domain
    allocate(v_loc(grid_pf_npts_local))
    v_loc = 1.d0
    call VecSetValues(face_ids, grid_pf_npts_local, grid_pf_cell_ids_nindex, &
                      v_loc, INSERT_VALUES, ierr); CHKERRQ(ierr)
    deallocate(v_loc)

    call VecAssemblyBegin(face_ids, ierr); CHKERRQ(ierr)
    call VecAssemblyEnd(face_ids, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(face_ids, v_loc, ierr); CHKERRQ(ierr)
    count = 0
    do local_id=1,grid%nlmax
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr); CHKERRQ(ierr)

    count = 0
    do local_id=1,grid%nlmax
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif

    enddo
    call VecRestoreArrayF90(face_ids, v_loc, ierr); CHKERRQ(ierr)

    !
    allocate(int_array(grid_pf_npts_local))
    do iconn = 1, grid_pf_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, grid_pf_npts_local, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr); CHKERRQ(ierr)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm, grid_pf_npts_local, grid_pf_cell_ids_nindex, &
                         PETSC_COPY_VALUES, is_from, ierr); CHKERRQ(ierr)

    ! create scatter context
    call VecCreateSeq(PETSC_COMM_SELF, grid_pf_npts_local, face_ids_loc, &
      ierr); CHKERRQ(ierr)

    call VecScatterCreate(face_ids, is_from, face_ids_loc, is_to, vec_scat, &
                          ierr); CHKERRQ(ierr)
    call ISDestroy(is_from, ierr); CHKERRQ(ierr)
    call ISDestroy(is_to, ierr); CHKERRQ(ierr)

    call VecScatterBegin(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterEnd(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(face_ids_loc, v_loc, ierr); CHKERRQ(ierr)
    count = 0
    do iconn = 1, grid_pf_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_pf_cell_ids_nindex(count) = INT(v_loc(iconn))
      endif

    enddo
    call VecRestoreArrayF90(face_ids_loc, v_loc, ierr); CHKERRQ(ierr)
    call VecDestroy(face_ids_loc, ierr); CHKERRQ(ierr)

    !
    ! Step-2: Recompute 'map%s2d_i/jscr' for pf mesh
    !
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_nwts, face_ids_loc, ierr)
    CHKERRQ(ierr)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                     PETSC_COPY_VALUES, is_to, ierr); CHKERRQ(ierr)

    do iconn = 1, map%s2d_nwts
      if (dest_mesh_id == PF_FACE_MESH) then
         int_array(iconn) = map%s2d_icsr(iconn)
      elseif (source_mesh_id == PF_FACE_MESH) then
         int_array(iconn) = map%s2d_jcsr(iconn)
      endif
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                 PETSC_COPY_VALUES, is_from, ierr); CHKERRQ(ierr)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(face_ids, is_from, face_ids_loc, is_to, vec_scat, &
                          ierr); CHKERRQ(ierr)
    call ISDestroy(is_from, ierr); CHKERRQ(ierr)
    call ISDestroy(is_to, ierr); CHKERRQ(ierr)

    call VecScatterBegin(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr); CHKERRQ(ierr)

    call VecScatterEnd(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(face_ids_loc, v_loc, ierr); CHKERRQ(ierr)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        if (dest_mesh_id == PF_FACE_MESH) then
          map%s2d_icsr(count) = INT(v_loc(iconn))
        elseif (source_mesh_id == PF_FACE_MESH) then
          map%s2d_jcsr(count) = INT(v_loc(iconn))
        endif
      endif

    enddo
    call VecRestoreArrayF90(face_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    call VecDestroy(face_ids_loc, ierr)
    CHKERRQ(ierr)

    if(count /= map%s2d_nwts) then
      option%io_buffer='No. of face cells in mapping dataset does not ' // &
        'match face cells on which BC is applied - PFLOTRAN.'
      call printErrMsg(option)
    endif
    call VecDestroy(face_ids, ierr); CHKERRQ(ierr)

    !
    ! Step-3: Find face cells-ids of CLM soil/below-ground domain
    !
    allocate(v_loc(grid_clm_npts_local))
    v_loc = 1.d0
    call VecCreateSeq(PETSC_COMM_SELF, grid_clm_npts_local, face_ids_loc, ierr)
    CHKERRQ(ierr)
    call VecCreateMPI(option%mycomm, clm_pf_idata%nlclm_sub, PETSC_DECIDE, face_ids, ierr)
    CHKERRQ(ierr)
    call VecSet(face_ids, -1.d0, ierr)
    CHKERRQ(ierr)

    ! Set 1.0 to all cells that make up surface of CLM subsurface domain
    call VecSetValues(face_ids, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                      v_loc, INSERT_VALUES, ierr)
    CHKERRQ(ierr)

    deallocate(v_loc)
    call VecAssemblyBegin(face_ids, ierr)
    CHKERRQ(ierr)
    call VecAssemblyEnd(face_ids, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(face_ids, v_loc, ierr)
    CHKERRQ(ierr)
    count = 0
    do local_id=1,clm_pf_idata%nlclm_sub
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr)
    CHKERRQ(ierr)

    count = 0
    do local_id=1,clm_pf_idata%nlclm_sub
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif
    enddo
    call VecRestoreArrayF90(face_ids, v_loc, ierr)
    CHKERRQ(ierr)

    !
    allocate(int_array(grid_clm_npts_local))
    do iconn = 1, grid_clm_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, grid_clm_npts_local, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    CHKERRQ(ierr)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                         PETSC_COPY_VALUES, is_from, ierr)
    CHKERRQ(ierr)


    ! create scatter context
    call VecScatterCreate(face_ids, is_from, face_ids_loc, is_to, vec_scat, &
                          ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_from, ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_to, ierr)
    CHKERRQ(ierr)

    call VecScatterBegin(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterEnd(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(face_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    count = 0
    do iconn = 1, grid_clm_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_clm_cell_ids_nindex_copy(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(face_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    call VecDestroy(face_ids_loc, ierr)
    CHKERRQ(ierr)

    !
    ! Step-4: Recompute 'map%s2d_i/jscr' for clm mesh
    !
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_nwts, face_ids_loc, ierr)
    CHKERRQ(ierr)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    CHKERRQ(ierr)


    do iconn = 1, map%s2d_nwts
      if (source_mesh_id == CLM_FACE_MESH) then
         int_array(iconn) = map%s2d_jcsr(iconn)
      elseif (dest_mesh_id == CLM_FACE_MESH) then
         int_array(iconn) = map%s2d_icsr(iconn)
      endif
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_from, ierr)
    CHKERRQ(ierr)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(face_ids, is_from, face_ids_loc, is_to, vec_scat, &
                          ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_from, ierr)
    CHKERRQ(ierr)
    call ISDestroy(is_to, ierr)
    CHKERRQ(ierr)

    call VecScatterBegin(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterEnd(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(face_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        if (source_mesh_id == CLM_FACE_MESH) then
           map%s2d_jcsr(count) = INT(v_loc(iconn))
        elseif (dest_mesh_id == CLM_FACE_MESH) then
           map%s2d_icsr(count) = INT(v_loc(iconn))
        endif
      endif
    enddo
    call VecRestoreArrayF90(face_ids_loc, v_loc, ierr)
    CHKERRQ(ierr)
    call VecDestroy(face_ids_loc, ierr)
    CHKERRQ(ierr)

    if(count /= map%s2d_nwts) then
      option%io_buffer='No. of face cells in mapping dataset does not ' // &
        'match face cells on which BC is applied - CLM.'
      call printErrMsg(option)
    endif

    call VecDestroy(face_ids, ierr)
    CHKERRQ(ierr)

    !
    select case(source_mesh_id)
      case(CLM_FACE_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                         grid_clm_npts_ghost, &
                                         grid_clm_cell_ids_nindex_copy, &
                                         grid_clm_local_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_FACE_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_pf_npts_local, &
                                         grid_pf_npts_ghost, &
                                         grid_pf_cell_ids_nindex, &
                                         grid_pf_local_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_clm_npts_local, &
                                              grid_clm_npts_ghost, &
                                              grid_clm_cell_ids_nindex_copy, &
                                              grid_clm_local_nindex)
      case default
        option%io_buffer = 'Invalid argument source_mesh_id passed to ' // &
          'pflotranModelInitMappingFaceToFace'
        call printErrMsg(option)
    end select

    deallocate(grid_pf_cell_ids_nindex)
    deallocate(grid_pf_local_nindex)
    deallocate(grid_clm_cell_ids_nindex_copy)
    deallocate(grid_clm_local_nindex)

    call MappingDecompose(map, option)
    call MappingFindDistinctSourceMeshCellIds(map, option)
    call MappingCreateWeightMatrix(map, option)
    call MappingCreateScatterOfSourceMesh(map, option)
    call MappingFreeNotNeeded(map)

    ! Setting the number of cells constituting the face of the 3D
    ! subsurface domain for each model.
    select case(map_id)
      case(CLM_2DTOP_TO_PF_2DTOP, PF_2DTOP_TO_CLM_2DTOP)
        clm_pf_idata%nlclm_2dtop = grid_clm_npts_local
        clm_pf_idata%ngclm_2dtop = grid_clm_npts_local
        clm_pf_idata%nlpf_2dtop  = grid_pf_npts_local
        clm_pf_idata%ngpf_2dtop  = grid_pf_npts_local
      case(CLM_2DBOT_TO_PF_2DBOT, PF_2DBOT_TO_CLM_2DBOT)
        clm_pf_idata%nlclm_2dbot = grid_clm_npts_local
        clm_pf_idata%ngclm_2dbot = grid_clm_npts_local
        clm_pf_idata%nlpf_2dbot  = grid_pf_npts_local
        clm_pf_idata%ngpf_2dbot  = grid_pf_npts_local
      case default
        option%io_buffer = 'map_id argument NOT yet supported in ' // &
                        'pflotranModelInitMappingFaceToFace'
        call printErrMsg(option)
    end select

  end subroutine pflotranModelInitMapFaceToFace

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  
end module pflotran_clm_setmapping_module

