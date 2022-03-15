module PM_Material_Transform_class

! MODULE DESCRIPTION:
! ===========================================================================
! This process model incorporates material transformations as surrogate models
! for physical phenomena given evolving system conditions
! ===========================================================================

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
  use petscsys
  use PM_Base_class
  use Realization_Subsurface_class
  use Option_module
  use PFLOTRAN_Constants_module
  use Material_Transform_module

  implicit none

  private

! OBJECT pm_material_transform_type:
! ==================================
! ---------------------------------------------------------------------------
! Description:  This is the material transform process model object. It contains
! a list of material transformation objects that may incorporate one or more
! surrogate models for physical phenomena.
! ---------------------------------------------------------------------------
! realization: pointer to subsurface realization object
! mtl: pointer to linked list of material transforms in the pm block
! --------------------------------------------------------------------------  
  type, public, extends(pm_base_type) :: pm_material_transform_type
    class(realization_subsurface_type), pointer :: realization
    type(material_transform_type), pointer :: mtl
  contains
    procedure, public :: Setup => PMMaterialTransformSetup
    procedure, public :: ReadPMBlock => PMMaterialTransformReadPMBlock
    procedure, public :: SetRealization => PMMaterialTransformSetRealization
    procedure, public :: InitializeRun => PMMaterialTransformInitializeRun
    procedure, public :: FinalizeRun => PMMaterialTransformFinalizeRun
    procedure, public :: InitializeTimestep => PMMaterialTransformInitializeTS
    procedure, public :: FinalizeTimestep => PMMaterialTransformFinalizeTS
    procedure, public :: UpdateSolution => PMMaterialTransformUpdateSolution
    procedure, public :: Solve => PMMaterialTransformSolve
    procedure, public :: TimeCut => PMMaterialTransformTimeCut
    procedure, public :: UpdateAuxVars => PMMaterialTransformUpdateAuxVars
    procedure, public :: CheckpointHDF5 => PMMaterialTransformCheckpointHDF5
    procedure, public :: CheckpointBinary => PMMaterialTransformCheckpointBinary
    procedure, public :: RestartHDF5 => PMMaterialTransformRestartHDF5
    procedure, public :: RestartBinary => PMMaterialTransformRestartBinary
    procedure, public :: InputRecord => PMMaterialTransformInputRecord
    procedure, public :: Destroy => PMMaterialTransformDestroy
  end type pm_material_transform_type

  public :: PMMaterialTransformCreate

contains

! ************************************************************************** !

function PMMaterialTransformCreate()
  !
  ! Creates process model
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022
  !

  implicit none

  class(pm_material_transform_type), pointer :: PMMaterialTransformCreate

  class(pm_material_transform_type), pointer :: pm

  allocate(pm)
  call PMBaseInit(pm)

  pm%header = 'MATERIAL TRANSFORM'

  nullify(pm%realization)
  nullify(pm%mtl)

  PMMaterialTransformCreate => pm

end function PMMaterialTransformCreate

! ************************************************************************** !

subroutine PMMaterialTransformSetRealization(this, realization)
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! realization (input): pointer to subsurface realization object
  ! ----------------------------------
  class(pm_material_transform_type) :: this
  class(realization_subsurface_type), pointer :: realization
  ! ----------------------------------

  this%realization => realization
  this%realization_base => realization

end subroutine PMMaterialTransformSetRealization

! ************************************************************************** !

subroutine PMMaterialTransformSetup(this)
  !
  ! Sets up auxiliary process model
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Material_module
  use Grid_module

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! ----------------------------------
  class(pm_material_transform_type) :: this
  ! ----------------------------------
  ! LOCAL VARIABLES:
  ! ================
  ! patch: pointer to patch object within realization
  ! option: pointer to option object within realization
  ! grid: pointer to grid object within realization
  ! mtf: pointer to material transform object within patch
  ! MT_auxvars: pointer to array of material transform auxiliary variables
  ! cur_material_property: pointer to material property within realization
  ! null_material_property: null pointer for regions without materials
  ! local_id: grid cell id number
  ! ghosted_id: ghosted grid cell id number
  ! material_id: id number of material
  ! ----------------------------------
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(material_transform_type), pointer :: mtf
  type(material_transform_auxvar_type), pointer :: MT_auxvars(:)
  type(material_property_type), pointer :: cur_material_property
  type(material_property_type), pointer :: null_material_property
  PetscInt :: local_id, ghosted_id, material_id
  ! ----------------------------------

  patch => this%realization%patch
  option => this%realization%option
  grid => patch%grid

  ! pass material transform list from PM to realization
  if (associated(this%mtl)) then
    call MaterialTransformAddToList(this%mtl, &
                                    this%realization%material_transform)
  endif

  ! set up mapping for material transform functions
  patch%material_transform => this%realization%material_transform
  call MaterialTransformConvertListToArray(patch%material_transform, &
                                           patch%material_transform_array, &
                                           option)

  ! material property mapping to PM Material Transform
  cur_material_property => this%realization%material_properties

  do
    if (.not. associated(cur_material_property)) exit

    ! material transform function id 
    if (associated(patch%material_transform_array)) then
      if (cur_material_property%mtf) then
        ! find ID
        cur_material_property%material_transform_id = &
          MaterialTransformGetID( &
            patch%material_transform_array, &
            cur_material_property%material_transform_name, &
            cur_material_property%name,option)
      endif
    endif

    ! check for errors (0 = not found, -999 = not applicable)
    if (cur_material_property%material_transform_id == 0) then
      option%io_buffer = 'Material transform function "' // &
        trim(cur_material_property%material_transform_name) // &
        '" not found for material "'//trim(cur_material_property%name) // '".'
      call PrintErrMsg(option)
    endif

    cur_material_property => cur_material_property%next

  enddo

  ! create null material property for inactive cells
  null_material_property => MaterialPropertyCreate(option)
  ! cell mapping for material transform id and initialize auxilary variables
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    material_id = patch%imat(ghosted_id)
    if (material_id <= 0) cycle

    ! get material property from id
    if (material_id == 0) then
      cur_material_property => null_material_property
    else if (abs(material_id) <= size(patch%material_property_array)) then
      ! error conditons already checked in InitSubsurfAssignMatProperties
      if (material_id < 0) then
        cur_material_property => null_material_property
      else
        cur_material_property => &
          patch%material_property_array(material_id)%ptr
      endif
    endif

    if (option%nflowdof > 0) then
      if (associated(patch%mtf_id)) then
        patch%mtf_id(ghosted_id) = &  
          cur_material_property%material_transform_id

        call RealLocalToLocalWithArray(this%realization,MTF_ID_ARRAY)

      endif
    endif

  enddo

  patch%aux%MT => MaterialTransformCreate()
  allocate(MT_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    material_id = patch%imat(ghosted_id)
    if (material_id <= 0) cycle

    call MaterialTransformAuxVarInit(MT_auxvars(ghosted_id))
    
    if (Initialized(patch%mtf_id(ghosted_id))) then
      ! pointer to material transform in patch ghosted id
      mtf => patch%material_transform_array(patch%mtf_id(ghosted_id))%ptr

      if (associated(mtf)) then
        if (associated(mtf%illitization)) then
          MT_auxvars(ghosted_id)%il_aux => IllitizationAuxVarInit(option)
        endif
        if (associated(mtf%buffer_erosion)) then
          MT_auxvars(ghosted_id)%be_aux => BufferErosionAuxVarInit()
        endif
      endif

    endif
    
  enddo
  patch%aux%MT%auxvars => MT_auxvars
  patch%aux%MT%num_aux = grid%ngmax

end subroutine PMMaterialTransformSetup

! ************************************************************************** !

subroutine PMMaterialTransformReadPMBlock(this,input)
  !
  ! Reads input file parameters associated with the material transform
  !   process model
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022

  use Input_Aux_module
  use Option_module
  use String_module
  use Material_Transform_module

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! input (input/output): pointer to input object
  ! ----------------------------------
  class(pm_material_transform_type) :: this
  type(input_type), pointer :: input
  ! ----------------------------------
  ! LOCAL VARIABLES:
  ! ================
  ! option: pointer to option object
  ! word: temporary string
  ! error_string: error message string
  ! material_transform: pointer to material transform object
  ! prev_material_transform: pointer for linked list creation
  ! ----------------------------------
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  class(material_transform_type), pointer :: material_transform
  class(material_transform_type), pointer :: prev_material_transform
  ! ----------------------------------

  option => this%option
  
  option%io_buffer = 'pflotran card:: MATERIAL_TRANSFORM_GENERAL'
  call PrintMsg(option)
  
  input%ierr = 0
  
  nullify(prev_material_transform)
  call InputPushBlock(input,option)
  do
    error_string = 'MATERIAL_TRANSFORM_GENERAL'
    
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    select case(trim(word))
    !-------------------------------------
      case('MATERIAL_TRANSFORM')
        error_string = 'MATERIAL_TRANSFORM_GENERAL, MATERIAL_TRANSFORM'
        material_transform => MaterialTransformCreate()
        call InputReadWord(input,option,material_transform%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        
        option%io_buffer = '  MATERIAL_TRANSFORM :: ' // &
                             trim(material_transform%name)
        call PrintMsg(option)
        
        error_string = 'MATERIAL_TRANSFORM, ' // trim(material_transform%name)
        
        call MaterialTransformRead(material_transform,input,option)
        
        if (associated(prev_material_transform)) then
          prev_material_transform%next => material_transform
        else
          this%mtl => material_transform
        endif
        prev_material_transform => material_transform
        
        nullify(material_transform)
    !-------------------------------------
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    !-------------------------------------
    end select
  enddo
  
  call InputPopBlock(input,option)

end subroutine PMMaterialTransformReadPMBlock

! ************************************************************************** !

recursive subroutine PMMaterialTransformInitializeRun(this)
  !
  ! Initializes the time stepping
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022
  use Realization_Subsurface_class
  use Patch_module
  use Material_module
  use Material_Aux_module
  use Grid_module
  use Option_module

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! --------------------------------
  class(pm_material_transform_type) :: this
  ! --------------------------------
  ! LOCAL VARIABLES:
  ! ================
  ! patch: pointer to patch object within realization
  ! option: pointer to option object within realization
  ! grid: pointer to grid object within realization
  ! mtf: pointer to material transform object within patch
  ! material_auxvars: pointer to array of material auxiliary variables
  ! material aux: pointer to material auxiliary variable object in list
  ! MT_auxvars: pointer to array of material transform auxiliary variables
  ! MT_au: pointer to material transform auxiliary variable object in list
  ! cur_material_property: pointer to material property within realization
  ! null_material_property: null pointer for regions without materials
  ! local_id: grid cell id number
  ! ghosted_id: ghosted grid cell id number
  ! material_id: id number of material
  ! --------------------------------
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(material_transform_type), pointer :: mtf
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_auxvar_type), pointer :: material_aux
  type(material_transform_auxvar_type), pointer :: MT_auxvars(:)
  type(material_transform_auxvar_type), pointer :: MT_aux
  PetscInt :: local_id, ghosted_id, material_id
  ! ----------------------------------

  patch => this%realization%patch
  option => this%realization%option
  grid => patch%grid

  material_auxvars => patch%aux%Material%auxvars
  MT_auxvars => patch%aux%MT%auxvars

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    material_id = patch%imat(ghosted_id)
    if (material_id <= 0) cycle
    
    material_aux => material_auxvars(ghosted_id)
    MT_aux => MT_auxvars(ghosted_id)
    
    if (Initialized(patch%mtf_id(ghosted_id)) .and. &
        .not. option%restart_flag) then
      ! pointer to material transform in patch ghosted id
      if (associated(patch%material_transform_array)) then
        allocate(mtf)
        mtf => patch%material_transform_array(patch%mtf_id(ghosted_id))%ptr
        if (associated(MT_aux%il_aux) .and. associated(mtf)) then
          MT_aux%il_aux%fs0 = &
            mtf%illitization%illitization_function%fs0
          MT_aux%il_aux%fs = &
            mtf%illitization%illitization_function%fs0
        endif
        nullify(mtf)
      endif
      
    endif
    
  enddo

end subroutine PMMaterialTransformInitializeRun

! ************************************************************************** !

subroutine PMMaterialTransformInitializeTS(this)
  !
  ! Initializes the time step
  !
  ! Author: Alex Salazar III
  ! Date: 01/20/2022

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! --------------------------------
  class(pm_material_transform_type) :: this
  ! --------------------------------

end subroutine PMMaterialTransformInitializeTS

! ************************************************************************** !

subroutine PMMaterialTransformFinalizeTS(this)
  !
  ! Finalizes the time step
  !
  ! Author: Alex Salazar III
  ! Date: 01/20/2022

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! --------------------------------
  class(pm_material_transform_type) :: this
  ! --------------------------------

  call RealizationUpdatePropertiesTS(this%realization)

end subroutine PMMaterialTransformFinalizeTS

! ************************************************************************** !

recursive subroutine PMMaterialTransformFinalizeRun(this)
  !
  ! Finalizes the run
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022
  !

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! --------------------------------
  class(pm_material_transform_type) :: this
  ! --------------------------------

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  

end subroutine PMMaterialTransformFinalizeRun

! ************************************************************************** !

subroutine PMMaterialTransformUpdateSolution(this)
  !
  ! Updates data in process model after a successful time step
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! ---------------------------------
  class(pm_material_transform_type) :: this
  ! ---------------------------------

end subroutine PMMaterialTransformUpdateSolution

! ************************************************************************** !

subroutine PMMaterialTransformUpdateAuxVars(this)
  !
  ! Updates the auxiliary variables associated with the process model
  !
  ! Author: Alex Salazar III
  ! Date: 03/03/2022

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! ----------------------------------
  class(pm_material_transform_type) :: this
  ! ----------------------------------

end subroutine PMMaterialTransformUpdateAuxVars

! ************************************************************************** !

subroutine PMMaterialTransformTimeCut(this)
  !
  ! Resets arrays for time step cut
  !
  ! Author: Alex Salazar III
  ! Date: 03/03/2022

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! ----------------------------------
  class(pm_material_transform_type) :: this
  ! ----------------------------------

end subroutine PMMaterialTransformTimeCut

! ************************************************************************** !

subroutine PMMaterialTransformSolve(this, time, ierr)
  !
  ! Updates materials based on the active models
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022
  !
  use Realization_Subsurface_class
  use Patch_module
  use Material_module
  use Material_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! time (input): [sec] simulation time
  ! ierr (input/output): [-] PETSc error integer
  ! ---------------------------------
  class(pm_material_transform_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  ! ---------------------------------
  ! LOCAL VARIABLES:
  ! ================
  ! patch: pointer to patch object within realization
  ! option: pointer to option object within realization
  ! grid: pointer to grid object within realization
  ! mtf: pointer to material transform object within patch
  ! material_auxvars: pointer to array of material auxiliary variables
  ! material aux: pointer to material auxiliary variable object in list
  ! global_auxvars: pointer to array of global auxiliary variables
  ! global aux: pointer to global auxiliary variable object in list
  ! MT_auxvars: pointer to array of material transform auxiliary variables
  ! MT_aux: pointer to material transform auxiliary variable object in list
  ! local_id: grid cell id number
  ! ghosted_id: ghosted grid cell id number
  ! material_id: id number of material
  ! ---------------------------------
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  class(material_transform_type), pointer :: mtf
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_auxvar_type), pointer :: material_aux
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_aux
  type(material_transform_auxvar_type), pointer :: MT_auxvars(:)
  type(material_transform_auxvar_type), pointer :: MT_aux
  PetscInt :: local_id, ghosted_id, material_id
  ! ----------------------------------

  patch => this%realization%patch
  option => this%realization%option
  grid => patch%grid

  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  MT_auxvars => patch%aux%MT%auxvars

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    material_id = patch%imat(ghosted_id)
    if (material_id <= 0) cycle

    global_aux => global_auxvars(ghosted_id)
    material_aux => material_auxvars(ghosted_id)
    MT_aux => MT_auxvars(ghosted_id)

    if (Initialized(patch%mtf_id(ghosted_id))) then
      ! pointer to material transform in patch ghosted id
      mtf => patch%material_transform_array(patch%mtf_id(ghosted_id))%ptr

      if (associated(mtf)) then
        if (associated(mtf%illitization)) then
          call mtf%illitization%illitization_function%CalculateILT( &
                 MT_aux%il_aux%fs, &
                 global_aux%temp, &
                 option%dt, &
                 MT_aux%il_aux%fi, &
                 MT_aux%il_aux%scale, &
                 option)
          call mtf%illitization%illitization_function%ShiftPerm( &
                 material_aux, &
                 MT_aux%il_aux, &
                 option)
        endif
        ! if (associated(mtf%buffer_erosion)) then
        ! endif
      endif
    endif

  enddo

  ierr = 0

end subroutine PMMaterialTransformSolve

! ************************************************************************** !

subroutine PMMaterialTransformCheckpointHDF5(this, pm_grp_id)
  ! 
  ! Checkpoints data associated with the material transform process model
  !
  ! Author: Alex Salazar III
  ! Date: 01/20/2022
  !

  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec
  use Field_module
  use Discretization_module
  use Variables_module, only: SMECTITE

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! pm_grp_id: file id number
  ! ----------------------------------
  class(pm_material_transform_type) :: this
  integer(HID_T) :: pm_grp_id
  ! ----------------------------------
  ! LOCAL VARIABLES:
  ! ================
  ! is: PETSc index set
  ! scatter_ctx: copies an MPI vector to sequential vectors on all MPI ranks
  ! local_mt_vec: local vector of material transform checkpoint values
  ! global_mt_vec: global vector of material transform checkpoint values
  ! ierr: I/O status indicator
  ! local_stride: local number of vector elements between start and end of block
  ! local_stride_tmp: temporary counter for local stride
  ! stride: global number of vector elements between start and end of block
  ! n_mt_local: number of material transforms on process
  ! n_mt_global: number of material transforms globally
  ! n_check_vars: number of values to checkpoint
  ! i, j: iterators
  ! indices: indices of the local material transform vector
  ! int_array: keeps track of the material transform number
  ! check_vars: array of checkpointed values
  ! cur_mt: material transform object
  ! dataset_name: descriptor of the material transform checkpoint data
  ! global_vec: global discretization PETSc vector 
  ! natural_vec: local discretization PETSc vector 
  ! check_il: logical check for presence of illitization functions in the
  !   material transform objects so auxiliary variables can be checkpointed
  ! check_be: logical check for presence of buffer erosion models in the
  !   material transform objects so auxiliary variables can be checkpointed
  ! option: option object
  ! field: field object
  ! discretization: discretization object
  ! ----------------------------------
  IS :: is
  VecScatter :: scatter_ctx
  Vec :: local_mt_vec
  Vec :: global_mt_vec
  PetscErrorCode :: ierr
  PetscInt :: local_stride
  PetscInt :: local_stride_tmp
  PetscInt :: stride
  PetscInt :: n_mt_local
  PetscInt :: n_mt_global
  PetscInt :: n_check_vars
  PetscInt :: i, j
  PetscInt, allocatable :: indices(:)
  PetscInt, allocatable :: int_array(:)
  PetscReal, allocatable :: check_vars(:)
  class(material_transform_type), pointer :: cur_mt
  character(len=MAXSTRINGLENGTH) :: dataset_name
  Vec :: global_vec
  Vec :: natural_vec
  PetscBool :: check_il
  PetscBool :: check_be
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  ! ----------------------------------

  option => this%realization%option
  field => this%realization%field
  discretization => this%realization%discretization

  local_stride = 0
  local_stride_tmp = 0
  n_mt_local = 0
  n_mt_global = 0

  ! current information in PM Material Transform that is checkpointed:
  !   (1) num_aux
  n_check_vars = 1 !number of scalar checkpoint variables
  
  cur_mt => this%mtl
  do 
    if (.not. associated(cur_mt)) exit
    n_mt_local = n_mt_local + 1
    local_stride_tmp = local_stride_tmp + n_check_vars
    cur_mt => cur_mt%next
    if (local_stride_tmp > local_stride) then
      local_stride = local_stride_tmp
    endif
    local_stride_tmp = 0
  enddo

  allocate(int_array(n_mt_local))
  cur_mt => this%mtl
  i = 1
  do
    if (.not. associated(cur_mt)) exit
    int_array(i) = i - 1
    i = i + 1
    cur_mt => cur_mt%next
  enddo

  ! gather relevant information from all processes
  call MPI_Allreduce(local_stride, stride, ONE_INTEGER_MPI, &
                     MPI_INTEGER, MPI_MAX, this%option%mycomm, ierr)
  call MPI_Allreduce(n_mt_local, n_mt_global, ONE_INTEGER_MPI, &
                     MPI_INTEGER, MPI_SUM, this%option%mycomm, ierr)   

  ! create MPI vector and sequential vector for mapping
  call VecCreateMPI(this%option%mycomm, n_mt_local*stride, n_mt_global*stride, & 
                    global_mt_vec,ierr); CHKERRQ(ierr)
  call VecCreateSeq(PETSC_COMM_SELF, n_mt_local*stride, local_mt_vec, ierr); &
         CHKERRQ(ierr)
  call VecSetBlockSize(global_mt_vec, stride, ierr); CHKERRQ(ierr)
  call VecSetBlockSize(local_mt_vec, stride, ierr); CHKERRQ(ierr)

  allocate(check_vars(stride))
  allocate(indices(stride))

  ! collect data for checkpointing
  j = 1
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit

    check_vars(1) = cur_mt%num_aux ! checkpoint #1

    i = n_check_vars + 1
    do
      if (i > stride) exit
      check_vars(i) = -9999
      i = i + 1
    enddo

    do i = 1, stride
      indices(i) = (j - 1)*stride + i - 1
    enddo
    j = j + 1

    call VecSetValues(local_mt_vec, stride, indices, check_vars, &
                     INSERT_VALUES, ierr); CHKERRQ(ierr)

    cur_mt => cur_mt%next

  enddo

  !Create map and add values from the sequential vector to the global 
  call ISCreateBlock(this%option%mycomm, stride, n_mt_local, int_array, &
                     PETSC_COPY_VALUES, is, ierr); CHKERRQ(ierr)
  call VecScatterCreate(local_mt_vec, PETSC_NULL_IS, global_mt_vec, &
                        is, scatter_ctx, ierr); CHKERRQ(ierr)
  call VecScatterBegin(scatter_ctx, local_mt_vec, global_mt_vec, &  
                       INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx, local_mt_vec, global_mt_vec, & 
                     INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)

  ! write the checkpoint file
  dataset_name='material transform model info'
  call HDF5WriteDataSetFromVec(dataset_name, this%option, global_mt_vec,&
                               pm_grp_id, H5T_NATIVE_DOUBLE)
  call VecScatterDestroy(scatter_ctx, ierr); CHKERRQ(ierr)
  call ISDestroy(is, ierr); CHKERRQ(ierr)
  call VecDestroy(global_mt_vec, ierr); CHKERRQ(ierr)
  call VecDestroy(local_mt_vec, ierr); CHKERRQ(ierr)

  ! checkpoint the auxiliary variables
  check_il = PETSC_FALSE
  check_be = PETSC_FALSE
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit
    if (associated(cur_mt%illitization)) then
      check_il = PETSC_TRUE
    endif
    if (associated(cur_mt%buffer_erosion)) then
      check_be = PETSC_TRUE
    endif
    cur_mt => cur_mt%next
  enddo

  if (check_il) then
    global_vec = PETSC_NULL_VEC
    call DiscretizationCreateVector(this%realization%discretization, ONEDOF, &
                                    global_vec, GLOBAL, option)
    call DiscretizationCreateVector(this%realization%discretization, ONEDOF, &
                                    natural_vec, NATURAL, option)
    call MaterialTransformGetAuxVarVecLoc(this%realization%patch%aux%MT, &
                                          field%work_loc, SMECTITE, &
                                          ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization, field%work_loc, &
                                     global_vec, ONEDOF)
    call DiscretizationGlobalToNatural(discretization, global_vec, &
                                       natural_vec, ONEDOF)
    dataset_name = "Smectite" // CHAR(0)
    call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
                                 pm_grp_id, H5T_NATIVE_DOUBLE)
    call VecDestroy(global_vec, ierr); CHKERRQ(ierr)
    call VecDestroy(natural_vec, ierr); CHKERRQ(ierr)
  endif
  ! if (check_be) then
  ! endif

end subroutine PMMaterialTransformCheckpointHDF5

! ************************************************************************** !

subroutine PMMaterialTransformRestartHDF5(this, pm_grp_id)
  ! 
  ! Restarts data associated with material transform process model
  ! 
  ! Author: Alex Salazar III
  ! Date: 01/20/2022
  !

  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use hdf5
  use HDF5_module, only : HDF5ReadDataSetInVec
  use Field_module
  use Discretization_module
  use Variables_module, only: SMECTITE

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! pm_grp_id: file id number
  ! ----------------------------------
  class(pm_material_transform_type) :: this
  integer(HID_T) :: pm_grp_id
  ! ----------------------------------
  ! LOCAL VARIABLES:
  ! ================
  ! is: PETSc index set
  ! scatter_ctx: copies an MPI vector to sequential vectors on all MPI ranks
  ! local_mt_vec: local vector of material transform checkpoint values
  ! global_mt_vec: global vector of material transform checkpoint values
  ! ierr: I/O status indicator
  ! local_stride: local number of vector elements between start and end of block
  ! local_stride_tmp: temporary counter for local stride
  ! stride: global number of vector elements between start and end of block
  ! n_mt_local: number of material transforms on process
  ! n_mt_global: number of material transforms globally
  ! n_check_vars: number of values to checkpoint
  ! i: iterator
  ! indices: indices of the local material transform vector
  ! int_array: keeps track of the material transform number
  ! check_vars: array of checkpointed values
  ! local_mt_array: data converted into a Fortran array
  ! cur_mt: material transform object
  ! dataset_name: descriptor of the material transform checkpoint data
  ! global_vec: global discretization PETSc vector 
  ! natural_vec: local discretization PETSc vector 
  ! check_il: logical check for presence of illitization functions in the
  !   material transform objects so auxiliary variables can be checkpointed
  ! check_be: logical check for presence of buffer erosion models in the
  !   material transform objects so auxiliary variables can be checkpointed
  ! option: option object
  ! field: field object
  ! discretization: discretization object
  ! ----------------------------------
  IS :: is
  VecScatter :: scatter_ctx
  Vec :: local_mt_vec
  Vec :: global_mt_vec
  PetscErrorCode :: ierr
  PetscInt :: local_stride
  PetscInt :: n_mt_local
  PetscInt :: n_mt_global
  PetscInt :: n_check_vars
  PetscInt :: local_stride_tmp
  PetscInt :: i
  PetscInt :: stride
  PetscInt, allocatable :: indices(:)
  PetscInt, allocatable :: int_array(:)
  PetscReal, allocatable :: check_vars(:)
  PetscReal, pointer :: local_mt_array(:)
  class(material_transform_type), pointer :: cur_mt
  character(len=MAXSTRINGLENGTH) :: dataset_name
  Vec :: global_vec
  Vec :: natural_vec
  PetscBool :: check_il
  PetscBool :: check_be
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  ! ----------------------------------

  option => this%realization%option
  field => this%realization%field
  discretization => this%realization%discretization

  local_stride = 0
  local_stride_tmp = 0
  n_mt_local = 0
  n_mt_global = 0

  ! current information in PM Material Transform that is checkpointed:
  !   (1) num_aux
  n_check_vars = 1 !number of scalar checkpoint variables

  cur_mt => this%mtl
  do 
    if (.not. associated(cur_mt)) exit
    n_mt_local = n_mt_local + 1
    local_stride_tmp = local_stride_tmp + n_check_vars
    cur_mt => cur_mt%next
    if (local_stride_tmp > local_stride) then
      local_stride = local_stride_tmp
    endif
    local_stride_tmp = 0
  enddo

  allocate(int_array(n_mt_local))
  i = 1
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit
    int_array(i) = i - 1
    i = i + 1
    cur_mt => cur_mt%next
  enddo

  ! gather relevant information from all processes
  call MPI_Allreduce(local_stride, stride, ONE_INTEGER_MPI, &
                     MPI_INTEGER, MPI_MAX, this%option%mycomm, ierr)
  call MPI_Allreduce(n_mt_local, n_mt_global, ONE_INTEGER_MPI, &
                     MPI_INTEGER, MPI_SUM, this%option%mycomm, ierr)   

  ! create MPI vector for HDF5 reading and sequential vector for mt information
  !   stored in the process
  call VecCreateMPI(this%option%mycomm, n_mt_local*stride, n_mt_global*stride, & 
                    global_mt_vec,ierr); CHKERRQ(ierr)
  call VecCreateSeq(PETSC_COMM_SELF, n_mt_local*stride, local_mt_vec, ierr); &
         CHKERRQ(ierr)
  call VecSetBlockSize(global_mt_vec, stride, ierr); CHKERRQ(ierr)
  call VecSetBlockSize(local_mt_vec, stride, ierr); CHKERRQ(ierr)

  ! read data from HDF5
  dataset_name='material transform model info'
  call HDF5ReadDataSetInVec(dataset_name, this%option, global_mt_vec, &
                            pm_grp_id, H5T_NATIVE_DOUBLE)

  ! create mapping between MPI and sequential vectors
  call ISCreateBlock(this%option%mycomm, stride, n_mt_local, int_array, &
                     PETSC_COPY_VALUES, is, ierr); CHKERRQ(ierr)
  call VecScatterCreate(global_mt_vec, is, local_mt_vec, &
                        PETSC_NULL_IS, scatter_ctx, ierr); CHKERRQ(ierr)

  ! obtain data from the MPI vector
  call VecScatterBegin(scatter_ctx, global_mt_vec, local_mt_vec, &  
                       INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx, global_mt_vec, local_mt_vec, & 
                     INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)

  ! convert the data into a Fortran array
  call VecGetArrayF90(local_mt_vec, local_mt_array, ierr); CHKERRQ(ierr)

  ! assign checkpointed material transform information
  i = 1
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit

    cur_mt%num_aux = local_mt_array(i) ! checkpoint #1

    cur_mt => cur_mt%next
    i = i + stride
  enddo

  call VecRestoreArrayF90(local_mt_vec, local_mt_array, ierr); CHKERRQ(ierr)
  call VecScatterDestroy(scatter_ctx, ierr); CHKERRQ(ierr)
  call ISDestroy(is, ierr); CHKERRQ(ierr)
  call VecDestroy(global_mt_vec, ierr); CHKERRQ(ierr)
  call VecDestroy(local_mt_vec, ierr); CHKERRQ(ierr)

  ! retrieve the auxiliary variables
  check_il = PETSC_FALSE
  check_be = PETSC_FALSE
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit
    if (associated(cur_mt%illitization)) then
      check_il = PETSC_TRUE
    endif
    if (associated(cur_mt%buffer_erosion)) then
      check_be = PETSC_TRUE
    endif
    cur_mt => cur_mt%next
  enddo

  if (check_il) then
    global_vec = PETSC_NULL_VEC
    call DiscretizationCreateVector(this%realization%discretization, ONEDOF, &
                                    global_vec, GLOBAL, option)
    call DiscretizationCreateVector(this%realization%discretization, ONEDOF, &
                                    natural_vec, NATURAL, option)
    dataset_name = "Smectite" // CHAR(0)
    call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
                              pm_grp_id, H5T_NATIVE_DOUBLE)
    call DiscretizationNaturalToGlobal(discretization, natural_vec, &
                                       global_vec, ONEDOF)
    call DiscretizationGlobalToLocal(discretization, global_vec, &
                                     field%work_loc, ONEDOF)
    call MaterialTransformSetAuxVarVecLoc(this%realization%patch%aux%MT, &
                                          field%work_loc, SMECTITE, &
                                          ZERO_INTEGER)
    call VecDestroy(global_vec, ierr); CHKERRQ(ierr)
    call VecDestroy(natural_vec, ierr); CHKERRQ(ierr)
  endif
  ! if (check_be) then
  ! endif

end subroutine PMMaterialTransformRestartHDF5

! ************************************************************************** !

subroutine PMMaterialTransformCheckpointBinary(this, viewer)
  ! 
  ! Checkpoints data associated with the material transform process model
  !
  ! Author: Alex Salazar III
  ! Date: 01/20/2022
  ! 

  use petscvec
  use Option_module
  use Discretization_module
  use Grid_module
  use Field_module
  use Discretization_module
  use Variables_module, only: SMECTITE
  
  implicit none
  
  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! viewer: petsc viewer variable
  ! ----------------------------------
  PetscViewer :: viewer
  class(pm_material_transform_type) :: this
  ! ----------------------------------
  ! LOCAL VARIABLES:
  ! ================
  ! is: PETSc index set
  ! scatter_ctx: copies an MPI vector to sequential vectors on all MPI ranks
  ! local_mt_vec: local vector of material transform checkpoint values
  ! global_mt_vec: global vector of material transform checkpoint values
  ! ierr: I/O status indicator
  ! local_stride: local number of vector elements between start and end of block
  ! local_stride_tmp: temporary counter for local stride
  ! stride: global number of vector elements between start and end of block
  ! n_mt_local: number of material transforms on process
  ! n_mt_global: number of material transforms globally
  ! n_check_vars: number of values to checkpoint
  ! i, j: iterators
  ! indices: indices of the local material transform vector
  ! int_array: keeps track of the material transform number
  ! check_vars: array of checkpointed values
  ! cur_mt: material transform object
  ! dataset_name: descriptor of the material transform checkpoint data
  ! global_vec: global discretization PETSc vector 
  ! natural_vec: local discretization PETSc vector 
  ! check_il: logical check for presence of illitization functions in the
  !   material transform objects so auxiliary variables can be checkpointed
  ! check_be: logical check for presence of buffer erosion models in the
  !   material transform objects so auxiliary variables can be checkpointed
  ! option: option object
  ! field: field object
  ! discretization: discretization object
  ! ----------------------------------
  IS :: is
  VecScatter :: scatter_ctx
  Vec :: local_mt_vec
  Vec :: global_mt_vec
  PetscErrorCode :: ierr
  PetscInt :: local_stride
  PetscInt :: local_stride_tmp
  PetscInt :: stride
  PetscInt :: n_mt_local
  PetscInt :: n_mt_global
  PetscInt :: n_check_vars
  PetscInt :: i, j
  PetscInt, allocatable :: indices(:)
  PetscInt, allocatable :: int_array(:)
  PetscReal, allocatable :: check_vars(:)
  class(material_transform_type), pointer :: cur_mt
  character(len=MAXSTRINGLENGTH) :: dataset_name
  Vec :: global_vec
  Vec :: natural_vec
  PetscBool :: check_il
  PetscBool :: check_be
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  ! ----------------------------------

  option => this%realization%option
  field => this%realization%field
  discretization => this%realization%discretization

  local_stride = 0
  local_stride_tmp = 0
  n_mt_local = 0
  n_mt_global = 0

  ! current information in PM Material Transform that is checkpointed:
  !   (1) num_aux
  n_check_vars = 1 !number of scalar checkpoint variables
  
  cur_mt => this%mtl
  do 
    if (.not. associated(cur_mt)) exit
    n_mt_local = n_mt_local + 1
    local_stride_tmp = local_stride_tmp + n_check_vars
    cur_mt => cur_mt%next
    if (local_stride_tmp > local_stride) then
      local_stride = local_stride_tmp
    endif
    local_stride_tmp = 0
  enddo

  allocate(int_array(n_mt_local))
  cur_mt => this%mtl
  i = 1
  do
    if (.not. associated(cur_mt)) exit
    int_array(i) = i - 1
    i = i + 1
    cur_mt => cur_mt%next
  enddo

  ! gather relevant information from all processes
  call MPI_Allreduce(local_stride, stride, ONE_INTEGER_MPI, &
                     MPI_INTEGER, MPI_MAX, this%option%mycomm, ierr)
  call MPI_Allreduce(n_mt_local, n_mt_global, ONE_INTEGER_MPI, &
                     MPI_INTEGER, MPI_SUM, this%option%mycomm, ierr)   

  ! create MPI vector and sequential vector for mapping
  call VecCreateMPI(this%option%mycomm, n_mt_local*stride, n_mt_global*stride, & 
                    global_mt_vec,ierr); CHKERRQ(ierr)
  call VecCreateSeq(PETSC_COMM_SELF, n_mt_local*stride, local_mt_vec, ierr); &
         CHKERRQ(ierr)
  call VecSetBlockSize(global_mt_vec, stride, ierr); CHKERRQ(ierr)
  call VecSetBlockSize(local_mt_vec, stride, ierr); CHKERRQ(ierr)

  allocate(check_vars(stride))
  allocate(indices(stride))

  ! collect data for checkpointing
  j = 1
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit

    check_vars(1) = cur_mt%num_aux ! checkpoint #1

    i = n_check_vars + 1
    do
      if (i > stride) exit
      check_vars(i) = -9999
      i = i + 1
    enddo

    do i = 1, stride
      indices(i) = (j - 1)*stride + i - 1
    enddo
    j = j + 1

    call VecSetValues(local_mt_vec, stride, indices, check_vars, &
                     INSERT_VALUES, ierr); CHKERRQ(ierr)

    cur_mt => cur_mt%next

  enddo

  !Create map and add values from the sequential vector to the global 
  call ISCreateBlock(this%option%mycomm, stride, n_mt_local, int_array, &
                     PETSC_COPY_VALUES, is, ierr); CHKERRQ(ierr)
  call VecScatterCreate(local_mt_vec, PETSC_NULL_IS, global_mt_vec, &
                        is, scatter_ctx, ierr); CHKERRQ(ierr)
  call VecScatterBegin(scatter_ctx, local_mt_vec, global_mt_vec, &  
                       INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx, local_mt_vec, global_mt_vec, & 
                     INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)

  ! write the checkpoint file
  dataset_name='material transform model info'
  call VecView(global_mt_vec,viewer,ierr); CHKERRQ(ierr)   
  call VecScatterDestroy(scatter_ctx, ierr); CHKERRQ(ierr)
  call ISDestroy(is, ierr); CHKERRQ(ierr)
  call VecDestroy(global_mt_vec, ierr); CHKERRQ(ierr)
  call VecDestroy(local_mt_vec, ierr); CHKERRQ(ierr)

  ! checkpoint the auxiliary variables
  check_il = PETSC_FALSE
  check_be = PETSC_FALSE
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit
    if (associated(cur_mt%illitization)) then
      check_il = PETSC_TRUE
    endif
    if (associated(cur_mt%buffer_erosion)) then
      check_be = PETSC_TRUE
    endif
    cur_mt => cur_mt%next
  enddo

  if (check_il) then
    global_vec = PETSC_NULL_VEC
    call DiscretizationCreateVector(this%realization%discretization, ONEDOF, &
                                    global_vec, GLOBAL, option)
    call MaterialTransformGetAuxVarVecLoc(this%realization%patch%aux%MT, &
                                          field%work_loc, SMECTITE, &
                                          ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization, field%work_loc, &
                                     global_vec, ONEDOF)
    call VecView(global_vec, viewer, ierr); CHKERRQ(ierr)


    call VecDestroy(global_vec, ierr); CHKERRQ(ierr)
  endif
  ! if (check_be) then
  ! endif

end subroutine PMMaterialTransformCheckpointBinary

! ************************************************************************** !

subroutine PMMaterialTransformRestartBinary(this, viewer)
  !
  ! Restarts data associated with material transform process model
  ! 
  ! Author: Alex Salazar III
  ! Date: 01/20/2022

  use petscvec
  use Option_module
  use Field_module
  use Discretization_module
  use Variables_module, only: SMECTITE

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! viewer: petsc viewer variable
  ! ----------------------------------
  PetscViewer :: viewer
  class(pm_material_transform_type) :: this
  ! ----------------------------------
  ! LOCAL VARIABLES:
  ! ================
  ! is: PETSc index set
  ! scatter_ctx: copies an MPI vector to sequential vectors on all MPI ranks
  ! local_mt_vec: local vector of material transform checkpoint values
  ! global_mt_vec: global vector of material transform checkpoint values
  ! ierr: I/O status indicator
  ! local_stride: local number of vector elements between start and end of block
  ! local_stride_tmp: temporary counter for local stride
  ! stride: global number of vector elements between start and end of block
  ! n_mt_local: number of material transforms on process
  ! n_mt_global: number of material transforms globally
  ! n_check_vars: number of values to checkpoint
  ! i: iterator
  ! indices: indices of the local material transform vector
  ! int_array: keeps track of the material transform number
  ! check_vars: array of checkpointed values
  ! local_mt_array: data converted into a Fortran array
  ! cur_mt: material transform object
  ! dataset_name: descriptor of the material transform checkpoint data
  ! global_vec: global discretization PETSc vector 
  ! natural_vec: local discretization PETSc vector 
  ! check_il: logical check for presence of illitization functions in the
  !   material transform objects so auxiliary variables can be checkpointed
  ! check_be: logical check for presence of buffer erosion models in the
  !   material transform objects so auxiliary variables can be checkpointed
  ! option: option object
  ! field: field object
  ! discretization: discretization object
  ! ----------------------------------
  IS :: is
  VecScatter :: scatter_ctx
  Vec :: local_mt_vec
  Vec :: global_mt_vec
  PetscErrorCode :: ierr
  PetscInt :: local_stride
  PetscInt :: n_mt_local
  PetscInt :: n_mt_global
  PetscInt :: n_check_vars
  PetscInt :: local_stride_tmp
  PetscInt :: i
  PetscInt :: stride
  PetscInt, allocatable :: indices(:)
  PetscInt, allocatable :: int_array(:)
  PetscReal, allocatable :: check_vars(:)
  PetscReal, pointer :: local_mt_array(:)
  class(material_transform_type), pointer :: cur_mt
  character(len=MAXSTRINGLENGTH) :: dataset_name
  Vec :: global_vec
  Vec :: natural_vec
  PetscBool :: check_il
  PetscBool :: check_be
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  ! ----------------------------------

  option => this%realization%option
  field => this%realization%field
  discretization => this%realization%discretization

  local_stride = 0
  local_stride_tmp = 0
  n_mt_local = 0
  n_mt_global = 0

  ! current information in PM Material Transform that is checkpointed:
  !   (1) num_aux
  n_check_vars = 1 !number of scalar checkpoint variables

  cur_mt => this%mtl
  do 
    if (.not. associated(cur_mt)) exit
    n_mt_local = n_mt_local + 1
    local_stride_tmp = local_stride_tmp + n_check_vars
    cur_mt => cur_mt%next
    if (local_stride_tmp > local_stride) then
      local_stride = local_stride_tmp
    endif
    local_stride_tmp = 0
  enddo

  allocate(int_array(n_mt_local))
  i = 1
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit
    int_array(i) = i - 1
    i = i + 1
    cur_mt => cur_mt%next
  enddo

  ! gather relevant information from all processes
  call MPI_Allreduce(local_stride, stride, ONE_INTEGER_MPI, &
                     MPI_INTEGER, MPI_MAX, this%option%mycomm, ierr)
  call MPI_Allreduce(n_mt_local, n_mt_global, ONE_INTEGER_MPI, &
                     MPI_INTEGER, MPI_SUM, this%option%mycomm, ierr)   

  ! create MPI vector for HDF5 reading and sequential vector for mt information
  !   stored in the process
  call VecCreateMPI(this%option%mycomm, n_mt_local*stride, n_mt_global*stride, & 
                    global_mt_vec,ierr); CHKERRQ(ierr)
  call VecCreateSeq(PETSC_COMM_SELF, n_mt_local*stride, local_mt_vec, ierr); &
         CHKERRQ(ierr)
  call VecSetBlockSize(global_mt_vec, stride, ierr); CHKERRQ(ierr)
  call VecSetBlockSize(local_mt_vec, stride, ierr); CHKERRQ(ierr)

  ! read data from HDF5
  call VecLoad(global_mt_vec, viewer, ierr);CHKERRQ(ierr)

  ! create mapping between MPI and sequential vectors
  call ISCreateBlock(this%option%mycomm, stride, n_mt_local, int_array, &
                     PETSC_COPY_VALUES, is, ierr); CHKERRQ(ierr)
  call VecScatterCreate(global_mt_vec, is, local_mt_vec, &
                        PETSC_NULL_IS, scatter_ctx, ierr); CHKERRQ(ierr)

  ! obtain data from the MPI vector
  call VecScatterBegin(scatter_ctx, global_mt_vec, local_mt_vec, &  
                       INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
  call VecScatterEnd(scatter_ctx, global_mt_vec, local_mt_vec, & 
                     INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)

  ! convert the data into a Fortran array
  call VecGetArrayF90(local_mt_vec, local_mt_array, ierr); CHKERRQ(ierr)

  ! assign checkpointed material transform information
  i = 1
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit

    cur_mt%num_aux = local_mt_array(i) ! checkpoint #1

    cur_mt => cur_mt%next
    i = i + stride
  enddo

  call VecRestoreArrayF90(local_mt_vec, local_mt_array, ierr); CHKERRQ(ierr)
  call VecScatterDestroy(scatter_ctx, ierr); CHKERRQ(ierr)
  call ISDestroy(is, ierr); CHKERRQ(ierr)
  call VecDestroy(global_mt_vec, ierr); CHKERRQ(ierr)
  call VecDestroy(local_mt_vec, ierr); CHKERRQ(ierr)

  ! retrieve the auxiliary variables
  check_il = PETSC_FALSE
  check_be = PETSC_FALSE
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit
    if (associated(cur_mt%illitization)) then
      check_il = PETSC_TRUE
    endif
    if (associated(cur_mt%buffer_erosion)) then
      check_be = PETSC_TRUE
    endif
    cur_mt => cur_mt%next
  enddo

  if (check_il) then
    global_vec = PETSC_NULL_VEC
    call DiscretizationCreateVector(this%realization%discretization, ONEDOF, &
                                    global_vec, GLOBAL, option)
    call VecLoad(global_vec, viewer, ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization, global_vec, &
                                     field%work_loc, ONEDOF)
    call MaterialTransformSetAuxVarVecLoc(this%realization%patch%aux%MT, &
                                          field%work_loc, SMECTITE, &
                                          ZERO_INTEGER)
    call VecDestroy(global_vec, ierr); CHKERRQ(ierr)
  endif
  ! if (check_be) then
  ! endif

end subroutine PMMaterialTransformRestartBinary

! ************************************************************************** !

subroutine PMMaterialTransformInputRecord(this)
  !
  ! Writes process model information to the input record file.
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022
  !

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! --------------------------------
  class(pm_material_transform_type) :: this
  ! --------------------------------
  ! LOCAL VARIABLES:
  ! ================
  ! id: number of output unit
  ! --------------------------------
  PetscInt :: id
  ! --------------------------------

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  if (associated(this%mtl)) then
    ! print material transform model information
    call MaterialTransformInputRecord(this%mtl)
  endif

end subroutine PMMaterialTransformInputRecord

! ************************************************************************** !

subroutine PMMaterialTransformStrip(this)
  !
  ! Strips the material transform process model
  !
  ! Author: Alex Salazar III
  ! Date: 03/03/2022

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! --------------------------------
  class(pm_material_transform_type) :: this
  ! --------------------------------
  ! LOCAL VARIABLES:
  ! ================
  ! cur_mt: pointer to current material transform object
  ! prev_mt: pointer to previous material transform object
  ! --------------------------------
  type(material_transform_type), pointer :: cur_mt, prev_mt
  ! --------------------------------

  nullify(this%realization)
  cur_mt => this%mtl
  do
    if (.not. associated(cur_mt)) exit
    prev_mt => cur_mt
    cur_mt => cur_mt%next
    call MaterialTransformDestroy(prev_mt)
  enddo
  deallocate(this%mtl)
  nullify(this%mtl)

end subroutine PMMaterialTransformStrip

! ************************************************************************** !

subroutine PMMaterialTransformDestroy(this)
  !
  ! Destroys auxiliary process model
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022

  implicit none

  ! INPUT ARGUMENTS:
  ! ================
  ! this (input/output): material transform process model object
  ! --------------------------------
  class(pm_material_transform_type) :: this
  ! --------------------------------

  call PMBaseDestroy(this)
  call PMMaterialTransformStrip(this)

end subroutine PMMaterialTransformDestroy

end module PM_Material_Transform_class
