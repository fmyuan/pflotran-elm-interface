module PM_Material_Transform_class

! MODULE DESCRIPTION:
! ===========================================================================
! This process model incorporates material transformations as surrogate models
! for physical phenomena given evolving system conditions
! ===========================================================================

#include "petsc/finclude/petscsys.h"
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
    class(material_transform_type), pointer :: mtl
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
    ! procedure, public :: CheckpointHDF5 => PMMaterialTransformCheckpointHDF5
    ! procedure, public :: CheckpointBinary => PMMaterialTransformCheckpointBinary
    ! procedure, public :: RestartHDF5 => PMMaterialTransformRestartHDF5
    ! procedure, public :: RestartBinary => PMMaterialTransformRestartBinary
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

subroutine PMMaterialTransformSetRealization(this,realization)
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022

  implicit none

! INPUT ARGUMENTS:
! ================
! this (input/output): material transform process model object
! realization (input): subsurface realization object
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

    allocate(mtf)

    mtf => this%mtl

    call MaterialTransformAddToList(mtf,this%realization%material_transform)

    nullify(mtf)

  endif

  ! set up mapping for material transform functions
  patch%material_transform => this%realization%material_transform
  call MaterialTransformConvertListToArray(patch%material_transform, &
                                           patch%material_transform_array, &
                                           option)

  ! material property mapping to PM Material Transform
  cur_material_property => this%realization%material_properties

  do
    if (.not.associated(cur_material_property)) exit
    
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
  allocate(MT_auxvars(grid%nlmax))
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    material_id = patch%imat(ghosted_id)
    
    call MaterialTransformAuxVarInit(MT_auxvars(ghosted_id))
    
    if (Initialized(patch%mtf_id(ghosted_id))) then
      ! pointer to material transform in patch ghosted id
      mtf => patch%material_transform_array(patch%mtf_id(ghosted_id))%ptr
      
      if (associated(mtf%illitization)) then
        MT_auxvars(ghosted_id)%il_aux => IllitizationAuxVarInit(option)
      endif

      if (associated(mtf%buffer_erosion)) then
        MT_auxvars(ghosted_id)%be_aux => BufferErosionAuxVarInit()
      endif
      
    endif
    
  enddo
  patch%aux%MT%auxvars => MT_auxvars
  patch%aux%MT%num_aux = grid%nlmax

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
    
    material_aux => material_auxvars(ghosted_id)
    MT_aux => MT_auxvars(ghosted_id)
    
    if (Initialized(patch%mtf_id(ghosted_id))) then
      ! pointer to material transform in patch ghosted id
      if (associated(patch%material_transform_array)) then
        allocate(mtf)
        mtf => patch%material_transform_array(patch%mtf_id(ghosted_id))%ptr
        if (associated(MT_aux%il_aux)) then
          MT_aux%il_aux%fs0 = &
            mtf%illitization%illitization_function%ilt_fs0
          MT_aux%il_aux%fs = &
            mtf%illitization%illitization_function%ilt_fs0
        endif
        nullify(mtf)
      endif
      
    endif
    
  enddo

end subroutine PMMaterialTransformInitializeRun

! ************************************************************************** !

subroutine PMMaterialTransformInitializeTS(this)
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

subroutine PMMaterialTransformSolve(this,time,ierr)
  !
  ! Updates the soil composition based on the model chosen
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
    
    global_aux => global_auxvars(ghosted_id)
    material_aux => material_auxvars(ghosted_id)
    MT_aux => MT_auxvars(ghosted_id)
    
    if (Initialized(patch%mtf_id(ghosted_id))) then
      ! pointer to material transform in patch ghosted id
      mtf => patch%material_transform_array(patch%mtf_id(ghosted_id))%ptr
      
      if (associated(mtf%illitization)) then
        call mtf%illitization%illitization_function%CalculateILT( &
               MT_aux%il_aux%fs, &
               global_aux%temp, &
               time, &
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
    
  enddo


  ierr = 0

end subroutine PMMaterialTransformSolve

! ***************************************************************************** !

subroutine PMMaterialTransformCheckpointHDF5(this,pm_grp_id)
  ! 
  ! Checkpoints data associated with the material transform process model
  !
  ! Author: Alex Salazar III
  ! Date: 01/20/2022
  ! 

  use Option_module
  use Realization_Subsurface_class
  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec

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
!
! ----------------------------------

! ----------------------------------
  
end subroutine PMMaterialTransformCheckpointHDF5

! ***************************************************************************** !


subroutine PMMaterialTransformRestartHDF5(this,pm_grp_id)
  ! 
  ! Restarts data associated with material transform process model
  ! 
  ! Author: Alex Salazar III
  ! Date: 01/20/2022

  use Option_module
  use Realization_Subsurface_class
  use hdf5
  use HDF5_module, only : HDF5ReadDataSetInVec

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
!
! ----------------------------------

! ----------------------------------

  
end subroutine PMMaterialTransformRestartHDF5

! ************************************************************************** !

subroutine PMMaterialTransformCheckpointBinary(this,viewer)
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
!
! ----------------------------------

! ----------------------------------

end subroutine PMMaterialTransformCheckpointBinary

! ***************************************************************************** !

subroutine PMMaterialTransformRestartBinary(this,viewer)
  !
  ! Restarts data associated with material transform process model
  ! 
  ! Author: Alex Salazar III
  ! Date: 01/20/2022

  use Option_module
  use hdf5
  use HDF5_module, only : HDF5ReadDataSetInVec

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
!
! ----------------------------------

! ----------------------------------

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
    if (.not.associated(cur_mt)) exit
    prev_mt => cur_mt
    cur_mt => cur_mt%next
    call MaterialTransformDestroy(prev_mt)
  enddo
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
