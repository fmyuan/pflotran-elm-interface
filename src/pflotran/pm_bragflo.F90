module PM_Bragflo_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_WIPP_Flow_class 
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_wippflo_type) :: pm_bragflo_type
    character(len=MAXWORDLENGTH) :: alpha_dataset_name
    character(len=MAXWORDLENGTH) :: elevation_dataset_name
  contains
    procedure, public :: Read => PMBragfloRead
    procedure, public :: InitializeRun => PMBragfloInitializeRun
    procedure, public :: Residual => PMBragfloResidual
    procedure, public :: Jacobian => PMBragfloJacobian
    procedure, public :: InputRecord => PMBragfloInputRecord
  end type pm_bragflo_type
  
  public :: PMBragfloCreate
  
contains

! ************************************************************************** !

function PMBragfloCreate()
  ! 
  ! Creates Bragflo process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/17
  ! 
  use PM_WIPP_Flow_class
  use WIPP_Flow_Aux_module
  use Bragflo_module

  implicit none
  
  class(pm_bragflo_type), pointer :: PMBragfloCreate

  class(pm_bragflo_type), pointer :: bragflo_pm

  wippflo_use_bragflo_flux = PETSC_TRUE
  
  allocate(bragflo_pm)
  call PMWIPPFloInitObject(bragflo_pm)

  bragflo_pm%name = 'BRAGFLO'
  bragflo_pm%alpha_dataset_name = ''
  bragflo_pm%elevation_dataset_name = ''

  PMBragfloCreate => bragflo_pm
  
end function PMBragfloCreate

! ************************************************************************** !

subroutine PMBragfloRead(this,input)
  ! 
  ! Reads BRAGFLO PM options block
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  !
  use PM_Subsurface_Flow_class
  use WIPP_Flow_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword, word
  class(pm_bragflo_type) :: this
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'BRAGFLO Options'

  input%ierr = 0
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMWIPPFloReadSelectCase(this,input,keyword,found, &
                                 error_string,option)
    if (found) cycle

    select case(keyword)
      case('ALPHA_DATASET')
        call InputReadNChars(input,option,this%alpha_dataset_name,&
                             MAXWORDLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'ALPHA DATASET,NAME',error_string)
      case('ELEVATION_DATASET')
        call InputReadNChars(input,option,this%elevation_dataset_name,&
                             MAXWORDLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'ELEVATION DATASET,NAME',error_string)
      case default
        call InputKeywordUnrecognized(keyword,'BRAGFLO Mode',option)
    end select
  enddo

end subroutine PMBragfloRead

! ************************************************************************** !

subroutine PMBragfloInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/17
  ! 
  
  implicit none
  
  class(pm_bragflo_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'BRAGFLO'

end subroutine PMBragfloInputRecord

! ************************************************************************** !

recursive subroutine PMBragfloInitializeRun(this)
  ! 
  ! Initializes the Bragflo mode run.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17

  use Realization_Base_class
  use WIPP_Flow_Aux_module
  use Input_Aux_module
  use Field_module
  use Dataset_Base_class
  use Dataset_Common_HDF5_class
  use Dataset_module
  use HDF5_module
  use Option_module

  implicit none

  class(pm_bragflo_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr
  type(input_type), pointer :: input
  class(dataset_base_type), pointer :: dataset
  type(field_type), pointer :: field
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  character(len=MAXSTRINGLENGTH) :: block_string
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscReal, pointer :: work_loc_p(:)
  PetscInt :: idof, ghosted_id


  call PMWIPPFloInitializeRun(this)

  ! read in alphas
  if (len_trim(this%alpha_dataset_name) > 0) then
    field => this%realization%field
    string = 'BRAGFLO ALPHA Dataset'
    dataset => DatasetBaseGetPointer(this%realization%datasets, &
                                     this%alpha_dataset_name, &
                                     string,this%option)
    select type(d => dataset)
      class is(dataset_common_hdf5_type)
        string2 = ''
        string = d%hdf5_dataset_name
        call HDF5ReadCellIndexedRealArray(this%realization, &
                                          field%work, &
                                          d%filename,&
                                          string2, &
                                          string, &
                                          d%realization_dependent)
        wippflo_auxvars => this%realization%patch%aux%WIPPFlo%auxvars
        call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
        call VecGetArrayReadF90(this%realization%field%work_loc,work_loc_p, &
                                ierr);CHKERRQ(ierr)
        do ghosted_id = 1, this%realization%patch%grid%ngmax
          do idof = 0, this%option%nflowdof
            wippflo_auxvars(idof,ghosted_id)%alpha = work_loc_p(ghosted_id)
          enddo
        enddo
        call VecRestoreArrayReadF90(this%realization%field%work_loc, &
                                    work_loc_p,ierr);CHKERRQ(ierr)
      class default
        this%option%io_buffer = 'Unsupported dataset type for BRAGFLO ALPHA.'
        call printErrMsg(this%option)
    end select
  else
    this%option%io_buffer = 'ALPHA should have been read from a dataset.'
    call printErrMsg(this%option)
  endif

  ! read in elevations
  if (len_trim(this%elevation_dataset_name) > 0) then
    field => this%realization%field
    string = 'BRAGFLO Elevation Dataset'
    dataset => DatasetBaseGetPointer(this%realization%datasets, &
                                     this%elevation_dataset_name, &
                                     string,this%option)
    select type(d => dataset)
      class is(dataset_common_hdf5_type)
        string2 = ''
        string = d%hdf5_dataset_name
        call HDF5ReadCellIndexedRealArray(this%realization, &
                                          field%work, &
                                          d%filename,&
                                          string2, &
                                          string, &
                                          d%realization_dependent)
        wippflo_auxvars => this%realization%patch%aux%WIPPFlo%auxvars
        call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
        call VecGetArrayReadF90(this%realization%field%work_loc,work_loc_p, &
                                ierr);CHKERRQ(ierr)
        do ghosted_id = 1, this%realization%patch%grid%ngmax
          do idof = 0, this%option%nflowdof
            wippflo_auxvars(idof,ghosted_id)%elevation = work_loc_p(ghosted_id)
          enddo
        enddo
        call VecRestoreArrayReadF90(this%realization%field%work_loc, &
                                    work_loc_p,ierr);CHKERRQ(ierr)
      class default
        this%option%io_buffer = 'Unsupported dataset type for BRAGFLO &
          &Elevation.'
        call printErrMsg(this%option)
    end select
  endif

end subroutine PMBragfloInitializeRun

! ************************************************************************** !

subroutine PMBragfloResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use PM_Subsurface_Flow_class
  use Bragflo_module, only : BragfloResidual

  implicit none

  class(pm_bragflo_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  call PMSubsurfaceFlowUpdatePropertiesNI(this)

  ! calculate residual
  call BragfloResidual(snes,xx,r,this%realization,this%pmwss_ptr,ierr)

  call this%PostSolve()

end subroutine PMBragfloResidual

! ************************************************************************** !

subroutine PMBragfloJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Bragflo_module, only : BragfloJacobian

  implicit none

  class(pm_bragflo_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  call BragfloJacobian(snes,xx,A,B,this%realization,this%pmwss_ptr,ierr)

end subroutine PMBragfloJacobian

end module PM_Bragflo_class
