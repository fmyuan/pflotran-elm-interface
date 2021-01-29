module survey_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: survey_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: filename

    PetscInt :: num_electrode                   ! Number of electrodes
    PetscInt, pointer :: ipos_electrode(:)      ! cell id of electrode pos
    PetscInt, pointer :: flag_electrode(:)      ! 0-> below/1-> on surface 
    PetscReal, pointer :: pos_electrode(:,:)    ! electrode positions

    
    PetscInt :: num_measurement                 ! number of data
    PetscInt, pointer :: config(:,:)            ! survey configuration
    PetscReal, pointer :: dsim(:)               ! Simulated data
    PetscReal, pointer :: dobs(:)               ! Observed data
    PetscReal, pointer :: Wd(:)                 ! data weight

  end type survey_type


  public :: SurveyCreate, &
            SurveyReadFile, &
            SurveyDestroy
 
contains

! ************************************************************************** !

function SurveyCreate()
  ! 
  ! Creates survey type
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/27/2021
  ! 
  implicit none
  
  type(survey_type), pointer :: SurveyCreate
  type(survey_type), pointer :: survey
  
  allocate(survey)

  survey%name = ''
  survey%filename = ''  
  survey%num_electrode = 0
  survey%num_measurement = 0

  nullify(survey%ipos_electrode)
  nullify(survey%pos_electrode)
  nullify(survey%flag_electrode)
  nullify(survey%config)
  nullify(survey%dsim)
  nullify(survey%dobs)
  nullify(survey%Wd)
   
  SurveyCreate => survey 
  
end function SurveyCreate

! ************************************************************************** !

subroutine SurveyReadFile(survey)
  ! 
  ! Read a geophysics survey 
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/29/2021

  use Input_Aux_module
  use Option_module

  implicit none

  type(survey_type), pointer :: survey
  type(option_type) :: option

  type(input_type), pointer :: input
  
  input => InputCreate(IUNIT_TEMP,survey%filename,option)

  if (option%igeopmode==ERT_MODE) then
    call SurveyReadERT(survey,input,option)
  else
    option%io_buffer = 'Only ERT mode is supported for SurveyReadFile.'
    call PrintErrMsg(option)
  endif  

  call InputDestroy(input)

end subroutine SurveyReadFile

! ************************************************************************** !

subroutine SurveyReadERT(survey,input,option)
  ! 
  ! Read an ERT survey 
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/29/2021

  use Option_module
  use Input_Aux_module

  implicit none

  type(survey_type), pointer :: survey
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: ielec, idata
  PetscInt :: itemp

  call InputReadPflotranString(input,option)
  call InputReadInt(input,option,survey%num_electrode)
  call InputErrorMsg(input,option,'num_electrode','ERT SURVEY FILE')

  allocate(survey%pos_electrode(THREE_INTEGER,survey%num_electrode))
  survey%pos_electrode = 0.d0
  allocate(survey%flag_electrode(survey%num_electrode))
  survey%flag_electrode = 1

  do ielec=1,survey%num_electrode
    call InputReadPflotranString(input,option)
    call InputReadInt(input,option,itemp)
    call InputReadDouble(input,option,survey%pos_electrode(X_DIRECTION,ielec))
    call InputErrorMsg(input,option,'x-position of electrode','ERT SURVEY FILE')
    call InputReadDouble(input,option,survey%pos_electrode(Y_DIRECTION,ielec))
    call InputErrorMsg(input,option,'y-position of electrode','ERT SURVEY FILE')
    call InputReadDouble(input,option,survey%pos_electrode(Z_DIRECTION,ielec))
    call InputErrorMsg(input,option,'z-position of electrode','ERT SURVEY FILE')
    call InputReadInt(input,option,survey%flag_electrode(ielec))
    call InputErrorMsg(input,option,'flag for electrode','ERT SURVEY FILE')
  enddo

  call InputReadPflotranString(input,option)
  call InputReadInt(input,option,survey%num_measurement)
  call InputErrorMsg(input,option,'num_measurement','ERT SURVEY FILE')

  allocate(survey%config(FOUR_INTEGER,survey%num_measurement))
  allocate(survey%dsim(survey%num_measurement))
  allocate(survey%Wd(survey%num_measurement))
  survey%config = 0
  survey%dsim = 0.d0
  survey%Wd = 0.d0

  do idata=1,survey%num_measurement
    call InputReadPflotranString(input,option)
    call InputReadInt(input,option,itemp)
    call InputReadInt(input,option,survey%config(ONE_INTEGER,idata))
    call InputErrorMsg(input,option,'A electrode configuration','ERT SURVEY FILE')
    call InputReadInt(input,option,survey%config(TWO_INTEGER,idata))
    call InputErrorMsg(input,option,'B electrode configuration','ERT SURVEY FILE')
    call InputReadInt(input,option,survey%config(THREE_INTEGER,idata))
    call InputErrorMsg(input,option,'M electrode configuration','ERT SURVEY FILE')
    call InputReadInt(input,option,survey%config(FOUR_INTEGER,idata))
    call InputErrorMsg(input,option,'N electrode configuration','ERT SURVEY FILE')
    call InputReadDouble(input,option,survey%dsim(idata))
    call InputErrorMsg(input,option,'data dobs','ERT SURVEY FILE')
    call InputReadDouble(input,option,survey%Wd(idata))
    call InputErrorMsg(input,option,'weight Wd','ERT SURVEY FILE')
  enddo

end subroutine SurveyReadERT

! ************************************************************************** !

subroutine SurveyDestroy(survey)
  ! 
  ! Deallocates survey
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/28/2021
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  type(survey_type), pointer :: survey
    
  if (.not.associated(survey)) return
  
  call DeallocateArray(survey%ipos_electrode)
  call DeallocateArray(survey%pos_electrode)
  call DeallocateArray(survey%flag_electrode)
  call DeallocateArray(survey%config)
  
  call DeallocateArray(survey%dsim)
  call DeallocateArray(survey%dobs)
  call DeallocateArray(survey%Wd)
  
  deallocate(survey)
  nullify(survey)

end subroutine SurveyDestroy

end module survey_module