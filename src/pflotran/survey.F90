module Survey_module

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
            SurveyRead, &
            SurveyReadERT, &
            SurveyDestroy
 
contains

! ************************************************************************** !

function SurveyCreate()
  ! 
  ! Creates survey type
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/27/21
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

subroutine SurveyRead(survey,input,option)
  ! 
  ! Read a geophysics survey 
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/29/21

  use Input_Aux_module
  use Option_module
  use String_module  

  implicit none

  type(survey_type), pointer :: survey
  type(input_type), pointer :: input
  type(option_type) :: option

  type(input_type), pointer :: input_tmp
  character(len=MAXWORDLENGTH) :: word

  ! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (input%ierr /= 0) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','SURVEY')
    call StringToUpper(word)
    select case(trim(word))
    case('FILE_NAME')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'FILENAME','SURVEY')
      survey%filename = word 
    case('FORMAT')
      call InputReadCard(input,option,word)
      call InputErrorMsg(input,option,'keyword','GRID')
      call StringToUpper(word)
      select case (trim(word))
      case ('E4D_SRV')
        !survey%format = E4D_SRV  
      case default
        option%io_buffer = 'Only E4D_SRV FORMAT can be ' // &
                       'provided under SURVEY.'
        call PrintErrMsg(option)
      end select
    case default
      call InputKeywordUnrecognized(input,word,'SURVEY',option)
    end select
  enddo
  call InputPopBlock(input,option)  
  
end subroutine SurveyRead

! ************************************************************************** !

subroutine SurveyReadERT(survey,grid,input,option)
  ! 
  ! Read an ERT survey 
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/29/21
  !

  use Option_module
  use Input_Aux_module
  use Grid_module
  use String_module 

  implicit none

  type(survey_type), pointer :: survey
  type(grid_type), pointer :: grid  
  type(input_type), pointer :: input
  type(option_type), pointer :: option

  PetscInt :: ielec, idata
  PetscInt :: itemp
  character(len=MAXWORDLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: string_ielec,string_idata

  error_string = "SURVEY from file "//trim(survey%filename)

  call InputReadPflotranString(input,option)
  call InputReadInt(input,option,survey%num_electrode)
  call InputErrorMsg(input,option,'num_electrode',error_string)

  allocate(survey%pos_electrode(THREE_INTEGER,survey%num_electrode))
  survey%pos_electrode = 0.d0
  allocate(survey%flag_electrode(survey%num_electrode))
  survey%flag_electrode = 1

  do ielec=1,survey%num_electrode
    write(string_ielec,*) ielec
    string_ielec = trim(adjustl(string_ielec)) 
    call InputReadPflotranString(input,option)
    call InputReadInt(input,option,itemp)
    call InputReadDouble(input,option,survey%pos_electrode(X_DIRECTION,ielec))
    call InputErrorMsg(input,option,'x-position of electrode ' &
                                     //string_ielec,error_string)
    call InputReadDouble(input,option,survey%pos_electrode(Y_DIRECTION,ielec))
    call InputErrorMsg(input,option,'y-position of electrode ' &
                                     //string_ielec,error_string)
    call InputReadDouble(input,option,survey%pos_electrode(Z_DIRECTION,ielec))
    call InputErrorMsg(input,option,'z-position of electrode ' &
                                     //string_ielec,error_string)
    call InputReadInt(input,option,survey%flag_electrode(ielec))
    call InputErrorMsg(input,option,'flag for electrode ' &
                                     //string_ielec,error_string)
  enddo

  call InputReadPflotranString(input,option)
  call InputReadInt(input,option,survey%num_measurement)
  call InputErrorMsg(input,option,'num_measurement',error_string)

  allocate(survey%config(FOUR_INTEGER,survey%num_measurement))
  allocate(survey%dsim(survey%num_measurement))
  allocate(survey%Wd(survey%num_measurement))
  survey%config = 0
  survey%dsim = 0.d0
  survey%Wd = 0.d0

  do idata=1,survey%num_measurement
    write(string_idata,*) idata 
    string_idata = trim(adjustl(string_idata)) 
    call InputReadPflotranString(input,option)
    call InputReadInt(input,option,itemp)
    call InputReadInt(input,option,survey%config(ONE_INTEGER,idata))
    call InputErrorMsg(input,option,'A electrode configuration for data ' &
                                    //string_idata,error_string)
    call InputReadInt(input,option,survey%config(TWO_INTEGER,idata))
    call InputErrorMsg(input,option,'B electrode configuration for data ' &
                                    //string_idata,error_string)
    call InputReadInt(input,option,survey%config(THREE_INTEGER,idata))
    call InputErrorMsg(input,option,'M electrode configuration for data ' &
                                    //string_idata,error_string)
    call InputReadInt(input,option,survey%config(FOUR_INTEGER,idata))
    call InputErrorMsg(input,option,'N electrode configuration for data ' &
                                    //string_idata,error_string)
    call InputReadDouble(input,option,survey%dsim(idata))
    call InputErrorMsg(input,option,'data dobs for data ' &
                                    //string_idata,error_string)
    call InputReadDouble(input,option,survey%Wd(idata))
    call InputErrorMsg(input,option,'weight Wd for data ' &
                                    //string_idata,error_string)
  enddo

  ! Get cell ids corrsponding to electrode positions
  allocate(survey%ipos_electrode(survey%num_electrode))
  survey%ipos_electrode = 1

  call SurveyGetElectrodeIndexFromPos(survey,grid,option)

end subroutine SurveyReadERT

! ************************************************************************** !

subroutine SurveyGetElectrodeIndexFromPos(survey,grid,option)
  ! 
  ! Get the local-ids corrspodning to each elctrode location
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 02/01/21
  !

  use Grid_module
  use Geometry_module  
  use Option_module

  implicit none

  type(survey_type), pointer :: survey
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option

  type(point3d_type) :: coordinate  
  PetscInt :: ielec, local_id

  do ielec=1,survey%num_electrode
    coordinate%x = survey%pos_electrode(X_DIRECTION,ielec)
    coordinate%y = survey%pos_electrode(Y_DIRECTION,ielec)
    coordinate%z = survey%pos_electrode(Z_DIRECTION,ielec)

    call GridGetLocalIDFromCoordinate(grid,coordinate,option,local_id)
    survey%ipos_electrode(ielec) = local_id
  enddo

end subroutine SurveyGetElectrodeIndexFromPos

! ************************************************************************** !

subroutine SurveyDestroy(survey)
  ! 
  ! Deallocates survey
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/28/21
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

end module Survey_module