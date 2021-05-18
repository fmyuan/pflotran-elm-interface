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

    PetscReal :: apparent_conductivity          ! app cond for an ERT survey
    PetscReal :: average_conductivity           ! avg cond of given cond model

  end type survey_type


  public :: SurveyCreate, &
            SurveyRead, &
            SurveyReadERT, &
            SurveyWriteERT, &
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
  survey%num_electrode = UNINITIALIZED_INTEGER
  survey%num_measurement = UNINITIALIZED_INTEGER
  survey%apparent_conductivity = UNINITIALIZED_DOUBLE
  survey%average_conductivity = UNINITIALIZED_DOUBLE

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
  use Option_Geophysics_module
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
  option%geophysics%num_electrodes = survey%num_electrode

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
  allocate(survey%dobs(survey%num_measurement))
  allocate(survey%Wd(survey%num_measurement))
  survey%config = 0
  survey%dsim = 0.d0
  survey%dobs = 0.d0
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
    call InputReadDouble(input,option,survey%dobs(idata))
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

subroutine SurveyWriteERT(survey,time_suffix,option)
  !
  ! writes simulated ERT data in a .srv file with
  ! filename = prefix(survey%filename)//-simulated.srv
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/08/21
  !
  use Option_module

  implicit none

  type(survey_type), pointer :: survey
  character(len=MAXWORDLENGTH) :: time_suffix
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: iprefix,i
  PetscInt :: fid

  filename = trim(option%global_prefix) // trim(option%group_prefix) // &
             '-ert-' // trim(time_suffix) // '.srv'
  fid = IUNIT_TEMP
  open(fid,file=filename,status='replace',action='write')
  write(word,*) survey%num_electrode
  write(fid,'(a)',advance='no') trim(adjustl(word))
  write(fid,'(15x,a)',advance="yes") "number of electrodes"
  do i=1,survey%num_electrode
     write(fid,"(i10,3f15.5,i10)") i,survey%pos_electrode(1:3,i), &
                                  survey%flag_electrode(i)
  end do

  write(fid,*)
  write(word,*) survey%num_measurement
  write(fid,'(a)',advance='no') trim(adjustl(word))
  write(fid,'(15x,a)',advance="yes") "number of measurements"
  do i=1,survey%num_measurement
    write(fid,"(i10,4i10,2es15.5)") i,survey%config(1:4,i), &
                  survey%dsim(i),0.05*abs(survey%dsim(i))+0.01
  end do

  close(fid)

end subroutine SurveyWriteERT

! ************************************************************************** !

subroutine SurveyCalculateApparentCond(survey)
  !
  ! Computes apparent conductivity for an ERT survey data
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/04/21
  !

  implicit none

  type(survey_type), pointer :: survey

  PetscReal :: cond,max_cond,min_cond
  PetscInt :: idata,ia,ib,im,in
  ! Geometric Factors
  PetscReal :: GFam,GFan,GFbm,GFbn,GF
  PetscReal :: A_pos(3),B_pos(3),M_pos(3),N_pos(3)
  PetscReal :: sum_cond,sum_weight

  max_cond = 0.d0
  min_cond = 1.d30
  sum_cond = 0.d0
  sum_weight = 0.d0

  do idata=1,survey%num_measurement
    ! for A and B electrodes
    ia = survey%config(1,idata)
    ib = survey%config(2,idata)
    im = survey%config(3,idata)
    in = survey%config(4,idata)

    if (ia/=0) A_pos = survey%pos_electrode(:,ia)
    if (ib/=0) B_pos = survey%pos_electrode(:,ib)
    if (im/=0) M_pos = survey%pos_electrode(:,im)
    if (in/=0) N_pos = survey%pos_electrode(:,in)

    GFam = 0.d0; GFan = 0.d0; GFbm = 0.d0; GFbn = 0.d0

    if (ia/=0 .and. im/=0) GFam = GeometricFactor(ia,im,A_pos,M_pos)
    if (ia/=0 .and. in/=0) GFan = GeometricFactor(ia,in,A_pos,N_pos)
    if (ib/=0 .and. im/=0) GFbm = GeometricFactor(ib,im,B_pos,M_pos)
    if (ib/=0 .and. in/=0) GFbn = GeometricFactor(ib,in,B_pos,N_pos)

    ! Total Geometric Factor
    GF = (GFam - GFan - GFbm + GFbn) / 2*PI

    ! TODO: CHANGE dsim to dobs
    if (survey%dobs(idata) /= 0) cond = GF/survey%dobs(idata)
    if (cond > 0) then
      ! finds max and min conductivity
      if (cond > max_cond) max_cond = cond
      if (cond < min_cond) min_cond = cond
      if (survey%Wd(idata) /= 0) then
        sum_cond = sum_cond + cond/survey%Wd(idata)
        sum_weight = sum_weight + 1/survey%Wd(idata)
      endif
    endif

  enddo

  survey%apparent_conductivity = sum_cond / sum_weight

  ! TODO (pj): write max and min apparent cond info
  ! and

  contains

  function GeometricFactor(ip,iq,p_pos,q_pos)
    implicit none

    PetscInt :: ip,iq
    PetscReal :: p_pos(3),q_pos(3)
    PetscReal :: GeometricFactor
    PetscReal :: dist

    GeometricFactor = 0.d0
    if (ip==0 .or. iq==0) return

    ! add small number to avoid overshooting for dist=0
    dist = NORM2(p_pos - q_pos) + 1.d-15

    GeometricFactor = 1 / dist

  end function GeometricFactor

end subroutine SurveyCalculateApparentCond

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
