module Output_Observation_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Logging_module 
  use Output_Aux_module
  use Output_Common_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  ! flags signifying the first time a routine is called during a given
  ! simulation
  PetscBool :: observation_first
  PetscBool :: observation_aggregate_first
  PetscBool :: check_for_obs_points
  PetscBool :: calculate_velocities ! true if any obs. pt. prints velocity
  PetscBool :: secondary_observation_first
  PetscBool :: secondary_check_for_obs_points
  PetscBool :: mass_balance_first
  PetscBool :: integral_flux_first
  PetscInt  :: ewriter_summ_count
  PetscInt  :: ewriter_rest_count
  PetscInt  :: linerept_count

  public :: OutputObservation, &
            OutputObservationInit, &
            OutputMassBalance, &
            OutputEclipseFiles, &
            OutputLineRept, &
            OutputIntegralFlux
            
contains

! ************************************************************************** !

subroutine OutputObservationInit(num_steps)
  ! 
  ! Initializes module variables for observation
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/13
  ! 

  use Option_module

  implicit none
  
  PetscInt :: num_steps
  
  check_for_obs_points = PETSC_TRUE
  calculate_velocities = PETSC_FALSE
  secondary_check_for_obs_points = PETSC_TRUE
  if (num_steps == 0) then
    observation_first = PETSC_TRUE
    secondary_observation_first = PETSC_TRUE
    mass_balance_first = PETSC_TRUE
    integral_flux_first = PETSC_TRUE
    observation_aggregate_first = PETSC_TRUE
  else
    observation_first = PETSC_FALSE
    secondary_observation_first = PETSC_FALSE
    mass_balance_first = PETSC_FALSE
    integral_flux_first = PETSC_FALSE
    observation_aggregate_first = PETSC_FALSE
  endif

  ewriter_summ_count = 0
  ewriter_rest_count = 0
  linerept_count     = 0

end subroutine OutputObservationInit

! ************************************************************************** !

subroutine OutputObservation(realization_base)
  ! 
  ! Main driver for all observation output subroutines
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module
  
  implicit none
  
  class(realization_base_type) :: realization_base

  if (realization_base%output_option%print_observation) then
    call OutputObservationTecplotColumnTXT(realization_base)
    call OutputAggregateToFile(realization_base)
    call OutputIntegralFlux(realization_base)
    if (realization_base%option%use_mc) then
      call OutputObservationTecplotSecTXT(realization_base)
    endif
  endif

end subroutine OutputObservation

! ************************************************************************** !

subroutine OutputObservationTecplotColumnTXT(realization_base)
  ! 
  ! Print to observation data to text file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Patch_module
  use Observation_module
  use Utility_module
  use String_module
 
  implicit none

  class(realization_base_type) :: realization_base
  
  PetscInt :: fid, icell
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  type(observation_type), pointer :: observation
  type(observation_aggregate_type), pointer :: aggregate
  type(output_variable_type), pointer :: cur_variable
  PetscBool, save :: open_file = PETSC_FALSE
  PetscReal, allocatable :: velocities(:,:,:)
  PetscInt :: local_id
  PetscInt :: icolumn
  PetscInt :: nphase
  PetscInt :: i
  PetscReal :: temp_real_comp, temp_real
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_output_observation,ierr);CHKERRQ(ierr)
  
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option
  
  if (check_for_obs_points) then
    open_file = PETSC_FALSE
    observation => patch%observation_list%first
    do
      if (.not.associated(observation)) exit
      if (observation%print_velocities) calculate_velocities = PETSC_TRUE
      if (observation%itype == OBSERVATION_SCALAR .or. &
          (observation%itype == OBSERVATION_FLUX .and. &
           option%myrank == option%io_rank)) then
        open_file = PETSC_TRUE
        exit
      endif
      observation => observation%next
    enddo
    check_for_obs_points = PETSC_FALSE
  endif

  if (calculate_velocities) then
    nphase = max(option%nphase,option%transport%nphase)
    allocate(velocities(3,realization_base%patch%grid%nlmax,nphase))
    call PatchGetCellCenteredVelocities(realization_base%patch, &
                                        ONE_INTEGER,velocities(:,:,1))
    if (nphase > 1) then
      call PatchGetCellCenteredVelocities(realization_base%patch, &
                                          TWO_INTEGER,velocities(:,:,2))
    endif
  endif
    
  if (open_file) then
    write(string,'(i6)') option%myrank
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-obs-' // trim(adjustl(string)) // '.pft'
  
    ! open file
    fid = 86
    if (observation_first .or. .not.FileExists(filename)) then
      open(unit=fid,file=filename,action="write",status="replace")
      ! write header
      ! write title
      write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // &
        ']"'
      observation => patch%observation_list%first

      ! must initialize icolumn here so that icolumn does not restart with
      ! each observation point
      if (output_option%print_column_ids) then
        icolumn = 1
      else
        icolumn = -1
      endif

      do 
        if (.not.associated(observation)) exit
        
        select case(observation%itype)
          case(OBSERVATION_SCALAR)
            if (associated(observation%region%coordinates) .and. &
                .not.observation%at_cell_center) then
 !             option%io_buffer = 'Writing of data at coordinates not ' // &
 !               'functioning properly for minerals.  Perhaps due to ' // &
 !               'non-ghosting of vol frac....>? - geh'
 !             call PrintErrMsg(option)
              call WriteObservationHeaderForCoord(fid,realization_base, &
                                                  observation%region, &
                                                 observation%print_velocities, &
                                                  icolumn)
            else
              do icell=1,observation%region%num_cells
                call WriteObservationHeaderForCell(fid,realization_base, &
                                                   observation%region,icell, &
                                                 observation%print_velocities, &
                                                   icolumn)
              enddo
            endif
          case(OBSERVATION_FLUX)
            if (option%myrank == option%io_rank) then
              call WriteObservationHeaderForBC(fid,realization_base, &
                                                observation%linkage_name)
            endif
        end select
        observation => observation%next
      enddo
      write(fid,'(a)',advance="yes") ""
    else
      open(unit=fid,file=filename,action="write",status="old", &
           position="append")
    endif
  
    observation => patch%observation_list%first
    write(fid,'(1es14.6)',advance="no") option%time/output_option%tconv
    do 
      if (.not.associated(observation)) exit
        select case(observation%itype)
          case(OBSERVATION_SCALAR)
            if (associated(observation%region%coordinates) .and. &
                .not.observation%at_cell_center) then
              call WriteObservationDataForCoord(fid,realization_base, &
                                                 observation%region)
              if (observation%print_velocities) then
                call WriteVelocityAtCoord(fid,realization_base, &
                                          observation%region)
              endif
            else
              do icell=1,observation%region%num_cells
                local_id = observation%region%cell_ids(icell)
                call WriteObservationDataForCell(fid,realization_base,local_id)
                if (observation%print_velocities) then
                  call WriteVelocityAtCell2(fid,realization_base,local_id, &
                                            velocities)
                endif
              enddo
            endif
          case(OBSERVATION_FLUX)
            call WriteObservationDataForBC(fid,realization_base, &
                                            patch, &
                                            observation%connection_set)
        end select
      observation => observation%next
    enddo
    write(fid,'(a)',advance="yes") ""
    close(fid)

  endif

  observation_first = PETSC_FALSE
  if (allocated(velocities)) deallocate(velocities)
  
  call PetscLogEventEnd(logging%event_output_observation,ierr);CHKERRQ(ierr)
      
end subroutine OutputObservationTecplotColumnTXT

! ************************************************************************** !

subroutine OutputAggregateToFile(realization_base)
  ! 
  ! Print to observation aggregate data to text file
  ! 
  ! Author: Michael Nole
  ! Date: 05/07/20
  !

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Patch_module
  use Observation_module
  use Utility_module
  use String_module

  implicit none

  class(realization_base_type) :: realization_base

  PetscInt :: fid, icell
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(observation_type), pointer :: observation
  type(observation_aggregate_type), pointer :: aggregate
  PetscInt :: icolumn
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_output_observation_agg,ierr);CHKERRQ(ierr)

  patch => realization_base%patch
  option => realization_base%option
  output_option => realization_base%output_option

  ! Write aggregate output separately because any rank could write.
  observation => patch%observation_list%first
  do
    if (.not. associated(observation)) exit
    if (observation%itype == OBSERVATION_AGGREGATE) then
      aggregate => observation%aggregate
      if (.not. associated(aggregate%output_variable)) then
        do
          if (.not. associated(aggregate)) exit
        ! Link aggregates with their output variables if it hasn't already
        ! been done.
          call ObservationAggregateLinkToVar(aggregate%output_variable, &
                                       output_option%output_obs_variable_list, &
                                       aggregate%var_name, option)
          aggregate => aggregate%next
        enddo
        aggregate => observation%aggregate
      endif

      do
        if (.not. associated(aggregate)) exit

        write(string,'(i6)') observation%id
        write(string2,'(i6)') aggregate%id
        filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-obs-' // trim(adjustl(string)) // '-agg-' // &
               trim(adjustl(string2)) // '.pft'
        fid = 86

        if (observation_aggregate_first .or. .not.FileExists(filename)) then
          if (option%myrank == option%io_rank) then
            open(unit=fid,file=filename,action="write",status="replace")
            ! write header
            ! write title
            write(fid,'(a)',advance="no") ' Aggregate Metric: '
            select case(aggregate%metric)
              case(OBSERVATION_AGGREGATE_MAX)
                 write(fid,'(a)',advance="no") 'Max '
              case(OBSERVATION_AGGREGATE_MIN)
                 write(fid,'(a)',advance="no") 'Min '
              case(OBSERVATION_AGGREGATE_AVG)
                 write(fid,'(a)',advance="no") 'Average '
            end select
            write(fid,'(a)',advance="yes") trim(aggregate%var_name)

            write(fid,'(a)',advance="no") ' "Time [' // &
                  trim(output_option%tunit) // ']"'

            ! must initialize icolumn here so that icolumn does not restart with
            ! each observation point
            if (output_option%print_column_ids) then
              icolumn = 1
            else
              icolumn = -1
            endif

            call WriteObservationHeaderAgg(fid,realization_base, &
                                           observation%region,icell, &
                                           observation%print_velocities, &
                                           icolumn)

            write(fid,'(a)',advance="yes") ""
            close(fid)
          endif
        endif
        aggregate => aggregate%next
      enddo

      aggregate => observation%aggregate
      do
        if (.not. associated(aggregate)) exit

        write(string,'(i6)') observation%id
        write(string2,'(i6)') aggregate%id
        filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-obs-' // trim(adjustl(string)) // '-agg-' // &
               trim(adjustl(string2)) // '.pft'
        fid = 86

        ! Compute the aggregate metric on each process
        call ObservationAggComputeMetric(realization_base, aggregate, &
                                         observation%region, option)
        ! Do the reduction and write
        call WriteObservationAggData(aggregate,realization_base, &
                                     string,filename,fid,output_option,option)
        aggregate => aggregate%next
      enddo

    endif
    observation => observation%next
  enddo

  observation_aggregate_first = PETSC_FALSE
  call PetscLogEventEnd(logging%event_output_observation_agg,ierr);CHKERRQ(ierr)

end subroutine OutputAggregateToFile

! ************************************************************************** !

subroutine WriteObservationHeaderAgg(fid,realization_base,region,icell, &
                                     print_velocities,icolumn)
  ! Print a header for aggregated data
  ! 
  ! Author: Michael Nole
  ! Date: 04/16/20
  ! 
  
  use Realization_Base_class, only : realization_base_type
  use Output_Aux_module
  use Region_module

  implicit none

  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(region_type) :: region
  PetscInt :: icell
  PetscBool :: print_velocities
  PetscInt :: icolumn

  character(len=MAXSTRINGLENGTH) :: cell_string

  cell_string = ' '

  call OutputWriteToHeader(fid,'x(m)',' ',cell_string,icolumn)
  call OutputWriteToHeader(fid,'y(m)',' ',cell_string,icolumn)
  call OutputWriteToHeader(fid,'z(m)',' ',cell_string,icolumn)

  call WriteObservationHeader(fid,realization_base,cell_string, &
                              print_velocities,icolumn)


end subroutine WriteObservationHeaderAgg

! ************************************************************************** !

subroutine WriteObservationHeaderForCell(fid,realization_base,region,icell, &
                                         print_velocities, &
                                         icolumn)
  ! 
  ! Print a header for data at a cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Output_Aux_module
  use Patch_module
  use Region_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(region_type) :: region
  PetscInt :: icell
  PetscBool :: print_velocities
  PetscInt :: icolumn
  
  PetscInt :: local_id
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  type(grid_type), pointer :: grid

  grid => realization_base%patch%grid
  
  local_id = region%cell_ids(icell)
  write(cell_string,*) grid%nG2A(grid%nL2G(region%cell_ids(icell)))
  cell_string = trim(region%name) // ' (' // trim(adjustl(cell_string)) // ')'

  ! add coordinate of cell center
  x_string = BestFloat(grid%x(grid%nL2G(local_id)),1.d4,1.d-2)
  y_string = BestFloat(grid%y(grid%nL2G(local_id)),1.d4,1.d-2)
  z_string = BestFloat(grid%z(grid%nL2G(local_id)),1.d4,1.d-2)
  cell_string = trim(cell_string) // ' (' // trim(adjustl(x_string)) // &
                ' ' // trim(adjustl(y_string)) // &
                ' ' // trim(adjustl(z_string)) // ')'
  
  call WriteObservationHeader(fid,realization_base,cell_string, &
                              print_velocities,icolumn)

end subroutine WriteObservationHeaderForCell

! ************************************************************************** !

subroutine WriteObservationHeaderForCoord(fid,realization_base,region, &
                                          print_velocities, &
                                          icolumn)
  ! 
  ! Print a header for data at a coordinate
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Patch_module
  use Region_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(region_type) :: region
  PetscBool :: print_velocities
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  
  cell_string = trim(region%name)
  
  x_string = BestFloat(region%coordinates(ONE_INTEGER)%x,1.d4,1.d-2)
  y_string = BestFloat(region%coordinates(ONE_INTEGER)%y,1.d4,1.d-2)
  z_string = BestFloat(region%coordinates(ONE_INTEGER)%z,1.d4,1.d-2)
  cell_string = trim(cell_string) // ' (' // trim(adjustl(x_string)) // ' ' // &
                trim(adjustl(y_string)) // ' ' // &
                trim(adjustl(z_string)) // ')'

  call WriteObservationHeader(fid,realization_base,cell_string, &
                              print_velocities,icolumn)

end subroutine WriteObservationHeaderForCoord

! ************************************************************************** !

subroutine WriteObservationHeader(fid,realization_base,cell_string, &
                                  print_velocities,icolumn)
  ! 
  ! Print a header for data
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/27/11
  ! 
                                  
  use Realization_Base_class, only : realization_base_type
  use Option_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  PetscBool :: print_velocities
  character(len=MAXSTRINGLENGTH) :: cell_string
  PetscInt :: icolumn
  
  PetscInt :: variable_count
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option  
  
  option => realization_base%option
  output_option => realization_base%output_option
  
  call OutputWriteVariableListToHeader(fid, &
                                       output_option%output_obs_variable_list, &
                                       cell_string,icolumn,PETSC_FALSE, &
                                       variable_count)

  if (print_velocities) then
    write(string,'(''m/'',a,'' '')') trim(realization_base%output_option%tunit)
    call OutputWriteToHeader(fid,'qlx',string,cell_string,icolumn)
    call OutputWriteToHeader(fid,'qly',string,cell_string,icolumn)
    call OutputWriteToHeader(fid,'qlz',string,cell_string,icolumn)

    if (max(option%nphase,option%transport%nphase) > 1) then
      call OutputWriteToHeader(fid,'qgx',string,cell_string,icolumn)
      call OutputWriteToHeader(fid,'qgy',string,cell_string,icolumn)
      call OutputWriteToHeader(fid,'qgz',string,cell_string,icolumn)
    endif
  endif
    
end subroutine WriteObservationHeader

! ************************************************************************** !

subroutine OutputObservationTecplotSecTXT(realization_base)
  ! 
  ! Print to secondary continuum observation
  ! data to text file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 04/08/13
  ! 

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Patch_module
  use Observation_module
  use Utility_module
 
  implicit none

  class(realization_base_type) :: realization_base
  
  PetscInt :: fid, icell
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  type(observation_type), pointer :: observation
  PetscBool, save :: open_file = PETSC_FALSE
  PetscInt :: local_id
  PetscInt :: icolumn
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_output_observation,ierr);CHKERRQ(ierr)
  
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option
  
  if (secondary_check_for_obs_points) then
    open_file = PETSC_FALSE
    observation => patch%observation_list%first
    do
      if (.not.associated(observation)) exit
      if (observation%itype == OBSERVATION_SCALAR .or. &
          (observation%itype == OBSERVATION_FLUX .and. &
           option%myrank == option%io_rank)) then
        open_file = PETSC_TRUE
        exit
      endif
      observation => observation%next
    enddo
    secondary_check_for_obs_points = PETSC_FALSE
  endif
  
  
  if (open_file) then
    write(string,'(i6)') option%myrank
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-obs-sec-' // trim(adjustl(string)) // '.pft'
  
    ! open file
    fid = 86
    if (secondary_observation_first .or. .not.FileExists(filename)) then
      open(unit=fid,file=filename,action="write",status="replace")
      ! write header
      ! write title
      write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // &
        ']"'
      observation => patch%observation_list%first

      ! must initialize icolumn here so that icolumn does not restart with
      ! each observation point
      if (output_option%print_column_ids) then
        icolumn = 1
      else
        icolumn = -1
      endif

      do 
        if (.not.associated(observation)) exit
        
        select case(observation%itype)
          case(OBSERVATION_SCALAR)
            if (associated(observation%region%coordinates) .and. &
                .not.observation%at_cell_center) then
              option%io_buffer = 'Writing of data at coordinates not &
                &functioning properly for minerals.  Perhaps due to &
                &non-ghosting of vol frac....>? - geh'
              call PrintErrMsg(option)
              call WriteObservationHeaderForCoordSec(fid,realization_base, &
                                                  observation%region, &
                                                  observation% &
                                                  print_secondary_data, &
                                                  icolumn)
            else
              do icell=1,observation%region%num_cells
                call WriteObservationHeaderForCellSec(fid,realization_base, &
                                                   observation%region,icell, &
                                                   observation% &
                                                   print_secondary_data, &
                                                   icolumn)
              enddo
            endif
        end select
        observation => observation%next
      enddo
      write(fid,'(a)',advance="yes") ""
    else
      open(unit=fid,file=filename,action="write",status="old", &
           position="append")
    endif
  
    observation => patch%observation_list%first
    write(fid,'(1es14.6)',advance="no") option%time/output_option%tconv
    do 
      if (.not.associated(observation)) exit
        select case(observation%itype)
          case(OBSERVATION_SCALAR)
              do icell=1,observation%region%num_cells
                local_id = observation%region%cell_ids(icell)
                if (observation%print_secondary_data(1)) then
                  call WriteObservationSecondaryDataAtCell(fid, &
                                                           realization_base, &
                                                           local_id, &
                                                           PRINT_SEC_TEMP)
                endif
                if (observation%print_secondary_data(2)) then
                  call WriteObservationSecondaryDataAtCell(fid, &
                                                           realization_base, &
                                                           local_id, &
                                                           PRINT_SEC_CONC)
                endif
                if (observation%print_secondary_data(3)) then
                  call WriteObservationSecondaryDataAtCell(fid, &
                                                           realization_base, &
                                                           local_id, &
                                                          PRINT_SEC_MIN_VOLFRAC)
                endif
                if (observation%print_secondary_data(4)) then
                  call WriteObservationSecondaryDataAtCell(fid, &
                                                           realization_base, &
                                                           local_id, &
                                                           PRINT_SEC_MIN_RATE)
                endif
                if (observation%print_secondary_data(5)) then
                  call WriteObservationSecondaryDataAtCell(fid, &
                                                           realization_base, &
                                                           local_id, &
                                                           PRINT_SEC_MIN_SI)
                endif
              enddo
      end select
      observation => observation%next
    enddo
    write(fid,'(a)',advance="yes") ""
    close(fid)

  endif

  secondary_observation_first = PETSC_FALSE
  
  call PetscLogEventEnd(logging%event_output_observation,ierr);CHKERRQ(ierr)
      
end subroutine OutputObservationTecplotSecTXT

! ************************************************************************** !

subroutine WriteObservationHeaderForCellSec(fid,realization_base,region,icell, &
                                            print_secondary_data, &
                                            icolumn)
  ! 
  ! Print a header for data at a cell for
  ! secondary continuum
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 04/08/13
  ! 

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Output_Aux_module
  use Patch_module
  use Region_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(region_type) :: region
  PetscInt :: icell
  PetscBool :: print_secondary_data(5)
  PetscInt :: icolumn
  
  PetscInt :: local_id
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  type(grid_type), pointer :: grid

  grid => realization_base%patch%grid
  
  local_id = region%cell_ids(icell)
  write(cell_string,*) grid%nG2A(grid%nL2G(region%cell_ids(icell)))
  cell_string = trim(region%name) // ' (' // trim(adjustl(cell_string)) // ')'

  ! add coordinate of cell center
  x_string = BestFloat(grid%x(grid%nL2G(local_id)),1.d4,1.d-2)
  y_string = BestFloat(grid%y(grid%nL2G(local_id)),1.d4,1.d-2)
  z_string = BestFloat(grid%z(grid%nL2G(local_id)),1.d4,1.d-2)
  cell_string = trim(cell_string) // ' (' // trim(adjustl(x_string)) // &
                ' ' // trim(adjustl(y_string)) // &
                ' ' // trim(adjustl(z_string)) // ')'
  
  call WriteObservationHeaderSec(fid,realization_base,cell_string, &
                                 print_secondary_data,icolumn)

end subroutine WriteObservationHeaderForCellSec

! ************************************************************************** !

subroutine WriteObservationHeaderForCoordSec(fid,realization_base,region, &
                                             print_secondary_data, &
                                             icolumn)
  ! 
  ! Print a header for data at a coordinate
  ! for secondary continuum
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 04/08/13
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Patch_module
  use Region_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(region_type) :: region
  PetscBool :: print_secondary_data(5)
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  
  cell_string = trim(region%name)
  
  x_string = BestFloat(region%coordinates(ONE_INTEGER)%x,1.d4,1.d-2)
  y_string = BestFloat(region%coordinates(ONE_INTEGER)%y,1.d4,1.d-2)
  z_string = BestFloat(region%coordinates(ONE_INTEGER)%z,1.d4,1.d-2)
  cell_string = trim(cell_string) // ' (' // trim(adjustl(x_string)) // ' ' // &
                trim(adjustl(y_string)) // ' ' // &
                trim(adjustl(z_string)) // ')'

  call WriteObservationHeaderSec(fid,realization_base,cell_string, &
                                 print_secondary_data,icolumn)

end subroutine WriteObservationHeaderForCoordSec

! ************************************************************************** !

subroutine WriteObservationHeaderSec(fid,realization_base,cell_string, &
                                     print_secondary_data,icolumn)
  ! 
  ! Print a header for secondary continuum data
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 10/27/13
  ! 
                                     
  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Reaction_Aux_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  class(reaction_rt_type), pointer :: reaction 
  PetscBool :: print_secondary_data(5)
  character(len=MAXSTRINGLENGTH) :: cell_string
  PetscInt :: icolumn
  
  PetscInt :: i,j
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option  
  
  option => realization_base%option
  output_option => realization_base%output_option
  
  ! add secondary temperature to header
  if (print_secondary_data(1)) then
    select case (option%iflowmode) 
      case (TH_MODE,TH_TS_MODE,MPH_MODE)
        do i = 1, option%nsec_cells
          write(string,'(i2)') i
          string = 'T(' // trim(adjustl(string)) // ')'
          call OutputWriteToHeader(fid,string,'C',cell_string,icolumn)
        enddo
      case default
    end select
  endif
  
  ! add secondary concentrations to header
  if (option%ntrandof > 0) then 
    select case(option%itranmode)
      case(RT_MODE)
        reaction => ReactionCast(realization_base%reaction_base)
        if (print_secondary_data(2)) then
          do j = 1, reaction%naqcomp
            do i = 1, option%nsec_cells
              write(string,'(i2)') i
              string = 'C(' // trim(adjustl(string)) // ') ' &
                         // trim(reaction%primary_species_names(j))
              call OutputWriteToHeader(fid,string,'molal',cell_string, &
                                       icolumn)
            enddo
          enddo
        endif
      
      ! add secondary mineral volume fractions to header
        if (print_secondary_data(3)) then
          do j = 1, reaction%mineral%nkinmnrl
            do i = 1, option%nsec_cells
              write(string,'(i2)') i
              string = 'VF(' // trim(adjustl(string)) // ') ' &
                       // trim(reaction%mineral%mineral_names(j))
              call OutputWriteToHeader(fid,string,'',cell_string,icolumn)
            enddo
          enddo
        endif  
        
      ! add secondary mineral rates to header
        if (print_secondary_data(4)) then
          do j = 1, reaction%mineral%nkinmnrl
            do i = 1, option%nsec_cells
              write(string,'(i2)') i
              string = 'Rate(' // trim(adjustl(string)) // ') ' &
                       // trim(reaction%mineral%mineral_names(j))
              call OutputWriteToHeader(fid,string,'',cell_string,icolumn)
            enddo
          enddo
        endif    
        
      ! add secondary mineral volume fractions to header
        if (print_secondary_data(5)) then
          do j = 1, reaction%mineral%nkinmnrl
            do i = 1, option%nsec_cells
              write(string,'(i2)') i
              string = 'SI(' // trim(adjustl(string)) // ') ' &
                       // trim(reaction%mineral%mineral_names(j))
              call OutputWriteToHeader(fid,string,'',cell_string,icolumn)
            enddo
          enddo
        endif    
      case(NWT_MODE)
        option%io_buffer = 'WriteObservationHeaderSec has not been setup &
          &for NW Transport.'
        call PrintErrMsg(option)
    end select
  endif 
  
end subroutine WriteObservationHeaderSec

! ************************************************************************** !

subroutine WriteObservationHeaderForBC(fid,realization_base,coupler_name)
  ! 
  ! Print a header for data over a region
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/18/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Reaction_Aux_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  character(len=MAXWORDLENGTH) :: coupler_name
  
  PetscInt :: i
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  class(reaction_rt_type), pointer :: reaction 
  
  option => realization_base%option
  reaction => ReactionCast(realization_base%reaction_base)
  
  select case(option%iflowmode)
    case(FLASH2_MODE)
    case(MPH_MODE)
    case(IMS_MODE)
    case(TH_MODE,TH_TS_MODE)
    case(MIS_MODE)
    case(RICHARDS_MODE,RICHARDS_TS_MODE)
      string = ',"Darcy flux ' // trim(coupler_name) // &
               ' [m^3/' // trim(realization_base%output_option%tunit) // ']"'
    case default
  end select
  write(fid,'(a)',advance="no") trim(string)

  if (associated(reaction)) then
    do i=1, reaction%naqcomp 
      ! may need to modify for molality vs molarity, but I believe molarity is correct
      write(fid,'(a)',advance="no") ',"' // &
        trim(reaction%primary_species_names(i)) // ' ' // &
        trim(coupler_name) // &
        ' [mol/' // trim(realization_base%output_option%tunit) // ']"'
    enddo
  endif

end subroutine WriteObservationHeaderForBC

! ************************************************************************** !

subroutine WriteObservationAggData(aggregate,realization_base,string,&
                                   filename,fid,output_option,option)
  !
  ! Print aggregate data of interest
  !
  ! Author: Michael Nole
  ! Date: 04/16/20

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Observation_module

  implicit none

  type(observation_aggregate_type), pointer :: aggregate
  class(realization_base_type) :: realization_base
  type(output_option_type), pointer :: output_option
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string, filename
  PetscInt :: fid  

  PetscErrorCode :: ierr
  PetscInt :: agg_rank, local_id
  PetscReal :: local_metric(2), global_metric(2)

110 format(es14.6)

  local_metric(1) = aggregate%metric_value
  local_metric(2) = option%myrank

  select case(aggregate%metric)
    case(OBSERVATION_AGGREGATE_MAX)
      call MPI_Allreduce(local_metric,global_metric,ONE_INTEGER, &
                         MPI_2DOUBLE_PRECISION,MPI_MAXLOC,option%mycomm,ierr) 
  end select

  agg_rank = global_metric(2)

  if (option%myrank == agg_rank) then

    local_id = aggregate%local_id

    open(unit=fid,file=filename,action="write",status="old",position="append")
    write(fid,'(1es14.6)',advance="no") option%time/output_option%tconv

    write(fid,110,advance="no") realization_base%patch%grid%x(realization_base%&
                                patch%grid%nL2G(local_id))
    write(fid,110,advance="no") realization_base%patch%grid%y(realization_base%&
                                patch%grid%nL2G(local_id))
    write(fid,110,advance="no") realization_base%patch%grid%z(realization_base%&
                                patch%grid%nL2G(local_id))

    call WriteObservationDataForCell(fid,realization_base,local_id)

    write(fid,'(a)',advance="yes") ""
    close(fid)    
  endif

end subroutine WriteObservationAggData

! ************************************************************************** !

subroutine WriteObservationDataForCell(fid,realization_base,local_id)
  ! 
  ! Print data for data at a cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  use Realization_Base_class, only : realization_base_type, &
                                     RealizGetVariableValueAtCell
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Variables_module
  
  implicit none
  
  PetscInt :: fid, i
  class(realization_base_type) :: realization_base
  PetscInt :: local_id
  PetscReal :: temp_real
  PetscInt :: ghosted_id
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option  
  type(output_variable_type), pointer :: cur_variable
  
  option => realization_base%option
  patch => realization_base%patch
  grid => patch%grid
  field => realization_base%field
  output_option => realization_base%output_option

110 format(es14.6)
111 format(i2)

  ghosted_id = grid%nL2G(local_id)

  ! loop over observation variables and write to file
  cur_variable => output_option%output_obs_variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    if (cur_variable%plot_only) then
      cur_variable => cur_variable%next
      cycle
    endif     
    temp_real = OutputGetVariableAtCell(realization_base,ghosted_id, &
                                        cur_variable)
    if (cur_variable%iformat == 0) then ! real
      write(fid,110,advance="no") temp_real
    else ! integer
      write(fid,111,advance="no") int(temp_real + 1.d-5)
    endif
    cur_variable => cur_variable%next
  enddo  

end subroutine WriteObservationDataForCell

! ************************************************************************** !

subroutine WriteObservationDataForCoord(fid,realization_base,region)
  ! 
  ! Print data for data at a coordinate
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Region_module  
  use Grid_module
  use Field_module
  use Patch_module
  use Variables_module
  
  use Grid_Structured_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(region_type) :: region

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: temp_real
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
    
  PetscInt :: ghosted_ids(8)
  PetscInt :: count
  PetscInt :: i, j, k
  PetscInt :: istart, iend, jstart, jend, kstart, kend
  
  option => realization_base%option
  patch => realization_base%patch
  grid => patch%grid
  field => realization_base%field
  output_option => realization_base%output_option

110 format(es14.6)
111 format(i2)

  count = 0
  local_id = region%cell_ids(1)
  ghosted_id = grid%nL2G(local_id)
  call StructGridGetIJKFromGhostedID(grid%structured_grid,ghosted_id,i,j,k)
  istart = i
  iend = i
  jstart = j
  jend = j
  kstart = k
  kend = k
  ! find the neighboring cells, between which to interpolate
  if (grid%x(ghosted_id) > region%coordinates(ONE_INTEGER)%x) then
    if (i > 1) then
      istart = i-1
    endif
  else
    if (i < grid%structured_grid%ngx) then
      iend = i+1
    endif
  endif
  if (grid%y(ghosted_id) > region%coordinates(ONE_INTEGER)%y) then
    if (j > 1) then
      jstart = j-1
    endif
  else
    if (j < grid%structured_grid%ngy) then
      jend = j+1
    endif
  endif
  if (grid%z(ghosted_id) > region%coordinates(ONE_INTEGER)%z) then
    if (k > 1) then
      kstart = k-1
    endif
  else
    if (k < grid%structured_grid%ngz) then
      kend = k+1
    endif
  endif
  count = 0
  do k=kstart,kend
    do j=jstart,jend
      do i=istart,iend
        count = count + 1
        ghosted_ids(count) = i + (j-1)*grid%structured_grid%ngx + &
                             (k-1)*grid%structured_grid%ngxy
      enddo
    enddo
  enddo
  
  ! loop over observation variables and write to file
  cur_variable => output_option%output_obs_variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    if (cur_variable%plot_only) then
      cur_variable => cur_variable%next
      cycle
    endif    
    temp_real = OutputGetVariableAtCoord(realization_base, &
                                         cur_variable, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
    if (cur_variable%iformat == 0) then ! real
      write(fid,110,advance="no") temp_real
    else ! integer
      write(fid,111,advance="no") int(temp_real + 1.d-5)
    endif
    cur_variable => cur_variable%next
  enddo

end subroutine WriteObservationDataForCoord

! ************************************************************************** !

subroutine WriteObservationDataForBC(fid,realization_base,patch,connection_set)
  ! 
  ! Print flux data for a boundary condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/18/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Connection_module  
  use Patch_module
  use Reaction_Aux_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(patch_type), pointer :: patch
  type(connection_set_type), pointer :: connection_set

  PetscInt :: i
  PetscInt :: iconn
  PetscInt :: offset
  PetscInt :: iphase
  PetscMPIInt :: int_mpi
  PetscReal :: sum_volumetric_flux(realization_base%option%nphase)
  PetscReal :: sum_volumetric_flux_global(realization_base%option%nphase)
  PetscReal :: sum_solute_flux(realization_base%option%ntrandof)
  PetscReal :: sum_solute_flux_global(realization_base%option%ntrandof)
  type(option_type), pointer :: option
  class(reaction_rt_type), pointer :: reaction
  PetscErrorCode :: ierr
  
  option => realization_base%option
  reaction => ReactionCast(realization_base%reaction_base)

110 format(es14.6)
 
  iphase = 1

  ! sum up fluxes across region
  if (associated(connection_set)) then
    offset = connection_set%offset
    select case(option%iflowmode)
      case(MPH_MODE,TH_MODE,TH_TS_MODE,IMS_MODE,FLASH2_MODE,G_MODE,H_MODE)
      case(WF_MODE)
        option%io_buffer = 'WriteObservationDataForBC() needs to be set up &
          & for WIPP Flow, and perhaps the other multiphase flow modes.'
        call PrintErrMsg(option)
      case(MIS_MODE)
      case(RICHARDS_MODE,RICHARDS_TS_MODE)
        sum_volumetric_flux = 0.d0
        if (associated(connection_set)) then
          do iconn = 1, connection_set%num_connections
            sum_volumetric_flux(:) = sum_volumetric_flux(:) + &
                            patch%boundary_velocities(iphase,offset+iconn)* &
                            connection_set%area(iconn)
          enddo
        endif
        int_mpi = option%nphase
        call MPI_Reduce(sum_volumetric_flux,sum_volumetric_flux_global, &
                        int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        option%io_rank,option%mycomm,ierr)
        if (option%myrank == option%io_rank) then
          do i = 1, option%nphase
            write(fid,110,advance="no") sum_volumetric_flux_global(i)
          enddo
        endif
    end select

    if (associated(reaction)) then
      sum_solute_flux = 0.d0
      if (associated(connection_set)) then
        do iconn = 1, connection_set%num_connections
          sum_solute_flux(:) = sum_solute_flux(:) + &
                               patch%boundary_tran_fluxes(:,offset+iconn)* &
                               connection_set%area(iconn)
        enddo
      endif
      int_mpi = option%ntrandof
      call MPI_Reduce(sum_solute_flux,sum_solute_flux_global, &
                      int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)
      if (option%myrank == option%io_rank) then
        !we currently only print the aqueous components
        do i = 1, reaction%naqcomp
          write(fid,110,advance="no") sum_solute_flux_global(i)
        enddo
      endif
    endif

  endif

end subroutine WriteObservationDataForBC

! ************************************************************************** !

subroutine WriteVelocityAtCell(fid,realization_base,local_id)
  ! 
  ! Computes velocities at a grid cell
  ! note: limited to structured grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  PetscInt :: local_id
  PetscInt :: iphase

  PetscReal :: velocity(1:3)
  option => realization_base%option
  
200 format(3(es14.6))

  iphase = 1
  velocity = GetVelocityAtCell(fid,realization_base,local_id,iphase)
  write(fid,200,advance="no") velocity(1:3)* &
                              realization_base%output_option%tconv

  if (max(option%nphase,option%transport%nphase) > 1) then
    iphase = 2
    velocity = GetVelocityAtCell(fid,realization_base,local_id,iphase)
    write(fid,200,advance="no") velocity(1:3)* &
                                realization_base%output_option%tconv
  endif

end subroutine WriteVelocityAtCell

! ************************************************************************** !

subroutine WriteVelocityAtCell2(fid,realization_base,local_id,velocities)
  ! 
  ! Writes the velocity previoiusly calculated and stored in vecs at cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/08/16
  ! 
  use Realization_Base_class, only : realization_base_type
  use Option_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  PetscInt :: local_id
  PetscReal :: velocities(:,:,:)

  PetscReal :: velocity(3)
  option => realization_base%option
  
200 format(3(es14.6))

  velocity = velocities(:,local_id,1)
  write(fid,200,advance="no") velocity(1:3)* &
                              realization_base%output_option%tconv

  if (max(option%nphase,option%transport%nphase) > 1) then
    velocity = velocities(:,local_id,2)
    write(fid,200,advance="no") velocity(1:3)* &
                                realization_base%output_option%tconv
  endif

end subroutine WriteVelocityAtCell2

! ************************************************************************** !

function GetVelocityAtCell(fid,realization_base,local_id,iphase)
  ! 
  ! Computes velocities at a grid cell
  ! note: limited to structured grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Connection_module
  use Coupler_module

  implicit none
  
  PetscReal :: GetVelocityAtCell(3)
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  PetscInt :: local_id

  PetscInt :: ghosted_id
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(connection_set_list_type), pointer :: connection_set_list
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: direction, iphase
  PetscReal :: area
  PetscReal :: sum_velocity(1:3), sum_area(1:3), velocity(1:3)
  
  option => realization_base%option
  patch => realization_base%patch
  grid => patch%grid
  field => realization_base%field

  sum_velocity = 0.d0
  sum_area = 0.d0
! iphase = 1

  ! interior velocities  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id_up = grid%nG2L(cur_connection_set%id_up(iconn)) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(cur_connection_set%id_dn(iconn)) ! = zero for ghost nodes
      if (local_id_up == local_id .or. local_id_dn == local_id) then
        do direction=1,3        
          area = cur_connection_set%area(iconn)* &
                 !geh: no dabs() here
                 cur_connection_set%dist(direction,iconn)
          sum_velocity(direction) = sum_velocity(direction) + &
                                    patch%internal_velocities(iphase,sum_connection)* &
                                    area
          sum_area(direction) = sum_area(direction) + dabs(area)
        enddo
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! boundary velocities
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      if (cur_connection_set%id_dn(iconn) == local_id) then
        do direction=1,3        
          area = cur_connection_set%area(iconn)* &
                 !geh: no dabs() here
                 cur_connection_set%dist(direction,iconn)
          sum_velocity(direction) = sum_velocity(direction) + &
                                    patch%boundary_velocities(iphase,sum_connection)* &
                                    area
          sum_area(direction) = sum_area(direction) + dabs(area)
        enddo
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  velocity = 0.d0
  do direction = 1,3
    if (abs(sum_area(direction)) > 1.d-40) &
      velocity(direction) = sum_velocity(direction)/sum_area(direction)
  enddo

  GetVelocityAtCell = velocity  

end function GetVelocityAtCell

! ************************************************************************** !

subroutine WriteVelocityAtCoord(fid,realization_base,region)
  ! 
  ! Computes velocities at a coordinate
  ! note: limited to structured grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Region_module
  use Option_module

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  type(region_type) :: region
  type(option_type), pointer :: option
  PetscInt :: local_id
  PetscInt :: iphase
  PetscReal :: coordinate(3)

  PetscReal :: velocity(1:3)

  option => realization_base%option
  
200 format(3(es14.6))

  iphase = 1
  velocity = GetVelocityAtCoord(fid,realization_base,region%cell_ids(1), &
                                region%coordinates(ONE_INTEGER)%x, &
                                region%coordinates(ONE_INTEGER)%y, &
                                region%coordinates(ONE_INTEGER)%z,iphase)
  write(fid,200,advance="no") velocity(1:3)*realization_base%output_option%tconv   

  if (max(option%nphase,option%transport%nphase) > 1) then
    iphase = 2
    velocity = GetVelocityAtCoord(fid,realization_base,region%cell_ids(1), &
                                region%coordinates(ONE_INTEGER)%x, &
                                region%coordinates(ONE_INTEGER)%y, &
                                region%coordinates(ONE_INTEGER)%z,iphase)
    write(fid,200,advance="no") velocity(1:3)*realization_base%output_option%tconv   
  endif

end subroutine WriteVelocityAtCoord

! ************************************************************************** !

function GetVelocityAtCoord(fid,realization_base,local_id,x,y,z,iphase)
  ! 
  ! Computes velocities at a coordinate
  ! note: limited to structured grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/08
  ! 
  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Connection_module
  use Coupler_module

  implicit none
  
  PetscReal :: GetVelocityAtCoord(3)
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  PetscInt :: local_id
  PetscReal :: x, y, z
  
  PetscInt :: ghosted_id
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(connection_set_list_type), pointer :: connection_set_list
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscInt :: local_id_up, local_id_dn
  PetscReal :: cell_coord(3), face_coord
  PetscReal :: coordinate(3)
  PetscInt :: direction, iphase
  PetscReal :: area, weight, distance
  PetscReal :: sum_velocity(1:3), velocity(1:3)
  PetscReal :: sum_weight(1:3)
  
  option => realization_base%option
  patch => realization_base%patch
  grid => patch%grid
  field => realization_base%field

  sum_velocity = 0.d0
  sum_weight = 0.d0
! iphase = 1

  ghosted_id = grid%nL2G(local_id)
  
  coordinate(X_DIRECTION) = x
  coordinate(Y_DIRECTION) = y
  coordinate(Z_DIRECTION) = z

  cell_coord(X_DIRECTION) = grid%x(ghosted_id)
  cell_coord(Y_DIRECTION) = grid%y(ghosted_id)
  cell_coord(Z_DIRECTION) = grid%z(ghosted_id)

  ! interior velocities  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id_up = grid%nG2L(cur_connection_set%id_up(iconn)) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(cur_connection_set%id_dn(iconn)) ! = zero for ghost nodes
      if (local_id_up == local_id .or. local_id_dn == local_id) then
        do direction=1,3
          if (local_id_up == local_id) then
            face_coord = cell_coord(direction) + &
                         cur_connection_set%dist(-1,iconn)* &
                         cur_connection_set%dist(0,iconn)* &
                         cur_connection_set%dist(direction,iconn)
          else
            face_coord = cell_coord(direction) - &
                         (1.d0-cur_connection_set%dist(-1,iconn))* &
                         cur_connection_set%dist(0,iconn)* &
                         cur_connection_set%dist(direction,iconn)
          endif
          distance = dabs(face_coord-coordinate(direction))
          if (distance < 1.d-40) distance = 1.d-40
          weight = cur_connection_set%area(iconn)* &
                 dabs(cur_connection_set%dist(direction,iconn))/ &
                 distance
 
          sum_velocity(direction) = sum_velocity(direction) + &
                                    cur_connection_set%dist(direction,iconn)* &
                                    patch%internal_velocities(iphase,sum_connection)* &
                                    weight
          sum_weight(direction) = sum_weight(direction) + weight
       enddo
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! boundary velocities
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      if (cur_connection_set%id_dn(iconn) == local_id) then
        do direction=1,3        
          face_coord = cell_coord(direction) - &
                    !   (1.d0-cur_connection_set%dist(-1,iconn))* & ! fraction upwind is always 0.d0
                       cur_connection_set%dist(0,iconn)* &
                       cur_connection_set%dist(direction,iconn)
          distance = dabs(face_coord-coordinate(direction))
          if (distance < 1.d-40) distance = 1.d-40
          weight = cur_connection_set%area(iconn)* &
                   dabs(cur_connection_set%dist(direction,iconn))/ &
                   distance
          sum_velocity(direction) = sum_velocity(direction) + &
                                    cur_connection_set%dist(direction,iconn)* &
                                    patch%boundary_velocities(iphase,sum_connection)* &
                                    weight
          sum_weight(direction) = sum_weight(direction) + weight
        enddo
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  velocity = 0.d0
  do direction = 1,3
    if (abs(sum_weight(direction)) > 1.d-40) &
      velocity(direction) = sum_velocity(direction)/sum_weight(direction)
  enddo

  GetVelocityAtCoord = velocity  

end function GetVelocityAtCoord

! ************************************************************************** !

subroutine WriteObservationSecondaryDataAtCell(fid,realization_base,local_id,ivar)
  ! 
  ! Print data for data at a cell
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 10/4/12
  ! 

  use Realization_Base_class, only : realization_base_type, &
                                     RealizGetVariableValueAtCell
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Variables_module
  use Reaction_Aux_module

  implicit none
  
  PetscInt :: fid,i,naqcomp,nkinmnrl
  class(realization_base_type) :: realization_base
  PetscInt :: local_id
  PetscInt :: ghosted_id
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option 
  class(reaction_rt_type), pointer :: reaction   
  PetscInt :: ivar
  
  option => realization_base%option
  patch => realization_base%patch
  grid => patch%grid
  field => realization_base%field
  output_option => realization_base%output_option

110 format(es14.6)

  ghosted_id = grid%nL2G(local_id)

  if (option%nsec_cells > 0) then
    if (ivar == PRINT_SEC_TEMP) then
      select case(option%iflowmode)
        case(MPH_MODE,TH_MODE,TH_TS_MODE)
          do i = 1, option%nsec_cells 
            write(fid,110,advance="no") &
              RealizGetVariableValueAtCell(realization_base,ghosted_id, &
                                           SECONDARY_TEMPERATURE,i)
          enddo
        end select
     endif
    if (option%ntrandof > 0) then
      select case(option%itranmode)
        case(RT_MODE)
          reaction => ReactionCast(realization_base%reaction_base)
          if (ivar == PRINT_SEC_CONC) then
            do naqcomp = 1, reaction%naqcomp
              do i = 1, option%nsec_cells 
                write(fid,110,advance="no") &
                  RealizGetVariableValueAtCell(realization_base,ghosted_id, &
                                             SECONDARY_CONCENTRATION,i,naqcomp)
              enddo
            enddo 
          endif
          if (ivar == PRINT_SEC_MIN_VOLFRAC) then
            do nkinmnrl = 1, reaction%mineral%nkinmnrl
              do i = 1, option%nsec_cells 
                write(fid,110,advance="no") &
                  RealizGetVariableValueAtCell(realization_base,ghosted_id, &
                                               SEC_MIN_VOLFRAC,i,nkinmnrl)
              enddo
            enddo
          endif
           if (ivar == PRINT_SEC_MIN_RATE) then
            do nkinmnrl = 1, reaction%mineral%nkinmnrl
              do i = 1, option%nsec_cells 
                write(fid,110,advance="no") &
                  RealizGetVariableValueAtCell(realization_base,ghosted_id, &
                                               SEC_MIN_RATE,i,nkinmnrl)
              enddo
            enddo
          endif
          if (ivar == PRINT_SEC_MIN_SI) then
            do nkinmnrl = 1, reaction%mineral%nkinmnrl
              do i = 1, option%nsec_cells 
                write(fid,110,advance="no") &
                  RealizGetVariableValueAtCell(realization_base,ghosted_id, &
                                               SEC_MIN_SI,i,nkinmnrl)
              enddo
            enddo
          endif           
        case(NWT_MODE)
          option%io_buffer = 'WriteObservationSecondaryDataAtCell has not &
            &been setup for NW Transport.'
          call PrintErrMsg(option)
      end select
    endif 
  endif 
   
end subroutine WriteObservationSecondaryDataAtCell

! ************************************************************************** !

subroutine ObservationAggregateLinkToVar(aggregate_var,output_var_list, &
                                         var_name,option)
  !
  ! Links aggregator variable to output variable
  !
  ! Author: Michael Nole
  ! Date: 04/17/20
 
  use Option_module 
  use String_module

  implicit none

  type(output_variable_type), pointer :: aggregate_var
  type(output_variable_list_type), pointer :: output_var_list
  character(len=MAXWORDLENGTH) :: var_name
  type(option_type), pointer :: option

  type(output_variable_type), pointer :: cur_variable

  cur_variable => output_var_list%first
  do
    if (.not. associated(cur_variable)) then
      option%io_buffer = 'Variable requested for aggregate metric ' //&
                         'does not match any output variables.'
      call PrintErrMsg(option)
    elseif (StringCompareIgnoreCase(cur_variable%name,var_name)) then
      aggregate_var => cur_variable
      exit
    endif
      cur_variable => cur_variable%next
  enddo

end subroutine ObservationAggregateLinkToVar

! ************************************************************************** !

subroutine ObservationAggComputeMetric(realization_base,aggregate,region,option)
  !
  ! Computes the user-specified aggregate metric
  !
  ! Author: Michael Nole
  ! Date: 04/17/20

  use Realization_Base_class, only : realization_base_type
  use Observation_module
  use Option_module
  use Region_module

  implicit none

  class(realization_base_type) :: realization_base
  type(observation_aggregate_type), pointer :: aggregate
  type(region_type), pointer :: region
  type(option_type), pointer :: option

  type(output_variable_type), pointer :: cur_variable
  PetscReal :: temp_real, temp_real_comp
  PetscInt :: icell, local_id, ghosted_id

  cur_variable => aggregate%output_variable
  temp_real = UNINITIALIZED_DOUBLE
  
  if (region%num_cells == 0) then
    select case(aggregate%metric)
      case(OBSERVATION_AGGREGATE_MAX)
        aggregate%metric_value = -1.d20
        aggregate%local_id = 0
      case(OBSERVATION_AGGREGATE_MIN)
        aggregate%metric_value = 1.d20
        aggregate%local_id = 0
    end select
  else
    do icell= 1,region%num_cells
      local_id = region%cell_ids(icell)
      ghosted_id = realization_base%patch%grid%nL2G(local_id)
      temp_real_comp = OutputGetVariableAtCell(realization_base, &
                                               ghosted_id,cur_variable)
      select case(aggregate%metric)
        case(OBSERVATION_AGGREGATE_MAX)
          if (temp_real_comp > temp_real .or. temp_real == &
              UNINITIALIZED_DOUBLE) then
            aggregate%metric_value = temp_real_comp
            aggregate%local_id = local_id
            temp_real = temp_real_comp
          endif
        case(OBSERVATION_AGGREGATE_MIN)
          if (temp_real_comp < temp_real .or. temp_real == &
              UNINITIALIZED_DOUBLE) then
            aggregate%metric_value = temp_real_comp
            aggregate%local_id = local_id
            temp_real = temp_real_comp
          endif
        case default
          option%io_buffer = 'Aggregate metric not assigned.'
          call PrintErrMsg(option)
      end select
    enddo
  endif

end subroutine ObservationAggComputeMetric

! ************************************************************************** !

subroutine OutputIntegralFlux(realization_base)
  ! 
  ! Print integral fluxes to Tecplot POINT format
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/14
  ! 

  use Realization_Subsurface_class, only : realization_subsurface_type
  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Grid_module
  use Patch_module
  use Output_Aux_module
  use Reaction_Aux_module
  use Integral_Flux_module
  use Utility_module
  use General_Aux_module, only : general_fmw => fmw_comp
  use WIPP_Flow_Aux_module, only : wipp_flow_fmw => fmw_comp

  implicit none

  class(realization_base_type), target :: realization_base

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(output_option_type), pointer :: output_option
  class(reaction_rt_type), pointer :: reaction

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: word, units
  character(len=MAXSTRINGLENGTH) :: string
  type(integral_flux_type), pointer :: integral_flux
  PetscReal :: flow_dof_scale(10)
  PetscReal, allocatable :: array(:,:)
  PetscReal, allocatable :: array_global(:,:)
  PetscReal, allocatable :: instantaneous_array(:)
  PetscInt, parameter :: fid = 86
  PetscInt :: i, j
  PetscInt :: istart, iend
  PetscInt :: icol
  PetscMPIInt :: int_mpi
  PetscReal :: tempreal
  PetscErrorCode :: ierr

  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  output_option => realization_base%output_option
  reaction => ReactionCast(realization_base%reaction_base)

  if (.not.associated(patch%integral_flux_list%first)) return

  flow_dof_scale = 1.d0
  select case(option%iflowmode)
    case(RICHARDS_MODE,RICHARDS_TS_MODE)
      flow_dof_scale(1) = FMWH2O
    case(TH_MODE,TH_TS_MODE)
      flow_dof_scale(1) = FMWH2O
    case(MIS_MODE)
      flow_dof_scale(1) = FMWH2O
      flow_dof_scale(2) = FMWGLYC
    case(G_MODE,H_MODE)
      flow_dof_scale(1) = FMWH2O
      flow_dof_scale(2) = general_fmw(2)
    case(WF_MODE)
      flow_dof_scale(1) = FMWH2O
      flow_dof_scale(2) = wipp_flow_fmw(2)
    case(MPH_MODE,FLASH2_MODE,IMS_MODE)
      flow_dof_scale(1) = FMWH2O
      flow_dof_scale(2) = FMWCO2
  end select

  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '-int.dat'
  else
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-int.dat'
  endif
  
  ! open file
  if (option%myrank == option%io_rank) then

    if (output_option%print_column_ids) then
      icol = 1
    else
      icol = -1
    endif
  
    if (integral_flux_first .or. .not.FileExists(filename)) then
      open(unit=fid,file=filename,action="write",status="replace")

      ! write header
      write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // &
        ']"'  
      
      if (option%iflowmode > 0) then
        call OutputWriteToHeader(fid,'dt_flow',output_option%tunit,'',icol)
      endif
      
      if (option%ntrandof > 0) then
        call OutputWriteToHeader(fid,'dt_tran',output_option%tunit,'',icol)
      endif
      
      integral_flux => patch%integral_flux_list%first
      do
        if (.not.associated(integral_flux)) exit
        select case(option%iflowmode)
          case(RICHARDS_MODE,RICHARDS_TS_MODE, &
               TH_MODE,TH_TS_MODE,MIS_MODE,G_MODE,H_MODE,MPH_MODE,FLASH2_MODE, &
               IMS_MODE,WF_MODE)
            string = trim(integral_flux%name) // ' Water'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(integral_flux%name) // ' Water'
            call OutputWriteToHeader(fid,string,units,'',icol)
        end select
        select case(option%iflowmode)
          case(MIS_MODE)
            string = trim(integral_flux%name) // ' Glycol'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(integral_flux%name) // ' Glycol'
            call OutputWriteToHeader(fid,string,units,'',icol)
          case(G_MODE,H_MODE)
            string = trim(integral_flux%name) // ' Air'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(integral_flux%name) // ' Air'
            call OutputWriteToHeader(fid,string,units,'',icol)
          case(WF_MODE)
            string = trim(integral_flux%name) // ' Gas Component'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(integral_flux%name) // ' Gas Component'
            call OutputWriteToHeader(fid,string,units,'',icol)
          case(MPH_MODE,FLASH2_MODE,IMS_MODE)
            string = trim(integral_flux%name) // ' CO2'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(integral_flux%name) // ' CO2'
            call OutputWriteToHeader(fid,string,units,'',icol)
        end select
        select case(option%iflowmode)
          case(TH_MODE,TH_TS_MODE,MIS_MODE,G_MODE,H_MODE,MPH_MODE,FLASH2_MODE, &
               IMS_MODE)
            string = trim(integral_flux%name) // ' Energy'
            call OutputWriteToHeader(fid,string,'MJ','',icol)
            units = 'MJ/' // trim(output_option%tunit) // ''
            string = trim(integral_flux%name) // ' Energy'
            call OutputWriteToHeader(fid,string,units,'',icol)
        end select
        
        if (option%ntrandof > 0) then
          select case(option%itranmode)
            case(RT_MODE)
              units = 'mol/' // trim(output_option%tunit) // ''
              do i=1,reaction%naqcomp
                if (reaction%primary_species_print(i)) then
                  string = trim(integral_flux%name) // ' ' // &
                           trim(reaction%primary_species_names(i))
                  call OutputWriteToHeader(fid,string,'mol','',icol)
                  string = trim(integral_flux%name) // ' ' // &
                           trim(reaction%primary_species_names(i))
                  call OutputWriteToHeader(fid,string,units,'',icol)
                endif
              enddo
            case(NWT_MODE)
              option%io_buffer = 'OutputIntegralFlux has not &
                &been setup for NW Transport.'
              call PrintErrMsg(option)
          end select
        endif
        integral_flux => integral_flux%next
      enddo
      write(fid,'(a)') '' 
    else
      open(unit=fid,file=filename,action="write",status="old",position="append")
    endif 
  endif     

100 format(100es17.8)
110 format(100es17.8)
120 format(100es17.8e3)

  ! write time
  if (option%myrank == option%io_rank) then
    write(fid,100,advance="no") option%time/output_option%tconv
  endif
  
  if (option%nflowdof > 0) then
    if (option%myrank == option%io_rank) &
      write(fid,100,advance="no") option%flow_dt/output_option%tconv
  endif
  if (option%ntrandof > 0) then
    if (option%myrank == option%io_rank) &
      write(fid,100,advance="no") option%tran_dt/output_option%tconv
  endif
  
  allocate(array(option%nflowdof + option%ntrandof,2))
  allocate(array_global(option%nflowdof + option%ntrandof,2))
  allocate(instantaneous_array(max(option%nflowdof,option%ntrandof)))
  integral_flux => patch%integral_flux_list%first
  do
    if (.not.associated(integral_flux)) exit
    array = 0.d0
    array_global = 0.d0
    if (option%nflowdof > 0) then
      istart = 1
      iend = option%nflowdof
      instantaneous_array = 0.d0
      call IntegralFluxGetInstantaneous(integral_flux, &
                                        patch%internal_flow_fluxes, &
                                        patch%boundary_flow_fluxes, &
                                        option%nflowdof, &
                                        instantaneous_array,option)
      array(istart:iend,1) = &
        integral_flux%integral_value(istart:iend)
      array(istart:iend,2) = &
        instantaneous_array(1:option%nflowdof)
    endif
    if (option%ntrandof > 0) then
      istart = option%nflowdof+1
      iend = option%nflowdof+option%ntrandof
      instantaneous_array = 0.d0
      call IntegralFluxGetInstantaneous(integral_flux, &
                                        patch%internal_tran_fluxes, &
                                        patch%boundary_tran_fluxes, &
                                        option%ntrandof, &
                                        instantaneous_array,option)
      array(istart:iend,1) = &
        integral_flux%integral_value(istart:iend)
      array(istart:iend,2) = &
        instantaneous_array(1:option%ntrandof)
    endif
    int_mpi = size(array)
    call MPI_Reduce(array,array_global, &
                    int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)
    ! time units conversion
    array_global(:,2) = array_global(:,2) * output_option%tconv
    if (option%myrank == option%io_rank) then
      if (option%nflowdof > 0) then
        do i = 1, option%nflowdof
          do j = 1, 2  ! 1 = integral, 2 = instantaneous
            tempreal = array_global(i,j)*flow_dof_scale(i)
            if (dabs(tempreal) > 0.d0 .and. dabs(tempreal) < 1.d-99) then
              write(fid,120,advance="no") tempreal
            else
              write(fid,110,advance="no") tempreal
            endif
          enddo
        enddo
      endif
      if (option%ntrandof > 0) then
        istart = option%nflowdof
        select case(option%itranmode)
          case(RT_MODE)
            do i=1,reaction%naqcomp
              do j = 1, 2  ! 1 = integral, 2 = instantaneous
                if (reaction%primary_species_print(i)) then
                  tempreal = array_global(istart+i,j)
                  if (dabs(tempreal) > 0.d0 .and. dabs(tempreal) < 1.d-99) then
                    write(fid,120,advance="no") tempreal
                  else
                    write(fid,110,advance="no") tempreal
                  endif
                endif
              enddo
            enddo
          case(NWT_MODE)
        end select
      endif
    endif
    integral_flux => integral_flux%next
  enddo
  deallocate(array)
  deallocate(array_global)
  deallocate(instantaneous_array)
  
  if (option%myrank == option%io_rank) then
    write(fid,'(a)') ''
    close(fid)
  endif
  
  integral_flux_first = PETSC_FALSE

end subroutine OutputIntegralFlux

! ************************************************************************** !

subroutine OutputMassBalance(realization_base)
  ! 
  ! Print to Tecplot POINT format
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/18/08
  ! 

  use Realization_Subsurface_class, only : realization_subsurface_type
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Utility_module
  use Output_Aux_module
  
  use Richards_module, only : RichardsComputeMassBalance
  use Mphase_module, only : MphaseComputeMassBalance
  use Flash2_module, only : Flash2ComputeMassBalance
  use Immis_module, only : ImmisComputeMassBalance
  use Miscible_module, only : MiscibleComputeMassBalance
  use TH_module, only : THComputeMassBalance
  use Reactive_Transport_module, only : RTComputeMassBalance
  use General_module, only : GeneralComputeMassBalance
  use Hydrate_module, only : HydrateComputeMassBalance
  use WIPP_Flow_module, only : WIPPFloComputeMassBalance
  use TOilIms_module, only : TOilImsComputeMassBalance
  use TOWG_module, only : TOWGComputeMassBalance
  use PM_TOilIms_Aux_module
  use PM_TOWG_Aux_module

  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use Material_Aux_class
  use General_Aux_module, only : general_fmw => fmw_comp
  use WIPP_Flow_Aux_module, only : wipp_flow_fmw => fmw_comp
  use Well_Data_class

  implicit none

  class(realization_base_type), target :: realization_base

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(output_option_type), pointer :: output_option
  type(coupler_type), pointer :: coupler
  type(mass_balance_region_type), pointer :: cur_mbr
  type(global_auxvar_type), pointer :: global_auxvars_bc_or_ss(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_bc_or_ss(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)

  class(reaction_rt_type), pointer :: reaction

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: word, units
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: fid = 86
  PetscInt :: ios
  PetscInt :: i,icol
  PetscInt :: k, j
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iconn
  PetscInt :: offset
  PetscInt :: iphase, ispec
  PetscInt :: icomp, nmobilecomp
  PetscInt :: max_tran_size
  PetscReal :: sum_area(4)
  PetscReal :: sum_area_global(4)
  PetscReal :: sum_kg(realization_base%option%nflowspec, &
               realization_base%option%nphase)
  PetscReal :: sum_kg_global(realization_base%option%nflowspec, &
               realization_base%option%nphase)
  PetscReal, allocatable :: sum_mol(:,:), sum_mol_global(:,:)
  
  PetscReal :: global_total_mass, global_water_mass

  PetscReal :: sum_trapped(realization_base%option%nphase)
  PetscReal :: sum_trapped_global(realization_base%option%nphase)
  PetscReal :: sum_mol_ye(3), sum_mol_global_ye(3)
  
  PetscMPIInt :: int_mpi
  PetscBool :: bcs_done
  PetscErrorCode :: ierr
  PetscBool,parameter :: wecl=PETSC_FALSE
  
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  reaction => ReactionCast(realization_base%reaction_base)
  output_option => realization_base%output_option

  if (option%ntrandof > 0) then
    select case(option%itranmode)
      case(RT_MODE)
        rt_auxvars => patch%aux%RT%auxvars
      case(NWT_MODE)
        option%io_buffer = 'OutputMassBalance has not been setup &
          &for NW Transport.'
        call PrintErrMsg(option)
    end select
  endif

  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '-mas.dat'
  else
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-mas.dat'
  endif
  
  ! open file
  if (option%myrank == option%io_rank) then

    if (output_option%print_column_ids) then
      icol = 1
    else
      icol = -1
    endif
  
    if (mass_balance_first .or. .not.FileExists(filename)) then
      open(unit=fid,file=filename,action="write",status="replace")

      ! write header
      write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // &
        ']"'  
      
      if (option%iflowmode > 0) then
        call OutputWriteToHeader(fid,'dt_flow',output_option%tunit,'',icol)
      endif
      
      if (option%ntrandof > 0) then
        call OutputWriteToHeader(fid,'dt_tran',output_option%tunit,'',icol)
      endif
      
      select case(option%iflowmode)
        case(RICHARDS_MODE,RICHARDS_TS_MODE)
          call OutputWriteToHeader(fid,'Global Water Mass','kg','',icol)
          
        case(TH_MODE,TH_TS_MODE)
          call OutputWriteToHeader(fid,'Global Water Mass in Liquid Phase', &
                                    'kg','',icol)
        case(G_MODE,H_MODE)
          call OutputWriteToHeader(fid,'Global Water Mass in Liquid Phase', &
                                    'kg','',icol)
          call OutputWriteToHeader(fid,'Global Air Mass in Liquid Phase', &
                                    'kg','',icol)
          call OutputWriteToHeader(fid,'Global Water Mass in Gas Phase', &
                                    'kg','',icol)
          call OutputWriteToHeader(fid,'Global Air Mass in Gas Phase', &
                                    'kg','',icol)
        case(WF_MODE)
          call OutputWriteToHeader(fid,'Global Water Mass in Liquid Phase', &
                                    'kg','',icol)
          call OutputWriteToHeader(fid,'Global Gas Component Mass in Gas &
                                   &Phase', 'kg','',icol)
        case(TOIL_IMS_MODE)
          call OutputWriteToHeader(fid,'Global Water Mass', &
                                    'kg','',icol)
          call OutputWriteToHeader(fid,'Global Oil Mass', &
                                    'kg','',icol)
        case(TOWG_MODE)
          call OutputWriteToHeader(fid,'Global Water Mass', &
                                    'kg','',icol)
          select case(towg_miscibility_model)
            case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF,TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
              call OutputWriteToHeader(fid,'Global Oil Mass', &
                                       'kg','',icol)
              call OutputWriteToHeader(fid,'Global Gas Mass', &
                                       'kg','',icol)
              if (towg_miscibility_model == TOWG_SOLVENT_TL) then
                call OutputWriteToHeader(fid,'Global Solvent Mass', &
                                         'kg','',icol)
              endif
          end select
        case(MPH_MODE,FLASH2_MODE)
          call OutputWriteToHeader(fid,'Global Water Mass in Water Phase', &
                                    'kmol','',icol)
          call OutputWriteToHeader(fid,'Global CO2 Mass in Water Phase', &
                                    'kmol','',icol)
          call OutputWriteToHeader(fid,'Trapped CO2 Mass in Water Phase', &
                                    'kmol','',icol)
          call OutputWriteToHeader(fid,'Global Water Mass in Gas Phase', &
                                    'kmol','',icol)
          call OutputWriteToHeader(fid,'Global CO2 Mass in Gas Phase', &
                                    'kmol','',icol)
          call OutputWriteToHeader(fid,'Trapped CO2 Mass in Gas Phase', &
                                    'kmol','',icol)
        case(IMS_MODE)
          call OutputWriteToHeader(fid,'Global Water Mass in Water Phase', &
                                    'kmol','',icol)
          call OutputWriteToHeader(fid,'Global CO2 Mass in Gas Phase', &
                                    'kmol','',icol)
        case(MIS_MODE)
          call OutputWriteToHeader(fid,'Global Water Mass in Liquid Phase', &
                                    'kg','',icol)
          call OutputWriteToHeader(fid,'Global Glycol Mass in Liquid Phase', &
                                    'kg','',icol)
      end select

      if (option%ntrandof > 0) then
        select case(option%itranmode)
          case(RT_MODE)
            do i=1,reaction%naqcomp
              if (reaction%primary_species_print(i)) then
                string = 'Global ' // trim(reaction%primary_species_names(i))
                call OutputWriteToHeader(fid,string,'mol','',icol)
              endif
            enddo

            do i=1,reaction%immobile%nimmobile
              if (reaction%immobile%print_me(i)) then
                string = 'Global ' // trim(reaction%immobile%names(i))
                call OutputWriteToHeader(fid,string,'mol','',icol)
              endif
            enddo

            do i=1,reaction%gas%nactive_gas
              if (reaction%gas%active_print_me(i)) then
                string = 'Global ' // trim(reaction%gas%active_names(i))
                call OutputWriteToHeader(fid,string,'mol','',icol)
              endif
            enddo

            if (option%mass_bal_detailed) then
              do i=1,reaction%mineral%nkinmnrl
                if (reaction%mineral%kinmnrl_print(i)) then
                  string = 'Global ' // trim(reaction%mineral%kinmnrl_names(i))
                  call OutputWriteToHeader(fid,string,'mol','',icol)
                endif
              enddo
            endif
          case(NWT_MODE)
        end select
      endif
      
      coupler => patch%boundary_condition_list%first
      bcs_done = PETSC_FALSE
      do
        if (.not.associated(coupler)) then
          if (bcs_done) then
            exit
          else
            bcs_done = PETSC_TRUE
            if (associated(patch%source_sink_list)) then
              coupler => patch%source_sink_list%first
              if (.not.associated(coupler)) exit
            else
              exit
            endif
          endif
        endif

        select case(option%iflowmode)
          case(RICHARDS_MODE,RICHARDS_TS_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            
            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
          case(TH_MODE,TH_TS_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            
            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
          case(MIS_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            string = trim(coupler%name) // ' Glycol Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            
            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
            string = trim(coupler%name) // ' Glycol Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
          case(G_MODE,H_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            string = trim(coupler%name) // ' Air Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)

            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
            string = trim(coupler%name) // ' Air Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
          case(WF_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            string = trim(coupler%name) // ' Gas Component Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)

            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
            string = trim(coupler%name) // ' Gas Component Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
          case(TOIL_IMS_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            string = trim(coupler%name) // ' Oil Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)

            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
            string = trim(coupler%name) // ' Oil Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
          case(TOWG_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            string = trim(coupler%name) // ' Oil Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            string = trim(coupler%name) // ' Gas Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            if (towg_miscibility_model == TOWG_SOLVENT_TL) then
              string = trim(coupler%name) // ' Solvent Mass'
              call OutputWriteToHeader(fid,string,'kg','',icol)          
            end if

            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
            string = trim(coupler%name) // ' Oil Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
            string = trim(coupler%name) // ' Gas Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
            if (towg_miscibility_model == TOWG_SOLVENT_TL) then
              string = trim(coupler%name) // ' Solvent Mass'
              call OutputWriteToHeader(fid,string,units,'',icol)
            end if
          case(MPH_MODE,FLASH2_MODE,IMS_MODE)
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,'kmol','',icol)
            string = trim(coupler%name) // ' CO2 Mass'
            call OutputWriteToHeader(fid,string,'kmol','',icol)
            
            units = 'kmol/' // trim(output_option%tunit) // ''
            string = trim(coupler%name) // ' Water Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
            string = trim(coupler%name) // ' CO2 Mass'
            call OutputWriteToHeader(fid,string,units,'',icol)
        end select
        
        if (option%ntrandof > 0) then
          select case(option%itranmode)
            case(RT_MODE)
              do i=1,reaction%naqcomp
                if (reaction%primary_species_print(i)) then
                  string = trim(coupler%name) // ' ' // &
                           trim(reaction%primary_species_names(i))
                  call OutputWriteToHeader(fid,string,'mol','',icol)
                endif
              enddo
  
              units = 'mol/' // trim(output_option%tunit) // ''
              do i=1,reaction%naqcomp
                if (reaction%primary_species_print(i)) then
                  string = trim(coupler%name) // ' ' // &
                           trim(reaction%primary_species_names(i))
                  call OutputWriteToHeader(fid,string,units,'',icol)
                endif
              enddo
            case(NWT_MODE)
          end select
        endif
        coupler => coupler%next
      
      enddo
      
      ! Print the water mass [kg] and species mass [mol] in the specified regions (header)
      if (associated(output_option%mass_balance_region_list)) then
        cur_mbr => output_option%mass_balance_region_list
        do
          if (.not.associated(cur_mbr)) exit
          string = 'Region ' // trim(cur_mbr%region_name) // ' Water Mass'
          call OutputWriteToHeader(fid,string,'kg','',icol)
          if (option%ntrandof > 0) then
            string = 'Region ' // trim(cur_mbr%region_name) // ' Total Mass'
            call OutputWriteToHeader(fid,string,'kg','',icol)
          endif
          cur_mbr => cur_mbr%next
        enddo
      endif

!  Write out well rates and total headers if required

      if (WellDataGetFlag()) then
        if (option%iflowmode == TOIL_IMS_MODE &
            .or. option%iflowmode == TOWG_MODE) then
          select type(realization_base)
           class is(realization_subsurface_type)
             call WriteWellHeaders(fid, icol, &
                                   realization_base, towg_miscibility_model, &
                                   option, wecl)
             if (output_option%write_masses) then
               call WriteWellMassHeaders(fid, icol, &
                                         realization_base, towg_miscibility_model)
             endif
          end select
        endif
      endif
      write(fid,'(a)') '' 
    else
      open(unit=fid,file=filename,action="write",status="old",position="append")
    endif 
    
  endif     

100 format(100es16.8)
110 format(100es16.8)

  ! write time
  if (option%myrank == option%io_rank) then
    write(fid,100,advance="no") option%time/output_option%tconv
  endif
  
  if (option%nflowdof > 0) then
    if (option%myrank == option%io_rank) &
      write(fid,100,advance="no") option%flow_dt/output_option%tconv
  endif
  if (option%ntrandof > 0) then
    if (option%myrank == option%io_rank) &
      write(fid,100,advance="no") option%tran_dt/output_option%tconv
  endif
  
! print out global mass balance

  if (option%nflowdof > 0) then
    sum_kg = 0.d0
    sum_trapped = 0.d0
    select type(realization_base)
      class is(realization_subsurface_type)
        select case(option%iflowmode)
          case(RICHARDS_MODE,RICHARDS_TS_MODE)
            call RichardsComputeMassBalance(realization_base,sum_kg(1,:))
          case(TH_MODE,TH_TS_MODE)
            call THComputeMassBalance(realization_base,sum_kg(1,:))
          case(MIS_MODE)
            call MiscibleComputeMassBalance(realization_base,sum_kg(:,1))
          case(MPH_MODE)
            call MphaseComputeMassBalance(realization_base,sum_kg(:,:),sum_trapped(:))
          case(FLASH2_MODE)
            call Flash2ComputeMassBalance(realization_base,sum_kg(:,:),sum_trapped(:))
          case(IMS_MODE)
            call ImmisComputeMassBalance(realization_base,sum_kg(:,1))
          case(G_MODE)
            call GeneralComputeMassBalance(realization_base,sum_kg(:,:))
          case(H_MODE)
            call HydrateComputeMassBalance(realization_base,sum_kg(:,:))
          case(WF_MODE)
            call WIPPFloComputeMassBalance(realization_base,sum_kg(:,1))
          case(TOIL_IMS_MODE)
            call TOilImsComputeMassBalance(realization_base,sum_kg(:,:))
          case(TOWG_MODE)
            call TOWGComputeMassBalance(realization_base,sum_kg(:,:))
        end select
      class default
        option%io_buffer = 'Unrecognized realization class in MassBalance().'
        call PrintErrMsg(option)
    end select

    int_mpi = option%nflowspec*option%nphase
    call MPI_Reduce(sum_kg,sum_kg_global, &
                    int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)

    if (option%iflowmode == MPH_MODE .or. option%iflowmode == FLASH2_MODE) then
!     call MPI_Barrier(option%mycomm,ierr)
      int_mpi = option%nphase
      call MPI_Reduce(sum_trapped,sum_trapped_global, &
                    int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)
    endif

    if (option%myrank == option%io_rank) then
      select case(option%iflowmode)
        case(RICHARDS_MODE,RICHARDS_TS_MODE,IMS_MODE,MIS_MODE,G_MODE,H_MODE, &
             TH_MODE,TH_TS_MODE)
          do iphase = 1, option%nphase
            do ispec = 1, option%nflowspec
              write(fid,110,advance="no") sum_kg_global(ispec,iphase)
            enddo
          enddo
        case(WF_MODE)
          do iphase = 1, option%nphase
            write(fid,110,advance="no") sum_kg_global(iphase,1)
          enddo
        case(TOIL_IMS_MODE)
          do iphase = 1, option%nphase
              write(fid,110,advance="no") sum_kg_global(iphase,1)
          enddo
        case(TOWG_MODE)
          select case(towg_miscibility_model)
            case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF,TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
              do iphase = 1, option%nphase
                write(fid,110,advance="no") sum_kg_global(iphase,1)
              enddo
          end select
        case(MPH_MODE,FLASH2_MODE)
          do iphase = 1, option%nphase
            do ispec = 1, option%nflowspec
              write(fid,110,advance="no") sum_kg_global(ispec,iphase)
            enddo
            write(fid,110,advance="no") sum_trapped_global(iphase)
          enddo
      end select
    endif
  endif
  
  if (option%ntrandof > 0) then
    if (option%transport%nphase > 1) then
      !TODO(geh): Within RTComputeMassBalance() all the mass is lumped into the
      !           liquid phase.  This need to be split out.  Also, the below
      !           where mass is summed across minerals needs to be moved
      !           to reactive_transport.F90.
      option%io_buffer = 'OutputMassBalance() needs to be refactored to &
        &consider species in the gas phase.'
!      call PrintErrMsg(option)
    endif
    max_tran_size = max(reaction%naqcomp,reaction%mineral%nkinmnrl, &
                        reaction%immobile%nimmobile,reaction%gas%nactive_gas)
    ! see RTComputeMassBalance for indexing used below
    allocate(sum_mol(max_tran_size,8))
    allocate(sum_mol_global(max_tran_size,8))
    sum_mol = 0.d0
    select type(realization_base)
      class is(realization_subsurface_type)
        call RTComputeMassBalance(realization_base,max_tran_size,sum_mol)
      class default
        option%io_buffer = 'Unrecognized realization class in MassBalance().'
        call PrintErrMsg(option)
    end select
    int_mpi = max_tran_size*8
    call MPI_Reduce(sum_mol,sum_mol_global,int_mpi, &
                    MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)

    if (option%myrank == option%io_rank) then
      do icomp = 1, reaction%naqcomp
        if (reaction%primary_species_print(icomp)) then
          write(fid,110,advance="no") sum_mol_global(icomp,1)
        endif
      enddo
      ! immobile species
      do i = 1, reaction%immobile%nimmobile
        if (reaction%immobile%print_me(i)) then
          write(fid,110,advance="no") &
            sum_mol_global(i,7)
        endif
      enddo
      ! gas species
      do i = 1, reaction%gas%nactive_gas
        if (reaction%gas%active_print_me(i)) then
          write(fid,110,advance="no") &
            sum_mol_global(i,8)
        endif
      enddo
    endif

!   print out mineral contribution to mass balance
    if (option%mass_bal_detailed) then
      if (option%myrank == option%io_rank) then
        do i = 1, reaction%mineral%nkinmnrl
          if (reaction%mineral%kinmnrl_print(i)) then
            write(fid,110,advance="no") sum_mol_global(i,6)
          endif
        enddo
      endif
    endif
    deallocate(sum_mol,sum_mol_global)
  endif

  coupler => patch%boundary_condition_list%first
  global_auxvars_bc_or_ss => patch%aux%Global%auxvars_bc
  if (option%ntrandof > 0) then
    select case(option%itranmode)
      case(RT_MODE)
        rt_auxvars_bc_or_ss => patch%aux%RT%auxvars_bc
      case(NWT_MODE)
    end select
  endif    
  bcs_done = PETSC_FALSE
  do 
    if (.not.associated(coupler)) then
      if (bcs_done) then
        exit
      else
        bcs_done = PETSC_TRUE
        if (associated(patch%source_sink_list)) then
          coupler => patch%source_sink_list%first
          if (.not.associated(coupler)) exit
          global_auxvars_bc_or_ss => patch%aux%Global%auxvars_ss
          if (option%ntrandof > 0) then
            select case(option%itranmode)
              case(RT_MODE)
                rt_auxvars_bc_or_ss => patch%aux%RT%auxvars_ss
              case(NWT_MODE)
            end select
          endif    
        else
          exit
        endif
      endif
    endif

    offset = coupler%connection_set%offset
    
    if (option%nflowdof > 0) then

#if 0
! compute the total area of the boundary condition
      if (.not.bcs_done) then
        sum_area = 0.d0
        do iconn = 1, coupler%connection_set%num_connections
          sum_area(1) = sum_area(1) + &
            coupler%connection_set%area(iconn)
          if (global_auxvars_bc_or_ss(offset+iconn)%sat(1) >= 0.5d0) then
            sum_area(2) = sum_area(2) + &
              coupler%connection_set%area(iconn)
          endif
          if (global_auxvars_bc_or_ss(offset+iconn)%sat(1) > 0.99d0) then
            sum_area(3) = sum_area(3) + &
              coupler%connection_set%area(iconn)
          endif
          sum_area(4) = sum_area(4) + &
            coupler%connection_set%area(iconn)* &
            global_auxvars_bc_or_ss(offset+iconn)%sat(1)
        enddo

        call MPI_Reduce(sum_area,sum_area_global, &
                        FOUR_INTEGER_MPI,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        option%io_rank,option%mycomm,ierr)
                          
        if (option%myrank == option%io_rank) then
          print *
          write(word,'(es16.6)') sum_area_global(1)
          print *, 'Total area in ' // trim(coupler%name) // &
                   ' boundary condition: ' // trim(adjustl(word)) // ' m^2'
          write(word,'(es16.6)') sum_area_global(2)
          print *, 'Total half-saturated area in '// &
                   trim(coupler%name) // &
                   ' boundary condition: ' // trim(adjustl(word)) // ' m^2'
          write(word,'(es16.6)') sum_area_global(3)
          print *, 'Total saturated area in '// trim(coupler%name) // &
                   ' boundary condition: ' // trim(adjustl(word)) // ' m^2'
          write(word,'(es16.6)') sum_area_global(4)
          print *, 'Total saturation-weighted area [=sum(saturation*area)] in '//&
                     trim(coupler%name) // &
                   ' boundary condition: ' // trim(adjustl(word)) // ' m^2'
          print *
        endif
      endif
#endif

      select case(option%iflowmode)
        case(RICHARDS_MODE,RICHARDS_TS_MODE)
          ! print out cumulative H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_auxvars_bc_or_ss(offset+iconn)%mass_balance
          enddo

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg,sum_kg_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global
          endif

          ! print out H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta
          enddo
          
          ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
          sum_kg = sum_kg*FMWH2O

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg,sum_kg_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global*output_option%tconv
          endif

        case(TH_MODE,TH_TS_MODE)
          ! print out cumulative H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_auxvars_bc_or_ss(offset+iconn)%mass_balance
          enddo

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg,sum_kg_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global
          endif

          ! print out H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta
          enddo
          ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
          sum_kg = sum_kg*FMWH2O

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg,sum_kg_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global*output_option%tconv
          endif

        case(MIS_MODE)
          ! print out cumulative mixture flux
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_auxvars_bc_or_ss(offset+iconn)%mass_balance(icomp,1)
            enddo
            
            if (icomp == 1) then
              sum_kg(icomp,1) = sum_kg(icomp,1)*FMWH2O
            else
              sum_kg(icomp,1) = sum_kg(icomp,1)*FMWGLYC
            endif
            
            int_mpi = option%nphase
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
          
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)
            endif
          enddo

          ! print out mixture flux
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta(icomp,1)
            enddo
            
        !   mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o/glycol
            if (icomp == 1) then
              sum_kg(icomp,1) = sum_kg(icomp,1)*FMWH2O
            else
              sum_kg(icomp,1) = sum_kg(icomp,1)*FMWGLYC
            endif

            int_mpi = option%nphase
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)*output_option%tconv
            endif
          enddo

        case(MPH_MODE,FLASH2_MODE)
        ! print out cumulative H2O & CO2 fluxes in kmol and kmol/time
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_auxvars_bc_or_ss(offset+iconn)%mass_balance(icomp,1)
            enddo
!geh            int_mpi = option%nphase
            int_mpi = 1
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)
            endif
          enddo
          
        ! print out H2O & CO2 fluxes in kmol and kmol/time
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta(icomp,1)
            enddo

          ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
!           sum_kg(icomp,1) = sum_kg(icomp,1)*FMWH2O ! <<---fix for multiphase!

!geh            int_mpi = option%nphase
            int_mpi = 1
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)*output_option%tconv
            endif
          enddo

        case(IMS_MODE)
        ! print out cumulative H2O & CO2 fluxes
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_auxvars_bc_or_ss(offset+iconn)%mass_balance(icomp,1)
            enddo
!geh            int_mpi = option%nphase
            int_mpi = 1
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)
            endif
          enddo
          
        ! print out H2O & CO2 fluxes
          sum_kg = 0.d0
          do icomp = 1, option%nflowspec
            do iconn = 1, coupler%connection_set%num_connections
              sum_kg(icomp,1) = sum_kg(icomp,1) + &
                global_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta(icomp,1)
            enddo

          ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
!           sum_kg(icomp,1) = sum_kg(icomp,1)*FMWH2O ! <<---fix for multiphase!

!geh            int_mpi = option%nphase
            int_mpi = 1
            call MPI_Reduce(sum_kg(icomp,1),sum_kg_global(icomp,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
            if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
              write(fid,110,advance="no") -sum_kg_global(icomp,1)*output_option%tconv
            endif
          enddo
        case(G_MODE,H_MODE)
          ! print out cumulative H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_auxvars_bc_or_ss(offset+iconn)%mass_balance
          enddo

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg(:,1),sum_kg_global(:,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global(:,1)
          endif

          ! print out H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg(:,1) = sum_kg(:,1) + &
              global_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta(:,1)
          enddo
          sum_kg(1,1) = sum_kg(1,1)*FMWH2O
          sum_kg(2,1) = sum_kg(2,1)*general_fmw(2)
          
          int_mpi = option%nphase
          call MPI_Reduce(sum_kg(:,1),sum_kg_global(:,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global(:,1)*output_option%tconv
          endif
        case(WF_MODE)
          ! print out cumulative H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_auxvars_bc_or_ss(offset+iconn)%mass_balance
          enddo

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg(:,1),sum_kg_global(:,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global(:,1)
          endif

          ! print out H2O flux
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg(:,1) = sum_kg(:,1) + &
              global_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta(:,1)
          enddo
          sum_kg(1,1) = sum_kg(1,1)*FMWH2O
          sum_kg(2,1) = sum_kg(2,1)*wipp_flow_fmw(2)
          
          int_mpi = option%nphase
          call MPI_Reduce(sum_kg(:,1),sum_kg_global(:,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global(:,1)*output_option%tconv
          endif
        case(TOIL_IMS_MODE)
          ! print out cumulative H2O and Oil fluxes
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_auxvars_bc_or_ss(offset+iconn)%mass_balance
          enddo

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg(:,1),sum_kg_global(:,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global(:,1)
          endif

          ! print out H2O and oil fluxes
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg(:,1) = sum_kg(:,1) + &
              global_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta(:,1)
          enddo
          sum_kg(1,1) = sum_kg(1,1)*toil_ims_fmw_comp(1) 
          sum_kg(2,1) = sum_kg(2,1)*toil_ims_fmw_comp(2)
          
          int_mpi = option%nphase
          call MPI_Reduce(sum_kg(:,1),sum_kg_global(:,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global(:,1)*output_option%tconv
          endif
        case(TOWG_MODE)
          ! print out cumulative Water, Oil, Gas and Solvent fluxes
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg = sum_kg + global_auxvars_bc_or_ss(offset+iconn)%mass_balance
          enddo

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg(:,1),sum_kg_global(:,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global(:,1)
          endif

          ! print out H2O, oil and gas fluxes
          sum_kg = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_kg(:,1) = sum_kg(:,1) + &
              global_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta(:,1)
          enddo
          sum_kg(1,1) = sum_kg(1,1)*towg_fmw_comp(1) 
          sum_kg(2,1) = sum_kg(2,1)*towg_fmw_comp(2)
          sum_kg(3,1) = sum_kg(3,1)*towg_fmw_comp(3)
          if (towg_miscibility_model == TOWG_SOLVENT_TL) then
            sum_kg(4,1) = sum_kg(4,1)*towg_fmw_comp(4)
          endif

          int_mpi = option%nphase
          call MPI_Reduce(sum_kg(:,1),sum_kg_global(:,1), &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)

          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            write(fid,110,advance="no") -sum_kg_global(:,1)*output_option%tconv
          endif
      end select
    endif
    
    if (option%ntrandof > 0) then
      select case(option%itranmode)
        case(RT_MODE)
          nmobilecomp = reaction%naqcomp
          allocate(sum_mol(nmobilecomp,option%transport%nphase))
          allocate(sum_mol_global(nmobilecomp,option%transport%nphase))
          ! print out cumulative boundary flux
          sum_mol = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_mol = sum_mol + &
              rt_auxvars_bc_or_ss(offset+iconn)%mass_balance(1:nmobilecomp,:)
          enddo

          int_mpi = nmobilecomp
          call MPI_Reduce(sum_mol,sum_mol_global,int_mpi, &
                          MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)

          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            do icomp = 1, reaction%naqcomp
              if (reaction%primary_species_print(icomp)) then
                write(fid,110,advance="no") -sum_mol_global(icomp,1)
              endif
            enddo
          endif

          ! print out boundary flux
          sum_mol = 0.d0
          do iconn = 1, coupler%connection_set%num_connections
            sum_mol = sum_mol + rt_auxvars_bc_or_ss(offset+iconn)% &
                                  mass_balance_delta(1:nmobilecomp,:) 
          enddo

          int_mpi = nmobilecomp
          call MPI_Reduce(sum_mol,sum_mol_global,int_mpi, &
                          MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)

          if (option%myrank == option%io_rank) then
            ! change sign for positive in / negative out
            do icomp = 1, reaction%naqcomp
              if (reaction%primary_species_print(icomp)) then
                write(fid,110,advance="no") -sum_mol_global(icomp,1)* &
                                              output_option%tconv
              endif
            enddo
          endif
        case(NWT_MODE)
      end select
      deallocate(sum_mol,sum_mol_global)
    endif

    coupler => coupler%next 
  enddo
  
  ! Print the total water and component mass in the specified regions (data)
  if (associated(output_option%mass_balance_region_list)) then
    cur_mbr => output_option%mass_balance_region_list
    do
      if (.not.associated(cur_mbr)) exit
      call PatchGetWaterMassInRegion(cur_mbr%region_cell_ids, &
                                     cur_mbr%num_cells,patch,option, &
                                     global_water_mass)
      write(fid,110,advance="no") global_water_mass
      if (option%ntrandof > 0) then
        call PatchGetCompMassInRegion(cur_mbr%region_cell_ids, &
             cur_mbr%num_cells,patch,option,global_total_mass)
        write(fid,110,advance="no") global_total_mass
      endif
      cur_mbr => cur_mbr%next
    enddo
  endif

!  Write out well rates and totals if required

  if (WellDataGetFlag()) then
    if (     option%iflowmode == TOIL_IMS_MODE &
        .or. option%iflowmode == TOWG_MODE       ) then
    if (option%myrank == option%io_rank) then
        select type(realization_base)
         class is(realization_subsurface_type)
          call WriteWellValues(fid,realization_base, &
                               towg_miscibility_model, &
                               option, wecl, sum_kg_global)
          if (output_option%write_masses) then
            call WriteWellMassValues(fid,realization_base, &
                                     towg_miscibility_model)
          endif
        end select
      endif
    endif
  endif

  if (option%myrank == option%io_rank) then
    write(fid,'(a)') ''
    close(fid)
  endif
  
  mass_balance_first = PETSC_FALSE

end subroutine OutputMassBalance

! ************************************************************************** !

subroutine OutputEclipseFiles(realization_base)
  !
  ! Write out Eclipse spec, summary and restart files
  !
  ! Author: Dave Ponting
  ! Date: 01/10/19
  !

  use Realization_Subsurface_class, only : realization_subsurface_type
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Utility_module
  use Output_Aux_module
  use PM_TOWG_Aux_module, only: towg_miscibility_model
  use Grid_Grdecl_module, only : GetIsGrdecl
  use Well_Data_class
  use TOilIms_module, only : TOilImsComputeMassBalance
  use TOWG_module, only : TOWGComputeMassBalance

  implicit none

  class(realization_base_type), target :: realization_base
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(output_option_type), pointer :: output_option

  PetscInt :: fid, icol
  PetscBool, parameter :: wecl = PETSC_TRUE
  PetscBool :: write_summ, write_rest,is_grdecl
  PetscInt  :: sum_ds, rst_ds, sum_ls, rst_ls
  PetscReal :: sum_dt, rst_dt, sum_lt, rst_lt, time
  PetscReal, parameter :: eps = 0.001

  PetscReal :: sum_kg(realization_base%option%nflowspec, &
               realization_base%option%nphase)
  PetscReal :: sum_kg_global(realization_base%option%nflowspec, &
               realization_base%option%nphase)
  PetscMPIInt :: int_mpi
  PetscErrorCode :: ierr

  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  output_option => realization_base%output_option

  !  Check that we have grid locations

  is_grdecl = GetIsGrdecl()

  if (.not.is_grdecl) then
    option%io_buffer = 'Eclipse file output requires grdecl type input'
    call PrintErrMsg(option)
  endif

  !  Find mass balance

  sum_kg = 0.d0
  select type(realization_base)
    class is(realization_subsurface_type)
      select case(option%iflowmode)
        case(TOIL_IMS_MODE)
          call TOilImsComputeMassBalance(realization_base,sum_kg(:,:))
        case(TOWG_MODE)
          call TOWGComputeMassBalance(realization_base,sum_kg(:,:))
      end select
    class default
      option%io_buffer = &
        'Unrecognized realization class in OutputEclipseFiles().'
      call PrintErrMsg(option)
  end select

  int_mpi = option%nflowspec*option%nphase
  call MPI_Reduce(sum_kg,sum_kg_global, &
                  int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                  option%io_rank,option%mycomm,ierr)

  !  Set useful scalars (negative fid prevents writing to -mas files)

  icol =  1
  fid  = -1
  time = option%time

  sum_dt = output_option%eclipse_options%write_ecl_sum_deltat
  rst_dt = output_option%eclipse_options%write_ecl_rst_deltat
  sum_ds = output_option%eclipse_options%write_ecl_sum_deltas
  rst_ds = output_option%eclipse_options%write_ecl_rst_deltas

  sum_lt  = output_option%eclipse_options%write_ecl_sum_lastt
  rst_lt  = output_option%eclipse_options%write_ecl_rst_lastt
  sum_ls  = output_option%eclipse_options%write_ecl_sum_lasts
  rst_ls  = output_option%eclipse_options%write_ecl_rst_lasts

  write_summ=GetEclWrtFlg(ewriter_summ_count, &
                          time, sum_dt, sum_ds, sum_lt, sum_ls)
  write_rest=GetEclWrtFlg(ewriter_rest_count, &
                          time, rst_dt, rst_ds, rst_lt, rst_ls)

  output_option%eclipse_options%write_ecl_sum_lastt = sum_lt
  output_option%eclipse_options%write_ecl_rst_lastt = rst_lt
  output_option%eclipse_options%write_ecl_sum_lasts = sum_ls
  output_option%eclipse_options%write_ecl_rst_lasts = rst_ls

  ! Summary files - just needs the io rank

  if (write_summ) then
    if (option%myrank == option%io_rank) then
      if (     option%iflowmode == TOIL_IMS_MODE &
          .or. option%iflowmode == TOWG_MODE      ) then

        if (ewriter_summ_count == 0) then

  !  Write out well and field headers if required

          select type(realization_base)
           class is(realization_subsurface_type)
             call WriteWellHeaders(fid, icol, &
                                   realization_base, &
                                   towg_miscibility_model, &
                                   option, wecl)
          end select
        endif

  !  Write out well rates and totals if required

        select type(realization_base)
         class is(realization_subsurface_type)
           call WriteWellValues(fid, realization_base, &
                                towg_miscibility_model, &
                                option, wecl, sum_kg_global)
        end select

      endif ! IMS or TOWG
    endif ! IO rank
  endif

  ! Restart files - needs all the ranks

  if (write_rest) then
    if (ewriter_rest_count == 0) then
      call setupEwriterRestMaps(patch, grid, option)
    endif
    select type(realization_base)
     class is(realization_subsurface_type)
     call WriteRestValues(realization_base, option)
    end select
  endif

  !  Set flags indicating first write operations done

  ewriter_summ_count = ewriter_summ_count+1
  ewriter_rest_count = ewriter_rest_count+1

end subroutine OutputEclipseFiles

! *************************************************************************** !

subroutine WriteWellHeaders(fid, icol, realization, &
                            towg_miscibility_model, &
                            option, wecl)
  !
  ! Used to write out file headers specific to TOIL and TOWG modes
  ! This routine must match the headers written by write_well_values

  ! Author: Dave Ponting
  ! Date  : 09/15/18

  use Realization_Subsurface_class
  use Well_Data_class
  use Option_module
  use Output_Eclipse_module, only:WriteEclipseFilesSpec

  implicit none

  PetscInt, intent(in   ) :: fid
  PetscInt, intent(inout) :: icol
  type(realization_subsurface_type) :: realization
  PetscInt, intent(in   ) :: towg_miscibility_model
  type(option_type), intent(in), pointer :: option
  PetscBool, intent(in) :: wecl

  type(well_data_list_type), pointer :: well_data_list

  PetscInt :: iwell, nwell, ni, mi, iword, ichar, irfn, lrfn

  character(len=MAXSTRINGLENGTH) :: name
  PetscBool :: is_restart
  character(len=8) :: restart_filename(9)
  character(len=1) :: z1

  character(len=8), allocatable :: zm(:)
  character(len=8), allocatable :: zn(:)
  character(len=8), allocatable :: zu(:)

  ! Write out Eclipse files if required

  ni = 0
  mi = 0
  if (wecl) then
    mi = 1
    ni = 1
    allocate(zm(mi))
    allocate(zn(mi))
    allocate(zu(mi))
    zm = ' '
    zn = ' '
    zu = ' '
    zm(1) = 'TIME'
    zn(1) = ' '
    zu(1) = 'DAYS'
  endif

  !  Find well list and loop over wells

  well_data_list => realization%well_data
  nwell = getnwell(well_data_list)
  do iwell = 1, nwell

  ! Get name and type of this well

    call getWellNameI(iwell, well_data_list, name)

  ! Oil rates and totals

    call WrtHrd(fid, 'wopr', name, 'm^3/d' , icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wopt', name, 'm^3'   , icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'woir', name, 'm^3/d' , icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'woit', name, 'm^3'   , icol, zm, zn, zu, ni, mi, wecl)

  ! Gas rates and totals

    call WrtHrd(fid, 'wgpr', name, 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wgpt', name, 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wgir', name, 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wgit', name, 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)

  ! Water rates and totals

    call WrtHrd(fid, 'wwpr', name, 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wwpt', name, 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wwir', name, 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wwit', name, 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)

  ! Solvent rates and totals if required

    if (towg_miscibility_model == TOWG_SOLVENT_TL) then
      if (wecl) then
        ! Eclipse uses n for solvent as s is salt
        call WrtHrd(fid, 'wnpr', name, 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
        call WrtHrd(fid, 'wnpt', name, 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
        call WrtHrd(fid, 'wnir', name, 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
        call WrtHrd(fid, 'wnit', name, 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
      else
        call WrtHrd(fid, 'wspr', name, 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
        call WrtHrd(fid, 'wspt', name, 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
        call WrtHrd(fid, 'wsir', name, 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
        call WrtHrd(fid, 'wsit', name, 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
      endif
    endif

  ! Liquid rates and totals

    call WrtHrd(fid, 'wlpr', name, 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wlpt', name, 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wgor', name, 'bar'  , icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wwct', name, 'bar'  , icol, zm, zn, zu, ni, mi, wecl)
    call WrtHrd(fid, 'wbhp', name, 'bar'  , icol, zm, zn, zu, ni, mi, wecl)

  enddo

  call WrtHrd(fid, 'fopr', 'field', 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fopt', 'field', 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'foir', 'field', 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'foit', 'field', 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)

  call WrtHrd(fid, 'fgpr', 'field', 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fgpt', 'field', 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fgir', 'field', 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fgit', 'field', 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)

  call WrtHrd(fid, 'fwpr', 'field', 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fwpt', 'field', 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fwir', 'field', 'm^3/d', icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fwit', 'field', 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)

  if (towg_miscibility_model == TOWG_SOLVENT_TL) then
    if (wecl) then
      ! Eclipse uses n for solvent as s is salt
      call WrtHrd(fid, 'fnpr', 'field', 'm^3/d', &
                  icol, zm, zn, zu, ni, mi, wecl)
      call WrtHrd(fid, 'fnpt', 'field', 'm^3'  , &
                  icol, zm, zn, zu, ni, mi, wecl)
      call WrtHrd(fid, 'fnir', 'field', 'm^3/d', &
                  icol, zm, zn, zu, ni, mi, wecl)
      call WrtHrd(fid, 'fnit', 'field', 'm^3'  , &
                  icol, zm, zn, zu, ni, mi, wecl)
    else
      call WrtHrd(fid, 'fspr', 'field', 'm^3/d', &
                  icol, zm, zn, zu, ni, mi, wecl)
      call WrtHrd(fid, 'fspt', 'field', 'm^3'  , &
                  icol, zm, zn, zu, ni, mi, wecl)
      call WrtHrd(fid, 'fsir', 'field', 'm^3/d', &
                  icol, zm, zn, zu, ni, mi, wecl)
      call WrtHrd(fid, 'fsit', 'field', 'm^3'  , &
                 icol, zm, zn, zu, ni, mi, wecl)
    endif
  endif

  call WrtHrd(fid, 'flpr', 'field', 'm^3/d'  , icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'flpt', 'field', 'm^3'    , icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fgor', 'field', 'm^3/m^3', icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fwct', 'field', 'm^3/m^3', icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fpav', 'field', 'Bar    ', icol, zm, zn, zu, ni, mi, wecl)

  call WrtHrd(fid, 'foip', 'field', 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fgip', 'field', 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
  call WrtHrd(fid, 'fwip', 'field', 'm^3'  , icol, zm, zn, zu, ni, mi, wecl)
  if (towg_miscibility_model == TOWG_SOLVENT_TL) then
    if (wecl) then
      ! Eclipse uses n for solvent as s is salt
      call WrtHrd(fid, 'fnip', 'field', 'm^3', icol, zm, zn, zu, ni, mi, wecl)
    else
      call WrtHrd(fid, 'fsip', 'field', 'm^3', icol, zm, zn, zu, ni, mi, wecl)
    endif
  endif

  ! Write out Eclipse files if required

  if (wecl) then
    is_restart = PETSC_FALSE
    restart_filename = ' '
    lrfn = len(trim(option%restart_filename))
    if (lrfn > 0) then
      is_restart = PETSC_TRUE
      do irfn = 1, lrfn
        z1=option%restart_filename(irfn:irfn)
        if (z1 == '-') exit
        iword = (irfn-1)/8 + 1
        ichar = irfn - 8*(iword-1)
        restart_filename(iword)(ichar:ichar) = z1
      enddo
    endif
    call WriteEclipseFilesSpec(zm, zn, zu, ni, is_restart, restart_filename)
    deallocate(zm)
    deallocate(zn)
    deallocate(zu)
  endif

end subroutine WriteWellHeaders

! *************************************************************************** !

subroutine WriteWellMassHeaders(fid, icol, realization, towg_miscibility_model)
  !
  ! Used to write out file headers specific to TOIL and TOWG modes
  ! This routine must match the headers written by write_well_values

  ! Author: Dave Ponting
  ! Date  : 09/15/18

  use Realization_Subsurface_class
  use Well_Data_class

  implicit none

  PetscInt, intent(in   ) :: fid
  PetscInt, intent(inout) :: icol
  type(realization_subsurface_type) :: realization
  PetscInt, intent(in   ) :: towg_miscibility_model

  type(well_data_list_type), pointer :: well_data_list

  PetscInt :: iwell, nwell
  character(len=MAXSTRINGLENGTH) :: name

  !  Find well list and loop over wells

  well_data_list => realization%well_data
  nwell = getnwell(well_data_list)
  do iwell = 1, nwell

  ! Get name and type of this well

    call getWellNameI(iwell, well_data_list, name)

  ! Oil mass rates and totals

    call WrtHrdMO(fid, 'wompr', name, 'kg/d', icol)
    call WrtHrdMO(fid, 'wompt', name, 'kg'  , icol)
    call WrtHrdMO(fid, 'womir', name, 'kg/d', icol)
    call WrtHrdMO(fid, 'womit', name, 'kg'  , icol)

  ! Gas mass rates and totals

    call WrtHrdMO(fid, 'wgmpr' , name, 'kg/d', icol)
    call WrtHrdMO(fid, 'wgmpt' , name, 'kg'  , icol)
    call WrtHrdMO(fid, 'wgmir' , name, 'kg/d', icol)
    call WrtHrdMO(fid, 'wgmit' , name, 'kg'  , icol)

  ! Water mass rates and totals

    call WrtHrdMO(fid, 'wwmpr', name, 'kg/d', icol)
    call WrtHrdMO(fid, 'wwmpt', name, 'kg'  , icol)
    call WrtHrdMO(fid, 'wwmir', name, 'kg/d', icol)
    call WrtHrdMO(fid, 'wwmit', name, 'kg'  , icol)

  ! Solvent mass rates and totals if required

    if (towg_miscibility_model == TOWG_SOLVENT_TL) then
      call WrtHrdMO(fid, 'wsmpr', name, 'kg/d', icol)
      call WrtHrdMO(fid, 'wsmpt', name, 'kg'  , icol)
      call WrtHrdMO(fid, 'wsmir', name, 'kg/d', icol)
      call WrtHrdMO(fid, 'wsmit', name, 'kg'  , icol)
    endif

  enddo

  call WrtHrdMO(fid, 'fompr', 'field', 'kg/d', icol)
  call WrtHrdMO(fid, 'fompt', 'field', 'kg'  , icol)
  call WrtHrdMO(fid, 'fomir', 'field', 'kg/d', icol)
  call WrtHrdMO(fid, 'fomit', 'field', 'kg'  , icol)

  call WrtHrdMO(fid, 'fgmpr', 'field', 'kg/d', icol)
  call WrtHrdMO(fid, 'fgmpt', 'field', 'kg'  , icol)
  call WrtHrdMo(fid, 'fgmir', 'field', 'kg/d', icol)
  call WrtHrdMO(fid, 'fgmit', 'field', 'kg'  , icol)

  call WrtHrdMO(fid, 'fwmpr', 'field', 'kg/d', icol)
  call WrtHrdMO(fid, 'fwmpt', 'field', 'kg'  , icol)
  call WrtHrdMO(fid, 'fwmir', 'field', 'kg/d', icol)
  call WrtHrdMO(fid, 'fwmit', 'field', 'kg'  , icol)

  if (towg_miscibility_model == TOWG_SOLVENT_TL) then
    call WrtHrdMO(fid, 'fsmpr', 'field', 'kg/d', icol)
    call WrtHrdMO(fid, 'fsmpt', 'field', 'kg'  , icol)
    call WrtHrdMO(fid, 'fsmir', 'field', 'kg/d', icol)
    call WrtHrdMO(fid, 'fsmit', 'field', 'kg'  , icol)
  endif

end subroutine WriteWellMassHeaders

! *************************************************************************** !

subroutine WriteWellValues(fid, realization, towg_miscibility_model, &
                           option, wecl, sum_kg_global)
  !
  ! Used to write out mas file values specific to TOIL and TOWG modes
  ! This routine must match the headers written by WriteWellHeaders
  !
  ! Author: Dave Ponting
  ! Date  : 09/15/18

  use Realization_Subsurface_class
  use Realization_Base_class, only : realization_base_type
  use Well_Type_class
  use Well_Data_class
  use Output_Eclipse_module, only:WriteEclipseFilesSumm
  use Option_module
  use EOS_Oil_module,  only : EOSOilGetSurfaceDensity
  use EOS_Gas_module,  only : EOSGasGetSurfaceDensity
  use EOS_Water_module,only : EOSWaterGetSurfaceDensity
  use EOS_Slv_module,  only : EOSSlvGetSurfaceDensity

  implicit none

  PetscInt, intent(in   ) :: fid
  type(realization_subsurface_type) :: realization
  type(option_type), intent(in), pointer :: option
  PetscBool, intent(in) :: wecl
  PetscReal :: sum_kg_global(:,:)

  PetscReal :: tconv, sign
  PetscInt :: towg_miscibility_model
  type(well_data_list_type), pointer :: well_data_list

  PetscInt :: iwell, nwell, welltype

  PetscReal :: wopriu, wgpriu, wwpriu, wspriu, &
               wopr, wopt, woir, woit, &
               wgpr, wgpt, wgir, wgit, &
               wwpr, wwpt, wwir, wwit, &
               wspr, wspt, wsir, wsit, &
               wlpr, wlpt, wgor, wwct, wbhp

  PetscReal :: fopr, fopt, foir, foit, &
               fgpr, fgpt, fgir, fgit, &
               fwpr, fwpt, fwir, fwit, &
               fspr, fspt, fsir, fsit, &
               flpr, flpt, fgor, fwct, fpav, &
               foip, fgip, fwip, fsip, sd

  PetscReal, pointer :: vd(:)
  PetscInt           :: nd, md, iphase

  tconv = 3600.0*24.0

  nd = 0

  !  Find well list and loop over wells

  well_data_list => realization%well_data
  nwell = getnwell(well_data_list)

  !  Set up array of Eclipse summary file data if required

  if (wecl) then
    nd = 1
    md = 1
    allocate(vd(md))
    vd = 0.0
    vd(nd) = option%time/tconv
  endif

  !  Loop over wells

  do iwell = 1, nwell

  !  Set up well flow sign (+ ve producers and -ve injectors)

    welltype = getWellTypeI(iwell, well_data_list)
    if (wellType == PROD_WELL_TYPE) then
      sign = 1.0
    else
      sign =-1.0
    endif

  !  Get values in internal units

    wopriu = getWellTTValI(iwell, W_TARG_OSV, VALTYPE_ACTUAL, well_data_list)
    wopt   = getWellTTValI(iwell, W_TARG_OSV, VALTYPE_TOTALP, well_data_list)
    woit   = getWellTTValI(iwell, W_TARG_OSV, VALTYPE_TOTALI, well_data_list)

    wgpriu = getWellTTValI(iwell, W_TARG_GSV, VALTYPE_ACTUAL, well_data_list)
    wgpt   = getWellTTValI(iwell, W_TARG_GSV, VALTYPE_TOTALP, well_data_list)
    wgit   = getWellTTValI(iwell, W_TARG_GSV, VALTYPE_TOTALI, well_data_list)

    wwpriu = getWellTTValI(iwell, W_TARG_WSV, VALTYPE_ACTUAL, well_data_list)
    wwpt   = getWellTTValI(iwell, W_TARG_WSV, VALTYPE_TOTALP, well_data_list)
    wwit   = getWellTTValI(iwell, W_TARG_WSV, VALTYPE_TOTALI, well_data_list)

    if (towg_miscibility_model == TOWG_SOLVENT_TL) then
      wspriu = getWellTTValI(iwell, W_TARG_SSV, VALTYPE_ACTUAL, well_data_list)
      wspt   = getWellTTValI(iwell, W_TARG_SSV, VALTYPE_TOTALP, well_data_list)
      wsit   = getWellTTValI(iwell, W_TARG_SSV, VALTYPE_TOTALI, well_data_list)
    else
      wspriu = 0.0
      wspt   = 0.0
      wsit   = 0.0
    endif

    wbhp = getWellTTValI(iwell, W_BHP_LIMIT, VALTYPE_ACTUAL, well_data_list)

  !  Convert rates to user units (per day not per sec) and sign convention

    if (wellType == PROD_WELL_TYPE) then
      wopr = wopriu*tconv*sign
      wgpr = wgpriu*tconv*sign
      wwpr = wwpriu*tconv*sign
      wspr = wspriu*tconv*sign
      woir = 0.0
      wgir = 0.0
      wwir = 0.0
      wsir = 0.0
    else
      wopr = 0.0
      wgpr = 0.0
      wwpr = 0.0
      wspr = 0.0
      woir = wopriu*tconv*sign
      wgir = wgpriu*tconv*sign
      wwir = wwpriu*tconv*sign
      wsir = wspriu*tconv*sign
    endif

  !  Convert BHP to Bars

    wbhp = wbhp*1.0D-5

  !  Set up dependent values (liquid rates and ratios)

    wlpr = wopr+wwpr
    wlpt = wopt+wwpt

    wwct = 0.0
    if (wlpr > 0.0) wwct = wwpr/wlpr
    wgor = 0.0
    if (wopr > 0.0) wgor = wgpr/wopr

  !  Write out well values

    call wrtToTableAndSumm(fid, wopr, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wopt, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, woir, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, woit, vd, nd, md, wecl)

    call wrtToTableAndSumm(fid, wgpr, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wgpt, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wgir, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wgit, vd, nd, md, wecl)

    call wrtToTableAndSumm(fid, wwpr, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wwpt, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wwir, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wwit, vd, nd, md, wecl)

    if (towg_miscibility_model == TOWG_SOLVENT_TL) then
      call wrtToTableAndSumm(fid, wspr, vd, nd, md, wecl)
      call wrtToTableAndSumm(fid, wspt, vd, nd, md, wecl)
      call wrtToTableAndSumm(fid, wsir, vd, nd, md, wecl)
      call wrtToTableAndSumm(fid, wsit, vd, nd, md, wecl)
    endif

    call wrtToTableAndSumm(fid, wlpr, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wlpt, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wgor, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wwct, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, wbhp, vd, nd, md, wecl)

  enddo

  !  Now the field values

  fopr = GetFieldTTVal(W_TARG_OSV, VALTYPE_ACTUALP, well_data_list)
  foir = GetFieldTTVal(W_TARG_OSV, VALTYPE_ACTUALI, well_data_list)
  fopt = GetFieldTTVal(W_TARG_OSV, VALTYPE_TOTALP , well_data_list)
  foit = GetFieldTTVal(W_TARG_OSV, VALTYPE_TOTALI , well_data_list)

  fgpr = GetFieldTTVal(W_TARG_GSV, VALTYPE_ACTUALP, well_data_list)
  fgir = GetFieldTTVal(W_TARG_GSV, VALTYPE_ACTUALI, well_data_list)
  fgpt = GetFieldTTVal(W_TARG_GSV, VALTYPE_TOTALP , well_data_list)
  fgit = GetFieldTTVal(W_TARG_GSV, VALTYPE_TOTALI , well_data_list)

  fwpr = GetFieldTTVal(W_TARG_WSV, VALTYPE_ACTUALP, well_data_list)
  fwir = GetFieldTTVal(W_TARG_WSV, VALTYPE_ACTUALI, well_data_list)
  fwpt = GetFieldTTVal(W_TARG_WSV, VALTYPE_TOTALP , well_data_list)
  fwit = GetFieldTTVal(W_TARG_WSV, VALTYPE_TOTALI , well_data_list)

  if (towg_miscibility_model == TOWG_SOLVENT_TL) then
    fspr = GetFieldTTVal(W_TARG_SSV, VALTYPE_ACTUALP, well_data_list)
    fsir = GetFieldTTVal(W_TARG_SSV, VALTYPE_ACTUALI, well_data_list)
    fspt = GetFieldTTVal(W_TARG_SSV, VALTYPE_TOTALP , well_data_list)
    fsit = GetFieldTTVal(W_TARG_SSV, VALTYPE_TOTALI , well_data_list)
  else
    fspr = 0.0
    fsir = 0.0
    fspt = 0.0
    fsit = 0.0
  endif

!  Convert to rates per day

  fopr = fopr*tconv
  foir = foir*tconv

  fgpr = fgpr*tconv
  fgir = fgir*tconv

  fwpr = fwpr*tconv
  fwir = fwir*tconv

  fspr = fspr*tconv
  fsir = fsir*tconv

!  Find liquid rates, water cut and gor

  flpr = fopr+fwpr
  flpt = fopt+fwpt

  fwct = 0.0
  fgor = 0.0
  if (flpr > 0.0) fwct = fwpr/flpr
  if (fopr > 0.0) fgor = fgpr/fopr

  !  Convert field pressure to Bars

  call GetFieldData(fpav)
  fpav = fpav*1.0D-5

  !  Get field mass-in-place values

  foip = 0.0
  fgip = 0.0
  fwip = 0.0
  fsip = 0.0

  do iphase = 1, option%nphase
    ! Oil phase
    if (iphase == option%oil_phase) then
      sd = EOSOilGetSurfaceDensity()
      if (sd>0.0) then
        foip = sum_kg_global(iphase,1)/sd
      endif
    endif
    ! Gas phase
    if (iphase == option%gas_phase) then
      sd = EOSGasGetSurfaceDensity()
      if (sd>0.0) then
        fgip = sum_kg_global(iphase,1)/sd
      endif
    endif
    !  Water (called liquid) phase
    if (iphase == option%liquid_phase) then
      sd = EOSWaterGetSurfaceDensity()
      if (sd>0.0) then
        fwip = sum_kg_global(iphase,1)/sd
      endif
    endif
    !  Solvent phase
    if (iphase == option%solvent_phase) then
      sd = EOSSlvGetSurfaceDensity()
      if (sd>0.0) then
        fsip = sum_kg_global(iphase,1)/sd
      endif
    endif
  enddo

  !  Write out field values

  call wrtToTableAndSumm(fid, fopr, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fopt, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, foir, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, foit, vd, nd, md, wecl)

  call wrtToTableAndSumm(fid, fgpr, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fgpt, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fgir, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fgit, vd, nd, md, wecl)

  call wrtToTableAndSumm(fid, fwpr, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fwpt, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fwir, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fwit, vd, nd, md, wecl)

  if (towg_miscibility_model == TOWG_SOLVENT_TL) then
    call wrtToTableAndSumm(fid, fspr, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, fspt, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, fsir, vd, nd, md, wecl)
    call wrtToTableAndSumm(fid, fsit, vd, nd, md, wecl)
  endif

  call wrtToTableAndSumm(fid, flpr, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, flpt, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fgor, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fwct, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fpav, vd, nd, md, wecl)

  call wrtToTableAndSumm(fid, foip, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fgip, vd, nd, md, wecl)
  call wrtToTableAndSumm(fid, fwip, vd, nd, md, wecl)
  if (towg_miscibility_model == TOWG_SOLVENT_TL) then
    call wrtToTableAndSumm(fid, fsip, vd, nd, md, wecl)
  endif

  ! Write out Eclipse files if required

  if (wecl) then
    call WriteEclipseFilesSumm(vd, nd)
    deallocate(vd)
   endif

end subroutine WriteWellValues

! *************************************************************************** !

subroutine WriteWellMassValues(fid, realization, towg_miscibility_model)
  !
  ! Used to write out mass rates to the mas file values
  ! Specific to TOIL and TOWG modes
  ! This routine must match the headers written by WriteWellMassHeaders
  !
  ! Author: Dave Ponting
  ! Date  : 09/15/18

  use Realization_Subsurface_class
  use Well_Type_class
  use Well_Data_class

  implicit none

  PetscInt, intent(in   ) :: fid
  type(realization_subsurface_type) :: realization
  PetscReal :: tconv, sign
  PetscInt :: towg_miscibility_model
  type(well_data_list_type), pointer :: well_data_list

  PetscInt :: iwell, nwell, welltype

  PetscReal :: wompriu, wgmpriu, wwmpriu, wsmpriu, &
               wompr, wompt, womir, womit, &
               wgmpr, wgmpt, wgmir, wgmit, &
               wwmpr, wwmpt, wwmir, wwmit, &
               wsmpr, wsmpt, wsmir, wsmit

  PetscReal :: fompr, fompt, fomir, fomit, &
               fgmpr, fgmpt, fgmir, fgmit, &
               fwmpr, fwmpt, fwmir, fwmit, &
               fsmpr, fsmpt, fsmir, fsmit

  tconv = 3600.0*24.0

  !  Find well list and loop over wells

  well_data_list => realization%well_data
  nwell = getnwell(well_data_list)

  !  Loop over wells

  do iwell = 1, nwell

  !  Set up well flow sign (+ ve producer and -ve injectors)

    welltype = getWellTypeI(iwell, well_data_list)
    if (wellType == PROD_WELL_TYPE) then
      sign = 1.0
    else
      sign =-1.0
    endif

  !  Get the well mass values

    wompriu = 0.0
    wgmpriu = 0.0
    wwmpriu = 0.0
    wsmpriu = 0.0

    wompr   = 0.0
    wgmpr   = 0.0
    wwmpr   = 0.0
    wsmpr   = 0.0

    womir   = 0.0
    wgmir   = 0.0
    wwmir   = 0.0
    wsmir   = 0.0

    wompriu = getWellTTValI(iwell, W_TARG_OM, VALTYPE_ACTUAL, well_data_list)
    wompt   = getWellTTValI(iwell, W_TARG_OM, VALTYPE_TOTALP, well_data_list)
    womit   = getWellTTValI(iwell, W_TARG_OM, VALTYPE_TOTALI, well_data_list)

    wgmpriu = getWellTTValI(iwell, W_TARG_GM, VALTYPE_ACTUAL, well_data_list)
    wgmpt   = getWellTTValI(iwell, W_TARG_GM, VALTYPE_TOTALP, well_data_list)
    wgmit   = getWellTTValI(iwell, W_TARG_GM, VALTYPE_TOTALI, well_data_list)

    wwmpriu = getWellTTValI(iwell, W_TARG_WM, VALTYPE_ACTUAL, well_data_list)
    wwmpt   = getWellTTValI(iwell, W_TARG_WM, VALTYPE_TOTALP, well_data_list)
    wwmit   = getWellTTValI(iwell, W_TARG_WM, VALTYPE_TOTALI, well_data_list)

    if (towg_miscibility_model == TOWG_SOLVENT_TL) then
      wsmpriu = getWellTTValI(iwell, W_TARG_SM, VALTYPE_ACTUAL, well_data_list)
      wsmpt   = getWellTTValI(iwell, W_TARG_SM, VALTYPE_TOTALP, well_data_list)
      wsmit   = getWellTTValI(iwell, W_TARG_SM, VALTYPE_TOTALI, well_data_list)
    endif

    if (wellType == PROD_WELL_TYPE) then
      wompr = wompriu * tconv * sign
      wgmpr = wgmpriu * tconv * sign
      wwmpr = wwmpriu * tconv * sign
      wsmpr = wsmpriu * tconv * sign
    else
      womir = wompriu * tconv * sign
      wgmir = wgmpriu * tconv * sign
      wwmir = wwmpriu * tconv * sign
      wsmir = wsmpriu * tconv * sign
    endif

  !  Write out well mass values

    call wrtToTable(fid, wompr)
    call wrtToTable(fid, wompt)
    call wrtToTable(fid, womir)
    call wrtToTable(fid, womit)

    call wrtToTable(fid, wgmpr)
    call wrtToTable(fid, wgmpt)
    call wrtToTable(fid, wgmir)
    call wrtToTable(fid, wgmit)

    call wrtToTable(fid, wwmpr)
    call wrtToTable(fid, wwmpt)
    call wrtToTable(fid, wwmir)
    call wrtToTable(fid, wwmit)

    if (towg_miscibility_model == TOWG_SOLVENT_TL) then
      call wrtToTable(fid, wsmpr)
      call wrtToTable(fid, wsmpt)
      call wrtToTable(fid, wsmir)
      call wrtToTable(fid, wsmit)
    endif

  enddo

  !  Now the field mass values

  fompr = GetFieldTTVal(W_TARG_OM, VALTYPE_ACTUALP, well_data_list)
  fomir = GetFieldTTVal(W_TARG_OM, VALTYPE_ACTUALI, well_data_list)
  fompt = GetFieldTTVal(W_TARG_OM, VALTYPE_TOTALP , well_data_list)
  fomit = GetFieldTTVal(W_TARG_OM, VALTYPE_TOTALI , well_data_list)

  fgmpr = GetFieldTTVal(W_TARG_GM, VALTYPE_ACTUALP, well_data_list)
  fgmir = GetFieldTTVal(W_TARG_GM, VALTYPE_ACTUALI, well_data_list)
  fgmpt = GetFieldTTVal(W_TARG_GM, VALTYPE_TOTALP , well_data_list)
  fgmit = GetFieldTTVal(W_TARG_GM, VALTYPE_TOTALI , well_data_list)

  fwmpr = GetFieldTTVal(W_TARG_WM, VALTYPE_ACTUALP, well_data_list)
  fwmir = GetFieldTTVal(W_TARG_WM, VALTYPE_ACTUALI, well_data_list)
  fwmpt = GetFieldTTVal(W_TARG_WM, VALTYPE_TOTALP , well_data_list)
  fwmit = GetFieldTTVal(W_TARG_WM, VALTYPE_TOTALI , well_data_list)

  if (towg_miscibility_model == TOWG_SOLVENT_TL) then
    fsmpr = GetFieldTTVal(W_TARG_SM, VALTYPE_ACTUALP, well_data_list)
    fsmir = GetFieldTTVal(W_TARG_SM, VALTYPE_ACTUALI, well_data_list)
    fsmpt = GetFieldTTVal(W_TARG_SM, VALTYPE_TOTALP , well_data_list)
    fsmit = GetFieldTTVal(W_TARG_SM, VALTYPE_TOTALI , well_data_list)
  else
    fsmpr = 0.0
    fsmir = 0.0
    fsmpt = 0.0
    fsmit = 0.0
  endif

!  Convert to rates per day

  fompr = fompr*tconv
  fomir = fomir*tconv

  fgmpr = fgmpr*tconv
  fgmir = fgmir*tconv

  fwmpr = fwmpr*tconv
  fwmir = fwmir*tconv

  fsmpr = fsmpr*tconv
  fsmir = fsmir*tconv

  !  Write out field mass values

   call wrtToTable(fid, fompr)
   call wrtToTable(fid, fompt)
   call wrtToTable(fid, fomir)
   call wrtToTable(fid, fomit)

   call wrtToTable(fid, fgmpr)
   call wrtToTable(fid, fgmpt)
   call wrtToTable(fid, fgmir)
   call wrtToTable(fid, fgmit)

   call wrtToTable(fid, fwmpr)
   call wrtToTable(fid, fwmpt)
   call wrtToTable(fid, fwmir)
   call wrtToTable(fid, fwmit)

   if (towg_miscibility_model == TOWG_SOLVENT_TL) then
     call wrtToTable(fid, fsmpr)
     call wrtToTable(fid, fsmpt)
     call wrtToTable(fid, fsmir)
     call wrtToTable(fid, fsmit)
   endif

end subroutine WriteWellMassValues

! *************************************************************************** !

subroutine WriteRestValues(realization, option)
  !
  ! Used to write restart file values
  !
  ! Author: Dave Ponting
  ! Date  : 12/15/18

  use Realization_Subsurface_class
  use Well_Data_class
  use Output_Eclipse_module, only:WriteEclipseFilesRest, GetMlmax
  use Option_module
  use Grid_module
  use Patch_module

  implicit none

  type(realization_subsurface_type) :: realization
  type(option_type), intent(in), pointer :: option
  type(well_data_list_type), pointer :: well_data_list
  PetscInt  :: ierr, nlmax, mlmax
  PetscBool :: is_ioproc
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch

  PetscReal, pointer        :: vsoll(:,:)
  PetscInt                  :: nsol
  character(len=8), pointer :: zsol(:)

  character(len=8), pointer :: wname(:)
  PetscInt, pointer :: wtype(:), wncmpl(:), &
                       ixcmpl(:), iycmpl(:), izcmpl(:), idcmpl(:)

  PetscReal :: time

  time = option%time

  ierr = 0

  is_ioproc = PETSC_FALSE
  if (option%myrank == option%io_rank) then
    is_ioproc = PETSC_TRUE
  endif

  grid => realization%patch%grid
  patch => realization%patch

  nlmax = grid%nlmax
  mlmax = GetMlmax()

  call allocateLocalSolution(vsoll, nsol, zsol, nlmax)
  call loadLocalSolution    (vsoll, nsol, zsol, patch, grid, option)
  well_data_list => realization%well_data
  call setupWellData(wname, wtype, wncmpl, ixcmpl, iycmpl, izcmpl, idcmpl, &
                     well_data_list)
  call WriteEclipseFilesRest(vsoll, nsol, zsol, time, is_ioproc, &
                             wname, wtype, wncmpl, ixcmpl, iycmpl, izcmpl, &
                             idcmpl, option)
  call deleteLocalSolution  (vsoll, zsol)
  call deleteWellData(wname, wtype, wncmpl, ixcmpl, iycmpl, izcmpl, idcmpl)

end subroutine WriteRestValues

! *************************************************************************** !

subroutine OutputLineRept(realization_base, option)
  !
  ! Write out single line per step progress reports in TOI_IMS and TOWG mode
  !
  ! Author: Dave Ponting
  ! Date  : 01/15/19

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Realization_Subsurface_class

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option

  if (     option%iflowmode == TOIL_IMS_MODE &
      .or. option%iflowmode == TOWG_MODE       ) then
  if (option%myrank == option%io_rank) then
      select type(realization_base)
       class is(realization_subsurface_type)
        call WriteLineRept(realization_base, option)
      end select
    endif
  endif

end subroutine OutputLineRept

! ************************************************************************** !

subroutine WriteLineRept(realization, option)
  !
  ! Used to write out single line progress reports
  !
  ! Author: Dave Ponting
  ! Date  : 01/15/19

  use Realization_Subsurface_class
  use Well_Type_class
  use Well_Data_class
  use Option_module

  implicit none

  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option

  PetscReal :: time, tconv, dt
  type(well_data_list_type), pointer :: well_data_list

  PetscInt :: icountp

  PetscReal :: fopt, fopr, &
               fgpr, fgir, &
               fwpr, fwir, &
               flpr, fgor, fwct, fpav
  character(len=8) :: stime, sdt, sopt, sopr, swpr, sgpr, swir, sgir, sgor, spav

100 format('Step   time     tstep    fopt     fopr     fwpr     ', &
           'fgpr     fwir     fgir     fwct    fgor     fpav     NL LI Ch')
101 format('       days     days     ksm3     sm3/d    ksm3/d   ', &
           'ksm3/d   ksm3/d   ksm3/d           ksm3/sm3 Bar')
102 format('------ -------- -------- -------- -------- -------- ', &
           '-------- -------- -------- ------- -------- -------- -- -- --')
103 format(I6, 1X, A8, 1X, A8, 1X, A8, 1X, A8, 1X, A8, 1X, A8, &
           1X, A8, 1X, A8, 1X, F7.5, 1X, A8, 1X, A8, &
           1X, I2, 1X, I2, 1X, I2)

  well_data_list => realization%well_data
  option => realization%option

  tconv = 3600.0*24.0

  time = option%time/tconv ! Convert to days
  dt   = option%dt/tconv
  if (time>0.0) then

  !  Get the field values

    fopt = GetFieldTTVal(W_TARG_OSV, VALTYPE_TOTALP , well_data_list)
    fopr = GetFieldTTVal(W_TARG_OSV, VALTYPE_ACTUALP, well_data_list)
    fgpr = GetFieldTTVal(W_TARG_GSV, VALTYPE_ACTUALP, well_data_list)
    fgir = GetFieldTTVal(W_TARG_GSV, VALTYPE_ACTUALI, well_data_list)
    fwpr = GetFieldTTVal(W_TARG_WSV, VALTYPE_ACTUALP, well_data_list)
    fwir = GetFieldTTVal(W_TARG_WSV, VALTYPE_ACTUALI, well_data_list)

    fopt = fopt*0.001
    fopr = fopr*tconv

    fgpr = fgpr*tconv*0.001
    fgir = fgir*tconv*0.001

    fwpr = fwpr*tconv
    fwir = fwir*tconv

    flpr = fopr+fwpr

    fwct = 0.0
    fgor = 0.0
    if (flpr > 0.0) fwct = fwpr/flpr
    if (fopr > 0.0) fgor = fgpr/fopr

  !  Convert field pressure to Bars

    call GetFieldData(fpav)
    fpav = fpav*1.0D-5

    if (mod(linerept_count, 20) == 0) then
      write(*, 100)
      write(*, 101)
      write(*, 102)
    endif

!   Increment linecount

    icountp = linerept_count + 1

!  Format the time field (prefer F8.1 if possible)

    if      (time<999999.0) then
      write(stime,'(F8.1)') time    ! Write as XXXXXX.X days
    else if (time<99999999.0 ) then
      write(stime,'(I8)') int(time) ! Write as XXXXXXXX days
    else
      call PrintTidy8(time,stime)   ! Write as general real days
    endif

!  Format the steplength (prefer F8.3 if possible)

    if      (dt<9999.0) then
      write(sdt,'(F8.3)') dt    ! Write as XXXX.XXX days
    else
      call PrintTidy8(dt,sdt)   ! Write as general real days
    endif

!  Format other real values

    call PrintTidy8(fopt, sopt)
    call PrintTidy8(fopr, sopr)
    call PrintTidy8(fwpr, swpr)
    call PrintTidy8(fgpr, sgpr)
    call PrintTidy8(fwir, swir)
    call PrintTidy8(fgir, sgir)
    call PrintTidy8(fgor, sgor)
    call PrintTidy8(fpav, spav)

! Print line report

    write(*, 103) icountp, stime, sdt, &
                  sopt, sopr, swpr, sgpr, &
                  swir, sgir, fwct, sgor, spav, &
                  option%nnl, option%linpernl, option%nchperst

    linerept_count = icountp

  endif

end subroutine WriteLineRept

!*****************************************************************************!

subroutine setupWellData(wname, wtype, wncmpl, &
                         ixcmpl, iycmpl, izcmpl, idcmpl, &
                         well_data_list)
  !
  ! Setup structures holding well locations for Output_Eclipse_module
  !
  ! Author: Dave Ponting
  ! Date  : 01/15/19

  use Well_Data_class

  implicit none

  character(len=8), pointer :: wname(:)
  PetscInt, pointer :: wtype(:), wncmpl(:), &
                       ixcmpl(:), iycmpl(:), izcmpl(:), idcmpl(:)
  type(well_data_list_type), pointer :: well_data_list

  PetscInt :: iw, nw, mw, ic, ncg, nct, mct, ci, cj, ck, cdd, welltype
  character(len=MAXSTRINGLENGTH) :: name

  nw = getnwell(well_data_list)
  mw = max(1, nw)

  allocate(wname (mw))
  allocate(wtype (mw))
  allocate(wncmpl(mw))

  nct = 0

  do iw=1, nw
    call getWellNameI(iw, well_data_list, name)
    welltype = getWellTypeI(iw, well_data_list)
    ncg = GetWellNCmplGI(iw, well_data_list)
    nct = nct+ncg
    wname (iw) = name(1:8)
    wtype (iw) = welltype
    wncmpl(iw) = ncg
  enddo

  mct = max(1, nct)

  allocate(ixcmpl(mct))
  allocate(iycmpl(mct))
  allocate(izcmpl(mct))
  allocate(idcmpl(mct))

  nct = 0

  do iw = 1, nw

    ncg = wncmpl(iw)

    do ic = 1, ncg
      call GetCmplGlobalLocI(iw, ic, ci, cj, ck, cdd, well_data_list)
      ixcmpl(nct+ic) = ci
      iycmpl(nct+ic) = cj
      izcmpl(nct+ic) = ck
      idcmpl(nct+ic) = cdd
    enddo

    nct = nct+ncg

  enddo

end subroutine setupWellData

!*****************************************************************************!

subroutine DeleteWellData(wname, wtype, wncmpl, ixcmpl, iycmpl, izcmpl, idcmpl)
  !
  ! Delete structures holding well locations for Output_Eclipse_module
  !
  ! Author: Dave Ponting
  ! Date  : 01/15/19

  implicit none

  character(len=8), pointer :: wname(:)
  PetscInt, pointer :: wtype(:), wncmpl(:), &
                       ixcmpl(:), iycmpl(:), izcmpl(:), idcmpl(:)

  deallocate(wname )
  deallocate(wtype )
  deallocate(wncmpl)

  deallocate(ixcmpl)
  deallocate(iycmpl)
  deallocate(izcmpl)
  deallocate(idcmpl)

end subroutine DeleteWellData

!*****************************************************************************!

subroutine allocateLocalSolution(vsoll, nsol, zsol, nlmax)
  !
  ! Allocate structures holding well locations for Output_Eclipse_module
  !
  ! Author: Dave Ponting
  ! Date  : 01/15/19

  implicit none

  PetscReal, pointer        :: vsoll(:,:)
  PetscInt                  :: nsol
  character(len=8), pointer :: zsol(:)
  PetscInt, intent(in)      :: nlmax

  nsol = 4
  allocate(vsoll(nlmax, nsol))
  allocate(zsol (       nsol))
  zsol(1) = 'PRESSURE'
  zsol(2) = 'SOIL'
  zsol(3) = 'SGAS'
  zsol(4) = 'SWAT'

end subroutine allocateLocalSolution

! ************************************************************************** !

subroutine loadLocalSolution(vsoll, nsol, zsol, patch, grid, option)
  !
  ! Load structures holding well locations for Output_Eclipse_module
  !
  ! Author: Dave Ponting
  ! Date  : 01/15/19

  use Patch_module
  use Grid_module
  use Option_module
  use PM_TOWG_Aux_module
  use PM_TOilIms_Aux_module

  implicit none

  PetscReal, pointer        :: vsoll(:,:)
  PetscInt                  :: nsol
  character(len=8), pointer :: zsol(:)

  type(patch_type), pointer :: patch
  type(grid_type ), pointer :: grid
  type(option_type) :: option

  PetscInt :: isol

  do isol = 1, nsol
    select case(option%iflowmode)
      case(TOIL_IMS_MODE)
        call patch%aux%TOil_ims%GetLocalSol(grid, patch%aux%material, &
                                            patch%imat, option, &
                                            vsoll, isol, zsol(isol))
      case(TOWG_MODE)
        call patch%aux%TOWG%GetLocalSol(grid, patch%aux%material, &
                                        patch%imat, option, &
                                        vsoll, isol, zsol(isol))
    end select
  enddo

end subroutine loadLocalSolution

! *************************************************************************** !

subroutine deleteLocalSolution(vsoll, zsol)
  !
  ! Delete structures holding well locations for Output_Eclipse_module
  !
  ! Author: Dave Ponting
  ! Date  : 01/15/19

  implicit none

  PetscReal, pointer        :: vsoll(:,:)
  character(len=8), pointer :: zsol(:)

  deallocate(vsoll)
  deallocate(zsol )

end subroutine deleteLocalSolution

! *************************************************************************** !

subroutine WrtHrd(fid, mnem, name, units, icolumn, zm, zn, zu, ni, mi, wecl)
  !
  ! Write out mnemonic, name and units values on stream fid
  ! and/or store in the zm/zn/zu buffers
  !
  ! Author: Dave Ponting
  ! Date  : 12/15/18

  use String_module

  implicit none

  PetscInt :: fid
  character(len=*), intent(in) :: mnem, name, units
  PetscInt, intent(in) :: icolumn
  character(len=8), allocatable, intent(inout) :: zm(:), zn(:), zu(:)
  PetscInt, intent(inout) :: ni
  PetscInt, intent(inout) :: mi
  PetscBool, intent(in) :: wecl

  character(len=8) :: mnemu, nameu, unitsu

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: cell

  !  Positive fid indicates -mas file output required

  if (fid > 0) then
    string = trim(name) // ' ' // trim(mnem)
    cell   = ' '
    call OutputWriteToHeader(fid, string, units, cell, icolumn)
  endif

  !  wecl indicates value storage for Eclipse output required

  if (wecl) then
    call checkHeaderBufferSize(zm, zn, zu, ni, mi)
    ni = ni+1
    mnemu  = mnem
    nameu  = name
    unitsu = units
    call StringToUpper(mnemu )
    call StringToUpper(nameu )
    call StringToUpper(unitsu)
    zm(ni) = mnemu
    zn(ni) = nameu
    zu(ni) = unitsu
  endif

end subroutine WrtHrd

! ************************************************************************** !

subroutine WrtHrdMO(fid, mnem, name, units, icolumn)
  !
  ! Write out mnemonic, name and units values on stream fid
  !
  ! Author: Dave Ponting
  ! Date  : 12/15/18

  use String_module

  implicit none

  PetscInt :: fid
  character(len=*), intent(in) :: mnem, name, units
  PetscInt, intent(in) :: icolumn

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: cell

  string = trim(name) // ' ' // trim(mnem)
  cell   = ' '
  call OutputWriteToHeader(fid, string, units, cell, icolumn)

end subroutine WrtHrdMO

! *************************************************************************** !

subroutine wrtToTableAndSumm(fid, val, vd, nd, md, wecl)
  !
  ! Write out values on stream fid and/or store in the vd value buffer
  !
  ! Author: Dave Ponting
  ! Date  : 12/15/18

  use Utility_module, only : ReallocateArray

  implicit none

110 format(es14.6)

  PetscInt  :: fid
  PetscReal :: val
  PetscReal, pointer :: vd(:)
  PetscInt, intent(inout) :: nd
  PetscInt, intent(in) :: md
  PetscBool, intent(in) :: wecl

  !  Positive fid indicates -mas file output required

  if (fid > 0) then
    write(fid, 110, advance="no") val
  endif

  !  wecl indicates value storage for Eclipse output required

  if (wecl) then
  !  Check if the buffer has space for another value: reallocate if not
    if (nd >= (md-1)) call reallocateArray(vd, md)
  !  Store value
    nd = nd+1
    vd(nd) = val
  endif

end subroutine wrtToTableAndSumm

! *************************************************************************** !

subroutine wrtToTable(fid, val)
  !
  ! Write out values on stream fid
  !
  ! Author: Dave Ponting
  ! Date  : 12/15/18

  use Utility_module, only : ReallocateArray

  implicit none

110 format(es14.6)

  PetscInt  :: fid
  PetscReal :: val

  write(fid, 110, advance="no") val

end subroutine wrtToTable

! ************************************************************************** !

subroutine checkHeaderBufferSize(zm, zn, zu, ni, mi)
  !
  ! Check header buffer is large enough and extend if required
  !
  ! Author: Dave Ponting
  ! Date  : 01/15/19

  implicit none

  character(len=8), allocatable, intent(inout) :: zm(:), zn(:), zu(:)
  PetscInt, intent(inout) :: ni
  PetscInt, intent(inout) :: mi

  character(len=8), allocatable :: zmt(:), znt(:), zut(:)
  PetscInt :: i

  ! Check if buffer as space for anaother set of values; reallocate if not

  if (ni >= (mi-1)) then

  ! Allocate temporary stores

   allocate(zmt(ni))
   allocate(znt(ni))
   allocate(zut(ni))

  ! Copy to temporary stores

   do i = 1, ni
     zmt(i) = zm(i)
     znt(i) = zn(i)
     zut(i) = zu(i)
   enddo

  !  Deallocate, extend and reallocate actual stores

   deallocate(zm)
   deallocate(zn)
   deallocate(zu)

   mi = 2*mi

   allocate(zm(mi))
   allocate(zn(mi))
   allocate(zu(mi))

   zm = ' '
   zn = ' '
   zu = ' '

  !  Copy values back

   do i = 1, ni
     zm(i) = zmt(i)
     zn(i) = znt(i)
     zu(i) = zut(i)
   enddo

  !  Deallocate the temporary stores

   deallocate(zmt)
   deallocate(znt)
   deallocate(zut)

  endif

end subroutine checkHeaderBufferSize

! ************************************************************************** !

subroutine setupEwriterRestMaps(patch, grid, option)
  !
  ! Set up maps required to convert from Pflotran order to that
  ! required for Eclipse restart files
  !
  ! Author: Dave Ponting
  ! Date  : 01/15/19

  use Patch_module
  use Grid_module
  use Option_module
  use Output_Eclipse_module, only: setupRestMaps

  implicit none

  type(patch_type), pointer :: patch
  type(grid_type ), pointer :: grid
  type(option_type) :: option

  PetscInt :: nlmax, lid, gid, nid, ierr

  PetscInt :: mlmax

  PetscInt, allocatable :: ltoa(:)

  !  First, find the maximum value of nlmax over all procs

  nlmax = grid%nlmax
  ierr  = 0
  call MPI_AllReduce(nlmax, mlmax, ONE_INTEGER_MPI, MPI_INTEGER, MPI_MAX, &
                     option%mycomm, ierr)

  !  Allocate local to active map and fill on this proc

  allocate(ltoa(mlmax))
  do lid = 1, nlmax
    gid        = grid%nL2G(lid  )
    nid = grid%nG2A(gid)
    if (patch%imat(gid) <= 0) cycle
      ltoa(lid) = nid
  enddo

  !  Setup the restart maps in Output_Eclipse_module

  call setupRestMaps(ltoa, option, nlmax, mlmax)

  !  Delete the ltoa work array

  deallocate(ltoa)

end subroutine setupEwriterRestMaps

! ************************************************************************** !

function GetEclWrtFlg(count, time, deltat, deltas, lastt, lasts)
  !
  ! Check if Eclipse output required for this step and time
  !
  ! Author: Dave Ponting
  ! Date  : 01/15/19

  implicit none

  PetscBool :: GetEclWrtFlg

  PetscInt , intent(in)    :: count
  PetscReal, intent(in)    :: time
  PetscReal, intent(in)    :: deltat
  PetscInt , intent(in)    :: deltas
  PetscReal, intent(inout) :: lastt
  PetscInt , intent(inout) :: lasts

  GetEclWrtFlg = PETSC_TRUE

  if (count == 0) then

  !  First call: will write, so set the last-write values to now

    lastt  = time
    lasts  = 0
  else

  !  Not first call: assume not writing unless a case qualifies

    GetEclWrtFlg = PETSC_FALSE

  ! delta-time has been set

    if (deltat > 0.0) then
  !  If deltat has elapsed since last write, write and reset last time
      if ((time-lastt) >= deltat) then
        GetEclWrtFlg = PETSC_TRUE
        if (abs(mod(time, deltat)) == 0.0) then
          lastt = time
        else
          lastt = deltat*int(time/deltat)
        endif
      endif
    endif

  !  delta-step has been set

    if (deltas > 0) then
  !  If deltas steps since last write, write and reset last step
      if ((count-lasts) >= deltas) then
        GetEclWrtFlg = PETSC_TRUE
        lasts = count
      endif
    endif
  endif

end function GetEclWrtFlg

! ************************************************************************** !

subroutine PrintTidy8(va, s)
  !
  ! Write out real number in neat 8-character format
  !
  ! Author: Dave Ponting
  ! Date  : 03/08/19

  implicit none
  PetscReal, intent(in) :: va
  character(len=8), intent(out) :: s
  character(len=8) :: zn
  PetscBool :: isNeg
  PetscInt  :: exponent, ilog10
  PetscReal :: v, vl, vu

  v  = va
  vu = 9999.99
  vl = 1.0e-3
  isNeg    = PETSC_FALSE
  exponent = 0

  s = '     0.0' ! Default

  if (v<0.0) then
    isNeg = PETSC_TRUE;v=-v
  endif

  if (v>1.0e-10) then

    iLog10 = int(log10(1.0000001*v))
    if ((v >= vu) .or. (v < vl)) then
      exponent = iLog10
      v = v/10.0**iLog10
      if (isNeg) then
        write(zn, '(F4.2)') v ! Write value as XX.X
      else
        write(zn, '(F5.3)') v ! Write value as XX.XX
      endif
      s = trim(adjustl(zn))//'E' ! Set to 4 or 5 chars+exp => 5 or 6 chars
      write(zn, '(I2)') exponent
      s = trim(adjustl(s))//trim(adjustl(zn)) ! Add two chars exponent
    else
      if (isNeg) then
        write(s, '(F7.2)') v ! Allow space for -9999.99
      else
        write(s, '(F8.3)') v ! Use all the chars up to 9999.999
      endif
    endif

    if (isNeg) then
     s = '-' // trim(adjustl(s))
    endif

    s = trim(adjustr(s))

  else
    s = '     0.0'
  endif

end subroutine PrintTidy8

end module Output_Observation_module
