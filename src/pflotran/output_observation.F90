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
  PetscBool :: check_for_obs_points
  PetscBool :: calculate_velocities ! true if any obs. pt. prints velocity
  PetscBool :: mass_balance_first
  PetscBool :: integral_flux_first

  public :: OutputObservation, &
            OutputObservationInit, &
            OutputMassBalance, &
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
  if (num_steps == 0) then
    observation_first = PETSC_TRUE
    mass_balance_first = PETSC_TRUE
    integral_flux_first = PETSC_TRUE
  else
    observation_first = PETSC_FALSE
    mass_balance_first = PETSC_FALSE
    integral_flux_first = PETSC_FALSE
  endif

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

!  if (realization_base%output_option%print_hdf5) then
!    call OutputObservationHDF5(realization)
!    call OutputObservationTecplot(realization)
!  endif
 
!  if (realization_base%output_option%print_tecplot .or. &
!      realization_base%output_option%print_hdf5) then
  if (realization_base%output_option%print_observation) then
    call OutputObservationTecplotColumnTXT(realization_base)
    call OutputIntegralFlux(realization_base)
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
 
  implicit none

  class(realization_base_type) :: realization_base
  
  PetscInt :: fid, icell
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string!, string2
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  type(observation_type), pointer :: observation
  PetscBool, save :: open_file = PETSC_FALSE
  PetscReal, allocatable :: velocities(:,:,:)
  PetscInt :: local_id
  PetscInt :: icolumn
  PetscInt :: nphase
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
    nphase = option%nphase
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
               '-obs-' // trim(adjustl(string)) // '.tec'
  
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
 !             call printErrMsg(option)
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

    if (option%nphase > 1) then
      call OutputWriteToHeader(fid,'qgx',string,cell_string,icolumn)
      call OutputWriteToHeader(fid,'qgy',string,cell_string,icolumn)
      call OutputWriteToHeader(fid,'qgz',string,cell_string,icolumn)
    endif
  endif
    
end subroutine WriteObservationHeader

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

  implicit none
  
  PetscInt :: fid
  class(realization_base_type) :: realization_base
  character(len=MAXWORDLENGTH) :: coupler_name
  
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  
  option => realization_base%option
  
  select case(option%iflowmode)
    case(TH_MODE)
      string = ',"Darcy flux ' // trim(coupler_name) // &
               ' [m^3/' // trim(realization_base%output_option%tunit) // ']"'
    case default
  end select
  write(fid,'(a)',advance="no") trim(string)


end subroutine WriteObservationHeaderForBC

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
  
  PetscInt :: fid
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
  type(option_type), pointer :: option
  PetscErrorCode :: ierr
  
  option => realization_base%option

110 format(es14.6)
 
  iphase = 1

  ! sum up fluxes across region
  if (associated(connection_set)) then
    offset = connection_set%offset
    select case(option%iflowmode)
      case(TH_MODE)
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

  endif

end subroutine WriteObservationDataForBC

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

  if (option%nphase > 1) then
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
  PetscInt :: iphase

  PetscReal :: velocity(1:3)

  option => realization_base%option
  
200 format(3(es14.6))

  iphase = 1
  velocity = GetVelocityAtCoord(fid,realization_base,region%cell_ids(1), &
                                region%coordinates(ONE_INTEGER)%x, &
                                region%coordinates(ONE_INTEGER)%y, &
                                region%coordinates(ONE_INTEGER)%z,iphase)
  write(fid,200,advance="no") velocity(1:3)*realization_base%output_option%tconv   

  if (option%nphase > 1) then
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
  PetscReal :: weight, distance
  PetscReal :: sum_velocity(1:3), velocity(1:3)
  PetscReal :: sum_weight(1:3)
  
  option => realization_base%option
  patch => realization_base%patch
  grid => patch%grid
  field => realization_base%field

  sum_velocity = 0.d0
  sum_weight = 0.d0

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
  use Integral_Flux_module
  use Utility_module

  implicit none

  class(realization_base_type), target :: realization_base

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(output_option_type), pointer :: output_option


  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: units
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


  if (.not.associated(patch%integral_flux_list%first)) return

  flow_dof_scale = 1.d0
  select case(option%iflowmode)
    case(TH_MODE)
      flow_dof_scale(1) = FMWH2O
  end select

  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '-int.dat'
  else
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-int.dat'
  endif
  
  ! open file
  if (option%myrank == option%io_rank) then

!geh    option%io_buffer = '--> write tecplot mass balance file: ' // trim(filename)
!geh    call printMsg(option)    

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
      

      integral_flux => patch%integral_flux_list%first
      do
        if (.not.associated(integral_flux)) exit
        select case(option%iflowmode)
          case(TH_MODE)
            string = trim(integral_flux%name) // ' Water'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            units = 'kg/' // trim(output_option%tunit) // ''
            string = trim(integral_flux%name) // ' Water'
            call OutputWriteToHeader(fid,string,units,'',icol)
        end select

        select case(option%iflowmode)
          case(TH_MODE)
            string = trim(integral_flux%name) // ' Energy'
            call OutputWriteToHeader(fid,string,'MJ','',icol)
            units = 'MJ/' // trim(output_option%tunit) // ''
            string = trim(integral_flux%name) // ' Energy'
            call OutputWriteToHeader(fid,string,units,'',icol)
        end select
        
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
  
  allocate(array(option%nflowdof,2))
  allocate(array_global(option%nflowdof,2))
  allocate(instantaneous_array(option%nflowdof))
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
  
  use Flowmode_module, only : FlowmodeComputeMassBalance

  use Global_Aux_module
  use Material_Aux_class
  use String_module

  implicit none

  class(realization_base_type), target :: realization_base

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(output_option_type), pointer :: output_option
  type(coupler_type), pointer :: coupler
  type(mass_balance_region_type), pointer :: cur_mbr
  type(global_auxvar_type), pointer :: global_auxvars_bc_or_ss(:)

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: units
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: fid = 86
  PetscInt :: icol
  PetscInt :: iconn
  PetscInt :: offset

  PetscInt :: iphase, ispec
  PetscReal :: mass_balance(realization_base%option%nflowspec,realization_base%option%nphase)
  PetscReal :: sum_local(realization_base%option%nphase)
  PetscReal :: sum_global(realization_base%option%nphase)
  
  PetscReal :: global_water_mass
  PetscReal :: sum_allwater  ! sum of global water mass and fluxes
  
  PetscMPIInt :: int_mpi
  PetscBool :: bcs_done
  PetscErrorCode :: ierr

  !-----------------------------------------------------------------------------
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option

  output_option => realization_base%output_option

  if (option%nflowspec<=0) return

  do ispec = 1, option%nflowspec  ! not yet output those for energy (TODO)

    if (len_trim(output_option%plot_name) > 2) then
      filename = trim(output_option%plot_name) // &
               '_' // trim(StringFormatInt(ispec)) // &
               '-mas.dat'
    else
      filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '_' // trim(StringFormatInt(ispec)) // &
               '-mas.dat'
    endif

    !-------------------------
    ! open file and write headers if first time
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
      
        ! total mass
        select case(option%iflowmode)
          case(TH_MODE)
            call OutputWriteToHeader(fid,'Global: Water(L)', &
                                    '-kg','',icol)
            call OutputWriteToHeader(fid,'Global: Water(G)', &
                                    '-kg','',icol)
            call OutputWriteToHeader(fid,'Global: Water(S)', &
                                    '-kg','',icol)
          case default
            ! do nothing
        end select

        ! BC/SS
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
        
          string = trim(coupler%name) // ': Water(L)'
          call OutputWriteToHeader(fid,string,'-kg','',icol)
          string = trim(coupler%name) // ': Water(G)'
          call OutputWriteToHeader(fid,string,'-kg','',icol)
          string = trim(coupler%name) // ': Water(S)'
          call OutputWriteToHeader(fid,string,'-kg','',icol)

          units = 'kg/' // trim(output_option%tunit) // ''
          string = trim(coupler%name) // ': Water(L)'
          call OutputWriteToHeader(fid,string,units,'',icol)
          string = trim(coupler%name) // ': Water(G)'
          call OutputWriteToHeader(fid,string,units,'',icol)
          string = trim(coupler%name) // ': Water(S)'
          call OutputWriteToHeader(fid,string,units,'',icol)

          coupler => coupler%next
        enddo
        ! sum of total mass and bc/ss
        if (option%nflowdof > 0) then
          string = ' SUM of All Water Mass and Flux'
          call OutputWriteToHeader(fid,string,'kg','',icol)
        endif

        ! Print the water mass [kg] and species mass [mol] in the specified regions (header)
        if (associated(output_option%mass_balance_region_list)) then
          cur_mbr => output_option%mass_balance_region_list
          do
            if (.not.associated(cur_mbr)) exit
            string = 'Region ' // trim(cur_mbr%region_name) // ' Water(L)'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            string = 'Region ' // trim(cur_mbr%region_name) // ' Water(G)'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            string = 'Region ' // trim(cur_mbr%region_name) // ' Water(S)'
            call OutputWriteToHeader(fid,string,'kg','',icol)


            cur_mbr => cur_mbr%next
          enddo
        endif
      
        write(fid,'(a)') ''
      
      else
        open(unit=fid,file=filename,action="write",status="old",position="append")
      endif
    
    endif

100 format(100es16.8)
110 format(100es16.8)

    ! ----------------------------------
    ! write time
    if (option%myrank == option%io_rank) then
      write(fid,100,advance="no") option%time/output_option%tconv
    endif
  
    if (option%nflowdof > 0) then
      if (option%myrank == option%io_rank) &
        write(fid,100,advance="no") option%flow_dt/output_option%tconv
    endif
  
    ! ----------------------------------
    ! print out mass balance data

    if (option%nflowdof > 0) then
      select type(realization_base)
        class is(realization_subsurface_type)
          select case(option%iflowmode)
            case(TH_MODE)
              call FlowmodeComputeMassBalance(realization_base,mass_balance)
          end select
        class default
          option%io_buffer = 'Unrecognized realization class in MassBalance().'
          call printErrMsg(option)
      end select

      sum_local = mass_balance(ispec,:)
      sum_allwater = 0.d0

      int_mpi = option%nflowspec*option%nphase
      call MPI_Reduce(sum_local,sum_global, &
                    int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)

      if (option%myrank == option%io_rank) then
        do iphase = 1, option%nphase
          write(fid,110,advance="no") sum_global(iphase)
            !
          sum_allwater = sum_allwater + sum_global(iphase)
        enddo
      endif
    endif
  
    coupler => patch%boundary_condition_list%first
    global_auxvars_bc_or_ss => patch%aux%Global%auxvars_bc
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
          else
            exit
          endif
        endif
      endif

      offset = coupler%connection_set%offset
      if (option%nflowdof > 0 &
         .and. option%compute_mass_balance_new) then
        ! print out cumulative H2O flux
        sum_local = 0.d0
        do iconn = 1, coupler%connection_set%num_connections
          sum_local = sum_local + global_auxvars_bc_or_ss(offset+iconn)%mass_balance(ispec,:)
        enddo
        int_mpi = option%nphase
        call MPI_Reduce(sum_local,sum_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
        if (option%myrank == option%io_rank) then
          ! change sign for positive in / negative out
          write(fid,110,advance="no") -sum_global
          !
          do iphase = 1, option%nphase
              sum_allwater = sum_allwater+sum_global(iphase)
          enddo
        endif

        ! print out H2O flux (mass)
        sum_local = 0.d0
        do iconn = 1, coupler%connection_set%num_connections
          sum_local = sum_local + global_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta(ispec,:)
        enddo
        ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
        sum_local = sum_local*FMWH2O

        int_mpi = option%nphase
        call MPI_Reduce(sum_local,sum_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
        if (option%myrank == option%io_rank) then
          ! change sign for positive in / negative out
          write(fid,110,advance="no") -sum_global*output_option%tconv   ! ? need checking here
        endif

      endif

      coupler => coupler%next
    enddo
    ! sum of all mass-balance items
    ! print out sum of global water mass and bc/ss fluxes
    ! (NOTE: THIS amount MUST BE constant, otherwise error in mass-conservation)
    if (option%myrank == option%io_rank) then
      write(fid,110,advance="no") sum_allwater
    endif
  
    ! Print the total water and component mass in the specified regions (data)
    if (associated(output_option%mass_balance_region_list)) then
      cur_mbr => output_option%mass_balance_region_list
      do
        if (.not.associated(cur_mbr)) exit
        call PatchGetWaterMassInRegion(cur_mbr%region_cell_ids, &
                                     cur_mbr%num_cells,patch,option, &
                                     global_water_mass)
        write(fid,110,advance="no") global_water_mass
        cur_mbr => cur_mbr%next
      enddo
    endif
  
    if (option%myrank == option%io_rank) then
      write(fid,'(a)') ''
      close(fid)
    endif
  
  end do  ! ispec =1, nflowspec-1

  mass_balance_first = PETSC_FALSE

end subroutine OutputMassBalance

end module Output_Observation_module
