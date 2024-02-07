module Output_HDF5_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Logging_module
  use Output_Aux_module
  use Output_Common_module

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal

  implicit none

  private

  PetscMPIInt, private, parameter :: ON=1, OFF=0

  ! flags signifying the first time a routine is called during a given
  ! simulation
  PetscBool :: hdf5_first

  public :: OutputHDF5Init, &
            OutputHDF5, &
            OutputHDF5UGridXDMF, &
            OutputHDF5FilenameID, &
            OutputHDF5UGridXDMFExplicit, &
            OutputHDF5DatasetStringArray, &
            OutputHDF5AttributeStringArray, &
            OutputHDF5OpenFile, &
            OutputHDF5CloseFile

  public :: OutputH5OpenFile, &
            OutputH5CloseFile, &
            OutputH5OpenGroup, &
            OutputH5CloseGroup, &
            OutputXMFOpenFile, &
            DetermineNumVertices, &
            OutputHDF5WriteStructCoordGroup, &
            WriteHDF5CoordinatesUGridXDMF

contains

! ************************************************************************** !

subroutine OutputHDF5Init(num_steps)
  !
  ! Initializes module variables for HDF5 output
  !
  ! Author: Glenn Hammond
  ! Date: 01/16/13
  !

  use Option_module

  implicit none

  PetscInt :: num_steps

  if (num_steps == 0) then
    hdf5_first = PETSC_TRUE
  else
    hdf5_first = PETSC_FALSE
  endif

end subroutine OutputHDF5Init

! ************************************************************************** !

subroutine OutputHDF5(realization_base,var_list_type)
  !
  ! Print to HDF5 file
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  use String_module

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use HDF5_module
  use HDF5_Aux_module

  implicit none

  class(realization_base_type) :: realization_base
  PetscInt :: var_list_type

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable

  Vec :: global_vec
  Vec :: global_vec_vx
  Vec :: global_vec_vy
  Vec :: global_vec_vz

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscBool :: first
  PetscInt :: ivar
  PetscBool :: include_gas_phase
  PetscErrorCode :: ierr

  discretization => realization_base%discretization
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option

  call OutputHDF5OpenFile(option, output_option, var_list_type, file_id, first)

  grid => patch%grid
  if (first) then
    call OutputHDF5Provenance(option, output_option, file_id)
    call OutputHDF5WriteStructCoordGroup(file_id,discretization,grid,option)
  endif

  ! create a group for the data set
  write(string,'(''Time:'',es13.5,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  if (len_trim(output_option%plot_name) > 2) then
    string = trim(string) // ' ' // output_option%plot_name
  endif
  !string = trim(string3) // ' ' // trim(string)
  call HDF5GroupOpenOrCreate(file_id,string,grp_id,option)

  ! write group attributes
  call OutputHDF5WriteSnapShotAtts(grp_id,option)

  ! write out data sets
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vx)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vy)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vz)

  select case (var_list_type)

    case (INSTANTANEOUS_VARS)
      ! loop over snapshot variables and write to file
      cur_variable => output_option%output_snap_variable_list%first
      do
        if (.not.associated(cur_variable)) exit
        call OutputGetVariableArray(realization_base,global_vec,cur_variable)
        string = OutputVariableGetName(cur_variable)
        call StringSwapChar(string," ","_")
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif
        if (cur_variable%iformat == 0) then
          call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                             global_vec,grp_id, &
                                             H5T_NATIVE_DOUBLE)
        else
          call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                             global_vec,grp_id, &
                                             H5T_NATIVE_INTEGER)
        endif
        cur_variable => cur_variable%next
      enddo

    case (AVERAGED_VARS)
      cur_variable => output_option%aveg_output_variable_list%first
      do ivar = 1,output_option%aveg_output_variable_list%nvars
        string = 'Aveg. ' // OutputVariableGetName(cur_variable)
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif

        call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                           field%avg_vars_vec(ivar),grp_id, &
                                           H5T_NATIVE_DOUBLE)

        cur_variable => cur_variable%next
      enddo

  end select

  include_gas_phase = PETSC_FALSE
  if (option%nphase > 1 .or. option%transport%nphase > 1) then
    include_gas_phase = PETSC_TRUE
  endif
  if (output_option%print_hdf5_vel_cent .and. &
      (var_list_type==INSTANTANEOUS_VARS)) then

    ! velocities
    call OutputGetCellCenteredVelocities(realization_base, global_vec_vx, &
                                         global_vec_vy,global_vec_vz, &
                                         LIQUID_PHASE)

    string = "Liquid X-Velocity [m_per_" // trim(output_option%tunit) // "]"
    call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                       global_vec_vx,grp_id,H5T_NATIVE_DOUBLE)

    string = "Liquid Y-Velocity [m_per_" // trim(output_option%tunit) // "]"
    call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                       global_vec_vy,grp_id,H5T_NATIVE_DOUBLE)

    string = "Liquid Z-Velocity [m_per_" // trim(output_option%tunit) // "]"
    call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                       global_vec_vz,grp_id,H5T_NATIVE_DOUBLE)

    if (include_gas_phase) then
        call OutputGetCellCenteredVelocities(realization_base,global_vec_vx, &
                                             global_vec_vy,global_vec_vz, &
                                             GAS_PHASE)
        string = "Gas X-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                           global_vec_vx,grp_id, &
                                           H5T_NATIVE_DOUBLE)

        string = "Gas Y-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                           global_vec_vy,grp_id, &
                                           H5T_NATIVE_DOUBLE)

        string = "Gas Z-Velocity"
        call HDF5WriteStructDataSetFromVec(string,realization_base, &
                                           global_vec_vz,grp_id, &
                                           H5T_NATIVE_DOUBLE)
    endif
  endif

  if (output_option%print_hdf5_vel_face .and. &
     (var_list_type==INSTANTANEOUS_VARS)) then

    ! internal flux velocities
    if (grid%structured_grid%nx > 1) then
        string = "Liquid X-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE, &
                                     X_DIRECTION,grp_id)
        if (include_gas_phase) then
          string = "Gas X-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE, &
                                       X_DIRECTION,grp_id)
        endif
    endif

    if (grid%structured_grid%ny > 1) then
        string = "Liquid Y-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE, &
                                     Y_DIRECTION,grp_id)
        if (include_gas_phase) then
          string = "Gas Y-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE, &
                                       Y_DIRECTION,grp_id)
        endif
    endif

    if (grid%structured_grid%nz > 1) then
        string = "Liquid Z-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization_base,LIQUID_PHASE, &
                                     Z_DIRECTION,grp_id)
        if (include_gas_phase) then
          string = "Gas Z-Flux Velocities"
          call WriteHDF5FluxVelocities(string,realization_base,GAS_PHASE, &
                                       Z_DIRECTION,grp_id)
        endif
    endif

  endif

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vx,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vy,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vz,ierr);CHKERRQ(ierr)

  call HDF5GroupClose(grp_id,option)

  call OutputHDF5CloseFile(option, file_id)

  hdf5_first = PETSC_FALSE

end subroutine OutputHDF5

! ************************************************************************** !

subroutine OutputHDF5OpenFile(option, output_option, var_list_type, file_id, &
                              first)
  !
  ! Determine the propper hdf5 output file name and open it.
  !
  ! Return the file handle and 'first' flag indicating if this is the
  ! first time the file has been opened.
  !
  use Option_module
  use hdf5
  use HDF5_Aux_module

  implicit none

  type(option_type), intent(inout) :: option
  type(output_option_type), intent(in) :: output_option
  PetscInt, intent(in) :: var_list_type
  character(len=MAXSTRINGLENGTH) :: filename
  PetscBool, intent(out) :: first

  integer(HID_T), intent(out) :: file_id

  character(len=MAXSTRINGLENGTH) :: string,string2,string3

  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      string2=''
      write(string3,'(i4)') output_option%plot_number
    case (AVERAGED_VARS)
      string2='-aveg'
      write(string3,'(i4)') &
        int(option%time/output_option%periodic_snap_output_time_incr)
  end select

  if (output_option%print_single_h5_file) then
    first = hdf5_first
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               trim(string2) // '.h5'
  else
    string = OutputHDF5FilenameID(output_option,option,var_list_type)
    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        if (mod(output_option%plot_number, &
                output_option%times_per_h5_file) == 0) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
      case (AVERAGED_VARS)
        if (Equal(mod((option%time-output_option%periodic_snap_output_time_incr)/ &
             output_option%periodic_snap_output_time_incr, &
             dble(output_option%times_per_h5_file)),0.d0)) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
    end select

    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                '-' // trim(string) // trim(string2) // '.h5'
  endif

  if (.not.first) then
    call HDF5FileTryOpen(filename,file_id,first,option%comm)
  endif
  if (first) then
    call HDF5FileOpen(filename,file_id,PETSC_TRUE,option)
  endif

  if (first) then
    option%io_buffer = ' --> creating hdf5 output file: ' // trim(filename)
  else
    option%io_buffer = ' --> appending to hdf5 output file: ' // trim(filename)
  endif
  call PrintMsg(option)

end subroutine OutputHDF5OpenFile

! ************************************************************************** !

subroutine OutputHDF5CloseFile(option, file_id)

  use Option_module
  use hdf5
  use HDF5_Aux_module

  implicit none

  type(option_type), intent(in) :: option
  integer(HID_T), intent(in) :: file_id

  call HDF5FileClose(file_id,option)

end subroutine OutputHDF5CloseFile

! ************************************************************************** !

subroutine OutputHDF5UGridXDMF(realization_base,var_list_type)
  !
  ! This routine writes unstructured grid data in HDF5 XDMF format.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 10/29/2012
  !

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  use String_module

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec
  use HDF5_Aux_module

  implicit none

  class(realization_base_type) :: realization_base
  PetscInt :: var_list_type

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable

  Vec :: global_vec
  Vec :: global_vec_vx,global_vec_vy,global_vec_vz
  Vec :: natural_vec

  character(len=MAXSTRINGLENGTH) :: filename_path, filename_header
  character(len=MAXSTRINGLENGTH) :: xmf_filename, att_datasetname, group_name
  character(len=MAXSTRINGLENGTH) :: string, string2,string3
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscBool :: first
  PetscInt :: ivar
  Vec :: ivec
  PetscErrorCode :: ierr

  discretization => realization_base%discretization
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option

  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      string2=''
      write(string3,'(i4)') output_option%plot_number
      xmf_filename = OutputFilename(output_option,option,'xmf','')
    case (AVERAGED_VARS)
      string2='-aveg'
      write(string3,'(i4)') &
        int(option%time/output_option%periodic_snap_output_time_incr)
      xmf_filename = OutputFilename(output_option,option,'xmf','aveg')
  end select

  if (output_option%print_single_h5_file) then
    first = hdf5_first
    !filename = trim(option%global_prefix) // trim(string2) // &
    !           trim(option%group_prefix) // '.h5'
    filename_path = trim(option%global_prefix) // trim(string2) // &
               trim(option%group_prefix) // '.h5'
    filename_header = trim(StringGetFilename(option%global_prefix)) // &
                      trim(string2) // trim(option%group_prefix) // '.h5'
  else
    string = OutputHDF5FilenameID(output_option,option,var_list_type)
    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        if (mod(output_option%plot_number, &
                output_option%times_per_h5_file) == 0) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
      case (AVERAGED_VARS)
        if (Equal(mod((option%time-output_option%periodic_snap_output_time_incr)/ &
             output_option%periodic_snap_output_time_incr, &
             dble(output_option%times_per_h5_file)),0.d0)) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
    end select

    !filename = trim(option%global_prefix) // trim(option%group_prefix) // &
    !           trim(string2) // '-' // trim(string) // '.h5'
    filename_path = trim(option%global_prefix) // &
                    trim(option%group_prefix) // &
                    trim(string2) // '-' // trim(string) // '.h5'
    filename_header = trim(StringGetFilename(option%global_prefix)) // &
                    trim(option%group_prefix) // &
                    trim(string2) // '-' // trim(string) // '.h5'
  endif

  grid => patch%grid

  if (.not.first) then
    call HDF5FileTryOpen(filename_path,file_id,first,option%comm)
  endif
  if (first) then
    call HDF5FileOpen(filename_path,file_id,PETSC_TRUE,option)
  else if (Uninitialized(realization_base%output_option%xmf_vert_len)) then
    call DetermineNumVertices(realization_base,option)
  endif

  if (first) then
    option%io_buffer = ' --> creating hdf5 output file: ' // trim(filename_path)
  else
    option%io_buffer = ' --> appending to hdf5 output file: ' // trim(filename_path)
  endif
  call PrintMsg(option)

  if (first) then
    ! create a group for the coordinates data set
    string = "Domain"
    call HDF5GroupCreate(file_id,string,grp_id,option)
    call WriteHDF5CoordinatesUGridXDMF(realization_base,option,grp_id)
    call HDF5GroupClose(grp_id,option)
  endif

  if (OptionIsIORank(option)) then
    option%io_buffer = ' --> write xmf output file: ' // trim(filename_path)
    call PrintMsg(option)
    open(unit=OUTPUT_UNIT,file=xmf_filename,action="write")
    call OutputXMFHeader(OUTPUT_UNIT, &
                         option%time/output_option%tconv, &
                         grid%nmax, &
                         realization_base%output_option%xmf_vert_len, &
                         grid%unstructured_grid%num_vertices_global,&
                         filename_header,PETSC_TRUE)
  endif

  ! create a group for the data set
  write(string,'(''Time'',es13.5,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  if (len_trim(output_option%plot_name) > 2) then
    string = trim(string) // ' ' // output_option%plot_name
  endif
  string = trim(string3) // ' ' // trim(string)

  call HDF5GroupOpenOrCreate(file_id,string,grp_id,option)
  group_name=string

  ! write group attributes
  call OutputHDF5WriteSnapShotAtts(grp_id,option)

  ! write out data sets
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vx)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vy)
  call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vz)

  select case (var_list_type)

    case (INSTANTANEOUS_VARS)
      ! loop over snapshot variables and write to file
      cur_variable => output_option%output_snap_variable_list%first
      do
        if (.not.associated(cur_variable)) exit
        call OutputGetVariableArray(realization_base,global_vec,cur_variable)
        call DiscretizationGlobalToNatural(discretization,global_vec, &
                                           natural_vec,ONEDOF)
        string = OutputVariableGetName(cur_variable)
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif
        if (cur_variable%iformat == 0) then
          call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                       H5T_NATIVE_DOUBLE)
        else
          call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                       H5T_NATIVE_INTEGER)
        endif
        att_datasetname = trim(filename_header) // ":/" // trim(group_name) // &
                          "/" // trim(string)
        if (OptionIsIORank(option)) then
          call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                  att_datasetname,CELL_CENTERED_OUTPUT_MESH)
        endif
        cur_variable => cur_variable%next
      enddo

    case (AVERAGED_VARS)
      if (associated(output_option%aveg_output_variable_list%first)) then
        cur_variable => output_option%aveg_output_variable_list%first
        do ivar = 1,output_option%aveg_output_variable_list%nvars
          string = 'Aveg. ' // OutputVariableGetName(cur_variable)
          if (len_trim(cur_variable%units) > 0) then
            word = cur_variable%units
            call HDF5MakeStringCompatible(word)
            string = trim(string) // ' [' // trim(word) // ']'
          endif

          call DiscretizationGlobalToNatural(discretization, &
                                             field%avg_vars_vec(ivar), &
                                             natural_vec,ONEDOF)
          call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                       H5T_NATIVE_DOUBLE)
          att_datasetname = trim(filename_header) // ":/" // trim(group_name) // &
                            "/" // trim(string)
          if (OptionIsIORank(option)) then
            call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                    att_datasetname,CELL_CENTERED_OUTPUT_MESH)
          endif
          cur_variable => cur_variable%next
        enddo
      endif

  end select

  !Output flowrates
  if (output_option%print_hdf5_mass_flowrate.or. &
     output_option%print_hdf5_energy_flowrate.or. &
     output_option%print_hdf5_aveg_mass_flowrate.or. &
     output_option%print_hdf5_aveg_energy_flowrate) then

    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        call OutputGetFaceFlowrateUGrid(realization_base)
        if (output_option%print_hdf5_mass_flowrate.or.&
           output_option%print_hdf5_energy_flowrate) then
          call WriteHDF5FlowratesUGrid(realization_base,option,grp_id, &
                                       var_list_type)
        endif
      case (AVERAGED_VARS)
        if (output_option%print_hdf5_aveg_mass_flowrate.or.&
           output_option%print_hdf5_aveg_energy_flowrate) then
          call WriteHDF5FlowratesUGrid(realization_base,option,grp_id, &
                                       var_list_type)
        endif
    end select
  endif

  if (output_option%print_hdf5_vel_cent .and. &
      (var_list_type==INSTANTANEOUS_VARS)) then

    ! velocities
    call OutputGetCellCenteredVelocities(realization_base,global_vec_vx, &
                                         global_vec_vy,global_vec_vz, &
                                         LIQUID_PHASE)
    do i = 1, 3
      select case(i)
        case(1)
          word = 'X'
          ivec = global_vec_vx
        case(2)
          word = 'Y'
          ivec = global_vec_vy
        case(3)
          word = 'Z'
          ivec = global_vec_vz
      end select
      string = 'Liquid ' // trim(word) // '-Velocity [m_per_' // &
               trim(output_option%tunit) // ']'
      call DiscretizationGlobalToNatural(discretization,ivec, &
                                         natural_vec,ONEDOF)
      call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                   H5T_NATIVE_DOUBLE)
      att_datasetname = trim(filename_header) // ":/" // &
                        trim(group_name) // "/" // trim(string)
      if (OptionIsIORank(option)) then
      call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string,att_datasetname, &
                              CELL_CENTERED_OUTPUT_MESH)
      endif
    enddo

    if (option%nphase > 1) then
        call OutputGetCellCenteredVelocities(realization_base,global_vec_vx, &
                                             global_vec_vy,global_vec_vz, &
                                             GAS_PHASE)
      do i = 1, 3
        select case(i)
          case(1)
            word = 'X'
            ivec = global_vec_vx
          case(2)
            word = 'Y'
            ivec = global_vec_vy
          case(3)
            word = 'Z'
            ivec = global_vec_vz
        end select
        string = 'Gas ' // trim(word) // '-Velocity [m_per_' // &
                 trim(output_option%tunit) // ']'
        call DiscretizationGlobalToNatural(discretization,ivec, &
                                         natural_vec,ONEDOF)
        call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                     H5T_NATIVE_DOUBLE)
        att_datasetname = trim(filename_header) // ":/" // &
                          trim(group_name) // "/" // trim(string)
        if (OptionIsIORank(option)) then
          call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                  att_datasetname,CELL_CENTERED_OUTPUT_MESH)
        endif
      enddo
    endif
  endif

  ! Output velocity at cell-face
  if (output_option%print_hdf5_vel_face) then

    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        call OutputGetFaceVelUGrid(realization_base)
        if (output_option%print_hdf5_vel_face) then
          call WriteHDF5FaceVelUGrid(realization_base,option,grp_id, &
                                     var_list_type)
        endif
      case (AVERAGED_VARS)
    end select
  endif

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vx,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vy,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vz,ierr);CHKERRQ(ierr)

  call HDF5GroupClose(grp_id,option)

  call HDF5FileClose(file_id,option)

  if (OptionIsIORank(option)) then
    call OutputXMFFooter(OUTPUT_UNIT)
    close(OUTPUT_UNIT)
  endif

  hdf5_first = PETSC_FALSE

end subroutine OutputHDF5UGridXDMF

! ************************************************************************** !

subroutine OutputHDF5UGridXDMFExplicit(realization_base,var_list_type)
  !
  ! This subroutine prints the explicit
  ! unstructured grid information in xdmf format
  !
  ! Author: Satish Karra, LANL
  ! Date: 07/17/2013
  !

  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  use String_module

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec
  use HDF5_Aux_module

  implicit none

  class(realization_base_type) :: realization_base
  PetscInt :: var_list_type

  integer(HID_T) :: file_id, new_file_id, file_id2
  integer(HID_T) :: grp_id, new_grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: data_set_id
  integer(HSIZE_T) :: dims(3), max_dims(3)

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable

  Vec :: global_vec
  Vec :: global_vec_vx, global_vec_vy, global_vec_vz
  Vec :: natural_vec
  PetscBool :: include_gas_phase

  character(len=MAXSTRINGLENGTH) :: filename_path, filename_header
  character(len=MAXSTRINGLENGTH) :: xmf_filename, att_datasetname, group_name
  character(len=MAXSTRINGLENGTH) :: domain_filename_path, domain_filename_header
  character(len=MAXSTRINGLENGTH) :: string, string2,string3
  character(len=MAXWORDLENGTH) :: word
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscMPIInt :: hdf5_err
  PetscBool :: first
  PetscInt :: ivar
  PetscBool :: write_xdmf
  PetscBool :: include_cell_centers
  PetscInt :: num_vertices, num_cells
  PetscInt :: mesh_type
  PetscErrorCode :: ierr

  discretization => realization_base%discretization
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option

  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      string2=''
      write(string3,'(i4)') output_option%plot_number
      xmf_filename = OutputFilename(output_option,option,'xmf','')
    case (AVERAGED_VARS)
      string2='-aveg'
      write(string3,'(i4)') &
        int(option%time/output_option%periodic_snap_output_time_incr)
      xmf_filename = OutputFilename(output_option,option,'xmf','aveg')
  end select

  if (output_option%print_single_h5_file) then
    first = hdf5_first
    !filename = trim(option%global_prefix) // trim(string2) // &
    !           trim(option%group_prefix) // '.h5'
    filename_path = trim(option%global_prefix) // trim(string2) // &
               trim(option%group_prefix) // '.h5'
    filename_header = trim(StringGetFilename(option%global_prefix)) // &
               trim(string2) // trim(option%group_prefix) // '.h5'
  else
    string = OutputHDF5FilenameID(output_option,option,var_list_type)
    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        if (mod(output_option%plot_number, &
                output_option%times_per_h5_file) == 0) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
      case (AVERAGED_VARS)
        if (Equal(mod((option%time-output_option%periodic_snap_output_time_incr)/ &
             output_option%periodic_snap_output_time_incr, &
             dble(output_option%times_per_h5_file)),0.d0)) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
    end select

    !filename = trim(option%global_prefix) // trim(option%group_prefix) // &
    !           trim(string2) // '-' // trim(string) // '.h5'
    filename_path = trim(option%global_prefix) // &
                    trim(option%group_prefix) // &
                    trim(string2) // '-' // trim(string) // '.h5'
    filename_header = trim(StringGetFilename(option%global_prefix)) // &
                    trim(option%group_prefix) // &
                    trim(string2) // '-' // trim(string) // '.h5'
  endif

  grid => patch%grid

  if (.not.first) then
    call HDF5FileTryOpen(filename_path,file_id,first,option%comm)
  endif
  if (first) then
    call HDF5FileOpen(filename_path,file_id,PETSC_TRUE,option)
  endif

  if (first) then
    option%io_buffer = ' --> creating hdf5 output file: ' // trim(filename_path)
  else
    option%io_buffer = ' --> appending to hdf5 output file: ' // &
                       trim(filename_path)
  endif
  call PrintMsg(option)

  domain_filename_path = trim(option%global_prefix) // '-domain.h5'
  domain_filename_header = &
    trim(StringGetFilename(option%global_prefix)) // '-domain.h5'
  write_xdmf = PETSC_FALSE
  include_cell_centers = PETSC_FALSE
  mesh_type = grid%unstructured_grid%explicit_grid%output_mesh_type
  if (len_trim(grid%unstructured_grid%explicit_grid%domain_filename) > 0 &
      .and. output_option%print_explicit_primal_grid) then
    option%io_buffer = 'PRINT_PRIMAL_GRID under OUTPUT may not be used &
      &when DOMAIN_FILENAME is defined under GRID. Please remove &
      &the DOMAIN_FILENAME card.'
    call PrintErrMsg(option)
  endif
  if (OptionIsIORank(option) .and. &
      (output_option%print_explicit_primal_grid .or. &
       len_trim(grid%unstructured_grid%explicit_grid% &
                  domain_filename) > 0)) then
    if (output_option%print_explicit_primal_grid) then
      ! for primal grid output, num_cells is set in the call to
      ! WriteHDF5CoordinatesUGridXDMFExplicit() below.  Therefore, this value
      ! for num_cells will be overwritten the first time called.
      num_cells = realization_base%output_option%xmf_vert_len
      num_vertices = grid%unstructured_grid%explicit_grid%num_vertices
    else
      ! have to open up domain file read the size of the vertex array
      domain_filename_path = &
        grid%unstructured_grid%explicit_grid%domain_filename

      domain_filename_header = domain_filename_path
      option%io_buffer = 'Opening HDF5 primary grid "' // &
        trim(domain_filename_path) // &
        '" for referencing in file "' // trim(xmf_filename) // '".'
      call PrintMsg(option)
      call HDF5FileOpenReadOnly(domain_filename_path,file_id2, &
                                PETSC_FALSE,'',option)
      string = 'Domain/Cells'
      call h5dopen_f(file_id2,string,data_set_id,hdf5_err)
      if (hdf5_err /= 0) then
        option%io_buffer = 'HDF5 dataset "' // trim(string) // '" not found &
          &in file "' // trim(domain_filename_path) // '".'
        call PrintErrMsg(option)
      endif
      call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
      ! should be a rank=2 data space
      call h5sget_simple_extent_dims_f(file_space_id,dims, &
                                       max_dims,hdf5_err)
      num_cells = int(dims(1))
      call h5sclose_f(file_space_id,hdf5_err)
      call HDF5DatasetClose(data_set_id,option)
      string = 'Domain/Vertices'
      call h5dopen_f(file_id2,string,data_set_id,hdf5_err)
      if (hdf5_err /= 0) then
        option%io_buffer = 'HDF5 dataset "' // trim(string) // '" not found &
          &in file "' // trim(domain_filename_path) // '".'
        call PrintErrMsg(option)
      endif
      call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
      ! should be a rank=2 data space
      call h5sget_simple_extent_dims_f(file_space_id,dims, &
                                       max_dims,hdf5_err)
      num_vertices = int(dims(2))
      call h5sclose_f(file_space_id,hdf5_err)
      call HDF5DatasetClose(data_set_id,option)
      call HDF5FileClose(file_id2,option)
      include_cell_centers = PETSC_TRUE
    endif
    write_xdmf = PETSC_TRUE
  endif

  if (OptionIsIORank(option) .and. &
      first .and. output_option%print_explicit_primal_grid) then
    call HDF5FileOpen(domain_filename_path,new_file_id,PETSC_TRUE,&
                      PETSC_FALSE,option)
    ! create a group for the coordinates data set
    string = "Domain"
    call HDF5GroupCreate(new_file_id,string,new_grp_id,option)
    call WriteHDF5CoordinatesUGridXDMFExplicit(realization_base,option, &
                                               new_grp_id)
    num_cells = realization_base%output_option%xmf_vert_len
    call HDF5GroupClose(new_grp_id,option)
    call HDF5FileClose(new_file_id,option)
  endif

  if (write_xdmf) then
    option%io_buffer = ' --> write xmf output file: ' // trim(xmf_filename)
    call PrintMsg(option)
    open(unit=OUTPUT_UNIT,file=xmf_filename,action="write")
    call OutputXMFHeader(OUTPUT_UNIT, &
                       option%time/output_option%tconv, &
                       grid%unstructured_grid%explicit_grid%num_elems, &
                       num_cells, &
                       num_vertices, &
                       domain_filename_header,include_cell_centers)
  endif

  ! create a group for the data set
  if (output_option%extend_hdf5_time_format) then
    write(string,'(''Time'',es20.12,x,a1)') &
          option%time/output_option%tconv,output_option%tunit
  else
    write(string,'(''Time'',es13.5,x,a1)') &
          option%time/output_option%tconv,output_option%tunit
  endif
  if (len_trim(output_option%plot_name) > 2) then
    string = trim(string) // ' ' // output_option%plot_name
  endif
  string = trim(string3) // ' ' // trim(string)

  call HDF5GroupOpenOrCreate(file_id,string,grp_id,option)
  group_name=string

  ! write group attributes
  call OutputHDF5WriteSnapShotAtts(grp_id,option)

  ! write out data sets
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)

  select case (var_list_type)

    case (INSTANTANEOUS_VARS)
      ! loop over snapshot variables and write to file
      cur_variable => output_option%output_snap_variable_list%first
      do
        if (.not.associated(cur_variable)) exit
        call OutputGetVariableArray(realization_base,global_vec,cur_variable)
        call DiscretizationGlobalToNatural(discretization,global_vec, &
                                           natural_vec,ONEDOF)
        string = OutputVariableGetName(cur_variable)
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif
        if (cur_variable%iformat == 0) then
          call HDF5WriteDataSetFromVec(string,option,natural_vec, &
                                       grp_id,H5T_NATIVE_DOUBLE)
        else
          call HDF5WriteDataSetFromVec(string,option,natural_vec, &
                                       grp_id,H5T_NATIVE_INTEGER)
        endif
        att_datasetname = trim(filename_header) // ":/" // &
                          trim(group_name) // "/" // trim(string)
        if (write_xdmf) then
          !call OutputXMFAttributeExplicit(OUTPUT_UNIT,grid%nmax,string, &
          !                                att_datasetname)
          call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                  att_datasetname,mesh_type)
        endif
        cur_variable => cur_variable%next
      enddo

    case (AVERAGED_VARS)
      if (associated(output_option%aveg_output_variable_list%first)) then
        cur_variable => output_option%aveg_output_variable_list%first
        do ivar = 1,output_option%aveg_output_variable_list%nvars
          string = 'Aveg. ' // OutputVariableGetName(cur_variable)
          if (len_trim(cur_variable%units) > 0) then
            word = cur_variable%units
            call HDF5MakeStringCompatible(word)
            string = trim(string) // ' [' // trim(word) // ']'
          endif

          call DiscretizationGlobalToNatural(discretization, &
                                             field%avg_vars_vec(ivar), &
                                             natural_vec,ONEDOF)
          call HDF5WriteDataSetFromVec(string,option,natural_vec, &
                                       grp_id,H5T_NATIVE_DOUBLE)
          att_datasetname = trim(filename_header) // ":/" // &
                            trim(group_name) // "/" // trim(string)
          if (write_xdmf) then
            call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                    att_datasetname,mesh_type)
          endif
          cur_variable => cur_variable%next
        enddo
      endif

  end select

  ! output cell-centered velocity
  include_gas_phase = PETSC_FALSE
  if (option%nphase > 1 .or. option%transport%nphase > 1) then
     include_gas_phase = PETSC_TRUE
  endif
  if (output_option%print_hdf5_vel_cent .and. &
       (var_list_type==INSTANTANEOUS_VARS)) then
     call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vx)
     call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vy)
     call DiscretizationDuplicateVector(discretization,global_vec,global_vec_vz)

     call OutputGetCellCenteredVelocities(realization_base, global_vec_vx, &
                                          global_vec_vy,global_vec_vz, &
                                          LIQUID_PHASE)

     string = "Liquid X-Velocity [m_per_" // trim(output_option%tunit) // "]"
     call DiscretizationGlobalToNatural(discretization,global_vec_vx, &
                                        natural_vec,ONEDOF)
     call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                  H5T_NATIVE_DOUBLE)
     att_datasetname = trim(filename_header) // ":/" // &
                       trim(group_name) // "/" // trim(string)
     if (write_xdmf) then
        call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                att_datasetname,mesh_type)
     endif

     string = "Liquid Y-Velocity [m_per_" // trim(output_option%tunit) // "]"
     call DiscretizationGlobalToNatural(discretization,global_vec_vy, &
                                        natural_vec,ONEDOF)
     call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                  H5T_NATIVE_DOUBLE)
     att_datasetname = trim(filename_header) // ":/" // &
                       trim(group_name) // "/" // trim(string)
     if (write_xdmf) then
        call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                att_datasetname,mesh_type)
     endif

     string = "Liquid Z-Velocity [m_per_" // trim(output_option%tunit) // "]"
     call DiscretizationGlobalToNatural(discretization,global_vec_vz, &
                                        natural_vec,ONEDOF)
     call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                  H5T_NATIVE_DOUBLE)

     att_datasetname = trim(filename_header) // ":/" // &
                       trim(group_name) // "/" // trim(string)
     if (write_xdmf) then
        call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                att_datasetname,mesh_type)
     endif

     if (include_gas_phase) then
        call OutputGetCellCenteredVelocities(realization_base,global_vec_vx, &
                                             global_vec_vy,global_vec_vz, &
                                             GAS_PHASE)

        string = "Gas X-Velocity"
        call DiscretizationGlobalToNatural(discretization,global_vec_vx, &
                                           natural_vec,ONEDOF)
        call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                    H5T_NATIVE_DOUBLE)
        att_datasetname = trim(filename_header) // ":/" // &
                          trim(group_name) // "/" // trim(string)
        if (write_xdmf) then
           call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                   att_datasetname,mesh_type)
        endif

        string = "Gas Y-Velocity"
        call DiscretizationGlobalToNatural(discretization,global_vec_vy, &
                                           natural_vec,ONEDOF)
        call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                     H5T_NATIVE_DOUBLE)
        att_datasetname = trim(filename_header) // ":/" // &
                          trim(group_name) // "/" // trim(string)
        if (write_xdmf) then
           call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                   att_datasetname,mesh_type)
        endif

        string = "Gas Z-Velocity"
        call DiscretizationGlobalToNatural(discretization,global_vec_vz, &
                                           natural_vec,ONEDOF)
        call HDF5WriteDataSetFromVec(string,option,natural_vec,grp_id, &
                                     H5T_NATIVE_DOUBLE)
        att_datasetname = trim(filename_header) // ":/" // &
                          trim(group_name) // "/" // trim(string)
        if (write_xdmf) then
           call OutputXMFAttribute(OUTPUT_UNIT,grid%nmax,string, &
                                   att_datasetname,mesh_type)
        endif
     endif
     call VecDestroy(global_vec_vx,ierr);CHKERRQ(ierr)
     call VecDestroy(global_vec_vy,ierr);CHKERRQ(ierr)
     call VecDestroy(global_vec_vz,ierr);CHKERRQ(ierr)
  endif

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)

  call HDF5GroupClose(grp_id,option)
  call HDF5FileClose(file_id,option)

  if (write_xdmf) then
    call OutputXMFFooter(OUTPUT_UNIT)
    close(OUTPUT_UNIT)
  endif

  hdf5_first = PETSC_FALSE

end subroutine OutputHDF5UGridXDMFExplicit

! ************************************************************************** !

function OutputHDF5FilenameID(output_option,option,var_list_type)
  !
  ! This subroutine creates an ID for HDF5 filename for:
  ! - Instantaneous, or
  ! - Temporally averaged variables.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/10/13
  !
  use Option_module

  implicit none

  type(option_type) :: option
  type(output_option_type) :: output_option
  PetscInt :: var_list_type

  character(len=MAXWORDLENGTH) :: OutputHDF5FilenameID
  PetscInt :: file_number

  select case(var_list_type)
    case (INSTANTANEOUS_VARS)
      file_number = floor(real(output_option%plot_number)/ &
                               output_option%times_per_h5_file)
    case (AVERAGED_VARS)
      file_number = floor((option%time - &
                           output_option%periodic_snap_output_time_incr)/ &
                           output_option%periodic_snap_output_time_incr/ &
                           output_option%times_per_h5_file)
  end select

  if (file_number < 10) then
    write(OutputHDF5FilenameID,'("00",i1)') file_number
  else if (output_option%plot_number < 100) then
    write(OutputHDF5FilenameID,'("0",i2)') file_number
  else if (output_option%plot_number < 1000) then
    write(OutputHDF5FilenameID,'(i3)') file_number
  else if (output_option%plot_number < 10000) then
    write(OutputHDF5FilenameID,'(i4)') file_number
  else if (output_option%plot_number < 100000) then
    write(OutputHDF5FilenameID,'(i5)') file_number
  else
    option%io_buffer = 'Plot number exceeds current maximum of 10^5.'
    call PrintErrMsgToDev(option,'ask for a higher maximum')
  endif

  OutputHDF5FilenameID = adjustl(OutputHDF5FilenameID)

end function OutputHDF5FilenameID

! ************************************************************************** !

subroutine WriteHDF5FluxVelocities(name,realization_base,iphase,direction, &
                                   file_id)
  !
  ! Print flux velocities to HDF5 file
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Connection_module
  use Patch_module
  use hdf5
  use HDF5_module, only : HDF5WriteStructuredDataSet, trick_hdf5

  implicit none

  character(len=32) :: name
  class(realization_base_type) :: realization_base
  PetscInt :: iphase
  PetscInt :: direction
  integer(HID_T) :: file_id

  PetscInt :: i, j, k
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option

  PetscReal, allocatable :: array(:)

  PetscBool, save :: trick_flux_vel_x = PETSC_FALSE
  PetscBool, save :: trick_flux_vel_y = PETSC_FALSE
  PetscBool, save :: trick_flux_vel_z = PETSC_FALSE

  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option

  ! in a few cases (i.e. for small test problems), some processors may
  ! have no velocities to print.  This results in zero-length arrays
  ! in collective H5Dwrite().  To avoid, we switch to independent
  ! H5Dwrite() and don't write from the zero-length procs.
!GEH - Structured Grid Dependence - Begin
  if (hdf5_first) then
    trick_flux_vel_x = PETSC_FALSE
    trick_flux_vel_y = PETSC_FALSE
    trick_flux_vel_z = PETSC_FALSE

    nx_local = grid%structured_grid%nlx
    ny_local = grid%structured_grid%nly
    nz_local = grid%structured_grid%nlz
    if (grid%structured_grid%gxe-grid%structured_grid%lxe == 0) then
      nx_local = grid%structured_grid%nlx-1
    endif
    call MPI_Allreduce(nx_local,i,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr);CHKERRQ(ierr)
    if (i == 0) trick_flux_vel_x = PETSC_TRUE
    if (grid%structured_grid%gye-grid%structured_grid%lye == 0) then
      ny_local = grid%structured_grid%nly-1
    endif
    call MPI_Allreduce(ny_local,j,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr);CHKERRQ(ierr)
    if (j == 0) trick_flux_vel_y = PETSC_TRUE
    if (grid%structured_grid%gze-grid%structured_grid%lze == 0) then
      nz_local = grid%structured_grid%nlz-1
    endif
    call MPI_Allreduce(nz_local,k,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr);CHKERRQ(ierr)
    if (k == 0) trick_flux_vel_z = PETSC_TRUE
  endif

  nx_local = grid%structured_grid%nlx
  ny_local = grid%structured_grid%nly
  nz_local = grid%structured_grid%nlz
  nx_global = grid%structured_grid%nx
  ny_global = grid%structured_grid%ny
  nz_global = grid%structured_grid%nz

  select case(direction)
    case(X_DIRECTION)
      nx_global = grid%structured_grid%nx-1
      if (grid%structured_grid%gxe-grid%structured_grid%lxe == 0) then
        nx_local = grid%structured_grid%nlx-1
      endif
      if (trick_flux_vel_x) trick_hdf5 = PETSC_TRUE
    case(Y_DIRECTION)
      ny_global = grid%structured_grid%ny-1
      if (grid%structured_grid%gye-grid%structured_grid%lye == 0) then
        ny_local = grid%structured_grid%nly-1
      endif
      if (trick_flux_vel_y) trick_hdf5 = PETSC_TRUE
    case(Z_DIRECTION)
      nz_global = grid%structured_grid%nz-1
      if (grid%structured_grid%gze-grid%structured_grid%lze == 0) then
        nz_local = grid%structured_grid%nlz-1
      endif
      if (trick_flux_vel_z) trick_hdf5 = PETSC_TRUE
  end select

  allocate(array(nx_local*ny_local*nz_local))
  call OutputCollectVelocityOrFlux(realization_base, iphase, direction, &
                                   PETSC_FALSE, array)

  array(1:nx_local*ny_local*nz_local) = &  ! convert time units
    array(1:nx_local*ny_local*nz_local) * output_option%tconv

  call HDF5WriteStructuredDataSet(name,array,file_id,H5T_NATIVE_DOUBLE, &
                                  option,nx_global,ny_global,nz_global, &
                                  nx_local,ny_local,nz_local, &
                                  grid%structured_grid%lxs, &
                                  grid%structured_grid%lys, &
                                  grid%structured_grid%lzs)
!GEH - Structured Grid Dependence - End

  deallocate(array)
  trick_hdf5 = PETSC_FALSE

end subroutine WriteHDF5FluxVelocities

! ************************************************************************** !

subroutine OutputHDF5WriteStructCoordGroup(file_id,discretization,grid,option)
  !
  ! Writes the Coordinates group to an structured HDF5 output file
  !
  ! Author: Glenn Hammond
  ! Date: 10/12/21
  !
  use hdf5
  use HDF5_Aux_module
  use Discretization_module
  use Option_module
  use Grid_module
  use String_module

  implicit none

  integer(HID_T) :: file_id
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option

  integer(HID_T) :: grp_id
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: array(:)
  PetscInt :: i

  ! create a group for the coordinates data set
  string = "Coordinates"
  call HDF5GroupCreate(file_id,string,grp_id,option)

  !GEH - Structured Grid Dependence - Begin
  ! write out coordinates in x, y, and z directions
  string = "X [m]"
  allocate(array(grid%structured_grid%nx+1))
  array(1) = discretization%origin_global(X_DIRECTION)
  do i=2,grid%structured_grid%nx+1
    array(i) = array(i-1) + grid%structured_grid%dx_global(i-1)
  enddo
  call WriteHDF5Coordinates(string,option,grid%structured_grid%nx+1, &
                            array,grp_id)
  deallocate(array)

  string = "Y [m]"
  allocate(array(grid%structured_grid%ny+1))
  array(1) = discretization%origin_global(Y_DIRECTION)
  do i=2,grid%structured_grid%ny+1
    array(i) = array(i-1) + grid%structured_grid%dy_global(i-1)
  enddo
  call WriteHDF5Coordinates(string,option,grid%structured_grid%ny+1, &
                            array,grp_id)
  deallocate(array)

  string = "Z [m]"
  allocate(array(grid%structured_grid%nz+1))
  array(1) = discretization%origin_global(Z_DIRECTION)
  do i=2,grid%structured_grid%nz+1
    array(i) = array(i-1) + grid%structured_grid%dz_global(i-1)
  enddo
  call WriteHDF5Coordinates(string,option,grid%structured_grid%nz+1, &
                            array,grp_id)
  deallocate(array)
  !GEH - Structured Grid Dependence - End

  call HDF5GroupClose(grp_id,option)

end subroutine OutputHDF5WriteStructCoordGroup

! ************************************************************************** !

subroutine WriteHDF5Coordinates(name,option,length,array,file_id)
  !
  ! Writes structured coordinates to HDF5 file
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

  use hdf5
  use HDF5_Aux_module
  use Option_module

  implicit none

  character(len=32) :: name
  type(option_type) :: option
  PetscInt :: length
  PetscReal :: array(:)
  integer(HID_T) :: file_id

  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  PetscMPIInt :: rank
  integer :: hdf5_err
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_output_coordinates_hdf5, &
                          ierr);CHKERRQ(ierr)

  ! write out grid structure
  rank = 1
  dims = 0
  ! x-direction
  dims(1) = length
  call h5screate_simple_f(rank,dims,file_space_id,hdf5_err,dims)
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,file_space_id, &
                   data_set_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err) ! must be independent and only from p0
#endif
  if (OptionIsIORank(option)) then
     call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
     ! this is due to a bug in hdf5-1.8.18 hwere H5S_ALL_F is an INTEGER.  It
     ! should be INTEGER(HID_T)
     memory_space_id = H5S_ALL_F
     call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,array,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
     call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  endif
  call h5pclose_f(prop_id,hdf5_err)
  call HDF5DatasetClose(data_set_id,option)
  call h5sclose_f(file_space_id,hdf5_err)

  call PetscLogEventEnd(logging%event_output_coordinates_hdf5, &
                        ierr);CHKERRQ(ierr)

end subroutine WriteHDF5Coordinates

! ************************************************************************** !

subroutine WriteHDF5CoordinatesUGrid(grid,option,file_id)
  !
  ! This subroutine writes unstructured coordinates to HDF5 file
  !
  ! Author: Gautam Bisht, ORNL
  ! Date: 05/31/12
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use HDF5_Aux_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use HDF5_module, only : trick_hdf5
  use Variables_module

  implicit none

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option

  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt :: hdf5_err

  PetscInt :: local_size
  PetscInt :: istart
  PetscInt :: i,j
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)
  Vec :: global_x_vertex_vec,global_y_vertex_vec,global_z_vertex_vec

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  type(ugdm_type),pointer :: ugdm_element
  PetscErrorCode :: ierr

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_x_vertex_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_y_vertex_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_z_vertex_vec,ierr);CHKERRQ(ierr)

  call VecGetLocalSize(global_x_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_y_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_z_vertex_vec,local_size,ierr);CHKERRQ(ierr)

  call OutputGetVertexCoordinates(grid, global_x_vertex_vec,X_COORDINATE,option)
  call OutputGetVertexCoordinates(grid, global_y_vertex_vec,Y_COORDINATE,option)
  call OutputGetVertexCoordinates(grid, global_z_vertex_vec,Z_COORDINATE,option)

  call VecGetArrayF90(global_x_vertex_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_y_vertex_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_z_vertex_vec,vec_z_ptr,ierr);CHKERRQ(ierr)

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size * 3
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 3D block
  rank_mpi = 2
  dims = 0
  dims(2) = grid%unstructured_grid%num_vertices_global
  dims(1) = 3
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Vertices" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                  option%mycomm,ierr);CHKERRQ(ierr)

  start(2) = istart
  start(1) = 0

  length(2) = local_size
  length(1) = 3

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif

  allocate(double_array(local_size*3))

  do i=1,local_size
    double_array((i-1)*3+1) = vec_x_ptr(i)
    double_array((i-1)*3+2) = vec_y_ptr(i)
    double_array((i-1)*3+3) = vec_z_ptr(i)
  enddo

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  deallocate(double_array)
  call h5pclose_f(prop_id,hdf5_err)


  call HDF5DatasetClose(data_set_id,option)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(global_x_vertex_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_y_vertex_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_z_vertex_vec,vec_z_ptr,ierr);CHKERRQ(ierr)


  call VecDestroy(global_x_vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_y_vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_z_vertex_vec,ierr);CHKERRQ(ierr)


  !
  !  Write elements
  !
  call UGridCreateUGDM(grid%unstructured_grid,ugdm_element,EIGHT_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,natural_vec, &
                           NATURAL,option)
  call OutputGetCellVertices(grid,global_vec)
  call VecScatterBegin(ugdm_element%scatter_gton,global_vec,natural_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_element%scatter_gton,global_vec,natural_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)

  local_size = grid%unstructured_grid%nlmax

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size*NINE_INTEGER
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 3D block
  rank_mpi = 2
  dims = 0
  dims(2) = grid%unstructured_grid%nmax
  dims(1) = NINE_INTEGER
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Cells" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_INTEGER,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                  option%mycomm,ierr);CHKERRQ(ierr)

  start(2) = istart
  start(1) = 0

  length(2) = local_size
  length(1) = NINE_INTEGER

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif

  allocate(int_array(local_size*NINE_INTEGER))

  do i=1,local_size
    int_array((i-1)*9 + 1) = 0
    int_array((i-1)*9 + 2) = INT(vec_ptr((i-1)*8+1))
    int_array((i-1)*9 + 3) = INT(vec_ptr((i-1)*8+2))
    int_array((i-1)*9 + 4) = INT(vec_ptr((i-1)*8+3))
    int_array((i-1)*9 + 5) = INT(vec_ptr((i-1)*8+4))
    int_array((i-1)*9 + 6) = INT(vec_ptr((i-1)*8+5))
    int_array((i-1)*9 + 7) = INT(vec_ptr((i-1)*8+6))
    int_array((i-1)*9 + 8) = INT(vec_ptr((i-1)*8+7))
    int_array((i-1)*9 + 9) = INT(vec_ptr((i-1)*8+8))
    do j=2,9
      if (int_array((i-1)*9 + j)>0) int_array((i-1)*9 + 1)= int_array((i-1)*9 + 1) +1
    enddo
  enddo

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_INTEGER,int_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  deallocate(int_array)
  call h5pclose_f(prop_id,hdf5_err)


  call HDF5DatasetClose(data_set_id,option)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call UGridDMDestroy(ugdm_element)

end subroutine WriteHDF5CoordinatesUGrid

! ************************************************************************** !

subroutine WriteHDF5CoordinatesUGridXDMF(realization_base,option,file_id)
  !
  ! This routine writes unstructured coordinates to HDF5 file in XDMF format
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 10/29/2012
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use HDF5_Aux_module
  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Variables_module

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option

  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  type(grid_type), pointer :: grid
  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt :: hdf5_err

  PetscInt :: local_size,vert_count,nverts
  PetscInt :: i,j
  PetscInt :: temp_int, istart
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)
  Vec :: global_x_vertex_vec,global_y_vertex_vec,global_z_vertex_vec
  Vec :: global_x_cell_vec,global_y_cell_vec,global_z_cell_vec
  Vec :: natural_x_cell_vec,natural_y_cell_vec,natural_z_cell_vec

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  type(ugdm_type),pointer :: ugdm_element, ugdm_cell
  PetscErrorCode :: ierr

  PetscInt :: TET_ID_XDMF = 6
  PetscInt :: PYR_ID_XDMF = 7
  PetscInt :: WED_ID_XDMF = 8
  PetscInt :: HEX_ID_XDMF = 9

  grid => realization_base%patch%grid

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_x_vertex_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_y_vertex_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%unstructured_grid%num_vertices_global, &
                    global_z_vertex_vec,ierr);CHKERRQ(ierr)

  call VecGetLocalSize(global_x_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_y_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_z_vertex_vec,local_size,ierr);CHKERRQ(ierr)

  call OutputGetVertexCoordinates(grid, global_x_vertex_vec,X_COORDINATE,option)
  call OutputGetVertexCoordinates(grid, global_y_vertex_vec,Y_COORDINATE,option)
  call OutputGetVertexCoordinates(grid, global_z_vertex_vec,Z_COORDINATE,option)

  call VecGetArrayF90(global_x_vertex_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_y_vertex_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_z_vertex_vec,vec_z_ptr,ierr);CHKERRQ(ierr)

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size * 3
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 2D block
  rank_mpi = 2
  dims = 0
  dims(2) = grid%unstructured_grid%num_vertices_global
  dims(1) = 3
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Vertices" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  !geh: cannot use dims(1) in MPI_Allreduce as it causes errors on
  !     Juqueen
  call MPI_Exscan(local_size,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                  option%mycomm,ierr);CHKERRQ(ierr)

  start(2) = istart
  start(1) = 0

  length(2) = local_size
  length(1) = 3

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  allocate(double_array(local_size*3))
  do i=1,local_size
    double_array((i-1)*3+1) = vec_x_ptr(i)
    double_array((i-1)*3+2) = vec_y_ptr(i)
    double_array((i-1)*3+3) = vec_z_ptr(i)
  enddo

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  deallocate(double_array)
  call h5pclose_f(prop_id,hdf5_err)

  call HDF5DatasetClose(data_set_id,option)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(global_x_vertex_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_y_vertex_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_z_vertex_vec,vec_z_ptr,ierr);CHKERRQ(ierr)


  call VecDestroy(global_x_vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_y_vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_z_vertex_vec,ierr);CHKERRQ(ierr)

  !
  !  Write elements
  !
  call UGridCreateUGDM(grid%unstructured_grid,ugdm_element,EIGHT_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,natural_vec, &
                           NATURAL,option)
  call OutputGetCellVertices(grid,global_vec)
  call VecScatterBegin(ugdm_element%scatter_gton,global_vec,natural_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_element%scatter_gton,global_vec,natural_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)

  local_size = grid%unstructured_grid%nlmax

  vert_count=0
  do i=1,local_size*EIGHT_INTEGER
    if (int(vec_ptr(i)) >0 ) vert_count=vert_count+1
  enddo
  vert_count=vert_count+grid%nlmax

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims(1) = vert_count
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  !geh: cannot use dims(1) in MPI_Allreduce as it causes errors on
  !     Juqueen
  call MPI_Allreduce(vert_count,temp_int,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                     option%mycomm,ierr);CHKERRQ(ierr)
  dims(1) = temp_int
  realization_base%output_option%xmf_vert_len=int(dims(1))

  ! file space which is a 2D block
  rank_mpi = 1
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Cells" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_INTEGER,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(vert_count,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                  option%mycomm,ierr);CHKERRQ(ierr)

  start(1) = istart
  length(1) = vert_count
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  allocate(int_array(vert_count))

  vert_count=0
  do i=1,local_size
    nverts=0
    do j=1,8
      if (vec_ptr((i-1)*8+j)>0) nverts=nverts+1
    enddo
    vert_count=vert_count+1
    select case (nverts)
      case (4) ! Tetrahedron
        int_array(vert_count) = TET_ID_XDMF
      case (5) ! Pyramid
        int_array(vert_count) = PYR_ID_XDMF
      case (6) ! Wedge
        int_array(vert_count) = WED_ID_XDMF
      case (8) ! Hexahedron
        int_array(vert_count) = HEX_ID_XDMF
    end select

    do j=1,8
      if (vec_ptr((i-1)*8+j)>0) then
        vert_count=vert_count+1
        int_array(vert_count) = INT(vec_ptr((i-1)*8+j))-1
      endif
    enddo
  enddo

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_INTEGER,int_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  deallocate(int_array)
  call h5pclose_f(prop_id,hdf5_err)

  call HDF5DatasetClose(data_set_id,option)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call UGridDMDestroy(ugdm_element)

  ! Cell center X/Y/Z
  call VecCreateMPI(option%mycomm,grid%nlmax,PETSC_DETERMINE, &
                    global_x_cell_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,grid%nlmax,PETSC_DETERMINE, &
                    global_y_cell_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,grid%nlmax,PETSC_DETERMINE, &
                    global_z_cell_vec,ierr);CHKERRQ(ierr)

  call OutputGetCellCoordinates(grid, global_x_cell_vec,X_COORDINATE)
  call OutputGetCellCoordinates(grid, global_y_cell_vec,Y_COORDINATE)
  call OutputGetCellCoordinates(grid, global_z_cell_vec,Z_COORDINATE)


  call UGridCreateUGDM(grid%unstructured_grid,ugdm_cell,ONE_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_cell, &
                           natural_x_cell_vec,NATURAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_cell, &
                           natural_y_cell_vec,NATURAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_cell, &
                           natural_z_cell_vec,NATURAL,option)

  call VecScatterBegin(ugdm_cell%scatter_gton,global_x_cell_vec, &
                       natural_x_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_cell%scatter_gton,global_x_cell_vec, &
                     natural_x_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                     ierr);CHKERRQ(ierr)

  call VecScatterBegin(ugdm_cell%scatter_gton,global_y_cell_vec, &
                       natural_y_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_cell%scatter_gton,global_y_cell_vec, &
                     natural_y_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                     ierr);CHKERRQ(ierr)

  call VecScatterBegin(ugdm_cell%scatter_gton,global_z_cell_vec, &
                       natural_z_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_cell%scatter_gton,global_z_cell_vec, &
                     natural_z_cell_vec,INSERT_VALUES,SCATTER_FORWARD, &
                     ierr);CHKERRQ(ierr)

  call VecGetArrayF90(natural_x_cell_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(natural_y_cell_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(natural_z_cell_vec,vec_z_ptr,ierr);CHKERRQ(ierr)
  local_size = grid%unstructured_grid%nlmax

  ! XC
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 2D block
  rank_mpi = 1
  dims = 0
  dims(1) = grid%nmax
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "XC" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                  option%mycomm,ierr);CHKERRQ(ierr)
  start(1) = istart
  length(1) = local_size

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,vec_x_ptr,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  call h5pclose_f(prop_id,hdf5_err)

  call HDF5DatasetClose(data_set_id,option)
  call h5sclose_f(file_space_id,hdf5_err)

  ! YC
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 2D block
  rank_mpi = 1
  dims = 0
  dims(1) = grid%nmax
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "YC" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                  option%mycomm,ierr);CHKERRQ(ierr)
  start(1) = istart
  length(1) = local_size

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,vec_y_ptr,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  call h5pclose_f(prop_id,hdf5_err)

  call HDF5DatasetClose(data_set_id,option)
  call h5sclose_f(file_space_id,hdf5_err)

  ! ZC
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 2D block
  rank_mpi = 1
  dims = 0
  dims(1) = grid%nmax
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "ZC" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                  option%mycomm,ierr);CHKERRQ(ierr)
  start(1) = istart
  length(1) = local_size

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,vec_z_ptr,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  call h5pclose_f(prop_id,hdf5_err)

  call HDF5DatasetClose(data_set_id,option)
  call h5sclose_f(file_space_id,hdf5_err)


  call VecRestoreArrayF90(natural_x_cell_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(natural_y_cell_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(natural_z_cell_vec,vec_z_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(global_x_cell_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_y_cell_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_z_cell_vec,ierr);CHKERRQ(ierr)

  call VecDestroy(natural_x_cell_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_y_cell_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_z_cell_vec,ierr);CHKERRQ(ierr)

  call UGridDMDestroy(ugdm_cell)

end subroutine WriteHDF5CoordinatesUGridXDMF

! ************************************************************************** !

subroutine DetermineNumVertices(realization_base,option)
  !
  ! Determine the number of vertices written out in the output HDF5 file
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/13/2015
  !
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Variables_module

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option

  type(grid_type), pointer :: grid
  PetscInt :: local_size,vert_count
  PetscInt :: i
  PetscInt :: temp_int

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  type(ugdm_type),pointer :: ugdm_element
  PetscErrorCode :: ierr

  grid => realization_base%patch%grid

  call UGridCreateUGDM(grid%unstructured_grid,ugdm_element,EIGHT_INTEGER, &
                       option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm_element,natural_vec, &
                           NATURAL,option)
  call OutputGetCellVertices(grid,global_vec)
  call VecScatterBegin(ugdm_element%scatter_gton,global_vec,natural_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm_element%scatter_gton,global_vec,natural_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  local_size = grid%unstructured_grid%nlmax

  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  vert_count=0
  do i=1,local_size*EIGHT_INTEGER
    if (int(vec_ptr(i)) >0 ) vert_count=vert_count+1
  enddo
  vert_count=vert_count+grid%nlmax
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)

  call MPI_Allreduce(vert_count,temp_int,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                     option%mycomm,ierr);CHKERRQ(ierr)
  realization_base%output_option%xmf_vert_len=temp_int

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call UGridDMDestroy(ugdm_element)

end subroutine DetermineNumVertices

! ************************************************************************** !

subroutine WriteHDF5CoordinatesUGridXDMFExplicit(realization_base,option, &
                                                 file_id)
  !
  ! Writes the coordinates of
  ! explicit grid to HDF5 file
  !
  ! Author: Satish Karra, LANL
  ! Date: 07/17/2013
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use HDF5_Aux_module
  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Variables_module

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option

  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  type(grid_type), pointer :: grid
  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt :: hdf5_err

  PetscInt :: local_size,vert_count,nverts
  PetscInt :: i,j
  PetscInt :: istart
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)

  PetscReal, pointer :: vec_ptr(:)
  Vec :: natural_vec
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  PetscErrorCode :: ierr

  PetscInt :: TRI_ID_XDMF = 4
  PetscInt :: QUAD_ID_XDMF = 5
  PetscInt :: TET_ID_XDMF = 6
  PetscInt :: PYR_ID_XDMF = 7
  PetscInt :: WED_ID_XDMF = 8
  PetscInt :: HEX_ID_XDMF = 9

  grid => realization_base%patch%grid

  allocate(vec_x_ptr(grid%unstructured_grid%explicit_grid%num_vertices))
  allocate(vec_y_ptr(grid%unstructured_grid%explicit_grid%num_vertices))
  allocate(vec_z_ptr(grid%unstructured_grid%explicit_grid%num_vertices))

  do i = 1, grid%unstructured_grid%explicit_grid%num_vertices
    vec_x_ptr(i) = grid%unstructured_grid%explicit_grid%vertex_coordinates(i)%x
    vec_y_ptr(i) = grid%unstructured_grid%explicit_grid%vertex_coordinates(i)%y
    vec_z_ptr(i) = grid%unstructured_grid%explicit_grid%vertex_coordinates(i)%z
  enddo

  !local_size = grid%unstructured_grid%explicit_grid%num_cells_global
  local_size = grid%unstructured_grid%explicit_grid%num_vertices
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size * 3
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 2D block
  rank_mpi = 2
  dims = 0
  !POc
  !dims(2) = grid%unstructured_grid%explicit_grid%num_cells_global
  dims(2) = grid%unstructured_grid%explicit_grid%num_vertices
  dims(1) = 3
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Vertices" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  start(2) = istart
  start(1) = 0

  length(2) = local_size
  length(1) = 3

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)


  allocate(double_array(local_size*3))
  do i=1,local_size
    double_array((i-1)*3+1) = vec_x_ptr(i)
    double_array((i-1)*3+2) = vec_y_ptr(i)
    double_array((i-1)*3+3) = vec_z_ptr(i)
  enddo

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  deallocate(double_array)

  call h5pclose_f(prop_id,hdf5_err)

  call HDF5DatasetClose(data_set_id,option)
  call h5sclose_f(file_space_id,hdf5_err)

  deallocate(vec_x_ptr)
  deallocate(vec_y_ptr)
  deallocate(vec_z_ptr)

  !
  !  Write elements
  !
  local_size = grid%unstructured_grid%explicit_grid%num_elems

  call VecCreate(PETSC_COMM_SELF,natural_vec,ierr);CHKERRQ(ierr)
  call VecSetSizes(natural_vec,local_size*EIGHT_INTEGER,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetBlockSize(natural_vec,EIGHT_INTEGER,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(natural_vec,ierr);CHKERRQ(ierr)

  call OutputGetCellVerticesExplicit(grid,natural_vec)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)

  vert_count=0

  do i=1,local_size*EIGHT_INTEGER
    if (int(vec_ptr(i)) >0 ) vert_count=vert_count+1
  enddo
  vert_count=vert_count+grid%unstructured_grid%explicit_grid%num_elems

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims(1) = vert_count
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  realization_base%output_option%xmf_vert_len=int(dims(1))

  ! file space which is a 2D block
  rank_mpi = 1
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Cells" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_INTEGER,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  start(1) = istart
  length(1) = vert_count
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)

  allocate(int_array(vert_count))

  vert_count=0
  do i=1,local_size
    nverts=0
    do j=1,8
      if (vec_ptr((i-1)*8+j)>0) nverts=nverts+1
    enddo
    vert_count=vert_count+1
    select case (nverts)
      case (3)
        int_array(vert_count) = TRI_ID_XDMF
      case (4)
        if (grid%unstructured_grid%grid_type /= TWO_DIM_GRID) then
        ! Tetrahedron
          int_array(vert_count) = TET_ID_XDMF
        else
          int_array(vert_count) = QUAD_ID_XDMF
        endif
      case (5) ! Pyramid
        int_array(vert_count) = PYR_ID_XDMF
      case (6) ! Wedge
        int_array(vert_count) = WED_ID_XDMF
      case (8) ! Hexahedron
        int_array(vert_count) = HEX_ID_XDMF
    end select

    do j=1,8
      if (vec_ptr((i-1)*8+j)>0) then
        vert_count=vert_count+1
        int_array(vert_count) = INT(vec_ptr((i-1)*8+j))-1
      endif
    enddo
  enddo

  call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_INTEGER,int_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

  deallocate(int_array)

  call h5pclose_f(prop_id,hdf5_err)

  call HDF5DatasetClose(data_set_id,option)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)

end subroutine WriteHDF5CoordinatesUGridXDMFExplicit

! ************************************************************************** !

subroutine WriteHDF5FlowratesUGrid(realization_base,option,file_id, &
                                   var_list_type)
  !
  ! This routine writes (mass/energy) flowrate for unstructured grid.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/19/2013
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  use Variables_module
  use Connection_module
  use Coupler_module
  use HDF5_Aux_module
  use Output_Aux_module
  use Field_module

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  PetscInt :: var_list_type

  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field

  PetscInt :: dof
  PetscInt :: offset
  PetscInt :: local_size
  PetscInt :: i
  PetscInt :: iface
  PetscInt :: ndof

  PetscReal, pointer :: vec_ptr1(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscReal, pointer :: double_array(:)
  PetscInt :: istart

  PetscBool :: mass_flowrate
  PetscBool :: energy_flowrate

  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string

  patch => realization_base%patch
  grid => patch%grid
  ugrid => grid%unstructured_grid
  output_option =>realization_base%output_option
  field => realization_base%field

  select case(option%iflowmode)
    case (RICHARDS_MODE,RICHARDS_TS_MODE,ZFLOW_MODE,PNF_MODE)
      ndof=1
    case (TH_MODE,TH_TS_MODE)
      ndof=1
      if (output_option%print_hdf5_mass_flowrate .and. &
          output_option%print_hdf5_energy_flowrate) ndof = 2
    case default
      option%io_buffer='FLOWRATE output not supported in this mode'
      call PrintErrMsg(option)
  end select

  call VecGetLocalSize(field%flowrate_inst,local_size,ierr);CHKERRQ(ierr)
  local_size = local_size/(option%nflowdof*MAX_FACE_PER_CELL + 1)

  allocate(double_array(local_size*(MAX_FACE_PER_CELL+1)))
  double_array = 0.d0

  offset = option%nflowdof*MAX_FACE_PER_CELL+1

  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      call VecGetArrayF90(field%flowrate_inst,vec_ptr1,ierr);CHKERRQ(ierr)
      mass_flowrate = output_option%print_hdf5_mass_flowrate
      energy_flowrate = output_option%print_hdf5_energy_flowrate
    case (AVERAGED_VARS)
      call VecGetArrayF90(field%flowrate_inst,vec_ptr1,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(field%flowrate_aveg,vec_ptr2,ierr);CHKERRQ(ierr)
      mass_flowrate = output_option%print_hdf5_aveg_mass_flowrate
      energy_flowrate = output_option%print_hdf5_aveg_energy_flowrate
  end select


  do dof = 1,option%nflowdof

    if (dof==1 .and. (.not.mass_flowrate)) exit
    if (dof==2 .and. (.not.energy_flowrate)) exit

    select case(option%iflowmode)
      case(RICHARDS_MODE,RICHARDS_TS_MODE,PNF_MODE)
        string = "Mass_Flowrate [kg_per_s]" // CHAR(0)
      case(ZFLOW_MODE)
        string = "Mass_Flowrate [m^3_per_s]" // CHAR(0)
        option%io_buffer = 'Fix mass flow rate for zflow in output_hdf5.F90'
        call PrintErrMsg(option)
      case(TH_MODE,TH_TS_MODE)
        if (dof==1) then
          string = "Mass_Flowrate [kg_per_s]" // CHAR(0)
        else
          string = "Energy_Flowrate [MJ_per_s]" // CHAR(0)
        endif
      case default
        option%io_buffer='FLOWRATE output not implemented in this mode.'
        call PrintErrMsg(option)
    end select

    if (var_list_type==AVERAGED_VARS) string = 'Aveg_' // trim(string) // &
                                               char(0)

    ! memory space which is a 1D vector
    rank_mpi = 1
    dims = 0
    dims(1) = local_size*(MAX_FACE_PER_CELL+1)
    call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

    ! file space which is a 2D block
    rank_mpi = 2
    dims = 0
    dims(2) = ugrid%nmax
    dims(1) = MAX_FACE_PER_CELL + 1
    call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

    call h5eset_auto_f(OFF,hdf5_err)
    call h5dopen_f(file_id,trim(string),data_set_id,hdf5_err)
    hdf5_flag = hdf5_err
    call h5eset_auto_f(ON,hdf5_err)
    if (hdf5_flag < 0) then
    ! if the dataset does not exist, create it
      call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
      call h5dcreate_f(file_id,trim(string),H5T_NATIVE_DOUBLE,file_space_id, &
                      data_set_id,hdf5_err,prop_id)
    else
      call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
    endif

    call h5pclose_f(prop_id,hdf5_err)

    istart = 0
    call MPI_Exscan(local_size,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                    option%mycomm,ierr);CHKERRQ(ierr)

    start(2) = istart
    start(1) = 0

    length(2) = local_size
    length(1) = MAX_FACE_PER_CELL + 1

    stride = 1
    call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                              hdf5_err,stride,stride)
    ! write the data
    call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

    select case (var_list_type)
      case (INSTANTANEOUS_VARS)

        do i=1,local_size
          ! Num. of faces for each cell (Note: Use vec_ptr1 not vec_ptr2)
          double_array((i-1)*(MAX_FACE_PER_CELL+1)+1) = &
            vec_ptr1((i-1)*offset+1)
          ! Flowrate values for each face
          do iface = 1,MAX_FACE_PER_CELL
            double_array((i-1)*(MAX_FACE_PER_CELL+1)+iface+1) = &
            vec_ptr1((i-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface + 1)
          enddo
        enddo

      case (AVERAGED_VARS)

        do i=1,local_size
          ! Num. of faces for each cell (Note: Use vec_ptr1 not vec_ptr2)
          double_array((i-1)*(MAX_FACE_PER_CELL+1)+1) = &
            vec_ptr1((i-1)*offset+1)
          ! Divide the flowrate values by integration 'time'
          do iface = 1,MAX_FACE_PER_CELL
            double_array((i-1)*(MAX_FACE_PER_CELL+1)+iface+1) = &
            vec_ptr2((i-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface + 1)/ &
            output_option%periodic_snap_output_time_incr
          enddo
        enddo
    end select

    call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
    call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                    hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

    call h5pclose_f(prop_id,hdf5_err)

    call HDF5DatasetClose(data_set_id,option)
    call h5sclose_f(file_space_id,hdf5_err)

  enddo

  ! Free up memory
  deallocate(double_array)

end subroutine WriteHDF5FlowratesUGrid

! ************************************************************************** !

subroutine WriteHDF5FaceVelUGrid(realization_base,option,file_id, &
                                 var_list_type)
  !
  ! This routine writes velocity at cell faces for unstructured grid.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 05/25/2014
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use hdf5
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  use Variables_module
  use Connection_module
  use Coupler_module
  use HDF5_Aux_module
  use Output_Aux_module
  use Field_module

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  PetscInt :: var_list_type

  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field

  PetscInt :: iphase
  PetscInt :: offset
  PetscInt :: local_size
  PetscInt :: i
  PetscInt :: idir
  PetscInt :: istart
  PetscInt :: iface

  PetscReal, pointer :: vec_ptr1(:)
  PetscReal, pointer :: double_array(:)

  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string
  character(len=1) :: string_dir

  patch => realization_base%patch
  grid => patch%grid
  ugrid => grid%unstructured_grid
  output_option =>realization_base%output_option
  field => realization_base%field

  if (option%nphase == 1 .and. option%transport%nphase > 1) then
    option%io_buffer = 'WriteHDF5FaceVelUGrid not supported for gas &
      &transport without flow in the gas phase.'
    call PrintErrMsg(option)
  endif
  call VecGetLocalSize(field%vx_face_inst,local_size,ierr);CHKERRQ(ierr)
  local_size = local_size/(option%nphase*MAX_FACE_PER_CELL + 1)

  allocate(double_array(local_size*(MAX_FACE_PER_CELL+1)))
  double_array = 0.d0

  offset = option%nphase*MAX_FACE_PER_CELL+1

  do idir = 1,3

    select case (idir)
      case (1)
        string_dir = 'X'
        call VecGetArrayF90(field%vx_face_inst,vec_ptr1,ierr);CHKERRQ(ierr)
      case (2)
        string_dir = 'Y'
        call VecGetArrayF90(field%vy_face_inst,vec_ptr1,ierr);CHKERRQ(ierr)
      case (3)
        string_dir = 'Z'
        call VecGetArrayF90(field%vz_face_inst,vec_ptr1,ierr);CHKERRQ(ierr)
    end select

    do iphase = 1,option%nphase

      select case (iphase)
        case (LIQUID_PHASE)
          string = "Liquid " // string_dir // "-Velocity at cell face [m_per_" &
            // trim(output_option%tunit) // "]"
        case (GAS_PHASE)
          string = "Gas " // string_dir // "-Velocity at cell face [m_per_" // &
            trim(output_option%tunit) // "]"
      end select

      ! memory space which is a 1D vector
      rank_mpi = 1
      dims = 0
      dims(1) = local_size*(MAX_FACE_PER_CELL+1)
      call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

      ! file space which is a 2D block
      rank_mpi = 2
      dims = 0
      dims(2) = ugrid%nmax
      dims(1) = MAX_FACE_PER_CELL + 1
      call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

      call h5eset_auto_f(OFF,hdf5_err)
      call h5dopen_f(file_id,trim(string),data_set_id,hdf5_err)
      hdf5_flag = hdf5_err
      call h5eset_auto_f(ON,hdf5_err)
      if (hdf5_flag < 0) then
        ! if the dataset does not exist, create it
        call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
        call h5dcreate_f(file_id,trim(string),H5T_NATIVE_DOUBLE,file_space_id, &
                        data_set_id,hdf5_err,prop_id)
      else
        call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
      endif

      call h5pclose_f(prop_id,hdf5_err)

      istart = 0
      call MPI_Exscan(local_size,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                      option%mycomm,ierr);CHKERRQ(ierr)

      start(2) = istart
      start(1) = 0

      length(2) = local_size
      length(1) = MAX_FACE_PER_CELL + 1

      stride = 1
      call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                                hdf5_err,stride,stride)
      ! write the data
      call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
      call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                              hdf5_err)
#endif

      select case (var_list_type)
        case (INSTANTANEOUS_VARS)

          do i=1,local_size
            ! Num. of faces for each cell
            double_array((i-1)*(MAX_FACE_PER_CELL+1)+1) = &
              vec_ptr1((i-1)*offset+1)
            ! Flowrate values for each face
            do iface = 1,MAX_FACE_PER_CELL
              double_array((i-1)*(MAX_FACE_PER_CELL+1)+iface+1) = &
              vec_ptr1((i-1)*offset + (iphase-1)*MAX_FACE_PER_CELL + iface + 1)* &
              realization_base%output_option%tconv
            enddo
          enddo

        case (AVERAGED_VARS)

      end select

      call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
      call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)

      call h5pclose_f(prop_id,hdf5_err)

      call HDF5DatasetClose(data_set_id,option)
      call h5sclose_f(file_space_id,hdf5_err)

    enddo

    select case (idir)
      case (1)
        call VecRestoreArrayF90(field%vx_face_inst,vec_ptr1, &
                                ierr);CHKERRQ(ierr)
      case (2)
        call VecRestoreArrayF90(field%vy_face_inst,vec_ptr1, &
                                ierr);CHKERRQ(ierr)
      case (3)
        call VecRestoreArrayF90(field%vz_face_inst,vec_ptr1, &
                                ierr);CHKERRQ(ierr)
    end select

  enddo

  ! Free up memory
  deallocate(double_array)

end subroutine WriteHDF5FaceVelUGrid

! ************************************************************************** !

subroutine OutputHDF5Provenance(option, output_option, file_id)
  !
  ! write pflotran and petsc provenance information including a copy
  ! of the inputfile
  !

  use Option_module, only : option_type
  use Output_Aux_module, only : output_option_type
  use PFLOTRAN_Provenance_module, only : provenance_max_str_len

  use hdf5
  use HDF5_Aux_module

  implicit none

  type(option_type), intent(in) :: option
  type(output_option_type), intent(in) :: output_option
  integer(HID_T), intent(in) :: file_id

  character(len=32) :: name
  integer(HID_T) :: provenance_id, string_type
  PetscMPIInt :: hdf5_err
  integer(SIZE_T) :: size_t_int

  ! create the provenance group
  name = "Provenance"
  call HDF5GroupCreate(file_id,name,provenance_id,option)

  ! create fixed length string datatype
  call h5tcopy_f(H5T_FORTRAN_S1, string_type, hdf5_err)
  size_t_int = provenance_max_str_len
  call h5tset_size_f(string_type, size_t_int, hdf5_err)

  call OutputHDF5Provenance_PFLOTRAN(option, provenance_id, string_type)
  call OutputHDF5Provenance_PETSc(provenance_id, string_type, option)

  ! close the provenance group
  call h5tclose_f(string_type, hdf5_err)
  call HDF5GroupClose(provenance_id,option)

end subroutine OutputHDF5Provenance

! ************************************************************************** !

subroutine OutputHDF5Provenance_PFLOTRAN(option, provenance_id, string_type)
  !
  ! write the pflotran provenance data as attributes (small) or
  ! datasets (big details)
  !

  use Option_module, only : option_type
  use PFLOTRAN_Provenance_module

  use hdf5
  use HDF5_Aux_module

  implicit none

  type(option_type), intent(in) :: option
  integer(HID_T), intent(in) :: provenance_id
  integer(HID_T), intent(in) :: string_type

  character(len=32) :: name
  integer(HID_T) :: pflotran_id

  ! Create the pflotran group under provenance
  name = "PFLOTRAN"
  call HDF5GroupCreate(provenance_id, name, pflotran_id, option)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "pflotran_compile_date_time", &
                                    ONE_INTEGER, pflotran_compile_date_time)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "pflotran_compile_user", &
                                    ONE_INTEGER, pflotran_compile_user)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "pflotran_compile_hostname", &
                                    ONE_INTEGER, pflotran_compile_hostname)

  call OutputHDF5AttributeStringArray(pflotran_id, string_type, &
                                      "pflotran_status", &
                                      ONE_INTEGER, pflotran_status)

  call OutputHDF5AttributeStringArray(pflotran_id, string_type, &
                                      "pflotran_changeset", &
                                      ONE_INTEGER, pflotran_changeset)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "detail_pflotran_fflags", &
                                    detail_pflotran_fflags_len, &
                                    detail_pflotran_fflags)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "detail_pflotran_status", &
                                    detail_pflotran_status_len, &
                                    detail_pflotran_status)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "detail_pflotran_parent", &
                                    detail_pflotran_parent_len, &
                                    detail_pflotran_parent)

  ! FIXME(bja, 2013-11-25): break gcc when diffs are present
  call OutputHDF5DatasetStringArray(pflotran_id, string_type, &
                                    "detail_pflotran_diff", &
                                    detail_pflotran_diff_len, &
                                    detail_pflotran_diff)

  call OutputHDF5Provenance_input(option, pflotran_id)

  ! close pflotran group
  call HDF5GroupClose(pflotran_id,option)

end subroutine OutputHDF5Provenance_PFLOTRAN

! ************************************************************************** !

subroutine OutputHDF5Provenance_input(option, pflotran_id)
  !
  ! open the pflotran input file, figure out how long it is, read it
  ! into a buffer, then write the buffer as a pflotran provenance
  ! group dataset.
  !
  use hdf5
  use Input_Aux_module, only : input_type, InputCreate, InputDestroy, &
       InputGetLineCount, InputReadToBuffer
  use Option_module, only : option_type
  use PFLOTRAN_Constants_module, only : IN_UNIT, MAXSTRINGLENGTH

  implicit none

  type(option_type), intent(in) :: option
  integer(HID_T), intent(in) :: pflotran_id

  integer(HID_T) :: input_string_type
  type(input_type), pointer :: input
  PetscInt :: input_line_count
  character(len=MAXSTRINGLENGTH), allocatable :: input_buffer(:)
  PetscMPIInt :: hdf5_err
  integer(SIZE_T) :: size_t_int

  input => InputCreate(IN_UNIT, option%input_filename, option)
  input_line_count = InputGetLineCount(input,option)
  allocate(input_buffer(input_line_count))
  call InputReadToBuffer(input, input_buffer, option)
  call h5tcopy_f(H5T_FORTRAN_S1, input_string_type, hdf5_err)
  size_t_int = MAXSTRINGLENGTH
  call h5tset_size_f(input_string_type, size_t_int, hdf5_err)
  call OutputHDF5DatasetStringArray(pflotran_id, input_string_type, &
                                    "pflotran_input_file", &
                                    input_line_count, input_buffer)
  call h5tclose_f(input_string_type, hdf5_err)
  deallocate(input_buffer)
  call InputDestroy(input)

end subroutine OutputHDF5Provenance_input

! ************************************************************************** !

subroutine OutputHDF5Provenance_PETSc(provenance_id, string_type, option)
  !
  ! write the petsc provenance data as attributes (small) or datasets
  ! (big details)
  !

  use PFLOTRAN_Provenance_module
  use hdf5
  use HDF5_Aux_module
  use Option_module

  implicit none

  integer(HID_T), intent(in) :: provenance_id
  integer(HID_T), intent(in) :: string_type
  type(option_type), intent(in) :: option

  character(len=32) :: name
  integer(HID_T) :: petsc_id

  ! create the petsc group under provenance
  name = "PETSc"
  call HDF5GroupCreate(provenance_id, name, petsc_id, option)

  call OutputHDF5AttributeStringArray(petsc_id, string_type, "petsc_status", &
                                      ONE_INTEGER, petsc_status)

  call OutputHDF5AttributeStringArray(petsc_id, string_type, &
                                      "petsc_changeset", &
                                      ONE_INTEGER, petsc_changeset)

  call OutputHDF5DatasetStringArray(petsc_id, string_type, &
                                    "detail_petsc_status", &
                                    detail_petsc_status_len, &
                                    detail_petsc_status)

  call OutputHDF5DatasetStringArray(petsc_id, string_type, &
                                    "detail_petsc_parent", &
                                    detail_petsc_parent_len, &
                                    detail_petsc_parent)

  call OutputHDF5DatasetStringArray(petsc_id, string_type, &
                                    "detail_petsc_config", &
                                    detail_petsc_config_len, &
                                    detail_petsc_config)

  ! close the petsc group
  call HDF5GroupClose(petsc_id,option)

end subroutine OutputHDF5Provenance_PETSc

! ************************************************************************** !

subroutine OutputHDF5AttributeStringArray(parent_id, type, name, length, data)
  ! create the dataspaces and attributes consisting of an array of
  ! strings, then write the data and cleanup

  use hdf5

  implicit none

  integer(HID_T), intent(in) ::  parent_id, type
  character(len=*), intent(in) :: name
  PetscInt, intent(in) :: length
  character(len=*), intent(in) :: data(length)

  integer(HID_T) :: dataspace_id, attribute_id
  integer(HSIZE_T), dimension(1:1) :: dims
  PetscMPIInt :: hdf5_err

  dims = length
  call h5screate_simple_f(1, dims, dataspace_id, hdf5_err)
  call h5eset_auto_f(OFF,hdf5_err)
  call h5aopen_f(parent_id, name, attribute_id, hdf5_err)
  if (hdf5_err /= 0) then
    call h5acreate_f(parent_id, name, type, dataspace_id, &
                     attribute_id, hdf5_err)
  endif
  call h5eset_auto_f(ON,hdf5_err)
  call h5awrite_f(attribute_id, type, data, dims, hdf5_err)
  call h5aclose_f(attribute_id, hdf5_err)
  call h5sclose_f(dataspace_id, hdf5_err)

end subroutine OutputHDF5AttributeStringArray

! ************************************************************************** !

subroutine OutputHDF5DatasetStringArray(parent_id, type, name, length, data)
  ! create the dataspaces and dataset consisting of an array of
  ! strings, then write the data and cleanup

  use hdf5

  implicit none

  integer(HID_T), intent(in) ::  parent_id, type
  character(len=*), intent(in) :: name
  PetscInt, intent(in) :: length
  character(len=*), intent(in) :: data(length)

  integer(HID_T) :: dataspace_id, attribute_id
  integer(HSIZE_T), dimension(1:1) :: dims
  PetscMPIInt :: hdf5_err

  dims = length
  call h5screate_simple_f(1, dims, dataspace_id, hdf5_err)
  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(parent_id, name, attribute_id, hdf5_err)
  if (hdf5_err /= 0) then
    call h5dcreate_f(parent_id, name, type, dataspace_id, &
                     attribute_id, hdf5_err)
  endif
  call h5eset_auto_f(ON,hdf5_err)
  call h5dwrite_f(attribute_id, type, data, dims, hdf5_err)
  call h5dclose_f(attribute_id, hdf5_err)
  call h5sclose_f(dataspace_id, hdf5_err)

end subroutine OutputHDF5DatasetStringArray

! ************************************************************************** !

subroutine OutputHDF5WriteSnapShotAtts(parent_id,option)
  !
  ! Writes attributes associated with a snapshot time in the output file.
  !
  ! Author: Glenn Hammond
  ! Date: 07/31/19
  !
  use hdf5
  use Option_module

  implicit none

  integer(HID_T) :: parent_id
  type(option_type) :: option

  integer(HID_T) :: attribute_id
  integer(HID_T) :: dataspace_id
  character(len=MAXSTRINGLENGTH) :: string
  integer(HSIZE_T) :: dims(1)
  PetscMPIInt :: hdf5_err

  dims = 1
  call h5screate_simple_f(1,dims,dataspace_id,hdf5_err)
  string = 'Time (s)'
  call h5eset_auto_f(OFF,hdf5_err)
  call h5aopen_f(parent_id, string, attribute_id, hdf5_err)
  if (hdf5_err /= 0) then
    call h5acreate_f(parent_id,string,H5T_NATIVE_DOUBLE,dataspace_id, &
                     attribute_id,hdf5_err)
  endif
  call h5eset_auto_f(ON,hdf5_err)
  call h5awrite_f(attribute_id,H5T_NATIVE_DOUBLE,option%time,dims,hdf5_err)
  call h5aclose_f(attribute_id, hdf5_err)
  call h5sclose_f(dataspace_id, hdf5_err)

end subroutine OutputHDF5WriteSnapShotAtts

! ************************************************************************** !

subroutine OutputH5OpenFile(option, h5obj, filename, file_id)
  !
  ! Opens an HDF5 file, creating it if the first time
  !
  ! Author: Glenn Hammond
  ! Date: 10/19/19
  !
  use hdf5
  use HDF5_Aux_module
  use Option_module

  implicit none

  type(option_type) :: option
  type(output_h5_type) :: h5obj
  character(len=MAXSTRINGLENGTH) :: filename
  integer(HID_T), intent(out) :: file_id

  if (.not.h5obj%first_write) then
    call HDF5FileTryOpen(filename,file_id,h5obj%first_write,option%comm)
  endif
  if (h5obj%first_write) then
    call HDF5FileOpen(filename,file_id,PETSC_TRUE,option)
  endif

  if (h5obj%first_write) then
    option%io_buffer = ' --> creating hdf5 output file: ' // trim(filename)
  else
    option%io_buffer = ' --> appending to hdf5 output file: ' // trim(filename)
  endif
  call PrintMsg(option)

end subroutine OutputH5OpenFile

! ************************************************************************** !

subroutine OutputH5CloseFile(option, h5file, file_id)
  !
  ! Closes an HDF5 file
  !
  ! Author: Glenn Hammond
  ! Date: 10/19/19
  !
  use hdf5
  use HDF5_Aux_module
  use Option_module

  implicit none

  type(option_type), intent(in), pointer :: option
  type(output_h5_type) :: h5file
  integer(HID_T), intent(in) :: file_id

  call HDF5FileClose(file_id,option)
  h5file%first_write = PETSC_FALSE

end subroutine OutputH5CloseFile

! ************************************************************************** !

subroutine OutputXMFOpenFile(option, filename, fid)
  !
  ! Opens an XMF file
  !
  ! Author: Glenn Hammond
  ! Date: 10/19/19
  !
  use Option_module

  implicit none

  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid

  if (OptionIsIORank(option)) then
    option%io_buffer = ' --> write xmf output file: ' // trim(filename)
    call PrintMsg(option)
    open(unit=fid,file=filename,action="write")
  endif

end subroutine OutputXMFOpenFile

! ************************************************************************** !

subroutine OutputH5OpenGroup(option, group_name, file_id, grp_id)
  !
  ! Opens an HDF5 group, creating it if it does not exist
  !
  ! Author: Glenn Hammond
  ! Date: 10/19/19
  !
  use hdf5
  use HDF5_Aux_module
  use Option_module

  implicit none

  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: group_name
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id

  call HDF5GroupOpenOrCreate(file_id,group_name,grp_id,option)

end subroutine OutputH5OpenGroup

! ************************************************************************** !

subroutine OutputH5CloseGroup(option,grp_id)
  !
  ! Closes an HDF5 group
  !
  ! Author: Glenn Hammond
  ! Date: 10/19/19
  !
  use hdf5
  use HDF5_Aux_module
  use Option_module

  implicit none

  type(option_type) :: option
  integer(HID_T) :: grp_id

  call HDF5GroupClose(grp_id,option)

end subroutine OutputH5CloseGroup

end module Output_HDF5_module
