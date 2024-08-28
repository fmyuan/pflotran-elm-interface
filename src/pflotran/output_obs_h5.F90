module Output_Obs_H5_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use hdf5
  use Logging_module
  use PFLOTRAN_Constants_module
  use Output_Aux_module
  use Output_Common_module
  use Region_module

  implicit none

  private

  PetscMPIInt, private, parameter :: ON=1, OFF=0

  ! flag signifying the first time output routine is called
  PetscBool :: observation_hdf5_first

  ! communicator information
  PetscMPIInt :: obs_h5_comm
  PetscMPIInt :: comm_size
  PetscMPIInt :: comm_rank

  character(len=MAXSTRINGLENGTH) :: h5_filename

  ! attribute information
  character(len=MAXWORDLENGTH) :: sim_time
  character(len=MAXWORDLENGTH) :: version

  ! information used to describe data
  PetscInt :: num_regions  ! Total number of defined regions

  ! arrays of dimension num_regions used for metadata
  PetscInt, pointer :: obs_region_ids(:)
  PetscInt, pointer :: obs_region_total_numcells(:)
  character(len=MAXWORDLENGTH), pointer :: region_names(:)

  ! arrays of dimension comm_size x num_regions
  PetscInt, pointer :: obs_region_numcells(:,:) ! extent in writes
  PetscInt, pointer :: obs_region_offset(:,:)   ! start in writes

  ! generic interface for writing data
  interface WriteH5RegionDataset
    module procedure WriteH5RegionRealDataset, &
                     WriteH5RegionIntegerDataset
  end interface WriteH5RegionDataset


  public :: OutputObsH5,      &
            OutputObsH5Init

contains

  subroutine OutputObsH5Init(num_steps)
    !
    ! Initializes module variables for HDF5 output
    ! Gets date/time information for simulation as
    ! well as version information and stores these
    ! to be output as attributes.
    !
    ! Author: R. McKeown
    ! Date: 01/26/2022
    !

    implicit none

    PetscInt :: num_steps

    character(len=8) :: date_word
    character(len=10) :: time_word
    character(len=5) :: zone_word

    if (num_steps == 0) then
      observation_hdf5_first = PETSC_TRUE

      call date_and_time(date_word,time_word,zone_word)
      write(sim_time,'(a)') date_word(5:6) // '/'         &
                            // date_word(7:8) // '/'      &
                            // date_word(1:4) // ' '      &
                            // time_word(1:2) // ':'      &
                            // time_word(3:4) // ' ('     &
                            // zone_word(1:3) // ':'      &
                            // zone_word(4:5) // ' UTC)'

    ! retrieve version
    version = GetVersion()

    else
      observation_hdf5_first = PETSC_FALSE
    endif

  end subroutine OutputObsH5Init


  subroutine OutputObsH5(realization_base)
    !
    ! Output observation variables by region to HDF5 file.
    !
    ! Author: R. McKeown
    ! Date: 01/29/2022
    !

    use Realization_Base_class
    use Option_module

    implicit none

    class(realization_base_type) :: realization_base
    type(option_type), pointer   :: option

    PetscMPIInt :: mpi_err

    option => realization_base%option

    ! first time step only
    if (observation_hdf5_first) then

      ! duplicate and save the usual communicator
      call MPI_Comm_dup(option%mycomm,obs_h5_comm,mpi_err);CHKERRQ(mpi_err)
      call MPI_Comm_size(obs_h5_comm,comm_size,mpi_err);CHKERRQ(mpi_err)
      call MPI_Comm_rank(obs_h5_comm,comm_rank,mpi_err);CHKERRQ(mpi_err)

      ! create empty hdf5 file for output
      call CreateH5ObservationFile(realization_base)

      ! get information required for metadata
      call GetH5ObservationRegionMetadata(realization_base)

      ! set up file structure and write time invariant data
      call WriteH5ObservationRegionDomain(realization_base)

      observation_hdf5_first = PETSC_FALSE
    endif

    ! for all time steps
    call WriteH5ObservationRegionVariables(realization_base)

  end subroutine OutputObsH5


  subroutine CreateH5ObservationFile(realization_base)
    !
    ! Create an HDF5 file for observation variables by region.
    !
    ! Author: R. McKeown
    ! Date: 01/26/2022
    !

    use Realization_Base_class, only : realization_base_type
    use Option_module

    implicit none

    class(realization_base_type) :: realization_base
    type(option_type), pointer   :: option

    PetscMPIInt :: h5_err

    ! HDF5 handles
    integer(HID_T) :: file_id
    integer(HID_T) :: fapl_id


    option => realization_base%option

    h5_filename = trim(option%global_prefix) // '-obs-region.h5'

    ! create an hdf5 file access property list
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, h5_err)
#ifndef SERIAL_HDF5
    call h5pset_fapl_mpio_f(fapl_id, obs_h5_comm, MPI_INFO_NULL, h5_err)
#endif

    call h5fcreate_f(h5_filename, H5F_ACC_TRUNC_F, file_id, h5_err, &
                     access_prp=fapl_id)

    ! add version and date attributes
    call AddVersionAndDateAttributes(file_id)

    call h5pclose_f(fapl_id, h5_err)

    call h5fclose_f(file_id, h5_err)

  end subroutine CreateH5ObservationFile


  subroutine AddVersionAndDateAttributes(file_id)
    !
    ! Output the verion and build date of PFLOTRAN as
    ! attributes in the HDF5 observation output file.
    !
    !
    ! Author: R. McKeown
    ! Date: 10/28/2022
    !

    implicit none

    PetscMPIInt :: h5_err

    character(len=16), parameter :: vname = "version"
    character(len=16), parameter :: dname = "simulation date"

    integer        :: arank = 1   ! Attribute rank
    integer(HID_T) :: file_id     ! HDF5 file handle
    integer(HID_T) :: atype_id    ! Attribute type
    integer(HID_T) :: attr_id     ! Attribute identifier
    integer(HID_T) :: attr_space  ! Attribute space
    integer(HSIZE_T) :: attrlen   ! Length of attribute string
    integer(HSIZE_T), dimension(1) :: adims = (/1/)   ! Attribute dimension
    integer(HSIZE_T), dimension(1) :: data_dims


    attrlen = MAXWORDLENGTH

    ! create scalar data space for the version attribute
    call h5screate_simple_f(arank, adims, attr_space, h5_err)

    ! create datatype for the version attribute
    call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5_err)
    call h5tset_size_f(atype_id, attrlen, h5_err)

    ! create the version attribute for observation output file
    call h5acreate_f(file_id, vname, atype_id, attr_space, attr_id, h5_err)

    ! write the attribute data
    data_dims(1) = 2
    call h5awrite_f(attr_id, atype_id, version, data_dims, h5_err)

    ! close the attribute, type, and data space
    call h5aclose_f(attr_id, h5_err)
    call h5tclose_f(atype_id, h5_err)
    call h5sclose_f(attr_space, h5_err)

    ! create scalar data space for the simulation time attribute
    call h5screate_simple_f(arank, adims, attr_space, h5_err)

    ! create datatype for the simulation date attribute
    call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5_err)
    call h5tset_size_f(atype_id, attrlen, h5_err)

    ! create the simulation date attribute for observation output file
    call h5acreate_f(file_id, dname, atype_id, attr_space, attr_id, h5_err)

    ! write the attribute data
    data_dims(1) = 2
    call h5awrite_f(attr_id, atype_id, sim_time, data_dims, h5_err)

    ! close the attribute, type, and data space
    call h5aclose_f(attr_id, h5_err)
    call h5tclose_f(atype_id, h5_err)
    call h5sclose_f(attr_space, h5_err)


  end subroutine AddVersionAndDateAttributes

  subroutine GetH5ObservationRegionMetadata(realization_base)
    !
    ! For requested observation regions, collect and disseminate to
    ! all MPI processes the region ID, the total number of cells
    ! in the region, and the offset  for writing data to the
    ! HDF5 dataset. This is the information required to define
    ! the structure of the HDF5 observation output file.
    !
    ! This information is stored in vectors of length comm_size
    ! and arrays of size comm_size x num_regions, where comm_size
    ! is the number of processes associated with the communicator
    ! and num_regions is the total number of regions defined in
    ! the input file (not the number requested by observation).
    !
    !
    ! Author: R. McKeown
    ! Date: 01/26/2022
    !

    use Realization_Base_class, only : realization_base_type
    use Patch_module
    use Observation_module
    use Region_module

    implicit none

    class(realization_base_type)      :: realization_base
    type(patch_type), pointer         :: patch
    type(observation_type), pointer   :: observation
    type(region_type), pointer        :: region
    type(region_list_type), pointer   :: region_list

    PetscInt     :: i_obs
    PetscInt     :: i_rank, i_region
    PetscInt     :: num_obs
    PetscInt     :: arr_size
    PetscInt     :: tmp_region_id
    PetscMPIInt  :: mpi_err


    patch => realization_base%patch
    region_list => patch%region_list

    ! total number of observations and regions
    num_regions = region_list%num_regions
    num_obs = patch%observation_list%num_observations

    call MPI_Allreduce(MPI_IN_PLACE,num_obs,comm_size,MPI_INTEGER,MPI_MAX, &
                       obs_h5_comm,mpi_err);CHKERRQ(mpi_err)

    ! total number of elements in 2D arrays
    !    - required for MPI_Allreduce
    arr_size = comm_size * num_regions

    allocate(obs_region_ids(num_regions))
    allocate(region_names(num_regions))
    allocate(obs_region_total_numcells(num_regions))
    allocate(obs_region_numcells(comm_size, num_regions))
    allocate(obs_region_offset(comm_size, num_regions))

    ! initialize module arrays
    obs_region_ids(:) = 0
    obs_region_total_numcells(:) = 0
    region_names(:) = ' '
    obs_region_numcells(:,:) = 0
    obs_region_offset(:,:) = 0

    observation => patch%observation_list%first

    ! is observation local
    do i_obs = 1, num_obs

      tmp_region_id = 0

      if (associated(observation)) then
        if (associated(observation%region)) then

          tmp_region_id = observation%region%id

          ! fill some array elements
          obs_region_ids(tmp_region_id) = tmp_region_id
          obs_region_numcells(comm_rank+1, tmp_region_id) =   &
                              observation%region%num_cells
        endif

      endif

      if (associated(observation)) observation => observation%next

    enddo

    ! loop through possible observation regions
    i_region = 0
    region => region_list%first
    do
      if (.not. associated(region) ) exit

      i_region = i_region + 1
      region_names(i_region) = region%name

      region => region%next
    enddo

    ! reduce so all processes have the same information
    call MPI_Allreduce(MPI_IN_PLACE,obs_region_numcells,arr_size,MPI_INTEGER, &
                       MPI_SUM,obs_h5_comm,mpi_err);CHKERRQ(mpi_err)

    call MPI_Allreduce(MPI_IN_PLACE,obs_region_ids,num_regions,MPI_INTEGER, &
                       MPI_MAX,obs_h5_comm,mpi_err);CHKERRQ(mpi_err)


    ! after MPI_Allreduce, every MPI process that is part
    ! of the communicator should see the same values for
    ! obs_region_numcells and obs_region_ids

    ! partial sums for offsets - HDF5 offsets are zero based
    do i_region = 1, num_regions

      do i_rank = 2, comm_size
        obs_region_offset(i_rank, i_region) =  &
            SUM(obs_region_numcells(1:i_rank-1, i_region))
      enddo

      ! full sum for total number of cells for the region
      obs_region_total_numcells = SUM(obs_region_numcells, dim=1)
    enddo

  end subroutine GetH5ObservationRegionMetadata


  subroutine WriteH5ObservationRegionDomain(realization_base)
    !
    ! Add region domain group and time invariant data
    ! to HDF5 observation file.
    !
    ! Author: R. McKeown
    ! Date: 02/07/2022
    !

    use Realization_Base_class, only : realization_base_type

    implicit none

    class(realization_base_type) :: realization_base

    character(len=MAXWORDLENGTH) :: group_name
    character(len=MAXWORDLENGTH) :: subgroup_name
    character(len=MAXWORDLENGTH) :: dset_name

    PetscInt :: ireg
    PetscInt :: region_id
    PetscMPIInt :: h5_err

    ! HDF5 handles
    integer(HID_T) :: file_id
    integer(HID_T) :: fapl_id
    integer(HID_T) :: group_id
    integer(HID_T) :: subgroup_id
    integer(HID_T) :: data_type

    ! parallel hdf5 property list
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, h5_err)
#ifndef SERIAL_HDF5
    call h5pset_fapl_mpio_f(fapl_id, obs_h5_comm, MPI_INFO_NULL, h5_err)
#endif

    ! open HDF5 observation file
    call h5fopen_f(h5_filename, H5F_ACC_RDWR_F, file_id, h5_err, &
                   fapl_id)

    ! close property list
    call h5pclose_f(fapl_id, h5_err)

    do ireg = 1, num_regions

      if (obs_region_ids(ireg) > 0) then
        ! get region name
        group_name = region_names(ireg)

        !call h5gcreate_f(file_id, group_name, group_id, h5_err)
        call CreateOrOpenH5Group(file_id, group_id, group_name)

        ! create Domain group by region for time invariant information
        subgroup_name = 'Domain'
        !call h5gcreate_f(group_id, subgroup_name, subgroup_id, h5_err)
        call CreateOrOpenH5Group(group_id, subgroup_id, subgroup_name)

        ! create datasets for grid info
        region_id = obs_region_ids(ireg)
        data_type = H5T_NATIVE_DOUBLE

        dset_name = "XC"
        call CreateH5RegionDataset(subgroup_id, dset_name,   &
                                     data_type, region_id)

        dset_name = "YC"
        call CreateH5RegionDataset(subgroup_id, dset_name,   &
                                     data_type, region_id)

        dset_name = "ZC"
        call CreateH5RegionDataset(subgroup_id, dset_name,   &
                                     data_type, region_id)

        dset_name = "Volume"
        call CreateH5RegionDataset(subgroup_id, dset_name,   &
                                     data_type, region_id)

        call h5gclose_f(subgroup_id, h5_err)
        call h5gclose_f(group_id, h5_err)

      endif
    enddo

    call WriteH5RegionGridInfo(realization_base, file_id)

    call h5fclose_f(file_id, h5_err)


  end subroutine WriteH5ObservationRegionDomain


  subroutine WriteH5RegionGridInfo(realization_base, file_id)
    !
    ! Add grid info to an HDF5 observation region file.
    !
    ! Author: R. McKeown
    ! Date: 01/29/2022
    !

    use Realization_Base_class, only : realization_base_type
    use Material_Aux_module, only : material_auxvar_type
    use Patch_module
    use Observation_module
    use Grid_module
    use Region_module

    implicit none

    class(realization_base_type)         :: realization_base
    type(material_auxvar_type), pointer  :: material_auxvars(:)
    type(patch_type), pointer            :: patch
    type(observation_type), pointer      :: observation
    type(grid_type), pointer             :: grid
    type(region_type), pointer           :: region

    integer(HID_T)               :: file_id
    integer(HID_T)               :: group_id

    PetscInt   :: icell
    PetscInt   :: local_id, ghosted_id
    PetscMPIInt :: h5_err
    PetscReal, pointer :: x(:), y(:), z(:)
    PetscReal, pointer :: volume(:)

    character(len=MAXSTRINGLENGTH) :: group_name
    character(len=MAXWORDLENGTH) :: dset_name


    patch => realization_base%patch
    material_auxvars => patch%aux%Material%auxvars
    grid => patch%grid

    ! loop over local observations
    observation => patch%observation_list%first

    do
      if (.not. associated(observation)) exit

      region => observation%region

      allocate(x(region%num_cells))
      allocate(y(region%num_cells))
      allocate(z(region%num_cells))
      allocate(volume(region%num_cells))

      ! collect x, y, z, and volume
      do icell = 1, region%num_cells
        local_id = region%cell_ids(icell)
        ghosted_id = grid%nL2G(local_id)

        x(icell) = grid%x(ghosted_id)
        y(icell) = grid%y(ghosted_id)
        z(icell) = grid%z(ghosted_id)
        volume(icell) = material_auxvars(ghosted_id)%volume
      enddo

      ! open correct region/Domain group
      group_name = trim(region%name) // '/Domain'
      call h5gopen_f(file_id, group_name, group_id, h5_err)

      ! open previously generated data space and write data
      dset_name = "XC"
      call WriteH5RegionDataset(group_id, x, region%id,   &
                                  dset_name)

      ! open previously generated data space and write data
      dset_name = "YC"
      call WriteH5RegionDataset(group_id, y, region%id,   &
                                  dset_name)

      ! open previously generated data space and write data
      dset_name = "ZC"
      call WriteH5RegionDataset(group_id, z, region%id,   &
                                  dset_name)

      ! open previously generated data space and write data
      dset_name = "Volume"
      call WriteH5RegionDataset(group_id, volume, region%id,   &
                                  dset_name)

      deallocate(x)
      deallocate(y)
      deallocate(z)
      deallocate(volume)

      call h5gclose_f(group_id, h5_err)

      observation => observation%next
    enddo

  end subroutine WriteH5RegionGridInfo


  subroutine WriteH5ObservationRegionVariables(realization_base)
    !
    ! Write observation variable values at current time to an HDF5
    ! observation region file.
    !
    ! Author: R. McKeown
    ! Date: 02/10/2022
    !

    use Realization_Base_class, only : realization_base_type
    use Option_module
    use Grid_module
    use Patch_module
    use Observation_module
    use Region_module

    implicit none

    class(realization_base_type)        :: realization_base
    type(option_type), pointer          :: option
    type(patch_type), pointer           :: patch
    type(observation_type), pointer     :: observation
    type(grid_type), pointer            :: grid
    type(region_type), pointer          :: region
    type(output_variable_type), pointer :: cur_variable
    type(output_option_type), pointer   :: output_option

    integer(HID_T) :: file_id
    integer(HID_T) :: fapl_id
    integer(HID_T) :: tgroup_id
    integer(HID_T) :: rgroup_id

    PetscMPIInt :: h5_err
    PetscInt    :: icell
    PetscInt    :: local_id, ghosted_id
    PetscInt    :: region_id
    PetscReal   :: tmp_real

    PetscInt, pointer  :: idata(:)
    PetscReal, pointer :: rdata(:)

    character(len=MAXWORDLENGTH) :: v_name
    character(len=MAXSTRINGLENGTH) :: t_name
    character(len=MAXSTRINGLENGTH) :: r_name

    option => realization_base%option
    patch => realization_base%patch
    grid => patch%grid
    output_option => realization_base%output_option

    ! parallel hdf5 property list
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, h5_err)
#ifndef SERIAL_HDF5
    call h5pset_fapl_mpio_f(fapl_id, obs_h5_comm, MPI_INFO_NULL, h5_err)
#endif

    ! open HDF5 observation file
    call h5fopen_f(h5_filename, H5F_ACC_RDWR_F, file_id, h5_err, &
                   fapl_id)

    ! close property list
    call h5pclose_f(fapl_id, h5_err)

    ! create time group collectively
    write(t_name,'(''Time:'',es13.5,x,a1)') &
          option%time/output_option%tconv, output_option%tunit

    ! create time group and variable data spaces
    call WriteH5ObservationRegionVariableMetadata   &
                (realization_base, file_id, t_name)

    ! loop over local observations - independent processes
    observation => patch%observation_list%first

    do
      if (.not. associated(observation)) exit

      region => observation%region
      r_name = trim(region%name)

      ! open region group
      call h5gopen_f(file_id, r_name, rgroup_id, h5_err)

      ! open time group
      call h5gopen_f(rgroup_id, t_name, tgroup_id, h5_err)

      allocate(idata(region%num_cells))
      allocate(rdata(region%num_cells))

      idata = 0
      rdata = 0.0

      ! loop over observation variables and write to file
      cur_variable => output_option%output_obs_variable_list%first

      do
        if (.not. associated(cur_variable)) exit
        if (cur_variable%plot_only) then
          cur_variable => cur_variable%next
          cycle
        endif

        idata = 0
        rdata = 0.0
        v_name = ' '

        region_id = observation%region%id

        do icell = 1, region%num_cells

          local_id = region%cell_ids(icell)
          ghosted_id = grid%nL2G(local_id)

          tmp_real = OutputGetVariableAtCell(realization_base,   &
                                              ghosted_id, &
                                              cur_variable)

          if (cur_variable%iformat == 0) then ! real
            rdata(icell) = tmp_real
          else                                ! integer
            idata(icell) = int(tmp_real + 1.d-5)
          endif

        enddo

        ! get variable name
        v_name = OutputVariableGetName(cur_variable)

        if (cur_variable%iformat == 0) then
          call WriteH5RegionDataset(tgroup_id, rdata,    &
                                      region_id, v_name)
        else
          call WriteH5RegionDataset(tgroup_id, idata,    &
                                      region_id, v_name)
        endif

        cur_variable => cur_variable%next
      enddo

      deallocate(idata)
      deallocate(rdata)

      ! close time and region groups
      call h5gclose_f(tgroup_id, h5_err)
      call h5gclose_f(rgroup_id, h5_err)

      observation => observation%next
    enddo

    ! close file
    call h5fclose_f(file_id, h5_err)

  end subroutine WriteH5ObservationRegionVariables


  subroutine WriteH5ObservationRegionVariableMetadata   &
                    (realization_base, file_id, t_name)
    !
    ! Add time group and variable datasets to region groups
    !
    ! Author: R. McKeown
    ! Date: 02/19/2022
    !

    use Realization_Base_class, only : realization_base_type
    use Output_Aux_module
    use Patch_module
    use Region_module

    implicit none

    type(output_variable_type), pointer :: cur_variable
    type(output_option_type), pointer   :: output_option

    class(realization_base_type) :: realization_base

    character(len=MAXWORDLENGTH) :: dset_name
    character(len=MAXSTRINGLENGTH) :: t_name
    character(len=MAXSTRINGLENGTH) :: r_name

    PetscInt    :: ireg
    PetscInt    :: region_id
    PetscMPIInt :: h5_err

    ! HDF5 handles
    integer(HID_T) :: file_id
    integer(HID_T) :: rgroup_id
    integer(HID_T) :: tgroup_id
    integer(HID_T) :: data_type

    output_option => realization_base%output_option

    ! loop over regions for collective metadata calls
    do ireg = 1, num_regions
      if (obs_region_ids(ireg) > 0) then
        ! create groups and datasets for grid info
        region_id = obs_region_ids(ireg)

        ! open correct region group and create time group
        r_name = trim(region_names(ireg))
        call h5gopen_f(file_id, r_name, rgroup_id, h5_err)

        ! create time group by region
        ! call h5gcreate_f(rgroup_id, t_name, tgroup_id, h5_err)
        call CreateOrOpenH5Group(rgroup_id, tgroup_id, t_name)

        ! loop over observation variables and write to file
        cur_variable => output_option%output_obs_variable_list%first

        do
          if (.not.associated(cur_variable)) exit
          if (cur_variable%plot_only) then
            cur_variable => cur_variable%next
            cycle
          endif

          if (cur_variable%iformat == 0) then
            data_type = H5T_NATIVE_DOUBLE
          else
            data_type = H5T_NATIVE_INTEGER
          endif

          dset_name = OutputVariableGetName(cur_variable)
          call CreateH5RegionDataset(tgroup_id, dset_name,   &
                                       data_type, region_id)

          cur_variable => cur_variable%next
        enddo

        ! close groups
        call h5gclose_f(tgroup_id, h5_err)
        call h5gclose_f(rgroup_id, h5_err)

      endif
    enddo

  end subroutine WriteH5ObservationRegionVariableMetadata


  subroutine CreateH5RegionDataset(group_id, dset_name, data_type, region_id)
    !
    ! Create a dataset in the observation region file.
    !
    ! Author: R. McKeown
    ! Date: 02/16/2022
    !

    implicit none

    character(len=MAXWORDLENGTH) :: dset_name

    PetscInt :: region_id
    PetscInt :: ds_rank = 1  ! data is always a 1D array

    PetscMPIInt :: h5_err

    integer(HID_T) :: data_type
    integer(HID_T) :: group_id
    integer(HID_T) :: space_id
    integer(HID_T) :: dset_id
    integer(HSIZE_T) :: dset_size(1)

    ! create data space entire dataset
    dset_size(1) = obs_region_total_numcells(region_id)
    call h5screate_simple_f(ds_rank, dset_size, space_id, h5_err)

    ! create the dataset
    call h5dcreate_f(group_id, dset_name, data_type, space_id, &
             dset_id, h5_err)

    call h5sclose_f(space_id, h5_err)
    call h5dclose_f(dset_id, h5_err)

  end subroutine CreateH5RegionDataset


  subroutine WriteH5RegionRealDataset(group_id,  v_data,     &
                                        region_id, dset_name)
    !
    ! Write a dataset to the HDF5 observation region file.
    !
    ! Author: R. McKeown
    ! Date: 01/26/2022
    !

    implicit none

    ! arguments
    integer(HID_T), intent(in)               :: group_id
    PetscReal, intent(in)                    :: v_data(:)
    PetscInt, intent(in)                     :: region_id
    character(len=MAXWORDLENGTH), intent(in) :: dset_name

    ! local variables
    PetscMPIInt :: h5_err

    integer(HID_T) :: fapl_id
    integer(HID_T) :: mspace_id
    integer(HID_T) :: space_id
    integer(HID_T) :: dset_id
    integer(HSIZE_T), dimension(1) :: dset_size
    integer(HSIZE_T), dimension(1) :: start
    integer(HSIZE_T), dimension(1) :: extent

    PetscInt :: ds_rank = 1  ! data is always a 1D array


    ! to set up hyperslab use offset and num_cells
    start(1) = obs_region_offset(comm_rank+1, region_id)
    extent(1) = obs_region_numcells(comm_rank+1, region_id)

    ! size of whole dataset
    dset_size(1) = obs_region_total_numcells(region_id)

    ! open previously created data space
    call h5dopen_f(group_id, dset_name, dset_id, h5_err)

    ! create the memory data space for each process
    call h5screate_simple_f(ds_rank, extent, mspace_id, h5_err)

    ! select hyperslab
    call h5dget_space_f(dset_id, space_id, h5_err)
    call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, start, &
                   extent, h5_err)
    ! create property list for write
    call h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, h5_err)
    call h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_INDEPENDENT_F, h5_err)

    ! write the dataset
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, v_data, dset_size, &
                      h5_err, file_space_id = space_id, &
                      mem_space_id = mspace_id, xfer_prp = fapl_id)

    ! close some resources
    call h5pclose_f(fapl_id, h5_err)
    call h5sclose_f(mspace_id, h5_err)
    call h5sclose_f(space_id, h5_err)
    call h5dclose_f(dset_id, h5_err)


  end subroutine WriteH5RegionRealDataset


  subroutine WriteH5RegionIntegerDataset(group_id, v_data,    &
                                           region_id, dset_name)
    !
    ! Write a dataset to the HDF5 observation region file.
    !
    ! Author: R. McKeown
    ! Date: 01/26/2022
    !

    implicit none


    ! arguments
    integer(HID_T), intent(in)               :: group_id
    PetscInt, intent(in)                     :: v_data(:)
    PetscInt, intent(in)                     :: region_id
    character(len=MAXWORDLENGTH), intent(in) :: dset_name

    ! local variables
    PetscMPIInt :: h5_err

    integer(HID_T) :: fapl_id
    integer(HID_T) :: mspace_id
    integer(HID_T) :: space_id
    integer(HID_T) :: dset_id
    integer(HSIZE_T), dimension(1) :: start
    integer(HSIZE_T), dimension(1) :: extent
    integer(HSIZE_T), dimension(1) :: dset_size

    PetscInt :: ds_rank = 1  ! data is always a 1D array


    ! to set up hyperslab use offset and num_cells
    start(1) = obs_region_offset(comm_rank+1, region_id)
    extent(1) = obs_region_numcells(comm_rank+1, region_id)

    ! size of whole dataset
    dset_size(1) = obs_region_total_numcells(region_id)

    ! open previously created data space
    call h5dopen_f(group_id, dset_name, dset_id, h5_err)

    ! create the memory data space for each process
    call h5screate_simple_f(ds_rank, extent, mspace_id, h5_err)

    ! select hyperslab
    call h5dget_space_f(dset_id, space_id, h5_err)
    call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, start, &
                   extent, h5_err)
    ! create property list for write
    call h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, h5_err)
    call h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_INDEPENDENT_F, h5_err)

    ! write the dataset
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, v_data, dset_size, &
                      h5_err, file_space_id = space_id, &
                      mem_space_id = mspace_id, xfer_prp = fapl_id)

    ! close some resources
    call h5pclose_f(fapl_id, h5_err)
    call h5sclose_f(mspace_id, h5_err)
    call h5sclose_f(space_id, h5_err)
    call h5dclose_f(dset_id, h5_err)


  end subroutine WriteH5RegionIntegerDataset

  subroutine CreateOrOpenH5Group(main_id, group_id, group_name)
    !
    ! In case of restart, the file and group may already exist.
    ! Check first and then create the group if necessary.
    !
    ! Author: R. McKeown
    ! Date: 02/07/2022
    !

    implicit none

    character(len=MAXWORDLENGTH) :: group_name

    integer(HID_T) :: main_id
    integer(HID_T) :: group_id
    PetscMPIInt :: h5_err

    call h5eset_auto_f(OFF,h5_err)
    call h5gopen_f(main_id, group_name, group_id, h5_err)
    if (h5_err /= 0) then
      call h5gcreate_f(main_id, group_name, group_id,   &
                       h5_err, OBJECT_NAMELEN_DEFAULT_F)
    endif
    call h5eset_auto_f(ON,h5_err)


  end subroutine CreateOrOpenH5Group

end module Output_Obs_H5_module
