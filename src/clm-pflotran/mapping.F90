module Mapping_module

  use PFLOTRAN_Constants_module

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscmat.h"

  use petscsys
  use petscvec
  use petscmat
  use petscis
  
  implicit none
  private

  type, public  :: mapping_type

    !
    ! Linear Mapping from Source mesh to Destination mesh is matrix-vector product and
    ! can be written as:
    !
    !                W * s = d                                 Eq[1]
    !
    ! where W - Weight matrix      (nd x ns)
    !       s - Source vector      (ns x 1)
    !       d - Destination vector (ns x 1)
    !
    ! In CLM-PFLOTRAN coupling, s and d vectors are decomposed over multiple processors.
    ! The decomposition of vectors need not be in a contiguous order.
    !
    ! Each processor perfoms a local matrix-vector product:
    !
    !                W_loc * s_dloc = d_loc                    Eq[2]
    !
    !   - Obtains from Global Vec 's', a subset of local 's_dloc'
    !   - Performs a local matrix-vector product 
    !

    character(len=MAXSTRINGLENGTH) :: filename

    ! Note: IDs of source/destination mesh are 0-based

    ! Source mesh
    PetscInt           :: s_ncells_loc              ! # of local source mesh cells present
    PetscInt,pointer   :: s_ids_loc_nidx(:)         ! IDs of source mesh cells

    ! Destination mesh
    PetscInt           :: d_ncells_loc              ! # of local destination mesh cells present
    PetscInt           :: d_ncells_gh               ! # of ghost destination mesh cells present
    PetscInt           :: d_ncells_ghd              ! local+ghost=ghosted destination mesh cells

    ! natuaral-index starting with 0
    PetscInt,pointer   :: d_ids_ghd_nidx(:)         ! IDs of ghosted destination mesh cells present
    PetscInt,pointer   :: d_ids_nidx_sor(:)         ! Sorted Ghosted IDs of destination mesh cells
    PetscInt,pointer   :: d_nGhd2Sor(:)             ! Ghosted to Sorted
    PetscInt,pointer   :: d_nSor2Ghd(:)             ! Sorted to Ghosted
    PetscInt,pointer   :: d_loc_or_gh(:)            ! To flag if a cell is local(1) or ghost(0)

    ! Mapping from Source-to-Destination mesh
    PetscInt           :: s2d_s_ncells              ! # of source cells for mapping
    PetscInt           :: s2d_s_ncells_dis          ! # of "distinct" source cells for mapping
    PetscInt,pointer   :: s2d_s_ids_nidx(:)         ! IDs of source cells for mapping
    PetscInt,pointer   :: s2d_s_ids_nidx_dis(:)     ! IDs of "distinct" source cells for mapping

    ! Compressed Sparse Row (CSR) Matrix
    PetscReal,pointer  :: s2d_wts(:)                ! Wts for mapping
    PetscInt           :: s2d_nwts                  ! Number of wts
    PetscInt,pointer   :: s2d_jcsr(:)               ! J-th entry for CSR
    PetscInt,pointer   :: s2d_icsr(:)               ! I-th entry for CSR
    PetscInt,pointer   :: s2d_nonzero_rcount_csr(:) ! Non-Zero entries within a row

    Mat                :: wts_mat                   ! Sparse matrix for linear mapping
    VecScatter         :: s2d_scat_s_gb2disloc      ! Vec-Scatter of source mesh:Global to "distinct" local
    
    Vec                :: s_disloc_vec              ! Sequential vector to save "distinct" local
                                                    ! component of source vector

    ! Header information about number of layers mapped
    PetscInt           :: clm_nlevsoi               ! Number of CLM nlevsoi
    PetscInt           :: clm_nlevgrnd              ! Number of CLM nlevgrnd
    PetscInt           :: clm_nlev_mapped           ! Number of CLM layers mapped
    PetscInt           :: pflotran_nlev             ! Number of PFLOTRAN layers
    PetscInt           :: pflotran_nlev_mapped      ! Number of PFLOTRAN layers mapped

    type(mapping_type), pointer :: next

  end type mapping_type

  type, public :: mapping_list_type
    PetscInt                       :: nmap
    type(mapping_type), pointer    :: first
    type(mapping_type), pointer    :: last
  end type mapping_list_type

  public :: MappingCreate, &
            MappingSetSourceMeshCellIds, &
            MappingSetDestinationMeshCellIds, &
            MappingReadTxtFile, &
            MappingReadHDF5, &
            MappingDecompose, &
            MappingFindDistinctSourceMeshCellIds, &
            MappingCreateWeightMatrix, &
            MappingCreateScatterOfSourceMesh, &
            MappingSourceToDestination, &
            MappingListCreate, &
            MappingListAddToList, &
            MappingDestroy
contains

! ************************************************************************** !

  function MappingCreate()
  ! 
  ! This routine creates a mapping.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 

    implicit none

    type(mapping_type), pointer :: MappingCreate
    type(mapping_type), pointer :: map

    allocate(map)

    map%s_ncells_loc = 0
    nullify(map%s_ids_loc_nidx)

    ! Destination mesh
    map%d_ncells_loc = 0
    map%d_ncells_gh = 0
    map%d_ncells_ghd = 0

    ! natuaral-index starting with 0
    nullify(map%d_ids_ghd_nidx)
    nullify(map%d_ids_nidx_sor)
    nullify(map%d_nGhd2Sor)
    nullify(map%d_nSor2Ghd)
    nullify(map%d_loc_or_gh)

    ! Mapping from Source-to-Destination mesh
    map%s2d_s_ncells = 0
    map%s2d_s_ncells_dis = 0
    nullify(map%s2d_s_ids_nidx)
    nullify(map%s2d_s_ids_nidx_dis)

    ! Compressed Sparse Row (CSR) Matrix
    nullify(map%s2d_wts)
    map%s2d_nwts = 0
    nullify(map%s2d_jcsr)
    nullify(map%s2d_icsr)
    nullify(map%s2d_nonzero_rcount_csr)

    map%wts_mat = PETSC_NULL_MAT
    map%s2d_scat_s_gb2disloc = PETSC_NULL_VECSCATTER
    map%s_disloc_vec = PETSC_NULL_VEC

    map%clm_nlevsoi = 0
    map%clm_nlevgrnd = 0
    map%clm_nlev_mapped = 0
    map%pflotran_nlev = 0
    map%pflotran_nlev_mapped = 0
    
    nullify(map%next)

    MappingCreate => map

  end function MappingCreate

! ************************************************************************** !

  subroutine MappingSetSourceMeshCellIds( map, &
                                          ncells, &
                                          cell_ids &
                                        )
  ! 
  ! This routine sets cell ids source mesh.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 

    implicit none

    type(mapping_type), pointer :: map

    PetscInt                    :: ncells,ii
    PetscInt,pointer            :: cell_ids(:)

    map%s_ncells_loc = ncells
    allocate(map%s_ids_loc_nidx(map%s_ncells_loc))

    do ii = 1,ncells
      map%s_ids_loc_nidx(ii) = cell_ids(ii)
    enddo

  end subroutine MappingSetSourceMeshCellIds

! ************************************************************************** !

  subroutine MappingSetDestinationMeshCellIds(map, &
                                              ncells_loc, &
                                              ncells_gh, &
                                              cell_ids_ghd, &
                                              loc_or_gh &
                                              )
  ! 
  ! This routine sets the cell ids of destination mesh.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 

    implicit none

    type(mapping_type), pointer :: map

    PetscInt                    :: ncells_loc, ncells_gh
    PetscInt,pointer            :: cell_ids_ghd(:), loc_or_gh(:)

    PetscInt                    :: ii
    PetscInt, pointer           :: index(:),rev_index(:)
    PetscErrorCode               :: ierr

    ! Initialize
    map%d_ncells_loc = ncells_loc
    map%d_ncells_gh  = ncells_gh
    map%d_ncells_ghd = ncells_loc + ncells_gh

    ! Allocate memory
    allocate(map%d_ids_ghd_nidx(map%d_ncells_ghd))
    allocate(map%d_ids_nidx_sor(map%d_ncells_ghd))
    allocate(map%d_loc_or_gh(   map%d_ncells_ghd))
    allocate(map%d_nGhd2Sor(    map%d_ncells_ghd))
    allocate(map%d_nSor2Ghd(    map%d_ncells_ghd))
    allocate(index(             map%d_ncells_ghd))
    allocate(rev_index(         map%d_ncells_ghd))

    do ii=1,ncells_loc+ncells_gh
      map%d_ids_ghd_nidx(ii) = cell_ids_ghd(ii)
      map%d_loc_or_gh(ii)    = loc_or_gh(ii)
      index(ii)              = ii
      rev_index(ii)          = ii
    enddo

    ! Sort cell_ids_ghd
    index = index - 1 ! Needs to be 0-based
    call PetscSortIntWithPermutation(map%d_ncells_ghd,cell_ids_ghd,index,ierr)
    index = index + 1

    do ii=1,ncells_loc+ncells_gh
      map%d_ids_nidx_sor(ii) = cell_ids_ghd(index(ii))
      map%d_nGhd2Sor(ii)     = index(ii)
    enddo

    ! Sort the index
    rev_index = rev_index - 1
    call PetscSortIntWithPermutation(map%d_ncells_ghd,index,rev_index,ierr)
    map%d_nSor2Ghd = rev_index

    ! Free memory
    deallocate(index)
    deallocate(rev_index)


  end subroutine MappingSetDestinationMeshCellIds

! ************************************************************************** !

  subroutine MappingReadTxtFile(map,map_filename,option)
  ! 
  ! This routine reads a ASCII mapping file.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    use Input_Aux_module
    use Option_module
    use String_module
    use PFLOTRAN_Constants_module
    
    implicit none
    
    ! argument
    type(mapping_type), pointer     :: map
    character(len=MAXSTRINGLENGTH)  :: map_filename
    type(option_type), pointer      :: option
    
    ! local variables
    type(input_type), pointer       :: input
    character(len=MAXSTRINGLENGTH)  :: card
    character(len=MAXSTRINGLENGTH)  :: string
    character(len=MAXWORDLENGTH)    :: word
    character(len=MAXSTRINGLENGTH)  :: hint
    PetscInt                        :: temp_int
    PetscInt                        :: temp_int_array(5)
    PetscInt                        :: nwts
    PetscInt                        :: nread
    PetscInt                        :: remainder, prev_row
    PetscInt                        :: irank, ii
    
    ! Only used by io_rank
    PetscInt                        :: nwts_tmp
    PetscInt,pointer                :: wts_row_tmp(:), wts_col_tmp(:)
    PetscReal,pointer               :: wts_tmp(:)
    PetscInt                        :: nheader
    
    PetscMPIInt                     :: status_mpi(MPI_STATUS_SIZE)
    PetscErrorCode                  :: ierr

    card     = 'MppingReadTxtFile'

    ! Read ASCII file through io_rank and communicate to other ranks
    if(option%myrank == option%io_rank) then

      input => InputCreate(20,map_filename,option)

      nwts     = -1
      prev_row = -1

      ! Read first six entries in the mapping file.
      do nheader = 1, 6
        call InputReadPflotranString(input,option)
        call InputReadWord(input,option,card,PETSC_TRUE)
        call StringToLower(card)

        select case (trim(card))
          case('clm_nlevsoi')
            hint = 'clm_nlevsoi'
            call InputReadInt(input,option,map%clm_nlevsoi)
            call InputErrorMsg(input,option,'CLM nlevsoi',hint)
          case('clm_nlevgrnd')
            hint = 'clm_nlevgrnd'
            call InputReadInt(input,option,map%clm_nlevgrnd)
            call InputErrorMsg(input,option,'CLM nlevgrnd',hint)
          case('clm_nlev_mapped')
            hint = 'clm_nlev_mapped'
            call InputReadInt(input,option,map%clm_nlev_mapped)
            call InputErrorMsg(input,option,'CLM nlev mapped',hint)
          case('pflotran_nlev')
            hint = 'pflotran_nlev'
            call InputReadInt(input,option,map%pflotran_nlev)
            call InputErrorMsg(input,option,'PFLOTRAN nlev',hint)
          case('pflotran_nlev_mapped')
            hint = 'pflotran_nlev_mapped'
            call InputReadInt(input,option,map%pflotran_nlev_mapped)
            call InputErrorMsg(input,option,'PFLOTRAN nlev mapped',hint)
          case('num_weights')
            hint = 'num_weights'
            call InputReadInt(input,option,nwts)
            call InputErrorMsg(input,option,'Number of weights',hint)
            write(*,*),'nwts = ',nwts
          case default
            option%io_buffer = 'Unrecognized keyword "' // trim(card) // &
              '" in explicit grid file.'
            call PrintErrMsgByRank(option)
        end select
      enddo
      
      nwts_tmp = nwts/option%mycommsize
      remainder= nwts - nwts_tmp*option%mycommsize
      
      allocate(wts_row_tmp(nwts_tmp + 1))
      allocate(wts_col_tmp(nwts_tmp + 1))
      allocate(wts_tmp(    nwts_tmp + 1))

      do irank = 0,option%mycommsize-1
        
        ! Determine the number of row to be read
        nread = nwts_tmp
        if(irank<remainder) nread = nread+1
        
        ! Read the data
        do ii = 1,nread
          call InputReadPflotranString(input,option)
          call InputReadInt(input,option,wts_row_tmp(ii))
          call InputReadInt(input,option,wts_col_tmp(ii))
          call InputReadDouble(input,option,wts_tmp(ii))
          
          !Perform checks on the data read
          if(wts_row_tmp(ii) < 1) then
            write(*,string),'Row entry for ii = ',ii,' less than 1'
            option%io_buffer = string
            call PrintErrMsg(option)
          endif
          
          if(wts_col_tmp(ii) < 1) then
            write(*,string),'Col entry for ii = ',ii,' less than 1'
            option%io_buffer = string
            call PrintErrMsg(option)
          endif
          
          if((wts_tmp(ii) < 0.d0).or.(wts_tmp(ii) > 1.d0)) then
            write(*,string),'Invalid wt value for ii = ',ii
            option%io_buffer = string
            call PrintErrMsg(option)
          endif
          
          ! ensure that row values in the data are stored in ascending order
          if(wts_row_tmp(ii) < prev_row) then
            write(*,string),'Row value in the mapping data not store in ascending order: ii ',ii
            option%io_buffer = string
            call PrintErrMsg(option)
          endif
          prev_row = wts_row_tmp(ii)
          
          ! Convert row/col values to 0-based
          wts_row_tmp(ii) = wts_row_tmp(ii) - 1
          wts_col_tmp(ii) = wts_col_tmp(ii) - 1
          
        enddo
        
        ! Save data locally
        if (irank == option%myrank) then
          
            map%s2d_nwts = nread
            allocate(map%s2d_icsr(map%s2d_nwts))
            allocate(map%s2d_jcsr(map%s2d_nwts))
            allocate(map%s2d_wts( map%s2d_nwts))
            
            do ii = 1,map%s2d_nwts
              map%s2d_icsr(ii) = wts_row_tmp(ii)
              map%s2d_jcsr(ii) = wts_col_tmp(ii)
              map%s2d_wts(ii)  = wts_tmp(ii)
            enddo
        
        else

          ! Otherwise communicate data to other ranks
          call MPI_Send(nread,1,MPI_INTEGER,irank,option%myrank,option%mycomm, &
                        ierr)
          call MPI_Send(wts_row_tmp,nread,MPI_INTEGER,irank,option%myrank, &
                        option%mycomm,ierr)
          call MPI_Send(wts_col_tmp,nread,MPI_INTEGER,irank,option%myrank, &
                        option%mycomm,ierr)
          call MPI_Send(wts_tmp,nread,MPI_DOUBLE_PRECISION,irank,option%myrank, &
                        option%mycomm,ierr)
          
        endif
        
      enddo

      deallocate(wts_row_tmp)
      deallocate(wts_col_tmp)
      deallocate(wts_tmp)
      call InputDestroy(input)
      
    else
      ! Other ranks receive data from io_rank
      
      ! Get the number of data
      call MPI_Recv(map%s2d_nwts,1,MPI_INTEGER,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)

      ! Allocate memory
      allocate(map%s2d_icsr(map%s2d_nwts))
      allocate(map%s2d_jcsr(map%s2d_nwts))
      allocate(map%s2d_wts( map%s2d_nwts))
      
      call MPI_Recv(map%s2d_icsr,map%s2d_nwts,MPI_INTEGER,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)
      call MPI_Recv(map%s2d_jcsr,map%s2d_nwts,MPI_INTEGER,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)
      call MPI_Recv(map%s2d_wts,map%s2d_nwts,MPI_DOUBLE_PRECISION,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)
                    
    endif

  ! Broadcast from root information regarding CLM/PFLOTRAN num soil layers
  temp_int_array(1) = map%clm_nlevsoi
  temp_int_array(2) = map%clm_nlevgrnd
  temp_int_array(3) = map%clm_nlev_mapped
  temp_int_array(4) = map%pflotran_nlev
  temp_int_array(5) = map%pflotran_nlev_mapped

  call MPI_Bcast(temp_int_array,FIVE_INTEGER,MPI_INTEGER,option%io_rank, &
                 option%mycomm,ierr)
    
  map%clm_nlevsoi = temp_int_array(1)
  map%clm_nlevgrnd = temp_int_array(2)
  map%clm_nlev_mapped = temp_int_array(3)
  map%pflotran_nlev = temp_int_array(4)
  map%pflotran_nlev_mapped = temp_int_array(5)

  end subroutine MappingReadTxtFile

! ************************************************************************** !

  subroutine MappingReadHDF5(map,map_filename,option)
  ! 
  ! This routine reads a mapping file in HDF5 format.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/8/2013
  ! 
  
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

    use Input_Aux_module
    use Option_module
    
    implicit none
    
    ! argument
    type(mapping_type), pointer :: map
    character(len=MAXSTRINGLENGTH) :: map_filename
    type(option_type), pointer :: option
    
    ! local
    PetscMPIInt       :: hdf5_err
    PetscMPIInt       :: rank_mpi
    PetscInt          :: ndims
    PetscInt          :: istart, iend, ii, jj
    PetscInt          :: num_cells_local
    PetscInt          :: num_cells_local_save
    PetscInt          :: num_vertices_local
    PetscInt          :: num_vertices_local_save
    PetscInt          :: remainder
    PetscInt,pointer  :: int_buffer(:)
    PetscReal,pointer :: double_buffer(:)
    PetscInt, parameter :: max_nvert_per_cell = 8  
    PetscErrorCode    :: ierr

    character(len=MAXSTRINGLENGTH) :: group_name
    character(len=MAXSTRINGLENGTH) :: dataset_name

#if defined(PETSC_HAVE_HDF5)
    integer(HID_T) :: file_id
    integer(HID_T) :: grp_id, grp_id2
    integer(HID_T) :: prop_id
    integer(HID_T) :: data_set_id
    integer(HID_T) :: file_space_id
    integer(HID_T) :: data_space_id
    integer(HID_T) :: memory_space_id
    integer(HSIZE_T) :: num_data_in_file
    integer(HSIZE_T) :: dims_h5(1), max_dims_h5(1)
    integer(HSIZE_T) :: offset(1), length(1), stride(1), block(1), dims(1)
#endif

    ! Initialize FORTRAN predefined datatypes
    call h5open_f(hdf5_err)

    ! Setup file access property with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)

#ifndef SERIAL_HDF5
    call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif

    ! Open the file collectively
    call h5fopen_f(map_filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
    call h5pclose_f(prop_id,hdf5_err)
    
    !
    ! /row
    !
    
    ! Open group
    group_name = "/row"
    option%io_buffer = 'Opening group: ' // trim(group_name)
    call PrintMsg(option)

    ! Open dataset
    call h5dopen_f(file_id,"row",data_set_id,hdf5_err)

    ! Get dataset's dataspace
    call h5dget_space_f(data_set_id,data_space_id,hdf5_err)
    
    ! Get number of dimensions and check
    call h5sget_simple_extent_ndims_f(data_space_id,ndims,hdf5_err)
    if (ndims /= 1) then
      option%io_buffer='Dimension of row dataset in ' // trim(map_filename) // &
            ' is not equal to 1.'
      call PrintErrMsg(option)
    endif

    ! Get dimensions of dataset
    call h5sget_simple_extent_dims_f(data_space_id,dims_h5,max_dims_h5, &
                                     hdf5_err)
    
    ! Determine the number of cells each that will be saved on each processor
    map%s2d_nwts=INT(dims_h5(1))/option%mycommsize
    remainder=INT(dims_h5(1))-map%s2d_nwts*option%mycommsize
    if (option%myrank < remainder) map%s2d_nwts=map%s2d_nwts + 1
    
    ! allocate array to store vertices for each cell
    allocate(map%s2d_icsr(map%s2d_nwts))
    allocate(map%s2d_jcsr(map%s2d_nwts))
    allocate(map%s2d_wts( map%s2d_nwts))

    ! Find istart and iend
    istart = 0
    iend   = 0
    call MPI_Exscan(map%s2d_nwts,istart,ONE_INTEGER_MPI, &
                    MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(map%s2d_nwts,iend,ONE_INTEGER_MPI, &
                  MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    
    ! Determine the length and offset of data to be read by each processor
    length(1) = iend-istart
    offset(1) = istart
    
    !
    rank_mpi = 1
    memory_space_id = -1
    
    ! Create data space for dataset
    call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err)

    ! Select hyperslab
    call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
    call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                               hdf5_err)
    
    ! Initialize data buffer
    allocate(int_buffer(length(1)))
    
    ! Create property list
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
    
    ! Read the dataset collectively
    call h5dread_f(data_set_id, H5T_NATIVE_INTEGER, int_buffer, &
                   dims_h5, hdf5_err, memory_space_id, data_space_id)
    
    ! Convert 1-based to 0-based
    map%s2d_icsr = int_buffer-1

    call h5dclose_f(data_set_id, hdf5_err)

    !
    ! /col
    !
    
    ! Open group
    group_name = "/col"
    option%io_buffer = 'Opening group: ' // trim(group_name)
    call PrintMsg(option)

    ! Open dataset
    call h5dopen_f(file_id,"col",data_set_id,hdf5_err)

    ! Get dataset's dataspace
    call h5dget_space_f(data_set_id,data_space_id,hdf5_err)
    
    ! Get number of dimensions and check
    call h5sget_simple_extent_ndims_f(data_space_id,ndims,hdf5_err)
    if (ndims /= 1) then
      option%io_buffer='Dimension of row dataset in ' // trim(map_filename) // &
            ' is not equal to 1.'
      call PrintErrMsg(option)
    endif

    ! Get dimensions of dataset
    call h5sget_simple_extent_dims_f(data_space_id,dims_h5,max_dims_h5, &
                                     hdf5_err)
    
    ! Determine the number of cells each that will be saved on each processor
    map%s2d_nwts=INT(dims_h5(1))/option%mycommsize
    remainder=INT(dims_h5(1))-map%s2d_nwts*option%mycommsize
    if (option%myrank < remainder) map%s2d_nwts=map%s2d_nwts + 1
    
    ! Find istart and iend
    istart = 0
    iend   = 0
    call MPI_Exscan(map%s2d_nwts,istart,ONE_INTEGER_MPI, &
                    MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(map%s2d_nwts,iend,ONE_INTEGER_MPI, &
                  MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    
    ! Determine the length and offset of data to be read by each processor
    length(1) = iend-istart
    offset(1) = istart
    
    !
    rank_mpi = 1
    memory_space_id = -1
    
    ! Create data space for dataset
    call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err)

    ! Select hyperslab
    call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
    call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                               hdf5_err)
    
    ! Create property list
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
    
    ! Read the dataset collectively
    call h5dread_f(data_set_id, H5T_NATIVE_INTEGER, int_buffer, &
                   dims_h5, hdf5_err, memory_space_id, data_space_id)
    
    ! Convert 1-based to 0-based
    map%s2d_jcsr = int_buffer-1

    call h5dclose_f(data_set_id, hdf5_err)

    !
    ! /S
    !
    
    ! Open group
    group_name = "/S"
    option%io_buffer = 'Opening group: ' // trim(group_name)
    call PrintMsg(option)

    ! Open dataset
    call h5dopen_f(file_id,"S",data_set_id,hdf5_err)

    ! Get dataset's dataspace
    call h5dget_space_f(data_set_id,data_space_id,hdf5_err)
    
    ! Get number of dimensions and check
    call h5sget_simple_extent_ndims_f(data_space_id,ndims,hdf5_err)
    if (ndims /= 1) then
      option%io_buffer='Dimension of row dataset in ' // trim(map_filename) // &
            ' is not equal to 1.'
      call PrintErrMsg(option)
    endif

    ! Get dimensions of dataset
    call h5sget_simple_extent_dims_f(data_space_id,dims_h5,max_dims_h5, &
                                     hdf5_err)
    
    ! Determine the number of cells each that will be saved on each processor
    map%s2d_nwts=INT(dims_h5(1))/option%mycommsize
    remainder=INT(dims_h5(1))-map%s2d_nwts*option%mycommsize
    if (option%myrank < remainder) map%s2d_nwts=map%s2d_nwts + 1
    
    ! Find istart and iend
    istart = 0
    iend   = 0
    call MPI_Exscan(map%s2d_nwts,istart,ONE_INTEGER_MPI, &
                    MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(map%s2d_nwts,iend,ONE_INTEGER_MPI, &
                  MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    
    ! Determine the length and offset of data to be read by each processor
    length(1) = iend-istart
    offset(1) = istart
    
    !
    rank_mpi = 1
    memory_space_id = -1
    
    ! Create data space for dataset
    call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err)

    ! Select hyperslab
    call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
    call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                               hdf5_err)
    
    ! Initialize data buffer
    allocate(double_buffer(length(1)))
    
    ! Create property list
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
    
    ! Read the dataset collectively
    call h5dread_f(data_set_id, H5T_NATIVE_DOUBLE, double_buffer, &
                   dims_h5, hdf5_err, memory_space_id, data_space_id)
    
    map%s2d_wts = double_buffer

    call h5dclose_f(data_set_id, hdf5_err)

    ! Close file
    call h5fclose_f(file_id, hdf5_err)
    call h5close_f(hdf5_err)

    ! Free memory
    deallocate(int_buffer)
    deallocate(double_buffer)

  end subroutine MappingReadHDF5

! ************************************************************************** !

  subroutine MappingDecompose(map,mycomm)
  ! 
  ! This routine decomposes the mapping when running on more than processor,
  ! while accounting for different domain decomposition of source and
  ! destination grid.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    use Option_module
    
    implicit none
    
    ! argument
    type(mapping_type), pointer   :: map
    PetscMPIInt :: mycomm
    
    ! local variables
    Vec                           :: nonzero_row_count_vec           ! MPI
    Vec                           :: nonzero_row_count_loc_vec       ! Seq
    Vec                           :: cumsum_nonzero_row_count_vec    ! MPI
    Vec                           :: cumsum_nonzero_row_count_loc_vec! Seq
    
    Vec                           :: row_vec, col_vec, wts_vec             ! MPI
    Vec                           :: row_loc_vec, col_loc_vec, wts_loc_vec ! Seq
    
    IS                            :: is_from, is_to
    
    VecScatter                    :: vec_scat
    
    PetscViewer                   :: viewer

    PetscInt, pointer             :: row(:)
    PetscScalar,pointer           :: row_count(:)
    PetscInt, pointer             :: tmp_int_array(:)
    PetscInt                      :: ii,jj,kk
    PetscInt                      :: nrow,cumsum_start,count
    
    PetscScalar,pointer           :: vloc1(:),vloc2(:),vloc3(:),vloc4(:)
    PetscErrorCode                :: ierr
    character(len=MAXSTRINGLENGTH):: string

    ! 0) Save dataset related to wts in MPI-Vecs
    call VecCreateMPI(mycomm,map%s2d_nwts,PETSC_DECIDE,row_vec,ierr)
    call VecCreateMPI(mycomm,map%s2d_nwts,PETSC_DECIDE,col_vec,ierr)
    call VecCreateMPI(mycomm,map%s2d_nwts,PETSC_DECIDE,wts_vec,ierr)

    call VecGetArrayF90(row_vec,vloc1,ierr)
    call VecGetArrayF90(col_vec,vloc2,ierr)
    call VecGetArrayF90(wts_vec,vloc3,ierr)

    do ii = 1,map%s2d_nwts
       vloc1(ii) = map%s2d_icsr(ii)
       vloc2(ii) = map%s2d_jcsr(ii)
       vloc3(ii) = map%s2d_wts(ii)
    enddo

    call VecRestoreArrayF90(row_vec,vloc1,ierr)
    call VecRestoreArrayF90(col_vec,vloc2,ierr)
    call VecRestoreArrayF90(wts_vec,vloc3,ierr)

    ! 1) For each cell of destination mesh, find the number of source mesh cells
    !    overlapped.
    !                             OR
    !    Number of non-zero entries for each row of the global W matrix

    ! Create a MPI vector
    call VecCreateMPI(mycomm,map%d_ncells_loc,PETSC_DECIDE, &
                      nonzero_row_count_vec,ierr)

    call VecSet(nonzero_row_count_vec,0.d0,ierr)
    
    ! Find non-zero entries for each of the W matrix with the locally saved data
    allocate(row_count(map%s2d_nwts))
    allocate(row(map%s2d_nwts))

    ii = 1
    nrow = 1
    row(nrow)       = map%s2d_icsr(ii) 
    row_count(nrow) = 1.d0
    
    do ii = 2,map%s2d_nwts
      if (map%s2d_icsr(ii) == row(nrow)) then
        row_count(nrow) = row_count(nrow) + 1
      else
        nrow = nrow + 1
        row(nrow)       = map%s2d_icsr(ii) 
        row_count(nrow) = 1
      endif
    enddo

    ! Save values in the MPI vector
    call VecSetValues(nonzero_row_count_vec,nrow,row,row_count,ADD_VALUES,ierr)
    call VecAssemblyBegin(nonzero_row_count_vec,ierr)
    call VecAssemblyEnd(nonzero_row_count_vec,ierr)
    deallocate(row)
    deallocate(row_count)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'nonzero_row_count_vec.out', viewer, ierr)
    call VecView(nonzero_row_count_vec, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
    
    ! 2) Find cummulative sum of the nonzero_row_count_vec

    ! Create a MPI vector
    call VecCreateMPI(mycomm,map%d_ncells_loc,PETSC_DECIDE, &
                      cumsum_nonzero_row_count_vec,ierr)
    call VecGetArrayF90(nonzero_row_count_vec,vloc1,ierr)
    call VecGetArrayF90(cumsum_nonzero_row_count_vec,vloc2,ierr)
    
    ii = 1
    vloc2(ii) = vloc1(ii)
    do ii = 2,map%d_ncells_loc
      vloc2(ii) = vloc2(ii-1) + vloc1(ii)
    enddo

    cumsum_start = 0
    call MPI_Exscan(INT(vloc2(map%d_ncells_loc)),cumsum_start,ONE_INTEGER_MPI, &
                    MPIU_INTEGER,MPI_SUM,mycomm,ierr)

    do ii = 1,map%d_ncells_loc
      vloc2(ii) = vloc2(ii) + cumsum_start
    enddo
    
    call VecRestoreArrayF90(nonzero_row_count_vec,vloc1,ierr)
    call VecRestoreArrayF90(cumsum_nonzero_row_count_vec,vloc2,ierr)

#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'cumsum_nonzero_row_count_vec.out', viewer, ierr)
    call VecView(cumsum_nonzero_row_count_vec, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
    
    ! 3) On a given processor, in order to map a variable from source mesh to
    !    destination mesh, find 
    !    - the number source mesh cells required
    !    - cell ids of source mesh
    !
    !    Use VecScatter() to save portion of nonzero_row_count_vec and
    !    cumsum_nonzero_row_count_vec corresponding which correspond to
    !    ghosted (local+ghost) cell ids destination mesh on a given proc.
    
    !
    call VecCreateSeq(PETSC_COMM_SELF,map%d_ncells_ghd, &
                      nonzero_row_count_loc_vec,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,map%d_ncells_ghd, &
                      cumsum_nonzero_row_count_loc_vec,ierr)

    ! Create index sets (IS) for VecScatter()
    allocate(tmp_int_array(map%d_ncells_ghd))
    do ii = 1,map%d_ncells_ghd
      tmp_int_array(ii) = ii - 1
    enddo
    call ISCreateBlock(mycomm,1,map%d_ncells_ghd,tmp_int_array, &
                       PETSC_COPY_VALUES,is_to,ierr)
    deallocate(tmp_int_array)
    
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'is_to.out', viewer, ierr)
    call ISView(is_to,viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
    
    call ISCreateBlock(mycomm,1,map%d_ncells_ghd, &
                      map%d_ids_ghd_nidx, PETSC_COPY_VALUES,is_from,ierr)

#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'is_from.out', viewer, ierr)
    call ISView(is_from, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Create VecScatter()
    call VecScatterCreate(nonzero_row_count_vec,is_from, &
                          nonzero_row_count_loc_vec,is_to, &
                          vec_scat,ierr)
    call ISDestroy(is_from,ierr)
    call ISDestroy(is_to,ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'vec_scat.out', viewer, ierr)
    call VecScatterView(vec_scat, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
    
    ! Scatter the data
    call VecScatterBegin(vec_scat,nonzero_row_count_vec, &
                        nonzero_row_count_loc_vec,INSERT_VALUES,SCATTER_FORWARD, &
                        ierr)
    call VecScatterEnd(vec_scat,nonzero_row_count_vec, &
                      nonzero_row_count_loc_vec,INSERT_VALUES,SCATTER_FORWARD, &
                      ierr)
    call VecScatterBegin(vec_scat,cumsum_nonzero_row_count_vec, &
                        cumsum_nonzero_row_count_loc_vec,INSERT_VALUES,SCATTER_FORWARD, &
                        ierr)
    call VecScatterEnd(vec_scat,cumsum_nonzero_row_count_vec, &
                      cumsum_nonzero_row_count_loc_vec,INSERT_VALUES,SCATTER_FORWARD, &
                      ierr)
    call VecScatterDestroy(vec_scat,ierr)
    
    call VecGetArrayF90(nonzero_row_count_loc_vec,vloc1,ierr)
    call VecGetArrayF90(cumsum_nonzero_row_count_loc_vec,vloc2,ierr)
    
    allocate(map%s2d_nonzero_rcount_csr(map%d_ncells_ghd))
    
    ! For each destination, save the number of overlapped source cells
    count = 0
    do ii = 1,map%d_ncells_ghd
      count = count + INT(vloc1(ii))
      map%s2d_nonzero_rcount_csr(ii) = INT(vloc1(ii))
    enddo
    map%s2d_s_ncells = count

    ! Allocate memory
    deallocate(map%s2d_wts)
      
    allocate(map%s2d_s_ids_nidx(map%s2d_s_ncells))
    allocate(map%s2d_wts(map%s2d_s_ncells))
    allocate(tmp_int_array(map%s2d_s_ncells))
      
    ! For each cell in destination mesh, save indices of MPI Vectors, which
    ! contain data read from mapping file, for all overlapped cells of 
    ! of source mesh
    kk = 0
    do ii = 1,map%d_ncells_ghd
       do jj = 1,INT(vloc1(ii))
          kk = kk + 1
          tmp_int_array(kk) = INT(vloc2(ii)) - INT(vloc1(ii)) + jj - 1
       enddo
    enddo

    ! Create an index set to scatter from
    call ISCreateBlock(mycomm,1,map%s2d_s_ncells,tmp_int_array, &
         PETSC_COPY_VALUES,is_from,ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'is_from2.out', viewer, ierr)
    call ISView(is_from,viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    do ii=1,map%s2d_s_ncells
       tmp_int_array(ii) = ii-1
    enddo

    call ISCreateBlock(mycomm,1,map%s2d_s_ncells,tmp_int_array, &
         PETSC_COPY_VALUES,is_to,ierr)
    deallocate(tmp_int_array)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'is_to2.out', viewer, ierr)
    call ISView(is_to,viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Allocate memory
    call VecCreateSeq(PETSC_COMM_SELF,map%s2d_s_ncells,row_loc_vec,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,map%s2d_s_ncells,col_loc_vec,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,map%s2d_s_ncells,wts_loc_vec,ierr)

    ! Create scatter context
    call VecScatterCreate(row_vec,is_from,row_loc_vec,is_to,vec_scat,ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'scatter_wts_data.out', viewer, ierr)
    call VecScatterView(vec_scat, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Scatter the data
    call VecScatterBegin(vec_scat,col_vec,col_loc_vec, INSERT_VALUES, &
         SCATTER_FORWARD,ierr)
    call VecScatterEnd(vec_scat,col_vec,col_loc_vec, INSERT_VALUES, &
         SCATTER_FORWARD,ierr)
    call VecScatterBegin(vec_scat,wts_vec,wts_loc_vec, INSERT_VALUES, &
         SCATTER_FORWARD,ierr)
    call VecScatterEnd(vec_scat,wts_vec,wts_loc_vec, INSERT_VALUES, &
         SCATTER_FORWARD,ierr)

    ! Attach to the local copy of the scatterd data
    call VecGetArrayF90(col_loc_vec,vloc3,ierr)
    call VecGetArrayF90(wts_loc_vec,vloc4,ierr)

    ! Save the scattered data
    do ii = 1,map%s2d_s_ncells
       map%s2d_s_ids_nidx(ii) = INT(vloc3(ii))
       map%s2d_wts(ii)        = vloc4(ii)
    enddo

    ! Restore data
    call VecRestoreArrayF90(col_loc_vec,vloc3,ierr)
    call VecRestoreArrayF90(wts_loc_vec,vloc4,ierr)

    ! Free memory
    call VecDestroy(row_loc_vec,ierr)
    call VecDestroy(col_loc_vec,ierr)
    call VecDestroy(wts_loc_vec,ierr)
    
    ! Restore data
    call VecRestoreArrayF90(nonzero_row_count_loc_vec,vloc1,ierr)
    call VecRestoreArrayF90(cumsum_nonzero_row_count_loc_vec,vloc2,ierr)

    ! Free memory
    call VecDestroy(nonzero_row_count_vec,ierr)
    call VecDestroy(cumsum_nonzero_row_count_vec,ierr)
    call VecDestroy(nonzero_row_count_loc_vec,ierr)
    call VecDestroy(cumsum_nonzero_row_count_loc_vec,ierr)
    
  end subroutine MappingDecompose

! ************************************************************************** !

  subroutine MappingFindDistinctSourceMeshCellIds(map)
  ! 
  ! This routine finds distinct cell ids of source mesh
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    implicit none
    
    ! argument
    type(mapping_type), pointer   :: map
    
    ! local variables
    PetscErrorCode               :: ierr
    PetscInt                     :: ii,jj,kk,count
    PetscInt, pointer            :: index(:)
    PetscInt,pointer             :: int_array(:),int_array2(:)
    PetscInt,pointer             :: int_array3(:),int_array4(:)

    Vec :: xx, yy
    
    ! No overlapped cells with Source Mesh, then return
    if(map%s2d_s_ncells == 0) then
       map%s2d_s_ncells_dis = 0
       call VecCreateSeq(PETSC_COMM_SELF, map%s2d_s_ncells_dis, map%s_disloc_vec, ierr)
       return
    end if

    ! Allocate memory
    allocate(map%s2d_jcsr( map%s2d_s_ncells))
    allocate(int_array (map%s2d_s_ncells))
    allocate(int_array2(map%s2d_s_ncells))
    allocate(int_array3(map%s2d_s_ncells))
    allocate(int_array4(map%s2d_s_ncells))

    !
    ! Follows Glenn's approach in unstructured code to remove duplicate
    !   vertices.
    !
    ! map%s2d_s_ids_loc_nidx - Contains source mesh cell ids with duplicate entries
    ! The algo is explained below: 
    !     int_array  : Cell-ids
    !     int_array2 : Sorted index
    !     int_array3 : Distinct values
    !     int_array4 : Indices w.r.t. new distinct value vector
    !
    !  ii  int_array  int_array2  int_array3  int_array4
    !   1     90         6           70          3
    !   2    100         3           80          4
    !   3     80         1           90          2
    !   4    100         2          100          4
    !   5    101         4          101          5
    !   6     70         5                       1
    !
  
    do ii = 1,map%s2d_s_ncells
      int_array(ii)  = map%s2d_s_ids_nidx(ii)
      int_array2(ii) = ii
    enddo
    
    int_array2 = int_array2 - 1
    call PetscSortIntWithPermutation(map%s2d_s_ncells,int_array,int_array2,ierr)
    int_array2 = int_array2 + 1
    
    int_array3 = 0
    int_array4 = 0
    count = 1
    int_array3(1)             = int_array(int_array2(1))
    int_array4(int_array2(1)) = count
    
    do ii=2,map%s2d_s_ncells
      jj = int_array(int_array2(ii))
      if (jj > int_array3(count)) then
        count = count + 1
        int_array3(count) = jj
      endif
      int_array4(int_array2(ii)) = count 
    enddo
    
    ! Change 1-based index to 0-based index
    int_array4 = int_array4 - 1
    
    map%s2d_s_ncells_dis = count
    ! Save the distinct ids
    allocate(map%s2d_s_ids_nidx_dis(map%s2d_s_ncells_dis))
    
    map%s2d_s_ids_nidx_dis(1:count) = int_array3(1:count)
    map%s2d_jcsr = int_array4

    ! Create a sequential vector
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_s_ncells_dis, map%s_disloc_vec, ierr)

    ! Free memory
    deallocate(int_array)
    deallocate(int_array2)
    deallocate(int_array3)
    deallocate(int_array4)

  end subroutine MappingFindDistinctSourceMeshCellIds

! ************************************************************************** !

  subroutine MappingCreateWeightMatrix(map,myrank)
  ! 
  ! This routine creates a weight matrix to map data from source to destination
  ! grid.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    use Option_module
    
    implicit none
    
    ! argument
    type(mapping_type), pointer   :: map
    PetscMPIInt :: myrank

    ! local variables
    PetscInt, pointer            :: index(:)
    PetscInt                     :: ii,jj,kk
    PetscErrorCode               :: ierr
    character(len=MAXSTRINGLENGTH)     :: string  
    PetscViewer :: viewer

    allocate(index(map%s2d_s_ncells))

    kk = 0
    do ii = 1,map%d_ncells_ghd
       do jj = 1,map%s2d_nonzero_rcount_csr(ii)
          kk = kk + 1
          index(kk) = ii -1
       enddo
    enddo

    !
    ! size(wts_mat) = [d_ncells_ghd x s2d_s_ncells_dis]
    !
    call MatCreateSeqAIJ(PETSC_COMM_SELF, &
         map%d_ncells_ghd,                & ! m
         map%s2d_s_ncells_dis,            & ! n
         PETSC_DEFAULT_INTEGER,           & ! nz
         map%s2d_nonzero_rcount_csr,      & ! nnz
         map%wts_mat, ierr)

    do ii = 1,map%s2d_s_ncells
       call MatSetValues(map%wts_mat,1,index(ii),1,map%s2d_jcsr(ii), &
          map%s2d_wts(ii),INSERT_VALUES,ierr)
    enddo
    call MatAssemblyBegin(map%wts_mat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(  map%wts_mat,MAT_FINAL_ASSEMBLY,ierr)
  
    deallocate(index)

#ifdef MAP_DEBUG
    write(string,*) myrank
    string = 'mat_wts' // trim(adjustl(string)) // '.out'
    call PetscViewerASCIIOpen(PETSC_COMM_SELF, trim(string), viewer, ierr)
    call MatView(map%wts_mat, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

  end subroutine MappingCreateWeightMatrix

! ************************************************************************** !

  subroutine MappingCreateScatterOfSourceMesh(map,mycomm)
  ! 
  ! This routine screates a vector scatter context from source to destination
  ! grid.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    use Option_module
    implicit none
#include "petsc/finclude/petscviewer.h"

    ! argument
    type(mapping_type), pointer  :: map
    PetscMPIInt :: mycomm

    ! local variables
    Vec :: pindex       ! PETSc index
    Vec :: nindex       ! Natural index
    Vec :: N2P          ! Natural to PETSc index
    Vec :: pindex_req   ! Required PETSc indices
    IS  :: is_from, is_to
    VecScatter :: vscat
    PetscViewer :: viewer

    PetscInt                     :: ii,jj,kk
    PetscInt                     :: istart, iend
    PetscInt,pointer             :: tmp_int_array(:)
    PetscScalar,pointer          :: v_loc(:)
    PetscErrorCode               :: ierr

    !
    ! Example:
    !
    ! GIVEN:
    ! source vector (MPI):
    !                       p0        |    p1         : processor id
    !                [ s3 s4 s5 s6 s7 | s3 s1 s0 ]    : natural index
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    !
    ! required source vector (Seq):
    !                     p0        |             p1         : processor id
    !                [ s0 s4 s6 ]   | [s0 s1 s2 s3 s4 s5 s6] : natural index
    !                   0  1  2     |   0  1  2  3  4  5  6  : PETSc order
    !
    ! End product of this subroutine is construction of a vector-scatter, which
    ! has indices of MPI vector needed for creation of sequential vector.
    !
    !                     p0        |             p1         : processor id
    !                [  7  1  3 ]   | [7  6  5  0  1  2  4]  : pindex_req 

    ! Create the vectors
    call VecCreateMPI(mycomm, map%s_ncells_loc, PETSC_DECIDE, pindex, ierr)
    call VecCreateMPI(mycomm, map%s_ncells_loc, PETSC_DECIDE, nindex, ierr)
    call VecCreateMPI(mycomm, map%s_ncells_loc, PETSC_DECIDE, N2P  , ierr)
    call VecCreateMPI(mycomm, map%s2d_s_ncells_dis, PETSC_DECIDE, pindex_req, ierr)


    ! STEP-1 -
    !
    ! GIVEN:
    ! source vector (MPI):
    !                       p0        |    p1         : processor id
    !                [ s3 s4 s5 s6 s7 | s3 s1 s0 ]    : natural index
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    ! source vector (MPI) sorted in asceding order:
    !                       p0        |    p1         : processor id
    !                [ s0 s1 s2 s3 s4 | s5 s6 s7 ]    : natural index
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    ! RESULT:
    ! Indices of MPI vector required to do the sorting (N2P)
    !                       p0        |    p1         : processor id
    !                [  7  6  5  0  1 |  2  3  4  ]   : N2P
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    
    ! Initialize 'nindex' vector
    call VecGetArrayF90(nindex,v_loc,ierr)
    do ii=1,map%s_ncells_loc
       v_loc(ii) = map%s_ids_loc_nidx(ii)
    enddo
    call VecRestoreArrayF90(nindex,v_loc,ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm,'nindex.out',viewer,ierr)
    call VecView(nindex,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Initialize 'pindex' vector
    call VecGetOwnershipRange(pindex,istart,iend,ierr)
    call VecGetArrayF90(pindex,v_loc,ierr)
    do ii=1,map%s_ncells_loc
       v_loc(ii) = istart + ii - 1
    enddo
    call VecRestoreArrayF90(pindex,v_loc,ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm,'pindex.out',viewer,ierr)
    call VecView(pindex,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Create 'is_from'
    allocate(tmp_int_array(map%s_ncells_loc))
    do ii=1,map%s_ncells_loc
       tmp_int_array(ii) = istart + ii - 1
    enddo
    call ISCreateBlock(mycomm, 1, map%s_ncells_loc, tmp_int_array, PETSC_COPY_VALUES, &
         is_from, ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm,'is_from.out', &
                              viewer,ierr)
    call ISView(is_from,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Create 'is_to'
    do ii=1,map%s_ncells_loc
       tmp_int_array(ii) = map%s_ids_loc_nidx(ii)
    enddo
    call ISCreateBlock(mycomm, 1, map%s_ncells_loc, tmp_int_array, PETSC_COPY_VALUES, &
         is_to, ierr)
    deallocate(tmp_int_array)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm,'is_to.out', &
                              viewer,ierr)
    call ISView(is_to,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Create 'vscat'
    call VecScatterCreate(nindex, is_from, N2P, is_to, vscat, ierr)
    call ISDestroy(is_to,ierr)
    call ISDestroy(is_from,ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'vscat.out', viewer, ierr)
    call VecScatterView(vscat, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Scatter data
    call VecScatterBegin(vscat, pindex, N2P, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  vscat, pindex, N2P, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vscat,ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm,'N2P.out',viewer,ierr)
    call VecView(N2P,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! STEP-2 -
    !
    ! GIVEN:
    ! source vector (MPI):
    !                       p0        |    p1         : processor id
    !                [ s3 s4 s5 s6 s7 | s3 s1 s0 ]    : natural index
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    ! Indices of MPI vector required to do the sorting (N2P)
    !                       p0        |    p1         : processor id
    !                [  7  6  5  0  1 |  2  3  4  ]   : N2P
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    ! required source vector (Seq):
    !                     p0        |             p1         : processor id
    !                [ s0 s4 s6 ]   | [s0 s1 s2 s3 s4 s5 s6] : natural index
    !                   0  1  2     |   0  1  2  3  4  5  6  : PETSc order
    !
    ! RESULT:
    ! Indices of MPI vector needed for creation of sequential vector.
    !                     p0        |             p1         : processor id
    !                [  7  1  3 ]   | [7  6  5  0  1  2  4]  : pindex_req
    !

    ! Create 'is_to'
    call ISCreateBlock(mycomm, 1, map%s2d_s_ncells_dis, map%s2d_s_ids_nidx_dis, &
         PETSC_COPY_VALUES, is_to, ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'is_to1.out', viewer, ierr)
    call ISView(is_to, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Create 'is_from'
    call VecGetOwnershipRange(pindex_req, istart, iend, ierr)
    allocate(tmp_int_array(map%s2d_s_ncells_dis))
    do ii = 1,map%s2d_s_ncells_dis
       tmp_int_array(ii) = istart + ii - 1
    enddo
    call ISCreateBlock(mycomm, 1, map%s2d_s_ncells_dis, tmp_int_array, &
         PETSC_COPY_VALUES, is_from, ierr)
    deallocate(tmp_int_array)

#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm, 'is_from1.out', viewer, ierr)
    call ISView(is_from, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Create vector scatter
    call VecScatterCreate(N2P, is_to, pindex_req, is_from, vscat, ierr)
    call ISDestroy(is_to,ierr)
    call ISDestroy(is_from, ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm,'vscat.out',viewer,ierr)
    call VecScatterView(vscat,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Scatter data
    call VecScatterBegin(vscat, N2P, pindex_req, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  vscat, N2P, pindex_req, INSERT_VALUES, SCATTER_FORWARD, ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm,'pindex_req.out',viewer,ierr)
    call VecView(pindex_req,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    !
    ! Step-3 -
    ! Using 'pindex_req' create and save a vector-scatter from MPI source vector
    ! to sequential source vectors

    call VecGetArrayF90(pindex_req,v_loc,ierr)
    call ISCreateBlock(mycomm, 1, map%s2d_s_ncells_dis, INT(v_loc),&
         PETSC_COPY_VALUES, is_to, ierr)
    call VecRestoreArrayF90(pindex_req,v_loc,ierr)

    call VecScatterCreate(N2P, is_to, pindex_req, is_from, map%s2d_scat_s_gb2disloc, ierr)
#ifdef MAP_DEBUG
    call PetscViewerASCIIOpen(mycomm,'s2d_scat_s_gb2disloc.out',viewer,ierr)
    call VecScatterView(map%s2d_scat_s_gb2disloc,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Free-memory
    call VecDestroy(pindex, ierr)
    call VecDestroy(nindex, ierr)
    call VecDestroy(N2P  , ierr)
    call VecScatterDestroy(vscat,ierr)

  end subroutine MappingCreateScatterOfSourceMesh

! ************************************************************************** !

  subroutine MappingSourceToDestination(map,s_vec,d_vec)
  ! 
  ! This routine maps the data from source to destination grid.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    implicit none
    
    ! argument
    type(mapping_type), pointer :: map
    Vec                         :: s_vec ! MPI
    Vec                         :: d_vec ! Seq
    
    ! local variables
    PetscErrorCode              :: ierr
    
    if (map%s2d_s_ncells > 0) then  
       ! Initialize local vector
       call VecSet(map%s_disloc_vec, 0.d0, ierr)
    end if

    ! Scatter the source vector
    call VecScatterBegin(map%s2d_scat_s_gb2disloc, s_vec, map%s_disloc_vec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(map%s2d_scat_s_gb2disloc, s_vec, map%s_disloc_vec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr)
    
    if (map%s2d_s_ncells > 0) then  
       ! Perform Matrix-Vector product
       call MatMult(map%wts_mat, map%s_disloc_vec, d_vec, ierr)
    end if
  end subroutine

! ************************************************************************** !

function MappingListCreate()
  !
  ! This routine creates an empty map-list
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 11/09/2014
  !

  implicit none

  type(mapping_list_type), pointer :: MappingListCreate

  type(mapping_list_type), pointer :: map_list

  allocate(map_list)
  nullify(map_list%first)
  nullify(map_list%last)
  map_list%nmap = 0

  MappingListCreate => map_list

end function MappingListCreate

! ************************************************************************** !

subroutine MappingListAddToList(list, map)
  !
  ! This routine adds a map to map-list
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 11/09/2014
  !

  implicit none

  type(mapping_list_type) :: list
  type(mapping_type), pointer :: map

  if (.not. associated(list%first)) then
    list%first => map
  else
    list%last%next => map
  endif
  list%last => map

  list%nmap = list%nmap + 1

end subroutine MappingListAddToList

! ************************************************************************** !

subroutine MappingListDestroy(list)

  implicit none

  type(mapping_list_type), pointer :: list

  nullify(list%last)
  call MappingDestroy(list%first)

  deallocate(list)
  nullify(list)

end subroutine MappingListDestroy

! ************************************************************************** !

recursive subroutine MappingDestroy(map)
  !
  ! This routine frees up memoery
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/22/2014
  !

    implicit none

    ! argument
    type(mapping_type), pointer :: map

    if (.not.associated(map)) return

    if (associated(map%s_ids_loc_nidx)) deallocate(map%s_ids_loc_nidx)
    if (associated(map%d_ids_ghd_nidx)) deallocate(map%d_ids_ghd_nidx)
    if (associated(map%d_ids_nidx_sor)) deallocate(map%d_ids_nidx_sor)
    if (associated(map%d_nGhd2Sor)) deallocate(map%d_nGhd2Sor)
    if (associated(map%d_nSor2Ghd)) deallocate(map%d_nSor2Ghd)
    if (associated(map%d_loc_or_gh)) deallocate(map%d_loc_or_gh)
    if (associated(map%s2d_s_ids_nidx)) deallocate(map%s2d_s_ids_nidx)
    if (associated(map%s2d_s_ids_nidx_dis)) deallocate(map%s2d_s_ids_nidx_dis)
    if (associated(map%s2d_wts)) deallocate(map%s2d_wts)
    if (associated(map%s2d_jcsr)) deallocate(map%s2d_jcsr)
    if (associated(map%s2d_icsr)) deallocate(map%s2d_icsr)
    if (associated(map%s2d_nonzero_rcount_csr)) deallocate(map%s2d_nonzero_rcount_csr)

    nullify(map%s_ids_loc_nidx)
    nullify(map%d_ids_ghd_nidx)
    nullify(map%d_ids_nidx_sor)
    nullify(map%d_nGhd2Sor)
    nullify(map%d_nSor2Ghd)
    nullify(map%d_loc_or_gh)
    nullify(map%s2d_s_ids_nidx)
    nullify(map%s2d_s_ids_nidx_dis)
    nullify(map%s2d_wts)
    nullify(map%s2d_jcsr)
    nullify(map%s2d_icsr)
    nullify(map%s2d_nonzero_rcount_csr)

    call MappingDestroy(map%next)

  end subroutine MappingDestroy

end module Mapping_module

