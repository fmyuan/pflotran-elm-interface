module Mapping_module

  !  use clm_pflotran_interface_type

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"
#include "finclude/petscviewer.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscmat.h"

  private

  type, public :: inside_each_pflotran_cell

     PetscInt           :: num_clm_cells
     PetscInt,  pointer :: id_clm_cells(:)
     PetscReal, pointer :: perc_vol_overlap(:)
     PetscReal          :: total_vol_overlap
  end type inside_each_pflotran_cell

  type, public :: inside_each_clm_cell

     PetscInt           :: num_pflotran_cells
     PetscInt,  pointer :: id_pflotran_cells(:)
     PetscReal, pointer :: perc_vol_overlap(:)
     PetscReal          :: total_vol_overlap
  end type inside_each_clm_cell

  type, public  :: mapping_type

     type(inside_each_pflotran_cell), pointer :: pf2clm(:)
     type(inside_each_clm_cell),      pointer :: clm2pf(:)

     !
     ! Linear Mapping from Source mest to Destination mesh is matrix-vector product and
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

     ! Source mesh
     PetscInt           :: s_ncells_local      ! # of local source mesh cells present
     PetscInt,pointer   :: s_cell_ids(:)       ! IDs of source mesh cells

     ! Destination mesh
     PetscInt           :: d_ncells_local      ! # of local destination mesh cells present
     PetscInt           :: d_ncells_ghost      ! # of ghost destination mesh cells present
     PetscInt           :: d_ncells_ghosted    ! local+ghost=ghosted destination mesh cells

     ! natuaral-index starting with 0
     PetscInt,pointer   :: d_cell_ids(:)       ! IDs of ghosted destination mesh cells present
     PetscInt,pointer   :: d_cell_ids_sort(:)  ! Sorted Ghosted IDs of destination mesh cells
     PetscInt,pointer   :: d_nG2S(:)           ! Ghosted to Sorted
     PetscInt,pointer   :: d_nS2G(:)           ! Sorted to Ghosted
     PetscInt,pointer   :: d_local_or_ghost(:) ! To flag if a cell is local(1) or ghost(0)


     ! Mapping from Source-to-Destination mesh
     PetscInt           :: s2d_s_ncells               ! # of source cells for mapping
     PetscInt           :: s2d_s_ncells_distinct      ! # of "distinct" source cells for mapping
     PetscInt,pointer   :: s2d_s_cell_ids(:)          ! IDs of source cells for mapping
     PetscInt,pointer   :: s2d_s_cell_ids_distinct(:) ! IDs of "distinct" source cells for mapping

     ! Compressed Sparse Row (CSR) Matrix
     PetscReal,pointer  :: s2d_wts(:)                    ! Wts for mapping
     PetscInt,pointer   :: s2d_wts_nonzero_rcount_csr(:) ! Non-Zero entries within a row
     PetscInt,pointer   :: s2d_wts_jcsr(:)               ! J-th entry for CSR

     Mat                :: mat_wts            ! Sparse matrix for linear mapping

     VecScatter         :: s2d_scat_s_g2dl   ! Vec-Scatter of source mesh:Global to "distinct" local
     !                                          source mesh cells needed for mapping


     VecScatter         :: scatter_d_l2n     ! Vec-Scatter of destination mesh:Local to Natural order
     !

  end type mapping_type

  public :: MappingCreate, &
       MappingFindDistinctSourceMeshCellIds, &
       MappingSetDestinationMeshCellIds, &
       MappingSetSourceMeshCellIds, &
       MappingCreateScatterOfSourceMesh, &
       MappingCreateWeightMatrix, &
       MappingSourceToDestination, &
       MappingReadTxtFile, &
       MappingReadTxtFileMPI, &
       MappingDestroy

contains

  ! ************************************************************************** !
  ! 
  ! MappingCreate: Creates a Mapping Object
  ! author: Gautam Bisht
  ! date: 04/27/2011
  !
  ! ************************************************************************** !
  function MappingCreate()

    implicit none

    type(mapping_type), pointer :: MappingCreate
    type(mapping_type), pointer :: map

    allocate(map)
    
    map%filename = ''
    
    map%s_ncells_local        = 0
    map%d_ncells_local        = 0
    map%d_ncells_ghost        = 0
    map%d_ncells_ghosted      = 0
    map%s2d_s_ncells          = 0
    map%s2d_s_ncells_distinct = 0
    
    nullify(map%pf2clm)
    nullify(map%clm2pf)

    nullify(map%s_cell_ids)
    nullify(map%d_cell_ids)
    nullify(map%d_cell_ids_sort)
    nullify(map%d_nG2S)
    nullify(map%d_nS2G)
    nullify(map%d_local_or_ghost)
    nullify(map%s2d_s_cell_ids)
    nullify(map%s2d_s_cell_ids_distinct)
    nullify(map%s2d_wts)
    nullify(map%s2d_wts_nonzero_rcount_csr)
    nullify(map%s2d_wts_jcsr)

    MappingCreate => map

  end function MappingCreate

  ! ************************************************************************** !
  !
  ! MappingDestory: Deallocates any allocated pointers in Mapping object
  ! author: Gautam Bisht
  ! date: 04/27/2011
  !
  ! ************************************************************************** !
  subroutine MappingDestroy (map)

    implicit none

    type (mapping_type), pointer :: map
    PetscErrorCode               :: ierr

    if(.not.associated(map)) return

    if(associated(map%s_cell_ids)) deallocate(map%s_cell_ids)

    if(associated(map%d_cell_ids      )) deallocate(map%d_cell_ids)
    if(associated(map%d_cell_ids_sort )) deallocate(map%d_cell_ids_sort)
    if(associated(map%d_nG2S          )) deallocate(map%d_nG2S)
    if(associated(map%d_nS2G          )) deallocate(map%d_nS2G)
    if(associated(map%d_local_or_ghost)) deallocate(map%d_local_or_ghost)

    if(associated(map%s2d_s_cell_ids         )) deallocate(map%s2d_s_cell_ids)
    if(associated(map%s2d_s_cell_ids_distinct)) deallocate(map%s2d_s_cell_ids_distinct)
    if(associated(map%s2d_wts                )) deallocate(map%s2d_wts)
    if(associated(map%s2d_wts_nonzero_rcount_csr           )) deallocate(map%s2d_wts_nonzero_rcount_csr)
    if(associated(map%s2d_wts_jcsr           )) deallocate(map%s2d_wts_jcsr)

    call MatDestroy(map%mat_wts,ierr)

    call VecScatterDestroy(map%s2d_scat_s_g2dl,ierr)
    call VecScatterDestroy(map%scatter_d_l2n    ,ierr)

  end subroutine MappingDestroy

  ! ************************************************************************** !
  !
  ! MappingFindDistinctSourceMeshCellIds: Multiple cells in the destination mesh
  !   can intersect with a same cell in the source mesh. This subroutines finds
  !   distinct cells in the source mesh that overlap with destination mesh.
  !
  ! Assumption: MappingSetSourceMeshCellIds,
  !             MappingSetDestinationMeshCellsIds, and
  !             MappingReadTxtFileMPI/MappingReadTxtFile have been called 
  !             earlier
  !
  ! author: Gautam Bisht
  ! date: 04/27/2011
  !
  ! ************************************************************************** !
  subroutine MappingFindDistinctSourceMeshCellIds(map, option)

    use Option_module

    implicit none

    type (mapping_type), pointer :: map
    type(option_type), pointer   :: option
    PetscErrorCode               :: ierr
    PetscInt                     :: ii,jj,kk,count
    PetscInt, pointer            :: index(:)
    PetscInt,pointer             :: int_array(:),int_array2(:)
    PetscInt,pointer             :: int_array3(:),int_array4(:)

    Vec :: xx, yy
    
    ! No overlapped cells with Source Mesh, then return
    if(map%s2d_s_ncells.eq.0) return
    
    ! Allocate memory
    allocate(map%s2d_wts_jcsr( map%s2d_s_ncells))
    allocate(int_array (map%s2d_s_ncells))
    allocate(int_array2(map%s2d_s_ncells))
    allocate(int_array3(map%s2d_s_ncells))
    allocate(int_array4(map%s2d_s_ncells))
  
  !
  ! Follows Glenn's approach in unstructured code to remove duplicate
  ! vertices.
  !
  ! map%s2d_s_cell_ids - Contains source mesh cell ids with duplicate entries
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
      int_array(ii)  = map%s2d_s_cell_ids(ii)
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
    
    map%s2d_s_ncells_distinct = count
    ! Save the distinct ids
    allocate(map%s2d_s_cell_ids_distinct(map%s2d_s_ncells_distinct))
    
    map%s2d_s_cell_ids_distinct(1:count) = int_array3(1:count)
    map%s2d_wts_jcsr = int_array4

    ! Free memory
    deallocate(int_array)
    deallocate(int_array2)
    deallocate(int_array3)
    deallocate(int_array4)

  end subroutine MappingFindDistinctSourceMeshCellIds


  ! ************************************************************************** !
  !
  ! MappingReadTxtFile: Reads a mapping file from txt via io_rank. 
  !    - I-th processor determines the entries of Vec 'd' (refer to Eq [1] at the
  !      top) it owns, and sends them to io_rank.
  !    - io_rank determines the rows of W matrix that I-th processor would need
  !      and sends them to it.
  !
  ! ASSUMPTION:
  !    - The format of txt file to be read is assumed to be a  "Coordinate list"
  !       format i.e. (row, col, wts).
  !    - The row and col values are 1-based index.
  !    - The data is stored with ascending rows.
  !
  ! author: Gautam Bisht
  ! date: 04/27/2011
  !
  ! ************************************************************************** !
  subroutine MappingReadTxtFile(map, filename, option)

    use Input_module
    use Option_module

    implicit none

    type (mapping_type), pointer :: map
    character(len=MAXSTRINGLENGTH)     :: filename
    character(len=MAXWORDLENGTH)       :: card
    type(input_type), pointer          :: input
    type(option_type), pointer         :: option
    PetscInt                           :: num_wts
    PetscInt                           :: fileid,ii,jj,kk,irank
    PetscInt                           :: d_ncells_tmp
    PetscInt                           :: s2d_s_ncells
    PetscInt,pointer                   :: d_cells_ids_tmp(:)
    PetscInt,pointer                   :: s2d_wts_nonzero_rcount_csr(:)
    PetscInt,pointer                   :: s2d_s_cell_ids(:)
    PetscInt,pointer                   :: wts_mat_row(:), wts_mat_col(:)
    PetscReal,pointer                  :: wts_mat(:), wt_ovp_tmp(:)
    PetscErrorCode                     :: ierr
    character(len=MAXSTRINGLENGTH)     :: string  

    PetscMPIInt                        :: status_mpi(MPI_STATUS_SIZE)

    fileid   = 20
    card     = 'MppingReadTxtFile'

    if(option%myrank == option%io_rank) then
       
       ! Open the file
       input => InputCreate(fileid,filename)
       call InputReadFlotranString(input,option)
       call InputErrorMsg(input,option,'number of cells',card)
       
       ! Read the number of non-zero entries in the W matrix
       num_wts = -1
       call InputReadInt(input,option,num_wts)
       
       ! Allocate memory
       allocate(wts_mat_row(num_wts))
       allocate(wts_mat_col(num_wts))
       allocate(wts_mat(    num_wts))
       
       ! Read the entries
       do ii = 1,num_wts
          call InputReadFlotranString(input,option)
          call InputReadInt(input,option,wts_mat_row(ii))
          call InputReadInt(input,option,wts_mat_col(ii))
          call InputReadDouble(input,option,wts_mat(ii))
          
          ! Perform checks on data read
          if((wts_mat_row(ii).lt.1).or.(wts_mat_col(ii).lt.1)) then
             write(*,string),'Row/Column entry invalid for ii = ',ii
             option%io_buffer = string
             call printErrMsg(option)
          endif
          if((wts_mat(ii).lt.0d0).or.(wts_mat(ii).gt.1d0)) then
             write(*,string),'Invalid Weight value for ii = ',ii
             option%io_buffer = string
             call printErrMsg(option)
          endif
          
          ! Row/Col in the txt file are in 1-based indexing,
          ! converting it into 0-based indexing
          wts_mat_row(ii) = wts_mat_row(ii) - 1
          wts_mat_col(ii) = wts_mat_col(ii) - 1

       enddo

       do irank = 0,option%mycommsize - 1

          if(irank.ne.option%io_rank) then

             ! 0) Get from irank-th processor information regarding destination 
             !    mesh:
             !    - Number of cells present
             !    - Cell ids
             call MPI_Recv(d_ncells_tmp, 1, MPI_INTEGER, irank, MPI_ANY_TAG, &
                option%mycomm, status_mpi, ierr)

             allocate(d_cells_ids_tmp(d_ncells_tmp))
             allocate(s2d_wts_nonzero_rcount_csr(   d_ncells_tmp))

             call MPI_Recv(d_cells_ids_tmp, d_ncells_tmp, MPI_INTEGER, irank, &
                MPI_ANY_TAG, option%mycomm, status_mpi, ierr)
          else
             ! 0) Get local data
             d_ncells_tmp = map%d_ncells_ghosted
             allocate(d_cells_ids_tmp(d_ncells_tmp))
             allocate(s2d_wts_nonzero_rcount_csr(   d_ncells_tmp))
             
             do ii=1,map%d_ncells_ghosted
                d_cells_ids_tmp(ii) = map%d_cell_ids_sort(ii)
             enddo
          endif

          ! 1) Find the number of overlapped cells
          s2d_s_ncells = 0
          do ii = 1,d_ncells_tmp
             do jj = 1,num_wts
                if(d_cells_ids_tmp(ii).eq.wts_mat_row(jj)) then
                   s2d_s_ncells = s2d_s_ncells + 1
                endif
             enddo
          enddo
      
          ! 2) Save the cell-ids of overlapped cells
          s2d_wts_nonzero_rcount_csr = 0

          if (s2d_s_ncells.gt.0) then
             allocate(s2d_s_cell_ids(s2d_s_ncells))
             allocate(wt_ovp_tmp(    s2d_s_ncells))

             kk = 0
             do ii = 1,d_ncells_tmp
                do jj = 1,num_wts
                   if(d_cells_ids_tmp(ii).eq.wts_mat_row(jj)) then
                      kk                 = kk + 1
                      s2d_s_cell_ids(kk) = wts_mat_col(jj)
                      wt_ovp_tmp(kk)     = wts_mat(jj)
                      s2d_wts_nonzero_rcount_csr(ii) = &
                         s2d_wts_nonzero_rcount_csr(ii) + 1
                   endif
                enddo
             enddo
          endif

          if(irank.ne.option%io_rank) then
             ! 3) Send information back to irank-th processor
             call MPI_Send(s2d_s_ncells, 1         , MPI_INTEGER, irank, &
                option%myrank, option%mycomm, ierr)
             call MPI_Send(s2d_wts_nonzero_rcount_csr,d_ncells_tmp, MPI_INTEGER, irank, &
                option%myrank, option%mycomm, ierr)

             if(s2d_s_ncells.gt.0) then
                call MPI_Send(wt_ovp_tmp    , s2d_s_ncells, MPI_DOUBLE_PRECISION,&
                   irank, option%myrank, option%mycomm, ierr)
                call MPI_Send(s2d_s_cell_ids, s2d_s_ncells, MPI_INTEGER         ,&
                   irank, option%myrank, option%mycomm, ierr)
                do ii = 1,s2d_s_ncells
           !write(*,*), ii, s2d_s_cell_ids(ii), wt_ovp_tmp(ii)
        enddo

                ! Free memory
                deallocate(s2d_s_cell_ids)
                deallocate(wt_ovp_tmp)
             endif
          else
             ! 3) Save local information
             allocate(map%s2d_wts_nonzero_rcount_csr(d_ncells_tmp))
             do ii = 1,d_ncells_tmp
                map%s2d_wts_nonzero_rcount_csr(ii) = s2d_wts_nonzero_rcount_csr(ii)
             enddo


             map%s2d_s_ncells = s2d_s_ncells
             if (map%s2d_s_ncells.gt.0) then
                allocate(map%s2d_s_cell_ids(map%s2d_s_ncells))
                allocate(map%s2d_wts(       map%s2d_s_ncells))

                do ii = 1,map%s2d_s_ncells
                   map%s2d_s_cell_ids(ii)  = s2d_s_cell_ids(ii)
                   map%s2d_wts(ii)         = wt_ovp_tmp(ii)
           !write(*,*), ii, s2d_s_cell_ids(ii), wt_ovp_tmp(ii)
                enddo

                ! Free memory
                deallocate(s2d_s_cell_ids)
                deallocate(wt_ovp_tmp)
             endif

          endif

          ! Free memory
          deallocate(d_cells_ids_tmp)
          deallocate(s2d_wts_nonzero_rcount_csr)

       enddo

       ! Free memory
       deallocate(wts_mat_row)
       deallocate(wts_mat_col)
       deallocate(wts_mat)

    else
       !
       ! If not io_rank
       !
       call MPI_Send(map%d_ncells_ghosted, 1, MPI_INTEGER, option%io_rank, &
            option%myrank, option%mycomm, ierr)
       call MPI_Send(map%d_cell_ids_sort, map%d_ncells_ghosted,&
            MPI_INTEGER, option%io_rank, option%myrank, option%mycomm, ierr)

       call MPI_Recv(map%s2d_s_ncells, 1, MPI_INTEGER, option%io_rank, MPI_ANY_TAG, &
            option%mycomm,status_mpi, ierr)

       allocate(map%s2d_wts_nonzero_rcount_csr(map%d_ncells_ghosted))
       call MPI_Recv(map%s2d_wts_nonzero_rcount_csr, map%d_ncells_ghosted, MPI_INTEGER,&
            option%io_rank,MPI_ANY_TAG, option%mycomm, status_mpi, ierr)

       if(map%s2d_s_ncells.gt.0) then
          allocate(map%s2d_s_cell_ids(map%s2d_s_ncells))
          allocate(map%s2d_wts(       map%s2d_s_ncells))

          call MPI_Recv(map%s2d_wts       , map%s2d_s_ncells, MPI_DOUBLE_PRECISION,option%io_rank,&
               MPI_ANY_TAG, option%mycomm,status_mpi,ierr)
          call MPI_Recv(map%s2d_s_cell_ids, map%s2d_s_ncells, MPI_INTEGER         ,option%io_rank,&
               MPI_ANY_TAG, option%mycomm,status_mpi,ierr)

       endif

    endif
  
  end subroutine MappingReadTxtFile

  ! ************************************************************************** !
  !
  ! MappingReadTxtFileMPI: Reads a mapping file from txt via io_rank. The io_rank
  !    reads the file in chunks and distributes the chuncks to all processors. 
  !
  ! ASSUMPTION:
  !    - The format of txt file to be read is assumed to be a  "Coordinate list"
  !       format i.e. (row, col, wt).
  !    - The row and col values are 1-based index.
  !    - The data is stored with ascending rows.
  !
  ! author: Gautam Bisht
  ! date: 04/27/2011
  !
  ! ************************************************************************** !
  subroutine MappingReadTxtFileMPI(map, filename, option)

    use Input_module
    use Option_module

    implicit none

    type (mapping_type), pointer :: map
    character(len=MAXSTRINGLENGTH)     :: filename
    character(len=MAXWORDLENGTH)       :: card
    type(input_type), pointer          :: input
    type(option_type), pointer         :: option
    PetscInt                           :: num_wts, num_wts_loc, num_wts_tmp
    PetscInt                           :: cumsum_start,num_to_read
    PetscInt                           :: fileid,ii,jj,kk,irank,count
    PetscInt                           :: d_ncells_tmp
    PetscInt                           :: s2d_s_ncells
    PetscInt,pointer                   :: d_cells_ids_tmp(:)
    PetscInt,pointer                   :: s2d_wts_nonzero_rcount_csr(:)
    PetscInt,pointer                   :: s2d_s_cell_ids(:)
    PetscInt,pointer                   :: wts_mat_row_loc(:), wts_mat_row_tmp(:)
    PetscInt,pointer                   :: wts_mat_col_loc(:), wts_mat_col_tmp(:)
    PetscReal,pointer                  :: wts_mat_loc(:)    , wts_mat_tmp(:)
    PetscReal,pointer                  :: wt_ovp_tmp(:)
    PetscScalar,pointer                :: wts_mat_row_count(:),values(:)

    PetscViewer                        :: viewer
    PetscReal, pointer                 :: mat_values(:)
    PetscInt                           :: istart,iend
    PetscInt,pointer                   :: r_idx(:),c_idx(:),num_col(:)

    Vec :: wts_mat_row_vec, wts_mat_col_vec, wts_mat_vec
    Vec :: wts_mat_row_vec_loc, wts_mat_col_vec_loc, wts_mat_vec_loc

    Vec :: nonzero_rcount, nonzero_rcount_cumsum
    Vec :: nonzero_rcount_loc, nonzero_rcount_cumsum_loc
    Mat :: low_tri_mat

    PetscScalar,pointer          :: v_loc_1(:),v_loc_2(:),v_loc_3(:)
    PetscInt, pointer            :: tmp_int_array(:)
    IS                           :: is_from, is_to
    VecScatter                   :: vec_scat

    PetscInt                           :: remainder,buffer,row_prev
    PetscErrorCode                     :: ierr

    PetscMPIInt                        :: status_mpi(MPI_STATUS_SIZE)
    character(len=MAXSTRINGLENGTH)     :: string  
    
    PetscInt :: PRINT_RANK
    
    PRINT_RANK = -1

    row_prev = -1
    fileid   = 20
    card     = 'MppingReadTxtFileMPI'

    ! Read ASCII file through io_rank and communicate to other ranks
    if(option%myrank == option%io_rank) then

       input => InputCreate(fileid,filename)
       call InputReadFlotranString(input,option)
       call InputErrorMsg(input,option,'number of cells',card)

       num_wts = -1
       call InputReadInt(input,option,num_wts)

       num_wts_tmp = num_wts/option%mycommsize
       remainder   = num_wts - num_wts_tmp*option%mycommsize

       allocate(wts_mat_row_tmp(num_wts_tmp + 1))
       allocate(wts_mat_col_tmp(num_wts_tmp + 1))
       allocate(wts_mat_tmp(    num_wts_tmp + 1))

       do irank = 0,option%mycommsize - 1

          num_to_read = num_wts_tmp
          if(irank<remainder) num_to_read = num_to_read + 1

          ! Read the data
          do ii = 1,num_to_read
             call InputReadFlotranString(input,option)
             call InputReadInt(input,option,wts_mat_row_tmp(ii))
             call InputReadInt(input,option,wts_mat_col_tmp(ii))
             call InputReadDouble(input,option,wts_mat_tmp(ii))

            ! Perform checks on data read
            if((wts_mat_row_tmp(ii).lt.1).or.(wts_mat_col_tmp(ii).lt.1)) then
                write(*,string),'Row/Column entry invalid for ii = ',ii
                option%io_buffer = string
                call printErrMsg(option)
             endif
             if((wts_mat_tmp(ii).lt.0d0).or.(wts_mat_tmp(ii).gt.1d0)) then
                write(*,string),'Invalid Weight value for ii = ',ii
                option%io_buffer = string
                call printErrMsg(option)
             endif
             ! Ensure that input data is stored in ascending row
             if(wts_mat_row_tmp(ii).lt.row_prev) then
                write(*,string),'Data is not stored in ascending row: ii = ',ii
                option%io_buffer = string
                call printErrMsg(option)
             endif

             ! Row/Col in the txt file are in 1-based indexing,
             ! converting it into 0-based indexing
             wts_mat_row_tmp(ii) = wts_mat_row_tmp(ii) - 1
             wts_mat_col_tmp(ii) = wts_mat_col_tmp(ii) - 1

          enddo

          ! If wts reside on io_rank
          if(irank.eq.option%io_rank) then
             num_wts_loc = num_to_read
             allocate(wts_mat_row_loc(num_wts_loc))
             allocate(wts_mat_col_loc(num_wts_loc))
             allocate(wts_mat_loc(    num_wts_loc))

             do ii=1,num_wts_loc
                wts_mat_row_loc(ii) = wts_mat_row_tmp(ii)
                wts_mat_col_loc(ii) = wts_mat_col_tmp(ii)
                wts_mat_loc(ii)     = wts_mat_tmp(ii)
                !write(*,*), ii, wts_mat_row_loc(ii), wts_mat_col_loc(ii), wts_mat_loc(ii)
             enddo
          else
             ! Otherwise communicate to other ranks
             call MPI_Send(num_to_read,1,MPI_INTEGER,irank,option%myrank, &
          option%mycomm,ierr)
             call MPI_Send(wts_mat_row_tmp,num_to_read,MPI_INTEGER,irank, &
          option%myrank,option%mycomm,ierr)
             call MPI_Send(wts_mat_col_tmp,num_to_read,MPI_INTEGER,irank, &
          option%myrank,option%mycomm,ierr)
             call MPI_Send(wts_mat_tmp,num_to_read,MPI_DOUBLE_PRECISION, &
          irank,option%myrank,option%mycomm,ierr)
          endif

       enddo

       deallocate(wts_mat_row_tmp)
       deallocate(wts_mat_col_tmp)
       deallocate(wts_mat_tmp)

    else
       ! Other ranks receive data from io_rank
       call MPI_Recv(num_wts_loc,1,MPI_INTEGER,option%io_rank,MPI_ANY_TAG, &
        option%mycomm,status_mpi,ierr)
       allocate(wts_mat_row_loc(num_wts_loc))
       allocate(wts_mat_col_loc(num_wts_loc))
       allocate(wts_mat_loc(    num_wts_loc))
       call MPI_Recv(wts_mat_row_loc,num_wts_loc,MPI_INTEGER,option%io_rank, &
        MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
       call MPI_Recv(wts_mat_col_loc,num_wts_loc,MPI_INTEGER,option%io_rank, &
        MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
       call MPI_Recv(wts_mat_loc,num_wts_loc,MPI_DOUBLE_PRECISION,option%io_rank, &
        MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
    endif

    ! Save dataset related to wts in MPI-Vecs
    call VecCreateMPI(option%mycomm, num_wts_loc, PETSC_DECIDE, wts_mat_row_vec, ierr)
    call VecCreateMPI(option%mycomm, num_wts_loc, PETSC_DECIDE, wts_mat_col_vec, ierr)
    call VecCreateMPI(option%mycomm, num_wts_loc, PETSC_DECIDE, wts_mat_vec    , ierr)

    call VecGetArrayF90(wts_mat_row_vec, v_loc_1, ierr)
    call VecGetArrayF90(wts_mat_col_vec, v_loc_2, ierr)
    call VecGetArrayF90(wts_mat_vec    , v_loc_3, ierr)

    do ii = 1,num_wts_loc
       v_loc_1(ii) = wts_mat_row_loc(ii)
       v_loc_2(ii) = wts_mat_col_loc(ii)
       v_loc_3(ii) = wts_mat_loc(ii)
    enddo

    call VecRestoreArrayF90(wts_mat_row_vec, v_loc_1, ierr)
    call VecRestoreArrayF90(wts_mat_col_vec, v_loc_2, ierr)
    call VecRestoreArrayF90(wts_mat_vec    , v_loc_3, ierr)

    !
    ! After reading the Mapping file in parallel, reshuffling of the data
    ! needs to occur. 
    !
    !        W_loc * s_loc  = d_loc
    !        W_loc * s_dloc = d_loc
    !
    ! For all points of destination mesh that are active on a given
    ! processor, find:
    !
    ! - The points of source meshes that are required
    ! -
    !


    ! For each cell of destination mesh, find the count of cells overlaped in
    ! source mesh.
    !                              OR
    ! Num of non-zero entries in each row of the global W matrix

    call VecCreateMPI(option%mycomm, map%d_ncells_local, PETSC_DECIDE, nonzero_rcount       , ierr)

    call VecGetSize(nonzero_rcount,ii,ierr)

    allocate(wts_mat_row_tmp(num_wts_loc))
    allocate(wts_mat_row_count(num_wts_loc))
    wts_mat_row_tmp   = 0
    wts_mat_row_count = 0

    jj = 1
    ii = 1
    wts_mat_row_tmp(jj)   = wts_mat_row_loc(ii)
    wts_mat_row_count(jj) = 1

    do ii = 2,num_wts_loc
       if( wts_mat_row_tmp(jj).eq.wts_mat_row_loc(ii)) then
          wts_mat_row_count(jj) = wts_mat_row_count(jj) + 1
       else
          jj = jj + 1
          wts_mat_row_tmp(jj)   = wts_mat_row_loc(ii)
          wts_mat_row_count(jj) = 1
       endif
    enddo
    
    if(option%myrank.eq.PRINT_RANK) then
      write(*,*), 'wts_mat_row_count: '
      do ii = 1,jj
        write(*,*), ii, wts_mat_row_tmp(ii),wts_mat_row_count(ii)
      enddo
    endif

    call VecSetValues(nonzero_rcount, jj, wts_mat_row_tmp, wts_mat_row_count, ADD_VALUES, ierr);
    call VecAssemblyBegin(nonzero_rcount,ierr)
    call VecAssemblyEnd(nonzero_rcount,ierr)
    deallocate(wts_mat_row_tmp)
    deallocate(wts_mat_row_count)

    ! Find cummulative sum of the nonzero_rcount vector
    call VecCreateMPI(option%mycomm, map%d_ncells_local, PETSC_DECIDE, nonzero_rcount_cumsum, ierr)
    call VecGetArrayF90(nonzero_rcount       ,v_loc_1,ierr)
    call VecGetArrayF90(nonzero_rcount_cumsum,v_loc_2,ierr)

    ii = 1
    v_loc_2(ii) = v_loc_1(ii)
    if(option%myrank.eq.PRINT_RANK) write (*,*), 'v_loc_1: ',ii, v_loc_1(ii), v_loc_2(ii)
    do ii = 2,map%d_ncells_local
       v_loc_2(ii) = v_loc_1(ii) + v_loc_2(ii-1)
       if(option%myrank.eq.PRINT_RANK) write (*,*), 'v_loc_1: ',ii, v_loc_1(ii), v_loc_2(ii)
    enddo

    cumsum_start = 0
    call MPI_Exscan(INT(v_loc_2(map%d_ncells_local)),cumsum_start,ONE_INTEGER_MPI,MPIU_INTEGER, &
         MPI_SUM,option%mycomm,ierr)

    if(option%myrank.eq.PRINT_RANK) write (*,*), 'cum_start: ',cumsum_start
    do ii = 1,map%d_ncells_local
       v_loc_2(ii) = v_loc_2(ii) + cumsum_start
       if(option%myrank.eq.PRINT_RANK) write (*,*), 'v_loc_2: ',ii, v_loc_2(ii)
    enddo

    call VecRestoreArrayF90(nonzero_rcount       , v_loc_1,ierr)
    call VecRestoreArrayF90(nonzero_rcount_cumsum, v_loc_2,ierr)

    ! On a given processor, find the number of source mesh cells required to 
    ! reconstruct variables for the destination mesh.
    !
    ! From the MPI-Vec nonzero_rcount, get entries corresponding to destination
    ! cells present on a given processor
    !
    
    ! Create local vectors to hold the scattered data
    call VecCreateSeq(PETSC_COMM_SELF, map%d_ncells_ghosted, &
                      nonzero_rcount_loc, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, map%d_ncells_ghosted, &
                      nonzero_rcount_cumsum_loc, ierr)

    ! Allocate and initialize a temporary array to be used for creation of index
    ! set
    allocate(tmp_int_array(map%d_ncells_ghosted))
    do ii=1,map%d_ncells_ghosted
       tmp_int_array(ii) = ii-1
    enddo

    ! Create an index set to scatter to local vector
    call ISCreateBlock(option%mycomm, 1, map%d_ncells_ghosted, tmp_int_array, &
                       PETSC_COPY_VALUES, is_to, ierr)
    deallocate(tmp_int_array)

    ! Create an index set to scatter from MPI vector
    call ISCreateBlock(option%mycomm, 1, map%d_ncells_ghosted, &
                      map%d_cell_ids_sort, PETSC_COPY_VALUES, is_from, ierr)

    if(option%myrank.eq.PRINT_RANK) then
      write(*,*), 'd_cell_ids_sort: ', option%myrank
      do ii = 1,map%d_ncells_ghosted
        write(*,*), ii, map%d_cell_ids_sort(ii)
      enddo
    endif

    ! Create a vector scatter context
    call VecScatterCreate(nonzero_rcount, is_from, nonzero_rcount_loc, is_to, &
                          vec_scat, ierr)
    call ISDestroy(is_from, ierr)
    call ISDestroy(is_to  , ierr)

    ! Scatter the data (i) nonzero_rcount; (ii) nonzero_rcount_cumsum
    call VecScatterBegin(vec_scat, nonzero_rcount, nonzero_rcount_loc, &
                         INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(vec_scat, nonzero_rcount, nonzero_rcount_loc, &
                       INSERT_VALUES, SCATTER_FORWARD, ierr)

    call VecScatterBegin(vec_scat, nonzero_rcount_cumsum, &
                         nonzero_rcount_cumsum_loc, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(vec_scat, nonzero_rcount_cumsum, &
                       nonzero_rcount_cumsum_loc, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)

    ! Destroy vector scatter
    call VecScatterDestroy(vec_scat, ierr)

    call VecGetArrayF90(nonzero_rcount_loc       , v_loc_1, ierr)
    call VecGetArrayF90(nonzero_rcount_cumsum_loc, v_loc_2, ierr)

    allocate(map%s2d_wts_nonzero_rcount_csr(map%d_ncells_ghosted))

    ! For each destination, save the number of overlapped source mesh cells
    count = 0 ! cummulative count of number of overlapped source mesh cells
    if(option%myrank.eq.PRINT_RANK) write(*,*), 'local rcount and rcount_cumsum: rank ',option%myrank
    do ii=1,map%d_ncells_ghosted
       count = count + INT(v_loc_1(ii))
       map%s2d_wts_nonzero_rcount_csr(ii) = INT(v_loc_1(ii))
       if(option%myrank.eq.PRINT_RANK) write(*,*), ii, count, map%s2d_wts_nonzero_rcount_csr(ii)
    enddo

    map%s2d_s_ncells = count

    !if (map%s2d_s_ncells > 0) then
      
      ! Allocate memory to save cell ids of source mesh that overlap with 
      ! cells in destination mesh and corresponding weights
      allocate(map%s2d_s_cell_ids(map%s2d_s_ncells))
      allocate(map%s2d_wts(       map%s2d_s_ncells))

      ! Allocate memory
      allocate(tmp_int_array(map%s2d_s_ncells))

      ! For each cell in destination mesh, save indices of MPI Vectors, which
      ! contain data read from mapping file, for all overlapped cells of 
      ! of source mesh
      if(option%myrank.eq.PRINT_RANK) write(*,*), 'tmp_int_array: rank ',option%myrank      
      kk = 0
      do ii = 1,map%d_ncells_ghosted
        do jj = 1,INT(v_loc_1(ii))
          kk = kk + 1
          tmp_int_array(kk) = INT(v_loc_2(ii)) - INT(v_loc_1(ii)) + jj - 1
          if(option%myrank.eq.PRINT_RANK) write(*,*), kk, tmp_int_array(kk), map%d_cell_ids_sort(ii)
        enddo
      enddo

      ! Create an index set to scatter from
      call ISCreateBlock(option%mycomm, 1, map%s2d_s_ncells, tmp_int_array, &
                         PETSC_COPY_VALUES, is_from, ierr)

      do ii=1,map%s2d_s_ncells
        tmp_int_array(ii) = ii-1
      enddo

      call ISCreateBlock(option%mycomm, 1, map%s2d_s_ncells, tmp_int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
      deallocate(tmp_int_array)

      ! Allocate memory
      call VecCreateSeq(PETSC_COMM_SELF, map%s2d_s_ncells, wts_mat_row_vec_loc, &
                        ierr)
      call VecCreateSeq(PETSC_COMM_SELF, map%s2d_s_ncells, wts_mat_col_vec_loc, &
                        ierr)
      call VecCreateSeq(PETSC_COMM_SELF, map%s2d_s_ncells, wts_mat_vec_loc, ierr)

      ! Create scatter context
      call VecScatterCreate(wts_mat_row_vec, is_from, wts_mat_row_vec_loc, &
                            is_to, vec_scat, ierr)

      ! Scatter the data
      call VecScatterBegin(vec_scat, wts_mat_col_vec, wts_mat_col_vec_loc, &
                           INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(vec_scat, wts_mat_col_vec, wts_mat_col_vec_loc, &
                         INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterBegin(vec_scat, wts_mat_vec, wts_mat_vec_loc, &
                           INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(vec_scat, wts_mat_vec, wts_mat_vec_loc, &
                         INSERT_VALUES, SCATTER_FORWARD, ierr)
      
      ! Attach to the local copy of the scatterd data
      call VecGetArrayF90(wts_mat_col_vec_loc,v_loc_1,ierr)
      call VecGetArrayF90(wts_mat_vec_loc    ,v_loc_2,ierr)

      ! Save the scattered data
      if(option%myrank.eq.PRINT_RANK) write(*,*), 'd_s_cell_ids: rank ',option%myrank
      do ii = 1,map%s2d_s_ncells
        map%s2d_s_cell_ids(ii) = INT(v_loc_1(ii))
        map%s2d_wts(ii)        = v_loc_2(ii)
        if(option%myrank.eq.PRINT_RANK) write(*,*), ii, map%s2d_s_cell_ids(ii)
      enddo 
      
      ! Restore vectors
      call VecRestoreArrayF90(wts_mat_col_vec_loc, v_loc_1, ierr)
      call VecRestoreArrayF90(wts_mat_vec_loc    , v_loc_2, ierr)
      
      ! Free memory
      call VecDestroy(wts_mat_col_vec_loc, ierr)
      call VecDestroy(wts_mat_vec_loc, ierr)

    !endif

    ! Restore vectors
    call VecRestoreArrayF90(nonzero_rcount_cumsum_loc, v_loc_1, ierr)
    call VecRestoreArrayF90(nonzero_rcount_loc       , v_loc_2, ierr)


    ! Free Memory
    call VecDestroy(nonzero_rcount_loc,ierr)
    call VecDestroy(nonzero_rcount_cumsum_loc,ierr)
    call VecDestroy(nonzero_rcount,ierr)
    call VecDestroy(nonzero_rcount_cumsum,ierr)
    call VecDestroy(wts_mat_row_vec,ierr)
    call VecDestroy(wts_mat_col_vec,ierr)
    call VecDestroy(wts_mat_vec,ierr)

  end subroutine MappingReadTxtFileMPI

  ! ************************************************************************** !
  !
  ! MappingSetSourceMeshCellIds: Sets the cell-ids of the source mesh owned by
  !   the processor.
  !
  ! author: Gautam Bisht
  ! date: 04/27/2011
  !
  ! ************************************************************************** !
  subroutine MappingSetSourceMeshCellIds(map, option, num_cells, cell_ids)

    use Option_module
    implicit none

    type(mapping_type), pointer  :: map
    type(option_type), pointer   :: option

    PetscInt                     :: num_cells, ii
    PetscInt, intent(in),pointer :: cell_ids(:)

    map%s_ncells_local = num_cells
    allocate(map%s_cell_ids(num_cells))

    do ii = 1,num_cells
       map%s_cell_ids(ii) = cell_ids(ii)
    enddo

  end subroutine MappingSetSourceMeshCellIds


  ! ************************************************************************** !
  !
  ! MappingSetDestinationMeshCellIds: 
  !   - Sets the cell-ids of the destination mesh owned by the processor.
  !   - Also, sets a flag to determine if a cell is local or ghost.
  !   - Sorts the destination mesh cell-ids in ascending order.
  !   - Saves indicies for: Ghosted-to-Sort and Sort-to-Ghosted order.
  !
  ! author: Gautam Bisht
  ! date: 04/27/2011
  !
  ! ************************************************************************** !
  subroutine MappingSetDestinationMeshCellIds(map, option, num_cells_local, &
       num_cells_ghost, cell_ids_ghosted, local_or_ghost)

    use Option_module
    implicit none

    type(mapping_type), pointer  :: map
    type(option_type), pointer   :: option

    PetscInt                     :: num_cells_local, num_cells_ghost, ii
    PetscInt, intent(in),pointer :: cell_ids_ghosted(:), local_or_ghost(:)
    PetscInt, pointer            :: index(:), rev_index(:)
    PetscInt,pointer             :: tmp_int_array(:)
    PetscErrorCode               :: ierr
    IS                           :: is_from, is_to
    Vec                          :: vec_tmp_1, vec_tmp_2
    PetscViewer :: viewer
    PetscInt                     :: istart, iend,count
  
  ! Initialize
    map%d_ncells_local   = num_cells_local
    map%d_ncells_ghost   = num_cells_ghost
    map%d_ncells_ghosted = num_cells_local + num_cells_ghost
  
  ! Assign memory
    allocate(map%d_cell_ids(        map%d_ncells_ghosted))
    allocate(map%d_cell_ids_sort(   map%d_ncells_ghosted))
    allocate(map%d_local_or_ghost(  map%d_ncells_ghosted))
    allocate(map%d_nG2S(            map%d_ncells_ghosted))
    allocate(map%d_nS2G(            map%d_ncells_ghosted))
    allocate(index(                 map%d_ncells_ghosted))
    allocate(rev_index(             map%d_ncells_ghosted))
  
  ! Save cell-ids, local/ghost flag, and initialize index and rev-index
    do ii = 1,num_cells_local+num_cells_ghost
       map%d_cell_ids(ii)       = cell_ids_ghosted(ii)
       map%d_local_or_ghost(ii) = local_or_ghost(ii)
       index(ii)     = ii
       rev_index(ii) = ii
    enddo

    ! Sort the d_cell_ids
    index = index - 1
    call PetscSortIntWithPermutation( map%d_ncells_ghosted, cell_ids_ghosted, index, ierr)
    index = index + 1
  
  ! Save the d_cells_ids in ascending sort order. Also save the index of
  ! Ghosted-to-Sorted order 
    do ii = 1,num_cells_local+num_cells_ghost
       map%d_cell_ids_sort(ii) = cell_ids_ghosted(index(ii))
       map%d_nG2S(ii)          = index(ii)
    enddo

    ! Sort the index
    rev_index = rev_index - 1
    call PetscSortIntWithPermutation( map%d_ncells_ghosted, index, rev_index, ierr)
    rev_index = rev_index + 1
  
  ! Save the reverse sort index: Sorted-to-Ghosted
    do ii = 1,num_cells_local+num_cells_ghost
       map%d_nS2G(ii)             = rev_index(ii)
    enddo
  
  ! Create vectors
    !write (*,*), 'map%d_ncells_ghosted = ', map%d_ncells_ghosted
    !write (*,*), 'map%d_ncells_local   = ', map%d_ncells_local
    call VecCreateMPI(option%mycomm, map%d_ncells_ghosted, PETSC_DECIDE, vec_tmp_1, ierr)
    call VecCreateMPI(option%mycomm, map%d_ncells_local  , PETSC_DECIDE, vec_tmp_2, ierr)
    ! 
    call VecGetOwnershipRange(vec_tmp_1,istart,iend,ierr)
    allocate(tmp_int_array(map%d_ncells_local))
    count = 0
    do ii=1,map%d_ncells_ghosted
       if(map%d_local_or_ghost(ii).eq.1) then
          count = count + 1
          tmp_int_array(count) = istart+ii-1
       endif
    enddo
    call ISCreateBlock(option%mycomm,1,map%d_ncells_local, tmp_int_array, PETSC_COPY_VALUES, &
         is_from, ierr)
    deallocate(tmp_int_array)

    allocate(tmp_int_array(map%d_ncells_local))
    count = 0
    do ii=1,map%d_ncells_ghosted
     if(map%d_local_or_ghost(ii).eq.1) then
        count = count + 1
          tmp_int_array(count) = map%d_cell_ids(ii)
    endif
    enddo
    
    call ISCreateBlock(option%mycomm,1,map%d_ncells_local, tmp_int_array, PETSC_COPY_VALUES, &
         is_to, ierr)
    deallocate(tmp_int_array)
  
    !call PetscViewerASCIIOpen(option%mycomm, 'is_from_1.out', viewer, ierr)
    !call ISView(is_from, viewer,ierr)
    !call PetscViewerDestroy(viewer, ierr)

    !call PetscViewerASCIIOpen(option%mycomm, 'is_to_1.out', viewer, ierr)
    !call ISView(is_to, viewer,ierr)
    !call PetscViewerDestroy(viewer, ierr)
  
    !write(*,*), 'call VecScatterCreate()'
    call VecScatterCreate(vec_tmp_1, is_from, vec_tmp_2, is_to, map%scatter_d_l2n, ierr)

    !call PetscViewerASCIIOpen(option%mycomm, 'scatter_d_l2n.out', viewer, ierr)
    !call VecScatterView(map%scatter_d_l2n, viewer,ierr)
    !call PetscViewerDestroy(viewer, ierr)

    ! Free memory
    !deallocate(index)
    !deallocate(rev_index)
    !call VecDestroy(vec_tmp_1)
    !call VecDestroy(vec_tmp_2)
    !call ISDestroy(is_from)
    !call ISDestroy(is_to)

    !write(*,*), 'done in MappingSetDestinationMeshCellIds....'
  end subroutine MappingSetDestinationMeshCellIds


  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingCreateScatterOfSourceMesh (map, option)

    use Option_module
    implicit none

    type(mapping_type), pointer  :: map
    type(option_type), pointer   :: option

    Vec :: porder ! PETSc Order
    Vec :: aorder ! Application/Natural Order
    Vec :: nA2P   ! Application to PETSc Order
    Vec :: porder_req
    IS  :: is_ltol, is_ltog, is_gtol
    VecScatter :: vscat_1,vscat_2
    PetscViewer :: viewer

    PetscInt                     :: ii,jj,kk
    PetscInt                     :: istart, iend
    PetscInt,pointer             :: tmp_int_array(:)
    PetscScalar,pointer          :: v_loc(:)
    PetscErrorCode               :: ierr

    ! Create the vectors
    call VecCreateMPI(option%mycomm, map%s_ncells_local, PETSC_DECIDE, porder, ierr)
    call VecCreateMPI(option%mycomm, map%s_ncells_local, PETSC_DECIDE, aorder, ierr)
    call VecCreateMPI(option%mycomm, map%s_ncells_local, PETSC_DECIDE, nA2P  , ierr)
    call VecCreateMPI(option%mycomm, map%s2d_s_ncells_distinct, PETSC_DECIDE, porder_req, ierr)

    ! Initialize 'aorder' vector
    call VecGetArrayF90(aorder,v_loc,ierr)
    do ii=1,map%s_ncells_local
       v_loc(ii) = map%s_cell_ids(ii)
    enddo
    call VecRestoreArrayF90(aorder,v_loc,ierr)

    ! Initialize 'porder' vector
    call VecGetOwnershipRange(porder,istart,iend,ierr)
    call VecGetArrayF90(porder,v_loc,ierr)
    do ii=1,map%s_ncells_local
       v_loc(ii) = istart + ii - 1
    enddo
    call VecRestoreArrayF90(porder,v_loc,ierr)

    ! Create 'is_ltol'
    allocate(tmp_int_array(map%s_ncells_local))
    do ii=1,map%s_ncells_local
       tmp_int_array(ii) = istart + ii - 1
    enddo
    call ISCreateBlock(option%mycomm, 1, map%s_ncells_local, tmp_int_array, PETSC_COPY_VALUES, &
         is_ltol, ierr)

    ! Create 'is_ltog'
    do ii=1,map%s_ncells_local
       tmp_int_array(ii) = map%s_cell_ids(ii)
    enddo
    call ISCreateBlock(option%mycomm, 1, map%s_ncells_local, tmp_int_array, PETSC_COPY_VALUES, &
         is_ltog, ierr)
    deallocate(tmp_int_array)

    ! Create 'vscat'
    call VecScatterCreate(aorder, is_ltol, nA2P, is_ltog, vscat_1, ierr)
    call ISDestroy(is_ltog,ierr)
    call ISDestroy(is_ltol,ierr)
  
    !call PetscViewerASCIIOpen(option%mycomm, 'vscat_1.out', viewer, ierr)
    !call VecScatterView(vscat_1, viewer,ierr)
    !call PetscViewerDestroy(viewer, ierr)

    ! Scatter data
    call VecScatterBegin(vscat_1, porder, nA2P, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  vscat_1, porder, nA2P, INSERT_VALUES, SCATTER_FORWARD, ierr)

    !
    call ISCreateBlock(option%mycomm, 1, map%s2d_s_ncells_distinct, map%s2d_s_cell_ids_distinct, &
         PETSC_COPY_VALUES, is_gtol, ierr)

    call VecGetOwnershipRange(porder_req, istart, iend, ierr)

    allocate(tmp_int_array(map%s2d_s_ncells_distinct))
    do ii = 1,map%s2d_s_ncells_distinct
       tmp_int_array(ii) = istart + ii - 1
    enddo
    call ISCreateBlock(option%mycomm, 1, map%s2d_s_ncells_distinct, tmp_int_array, &
         PETSC_COPY_VALUES, is_ltol, ierr)

    !call PetscViewerASCIIOpen(option%mycomm, 'is_ltol.out', viewer, ierr)
    !call ISView(is_ltol, viewer,ierr)
    !call PetscViewerDestroy(viewer, ierr)

    !call PetscViewerASCIIOpen(option%mycomm, 'is_gtol.out', viewer, ierr)
    !call ISView(is_gtol, viewer,ierr)
    !call PetscViewerDestroy(viewer, ierr)

    !if(option%myrank == option%io_rank) write(*,*), 'attempting to create vscat_2 ....'
    call VecScatterCreate(nA2P, is_gtol, porder_req, is_ltol, vscat_2, ierr)
    !if(option%myrank == option%io_rank) write(*,*), 'vscat_2 created'

    call VecScatterBegin(vscat_2, nA2P, porder_req, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  vscat_2, nA2P, porder_req, INSERT_VALUES, SCATTER_FORWARD, ierr)

    call ISDestroy(is_gtol,ierr)
    call VecGetArrayF90(porder_req,v_loc,ierr)
    call ISCreateBlock(option%mycomm, 1, map%s2d_s_ncells_distinct, INT(v_loc),&
         PETSC_COPY_VALUES, is_gtol, ierr)
    call VecRestoreArrayF90(porder_req,v_loc,ierr)


    call VecScatterCreate(nA2P, is_gtol, porder_req, is_ltol, map%s2d_scat_s_g2dl, ierr)

    ! Free-memory
    call VecDestroy(porder, ierr)
    call VecDestroy(aorder, ierr)
    call VecDestroy(nA2P  , ierr)
    call VecScatterDestroy(vscat_1,ierr)
    call VecScatterDestroy(vscat_2,ierr)
    !call ISDestroy(is_gtol,ierr)
    !call ISDestroy(is_ltol, ierr)


  end subroutine MappingCreateScatterOfSourceMesh


  ! ************************************************************************** !
  !
  ! MappingCreateWeightMatrix: Creates a sequential weight matrix
  ! author: Gautam Bisht
  ! date: 04/27/2011
  !
  ! ************************************************************************** !
  subroutine MappingCreateWeightMatrix(map, option)

    use Option_module
    implicit none

    type(mapping_type), pointer  :: map
    type(option_type), pointer   :: option
    PetscInt, pointer            :: index(:)
    PetscInt                     :: ii,jj,kk
    PetscErrorCode               :: ierr
    character(len=MAXSTRINGLENGTH)     :: string  
    PetscViewer :: viewer

    allocate(index(map%s2d_s_ncells))

    kk = 0
    do ii = 1,map%d_ncells_ghosted
       do jj = 1,map%s2d_wts_nonzero_rcount_csr(ii)
          kk = kk + 1
          index(kk) = ii -1
       enddo
    enddo

    !
    ! size(mat_wts) = [d_ncells_ghosted x s2d_s_ncells_distinct]
    !
    call MatCreateSeqAIJ(PETSC_COMM_SELF, map%d_ncells_ghosted, &
         map%s2d_s_ncells_distinct, PETSC_NULL, &
         map%s2d_wts_nonzero_rcount_csr, map%mat_wts, ierr)

    do ii = 1,map%s2d_s_ncells
       call MatSetValues(map%mat_wts,1,index(ii),1,map%s2d_wts_jcsr(ii), &
          map%s2d_wts(ii),INSERT_VALUES,ierr)
    enddo
    call MatAssemblyBegin(map%mat_wts,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(  map%mat_wts,MAT_FINAL_ASSEMBLY,ierr)
  
    deallocate(index)
    
    if(option%myrank.eq.0) then
       !call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'mat_wts_0.out', viewer, ierr)
    else
       !call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'mat_wts_1.out', viewer, ierr)
    endif
    !call MatView(map%mat_wts, viewer,ierr)
    !call PetscViewerDestroy(viewer, ierr)
  
  
  end subroutine MappingCreateWeightMatrix

  ! ************************************************************************** !
  !
  ! MappingSourceToDestination: 
  ! author: Gautam Bisht
  ! date: 04/27/2011
  !
  ! ************************************************************************** !
  subroutine MappingSourceToDestination(map, option, vec_s, vec_d)

    use Option_module
    implicit none

    type(mapping_type), pointer  :: map
    type(option_type), pointer   :: option
    Vec :: vec_s      , vec_d
    Vec :: vec_s_local, vec_d_local
    Vec :: vec_s_loc  , vec_d_loc
    Vec :: vec_d_nat
    PetscScalar,pointer :: v_loc_1(:), v_loc_2(:)
    PetscInt :: ii, rsize, csize
    PetscErrorCode               :: ierr
    PetscViewer :: viewer

    call VecCreateMPI(option%mycomm  , map%s2d_s_ncells_distinct, PETSC_DECIDE, vec_s_local, ierr)
    call VecCreateMPI(option%mycomm  , map%d_ncells_local       , PETSC_DECIDE, vec_d_nat, ierr)

    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_s_ncells_distinct, vec_s_loc, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, map%d_ncells_ghosted     , vec_d_loc, ierr)

    call VecScatterBegin(map%s2d_scat_s_g2dl, vec_s, vec_s_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  map%s2d_scat_s_g2dl, vec_s, vec_s_local, INSERT_VALUES, SCATTER_FORWARD, ierr)

    call VecGetArrayF90(vec_s_local,v_loc_1,ierr)
    call VecGetArrayF90(vec_s_loc  ,v_loc_2,ierr)
    v_loc_2 = v_loc_1
    call VecRestoreArrayF90(vec_s_local,v_loc_1,ierr)
    call VecRestoreArrayF90(vec_s_loc  ,v_loc_2,ierr)

    call MatMult(map%mat_wts, vec_s_loc, vec_d_loc, ierr)


    call VecGetArrayF90(vec_d    ,v_loc_1,ierr)
    call VecGetArrayF90(vec_d_loc,v_loc_2,ierr)
    v_loc_1 = v_loc_2
    call VecRestoreArrayF90(vec_d    ,v_loc_1,ierr)
    call VecRestoreArrayF90(vec_d_loc,v_loc_2,ierr)
  
    call VecScatterBegin(map%scatter_d_l2n, vec_d, vec_d_nat, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  map%scatter_d_l2n, vec_d, vec_d_nat, INSERT_VALUES, SCATTER_FORWARD, ierr)

    !call PetscViewerASCIIOpen(option%mycomm, 'u_b.out', viewer, ierr)
    !call VecView(vec_d_nat, viewer,ierr)
    !call PetscViewerDestroy(viewer, ierr)

    !if(option%myrank.eq.0) then
    !   call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'vec_d_loc_0.out', viewer, ierr)
    !else
    !   call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'vec_d_loc_1.out', viewer, ierr)
    !endif
    !call VecView(vec_d_loc, viewer,ierr)
    !call PetscViewerDestroy(viewer, ierr)
    
    !
    call VecDestroy(vec_s_local,ierr)
    call VecDestroy(vec_s_loc,ierr)
    call VecDestroy(vec_d_loc,ierr)
    call VecDestroy(vec_d_nat,ierr)
    
    

  end subroutine MappingSourceToDestination


end module Mapping_module
