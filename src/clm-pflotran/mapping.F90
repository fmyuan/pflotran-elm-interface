module Mapping_module

  !  use clm_pflotran_interface_type

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"
!#include "finclude/petscsys.h"
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

     ! This module does the mapping between two meshs:
     !  'fmesh'   - The mesh from which data is to be mapped.
     !  'tmesh'   - The mesh to which data is being mapped.
     !  'ocell'   - Cells of 'fmesh' which overlapped with 'tmesh'
     !  'docell'  - Since a cell of 'fmesh' can overlap with multiple cells
     !              of 'tmesh', this refers to distinct (i.e. without duplications)
     !              overlapped cells of 'fmesh'

     PetscInt :: map_id                            !
#if 0
     ! Size of 'tmesh'
     PetscInt :: tmesh_num_cells_local                    ! Number of local cells
     PetscInt :: tmesh_num_cells_ghost                    ! Number of ghost cells
     PetscInt :: tmesh_num_cells_ghosted                  ! num_cells_local + num_cells_ghost

     PetscInt :: num_ocells_with_fmesh                    ! Total number of cells 'fmesh'
     !                                                      which overlap with 'tmesh'
     !PetscInt :: num_docells_with_fmesh                   ! Number of distinct cells 'fmesh'
     !                                                      which overlap with 'tmesh'

     ! 'tmesh' Cells:
     PetscInt, pointer  :: tmesh_cell_ids_ghosted_nindex(:)! Cell-IDs in Natural indexing             [tmesh_num_cells_ghosted]


     !                                                     'fmesh' Overlapped Cells:
     PetscInt, pointer  :: ocell_cnt_ghosted(:)            ! For each 'tmesh' cell, the count of
     !                                                       cells within 'fmesh' overlapped with     [tmesh_num_cells_ghosted]
     PetscInt, pointer  :: ocell_cnt_cumsum_ghosted(:)     ! Cummulative sum of 'ocell_cnt_ghosted'   [tmesh_num_cells_ghosted]
     PetscInt, pointer  :: ocell_ids_nindex (:)            ! Overlapped Cell-IDs in Natural indexing  [num_ocells_with_fmesh]
     PetscReal, pointer :: ocell_vol_nindex(:)             ! Overlapped volume in Natural indexing    [num_ocells_with_fmesh]
     
     Mat :: weights_f2tmesh
	 
     PetscInt, pointer  :: docell_ids_sorted_nindex(:)     ! Distinct overlapped Cell-IDs sorted and
     !                                                       in natural indexing                      [num_docells_with_fmesh]
     PetscInt, pointer   :: hash(:)                        ! Lookup table for obtaining position in
     !                                                       the sorted distinct cell ids
     !                                                       (docell_ids_sorted) list for every
     !                                                       overlapped cell-id (ocell)               [num_ocells_with_fmesh]

     !
     ! This class allows 'scatter' of local data on fmesh to a
     ! global vector; from which local data is 'sactter' for
     ! tmesh
     !
     !
     !               [scatter]             [scatter]
     ! fmeshlocal --------------> global ------------> tmeshlocal
     !
     !

     ! Index set
     IS         :: is_fmeshlocal_to_global
     IS         :: is_global_to_tmeshlocal

     ! Vector scatter
     VecScatter :: scatter_fmeshlocal_to_global
     VecScatter :: scatter_global_to_tmeshlocal
#endif




     ! Source mesh
	 PetscInt           :: s_ncells_local
	 PetscInt,pointer	:: s_cell_ids(:)
	 
	 ! Destination mesh
	 PetscInt           :: d_ncells_local
	 PetscInt           :: d_ncells_ghost
	 PetscInt           :: d_ncells_ghosted
	 PetscInt,pointer   :: d_cell_ids(:)
	 PetscInt,pointer   :: d_cell_ids_sort(:)
	 PetscInt,pointer   :: d_nG2S(:)          ! Ghosted to Sorted
	 PetscInt,pointer   :: d_nS2G(:)          ! Sorted to Ghosted
	 PetscInt,pointer   :: d_local_or_ghost(:)
	 
	 
	 !
	 PetscInt           :: s2d_s_ncells       ! Num of source cells required for mapping
     PetscInt           :: s2d_s_ncells_distinct
	 PetscInt,pointer   :: s2d_s_cell_ids(:)  ! Ids of source cells required for mapping
     PetscInt,pointer   :: s2d_s_cell_ids_distinct(:)
	 PetscReal,pointer  :: s2d_wts(:)         ! Wts for mapping
	 PetscInt,pointer   :: s2d_wts_icsr(:)
	 PetscInt,pointer   :: s2d_wts_jcsr(:)    !
	 !PetscInt,pointer   :: s2d_nonzerowt_count_for_each_d(:) ! s2d_wts_icsr()

     VecScatter         :: s2d_vscat
	 
	 Mat                :: mat_wts

     VecScatter         :: d_local_to_natural_vscat













  end type mapping_type

  public :: MappingCreate, &
       MappingAllocateMemory, &
       MappingCreateIS, &
       MappingCreateVecScatter, &
	   MappingFindDistinctSourceMeshCellIds, &
       MappingScatterLocal2Global, &
       MappingScatterGlobal2Local, &
       MappingSetOcells, &
	   MappingSetNumCells, &
	   MappingSetTmeshCellIds, &
	   MappingSetDestinationMeshCellIds, &
	   MappingSetSourceMeshCellIds, &
	   MappingCreateScatterOfSourceMesh, &
	   MappingCreateWeightMatrix, &
	   MappingSourceToDestination, &
	   MappingReadTxtFile, &
       MappingDestroy

contains

  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  function MappingCreate()

    implicit none

    type(mapping_type), pointer :: MappingCreate
    type(mapping_type), pointer :: map

    allocate(map)
    nullify(map%pf2clm)
    nullify(map%clm2pf)

    map%map_id                    = -1
#if 0	
	map%tmesh_num_cells_local     = -1
    map%tmesh_num_cells_ghost     = -1
    map%tmesh_num_cells_ghosted   = -1
    map%num_ocells_with_fmesh     = -1

    nullify (map%tmesh_cell_ids_ghosted_nindex)
    nullify (map%ocell_cnt_ghosted)
    nullify (map%ocell_cnt_cumsum_ghosted)
    nullify (map%ocell_ids_nindex)
    nullify (map%ocell_vol_nindex)
    nullify (map%docell_ids_sorted_nindex)
    nullify (map%hash)

    map%is_fmeshlocal_to_global      = 0
    map%is_global_to_tmeshlocal      = 0
    map%scatter_fmeshlocal_to_global = 0
    map%scatter_global_to_tmeshlocal = 0
#endif
    MappingCreate => map

  end function MappingCreate

  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingAllocateMemory( map )

    implicit none

    type(mapping_type) :: map
    PetscErrorCode     :: ierr

#if 0
    allocate( map%tmesh_cell_ids_ghosted_nindex  ( map%tmesh_num_cells_ghosted) )
    allocate( map%ocell_cnt_ghosted              ( map%tmesh_num_cells_ghosted) )
    allocate( map%ocell_cnt_cumsum_ghosted       ( map%tmesh_num_cells_ghosted) )

    allocate( map%ocell_ids_nindex ( map%num_ocells_with_fmesh) )
    allocate( map%ocell_vol_nindex ( map%num_ocells_with_fmesh) )
#endif
  end subroutine MappingAllocateMemory


  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingSetOcells ( map, ocell_ids, ocell_vol )

    implicit none

    type(mapping_type) :: map
    PetscInt           :: ocell_ids(:)
    PetscReal          :: ocell_vol(:)
    PetscInt           :: ii, jj, count
    PetscErrorCode     :: ierr
    PetscInt, pointer  :: index(:)
    PetscInt           :: rank

#if 0
    call MPI_Comm_rank(MPI_COMM_WORLD,rank, ierr)

    ! Check if number of overlapped cells > 0
	if (map%num_ocells_with_fmesh.eq.0) then
		map%num_docells_with_fmesh = 0
		return
	endif
	
    allocate (index(map%num_ocells_with_fmesh))

    do ii = 1,map%num_ocells_with_fmesh
       map%ocell_ids_nindex(ii) = ocell_ids(ii)
       map%ocell_vol_nindex(ii) = ocell_vol(ii)
       index(ii) = ii
    enddo

    ! 'Distinct' overlapped cells computation:
    ! First, sort the ocell_ids
    index = index - 1
    call PetscSortIntWithPermutation( map%num_ocells_with_fmesh, ocell_ids, index, ierr)
    index = index + 1

    ! Second, compute the number of distinct overlapped cells
    map%num_docells_with_fmesh = 1
    do ii = 2,map%num_ocells_with_fmesh
       if ( ocell_ids(index(ii)).gt.ocell_ids(index(ii-1)) ) then
          map%num_docells_with_fmesh = map%num_docells_with_fmesh + 1
       endif
    enddo

    ! Allocate memory
    allocate( map%docell_ids_sorted_nindex ( map%num_docells_with_fmesh ) )
    allocate( map%hash                     ( map%num_ocells_with_fmesh  ) )

    ! Save the ids of distinct overlapped cells
    count = 1
    map%docell_ids_sorted_nindex( count ) = ocell_ids(index(1))

    do ii = 2,map%num_ocells_with_fmesh
       if ( ocell_ids(index(ii)).gt.ocell_ids(index(ii-1)) ) then
          count = count + 1
          map%docell_ids_sorted_nindex( count ) = &
               ocell_ids(index(ii))
       endif
    enddo

    ! Populate the lookup table
    do ii = 1,map%num_ocells_with_fmesh
       do jj = 1,map%num_docells_with_fmesh
          if( ocell_ids(ii).eq.map%docell_ids_sorted_nindex(jj)) then
             map%hash(ii) = jj
          endif
       enddo
    enddo

    deallocate (index)
#endif
  end subroutine MappingSetOcells

  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingSetNumCells ( map, t_num_local, t_num_ghost, &
       t_num_ghosted, f_num_ocells )
	   
    implicit none
	
	type(mapping_type), pointer :: map
	PetscInt                    :: t_num_local, t_num_ghost
	PetscInt                    :: t_num_ghosted, f_num_ocells
#if 0	
    map%tmesh_num_cells_local   = t_num_local
    map%tmesh_num_cells_ghost   = t_num_ghost
    map%tmesh_num_cells_ghosted = t_num_ghosted
    map%num_ocells_with_fmesh   = f_num_ocells
#endif

  end subroutine MappingSetNumCells
  
  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingSetTmeshCellIds ( map, cell_ids )
  
    implicit none 
	
	type(mapping_type), pointer  :: map
	PetscInt, intent(in),pointer :: cell_ids(:)
	PetscInt                     :: local_id

#if 0	
	do local_id = 1,map%tmesh_num_cells_ghosted
	  map%tmesh_cell_ids_ghosted_nindex(local_id) = cell_ids(local_id)
	end do
#endif

  end subroutine MappingSetTmeshCellIds

  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingCreateIS(map, mycomm, fmeshlocal_npts, fmesh_cell_ids)

    implicit none

    type(mapping_type), pointer :: map
    PetscInt                    :: fmeshlocal_npts
    PetscInt, pointer           :: fmesh_cell_ids(:)
    PetscMPIInt                 :: mycomm
    PetscErrorCode               :: ierr

#if 0
    call ISCreateBlock( mycomm, 1, fmeshlocal_npts, &
         fmesh_cell_ids, PETSC_COPY_VALUES, &
         map%is_fmeshlocal_to_global, ierr)


    call ISCreateBlock( mycomm, 1, map%num_docells_with_fmesh, &
         map%docell_ids_sorted_nindex, PETSC_COPY_VALUES, &
         map%is_global_to_tmeshlocal, ierr)
#endif
  end subroutine MappingCreateIS

  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingCreateVecScatter(map, mycomm, fmesh_npts, vec_fmeshlocal, &
       vec_global, vec_tmeshlocal)

    implicit none

    type(mapping_type), pointer :: map
    PetscInt                    :: fmesh_npts
    PetscInt, allocatable       :: tmp_int_array(:)
    PetscInt                    :: ii
    PetscMPIInt                 :: mycomm
    PetscErrorCode              :: ierr

    PetscViewer :: viewer
    Vec :: vec_fmeshlocal, vec_global, vec_tmeshlocal

    IS :: is_fmeshlocal_to_fmeshlocal
    IS :: is_tmeshlocal_to_tmeshlocal
    PetscInt           :: rank

#if 0
    call MPI_Comm_rank(MPI_COMM_WORLD,rank, ierr)

    ! Create a temporary IS set 'tmeshlocal_to_tmeshlocal'
    allocate(tmp_int_array(map%num_docells_with_fmesh))
    do ii = 1,map%num_docells_with_fmesh
       tmp_int_array(ii) = ii - 1
    enddo
    call ISCreateBlock( mycomm, 1, map%num_docells_with_fmesh, tmp_int_array, PETSC_COPY_VALUES, &
         is_tmeshlocal_to_tmeshlocal, ierr)
    deallocate(tmp_int_array)
    call PetscViewerASCIIOpen(mycomm, 'is_tmeshlocal_to_tmeshlocal.out',viewer,ierr)
    call ISView(is_tmeshlocal_to_tmeshlocal,viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)

    ! Create a temporary IS set 'fmeshlocal_to_fmeshlocal'
    allocate(tmp_int_array(fmesh_npts))
    do ii = 1,fmesh_npts
       tmp_int_array(ii) = ii - 1
    enddo
    call ISCreateBlock( mycomm, 1, fmesh_npts, tmp_int_array, PETSC_COPY_VALUES, &
         is_fmeshlocal_to_fmeshlocal, ierr)
    deallocate(tmp_int_array)
    call PetscViewerASCIIOpen(mycomm, 'is_fmeshlocal_to_fmeshlocal.out',viewer,ierr)
    call ISView(is_fmeshlocal_to_fmeshlocal,viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)

    call VecScatterCreate(vec_fmeshlocal, is_fmeshlocal_to_fmeshlocal, &
         vec_global, map%is_fmeshlocal_to_global, &
         map%scatter_fmeshlocal_to_global, ierr)

    call PetscViewerASCIIOpen(mycomm, 'scatter_local_to_global.out', viewer, ierr)
    call VecScatterView(map%scatter_fmeshlocal_to_global, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)

    call VecScatterCreate(vec_global, map%is_global_to_tmeshlocal, &
         vec_tmeshlocal, is_tmeshlocal_to_tmeshlocal, &
         map%scatter_global_to_tmeshlocal, ierr)

    call PetscViewerASCIIOpen(mycomm, 'scatter_global_to_local.out', viewer, ierr)
    call VecScatterView(map%scatter_global_to_tmeshlocal, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)

    !call ISDestroy(is_fmeshlocal_to_fmeshlocal)
    !call ISDestroy(is_tmeshlocal_to_tmeshlocal)
#endif
  end subroutine MappingCreateVecScatter

  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingScatterLocal2Global(map, vec_local, vec_global)

    implicit none

    type (mapping_type), pointer :: map
    PetscErrorCode               :: ierr

    Vec :: vec_local, vec_global
#if 0
    call VecScatterBegin( map%scatter_fmeshlocal_to_global, &
         vec_local,vec_global, INSERT_VALUES,SCATTER_FORWARD, ierr)
    call VecScatterEnd( map%scatter_fmeshlocal_to_global, &
         vec_local,vec_global, INSERT_VALUES,SCATTER_FORWARD, ierr)
#endif
  end subroutine MappingScatterLocal2Global


  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingScatterGlobal2Local(map, vec_local, vec_global)

    implicit none

    type (mapping_type), pointer :: map
    PetscErrorCode               :: ierr

    Vec :: vec_local, vec_global
#if 0
    call VecScatterBegin( map%scatter_global_to_tmeshlocal, &
         vec_global, vec_local, INSERT_VALUES,SCATTER_FORWARD, ierr)
    call VecScatterEnd( map%scatter_global_to_tmeshlocal, &
         vec_global,vec_local, INSERT_VALUES,SCATTER_FORWARD, ierr)
#endif
  end subroutine MappingScatterGlobal2Local

  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingFindDistinctSourceMeshCellIds(map, option)

	use Option_module
	
    implicit none
	
	type (mapping_type), pointer :: map
    type(option_type), pointer         :: option
	PetscErrorCode               :: ierr
	PetscInt                     :: ii,jj,kk,count
	PetscInt, pointer            :: index(:), rev_index(:), index_asc(:)
	PetscInt, pointer            :: org_or_dup(:), org_or_dup_cumsum(:)
	PetscScalar, pointer         :: v_loc(:)

	Vec :: xx, yy
	
	allocate(index(            map%s2d_s_ncells))
	allocate(index_asc(        map%s2d_s_ncells))
	allocate(rev_index(        map%s2d_s_ncells))
	allocate(org_or_dup(       map%s2d_s_ncells))
	allocate(org_or_dup_cumsum(map%s2d_s_ncells))
	allocate(map%s2d_wts_jcsr( map%s2d_s_ncells))
		

	do ii = 1,map%s2d_s_ncells
	  index(ii)     = ii
	  rev_index(ii) = ii
	enddo
	
    ! Sort the s2d_s_cell_ids
    index = index - 1
    call PetscSortIntWithPermutation( map%s2d_s_ncells, map%s2d_s_cell_ids, index, ierr)
    index = index + 1

    ! Find the number of distinct s2d_s_cell_ids
    map%s2d_s_ncells_distinct = 1
	ii = 1
	kk = 1
	index_asc(kk) = index(ii)
	
    do ii = 2,map%s2d_s_ncells
       if ( map%s2d_s_cell_ids(index(ii)).gt.map%s2d_s_cell_ids(index(ii-1)) ) then
          map%s2d_s_ncells_distinct = map%s2d_s_ncells_distinct + 1
		  
		  ! Sort the 'index_asc' to ensure for the duplicate entries in the
		  ! s2d_s_cell_ids, the index values are in asceding order
		  call PetscSortInt( kk, index_asc, ierr)
		  
		  ! Update the index values
		  count = 1
		  do jj = ii-kk,ii-1
		    index(jj) = index_asc(count)
			count = count + 1
		  enddo
		  
		  ! Reset the index_asc
		  index_asc = 0
		  kk    = 1
		  index_asc(kk) = index(ii)
	   else
	      ! A duplicate entry is found. Save it in index_asc
	      kk = kk + 1
		  index_asc(kk) = index(ii)
       endif
    enddo
	
	! Sort the index_asc
	call PetscSortInt( kk, index_asc, ierr)
	count = 1 		  
	do jj = ii-kk,ii-1		  
	  index(jj) = index_asc(count)
	  count = count + 1
	enddo
	
	! Generate reverse index mapping of sorting of s2d_s_ncells
    rev_index = rev_index - 1
    call PetscSortIntWithPermutation( map%s2d_s_ncells, index, rev_index, ierr)
    rev_index = rev_index + 1

	! Save the distinct ids
    allocate(map%s2d_s_cell_ids_distinct(map%s2d_s_ncells_distinct))
	
	kk = 1
	ii = 1
	map%s2d_s_cell_ids_distinct(kk) = map%s2d_s_cell_ids(index(ii))
	org_or_dup(ii) = 1
	org_or_dup_cumsum(ii) = org_or_dup(ii)
    do ii = 2,map%s2d_s_ncells
       if ( map%s2d_s_cell_ids(index(ii)).gt.map%s2d_s_cell_ids(index(ii-1)) ) then
	     kk = kk + 1
         map%s2d_s_cell_ids_distinct(kk) = map%s2d_s_cell_ids(index(ii))
		 org_or_dup(ii) = 1
	   else
	     org_or_dup(ii) = 0 
       endif
	   org_or_dup_cumsum(ii) = org_or_dup(ii) + org_or_dup_cumsum(ii-1)
    enddo
	
	do ii = 1,map%s2d_s_ncells
	  map%s2d_wts_jcsr(ii) = rev_index(ii) - (rev_index(ii) - org_or_dup_cumsum(rev_index(ii))) - 1
	enddo
		
	deallocate(index)
	allocate(index(map%s2d_s_ncells))
	
	kk = 1
	do ii = 1,map%d_ncells_ghosted
	  do jj = 1,map%s2d_wts_icsr(ii)
	    index(kk) = ii -1
		kk = kk + 1
	  enddo
	enddo

    !
    ! size(mat_wts) = [d_ncells_ghosted x s2d_s_ncells_distinct]
	!
	call MatCreateSeqAIJ(PETSC_COMM_SELF, map%d_ncells_ghosted, &
	     map%s2d_s_ncells_distinct, PETSC_NULL, &             
		 map%s2d_wts_icsr, map%mat_wts, ierr)
	
    do ii = 1,map%s2d_s_ncells
	  call MatSetValues(map%mat_wts,1,index(ii),1,map%s2d_wts_jcsr(ii),map%s2d_wts(ii),INSERT_VALUES,ierr)
	enddo
	call MatAssemblyBegin(map%mat_wts,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(map%mat_wts,MAT_FINAL_ASSEMBLY,ierr)
	
	if(option%myrank.eq.1) then
	  !call MatView(map%mat_wts,PETSC_VIEWER_STDOUT_SELF )
	  call VecCreateSeq(PETSC_COMM_SELF, map%s2d_s_ncells_distinct, xx,ierr)  
	  call VecCreateSeq(PETSC_COMM_SELF, map%d_ncells_ghosted     , yy,ierr)  
	  ii = 0
	  
	  deallocate(index)
	  allocate(v_loc(map%s2d_s_ncells_distinct))
	  allocate(index(map%s2d_s_ncells_distinct))
	  
	  do ii=1,map%s2d_s_ncells_distinct
	    v_loc(ii) = 1.0d0
		index(ii) = ii-1
	  enddo
	    call VecSetValues(xx,map%s2d_s_ncells_distinct,index,v_loc,INSERT_VALUES,ierr)
		!call VecView(xx,PETSC_VIEWER_STDOUT_SELF)
	  call MatMult( map%mat_wts, xx, yy, ierr)
	  !call VecView(yy,PETSC_VIEWER_STDOUT_SELF)
	endif
	

    deallocate(index)
    deallocate(index_asc)
	deallocate(rev_index)
	deallocate(org_or_dup)
	deallocate(org_or_dup_cumsum)
	
  end subroutine MappingFindDistinctSourceMeshCellIds
  
  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingDestroy (map)

    implicit none

    type (mapping_type), pointer :: map
    PetscErrorCode               :: ierr

#if 0
    if(.not.associated(map)) return

    if(associated(map%tmesh_cell_ids_ghosted_nindex )) deallocate(map%tmesh_cell_ids_ghosted_nindex)
    if(associated(map%ocell_cnt_ghosted       )) deallocate(map%ocell_cnt_ghosted)
    if(associated(map%ocell_cnt_cumsum_ghosted)) deallocate(map%ocell_cnt_cumsum_ghosted)
    if(associated(map%ocell_ids_nindex        )) deallocate(map%ocell_ids_nindex)
    if(associated(map%ocell_vol_nindex        )) deallocate(map%ocell_vol_nindex)
    if(associated(map%docell_ids_sorted_nindex)) deallocate(map%docell_ids_sorted_nindex)
    if(associated(map%hash                    )) deallocate(map%hash)

    call ISDestroy(map%is_fmeshlocal_to_global,ierr)
    call ISDestroy(map%is_global_to_tmeshlocal,ierr)
    call VecScatterDestroy(map%scatter_fmeshlocal_to_global,ierr)
    call VecScatterDestroy(map%scatter_global_to_tmeshlocal,ierr)
#endif

  end subroutine MappingDestroy

  ! ************************************************************************** !
  !
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
	PetscInt,pointer                   :: s2d_wts_icsr(:)
	PetscInt,pointer                   :: s2d_s_cell_ids(:)
	PetscInt,pointer                   :: wts_mat_row(:), wts_mat_col(:)
	PetscReal,pointer                  :: wts_mat(:), wt_ovp_tmp(:)
    PetscErrorCode                     :: ierr

    PetscMPIInt                        :: status_mpi(MPI_STATUS_SIZE)
    
    fileid   = 20
    card     = 'pflotran_interface_main'


    if(option%myrank == option%io_rank) then
	
      input => InputCreate(fileid,filename)
      call InputReadFlotranString(input,option)
      call InputErrorMsg(input,option,'number of cells',card)

      num_wts = -1
      call InputReadInt(input,option,num_wts)

       allocate(wts_mat_row(num_wts))
	   allocate(wts_mat_col(num_wts))
	   allocate(wts_mat(    num_wts))
	   
	   do ii = 1,num_wts
	     call InputReadFlotranString(input,option)
	     call InputReadInt(input,option,wts_mat_row(ii))
		 call InputReadInt(input,option,wts_mat_col(ii))
		 call InputReadDouble(input,option,wts_mat(ii))
		 
		 ! Row/Col in the txt file are in 1-based indexing,
		 ! converting it into 0-based indexing
		 wts_mat_row(ii) = wts_mat_row(ii) - 1d0
		 wts_mat_col(ii) = wts_mat_col(ii) - 1d0
		 
	   enddo
	   
	   do irank = 0,option%mycommsize - 1

	     if(irank.ne.option%io_rank) then	   

   	       ! 0) Get from irank-th processor information regarding destination mesh:
		   !    - Number of cells present
		   !    - Cell ids
	       call MPI_Recv(d_ncells_tmp, 1, MPI_INTEGER, irank, MPI_ANY_TAG, &
		        option%mycomm, status_mpi, ierr)
		 
		   allocate(d_cells_ids_tmp(               d_ncells_tmp))
		   allocate(s2d_wts_icsr(d_ncells_tmp))
		   call MPI_Recv(d_cells_ids_tmp, d_ncells_tmp, MPI_INTEGER, irank, MPI_ANY_TAG, &
		        option%mycomm, status_mpi, ierr)
		 else
		   ! 0) Get local data
		   d_ncells_tmp = map%d_ncells_ghosted
		   allocate(d_cells_ids_tmp(               d_ncells_tmp))
		   allocate(s2d_wts_icsr(d_ncells_tmp))
		   
		   do ii=1,map%d_ncells_ghosted
		     d_cells_ids_tmp(ii) = map%d_cell_ids_sort(ii)
		   enddo
		 endif
		 
		 ! 1) Find the number of overlapped cells
		 s2d_s_ncells = 0d0
		 do ii = 1,d_ncells_tmp
		   do jj = 1,num_wts
		     if(d_cells_ids_tmp(ii).eq.wts_mat_row(jj)) then
			   s2d_s_ncells = s2d_s_ncells + 1d0
			 endif
		   enddo 
		 enddo
		 
		 ! 2) Save the cell-ids of overlapped cells
		 s2d_wts_icsr = 0
		 
		 if (s2d_s_ncells.gt.0) then
  		   allocate(s2d_s_cell_ids(s2d_s_ncells))
		   allocate(wt_ovp_tmp(    s2d_s_ncells))
		 
		   kk = 0d0
		   do ii = 1,d_ncells_tmp
		     do jj = 1,num_wts
		       if(d_cells_ids_tmp(ii).eq.wts_mat_row(jj)) then
			     kk                 = kk + 1d0
			     s2d_s_cell_ids(kk) = wts_mat_col(jj)
		  	     wt_ovp_tmp(kk)     = wts_mat(jj)
				 s2d_wts_icsr(ii) = s2d_wts_icsr(ii) + 1d0
			   endif
		     enddo
		   enddo
         endif
		
		 if(irank.ne.option%io_rank) then
  		   ! 3) Send information back to irank-th processor
 		   call MPI_Send(s2d_s_ncells  , 1         , MPI_INTEGER, irank, option%myrank, &
		        option%mycomm, ierr)
  	       call MPI_Send(s2d_wts_icsr,d_ncells_tmp, MPI_INTEGER, irank, option%myrank, &
		        option%mycomm, ierr)

		   if(s2d_s_ncells.gt.0) then   
		      call MPI_Send(wt_ovp_tmp    , s2d_s_ncells, MPI_DOUBLE_PRECISION   , irank, option%myrank, option%mycomm, ierr)
     		  call MPI_Send(s2d_s_cell_ids, s2d_s_ncells, MPI_INTEGER, irank, option%myrank, option%mycomm, ierr)

              ! Free memory
			  deallocate(s2d_s_cell_ids)
			  deallocate(wt_ovp_tmp)
		    endif	
		  else
		    ! 3) Save local information
		    allocate(map%s2d_wts_icsr(d_ncells_tmp))
			do ii = 1,d_ncells_tmp
			  map%s2d_wts_icsr(ii) = s2d_wts_icsr(ii)
			enddo
			
			
		    map%s2d_s_ncells = s2d_s_ncells
			if (map%s2d_s_ncells.gt.0) then
			  allocate(map%s2d_s_cell_ids(map%s2d_s_ncells))
			  allocate(map%s2d_wts(       map%s2d_s_ncells))
			  
			  do ii = 1,map%s2d_s_ncells
			    map%s2d_s_cell_ids(ii)  = s2d_s_cell_ids(ii)
				map%s2d_wts(ii)         = wt_ovp_tmp(ii)
			  enddo
			  			  
			  ! Free memory
			  deallocate(s2d_s_cell_ids)
			  deallocate(wt_ovp_tmp)
			endif
			
		  endif
		  		  
		  ! Free memory
		  deallocate(d_cells_ids_tmp)
		  deallocate(s2d_wts_icsr)

	   enddo
	   
	   ! Free memory
       deallocate(wts_mat_row)
	   deallocate(wts_mat_col)
	   deallocate(wts_mat)

	else
	
	   call MPI_Send(map%d_ncells_ghosted, 1, MPI_INTEGER, option%io_rank, &
	        option%myrank, option%mycomm, ierr)
	   call MPI_Send(map%d_cell_ids_sort, map%d_ncells_ghosted,&
	        MPI_INTEGER, option%io_rank, option%myrank, option%mycomm, ierr)

	   call MPI_Recv(map%s2d_s_ncells, 1, MPI_INTEGER, option%io_rank, MPI_ANY_TAG, &
	        option%mycomm,status_mpi, ierr)

       allocate(map%s2d_wts_icsr(map%d_ncells_ghosted))
       call MPI_Recv(map%s2d_wts_icsr, map%d_ncells_ghosted, MPI_INTEGER,&
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

    write(*,*), 'map%s_ncells_local = ', map%s_ncells_local
  end subroutine MappingSetSourceMeshCellIds


  ! ************************************************************************** !
  !
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
	IS                           :: is_ltol, is_ltog
	Vec                          :: vec_tmp_1, vec_tmp_2
    PetscViewer :: viewer
	PetscInt                     :: istart, iend,count
	
	map%d_ncells_local   = num_cells_local
	map%d_ncells_ghost   = num_cells_ghost
	map%d_ncells_ghosted = num_cells_local + num_cells_ghost
	
    call VecCreateMPI(option%mycomm, map%d_ncells_local, PETSC_DECIDE, vec_tmp_1, ierr)
    call VecCreateMPI(option%mycomm, map%d_ncells_local, PETSC_DECIDE, vec_tmp_2, ierr)

	allocate(map%d_cell_ids(        map%d_ncells_ghosted))
	allocate(map%d_cell_ids_sort(   map%d_ncells_ghosted))
	allocate(map%d_local_or_ghost(  map%d_ncells_ghosted))
	allocate(map%d_nG2S(            map%d_ncells_ghosted))
	allocate(map%d_nS2G(            map%d_ncells_ghosted))
	allocate(index(                 map%d_ncells_ghosted))
	allocate(rev_index(             map%d_ncells_ghosted))
			
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

	do ii = 1,num_cells_local+num_cells_ghost
	  map%d_cell_ids_sort(ii) = cell_ids_ghosted(index(ii))
	  map%d_nG2S(ii)             = index(ii)
	enddo

    ! Sort the d_cell_ids
    rev_index = rev_index - 1
    call PetscSortIntWithPermutation( map%d_ncells_ghosted, index, rev_index, ierr)
    rev_index = rev_index + 1
	
	do ii = 1,num_cells_local+num_cells_ghost
	  map%d_nS2G(ii)             = rev_index(ii)
	enddo
	
    ! is_d_local_to_natural
    call VecGetOwnershipRange(vec_tmp_1,istart,iend,ierr)
	allocate(tmp_int_array(map%d_ncells_local))
	do ii=1,map%d_ncells_local
	  tmp_int_array(ii) = istart+ii-1
	  !if(option%myrank.eq.1) write(*,*),ii,tmp_int_array(ii)
	enddo
	write(*,*)
	call ISCreateBlock(option%mycomm,1,map%d_ncells_local, tmp_int_array, PETSC_COPY_VALUES, &
	     is_ltol, ierr)
    
	!call ISView(is_ltol,PETSC_VIEWER_STDOUT_WORLD,ierr)
	
	count = count + 1
	do ii=1,map%d_ncells_local
	  if(map%d_local_or_ghost(ii).eq.1) then
	    count = count + 1
		tmp_int_array(count) = map%d_cell_ids(ii)
	  endif 
	  !if(option%myrank.eq.1) write(*,*),ii,tmp_int_array(ii)
	enddo
	call ISCreateBlock(option%mycomm,1,map%d_ncells_local, tmp_int_array, PETSC_COPY_VALUES, &
	     is_ltog, ierr)
    deallocate(tmp_int_array)
	
	!call ISView(is_ltog,PETSC_VIEWER_STDOUT_WORLD,ierr)
	
	call VecScatterCreate(vec_tmp_1, is_ltol, vec_tmp_2, is_ltog, map%d_local_to_natural_vscat, ierr)

    call PetscViewerASCIIOpen(option%mycomm, 'd_local_to_natural_vscat.out', viewer, ierr)
	call VecScatterView(map%d_local_to_natural_vscat, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
   

	deallocate(index)
    deallocate(rev_index)
	call VecDestroy(vec_tmp_1)
	call VecDestroy(vec_tmp_2)
	call ISDestroy(is_ltol)
	call ISDestroy(is_ltog)
	
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

    call VecScatterCreate(nA2P, is_gtol, porder_req, is_ltol, vscat_2, ierr)
		
	call VecScatterBegin(vscat_2, nA2P, porder_req, INSERT_VALUES, SCATTER_FORWARD, ierr)
	call VecScatterEnd(  vscat_2, nA2P, porder_req, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call PetscViewerASCIIOpen(option%mycomm, 's2d_vscat.out', viewer, ierr)
	call VecScatterView(vscat_2, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
	
	call ISDestroy(is_gtol,ierr)
	call VecGetArrayF90(porder_req,v_loc,ierr)
	call ISCreateBlock(option%mycomm, 1, map%s2d_s_ncells_distinct, INT(v_loc),&
		  PETSC_COPY_VALUES, is_gtol, ierr)
	call VecRestoreArrayF90(porder_req,v_loc,ierr)	
	!call ISView(is_gtol, PETSC_VIEWER_STDOUT_WORLD,ierr)
	
	
	call VecScatterCreate(nA2P, is_gtol, porder_req, is_ltol, map%s2d_vscat, ierr)
	
	
	call VecDestroy(porder, ierr)
	call VecDestroy(aorder, ierr)
	call VecDestroy(nA2P  , ierr)
	call VecScatterDestroy(vscat_1,ierr)
	call VecScatterDestroy(vscat_2,ierr)
    call ISDestroy(is_gtol,ierr)
	call ISDestroy(is_ltol, ierr)
	
	
  end subroutine MappingCreateScatterOfSourceMesh


  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingCreateWeightMatrix(map, option)

    use Option_module
    implicit none
  
	type(mapping_type), pointer  :: map
    type(option_type), pointer   :: option

  end subroutine MappingCreateWeightMatrix

  ! ************************************************************************** !
  !
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
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_s_ncells_distinct, vec_s_loc, ierr)
	call VecCreateSeq(PETSC_COMM_SELF, map%d_ncells_ghosted     , vec_d_loc, ierr)
	call VecCreateMPI(option%mycomm  , map%d_ncells_local       , PETSC_DECIDE, vec_d_nat, ierr)
	
	
    !call VecView(vec_s      , PETSC_VIEWER_STDOUT_WORLD, ierr)
	call VecScatterBegin(map%s2d_vscat, vec_s, vec_s_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
	call VecScatterEnd(  map%s2d_vscat, vec_s, vec_s_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
   	!call VecView(vec_s_local, PETSC_VIEWER_STDOUT_WORLD, ierr)
   
	call VecGetArrayF90(vec_s_local,v_loc_1,ierr)
	call VecGetArrayF90(vec_s_loc  ,v_loc_2,ierr)
    do ii = 1, map%s2d_s_ncells_distinct
	  v_loc_2(ii) = v_loc_1(ii)
	enddo
	call VecRestoreArrayF90(vec_s_local,v_loc_1,ierr)
	call VecRestoreArrayF90(vec_s_loc  ,v_loc_2,ierr)

	
	  
	!call VecView(vec_s_loc  , PETSC_VIEWER_STDOUT_SELF, ierr)
	!call MatView(map%mat_wts, PETSC_VIEWER_STDOUT_SELF, ierr)
    call MatMult(map%mat_wts, vec_s_loc, vec_d_loc, ierr)
	
	if(option%myrank.eq.1) then
	  !call MatView(map%mat_wts,PETSC_VIEWER_STDOUT_SELF, ierr)
	  !call VecView(vec_d_loc, PETSC_VIEWER_STDOUT_SELF, ierr)
      call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'u_b_2.txt',viewer,ierr)
      call VecView(vec_d_loc,viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)
	else
      call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'u_b_1.txt',viewer,ierr)
      call VecView(vec_d_loc,viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)
	endif
	
	call MPI_Barrier( PETSC_COMM_WORLD, ierr)  
	
	
	call VecGetArrayF90(vec_d    ,v_loc_1,ierr)
	call VecGetArrayF90(vec_d_loc,v_loc_2,ierr)
	do ii = 1,map%d_ncells_local
	  v_loc_1(ii) = v_loc_2(ii)
	enddo
	call VecRestoreArrayF90(vec_d    ,v_loc_1,ierr)
	call VecRestoreArrayF90(vec_d_loc,v_loc_2,ierr)
	
	call VecScatterBegin(map%d_local_to_natural_vscat, vec_d, vec_d_nat, INSERT_VALUES, SCATTER_FORWARD, ierr)
	call VecScatterEnd(  map%d_local_to_natural_vscat, vec_d, vec_d_nat, INSERT_VALUES, SCATTER_FORWARD, ierr)
	
	!call VecView(vec_d_nat, PETSC_VIEWER_STDOUT_WORLD, ierr)
	
   end subroutine MappingSourceToDestination


end module Mapping_module
