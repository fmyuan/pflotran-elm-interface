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

     ! Size of 'tmesh'
     PetscInt :: tmesh_num_cells_local                    ! Number of local cells
     PetscInt :: tmesh_num_cells_ghost                    ! Number of ghost cells
     PetscInt :: tmesh_num_cells_ghosted                  ! num_cells_local + num_cells_ghost

     PetscInt :: num_ocells_with_fmesh                    ! Total number of cells 'fmesh'
     !                                                      which overlap with 'tmesh'
     PetscInt :: num_docells_with_fmesh                   ! Number of distinct cells 'fmesh'
     !                                                      which overlap with 'tmesh'

     ! 'tmesh' Cells:
     PetscInt, pointer  :: tmesh_cell_ids_ghosted_nindex(:)! Cell-IDs in Natural indexing             [tmesh_num_cells_ghosted]

     !                                                     'fmesh' Overlapped Cells:
     PetscInt, pointer  :: ocell_cnt_ghosted(:)            ! For each 'tmesh' cell, the count of
     !                                                       cells within 'fmesh' overlapped with     [tmesh_num_cells_ghosted]
     PetscInt, pointer  :: ocell_cnt_cumsum_ghosted(:)     ! Cummulative sum of 'ocell_cnt_ghosted'   [tmesh_num_cells_ghosted]
     PetscInt, pointer  :: ocell_ids_nindex (:)            ! Overlapped Cell-IDs in Natural indexing  [num_ocells_with_fmesh]
     PetscReal, pointer :: ocell_vol_nindex(:)             ! Overlapped volume in Natural indexing    [num_ocells_with_fmesh]

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


  end type mapping_type

  public :: MappingCreate, &
       MappingAllocateMemory, &
       MappingSetOcells, &
       MappingCreateIS, &
       MappingCreateVecScatter, &
       MappingScatterLocal2Global, &
       MappingScatterGlobal2Local, &
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

    allocate( map%tmesh_cell_ids_ghosted_nindex  ( map%tmesh_num_cells_ghosted) )
    allocate( map%ocell_cnt_ghosted        ( map%tmesh_num_cells_ghosted) )
    allocate( map%ocell_cnt_cumsum_ghosted ( map%tmesh_num_cells_ghosted) )

    allocate( map%ocell_ids_nindex ( map%num_ocells_with_fmesh) )
    allocate( map%ocell_vol_nindex ( map%num_ocells_with_fmesh) )

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

  end subroutine MappingSetOcells

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

    call ISCreateBlock( mycomm, 1, fmeshlocal_npts, &
         fmesh_cell_ids, PETSC_COPY_VALUES, &
         map%is_fmeshlocal_to_global, ierr)


    call ISCreateBlock( mycomm, 1, map%num_docells_with_fmesh, &
         map%docell_ids_sorted_nindex, PETSC_COPY_VALUES, &
         map%is_global_to_tmeshlocal, ierr)

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

    call PetscViewerASCIIOpen(mycomm, 'scatter_clmlocal_to_global.out', viewer, ierr)
    call VecScatterView(map%scatter_fmeshlocal_to_global, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)

    call VecScatterCreate(vec_global, map%is_global_to_tmeshlocal, &
         vec_tmeshlocal, is_tmeshlocal_to_tmeshlocal, &
         map%scatter_global_to_tmeshlocal, ierr)

    call PetscViewerASCIIOpen(mycomm, 'scatter_global_to_pflocal.out', viewer, ierr)
    call VecScatterView(map%scatter_global_to_tmeshlocal, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)

    !call ISDestroy(is_fmeshlocal_to_fmeshlocal)
    !call ISDestroy(is_tmeshlocal_to_tmeshlocal)

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

    call VecScatterBegin( map%scatter_fmeshlocal_to_global, &
         vec_local,vec_global, INSERT_VALUES,SCATTER_FORWARD, ierr)
    call VecScatterEnd( map%scatter_fmeshlocal_to_global, &
         vec_local,vec_global, INSERT_VALUES,SCATTER_FORWARD, ierr)

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

    call VecScatterBegin( map%scatter_global_to_tmeshlocal, &
         vec_global, vec_local, INSERT_VALUES,SCATTER_FORWARD, ierr)
    call VecScatterEnd( map%scatter_global_to_tmeshlocal, &
         vec_global,vec_local, INSERT_VALUES,SCATTER_FORWARD, ierr)

  end subroutine MappingScatterGlobal2Local

  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingDestroy (map)

    implicit none

    type (mapping_type), pointer :: map
    PetscErrorCode               :: ierr

    if(.not.associated(map)) return

    if(associated(map%tmesh_cell_ids_ghosted_nindex )) deallocate(map%tmesh_cell_ids_ghosted_nindex)
    if(associated(map%ocell_cnt_ghosted       )) deallocate(map%ocell_cnt_ghosted)
    if(associated(map%ocell_cnt_cumsum_ghosted)) deallocate(map%ocell_cnt_cumsum_ghosted)
    if(associated(map%ocell_ids_nindex        )) deallocate(map%ocell_ids_nindex)
    if(associated(map%ocell_vol_nindex        )) deallocate(map%ocell_vol_nindex)
    if(associated(map%docell_ids_sorted_nindex)) deallocate(map%docell_ids_sorted_nindex)
    if(associated(map%hash                    )) deallocate(map%hash)

    call ISDestroy(map%is_fmeshlocal_to_global)
    call ISDestroy(map%is_global_to_tmeshlocal)
    call VecScatterDestroy(map%scatter_fmeshlocal_to_global)
    call VecScatterDestroy(map%scatter_global_to_tmeshlocal)

  end subroutine MappingDestroy

end module Mapping_module
