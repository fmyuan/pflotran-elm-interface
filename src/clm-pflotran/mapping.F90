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

     ! This module does mapping between two meshs:
     !  'from_mesh'      - The mesh from which data is to be mapped.
     !  'to_mesh'        - The mesh to which data is being mapped.
     !  'ocell'          - Cells of 'from_mesh' which overlapped with 'to_mesh'
     !  'docell'         - Since a cell 'from_mesh' can overlap with multiple cells
     !                     of 'to_mesh', this refers to distinct ocells (i.e. without duplications)

     PetscInt :: map_id                            !

     ! Size of 'to_mesh'
     PetscInt :: num_cells_local                    ! Number of local cells
     PetscInt :: num_cells_ghost                    ! Number of ghost cells
     PetscInt :: num_cells_ghosted                  ! num_cells_local + num_cells_ghost

     PetscInt :: num_ocells                         ! Total number of cells 'from_mesh'
     !                                                which overlap with 'to_mesh'
     PetscInt :: num_docells                        ! Number of distinct cells 'from_mesh'
     !                                                which overlap with 'to_mesh'

     ! 'to_mesh' Cells:
     PetscInt, pointer  :: cell_ids_ghosted_nindex(:)      ! Cell-IDs in Natural indexing             [num_cells_ghosted]

     !                                                    'from_mesh' Overlapped Cells:
     PetscInt, pointer  :: ocell_cnt_ghosted(:)            ! For each 'to_mesh' cell, the count of
     !                                                       cells within 'from_mesh' overlapped with [num_cells_ghosted]
     PetscInt, pointer  :: ocell_cnt_cumsum_ghosted(:)     ! Cummulative sum of 'ocell_cnt_ghosted'   [num_cells_ghosted]
     PetscInt, pointer  :: ocell_ids_nindex (:)            ! Overlapped Cell-IDs in Natural indexing  [num_ocells]
     PetscReal, pointer :: ocell_vol_nindex(:)             ! Overlapped volume in Natural indexing    [num_ocells]

     PetscInt, pointer  :: docell_ids_sorted_nindex(:)     ! Distinct overlapped Cell-IDs sorted and
     !                                                       in natural indexing                      [num_docells]
     PetscInt,pointer   :: hash(:)                         ! Lookup table for obtaining position in 
     !                                                       the sorted distinct cell ids 
     !                                                       (docell_ids_sorted) list for every 
	 !                                                       overlapped cell-id (ocell)               [num_ocells]

  end type mapping_type

  public :: MappingCreate, &
       MappingAllocateMemory, &
       MappingSetOcells, &
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

    map%map_id              = -1
    map%num_cells_local     = -1
    map%num_cells_ghost     = -1
    map%num_cells_ghosted   = -1
    map%num_ocells          = -1

    nullify (map%cell_ids_ghosted_nindex)
    nullify (map%ocell_cnt_ghosted)
    nullify (map%ocell_cnt_cumsum_ghosted)
    nullify (map%ocell_ids_nindex)
    nullify (map%ocell_vol_nindex)
    nullify (map%docell_ids_sorted_nindex)
    nullify (map%hash)

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

    allocate( map%cell_ids_ghosted_nindex  ( map%num_cells_ghosted) )
    allocate( map%ocell_cnt_ghosted        ( map%num_cells_ghosted) )
    allocate( map%ocell_cnt_cumsum_ghosted ( map%num_cells_ghosted) )

    allocate( map%ocell_ids_nindex ( map%num_ocells) )
    allocate( map%ocell_vol_nindex ( map%num_ocells) )

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

    allocate (index(map%num_ocells))

    do ii = 1,map%num_ocells
       map%ocell_ids_nindex(ii) = ocell_ids(ii)
       map%ocell_vol_nindex(ii) = ocell_vol(ii)
       index(ii) = ii
    enddo

    ! 'Distinct' overlapped cells computation:
    ! First, sort the ocell_ids
    index = index - 1
    call PetscSortIntWithPermutation( map%num_ocells, ocell_ids, index, ierr)
    index = index + 1

    ! Second, compute the number of distinct overlapped cells
    map%num_docells = 1
    do ii = 2,map%num_ocells
       if ( ocell_ids(index(ii)).gt.ocell_ids(index(ii-1)) ) then
          map%num_docells = map%num_docells + 1
       endif
    enddo

    !print *, 'num_docells: ', map%num_docells, 'num_ocells',map%num_ocells
    ! Allocate memory
    allocate( map%docell_ids_sorted_nindex ( map%num_docells ) )
    allocate( map%hash                     ( map%num_ocells  ) )

    ! Save the ids of distinct overlapped cells
    count = 1
    map%docell_ids_sorted_nindex( count ) = ocell_ids(index(1))
    !print *, '1', ocell_ids(index(1))

    do ii = 2,map%num_ocells
       if ( ocell_ids(index(ii)).gt.ocell_ids(index(ii-1)) ) then
          count = count + 1
          map%docell_ids_sorted_nindex( count ) = &
               ocell_ids(index(ii))
          !print *, count, ocell_ids(index(ii))
       endif
    enddo

    ! Populate the lookup table
    print *, 'Hash table: '
    do ii = 1,map%num_ocells
       do jj = 1,map%num_docells
          if( ocell_ids(ii).eq.map%docell_ids_sorted_nindex(jj)) then
             map%hash(ii) = jj
             !print *, ocell_ids(ii), map%docell_ids_sorted_nindex(jj), &
             !     jj
          endif
       enddo
    enddo

    deallocate (index)

  end subroutine MappingSetOcells


  ! ************************************************************************** !
  !
  !
  ! ************************************************************************** !
  subroutine MappingDestroy (map)

    implicit none 

    type (mapping_type), pointer :: map
    PetscErrorCode               :: ierr

    if(.not.associated(map)) return

    if(associated(map%cell_ids_ghosted_nindex )) deallocate(map%cell_ids_ghosted_nindex)
    if(associated(map%ocell_cnt_ghosted       )) deallocate(map%ocell_cnt_ghosted)
    if(associated(map%ocell_cnt_cumsum_ghosted)) deallocate(map%ocell_cnt_cumsum_ghosted)
    if(associated(map%ocell_ids_nindex        )) deallocate(map%ocell_ids_nindex)
    if(associated(map%ocell_vol_nindex        )) deallocate(map%ocell_vol_nindex)
    if(associated(map%docell_ids_sorted_nindex)) deallocate(map%docell_ids_sorted_nindex)
    if(associated(map%hash                    )) deallocate(map%hash)

  end subroutine MappingDestroy

end module Mapping_module
