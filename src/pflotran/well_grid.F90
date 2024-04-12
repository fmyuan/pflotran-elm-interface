module Well_Grid_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Geometry_module
  use PFLOTRAN_Constants_module

implicit none

private

type, public :: well_grid_type
    ! number of well segments
    PetscInt :: nsegments
    ! number of well connections
    PetscInt :: nconnections
    ! well index of each well segment [0,1]  0 = cased; 1 = open
    PetscReal, pointer :: casing(:)
    ! delta h discretization of each segment center [m]
    PetscReal, pointer :: dh(:)
    ! x-, y-, and z- deltas of each well segment
    PetscReal, pointer :: dx(:)
    PetscReal, pointer :: dy(:)
    PetscReal, pointer :: dz(:)
    ! reservoir dz
    PetscReal, pointer :: res_dz(:)
    ! h coordinate of each segment center [m]
    type(point3d_type), pointer :: h(:)
    ! the local id of the reservoir grid cell within which each segment
    ! center resides
    PetscInt, pointer :: h_local_id(:)
    ! the ghosted id of the reservoir grid cell within which each segment
    ! center resides
    PetscInt, pointer :: h_ghosted_id(:)
    ! the global id of the reservoir grid cell within which each segment
    ! center resides
    PetscInt, pointer :: h_global_id(:)
    ! the rank id of the reservoir grid cell within which each segment
    ! center resides
    PetscInt, pointer :: h_rank_id(:)
    ! the strata id number associated with each well segment
    PetscInt, pointer :: strata_id(:)
    ! coordinate of the top/bottom of the well [m]
    PetscReal :: tophole(3)
    PetscReal :: bottomhole(3)
    ! x-y span search distance multipier
    PetscInt :: xy_span_multiplier
    ! flag to match well grid discretization to reservoir grid
    PetscBool :: match_reservoir
    ! dz for "match reservoir" search and peck method [m]
    PetscReal :: dz_peck
    ! minimum dz for "match reservoir" method [m]
    PetscReal :: min_dz
    ! ratio for # well cells per reservoir cell
    PetscInt :: well_res_ratio
    ! list of segment center z values [m]
    PetscReal, pointer :: z_list(:)
    ! list of segment length values [m]
    PetscReal, pointer :: l_list(:)
    type(deviated_well_type), pointer :: deviated_well_segment_list
end type well_grid_type

type, public :: deviated_well_type
    PetscReal :: surface_origin(3)
    PetscReal :: dxyz(3)
    PetscReal :: radius_to_horizontal_x
    PetscReal :: radius_to_horizontal_y
    PetscReal :: radius_to_horizontal_angle(2)
    PetscReal :: radius_to_vertical
    PetscBool :: cased
    type(deviated_well_type), pointer :: next
end type deviated_well_type

public :: WellGridCreate, &
          WellSegmentInit, &
          WellGridAddConnectionsExplicit, &
          WellGridDestroy

contains

! ************************************************************************** !

function WellGridCreate()
    !
    ! Creates grid variables.
    !
    ! Author: Michael Nole
    ! Date: 03/04/2023
  
    implicit none
  
    type(well_grid_type), pointer :: WellGridCreate
    type(well_grid_type), pointer :: well_grid
  
    ! create the well grid object:
    allocate(well_grid)
    well_grid%nsegments = UNINITIALIZED_INTEGER
    well_grid%nconnections = UNINITIALIZED_INTEGER
    nullify(well_grid%casing)
    nullify(well_grid%dh)
    nullify(well_grid%dx)
    nullify(well_grid%dy)
    nullify(well_grid%dz)
    nullify(well_grid%res_dz)
    nullify(well_grid%h)
    nullify(well_grid%h_local_id)
    nullify(well_grid%h_ghosted_id)
    nullify(well_grid%h_global_id)
    nullify(well_grid%h_rank_id)
    nullify(well_grid%strata_id)
    well_grid%tophole(:) = UNINITIALIZED_DOUBLE
    well_grid%bottomhole(:) = UNINITIALIZED_DOUBLE
    well_grid%xy_span_multiplier = 10
    well_grid%match_reservoir = PETSC_FALSE
    well_grid%dz_peck = 1.0d-2
    well_grid%min_dz = UNINITIALIZED_DOUBLE
    well_grid%well_res_ratio = 1
    nullify(well_grid%z_list)
    nullify(well_grid%l_list)
    nullify(well_grid%deviated_well_segment_list)
  
    WellGridCreate => well_grid
  
end function WellGridCreate

! ************************************************************************** !

subroutine WellSegmentInit(well_segment)
  !
  ! Initialize values for a well segment defined by a specific trajectory.
  !
  ! Author: Michael Nole
  ! Date: 04/11/2024
  !

  implicit none

  type(deviated_well_type), pointer :: well_segment

  well_segment%surface_origin(:) = UNINITIALIZED_DOUBLE
  well_segment%dxyz(:) = UNINITIALIZED_DOUBLE
  well_segment%radius_to_horizontal_x = UNINITIALIZED_DOUBLE
  well_segment%radius_to_horizontal_y = UNINITIALIZED_DOUBLE
  well_segment%radius_to_horizontal_angle(:) = UNINITIALIZED_DOUBLE
  well_segment%radius_to_vertical = UNINITIALIZED_DOUBLE
  well_segment%cased = PETSC_FALSE

end subroutine WellSegmentInit

! ************************************************************************** !

subroutine WellGridDestroy(well_grid)
    !
    ! Destroys a well grid.
    !
    ! Author: Michael Nole
    ! Date: 04/02/2024
    !
    
    use Utility_module, only : DeallocateArray
  
    implicit none
  
    type(well_grid_type), pointer :: well_grid
  
  
    call DeallocateArray(well_grid%h_local_id)
    call DeallocateArray(well_grid%h_ghosted_id)
    call DeallocateArray(well_grid%dh)
    nullify(well_grid)
  
end subroutine WellGridDestroy
  
! ************************************************************************** !

subroutine WellGridAddConnectionsExplicit(cell_centroids,connections,&
                                          face_areas, face_centroids,&
                                          well_grid,option)
  !
  ! Modifies grid connectivity to include well connections.
  !
  ! Author: Michael Nole
  ! Date: 04/03/2024
  !

  use Option_module

  implicit none

  type(point3d_type), pointer :: cell_centroids(:)
  PetscInt, pointer :: connections(:,:)
  PetscReal, pointer :: face_areas(:)
  type(point3d_type), pointer :: face_centroids(:)
  type(well_grid_type), pointer :: well_grid
  type(option_type), pointer :: option

  PetscInt :: n_well_conn, total_connections, num_connections, iconn
  PetscInt :: isegment, dual_segment, local_id, dual_id
  PetscBool, allocatable :: well_connections(:,:)
  PetscBool, allocatable :: mask(:)
  PetscInt, pointer :: new_connections(:,:)
  PetscReal, pointer :: new_face_areas(:)
  type(point3d_type), pointer :: new_face_centroids(:)

  num_connections = size(connections,2)
  allocate(mask(num_connections))
  n_well_conn = 0
  
  ! Add embedded well connectivity. This just connects the bottom-hole
  ! cell to all other reservoir cells connected to the well. Does not
  ! connect those other reservoir cells to each other. This is therefore
  ! only applicable to a HYDROSTATIC well model type.
  allocate(well_connections(well_grid%nsegments,well_grid%nsegments))
  well_connections(:,:) = PETSC_FALSE
  dual_segment = 1
  do isegment = 1,well_grid%nsegments
    if (well_grid%h_rank_id(isegment) /= option%myrank) cycle
    if (well_grid%casing(isegment) <= 0.d0) cycle
      local_id = well_grid%h_ghosted_id(isegment)
    ! For each local cell, check if a connection to a well cell exists. If 
    ! not, add it, and keep track of connections that are only well-related.
    dual_id = well_grid%h_ghosted_id(dual_segment)
    if (dual_id == local_id) cycle
    if (well_connections(dual_segment,isegment)) cycle
    mask = PETSC_FALSE
    where (connections(1,:) == local_id)
      where (connections(2,:) == dual_id)
        mask = PETSC_TRUE
      endwhere
    endwhere
    if (any(mask)) cycle
    where (connections(1,:) == dual_id)
      where (connections(2,:) == local_id)
        mask = PETSC_TRUE
      endwhere
    endwhere
    if (any(mask)) cycle
    n_well_conn = n_well_conn+1
    well_connections(isegment,dual_segment) = PETSC_TRUE
  enddo

  total_connections = num_connections + n_well_conn
  allocate(new_connections(2,total_connections))
  allocate(new_face_areas(total_connections))
  allocate(new_face_centroids(total_connections))
  new_connections(:,1:num_connections) = connections(:,:)
  new_face_areas(1:num_connections) = face_areas(:)
  new_face_centroids(1:num_connections) = face_centroids(:)
  iconn = total_connections - n_well_conn + 1
  do isegment = 1,well_grid%nsegments
    dual_segment = 1
    if (well_connections(isegment,dual_segment)) then
      new_connections(1,iconn) = isegment
      new_connections(2,iconn) = dual_segment
      ! Create fake faces
      new_face_areas(iconn) = 0.d0
      new_face_centroids(iconn) = &
            cell_centroids(well_grid%h_ghosted_id(dual_segment))
      iconn = iconn + 1
    endif
  enddo

  nullify(connections)
  nullify(face_areas)
  nullify(face_centroids)
  connections => new_connections
  face_areas => new_face_areas
  face_centroids => new_face_centroids

  deallocate(mask)
  deallocate(well_connections)

end subroutine WellGridAddConnectionsExplicit

! ************************************************************************** !

end module Well_Grid_module
