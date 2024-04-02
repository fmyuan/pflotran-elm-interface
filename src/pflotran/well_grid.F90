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
    ! delta h discretization of each segment center [m]
    PetscReal, pointer :: dh(:)
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
end type well_grid_type

public :: WellGridCreate, &
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
    nullify(well_grid%dh)
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
  
    WellGridCreate => well_grid
  
end function WellGridCreate

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

end module Well_Grid_module