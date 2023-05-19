module Realization_Common_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Patch_module
  use Region_module

  implicit none

  public :: RealizationLocalizeRegions, &
            RealizationAddCoupler, &
            RealizationAddStrata


contains

! ************************************************************************** !

subroutine RealizationLocalizeRegions(patch,region_list,option)
  !
  ! Localizes regions within each patch
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use String_module
  use Grid_module

  implicit none

  type(patch_type), pointer :: patch
  type(region_list_type), pointer :: region_list
  type(option_type), pointer :: option

  type (region_type), pointer :: cur_region, cur_region2
  type(region_type), pointer :: region

  ! check to ensure that region names are not duplicated
  cur_region => region_list%first
  do
    if (.not.associated(cur_region)) exit
    cur_region2 => cur_region%next
    do
      if (.not.associated(cur_region2)) exit
      if (StringCompare(cur_region%name,cur_region2%name,MAXWORDLENGTH)) then
        option%io_buffer = 'Duplicate region names: ' // trim(cur_region%name)
        call PrintErrMsg(option)
      endif
      cur_region2 => cur_region2%next
    enddo
    cur_region => cur_region%next
  enddo

  call PatchLocalizeRegions(patch,region_list,option)
  ! destroy realization's copy of region list as it can be confused with the
  ! localized patch regions later in teh simulation.
  call RegionDestroyList(region_list)

  ! compute regional connections for inline surface flow
  if (option%flow%inline_surface_flow) then
     region => RegionGetPtrFromList(option%flow%inline_surface_region_name, &
          patch%region_list)
     if (.not.associated(region)) then
        option%io_buffer = 'realization_subsurface.F90:RealizationLocalize&
             &Regions() --> Could not find a required region named "' // &
             trim(option%flow%inline_surface_region_name) // &
             '" from the list of regions.'
        call PrintErrMsg(option)
     endif
     call GridRestrictRegionalConnect(patch%grid,region)
   endif

end subroutine RealizationLocalizeRegions


! ************************************************************************** !

subroutine RealizationAddCoupler(patch,coupler)
  !
  ! Adds a copy of a coupler to a list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Coupler_module

  implicit none

  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: coupler

  type(coupler_type), pointer :: new_coupler

  ! only add to flow list for now, since they will be split out later
  new_coupler => CouplerCreate(coupler)
  select case(coupler%itype)
    case(BOUNDARY_COUPLER_TYPE)
      call CouplerAddToList(new_coupler,patch%boundary_condition_list)
    case(INITIAL_COUPLER_TYPE)
      call CouplerAddToList(new_coupler,patch%initial_condition_list)
    case(SRC_SINK_COUPLER_TYPE)
      call CouplerAddToList(new_coupler,patch%source_sink_list)
  end select
  nullify(new_coupler)

  call CouplerDestroy(coupler)

end subroutine RealizationAddCoupler

! ************************************************************************** !

subroutine RealizationAddStrata(patch,strata)
  !
  ! Adds a copy of a strata to a list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Strata_module

  implicit none

  type(patch_type), pointer :: patch
  type(strata_type), pointer :: strata

  type(strata_type), pointer :: new_strata

  new_strata => StrataCreate(strata)
  call StrataAddToList(new_strata,patch%strata_list)
  nullify(new_strata)

  call StrataDestroy(strata)

end subroutine RealizationAddStrata

end module Realization_Common_module