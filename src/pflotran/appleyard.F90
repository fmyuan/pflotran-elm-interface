module Appleyard_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  implicit none

  private

  public :: TOilAppleyard, &
            TOWGAppleyard

contains

subroutine  TOWGAppleyard(saturation0, del_saturation, ghosted_id, realization, oid, changed)

  use Realization_Subsurface_class
  implicit none
  PetscReal :: saturation0, del_saturation
  PetscInt :: ghosted_id, oid
  !class(realization_subsurface_type), pointer :: realization
  !type(realization_subsurface_type), pointer :: realization
  type(realization_subsurface_type) :: realization

  PetscReal ::  ds_out_o, slc, soc
  PetscReal :: saturation0_oil, del_saturation_oil

  PetscBool :: changed

  changed = PETSC_FALSE


      !! 1) get residual (critical) saturations
      call GetCriticalSaturation(soc, oid, realization, ghosted_id)

      !! 2) perform the chop on each phase
      !!    oil:
      saturation0_oil = saturation0
      del_saturation_oil = del_saturation
      call AppleyardChopSuggest(saturation0_oil, del_saturation_oil, soc, ds_out_o)

      !! 3) update the saturation change if it needs to be changed:
      if (del_saturation_oil /= ds_out_o) then
        print *, "oil appleyard ", del_saturation_oil, " ", ds_out_o, " ", &
                 saturation0_oil - ds_out_o
      !call AppleyardChop(saturation0_oil, del_saturation_oil, soc, ds_out_o)

        del_saturation = ds_out_o
        changed = PETSC_TRUE
      endif

      !print *, "oil crit is ", soc
      !!! /end of Applyard chop
end subroutine TOWGAppleyard

! ************************************************************************** !


subroutine  TOilAppleyard(saturation0, del_saturation, ghosted_id, realization, lid, oid)

  use Realization_Subsurface_class
  implicit none
  PetscReal :: saturation0, del_saturation
  PetscInt :: ghosted_id, lid, oid
  class(realization_subsurface_type), pointer :: realization

  PetscReal :: ds_out_l, ds_out_o, slc, soc
  PetscReal :: saturation0_oil, del_saturation_oil, saturation0_liq, del_saturation_liq


      !! 1) get residual (critical) saturations
      call GetCriticalSaturation(slc, lid, realization, ghosted_id)
      call GetCriticalSaturation(soc, oid, realization, ghosted_id)

      !! 2) perform the chop on each phase
      !!    liquid:
      saturation0_liq = 1.d0 - saturation0
      del_saturation_liq = -1.d0 * del_saturation
      call AppleyardChopSuggest(saturation0_liq, del_saturation_liq, slc, ds_out_l)
      !!    oil:
      saturation0_oil = saturation0
      del_saturation_oil = del_saturation
      call AppleyardChopSuggest(saturation0_oil, del_saturation_oil, soc, ds_out_o)

      !! 3) update the saturation change if it needs to be changed:
      if (del_saturation_liq /= ds_out_l .AND. del_saturation_oil /= ds_out_o) then
        print *, "Appleyard chop has trucated both oil and liquid saturation. How?"
        !! do the biggest one:
        if (abs(ds_out_l) > abs(ds_out_o)) then
          del_saturation = -1.d0*ds_out_l
        else
          del_saturation = ds_out_o
        endif
      !! check liquid:
      elseif (del_saturation_liq /= ds_out_l) then
        !print *, "liquid appleyard ", del_saturation_liq, " ", ds_out_l, " ", &
                 !saturation0_liq - ds_out_l
      !call AppleyardChop(saturation0, del_saturation, slc, ds_out_l)

        del_saturation = -1.d0*ds_out_l
      !! check oil:
      elseif (del_saturation_oil /= ds_out_o) then
        !print *, "oil appleyard ", del_saturation_oil, " ", ds_out_o, " ", &
                 !saturation0_oil - ds_out_o
      !call AppleyardChop(saturation0_oil, del_saturation_oil, soc, ds_out_o)

        del_saturation = ds_out_o
      endif
      !!! /end of Applyard chop


      !!!!! what about default case?
end subroutine TOilAppleyard

! ************************************************************************** !

subroutine AppleyardChopSuggest(s, ds, sc, ds_out)

!!! note based on assumption that ds is NEGATIVE increment

  implicit none
  PetscReal :: s, ds, sc
  PetscReal :: ds_out

  PetscReal :: s_new
  PetscReal :: margin, scl, scu, passby

  !! increment is a subtraction:
  s_new = s - ds

  !margin = 2.5d-2
  !margin = 2.0d-3
  margin = 1.0d-6
  !margin = 0.d0 !! terrible idea

  scl = sc - margin
  scu = sc + margin

  ds_out = ds

  !passby = 0.d0
  passby = 0.5d0 * margin
  !passby = 1.1d0*margin

  if (abs(ds) <= 2.d0*margin)  then
    !print *, "not doing tiny appleyard"
    return
  endif

  if (s < scl .AND. s_new > scu) then
    !! crossed envelope going UP, now we want to 
    !! override such that:
    ! s_new = sc + margin/2.0
    ! i.e.,
    ! s - ds = sc + margin/2.0
    ! ds = s - sc - margin/2.0
    !ds_out = s - sc - margin/2.0
    ds_out = s - sc - passby
  endif

  if (s > scu .AND. s_new < scl) then
    !! crossed envelope going DOWN, now we want to 
    !! override such that:
    ! s_new = sc - margin/2.0
    ! i.e.,
    ! s - ds = sc - margin/2.0
    ! ds = s - sc + margin/2.0
    !ds_out = s - sc + margin/2.0
    ds_out = s - sc + passby
  endif

end subroutine AppleyardChopSuggest

! ************************************************************************** !

subroutine GetCriticalSaturation(sc, phase_id, realization, cell_id)

! pull out the residual situation for phase phase_id from the realization
! object.
!
! Dan Stone April 2018

  use Realization_Subsurface_class

  implicit none
  PetscReal :: sc
  PetscInt :: phase_id, cell_id
  !class(realization_subsurface_type), pointer :: realization
  !type(realization_subsurface_type), pointer :: realization
  type(realization_subsurface_type) :: realization


  ! we're only doing three things here:
  ! 1) drill down into:
  !    realization%patch%aux%Material%material_parameter%soil_residual_saturation(:,:)
  !    which is a 2d array: phase x (saturation function id)
  ! 2) To get the sat func id, drill down into:
  !    realization%patch%sat_func_id(cell_id)
  ! 3) address the correct entry of the soil residual saturation:
  !    soil_residual_saturation(phase_id, (soil_residual_saturation)) 

  sc = realization%patch%aux%Material%material_parameter%soil_residual_saturation( &
                                     phase_id, realization%patch%sat_func_id(cell_id))

end subroutine GetCriticalSaturation
! ************************************************************************** !

end module Appleyard_module
