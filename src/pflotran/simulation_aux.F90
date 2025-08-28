module Simulation_Aux_module

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PFLOTRAN_Constants_module

  implicit none

  private

  type,public :: simulation_aux_type

    ! Note: These are GLOBAL vectors (i.e. they do not contain ghost control
    !       volumes)

    ! Size of entire subsurface domain
    Vec :: subsurf_pres
    Vec :: subsurf_temp
    Vec :: subsurf_fluid_den

    Vec :: subsurf_por0
    Vec :: subsurf_por
    Vec :: subsurf_strain
    Vec :: subsurf_stress
    Vec :: subsurf_perm0
    Vec :: subsurf_perm

    VecScatter :: geomechanics_to_subsurf
    VecScatter :: subsurf_to_geomechanics

  end type simulation_aux_type

  public :: SimAuxCreate, &
            SimAuxCopyVecScatter, &
            SimAuxCopySubsurfVec, &
            SimAuxCopySubsurfGeomechVec, &
            SimAuxDestroy

contains

! ************************************************************************** !

function SimAuxCreate()
  !
  ! This routine allocates auxillary object.
  !
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/20/13
  !

  use Option_module

  implicit none

  type (simulation_aux_type),pointer :: SimAuxCreate

  type (simulation_aux_type),pointer :: aux

  allocate(aux)
  PetscObjectNullify(aux%subsurf_pres)
  PetscObjectNullify(aux%subsurf_temp)
  PetscObjectNullify(aux%subsurf_fluid_den)
  PetscObjectNullify(aux%subsurf_por0)
  PetscObjectNullify(aux%subsurf_por)
  PetscObjectNullify(aux%subsurf_perm0)
  PetscObjectNullify(aux%subsurf_perm)
  PetscObjectNullify(aux%subsurf_strain)
  PetscObjectNullify(aux%subsurf_stress)

  PetscObjectNullify(aux%subsurf_to_geomechanics)
  PetscObjectNullify(aux%geomechanics_to_subsurf)

  SimAuxCreate => aux

end function SimAuxCreate

! ************************************************************************** !

subroutine SimAuxCopyVecScatter(aux, vscat, vscat_index)
  !
  ! This routine copies VectorScatter to an appropriate context.
  !
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  !

  implicit none

  type (simulation_aux_type),pointer :: aux
  VecScatter :: vscat
  PetscInt :: vscat_index

  PetscErrorCode :: ierr

  select case (vscat_index)
    case(SUBSURF_TO_GEOMECHANICS)
      call VecScatterCopy(vscat,aux%subsurf_to_geomechanics, &
                          ierr);CHKERRQ(ierr)
    case(GEOMECHANICS_TO_SUBSURF)
      call VecScatterCopy(vscat,aux%geomechanics_to_subsurf, &
                          ierr);CHKERRQ(ierr)
  end select

end subroutine SimAuxCopyVecScatter

! ************************************************************************** !

subroutine SimAuxCopySubsurfVec(aux, subsurf_vec)
  !
  ! This routine creates 3D vectors related with subsurface-flow.
  !
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  !

  implicit none

  type (simulation_aux_type),pointer :: aux
  Vec :: subsurf_vec

  PetscErrorCode :: ierr

  call VecDuplicate(subsurf_vec,aux%subsurf_pres,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_temp,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_fluid_den,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_por0,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_por,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_perm0,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_perm,ierr);CHKERRQ(ierr)

end subroutine SimAuxCopySubsurfVec

! ************************************************************************** !

subroutine SimAuxCopySubsurfGeomechVec(aux, subsurf_geomech_vec)
  !
  ! This routine creates vectors associated with geomechanics.
  !
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  !
  implicit none

  type (simulation_aux_type),pointer :: aux
  Vec :: subsurf_geomech_vec

  PetscErrorCode :: ierr

  call VecDuplicate(subsurf_geomech_vec,aux%subsurf_stress, &
                    ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_geomech_vec,aux%subsurf_strain, &
                    ierr);CHKERRQ(ierr)

end subroutine SimAuxCopySubsurfGeomechVec

! ************************************************************************** !

subroutine SimAuxDestroy(aux)
  !
  ! This routine deallocates auxillary object.
  !
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/20/13
  !

  implicit none

  type(simulation_aux_type), pointer :: aux

  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

  if (.not.PetscObjectIsNull(aux%subsurf_pres)) then
    call VecDestroy(aux%subsurf_pres,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(aux%subsurf_temp)) then
    call VecDestroy(aux%subsurf_temp,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(aux%subsurf_fluid_den)) then
    call VecDestroy(aux%subsurf_fluid_den,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(aux%subsurf_por0)) then
    call VecDestroy(aux%subsurf_por0,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(aux%subsurf_por)) then
    call VecDestroy(aux%subsurf_por,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(aux%subsurf_perm0)) then !DANNY
    call VecDestroy(aux%subsurf_perm0,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(aux%subsurf_perm)) then
    call VecDestroy(aux%subsurf_perm,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(aux%subsurf_stress)) then
    call VecDestroy(aux%subsurf_stress,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(aux%subsurf_strain)) then
    call VecDestroy(aux%subsurf_strain,ierr);CHKERRQ(ierr)
  endif

  if (.not.PetscObjectIsNull(aux%subsurf_to_geomechanics)) then
    call VecScatterDestroy(aux%subsurf_to_geomechanics,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(aux%geomechanics_to_subsurf)) then
    call VecScatterDestroy(aux%geomechanics_to_subsurf,ierr);CHKERRQ(ierr)
  endif

  deallocate(aux)
  nullify(aux)

end subroutine SimAuxDestroy

end module Simulation_Aux_module
