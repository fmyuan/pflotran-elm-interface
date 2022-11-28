module Inversion_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module
  use Inversion_Coupled_Aux_module
  use Inversion_Measurement_Aux_module
  use Inversion_Parameter_module
  use Inversion_TS_Aux_module
  use Option_Inversion_module

  implicit none

  private

  type, public :: inversion_aux_type
    Vec :: solution ! solely a pointer
    PetscInt :: isync_time              ! current index of sync_times
    PetscReal, pointer :: sync_times(:) ! an array with all measurement times
    type(inversion_coupled_aux_type), pointer :: coupled_aux
    type(inversion_measurement_aux_type), pointer :: measurements(:)
    type(inversion_parameter_type), pointer :: parameters(:)
    Vec :: measurement_vec
    Vec :: dist_measurement_vec
    Vec :: parameter_vec
    Vec :: dist_parameter_vec
    VecScatter :: scatter_measure_to_dist_measure
    VecScatter :: scatter_param_to_dist_param
    VecScatter :: scatter_global_to_dist_param
    Mat :: JsensitivityT
    ! adjoint data structures
    PetscBool :: store_adjoint
    Mat :: M_ptr
    type(inversion_forward_ts_aux_type), pointer :: first_forward_ts_aux
    type(inversion_forward_ts_aux_type), pointer :: last_forward_ts_aux
    PetscReal, pointer :: local_measurement_values_ptr(:)
    PetscReal, pointer :: local_dobs_dunknown_values_ptr(:)
    PetscReal, pointer :: local_dobs_dparam_values_ptr(:)
  end type inversion_aux_type

  public :: InversionAuxCreate, &
            InversionAuxResetMeasurements, &
            InversionAuxAdjointRecordTS, &
            InvAuxAdjCleanupAfterForwardRun, &
            InversionAuxDestroy

contains

! ************************************************************************** !

function InversionAuxCreate()
  !
  ! Allocate and initialize auxiliary inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  type(inversion_aux_type), pointer :: InversionAuxCreate

  type(inversion_aux_type), pointer :: aux

  allocate(aux)

  aux%solution = PETSC_NULL_VEC

  aux%isync_time = 1
  nullify(aux%sync_times)
  nullify(aux%coupled_aux)
  nullify(aux%measurements)
  nullify(aux%parameters)
  aux%measurement_vec = PETSC_NULL_VEC
  aux%dist_measurement_vec = PETSC_NULL_VEC
  aux%parameter_vec = PETSC_NULL_VEC
  aux%dist_parameter_vec = PETSC_NULL_VEC
  aux%scatter_measure_to_dist_measure = PETSC_NULL_VECSCATTER
  aux%scatter_param_to_dist_param = PETSC_NULL_VECSCATTER
  aux%scatter_global_to_dist_param = PETSC_NULL_VECSCATTER
  aux%JsensitivityT = PETSC_NULL_MAT
  nullify(aux%local_measurement_values_ptr)
  nullify(aux%local_dobs_dunknown_values_ptr)
  nullify(aux%local_dobs_dparam_values_ptr)
  ! adjoint
  call InversionAuxInitAdjoint(aux)

  InversionAuxCreate => aux

end function InversionAuxCreate

! ************************************************************************** !

subroutine InversionAuxInitAdjoint(aux)
  !
  ! Initializes adjoint portion of object
  !
  ! Author: Glenn Hammond
  ! Date: 11/28/22

  type(inversion_aux_type) :: aux

  aux%store_adjoint = PETSC_TRUE
  aux%M_ptr = PETSC_NULL_MAT
  nullify(aux%first_forward_ts_aux)
  nullify(aux%last_forward_ts_aux)

end subroutine InversionAuxInitAdjoint

! ************************************************************************** !

subroutine InversionAuxResetMeasurements(aux)
  !
  ! Resets flags for forward run back to original settings.
  !
  ! Author: Glenn Hammond
  ! Date: 02/21/22

  type(inversion_aux_type), pointer :: aux

  PetscInt :: imeasurement

  aux%isync_time = 1
  do imeasurement = 1, size(aux%measurements)
    call InversionMeasurementAuxReset(aux%measurements(imeasurement))
  enddo

end subroutine InversionAuxResetMeasurements

! ************************************************************************** !

subroutine InversionAuxAdjointRecordTS(aux,time)
  !
  ! Appends a time step to the linked list
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/22, 11/28/22

  use Utility_module

  implicit none

  type(inversion_aux_type), pointer :: aux
  PetscReal :: time

  if (associated(aux%last_forward_ts_aux)) then
    aux%last_forward_ts_aux%time = time
    ! store the solution
    call InvForTSAuxDupForwardJacobian(aux%M_ptr,aux%last_forward_ts_aux)
    ! append next time step
    aux%last_forward_ts_aux => &
      InversionForwardTSAuxCreate(aux%last_forward_ts_aux)
  endif

end subroutine InversionAuxAdjointRecordTS

! ************************************************************************** !

subroutine InvAuxAdjCleanupAfterForwardRun(aux)
  !
  ! Destroys the linked list of adjoint objects
  !
  ! Author: Glenn Hammond
  ! Date: 11/28/22

  implicit none

  type(inversion_aux_type) :: aux

  ! destroy the lists
  call InvForwardTSAuxDestroyList(aux%first_forward_ts_aux,PETSC_FALSE)
  ! initialize everything else
  call InversionAuxInitAdjoint(aux)

end subroutine InvAuxAdjCleanupAfterForwardRun

! ************************************************************************** !

subroutine InversionAuxDestroy(aux)
  !
  ! Deallocates a inversion auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  use Utility_module, only : DeallocateArray

  type(inversion_aux_type), pointer :: aux

  PetscInt :: i
  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

  call DeallocateArray(aux%sync_times)

  ! these are owned and must be destroyed
  if (associated(aux%coupled_aux)) then
    call InversionCoupledAuxDestroy(aux%coupled_aux)
  endif
  if (associated(aux%measurements)) then
    do i = 1, size(aux%measurements)
      call InversionMeasurementAuxStrip(aux%measurements(i))
    enddo
    deallocate(aux%measurements)
  endif
  nullify(aux%measurements)
  if (associated(aux%parameters)) then
    do i = 1, size(aux%parameters)
      call InversionParameterStrip(aux%parameters(i))
    enddo
    deallocate(aux%parameters)
  endif
  nullify(aux%parameters)
  if (aux%measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(aux%measurement_vec,ierr);CHKERRQ(ierr)
  endif
  if (aux%dist_measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(aux%dist_measurement_vec,ierr);CHKERRQ(ierr)
  endif
  if (aux%parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(aux%parameter_vec,ierr);CHKERRQ(ierr)
  endif
  if (aux%dist_parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(aux%dist_parameter_vec,ierr);CHKERRQ(ierr)
  endif
  if (aux%scatter_measure_to_dist_measure /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(aux%scatter_measure_to_dist_measure, &
                           ierr);CHKERRQ(ierr)
  endif
  if (aux%scatter_param_to_dist_param /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(aux%scatter_measure_to_dist_measure, &
                           ierr);CHKERRQ(ierr)
  endif
  if (aux%scatter_global_to_dist_param /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(aux%scatter_measure_to_dist_measure, &
                           ierr);CHKERRQ(ierr)
  endif
  if (aux%JsensitivityT /= PETSC_NULL_MAT) then
    call MatDestroy(aux%JsensitivityT,ierr);CHKERRQ(ierr)
  endif

  ! nullify objects owned by other objects
  aux%solution = PETSC_NULL_VEC
  nullify(aux%local_measurement_values_ptr)
  nullify(aux%local_dobs_dunknown_values_ptr)
  nullify(aux%local_dobs_dparam_values_ptr)
  ! adjoints
  call InvAuxAdjCleanupAfterForwardRun(aux)

  deallocate(aux)
  nullify(aux)

end subroutine InversionAuxDestroy

end module Inversion_Aux_module
