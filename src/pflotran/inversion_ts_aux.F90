module Inversion_TS_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module
  use Inversion_Measurement_Aux_module

  implicit none

  private

  type, public :: inversion_forward_aux_type
    PetscBool :: store_adjoint
    PetscInt :: num_timesteps
    PetscInt :: iobs_var
    PetscInt :: isync_time
    PetscReal, pointer :: sync_times(:)
    Mat :: M_ptr
    Vec :: solution_ptr
    type(inversion_forward_ts_aux_type), pointer :: first
    type(inversion_forward_ts_aux_type), pointer :: last
    type(inversion_forward_ts_aux_type), pointer :: current
    type(inversion_measurement_aux_type), pointer :: measurements(:)
    Vec :: measurement_vec
    PetscReal, pointer :: local_measurement_values_ptr(:)
    PetscReal, pointer :: local_derivative_values_ptr(:)
  end type inversion_forward_aux_type

  type, public :: inversion_forward_ts_aux_type
    PetscInt :: timestep
    PetscReal :: time
    Mat :: dResdu             ! copy of Jacobian: df(u)/du
    Mat :: dResdparam         ! matrix storing df(u)/dparameter
    Vec, pointer :: lambda(:) ! arrays storing dg(u)/df(u)
    ! derivative of residual wrt unknown at old time level (k)
    PetscReal, pointer :: dRes_du_k(:,:,:)  ! array storing df(u)^k+1/du^k
    PetscReal, pointer :: dFluxdIntConn(:,:)
    PetscReal, pointer :: dFluxdBCConn(:,:)
    Vec :: solution
    type(inversion_forward_ts_aux_type), pointer :: prev
    type(inversion_forward_ts_aux_type), pointer :: next
  end type inversion_forward_ts_aux_type


  public :: InversionForwardAuxCreate, &
            InvForwardAuxResetMeasurements, &
            InversionForwardAuxStep, &
            InvForwardAuxDestroyList, &
            InversionForwardAuxDestroy, &
            InversionTSAuxCreate, &
            InvTSAuxAllocate, &
            InvTSAuxAllocateFluxCoefArrays, &
            InvTSAuxStoreCopyGlobalMatVecs, &
            InversionTSAuxDestroy

contains

! ************************************************************************** !

function InversionForwardAuxCreate()
  !
  ! Allocate and initialize auxiliary inversion forward object
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/22
  !
  implicit none

  type(inversion_forward_aux_type), pointer :: InversionForwardAuxCreate

  type(inversion_forward_aux_type), pointer :: aux

  allocate(aux)

  aux%store_adjoint = PETSC_TRUE
  aux%num_timesteps = 0
  aux%iobs_var = UNINITIALIZED_INTEGER
  aux%isync_time = 1
  nullify(aux%sync_times)
  aux%M_ptr = PETSC_NULL_MAT
  aux%solution_ptr = PETSC_NULL_VEC
  nullify(aux%first)
  nullify(aux%last)
  nullify(aux%current)
  nullify(aux%measurements)
  aux%measurement_vec = PETSC_NULL_VEC
  nullify(aux%local_measurement_values_ptr)
  nullify(aux%local_derivative_values_ptr)

  InversionForwardAuxCreate => aux

end function InversionForwardAuxCreate

! ************************************************************************** !

subroutine InvForwardAuxResetMeasurements(aux)
  !
  ! Resets flags for forward run back to original settings.
  !
  ! Author: Glenn Hammond
  ! Date: 02/21/22

  type(inversion_forward_aux_type), pointer :: aux

  PetscInt :: imeasurement

  aux%isync_time = 1
  do imeasurement = 1, size(aux%measurements)
    call InversionMeasurementAuxReset(aux%measurements(imeasurement))
  enddo

end subroutine InvForwardAuxResetMeasurements

! ************************************************************************** !

subroutine InversionForwardAuxStep(aux,time)
  !
  ! Appends a time step to the linked list
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/22

  use Utility_module

  implicit none

  type(inversion_forward_aux_type), pointer :: aux
  PetscReal :: time

  if (associated(aux%current)) then
    aux%current%time = time
    ! store the solution
    call InvTSAuxStoreCopyGlobalMatVecs(aux,aux%current)
    ! append next time step
    aux%current => InversionTSAuxCreate(aux%current,aux%M_ptr)
    aux%last => aux%current
  endif

end subroutine InversionForwardAuxStep

! ************************************************************************** !

function InversionTSAuxCreate(prev_ts_aux,M_ptr)
  !
  ! Allocate and initialize auxiliary inversion time step object
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21
  !
  implicit none

  type(inversion_forward_ts_aux_type), pointer :: prev_ts_aux
  Mat :: M_ptr

  type(inversion_forward_ts_aux_type), pointer :: InversionTSAuxCreate

  type(inversion_forward_ts_aux_type), pointer :: aux

  allocate(aux)

  aux%time = UNINITIALIZED_DOUBLE
  aux%dResdu = PETSC_NULL_MAT
  aux%dResdparam = PETSC_NULL_MAT
  aux%solution = PETSC_NULL_VEC
  nullify(aux%lambda)

  nullify(aux%dRes_du_k)
  nullify(aux%dFluxdIntConn)
  nullify(aux%dFluxdBCConn)
  nullify(aux%prev)
  nullify(aux%next)

  if (associated(prev_ts_aux)) then
    aux%timestep = prev_ts_aux%timestep + 1
    prev_ts_aux%next => aux
    aux%prev => prev_ts_aux
    call InvTSAuxAllocate(aux,M_ptr,size(prev_ts_aux%dRes_du_k,1), &
                          size(prev_ts_aux%dRes_du_k,3))
    if (associated(prev_ts_aux%dFluxdIntConn)) then
      call InvTSAuxAllocateFluxCoefArrays(aux, &
                                          size(prev_ts_aux%dFluxdIntConn,2), &
                                          size(prev_ts_aux%dFluxdBCConn,2))
    endif
  else
    aux%timestep = 1
  endif

  InversionTSAuxCreate => aux

end function InversionTSAuxCreate

! ************************************************************************** !

subroutine InvTSAuxAllocate(aux,M_ptr,ndof,ncell)
  !
  ! Allocated array holding Jacobian and dResdk matrix coefficients
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21

  type(inversion_forward_ts_aux_type), pointer :: aux
  Mat :: M_ptr
  PetscInt :: ndof
  PetscInt :: ncell
  PetscErrorCode :: ierr

  call MatDuplicate(M_ptr,MAT_SHARE_NONZERO_PATTERN,aux%dResdparam, &
                    ierr);CHKERRQ(ierr)
  allocate(aux%dRes_du_k(ndof,ndof,ncell))
  aux%dRes_du_k = 0.d0

end subroutine InvTSAuxAllocate

! ************************************************************************** !

subroutine InvTSAuxAllocateFluxCoefArrays(aux,num_internal,num_boundary)
  !
  ! Allocated array holding flux coefficients
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21

  type(inversion_forward_ts_aux_type), pointer :: aux
  PetscInt :: num_internal
  PetscInt :: num_boundary

  allocate(aux%dFluxdIntConn(6,num_internal))
  aux%dFluxdIntConn = 0.d0
  allocate(aux%dFluxdBCConn(2,num_boundary))
  aux%dFluxdBCConn = 0.d0

end subroutine InvTSAuxAllocateFluxCoefArrays

! ************************************************************************** !

subroutine InvTSAuxStoreCopyGlobalMatVecs(forward_aux,ts_aux)
  !
  ! Copies Jacobian matrix and solution vector
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/22

  type(inversion_forward_aux_type), pointer :: forward_aux
  type(inversion_forward_ts_aux_type), pointer :: ts_aux

  PetscErrorCode :: ierr

  call MatDuplicate(forward_aux%M_ptr,MAT_COPY_VALUES,ts_aux%dResdu, &
                    ierr);CHKERRQ(ierr)
  if (forward_aux%solution_ptr /= PETSC_NULL_VEC) then
    call VecDuplicate(forward_aux%solution_ptr,ts_aux%solution, &
                      ierr);CHKERRQ(ierr)
    call VecCopy(forward_aux%solution_ptr,ts_aux%solution,ierr);CHKERRQ(ierr)
  endif

end subroutine InvTSAuxStoreCopyGlobalMatVecs

! ************************************************************************** !

subroutine InversionForwardAuxDestroy(aux)
  !
  ! Deallocates a inversion forward auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/22
  !
  use Utility_module, only : DeallocateArray

  type(inversion_forward_aux_type), pointer :: aux

  call InvForwardAuxDestroyList(aux,PETSC_FALSE)

  call DeallocateArray(aux%sync_times)

  ! simply nullify
  nullify(aux%last)
  nullify(aux%current)
  aux%M_ptr = PETSC_NULL_MAT
  aux%solution_ptr = PETSC_NULL_VEC
  nullify(aux%measurements)
  aux%measurement_vec = PETSC_NULL_VEC
  nullify(aux%local_measurement_values_ptr)
  nullify(aux%local_derivative_values_ptr)

  deallocate(aux)
  nullify(aux)

end subroutine InversionForwardAuxDestroy

! ************************************************************************** !

subroutine InvForwardAuxDestroyList(aux,print_msg)
  !
  ! Deallocates a inversion auxiliary timestep list object
  !
  ! Author: Glenn Hammond
  ! Date: 12/03/21
  !
  type(inversion_forward_aux_type) :: aux
  PetscBool :: print_msg

  type(inversion_forward_ts_aux_type), pointer :: cur_inversion_ts_aux
  type(inversion_forward_ts_aux_type), pointer :: next_inversion_ts_aux

  cur_inversion_ts_aux => aux%first
  do
    if (.not.associated(cur_inversion_ts_aux)) exit
    next_inversion_ts_aux => cur_inversion_ts_aux%next
#if 0
    if (.not.associated(next_inversion_ts_aux)) then
      if (print_msg) then
        print *
        print *, '  Last Inversion TS Aux timestep and time: ', &
          cur_inversion_ts_aux%timestep, cur_inversion_ts_aux%time
        print *
      endif
    endif
#endif
    call InversionTSAuxDestroy(cur_inversion_ts_aux)
    cur_inversion_ts_aux => next_inversion_ts_aux
  enddo

  nullify(aux%first)
  nullify(aux%current)
  nullify(aux%last)

end subroutine InvForwardAuxDestroyList

! ************************************************************************** !

subroutine InversionTSAuxDestroy(aux)
  !
  ! Deallocates a inversion auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21
  !
  use Utility_module, only : DeallocateArray

  type(inversion_forward_ts_aux_type), pointer :: aux

  PetscInt :: iaux
  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

  nullify(aux%prev)
  nullify(aux%next)

  if (aux%dResdu /= PETSC_NULL_MAT) then
    call MatDestroy(aux%dResdu,ierr);CHKERRQ(ierr)
  endif
  if (aux%dResdparam /= PETSC_NULL_MAT) then
    call MatDestroy(aux%dResdparam,ierr);CHKERRQ(ierr)
  endif
  if (aux%solution /= PETSC_NULL_VEC) then
    call VecDestroy(aux%solution,ierr);CHKERRQ(ierr)
  endif
  if (associated(aux%lambda)) then
    call VecDestroyVecs(size(aux%lambda),aux%lambda,ierr);CHKERRQ(ierr)
    nullify(aux%lambda)
  endif
  call DeallocateArray(aux%dRes_du_k)
  call DeallocateArray(aux%dFluxdIntConn)
  call DeallocateArray(aux%dFluxdBCConn)

  deallocate(aux)
  nullify(aux)

end subroutine InversionTSAuxDestroy

end module Inversion_TS_Aux_module
