module Inversion_Coupled_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module
  use Inversion_Parameter_module

  implicit none

  private

  type, public :: inversion_coupled_aux_type
    type(inversion_coupled_soln_type), pointer :: solutions(:)
  end type inversion_coupled_aux_type

  type, public :: inversion_coupled_soln_type
    PetscReal :: time
    PetscBool :: measured
    Vec :: original_saturation_solution
    Vec :: perturbed_saturation_solution
    Vec :: original_solute_solution
    Vec :: perturbed_solute_solution
    Vec, pointer :: dsaturation_dparameter(:)
    Vec, pointer :: dsolute_dparameter(:)
  end type inversion_coupled_soln_type

  public :: InversionCoupledAuxCreate, &
            InversionCoupledSolutionCreate, &
            InversionCoupledSolutionInit, &
            InvCoupledAllocateSolnVecs, &
            InvCoupledUpdateSolnVecs, &
            InversionCoupledAuxReset, &
            InversionCoupledAuxDestroy

contains

! ************************************************************************** !

function InversionCoupledAuxCreate()
  !
  ! Allocate and initialize auxiliary object for coupled flow and ert
  !
  ! Author: Glenn Hammond
  ! Date: 09/28/22
  !
  implicit none

  type(inversion_coupled_aux_type), pointer :: InversionCoupledAuxCreate

  type(inversion_coupled_aux_type), pointer :: aux

  allocate(aux)

  nullify(aux%solutions)

  InversionCoupledAuxCreate => aux

end function InversionCoupledAuxCreate

! ************************************************************************** !

function InversionCoupledSolutionCreate()
  !
  ! Allocate and initialize solution object for coupled flow and ert
  !
  ! Author: Glenn Hammond
  ! Date: 09/28/22
  !
  implicit none

  type(inversion_coupled_soln_type), pointer :: InversionCoupledSolutionCreate

  type(inversion_coupled_soln_type), pointer :: aux

  allocate(aux)
  call InversionCoupledSolutionInit(aux)

  InversionCoupledSolutionCreate => aux

end function InversionCoupledSolutionCreate

! ************************************************************************** !

subroutine InversionCoupledSolutionInit(aux)
  !
  ! Allocate and initialize solution object for coupled flow and ert
  !
  ! Author: Glenn Hammond
  ! Date: 09/28/22
  !
  implicit none

  type(inversion_coupled_soln_type) :: aux

  aux%time = UNINITIALIZED_DOUBLE
  aux%measured = PETSC_FALSE
  aux%original_saturation_solution = PETSC_NULL_VEC
  aux%perturbed_saturation_solution = PETSC_NULL_VEC
  aux%original_solute_solution = PETSC_NULL_VEC
  aux%perturbed_solute_solution = PETSC_NULL_VEC
  nullify(aux%dsaturation_dparameter)
  nullify(aux%dsolute_dparameter)

end subroutine InversionCoupledSolutionInit

! ************************************************************************** !

subroutine InvCoupledAllocateSolnVecs(aux,onedof_vec,num_parameters)
  !
  ! Allocate and initialize solution object for coupled flow and ert
  !
  ! Author: Glenn Hammond
  ! Date: 09/28/22
  !
  use ZFlow_Aux_module

  implicit none

  type(inversion_coupled_aux_type) :: aux
  Vec :: onedof_vec
  PetscInt :: num_parameters

  PetscInt :: i, j
  PetscErrorCode :: ierr

  do i = 1, size(aux%solutions)
    call VecDuplicate(onedof_vec,aux%solutions(i)%original_saturation_solution, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(onedof_vec,aux%solutions(i)%perturbed_saturation_solution, &
                      ierr);CHKERRQ(ierr)
    if (Initialized(zflow_sol_tran_eq)) then
      call VecDuplicate(onedof_vec,aux%solutions(i)%original_solute_solution, &
                        ierr);CHKERRQ(ierr)
      call VecDuplicate(onedof_vec,aux%solutions(i)%perturbed_solute_solution, &
                        ierr);CHKERRQ(ierr)
    endif
    allocate(aux%solutions(i)%dsaturation_dparameter(num_parameters))
    aux%solutions(i)%dsaturation_dparameter(:) = PETSC_NULL_VEC
    call VecDuplicateVecsF90(onedof_vec,num_parameters, &
                             aux%solutions(i)%dsaturation_dparameter, &
                             ierr);CHKERRQ(ierr)
    if (Initialized(zflow_sol_tran_eq)) then
      allocate(aux%solutions(i)%dsolute_dparameter(num_parameters))
      aux%solutions(i)%dsolute_dparameter(:) = PETSC_NULL_VEC
      call VecDuplicateVecsF90(onedof_vec,num_parameters, &
                              aux%solutions(i)%dsolute_dparameter, &
                              ierr);CHKERRQ(ierr)
    endif
    do j = 1, num_parameters
      call VecSet(aux%solutions(i)%dsaturation_dparameter(j), &
                  UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
      if (Initialized(zflow_sol_tran_eq)) then
        call VecSet(aux%solutions(i)%dsolute_dparameter(j), &
                    UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
      endif
    enddo
  enddo

end subroutine InvCoupledAllocateSolnVecs

! ************************************************************************** !

subroutine InvCoupledUpdateSolnVecs(iparameter,perturbed_solution, &
                                    original_solution,solution_derivatives, &
                                    pert)
  !
  ! Allocate and initialize solution object for coupled flow and ert
  !
  ! Author: Glenn Hammond
  ! Date: 09/28/22
  !
  implicit none

  PetscInt :: iparameter
  Vec :: perturbed_solution
  Vec :: original_solution
  Vec :: solution_derivatives(:)
  PetscReal :: pert

  PetscErrorCode :: ierr

  call VecWAXPY(solution_derivatives(iparameter),-1.d0, &
                original_solution,perturbed_solution,ierr);CHKERRQ(ierr)
  call VecScale(solution_derivatives(iparameter),1.d0/pert, &
                ierr);CHKERRQ(ierr)
  call VecSet(perturbed_solution,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)

end subroutine InvCoupledUpdateSolnVecs

! ************************************************************************** !

subroutine InversionCoupledAuxReset(aux)
  !
  ! Resets flags for auxiliary object for coupled flow and ert
  !
  ! Author: Glenn Hammond
  ! Date: 09/28/22

  type(inversion_coupled_aux_type), pointer :: aux

  PetscInt :: i

!  aux%cur_solution => aux%solution_list

  do i = 1, size(aux%solutions)
    aux%solutions(i)%measured = PETSC_FALSE
  enddo

end subroutine InversionCoupledAuxReset

! ************************************************************************** !

subroutine InversionCoupledAuxDestroy(aux)
  !
  ! Deallocates a auxiliary object for coupled flow and ert
  !
  ! Author: Glenn Hammond
  ! Date: 09/28/22
  !
  type(inversion_coupled_aux_type), pointer :: aux

  PetscInt :: i

  if (.not.associated(aux)) return

  if (associated(aux%solutions)) then
    do i = 1, size(aux%solutions)
      call InversionCoupledSolutionDestroy(aux%solutions(i))
    enddo
    deallocate(aux%solutions)
    nullify(aux%solutions)
  endif


  deallocate(aux)
  nullify(aux)

end subroutine InversionCoupledAuxDestroy

! ************************************************************************** !

subroutine InversionCoupledSolutionDestroy(solution)
  !
  ! Deallocates a auxiliary object for coupled flow and ert
  !
  ! Author: Glenn Hammond
  ! Date: 09/28/22
  !
  type(inversion_coupled_soln_type) :: solution

  PetscBool :: deallocate_solute
  PetscInt :: i
  PetscErrorCode :: ierr

  deallocate_solute = solution%original_solute_solution /= PETSC_NULL_VEC

  call VecDestroy(solution%original_saturation_solution,ierr);CHKERRQ(ierr)
  call VecDestroy(solution%perturbed_saturation_solution,ierr);CHKERRQ(ierr)
  ! solute is optional
  if (deallocate_solute) then
    call VecDestroy(solution%original_solute_solution,ierr);CHKERRQ(ierr)
    call VecDestroy(solution%perturbed_solute_solution,ierr);CHKERRQ(ierr)
  endif

  do i = 1, size(solution%dsaturation_dparameter)
    call VecDestroy(solution%dsaturation_dparameter(i),ierr);CHKERRQ(ierr)
    if (deallocate_solute) then
      call VecDestroy(solution%dsolute_dparameter(i),ierr);CHKERRQ(ierr)
    endif
  enddo
  nullify(solution%dsaturation_dparameter)
  nullify(solution%dsolute_dparameter)

end subroutine InversionCoupledSolutionDestroy

end module Inversion_Coupled_Aux_module
