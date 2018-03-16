module Solver_cpr_module

#include "petsc/finclude/petscts.h"
  use petscts
  !use Solver_module
  use CPR_Precondititioner_module 

  implicit none

  private

  public :: SolverCPRInit


contains

subroutine SolverCPRInit(J, stash, pcin, ierr, option)
  use Option_module
  implicit none
  !type(solver_type) :: solver
  Mat :: J
  type(cpr_pc_type) :: stash
  PC :: pcin
  MPI_Comm :: C
  PetscErrorCode :: ierr
  type(option_type) :: option

  call PetscObjectGetComm(pcin, C, ierr); CHKERRQ(ierr)

  call CPRmake(pcin, stash, C, ierr, option)

  !! set the A matrix in the stash to J:
  stash%A = J

end subroutine SolverCPRInit

end module Solver_cpr_module
