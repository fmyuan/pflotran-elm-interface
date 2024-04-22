module Init_Subsurface_Flow_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: InitSubsurfFlowReadInitCond

contains

! ************************************************************************** !

subroutine InitSubsurfFlowReadInitCond(realization,filename)
  !
  ! Assigns flow initial condition from HDF5 file
  !
  ! Author: Glenn Hammond
  ! Date: 03/05/10, 12/04/14
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
  use Discretization_module
  use HDF5_module

  implicit none

  class(realization_subsurface_type) :: realization
  character(len=MAXSTRINGLENGTH) :: filename

  PetscInt :: local_id, idx, offset
  PetscReal, pointer :: xx_p(:)
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscReal, pointer :: vec_p(:)
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  if (option%iflowmode /= RICHARDS_MODE .and. &
      option%iflowmode /= RICHARDS_TS_MODE) then
    option%io_buffer = 'Reading of flow initial conditions from HDF5 ' // &
                       'file (' // trim(filename) // &
                       'not currently not supported for mode: ' // &
                       trim(option%flowmode)
  endif

    ! assign initial conditions values to domain
  call VecGetArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

  ! Pressure for all modes
  offset = 1
  group_name = ''
  dataset_name = 'Pressure'
  call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                    filename,group_name, &
                                    dataset_name,option%id>0)
  call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
  do local_id=1, grid%nlmax
    if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    if (dabs(vec_p(local_id)) < 1.d-40) then
      print *,  option%myrank, grid%nG2A(grid%nL2G(local_id)), &
            ': Potential error - zero pressure in Initial Condition read from file.'
    endif
    idx = (local_id-1)*option%nflowdof + offset
    xx_p(idx) = vec_p(local_id)
  enddo
  call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                   field%flow_xx_loc,NFLOWDOF)
  call VecCopy(field%flow_xx,field%flow_yy,ierr);CHKERRQ(ierr)

end subroutine InitSubsurfFlowReadInitCond

end module Init_Subsurface_Flow_module
