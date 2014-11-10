module EPICMod

  use EPICConstants
  use EPICData
  use EPICMesh
  use Mapping_module

 implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type, public :: epic_type

    type(epic_mesh_list_type), pointer :: clm_mesh_list
    type(epic_mesh_list_type), pointer :: pflotran_mesh_list

    type(epic_data_list_type), pointer :: clm_data_list
    type(epic_data_list_type), pointer :: pflotran_data_list

    type(mapping_list_type), pointer :: map_list

    ! Number of cells for the 3D subsurface domain
    !PetscInt :: nlclm_sub ! num of local clm cells
    !PetscInt :: ngclm_sub ! num of ghosted clm cells (ghosted = local+ghosts)
    !PetscInt :: nlpf_sub  ! num of local pflotran cells
    !PetscInt :: ngpf_sub  ! num of ghosted pflotran cells (ghosted = local+ghosts)

    ! Number of cells for the surface of the 3D subsurface domain
    !PetscInt :: nlclm_2dsub  ! num of local clm cells
    !PetscInt :: ngclm_2dsub  ! num of ghosted clm cells (ghosted = local+ghosts)
    !PetscInt :: nlpf_2dsub   ! num of local pflotran cells
    !PetscInt :: ngpf_2dsub   ! num of ghosted pflotran cells (ghosted = local+ghosts)

    ! Number of cells for the 2D surface domain
    !PetscInt :: nlclm_srf  ! num of local clm cells
    !PetscInt :: ngclm_srf  ! num of ghosted clm cells (ghosted = local+ghosts)
    !PetscInt :: nlpf_srf   ! num of local pflotran cells
    !PetscInt :: ngpf_srf   ! num of ghosted pflotran cells (ghosted = local+ghosts)

    PetscInt :: nzclm_mapped ! num of CLM soil layers that are mapped

  end type epic_type

  type(epic_type) , public, target , save :: epic

  public :: EPCIInit, &
            EPICAddData

contains

! ************************************************************************** !

subroutine EPCIInit()

  implicit none

  epic%clm_data_list => EPICDataListCreate()
  epic%pflotran_data_list => EPICDataListCreate()

  epic%clm_mesh_list => EPICMeshListCreate()
  epic%pflotran_mesh_list => EPICMeshListCreate()

  !epic%nlclm_sub = 0
  !epic%ngclm_sub = 0
  !epic%nlpf_sub = 0
  !epic%ngpf_sub = 0

  !epic%nlclm_2dsub = 0
  !epic%ngclm_2dsub = 0
  !epic%nlpf_2dsub = 0
  !epic%ngpf_2dsub = 0

  !epic%nlclm_srf = 0
  !epic%ngclm_srf = 0
  !epic%nlpf_srf = 0
  !epic%ngpf_srf = 0

  epic%nzclm_mapped = 0

end subroutine EPCIInit

! ************************************************************************** !

function EPICGetPointerToModelList(model_type)

  implicit none

  PetscInt :: model_type

  type(epic_data_list_type), pointer :: model_data_list
  type(epic_data_list_type), pointer :: EPICGetPointerToModelList

  select case(model_type)
  case (EPIC_MODEL_CLM)
    model_data_list => epic%clm_data_list
  case (EPIC_MODEL_PFLOTRAN)
    model_data_list => epic%pflotran_data_list
  case default
    write(*,*) 'EPICGetPointerToModelList: Unknown model_type'
    stop
  end select

  EPICGetPointerToModelList => model_data_list

end function EPICGetPointerToModelList

! ************************************************************************** !

function EPICGetPointerToMeshList(model_type)

  implicit none

  PetscInt :: model_type

  type(epic_mesh_list_type), pointer :: model_mesh_list
  type(epic_mesh_list_type), pointer :: EPICGetPointerToMeshList

  select case(model_type)
  case (EPIC_MODEL_CLM)
    model_mesh_list => epic%clm_mesh_list
  case (EPIC_MODEL_PFLOTRAN)
    model_mesh_list => epic%pflotran_mesh_list
  case default
    write(*,*) 'EPICGetPointerToMeshList: Unknown model_type'
    stop
  end select

  EPICGetPointerToMeshList => model_mesh_list

end function EPICGetPointerToMeshList

! ************************************************************************** !

subroutine EPICAddData(model_type, comm, name, units, localsize, var_type, is_seq)

  implicit none

  PetscInt :: model_type
  PetscMPIInt :: comm
  character(len=EPIC_MAXWORDLENGTH)    :: name
  character(len=EPIC_MAXWORDLENGTH)    :: units
  PetscInt :: localsize
  PetscInt :: var_type
  PetscBool :: is_seq

  type(epic_data_list_type), pointer :: model_data_list

  model_data_list => EPICGetPointerToModelList(model_type)

  call EPICDataAddToList(model_data_list, comm, name, units, localsize, &
                         var_type, is_seq)

end subroutine EPICAddData

! ************************************************************************** !

subroutine EPICAddMesh(model_type, name, ncell, nlcell, ngcell, nidx_nlcell, &
                       nidx_ngcell)

  implicit none

  PetscInt :: model_type
  character(len=EPIC_MAXWORDLENGTH)    :: name
  PetscInt :: ncell
  PetscInt :: nlcell
  PetscInt :: ngcell
  PetscInt, pointer :: nidx_nlcell(:)
  PetscInt, pointer :: nidx_ngcell(:)

  type(epic_mesh_list_type), pointer :: model_mesh_list

  model_mesh_list => EPICGetPointerToMeshList(model_type)

  call EPICMeshAddToList(model_mesh_list, name, ncell, nlcell, ngcell, nidx_nlcell, &
                         nidx_ngcell)

end subroutine EPICAddMesh

! ************************************************************************** !

function EPICGetPointerToModeData(model_type, data_type)

  implicit none

  PetscInt :: model_type
  PetscInt :: data_type

  Vec, pointer :: EPICGetPointerToModeData

  type(epic_data_list_type), pointer :: model_data_list

  select case(model_type)
  case (EPIC_MODEL_CLM)
    model_data_list => epic%clm_data_list
  case (EPIC_MODEL_PFLOTRAN)
    model_data_list => epic%pflotran_data_list
  case default
    write(*,*) 'EPICGetPointerToModeData: Unknown model_type'
    stop
  end select

  EPICGetPointerToModeData => EPICDataGetPointer(model_data_list, data_type)

end function EPICGetPointerToModeData

end module EPICMod