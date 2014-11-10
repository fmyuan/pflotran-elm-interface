module EPICMesh

  use EPICConstants

  implicit none

#include "finclude/petscsys.h"

  private

  type, public :: epic_mesh_type

     character(len=EPIC_MAXWORDLENGTH)  :: name              ! Data name
     PetscInt                           :: ncell             ! Total number of cells in the mesh
     PetscInt                           :: nlcell            ! Number of non-ghosted ("local") cells
     PetscInt                           :: ngcell            ! Number of ghosted + non-ghosted cells

     PetscInt,pointer                   :: nidx_nlcell(:)    ! Natural ID of non-ghosted cells
     PetscInt,pointer                   :: nidx_ngcell(:)    ! Natural ID of ghosted + non-ghosted cells

    type(epic_mesh_type), pointer :: next

  end type epic_mesh_type

  type, public :: epic_mesh_list_type
    PetscInt                         :: nmesh
    type(epic_mesh_type), pointer    :: first
    type(epic_mesh_type), pointer    :: last
  end type epic_mesh_list_type

  interface EPICMeshCreate
    module procedure EPICMeshCreate1
    module procedure EPICMeshCreate2
  end interface EPICMeshCreate

  interface EPICMeshAddToList
    module procedure EPICMeshAddToList1
    module procedure EPICMeshAddToList2
  end interface EPICMeshAddToList

  public :: EPICMeshCreate
  public :: EPICMeshDestroy
  public :: EPICMeshListCreate
  public :: EPICMeshAddToList
!  public :: EPICMeshGetPointer
  public :: EPICMeshListDestroy

contains

! ************************************************************************** !

function EPICMeshCreate1()

  implicit none
  
  type(epic_mesh_type),pointer :: EPICMeshCreate1

  type(epic_mesh_type),pointer :: mesh

  allocate(mesh)
  mesh%name = ''

  mesh%ncell = 0
  mesh%nlcell = 0
  mesh%ngcell = 0
  nullify(mesh%nidx_nlcell)
  nullify(mesh%nidx_ngcell)
  nullify(mesh%next)

  EPICMeshCreate1 => mesh

end function EPICMeshCreate1

! ************************************************************************** !

function EPICMeshCreate2(name, ncell, nlcell, ngcell, nidx_nlcell, nidx_ngcell)

  implicit none
  
  character(len=EPIC_MAXWORDLENGTH)    :: name
  PetscInt :: ncell
  PetscInt :: nlcell
  PetscInt :: ngcell
  PetscInt, pointer :: nidx_nlcell(:)
  PetscInt, pointer :: nidx_ngcell(:)

  type(epic_mesh_type),pointer :: EPICMeshCreate2

  type(epic_mesh_type),pointer :: mesh

  allocate(mesh)

  mesh%name = trim(adjustl(name))
  mesh%ncell = ncell
  mesh%nlcell = nlcell
  mesh%ngcell = ngcell

  allocate(mesh%nidx_nlcell(nlcell))
  allocate(mesh%nidx_ngcell(ngcell))
  mesh%nidx_nlcell(:) = nidx_nlcell(:)
  mesh%nidx_ngcell(:) = nidx_ngcell(:)

  nullify(mesh%next)

  EPICMeshCreate2 => mesh

end function EPICMeshCreate2

! ************************************************************************** !

function EPICMeshListCreate()

  implicit none
  
  type(epic_mesh_list_type), pointer :: EPICMeshListCreate
  
  type(epic_mesh_list_type), pointer :: mesh_list
  
  allocate(mesh_list)
  nullify(mesh_list%first)
  nullify(mesh_list%last)
  mesh_list%nmesh = 0
  
  EPICMeshListCreate => mesh_list
  
end function EPICMeshListCreate

! ************************************************************************** !

subroutine EPICMeshAddToList1(list, mesh)

  implicit none
  
  type(epic_mesh_list_type) :: list
  type(epic_mesh_type), pointer :: mesh
  
  if (.not. associated(list%first)) then
    list%first => mesh
  else
    list%last%next => mesh
  endif
  list%last => mesh
  
  list%nmesh = list%nmesh + 1
  
end subroutine EPICMeshAddToList1

! ************************************************************************** !

subroutine EPICMeshAddToList2(list, name, ncell, nlcell, ngcell, nidx_nlcell, &
                              nidx_ngcell)

  implicit none
  
  type(epic_mesh_list_type) :: list
  character(len=EPIC_MAXWORDLENGTH)    :: name
  PetscInt :: ncell
  PetscInt :: nlcell
  PetscInt :: ngcell
  PetscInt, pointer :: nidx_nlcell(:)
  PetscInt, pointer :: nidx_ngcell(:)

  type(epic_mesh_type), pointer :: mesh

  mesh => EPICMeshCreate(name, ncell, nlcell, ngcell, nidx_nlcell, nidx_ngcell)
  call EPICMeshAddToList1(list, mesh)

end subroutine EPICMeshAddToList2

! ************************************************************************** !

subroutine EPICMeshListDestroy(list)

  implicit none
  
  type(epic_mesh_list_type), pointer :: list
  
  nullify(list%last)
  call OutputVariableDestroy(list%first)
  
  deallocate(list)
  nullify(list)
  
end subroutine EPICMeshListDestroy

! ************************************************************************** !

recursive subroutine EPICMeshDestroy(data)

  implicit none
  
  type(epic_mesh_type), pointer :: data
  
  if (.not.associated(data)) return
  
  call EPICMeshDestroy(data%next)
  
  deallocate(data)
  nullify(data)
  
end subroutine EPICMeshDestroy

end module EPICMesh