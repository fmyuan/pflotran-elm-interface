module ERT_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module
  
  implicit none
  
  private 

  type, public :: ert_auxvar_type
  
    PetscReal :: bulk_conductivity

  end type ert_auxvar_type
  
  type, public :: ert_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(ert_auxvar_type), pointer :: auxvars(:)
    type(ert_auxvar_type), pointer :: auxvars_bc(:)
    type(ert_auxvar_type), pointer :: auxvars_ss(:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type ert_type

  public :: ERTAuxCreate, ERTAuxDestroy, &
            ERTAuxVarCompute, ERTAuxVarInit, &
            ERTAuxVarCopy

contains

! ************************************************************************** !

function ERTAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  ! 

  use Option_module

  implicit none
  
  type(ert_type), pointer :: ERTAuxCreate
  
  type(ert_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%matrix_zeroing)

  ERTAuxCreate => aux
  
end function ERTAuxCreate

! ************************************************************************** !

subroutine ERTAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  ! 

  use Option_module

  implicit none
  
  type(ert_auxvar_type) :: auxvar
  type(option_type) :: option
  
  auxvar%bulk_conductivity = 0.d0

end subroutine ERTAuxVarInit

! ************************************************************************** !

subroutine ERTAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  ! 

  use Option_module

  implicit none
  
  type(ert_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%bulk_conductivity = auxvar%bulk_conductivity

end subroutine ERTAuxVarCopy

! ************************************************************************** !

subroutine ERTAuxVarCompute(x,ert_auxvar,global_auxvar,rt_auxvar, &
                            material_auxvar,natural_id,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  ! 

  use Option_module
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  type(ert_auxvar_type) :: ert_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reactive_transport_auxvar_type) :: rt_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  
  ! calculate ert_auxvar%bulk_conductivity = f(global,rt,material-auxars)
  
end subroutine ERTAuxVarCompute

! ************************************************************************** !

subroutine AuxVarDestroy(auxvar)
  ! 
  ! Deallocates a ert auxiliary object
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  ! 

  implicit none

  type(ert_auxvar_type) :: auxvar
  
end subroutine AuxVarDestroy

! ************************************************************************** !

subroutine ERTAuxDestroy(aux)
  ! 
  ! Deallocates a ert auxiliary object
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/11/21
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none

  type(ert_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call AuxVarDestroy(aux%auxvars(iaux))
    enddo  
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call AuxVarDestroy(aux%auxvars_bc(iaux))
    enddo  
    deallocate(aux%auxvars_bc)
  endif
  nullify(aux%auxvars_bc)
  if (associated(aux%auxvars_ss)) then
    do iaux = 1, aux%num_aux_ss
      call AuxVarDestroy(aux%auxvars_ss(iaux))
    enddo  
    deallocate(aux%auxvars_ss)
  endif
  nullify(aux%auxvars_ss)
  
  call MatrixZeroingDestroy(aux%matrix_zeroing)

  deallocate(aux)
  nullify(aux)
    
end subroutine ERTAuxDestroy

end module ERT_Aux_module
