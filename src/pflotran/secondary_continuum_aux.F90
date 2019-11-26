! added by S. Karra 07/11/12

module Secondary_Continuum_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: slab_type
    PetscReal :: length                       ! input - length of slab
    PetscReal :: area                         ! input - surface area
  end type slab_type
  
  type, public :: nested_cube_type
    PetscReal :: matrix_block_size            ! input - side of cube
    PetscReal :: fracture_spacing             ! input - fracture spacing
  end type nested_cube_type
  
  type, public :: nested_sphere_type
    PetscReal :: radius                       ! input - radius of sphere
  end type nested_sphere_type
  
  type, public :: sec_continuum_type
    PetscInt :: itype                         ! input - type of sec. continuum (slab, nested_cube, nested_sphere,....) 
    type(slab_type) :: slab
    type(nested_cube_type) :: nested_cube
    type(nested_sphere_type) :: nested_sphere 
    PetscReal, pointer :: distance(:)         ! This is the array of positions of cells centers from the center
  end type sec_continuum_type                 ! farthest is the cell center closest to interface between prim. and sec. continuua

  type, public :: sec_heat_type  
    PetscInt :: ncells                         ! number of secondary grid cells
    PetscReal :: aperture                      ! fracture aperture
    PetscReal :: epsilon                       ! vol. frac. of primary continuum
    type(sec_continuum_type) :: sec_continuum
    PetscReal, pointer :: sec_temp(:)          ! array of temp. at secondary grid cells
    PetscReal, pointer :: area(:)              ! surface area
    PetscReal, pointer :: vol(:)               ! volume     face      node       face
    PetscReal, pointer :: dm_plus(:)           ! see fig.    |----------o----------|
    PetscReal, pointer :: dm_minus(:)          ! see fig.      <dm_minus> <dm_plus>
    PetscReal :: interfacial_area              ! interfacial area between prim. and sec. per unit volume of prim.+sec.
    PetscBool :: log_spacing                   ! flag to check if log spacing is set
    PetscReal :: outer_spacing                 ! value of the outer most grid cell spacing
  end type sec_heat_type  
 
  type, public :: sc_heat_type
    type(sec_heat_type), pointer :: sec_heat_vars(:)
  end type sc_heat_type

  public :: SecondaryAuxHeatCreate, SecondaryAuxHeatDestroy
            
contains

! ************************************************************************** !

function SecondaryAuxHeatCreate(option)
  ! 
  ! Allocate and initialize secondary continuum heat
  ! auxiliary object
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 01/10/13
  ! 

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(sc_heat_type), pointer :: SecondaryAuxHeatCreate
  
  type(sc_heat_type), pointer :: aux

  allocate(aux) 
  nullify(aux%sec_heat_vars)
  
  SecondaryAuxHeatCreate => aux
  
end function SecondaryAuxHeatCreate  

! ************************************************************************** !

subroutine SecondaryAuxHeatDestroy(aux)
  ! 
  ! Deallocates a secondary continuum heat
  ! auxiliary object
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 01/10/13
  ! 

  implicit none

  type(sc_heat_type), pointer :: aux
   
  if (.not.associated(aux)) return
  
  deallocate(aux)
  nullify(aux)  

end subroutine SecondaryAuxHeatDestroy

! ************************************************************************** !


end module Secondary_Continuum_Aux_module
            
