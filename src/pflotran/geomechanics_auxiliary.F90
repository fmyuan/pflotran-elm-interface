module Geomechanics_Auxiliary_module

#include "petsc/finclude/petscvec.h"
  use petscvec

  use Geomechanics_Global_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: geomech_auxiliary_type
    type(geomech_global_type), pointer :: GeomechGlobal
    type(geomech_parameter_type), pointer :: GeomechParam
  end type geomech_auxiliary_type

  type, public :: geomech_parameter_type
    PetscReal, pointer :: youngs_modulus(:)
    PetscReal, pointer :: poissons_ratio(:)
    PetscReal, pointer :: biot_coeff(:)
    PetscReal, pointer :: thermal_exp_coeff(:)
    PetscReal, pointer :: density(:)

    PetscBool :: youngs_modulus_spatially_varying
    PetscBool :: poissons_ratio_spatially_varying
    PetscBool :: density_spatially_varying
    PetscBool :: biot_coeff_spatially_varying
    PetscBool :: thermal_exp_coeff_spatially_varying

    Vec :: youngs_modulus_vec
    Vec :: poissons_ratio_vec
    Vec :: density_vec
    Vec :: biot_coeff_vec
    Vec :: thermal_exp_coeff_vec
  end type geomech_parameter_type

  public :: GeomechAuxInit, &
            GeomechAuxDestroy

contains

! ************************************************************************** !

subroutine GeomechAuxInit(geomech_aux)
  !
  ! Nullifies pointers in geomech auxiliary type
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13
  !

  implicit none

  type(geomech_auxiliary_type) :: geomech_aux

  nullify(geomech_aux%GeomechGlobal)
  allocate(geomech_aux%GeomechParam)
  nullify(geomech_aux%GeomechParam%youngs_modulus)
  nullify(geomech_aux%GeomechParam%poissons_ratio)
  nullify(geomech_aux%GeomechParam%biot_coeff)
  nullify(geomech_aux%GeomechParam%thermal_exp_coeff)
  nullify(geomech_aux%GeomechParam%density)

end subroutine GeomechAuxInit

! ************************************************************************** !

subroutine GeomechAuxDestroy(geomech_aux)
  !
  ! Strips a geomech auxiliary type
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13
  !

  implicit none

  type(geomech_auxiliary_type) :: geomech_aux

  call GeomechGlobalAuxDestroy(geomech_aux%GeomechGlobal)

  nullify(geomech_aux%GeomechGlobal)

  if (associated(geomech_aux%GeomechParam)) then
    if (associated(geomech_aux%GeomechParam%youngs_modulus)) &
      deallocate(geomech_aux%GeomechParam%youngs_modulus)
    nullify(geomech_aux%GeomechParam%youngs_modulus)
    if (associated(geomech_aux%GeomechParam%poissons_ratio)) &
      deallocate(geomech_aux%GeomechParam%poissons_ratio)
    nullify(geomech_aux%GeomechParam%poissons_ratio)
    if (associated(geomech_aux%GeomechParam%biot_coeff)) &
      deallocate(geomech_aux%GeomechParam%biot_coeff)
    nullify(geomech_aux%GeomechParam%biot_coeff)
    if (associated(geomech_aux%GeomechParam%thermal_exp_coeff)) &
      deallocate(geomech_aux%GeomechParam%thermal_exp_coeff)
    nullify(geomech_aux%GeomechParam%thermal_exp_coeff)
    if (associated(geomech_aux%GeomechParam%density)) &
      deallocate(geomech_aux%GeomechParam%density)
    nullify(geomech_aux%GeomechParam%density)
  endif

  if (associated(geomech_aux%GeomechParam)) then
    deallocate(geomech_aux%GeomechParam)
  endif
  nullify(geomech_aux%GeomechParam)

end subroutine GeomechAuxDestroy

end module Geomechanics_Auxiliary_module
