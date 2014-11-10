module EPICConstants

  implicit none

  private

#include "finclude/petscsys.h"

  ! Indentifiers for various models (eg. CLM, PFLOTRAN, etc.,)
  PetscInt, parameter, public :: EPIC_MODEL_NULL      = 0
  PetscInt, parameter, public :: EPIC_MODEL_CLM       = 1
  PetscInt, parameter, public :: EPIC_MODEL_PFLOTRAN  = 2

  ! Indentifiers for various types of data
  PetscInt, parameter, public :: EPIC_VAR_NULL                     = 0

  ! 3D subsurface parameters
  PetscInt, parameter, public :: EPIC_VAR_SOIL_HKSAT_X             = 1
  PetscInt, parameter, public :: EPIC_VAR_SOIL_HKSAT_Y             = 2
  PetscInt, parameter, public :: EPIC_VAR_SOIL_HKSAT_Z             = 3
  PetscInt, parameter, public :: EPIC_VAR_SOIL_SUCSAT              = 4
  PetscInt, parameter, public :: EPIC_VAR_SOIL_WATSAT              = 5
  PetscInt, parameter, public :: EPIC_VAR_SOIL_BSW                 = 6

  PetscInt, parameter, public :: EPIC_VAR_MESH_FACE_AREA           = 101

  PetscInt, parameter, public :: EPIC_VAR_SUBSURF_PRESS            = 201
  PetscInt, parameter, public :: EPIC_VAR_SUBSURF_TMP              = 202
  PetscInt, parameter, public :: EPIC_VAR_SUBSURF_LIQ_SAT          = 203
  PetscInt, parameter, public :: EPIC_VAR_SUBSURF_ICE_SAT          = 204
  PetscInt, parameter, public :: EPIC_VAR_SUBSURF_THERMAL_COND     = 205

  PetscInt, parameter, public :: EPIC_VAR_SURF_HEAD                = 301
  PetscInt, parameter, public :: EPIC_VAR_SURF_TMP                 = 302

  PetscInt, parameter, public :: EPIC_VAR_RAIN_FLUX                = 401
  PetscInt, parameter, public :: EPIC_VAR_RAIN_TMP                 = 402
  PetscInt, parameter, public :: EPIC_VAR_ET_FLUX                  = 403
  PetscInt, parameter, public :: EPIC_VAR_GFLUX                    = 404

  PetscInt, parameter, public :: EPIC_MAXWORDLENGTH = 512

end module EPICConstants