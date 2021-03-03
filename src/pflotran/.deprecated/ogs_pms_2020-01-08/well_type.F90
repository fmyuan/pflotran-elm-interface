module Well_Type_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  public

! Well types (producer is always multi-phase, different types of injectors)

  PetscInt, parameter, public :: PROD_WELL_TYPE      = 1
  PetscInt, parameter, public :: WAT_INJ_WELL_TYPE   = 2
  PetscInt, parameter, public :: OIL_INJ_WELL_TYPE   = 3
  PetscInt, parameter, public :: GAS_INJ_WELL_TYPE   = 4
  PetscInt, parameter, public :: SLV_INJ_WELL_TYPE   = 5

  PetscBool :: wd_isothermal=PETSC_TRUE

! Well target types (BHP=well pressure, SV=surface volume rate, M= mass rate)

  PetscInt, parameter, public :: W_BHP_LIMIT         =  1
  PetscInt, parameter, public :: W_TARG_OSV          =  2
  PetscInt, parameter, public :: W_TARG_GSV          =  3
  PetscInt, parameter, public :: W_TARG_WSV          =  4
  PetscInt, parameter, public :: W_TARG_SSV          =  5
  PetscInt, parameter, public :: W_TARG_LSV          =  6
  PetscInt, parameter, public :: W_TARG_OM           =  7
  PetscInt, parameter, public :: W_TARG_GM           =  8
  PetscInt, parameter, public :: W_TARG_WM           =  9
  PetscInt, parameter, public :: W_TARG_SM           = 10
  PetscInt, parameter, public :: W_TARG_RV           = 11
! Next value should be equal to number of target types
  PetscInt, parameter, public :: N_WELL_TT           = 11

! Well target type flag values: pressure, volume or mass

  PetscInt, parameter :: TT_P                        = 1
  PetscInt, parameter :: TT_V                        = 2
  PetscInt, parameter :: TT_M                        = 3

! Well completion factor types

  PetscInt, parameter, public :: WELL_FACTOR_CONST    = 1
  PetscInt, parameter, public :: WELL_FACTOR_PEACEMAN = 2

! Well open/shut status flags

  PetscInt, parameter, public :: W_STATUS_OPEN = 1
  PetscInt, parameter, public :: W_STATUS_SHUT = 2

! Well completion open/closed flags

  PetscInt, parameter, public :: CONN_STATUS_OPEN  = 1
  PetscInt, parameter, public :: CONN_STATUS_CLOSE = 2

end module Well_Type_class
