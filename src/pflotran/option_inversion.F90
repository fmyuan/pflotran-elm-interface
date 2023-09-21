module Option_Inversion_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Communicator_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: inversion_option_type
    type(comm_type), pointer :: invcomm
    type(comm_type), pointer :: forcomm
    type(comm_type), pointer :: forcomm_i
    PetscInt :: num_process_groups
    PetscBool :: use_perturbation
    PetscBool :: perturbation_run
    PetscBool :: coupled_flow_ert
    PetscBool :: record_measurements
    PetscBool :: calculate_ert
    PetscBool :: calculate_ert_jacobian
    ! parameter flags
    PetscBool :: invert_for_elec_cond
    PetscBool :: invert_for_permeability
    PetscBool :: invert_for_porosity
    PetscBool :: invert_for_vg_alpha
    PetscBool :: invert_for_vg_m
    PetscBool :: invert_for_vg_sr
    PetscBool :: invert_for_arch_cement_exp
    PetscBool :: invert_for_arch_sat_exp
    PetscBool :: invert_for_arch_tort_const
    character(len=MAXWORDLENGTH) :: iteration_prefix
    character(len=MAXSTRINGLENGTH) :: restart_filename
  end type inversion_option_type

  public :: OptionInversionCreate, &
            OptionInversionInit, &
            OptionInversionDestroy

contains

! ************************************************************************** !

function OptionInversionCreate()
  !
  ! Allocates and initializes a new OptionInversion object
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/21
  !

  implicit none

  type(inversion_option_type), pointer :: OptionInversionCreate

  type(inversion_option_type), pointer :: option

  allocate(option)

  call OptionInversionInit(option)
  OptionInversionCreate => option

end function OptionInversionCreate

! ************************************************************************** !

subroutine OptionInversionInit(option)
  !
  ! Initializes all option variables
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/21
  !

  implicit none

  type(inversion_option_type) :: option

  nullify(option%invcomm)
  nullify(option%forcomm)
  nullify(option%forcomm_i)
  option%num_process_groups = 1

  option%use_perturbation = PETSC_FALSE
  option%perturbation_run = PETSC_FALSE
  option%coupled_flow_ert = PETSC_FALSE
  option%record_measurements = PETSC_TRUE
  option%calculate_ert = PETSC_FALSE
  option%calculate_ert_jacobian = PETSC_FALSE
  option%iteration_prefix = ''
  option%restart_filename = ''

  option%invert_for_elec_cond = PETSC_FALSE
  option%invert_for_permeability = PETSC_FALSE
  option%invert_for_porosity = PETSC_FALSE
  option%invert_for_vg_alpha = PETSC_FALSE
  option%invert_for_vg_m = PETSC_FALSE
  option%invert_for_vg_sr = PETSC_FALSE
  option%invert_for_arch_cement_exp = PETSC_FALSE
  option%invert_for_arch_sat_exp = PETSC_FALSE
  option%invert_for_arch_tort_const = PETSC_FALSE

end subroutine OptionInversionInit

! ************************************************************************** !

subroutine OptionInversionDestroy(option)
  !
  ! Deallocates an option
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/21
  !

  implicit none

  type(inversion_option_type), pointer :: option

  if (.not.associated(option)) return

  call CommDestroy(option%invcomm)
  call CommDestroy(option%forcomm)
  call CommDestroy(option%forcomm_i)

  deallocate(option)
  nullify(option)

end subroutine OptionInversionDestroy

end module Option_Inversion_module
