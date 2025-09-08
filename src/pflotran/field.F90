module Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: field_type

    !get material id
    ! 1 degree of freedom
    Vec :: porosity0
    Vec :: porosity_base_store
    Vec :: porosity_t
    Vec :: porosity_tpdt
    Vec :: porosity_geomech_store
    Vec :: tortuosity0

    Vec :: perm0_xx, perm0_yy, perm0_zz
    Vec :: perm0_xz, perm0_xy, perm0_yz

    Vec :: work, work_loc

    Vec :: volume0
    Vec :: compressibility0

    !TODO(geh): move these Vecs into their respective pms
    ! residual vectors
    Vec :: flow_r
    Vec :: tran_r

    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: flow_xx, flow_xx_loc, flow_dxx, flow_yy, flow_accum, flow_accum2
    Vec :: tran_xx, tran_xx_loc, tran_dxx, tran_yy, tran_accum
    Vec :: flow_xxdot, flow_xxdot_loc
    Vec :: flow_rhs

    ! vectors for advanced nonlinear solvers other than Newton - Heeho
    Vec :: flow_scaled_xx, flow_work_loc

    Vec :: tran_log_xx, tran_work_loc

    ! mass transfer
    Vec :: flow_mass_transfer
    Vec :: tran_mass_transfer

    Vec :: flow_ts_mass_balance, flow_total_mass_balance
    Vec :: tran_ts_mass_balance, tran_total_mass_balance

    ! vector that holds the second layer of ghost cells for tvd
    Vec :: tvd_ghosts

    ! vectors to save temporally average quantities
    Vec, pointer :: avg_vars_vec(:)
    PetscInt :: nvars

    ! vectors to save temporally average flowrates
    Vec :: flowrate_inst
    Vec :: flowrate_aveg

    ! vectors to save velocity at face
    Vec :: vx_face_inst
    Vec :: vy_face_inst
    Vec :: vz_face_inst

    Vec, pointer :: max_change_vecs(:)

  end type field_type

  public :: FieldCreate, &
            FieldDestroy

contains

! ************************************************************************** !

function FieldCreate()
  !
  ! Allocates and initializes a new Field object
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

  implicit none

  type(field_type), pointer :: FieldCreate

  type(field_type), pointer :: field

  allocate(field)

  ! nullify PetscVecs
  PetscObjectNullify(field%porosity0)
  PetscObjectNullify(field%porosity_base_store)
  PetscObjectNullify(field%porosity_geomech_store)
  PetscObjectNullify(field%porosity_t)
  PetscObjectNullify(field%porosity_tpdt)
  PetscObjectNullify(field%tortuosity0)

  PetscObjectNullify(field%perm0_xx)
  PetscObjectNullify(field%perm0_yy)
  PetscObjectNullify(field%perm0_zz)
  PetscObjectNullify(field%perm0_xy)
  PetscObjectNullify(field%perm0_xz)
  PetscObjectNullify(field%perm0_yz)

  PetscObjectNullify(field%work)
  PetscObjectNullify(field%work_loc)

  PetscObjectNullify(field%volume0)
  PetscObjectNullify(field%compressibility0)

  PetscObjectNullify(field%flow_r)
  PetscObjectNullify(field%flow_xx)
  PetscObjectNullify(field%flow_xx_loc)
  PetscObjectNullify(field%flow_scaled_xx)
  PetscObjectNullify(field%flow_work_loc)
  PetscObjectNullify(field%flow_dxx)
  PetscObjectNullify(field%flow_yy)
  PetscObjectNullify(field%flow_accum)
  PetscObjectNullify(field%flow_accum2)
  PetscObjectNullify(field%flow_xxdot)
  PetscObjectNullify(field%flow_xxdot_loc)
  PetscObjectNullify(field%flow_rhs)

  PetscObjectNullify(field%tran_r)
  PetscObjectNullify(field%tran_log_xx)
  PetscObjectNullify(field%tran_xx)
  PetscObjectNullify(field%tran_xx_loc)
  PetscObjectNullify(field%tran_dxx)
  PetscObjectNullify(field%tran_yy)
  PetscObjectNullify(field%tran_accum)
  PetscObjectNullify(field%tran_work_loc)

  PetscObjectNullify(field%tvd_ghosts)

  PetscObjectNullify(field%flow_mass_transfer)
  PetscObjectNullify(field%tran_mass_transfer)

  PetscObjectNullify(field%flow_ts_mass_balance)
  PetscObjectNullify(field%flow_total_mass_balance)
  PetscObjectNullify(field%tran_ts_mass_balance)
  PetscObjectNullify(field%tran_total_mass_balance)

  nullify(field%avg_vars_vec)
  field%nvars = 0

  PetscObjectNullify(field%flowrate_inst)
  PetscObjectNullify(field%flowrate_aveg)

  PetscObjectNullify(field%vx_face_inst)
  PetscObjectNullify(field%vy_face_inst)
  PetscObjectNullify(field%vz_face_inst)

  nullify(field%max_change_vecs)

  FieldCreate => field

end function FieldCreate

! ************************************************************************** !

subroutine FieldDestroy(field)
  !
  ! Deallocates a field object
  !
  ! Author: Glenn Hammond
  ! Date: 11/15/07
  !

  implicit none

  type(field_type), pointer :: field

  PetscErrorCode :: ierr
  PetscInt :: ivar
  PetscInt :: num_vecs

  if (.not.associated(field)) return

  ! Destroy PetscVecs
  if (.not.PetscObjectIsNull(field%porosity0)) then
    call VecDestroy(field%porosity0,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%porosity_base_store)) then
    call VecDestroy(field%porosity_base_store,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%porosity_geomech_store)) then
    call VecDestroy(field%porosity_geomech_store,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%porosity_t)) then
    call VecDestroy(field%porosity_t,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%porosity_tpdt)) then
    call VecDestroy(field%porosity_tpdt,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tortuosity0)) then
    call VecDestroy(field%tortuosity0,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%perm0_xx)) then
    call VecDestroy(field%perm0_xx,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%perm0_yy)) then
    call VecDestroy(field%perm0_yy,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%perm0_zz)) then
    call VecDestroy(field%perm0_zz,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%perm0_xy)) then
    call VecDestroy(field%perm0_xy,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%perm0_xz)) then
    call VecDestroy(field%perm0_xz,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%perm0_yz)) then
    call VecDestroy(field%perm0_yz,ierr);CHKERRQ(ierr)
  endif

  if (.not.PetscObjectIsNull(field%work)) then
    call VecDestroy(field%work,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%work_loc)) then
    call VecDestroy(field%work_loc,ierr);CHKERRQ(ierr)
  endif

  if (.not.PetscObjectIsNull(field%volume0)) then
    call VecDestroy(field%volume0,ierr);CHKERRQ(ierr)
  endif

  if (.not.PetscObjectIsNull(field%compressibility0)) then
    call VecDestroy(field%compressibility0,ierr);CHKERRQ(ierr)
  endif

  if (.not.PetscObjectIsNull(field%flow_r)) then
    call VecDestroy(field%flow_r,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_scaled_xx)) then
    call VecDestroy(field%flow_scaled_xx,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_work_loc)) then
    call VecDestroy(field%flow_work_loc,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_xx)) then
    call VecDestroy(field%flow_xx,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_xx_loc)) then
    call VecDestroy(field%flow_xx_loc,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_dxx)) then
    call VecDestroy(field%flow_dxx,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_yy)) then
    call VecDestroy(field%flow_yy,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_accum)) then
    call VecDestroy(field%flow_accum,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_accum2)) then
    call VecDestroy(field%flow_accum2,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_xxdot)) then
    call VecDestroy(field%flow_xxdot,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_xxdot_loc)) then
    call VecDestroy(field%flow_xxdot_loc,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_rhs)) then
    call VecDestroy(field%flow_rhs,ierr);CHKERRQ(ierr)
  endif

  if (.not.PetscObjectIsNull(field%tran_r)) then
    call VecDestroy(field%tran_r,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tran_log_xx)) then
    call VecDestroy(field%tran_log_xx,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tran_xx)) then
    call VecDestroy(field%tran_xx,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tran_xx_loc)) then
    call VecDestroy(field%tran_xx_loc,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tran_dxx)) then
    call VecDestroy(field%tran_dxx,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tran_yy)) then
    call VecDestroy(field%tran_yy,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tran_accum)) then
    call VecDestroy(field%tran_accum,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tran_work_loc)) then
    call VecDestroy(field%tran_work_loc,ierr);CHKERRQ(ierr)
  endif

  if (.not.PetscObjectIsNull(field%flow_mass_transfer)) then
    call VecDestroy(field%flow_mass_transfer,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tran_mass_transfer)) then
    call VecDestroy(field%tran_mass_transfer,ierr);CHKERRQ(ierr)
  endif

  if (.not.PetscObjectIsNull(field%flow_ts_mass_balance)) then
    call VecDestroy(field%flow_ts_mass_balance,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flow_total_mass_balance)) then
    call VecDestroy(field%flow_total_mass_balance,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tran_ts_mass_balance)) then
    call VecDestroy(field%tran_ts_mass_balance,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%tran_total_mass_balance)) then
    call VecDestroy(field%tran_total_mass_balance,ierr);CHKERRQ(ierr)
  endif

  if (.not.PetscObjectIsNull(field%tvd_ghosts)) then
    call VecDestroy(field%tvd_ghosts,ierr);CHKERRQ(ierr)
  endif

  do ivar = 1,field%nvars
    call VecDestroy(field%avg_vars_vec(ivar),ierr);CHKERRQ(ierr)
  enddo

  if (.not.PetscObjectIsNull(field%flowrate_inst)) then
    call VecDestroy(field%flowrate_inst,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%flowrate_aveg)) then
    call VecDestroy(field%flowrate_aveg,ierr);CHKERRQ(ierr)
  endif

  if (.not.PetscObjectIsNull(field%vx_face_inst)) then
    call VecDestroy(field%vx_face_inst,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%vy_face_inst)) then
    call VecDestroy(field%vy_face_inst,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(field%vz_face_inst)) then
    call VecDestroy(field%vz_face_inst,ierr);CHKERRQ(ierr)
  endif

  if (associated(field%max_change_vecs)) then
    !geh: kludge as the compiler returns i4 in 64-bit
    num_vecs = size(field%max_change_vecs)
    call VecDestroyVecs(num_vecs,field%max_change_vecs,ierr);CHKERRQ(ierr)
  endif

  deallocate(field)
  nullify(field)

end subroutine FieldDestroy

end module Field_module
