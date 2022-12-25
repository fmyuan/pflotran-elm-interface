module Global_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  public GlobalSetup, &
         GlobalSetAuxVarScalar, &
         GlobalSetAuxVarVecLoc, &
         GlobalGetAuxVarVecLoc, &
         GlobalWeightAuxVars, &
         GlobalUpdateState, &
         GlobalUpdateAuxVars

contains

! ************************************************************************** !

subroutine GlobalSetup(realization)
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink

  PetscInt :: ghosted_id, iconn, sum_connection
  type(global_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: auxvars_bc(:)
  type(global_auxvar_type), pointer :: auxvars_ss(:)

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%Global => GlobalAuxCreate()

  ! allocate auxvar data structures for all grid cells
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  option%iflag = 1 ! allocate mass_balance array
#else
  option%iflag = 0 ! be sure not to allocate mass_balance array
#endif
  allocate(auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call GlobalAuxVarInit(auxvars(ghosted_id),option)
  enddo
  patch%aux%Global%auxvars => auxvars
  patch%aux%Global%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo

  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array
    allocate(auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call GlobalAuxVarInit(auxvars_bc(iconn),option)
    enddo
    patch%aux%Global%auxvars_bc => auxvars_bc
  endif
  patch%aux%Global%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    sum_connection = sum_connection + &
                     source_sink%connection_set%num_connections
    source_sink => source_sink%next
  enddo

  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array
    allocate(auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call GlobalAuxVarInit(auxvars_ss(iconn),option)
    enddo
    patch%aux%Global%auxvars_ss => auxvars_ss
  endif
  patch%aux%Global%num_aux_ss = sum_connection

  option%iflag = 0

end subroutine GlobalSetup

! ************************************************************************** !

subroutine GlobalSetAuxVarScalar(realization,value,ivar)
  !
  ! Sets values of auxvar data using a scalar value.
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/08
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_SATURATION, &
                               LIQUID_DENSITY, &
                               GAS_DENSITY, GAS_SATURATION, &
                               TEMPERATURE, LIQUID_DENSITY_MOL, &
                               GAS_DENSITY_MOL

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: value
  PetscInt :: ivar

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: auxvars_bc(:)

  PetscInt :: i
  PetscInt :: iphase
  PetscInt :: num_aux
  PetscInt :: num_aux_bc

  patch => realization%patch
  option => realization%option

  auxvars => patch%aux%Global%auxvars
  auxvars_bc => patch%aux%Global%auxvars_bc
  num_aux = patch%aux%Global%num_aux
  num_aux_bc = patch%aux%Global%num_aux_bc

  select case(ivar)
    case(LIQUID_PRESSURE)
      iphase = option%liquid_phase
      do i=1, num_aux
        auxvars(i)%pres(iphase) = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%pres(iphase) = value
      enddo
    case(TEMPERATURE)
      do i=1, num_aux
        auxvars(i)%temp = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%temp = value
      enddo
    case(LIQUID_DENSITY)
      iphase = option%liquid_phase
      do i=1, num_aux
        auxvars(i)%den_kg(iphase) = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%den_kg(iphase) = value
      enddo
    case(LIQUID_DENSITY_MOL)
      iphase = option%liquid_phase
      do i=1, num_aux
        auxvars(i)%den(iphase) = value/FMWH2O
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%den(iphase) = value/FMWH2O
      enddo
    case(LIQUID_SATURATION)
      iphase = option%liquid_phase
      do i=1, num_aux
        auxvars(i)%sat(iphase) = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%sat(iphase) = value
      enddo
    case(GAS_DENSITY)
      iphase = option%gas_phase
      do i=1, num_aux
        auxvars(i)%den_kg(iphase) = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%den_kg(iphase) = value
      enddo
    case(GAS_DENSITY_MOL)
      iphase = option%gas_phase
      do i=1, num_aux
        auxvars(i)%den(iphase) = value/FMWAIR
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%den(iphase) = value/FMWAIR
      enddo
    case(GAS_SATURATION)
      iphase = option%gas_phase
      do i=1, num_aux
        auxvars(i)%sat(iphase) = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%sat(iphase) = value
      enddo
    case default
      print *, 'Case(', ivar, ') not supported in GlobalSetAuxVarScalar'
      stop
  end select

end subroutine GlobalSetAuxVarScalar

! ************************************************************************** !

subroutine GlobalSetAuxVarVecLoc(realization,vec_loc,ivar,time_level)
  !
  ! Sets values of auxvar data using a vector.
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/08
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_SATURATION, &
                               LIQUID_DENSITY, GAS_PRESSURE, &
                               GAS_DENSITY, GAS_SATURATION, &
                               TEMPERATURE, SC_FUGA_COEFF, GAS_DENSITY_MOL, &
                               STATE, DARCY_VELOCITY

  implicit none

  class(realization_subsurface_type) :: realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: time_level

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: auxvars(:)

  PetscInt :: ghosted_id
  PetscInt :: iphase
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr

  patch => realization%patch
  grid => patch%grid
  option => realization%option

  auxvars => patch%aux%Global%auxvars

  call VecGetArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)

  select case(ivar)
    case(LIQUID_PRESSURE)
      iphase = option%liquid_phase
      select case(time_level)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres_store(iphase,TIME_T) = &
                vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres_store(iphase,TIME_TpDT) = &
                vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres(iphase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(GAS_PRESSURE)
      iphase = option%gas_phase
      select case(time_level)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres_store(iphase,TIME_T) = &
                vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres_store(iphase,TIME_TpDT) = &
                vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres(iphase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(TEMPERATURE)
      select case(time_level)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%temp_store(TIME_T) = vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%temp_store(TIME_TpDT) = vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%temp = vec_loc_p(ghosted_id)
          enddo
      end select
    case(LIQUID_DENSITY)
      iphase = option%liquid_phase
      select case(time_level)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_kg_store(iphase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_kg_store(iphase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_kg(iphase) = vec_loc_p(ghosted_id)
            auxvars(ghosted_id)%den(iphase) = vec_loc_p(ghosted_id)/FMWH2O
          enddo
      end select
    case(GAS_SATURATION)
      iphase = option%gas_phase
      select case(time_level)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat_store(iphase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat_store(iphase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat(iphase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(GAS_DENSITY)
      iphase = option%gas_phase
      select case(time_level)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_kg_store(iphase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_kg_store(iphase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_kg(iphase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(GAS_DENSITY_MOL)
      iphase = option%gas_phase
      select case(time_level)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_store(iphase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_store(iphase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den(iphase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(LIQUID_SATURATION)
      iphase = option%liquid_phase
      select case(time_level)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat_store(iphase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat_store(iphase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat(iphase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(SC_FUGA_COEFF)
      select case(time_level)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%fugacoeff_store(1,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%fugacoeff_store(1,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%fugacoeff(1) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(STATE)
      do ghosted_id=1, grid%ngmax
        auxvars(ghosted_id)%istate = int(vec_loc_p(ghosted_id)+1.d-10)
      enddo
    case(DARCY_VELOCITY)
      do ghosted_id=1, grid%ngmax
        patch%aux%Global%auxvars(ghosted_id)%darcy_vel(option%liquid_phase) = vec_loc_p(ghosted_id)
      enddo
      !end select
    case default
      print *, 'Case(', ivar, ') not supported in GlobalSetAuxVarVecLoc'
      stop
  end select

  call VecRestoreArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)

end subroutine GlobalSetAuxVarVecLoc

! ************************************************************************** !

subroutine GlobalGetAuxVarVecLoc(realization,vec_loc,ivar)
  !
  ! Sets values of auxvar data using a vector.
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/08
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Variables_module, only : STATE

  implicit none

  class(realization_subsurface_type) :: realization
  Vec :: vec_loc
  PetscInt :: ivar

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid

  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr

  patch => realization%patch
  grid => patch%grid
  option => realization%option

  call VecGetArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)

  select case(ivar)
    case(STATE)
      do ghosted_id=1, grid%ngmax
        vec_loc_p(ghosted_id) = &
          dble(patch%aux%Global%auxvars(ghosted_id)%istate)
      enddo
    case default
      option%io_buffer = 'Variable unrecognized in GlobalGetAuxVarVecLoc.'
      call PrintErrMsg(option)
  end select

  call VecRestoreArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)

end subroutine GlobalGetAuxVarVecLoc

! ************************************************************************** !

subroutine GlobalWeightAuxVars(realization,weight)
  !
  ! Updates the densities and saturations in auxiliary
  ! variables associated with reactive transport
  !
  ! Author: Glenn Hammond
  ! Date: 11/03/08
  !

  use Realization_Subsurface_class
  use Option_module
  use Material_module, only : MaterialWeightAuxVars

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: weight

  type(option_type), pointer :: option
  type(global_auxvar_type), pointer :: auxvars(:)
  PetscInt :: ghosted_id

  option => realization%option
  auxvars => realization%patch%aux%Global%auxvars

  ! interpolate variables based on weight
  do ghosted_id = 1, realization%patch%aux%Global%num_aux
    auxvars(ghosted_id)%sat(:) = &
      (weight*auxvars(ghosted_id)%sat_store(:,TIME_TpDT)+ &
       (1.d0-weight)*auxvars(ghosted_id)%sat_store(:,TIME_T))
  enddo

  select case(option%iflowmode)
    case(ZFLOW_MODE)
    case default
      do ghosted_id = 1, realization%patch%aux%Global%num_aux
        ! interpolate density and saturation based on weight
        auxvars(ghosted_id)%den_kg(:) = &
          (weight*auxvars(ghosted_id)%den_kg_store(:,TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%den_kg_store(:,TIME_T))
      enddo
  end select

  select case(option%iflowmode)
    case(G_MODE,H_MODE)
      do ghosted_id = 1, realization%patch%aux%Global%num_aux
        auxvars(ghosted_id)%pres(:) = &
          (weight*auxvars(ghosted_id)%pres_store(:,TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%pres_store(:,TIME_T))
        auxvars(ghosted_id)%temp = &
          (weight*auxvars(ghosted_id)%temp_store(TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%temp_store(TIME_T))
      enddo
    case(WF_MODE)
      do ghosted_id = 1, realization%patch%aux%Global%num_aux
        auxvars(ghosted_id)%pres(:) = &
          (weight*auxvars(ghosted_id)%pres_store(:,TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%pres_store(:,TIME_T))
      enddo
    case(MPH_MODE)
      ! need future implementation for ims_mode too
      do ghosted_id = 1, realization%patch%aux%Global%num_aux
        auxvars(ghosted_id)%pres(:) = &
          (weight*auxvars(ghosted_id)%pres_store(:,TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%pres_store(:,TIME_T))
        auxvars(ghosted_id)%temp = &
          (weight*auxvars(ghosted_id)%temp_store(TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%temp_store(TIME_T))
        auxvars(ghosted_id)%fugacoeff(:) = &
          (weight*auxvars(ghosted_id)%fugacoeff_store(:,TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%fugacoeff_store(:,TIME_T))
        if (weight<1D-12) auxvars(ghosted_id)%reaction_rate(:)=0D0
  !      auxvars(ghosted_id)%den(:) = &
  !        (weight*auxvars(ghosted_id)%den_store(:,TIME_TpDT)+ &
  !         (1.d0-weight)*auxvars(ghosted_id)%den_store(:,TIME_T))
      enddo
  end select

end subroutine GlobalWeightAuxVars

! ************************************************************************** !

subroutine GlobalUpdateState(realization)
  !
  ! Updates global aux var variables for use in
  ! reactive transport
  !
  ! Author: Glenn Hammond
  ! Date: 01/14/09
  !

  use Realization_Subsurface_class
  use Realization_Base_class, only : RealizationGetVariable
  use Communicator_Base_class
  use Variables_module, only : STATE

  class(realization_subsurface_type) :: realization

  call RealizationGetVariable(realization,realization%field%work,STATE, &
                              ZERO_INTEGER)
  call realization%comm1%GlobalToLocal(realization%field%work, &
                                       realization%field%work_loc)
  call GlobalSetAuxVarVecLoc(realization,realization%field%work_loc,STATE, &
                             ZERO_INTEGER)

end subroutine GlobalUpdateState

! ************************************************************************** !

subroutine GlobalUpdateAuxVars(realization,time_level,time)
  !
  ! Updates global aux var variables for use in
  ! reactive transport
  !
  ! Author: Glenn Hammond
  ! Date: 01/14/09
  !

  use Realization_Subsurface_class
  use Realization_Base_class, only : RealizationGetVariable
  use Field_module
  use Option_module
  use Communicator_Base_class
  use Material_module, only : MaterialUpdateAuxVars
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_SATURATION, &
                               LIQUID_DENSITY, GAS_PRESSURE, &
                               GAS_DENSITY, GAS_SATURATION, &
                               TEMPERATURE, SC_FUGA_COEFF, GAS_DENSITY_MOL, &
                               DARCY_VELOCITY
  use ZFlow_Aux_module, only : zflow_density_kg

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Discretization_module
  use Output_Common_module

  class(realization_subsurface_type) :: realization
  PetscReal :: time
  PetscInt  :: time_level

  PetscInt :: local_id, local_id_max
  Vec :: vec_x,vec_y,vec_z,global_vec
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:), vec_calc_ptr(:)
  Vec :: vec_calc
  PetscBool :: flag
  PetscErrorCode :: ierr
  type(discretization_type), pointer :: discretization

  type(field_type), pointer :: field
  type(option_type), pointer :: option

  option => realization%option
  field => realization%field

  select case(time_level)
    case(TIME_T)
      realization%patch%aux%Global%time_t = time
    case(TIME_TpDT)
      realization%patch%aux%Global%time_tpdt = time
  end select

  ! liquid density
  select case(option%iflowmode)
    case(ZFLOW_MODE)
      call GlobalSetAuxVarScalar(realization,zflow_density_kg,LIQUID_DENSITY)
    case default
      ! liquid density
      call GlobalUpdateSingleAuxVar(realization,LIQUID_DENSITY,time_level)
  end select

  ! liquid saturation
  call GlobalUpdateSingleAuxVar(realization,LIQUID_SATURATION,time_level)

  ! gas saturation
  flag = PETSC_FALSE
  select case(option%iflowmode)
    case(RICHARDS_MODE,TH_MODE,TH_TS_MODE,ZFLOW_MODE)
      if (option%transport%nphase > 1) flag = PETSC_TRUE
    case default
      flag = PETSC_TRUE
  end select
  if (flag) then
    call GlobalUpdateSingleAuxVar(realization,GAS_SATURATION,time_level)
  endif


  select case(option%iflowmode)
    case(MPH_MODE)
      ! Gas density
      call GlobalUpdateSingleAuxVar(realization,GAS_DENSITY,time_level)
      call GlobalUpdateSingleAuxVar(realization,GAS_DENSITY_MOL,time_level)

      ! liquid pressure
      call GlobalUpdateSingleAuxVar(realization,LIQUID_PRESSURE,time_level)

      ! gas pressure
      call GlobalUpdateSingleAuxVar(realization,GAS_PRESSURE,time_level)

      ! temperature
      call GlobalUpdateSingleAuxVar(realization,TEMPERATURE,time_level)

      ! fugacity coeff
      call GlobalUpdateSingleAuxVar(realization,SC_FUGA_COEFF,time_level)
    case(TH_MODE,TH_TS_MODE)
      ! pressure
      call GlobalUpdateSingleAuxVar(realization,LIQUID_DENSITY,time_level)

      ! temperature
      call GlobalUpdateSingleAuxVar(realization,TEMPERATURE,time_level)
    case(G_MODE,WF_MODE)
      ! pressure
      call GlobalUpdateSingleAuxVar(realization,LIQUID_DENSITY,time_level)
      if (option%iflowmode == G_MODE) then
        ! temperature
        call GlobalUpdateSingleAuxVar(realization,TEMPERATURE,time_level)
      endif
      ! Gas pressure
      call GlobalUpdateSingleAuxVar(realization,GAS_DENSITY,time_level)
      ! Gas density
      call GlobalUpdateSingleAuxVar(realization,GAS_DENSITY,time_level)
  end select

  ! darcy velocity (start)
  if (option%flow%store_darcy_vel) then

    !Create vectors of approapriate size
    discretization => realization%discretization
    call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                    option)
    call DiscretizationDuplicateVector(discretization,global_vec,vec_x)
    call DiscretizationDuplicateVector(discretization,global_vec,vec_y)
    call DiscretizationDuplicateVector(discretization,global_vec,vec_z)
    call DiscretizationDuplicateVector(discretization,global_vec,vec_calc)
    call OutputGetCellCenteredVelocities(realization,vec_x,vec_y,vec_z, &
                                         LIQUID_PHASE)

    ! open the vectors
    call VecGetArrayF90(vec_x,vec_x_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(vec_y,vec_y_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(vec_z,vec_z_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(vec_calc,vec_calc_ptr,ierr);CHKERRQ(ierr)

    ! the local size of the velocity vector
    ! local size = the number of cells calculated on the processor
    call VecGetLocalSize(vec_x,local_id_max,ierr);CHKERRQ(ierr)

    ! for each local(!) calculation calculate the velocity and store it
    do local_id=1, local_id_max
      vec_calc_ptr(local_id) = sqrt(vec_x_ptr(local_id)**2 &
                                   + vec_y_ptr(local_id)**2 &
                                   + vec_z_ptr(local_id)**2 ) &
                                   / realization%output_option%tconv
    enddo

    ! close the vectors
    call VecRestoreArrayF90(vec_calc,vec_calc_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(vec_x,vec_x_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(vec_y,vec_y_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(vec_z,vec_z_ptr,ierr);CHKERRQ(ierr)

    ! Set the auxvar variable for the DARCY_VELOCITY
    call realization%comm1%GlobalToLocal(vec_calc,field%work_loc)
    call GlobalSetAuxVarVecLoc(realization,field%work_loc,DARCY_VELOCITY, &
                              time_level)

    ! Destroy all vectors which were used for calculations
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
    call VecDestroy(vec_x,ierr);CHKERRQ(ierr)
    call VecDestroy(vec_y,ierr);CHKERRQ(ierr)
    call VecDestroy(vec_z,ierr);CHKERRQ(ierr)
    call VecDestroy(vec_calc,ierr);CHKERRQ(ierr)

  end if
  ! darcy velocity (end)

end subroutine GlobalUpdateAuxVars

! ************************************************************************** !

subroutine GlobalUpdateSingleAuxVar(realization,ivar,time_level)
  !
  ! Updates a single variable in global auxvar
  !
  ! Author: Glenn Hammond
  ! Date: 08/23/21
  !
  use Field_module
  use Realization_Subsurface_class
  use Realization_Base_class, only : RealizationGetVariable

  implicit none

  class(realization_subsurface_type) :: realization
  PetscInt :: ivar
  PetscInt :: time_level

  type(field_type), pointer :: field

  field => realization%patch%field
  call RealizationGetVariable(realization,field%work,ivar,ZERO_INTEGER)
  call realization%comm1%GlobalToLocal(field%work,field%work_loc)
  call GlobalSetAuxVarVecLoc(realization,field%work_loc,ivar,time_level)

end subroutine GlobalUpdateSingleAuxVar

end module Global_module
