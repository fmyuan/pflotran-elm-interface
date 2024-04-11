module Global_module

#include "petsc/finclude/petscsys.h"
  use petscsys  
  use Global_Aux_module
  
  use PFLOTRAN_Constants_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_SATURATION, LIQUID_DENSITY, &
                               GAS_PRESSURE, GAS_SATURATION, GAS_DENSITY, &
                               TEMPERATURE

  implicit none
  
  private 

  public GlobalSetup,           &
         GlobalSetAuxVarScalar, &
         GlobalSetAuxVarVecLoc, &
         GlobalWeightAuxVars,   &
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
    allocate(auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call GlobalAuxVarInit(auxvars_ss(iconn),option)
    enddo
    patch%aux%Global%auxvars_ss => auxvars_ss
  endif
  patch%aux%Global%num_aux_ss = sum_connection
  
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

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: value
  PetscInt :: ivar

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: auxvars_bc(:)

  PetscInt :: i
  PetscInt :: num_aux
  PetscInt :: num_aux_bc

  patch => realization%patch
  option => realization%option

  auxvars => patch%aux%Global%auxvars
  auxvars_bc => patch%aux%Global%auxvars_bc
  num_aux = patch%aux%Global%num_aux
  num_aux_bc = patch%aux%Global%num_aux_bc

  select case(ivar)
    case(TEMPERATURE)
      do i=1, num_aux
        auxvars(i)%tc = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%tc = value
      enddo
    case(LIQUID_DENSITY)
      do i=1, num_aux
        auxvars(i)%den(LIQ_FLUID) = value
      enddo
    case(LIQUID_SATURATION)
      do i=1, num_aux
        auxvars(i)%sat(LIQ_FLUID) = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%sat(LIQ_FLUID) = value
      enddo
    case(LIQUID_PRESSURE)
      do i=1, num_aux
        auxvars(i)%pres(LIQ_FLUID) = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%pres(LIQ_FLUID) = value
      enddo
    case(GAS_DENSITY)
      do i=1, num_aux
        auxvars(i)%den(AIR_FLUID) = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%den(AIR_FLUID) = value
      enddo
    case(GAS_SATURATION)
      do i=1, num_aux
        auxvars(i)%sat(AIR_FLUID) = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%sat(AIR_FLUID) = value
      enddo
    case(GAS_PRESSURE)
      do i=1, num_aux
        auxvars(i)%pres(AIR_FLUID) = value
      enddo
      do i=1, num_aux_bc
        auxvars_bc(i)%pres(AIR_FLUID) = value
      enddo
    case default
      print *, 'Case(', ivar, ') not supported in GlobalSetAuxVarScalar'
      stop
  end select

end subroutine GlobalSetAuxVarScalar

! ************************************************************************** !

subroutine GlobalSetAuxVarVecLoc(realization,vec_loc,ivar,isubvar)
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
  
  implicit none

  class(realization_subsurface_type) :: realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  
  auxvars => patch%aux%Global%auxvars

  call VecGetArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)
  
  select case(ivar)
    case(TEMPERATURE)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%tc_store(TIME_T) = vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%tc_store(TIME_TpDT) = vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%tc = vec_loc_p(ghosted_id)
          enddo
      end select
    !
    case(LIQUID_SATURATION)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat_store(LIQ_FLUID,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat_store(LIQ_FLUID,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat(LIQ_FLUID) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(LIQUID_DENSITY)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_store(LIQ_FLUID,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_store(LIQ_FLUID,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den(LIQ_FLUID) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(LIQUID_PRESSURE)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres_store(LIQ_FLUID,TIME_T) = &
                vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres_store(LIQ_FLUID,TIME_TpDT) = &
                vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres(LIQ_FLUID) = vec_loc_p(ghosted_id)
          enddo
      end select
    !
    case(GAS_SATURATION)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat_store(AIR_FLUID,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat_store(AIR_FLUID,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%sat(AIR_FLUID) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(GAS_DENSITY)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_store(AIR_FLUID,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den_store(AIR_FLUID,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%den(AIR_FLUID) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(GAS_PRESSURE)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres_store(AIR_FLUID,TIME_T) = &
                vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres_store(AIR_FLUID,TIME_TpDT) = &
                vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            auxvars(ghosted_id)%pres(AIR_FLUID) = vec_loc_p(ghosted_id)
          enddo
      end select
    !
    case default
      print *, 'Case(', ivar, ') not supported in GlobalSetAuxVarVecLoc'
      stop
  end select

  call VecRestoreArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)

end subroutine GlobalSetAuxVarVecLoc

! ************************************************************************** !
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
  
  do ghosted_id = 1, realization%patch%aux%Global%num_aux
    ! interpolate density and saturation based on weight
    auxvars(ghosted_id)%tc = &
      (weight*auxvars(ghosted_id)%tc_store(TIME_TpDT)+ &
        (1.d0-weight)*auxvars(ghosted_id)%tc_store(TIME_T))
    auxvars(ghosted_id)%den(:) = &
      (weight*auxvars(ghosted_id)%den_store(:,TIME_TpDT)+ &
        (1.d0-weight)*auxvars(ghosted_id)%den_store(:,TIME_T))
    auxvars(ghosted_id)%sat(:) = &
      (weight*auxvars(ghosted_id)%sat_store(:,TIME_TpDT)+ &
        (1.d0-weight)*auxvars(ghosted_id)%sat_store(:,TIME_T))
  enddo
  
  select case(option%iflowmode) 
    case(MPFLOW_MODE)
      do ghosted_id = 1, realization%patch%aux%Global%num_aux
        auxvars(ghosted_id)%pres(:) = &
          (weight*auxvars(ghosted_id)%pres_store(:,TIME_TpDT)+ &
            (1.d0-weight)*auxvars(ghosted_id)%pres_store(:,TIME_T))
      enddo  

  end select
  
end subroutine GlobalWeightAuxVars

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
  use Communicator_Base_module
  use Material_module, only : MaterialUpdateAuxVars
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_SATURATION, &
                               LIQUID_DENSITY, GAS_PRESSURE, &
                               GAS_DENSITY, GAS_SATURATION, &
                               TEMPERATURE
  
  class(realization_subsurface_type) :: realization
  PetscReal :: time
  PetscInt :: time_level
  
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
  call RealizationGetVariable(realization,field%work,LIQUID_DENSITY, &
                             ZERO_INTEGER)
  call realization%comm1%GlobalToLocal(field%work,field%work_loc)
  call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_DENSITY, &
                             time_level)

  ! liquid saturation
  call RealizationGetVariable(realization,field%work,LIQUID_SATURATION, &
                             ZERO_INTEGER)
  call realization%comm1%GlobalToLocal(field%work,field%work_loc)
  call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_SATURATION, &
                             time_level)
  select case(option%iflowmode)

    case(MPFLOW_MODE)
      ! pressure
      call RealizationGetVariable(realization,field%work,LIQUID_PRESSURE, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_PRESSURE, &
                                 time_level)
      ! temperature
      call RealizationGetVariable(realization,field%work,TEMPERATURE, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,TEMPERATURE, &
                                 time_level)
      ! Gas density
      call RealizationGetVariable(realization,field%work,GAS_DENSITY, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_DENSITY, &
                                 time_level)
      ! Gas saturation
      call RealizationGetVariable(realization,field%work,GAS_SATURATION, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_SATURATION, &
                                 time_level)                         

  end select

end subroutine GlobalUpdateAuxVars

end module Global_module
