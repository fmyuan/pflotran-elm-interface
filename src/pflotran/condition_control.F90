module Condition_Control_module
  
  ! This module store routines that operate on conditions from a level above 
  ! that of the realization_module.  This is necessary to access capability
  ! such as HDF5 which is unavailable from within the realization object
  ! and below.  Routines in this module will loop over realization, levels,
  ! and patches without calling underlying level/patch versions of the
  ! subroutines, which is common in realization.F90 - GEH
#include "petsc/finclude/petscvec.h"
  use petscvec

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: CondControlAssignFlowInitCond, &
            !CondControlAssignTranInitCond, &
            CondControlScaleSourceSink
 
contains

! ************************************************************************** !

subroutine CondControlAssignFlowInitCond(realization)
  ! 
  ! Assigns flow initial conditions to model
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07, 10/18/11
  ! 
  use Realization_Subsurface_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use FlowCondition_module
  use Dataset_Base_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Common_HDF5_class
  use Dataset_module
  use Grid_module
  use Patch_module
  use EOS_Water_module

  use Global_module
  use Global_Aux_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscInt :: iconn, idof
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscReal, pointer :: xx_p(:)
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(patch_type), pointer :: cur_patch
  class(dataset_base_type), pointer :: dataset
  PetscBool :: dataset_flag(realization%option%nflowdof)
  PetscInt :: num_connections
  PetscInt, pointer :: conn_id_ptr(:)

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    grid => cur_patch%grid

    select case(option%iflowmode)

      case default
        ! assign initial conditions values to domain
        call VecGetArrayF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
      
        xx_p = UNINITIALIZED_DOUBLE
      
        initial_condition => cur_patch%initial_condition_list%first
        do
      
          if (.not.associated(initial_condition)) exit

          dataset_flag = PETSC_FALSE
          do idof = 1, option%nflowdof
            dataset =>  initial_condition%flow_condition% &
                              sub_condition_ptr(idof)%ptr%dataset
            select type(dataset_ptr => dataset)
              class is(dataset_gridded_hdf5_type)
                ! already mapped to flow_aux_real_var
                if (.not.associated(initial_condition%flow_aux_real_var)) then
                  option%io_buffer = 'A gridded dataset is being &
                    &used with ' // trim(option%flowmode) // &
                    ', yet flow_aux_real_var is not allocated.'
                  call PrintErrMsgToDev(option,'')
                endif
              class is(dataset_common_hdf5_type)
                dataset_flag(idof) = PETSC_TRUE
                call VecRestoreArrayF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
                call ConditionControlMapDatasetToVec(realization, &
                        initial_condition%flow_condition% &
                          sub_condition_ptr(idof)%ptr%dataset,idof, &
                        field%flow_xx,GLOBAL)
                call VecGetArrayF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
            end select
          enddo            
          if (.not.associated(initial_condition%flow_aux_real_var) .and. &
              .not.associated(initial_condition%flow_condition)) then
            option%io_buffer = 'Flow condition is NULL in initial condition'
            call PrintErrMsg(option)
          endif
          if (associated(initial_condition%flow_aux_real_var)) then
            num_connections = &
              initial_condition%connection_set%num_connections
            conn_id_ptr => initial_condition%connection_set%id_dn
          else
            num_connections = initial_condition%region%num_cells
            conn_id_ptr => initial_condition%region%cell_ids
          endif
          do iconn=1, num_connections
            local_id = conn_id_ptr(iconn)
            ghosted_id = grid%nL2G(local_id)
            iend = local_id*option%nflowdof
            ibegin = iend-option%nflowdof+1
            if (cur_patch%imat(ghosted_id) <= 0) then
              xx_p(ibegin:iend) = 0.d0
            endif
            if (associated(initial_condition%flow_aux_real_var)) then
              do idof = 1, option%nflowdof
                if (.not.dataset_flag(idof)) then
                  xx_p(ibegin+idof-1) =  &
                    initial_condition%flow_aux_real_var(idof,iconn)
                endif
              enddo
            else
              do idof = 1, option%nflowdof
                if (.not.dataset_flag(idof)) then
                  xx_p(ibegin+idof-1) = &
                    initial_condition%flow_condition% &
                      sub_condition_ptr(idof)%ptr%dataset%rarray(1)
                endif
              enddo
            endif
          enddo
          initial_condition => initial_condition%next
        enddo
        call VecRestoreArrayF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

    end select 
   
    cur_patch => cur_patch%next
  enddo

  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                   field%flow_xx_loc,NFLOWDOF)  

  call VecCopy(field%flow_xx, field%flow_yy, ierr);CHKERRQ(ierr)

end subroutine CondControlAssignFlowInitCond

! ************************************************************************** !

! ************************************************************************** !

subroutine ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                           mdof_vec,vec_type)
  ! 
  ! maps an external dataset to a PETSc vec
  ! representing values at each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/23/12
  ! 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Dataset_Common_HDF5_class
  use Dataset_Base_class
  use HDF5_module
  use Discretization_module

  implicit none
  

  class(realization_subsurface_type) :: realization
  class(dataset_base_type), pointer :: dataset
  PetscInt :: idof
  Vec :: mdof_vec
  PetscInt :: vec_type

  type(field_type), pointer :: field
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  field => realization%field
  option => realization%option
  
  call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
  if (associated(dataset)) then
    select type(dataset)
      class is (dataset_common_hdf5_type)
        string = '' ! group name
        ! have to copy to string2 due to mismatch in string size
        string2 = dataset%hdf5_dataset_name
        call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                          dataset%filename, &
                                          string,string2, &
                                          dataset%realization_dependent)
        if (vec_type == GLOBAL) then
          call VecStrideScatter(field%work,idof-1,mdof_vec, &
                                INSERT_VALUES,ierr);CHKERRQ(ierr)
        else
          call DiscretizationGlobalToLocal(realization%discretization, &
                                           field%work, &
                                           field%work_loc,ONEDOF)
          call VecStrideScatter(field%work_loc,idof-1,mdof_vec, &
                                INSERT_VALUES,ierr);CHKERRQ(ierr)
        endif
      class default
        option%io_buffer = 'Dataset "' // trim(dataset%name) // &
          '" not supported in ConditionControlMapDatasetToVec.'
        call PrintErrMsg(option)
    end select
  endif

end subroutine ConditionControlMapDatasetToVec

! ************************************************************************** !

subroutine CondControlScaleSourceSink(realization)
  ! 
  ! Scales select source/sinks based on perms
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/03/08, 10/18/11
  ! 
#include "petsc/finclude/petscdmda.h"
  use petscdmda
      
  use Realization_Subsurface_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use FlowCondition_module
  use Grid_module
  use Patch_module
  use Material_Aux_class
  use Variables_module, only : PERMEABILITY_X

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: cur_source_sink
  type(connection_set_type), pointer :: cur_connection_set
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(patch_type), pointer :: cur_patch
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id, neighbor_ghosted_id
  PetscInt :: iconn
  PetscReal :: scale, sum
  PetscInt :: icount
  PetscInt :: x_count, y_count, z_count
  PetscInt, parameter :: x_width = 1, y_width = 1, z_width = 0
  
  PetscInt :: ghosted_neighbors(0:27)
  
  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  material_auxvars => realization%patch%aux%Material%auxvars
  
  ! GB: grid was uninitialized
  grid => patch%grid

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    ! BIG-TIME warning here.  I assume that all source/sink cells are within 
    ! a single patch - geh

    grid => cur_patch%grid

    cur_source_sink => cur_patch%source_sink_list%first
    do
      if (.not.associated(cur_source_sink)) exit

      call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

      cur_connection_set => cur_source_sink%connection_set
    
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)

        select case(option%iflowmode)
          case(MPFLOW_MODE)
              call GridGetGhostedNeighbors(grid,ghosted_id,DMDA_STENCIL_STAR, &
                                          x_width,y_width,z_width, &
                                          x_count,y_count,z_count, &
                                          ghosted_neighbors,option)
              ! ghosted neighbors is ordered first in x, then, y, then z
              icount = 0
              sum = 0.d0
              ! x-direction
              do while (icount < x_count)
                icount = icount + 1
                neighbor_ghosted_id = ghosted_neighbors(icount)
                sum = sum + MaterialAuxVarGetValue(material_auxvars( &
                              neighbor_ghosted_id),PERMEABILITY_X) * &
                            grid%structured_grid%dy(neighbor_ghosted_id)* &
                            grid%structured_grid%dz(neighbor_ghosted_id)
                 
              enddo
              ! y-direction
              do while (icount < x_count + y_count)
                icount = icount + 1
                neighbor_ghosted_id = ghosted_neighbors(icount)                 
                sum = sum + MaterialAuxVarGetValue(material_auxvars( &
                              neighbor_ghosted_id),PERMEABILITY_X) * &
                            grid%structured_grid%dx(neighbor_ghosted_id)* &
                            grid%structured_grid%dz(neighbor_ghosted_id)
                 
              enddo
              ! z-direction
              do while (icount < x_count + y_count + z_count)
                icount = icount + 1
                neighbor_ghosted_id = ghosted_neighbors(icount)                 
                sum = sum + MaterialAuxVarGetValue(material_auxvars( &
                              neighbor_ghosted_id),PERMEABILITY_X) * &
                            grid%structured_grid%dx(neighbor_ghosted_id)* &
                            grid%structured_grid%dy(neighbor_ghosted_id)
              enddo
              vec_ptr(local_id) = vec_ptr(local_id) + sum
        end select

      enddo
        
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      call VecNorm(field%work,NORM_1,scale,ierr);CHKERRQ(ierr)
      scale = 1.d0/scale
      call VecScale(field%work,scale,ierr);CHKERRQ(ierr)

      call VecGetArrayF90(field%work,vec_ptr, ierr);CHKERRQ(ierr)
      do iconn = 1, cur_connection_set%num_connections      
        local_id = cur_connection_set%id_dn(iconn)
        select case(option%iflowmode)
          case(MPFLOW_MODE)
            cur_source_sink%flow_aux_real_var(ONE_INTEGER,iconn) = &
              vec_ptr(local_id)
        end select

      enddo
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
        
      cur_source_sink => cur_source_sink%next
    enddo
    cur_patch => cur_patch%next
  enddo

end subroutine CondControlScaleSourceSink

! ************************************************************************** !

end module Condition_Control_module
