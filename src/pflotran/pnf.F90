module PNF_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PNF_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: PNFSetup, &
            PNFInitializeTimestep, &
            PNFUpdateSolution, &
            PNFTimeCut,&
            PNFUpdateAuxVars, &
            PNFComputeMassBalance, &
            PNFZeroMassBalanceDelta, &
            PNFSetPlotVariables, &
            PNFDestroy

contains

! ************************************************************************** !

subroutine PNFSetup(realization)
  !
  ! Creates arrays for auxiliary variables
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Material_Aux_module
  use Output_Aux_module
  use Characteristic_Curves_module
  use Matrix_Zeroing_module
  use EOS_Water_module

  implicit none

  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(output_variable_list_type), pointer :: list
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscErrorCode :: ierr
                                                ! extra index for derivatives
  type(pnf_auxvar_type), pointer :: pnf_auxvars(:)
  type(pnf_auxvar_type), pointer :: pnf_auxvars_bc(:)
  type(pnf_auxvar_type), pointer :: pnf_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%PNF => PNFAuxCreate(option)

  allocate(pnf_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call PNFAuxVarInit(pnf_auxvars(ghosted_id),option)
  enddo
  patch%aux%PNF%auxvars => pnf_auxvars
  patch%aux%PNF%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(pnf_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call PNFAuxVarInit(pnf_auxvars_bc(iconn),option)
    enddo
    patch%aux%PNF%auxvars_bc => pnf_auxvars_bc
  endif
  patch%aux%PNF%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(pnf_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call PNFAuxVarInit(pnf_auxvars_ss(iconn),option)
    enddo
    patch%aux%PNF%auxvars_ss => pnf_auxvars_ss
  endif
  patch%aux%PNF%num_aux_ss = sum_connection

  list => realization%output_option%output_snap_variable_list
  call PNFSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call PNFSetPlotVariables(realization,list)

  PNF_ts_count = 0
  PNF_ts_cut_count = 0

end subroutine PNFSetup

! ************************************************************************** !

subroutine PNFInitializeTimestep(realization)
  !
  ! Update data in module prior to time step
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Realization_Subsurface_class
  use Upwind_Direction_module

  implicit none

  type(realization_subsurface_type) :: realization


end subroutine PNFInitializeTimestep

! ************************************************************************** !

subroutine PNFUpdateSolution(realization)
  !
  ! Updates data in module after a successful time
  ! step
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization

  call PNFUpdateAuxVars(realization)
  if (realization%option%compute_mass_balance_new) then
    call PNFUpdateMassBalance(realization)
  endif

  PNF_ts_count = PNF_ts_count + 1
  PNF_ts_cut_count = 0

end subroutine PNFUpdateSolution

! ************************************************************************** !

subroutine PNFTimeCut(realization)
  !
  ! Resets arrays for time step cut
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization

  PNF_ts_cut_count = PNF_ts_cut_count + 1

  call PNFInitializeTimestep(realization)

end subroutine PNFTimeCut

! ************************************************************************** !

subroutine PNFComputeMassBalance(realization,mass_balance)
  !
  ! Initializes mass balance
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_module

  implicit none

  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(1)

end subroutine PNFComputeMassBalance

! ************************************************************************** !

subroutine PNFZeroMassBalanceDelta(realization)
  !
  ! Zeros mass balance delta array
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module

  implicit none

  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%PNF%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%PNF%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine PNFZeroMassBalanceDelta

! ************************************************************************** !

subroutine PNFUpdateMassBalance(realization)
  !
  ! Updates mass balance
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module

  implicit none

  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%PNF%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance(1,1) = &
      global_auxvars_bc(iconn)%mass_balance(1,1) + &
      global_auxvars_bc(iconn)%mass_balance_delta(1,1)*option%flow_dt
  enddo
  do iconn = 1, patch%aux%PNF%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance(1,1) = &
      global_auxvars_ss(iconn)%mass_balance(1,1) + &
      global_auxvars_ss(iconn)%mass_balance_delta(1,1)*option%flow_dt
  enddo

end subroutine PNFUpdateMassBalance

! ************************************************************************** !

subroutine PNFUpdateAuxVars(realization)
  !
  ! Updates the auxiliary variables associated with the PNF problem
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Discretization_module
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Material_Aux_module
  use General_Aux_module, only : ANY_STATE, TWO_PHASE_STATE

  implicit none

  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(pnf_auxvar_type), pointer :: pnf_auxvars(:)
  type(pnf_auxvar_type), pointer :: pnf_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, iconn, natural_id
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  pnf_auxvars => patch%aux%PNF%auxvars
  pnf_auxvars_bc => patch%aux%PNF%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  call DiscretizationGlobalToLocal(realization%discretization,field%flow_xx, &
                                   field%flow_xx_loc,NFLOWDOF)
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! PNF_UPDATE_FOR_ACCUM indicates call from non-perturbation
    natural_id = grid%nG2A(ghosted_id)
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
    call PNFAuxVarCompute(xx_loc_p(ghosted_id:ghosted_id), &
                            pnf_auxvars(ghosted_id), &
                            global_auxvars(ghosted_id), &
                            natural_id, &
                            option)
  enddo

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      !geh: negate to indicate boundary connection, not actual cell
      natural_id = -grid%nG2A(ghosted_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(boundary_condition%flow_condition% &
                    itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC)
          xxbc(1) = boundary_condition% &
                      flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
        case(NEUMANN_BC)
          xxbc(1) = xx_loc_p(local_id)
        case default
          option%io_buffer = 'boundary itype not set up in PNFUpdateAuxVars'
          call PrintErrMsg(option)
      end select

      ! PNF_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      call PNFAuxVarCompute(xxbc,pnf_auxvars_bc(sum_connection), &
                            global_auxvars_bc(sum_connection), &
                            natural_id, &
                            option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%PNF%auxvars_up_to_date = PETSC_TRUE

end subroutine PNFUpdateAuxVars

! ************************************************************************** !

subroutine PNFSetPlotVariables(realization,list)
  !
  ! Adds variables to be printed to list
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module

  implicit none

  type(realization_subsurface_type) :: realization
  type(output_variable_list_type), pointer :: list

  character(len=MAXWORDLENGTH) :: name, units

  if (associated(list%first)) then
    return
  endif

  if (list%flow_vars) then

    name = 'Liquid Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                LIQUID_PRESSURE)

  endif

end subroutine PNFSetPlotVariables

! ************************************************************************** !

subroutine PNFDestroy(realization)
  !
  ! Deallocates variables associated with Richard
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Option_module

  implicit none

  type(realization_subsurface_type) :: realization

  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary.F90.

end subroutine PNFDestroy

end module PNF_module
