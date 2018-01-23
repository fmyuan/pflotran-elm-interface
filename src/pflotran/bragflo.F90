module Bragflo_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use WIPP_Flow_module
  use WIPP_Flow_Common_module
  use WIPP_Flow_Aux_module
  use Bragflo_Common_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  public :: BragfloResidual, &
            BragfloJacobian

contains

! ************************************************************************** !

subroutine BragfloResidual(snes,xx,r,realization,pmwss_ptr,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module

  use Connection_module
  use Grid_module
  use Coupler_module  
  use Debug_module
  use Material_Aux_class
  use PM_WIPP_SrcSink_class

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  class(pm_wipp_srcsink_type), pointer :: pmwss_ptr
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  Mat, parameter :: null_mat = tMat(-1)
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(material_parameter_type), pointer :: material_parameter
  type(wippflo_parameter_type), pointer :: wippflo_parameter
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscReal :: scale
  PetscReal :: ss_flow_vol_flux(realization%option%nphase)
  PetscInt :: sum_connection
  PetscInt :: local_start, local_end
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, imat, imat_up, imat_dn
  PetscInt :: iteration_number
  PetscInt, pointer :: upwind_direction(:,:), upwind_direction_bc(:,:)

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  PetscReal, pointer :: vec_p(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k

  PetscInt :: icap_up, icap_dn
  PetscReal :: Res(realization%option%nflowdof)
  PetscReal :: Jac_dummy(realization%option%nflowdof, &
                         realization%option%nflowdof)
  PetscReal :: v_darcy(realization%option%nphase)
  
  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  wippflo_auxvars_bc => patch%aux%WIPPFlo%auxvars_bc
  wippflo_parameter => patch%aux%WIPPFlo%wippflo_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  upwind_direction => patch%aux%WIPPFlo%upwind_direction
  upwind_direction_bc => patch%aux%WIPPFlo%upwind_direction_bc

  wippflo_newton_iteration_number = wippflo_newton_iteration_number + 1
  ! bragflo uses the following logic, update when
  !   it == 1, before entering iteration loop
  !   it > 1 and mod(it-1,frequency) == 0
  ! the first is set in WIPPFloInitializeTimestep, the second is set here
  if (wippflo_newton_iteration_number > 1 .and. &
      mod(wippflo_newton_iteration_number-1, &
          wippflo_upwind_dir_update_freq) == 0) then
    wippflo_update_upwind_direction = PETSC_TRUE
  endif

  ! Communication -----------------------------------------
  ! must be called before WIPPFloUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call WIPPFloUpdateAuxVars(realization)

  ! override flags since they will soon be out of date
  patch%aux%WIPPFlo%auxvars_up_to_date = PETSC_FALSE 

  if (option%compute_mass_balance_new) then
    call WIPPFloZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  r_p = -accum_p
  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  
  ! accumulation at t(k+1)
  call VecGetArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    call WIPPFloAccumulation(wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,Res,Jac_dummy, &
                             PETSC_FALSE)
    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    accum_p2(local_start:local_end) = Res(:)
  enddo
  call VecRestoreArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      imat_up = patch%imat(ghosted_id_up) 
      imat_dn = patch%imat(ghosted_id_dn) 
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)

      call BragfloFlux(wippflo_auxvars(ZERO_INTEGER,ghosted_id_up), &
                       global_auxvars(ghosted_id_up), &
                       material_auxvars(ghosted_id_up), &
                       wippflo_auxvars(ZERO_INTEGER,ghosted_id_dn), &
                       global_auxvars(ghosted_id_dn), &
                       material_auxvars(ghosted_id_dn), &
                       cur_connection_set%area(iconn), &
                       cur_connection_set%dist(:,iconn), &
                       upwind_direction(:,sum_connection), &
                       wippflo_parameter,option,v_darcy,Res, &
                       Jac_dummy,Jac_dummy, &
                       PETSC_FALSE, & ! derivative call
                       wippflo_fix_upwind_direction, &
                       wippflo_update_upwind_direction, &
                       wippflo_count_upwind_dir_flip, & 
                       PETSC_FALSE)

      patch%internal_velocities(:,sum_connection) = v_darcy
      if (associated(patch%internal_flow_fluxes)) then
        patch%internal_flow_fluxes(:,sum_connection) = Res(:)
      endif
      
      if (local_id_up > 0) then
        local_end = local_id_up * option%nflowdof
        local_start = local_end - option%nflowdof + 1
        r_p(local_start:local_end) = r_p(local_start:local_end) + Res(:)
      endif
         
      if (local_id_dn > 0) then
        local_end = local_id_dn * option%nflowdof
        local_start = local_end - option%nflowdof + 1
        r_p(local_start:local_end) = r_p(local_start:local_end) - Res(:)
      endif
    enddo

    cur_connection_set => cur_connection_set%next
  enddo    

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      imat_dn = patch%imat(ghosted_id)
      if (imat_dn <= 0) cycle

      icap_dn = patch%sat_func_id(ghosted_id)

      call BragfloBCFlux(boundary_condition%flow_bc_type, &
                     boundary_condition%flow_aux_mapping, &
                     boundary_condition%flow_aux_real_var(:,iconn), &
                     wippflo_auxvars_bc(sum_connection), &
                     global_auxvars_bc(sum_connection), &
                     wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     upwind_direction_bc(:,sum_connection), &
                     wippflo_parameter,option, &
                     v_darcy,Res,Jac_dummy, &
                     PETSC_FALSE, & ! derivative call
                     wippflo_fix_upwind_direction, &
                     wippflo_update_upwind_direction, &
                     wippflo_count_upwind_dir_flip, & 
                     PETSC_FALSE)
      patch%boundary_velocities(:,sum_connection) = v_darcy
      if (associated(patch%boundary_flow_fluxes)) then
        patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2)
      endif

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1
      r_p(local_start:local_end)= r_p(local_start:local_end) - Res(:)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif
      
      call WIPPFloSrcSink(option,source_sink%flow_condition%general%rate% &
                                  dataset%rarray(:), &
                          source_sink%flow_condition%general%rate%itype, &
                          wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          ss_flow_vol_flux, &
                          scale,Res,Jac_dummy, &
                          PETSC_FALSE)

      r_p(local_start:local_end) =  r_p(local_start:local_end) - Res(:)

      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif      
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)
      endif      
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2)
      endif

    enddo
    source_sink => source_sink%next
  enddo
  
  ! WIPP gas/brine generation process model source/sinks
  if (associated(pmwss_ptr)) then
    call pmwss_ptr%Solve(option%time,ierr)
    call PMWSSCalcResidualValues(pmwss_ptr,local_start,local_end, &
                                 r_p,ss_flow_vol_flux)    
  endif

  if (patch%aux%WIPPFlo%inactive_cells_exist) then
    do i=1,patch%aux%WIPPFlo%n_inactive_rows
      r_p(patch%aux%WIPPFlo%inactive_rows_local(i)) = 0.d0
    enddo
  endif
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  
  call WIPPFloSSSandbox(r,null_mat,PETSC_FALSE,grid,material_auxvars, &
                        wippflo_auxvars,option)

  ! Mass Transfer
  if (field%flow_mass_transfer /= PETSC_NULL_VEC) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    !call VecGetArrayF90(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    !call VecRestoreArrayF90(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    call VecAXPY(r,-1.d0,field%flow_mass_transfer,ierr);CHKERRQ(ierr)
  endif                      
                        
  if (realization%debug%vecview_residual) then
    string = 'BFresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    string = 'BFxx'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  wippflo_update_upwind_direction = PETSC_FALSE

end subroutine BragfloResidual

! ************************************************************************** !

subroutine BragfloJacobian(snes,xx,A,B,realization,pmwss_ptr,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_class
  use PM_WIPP_SrcSink_class

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization
  class(pm_wipp_srcsink_type), pointer :: pmwss_ptr
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscReal :: norm
  PetscViewer :: viewer

  PetscInt :: icap_up,icap_dn
  PetscReal :: qsrc, scale
  PetscInt :: imat, imat_up, imat_dn
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  Vec, parameter :: null_vec = tVec(-1)
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscInt, pointer :: zeros(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(material_parameter_type), pointer :: material_parameter
  type(wippflo_parameter_type), pointer :: wippflo_parameter
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt, pointer :: upwind_direction(:,:), upwind_direction_bc(:,:)
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  wippflo_auxvars_bc => patch%aux%WIPPFlo%auxvars_bc
  wippflo_parameter => patch%aux%WIPPFlo%wippflo_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  upwind_direction => patch%aux%WIPPFlo%upwind_direction
  upwind_direction_bc => patch%aux%WIPPFlo%upwind_direction_bc

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  ! Perturb aux vars
  do ghosted_id = 1, grid%ngmax  ! For each local node do...
    if (patch%imat(ghosted_id) <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    call WIPPFloAuxVarPerturb(wippflo_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              natural_id,option)
  enddo
  
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    call WIPPFloAccumDerivative(wippflo_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              material_parameter%soil_heat_capacity(imat), &
                              option, &
                              Jup) 
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo

  ! Interior Flux Terms -----------------------------------  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0    
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      imat_up = patch%imat(ghosted_id_up)
      imat_dn = patch%imat(ghosted_id_dn)
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
   
      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)
                              
      call BragfloFluxDerivative(wippflo_auxvars(:,ghosted_id_up), &
                     global_auxvars(ghosted_id_up), &
                     material_auxvars(ghosted_id_up), &
                     wippflo_auxvars(:,ghosted_id_dn), &
                     global_auxvars(ghosted_id_dn), &
                     material_auxvars(ghosted_id_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     upwind_direction(:,sum_connection), &
                     wippflo_parameter,option, &
                     Jup,Jdn)
      if (local_id_up > 0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      imat_dn = patch%imat(ghosted_id)
      if (imat_dn <= 0) cycle

      icap_dn = patch%sat_func_id(ghosted_id)

      call BragfloBCFluxDerivative(boundary_condition%flow_bc_type, &
                      boundary_condition%flow_aux_mapping, &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      wippflo_auxvars_bc(sum_connection), &
                      global_auxvars_bc(sum_connection), &
                      wippflo_auxvars(:,ghosted_id), &
                      global_auxvars(ghosted_id), &
                      material_auxvars(ghosted_id), &
                      cur_connection_set%area(iconn), &
                      cur_connection_set%dist(:,iconn), &
                      upwind_direction_bc(:,sum_connection), &
                      wippflo_parameter,option, &
                      Jdn)

      Jdn = -Jdn
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sinks
  source_sink => patch%source_sink_list%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif
      
      Jup = 0.d0
      call WIPPFloSrcSinkDerivative(option, &
                        source_sink%flow_condition%general%rate% &
                                  dataset%rarray(:), &
                        source_sink%flow_condition%general%rate%itype, &
                        wippflo_auxvars(:,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        scale,Jup)

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    source_sink => source_sink%next
  enddo
  
  
  ! WIPP gas/brine generation process model source/sinks
  if (associated(pmwss_ptr)) then
    call PMWSSCalcJacobianValues(pmwss_ptr,A,ierr)
  endif
  
  call WIPPFloSSSandbox(null_vec,A,PETSC_TRUE,grid,material_auxvars, &
                        wippflo_auxvars,option)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! zero out inactive cells
  if (patch%aux%WIPPFlo%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,patch%aux%WIPPFlo%n_inactive_rows, &
                          patch%aux%WIPPFlo%inactive_rows_local_ghosted, &
                          qsrc,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif

  if (realization%debug%matview_Jacobian) then
    string = 'BFjacobian'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option) 
  endif

end subroutine BragfloJacobian

end module Bragflo_module
