module Richards_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Richards_Aux_module
  use Richards_Common_module
  use Global_Aux_module
  use Material_Aux_class
  use InlineSurface_Aux_module
  use InlineSurface_module
#ifdef BUFFER_MATRIX
  use Matrix_Buffer_module
#endif
  
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none
  
  private 

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-6
  PetscReal, parameter :: unit_z(3) = [0.d0,0.d0,1.d0]
  PetscInt, parameter :: NO_CONN = 1
  PetscInt, parameter :: VERT_CONN = 2
  PetscInt, parameter :: HORZ_CONN = 3

  public RichardsResidual, &
         RichardsJacobian, &
         RichardsTimeCut,&
         RichardsSetup, &
         RichardsInitializeTimestep, &
         RichardsUpdateAuxVars, &
         RichardsMaxChange, &
         RichardsUpdateSolution, &
         RichardsComputeMassBalance, &
         RichardsDestroy, &
         RichardsResidualInternalConn, &
         RichardsResidualBoundaryConn, &
         RichardsResidualSourceSink, &
         RichardsJacobianInternalConn, &
         RichardsJacobianBoundaryConn, &
         RichardsJacobianSourceSink

contains

! ************************************************************************** !

subroutine RichardsTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  
  richards_ts_cut_count = richards_ts_cut_count + 1
  call RichardsInitializeTimestep(realization)  
 
end subroutine RichardsTimeCut

! ************************************************************************** !

subroutine RichardsSetup(realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Output_Aux_module

  type(realization_subsurface_type) :: realization

  type(output_variable_list_type), pointer :: list
  
  call RichardsSetupPatch(realization)

  list => realization%output_option%output_snap_variable_list
  call RichardsSetPlotVariables(list)
  list => realization%output_option%output_obs_variable_list
  call RichardsSetPlotVariables(list)
  
  richards_ts_count = 0
  richards_ts_cut_count = 0
  richards_ni_count = 0

end subroutine RichardsSetup

! ************************************************************************** !

subroutine RichardsSetupPatch(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Region_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink

  PetscInt :: local_id, ghosted_id, iconn, sum_connection
  PetscInt :: ivertex, nvert, region_id, vertex_id
  PetscInt :: i
  PetscBool :: error_found, found
  PetscInt :: flag(10)
  PetscReal :: minz, maxz, zcenter, zface
  PetscReal :: tempreal
  PetscErrorCode :: ierr
  type(material_parameter_type), pointer :: material_parameter
  class(material_auxvar_type), pointer :: material_auxvars(:)  
  type(richards_auxvar_type), pointer :: rich_auxvars(:)  
  type(richards_auxvar_type), pointer :: rich_auxvars_bc(:)  
  type(richards_auxvar_type), pointer :: rich_auxvars_ss(:)  
  type(region_type), pointer :: region
  type(coupler_type), pointer :: coupler
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%Richards => RichardsAuxCreate()

  ! ensure that material properties specific to this module are properly
  ! initialized
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE
  material_auxvars => patch%aux%Material%auxvars
  flag = 0
  !TODO(geh): change to looping over ghosted ids once the legacy code is 
  !           history and the communicator can be passed down.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (material_auxvars(ghosted_id)%volume < 0.d0 .and. flag(1) == 0) then
      flag(1) = 1
      option%io_buffer = 'ERROR: Non-initialized cell volume.'
      call PrintMsgByRank(option)
    endif
    if (material_auxvars(ghosted_id)%porosity_base < 0.d0 .and. &
        flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'ERROR: Non-initialized porosity.'
      call PrintMsgByRank(option)
    endif
    tempreal = minval(material_auxvars(ghosted_id)%permeability)
    if (Uninitialized(tempreal) .and. flag(5) == 0) then
      option%io_buffer = 'ERROR: Non-initialized permeability.'
      call PrintMsgByRank(option)
      flag(5) = 1
    endif
  enddo

  error_found = error_found .or. (maxval(flag) > 0)
  call MPI_Allreduce(MPI_IN_PLACE,error_found,ONE_INTEGER_MPI,MPI_LOGICAL, &
                     MPI_LOR,option%mycomm,ierr)
  if (error_found) then
    option%io_buffer = 'Material property errors found in RichardsSetup.'
    call PrintErrMsg(option)
  endif
  
  ! allocate auxvar data structures for all grid cells  
  allocate(rich_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call RichardsAuxVarInit(rich_auxvars(ghosted_id),option)
  enddo
  patch%aux%Richards%auxvars => rich_auxvars
  patch%aux%Richards%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(rich_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call RichardsAuxVarInit(rich_auxvars_bc(iconn),option)
    enddo
    patch%aux%Richards%auxvars_bc => rich_auxvars_bc
  endif
  patch%aux%Richards%num_aux_bc = sum_connection
  
  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(rich_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call RichardsAuxVarInit(rich_auxvars_ss(iconn),option)
    enddo
    patch%aux%Richards%auxvars_ss => rich_auxvars_ss
  endif
  patch%aux%Richards%num_aux_ss = sum_connection

  ! if we are doing inline surface flow, build those auxillary
  ! variables for grid cells as well as boundary conditions
  if (option%flow%inline_surface_flow) then

    ! point to the top cell region
    region => RegionGetPtrFromList(option%flow%inline_surface_region_name, &
                                   patch%region_list)

    ! create the auxvar structure
    patch%aux%InlineSurface => InlineSurfaceAuxCreate()

    ! allocate auxvar data structures for the region grid cells, we
    ! also need to compute a few extra quantities that aren't part of
    ! the standard structures.
    allocate(patch%aux%InlineSurface%auxvars(region%num_cells))
    patch%aux%InlineSurface%num_aux = region%num_cells
    do region_id = 1, patch%aux%InlineSurface%num_aux

      call InlineSurfaceAuxVarInit(patch%aux%InlineSurface%auxvars(region_id), &
                                   option)

      ! find the half cell height, this is used in the approximation
      ghosted_id = region%cell_ids(region_id)
      if (associated(grid%unstructured_grid)) then
        if (associated(grid%unstructured_grid%explicit_grid)) then
          zcenter = &
            grid%unstructured_grid%explicit_grid%cell_centroids(ghosted_id)%z 
          zface = region%explicit_faceset%face_centroids(region_id)%z  
          patch%aux%InlineSurface%auxvars(region_id)%half_cell_height = &
            zface-zcenter
        else if( associated(grid%unstructured_grid%polyhedra_grid) ) then
          option%io_buffer = 'richards.F90:RichardsSetupPatch() --> &
             &unsupported grid type, could not compute top cell half heights.'
          call PrintErrMsg(option)
        else !implicit unstructured 
          minz  = 0.0d0
          maxz  = 0.0d0
          nvert = grid%unstructured_grid%cell_vertices(0,ghosted_id)/2
          do ivertex = 1, nvert
            vertex_id = grid%unstructured_grid%cell_vertices(ivertex,ghosted_id)
            minz = minz + grid%unstructured_grid%vertices(vertex_id)%z
            vertex_id = grid%unstructured_grid%cell_vertices(ivertex+nvert, &
                                                             ghosted_id)
            maxz = maxz + grid%unstructured_grid%vertices(vertex_id)%z
          enddo
          minz = minz/(DBLE(nvert))
          maxz = maxz/(DBLE(nvert))
          patch%aux%InlineSurface%auxvars(region_id)%half_cell_height = &
            abs(0.5d0*(maxz-minz))
        end if
      else if (associated(grid%structured_grid)) then
        patch%aux%InlineSurface%auxvars(region_id)%half_cell_height = &
          0.5d0*grid%structured_grid%dz(ghosted_id)
      else
        option%io_buffer = 'richards.F90:RichardsSetupPatch() --> &
          &unsupported grid type, could not compute top cell half heights.'
        call PrintErrMsg(option)
      endif

      ! set Manning's coefficient
      patch%aux%InlineSurface%auxvars(region_id)%Mannings_coeff = &
        option%flow%inline_surface_Mannings_coeff

    enddo

    ! loop over bc's, if a surface bc, repeat the above     
    sum_connection = 0
    coupler => patch%boundary_condition_list%first
    do
      if (.not.associated(coupler)) exit
      if (coupler%flow_condition%pressure%itype == SURFACE_DIRICHLET .or. &
          coupler%flow_condition%pressure%itype == &
                                             SURFACE_ZERO_GRADHEIGHT .or. &
          coupler%flow_condition%pressure%itype == SURFACE_SPILLOVER ) then
        do iconn = 1,coupler%connection_set%num_connections

          ! the connection down cell must be in the surface region
          found = PETSC_FALSE
          do region_id = 1, region%num_cells
            if (coupler%connection_set%id_dn(iconn) == &
                  region%cell_ids(region_id)) then
              found = PETSC_TRUE
            endif
          enddo
          if (found .eqv. PETSC_FALSE) then
            option%io_buffer = 'richards.F90:RichardsSetupPatch() --> &
              &surface boundary condition, assigned to a region which is &
              &not the boundary of the inline surface region.'
            call PrintErrMsg(option)
          endif
        enddo
        sum_connection = sum_connection + coupler%connection_set%num_connections
      endif
      coupler => coupler%next
    enddo

    ! if we have something to allocate, then allocate it!
    if (sum_connection > 0) then
      allocate(patch%aux%InlineSurface%auxvars_bc(sum_connection))
      patch%aux%InlineSurface%num_aux_bc = sum_connection
      do iconn = 1, sum_connection
        call InlineSurfaceAuxVarInit(patch%aux%InlineSurface%&
                                       auxvars_bc(iconn),option)
        patch%aux%InlineSurface%auxvars_bc(iconn)%Mannings_coeff = &
          option%flow%inline_surface_Mannings_coeff
      enddo
    endif

    ! I will just get the half cell height from the interior auxvar
    sum_connection = 0
    coupler => patch%boundary_condition_list%first
    do
      if (.not.associated(coupler)) exit
      if (coupler%flow_condition%pressure%itype == SURFACE_DIRICHLET .or. &
          coupler%flow_condition%pressure%itype == &
                                             SURFACE_ZERO_GRADHEIGHT .or. &
          coupler%flow_condition%pressure%itype == SURFACE_SPILLOVER) then
        do iconn = 1,coupler%connection_set%num_connections
          sum_connection = sum_connection + 1
          do region_id = 1, region%num_cells
            if (coupler%connection_set%id_dn(iconn) == &
                region%cell_ids(region_id)) then
              patch%aux%InlineSurface%auxvars_bc(sum_connection)% &
                half_cell_height = &
                patch%aux%InlineSurface%auxvars(region_id)%half_cell_height
              cycle
            endif
          enddo
        enddo
      endif
      coupler => coupler%next
    enddo

  endif
  
end subroutine RichardsSetupPatch

! ************************************************************************** !

subroutine RichardsComputeMassBalance(realization,mass_balance)
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Realization_Subsurface_class

  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)
  
  mass_balance = 0.d0
  
  call RichardsComputeMassBalancePatch(realization,mass_balance)

end subroutine RichardsComputeMassBalance

! ************************************************************************** !

subroutine RichardsComputeMassBalancePatch(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Region_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(inlinesurface_auxvar_type), pointer :: inlinesurface_auxvars(:)
  type(region_type), pointer :: region
  
  PetscErrorCode :: ierr
  PetscInt :: local_id, region_id
  PetscInt :: ghosted_id
  PetscReal :: mass
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! mass = volume*saturation*density
    mass_balance = mass_balance + &
      global_auxvars(ghosted_id)%den_kg* &
      global_auxvars(ghosted_id)%sat* &
      material_auxvars(ghosted_id)%porosity* &
      material_auxvars(ghosted_id)%volume
  enddo

  if (option%flow%inline_surface_flow) then
    inlinesurface_auxvars => patch%aux%InlineSurface%auxvars
    region => RegionGetPtrFromList(option%flow%inline_surface_region_name, &
                                   patch%region_list)
    do region_id = 1, patch%aux%InlineSurface%num_aux
      ghosted_id = region%cell_ids(region_id)
      mass = inlinesurface_auxvars(region_id)%surface_water_depth
      mass = mass * global_auxvars(ghosted_id)%den_kg(1)
      mass = mass * material_auxvars(ghosted_id)%volume
      mass = mass / (2.0d0*inlinesurface_auxvars(region_id)%half_cell_height)
      mass_balance = mass_balance + mass
    enddo
  endif
  
end subroutine RichardsComputeMassBalancePatch

! ************************************************************************** !

subroutine RichardsZeroMassBalDeltaPatch(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
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

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%Richards%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%Richards%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
  if (patch%aux%Richards%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_ss
      global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
    enddo
  endif

end subroutine RichardsZeroMassBalDeltaPatch

! ************************************************************************** !

subroutine RichardsUpdateMassBalancePatch(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
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

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%Richards%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance = &
      patch%aux%Global%auxvars(iconn)%mass_balance + &
      patch%aux%Global%auxvars(iconn)%mass_balance_delta*FMWH2O* &
      option%flow_dt
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%Richards%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance = &
        global_auxvars_bc(iconn)%mass_balance + &
        global_auxvars_bc(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif

  if (patch%aux%Richards%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_ss
      global_auxvars_ss(iconn)%mass_balance = &
        global_auxvars_ss(iconn)%mass_balance + &
        global_auxvars_ss(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif

end subroutine RichardsUpdateMassBalancePatch

! ************************************************************************** !

subroutine RichardsUpdatePermPatch(realization)
  ! 
  ! Updates the permeability based on pressure
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 01/09/12
  ! 

  use Grid_module
  use Realization_Subsurface_class
  use Option_module
  use Discretization_module
  use Patch_module
  use Field_module
  use Material_module
  use Material_Aux_class
  use Variables_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(material_property_ptr_type), pointer :: material_property_array(:)
  type(discretization_type), pointer :: discretization
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscReal :: scale
  PetscReal :: p_min, p_max, permfactor_max
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscReal, pointer :: perm_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  field => realization%field
  grid => patch%grid
  material_property_array => patch%material_property_array
  material_auxvars => patch%aux%Material%auxvars

  if (.not.associated(patch%imat)) then
    option%io_buffer = 'Materials IDs not present in run.  Material ' // &
      ' properties cannot be updated without material ids'
    call PrintErrMsg(option)
  endif
  
  call VecGetArrayF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    p_min = material_property_array(patch%imat(ghosted_id))%ptr%min_pressure
    p_max = material_property_array(patch%imat(ghosted_id))%ptr%max_pressure
    permfactor_max = material_property_array(patch%imat(ghosted_id))%ptr% &
                     max_permfactor
    if (xx_loc_p(local_id) < p_min) then
      scale = 1
    else 
      if (xx_loc_p(local_id) < p_max) then
        scale = (xx_loc_p(local_id) - p_min)/ &
                (p_max - p_min)*(permfactor_max - 1.d0) + 1.d0
      else
        scale = permfactor_max
      endif
    endif
    !geh: this is a kludge for gfortran.  the code reports errors when 
    !     material_auxvars(ghosted_id)%permeability is used.
    ! Not an issue with Intel
    perm_ptr => material_auxvars(ghosted_id)%permeability
    perm_ptr(perm_xx_index) = perm0_xx_p(local_id)*scale
    perm_ptr(perm_yy_index) = perm0_yy_p(local_id)*scale
    perm_ptr(perm_zz_index) = perm0_zz_p(local_id)*scale
!    material_auxvars(ghosted_id)%permeability(perm_xx_index) = &
!      perm0_xx_p(local_id)*scale
!    material_auxvars(ghosted_id)%permeability(perm_yy_index) = &
!      perm0_yy_p(local_id)*scale
!    material_auxvars(ghosted_id)%permeability(perm_zz_index) = &
!      perm0_zz_p(local_id)*scale
  enddo
  
  call VecRestoreArrayF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)
  call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)
  call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)

  
end subroutine RichardsUpdatePermPatch

! ************************************************************************** !

subroutine RichardsUpdateAuxVars(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the Richards problem
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  type(realization_subsurface_type) :: realization
  
  call RichardsUpdateAuxVarsPatch(realization)

end subroutine RichardsUpdateAuxVars

! ************************************************************************** !

subroutine RichardsUpdateAuxVarsPatch(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the Richards problem
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Logging_module
  use Region_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(region_type), pointer :: region
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_auxvars(:) 
  type(richards_auxvar_type), pointer :: rich_auxvars_bc(:)
  type(richards_auxvar_type), pointer :: rich_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, region_id
  PetscInt :: iphasebc, iphase, i, istart, iend
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof), Pl
  PetscErrorCode :: ierr
  Vec :: phi
  
  call PetscLogEventBegin(logging%event_r_auxvars,ierr);CHKERRQ(ierr)

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  rich_auxvars_ss => patch%aux%Richards%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
    
  call VecGetArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle

    iend = ghosted_id*option%nflowdof
    istart = iend - option%nflowdof + 1

    call RichardsAuxVarCompute(xx_loc_p(istart:iend), &
                               rich_auxvars(ghosted_id), &
                               global_auxvars(ghosted_id), &
                               material_auxvars(ghosted_id), &
                               patch%characteristic_curves_array( &
                                 patch%cc_id(ghosted_id))%ptr, &
                               grid%nG2A(ghosted_id), &
                               PETSC_TRUE,option)   
  enddo

  if (option%flow%inline_surface_flow) then
    region => RegionGetPtrFromList(option%flow%inline_surface_region_name, &
                                   patch%region_list)
     do region_id = 1, patch%aux%InlineSurface%num_aux
        ghosted_id = region%cell_ids(region_id)
        call InlineSurfaceAuxVarCompute(patch%aux%InlineSurface% &
                                          auxvars(region_id), &
                                        global_auxvars(ghosted_id),option)
     enddo
  endif
  
  call PetscLogEventEnd(logging%event_r_auxvars,ierr);CHKERRQ(ierr)

  call PetscLogEventBegin(logging%event_r_auxvars_bc,ierr);CHKERRQ(ierr)

  ! boundary conditions
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      istart = (ghosted_id-1)*option%nflowdof + 1

      if (patch%imat(ghosted_id) <= 0) cycle

      select case(boundary_condition%flow_condition% &
                    itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC, DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
             HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
             HET_SURF_HYDROSTATIC_SEEPAGE_BC, &
             HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
             HET_HYDROSTATIC_CONDUCTANCE_BC, &
             SURFACE_DIRICHLET,SURFACE_SPILLOVER)
          xxbc(1) = boundary_condition% &
                      flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC,UNIT_GRADIENT_BC, &
             SURFACE_ZERO_GRADHEIGHT)
          xxbc(1) = xx_loc_p(istart)
      end select
     
      call RichardsAuxVarCompute(xxbc(1),rich_auxvars_bc(sum_connection), &
                                 global_auxvars_bc(sum_connection), &
                                 material_auxvars(ghosted_id), &
                                 patch%characteristic_curves_array( &
                                   patch%cc_id(ghosted_id))%ptr, &
                                 -grid%nG2A(ghosted_id), &
                                 PETSC_FALSE,option)
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  ! inline surface boundary conditions
  if (option%flow%inline_surface_flow) then  
    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0
    do 
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        if (patch%imat(ghosted_id) <= 0) cycle

        select case(boundary_condition%flow_condition% &
                      itype(RICHARDS_PRESSURE_DOF))
          case(SURFACE_DIRICHLET,SURFACE_SPILLOVER)
            ! bc is specified as a water height, but we need pressure
            Pl = global_auxvars_bc(sum_connection)%pres(1)
            if (Pl < 100.0d0) then
              Pl = (Pl + &
                patch%aux%InlineSurface%auxvars_bc(iconn)%half_cell_height)* &
                patch%aux%InlineSurface%auxvars_bc(iconn)%density* &
                FMWH2O*ABS(option%gravity(3)) &
                + option%flow%reference_pressure
              global_auxvars_bc(sum_connection)%pres(1) = Pl
            endif

            call InlineSurfaceAuxVarCompute(patch%aux%InlineSurface% &
                                              auxvars_bc(iconn), &
                                            global_auxvars_bc(sum_connection), &
                                            option)

          case(SURFACE_ZERO_GRADHEIGHT)
            call InlineSurfaceAuxVarCompute(patch%aux%InlineSurface% &
                                            auxvars_bc(iconn), &
                                            global_auxvars_bc(sum_connection), &
                                            option)
        end select

      enddo
      boundary_condition => boundary_condition%next
    enddo
  endif

  ! source/sinks
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

      call RichardsAuxVarCopy(rich_auxvars(ghosted_id), &
                              rich_auxvars_ss(sum_connection),option)
      call GlobalAuxVarCopy(global_auxvars(ghosted_id), &
                            global_auxvars_ss(sum_connection),option)

    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%Richards%auxvars_up_to_date = PETSC_TRUE
  if (option%flow%inline_surface_flow) then
     patch%aux%InlineSurface%auxvars_up_to_date = PETSC_TRUE
  endif

  call PetscLogEventEnd(logging%event_r_auxvars_bc,ierr);CHKERRQ(ierr)

end subroutine RichardsUpdateAuxVarsPatch

! ************************************************************************** !

subroutine RichardsInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/20/08
  ! 

  use Realization_Subsurface_class
  use Field_module 
  
  implicit none
  


  type(realization_subsurface_type) :: realization

  PetscViewer :: viewer
  PetscErrorCode :: ierr

  type(field_type), pointer :: field


  field => realization%field


  call RichardsUpdateFixedAccum(realization)

  richards_ni_count = 0

  if (realization%option%flow%quasi_3d) &
    call RichardsComputeLateralMassFlux(realization)

!   call PetscViewerASCIIOpen(realization%option%mycomm,'flow_yy.out', &
!                              viewer,ierr)
!    call VecView(field%flow_xx_faces, viewer, ierr)
!    call VecView(field%flow_yy, viewer, ierr)
!
!    call PetscViewerDestroy(viewer,ierr)
!    write(*,*) "Flow_yy" 
!    read(*,*)    

end subroutine RichardsInitializeTimestep

! ************************************************************************** !

subroutine RichardsUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/13/08
  ! 

  use Realization_Subsurface_class
  use Field_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  call RichardsUpdateSolutionPatch(realization)

end subroutine RichardsUpdateSolution

! ************************************************************************** !

subroutine RichardsUpdateSolutionPatch(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/13/08
  ! 

  use Realization_Subsurface_class
    
  implicit none
  
  type(realization_subsurface_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call RichardsUpdateMassBalancePatch(realization)
  endif
  
  richards_ts_count = richards_ts_count + 1
  richards_ts_cut_count = 0
  richards_ni_count = 0

  if (realization%option%flow%update_flow_perm) then
!TODO(geh): this is in the wrong place  
    call RichardsUpdatePermPatch(realization)
  endif

end subroutine RichardsUpdateSolutionPatch

! ************************************************************************** !

subroutine RichardsUpdateFixedAccum(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class

  type(realization_subsurface_type) :: realization
  
  call RichardsUpdateFixedAccumPatch(realization)

end subroutine RichardsUpdateFixedAccum

! ************************************************************************** !

subroutine RichardsUpdateFixedAccumPatch(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Region_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(region_type), pointer :: region
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id
  PetscInt :: numfaces, jface, ghost_face_id, j, region_id
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: accum_p(:)
  PetscReal :: Res(1)
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
    
  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

!  numfaces = 6     ! hex only
!  allocate(sq_faces(numfaces))
!  allocate(faces_pr(numfaces))

  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    call RichardsAuxVarCompute(xx_p(local_id:local_id), &
                   rich_auxvars(ghosted_id),global_auxvars(ghosted_id), &
                   material_auxvars(ghosted_id), &
                   patch%characteristic_curves_array( &
                         patch%cc_id(ghosted_id))%ptr, &
                   grid%nG2A(ghosted_id), &
                   PETSC_TRUE,option)
    call RichardsAccumulation(rich_auxvars(ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              option,accum_p(local_id:local_id))
  enddo

  
  if (option%flow%inline_surface_flow) then
    region => RegionGetPtrFromList(option%flow%inline_surface_region_name, &
                                   patch%region_list)
    do region_id = 1, patch%aux%InlineSurface%num_aux
      local_id = region%cell_ids(region_id)
      ghosted_id = grid%nL2G(local_id)
      call InlineSurfaceAuxVarCompute(patch%aux%InlineSurface% &
                                        auxvars(region_id), &
                                      global_auxvars(ghosted_id),option)
      if (patch%imat(ghosted_id) <= 0) cycle
      call InlineSurfaceAccumulation(patch%aux%InlineSurface% &
                                       auxvars(region_id), &
                                     material_auxvars(ghosted_id),option,Res)
      accum_p(local_id:local_id) = accum_p(local_id:local_id) + Res(1)
    enddo
  endif
   
  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

end subroutine RichardsUpdateFixedAccumPatch

! ************************************************************************** !

subroutine RichardsNumericalJacTest(xx,realization)
  ! 
  ! Computes the a test numerical jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none

  Vec :: xx
  type(realization_subsurface_type) :: realization

  Vec :: xx_pert
  Vec :: res
  Vec :: res_pert
  Mat :: A
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  PetscReal :: derivative, perturbation
  
  PetscReal, pointer :: vec_p(:), vec2_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  
  PetscInt :: idof, idof2, icell

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  
  call VecDuplicate(xx,xx_pert,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res_pert,ierr);CHKERRQ(ierr)
  
  call MatCreate(option%mycomm,A,ierr);CHKERRQ(ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof, &
                   grid%nlmax*option%nflowdof, &
                   ierr);CHKERRQ(ierr)
  call MatSetType(A,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);CHKERRQ(ierr)
    
  call RichardsResidual(PETSC_NULL_SNES,xx,res,realization,ierr)
  call VecGetArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)
  do icell = 1,grid%nlmax
    if (patch%imat(grid%nL2G(icell)) <= 0) cycle
     idof = icell
!    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof 
      call veccopy(xx,xx_pert,ierr);CHKERRQ(ierr)
      call vecgetarrayf90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call vecrestorearrayf90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      call RichardsResidual(PETSC_NULL_SNES,xx_pert,res_pert,realization,ierr)
      call vecgetarrayf90(res_pert,vec_p,ierr);CHKERRQ(ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call matsetvalue(a,idof2-1,idof-1,derivative,insert_values, &
                           ierr);CHKERRQ(ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
!    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call PetscViewerASCIIOpen(option%mycomm,'numerical_jacobian.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(A,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  call MatDestroy(A,ierr);CHKERRQ(ierr)
  
  call VecDestroy(xx_pert,ierr);CHKERRQ(ierr)
  call VecDestroy(res,ierr);CHKERRQ(ierr)
  call VecDestroy(res_pert,ierr);CHKERRQ(ierr)
  
end subroutine RichardsNumericalJacTest

! ************************************************************************** !

subroutine RichardsResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Option_module
  use Logging_module
  use Material_module
  use Material_Aux_class
  use Variables_module
  use Debug_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscViewer :: viewer
  PetscInt :: skip_conn_type
  PetscErrorCode :: ierr

  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string

  call PetscLogEventBegin(logging%event_r_residual,ierr);CHKERRQ(ierr)
  
  field => realization%field
  option => realization%option

  call RichardsResidualPreliminaries(xx,r,realization,ierr)

  skip_conn_type = NO_CONN
  if (option%flow%only_vertical_flow) skip_conn_type = HORZ_CONN

  call RichardsResidualInternalConn(r,realization,skip_conn_type,ierr)
  call RichardsResidualBoundaryConn(r,realization,ierr)
  call RichardsResidualAccumulation(r,realization,ierr)
  call RichardsResidualSourceSink(r,realization,ierr)

  if (realization%debug%vecview_residual) then
    call DebugWriteFilename(realization%debug,string,'Rresidual','', &
                            richards_ts_count,richards_ts_cut_count, &
                            richards_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  if (realization%debug%vecview_solution) then
    call DebugWriteFilename(realization%debug,string,'Rxx','', &
                            richards_ts_count,richards_ts_cut_count, &
                            richards_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

  call PetscLogEventEnd(logging%event_r_residual,ierr);CHKERRQ(ierr)

end subroutine RichardsResidual

! ************************************************************************** !

subroutine RichardsResidualPreliminaries(xx,r,realization,ierr)
  ! 
  ! Perform preliminary work prior to residual computation
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/09/2016
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
  implicit none

  Vec, intent(inout) :: xx
  Vec, intent(inout) :: r
  type(realization_subsurface_type) :: realization

  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscErrorCode :: ierr

  patch => realization%patch
  option => realization%option

  call VecZeroEntries(r, ierr); CHKERRQ(ierr)

  call RichardsUpdateLocalVecs(xx,realization,ierr)

  call RichardsUpdateAuxVarsPatch(realization)

  ! override flags since they will soon be out of date
  patch%aux%Richards%auxvars_up_to_date = PETSC_FALSE 
  ! override flags since they will soon be out of date
  patch%aux%Richards%auxvars_cell_pressures_up_to_date = PETSC_FALSE 

  if (option%compute_mass_balance_new) &
    call RichardsZeroMassBalDeltaPatch(realization)

end subroutine RichardsResidualPreliminaries

! ************************************************************************** !

subroutine RichardsUpdateLocalVecs(xx,realization,ierr)
  ! 
  ! Updates local vectors needed for residual computation
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/09/2016
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Option_module
  use Logging_module
  use Material_module
  use Material_Aux_class
  use Variables_module
  use Debug_module

  implicit none

  Vec :: xx
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string

  field => realization%field
  discretization => realization%discretization
  option => realization%option

  ! Communication -----------------------------------------
  ! These 3 must be called before RichardsUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)

  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)
  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)
  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)

end subroutine RichardsUpdateLocalVecs

! ************************************************************************** !

subroutine RichardsResidualInternalConn(r,realization,skip_conn_type,ierr)
  ! 
  ! Computes the interior flux terms of the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Debug_module
  use Region_module
  
  implicit none

  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscInt :: skip_conn_type
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(region_type), pointer :: region
  type(material_parameter_type), pointer :: material_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(inlinesurface_auxvar_type), pointer :: insurf_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: istart
  PetscInt :: local_id_up
  PetscInt :: local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: region_id_up, region_id_dn
  PetscInt :: icc_up
  PetscInt :: icc_dn
  PetscInt :: iconn
  PetscInt :: sum_connection

  PetscReal :: Res(realization%option%nflowdof)
  PetscReal :: v_darcy
  PetscReal, pointer :: r_p(:)

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  material_parameter => patch%aux%Material%material_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  if (option%flow%inline_surface_flow) then
    insurf_auxvars => patch%aux%InlineSurface%auxvars
  endif
   
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

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

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      if (.not.(skip_conn_type == NO_CONN)) then
        if (skip_conn(cur_connection_set%dist(1:3,iconn), skip_conn_type)) cycle
      endif

      icc_up = patch%cc_id(ghosted_id_up)
      icc_dn = patch%cc_id(ghosted_id_dn)

      call RichardsFlux(rich_auxvars(ghosted_id_up), &
                      global_auxvars(ghosted_id_up), &
                      material_auxvars(ghosted_id_up), &
                      rich_auxvars(ghosted_id_dn), &
                      global_auxvars(ghosted_id_dn), &
                      material_auxvars(ghosted_id_dn), &
                      cur_connection_set%area(iconn), &
                      cur_connection_set%dist(:,iconn), &
                      option,v_darcy,Res)

      patch%internal_velocities(1,sum_connection) = v_darcy
      if (associated(patch%internal_flow_fluxes)) then
        patch%internal_flow_fluxes(1,sum_connection) = Res(1)
      endif
      if (local_id_up>0) then
        istart = (local_id_up-1)*option%nflowdof + 1
        r_p(istart) = r_p(istart) + Res(1)
      endif

      if (local_id_dn>0) then
        istart = (local_id_dn-1)*option%nflowdof + 1
        r_p(istart) = r_p(istart) - Res(1)
      endif

    enddo

    cur_connection_set => cur_connection_set%next
  enddo

  ! Regional Interior Flux Terms -----------------------------------
  if (option%flow%inline_surface_flow) then
    connection_set_list => grid%reg_internal_connection_set_list
    cur_connection_set => connection_set_list%first
    region => RegionGetPtrFromList(option%flow%inline_surface_region_name, &
                                   patch%region_list)
  else
    nullify(cur_connection_set)
  endif
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      
      region_id_up = cur_connection_set%id_up(iconn)
      region_id_dn = cur_connection_set%id_dn(iconn)
      
      ghosted_id_up = region%cell_ids(region_id_up)
      ghosted_id_dn = region%cell_ids(region_id_dn)
      
      local_id_up = grid%nG2L(ghosted_id_up) 
      local_id_dn = grid%nG2L(ghosted_id_dn) 
      
      if (patch%imat(ghosted_id_up) <= 0 .or.  &
           patch%imat(ghosted_id_dn) <= 0) cycle
      
      call InlineSurfaceFlux(insurf_auxvars(region_id_up), &
                             insurf_auxvars(region_id_dn), &
                             cur_connection_set%area(iconn), &
                             cur_connection_set%dist(:,iconn), &
                             Res)
      
      if (local_id_up>0) then
        istart = (local_id_up-1)*option%nflowdof + 1
        r_p(istart) = r_p(istart) + Res(1)
      endif
      
      if (local_id_dn>0) then
        istart = (local_id_dn-1)*option%nflowdof + 1
        r_p(istart) = r_p(istart) - Res(1)
      endif
      
    enddo

    cur_connection_set => cur_connection_set%next
  enddo
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)

end subroutine RichardsResidualInternalConn

! ************************************************************************** !

subroutine RichardsResidualBoundaryConn(r,realization,ierr)
  ! 
  ! Computes the boundary flux terms of the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Debug_module
  use Region_module
  
  implicit none

  Vec :: r
  type(realization_subsurface_type) :: realization

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(region_type), pointer :: region
  type(coupler_type), pointer :: boundary_condition
  type(material_parameter_type), pointer :: material_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:), rich_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: region_id, i
  PetscInt :: istart
  PetscInt :: icc_up
  PetscInt :: icc_dn
  PetscInt :: iconn
  PetscInt :: sum_connection

  PetscReal :: Res(realization%option%nflowdof), v_darcy
  PetscReal, pointer :: r_p(:)

  PetscErrorCode :: ierr

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  material_parameter => patch%aux%Material%material_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

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

      if (patch%imat(ghosted_id) <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icc_dn = patch%cc_id(ghosted_id)

      call RichardsBCFlux(boundary_condition%flow_condition%itype, &
                       boundary_condition%flow_aux_real_var(:,iconn), &
                       rich_auxvars_bc(sum_connection), &
                       global_auxvars_bc(sum_connection), &
                       rich_auxvars(ghosted_id), &
                       global_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       cur_connection_set%area(iconn), &
                       cur_connection_set%dist(:,iconn), &
                       option, &
                       v_darcy,Res)
      patch%boundary_velocities(1,sum_connection) = v_darcy
      if (associated(patch%boundary_flow_fluxes)) then
        patch%boundary_flow_fluxes(1,sum_connection) = Res(1)
      endif

      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)%mass_balance_delta(1,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(1,1) - Res(1)
      endif

      istart = (local_id-1)*option%nflowdof + 1
      r_p(istart)= r_p(istart) - Res(1)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Inline surface BCs
  if (option%flow%inline_surface_flow) then
    sum_connection = 0
    boundary_condition => patch%boundary_condition_list%first
    region => RegionGetPtrFromList(option%flow%inline_surface_region_name, &
                                   patch%region_list)    
    do
      if (.not.associated(boundary_condition)) exit
      if (boundary_condition%flow_condition%pressure%itype == &
                                                 SURFACE_DIRICHLET       .or. &
          boundary_condition%flow_condition%pressure%itype == &
                                                 SURFACE_ZERO_GRADHEIGHT .or. &
          boundary_condition%flow_condition%pressure%itype == &
                                                 SURFACE_SPILLOVER ) then
        do iconn = 1,boundary_condition%connection_set%num_connections
          sum_connection = sum_connection + 1
          
          local_id   = boundary_condition%connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          
          region_id  = -1
          do i = 1,region%num_cells
            if (region%cell_ids(i) == local_id) then
              region_id = i
              exit
            endif
          enddo
          
          Res(1) = 0.0d0
          call InlineSurfaceBCFlux(boundary_condition%flow_condition%itype, &
               patch%aux%InlineSurface%auxvars_bc(sum_connection), &
               patch%aux%InlineSurface%auxvars(region_id), &
               boundary_condition%connection_set%area(iconn), &
               boundary_condition%connection_set%dist(:,iconn), &
               Res)
          
          istart = (local_id-1)*option%nflowdof + 1
          r_p(istart)= r_p(istart) - Res(1)
          
        enddo
      endif
      boundary_condition => boundary_condition%next
    enddo
     
  endif
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)

end subroutine RichardsResidualBoundaryConn

! ************************************************************************** !

subroutine RichardsResidualSourceSink(r,realization,ierr)
  ! 
  ! Computes the accumulation and source/sink terms of
  ! the residual equation on a single patch
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Field_module
  use Debug_module
  
  implicit none

  Vec :: r
  type(realization_subsurface_type) :: realization

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: rich_auxvars(:), rich_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: i
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: qsrc, qsrc_mol
  PetscReal :: Res(realization%option%nflowdof)
  PetscReal, pointer :: r_p(:), accum_p(:)
  PetscReal, pointer :: mmsrc(:)
  PetscReal, allocatable :: msrc(:)
  PetscReal :: well_status
  PetscReal :: well_factor
  PetscReal :: pressure_bh
  PetscReal :: pressure_max
  PetscReal :: pressure_min
  PetscReal :: well_inj_water
  PetscReal :: Dq, dphi, v_darcy, ukvr

  Mat, parameter :: null_mat = tMat(0)

  PetscErrorCode :: ierr

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_ss => patch%aux%Richards%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars

  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

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

      if (source_sink%flow_condition%itype(1)/=HET_VOL_RATE_SS .and. &
          source_sink%flow_condition%itype(1)/=HET_MASS_RATE_SS .and. &
          source_sink%flow_condition%itype(1)/=WELL_SS) &
        qsrc = source_sink%flow_condition%rate%dataset%rarray(1)

      select case(source_sink%flow_condition%itype(1))
        case(MASS_RATE_SS)
          qsrc_mol = qsrc/FMWH2O ! kg/sec -> kmol/sec
        case(SCALED_MASS_RATE_SS)
          qsrc_mol = qsrc/FMWH2O* & ! kg/sec -> kmol/sec
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc_mol = qsrc*global_auxvars(ghosted_id)%den(1) ! den = kmol/m^3
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc_mol = qsrc*global_auxvars(ghosted_id)%den(1)* & ! den = kmol/m^3
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_VOL_RATE_SS)
          ! qsrc1 = m^3/sec
                                 ! flow = m^3/s
          qsrc_mol = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)* & 
                     global_auxvars(ghosted_id)%den(1)    ! den  = kmol/m^3
        case(HET_MASS_RATE_SS)
                                                           ! kg/sec -> kmol/sec
          qsrc_mol = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)/FMWH2O 
      
        case(WELL_SS) ! production well, SK 12/19/13
          ! if node pessure is lower than the given extraction pressure, 
          ! shut it down
          !  well parameter explanation
          !   1. well status. 1 injection; -1 production; 0 shut in!
          !   2. well factor [m^3],  the effective permeability [m^2/s]
          !   3. bottomhole pressure:  [Pa]
          !   4. max pressure: [Pa]
          !   5. min pressure: [Pa]   
          mmsrc => source_sink%flow_condition%well%dataset%rarray

          well_status = mmsrc(1)
          well_factor = mmsrc(2)
          pressure_bh = mmsrc(3)
          pressure_max = mmsrc(4)
          pressure_min = mmsrc(5)
    
          ! production well (well status = -1)
          if (dabs(well_status + 1.D0) < 1.D-1) then
            if (global_auxvars(ghosted_id)%pres(1) > pressure_min) then
              Dq = well_factor 
              dphi = global_auxvars(ghosted_id)%pres(1) - pressure_bh
              if (dphi >= 0.D0) then ! outflow only
                ukvr = rich_auxvars(ghosted_id)%kvr
                if (ukvr < 1.e-20) ukvr = 0.D0
                v_darcy = 0.D0
                if (ukvr*Dq > floweps) then
                  v_darcy = Dq * ukvr * dphi
                  ! store volumetric rate for ss_fluid_fluxes()
                  qsrc_mol = -1.d0*v_darcy*global_auxvars(ghosted_id)%den(1)
                endif
              endif
            endif
          endif 
      end select

      if (option%compute_mass_balance_new) then
        ! need to added global auxvar for src/sink
        global_auxvars_ss(sum_connection)%mass_balance_delta(1,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1,1) - &
          qsrc_mol
      endif

      istart = (local_id-1)*option%nflowdof + 1
      r_p(istart) = r_p(istart) - qsrc_mol

      if (associated(patch%ss_flow_vol_fluxes)) then
        ! fluid flux [m^3/sec] = qsrc_mol [kmol/sec] / den [kmol/m^3]
        patch%ss_flow_vol_fluxes(1,sum_connection) = qsrc_mol / &
                                           global_auxvars(ghosted_id)%den(1)
      endif
      if (associated(patch%ss_flow_fluxes)) then
        ! fluid flux [m^3/sec] = qsrc_mol [kmol/sec] / den [kmol/m^3]
        patch%ss_flow_fluxes(1,sum_connection) = qsrc_mol
      endif
    enddo
    source_sink => source_sink%next
  enddo

  if (patch%aux%Richards%inactive_cells_exist) then
    do i=1,patch%aux%Richards%matrix_zeroing%n_inactive_rows
      r_p(patch%aux%Richards%matrix_zeroing%inactive_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  call RichardsSSSandbox(r,null_mat,PETSC_FALSE,grid,material_auxvars, &
                         global_auxvars,rich_auxvars,option)
  
  ! Mass Transfer
  if (field%flow_mass_transfer /= PETSC_NULL_VEC) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    call VecAXPY(r,-1.d0,field%flow_mass_transfer,ierr);CHKERRQ(ierr)
  endif

end subroutine RichardsResidualSourceSink

! ************************************************************************** !

subroutine RichardsResidualAccumulation(r,realization,ierr)
  ! 
  ! Computes the accumulation terms of the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Region_module

  implicit none

  Vec :: r
  type(realization_subsurface_type) :: realization

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(region_type), pointer :: region
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(inlinesurface_auxvar_type), pointer :: inlinesurface_auxvars(:)
  
  PetscInt :: local_id, ghosted_id, region_id
  PetscInt :: istart

  PetscReal, pointer :: r_p(:), accum_p(:), accum2_p(:)
  PetscReal :: Res(realization%option%nflowdof)

  PetscErrorCode :: ierr

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  if (option%flow%inline_surface_flow) then
    region => RegionGetPtrFromList(option%flow%inline_surface_region_name, &
         patch%region_list)
    inlinesurface_auxvars => patch%aux%InlineSurface%auxvars
  endif

  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum2, accum2_p, ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  if (.not.option%flow%steady_state) then
    r_p(1:grid%nlmax) = r_p(1:grid%nlmax) - accum_p(1:grid%nlmax)

    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      call RichardsAccumulation(rich_auxvars(ghosted_id), &
           global_auxvars(ghosted_id), &
           material_auxvars(ghosted_id), &
           option,Res)
      istart = (local_id-1)*option%nflowdof + 1
      r_p(istart) = r_p(istart) + Res(1)
      accum2_p(istart) = Res(1)
    enddo

    if (option%flow%inline_surface_flow) then
      ! Loop through cells in the defined region
      do region_id = 1, region%num_cells 
        local_id = region%cell_ids(region_id)
        ghosted_id = grid%nL2G(local_id)         
        if (patch%imat(ghosted_id) <= 0) cycle
        call InlineSurfaceAccumulation(inlinesurface_auxvars(region_id), &
             material_auxvars(ghosted_id),option,Res)
        istart = (local_id-1)*option%nflowdof + 1
        r_p(istart) = r_p(istart) + Res(1)
        accum2_p(istart) = accum2_p(istart) + Res(1)
      enddo
    endif

  endif

  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum2, accum2_p, ierr);CHKERRQ(ierr)

end subroutine RichardsResidualAccumulation

! ************************************************************************** !

subroutine RichardsJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Logging_module
  use Debug_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type_A, mat_type_B
  PetscViewer :: viewer
  type(grid_type),  pointer :: grid
  type(option_type), pointer :: option
  PetscReal :: norm
  character(len=MAXSTRINGLENGTH) :: string

  call PetscLogEventBegin(logging%event_r_jacobian,ierr);CHKERRQ(ierr)

  option => realization%option

  call MatGetType(A,mat_type_A,ierr);CHKERRQ(ierr)
  call MatGetType(B,mat_type_B,ierr);CHKERRQ(ierr)
  if (mat_type_A == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  call RichardsJacobianInternalConn(J,realization,ierr)
  call RichardsJacobianBoundaryConn(J,realization,ierr)
  call RichardsJacobianAccumulation(J,realization,ierr)
  call RichardsJacobianSourceSink(J,realization,ierr)

  if (A /= B .and. mat_type_A /= MATMFFD) then
    ! If the Jacobian and preconditioner matrices are different (and not 
    ! because we are using a "matrix free" Jacobian), then we need to 
    ! copy the computed Jacobian into the preconditioner matrix.
    call MatConvert(J,mat_type_B,MAT_REUSE_MATRIX,B,ierr);CHKERRQ(ierr);
  endif

  if (realization%debug%matview_Jacobian) then
    call DebugWriteFilename(realization%debug,string,'Rjacobian','', &
                            richards_ts_count,richards_ts_cut_count, &
                            richards_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call PrintMsg(option)
    call MatNorm(J,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call PrintMsg(option)
    call MatNorm(J,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call PrintMsg(option)
  endif

#if 0
    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_dxx.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(realization%field%flow_dxx,viewer,ierr);CHKERRQ(ierr)

    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
 

    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_yy.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(realization%field%flow_yy,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call PetscLogEventEnd(logging%event_r_jacobian,ierr);CHKERRQ(ierr)
!  call PrintErrMsg(option)

  richards_ni_count = richards_ni_count + 1
  
end subroutine RichardsJacobian

! ************************************************************************** !

subroutine RichardsJacobianInternalConn(A,realization,ierr)
  ! 
  ! Computes the interior flux terms of the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
       
  use Connection_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_class
  use Region_module
  
  implicit none

  Mat, intent(inout) :: A
  type(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr

  PetscInt :: icc_up,icc_dn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: region_id_up, region_id_dn
  PetscInt :: istart_up, istart_dn, istart

  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)

  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(region_type), pointer :: region
  type(material_parameter_type), pointer :: material_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(inlinesurface_auxvar_type), pointer :: insurf_auxvars(:)
  
  character(len=MAXSTRINGLENGTH) :: string

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  if (option%flow%inline_surface_flow) then
    insurf_auxvars => patch%aux%InlineSurface%auxvars
  endif
   
#ifdef BUFFER_MATRIX
  if (option%use_matrix_buffer) then
    if (associated(patch%aux%Richards%matrix_buffer)) then
      call MatrixBufferZero(patch%aux%Richards%matrix_buffer)
    else
      patch%aux%Richards%matrix_buffer => MatrixBufferCreate()
      call MatrixBufferInit(A,patch%aux%Richards%matrix_buffer,grid)
    endif
  endif
#endif

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

      if (patch%imat(ghosted_id_up) <= 0 .or. &
          patch%imat(ghosted_id_dn) <= 0) cycle

      if (option%flow%only_vertical_flow) then
        !geh: place second conditional within first to avoid excessive
        !     dot products when .not. option%flow%only_vertical_flow
        if (abs(dot_product(cur_connection_set%dist(1:3,iconn),unit_z)) < &
            0.99d0) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

      icc_up = patch%cc_id(ghosted_id_up)
      icc_dn = patch%cc_id(ghosted_id_dn)

      call RichardsFluxDerivative(rich_auxvars(ghosted_id_up), &
                     global_auxvars(ghosted_id_up), &
                     material_auxvars(ghosted_id_up), &
                     rich_auxvars(ghosted_id_dn), &
                     global_auxvars(ghosted_id_dn), &
                     material_auxvars(ghosted_id_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(-1:3,iconn),&
                     option,&
                     patch%characteristic_curves_array(icc_up)%ptr, &
                     patch%characteristic_curves_array(icc_dn)%ptr, &
                     Jup,Jdn)

      if (local_id_up > 0) then

#ifdef BUFFER_MATRIX
        if (option%use_matrix_buffer) then
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer, &
                               ghosted_id_up,ghosted_id_up,Jup(1,1))
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer, &
                               ghosted_id_up,ghosted_id_dn,Jdn(1,1))
        else
#endif
          istart_up = (ghosted_id_up-1)*option%nflowdof + 1
          istart_dn = (ghosted_id_dn-1)*option%nflowdof + 1

          call MatSetValuesLocal(A,1,istart_up-1,1,istart_up-1, &
                                        Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValuesLocal(A,1,istart_up-1,1,istart_dn-1, &
                                        Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
#ifdef BUFFER_MATRIX
        endif
#endif
      endif

      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
#ifdef BUFFER_MATRIX
        if (option%use_matrix_buffer) then
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer, &
                               ghosted_id_dn,ghosted_id_dn,Jdn(1,1))
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer, &
                               ghosted_id_dn,ghosted_id_up,Jup(1,1))
        else
#endif
          istart_up = (ghosted_id_up-1)*option%nflowdof + 1
          istart_dn = (ghosted_id_dn-1)*option%nflowdof + 1

          call MatSetValuesLocal(A,1,istart_dn-1,1,istart_dn-1, &
                                        Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValuesLocal(A,1,istart_dn-1,1,istart_up-1, &
                                        Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
#ifdef BUFFER_MATRIX
        endif
#endif
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Regional Interior Flux Terms -----------------------------------
  if (option%flow%inline_surface_flow) then
    connection_set_list => grid%reg_internal_connection_set_list
    cur_connection_set => connection_set_list%first
    region => RegionGetPtrFromList(option%flow%inline_surface_region_name, &
         patch%region_list)
  else
    nullify(cur_connection_set)
  endif
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      region_id_up = cur_connection_set%id_up(iconn)
      region_id_dn = cur_connection_set%id_dn(iconn)

      ghosted_id_up = region%cell_ids(region_id_up)
      ghosted_id_dn = region%cell_ids(region_id_dn)

      local_id_up = grid%nG2L(ghosted_id_up) 
      local_id_dn = grid%nG2L(ghosted_id_dn) 

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
           patch%imat(ghosted_id_dn) <= 0) cycle

      call InlineSurfaceFluxJac(insurf_auxvars(region_id_up), &
           insurf_auxvars(region_id_dn), &
           cur_connection_set%area(iconn), &
           cur_connection_set%dist(:,iconn), &
           option, Jup,Jdn)

      if (local_id_up>0) then
        istart_up = (ghosted_id_up-1)*option%nflowdof + 1
        istart_dn = (ghosted_id_dn-1)*option%nflowdof + 1
        call MatSetValuesLocal(A,1,istart_up-1,1,istart_up-1, &
             Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesLocal(A,1,istart_up-1,1,istart_dn-1, &
             Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif

      if (local_id_dn>0) then
        istart_up = (ghosted_id_up-1)*option%nflowdof + 1
        istart_dn = (ghosted_id_dn-1)*option%nflowdof + 1
        call MatSetValuesLocal(A,1,istart_dn-1,1,istart_dn-1, &
             -Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesLocal(A,1,istart_dn-1,1,istart_up-1, &
             -Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Rjacobian_flux','', &
                            richards_ts_count,richards_ts_cut_count, &
                            richards_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

end subroutine RichardsJacobianInternalConn

! ************************************************************************** !

subroutine RichardsJacobianBoundaryConn(A,realization,ierr)
  ! 
  ! Computes the boundary flux terms of the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_class
  use Region_module
  
  implicit none

  Mat, intent(inout) :: A
  type(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr

  PetscInt :: icc_up,icc_dn
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, region_id, i
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: istart_up, istart_dn, istart
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(region_type), pointer :: region
  type(material_parameter_type), pointer :: material_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:), rich_auxvars_bc(:) 
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  
  character(len=MAXSTRINGLENGTH) :: string

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  
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

      if (patch%imat(ghosted_id) <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icc_dn = patch%cc_id(ghosted_id) 

      call RichardsBCFluxDerivative(boundary_condition%flow_condition%itype, &
                     boundary_condition%flow_aux_real_var(:,iconn), &
                     rich_auxvars_bc(sum_connection), &
                     global_auxvars_bc(sum_connection), &
                     rich_auxvars(ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     option, &
                     patch%characteristic_curves_array(icc_dn)%ptr, &
                     Jdn)
      Jdn = -Jdn

#ifdef BUFFER_MATRIX
      if (option%use_matrix_buffer) then
        call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
                             ghosted_id,Jdn(1,1))
      else
#endif
        istart = (ghosted_id-1)*option%nflowdof + 1

        call MatSetValuesLocal(A,1,istart-1,1,istart-1,Jdn, &
                               ADD_VALUES,ierr);CHKERRQ(ierr)
#ifdef BUFFER_MATRIX
      endif
#endif
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  ! Inline surface BCs
  if (option%flow%inline_surface_flow) then
    sum_connection = 0
    boundary_condition => patch%boundary_condition_list%first
    region => RegionGetPtrFromList(option%flow%inline_surface_region_name, &
         patch%region_list)
    do
      if (.not.associated(boundary_condition)) exit
      if (boundary_condition%flow_condition%pressure%itype == &
                                                SURFACE_DIRICHLET .or. &
          boundary_condition%flow_condition%pressure%itype == &
                                                SURFACE_ZERO_GRADHEIGHT .or. &
          boundary_condition%flow_condition%pressure%itype == &
                                                SURFACE_SPILLOVER ) then
        do iconn = 1,boundary_condition%connection_set%num_connections
          sum_connection = sum_connection + 1

          local_id   = boundary_condition%connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)

          region_id  = -1
          do i = 1,region%num_cells
            if (region%cell_ids(i) == local_id) then
              region_id = i
              exit
            endif
          enddo

          Jup = 0.0d0
          Jdn = 0.0d0
          call InlineSurfaceBCFluxJac(boundary_condition%flow_condition%itype, &
               patch%aux%InlineSurface%auxvars_bc(sum_connection), &
               patch%aux%InlineSurface%auxvars(region_id), &
               boundary_condition%connection_set%area(iconn), &
               boundary_condition%connection_set%dist(:,iconn), &
               option,Jdn)
          Jdn = -Jdn

          istart = (ghosted_id-1)*option%nflowdof + 1
          call MatSetValuesLocal(A,1,istart-1,1,istart-1,Jdn, &
               ADD_VALUES,ierr);CHKERRQ(ierr)

        enddo
      endif
      boundary_condition => boundary_condition%next
    enddo

  endif

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Rjacobian_bcflux','', &
                            richards_ts_count,richards_ts_cut_count, &
                            richards_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)    
  endif
  
end subroutine RichardsJacobianBoundaryConn

! ************************************************************************** !

subroutine RichardsJacobianAccumulation(A,realization,ierr)
  ! 
  ! Computes the accumulation terms of the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Debug_module
  use Region_module
  
  implicit none

  Mat, intent(inout) :: A
  type(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr

  PetscInt :: local_id, ghosted_id, region_id
  PetscInt :: istart

  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof)

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(region_type), pointer :: region
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(inlinesurface_auxvar_type), pointer :: inlinesurface_auxvars(:)
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  if (option%flow%inline_surface_flow) then
    region => &
      RegionGetPtrFromList(option%flow%inline_surface_region_name, &
                           patch%region_list)
    inlinesurface_auxvars => patch%aux%InlineSurface%auxvars
  endif

  if (.not.option%flow%steady_state) then

    ! Accumulation terms ------------------------------------
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      call RichardsAccumDerivative(rich_auxvars(ghosted_id), &
           global_auxvars(ghosted_id), &
           material_auxvars(ghosted_id), &
           option, &
           patch%characteristic_curves_array( &
           patch%cc_id(ghosted_id))%ptr, &
           Jup)

#ifdef BUFFER_MATRIX
      if (option%use_matrix_buffer) then
        call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
             ghosted_id,Jup(1,1))
      else
#endif
        istart = (ghosted_id-1)*option%nflowdof + 1

        call MatSetValuesLocal(A,1,istart-1,1,istart-1,Jup, &
             ADD_VALUES,ierr);CHKERRQ(ierr)
#ifdef BUFFER_MATRIX
      endif
#endif
    enddo

    if (option%flow%inline_surface_flow) then
      do region_id = 1, region%num_cells
        local_id = region%cell_ids(region_id)
        ghosted_id = grid%nL2G(local_id)         
        if (patch%imat(ghosted_id) <= 0) cycle
        call InlineSurfaceAccumulationJac(inlinesurface_auxvars(region_id), &
             material_auxvars(ghosted_id),option,Jup)
        istart = (ghosted_id-1)*option%nflowdof + 1
        call MatSetValuesLocal(A,1,istart-1,1,istart-1,Jup, &
             ADD_VALUES,ierr);CHKERRQ(ierr)
      enddo
    endif

  endif

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Rjacobian_accum','', &
                            richards_ts_count,richards_ts_cut_count, &
                            richards_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

end subroutine RichardsJacobianAccumulation

! ************************************************************************** !

subroutine RichardsJacobianSourceSink(A,realization,ierr)
  ! 
  ! Computes the accumulation and source/sink terms of
  ! the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
    
  implicit none

  Mat, intent(inout) :: A
  type(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr

  PetscReal :: qsrc
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: flow_pc
  PetscViewer :: viewer
  PetscReal, pointer :: mmsrc(:)
  PetscReal :: well_status
  PetscReal :: well_factor
  PetscReal :: pressure_bh
  PetscReal :: pressure_max
  PetscReal :: pressure_min
  PetscReal :: ukvr, Dq, dphi, v_darcy
  Vec, parameter :: null_vec = tVec(0)
  character(len=MAXSTRINGLENGTH) :: string

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first 
  do 
    if (.not.associated(source_sink)) exit
    
    if (source_sink%flow_condition%itype(1)/=HET_VOL_RATE_SS.and. &
       source_sink%flow_condition%itype(1)/=HET_MASS_RATE_SS .and. &
       source_sink%flow_condition%itype(1)/=WELL_SS) &
      qsrc = source_sink%flow_condition%rate%dataset%rarray(1)

    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle
      
      Jup = 0.d0
      select case(source_sink%flow_condition%itype(1))
        case(MASS_RATE_SS,SCALED_MASS_RATE_SS,HET_MASS_RATE_SS)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          Jup(1,1) = -qsrc*rich_auxvars(ghosted_id)%dden_dp*FMWH2O
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          Jup(1,1) = -qsrc*rich_auxvars(ghosted_id)%dden_dp*FMWH2O* &
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_VOL_RATE_SS)
          Jup(1,1) = -source_sink%flow_aux_real_var(ONE_INTEGER,iconn)* &
                    rich_auxvars(ghosted_id)%dden_dp*FMWH2O
        case(WELL_SS) ! production well, SK 12/19/13
          ! if node pessure is lower than the given extraction pressure, 
          ! shut it down
          !  well parameter explanation
          !   1. well status. 1 injection; -1 production; 0 shut in!
          !   2. well factor [m^3],  the effective permeability [m^2/s]
          !   3. bottomhole pressure:  [Pa]
          !   4. max pressure: [Pa]
          !   5. min pressure: [Pa]
          mmsrc => source_sink%flow_condition%well%dataset%rarray

          well_status = mmsrc(1)
          well_factor = mmsrc(2)
          pressure_bh = mmsrc(3)
          pressure_max = mmsrc(4)
          pressure_min = mmsrc(5)
    
          ! production well (well status = -1)
          if (dabs(well_status + 1.D0) < 1.D-1) then
            if (global_auxvars(ghosted_id)%pres(1) > pressure_min) then
              Dq = well_factor 
              dphi = global_auxvars(ghosted_id)%pres(1) - pressure_bh
              if (dphi >= 0.D0) then ! outflow only
                ukvr = rich_auxvars(ghosted_id)%kvr
                if (ukvr < 1.e-20) ukvr = 0.D0
                v_darcy = 0.D0
                if (ukvr*Dq > floweps) then
                  v_darcy = Dq * ukvr * dphi
                  ! store volumetric rate for ss_fluid_fluxes()
                  Jup(1,1) = -1.d0*(-Dq*rich_auxvars(ghosted_id)%dkvr_dp*dphi* &
                             global_auxvars(ghosted_id)%den(1) &
                             -Dq*ukvr*1.d0*global_auxvars(ghosted_id)%den(1) &
                             -Dq*ukvr*dphi*rich_auxvars(ghosted_id)%dden_dp)
                endif
              endif
            endif
          endif 
      end select
#ifdef BUFFER_MATRIX
      if (option%use_matrix_buffer) then
        call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
                             ghosted_id,Jup(1,1))
      else
#endif
        istart = (ghosted_id-1)*option%nflowdof + 1

        call MatSetValuesLocal(A,1,istart-1,1,istart-1,Jup,ADD_VALUES, &
                               ierr);CHKERRQ(ierr)
#ifdef BUFFER_MATRIX
      endif
#endif
    enddo
    source_sink => source_sink%next
  enddo

  call RichardsSSSandbox(null_vec,A,PETSC_TRUE,grid,material_auxvars, &
                         global_auxvars,rich_auxvars,option)

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Rjacobian_srcsink','', &
                            richards_ts_count,richards_ts_cut_count, &
                            richards_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  
#ifdef BUFFER_MATRIX
  if (option%use_matrix_buffer) then
    if (patch%aux%Richards%inactive_cells_exist) then
      call MatrixBufferZeroRows(patch%aux%Richards%matrix_buffer, &
                                patch%aux%Richards%matrix_zeroing% &
                                  n_inactive_rows, &
                                patch%aux%Richards%matrix_zeroing% &
                                  inactive_rows_local_ghosted)
    endif
    call MatrixBufferSetValues(A,patch%aux%Richards%matrix_buffer)
  endif
#endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

! zero out isothermal and inactive cells
#ifdef BUFFER_MATRIX
  if (.not.option%use_matrix_buffer) then
#endif
    if (patch%aux%Richards%inactive_cells_exist) then
      qsrc = 1.d0 ! solely a temporary variable in this conditional
      call MatZeroRowsLocal(A,patch%aux%Richards%matrix_zeroing% &
                              n_inactive_rows, &
                            patch%aux%Richards%matrix_zeroing% &
                              inactive_rows_local_ghosted, &
                            qsrc,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                            ierr);CHKERRQ(ierr)
    endif
#ifdef BUFFER_MATRIX
  endif
#endif

end subroutine RichardsJacobianSourceSink

! ************************************************************************** !

subroutine RichardsMaxChange(realization,dpmax)
  ! 
  ! Computes the maximum change in the solution vector
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/08
  ! 

  use Realization_Base_class
  use Option_module
  use Field_module
  
  implicit none
  
  class(realization_base_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  PetscReal :: dpmax
  
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  option => realization%option
  field => realization%field

  dpmax = 0.d0
  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy, &
                ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY, &
                     dpmax,ierr);CHKERRQ(ierr)

end subroutine RichardsMaxChange

! ************************************************************************** !

subroutine RichardsSetPlotVariables(list)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 
  
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
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

    name = 'Liquid Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                LIQUID_SATURATION)
                                
  endif
  
end subroutine RichardsSetPlotVariables

! ************************************************************************** !

subroutine RichardsPrintAuxVars(richards_auxvar,global_auxvar,cell_id)
  ! 
  ! Prints out the contents of an auxvar
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/21/12
  ! 

  use Global_Aux_module

  implicit none

  type(richards_auxvar_type) :: richards_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: cell_id

  print *, '      cell: ', cell_id
  print *, '  pressure: ', global_auxvar%pres(1)
  print *, 'saturation: ', global_auxvar%sat(1)
  print *, '   density: ', global_auxvar%den_kg(1)
  print *, '        pc: ', richards_auxvar%pc
  print *, '       kvr: ', richards_auxvar%kvr
  print *, '   dkvr_dp: ', richards_auxvar%dkvr_dp
  print *, '   dsat_dp: ', richards_auxvar%dsat_dp
  print *, '   dden_dp: ', richards_auxvar%dden_dp

end subroutine RichardsPrintAuxVars

! ************************************************************************** !

subroutine RichardsSSSandbox(residual,Jacobian,compute_derivative, &
                             grid,material_auxvars,global_auxvars, &
                             rich_auxvars,option)
  ! 
  ! Evaluates source/sink term storing residual and/or Jacobian
  ! 
  ! Author: Guoping Tang
  ! Date: 06/03/14
  ! 
  ! Modified by: Ayman Alzraiee on 04/05/2016 

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module
  use Grid_module
  use Material_Aux_class, only: material_auxvar_type
  use SrcSink_Sandbox_module
  use SrcSink_Sandbox_Base_class
  
  implicit none

  PetscBool :: compute_derivative
  Vec :: residual
  Mat :: Jacobian
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(grid_type) :: grid
  type(option_type) :: option
  
  PetscReal, pointer :: r_p(:)
  PetscReal :: res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  PetscInt :: local_id, ghosted_id, istart, iend
  PetscReal :: aux_real(10)
  PetscErrorCode :: ierr
  
  if (.not.compute_derivative) then
    call VecGetArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif
  
  cur_srcsink => ss_sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
    aux_real = 0.d0
    local_id = cur_srcsink%local_cell_id
    ghosted_id = grid%nL2G(local_id)
    res = 0.d0
    Jac = 0.d0
    call RichardsSSSandboxLoadAuxReal(cur_srcsink,aux_real, &
                                      global_auxvars(ghosted_id), &
                                      rich_auxvars(ghosted_id),option)
    call cur_srcsink%Evaluate(res,Jac,PETSC_FALSE, &
                              material_auxvars(ghosted_id), &
                              aux_real,option)
    if (compute_derivative) then
      call RichardsSSSandboxLoadAuxReal(cur_srcsink,aux_real, &
                                        global_auxvars(ghosted_id), &
                                        rich_auxvars(ghosted_id),option)
      call cur_srcsink%Evaluate(res,Jac,PETSC_TRUE, &
                                material_auxvars(ghosted_id), &
                                aux_real,option)
      call MatSetValuesBlockedLocal(Jacobian,1,ghosted_id-1,1, &
                                    ghosted_id-1,Jac,ADD_VALUES, &
                                    ierr);CHKERRQ(ierr)
    else
      iend = local_id*option%nflowdof
      istart = iend - option%nflowdof + 1
      r_p(istart:iend) = r_p(istart:iend) - res
    endif
    cur_srcsink => cur_srcsink%next
  enddo
  
  if (.not.compute_derivative) then
    call VecRestoreArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine RichardsSSSandbox

! ************************************************************************** !

subroutine RichardsSSSandboxLoadAuxReal(srcsink,aux_real,global_auxvar, &
                                        rich_auxvars,option)
! Modified by: Ayman Alzraiee on 04/05/2016 
  use Option_module
  use SrcSink_Sandbox_Base_class
  use SrcSink_Sandbox_Downreg_class

  implicit none

  class(srcsink_sandbox_base_type) :: srcsink
  PetscReal :: aux_real(:)
  type(global_auxvar_type) :: global_auxvar
  type(richards_auxvar_type) :: rich_auxvars
  type(option_type) :: option
  
  aux_real = 0.d0

  !select type(srcsink)
  !  class is(srcsink_sandbox_downreg_type)
      aux_real(1) = rich_auxvars%kvr ! fluid mobility
      aux_real(3) = global_auxvar%pres(1)
      aux_real(9) = global_auxvar%den(1)
  !end select
  
end subroutine RichardsSSSandboxLoadAuxReal

! ************************************************************************** !

subroutine RichardsComputeLateralMassFlux(realization)
  !
  ! Computes lateral flux source/sink term when QUASI_3D is true
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/09/2016
  !

  use Connection_module
  use Realization_Subsurface_class
  use Field_module

  implicit none

  type(realization_subsurface_type) :: realization
  type(field_type), pointer :: field
  PetscErrorCode :: ierr

  field => realization%field

  call RichardsUpdateLocalVecs(field%flow_xx, realization, ierr)

  call RichardsUpdateAuxVarsPatch(realization)

  call VecZeroEntries(field%flow_mass_transfer, ierr); CHKERRQ(ierr)

  call RichardsResidualInternalConn(field%flow_mass_transfer, &
                                    realization, VERT_CONN, ierr)

  call VecScale(field%flow_mass_transfer, -1.d0, ierr); CHKERRQ(ierr)

end subroutine RichardsComputeLateralMassFlux

! ************************************************************************** !

function skip_conn(dist,skip_conn_type)
  !
  ! Returns if a connection should be skipped depending on the skip_conn_type
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/10/2016
  !

  implicit none

  PetscReal :: dist(1:3)
  PetscInt  :: skip_conn_type

  PetscBool :: skip_conn
  PetscBool :: is_conn_vertical

  skip_conn = PETSC_FALSE

  is_conn_vertical = (abs(dot_product(dist(1:3),unit_z)) < 0.99d0)

  select case(skip_conn_type)
    case (HORZ_CONN)
      if (is_conn_vertical) skip_conn = PETSC_TRUE
    case (VERT_CONN)
      if (.not.is_conn_vertical) skip_conn = PETSC_TRUE
  end select

end function skip_conn

! ************************************************************************** !

subroutine RichardsDestroy(realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Realization_Subsurface_class
  
  implicit none

  type(realization_subsurface_type) :: realization
  
  call RichardsDestroyPatch(realization)

end subroutine RichardsDestroy

! ************************************************************************** !

subroutine RichardsDestroyPatch(realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/09
  ! 

  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization
  
  ! taken care of in auxiliary.F90

end subroutine RichardsDestroyPatch

end module Richards_module
