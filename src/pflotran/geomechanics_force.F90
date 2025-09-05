module Geomechanics_Force_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Geomechanics_Global_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public :: GeomechForceSetup, &
            GeomechForceUpdateAuxVars, &
            GeomechanicsForceInitialGuess, &
            GeomechUpdateFromSubsurf, &
            GeomechUpdateSubsurfFromGeomech, &
            GeomechCreateGeomechSubsurfVec, &
            GeomechCreateSubsurfStressStrainVec, &
            GeomechUpdateSolution, &
            GeomechStoreInitialPressTemp, &
            GeomechStoreInitialDisp, &
            GeomechStoreInitialPorosity, &
            GeomechForceSetupLinearSystem, &
            GeomechForceAssembleCoeffMatrix

contains

! ************************************************************************** !

subroutine GeomechForceSetup(geomech_realization)
  !
  ! Sets up the geomechanics calculations
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  !

  use Geomechanics_Realization_class
  use Output_Aux_module

  class(realization_geomech_type) :: geomech_realization
  type(output_variable_list_type), pointer :: list

  call GeomechForceSetupPatch(geomech_realization)

  list => geomech_realization%output_option%output_snap_variable_list
  call GeomechForceSetPlotVariables(list)
  list => geomech_realization%output_option%output_obs_variable_list
  call GeomechForceSetPlotVariables(list)

end subroutine GeomechForceSetup

! ************************************************************************** !

subroutine GeomechForceSetupPatch(geomech_realization)
  !
  ! Sets up the arrays for geomech parameters
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/11/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Option_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch

  PetscInt :: i

  option => geomech_realization%option
  patch => geomech_realization%geomech_patch

  patch%geomech_aux%GeomechParam%youngs_modulus_spatially_varying = .false.
  patch%geomech_aux%GeomechParam%poissons_ratio_spatially_varying = .false.
  patch%geomech_aux%GeomechParam%density_spatially_varying = .false.
  patch%geomech_aux%GeomechParam%biot_coeff_spatially_varying = .false.
  patch%geomech_aux%GeomechParam%thermal_exp_coeff_spatially_varying = .false.
  do i = 1, size(geomech_realization%geomech_material_property_array)
    if (associated(geomech_realization%geomech_material_property_array(i)%ptr%youngs_modulus_dataset)) then
      patch%geomech_aux%GeomechParam%youngs_modulus_spatially_varying = .true.
    endif
    if (associated(geomech_realization%geomech_material_property_array(i)%ptr%poissons_ratio_dataset)) then
      patch%geomech_aux%GeomechParam%poissons_ratio_spatially_varying = .true.
    endif
    if (associated(geomech_realization%geomech_material_property_array(i)%ptr%density_dataset)) then
      patch%geomech_aux%GeomechParam%density_spatially_varying = .true.
    endif
    if (associated(geomech_realization%geomech_material_property_array(i)%ptr%biot_coeff_dataset)) then
      patch%geomech_aux%GeomechParam%biot_coeff_spatially_varying = .true.
    endif
    if (associated(geomech_realization%geomech_material_property_array(i)%ptr%thermal_exp_coeff_dataset)) then
      patch%geomech_aux%GeomechParam%thermal_exp_coeff_spatially_varying = .true.
    endif
  enddo

  ! Young's modulus
  if (patch%geomech_aux%GeomechParam%youngs_modulus_spatially_varying) then
    patch%geomech_aux%GeomechParam%youngs_modulus_vec = patch%geomech_field%youngs_modulus
  else
    allocate(patch%geomech_aux%GeomechParam%youngs_modulus &
      (size(geomech_realization%geomech_material_property_array)))
    do i = 1, size(geomech_realization%geomech_material_property_array)
      patch%geomech_aux%GeomechParam%youngs_modulus(geomech_realization% &
        geomech_material_property_array(i)%ptr%id) = geomech_realization% &
        geomech_material_property_array(i)%ptr%youngs_modulus
    enddo
  endif
  ! Poisson's ratio
  if (patch%geomech_aux%GeomechParam%poissons_ratio_spatially_varying) then
    patch%geomech_aux%GeomechParam%poissons_ratio_vec = patch%geomech_field%poissons_ratio
  else
    allocate(patch%geomech_aux%GeomechParam%poissons_ratio &
      (size(geomech_realization%geomech_material_property_array)))
    do i = 1, size(geomech_realization%geomech_material_property_array)
      patch%geomech_aux%GeomechParam%poissons_ratio(geomech_realization% &
        geomech_material_property_array(i)%ptr%id) = geomech_realization% &
        geomech_material_property_array(i)%ptr%poissons_ratio
    enddo
  endif
  ! Density
  if (patch%geomech_aux%GeomechParam%density_spatially_varying) then
    patch%geomech_aux%GeomechParam%density_vec = patch%geomech_field%density
  else
    allocate(patch%geomech_aux%GeomechParam%density &
      (size(geomech_realization%geomech_material_property_array)))
    do i = 1, size(geomech_realization%geomech_material_property_array)
      patch%geomech_aux%GeomechParam%density(geomech_realization% &
        geomech_material_property_array(i)%ptr%id) = geomech_realization% &
        geomech_material_property_array(i)%ptr%density
    enddo
  endif
  ! Biot's coefficient
  if (patch%geomech_aux%GeomechParam%biot_coeff_spatially_varying) then
    patch%geomech_aux%GeomechParam%biot_coeff_vec = patch%geomech_field%biot_coeff
  else
    allocate(patch%geomech_aux%GeomechParam%biot_coeff &
      (size(geomech_realization%geomech_material_property_array)))
    do i = 1, size(geomech_realization%geomech_material_property_array)
      patch%geomech_aux%GeomechParam%biot_coeff(geomech_realization% &
        geomech_material_property_array(i)%ptr%id) = geomech_realization% &
        geomech_material_property_array(i)%ptr%biot_coeff
    enddo
  endif
  ! Thermal expansion coefficient
  if (patch%geomech_aux%GeomechParam%thermal_exp_coeff_spatially_varying) then
    patch%geomech_aux%GeomechParam%thermal_exp_coeff_vec = patch%geomech_field%thermal_exp_coeff
  else
    allocate(patch%geomech_aux%GeomechParam%thermal_exp_coeff &
      (size(geomech_realization%geomech_material_property_array)))
    do i = 1, size(geomech_realization%geomech_material_property_array)
      patch%geomech_aux%GeomechParam%thermal_exp_coeff(geomech_realization% &
        geomech_material_property_array(i)%ptr%id) = geomech_realization% &
        geomech_material_property_array(i)%ptr%thermal_exp_coeff
    enddo
  endif

end subroutine GeomechForceSetupPatch

! ************************************************************************** !

subroutine GeomechForceSetPlotVariables(list)
  !
  ! Set up of geomechanics plot variables
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  !

  use Output_Aux_module
  use Variables_module

  implicit none

  type(output_variable_list_type), pointer :: list
  type(output_variable_type), pointer :: output_variable

  character(len=MAXWORDLENGTH) :: name, units

  if (associated(list%first)) then
    return
  endif

  name = 'displacement_x'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_DISP_X)

  name = 'displacement_y'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_DISP_Y)

  name = 'displacement_z'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_DISP_Z)

  units = ''
  name = 'Material ID'
  output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                          units,GEOMECH_MATERIAL_ID)
  output_variable%iformat = 1 ! integer
  call OutputVariableAddToList(list,output_variable)

  name = 'volumetric strain'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               GEOMECH_VOLUMETRIC_STRAIN)

  name = 'strain_xx'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_XX)

  name = 'strain_yy'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_YY)

  name = 'strain_zz'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_ZZ)

  name = 'strain_xy'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_XY)

  name = 'strain_yz'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_YZ)

  name = 'strain_zx'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_ZX)

  name = 'stress_xx'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_XX)

  name = 'stress_yy'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_YY)

  name = 'stress_zz'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_ZZ)

  name = 'stress_xy'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_XY)

  name = 'stress_yz'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_YZ)

  name = 'stress_zx'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_ZX)

  name = 'stress_total_xx'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_XX)

  name = 'stress_total_yy'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_YY)

  name = 'stress_total_zz'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_ZZ)

  name = 'stress_total_xy'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_XY)

  name = 'stress_total_yz'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_YZ)

  name = 'stress_total_zx'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_ZX)

  name = 'relative_displacement_x'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_REL_DISP_X)

  name = 'relative_displacement_y'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_REL_DISP_Y)

  name = 'relative_displacement_z'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_REL_DISP_Z)


end subroutine GeomechForceSetPlotVariables

! ************************************************************************** !

subroutine GeomechanicsForceInitialGuess(geomech_realization)
  !
  ! Sets up the inital guess for the solution
  ! The boundary conditions are set here
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/19/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Option_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Patch_module
  use Geomechanics_Coupler_module
  use Geomechanics_Region_module

  implicit none

  class(realization_geomech_type) :: geomech_realization

  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: field
  type(geomech_patch_type), pointer :: patch
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_grid_type), pointer :: grid
  type(gm_region_type), pointer :: region

  PetscInt :: ghosted_id,local_id,total_verts,ivertex
  PetscReal, pointer :: xx_p(:)
  PetscErrorCode :: ierr

  option => geomech_realization%option
  field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid

  call VecGetArrayf90(field%disp_xx,xx_p,ierr);CHKERRQ(ierr)

  boundary_condition => patch%geomech_boundary_condition_list%first
  total_verts = 0
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      total_verts = total_verts + 1
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      ! X displacement
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_X_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_X_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

      ! Y displacement
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_Y_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_Y_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

      ! Z displacement
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_Z_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_Z_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayf90(field%disp_xx,xx_p,ierr);CHKERRQ(ierr)

end subroutine GeomechanicsForceInitialGuess

! ************************************************************************** !

subroutine GeomechForceUpdateAuxVars(geomech_realization)
  !
  ! Updates the geomechanics variables
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/18/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Option_module
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Coupler_module
  use Geomechanics_Material_module
  use Geomechanics_Global_Aux_module
  use Geomechanics_Region_module

  implicit none

  class(realization_geomech_type) :: geomech_realization

  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch
  type(geomech_grid_type), pointer :: grid
  type(geomech_field_type), pointer :: geomech_field
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)

  PetscInt :: ghosted_id
  PetscReal, pointer :: xx_loc_p(:), xx_init_loc_p(:)
  PetscErrorCode :: ierr

  option => geomech_realization%option
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  geomech_field => geomech_realization%geomech_field

  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars

  ! Communication -----------------------------------------
  call GeomechDiscretizationGlobalToLocal(geomech_realization% &
                                          geomech_discretization, &
                                          geomech_realization% &
                                          geomech_field%disp_xx, &
                                          geomech_realization% &
                                          geomech_field%disp_xx_loc,NGEODOF)

  call VecGetArrayF90(geomech_field%disp_xx_loc,xx_loc_p,ierr)
  call VecGetArrayF90(geomech_field%disp_xx_init_loc,xx_init_loc_p,ierr)

  ! Internal aux vars
  do ghosted_id = 1, grid%ngmax_node
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_X_DOF) = &
      xx_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_Y_DOF) = &
      xx_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_Z_DOF) = &
      xx_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER)

    geomech_global_aux_vars(ghosted_id)%rel_disp_vector(GEOMECH_DISP_X_DOF) = &
      xx_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER) - &
      xx_init_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%rel_disp_vector(GEOMECH_DISP_Y_DOF) = &
      xx_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER) - &
      xx_init_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%rel_disp_vector(GEOMECH_DISP_Z_DOF) = &
      xx_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER) - &
      xx_init_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER)
 enddo

  call VecRestoreArrayF90(geomech_field%disp_xx_loc,xx_loc_p,ierr)
  call VecRestoreArrayF90(geomech_field%disp_xx_init_loc,xx_init_loc_p,ierr)


end subroutine GeomechForceUpdateAuxVars

! ************************************************************************** !

subroutine GeomechForceResidual(snes,xx,r,geomech_realization,ierr)
  !
  ! Computes the residual equation
  !
  ! Author: Satish Karra
  ! Date: 06/21/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  class(realization_geomech_type) :: geomech_realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_field_type), pointer :: field
  type(option_type), pointer :: option

  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  option => geomech_realization%option

  ! Communication -----------------------------------------
  call GeomechDiscretizationGlobalToLocal(geomech_discretization,xx, &
                                          field%disp_xx_loc,NGEODOF)

  call GeomechForceResidualPatch(snes,xx,r,geomech_realization,ierr)

  if (geomech_realization%geomech_debug%vecview_residual) then
    call PetscViewerASCIIOpen(geomech_realization%option%mycomm, &
                              'Geomech_residual.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  endif

  if (geomech_realization%geomech_debug%vecview_solution) then
    call PetscViewerASCIIOpen(geomech_realization%option%mycomm, &
                              'Geomech_xx.out',viewer,ierr);CHKERRQ(ierr)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

end subroutine GeomechForceResidual

! ************************************************************************** !

subroutine GeomechForceResidualPatch(snes,xx,r,geomech_realization,ierr)
  !
  ! Computes the residual equation on a patch
  !
  ! Author: Satish Karra
  ! Date: 06/24/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Grid_Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Coupler_module
  use Option_module
  use Geomechanics_Auxiliary_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  class(realization_geomech_type) :: geomech_realization
  PetscErrorCode :: ierr

  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(gm_region_type), pointer :: region
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_parameter_type), pointer :: GeomechParam

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscReal, allocatable :: local_press(:), local_temp(:)
  PetscInt, allocatable :: petsc_ids(:)
  PetscInt, allocatable :: ids(:)
  PetscReal, allocatable :: res_vec(:)
  PetscReal, pointer :: press(:), temp(:)
  PetscReal, pointer :: fluid_density(:), porosity(:)
  PetscReal, pointer :: press_init(:), temp_init(:)
  PetscReal, pointer :: fluid_density_init(:)
  PetscReal, allocatable :: beta_vec(:), alpha_vec(:)
  PetscReal, allocatable :: density_rock_vec(:), density_fluid_vec(:)
  PetscReal, allocatable :: density_bulk_vec(:)
  PetscReal, allocatable :: youngs_vec(:), poissons_vec(:)
  PetscReal, allocatable :: porosity_vec(:)
  PetscInt :: ielem, ivertex
  PetscInt :: ghosted_id
  PetscInt :: eletype, idof
  PetscInt :: petsc_id, local_id
  PetscReal, pointer :: imech_loc_p(:)
  PetscInt :: size_elenodes
  PetscInt :: facetype, nfaces

  PetscInt :: iface, num_vertices
  PetscReal :: stress_bc(SIX_INTEGER)
  PetscInt, allocatable :: face_vertices(:)

  PetscReal, pointer :: temp_youngs_modulus_p(:)
  PetscReal, pointer :: temp_poissons_ratio_p(:)
  PetscReal, pointer :: temp_density_p(:)
  PetscReal, pointer :: temp_biot_coeff_p(:)
  PetscReal, pointer :: temp_thermal_exp_coeff_p(:)

  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars
  GeomechParam => patch%geomech_aux%GeomechParam

  call GeomechForceUpdateAuxVars(geomech_realization)
  ! Add flag for the update

  call VecSet(r,0.d0,ierr);CHKERRQ(ierr)

#if 0
  error_H1_global = 0.d0
  error_L2_global = 0.d0
#endif

  ! Get pressure and temperature from subsurface
  call VecGetArrayF90(field%press_loc,press,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%temp_loc,temp,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%porosity_loc,porosity,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%fluid_density_loc,fluid_density,ierr);CHKERRQ(ierr)

  ! Get initial pressure and temperature
  call VecGetArrayF90(field%press_init_loc,press_init,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%temp_init_loc,temp_init,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%fluid_density_init_loc,fluid_density_init, &
                      ierr);CHKERRQ(ierr)

  if (GeomechParam%youngs_modulus_spatially_varying) then
    call VecGetArrayF90(field%youngs_modulus,temp_youngs_modulus_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%poissons_ratio_spatially_varying) then
    call VecGetArrayF90(field%poissons_ratio,temp_poissons_ratio_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%density_spatially_varying) then
    call VecGetArrayF90(field%density,temp_density_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%biot_coeff_spatially_varying) then
    call VecGetArrayF90(field%biot_coeff,temp_biot_coeff_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%thermal_exp_coeff_spatially_varying) then
    call VecGetArrayF90(field%thermal_exp_coeff,temp_thermal_exp_coeff_p,ierr);CHKERRQ(ierr)
  endif

  ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    allocate(elenodes(grid%elem_nodes(0,ielem)))
    allocate(local_coordinates(size(elenodes),THREE_INTEGER))
    allocate(local_disp(size(elenodes),option%ngeomechdof))
    allocate(local_press(size(elenodes)))
    allocate(local_temp(size(elenodes)))
    allocate(petsc_ids(size(elenodes)))
    allocate(ids(size(elenodes)*option%ngeomechdof))
    allocate(res_vec(size(elenodes)*option%ngeomechdof))
    allocate(beta_vec(size(elenodes)))
    allocate(alpha_vec(size(elenodes)))
    allocate(density_rock_vec(size(elenodes)))
    allocate(density_fluid_vec(size(elenodes)))
    allocate(youngs_vec(size(elenodes)))
    allocate(poissons_vec(size(elenodes)))
    allocate(porosity_vec(size(elenodes)))
    allocate(density_bulk_vec(size(elenodes)))
    elenodes = grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    eletype = grid%gauss_node(ielem)%entity_type
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
    enddo
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        local_disp(ivertex,idof) = &
          geomech_global_aux_vars(ghosted_id)%disp_vector(idof)
        ids(idof + (ivertex-1)*option%ngeomechdof) = &
          (petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
      enddo
      local_press(ivertex) = press(ghosted_id) - press_init(ghosted_id)  ! p - p_0
      local_temp(ivertex) = temp(ghosted_id) - temp_init(ghosted_id)     ! T - T_0
      if (GeomechParam%thermal_exp_coeff_spatially_varying) then
        alpha_vec(ivertex) = temp_thermal_exp_coeff_p(grid%nG2L(ghosted_id))
      else
        alpha_vec(ivertex) = &
          GeomechParam%thermal_exp_coeff(nint(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%biot_coeff_spatially_varying) then
        beta_vec(ivertex) = temp_biot_coeff_p(grid%nG2L(ghosted_id))
      else
        beta_vec(ivertex) = &
          GeomechParam%biot_coeff(nint(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%density_spatially_varying) then
        density_rock_vec(ivertex) = temp_density_p(grid%nG2L(ghosted_id))
      else
        density_rock_vec(ivertex) = &
          GeomechParam%density(nint(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%youngs_modulus_spatially_varying) then
        youngs_vec(ivertex) = temp_youngs_modulus_p(grid%nG2L(ghosted_id))
      else
        youngs_vec(ivertex) = &
          GeomechParam%youngs_modulus(nint(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%poissons_ratio_spatially_varying) then
        poissons_vec(ivertex) = temp_poissons_ratio_p(grid%nG2L(ghosted_id))
      else
        poissons_vec(ivertex) = &
          GeomechParam%poissons_ratio(nint(imech_loc_p(ghosted_id)))
      endif
      density_fluid_vec(ivertex) = fluid_density(ghosted_id)
      porosity_vec(ivertex) = porosity(ghosted_id)
      density_bulk_vec(ivertex) = (porosity_vec(ivertex) * &
                                   density_fluid_vec(ivertex)) + &
                                  ((1.d0 - porosity_vec(ivertex)) * &
                                   density_rock_vec(ivertex))
    enddo
    size_elenodes = size(elenodes)
    call GeomechForceLocalElemResidual(size_elenodes,local_coordinates, &
       local_disp,local_press,local_temp,youngs_vec,poissons_vec, &
       density_bulk_vec,beta_vec,alpha_vec,eletype, &
       grid%gauss_node(ielem)%dim,grid%gauss_node(ielem)%r, &
       grid%gauss_node(ielem)%w,res_vec,option)
    call VecSetValues(r,size(ids),ids,res_vec,ADD_VALUES,ierr);CHKERRQ(ierr)
#if 0
    call GeomechForceLocalElemError(size_elenodes,local_coordinates, &
                                    local_disp, &
                                    eletype,grid%gauss_node(ielem)%dim, &
                                    grid%gauss_node(ielem)%r, &
                                    grid%gauss_node(ielem)%w,error_L2, &
                                    error_H1,option)
    error_H1_global = error_H1_global + error_H1
    error_L2_global = error_L2_global + error_L2
#endif
    deallocate(elenodes)
    deallocate(local_coordinates)
    deallocate(local_disp)
    deallocate(petsc_ids)
    deallocate(ids)
    deallocate(res_vec)
    deallocate(local_press)
    deallocate(local_temp)
    deallocate(beta_vec)
    deallocate(alpha_vec)
    deallocate(density_rock_vec)
    deallocate(density_fluid_vec)
    deallocate(youngs_vec)
    deallocate(poissons_vec)
    deallocate(porosity_vec)
    deallocate(density_bulk_vec)
  enddo

  call VecRestoreArrayF90(field%press_loc,press,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%temp_loc,temp,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%fluid_density_loc,fluid_density, &
                          ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%press_init_loc,press_init,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%temp_init_loc,temp_init,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%fluid_density_init_loc,fluid_density_init, &
                          ierr);CHKERRQ(ierr)

  if (GeomechParam%youngs_modulus_spatially_varying) then
    call VecRestoreArrayF90(field%youngs_modulus,temp_youngs_modulus_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%poissons_ratio_spatially_varying) then
    call VecRestoreArrayF90(field%poissons_ratio,temp_poissons_ratio_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%density_spatially_varying) then
    call VecRestoreArrayF90(field%density,temp_density_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%biot_coeff_spatially_varying) then
    call VecRestoreArrayF90(field%biot_coeff,temp_biot_coeff_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%thermal_exp_coeff_spatially_varying) then
    call VecRestoreArrayF90(field%thermal_exp_coeff,temp_thermal_exp_coeff_p,ierr);CHKERRQ(ierr)
  endif

#if 0
  call MPI_Allreduce(error_H1_global,error_H1_global,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                     ierr);CHKERRQ(ierr)
  call MPI_Allreduce(error_L2_global,error_L2_global,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                     ierr);CHKERRQ(ierr)

  if (OptionIsIORank(option)) then
    print *, 'L2 error:', sqrt(error_L2_global)
    print *, 'H1 error:', sqrt(error_H1_global)
  endif
#endif


  call VecAssemblyBegin(r,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(r,ierr);CHKERRQ(ierr)

  ! Find the boundary nodes with dirichlet and set the residual at those nodes
  ! to zero, later set the Jacobian to 1

  ! displacement boundary conditions
  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      ! X displacement
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_X_DOF-1, &
                             0.d0,INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for displacement not available.'
            call PrintErrMsg(option)
        end select
      endif

      ! Y displacement
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_Y_DOF-1, &
                             0.d0,INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for displacement not available.'
            call PrintErrMsg(option)
        end select
      endif

      ! Z displacement
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_Z_DOF-1, &
                             0.d0,INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for displacement not available.'
            call PrintErrMsg(option)
        end select
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Need to assemby here since one cannot mix INSERT_VALUES
  ! and ADD_VALUES
  call VecAssemblyBegin(r,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(r,ierr);CHKERRQ(ierr)

  ! Force boundary conditions
  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      ! X force
      if (associated(boundary_condition%geomech_condition%force_x)) then
        select case(boundary_condition%geomech_condition%force_x%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_X_DOF-1, &
                             -boundary_condition% &
                               geomech_aux_real_var(GEOMECH_DISP_X_DOF,ivertex), &
                             ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for force not available.'
            call PrintErrMsg(option)
        end select
      endif

       ! Y force
      if (associated(boundary_condition%geomech_condition%force_y)) then
        select case(boundary_condition%geomech_condition%force_y%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_Y_DOF-1, &
                             -boundary_condition% &
                               geomech_aux_real_var(GEOMECH_DISP_Y_DOF,ivertex), &
                             ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for force not available.'
            call PrintErrMsg(option)

        end select
      endif

       ! Z force
      if (associated(boundary_condition%geomech_condition%force_z)) then
        select case(boundary_condition%geomech_condition%force_z%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_Z_DOF-1, &
                             -boundary_condition% &
                               geomech_aux_real_var(GEOMECH_DISP_Z_DOF,ivertex), &
                             ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for force not available.'
            call PrintErrMsg(option)
        end select
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecAssemblyBegin(r,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(r,ierr);CHKERRQ(ierr)

  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    ! Traction
    if (associated(boundary_condition%geomech_condition%traction)) then
      select case(boundary_condition%geomech_condition%traction%itype)
        case(DIRICHLET_BC)
          option%io_buffer = 'Dirichlet BC for traction not available.'
          call PrintErrMsg(option)
        case(ZERO_GRADIENT_BC)
         ! do nothing
        case(NEUMANN_BC)
          stress_bc = boundary_condition%geomech_condition%traction% &
                        dataset%rarray
          nfaces = boundary_condition%region%sideset%nfaces
          do iface = 1, nfaces
            num_vertices = size(boundary_condition%region%sideset% &
                             face_vertices(:,iface))
            if (num_vertices == THREE_INTEGER) then
              facetype = TRI_FACE_TYPE
            else
              facetype = QUAD_FACE_TYPE
            endif
            allocate(face_vertices(num_vertices))
            face_vertices = boundary_condition%region%sideset% &
                              face_vertices(1:num_vertices,iface)
            allocate(local_coordinates(num_vertices,THREE_INTEGER))
            allocate(petsc_ids(num_vertices))
            allocate(ids(num_vertices*option%ngeomechdof))
            allocate(res_vec(num_vertices*option%ngeomechdof))
            do ivertex = 1, num_vertices
              ghosted_id = face_vertices(ivertex)
              local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = &
                grid%nodes(ghosted_id)%x
              local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = &
                grid%nodes(ghosted_id)%y
              local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = &
                grid%nodes(ghosted_id)%z
              petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
              do idof = 1, option%ngeomechdof
                ids(idof + (ivertex-1)*option%ngeomechdof) = &
                  (petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
              enddo
            enddo
            call GeomechForceApplyTractionBCtoResidual( &
                   local_coordinates, &
                   facetype, &
                   stress_bc, &
                   grid%gauss_surf_node(facetype)%r, &
                   grid%gauss_surf_node(facetype)%w, &
                   res_vec,option)
            call VecSetValues(r,size(ids),ids,res_vec,ADD_VALUES, &
                   ierr);CHKERRQ(ierr)
            deallocate(local_coordinates)
            deallocate(petsc_ids)
            deallocate(ids)
            deallocate(res_vec)
          enddo
      end select
    endif
    boundary_condition => boundary_condition%next
  enddo
  call VecAssemblyBegin(r,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(r,ierr);CHKERRQ(ierr)

end subroutine GeomechForceResidualPatch

! ************************************************************************** !

subroutine ComputeTetVolAtVertex(vert_0, &
                                 vert_1, &
                                 vert_2, &
                                 vert_3, &
                                 volume)
  !
  ! Returns the volume of the specified tetrahedron after
  ! being clipped on the positive side of the planes computed
  ! in this function.
  ! This works by calling ClippedVolume() recursively on one plane
  ! at a time.
  !
  ! Author: Joe Eyles, WSP
  ! Date: 29/04/2025
  !


  PetscReal :: vert_0(3)
  PetscReal :: vert_1(3)
  PetscReal :: vert_2(3)
  PetscReal :: vert_3(3)

  PetscReal :: midPoints(3, 3), normals(3, 3)
  PetscReal :: volume

  midPoints(1, :) = (vert_0(:) + vert_1(:))*0.5
  midPoints(2, :) = (vert_0(:) + vert_2(:))*0.5
  midPoints(3, :) = (vert_0(:) + vert_3(:))*0.5

  normals(1, :) = vert_0(:) - vert_1(:)
  normals(2, :) = vert_0(:) - vert_2(:)
  normals(3, :) = vert_0(:) - vert_3(:)

  normals(1, :) = normals(1, :) / sqrt(dot_product(normals(1, :), normals(1, :)))
  normals(2, :) = normals(2, :) / sqrt(dot_product(normals(2, :), normals(2, :)))
  normals(3, :) = normals(3, :) / sqrt(dot_product(normals(3, :), normals(3, :)))

  volume = ClippedVolume(vert_0, vert_1, vert_2, vert_3, midPoints, normals, 3)
end subroutine ComputeTetVolAtVertex

! ************************************************************************** !

recursive function ClippedVolume(vert_0, vert_1, vert_2, vert_3, midPoints, normals, nSize_in) result(nRet)
  !
  ! Clips on the given plane, splits the resulting
  ! polyhedron into tetrahedra, and then calls itself recursively,
  ! ultimately returning the volume.
  !
  ! Author: Joe Eyles, WSP
  ! Date: 29/04/2025
  !

  PetscReal :: vert_0(3), vert_1(3), vert_2(3), vert_3(3)
  PetscReal :: midPoints(3, 3), normals(3, 3)
  PetscInt :: nSize_in
  PetscInt :: nSize
  PetscReal :: nRet
  PetscReal :: vP(3), vN(3), vvert_0(3), vvert_1(3), vvert_2(3), vvert_3(3), vvert_4(3)
  PetscReal :: V(4, 3)
  PetscReal :: d(4)

  nRet = 0
  nSize = nSize_in

  if (nSize == 0) then
    nRet = abs(dot_product(cross_product(vert_1 - vert_0, vert_2 - vert_0), vert_3 - vert_0) / 6.0)
  else
    nSize = nSize - 1
    vP = midPoints(nSize + 1, :)
    vN = normals(nSize + 1, :)
    V(1, :) = vert_0
    d(1) = dot_product(vert_0 - vP, vN)
    V(2, :) = vert_1
    d(2) = dot_product(vert_1 - vP, vN)
    V(3, :) = vert_2
    d(3) = dot_product(vert_2 - vP, vN)
    V(4, :) = vert_3
    d(4) = dot_product(vert_3 - vP, vN)

    call sort(V, d)

    if (d(1) <= 0.0) then
      nRet = 0.0
    else if (d(4) >= 0.0) then
      nRet = nRet + ClippedVolume(vert_0, vert_1, vert_2, vert_3, midPoints, normals, nSize)
    else if (d(3) > 0.0) then
      vvert_0 = Intersect(V(1, :), V(4, :), d(1), d(4))
      vvert_1 = Intersect(V(2, :), V(4, :), d(2), d(4))
      vvert_2 = Intersect(V(3, :), V(4, :), d(3), d(4))
      nRet = nRet + ClippedVolume(V(1, :), V(2, :), V(3, :), vvert_2, midPoints, normals, nSize) + &
             ClippedVolume(V(1, :), vvert_0, V(2, :), vvert_2, midPoints, normals, nSize) + &
             ClippedVolume(V(2, :), vvert_0, vvert_1, vvert_2, midPoints, normals, nSize)
    else if (d(3) == 0.0) then
      vvert_0 = Intersect(V(1, :), V(4, :), d(1), d(4))
      vvert_1 = Intersect(V(2, :), V(4, :), d(2), d(4))
      nRet = nRet + ClippedVolume(V(1, :), vvert_0, V(2, :), V(3, :), midPoints, normals, nSize) + &
             ClippedVolume(V(2, :), vvert_0, vvert_1, V(3, :), midPoints, normals, nSize)
    else if (d(2) > 0.0) then
      vvert_1 = Intersect(V(1, :), V(3, :), d(1), d(3))
      vvert_2 = Intersect(V(1, :), V(4, :), d(1), d(4))
      vvert_3 = Intersect(V(2, :), V(3, :), d(2), d(3))
      vvert_4 = Intersect(V(2, :), V(4, :), d(2), d(4))
      nRet = nRet + ClippedVolume(V(1, :), vvert_1, vvert_2, vvert_3, midPoints, normals, nSize) + &
             ClippedVolume(V(1, :), vvert_2, vvert_3, vvert_4, midPoints, normals, nSize) + &
             ClippedVolume(V(1, :), V(2, :), vvert_3, vvert_4, midPoints, normals, nSize)
    else
      vvert_1 = Intersect(V(1, :), V(2, :), d(1), d(2))
      vvert_2 = Intersect(V(1, :), V(3, :), d(1), d(3))
      vvert_3 = Intersect(V(1, :), V(4, :), d(1), d(4))
      nRet = nRet + ClippedVolume(V(1, :), vvert_1, vvert_2, vvert_3, midPoints, normals, nSize)
    end if
  end if
end function ClippedVolume

! ************************************************************************** !

function Intersect(v1, v2, d1, d2) result(v)
  !
  ! Computes the intesection
  !
  ! Author: Joe Eyles, WSP
  ! Date: 29/04/2025
  !

  PetscReal :: v1(3), v2(3)
  PetscReal :: d1, d2
  PetscReal :: v(3)

  v = v1 * (-d2 / (d1 - d2)) + v2 * (d1 / (d1 - d2))
end function Intersect

! ************************************************************************** !

function cross_product(v1, v2) result(v)
  !
  ! Computes the cross product between v1 and v2
  !
  ! Author: Joe Eyles, WSP
  ! Date: 29/04/2025
  !

  PetscReal :: v1(3), v2(3)
  PetscReal :: v(3)

  v(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v(3) = v1(1) * v2(2) - v1(2) * v2(1)
end function cross_product

! ************************************************************************** !

subroutine sort(V, d)
  !
  ! Sorts vectors V and d based on d
  !
  ! Author: Joe Eyles, WSP
  ! Date: 29/04/2025
  !

  PetscReal :: V(4, 3)
  PetscReal :: d(4)
  integer :: i, j
  PetscReal :: temp_d
  PetscReal :: temp_V(3)

  do i = 1, 3
    do j = i + 1, 4
      if (d(i) < d(j)) then
        temp_d = d(i)
        d(i) = d(j)
        d(j) = temp_d
        temp_V = V(i, :)
        V(i, :) = V(j, :)
        V(j, :) = temp_V
      end if
    end do
  end do
end subroutine sort

! ************************************************************************** !

function face_unitnormal(v1, v2, v3) result(n)
  !
  ! calculates the outward pointing normal vector
  ! for tri face, pass coordinates in the same order
  ! for quad face, pass (v1, v2, v4)
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 09/04/2025
  !

  PetscReal :: v1(THREE_INTEGER), v2(THREE_INTEGER), v3(THREE_INTEGER)
  PetscReal :: n(THREE_INTEGER)

  PetscReal :: e1(THREE_INTEGER), e2(THREE_INTEGER)

  e1 = v2 - v1
  e2 = v3 - v1
  n = cross_product(e1,e2)
  n = n/sqrt(dot_product(n,n))

end function face_unitnormal

! ************************************************************************** !

subroutine GeomechForceLocalElemResidual(size_elenodes,local_coordinates, &
                                         local_disp, &
                                         local_press,local_temp, &
                                         local_youngs,local_poissons, &
                                         local_density,local_beta, &
                                         local_alpha, &
                                         eletype,dim,r,w,res_vec,option)
  !
  ! Computes the residual for a local element
  !
  ! Author: Satish Karra
  ! Date: 06/24/13
  !

  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module

  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: B(:,:), Kmat(:,:)
  PetscReal, allocatable :: res_vec(:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscReal, allocatable :: local_press(:)
  PetscReal, allocatable :: local_temp(:)
  PetscReal, allocatable :: local_youngs(:)
  PetscReal, allocatable :: local_poissons(:)
  PetscReal, allocatable :: local_density(:)
  PetscReal, allocatable :: local_beta(:)
  PetscReal, allocatable :: local_alpha(:)

  PetscReal, pointer :: r(:,:), w(:)
  PetscInt :: igpt
  PetscInt :: len_w
  PetscInt :: eletype
  PetscReal :: x(THREE_INTEGER), J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: inv_J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: detJ_map
  PetscInt :: i,j
  PetscInt :: dim
  PetscReal :: lambda, mu, beta, alpha
  PetscReal :: density, youngs_mod, poissons_ratio
  PetscInt :: load_type
  PetscReal :: bf(THREE_INTEGER)
  PetscReal :: identity(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: gauss_tot_weight
  PetscReal, allocatable :: N(:,:)
  PetscReal, allocatable :: gauss_tet_vol_weight(:)
  PetscReal, allocatable :: vecB_transpose(:,:)
  PetscReal, allocatable :: kron_B_eye(:,:)
  PetscReal, allocatable :: kron_B_transpose_eye(:,:)
  PetscReal, allocatable :: Trans(:,:)
  PetscReal, allocatable :: kron_eye_B_transpose(:,:)
  PetscReal, allocatable :: kron_N_eye(:,:)
  PetscReal, allocatable :: vec_local_disp(:,:)
  PetscReal, allocatable :: force(:), res_vec_mat(:,:)
  PetscInt :: size_elenodes

  allocate(B(size_elenodes,dim))
  allocate(Kmat(size_elenodes*option%ngeomechdof, &
                size_elenodes*option%ngeomechdof))
  allocate(force(size_elenodes*option%ngeomechdof))
  allocate(res_vec_mat(size_elenodes*option%ngeomechdof,1))

  res_vec = 0.d0
  res_vec_mat = 0.d0
  Kmat = 0.d0
  force = 0.d0
  len_w = size(w)

  identity = 0.d0
  do i = 1, THREE_INTEGER
    do j = 1, THREE_INTEGER
      if (i == j) identity(i,j) = 1.d0
    enddo
  enddo

  if(option%geomechanics%improve_tet_weighting .and. &
    eletype == TET_TYPE .and. len_w == 4) then
    allocate(gauss_tet_vol_weight(len_w))
    call ComputeTetVolAtVertex(local_coordinates(1, :), &
                               local_coordinates(2, :), &
                               local_coordinates(3, :), &
                               local_coordinates(4, :), &
                               gauss_tet_vol_weight(1))

    call ComputeTetVolAtVertex(local_coordinates(2, :), &
                               local_coordinates(3, :), &
                               local_coordinates(4, :), &
                               local_coordinates(1, :), &
                               gauss_tet_vol_weight(2))

    call ComputeTetVolAtVertex(local_coordinates(3, :), &
                               local_coordinates(4, :), &
                               local_coordinates(1, :), &
                               local_coordinates(2, :), &
                               gauss_tet_vol_weight(3))

    call ComputeTetVolAtVertex(local_coordinates(4, :), &
                               local_coordinates(1, :), &
                               local_coordinates(2, :), &
                               local_coordinates(3, :), &
                               gauss_tet_vol_weight(4))

    gauss_tot_weight = gauss_tet_vol_weight(1) + &
                       gauss_tet_vol_weight(2) + &
                       gauss_tet_vol_weight(3) + &
                       gauss_tet_vol_weight(4)
    do i = 1, 4
      gauss_tet_vol_weight(i) = gauss_tet_vol_weight(i) / gauss_tot_weight
    enddo
  endif

  call Transposer(option%ngeomechdof,size_elenodes,Trans)

  do igpt = 1, len_w
    shapefunction%element_type = eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction)
    x = matmul(transpose(local_coordinates),shapefunction%N)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)
    allocate(N(size(shapefunction%N),ONE_INTEGER))
    call Determinant(J_map,detJ_map)
    if (detJ_map <= 0.d0) then
      option%io_buffer = 'GEOMECHANICS: Determinant of J_map has' // &
                         ' to be positive!'
      call PrintErrMsg(option)
    endif
    ! Find the inverse of J_map
    call MatInv3(J_map,inv_J_map)
    B = matmul(shapefunction%DN,inv_J_map)
    youngs_mod = dot_product(shapefunction%N,local_youngs)
    poissons_ratio = dot_product(shapefunction%N,local_poissons)
    alpha = dot_product(shapefunction%N,local_alpha)
    beta = dot_product(shapefunction%N,local_beta)
    density = dot_product(shapefunction%N,local_density)
    call GeomechGetLambdaMu(lambda,mu,youngs_mod,poissons_ratio)
    call GeomechGetBodyForce(load_type,lambda,mu,x,bf,option)
    call ConvertMatrixToVector(transpose(B),vecB_transpose)
    Kmat = Kmat + w(igpt)*lambda* &
      matmul(vecB_transpose,transpose(vecB_transpose))*detJ_map
    call Kron(B,identity,kron_B_eye)
    call Kron(transpose(B),identity,kron_B_transpose_eye)
    call Kron(identity,transpose(B),kron_eye_B_transpose)
    N(:,1)= shapefunction%N
    call Kron(N,identity,kron_N_eye)
    Kmat = Kmat + w(igpt)*mu*matmul(kron_B_eye,kron_B_transpose_eye)*detJ_map
    Kmat = Kmat + w(igpt)*mu* &
      matmul(matmul(kron_B_eye,kron_eye_B_transpose),Trans)*detJ_map
    if(option%geomechanics%improve_tet_weighting .and. &
       eletype == TET_TYPE .and. len_w == 4) then
      ! w(igpt) = 1/4 * 1/6 = 1/n_gausspoints * reference_tet_volume
      ! and detJ_map * reference_tet_volume = current_tet_volume
      ! so w(igpt)*detJ_map = 1/4 * current_tet_volume = proportion_of_current_tet_assigned_to_this_gauss_point
      ! However we can do better than just using 1/4 (i.e. assigning it evenly).
      ! This is what gauss_tet_vol_weight is
      ! (a more careful assigning of the volume of the tetrahedron amongst the gauss points).
      ! This approach only works for tetrahedrons with 4 gauss points,
      ! as each guass point sits near a unique vertex.
      force = force + w(igpt)*4.d0*gauss_tet_vol_weight(igpt)*density* &
                      matmul(kron_N_eye,bf)*detJ_map
    else
      force = force + w(igpt)*density*matmul(kron_N_eye,bf)*detJ_map
    endif
    force = force + w(igpt)*beta*dot_product(N(:,1),local_press)* &
      vecB_transpose(:,1)*detJ_map
    force = force + w(igpt)*alpha*(3.d0*lambda+2.d0*mu)* &
      dot_product(N(:,1),local_temp)*vecB_transpose(:,1)*detJ_map
    call ShapeFunctionDestroy(shapefunction)
    deallocate(N)
    deallocate(vecB_transpose)
    deallocate(kron_B_eye)
    deallocate(kron_B_transpose_eye)
    deallocate(kron_eye_B_transpose)
    deallocate(kron_N_eye)
  enddo

  if(option%geomechanics%improve_tet_weighting .and. &
     eletype.eq.2 .and. len_w.eq.4) then
    deallocate(gauss_tet_vol_weight)
  endif

  call ConvertMatrixToVector(transpose(local_disp),vec_local_disp)
  res_vec_mat = matmul(Kmat,vec_local_disp)
  res_vec = res_vec + res_vec_mat(:,1)
  res_vec = res_vec - force

  deallocate(B)
  deallocate(force)
  deallocate(Kmat)
  deallocate(res_vec_mat)
  deallocate(vec_local_disp)
  deallocate(Trans)

end subroutine GeomechForceLocalElemResidual

! ************************************************************************** !

subroutine GeomechForceApplyTractionBCtoResidual(local_coordinates, &
                                        facetype, &
                                        stress_bc, &
                                        r,w,res_vec,option)

  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module

  PetscInt :: size_facenodes
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscInt :: facetype
  PetscReal, pointer :: r(:,:), w(:)
  PetscReal, allocatable :: res_vec(:)
  type(option_type) :: option
  PetscReal :: stress_bc(SIX_INTEGER)
  type(shapefunction_type) :: shapefunction
  PetscInt :: igpt
  PetscInt :: len_w
  PetscReal, allocatable :: x(:), J_map(:,:)
  PetscReal :: xp_J(THREE_INTEGER)
  PetscReal :: boundary_stress(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: traction(THREE_INTEGER,ONE_INTEGER)
  PetscReal :: surf_J
  PetscReal, allocatable :: N(:,:), force(:), kron_N_traction(:,:)
  PetscReal :: normal_vec(THREE_INTEGER)
  PetscInt :: shapefunc_eletype

  res_vec = 0.d0
  len_w = size(w)

  if (facetype == TRI_FACE_TYPE) then
    size_facenodes = THREE_INTEGER
    shapefunc_eletype = TRI_TYPE
  endif
  if (facetype == QUAD_FACE_TYPE) then
    size_facenodes = FOUR_INTEGER
    shapefunc_eletype = QUAD_TYPE
  endif

  allocate(force(size_facenodes*option%ngeomechdof))
  allocate(x(size_facenodes))
  allocate(J_map(size_facenodes,TWO_INTEGER))

  force = 0.d0
  boundary_stress = 0.d0

  boundary_stress(1,1) = stress_bc(1) ! sigma_xx
  boundary_stress(2,2) = stress_bc(2) ! sigma_yy
  boundary_stress(3,3) = stress_bc(3) ! sigma_zz
  boundary_stress(1,2) = stress_bc(4) ! sigma_xy
  boundary_stress(2,3) = stress_bc(5) ! sigma_yz
  boundary_stress(3,1) = stress_bc(6) ! sigma_zx

  ! symm
  boundary_stress(1,3) = boundary_stress(3,1) ! sigma_xz
  boundary_stress(3,2) = boundary_stress(2,3) ! sigma_zy
  boundary_stress(2,1) = boundary_stress(1,2) ! sigma_yx

  if (facetype == TRI_FACE_TYPE) &
    normal_vec = face_unitnormal(local_coordinates(1,:), &
                                 local_coordinates(2,:), &
                                 local_coordinates(3,:))
  if (facetype == QUAD_FACE_TYPE) &
    normal_vec = face_unitnormal(local_coordinates(1,:), &
                                 local_coordinates(2,:), &
                                 local_coordinates(4,:))

  do igpt = 1, len_w

    shapefunction%element_type = shapefunc_eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction)
    x = matmul(transpose(local_coordinates),shapefunction%N)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)

    allocate(N(size(shapefunction%N),ONE_INTEGER))

    N(:,1)= shapefunction%N
    xp_J = cross_product(J_map(:,1), J_map(:,2))
    surf_J = sqrt(dot_product(xp_J,xp_J))

    if (surf_J <= 0.d0) then
      option%io_buffer = 'GEOMECHANICS: The surface jacobian has' // &
                         ' to be positive!'
      call PrintErrMsg(option)
    endif

    traction(:,1) = matmul(boundary_stress,normal_vec)
    call Kron(N,traction,kron_N_traction)
    force = force + w(igpt)*kron_N_traction(:,1)*surf_J
    call ShapeFunctionDestroy(shapefunction)

    deallocate(N)

  enddo

  res_vec = res_vec - force

  deallocate(force)

end subroutine GeomechForceApplyTractionBCtoResidual

! ************************************************************************** !

subroutine GeomechForceApplyTractionBCtoRHS(local_coordinates, &
                                            facetype, &
                                            stress_bc, &
                                            r,w,rhs_vec,option)

  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module

  PetscInt :: size_facenodes
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscInt :: facetype
  PetscReal, pointer :: r(:,:), w(:)
  PetscReal, allocatable :: rhs_vec(:)
  type(option_type) :: option
  PetscReal :: stress_bc(SIX_INTEGER)
  type(shapefunction_type) :: shapefunction
  PetscInt :: igpt
  PetscInt :: len_w
  PetscReal, allocatable :: x(:), J_map(:,:)
  PetscReal :: xp_J(THREE_INTEGER)
  PetscReal :: boundary_stress(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: traction(THREE_INTEGER,ONE_INTEGER)
  PetscReal :: surf_J
  PetscReal, allocatable :: N(:,:), force(:), kron_N_traction(:,:)
  PetscReal :: normal_vec(THREE_INTEGER)
  PetscInt :: shapefunc_eletype

  rhs_vec = 0.d0
  len_w = size(w)

  if (facetype == TRI_FACE_TYPE) then
    size_facenodes = THREE_INTEGER
    shapefunc_eletype = TRI_TYPE
  endif
  if (facetype == QUAD_FACE_TYPE) then
    size_facenodes = FOUR_INTEGER
    shapefunc_eletype = QUAD_TYPE
  endif

  allocate(force(size_facenodes*option%ngeomechdof))
  allocate(x(size_facenodes))
  allocate(J_map(size_facenodes,TWO_INTEGER))

  force = 0.d0
  boundary_stress = 0.d0

  boundary_stress(1,1) = stress_bc(1) ! sigma_xx
  boundary_stress(2,2) = stress_bc(2) ! sigma_yy
  boundary_stress(3,3) = stress_bc(3) ! sigma_zz
  boundary_stress(1,2) = stress_bc(4) ! sigma_xy
  boundary_stress(2,3) = stress_bc(5) ! sigma_yz
  boundary_stress(3,1) = stress_bc(6) ! sigma_zx

  ! symm
  boundary_stress(1,3) = boundary_stress(3,1) ! sigma_xz
  boundary_stress(3,2) = boundary_stress(2,3) ! sigma_zy
  boundary_stress(2,1) = boundary_stress(1,2) ! sigma_yx

  if (facetype == TRI_FACE_TYPE) &
    normal_vec = face_unitnormal(local_coordinates(1,:), &
                                 local_coordinates(2,:), &
                                 local_coordinates(3,:))
  if (facetype == QUAD_FACE_TYPE) &
    normal_vec = face_unitnormal(local_coordinates(1,:), &
                                 local_coordinates(2,:), &
                                 local_coordinates(4,:))

  do igpt = 1, len_w

    shapefunction%element_type = shapefunc_eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction)
    x = matmul(transpose(local_coordinates),shapefunction%N)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)

    allocate(N(size(shapefunction%N),ONE_INTEGER))

    N(:,1)= shapefunction%N
    xp_J = cross_product(J_map(:,1), J_map(:,2))
    surf_J = sqrt(dot_product(xp_J,xp_J))

    if (surf_J <= 0.d0) then
      option%io_buffer = 'GEOMECHANICS: The surface jacobian has' // &
                         ' to be positive!'
      call PrintErrMsg(option)
    endif

    traction(:,1) = matmul(boundary_stress,normal_vec)
    call Kron(N,traction,kron_N_traction)
    force = force + w(igpt)*kron_N_traction(:,1)*surf_J
    call ShapeFunctionDestroy(shapefunction)

    deallocate(N)

  enddo

  rhs_vec = rhs_vec + force

  deallocate(force)

end subroutine GeomechForceApplyTractionBCtoRHS

! ************************************************************************** !

subroutine GeomechForceSetupLinearSystem(A,solution,rhs,geomech_realization, &
                                         ierr)
  !
  ! Computes the Coefficient matrix and the right hand side of the system
  !
  ! Author: Satish Karra
  ! Date: 06/04/2025
  !

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Grid_Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Coupler_module
  use Option_module
  use Geomechanics_Auxiliary_module

  implicit none

  Vec :: solution
  Vec :: rhs
  Mat :: A
  class(realization_geomech_type) :: geomech_realization
  PetscErrorCode :: ierr

  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(gm_region_type), pointer :: region
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_parameter_type), pointer :: GeomechParam

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_press(:), local_temp(:)
  PetscInt, allocatable :: petsc_ids(:)
  PetscInt, allocatable :: ids(:)
  PetscReal, allocatable :: rhs_local_vec(:)
  PetscReal, pointer :: press(:), temp(:)
  PetscReal, pointer :: fluid_density(:), porosity(:)
  PetscReal, pointer :: press_init(:), temp_init(:)
  PetscReal, pointer :: fluid_density_init(:)
  PetscReal, allocatable :: beta_vec(:), alpha_vec(:)
  PetscReal, allocatable :: density_rock_vec(:), density_fluid_vec(:)
  PetscReal, allocatable :: density_bulk_vec(:)
  PetscReal, allocatable :: youngs_vec(:), poissons_vec(:)
  PetscReal, allocatable :: porosity_vec(:)
  PetscInt :: ielem, ivertex
  PetscInt :: ghosted_id
  PetscInt :: eletype
  PetscInt :: petsc_id, local_id
  PetscReal, pointer :: imech_loc_p(:)
  PetscInt :: size_elenodes, idof

  PetscInt :: facetype, nfaces
  PetscInt :: iface, num_vertices
  PetscReal :: stress_bc(SIX_INTEGER)
  PetscInt, allocatable :: face_vertices(:)

  PetscReal, pointer :: temp_youngs_modulus_p(:)
  PetscReal, pointer :: temp_poissons_ratio_p(:)
  PetscReal, pointer :: temp_density_p(:)
  PetscReal, pointer :: temp_biot_coeff_p(:)
  PetscReal, pointer :: temp_thermal_exp_coeff_p(:)

  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars
  GeomechParam => patch%geomech_aux%GeomechParam


  solution = field%disp_xx

  ! use the stored matrix
  A = field%A
  rhs = field%rhs

  call VecZeroEntries(rhs,ierr);CHKERRQ(ierr)

  ! Get pressure and temperature from subsurface
  call VecGetArrayF90(field%press_loc,press,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%temp_loc,temp,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%porosity_loc,porosity,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%fluid_density_loc,fluid_density,ierr);CHKERRQ(ierr)

  ! Get geomech properties
  if (GeomechParam%youngs_modulus_spatially_varying) then
    call VecGetArrayF90(field%youngs_modulus,temp_youngs_modulus_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%poissons_ratio_spatially_varying) then
    call VecGetArrayF90(field%poissons_ratio,temp_poissons_ratio_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%density_spatially_varying) then
    call VecGetArrayF90(field%density,temp_density_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%biot_coeff_spatially_varying) then
    call VecGetArrayF90(field%biot_coeff,temp_biot_coeff_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%thermal_exp_coeff_spatially_varying) then
    call VecGetArrayF90(field%thermal_exp_coeff,temp_thermal_exp_coeff_p,ierr);CHKERRQ(ierr)
  endif

  ! Get initial pressure and temperature
  call VecGetArrayF90(field%press_init_loc,press_init,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%temp_init_loc,temp_init,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%fluid_density_init_loc,fluid_density_init, &
                      ierr);CHKERRQ(ierr)

  ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    allocate(elenodes(grid%elem_nodes(0,ielem)))
    allocate(local_coordinates(size(elenodes),THREE_INTEGER))
    allocate(local_press(size(elenodes)))
    allocate(local_temp(size(elenodes)))
    allocate(petsc_ids(size(elenodes)))
    allocate(ids(size(elenodes)*option%ngeomechdof))
    allocate(rhs_local_vec(size(elenodes)*option%ngeomechdof))
    allocate(beta_vec(size(elenodes)))
    allocate(alpha_vec(size(elenodes)))
    allocate(density_rock_vec(size(elenodes)))
    allocate(density_fluid_vec(size(elenodes)))
    allocate(youngs_vec(size(elenodes)))
    allocate(poissons_vec(size(elenodes)))
    allocate(porosity_vec(size(elenodes)))
    allocate(density_bulk_vec(size(elenodes)))
    elenodes = grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    eletype = grid%gauss_node(ielem)%entity_type
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
    enddo
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        ids(idof + (ivertex-1)*option%ngeomechdof) = &
          (petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
      enddo
      local_press(ivertex) = press(ghosted_id) - press_init(ghosted_id)  ! p - p_0
      local_temp(ivertex) = temp(ghosted_id) - temp_init(ghosted_id)     ! T - T_0
      if (GeomechParam%thermal_exp_coeff_spatially_varying) then
        alpha_vec(ivertex) = temp_thermal_exp_coeff_p(grid%nG2L(ghosted_id))
      else
        alpha_vec(ivertex) = &
          GeomechParam%thermal_exp_coeff(int(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%biot_coeff_spatially_varying) then
        beta_vec(ivertex) = temp_biot_coeff_p(grid%nG2L(ghosted_id))
      else
        beta_vec(ivertex) = &
          GeomechParam%biot_coeff(int(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%density_spatially_varying) then
        density_rock_vec(ivertex) = temp_density_p(grid%nG2L(ghosted_id))
      else
        density_rock_vec(ivertex) = &
          GeomechParam%density(int(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%youngs_modulus_spatially_varying) then
        youngs_vec(ivertex) = temp_youngs_modulus_p(grid%nG2L(ghosted_id))
      else
        youngs_vec(ivertex) = &
          GeomechParam%youngs_modulus(int(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%poissons_ratio_spatially_varying) then
        poissons_vec(ivertex) = temp_poissons_ratio_p(grid%nG2L(ghosted_id))
      else
        poissons_vec(ivertex) = &
          GeomechParam%poissons_ratio(int(imech_loc_p(ghosted_id)))
      endif
      density_fluid_vec(ivertex) = fluid_density(ghosted_id)
      porosity_vec(ivertex) = porosity(ghosted_id)
      density_bulk_vec(ivertex) = (porosity_vec(ivertex) * &
                                   density_fluid_vec(ivertex)) + &
                                  ((1.d0 - porosity_vec(ivertex)) * &
                                   density_rock_vec(ivertex))
    enddo
    size_elenodes = size(elenodes)
    call GeomechForceLocalElemRHS(size_elenodes,local_coordinates, &
       local_press,local_temp,youngs_vec,poissons_vec, &
       density_bulk_vec,beta_vec,alpha_vec,eletype, &
       grid%gauss_node(ielem)%dim,grid%gauss_node(ielem)%r, &
       grid%gauss_node(ielem)%w,rhs_local_vec,option)
    call VecSetValues(rhs,size(ids),ids,rhs_local_vec,ADD_VALUES, &
                      ierr);CHKERRQ(ierr)

    deallocate(elenodes)
    deallocate(local_coordinates)
    deallocate(petsc_ids)
    deallocate(ids)
    deallocate(rhs_local_vec)
    deallocate(local_press)
    deallocate(local_temp)
    deallocate(beta_vec)
    deallocate(alpha_vec)
    deallocate(density_rock_vec)
    deallocate(density_fluid_vec)
    deallocate(youngs_vec)
    deallocate(poissons_vec)
    deallocate(porosity_vec)
    deallocate(density_bulk_vec)
  enddo

  call VecRestoreArrayF90(field%press_loc,press,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%temp_loc,temp,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%fluid_density_loc,fluid_density, &
                          ierr);CHKERRQ(ierr)

  if (GeomechParam%youngs_modulus_spatially_varying) then
    call VecRestoreArrayF90(field%youngs_modulus,temp_youngs_modulus_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%poissons_ratio_spatially_varying) then
    call VecRestoreArrayF90(field%poissons_ratio,temp_poissons_ratio_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%density_spatially_varying) then
    call VecRestoreArrayF90(field%density,temp_density_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%biot_coeff_spatially_varying) then
    call VecRestoreArrayF90(field%biot_coeff,temp_biot_coeff_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%thermal_exp_coeff_spatially_varying) then
    call VecRestoreArrayF90(field%thermal_exp_coeff,temp_thermal_exp_coeff_p,ierr);CHKERRQ(ierr)
  endif

  call VecRestoreArrayF90(field%press_init_loc,press_init,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%temp_init_loc,temp_init,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%fluid_density_init_loc,fluid_density_init, &
                          ierr);CHKERRQ(ierr)


  call VecAssemblyBegin(rhs,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(rhs,ierr);CHKERRQ(ierr)

  ! Find the boundary nodes with dirichlet and set the residual at those nodes
  ! to zero, later set the Jacobian to 1

  ! displacement boundary conditions
  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      ! X displacement
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_X_DOF-1, &
                             boundary_condition% &
                             geomech_aux_real_var(GEOMECH_DISP_X_DOF, &
                             ivertex),INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for displacement not available.'
            call PrintErrMsg(option)
        end select
      endif

      ! Y displacement
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_Y_DOF-1, &
                             boundary_condition% &
                             geomech_aux_real_var(GEOMECH_DISP_Y_DOF, &
                             ivertex),INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for displacement not available.'
            call PrintErrMsg(option)
        end select
      endif

      ! Z displacement
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_Z_DOF-1, &
                             boundary_condition% &
                             geomech_aux_real_var(GEOMECH_DISP_Z_DOF, &
                             ivertex),INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for displacement not available.'
            call PrintErrMsg(option)
        end select
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Need to assemby here since one cannot mix INSERT_VALUES
  ! and ADD_VALUES
  call VecAssemblyBegin(rhs,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(rhs,ierr);CHKERRQ(ierr)

  ! Force boundary conditions
  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      ! X force
      if (associated(boundary_condition%geomech_condition%force_x)) then
        select case(boundary_condition%geomech_condition%force_x%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                         (petsc_id-1)*option% &
                           ngeomechdof+GEOMECH_DISP_X_DOF-1, &
                         boundary_condition% &
                           geomech_aux_real_var(GEOMECH_DISP_X_DOF,ivertex), &
                         ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for force not available.'
            call PrintErrMsg(option)
        end select
      endif

       ! Y force
      if (associated(boundary_condition%geomech_condition%force_y)) then
        select case(boundary_condition%geomech_condition%force_y%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                         (petsc_id-1)*option% &
                           ngeomechdof+GEOMECH_DISP_Y_DOF-1, &
                         boundary_condition% &
                           geomech_aux_real_var(GEOMECH_DISP_Y_DOF,ivertex), &
                         ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for force not available.'
            call PrintErrMsg(option)

        end select
      endif

       ! Z force
      if (associated(boundary_condition%geomech_condition%force_z)) then
        select case(boundary_condition%geomech_condition%force_z%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                        (petsc_id-1)*option% &
                          ngeomechdof+GEOMECH_DISP_Z_DOF-1, &
                         boundary_condition% &
                           geomech_aux_real_var(GEOMECH_DISP_Z_DOF,ivertex), &
                         ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for force not available.'
            call PrintErrMsg(option)
        end select
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecAssemblyBegin(rhs,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(rhs,ierr);CHKERRQ(ierr)

  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    ! Traction
    if (associated(boundary_condition%geomech_condition%traction)) then
      select case(boundary_condition%geomech_condition%traction%itype)
        case(DIRICHLET_BC)
          option%io_buffer = 'Dirichlet BC for traction not available.'
          call PrintErrMsg(option)
        case(ZERO_GRADIENT_BC)
         ! do nothing
        case(NEUMANN_BC)
          stress_bc = boundary_condition%geomech_condition%traction% &
                        dataset%rarray
          nfaces = boundary_condition%region%sideset%nfaces
          do iface = 1, nfaces
            num_vertices = size(boundary_condition%region%sideset% &
              face_vertices(:,iface))
            if (num_vertices == THREE_INTEGER) then
              facetype = TRI_FACE_TYPE
            else
              facetype = QUAD_FACE_TYPE
            endif
            allocate(face_vertices(num_vertices))
            face_vertices = boundary_condition%region%sideset% &
              face_vertices(1:num_vertices,iface)
            allocate(local_coordinates(num_vertices,THREE_INTEGER))
            allocate(petsc_ids(num_vertices))
            allocate(ids(num_vertices*option%ngeomechdof))
            allocate(rhs_local_vec(num_vertices*option%ngeomechdof))
            do ivertex = 1, num_vertices
              ghosted_id = face_vertices(ivertex)
              local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = &
                grid%nodes(ghosted_id)%x
              local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = &
                grid%nodes(ghosted_id)%y
              local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = &
                grid%nodes(ghosted_id)%z
              petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
              do idof = 1, option%ngeomechdof
                ids(idof + (ivertex-1)*option%ngeomechdof) = &
                  (petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
              enddo
            enddo
            call GeomechForceApplyTractionBCtoRHS( &
                   local_coordinates, &
                   facetype, &
                   stress_bc, &
                   grid%gauss_surf_node(facetype)%r, &
                   grid%gauss_surf_node(facetype)%w, &
                   rhs_local_vec,option)
            call VecSetValues(rhs,size(ids),ids,rhs_local_vec,ADD_VALUES, &
                   ierr);CHKERRQ(ierr)
            deallocate(face_vertices)
            deallocate(local_coordinates)
            deallocate(petsc_ids)
            deallocate(ids)
            deallocate(rhs_local_vec)
          enddo
      end select
    endif
    boundary_condition => boundary_condition%next
  enddo
  call VecAssemblyBegin(rhs,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(rhs,ierr);CHKERRQ(ierr)

end subroutine GeomechForceSetupLinearSystem

! ************************************************************************** !

subroutine GeomechForceLocalElemRHS(size_elenodes,local_coordinates, &
                                         local_press,local_temp, &
                                         local_youngs,local_poissons, &
                                         local_density,local_beta, &
                                         local_alpha, &
                                         eletype,dim,r,w,rhs_vec,option)
  !
  ! Computes the RHS for a local element
  !
  ! Author: Satish Karra
  ! Date: 06/24/13
  !

  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module

  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: B(:,:)
  PetscReal, allocatable :: rhs_vec(:)
  PetscReal, allocatable :: local_press(:)
  PetscReal, allocatable :: local_temp(:)
  PetscReal, allocatable :: local_youngs(:)
  PetscReal, allocatable :: local_poissons(:)
  PetscReal, allocatable :: local_density(:)
  PetscReal, allocatable :: local_beta(:)
  PetscReal, allocatable :: local_alpha(:)

  PetscReal, pointer :: r(:,:), w(:)
  PetscInt :: igpt
  PetscInt :: len_w
  PetscInt :: eletype
  PetscReal :: x(THREE_INTEGER), J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: inv_J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: detJ_map
  PetscInt :: i,j
  PetscInt :: dim
  PetscReal :: lambda, mu, beta, alpha
  PetscReal :: density, youngs_mod, poissons_ratio
  PetscInt :: load_type
  PetscReal :: bf(THREE_INTEGER)
  PetscReal :: identity(THREE_INTEGER,THREE_INTEGER)
  PetscReal, allocatable :: N(:,:)
  PetscReal, allocatable :: vecB_transpose(:,:)
  PetscReal, allocatable :: kron_N_eye(:,:)
  PetscInt :: size_elenodes

  allocate(B(size_elenodes,dim))

  rhs_vec = 0.d0
  len_w = size(w)

  identity = 0.d0
  do i = 1, THREE_INTEGER
    do j = 1, THREE_INTEGER
      if (i == j) identity(i,j) = 1.d0
    enddo
  enddo

  do igpt = 1, len_w
    shapefunction%element_type = eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction)
    x = matmul(transpose(local_coordinates),shapefunction%N)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)
    allocate(N(size(shapefunction%N),ONE_INTEGER))
    call Determinant(J_map,detJ_map)
    if (detJ_map <= 0.d0) then
      option%io_buffer = 'GEOMECHANICS: Determinant of J_map has' // &
                         ' to be positive!'
      call PrintErrMsg(option)
    endif
    ! Find the inverse of J_map
    call MatInv3(J_map,inv_J_map)
    B = matmul(shapefunction%DN,inv_J_map)
    youngs_mod = dot_product(shapefunction%N,local_youngs)
    poissons_ratio = dot_product(shapefunction%N,local_poissons)
    alpha = dot_product(shapefunction%N,local_alpha)
    beta = dot_product(shapefunction%N,local_beta)
    density = dot_product(shapefunction%N,local_density)
    call GeomechGetLambdaMu(lambda,mu,youngs_mod,poissons_ratio)
    call GeomechGetBodyForce(load_type,lambda,mu,x,bf,option)
    call ConvertMatrixToVector(transpose(B),vecB_transpose)
    N(:,1)= shapefunction%N
    call Kron(N,identity,kron_N_eye)
    rhs_vec = rhs_vec + w(igpt)*density*matmul(kron_N_eye,bf)*detJ_map
    rhs_vec = rhs_vec + w(igpt)*beta*dot_product(N(:,1),local_press)* &
      vecB_transpose(:,1)*detJ_map
    rhs_vec = rhs_vec + w(igpt)*alpha*(3*lambda+2*mu)* &
      dot_product(N(:,1),local_temp)*vecB_transpose(:,1)*detJ_map
    call ShapeFunctionDestroy(shapefunction)
    deallocate(N)
    deallocate(vecB_transpose)
    deallocate(kron_N_eye)
  enddo

  deallocate(B)

end subroutine GeomechForceLocalElemRHS

! ************************************************************************** !

subroutine GeomechForceAssembleCoeffMatrixLocal(size_elenodes,local_coordinates, &
                                                local_disp, &
                                                local_youngs,local_poissons, &
                                                eletype,dim,r,w,Kmat,option)
  !
  ! Computes the Coefficient matrix locally of an element
  !
  ! Author: Satish Karra
  ! Date: 06/24/13
  ! Modified: 06/12/25
  !

  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module

  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: B(:,:), Kmat(:,:)
  PetscReal, allocatable :: local_disp(:)
  PetscReal, pointer :: r(:,:), w(:)
  PetscInt :: igpt
  PetscInt :: len_w
  PetscInt :: eletype
  PetscReal :: x(THREE_INTEGER), J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: inv_J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: detJ_map
  PetscInt :: i,j
  PetscInt :: dim
  PetscReal :: lambda, mu
  PetscReal :: youngs_mod, poissons_ratio
  PetscReal :: identity(THREE_INTEGER,THREE_INTEGER)
  PetscReal, allocatable :: N(:,:)
  PetscReal, allocatable :: vecB_transpose(:,:)
  PetscReal, allocatable :: kron_B_eye(:,:)
  PetscReal, allocatable :: kron_B_transpose_eye(:,:)
  PetscReal, allocatable :: Trans(:,:)
  PetscReal, allocatable :: kron_eye_B_transpose(:,:)
  PetscReal, allocatable :: kron_N_eye(:,:)
  PetscReal, allocatable :: local_youngs(:)
  PetscReal, allocatable :: local_poissons(:)
  PetscInt :: size_elenodes

  allocate(B(size_elenodes,dim))

  Kmat = 0.d0
  len_w = size(w)

  identity = 0.d0
  do i = 1, THREE_INTEGER
    do j = 1, THREE_INTEGER
      if (i == j) identity(i,j) = 1.d0
    enddo
  enddo

  call Transposer(option%ngeomechdof,size_elenodes,Trans)

  do igpt = 1, len_w
    shapefunction%element_type = eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction)
    x = matmul(transpose(local_coordinates),shapefunction%N)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)
    allocate(N(size(shapefunction%N),ONE_INTEGER))
    call Determinant(J_map,detJ_map)
    if (detJ_map <= 0.d0) then
      option%io_buffer = 'GEOMECHANICS: Determinant of J_map has' // &
                         ' to be positive!'
      call PrintErrMsg(option)
    endif
    ! Find the inverse of J_map
    call MatInv3(J_map,inv_J_map)
    B = matmul(shapefunction%DN,inv_J_map)
    youngs_mod = dot_product(shapefunction%N,local_youngs)
    poissons_ratio = dot_product(shapefunction%N,local_poissons)
    call GeomechGetLambdaMu(lambda,mu,youngs_mod,poissons_ratio)
    call ConvertMatrixToVector(transpose(B),vecB_transpose)
    Kmat = Kmat + w(igpt)*lambda* &
      matmul(vecB_transpose,transpose(vecB_transpose))*detJ_map
    call Kron(B,identity,kron_B_eye)
    call Kron(transpose(B),identity,kron_B_transpose_eye)
    call Kron(identity,transpose(B),kron_eye_B_transpose)
    N(:,1)= shapefunction%N
    call Kron(N,identity,kron_N_eye)
    Kmat = Kmat + w(igpt)*mu* &
      matmul(kron_B_eye,kron_B_transpose_eye)*detJ_map
    Kmat = Kmat + w(igpt)*mu* &
      matmul(matmul(kron_B_eye,kron_eye_B_transpose),Trans)*detJ_map
    call ShapeFunctionDestroy(shapefunction)
    deallocate(N)
    deallocate(vecB_transpose)
    deallocate(kron_B_eye)
    deallocate(kron_B_transpose_eye)
    deallocate(kron_eye_B_transpose)
    deallocate(kron_N_eye)
  enddo

  deallocate(B)
  deallocate(Trans)

end subroutine GeomechForceAssembleCoeffMatrixLocal

! ************************************************************************** !

subroutine GeomechGetLambdaMu(lambda,mu,E,nu)
  !
  ! Gets the material properties given the position
  ! of the point
  !
  ! Author: Satish Karra
  ! Date: 06/24/13
  !

  PetscReal :: lambda, mu
  PetscReal :: E, nu

  lambda = E*nu/(1.d0+nu)/(1.d0-2.d0*nu)
  mu = E/2.d0/(1.d0+nu)


end subroutine GeomechGetLambdaMu

! ************************************************************************** !

subroutine GeomechGetBodyForce(load_type,lambda,mu,coord,bf,option)
  !
  ! Gets the body force at a given position
  ! of the point
  !
  ! Author: Satish Karra
  ! Date: 06/24/13
  !

  use Option_module

  type(option_type) :: option

  PetscInt :: load_type
  PetscReal :: lambda, mu
  PetscReal :: coord(THREE_INTEGER)
  PetscReal :: bf(THREE_INTEGER)
  PetscReal :: x, y, z

  bf = 0.d0

  x = coord(1)
  y = coord(2)
  z = coord(3)

  ! This subroutine needs major changes. For given position, it needs to give
  ! out lambda, mu, also need to add density of rock


  select case(load_type)
    case default
      bf(GEOMECH_DISP_X_DOF) = option%geomechanics%gravity(X_DIRECTION)
      bf(GEOMECH_DISP_Y_DOF) = option%geomechanics%gravity(Y_DIRECTION)
      bf(GEOMECH_DISP_Z_DOF) = option%geomechanics%gravity(Z_DIRECTION)
  end select

end subroutine GeomechGetBodyForce

! ************************************************************************** !

subroutine GeomechForceAssembleCoeffMatrix(A,geomech_realization)
  !
  ! Computes the Assembled Coefficient Matrix
  !
  ! Author: Satish Karra
  ! Date: 06/21/13
  ! Modified: 07/12/16, 06/12/25
  !

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Coupler_module
  use Geomechanics_Field_module
  use Geomechanics_Debug_module
  use Geomechanics_Discretization_module
  use Option_module
  use Grid_Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Auxiliary_module

  implicit none

  Mat :: A

  PetscErrorCode :: ierr

  class(realization_geomech_type) :: geomech_realization
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(gm_region_type), pointer :: region
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_parameter_type), pointer :: GeomechParam

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_disp(:)
  PetscInt, allocatable :: ghosted_ids(:)
  PetscReal, allocatable :: Jac_full(:,:)
  ! due to PETSc explicit interface declaring the sub matrix as pointer,
  ! must be pointer and not allocatable
  PetscReal, pointer :: Jac_sub_mat(:,:)
  PetscInt, allocatable :: rows(:)
  PetscReal, allocatable :: youngs_vec(:), poissons_vec(:)
  PetscInt :: ielem,ivertex
  PetscInt :: ghosted_id
  PetscInt :: eletype, idof
  PetscInt :: local_id, petsc_id
  PetscInt :: ghosted_id1, ghosted_id2
  PetscInt :: petsc_id1, petsc_id2
  PetscInt :: id1, id2, vertex_count, count
  PetscReal, pointer :: imech_loc_p(:)
  PetscInt :: size_elenodes

  PetscReal, pointer :: temp_youngs_modulus_p(:)
  PetscReal, pointer :: temp_poissons_ratio_p(:)

  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars
  GeomechParam => patch%geomech_aux%GeomechParam

  call MatZeroEntries(A,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)

  if (GeomechParam%youngs_modulus_spatially_varying) then
    call VecGetArrayF90(field%youngs_modulus,temp_youngs_modulus_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%poissons_ratio_spatially_varying) then
    call VecGetArrayF90(field%poissons_ratio,temp_poissons_ratio_p,ierr);CHKERRQ(ierr)
  endif

  ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    allocate(elenodes(grid%elem_nodes(0,ielem)))
    allocate(local_coordinates(size(elenodes),THREE_INTEGER))
    allocate(local_disp(size(elenodes)*option%ngeomechdof))
    allocate(ghosted_ids(size(elenodes)))
    allocate(Jac_full(size(elenodes)*option%ngeomechdof, &
                      size(elenodes)*option%ngeomechdof))
    allocate(Jac_sub_mat(option%ngeomechdof,option%ngeomechdof))
    allocate(youngs_vec(size(elenodes)))
    allocate(poissons_vec(size(elenodes)))
    elenodes = grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    eletype = grid%gauss_node(ielem)%entity_type
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      ghosted_ids(ivertex) = ghosted_id
    enddo
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        local_disp(idof + (ivertex-1)*option%ngeomechdof) = &
          geomech_global_aux_vars(ghosted_id)%disp_vector(idof)
      enddo
      if (GeomechParam%youngs_modulus_spatially_varying) then
        youngs_vec(ivertex) = temp_youngs_modulus_p(grid%nG2L(ghosted_id))
      else
        youngs_vec(ivertex) = &
          GeomechParam%youngs_modulus(nint(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%poissons_ratio_spatially_varying) then
        poissons_vec(ivertex) = temp_poissons_ratio_p(grid%nG2L(ghosted_id))
      else
        poissons_vec(ivertex) = &
          GeomechParam%poissons_ratio(nint(imech_loc_p(ghosted_id)))
      endif
    enddo
    size_elenodes = size(elenodes)
    call GeomechForceAssembleCoeffMatrixLocal(size_elenodes,local_coordinates, &
       local_disp,youngs_vec,poissons_vec,eletype, &
       grid%gauss_node(ielem)%dim,grid%gauss_node(ielem)%r, &
       grid%gauss_node(ielem)%w,Jac_full,option)
    do id1 = 1, size(ghosted_ids)
      ghosted_id1 = ghosted_ids(id1)
      petsc_id1 = grid%node_ids_ghosted_petsc(ghosted_id1)
      do id2 = 1, size(ghosted_ids)
        ghosted_id2 = ghosted_ids(id2)
        petsc_id2 = grid%node_ids_ghosted_petsc(ghosted_id2)
        Jac_sub_mat = 0.d0
        Jac_sub_mat =  &
          Jac_full(option%ngeomechdof*(id1-1)+GEOMECH_DISP_X_DOF: &
                   option%ngeomechdof*(id1-1)+GEOMECH_DISP_Z_DOF, &
                   option%ngeomechdof*(id2-1)+GEOMECH_DISP_X_DOF: &
                   option%ngeomechdof*(id2-1)+GEOMECH_DISP_Z_DOF)

        call MatSetValuesBlocked(A,1,petsc_id1-1,1,petsc_id2-1,Jac_sub_mat, &
                                 ADD_VALUES,ierr);CHKERRQ(ierr)
      enddo
    enddo

    deallocate(elenodes)
    deallocate(local_coordinates)
    deallocate(local_disp)
    deallocate(ghosted_ids)
    deallocate(Jac_full)
    deallocate(Jac_sub_mat)
    deallocate(youngs_vec)
    deallocate(poissons_vec)
  enddo

  call VecRestoreArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)

  if (GeomechParam%youngs_modulus_spatially_varying) then
    call VecRestoreArrayF90(field%youngs_modulus,temp_youngs_modulus_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%poissons_ratio_spatially_varying) then
    call VecRestoreArrayF90(field%poissons_ratio,temp_poissons_ratio_p,ierr);CHKERRQ(ierr)
  endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! Find the boundary nodes with dirichlet and set the residual at those nodes
  ! to zero, later set the Jacobian to 1

  ! Find the number of boundary vertices
  vertex_count = 0
  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    vertex_count = vertex_count + region%num_verts
    boundary_condition => boundary_condition%next
  enddo

  allocate(rows(vertex_count*option%ngeomechdof))
  count = 0

  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      ! X displacement
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            count = count + 1
            rows(count) = (ghosted_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_X_DOF-1
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

      ! Y displacement
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            count = count + 1
            rows(count) = (ghosted_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Y_DOF-1
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

      ! Z displacement
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            count = count + 1
            rows(count) = (ghosted_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Z_DOF-1
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call MatZeroRowsLocal(A,count,rows,1.d0,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                        ierr);CHKERRQ(ierr)
  call MatSetOption(A,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE, &
                    ierr);CHKERRQ(ierr)
  call MatStoreValues(A,ierr);CHKERRQ(ierr)

  deallocate(rows)

end subroutine GeomechForceAssembleCoeffMatrix

! ************************************************************************** !

subroutine GeomechUpdateFromSubsurf(realization,geomech_realization)
  !
  ! The pressure/temperature from subsurface are
  ! mapped to geomech
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/10/13
  !

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization
  type(grid_type), pointer :: grid
  type(geomech_grid_type), pointer :: geomech_grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(geomech_field_type), pointer :: geomech_field
  type(gmdm_ptr_type), pointer :: dm_ptr

  PetscErrorCode :: ierr
  PetscReal, pointer :: vec_p(:), xx_loc_p(:)
  PetscInt :: local_id, ghosted_id

  option        => realization%option
  grid          => realization%discretization%grid
  field         => realization%field
  geomech_grid  => geomech_realization%geomech_discretization%grid
  geomech_field => geomech_realization%geomech_field

  ! use the subsurface output option parameters for geomechanics as well
  geomech_realization%output_option%tunit = realization%output_option%tunit
  geomech_realization%output_option%tconv = realization%output_option%tconv

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_realization% &
                                                   geomech_discretization, &
                                                   ONEDOF)


  ! pressure
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(geomech_field%subsurf_vec_1dof,vec_p,ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    vec_p(local_id) = xx_loc_p(option%nflowdof*(ghosted_id-1)+1)
  enddo
  call VecRestoreArrayF90(geomech_field%subsurf_vec_1dof,vec_p,ierr)
  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  ! Scatter the data
  call VecScatterBegin(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                       geomech_field%subsurf_vec_1dof,geomech_field%press, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                     geomech_field%subsurf_vec_1dof,geomech_field%press, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  ! temperature
  if (option%nflowdof > 1) then
    call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(geomech_field%subsurf_vec_1dof,vec_p,ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      vec_p(local_id) = xx_loc_p(option%nflowdof*(ghosted_id-1)+2)
    enddo
    call VecRestoreArrayF90(geomech_field%subsurf_vec_1dof,vec_p,ierr)
    call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

    ! Scatter the data
    call VecScatterBegin(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                         geomech_field%subsurf_vec_1dof,geomech_field%temp, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                       geomech_field%subsurf_vec_1dof,geomech_field%temp, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  endif

  call GeomechDiscretizationGlobalToLocal(&
                                geomech_realization%geomech_discretization, &
                                geomech_field%press, &
                                geomech_field%press_loc,ONEDOF)

  if (option%nflowdof > 1) &
    call GeomechDiscretizationGlobalToLocal(&
                                geomech_realization%geomech_discretization, &
                                geomech_field%temp, &
                                geomech_field%temp_loc,ONEDOF)

end subroutine GeomechUpdateFromSubsurf

! ************************************************************************** !

subroutine GeomechUpdateSubsurfFromGeomech(realization,geomech_realization)
  !
  ! The stresses and strains from geomech
  ! are mapped to subsurf.
  !
  ! Author: Satish Karra, LANL
  ! Date: 10/10/13
  !

  use Realization_Subsurface_class
  use Discretization_module
  use Grid_module
  use Field_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization
  type(grid_type), pointer :: grid
  type(geomech_grid_type), pointer :: geomech_grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(geomech_field_type), pointer :: geomech_field
  type(gmdm_ptr_type), pointer :: dm_ptr

  PetscErrorCode :: ierr

  option        => realization%option
  grid          => realization%discretization%grid
  field         => realization%field
  geomech_grid  => geomech_realization%geomech_discretization%grid
  geomech_field => geomech_realization%geomech_field

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_realization% &
                                                   geomech_discretization, &
                                                   ONEDOF)

  ! Scatter the strains
  call VecScatterBegin(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                       geomech_field%strain,geomech_field%strain_subsurf, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                     geomech_field%strain,geomech_field%strain_subsurf, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  ! Scatter the stresses
  call VecScatterBegin(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                       geomech_field%stress,geomech_field%stress_subsurf, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                     geomech_field%stress,geomech_field%stress_subsurf, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  ! Scatter from global to local vectors
  call DiscretizationGlobalToLocal(realization%discretization, &
                                   geomech_field%strain_subsurf, &
                                   geomech_field%strain_subsurf_loc, &
                                   NGEODOF)
  call DiscretizationGlobalToLocal(realization%discretization, &
                                   geomech_field%stress_subsurf, &
                                   geomech_field%stress_subsurf_loc, &
                                   NGEODOF)

end subroutine GeomechUpdateSubsurfFromGeomech

! ************************************************************************** !

subroutine GeomechCreateGeomechSubsurfVec(realization,geomech_realization)
  !
  ! Creates the MPI vector that stores the
  ! variables from subsurface
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/10/13
  !

  use Grid_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Field_module
  use String_module
  use Realization_Subsurface_class
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: geomech_field

  PetscErrorCode :: ierr

  option     => realization%option
  grid       => realization%discretization%grid
  geomech_field => geomech_realization%geomech_field

  call VecCreate(option%mycomm,geomech_field%subsurf_vec_1dof, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%subsurf_vec_1dof,grid%nlmax,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%subsurf_vec_1dof,ierr);CHKERRQ(ierr)

end subroutine GeomechCreateGeomechSubsurfVec

! ************************************************************************** !

subroutine GeomechCreateSubsurfStressStrainVec(realization,geomech_realization)
  !
  ! Creates the subsurface stress and strain
  ! MPI vectors to store information from geomechanics
  !
  ! Author: Satish Karra, LANL
  ! Date: 10/10/13
  !

  use Grid_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Field_module
  use String_module
  use Realization_Subsurface_class
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: geomech_field

  PetscErrorCode :: ierr

  option     => realization%option
  grid       => realization%discretization%grid
  geomech_field => geomech_realization%geomech_field

  ! strain
  call VecCreate(option%mycomm,geomech_field%strain_subsurf, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%strain_subsurf,grid%nlmax*SIX_INTEGER, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%strain_subsurf,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%strain_subsurf,ierr);CHKERRQ(ierr)

  ! stress
  call VecCreate(option%mycomm,geomech_field%stress_subsurf, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%stress_subsurf,grid%nlmax*SIX_INTEGER, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%stress_subsurf,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%stress_subsurf,ierr);CHKERRQ(ierr)

  ! strain_loc
  call VecCreate(PETSC_COMM_SELF,geomech_field%strain_subsurf_loc, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%strain_subsurf_loc,grid%ngmax*SIX_INTEGER, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%strain_subsurf_loc,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%strain_subsurf_loc,ierr);CHKERRQ(ierr)

  ! stress_loc
  call VecCreate(PETSC_COMM_SELF,geomech_field%stress_subsurf_loc, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%stress_subsurf_loc,grid%ngmax*SIX_INTEGER, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%stress_subsurf_loc,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%stress_subsurf_loc,ierr);CHKERRQ(ierr)

end subroutine GeomechCreateSubsurfStressStrainVec

! ************************************************************************** !

subroutine GeomechForceStressStrain(geomech_realization)
  !
  ! Computes the stress strain on a patch
  !
  ! Author: Satish Karra
  ! Date: 09/17/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Grid_Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Coupler_module
  use Option_module
  use Geomechanics_Auxiliary_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(geomech_parameter_type), pointer :: GeomechParam

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscInt, allocatable :: petsc_ids(:)
  PetscInt, allocatable :: ids(:)
  PetscReal, allocatable :: youngs_vec(:), poissons_vec(:)
  PetscReal, allocatable :: alpha_vec(:)  ! thermal expansion coefficient
  PetscReal, allocatable :: beta_vec(:)  ! Biot's coefficient
  PetscReal, allocatable :: local_temp(:)
  PetscReal, allocatable :: local_press(:)
  PetscReal, allocatable :: strain(:,:), stress(:,:), stress_total(:,:)
  PetscInt :: ielem, ivertex
  PetscInt :: ghosted_id
  PetscInt :: eletype, idof
  PetscInt :: local_id
  PetscInt :: size_elenodes
  PetscReal, pointer :: imech_loc_p(:)
  PetscReal, pointer :: strain_loc_p(:)
  PetscReal, pointer :: stress_loc_p(:)
  PetscReal, pointer :: stress_total_loc_p(:)
  PetscReal, pointer :: temp_loc_p(:), temp_init_loc_p(:)
  PetscReal, pointer :: press_loc_p(:), press_init_loc_p(:)
  PetscReal, pointer :: strain_p(:), stress_p(:), stress_total_p(:)
  PetscReal, pointer :: no_elems_p(:)

  PetscReal, pointer :: temp_youngs_modulus_p(:)
  PetscReal, pointer :: temp_poissons_ratio_p(:)
  PetscReal, pointer :: temp_biot_coeff_p(:)
  PetscReal, pointer :: temp_thermal_exp_coeff_p(:)

  PetscErrorCode :: ierr

  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars
  GeomechParam => patch%geomech_aux%GeomechParam

  call VecSet(field%strain,0.d0,ierr);CHKERRQ(ierr)
  call VecSet(field%stress,0.d0,ierr);CHKERRQ(ierr)
  call VecSet(field%stress_total,0.d0,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%stress_total_loc,stress_total_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%temp_loc,temp_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%temp_init_loc,temp_init_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%press_loc,press_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%press_init_loc,press_init_loc_p,ierr);CHKERRQ(ierr)

  if (GeomechParam%youngs_modulus_spatially_varying) then
    call VecGetArrayF90(field%youngs_modulus,temp_youngs_modulus_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%poissons_ratio_spatially_varying) then
    call VecGetArrayF90(field%poissons_ratio,temp_poissons_ratio_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%biot_coeff_spatially_varying) then
    call VecGetArrayF90(field%biot_coeff,temp_biot_coeff_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%thermal_exp_coeff_spatially_varying) then
    call VecGetArrayF90(field%thermal_exp_coeff,temp_thermal_exp_coeff_p,ierr);CHKERRQ(ierr)
  endif

  strain_loc_p = 0.d0
  stress_loc_p = 0.d0
  stress_total_loc_p = 0.d0

   ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    allocate(elenodes(grid%elem_nodes(0,ielem)))
    allocate(local_coordinates(size(elenodes),THREE_INTEGER))
    allocate(local_disp(size(elenodes),option%ngeomechdof))
    allocate(petsc_ids(size(elenodes)))
    allocate(local_temp(size(elenodes)))
    allocate(local_press(size(elenodes)))
    allocate(ids(size(elenodes)*option%ngeomechdof))
    allocate(youngs_vec(size(elenodes)))
    allocate(poissons_vec(size(elenodes)))
    allocate(alpha_vec(size(elenodes)))
    allocate(beta_vec(size(elenodes)))
    allocate(strain(size(elenodes),SIX_INTEGER))
    allocate(stress(size(elenodes),SIX_INTEGER))
    allocate(stress_total(size(elenodes),SIX_INTEGER))
    elenodes = grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    eletype = grid%gauss_node(ielem)%entity_type
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
    enddo
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        local_disp(ivertex,idof) = &
          geomech_global_aux_vars(ghosted_id)%disp_vector(idof)
        ids(idof + (ivertex-1)*option%ngeomechdof) = &
          (petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
      enddo
      local_temp(ivertex) = &
        temp_loc_p(ghosted_id) - temp_init_loc_p(ghosted_id)
      local_press(ivertex) = &
        press_loc_p(ghosted_id) - press_init_loc_p(ghosted_id)
      if (GeomechParam%youngs_modulus_spatially_varying) then
        youngs_vec(ivertex) = temp_youngs_modulus_p(grid%nG2L(ghosted_id))
      else
        youngs_vec(ivertex) = &
          GeomechParam%youngs_modulus(nint(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%poissons_ratio_spatially_varying) then
        poissons_vec(ivertex) = temp_poissons_ratio_p(grid%nG2L(ghosted_id))
      else
        poissons_vec(ivertex) = &
          GeomechParam%poissons_ratio(nint(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%biot_coeff_spatially_varying) then
        beta_vec(ivertex) = temp_biot_coeff_p(grid%nG2L(ghosted_id))
      else
        beta_vec(ivertex) = &
          GeomechParam%biot_coeff(nint(imech_loc_p(ghosted_id)))
      endif
      if (GeomechParam%thermal_exp_coeff_spatially_varying) then
        alpha_vec(ivertex) = temp_thermal_exp_coeff_p(grid%nG2L(ghosted_id))
      else
        alpha_vec(ivertex) = &
          GeomechParam%thermal_exp_coeff(nint(imech_loc_p(ghosted_id)))
      endif
    enddo
    size_elenodes = size(elenodes)
    call GeomechForceLocalElemStressStrain(size_elenodes,local_coordinates, &
       local_disp,local_temp,local_press,youngs_vec, &
       poissons_vec,alpha_vec,beta_vec,eletype,grid%gauss_node(ielem)%dim, &
       strain,stress,PETSC_FALSE,option)
    call GeomechForceLocalElemStressStrain(size_elenodes,local_coordinates, &
       local_disp,local_temp,local_press,youngs_vec, &
       poissons_vec,alpha_vec,beta_vec,eletype,grid%gauss_node(ielem)%dim, &
       strain,stress_total,PETSC_TRUE,option)

    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, SIX_INTEGER
        strain_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) = &
          strain_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) + &
          strain(ivertex,idof)
        stress_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) = &
          stress_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) + &
          stress(ivertex,idof)
        stress_total_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) = &
          stress_total_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) + &
          stress_total(ivertex,idof)
      enddo
    enddo

    deallocate(elenodes)
    deallocate(local_coordinates)
    deallocate(local_disp)
    deallocate(local_temp)
    deallocate(local_press)
    deallocate(petsc_ids)
    deallocate(ids)
    deallocate(youngs_vec)
    deallocate(poissons_vec)
    deallocate(alpha_vec)
    deallocate(beta_vec)
    deallocate(strain)
    deallocate(stress)
    deallocate(stress_total)
  enddo

  call VecRestoreArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%stress_total_loc,stress_total_loc_p,ierr);CHKERRQ(ierr)

  if (GeomechParam%youngs_modulus_spatially_varying) then
    call VecRestoreArrayF90(field%youngs_modulus,temp_youngs_modulus_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%poissons_ratio_spatially_varying) then
    call VecRestoreArrayF90(field%poissons_ratio,temp_poissons_ratio_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%biot_coeff_spatially_varying) then
    call VecRestoreArrayF90(field%biot_coeff,temp_biot_coeff_p,ierr);CHKERRQ(ierr)
  endif
  if (GeomechParam%thermal_exp_coeff_spatially_varying) then
    call VecRestoreArrayF90(field%thermal_exp_coeff,temp_thermal_exp_coeff_p,ierr);CHKERRQ(ierr)
  endif

  call GeomechDiscretizationLocalToGlobalAdd(geomech_discretization, &
                                             field%strain_loc,field%strain, &
                                             SIX_INTEGER)
  call GeomechDiscretizationLocalToGlobalAdd(geomech_discretization, &
                                             field%stress_loc,field%stress, &
                                             SIX_INTEGER)
  call GeomechDiscretizationLocalToGlobalAdd(geomech_discretization, &
                                             field%stress_total_loc, &
                                             field%stress_total, &
                                             SIX_INTEGER)

! Now take the average at each node for elements sharing the node
  call VecGetArrayF90(grid%no_elems_sharing_node,no_elems_p, &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%strain,strain_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%stress,stress_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%stress_total,stress_total_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax_node
    do idof = 1, SIX_INTEGER
      strain_p(idof + (local_id-1)*SIX_INTEGER) = &
        strain_p(idof + (local_id-1)*SIX_INTEGER)/nint(no_elems_p(local_id))
      stress_p(idof + (local_id-1)*SIX_INTEGER) = &
        stress_p(idof + (local_id-1)*SIX_INTEGER)/nint(no_elems_p(local_id))
      stress_total_p(idof + (local_id-1)*SIX_INTEGER) = &
        stress_total_p(idof + (local_id-1)*SIX_INTEGER)/nint(no_elems_p(local_id))
    enddo
  enddo
  call VecRestoreArrayF90(field%strain,strain_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%stress,stress_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%stress_total,stress_total_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%no_elems_sharing_node,no_elems_p, &
                          ierr);CHKERRQ(ierr)

! Now scatter back to local domains
  call GeomechDiscretizationGlobalToLocal(geomech_discretization, &
                                          field%strain,field%strain_loc, &
                                          SIX_INTEGER)
  call GeomechDiscretizationGlobalToLocal(geomech_discretization, &
                                          field%stress,field%stress_loc, &
                                          SIX_INTEGER)
  call GeomechDiscretizationGlobalToLocal(geomech_discretization, &
                                          field%stress_total, &
                                          field%stress_total_loc, &
                                          SIX_INTEGER)

  call VecGetArrayF90(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%stress_total_loc,stress_total_loc_p,ierr);CHKERRQ(ierr)
! Copy them to global_aux_vars
  do ghosted_id = 1, grid%ngmax_node
    do idof = 1, SIX_INTEGER
      geomech_global_aux_vars(ghosted_id)%strain(idof) = &
        strain_loc_p(idof + (ghosted_id-1)*SIX_INTEGER)
      geomech_global_aux_vars(ghosted_id)%stress(idof) = &
        stress_loc_p(idof + (ghosted_id-1)*SIX_INTEGER)
      geomech_global_aux_vars(ghosted_id)%stress_total(idof) = &
        stress_total_loc_p(idof + (ghosted_id-1)*SIX_INTEGER)
    enddo
  enddo
  call VecRestoreArrayF90(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%stress_total_loc,stress_total_loc_p,ierr);CHKERRQ(ierr)

end subroutine GeomechForceStressStrain

! ************************************************************************** !

subroutine GeomechForceLocalElemStressStrain(size_elenodes,local_coordinates, &
                                             local_disp,local_temp, &
                                             local_press, &
                                             local_youngs,local_poissons, &
                                             local_alpha, local_beta, &
                                             eletype,dim,strain,stress, &
                                             compute_stress_total,option)
  !
  ! Computes the stress-strain for a local
  ! element
  !
  ! Author: Satish Karra
  ! Date: 09/17/13
  !

  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module

  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: B(:,:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscReal, allocatable :: local_temp(:)
  PetscReal, allocatable :: local_press(:)
  PetscReal, allocatable :: local_youngs(:)
  PetscReal, allocatable :: local_poissons(:)
  PetscReal, allocatable :: local_alpha(:)  ! Thermal expansion coefficient
  PetscReal, allocatable :: local_beta(:)  ! Biot's coefficient
  PetscReal, allocatable :: strain(:,:)
  PetscReal, allocatable :: stress(:,:)
  PetscReal :: strain_local(NINE_INTEGER,ONE_INTEGER)
  PetscReal :: stress_local(NINE_INTEGER,ONE_INTEGER)
  PetscBool :: compute_stress_total

  PetscInt :: ivertex
  PetscInt :: eletype
  PetscReal :: identity(THREE_INTEGER,THREE_INTEGER)
  PetscInt :: dim
  PetscInt :: i, j
  PetscReal :: lambda, mu
  PetscReal :: youngs_mod, poissons_ratio, alpha, beta
  PetscReal, allocatable :: kron_B_eye(:,:)
  PetscReal, allocatable :: kron_B_transpose_eye(:,:)
  PetscReal, allocatable :: Trans(:,:)
  PetscReal, allocatable :: kron_eye_B_transpose(:,:)
  PetscReal, allocatable :: vec_local_disp(:,:)
  PetscInt :: size_elenodes
  PetscReal :: J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: inv_J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: eye_vec(NINE_INTEGER,ONE_INTEGER)
  PetscReal :: dT, dP

  allocate(B(size_elenodes,dim))

  call Transposer(option%ngeomechdof,size_elenodes,Trans)
  strain = 0.d0
  stress = 0.d0

  call ConvertMatrixToVector(transpose(local_disp),vec_local_disp)

  identity = 0.d0
  do i = 1, THREE_INTEGER
    do j = 1, THREE_INTEGER
      if (i == j) identity(i,j) = 1.d0
    enddo
  enddo

  eye_vec = 0.d0
  eye_vec(1,1) = 1.d0
  eye_vec(5,1) = 1.d0
  eye_vec(9,1) = 1.d0

  do ivertex = 1, size_elenodes
    strain_local = 0.d0
    stress_local = 0.d0
    shapefunction%element_type = eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = shapefunction%coord(ivertex,:)
    call ShapeFunctionCalculate(shapefunction)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)
    call MatInv3(J_map,inv_J_map)
    B = matmul(shapefunction%DN,inv_J_map)
    youngs_mod = dot_product(shapefunction%N,local_youngs)
    poissons_ratio = dot_product(shapefunction%N,local_poissons)
    alpha = dot_product(shapefunction%N,local_alpha)
    beta = dot_product(shapefunction%N,local_beta)
    call GeomechGetLambdaMu(lambda,mu,youngs_mod,poissons_ratio)
    call Kron(B,identity,kron_B_eye)
    call Kron(transpose(B),identity,kron_B_transpose_eye)
    call Kron(identity,transpose(B),kron_eye_B_transpose)
    strain_local =  0.5d0*matmul((kron_B_transpose_eye + &
      matmul(kron_eye_B_transpose,Trans)),vec_local_disp)
    stress_local = lambda*(strain_local(1,1)+ &
                   strain_local(5,1)+strain_local(9,1))*eye_vec + &
                   2.d0*mu*strain_local
    if (compute_stress_total) then
      dT = local_temp(ivertex)
      dP = local_press(ivertex)
      ! Add the thermal stress
      stress_local = stress_local - alpha * dT * eye_vec * (3.d0*lambda +2.d0*mu)
      ! Add Biot's contribution
      stress_local = stress_local - beta * dP * eye_vec
    endif
    call ShapeFunctionDestroy(shapefunction)
    deallocate(kron_B_eye)
    deallocate(kron_B_transpose_eye)
    deallocate(kron_eye_B_transpose)
    strain(ivertex,1) = strain_local(1,1)
    strain(ivertex,2) = strain_local(5,1)
    strain(ivertex,3) = strain_local(9,1)
    strain(ivertex,4) = strain_local(2,1)
    strain(ivertex,5) = strain_local(6,1)
    strain(ivertex,6) = strain_local(3,1)
    stress(ivertex,1) = stress_local(1,1)
    stress(ivertex,2) = stress_local(5,1)
    stress(ivertex,3) = stress_local(9,1)
    stress(ivertex,4) = stress_local(2,1)
    stress(ivertex,5) = stress_local(6,1)
    stress(ivertex,6) = stress_local(3,1)
  enddo

  deallocate(B)
  deallocate(vec_local_disp)
  deallocate(Trans)

end subroutine GeomechForceLocalElemStressStrain

! ************************************************************************** !

subroutine GeomechUpdateSolution(geomech_realization)
  !
  ! Updates data in module after a successful time
  ! step
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/17/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Field_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(geomech_field_type), pointer :: field

  field => geomech_realization%geomech_field

  call GeomechForceUpdateAuxVars(geomech_realization)
  call GeomechUpdateSolutionPatch(geomech_realization)

end subroutine GeomechUpdateSolution

! ************************************************************************** !

subroutine GeomechUpdateSolutionPatch(geomech_realization)
  !
  ! updates data in module after a successful time
  ! step
  !
  ! Author: satish karra, lanl
  ! Date: 09/17/13
  !

  use Geomechanics_Realization_class

  implicit none

  class(realization_geomech_type) :: geomech_realization

  call GeomechForceStressStrain(geomech_realization)

end subroutine GeomechUpdateSolutionPatch

! ************************************************************************** !

subroutine GeomechStoreInitialPressTemp(geomech_realization)
  !
  ! Stores initial pressure and temperature from
  ! subsurface
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/24/13
  !

  use Geomechanics_Realization_class

  implicit none

  class(realization_geomech_type) :: geomech_realization

  PetscErrorCode :: ierr

  call VecCopy(geomech_realization%geomech_field%press_loc, &
               geomech_realization%geomech_field%press_init_loc, &
               ierr);CHKERRQ(ierr)

  call VecCopy(geomech_realization%geomech_field%temp_loc, &
               geomech_realization%geomech_field%temp_init_loc, &
               ierr);CHKERRQ(ierr)

  call VecCopy(geomech_realization%geomech_field%fluid_density_loc, &
               geomech_realization%geomech_field%fluid_density_init_loc, &
               ierr);CHKERRQ(ierr)

end subroutine GeomechStoreInitialPressTemp

! ************************************************************************** !

subroutine GeomechStoreInitialPorosity(realization,geomech_realization)
  !
  ! Stores initial porosity from
  ! subsurface
  !
  ! Author: Satish Karra, LANL
  ! Date: 10/22/13
  !

  use Geomechanics_Realization_class
  use Realization_Subsurface_class
  use Discretization_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  class(realization_subsurface_type) :: realization
  type(discretization_type) :: discretization

  call DiscretizationDuplicateVector(discretization, &
                                     realization%field%work_loc, &
                                     geomech_realization%geomech_field% &
                                     porosity_init_loc)

end subroutine GeomechStoreInitialPorosity

! ************************************************************************** !

subroutine GeomechStoreInitialDisp(geomech_realization)
  !
  ! Stores initial displacement for calculating
  ! relative displacements
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/30/13
  !

  use Geomechanics_Realization_class

  implicit none

  class(realization_geomech_type) :: geomech_realization

  PetscErrorCode :: ierr

  call VecCopy(geomech_realization%geomech_field%disp_xx_loc, &
               geomech_realization%geomech_field%disp_xx_init_loc, &
               ierr);CHKERRQ(ierr)

end subroutine GeomechStoreInitialDisp

end module Geomechanics_Force_module
