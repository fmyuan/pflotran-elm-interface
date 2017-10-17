module Init_Subsurface_module

#include "petsc/finclude/petscvec.h"
  use petscvec
  !geh: there can be no dependencies on simulation object in this file
  use PFLOTRAN_Constants_module

  implicit none

  private


  public :: InitSubsurfAssignMatIDsToRegns, &
            InitSubsurfAssignMatProperties, &
            SubsurfInitMaterialProperties, &
            SubsurfAssignVolsToMatAuxVars, &
            SubsurfSandboxesSetup, &
            InitSubsurfaceSetupZeroArrays
  
contains

! ************************************************************************** !

subroutine SubsurfInitMaterialProperties(realization)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/14
  ! 
  use Realization_Subsurface_class
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  call SubsurfAllocMatPropDataStructs(realization)
  call InitSubsurfAssignMatIDsToRegns(realization)
  call InitSubsurfAssignMatProperties(realization)
  
end subroutine SubsurfInitMaterialProperties

! ************************************************************************** !

subroutine SubsurfAllocMatPropDataStructs(realization)
  ! 
  ! Allocates data structures associated with storage of material properties
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/14
  ! 
  use Realization_Subsurface_class
  use Material_module
  use Option_module
  use Discretization_module
  use Grid_module
  use Patch_module
  use Material_Aux_class
  use Fracture_module, only : FractureAuxVarInit
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscInt :: ghosted_id
  PetscInt :: istart, iend
  PetscInt :: i
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: cur_patch
  class(material_auxvar_type), pointer :: material_auxvars(:)

  option => realization%option

  ! initialize material auxiliary indices.  this does not have to be done
  ! for each patch, just once.
  call MaterialInitAuxIndices(realization%patch%material_property_array,option)

  ! create mappinging
  ! loop over all patches and allocation material id arrays
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    grid => cur_patch%grid
    if (.not.associated(cur_patch%imat)) then
      allocate(cur_patch%imat(grid%ngmax))
      ! initialize to "unset"
      cur_patch%imat = UNINITIALIZED_INTEGER
      select case(option%iflowmode)
        case(NULL_MODE)
        case(RICHARDS_MODE,WF_MODE)
          allocate(cur_patch%cc_id(grid%ngmax))
          cur_patch%cc_id = UNINITIALIZED_INTEGER
        case default
          allocate(cur_patch%cc_id(grid%ngmax))
          cur_patch%cc_id = UNINITIALIZED_INTEGER
          allocate(cur_patch%cct_id(grid%ngmax)) 
          cur_patch%cct_id = UNINITIALIZED_INTEGER
      end select
    endif
    
    cur_patch%aux%Material => MaterialAuxCreate()
    allocate(material_auxvars(grid%ngmax))
    do ghosted_id = 1, grid%ngmax
      call MaterialAuxVarInit(material_auxvars(ghosted_id),option)
      if (option%flow%fracture_on) then
        call FractureAuxVarInit(material_auxvars(ghosted_id))
      endif
    enddo
    cur_patch%aux%Material%num_aux = grid%ngmax
    cur_patch%aux%Material%auxvars => material_auxvars
    nullify(material_auxvars)
    
    cur_patch => cur_patch%next
  enddo

  ! Create Vec that holds compressibility
  if (soil_compressibility_index > 0) then
    call DiscretizationDuplicateVector(realization%discretization, &
                                       realization%field%work, &
                                       realization%field%compressibility0)
  endif

end subroutine SubsurfAllocMatPropDataStructs

! ************************************************************************** !

subroutine InitSubsurfAssignMatIDsToRegns(realization)
  ! 
  ! Assigns material properties to associated regions in the model
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 
  use Realization_Subsurface_class
  use Strata_module
  use Region_module
  use Option_module
  use Grid_module
  use Patch_module
  use Field_module
  use Material_module
  use Material_Aux_class
  use String_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscInt :: icell, local_id, ghosted_id
  PetscInt :: istart, iend
  PetscInt :: local_min, global_min
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(strata_type), pointer :: strata
  type(patch_type), pointer :: cur_patch

  type(material_property_type), pointer :: material_property
  type(region_type), pointer :: region
  class(material_auxvar_type), pointer :: material_auxvars(:)
#ifdef CLM_PFLOTRAN
  type(material_property_type), pointer :: material_property_default
  type(strata_type), pointer :: pre_strata
  character(len=MAXSTRINGLENGTH) :: string
#endif
  
  option => realization%option

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    ! set material ids to uninitialized
    material_auxvars => cur_patch%aux%Material%auxvars
    cur_patch%imat = UNINITIALIZED_INTEGER
    grid => cur_patch%grid

#ifdef CLM_PFLOTRAN
    ! strata(_list) and material_property are unique for each cell, if CLM coupling with PFLOTRAN
    if(cur_patch%strata_list%num_strata /= grid%nlmax) then
      write(string,*) 'num_strata -', cur_patch%strata_list%num_strata, &
        '; grid number -', grid%nlmax
      option%io_buffer = 'CLM-PFLOTRAN: strata_list number is ' // &
              ' not same as grid numbers: ' // &
              trim(string)
      call printErrMsg(option)
    end if
    strata => cur_patch%strata_list%first
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)

      ! pointer to material
      if (cur_patch%surf_or_subsurf_flag == SUBSURFACE) then
        strata%material_property => &
            MaterialPropGetPtrFromArray(strata%material_property_name, &
                                        cur_patch%material_property_array)
        if (.not.associated(strata%material_property)) then
          option%io_buffer = 'Material "' // &
                              trim(strata%material_property_name) // &
                              '" not found in material list'
          call printErrMsg(option)
        endif
      endif

      cur_patch%imat(ghosted_id) = strata%material_property%internal_id

      strata => strata%next

    end do
#else

    strata => cur_patch%strata_list%first
    do
      if (.not.associated(strata)) exit
      ! if not within time period specified, skip the strata.
      ! use a one second tolerance on the start time and end time
      if (StrataWithinTimePeriod(strata,option%time)) then
        ! Read in cell by cell material ids if they exist
        if (.not.associated(strata%region) .and. strata%active) then
          call SubsurfReadMaterialIDsFromFile(realization, &
                                              strata%realization_dependent, &
                                              strata%material_property_filename)
        ! Otherwise, set based on region
        else
          region => strata%region
          material_property => strata%material_property
          if (associated(region)) then
            istart = 1
            iend = region%num_cells
          else
            istart = 1
            iend = grid%nlmax
          endif
          do icell=istart, iend
            if (associated(region)) then
              local_id = region%cell_ids(icell)
            else
              local_id = icell
            endif
            ghosted_id = grid%nL2G(local_id)
            if (strata%active) then
              cur_patch%imat(ghosted_id) = material_property%internal_id
            else
              ! if not active, set material id to zero
              cur_patch%imat(ghosted_id) = 0
            endif
          enddo
        endif
      endif
      strata => strata%next
    enddo
#endif

    cur_patch => cur_patch%next
  enddo
  
  ! ensure that ghosted values for material ids are up to date
  call RealLocalToLocalWithArray(realization,MATERIAL_ID_ARRAY)
  
  ! set material ids in material auxvar.  this must come after the update of 
  ! ghost values.
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    ! set material ids to uninitialized
    material_auxvars => cur_patch%aux%Material%auxvars
    grid => cur_patch%grid
    do ghosted_id = 1, grid%ngmax
      material_auxvars(ghosted_id)%id = cur_patch%imat(ghosted_id)
    enddo
    cur_patch => cur_patch%next
  enddo

end subroutine InitSubsurfAssignMatIDsToRegns

! ************************************************************************** !

subroutine InitSubsurfAssignMatProperties(realization)
  ! 
  ! Assigns material properties based on material ids
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/14
  ! 
  use Realization_Subsurface_class
  use Grid_module
  use Discretization_module
  use Field_module
  use Patch_module
  use Material_Aux_class
  use Material_module
  use Option_module
  use WIPP_module
  use Creep_Closure_module
  use Fracture_module
  use Geomechanics_Subsurface_Properties_module
  use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, PERMEABILITY_XY, &
                               PERMEABILITY_YZ, PERMEABILITY_XZ, &
                               TORTUOSITY, POROSITY, SOIL_COMPRESSIBILITY
  use HDF5_module
  use Grid_Grdecl_module, only : GetPoroPermValues, &
                                 WriteStaticDataAndCleanup, &
                                 DeallocatePoroPermArrays, &
                                 PermPoroExchangeAndSet,SatnumExchangeAndSet, &
                                 GetIsGrdecl,GetSatnumSet,GetSatnumValue
  use Utility_module, only : DeallocateArray
  
  implicit none

  class(realization_subsurface_type) :: realization
  
  PetscReal, pointer :: por0_p(:)
  PetscReal, pointer :: tor0_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: perm_xy_p(:)
  PetscReal, pointer :: perm_xz_p(:)
  PetscReal, pointer :: perm_yz_p(:)
  PetscReal, pointer :: perm_pow_p(:)
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: compress_p(:)
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(material_property_type), pointer :: material_property
  type(material_property_type), pointer :: null_material_property
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(material_type), pointer :: Material
  PetscInt :: local_id, ghosted_id, material_id,natural_id,i
  PetscReal :: tempreal, poro, permx, permy, permz
  PetscInt :: tempint, isatnum, maxsatn
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  PetscInt ,pointer,dimension(:)::inatsend
  PetscBool :: write_ecl, satnum_set

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  
  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_material_property => MaterialPropertyCreate()
  if (option%nflowdof > 0) then
    call VecGetArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
    if (option%flow%full_perm_tensor) then
      call VecGetArrayF90(field%perm0_xy,perm_xy_p,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(field%perm0_xz,perm_xz_p,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(field%perm0_yz,perm_yz_p,ierr);CHKERRQ(ierr)
    endif
    if (soil_compressibility_index > 0) then
      call VecGetArrayF90(field%compressibility0,compress_p,ierr);CHKERRQ(ierr)
    endif
  endif
  call VecGetArrayF90(field%porosity0,por0_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%tortuosity0,tor0_p,ierr);CHKERRQ(ierr)
        
  ! have to use Material%auxvars() and not material_auxvars() due to memory
  ! errors in gfortran
  Material => patch%aux%Material
  
  !if material is associated with fracture, then allocate memory.
  wipp => WIPPGetPtr()
  call CreepClosureConvertListToArray(wipp%creep_closure_tables, &
                                      wipp%creep_closure_tables_array, &
                                      option)
  
  do ghosted_id = 1, grid%ngmax
    material_id = patch%imat(ghosted_id)
    if (material_id > 0) then
      material_property => &
        patch%material_property_array(material_id)%ptr
      
    call GeomechanicsSubsurfacePropsAuxvarInit( &
          material_property%geomechanics_subsurface_properties, &
          patch%aux%Material%auxvars(ghosted_id))
        
      ! lookup creep closure table id from creep closure table name
      if (option%flow%creep_closure_on) then
        material_property%creep_closure_id = &
          CreepClosureGetID(wipp%creep_closure_tables_array, &
                             material_property%creep_closure_name, &
                             material_property%name,option)
      ! copy creep closure table id from material to material_aux
      patch%aux%Material%auxvars(ghosted_id)%creep_closure_id = &
        material_property%creep_closure_id
      endif
    endif
  enddo
  
  !  Prepare for exchange of cell indices and check if satnum set

  satnum_set = PETSC_FALSE
  if (GetIsGrdecl()) then
    if (option%myrank /= option%io_rank) then
      allocate(inatsend(grid%nlmax))
    endif
    satnum_set = GetSatnumSet(maxsatn)
  endif
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    material_id = patch%imat(ghosted_id)
    if (material_id == 0) then
      material_property => null_material_property
    else if (abs(material_id) <= &
             size(patch%material_property_array)) then
      if (material_id < 0) then
        material_property => null_material_property
      else
        material_property => &
          patch%material_property_array(material_id)%ptr
        if (.not.associated(material_property)) then
          write(string,*) &
            patch%imat_internal_to_external(material_id)
          option%io_buffer = 'No material property for material id ' // &
                              trim(adjustl(string)) &
                              //  ' defined in input file.'
          call PrintErrMsgByRank(option)
        endif
      endif
    else if (Uninitialized(material_id)) then 
      write(string,*) grid%nG2A(ghosted_id)
      option%io_buffer = 'Uninitialized material id in patch at cell ' // &
                          trim(adjustl(string))
      call PrintErrMsgByRank(option)
    else if (material_id > size(patch%material_property_array)) then
      write(option%io_buffer,*) patch%imat_internal_to_external(material_id)
      option%io_buffer = 'Unmatched material id in patch: ' // &
        adjustl(trim(option%io_buffer))
      call PrintErrMsgByRank(option)
    else if (material_id == UNMAPPED_MATERIAL_ID) then
      option%io_buffer = 'A material ID has been assigned to a cell that &
        &does not match the IDs of any MATERIAL_PROPERTIES.'
      call PrintErrMsgByRank(option)
    else
      option%io_buffer = 'Something messed up with material ids. Possibly &
        &material ids not assigned to all grid cells. Contact Glenn!'
      call PrintErrMsgByRank(option)
    endif
    if (option%nflowdof > 0) then
      patch%cc_id(ghosted_id) = &
        material_property%saturation_function_id
      if (associated(patch%cct_id)) then
        patch%cct_id(ghosted_id) = &  
          material_property%thermal_conductivity_function_id
      endif
      perm_xx_p(local_id) = material_property%permeability(1,1)
      perm_yy_p(local_id) = material_property%permeability(2,2)
      perm_zz_p(local_id) = material_property%permeability(3,3)
      if (option%flow%full_perm_tensor) then
        perm_xy_p(local_id) = material_property%permeability(1,2)
        perm_xz_p(local_id) = material_property%permeability(1,3)
        perm_yz_p(local_id) = material_property%permeability(2,3)
      endif
      if (soil_compressibility_index > 0) then
        compress_p(local_id) = material_property%soil_compressibility
      endif
    endif
    if (associated(Material%auxvars)) then
      call MaterialAssignPropertyToAux(Material%auxvars(ghosted_id), &
                                        material_property,option)
    endif
    por0_p(local_id) = material_property%porosity
    tor0_p(local_id) = material_property%tortuosity

    if (GetIsGrdecl()) then

      natural_id = grid%nG2A(ghosted_id)

      if (option%myrank == option%io_rank) then
  !  Simply set up the values on the I/O proc
        call GetPoroPermValues(natural_id,poro,permx,permy,permz)
        por0_p(local_id)    = poro
        perm_xx_p(local_id) = permx
        perm_yy_p(local_id) = permy
        perm_zz_p(local_id) = permz
        if (option%flow%full_perm_tensor) then
          perm_xy_p(local_id) = 0.d0
          perm_xz_p(local_id) = 0.d0
          perm_yz_p(local_id) = 0.d0
        endif
        if( satnum_set ) then
  !  Set satnums on this proc
          isatnum = GetSatnumValue(natural_id)
          if (option%nflowdof > 0) then
             patch%cc_id(ghosted_id) = isatnum
          endif
        endif
                
      else
  !  Add to the request list on other procs
        inatsend(local_id)=natural_id
      endif

    endif
  enddo

  if (GetIsGrdecl()) then
    call PermPoroExchangeAndSet(por0_p,perm_xx_p,perm_yy_p,perm_zz_p, &
                                perm_xy_p,perm_xz_p,perm_yz_p, &
                                inatsend,grid%nlmax,option)
    if( satnum_set ) then
      call SatnumExchangeAndSet(patch%cc_id, &
                                inatsend, grid%nlmax, grid%nL2G, option)
    endif  
    if (option%myrank .ne. option%io_rank) then
      call DeallocateArray(inatsend)
    endif
    write_ecl = realization%output_option%write_ecl
    call WriteStaticDataAndCleanup(write_ecl, &
                                realization%output_option%eclipse_options, &
                                option)
    call DeallocatePoroPermArrays(option)
  endif

  call MaterialPropertyDestroy(null_material_property)

  if (option%nflowdof > 0) then
    call VecRestoreArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
    if (option%flow%full_perm_tensor) then
      call VecRestoreArrayF90(field%perm0_xy,perm_xy_p,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(field%perm0_xz,perm_xz_p,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(field%perm0_yz,perm_yz_p,ierr);CHKERRQ(ierr)
    endif
    if (soil_compressibility_index > 0) then
      call VecRestoreArrayF90(field%compressibility0,compress_p, &
                              ierr);CHKERRQ(ierr)
    endif
  endif
  call VecRestoreArrayF90(field%porosity0,por0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%tortuosity0,tor0_p,ierr);CHKERRQ(ierr)
        
  ! read in any user-defined property fields
  do material_id = 1, size(patch%material_property_array)
    material_property => &
            patch%material_property_array(material_id)%ptr
    if (.not.associated(material_property)) cycle
    if (material_property%active) then
      if (associated(material_property%permeability_dataset)) then
        call SubsurfReadPermsFromFile(realization,material_property)
      endif
      if (associated(material_property%compressibility_dataset)) then
        call SubsurfReadDatasetToVecWithMask(realization, &
               material_property%compressibility_dataset, &
               material_property%internal_id,PETSC_FALSE,field%compressibility0)
      endif
      if (associated(material_property%porosity_dataset)) then
        call SubsurfReadDatasetToVecWithMask(realization, &
               material_property%porosity_dataset, &
               material_property%internal_id,PETSC_FALSE,field%porosity0)
        ! if tortuosity is a function of porosity, we must calculate the
        ! the tortuosity on a cell to cell basis.
        if (field%tortuosity0 /= PETSC_NULL_VEC .and. &
            material_property%tortuosity_function_of_porosity) then
          call VecGetArrayF90(field%porosity0,por0_p,ierr);CHKERRQ(ierr)
          call VecGetArrayF90(field%tortuosity0,tor0_p,ierr);CHKERRQ(ierr)
          do local_id = 1, grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%imat(ghosted_id) == material_property%internal_id) then
              tor0_p(local_id) = por0_p(local_id)** &
                material_property%tortuosity_func_porosity_pwr
            endif
          enddo
          call VecRestoreArrayF90(field%porosity0,por0_p,ierr);CHKERRQ(ierr)
          call VecRestoreArrayF90(field%tortuosity0,tor0_p,ierr);CHKERRQ(ierr)
        endif
      endif
      if (associated(material_property%tortuosity_dataset)) then
        call SubsurfReadDatasetToVecWithMask(realization, &
               material_property%tortuosity_dataset, &
               material_property%internal_id,PETSC_FALSE,field%tortuosity0)
      endif
    endif
  enddo
      
  ! update ghosted values
  if (option%nflowdof > 0) then
    call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                                     field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_X,ZERO_INTEGER)
    call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                                     field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Y,ZERO_INTEGER)
    call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                                     field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Z,ZERO_INTEGER)
    if (option%flow%full_perm_tensor) then
      call DiscretizationGlobalToLocal(discretization,field%perm0_xy, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_XY,ZERO_INTEGER)
      call DiscretizationGlobalToLocal(discretization,field%perm0_xz, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_XZ,ZERO_INTEGER)
      call DiscretizationGlobalToLocal(discretization,field%perm0_yz, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_YZ,ZERO_INTEGER)
    endif

    call RealLocalToLocalWithArray(realization,CC_ID_ARRAY)
    if (associated(patch%cct_id)) then
      call RealLocalToLocalWithArray(realization,CCT_ID_ARRAY)
    endif
    
    if (soil_compressibility_index > 0) then
      call DiscretizationGlobalToLocal(discretization,field%compressibility0, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   SOIL_COMPRESSIBILITY,ZERO_INTEGER)
    endif
  endif
  
  call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                                   field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_INITIAL)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_BASE)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_CURRENT)
  call DiscretizationGlobalToLocal(discretization,field%tortuosity0, &
                                    field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               TORTUOSITY,ZERO_INTEGER)

  ! copy rock properties to neighboring ghost cells
  do i = 1, max_material_index
    call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    do local_id = 1, patch%grid%nlmax
      vec_p(local_id) = &
          Material%auxvars(patch%grid%nL2G(local_id))%soil_properties(i)
    enddo
    call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,field%work, &
                                     field%work_loc,ONEDOF)
    call VecGetArrayF90(field%work_loc,vec_p,ierr);CHKERRQ(ierr)
    do ghosted_id = 1, patch%grid%ngmax
      Material%auxvars(ghosted_id)%soil_properties(i) = &
         vec_p(ghosted_id)
    enddo
    call VecRestoreArrayF90(field%work_loc,vec_p,ierr);CHKERRQ(ierr)
  enddo

  if (option%geomech_on) then
    call VecCopy(field%porosity0,field%porosity_geomech_store,ierr);CHKERRQ(ierr)
#ifdef GEOMECH_DEBUG
    print *, 'InitSubsurfAssignMatProperties'
    call PetscViewerASCIIOpen(realization%option%mycomm, &
                              'porosity_geomech_store_por0.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(field%porosity_geomech_store,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  endif

end subroutine InitSubsurfAssignMatProperties

! ************************************************************************** !

subroutine SubsurfReadMaterialIDsFromFile(realization,realization_dependent, &
                                          filename)
  ! 
  ! Reads in grid cell materials
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/03/08
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
  use Discretization_module
  use Logging_module
  use Input_Aux_module
  use Material_module

  use HDF5_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  PetscBool :: realization_dependent
  character(len=MAXSTRINGLENGTH) :: filename
  
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch   
  type(input_type), pointer :: input
  type(discretization_type), pointer :: discretization
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscBool :: append_realization_id
  PetscInt :: ghosted_id, natural_id, material_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscInt, pointer :: external_to_internal_mapping(:)
  Vec :: global_vec
  Vec :: local_vec
  PetscErrorCode :: ierr

  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  discretization => realization%discretization
  
  if (index(filename,'.h5') > 0) then
    group_name = 'Materials'
    dataset_name = 'Material Ids'
    call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                    option)
    call DiscretizationCreateVector(discretization,ONEDOF,local_vec,LOCAL, &
                                    option)
    call HDF5ReadCellIndexedIntegerArray(realization,global_vec, &
                                         filename,group_name, &
                                         dataset_name,realization_dependent)
    call DiscretizationGlobalToLocal(discretization,global_vec,local_vec,ONEDOF)
    call GridCopyVecToIntegerArray(grid,patch%imat,local_vec,grid%ngmax)
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
    call VecDestroy(local_vec,ierr);CHKERRQ(ierr)
  else
    call PetscLogEventBegin(logging%event_hash_map,ierr);CHKERRQ(ierr)
    call GridCreateNaturalToGhostedHash(grid,option)
    input => InputCreate(IUNIT_TEMP,filename,option)
    do
      call InputReadPflotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,natural_id)
      call InputErrorMsg(input,option,'natural id','STRATA')
      ! natural ids in hash are zero-based
      ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (ghosted_id > 0) then
        call InputReadInt(input,option,material_id)
        call InputErrorMsg(input,option,'material id','STRATA')
        patch%imat(ghosted_id) = material_id
      endif
    enddo
    call InputDestroy(input)
    call GridDestroyHashTable(grid)
    call PetscLogEventEnd(logging%event_hash_map,ierr);CHKERRQ(ierr)
  endif
  
  call MaterialCreateExtToIntMapping(patch%material_property_array, &
                                     external_to_internal_mapping)
  call MaterialApplyMapping(external_to_internal_mapping,patch%imat)
  deallocate(external_to_internal_mapping)
  nullify(external_to_internal_mapping)
  
end subroutine SubsurfReadMaterialIDsFromFile

! ************************************************************************** !

subroutine SubsurfReadPermsFromFile(realization,material_property)
  ! 
  ! Reads in grid cell permeabilities
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/19/09
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
  use Discretization_module
  use Logging_module
  use Input_Aux_module
  use Material_module
  use HDF5_module
  use Dataset_Common_HDF5_class
  
  implicit none

  class(realization_subsurface_type) :: realization
  type(material_property_type) :: material_property

  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(discretization_type), pointer :: discretization
  character(len=MAXWORDLENGTH) :: word
  class(dataset_common_hdf5_type), pointer :: dataset_common_hdf5_ptr
  PetscInt :: local_id
  PetscInt :: idirection, temp_int
  PetscReal :: ratio, scale
  Vec :: global_vec
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: perm_xy_p(:)
  PetscReal, pointer :: perm_yz_p(:)
  PetscReal, pointer :: perm_xz_p(:)
  PetscReal, pointer :: perm_ptr(:)

  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  discretization => realization%discretization

  call VecGetArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
  if (option%flow%full_perm_tensor) then
    call VecGetArrayF90(field%perm0_xy,perm_xy_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%perm0_xz,perm_xz_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%perm0_yz,perm_yz_p,ierr);CHKERRQ(ierr)
  endif
  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  if (material_property%isotropic_permeability .or. &
      (.not.material_property%isotropic_permeability .and. &
       Initialized(material_property%vertical_anisotropy_ratio))) then
    ! Although the mask of material ID is applied below, we must only read
    ! in the permeabilities that apply to this material so that small, 
    ! localized gridded datasets (that only apply to a subset of the domain)
    ! can be used.
    call VecZeroEntries(global_vec,ierr);CHKERRQ(ierr)
    call SubsurfReadDatasetToVecWithMask(realization, &
                                    material_property%permeability_dataset, &
                                    material_property%internal_id, &
                                    PETSC_FALSE,global_vec)
    call VecGetArrayF90(global_vec,vec_p,ierr);CHKERRQ(ierr)
    ratio = 1.d0
    scale = 1.d0
    !TODO(geh): fix so that ratio and scale work for perms outside
    ! of dataset
    if (Initialized(material_property%vertical_anisotropy_ratio)) then
      ratio = material_property%vertical_anisotropy_ratio
    endif
    if (material_property%permeability_scaling_factor > 0.d0) then
      scale = material_property%permeability_scaling_factor
    endif
    do local_id = 1, grid%nlmax
      if (patch%imat(grid%nL2G(local_id)) == &
          material_property%internal_id) then
        perm_xx_p(local_id) = vec_p(local_id)*scale
        perm_yy_p(local_id) = vec_p(local_id)*scale
        perm_zz_p(local_id) = vec_p(local_id)*ratio*scale
        if (option%flow%full_perm_tensor) then 
          perm_xy_p(local_id) = 0.d0
          perm_xz_p(local_id) = 0.d0
          perm_yz_p(local_id) = 0.d0
        endif
      endif
    enddo
    call VecRestoreArrayF90(global_vec,vec_p,ierr);CHKERRQ(ierr)
  else
    if (material_property%full_permeability_tensor) then
      temp_int = YZ_DIRECTION
    else
      temp_int = Z_DIRECTION
    endif
    do idirection = X_DIRECTION, temp_int
      select case(idirection)
        case(X_DIRECTION)
          dataset_common_hdf5_ptr => &
             DatasetCommonHDF5Cast(material_property%permeability_dataset)
          perm_ptr => perm_xx_p
        case(Y_DIRECTION)
          dataset_common_hdf5_ptr => &
             DatasetCommonHDF5Cast(material_property%permeability_dataset_y)
          perm_ptr => perm_yy_p
        case(Z_DIRECTION)
          dataset_common_hdf5_ptr => &
             DatasetCommonHDF5Cast(material_property%permeability_dataset_z)
          perm_ptr => perm_zz_p
        case(XY_DIRECTION)
          dataset_common_hdf5_ptr => &
             DatasetCommonHDF5Cast(material_property%permeability_dataset_xy)
          perm_ptr => perm_xy_p
        case(XZ_DIRECTION)
          dataset_common_hdf5_ptr => &
             DatasetCommonHDF5Cast(material_property%permeability_dataset_xz)
          perm_ptr => perm_xz_p
        case(YZ_DIRECTION)
          dataset_common_hdf5_ptr => &
             DatasetCommonHDF5Cast(material_property%permeability_dataset_yz)
          perm_ptr => perm_yz_p
      end select
      ! Although the mask of material ID is applied below, we must only read
      ! in the permeabilities that apply to this material so that small, 
      ! localized gridded datasets (that only apply to a subset of the domain)
      ! can be used.
      call VecZeroEntries(global_vec,ierr);CHKERRQ(ierr)
      call SubsurfReadDatasetToVecWithMask(realization, &
                                           dataset_common_hdf5_ptr,&
                                           material_property%internal_id, &
                                           PETSC_FALSE,global_vec)
      call VecGetArrayF90(global_vec,vec_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax
        if (patch%imat(grid%nL2G(local_id)) == &
            material_property%internal_id) then
          perm_ptr(local_id) = vec_p(local_id)
        endif
      enddo
      call VecRestoreArrayF90(global_vec,vec_p,ierr);CHKERRQ(ierr)
    enddo
  endif
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  
  call VecRestoreArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
  if (option%flow%full_perm_tensor) then
    call VecRestoreArrayF90(field%perm0_xy,perm_xy_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_xz,perm_xz_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_yz,perm_yz_p,ierr);CHKERRQ(ierr)
  endif
  
end subroutine SubsurfReadPermsFromFile

! ************************************************************************** !

subroutine SubsurfReadDatasetToVecWithMask(realization,dataset, &
                                           material_id,read_all_values,vec)
  ! 
  ! Reads a dataset into a PETSc Vec using the material id as a mask
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/19/2016
  ! 
  use Realization_Subsurface_class
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
  use Logging_module
  use Input_Aux_module
  use Material_module
  use HDF5_module
  use Dataset_Base_class
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  
  implicit none

  class(realization_subsurface_type) :: realization
  class(dataset_base_type) :: dataset
  PetscInt :: material_id
  PetscBool :: read_all_values
  Vec :: vec

  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: local_id, ghosted_id, natural_id
  PetscReal :: tempreal
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: work_p(:)
  
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  call VecGetArrayF90(vec,vec_p,ierr);CHKERRQ(ierr)

  if (index(dataset%filename,'.h5') > 0) then
    group_name = ''
    dataset_name = dataset%name
    select type(dataset)
      class is(dataset_gridded_hdf5_type)
        call DatasetGriddedHDF5Load(dataset,option)
        do local_id = 1, grid%nlmax
          ghosted_id = grid%nL2G(local_id)
          if (read_all_values .or. &
              patch%imat(ghosted_id) == material_id) then
            call DatasetGriddedHDF5InterpolateReal(dataset, &
                   grid%x(ghosted_id),grid%y(ghosted_id),grid%z(ghosted_id), &
                   vec_p(local_id),option)
          endif
        enddo
        ! now we strip the dataset to save storage, saving only the name
        ! and filename incase it must be read later
        filename = dataset%filename
        dataset_name = dataset%name
        call DatasetGriddedHDF5Strip(dataset)
        call DatasetGriddedHDF5Init(dataset)
        dataset%filename = filename
        dataset%name = trim(dataset_name)
      class is(dataset_common_hdf5_type)
        dataset_name = dataset%hdf5_dataset_name
        call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                          dataset%filename, &
                                          group_name,dataset_name, &
                                          dataset%realization_dependent)
        call VecGetArrayF90(field%work,work_p,ierr);CHKERRQ(ierr)
        if (read_all_values) then
          vec_p(:) = work_p(:)
        else
          do local_id = 1, grid%nlmax
            if (patch%imat(grid%nL2G(local_id)) == material_id) then
              vec_p(local_id) = work_p(local_id)
            endif
          enddo
        endif
        call VecRestoreArrayF90(field%work,work_p,ierr);CHKERRQ(ierr)
      class default
        option%io_buffer = 'Dataset "' // trim(dataset%name) // '" is of the &
          &wrong type for SubsurfReadDatasetToVecWithMask()'
        call PrintErrMsg(option)
    end select
  else
    call PetscLogEventBegin(logging%event_hash_map,ierr);CHKERRQ(ierr)
    call GridCreateNaturalToGhostedHash(grid,option)
    input => InputCreate(IUNIT_TEMP,dataset%filename,option)
    do
      call InputReadPflotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,natural_id)
      call InputErrorMsg(input,option,'ASCII natural id', &
                         'SubsurfReadDatasetToVecWithMask')
      ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (ghosted_id > 0) then
        if (read_all_values .or. &
            patch%imat(ghosted_id) == material_id) then
          local_id = grid%nG2L(ghosted_id)
          if (local_id > 0) then
            call InputReadDouble(input,option,tempreal)
            call InputErrorMsg(input,option,'dataset value', &
                               'SubsurfReadDatasetToVecWithMask')
            vec_p(local_id) = tempreal
          endif
        endif
      endif
    enddo
    call InputDestroy(input)
    call GridDestroyHashTable(grid)
    call PetscLogEventEnd(logging%event_hash_map,ierr);CHKERRQ(ierr)
  endif
  
  call VecRestoreArrayF90(vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine SubsurfReadDatasetToVecWithMask

! ************************************************************************** !

subroutine SubsurfAssignVolsToMatAuxVars(realization)
  ! 
  ! Assigns the cell volumes currently stored in field%volume0 to the 
  ! material auxiliary variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/14, 12/04/14
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Material_module
  use Discretization_module
  use Field_module
  use Variables_module, only : VOLUME
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field

  option => realization%option
  field => realization%field

  call DiscretizationGlobalToLocal(realization%discretization,field%volume0, &
                                   field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                               field%work_loc,VOLUME,ZERO_INTEGER)

end subroutine SubsurfAssignVolsToMatAuxVars

! ************************************************************************** !

subroutine SubsurfSandboxesSetup(realization)
  ! 
  ! Initializes sandbox objects.
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/06/14, 12/04/14

  use Realization_Subsurface_class
  use SrcSink_Sandbox_module
  
  class(realization_subsurface_type) :: realization
  
  call SSSandboxSetup(realization%patch%grid,realization%option, &
                      realization%output_option)
  
end subroutine SubsurfSandboxesSetup

! ************************************************************************** !

subroutine InitSubsurfaceSetupZeroArrays(realization)
  ! 
  ! Initializes sandbox objects.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/11/16

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Matrix_Zeroing_module

  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  PetscBool, allocatable :: dof_is_active(:)
  PetscInt :: ndof

  type(matrix_zeroing_type), pointer :: matrix_zeroing
  PetscBool :: inactive_cells_exist

  inactive_cells_exist = PETSC_FALSE
  nullify(matrix_zeroing)
  
  option => realization%option
  patch => realization%patch
  
  if (option%nflowdof > 0) then
    allocate(dof_is_active(option%nflowdof))
    dof_is_active = PETSC_TRUE
#if defined(ISOTHERMAL)
    select case(option%iflowmode)
      case(TH_MODE,TH_TS_MODE)
        ! second equation is energy
        dof_is_active(TWO_INTEGER) = PETSC_FALSE
      case(MPH_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE)
        ! third equation is energy
        dof_is_active(THREE_INTEGER) = PETSC_FALSE
    end select
#endif
    select case(option%iflowmode)
      case(RICHARDS_MODE,RICHARDS_TS_MODE)
        matrix_zeroing => patch%aux%Richards%matrix_zeroing
      case(TH_MODE,TH_TS_MODE)
        matrix_zeroing => patch%aux%TH%matrix_zeroing
      case(MPH_MODE)
        matrix_zeroing => patch%aux%Mphase%matrix_zeroing
      case(G_MODE)
        matrix_zeroing => patch%aux%General%matrix_zeroing
      case(H_MODE)
        matrix_zeroing => patch%aux%Hydrate%matrix_zeroing
      case(WF_MODE)
        matrix_zeroing => patch%aux%WIPPFlo%matrix_zeroing
      case(TOIL_IMS_MODE)
        matrix_zeroing => patch%aux%TOil_ims%matrix_zeroing
      case(TOWG_MODE)
        matrix_zeroing => patch%aux%TOWG%matrix_zeroing
      case(IMS_MODE)
        matrix_zeroing => patch%aux%Immis%matrix_zeroing
      case(MIS_MODE)
        matrix_zeroing => patch%aux%Miscible%matrix_zeroing
      case(FLASH2_MODE)
        matrix_zeroing => patch%aux%Flash2%matrix_zeroing
    end select
    call InitSubsurfaceCreateZeroArray(patch,dof_is_active,matrix_zeroing, &
                                       inactive_cells_exist,option)
    select case(option%iflowmode)
      case(RICHARDS_MODE,RICHARDS_TS_MODE)
        patch%aux%Richards%matrix_zeroing => matrix_zeroing
        patch%aux%Richards%inactive_cells_exist = inactive_cells_exist
      case(TH_MODE,TH_TS_MODE)
        patch%aux%TH%matrix_zeroing => matrix_zeroing
        patch%aux%TH%inactive_cells_exist = inactive_cells_exist
      case(MPH_MODE)
        patch%aux%Mphase%matrix_zeroing => matrix_zeroing
        patch%aux%Mphase%inactive_cells_exist = inactive_cells_exist
      case(G_MODE)
        patch%aux%General%matrix_zeroing => matrix_zeroing
        patch%aux%General%inactive_cells_exist = inactive_cells_exist
      case(H_MODE)
        patch%aux%Hydrate%matrix_zeroing => matrix_zeroing
        patch%aux%Hydrate%inactive_cells_exist = inactive_cells_exist
      case(WF_MODE)
        patch%aux%WIPPFlo%matrix_zeroing => matrix_zeroing
        patch%aux%WIPPFlo%inactive_cells_exist = inactive_cells_exist
      case(TOIL_IMS_MODE)
        patch%aux%TOil_ims%matrix_zeroing => matrix_zeroing
        patch%aux%TOil_ims%inactive_cells_exist = inactive_cells_exist
      case(TOWG_MODE)
        !PO: same for all pm_XXX_aux - can be defined in PM_Base_Aux_module
        patch%aux%TOWG%matrix_zeroing => matrix_zeroing
        patch%aux%TOWG%inactive_cells_exist = inactive_cells_exist
      case(IMS_MODE)
        patch%aux%Immis%matrix_zeroing => matrix_zeroing
        patch%aux%Immis%inactive_cells_exist = inactive_cells_exist
      case(MIS_MODE)
        patch%aux%Miscible%matrix_zeroing => matrix_zeroing
        patch%aux%Miscible%inactive_cells_exist = inactive_cells_exist
      case(FLASH2_MODE)
        patch%aux%Flash2%matrix_zeroing => matrix_zeroing
        patch%aux%Flash2%inactive_cells_exist = inactive_cells_exist
    end select
    deallocate(dof_is_active)
  endif

  nullify(matrix_zeroing)
  inactive_cells_exist = PETSC_FALSE

  if (option%ntrandof > 0) then
    select case(option%itranmode)
      case(RT_MODE,EXPLICIT_ADVECTION)
        ! remove ndof above if this is moved
        if (option%transport%reactive_transport_coupling == &
            GLOBAL_IMPLICIT) then
          ndof = realization%reaction%ncomp
        else
          ndof = 1
        endif
        allocate(dof_is_active(ndof))
        dof_is_active = PETSC_TRUE  
      case(NWT_MODE)
        if (option%transport%nw_transport_coupling == GLOBAL_IMPLICIT) then
          ndof = realization%reaction_nw%params%nspecies
        else
          ndof = 1
        endif
        allocate(dof_is_active(ndof))
        dof_is_active = PETSC_TRUE 
    end select
    select case(option%itranmode)
      case(RT_MODE,EXPLICIT_ADVECTION)
        matrix_zeroing => patch%aux%RT%matrix_zeroing
      case(NWT_MODE)
        matrix_zeroing => patch%aux%NWT%matrix_zeroing
    end select
    call InitSubsurfaceCreateZeroArray(patch,dof_is_active,matrix_zeroing, &
                                       inactive_cells_exist,option)
    deallocate(dof_is_active)
    select case(option%itranmode)
      case(RT_MODE,EXPLICIT_ADVECTION)
        patch%aux%RT%matrix_zeroing => matrix_zeroing
        patch%aux%RT%inactive_cells_exist = inactive_cells_exist
      case(NWT_MODE)
        patch%aux%NWT%matrix_zeroing => matrix_zeroing
        patch%aux%NWT%inactive_cells_exist = inactive_cells_exist
    end select
  endif  

end subroutine InitSubsurfaceSetupZeroArrays

! ************************************************************************** !

subroutine InitSubsurfaceCreateZeroArray(patch,dof_is_active, &
                                         matrix_zeroing, &
                                         inactive_cells_exist, &
                                         option)
  ! 
  ! Computes the zeroed rows for inactive grid cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07, 03/02/16, 12/04/19
  ! 
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use String_module
  use Matrix_Zeroing_module
  use Utility_module, only : DeallocateArray
  
  implicit none

  type(patch_type) :: patch
  PetscBool :: dof_is_active(:)
  type(matrix_zeroing_type), pointer :: matrix_zeroing
  PetscBool :: inactive_cells_exist
  type(option_type) :: option
  
  PetscInt :: ncount, idof
  PetscInt :: local_id, ghosted_id
  PetscInt :: ndof, n_active_dof
  PetscInt :: n_inactive_rows

  type(grid_type), pointer :: grid
  PetscInt :: flag
  PetscErrorCode :: ierr
    
  flag = 0
  grid => patch%grid
  ndof = size(dof_is_active)
  n_active_dof = 0
  do idof = 1, ndof
    if (dof_is_active(idof)) n_active_dof = n_active_dof + 1
  enddo
  
  n_inactive_rows = 0
  inactive_cells_exist = PETSC_FALSE

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      n_inactive_rows = n_inactive_rows + ndof
    else if (n_active_dof < ndof) then
      n_inactive_rows = n_inactive_rows + (ndof-n_active_dof)
    endif
  enddo

  call MatrixZeroingInitInactive(matrix_zeroing,n_inactive_rows)

  ncount = 0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      do idof = 1, ndof
        ncount = ncount + 1
        ! 1-based indexing
        matrix_zeroing%inactive_rows_local(ncount) = (local_id-1)*ndof+idof
        ! 0-based indexing
        matrix_zeroing%inactive_rows_local_ghosted(ncount) = &
          (ghosted_id-1)*ndof+idof-1
      enddo
    else if (n_active_dof < ndof) then
      do idof = 1, ndof
        if (dof_is_active(idof)) cycle
        ncount = ncount + 1
        ! 1-based indexing
        matrix_zeroing%inactive_rows_local(ncount) = (local_id-1)*ndof+idof
        ! 0-based indexing
        matrix_zeroing%inactive_rows_local_ghosted(ncount) = &
          (ghosted_id-1)*ndof+idof-1
      enddo
    endif
  enddo

  call MPI_Allreduce(n_inactive_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr)
  if (flag > 0) then
    inactive_cells_exist = PETSC_TRUE
  endif
     
  if (ncount /= n_inactive_rows) then
    option%io_buffer = 'Error:  Mismatch in non-zero row count! ' // &
      StringWrite(ncount) // ' ' // StringWrite(n_inactive_rows)
    call PrintErrMsgByRank(option)
  endif

end subroutine InitSubsurfaceCreateZeroArray

end module Init_Subsurface_module
