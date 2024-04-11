module MpFlow_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use MpFlow_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none
  
  private 

  ! Cutoff parameters
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: unit_z(3) = [0.d0,0.d0,1.d0]

  public :: MpFlowTimeCut,  &
            MpFlowSetup,    &
            MpFlowResidual, &
            MpFlowJacobian, &
            MpFlowMaxChange,          &
            MpFlowUpdateSolution,     &
            MpFlowInitializeTimestep, &
            MpFlowComputeMassBalance, &
            MpFlowUpdateAuxVars,      &
            MpFlowDestroy,            &
            MpFlowAccumulation
         
contains

! ************************************************************************** !

subroutine MpFlowTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  option => realization%option
  field => realization%field
  call MpFlowInitializeTimestep(realization)
 
end subroutine MpFlowTimeCut

! ************************************************************************** !

subroutine MpFlowSetup(realization)
  ! 
  ! Author: ???
  ! Date: 02/22/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Output_Aux_module

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  type(output_variable_list_type), pointer :: list
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MpFlowSetupPatch(realization)
    cur_patch => cur_patch%next
  enddo

  list => realization%output_option%output_snap_variable_list
  call MpFlowSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call MpFlowSetPlotVariables(realization,list)

end subroutine MpFlowSetup

! ************************************************************************** !

subroutine MpFlowSetupPatch(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: ???
  ! Date: 02/22/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Connection_module
  !use Fluid_module
  use Characteristic_Curves_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  !type(fluid_property_type), pointer :: cur_fluid_property
  character(len=MAXWORDLENGTH) :: word


  PetscInt :: ghosted_id, sum_connection
  PetscInt :: i, material_id, sat_func_id
  PetscBool :: error_found
  PetscInt :: iflag
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
    
  patch%aux%MpFlow => MpFlowAuxCreate(option)

  ! mpflow_parameters of materials
  allocate(patch%aux%MpFlow%mpflow_parameters(size(patch%material_property_array)))
  error_found = PETSC_FALSE

  do i = 1, size(patch%material_property_array)

    word = patch%material_property_array(i)%ptr%name 
    if (Uninitialized(patch%material_property_array(i)%ptr%specific_heat)) then
      option%io_buffer = 'ERROR: Non-initialized HEAT_CAPACITY in material ' &
                         // trim(word)
      call printMsg(option)
      error_found = PETSC_TRUE
    endif
    if (Uninitialized(patch%material_property_array(i)%ptr% &
                      thermal_conductivity_wet)) then
      option%io_buffer = 'ERROR: Non-initialized THERMAL_CONDUCTIVITY_WET in ' // &
                         'material ' // trim(word)
      call printMsg(option)
      error_found = PETSC_TRUE
    endif
    if (Uninitialized(patch%material_property_array(i)%ptr% &
                      thermal_conductivity_dry)) then
      option%io_buffer = 'ERROR: Non-initialized THERMAL_CONDUCTIVITY_DRY in' // &
                         'material ' // trim(word)
      call printMsg(option)
      error_found = PETSC_TRUE
    endif

    if (Uninitialized(patch%material_property_array(i)%ptr% &
                        thermal_conductivity_frozen)) then
        option%io_buffer = 'ERROR: Non-initialized THERMAL_CONDUCTIVITY_FROZEN' // &
                           ' in material ' // trim(word)
        call printMsg(option)
        error_found = PETSC_TRUE
    endif
    if (Uninitialized(patch%material_property_array(i)%ptr% &
                        alpha_fr)) then
        option%io_buffer = 'ERROR: Non-initialized THERMAL_COND_EXPONENT_FROZEN' // &
                           ' in material ' // trim(word)
        call printMsg(option)
        error_found = PETSC_TRUE
    endif
    !

    material_id = abs(patch%material_property_array(i)%ptr%internal_id)

    patch%aux%MpFlow%mpflow_parameters(material_id)%dencpr = &
      patch%material_property_array(i)%ptr%matrix_density*1.d-6* &
        patch%material_property_array(i)%ptr%specific_heat

    patch%aux%MpFlow%mpflow_parameters(material_id)%ckwet = &
      patch%material_property_array(i)%ptr%thermal_conductivity_wet*1.d-6

    patch%aux%MpFlow%mpflow_parameters(material_id)%ckdry = &
      patch%material_property_array(i)%ptr%thermal_conductivity_dry*1.d-6

    patch%aux%MpFlow%mpflow_parameters(material_id)%alpha = &
      patch%material_property_array(i)%ptr%alpha

    patch%aux%MpFlow%mpflow_parameters(material_id)%ckfrozen = &
      patch%material_property_array(i)%ptr%thermal_conductivity_frozen*1.d-6

    patch%aux%MpFlow%mpflow_parameters(material_id)%alpha_fr = &
      patch%material_property_array(i)%ptr%alpha_fr
    !

    allocate(patch%aux%MpFlow%mpflow_parameters(material_id)%sr(option%nfluid))
    sat_func_id = patch%material_property_array(i)%ptr%saturation_function_id
    patch%aux%MpFlow%mpflow_parameters(material_id)%sr(:) = &
      CharCurvesGetGetResidualSats(patch%characteristic_curves_array(sat_func_id)%ptr,option)

    allocate(patch%aux%MpFlow%mpflow_parameters(material_id)%diffusion_coefficient(option%nfluid))
    allocate(patch%aux%MpFlow%mpflow_parameters(material_id)%diffusion_activation_energy(option%nfluid))

  enddo

  ! allocate auxvar data structures for all grid cells
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  option%iflag = 1 ! allocate mass_balance array
#else
  option%iflag = 0 ! be sure not to allocate mass_balance array
#endif
  allocate(patch%aux%MpFlow%auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call MpFlowAuxVarInit(patch%aux%MpFlow%auxvars(ghosted_id),option)
  enddo
  patch%aux%MpFlow%num_aux = grid%ngmax

  ! count the number of boundary connections
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo
  patch%aux%MpFlow%num_aux_bc = sum_connection
  if (sum_connection>0) then
    iflag = option%iflag    ! save the originally-set flag
    option%iflag = 1        ! for BC, always has mass_balance data arrays
    allocate(patch%aux%MpFlow%auxvars_bc(sum_connection))
    do ghosted_id = 1, sum_connection
        call MpFlowAuxVarInit(patch%aux%MpFlow%auxvars_bc(ghosted_id),option)
    enddo
    option%iflag = iflag    ! recover the original flag
  endif

  ! Create aux vars for source/sink
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  patch%aux%MpFlow%num_aux_ss = sum_connection
  if (sum_connection>0) then
    iflag = option%iflag    ! save the originally-set flag
    option%iflag = 1        ! for SrcSink, always has mass_balance data arrays
    allocate(patch%aux%MpFlow%auxvars_ss(sum_connection))
    do ghosted_id = 1, sum_connection
        call MpFlowAuxVarInit(patch%aux%MpFlow%auxvars_ss(ghosted_id),option)
    enddo
    option%iflag = iflag    ! recover the original flag
  endif
end subroutine MpFlowSetupPatch


! ************************************************************************** !
! - Residual --------------------------------------------------------
! ************************************************************************** !

subroutine MpFlowResidual(snes,xx,r,realization,ierr)
  !
  ! Computes the residual equation
  !
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
  !

  use Realization_Subsurface_class
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Variables_module
  use Material_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  type(option_type), pointer :: option

  field => realization%field
  discretization => realization%discretization
  option => realization%option

 ! check initial guess -----------------------------------------------
  ierr = MpFlowInitGuessCheck(xx,option)
  if (ierr<0) then
    call SNESSetFunctionDomainError(snes,ierr);CHKERRQ(ierr)
    return
  endif

  ! Communication -----------------------------------------
  ! These 3 must be called before MpFlowUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)

  call DiscretizationLocalToLocal(discretization,field%icap_loc,field%icap_loc,ONEDOF)

  call DiscretizationLocalToLocal(discretization,field%ithrm_loc,field%ithrm_loc,ONEDOF)

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

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    ! flag of PETSC_FALSE not to update patch global_auxvars
    ! because this is an iteration
    cur_patch%aux%MpFlow%auxvars_up_to_date = PETSC_FALSE

    realization%patch => cur_patch

    call MpFlowResidualPatch(r,realization,ierr)
    cur_patch => cur_patch%next
  enddo

end subroutine MpFlowResidual

! ************************************************************************** !

subroutine MpFlowResidualPatch(r,realization,ierr)
  !
  ! Computes the residual equation at patch level
  !
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
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

  Vec, intent(inout) :: r
  type(realization_subsurface_type) :: realization

  ! local variables
  PetscErrorCode :: ierr
  PetscInt :: i
  PetscInt :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: accum_p(:)

  PetscReal, pointer :: r_p(:), xx_loc_p(:)
  PetscReal, pointer :: icap_loc_p(:), ithrm_loc_p(:)

  PetscInt :: ithrm_up, ithrm_dn
  PetscReal :: qsrc1, esrc1

  PetscReal :: Res(realization%option%nflowdof)
  PetscReal :: Res_src(realization%option%nflowdof)
  PetscViewer :: viewer

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field

  type(mpflow_parameter_type), pointer :: mpflow_parameters(:)
  type(mpflow_auxvar_type), pointer :: auxvars(:)
  type(mpflow_auxvar_type), pointer :: auxvars_bc(:)
  type(mpflow_auxvar_type), pointer :: auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: v_darcy
  PetscInt :: iconn, istart, iend, itherm, ii
  PetscInt :: sum_connection
  PetscReal :: fluxe_bulk, fluxe_cond

  PetscReal :: sum_mass_flux(realization%patch%grid%ngmax)

  PetscInt :: ifluid

  !---------------------------------------------------------------------------------------------
  ifluid = IFLUID1   ! (TODO: adding IFLUID2 for air flow)

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  mpflow_parameters => patch%aux%MpFlow%mpflow_parameters
  auxvars => patch%aux%MpFlow%auxvars
  auxvars_bc => patch%aux%MpFlow%auxvars_bc
  auxvars_ss => patch%aux%MpFlow%auxvars_ss

  material_auxvars => patch%aux%Material%auxvars

  call MpFlowUpdateAuxVarsPatch(realization)

  if (option%compute_mass_balance_new) then
    call MpFlowZeroMassBalanceDeltaPatch(realization)

    sum_mass_flux=0.d0
  endif


! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90( r, r_p, ierr);CHKERRQ(ierr)

  r_p = - accum_p

  ! Accumulation terms ------------------------------------

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1

    itherm = int(ithrm_loc_p(ghosted_id))

    call MpFlowAccumulation(auxvars(ghosted_id),      &
                        material_auxvars(ghosted_id), &
                        mpflow_parameters(itherm),    &
                        option,                       &
                        Res)

    r_p(istart:iend) = r_p(istart:iend) + Res

    do ii=1,option%nflowdof
      if(Res(ii) /= Res(ii) .or. abs(Res(ii))>huge(Res(ii)) ) then
        write(string,*) 'local_id: ', local_id, 'Res: ', ii, Res(ii)
        option%io_buffer = ' NaN or INF of Residuals @ MpFlow.F90: FlowResidualPatch - Accumulation of ' // &
          trim(string)
        call printErrMsg(option)
      endif
    enddo

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
      iend = local_id * option%nflowdof
      istart = iend - option%nflowdof + 1
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (source_sink%flow_condition%rate%itype /= HET_MASS_RATE_SS) &
        qsrc1 = source_sink%flow_condition%rate%dataset%rarray(1)

      Res_src = 0.d0
      select case (source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS)
          qsrc1 = qsrc1 / FMW_FLUIDS(ifluid) ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
        case(SCALED_MASS_RATE_SS)
          qsrc1 = qsrc1 / FMW_FLUIDS(ifluid) * &
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn) ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
        case(VOLUMETRIC_RATE_SS)
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*auxvars(ghosted_id)%den(1) ! den = kmol/m^3
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*auxvars(ghosted_id)%den(1)* & ! den = kmol/m^3
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_MASS_RATE_SS)
          qsrc1 = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)/FMW_FLUIDS(ifluid)

        case default
          write(string,*) source_sink%flow_condition%rate%itype
          option%io_buffer='Flow mode source_sink%flow_condition%rate%itype = ' // &
          trim(adjustl(string)) // ', not implemented.'
      end select

      Res_src(IFDOF1) = qsrc1

      esrc1 = 0.d0
      select case(source_sink%flow_condition%itype(IFDOF2))
        case (ENERGY_RATE_SS)
          esrc1 = source_sink%flow_condition%energy_rate%dataset%rarray(1)
        case (SCALED_ENERGY_RATE_SS)
          esrc1 = source_sink%flow_condition%energy_rate%dataset%rarray(1) * &
                  source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case (HET_ENERGY_RATE_SS)
          esrc1 = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
      end select

      ! convert J/s --> MJ/s
      ! default internal energy units are MJ
      Res_src(IFDOF2) = esrc1*1.d-6

      ! Update residual term associated with T
      Res_src(IFDOF2) = Res_src(IFDOF2) + &
          qsrc1*auxvars(ghosted_id)%H(LIQ_FLUID)

      r_p(istart:iend) = r_p(istart:iend) - Res_src

      if (option%compute_mass_balance_new) then
        auxvars_ss(sum_connection)%mass_balance_delta(LIQ_FLUID, IFLOWSPEC1) = &
          auxvars_ss(sum_connection)%mass_balance_delta(LIQ_FLUID, IFLOWSPEC1) - Res_src(IFDOF1)

        ! cumulative mass_balance_delta
        sum_mass_flux(ghosted_id) = sum_mass_flux(ghosted_id) - Res_src(IFDOF1)
      endif
      if (associated(patch%ss_flow_vol_fluxes)) then
        ! fluid flux [m^3/sec] = qsrc_mol [kmol/sec] / den [kmol/m^3]
        patch%ss_flow_vol_fluxes(1,sum_connection) = qsrc1 / &
                                           auxvars(ghosted_id)%den(LIQ_FLUID)
      endif
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(1,sum_connection) = qsrc1
      endif

      Res = Res_src

      do ii=1,option%nflowdof
        if(Res(ii) /= Res(ii) .or. abs(Res(ii))>huge(Res(ii)) ) then
          write(string, *) ' name -', source_sink%name, ' @local_id -', local_id, 'with Res -', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ MpFlow.F90: FlowResidualPatch - source_sink of ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

    enddo
    source_sink => source_sink%next
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

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      if (option_flow%only_vertical_flow) then
        !geh: place second conditional within first to avoid excessive
        !     dot products when .not. option_flow%only_vertical_flow
        if (dot_product(cur_connection_set%dist(1:3,iconn),unit_z) < &
            1.d-10) cycle
      endif

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))

      call MpFlowInternalFlux(auxvars(ghosted_id_up),  &
                  material_auxvars(ghosted_id_up),     &
                  mpflow_parameters(ithrm_up),         &
                  auxvars(ghosted_id_dn),              &
                  material_auxvars(ghosted_id_dn),     &
                  mpflow_parameters(ithrm_dn),         &
                  cur_connection_set%area(iconn),      &
                  cur_connection_set%dist(:,iconn),    &
                  option,                              &
                  v_darcy, fluxe_bulk, fluxe_cond, Res)

      patch%internal_velocities(1,sum_connection) = v_darcy
      patch%internal_flow_fluxes(:,sum_connection) = Res(:)

#ifdef COMPUTE_INTERNAL_MASS_FLUX
      if (option%compute_mass_balance_new) then
        ! contribution to up cells
        auxvars(ghosted_id_up)%mass_balance_delta(LIQ_FLUID, IFLOWSPEC1) = &
          auxvars(ghosted_id_up)%mass_balance_delta(LIQ_FLUID, IFLOWSPEC1) - Res(1)
        ! contribution to dn cells
        auxvars(ghosted_id_dn)%mass_balance_delta(LIQ_FLUID, IFLOWSPEC1) = &
          auxvars(ghosted_id_dn)%mass_balance_delta(LIQ_FLUID, IFLOWSPEC1) + Res(1)

        ! cumulative mass_balance_delta for scaling
        sum_mass_flux(ghosted_id_up) = sum_mass_flux(ghosted_id_up) - Res(1)
        sum_mass_flux(ghosted_id_dn) = sum_mass_flux(ghosted_id_dn) + Res(1)

      endif
#endif

      if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
      endif

      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
      endif

      do ii=1,option%nflowdof
        if(Res(ii) /= Res(ii) .or. abs(Res(ii))>huge(Res(ii)) ) then
          write(string, *) ' @local_id -', local_id_up, local_id_dn, ' with Res -', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ MpFlow.F90: FlowResidualPatch - interior flux between ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

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

      if (patch%imat(ghosted_id) <= 0) cycle

      ithrm_dn = int(ithrm_loc_p(ghosted_id))

      call MpFlowBCFlux(boundary_condition%flow_condition%itype,     &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      auxvars(ghosted_id),                           &
                      material_auxvars(ghosted_id),                  &
                      mpflow_parameters(ithrm_dn),                   &
                      cur_connection_set%area(iconn),                &
                      cur_connection_set%dist(-1:3,iconn),           &
                      option,                                        &
                      v_darcy,                                       &
                      fluxe_bulk, fluxe_cond,                        &
                      Res)


      patch%boundary_velocities(1,sum_connection) = v_darcy
      patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      patch%boundary_energy_flux(1,sum_connection) = fluxe_bulk
      patch%boundary_energy_flux(2,sum_connection) = fluxe_cond

      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        auxvars_bc(sum_connection)%mass_balance_delta(1,1) = &
          auxvars_bc(sum_connection)%mass_balance_delta(1,1) - Res(1)

        ! cumulative mass_balance_delta for scaling
        sum_mass_flux(ghosted_id) = sum_mass_flux(ghosted_id) - Res(1)

      endif

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:option%nflowdof)

      do ii=1,option%nflowdof
        if(Res(ii) /= Res(ii) .or. abs(Res(ii))>huge(Res(ii)) ) then
          write(string, *) ' name -', boundary_condition%name, ' @local_id -', local_id, 'with Res -', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ MpFlow.F90: FlowResidualPatch - boundary_condition of ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! 'mass balance delta' for BC and/or SrcSink conditions scaled by non-zero Residuals
  ! NOTE: if finally solution convergence NOT reached by near-zero Resdiduals,
  !       BC or SrcSink conditions produced mass balance delta have to be scaled,
  !       otherwise mass-balance issue (because mass-balance-delta not recalculated after convergence)
  if (option%compute_mass_balance_new) then

#ifdef COMPUTE_INTERNAL_MASS_FLUX
    connection_set_list => grid%internal_connection_set_list
    cur_connection_set => connection_set_list%first
    sum_connection = 0
    do
      if (.not.associated(cur_connection_set)) exit
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        ghosted_id_up = cur_connection_set%id_up(iconn)
        ghosted_id_dn = cur_connection_set%id_dn(iconn)
        local_id_up = grid%nG2L(ghosted_id_up)
        local_id_dn = grid%nG2L(ghosted_id_dn)
        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle

        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        if(abs(r_p(istart))>1.d-20 .and. abs(sum_mass_flux(ghosted_id_up))>1.d-20) then
          auxvars(ghosted_id_up)%mass_balance_delta = &
            auxvars(ghosted_id_up)%mass_balance_delta * &
            (1.0d0-r_p(istart)/sum_mass_flux(ghosted_id_up))
        endif

        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        if(abs(r_p(istart))>1.d-20 .and. abs(sum_mass_flux(ghosted_id_dn))>1.d-20) then
          auxvars(ghosted_id_dn)%mass_balance_delta = &
            auxvars(ghosted_id_dn)%mass_balance_delta * &
            (1.0d0-r_p(istart)/sum_mass_flux(ghosted_id_dn))
        endif

      enddo
      cur_connection_set => cur_connection_set%next

    enddo
#endif

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
        iend = local_id*option%nflowdof
        istart = iend-option%nflowdof+1

        if(abs(r_p(istart))>1.d-20 .and. abs(sum_mass_flux(ghosted_id))>1.d-20) then
          auxvars_bc(sum_connection)%mass_balance_delta =   &
            auxvars_bc(sum_connection)%mass_balance_delta* &
            (1.0d0-r_p(istart)/sum_mass_flux(ghosted_id))
        endif
      enddo
      boundary_condition => boundary_condition%next
    enddo

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
        iend = local_id * option%nflowdof
        istart = iend - option%nflowdof + 1

        if(abs(r_p(istart))>1.d-20 .and. abs(sum_mass_flux(ghosted_id))>1.d-20) then
          auxvars_ss(sum_connection)%mass_balance_delta =  &
            auxvars_ss(sum_connection)%mass_balance_delta* &
            (1.d0-r_p(istart)/sum_mass_flux(ghosted_id))
        endif
      enddo
      source_sink => source_sink%next
    enddo

  endif !option%compute_mass_balance_new

  ! scale the residual by the volume
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    r_p (istart:iend)= r_p(istart:iend)/material_auxvars(ghosted_id)%volume
  enddo

  if (patch%aux%MpFlow%inactive_cells_exist) then
    do i=1,patch%aux%MpFlow%n_zero_rows
      r_p(patch%aux%MpFlow%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)


  if (realization%debug%vecview_residual) then
    string = 'Flowresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

end subroutine MpFlowResidualPatch



! ************************************************************************** !
! - Jacobian --------------------------------------------------------
! ************************************************************************** !

subroutine MpFlowJacobian(A,B,realization,ierr)
  !
  ! Computes the Jacobian
  !
  ! Author: F-M. Yuan
  ! Date: 02/15/2019
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Debug_module

  implicit none

  Mat :: A, B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(patch_type), pointer :: cur_patch
  type(option_type),  pointer :: option
  PetscReal :: norm

  character(len=MAXSTRINGLENGTH) :: string

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MpFlowJacobianPatch(J, realization,ierr)
    cur_patch => cur_patch%next
  enddo

  if (realization%debug%matview_Jacobian) then
    string = 'Flowjacobian'
    call DebugCreateViewer(realization%debug,string,realization%option,viewer)
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

end subroutine MpFlowJacobian

! ************************************************************************** !

subroutine MpFlowJacobianPatch(AB,realization,ierr)
  !
  ! Computes the Jacobian
  !
  ! Author: F-M. Yuan
  ! Date: 02/1/2019
  !

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_Subsurface_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Characteristic_Curves_module

  Mat :: AB
  type(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr

  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: icap_loc_p(:), ithrm_loc_p(:)

  PetscInt :: ii, jj
  PetscInt :: istart, iend
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: ithrm_up, ithrm_dn

  PetscReal :: qsrc1
  PetscReal :: f_up

  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof), &
               Jsrc(realization%option%nflowdof,realization%option%nflowdof)

  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(mpflow_parameter_type), pointer :: mpflow_parameters(:)
  type(mpflow_auxvar_type), pointer :: auxvars(:)
  !type(global_auxvar_type), pointer :: global_auxvars(:)
  type(mpflow_auxvar_type), pointer :: auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: Rdummy(realization%option%nflowdof)
  PetscReal :: v_darcy, fluxe_bulk, fluxe_cond
  PetscViewer :: viewer
  PetscInt :: ifluid, idof

  ! for 'NUMERICAL_JACOBIAN'
  PetscInt :: icap
  PetscReal :: pert
  PetscReal :: Rpert(realization%option%nflowdof)
  PetscReal :: Jpert(realization%option%nflowdof,realization%option%nflowdof)
  PetscReal :: xx(realization%option%nflowdof),xx_pert(realization%option%nflowdof)
  type(mpflow_auxvar_type) :: auxvar_pert
  class(characteristic_curves_type), pointer :: characteristic_curves

  !---------------------------------------------------------------------------------
  ifluid = IFLUID1

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  mpflow_parameters => patch%aux%MpFlow%mpflow_parameters
  auxvars => patch%aux%MpFlow%auxvars
  auxvars_bc => patch%aux%MpFlow%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)


  ! Accumulation terms ------------------------------------
  do local_id_up = 1, grid%nlmax  ! For each local node do...
    ghosted_id_up = grid%nL2G(local_id_up)
    ! Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id_up) <= 0) cycle
    iend = local_id_up*option%nflowdof
    istart = iend-option%nflowdof+1
    ithrm_up = int(ithrm_loc_p(ghosted_id_up))

    if (option_flow%numerical_derivatives) then
      icap = int(icap_loc_p(ghosted_id_up))
      characteristic_curves => patch%characteristic_curves_array(icap)%ptr

      ! current 'res'
      call MpFlowAccumulation(auxvars(ghosted_id_up),         &
                            material_auxvars(ghosted_id_up),  &
                            mpflow_parameters(ithrm_up),      &
                            option,                           &
                            Rdummy)

      ! for 'NUMERICAL_JACOBIAN', set 'perturbance'
      call MpFlowAuxVarInit(auxvar_pert,option)
      call MpFlowAuxVarCopy(auxvars(ghosted_id_up), auxvar_pert)

      xx(IFDOF1) = xx_loc_p(istart)
      xx(IFDOF2) = xx_loc_p(iend)
      do idof = 1, option%nflowdof
        xx_pert = xx

        pert = xx(idof)*perturbation_tolerance
        if (idof == IFDOF1 .and. xx(IFDOF1)<option%reference_pressure) pert = -pert
        if (idof == IFDOF2 .and. xx(IFDOF2)<0.d0) pert = -pert
        xx_pert(idof) = xx(idof) + pert

        call MpFlowAuxVarCompute(xx_pert,            &
                      characteristic_curves,         &
                      mpflow_parameters(ithrm_up),   &
                      auxvar_pert,                   &
                      option)

        call MpFlowAccumulation(auxvar_pert,                 &
                            material_auxvars(ghosted_id_up), &
                            mpflow_parameters(ithrm_up),     &
                            option,                          &
                            Rpert)

        Jpert(:, idof) = (Rpert(:)-Rdummy(:))/pert
      enddo
      Jup = Jpert
      call MpFlowAuxVarStrip(auxvar_pert)

    else
      ! analytical jacobian
      call MpFlowAccumDerivative(auxvars(ghosted_id_up),     &
                            material_auxvars(ghosted_id_up), &
                            mpflow_parameters(ithrm_up),     &
                            option,                          &
                            Rdummy, Jup, PETSC_TRUE)

    endif

    do ii=1,option%nflowdof
      do jj=1,option%nflowdof
        if(Jup(ii,jj) /= Jup(ii,jj) &
           .or. abs(Jup(ii,jj))>huge(Jup(ii,jj)) ) then
          write(string, *) ' @local_id -', local_id_up, 'with Jacobin -', ii,jj,Jup(ii,jj)
          option%io_buffer = ' NaN or INF of Jacobians @ MpFlow.F90: FlowJacobinPatch - Accumulation ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo
    enddo


    ! scale by the volume of the cell
    Jup = Jup/material_auxvars(ghosted_id_up)%volume

    call MatSetValuesBlockedLocal(AB,1,ghosted_id_up-1,1,ghosted_id_up-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo


  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(AB,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(AB,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_accum'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(AB,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit

    cur_connection_set => source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections

      sum_connection = sum_connection + 1
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)

      if (patch%imat(ghosted_id_dn) <= 0) cycle

      if (source_sink%flow_condition%rate%itype /= HET_MASS_RATE_SS) &
        qsrc1 = source_sink%flow_condition%rate%dataset%rarray(1)

      select case (source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS)
          qsrc1 = qsrc1 / FMW_FLUIDS(ifluid)                  ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
        case(SCALED_MASS_RATE_SS)
          qsrc1 = qsrc1 / FMW_FLUIDS(ifluid) * &
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)  ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
        case(VOLUMETRIC_RATE_SS)
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*auxvars(ghosted_id_dn)%den(LIQ_FLUID)    ! den = kmol/m^3
        case(SCALED_VOLUMETRIC_RATE_SS)
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*auxvars(ghosted_id_dn)%den(LIQ_FLUID)* & ! den = kmol/m^3
                   source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_MASS_RATE_SS)
          qsrc1 = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)/FMW_FLUIDS(ifluid)
        case default
          write(string,*) source_sink%flow_condition%rate%itype
          option%io_buffer='Flow mode source_sink%flow_condition%rate%itype = ' // &
          trim(adjustl(string)) // ', not implemented.'
      end select

      Jsrc = 0.d0

      if (qsrc1 > 0.d0) then ! injection
        Jsrc(IFDOF2,:) = 0.d0
      else
        ! extraction
        Jsrc(IFDOF2,IFDOF1) = &
          -qsrc1*auxvars(ghosted_id_dn)%D_H(LIQ_FLUID, IFDOF1)
        Jsrc(IFDOF2,IFDOF2) = &
          -qsrc1*auxvars(ghosted_id_dn)%D_H(LIQ_FLUID, IFDOF2)
      endif

      do ii=1,option%nflowdof
        do jj=1,option%nflowdof
          if(Jsrc(ii,jj) /= Jsrc(ii,jj) &
             .or. abs(Jsrc(ii,jj))>huge(Jsrc(ii,jj)) ) then
            write(string, *) ' name -', source_sink%name, ' @local_id -', local_id_dn, 'with Jacobin -', ii,jj, Jsrc
            option%io_buffer = ' NaN or INF of Jacobians @ MpFlow.F90: FlowJacobinPatch - Source_Sink of ' // &
              trim(string)
            call printErrMsg(option)
          endif
        enddo
      enddo

      ! scale by the volume of the cell
      Jsrc = Jsrc/material_auxvars(ghosted_id_dn)%volume

      call MatSetValuesBlockedLocal(AB,1,ghosted_id_dn-1,1,ghosted_id_dn-1,Jsrc, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)


    enddo
    source_sink => source_sink%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(AB,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(AB,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_srcsink'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(AB,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

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

      if (option_flow%only_vertical_flow) then
        !geh: place second conditional within first to avoid excessive
        !     dot products when .not. option_flow%only_vertical_flow
        if (dot_product(cur_connection_set%dist(1:3,iconn),unit_z) < &
            1.d-10) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))

      if (option_flow%numerical_derivatives) then

        ! current 'res'
        call MpFlowInternalFlux(auxvars(ghosted_id_up), &
                   material_auxvars(ghosted_id_up),     &
                   mpflow_parameters(ithrm_up),         &
                   auxvars(ghosted_id_dn),              &
                   material_auxvars(ghosted_id_dn),     &
                   mpflow_parameters(ithrm_dn),         &
                   cur_connection_set%area(iconn),      &
                   cur_connection_set%dist(:,iconn),    &
                   option,                              &
                   v_darcy, fluxe_bulk, fluxe_cond,     &
                   Rdummy)


        ! for 'NUMERICAL_JACOBIAN' of 'up'
        icap = int(icap_loc_p(ghosted_id_up))
        characteristic_curves => patch%characteristic_curves_array(icap)%ptr

        call MpFlowAuxVarInit(auxvar_pert,option)
        call MpFlowAuxVarCopy(auxvars(ghosted_id_up), auxvar_pert)

        iend   = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        xx(IFDOF1) = xx_loc_p(istart)
        xx(IFDOF2) = xx_loc_p(iend)
        do idof = 1, option%nflowdof
          xx_pert = xx

          pert = xx(idof)*perturbation_tolerance
          if (idof == IFDOF1 .and. xx(IFDOF1)<option%reference_pressure) pert = -pert
          if (idof == IFDOF2 .and. xx(IFDOF2)<0.d0) pert = -pert
          xx_pert(idof) = xx(idof) + pert

          call MpFlowAuxVarCompute(xx_pert,          &
                        characteristic_curves,       &
                        mpflow_parameters(ithrm_up), &
                        auxvar_pert,                 &
                        option)
          call MpFlowInternalFlux(                        &
                     auxvar_pert,                         &  ! 'auxvar_up' perturbation
                     material_auxvars(ghosted_id_up),     &
                     mpflow_parameters(ithrm_up),         &
                     auxvars(ghosted_id_dn),              &
                     material_auxvars(ghosted_id_dn),     &
                     mpflow_parameters(ithrm_dn),         &
                     cur_connection_set%area(iconn),      &
                     cur_connection_set%dist(:,iconn),    &
                     option,                              &
                     v_darcy, fluxe_bulk, fluxe_cond,     &
                     Rpert)


          Jpert(:, idof) = (Rpert(:)-Rdummy(:))/pert
        enddo
        Jup = Jpert
        call MpFlowAuxVarStrip(auxvar_pert)

        ! for 'NUMERICAL_JACOBIAN' of 'dn'
        icap = int(icap_loc_p(ghosted_id_dn))
        characteristic_curves => patch%characteristic_curves_array(icap)%ptr

        call MpFlowAuxVarInit(auxvar_pert,option)
        call MpFlowAuxVarCopy(auxvars(ghosted_id_dn), auxvar_pert)

        iend   = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        xx(IFDOF1) = xx_loc_p(istart)
        xx(IFDOF2) = xx_loc_p(iend)
        do idof = 1, option%nflowdof
          xx_pert = xx

          pert = xx(idof)*perturbation_tolerance
          if (idof == IFDOF1 .and. xx(IFDOF1)<option%reference_pressure) pert = -pert
          if (idof == IFDOF2 .and. xx(IFDOF2)<0.d0) pert = -pert
          xx_pert(idof) = xx(idof) + pert

          call MpFlowAuxVarCompute(xx_pert,          &
                        characteristic_curves,       &
                        mpflow_parameters(ithrm_up), &
                        auxvar_pert,                 &
                        option)
          call MpFlowInternalFlux(                        &
                     auxvars(ghosted_id_up),              &
                     material_auxvars(ghosted_id_up),     &
                     mpflow_parameters(ithrm_up),         &
                     auxvar_pert,                         &  ! 'auxvar_dn' perturbation
                     material_auxvars(ghosted_id_dn),     &
                     mpflow_parameters(ithrm_dn),         &
                     cur_connection_set%area(iconn),      &
                     cur_connection_set%dist(:,iconn),    &
                     option,                              &
                     v_darcy, fluxe_bulk, fluxe_cond,     &
                     Rpert)


          Jpert(:, idof) = (Rpert(:)-Rdummy(:))/pert
        enddo
        Jdn = Jpert
        call MpFlowAuxVarStrip(auxvar_pert)
        !

      else
        ! analytical jacobian

        call MpFlowInternalFluxDerivative(auxvars(ghosted_id_up),          &
                                      material_auxvars(ghosted_id_up),     &
                                      mpflow_parameters(ithrm_up),         &
                                      auxvars(ghosted_id_dn),              &
                                      material_auxvars(ghosted_id_dn),     &
                                      mpflow_parameters(ithrm_dn),         &
                                      cur_connection_set%area(iconn),      &
                                      cur_connection_set%dist(-1:3,iconn), &
                                      option,                              &
                                      v_darcy, fluxe_bulk, fluxe_cond,     &
                                      Rdummy, Jup, Jdn, PETSC_TRUE)

      endif  ! end if (option_flow%numerical_derivatives)



      do ii=1,option%nflowdof
        do jj=1,option%nflowdof
          if(Jup(ii,jj) /= Jup(ii,jj) .or. Jdn(ii,jj) /= Jdn(ii,jj) &
             .or. abs(Jup(ii,jj))>huge(Jup(ii,jj)) .or. abs(Jdn(ii,jj))>huge(Jdn(ii,jj)) ) then
            write(string, *) ' between local_id up/dn -', local_id_up, local_id_dn, &
                'with Jacobin -', ii,jj, Jup(ii,jj), Jdn(ii,jj)
            option%io_buffer = ' NaN or INF of Jacobians @ MpFlow.F90: FlowJacobinPatch - Interior flux ' // &
                trim(string)
            call printErrMsg(option)
          endif
        enddo
      enddo


!  scale by the volume of the cell

      if (local_id_up > 0) then
        call MatSetValuesBlockedLocal(AB,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup/material_auxvars(ghosted_id_up)%volume,ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(AB,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn/material_auxvars(ghosted_id_up)%volume,ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn

        call MatSetValuesBlockedLocal(AB,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn/material_auxvars(ghosted_id_dn)%volume,ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(AB,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup/material_auxvars(ghosted_id_dn)%volume,ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(AB,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(AB,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_flux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(AB,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit

    cur_connection_set => boundary_condition%connection_set

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)
      if (patch%imat(ghosted_id_dn) <= 0) cycle

      ithrm_dn  = int(ithrm_loc_p(ghosted_id_dn))

      if (option_flow%numerical_derivatives) then

        ! current 'res'
        call MpFlowBCFlux(boundary_condition%flow_condition%itype,   &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      auxvars(ghosted_id_dn),                        &
                      material_auxvars(ghosted_id_dn),               &
                      mpflow_parameters(ithrm_dn),                   &
                      cur_connection_set%area(iconn),                &
                      cur_connection_set%dist(-1:3,iconn),           &
                      option,                                        &
                      v_darcy,                                       &
                      fluxe_bulk, fluxe_cond,                        &
                      Rdummy)

        ! for 'NUMERICAL_JACOBIAN' of 'dn'
        icap = int(icap_loc_p(ghosted_id_dn))
        characteristic_curves => patch%characteristic_curves_array(icap)%ptr

        call MpFlowAuxVarInit(auxvar_pert,option)
        call MpFlowAuxVarCopy(auxvars(ghosted_id_dn), auxvar_pert)

        iend   = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        xx(IFDOF1) = xx_loc_p(istart)
        xx(IFDOF2) = xx_loc_p(iend)
        do idof = 1, option%nflowdof
          xx_pert = xx

          pert = xx(idof)*perturbation_tolerance
          if (idof == IFDOF1 .and. xx(IFDOF1)<option%reference_pressure) pert = -pert
          if (idof == IFDOF2 .and. xx(IFDOF2)<0.d0) pert = -pert
          xx_pert(idof) = xx(idof) + pert

          call MpFlowAuxVarCompute(xx_pert,          &
                        characteristic_curves,       &
                        mpflow_parameters(ithrm_up), &
                        auxvar_pert,                 &
                        option)
          call MpFlowBCFlux(boundary_condition%flow_condition%itype, &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      auxvar_pert,                                   &  ! 'auxvar_dn' perturbalated
                      material_auxvars(ghosted_id_dn),               &
                      mpflow_parameters(ithrm_dn),                   &
                      cur_connection_set%area(iconn),                &
                      cur_connection_set%dist(-1:3,iconn),           &
                      option,                                        &
                      v_darcy,                                       &
                      fluxe_bulk, fluxe_cond,                        &
                      Rdummy)

          Jpert(:, idof) = (Rpert(:)-Rdummy(:))/pert
        enddo
        Jdn = Jpert
        call MpFlowAuxVarStrip(auxvar_pert)
        !

      else
        ! analytical jacobian
        call MpFlowBCFluxDerivative(boundary_condition%flow_condition%itype,   &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                auxvars(ghosted_id_dn),                        &
                                material_auxvars(ghosted_id_dn),               &
                                mpflow_parameters(ithrm_dn),                   &
                                cur_connection_set%area(iconn),                &
                                cur_connection_set%dist(-1:3,iconn),           &
                                option,                                        &
                                v_darcy, fluxe_bulk, fluxe_cond,               &
                                Rdummy, Jdn, PETSC_TRUE)

      endif ! end of if(option_flow%numerical_derivatives)

      !
      Jdn = -Jdn

      do ii=1,option%nflowdof
        do jj=1,option%nflowdof
          if(Jdn(ii,jj) /= Jdn(ii,jj) .or. &
             abs(Jdn(ii,jj))>huge(Jdn(ii,jj)) ) then
            write(string, *) ' name -', boundary_condition%name, ' @local_id -', local_id_dn, 'with Jacobin -', ii,jj, Jdn(ii,jj)
            option%io_buffer = ' NaN or INF of Jacobians @ MpFlow.F90: FlowJacobinPatch - Boundary_Condition of ' // &
                trim(string)
            call printErrMsg(option)
          endif
        enddo
      enddo

      !  scale by the volume of the cell
      Jdn = Jdn/material_auxvars(ghosted_id_dn)%volume

      call MatSetValuesBlockedLocal(AB,1,ghosted_id_dn-1,1,ghosted_id_dn-1,Jdn,ADD_VALUES, &
                                    ierr);CHKERRQ(ierr)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(AB,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(AB,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_bcflux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(AB,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(AB,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(AB,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  if (patch%aux%MpFlow%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(AB,patch%aux%MpFlow%n_zero_rows, &
                          patch%aux%MpFlow%zero_rows_local_ghosted,f_up, &
                          PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif

end subroutine MpFlowJacobianPatch

! ************************************************************************** !
! - Others --------------------------------------------------------
! ************************************************************************** !

subroutine MpFlowMaxChange(realization,dpmax,dtmpmax)
  !
  ! Computes the maximum change in the solution vector
  !
  ! Author: ???
  ! Date: 01/15/08
  !

  use Realization_Subsurface_class
  use Option_module
  use Field_module

  implicit none

  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(field_type), pointer :: field

  PetscReal :: dpmax, dtmpmax
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field

  dpmax = 0.d0
  dtmpmax = 0.d0

  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy, &
                ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,dpmax, &
                     ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,dtmpmax, &
                     ierr);CHKERRQ(ierr)

end subroutine MpFlowMaxChange

! ************************************************************************** !

subroutine MpFlowComputeMassBalance(realization, mass_balance)
  ! 
  ! MpFlowomputeMassBalance:
  ! 

  use Realization_Subsurface_class
  use Patch_module

  type(realization_subsurface_type) :: realization
  PetscReal, pointer :: mass_balance(:,:)
   
  type(patch_type), pointer :: cur_patch

  ! for all patch(es)
  mass_balance = 0.d0

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MpFlowComputeMassBalancePatch(realization, mass_balance)
    cur_patch => cur_patch%next
  enddo

end subroutine MpFlowComputeMassBalance

! ************************************************************************** !

subroutine MpFlowComputeMassBalancePatch(realization,mass_balance)
  ! 
  ! MpFlowomputeMassBalancePatch:
  !  Compute Sum of Mass_Balance for one Patch
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class, only : material_auxvar_type, &
                                 soil_compressibility_index, &
                                 MaterialCompressSoil
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal,pointer :: mass_balance(:,:)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(mpflow_auxvar_type),pointer :: mpflow_auxvars(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iflowsp
  PetscReal :: compressed_porosity, pres
  PetscReal :: por
  PetscReal :: tempreal

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  !global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  mpflow_auxvars => patch%aux%MpFlow%auxvars

  ! note: don't zero mass_balance here,
  !       which already done prior to first patch of realization
  do iflowsp = 1, option_flow%nspecflow
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (soil_compressibility_index > 0) then
        pres = max(mpflow_auxvars(ghosted_id)%pres(LIQ_FLUID), &
                 mpflow_auxvars(ghosted_id)%pres(AIR_FLUID))
        call MaterialCompressSoil(material_auxvars(ghosted_id), &
                                 pres, &
                                 compressed_porosity,tempreal)
        por = compressed_porosity
      else
        por = material_auxvars(ghosted_id)%porosity
      endif

      ! mass = por*volume*saturation*den_kg*fraction_in_fluid
      tempreal = por*material_auxvars(ghosted_id)%volume

      ! in 'liq_fluid'
      if (option_flow%nspecliq>0) then
        mass_balance(LIQ_FLUID, iflowsp) = mass_balance(LIQ_FLUID, iflowsp) + &
          mpflow_auxvars(ghosted_id)%sat(LIQ_FLUID)* &
          mpflow_auxvars(ghosted_id)%den_kg(LIQ_FLUID)* &
          tempreal
      endif

      ! in 'air_fluid'
      if (option_flow%nspecgas>0) then
        mass_balance(AIR_FLUID, iflowsp) = mass_balance(AIR_FLUID, iflowsp) + &
          mpflow_auxvars(ghosted_id)%sat(AIR_FLUID)* &
          mpflow_auxvars(ghosted_id)%den_kg(AIR_FLUID)* &
          tempreal
      endif

      ! immobilized (e.g. frozen) 'fluid'
      tempreal = tempreal * (1.d0 - mpflow_auxvars(ghosted_id)%sat(LIQ_FLUID) &
                                  - mpflow_auxvars(ghosted_id)%sat(AIR_FLUID) )
      mass_balance(option%nfluid+1,iflowsp) = &
        mass_balance(option%nfluid+1,iflowsp) + &
        mpflow_auxvars(ghosted_id)%den_kg(option%nfluid+1)* &
        tempreal

    enddo
  enddo

end subroutine MpFlowComputeMassBalancePatch

! ************************************************************************** !

subroutine MpFlowZeroMassBalanceDeltaPatch(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 12/13/11
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  type(MpFlow_auxvar_type), pointer :: auxvars(:)
#endif
  type(MpFlow_auxvar_type), pointer :: auxvars_bc(:)
  type(MpFlow_auxvar_type), pointer :: auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  auxvars    => patch%aux%MpFlow%auxvars
#endif
  auxvars_bc => patch%aux%MpFlow%auxvars_bc
  auxvars_ss => patch%aux%MpFlow%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%MpFlow%num_aux
    auxvars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%MpFlow%num_aux_bc > 0) then
    do iconn = 1, patch%aux%MpFlow%num_aux_bc
      auxvars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
  if (patch%aux%MpFlow%num_aux_ss > 0) then
    do iconn = 1, patch%aux%MpFlow%num_aux_ss
      auxvars_ss(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
 
end subroutine MpFlowZeroMassBalanceDeltaPatch

! ************************************************************************** !

subroutine MpFlowUpdateMassBalanceDeltaPatch(realization)
  ! 
  ! Updates mass balance
  !  Compute Mass_Balance Delta, since 'ZeroMassBalanceDelta' called, in a Patch for each gridcell
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  type(MpFlow_auxvar_type), pointer :: auxvars(:)
#endif
  type(MpFlow_auxvar_type), pointer :: auxvars_bc(:)
  type(MpFlow_auxvar_type), pointer :: auxvars_ss(:)

  PetscInt :: iconn, ifluid, iflowsp

  !---------------------------------------------------------------------

  option => realization%option
  patch => realization%patch

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  auxvars    => patch%aux%MpFlow%auxvars
#endif
  auxvars_bc => patch%aux%MpFlow%auxvars_bc
  auxvars_ss => patch%aux%MpFlow%auxvars_ss

  ! note: we don't zero mass_balance here, so it's upon where this subroutine called for 'balance' as accumulative or not.

  iflowsp= IFLOWSPEC1
  do ifluid = 1, option%nfluid+1 ! extra '1' for immoblized mass of a fluid (eq. water->ice, vapor->ice, air->frozen air?)
    ! here now water (indexed as 'IFLOW1') only.
    ! NOTE- if no gas-fluid considered, vapor will not be accounted; if no liq-fluid, vapor to either ice or water will no be accounted
#ifdef COMPUTE_INTERNAL_MASS_FLUX
    do iconn = 1, patch%aux%MpFlow%num_aux
      auxvars(iconn)%mass_balance(ifluid, iflowsp) = &
        auxvars(iconn)%mass_balance(ifluid, iflowsp) + &
        auxvars(iconn)%mass_balance_delta(ifluid, iflowsp)* &
        FMW_FLUIDS(iflowsp)* &
        option_flow%dt
    enddo
#endif

    if (patch%aux%MpFlow%num_aux_bc > 0) then
      do iconn = 1, patch%aux%MpFlow%num_aux_bc
        auxvars_bc(iconn)%mass_balance(ifluid, iflowsp) = &
          auxvars_bc(iconn)%mass_balance(ifluid, iflowsp)+ &
          auxvars_bc(iconn)%mass_balance_delta(ifluid, iflowsp)* &
          FMW_FLUIDS(iflowsp)* &
          option_flow%dt
      enddo
    endif

    if (patch%aux%MpFlow%num_aux_ss > 0) then
      do iconn = 1, patch%aux%MpFlow%num_aux_ss
        auxvars_ss(iconn)%mass_balance(ifluid,iflowsp) = &
          auxvars_ss(iconn)%mass_balance(ifluid, iflowsp)+ &
          auxvars_ss(iconn)%mass_balance_delta(ifluid, iflowsp)* &
          FMW_FLUIDS(IFLOWSPEC1)* &
          option_flow%dt
      enddo
    endif

  enddo

end subroutine MpFlowUpdateMassBalanceDeltaPatch

! ************************************************************************** !

subroutine MpFlowUpdateAuxVars(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the Flow problem
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch

    ! the flag to update global_auxvars
    ! because this is auxvar computing after solution
    cur_patch%aux%MpFlow%auxvars_up_to_date = PETSC_TRUE

    call MpFlowUpdateAuxVarsPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine MpFlowUpdateAuxVars

! ************************************************************************** !

subroutine MpFlowUpdateAuxVarsPatch(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the Flow problem
  ! 
  ! Author: ???
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
  use Characteristic_Curves_module
   
  implicit none

  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(mpflow_auxvar_type), pointer :: mpflow_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(mpflow_parameter_type), pointer :: mpflow_parameters(:)
  class(characteristic_curves_type), pointer :: characteristic_curves

  PetscInt :: ghosted_id, local_id, istart, iend, iconn, sum_connection, idof
  PetscInt :: ithrm, icap
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:), icap_loc_p(:)

  PetscReal :: xx(realization%option%nflowdof)

  PetscErrorCode :: ierr
  !------------------------------------------------------------------------------------
  
  option => realization%option

  field => realization%field
  patch => realization%patch
  grid => patch%grid

  material_auxvars  => patch%aux%Material%auxvars
  mpflow_auxvars    => patch%aux%MpFlow%auxvars
  mpflow_parameters => patch%aux%MpFlow%mpflow_parameters

  global_auxvars      => patch%aux%Global%auxvars
  global_auxvars_bc   => patch%aux%Global%auxvars_bc
  global_auxvars_ss   => patch%aux%Global%auxvars_ss

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    if (patch%imat(ghosted_id) <= 0) cycle
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    ithrm = int(ithrm_loc_p(ghosted_id))

    xx(1) = xx_loc_p(istart)
    xx(2) = xx_loc_p(iend)

    icap = int(icap_loc_p(ghosted_id))
    characteristic_curves => patch%characteristic_curves_array(icap)%ptr

    call MpFlowAuxVarCompute(xx,                 &
            characteristic_curves,               &
            mpflow_parameters(ithrm),            &
            mpflow_auxvars(ghosted_id),          &
            option)

    ! Assign values to GLOBAL_auxvar, which shared with other MODEs in PFLOTRAN
    ! when instructed so (e.g. when iteration is done and thereafter passing 'mpflow_auxvar' to 'global_auxvar')
    if (patch%aux%MpFlow%auxvars_up_to_date) then
      global_auxvars(ghosted_id)%tc   = mpflow_auxvars(ghosted_id)%tc
      global_auxvars(ghosted_id)%pres = mpflow_auxvars(ghosted_id)%pres
      global_auxvars(ghosted_id)%sat  = mpflow_auxvars(ghosted_id)%sat
      global_auxvars(ghosted_id)%den  = mpflow_auxvars(ghosted_id)%den
    endif

  enddo

  ! BCs
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
      ithrm = int(ithrm_loc_p(ghosted_id))

      ! set global_auxvars to that connected cell at first,
      ! then, reset relevant vars to BC's values
      call GlobalAuxVarCopy(global_auxvars(ghosted_id), &
                            global_auxvars_bc(sum_connection))

      do idof=1,option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC, &
               HET_SURF_HYDROSTATIC_SEEPAGE_BC, &
               HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC)
            if (idof==IFDOF1) then
              global_auxvars_bc(sum_connection)%pres(LIQ_FLUID) = &
                boundary_condition%flow_aux_real_var(IFDOF1,iconn)
            elseif(idof==IFDOF2) then
              global_auxvars_bc(sum_connection)%tc = &
                boundary_condition%flow_aux_real_var(IFDOF2,iconn)
            endif

          case(HYDROSTATIC_CONDUCTANCE_BC, HET_HYDROSTATIC_CONDUCTANCE_BC)
            !(TODO) the following is incorrect?
            if(idof==IFDOF2) then
              global_auxvars_bc(sum_connection)%tc = &
                boundary_condition%flow_aux_real_var(IFDOF3,iconn)
            endif

          case(NEUMANN_BC,ZERO_GRADIENT_BC)
             ! nothing to do here

        end select
      enddo
      
    enddo
    boundary_condition => boundary_condition%next
  enddo

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
      ithrm = int(ithrm_loc_p(ghosted_id))

      iend = ghosted_id*option%nflowdof
      istart = iend-option%nflowdof+1

      ! set global_auxvars to that connected cell at first,
      ! then reset its relevant vars to SrcSink's
      call GlobalAuxVarCopy(global_auxvars(ghosted_id), &
                            global_auxvars_bc(sum_connection))

#if 0
      ! (TODO)
      select case(source_sink%flow_condition%itype(IFDOF2))
        case (HET_DIRICHLET_BC)
          tsrc1 = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
        case (DIRICHLET_BC)
          tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
        case (ENERGY_RATE_SS,SCALED_ENERGY_RATE_SS,HET_ENERGY_RATE_SS, &
              ZERO_GRADIENT_BC)
          tsrc1 = xx_loc_p((ghosted_id-1)*option%nflowdof+2)
        case default
          option%io_buffer='Unsupported temperature flow condtion for ' // &
            'a source-sink in Flow mode: ' // trim(source_sink%name)
          call printErrMsg(option)
      end select
#endif

    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)


end subroutine MpFlowUpdateAuxVarsPatch

! ************************************************************************** !

subroutine MpFlowInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: ???
  ! Date: 02/20/08
  ! 

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  call MpFlowUpdateFixedAccumulation(realization)

end subroutine MpFlowInitializeTimestep

! ************************************************************************** !

subroutine MpFlowUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time step
  ! 
  ! Author: ???
  ! Date: 02/13/08
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  
  field => realization%field
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MpFlowUpdateSolutionPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine MpFlowUpdateSolution

! ************************************************************************** !

subroutine MpFlowUpdateSolutionPatch(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 12/13/11, 02/28/14
  ! 


  use Realization_Subsurface_class
  use Patch_module
    
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(patch_type), pointer :: patch
  type(mpflow_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)

  patch => realization%patch
  auxvars => patch%aux%MpFlow%auxvars
  global_auxvars => patch%aux%Global%auxvars
  
  ! (TODO checking here? - ONLY thing to do ?)
  !if (realization%option%compute_mass_balance_new) then
  !  call MpFlowUpdateMassBalancePatch(realization)
  !endif

end subroutine MpFlowUpdateSolutionPatch

! ************************************************************************** !

subroutine MpFlowUpdateFixedAccumulation(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    ! flag of PETSC_FALSE not to update patch global_auxvars
    ! because this is an iteration
    cur_patch%aux%MpFlow%auxvars_up_to_date = PETSC_FALSE

    realization%patch => cur_patch
    call MpFlowUpdateFixedAccumPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine MpFlowUpdateFixedAccumulation

! ************************************************************************** !

subroutine MpFlowUpdateFixedAccumPatch(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Characteristic_Curves_module

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  ! locals
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(mpflow_auxvar_type), pointer :: mpflow_auxvars(:)
  type(mpflow_parameter_type), pointer :: mpflow_parameters(:)
  class(characteristic_curves_type), pointer :: characteristic_curves

  !class(material_auxvar_type), pointer :: material_auxvars(:)  ! this unknownly would likely mess up %soil_properties(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, istart, iend
  PetscReal, pointer :: xx_p(:), icap_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:), accum_p(:)
  PetscInt :: ithrm, icap
  PetscReal :: xx(realization%option%nflowdof)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  mpflow_parameters => patch%aux%MpFlow%mpflow_parameters
  mpflow_auxvars    => patch%aux%MpFlow%auxvars
  global_auxvars  => patch%aux%Global%auxvars
  material_auxvars=> patch%aux%Material%auxvars

  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1

    ! material property id
    ithrm = int(ithrm_loc_p(ghosted_id))

    ! charateristic curve id
    icap = int(icap_loc_p(ghosted_id))
    characteristic_curves => patch%characteristic_curves_array(icap)%ptr

    xx(1) = xx_p(istart)
    xx(2) = xx_p(iend)

    call MpFlowAuxVarCompute(xx,                 &
            characteristic_curves,               &
            mpflow_parameters(ithrm),            &
            mpflow_auxvars(ghosted_id),          &
            option)
    
    call MpFlowAccumulation(mpflow_auxvars(ghosted_id),  &
                          material_auxvars(ghosted_id),  &
                          mpflow_parameters(ithrm),        &
                          option,accum_p(istart:iend))

  enddo

  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

end subroutine MpFlowUpdateFixedAccumPatch

! ************************************************************************** !

subroutine MpFlowAccumDerivative(mpflow_auxvar,   &
                             material_auxvar,     &
                             mpflow_parameter,    &
                             option,              &
                             R, J, ifderivative)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: F-M Yuan
  ! Date: 2/13/19
  ! 

  use Option_module
  use Material_Aux_class, only : material_auxvar_type,       &
                                 soil_compressibility_index, &
                                 MaterialCompressSoil
  
  implicit none

  type(mpflow_auxvar_type) :: mpflow_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(mpflow_parameter_type) :: mpflow_parameter
  type(option_type) :: option
  PetscBool, intent(in)  :: ifderivative
  PetscReal, intent(out) :: R(option%nflowdof)
  PetscReal, intent(out) :: J(option%nflowdof,option%nflowdof)

  PetscReal :: compressed_porosity, dcompressed_porosity_dp
  PetscReal :: vol, por, porXvol, dporXvol(option%nflowdof)

  PetscReal :: presl, tc
  PetscReal :: mass, eng
  PetscReal :: dmass, deng

  PetscReal :: u_rock, du_rock(option%nflowdof)
  PetscInt :: ifluid, idof
  PetscReal :: sat, dsat
  PetscReal :: den, dden
  PetscReal :: uh,  duh

  !--------------------------------------------------------------------------------------------
  R(:)   = 0.d0
  J(:,:) = 0.d0

  presl = mpflow_auxvar%pres(LIQ_FLUID)
  tc    = mpflow_auxvar%tc

  vol   = material_auxvar%volume
  if (soil_compressibility_index > 0) then
    call MaterialCompressSoil(material_auxvar,presl, &
                              compressed_porosity,dcompressed_porosity_dp)
    por = compressed_porosity
  else
    por = material_auxvar%porosity_base
    dcompressed_porosity_dp = 0.d0
  endif
  porXvol         = por*vol
  dporXvol(IFDOF1) = dcompressed_porosity_dp*vol
  dporXvol(IFDOF2) = 0.d0

  !------------------------------------------------------------------
  ! NOTE: 'u' or 'h' for soil particles not included in 'mpflow_auxvar' calculation
  !      AND, here assuming that no dependence on 'P' for both mass and energy
  u_rock            = mpflow_parameter%dencpr*(tc + TC2TK)
  du_rock(IFDOF1)   = 0.d0
  du_rock(IFDOF2)   = mpflow_parameter%dencpr
  if (option_flow%isothermal) then
    u_rock         = 0.d0
    du_rock(IFDOF2) = 0.d0
  endif

  !------------------------------------------------------------------
  mass = 0.d0
  eng  = u_rock * (1.d0 - por) * vol
  !
  do ifluid = 1, option%nfluid+1  ! '+1' for immobilized fluid
    if(ifluid>option%nfluid) then
      sat = 1.d0 - mpflow_auxvar%sat(AIR_FLUID) &
                 - mpflow_auxvar%sat(LIQ_FLUID)
    else
      sat = mpflow_auxvar%sat(ifluid)
    endif
    den = mpflow_auxvar%den_kg(ifluid)/FMW_FLUIDS(IFLOWSPEC1)
    uh  = mpflow_auxvar%U(ifluid)
    if (ifluid==LIQ_FLUID) uh = mpflow_auxvar%H(ifluid)    ! needs more thinking here (TODO)

    mass = mass + porXvol*sat*den
    eng  = eng + porXvol*sat*den*uh
  end do

  !residual
  R(IFDOF1) = mass/option_flow%dt
  R(IFDOF2) = eng/option_flow%dt

  !------------------------------------------------------------------
  if (ifderivative .and. .not.option_flow%numerical_derivatives) then
    ! soil matrix (particles): u_rock * vol - u_rock * porXvol
    J(IFDOF2,IFDOF1) = vol*du_rock(IFDOF1) - &
      (u_rock*dporXvol(IFDOF1) + du_rock(IFDOF1)*porXvol)

    J(IFDOF2,IFDOF2) = vol*du_rock(IFDOF2) - &
      (u_rock*dporXvol(IFDOF2) + du_rock(IFDOF2)*porXvol)


    do idof = 1, option%nflowdof
      dmass = 0.d0
      deng  = 0.d0
      do ifluid = 1, option%nfluid+1  !'+1' for immobilized fluid

        if(ifluid>option%nfluid) then
          sat = 1.d0 - mpflow_auxvar%sat(AIR_FLUID) &
                     - mpflow_auxvar%sat(LIQ_FLUID)
          dsat = 0.d0- mpflow_auxvar%D_sat(AIR_FLUID, idof) &
                     - mpflow_auxvar%D_sat(LIQ_FLUID, idof)
        else
          sat = mpflow_auxvar%sat(ifluid)
          dsat = mpflow_auxvar%D_sat(ifluid, idof)
        endif

        den = mpflow_auxvar%den_kg(ifluid)/FMW_FLUIDS(IFLOWSPEC1)
        dden = mpflow_auxvar%D_den_kg(ifluid, idof)/FMW_FLUIDS(IFLOWSPEC1)

        uh  = mpflow_auxvar%U(ifluid)
        duh  = mpflow_auxvar%D_U(ifluid, idof)    ! no 'H' for solid/gas phases (TODO)
        if (ifluid==LIQ_FLUID) then
          uh = mpflow_auxvar%H(ifluid)            ! needs more thinking here (TODO)
          duh = mpflow_auxvar%D_H(ifluid,idof)    ! needs more thinking here (TODO)
        endif

        !mass = porXvol*sat*den
        dmass = dmass + &
                porXvol       *sat  *dden + &
                porXvol       *dsat *den  + &
                dporXvol(idof)*sat  *den

        !eng  = porXvol*sat*den*uh
        deng  = deng + &
                porXvol       *sat  *den  *duh + &
                porXvol       *sat  *dden *uh  + &
                porXvol       *dsat *den  *uh  + &
                dporXvol(idof)*sat  *den  *uh

      end do

      !
      J(IFDOF1, idof) = J(IFDOF1, idof) + dmass
      J(IFDOF2, idof) = J(IFDOF2, idof) + deng
    end do
    !
    J = J/option_flow%dt
  endif

end subroutine MpFlowAccumDerivative

! ************************************************************************** !

subroutine MpFlowAccumulation(mpflow_auxvar,     &
                            material_auxvar,     &
                            mpflow_parameter,    &
                            option,              &
                            Res)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: F-M Yuan
  ! Date: 2/13/19
  ! 

  use Option_module
  use Material_Aux_class, only : material_auxvar_type
  
  implicit none

  type(mpflow_auxvar_type) :: mpflow_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(mpflow_parameter_type) :: mpflow_parameter
  type(option_type) :: option
  PetscReal, intent(out) :: Res(option%nflowdof)

  PetscReal :: J(option%nflowdof,option%nflowdof)
  !
  call MpFlowAccumDerivative(mpflow_auxvar,     &
                           material_auxvar,     &
                           mpflow_parameter,    &
                           option,              &
                           Res, J, PETSC_FALSE)

end subroutine MpFlowAccumulation

! ************************************************************************** !

subroutine MpFlowInternalFluxDerivative(                               &
                              mpflow_auxvar_up,                        &
                              material_auxvar_up, mpflow_parameter_up, &
                              mpflow_auxvar_dn,                        &
                              material_auxvar_dn, mpflow_parameter_dn, &
                              area, dist,                              &
                              option,                                  &
                              v_darcy, fluxe_bulk, fluxe_cond,  Res,   &
                              Jup, Jdn, ifderivative)
  !
  ! Computes the derivatives of fluxes btw two adjacent cells
  !
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
  !
  !
  use Option_module
  use FlowCondition_module

  implicit none

  type(mpflow_auxvar_type) :: mpflow_auxvar_up, mpflow_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(mpflow_parameter_type) :: mpflow_parameter_up, mpflow_parameter_dn
  PetscReal, intent(in) :: area, dist(-1:3)
  type(option_type) :: option
  PetscBool, intent(in)  :: ifderivative
  PetscReal, intent(out) :: v_darcy                                ! unit: m/sec (default)
  PetscReal, intent(out) :: fluxe_bulk, fluxe_cond                 ! unit: MJ/sec (for output)
  PetscReal, intent(out) :: Res(option%nflowdof)
  PetscReal, intent(out) :: Jdn(option%nflowdof, option%nflowdof)
  PetscReal, intent(out) :: Jup(option%nflowdof, option%nflowdof)

  ! locals
  PetscInt :: idof
  PetscReal :: q, dq_up(option%nflowdof), dq_dn(option%nflowdof)
  PetscReal :: qe, dqe_up(option%nflowdof), dqe_dn(option%nflowdof)

!----------------------------------------------------------------------------------------
  Res = 0.d0
  Jup = 0.d0
  Jdn = 0.d0

  !----------------
  q = 0.d0
  dq_up = 0.d0
  dq_dn = 0.d0

  call DarcyFlowDerivative(mpflow_auxvar_up,    mpflow_auxvar_dn,     &
                           material_auxvar_up,  material_auxvar_dn,   &
                           mpflow_parameter_up, mpflow_parameter_dn,  &
                           area, dist,                                &
                           option,                                    &
                           v_darcy,                                   &
                           q,                                         &
                           dq_up, dq_dn, ifderivative)

  !
  Res(IFDOF1)   = q
  if (ifderivative) then
    do idof = 1, option%nflowdof
      Jup(IFDOF1,idof) = dq_up(idof)
      Jdn(IFDOF1,idof) = dq_dn(idof)
    end do
  endif

  !------------------------------------------------
  ! (TODO) vapor transport (air advection, vapor diffusion, and so on)

  !-------------------------------------------------------------------------------------------
  ! thermal conduction
  qe     = 0.d0
  dqe_up = 0.d0
  dqe_dn = 0.d0
  call ConductionDerivative(mpflow_auxvar_up,   mpflow_auxvar_dn,    &
                            material_auxvar_up, material_auxvar_dn,  &
                            area, dist,                              &
                            option,                                  &
                            qe,                                      &
                            dqe_up, dqe_dn, ifderivative)

  !
  fluxe_cond    = qe
  Res(IFDOF2)   = qe
  if (ifderivative) then
    do idof = 1, option%nflowdof
      Jup(IFDOF2,idof) = dqe_up(idof)
      Jdn(IFDOF2,idof) = dqe_dn(idof)
    enddo
  endif

  !--------------------
  ! thermal advection
  qe     = 0.d0
  dqe_up = 0.d0
  dqe_dn = 0.d0
  call AdvectionDerivative(mpflow_auxvar_up,   mpflow_auxvar_dn,    &
                           q, dq_up, dq_dn,                     &
                           option,                              &
                           qe,                                  &
                           dqe_up, dqe_dn, ifderivative)
  !
  fluxe_bulk    = qe
  Res(IFDOF2)   = Res(IFDOF2) + qe
  if (ifderivative) then
    do idof = 1, option%nflowdof
      Jup(IFDOF2,idof) = Jup(IFDOF2,idof) + dqe_up(idof)
      Jdn(IFDOF2,idof) = Jdn(IFDOF2,idof) + dqe_dn(idof)
    end do
  endif

end subroutine MpFlowInternalFluxDerivative

! ************************************************************************** !

subroutine MpFlowInternalFlux(mpflow_auxvar_up,                        &
                              material_auxvar_up, mpflow_parameter_up, &
                              mpflow_auxvar_dn,                        &
                              material_auxvar_dn, mpflow_parameter_dn, &
                              area, dist,                              &
                              option,                                  &
                              v_darcy, fluxe_bulk, fluxe_cond,  Res)
  !
  ! Computes internal fluxes of two-adjacent cells
  !
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
  !
  use Option_module

  implicit none

  type(mpflow_auxvar_type) :: mpflow_auxvar_up, mpflow_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(mpflow_parameter_type) :: mpflow_parameter_up, mpflow_parameter_dn
  type(option_type) :: option
  PetscReal :: dist(-1:3), area
  PetscReal, intent(out) :: Res(option%nflowdof)
  PetscReal, intent(out) :: v_darcy
  PetscReal, intent(out) :: fluxe_bulk, fluxe_cond

  ! locals
  PetscReal :: Jup(option%nflowdof, option%nflowdof)
  PetscReal :: Jdn(option%nflowdof, option%nflowdof)

  !-----------------------------------------------------------------------------
  call MpFlowInternalFluxDerivative(mpflow_auxvar_up,                  &
                              material_auxvar_up, mpflow_parameter_up, &
                              mpflow_auxvar_dn,                        &
                              material_auxvar_dn, mpflow_parameter_dn, &
                              area, dist,                              &
                              option,                                  &
                              v_darcy, fluxe_bulk, fluxe_cond,  Res,   &
                              Jup, Jdn, PETSC_FALSE)

end subroutine MpFlowInternalFlux


! ************************************************************************************************** !

subroutine MpFlowBCFluxDerivative(ibndtype, bc_aux_real_var,           &
                              mpflow_auxvar_dn,                        &
                              material_auxvar_dn, mpflow_parameter_dn, &
                              area, dist,                              &
                              option,                                  &
                              v_darcy, fluxe_bulk, fluxe_cond,  Rdn,   &
                              Jdn, ifderivative)
  !
  ! Computes the derivatives of BCs (Jacobians)
  !
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
  !
  use Option_module
  use FlowCondition_module

  implicit none
  
  PetscInt :: ibndtype(:)
  PetscReal :: bc_aux_real_var(:)
  type(mpflow_auxvar_type) :: mpflow_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(mpflow_parameter_type) :: mpflow_parameter_dn
  PetscReal, intent(in) :: area, dist(-1:3)
  type(option_type) :: option
  PetscBool, intent(in)  :: ifderivative
  PetscReal, intent(out) :: v_darcy                                ! unit: m/sec (default)
  PetscReal, intent(out) :: fluxe_bulk, fluxe_cond                 ! unit: MJ/sec (for output)
  PetscReal, intent(out) :: Rdn(option%nflowdof)
  PetscReal, intent(out) :: Jdn(option%nflowdof, option%nflowdof)

  ! locals
  type(mpflow_auxvar_type) :: mpflow_auxvar_bc

  PetscInt :: idof
  PetscReal :: q, dq_bc(option%nflowdof), dq_dn(option%nflowdof)
  PetscReal :: qe, dqe_bc(option%nflowdof), dqe_dn(option%nflowdof)

  PetscBool :: skip_thermal_conduction
  PetscBool :: skip_mass_flow

!----------------------------------------------------------------------------------------
  skip_thermal_conduction = PETSC_FALSE
  skip_mass_flow = PETSC_FALSE

  Rdn = 0.d0
  Jdn = 0.d0 
  
  ! duplicate BC-connected cell auxvar at first, and later on re-assign specific var accordingly
  !
  call MpFlowAuxVarInit(mpflow_auxvar_bc,option)
  call MpFlowAuxVarCopy(mpflow_auxvar_dn,mpflow_auxvar_bc)

  select case(ibndtype(IFDOF1))
    case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC)
      mpflow_auxvar_bc%pres(LIQ_FLUID) = bc_aux_real_var(IFDOF1) ! all water pressure on soil pore
      mpflow_auxvar_bc%pres(SOLID_PHASE)  = bc_aux_real_var(IFDOF1) ! liq. water pressure on ICE/soil pore
      mpflow_auxvar_bc%pres_fh2o = bc_aux_real_var(IFDOF1)          ! currently this is what to calculated dPressure for flow
      !TODO: needs pressure of air somehow

      ! the following may not be needed
      mpflow_auxvar_bc%D_pres(:,:) = 0.d0
      mpflow_auxvar_bc%D_pc(:) = 0.d0
      mpflow_auxvar_bc%D_sat(:,:) = 0.d0
      mpflow_auxvar_bc%D_den(:,:) = 0.d0
      mpflow_auxvar_bc%D_den_kg(:,:) = 0.d0
      mpflow_auxvar_bc%D_mobility(:,:) = 0.d0
      mpflow_auxvar_bc%D_por(:) = 0.d0
      mpflow_auxvar_bc%D_H(:,:) = 0.d0
      mpflow_auxvar_bc%D_U(:,:) = 0.d0

  end select

  select case(ibndtype(IFDOF2))
    case(DIRICHLET_BC)
      mpflow_auxvar_bc%tc = bc_aux_real_var(IFDOF2)
  end select

  q = 0.d0
  dq_bc = 0.d0
  dq_dn = 0.d0
  select case(ibndtype(IFDOF1))
    case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC)

      call DarcyFlowDerivative(mpflow_auxvar_bc,    mpflow_auxvar_dn,    &
                               material_auxvar_dn,  material_auxvar_dn,  &
                               mpflow_parameter_dn, mpflow_parameter_dn, &
                               area, dist,                               &
                               option,                                   &
                               v_darcy,                                  &
                               q,                                        &
                               dq_bc, dq_dn, ifderivative)

      if (ibndtype(IFDOF1) == HYDROSTATIC_SEEPAGE_BC) then
        ! ONE-way flow, i.e. flow outward ONLY
        ! + from up (bc) -> dn; - from dn -> bc (up)
        if (q > 0.d0) then
          q     = 0.d0
          dq_bc = 0.d0
          dq_dn = 0.d0
        endif
        skip_thermal_conduction = PETSC_TRUE
      endif

      
    case(NEUMANN_BC)
      v_darcy = bc_aux_real_var(IFDOF1)
      q = v_darcy * area
      dq_bc = 0.d0
      dq_dn = 0.d0

    case(ZERO_GRADIENT_BC)
      q     = 0.d0
      dq_bc = 0.d0
      dq_dn = 0.d0
      ! if flux-type BC for T, the fluid is totally energy form without mass
      !(may not be needed, but just in case)
      if(ibndtype(IFDOF2) == NEUMANN_BC) skip_mass_flow = PETSC_TRUE

    case default
      option%io_buffer = 'BC type for Liq. water flow: "' // trim(GetSubConditionName(ibndtype(IFDOF1))) // &
        '" not implemented in Flow mode.'
      call printErrMsg(option)

  end select

  !
  Rdn(IFDOF1)   = q
  if (ifderivative) then
    do idof = 1, option%nflowdof
      Jdn(IFDOF1,idof) = dq_dn(idof)
    end do
  endif

  !------------------------------------------------
  ! (TODO) vapor transport (air advection, vapor diffusion, and so on)

  !-------------------------------------------------------------------------------------------
  ! BC Thermal terms
  qe     = 0.d0
  dqe_bc = 0.d0
  dqe_dn = 0.d0
  select case(ibndtype(IFDOF2))
    ! conduction
    case(DIRICHLET_BC)

      if (skip_thermal_conduction) then
        ! skip thermal conducton for some specific boundary condition
        ! (e.g. BC pressure is below reference pressure such as river stage is below cell center).

        qe     = 0.d0
        dqe_bc = 0.d0
        dqe_dn = 0.d0
      else

        call ConductionDerivative(mpflow_auxvar_bc, mpflow_auxvar_dn,  &
                              material_auxvar_dn, material_auxvar_dn,  &
                              area, dist,                              &
                              option,                                  &
                              qe,                                      &
                              dqe_bc, dqe_dn, ifderivative)

      endif
      !
    case(NEUMANN_BC)
      qe = bc_aux_real_var(IFDOF2)*area
      !
    case(ZERO_GRADIENT_BC)
      qe = 0.d0
      !
    case default
      option%io_buffer = 'BC type for Thermal condition: "' // trim(GetSubConditionName(ibndtype(IFDOF2))) // &
        '" not implemented in Flow mode.'
      call printErrMsg(option)

  end select
  !
  fluxe_cond    = qe
  Rdn(IFDOF2)   = Rdn(IFDOF2) + qe
  if (ifderivative) then
    do idof = 1, option%nflowdof
      Jdn(IFDOF2,idof) = dqe_dn(idof)
    end do
  endif

  !--------------------
  ! thermal advection (after conduction term, so that having correct BC vars)
  call AdvectionDerivative(mpflow_auxvar_bc, mpflow_auxvar_dn,  &
                           q, dq_bc, dq_dn,                     &
                           option,                              &
                           qe,                                  &
                           dqe_bc, dqe_dn, ifderivative)
  !
  fluxe_bulk    = qe
  Rdn(IFDOF2)   = Rdn(IFDOF2) + qe
  if (ifderivative) then
    do idof = 1, option%nflowdof
      Jdn(IFDOF2,idof) = Jdn(IFDOF2,idof) + dqe_dn(idof)
    end do
  endif

  call MpFlowAuxVarStrip(mpflow_auxvar_bc)

end subroutine MpFlowBCFluxDerivative

! ************************************************************************** !

subroutine MpFlowBCFlux(ibndtype, bc_aux_real_var,           &
                      mpflow_auxvar_dn,                      &
                      material_auxvar_dn,                    &
                      mpflow_parameter_dn,                   &
                      area, dist,                            &
                      option,                                &
                      v_darcy, fluxe_bulk, fluxe_cond, Res)
  !
  ! Computes the  boundary flux terms for the residual
  ! 
  ! Author: F-M. Yuan
  ! Date: 02/13/2019
  ! 
  use Option_module
 
  implicit none
  
  PetscInt, intent(in) :: ibndtype(:)
  PetscReal, intent(in):: bc_aux_real_var(:) ! from aux_real_var array
  type(mpflow_auxvar_type), intent(in) :: mpflow_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(mpflow_parameter_type)   :: mpflow_parameter_dn
  type(option_type) :: option
  PetscReal :: dist(-1:3), area
  PetscReal, intent(out) :: Res(option%nflowdof)
  PetscReal, intent(out) :: v_darcy
  PetscReal, intent(out) :: fluxe_bulk, fluxe_cond
  
  ! locals
  PetscReal :: Jdn(option%nflowdof, option%nflowdof)

  !-----------------------------------------------------------------------------
  call MpFlowBCFluxDerivative(ibndtype, bc_aux_real_var,               &
                              mpflow_auxvar_dn,                        &
                              material_auxvar_dn,                      &
                              mpflow_parameter_dn,                     &
                              area, dist,                              &
                              option,                                  &
                              v_darcy, fluxe_bulk, fluxe_cond, Res,    &
                              Jdn, PETSC_FALSE)

end subroutine MpFlowBCFlux

! ********************************************************************************************************* !

subroutine MpFlowSetPlotVariables(realization,list)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 
  
  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module
  use Material_Aux_class
  use Option_module

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

    name = 'Liquid Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                LIQUID_SATURATION)
                                
    name = 'Liquid Density'
    units = 'kg/m^3'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_DENSITY)

    name = 'Liquid Energy'
    units = 'kJ/mol'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_ENERGY)

    name = 'Liquid Viscosity'
    units = 'Pa.s'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_VISCOSITY)

    name = 'Liquid Mobility'
    units = '1/Pa.s'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_MOBILITY)
  
  endif
  
  if (list%energy_vars) then

    name = 'Temperature'
    units = 'C'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                TEMPERATURE)
  
  endif

  name = 'Gas Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
          GAS_PRESSURE)

  name = 'Gas Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
          GAS_SATURATION)

  name = 'Gas Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
          GAS_DENSITY)

  name = 'Ice Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
        ICE_SATURATION)

  name = 'Ice Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
        ICE_DENSITY)
        

  if (soil_compressibility_index > 0) then
  
    name = 'Transient Porosity'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 POROSITY_CURRENT)
                                 
  endif

end subroutine MpFlowSetPlotVariables


! ************************************************************************** !
function MpFlowInitGuessCheck(xx, option)
  !
  ! Checks if the initial guess is valid.
  ! Note: Only implemented for DALL_AMICO formulation.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/04/2014
  !
  use Option_module

  Vec :: xx
  type(option_type), pointer :: option

  PetscInt :: MpFlowInitGuessCheck
  PetscInt :: idx
  PetscReal :: pres_min, pres_max
  PetscReal :: temp_min, temp_max
  PetscInt :: ipass, ipass0
  PetscErrorCode :: ierr

  ipass = 1

  call VecStrideMin(xx,ZERO_INTEGER,idx,pres_min,ierr)
  call VecStrideMin(xx,ONE_INTEGER ,idx,temp_min,ierr)
  call VecStrideMax(xx,ZERO_INTEGER,idx,pres_max,ierr)
  call VecStrideMax(xx,ONE_INTEGER ,idx,temp_max,ierr)

  if (pres_min < -1.d10 .or. pres_min > 1.d10 .or. &
      temp_min < -100.d0 .or. temp_max > 100.d0) then
      ipass = -1
  endif

   call MPI_Barrier(option%mycomm,ierr)
   if (option%mycommsize>1)then
      call MPI_Allreduce(ipass,ipass0,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                         option%mycomm,ierr)
      if (ipass0 < option%mycommsize) ipass=-1
   endif
   MpFlowInitGuessCheck = ipass

end function MpFlowInitGuessCheck


! ************************************************************************** !

subroutine MpFlowDestroy(patch)
  ! Author: ???
  ! Date: 02/14/08
  ! 

  use Patch_module

  implicit none
  
  type(patch_type) :: patch
  
  ! need to free array in aux vars
  call MpFlowAuxDestroy(patch%aux%MpFlow)

end subroutine MpFlowDestroy

! ************************************************************************** !

end module MpFlow_module
