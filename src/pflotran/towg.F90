module TOWG_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_TOWG_Aux_module
  use AuxVars_TOWG_module
  use Global_Aux_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#define CONVECTION
#define DIFFUSION
#define LIQUID_DIFFUSION
#define CONDUCTION

#define GLOBALWORKERS

!#define DEBUG_TOWG_FILEOUTPUT
!#define DEBUG_TOWG_FLUXES  

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

#ifdef DEBUG_TOWG_FILEOUTPUT
  PetscInt, parameter :: debug_unit = 87
  PetscInt, parameter :: debug_info_unit = 86
  character(len=MAXWORDLENGTH) :: debug_filename
  PetscInt :: debug_flag = 0
  PetscInt :: debug_iteration_count
  PetscInt :: debug_timestep_cut_count
  PetscInt :: debug_timestep_count
#endif

#ifdef GLOBALWORKERS
  PetscReal, allocatable, dimension(:) :: D_den_kg_ave_up,D_den_kg_ave_dn
  PetscReal, allocatable, dimension(:) :: D_den_ave_up,D_den_ave_dn
  PetscReal, allocatable, dimension(:) :: D_delta_presure_up,D_delta_presure_dn
  PetscReal, allocatable, dimension(:) :: D_mobility_up,D_mobility_dn
  PetscReal, allocatable, dimension(:) :: D_uH_up,D_uH_dn
  PetscReal, allocatable, dimension(:) :: D_v_darcy_up,D_v_darcy_dn
  PetscReal, allocatable, dimension(:) :: D_q_up,D_q_dn
  PetscReal, allocatable, dimension(:) :: D_mole_flux_up,D_mole_flux_dn
  PetscReal, allocatable, dimension(:) :: D_xmf_up,D_xmf_dn

  PetscReal, allocatable, dimension(:) :: D_sat_liquid_up,D_sat_liquid_dn,D_k_eff_up,D_k_eff_dn
  PetscReal, allocatable, dimension(:) :: D_k_eff_ave_up,D_k_eff_ave_dn,D_delta_temp_up,D_delta_temp_dn
  PetscReal, allocatable, dimension(:) :: D_worker1,D_worker2
#endif


  !pointing to null() function
  procedure(TOWGUpdateAuxVarsDummy), pointer :: TOWGUpdateAuxVars => null()
  procedure(TOWGAccumulationDummy), pointer :: TOWGAccumulation => null()
  procedure(TOWGComputeMassBalanceDummy), pointer :: &
                                           TOWGComputeMassBalance => null()
  procedure(TOWGFluxDummy), pointer :: TOWGFlux => null()
  procedure(TOWGBCFluxDummy), pointer :: TOWGBCFlux => null()
  procedure(TOWGSrcSinkDummy), pointer :: TOWGSrcSink => null()
  procedure(TOWGCheckUpdatePreDummy), pointer :: TOWGCheckUpdatePre => null()
  procedure(TOWGMaxChangeDummy), pointer :: TOWGMaxChange => null()

  abstract interface

    subroutine TOWGUpdateAuxVarsDummy(realization,update_state)
      use Realization_Subsurface_class  
      implicit none
      type(realization_subsurface_type) :: realization
      PetscBool :: update_state
    end subroutine TOWGUpdateAuxVarsDummy

    subroutine TOWGAccumulationDummy(auxvar,global_auxvar,material_auxvar, &
                                     soil_heat_capacity,option,Res,debug_cell,&
                                     j,analytical_derivatives)
      use AuxVars_TOWG_module
      use Global_Aux_module
      use Option_module
      use Material_module
      use Material_Aux_class
      implicit none
      class(auxvar_towg_type) :: auxvar
      type(global_auxvar_type) :: global_auxvar
      class(material_auxvar_type) :: material_auxvar
      PetscReal :: soil_heat_capacity
      type(option_type) :: option
      PetscReal :: Res(option%nflowdof) 
      PetscBool :: debug_cell
      PetscReal :: j(1:option%nflowspec,1:option%nflowdof) 
      PetscBool :: analytical_derivatives
    end subroutine TOWGAccumulationDummy

    subroutine TOWGComputeMassBalanceDummy(realization,mass_balance)
      use Realization_Subsurface_class 
      implicit none
      type(realization_subsurface_type) :: realization
      PetscReal :: mass_balance(realization%option%nflowspec, &
                            realization%option%nphase)
    end subroutine TOWGComputeMassBalanceDummy

    subroutine TOWGFluxDummy(auxvar_up,global_auxvar_up, &
                             material_auxvar_up, &
                             thermal_conductivity_up, &
                             auxvar_dn,global_auxvar_dn, &
                             material_auxvar_dn, &
                             thermal_conductivity_dn, &
                             area, dist, towg_parameter, &
                             option,v_darcy,Res, &
                             debug_connection, &
                             jup,jdn,analytical_derivatives)
      use PM_TOWG_Aux_module
      use AuxVars_TOWG_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_class
      implicit none
      class(auxvar_towg_type) :: auxvar_up, auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
      type(option_type) :: option
      PetscReal :: v_darcy(option%nphase)
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(towg_parameter_type) :: towg_parameter
      PetscReal :: thermal_conductivity_dn(2)
      PetscReal :: thermal_conductivity_up(2)
      PetscReal :: Res(option%nflowdof)
      PetscBool :: debug_connection
      PetscReal, dimension(1:option%nflowdof,1:option%nflowdof) :: jup,jdn
      PetscBool :: analytical_derivatives
    end subroutine TOWGFluxDummy

    subroutine TOWGBCFluxDummy(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                               auxvar_up,global_auxvar_up, &
                               auxvar_dn,global_auxvar_dn, &
                               material_auxvar_dn, &
                               thermal_conductivity_dn, &
                               area,dist,towg_parameter, &
                               option,v_darcy,Res,debug_connection,&
                               Jdn,analytical_derivatives)
      use PM_TOWG_Aux_module
      use AuxVars_TOWG_module
      use Global_Aux_module
      use Option_module                              
      use Material_Aux_class
      implicit none
      class(auxvar_towg_type) :: auxvar_up, auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      class(material_auxvar_type) :: material_auxvar_dn
      type(option_type) :: option
      PetscReal :: bc_auxvars(:)
      PetscReal :: v_darcy(option%nphase), area
      type(towg_parameter_type) :: towg_parameter
      PetscReal :: dist(-1:3)
      PetscReal :: Res(1:option%nflowdof)
      PetscInt :: ibndtype(1:option%nflowdof)
      PetscInt :: bc_auxvar_mapping(TOWG_MAX_INDEX)
      PetscReal :: thermal_conductivity_dn(2)
      PetscBool :: debug_connection
      PetscReal, dimension(1:option%nflowdof,1:option%nflowdof) :: Jdn
      PetscBool :: analytical_derivatives
    end subroutine TOWGBCFluxDummy

    subroutine TOWGSrcSinkDummy(option,src_sink_condition, auxvar, &
                            global_auxvar,ss_flow_vol_flux,scale,Res,&
                            j,analytical_derivatives)
      use Option_module
      use AuxVars_TOWG_module
      use Global_Aux_module
      use Condition_module  
      implicit none
      type(option_type) :: option
      type(flow_towg_condition_type), pointer :: src_sink_condition
      class(auxvar_towg_type) :: auxvar
      type(global_auxvar_type) :: global_auxvar
      PetscReal :: ss_flow_vol_flux(option%nphase)
      PetscReal :: scale  
      PetscReal :: Res(option%nflowdof)
      PetscReal, dimension(1:option%nflowdof,1:option%nflowdof) :: j
      PetscBool :: analytical_derivatives
    end subroutine TOWGSrcSinkDummy

    subroutine TOWGCheckUpdatePreDummy(line_search,X,dX,changed,realization, &
                                       max_it_before_damping,damping_factor, &
                                       max_pressure_change,ierr)
      use Realization_Subsurface_class
#include "petsc/finclude/petscsnes.h"
      use petscsnes
      implicit none
      SNESLineSearch :: line_search
      Vec :: X
      Vec :: dX
      PetscBool :: changed
      PetscErrorCode :: ierr
      type(realization_subsurface_type) :: realization
      PetscInt :: max_it_before_damping
      PetscReal :: damping_factor
      PetscReal :: max_pressure_change
    end subroutine TOWGCheckUpdatePreDummy

    subroutine TOWGMaxChangeDummy(realization,max_change_ivar, &
                                  max_change_isubvar,max_pressure_change, &
                                  max_xmol_change,max_saturation_change, &
                                  max_temperature_change)
      use Realization_Subsurface_class
      implicit none
      class(realization_subsurface_type), pointer :: realization
      PetscInt :: max_change_ivar(:)
      PetscInt :: max_change_isubvar(:)
      PetscReal :: max_pressure_change
      PetscReal :: max_xmol_change
      PetscReal :: max_saturation_change
      PetscReal :: max_temperature_change
    end subroutine TOWGMaxChangeDummy

  end interface 

  public :: TOWGSetup,&
            TOWGUpdateAuxVars, &
            TOWGUpdateSolution, &
            TOWGMapBCAuxVarsToGlobal, &
            TOWGInitializeTimestep, &
            TOWGComputeMassBalance, &
            TOWGResidual, &
            TOWGJacobian, &
            TOWGCheckUpdatePre, &
            TOWGTimeCut, &
            TOWGMaxChange, &
            TOWGDestroy

contains

! ************************************************************************** !

function checkBlackOilCIP(iphase,icomp,is_oil_in_oil,is_gas_in_oil,option)

!------------------------------------------------------------------------------
! Used in black oil, returns .true. if phase iphase contains component icomp
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Oct 2017
!------------------------------------------------------------------------------

  use Option_module, only : option_type

  PetscBool                    :: checkBlackOilCIP
  PetscInt,  intent(in )       :: iphase
  PetscInt,  intent(in )       :: icomp
  PetscBool, intent(out)       :: is_oil_in_oil,is_gas_in_oil
  type(option_type),intent(in) :: option
  PetscBool                    :: componentInPhase

!  Default return values

  componentInPhase = PETSC_FALSE
  is_oil_in_oil    = PETSC_FALSE
  is_gas_in_oil    = PETSC_FALSE

! disgas is special case

  if( iphase == option%oil_phase ) then
    if( icomp == option%oil_phase ) is_oil_in_oil = PETSC_TRUE
    if( icomp == option%gas_phase ) is_gas_in_oil = PETSC_TRUE
  endif

! OK if phase and component match or dissolved gas

  if( iphase == icomp .or. is_gas_in_oil )  componentInPhase = PETSC_TRUE

  checkBlackOilCIP = componentInPhase

end function checkBlackOilCIP

! ************************************************************************** !

subroutine TOWGSetup(realization)
  ! 
  ! Creates arrays for TOWG auxiliary variables
  ! 
  ! Author: Paolo Orsini and Dave Ponting (OGS)
  ! Date: 11/07/16
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  !use Fluid_module
  use Material_Aux_class
  use Output_Aux_module
  use AuxVars_Flow_module

  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(material_parameter_type), pointer :: material_parameter
  type(output_variable_list_type), pointer :: list

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscInt :: num_bc_connection, num_ss_connection
  PetscInt :: i, idof, count
  PetscBool :: error_found
  PetscInt :: flag(10)

  class(material_auxvar_type), pointer :: material_auxvars(:)
  !type(fluid_property_type), pointer :: cur_fluid_property
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  
  patch%aux%TOWG => TOWGAuxCreate(option)

  towg_analytical_derivatives = .not. option%flow%numerical_derivatives
  towg_analytical_derivatives_compare = option%flow%numerical_derivatives_compare

  towg_dcomp_tol = flow_aux_debug_tol
  towg_dcomp_reltol = flow_aux_debug_reltol

#ifdef GLOBALWORKERS
  if (towg_analytical_derivatives) then
    allocate(D_den_kg_ave_up (1:option%nflowdof))
    allocate(D_den_kg_ave_dn (1:option%nflowdof))
    allocate(D_den_ave_up (1:option%nflowdof))
    allocate(D_den_ave_dn (1:option%nflowdof))
    allocate(D_delta_presure_up (1:option%nflowdof))
    allocate(D_delta_presure_dn (1:option%nflowdof))
    allocate( D_mobility_up (1:option%nflowdof))
    allocate(D_mobility_dn (1:option%nflowdof))
    allocate(D_uH_up (1:option%nflowdof))
    allocate(D_uH_dn (1:option%nflowdof))
    allocate(D_v_darcy_up (1:option%nflowdof))
    allocate(D_v_darcy_dn (1:option%nflowdof))
    allocate(D_q_up (1:option%nflowdof))
    allocate(D_q_dn (1:option%nflowdof))
    allocate(D_mole_flux_up (1:option%nflowdof))
    allocate(D_mole_flux_dn (1:option%nflowdof))
    allocate(D_xmf_up (1:option%nflowdof))
    allocate(D_xmf_dn (1:option%nflowdof))

    allocate( D_sat_liquid_up  (1:option%nflowdof))
    allocate( D_sat_liquid_dn (1:option%nflowdof))
    allocate(D_k_eff_up  (1:option%nflowdof))
    allocate(D_k_eff_dn  (1:option%nflowdof))
    allocate(D_k_eff_ave_up (1:option%nflowdof))
    allocate(D_k_eff_ave_dn (1:option%nflowdof))
    allocate(D_delta_temp_up (1:option%nflowdof))
    allocate(  D_delta_temp_dn (1:option%nflowdof))
    allocate( D_worker1 (1:option%nflowdof))
    allocate(D_worker2 (1:option%nflowdof))

  endif
#endif

  ! ensure that material properties specific to this module are properly
  ! initialized
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE

  if (minval(material_parameter%soil_heat_capacity(:)) < 0.d0) then
    option%io_buffer = 'Non-initialized soil heat capacity.'
    call PrintMsg(option)
    error_found = PETSC_TRUE
  endif
  if (minval(material_parameter%soil_thermal_conductivity(:,:)) < 0.d0) then
    option%io_buffer = 'Non-initialized soil thermal conductivity.'
    call PrintMsg(option)
    error_found = PETSC_TRUE
  endif
  
  material_auxvars => patch%aux%Material%auxvars
  flag = 0

  !TODO(geh): change to looping over ghosted ids once the legacy code is 
  !           history and the communicator can be passed down.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (material_auxvars(ghosted_id)%volume < 0.d0 .and. flag(1) == 0) then
      flag(1) = 1
      option%io_buffer = 'Non-initialized cell volume.'
      call PrintMsg(option)
    endif
    if (material_auxvars(ghosted_id)%porosity_base < 0.d0 .and. &
        flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'Non-initialized porosity.'
      call PrintMsg(option)
    endif
    if (material_auxvars(ghosted_id)%tortuosity < 0.d0 .and. flag(3) == 0) then
      flag(3) = 1
      option%io_buffer = 'Non-initialized tortuosity.'
      call PrintMsg(option)
    endif
    if (material_auxvars(ghosted_id)%soil_particle_density < 0.d0 .and. &
        flag(4) == 0) then
      flag(4) = 1
      option%io_buffer = 'Non-initialized soil particle density.'
      call PrintMsg(option)
    endif
    if (minval(material_auxvars(ghosted_id)%permeability) < 0.d0 .and. &
        flag(5) == 0) then
      option%io_buffer = 'Non-initialized permeability.'
      call PrintMsg(option)
      flag(5) = 1
    endif
  enddo

  if (error_found .or. maxval(flag) > 0) then
    option%io_buffer = 'Material property errors found in TOWGSetup.'
    call PrintErrMsg(option)
  endif

  num_bc_connection = &
               CouplerGetNumConnectionsInList(patch%boundary_condition_list)

  num_ss_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)

  call patch%aux%TOWG%Init(grid,num_bc_connection,num_ss_connection,option)


  !initialise here diffusion coefficient when supporting the oil and gas
  ! diffusione terms
  ! initialize parameters
  !cur_fluid_property => realization%fluid_properties
  !do 
  !  if (.not.associated(cur_fluid_property)) exit
  !  patch%aux%TOWGl%parameter% &
  !    diffusion_coefficient(cur_fluid_property%phase_id) = &
  !      cur_fluid_property%diffusion_coefficient
  !  cur_fluid_property => cur_fluid_property%next
  !enddo  
  ! check whether diffusion coefficients are initialized.
  ! check initialisation of diffusiont coefficient phase by phase as for 
  ! example below
  !if (Uninitialized(patch%aux%TOWG%parameter% &
  !    diffusion_coefficient(option%liquid_phase))) then
  !  option%io_buffer = &
  !    UninitializedMessage('Liquid phase diffusion coefficient','')
  !  call PrintErrMsg(option)
  !endif
  !

  list => realization%output_option%output_snap_variable_list
  call TOWGSetPlotVariables(list)
  list => realization%output_option%output_obs_variable_list
  call TOWGSetPlotVariables(list) 

  ! convergence criteria to be chosen (can use TOUGH or TOWG type)
  ! set up here tough convergence creteria if needed.

!------------------------------------------------------------------------------
!  TOWG functions that can vary with the miscibility model
!------------------------------------------------------------------------------

  select case(towg_miscibility_model)
    case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF)
      TOWGCheckUpdatePre => TOWGImsTLCheckUpdatePre
      TOWGMaxChange => TOWGImsTLMaxChange
      TOWGSrcSink   => TOWGImsTLSrcSink
      select case(towg_miscibility_model)
        case(TOWG_IMMISCIBLE)
          call TOWGImsAuxVarComputeSetup()
        case(TOWG_TODD_LONGSTAFF)
          call TOWGTLAuxVarComputeSetup()
      end select
    case(TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
! No trap on neg Sg, avoid dampling Pb as Sg, so special black oil/TL4P version
       TOWGCheckUpdatePre => TOWGBlackOilCheckUpdatePre
! Detailed lookup needed for the oil source case, so special black oil/TL4P version
       TOWGSrcSink => TOWGBOSrcSink
! Cases in which black oil and TL4P are different
       select case(towg_miscibility_model)
          case(TOWG_BLACK_OIL)
! Must convert Pbub changes to eff. satn. change, so special black oil version
            TOWGMaxChange => TOWGBOMaxChange
            call TOWGBlackOilAuxVarComputeSetup()
          case(TOWG_SOLVENT_TL)
! Pbub changes and extra saturation, so special solvent version
            TOWGMaxChange => TOWGTL4PMaxChange
            call TL4PAuxVarComputeSetup()
       end select
    case default
       option%io_buffer = 'TOWGSetup: mode not supported.'
       call PrintErrMsg(option)
  end select

!------------------------------------------------------------------------------
! TOWG functions that do not vary with the miscibility model
!------------------------------------------------------------------------------

  select case(towg_miscibility_model)
    case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF,TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
       TOWGAccumulation => TOWGImsTLBOAccumulation
       TOWGUpdateAuxVars => TOWGImsTLBOUpdateAuxVars
       TOWGComputeMassBalance => TOWGImsTLBOComputeMassBalance
       TOWGFlux => TOWGImsTLBOFlux
       TOWGBCFlux => TOWGImsTLBOBCFlux
    case default
      option%io_buffer = 'TOWGSetup: mode not supported.'
      call PrintErrMsg(option)
  end select

#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_flag = 0
  debug_iteration_count = 0
  debug_timestep_cut_count = 0
  debug_timestep_count = 0
  ! create new file
  open(debug_info_unit, file='debug_towg_info.txt', action="write", &
       status="unknown")
  write(debug_info_unit,*) 'type timestep cut iteration'
  close(debug_info_unit)
#endif  

end subroutine TOWGSetup

! ************************************************************************** !

#ifdef GLOBALWORKERS
function CheckWorkersAllocated()

  implicit none

  PetscBool :: CheckWorkersAllocated

  CheckWorkersAllocated = PETSC_TRUE

  if (.NOT. allocated(D_den_kg_ave_up))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_den_kg_ave_dn))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_den_ave_up))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_den_ave_dn))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_delta_presure_up))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_delta_presure_dn))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_mobility_up))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_mobility_dn))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_uH_up))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_uH_dn))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_v_darcy_up))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_v_darcy_dn))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_q_up))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_q_dn))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_mole_flux_up))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_mole_flux_dn))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_xmf_up))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif
  if (.NOT. allocated(D_xmf_dn))  then
    CheckWorkersAllocated = PETSC_FALSE; return
  endif

end function CheckWorkersAllocated
#endif

! ************************************************************************** !

subroutine TOWGImsTLBOUpdateAuxVars(realization,update_state)
  ! 
  ! Updates the auxiliary variables associated with the TOWG_IMS problem
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/30/16
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Material_Aux_class
  use EOS_Water_module
  use Saturation_Function_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  PetscBool :: update_state
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  class(pm_towg_aux_type), pointer :: towg
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)  

  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, natural_id
  PetscInt :: ghosted_start, ghosted_end
  PetscInt :: iphasebc, iphase
  PetscInt :: offset
  PetscInt :: istate
  PetscReal :: gas_pressure, capillary_pressure, liquid_saturation
  PetscReal :: saturation_pressure, temperature
  PetscInt :: real_index, variable
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)

!#define DEBUG_AUXVARS
#ifdef DEBUG_AUXVARS
  character(len=MAXWORDLENGTH) :: word
  PetscInt, save :: icall = 0
#endif

  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  towg => patch%aux%TOWG
  !gen_auxvars => patch%aux%General%auxvars
  !gen_auxvars_bc => patch%aux%General%auxvars_bc

  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
    
!------------------------------------------------------------------------------
! Extract solution
!------------------------------------------------------------------------------

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

#ifdef DEBUG_AUXVARS
  icall = icall + 1
  write(word,*) icall
  word = 'towgaux' // trim(adjustl(word))
#endif

!------------------------------------------------------------------------------
! Loop over cells
!------------------------------------------------------------------------------

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ghosted_end = ghosted_id*option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    !TOWG_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = TOWG_UPDATE_FOR_FIXED_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
 
    !to replace with TOWGImsAuxVarCompute (there is no need for pointer here)
    call TOWGAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                        towg%auxvars(ZERO_INTEGER,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        material_auxvars(ghosted_id), &
                        patch%characteristic_curves_array( &
                        patch%sat_func_id(ghosted_id))%ptr, &
                       natural_id, &
                       option)
  enddo

!------------------------------------------------------------------------------
! Now boundary conditions
!------------------------------------------------------------------------------

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
      offset = (ghosted_id-1)*option%nflowdof
      if (patch%imat(ghosted_id) <= 0) cycle

      xxbc(:) = xx_loc_p(offset+1:offset+option%nflowdof)
      istate = boundary_condition%flow_aux_int_var(TOWG_STATE_INDEX,iconn)
      if (istate == TOWG_NULL_STATE) then !this is applied to flux (Neumann) conditions
        istate = global_auxvars(ghosted_id)%istate !look into the state of down cell
        select case(istate)
! Two states possible: TOWG_THREE_PHASE_STATE,TOWG_LIQ_OIL_STATE
          case(TOWG_THREE_PHASE_STATE,TOWG_LIQ_OIL_STATE)
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(HYDROSTATIC_BC)
                  real_index = &
                    boundary_condition% &
                    flow_aux_mapping(towg_dof_to_primary_variable(idof,istate))
                  xxbc(idof) = &
                    boundary_condition%flow_aux_real_var(real_index,iconn)
                case(DIRICHLET_BC)
                  variable = towg_dof_to_primary_variable(idof,istate)
                  select case(variable)
                    ! for oil pressure dof
                    case(TOWG_OIL_PRESSURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        option%io_buffer = 'TOWG Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" needs oil pressure defined.'
                        call PrintErrMsg(option)
                      endif
                    ! for oil saturation dof
                    case(TOWG_OIL_SATURATION_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        ! should be able to use oil saturation in DN cell
                        !option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                        !  trim(boundary_condition%flow_condition%name) // &
                        !  '" needs oil saturation defined.'
                        !call PrintErrMsg(option)
                      endif
                    case(TOWG_GAS_SATURATION_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        ! should be able to use oil saturation in DN cell
                        !option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                        !  trim(boundary_condition%flow_condition%name) // &
                        !  '" needs gas saturation defined.'
                        !call PrintErrMsg(option)
                      endif
                    case(TOWG_TEMPERATURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        option%io_buffer = 'TOWG Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" needs temperature defined.'
                        call PrintErrMsg(option)
                      endif
                  end select
                case(NEUMANN_BC)
                case default
                  option%io_buffer = 'Unknown BC type in TOWGUpdateAuxVars().'
                  call PrintErrMsg(option)
              end select
            enddo  
        end select
      else
        ! we do this for all BCs; Neumann bcs will be set later
        do idof = 1, option%nflowdof
          real_index = boundary_condition% &
              flow_aux_mapping(towg_dof_to_primary_variable(idof,istate))
          if (real_index > 0) then
            xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
          else
            option%io_buffer = 'Error setting up boundary ' // &
                                'condition in TOWGUpdateAuxVars'
            call PrintErrMsg(option)
          endif
        enddo
      endif
          
!--Now do the auxvar setup for the boundary condition--------------------------

      ! set this based on data given 
      global_auxvars_bc(sum_connection)%istate = istate
      ! TOWG_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = TOWG_UPDATE_FOR_BOUNDARY
      call TOWGAuxVarCompute(xxbc,towg%auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%sat_func_id(ghosted_id))%ptr, &
                                natural_id, &
                                option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

!------------------------------------------------------------------------------
! Restore state and solution
!------------------------------------------------------------------------------

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%TOWG%auxvars_up_to_date = PETSC_TRUE

end subroutine TOWGImsTLBOUpdateAuxVars

! ************************************************************************** !

subroutine TOWGInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/07/16 
  ! 

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  call TOWGUpdateFixedAccum(realization)
  
#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_flag = 0
  if (.true.) then
    debug_iteration_count = 0
    debug_flag = 1
  endif
  debug_iteration_count = 0
#endif

end subroutine TOWGInitializeTimestep

! ************************************************************************** !

subroutine TOWGUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/06/16 
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Well_Data_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  class(well_data_type), pointer :: well_data
  type(well_data_list_type),pointer :: well_data_list

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  class(pm_towg_aux_type), pointer :: towg
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr
  PetscReal :: dt
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  towg => patch%aux%TOWG
  global_auxvars => patch%aux%Global%auxvars

! Loop over well_data wells if present

  dt = option%flow_dt

  if (WellDataGetFlag()) then
    well_data_list => realization%well_data
    well_data => well_data_list%first

    do
      if (.not.associated(well_data)) exit
      call well_data%DoUpdate(dt,option)
      well_data => well_data%next
    enddo

    call FindGroupRates (well_data_list)
    call FindGroupTotals(well_data_list)

  endif

  if (realization%option%compute_mass_balance_new) then
    call TOWGUpdateMassBalance(realization)
  endif
  
  do ghosted_id = 1, grid%ngmax
    towg%auxvars(ZERO_INTEGER,ghosted_id)%istate_store(TOWG_PREV_TS) = &
      global_auxvars(ghosted_id)%istate
  enddo

#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_iteration_count = 0
  debug_timestep_cut_count = 0
  debug_timestep_count = debug_timestep_count + 1
#endif 
    
end subroutine TOWGUpdateSolution

! ************************************************************************** !

subroutine TOWGTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/30/16
  ! 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Patch_module
  use Discretization_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(pm_towg_aux_type), pointer :: towg
  
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars
  towg => patch%aux%TOWG

  ! restore stored state
  do ghosted_id = 1, grid%ngmax
    global_auxvars(ghosted_id)%istate = &
      towg%auxvars(ZERO_INTEGER,ghosted_id)%istate_store(TOWG_PREV_TS)
  enddo

#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_timestep_cut_count = debug_timestep_cut_count + 1
#endif 

  call TOWGInitializeTimestep(realization)  

end subroutine TOWGTimeCut

! ************************************************************************** !

subroutine TOWGZeroMassBalanceDelta(realization)
  ! 
  ! Zeros mass balance delta array  
  ! PO: identical for many flow modes (Genral, Toil_Ims, TOWG), where can it 
  !     be located to be shared?? flow_mode_common.F90 ?
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/08/16 
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

  do iconn = 1, patch%aux%TOWG%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%TOWG%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine TOWGZeroMassBalanceDelta

! ************************************************************************** !

subroutine TOWGUpdateMassBalance(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/06/16 
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
  PetscInt :: icomp

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%TOWG%num_aux_bc
    do icomp = 1, option%nflowspec
      global_auxvars_bc(iconn)%mass_balance(icomp,:) = &
        global_auxvars_bc(iconn)%mass_balance(icomp,:) + &
        global_auxvars_bc(iconn)%mass_balance_delta(icomp,:)* &
        towg_fmw_comp(icomp)*option%flow_dt
    enddo
  enddo
  do iconn = 1, patch%aux%TOWG%num_aux_ss
    do icomp = 1, option%nflowspec
      global_auxvars_ss(iconn)%mass_balance(icomp,:) = &
        global_auxvars_ss(iconn)%mass_balance(icomp,:) + &
        global_auxvars_ss(iconn)%mass_balance_delta(icomp,:)* &
        towg_fmw_comp(icomp)*option%flow_dt
    enddo
  enddo

end subroutine TOWGUpdateMassBalance

! ************************************************************************** !

subroutine TOWGImsTLBOComputeMassBalance(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/21/16
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec, &
                            realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  class(pm_towg_aux_type), pointer :: towg
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase, icomp
  PetscReal :: vol_phase,xmf,molar_density,mole_wt

  PetscBool :: is_black_oil
  PetscBool :: componentInPhase,is_oil_in_oil,is_gas_in_oil

  is_black_oil=PETSC_FALSE
  if(     ( towg_miscibility_model == TOWG_BLACK_OIL  ) &
     .or. ( towg_miscibility_model == TOWG_SOLVENT_TL ) ) is_black_oil=PETSC_TRUE

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  towg => patch%aux%TOWG
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    !not for both TOWG_IMS and TL phases and components coincides
    do iphase = 1, option%nphase
      ! volume_phase = saturation*porosity*volume
      vol_phase = &
        towg%auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase)* &
        towg%auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
        material_auxvars(ghosted_id)%volume
      ! mass = volume_phase*density 
      if( is_black_oil ) then
        do icomp=1,option%nphase
          componentInPhase=checkBlackOilCIP(iphase,icomp,is_oil_in_oil,is_gas_in_oil,option)
          if( componentInPhase ) then
            ! mole fraction of component in phase
            xmf=1.0
            if( is_oil_in_oil ) xmf=towg%auxvars(ZERO_INTEGER,ghosted_id)%bo%xo
            if( is_gas_in_oil ) xmf=towg%auxvars(ZERO_INTEGER,ghosted_id)%bo%xg
          ! mass = xmf*volume_phase*density
            molar_density=towg%auxvars(ZERO_INTEGER,ghosted_id)%den(iphase)
            mole_wt      =towg_fmw_comp(icomp)
            mass_balance(icomp,1) = mass_balance(icomp,1)+molar_density*mole_wt*vol_phase*xmf
          endif
        enddo
      else
      mass_balance(iphase,1) = mass_balance(iphase,1) + &
        towg%auxvars(ZERO_INTEGER,ghosted_id)%den(iphase)* &
        towg_fmw_comp(iphase) * vol_phase
      endif
    enddo
  enddo

end subroutine TOWGImsTLBOComputeMassBalance

! ************************************************************************** !

subroutine TOWGImsTLMaxChange(realization,max_change_ivar,max_change_isubvar,&
                              max_pressure_change,max_xmol_change, &
                              max_saturation_change,max_temperature_change)
  ! 
  ! Compute primary variable max changes
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/30/16
  ! 

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  !use Global_Aux_module
  !use General_Aux_module
  !use Variables_module, only : LIQUID_PRESSURE, LIQUID_MOLE_FRACTION, &
  !                             TEMPERATURE, GAS_PRESSURE, AIR_PRESSURE, &
  !                             GAS_SATURATION
  implicit none
  
  class(realization_subsurface_type), pointer :: realization
  PetscInt :: max_change_ivar(:)
  PetscInt :: max_change_isubvar(:)
  PetscReal :: max_pressure_change
  PetscReal :: max_xmol_change
  PetscReal :: max_saturation_change
  PetscReal :: max_temperature_change

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: max_change_local(4)
  PetscReal :: max_change_global(4)
  PetscReal :: max_change
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id

  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  max_change_global = 0.d0
  max_change_local = 0.d0
  
  ! 'TOWG_IMMISCIBLE','TODD_LONGOSTAFF'
  ! max change variables = [OIL_PRESSURE, OIL_SATURATION, &
  !                         GAS_SATURATION,TEMPERATURE]
  do i = 1, 4
    call RealizationGetVariable(realization,field%work, &
                                max_change_ivar(i),max_change_isubvar(i))
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_ptr2,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      ! have to weed out cells that changed state
      if (dabs(vec_ptr(j)) > 1.d-40 .and. dabs(vec_ptr2(j)) > 1.d-40) then
        max_change = max(max_change,dabs(vec_ptr(j)-vec_ptr2(j)))
      endif
    enddo
    max_change_local(i) = max_change
    call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%max_change_vecs(i),vec_ptr2, &
                            ierr);CHKERRQ(ierr)
    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)
  enddo
  call MPI_Allreduce(max_change_local,max_change_global,FOUR_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpo= ",1pe12.4, " dso= ",1pe12.4, &
      & " dsg= ",1pe12.4," dt= ",1pe12.4)') &
      max_change_global(1:4)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpo= ",1pe12.4, " dso= ",1pe12.4,&
      & " dsg= ",1pe12.4,/,15x,"  dt= ",1pe12.4)') &
      max_change_global(1:4)
  endif
 
  ! 'TOWG_IMMISCIBLE','TODD_LONGOSTAFF'
  ! max change variables = [OIL_PRESSURE, OIL_SATURATION, &
  !                         GAS_SATURATION,TEMPERATURE]
  max_pressure_change = max_change_global(1)
  max_xmol_change = 0.0d0
  max_saturation_change = maxval(max_change_global(2:3))
  max_temperature_change = max_change_global(4)
  
end subroutine TOWGImsTLMaxChange

! ************************************************************************** !

subroutine TOWGBOMaxChange(realization,max_change_ivar,max_change_isubvar,&
                           max_pressure_change,max_xmol_change, &
                           max_saturation_change,max_temperature_change)

!------------------------------------------------------------------------------
! Used in TOWG_BLACK_OIL, to compute primary variable max changes
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Dec 2017
!------------------------------------------------------------------------------

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Variables_module, only : OIL_PRESSURE

  implicit none

  class(realization_subsurface_type), pointer :: realization
  PetscInt :: max_change_ivar(:)
  PetscInt :: max_change_isubvar(:)
  PetscReal :: max_pressure_change
  PetscReal :: max_xmol_change
  PetscReal :: max_saturation_change
  PetscReal :: max_temperature_change

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: max_change_local(5)
  PetscReal :: max_change_global(5)
  PetscReal :: sum_values_local (2)
  PetscReal :: sum_values_global(2)
  PetscReal :: max_change,sum_pressure,normed_max_pb_change,average_pressure
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id

  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  max_change_global   = 0.d0
  max_change_local    = 0.d0
  sum_pressure        = 0.d0
  normed_max_pb_change= 0.d0
  average_pressure    = 0.d0

!------------------------------------------------------------------------------
! Look at the black oil variables:
! {OIL_PRESSURE,OIL_SATURATION,GAS_SATURATION,TEMPERATURE,BUBBLE_POINT}
!------------------------------------------------------------------------------

  do i = 1, 5

!--Get values------------------------------------------------------------------

    call RealizationGetVariable(realization,field%work, &
                                max_change_ivar(i),max_change_isubvar(i))
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_ptr2,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      ! have to weed out cells that changed state
      if (dabs(vec_ptr(j)) > 1.d-40 .and. dabs(vec_ptr2(j)) > 1.d-40) then
        max_change = max(max_change,dabs(vec_ptr(j)-vec_ptr2(j)))
      endif
    enddo
    max_change_local(i) = max_change

!--Case of pressure to provide norm for bubble point change--------------------

    if( max_change_ivar(i)==OIL_PRESSURE ) then
      do j = 1, grid%nlmax
        sum_values_local(1)=sum_pressure
        sum_values_local(2)=real(grid%nlmax)
      enddo
    endif

!--Restore values--------------------------------------------------------------

    call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%max_change_vecs(i),vec_ptr2, &
                            ierr);CHKERRQ(ierr)

    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)

  enddo

!------------------------------------------------------------------------------
! Do global reductions
!------------------------------------------------------------------------------

  call MPI_Allreduce(max_change_local,max_change_global,FIVE_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  call MPI_Allreduce(sum_values_local,sum_values_global,TWO_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

  if( sum_values_global(2)>0.0 ) then
    average_pressure=sum_values_global(1)/sum_values_global(2)
  endif

  if( average_pressure>0.0 ) then
   normed_max_pb_change=max_change_global(5)/average_pressure
  endif

!------------------------------------------------------------------------------
! Report the changes
!------------------------------------------------------------------------------

  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpo= ",1pe12.4, " dso= ",1pe12.4, &
      & " dsg= ",1pe12.4," dt= ",1pe12.4," dpb= ",1pe12.4)') &
      max_change_global(1:5)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpo= ",1pe12.4, " dso= ",1pe12.4,&
      & " dsg= ",1pe12.4,/,15x,"  dt= ",1pe12.4,"  dpb= ",1pe12.4)') &
      max_change_global(1:5)
  endif

!------------------------------------------------------------------------------
! Assemble the changes that control the simulation step
! Black oil model should be independent of mole weights,
! so do not set mole fraction max change
! Include normed maximum bubble point as an effective saturation change
!------------------------------------------------------------------------------

  max_pressure_change    = max_change_global(1)
  max_xmol_change        = 0.0d0
  max_saturation_change  = maxval(max_change_global(2:3))
  max_saturation_change  = max(max_saturation_change,normed_max_pb_change)
  max_temperature_change = max_change_global(4)

end subroutine TOWGBOMaxChange

subroutine TOWGTL4PMaxChange(realization,max_change_ivar,max_change_isubvar,&
                             max_pressure_change,max_xmol_change, &
                             max_saturation_change,max_temperature_change)

!------------------------------------------------------------------------------
! Used in TOWG_SOLVENT_TL, to compute primary variable max changes
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Variables_module, only : OIL_PRESSURE

  implicit none

  class(realization_subsurface_type), pointer :: realization
  PetscInt :: max_change_ivar(:)
  PetscInt :: max_change_isubvar(:)
  PetscReal :: max_pressure_change
  PetscReal :: max_xmol_change
  PetscReal :: max_saturation_change
  PetscReal :: max_temperature_change

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: max_change_local(6)
  PetscReal :: max_change_global(6)
  PetscReal :: sum_values_local (2)
  PetscReal :: sum_values_global(2)
  PetscReal :: max_change,sum_pressure,normed_max_pb_change,average_pressure
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id

  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  max_change_global   = 0.d0
  max_change_local    = 0.d0
  sum_pressure        = 0.d0
  normed_max_pb_change= 0.d0
  average_pressure    = 0.d0

!------------------------------------------------------------------------------
! Look at the TL4P oil variables (order defined by max_change_ivar)
! {OIL_PRESSURE,OIL_SATURATION,GAS_SATURATION,SOLVENT_SATURATION,TEMPERATURE,BUBBLE_POINT}
!------------------------------------------------------------------------------

  do i = 1, 6

!--Get values------------------------------------------------------------------

    call RealizationGetVariable(realization,field%work, &
                                max_change_ivar(i),max_change_isubvar(i))
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_ptr2,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      ! have to weed out cells that changed state
      if (dabs(vec_ptr(j)) > 1.d-40 .and. dabs(vec_ptr2(j)) > 1.d-40) then
        max_change = max(max_change,dabs(vec_ptr(j)-vec_ptr2(j)))
      endif
    enddo
    max_change_local(i) = max_change

!--Case of pressure to provide norm for bubble point change--------------------

    if( max_change_ivar(i)==OIL_PRESSURE ) then
      do j = 1, grid%nlmax
        sum_values_local(1)=sum_pressure
        sum_values_local(2)=real(grid%nlmax)
      enddo
    endif

!--Restore values--------------------------------------------------------------

    call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%max_change_vecs(i),vec_ptr2, &
                            ierr);CHKERRQ(ierr)

    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)

  enddo

!------------------------------------------------------------------------------
! Do global reductions
!------------------------------------------------------------------------------

  call MPI_Allreduce(max_change_local,max_change_global,SIX_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  call MPI_Allreduce(sum_values_local,sum_values_global,TWO_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

  if( sum_values_global(2)>0.0 ) then
    average_pressure=sum_values_global(1)/sum_values_global(2)
  endif

  if( average_pressure>0.0 ) then
   normed_max_pb_change=max_change_global(6)/average_pressure
  endif

!------------------------------------------------------------------------------
! Report the changes
!------------------------------------------------------------------------------

  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpo= ",1pe12.4, " dso= ",1pe12.4, &
      & " dsg= ",1pe12.4," dss= ",1pe12.4," dt= ",1pe12.4," dpb= ",1pe12.4)') &
      max_change_global(1:6)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpo= ",1pe12.4, " dso= ",1pe12.4,&
      & " dsg= ",1pe12.4," dss= ",1pe12.4,/,15x,"  dt= ",1pe12.4,"  dpb= ",1pe12.4)') &
      max_change_global(1:6)
  endif

!------------------------------------------------------------------------------
! Assemble the changes that control the simulation step
! TL4P model should be independent of mole weights,
! so do not set mole fraction max change
! Include normed maximum bubble point as an effective saturation change
!------------------------------------------------------------------------------

  max_pressure_change    = max_change_global(1)
  max_xmol_change        = 0.0d0
  max_saturation_change  = maxval(max_change_global(2:4))
  max_saturation_change  = max(max_saturation_change,normed_max_pb_change)
  max_temperature_change = max_change_global(5)

end subroutine TOWGTL4PMaxChange

! ************************************************************************** !

subroutine TOWGMapBCAuxVarsToGlobal(realization)
  !
  ! Map BC auxvars for TOWG problem coupled with reactive transport
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/06/16 
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Connection_module

  implicit none

  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  class(pm_towg_aux_type), pointer :: towg
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  

  PetscInt :: sum_connection, iconn
  
  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update
  
  towg => patch%aux%TOWG
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        towg%auxvars_bc(sum_connection)%sat
      global_auxvars_bc(sum_connection)%den_kg = &
        towg%auxvars_bc(sum_connection)%den_kg
      global_auxvars_bc(sum_connection)%temp = &
        towg%auxvars_bc(sum_connection)%temp
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine TOWGMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine TOWGSetPlotVariables(list)
  ! 
  ! Adds variables to be printed to list for TOWG module
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/29/16
  ! 
  
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(output_variable_list_type), pointer :: list
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_type), pointer :: output_variable

  if (associated(list%first)) then
    return
  endif
  
  name = 'Temperature'
  units = 'C'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               TEMPERATURE)

  name = 'Liquid Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               LIQUID_PRESSURE)

  name = 'Oil Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               OIL_PRESSURE)

  name = 'Gas Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               GAS_PRESSURE)

  name = 'Liquid Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)
  
  name = 'Oil Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               OIL_SATURATION)

  name = 'Gas Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               GAS_SATURATION)

  name = 'Liquid Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)
  
  name = 'Oil Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               OIL_DENSITY)

  name = 'Gas Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_DENSITY)
  
  name = 'Liquid Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_ENERGY)

  name = 'Oil Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               OIL_ENERGY)

  name = 'Gas Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_ENERGY)

! Bubble point output option for black oil and full solvent
  if ( towg_miscibility_model == TOWG_BLACK_OIL .or. &
       towg_miscibility_model == TOWG_SOLVENT_TL ) then
    name = 'Bubble Point'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               BUBBLE_POINT)
  end if

! Solvent model
  if ( towg_miscibility_model == TOWG_SOLVENT_TL ) then
    name = 'Solvent Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                 SOLVENT_SATURATION)
  end if

end subroutine TOWGSetPlotVariables

! ************************************************************************** !

subroutine TOWGUpdateFixedAccum(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/07/16 
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_class

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  class(pm_towg_aux_type), pointer :: towg
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end, natural_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
                          
  PetscErrorCode :: ierr

  PetscReal :: jdum(realization%option%nflowdof,realization%option%nflowdof)
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  towg => patch%aux%TOWG
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter
    
  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  !IF towg convergence criteria required: 
  !   initialize dynamic accumulation term for every p iteration step
  !if (towg_tough2_conv_criteria) then
  !  call VecGetArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  !endif
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    local_end = local_id*option%nflowdof
    local_start = local_end - option%nflowdof + 1
    ! TOWG_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = TOWG_UPDATE_FOR_FIXED_ACCUM
    call TOWGAuxVarCompute(xx_p(local_start:local_end), &
                           towg%auxvars(ZERO_INTEGER,ghosted_id), &
                           global_auxvars(ghosted_id), &
                           material_auxvars(ghosted_id), &
                           patch%characteristic_curves_array( &
                           patch%sat_func_id(ghosted_id))%ptr, &
                           natural_id, &
                           option)
    call TOWGAccumulation(towg%auxvars(ZERO_INTEGER,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          material_auxvars(ghosted_id), &
                          material_parameter%soil_heat_capacity(imat), &
                          option,accum_p(local_start:local_end), &
                          local_id == towg_debug_cell_id,jdum,PETSC_FALSE) 
  enddo
  
  !for tough2 convergence criteria
  !if (towg_tough2_conv_criteria) then
  !  accum_p2 = accum_p
  !endif
  
  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  
  !tough2 convergence criteria:
  ! initialize dynamic accumulation term for every p iteration step
  !if (towg_tough2_conv_criteria) then
  !  call VecRestoreArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  !endif
  
end subroutine TOWGUpdateFixedAccum

! ************************************************************************** !

subroutine TOWGImsTLBOAccumulation(auxvar,global_auxvar,material_auxvar, &
                                 soil_heat_capacity,option,Res,debug_cell,&
                                 j,analytical_derivatives)
  ! 
  ! Computes the non-fixed portion of the accumulation term for the residual
  ! Used for TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF,TOWG_BLACK_OIL,TOWG_SOLVENT_TL
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/07/16 
  ! 

  use Option_module
  use Material_module
  use Material_Aux_class
  use Derivatives_utilities_module 
  
  implicit none

  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof) 
  PetscBool :: debug_cell
  
  PetscInt :: iphase,icomp, energy_id
  
  PetscReal :: porosity
  PetscReal :: volume_over_dt,xmf

  PetscBool :: is_black_oil
  PetscBool :: componentInPhase,is_oil_in_oil,is_gas_in_oil

  PetscReal :: j(1:option%nflowdof,1:option%nflowdof) 
  PetscReal :: D_xmf(1:option%nflowdof) 
  PetscReal :: D_temp(1:option%nflowdof) 
  PetscBool :: analytical_derivatives
  PetscInt :: ndof 
  
  ndof = option%nflowdof

  if (analytical_derivatives) then
    j = 0.d0

    ! a bit silly but helps elegance later on:
    D_temp = 0.d0
    D_temp(towg_energy_dof) = 1.d0
  endif

! Set flag indicating a black oil run

  is_black_oil=PETSC_FALSE
  if(     ( towg_miscibility_model == TOWG_BLACK_OIL  ) &
     .or. ( towg_miscibility_model == TOWG_SOLVENT_TL ) ) is_black_oil=PETSC_TRUE
 
  energy_id = option%energy_id
  
  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use gen_auxvar%effective porosity here as it enables numerical 
  ! derivatives to be employed 
  porosity = auxvar%effective_porosity
  
  ! accumulation term units = kmol/s 
  ! not for TOWG IMS and TL (kmol phase) = (kmol comp)
  ! and nphase = nflowspec
  Res = 0.d0

  if( is_black_oil ) then
    do iphase = 1, option%nphase
      do icomp = 1, option%nphase ! In black oil, same number of components as phases
        ! Res[kmol comp/m^3 void] =  moleFraction
        !                           *sat[m^3 phase /m^3 void ]
        !                           *den[kmol phase/m^3 phase]
        componentInPhase=checkBlackOilCIP(iphase,icomp,is_oil_in_oil,is_gas_in_oil,option)

        if( componentInPhase ) then
          xmf=1.0d0
          if( is_oil_in_oil ) xmf=auxvar%bo%xo
          if( is_gas_in_oil ) xmf=auxvar%bo%xg

          if (analytical_derivatives) then
            D_xmf = 0.d0
            if( is_oil_in_oil ) D_xmf=auxvar%bo%D_xo
            if( is_gas_in_oil ) D_xmf=auxvar%bo%D_xg
          endif

          Res(icomp)=Res(icomp)+xmf*auxvar%sat(iphase)*auxvar%den(iphase)

          if (analytical_derivatives) then
            ! This is first instance of a XRule() type routine being used; see 
            ! derivatives_utilites.F90 for definition and description.
            J(icomp,:) = J(icomp,:) + ProdRule3(xmf,D_xmf,                                  &
                                            auxvar%sat(iphase),auxvar%D_sat(iphase,:),      &
                                            auxvar%den(iphase),auxvar%D_den(iphase,:),ndof )
          endif

        endif

      enddo
    enddo
  else
  do iphase = 1, option%nphase
    ! Res[kmol phase/m^3 void] = sat[m^3 phase/m^3 void] * 
    !                           den[kmol phase/m^3 phase] * 
    Res(iphase) = Res(iphase) + auxvar%sat(iphase) * &
                                auxvar%den(iphase) 

    if (analytical_derivatives) then
      J(iphase,:) = J(iphase,:) + ProdRule(auxvar%sat(iphase),auxvar%D_sat(iphase,:),      &
                                           auxvar%den(iphase),auxvar%D_den(iphase,:),ndof )
    endif

  enddo
  endif

  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] * 
  !---------  scale by porosity * volume / dt ----------------
  !                 vol[m^3 bulk] / dt[sec]

  ! do derivs first because it depends on unmodified value of Res.
  if (analytical_derivatives) then
    do icomp = 1,option%nflowspec
      J(icomp,:) = ProdRule(Res(icomp),J(icomp,:), &
                            porosity,auxvar%D_por,ndof) &
                 * volume_over_dt
    enddo
  endif

  Res(1:option%nflowspec) = Res(1:option%nflowspec) * &
                            porosity * volume_over_dt

  !---------  scale by porosity * volume / dt ----------------


  do iphase = 1, option%nphase
    ! Res[MJ/m^3 void] = sat[m^3 phase/m^3 void] *
    !                    den[kmol phase/m^3 phase] * U[MJ/kmol phase]
    Res(energy_id) = Res(energy_id) + auxvar%sat(iphase) * &
                                      auxvar%den(iphase) * &
                                      auxvar%U(iphase)
    if (analytical_derivatives) then
      J(energy_id,:) = J(energy_id,:) + ProdRule3(auxvar%sat(iphase),auxvar%D_sat(iphase,:),   &
                                                  auxvar%den(iphase),auxvar%D_den(iphase,:),   &
                                                  auxvar%U(iphase),auxvar%D_U(iphase,:),ndof    )
    endif

  enddo

  
  if (analytical_derivatives) then
    J(energy_id,:) = (  ProdRule(Res(energy_id),J(energy_id,:),     &
                                porosity,auxvar%D_por,ndof    )     &
                      + material_auxvar%soil_particle_density       &
                      * soil_heat_capacity                          &
                      * ProdRule(1.d0-porosity,-1.d0*auxvar%D_por,  &
                                 auxvar%temp,D_temp,ndof         )  &
                     ) * volume_over_dt
  endif
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] + 
  !                (1-por)[m^3 rock/m^3 bulk] * 
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * porosity + &
                    (1.d0 - porosity) * &
                    material_auxvar%soil_particle_density * &
                    soil_heat_capacity * auxvar%temp) * volume_over_dt

  
#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,7es24.15)') 'accum:', Res
  endif
#endif                    

end subroutine TOWGImsTLBOAccumulation

! ************************************************************************** !

subroutine TOWGImsTLBOFlux(auxvar_up,global_auxvar_up, &
                         material_auxvar_up, &
                         thermal_conductivity_up, &
                         auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         thermal_conductivity_dn, &
                         area, dist, towg_parameter, &
                         option,v_darcy,Res, &
                         debug_connection, &
                         jup,jdn,analytical_derivatives)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/08/16 
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  use Derivatives_utilities_module
  !use Fracture_module
  !use Klinkenberg_module
  
  implicit none
  
  class(auxvar_towg_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(towg_parameter_type) :: towg_parameter
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: Res(option%nflowdof)
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscInt :: energy_id
  PetscInt :: iphase,icomp
  
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure, delta_temp
  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: sat_liquid,sat_liquid_pos
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  !for debugging
  PetscReal :: adv_flux(option%nflowdof)
  PetscReal :: debug_flux(3), debug_dphi(3)
  
  PetscReal :: dummy_dperm_up, dummy_dperm_dn
  PetscReal :: temp_perm_up, temp_perm_dn
  PetscReal :: denup,dendn,xmf

  PetscBool :: is_black_oil,istl,isoil,isgas
  PetscBool :: componentInPhase,is_oil_in_oil,is_gas_in_oil

  PetscReal, dimension(1:option%nflowdof,1:option%nflowdof) :: jup, jdn
  PetscBool :: analytical_derivatives
  PetscInt :: ndof 
  

  PetscReal :: d_den_kg_ave_dden_up,d_den_kg_ave_dden_dn
  PetscReal :: d_den_ave_dden_up,d_den_ave_dden_dn
  PetscReal :: d_delta_temp_dt_up,d_delta_temp_dt_dn,dheat_flux_ddelta_temp

  PetscReal :: D_denup(option%nflowdof),D_dendn(option%nflowdof)
#ifndef GLOBALWORKERS
  PetscReal :: D_den_kg_ave_up(option%nflowdof),D_den_kg_ave_dn(option%nflowdof)
  PetscReal :: D_den_ave_up(option%nflowdof),D_den_ave_dn(option%nflowdof)
  PetscReal :: D_delta_presure_up(option%nflowdof),D_delta_presure_dn(option%nflowdof)
  PetscReal :: D_mobility_up(option%nflowdof),D_mobility_dn(option%nflowdof)
  PetscReal :: D_uH_up(option%nflowdof),D_uH_dn(option%nflowdof)
  PetscReal :: D_v_darcy_up(option%nflowdof),D_v_darcy_dn(option%nflowdof)
  PetscReal :: D_q_up(option%nflowdof),D_q_dn(option%nflowdof)
  PetscReal :: D_mole_flux_up(option%nflowdof),D_mole_flux_dn(option%nflowdof)
  PetscReal :: D_xmf_up(option%nflowdof),D_xmf_dn(option%nflowdof)
  PetscReal, dimension(1:option%nflowdof) :: D_sat_liquid_up,D_sat_liquid_dn,D_k_eff_up,D_k_eff_dn
  PetscReal, dimension(1:option%nflowdof) :: D_k_eff_ave_up,D_k_eff_ave_dn,D_delta_temp_up,D_delta_temp_dn
  PetscReal, dimension(1:option%nflowdof) :: D_worker1,D_worker2
#endif

  PetscReal :: worker1,worker2
  PetscInt :: i

#ifdef GLOBALWORKERS
  if (analytical_derivatives) then
    if (.NOT. CheckWorkersAllocated()) then
      ! something has gone horribly wrong here
      option%io_buffer = 'TOWGImsTLBOFlux: analytical derivatives mode is ON but &
                          intermediate workers are not allocated.'
     call PrintErrMsg(option)
    endif
  endif
#endif



  ndof = option%nflowdof
  if (analytical_derivatives) then
    jup = 0.d0; jdn = 0.d0
  endif

! Set flag indicating a black oil run

  is_black_oil=PETSC_FALSE
  if(     ( towg_miscibility_model == TOWG_BLACK_OIL  ) &
     .or. ( towg_miscibility_model == TOWG_SOLVENT_TL ) ) is_black_oil=PETSC_TRUE

  energy_id = option%energy_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalarSafe(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalarSafe(dist,perm_dn)
  
  ! Fracture permeability change only available for structured grid (Heeho)
  !if (associated(material_auxvar_up%fracture)) then
  !  call FracturePermEvaluate(material_auxvar_up,perm_up,temp_perm_up, &
  !                            dummy_dperm_up,dist)
  !  perm_up = temp_perm_up
  !endif
  !if (associated(material_auxvar_dn%fracture)) then
  !  call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
  !                            dummy_dperm_dn,dist)
  !  perm_dn = temp_perm_dn
  !endif
  
  !if (associated(klinkenberg)) then
  !  perm_ave_over_dist(1) = (perm_up * perm_dn) / &
  !                          (dist_up*perm_dn + dist_dn*perm_up)
  !  temp_perm_up = klinkenberg%Evaluate(perm_up, &
  !                                       auxvar_up%pres(option%gas_phase))
  !  temp_perm_dn = klinkenberg%Evaluate(perm_dn, &
  !                                       auxvar_dn%pres(option%gas_phase))
  !  perm_ave_over_dist(2) = (temp_perm_up * temp_perm_dn) / &
  !                          (dist_up*temp_perm_dn + dist_dn*temp_perm_up)
  !else
    if ((perm_up>0.0) .and. (perm_dn>0.0) ) then
      perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                              (dist_up*perm_dn + dist_dn*perm_up)
    else
      perm_ave_over_dist(:) = 0.0
    endif
  !endif
      
  Res = 0.d0
  
  v_darcy = 0.d0
#ifdef DEBUG_FLUXES  
  adv_flux = 0.d0
#endif
#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_flux = 0.d0
  debug_dphi = 0.d0
#endif

#ifdef CONVECTION

  do iphase = 1, option%nphase

    if (analytical_derivatives) then
      ! ensure all the worker variables are zero at the start of each loop
      d_den_kg_ave_dden_up=0.d0; d_den_kg_ave_dden_dn=0.d0
      d_den_ave_dden_up=0.d0;       d_den_ave_dden_dn=0.d0
      d_delta_temp_dt_up=0.d0;     d_delta_temp_dt_dn=0.d0; dheat_flux_ddelta_temp=0.d0
      D_den_kg_ave_up=0.d0;           D_den_kg_ave_dn=0.d0
      D_den_ave_up=0.d0;                 D_den_ave_dn=0.d0
      D_delta_presure_up=0.d0;     D_delta_presure_dn=0.d0
      D_mobility_up=0.d0;               D_mobility_dn=0.d0
      D_uH_up=0.d0;                           D_uH_dn=0.d0
      D_v_darcy_up=0.d0;                 D_v_darcy_dn=0.d0
      D_q_up=0.d0;                             D_q_dn=0.d0
      D_mole_flux_up=0.d0;             D_mole_flux_dn=0.d0
      D_xmf_up=0.d0;                         D_xmf_dn=0.d0
    endif
 
 !!! EXPERIMENTAL - turning off the cycle command remains a bad idea
 !!! update - this seems to be fixed now
 if (.NOT. analytical_derivatives) then
    if (auxvar_up%mobility(iphase) + &
        auxvar_dn%mobility(iphase) < eps) then
      cycle
    endif
 endif

    istl =(towg_miscibility_model == TOWG_TODD_LONGSTAFF)
    isoil=(option%phase_map(iphase) == OIL_PHASE)
    isgas=(option%phase_map(iphase) == GAS_PHASE)
    if( istl .and. (isoil.or.isgas) ) then
      if( isoil ) then
        denup=auxvar_up%tl%den_oil_eff_kg
        dendn=auxvar_dn%tl%den_oil_eff_kg
        if (analytical_derivatives) then
          D_denup = auxvar_up%tl%D_den_oil_eff_kg
          D_dendn = auxvar_dn%tl%D_den_oil_eff_kg
        endif
      else
        denup=auxvar_up%tl%den_gas_eff_kg
        dendn=auxvar_dn%tl%den_gas_eff_kg
        if (analytical_derivatives) then
          D_denup = auxvar_up%tl%D_den_gas_eff_kg
          D_dendn = auxvar_dn%tl%D_den_gas_eff_kg
        endif
      endif

      if (analytical_derivatives) then
        density_kg_ave = TOWGImsTLAverageDensity( auxvar_up%sat(iphase), &
                                                  auxvar_dn%sat(iphase), &
                                                  denup                , &
                                                  dendn                , &
                                                  d_den_kg_ave_dden_up , &
                                                  d_den_kg_ave_dden_dn    )
        !D_den_kg_ave_up = auxvar_up%D_den_kg(iphase,:)*d_den_kg_ave_dden_up
        !D_den_kg_ave_dn = auxvar_dn%D_den_kg(iphase,:)*d_den_kg_ave_dden_dn
        D_den_kg_ave_up = D_denup*d_den_kg_ave_dden_up
        D_den_kg_ave_dn = D_dendn*d_den_kg_ave_dden_dn
      else
        density_kg_ave = TOWGImsTLAverageDensity( auxvar_up%sat(iphase), &
                                                  auxvar_dn%sat(iphase), &
                                                  denup                , &
                                                  dendn )
      endif

    else
      if (analytical_derivatives) then
        density_kg_ave = TOWGImsTLAverageDensity( auxvar_up%sat(iphase)    , &
                                                  auxvar_dn%sat(iphase)    , &
                                                  auxvar_up%den_kg(iphase) , &
                                                  auxvar_dn%den_kg(iphase) , &
                                                  d_den_kg_ave_dden_up     , &
                                                  d_den_kg_ave_dden_dn         )
        ! d_den_kg_ave_dden_up comes from average density calc
        ! ( d (ave den) / d (den up) )
        D_den_kg_ave_up = auxvar_up%D_den_kg(iphase,:)*d_den_kg_ave_dden_up
        D_den_kg_ave_dn = auxvar_dn%D_den_kg(iphase,:)*d_den_kg_ave_dden_dn
      else
        density_kg_ave = TOWGImsTLAverageDensity( auxvar_up%sat(iphase)   , &
                                                  auxvar_dn%sat(iphase)   , &
                                                  auxvar_up%den_kg(iphase), &
                                                  auxvar_dn%den_kg(iphase) )
      endif

    endif

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = auxvar_up%pres(iphase) - &
                     auxvar_dn%pres(iphase) + &
                     gravity_term

    if (analytical_derivatives) then
      D_delta_presure_up = auxvar_up%D_pres(iphase,:) &
                         + D_den_kg_ave_up*dist_gravity
      D_delta_presure_dn = -auxvar_dn%D_pres(iphase,:) &
                           + D_den_kg_ave_dn*dist_gravity
    endif

#ifdef TOWG_DEBUG
      debug_dphi(iphase) = delta_pressure
#endif

    if (delta_pressure >= 0.D0) then
      mobility = auxvar_up%mobility(iphase)
      H_ave = auxvar_up%H(iphase)
      uH = H_ave

      if (analytical_derivatives) then
        D_mobility_up = auxvar_up%D_mobility(iphase,:)
        D_mobility_dn = 0.d0
        D_uH_up = auxvar_up%D_H(iphase,:)
        D_uH_dn = 0.d0
      endif

    else
      mobility = auxvar_dn%mobility(iphase)
      H_ave = auxvar_dn%H(iphase)
      uH = H_ave

      if (analytical_derivatives) then
        D_mobility_up = 0.d0
        D_mobility_dn = auxvar_dn%D_mobility(iphase,:)
        D_uH_up = 0.d0
        D_uH_dn = auxvar_dn%D_H(iphase,:)
      endif

    endif      

    !if (mobility > floweps) then  
    if (mobility > floweps .OR. analytical_derivatives) then ! not clear how much differnce this makes
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure

      if (analytical_derivatives) then
        D_v_darcy_up(:) = perm_ave_over_dist(iphase)                             &
                               * ProdRule(mobility,D_mobility_up,                &
                                          delta_pressure,D_delta_presure_up,ndof  )
        D_v_darcy_dn(:) = perm_ave_over_dist(iphase)                             &
                               * ProdRule(mobility,D_mobility_dn,                &
                                          delta_pressure,D_delta_presure_dn,ndof  )

        density_ave = TOWGImsTLAverageDensity(auxvar_up%sat(iphase), &
                                              auxvar_dn%sat(iphase), &
                                              auxvar_up%den(iphase), &
                                              auxvar_dn%den(iphase), &
                                              d_den_ave_dden_up    , &
                                              d_den_ave_dden_dn       )
        ! d_den_ave_dden_up comes from average density calc
        ! ( d (ave den) / d (den up) )
        D_den_ave_up = auxvar_up%D_den(iphase,:)*d_den_ave_dden_up
        D_den_ave_dn = auxvar_dn%D_den(iphase,:)*d_den_ave_dden_dn
      else
        density_ave = TOWGImsTLAverageDensity(auxvar_up%sat(iphase), &
                                              auxvar_dn%sat(iphase), &
                                              auxvar_up%den(iphase), &
                                              auxvar_dn%den(iphase) )
      endif

      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      if (analytical_derivatives) then
        D_q_up = D_v_darcy_up * area
        D_q_dn = D_v_darcy_dn * area
      endif

      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      mole_flux = q*density_ave

      if (analytical_derivatives) then
        D_mole_flux_up = ProdRule(q,D_q_up,                     &
                                  density_ave,D_den_ave_up,ndof  )
        D_mole_flux_dn = ProdRule(q,D_q_dn,                     &
                                  density_ave,D_den_ave_dn,ndof  )
      endif

      ! Res[kmol total/sec] = mole_flux[kmol phase/sec]

      if( is_black_oil ) then
! Loop over components in phase
        do icomp = 1, option%nphase
          componentInPhase=checkBlackOilCIP(iphase,icomp,is_oil_in_oil,is_gas_in_oil,option)
          if( componentInPhase ) then

            xmf=1.0d0
            if (delta_pressure >= 0.D0) then
              if( is_oil_in_oil ) xmf=auxvar_up%bo%xo
              if( is_gas_in_oil ) xmf=auxvar_up%bo%xg
            else
              if( is_oil_in_oil ) xmf=auxvar_dn%bo%xo
              if( is_gas_in_oil ) xmf=auxvar_dn%bo%xg
            endif

            if (analytical_derivatives) then
              D_xmf_up=0.0d0
              D_xmf_dn=0.0d0
              if (delta_pressure >= 0.D0) then
                if( is_oil_in_oil ) D_xmf_up=auxvar_up%bo%D_xo
                if( is_gas_in_oil ) D_xmf_up=auxvar_up%bo%D_xg
              else
                if( is_oil_in_oil ) D_xmf_dn=auxvar_dn%bo%D_xo
                if( is_gas_in_oil ) D_xmf_dn=auxvar_dn%bo%D_xg
              endif
            endif

            if (analytical_derivatives) then
              Jup(icomp,:) = Jup(icomp,:) + ProdRule(xmf,D_xmf_up,                &
                                                     mole_flux,D_mole_flux_up,ndof )
              Jdn(icomp,:) = Jdn(icomp,:) + ProdRule(xmf,D_xmf_dn,                &
                                                     mole_flux,D_mole_flux_dn,ndof )
            endif

            Res(icomp)=Res(icomp)+xmf*mole_flux

          endif
        enddo
      else
! One component in phase
        Res(iphase) = mole_flux
        icomp = iphase

        if (analytical_derivatives) then
          Jup(icomp,:) = Jup(icomp,:) + D_mole_flux_up(:)
          Jdn(icomp,:) = Jdn(icomp,:) + D_mole_flux_dn(:)
        endif

      endif

#ifdef DEBUG_FLUXES  
      adv_flux(iphase) = mole_flux
#endif

      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
      Res(energy_id) = Res(energy_id) + mole_flux * uH


      if (analytical_derivatives) then
        Jup(energy_id,:) = Jup(energy_id,:) + ProdRule(mole_flux,D_mole_flux_up, &
                                               uH,D_uH_up,ndof             )
        Jdn(energy_id,:) = Jdn(energy_id,:) + ProdRule(mole_flux,D_mole_flux_dn, &
                                               uH,D_uH_dn,ndof             )
      endif

#ifdef DEBUG_FLUXES  
      adv_flux(energy_id) = adv_flux(energy_id) + mole_flux * uH
#endif

#ifdef DEBUG_TOWG_FILEOUTPUT
      debug_dphi(iphase) = delta_pressure
      debug_flux(iphase) = mole_flux * uH
#endif
    endif                   
  enddo
! CONVECTION
#endif

#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'delta pressure :', debug_dphi(:)
    write(debug_unit,'(a,7es24.15)') 'adv flux (energy):', debug_flux(:)
  endif
  debug_flux = 0.d0
#endif                    


!!! Note energy contributions in tl4p mode analytical derivatives development was often questionable. 
!!! May be bug or potential for improvment here.
#ifdef CONDUCTION
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry) 
  ! Assuming that oil and water have same conductivity:
  ! s_l = S_water + S_oil 
  sat_liquid = auxvar_up%sat(option%liquid_phase) + &
               auxvar_up%sat(option%oil_phase)
  sat_liquid_pos=max(sat_liquid,0.0d0)
  k_eff_up = thermal_conductivity_up(1) + &
             sqrt(sat_liquid_pos) * &
             (thermal_conductivity_up(2) - thermal_conductivity_up(1))

  ! get corresponding upwind derivative contributions
  if (analytical_derivatives) then
    D_sat_liquid_up(:) = auxvar_up%D_sat(option%liquid_phase,:) + &
                         auxvar_up%D_sat(option%oil_phase,:)

    D_k_eff_up = 0.d0
    do i = 1,ndof
      !if (abs(D_sat_liquid_up(i)) < epsilon(sat_liquid) .OR.  sat_liquid_pos == 0.d0) then 
      if (sat_liquid_pos == 0.d0) then 
        D_k_eff_up(i) = 0.d0 
      else
        D_k_eff_up(i) = 0.5d0*D_sat_liquid_up(i)/sqrt(sat_liquid_pos)
      endif
    enddo
    D_k_eff_up = D_k_eff_up * &
                 (thermal_conductivity_up(2) - thermal_conductivity_up(1))

  endif


  sat_liquid = auxvar_dn%sat(option%liquid_phase) + &
               auxvar_dn%sat(option%oil_phase)
  sat_liquid_pos=max(sat_liquid,0.0d0)
  k_eff_dn = thermal_conductivity_dn(1) + &
             sqrt(sat_liquid_pos) * &
             (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))

  ! get corresponding downwind derivative contributions
  if (analytical_derivatives) then
    D_sat_liquid_dn = auxvar_dn%D_sat(option%liquid_phase,:) + &
                      auxvar_dn%D_sat(option%oil_phase,:)
    D_k_eff_dn = 0.d0
    do i = 1,ndof
      !if (abs(D_sat_liquid_dn(i)) < epsilon(sat_liquid) .OR.  sat_liquid_pos == 0.d0) then 
      if (sat_liquid_pos == 0.d0) then
        D_k_eff_dn(i) = 0.d0 
      else
        D_k_eff_dn(i) = 0.5d0*D_sat_liquid_dn(i)/sqrt(sat_liquid_pos)
      endif
    enddo
    D_k_eff_dn = D_k_eff_dn * &
                 (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
  endif

  if (k_eff_up > 0.d0 .or. k_eff_dn > 0.d0) then
    k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)


    if (analytical_derivatives) then
      ! complete the upwind derivatives
      D_k_eff_ave_up = 0.d0

      worker1 = k_eff_up*k_eff_dn
      ! derivative w.r.t. upwind variables:
      D_worker1 = k_eff_dn*D_k_eff_up 

      worker2 = k_eff_up*dist_dn+k_eff_dn*dist_up
      ! derivative w.r.t. upwind variables:
      D_worker2 = D_k_eff_up*dist_dn 

      D_k_eff_ave_up = DivRule(worker1,D_worker1, &
                               worker2,D_worker2,ndof)

      !complete the downwind derivatives
      D_k_eff_ave_dn = 0.d0

      worker1 = k_eff_up*k_eff_dn
      ! derivative w.r.t. downwind variables:
      D_worker1 = k_eff_up*D_k_eff_dn

      worker2 = k_eff_up*dist_dn+k_eff_dn*dist_up
      ! derivative w.r.t. downwind variables:
      D_worker2 = D_k_eff_dn*dist_up

      D_k_eff_ave_dn = DivRule(worker1,D_worker1, &
                               worker2,D_worker2,ndof)
    endif
  else
    k_eff_ave = 0.d0
    if (analytical_derivatives) then
      D_k_eff_ave_up = 0.d0
      D_k_eff_ave_dn = 0.d0
    endif
  endif
  ! units:
  ! k_eff = W/K-m = J/s/K-m
  ! delta_temp = K
  ! area = m^2
  ! heat_flux = k_eff * delta_temp * area = J/s
  delta_temp = auxvar_up%temp - auxvar_dn%temp
  heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! J/s -> MJ/s
  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux

  if (analytical_derivatives) then
    !!!! don't need these two anymore
    d_delta_temp_dt_up = 1.d0
    d_delta_temp_dt_dn = - 1.d0

    D_delta_temp_up = 0.d0
    D_delta_temp_up(towg_energy_dof) = 1.d0
    D_delta_temp_dn = 0.d0
    D_delta_temp_dn(towg_energy_dof) = -1.d0


    jup(energy_id,:) = jup(energy_id,:)              &
                                   + ProdRule(k_eff_ave,D_k_eff_ave_up,          &
                                              delta_temp,D_delta_temp_up,ndof )  &
                                   * area * 1.d-6

    jdn(energy_id,:) = jdn(energy_id,:)              &
                                   + ProdRule(k_eff_ave,D_k_eff_ave_dn,          &
                                              delta_temp,D_delta_temp_dn,ndof )  &
                                   * area * 1.d-6
  endif

! CONDUCTION
#endif
  
#ifdef DEBUG_FLUXES  
  if (debug_connection) then  
!    write(*,'(a,7es12.4)') 'in: ', adv_flux(:)*dist(1), diff_flux(:)*dist(1)
    write(*,'('' phase: liquid'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(1), auxvar_dn%pres(1)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(1), auxvar_dn%sat(1)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(1) * 1.d6
    write(*,'('' phase: oil'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(2), auxvar_dn%pres(2)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(2), auxvar_dn%sat(2)
    write(*,'(''  oil --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(2) * 1.d6
    write(*,'('' phase: gas'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(3), auxvar_dn%pres(3)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(3), auxvar_dn%sat(3)
    write(*,'(''  gas --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(3)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(3) * 1.d6
    write(*,'(''  energy --'')')
    write(*,'(''   advective heat flux:'',es12.4)') adv_flux(4) * 1.d6
    write(*,'(''   conductive heat flux:'',es12.4)') heat_flux * 1.d6
    write(*,'(''   total heat flux:'',es12.4)') (heat_flux + adv_flux(4))*1.d6

  endif
#endif

end subroutine TOWGImsTLBOFlux

! ************************************************************************** !

subroutine TOWGImsTLBOBCFlux(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                           auxvar_up,global_auxvar_up, &
                           auxvar_dn,global_auxvar_dn, &
                           material_auxvar_dn, &
                           thermal_conductivity_dn, &
                           area,dist,towg_parameter, &
                           option,v_darcy,Res,debug_connection,&
                           jdn,analytical_derivatives)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/22/16
  ! 
  use Option_module                              
  use Material_Aux_class
  !use Fracture_module
  !use Klinkenberg_module
  use Derivatives_utilities_module
  
  implicit none
  
  class(auxvar_towg_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: bc_auxvars(:) ! from aux_real_var array
  PetscReal :: v_darcy(option%nphase), area
  type(towg_parameter_type) :: towg_parameter
  PetscReal :: dist(-1:3)
  PetscReal :: Res(1:option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: bc_auxvar_mapping(TOWG_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)
  PetscBool :: debug_connection
  
  PetscInt :: energy_id
  PetscInt :: iphase,icomp
  PetscInt :: bc_type
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: H_ave, uH
  PetscReal :: perm_dn_adj(option%nphase)
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure, delta_temp !, delta_xmol
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: sat_dn, perm_dn, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, pres_ave
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux

  PetscReal :: adv_flux(option%nflowdof)
  PetscReal :: debug_flux(3), debug_dphi(3)

  PetscReal :: boundary_pressure
  PetscReal :: sat_liquid,sat_liquid_pos

  PetscReal :: dden_dn, dden_up

  PetscInt :: idof
  PetscBool :: neumann_bc_present,is_black_oil
  
  PetscReal :: temp_perm_dn
  PetscReal :: dummy_dperm_dn

  PetscReal :: denup,dendn,xmf

  PetscBool :: istl,isoil,isgas,is_oil_in_oil,is_gas_in_oil,componentInPhase

  PetscReal, dimension(1:option%nflowdof,1:option%nflowdof) :: jdn
  PetscBool :: analytical_derivatives
  PetscInt :: ndof 

  PetscReal, dimension(1:option%nflowdof) :: D_dendn
  PetscReal, dimension(1:option%nflowdof) :: D_den_kg_ave_dn
  PetscReal, dimension(1:option%nflowdof) :: D_den_ave_dn
  PetscReal, dimension(1:option%nflowdof) :: D_delta_presure_dn
  PetscReal, dimension(1:option%nflowdof) :: D_mobility_dn
  PetscReal, dimension(1:option%nflowdof) :: D_uH_dn
  PetscReal, dimension(1:option%nflowdof) :: D_v_darcy_dn
  PetscReal, dimension(1:option%nflowdof) :: D_q_dn
  PetscReal, dimension(1:option%nflowdof) :: D_mole_flux_dn
  PetscReal, dimension(1:option%nflowdof) :: D_xmf_dn
  PetscReal :: d_den_kg_ave_dden_dn
  PetscReal :: d_den_ave_dden_dn
  PetscReal :: d_delta_temp_dt_dn
  PetscReal :: dheat_flux_ddelta_temp
  PetscReal :: dummy

  PetscReal, dimension(1:option%nflowdof) :: D_sat_liquid_dn,D_k_eff_dn
  PetscReal, dimension(1:option%nflowdof) :: D_k_eff_ave_dn,D_delta_temp_dn
  PetscInt :: i

  ndof = option%nflowdof

  if (analytical_derivatives) then
    jdn = 0.d0
  endif

! Set flag indicating a black oil run

  is_black_oil=PETSC_FALSE
  if(     ( towg_miscibility_model == TOWG_BLACK_OIL  ) &
     .or. ( towg_miscibility_model == TOWG_SOLVENT_TL ) ) is_black_oil=PETSC_TRUE
  
  energy_id = option%energy_id

  Res = 0.d0
  v_darcy = 0.d0  

#ifdef DEBUG_FLUXES  
  adv_flux = 0.d0
#endif
#ifdef DEBUG_TOWG_FILEOUTPUT
  debug_flux = 0.d0
  debug_dphi = 0.d0
#endif

  neumann_bc_present = PETSC_FALSE
  
  call material_auxvar_dn%PermeabilityTensorToScalarSafe(dist,perm_dn)

  ! Fracture permeability change only available for structured grid (Heeho)
  !if (associated(material_auxvar_dn%fracture)) then
  !  call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
  !                            dummy_dperm_dn,dist)
  !  perm_dn = temp_perm_dn
  !endif  
  
  !if (associated(klinkenberg)) then
  !  perm_dn_adj(1) = perm_dn
  !                                        
  !  perm_dn_adj(2) = klinkenberg%Evaluate(perm_dn, &
  !                                        gen_auxvar_dn%pres(option%gas_phase))
  !else
    perm_dn_adj(:) = perm_dn
  !endif
  
#ifdef CONVECTION  
  do iphase = 1, option%nphase

    if (analytical_derivatives) then
      ! zero out all the workers 
      D_den_kg_ave_dn=0.d0
      D_den_ave_dn=0.d0
      D_delta_presure_dn=0.d0
      D_mobility_dn=0.d0
      D_uH_dn=0.d0
      D_v_darcy_dn=0.d0
      D_q_dn=0.d0
      D_mole_flux_dn=0.d0
      D_xmf_dn=0.d0
      d_den_kg_ave_dden_dn=0.d0
      d_den_ave_dden_dn=0.d0
      d_delta_temp_dt_dn=0.d0
      dheat_flux_ddelta_temp=0.d0
   endif
   ! zero this out
   density_ave = 0.d0
 
    bc_type = ibndtype(iphase)
    select case(bc_type)
      ! figure out the direction of flow
      case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,CONDUCTANCE_BC)

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

        if (bc_type == CONDUCTANCE_BC) then

! The values at TOWG_LIQ_CONDUCTANCE_INDEX etc are not actually set
         option%io_buffer = 'Boundary conductances are not available'
         call PrintErrMsg(option)

          select case(option%phase_map(iphase)) 
            case(LIQUID_PHASE)
              idof = bc_auxvar_mapping(TOWG_LIQ_CONDUCTANCE_INDEX)
            case(OIL_PHASE)
              idof = bc_auxvar_mapping(TOWG_OIL_CONDUCTANCE_INDEX)
            case(GAS_PHASE)
              idof = bc_auxvar_mapping(TOWG_GAS_CONDUCTANCE_INDEX)
            case(SOLVENT_PHASE)
              idof = bc_auxvar_mapping(TOWG_SOLVENT_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = bc_auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
          
        ! using residual saturation cannot be correct! - geh
        ! reusing sir_dn for bounary auxvar
#define BAD_MOVE1 ! this works
#ifndef BAD_MOVE1       
        if (auxvar_up%mobility(iphase) > eps .or. &
            auxvar_dn%mobility(iphase) > eps) then
#endif
          boundary_pressure = auxvar_up%pres(iphase)

          !PO: no free surfce boundaries considered           
          !if (iphase == LIQUID_PHASE .and. &
          !    global_auxvar_up%istate == GAS_STATE) then
          !  ! the idea here is to accommodate a free surface boundary
          !  ! face.  this will not work for an interior grid cell as
          !  ! there should be capillary pressure in force.
          !  boundary_pressure = gen_auxvar_up%pres(option%gas_phase)
          !endif

          istl =(towg_miscibility_model == TOWG_TODD_LONGSTAFF)
          isoil=(option%phase_map(iphase) == OIL_PHASE)
          isgas=(option%phase_map(iphase) == GAS_PHASE)
          if( istl .and. (isoil.or.isgas) ) then
            if( isoil ) then
              dendn=auxvar_dn%tl%den_oil_eff_kg
              if (analytical_derivatives) then
                D_dendn = auxvar_dn%tl%D_den_oil_eff_kg
              endif
            else
              dendn=auxvar_dn%tl%den_gas_eff_kg
              if (analytical_derivatives) then
                D_dendn = auxvar_dn%tl%D_den_gas_eff_kg
              endif
            endif
            if (analytical_derivatives) then
              density_kg_ave = TOWGImsTLAverageDensity( auxvar_up%sat(iphase), &
                                                        auxvar_dn%sat(iphase), &
                                                        denup                , &
                                                        dendn                , &
                                                        dummy                , &
                                                        d_den_kg_ave_dden_dn    )

              ! d_den_kg_ave_dden_up comes from average density calc
              ! ( d (ave den) / d (den up) )
              !D_den_kg_ave_up = auxvar_up%D_den_kg(iphase,:)*d_den_kg_ave_dden_up
              !D_den_kg_ave_dn = auxvar_dn%D_den_kg(iphase,:)*d_den_kg_ave_dden_dn
              D_den_kg_ave_dn = D_dendn*d_den_kg_ave_dden_dn
              else
                density_kg_ave = TOWGImsTLAverageDensity(auxvar_up%sat(iphase), &
                                                         auxvar_dn%sat(iphase), &
                                                         denup                , &
                                                         dendn )
            endif
                else

            if (analytical_derivatives) then
              density_kg_ave = TOWGImsTLAverageDensity( auxvar_up%sat(iphase), &
                                                        auxvar_dn%sat(iphase), &
                                                     auxvar_up%den_kg(iphase), &
                                                     auxvar_dn%den_kg(iphase), &
                                                        dummy , &
                                                        d_den_kg_ave_dden_dn    )
              ! d_den_kg_ave_dden_up comes from average density calc
              ! ( d (ave den) / d (den up) )
              D_den_kg_ave_dn = auxvar_dn%D_den_kg(iphase,:)*d_den_kg_ave_dden_dn
            else
              density_kg_ave = TOWGImsTLAverageDensity(auxvar_up%sat(iphase)   , &
                                                       auxvar_dn%sat(iphase)   , &
                                                       auxvar_up%den_kg(iphase), &
                                                       auxvar_dn%den_kg(iphase) )
            endif

          endif


          gravity_term = density_kg_ave * dist_gravity
          delta_pressure = boundary_pressure - &
                           auxvar_dn%pres(iphase) + &
                           gravity_term

          if (analytical_derivatives) then
            D_delta_presure_dn = -auxvar_dn%D_pres(iphase,:) &
                                 + D_den_kg_ave_dn*dist_gravity
          endif

#ifdef DEBUG_TOWG_FILEOUTPUT
          debug_dphi(iphase) = delta_pressure
#endif

          ! PO CONDUCTANCE_BC and HYDROSTATIC_SEEPAGE_BC to be implemented/tested
          if (bc_type == HYDROSTATIC_SEEPAGE_BC .or. &
              bc_type == CONDUCTANCE_BC) then
                ! flow in         ! boundary cell is <= pref
            if (delta_pressure > 0.d0 .and. &
                auxvar_up%pres(iphase) - &
                 option%reference_pressure < eps) then
              delta_pressure = 0.d0
              if (analytical_derivatives) then !!!! CHECK is this correct?
                D_delta_presure_up = 0.d0
                D_delta_presure_dn = 0.d0
              endif
            endif
          endif
          
          !upwinding mobility and enthalpy  
          if (delta_pressure >= 0.D0) then
            mobility = auxvar_up%mobility(iphase)
            uH = auxvar_up%H(iphase)
            if (analytical_derivatives) then
              D_mobility_dn = 0.d0
              D_uH_dn = 0.d0
            endif
          else
            mobility = auxvar_dn%mobility(iphase)
            uH = auxvar_dn%H(iphase)
            if (analytical_derivatives) then
              D_mobility_dn = auxvar_dn%D_mobility(iphase,:)
              D_uH_dn = auxvar_dn%D_H(iphase,:)
            endif
          endif      

          if (mobility > floweps .or. analytical_derivatives) then ! not clear how much difference this makes
          !if (mobility > floweps) then
            ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
            !                    dP[Pa]]
            v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure

            if (analytical_derivatives) then
              D_v_darcy_dn = perm_ave_over_dist                                        &
                                     * ProdRule(mobility,D_mobility_dn,                &
                                                delta_pressure,D_delta_presure_dn,ndof  )
            endif

            ! only need average density if velocity > 0.
            if (analytical_derivatives) then

              density_ave = TOWGImsTLAverageDensity(auxvar_up%sat(iphase), &
                                                    auxvar_dn%sat(iphase), &
                                                    auxvar_up%den(iphase), &
                                                    auxvar_dn%den(iphase), &
                                                    dummy    , &
                                                    d_den_ave_dden_dn       )
              ! d_den_ave_dden_up comes from average density calc
              ! ( d (ave den) / d (den up) )
              D_den_ave_dn = auxvar_dn%D_den(iphase,:)*d_den_ave_dden_dn
            else
              density_ave = TOWGImsTLAverageDensity(auxvar_up%sat(iphase), &
                                                    auxvar_dn%sat(iphase), &
                                                    auxvar_up%den(iphase), &
                                                    auxvar_dn%den(iphase) )
            endif
          endif
#ifndef BAD_MOVE1        
        endif ! sat > eps
#endif

      case(NEUMANN_BC)
        select case(option%phase_map(iphase))
          case(LIQUID_PHASE)
            idof = bc_auxvar_mapping(TOWG_LIQUID_FLUX_INDEX)
          case(OIL_PHASE)
            idof = bc_auxvar_mapping(TOWG_OIL_FLUX_INDEX)
          case(GAS_PHASE)
            idof = bc_auxvar_mapping(TOWG_GAS_FLUX_INDEX)
          case(SOLVENT_PHASE)
            idof = bc_auxvar_mapping(TOWG_SOLV_FLUX_INDEX)
        end select
        
        neumann_bc_present = PETSC_TRUE
        !xmol = 0.d0
        !xmol(iphase) = 1.d0
        if (dabs(bc_auxvars(idof)) > floweps) then
          v_darcy(iphase) = bc_auxvars(idof)
          !upwinding based on given BC flux sign
          if (analytical_derivatives) then
            D_v_darcy_dn(:) = 0.d0
          endif
          if (v_darcy(iphase) > 0.d0) then 
            density_ave = auxvar_up%den(iphase)
            uH = auxvar_up%H(iphase)
            if (analytical_derivatives) then
              D_den_ave_dn = 0.d0
              D_uH_dn = 0.d0
            endif
          else 
            density_ave = auxvar_dn%den(iphase)
            uH = auxvar_dn%H(iphase)
            if (analytical_derivatives) then
              D_den_ave_dn = auxvar_dn%D_den(iphase,:)
              D_uH_dn = auxvar_dn%D_H(iphase,:)
            endif
          endif 
        endif
      case default
        option%io_buffer = &
         'Boundary condition type not recognized in TOWGImsTLBOBCFlux phase loop'
        call PrintErrMsg(option)
    end select

    if (dabs(v_darcy(iphase)) > 0.d0 .OR. analytical_derivatives) then   !!!!! CHECK - need exception for aderivs here too?  -YES
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area
      if (analytical_derivatives) then
        D_q_dn = D_v_darcy_dn * area
      endif
      if (density_ave < 1.d-40) then
        option%io_buffer = 'Zero density in TOWGImsTLBOBCFlux()'
        call PrintErrMsgByRank(option)
      endif
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                              density_ave[kmol phase/m^3 phase]
      mole_flux = q*density_ave
      ! Res[kmol phase/sec] 
      if (analytical_derivatives) then
        D_mole_flux_dn = ProdRule(q,D_q_dn,                     &
                                  density_ave,D_den_ave_dn,ndof  )
      endif

      if( is_black_oil ) then

! Loop over components in this phase

        do icomp = 1, option%nphase
          componentInPhase=checkBlackOilCIP(iphase,icomp,is_oil_in_oil,is_gas_in_oil,option)
          if( componentInPhase ) then
            xmf=1.0d0
            if (v_darcy(iphase) > 0.d0) then
              if( is_oil_in_oil ) xmf=auxvar_up%bo%xo
              if( is_gas_in_oil ) xmf=auxvar_up%bo%xg
            else
              if( is_oil_in_oil ) xmf=auxvar_dn%bo%xo
              if( is_gas_in_oil ) xmf=auxvar_dn%bo%xg
            endif

            if (analytical_derivatives) then
              D_xmf_dn=0.0d0
              if (delta_pressure >= 0.D0) then
                ! we'd set upstream derivatives here but we won't use them
              else
                if( is_oil_in_oil ) D_xmf_dn=auxvar_dn%bo%D_xo
                if( is_gas_in_oil ) D_xmf_dn=auxvar_dn%bo%D_xg
              endif
            endif

            if (analytical_derivatives) then
              Jdn(icomp,:) = Jdn(icomp,:) + ProdRule(xmf,D_xmf_dn,                &
                                                     mole_flux,D_mole_flux_dn,ndof )
            endif

            Res(icomp)=Res(icomp)+xmf*mole_flux

          endif
        enddo

      else

! Just one component in this phase

      Res(iphase) = mole_flux 

      if (analytical_derivatives) then
        Jdn(iphase,:) = Jdn(iphase,:) + D_mole_flux_dn(:)
      endif

      endif

#ifdef DEBUG_FLUXES  
      adv_flux(iphase) = mole_flux 
#endif
      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
      Res(energy_id) = Res(energy_id) + mole_flux * uH ! H_ave

      if (analytical_derivatives) then
        Jdn(energy_id,:) = Jdn(energy_id,:) + ProdRule(mole_flux,D_mole_flux_dn, &
                                               uH,D_uH_dn,ndof             )
      endif

#ifdef DEBUG_FLUXES  
      adv_flux(energy_id) = adv_flux(energy_id) + mole_flux * uH
#endif
#ifdef DEBUG_TOWG_FILEOUTPUT
      debug_flux(iphase) = mole_flux * uH
#endif
    endif

  enddo
! CONVECTION
#endif

#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then  
    write(debug_unit,'(a,7es24.15)') 'bc delta pressure :', debug_dphi(:)
    write(debug_unit,'(a,7es24.15)') 'bc adv flux (energy):', debug_flux(:)
  endif
  debug_flux = 0.d0
#endif                    

#ifdef CONDUCTION
  ! add heat conduction flux
  heat_flux = 0.d0
  select case (ibndtype(towg_energy_eq_idx))
    case (DIRICHLET_BC)
      ! based on Somerton et al., 1974:
      ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
      ! Assuming that oil and water have same conductivity:
      ! s_l = S_water + S_oil 
      sat_liquid = auxvar_dn%sat(option%liquid_phase) + &
                   auxvar_dn%sat(option%oil_phase)
      sat_liquid_pos=max(sat_liquid,0.0d0)

      k_eff_dn = thermal_conductivity_dn(1) + &
                 sqrt(sat_liquid_pos) * &
                 (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))

  ! corresponding downwind derivatives:
  if (analytical_derivatives) then
    D_sat_liquid_dn = auxvar_dn%D_sat(option%liquid_phase,:) + &
                      auxvar_dn%D_sat(option%oil_phase,:)
    D_k_eff_dn = 0.d0
    do i = 1,ndof
      !if (abs(D_sat_liquid_dn(i)) < epsilon(sat_liquid) .OR.  sat_liquid_pos == 0.d0) then  ! probably don't need this
      if (sat_liquid_pos == 0.d0) then
        D_k_eff_dn(i) = 0.d0 
      else
        D_k_eff_dn(i) = 0.5d0*D_sat_liquid_dn(i)/sqrt(sat_liquid_pos)
      endif
    enddo
    D_k_eff_dn = D_k_eff_dn * &
                 (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
  endif

      ! units:
      ! k_eff = W/K/m/m = J/s/K/m/m
      ! delta_temp = K
      ! area = m^2
      ! heat_flux = J/s
      k_eff_ave = k_eff_dn / dist(0)
      delta_temp = auxvar_up%temp - auxvar_dn%temp
      heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! convert W -> MW

  ! corresponding downwind derivatives:
  if (analytical_derivatives) then

      D_k_eff_ave_dn = D_k_eff_dn / dist(0)
      D_delta_temp_dn = 0.d0
      D_delta_temp_dn(towg_energy_dof) = -1.d0

      jdn(energy_id,:) = jdn(energy_id,:)              &
                                     + ProdRule(k_eff_ave,D_k_eff_ave_dn,          &
                                                delta_temp,D_delta_temp_dn,ndof )  &
                                     * area * 1.d-6

  endif

    case(NEUMANN_BC)
                  ! flux prescribed as MW/m^2
      heat_flux = bc_auxvars(bc_auxvar_mapping(TOWG_ENERGY_FLUX_INDEX)) * area

    case default
      option%io_buffer = 'Boundary condition type not recognized in ' // &
        'TOWGImsTLBOBCFlux heat conduction loop.'
      call PrintErrMsg(option)
  end select
  Res(energy_id) = Res(energy_id) + heat_flux ! MW
! CONDUCTION
#endif


#ifdef DEBUG_FLUXES  
  if (debug_connection) then  
    write(*,'('' bc phase: liquid'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(1), auxvar_dn%pres(1)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(1), auxvar_dn%sat(1)
    write(*,'(''  water --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(1)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(1) * 1.d6
    write(*,'('' bc phase: oil'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(2), auxvar_dn%pres(2)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(2), auxvar_dn%sat(2)
    write(*,'(''  oil --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(2)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(2) * 1.d6
    write(*,'('' bc phase: gas'')')
    write(*,'(''  pressure   :'',2es12.4)') auxvar_up%pres(3), auxvar_dn%pres(3)
    write(*,'(''  saturation :'',2es12.4)') auxvar_up%sat(3), auxvar_dn%sat(3)
    write(*,'(''  gas --'')')
    write(*,'(''   darcy flux:'',es12.4)') adv_flux(3)
    write(*,'(''   heat adv. flux:'',es12.4)') debug_flux(3) * 1.d6
    write(*,'(''  bc energy --'')')
    write(*,'(''   advective heat flux:'',es12.4)') adv_flux(4) * 1.d6
    write(*,'(''   conductive heat flux:'',es12.4)') heat_flux * 1.d6
    write(*,'(''   total heat flux:'',es12.4)') (heat_flux + adv_flux(4))*1.d6

  endif
#endif
  
end subroutine TOWGImsTLBOBCFlux

! ************************************************************************** !

subroutine TOWGImsTLSrcSink(option,src_sink_condition, auxvar, &
                            global_auxvar,ss_flow_vol_flux,scale,Res,&
                            j,analytical_derivatives)
  ! 
  ! Computes the source/sink terms for the residual 
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/27/16
  ! 

  use Option_module
  use Condition_module  

  use EOS_Water_module
  use EOS_Oil_module
  use EOS_Gas_module
  use EOS_Slv_module
  use Derivatives_utilities_module

  implicit none

  type(option_type) :: option
  type(flow_towg_condition_type), pointer :: src_sink_condition
  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar !keep global_auxvar for salinity
  PetscReal :: ss_flow_vol_flux(option%nphase)
  PetscReal :: scale  
  PetscReal :: Res(option%nflowdof)

  ! local parameter
  PetscInt, parameter :: SRC_TEMPERATURE = 1
  PetscInt, parameter :: SRC_ENTHALPY = 2 
  ! local variables
  PetscReal, pointer :: qsrc(:)
  PetscInt :: flow_src_sink_type    
  PetscReal :: qsrc_mol
  PetscReal :: den, den_kg, enthalpy, internal_energy_dummy, temperature
  PetscReal :: cell_pressure
  PetscInt :: iphase
  PetscInt :: energy_var
  PetscErrorCode :: ierr

  PetscReal, dimension(1:option%nflowdof,1:option%nflowdof) :: j
  PetscBool :: analytical_derivatives
  PetscReal, dimension(1:option%nflowdof) :: D_den,D_qsrc_mol,D_enthalpy,D_cpres

  PetscReal :: dx_dcpres, dx_dtcell, dT_dTcell
  PetscInt :: ndof, cp_loc, idex
  PetscInt :: dof_op,dof_osat,dof_gsat,dof_temp
  PetscReal :: dum1,dum2



#if 0
  if (analytical_derivatives) then
    option%io_buffer = 'TOWGImsTLSrcSink: analytical derivatives are not yet available.'
    call PrintErrMsg(option)
  endif
#endif

  dof_op = TOWG_OIL_PRESSURE_DOF
  dof_osat = TOWG_OIL_SATURATION_DOF
  dof_gsat = TOWG_GAS_SATURATION_3PH_DOF
  dof_temp = towg_energy_dof

  if (analytical_derivatives) then
    j = 0.d0
    ndof = option%nflowdof
    ! zero out the intermediates:
    D_den = 0.d0
    D_qsrc_mol = 0.d0; D_enthalpy = 0.d0
  endif

  ! this can be removed if extending to pressure condition
  if (.not.associated(src_sink_condition%rate) ) then
    option%io_buffer = 'TOWGImsTLSrcSink fow condition rate not defined ' // &
    'rate is needed for a valid src/sink term'
    call PrintErrMsg(option)
  end if

  qsrc => src_sink_condition%rate%dataset%rarray

  energy_var = 0
  if ( associated(src_sink_condition%temperature) ) then
    energy_var = SRC_TEMPERATURE 
  else if ( associated(src_sink_condition%enthalpy) ) then
    energy_var = SRC_ENTHALPY
  end if

  flow_src_sink_type = src_sink_condition%rate%itype

  ! checks that qsrc(liquid_phase), qsrc(oil_phase), qsrc(gas_phase) 
  ! do not have different signs
  ! if ( (qsrc(option%liquid_phase)>0.0d0 .and. qsrc(option%oil_phase)<0.d0).or.&
  !     (qsrc(option%liquid_phase)<0.0d0 .and. qsrc(option%oil_phase)>0.d0)  & 
  !   ) then
  !   option%io_buffer = "TOilImsSrcSink error: " // &
  !     "src(wat) and src(oil) with opposite sign"
  !   call PrintErrMsg(option)
  ! end if

  ! if not given, approximates BHP with pressure of perforated cell
  if ( associated(src_sink_condition%bhp_pressure) ) then
    cell_pressure = src_sink_condition%bhp_pressure%dataset%rarray(1)
      if (analytical_derivatives) then
        D_cpres = 0.d0
      endif
  else
    cell_pressure = &
        maxval(auxvar%pres(option%liquid_phase:option%gas_phase))
        if (analytical_derivatives) then
          cp_loc = option%liquid_phase
          do idex = option%liquid_phase,option%gas_phase
            if (auxvar%pres(idex) > auxvar%pres(cp_loc)) then
            cp_loc = idex
            endif
          end do
          D_cpres = auxvar%D_pres(cp_loc,:)
        endif
  end if

  ! if enthalpy is used to define enthalpy or energy rate is used  
  ! approximate bottom hole temperature (BHT) with local temp
  if ( associated(src_sink_condition%temperature) ) then
    temperature = src_sink_condition%temperature%dataset%rarray(1)
    if (analytical_derivatives) then
      dT_dTcell = 0.d0
    endif
  else   
    temperature = auxvar%temp
    if (analytical_derivatives) then
      dT_dTcell = 1.d0
    endif
  end if

  Res = 0.d0
  do iphase = 1, option%nphase
    qsrc_mol = 0.d0
    if ( qsrc(iphase) > 0.d0) then 
      select case(option%phase_map(iphase))
        case(LIQUID_PHASE)
          if (analytical_derivatives) then
            call EOSWaterDensity(temperature,cell_pressure, &
                                  den_kg,den, &
                                  dx_dcpres, &
                                  dx_dtcell, ierr, auxvar%table_idx)
            D_den = D_cpres*dx_dcpres
            D_den(dof_temp) = D_den(dof_temp) + dx_dtcell*dT_dTcell
          else
            call EOSWaterDensity(temperature,cell_pressure,den_kg,den,ierr, &
                                 auxvar%table_idx)
          endif
        case(OIL_PHASE)
          if (analytical_derivatives) then
             call EOSOilDensity(temperature,cell_pressure,den,dx_dtcell,dx_dcpres,&
                     ierr,auxvar%table_idx)
            D_den = D_cpres*dx_dcpres
            D_den(dof_temp) = D_den(dof_temp) + dx_dtcell*dT_dTcell
          else
            call EOSOilDensity(temperature,cell_pressure,den,ierr, &
                               auxvar%table_idx)
          endif
        case(GAS_PHASE)
          if (analytical_derivatives) then
            call EOSGasDensity(temperature,cell_pressure,den,dx_dtcell,dx_dcpres, &
                                ierr,auxvar%table_idx)
            D_den = D_cpres*dx_dcpres
            D_den(dof_temp) = D_den(dof_temp) + dx_dtcell*dT_dTcell
          else
            call EOSGasDensity(temperature,cell_pressure,den,ierr, &
                               auxvar%table_idx)
          endif
        case(SOLVENT_PHASE)
          if (analytical_derivatives) then
           call EOSSlvDensity(temperature,cell_pressure,den,dx_dcpres,dx_dtcell, &
                               ierr,auxvar%table_idx)
            D_den = D_cpres*dx_dcpres
            D_den(dof_temp) = D_den(dof_temp) + dx_dtcell*dT_dTcell
          else
            call EOSSlvDensity(temperature,cell_pressure,den,ierr, &
                               auxvar%table_idx)
          endif
      end select 
    else
      if (analytical_derivatives) then
        D_den = auxvar%D_den(iphase,:)
      endif
      den = auxvar%den(iphase)
    end if

    select case(flow_src_sink_type)
      ! injection and production 
      case(MASS_RATE_SS)
        qsrc_mol = qsrc(iphase)/towg_fmw_comp(iphase) ! kg/sec -> kmol/sec
        if (analytical_derivatives) then
          !!! qsrc and fmw both 0
          D_qsrc_mol = 0.d0
        endif
      case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
        qsrc_mol = qsrc(iphase)/towg_fmw_comp(iphase)*scale 
        if (analytical_derivatives) then
          !!! qsrc and fmw both 0
          D_qsrc_mol = 0.d0
        endif
      case(VOLUMETRIC_RATE_SS)  ! assume local density for now 
                  ! qsrc(iphase) = m^3/sec  
        qsrc_mol = qsrc(iphase)*den ! den = kmol/m^3 
        if (analytical_derivatives) then
          D_qsrc_mol = qsrc(iphase)*D_den
        endif
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec             ! den = kmol/m^3
        qsrc_mol = qsrc(iphase)* den * scale
        if (analytical_derivatives) then
          D_qsrc_mol = qsrc(iphase)*D_den*scale
        endif
    end select
    ss_flow_vol_flux(iphase) = qsrc_mol/ den
    !!! note we don't provide an ss_flow_vol_flux deriv
    if (analytical_derivatives) then
      J(iphase,:) = D_qsrc_mol(:)
    endif
    Res(iphase) = qsrc_mol
  enddo

  ! when using scaled src/sinks, the rates (mass or vol) scaling 
  ! at this point the scale factor is already included in Res(iphase)

  ! Res(option%energy_id), energy units: MJ/sec

  !if ( associated(src_sink_condition%temperature) .or. &
  !    associated(src_sink_condition%enthalpy) &
  !   ) then
  !if the energy rate is not given, use either temperature or enthalpy
  if (analytical_derivatives) then 
    D_enthalpy= 0.d0
  endif
  if ( dabs(qsrc(option%energy_id)) < 1.0d-40 ) then
    ! water injection 
    if (qsrc(option%liquid_phase) > 0.d0) then !implies qsrc(option%oil_phase)>=0
      if ( energy_var == SRC_ENTHALPY ) then
        !input as J/kg
        enthalpy = src_sink_condition%enthalpy% &
                       dataset%rarray(option%liquid_phase)
                     ! J/kg * kg/kmol = J/kmol  
        enthalpy = enthalpy * towg_fmw_comp(option%liquid_phase)
        !!! no analytical derivatives here
      else !note: temp can either be input or taken as the one of perf. block
      !else if ( energy_var == SRC_TEMPERATURE ) then
        if (analytical_derivatives) then
          call EOSWaterEnthalpy(temperature, cell_pressure,enthalpy,dx_dcpres,dx_dtcell,ierr)
          D_enthalpy = D_cpres*dx_dcpres
          D_enthalpy(dof_temp) = D_enthalpy(dof_temp) + dx_dtcell*dT_dTcell
        else
          call EOSWaterEnthalpy(temperature, cell_pressure,enthalpy,ierr)
          ! enthalpy = [J/kmol]
        endif
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      if (analytical_derivatives) then
        D_enthalpy = D_enthalpy * 1.d-6 ! J/kmol -> whatever units
        J(option%energy_id,:) = J(option%energy_id,:) &
                              + ProdRule(Res(option%liquid_phase),J(option%liquid_phase,:),&
                                             enthalpy,D_enthalpy,ndof)
      endif
      ! enthalpy units: MJ/kmol
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%liquid_phase) * enthalpy
    end if
    ! oil injection 
    if (analytical_derivatives) then 
      D_enthalpy= 0.d0
    endif
    if (qsrc(option%oil_phase) > 0.d0) then !implies qsrc(option%liquid_phase)>=0
      if ( energy_var == SRC_ENTHALPY ) then
        enthalpy = src_sink_condition%enthalpy% &
                     dataset%rarray(option%oil_phase)
                      !J/kg * kg/kmol = J/kmol  
        enthalpy = enthalpy * towg_fmw_comp(option%oil_phase)
        !!! no analytical derivatives here
      else !note: temp can either be input or taken as the one of perf. block
      !if ( energy_var == SRC_TEMPERATURE ) then
        if (analytical_derivatives) then
          call EOSOilEnthalpy(temperature, cell_pressure,enthalpy,dx_dcpres,dx_dtcell,ierr)
          D_enthalpy = D_cpres*dx_dcpres
          D_enthalpy(dof_temp) = D_enthalpy(dof_temp) + dx_dtcell*dT_dTcell
        else
          call EOSOilEnthalpy(temperature,cell_pressure,enthalpy,ierr)
        endif
        ! enthalpy = [J/kmol]
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      if (analytical_derivatives) then
        D_enthalpy = D_enthalpy * 1.d-6 ! J/kmol -> whatever units
        J(option%energy_id,:) = J(option%energy_id,:) &
                              + ProdRule(Res(option%oil_phase),J(option%oil_phase,:),&
                                             enthalpy,D_enthalpy,ndof)
      endif
      ! enthalpy units: MJ/kmol
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%oil_phase) * enthalpy
    end if
    ! gas injection 
    if (qsrc(option%gas_phase) > 0.d0) then
      if ( energy_var == SRC_ENTHALPY ) then
        enthalpy = src_sink_condition%enthalpy% &
                     dataset%rarray(option%gas_phase)
                      !J/kg * kg/kmol = J/kmol
        enthalpy = enthalpy * towg_fmw_comp(option%gas_phase)
      else !note: temp can either be input or taken as the one of perf. block
      !if ( energy_var == SRC_TEMPERATURE ) then
        if (analytical_derivatives) then
          !subroutine EOSGasEnergyDerive(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
          call EOSGasEnergy(temperature, cell_pressure,enthalpy,dx_dtcell,dx_dcpres,&
                            internal_energy_dummy,dum1,dum2,ierr)
          D_enthalpy = D_cpres*dx_dcpres
          D_enthalpy(dof_temp) = D_enthalpy(dof_temp) + dx_dtcell*dT_dTcell
        else
          call EOSGasEnergy(temperature,cell_pressure,enthalpy, &
                              internal_energy_dummy,ierr)
        endif
        ! enthalpy = [J/kmol]
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      if (analytical_derivatives) then
        D_enthalpy = D_enthalpy * 1.d-6 ! J/kmol -> whatever units
        J(option%energy_id,:) = J(option%energy_id,:) &
                              + ProdRule(Res(option%gas_phase),J(option%gas_phase,:),&
                                             enthalpy,D_enthalpy,ndof)
      endif
      ! enthalpy units: MJ/kmol
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%gas_phase) * enthalpy
    end if
    ! water energy extraction due to water production
    if (qsrc(option%liquid_phase) < 0.d0) then !implies qsrc(option%oil_phase)<=0
      ! auxvar enthalpy units: MJ/kmol
      if (analytical_derivatives) then
        J(option%energy_id,:) = J(option%energy_id,:)                                               &
                              + ProdRule(Res(option%liquid_phase),J(option%liquid_phase,:),               &
                                         auxvar%H(option%liquid_phase),auxvar%D_H(option%liquid_phase,:), &
                                         ndof)
      endif
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%liquid_phase) * &
                              auxvar%H(option%liquid_phase)
    end if
    !oil energy extraction due to oil production
    if (qsrc(option%oil_phase) < 0.d0) then !implies qsrc(option%liquid_phase)<=0
      ! auxvar enthalpy units: MJ/kmol
      if (analytical_derivatives) then
        J(option%energy_id,:) = J(option%energy_id,:)                                               &
                              + ProdRule(Res(option%oil_phase),J(option%oil_phase,:),               &
                                         auxvar%H(option%oil_phase),auxvar%D_H(option%oil_phase,:), &
                                         ndof)
      endif
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%oil_phase) * &
                              auxvar%H(option%oil_phase)
    end if
    if (qsrc(option%gas_phase) < 0.d0) then !implies qsrc(option%liquid_phase)<=0
      ! auxvar enthalpy units: MJ/kmol
      if (analytical_derivatives) then
        J(option%energy_id,:) = J(option%energy_id,:)                                               &
                              + ProdRule(Res(option%gas_phase),J(option%gas_phase,:),               &
                                         auxvar%H(option%gas_phase),auxvar%D_H(option%gas_phase,:), &
                                         ndof)
      endif
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%gas_phase) * &
                              auxvar%H(option%gas_phase)
    end if
    if( towg_miscibility_model == TOWG_SOLVENT_TL ) then
      if (qsrc(option%solvent_phase) < 0.d0) then !implies qsrc(option%solvent_phase)<=0
        ! auxvar enthalpy units: MJ/kmol
      if (analytical_derivatives) then
        J(option%energy_id,:) = J(option%energy_id,:)                                               &
                              + ProdRule(Res(option%solvent_phase),J(option%solvent_phase,:),               &
                                         auxvar%H(option%solvent_phase),auxvar%D_H(option%solvent_phase,:), &
                                         ndof)
      endif
        Res(option%energy_id) = Res(option%energy_id) + &
                                Res(option%solvent_phase) * &
                                auxvar%H(option%solvent_phase)
      end if
    end if
  else !if the energy rate is given, it overwrites both temp and enthalpy
    ! if energy rate is given, loaded in qsrc(4) in MJ/sec 
    Res(option%energy_id) = qsrc(option%energy_id)* scale ! MJ/s
  end if

  nullify(qsrc)      
  
end subroutine TOWGImsTLSrcSink

! ************************************************************************** !

subroutine TOWGBOSrcSink(option,src_sink_condition, auxvar, &
                         global_auxvar,ss_flow_vol_flux,scale,Res,&
                         j, analytical_derivatives)
  ! 
  ! Computes the source/sink terms for the residual in the black oil case
  ! 
  ! Author: Dave Ponting
  ! Date: Dec 2017
  ! 

  use Option_module
  use Condition_module

  use EOS_Water_module
  use EOS_Oil_module
  use EOS_Gas_module
  use EOS_Slv_module
  use Derivatives_utilities_module
  
  implicit none

  type(option_type) :: option
  type(flow_towg_condition_type), pointer :: src_sink_condition
  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar !keep global_auxvar for salinity
  PetscReal :: ss_flow_vol_flux(option%nphase)
  PetscReal :: scale
  PetscReal :: Res(option%nflowdof)

  ! local parameter
  PetscInt, parameter :: SRC_TEMPERATURE = 1
  PetscInt, parameter :: SRC_ENTHALPY = 2
  ! local variables
  PetscReal, pointer :: qsrc(:)
  PetscInt :: flow_src_sink_type
  PetscReal :: qsrc_mol,mole_wt,ref_pressure
  PetscReal :: den, den_kg, enthalpy, internal_energy_dummy, temperature &
              ,xmf,xmfo,xmfg,po,pb,cr,crusp
  PetscReal :: cell_pressure
  PetscInt :: iphase,icomp
  PetscInt :: energy_var
  PetscBool :: componentInPhase,is_oil_in_oil,is_gas_in_oil
  PetscErrorCode :: ierr

  PetscReal, dimension(1:option%nflowdof,1:option%nflowdof) :: j
  PetscBool :: analytical_derivatives

  PetscReal :: dp_dpo,dT_dTcell
  PetscReal, dimension(1:option%nflowdof) :: D_xmfo,D_xmfg,D_mole_wt,D_den
  PetscReal, dimension(1:option%nflowdof) :: D_qsrc_mol,D_enthalpy,D_xmf
  PetscReal, dimension(1:option%nflowdof) :: D_cpres
  PetscReal, dimension(1:option%nflowdof) :: D_crusp, D_cor

  PetscInt :: ndof, cp_loc, idex
  PetscInt :: dof_op,dof_osat,dof_gsat,dof_temp
  PetscReal :: cor,one_p_crusp,dcor_dpo,dcor_dpb,dcor_dt
  PetscReal :: dcr_dt,dcr_dpb,dum1,dum2,mxpcand
  PetscReal :: dcrusp_dpo,dcrusp_dpb,dcrusp_dt
  PetscReal :: dx_dcpres, dx_dtcell
  PetscBool :: isSat

  PetscReal :: dden_dt, dummy


  dof_op = TOWG_OIL_PRESSURE_DOF
  dof_osat = TOWG_OIL_SATURATION_DOF
  dof_gsat = TOWG_GAS_SATURATION_3PH_DOF
  dof_temp = towg_energy_dof


  !ref_pressure=option%reference_pressure
  ref_pressure = 1.0D5 !used as Pb for dead oil injection

  if (analytical_derivatives) then
    j = 0.d0
    ndof = option%nflowdof
    ! let's be thorough and zero out the intermediates:
    D_xmfo = 0.d0; D_xmfg = 0.d0; D_mole_wt = 0.d0; D_den = 0.d0
    D_qsrc_mol = 0.d0; D_enthalpy = 0.d0; D_xmf = 0.d0
  endif

  ! this can be removed if extending to pressure condition
  if (.not.associated(src_sink_condition%rate) ) then
    option%io_buffer = 'TOWGBOSrcSink flow condition rate not defined ' // &
    'rate is needed for a valid src/sink term'
    call PrintErrMsg(option)
  endif

  qsrc => src_sink_condition%rate%dataset%rarray

  energy_var = 0
  if ( associated(src_sink_condition%temperature) ) then
    energy_var = SRC_TEMPERATURE
  else if ( associated(src_sink_condition%enthalpy) ) then
    energy_var = SRC_ENTHALPY
  endif
  
  flow_src_sink_type = src_sink_condition%rate%itype

  ! if not given, approximates BHP with pressure of perforated cell
  if ( associated(src_sink_condition%bhp_pressure) ) then
    cell_pressure = src_sink_condition%bhp_pressure%dataset%rarray(1)
    ! there should be no cell (oil) pressure dependence - zero out later on pressure derivs
    if (analytical_derivatives) then
      dp_dpo = 0.d0
    endif
  else
    cell_pressure = &
        maxval(auxvar%pres(option%liquid_phase:option%gas_phase))
    if (analytical_derivatives) then
      dp_dpo = 1.d0
      ! find dex corresponding to cell pres
      cp_loc = option%liquid_phase
      do idex = option%liquid_phase,option%gas_phase
        if (auxvar%pres(idex) > auxvar%pres(cp_loc)) then
          cp_loc = idex
        endif
      end do
      !cp_loc =  maxloc(auxvar%pres(option%liquid_phase:option%gas_phase))
      D_cpres = auxvar%D_pres(cp_loc,:)
    endif
  end if


  ! if enthalpy is used to define enthalpy or energy rate is used
  ! approximate bottom hole temperature (BHT) with local temp
  if ( associated(src_sink_condition%temperature) ) then
    temperature = src_sink_condition%temperature%dataset%rarray(1)
    ! there should be no cell temp dependence - zero out later on temp derivs
    if (analytical_derivatives) then
      dT_dTcell = 0.d0
    endif
  else
    temperature = auxvar%temp
    if (analytical_derivatives) then
      dT_dTcell = 1.d0
    endif
  end if

  Res = 0.d0
  do iphase = 1, option%nphase
    qsrc_mol = 0.d0
    mole_wt  = towg_fmw_comp(iphase)
    xmfo=auxvar%bo%xo
    xmfg=auxvar%bo%xg
    if (analytical_derivatives) then
      ! get xo and xg derivs from auxvars,
      D_xmfo = auxvar%bo%D_xo
      D_xmfg = auxvar%bo%D_xg
      ! mol wt has 0 derivs
      D_mole_wt = 0.d0
      ! zero out other derivatives
      D_den = 0.d0
    endif
    if ( qsrc(iphase) > 0.d0) then
      select case(option%phase_map(iphase))
        case(LIQUID_PHASE)
! Case of water
          if (.NOT. analytical_derivatives) then
            call EOSWaterDensity(temperature,cell_pressure,den_kg,den,ierr, &
                                                           auxvar%table_idx)
          else
            ! there is a density deriv w.r.t. pressure and temp
            call EOSWaterDensity(temperature,cell_pressure, &
                                 den_kg,den, &
                                 dx_dcpres, &
                                 dx_dtcell, ierr, auxvar%table_idx)
            D_den = D_cpres*dx_dcpres


          ! dx_dtcell may or may not be derivative by true 
          ! temperature, so we filter:
          dx_dtcell = dx_dtcell * dT_dTcell

          ! technically we might have to deal with a case where
          ! cell pressure is depenent on temp; this has been assumed
          ! possible since we use the full D_cpres array or derivatives,
          ! and is taken care of correctly above.
          ! However the density routine also introduces a seperate 
          ! dependence on temperature (though it may be fixed temp case thus
          ! zero etc etc)

          ! The cell pressure temp derivative is going to be 0 in any model we're dealing 
          ! with now, but for generality and correctness we should treat this
          ! properly.
          ! Let the density routine above be 
          ! den = den(a,b) ; it takes two arguments.
          ! We have selected a = cell pressure, b = temp:
          ! den = den(cp,t).
          ! Now we assume cp is a function of all nflowdof solution variables.
          ! Denote the cell variables as
          ! x;t
          ! where t is temperature and x are the non temperature variables
          ! we don't care about.
          ! the d den / d x_i are taken care of properly above. For derivative
          ! w.r.t. t we have a part (potentially) missing, since:
          ! d den / d t = (d den / d a) * (d cp / d t) + (d t / d t = 1) * (d den / d b)
          ! The part (d den / d cp) * (d cp / d t) was correctly taken care of
          ! by the line " D_den = D_cpres*dx_dcpres " above, so we just add in the
          ! missing part:
          D_den(dof_temp) = D_den(dof_temp) + dx_dtcell

          endif
        case(OIL_PHASE)
! Note this is dead oil, so take bubble point as reference pressure
          po=cell_pressure
          pb=ref_pressure
! Density and compressibility lookup
          if (.NOT. analytical_derivatives) then
            call EOSOilDensity        (temperature,pb,den,ierr,auxvar%table_idx)
            call EOSOilCompressibility(temperature,pb,cr ,ierr,auxvar%table_idx)
          else
            ! differentiate following correction 
            !crusp=cr*(po-pb)
            !den=den*(1.0+crusp*(1.0+0.5*crusp))
            ! given:
            !       - po is cell pres which should be treated as having full derivs
            !       - pb is ref pressure which is constant
            ! should care about the temp derivatives in the density and compressibility calls
            ! but not the pb ones:
            call EOSOilDensity(temperature,pb,den,dden_dt,dummy,ierr,auxvar%table_idx)
            call EOSOilCompressibility(temperature,pb,cr,dcr_dt,dummy,ierr,auxvar%table_idx)


            ! should also make sure dcr_dt is acted on by dT_dTcell:
            dcr_dt = dcr_dt * dT_dTcell
            dden_dt = dden_dt * dT_dTcell

            ! fully, we have:
            ! crusp(X) = cr(X) * (po(X) - pb)
            ! =>
            ! D_crusp(X) = D_cr(X) * (po(X) - pb) + cr(X) * D_po(X)
            ! The only nonzero term in D_cr would be the dof_temp one, so we 
            ! won't do full arrays for the first term of the above.
            ! the second term is:
            D_crusp = cr * D_cpres
            ! correction for the nonzero temp deriv of cr:
            D_crusp(dof_temp) = D_crusp(dof_temp) + dcr_dt * (po - pb)

            ! then we consider
            cor = 1 + crusp + 0.5*crusp*crusp
            ! dcor/dx = dcrusp/dx + crusp*dcrusp/dx  = (1 + crusp)*dcrusp/dx
            ! or
            D_cor = (1.d0 + crusp) * D_crusp
            ! Next:
            ! den = cor * den
            ! In principle full prod rule applies:
            ! D_den = cor * D_den + D_cor * den
            ! but only nonzero part of D_den is temp again so like before:
            D_den = den * D_cor
            ! correct temp part:
            D_cor(dof_temp) = D_cor(dof_temp) + dden_dt * cor
          endif
! Correct for undersaturation: correction not yet available for energy
          crusp=cr*(po-pb)
          den=den*(1.0+crusp*(1.0+0.5*crusp))
! Set to pure dead oil composition and molecular weight
          xmfo=1.0D0
          xmfg=0.0D0
          mole_wt=EOSOilGetFMW()

          if (analytical_derivatives) then
            ! now constants:
            D_xmfo = 0.d0
            D_xmfg = 0.d0
            D_mole_wt = 0.d0
          endif

        case(GAS_PHASE)
          if (.NOT. analytical_derivatives) then
            call EOSGasDensity(temperature,cell_pressure,den,ierr, &
                               auxvar%table_idx)
          else
            ! there is a density deriv w.r.t. pressure and temp at least
            call EOSGasDensity(temperature,cell_pressure,den,dx_dtcell,dx_dcpres, &
                               ierr,auxvar%table_idx)
            ! refer to explaination after analgous water density
            ! calls above
            dx_dtcell = dx_dtcell * dT_dTcell
            D_den = D_cpres * dx_dcpres
            D_den(dof_temp) = D_den(dof_temp) + dx_dtcell

          endif
        case(SOLVENT_PHASE)
          if (.NOT. analytical_derivatives) then
            call EOSSlvDensity(temperature,cell_pressure,den,ierr, &
                               auxvar%table_idx)
          else
          !!!EOSSlvDensityDerive(T,P,Rho_slv,dRho_dT,dRho_dP,ierr,table_idxs)
#if 0
            call EOSSlvDensity(temperature,cell_pressure,den,dx_dcpres,dx_dtcell, &
                               ierr,auxvar%table_idx)
#endif
            call EOSSlvDensity(temperature,cell_pressure,den,dx_dtcell,dx_dcpres, &
                               ierr,auxvar%table_idx)
            ! refer to explaination after analgous water density
            ! calls above
            dx_dtcell = dx_dtcell * dT_dTcell
            D_den = D_cpres * dx_dcpres
            D_den(dof_temp) = D_den(dof_temp) + dx_dtcell
          endif
      end select
    else
      den = auxvar%den(iphase)
      if (analytical_derivatives) then
        ! get den derivs from auxvar
        D_Den = auxvar%D_den(iphase,:)
      endif
      select case(option%phase_map(iphase))
        case(OIL_PHASE)
          mole_wt= xmfo*EOSOilGetFMW() &
                  +xmfg*EOSGasGetFMW()
          if (analytical_derivatives) then
            ! xmf derivs have been taken as the ones from 
            ! the auxvar, so just add and multiply by ...FMW()
            D_mole_wt = D_xmfo*EOSOilGetFMW() &
                      + D_xmfg*EOSGasGetFMW()
          endif
      end select
    end if

    select case(flow_src_sink_type)
      ! injection and production
      case(MASS_RATE_SS)
        qsrc_mol = qsrc(iphase)/mole_wt          ! kg/sec -> kmol/sec
        if (analytical_derivatives) then
          ! qsrc is constant, mole_wt might not be
          D_qsrc_mol = qsrc(iphase)*DivRule1(mole_wt,D_mole_wt,ndof)

        endif
      case(SCALED_MASS_RATE_SS)                  ! kg/sec -> kmol/sec
        qsrc_mol = qsrc(iphase)/mole_wt*scale
        ! 
        if (analytical_derivatives) then
          D_qsrc_mol = scale*qsrc(iphase)*DivRule1(mole_wt,D_mole_wt,ndof)
        endif
      case(VOLUMETRIC_RATE_SS)  ! assume local density for now
                                ! qsrc(iphase) = m^3/sec
        qsrc_mol = qsrc(iphase)*den ! den = kmol/m^3
        ! 
        if (analytical_derivatives) then
          !D_qsrc_mol = qsrc(iphase)*DivRule1(den,D_den,ndof)
          D_qsrc_mol = qsrc(iphase)*D_den
        endif
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec              ! den = kmol/m^3
        qsrc_mol = qsrc(iphase)* den * scale
        ! 
        if (analytical_derivatives) then
          !D_qsrc_mol = scale*qsrc(iphase)*DivRule1(den,D_den,ndof)
          D_qsrc_mol = scale*qsrc(iphase)*D_den
        endif
    end select

    ss_flow_vol_flux(iphase) = qsrc_mol/ den

! Loop over components (note number of comps=number of phases in black oil)
! Increment as given icomp value can occur more than once

    do icomp=1, option%nphase
      componentInPhase=checkBlackOilCIP(iphase,icomp,is_oil_in_oil,is_gas_in_oil,option)
      if( componentInPhase ) then
        xmf=1.0d0
        if( is_oil_in_oil ) xmf=xmfo
        if( is_gas_in_oil ) xmf=xmfg

        if (analytical_derivatives) then
          D_xmf = 0.d0
          if( is_oil_in_oil ) D_xmf=D_xmfo
          if( is_gas_in_oil ) D_xmf=D_xmfg
        endif

        Res(icomp) = Res(icomp)+xmf*qsrc_mol
        ! jac contribution here
        if (analytical_derivatives) then
          J(icomp,:) = J(icomp,:) + ProdRule(xmf,D_xmf,              &
                                             qsrc_mol,D_qsrc_mol,ndof )
        endif
      endif
    enddo

  enddo

  ! when using scaled src/sinks, the rates (mass or vol) scaling
  ! at this point the scale factor is already included in Res(iphase)

  ! Res(option%energy_id), energy units: MJ/sec

  !if ( associated(src_sink_condition%temperature) .or. &
  !    associated(src_sink_condition%enthalpy) &
  !   ) then
  !if the energy rate is not given, use either temperature or enthalpy
  if ( dabs(qsrc(option%energy_id)) < 1.0d-40 ) then
    ! water injection
    if (analytical_derivatives) then
      ! zero out enthalpy derivs
      D_enthalpy = 0.d0
    endif
    if (qsrc(option%liquid_phase) > 0.d0) then !implies qsrc(option%oil_phase)>=0
      if ( energy_var == SRC_ENTHALPY ) then
        !input as J/kg
        enthalpy = src_sink_condition%enthalpy% &
                       dataset%rarray(option%liquid_phase)
                     ! J/kg * kg/kmol = J/kmol
        enthalpy = enthalpy * towg_fmw_comp(option%liquid_phase)
      else !note: temp can either be input or taken as the one of perf. block
      !else if ( energy_var == SRC_TEMPERATURE ) then
        if (analytical_derivatives) then
          ! enthalpy derivatives from eos
          ! all anagous to density calls above
          call EOSWaterEnthalpy(temperature, cell_pressure,enthalpy,dx_dcpres,dx_dtcell,ierr)

          dx_dtcell = dx_dtcell * dT_dTcell
          D_enthalpy = dx_dcpres * D_cpres
          D_enthalpy(dof_temp) = D_enthalpy(dof_temp) + dx_dtcell
          D_enthalpy = D_enthalpy * 1.d-6
        else
          call EOSWaterEnthalpy(temperature, cell_pressure,enthalpy,ierr)
        endif
        ! enthalpy = [J/kmol]
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      ! enthalpy units: MJ/kmol

      if (analytical_derivatives) then
        ! jac contribution involves using previous residual values, so do it before those
        ! values change
        J(option%energy_id,:) = J(option%energy_id,:) + ProdRule(Res(option%liquid_phase),J(option%liquid_phase,:), &
                                                   enthalpy,D_enthalpy,ndof                           )
      endif

      ! enthalpy units: MJ/kmol ! water component mass
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%liquid_phase) * enthalpy
    end if
    ! oil injection (is assumed dead oil, so can use simple oil molecular weight)
    if (qsrc(option%oil_phase) > 0.d0) then !implies qsrc(option%liquid_phase)>=0
      if (analytical_derivatives) then
        ! zero out enthalpy derivs
        D_enthalpy = 0.d0
      endif
      if ( energy_var == SRC_ENTHALPY ) then
        enthalpy = src_sink_condition%enthalpy% &
                     dataset%rarray(option%oil_phase)
                      !J/kg * kg/kmol = J/kmol
        enthalpy = enthalpy * towg_fmw_comp(option%oil_phase)

        ! derivs should be 0 here
      else !note: temp can either be input or taken as the one of perf. block
      !if ( energy_var == SRC_TEMPERATURE ) then
        ! enthalpy = [J/kmol]
        if (analytical_derivatives) then
          ! enthalpy derivatives from eos
          ! all anagous to density calls above
          call EOSOilEnthalpy(temperature, cell_pressure,enthalpy,dx_dcpres,dx_dtcell,ierr)

          dx_dtcell = dx_dtcell * dT_dTcell
          D_enthalpy = dx_dcpres * D_cpres
          D_enthalpy(dof_temp) = D_enthalpy(dof_temp) + dx_dtcell
          D_enthalpy = D_enthalpy * 1.d-6
        else
          call EOSOilEnthalpy(temperature,cell_pressure,enthalpy,ierr)
        endif
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      ! enthalpy units: MJ/kmol
      if (analytical_derivatives) then
        ! jac contribution involves using previous residual values, so do it before those
        ! values change
        J(option%energy_id,:) = J(option%energy_id,:) + ProdRule(Res(option%oil_phase),J(option%oil_phase,:), &
                                                   enthalpy,D_enthalpy,ndof                           )
      endif
      ! enthalpy units: MJ/kmol ! oil component mass
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%oil_phase) * enthalpy
    end if
    ! gas injection
    if (qsrc(option%gas_phase) > 0.d0) then
      if (analytical_derivatives) then
        ! zero out enthalpy derivs
        D_enthalpy = 0.d0
      endif
      if ( energy_var == SRC_ENTHALPY ) then
        enthalpy = src_sink_condition%enthalpy% &
                     dataset%rarray(option%gas_phase)
                      !J/kg * kg/kmol = J/kmol
        enthalpy = enthalpy * towg_fmw_comp(option%gas_phase)
      else !note: temp can either be input or taken as the one of perf. block
      !if ( energy_var == SRC_TEMPERATURE ) then
        ! enthalpy = [J/kmol]
        if (analytical_derivatives) then
          ! enthalpy derivatives from eos
          ! all anagous to density calls above
          call EOSGasEnergy(temperature, cell_pressure,enthalpy,dx_dtcell,dx_dcpres,&
                            internal_energy_dummy,dum1,dum2,ierr)

          dx_dtcell = dx_dtcell * dT_dTcell
          D_enthalpy = dx_dcpres * D_cpres
          D_enthalpy(dof_temp) = D_enthalpy(dof_temp) + dx_dtcell
          D_enthalpy = D_enthalpy * 1.d-6

        else
          call EOSGasEnergy(temperature,cell_pressure,enthalpy, &
                                internal_energy_dummy,ierr)
        endif
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      ! enthalpy units: MJ/kmol ! oil component mass
      if (analytical_derivatives) then
        ! jac contribution involves using previous residual values, so do it before those
        ! values change
        J(option%energy_id,:) = J(option%energy_id,:) + ProdRule(Res(option%gas_phase),J(option%gas_phase,:), &
                                                   enthalpy,D_enthalpy,ndof                           )
      endif
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%gas_phase) * enthalpy
    end if

!--Solvent injection-----------------------------------------------------------

    if( towg_miscibility_model == TOWG_SOLVENT_TL ) then
      if (qsrc(option%solvent_phase) > 0.d0) then
        if ( energy_var == SRC_ENTHALPY ) then
          enthalpy = src_sink_condition%enthalpy% &
                       dataset%rarray(option%solvent_phase)
                        !J/kg * kg/kmol = J/kmol
          enthalpy = enthalpy * towg_fmw_comp(option%solvent_phase)
        else !note: temp can either be input or taken as the one of perf. block
        !if ( energy_var == SRC_TEMPERATURE ) then
          if (analytical_derivatives) then
!subroutine EOSSlvEnergyDerive(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
            call EOSSlvEnergy(temperature,cell_pressure,enthalpy,dx_dtcell,dx_dcpres, &
                              internal_energy_dummy,dummy,dummy,ierr)

            dx_dtcell = dx_dtcell * dT_dTcell
            D_enthalpy = dx_dcpres * D_cpres
            D_enthalpy(dof_temp) = D_enthalpy(dof_temp) + dx_dtcell
            D_enthalpy = D_enthalpy * 1.d-6
          else
            call EOSSlvEnergy(temperature,cell_pressure,enthalpy, &
                              internal_energy_dummy,ierr)
            ! enthalpy = [J/kmol]
          endif
        end if
        enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
        ! enthalpy units: MJ/kmol
        if (analytical_derivatives) then
          ! jac contribution involves using previous residual values, so do it before those
          ! values change
          J(option%energy_id,:) = J(option%energy_id,:) + ProdRule(Res(option%solvent_phase),J(option%solvent_phase,:), &
                                                     enthalpy,D_enthalpy,ndof                           )
        endif
        Res(option%energy_id) = Res(option%energy_id) + &
                                Res(option%solvent_phase) * enthalpy
      end if
    endif

    ! water energy extraction due to water production
    if (qsrc(option%liquid_phase) < 0.d0) then !implies qsrc(option%oil_phase)<=0
      ! auxvar enthalpy units: MJ/kmol
      if (analytical_derivatives) then
        ! jac contribution involves using previous residual values, so do it before those
        ! values change
        J(option%energy_id,:) = J(option%energy_id,:)                                                           &
                              + ProdRule(Res(option%liquid_phase),J(option%liquid_phase,:),                     &
                                         auxvar%H(option%liquid_phase),auxvar%D_H(option%liquid_phase,:),ndof )
      endif
      ! auxvar enthalpy units: MJ/kmol ! water component mass
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%liquid_phase) * &
                              auxvar%H(option%liquid_phase)
    end if
    !oil energy extraction due to oil production
    if (qsrc(option%oil_phase) < 0.d0) then !implies qsrc(option%liquid_phase)<=0
      ! auxvar enthalpy units: MJ/kmol
      if (analytical_derivatives) then
        ! jac contribution involves using previous residual values, so do it before those
        ! values change
        J(option%energy_id,:) = J(option%energy_id,:) + ProdRule(Res(option%oil_phase),J(option%oil_phase,:),                   &
                                                                 auxvar%H(option%oil_phase),auxvar%D_H(option%oil_phase,:),ndof )
      endif
      ! auxvar enthalpy units: MJ/kmol ! water component mass
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%oil_phase) * &
                              auxvar%H(option%oil_phase)
    end if
    if (qsrc(option%gas_phase) < 0.d0) then !implies qsrc(option%liquid_phase)<=0
      ! auxvar enthalpy units: MJ/kmol
      if (analytical_derivatives) then
        ! jac contribution involves using previous residual values, so do it before those
        ! values change
        J(option%energy_id,:) = J(option%energy_id,:) + ProdRule(Res(option%gas_phase),J(option%gas_phase,:),                   &
                                                                 auxvar%H(option%gas_phase),auxvar%D_H(option%gas_phase,:),ndof )
      endif
      ! auxvar enthalpy units: MJ/kmol ! water component mass
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%gas_phase) * &
                              auxvar%H(option%gas_phase)
    end if

    !--Solvent energy extraction due to solvent production-------------------------
    if( towg_miscibility_model == TOWG_SOLVENT_TL ) then

      if (qsrc(option%solvent_phase) < 0.d0) then !implies qsrc(option%gas_phase)<=0
        ! auxvar enthalpy units: MJ/kmol

        ! auxvar enthalpy units: MJ/kmol
        if (analytical_derivatives) then
          ! jac contribution involves using previous residual values, so do it before those
          ! values change
          J(option%energy_id,:) = J(option%energy_id,:)                                             &
                                + ProdRule(Res(option%solvent_phase),J(option%solvent_phase,:),     &
                                           auxvar%H(option%solvent_phase),auxvar%D_H(option%solvent_phase,:),ndof )
        endif

        Res(option%energy_id) = Res(option%energy_id) + &
                                Res(option%solvent_phase) * &
                                auxvar%H(option%solvent_phase)
      end if
    endif

  else !if the energy rate is given, it overwrites both temp and enthalpy
    ! if energy rate is given, loaded in qsrc(option%energy_id) in MJ/sec
    Res(option%energy_id) = qsrc(option%energy_id)* scale ! MJ/s
  end if

  nullify(qsrc)

end subroutine TOWGBOSrcSink

! ************************************************************************** !

subroutine TOWGAccumDerivative(auxvar,global_auxvar,material_auxvar, &
                               soil_heat_capacity,option,J,ghosted_id)
  !
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  !
  ! Author: Paolo Orsini
  ! Date: 12/27/16
  !

  use Option_module
  use Saturation_Function_module
  use Material_Aux_class
  use Utility_module

  implicit none

  type(auxvar_towg_type) :: auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  PetscReal :: J(option%nflowdof,option%nflowdof)

  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: jac(option%nflowdof,option%nflowdof)
  PetscReal :: jac_pert(option%nflowdof,option%nflowdof)
  PetscInt :: idof, irow, idex, jdex

  PetscReal :: jdum(option%nflowdof,option%nflowdof)
  PetscReal :: jalyt(option%nflowdof,option%nflowdof)
  PetscBool :: flagged
  PetscInt :: ghosted_id

  if (.NOT. towg_analytical_derivatives .OR. towg_analytical_derivatives_compare) then
    call TOWGAccumulation(auxvar(ZERO_INTEGER),global_auxvar,material_auxvar, &
                          soil_heat_capacity,option,res,PETSC_FALSE,jdum,PETSC_FALSE)

    do idof = 1, option%nflowdof
      call TOWGAccumulation(auxvar(idof),global_auxvar,material_auxvar, &
                            soil_heat_capacity,option,res_pert,PETSC_FALSE,jdum,PETSC_FALSE)
      do irow = 1, option%nflowdof
        J(irow,idof) = (res_pert(irow)-res(irow))/auxvar(idof)%pert
      enddo !irow
    enddo ! idof

    if (towg_isothermal) then
      J(towg_energy_eq_idx,:) = 0.d0
      J(:,towg_energy_eq_idx) = 0.d0
    endif

    if (towg_no_oil) then
      J(TOWG_OIL_EQ_IDX,:) = 0.d0
      J(:,TOWG_OIL_EQ_IDX) = 0.d0
    endif

    if (towg_no_gas) then
      J(TOWG_GAS_EQ_IDX,:) = 0.d0
      J(:,TOWG_GAS_EQ_IDX) = 0.d0
    endif
  endif

  if (towg_analytical_derivatives) then
    call TOWGAccumulation(auxvar(ZERO_INTEGER),global_auxvar,material_auxvar, &
                          soil_heat_capacity,option,res,PETSC_FALSE,jalyt,PETSC_TRUE)

#if 0
!!! we may wish to have an option for a robust checker for this at some point, it can happen surpringly 
!!! easilly with marginal saturations
    do idex = 1,option%nflowdof
      do jdex = 1,option%nflowdof
        if (isnan(Jalyt(idex,jdex))) then
          print *, "NAN HERE! accum", Jalyt(idex,jdex), idex, jdex
        endif
      enddo
    enddo
#endif


    if (towg_isothermal) then
      jalyt(towg_energy_eq_idx,:) = 0.d0
      jalyt(:,towg_energy_eq_idx) = 0.d0
    endif

    if (towg_no_oil) then
      jalyt(TOWG_OIL_EQ_IDX,:) = 0.d0
      jalyt(:,TOWG_OIL_EQ_IDX) = 0.d0
    endif

    if (towg_no_gas) then
      jalyt(TOWG_GAS_EQ_IDX,:) = 0.d0
      jalyt(:,TOWG_GAS_EQ_IDX) = 0.d0
    endif

    if (towg_analytical_derivatives_compare) then
      flagged = PETSC_FALSE
      call MatCompare(J, jalyt, option%nflowdof, option%nflowdof, towg_dcomp_tol, towg_dcomp_reltol,flagged)
      if (flagged) then
        print *, "this is accum derivative, ghosted_id: ", ghosted_id
        print *, auxvar(ZERO_INTEGER)%sat
      endif
    endif

    j = jalyt


  endif


end subroutine TOWGAccumDerivative

! ************************************************************************** !

subroutine TOWGFluxDerivative(auxvar_up,global_auxvar_up, &
                              material_auxvar_up, &
                              thermal_conductivity_up, &
                              auxvar_dn,global_auxvar_dn, &
                              material_auxvar_dn, &
                              thermal_conductivity_dn, &
                              area, dist, &
                              towg_parameter, &
                              option,Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/27/16
  ! 
  use Option_module
  use Material_Aux_class
  use Utility_module
  
  implicit none
  
  type(auxvar_towg_type) :: auxvar_up(0:), auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(towg_parameter_type) :: towg_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  PetscReal :: Jalyt_up(option%nflowdof,option%nflowdof)
  PetscReal :: Jalyt_dn(option%nflowdof,option%nflowdof)

  PetscReal :: Jdum_up(option%nflowdof,option%nflowdof)
  PetscReal :: Jdum_dn(option%nflowdof,option%nflowdof)

  PetscBool :: flagged

  PetscInt :: i,j

  Jup = 0.d0
  Jdn = 0.d0


  if (.NOT. towg_analytical_derivatives .OR. towg_analytical_derivatives_compare) then
  
    !print *, 'TOWGFluxDerivative'
    option%iflag = -2
    call TOWGFlux(auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                  material_auxvar_up, &
                  thermal_conductivity_up, &
                  auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn, &
                  thermal_conductivity_dn, &
                  area,dist,towg_parameter, &
                  option,v_darcy,res,PETSC_FALSE,Jdum_up,Jdum_dn,PETSC_FALSE)
                             
    ! upgradient derivatives
    do idof = 1, option%nflowdof
      call TOWGFlux(auxvar_up(idof),global_auxvar_up, &
                    material_auxvar_up, &
                    thermal_conductivity_up, &
                    auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                    material_auxvar_dn, &
                    thermal_conductivity_dn, &
                    area,dist,towg_parameter, &
                    option,v_darcy,res_pert,PETSC_FALSE,Jdum_up,Jdum_dn,PETSC_FALSE)
      do irow = 1, option%nflowdof
        Jup(irow,idof) = (res_pert(irow)-res(irow))/auxvar_up(idof)%pert
      enddo !irow
    enddo ! idof

    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call TOWGFlux(auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                    material_auxvar_up, &
                    thermal_conductivity_up, &
                    auxvar_dn(idof),global_auxvar_dn, &
                    material_auxvar_dn, &
                    thermal_conductivity_dn, &
                    area,dist,towg_parameter, &
                    option,v_darcy,res_pert,PETSC_FALSE,Jdum_up,Jdum_dn,PETSC_FALSE)
      do irow = 1, option%nflowdof
        Jdn(irow,idof) = (res_pert(irow)-res(irow))/auxvar_dn(idof)%pert
      enddo !irow
    enddo ! idof

    if (towg_isothermal) then
      Jup(towg_energy_eq_idx,:) = 0.d0
      Jup(:,towg_energy_eq_idx) = 0.d0
      Jdn(towg_energy_eq_idx,:) = 0.d0
      Jdn(:,towg_energy_eq_idx) = 0.d0
    endif
    
    if (towg_no_oil) then
      Jup(TOWG_OIL_EQ_IDX,:) = 0.d0
      Jup(:,TOWG_OIL_EQ_IDX) = 0.d0
      Jdn(TOWG_OIL_EQ_IDX,:) = 0.d0
      Jdn(:,TOWG_OIL_EQ_IDX) = 0.d0
    endif  

    if (towg_no_gas) then
      Jup(TOWG_GAS_EQ_IDX,:) = 0.d0
      Jup(:,TOWG_GAS_EQ_IDX) = 0.d0
      Jdn(TOWG_GAS_EQ_IDX,:) = 0.d0
      Jdn(:,TOWG_GAS_EQ_IDX) = 0.d0
    endif  
  endif

  if (towg_analytical_derivatives) then
    call TOWGFlux(auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                  material_auxvar_up, &
                  thermal_conductivity_up, &
                  auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn, &
                  thermal_conductivity_dn, &
                  area,dist,towg_parameter, &
                  option,v_darcy,res,PETSC_FALSE,Jalyt_up,Jalyt_dn,PETSC_TRUE)

#if 0
!!! we may wish to have an option for a robust checker for this at some point, it can happen surpringly 
!!! easilly with marginal saturations
    do i = 1,option%nflowdof
      do j = 1,option%nflowdof
        if (isnan(Jalyt_up(i,j))) then
          print *, "NAN HERE! (flux up)", Jalyt_up(i,j), i, j
          !Jalyt_up(i,j) = 0.d0
        endif
        if (isnan(Jalyt_dn(i,j))) then
          print *, "NAN HERE! (flux dn)", Jalyt_dn(i,j), i, j
          !Jalyt_dn(i,j) = 0.d0
        endif
      enddo
    enddo
#endif

    if (towg_isothermal) then
      Jalyt_up(towg_energy_eq_idx,:) = 0.d0
      Jalyt_up(:,towg_energy_eq_idx) = 0.d0
      Jalyt_dn(towg_energy_eq_idx,:) = 0.d0
      Jalyt_dn(:,towg_energy_eq_idx) = 0.d0
    endif
    
    if (towg_no_oil) then
      Jalyt_up(TOWG_OIL_EQ_IDX,:) = 0.d0
      Jalyt_up(:,TOWG_OIL_EQ_IDX) = 0.d0
      Jalyt_dn(TOWG_OIL_EQ_IDX,:) = 0.d0
      Jalyt_dn(:,TOWG_OIL_EQ_IDX) = 0.d0
    endif  

    if (towg_no_gas) then
      Jalyt_up(TOWG_GAS_EQ_IDX,:) = 0.d0
      Jalyt_up(:,TOWG_GAS_EQ_IDX) = 0.d0
      Jalyt_dn(TOWG_GAS_EQ_IDX,:) = 0.d0
      Jalyt_dn(:,TOWG_GAS_EQ_IDX) = 0.d0
    endif  

    if (towg_analytical_derivatives_compare) then
      flagged = PETSC_FALSE
      call MatCompare(Jup, Jalyt_up, option%nflowdof,option%nflowdof, towg_dcomp_tol, towg_dcomp_reltol,flagged)
      if (flagged) then
        print *, "this is flux derivative, that was matrix up"
      endif
      flagged = PETSC_FALSE
      call MatCompare(Jdn, Jalyt_dn, option%nflowdof, option%nflowdof, towg_dcomp_tol, towg_dcomp_reltol,flagged)
      if (flagged) then
        print *, "this is flux derivative, that was matrix dn"
      endif

    call TOWGFlux(auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                  material_auxvar_up, &
                  thermal_conductivity_up, &
                  auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn, &
                  thermal_conductivity_dn, &
                  area,dist,towg_parameter, &
                  option,v_darcy,res,PETSC_FALSE,Jalyt_up,Jalyt_dn,PETSC_TRUE)

    call TOWGFlux(auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                  material_auxvar_up, &
                  thermal_conductivity_up, &
                  auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn, &
                  thermal_conductivity_dn, &
                  area,dist,towg_parameter, &
                  option,v_darcy,res,PETSC_FALSE,Jdum_up,Jdum_dn,PETSC_FALSE)

      idof = 4
      call TOWGFlux(auxvar_up(idof),global_auxvar_up, &
                    material_auxvar_up, &
                    thermal_conductivity_up, &
                    auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                    material_auxvar_dn, &
                    thermal_conductivity_dn, &
                    area,dist,towg_parameter, &
                    option,v_darcy,res_pert,PETSC_FALSE,Jdum_up,Jdum_dn,PETSC_FALSE)

    endif

    jup = jalyt_up
    jdn = jalyt_dn
  endif


end subroutine TOWGFluxDerivative

! ************************************************************************** !

subroutine TOWGBCFluxDerivative(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                                auxvar_up, global_auxvar_up, &
                                auxvar_dn,global_auxvar_dn, &
                                material_auxvar_dn, &
                                thermal_conductivity_dn, &
                                area,dist,towg_parameter, &
                                option,Jdn)
  ! 
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/27/16
  ! 

  use Option_module 
  use Material_Aux_class
  use Utility_module
  
  implicit none

  PetscReal :: bc_auxvars(:) ! from aux_real_var array
  type(auxvar_towg_type) :: auxvar_up, auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(towg_parameter_type) :: towg_parameter
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: bc_auxvar_mapping(TOWG_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  PetscReal :: jdum(option%nflowdof,option%nflowdof)
  PetscReal :: Jalyt_dn(option%nflowdof,option%nflowdof)

  PetscBool :: flagged

  Jdn = 0.d0
  !print *, 'TOWGBCFluxDerivative'

  if (.NOT. towg_analytical_derivatives .OR. towg_analytical_derivatives_compare) then

    option%iflag = -2
    call TOWGBCFlux(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                    auxvar_up,global_auxvar_up, &
                    auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                    material_auxvar_dn, &
                    thermal_conductivity_dn, &
                    area,dist,towg_parameter, &
                    option,v_darcy,res,PETSC_FALSE,jdum,PETSC_FALSE)
                      
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call TOWGBCFlux(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                      auxvar_up,global_auxvar_up, &
                      auxvar_dn(idof),global_auxvar_dn, &
                      material_auxvar_dn, &
                      thermal_conductivity_dn, &
                      area,dist,towg_parameter, &
                      option,v_darcy,res_pert,PETSC_FALSE,jdum,PETSC_FALSE)
      do irow = 1, option%nflowdof
        Jdn(irow,idof) = (res_pert(irow)-res(irow))/auxvar_dn(idof)%pert
      enddo !irow
    enddo ! idof

    if (towg_isothermal) then
      Jdn(towg_energy_eq_idx,:) = 0.d0
      Jdn(:,towg_energy_eq_idx) = 0.d0
    endif
    
    if (towg_no_oil) then
      Jdn(TOWG_OIL_EQ_IDX,:) = 0.d0
      Jdn(:,TOWG_OIL_EQ_IDX) = 0.d0
    endif  

    if (towg_no_gas) then
      Jdn(TOWG_GAS_EQ_IDX,:) = 0.d0
      Jdn(:,TOWG_GAS_EQ_IDX) = 0.d0
    endif  

  endif

  if (towg_analytical_derivatives) then
    call TOWGBCFlux(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                    auxvar_up,global_auxvar_up, &
                    auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                    material_auxvar_dn, &
                    thermal_conductivity_dn, &
                    area,dist,towg_parameter, &
                    option,v_darcy,res,PETSC_FALSE,jalyt_dn,PETSC_TRUE)

    if (towg_isothermal) then
     jalyt_dn(towg_energy_eq_idx,:) = 0.d0
     jalyt_dn(:,towg_energy_eq_idx) = 0.d0
    endif
    
    if (towg_no_oil) then
     jalyt_dn(TOWG_OIL_EQ_IDX,:) = 0.d0
     jalyt_dn(:,TOWG_OIL_EQ_IDX) = 0.d0
    endif  

    if (towg_no_gas) then
      jalyt_dn(TOWG_GAS_EQ_IDX,:) = 0.d0
      jalyt_dn(:,TOWG_GAS_EQ_IDX) = 0.d0
    endif  

    if (towg_analytical_derivatives_compare) then
      flagged = PETSC_FALSE
      call MatCompare(Jdn, Jalyt_dn, option%nflowdof, option%nflowdof, towg_dcomp_tol, towg_dcomp_reltol,flagged)
      if (flagged) then
        print *, "this is bc flux derivative"
      endif

      call TOWGBCFlux(ibndtype,bc_auxvar_mapping,bc_auxvars, &
                      auxvar_up,global_auxvar_up, &
                      auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                      material_auxvar_dn, &
                      thermal_conductivity_dn, &
                      area,dist,towg_parameter, &
                      option,v_darcy,res,PETSC_FALSE,jalyt_dn,PETSC_TRUE)
    endif

    jdn = jalyt_dn
  endif
    
  end subroutine TOWGBCFluxDerivative

! ************************************************************************** !

subroutine TOWGSrcSinkDerivative(option,src_sink_condition,auxvars, &
                                 global_auxvar,scale,Jac)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 12/28/16
  ! 

  use Option_module
  use Condition_module
  use Utility_module

  implicit none

  type(option_type) :: option
  type(flow_towg_condition_type), pointer :: src_sink_condition
  type(auxvar_towg_type) :: auxvars(0:)
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: scale
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow

  PetscReal :: jdum(option%nflowdof,option%nflowdof)
  PetscReal :: Jalyt(option%nflowdof,option%nflowdof)

  PetscBool :: flagged


  if (.NOT. towg_analytical_derivatives .OR. towg_analytical_derivatives_compare) then
    option%iflag = -3
    call TOWGSrcSink(option,src_sink_condition,auxvars(ZERO_INTEGER), &
                     global_auxvar,dummy_real,scale,Res,jdum,PETSC_FALSE)

    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call TOWGSrcSink(option,src_sink_condition,auxvars(idof), &
                       global_auxvar,dummy_real,scale,res_pert,jdum,PETSC_FALSE)
      do irow = 1, option%nflowdof
        Jac(irow,idof) = (res_pert(irow)-res(irow))/auxvars(idof)%pert
      enddo !irow
    enddo ! idof
   
    if (towg_isothermal) then
      Jac(towg_energy_eq_idx,:) = 0.d0
      Jac(:,towg_energy_eq_idx) = 0.d0
    endif
    
    if (towg_no_oil) then
      Jac(TOWG_OIL_EQ_IDX,:) = 0.d0
      Jac(:,TOWG_OIL_EQ_IDX) = 0.d0
    endif  

    if (towg_no_gas) then
      Jac(TOWG_GAS_EQ_IDX,:) = 0.d0
      Jac(:,TOWG_GAS_EQ_IDX) = 0.d0
    endif  
  endif

  if (towg_analytical_derivatives) then

    call TOWGSrcSink(option,src_sink_condition,auxvars(ZERO_INTEGER), &
                     global_auxvar,dummy_real,scale,res_pert,Jalyt,PETSC_TRUE)

    if (towg_isothermal) then
      jalyt(towg_energy_eq_idx,:) = 0.d0
      jalyt(:,towg_energy_eq_idx) = 0.d0
    endif
    
    if (towg_no_oil) then
     jalyt(TOWG_OIL_EQ_IDX,:) = 0.d0
     jalyt(:,TOWG_OIL_EQ_IDX) = 0.d0
    endif  

    if (towg_no_gas) then
      jalyt(TOWG_GAS_EQ_IDX,:) = 0.d0
      jalyt(:,TOWG_GAS_EQ_IDX) = 0.d0
    endif  

    if (towg_analytical_derivatives_compare) then
      flagged = PETSC_FALSE
      call MatCompare(Jac, Jalyt,option%nflowdof,option%nflowdof, towg_dcomp_tol, towg_dcomp_reltol,flagged)
      if (flagged) then 
        print *, "this is src sink derivative"
      endif
      call TOWGSrcSink(option,src_sink_condition,auxvars(ZERO_INTEGER), &
                       global_auxvar,dummy_real,scale,res_pert,Jalyt,PETSC_TRUE)
    endif

    jac = jalyt
  endif

end subroutine TOWGSrcSinkDerivative

! ************************************************************************** !

subroutine TOWGResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 12/28/16
  ! 
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscmat.h"
  use petscsnes
  use petscmat
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
  use Well_Data_class
  use Well_Solver_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr,jerr
  
  Mat, parameter :: null_mat = PETSC_NULL_MAT
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(material_parameter_type), pointer :: material_parameter
  class(pm_towg_aux_type), pointer :: towg
  type(towg_parameter_type), pointer :: towg_parameter
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscInt :: iphase
  PetscReal :: scale
  PetscReal :: ss_flow_vol_flux(realization%option%nphase)
  PetscInt :: sum_connection
  PetscInt :: local_start, local_end
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, imat, imat_up, imat_dn
  PetscInt, save :: iplot = 0

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word

  PetscInt :: icap_up, icap_dn
  PetscReal :: Res(realization%option%nflowdof)
  !PetscReal :: Jac_dummy(realization%option%nflowdof, &
  !                       realization%option%nflowdof)
  PetscReal :: v_darcy(realization%option%nphase)

  type(well_data_list_type),pointer :: well_data_list
  class(well_data_type), pointer :: well_data

  PetscReal :: jdum(realization%option%nflowdof,realization%option%nflowdof)
  PetscReal :: jdum2(realization%option%nflowdof,realization%option%nflowdof)
  
  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  towg => patch%aux%TOWG
  towg_parameter => towg%parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  
#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then
    debug_iteration_count = debug_iteration_count + 1
    write(word,*) debug_timestep_count
    string = 'residual_debug_data_' // trim(adjustl(word))
    write(word,*) debug_timestep_cut_count
    string = trim(string) // '_' // trim(adjustl(word))
    write(word,*) debug_iteration_count
    debug_filename = trim(string) // '_' // trim(adjustl(word)) // '.txt'
    open(debug_unit, file=debug_filename, action="write", status="unknown")
    open(debug_info_unit, file='debug_info.txt', action="write", &
         position="append", status="unknown")
    write(debug_info_unit,*) 'residual ', debug_timestep_count, &
      debug_timestep_cut_count, debug_iteration_count
    close(debug_info_unit)
  endif
#endif

  ! Communication -----------------------------------------
  ! These 3 must be called before TOWGUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  
                                     ! do update state
  call TOWGUpdateAuxVars(realization,PETSC_TRUE)

  ! override flags since they will soon be out of date
  towg%auxvars_up_to_date = PETSC_FALSE 

  ! always assume variables have been swapped; therefore, must copy back
  call VecLockPop(xx,ierr); CHKERRQ(ierr)
  call DiscretizationLocalToGlobal(discretization,field%flow_xx_loc,xx, &
                                   NFLOWDOF)
  call VecLockPush(xx,ierr); CHKERRQ(ierr)

  if (option%compute_mass_balance_new) then
    call TOWGZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  r_p = -accum_p

  
  !Heeho dynamically update p+1 accumulation term
  !if (towg_tough2_conv_criteria) then
  !  call VecGetArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  !endif
  
  ! accumulation at t(k+1)
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    call TOWGAccumulation(towg%auxvars(ZERO_INTEGER,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          material_auxvars(ghosted_id), &
                          material_parameter%soil_heat_capacity(imat), &
                          option,Res,local_id == towg_debug_cell_id,jdum,PETSC_FALSE)


    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    
    !Heeho dynamically update p+1 accumulation term
    !if (towg_tough2_conv_criteria) then
    !  accum_p2(local_start:local_end) = Res(:)
    !endif
    
  enddo

  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  !Heeho dynamically update p+1 accumulation term
  !if (towg_tough2_conv_criteria) then
  !  call VecRestoreArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  !endif

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

      call TOWGFlux(towg%auxvars(ZERO_INTEGER,ghosted_id_up), &
                    global_auxvars(ghosted_id_up), &
                    material_auxvars(ghosted_id_up), & 
                    material_parameter%soil_thermal_conductivity(:,imat_up), &
                    towg%auxvars(ZERO_INTEGER,ghosted_id_dn), &
                    global_auxvars(ghosted_id_dn), &
                    material_auxvars(ghosted_id_dn), &
                    material_parameter%soil_thermal_conductivity(:,imat_dn), &
                    cur_connection_set%area(iconn), &
                    cur_connection_set%dist(:,iconn), &
                    towg_parameter,option,v_darcy,Res, &
                    (local_id_up == towg_debug_cell_id .or. &
                     local_id_dn == towg_debug_cell_id), &
                     jdum,jdum2,PETSC_FALSE)

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

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icap_dn = patch%sat_func_id(ghosted_id)

      call TOWGBCFlux(boundary_condition%flow_bc_type, &
                    boundary_condition%flow_aux_mapping, & 
                    boundary_condition%flow_aux_real_var(:,iconn), &
                    towg%auxvars_bc(sum_connection), &
                    global_auxvars_bc(sum_connection), &
                    towg%auxvars(ZERO_INTEGER,ghosted_id), &
                    global_auxvars(ghosted_id), &
                    material_auxvars(ghosted_id), &
                    material_parameter%soil_thermal_conductivity(:,imat_dn), &
                    cur_connection_set%area(iconn), &
                    cur_connection_set%dist(:,iconn), &
                    towg_parameter,option,v_darcy,Res, &
                    local_id == towg_debug_cell_id,jdum,PETSC_FALSE)

      patch%boundary_velocities(:,sum_connection) = v_darcy
      if (associated(patch%boundary_flow_fluxes)) then
        patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)% &
            mass_balance_delta(1:option%nflowspec,1) = &
          global_auxvars_bc(sum_connection)% &
            mass_balance_delta(1:option%nflowspec,1) - &
          Res(1:option%nflowspec)
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

      call TOWGSrcSink(option,source_sink%flow_condition%towg, &
                       towg%auxvars(ZERO_INTEGER,ghosted_id), &
                       global_auxvars(ghosted_id), ss_flow_vol_flux, &
                       scale,Res,jdum,PETSC_FALSE)

      r_p(local_start:local_end) =  r_p(local_start:local_end) - Res(:)

      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif      
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)
      endif      
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_ss(sum_connection)% &
            mass_balance_delta(1:option%nflowspec,1) = &
          global_auxvars_ss(sum_connection)% &
            mass_balance_delta(1:option%nflowspec,1) - &
          Res(1:option%nflowspec)
      endif

    enddo
    source_sink => source_sink%next
  enddo

  ! Set up the average pressure (may be needed for voidage calculations)

  call patch%aux%TOWG%FieldVolRefAve(grid,patch%aux%material, &
                                          patch%imat,option)

  ! Loop over well_data wells if present

  if (WellDataGetFlag()) then
    jerr = 0
    well_data_list => realization%well_data
    well_data => well_data_list%first

    do
      if (.not.associated(well_data)) exit
        call SolveWell(patch%aux,option,well_data,r_p)
        call MPI_Barrier(option%mycomm,jerr)
      well_data => well_data%next
    enddo
  endif

  if (towg%inactive_cells_exist) then
    do i = 1,towg%n_inactive_rows
      r_p(towg%inactive_rows_local(i)) = 0.d0
    enddo
  endif
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  
  if (Initialized(towg_debug_cell_id)) then
    call VecGetArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
    do local_id = towg_debug_cell_id-1, towg_debug_cell_id+1
      write(*,'(''  residual   : '',i2,10es12.4)') local_id, &
        r_p((local_id-1)*option%nflowdof+1:(local_id-1)*option%nflowdof+2), &
        r_p(local_id*option%nflowdof)*1.d6
    enddo
    call VecRestoreArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
  endif
  
  if (towg_isothermal) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+towg_energy_eq_idx) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif
  if (towg_no_oil) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+TOWG_OIL_EQ_IDX) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif  
  if (towg_no_gas) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+TOWG_GAS_EQ_IDX) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif  

#ifdef DEBUG_TOWG_FILEOUTPUT
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    write(debug_unit,'(a,i5,7es24.15)') 'fixed residual:', local_id, &
      accum_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof)
  enddo
  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    write(debug_unit,'(a,i5,7es24.15)') 'residual:', local_id, &
      r_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof)
  enddo
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
#endif
  
  if (realization%debug%vecview_residual) then
    string = 'TOWGresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    string = 'TOWGxx'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then
    close(debug_unit)
  endif
#endif
  
end subroutine TOWGResidual

! ************************************************************************** !

subroutine TOWGJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/28/16
  ! 

#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscmat.h"
  use petscsnes
  use petscmat
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_class
  use Well_Data_class

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr,jerr,nflowdof

  Mat :: J
  MatType :: mat_type
  PetscReal :: norm
  PetscViewer :: viewer

  PetscInt :: icap_up,icap_dn
  PetscReal :: qsrc, scale
  PetscInt :: imat, imat_up, imat_dn
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: irow
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  Vec, parameter :: null_vec = PETSC_NULL_VEC

  class(well_data_type), pointer :: well_data
  type(well_data_list_type),pointer :: well_data_list

  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  PetscInt, pointer :: zeros(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(material_parameter_type), pointer :: material_parameter
  class(pm_towg_aux_type), pointer :: towg
  type(towg_parameter_type), pointer :: towg_parameter
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscBool::analytical_derivatives
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  towg => patch%aux%TOWG
  towg_parameter => towg%parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then
    write(word,*) debug_timestep_count
    string = 'jacobian_debug_data_' // trim(adjustl(word))
    write(word,*) debug_timestep_cut_count
    string = trim(string) // '_' // trim(adjustl(word))
    write(word,*) debug_iteration_count
    debug_filename = trim(string) // '_' // trim(adjustl(word)) // '.txt'
    open(debug_unit, file=debug_filename, action="write", status="unknown")
    open(debug_info_unit, file='debug_info.txt', action="write", &
         position="append", status="unknown")
    write(debug_info_unit,*) 'jacobian ', debug_timestep_count, &
      debug_timestep_cut_count, debug_iteration_count
    close(debug_info_unit)
  endif
#endif

  ! Perturb aux vars
  do ghosted_id = 1, grid%ngmax  ! For each local node do...
    if (patch%imat(ghosted_id) <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    if (.NOT. towg_analytical_derivatives .OR. towg_analytical_derivatives_compare &
        .OR. option%flow%num_as_alyt_derivs) then
      call TOWGAuxVarPerturb(towg%auxvars(:,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             patch%characteristic_curves_array( &
                             patch%sat_func_id(ghosted_id))%ptr, &
                             natural_id,option)
    endif

  enddo

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    call TOWGAccumDerivative(towg%auxvars(:,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option, &
                             Jup,ghosted_id) 
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)

  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_accum'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
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

      imat_up = patch%imat(ghosted_id_up)
      imat_dn = patch%imat(ghosted_id_dn)
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
   
      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)
                              
      call TOWGFluxDerivative(towg%auxvars(:,ghosted_id_up), &
                     global_auxvars(ghosted_id_up), &
                     material_auxvars(ghosted_id_up), &
                     material_parameter%soil_thermal_conductivity(:,imat_up), &
                     towg%auxvars(:,ghosted_id_dn), &
                     global_auxvars(ghosted_id_dn), &
                     material_auxvars(ghosted_id_dn), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     towg_parameter,option,&
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

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_flux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
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
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      imat_dn = patch%imat(ghosted_id)
      if (imat_dn <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icap_dn = patch%sat_func_id(ghosted_id)

      call TOWGBCFluxDerivative(boundary_condition%flow_bc_type, &
                      boundary_condition%flow_aux_mapping, &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      towg%auxvars_bc(sum_connection), &
                      global_auxvars_bc(sum_connection), &
                      towg%auxvars(:,ghosted_id), &
                      global_auxvars(ghosted_id), &
                      material_auxvars(ghosted_id), &
                      material_parameter%soil_thermal_conductivity(:,imat_dn), &
                      cur_connection_set%area(iconn), &
                      cur_connection_set%dist(:,iconn), &
                      towg_parameter,option, &
                      Jdn)

      Jdn = -Jdn
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_bcflux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

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
      call TOWGSrcSinkDerivative(option, &
                        source_sink%flow_condition%towg, &
                        towg%auxvars(:,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        scale,Jup)

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)


    enddo
    source_sink => source_sink%next
  enddo

! Loop over well_data wells if present

  analytical_derivatives = .not. option%flow%numerical_derivatives
  if( analytical_derivatives ) then
    if (WellDataGetFlag()) then
      nflowdof=realization%option%nflowdof
      jerr = 0
      well_data_list => realization%well_data
      well_data => well_data_list%first

      do
        if (.not.associated(well_data)) exit
          call well_data%DoIncrJac(option,nflowdof,Jup,A)
          well_data => well_data%next
      enddo
    endif
  endif

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_srcsink'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! zero out isothermal and inactive cells
  if (towg%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,towg%n_inactive_rows, &
                          towg%inactive_rows_local_ghosted, &
                          qsrc,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif

  if (towg_isothermal) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => towg%row_zeroing_array
    ! zero energy residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        towg_energy_eq_idx - 1 ! zero-based
    enddo
    call MatZeroRowsLocal(A,grid%nlmax,zeros,qsrc,PETSC_NULL_VEC, &
                          PETSC_NULL_VEC,ierr);CHKERRQ(ierr)
  endif

  if (towg_no_oil) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => towg%row_zeroing_array
    ! zero gas component mass balance residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        TOWG_OIL_EQ_IDX - 1 ! zero-based
    enddo
    call MatZeroRowsLocal(A,grid%nlmax,zeros,qsrc,PETSC_NULL_VEC, &
                          PETSC_NULL_VEC,ierr);CHKERRQ(ierr)
  endif

  if (towg_no_gas) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => towg%row_zeroing_array
    ! zero gas component mass balance residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        TOWG_GAS_EQ_IDX - 1 ! zero-based
    enddo
    call MatZeroRowsLocal(A,grid%nlmax,zeros,qsrc,PETSC_NULL_VEC, &
                          PETSC_NULL_VEC,ierr);CHKERRQ(ierr)
  endif
  
  if (realization%debug%matview_Jacobian) then
    string = 'TOWGjacobian'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
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

!  call MatView(J,PETSC_VIEWER_STDOUT_WORLD,ierr)

#ifdef DEBUG_TOWG_FILEOUTPUT
  if (debug_flag > 0) then
    write(word,*) debug_timestep_count
    string = 'jacobian_' // trim(adjustl(word))
    write(word,*) debug_timestep_cut_count
    string = trim(string) // '_' // trim(adjustl(word))
    write(word,*) debug_iteration_count
    string = trim(string) // '_' // trim(adjustl(word)) // '.out'
    call PetscViewerASCIIOpen(realization%option%mycomm,trim(string), &
                              viewer,ierr);CHKERRQ(ierr)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    close(debug_unit)
  endif
#endif

end subroutine TOWGJacobian

! ************************************************************************** !

subroutine TOWGImsTLCheckUpdatePre(line_search,X,dX,changed,realization, &
                                   max_it_before_damping,damping_factor, &
                                   max_pressure_change,ierr)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 12/30/16
  ! 
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Patch_module

  implicit none

  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  type(realization_subsurface_type) :: realization
  PetscInt :: max_it_before_damping
  PetscReal :: damping_factor
  PetscReal :: max_pressure_change

  PetscReal, pointer :: X_p(:), dX_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field

  !type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:)
  !type(global_auxvar_type), pointer :: global_auxvars(:)  

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset

  PetscInt :: pressure_index, saturation_index, temperature_index

  PetscReal :: pressure0, pressure1, del_pressure
  PetscReal :: temperature0, temperature1, del_temperature
  PetscReal :: saturation0, saturation1, del_saturation

  PetscReal, parameter :: max_saturation_change = 0.125d0
  PetscReal, parameter :: max_temperature_change = 10.d0
  PetscReal :: scale, temp_scale, temp_real
  PetscReal, parameter :: tolerance = 0.99d0
  PetscReal, parameter :: initial_scale = 1.d0
  SNES :: snes
  PetscInt :: newton_iteration

  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  patch => realization%patch

  call SNESLineSearchGetSNES(line_search,snes,ierr)
  call SNESGetIterationNumber(snes,newton_iteration,ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  changed = PETSC_TRUE

  !print *, "TOWGImsTLCheckUpdatePre"
  ! truncation
  ! Oil and Gas Saturations must be truncated.  We do not use scaling
  ! here because of the very small values.  just truncation.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    offset = (local_id-1)*option%nflowdof
    saturation_index = offset + TOWG_OIL_SATURATION_DOF
    if ( (X_p(saturation_index) - dX_p(saturation_index)) < 0.d0 ) then
      ! we use 1.d-6 since cancelation can occur with smaller values
      ! this threshold is imposed in the initial condition
      dX_p(saturation_index) = X_p(saturation_index)
    end if
    saturation_index = offset + TOWG_GAS_SATURATION_3PH_DOF
    if ( (X_p(saturation_index) - dX_p(saturation_index)) < 0.d0 ) then
      ! we use 1.d-6 since cancelation can occur with smaller values
      ! this threshold is imposed in the initial condition
      dX_p(saturation_index) = X_p(saturation_index)
    end if
  enddo

  scale = initial_scale
  if (max_it_before_damping > 0 .and. &
      newton_iteration > max_it_before_damping) then
    scale = damping_factor
  endif

#define LIMIT_MAX_PRESSURE_CHANGE
#define LIMIT_MAX_SATURATION_CHANGE
!!#define LIMIT_MAX_TEMPERATURE_CHANGE
!! TRUNCATE_PRESSURE is needed for times when the solve wants
!! to pull them negative.
!#define TRUNCATE_PRESSURE

  ! scaling
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    offset = (local_id-1)*option%nflowdof
    temp_scale = 1.d0
    pressure_index = offset + TOWG_OIL_PRESSURE_DOF
    dX_p(pressure_index) = dX_p(pressure_index) * towg_pressure_scale
    del_pressure = dX_p(pressure_index)
    pressure0 = X_p(pressure_index)
    pressure1 = pressure0 - del_pressure
#ifdef LIMIT_MAX_PRESSURE_CHANGE
    if (dabs(del_pressure) > max_pressure_change) then
      temp_real = dabs(max_pressure_change/del_pressure)
      temp_scale = min(temp_scale,temp_real)
     endif
#endif
#ifdef TRUNCATE_PRESSURE
    if (pressure1 <= 0.d0) then
      if (dabs(del_pressure) > 1.d-40) then
        temp_real = tolerance * dabs(pressure0 / del_pressure)
        temp_scale = min(temp_scale,temp_real)
      endif
    endif
#endif 
!TRUNCATE_PRESSURE
#ifdef LIMIT_MAX_SATURATION_CHANGE
    !oil saturation
    saturation_index = offset + TOWG_OIL_SATURATION_DOF
    del_saturation = dX_p(saturation_index)
    !saturation0 = X_p(saturation_index)
    !saturation1 = saturation0 - del_saturation
    if (dabs(del_saturation) > max_saturation_change) then
       temp_real = dabs(max_saturation_change/del_saturation)
       temp_scale = min(temp_scale,temp_real)
    endif
    !gas saturation
    saturation_index = offset + TOWG_GAS_SATURATION_3PH_DOF
    del_saturation = dX_p(saturation_index)
    !saturation0 = X_p(saturation_index)
    !saturation1 = saturation0 - del_saturation
    if (dabs(del_saturation) > max_saturation_change) then
       temp_real = dabs(max_saturation_change/del_saturation)
       temp_scale = min(temp_scale,temp_real)
    endif
#endif 
!LIMIT_MAX_SATURATION_CHANGE
#ifdef LIMIT_MAX_TEMPERATURE_CHANGE        
    temperature_index  = offset + towg_energy_dof
    del_temperature = dX_p(temperature_index)
    if (dabs(del_temperature) > max_temperature_change) then
       temp_real = dabs(max_temperature_change/del_temperature)
       temp_scale = min(temp_scale,temp_real)
    endif
#endif 
!LIMIT_MAX_TEMPERATURE_CHANGE
    scale = min(scale,temp_scale) 
  enddo

  temp_scale = scale
  call MPI_Allreduce(temp_scale,scale,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_MIN,option%mycomm,ierr)

  ! it performs an homogenous scaling using the smallest scaling factor
  ! over all subdomains domains
  if (scale < 0.9999d0) then
    dX_p = scale*dX_p
  endif

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

end subroutine TOWGImsTLCheckUpdatePre

! ************************************************************************** !

subroutine TOWGBlackOilCheckUpdatePre(line_search,X,dX,changed,realization, &
                                      max_it_before_damping,damping_factor, &
                                      max_pressure_change,ierr)
!------------------------------------------------------------------------------
! Used in TOWG_BLACK_OIL and TOWG_SOLVENT_TL, prepares update for solver
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Oct 2017
!------------------------------------------------------------------------------

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Patch_module

  implicit none

  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  type(realization_subsurface_type) :: realization
  PetscInt :: max_it_before_damping
  PetscReal :: damping_factor
  PetscReal :: max_pressure_change

  PetscReal, pointer :: X_p(:), dX_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset

  PetscInt :: pressure_index, saturation_index, temperature_index

  PetscReal :: pressure0, pressure1, del_pressure
  PetscReal :: temperature0, temperature1, del_temperature
  PetscReal :: saturation0, saturation1, del_saturation, del_sat_cand

  PetscReal,parameter :: max_saturation_change = 0.125d0
  PetscReal,parameter :: max_temperature_change = 10.d0
  PetscReal :: max_pb_change
  PetscReal :: scale, temp_scale, temp_real
  PetscReal, parameter :: tolerance = 0.99d0
  PetscReal, parameter :: initial_scale = 1.d0
  SNES :: snes
  PetscInt :: newton_iteration,istate

  PetscReal :: scand
  PetscBool :: slv_sat_truncate

  type(global_auxvar_type), pointer :: global_auxvars(:)

  slv_sat_truncate = TL4P_slv_sat_truncate

  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  patch => realization%patch

  global_auxvars => patch%aux%Global%auxvars

  call SNESLineSearchGetSNES(line_search,snes,ierr)
  call SNESGetIterationNumber(snes,newton_iteration,ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  max_pb_change=2.0*max_pressure_change
  !max_pb_change= 1.1D5

  ! truncation
  ! Oil saturation must be truncated.  We do not use scaling
  ! here because of the very small values.  just truncation.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    offset = (local_id-1)*option%nflowdof
! Stop oil saturation going negative, not gas - need Sg<0 to get back to Pbub
    saturation_index = offset + TOWG_OIL_SATURATION_DOF
    if ( (X_p(saturation_index) - dX_p(saturation_index)) < 0.d0 ) then
      dX_p(saturation_index) = X_p(saturation_index)
    end if

    ! Stop SOLVENT saturation going negative (if that's even a good idea)
    if( slv_sat_truncate .AND. towg_miscibility_model == TOWG_SOLVENT_TL ) then
      saturation_index = offset + TOWG_SOLV_SATURATION_DOF 
      if ( (X_p(saturation_index) - dX_p(saturation_index)) < 0.d0 ) then
        dX_p(saturation_index) = X_p(saturation_index)
      end if
    endif
  enddo

  scale = initial_scale
  if (max_it_before_damping > 0 .and. &
      newton_iteration > max_it_before_damping) then
    scale = damping_factor
  endif

#define LIMIT_MAX_PRESSURE_CHANGE
#define LIMIT_MAX_SATURATION_CHANGE
!!#define LIMIT_MAX_TEMPERATURE_CHANGE
!! TRUNCATE_PRESSURE is needed for times when the solve wants
!! to pull them negative.
!#define TRUNCATE_PRESSURE

  ! scaling

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    offset = (local_id-1)*option%nflowdof

    istate=global_auxvars(ghosted_id)%istate

    temp_scale = 1.d0
    pressure_index = offset + TOWG_OIL_PRESSURE_DOF
    dX_p(pressure_index) = dX_p(pressure_index) * towg_pressure_scale
    del_pressure = dX_p(pressure_index)
    pressure0 = X_p(pressure_index)
    pressure1 = pressure0 - del_pressure
#ifdef LIMIT_MAX_PRESSURE_CHANGE
    if (dabs(del_pressure) > max_pressure_change) then
      temp_real = dabs(max_pressure_change/del_pressure)
      temp_scale = min(temp_scale,temp_real)
     endif
#endif
#ifdef TRUNCATE_PRESSURE
    if (pressure1 <= 0.d0) then
      if (dabs(del_pressure) > 1.d-40) then
        temp_real = tolerance * dabs(pressure0 / del_pressure)
        temp_scale = min(temp_scale,temp_real)
      endif
    endif
#endif
!TRUNCATE_PRESSURE
#ifdef LIMIT_MAX_SATURATION_CHANGE
    !oil saturation
    saturation_index = offset + TOWG_OIL_SATURATION_DOF
    del_saturation = dX_p(saturation_index)
    !saturation0 = X_p(saturation_index)
    !saturation1 = saturation0 - del_saturation
    if (dabs(del_saturation) > max_saturation_change) then
       temp_real = dabs(max_saturation_change/del_saturation)
       temp_scale = min(temp_scale,temp_real)
    endif
    !gas saturation location
    saturation_index = offset + TOWG_GAS_SATURATION_3PH_DOF
    del_saturation = dX_p(saturation_index)
    !saturation0 = X_p(saturation_index)
    !saturation1 = saturation0 - del_saturation
    if( istate == TOWG_THREE_PHASE_STATE ) then ! Is gas saturation variable
      if (dabs(del_saturation) > max_saturation_change) then
         temp_real = dabs(max_saturation_change/del_saturation)
         temp_scale = min(temp_scale,temp_real)
      endif
    else if( istate == TOWG_LIQ_OIL_STATE ) then ! Is bubble point variable
      if (dabs(del_saturation) > max_pb_change) then
         temp_real = dabs(max_pb_change/del_saturation)
         temp_scale = min(temp_scale,temp_real)
      endif
      ! let's try avoiding negative pb values:
      if ((X_p(saturation_index) - dX_p(saturation_index))  < 0.d0) then
      dX_p(saturation_index) = X_p(saturation_index) - 0.5
      endif
    endif
#endif
!LIMIT_MAX_SATURATION_CHANGE
#ifdef LIMIT_MAX_TEMPERATURE_CHANGE
    temperature_index  = offset + towg_energy_dof
    del_temperature = dX_p(temperature_index)
    if (dabs(del_temperature) > max_temperature_change) then
       temp_real = dabs(max_temperature_change/del_temperature)
       temp_scale = min(temp_scale,temp_real)
    endif
#endif
!LIMIT_MAX_TEMPERATURE_CHANGE
    scale = min(scale,temp_scale)
  enddo

  temp_scale = scale
  call MPI_Allreduce(temp_scale,scale,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_MIN,option%mycomm,ierr)

  ! it performs an homogenous scaling using the smallest scaling factor
  ! over all subdomains domains
  if (scale < 0.9999d0) then
    dX_p = scale*dX_p
  endif

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

end subroutine TOWGBlackOilCheckUpdatePre

! ************************************************************************** !

function TOWGImsTLAverageDensity(sat_up,sat_dn,density_up,density_dn, &
                                 d_denave_den_up,d_denave_den_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/18/16
  ! 

  implicit none

  PetscReal :: sat_up, sat_dn
  PetscReal :: density_up, density_dn

  PetscReal, optional, intent(OUT) :: d_denave_den_up,d_denave_den_dn

  PetscBool :: getDerivs

  PetscReal :: TOWGImsTLAverageDensity

  getDerivs = PETSC_FALSE
  if (present(d_denave_den_up) .AND. present(d_denave_den_dn)) then
    getDerivs = PETSC_TRUE
  endif

!  if ( (towg_miscibility_model == TOWG_IMMISCIBLE) .or. &
!       (towg_miscibility_model == TOWG_TODD_LONGSTAFF) &
!     ) then
!!! EXPERIMENTAL
#if 0
  if ( towg_miscibility_model == TOWG_TODD_LONGSTAFF  .OR. &
       towg_miscibility_model == TOWG_SOLVENT_TL) then
#endif
  if ( towg_miscibility_model == TOWG_TODD_LONGSTAFF  .OR. &
       (towg_miscibility_model == TOWG_SOLVENT_TL .AND. TL4P_altDensity)) then
  !if ( towg_miscibility_model == towg_todd_longstaff ) then
    TOWGImsTLAverageDensity = 0.5d0*(density_up+density_dn)
    if (getDerivs) then
      d_denave_den_up = 0.5d0;d_denave_den_dn = 0.5d0
    endif
  else
    if (sat_up < eps ) then
      TOWGImsTLAverageDensity = density_dn
      if (getDerivs) then
        d_denave_den_up = 0.0d0;d_denave_den_dn = 1.0d0
      endif
    else if (sat_dn < eps ) then 
      TOWGImsTLAverageDensity = density_up
      if (getDerivs) then
        d_denave_den_up = 1.0d0;d_denave_den_dn = 0.0d0
      endif
    else ! in here we could use an armonic average, 
         ! other idea sat weighted average but it needs truncation
      TOWGImsTLAverageDensity = 0.5d0*(density_up+density_dn)
      if (getDerivs) then 
        d_denave_den_up = 0.5d0;d_denave_den_dn = 0.5d0
      endif
    end if
  end if

end function TOWGImsTLAverageDensity

! ************************************************************************** !

subroutine TOWGDestroy(realization)
  ! 
  ! Deallocates variables associated with TOWG
  ! 
  ! Author: Paolo Orsini
  ! Date: 01/19/17
  ! 

  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization
  
  ! place anything that needs to be freed here.
#ifdef GLOBALWORKERS
  if (towg_analytical_derivatives) then
    deallocate(D_den_kg_ave_up)
    deallocate(D_den_kg_ave_dn)
    deallocate(D_den_ave_up)
    deallocate(D_den_ave_dn)
    deallocate(D_delta_presure_up)
    deallocate(D_delta_presure_dn)
    deallocate( D_mobility_up)
    deallocate(D_mobility_dn)
    deallocate(D_uH_up)
    deallocate(D_uH_dn)
    deallocate(D_v_darcy_up)
    deallocate(D_v_darcy_dn)
    deallocate(D_q_up)
    deallocate(D_q_dn)
    deallocate(D_mole_flux_up)
    deallocate(D_mole_flux_dn)
    deallocate(D_xmf_up)
    deallocate(D_xmf_dn)

    deallocate(D_sat_liquid_up)
    deallocate(D_sat_liquid_dn)
    deallocate(D_k_eff_up)
    deallocate(D_k_eff_dn)
    deallocate(D_k_eff_ave_up)
    deallocate(D_k_eff_ave_dn)
    deallocate(D_delta_temp_up)
    deallocate(D_delta_temp_dn)
    deallocate(D_worker1)
    deallocate(D_worker2)
  endif
#endif


end subroutine TOWGDestroy

! ************************************************************************** !

end module TOWG_module
