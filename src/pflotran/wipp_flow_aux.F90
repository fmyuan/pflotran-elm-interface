module WIPP_Flow_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#if 0
  PetscInt, public, parameter :: NULL_STATE = 0
  PetscInt, public, parameter :: LIQUID_STATE = 1
  PetscInt, public, parameter :: GAS_STATE = 2
  PetscInt, public, parameter :: TWO_PHASE_STATE = 3
  PetscInt, public, parameter :: ANY_STATE = 4
#endif
  
  PetscBool, public :: wippflo_analytical_derivatives = PETSC_FALSE
  PetscBool, public :: wippflo_fix_upwind_direction = PETSC_TRUE
  PetscBool, public :: wippflo_update_upwind_direction = PETSC_FALSE
  PetscBool, public :: wippflo_count_upwind_dir_flip = PETSC_FALSE
  PetscBool, public :: wippflo_use_fracture = PETSC_TRUE
  PetscBool, public :: wippflo_use_creep_closure = PETSC_TRUE
  !TODO(geh): hardwire gas to H2
  PetscReal, public :: fmw_comp(2) = [FMWH2O,FMWAIR]
  PetscReal, public :: wippflo_max_pressure_change = 5.d4
  PetscInt, public :: wippflo_max_it_before_damping = UNINITIALIZED_INTEGER
  PetscReal, public :: wippflo_damping_factor = 0.6d0
  

  PetscInt, parameter, public :: WIPPFLO_LIQUID_PRESSURE_DOF = 1
  PetscInt, parameter, public :: WIPPFLO_GAS_SATURATION_DOF = 2
  
  PetscInt, parameter, public :: WIPPFLO_LIQUID_EQUATION_INDEX = 1
  PetscInt, parameter, public :: WIPPFLO_GAS_EQUATION_INDEX = 2
  
  PetscInt, parameter, public :: WIPPFLO_STATE_INDEX = 1
  PetscInt, parameter, public :: WIPPFLO_LIQUID_PRESSURE_INDEX = 1
  PetscInt, parameter, public :: WIPPFLO_GAS_SATURATION_INDEX = 2
  PetscInt, parameter, public :: WIPPFLO_LIQUID_FLUX_INDEX = 3
  PetscInt, parameter, public :: WIPPFLO_GAS_FLUX_INDEX = 4
  PetscInt, parameter, public :: WIPPFLO_LIQUID_CONDUCTANCE_INDEX = 5
  PetscInt, parameter, public :: WIPPFLO_GAS_CONDUCTANCE_INDEX = 6
  PetscInt, parameter, public :: WIPPFLO_MAX_INDEX = 6
  
  PetscInt, parameter, public :: WIPPFLO_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: WIPPFLO_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: WIPPFLO_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: WIPPFLO_UPDATE_FOR_BOUNDARY = 2
  
  PetscReal, parameter, public :: WIPPFLO_PRESSURE_SCALE = 1.d0

  ! these variables, which are global to general, can be modified
  PetscInt, public :: dof_to_primary_variable(2)
  
  type, public :: wippflo_auxvar_type
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal, pointer :: den(:)    ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: den_kg(:) ! (iphase) kg/m^3 phase
    !geh: leave xmol in object until analytical derivatives have been fixed.
    PetscReal :: xmol(2,2)
    PetscReal :: temp
    PetscReal, pointer :: mobility(:) ! relative perm / kinematic viscosity
    PetscReal :: effective_porosity ! factors in compressibility
    PetscReal :: pert
    PetscReal :: fracture_perm_scaling_factor
    PetscReal :: klinkenberg_scaling_factor(3)
    type(wippflo_derivative_auxvar_type), pointer :: d
  end type wippflo_auxvar_type
  
  type, public :: wippflo_derivative_auxvar_type
    PetscReal :: pc_satg
    PetscReal :: por_p
    PetscReal :: denl_pl
    PetscReal :: deng_pg
    PetscReal :: deng_pa
    PetscReal :: dengkg_pg
    PetscReal :: denv
    PetscReal :: dena
    PetscReal :: denv_pg
    PetscReal :: dena_pg
    
    PetscReal :: psat_p
    PetscReal :: pv_p
    PetscReal :: pv_pa
    PetscReal :: mobilityl_pl
    PetscReal :: mobilityl_satg
    PetscReal :: mobilityg_pg
    PetscReal :: mobilityg_satg
    PetscReal :: mobilityg_pa
    PetscReal :: mug
    PetscReal :: mug_pg
  end type wippflo_derivative_auxvar_type
  
  type, public :: wippflo_parameter_type
    PetscBool :: check_post_converged
    PetscBool :: fix_upwind_direction
  end type wippflo_parameter_type
  
  type, public :: wippflo_type
    PetscInt :: n_inactive_rows
    PetscInt, pointer :: inactive_rows_local(:), inactive_rows_local_ghosted(:)
    PetscInt, pointer :: row_zeroing_array(:)

    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    PetscInt, pointer :: upwind_direction(:,:)
    PetscInt, pointer :: upwind_direction_bc(:,:)
    type(wippflo_parameter_type), pointer :: wippflo_parameter
    type(wippflo_auxvar_type), pointer :: auxvars(:,:)
    type(wippflo_auxvar_type), pointer :: auxvars_bc(:)
    type(wippflo_auxvar_type), pointer :: auxvars_ss(:)
  end type wippflo_type

  interface WIPPFloAuxVarDestroy
    module procedure WIPPFloAuxVarSingleDestroy
    module procedure WIPPFloAuxVarArray1Destroy
    module procedure WIPPFloAuxVarArray2Destroy
  end interface WIPPFloAuxVarDestroy
  
  interface WIPPFloOutputAuxVars
    module procedure WIPPFloOutputAuxVars1
    module procedure WIPPFloOutputAuxVars2
  end interface WIPPFloOutputAuxVars
  
  public :: WIPPFloAuxCreate, &
            WIPPFloAuxDestroy, &
            WIPPFloAuxVarCompute, &
            WIPPFloAuxVarInit, &
            WIPPFloAuxVarCopy, &
            WIPPFloScalePerm, &
            WIPPFloAuxVarDestroy, &
            WIPPFloAuxVarStrip, &
            WIPPFloAuxVarPerturb, &
            WIPPFloPrintAuxVars, &
            WIPPFloOutputAuxVars

contains

! ************************************************************************** !

function WIPPFloAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Option_module

  implicit none

  type(option_type) :: option
    
  type(wippflo_type), pointer :: WIPPFloAuxCreate
  
  type(wippflo_type), pointer :: aux

  dof_to_primary_variable(1:2) = &
    reshape([WIPPFLO_LIQUID_PRESSURE_INDEX, WIPPFLO_GAS_SATURATION_INDEX], &
             shape(dof_to_primary_variable))
  
  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%upwind_direction)
  nullify(aux%upwind_direction_bc)
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  aux%n_inactive_rows = 0
  nullify(aux%inactive_rows_local)
  nullify(aux%inactive_rows_local_ghosted)
  nullify(aux%row_zeroing_array)

  allocate(aux%wippflo_parameter)
  aux%wippflo_parameter%check_post_converged = PETSC_FALSE
  aux%wippflo_parameter%fix_upwind_direction = PETSC_TRUE
  
  WIPPFloAuxCreate => aux
  
end function WIPPFloAuxCreate

! ************************************************************************** !

subroutine WIPPFloAuxVarInit(auxvar,allocate_derivative,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Option_module

  implicit none
  
  type(wippflo_auxvar_type) :: auxvar
  PetscBool :: allocate_derivative
  type(option_type) :: option

  auxvar%temp = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%pert = 0.d0
  auxvar%fracture_perm_scaling_factor = 1.d0
  auxvar%klinkenberg_scaling_factor = 1.d0
  
  auxvar%xmol = 0.d0
  auxvar%xmol(1,1) = 1.d0
  auxvar%xmol(2,2) = 1.d0
  allocate(auxvar%pres(option%nphase+FOUR_INTEGER))
  auxvar%pres = 0.d0
  allocate(auxvar%sat(option%nphase))
  auxvar%sat = 0.d0
  allocate(auxvar%den(option%nphase))
  auxvar%den = 0.d0
  allocate(auxvar%den_kg(option%nphase))
  auxvar%den_kg = 0.d0
  allocate(auxvar%mobility(option%nphase))
  auxvar%mobility = 0.d0
  if (allocate_derivative) then
    allocate(auxvar%d)
    auxvar%d%pc_satg = 0.d0
    auxvar%d%por_p = 0.d0
    auxvar%d%denl_pl = 0.d0
    auxvar%d%deng_pg = 0.d0
    auxvar%d%deng_pa = 0.d0
    auxvar%d%dengkg_pg = 0.d0
    auxvar%d%denv = 0.d0
    auxvar%d%dena = 0.d0
    auxvar%d%denv_pg = 0.d0
    auxvar%d%dena_pg = 0.d0
        
    auxvar%d%psat_p = 0.d0
    auxvar%d%pv_p = 0.d0
    auxvar%d%pv_pa = 0.d0
    auxvar%d%mug = 0.d0
    auxvar%d%mug_pg = 0.d0
    auxvar%d%mobilityl_pl = 0.d0
    auxvar%d%mobilityl_satg = 0.d0
    auxvar%d%mobilityg_pg = 0.d0
    auxvar%d%mobilityg_satg = 0.d0
    auxvar%d%mobilityg_pa = 0.d0
  else
    nullify(auxvar%d)
  endif
  
end subroutine WIPPFloAuxVarInit

! ************************************************************************** !

subroutine WIPPFloAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Option_module

  implicit none
  
  type(wippflo_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%mobility = auxvar%mobility
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%pert = auxvar%pert

end subroutine WIPPFloAuxVarCopy

! ************************************************************************** !

subroutine WIPPFloAuxVarCompute(x,wippflo_auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use Creep_Closure_module
  use Fracture_module
  use Klinkenberg_module
  use WIPP_module
  use Variables_module, only : SOIL_REFERENCE_PRESSURE
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(wippflo_auxvar_type) :: wippflo_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(creep_closure_type), pointer :: creep_closure
  PetscInt :: natural_id

  PetscInt :: gid, lid
  PetscReal :: cell_pressure, water_vapor_pressure
  PetscReal :: den_water_vapor, den_kg_water_vapor
  PetscReal :: den_air
  PetscReal :: krl, visl, dvis_dp, dvis_dT, dvis_dpa
  PetscReal :: dkrl_dsatl, dkrl_dsatg
  PetscReal :: dkrg_dsatl, dkrg_dsatg
  PetscReal :: krg, visg
  PetscReal :: guess, dummy
  PetscInt :: cpid, spid
  PetscReal :: creep_closure_time
  PetscReal :: aux(1)
  PetscReal :: dpor_dp
  PetscReal :: tempreal
  PetscReal :: dpair_dpgas
  PetscReal :: dden_air_dpa, dden_air_dpg
  PetscReal :: dden_water_vapor_dpv
  PetscReal :: dpc_dsatl
  PetscReal :: perm_for_cc
  PetscReal :: prev_effective_porosity
  PetscErrorCode :: ierr

  ! from init.F90
!  option%nphase = 2
!  option%liquid_phase = 1  ! liquid_pressure
!  option%gas_phase = 2     ! gas_pressure

!  option%capillary_pressure_id = 3
!  option%saturation_pressure_id = 4

  lid = option%liquid_phase
  gid = option%gas_phase
  cpid = option%capillary_pressure_id
  spid = option%saturation_pressure_id

  ! Two Phase State Variables
  wippflo_auxvar%temp = option%reference_temperature
  wippflo_auxvar%pres(lid) = x(WIPPFLO_LIQUID_PRESSURE_DOF)
  wippflo_auxvar%sat(gid) = x(WIPPFLO_GAS_SATURATION_DOF)
  ! calculate saturation pressure as reference.
  call EOSWaterSaturationPressure(wippflo_auxvar%temp, &
                                  wippflo_auxvar%pres(spid),ierr)
  if (associated(wippflo_auxvar%d)) then
    dpair_dpgas = 1.d0
  endif
  wippflo_auxvar%sat(lid) = 1.d0 - wippflo_auxvar%sat(gid)

  cell_pressure = wippflo_auxvar%pres(lid)

  prev_effective_porosity = wippflo_auxvar%effective_porosity

  ! calculate effective porosity as a function of pressure
  if (option%iflag /= WIPPFLO_UPDATE_FOR_BOUNDARY) then
    dpor_dp = 0.d0
    wippflo_auxvar%effective_porosity = material_auxvar%porosity_base
    ! creep_closure, fracture, and soil_compressibility are mutually exclusive
    if (option%flow%creep_closure_on .and. wippflo_use_creep_closure) then
      creep_closure => wipp%creep_closure_tables_array( &
                         material_auxvar%creep_closure_id )%ptr
      if (associated(creep_closure)) then
        if (option%time > creep_closure%time_datamax .OR. & 
            option%time > creep_closure%time_closeoff) then
          material_auxvar%porosity_base = prev_effective_porosity
          call MaterialAuxVarSetValue(material_auxvar,SOIL_REFERENCE_PRESSURE, &
                                      cell_pressure)
          ! index 1 of wipp%creep_closure_tables_array is a null pointer
          material_auxvar%creep_closure_id = 1 
          nullify(creep_closure)
        else if (cell_pressure > creep_closure%shutdown_pressure) then
          ! fix to shutdown pressure and porosity at shutdown pressure
          wippflo_auxvar%effective_porosity = &
           creep_closure%Evaluate(option%time,creep_closure%shutdown_pressure)
          material_auxvar%porosity_base = wippflo_auxvar%effective_porosity
          call MaterialAuxVarSetValue(material_auxvar,SOIL_REFERENCE_PRESSURE, &
                                      creep_closure%shutdown_pressure)
          ! index 1 of wipp%creep_closure_tables_array is a null pointer
          material_auxvar%creep_closure_id = 1 
          nullify(creep_closure)
        else
          ! option%time here is the t time, not t + dt time.
          creep_closure_time = option%time
          if (option%iflag /= WIPPFLO_UPDATE_FOR_FIXED_ACCUM) then
            creep_closure_time = creep_closure_time + option%flow_dt
          endif
          
          wippflo_auxvar%effective_porosity = &
            creep_closure%Evaluate(creep_closure_time,cell_pressure)
          wippflo_auxvar%effective_porosity = & 
            max(wippflo_auxvar%effective_porosity,creep_closure%porosity_minimum)
        endif
      else if (associated(material_auxvar%fracture) .and. &
               wippflo_use_fracture) then
          call FracturePoroEvaluate(material_auxvar,cell_pressure, &
                                    wippflo_auxvar%effective_porosity,dpor_dp)
      else if (soil_compressibility_index > 0) then
          call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                    wippflo_auxvar%effective_porosity,dpor_dp)
      endif
    else if (associated(material_auxvar%fracture) .and. &
             wippflo_use_fracture) then
      call FracturePoroEvaluate(material_auxvar,cell_pressure, &
                                wippflo_auxvar%effective_porosity,dpor_dp)
    else if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                wippflo_auxvar%effective_porosity,dpor_dp)
    endif
    if (option%iflag /= WIPPFLO_UPDATE_FOR_DERIVATIVE) then
      ! this needs to be set for proper output
      material_auxvar%porosity = wippflo_auxvar%effective_porosity
    endif
  endif
  if (associated(wippflo_auxvar%d)) then
    wippflo_auxvar%d%por_p = dpor_dp
  endif

  wippflo_auxvar%fracture_perm_scaling_factor = 1.d0
  if (associated(material_auxvar%fracture) .and. &
      wippflo_use_fracture) then
    call FracturePermScale(material_auxvar,cell_pressure, &
                           wippflo_auxvar%effective_porosity, &
                           wippflo_auxvar%fracture_perm_scaling_factor)
  endif
  ! According to the order of operations (PTHRESH/RELPERM prior to ROCKCOMP) 
  ! in PROPS1 in BRAGFLO, fracture has no impact on PTHRESH perm. Thus, the
  ! permeability used in characteristic curves is unmodified.
  perm_for_cc = material_auxvar%permeability(perm_xx_index)
  select type(sf => characteristic_curves%saturation_function)
    class is(sat_func_WIPP_type)
      sf%pct = sf%pct_a * perm_for_cc ** sf%pct_exp
      option%pct_updated = PETSC_TRUE
    class default
      option%pct_updated = PETSC_FALSE
  end select
  call characteristic_curves%saturation_function% &
          CapillaryPressure(wippflo_auxvar%sat(lid),wippflo_auxvar%pres(cpid), &
                            dpc_dsatl,option)                             
  if (associated(wippflo_auxvar%d)) then
    wippflo_auxvar%d%pc_satg = -1.d0*dpc_dsatl
  endif
 
  wippflo_auxvar%pres(gid) = wippflo_auxvar%pres(lid) + wippflo_auxvar%pres(cpid)

  ! Klinkenberg effect is based on gas pressure. Therefore, it cannot be
  ! updated prior to this location and should be skipped if a fracture cell.
  wippflo_auxvar%klinkenberg_scaling_factor = 1.d0
  if (associated(klinkenberg)) then
    if (associated(material_auxvar%fracture) .and. wippflo_use_fracture) then
      if (.not.material_auxvar%fracture%fracture_is_on) then
        call klinkenberg%Scale(material_auxvar%permeability, &
                               wippflo_auxvar%pres(gid), &
                               wippflo_auxvar%klinkenberg_scaling_factor)
      endif
    else  
      call klinkenberg%Scale(material_auxvar%permeability, &
                             wippflo_auxvar%pres(gid), &
                             wippflo_auxvar%klinkenberg_scaling_factor)
    endif
  endif

  ! ALWAYS UPDATE THERMODYNAMIC PROPERTIES FOR BOTH PHASES!!!

  ! Liquid phase thermodynamic properties
  if (.not.option%flow%density_depends_on_salinity) then
    if (associated(wippflo_auxvar%d)) then
      call EOSWaterDensity(wippflo_auxvar%temp,wippflo_auxvar%pres(lid), &
                           wippflo_auxvar%den_kg(lid),wippflo_auxvar%den(lid), &
                           wippflo_auxvar%d%denl_pl,dummy,ierr)
    else
      call EOSWaterDensity(wippflo_auxvar%temp,wippflo_auxvar%pres(lid), &
                           wippflo_auxvar%den_kg(lid),wippflo_auxvar%den(lid),ierr)
    endif
  else
    aux(1) = global_auxvar%m_nacl(1)
    if (associated(wippflo_auxvar%d)) then
      call EOSWaterDensityExt(wippflo_auxvar%temp,wippflo_auxvar%pres(lid),aux, &
                              wippflo_auxvar%den_kg(lid),wippflo_auxvar%den(lid), &
                              wippflo_auxvar%d%denl_pl,dummy,ierr)
    else
      call EOSWaterDensityExt(wippflo_auxvar%temp,wippflo_auxvar%pres(lid),aux, &
                              wippflo_auxvar%den_kg(lid),wippflo_auxvar%den(lid),ierr)
    endif
  endif

  ! Gas phase thermodynamic properties
  if (associated(wippflo_auxvar%d)) then
    call EOSGasDensity(wippflo_auxvar%temp,wippflo_auxvar%pres(gid),den_air, &
                       dummy,dden_air_dpa,ierr)
    ! add in partial w/respec to pa_T
    dden_air_dpg = dden_air_dpa*dpair_dpgas
  else
    call EOSGasDensity(wippflo_auxvar%temp,wippflo_auxvar%pres(gid),den_air,ierr)
  endif
  den_water_vapor = 0.d0
  den_kg_water_vapor = 0.d0
  wippflo_auxvar%den(gid) = den_water_vapor + den_air
  wippflo_auxvar%den_kg(gid) = den_kg_water_vapor + den_air*fmw_comp(gid)
  den_water_vapor = 0.d0
  
  if (associated(wippflo_auxvar%d)) then
    dden_water_vapor_dpv = 0.d0
    wippflo_auxvar%d%pv_p = 0.d0
    wippflo_auxvar%d%pv_pa = 0.d0
    wippflo_auxvar%d%denv = den_water_vapor
    wippflo_auxvar%d%dena = den_air
    wippflo_auxvar%d%denv_pg = dden_water_vapor_dpv * wippflo_auxvar%d%pv_p 
    wippflo_auxvar%d%dena_pg = dden_air_dpa * dpair_dpgas
      
    wippflo_auxvar%d%deng_pg = dden_water_vapor_dpv * wippflo_auxvar%d%pv_p + &
                            dden_air_dpa * dpair_dpgas
    wippflo_auxvar%d%deng_pa = dden_water_vapor_dpv * wippflo_auxvar%d%pv_pa + &
                            dden_air_dpa
    wippflo_auxvar%d%dengkg_pg = dden_water_vapor_dpv * wippflo_auxvar%d%pv_p * FMWH2O + &
                              dden_air_dpa * dpair_dpgas * fmw_comp(2)
  endif
  
  ! Liquid Phase
  call characteristic_curves%liq_rel_perm_function% &
          RelativePermeability(wippflo_auxvar%sat(lid),krl,dkrl_dsatl,option)
  ! dkrl_sat is with respect to liquid pressure, but the primary dependent
  ! variable is gas pressure.  therefore, negate
  dkrl_dsatg = -1.d0 * dkrl_dsatl
  if (.not.option%flow%density_depends_on_salinity) then
    if (associated(wippflo_auxvar%d)) then
      call EOSWaterViscosity(wippflo_auxvar%temp,wippflo_auxvar%pres(lid), &
                              wippflo_auxvar%pres(spid), &
                              0.d0,visl, &
                              dummy,dvis_dp,ierr)
    else
      call EOSWaterViscosity(wippflo_auxvar%temp,wippflo_auxvar%pres(lid), &
                              wippflo_auxvar%pres(spid),visl,ierr)
    endif
  else
    aux(1) = global_auxvar%m_nacl(1)
    if (associated(wippflo_auxvar%d)) then
      call EOSWaterViscosityExt(wippflo_auxvar%temp,wippflo_auxvar%pres(lid), &
                                wippflo_auxvar%pres(spid), &
                                0.d0,aux,visl, &
                                dummy,dvis_dp,ierr)
    else
      call EOSWaterViscosityExt(wippflo_auxvar%temp,wippflo_auxvar%pres(lid), &
                                wippflo_auxvar%pres(spid),aux,visl,ierr)
    endif
  endif
  wippflo_auxvar%mobility(lid) = krl/visl
  if (associated(wippflo_auxvar%d)) then
    ! use chainrule for derivative
    tempreal = -1.d0*krl/(visl*visl)
    wippflo_auxvar%d%mobilityl_pl = tempreal*dvis_dp
    wippflo_auxvar%d%mobilityl_satg = dkrl_dsatg/visl
  endif

  ! Gas Phase
  call characteristic_curves%gas_rel_perm_function% &
          RelativePermeability(wippflo_auxvar%sat(lid),krg,dkrg_dsatl,option)                            
  ! dkrl_sat is with respect to liquid pressure, but the primary dependent
  ! variable is gas pressure.  therefore, negate
  dkrg_dsatg = -1.d0 * dkrg_dsatl
  ! STOMP uses separate functions for calculating viscosity of vapor and
  ! and air (WATGSV,AIRGSV) and then uses GASVIS to calculate mixture 
  ! viscosity.
  if (associated(wippflo_auxvar%d)) then
    call EOSGasViscosity(wippflo_auxvar%temp, wippflo_auxvar%pres(gid), &
                          wippflo_auxvar%pres(gid), den_air, &
  ! viscosity routine says it should be derivative of rho_comp wrt pg
!geh                           dden_air_dT, wippflo_auxvar%d%deng_pg, &
                          0.d0, dden_air_dpa, dden_air_dpg, &
                          0.d0, dpair_dpgas, &
                          visg,dummy,dvis_dpa,dvis_dp,ierr)      
  else    
    call EOSGasViscosity(wippflo_auxvar%temp,wippflo_auxvar%pres(gid), &
                          wippflo_auxvar%pres(gid),den_air,visg,ierr)
  endif
  wippflo_auxvar%mobility(gid) = krg/visg
  if (associated(wippflo_auxvar%d)) then
    ! use chainrule for derivative
    tempreal = -1.d0*krg/(visg*visg)
    wippflo_auxvar%d%mobilityg_pg = tempreal*dvis_dp
    wippflo_auxvar%d%mobilityg_satg = dkrg_dsatg/visg
    wippflo_auxvar%d%mobilityg_pa = tempreal*dvis_dpa
    wippflo_auxvar%d%mug = visg
    wippflo_auxvar%d%mug_pg = dvis_dp
  endif    

end subroutine WIPPFloAuxVarCompute

! ************************************************************************** !

subroutine WIPPFloAuxVarPerturb(wippflo_auxvar,global_auxvar, &
                                material_auxvar, &
                                characteristic_curves,natural_id, &
                                option)
  ! 
  ! Calculates auxiliary variables for perturbed system
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(wippflo_auxvar_type) :: wippflo_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
     
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), &
               pert(option%nflowdof), x_pert_save(option%nflowdof)

  PetscReal :: tempreal
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-10
  PetscInt :: idof

  x(WIPPFLO_LIQUID_PRESSURE_DOF) = &
    wippflo_auxvar(ZERO_INTEGER)%pres(option%liquid_phase)
  x(WIPPFLO_GAS_SATURATION_DOF) = &
    wippflo_auxvar(ZERO_INTEGER)%sat(option%gas_phase)
  pert(WIPPFLO_LIQUID_PRESSURE_DOF) = &
    perturbation_tolerance*x(WIPPFLO_LIQUID_PRESSURE_DOF)+min_perturbation
  if (x(WIPPFLO_GAS_SATURATION_DOF) > 0.5d0) then 
    pert(WIPPFLO_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
  else
    pert(WIPPFLO_GAS_SATURATION_DOF) = perturbation_tolerance
  endif
  
  ! WIPPFLO_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = WIPPFLO_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    wippflo_auxvar(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    x_pert_save = x_pert
    call WIPPFloAuxVarCompute(x_pert,wippflo_auxvar(idof),global_auxvar, &
                              material_auxvar, &
                              characteristic_curves,natural_id,option)
  enddo

  wippflo_auxvar(WIPPFLO_LIQUID_PRESSURE_DOF)%pert = &
    wippflo_auxvar(WIPPFLO_LIQUID_PRESSURE_DOF)%pert / WIPPFLO_PRESSURE_SCALE
  
end subroutine WIPPFloAuxVarPerturb

! ************************************************************************** !

subroutine WIPPFloScalePerm(wippflo_auxvar,material_auxvar,perm,ivar)
  ! 
  ! Scales the permeability by fracture and Klinkenbert (gas only)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/29/17
  ! 
  use Material_Aux_class
  use Variables_module, only : GAS_PERMEABILITY, GAS_PERMEABILITY_X, &
                               GAS_PERMEABILITY_Y, GAS_PERMEABILITY_Z

  implicit none

  type(wippflo_auxvar_type) :: wippflo_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: perm
  PetscInt :: ivar

  if (associated(material_auxvar%fracture)) then
    if (material_auxvar%fracture%fracture_is_on) then
      perm = perm * wippflo_auxvar%fracture_perm_scaling_factor
    endif
  endif
  select case(ivar)
    case(GAS_PERMEABILITY,GAS_PERMEABILITY_X)
      perm = perm * wippflo_auxvar%klinkenberg_scaling_factor(X_DIRECTION)
    case(GAS_PERMEABILITY_Y)
      perm = perm * wippflo_auxvar%klinkenberg_scaling_factor(Y_DIRECTION)
    case(GAS_PERMEABILITY_Z)
      perm = perm * wippflo_auxvar%klinkenberg_scaling_factor(Z_DIRECTION)
  end select

end subroutine WIPPFloScalePerm

! ************************************************************************** !

subroutine WIPPFloPrintAuxVars(wippflo_auxvar,global_auxvar,material_auxvar, &
                               natural_id,string,option)
  ! 
  ! Prints out the contents of an auxvar
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(wippflo_auxvar_type) :: wippflo_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  PetscInt :: cpid, spid
  PetscInt :: gid, lid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_saturation, gas_saturation

  lid = option%liquid_phase
  gid = option%gas_phase
  cpid = option%capillary_pressure_id
  spid = option%saturation_pressure_id

  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  liquid_density = wippflo_auxvar%den(lid)
  gas_density = wippflo_auxvar%den(gid)
  liquid_saturation = wippflo_auxvar%sat(lid)
  gas_saturation = wippflo_auxvar%sat(gid)
  liquid_mass = (liquid_density*liquid_saturation)* & 
                wippflo_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (gas_density*gas_saturation)* & 
              wippflo_auxvar%effective_porosity*material_auxvar%volume
  print *, 'tot liq comp mass [kmol]: ', liquid_mass
  print *, 'tot gas comp mass [kmol]: ', gas_mass
  print *, '         liquid pressure: ', wippflo_auxvar%pres(lid)
  print *, '            gas pressure: ', wippflo_auxvar%pres(gid)
  print *, '      capillary pressure: ', wippflo_auxvar%pres(cpid)
  print *, '     saturation pressure: ', wippflo_auxvar%pres(spid)
  print *, '       liquid saturation: ', wippflo_auxvar%sat(lid)
  print *, '          gas saturation: ', wippflo_auxvar%sat(gid)
  print *, '   liquid density [kmol]: ', wippflo_auxvar%den(lid)
  print *, '      gas density [kmol]: ', wippflo_auxvar%den(gid)
  print *, '     liquid density [kg]: ', wippflo_auxvar%den_kg(lid)
  print *, '        gas density [kg]: ', wippflo_auxvar%den_kg(gid)
  print *, '         temperature [C]: ', wippflo_auxvar%temp
  print *, '         liquid mobility: ', wippflo_auxvar%mobility(lid)
  print *, '            gas mobility: ', wippflo_auxvar%mobility(gid)
  print *, '      effective porosity: ', wippflo_auxvar%effective_porosity
  print *, '--------------------------------------------------------'

end subroutine WIPPFloPrintAuxVars

! ************************************************************************** !

subroutine WIPPFloOutputAuxVars1(wippflo_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,string,append,option)
  ! 
  ! Prints out the contents of an auxvar to a file
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(wippflo_auxvar_type) :: wippflo_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  PetscBool :: append
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string2
  PetscInt :: cpid, spid
  PetscInt :: gid, lid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_saturation, gas_saturation

  lid = option%liquid_phase
  gid = option%gas_phase
  cpid = option%capillary_pressure_id
  spid = option%saturation_pressure_id

  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0

  write(string2,*) natural_id
  string2 = trim(adjustl(string)) // '_' // trim(adjustl(string2)) // '.txt'
  if (append) then
    open(unit=86,file=string2,position='append')
  else
    open(unit=86,file=string2)
  endif

  write(86,*) '--------------------------------------------------------'
  write(86,*) trim(string)
  write(86,*) '             cell id: ', natural_id
  liquid_density = wippflo_auxvar%den(lid)
  gas_density = wippflo_auxvar%den(gid)
  liquid_saturation = wippflo_auxvar%sat(lid)
  gas_saturation = wippflo_auxvar%sat(gid)
  liquid_mass = (liquid_density*liquid_saturation)* & 
                 wippflo_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (gas_density*gas_saturation)* & 
              wippflo_auxvar%effective_porosity*material_auxvar%volume
  write(86,*) 'tot liq comp mass [kmol]: ', liquid_mass
  write(86,*) 'tot gas comp mass [kmol]: ', gas_mass
  write(86,*) '         liquid pressure: ', wippflo_auxvar%pres(lid)
  write(86,*) '            gas pressure: ', wippflo_auxvar%pres(gid)
  write(86,*) '      capillary pressure: ', wippflo_auxvar%pres(cpid)
  write(86,*) '     saturation pressure: ', wippflo_auxvar%pres(spid)
  write(86,*) '         temperature [C]: ', wippflo_auxvar%temp
  write(86,*) '       liquid saturation: ', wippflo_auxvar%sat(lid)
  write(86,*) '          gas saturation: ', wippflo_auxvar%sat(gid)
  write(86,*) '   liquid density [kmol]: ', wippflo_auxvar%den(lid)
  write(86,*) '     liquid density [kg]: ', wippflo_auxvar%den_kg(lid)
  write(86,*) '      gas density [kmol]: ', wippflo_auxvar%den(gid)
  write(86,*) '        gas density [kg]: ', wippflo_auxvar%den_kg(gid)
  write(86,*) '         liquid mobility: ', wippflo_auxvar%mobility(lid)
  write(86,*) '            gas mobility: ', wippflo_auxvar%mobility(gid)
  write(86,*) '      effective porosity: ', wippflo_auxvar%effective_porosity
  write(86,*) '...'
  write(86,*) liquid_mass
  write(86,*) gas_mass
  write(86,*) wippflo_auxvar%pres(lid)
  write(86,*) wippflo_auxvar%pres(gid)
  write(86,*) wippflo_auxvar%pres(cpid)
  write(86,*) wippflo_auxvar%pres(spid)
  write(86,*) wippflo_auxvar%temp
  write(86,*) wippflo_auxvar%sat(lid)
  write(86,*) wippflo_auxvar%sat(gid)
  write(86,*) wippflo_auxvar%den(lid)
  write(86,*) wippflo_auxvar%den_kg(lid)
  write(86,*) wippflo_auxvar%den(gid)
  write(86,*) wippflo_auxvar%den_kg(gid)
  write(86,*) ''
  write(86,*) wippflo_auxvar%mobility(lid)
  write(86,*) wippflo_auxvar%mobility(gid)
  write(86,*) wippflo_auxvar%effective_porosity
  write(86,*) '--------------------------------------------------------'
  
  close(86)

end subroutine WIPPFloOutputAuxVars1

! ************************************************************************** !

subroutine WIPPFloOutputAuxVars2(wippflo_auxvars,global_auxvars,option)
  ! 
  ! Prints out the contents of an auxvar to a file
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Global_Aux_module
  use Option_module

  implicit none

  type(wippflo_auxvar_type) :: wippflo_auxvars(0:,:)
  type(global_auxvar_type) :: global_auxvars(:)
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: cpid, spid
  PetscInt :: gid, lid
  PetscInt :: i, n, idof

  lid = option%liquid_phase
  gid = option%gas_phase
  cpid = option%capillary_pressure_id

  string = 'wippflo_auxvar.txt'
  open(unit=86,file=string)
  
  n = size(global_auxvars)

100 format(a,100('','',i9))
  
  write(86,'(a,100('','',i9))') '             cell id: ', &
    ((i,i=1,n),idof=0,2)
  write(86,'(a,100('','',i2))') '                idof: ', &
    ((idof,i=1,n),idof=0,2)
  write(86,100) '      liquid pressure: ', &
    ((wippflo_auxvars(idof,i)%pres(lid),i=1,n),idof=0,2)
  write(86,100) '         gas pressure: ', &
    ((wippflo_auxvars(idof,i)%pres(gid),i=1,n),idof=0,2)
  write(86,100) '   capillary pressure: ', &
    ((wippflo_auxvars(idof,i)%pres(cpid),i=1,n),idof=0,2)
  write(86,100) '    liquid saturation: ', &
    ((wippflo_auxvars(idof,i)%sat(lid),i=1,n),idof=0,2)
  write(86,100) '       gas saturation: ', &
    ((wippflo_auxvars(idof,i)%sat(gid),i=1,n),idof=0,2)
  write(86,100) 'liquid density [kmol]: ', &
    ((wippflo_auxvars(idof,i)%den(lid),i=1,n),idof=0,2)
  write(86,100) '  liquid density [kg]: ', &
    ((wippflo_auxvars(idof,i)%den_kg(lid),i=1,n),idof=0,2)
  write(86,100) '   gas density [kmol]: ', &
    ((wippflo_auxvars(idof,i)%den(gid),i=1,n),idof=0,2)
  write(86,100) '     gas density [kg]: ', &
    ((wippflo_auxvars(idof,i)%den_kg(gid),i=1,n),idof=0,2)
  write(86,*)
  write(86,100) '      liquid mobility: ', &
    ((wippflo_auxvars(idof,i)%mobility(lid),i=1,n),idof=0,2)
  write(86,100) '         gas mobility: ', &
    ((wippflo_auxvars(idof,i)%mobility(gid),i=1,n),idof=0,2)
  write(86,100) '   effective porosity: ', &
    ((wippflo_auxvars(idof,i)%effective_porosity,i=1,n),idof=0,2)
  
  close(86)

end subroutine WIPPFloOutputAuxVars2

! ************************************************************************** !

subroutine WIPPFloAuxVarSingleDestroy(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  implicit none

  type(wippflo_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call WIPPFloAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)  

end subroutine WIPPFloAuxVarSingleDestroy

! ************************************************************************** !

subroutine WIPPFloAuxVarArray1Destroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  implicit none

  type(wippflo_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call WIPPFloAuxVarStrip(auxvars(iaux))
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)  

end subroutine WIPPFloAuxVarArray1Destroy

! ************************************************************************** !

subroutine WIPPFloAuxVarArray2Destroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  implicit none

  type(wippflo_auxvar_type), pointer :: auxvars(:,:)
  
  PetscInt :: iaux, idof
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call WIPPFloAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)  

end subroutine WIPPFloAuxVarArray2Destroy

! ************************************************************************** !

subroutine WIPPFloAuxVarStrip(auxvar)
  ! 
  ! WIPPFloAuxVarDestroy: Deallocates a general auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  type(wippflo_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%pres)  
  call DeallocateArray(auxvar%sat)  
  call DeallocateArray(auxvar%den)  
  call DeallocateArray(auxvar%den_kg)  
  call DeallocateArray(auxvar%mobility)  
  if (associated(auxvar%d)) then
    deallocate(auxvar%d)
    nullify(auxvar%d)
  endif
  
end subroutine WIPPFloAuxVarStrip

! ************************************************************************** !

subroutine WIPPFloAuxDestroy(aux)
  ! 
  ! Deallocates a general auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  type(wippflo_type), pointer :: aux
  
  if (.not.associated(aux)) return
  
  call WIPPFloAuxVarDestroy(aux%auxvars)
  call WIPPFloAuxVarDestroy(aux%auxvars_bc)
  call WIPPFloAuxVarDestroy(aux%auxvars_ss)

  call DeallocateArray(aux%upwind_direction)
  call DeallocateArray(aux%upwind_direction_bc)
  call DeallocateArray(aux%inactive_rows_local)
  call DeallocateArray(aux%inactive_rows_local_ghosted)
  call DeallocateArray(aux%row_zeroing_array)

  if (associated(aux%wippflo_parameter)) then
  endif
  nullify(aux%wippflo_parameter)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine WIPPFloAuxDestroy

end module WIPP_Flow_Aux_module
