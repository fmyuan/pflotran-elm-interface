module PM_TOWG_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use PM_Base_Aux_module
  use AuxVars_TOWG_module

  implicit none

  private


  !global variable to TOWG
  PetscReal, public :: towg_window_epsilon = 1.d-4
  !PetscReal, public :: towg_fmw_comp(3) = & ! initialised after EOSread
  !           [UNINITIALIZED_DOUBLE,UNINITIALIZED_DOUBLE,UNINITIALIZED_DOUBLE]
  PetscReal, pointer, public :: towg_fmw_comp(:)
  PetscInt, pointer, public :: towg_dof_to_primary_variable(:,:)

  PetscInt, public :: towg_debug_cell_id = UNINITIALIZED_INTEGER
  PetscReal, parameter, public :: towg_pressure_scale = 1.d0
  
  PetscReal, public :: val_tl_omega = 0.0d0
  PetscReal, public :: fmis_sl      = 0.0d0
  PetscReal, public :: fmis_su      = 1.0d0
  PetscBool, public :: fmis_is_zero = .false.
  PetscBool, public :: fmis_is_unity= .false.

  PetscBool, public :: towg_isothermal = PETSC_FALSE
  PetscBool, public :: towg_no_gas = PETSC_FALSE
  PetscBool, public :: towg_no_oil = PETSC_FALSE
  PetscInt, public :: towg_miscibility_model = UNINITIALIZED_INTEGER

  ! list of TOWG paramters
  PetscInt, parameter, public :: TOWG_PREV_TS = 1
  PetscInt, parameter, public :: TOWG_PREV_IT = 2

  ! thermodynamic state of fluid ids - for BLACK OIL
  PetscInt, parameter, public :: TOWG_NULL_STATE = 0
  PetscInt, parameter, public :: TOWG_LIQ_OIL_STATE = 1
  PetscInt, parameter, public :: TOWG_LIQ_GAS_STATE = 2
  PetscInt, parameter, public :: TOWG_THREE_PHASE_STATE = 3
  PetscInt, parameter, public :: TOWG_ANY_STATE = 4

  ! Primary DOF indices -
  PetscInt, parameter, public :: TOWG_OIL_PRESSURE_DOF = 1
  PetscInt, parameter, public :: TOWG_OIL_SATURATION_DOF = 2
  PetscInt, parameter, public :: TOWG_GAS_SATURATION_2PH_DOF = 2
  PetscInt, parameter, public :: TOWG_GAS_SATURATION_3PH_DOF = 3
  PetscInt, parameter, public :: TOWG_BUBBLE_POINT_3PH_DOF   = 3 !Variable substitutes with gas saturation
  PetscInt, parameter, public :: TOWG_SOLV_SATURATION_DOF = 4
  PetscInt, parameter, public :: TOWG_3CMPS_ENERGY_DOF = 4
  PetscInt, parameter, public :: TOWG_SOLV_TL_ENERGY_DOF = 5
  !towg_energy_dof assigned TOWG_3CMPS_ENERGY_DOF or TOWG_SOLV_TL_ENERGY_DOF
  !in PMTOWGCreate
  PetscInt, public :: towg_energy_dof = UNINITIALIZED_INTEGER


  ! Equation indices -
  PetscInt, parameter, public :: TOWG_LIQ_EQ_IDX = 1
  PetscInt, parameter, public :: TOWG_OIL_EQ_IDX = 2
  PetscInt, parameter, public :: TOWG_GAS_EQ_IDX = 3
  PetscInt, parameter, public :: TOWG_SOLV_EQ_IDX = 4
  PetscInt, parameter, public :: TOWG_3CMPS_ENERGY_EQ_IDX = 4
  PetscInt, parameter, public :: TOWG_SOLV_TL_ENERGY_EQ_IDX = 5
  PetscInt, public :: towg_energy_eq_idx = UNINITIALIZED_INTEGER

  !Indices used to map aux_real for condition values
  PetscInt, parameter, public :: TOWG_OIL_PRESSURE_INDEX = 1
  PetscInt, parameter, public :: TOWG_OIL_SATURATION_INDEX = 2
  PetscInt, parameter, public :: TOWG_GAS_SATURATION_INDEX = 3
  PetscInt, parameter, public :: TOWG_BUBBLE_POINT_INDEX   = 3 !Variable substitutes with gas saturation
  PetscInt, parameter, public :: TOWG_SOLV_SATURATION_INDEX = 4
  PetscInt, parameter, public :: TOWG_X_GAS_IN_OIL_INDEX = 5
  PetscInt, parameter, public :: TOWG_TEMPERATURE_INDEX = 6
  PetscInt, parameter, public :: TOWG_LIQUID_FLUX_INDEX = 7
  PetscInt, parameter, public :: TOWG_OIL_FLUX_INDEX = 8
  PetscInt, parameter, public :: TOWG_GAS_FLUX_INDEX = 9
  PetscInt, parameter, public :: TOWG_SOLV_FLUX_INDEX = 10
  PetscInt, parameter, public :: TOWG_ENERGY_FLUX_INDEX = 11
  PetscInt, parameter, public :: TOWG_LIQ_CONDUCTANCE_INDEX = 12
  PetscInt, parameter, public :: TOWG_OIL_CONDUCTANCE_INDEX = 13
  PetscInt, parameter, public :: TOWG_GAS_CONDUCTANCE_INDEX = 14
  PetscInt, parameter, public :: TOWG_SOLVENT_CONDUCTANCE_INDEX = 15
  PetscInt, parameter, public :: TOWG_MAX_INDEX = 15

  !Indices used to map aux_int_var for condition values
  PetscInt, parameter, public :: TOWG_STATE_INDEX = 1

  !flags to identify type of auxvar update
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_BOUNDARY = 2

  !alternative density computation for tl4p, sometimes problematic?
  PetscBool, public :: TL4P_altDensity = PETSC_FALSE
  !switches for some minor changes
  PetscBool, public :: TL4P_slv_sat_truncate,TL4P_safemobs

  ! it might be required for thermal diffusion terms and tough conv criteria
  type, public :: towg_parameter_type
     !  PetscReal, pointer :: diffusion_coefficient(:) ! (iphase)
     !  PetscReal :: newton_inf_scaled_res_tol
     !  PetscBool :: check_post_converged
  end type towg_parameter_type

  !if required, could add other intermediate classes:
  type, public, extends(pm_base_aux_type) :: pm_towg_aux_type
    type(towg_parameter_type), pointer :: parameter
    class(auxvar_towg_type), pointer :: auxvars(:,:)
    class(auxvar_towg_type), pointer :: auxvars_bc(:)
    class(auxvar_towg_type), pointer :: auxvars_ss(:)
  contains
    !add bound-procedure
    procedure, public :: Init => InitTOWGAuxVars
    procedure, public :: FieldVolRefAve => TOWGAuxFieldVolRefAve
    procedure, public :: GetLocalSol    => TOWGGetLocalSol
    procedure, public :: IsSolventModel
    !procedure, public :: Perturb => PerturbTOilIms
  end type pm_towg_aux_type

  interface TOWGAuxVarStrip
    module procedure TOWGAuxVarArray1Strip
    module procedure TOWGAuxVarArray2Strip
  end interface TOWGAuxVarStrip


  !pointing to null() function
  procedure(TOWGAuxVarComputeDummy), pointer :: TOWGAuxVarCompute => null()
  procedure(TOWGAuxVarPerturbDummy), pointer :: TOWGAuxVarPerturb => null()

  abstract interface
    subroutine TOWGAuxVarComputeDummy(x,auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
      use Option_module
      use AuxVars_TOWG_module
      use Global_Aux_module
      use EOS_Water_module
      use EOS_Oil_module
      use Characteristic_Curves_module
      use Material_Aux_class
      implicit none
      type(option_type) :: option
      class(characteristic_curves_type) :: characteristic_curves
      PetscReal :: x(option%nflowdof)
      class(auxvar_towg_type) :: auxvar
      type(global_auxvar_type) :: global_auxvar
      class(material_auxvar_type) :: material_auxvar
      PetscInt :: natural_id
    end subroutine TOWGAuxVarComputeDummy

    subroutine TOWGAuxVarPerturbDummy(auxvar,global_auxvar, &
                                      material_auxvar, &
                                      characteristic_curves,natural_id, &
                                      option)
      use AuxVars_TOWG_module
      use Option_module
      use Characteristic_Curves_module
      use Global_Aux_module
      use Material_Aux_class
      implicit none
      type(option_type) :: option
      PetscInt :: natural_id
      type(auxvar_towg_type) :: auxvar(0:)
      type(global_auxvar_type) :: global_auxvar
      class(material_auxvar_type) :: material_auxvar
      class(characteristic_curves_type) :: characteristic_curves

    end subroutine

  end interface

  public :: TOWGAuxCreate, &
            TOWGAuxDestroy, &
            TOWGAuxVarStrip, &
            TOWGAuxVarCompute, &
            TOWGImsAuxVarComputeSetup, &
            TOWGTLAuxVarComputeSetup, &
            TOWGAuxVarPerturb, &
            TOWGBlackOilAuxVarComputeSetup, &
            TL4PAuxVarComputeSetup

contains

! ************************************************************************** !

function TOWGAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object for TOWG
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/05/16
  !

  use Option_module
  use EOS_Oil_module
  use EOS_Gas_module
  use EOS_Slv_module

  implicit none

  type(option_type) :: option

  class(pm_towg_aux_type), pointer :: TOWGAuxCreate

  class(pm_towg_aux_type), pointer :: aux

  allocate(towg_fmw_comp(option%nflowspec))

!--In TOWG case the gas FMW must be defined in the input deck

  if ( Uninitialized(EOSGasGetFMW()) ) then
    option%io_buffer = 'TOWG: gas FMW not initialised. ' // &
                       'Define its value in the the input deck' // &
                       ' or add EOS GAS card to default to FMWAIR'
    call printErrMsg(option)
  endif

  towg_fmw_comp(1) = FMWH2O
  towg_fmw_comp(2) = EOSOilGetFMW()
  towg_fmw_comp(3) = EOSGasGetFMW()

!--If solvent option, set up solvent molecular weight

  if ( towg_miscibility_model == TOWG_SOLVENT_TL ) then

    if ( Uninitialized(EOSSlvGetFMW()) ) then
      option%io_buffer = 'Solvent FMW not initialised. ' // &
                         'Define its value in the the input deck' // &
                         ' or add EOS SOLVENT card to default to FMWCO2'
      call printErrMsg(option)
    endif
    towg_fmw_comp(4) =EOSSlvGetFMW()
  end if

  allocate( towg_dof_to_primary_variable(1:option%nflowdof,1:3) )

  select case(towg_miscibility_model)
    case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF,TOWG_BLACK_OIL)
      towg_dof_to_primary_variable(1:option%nflowdof,1:3) = &
        reshape([TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, & ! LiQUID_STATE
                TOWG_BUBBLE_POINT_INDEX, TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_GAS_SATURATION_INDEX, &  ! GAS_STATE (not used yet,for Rv)
                TOWG_BUBBLE_POINT_INDEX, TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &  ! TWO_PHASE_STATE
                TOWG_GAS_SATURATION_INDEX,TOWG_TEMPERATURE_INDEX], &
                shape(towg_dof_to_primary_variable))
    case(TOWG_SOLVENT_TL)
      towg_dof_to_primary_variable(1:option%nflowdof,1:3) = &
        reshape([TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &
                TOWG_BUBBLE_POINT_INDEX, TOWG_SOLV_SATURATION_INDEX, &
                TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_GAS_SATURATION_INDEX, & ! GAS_STATE (not used yet,for Rv)
                TOWG_BUBBLE_POINT_INDEX, TOWG_SOLV_SATURATION_INDEX, &
                TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &
                TOWG_GAS_SATURATION_INDEX, TOWG_SOLV_SATURATION_INDEX, &
                TOWG_TEMPERATURE_INDEX], &
                shape(towg_dof_to_primary_variable))
  end select

  !allocate here to define this is a pm_toil_ims_aux_type
  allocate(aux)

  call PMBaseAuxInit(aux)

  !nullify here and not in the parent class, because auxvars are mode dependent
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)

  allocate(aux%parameter)

  !PO - to allocate when supporting diffusion the oil and gas phase
  !allocate(aux%parameter%diffusion_coefficient(option%nphase))
  !aux%parameter%diffusion_coefficient(option%liquid_phase) = &
  !                                                 UNINITIALIZED_DOUBLE
  !aux%parameter%diffusion_coefficient(option%oil_phase) = &
  !                                                 UNINITIALIZED_DOUBLE
  !aux%parameter%diffusion_coefficient(option%gas_phase) = 2.13d-5
  !aux%parameter%newton_inf_scaled_res_tol = 1.d-50
  !aux%parameter%check_post_converged = PETSC_FALSE


  TOWGAuxCreate => aux

end function TOWGAuxCreate

! ************************************************************************** !

subroutine InitTOWGAuxVars(this,grid,num_bc_connection, &
                              num_ss_connection,option)
  !
  ! Initialize pm_towg_auxvars
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/07/16
  !

  use Option_module
  use Grid_module

  implicit none

  class(pm_towg_aux_type) :: this
  PetscInt :: num_bc_connection
  PetscInt :: num_ss_connection
  type(grid_type) :: grid
  type(option_type) :: option

  PetscInt :: ghosted_id, iconn, local_id
  PetscInt :: idof

  if (option%flow%numerical_derivatives         .or. &
      option%flow%numerical_derivatives_compare .or. &
      option%flow%num_as_alyt_derivs)             then
    allocate(this%auxvars(0:option%nflowdof,grid%ngmax))
    do ghosted_id = 1, grid%ngmax
      do idof = 0, option%nflowdof
        call this%auxvars(idof,ghosted_id)%Init(option)
        if (towg_miscibility_model == TOWG_TODD_LONGSTAFF ) then
          call this%auxvars(idof,ghosted_id)%InitTL(option)
        end if
        if (    (towg_miscibility_model == TOWG_BLACK_OIL ) &
            .or.(towg_miscibility_model == TOWG_SOLVENT_TL) ) then
          call this%auxvars(idof,ghosted_id)%InitBO(option)
        end if
      end do
    end do
  else
    allocate(this%auxvars(0:0,grid%ngmax))
    do ghosted_id = 1, grid%ngmax
      !do idof = 0, option%nflowdof
        call this%auxvars(0,ghosted_id)%Init(option)
        if (towg_miscibility_model == TOWG_TODD_LONGSTAFF ) then
          call this%auxvars(0,ghosted_id)%InitTL(option)
        end if
        if (    (towg_miscibility_model == TOWG_BLACK_OIL ) &
            .or.(towg_miscibility_model == TOWG_SOLVENT_TL) ) then
          call this%auxvars(0,ghosted_id)%InitBO(option)
        end if
      !end do
    end do
  endif

  this%num_aux = grid%ngmax

  if (num_bc_connection > 0) then
    allocate(this%auxvars_bc(num_bc_connection))
    do iconn = 1, num_bc_connection
      call this%auxvars_bc(iconn)%Init(option)
      if (towg_miscibility_model == TOWG_TODD_LONGSTAFF ) then
        call this%auxvars_bc(iconn)%InitTL(option)
      end if
      if (    (towg_miscibility_model == TOWG_BLACK_OIL ) &
          .or.(towg_miscibility_model == TOWG_SOLVENT_TL) ) then
        call this%auxvars_bc(iconn)%InitBO(option)
      end if
    enddo
  endif
  this%num_aux_bc = num_bc_connection

  if (num_ss_connection > 0) then
    allocate(this%auxvars_ss(num_ss_connection))
    do iconn = 1, num_ss_connection
      call this%auxvars_ss(iconn)%Init(option)
      if (towg_miscibility_model == TOWG_TODD_LONGSTAFF ) then
        call this%auxvars_ss(iconn)%InitTL(option)
      end if
      if (    (towg_miscibility_model == TOWG_BLACK_OIL ) &
          .or.(towg_miscibility_model == TOWG_SOLVENT_TL) ) then
        call this%auxvars_ss(iconn)%InitBO(option)
      end if
    enddo
  endif
  this%num_aux_ss = num_ss_connection

  call PMBaseAuxSetup(this,grid,option)

end subroutine InitTOWGAuxVars

! ************************************************************************** !

subroutine TOWGImsAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell for TOWGIms
  !
  ! Author: Paolo Orsini
  ! Date: 12/02/16
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Oil_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar ! passing this for salt conc.
                                            ! not currenty used
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  PetscInt :: natural_id !only for debugging/print out - currently not used

  PetscInt :: wid, oid, gid
  PetscReal :: cell_pressure, wat_sat_pres
  PetscReal :: krw,dkrw_sato,dkrw_satg,dkrw_satw, visw
  PetscReal :: kro,dkro_sato,dkro_satg, viso
  PetscReal :: krg,dkrg_sato,dkrg_satg, visg
  PetscReal :: sat_liq_gas, sat_tot_liq
  PetscReal :: dummy, dummy2

  PetscErrorCode :: ierr

  ! from SubsurfaceSetFlowMode
  ! option%liquid_phase = 1           ! liquid_phase = water phase
  ! option%oil_phase = 2              ! oil_pressure
  ! option%gas_phase = 3              ! gas_pressure
  wid = option%liquid_phase
  oid = option%oil_phase
  gid = option%gas_phase

  auxvar%effective_porosity = 0.d0
  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%pc = 0.d0
  auxvar%den = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%mobility = 0.d0
  auxvar%H = 0.d0
  auxvar%U = 0.d0
  auxvar%temp = 0.0d0

  !passing auxvars given by the solution variables

  auxvar%pres(oid) = x(TOWG_OIL_PRESSURE_DOF)
  auxvar%sat(oid) = x(TOWG_OIL_SATURATION_DOF)
  auxvar%sat(gid) = x(TOWG_GAS_SATURATION_3PH_DOF)
  auxvar%temp = x(towg_energy_dof)

  auxvar%sat(wid) = 1.d0 - auxvar%sat(oid) - auxvar%sat(gid)

!  Get oil-water and oil-gas capilliary pressures
!  Note Pw=Po-Pcow, Pg=Po+Pcog

  call characteristic_curves%oil_wat_sat_func% &
            CapillaryPressure(auxvar%sat(wid), &
                              auxvar%pc(wid),dummy,option,auxvar%table_idx)
                              
  call characteristic_curves%oil_gas_sat_func% &
            CapillaryPressure(auxvar%sat(gid), &
                              auxvar%pc(oid),dummy,option,auxvar%table_idx)                            
                              
  auxvar%pres(wid) = auxvar%pres(oid) - auxvar%pc(wid)
  auxvar%pres(gid) = auxvar%pres(oid) + auxvar%pc(oid)

  cell_pressure = max(auxvar%pres(wid),auxvar%pres(oid),auxvar%pres(gid))

  ! calculate effective porosity as a function of pressure
  if (option%iflag /= TOWG_UPDATE_FOR_BOUNDARY) then
    auxvar%effective_porosity = material_auxvar%porosity_base

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                auxvar%effective_porosity,dummy)
    endif
    if (option%iflag /= TOWG_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = auxvar%effective_porosity
    endif
  endif

  ! UPDATE THERMODYNAMIC PROPERTIES FOR BOTH PHASES!!!

  ! Liquid phase thermodynamic properties
  ! using cell_pressure (which is the max press)? or %pres(wid)?
  call EOSWaterDensity(auxvar%temp,cell_pressure, &
                       auxvar%den_kg(wid),auxvar%den(wid),ierr)
  call EOSWaterEnthalpy(auxvar%temp,cell_pressure,auxvar%H(wid),ierr)
  auxvar%H(wid) = auxvar%H(wid) * 1.d-6 ! J/kmol -> MJ/kmol

  ! MJ/kmol comp                  ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
  auxvar%U(wid) = auxvar%H(wid) - (cell_pressure / auxvar%den(wid) * 1.d-6)

  ! ADD HERE BRINE dependency. Two options (see mphase)
  ! - salinity constant in space and time (passed in option%option%m_nacl)
  ! - salt can be trasnported by RT (sequential coupling) and passed
  !   and passed with global_auxvar%m_nacl
  !  ! Assign salinity
  !  m_na=option%m_nacl; m_cl=m_na; m_nacl=m_na
  !  if (option%ntrandof > 0) then
  !    m_na = global_auxvar%m_nacl(1)
  !    m_cl = global_auxvar%m_nacl(2)
  !    m_nacl = m_na
  !    if (m_cl > m_na) m_nacl = m_cl
  !  endif
  !
  !  ! calculate density for pure water
  !  call EOSWaterDensityEnthalpy(t,pw,dw_kg,dw_mol,hw,ierr)
  !  !..................
  !  xm_nacl = m_nacl*FMWNACL
  !  xm_nacl = xm_nacl/(1.D3 + xm_nacl)
  !  ! corrects water densit previously calculated as pure water
  !  call EOSWaterDensityNaCl(t,p,xm_nacl,dw_kg)
  !  ! water viscosity dependence on salt concetration, but no derivatives
  !  !  call EOSWaterViscosityNaCl(t,p,xm_nacl,visl)
  !  call EOSWaterViscosity(t,pw,sat_pressure,0.d0,visl,dvdt,dvdp,dvdps,ierr)

  call EOSOilDensityEnergy(auxvar%temp,auxvar%pres(oid),&
                           auxvar%den(oid),auxvar%H(oid), &
                           auxvar%U(oid),ierr)

  auxvar%den_kg(oid) = auxvar%den(oid) * EOSOilGetFMW()

  auxvar%H(oid) = auxvar%H(oid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(oid) = auxvar%U(oid) * 1.d-6 ! J/kmol -> MJ/kmol

  !compute gas properties (default is air - but methane can be set up)

  call EOSGasDensityEnergy(auxvar%temp,auxvar%pres(gid),auxvar%den(gid), &
                           auxvar%H(gid),auxvar%U(gid),ierr)

  auxvar%den_kg(gid) = auxvar%den(gid) * EOSGasGetFMW()
  auxvar%H(gid) = auxvar%H(gid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(gid) = auxvar%U(gid) * 1.d-6 ! J/kmol -> MJ/kmol

  ! compute water mobility (rel. perm / viscostiy)

  call characteristic_curves%wat_rel_perm_func_owg% &
                  RelativePermeability(auxvar%sat(wid),krw,dkrw_satw, &
                                            option,auxvar%table_idx)

  call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,ierr)

  ! use cell_pressure; cell_pressure - psat calculated internally
  call EOSWaterViscosity(auxvar%temp,cell_pressure,wat_sat_pres,visw,ierr)

  auxvar%mobility(wid) = krw/visw

  ! compute oil mobility (rel. perm / viscosity)
  !PO: For testing use krow - to be extended to kro (Eclipse model)
  call characteristic_curves%ow_rel_perm_func_owg% &
                       RelativePermeability(auxvar%sat(oid),kro,dkro_sato, &
                                                     option,auxvar%table_idx)

  call EOSOilViscosity(auxvar%temp,auxvar%pres(oid), &
                       auxvar%den(oid), viso, ierr)

  auxvar%mobility(oid) = kro/viso

  !compute gas mobility (rel. perm / viscosity)
  call characteristic_curves%gas_rel_perm_func_owg% &
                       RelativePermeability(auxvar%sat(gid),krg,dummy,&
                                                     option,auxvar%table_idx)

  !currently only a viscosity model for air or constant value

  call EOSGasViscosity(auxvar%temp,auxvar%pres(gid), &
                       auxvar%pres(gid),auxvar%den(gid),visg,ierr)

  auxvar%mobility(gid) = krg/visg


end subroutine TOWGImsAuxVarCompute

subroutine getBlackOilComposition(bubble_point,temperature,table_idxs,xo,xg,&
                                  dxo_dpb,dxg_dpb,dxo_dt,dxg_dt)

!------------------------------------------------------------------------------
! Find oil phase composition (xo,xg) as a function of bubble point
! Used in TOWG_BLACK_OIL and TOWG_SOLVENT_TL
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Dec 2017
!------------------------------------------------------------------------------

  use EOS_Oil_module
  use EOS_Gas_module

  PetscReal, intent(in ) :: bubble_point
  PetscReal, intent(in ) :: temperature
  PetscInt, pointer, intent(inout) :: table_idxs(:)
  PetscReal, intent(out) :: xo
  PetscReal, intent(out) :: xg
  PetscReal              :: rs_molar
  PetscInt               :: ierr

  PetscReal, optional, intent(out) :: dxo_dpb,dxg_dpb,dxo_dt,dxg_dt
  PetscReal              :: drs_molar_dt,drs_molar_dpb
  PetscBool              :: getDerivs

  getDerivs = PETSC_FALSE
  if (present(dxo_dpb).AND.present(dxg_dpb).AND.present(dxo_dt).AND.present(dxg_dt))then
    getDerivs = PETSC_TRUE
  endif

!--Get molar Rsm=(moles gas)/(moles oil)---------------------------------------

  if (getDerivs) then
    call EOSOilRS(temperature,bubble_point,rs_molar,drs_molar_dt,drs_molar_dpb,ierr,table_idxs)
  else
    call EOSOilRS(temperature,bubble_point,rs_molar,ierr,table_idxs)
  endif

!--Get oil and gas mole fractions----------------------------------------------
!
! xo=(moles oil)/(moles oil+moles gas)=1/(1+(moles gas)/(moles oil))=1/(1+Rsm)
! xg=1-xo=1-1/(1+Rsm)=(1+Rsm-1)/(1+Rsm)=Rsm/(1+Rsm)
!
!------------------------------------------------------------------------------

  xo=1.0d0   /(1.0d0+rs_molar)
  xg=rs_molar/(1.0d0+rs_molar)

  if (getDerivs) then
    dxo_dpb = -drs_molar_dpb/(1.0d0+rs_molar)/(1.0d0+rs_molar)
    dxg_dpb = -dxo_dpb
    ! (see above: xg = 1 - xo)
    dxo_dt = -drs_molar_dt/(1.0d0+rs_molar)/(1.0d0+rs_molar)
    dxg_dt = -dxo_dt
    if (dabs(dxo_dpb) < 1.d-10) then
      print *, "dxo dpb is ", dxo_dpb
    endif
  endif


end subroutine getBlackOilComposition

! ************************************************************************** !

subroutine TOWGBlackOilAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, &
                                     characteristic_curves,natural_id,option)
!------------------------------------------------------------------------------
! Auxillary variable calculation for black oil
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Jul 2017
!------------------------------------------------------------------------------


  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Oil_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use Derivatives_utilities_module

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  PetscInt :: natural_id !only for debugging/print out - currently not used

  PetscInt :: wid, oid, gid, istate
  PetscReal :: cell_pressure, wat_sat_pres
  PetscReal :: krw, visw
  PetscReal :: kro, viso
  PetscReal :: krg, visg
  PetscReal :: dummy,po,pb
  PetscReal :: dkrw_sato,dkrw_satg,dkrw_satw
  PetscReal :: dkro_sato,dkro_satg
  PetscReal :: dkrg_sato,dkrg_satg
  PetscErrorCode :: ierr
  PetscBool :: isSat
  PetscReal :: oneminuseps
  PetscReal :: deno,cr,cvisc,ho,uo,crusp,cvusp

  PetscReal, parameter :: eps_oil   = 1.0d-6
  PetscReal, parameter :: epss=1.0d-4
  PetscReal, parameter :: epsp=1.0d3
  PetscBool :: getDerivs
  
  PetscReal :: dpc_o_dsg, dpc_w_dsw
  PetscReal :: dcr_dt,dcr_dpb
  PetscReal :: dcrusp_dpo,dcrusp_dpb,dcrusp_dt
  PetscReal :: cor,one_p_crusp,dcor_dpo,dcor_dpb,dcor_dt
  PetscReal :: dps_dt
  PetscReal :: dvw_dt, dvw_dp
  PetscReal, dimension(1:option%nflowdof) :: dmobw,dmobo,dmobg
  PetscReal :: dviso_dpb,dcvisc_dT,dcvisc_dPb
  PetscReal :: dcvusp_dpo,dcvusp_dpb,dcvusp_dt
  PetscReal :: one_p_cvusp
  PetscReal :: dvo_dp,dvo_dpb,dvo_dt
  PetscReal :: dvg_dt,dvg_dpcomp,dvg_dpgas,d_deno_dpb,d_xg_dpb,d_xo_dpb
  PetscReal :: worker
  PetscReal, dimension(1:option%nflowdof) :: D_worker
  PetscReal, dimension(1:option%nflowdof) :: D_cell_pres
  PetscReal, dimension(1:option%nflowdof) :: D_visc,D_kr
  !PetscReal :: d_cellpres_dso,d_cellpres_dsg

  PetscInt :: dof_op,dof_osat,dof_gsat,dof_temp

  PetscReal :: stpt, drpt, dff, ndrv

  dof_op = TOWG_OIL_PRESSURE_DOF
  dof_osat = TOWG_OIL_SATURATION_DOF
  dof_gsat = TOWG_GAS_SATURATION_3PH_DOF
  dof_temp = towg_energy_dof



  if (towg_analytical_derivatives) then
    if (.NOT. auxvar%has_derivs) then
      ! how did this happen?
      option%io_buffer = 'towg bo auxvars: towg_analytical_derivatives is true, &
                          but auxvar%has_derivs is false, should both be true. &
                          How did this happen?'
      call printErrMsg(option)
    endif

    auxvar%D_pres = 0.d0
    auxvar%D_sat = 0.d0
    auxvar%D_pc = 0.d0
    auxvar%D_den = 0.d0
    auxvar%D_den_kg = 0.d0
    auxvar%D_mobility = 0.d0
    auxvar%D_por = 0.d0

    auxvar%D_H = 0.d0
    auxvar%D_U = 0.d0

    ! also zero out all black oil derivatives
    auxvar%bo%D_xo = 0.d0
    auxvar%bo%D_xg = 0.d0

    getDerivs = PETSC_TRUE
  else 
    getDerivs = PETSC_FALSE
  endif

!==============================================================================
!  Initialise
!==============================================================================

  oneminuseps=1.0d0-epss

! 'Liquid' is really the aqueous, water phase

  wid = option%liquid_phase
  oid = option%oil_phase
  gid = option%gas_phase

  auxvar%effective_porosity = 0.d0
  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%pc = 0.d0
  auxvar%den = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%mobility = 0.d0
  auxvar%H = 0.d0
  auxvar%U = 0.d0

!==============================================================================
!  Check state to see if saturated
!==============================================================================

  istate=global_auxvar%istate

  if( istate==TOWG_THREE_PHASE_STATE ) then
    isSat=PETSC_TRUE
  else
    isSat=PETSC_FALSE
  endif

!==============================================================================
! Getting auxvars as given by the solution variables,
! allowing for saturated/undersaturated switch
!==============================================================================

  auxvar%pres(oid) = x(TOWG_OIL_PRESSURE_DOF  )
  auxvar%sat (oid) = x(TOWG_OIL_SATURATION_DOF)
  if( isSat ) then
    auxvar%sat(gid)        = x(TOWG_GAS_SATURATION_3PH_DOF)
    auxvar%bo%bubble_point = auxvar%pres(oid)
  else
    auxvar%bo%bubble_point = x(TOWG_GAS_SATURATION_3PH_DOF)
    auxvar%sat(gid)        = 0.0
  endif
  auxvar%temp = x(towg_energy_dof)

!==============================================================================
! Check if this state still valid and flip if not (but not on diff call)
!==============================================================================

  if (option%iflag /= TOWG_UPDATE_FOR_DERIVATIVE) then
    if( isSat ) then
      if( (auxvar%sat(gid)<0.0d0) .and. (auxvar%sat(oid)>eps_oil) ) then
! Gas saturation has gone negative and significant oil in cell
        global_auxvar%istate          =TOWG_LIQ_OIL_STATE
        auxvar%sat(gid)               =0.0d0
        auxvar%bo%bubble_point        =auxvar%pres(oid)-epsp
        x(TOWG_BUBBLE_POINT_3PH_DOF)  =auxvar%pres(oid)-epsp
      endif
    else
      if(      (auxvar%bo%bubble_point > auxvar%pres(oid)) &
          .or. (auxvar%sat(oid)        < eps_oil         ) ) then
! Bubble point has exceeded oil pressure or no significant oil in cell
        global_auxvar%istate          =TOWG_THREE_PHASE_STATE
        auxvar%bo%bubble_point        =auxvar%pres(oid)
        auxvar%sat(gid)               =epss
        x(TOWG_GAS_SATURATION_3PH_DOF)=epss
! Make sure the extra gas does not push the water saturation negative
        if( auxvar%sat(oid) > oneminuseps ) then
          auxvar%sat(oid)             =oneminuseps
          x(TOWG_OIL_SATURATION_DOF)  =oneminuseps
        endif
      endif
    endif
  endif

  ! get state again in case it swapped. It's possible for 
  ! state to change and isSat still be same otherwise.
  ! Alternatively could flip isSat inside if stats above.
  istate=global_auxvar%istate

  if( istate==TOWG_THREE_PHASE_STATE ) then
    isSat=PETSC_TRUE
  else
    isSat=PETSC_FALSE
  endif

!==============================================================================
! Set up the final saturation value (water)
!==============================================================================

  auxvar%sat(wid) = 1.d0 - auxvar%sat(oid) - auxvar%sat(gid)



!==============================================================================
! Extract solution into local scalars
!==============================================================================

  po=auxvar%pres(oid)
  pb=auxvar%bo%bubble_point


!==============================================================================

  ! trivial saturation derivatives: 
  if (getDerivs) then
    auxvar%D_sat(oid,dof_osat) = 1.d0 ! diff oil sat by oil sat
    auxvar%D_sat(wid,dof_osat) = -1.d0 ! diff liquid sat by oil sat
    if (isSat) then
      auxvar%D_sat(wid,dof_gsat) = -1.d0 ! diff liquid sat by gas sat
      auxvar%D_sat(gid,dof_gsat) = 1.d0 ! diff gas sat by gas sat
    endif
  endif
  ! trivial pressure derivatives: 
  if (getDerivs) then
    auxvar%D_pres(oid,dof_op) = 1.d0 ! diff oil pres by oil pres
    ! these will be overwritten with cap pres derivs below
    auxvar%D_pres(wid,dof_op) = 1.d0 ! diff liquid pres by oil pres
    auxvar%D_pres(gid,dof_op) = 1.d0 ! diff gas pres by oil pres
  endif

!==============================================================================
! Look up the Rs value
!==============================================================================

  if (getDerivs) then
    call getBlackOilComposition(auxvar%bo%bubble_point,auxvar%temp, &
                                auxvar%table_idx,auxvar%bo%xo,auxvar%bo%xg, &
                                d_xo_dpb,d_xg_dpb,&
                                auxvar%bo%D_xo(dof_temp),auxvar%bo%D_xg(dof_temp))
    if (isSat) then
      ! bubble point is cell pressure
      auxvar%bo%D_xo(dof_op) = d_xo_dpb
      auxvar%bo%D_xg(dof_op) = d_xg_dpb
    else
      ! bubble point is a solution variable
      auxvar%bo%D_xo(dof_gsat) = d_xo_dpb
      auxvar%bo%D_xg(dof_gsat) = d_xg_dpb
    endif
  else
    call getBlackOilComposition(auxvar%bo%bubble_point,auxvar%temp, &
                                auxvar%table_idx,auxvar%bo%xo,auxvar%bo%xg )
  endif

  if (auxvar%bo%xo < 0.d0) then
    print *, "xo negative ", auxvar%bo%xo, " pb is ", auxvar%bo%bubble_point
    option%io_buffer = 'xo has gone negative; xo and bubble point are'
    call printMsg(option)
    write(option%io_buffer,*) auxvar%bo%xo
    call printMsg(option)
    write(option%io_buffer,*) auxvar%bo%bubble_point
    call printMsg(option)
  endif

!==============================================================================
!  Get the capillary pressures
!==============================================================================


! phase notation for the cap pressure can be confusing. There are (ndof -1)
! entries in the array as opposed to ndof. 
! pc(wid) = pc_{wo} : between oil and water
! pc(oid) = pc_{og} : between oil and gas

! and of course analagously for D_pc.
! Writing into or reading from D_pc(gid) is a mistake.
! Note this implies a hardcoding of the values of w/o/gid, since we must
! always have gid > wid,oid.

  call characteristic_curves%oil_wat_sat_func% &
            CapillaryPressure(auxvar%sat(wid), &
                              auxvar%pc(wid),dpc_w_dsw,option,auxvar%table_idx)

  auxvar%pc(oid) = 0.0d0
                
  call characteristic_curves%oil_gas_sat_func% &
            CapillaryPressure(auxvar%sat(gid), &
                              auxvar%pc(oid),dpc_o_dsg,option,auxvar%table_idx)

  if (getDerivs) then
    ! deriv of pc between oil and water, w.r.t. oil sat:
    auxvar%D_pc(wid,dof_osat) = -dpc_w_dsw
    ! deriv of pc between oil and water, w.r.t. gas sat:
    if (isSat) then
      auxvar%D_pc(wid,dof_gsat) = -dpc_w_dsw
    endif
    
    if (isSat) then
      ! deriv of pc between gas and oil w.r.t. oil sat:
      auxvar%D_pc(oid,dof_gsat) = dpc_o_dsg
      ! deriv of pc between gas and oil w.r.t. gas sat:
      ! 0 b/c g sat is only arg of the cp routine
      !auxvar%D_pc(oid,dof_osat) = -dpc_o_dsg
      ! gas saturation isn't a solution variable is not sat state
    endif

    ! pressure saved below
  endif


!==============================================================================
!  Get the phase pressures
!==============================================================================

  auxvar%pres(wid) = auxvar%pres(oid) - auxvar%pc(wid)
  auxvar%pres(gid) = auxvar%pres(oid) + auxvar%pc(oid)

  if (getDerivs) then
    auxvar%D_pres(wid,dof_op) = 1.d0
    ! pc(wid) has derivatives w.r.t. oil and gas sat
    auxvar%D_pres(wid,dof_osat) =  -auxvar%D_pc(wid,dof_osat)
    if (isSat) then
      auxvar%D_pres(wid,dof_gsat) = -auxvar%D_pc(wid,dof_gsat)
    endif
    auxvar%D_pres(gid,dof_op) = 1.d0
    ! pc(gid) has derivatives w.r.t. oil and gas sat
    auxvar%D_pres(gid,dof_osat) =  auxvar%D_pc(oid,dof_osat)
    if (isSat) then
      auxvar%D_pres(gid,dof_gsat) =  auxvar%D_pc(oid,dof_gsat)
    endif

  endif

  cell_pressure = max(auxvar%pres(wid),auxvar%pres(oid),auxvar%pres(gid))

  ! For analytical derivatives:
  ! In the case of nonzero cap pressure we may have dependencies on saturations,
  ! through the cap pressures, this is dealt with next:

  if (getDerivs) then
    D_cell_pres = 0.d0
    D_cell_pres(dof_op) = 1.d0
    if (cell_pressure == auxvar%pres(wid)) then
      D_cell_pres(:) = auxvar%D_pres(wid,:)
    elseif (cell_pressure == auxvar%pres(gid)) then
      D_cell_pres(:) = auxvar%D_pres(gid,:)
    !else ! cell pressure = oil pressure and we need to do nothing here
    endif
  endif


!==============================================================================
!  Get rock properties
!==============================================================================

!--Calculate effective porosity as a function of pressure----------------------

  if (option%iflag /= TOWG_UPDATE_FOR_BOUNDARY) then
    auxvar%effective_porosity = material_auxvar%porosity_base

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                auxvar%effective_porosity,dummy)
      if (getDerivs) then
        ! it's a cell pressure derivative:
        auxvar%D_por(dof_op) = dummy
        auxvar%D_por(dof_osat) = dummy*D_cell_pres(dof_osat)
        auxvar%D_por(dof_gsat) = dummy*D_cell_pres(dof_gsat)
        ! could also just be:
        ! D_por = dummy*D_cell_pres?

      endif
    endif
    if (option%iflag /= TOWG_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = auxvar%effective_porosity
    endif
  endif

!==============================================================================
!  Get fluid properties
!==============================================================================

!------------------------------------------------------------------------------
! Water phase thermodynamic properties
!------------------------------------------------------------------------------

  if (getDerivs) then
    call EOSWaterDensity(auxvar%temp,cell_pressure, &
                      auxvar%den_kg(wid),auxvar%den(wid), &
                      auxvar%D_den(wid,dof_op), &
                      auxvar%D_den(wid,dof_temp), ierr)

        ! pressure deriv is a a cell pressure derivative:
        auxvar%D_den(wid,dof_osat) = D_cell_pres(dof_osat)*auxvar%D_den(wid,dof_op)
        auxvar%D_den(wid,dof_gsat) = D_cell_pres(dof_gsat)*auxvar%D_den(wid,dof_op)
        ! (implicitly: D_cell_pres(dof_oid) = 1 and a multplication by 1 was not
        !  done here)

    call EOSWaterEnthalpy(auxvar%temp, &
                       cell_pressure, &
                       auxvar%H(wid), &
                       auxvar%D_H(wid,dof_op), &
                       auxvar%D_H(wid,dof_temp), &
                       ierr)

    ! pressure deriv is a a cell pressure derivative:
    auxvar%D_H(wid,dof_osat) = D_cell_pres(dof_osat)*auxvar%D_H(wid,dof_op)
    auxvar%D_H(wid,dof_gsat) = D_cell_pres(dof_gsat)*auxvar%D_H(wid,dof_op)
    ! (implicitly: D_cell_pres(dof_oid) = 1 and a multplication by 1 was not
    !  done here)
    
    auxvar%D_H(wid,:) = auxvar%D_H(wid,:) * 1.d-6 ! J/kmol -> MJ/kmol

    ! derivatives corresponding to computation of U(wid) below:
    auxvar%D_U(wid,:) = auxvar%D_H(wid,:)                                          &
                      - 1.d-6                                                      &
                      * DivRule(cell_pressure,D_cell_pres,                         &
                                auxvar%den(wid),auxvar%D_den(wid,:),option%nflowdof )
  else
    call EOSWaterDensity(auxvar%temp,cell_pressure, &
                         auxvar%den_kg(wid),auxvar%den(wid),ierr)
    call EOSWaterEnthalpy(auxvar%temp,cell_pressure,auxvar%H(wid),ierr)
  endif
  auxvar%H(wid) = auxvar%H(wid) * 1.d-6 ! J/kmol -> MJ/kmol

  ! MJ/kmol comp                  ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
  auxvar%U(wid) = auxvar%H(wid) - (cell_pressure / auxvar%den(wid) * 1.d-6)

  if (getDerivs) then
    auxvar%D_den_kg(oid,:) = auxvar%D_den(oid,:) * EOSOilGetFMW()
  endif



!------------------------------------------------------------------------------
! Gas phase thermodynamic properties. Can be ideal gas default or PVTG table
!------------------------------------------------------------------------------

  if (getDerivs) then
    !         note here derivative w.r.t. gas pressure saved as derivatives w.r.t. oil presures \/                                            
    call EOSGasDensityEnergy(auxvar%temp,auxvar%pres(gid),auxvar%den(gid),auxvar%D_den(gid,dof_temp),auxvar%D_den(gid,dof_op), &
                                     auxvar%H(gid),auxvar%D_H(gid,dof_temp),auxvar%D_H(gid,dof_op),auxvar%U(gid),&
                                     auxvar%D_U(gid,dof_temp),&
                                     auxvar%D_U(gid,dof_op),ierr,auxvar%table_idx)

    
    ! pressure derivatives:
    auxvar%D_den(gid,dof_osat) = auxvar%D_den(gid,dof_op)*D_cell_pres(dof_osat)
    if (isSat) then
      auxvar%D_den(gid,dof_gsat) = auxvar%D_den(gid,dof_op)*D_cell_pres(dof_gsat)
     endif
    auxvar%D_H(gid,dof_osat) = auxvar%D_H(gid,dof_op)*D_cell_pres(dof_osat)
    if (isSat) then
      auxvar%D_H(gid,dof_gsat) = auxvar%D_H(gid,dof_op)*D_cell_pres(dof_gsat)
     endif
    auxvar%D_U(gid,dof_osat) = auxvar%D_U(gid,dof_op)*D_cell_pres(dof_osat)
    if (isSat) then
      auxvar%D_U(gid,dof_gsat) = auxvar%D_U(gid,dof_op)*D_cell_pres(dof_gsat)
    endif

    ! scaling:
    auxvar%D_H(gid,:) = auxvar%D_H(gid,:) * 1.d-6 ! J/kmol -> MJ/kmol
    auxvar%D_U(gid,:) = auxvar%D_U(gid,:) * 1.d-6 ! J/kmol -> MJ/kmol

  else
    call EOSGasDensityEnergy(auxvar%temp,auxvar%pres(gid),auxvar%den(gid), &
                             auxvar%H(gid),auxvar%U(gid),ierr,auxvar%table_idx)
  endif

  auxvar%den_kg(gid) = auxvar%den(gid) * EOSGasGetFMW()

  if (getDerivs) then
    auxvar%D_den_kg(gid,:) = auxvar%D_den(gid,:) * EOSGasGetFMW()
  endif

  auxvar%H(gid) = auxvar%H(gid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(gid) = auxvar%U(gid) * 1.d-6 ! J/kmol -> MJ/kmol

!------------------------------------------------------------------------------
! Oil phase thermodynamic properties from PVCO
! Note that auxvar%table_idx for table lookup is passed after ierr
! because it is an optional argument (needed only if table lookup)
!------------------------------------------------------------------------------

  if (getDerivs) then
    call EOSOilDensityEnergy(auxvar%temp,po,deno, &
                             auxvar%D_den(oid,dof_temp),auxvar%D_den(oid,dof_op), &
                             ho,auxvar%D_H(oid,dof_temp),auxvar%D_H(oid,dof_op), &
                             uo,auxvar%D_U(oid,dof_temp),auxvar%D_U(oid,dof_op), &
                             ierr) 

    ! Density and compressibility look-up at bubble point
    call EOSOilDensity(auxvar%temp,pb,deno,auxvar%D_den(oid,dof_temp),d_deno_dpb,ierr,auxvar%table_idx)
    call EOSOilCompressibility(auxvar%temp,pb,cr,dcr_dt,dcr_dpb,ierr,auxvar%table_idx)

    if (isSat) then
      ! pb is pressure:
      auxvar%D_den(oid,dof_op) = d_deno_dpb
    else
      ! pb is a solution variable:
      auxvar%D_den(oid,dof_gsat) = d_deno_dpb
    endif

  else
    call EOSOilDensityEnergy  (auxvar%temp,po,deno,ho,uo,ierr,auxvar%table_idx)
  ! Density and compressibility look-up at bubble point
    call EOSOilDensity        (auxvar%temp,pb,deno      ,ierr,auxvar%table_idx)
    call EOSOilCompressibility(auxvar%temp,pb,cr        ,ierr,auxvar%table_idx)
  endif

!  Correct for undersaturation: correction not yet available for energy
! --------- Correct for undersaturation ---------------------------------------
  crusp=cr*(po-pb)

  if (getDerivs) then
    ! differentiate the correction for undersaturation (see below)
    ! NOTE: might have to adjust this when corrections are introduced for energy
    if (isSat) then
      ! leave density derivs w.r.t. pres, temp, and bubble point as set above since the correction
      ! is zero.
    else
      dcrusp_dpo =  cr !  +  dcr_dpo*(po-pb) but we know dcr_dpo is zero
      dcrusp_dpb = dcr_dpb*(po-pb) - cr
      dcrusp_dt =  dcr_dt*(po-pb)

      cor = 1 + crusp + 0.5*crusp*crusp
      ! dcor/dx = dcrusp/dx + crusp*dcrusp/dx  = (1 + crusp)*dcrusp/dx
      ! so
      one_p_crusp = 1.d0 + crusp
      dcor_dpo = one_p_crusp*dcrusp_dpo
      dcor_dpb = one_p_crusp*dcrusp_dpb
      dcor_dt = one_p_crusp*dcrusp_dt

      auxvar%D_den(oid,dof_op) = deno*dcor_dpo ! +  ddeno_dpo*cor but we know ddeno_dpo is zero
      auxvar%D_den(oid,dof_gsat) = auxvar%D_den(oid,dof_gsat)*cor + deno*dcor_dpb
      auxvar%D_den(oid,dof_temp) = auxvar%D_den(oid,dof_temp)*cor + deno*dcor_dt
    endif

  endif

  !crusp=cr*(po-pb) ! moved above
  deno=deno*(1.0+crusp*(1.0+0.5*crusp))
  auxvar%den(oid)=deno
  auxvar%H  (oid)=ho
  auxvar%U  (oid)=uo
! --------- /Correct for undersaturation --------------------------------------


!------------------------------------------------------------------------------
! Correct oil phase molar density and enthalpy for oil composition
! EOSOilDensity returns oil moles/volume in oil phase,
! but really have (1+Rsmolar) times as many moles when dissolved gas included
! 1+Rsmolar=1+xg/xo=(xo+xg)/xo=1/xo.
! Note xo=1/(1+Rsmolar) so cannot be zero for finite looked up Rsmolar
!------------------------------------------------------------------------------

!-----------Correct oil molar density------------------------------------------
  if (getDerivs) then
    ! assume we've already gotten D_den and D_xo
    ! correct (i.e., D_xo(dof_gsat) is correct regardless of state; it's
    ! nonzero if that index is for pb, zero if otherwise, sim. for 
    ! D_xo(dof_op);
    ! then jsut do this:
    auxvar%D_den(oid,:) = DivRule(auxvar%den(oid),auxvar%D_den(oid,:), &
                                  auxvar%bo%xo,auxvar%bo%D_xo,option%nflowdof      )
  endif

  auxvar%den(oid)=auxvar%den(oid)/auxvar%bo%xo
!----------/Correct oil molar density-----------------------------------------



! Get oil mass density as (mixture oil molar density).(mixture oil molecular weight)

!----------Oil mass density---------------------------------------------------
  if (getDerivs) then
    ! assume D_xo and D_xg have been previously set up
    ! correctly independently of state (so derivatives w.r.t. temp, and w.r.t.
    ! pres OR pb are in arrays correctly), then just differentiate 
    ! mechanically:

    worker = auxvar%bo%xo*EOSOilGetFMW()     &
           + auxvar%bo%xg*EOSGasGetFMW()
    D_worker = auxvar%bo%D_xo*EOSOilGetFMW() &
             + auxvar%bo%D_xg*EOSGasGetFMW()
    auxvar%D_den_kg(oid,:) = ProdRule(auxvar%den(oid),auxvar%D_den(oid,:), &
                                      worker,D_worker,option%nflowdof)
  endif

  auxvar%den_kg(oid) = auxvar%den(oid) * ( auxvar%bo%xo*EOSOilGetFMW() &
                                          +auxvar%bo%xg*EOSGasGetFMW() )
!----------/Oil mass density---------------------------------------------------


!----------H and U scaling-----------------------------------------------------
  auxvar%H(oid) = auxvar%H(oid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(oid) = auxvar%U(oid) * 1.d-6 ! J/kmol -> MJ/kmol

  if (getDerivs) then
    auxvar%D_H(oid,:) = auxvar%D_H(oid,:) * 1.d-6 ! J/kmol -> MJ/kmol
    auxvar%D_U(oid,:) = auxvar%D_U(oid,:) * 1.d-6 ! J/kmol -> MJ/kmol
  endif
!----------/H and U scaling----------------------------------------------------

!------------------------------------------------------------------------------
! Get oil enthalpy/oil mole in oil phase.
! Calculation was for pure oil, but really have total molar enthalpy of:
! oil hydrocarbon enthalpy/oil hydrocarbon mole= xo*(oil enthalpy/oil mole)
!                                               +xg*(gas enthalpy/gas mole)
!------------------------------------------------------------------------------

!----------Oil enthalpy/oil mole in oil phase----------------------------------
  if (getDerivs) then

    auxvar%D_H(oid,:) = prodrule(auxvar%bo%xo,auxvar%bo%D_xo,           &
                                 auxvar%H(oid),auxvar%D_H(oid,:),option%nflowdof)  &
                      + ProdRule(auxvar%bo%xg,auxvar%bo%D_xg,           &
                                 auxvar%H(gid),auxvar%D_H(gid,:),option%nflowdof) 

    auxvar%D_U(oid,:) = prodrule(auxvar%bo%xo,auxvar%bo%d_xo,           &
                                 auxvar%U(oid),auxvar%D_U(oid,:),option%nflowdof)  &
                      + ProdRule(auxvar%bo%xg,auxvar%bo%d_xg,           &
                                 auxvar%U(gid),auxvar%D_U(gid,:),option%nflowdof) 
  endif

  auxvar%H(oid) = auxvar%bo%xo*auxvar%H(oid)+auxvar%bo%xg*auxvar%H(gid)
  auxvar%U(oid) = auxvar%bo%xo*auxvar%U(oid)+auxvar%bo%xg*auxvar%U(gid)
!----------/Oil enthalpy/oil mole in oil phase---------------------------------

!===============================================================================
! Fluid mobility calculation
!===============================================================================

!-------------------------------------------------------------------------------
!  Water mobility (rel. perm / viscosity)
!-------------------------------------------------------------------------------

  call characteristic_curves%wat_rel_perm_func_owg% &
                    RelativePermeability(auxvar%sat(wid),krw,dkrw_satw,option, &
                    auxvar%table_idx)

  if (getDerivs) then
    call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,dps_dt,ierr)

    call EOSWaterViscosity(auxvar%temp, cell_pressure, &
                           wat_sat_pres, dps_dt, visw, &
                           dvw_dt,  dvw_dp, ierr)

        ! pressure deriv (dvw_dp) is a a cell pressure derivative:
        D_visc = 0.d0
        D_visc(dof_op) = dvw_dp

        D_visc(dof_osat)= D_cell_pres(dof_osat)*dvw_dp
        D_visc(dof_gsat)= D_cell_pres(dof_gsat)*dvw_dp
        D_visc(dof_temp)= dvw_dt
  else
    call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,ierr)

  ! use cell_pressure; cell_pressure - psat calculated internally
    call EOSWaterViscosity(auxvar%temp,cell_pressure,wat_sat_pres,visw,ierr)
  endif

  auxvar%mobility(wid) = krw/visw


  if (getDerivs) then

    ! since here we have to account that visc may have derivatives in all variables,
    ! might as well apply divrule to krw/visw instead of call the specific mobility
    ! derivatives routine
    D_kr = 0.d0
    ! char curve returns deriv w.r.t  water sat so account for that
    D_kr(dof_osat) = -dkrw_satw
    if (isSat) then 
      D_kr(dof_gsat) = -dkrw_satw 
    endif
    auxvar%D_mobility(wid, :) =  DivRule(krw,D_kr,                  &
                                         visw,D_visc,option%nflowdof )
  endif

!-------------------------------------------------------------------------------
! Oil mobility (rel. perm / viscosity)
!-------------------------------------------------------------------------------

  call characteristic_curves%oil_rel_perm_func_owg% &
                  RelativePermeability(auxvar%sat(oid),auxvar%sat(gid),kro, &
                                 dkro_sato,dkro_satg,option,auxvar%table_idx)

!--If PVCO defined in EOS OIL, the viscosities are extracted via table lookup--

  if (getDerivs) then
    call EOSOilViscosity(auxvar%temp,pb,auxvar%den(oid),viso,dvo_dt,dviso_dpb,ierr,auxvar%table_idx)
    call EOSOilViscosibility(auxvar%temp,pb,cvisc,dcvisc_dT,dcvisc_dPb,ierr, &
                             auxvar%table_idx)

  else
  ! Viscosity and viscosibility look-up at bubble point
    call EOSOilViscosity    (auxvar%temp,pb, &
                             auxvar%den(oid), viso, ierr,auxvar%table_idx)
    call EOSOilViscosibility(auxvar%temp,pb,cvisc,ierr,auxvar%table_idx)
  endif

!----------Correct oil viscosity-----------------------------------------------
  cvusp=cvisc*(po-pb)

  if (getDerivs) then
    ! get derivatives of corrected viso
    if (isSat) then
      ! 1) cvusp = 0; cor = 1 constants
      ! 2) pb is really oil pressure so pb derivs equal to  po derivs

      dvo_dp = dviso_dpb
    else ! unsat state
      ! 1) cvusp, cor now vary
      ! 2) pb derivs diferent from po derivs
      ! 3) have to apply corrector visc = cor*visc, so apply corresponding correctors
      !    to derivs too

      dcvusp_dpo = cvisc ! cvisc is independent of po here
      dcvusp_dpb = dcvisc_dpb*(po-pb) - cvisc
      dcvusp_dt = dcvisc_dt*(po-pb) 

      cor = (1.0+cvusp*(1.0+0.5*cvusp))
      ! dcor/dx = dcvusp/dx + cvusp*dcvusp so:
      one_p_cvusp = 1.d0 + cvusp
      dcor_dpo = dcvusp_dpo*one_p_cvusp 
      dcor_dpb = dcvusp_dpb*one_p_cvusp 
      dcor_dt = dcvusp_dt*one_p_cvusp 

      dvo_dp = viso*dcor_dpo ! viso is independent of po here
      dvo_dpb = viso*dcor_dpb + dviso_dpb*cor
      dvo_dt = viso*dcor_dt + dvo_dt*cor

    endif
  endif

  !cvusp=cvisc*(po-pb) !!! moved above
  viso=viso*(1.0+cvusp*(1.0+0.5*cvusp))
!----------/Correct oil viscosity-----------------------------------------------

  auxvar%mobility(oid) = kro/viso


  if (getDerivs) then
    ! mobility derivatives:
    D_kr = 0.d0
    D_kr(dof_osat) = dkro_sato
    if (isSat) then 
      D_kr(dof_gsat) = dkro_satg
    endif
    D_visc = 0.d0
    D_visc(dof_op) = dvo_dp
    if (.NOT. isSat) then 
      D_visc(dof_gsat) = dvo_dpb
    endif
    D_visc(dof_temp) = dvo_dt

    auxvar%D_mobility(oid, :) =  DivRule(kro,D_kr,                  &
                                         viso,D_visc,option%nflowdof )

  endif

!-------------------------------------------------------------------------------
!  Gas mobility (rel. perm / viscosity)
!-------------------------------------------------------------------------------

  call characteristic_curves%gas_rel_perm_func_owg% &
                   RelativePermeability(auxvar%sat(gid),krg, &
                                dkrg_satg,option,auxvar%table_idx)

  if (getDerivs) then
    ! seems to want derivatves of comp pressure w.r.t. pressure and temp, we use 
    ! comp pres = pres so set these to zero
    dummy = 0.d0

    call  EOSGasViscosity(auxvar%temp,auxvar%pres(gid),auxvar%pres(gid),&
                          auxvar%den(gid), &
                          auxvar%D_den(gid,dof_temp), auxvar%D_den(gid,dof_temp), auxvar%D_den(gid,dof_op), &
                          dummy,dummy, &                      ! dPcomp_dT, dPcomp_dPgas
                          visg, dvg_dt, dvg_dpcomp, dvg_dpgas, ierr, &
                          auxvar%table_idx)


        ! the pres deriv is a cell pres deriv so:
        D_visc = 0.d0
        D_visc(dof_op) = dvg_dpgas

        D_visc(dof_osat)= D_cell_pres(dof_osat)*dvg_dpgas
        D_visc(dof_gsat)= D_cell_pres(dof_gsat)*dvg_dpgas

        D_visc(dof_temp)= dvg_dt

  else
  !--If PVDG defined in EOS OIL, the viscosities are extracted via table lookup--
    call EOSGasViscosity(auxvar%temp,auxvar%pres(gid), &
                         auxvar%pres(gid),auxvar%den(gid),visg,ierr,&
                         auxvar%table_idx)
  endif

  auxvar%mobility(gid) = krg/visg
  if (getDerivs) then

    ! see comment for water mobility
    D_kr = 0.d0
    if (isSat) then 
      D_kr(dof_gsat) = dkrg_satg
    endif
    auxvar%D_mobility(gid, :) =  DivRule(krg,D_kr,                  &
                                         visg,D_visc,option%nflowdof )
  endif



end subroutine TOWGBlackOilAuxVarCompute

! ************************************************************************** !

! this is no longer used:
subroutine MobilityDerivs_TOWG_BO(dmob,kr,visc,dkr_dso,dkr_dsg,dvisc_dpo,dvisc_dpb,dvisc_dt,isSat,ndof)
  implicit none
  PetscInt :: ndof
  PetscReal, dimension(1:ndof) :: dmob
  PetscReal :: kr,visc 
  PetscReal :: dkr_dso, dkr_dsg
  PetscReal :: dvisc_dpo, dvisc_dpb, dvisc_dt
  PetscBool :: isSat


  PetscReal :: one_over_visc_sq

  dmob = 0.d0

  one_over_visc_sq = 1.d0/visc/visc

  ! oil pressure derivative
  dmob(TOWG_OIL_PRESSURE_DOF) = -1.d0*kr*dvisc_dpo*one_over_visc_sq

  ! oil saturation derivative 
  dmob(TOWG_OIL_SATURATION_DOF) = dkr_dso/visc

  ! temperature derivative
  dmob(towg_energy_dof) = -1.d0*kr*dvisc_dt*one_over_visc_sq

  ! gas saturation OR bubble point derivative
  if (isSat) then
    ! gas saturation
    dmob(TOWG_GAS_SATURATION_3PH_DOF) = dkr_dsg/visc
  else
    ! bubble point
    dmob(TOWG_GAS_SATURATION_3PH_DOF) = -1.d0*kr*dvisc_dpb*one_over_visc_sq
  endif


end subroutine MobilityDerivs_TOWG_BO


! ************************************************************************** !

subroutine TL4PAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, &
                             characteristic_curves,natural_id,option)
!------------------------------------------------------------------------------
! Auxillary variable calculation for 4-phase Todd-Longstaff
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------
! Modified: Daniel Stone
! Date  : Feb 2019
! Reason : Adding analytical derivatives code
!------------------------------------------------------------------------------

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Oil_module
  use EOS_Gas_module
  use EOS_Slv_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use Derivatives_utilities_module

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  PetscInt :: natural_id !only for debugging/print out - currently not used

  PetscInt :: wid, oid, gid, sid, istate
  PetscReal :: cell_pressure, wat_sat_pres
  PetscReal :: so,sg,sw,ss,sv,sh,t
  PetscReal,parameter::szero=0.0
  PetscReal :: uoil,ugas,uslv,uvap
  PetscReal :: krh,dkrh_sath,dkrh_satz
  PetscReal :: kroi,krgi,krsi
  PetscReal :: krom,krgm,krsm
  PetscReal :: krw ,visw,viswtl
  PetscReal :: kro ,viso,visotl
  PetscReal :: krg ,visg,visgtl
  PetscReal :: krs ,viss,visstl
  PetscReal :: dummy,po,pb
  PetscReal :: dkrw_sato,dkrw_satg,dkrw_satw
  PetscReal :: dkrg_sato,dkrg_satg
  PetscReal :: dkrv_sato,dkrv_satv,krvi
  PetscErrorCode :: ierr
  PetscBool isSat
  PetscReal :: oneminuseps
  PetscReal :: dummy2,cr,cvisc,ho,uo,crusp,cvusp
  PetscReal :: test
  PetscReal :: swcr,sgcr,sowcr,sogcr,swco,fm,fi,rat
  PetscReal :: deno,deng,dens,denotl,dengtl,denstl,krow,krog,sc,den,soa,sva,sumhydsat
  PetscReal :: swa
  PetscReal :: socrs,svcrs,krvm

  PetscInt :: dof_op,dof_osat,dof_gsat,dof_temp,dof_ssat,ndof,cploc

  PetscReal, parameter :: eps_oil=1.0d-5
  PetscReal, parameter :: eps_gas=1.0d-8
  PetscReal, parameter :: epss   =1.0d-4
  PetscReal, parameter :: epsp   =1.0d3
  PetscReal, parameter :: epseps =1.0d-10
  
!----------------- intermediate derivatives and related variables: ------------------------------------------------
  PetscBool :: getDerivs,mobSanityCheck                                                  ! utility flags
  PetscReal :: d_xo_dpb,d_xg_dpb                                                         ! xo and xg wrt bubble point
  PetscReal :: d_pcw_d_swa,d_pco_d_sva                                                   ! cap pres derivs, returned args
  PetscReal :: dx_dcell_pres,dx_dpres,dx_dt                                              ! used as return args for several routines
  PetscReal :: dh_dt,dh_dpres,du_dt,du_dpres                                             ! used as return args for several routines
  PetscReal :: ddeno_dpb,dcr_dt,dcr_dpb                                                  ! intermediates for oil den lookup and correction
  PetscReal :: dcrusp_dpo,dcrusp_dpb,dcrusp_dt,cor,one_p_crusp                           ! as above
  PetscReal :: dcor_dpo,dcor_dpb,dcor_dt,worker                                          ! as above
  PetscReal :: dps_dt                                                                    ! return arg from routine, water sat pres deriv
  PetscReal :: dkro_sato,dkrog_uoil                                                      ! kro and krog deriv routine return args
  PetscReal :: num                                                                       ! intermediate used for sevral krxx derivs
  PetscReal :: dvo_dt,dvo_dpb,dcvisc_dt,dcvisc_dpb                                       ! intermediates for oil visc lookup and correction
  PetscReal :: one_p_cvusp,dvo_dp,dviso_dpb                                              ! as above
  PetscReal :: dcvusp_dpo,dcvusp_dpb,dcvusp_dt                                           ! as above
  PetscReal :: denos,dengs,denogs,denos_pre,denog                                        ! for debugging, variables internal to tl visc/den
  ! arrays of derivatives, D_xx(i) is derivative of variable xx w.r.t. to variable index i
  ! (e.g. D_fm(dof_ssat) = deriv of fm w.r.t. solvent sat):
  PetscReal,dimension(1:option%nflowdof) :: D_fm,D_cell_pres                             ! fm and cell pressure 
  PetscReal,dimension(1:option%nflowdof) :: D_worker                                     ! utility used in a few places
  PetscReal,dimension(1:option%nflowdof) :: D_visc,D_kr                                  ! intermeds, for water mob derivs
  PetscReal,dimension(1:option%nflowdof) :: D_uoil,D_uvap                                ! derivs of uoil and uvap
  PetscReal,dimension(1:option%nflowdof) :: D_sv,D_sw,D_sh                               ! composite saturation derivatives
  PetscReal,dimension(1:option%nflowdof) :: D_num,D_den                                  ! intermediates used in several places
  PetscReal,dimension(1:option%nflowdof) :: D_krom,D_krvm,D_krgm,D_krsm                  ! various kr intermediates
  PetscReal,dimension(1:option%nflowdof) :: D_kroi,D_krh,D_krog,D_krow                   ! as above
  PetscReal,dimension(1:option%nflowdof) :: D_krgi,D_krsi                                ! as above
  PetscReal,dimension(1:option%nflowdof) :: D_viso,D_visg,D_viss                         ! viscosity derivs
  PetscReal,dimension(1:option%nflowdof) :: D_kro,D_krg,D_krs                            ! tl relperm derivs
  PetscReal,dimension(1:option%nflowdof) :: D_so,D_sg,D_ss                               ! for input to routine
  PetscReal,dimension(1:option%nflowdof) :: D_deno,D_deng,D_dens                         ! as above
  PetscReal,dimension(1:option%nflowdof) :: D_visotl,D_visgtl,D_visstl                   ! tl visc derivs
  PetscReal,dimension(1:option%nflowdof) :: D_denotl,D_dengtl,D_denstl                   ! tl den derivs
  PetscReal,dimension(1:option%nflowdof) :: D_denos,D_dengs,D_denogs,D_denos_pre,D_denog ! for debugging

  PetscInt :: idex,jdex
!-----------------/intermediate derivatives and related variables: -------------------------------------------------

  ! used for indexing into derivative arrays:
  dof_op = TOWG_OIL_PRESSURE_DOF
  dof_osat = TOWG_OIL_SATURATION_DOF
  dof_gsat = TOWG_GAS_SATURATION_3PH_DOF
  dof_temp = towg_energy_dof
  dof_ssat = TOWG_SOLV_SATURATION_DOF
  ndof = option%nflowdof
  mobSanityCheck = TL4P_safemobs

if (towg_analytical_derivatives) then
  if (.NOT. auxvar%has_derivs) then
    ! how did this happen?
    option%io_buffer = 'towg tl4p auxvars: towg_analytical_derivatives is true, &
                        but auxvar%has_derivs is false, should both be true. &
                        How did this happen?'
    call printErrMsg(option)
  endif

  auxvar%D_pres = 0.d0
  auxvar%D_sat = 0.d0
  auxvar%D_pc = 0.d0
  auxvar%D_den = 0.d0
  auxvar%D_den_kg = 0.d0
  auxvar%D_mobility = 0.d0
  auxvar%D_por = 0.d0

  auxvar%D_H = 0.d0
  auxvar%D_U = 0.d0

  ! also zero out all black oil derivatives
  auxvar%bo%D_xo = 0.d0
  auxvar%bo%D_xg = 0.d0

  getDerivs = PETSC_TRUE
else
  getDerivs = PETSC_FALSE
endif

!==============================================================================
!  Initialise
!==============================================================================

  krh      =0.0
  dkrh_sath=0.0
  dkrh_satz=0.0

  kroi=0.0
  krgi=0.0
  krsi=0.0

  krom=0.0
  krgm=0.0
  krsm=0.0

  krw =0.0
  kro =0.0
  krg =0.0
  krs =0.0

  krvi=0.0

  swcr =0.0
  sgcr =0.0
  sowcr=0.0
  sogcr=0.0
  swco =0.0

  fm =0.0
  fi =1.0

  oneminuseps=1.0d0-epss

! 'Liquid' is really the aqueous, water phase

  wid = option%liquid_phase
  oid = option%oil_phase
  gid = option%gas_phase
  sid = option%solvent_phase

  auxvar%effective_porosity = 0.d0
  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%pc = 0.d0
  auxvar%den = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%mobility = 0.d0
  auxvar%H = 0.d0
  auxvar%U = 0.d0

!==============================================================================
!  Check state to see if saturated
!==============================================================================

  istate=global_auxvar%istate

  if( istate==TOWG_THREE_PHASE_STATE ) then
    isSat=PETSC_TRUE
  else
    isSat=PETSC_FALSE
  endif

!==============================================================================
! Getting auxvars as given by the solution variables,
! allowing for saturated/undersaturated switch
!==============================================================================

  auxvar%pres(oid) = x(TOWG_OIL_PRESSURE_DOF  )
  auxvar%sat (oid) = x(TOWG_OIL_SATURATION_DOF)
  if( isSat ) then
    auxvar%sat(gid)        = x(TOWG_GAS_SATURATION_3PH_DOF)
    auxvar%bo%bubble_point = auxvar%pres(oid)
  else
    auxvar%bo%bubble_point = x(TOWG_BUBBLE_POINT_3PH_DOF)
    auxvar%sat(gid)        = 0.0
  endif
  auxvar%sat(sid)        = x(TOWG_SOLV_SATURATION_DOF)
  auxvar%temp = x(towg_energy_dof)

!==============================================================================
! Check if this state still valid and flip if not (but not on diff call)
!==============================================================================

  if (option%iflag /= TOWG_UPDATE_FOR_DERIVATIVE) then
    if( isSat ) then
      if( (auxvar%sat(gid)<(-epseps)) .and. (auxvar%sat(oid)>eps_oil) ) then
        ! warning: be aware that it can indeed happen that negative gas saturation will be accepted
        ! if oil sat is <= eps_oil.
! Gas saturation has gone negative and significant oil in cell
        global_auxvar%istate          =TOWG_LIQ_OIL_STATE
        auxvar%sat(gid)               =0.0d0
        auxvar%bo%bubble_point        =auxvar%pres(oid)-epsp
        x(TOWG_BUBBLE_POINT_3PH_DOF)  =auxvar%pres(oid)-epsp
        ! make sure this bool is still correct:
        isSat = PETSC_FALSE 
      endif
    else
      if(      (auxvar%bo%bubble_point > auxvar%pres(oid)) &
          .or. (auxvar%sat(oid)        < eps_oil         ) ) then
! Bubble point has exceeded oil pressure or no significant oil in cell
        global_auxvar%istate          =TOWG_THREE_PHASE_STATE
        auxvar%bo%bubble_point        =auxvar%pres(oid)
        if(auxvar%bo%bubble_point > auxvar%pres(oid)) then
          auxvar%sat(gid)=epss      ! Bubble point crossed
        else
          auxvar%sat(gid)=eps_gas   ! Run out of oil
        endif
        x(TOWG_GAS_SATURATION_3PH_DOF)=auxvar%sat(gid)
! Make sure the extra gas does not push the water saturation negative
        sumhydsat=auxvar%sat(oid)+auxvar%sat(sid)+auxvar%sat(gid)
        ! make sure this bool is still correct:
        isSat = PETSC_TRUE
        if( sumhydsat>1.0  ) then
          rat=1.0/sumhydsat
          auxvar%sat(oid)=auxvar%sat(oid)*rat
          auxvar%sat(gid)=auxvar%sat(gid)*rat
          auxvar%sat(sid)=auxvar%sat(sid)*rat
          x(TOWG_OIL_SATURATION_DOF    )=auxvar%sat(oid)
          x(TOWG_GAS_SATURATION_3PH_DOF)=auxvar%sat(gid)
          x(TOWG_SOLV_SATURATION_DOF   )=auxvar%sat(sid)
        endif
      endif
    endif
  endif

!==============================================================================
! Set up the final saturation value (water)
!==============================================================================

  auxvar%sat(wid) = 1.d0 - auxvar%sat(oid) - auxvar%sat(gid) - auxvar%sat(sid)

  ! trivial saturation derivatives:
  if (getDerivs) then
    auxvar%D_sat(oid,dof_osat) =  1.d0 ! diff oil sat by oil sat
    auxvar%D_sat(sid,dof_ssat) =  1.d0 ! diff solvent sat by solvent sat

    auxvar%D_sat(wid,dof_osat) = -1.d0 ! diff liquid sat by oil sat
    auxvar%D_sat(wid,dof_ssat) = -1.d0 ! diff liquid sat by solvent sat

    if (isSat) then
      auxvar%D_sat(gid,dof_gsat) = 1.d0 ! diff gas sat by gas sat
      auxvar%D_sat(wid,dof_gsat) = -1.d0 ! diff liquid sat by gas sat
    endif
  endif

!==============================================================================
! Extract solution into local scalars
!==============================================================================

  po=auxvar%pres(oid)
  pb=auxvar%bo%bubble_point

  so=auxvar%sat(oid)
  sg=auxvar%sat(gid)
  sw=auxvar%sat(wid)
  ss=auxvar%sat(sid)

  sv=sg+ss
  sh=so+sv

  t =auxvar%temp

!==============================================================================
! Look up the Rs value
!==============================================================================

  if (getDerivs) then
    call getBlackOilComposition(pb,t, &
                                auxvar%table_idx,auxvar%bo%xo,auxvar%bo%xg, &
                                d_xo_dpb,d_xg_dpb,&
                                auxvar%bo%D_xo(dof_temp),auxvar%bo%D_xg(dof_temp))
    if (isSat) then
      ! bubble point is cell pressure
      auxvar%bo%D_xo(dof_op) = d_xo_dpb
      auxvar%bo%D_xg(dof_op) = d_xg_dpb
    else
      ! bubble point is a solution variable
      auxvar%bo%D_xo(dof_gsat) = d_xo_dpb
      auxvar%bo%D_xg(dof_gsat) = d_xg_dpb
    endif
  else
    call getBlackOilComposition(pb,t, &
                                auxvar%table_idx,auxvar%bo%xo,auxvar%bo%xg )
  endif

!==============================================================================
!  Get the miscible-immiscible mixing fractions (functions of fs=Ss/(Sg+Ss))
!==============================================================================
  if (getDerivs) D_fm = 0.d0

  call TL4PMiscibilityFraction(sg,ss,fm,D_fm(dof_gsat),D_fm(dof_ssat))
  fi=1.0-fm

  if (getDerivs) then
    if (.NOT. isSat) D_fm(dof_gsat)= 0.d0
    ! store variables if debugging mode wants to:
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%fm = fm
      auxvar%tlT%D_fm = D_fm
    endif
  endif

!==============================================================================
!  Get the capillary pressures using three phase mode as functions of so,sv,sw
!==============================================================================

!--Prepare argument for lookups (these routines can fail if so+sv>1.0)

  soa=so
  sva=sv
!-- comment protection since this is now implemented within Pcow and Pcog
!  if( sh>=1.0d0 ) then
!    soa=(1.0d0-eps_oil)*so/sh
!    sva=(1.0d0-eps_oil)*sv/sh
!  endif
  swa = 1 - soa - sva
!--Pcow------------------------------------------------------------------------

  call characteristic_curves%oil_wat_sat_func% &
            CapillaryPressure(swa,auxvar%pc(wid),d_pcw_d_swa,option,auxvar%table_idx)

!--Pcog------------------------------------------------------------------------

  call characteristic_curves%oil_gas_sat_func% &
            CapillaryPressure(sva,auxvar%pc(oid),d_pco_d_sva,option,auxvar%table_idx)

  if (getDerivs) then
    ! recall these composite saturations
    ! swa = 1 - soa - sva
    ! sva = sv = sg + ss
    ! soa = so
    ! therefore
    ! swa = 1 - so - sg - ss
    ! derivs all -1
    !
    ! sva = sv = ss + sg
    ! d sva / dso = 0
    ! d sva / dsg = 1
    ! d sva / dss = 1

    ! deriv of pc between oil and water, w.r.t. oil sat:
    auxvar%D_pc(wid,dof_osat) =  -d_pcw_d_swa
    ! deriv of pc between oil and water, w.r.t. slv sat:
    auxvar%D_pc(wid,dof_ssat) =  -d_pcw_d_swa
    ! deriv of pc between oil and water, w.r.t. gas sat:
    if (isSat) then
      auxvar%D_pc(wid,dof_gsat) = -d_pcw_d_swa
    endif

    ! handle cap pres between vapour and oil:
    ! deriv of pc between vapour and oil w.r.t. slv sat:
    auxvar%D_pc(oid,dof_ssat) = d_pco_d_sva

    if (isSat) then
      ! deriv of pc between vapour and oil w.r.t. gas sat:
      auxvar%D_pc(oid,dof_gsat) = d_pco_d_sva
    endif

    !  Scale Pcog for degree of miscibility (derivs)
    ! note fi = 1 - fm so D_fi = -D_fm:
    auxvar%D_pc(oid,:) = ProdRule(fi,-D_fm,auxvar%pc(oid),auxvar%D_pc(oid,:),ndof)
  endif

!  Scale Pcog for degree of miscibility
  auxvar%pc(oid)=fi*auxvar%pc(oid)

!==============================================================================
!  Get the phase pressures
!==============================================================================

  auxvar%pres(wid) = po-auxvar%pc(wid)
  auxvar%pres(gid) = po+auxvar%pc(oid)

  auxvar%pres(sid) = auxvar%pres(gid)

  if (getDerivs) then
    ! first trivial presure derivatives:
    ! oil by oil:
    auxvar%D_pres(oid,dof_op) = 1.d0

    ! water pressure: 
    auxvar%D_pres(wid,:) = -auxvar%D_pc(wid,:) + auxvar%D_pres(oid,:)

    ! gas pressure:
    auxvar%D_pres(gid,:) = auxvar%D_pc(oid,:) + auxvar%D_pres(oid,:)

    ! solvent pressure: 
    auxvar%D_pres(sid,:) = auxvar%D_pres(gid,:)
  endif

  cell_pressure = max(auxvar%pres(wid),auxvar%pres(oid),auxvar%pres(gid))

  if (getDerivs) then
    ! find dex corresponding to cell pres:
    cploc = wid
    if (auxvar%pres(oid) > auxvar%pres(cploc)) then
      cploc = oid
    endif
    if (auxvar%pres(gid) > auxvar%pres(cploc)) then
      cploc = gid
    endif
    D_cell_pres = auxvar%D_pres(cploc,:)
  endif

  ! store variables if debugging mode wants to:
  if (auxvar%has_TL_test_object) then
    auxvar%tlT%cellpres= cell_pressure
    auxvar%tlT%D_cellpres= D_cell_pres
  endif

!==============================================================================
!  Get rock properties
!==============================================================================

!--Calculate effective porosity as a function of pressure----------------------

  if (option%iflag /= TOWG_UPDATE_FOR_BOUNDARY) then
    auxvar%effective_porosity = material_auxvar%porosity_base

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                auxvar%effective_porosity,dx_dcell_pres)
      if (getDerivs) then
        auxvar%D_por = dx_dcell_pres * D_cell_pres
      endif
    endif
    if (option%iflag /= TOWG_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = auxvar%effective_porosity
    endif
  endif

!==============================================================================
!  Get fluid properties
!==============================================================================

!------------------------------------------------------------------------------
! Water phase thermodynamic properties
!------------------------------------------------------------------------------

  if (getDerivs) then
    call EOSWaterDensity(t,cell_pressure, &
                         auxvar%den_kg(wid),auxvar%den(wid), &
                         dx_dcell_pres, &
                         dx_dt,ierr)

    auxvar%D_den(wid,:) =  dx_dcell_pres * D_cell_pres
    auxvar%D_den(wid,dof_temp) = auxvar%D_den(wid,dof_temp) + dx_dt

    ! kg density should be stored too for completeness:
    auxvar%D_den_kg(wid,:) = auxvar%D_den(wid,:) * FMWH2O

    call EOSWaterEnthalpy(auxvar%temp, &
                          cell_pressure, &
                          auxvar%H(wid), &
                          dx_dcell_pres, &
                          dx_dt, &
                          ierr)


    auxvar%D_H(wid,:) =  dx_dcell_pres * D_cell_pres
    auxvar%D_H(wid,dof_temp) = auxvar%D_H(wid,dof_temp) + dx_dt

    ! scaling:
    auxvar%D_H(wid,:) = auxvar%D_H(wid,:) * 1.d-6 ! J/kmol -> MJ/kmol

    ! derivatives corresponding to computation of U(wid) below:
    auxvar%D_U(wid,:) = auxvar%D_H(wid,:)                                          &
                      - 1.d-6                                                      &
                      * DivRule(cell_pressure,D_cell_pres,                         &
                                auxvar%den(wid),auxvar%D_den(wid,:),option%nflowdof )
  else
    call EOSWaterDensity(auxvar%temp,cell_pressure, &
                         auxvar%den_kg(wid),auxvar%den(wid),ierr)
    call EOSWaterEnthalpy(auxvar%temp,cell_pressure,auxvar%H(wid),ierr)
  endif

  auxvar%H(wid) = auxvar%H(wid) * 1.d-6 ! J/kmol -> MJ/kmol
  ! MJ/kmol comp                  ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
  auxvar%U(wid) = auxvar%H(wid) - (cell_pressure / auxvar%den(wid) * 1.d-6)



!------------------------------------------------------------------------------
! /end of Water phase thermodynamic properties
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Gas phase thermodynamic properties. Can be ideal gas default or PVTG table
!------------------------------------------------------------------------------

  if (getDerivs) then
    call EOSGasDensityEnergy(t,auxvar%pres(gid),auxvar%den(gid),dx_dt,dx_dpres, &
                                     auxvar%H(gid),dh_dt,dh_dpres,auxvar%U(gid),&
                                     du_dt,&
                                     du_dpres,ierr,auxvar%table_idx)

    ! pressure derivatives - den, h and u are functions of gas pressure and temp here
    auxvar%D_den(gid,:) = auxvar%D_pres(gid,:) * dx_dpres
    auxvar%D_den(gid,dof_temp) = auxvar%D_den(gid,dof_temp) + dx_dt

    auxvar%D_u(gid,:) = auxvar%D_pres(gid,:) * du_dpres
    auxvar%D_u(gid,dof_temp) = auxvar%D_u(gid,dof_temp) + du_dt

    auxvar%D_h(gid,:) = auxvar%D_pres(gid,:) * dh_dpres
    auxvar%D_h(gid,dof_temp) = auxvar%D_h(gid,dof_temp) + dh_dt


    auxvar%D_den_kg(gid,:) = auxvar%D_den(gid,:) * EOSGasGetFMW()
    ! scaling:
    auxvar%D_H(gid,:) = auxvar%D_H(gid,:) * 1.d-6 ! J/kmol -> MJ/kmol
    auxvar%D_U(gid,:) = auxvar%D_U(gid,:) * 1.d-6 ! J/kmol -> MJ/kmol

  else
    call EOSGasDensityEnergy(t,auxvar%pres(gid),auxvar%den(gid), &
                             auxvar%H(gid),auxvar%U(gid),ierr,auxvar%table_idx)
  endif

  auxvar%den_kg(gid) = auxvar%den(gid) * EOSGasGetFMW()

  auxvar%H(gid) = auxvar%H(gid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(gid) = auxvar%U(gid) * 1.d-6 ! J/kmol -> MJ/kmol


!------------------------------------------------------------------------------
! /end of Gas phase thermodynamic properties. 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Solvent phase thermodynamic properties. Can be ideal gas default or PVTS table
!------------------------------------------------------------------------------

  if (getDerivs) then

    call EOSSlvDensityEnergy(t,auxvar%pres(sid), &
                             auxvar%den(sid),dx_dt,dx_dpres,auxvar%H(sid),dh_dt,dh_dpres&
                             ,auxvar%U(sid),du_dt,du_dpres,ierr,&
                             auxvar%table_idx)

    ! pressure derivatives - den, h and u are functions of gas pressure and temp here
    auxvar%D_den(sid,:) = auxvar%D_pres(sid,:) * dx_dpres
    auxvar%D_den(sid,dof_temp) = auxvar%D_den(sid,dof_temp) + dx_dt

    auxvar%D_u(sid,:) = auxvar%D_pres(sid,:) * du_dpres
    auxvar%D_u(sid,dof_temp) = auxvar%D_u(sid,dof_temp) + du_dt

    auxvar%D_h(sid,:) = auxvar%D_pres(sid,:) * dh_dpres
    auxvar%D_h(sid,dof_temp) = auxvar%D_h(sid,dof_temp) + dh_dt

    auxvar%D_den_kg(sid,:) = auxvar%D_den(sid,:) * EOSSlvGetFMW()
    ! scaling:
    auxvar%D_H(sid,:) = auxvar%D_H(sid,:) * 1.d-6 ! J/kmol -> MJ/kmol
    auxvar%D_U(sid,:) = auxvar%D_U(sid,:) * 1.d-6 ! J/kmol -> MJ/kmol
  else
    call EOSSlvDensityEnergy(t,auxvar%pres(sid), &
                             auxvar%den(sid),auxvar%H(sid),auxvar%U(sid),ierr,&
                             auxvar%table_idx)

  endif

  auxvar%den_kg(sid) = auxvar%den(sid) * EOSSlvGetFMW()

  auxvar%H(sid) = auxvar%H(sid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(sid) = auxvar%U(sid) * 1.d-6 ! J/kmol -> MJ/kmol


!------------------------------------------------------------------------------
! /end of Solvent phase thermodynamic properties.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Oil phase thermodynamic properties from PVCO
! Note that auxvar%table_idx for table lookup is passed after ierr
! because it is an optional argument (needed only if table lookup)
!------------------------------------------------------------------------------

  if (getDerivs) then
   call EOSOilDensityEnergy(auxvar%temp,po,deno, &
                            auxvar%D_den(oid,dof_temp),auxvar%D_den(oid,dof_op), &
                            ho,auxvar%D_H(oid,dof_temp),auxvar%D_H(oid,dof_op), &
                            uo,auxvar%D_U(oid,dof_temp),auxvar%D_U(oid,dof_op), &
                            ierr,auxvar%table_idx)

   ! Density and compressibility look-up at bubble point
   call EOSOilDensity(auxvar%temp,pb,deno,auxvar%D_den(oid,dof_temp),ddeno_dpb,ierr,auxvar%table_idx)
   call EOSOilCompressibility(auxvar%temp,pb,cr,dcr_dt,dcr_dpb,ierr,auxvar%table_idx)

   if (isSat) then
     ! pb is pressure:
     auxvar%D_den(oid,dof_op) = ddeno_dpb
   else
     ! pb is a solution variable:
     auxvar%D_den(oid,dof_gsat) = ddeno_dpb
   endif

 else
   call EOSOilDensityEnergy  (auxvar%temp,po,deno,ho,uo,ierr,auxvar%table_idx)
   ! Density and compressibility look-up at bubble point
   call EOSOilDensity        (auxvar%temp,pb,deno      ,ierr,auxvar%table_idx)
   call EOSOilCompressibility(auxvar%temp,pb,cr        ,ierr,auxvar%table_idx)
 endif

!  Correct for undersaturation: correction not yet available for energy
! --------- Correct for undersaturation ---------------------------------------
  crusp=cr*(po-pb)


  if (getDerivs) then
    ! differentiate the correction for undersaturation (see below)
    ! NOTE: might have to adjust this when corrections are introduced for energy
    if (isSat) then
      ! leave density derivs w.r.t. pres, temp, and bubble point as set above since the correction
      ! is zero.
    else
      dcrusp_dpo =  cr !  +  dcr_dpo*(po-pb) but we know dcr_dpo is zero
      dcrusp_dpb = dcr_dpb*(po-pb) - cr
      dcrusp_dt =  dcr_dt*(po-pb)

      cor = 1 + crusp + 0.5*crusp*crusp
      ! dcor/dx = dcrusp/dx + crusp*dcrusp/dx  = (1 + crusp)*dcrusp/dx
      ! so
      one_p_crusp = 1.d0 + crusp
      dcor_dpo = one_p_crusp*dcrusp_dpo
      dcor_dpb = one_p_crusp*dcrusp_dpb
      dcor_dt = one_p_crusp*dcrusp_dt

      auxvar%D_den(oid,dof_op) = deno*dcor_dpo ! +  ddeno_dpo*cor but we know ddeno_dpo is zero
      auxvar%D_den(oid,dof_gsat) = auxvar%D_den(oid,dof_gsat)*cor + deno*dcor_dpb
      auxvar%D_den(oid,dof_temp) = auxvar%D_den(oid,dof_temp)*cor + deno*dcor_dt
    endif

  endif

  !crusp=cr*(po-pb) ! moved above
  deno=deno*(1.0+crusp*(1.0+0.5*crusp))
  auxvar%den(oid)=deno
  auxvar%H  (oid)=ho
  auxvar%U  (oid)=uo
! --------- /Correct for undersaturation --------------------------------------


!------------------------------------------------------------------------------
! /end of Oil phase thermodynamic properties from PVCO
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Correct oil phase molar density and enthalpy for oil composition
! EOSOilDensity returns oil moles/volume in oil phase,
! but really have (1+Rsmolar) times as many moles when dissolved gas included
! 1+Rsmolar=1+xg/xo=(xo+xg)/xo=1/xo.
! Note xo=1/(1+Rsmolar) so cannot be zero for finite looked up Rsmolar
!------------------------------------------------------------------------------

!--------------oil molar density correction------------------------------------
  if (getDerivs) then
    auxvar%D_den(oid,:) = DivRule(auxvar%den(oid),auxvar%D_den(oid,:), &
                                  auxvar%bo%xo,auxvar%bo%D_xo,option%nflowdof)
  endif
  auxvar%den(oid)=auxvar%den(oid)/auxvar%bo%xo

!--------------/oil molar density correction-----------------------------------


! Get oil mass density as (mixture oil molar density).(mixture oil molecular weight)

!--------------oil mass density ------------------------------------------------
  if (getDerivs) then
    worker = auxvar%bo%xo*EOSOilGetFMW()     &
           + auxvar%bo%xg*EOSGasGetFMW()
    D_worker = auxvar%bo%D_xo*EOSOilGetFMW() &
             + auxvar%bo%D_xg*EOSGasGetFMW()
    auxvar%D_den_kg(oid,:) = ProdRule(auxvar%den(oid),auxvar%D_den(oid,:), &
                                      worker,D_worker,option%nflowdof)
  endif
  auxvar%den_kg(oid) = auxvar%den(oid) * ( auxvar%bo%xo*EOSOilGetFMW() &
                                          +auxvar%bo%xg*EOSGasGetFMW() )
!--------------/oil mass density -----------------------------------------------

!--------------H and U scaling -------------------------------------------------
  if (getDerivs) then
    auxvar%D_H(oid,:) = auxvar%D_H(oid,:) * 1.d-6 ! J/kmol -> MJ/kmol
    auxvar%D_U(oid,:) = auxvar%D_U(oid,:) * 1.d-6 ! J/kmol -> MJ/kmol
  endif
  auxvar%H(oid) = auxvar%H(oid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(oid) = auxvar%U(oid) * 1.d-6 ! J/kmol -> MJ/kmol
!--------------/H and U scaling ------------------------------------------------


!------------------------------------------------------------------------------
! Get oil enthalpy/oil mole in oil phase.
! Calculation was for pure oil, but really have total molar enthalpy of:
! oil hydrocarbon enthalpy/oil hydrocarbon mole= xo*(oil enthalpy/oil mole)
!                                               +xg*(gas enthalpy/gas mole)
!------------------------------------------------------------------------------

!--------------H and U correction ---------------------------------------------
  if (getDerivs) then
    auxvar%D_H(oid,:)   = ProdRule(auxvar%bo%xo,auxvar%bo%D_xo,        &
                                 auxvar%H(oid),auxvar%D_H(oid,:),ndof) &
                        +  ProdRule(auxvar%bo%xg,auxvar%bo%D_xg,       &
                                 auxvar%H(gid),auxvar%D_H(gid,:),ndof) 

    auxvar%D_U(oid,:)   = ProdRule(auxvar%bo%xo,auxvar%bo%D_xo,        &
                                 auxvar%U(oid),auxvar%D_U(oid,:),ndof) &
                        +  ProdRule(auxvar%bo%xg,auxvar%bo%D_xg,       &
                                 auxvar%U(gid),auxvar%D_U(gid,:),ndof) 
    

  endif
  auxvar%H(oid) = auxvar%bo%xo*auxvar%H(oid)+auxvar%bo%xg*auxvar%H(gid)
  auxvar%U(oid) = auxvar%bo%xo*auxvar%U(oid)+auxvar%bo%xg*auxvar%U(gid)
!--------------/H and U correction --------------------------------------------

!===============================================================================
! Fluid mobility calculation
!===============================================================================

!-------------------------------------------------------------------------------
!  Water mobility (rel. perm / viscosity)
!-------------------------------------------------------------------------------

  call characteristic_curves%wat_rel_perm_func_owg% &
                RelativePermeability(sw,krw,dkrw_satw,option,auxvar%table_idx)

  if (getDerivs) then
    call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,dps_dt,ierr)

    call EOSWaterViscosity(auxvar%temp, cell_pressure, &
                           wat_sat_pres, dps_dt, visw, &
                           dx_dt,  dx_dcell_pres, ierr)

    ! pressure deriv (dvw_dp) is a a cell pressure derivative:
    D_visc = 0.d0
    D_visc = dx_dcell_pres * D_cell_pres
    D_visc(dof_temp) = D_visc(dof_temp) + dx_dt

  else
    call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,ierr)

    ! use cell_pressure; cell_pressure - psat calculated internally
    call EOSWaterViscosity(auxvar%temp,cell_pressure,wat_sat_pres,visw,ierr)
  endif

  auxvar%mobility(wid) = krw/visw
  if (getDerivs) then
    D_kr = 0.d0
    D_kr(dof_osat) = -dkrw_satw
    if (isSat) D_kr(dof_gsat) = -dkrw_satw
    D_kr(dof_ssat) = -dkrw_satw

    auxvar%D_mobility(wid,:) = DivRule(krw,D_kr,visw,D_visc,ndof)
  endif

!-------------------------------------------------------------------------------
!  Hydrocarbon mobilities
!-------------------------------------------------------------------------------

!--Get and scale the endpoints--------------------------------------------------

  call characteristic_curves%GetOWGCriticalAndConnateSats(swcr,sgcr,dummy, &
                      sowcr,sogcr,swco,option)

  call TL4PScaleCriticals(sgcr,sogcr,fm,so,sv,uoil,uvap, &
                          getDerivs,ndof,dof_osat,dof_gsat,dof_ssat,&
                          D_fm,D_uoil,D_uvap)

  if (.not. isSat) then
    D_uoil(dof_gsat) = 0.d0
    D_uvap(dof_gsat) = 0.d0
  endif

  ! store variables if debugging mode wants to:
  if (getDerivs) then
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%uoil = uoil
      auxvar%tlT%uvap = uvap

      auxvar%tlT%D_uoil = D_uoil
      auxvar%tlT%D_uvap = D_uvap
    endif
  endif


!-------------------------------------------------------------------------------
! Oil mobility (rel. perm / viscosity)
!-------------------------------------------------------------------------------

!--Oil relative permeabilities--------------------------------------------------*/

!--Immiscible lookup (Eclipse 3-phase)

! Get Krow(So) directly from the krow member function of kro
  call characteristic_curves%oil_rel_perm_func_owg% &
              RelPermOW(so,krow,dkro_sato,option,auxvar%table_idx)

! Get Krog(Uo) directly from the krow member function of kro
  call characteristic_curves%oil_rel_perm_func_owg% &
              RelPermOG(uoil,krog,dkrog_uoil,option,auxvar%table_idx)

!--Some simple intermediate derivs:-------------------------------------------
  D_krog = D_uoil*dkrog_uoil
  D_krow = 0.d0
  D_krow(dof_osat) = dkro_sato

  D_sv = 0.d0
  D_sv(dof_ssat) = 1.d0
  if (isSat) D_sv(dof_gsat) = 1.d0

  D_sw = 0.d0
  D_sw(dof_ssat) = -1.d0
  D_sw(dof_osat) = -1.d0
  if (isSat) D_sw(dof_gsat) = -1.d0
!--/Some simple intermediate derivs:------------------------------------------


! Form the Eclipse three-phase Kro expression


!--kroi and derivs:----------------------------------------------------------
  den=sv+sw-swco
  if( den>0.0 ) then
    kroi=(sv*krog+(sw-swco)*krow)/den
  else
    kroi=0.5*(krog+krow)
  endif
  
  if (getDerivs) then
    if (den>0.0) then
      num=(sv*krog+(sw-swco)*krow)
      D_num = ProdRule(sv,D_sv,krog,D_krog,ndof)     &
            + ProdRule(sw-swco,D_sw,krow,D_krow,ndof)
      D_den = 0.d0
      D_den(dof_osat) = -1.d0
      D_kroi = DivRule(num,D_num,den,D_den,ndof)
    else
      D_kroi = 0.5*(D_krog+D_krow)
    endif
  endif

  ! store variables if debugging mode wants to:
  if (getDerivs) then
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%kroi = kroi
      auxvar%tlT%D_kroi = D_kroi

      auxvar%tlT%krog = krog
      auxvar%tlT%D_krog = D_krog

      auxvar%tlT%krow = krow
      auxvar%tlT%D_krow = D_krow
    endif
  endif
!--/kroi and derivs:---------------------------------------------------------


!--Miscible lookup (Krow at Sh=So+Sg+Ss)
  call characteristic_curves%oil_rel_perm_func_owg% &
                RelPermOW(sh,krh,dkrh_sath,option,auxvar%table_idx)

!  Obtain all the miscible limit rel. perms. from Krh, allowing for criticals

!  Get the modified critical oil on gas and vapour

  socrs=fi*sogcr
  svcrs=fi*sgcr

!--sh and krh derivs:---------------------------------------------------------
  if (getDerivs) then
    D_sh(:) = auxvar%D_sat(dof_osat,:) + auxvar%D_sat(dof_gsat,:) + auxvar%D_sat(dof_ssat,:)
    D_krh = dkrh_sath*D_sh
  endif
!--/sh and krh derivs:--------------------------------------------------------

!--krom and derivs:----------------------------------------------------------
! For non-zero Krom, must have So and Sh > Socrs
  if( (sh > (socrs+epss)) .and. (so >= socrs) ) then
    krom=krh*(so-socrs)/(sh-socrs)
  else
    krom=0.0
  endif

 ! recall sh = so + sg + ss
  if (getDerivs) then
    D_krom = 0.d0
    if( (sh > (socrs+epss)) .and. (so >= socrs) ) then
      num = krh*(so-socrs)
      worker = so-socrs
      ! D_socrs = sogcr*D_fi = -sogcr*D_fm so:
      D_worker = sogcr*D_fm + auxvar%D_sat(dof_osat,:)

      D_num = ProdRule(krh,D_krh,worker,D_worker,ndof)

      den = sh - socrs
      D_den = D_sh + sogcr*D_fm

      D_krom = DivRule(num,D_num,den,D_den,ndof)

    else
      D_krom = 0.d0
    endif
  endif

 ! store variables if debugging mode wants to:
  if (getDerivs) then
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%krom = krom
      auxvar%tlT%D_krom = D_krom
   endif
  endif
!--/krom and derivs:---------------------------------------------------------

!--krvm and derivs:----------------------------------------------------------
! For non-zero Krvm, must have Sv and Sh > Svcrs
  if( (sh > (svcrs+epss)) .and. (sv >= svcrs) ) then
    krvm=krh*(sv-svcrs)/(sh-svcrs)
  else
    krvm=0.0
  endif

  if (getDerivs) then
    D_krvm = 0.d0
    if((sh > (svcrs+epss)) .and. (sv >= svcrs)) then
      num = krh*(sv-svcrs)
      den = sh-svcrs
      worker = sv-svcrs
      ! D_svcrs = sgcr*D_fi = -sgcr*D_fm so:
      D_worker = sgcr*D_fm + D_sv
      D_num= ProdRule(krh,D_krh,worker,D_worker,ndof)

      deno= sh - svcrs
      D_den= D_sh + sgcr*D_fm

      D_krvm = DivRule(num,D_num,den,D_den,ndof)
    else
      D_krvm = 0.d0
    endif
  endif
  ! store variables if debugging mode wants to:
  if (getDerivs) then
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%krvm = krvm
      auxvar%tlT%D_krvm = D_krvm
    endif
  endif
!--/krvm and derivs----------------------------------------------------------


! Now split vapour Krvm into Krgm and Krsm pro rata saturations

!--krgm, krsm and derivs:----------------------------------------------------
  if( sv > 0.0 ) then
    krgm=krvm*sg/sv
    krsm=krvm*ss/sv
  else
    krgm=0.0
    krsm=0.0
  endif

   if (getDerivs) then
    D_krgm = 0.d0
    D_krsm = 0.d0
    if( sv > 0.0) then

      num = krvm*sg;
      ! prod rule but we know there's only one nonzero deriv of sg so:
      D_num= D_krvm*sg
      if (isSat) D_num(dof_gsat) = D_num(dof_gsat) + krvm
      D_krgm = DivRule(num,D_num,sv,D_sv,ndof)

      num= krvm*ss;
      ! prod rule but we know there's only one nonzero deriv of ss so:
      D_num = D_krvm*ss; D_num(dof_ssat) = D_num(dof_ssat) + krvm
      ! worker and D_worker unchanged
      D_krsm = DivRule(num,D_num,sv,D_sv,ndof)
    endif
  endif
  ! store variables if debugging mode wants to:
  if (getDerivs) then
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%krgm = krgm
      auxvar%tlT%D_krgm = D_krgm

      auxvar%tlT%krsm = krsm
      auxvar%tlT%D_krsm = D_krsm
    endif
  endif
!--/krgm, krsm and derivs---------------------------------------------------


!--Oil viscosities--------------------------------------------------------------*/

!--If PVCO defined in EOS OIL, the viscosities are extracted via table lookup--

  ! Viscosity and viscosibility look-up at bubble point
  if (getDerivs) then
    call EOSOilViscosity(auxvar%temp,pb,auxvar%den(oid),viso,dvo_dt,dvo_dpb,ierr,auxvar%table_idx)
    call EOSOilViscosibility(auxvar%temp,pb,cvisc,dcvisc_dt,dcvisc_dpb,ierr, &
                             auxvar%table_idx)
  else
    call EOSOilViscosity    (auxvar%temp,pb, &
                             auxvar%den(oid), viso, ierr,auxvar%table_idx)
    call EOSOilViscosibility(auxvar%temp,pb,cvisc,ierr,auxvar%table_idx)
  endif


!----------Correct oil viscosity-----------------------------------------------
  cvusp=cvisc*(po-pb)

  if (getDerivs) then
    ! get derivatives of corrected viso

    ! if saturated:
      ! 1) cvusp = 0; cor = 1 constants
      ! 2) pb is really oil pressure so pb derivs equal to po derivs
    if (.NOT. isSat) then
      ! 1) cvusp, cor now vary
      ! 2) pb derivs diferent from po derivs
      ! 3) have to apply corrector visc = cor*visc, so apply corresponding correctors
      !    to derivs too

      dcvusp_dpo = cvisc ! cvisc is independent of po here
      dcvusp_dpb = dcvisc_dpb*(po-pb) - cvisc
      dcvusp_dt = dcvisc_dt*(po-pb)

      cor = (1.0+cvusp*(1.0+0.5*cvusp))
      ! dcor/dx = dcvusp/dx + cvusp*dcvusp so:
      one_p_cvusp = 1.d0 + cvusp
      dcor_dpo = dcvusp_dpo*one_p_cvusp
      dcor_dpb = dcvusp_dpb*one_p_cvusp
      dcor_dt = dcvusp_dt*one_p_cvusp

      dvo_dp = viso*dcor_dpo ! viso is independent of po here
      dvo_dpb = viso*dcor_dpb + dvo_dpb*cor
      dvo_dt = viso*dcor_dt + dvo_dt*cor

    endif
    D_viso = 0.d0
    if (isSat) then
      D_viso(dof_op) = dvo_dpb
    else
      D_viso(dof_op) = dvo_dp
      D_viso(dof_gsat) = dvo_dpb
    endif
    ! no ssat deriv
    D_viso(dof_temp) = dvo_dt

  endif

  !cvusp=cvisc*(po-pb) !!! moved above
  viso=viso*(1.0+cvusp*(1.0+0.5*cvusp))
!----------/Correct oil viscosity-----------------------------------------------

!-------------------------------------------------------------------------------
!  Vapour (gas and solvent) mobilities (rel. perm / viscosity)
!-------------------------------------------------------------------------------

  call characteristic_curves%gas_rel_perm_func_owg% &
               RelativePermeability(sv,krvi,dkrv_satv,option,auxvar%table_idx)


!  Obtain all the immiscible limit rel. perms. from Krvi

!-0=-krgi, krsi, and derivs:-----------------------------------------------------
  if( sv>0.0 ) then
    krgi=sg*krvi/sv
    krsi=ss*krvi/sv
  else
    krgi=0.0
    krsi=0.0
  endif

    if (getDerivs) then
    ! existence of dkrv_satv implies:
    ! dkrvi_sg = dkrv_satv
    ! dkrvi_ss = dkrv_satv
    ! b/c sv = ss + sg
    D_krgi = 0.d0
    D_krsi = 0.d0
    if( sv>0.0 ) then
      ! (would also be straightforward to replace this with a proddivrule() call)
      if (isSat) D_krgi(dof_gsat) = sg*dkrv_satv/sv + ss*krvi/sv/sv
      D_krgi(dof_ssat) =            sg*dkrv_satv/sv - sg*krvi/sv/sv

      ! (would also be straightforward to replace this with a proddivrule() call)
      if (isSat) D_krsi(dof_gsat) = ss*dkrv_satv/sv - ss*krvi/sv/sv
      D_krsi(dof_ssat) =            ss*dkrv_satv/sv + sg*krvi/sv/sv
    else
      ! nothing, or maybe constant slope?
    endif
  endif

  ! store variables if debugging mode wants to:
  if (getDerivs) then
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%krgi = krgi
      auxvar%tlT%D_krgi = D_krgi

      auxvar%tlT%krsi = krsi
      auxvar%tlT%D_krsi = D_krsi
    endif
  endif
!--/krgi, krsi, and derivs:-----------------------------------------------------

!--If PVDG defined in EOS GAS, the viscosities are extracted via table lookup--
  if (getDerivs) then
      D_visg = 0.d0
      call  EOSGasViscosity(auxvar%temp,auxvar%pres(gid),auxvar%pres(gid),&
                            auxvar%den(gid), &
                            auxvar%D_den(gid,dof_temp), auxvar%D_den(gid,dof_temp), auxvar%D_den(gid,dof_op), &
                            0.d0,1.d0, &                      ! dPcomp_dT, dPcomp_dPgas
                            visg, dx_dt, dummy, dx_dcell_pres, ierr, &
                            auxvar%table_idx)

    ! pressure deriv is a a gas pressure derivative:
    D_visg = dx_dcell_pres * auxvar%D_pres(gid,:)
    D_visg(dof_temp) = D_visg(dof_temp) + dx_dt

  else
    call EOSGasViscosity(auxvar%temp,auxvar%pres(gid), &
                         auxvar%pres(gid),auxvar%den(gid),visg,ierr,&
                         auxvar%table_idx)
  endif

!--If PVDS defined in EOS SLV, the viscosities are extracted via table lookup--

  if (getDerivs) then
    call EOSSlvViscosity(auxvar%temp,auxvar%pres(sid), &
                         viss,dx_dt,dx_dcell_pres,ierr,&
                         auxvar%table_idx)

    ! pressure deriv is a a solvent pressure derivative:
    !D_viss = dx_dcell_pres * D_cell_pres
    D_viss = dx_dcell_pres * auxvar%D_pres(sid,:)
    D_viss(dof_temp) = D_viss(dof_temp) + dx_dt
  else
    call EOSSlvViscosity(auxvar%temp,auxvar%pres(sid), &
                         viss,ierr,&
                         auxvar%table_idx)
  endif

  ! store variables if debugging mode wants to:
  if (getDerivs) then
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%viso = viso
      auxvar%tlT%visg = visg
      auxvar%tlT%viss = viss

      auxvar%tlT%D_viso = D_viso
      auxvar%tlT%D_visg = D_visg
      auxvar%tlT%D_visS = D_viss
    endif
  endif

!-------------------------------------------------------------------------------
!  Form the Todd-Longstaff rel perms
!-------------------------------------------------------------------------------

  call TL4PRelativePermeabilities( so,sg,sw,ss,sv,sh,fm,      &
                                 swcr,sgcr,sowcr,sogcr,swco,  &
                                 kroi,krgi,krsi,              &
                                 krom,krgm,krsm,              &
                                 kro ,krg ,krs,               &
                                 getDerivs,ndof,              &
                                 D_fm,                        & ! fm derivs (in)
                                 D_kroi,D_krgi,D_krsi,        & ! imiscible k derivs (in)
                                 D_krom,D_krgm,D_krsm,        & ! miskible k derivs  (in)
                                 D_kro ,D_krg ,D_krs   )        ! true k derivs     (out)

  ! store variables if debugging mode wants to:
  if (getDerivs) then
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%krotl = kro
      auxvar%tlT%krgtl = krg
      auxvar%tlT%krstl = krs
      auxvar%tlT%D_krotl = D_kro
      auxvar%tlT%D_krgtl = D_krg
      auxvar%tlT%D_krstl = D_krs
    endif
  endif

!--Form the Todd-Longstaff viscosities-----------------------------------------

  if (getDerivs) then
    ! can cause problems if we pass in these
    ! auxvar members directly if they are not allocated
    ! (i.e. when running numerical) and in parallel, crashes
    ! some tl np2 regtests.
    ! (alternative: optional arguments)
    D_so = auxvar%D_sat(oid,:)
    D_sg = auxvar%D_sat(gid,:)
    D_ss = auxvar%D_sat(sid,:)
  else
    ! they shouldn't be used in this case but for robustness:
    D_so = 0.d0
    D_sg = 0.d0
    D_ss = 0.d0
  endif

  call TL4PViscosity(so,sg,ss,viso,visg,viss,visotl,visgtl,visstl, &
                     getDerivs,ndof,                               &
                     D_so,D_sg,D_ss,                               & 
                     dof_osat,dof_gsat,dof_ssat,                   &
                     D_viso,D_visg,D_viss,                         & ! visc derivs in
                     D_visotl,D_visgtl,D_visstl,                   & ! tl visc derivs out
                     isSat,                                        &
                     denos,D_denos,dengs,D_dengs,denogs,D_denogs,  &
                     denos_pre,D_denos_pre)


  ! store variables if debugging mode wants to:
  if (getDerivs) then
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%viscotl = visotl
      auxvar%tlT%viscgtl = visgtl
      auxvar%tlT%viscstl = visstl

      auxvar%tlT%denos = denos
      auxvar%tlT%dengs = dengs
      auxvar%tlt%denogs = denogs

      auxvar%tlT%denos_pre = denos_pre
      auxvar%tlT%D_denos_pre = D_denos_pre

      auxvar%tlT%D_viscotl = D_visotl
      auxvar%tlT%D_viscgtl = D_visgtl
      auxvar%tlT%D_viscstl = D_visstl

      auxvar%tlT%D_denos = D_denos
      auxvar%tlT%D_dengs = D_dengs
      auxvar%tlT%D_denogs = D_denogs
    endif
  endif

!--Form the Todd-Longstaff densities (for Darcy flow only)---------------------

  deno=auxvar%den_kg(oid)
  deng=auxvar%den_kg(gid)
  dens=auxvar%den_kg(sid)

  if (getDerivs) then
    ! can cause problems if we pass in these
    ! auxvar members directly if they are not allocated
    ! (i.e. when running numerical) and in parallel, crashes
    ! some tl np2 regtests
    ! (alternative: optional arguments)
    D_deno = auxvar%D_den_kg(oid,:)
    D_deng = auxvar%D_den_kg(gid,:)
    D_dens = auxvar%D_den_kg(sid,:)
  else
    ! they shouldn't be used in this case but might as well
    D_deno = 0.d0
    D_deng = 0.d0
    D_dens = 0.d0
  endif

  call TL4PDensity( so,sg,ss,                             &
                    viso,visg,viss,visotl,visgtl,visstl,  &
                    deno,deng,dens,denotl,dengtl,denstl , &
                    getDerivs, ndof ,                     &
                    D_so,D_sg,D_ss,                       &
                    dof_osat,dof_gsat,dof_ssat,           &
                    D_deno,D_deng,D_dens,                 &
                    D_viso,D_visg,D_viss,                 &
                    D_visotl,D_visgtl,D_visstl,           &
                    D_denotl,D_dengtl,D_denstl,           &
                    isSat,                                &
                    denog,D_denog)

  ! store variables if debugging mode wants to:
  if (getDerivs) then
    if (auxvar%has_TL_test_object) then
      auxvar%tlT%denotl = denotl
      auxvar%tlT%dengtl = dengtl
      auxvar%tlT%denstl = denstl

      auxvar%tlT%denog = denog
      auxvar%tlT%D_denog = D_denog

      auxvar%tlT%D_denotl = D_denotl
      auxvar%tlT%D_dengtl = D_dengtl
      auxvar%tlT%D_denstl = D_denstl
    endif
  endif

!------------------------------------------------------------------------------
!  Make and load mobilities
!------------------------------------------------------------------------------


  if (getDerivs) then
    auxvar%D_mobility(oid,:) = DivRule(kro,D_kro,visotl,D_visotl,ndof)
    auxvar%D_mobility(gid,:) = DivRule(krg,D_krg,visgtl,D_visgtl,ndof)
    auxvar%D_mobility(sid,:) = DivRule(krs,D_krs,visstl,D_visstl,ndof)
  endif

  if (mobSanityCheck) then
    if (kro > 0.0) auxvar%mobility(oid) = kro/visotl
    if (krg > 0.0) auxvar%mobility(gid) = krg/visgtl
    if (krs > 0.0) auxvar%mobility(sid) = krs/visstl
  else
    auxvar%mobility(oid) = kro/visotl
    auxvar%mobility(gid) = krg/visgtl
    auxvar%mobility(sid) = krs/visstl
  endif

#if 0
!!! we may wish to have a robust optional error check for things like these, it can
!!! happen surprisingly often if saturations become marginal
  if (isnan(auxvar%mobility(oid))) then
    print *, "nan mob oil"
  endif
  if (isnan(auxvar%mobility(gid))) then
    print *, "nan mob gas"
  endif
  if (isnan(auxvar%mobility(sid))) then
    print *, "nan mob slv"
  endif
#endif

!------------------------------------------------------------------------------
!  Load densities
!------------------------------------------------------------------------------

  auxvar%den_kg(oid)=denotl
  auxvar%den_kg(gid)=dengtl
  auxvar%den_kg(sid)=denstl

  if (getDerivs) then
    auxvar%D_den_kg(oid,:)=D_denotl
    auxvar%D_den_kg(gid,:)=D_dengtl
    auxvar%D_den_kg(sid,:)=D_denstl
  endif

end subroutine TL4PAuxVarCompute


! ************************************************************************** !

subroutine TOWGTLAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, &
                               characteristic_curves,natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell for TOWGTL
  !
  ! Author: Paolo Orsini (OGS) - David Ponting (OGS)
  ! Date: 07/06/17
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Oil_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use Derivatives_utilities_module

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar ! passing this for salt conc.
                                            ! not currenty used
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  PetscInt :: natural_id !only for debugging/print out - currently not used

  PetscInt :: wid, oid, gid
  PetscReal :: cell_pressure, wat_sat_pres
  !PetscReal :: krl, dkrl_Se
  !PetscReal :: krh, dkrh_Se
  PetscReal :: sath
  PetscReal :: krw,dkrw_sato,dkrw_satg,dkrw_satw
  PetscReal :: krh,dkrh_sato,dkrh_satg, dkrh_sath
  PetscReal :: viso
  PetscReal :: visg
  !PetscReal :: visl
  PetscReal :: visw
  PetscReal :: sat_water
  PetscReal :: dummy,dummy2,deno,deng
  !PetscReal :: Uoil_J_kg, Hoil_J_kg
  PetscErrorCode :: ierr
  PetscReal :: krotl,krgtl,viscotl,viscgtl,denotl,dengtl

  PetscBool :: getDerivs
  PetscInt :: dof_op,dof_osat,dof_gsat,dof_temp,cploc,ndof
  PetscReal,dimension(1:option%nflowdof) :: D_cell_pres,D_kr,D_visc,D_visco,D_viscg
  PetscReal,dimension(1:option%nflowdof) :: D_deno,D_deng,D_sath
  PetscReal,dimension(1:option%nflowdof) :: D_krotl,D_krgtl,D_viscotl,D_viscgtl,D_denotl,D_dengtl
  PetscReal :: d_x_d_cp,d_x_d_T,d_ps_d_T

  dof_op = TOWG_OIL_PRESSURE_DOF
  dof_osat = TOWG_OIL_SATURATION_DOF
  dof_gsat = TOWG_GAS_SATURATION_3PH_DOF
  dof_temp = towg_energy_dof
  ndof = option%nflowdof



  if (towg_analytical_derivatives) then
    if (.NOT. auxvar%has_derivs) then
      ! how did this happen?
      option%io_buffer = 'towg tl auxvars: towg_analytical_derivatives is true, &
                          but auxvar%has_derivs is false, should both be true. &
                          How did this happen?'
      call printErrMsg(option)
    endif

    auxvar%D_pres = 0.d0
    auxvar%D_sat = 0.d0
    auxvar%D_pc = 0.d0
    auxvar%D_den = 0.d0
    auxvar%D_den_kg = 0.d0
    auxvar%D_mobility = 0.d0
    auxvar%D_por = 0.d0

    auxvar%D_H = 0.d0
    auxvar%D_U = 0.d0

    auxvar%tl%D_den_oil_eff_kg(:) = 0.d0
    auxvar%tl%D_den_gas_eff_kg(:) = 0.d0


    getDerivs = PETSC_TRUE
  else 
    getDerivs = PETSC_FALSE
  endif


!--Initialise-----------------------------------------------------------------

  krotl  =0.0
  krgtl  =0.0

  viscotl=0.0
  viscgtl=0.0

  denotl =0.0
  dengtl =0.0

!--Get phase pointers for water,oil and gas-----------------------------------

  ! from SubsurfaceSetFlowMode
  ! option%liquid_phase = 1           ! liquid_phase - water phase idx
  ! option%oil_phase = 2              ! oil_phase
  ! option%gas_phase = 3              ! gas_phase

  wid = option%liquid_phase
  oid = option%oil_phase
  gid = option%gas_phase

  auxvar%effective_porosity = 0.d0
  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%pc = 0.d0
  auxvar%den = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%mobility = 0.d0
  auxvar%H = 0.d0
  auxvar%U = 0.d0
  auxvar%temp = 0.0d0

!--Extract into auxvars from the solution variables----------------------------

  auxvar%pres(oid) = x(TOWG_OIL_PRESSURE_DOF)
  auxvar%sat(oid) = x(TOWG_OIL_SATURATION_DOF)
  auxvar%sat(gid) = x(TOWG_GAS_SATURATION_3PH_DOF)
  auxvar%temp = x(towg_energy_dof)

  auxvar%sat(wid) = 1.d0 - auxvar%sat(oid) - auxvar%sat(gid)

  ! trivial saturation derivatives:
  if (getDerivs) then
    auxvar%D_sat(oid,dof_osat) =  1.d0 ! diff oil sat by oil sat
    auxvar%D_sat(wid,dof_osat) = -1.d0 ! diff liquid sat by gas sat
    auxvar%D_sat(wid,dof_gsat) = -1.d0 ! diff liquid sat by gas sat
    auxvar%D_sat(gid,dof_gsat) =  1.d0 ! diff gas sat by gas sat
  endif


!--Set up phase pressures and capillary pressures-----------------------------

  ! In Todd-Longstaff mode assume no cap. pressure between oil and gas: pcog=0
  ! The cap. pressure is hydrocarbon/water (pchw) as function of water satn.

  auxvar%pc(oid) = 0.0d0

  sat_water = auxvar%sat(wid)

  call characteristic_curves%oil_wat_sat_func% &
              CapillaryPressure(sat_water,auxvar%pc(wid),dummy, &
                                option,auxvar%table_idx)

  auxvar%pres(gid) = auxvar%pres(oid)
  auxvar%pres(wid) = auxvar%pres(oid) - auxvar%pc(wid)

  if (getDerivs) then 
    auxvar%D_pres(oid,dof_op) = 1.d0
    auxvar%D_pres(gid,dof_op) = 1.d0
    auxvar%D_pres(wid,dof_op) = 1.d0 

    ! dummy is d cp / d s_w 
    ! sw = 1 - s_o - s_g s
    ! so d cp / d s_o = d cp / d s_w * d sw / d s_o = -1*d pc / d s_w
    ! so:
    !auxvar%D_pc(oid,dof_osat) = dummy
    !auxvar%D_pc(oid,dof_gsat) = dummy
    auxvar%D_pc(wid,dof_osat) = -dummy
    auxvar%D_pc(wid,dof_gsat) = -dummy
    auxvar%D_pres(wid,dof_osat) = -auxvar%D_pc(wid,dof_osat)
    auxvar%D_pres(wid,dof_gsat) = -auxvar%D_pc(wid,dof_gsat)

  endif 

  cell_pressure = max(auxvar%pres(wid),auxvar%pres(oid),auxvar%pres(gid))

  if (getDerivs) then
    ! a bit ugly but works:
    cploc = wid
    if (auxvar%pres(oid) > auxvar%pres(cploc)) then
      cploc = oid
    endif
    if (auxvar%pres(gid) > auxvar%pres(cploc)) then
      cploc = gid
    endif
    D_cell_pres = auxvar%D_pres(cploc,:)
  endif

!--Calculate effective porosity as a function of pressure---------------------

  if (option%iflag /= TOWG_UPDATE_FOR_BOUNDARY) then
    auxvar%effective_porosity = material_auxvar%porosity_base

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                auxvar%effective_porosity,dummy)
      if (getDerivs) then
        ! dummy here is d por / d cell pres so
        auxvar%D_por = dummy*D_cell_pres
      endif
    endif
    if (option%iflag /= TOWG_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = auxvar%effective_porosity
    endif
  endif

!--Update thermodynamic properties (density, enphalpy..) for all phases-------

  ! Liquid phase thermodynamic properties
  ! using cell_pressure (which is the max press)? or %pres(wid)?

! WATER THERMODYNAMIC PROPERTIES --------------------------------------------
  if (getDerivs) then
    call EOSWaterDensity(auxvar%temp,cell_pressure, &
                      auxvar%den_kg(wid),auxvar%den(wid), &
                      d_x_d_cp, &
                      d_x_d_T, ierr)

    ! den = den( cell_pressure, temp)
    ! chain rule for cell pressure part:
    auxvar%D_den(wid,:) = d_x_d_cp*D_cell_pres(:)
    ! then additional term for temperature derivative:
    auxvar%D_den(wid,dof_temp) = auxvar%D_den(wid,dof_temp) + d_x_d_T


    call EOSWaterEnthalpy(auxvar%temp, &
                       cell_pressure, &
                       auxvar%H(wid), &
                       d_x_d_cp, &
                       d_x_d_T, &
                       ierr)

    ! H = H H(  cell_pressure, temp)
    ! chain rule for cell pressure part:
    auxvar%D_H(wid,:) = d_x_d_cp*D_cell_pres(:)
    ! then additional term for temperature derivative:
    auxvar%D_H(wid,dof_temp) = auxvar%D_H(wid,dof_temp) + d_x_d_T

    ! derivatives corresponding to scaling of H then computation of U(wid) below:
    auxvar%D_H(wid,:) = auxvar%D_H(wid,:) * 1.d-6 ! J/kmol -> MJ/kmol
    auxvar%D_U(wid,:) = auxvar%D_H(wid,:)                                          &
                      - 1.d-6                                                      &
                      * DivRule(cell_pressure,D_cell_pres,                         &
                                auxvar%den(wid),auxvar%D_den(wid,:),option%nflowdof )

    ! don't forget den kg derivatives:
    auxvar%den_kg(wid) = auxvar%den(wid) * FMWH2O 

  else
    call EOSWaterDensity(auxvar%temp,cell_pressure, &
                         auxvar%den_kg(wid),auxvar%den(wid),ierr)
    call EOSWaterEnthalpy(auxvar%temp,cell_pressure,auxvar%H(wid),ierr)
  endif
  auxvar%H(wid) = auxvar%H(wid) * 1.d-6 ! J/kmol -> MJ/kmol

  ! MJ/kmol comp                  ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
  auxvar%U(wid) = auxvar%H(wid) - (cell_pressure / auxvar%den(wid) * 1.d-6)

! END OF WATER THERMODYNAMIC PROPERTIES -----------------------------------------
  

  ! ADD HERE BRINE dependency. Two options (see mphase)
  ! - salinity constant in space and time (passed in option%option%m_nacl)
  ! - salt can be transported by RT (sequential coupling) and passed
  !   and passed with global_auxvar%m_nacl


! OIL THERMODYNAMIC PROPERTIES ------------------------------------------------------

  if (getDerivs) then
    call EOSOilDensityEnergy(auxvar%temp,auxvar%pres(oid),auxvar%den(oid), &
                             auxvar%D_den(oid,dof_temp),auxvar%D_den(oid,dof_op), &
                             auxvar%H(oid),auxvar%D_H(oid,dof_temp),auxvar%D_H(oid,dof_op), &
                             auxvar%U(oid),auxvar%D_U(oid,dof_temp),auxvar%D_U(oid,dof_op), &
                             ierr,auxvar%table_idx) 

    auxvar%D_den_kg(oid,:) = auxvar%D_den(oid,:) * EOSOilGetFMW()

    auxvar%D_H(oid,:) = auxvar%D_H(oid,:) * 1.d-6 ! J/kmol -> MJ/kmol
    auxvar%D_U(oid,:) = auxvar%D_U(oid,:) * 1.d-6 ! J/kmol -> MJ/kmol
  else
    call EOSOilDensityEnergy(auxvar%temp,auxvar%pres(oid),&
                             auxvar%den(oid),auxvar%H(oid), &
                             auxvar%U(oid),ierr,auxvar%table_idx)
  endif

  auxvar%den_kg(oid) = auxvar%den(oid) * EOSOilGetFMW()

  auxvar%H(oid) = auxvar%H(oid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(oid) = auxvar%U(oid) * 1.d-6 ! J/kmol -> MJ/kmol

! END OF OIL THERMODYNAMIC PROPERTIES -----------------------------------------------

! GAS THERMODYNAMIC PROPERTIES ------------------------------------------------------

  !compute gas properties (default is air - but methane can be set up)
  if (getDerivs) then
    call EOSGasDensityEnergy(auxvar%temp,auxvar%pres(gid),auxvar%den(gid),auxvar%D_den(gid,dof_temp),auxvar%D_den(gid,dof_op), &
                                     auxvar%H(gid),auxvar%D_H(gid,dof_temp),auxvar%D_H(gid,dof_op),auxvar%U(gid),&
                                     auxvar%D_U(gid,dof_temp),&
                                     auxvar%D_U(gid,dof_op),ierr,auxvar%table_idx)
   ! Note some implicit chain rule here. The pressure derivatives this gives are with respect to gas pressure, not cell or oil pressure.
   ! Referencing above, we see that 
   ! "  auxvar%pres(gid) = auxvar%pres(oid)  "
   ! is the only dependence - on capilliary pres complications, gas pressure is just oil pressure. So we save some lines and just save
   ! the gas presure derivatives as oil pressure derivatives, since the correct chain rule is just to multpily by one.
   ! NOTE that if someone decides to change the cap pressure behaviour so that gas pressure is more complex, then need to revise this.

    auxvar%D_den_kg(gid,:) = auxvar%D_den(gid,:) * EOSGasGetFMW()
    auxvar%D_H(gid,:) = auxvar%D_H(gid,:) * 1.d-6 ! J/kmol -> MJ/kmol
    auxvar%D_U(gid,:) = auxvar%D_U(gid,:) * 1.d-6 ! J/kmol -> MJ/kmol

  else
    call EOSGasDensityEnergy(auxvar%temp,auxvar%pres(gid),auxvar%den(gid), &
                        auxvar%H(gid),auxvar%U(gid),ierr,auxvar%table_idx)
  endif

  auxvar%den_kg(gid) = auxvar%den(gid) * EOSGasGetFMW()
  auxvar%H(gid) = auxvar%H(gid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(gid) = auxvar%U(gid) * 1.d-6 ! J/kmol -> MJ/kmol

! END OF GAS THERMODYNAMIC PROPERTIES -----------------------------------------------

!--Set up water phase rel. perm. and unmixed viscosities-----------------------

! ------ WATER MOBILITY -----------------------------------------------
  ! compute water relative permability Krw(Sw)
  call characteristic_curves%wat_rel_perm_func_owg% &
                      RelativePermeability(auxvar%sat(wid),krw,dkrw_satw, &
                                                option,auxvar%table_idx)


  if (getDerivs) then
    call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,d_ps_d_T,ierr)
    ! use cell_pressure; cell_pressure - psat calculated internally
    call EOSWaterViscosity(auxvar%temp,cell_pressure,wat_sat_pres,d_ps_d_T,visw,d_x_d_T,d_x_d_cp,ierr)
    D_visc = d_x_d_T*D_cell_pres
    D_visc(dof_temp) = D_visc(dof_temp) + d_x_d_T
  else
    call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,ierr)
    ! use cell_pressure; cell_pressure - psat calculated internally
    call EOSWaterViscosity(auxvar%temp,cell_pressure,wat_sat_pres,visw,ierr)
  endif

  auxvar%mobility(wid) = krw/visw

  if (getDerivs) then
    D_kr = 0.d0
    D_kr(dof_osat) = -dkrw_satw
    D_kr(dof_gsat) = -dkrw_satw
    auxvar%D_mobility(wid,:) = DivRule(krw,D_kr,visw,D_visc,ndof)
  endif

! END OF WATER MOBILITY -----------------------------------------------

!--Set up hydrocarbon phase rel. perm. and unmixed viscosities----------------

  ! In Todd-Longstaff case, look up the hydrocarbon rel perm using the
  ! hydrocarbon saturation (Sh=So+Sg=1-Sw), with water saturation as argument

  sath = auxvar%sat(oid) + auxvar%sat(gid)
  dummy = 0.0

  if (getDerivs) then
    D_sath = auxvar%D_sat(oid,:) + auxvar%D_sat(gid,:)
  endif

  call characteristic_curves%ow_rel_perm_func_owg% &
               RelativePermeability(sath,krh,dkrh_sath,option,auxvar%table_idx)
  if (getDerivs) then
    D_kr = 0.d0
    D_kr = dkrh_sath*D_sath
  endif

  ! Oil viscosity
  if (getDerivs) then
    D_visco = 0.d0
    call EOSOilViscosity(auxvar%temp,auxvar%pres(oid), &
                         auxvar%den(oid), viso, D_visco(dof_temp),D_visco(dof_op), &
                         ierr,auxvar%table_idx)
  else
    call EOSOilViscosity(auxvar%temp,auxvar%pres(oid), &
                         auxvar%den(oid), viso, ierr,auxvar%table_idx)
  endif

  ! Gas viscosity : currently only viscosity model for air or constant value
  if (getDerivs) then
    D_viscg = 0.d0
    call EOSGasViscosity(auxvar%temp,auxvar%pres(gid),auxvar%pres(gid),auxvar%den(gid), &
                         auxvar%D_den(gid,dof_temp),auxvar%D_den(gid,dof_op),auxvar%D_den(gid,dof_op), &
                         auxvar%D_pres(gid,dof_temp),1.d0, &
                         visg,D_viscg(dof_temp),D_viscg(dof_op),D_viscg(dof_op), &
                         ierr,auxvar%table_idx)
  else
    call EOSGasViscosity(auxvar%temp,auxvar%pres(gid),auxvar%pres(gid), &
                         auxvar%den(gid),visg,ierr,auxvar%table_idx)
  endif


!--Set up the Todd-Longstaff rel perms, viscosities and densities-------------

  deno=auxvar%den_kg(oid)
  deng=auxvar%den_kg(gid)
  if (getDerivs) then
    D_deno=auxvar%D_den_kg(oid,:)
    D_deng=auxvar%D_den_kg(gid,:)
  endif

  if (getDerivs) then
    call vToddLongstaff( oid,gid,krh,viso,visg,deno,deng,auxvar    &
                        ,krotl,krgtl,viscotl,viscgtl,denotl,dengtl,&
                        PETSC_TRUE,ndof,option,                    &
                        D_visco,D_viscg,D_kr,                      &
                        D_krotl,D_krgtl,D_viscotl,D_viscgtl,       &
                        D_denotl,D_dengtl                           )
  else
    call vToddLongstaff( oid,gid,krh,viso,visg,deno,deng,auxvar &
                        ,krotl,krgtl,viscotl,viscgtl,denotl,dengtl,&
                        PETSC_FALSE,ndof,option)
  endif

  if (auxvar%hastl3p_test_object) then
    auxvar%tl3TEST%viscotl = viscotl
    auxvar%tl3TEST%D_viscotl = D_viscotl

    auxvar%tl3TEST%viscgtl = viscgtl
    auxvar%tl3TEST%D_viscgtl = D_viscgtl

    auxvar%tl3TEST%denotl = denotl
    auxvar%tl3TEST%D_denotl = D_denotl

    auxvar%tl3TEST%dengtl = dengtl
    auxvar%tl3TEST%D_dengtl = D_dengtl

    auxvar%tl3TEST%krotl = krotl
    auxvar%tl3TEST%D_krotl = D_krotl

    auxvar%tl3TEST%krgtl = krgtl
    auxvar%tl3TEST%D_krgtl = D_krgtl

    auxvar%tl3TEST%krh = krh
    auxvar%tl3TEST%D_krh = D_kr
  endif

!--Calculate and store oil and gas mobilities---------------------------------

  if( viscotl>0.0 ) then
    auxvar%mobility(oid) = krotl/viscotl
    if (getDerivs) then
      auxvar%D_mobility(oid,:) = DivRule(krotl,D_krotl,viscotl,D_viscotl,ndof)
    endif
  else
    auxvar%mobility(oid) = 0.0
    if (getDerivs) then
      auxvar%D_mobility(oid,:) = 0.d0
    endif
  endif

  if( viscgtl>0.0 ) then
    auxvar%mobility(gid) = krgtl/viscgtl
    if (getDerivs) then
      auxvar%D_mobility(gid,:) = DivRule(krgtl,D_krgtl,viscgtl,D_viscgtl,ndof)
    endif
  else
    auxvar%mobility(gid) = 0.0
    if (getDerivs) then
      auxvar%D_mobility(gid,:) = 0.d0
    endif
  endif

!--Store mixed density values-------------------------------------------------

  auxvar%tl%den_oil_eff_kg = denotl
  auxvar%tl%den_gas_eff_kg = dengtl
  if (getDerivs) then
    auxvar%tl%D_den_oil_eff_kg = D_denotl
    auxvar%tl%D_den_gas_eff_kg = D_dengtl
  endif

end subroutine TOWGTLAuxVarCompute

!==============================================================================

subroutine vToddLongstaff( oid,gid,krh,visco,viscg,deno,deng,auxvar   &
                          ,krotl,krgtl,viscotl,viscgtl,denotl,dengtl  &
                          ,getDerivs,ndof,option                      &
                          ,D_visco,D_viscg,D_krh                      &
                          ,D_krotl,D_krgtl,D_viscotl,D_viscgtl        &
                          ,D_denotl,D_dengtl                           )

!------------------------------------------------------------------------------
! Routine to calculate the Todd-Longstaff mobilities
! Based on the reference TL:
! 'The Development, Testing and Application of a Numerical Simulator for
!  Predicting Miscible Flood Performance', M.R.Todd and W.J.Longstaff,
!  Journal Pet. Tech., 1972
!  SPE 3484, downloadable from OnePetro on www.spe.org
!------------------------------------------------------------------------------
!
! input  : oid,gid         : Pointers to oil and gas saturations
! input  : krh             : krh contains hydrocarbon rel. perm. on input.
! input  : visco,viscg     : Unmixed phase viscosity values on input
! input  : deno ,deng      : Unmixed phase density   values on input
! input  : auxvar          : Used to obtain saturation values
! output : krh             : Contains oil and gas rel. perms. on output
! output : viscotl,viscgtl : Mixed phase viscosity values on output
! output : denotl ,dengtl  : Mixed phase density   values on output
!
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Jul 2017
!------------------------------------------------------------------------------

  use Derivatives_utilities_module
  use Option_module
  implicit none

  PetscInt ,intent(in )   ::oid,gid
  PetscReal,intent(in )   ::krh
  PetscReal,intent(in )   ::visco,viscg
  PetscReal,intent(in )   ::deno,deng
  PetscBool,intent(in)    :: getDerivs
  PetscInt,intent(in)     :: ndof

  PetscReal,dimension(1:ndof),optional,intent(in ):: D_visco,D_viscg,D_krh
  type(option_type) :: option
  PetscReal,dimension(1:ndof),optional,intent(out):: D_krotl,D_krgtl
  PetscReal,dimension(1:ndof),optional,intent(out):: D_viscotl,D_viscgtl
  PetscReal,dimension(1:ndof),optional,intent(out):: D_denotl,D_dengtl

  class(auxvar_towg_type) :: auxvar
  PetscReal,intent(out)   ::krotl,krgtl
  PetscReal,intent(out)   ::viscotl,viscgtl
  PetscReal,intent(out)   ::denotl,dengtl

  PetscReal::so,sg,den,deninv,fo,fg
  PetscReal,dimension(1:ndof) :: D_worker,D_fo,D_fg
  PetscReal,dimension(1:ndof):: D_deno,D_deng 


  if (getDerivs .AND. &
      .NOT.(present(D_visco).AND.present(D_viscg).AND.present(D_krh).AND.present(D_krotl).AND.present(D_krgtl).AND. &
            present(D_viscotl).AND.present(D_viscgtl).AND.present(D_denotl).AND.present(D_dengtl)                    )) then
      option%io_buffer = 'vToddLongstaff(): attempting to get analytical &
                          derivatives but not all optional arguments in place'
      call printErrMsg(option)
  endif
          


  if (getDerivs) then
    D_deno = auxvar%D_den_kg(oid,:)
    D_deng = auxvar%D_den_kg(gid,:)
  endif

!--Get oil and gas saturations and form fractions------------------------------

  so=auxvar%sat(oid)
  sg=auxvar%sat(gid)

  den=so+sg
  if( den>0.0 ) then
    deninv=1.0/den
    fo=so*deninv
    fg=sg*deninv
    if (getDerivs) then
      D_worker = auxvar%D_sat(oid,:)+auxvar%D_sat(gid,:)
      D_fo = DivRule(so,auxvar%D_sat(oid,:),den,D_worker,ndof)
      D_fg = DivRule(sg,auxvar%D_sat(gid,:),den,D_worker,ndof)
    endif
  else
    fo=0.5
    fg=0.5
    if (getDerivs) then
      D_fo = 0.d0
      D_fg = 0.d0
    endif
  endif

!--Split the hydrocarbon relative permeability using fractions (TL,eqn. 2a,2b)-

  krotl=fo*krh
  krgtl=fg*krh
  if (getDerivs) then
   D_krotl = ProdRule(fo,D_fo,krh,D_krh,ndof)
   D_krgtl = ProdRule(fg,D_fg,krh,D_krh,ndof)
  endif

!--Form the omega-weighted oil and gas viscosities-----------------------------
  
  if (getDerivs) then
    call vToddLongstaffViscosity(fo,fg,so,sg,visco,viscg,viscotl,viscgtl,PETSC_TRUE,ndof &
                                 ,option                                                 &
                                 ,D_fo,D_fg,D_visco,D_viscg,D_viscotl,D_viscgtl)
  else
    call vToddLongstaffViscosity(fo,fg,so,sg,visco,viscg,viscotl,viscgtl,PETSC_FALSE,ndof&
                                ,option)
  endif

!--Form the omega weighted oil and gas densities (for Darcy flow only)---------

  if (getDerivs) then
    call vToddLongstaffDensity( fo,fg,visco,viscg,viscotl,viscgtl &
                               ,deno,deng,denotl,dengtl,PETSC_TRUE,ndof &
                               ,option                                  &
                               ,D_deno,D_deng                           &
                               ,D_visco,D_viscg,D_fo,D_fg               &
                               ,D_viscotl,D_viscgtl, D_denotl,D_dengtl)
  else
    call vToddLongstaffDensity( fo,fg,visco,viscg,viscotl,viscgtl &
                               ,deno,deng,denotl,dengtl,PETSC_FALSE,NDOF,option)
  endif

end subroutine vToddLongstaff

!==============================================================================

subroutine vToddLongstaffViscosity(fo,fg,so,sg,visco,viscg,viscotl,viscgtl &
                                  ,getDerivs,ndof,option,D_fo,D_fg,D_visco,D_viscg &
                                  ,D_viscotl,D_viscgtl)

!------------------------------------------------------------------------------
! Routine to calculate the Todd-Longstaff mixed viscosities
! Based on the reference TL (see header of calling routine for details)
!------------------------------------------------------------------------------
! input : fo     ,fg     : Oil and gas fractions
! input : visco  ,viscg  : Unmixed phase viscosity values on input
! output: viscotl,viscgtl: Mixed   phase viscosity values on output
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Jul 2017
!------------------------------------------------------------------------------
  use Derivatives_utilities_module
  use Option_module
  implicit none

  PetscReal,intent(in ):: fo,fg,so,sg,visco,viscg

  PetscBool,intent(in)    :: getDerivs
  PetscInt,intent(in)     :: ndof
  PetscReal,dimension(1:ndof),optional,intent(in ):: D_fo,D_fg,D_visco,D_viscg
  type(option_type) :: option

  PetscReal,intent(out):: viscotl,viscgtl
  PetscReal,dimension(1:ndof),optional,intent(out):: D_viscotl,D_viscgtl
  PetscReal            :: tlomegac,sn,sni,viscqpo,viscqpg,wviscqp &
                         ,denom,denominv,viscm
  PetscReal            :: viscoimw,viscgimw,viscmw
  PetscReal            :: workerReal,numerator
  PetscReal,dimension(1:ndof)  :: D_viscqpo,D_viscqpg,D_wviscqp,D_denom,D_viscm
  PetscReal,dimension(1:ndof)  :: D_viscoimw,D_viscgimw,D_viscmw,D_numerator

  if (getDerivs .AND. &
      .NOT.(present(D_fo).AND.present(D_fg).AND.present(D_visco).AND.present(D_viscg).AND.&
            present(D_viscotl).AND.present(D_viscgtl)                   )) then
      option%io_buffer = 'vToddLongstaffViscosity(): attempting to get analytical &
                          derivatives but not all optional arguments in place'
      call printErrMsg(option)
  endif

!--Set up complement of the Todd Longstaff omega-------------------------------

  tlomegac=1.0-val_tl_omega

!--Form quarter-powers of basic viscosities------------------------------------

  viscqpo=visco**0.25
  viscqpg=viscg**0.25

  if (getDerivs) then
    workerReal = 0.25
    D_viscqpo=PowerRule(visco,D_visco,workerReal,ndof)
    D_viscqpg=PowerRule(viscg,D_viscg,workerReal,ndof)
  endif

!--Form weighted combination of the 1/4 powers & its 4th power (denom of TL 4a)
!  Note that:
!  fg multiplies the oil viscosity 1/4-power term
!  fo multiplies the gas viscosity 1/4-power term

  wviscqp=fg*viscqpo+fo*viscqpg
  denom=wviscqp**4.0
  if (getDerivs) then
    D_wviscqp=ProdRule(fg,D_fg,viscqpo,D_viscqpo,ndof) &
             +ProdRule(fo,D_fo,viscqpg,D_viscqpg,ndof)
    workerReal = 4.0
    D_denom=PowerRule(wviscqp,D_wviscqp,workerReal,ndof)

  endif

!--Obtain a safe denominator inverse-------------------------------------------

  if( denom>0.0 ) then
    denominv=1.0/denom
  else
    denominv=0.0
  endif

!--Form mixed viscosity--(TL equation 4a)--------------------------------------

  viscm=visco*viscg*denominv
  if (getDerivs) then
    numerator = visco*viscg
    D_numerator = ProdRule(visco,D_visco,viscg,D_viscg,ndof)
    D_viscm = DivRule(numerator,D_numerator,denom,D_denom,ndof)
    !D_viscm = ProdRule3(visco,D_visco,viscg,D_viscg,denom,D_denom,ndof)
  endif

!--Form omega-weighted replacement viscosities---------------------------------

! 1-omega and omega power contributions

  viscoimw=visco**tlomegac
  viscgimw=viscg**tlomegac
  viscmw  =viscm**val_tl_omega

  if (getDerivs) then
    workerReal = tlomegac
    D_viscoimw= PowerRule(visco,D_visco,workerReal,ndof)
    D_viscgimw= PowerRule(viscg,D_viscg,workerReal,ndof)
    workerReal = val_tl_omega
    D_viscmw= PowerRule(viscm,D_viscm,workerReal,ndof)
  endif

! Combine to get final value (TL 3a and 3b)

  viscotl=viscoimw*viscmw
  viscgtl=viscgimw*viscmw

  if (getDerivs) then
    D_viscotl=ProdRule(viscoimw,D_viscoimw,viscmw,D_viscmw,ndof)
    D_viscgtl=ProdRule(viscgimw,D_viscgimw,viscmw,D_viscmw,ndof)
  endif

end subroutine vToddLongstaffViscosity

!------------------------------------------------------------------------------

subroutine vToddLongstaffDensity( fo,fg,visco,viscg,viscotl,viscgtl &
                                 ,deno,deng,denotl,dengtl,getDerivs,ndof &
                                 ,option                                 &
                                 ,D_deno,D_deng,D_visco,D_viscg          &
                                 ,D_fo,D_fg                              &
                                 ,D_viscotl,D_viscgtl                    &
                                 ,D_denotl,D_dengtl)

!------------------------------------------------------------------------------
! Routine to calculate the Todd-Longstaff mixed densities
! Based on the reference TL (see header of calling routine for details)
!------------------------------------------------------------------------------
! input : fo     ,fg     : Oil and gas fractions
! input : visco  ,viscg  : unmixed phase viscosity values
! input : viscotl,viscgtl: mixed   phase viscosity values
! input : deno   ,deng   : unmixed phase density   values
! output: denotl ,dengtl : mixed   phase density   values
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Jul 2017
!------------------------------------------------------------------------------

  use Derivatives_utilities_module
  use Option_module
  implicit none

  PetscReal,intent(in) :: fo,fg,visco,viscg,viscotl,viscgtl,deno,deng

  PetscBool,intent(in)    :: getDerivs
  PetscInt,intent(in)     :: ndof
  PetscReal,dimension(1:ndof),optional,intent(in ):: D_deno,D_deng,D_visco,D_viscg
  PetscReal,dimension(1:ndof),optional,intent(in ):: D_viscotl,D_viscgtl,D_fo,D_fg

  type(option_type) :: option

  PetscReal,intent(out):: denotl,dengtl

  PetscReal,dimension(1:ndof),optional,intent(out):: D_denotl,D_dengtl

  PetscReal            :: tlomegac,viscginv,viscotlinv,viscgtlinv &
                         ,m,mooe,moge,mqp,mooeqp,mogeqp,mqpm,mqpm1,mqpm1inv &
                         ,fo_oe,fo_ge,fg_oe,fg_ge,denm,workerReal

  PetscReal,dimension(1:ndof) :: D_m,D_mqp,D_mooeqp,D_mogeqp,D_mooe,D_moge
  PetscReal,dimension(1:ndof) :: D_fo_oe,D_fo_ge,D_denm,P1,P2

  if (getDerivs .AND. &
      .NOT.(present(D_deno).AND.present(D_deng).AND.present(D_visco).AND.present(D_viscg).AND.&
            present(D_viscotl).AND.present(D_viscgtl).AND.present(D_fo).AND.present(D_fg).AND.&
            present(D_denotl).AND.present(D_dengtl)                                            )) then
      option%io_buffer = 'vToddLongstaffDensity(): attempting to get analytical &
                          derivatives but not all optional arguments in place'
      call printErrMsg(option)
  endif

!--Set up complement of omega--------------------------------------------------

  tlomegac=1.0-val_tl_omega

! Build quarter power of mobility ratio of:
!
!     m   =visco/viscg  , original oil and original     gas viscosity
!     mooe=visco/viscotl, original oil and effective TL oil viscosity
!     moge=visco/viscgtl, original oil and effective TL gas viscosity

  viscginv=0.0
  if( viscg>0.0 ) then
    viscginv=1.0/viscg
  endif

  viscotlinv=0.0
  if( viscotl>0.0 ) then
    viscotlinv=1.0/viscotl
  endif

  viscgtlinv=0.0
  if( viscgtl>0.0 ) then
    viscgtlinv=1.0/viscgtl
  endif

  m=visco*viscginv
  if (getDerivs) then
    D_m = DivRule(visco,D_visco,viscg,D_viscg,ndof)
  endif
  if( abs(m-1.0)>1.0E-6 .and. viscotl>0.0 .and. viscgtl>0.0 ) then

! Case of non-unit mobility ratio

! Form (visco(unmixed))/viscptl), p=o,g, as used in eqn. 8b and 8a in TL

    mooe=visco*viscotlinv
    moge=visco*viscgtlinv
    if (getDerivs) then 
      D_mooe=DivRule(visco,D_visco,viscotl,D_viscotl,ndof)
      D_moge=DivRule(visco,D_visco,viscgtl,D_viscgtl,ndof)
    endif

    mqp   =m   **0.25
    mooeqp=mooe**0.25
    mogeqp=moge**0.25

    if (getDerivs) then 
      workerReal = 0.25
      D_mqp   =PowerRule(m,D_m,workerReal,ndof)
      D_mooeqp=PowerRule(mooe,D_mooe,workerReal,ndof)
      D_mogeqp=PowerRule(moge,D_moge,workerReal,ndof)   
    endif

    mqpm1=mqp-1.0

    mqpm1inv=0.0
    if( abs(mqpm1)>0.0 ) then
      mqpm1inv=1.0/mqpm1
    endif

!  Form effective fractional saturations, eqn. 8b and 8a in TL

    fo_oe=(mqp-mooeqp)*mqpm1inv
    fo_ge=(mqp-mogeqp)*mqpm1inv
    if (getDerivs) then 
      D_fo_oe = DivRule(mqp,D_mqp,mqpm1,D_mqp,ndof) &
              - DivRule(mooeqp,D_mooeqp,mqpm1,D_mqp,ndof)
      D_fo_ge = DivRule(mqp,D_mqp,mqpm1,D_mqp,ndof) &
              - DivRule(mogeqp,D_mogeqp,mqpm1,D_mqp,ndof)
    endif

!  Complements of effective oil saturations, as used in eqn. 9a and 9b

    fg_oe=1.0-fo_oe
    fg_ge=1.0-fo_ge

!  Set up mixed densities, eqn. 9b and 9a

    denotl=deno*fo_oe+deng*fg_oe
    dengtl=deno*fo_ge+deng*fg_ge

    if (getDerivs) then 
      D_denotl = D_deno*  fo_oe &
               +   deno*D_fo_oe &
               + D_deng*  fg_oe & 
               -   deng*D_fo_oe ! D_fg_oe = -D_fo_oe

      D_dengtl = D_deno*  fo_ge &
               +   deno*D_fo_ge &
               + D_deng*  fg_ge & 
               -   deng*D_fo_ge ! D_fg_ge = -D_fo_ge
    endif

  else

! Case of unit mobility ratio, eqn. 10a and 10b

    denm=fo*deno+fg*deng
    denotl=tlomegac*deno+val_tl_omega*denm
    dengtl=tlomegac*deng+val_tl_omega*denm

    if (getDerivs) then 
      D_denm(:)=ProdRule(fo,D_fo,deno,D_deno,ndof) &
               +ProdRule(fg,D_fg,deng,D_deng,ndof)
      D_denotl=tlomegac*D_deno+val_tl_omega*D_denm
      D_dengtl=tlomegac*D_deng+val_tl_omega*D_denm
    endif

  endif

end subroutine vToddLongstaffDensity

subroutine TOWGImsTLAuxVarPerturb(auxvar,global_auxvar, &
                                  material_auxvar, &
                                  characteristic_curves,natural_id, &
                                  option)
  !
  ! Calculates auxiliary variables for perturbed system
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 12/27/16
  !

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class
  use Derivative_tests_module

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(auxvar_towg_type) :: auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), &
               pert(option%nflowdof), x_pert_save(option%nflowdof)

  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: tempreal
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-10
  PetscInt :: idof

  x(TOWG_OIL_PRESSURE_DOF) = auxvar(ZERO_INTEGER)%pres(option%oil_phase)
  x(TOWG_OIL_SATURATION_DOF) = auxvar(ZERO_INTEGER)%sat(option%oil_phase)
  x(TOWG_GAS_SATURATION_3PH_DOF) = auxvar(ZERO_INTEGER)%sat(option%gas_phase)
  x(TOWG_3CMPS_ENERGY_DOF) = auxvar(ZERO_INTEGER)%temp

  pert(TOWG_OIL_PRESSURE_DOF) = &
    perturbation_tolerance*x(TOWG_OIL_PRESSURE_DOF)+min_perturbation
  pert(TOWG_3CMPS_ENERGY_DOF) = &
    perturbation_tolerance*x(TOWG_3CMPS_ENERGY_DOF)+min_perturbation
  if (x(TOWG_OIL_SATURATION_DOF) > 0.5d0) then
    pert(TOWG_OIL_SATURATION_DOF) = -1.d0 * perturbation_tolerance
  else
    pert(TOWG_OIL_SATURATION_DOF) = perturbation_tolerance
  endif
  if (x(TOWG_GAS_SATURATION_3PH_DOF) > 0.5d0) then
    pert(TOWG_GAS_SATURATION_3PH_DOF) = -1.d0 * perturbation_tolerance
  else
    pert(TOWG_GAS_SATURATION_3PH_DOF) = perturbation_tolerance
  endif

  ! TOWG_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = TOWG_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    auxvar(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    x_pert_save = x_pert
    call TOWGAuxVarCompute(x_pert,auxvar(idof),global_auxvar, &
                           material_auxvar, &
                           characteristic_curves,natural_id,option)
  enddo

  auxvar(TOWG_OIL_PRESSURE_DOF)%pert = &
     auxvar(TOWG_OIL_PRESSURE_DOF)%pert / towg_pressure_scale


#if 0
  if (.NOT. option%flow%numerical_derivatives .AND. option%flow%num_as_alyt_derivs) then
    call NumCompare_tl3p(option%nphase,option%nflowdof,auxvar,option,&
                            TOWG_OIL_PRESSURE_DOF,TOWG_OIL_SATURATION_DOF,&
                            TOWG_GAS_SATURATION_3PH_DOF,towg_energy_dof)
  endif
#endif
  if (option%flow%numerical_derivatives_compare) then 
    call NumCompare_tl3p(option%nphase,option%nflowdof,auxvar,option,&
                            TOWG_OIL_PRESSURE_DOF,TOWG_OIL_SATURATION_DOF,&
                            TOWG_GAS_SATURATION_3PH_DOF,towg_energy_dof)
  endif

end subroutine TOWGImsTLAuxVarPerturb

! ************************************************************************** !

subroutine TOWGBlackOilAuxVarPerturb(auxvar,global_auxvar, &
                                     material_auxvar, &
                                     characteristic_curves,natural_id, &
                                     option)
!------------------------------------------------------------------------------
! Set up auxillary variables for perturbed system
! Used in TOWG_BLACK_OIL
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Oct 2017
!------------------------------------------------------------------------------

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class
  use Derivative_tests_module

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(auxvar_towg_type) :: auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), &
               pert(option%nflowdof), x_pert_save(option%nflowdof)

  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: tempreal
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-10
  PetscInt :: idof,saturated_state
  PetscBool::isSaturated

!--Look at the state and set saturated oil flag--------------------------------

  saturated_state=global_auxvar%istate
  if( saturated_state==TOWG_THREE_PHASE_STATE ) then
    isSaturated=.true.
  else
    isSaturated=.false.
  endif

!--Set up perturbed solution---------------------------------------------------

  x(TOWG_OIL_PRESSURE_DOF) = auxvar(ZERO_INTEGER)%pres(option%oil_phase)
  x(TOWG_OIL_SATURATION_DOF) = auxvar(ZERO_INTEGER)%sat(option%oil_phase)
  if( isSaturated ) then
    x(TOWG_GAS_SATURATION_3PH_DOF) = auxvar(ZERO_INTEGER)%sat(option%gas_phase)
  else
    x(TOWG_BUBBLE_POINT_3PH_DOF  ) = auxvar(ZERO_INTEGER)%bo%bubble_point
  endif
  x(TOWG_3CMPS_ENERGY_DOF) = auxvar(ZERO_INTEGER)%temp

  pert(TOWG_OIL_PRESSURE_DOF) = &
    perturbation_tolerance*x(TOWG_OIL_PRESSURE_DOF)+min_perturbation
  pert(TOWG_3CMPS_ENERGY_DOF) = &
    perturbation_tolerance*x(TOWG_3CMPS_ENERGY_DOF)+min_perturbation
  if (x(TOWG_OIL_SATURATION_DOF) > 0.5d0) then
    pert(TOWG_OIL_SATURATION_DOF) = -1.d0 * perturbation_tolerance
  else
    pert(TOWG_OIL_SATURATION_DOF) = perturbation_tolerance
  endif

  if( isSaturated ) then
    if (x(TOWG_GAS_SATURATION_3PH_DOF) > 0.5d0) then
      pert(TOWG_GAS_SATURATION_3PH_DOF) = -1.d0 * perturbation_tolerance
    else
      pert(TOWG_GAS_SATURATION_3PH_DOF) = perturbation_tolerance
    endif
  else
    pert(TOWG_BUBBLE_POINT_3PH_DOF) = &
    perturbation_tolerance*x(TOWG_BUBBLE_POINT_3PH_DOF)+min_perturbation
  endif

  ! TOWG_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = TOWG_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    auxvar(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    x_pert_save = x_pert
    call TOWGAuxVarCompute(x_pert,auxvar(idof),global_auxvar, &
                           material_auxvar, &
                           characteristic_curves,natural_id,option)
  enddo

  auxvar(TOWG_OIL_PRESSURE_DOF)%pert = &
     auxvar(TOWG_OIL_PRESSURE_DOF)%pert / towg_pressure_scale

  if (option%flow%numerical_derivatives_compare) then 
    call NumCompare_towg_bo(option%nphase,option%nflowdof,auxvar,option,&
                            TOWG_OIL_PRESSURE_DOF,TOWG_OIL_SATURATION_DOF,&
                            TOWG_GAS_SATURATION_3PH_DOF,towg_energy_dof,&
                            isSaturated)
  endif

end subroutine TOWGBlackOilAuxVarPerturb

! ************************************************************************** !

subroutine TL4PAuxVarPerturb(auxvar,global_auxvar, &
                             material_auxvar, &
                             characteristic_curves,natural_id, &
                             option)
!------------------------------------------------------------------------------
! Set up auxillary variables for perturbed system
! Used in TOWG_SOLVENT_TL
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class
  use Derivative_tests_module


  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(auxvar_towg_type) :: auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), &
               pert(option%nflowdof), x_pert_save(option%nflowdof)

  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: tempreal
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-10
  PetscInt :: idof,saturated_state
  PetscBool::isSaturated

!--Look at the state and set saturated oil flag--------------------------------

  saturated_state=global_auxvar%istate
  if( saturated_state==TOWG_THREE_PHASE_STATE ) then
    isSaturated=.true.
  else
    isSaturated=.false.
  endif

!--Set up perturbed solution---------------------------------------------------

! Oil pressure
  x(TOWG_OIL_PRESSURE_DOF) = auxvar(ZERO_INTEGER)%pres(option%oil_phase)

! Oil satn
  x(TOWG_OIL_SATURATION_DOF) = auxvar(ZERO_INTEGER)%sat(option%oil_phase)

! Gas saturation or bubble point
  if( isSaturated ) then
    x(TOWG_GAS_SATURATION_3PH_DOF) = auxvar(ZERO_INTEGER)%sat(option%gas_phase)
  else
    x(TOWG_BUBBLE_POINT_3PH_DOF  ) = auxvar(ZERO_INTEGER)%bo%bubble_point
  endif

! Solvent saturation
  x(TOWG_SOLV_SATURATION_DOF) =  auxvar(ZERO_INTEGER)%sat(option%solvent_phase)

! TL4P energy location
  x(TOWG_SOLV_TL_ENERGY_DOF) = auxvar(ZERO_INTEGER)%temp

!--Now the actual perturbation-------------------------------------------------

! Oil pressure
  pert(TOWG_OIL_PRESSURE_DOF) = &
    perturbation_tolerance*x(TOWG_OIL_PRESSURE_DOF)+min_perturbation

! Oil saturation
  if (x(TOWG_OIL_SATURATION_DOF) > 0.5d0) then
    pert(TOWG_OIL_SATURATION_DOF) = -1.d0 * perturbation_tolerance
  else
    pert(TOWG_OIL_SATURATION_DOF) = perturbation_tolerance
  endif

! Gas saturation or bubble point
  if( isSaturated ) then
    if (x(TOWG_GAS_SATURATION_3PH_DOF) > 0.5d0) then
      pert(TOWG_GAS_SATURATION_3PH_DOF) = -1.d0 * perturbation_tolerance
    else
      pert(TOWG_GAS_SATURATION_3PH_DOF) = perturbation_tolerance
    endif
  else
    pert(TOWG_BUBBLE_POINT_3PH_DOF) = &
    perturbation_tolerance*x(TOWG_BUBBLE_POINT_3PH_DOF)+min_perturbation
  endif

! Solvent saturation
  if (x(TOWG_SOLV_SATURATION_DOF) > 0.5d0) then
    pert(TOWG_SOLV_SATURATION_DOF) = -1.d0 * perturbation_tolerance
  else
    pert(TOWG_SOLV_SATURATION_DOF) = perturbation_tolerance
  endif

! TL4P energy location
  pert(TOWG_SOLV_TL_ENERGY_DOF) = &
    perturbation_tolerance*x(TOWG_SOLV_TL_ENERGY_DOF)+min_perturbation

!--TOWG_UPDATE_FOR_DERIVATIVE indicates call from perturbation-----------------

  option%iflag = TOWG_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    auxvar(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    x_pert_save = x_pert
    call TOWGAuxVarCompute(x_pert,auxvar(idof),global_auxvar, &
                           material_auxvar, &
                           characteristic_curves,natural_id,option)
  enddo

  auxvar(TOWG_OIL_PRESSURE_DOF)%pert = &
     auxvar(TOWG_OIL_PRESSURE_DOF)%pert / towg_pressure_scale


#if 0
!!! no longer functional, was used for development testing
  if (.NOT. option%flow%numerical_derivatives .AND. option%flow%num_as_alyt_derivs) then
    call Num_as_alyt_tl4p(option%nphase,option%nflowdof,auxvar,option,&
                            TOWG_OIL_PRESSURE_DOF,TOWG_OIL_SATURATION_DOF,&
                            TOWG_GAS_SATURATION_3PH_DOF,towg_energy_dof, &
                            isSaturated)
  endif
#endif

  if (.NOT. option%flow%numerical_derivatives .AND. option%flow%numerical_derivatives_compare) then 
    call NumCompare_tl4p(option%nphase,option%nflowdof,auxvar,option,&
                            TOWG_OIL_PRESSURE_DOF,TOWG_OIL_SATURATION_DOF,&
                            TOWG_GAS_SATURATION_3PH_DOF,towg_energy_dof, &
                            isSaturated)
  endif


end subroutine TL4PAuxVarPerturb

! ************************************************************************** !

subroutine TOWGImsAuxVarComputeSetup()

  implicit none

  TOWGAuxVarCompute => TOWGImsAuxVarCompute
  TOWGAuxVarPerturb => TOWGImsTLAuxVarPerturb

end subroutine TOWGImsAuxVarComputeSetup

! ************************************************************************** !

subroutine TOWGBlackOilAuxVarComputeSetup()

!------------------------------------------------------------------------------
! Used in TOWG_BLACK_OIL, set function pointers to required routines to
! do auxillary variable main setup and perturbation calculation
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Sep 2017
!------------------------------------------------------------------------------

  implicit none

  TOWGAuxVarCompute => TOWGBlackOilAuxVarCompute
  TOWGAuxVarPerturb => TOWGBlackOilAuxVarPerturb

end subroutine TOWGBlackOilAuxVarComputeSetup

! ************************************************************************** !

subroutine TL4PAuxVarComputeSetup()

!------------------------------------------------------------------------------
! Used in TOWG_SOLVENT_TL, set function pointers to required routines to
! do auxillary variable main setup and perturbation calculation
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  implicit none

  TOWGAuxVarCompute => TL4PAuxVarCompute
  TOWGAuxVarPerturb => TL4PAuxVarPerturb

end subroutine TL4PAuxVarComputeSetup

! ************************************************************************** !

subroutine TOWGTLAuxVarComputeSetup()

  implicit none

  TOWGAuxVarCompute => TOWGTLAuxVarCompute
  TOWGAuxVarPerturb => TOWGImsTLAuxVarPerturb

end subroutine TOWGTLAuxVarComputeSetup

! ************************************************************************** !

subroutine TOWGAuxDestroy(aux)
  !
  ! Deallocates a towg auxiliary object
  !
  ! Author: Paolo Orsini
  ! Date: 11/07/16
  !
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_towg_aux_type), pointer :: aux
  PetscInt :: iaux, idof

  if (.not.associated(aux)) return

  if (associated(aux%auxvars) ) then
    call TOWGAuxVarStrip(aux%auxvars)
    deallocate(aux%auxvars)
  end if
  nullify(aux%auxvars)

  if (associated(aux%auxvars_bc) ) then
    call TOWGAuxVarStrip(aux%auxvars_bc)
    deallocate(aux%auxvars_bc)
  end if
  nullify(aux%auxvars_bc)

  if ( associated(aux%auxvars_ss) ) then
    call TOWGAuxVarStrip(aux%auxvars_ss)
    deallocate(aux%auxvars_ss)
  end if
  nullify(aux%auxvars_ss)

  call PMBaseAuxStrip(aux)

  if (associated(aux%parameter)) then
    deallocate(aux%parameter)
    !to add paramter strip when introducing diff in oil and gas phases
  end if
  nullify(aux%parameter)

  deallocate(aux)
  nullify(aux)

end subroutine TOWGAuxDestroy

! ************************************************************************** !

! ************************************************************************** !

subroutine  TOWGAuxVarArray1Strip(auxvars)
  !
  ! Deallocates a mode auxiliary object
  ! this could be generalised for different modes
  ! using class(*) (unlimited polymorphic)
  !
  ! Author: Paolo Orsini
  ! Date: 11/07/16
  !

  use AuxVars_TOWG_module

  implicit none

  !can't use class due to gfortran (4.8.4) bug (values passed are not correct)
  !bug fixed in gfortran 6.2.0
  !here we can pass by pointer, we could destroy the array within the routine
  !but we don't to be consistent with TOilImsAuxVarArray2Strip
  !class(auxvar_towg_type), pointer :: auxvars(:)
  type(auxvar_towg_type) :: auxvars(:)

  PetscInt :: iaux

  !print *, "den oil bc/ss pass = ", auxvars(1)%den(2)

  do iaux = 1, size(auxvars)
    call auxvars(iaux)%Strip
  enddo

end subroutine TOWGAuxVarArray1Strip

! ************************************************************************** !

subroutine TOWGAuxVarArray2Strip(auxvars)
  !
  ! Deallocates a mode auxiliary object
  ! this could be generalised for different modes
  ! using class(*) (unlimited polymorphic)
  !
  ! Author: Paolo Orsini
  ! Date: 11/07/16
  !

  use AuxVars_TOWG_module

  implicit none

  !can't use class due to gfortran (4.8.4) bug (values passed are not correct)
  !bug fixed in gfortran 6.2.0
  !cannot use type(...) with pointer attribute, therefore we deallocate and
  !nullify pointer outide this routine
  !because the compiler does not allow to specify lower 0-bound in auxvar
  !type(auxvar_towg_type), pointer :: auxvars(:,:)
  !class(auxvar_towg_type) :: auxvars(0:,:)
  type(auxvar_towg_type) :: auxvars(0:,:)

  PetscInt :: iaux, idof

  do iaux = 1, size(auxvars,2)
    do idof = 1, size(auxvars,1)
      call auxvars(idof-1,iaux)%Strip()
    enddo
  enddo

end subroutine TOWGAuxVarArray2Strip

subroutine TL4PMiscibilityFraction(sg,ss,fm,dfmdsg,dfmdss)

!------------------------------------------------------------------------------
! Set up the TL miscibility fraction
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  implicit none

  PetscReal,intent(in ) :: sg,ss
  PetscReal,intent(out) :: fm,dfmdsg,dfmdss
  PetscReal             :: svi,fs,sv,dfmdfs,svsq,svsqi,sgi,ssi

  fm    =0.0

  dfmdsg=0.0
  dfmdss=0.0

  sv=sg+ss
  if( sv>0.0 ) then
    svi=1.0/sv
    fs=ss*svi

#if 0
    ! to be considered - the rounding errors this avoids come up often. Note this has indeed been observed to change
    ! the convergence behaviour of numerical and anlytical runs slightly
    if (sg ==0.d0) fs = 1.d0
#endif

    dfmdfs=0.d0
    call TL4PMiscibilityFractionFromSaturationFraction(fs,fm,dfmdfs)

    ! derivatives: note that the derivs of fs can be very problematic
    ! for near 0 saturations.
    dfmdsg = 0.0
    dfmdss = 0.0
    if (sg == 0.d0) then ! attempt to catch case of sg = 0, ss very small
      ! explicitly, subbing in sg = 0, we have:
      ! d fs / d sg = - ss / (ss+sg)^2 = -1/ss
      ! d fs / d ss =   sg / (ss+sg)^2 =  0
      ! so we do the following:
      ssi = 0.0
      if (ss > 0.0) ssi = 1.0/ss
      dfmdsg = -dfmdfs*ssi

      ! leave as default 0:
      !dfmdss =  0.0
    elseif (ss == 0.0) then
      ! could happen that ss = 0 and sg very small can give similar problems, so:
      ! d fs / d sg = - ss / (ss+sg)^2 = 0
      ! d fs / d ss =   sg / (ss+sg)^2 = 1/sg
      ! then:
      !dfmdsg = 0.0
      sgi = 0.d0
      if (sg > 0.d0) then 
        sgi = 1.d0/sg
        dfmdss =  dfmdfs*sgi
      endif
    else ! standard trap but can still fail
      svsq = sv*sv
      svsqi = 1.0/svsq
      if (svsq>0.0) then
        dfmdsg = -dfmdfs*ss*svsqi
        dfmdss =  dfmdfs*sg*svsqi
      endif
    endif

#if 0
!!! I'll leave these commented out but we may well want a nice robust error check for this at some point just in case
    if (isnan(dfmdsg)) then
      print *, "nan in fm computation derivs, dfmdsg"
    endif
    if (isnan(dfmdss)) then
      print *, "nan in fm computation derivs, dfmdss"
    endif
#endif

  endif

end subroutine TL4PMiscibilityFraction

! ************************************************************************** !

subroutine TL4PRelativePermeabilities( so  ,sg  ,sw   ,ss   ,sv  ,sh, fm , &
                                       swcr,sgcr,sowcr,sogcr,swco,         &
                                       kroi,krgi,krsi,                     &
                                       krom,krgm,krsm,                     &
                                       kro ,krg ,krs,                      &
                                       getDerivs,ndof,                     &
                                       D_fm ,                              & ! fm derivs (in)
                                       D_kroi,D_krgi,D_krsi,               & ! imiscible k derivs (in)
                                       D_krom,D_krgm,D_krsm,               & ! miskible k derivs  (in)
                                       D_kro ,D_krg ,D_krs   )               ! true k derivs     (out)

!------------------------------------------------------------------------------
! Set up the interpolated TL4P relative permeabilities
! Used in TOWG_BLACK_OIL and TOWG_SOLVENT_TL
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  use Derivatives_utilities_module
  implicit none

  PetscReal,intent(in )                    :: so  ,sg  ,sw   ,ss   ,sv  ,sh,fm 
  PetscReal,intent(in )                    :: swcr,sgcr,sowcr,sogcr,swco
  PetscReal,intent(in )                    :: kroi,krgi,krsi
  PetscReal,intent(in )                    :: krom,krgm,krsm
  PetscReal,dimension(1:ndof),intent(in)   :: D_fm
  PetscReal,dimension(1:ndof),intent(in)   :: D_kroi,D_krgi,D_krsi
  PetscReal,dimension(1:ndof),intent(in)   :: D_krom,D_krgm,D_krsm

  PetscReal,intent(out)                    :: kro ,krg ,krs
  PetscReal,dimension(1:ndof),intent(out)  :: D_kro,D_krg,D_krs

  PetscReal svi,fs,fi

  PetscBool :: getDerivs
  PetscInt :: ndof



  kro=0.0
  krg=0.0
  krs=0.0

!--Form new mixed results and derivatives--------------------------------------

  fi=1.0-fm

  kro=fm*krom+fi*kroi
  krg=fm*krgm+fi*krgi
  krs=fm*krsm+fi*krsi


  if (getDerivs) then

   D_kro = ProdRule(fm,D_fm,krom,D_krom,ndof) &
         + ProdRule(fi,-D_fm,kroi,D_kroi,ndof)   

   D_krg = ProdRule(fm,D_fm,krgm,D_krgm,ndof) &
         + ProdRule(fi,-D_fm,krgi,D_krgi,ndof)   

   D_krs = ProdRule(fm,D_fm,krsm,D_krsm,ndof) &
         + ProdRule(fi,-D_fm,krsi,D_krsi,ndof)   
  endif

end subroutine TL4PRelativePermeabilities

! ************************************************************************** !

subroutine TL4PViscosity( so   ,sg   ,ss                               , &
                           visco,viscg,viscs                           , &
                           viscotl,viscgtl,viscstl                     , &
                           getDerivs,ndof                              , &
                           D_so,D_sg,D_ss                              , &
                           dof_osat,dof_gsat,dof_ssat                  , &
                           D_visco,D_viscg,D_viscs                     , & ! visc derivs in
                           D_viscotl,D_viscgtl,D_viscstl               , & ! tl visc derivs out
                           isSat                                       , &
                           denos,D_denos,dengs,D_dengs,denogs,D_denogs , &
                           denos_pre,D_denos_pre)

!------------------------------------------------------------------------------
! Set up the interpolated TL4P viscosities
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  use Derivatives_utilities_module
  implicit none

  ! in :
  PetscReal,intent(in ) :: so     ,sg     ,ss
  PetscReal,intent(in ) :: visco  ,viscg  ,viscs
  PetscBool,intent(in ) :: getDerivs,isSat
  PetscInt ,intent(in ) :: ndof,dof_osat,dof_gsat,dof_ssat
  PetscReal,dimension(1:ndof),intent(in)   :: D_visco,D_viscg,D_viscs
  PetscReal,dimension(1:ndof),intent(in)   :: D_so,D_sg,D_ss


  ! out:
  PetscReal,dimension(1:ndof),intent(out)  :: D_viscotl,D_viscgtl,D_viscstl
  PetscReal,intent(out) :: viscotl,viscgtl,viscstl

  ! intermediates:
  PetscReal             :: viscom ,viscgm ,viscsm
  PetscReal             :: viscoimw,viscgimw,viscsimw
  PetscReal             :: sos,sosi
  PetscReal             :: sgs,sgsi
  PetscReal             :: sh ,shi
  PetscReal             :: voqp,vgqp,vsqp
  PetscReal             :: denos ,dengs ,denogs
  PetscReal             :: denosi,dengsi,denogsi
  PetscReal             :: tlomegac
  PetscReal             :: viscom_w,viscgm_w,viscsm_w
  PetscReal             :: visco_wc,viscg_wc,viscs_wc

  PetscReal,dimension(1:ndof)  :: D_voqp,D_vgqp,D_vsqp
  PetscReal,dimension(1:ndof)  :: D_denos,D_dengs,D_denogs
  PetscReal,dimension(1:ndof)  :: D_viscom,D_viscgm,D_viscsm
  PetscReal,dimension(1:ndof)  :: D_viscom_w,D_viscgm_w,D_viscsm_w
  PetscReal,dimension(1:ndof)  :: D_visco_wc,D_viscg_wc,D_viscs_wc
  PetscReal,dimension(1:ndof)  :: D_worker
  !PetscReal,dimension(1:ndof)  :: D_1,D_2,D_3,D_denos_alt
  PetscReal,dimension(1:ndof)  :: D_sos,D_sgs,D_sh
  PetscReal,dimension(1:ndof)  :: D_denos_pre
  PetscReal :: denos_pre
  PetscReal :: workerReal
  PetscReal :: sosi_sq,worker,sgsi_sq,shi_sq

  PetscInt :: i

!--Set up complement of the Todd Longstaff omega-------------------------------

  tlomegac=1.0-val_tl_omega

!--Saturation sums-------------------------------------------------------------

  sos=so+ss
  sosi=0.0
  if( sos>0.0 ) sosi=1.0/sos

  sgs=sg+ss
  sgsi=0.0
  if( sgs>0.0 ) sgsi=1.0/sgs

  sh=so+sg+ss
  shi=0.0
  if( sh>0.0  ) shi=1.0/sh

  if (getDerivs) then
    D_sos = D_so + D_ss
    D_sgs = D_sg + D_ss
    D_sh  = D_so + D_sg + D_ss
  endif

!--Use quarter-power mixing rule-----------------------------------------------

  voqp=visco**0.25
  vgqp=viscg**0.25
  vsqp=viscs**0.25

  if (getDerivs) then
    workerReal = 0.25
    D_voqp = PowerRule(visco,D_visco,workerReal,ndof)
    D_vgqp = PowerRule(viscg,D_viscg,workerReal,ndof)
    D_vsqp = PowerRule(viscs,D_viscs,workerReal,ndof)
  endif
!--/Use quarter-power mixing rule----------------------------------------------

!------------ denos and derivatives--------------------------------------------

  if( sosi>0.0 ) then
    denos= so*sosi*vsqp &
          +ss*sosi*voqp
  else
    denos= 0.5*(vsqp+voqp)
  endif

   if (getDerivs) then
    if (sosi>0.d0) then

      D_denos = ProdDivRule(so,D_so,vsqp,D_vsqp,sos,D_sos,ndof) &
              + ProdDivRule(ss,D_ss,voqp,D_voqp,sos,D_sos,ndof) 
    else
      D_denos= 0.5*(D_vsqp+D_voqp)
    endif
  endif

!------------ /denos and derivatives--------------------------------------------


!------------ dengs and derivatives--------------------------------------------
  if( sgsi>0.0 ) then
    dengs= sg*sgsi*vsqp &
          +ss*sgsi*vgqp
  else
    dengs= 0.5*(vsqp+vgqp)
  endif
   if (getDerivs) then
    if (sgsi>0.d0) then

      D_dengs = ProdDivRule(sg,D_sg,vsqp,D_vsqp,sgs,D_sgs,ndof) &
              + ProdDivRule(ss,D_ss,vgqp,D_vgqp,sgs,D_sgs,ndof) 
    else
      D_dengs= 0.5*(D_vsqp+D_vgqp)
    endif
  endif
!------------/dengs and derivatives--------------------------------------------

!------------ denogs and derivatives--------------------------------------------
  if( shi>0.0 ) then
    denogs = so*shi*vgqp*vsqp &
            +sg*shi*vsqp*voqp &
            +ss*shi*voqp*vgqp
  else
    denogs =( vgqp*vsqp &
             +vsqp*voqp &
             +voqp*vgqp )/3.0
  endif
   if (getDerivs) then
    if (shi>0.d0) then

     D_denogs = ProdProdDivRule(so,D_so,vgqp,D_vgqp,vsqp,D_vsqp,sh,D_sh,ndof) &
              + ProdProdDivRule(sg,D_sg,vsqp,D_vsqp,voqp,D_voqp,sh,D_sh,ndof) &
              + ProdProdDivRule(ss,D_ss,voqp,D_voqp,vgqp,D_vgqp,sh,D_sh,ndof)
    else
      D_denogs = ProdRule(vgqp,D_vgqp,vsqp,D_vsqp,ndof) &
         + ProdRule(vsqp,D_vsqp,voqp,D_voqp,ndof) &
         + ProdRule(voqp,D_voqp,vgqp,D_vgqp,ndof)
      D_denogs = D_denogs/3.0
    endif
   endif
!------------ /denogs and derivatives--------------------------------------------

  denos_pre = denos
  D_denos_pre = D_denos

  if (getDerivs) then
    workerReal = 4.d0
    D_denos  = PowerRule(denos,D_denos,workerReal,ndof)
    D_dengs  = PowerRule(dengs,D_dengs,workerReal,ndof)
    D_denogs = PowerRule(denogs,D_denogs,workerReal,ndof)
  endif

  denos =denos **4.0
  dengs =dengs **4.0
  denogs=denogs**4.0


  denosi =0.0
  if( denos >0.0 ) denosi =1.0/denos

  dengsi =0.0
  if( dengs >0.0 ) dengsi =1.0/dengs

  denogsi=0.0
  if( denogs>0.0 ) denogsi=1.0/denogs

!--Construct mixed forms-------------------------------------------------------

  viscom=visco      *viscs*denosi
  viscgm=viscg      *viscs*dengsi
  viscsm=visco*viscg*viscs*denogsi

  if (getDerivs) then
    D_viscom = ProdDivRule(visco,D_visco,viscs,D_viscs,denos,D_denos,ndof)
    D_viscgm = ProdDivRule(viscg,D_viscg,viscs,D_viscs,dengs,D_dengs,ndof)
    D_viscsm = ProdProdDivRule(visco,D_visco,viscg,D_viscg,viscs,D_viscs,denogs,D_denogs,ndof)
  endif

!--Construct mixed forms to power omega----------------------------------------

  viscom_w=viscom**val_tl_omega
  viscgm_w=viscgm**val_tl_omega
  viscsm_w=viscsm**val_tl_omega

  if (getDerivs) then
    D_viscom_w = PowerRule(viscom,D_viscom,val_tl_omega,ndof)
    D_viscgm_w = PowerRule(viscgm,D_viscgm,val_tl_omega,ndof)
    D_viscsm_w = PowerRule(viscsm,D_viscsm,val_tl_omega,ndof)
  endif

!--Construct original forms to power omega-------------------------------------

  visco_wc=visco**tlomegac
  viscg_wc=viscg**tlomegac
  viscs_wc=viscs**tlomegac

  if (getDerivs) then
    D_visco_wc = PowerRule(visco,D_visco,tlomegac,ndof)
    D_viscg_wc = PowerRule(viscg,D_viscg,tlomegac,ndof)
    D_viscs_wc = PowerRule(viscs,D_viscs,tlomegac,ndof)
  endif

!--Combine---------------------------------------------------------------------

  viscotl=viscom_w*visco_wc
  viscgtl=viscgm_w*viscg_wc
  viscstl=viscsm_w*viscs_wc

  if (getDerivs) then
    D_viscotl = ProdRule(viscom_w,D_viscom_w,visco_wc,D_visco_wc,ndof)
    D_viscgtl = ProdRule(viscgm_w,D_viscgm_w,viscg_wc,D_viscg_wc,ndof)
    D_viscstl = ProdRule(viscsm_w,D_viscsm_w,viscs_wc,D_viscs_wc,ndof)

! we may wish to uncomment and have as an optional error check at some point; 
! a lot can go wrong when the system gets to marginal saturations
#if 0
    do i = 1,ndof
      if (isnan(D_viscotl(i))) then
        print *, "viscotl nan deriv at ", i
      endif
      if (isnan(D_viscgtl(i))) then
        print *, "viscgtl nan deriv at ", i
      endif
      if (isnan(D_viscstl(i))) then
        print *, "viscstl nan deriv at ", i
      endif
    enddo
#endif
  endif


end subroutine TL4PViscosity

! ************************************************************************** !

subroutine TL4PDensity( so    ,sg    ,ss    ,       &
                        viso  ,visg  ,viss  ,       &
                        visotl,visgtl,visstl,       &
                        deno  ,deng  ,dens  ,       &
                        denotl,dengtl,denstl,       &
                        getDerivs, ndof     ,       &
                        D_so,D_sg,D_ss      ,       &
                        dof_osat,dof_gsat,dof_ssat, &
                        D_deno,D_deng,D_dens,       &
                        D_viso,D_visg,D_viss,       &
                        D_visotl,D_visgtl,D_visstl, &
                        D_denotl,D_dengtl,D_denstl, &
                        isSat               ,       &
                        denog,D_denog)

!------------------------------------------------------------------------------
! Set up the interpolated TL4P gravity densities
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  use Derivatives_utilities_module
  implicit none

  PetscReal,intent(in ) :: so    ,sg    ,ss
  PetscReal,intent(in ) :: viso  ,visg  ,viss
  PetscReal,intent(in ) :: visotl,visgtl,visstl
  PetscReal,intent(in ) :: deno  ,deng  ,dens


  PetscBool,intent(in ) :: getDerivs,isSat
  PetscInt ,intent(in ) :: ndof,dof_osat,dof_gsat,dof_ssat
  PetscReal,dimension(1:ndof),intent(in)  :: D_so,D_sg,D_ss
  PetscReal,dimension(1:ndof),intent(in)  :: D_deno,D_deng,D_dens
  PetscReal,dimension(1:ndof),intent(in)  :: D_viso,D_visg,D_viss
  PetscReal,dimension(1:ndof),intent(in)  :: D_visotl,D_visgtl,D_visstl

  PetscReal,intent(out) :: denotl,dengtl,denstl
  PetscReal,dimension(1:ndof),intent(out) :: D_denotl,D_dengtl,D_denstl

  PetscReal             :: tlomega,tlomegac,sh,shi,fo,fg,fs,denm,sog,sogi,go,gg
  PetscReal             :: visoqp,visgqp,d,d4,d4i,visog,denog,denmo,denmg,denms
  PetscReal,dimension(1:ndof) :: D_denm,D_denog
  PetscReal,dimension(1:ndof) :: D_d,D_d4,D_visoqp,D_visgqp,D_visog
  PetscReal,dimension(1:ndof) :: D_sog,D_sh,D_fo,D_fg,D_fs,D_go,D_gg
  PetscReal,dimension(1:ndof) :: D_denmo,D_denmg,D_denms
  PetscReal,dimension(1:ndof) :: D_worker
  PetscReal :: workerReal
  PetscReal :: sgoi_sq
  PetscInt :: i

!--Initialise------------------------------------------------------------------

  visoqp=0.0
  visgqp=0.0

!--Set up complement of the Todd Longstaff omega and hydrocarbon saturation----

  tlomega =    val_tl_omega
  tlomegac=1.0-val_tl_omega

  sh=so+sg+ss

  if (getDerivs) then
    D_sh  = D_so + D_sg + D_ss
  endif

!--Form the hydrocarbon saturation fractions-----------------------------------

  if( sh>0.0 ) then

    shi=1.0/sh

    fo=so*shi
    fg=sg*shi
    fs=ss*shi

    if (getDerivs) then
      D_fo = DivRule(so,D_so,sh,D_sh,ndof)
      D_fg = DivRule(sg,D_sg,sh,D_sh,ndof)
      D_fs = DivRule(ss,D_ss,sh,D_sh,ndof)
    endif

  else

    fo=1.0/3.0
    fg=1.0/3.0
    fs=1.0/3.0

    if (getDerivs) then
      D_fo = 0.0
      D_fg = 0.0
      D_fs = 0.0
    endif

  endif

!--Form the default mixed density----------------------------------------------

  denm= fo*deno &
       +fg*deng &
       +fs*dens
  ! these aren't needed?
  denmo=denm
  denmg=denm
  denms=denm


  if (getDerivs) then
    D_denm = ProdRule(fo,D_fo,deno,D_deno,ndof) &
           + ProdRule(fg,D_fg,deng,D_deng,ndof) &
           + ProdRule(fs,D_fs,dens,D_dens,ndof) 

    ! don't need this if denmo etc aren't needed
    D_denmo = D_denm
    D_denmg = D_denm
    D_denms = D_denm
  endif

!--Omega-dependent mixing to final density-------------------------------------

  denotl=tlomegac*deno+tlomega*denm
  dengtl=tlomegac*deng+tlomega*denm
  denstl=tlomegac*dens+tlomega*denm

  if (getDerivs) then
    D_denotl=tlomegac*D_deno+tlomega*D_denm
    D_dengtl=tlomegac*D_deng+tlomega*D_denm
    D_denstl=tlomegac*D_dens+tlomega*D_denm
  endif

!--Now the main calculation: start with forming an oil-gas viscosity estimate--

  sog=so+sg
  if( sog>0.0 ) then
    sogi=1.0/sog
    go=so*sogi
    gg=sg*sogi
    if (getDerivs) then
      D_sog=D_so+D_sg
      D_go = DivRule(so,D_so,sog,D_sog,ndof)
      D_gg = DivRule(sg,D_sg,sog,D_sog,ndof)
    endif
  else
    go=0.5
    gg=0.5
    if (getDerivs) then
      D_sog=0.0
      D_go = 0.0 
      D_gg = 0.0
    endif
  endif

  denog=go*deno+gg*deng
  if (getDerivs) then
    if( sog>0.0 ) then
      D_denog = ProdRule(go,D_go,deno,D_deno,ndof) &
              + ProdRule(gg,D_gg,deng,D_deng,ndof)
    else
      ! better handling of this case may be possible
      D_denog = 0.0
    endif
  endif

!--Form the quarter-powers of the two viscosities (unless zero value is OK)----

  if( viso>0.0 ) visoqp=viso**0.25
  if( visg>0.0 ) visgqp=visg**0.25

  if (getDerivs) then
    workerReal = 0.25
    if (viso>0.0) then
      D_visoqp = PowerRule(viso,D_viso,workerReal,ndof)
    else
      D_visoqp = 0.d0
    endif
    if (visg>0.0) then
      D_visgqp = PowerRule(visg,D_visg,workerReal,ndof)
    else
      D_visgqp = 0.d0
    endif
  endif

!--Form 1/4 power combination (note indices to get inverses)---------------- --

  d=go*visgqp+gg*visoqp
  d4=d**4.0
  if( d4>0.0 ) d4i=1.0/d4

  if (getDerivs) then 
    D_d = ProdRule(go,D_go,visgqp,D_visgqp,ndof) &
        + ProdRule(gg,D_gg,visoqp,D_visoqp,ndof)
    workerReal = 4.0
    D_d4 = PowerRule(d,D_d,workerReal,ndof)

  endif

! Form value
  visog=viso*visg*d4i
  if (getDerivs) then
    D_visog = ProdDivRule(viso,D_viso,visg,D_visg,d4,D_d4,ndof)
  endif

!--Now form the viscosity-weighted values--------------------------------------

  call formMixedDen2(viso,viss,visotl,deno,dens ,denotl,    &
                     getDerivs,ndof,                        &
                     D_viso,D_viss,D_visotl,D_deno,D_dens,  &
                     D_denotl)
  call formMixedDen2(visg,viss,visgtl,deng,dens ,dengtl,    &
                     getDerivs,ndof,                        &
                     D_visg,D_viss,D_visgtl,D_deng,D_dens,  &
                     D_dengtl)
  call formMixedDen2(viss,visog,visstl,dens,denog,denstl,   &
                     getDerivs,ndof,                        &
                     D_viss,D_visog,D_visstl,D_dens,D_denog,&
                     D_denstl)

#if 0
! we may wish to uncomment and have as an optional error check at some point; 
! a lot can go wrong when the system gets to marginal saturations
if (getDerivs) then
  do i = 1,ndof
    if (isnan(D_denotl(i))) then
      print *, "denotl nan deriv at ", i
    endif
    if (isnan(D_dengtl(i))) then
      print *, "dengtl nan deriv at ", i
    endif
    if (isnan(D_denstl(i))) then
      print *, "denstl nan deriv at ", i
    endif
  enddo
endif
#endif

end subroutine TL4PDensity

! ************************************************************************** !

subroutine TL4PMiscibilityFractionFromSaturationFraction(fs,fm,dfmdfs)

!------------------------------------------------------------------------------
! Obtain the miscibility fraction from the solvent in vapour fraction
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  implicit none

  PetscReal,intent(in) ::fs
  PetscReal,intent(out)::fm
  PetscReal,intent(out)::dfmdfs
  PetscReal::den,deni

!--Default to straight line----------------------------------------------------

  fm    =fs
  dfmdfs=1.0

!--Check for special cases or ramp function------------------------------------

  if( fmis_is_zero ) then
    fm    =0.0
    dfmdfs=0.0
  else if( fmis_is_unity ) then
    fm    =1.0
    dfmdfs=0.0
  else
    if( fs<=fmis_sl ) then
      fm    =0.0d0
      dfmdfs=0.0

      ! constant slope at endpoint:
      if (fs == fmis_sl) then
        den=fmis_su-fmis_sl
        deni=0.0
        if( abs(den)>0.0 ) deni=1.0/(fmis_su-fmis_sl)
        dfmdfs=deni
      endif

    else if( fs>=fmis_su ) then
      fm    =1.0d0
      dfmdfs=0.0

      ! constant slope at endpoint:
      if (fs == fmis_su) then
        den=fmis_su-fmis_sl
        deni=0.0
        if( abs(den)>0.0 ) deni=1.0/(fmis_su-fmis_sl)
        dfmdfs=deni
      endif

    else
      den=fmis_su-fmis_sl
      deni=0.0
      if( abs(den)>0.0 ) deni=1.0/(fmis_su-fmis_sl)
      fm    =(fs-fmis_sl)*deni
      dfmdfs=deni
    endif
  endif

end subroutine TL4PMiscibilityFractionFromSaturationFraction

! ************************************************************************** !

subroutine TL4PScaleCriticals(sgcr,sogcr,fm,soil,svap,uoil,uvap &
                              ,getDerivs,ndof &
                              ,dof_soil,dof_sgas,dof_ssol &
                              ,D_fm,D_uoil,D_uvap)

!------------------------------------------------------------------------------
! Obtain lookup args modified for scaling (to be used against unscaled data)
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  implicit none

  PetscReal,intent(in )                   :: sgcr,sogcr,fm,soil,svap
  PetscBool,intent(in)                    :: getDerivs
  PetscInt,intent(in)                     :: ndof,dof_soil,dof_sgas,dof_ssol
  PetscReal,dimension(1:ndof),intent(in)  :: D_fm


  PetscReal,intent(out)                   :: uoil,uvap
  PetscReal,dimension(1:ndof),intent(out) :: D_uoil,D_uvap

  PetscReal                               :: ugcr,uogcr,fi
  PetscReal                               :: duvap_dugcr,duoil_duogcr
  PetscReal                               :: duoil_soil,duvap_svap
  PetscReal,dimension(1:ndof) :: D_fi,D_ugcr,D_uogcr

  fi=1.0-fm

  ugcr =fi*sgcr
  uogcr=fi*sogcr

  if (getDerivs) then
    D_ugcr  = -1.d0*D_fm*sgcr
    D_uogcr = -1.d0*D_fm*sogcr
  endif

  call TL4PScaleLookupSaturation(sgcr ,ugcr ,svap,uvap,duvap_svap,duvap_dugcr)
  call TL4PScaleLookupSaturation(sogcr,uogcr,soil,uoil,duoil_soil,duoil_duogcr)

  if (getDerivs) then
    D_uoil = 0.d0
    D_uoil = duoil_duogcr*D_uogcr
    D_uoil(dof_soil) = D_uoil(dof_soil) + duoil_soil

    D_uvap = 0.d0
    D_uvap = duvap_dugcr*D_ugcr
    D_uvap(dof_sgas) = D_uvap(dof_sgas) + duvap_svap
    D_uvap(dof_ssol) = D_uvap(dof_ssol) + duvap_svap
  endif

end subroutine TL4PScaleCriticals

! ************************************************************************** !

!subroutine TL4PScaleLookupSaturation(sl,ul,s,u)
subroutine TL4PScaleLookupSaturation(sl,ul,s,u,du_s,du_ul)


!------------------------------------------------------------------------------
! Transform a true saturation to a value to be used with unscaled data
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  implicit none

  PetscReal,intent(in ) :: sl,ul,s
  PetscReal,intent(out) :: u
  PetscReal             :: slc
  PetscReal             :: sli,slci

  PetscReal,intent(out) :: du_s,du_ul

!--Form safe inverses----------------------------------------------------------

  sli=0.0
  if( sl>0.0 ) sli=1.0/sl

  slc=1.0-sl

  slci=0.0
  if( slc>0.0 ) slci=1.0/slc

  if( s<=sl ) then
! First interval: scale by ul/sl so s=sl maps to ul
    u=s*ul*sli

    du_s = ul*sli
    du_ul = s*sli
  else
! Second interval: scale so that s=sl maps to ul, s=1 maps to 1
    u=ul+(s-sl)*(1.0-ul)*slci

    du_s = (1.0-ul)*slci
    du_ul = 1.0 - (s-sl)*slci
  endif

end subroutine TL4PScaleLookupSaturation

! ************************************************************************** !

subroutine formMixedDen2(visa,visb,visatl,dena,denb,denatl, &
                         getDerivs,ndof,                    &
                         D_visa,D_visb,D_visatl,D_dena,D_denb,       &
                         D_denatl)

!------------------------------------------------------------------------------
! Obtain interpolated Todd-Longstaff density to match Todd-Longstaff viscosity
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Apr 2018
!------------------------------------------------------------------------------

  use Derivatives_utilities_module
  implicit none

  PetscReal,intent(in)   :: visa,visb,visatl,dena,denb
  PetscReal,intent(inout):: denatl
  PetscReal              :: visbi,visatli
  PetscReal              :: m,me,mqp,meqp,mqpm1,mqpm1i
  PetscReal              :: fa,fb

  PetscBool              :: getDerivs
  PetscInt               :: ndof

  PetscReal,dimension(1:ndof),intent(in) :: D_visa,D_visb,D_visatl,D_dena,D_denb
  PetscReal,dimension(1:ndof),intent(inout) :: D_denatl

  PetscReal,dimension(1:ndof) :: D_m,D_me,D_fa,D_mqp,D_meqp

  PetscReal :: workerReal

!--Form the viscosity ratio and check that it is not unity---------------------

  visbi=0.0
  if( visb  >0.0 ) visbi=1.0/visb

  visatli=0.0
  if( visatl>0.0 ) visatli=1.0/visatl

  m=visa*visbi

  if (getDerivs) then
    D_m = DivRule(visa,D_visa,visb,D_visb,ndof)
  endif

  if( (m.ne.1.0) .and. (visatl>0.0) ) then

! Build the interpolation coefficient, fa

    me=visa*visatli

    if (getDerivs) then
      D_me = DivRule(visa,D_visa,visatl,D_visatl,ndof)
    endif

    mqp =M**0.25
    meqp=Me**0.25

    if (getDerivs) then
      workerReal = 0.25
      D_mqp = PowerRule(m,D_m,workerReal,ndof)
      D_meqp = PowerRule(me,D_me,workerReal,ndof)
    endif

    mqpm1=mqp-1.0

    mqpm1i=0.0
    if( abs(mqpm1)>0.0 ) mqpm1i=1.0/mqpm1

    fa=(mqp-meqp)*mqpm1i

! Check for limiting to (0,1)

    if( fa<0.0 ) fa=0.0
    if( fa>1.0 ) fa=1.0

    if (getDerivs) then
      if (fa >= 0.0 .AND. fa <= 1.0) then
        D_fa = DivRule(mqp-meqp,D_mqp-D_meqp,mqpm1,D_mqp,ndof)
      else
        D_fa = 0.d0
      endif
    endif

    fb=1.0-fa

! Construct output density and over-write existing default value

    denatl=dena*fa+denb*fb

    if (getDerivs) then
      ! i.e. denatl = dena*fa+denb(1-fa)
      !             = fa*(dena-denb) + denb
      D_denatl = ProdRule(fa,D_fa,dena-denb,D_dena-D_denb,ndof) &
               + D_denb
    endif

  endif

end subroutine formMixedDen2

! ************************************************************************** !

subroutine TOWGAuxFieldVolRefAve(this,grid,material,imat,option)
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 08/24/18
  !

  use Option_module
  use Grid_module
  use Material_Aux_class
  use Well_Data_class,only : SetFieldData

  implicit none

  class(pm_towg_aux_type) :: this
  type(grid_type) :: grid
  type(material_type), pointer :: material
  PetscInt, intent(in) :: imat(:)
  type(option_type) :: option

  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: local_id, ghosted_id
  PetscReal :: rpv_sh_l, rpv_sh_g,rpv_l, rpv_g
  PetscReal :: pr,pv,sh
  PetscReal :: fhpav_l,fhpav_g,fpav_l,fpav_g,fpav
  PetscInt :: int_mpi
  PetscErrorCode :: ierr

!  Initial values for sums (with and without hydrocarbon fraction)

  rpv_sh_l = 0.0d0
  rpv_sh_g = 0.0d0
  rpv_l   = 0.0d0
  rpv_g   = 0.0d0

  fhpav_l = 0.0d0
  fhpav_g = 0.0d0
  fpav_l  = 0.0d0
  fpav_g  = 0.0d0

!  Do the summations on this rank

  material_auxvars => material%auxvars

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (imat(ghosted_id) <= 0) cycle

    pv = material_auxvars(ghosted_id)%porosity_base * &
         material_auxvars(ghosted_id)%volume
    sh = this%auxvars(ZERO_INTEGER,ghosted_id)%sat(option%oil_phase) + &
          this%auxvars(ZERO_INTEGER,ghosted_id)%sat(option%gas_phase)
    if (this%IsSolventModel()) then
       sh = sh + this%auxvars(ZERO_INTEGER,ghosted_id)% &
                                           sat(option%solvent_phase)
    end if
    pr=this%auxvars(ZERO_INTEGER,ghosted_id)%pres(option%oil_phase)

    rpv_sh_l = rpv_sh_l + sh * pv
    rpv_l    = rpv_l    +      pv
    fhpav_l  = fhpav_l  + sh * pv * pr
    fpav_l   = fpav_l   +      pv * pr

  end do

!  Do the summations over ranks

  int_mpi = ONE_INTEGER
  call MPI_AllReduce(rpv_sh_l,rpv_sh_g, &
                     int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
  call MPI_AllReduce(fhpav_l,fhpav_g, &
                     int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
  call MPI_AllReduce(rpv_l,rpv_g, &
                     int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
  call MPI_AllReduce(fpav_l,fpav_g, &
                     int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

! Form the average pressure: use simpler form if no hydrocarbon volume

  if ( rpv_sh_g>0.0 ) then
    fpav = fhpav_g/rpv_sh_g  !  Hydrocarbon volume exists
  else if ( rpv_g>0.0 ) then
    fpav = fpav_g/rpv_g      ! No hydrocarbon, but pore volume exists
  endif

!  Set value

  call SetFieldData(fpav)

  nullify(material_auxvars)

end subroutine TOWGAuxFieldVolRefAve

! ************************************************************************** !

subroutine TOWGGetLocalSol(this,grid,material,imat,option,vsoll,isol,zsol)
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 08/24/18

  use Grid_module
  use Material_Aux_class
  use Option_module
  use String_module,only : StringCompareIgnoreCase

  implicit none

  class(pm_towg_aux_type) :: this
  type(grid_type) :: grid
  type(material_type), pointer :: material
  PetscInt, intent(in) :: imat(:)
  type(option_type) :: option
  PetscReal :: vsoll(:,:)
  PetscInt,intent(in)::isol
  character(len=8)::zsol
  PetscInt::iphase
  PetscBool::ispres
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: local_id, ghosted_id

  material_auxvars => material%auxvars

  ispres=PETSC_FALSE
  iphase=0

  if( StringCompareIgnoreCase(zsol,'Pressure') ) ispres=PETSC_TRUE
  if( StringCompareIgnoreCase(zsol,'Soil')     ) iphase=option%oil_phase
  if( StringCompareIgnoreCase(zsol,'Sgas')     ) iphase=option%gas_phase
  if( StringCompareIgnoreCase(zsol,'Swat')     ) iphase=option%liquid_phase

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (imat(ghosted_id) <= 0) cycle
    if( isPres ) then
      vsoll(local_id,isol) = this%auxvars(ZERO_INTEGER,ghosted_id)%pres(option%oil_phase)
    else if ( iphase>0 ) then
      vsoll(local_id,isol) = this%auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase)
    else
      vsoll(local_id,isol)=0.0
    endif
  end do

end subroutine TOWGGetLocalSol

! ************************************************************************** !

function IsSolventModel(this)

  implicit none

  class(pm_towg_aux_type) :: this

  PetscBool :: IsSolventModel

  IsSolventModel = PETSC_FALSE

  if (towg_miscibility_model == TOWG_SOLVENT_TL) IsSolventModel = PETSC_TRUE

end function IsSolventModel

end module PM_TOWG_Aux_module
