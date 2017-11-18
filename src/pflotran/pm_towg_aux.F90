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

  PetscBool, public :: towg_isothermal = PETSC_FALSE
  !PO:needs to add input for towg_no_gas and towg_no_oil in pm_towg%Read
  !   towg_no_oil currently supported only for TOWG_IMMISCIBLE
  !   and TOWG_TODD_LONGSTAFF. To have it working also for BLACK_OIL and SOLV.
  !   must swap the orger of primary vars for TOWG_LIQ_GAS_STATE,
  !   Po,Sg,Xg^G -> Po,Xg^G,Sg
  PetscBool, public :: towg_no_gas = PETSC_FALSE
  PetscBool, public :: towg_no_oil = PETSC_FALSE
  PetscInt, public :: towg_miscibility_model = UNINITIALIZED_INTEGER

  ! list of TOWG paramters
  PetscInt, parameter, public :: TOWG_PREV_TS = 1
  PetscInt, parameter, public :: TOWG_PREV_IT = 2

  ! available miscibility models - now defined in PFLOTRAN_Constants_module
  ! PetscInt, parameter, public :: TOWG_IMMISCIBLE = 1
  ! PetscInt, parameter, public :: TOWG_TODD_LONGSTAFF = 2
  ! PetscInt, parameter, public :: TOWG_BLACK_OIL = 3
  ! PetscInt, parameter, public :: TOWG_SOLVENT_TL = 4

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
  PetscInt, parameter, public :: TOWG_BUBBLE_POINT_3PH_DOF   = 3 !Variable substitutes with gas saturation DKP
  PetscInt, parameter, public :: TOWG_X_GAS_IN_OIL_DOF = 3
  PetscInt, parameter, public :: TOWG_X_OIL_IN_GAS_DOF = 3
  PetscInt, parameter, public :: TOWG_SOLV_SATURATION = 4
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
! DKP Note twin values for Sg and Pb, TOWG_X_GAS_IN_GAS_INDEX removed
  PetscInt, parameter, public :: TOWG_OIL_PRESSURE_INDEX = 1
  PetscInt, parameter, public :: TOWG_OIL_SATURATION_INDEX = 2
  PetscInt, parameter, public :: TOWG_GAS_SATURATION_INDEX = 3
  PetscInt, parameter, public :: TOWG_BUBBLE_POINT_INDEX   = 3
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
  PetscInt, parameter, public :: TOWG_MAX_INDEX = 14

  !Indices used to map aux_int_var for condition values
  PetscInt, parameter, public :: TOWG_STATE_INDEX = 1

  !flags to identify type of auxvar update
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_BOUNDARY = 2

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
            TOWGBlackOilAuxVarComputeSetup
  !          TOilImsAuxVarPerturb, TOilImsAuxDestroy, &
  !          TOilImsAuxVarStrip

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

  implicit none

  type(option_type) :: option

  class(pm_towg_aux_type), pointer :: TOWGAuxCreate

  class(pm_towg_aux_type), pointer :: aux

  allocate(towg_fmw_comp(option%nflowspec))

  !in TOWG the gas FMW must be defined in the input deck
  if ( Uninitialized(EOSGasGetFMW()) ) then
    option%io_buffer = 'TOWG: gas FMW not initialised. ' // &
                       'Define its value in the the input deck' // &
                       ' or add EOS GAS card to default to FMWAIR'
    call printErrMsg(option)
  endif

  towg_fmw_comp(1) = FMWH2O
  towg_fmw_comp(2) = EOSOilGetFMW()
  towg_fmw_comp(3) = EOSGasGetFMW()

  !need to add an EOS for solvent, before adding solvent
  !if ( towg_miscibility_model == TOWG_SOLVENT_TL ) then
  !  towg_fmw_comp(4) =
  !end if

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

  allocate(this%auxvars(0:option%nflowdof,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, option%nflowdof
      call this%auxvars(idof,ghosted_id)%Init(option)
      if (towg_miscibility_model == TOWG_TODD_LONGSTAFF ) then
        call this%auxvars(idof,ghosted_id)%InitTL(option)
      end if
      if (towg_miscibility_model == TOWG_BLACK_OIL      ) then
        call this%auxvars(idof,ghosted_id)%InitBO(option)
      end if
    enddo
  enddo

  this%num_aux = grid%ngmax

  if (num_bc_connection > 0) then
    allocate(this%auxvars_bc(num_bc_connection))
    do iconn = 1, num_bc_connection
      call this%auxvars_bc(iconn)%Init(option)
      if (towg_miscibility_model == TOWG_TODD_LONGSTAFF ) then
        call this%auxvars_bc(iconn)%InitTL(option)
      end if
      if (towg_miscibility_model == TOWG_BLACK_OIL      ) then
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
      if (towg_miscibility_model == TOWG_BLACK_OIL      ) then
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

  PetscInt :: lid, oid, gid
  PetscReal :: cell_pressure, wat_sat_pres
  PetscReal :: krl, visl, dkrl_Se
  PetscReal :: kro, viso, dkro_Se
  PetscReal :: krg, visg, dkrg_Se
  PetscReal :: sat_liq_gas, sat_tot_liq
  PetscReal :: dummy, dummy2
  !PetscReal :: Uoil_J_kg, Hoil_J_kg
  PetscErrorCode :: ierr

  ! from SubsurfaceSetFlowMode
  ! option%liquid_phase = 1           ! liquid_pressure - water pressure
  ! option%oil_phase = 2              ! oil_pressure
  ! option%gas_phase = 3              ! gas_pressure
  lid = option%liquid_phase
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

  !assing auxvars given by the solution variables
  auxvar%pres(oid) = x(TOWG_OIL_PRESSURE_DOF)
  auxvar%sat(oid) = x(TOWG_OIL_SATURATION_DOF)
  auxvar%sat(gid) = x(TOWG_GAS_SATURATION_3PH_DOF)
  auxvar%temp = x(towg_energy_dof)

  auxvar%sat(lid) = 1.d0 - auxvar%sat(oid) - auxvar%sat(gid)

  !PO before characteristic function refactoring - to match TOUGH2-EOS8
  ! sat_tot_liq = auxvar%sat(lid) + auxvar%sat(oid)
  ! call characteristic_curves%saturation_function% &
  !      CapillaryPressure(sat_tot_liq,auxvar%pc(lid),dummy,option)
  !
  ! auxvar%pres(lid) = auxvar%pres(oid)
  ! auxvar%pres(gid) = auxvar%pres(lid) + auxvar%pc(lid)
  ! auxvar%pc(oid) = 0.0d0
  !end before refactoring

  !PO after characteristic function refactoring
  call characteristic_curves%oil_wat_sat_func% &
            CapillaryPressure(auxvar%sat(oid),auxvar%sat(gid),auxvar%pc(lid), &
                              dummy,dummy2,option)
  call characteristic_curves%oil_gas_sat_func% &
            CapillaryPressure(auxvar%sat(oid),auxvar%sat(gid),auxvar%pc(oid), &
                              dummy,dummy2,option)
  auxvar%pres(lid) = auxvar%pres(oid) - auxvar%pc(lid)
  auxvar%pres(gid) = auxvar%pres(oid) + auxvar%pc(oid)
  !end after refactoring

  cell_pressure = max(auxvar%pres(lid),auxvar%pres(oid),auxvar%pres(gid))


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
  ! using cell_pressure (which is the max press)? or %pres(lid)?
  call EOSWaterDensity(auxvar%temp,cell_pressure, &
                       auxvar%den_kg(lid),auxvar%den(lid),ierr)
  call EOSWaterEnthalpy(auxvar%temp,cell_pressure,auxvar%H(lid),ierr)
  auxvar%H(lid) = auxvar%H(lid) * 1.d-6 ! J/kmol -> MJ/kmol

  ! MJ/kmol comp                  ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
  auxvar%U(lid) = auxvar%H(lid) - (cell_pressure / auxvar%den(lid) * 1.d-6)

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
  call characteristic_curves%liq_rel_perm_function% &
         RelativePermeability(auxvar%sat(lid),krl,dkrl_Se,option)

  call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,ierr)

  ! use cell_pressure; cell_pressure - psat calculated internally
  call EOSWaterViscosity(auxvar%temp,cell_pressure,wat_sat_pres,visl,ierr)

  auxvar%mobility(lid) = krl/visl

  ! compute oil mobility (rel. perm / viscostiy)
  sat_liq_gas = auxvar%sat(lid) + auxvar%sat(gid)
  call characteristic_curves%oil_rel_perm_function% &
         RelativePermeability(sat_liq_gas,kro,dkro_Se,option)

  call EOSOilViscosity(auxvar%temp,auxvar%pres(oid), &
                       auxvar%den(oid), viso, ierr)

  auxvar%mobility(oid) = kro/viso

  !compute gas mobility (rel. perm / viscosity)
  sat_tot_liq = auxvar%sat(lid) + auxvar%sat(oid)
  call characteristic_curves%gas_rel_perm_function% &
         RelativePermeability(sat_tot_liq,krg,dkrg_Se,option)

  !currently only a viscosity model for air or constant value
  call EOSGasViscosity(auxvar%temp,auxvar%pres(gid), &
                       auxvar%pres(gid),auxvar%den(gid),visg,ierr)

  auxvar%mobility(gid) = krg/visg


end subroutine TOWGImsAuxVarCompute

! DKP New routines to find oil phase mole fractions as a function of (P,Pb)---

subroutine getRsVolume(bubble_point,rs_volume)

!------------------------------------------------------------------------------
! Obtain the Rs value (as surface vol gas)/(surface volume oil)
! This is a function for now - needs replacement with table
! Taken roughly from SPE1 (Rs~0.8*Pb(Bars)=0.8*1.0d-5*P(Pa)
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Jul 2017
!------------------------------------------------------------------------------

  PetscReal, intent(in ) :: bubble_point
  PetscReal, intent(out) :: rs_volume

  rs_volume=1.0+0.8*1.0d-5*bubble_point !DKPDKP This is dummy value for tests

end subroutine getRsVolume

subroutine getBlackOilComposition(bubble_point,temperature,table_idxs,&
                                  xo,xg,pref,tref)

  use EOS_Oil_module
  use EOS_Gas_module

  PetscReal, intent(in ) :: bubble_point
  PetscReal, intent(in ) :: temperature
  PetscInt, pointer, intent(inout) :: table_idxs(:)
  PetscReal, intent(out) :: xo
  PetscReal, intent(out) :: xg
  PetscReal              :: pref,tref
  PetscReal              :: rs_volume=0.0d0,rs_molar
  PetscReal              :: mdenrefo,mdenrefg,denrefo,denrefg,mwo,mwg,ideal_gas_molar_density

  PetscInt               :: ierr

!--Get molecular weights of oil and gas----------------------------------------

  mwo=EOSOilGetFMW()
  mwg=EOSGasGetFMW()

!--Get surface mass densities of oil and gas components------------------------

  !denrefo=EOSOilGetDenRef()
  !PO this has been changed to read a reference/surface density independently
  !   from the one defined within the InverseLiean model
  !   In the input deck the user can now enter one of the following
  !   REFERENCE_DENSITY 1000 kg/m^3
  !   SURFACE_DENSITY 1000 kg/m^3
  !   STANDARD_DENSITY 1000 kg/m^3
  ! EOSOilGetReferenceDensity() return this value
  ! EOSOilGetDenRef() has been replaced by EOSOilGetReferenceDensity()
  !... sorry for the long name
  denrefo=EOSOilGetReferenceDensity()

  ideal_gas_molar_density=pref/(IDEAL_GAS_CONSTANT*(273.15+tref))*1.d-3

  denrefg=mwg*ideal_gas_molar_density !DKPDKP ideal gas form needs work - ie PVTG table
  !PO Now can get the gas reference value using EOSGasGetReferenceDensity
  !   this function returns the value defined in the input deck by one of the following:
  !   REFERENCE_DENSITY 0.678366 kg/m^3
  !   SURFACE_DENSITY 0.678366 kg/m^3
  !   STANDARD_DENSITY 0.678366 kg/m^3
  denrefg = EOSGasGetReferenceDensity()
  ! PO if instead you need the density defined by the PVDG table,
  !    then call the EOSGasDensity.
  !    Note, if a call to gas density is needed here - consider passing
  !    gas_den to getBlackOilComposition, and call EOSGasDensityEnergy before
  !    getBlackOilComposition in TOWGBlackOilAuxVarCompute
  !call EOSGasDensity(temperature,bubble_point,denrefg,ierr,table_idxs)

!--Molar density = (mass density)/MW = (Kg/sm3)/(Kg/KG-mol) = KG-mol/sm3-------

  mdenrefo=denrefo/Mwo
  mdenrefg=denrefg/Mwg

!--Get GOR as a surface volume ratio Rsv=(vol gas)/(vol oil)-------------------

  call getRsVolume(bubble_point,rs_volume)
  !PO - lookup table now available - EOSOilRS returns mol_gas/mol_oil
  !call EOSOilRS(temperature,bubble_point,rs_molar,ierr,table_idxs)

!--Convert to molar GOR, Rsm---------------------------------------------------
!
!  Rsm=(moles gas)/(moles oil)=(mol den gas)*(vol gas)/((mol den oil)*(vol oil))
!                             =Rsv*(mol den gas)/(mol den oil)
!
!------------------------------------------------------------------------------

  rs_molar=rs_volume*(mdenrefg/mdenrefo)

!--Get oil and gas mole fractions----------------------------------------------
!
! xo=(moles oil)/(moles oil+moles gas)=1/(1+(moles gas)/(moles oil))=1/(1+Rsm)
! xg=1-xo=1-1/(1+Rsm)=(1+Rsm-1)/(1+Rsm)=Rsm/(1+Rsm)
!
!------------------------------------------------------------------------------

  xo=1.0d0   /(1.0d0+rs_molar)
  xg=rs_molar/(1.0d0+rs_molar)

end subroutine getBlackOilComposition

! DKP---------------------------------------------------------------------------

! DKP New black oil auxillary variable calculation------------------------------

subroutine TOWGBlackOilAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, & !This updates ths black oil variables
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
  !PO this must not be used because belonging to the general mode
  !   when the phase state integer name is needed one can use those specific to
  !   TOWG: TOWG_LIQ_OIL_STATE, TOWG_LIQ_GAS_STATE, TOWG_THREE_PHASE_STATE, etc
  !   This what is actually done here, so I commented the line below -to be deleted
  !use General_Aux_module, only : LIQUID_STATE, TWO_PHASE_STATE

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar ! passing this for salt conc.
                                            ! not currenty used
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  PetscInt :: natural_id !only for debugging/print out - currently not used

  PetscInt :: lid, oid, gid, istate
  PetscReal :: cell_pressure, wat_sat_pres
  PetscReal :: krl, visl, dkrl_Se
  PetscReal :: kro, viso, dkro_Se
  PetscReal :: krg, visg, dkrg_Se
  PetscReal :: sat_liq_gas, sat_tot_liq
  PetscReal :: dummy,So,Sg,Sw,Pb,So1,Sg1,Sw1,Pb1
  !PetscReal :: Uoil_J_kg, Hoil_J_kg
  PetscErrorCode :: ierr
  PetscBool isSat
  PetscReal :: epss=1.0d-4,epssc,epsp=1.0d3,pref,tref,pressure

  epssc=1.0d0-epss

  ! from SubsurfaceSetFlowMode
  ! option%liquid_phase = 1           ! liquid_pressure
  ! option%oil_phase = 2              ! oil_pressure
  ! option%gas_phase = 3              ! gas_pressure
  lid = option%liquid_phase
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

!--Check state to see if saturated----------------------------------------------

  istate=global_auxvar%istate

  if( istate==TOWG_THREE_PHASE_STATE ) then
    isSat=.true.
  else
    isSat=.false.
  endif

!-------------------------------------------------------------------------------
! Getting auxvars as given by the solution variables,
! allowing for saturated/undersaturated switch
!-------------------------------------------------------------------------------

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

  pressure=auxvar%pres(oid)
  Pb      =auxvar%bo%bubble_point

!------------------------------------------------------------------------------
! Check if this state still valid and flip if not (but not on diff call)
!------------------------------------------------------------------------------

  So=auxvar%sat(oid)
  Sg=auxvar%sat(gid)
  Sw=auxvar%sat(lid)
  Pb=auxvar%bo%bubble_point

  if (option%iflag /= TOWG_UPDATE_FOR_DERIVATIVE) then
    if( isSat ) then
      if( Sg<0.0d0 ) then
        global_auxvar%istate          =TOWG_LIQ_OIL_STATE
        auxvar%sat(gid)               =0.0d0
        auxvar%bo%bubble_point        =auxvar%pres(oid)-epsp
        x(TOWG_BUBBLE_POINT_3PH_DOF)  =auxvar%pres(oid)-epsp
      endif
    else
      if( auxvar%bo%bubble_point.gt.auxvar%pres(oid) ) then
        global_auxvar%istate          =TOWG_THREE_PHASE_STATE
        auxvar%bo%bubble_point        =auxvar%pres(oid)
        auxvar%sat(gid)               =epss
        x(TOWG_GAS_SATURATION_3PH_DOF)=epss
        if( auxvar%sat(oid).gt.epssc ) then
          auxvar%sat(oid)             =epssc
          x(TOWG_OIL_SATURATION_DOF)  =epssc
        endif
      endif
    endif
  endif

  So1=auxvar%sat(oid)
  Sg1=auxvar%sat(gid)
  Sw1=auxvar%sat(lid)
  Pb1=auxvar%bo%bubble_point

!------------------------------------------------------------------------------
! Set up the final saturation value (water)
!------------------------------------------------------------------------------

  auxvar%sat(lid) = 1.d0 - auxvar%sat(oid) - auxvar%sat(gid)

!--Look up the Rs value--------------------------------------------------------

  pref=option%reference_pressure
  tref=option%reference_temperature
  !PO either you pass the entire auxvar, or the you also need to pass:
  !    the table indices, and the temprature (on which RS might depend on)
  call getBlackOilComposition(auxvar%bo%bubble_point,auxvar%temp, &
                              auxvar%eos_table_idx,auxvar%bo%xo,auxvar%bo%xg, &
                              pref,tref)

!--Get the capillary pressures-------------------------------------------------

  !below pc_ow /= 0
  !compute capillary presssure water/oil (pc_ow)
  !call characteristic_curves%saturation_function% &
  !     CapillaryPressure(material_auxvar,auxvar%sat(lid),auxvar%pc(lid),option)
  !auxvar%pres(lid) = auxvar%pres(oid) - auxvar%pc(lid)

  !Assumptions below on capillary pressure for comparison with TOUGH2-EOS8:
  ! pc_ow = 0, only pc_gw /= 0 and computed considering water saturation only
  ! this results in pc_gw = pc_go

  !assuming no capillary pressure between oil and water: pc_ow = 0
  ! pc_ow = 0.0d0
  auxvar%pres(lid) = auxvar%pres(oid)

  !compute capillary pressure gas/water (pc_gw),
  ! pc_go = pc_gw when pc_ow = 0
  !call characteristic_curves%saturation_function% &
  !     CapillaryPressure(material_auxvar,auxvar%sat(lid),auxvar%pc(lid),option)
  !To match TOUGH-EOS8, consider water plu poil as wetting phase
  !for capillary press computation
  sat_tot_liq = auxvar%sat(lid) + auxvar%sat(oid)
  call characteristic_curves%saturation_function% &
       CapillaryPressure(sat_tot_liq,auxvar%pc(lid),dummy,option)

  auxvar%pres(gid) = auxvar%pres(lid) + auxvar%pc(lid)

  auxvar%pc(oid) = 0.0d0


  cell_pressure = max(auxvar%pres(lid),auxvar%pres(oid),auxvar%pres(gid))

!--Get rock and fluid properties-----------------------------------------------

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

!------------------------------------------------------------------------------
!DKPDKP Thermodynamic properties for all phases in black oil model
! Will need more work - eg total molar density as function of P and Pb
!------------------------------------------------------------------------------

! Water phase thermodynamic properties

  call EOSWaterDensity(auxvar%temp,cell_pressure, &
                       auxvar%den_kg(lid),auxvar%den(lid),ierr)
  call EOSWaterEnthalpy(auxvar%temp,cell_pressure,auxvar%H(lid),ierr)
  auxvar%H(lid) = auxvar%H(lid) * 1.d-6 ! J/kmol -> MJ/kmol

  ! MJ/kmol comp                  ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
  auxvar%U(lid) = auxvar%H(lid) - (cell_pressure / auxvar%den(lid) * 1.d-6)

! Gas phase thermodynamic properties (also needed for dissolved gas)

  !compute gas properties (default is air with an analytical function)
  !PO: if the PVDG table is defined for EOS GAS the molar density
  !   is computed from BG (i.e. FVF for the gas) and the SURFACE_DENSITY.
  !   The transormation in done once after the EOS GAS input has been read.
  !   Note that auxvar%eos_table_idx for table lookup is passed after ierr
  !   because it is an optional argument (needed only for table lookup).
  !   It might be worth calling this before getBlackOilComposition - so that
  !  the gas density is already avaiable
  call EOSGasDensityEnergy(auxvar%temp,auxvar%pres(gid),auxvar%den(gid), &
                          auxvar%H(gid),auxvar%U(gid),ierr,auxvar%eos_table_idx)

  auxvar%den_kg(gid) = auxvar%den(gid) * EOSGasGetFMW()

  auxvar%H(gid) = auxvar%H(gid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(gid) = auxvar%U(gid) * 1.d-6 ! J/kmol -> MJ/kmol

! Oil phase thermodynamic properties currently not including Pb dependence--

!  Note:see comment 'ADD HERE BRINE dependency...' elsewhere

  !PO: if either PVDO or PVCO tables is defined for EOS OIL, the molar density
  !    is computed from BO (i.e. FVF for the oil) and SURFAACE_DENSITY.
  !    The transormation in done once after the EOS OIL input has been read.
  !    Note that auxvar%eos_table_idx for table lookup is passed after ierr
  !    because it is an optional argument (needed only if table lookup)
  call EOSOilDensityEnergy(auxvar%temp,auxvar%pres(oid),&
                           auxvar%den(oid),auxvar%H(oid), &
                           auxvar%U(oid),ierr,auxvar%eos_table_idx)

! Correct oil phase molar density and enthalpy for oil composition----------

! EOSOilDensity returns oil moles/volume in oil phase,
! but really have (1+Rsmolar) times as many moles when dissolved gas included
! 1+Rsmolar=1+xg/xo=(xo+xg)/xo=1/xo

  auxvar%den(oid)=auxvar%den(oid)/auxvar%bo%xo

! Get oil mass density as (mixture oil molar density).(mixture oil molecular weight)
  auxvar%den_kg(oid) = auxvar%den(oid) * ( auxvar%bo%xo*EOSOilGetFMW() &
                                          +auxvar%bo%xg*EOSGasGetFMW() )

  auxvar%H(oid) = auxvar%H(oid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(oid) = auxvar%U(oid) * 1.d-6 ! J/kmol -> MJ/kmol

! We now have oil enthalpy/oil mole in oil phase (if calculation done for pure oil)
! but really have total molar enthalpy of:
! hydrocarbon enthalpy/hydrocarbon mole= xo.(oil enthalpy/oil mole+xg*(gas enthalpy/gas mole)

  auxvar%H(oid) = auxvar%bo%xo*auxvar%H(oid)+auxvar%bo%xg*auxvar%H(gid)
  auxvar%U(oid) = auxvar%bo%xo*auxvar%U(oid)+auxvar%bo%xg*auxvar%U(gid)

!--Remainder of the fluid mobility calculation is standard---------------------

  ! compute water mobility (rel. perm / viscosity)
  call characteristic_curves%liq_rel_perm_function% &
         RelativePermeability(auxvar%sat(lid),krl,dkrl_Se,option)

  call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,ierr)

  ! use cell_pressure; cell_pressure - psat calculated internally
  call EOSWaterViscosity(auxvar%temp,cell_pressure,wat_sat_pres,visl,ierr)

  auxvar%mobility(lid) = krl/visl

  ! compute oil mobility (rel. perm / viscosity)
  sat_liq_gas = auxvar%sat(lid) + auxvar%sat(gid)
  call characteristic_curves%oil_rel_perm_function% &
         RelativePermeability(sat_liq_gas,kro,dkro_Se,option)

  !PO - if PVDO or PVCO defined in EOS OIL, the viscosity values are extracted
  !     via table lookup.
  call EOSOilViscosity(auxvar%temp,auxvar%pres(oid), &
                       auxvar%den(oid), viso, ierr,auxvar%eos_table_idx)

  auxvar%mobility(oid) = kro/viso

  !compute gas mobility (rel. perm / viscosity)
  sat_tot_liq = auxvar%sat(lid) + auxvar%sat(oid)
  call characteristic_curves%gas_rel_perm_function% &
         RelativePermeability(sat_tot_liq,krg,dkrg_Se,option)

  !currently the only analytical viscosity model available is for air
  !PO - if PVDG is defined in EOS GAS, the viscosity values are extracted
  !     via table lookup.
  call EOSGasViscosity(auxvar%temp,auxvar%pres(gid), &
                       auxvar%pres(gid),auxvar%den(gid),visg,ierr,&
                       auxvar%eos_table_idx)

  auxvar%mobility(gid) = krg/visg

  !PO below example on how to call the fucntion computing the compressibility
  !   and the viscosibility
  !call EOSOilCompressibility(auxvar%temp,auxvar%pres(oid), &
  !                           dummy,ierr,auxvar%eos_table_idx)
  !
  !call EOSOilViscosibility(auxvar%temp,auxvar%pres(oid), &
  !                          dummy,ierr,auxvar%eos_table_idx)

end subroutine TOWGBlackOilAuxVarCompute

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

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar ! passing this for salt conc.
                                            ! not currenty used
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  PetscInt :: natural_id !only for debugging/print out - currently not used

  PetscInt :: lid, oid, gid
  PetscReal :: cell_pressure, wat_sat_pres
  PetscReal :: krl, dkrl_Se
  PetscReal :: krh, dkrh_Se
  PetscReal :: viso
  PetscReal :: visg
  PetscReal :: visl
  PetscReal :: sat_water
  PetscReal :: dummy,dummy2,deno,deng
  !PetscReal :: Uoil_J_kg, Hoil_J_kg
  PetscErrorCode :: ierr
  PetscReal :: krotl=0.0,krgtl=0.0,viscotl=0.0,viscgtl=0.0,denotl=0.0,dengtl=0.0
!--Get phase pointers for water,oil and gas-----------------------------------

  ! from SubsurfaceSetFlowMode
  ! option%liquid_phase = 1           ! liquid_pressure
  ! option%oil_phase = 2              ! oil_pressure
  ! option%gas_phase = 3              ! gas_pressure

  lid = option%liquid_phase
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

  auxvar%sat(lid) = 1.d0 - auxvar%sat(oid) - auxvar%sat(gid)

!--Set up phase pressures and capillary pressures-----------------------------

  ! In Todd-Longstaff mode assume no cap. pressure between oil and gas: pcog=0
  ! The cap. pressure is hydrocarbon/water (pchw) as function of water satn.

  auxvar%pc(oid) = 0.0d0

  sat_water = auxvar%sat(lid)

  ! PO before characteristic_curves refactoring
  !call characteristic_curves%saturation_function% &
  !     CapillaryPressure(sat_water,auxvar%pc(lid),dummy,option)

  ! PO after characteristic_curves refactoring
  !    dummy=dPcdSo,dummy2=dPcdSg
  call characteristic_curves%oil_wat_sat_func% &
            CapillaryPressure(auxvar%sat(oid),auxvar%sat(gid),auxvar%pc(lid), &
                              dummy,dummy2,option)

  auxvar%pres(gid) = auxvar%pres(oid)
  auxvar%pres(lid) = auxvar%pres(oid) - auxvar%pc(lid)

  cell_pressure = max(auxvar%pres(lid),auxvar%pres(oid),auxvar%pres(gid))

!--Calculate effective porosity as a function of pressure---------------------

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

!--Update thermodynamic properties (density, enphalpy..) for all phases-------

  ! Liquid phase thermodynamic properties
  ! using cell_pressure (which is the max press)? or %pres(lid)?

  call EOSWaterDensity(auxvar%temp,cell_pressure, &
                       auxvar%den_kg(lid),auxvar%den(lid),ierr)
  call EOSWaterEnthalpy(auxvar%temp,cell_pressure,auxvar%H(lid),ierr)
  auxvar%H(lid) = auxvar%H(lid) * 1.d-6 ! J/kmol -> MJ/kmol

  ! MJ/kmol comp                  ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
  auxvar%U(lid) = auxvar%H(lid) - (cell_pressure / auxvar%den(lid) * 1.d-6)

  ! ADD HERE BRINE dependency. Two options (see mphase)
  ! - salinity constant in space and time (passed in option%option%m_nacl)
  ! - salt can be transported by RT (sequential coupling) and passed
  !   and passed with global_auxvar%m_nacl

  call EOSOilDensityEnergy(auxvar%temp,auxvar%pres(oid),&
                           auxvar%den(oid),auxvar%H(oid), &
                           auxvar%U(oid),ierr,auxvar%eos_table_idx)

  auxvar%den_kg(oid) = auxvar%den(oid) * EOSOilGetFMW()

  auxvar%H(oid) = auxvar%H(oid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(oid) = auxvar%U(oid) * 1.d-6 ! J/kmol -> MJ/kmol

  !compute gas properties (default is air - but methane can be set up)
  call EOSGasDensityEnergy(auxvar%temp,auxvar%pres(gid),auxvar%den(gid), &
                      auxvar%H(gid),auxvar%U(gid),ierr,auxvar%eos_table_idx)

  auxvar%den_kg(gid) = auxvar%den(gid) * EOSGasGetFMW()
  auxvar%H(gid) = auxvar%H(gid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(gid) = auxvar%U(gid) * 1.d-6 ! J/kmol -> MJ/kmol

!--Set up water phase rel. perm. and unmixed viscosities-----------------------

  ! compute water relative permability Krw(Sw)
  call characteristic_curves%liq_rel_perm_function% &
         RelativePermeability(sat_water,krl,dkrl_Se,option)

  call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,ierr)

  ! use cell_pressure; cell_pressure - psat calculated internally
  call EOSWaterViscosity(auxvar%temp,cell_pressure,wat_sat_pres,visl,ierr)

  auxvar%mobility(lid) = krl/visl

!--Set up hydrocarbon phase rel. perm. and unmixed viscosities----------------

  ! In Todd-Longstaff case, look up the hydrocarbon rel perm using the
  ! hydrocarbon saturation (Sh=So+Sg=1-Sw), with water saturation as argument

  call characteristic_curves%oil_rel_perm_function% &
         RelativePermeability(sat_water,krh,dkrh_Se,option)

  ! Oil viscosity
  call EOSOilViscosity(auxvar%temp,auxvar%pres(oid), &
                       auxvar%den(oid), viso, ierr,auxvar%eos_table_idx)

  ! Gas viscosity : currently only viscosity model for air or constant value
  call EOSGasViscosity(auxvar%temp,auxvar%pres(gid),auxvar%pres(gid), &
                       auxvar%den(gid),visg,ierr,auxvar%eos_table_idx)

!--Set up the Todd-Longstaff rel perms, viscosities and densities-------------

  deno=auxvar%den_kg(oid)
  deng=auxvar%den_kg(gid)

  call vToddLongstaff( oid,gid,krh,viso,visg,deno,deng,auxvar &
                      ,krotl,krgtl,viscotl,viscgtl,denotl,dengtl )

!--Calculate and store oil and gas mobilities---------------------------------

  if( viscotl>0.0 ) then
    auxvar%mobility(oid) = krotl/viscotl
  else
    auxvar%mobility(oid) = 0.0
  endif

  if( viscgtl>0.0 ) then
    auxvar%mobility(gid) = krgtl/viscgtl
  else
    auxvar%mobility(gid) = 0.0
  endif

!--Store mixed density values-------------------------------------------------

  auxvar%tl%den_oil_eff_kg = denotl
  auxvar%tl%den_gas_eff_kg = dengtl

end subroutine TOWGTLAuxVarCompute

!==============================================================================

subroutine vToddLongstaff( oid,gid,krh,visco,viscg,deno,deng,auxvar   &
                          ,krotl,krgtl,viscotl,viscgtl,denotl,dengtl )

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

  implicit none

  PetscInt ,intent(in )   ::oid,gid
  PetscReal,intent(in )   ::krh
  PetscReal,intent(in )   ::visco,viscg
  PetscReal,intent(in )   ::deno,deng
  class(auxvar_towg_type) :: auxvar
  PetscReal,intent(out)   ::krotl,krgtl
  PetscReal,intent(out)   ::viscotl,viscgtl
  PetscReal,intent(out)   ::denotl,dengtl

  PetscReal::so,sg,den,deninv,fo,fg

!--Get oil and gas saturations and form fractions------------------------------

  so=auxvar%sat(oid)
  sg=auxvar%sat(gid)

  den=so+sg
  if( den>0.0 ) then
    deninv=1.0/den
    fo=so*deninv
    fg=sg*deninv
  else
    fo=0.5
    fg=0.5
  endif

!--Split the hydrocarbon relative permeability using fractions (TL,eqn. 2a,2b)-

  krotl=fo*krh
  krgtl=fg*krh

!--Form the omega-weighted oil and gas viscosities-----------------------------

  call vToddLongstaffViscosity(fo,fg,so,sg,visco,viscg,viscotl,viscgtl)

!--Form the omega weighted oil and gas densities (for Darcy flow only)---------

  call vToddLongstaffDensity( fo,fg,visco,viscg,viscotl,viscgtl &
                             ,deno,deng,denotl,dengtl)

end subroutine vToddLongstaff

!==============================================================================

subroutine vToddLongstaffViscosity(fo,fg,so,sg,visco,viscg,viscotl,viscgtl)

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
  implicit none

  PetscReal,intent(in ):: fo,fg,so,sg,visco,viscg
  PetscReal,intent(out):: viscotl,viscgtl
  PetscReal            :: tlomegac,sn,sni,viscqpo,viscqpg,wviscqp &
                         ,denom,denominv,viscm
  PetscReal            :: viscoimw,viscgimw,viscmw

!--Set up complement of the Todd Longstaff omega-------------------------------

  tlomegac=1.0-val_tl_omega

!--Form quarter-powers of basic viscosities------------------------------------

  viscqpo=visco**0.25
  viscqpg=viscg**0.25

!--Form weighted combination of the 1/4 powers & its 4th power (denom of TL 4a)
!  Note that:
!  fg multiplies the oil viscosity 1/4-power term
!  fo multiplies the gas viscosity 1/4-power term

  wviscqp=fg*viscqpo+fo*viscqpg
  denom=wviscqp**4.0

!--Obtain a safe denominator inverse-------------------------------------------

  if( denom>0.0 ) then
    denominv=1.0/denom
  else
    denominv=0.0
  endif

!--Form mixed viscosity--(TL equation 4a)--------------------------------------

  viscm=visco*viscg*denominv

!--Form omega-weighted replacement viscosities---------------------------------

! 1-omega and omega power contributions

  viscoimw=visco**tlomegac
  viscgimw=viscg**tlomegac
  viscmw  =viscm**val_tl_omega

! Combine to get final value (TL 3a and 3b)

  viscotl=viscoimw*viscmw
  viscgtl=viscgimw*viscmw

end subroutine vToddLongstaffViscosity

subroutine vToddLongstaffDensity( fo,fg,visco,viscg,viscotl,viscgtl &
                                 ,deno,deng,denotl,dengtl)

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

  implicit none

  PetscReal,intent(in) :: fo,fg,visco,viscg,viscotl,viscgtl,deno,deng
  PetscReal,intent(out):: denotl,dengtl
  PetscReal            :: tlomegac,viscginv,viscotlinv,viscgtlinv &
                         ,m,mooe,moge,mqp,mooeqp,mogeqp,mqpm,mqpm1,mqpm1inv &
                         ,fo_oe,fo_ge,fg_oe,fg_ge &
                         ,denm

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
  if( abs(m-1.0)>1.0E-6 .and. viscotl>0.0 .and. viscgtl>0.0 ) then

! Case of non-unit mobility ratio

! Form (visco(unmixed))/viscptl), p=o,g, as used in eqn. 8b and 8a in TL

    mooe=visco*viscotlinv
    moge=visco*viscgtlinv

    mqp   =m   **0.25
    mooeqp=mooe**0.25
    mogeqp=moge**0.25

    mqpm1=mqp-1.0

    mqpm1inv=0.0;
    if( abs(mqpm1)>0.0 ) then
      mqpm1inv=1.0/mqpm1
    endif

!  Form effective fractional saturations, eqn. 8b and 8a in TL

    fo_oe=(mqp-mooeqp)*mqpm1inv
    fo_ge=(mqp-mogeqp)*mqpm1inv

!  Complements of effective oil saturations, as used in eqn. 9a and 9b

    fg_oe=1.0-fo_oe
    fg_ge=1.0-fo_ge

!  Set up mixed densities, eqn. 9b and 9a

    denotl=deno*fo_oe+deng*fg_oe
    dengtl=deno*fo_ge+deng*fg_ge

  else

! Case of unit mobility ratio, eqn. 10a and 10b

    denm=fo*deno+fg*deng
    denotl=tlomegac*deno+val_tl_omega*denm
    dengtl=tlomegac*deng+val_tl_omega*denm

  endif

end subroutine vToddLongstaffDensity

! ************************************************************************** !

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

end subroutine TOWGImsTLAuxVarPerturb

! DKP Set up auxillary variables for perturbed black oil system----------------

subroutine TOWGBlackOilAuxVarPerturb(auxvar,global_auxvar, &
                                     material_auxvar, &
                                     characteristic_curves,natural_id, &
                                     option)
!------------------------------------------------------------------------------
! Used in TOWG_BLACK_OIL, auxillary variables for perturbed system
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Oct 2017
!------------------------------------------------------------------------------

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

end subroutine TOWGBlackOilAuxVarPerturb

! DKP--------------------------------------------------------------------------

subroutine TOWGImsAuxVarComputeSetup()

  implicit none

  TOWGAuxVarCompute => TOWGImsAuxVarCompute
  TOWGAuxVarPerturb => TOWGImsTLAuxVarPerturb

end subroutine TOWGImsAuxVarComputeSetup

! DKP Routine to set up black oil values and perturbed values------------------

subroutine TOWGBlackOilAuxVarComputeSetup()

  implicit none

  TOWGAuxVarCompute => TOWGBlackOilAuxVarCompute
  TOWGAuxVarPerturb => TOWGBlackOilAuxVarPerturb

end subroutine TOWGBlackOilAuxVarComputeSetup

! DKP--------------------------------------------------------------------------

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

! ************************************************************************** !

end module PM_TOWG_Aux_module
