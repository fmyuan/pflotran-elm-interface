module Hydrate_module

!
! Author: Michael Nole
! Date: 01/02/19
!
! MODULE DESCRIPTION:
! ***************************************************************************
! This module extends general mode to account for gas hydrate
! formation and dissociation
! ***************************************************************************

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Data_Mediator_Vec_class
  use Region_module

  use General_Aux_module
  use Global_Aux_module

  implicit none

  private

  PetscInt, parameter :: H_STATE = 5
  PetscInt, parameter :: ICE_STATE = 4
  PetscInt, parameter :: GA_STATE = 3
  PetscInt, parameter :: HG_STATE = 6
  PetscInt, parameter :: HA_STATE = 7
  PetscInt, parameter :: HI_STATE = 8
  PetscInt, parameter :: GI_STATE = 9
  PetscInt, parameter :: AI_STATE = 10
  PetscInt, parameter :: HGA_STATE = 11
  PetscInt, parameter :: HAI_STATE = 12
  PetscInt, parameter :: HGI_STATE = 13
  PetscInt, parameter :: GAI_STATE = 14
  PetscInt, parameter :: QUAD_STATE = 15
  
  PetscInt, parameter :: lid = 1
  PetscInt, parameter :: gid = 2
  PetscInt, parameter :: hid = 3
  PetscInt, parameter :: iid = 4

  !Structure 1 methane hydrate:
  PetscReal, parameter :: HYDRATE_DENSITY_KG = 910 !kg/m^3
  PetscReal, parameter :: HYDRATE_MW = 124.13 !957.04 !g/mol
  PetscReal, parameter :: HYDRATE_DENSITY = 51.79 !0.95 !mol/L
  PetscReal, parameter :: MOL_RATIO_METH = 0.1481d0
  PetscReal, parameter :: MOL_RATIO_H20 = 0.8519d0
  PetscReal, parameter :: MASS_RATIO_METH = 0.1341d0
  PetscReal, parameter :: MASS_RATIO_H20 = 0.8659d0
  
  !Ice: 
  PetscReal, parameter :: ICE_DENSITY_KG = 920 !kg/m^3
  PetscReal, parameter :: ICE_DENSITY = 50.86 !mol/L


  PetscReal, parameter :: lambda_hyd = 0.49 !W/m-K
  
!  type, public :: methanogenesis_type
!    character(len=MAXWORDLENGTH) :: source_name
!    PetscReal, parameter :: alpha 
!    PetscReal, parameter :: k_alpha
!    PetscReal, parameter :: lambda
!    PetscReal, parameter :: omega
!    PetscReal, parameter :: z_smt
!    type(methanogenesis_type), pointer :: next
!  end type methanogenesis_type
  
!  type, public :: methanogenesis_mediator_type
!    type(methanogenesis_type), pointer :: methanogenesis_list
!    class(data_mediator_vec_type), pointer :: data_mediator
!    PetscInt :: total_num_cells
!  end type criticality_mediator_type

  public :: HydrateSetFlowMode, &
            HydrateUpdateState, &
            HydrateAuxVarCompute, &
            HydrateAccumulation, &
            HydrateAuxVarPerturb, &
            HydratePE

contains

subroutine HydrateSetFlowMode(option)
!
! Sets the flow mode for equilibrium hydrate formation and dissociation
!
! Author: Michael Nole
! Date: 01/02/19
!

  use Option_module

  implicit none

  type(option_type) :: option

  option%iflowmode = G_MODE
  option%nphase = 4
  option%liquid_phase = 1  ! liquid_pressure
  option%gas_phase = 2     ! gas_pressure
  option%hydrate_phase = 3
  option%ice_phase = 4

  option%air_pressure_id = 3
  option%capillary_pressure_id = 4
  option%vapor_pressure_id = 5
  option%saturation_pressure_id = 6

  option%water_id = 1
  option%air_id = 2
  option%energy_id = 3

  option%nflowdof = 3
  option%nflowspec = 2
  option%use_isothermal = PETSC_FALSE

  option%hydrate_flag = PETSC_TRUE

end subroutine HydrateSetFlowMode

subroutine HydrateUpdateState(x,gen_auxvar,global_auxvar, material_auxvar, &
                              characteristic_curves,natural_id,option)
  ! 
  ! Decides on state changes and adds epsilons to new primary variables 
  ! accordingly. Primary variables for each phase state modeled 
  ! roughly after Sun, 2005
  ! 
  ! Author: Michael Nole
  ! Date: 01/28/18
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use General_Aux_module

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  class(characteristic_curves_type) :: characteristic_curves
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscReal, parameter :: epsilon = 0.d0
  PetscReal, parameter :: TQD = 0.d0 !degrees C 
  PetscReal :: liq_epsilon, gas_epsilon, hyd_epsilon, two_phase_epsilon
  PetscReal :: ga_epsilon, ha_epsilon
  PetscReal :: x(option%nflowdof)
  PetscReal :: PE_hyd, K_H, Tf_ice, dTf
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, hid, iid, acid, wid
  PetscBool :: istatechng
  PetscErrorCode :: ierr

  if (general_immiscible .or. global_auxvar%istatechng) return

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = 3
  iid = 4

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id

  gen_auxvar%istate_store(PREV_IT) = global_auxvar%istate
  istatechng = PETSC_FALSE

  gas_epsilon = 0.d0
  liq_epsilon = 0.d0
  hyd_epsilon = 0.d0
  two_phase_epsilon = 0.d0

  !man: need to implement ice once hydrate works
  !man: right now comparing hydrate equilib pressure to gas
  !pressure (assuming low water solubility in methane). 
  !Ideally would compare to partial pressure of methane.

  if (global_auxvar%hstate == ZERO_INTEGER .and. gen_auxvar%sat(gid) &
       < 0.d0) then
    global_auxvar%hstate = HA_STATE
    gen_auxvar%sat(hid) = 0.d0 !-1.d0 * gen_auxvar%sat(gid)
    gen_auxvar%sat(gid) = 0.d0
  endif

  if (global_auxvar%hstate == ZERO_INTEGER) global_auxvar% &
                              hstate = global_auxvar%istate 
  gen_auxvar%hstate_store(PREV_IT) = global_auxvar%hstate

  call HydratePE(gen_auxvar%temp,PE_hyd)
 
  call GibbsThomsonFreezing(1.d0-gen_auxvar%sat(iid),6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,option) 

  Tf_ice = TQD + dTf
  !Update State
  
  select case(global_auxvar%hstate)
    case(LIQUID_STATE)
      if (gen_auxvar%temp > Tf_ice) then
        if (gen_auxvar%pres(apid) >= gen_auxvar% &
             pres(lid)*(1.d0-window_epsilon)) then
          if (gen_auxvar%pres(apid) >= PE_hyd) then
            istatechng = PETSC_TRUE
            global_auxvar%hstate = HA_STATE
            global_auxvar%istate = TWO_PHASE_STATE
          else
            istatechng = PETSC_TRUE
            global_auxvar%hstate = GA_STATE
            global_auxvar%istate = TWO_PHASE_STATE
            liq_epsilon = option%phase_chng_epsilon
          endif
        else
          istatechng = PETSC_FALSE
        endif
      elseif (gen_auxvar%pres(apid) >= gen_auxvar%pres(lid)) then
        if (gen_auxvar%pres(apid) < PE_hyd) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GAI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          liq_epsilon = option%phase_chng_epsilon
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = QUAD_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          liq_epsilon = option%phase_chng_epsilon
        endif
      else
        istatechng = PETSC_TRUE
        global_auxvar%hstate = AI_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      endif

    case(GAS_STATE)
      if (gen_auxvar%pres(vpid) >= gen_auxvar%pres(spid)* &
         (1.d0-window_epsilon)) then
        if (gen_auxvar%pres(apid) < PE_hyd) then
          if (gen_auxvar%temp > Tf_ice) then
            istatechng = PETSC_TRUE
            global_auxvar%hstate = GA_STATE
            global_auxvar%istate = TWO_PHASE_STATE
            gas_epsilon = option%phase_chng_epsilon
          elseif (gen_auxvar%temp == Tf_ice) then
            istatechng = PETSC_TRUE             
            global_auxvar%hstate = GAI_STATE
            global_auxvar%istate =  TWO_PHASE_STATE
            gas_epsilon = option%phase_chng_epsilon           
          else
            istatechng = PETSC_TRUE
            global_auxvar%hstate = GI_STATE
            global_auxvar%istate = TWO_PHASE_STATE
            gas_epsilon = option%phase_chng_epsilon
          endif
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HG_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          gas_epsilon = option%phase_chng_epsilon
        endif
      else
        istatechng = PETSC_FALSE
      endif

    case(H_STATE)
      if (gen_auxvar%pres(apid) < PE_hyd) then
        if (gen_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HGA_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%temp == Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = QUAD_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HGI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        endif
      else
          istatechng = PETSC_FALSE
      endif

    case(ICE_STATE)
      if (gen_auxvar%temp > Tf_ice) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = QUAD_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      else
        istatechng = PETSC_FALSE
      endif 

    case(GA_STATE)
      if (gen_auxvar%pres(apid) < PE_hyd) then
        if (gen_auxvar%temp > Tf_ice) then
          if (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(lid) > 0.d0) then
            istatechng = PETSC_FALSE
          elseif (gen_auxvar%sat(gid) <= 0.d0) then
            istatechng = PETSC_TRUE
            global_auxvar%hstate = LIQUID_STATE
            global_auxvar%istate = LIQUID_STATE
            two_phase_epsilon = option%phase_chng_epsilon !*1.d-5
          elseif (gen_auxvar%sat(gid) >= 1.d0) then
            istatechng = PETSC_TRUE
            global_auxvar%hstate = GAS_STATE
            global_auxvar%istate = GAS_STATE
            two_phase_epsilon = option%phase_chng_epsilon !*1.d-5
          endif
        else 
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GAI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          two_phase_epsilon = option%phase_chng_epsilon
        endif
      else
        if (gen_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HGA_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          two_phase_epsilon = option%phase_chng_epsilon
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = QUAD_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          two_phase_epsilon = option%phase_chng_epsilon
        endif
      endif

    case(HG_STATE)
      if (gen_auxvar%pres(apid) > PE_hyd) then
        if (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (gen_auxvar%sat(hid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = H_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          two_phase_epsilon = option%phase_chng_epsilon
        elseif (gen_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GAS_STATE
          global_auxvar%istate = GAS_STATE
          two_phase_epsilon = option%phase_chng_epsilon
        endif
      else
        if (gen_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HGA_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          two_phase_epsilon = option%phase_chng_epsilon
        elseif (gen_auxvar%temp == Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = QUAD_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          two_phase_epsilon = option%phase_chng_epsilon
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HGI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          two_phase_epsilon = option%phase_chng_epsilon
        endif
      endif

    case(HA_STATE)
      if (gen_auxvar%pres(apid) > PE_hyd .and. gen_auxvar%temp > Tf_ice) then
        if (gen_auxvar%sat(hid) >0.d0 .and. gen_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (gen_auxvar%sat(hid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = H_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = LIQUID_STATE
          global_auxvar%istate = LIQUID_STATE
        endif
      elseif (gen_auxvar%temp > Tf_ice) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = HGA_STATE
        global_auxvar%istate = TWO_PHASE_STATE
        ha_epsilon = option%phase_chng_epsilon
      elseif (gen_auxvar%pres(apid) > PE_hyd) then
       istatechng = PETSC_TRUE
       global_auxvar%hstate = HAI_STATE
       global_auxvar%istate = TWO_PHASE_STATE
      else
       istatechng = PETSC_TRUE
       global_auxvar%hstate = QUAD_STATE
       global_auxvar%istate = TWO_PHASE_STATE
      endif

    case(HI_STATE)
      if (gen_auxvar%pres(apid) > PE_hyd) then
        if (gen_auxvar%temp < Tf_ice) then
          istatechng = PETSC_FALSE
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HAI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        endif
      else
        if (gen_auxvar%temp < Tf_ice) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HGI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = QUAD_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        endif
      endif

    case(GI_STATE)
      if (gen_auxvar%temp < Tf_ice .and. gen_auxvar%pres(apid) < PE_hyd) then
        if (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (gen_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GAS_STATE
          global_auxvar%istate = GAS_STATE
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = ICE_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        endif
      elseif (gen_auxvar%temp < Tf_ice) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = HGI_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%pres(apid) < PE_hyd) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = GAI_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      else
        istatechng = PETSC_TRUE
        global_auxvar%hstate = QUAD_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      endif

    case(AI_STATE)
      if (gen_auxvar%pres(apid) >= gen_auxvar% &
             pres(lid)*(1.d0-window_epsilon)) then 
        if (gen_auxvar%pres(apid) < PE_hyd) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GAI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HAI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        endif
      else
        if (gen_auxvar%sat(lid) > 0.d0 .and. gen_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (gen_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = LIQUID_STATE
          global_auxvar%istate = LIQUID_STATE
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = ICE_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        endif
      endif

    case(HGA_STATE)
      if (gen_auxvar%temp > Tf_ice) then
        if (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(gid) > 0.d0 &
            .and. gen_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(lid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HA_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(gid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HG_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(lid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GA_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = LIQUID_STATE
          global_auxvar%istate = LIQUID_STATE
        elseif (gen_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GAS_STATE
          global_auxvar%istate = GAS_STATE
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = H_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        endif
      else
        istatechng = PETSC_TRUE
        global_auxvar%hstate = QUAD_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      endif

    case(HAI_STATE)
      if (gen_auxvar%pres(apid) > PE_hyd) then
        if (gen_auxvar%sat(lid) > 0.d0 .and. gen_auxvar%sat(hid) > 0.d0 &
            .and. gen_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (gen_auxvar%sat(lid) > 0.d0 .and. gen_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = AI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(lid) > 0.d0 .and. gen_auxvar%sat(hid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HA_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = LIQUID_STATE
          global_auxvar%istate = LIQUID_STATE
        elseif (gen_auxvar%sat(hid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = H_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = ICE_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        endif
      else
        istatechng = PETSC_TRUE
        global_auxvar%hstate = QUAD_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      endif 
    case(HGI_STATE)
      if (gen_auxvar%temp < Tf_ice) then
        if (gen_auxvar%sat(iid) > 0.d0 .and. gen_auxvar%sat(hid) > 0.d0 &
            .and. gen_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(hid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HG_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = ICE_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GAS_STATE
          global_auxvar%istate = GAS_STATE
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = H_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        endif
      else
        istatechng = PETSC_TRUE
        global_auxvar%hstate = QUAD_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      endif
      
    case(GAI_STATE)
      if (gen_auxvar%pres(apid) < PE_hyd) then
        if (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(lid) > 0.d0 &
            .and. gen_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(lid) &
                > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GA_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(iid) &
                > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(lid) > 0.d0 .and. gen_auxvar%sat(iid) &
                > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = AI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        elseif (gen_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GAS_STATE
          global_auxvar%istate = GAS_STATE
        elseif (gen_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = LIQUID_STATE
          global_auxvar%istate = LIQUID_STATE
        else
          istatechng = PETSC_TRUE
          global_auxvar%hstate = ICE_STATE
          global_auxvar%istate = TWO_PHASE_STATE
        endif
      else
        istatechng = PETSC_TRUE
        global_auxvar%hstate = QUAD_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      endif
      
    case(QUAD_STATE)
!
      if (gen_auxvar%sat(lid) > 0.d0 .and. gen_auxvar%sat(gid) > 0.d0 &
          .and. gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(iid) &
          > 0.d0) then
        istatechng = PETSC_FALSE
      elseif (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(gid) &
               > 0.d0 .and. gen_auxvar%sat(lid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = HGA_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(lid) &
              > 0.d0 .and. gen_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = HAI_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(lid) &
              > 0.d0 .and. gen_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = GAI_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(gid) &
              > 0.d0 .and. gen_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = HGI_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(lid) > 0.d0 .and. gen_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = AI_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(lid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = GA_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(lid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = HA_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = HI_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = GI_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(hid) > 0.d0 .and. gen_auxvar%sat(hid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = HG_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(lid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = LIQUID_STATE
        global_auxvar%istate = LIQUID_STATE
      elseif (gen_auxvar%sat(hid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = H_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      elseif (gen_auxvar%sat(gid) > 0.d0) then
        istatechng = PETSC_TRUE
        global_auxvar%hstate = GAS_STATE
        global_auxvar%istate = GAS_STATE
      else
        istatechng = PETSC_TRUE
        global_auxvar%hstate = ICE_STATE
        global_auxvar%istate = TWO_PHASE_STATE
      endif
  end select


  !Update primary variables

  if (istatechng) then

    if (option%restrict_state_chng) global_auxvar%istatechng = PETSC_TRUE

!    MAN: need to figure out the epsilons
!
    select case(global_auxvar%hstate)

      case(LIQUID_STATE)
!     ********* Aqueous State (A) ********************************
!     Primary variables: Pa, Xma, T
!
        x(GENERAL_LIQUID_PRESSURE_DOF) = gen_auxvar%pres(gid)
        x(GENERAL_LIQUID_STATE_X_MOLE_DOF) = gen_auxvar%xmol(acid,lid)
        x(GENERAL_ENERGY_DOF) = x(GENERAL_ENERGY_DOF)

      case(GAS_STATE)
!     ********* Gas State (G) ********************************
!     Primary variables: Pg, Pa, T
!
        x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = gen_auxvar%pres(apid)
        x(GENERAL_ENERGY_DOF) = x(GENERAL_ENERGY_DOF)
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid)
      
      case(H_STATE)
!     ********* Hydrate State (H) ********************************
!     Primary variables: Pg, Xmh, T
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid)
        x(GENERAL_GAS_SATURATION_DOF) = MOL_RATIO_METH
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp
      
      case(ICE_STATE)
!     ********* Ice State (I) ********************************
!     Primary variables: Pg, Xmi, T
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid)
        x(GENERAL_GAS_SATURATION_DOF) = 0.d0
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp

      case(GA_STATE)
!     ********* Gas & Aqueous State (GA) ********************************
!     Primary variables: Pg, Sg, T
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%sat(gid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp
      
      case(HG_STATE)
!     ********* Hydrate & Gas State (HG) ********************************
!     Primary variables: Pg, Sg, T
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%sat(gid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp

      case(HA_STATE)
!     ********* Hydrate & Aqueous State (HA) ********************************
!     Primary variables: Pg, Sh, T
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%sat(hid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp
      
      case(HI_STATE)
!     ********* Hydrate & Ice State (HI) ********************************
!     Primary variables: Pg, Sh, T
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%sat(hid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp
        
      case(GI_STATE)
!     ********* Gas & Ice State (GI) ********************************
!     Primary variables: Pg, Si, T
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%sat(iid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp
      
      case(AI_STATE)
!     ********* Aqueous & Ice State (AI) ********************************
!     Primary variables: Pg, Xma, Sl
!
        x(GENERAL_LIQUID_PRESSURE_DOF) = gen_auxvar%pres(lid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%xmol(acid,lid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%sat(lid)
        
      case(HGA_STATE)
!     ********* Hydrate, Gas, & Aqueous State (HGA) **************************
!     Primary variables: Sl, Sh, T
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%sat(lid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%sat(hid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp
      
      case(HAI_STATE)
!     ********* Hydrate, Aqueous, & Ice State (HAI) **************************
!     Primary variables: Pg, Sl, Si
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%sat(lid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%sat(iid)
      
      case(HGI_STATE)
!     ********* Hydrate, Gas, & Ice State (HGI) ******************************
!     Primary variables: Sh, Si, T
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%sat(hid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%sat(iid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%temp
      
      case(GAI_STATE)
!     ********* Gas, Aqueous, & Ice State (GAI) ******************************
!     Primary variables: Pg, Sl, Si
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%pres(gid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%sat(lid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%sat(iid)
      
      case(QUAD_STATE)
!     ********* Quadruple Point (HGAI) ********************************
!     Primary variables: Sg, Sl, Si
!
        x(GENERAL_GAS_PRESSURE_DOF) = gen_auxvar%sat(gid)
        x(GENERAL_GAS_SATURATION_DOF) = gen_auxvar%sat(lid)
        x(GENERAL_ENERGY_DOF) = gen_auxvar%sat(iid)
      
      case default
        write(option%io_buffer,*) global_auxvar%hstate
        option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
          ') not recognized in HydrateUpdateState.'
        call printErrMsgByRank(option)

    end select


    call HydrateAuxVarCompute(x,gen_auxvar, global_auxvar,material_auxvar, &
          characteristic_curves,natural_id,option)

  endif

end subroutine HydrateUpdateState

subroutine HydrateAuxVarCompute(x,gen_auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell, with gas hydrate physics
  !
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use Creep_Closure_module
  use Fracture_module
  use WIPP_module

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(creep_closure_type), pointer :: creep_closure
  PetscInt :: natural_id

  PetscInt :: gid, lid, acid, wid, eid, hid, iid
  PetscReal :: cell_pressure, water_vapor_pressure
  PetscReal :: den_water_vapor, den_kg_water_vapor
  PetscReal :: u_water_vapor, h_water_vapor
  PetscReal :: den_air, h_air, u_air
  PetscReal :: xmol_air_in_gas, xmol_water_in_gas
  PetscReal :: krl, visl, dvis_dp, dvis_dT, dvis_dpa
  PetscReal :: dkrl_dsatl, dkrl_dsatg
  PetscReal :: dkrg_dsatl, dkrg_dsatg
  PetscReal :: krg, visg
  PetscReal :: K_H_tilde
  PetscReal :: guess, dummy
  PetscInt :: apid, cpid, vpid, spid
  PetscReal :: NaN
  PetscReal :: creep_closure_time
  PetscReal :: xmass_air_in_gas
  PetscReal :: Ugas_J_kg, Hgas_J_kg
  PetscReal :: Uair_J_kg, Hair_J_kg
  PetscReal :: Uvapor_J_kg, Hvapor_J_kg
  PetscReal :: Hg_mixture_fractioned
  PetscReal :: U_hyd, U_ice, PE_hyd
  PetscReal :: aux(1)
  PetscReal :: hw, hw_dp, hw_dT
  PetscReal :: dpor_dp
  PetscReal :: one_over_dw
  PetscReal :: tempreal, tempreal2, tempreal3
  PetscReal :: dpair_dT, dpair_dpgas
  PetscReal :: dden_air_dT, dden_air_dpa, dden_air_dpg
  PetscReal :: du_air_dT, dh_air_dT
  PetscReal :: du_air_dpa, dh_air_dpa
  PetscReal :: dden_water_vapor_dpv, dden_water_vapor_dT
  PetscReal :: dh_water_vapor_dpv, dh_water_vapor_dT
  PetscReal :: du_water_vapor_dpv, du_water_vapor_dT
  PetscReal :: dpc_dsatl
  PetscReal :: dden_ice_dT, dden_ice_dP
  character(len=8) :: state_char
  PetscErrorCode :: ierr
  PetscReal, parameter :: TQD = 0.d0 
  PetscReal :: dTf
  PetscReal :: sigma

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id


#ifdef DEBUG_GENERAL
  ! create a NaN
  NaN = InitToNan()

  gen_auxvar%H = NaN
  gen_auxvar%U = NaN
  gen_auxvar%pres = NaN
  gen_auxvar%sat = NaN
  gen_auxvar%den = NaN
  gen_auxvar%den_kg = NaN
  gen_auxvar%xmol = NaN
  gen_auxvar%effective_porosity = NaN
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      state_char = 'L'
    case(GAS_STATE)
      state_char = 'G'
    case(TWO_PHASE_STATE)
      state_char = '2P'
  end select
#else
  !geh: do not initialize gen_auxvar%temp as the previous value is used as the
  !     initial guess for two phase.
  gen_auxvar%H = 0.d0
  gen_auxvar%U = 0.d0
  gen_auxvar%pres = 0.d0
  gen_auxvar%sat = 0.d0
  gen_auxvar%den = 0.d0
  gen_auxvar%den_kg = 0.d0
  gen_auxvar%xmol = 0.d0
  gen_auxvar%effective_porosity = 0.d0
#endif
  gen_auxvar%mobility = 0.d0

#if 0
  if (option%iflag >= GENERAL_UPDATE_FOR_ACCUM) then
    if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
        natural_id, x(1:3), trim(state_char)
    else
!      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
!        -1*natural_id, x(1:3), trim(state_char)
    endif
  endif
#endif

  if (global_auxvar%hstate == ZERO_INTEGER) global_auxvar% &
                              hstate = global_auxvar%istate

  gen_auxvar%xmol(wid,hid) = MOL_RATIO_H20
  gen_auxvar%xmol(acid,hid) = MOL_RATIO_METH

  select case(global_auxvar%hstate)
    case(LIQUID_STATE)
!     ********* Aqueous State (A) ********************************
!     Primary variables: Pa, Xma, T
!
      gen_auxvar%pres(lid) = x(GENERAL_LIQUID_PRESSURE_DOF)
      gen_auxvar%xmol(acid,lid) = x(GENERAL_LIQUID_STATE_X_MOLE_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)

      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(wid,gid) = 0.d0
      gen_auxvar%xmol(acid,gid) = 0.d0
      gen_auxvar%sat(lid) = 1.d0
      gen_auxvar%sat(gid) = 0.d0
      gen_auxvar%sat(hid) = 0.d0
      gen_auxvar%sat(iid) = 0.d0

      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                        gen_auxvar%pres(spid),ierr)
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      
      gen_auxvar%pres(spid) = 1.d-6
      
      gen_auxvar%pres(gid) = max(gen_auxvar%pres(lid),gen_auxvar%pres(spid))
      gen_auxvar%pres(apid) = K_H_tilde*gen_auxvar%xmol(acid,lid)

      if (gen_auxvar%pres(gid) <= 0.d0) then
        write(option%io_buffer,'(''Negative gas pressure at cell '', &
          & i8,'' in HydrateAuxVarCompute(LIQUID_STATE).  Attempting bailout.'')') &
          natural_id
        call printMsgByRank(option)
        gen_auxvar%pres(vpid) = 0.5d0*gen_auxvar%pres(spid)
        gen_auxvar%pres(gid) = gen_auxvar%pres(vpid) + gen_auxvar%pres(apid)
      else
        gen_auxvar%pres(vpid) = gen_auxvar%pres(lid) - gen_auxvar%pres(apid)
      endif
      gen_auxvar%pres(cpid) = 0.d0
    case (GAS_STATE)
!     ********* Gas State (G) ********************************
!     Primary variables: Pg, Pa, T
!
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)

      gen_auxvar%pres(apid) = x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)

      gen_auxvar%sat(lid) = 0.d0
      gen_auxvar%sat(gid) = 1.d0
      gen_auxvar%sat(hid) = 0.d0
      gen_auxvar%sat(iid) = 0.d0

      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / &
                                   gen_auxvar%pres(gid)
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)

      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                    gen_auxvar%pres(spid),ierr)

      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)

      gen_auxvar%pres(spid) = 1.d-6
      
      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 0.d0

      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)

      call characteristic_curves%saturation_function% &
             CapillaryPressure(gen_auxvar%sat(lid), &
                               gen_auxvar%pres(cpid),dpc_dsatl,option)
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - &
                             gen_auxvar%pres(cpid)
    case (H_STATE)
!     ********* Hydrate State (H) ********************************
!     Primary variables: Pg, Xmh, T
!
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      x(GENERAL_GAS_SATURATION_DOF) = MOL_RATIO_METH
      gen_auxvar%xmol(acid,hid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%sat(lid) = 0.d0
      gen_auxvar%sat(gid) = 0.d0
      gen_auxvar%sat(hid) = 1.d0
      gen_auxvar%sat(iid) = 0.d0
      
      call HydratePE(gen_auxvar%temp,PE_hyd)
      gen_auxvar%pres(apid) = PE_hyd
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)
      
      gen_auxvar%pres(cpid) = 0.d0
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(wid,gid) = gen_auxvar%pres(vpid) / gen_auxvar%pres(gid)
      gen_auxvar%xmol(acid,gid) = 1.d0 - gen_auxvar%xmol(wid,gid)
      
    case(ICE_STATE)
!     ********* Ice State (I) ********************************
!     Primary variables: Pg, Xmi, T
!
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      x(GENERAL_GAS_SATURATION_DOF) = 0.d0
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%sat(lid) = 0.d0
      gen_auxvar%sat(gid) = 0.d0
      gen_auxvar%sat(hid) = 0.d0
      gen_auxvar%sat(iid) = 1.d0
      
      gen_auxvar%pres(cpid) = 0.d0
      gen_auxvar%pres(apid) = 0.d0
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)

      gen_auxvar%xmol(acid,lid) = 0.d0
      gen_auxvar%xmol(wid,lid) = 1.d0
      gen_auxvar%xmol(wid,gid) = 1.d0
      gen_auxvar%xmol(acid,gid) = 0.d0
      
    case(GA_STATE)
!     ********* Gas & Aqueous State (GA) ********************************
!     Primary variables: Pg, Sg, T
!
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(gid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%sat(lid) = 1.d0 - gen_auxvar%sat(gid)
      gen_auxvar%sat(hid) = 0.d0
      gen_auxvar%sat(iid) = 0.d0

      !man
!        gen_auxvar%sat(gid) = max(0.d0,gen_auxvar%sat(gid))
!        gen_auxvar%sat(gid) = min(1.d0,gen_auxvar%sat(gid))

      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                          gen_auxvar%pres(spid),ierr)
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
!        call EOSGasHenry(gen_auxvar%temp,gen_auxvar%pres(spid),K_H_tilde)

      gen_auxvar%pres(spid) = 1.d-6

      if (general_immiscible) then
        gen_auxvar%pres(spid) = GENERAL_IMMISCIBLE_VALUE
      endif
      gen_auxvar%pres(vpid) = gen_auxvar%pres(spid)
      gen_auxvar%pres(apid) = gen_auxvar%pres(gid) - gen_auxvar%pres(vpid)


      call characteristic_curves%saturation_function% &
             CapillaryPressure(gen_auxvar%sat(lid), &
                               gen_auxvar%pres(cpid),dpc_dsatl,option)

      !man: IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(gen_auxvar%temp,sigma)
      endif
      gen_auxvar%pres(cpid) = gen_auxvar%pres(cpid)*sigma

      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      if (general_immiscible) then
        gen_auxvar%xmol(acid,lid) = GENERAL_IMMISCIBLE_VALUE
      endif

      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / gen_auxvar%pres(gid)
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)
    
    case(HG_STATE)
!     ********* Hydrate & Gas State (HG) ********************************
!     Primary variables: Pg, Sg, T
!
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(gid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%sat(lid) = 0.d0
      gen_auxvar%sat(hid) = 1.d0 - gen_auxvar%sat(gid)
      gen_auxvar%sat(iid) = 0.d0
      
      call HydratePE(gen_auxvar%temp,PE_hyd)
      gen_auxvar%pres(apid) = PE_hyd
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)
      
      gen_auxvar%pres(cpid) = 0.d0
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(wid,gid) = gen_auxvar%pres(vpid) / gen_auxvar%pres(gid)
      gen_auxvar%xmol(acid,gid) = 1.d0 - gen_auxvar%xmol(wid,gid)
      
      
      
    case(HA_STATE)
!     ********* Hydrate & Aqueous State (HA) ********************************
!     Primary variables: Pg, Sh, T
!

      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(hid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)

!       gen_auxvar%sat(hid) = max(0.d0,gen_auxvar%sat(hid))
!       gen_auxvar%sat(hid) = min(1.d0,gen_auxvar%sat(hid))

      gen_auxvar%sat(lid) = 1.d0 - gen_auxvar%sat(hid)
      gen_auxvar%sat(gid) = 0.d0
      gen_auxvar%sat(iid) = 0.d0
      
      call HydratePE(gen_auxvar%temp,PE_hyd)
      gen_auxvar%pres(apid) = PE_hyd
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                          gen_auxvar%pres(spid),ierr)

      gen_auxvar%pres(cpid) = 0.d0
      ! Setting air pressure equal to gas pressure makes forming hydrate
      ! easier
      gen_auxvar%pres(apid) = gen_auxvar%pres(gid)
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(wid,gid) = 0.d0
      gen_auxvar%xmol(acid,gid) = 0.d0


    case(HI_STATE)
!     ********* Hydrate & Ice State (HI) ********************************
!     Primary variables: Pg, Sh, T
!
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(hid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%sat(lid) = 0.d0
      gen_auxvar%sat(gid) = 0.d0
      gen_auxvar%sat(iid) = 1.d0 - gen_auxvar%sat(hid)
      
      call HydratePE(gen_auxvar%temp,PE_hyd)
      gen_auxvar%pres(apid) = PE_hyd
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      
      gen_auxvar%pres(cpid) = 0.d0
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(wid,gid) = gen_auxvar%pres(vpid) / gen_auxvar%pres(gid)
      gen_auxvar%xmol(acid,gid) = 1.d0 - gen_auxvar%xmol(wid,gid)
      

    case(GI_STATE)
!     ********* Gas & Ice State (GI) ********************************
!     Primary variables: Pg, Si, T
!
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(iid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%sat(lid) = 0.d0
      gen_auxvar%sat(gid) = 1.d0 - gen_auxvar%sat(iid)
      gen_auxvar%sat(hid) = 0.d0
      
      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                          gen_auxvar%pres(spid),ierr)
                              
      gen_auxvar%pres(spid) = 1.d-6                        
                              
      gen_auxvar%pres(vpid) = gen_auxvar%pres(spid)
      gen_auxvar%pres(apid) = gen_auxvar%pres(gid) - gen_auxvar%pres(vpid)
      
      gen_auxvar%pres(cpid) = 0.d0
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)

      gen_auxvar%xmol(acid,lid) = 0.d0
      gen_auxvar%xmol(wid,lid) = 1.d0
      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / gen_auxvar%pres(gid)
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)

    case(AI_STATE)
!     ********* Aqueous & Ice State (AI) ********************************
!     Primary variables: Pl, Xma, Sl
!
      gen_auxvar%pres(lid) = x(GENERAL_LIQUID_PRESSURE_DOF)
      gen_auxvar%xmol(acid,lid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%sat(lid) = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%sat(lid) = max(0.d0,min(1.d0,gen_auxvar%sat(lid)))

      gen_auxvar%sat(gid) = 0.d0
      gen_auxvar%sat(hid) = 0.d0
      gen_auxvar%sat(iid) = 1.d0 - gen_auxvar%sat(lid)
      
      call GibbsThomsonFreezing(gen_auxvar%sat(lid),6017.1d0,ICE_DENSITY,&
                                TQD, dTf,characteristic_curves,option)      
      gen_auxvar%temp = TQD+dTf
      
      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                          gen_auxvar%pres(spid),ierr)
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)

      gen_auxvar%pres(spid) = 1.d-6
      gen_auxvar%pres(gid) = max(gen_auxvar%pres(lid),gen_auxvar%pres(spid))
      gen_auxvar%pres(cpid) = 0.d0
      gen_auxvar%pres(apid) = K_H_tilde*gen_auxvar%xmol(acid,lid)
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(acid,gid) = 0.d0
      gen_auxvar%xmol(wid,gid) = 0.d0
      

    case(HGA_STATE)
!     ********* Hydrate, Gas, & Aqueous State (HGA) **************************
!     Primary variables: Sl, Sh, T
!
      gen_auxvar%sat(lid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(hid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%sat(gid) = 1.d0 - gen_auxvar%sat(lid) - gen_auxvar%sat(hid)
      gen_auxvar%sat(gid) = max(gen_auxvar%sat(gid), 0.d0)
      gen_auxvar%sat(iid) = 0.d0
      
      
      call HydratePE(gen_auxvar%temp,PE_hyd)
      
      gen_auxvar%pres(apid) = PE_hyd

      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                    gen_auxvar%pres(spid),ierr)
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      
      gen_auxvar%pres(spid) = 1.d-6
      
      gen_auxvar%pres(vpid) = gen_auxvar%pres(spid)
      gen_auxvar%pres(gid) = gen_auxvar%pres(apid) + gen_auxvar%pres(vpid)

      call characteristic_curves%saturation_function% &
             CapillaryPressure(gen_auxvar%sat(lid), &
                               gen_auxvar%pres(cpid),dpc_dsatl,option)

      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(gen_auxvar%temp,sigma)
      endif
      gen_auxvar%pres(cpid) = gen_auxvar%pres(cpid)*sigma

      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / gen_auxvar%pres(gid)
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)


    case(HAI_STATE)
!     ********* Hydrate, Aqueous, & Ice State (HAI) **************************
!     Primary variables: Pg, Sl, Si
!
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(lid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%sat(iid) = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%sat(gid) = 0.d0
      gen_auxvar%sat(hid) = 1.d0 - gen_auxvar%sat(lid) - gen_auxvar%sat(iid)
     
      call GibbsThomsonFreezing(1.d0-gen_auxvar%sat(iid), 6017.1d0, &
                   ICE_DENSITY, TQD, dTf,characteristic_curves,option)
      gen_auxvar%temp = TQD+dTf 
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      
      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                          gen_auxvar%pres(spid),ierr)

      gen_auxvar%pres(cpid) = 0.d0
      gen_auxvar%pres(spid) = 1.d-6
      gen_auxvar%pres(apid) = gen_auxvar%pres(gid)
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)
      gen_auxvar%pres(vpid) = gen_auxvar%pres(gid) - gen_auxvar%pres(apid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(wid,gid) = 0.d0
      gen_auxvar%xmol(acid,gid) = 0.d0

    case(HGI_STATE)
!     ********* Hydrate, Gas, & Ice State (HGI) ******************************
!     Primary variables: Sh, Si, T
!
      gen_auxvar%sat(hid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(iid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%sat(lid) = 0.d0
      gen_auxvar%sat(gid) = 1.d0 - gen_auxvar%sat(hid) - gen_auxvar%sat(iid)
      
      call HydratePE(gen_auxvar%temp,PE_hyd)
      gen_auxvar%pres(apid) = PE_hyd
      
      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                    gen_auxvar%pres(spid),ierr)
      
      gen_auxvar%pres(spid) = 1.d-6
      
      gen_auxvar%pres(vpid) = gen_auxvar%pres(spid)
      gen_auxvar%pres(gid) = gen_auxvar%pres(apid) + gen_auxvar%pres(vpid)
      
      gen_auxvar%xmol(acid,lid) = 0.d0
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(wid,gid) = gen_auxvar%pres(vpid) / gen_auxvar%pres(gid)
      gen_auxvar%xmol(acid,gid) = 1.d0 - gen_auxvar%xmol(wid,gid)
    
    case(GAI_STATE)
!     ********* Gas, Aqueous, & Ice State (GAI) ******************************
!     Primary variables: Pg, Sl, Si
!
      gen_auxvar%pres(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(lid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%sat(iid) = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%temp = TQD
      
      gen_auxvar%sat(gid) = 1.d0 - gen_auxvar%sat(lid) - gen_auxvar%sat(iid)
      gen_auxvar%sat(hid) = 0.d0

      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                    gen_auxvar%pres(spid),ierr)
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      
      gen_auxvar%pres(spid) = 1.d-6
      
      gen_auxvar%pres(vpid) = gen_auxvar%pres(spid)
      gen_auxvar%pres(apid) = gen_auxvar%pres(gid) - gen_auxvar%pres(vpid)

      call characteristic_curves%saturation_function% &
             CapillaryPressure(gen_auxvar%sat(lid), &
                               gen_auxvar%pres(cpid),dpc_dsatl,option)

      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(gen_auxvar%temp,sigma)
      endif
      gen_auxvar%pres(cpid) = gen_auxvar%pres(cpid)*sigma

      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / gen_auxvar%pres(gid)
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)
      
    case(QUAD_STATE)
!     ********* Quadruple Point (HGAI) ********************************
!     Primary variables: Sg, Sl, Si
!
      gen_auxvar%sat(gid) = x(GENERAL_GAS_PRESSURE_DOF)
      gen_auxvar%sat(lid) = x(GENERAL_GAS_SATURATION_DOF)
      gen_auxvar%sat(iid) = x(GENERAL_ENERGY_DOF)
      
      gen_auxvar%temp = TQD
      
      call HydratePE(gen_auxvar%temp,PE_hyd)
      gen_auxvar%pres(apid) = PE_hyd
      
      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                    gen_auxvar%pres(spid),ierr)
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      
      gen_auxvar%pres(spid) = 1.d-6
      
      gen_auxvar%pres(vpid) = gen_auxvar%pres(spid)
      gen_auxvar%pres(gid) = gen_auxvar%pres(apid) + gen_auxvar%pres(vpid)

      call characteristic_curves%saturation_function% &
             CapillaryPressure(gen_auxvar%sat(lid), &
                               gen_auxvar%pres(cpid),dpc_dsatl,option)

      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(gen_auxvar%temp,sigma)
      endif
      gen_auxvar%pres(cpid) = gen_auxvar%pres(cpid)*sigma

      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)

      gen_auxvar%xmol(acid,lid) = gen_auxvar%pres(apid) / K_H_tilde
      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(acid,gid) = gen_auxvar%pres(apid) / gen_auxvar%pres(gid)
      gen_auxvar%xmol(wid,gid) = 1.d0 - gen_auxvar%xmol(acid,gid)
      
      
    case default
      write(option%io_buffer,*) global_auxvar%hstate
      option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
        ') not recognized in HydrateAuxVarCompute.'
      call printErrMsgByRank(option)

  end select

  cell_pressure = max(gen_auxvar%pres(lid),gen_auxvar%pres(gid), &
                      gen_auxvar%pres(spid))

  ! calculate effective porosity as a function of pressure
  if (option%iflag /= GENERAL_UPDATE_FOR_BOUNDARY) then
    dpor_dp = 0.d0
    gen_auxvar%effective_porosity = material_auxvar%porosity_base
#if 0
!geh this code is no longer valid
    if (associated(material_auxvar%fracture) .and. & 
      material_auxvar%fracture%setup) then
      ! The initiating pressure and maximum pressure must be calculated
      ! before fracture function applies - Heeho
      call FractureInitialSetup(material_auxvar,cell_pressure)
    endif
    if (soil_compressibility_index > 0 .and. &
      material_auxvar%setup_reference_pressure) then
      call MaterialReferencePressureSetup(material_auxvar,cell_pressure)
    endif
#endif
    ! creep_closure, fracture, and soil_compressibility are mutually exclusive
    if (option%flow%creep_closure_on) then
      creep_closure => wipp%creep_closure_tables_array( &
                         material_auxvar%creep_closure_id )%ptr
      if ( associated(creep_closure) ) then
        ! option%time here is the t time, not t + dt time.
        creep_closure_time = option%time
        if (option%iflag /= GENERAL_UPDATE_FOR_FIXED_ACCUM) then
          creep_closure_time = creep_closure_time + option%flow_dt
        endif
        
        gen_auxvar%effective_porosity = &
          creep_closure%Evaluate(creep_closure_time,cell_pressure)
      else if (associated(material_auxvar%fracture)) then
          call FracturePoroEvaluate(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dpor_dp)
      else if (soil_compressibility_index > 0) then
          call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dpor_dp)
      endif
    else if (associated(material_auxvar%fracture)) then
      call FracturePoroEvaluate(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dpor_dp)
    else if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dpor_dp)
    endif
    if (option%iflag /= GENERAL_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = gen_auxvar%effective_porosity
    endif
  endif
  if (associated(gen_auxvar%d)) then
    gen_auxvar%d%por_p = dpor_dp
  endif                      
                      
  !MAN: Need to add permeability change as function of hydrate saturation.

  ! ALWAYS UPDATE THERMODYNAMIC PROPERTIES

  ! Liquid phase thermodynamic properties
  ! must use cell_pressure as the pressure, not %pres(lid)
  if (.not.option%flow%density_depends_on_salinity) then
    if (associated(gen_auxvar%d)) then
      call EOSWaterDensity(gen_auxvar%temp,cell_pressure, &
                           gen_auxvar%den_kg(lid),gen_auxvar%den(lid), &
                           gen_auxvar%d%denl_pl,gen_auxvar%d%denl_T,ierr)
    else
      call EOSWaterDensity(gen_auxvar%temp,cell_pressure, &
                           gen_auxvar%den_kg(lid),gen_auxvar%den(lid),ierr)
    endif
  else
    aux(1) = global_auxvar%m_nacl(1)
    call EOSWaterDensityExt(gen_auxvar%temp,cell_pressure,aux, &
                              gen_auxvar%den_kg(lid),gen_auxvar%den(lid),ierr)
  endif

  call EOSWaterEnthalpy(gen_auxvar%temp,cell_pressure,hw,ierr)

  gen_auxvar%H(lid) = hw * 1.d-6 ! J/kmol -> MJ/kmol
  ! MJ/kmol comp
  gen_auxvar%U(lid) = gen_auxvar%H(lid) - &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        (cell_pressure / gen_auxvar%den(lid) * &
                        1.d-6)
  if (global_auxvar%hstate .ne. LIQUID_STATE) then 
    if (global_auxvar%hstate == GAS_STATE) then
      water_vapor_pressure = gen_auxvar%pres(vpid)
    else
      water_vapor_pressure = gen_auxvar%pres(spid)
    endif
    call EOSGasDensityEnergy(gen_auxvar%temp,gen_auxvar%pres(apid),den_air, &
                               h_air,u_air,ierr)
    h_air = h_air * 1.d-6 ! J/kmol -> MJ/kmol
    u_air = u_air * 1.d-6 ! J/kmol -> MJ/kmol

    call EOSWaterSteamDensityEnthalpy(gen_auxvar%temp,water_vapor_pressure, &
                                        den_kg_water_vapor,den_water_vapor, &
                                        h_water_vapor,ierr)
    u_water_vapor = h_water_vapor - &
                    ! Pa / kmol/m^3 = J/kmol
                    water_vapor_pressure / den_water_vapor
    h_water_vapor = h_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    u_water_vapor = u_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    gen_auxvar%den(gid) = den_water_vapor + den_air
    gen_auxvar%den_kg(gid) = den_kg_water_vapor + den_air*fmw_comp(gid)

    xmol_air_in_gas = gen_auxvar%xmol(acid,gid)
    xmol_water_in_gas = gen_auxvar%xmol(wid,gid)

#ifdef DEBUG_GENERAL
    xmass_air_in_gas = xmol_air_in_gas*fmw_comp(gid) / &
                       (xmol_water_in_gas*FMWH2O + &
                        xmol_air_in_gas*fmw_comp(gid))
    Hair_J_kg = h_air*1.d6/fmw_comp(gid)
    Uair_J_kg = u_air*1.d6/fmw_comp(gid)
    Hvapor_J_kg = h_water_vapor*1.d6/FMWH2O
    Uvapor_J_kg = u_water_vapor*1.d6/FMWH2O
    Ugas_J_kg = xmass_air_in_gas*Uair_J_kg + &
                (1.d0-xmass_air_in_gas)*Uvapor_J_kg
    Hgas_J_kg = Ugas_J_kg + &
                gen_auxvar%pres(gid)/gen_auxvar%den_kg(gid)
#endif

    ! MJ/kmol
    gen_auxvar%U(gid) = xmol_water_in_gas * u_water_vapor + &
                        xmol_air_in_gas * u_air
    Hg_mixture_fractioned = xmol_water_in_gas*h_water_vapor + &
                            xmol_air_in_gas*h_air
    gen_auxvar%H(gid) = gen_auxvar%U(gid) + &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        gen_auxvar%pres(gid)/gen_auxvar%den(gid) * 1.d-6

  endif ! istate /= LIQUID_STATE

  if (global_auxvar%hstate == LIQUID_STATE .or. &
      global_auxvar%istate == TWO_PHASE_STATE) then
    call characteristic_curves%liq_rel_perm_function% &
           RelativePermeability(gen_auxvar%sat(lid),krl,dkrl_dsatl,option)
    dkrl_dsatg = -1.d0 * dkrl_dsatl
    if (.not.option%flow%density_depends_on_salinity) then
      call EOSWaterViscosity(gen_auxvar%temp,cell_pressure, &
                               gen_auxvar%pres(spid),visl,ierr)
    else
      aux(1) = global_auxvar%m_nacl(1)
      call EOSWaterViscosityExt(gen_auxvar%temp,cell_pressure, &
                                  gen_auxvar%pres(spid),aux,visl,ierr)
    endif
    gen_auxvar%mobility(lid) = krl/visl
  endif

  if (global_auxvar%hstate == GA_STATE) then
    call characteristic_curves%gas_rel_perm_function% &
           RelativePermeability(gen_auxvar%sat(lid),krg,dkrg_dsatl,option)
    dkrg_dsatg = -1.d0 * dkrg_dsatl
    call EOSGasViscosity(gen_auxvar%temp,gen_auxvar%pres(apid), &
                           gen_auxvar%pres(gid),den_air,visg,ierr)
    gen_auxvar%mobility(gid) = krg/visg
  endif

  call EOSHydrateEnergy(gen_auxvar%temp, U_hyd)
  gen_auxvar%xmol(wid,hid) = MOL_RATIO_H20
  gen_auxvar%xmol(gid,hid) = MOL_RATIO_METH
  gen_auxvar%den(hid) = HYDRATE_DENSITY
  gen_auxvar%den_kg(hid) = HYDRATE_DENSITY_KG
  gen_auxvar%U(hid) = U_hyd
  gen_auxvar%H(hid) = U_hyd
  gen_auxvar%mobility(hid) = 0.d0
  
  call EOSIceEnergy(gen_auxvar%temp, U_ice)
  gen_auxvar%xmol(wid,iid) = 1.d0
  gen_auxvar%xmol(gid,iid) = 0.d0
  call EOSWaterDensityIcePainter(gen_auxvar%temp,gen_auxvar%pres(lid), &
                    PETSC_FALSE, gen_auxvar%den(iid), &
                    dden_ice_dT, dden_ice_dP, ierr)
  !gen_auxvar%den(iid) = ICE_DENSITY
  gen_auxvar%den_kg(iid) = ICE_DENSITY_KG
  gen_auxvar%U(iid) = U_ice
  gen_auxvar%H(iid) = U_ice
  gen_auxvar%mobility(iid) = 0.d0

end subroutine HydrateAuxVarCompute

subroutine HydrateAccumulation(gen_auxvar,global_auxvar,material_auxvar, &
                               soil_heat_capacity,option,Res,Jac, &
                               analytical_derivatives,debug_cell)
  !
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  !

  use Option_module
  use Material_Aux_class

  implicit none

  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_cell

  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase

  PetscReal :: porosity
  PetscReal :: volume_over_dt

  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use gen_auxvar%effective porosity here as it enables numerical
  ! derivatives to be employed
  porosity = gen_auxvar%effective_porosity

  ! accumulation term units = kmol/s
  Res = 0.d0
  do iphase = 1, option%nphase
    ! Res[kmol comp/m^3 void] = sat[m^3 phase/m^3 void] *
    !                           den[kmol phase/m^3 phase] *
    !                           xmol[kmol comp/kmol phase]
    do icomp = 1, option%nflowspec
#ifdef DEBUG_GENERAL
      ! for debug version, aux var entries are initialized to NaNs.  even if
      ! saturation is zero, density may be a NaN.  So the conditional prevents
      ! this calculation.  For non-debug, aux var entries are initialized to
      ! 0.d0
      if (gen_auxvar%sat(iphase) > 0.d0) then
#endif
      Res(icomp) = Res(icomp) + gen_auxvar%sat(iphase) * &
                                gen_auxvar%den(iphase) * &
                                gen_auxvar%xmol(icomp,iphase)
#ifdef DEBUG_GENERAL
      endif
#endif
    enddo
  enddo

  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] *
  !                 vol[m^3 bulk] / dt[sec]
  Res(1:option%nflowspec) = Res(1:option%nflowspec) * &
                            porosity * volume_over_dt

  do iphase = 1, option%nphase
    ! Res[MJ/m^3 void] = sat[m^3 phase/m^3 void] *
    !                    den[kmol phase/m^3 phase] * U[MJ/kmol phase]
#ifdef DEBUG_GENERAL
    ! for debug version, aux var entries are initialized to NaNs.  even if
    ! saturation is zero, density may be a NaN.  So the conditional prevents
    ! this calculation.  For non-debug, aux var entries are initialized to
    ! 0.d0
    if (gen_auxvar%sat(iphase) > 0.d0) then
#endif
    Res(energy_id) = Res(energy_id) + gen_auxvar%sat(iphase) * &
                                      gen_auxvar%den(iphase) * &
                                      gen_auxvar%U(iphase)
#ifdef DEBUG_GENERAL
    endif
#endif
  enddo
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] +
  !                (1-por)[m^3 rock/m^3 bulk] *
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * porosity + &
                    (1.d0 - porosity) * &
                    material_auxvar%soil_particle_density * &
                    soil_heat_capacity * gen_auxvar%temp) * volume_over_dt

#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,7es24.15)') 'accum:', Res
  endif
#endif

end subroutine HydrateAccumulation


subroutine HydrateAuxVarPerturb(gen_auxvar,global_auxvar, &
                                material_auxvar, &
                                characteristic_curves,natural_id, &
                                option)
  !
  ! Calculates auxiliary variables for perturbed system
  !
  ! Author: Michael Nole
  ! Date: 01/30/19
  !
  
  ! MAN: This subroutine needs work

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(general_auxvar_type) :: gen_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), &
               pert(option%nflowdof), x_pert_save(option%nflowdof)

  PetscReal :: tempreal
  PetscInt :: lid, gid, hid
!#define LEGACY_PERTURBATION
#ifdef LEGACY_PERTURBATION
  PetscReal, parameter :: perturbation_tolerance = 1.d-5
#else
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
! 1.d-11 works well for Emily's 1D nacl2
!  PetscReal, parameter :: perturbation_tolerance = 1.d-11
#endif
  PetscReal, parameter :: min_mole_fraction_pert = 1.d-12
  PetscReal, parameter :: min_perturbation = 1.d-10
  PetscInt :: idof

#ifdef DEBUG_GENERAL
  type(global_auxvar_type) :: global_auxvar_debug
  type(general_auxvar_type) :: general_auxvar_debug
  call GlobalAuxVarInit(global_auxvar_debug,option)
  call GeneralAuxVarInit(general_auxvar_debug,PETSC_FALSE,option)
#endif

  lid = 1
  gid = 2
  hid = 3

  select case(global_auxvar%hstate)
    case(LIQUID_STATE)
       x(GENERAL_LIQUID_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%liquid_phase)
       x(GENERAL_LIQUID_STATE_X_MOLE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%xmol(option%air_id,option%liquid_phase)
       x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%temp
#ifdef LEGACY_PERTURBATION
       ! if the liquid state, the liquid pressure will always be greater
       ! than zero.
       pert(GENERAL_LIQUID_PRESSURE_DOF) = &
         max(perturbation_tolerance*x(GENERAL_LIQUID_PRESSURE_DOF), &
             perturbation_tolerance)
       ! if the air mole fraction perturbation is too small, the derivatives
       ! can be poor.
       pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = &
         -1.d0*max(perturbation_tolerance*x(GENERAL_LIQUID_STATE_X_MOLE_DOF), &
                   min_mole_fraction_pert)
       pert(GENERAL_ENERGY_DOF) = &
         -1.d0*perturbation_tolerance*x(GENERAL_ENERGY_DOF)
#else
       pert(GENERAL_LIQUID_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_LIQUID_PRESSURE_DOF) + &
         min_perturbation
       if (x(GENERAL_LIQUID_STATE_X_MOLE_DOF) > &
           1.d3 * perturbation_tolerance) then
         pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = -1.d0 * perturbation_tolerance
       else
         pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = perturbation_tolerance
       endif
       pert(GENERAL_ENERGY_DOF) = -1.d0 * &
         (perturbation_tolerance*x(GENERAL_ENERGY_DOF) + min_perturbation)
#endif
    case(GAS_STATE)
       x(GENERAL_GAS_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
       x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
       x(GENERAL_ENERGY_DOF) = gen_auxvar(ZERO_INTEGER)%temp
#ifdef LEGACY_PERTURBATION
       ! gas pressure [p(g)] must always be perturbed down as p(v) = p(g) - p(a)
       ! and p(v) >= Psat (i.e. an increase in p(v)) results in two phase.
       pert(GENERAL_GAS_PRESSURE_DOF) = &
         -1.d0*perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)
       ! perturb air pressure towards gas pressure unless the perturbed
       ! air pressure exceeds the gas pressure
       if (x(GENERAL_GAS_PRESSURE_DOF) - &
           x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) > &
           perturbation_tolerance* &
           x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)) then
         pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
           perturbation_tolerance*x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
       else
         pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
           -1.d0*perturbation_tolerance*x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)
       endif
       pert(GENERAL_ENERGY_DOF) = &
         perturbation_tolerance*x(GENERAL_ENERGY_DOF)
#else
       ! gas pressure [p(g)] must always be perturbed down as p(v) = p(g) - p(a)
       ! and p(v) >= Psat (i.e. an increase in p(v)) results in two phase.
       pert(GENERAL_GAS_PRESSURE_DOF) = -1.d0 * &
         (perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF) + min_perturbation)
       ! perturb air pressure towards gas pressure unless the perturbed
       ! air pressure exceeds the gas pressure
       tempreal = perturbation_tolerance* &
                  x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) + min_perturbation
       if (x(GENERAL_GAS_PRESSURE_DOF) - &
           x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) > tempreal) then
         pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = tempreal
       else
         pert(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
       endif
       pert(GENERAL_ENERGY_DOF) = &
         perturbation_tolerance*x(GENERAL_ENERGY_DOF) + min_perturbation
#endif
    case(H_STATE)
   
    case(ICE_STATE)        
       x(GENERAL_GAS_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
       x(GENERAL_GAS_SATURATION_DOF) = 0.d0
       x(GENERAL_ENERGY_DOF) = gen_auxvar(ZERO_INTEGER)%temp

       pert(GENERAL_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)
       pert(GENERAL_GAS_SATURATION_DOF) = 0.d0
       pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)
       
    case(GA_STATE)
       x(GENERAL_GAS_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
!       x(GENERAL_AIR_PRESSURE_DOF) = &
!         gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
       x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(option%gas_phase)
#ifdef LEGACY_PERTURBATION
       if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%temp
         pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)
       else
         ! here GENERAL_2PH_STATE_AIR_PRESSURE_DOF = GENERAL_ENERGY_DOF
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
         ! perturb air pressure towards gas pressure unless the perturbed
         ! air pressure exceeds the gas pressure
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
             x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) > &
             perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)) then
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
             perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)
         else
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
             -1.d0*perturbation_tolerance*x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)
         endif
       endif
       pert(GENERAL_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)
       ! always perturb toward 0.5
       if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
         pert(GENERAL_GAS_SATURATION_DOF) = &
           -1.d0*perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF)
       else
         pert(GENERAL_GAS_SATURATION_DOF) = &
           perturbation_tolerance*x(GENERAL_GAS_SATURATION_DOF)
       endif
#else
       if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%temp
         pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation
       else
         ! here GENERAL_2PH_STATE_AIR_PRESSURE_DOF = GENERAL_ENERGY_DOF
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
         ! perturb air pressure towards gas pressure unless the perturbed
         ! air pressure exceeds the gas pressure
         tempreal = perturbation_tolerance* &
                    x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) + min_perturbation
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
             x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) > tempreal) then
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = tempreal
         else
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
         endif
       endif
       pert(GENERAL_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)+min_perturbation
       if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
         pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
       else
         pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance
       endif
#endif
    case(HG_STATE)

    case(HA_STATE)
      x(GENERAL_GAS_PRESSURE_DOF) = &
        gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(hid)
      if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%temp
         pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation
      else
         ! here GENERAL_2PH_STATE_AIR_PRESSURE_DOF = GENERAL_ENERGY_DOF
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
         ! perturb air pressure towards gas pressure unless the perturbed
         ! air pressure exceeds the gas pressure
         tempreal = perturbation_tolerance* &
                    x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) + min_perturbation
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
             x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) > tempreal) then
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = tempreal
         else
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
         endif
      endif
        pert(GENERAL_GAS_PRESSURE_DOF) = &
        perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)+min_perturbation
      if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
        pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance
      endif

    case(HI_STATE)

    case(GI_STATE)

    case(AI_STATE)
      x(GENERAL_LIQUID_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%liquid_phase)
      x(GENERAL_LIQUID_STATE_X_MOLE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%xmol(option%air_id,option%liquid_phase)
      x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(lid)
      pert(GENERAL_LIQUID_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_LIQUID_PRESSURE_DOF) + &
         min_perturbation
      if (x(GENERAL_LIQUID_STATE_X_MOLE_DOF) > &
           1.d3 * perturbation_tolerance) then
        pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_LIQUID_STATE_X_MOLE_DOF) = perturbation_tolerance
      endif
      if (x(GENERAL_ENERGY_DOF) > 0.5d0) then
        pert(GENERAL_ENERGY_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_ENERGY_DOF) = perturbation_tolerance
      endif

    case(HGA_STATE)
      x(GENERAL_GAS_PRESSURE_DOF) = &
        gen_auxvar(ZERO_INTEGER)%sat(lid)
      x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(hid)
      if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%temp
         pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation
      else
         ! here GENERAL_2PH_STATE_AIR_PRESSURE_DOF = GENERAL_ENERGY_DOF
         x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
         ! perturb air pressure towards gas pressure unless the perturbed
         ! air pressure exceeds the gas pressure
         tempreal = perturbation_tolerance* &
                    x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) + min_perturbation
         if (x(GENERAL_GAS_PRESSURE_DOF) - &
             x(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) > tempreal) then
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = tempreal
         else
           pert(GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
         endif
      endif
      
      if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
        pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      if (x(GENERAL_GAS_PRESSURE_DOF) > 0.5d0) then
        pert(GENERAL_GAS_PRESSURE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_PRESSURE_DOF) = perturbation_tolerance
      endif

    case(HAI_STATE)
      
      x(GENERAL_GAS_PRESSURE_DOF) = &
        gen_auxvar(ZERO_INTEGER)%pres(gid)
      x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(lid)   
      x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(iid)
      
      pert(GENERAL_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF) + &
         min_perturbation 
      if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
        pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      if (x(GENERAL_ENERGY_DOF) > 0.5d0) then
        pert(GENERAL_ENERGY_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_ENERGY_DOF) = perturbation_tolerance
      endif

    case(HGI_STATE)

    case(GAI_STATE)

    case(QUAD_STATE)

  end select

  ! GENERAL_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = GENERAL_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    gen_auxvar(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    x_pert_save = x_pert
    call HydrateAuxVarCompute(x_pert,gen_auxvar(idof),global_auxvar, &
                              material_auxvar, &
                              characteristic_curves,natural_id,option)
  enddo

  select case(global_auxvar%hstate)
    case(LIQUID_STATE)
      gen_auxvar(GENERAL_LIQUID_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_LIQUID_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(GAS_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
      gen_auxvar(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_STATE_AIR_PRESSURE_DOF)%pert / &
        GENERAL_PRESSURE_SCALE
    case(H_STATE)

    case(ICE_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(GA_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
      if (general_2ph_energy_dof == GENERAL_AIR_PRESSURE_INDEX) then
        gen_auxvar(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)%pert = &
          gen_auxvar(GENERAL_2PH_STATE_AIR_PRESSURE_DOF)%pert / &
          GENERAL_PRESSURE_SCALE
      endif
    case(HG_STATE)
    case(HA_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(HI_STATE)

    case(GI_STATE)

    case(AI_STATE)
      gen_auxvar(GENERAL_LIQUID_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_LIQUID_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(HGA_STATE)

    case(HAI_STATE)

    case(HGI_STATE)

    case(GAI_STATE)

    case(QUAD_STATE)
  end select


end subroutine HydrateAuxVarPerturb

subroutine HydratePE(T, PE)

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: PE

  PetscReal :: T_temp
  PetscReal, parameter :: TQD = 273.15 ! K

  T_temp = T + TQD

  if (T_temp < TQD) then
    PE = exp(-43.8921173434628 + 0.776302133739303 * T_temp &
          - 7.27291427030502d-3 * T_temp**2 + 3.85413985900724d-5 * T_temp**3 &
          - 1.03669656828834d-7 * T_temp**4 + 1.09882180475307d-10 * T_temp**5)
  else
    PE = exp(-1.9413850446456d5 + 3.31018213397926d3 * T_temp &
          - 22.5540264493806* T_temp**2 + 0.0767559117787059 * T_temp**3 &
          - 1.30465829788791d-4 * T_temp**4 + 8.86065316687571d-8 * T_temp**5)
  endif

  PE = PE * 1.d6

end subroutine HydratePE

subroutine HenrysConstantMethane(T,K_H)

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: K_H

  PetscReal, parameter :: TQD = 273.15 ! K
  PetscReal :: T_temp

  T_temp = T
  T_temp = T_temp + TQD

  ! Units: Pa/M
  K_H = exp(5.1345 + 7837/T_temp - 1.509d6/(T_temp**2) + 2.06d7/(T_temp**3)) *1000
 
end subroutine HenrysConstantMethane

subroutine EOSHydrateEnergy(T,U)

  implicit none

  PetscReal, intent(in):: T
  PetscReal, intent(out) :: U

  PetscReal, parameter :: TQD = 237.15d0 ! K
  PetscReal, parameter :: Hh0 = -54734d0 ! J/mol
  PetscReal, parameter :: MWH = 82.187d0 ! g/mol
  PetscReal :: Cph, T_temp

  T_temp = T + TQD
  
  ! Integral of Cph * dT ; Cph from Handa, 1998

  ! Units: J/mol
  U = Hh0 + 6.6d0 * (T_temp-TQD) + 7.269d-1 * (T_temp**2 - TQD **2) - 1.21333d-3 * &
        (T_temp**3 - TQD**3)  + 1.578d-6 * (T_temp**4 - TQD**4)
  ! Units: MJ/kmol
  U = U / 1.d3  

end subroutine EOSHydrateEnergy

subroutine EOSIceEnergy(T,U)

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: U
  PetscReal, parameter :: Lw = -6017.1d0 !Latent heat of fusion,  J/mol
  PetscReal, parameter :: TQD = 273.15d0
  PetscReal :: T_temp

  T_temp = T + TQD

  if (T_temp >= 90.d0) then
    U = Lw + 185.d0 * (T_temp-TQD) + 3.445 * (T_temp*T_temp - TQD*TQD) 
  else
    U = Lw + 4.475 * (T_temp*T_temp - TQD*TQD)
  endif

  ! J/mol to MJ/kmol
  U = U / 1.d3
       
end subroutine EOSIceEnergy

subroutine GibbsThomsonFreezing(sat,Hf,rho,Tb,dTf,characteristic_curves,option)

  use Characteristic_Curves_module
  use Option_module

  implicit none

  PetscReal, intent(in) :: sat
  PetscReal, intent(in) :: Hf
  PetscReal, intent(in) :: rho
  PetscReal, intent(in) :: Tb
  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal, intent(out) :: dTf

  PetscReal :: Pc,dpc_dsatl
  
  call characteristic_curves%saturation_function% &
             CapillaryPressure(sat,Pc,dpc_dsatl,option)

  dTf = -(Tb+273.15)*Pc/(Hf * rho * 1000.d0)
  

end subroutine GibbsThomsonFreezing

!subroutine Methanogenesis(depth, q_meth)
  ! A simple methanogenesis source parameterized as a function of depth
  
!#include "petsc/finclude/petscvec.h"
!  use petscvec
!  use Realization_Subsurface_class
  
!  implicit none
  
!  type(criticality_mediator_type), pointer :: this
!  class(realization_subsurface_type), pointer :: realization
!  PetscReal :: time
!  PetscErrorCode :: ierr
  
!  PetscInt :: i,j
!  type(criticality_type), pointer :: cur_criticality
!  PetscReal, pointer :: heat_source(:)
  
!  PetscReal, intent(in) :: depth
!  PetscReal, intent(out) :: q_meth
  
!  PetscReal, parameter :: alpha = 0.005
!  PetscReal, parameter :: k_alpha = 2241 ! Maliverno, 2010, corrected
!  PetscReal, parameter :: lambda = 1.d-13
!  PetscReal, parameter :: omega = 3.17d-11
!  PetscReal, parameter :: z_smt = 15.d0
  
  
!  if (depth > zsmt) then
!    q_meth = k_alpha * lamda * exp(-lambda/omega * (depth - zsmt))
!  endif
  
!  call VecGetArrayF90(this%data_mediator%vec,heat_source, &
!                      ierr);CHKERRQ(ierr)
  
!  cur_criticality => this%criticality_list
!  j = 0
!  do
!    if (.not. associated(cur_criticality)) exit
!    do i = 1, cur_criticality%region%num_cells
!      j = j + 1
!      heat_source(j) = -cur_criticality%crit_mech%heat_released
!    enddo
!    cur_criticality => cur_criticality%next
!  enddo
  
!  call VecRestoreArrayF90(this%data_mediator%vec,heat_source, &
!                          ierr);CHKERRQ(ierr)
  
  
!end subroutine Methanogenesis

end module Hydrate_module
