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

  PetscInt, parameter, public :: H_STATE = 17 !5 (4 and 5 conflict with 
  PetscInt, parameter, public :: ICE_STATE = 16 !4 ANY_STATE and MULTI_STATE)
  PetscInt, parameter, public :: GA_STATE = 3
  PetscInt, parameter, public :: HG_STATE = 6
  PetscInt, parameter, public :: HA_STATE = 7
  PetscInt, parameter, public :: HI_STATE = 8
  PetscInt, parameter, public :: GI_STATE = 9
  PetscInt, parameter, public :: AI_STATE = 10
  PetscInt, parameter, public :: HGA_STATE = 11
  PetscInt, parameter, public :: HAI_STATE = 12
  PetscInt, parameter, public :: HGI_STATE = 13
  PetscInt, parameter, public :: GAI_STATE = 14
  PetscInt, parameter, public :: QUAD_STATE = 15
  
  PetscInt, parameter :: lid = 1
  PetscInt, parameter :: gid = 2
  PetscInt, parameter :: hid = 3
  PetscInt, parameter :: iid = 4

  !Structure 1 methane hydrate:
  PetscReal, parameter :: Nhyd = 6.d0
  PetscReal, parameter :: HYDRATE_DENSITY_KG = 920.d0 !kg/m^3
  PetscReal, parameter :: HYDRATE_DENSITY = 52.15551276d0 !mol/L
  PetscReal, parameter :: MW_CH4 = 16.04d0
  PetscReal, parameter :: MW_H20 = 18.01d0 

  PetscReal, parameter :: MOL_RATIO_METH = 0.14285714285d0
  PetscReal, parameter :: MOL_RATIO_H20 = 1.d0 - MOL_RATIO_METH

  PetscReal, parameter :: TQD = 1.d-2 !0.d0 !1.0d-2 !Quad point temperature (C)
  
  !Ice: 
  PetscReal, parameter :: ICE_DENSITY_KG = 920.d0 !kg/m^3
  PetscReal, parameter :: ICE_DENSITY = 50.86d0 !mol/L


  PetscReal, parameter :: lambda_hyd = 0.49d0 !W/m-K

  PetscReal :: hydrate_perm_base(3) = -999.9d0
  PetscInt :: hydrate_perm_scaling_function = 0

  public :: HydrateSetFlowMode, &
            HydrateRead, &
            HydrateUpdateState, &
            HydrateAuxVarCompute, &
            HydrateAccumulation, &
            HydrateAuxVarPerturb, &
            HydrateCompositeThermalCond, &
            HydratePE, &
            Methanogenesis

contains

! ************************************************************************** !

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

end subroutine HydrateSetFlowMode

! ************************************************************************** !

subroutine HydrateRead(input,meth,option)
  
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none
   
  type(input_type), pointer :: input
  type(methanogenesis_type), pointer :: meth
  type(option_type), pointer :: option

  character(len=MAXSTRINGLENGTH) :: error_string 
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: temp_int

  do

    call InputReadPflotranString(input,option)
    if (input%ierr /= 0) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE) 
    call InputErrorMsg(input,option,'keyword','HYDRATE')
    call StringToUpper(word)

    select case(trim(word))
      case('METHANOGENESIS')
        if (.not. associated(meth)) then
          allocate(meth)
        endif
        do
          call InputReadPflotranString(input,option)
          if (input%ierr /= 0) exit
          if (InputCheckExit(input,option)) exit
          
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','HYDRATE')
          call StringToUpper(word)
          select case(trim(word))
            case('NAME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'methanogenesis source name', &
                                 error_string)
              call StringToUpper(word)
              meth%source_name = trim(word)
            case('ALPHA')
              call InputReadDouble(input,option,meth%alpha)
              call InputErrorMsg(input,option,'alpha',error_string) 
            case('LAMBDA')
              call InputReadDouble(input,option,meth%lambda)
              call InputErrorMsg(input,option,'lambda',error_string) 
            case('V_SED')
              call InputReadDouble(input,option,meth%omega)
              call InputErrorMsg(input,option,'v_sed',error_string) 
            case('SMT_DEPTH')
              call InputReadDouble(input,option,meth%z_smt)
              call InputErrorMsg(input,option,'smt_depth',error_string) 
            case('K_ALPHA')
              call InputReadDouble(input,option,meth%k_alpha)
              call InputErrorMsg(input,option,'k_alpha',error_string) 
          end select
        enddo
      case('PERM_SCALING_FUNCTION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','hyd_perm_scaling_function')
        call StringToUpper(word)
        select case(trim(word))
          case('DAI_AND_SEOL')
            temp_int = 1
            hydrate_perm_scaling_function = temp_int
        end select
    end select

  enddo 

end subroutine HydrateRead

! ************************************************************************** !

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
  PetscReal :: liq_epsilon, gas_epsilon, hyd_epsilon, two_phase_epsilon
  PetscReal :: ga_epsilon, ha_epsilon
  PetscReal :: x(option%nflowdof)
  PetscReal :: PE_hyd, K_H, Tf_ice, dTf, h_sat_eff, i_sat_eff
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

  !man: right now comparing hydrate equilib pressure to gas
  !pressure (assuming low water solubility in methane). 
  !Ideally would compare to partial pressure of methane.

  if (global_auxvar%hstate == ZERO_INTEGER .and. gen_auxvar%sat(gid) &
       < 0.d0) then
    global_auxvar%hstate = HA_STATE
    gen_auxvar%sat(hid) = -1.d0 * gen_auxvar%sat(gid)
    gen_auxvar%sat(gid) = 0.d0
  endif

  if (global_auxvar%hstate == ZERO_INTEGER) global_auxvar% &
                              hstate = global_auxvar%istate 
  gen_auxvar%hstate_store(PREV_IT) = global_auxvar%hstate
  
  if (gen_auxvar%sat(hid) > gen_auxvar%sat(iid)) then
    h_sat_eff = gen_auxvar%sat(hid)+gen_auxvar%sat(iid)
    i_sat_eff = 2.d0 * gen_auxvar%sat(iid)
  else
    h_sat_eff = 2.d0 * gen_auxvar%sat(hid)
    i_sat_eff = gen_auxvar%sat(hid) + gen_auxvar%sat(iid)
  endif

  call HydratePE(gen_auxvar%temp,h_sat_eff, PE_hyd, &
          characteristic_curves, option)

  call GibbsThomsonFreezing(1.d0-i_sat_eff,6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,option)

  Tf_ice = TQD + dTf
  !Update State
  
  select case(global_auxvar%hstate)
    case(LIQUID_STATE)
      if (gen_auxvar%temp > Tf_ice) then
        if (gen_auxvar%pres(apid) >= gen_auxvar% &
             pres(lid)*(1.d0-window_epsilon)) then
          !if (gen_auxvar%pres(apid) >= PE_hyd) then
          if (gen_auxvar%pres(gid) >= PE_hyd) then
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
      elseif (gen_auxvar%pres(apid) >= gen_auxvar%pres(lid)* &
              (1.d0-window_epsilon)) then
        if (gen_auxvar%pres(gid) < PE_hyd) then
      !elseif (gen_auxvar%pres(apid) >= gen_auxvar%pres(lid)) then
        !if (gen_auxvar%pres(apid) < PE_hyd) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = GAI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
          liq_epsilon = option%phase_chng_epsilon
        elseif (gen_auxvar%pres(gid) > PE_hyd) then
          istatechng = PETSC_TRUE
          global_auxvar%hstate = HAI_STATE
          global_auxvar%istate = TWO_PHASE_STATE
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
      !if (gen_auxvar%pres(gid) < PE_hyd) then
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
      !if (gen_auxvar%pres(apid) < PE_hyd) then
      if (gen_auxvar%pres(gid) < PE_hyd) then
        if (gen_auxvar%temp > Tf_ice) then
          if (gen_auxvar%sat(gid) > 0.d0 .and. gen_auxvar%sat(lid) > 0.d0) then
            istatechng = PETSC_FALSE
          elseif (gen_auxvar%sat(gid) <= 0.d0) then
            istatechng = PETSC_TRUE
            global_auxvar%hstate = LIQUID_STATE
            global_auxvar%istate = LIQUID_STATE
            two_phase_epsilon = option%phase_chng_epsilon
          elseif (gen_auxvar%sat(gid) >= 1.d0) then
            istatechng = PETSC_TRUE
            global_auxvar%hstate = GAS_STATE
            global_auxvar%istate = GAS_STATE
            two_phase_epsilon = option%phase_chng_epsilon
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
      !if (gen_auxvar%pres(apid) > PE_hyd) then
      if (gen_auxvar%pres(gid) > PE_hyd) then
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
      !if (gen_auxvar%pres(gid) > PE_hyd .and. gen_auxvar%temp > Tf_ice) then
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
      !if (gen_auxvar%pres(gid) > PE_hyd) then
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
      !if (gen_auxvar%temp < Tf_ice .and. gen_auxvar%pres(gid) < PE_hyd) then
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
      !if (gen_auxvar%pres(apid) >= gen_auxvar% &
      !       pres(lid)*(1.d0-window_epsilon)) then 
      !  if (gen_auxvar%pres(apid) < PE_hyd) then
      if (gen_auxvar%pres(apid) > PE_hyd*(1.d0-window_epsilon)) then
        if (gen_auxvar%pres(lid) < PE_hyd*(1.d0+window_epsilon)) then
      !  if (gen_auxvar%pres(lid) >  PE_hyd*(1.d0-window_epsilon)) then
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
      !if (gen_auxvar%pres(apid) > PE_hyd) then
      if (gen_auxvar%pres(gid) > PE_hyd*(1.d0-window_epsilon)) then
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
      !if (gen_auxvar%pres(gid) < PE_hyd) then
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
        call PrintErrMsgByRank(option)

    end select

    call HydrateAuxVarCompute(x,gen_auxvar, global_auxvar,material_auxvar, &
          characteristic_curves,natural_id,option)

  endif

end subroutine HydrateUpdateState

! ************************************************************************** !

subroutine HydrateAuxVarCompute(x,gen_auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell, with gas hydrate physics
  ! Author: Michael Nole
  ! Date: 04/04/19
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id

  PetscInt :: gid, lid, acid, wid, eid, hid, iid
  PetscReal :: cell_pressure, water_vapor_pressure
  PetscReal :: den_water_vapor, den_kg_water_vapor
  PetscReal :: u_water_vapor, h_water_vapor
  PetscReal :: den_air, h_air, u_air
  PetscReal :: xmol_air_in_gas, xmol_water_in_gas
  PetscReal :: krl, visl
  PetscReal :: dkrl_dsatl, dkrl_dsatg
  PetscReal :: dkrg_dsatl, dkrg_dsatg
  PetscReal :: krg, visg
  PetscReal :: K_H_tilde
  PetscInt :: apid, cpid, vpid, spid
  PetscReal :: xmass_air_in_gas
  PetscReal :: Ugas_J_kg, Hgas_J_kg
  PetscReal :: Uair_J_kg, Hair_J_kg
  PetscReal :: Uvapor_J_kg, Hvapor_J_kg
  PetscReal :: Hg_mixture_fractioned
  PetscReal :: H_hyd, U_ice, PE_hyd
  PetscReal :: aux(1)
  PetscReal :: hw, hw_dp, hw_dT
  PetscReal :: dpor_dp
  PetscReal :: dpc_dsatl
  PetscReal :: dden_ice_dT, dden_ice_dP
  character(len=8) :: state_char
  PetscErrorCode :: ierr
  PetscReal :: dTf, h_sat_eff, i_sat_eff, liq_sat_eff, g_sat_eff
  PetscReal :: solid_sat_eff
  PetscReal :: sigma

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id
  wid = option%water_id
  eid = option%energy_id


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
  gen_auxvar%mobility = 0.d0

#if 0
  if (option%iflag >= GENERAL_UPDATE_FOR_ACCUM) then
    if (option%iflag == GENERAL_UPDATE_FOR_ACCUM) then
      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
        natural_id, x(1:3), trim(state_char)
    else
    endif
  endif
#endif

  if (global_auxvar%hstate == ZERO_INTEGER) global_auxvar% &
                              hstate = global_auxvar%istate
  
  gen_auxvar%xmol(wid,hid) = MOL_RATIO_H20
  gen_auxvar%xmol(acid,hid) = MOL_RATIO_METH
  gen_auxvar%den(hid) = HYDRATE_DENSITY
  gen_auxvar%den_kg(hid) = HYDRATE_DENSITY_KG

  select case(global_auxvar%hstate)
    case(LIQUID_STATE)
!     ********* Aqueous State (A) ********************************
!     Primary variables: Pa, Xma, T
!
      gen_auxvar%pres(lid) = x(GENERAL_LIQUID_PRESSURE_DOF)
      gen_auxvar%xmol(acid,lid) = x(GENERAL_LIQUID_STATE_X_MOLE_DOF)
      gen_auxvar%temp = x(GENERAL_ENERGY_DOF)

      gen_auxvar%xmol(acid,lid) = max(0.d0,gen_auxvar%xmol(acid,lid))

      gen_auxvar%xmol(wid,lid) = 1.d0 - gen_auxvar%xmol(acid,lid)
      gen_auxvar%xmol(wid,gid) = 0.d0
      gen_auxvar%xmol(acid,gid) = 0.d0
      gen_auxvar%sat(lid) = 1.d0
      gen_auxvar%sat(gid) = 0.d0
      gen_auxvar%sat(hid) = 0.d0
      gen_auxvar%sat(iid) = 0.d0

      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                        gen_auxvar%pres(spid),ierr)
      call HydratePE(gen_auxvar%temp, 0.d0, PE_hyd, characteristic_curves, &
                     option)
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      call HydrateDavieBuffettCorrection(gen_auxvar%temp,gen_auxvar% &
                                           pres(lid),K_H_tilde)
      gen_auxvar%pres(spid) = 1.d-6 
      
      gen_auxvar%pres(gid) = max(gen_auxvar%pres(lid),gen_auxvar%pres(spid))
      gen_auxvar%pres(apid) = K_H_tilde*gen_auxvar%xmol(acid,lid)

      if (gen_auxvar%pres(gid) <= 0.d0) then
        write(option%io_buffer,'(''Negative gas pressure at cell '', &
          & i8,'' in HydrateAuxVarCompute(LIQUID_STATE).  Attempting bailout.'')') &
          natural_id
        call PrintMsgByRank(option)
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

      call HydratePE(gen_auxvar%temp,gen_auxvar%sat(hid),PE_hyd, &
              characteristic_curves, option)
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
     
      gen_auxvar%sat(gid) = max(0.d0,min(1.d0,gen_auxvar%sat(gid))) 
      gen_auxvar%sat(lid) = 1.d0 - gen_auxvar%sat(gid)
      gen_auxvar%sat(hid) = 0.d0
      gen_auxvar%sat(iid) = 0.d0

      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                          gen_auxvar%pres(spid),ierr)
      gen_auxvar%pres(spid) = 1.d-6 

      if (general_immiscible) then
        gen_auxvar%pres(spid) = GENERAL_IMMISCIBLE_VALUE
      endif
      gen_auxvar%pres(vpid) = gen_auxvar%pres(spid)
      gen_auxvar%pres(apid) = gen_auxvar%pres(gid) - gen_auxvar%pres(vpid)

      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)

      call characteristic_curves%saturation_function% &
             CapillaryPressure(gen_auxvar%sat(lid), gen_auxvar%pres(cpid), &
                               dpc_dsatl,option)

      !IFT calculation
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
     
      if (gen_auxvar%sat(hid) > gen_auxvar%sat(gid)) then
        call HydratePE(gen_auxvar%temp, gen_auxvar%sat(hid)+ &
                gen_auxvar%sat(gid), PE_hyd, characteristic_curves, option)
      else
        call HydratePE(gen_auxvar%temp, 2.d0 * gen_auxvar%sat(hid), &
                PE_hyd, characteristic_curves, option)
      endif

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

      gen_auxvar%sat(hid) = max(0.d0,min(1.d0,gen_auxvar%sat(hid)))

      gen_auxvar%sat(lid) = 1.d0 - gen_auxvar%sat(hid)
      gen_auxvar%sat(gid) = 0.d0
      gen_auxvar%sat(iid) = 0.d0
      
      call HydratePE(gen_auxvar%temp,gen_auxvar%sat(hid), PE_hyd, &
              characteristic_curves, option)
      gen_auxvar%pres(apid) = PE_hyd
      
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      call HydrateDavieBuffettCorrection(gen_auxvar%temp,gen_auxvar% &
                                           pres(gid),K_H_tilde)
      
      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                          gen_auxvar%pres(spid),ierr)


      gen_auxvar%pres(cpid) = 0.d0
      ! Setting air pressure equal to gas pressure makes forming hydrate
      ! easier
      gen_auxvar%pres(apid) = gen_auxvar%pres(gid)
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid)
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
      
      if (gen_auxvar%sat(hid) > gen_auxvar%sat(iid)) then
        call HydratePE(gen_auxvar%temp, gen_auxvar%sat(hid)+ &
                gen_auxvar%sat(iid), PE_hyd, characteristic_curves, option)
      else
        call HydratePE(gen_auxvar%temp, 2.d0 * gen_auxvar%sat(hid), PE_hyd, &
                characteristic_curves, option)
      endif

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
     
      gen_auxvar%xmol(acid,lid) = max(0.d0,gen_auxvar%xmol(acid,lid))

      if (global_auxvar%istatechng) then
        gen_auxvar%sat(lid) = max(0.d0,min(1.d0,gen_auxvar%sat(lid)))
      endif

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
      gen_auxvar%pres(vpid) = 0.d0 !gen_auxvar%pres(gid) - gen_auxvar%pres(apid)

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
      
      gen_auxvar%sat(lid) = max(0.d0,min(1.d0,gen_auxvar%sat(lid)))
      gen_auxvar%sat(hid) = max(0.d0,min(1.d0,gen_auxvar%sat(hid)))

      gen_auxvar%sat(gid) = 1.d0 - gen_auxvar%sat(lid) - gen_auxvar%sat(hid)
      gen_auxvar%sat(iid) = 0.d0
      
      !if (gen_auxvar%sat(hid) > gen_auxvar%sat(gid)) then
      !  h_sat_eff = gen_auxvar%sat(hid) + gen_auxvar%sat(gid) 
      !  g_sat_eff = 2.d0 * gen_auxvar%sat(gid)
      !else
      !  g_sat_eff = gen_auxvar%sat(hid) + gen_auxvar%sat(gid) 
      !  h_sat_eff = 2.d0 * gen_auxvar%sat(hid)
      !endif
      
      h_sat_eff = gen_auxvar%sat(hid)
      liq_sat_eff = gen_auxvar%sat(lid)/(gen_auxvar%sat(lid)+ &
                    gen_auxvar%sat(gid))
      call HydratePE(gen_auxvar%temp, h_sat_eff, PE_hyd, &
                      characteristic_curves, option)
      call characteristic_curves%saturation_function%CapillaryPressure( &
                liq_sat_eff, gen_auxvar%pres(cpid), &
                dpc_dsatl,option)
      
      gen_auxvar%pres(apid) = PE_hyd

      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                    gen_auxvar%pres(spid),ierr)
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)

      gen_auxvar%pres(spid) = 1.d-6 
      
      gen_auxvar%pres(vpid) = gen_auxvar%pres(spid)
      gen_auxvar%pres(gid) = gen_auxvar%pres(apid) + gen_auxvar%pres(vpid)

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

      if (gen_auxvar%sat(hid) > gen_auxvar%sat(iid)) then
        h_sat_eff = gen_auxvar%sat(hid)+gen_auxvar%sat(iid)
        i_sat_eff = 2.d0 * gen_auxvar%sat(iid)
      else
        h_sat_eff = 2.d0 * gen_auxvar%sat(hid)
        i_sat_eff = gen_auxvar%sat(hid) + gen_auxvar%sat(iid)
      endif

      call GibbsThomsonFreezing(1.d0-i_sat_eff,6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,option)

      gen_auxvar%temp = TQD+dTf
      call HenrysConstantMethane(gen_auxvar%temp,K_H_tilde)
      call HydratePE(gen_auxvar%temp,h_sat_eff, PE_hyd, &
          characteristic_curves, option)
      call EOSWaterSaturationPressure(gen_auxvar%temp, &
                                          gen_auxvar%pres(spid),ierr)

      gen_auxvar%pres(cpid) = 0.d0

      gen_auxvar%pres(spid) = 1.d-6 

      gen_auxvar%pres(apid) = PE_hyd !gen_auxvar%pres(gid)
      gen_auxvar%pres(lid) = gen_auxvar%pres(gid) - gen_auxvar%pres(cpid)
      gen_auxvar%pres(vpid) = 0.d0 !gen_auxvar%pres(gid) - gen_auxvar%pres(apid)

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
     
      if (gen_auxvar%sat(hid) > gen_auxvar%sat(iid)) then
        if (gen_auxvar%sat(hid) > gen_auxvar%sat(gid)) then
          call HydratePE(gen_auxvar%temp, 1.d0, PE_hyd, &
                  characteristic_curves, option)
        else
          call HydratePE(gen_auxvar%temp, 3.d0 * gen_auxvar%sat(iid) + &
                  2.d0 * (gen_auxvar%sat(hid)-gen_auxvar%sat(iid)), PE_hyd, &
                  characteristic_curves, option)
        endif
      elseif (gen_auxvar%sat(hid) > gen_auxvar%sat(gid)) then
        call HydratePE(gen_auxvar%temp, 3.d0 * gen_auxvar%sat(gid) + &
          2.d0 * (gen_auxvar%sat(hid) - gen_auxvar%sat(gid)), PE_hyd, &
          characteristic_curves, option)
      else
        call HydratePE(gen_auxvar%temp, 3.d0 * gen_auxvar%sat(hid), PE_hyd, &
              characteristic_curves, option)
      endif  
      
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
      
      gen_auxvar%sat(gid) = 1.d0 - gen_auxvar%sat(lid) - gen_auxvar%sat(iid)
      gen_auxvar%sat(hid) = 0.d0

      call GibbsThomsonFreezing(1.d0-gen_auxvar%sat(iid),6017.1d0, &
              ICE_DENSITY,TQD,dTf,characteristic_curves,option)

      gen_auxvar%temp = TQD+dTf

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
     
      if (gen_auxvar%sat(hid) > gen_auxvar%sat(iid)) then
        h_sat_eff = gen_auxvar%sat(hid)+gen_auxvar%sat(iid)
        i_sat_eff = 2.d0 * gen_auxvar%sat(iid)
      else
        h_sat_eff = 2.d0 * gen_auxvar%sat(hid)
        i_sat_eff = gen_auxvar%sat(hid) + gen_auxvar%sat(iid)
      endif

      call GibbsThomsonFreezing(1.d0-i_sat_eff,6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,option)

      !if (gen_auxvar%sat(hid) > gen_auxvar%sat(iid)) then
      !  if (gen_auxvar%sat(hid) > gen_auxvar%sat(gid)) then
      !    call HydratePE(gen_auxvar%temp, 1.d0 - gen_auxvar%sat(lid), &
      !            PE_hyd, characteristic_curves, option)
      !  else
      !    call HydratePE(gen_auxvar%temp, 3.d0 * gen_auxvar%sat(iid) + &
      !            2.d0 * (gen_auxvar%sat(hid) - gen_auxvar%sat(iid)), PE_hyd, &
      !            characteristic_curves, option)
      !    call GibbsThomsonFreezing(3.d0 * gen_auxvar%sat(iid), 6017.1d0, &
      !            ICE_DENSITY, TQD, dTf, characteristic_curves, option)
      !  endif
      !elseif (gen_auxvar%sat(hid) > gen_auxvar%sat(gid)) then
      !  call HydratePE(gen_auxvar%temp, 3.d0 * gen_auxvar%sat(gid) + &
      !          2.d0 * (gen_auxvar%sat(hid) - gen_auxvar%sat(gid)), &
      !          PE_hyd, characteristic_curves, option)
      !  call GibbsThomsonFreezing(1.d0-gen_auxvar%sat(lid), 6017.1d0, &
      !          ICE_DENSITY, TQD, dTf, characteristic_curves, option)
      !else
      !  call HydratePE(gen_auxvar%temp, 3.d0 *gen_auxvar%sat(hid), &
      !          PE_hyd, characteristic_curves, option)
      !  if (gen_auxvar%sat(iid) < gen_auxvar%sat(gid)) then
      !    call GibbsThomsonFreezing(3.d0 * gen_auxvar%sat(hid) + &
      !            2.d0 * (gen_auxvar%sat(iid) - gen_auxvar%sat(hid)), &
      !            6017.1d0, ICE_DENSITY, TQD, dTf, characteristic_curves, &
      !            option)
      !  endif
      !endif

      gen_auxvar%temp = TQD + dTf

      call HydratePE(gen_auxvar%temp,h_sat_eff, PE_hyd, &
          characteristic_curves, option)
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
      call PrintErrMsgByRank(option)

  end select

  cell_pressure = max(gen_auxvar%pres(lid),gen_auxvar%pres(gid), &
                      gen_auxvar%pres(spid))

  ! calculate effective porosity as a function of pressure
  if (option%iflag /= GENERAL_UPDATE_FOR_BOUNDARY) then
    dpor_dp = 0.d0
    gen_auxvar%effective_porosity = material_auxvar%porosity_base
    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                gen_auxvar%effective_porosity,dpor_dp)
    endif
    if (option%iflag /= GENERAL_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = gen_auxvar%effective_porosity
    endif
    solid_sat_eff = gen_auxvar%sat(hid) + gen_auxvar%sat(iid)
    
    select case (hydrate_perm_scaling_function)
      case(1) ! Dai and Seol, 2014
        if (hydrate_perm_base(1) < -999.d0) then
          hydrate_perm_base = material_auxvar%permeability
        endif
        material_auxvar%permeability = hydrate_perm_base * &
                    (1.d0-solid_sat_eff)**3/(1.d0+2.d0*solid_sat_eff)**2
      case default
    end select

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

    ! MJ/kmol
    gen_auxvar%U(gid) = xmol_water_in_gas * u_water_vapor + &
                        xmol_air_in_gas * u_air
    Hg_mixture_fractioned = xmol_water_in_gas*h_water_vapor + &
                            xmol_air_in_gas*h_air
    gen_auxvar%H(gid) = gen_auxvar%U(gid) + &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        gen_auxvar%pres(gid)/gen_auxvar%den(gid) * 1.d-6

  endif

  if (gen_auxvar%sat(lid) > 0.d0) then
    if (gen_auxvar%sat(lid) >= 1.d0) then
      krl = 1.d0
    else
      call characteristic_curves%liq_rel_perm_function% &
           RelativePermeability(gen_auxvar%sat(lid),krl,dkrl_dsatl,option)
      krl = max(0.d0,krl)
    endif
    call EOSWaterViscosity(gen_auxvar%temp,cell_pressure, &
                               gen_auxvar%pres(spid),visl,ierr)
    gen_auxvar%mobility(lid) = krl/visl
    gen_auxvar%kr(lid) = krl
  endif

  if (gen_auxvar%sat(gid) > 0.d0) then
    if (gen_auxvar%sat(gid) >=1.d0) then
      krg = 1.d0
    else
      call characteristic_curves%gas_rel_perm_function% &
           RelativePermeability(1.d0-gen_auxvar%sat(gid),krg,dkrg_dsatl,option)
      krg = max(0.d0,krg)
    endif
    call EOSGasViscosity(gen_auxvar%temp,gen_auxvar%pres(apid), &
                           gen_auxvar%pres(gid),den_air,visg,ierr)
    gen_auxvar%mobility(gid) = krg/visg
    gen_auxvar%kr(gid) = krg
  endif

  call EOSHydrateEnthalpy(gen_auxvar%temp, H_hyd)
  gen_auxvar%U(hid) = H_hyd !- cell_pressure/gen_auxvar%den(hid)*1.d-6
  gen_auxvar%H(hid) = H_hyd
  gen_auxvar%mobility(hid) = 0.d0
  
  call EOSIceEnergy(gen_auxvar%temp, U_ice)
  gen_auxvar%xmol(wid,iid) = 1.d0
  gen_auxvar%xmol(gid,iid) = 0.d0
  !call EOSWaterDensityIcePainter(gen_auxvar%temp,gen_auxvar%pres(lid), &
  !                  PETSC_FALSE, gen_auxvar%den(iid), &
  !                  dden_ice_dT, dden_ice_dP, ierr)
  gen_auxvar%den(iid) = ICE_DENSITY
  gen_auxvar%den_kg(iid) = ICE_DENSITY_KG
  gen_auxvar%U(iid) = U_ice
  gen_auxvar%H(iid) = U_ice
  gen_auxvar%mobility(iid) = 0.d0

end subroutine HydrateAuxVarCompute

! ************************************************************************** !

subroutine HydrateAccumulation(gen_auxvar,global_auxvar,material_auxvar, &
                               z,offset,meth,soil_heat_capacity, &
                               option,Res,Jac,analytical_derivatives,debug_cell)
  !
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual, for the hydrate sub-pm
  !
  ! Author: Michael Nole
  ! Date: 03/01/19
  !

  use Option_module
  use Material_Aux_class

  implicit none

  type(general_auxvar_type) :: gen_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: z, offset
  type(methanogenesis_type), pointer :: meth
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_cell

  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase

  PetscReal :: porosity, volume
  PetscReal :: volume_over_dt
  PetscReal :: q_meth

  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  volume_over_dt = material_auxvar%volume / option%flow_dt
  porosity = gen_auxvar%effective_porosity
  volume = material_auxvar%volume

  ! accumulation term units = kmol/s
  Res = 0.d0
  do iphase = 1, option%nphase
    ! Res[kmol comp/m^3 void] = sat[m^3 phase/m^3 void] *
    !                           den[kmol phase/m^3 phase] *
    !                           xmol[kmol comp/kmol phase]
    do icomp = 1, option%nflowspec
      Res(icomp) = Res(icomp) + gen_auxvar%sat(iphase) * &
                                gen_auxvar%den(iphase) * &
                                gen_auxvar%xmol(icomp,iphase)
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
    Res(energy_id) = Res(energy_id) + gen_auxvar%sat(iphase) * &
                                      gen_auxvar%den(iphase) * &
                                      gen_auxvar%U(iphase)
  enddo
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] +
  !                (1-por)[m^3 rock/m^3 bulk] *
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * porosity + &
                    (1.d0 - porosity) * &
                    material_auxvar%soil_particle_density * &
                    soil_heat_capacity * gen_auxvar%temp) * volume_over_dt

  q_meth = 0.d0
  if (associated(meth) .and. offset > 0.d0) then
    call Methanogenesis(z, offset, meth, q_meth)
    !kmol/m^3/s to kmol/s
    q_meth = q_meth*(1.d0 - porosity)*volume
    Res(air_comp_id) = Res(air_comp_id) + q_meth
  endif
end subroutine HydrateAccumulation

! ************************************************************************** !

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
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
!  PetscReal, parameter :: perturbation_tolerance = 1.d-11
  PetscReal, parameter :: min_mole_fraction_pert = 1.d-12
  PetscReal, parameter :: min_perturbation = 1.d-10
  PetscInt :: idof

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

    case(GAS_STATE)
      x(GENERAL_GAS_PRESSURE_DOF) = &
        gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
        gen_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
      x(GENERAL_ENERGY_DOF) = gen_auxvar(ZERO_INTEGER)%temp
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

    case(H_STATE)
      x(GENERAL_GAS_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(GENERAL_GAS_SATURATION_DOF) = MOL_RATIO_METH
      x(GENERAL_ENERGY_DOF) = gen_auxvar(ZERO_INTEGER)%temp

      pert(GENERAL_GAS_PRESSURE_DOF) = -1.d0 * &
         (perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF) + min_perturbation)
      pert(GENERAL_GAS_SATURATION_DOF) = x(GENERAL_GAS_SATURATION_DOF) 
      pert(GENERAL_ENERGY_DOF) = &
         perturbation_tolerance*x(GENERAL_ENERGY_DOF) + min_perturbation

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
      x(GENERAL_GAS_SATURATION_DOF) = &
        gen_auxvar(ZERO_INTEGER)%sat(option%gas_phase)
      x(GENERAL_ENERGY_DOF) = &
        gen_auxvar(ZERO_INTEGER)%temp
       
      pert(GENERAL_GAS_PRESSURE_DOF) = &
        perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)+min_perturbation
      if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
        pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance * &
                                         x(GENERAL_GAS_SATURATION_DOF)
      else
        pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance * &
                                         x(GENERAL_GAS_SATURATION_DOF)
      endif
      pert(GENERAL_ENERGY_DOF) = &
        perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation

    case(HG_STATE)
      x(GENERAL_GAS_PRESSURE_DOF) = &
         gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(option%gas_phase)
      x(GENERAL_ENERGY_DOF) = &
           gen_auxvar(ZERO_INTEGER)%temp
     
      pert(GENERAL_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)+min_perturbation
      if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
        pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation

    case(HA_STATE)
      x(GENERAL_GAS_PRESSURE_DOF) = &
        gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(hid)
      x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%temp

      pert(GENERAL_GAS_PRESSURE_DOF) = &
        perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)+min_perturbation
      if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
        pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(GENERAL_ENERGY_DOF) = &
         perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation

    case(HI_STATE)
      x(GENERAL_GAS_PRESSURE_DOF) = &
        gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(hid)
      x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%temp 

      pert(GENERAL_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)+min_perturbation
      if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
        pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation

    case(GI_STATE)
      x(GENERAL_GAS_PRESSURE_DOF) = &
        gen_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(iid)
      x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%temp

      pert(GENERAL_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(GENERAL_GAS_PRESSURE_DOF)+min_perturbation
      if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
        pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation

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
      x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%temp
    
      pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation
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
      x(GENERAL_GAS_PRESSURE_DOF) = &
        gen_auxvar(ZERO_INTEGER)%sat(hid)
      x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(iid)
      x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%temp

      if (x(GENERAL_GAS_PRESSURE_DOF) > 0.5d0) then
        pert(GENERAL_GAS_PRESSURE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_PRESSURE_DOF) = perturbation_tolerance
      endif
      if (x(GENERAL_GAS_SATURATION_DOF) > 0.5d0) then
        pert(GENERAL_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(GENERAL_ENERGY_DOF) = &
           perturbation_tolerance*x(GENERAL_ENERGY_DOF)+min_perturbation

    case(GAI_STATE)
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

    case(QUAD_STATE)
      x(GENERAL_GAS_PRESSURE_DOF) = &
        gen_auxvar(ZERO_INTEGER)%sat(gid)
      x(GENERAL_GAS_SATURATION_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(lid)
      x(GENERAL_ENERGY_DOF) = &
         gen_auxvar(ZERO_INTEGER)%sat(iid)

      if (x(GENERAL_GAS_PRESSURE_DOF) > 0.5d0) then
        pert(GENERAL_GAS_PRESSURE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(GENERAL_GAS_PRESSURE_DOF) = perturbation_tolerance
      endif
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
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
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
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(HA_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(HI_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(GI_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(AI_STATE)
      gen_auxvar(GENERAL_LIQUID_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_LIQUID_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(HGA_STATE)

    case(HAI_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(HGI_STATE)

    case(GAI_STATE)
      gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert = &
        gen_auxvar(GENERAL_GAS_PRESSURE_DOF)%pert / GENERAL_PRESSURE_SCALE
    case(QUAD_STATE)
  end select


end subroutine HydrateAuxVarPerturb

! ************************************************************************** !

subroutine HydratePE(T,sat, PE, characteristic_curves, option)

  !This subroutine calculates the 3-phase equilibrium pressure of methane
  !hydrate in pure water, from polynomial fit (Moridis, 2003)
  !
  !Author: Michael Nole
  !Date: 01/22/19
  !

  use Characteristic_Curves_module
  use Option_module

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: sat
  PetscReal, intent(out) :: PE

  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option

  PetscReal :: T_temp, dTf

  call GibbsThomsonFreezing(1.d0-sat, 54734.d0, HYDRATE_DENSITY, T, dTf, &
          characteristic_curves, option)

  !MAN: no phase boundary shift
  !dTf = 0.d0
  
  T_temp = T + 273.15d0 + dTf

  if (T < TQD) then
    !Moridis, 2003
    PE = exp(-43.8921173434628 + 0.776302133739303 * T_temp &
          - 7.27291427030502d-3 * T_temp**2 + 3.85413985900724d-5 * T_temp**3 &
          - 1.03669656828834d-7 * T_temp**4 + 1.09882180475307d-10 * T_temp**5)
    !Kamath, 1984
    !PE = exp(1.4717d1-1.88679d3/T_temp)*1.d-3
  else
    !Moridis, 2003
    PE = exp(-1.9413850446456d5 + 3.31018213397926d3 * T_temp &
          - 22.5540264493806* T_temp**2 + 0.0767559117787059 * T_temp**3 &
          - 1.30465829788791d-4 * T_temp**4 + 8.86065316687571d-8 * T_temp**5)
    !Kamath, 1984
    !PE = exp(3.898d1-8.533d3/T_temp)*1.d-3
  endif

  PE = PE * 1.d6

end subroutine HydratePE

! ************************************************************************** !

subroutine HenrysConstantMethane(T,K_H)

  !Calculates the Henry's constant of methane as a function of temperature
  !(Carroll and Mather, 1997)
  !
  !Author: Michael Nole
  !Date: 01/22/19
  !
  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: K_H

  PetscReal :: T_temp
  PetscReal, parameter :: R = 8.314 !J/mol-K

  T_temp = T + 273.15d0

  !Carroll & Mather  Units: Pa/mol frac
  K_H = exp(5.1345 + 7837.d0/T_temp - 1.509d6/(T_temp**2) + 2.06d7/ &
            (T_temp**3)) *1.d3
  
  !Cramer, 1982
  !K_H = 1.d5*(24582.4d0 + 6.71091d2*T + 6.87067d0*T **2 - &
  !            1.773079d-1*T**3 + 1.09652d-03*T**4 - &
  !            3.19599d-6*T**5 + 4.46172d-9*T**6 - &
  !            2.40294d-12*T**7)

end subroutine HenrysConstantMethane

! ************************************************************************** !
subroutine HydrateDavieBuffettCorrection(T,P,K_H)

  implicit none
  
  PetscReal, intent(in) :: T, P
  PetscReal :: K_H

  PetscReal, parameter :: C3_0 = 156.36d0 !mM
  PetscReal, parameter :: T_0 = 292.d0 !K
  PetscReal, parameter :: P_0 = 20.d0 !MPa
  PetscReal, parameter :: dC3_dT = 6.34d0 !mM
  PetscReal, parameter :: dC3_dP = 1.11d0 !mM
  PetscReal, parameter :: alpha = 14.4d0 !C
  PetscReal :: logP

  PetscReal :: T3

  ! Inverting the Moridis equation
  if (T > TQD) then
    !Lower-order
    T3 = 9.0622d0 * log(P*1.d-6) + 264.66d0 
    
    !Higher-order
    !logP = log(P*1.d-6)
    !T3 = -0.11d0*logP**6 + 0.1733d0*logP**5 - 0.9679d0*logP**4 + 2.3492d0*logP**3 &
    !     - 2.7715d0*logP**2 + 11.389d0*logP + 263.5d0
  else
    T3 = T + 273.15d0
  endif

  K_H = K_H / exp((T+273.15d0-T3)/alpha)
  
end subroutine HydrateDavieBuffettCorrection 
! ************************************************************************** !
subroutine EOSHydrateEnthalpy(T,H)

  !Enthalpy of gas hydrate as f(Temperature) (Handa, 1998)
  !
  !Author: Michael Nole
  !Date: 01/22/19
  !
  implicit none

  PetscReal, intent(in):: T
  PetscReal, intent(out) :: H

  PetscReal, parameter :: Hh0 = -54734.d0 ! J/mol
  PetscReal :: Cph, T_temp

  T_temp = T + 273.15d0
  
  ! Integral of Cph * dT ; Cph from Handa, 1998

  ! Units: J/mol
  !H = Hh0 + 6.6d0 * (T_temp-273.15d0) + 7.269d-1 * (T_temp**2 - 273.15d0**2) - 1.21333d-3 * &
  !      (T_temp**3 - 273.15d0**3)  + 1.578d-6 * (T_temp**4 - 273.15d0**4)
  ! Units: MJ/kmol
  !H = H / 1.d3

  !H = H / (Nhyd+1.d0) 

  !Constant Cp
  Cph = 1620.d0*(MW_H20*Nhyd + MW_CH4)/1.d3
  H = Cph * (T-TQD) + Hh0 / (Nhyd + 1.d0)
  H = H / 1.d3

end subroutine EOSHydrateEnthalpy

! ************************************************************************** !

subroutine EOSIceEnergy(T,U)
  
  !Internal energy of ice as f(Temperature) (Fukusako and Yamamoto, 1993)
  !
  !Author: Michael Nole
  !Date: 04/04/19
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: U
  PetscReal, parameter :: Lw = -6017.1d0 !Latent heat of fusion,  J/mol
  PetscReal :: T_temp

  T_temp = T + 273.15d0

  if (T_temp >= 90.d0) then
    U = Lw + 185.d0 * (T_temp-273.15d0) + 3.445 * (T_temp**2 - 273.15d0**2) 
  else
    U = Lw + 4.475 * (T_temp**2 - 273.15d0**2)
  endif

  ! J/mol to MJ/kmol
  U = U / 1.d3
       
end subroutine EOSIceEnergy

! ************************************************************************** !

subroutine GibbsThomsonFreezing(sat,Hf,rho,Tb,dTf,characteristic_curves,option)

  !This subroutine ties the capillary pressure function to a Gibbs-Thomson
  !subcooling required to precipitate a solid in pores.
  !
  !Author: Michael Nole
  !Date: 04/04/19
  !

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
  
  !MAN debugging
  !dTf = 0.d0
end subroutine GibbsThomsonFreezing

! ************************************************************************** !

subroutine HydrateCompositeThermalCond(phi,sat,kdry,kwet,keff)

  implicit none

  PetscReal :: phi, kdry, kwet
  PetscReal, pointer :: sat(:)
  
  PetscReal :: keff
  PetscReal :: k_h20,k_ch4,k_hyd,k_ice
  PetscInt :: lid, gid, hid, iid

  lid = 1
  gid = 2
  hid = 3
  iid = 4
 
  k_h20 = 0.49d0 !W/m-K
  k_ch4 = 30.d-3 !W/m-K
  k_hyd = 0.49d0 !W/m-K
  k_ice = 2.d0   !W/m-K

  keff = kdry + phi * (sat(lid)*kwet + sat(hid)*k_hyd + sat(iid) * k_ice &
          + sat(gid)*k_ch4)

end subroutine HydrateCompositeThermalCond

! ************************************************************************** !

subroutine Methanogenesis(z,offset,meth,q_meth)
  
  ! A simple methanogenesis source parameterized as a function of depth
  ! assuming top of domain is the seafloor
  ! Author: Michael Nole
  ! Date: 03/05/19
  !

  implicit none

  PetscReal :: z, offset
  type(methanogenesis_type), pointer :: meth
  PetscReal :: q_meth

  PetscReal :: alpha, k_alpha, lambda, omega, z_smt
 
  alpha = meth%alpha
  k_alpha = meth%k_alpha
  lambda = meth%lambda
  omega = meth%omega
  z_smt = meth%z_smt 
  
  if (offset - z > z_smt) then
    q_meth = k_alpha * lambda * alpha * exp(-lambda/omega * (offset - &
                                    z - z_smt))
  else
    q_meth = 0.d0
  endif

  !kg/m^3/s to kmol/s
  q_meth = q_meth / MW_CH4

end subroutine Methanogenesis

! ************************************************************************** !

end module Hydrate_module
