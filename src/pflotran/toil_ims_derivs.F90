module TOilIms_derivs_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use AuxVars_TOilIms_module 

  use Global_Aux_module

  use PFLOTRAN_Constants_module
  use PM_TOilIms_Aux_module

  implicit none
  
  private 
#define TOIL_CONVECTION
#define TOIL_CONDUCTION

! Cutoff parameters - no public
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

#if 0
  public :: toil_accum_derivs_alyt, &
            TOilImsFluxPFL_derivs
#endif


  public :: toil_accum_derivs_alyt, &
            TOilImsFluxPFL_derivs, &
            TOilImsSrcSink_derivs, &
            TOilImsBCFlux_derivs

contains


subroutine TOilImsBCFlux_derivs(ibndtype,auxvar_mapping,auxvars, &
                         toil_auxvar_up,global_auxvar_up, &
                         toil_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         sir_dn, &
                         thermal_conductivity_dn, &
                         area,dist,toil_parameter, &
                         option,v_darcy,Res, jdn)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/27/15
  ! 
  use Option_module                              
  use Material_Aux_class
  !use Fracture_module
  !use Klinkenberg_module
  
  implicit none
  
  !type(toil_ims_auxvar_type) :: toil_auxvar_up, toil_auxvar_dn
  class(auxvar_toil_ims_type) :: toil_auxvar_up, toil_auxvar_dn 
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn(:)
  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: v_darcy(option%nphase), area
  type(toil_ims_parameter_type) :: toil_parameter
  PetscReal :: dist(-1:3)
  PetscReal :: Res(1:option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(TOIL_IMS_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)
  !PetscBool :: debug_connection
  
  PetscInt :: energy_id
  PetscInt :: iphase
  PetscInt :: bc_type
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: H_ave, uH
  PetscReal :: perm_dn_adj(option%nphase)
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure, delta_temp
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: sat_dn, perm_dn, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, pres_ave
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  ! for debugging only
  PetscReal :: adv_flux(3,2), diff_flux(2,2)
  PetscReal :: debug_flux(3,3), debug_dphi(2)

  PetscReal :: boundary_pressure
  PetscReal :: tempreal
  
  PetscInt :: idof
  PetscBool :: neumann_bc_present
  
  PetscReal :: dummy_perm_dn

  !PetscReal :: jdn(1:option%nflowdof, 1:option%nflowspec) !! no!
  PetscReal :: jdn(1:option%nflowdof, 1:option%nflowdof)

  PetscReal :: up_scale, dn_scale

  PetscReal :: ddensity_kg_ave_dden_kg_up, ddensity_kg_ave_dden_kg_dn
  PetscReal :: d_delta_pres_dp_up, d_delta_pres_dp_dn
  PetscReal :: d_delta_pres_dT_up, d_delta_pres_dT_dn

  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: d_delta_temp_dt_up, d_delta_temp_dt_dn, dheat_flux_ddelta_temp
  PetscReal :: d_delta_pres_ds_up, d_delta_pres_ds_dn     

  PetscReal, dimension(1:3) :: d_v_darcy_up, d_v_darcy_dn
  PetscReal, dimension(1:3) :: d_q_up, d_q_dn
  PetscReal, dimension(1:3) :: d_mole_flux_up, d_mole_flux_dn
  PetscReal, dimension(1:3) :: d_energy_flux_up, d_energy_flux_dn
  
  energy_id = option%energy_id

  Res = 0.d0
  v_darcy = 0.d0 

  jdn = 0.d0
  d_v_darcy_up = 0.d0
  d_v_darcy_dn = 0.d0
  d_q_up = 0.d0
  d_q_dn = 0.d0
  d_mole_flux_up = 0.d0
  d_mole_flux_dn = 0.d0
  d_energy_flux_up = 0.d0
  d_energy_flux_dn = 0.d0

 
!#ifdef DEBUG_FLUXES    
!  adv_flux = 0.d0
!  diff_flux = 0.d0
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  debug_flux = 0.d0
!  debug_dphi = 0.d0
!#endif

  neumann_bc_present = PETSC_FALSE
  
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  ! currently no fractures considered 
  ! Fracture permeability change only available for structured grid (Heeho)
  !if (associated(material_auxvar_dn%fracture)) then
  !  call FracturePermEvaluate(material_auxvar_dn,perm_dn,perm_dn, &
  !                            dummy_perm_dn,dist)
  !endif  
  
  !if (associated(klinkenberg)) then
  !  perm_dn_adj(1) = perm_dn
  !                                        
  !  perm_dn_adj(2) = klinkenberg%Evaluate(perm_dn, &
  !                                        gen_auxvar_dn%pres(option%gas_phase))
  !else
    perm_dn_adj(:) = perm_dn
  !endif
  
#ifdef TOIL_CONVECTION  
  do iphase = 1, option%nphase
    d_v_darcy_up = 0.d0
    d_v_darcy_dn = 0.d0
    d_q_up = 0.d0
    d_q_dn = 0.d0
    d_mole_flux_up = 0.d0
    d_mole_flux_dn = 0.d0
    d_energy_flux_up = 0.d0
    d_energy_flux_dn = 0.d0
   
    bc_type = ibndtype(iphase) ! loop over equations 1.Liq and 2.Oil
    select case(bc_type)
      ! figure out the direction of flow
      case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == CONDUCTANCE_BC) then !not implemented yet
          select case(iphase)
            case(LIQUID_PHASE)
              idof = auxvar_mapping(TOIL_IMS_LIQ_CONDUCTANCE_INDEX)
            case(TOIL_IMS_OIL_PHASE)
              idof = auxvar_mapping(TOIL_IMS_OIL_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
        ! PO need to check what values of saturations are assigned to the BC ghost cells  
        ! using residual saturation cannot be correct! - geh
        ! reusing sir_dn for bounary auxvar
!#define BAD_MOVE1 ! this works
!#ifndef BAD_MOVE1       
        if (toil_auxvar_up%sat(iphase) > sir_dn(iphase) .or. &
            toil_auxvar_dn%sat(iphase) > sir_dn(iphase)) then
!#endif
          boundary_pressure = toil_auxvar_up%pres(iphase)

          ! PO no free surfce boundaries considered  
          !if (iphase == LIQUID_PHASE .and. &
          !    global_auxvar_up%istate == GAS_STATE) then
          !  ! the idea here is to accommodate a free surface boundary
          !  ! face.  this will not work for an interior grid cell as
          !  ! there should be capillary pressure in force.
          !  boundary_pressure = gen_auxvar_up%pres(option%gas_phase)
          !endif

          !density_kg_ave = 0.5d0 * (toil_auxvar_up%den_kg(iphase) + &
          !                          toil_auxvar_dn%den_kg(iphase) )

#if 0
          density_kg_ave = TOilImsAverageDensity(toil_auxvar_up%sat(iphase), &
                           toil_auxvar_dn%sat(iphase), &
                           toil_auxvar_up%den_kg(iphase), &
                           toil_auxvar_dn%den_kg(iphase))
#endif


          density_kg_ave = TOilImsAverageDensity_derivs(toil_auxvar_up%sat(iphase), &
                           toil_auxvar_dn%sat(iphase), &
                           toil_auxvar_up%den_kg(iphase), &
                           toil_auxvar_dn%den_kg(iphase), &
                           ddensity_kg_ave_dden_kg_up, &
                           ddensity_kg_ave_dden_kg_dn)

          gravity_term = density_kg_ave * dist_gravity
          delta_pressure = boundary_pressure - &
                           toil_auxvar_dn%pres(iphase) + &
                           gravity_term
          

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!          debug_dphi(iphase) = delta_pressure
!#endif
          ! PO CONDUCTANCE_BC and SEEPAGE_BC not implemented
          if (bc_type == SEEPAGE_BC .or. &
              bc_type == CONDUCTANCE_BC) then
                ! flow in         ! boundary cell is <= pref
            if (delta_pressure > 0.d0 .and. &
                toil_auxvar_up%pres(iphase) - &
                 option%reference_pressure < eps) then
              delta_pressure = 0.d0
            endif
          endif
          
          !upwinding mobilities and enthalpies   
          up_scale = 0.d0 !! ADDED
          dn_scale = 0.d0 !! ADDED
          if (delta_pressure >= 0.D0) then
            up_scale = 1.d0
            mobility = toil_auxvar_up%mobility(iphase)
            uH = toil_auxvar_up%H(iphase)
#ifdef TOIL_DEN_UPWIND
            density_ave = toil_auxvar_up%den(iphase)
#endif
          else
            dn_scale = 1.d0
            mobility = toil_auxvar_dn%mobility(iphase)
            uH = toil_auxvar_dn%H(iphase)
#ifdef TOIL_DEN_UPWIND
            density_ave = toil_auxvar_dn%den(iphase)
#endif
          endif      

          if (mobility > floweps) then
            ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
            !                    dP[Pa]]
            v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure
            ! only need average density if velocity > 0.

            ! when this is commented - using upwinding value
            !density_ave = 0.5d0 * (toil_auxvar_up%den(iphase) + &
            !                       toil_auxvar_dn%den(iphase) )
#ifndef TOIL_DEN_UPWIND
#if 0
            density_ave = TOilImsAverageDensity(toil_auxvar_up%sat(iphase), &
                           toil_auxvar_dn%sat(iphase), &
                           toil_auxvar_up%den(iphase), &
                           toil_auxvar_dn%den(iphase))
#endif

            density_ave = TOilImsAverageDensity_derivs(toil_auxvar_up%sat(iphase), &
                                 toil_auxvar_dn%sat(iphase), &
                                 toil_auxvar_up%den(iphase), &
                                 toil_auxvar_dn%den(iphase), &
                                 ddensity_ave_dden_up, &
                                 ddensity_ave_dden_dn)
            

      !! the upwind values computed by this are wrong but will not be used
#if 0
      call  DeltaPressureDerivs_up_and_down(d_delta_pres_dp_up, d_delta_pres_dp_dn,     &
                                            d_delta_pres_dT_up, d_delta_pres_dT_dn,     &
                                            dist_gravity,                               &
                                            ddensity_ave_dden_up, ddensity_ave_dden_dn, &
                                            toil_auxvar_up%d%dden_dp(iphase,1),                &
                                            toil_auxvar_dn%d%dden_dp(iphase,1),                &
                                            toil_auxvar_up%d%dden_dt(iphase),                 &
                                            toil_auxvar_dn%d%dden_dt(iphase)                 )
#endif

      call  DeltaPressureDerivs_up_and_down(d_delta_pres_dp_up, d_delta_pres_dp_dn,     &
                                            d_delta_pres_dT_up, d_delta_pres_dT_dn,     &
                                            dist_gravity,                               &
                                            ddensity_ave_dden_up, ddensity_ave_dden_dn, &
                                            toil_auxvar_up%d%dden_dp(iphase,1),                &
                                            toil_auxvar_dn%d%dden_dp(iphase,1),                &
                                            toil_auxvar_up%d%dden_dt(iphase),                 &
                                            toil_auxvar_dn%d%dden_dt(iphase),                 &
                                            toil_auxvar_up%d%dp_dsat(iphase),                 &
                                            toil_auxvar_dn%d%dp_dsat(iphase),                &
                                            d_delta_pres_ds_up, d_delta_pres_ds_dn,    &
                                            toil_ims_fmw_comp(iphase) )


      call v_darcy_derivs(d_v_darcy_dn, toil_auxvar_dn%d%dmobility(iphase, 1:3),     &
                          dn_scale, delta_pressure, mobility, perm_ave_over_dist,  &
                          d_delta_pres_dp_dn, d_delta_pres_dt_dn,            &
                          d_delta_pres_ds_dn)
                          

#endif
          endif
!#ifndef BAD_MOVE1        
        endif ! sat > eps
!#endif

      case(NEUMANN_BC)
        select case(iphase)
          case(LIQUID_PHASE)
            idof = auxvar_mapping(TOIL_IMS_LIQUID_FLUX_INDEX)
          case(TOIL_IMS_OIL_PHASE)
            idof = auxvar_mapping(TOIL_IMS_OIL_FLUX_INDEX)
        end select
        
        neumann_bc_present = PETSC_TRUE
        if (dabs(auxvars(idof)) > floweps) then
          up_scale = 0.d0 !! ADDED
          dn_scale = 0.d0 !! ADDED

          v_darcy(iphase) = auxvars(idof)
          !!! so d_v_darcy is 0 here

          if (v_darcy(iphase) > 0.d0) then 
            up_scale = 1.d0
            density_ave = toil_auxvar_up%den(iphase)
            uH = toil_auxvar_up%H(iphase)
          else 
            dn_scale = 1.d0
            density_ave = toil_auxvar_dn%den(iphase)
            uH = toil_auxvar_dn%H(iphase)
          endif 
        endif
      case default
        option%io_buffer = &
          'Boundary condition type not recognized in GeneralBCFlux phase loop.'
        call printErrMsg(option)
    end select

    !if (dabs(v_darcy(iphase)) > 0.d0) then !!! is this right? will get a lot of zero derivatives otherwise
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area

      d_q_up = d_v_darcy_up * area
      d_q_dn = d_v_darcy_dn * area

      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                              density_ave[kmol phase/m^3 phase]
      mole_flux = q*density_ave       

      call MoleFluxDerivs(d_mole_flux_dn, d_q_dn, q, density_ave, ddensity_ave_dden_dn, &
                          toil_auxvar_dn%d%dden_dp(iphase,1), toil_auxvar_dn%d%dden_dt(iphase))

      jdn(iphase, 1:3) = jdn(iphase, 1:3) + d_mole_flux_dn

  
      ! Res[kmol total/sec]
      Res(iphase) = Res(iphase) + mole_flux
 
      ! Res[kmol total/sec]
      !do icomp = 1, option%nflowspec
      !  ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
      !  !                      xmol[kmol comp/mol phase]
      !  Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      !enddo
!#ifdef DEBUG_FLUXES  
!      do icomp = 1, option%nflowspec
!        adv_flux(icomp,iphase) = adv_flux(icomp,iphase) + mole_flux * xmol(icomp)
!      enddo
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      do icomp = 1, option%nflowspec
!        debug_flux(icomp,iphase) = debug_flux(icomp,iphase) + mole_flux * xmol(icomp)
!      enddo
!#endif
      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
      Res(energy_id) = Res(energy_id) + mole_flux * uH ! H_ave

      !! by `energy flux' mean the term moleFlux * uH
      call EnergyFluxDerivs(d_energy_flux_dn, d_mole_flux_dn, uH, dn_scale, mole_flux, &
                            toil_auxvar_dn%d%dH_dp(iphase), toil_auxvar_dn%d%dH_dt(iphase))
      jdn(energy_id, 1:3) = jdn(energy_id, 1:3) + d_energy_flux_dn

!#ifdef DEBUG_FLUXES  
!      adv_flux(energy_id,iphase) = adv_flux(energy_id,iphase) + mole_flux * uH
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      debug_flux(energy_id,iphase) = debug_flux(energy_id,iphase) + mole_flux * uH
!#endif
    !endif
  enddo
#endif 
! end of TOIL_CONVECTION
  
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  if (debug_flag > 0) then 
!    write(debug_unit,'(a,7es24.15)') 'bc delta pressure :', debug_dphi(:)  
!    write(debug_unit,'(a,7es24.15)') 'bc adv flux (liquid):', debug_flux(:,1)
!    write(debug_unit,'(a,7es24.15)') 'bc adv flux (gas):', debug_flux(:,2)
!  endif
!  debug_flux = 0.d0
!#endif  


#ifdef TOIL_CONDUCTION
  ! add heat conduction flux
  heat_flux = 0.d0
  select case (ibndtype(TOIL_IMS_ENERGY_EQUATION_INDEX))
    case (DIRICHLET_BC)
      ! based on Somerton et al., 1974:
      ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
      !k_eff_dn = thermal_conductivity_dn(1) + &
      !           sqrt(gen_auxvar_dn%sat(option%liquid_phase)) * &
      !           (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
      ! considered the formation fully saturated in water for heat conduction
      k_eff_dn = thermal_conductivity_dn(1)
      ! units:
      ! k_eff = W/K/m/m = J/s/K/m/m
      ! delta_temp = K
      ! area = m^2
      ! heat_flux = J/s
      k_eff_ave = k_eff_dn / dist(0)
      delta_temp = toil_auxvar_up%temp - toil_auxvar_dn%temp

      d_delta_temp_dt_dn = - 1.d0

      heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! convert W -> MW

      dheat_flux_ddelta_temp = k_eff_ave * area * 1.d-6 ! J/s -> MJ/s
    case(NEUMANN_BC)
                  ! flux prescribed as MW/m^2
      heat_flux = auxvars(auxvar_mapping(TOIL_IMS_ENERGY_FLUX_INDEX)) * area

      dheat_flux_ddelta_temp = 0.d0 ! constant

    case(ZERO_GRADIENT_BC)
      ! No contribution to heat_flux

      dheat_flux_ddelta_temp = 0.d0 ! constant

    case default
      option%io_buffer = 'Boundary condition type not recognized in ' // &
        'TOilImsBCFlux heat conduction loop.'
      call printErrMsg(option)
  end select
  Res(energy_id) = Res(energy_id) + heat_flux ! MW

  jdn(energy_id, 3) = jdn(energy_id, 3) + d_delta_temp_dt_dn*dheat_flux_ddelta_temp
#endif 
! end of TOIL_CONDUCTION


end subroutine TOilImsBCFlux_derivs

! ************************************************************************** !

subroutine TOilImsSrcSink_derivs(option,src_sink_condition, toil_auxvar, &
                          global_auxvar,ss_flow_vol_flux,scale,Res, j)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/04/15
  ! 

  use Option_module
  use Condition_module  

  use EOS_Water_module
  use EOS_Oil_module

  implicit none

  type(option_type) :: option
  type(flow_toil_ims_condition_type), pointer :: src_sink_condition
  !type(toil_ims_auxvar_type) :: toil_auxvar
  class(auxvar_toil_ims_type) :: toil_auxvar
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
  PetscReal :: den, den_kg, enthalpy, internal_energy, temperature
  PetscReal :: cell_pressure, dummy_pressure
  PetscInt :: iphase
  PetscInt :: energy_var
  PetscErrorCode :: ierr

  PetscReal :: dden_bool, denth_bool
  PetscReal :: hw_dp, hw_dT, ho_dp, ho_dT
  PetscReal, dimension(1:3,1:3) :: j
  PetscReal, dimension(1:3) :: d_qsrc_mol, d_inj_en_part
  PetscReal :: scale_use

  j = 0.d0

  cell_pressure = 0.d0

  ! this can be removed when etxending to pressure condition
  if (.not.associated(src_sink_condition%rate) ) then
    option%io_buffer = 'TOilImsSrcSink fow condition rate not defined ' // &
    'rate is needed for a valid src/sink term'
    call printErrMsg(option)  
  end if

  !qsrc => src_sink_condition%rate%dataset%rarray(:)
  qsrc => src_sink_condition%rate%dataset%rarray

  energy_var = 0
  if ( associated(src_sink_condition%temperature) ) then
    energy_var = SRC_TEMPERATURE 
  else if ( associated(src_sink_condition%enthalpy) ) then
    energy_var = SRC_ENTHALPY
  end if

  flow_src_sink_type = src_sink_condition%rate%itype

 ! checks that qsrc(liquid_phase) and qsrc(oil_phase) 
 ! do not have different signs
  if ( (qsrc(option%liquid_phase)>0.0d0 .and. qsrc(option%oil_phase)<0.d0).or.&
      (qsrc(option%liquid_phase)<0.0d0 .and. qsrc(option%oil_phase)>0.d0)  & 
    ) then
    option%io_buffer = "TOilImsSrcSink error: " // &
      "src(wat) and src(oil) with opposite sign"
    call printErrMsg(option)
  end if

  ! approximates BHP with local pressure
  ! to compute BHP we need to solve an IPR equation
    cell_pressure = &
        maxval(toil_auxvar%pres(option%liquid_phase:option%oil_phase))

  ! if enthalpy is used to define enthelpy or energy rate is used  
  ! approximate bottom hole temperature (BHT) with local temp
  if ( energy_var == SRC_TEMPERATURE) then
    temperature = src_sink_condition%temperature%dataset%rarray(1)
  else   
    temperature = toil_auxvar%temp
  end if


  Res = 0.d0
  dden_bool = 0.d0
  do iphase = 1, option%nphase
    qsrc_mol = 0.d0
    if ( qsrc(iphase) > 0.d0) then 
      select case(iphase)
        case(LIQUID_PHASE)
          call EOSWaterDensity(temperature,cell_pressure,den_kg,den,ierr)
        case(TOIL_IMS_OIL_PHASE)
            call EOSOilDensity(temperature,cell_pressure,den,ierr)
      end select 
    else
      den = toil_auxvar%den(iphase)
    end if

    scale_use = 0.d0
    select case(flow_src_sink_type)
      ! injection and production 
      case(MASS_RATE_SS)
        qsrc_mol = qsrc(iphase)/toil_ims_fmw_comp(iphase) ! kg/sec -> kmol/sec
      case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
        qsrc_mol = qsrc(iphase)/toil_ims_fmw_comp(iphase)*scale 
        scale_use = scale
      case(VOLUMETRIC_RATE_SS)  ! assume local density for now 
                  ! qsrc(iphase) = m^3/sec  
        qsrc_mol = qsrc(iphase)*den ! den = kmol/m^3 
        dden_bool = 1.d0
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec             ! den = kmol/m^3
        qsrc_mol = qsrc(iphase)* den * scale
        dden_bool = 1.d0
        !qsrc_mol = qsrc(iphase)*gen_auxvar%den(iphase)*scale 
        scale_use = scale
    end select
    ss_flow_vol_flux(iphase) = qsrc_mol/ den
    Res(iphase) = qsrc_mol
    call Qsrc_mol_derivs(d_qsrc_mol,dden_bool, qsrc(iphase), toil_auxvar%d%dden_dp(iphase,1), &
                         toil_auxvar%d%dden_dt(iphase), scale_use)
    j(iphase, 1:3) = j(iphase, 1:3) + d_qsrc_mol(1:3)
  enddo

  ! when using scaled src/sinks, the rates (marr or vol) scaling 
  ! at this point the scale factor is already included in Res(iphase)

  ! Res(option%energy_id), energy units: MJ/sec

  if ( associated(src_sink_condition%temperature) .or. &
      associated(src_sink_condition%enthalpy) &
     ) then
    ! if injection compute local pressure that will be used as BHP
    ! approximation used to overcome the solution of an IPR
    !if ( qsrc(option%liquid_phase)>0.d0 .or. 
    !    qsrc(option%oil_phase)>0.d0 ) then
    !  cell_pressure = &
    !      maxval(toil_auxvar%pres(option%liquid_phase:option%oil_phase))
    !end if
    ! water injection 
    if (qsrc(option%liquid_phase) > 0.d0) then !implies qsrc(option%oil_phase)>=0
      denth_bool = 0.d0
      if ( energy_var == SRC_TEMPERATURE ) then
        call EOSWaterDensity(src_sink_condition%temperature% &
                             dataset%rarray(1), cell_pressure, &
                             den_kg,den,ierr)
        !!! changed to include derivatives
        call EOSWaterEnthalpy(src_sink_condition%temperature% &
                              dataset%rarray(1), cell_pressure, &
                              enthalpy, hw_dp, hw_dT,ierr)
          hw_dp = hw_dp * 1.d-6
          !hw_dT = hw_dT * 1.d-6
          hw_dT = 0.d0 !! because input temp is constant here
          denth_bool = 1.d0

        !call EOSWaterEnthalpy(gen_auxvar%temp,cell_pressure,enthalpy, &
                              !hw_dp,hw_dT,ierr)
        ! enthalpy = [J/kmol]
      else if ( energy_var == SRC_ENTHALPY ) then
        !input as J/kg
        enthalpy = src_sink_condition%enthalpy% &
                       dataset%rarray(option%liquid_phase)
                     ! J/kg * kg/kmol = J/kmol  
        enthalpy = enthalpy * toil_ims_fmw_comp(option%liquid_phase) 
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      ! enthalpy units: MJ/kmol ! water component mass                     
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%liquid_phase) * enthalpy

       !!! LIKELY NEEDS MORE - derivatives of Res(...) are part of J already
       !!! and should be used
      call InjectionEnergyPartDerivs(d_inj_en_part, denth_bool, &
                                     Res(option%liquid_phase), enthalpy, &
                                     hw_dp, hw_dT)
      j(option%energy_id, :) = j(option%energy_id, :) +  d_inj_en_part

    end if
    ! oil injection 
    if (qsrc(option%oil_phase) > 0.d0) then !implies qsrc(option%liquid_phase)>=0
      denth_bool = 0.d0
      if ( energy_var == SRC_TEMPERATURE ) then
        call EOSOilEnthalpy(src_sink_condition%temperature%dataset%rarray(1), &
                            cell_pressure, enthalpy, hw_dp, hw_dT,  ierr)
          ho_dp = ho_dp * 1.d-6
          !ho_dT = ho_dT * 1.d-6
          ho_dT = 0.d0 !! because input temp is constant here
          denth_bool = 1.d0

        ! enthalpy = [J/kmol] 
      else if ( energy_var == SRC_ENTHALPY ) then
        enthalpy = src_sink_condition%enthalpy% &
                     dataset%rarray(option%oil_phase)
                      !J/kg * kg/kmol = J/kmol  
        enthalpy = enthalpy * toil_ims_fmw_comp(option%oil_phase)        
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      ! enthalpy units: MJ/kmol ! oil component mass                     
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%oil_phase) * enthalpy

       !!! LIKELY NEEDS MORE - derivatives of Res(...) are part of J already
       !!! and should be used
      call InjectionEnergyPartDerivs(d_inj_en_part, denth_bool, &
                                     Res(option%oil_phase), enthalpy, &
                                     ho_dp, ho_dT)
      j(option%energy_id, :) = j(option%energy_id, :) +  d_inj_en_part

    end if
    ! water energy extraction due to water production
    if (qsrc(option%liquid_phase) < 0.d0) then !implies qsrc(option%oil_phase)<=0
      ! auxvar enthalpy units: MJ/kmol ! water component mass                     
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%liquid_phase) * &
                              toil_auxvar%H(option%liquid_phase)

       !!! LIKELY NEEDS MORE - derivatives of Res(...) are part of J already
       !!! and should be used

      !! bit of a misuse here since commenting in this routine says it's
      !! for r * enthalpy
      !! here use it for 
      !! r * H
      !! but it's fine if the arguments are changed appropriately.
      !! Really this is a general purpose simple derivative routine
      !! for const * x 
      !! where only dx_dp and dx_dT are nonzero
      !! and switched on or off by an upstreaming-like flag
      call InjectionEnergyPartDerivs(d_inj_en_part, denth_bool, &
                                     Res(option%liquid_phase), &
                                     toil_auxvar%H(option%liquid_phase) , &
                                     toil_auxvar%d%dH_dp(option%liquid_phase), &
                                     toil_auxvar%d%dH_dT(option%liquid_phase))
      j(option%energy_id, :) = j(option%energy_id, :) +  d_inj_en_part

    end if
    !oil energy extraction due to oil production 
    if (qsrc(option%oil_phase) < 0.d0) then !implies qsrc(option%liquid_phase)<=0
      ! auxvar enthalpy units: MJ/kmol ! water component mass                     
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%oil_phase) * &
                              toil_auxvar%H(option%oil_phase)

       !!! LIKELY NEEDS MORE - derivatives of Res(...) are part of J already
       !!! and should be used

      call InjectionEnergyPartDerivs(d_inj_en_part, denth_bool, &
                                     Res(option%oil_phase), &
                                     toil_auxvar%H(option%oil_phase) , &
                                     toil_auxvar%d%dH_dp(option%oil_phase), &
                                     toil_auxvar%d%dH_dT(option%oil_phase))
      j(option%energy_id, :) = j(option%energy_id, :) +  d_inj_en_part

    end if

  else !if not temp or enthalpy are given
    ! if energy rate is given, loaded in qsrc(3) in MJ/sec 
    Res(option%energy_id) = qsrc(THREE_INTEGER)* scale ! MJ/s

    !! no derivatives for this!
  end if


  nullify(qsrc)      
  
end subroutine TOilImsSrcSink_derivs
! ************************************************************************** !

!subroutine InjectionEnergyPartDerivs(d_inj_en_part, denth_bool, &
                               !r, enthalpy, &
                               !hw_dp, hw_dT)

subroutine InjectionEnergyPartDerivs(d_inj_en_part, denth_bool, r, &
                                      enthalpy, hw_dp, hw_dT)

  implicit none
  PetscReal, dimension(1:3) :: d_inj_en_part
  PetscReal :: denth_bool, r, enthalpy, hw_dp, hw_dT
     
  !r *  enthalpy

  d_inj_en_part = 0.d0

  !! w.r.t. oil pressure
  d_inj_en_part(1) = denth_bool * r * hw_dp


  !! w.r.t. temperature
  d_inj_en_part(3) = denth_bool * r * hw_dT

 

end subroutine InjectionEnergyPartDerivs

! ************************************************************************** !

subroutine Qsrc_mol_derivs(d_qsrc_mol,dden_bool, qsrc, dden_dp, dden_dt, sc)
  
  implicit none

  PetscReal, dimension(1:3) :: d_qsrc_mol
  PetscReal :: dden_bool
  PetscReal :: qsrc, dden_dp, dden_dt, sc

  !! qsrc_mol = sc*qsrc*den

  d_qsrc_mol = 0.d0

  !! w.r.t. pressure:
  d_qsrc_mol(1) = dden_bool * qsrc * dden_dp

  !! w.r.t. saturation
  !!  (is 0)

  !! w.r.t. temperature:
  d_qsrc_mol(3) = dden_bool * qsrc * dden_dT

  d_qsrc_mol = d_qsrc_mol * sc


end subroutine Qsrc_mol_derivs

! ************************************************************************** !

subroutine TOilImsFluxPFL_derivs(toil_auxvar_up,global_auxvar_up, &
                       material_auxvar_up, &
                       sir_up, &
                       thermal_conductivity_up, &
                       toil_auxvar_dn,global_auxvar_dn, &
                       material_auxvar_dn, &
                       sir_dn, &
                       thermal_conductivity_dn, &
                       area, dist, &
                       option,v_darcy,Res, &
                       jup, jdn)

  use Option_module
  use Material_Aux_class
  use Connection_module
 
  ! no fractures considered for now
  ! use Fracture_module
  !use Klinkenberg_module
  
  implicit none
  
  !type(toil_ims_auxvar_type) :: toil_auxvar_up, toil_auxvar_dn
  PetscReal, dimension(1:3,1:3) :: jup, jdn
  class(auxvar_toil_ims_type) :: toil_auxvar_up, toil_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  !type(toil_ims_parameter_type) :: parameter
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: Res(option%nflowdof)
  !PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscInt :: energy_id
  PetscInt :: iphase
 
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn           ! no mole fractions
  PetscReal :: delta_pressure, delta_temp !, delta_xmol,

  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: stpd_up, stpd_dn
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, tempreal
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux

  ! no diff fluxes - arrays used for debugging only
  PetscReal :: adv_flux(3,2), diff_flux(2,2)
  PetscReal :: debug_flux(3,3), debug_dphi(2)
  
  PetscReal :: dummy_perm_up, dummy_perm_dn

  PetscReal :: up_scale, dn_scale

  PetscReal :: ddensity_kg_ave_dden_kg_up, ddensity_kg_ave_dden_kg_dn
  PetscReal :: d_delta_pres_dp_up, d_delta_pres_dp_dn
  PetscReal :: d_delta_pres_dT_up, d_delta_pres_dT_dn

  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: d_delta_temp_dt_up, d_delta_temp_dt_dn, dheat_flux_ddelta_temp
  PetscReal :: d_delta_pres_ds_up, d_delta_pres_ds_dn   

  PetscReal, dimension(1:3) :: d_v_darcy_up, d_v_darcy_dn
  PetscReal, dimension(1:3) :: d_q_up, d_q_dn
  PetscReal, dimension(1:3) :: d_mole_flux_up, d_mole_flux_dn
  PetscReal, dimension(1:3) :: d_energy_flux_up, d_energy_flux_dn



  energy_id = option%energy_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  
    perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
      
  Res = 0.d0
  
  v_darcy = 0.d0

  jup = 0.d0
  jdn = 0.d0


#ifdef TOIL_CONVECTION
  do iphase = 1, option%nphase
    !print *, "iphase: ", iphase, ", phases: ", option%nphase
    d_v_darcy_up = 0.d0
    d_v_darcy_dn = 0.d0
    d_q_up = 0.d0
    d_q_dn = 0.d0
    d_mole_flux_up = 0.d0
    d_mole_flux_dn = 0.d0
    d_energy_flux_up = 0.d0
    d_energy_flux_dn = 0.d0
 
    if (toil_auxvar_up%mobility(iphase) + &
        toil_auxvar_dn%mobility(iphase) < eps) then
      cycle
    endif

    ! an alternative could be to avergae using oil_sat
    !density_kg_ave = 0.5d0* ( toil_auxvar_up%den_kg(iphase) + &
    !                          toil_auxvar_dn%den_kg(iphase) )
    density_kg_ave = TOilImsAverageDensity_derivs(toil_auxvar_up%sat(iphase), &
                     toil_auxvar_dn%sat(iphase), &
                     toil_auxvar_up%den_kg(iphase), &
                     toil_auxvar_dn%den_kg(iphase), &
                     ddensity_kg_ave_dden_kg_up, &
                     ddensity_kg_ave_dden_kg_dn)


    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = toil_auxvar_up%pres(iphase) - &
                     toil_auxvar_dn%pres(iphase) + &
                     gravity_term



!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      debug_dphi(iphase) = delta_pressure
!#endif

    ! upwinding the mobilities and enthalpies
    up_scale = 0.d0 !! ADDED
    dn_scale = 0.d0 !! ADDED
    if (delta_pressure >= 0.D0) then
      up_scale = 1.d0
      mobility = toil_auxvar_up%mobility(iphase)
      H_ave = toil_auxvar_up%H(iphase)
      uH = H_ave
#ifdef TOIL_DEN_UPWIND
      density_ave = toil_auxvar_up%den(iphase)
      ddensity_ave_dden_up = 1.d0
      ddensity_ave_dden_dn = 0.d0
#endif
    else
      dn_scale = 1.d0
      mobility = toil_auxvar_dn%mobility(iphase)
      H_ave = toil_auxvar_dn%H(iphase)
      uH = H_ave
#ifdef TOIL_DEN_UPWIND
      density_ave = toil_auxvar_dn%den(iphase)
      ddensity_ave_dden_up = 0.d0
      ddensity_ave_dden_dn = 1.d0
#endif
    endif      

    if (mobility > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure



      ! if comments below, use upwinding value
      !density_ave = 0.5d0*( toil_auxvar_up%den(iphase) + &
      !                      toil_auxvar_dn%den(iphase))
#ifndef TOIL_DEN_UPWIND
      density_ave = TOilImsAverageDensity_derivs(toil_auxvar_up%sat(iphase), &
                           toil_auxvar_dn%sat(iphase), &
                           toil_auxvar_up%den(iphase), &
                           toil_auxvar_dn%den(iphase), &
                           ddensity_ave_dden_up, &
                           ddensity_ave_dden_dn)
#endif 

      !! defer delta pressure derivtives to here because we know density average 
      !! derivatives will have been calculated by this point
      call  DeltaPressureDerivs_up_and_down(d_delta_pres_dp_up, d_delta_pres_dp_dn,     &
                                            d_delta_pres_dT_up, d_delta_pres_dT_dn,     &
                                            dist_gravity,                               &
                                            ddensity_ave_dden_up, ddensity_ave_dden_dn, &
                                            toil_auxvar_up%d%dden_dp(iphase,1),                &
                                            toil_auxvar_dn%d%dden_dp(iphase,1),                &
                                            toil_auxvar_up%d%dden_dt(iphase),                 &
                                            toil_auxvar_dn%d%dden_dt(iphase),                 &
                                            toil_auxvar_up%d%dp_dsat(iphase),                 &
                                            toil_auxvar_dn%d%dp_dsat(iphase),                 &
                                            d_delta_pres_ds_up, d_delta_pres_ds_dn,&
                                            toil_ims_fmw_comp(iphase))

                                           !dp_ds_up, dp_ds_dn,                         &
                                           !ddelta_pressure_dsatup, ddelta_pressure_dsatdn )


      call v_darcy_derivs(d_v_darcy_up, toil_auxvar_up%d%dmobility(iphase, 1:3),     &
                          up_scale, delta_pressure, mobility, perm_ave_over_dist(iphase),  &
                          d_delta_pres_dp_up, d_delta_pres_dt_up, d_delta_pres_ds_up)
      call v_darcy_derivs(d_v_darcy_dn, toil_auxvar_dn%d%dmobility(iphase, 1:3),     &
                          dn_scale, delta_pressure, mobility, perm_ave_over_dist(iphase),  &
                          d_delta_pres_dp_dn, d_delta_pres_dt_dn, d_delta_pres_ds_dn)

      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  

      d_q_up = d_v_darcy_up * area
      d_q_dn = d_v_darcy_dn * area

      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      mole_flux = q*density_ave

      call MoleFluxDerivs(d_mole_flux_up, d_q_up, q, density_ave, ddensity_ave_dden_up, &
                          toil_auxvar_up%d%dden_dp(iphase,1), toil_auxvar_up%d%dden_dt(iphase))
      call MoleFluxDerivs(d_mole_flux_dn, d_q_dn, q, density_ave, ddensity_ave_dden_dn, &
                          toil_auxvar_dn%d%dden_dp(iphase,1), toil_auxvar_dn%d%dden_dt(iphase))

      ! Res[kmol total/sec]

      ! Res[kmol phase/sec] = mole_flux[kmol phase/sec]  
      Res(iphase) = Res(iphase) + mole_flux 
      jup(iphase,1:3) = jup(iphase, 1:3) + d_mole_flux_up
      jdn(iphase, 1:3) = jdn(iphase, 1:3) + d_mole_flux_dn

      !do icomp = 1, option%nflowspec
      !  ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
      !  !                      xmol[kmol comp/kmol phase]
      !  Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      !enddo

!#ifdef DEBUG_FLUXES  
!      do icomp = 1, option%nflowspec
!        adv_flux(icomp) = adv_flux(icomp) + mole_flux * xmol(icomp)
!      enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      do icomp = 1, option%nflowspec
!        debug_flux(icomp,iphase) = debug_flux(icomp,iphase) + mole_flux * xmol(icomp)
!      enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
!#endif

      Res(energy_id) = Res(energy_id) + mole_flux * uH


      !! by `energy flux' mean the term moleFlux * uH
      call EnergyFluxDerivs(d_energy_flux_up, d_mole_flux_up, uH, up_scale, mole_flux, &
                            toil_auxvar_up%d%dH_dp(iphase), toil_auxvar_up%d%dH_dt(iphase))
      call EnergyFluxDerivs(d_energy_flux_dn, d_mole_flux_dn, uH, dn_scale, mole_flux, &
                            toil_auxvar_dn%d%dH_dp(iphase), toil_auxvar_dn%d%dH_dt(iphase))
      jup(energy_id, 1:3) = jup(energy_id, 1:3) + d_energy_flux_up
      jdn(energy_id, 1:3) = jdn(energy_id, 1:3) + d_energy_flux_dn


!#ifdef DEBUG_FLUXES  
!      adv_flux(energy_id) = adv_flux(energy_id) + mole_flux * uH
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      debug_dphi(iphase) = delta_pressure
!      debug_flux(energy_id,iphase) = debug_flux(energy_id,iphase) + mole_flux * uH
!#endif

    endif  ! if mobility larger than given tolerance                 

  enddo
#endif 
! TOIL_CONVECTION

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  if (debug_flag > 0) then  
!    write(debug_unit,'(a,7es24.15)') 'delta pressure :', debug_dphi(:)
!    write(debug_unit,'(a,7es24.15)') 'adv flux (liquid):', debug_flux(:,1)
!    write(debug_unit,'(a,7es24.15)') 'adv flux (gas):', debug_flux(:,2)
!  endif
!  debug_flux = 0.d0
!#endif                    

#ifdef TOIL_CONDUCTION
  ! model for liquid + gas
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
  !k_eff_up = thermal_conductivity_up(1) + &
  !           sqrt(gen_auxvar_up%sat(option%liquid_phase)) * &
  !           (thermal_conductivity_up(2) - thermal_conductivity_up(1))
  !k_eff_dn = thermal_conductivity_dn(1) + &
  !           sqrt(gen_auxvar_dn%sat(option%liquid_phase)) * &
  !           (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
  !if (k_eff_up > 0.d0 .or. k_eff_up > 0.d0) then
  !  k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  !else
  !  k_eff_ave = 0.d0
  !endif
  ! considered the formation fully saturated in water for heat conduction 
  k_eff_up = thermal_conductivity_up(1)
  k_eff_dn = thermal_conductivity_dn(1)
  if (k_eff_up > 0.d0 .or. k_eff_dn > 0.d0) then
    k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  else
    k_eff_ave = 0.d0
  endif

  ! units:
  ! k_eff = W/K-m = J/s/K-m
  ! delta_temp = K
  ! area = m^2
  ! heat_flux = k_eff * delta_temp * area = J/s
  delta_temp = toil_auxvar_up%temp - toil_auxvar_dn%temp

  heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! J/s -> MJ/s

  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux

  !!  analytical derivatives:
  d_delta_temp_dt_up = 1.d0
  d_delta_temp_dt_dn = - 1.d0

  dheat_flux_ddelta_temp = k_eff_ave * area * 1.d-6 ! J/s -> MJ/s

  jup(energy_id, 3) = jup(energy_id, 3) + d_delta_temp_dt_up*dheat_flux_ddelta_temp
  jdn(energy_id, 3) = jdn(energy_id, 3) + d_delta_temp_dt_dn*dheat_flux_ddelta_temp


! CONDUCTION
#endif

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  debug_flux(energy_id,1) = debug_flux(energy_id,1) + heat_flux
!  if (debug_flag > 0) then  
!    write(debug_unit,'(a,7es24.15)') 'dif flux (liquid):', debug_flux(:,1)
!    write(debug_unit,'(a,7es24.15)') 'dif flux (gas):', debug_flux(:,2)
!  endif
!#endif

end subroutine TOilImsFluxPFL_derivs

! ************************************************************************** !

function TOilImsAverageDensity_derivs(sat_up,sat_dn,density_up,density_dn, dden_up, dden_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/28/15
  ! 

  implicit none

  PetscReal :: sat_up, sat_dn
  PetscReal :: density_up, density_dn, dden_up, dden_dn

  PetscReal :: TOilImsAverageDensity_derivs

  dden_up = 0.d0
  dden_dn = 0.d0

  if (sat_up < eps ) then
    TOilImsAverageDensity_derivs = density_dn
    dden_dn = 1.d0
  else if (sat_dn < eps ) then 
    TOilImsAverageDensity_derivs = density_up
    dden_up = 1.d0
  else ! in here we could use an armonic average, 
       ! other idea sat weighted average but it needs truncation
    TOilImsAverageDensity_derivs = 0.5d0*(density_up+density_dn)
    dden_up = 0.5d0
    dden_dn = 0.5d0
    !print *, "split"
  end if

end function TOilImsAverageDensity_derivs
! ************************************************************************** !

subroutine DeltaPressureDerivs_up_and_down(ddelta_pressure_dpup, ddelta_pressure_dpdn, &
                                           ddelta_pressure_dTup, ddelta_pressure_dTdn, &
                                           dist_gravity,                               &
                                           ddensity_ave_dden_up, ddensity_ave_dden_dn, &
                                           dden_dp_up, dden_dp_dn,                     &
                                           dden_dt_up, dden_dt_dn,                     &
                                           dp_ds_up, dp_ds_dn,                         &
                                           ddelta_pressure_dsatup, ddelta_pressure_dsatdn, &
                                           fmw)


  implicit none

  PetscReal :: ddelta_pressure_dpup, ddelta_pressure_dpdn
  PetscReal :: ddelta_pressure_dTup, ddelta_pressure_dTdn
  PetscReal :: ddelta_pressure_dsatup, ddelta_pressure_dsatdn
  PetscReal :: dist_gravity
  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: dden_dp_up, dden_dp_dn
  PetscReal :: dden_dt_up, dden_dt_dn
  PetscReal :: dp_ds_up, dp_ds_dn
  PetscReal :: fmw
  PetscReal :: fmw_use


  !fmw_use = 0.d0
  fmw_use = fmw


  !gravity_term = density_kg_ave * dist_gravity
  !delta_pressure = toil_auxvar_up%pres(iphase) - &
                   !toil_auxvar_dn%pres(iphase) + &
                   !gravity_term

  ! w.r.t. pressure:
  ddelta_pressure_dpup = 1.d0 + dist_gravity * &
                         ddensity_ave_dden_up * &
                         dden_dp_up * fmw_use
  ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                         ddensity_ave_dden_dn * &
                         dden_dp_up * fmw_use

  ! w.r.t. saturation:
  ddelta_pressure_dsatup = dp_ds_up

  ddelta_pressure_dsatdn = -1.0*dp_ds_dn

#if 0
  ddelta_pressure_dsatup = 0.d0

  ddelta_pressure_dsatdn = 0.d0
#endif


  ! w.r.t. temperature:
  ddelta_pressure_dTup = dist_gravity * &
                         ddensity_ave_dden_up * &
                         dden_dt_up * fmw_use
  ddelta_pressure_dTdn = dist_gravity * &
                         ddensity_ave_dden_dn * &
                         dden_dt_dn * fmw_use


#if 0
    !! adapted from general_common:
    ddelta_pressure_dpup = 1.d0 + dist_gravity * &
                           ddensity_kg_ave_dden_kg_up * &
                           gen_auxvar_up%d%denl_pl * fmw_comp(iphase)
    ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                           ddensity_kg_ave_dden_kg_dn * &
                           gen_auxvar_dn%d%denl_pl * fmw_comp(iphase)
    ddelta_pressure_dTup = dist_gravity * ddensity_kg_ave_dden_kg_up * &
                           gen_auxvar_up%d%denl_T * fmw_comp(iphase)
    ddelta_pressure_dTdn = dist_gravity * ddensity_kg_ave_dden_kg_dn * &
                           gen_auxvar_dn%d%denl_T * fmw_comp(iphase)
#endif


end subroutine DeltaPressureDerivs_up_and_down

! ************************************************************************** !

subroutine EnergyFluxDerivs(d_energy_flux, d_mole_flux, uH, updn_scale, mole_flux, &
                            dH_dp, dH_dt)

  implicit none

  PetscReal, dimension(1:3) :: d_energy_flux
  PetscReal, dimension(1:3) :: d_mole_flux
  PetscReal :: uH, updn_scale, mole_flux
  PetscReal :: dH_dp, dH_dt

  !Res(energy_id) = Res(energy_id) + mole_flux * uH

  ! w.r.t. oil pressure
  d_energy_flux(1) = d_mole_flux(1)*uH + updn_scale*mole_flux*dH_dp

  ! w.r.t. oil saturation
  d_energy_flux(2) = d_mole_flux(2)*uH

  ! w.r.t. temperature
  d_energy_flux(3) = d_mole_flux(3)*uH + updn_scale*mole_flux*dH_dt

end subroutine EnergyFluxDerivs

! ************************************************************************** !

subroutine v_darcy_derivs(d_v_darcy, dmobility, updn_scale, delta_pressure, mobility, &
                          perm_ave_over_dist, d_delta_pres_dp, d_delta_pres_dt, &
                          d_delta_pres_ds)

  implicit none

  PetscReal, dimension(1:3) ::  d_v_darcy
  PetscReal, dimension(1:3) ::  dmobility
  PetscReal :: updn_scale, delta_pressure, mobility, perm_ave_over_dist
  PetscReal :: d_delta_pres_dp, d_delta_pres_dt
  PetscReal :: d_delta_pres_ds


  !v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure

   d_v_darcy = 0.d0

  ! w.r.t. oil pressure
  d_v_darcy(1) = updn_scale*dmobility(1)*delta_pressure + mobility*d_delta_pres_dp

  ! w.r.t. oil saturation
  d_v_darcy(2) = updn_scale*dmobility(2)*delta_pressure + mobility*d_delta_pres_ds

  ! w.r.t. temperature 
  d_v_darcy(3) = updn_scale*dmobility(3)*delta_pressure + mobility*d_delta_pres_dt


  !! scale by perm
  d_v_darcy = perm_ave_over_dist*d_v_darcy

end subroutine V_Darcy_Derivs

! ************************************************************************** !

subroutine MoleFluxDerivs(d_mole_flux, d_q, q, density_ave, ddensity_ave_dden, &
                          dden_dp, dden_dt)

  implicit none
  
  PetscReal, dimension(1:3) :: d_mole_flux
  PetscReal, dimension(1:3) :: d_q
  PetscReal :: q, density_ave
  PetscReal :: ddensity_ave_dden, dden_dp, dden_dt


  !mole_flux = q*density_ave
  ! d_mole_flux(i) = deriv of mole flux w.r.t. variable i
  ! 
  ! variables are:
  !
  ! pres
  ! sat 
  ! temp

  ! assume access to 
  !
  ! q, density ave
  ! dq: derivative array
  ! 
  ! ddensity_ave_dden
  ! so can chain rule with other density derivatives

  d_mole_flux = 0.d0

  ! w.r.t. oil pressure
  d_mole_flux(1) =  d_q(1)*density_ave + q*ddensity_ave_dden*dden_dp
  !                                       (  d (ave den) / d (p)    )

  ! w.r.t. oil sat
  d_mole_flux(2) = d_q(2)*density_ave

  ! w.r.t. temperature 
  d_mole_flux(3) =  d_q(3)*density_ave + q*ddensity_ave_dden*dden_dt
  !                                       (  d (ave den) / d (t)    )

end subroutine MoleFluxDerivs

! ************************************************************************** !

subroutine toil_accum_derivs_alyt(toil_auxvar,material_auxvar, option, j, soil_heat_capacity)

  use Material_Aux_class
  use Option_module

  implicit none

  !! Inputs:
  type(auxvar_toil_ims_type) :: toil_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  !! Outputs:
  PetscReal, dimension(1:3,1:3) :: j !! entry j is partial accum(j) / partial v_j
  !! workers:
  PetscReal :: porosity, volume
  PetscInt :: oid, lid, energy_id, iphase

  oid = option%oil_phase
  lid = option%liquid_phase
  energy_id = option%energy_id
  volume = material_auxvar%volume
  porosity = toil_auxvar%effective_porosity


  j = 0.d0


  !!    OIL EQUATION:

  !! w.r.t. pressure
  !print *, "ddendp: ", toil_auxvar%d%dden_dp(oid,1)
  j(oid, 1) = porosity*toil_auxvar%sat(oid)*toil_auxvar%d%dden_dp(oid,1) + &
              toil_auxvar%d%dpor_dp*(toil_auxvar%sat(oid)*toil_auxvar%den(oid))

  !! w.r.t. sat:
  j(oid, 2)  = porosity*toil_auxvar%den(oid)

  !! w.r.t. temp:
  j(oid,3) = porosity*toil_auxvar%sat(oid)*toil_auxvar%d%dden_dT(oid)


  !! END OIL EQUATION


  !!     LIQUID EQUATION

  !! w.r.t. pressure:
  !print *, "ddendp: ", toil_auxvar%d%dden_dp(lid,1)
  j(lid,1)  = porosity*toil_auxvar%sat(lid)*toil_auxvar%d%dden_dp(lid,1) + &
              toil_auxvar%d%dpor_dp*(toil_auxvar%sat(lid)*toil_auxvar%den(lid))

  !! w.r.t. sat:
  j(lid,2)  = -1.d0*toil_auxvar%den(lid)*porosity

  !! w.r.t. temp
  j(lid,3) =  porosity*toil_auxvar%sat(lid)*toil_auxvar%d%dden_dT(lid)

  !! END LIQUID EQUATION


!!! ADD D_por PARTS TO ALL THESE
  !!     ENERGY EQUATION

  !! first a sum over the two phases, with the term being
  !! 
  !!    (poro) (den) (U)

  !!  liquid phase
  !! 
  !! w.r.t pressure
  j(energy_id, 1) = j(energy_id, 1) + & 
                    toil_auxvar%sat(lid)* ( &
                    toil_auxvar%d%dden_dp(lid,1)*toil_auxvar%U(lid) + &
                    toil_auxvar%den(lid)*toil_auxvar%d%dU_dp(lid) )

  !! w.r.t oil sat:
  j(energy_id, 2) = j(energy_id, 2) - & !! note negative, next term is scaled by dsl/dso 
                    toil_auxvar%den(lid)*toil_auxvar%U(lid)

  !! w.r.t. temp
  j(energy_id,3) = j(energy_id,3) + &
                   toil_auxvar%sat(lid)* ( &
                   toil_auxvar%d%dden_dt(lid)*toil_auxvar%U(lid) + &
                   toil_auxvar%den(lid)*toil_auxvar%d%dU_dT(lid)  )

  !!  oil phase
  !!
  !! w.r.t pressure
  j(energy_id, 1) = j(energy_id, 1) + & 
                    toil_auxvar%sat(oid)* ( &
                    toil_auxvar%d%dden_dp(oid,1)*toil_auxvar%U(oid) + &
                    toil_auxvar%den(oid)*toil_auxvar%d%dU_dp(oid) )

  !! w.r.t oil sat:
  j(energy_id, 2) = j(energy_id, 2) + & 
                    toil_auxvar%den(oid)*toil_auxvar%U(oid)

  !! w.r.t. temp
  j(energy_id,3) = j(energy_id,3) + &
                   toil_auxvar%sat(oid)* ( &
                   toil_auxvar%d%dden_dt(oid)*toil_auxvar%U(oid) + &
                   toil_auxvar%den(oid)*toil_auxvar%d%dU_dT(oid)  )


  !! all should be scaled by poro:
  j(energy_id,1)  = porosity*j(energy_id,1)
  j(energy_id,2)  = porosity*j(energy_id,2)
  j(energy_id,3)  = porosity*j(energy_id,3)
                    
  !! also the (1-por) ... term
  j(energy_id,1) = j(energy_id,1) - toil_auxvar%d%dpor_dp*material_auxvar%soil_particle_density*soil_heat_capacity*toil_auxvar%temp
  j(energy_id,3) = j(energy_id,3) + (1.d0 - porosity)*material_auxvar%soil_particle_density * &
                                                             soil_heat_capacity


# if 0
  do iphase = 1, option%nphase
    j(energy_id,1) = j(energy_id,1) + &
                    toil_auxvar%sat(iphase)* ( &
                    toil_auxvar%d%dden_dp(iphase,1)*toil_auxvar%U(iphase) + &
                    toil_auxvar%den(iphase)*toil_auxvar%d%dU_dp(iphase) )

    j(energy_id,3) = j(energy_id,3) + &
                     toil_auxvar%sat(iphase)* ( &
                     toil_auxvar%d%dden_dt(iphase)*toil_auxvar%U(iphase) + &
                     toil_auxvar%den(iphase)*toil_auxvar%d%dU_dT(iphase)  )
  end do

  j(energy_id,1)  = porosity*j(energy_id,1)
  j(energy_id,2)  = porosity*j(energy_id,2)

  j(energy_id,3) = j(energy_id,3) + (1.d0 - porosity)*material_auxvar%soil_particle_density * &
                                                             soil_heat_capacity

  j(energy_id,3) = j(energy_id,3)

  j(energy_id,1) = j(energy_id,1)*porosity

# endif

  !! END ENERGY EQUATION


j = j*volume

end subroutine toil_accum_derivs_alyt

! ************************************************************************** !

#if 0
  do iphase = 1, option%nphase

    denomp = denomp + &
             toil_auxvar%sat(iphase)*(toil_auxvar%derivatives%dden_dp(iphase,1)*toil_auxvar%U(iphase) + &
                                      toil_auxvar%den(iphase)*toil_auxvar%derivatives%dU_dp(iphase))

    !! denoms = 0? Density and U should be independent of s - check again
    !! ComputeAuxVars for hidden dependencies. Should see agreement from numerical
    !! derivs.
    !denoms = denoms + &
             !!toil_auxvar%sat(iphase)*(toil_auxvar%derivatives%dden_ds(iphase)*toil_auxvar%U(iphase) + &
              !toil_auxvar%sat(iphase)*toil_auxvar%den(iphase)*toil_auxvar%derivatives%dU_ds(iphase)
              

    denomT = denomT + &
             toil_auxvar%sat(iphase)*(toil_auxvar%derivatives%dden_dt(iphase)*toil_auxvar%U(iphase) + &
                                      toil_auxvar%den(iphase)*toil_auxvar%derivatives%dU_dT(iphase))

  end do
  denomT = porosity*denomT
  denomT = denomT + (1.d0 - porosity)*material_auxvar%soil_particle_density * &
                                                             soil_heat_capacity
  denomT = denomT/volume

  porVol = porosity/volume
  denomp = denomp*porVol
  denoms = denoms*porVol

#endif

#if 0
   denomp = porosity*toil_auxvar%sat(lid)*toil_auxvar%derivatives%dden_dp(lid,1)/volume !! porosity * oil saturation * 
                                                                            !! d (oil mole frac)/ d(oil pres)
  ! denomp2 = denom*auxvar%dden_dp

   denoms = -1.d0*toil_auxvar%den(lid)  !! could be dden_ds part here to but that (should be) always zero
   denoms = porosity*denoms/volume

  ! denoms_sq = denoms**2

   denomT = porosity*toil_auxvar%sat(lid)*toil_auxvar%derivatives%dden_dT(lid)/volume
#endif

#if 0

subroutine u1(toil_auxvar,material_auxvar, option, for_accum, j_1) !! inputs: auxvar, outputs: two arrays j_1 j_2; include option to not compute j_2
  use Material_Aux_class
  use Option_module

  implicit none
  !! Inputs:
  type(auxvar_toil_ims_type) :: toil_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  !! Outputs:
  PetscReal, dimension(1:3) :: j_1 !! entry j is partial v_j / partial u
  ! PETSCREAL, dimension(1:3, 1:3) :: j_2 !! entry i,j is partial^2 v_i /  partial u partial v_j
  !! Internal:
  PetscReal :: denomp, denoms, denomT !! identify parts of computation that are reused
  PetscReal :: porosity, volume
  PetscInt :: oid 
  PetscBool :: for_accum !! if this is true return just the reciprocals of j1 - these are what are needed for analytical derivatives of the accumulation part.

  !! Pretty much just copy equations from paper/Maple and do them here one after
  !! the other

  oid = option%oil_phase
  volume = material_auxvar%volume

  porosity = toil_auxvar%effective_porosity
  
   denomp = porosity*toil_auxvar%sat(oid)*toil_auxvar%derivatives%dden_dp(oid,1)/volume !! porosity * oil saturation * 
                                                                            !! d (oil mole frac)/ d(oil pres)
  ! denomp2 = denom*auxvar%dden_dp

   denoms = porosity*toil_auxvar%den(oid)/volume
  ! denoms_sq = denoms**2

   denomT = porosity*toil_auxvar%sat(oid)*toil_auxvar%derivatives%dden_dT(oid)/volume
  ! denomT_sq = denomT**2

  if (for_accum) then
    j_1(1) = denomp
    j_1(2) = denoms
    j_1(3) = denomT
  else
   
    j_1(1) = 1.d0/denomp !! what to do when some denom are 0? Likely to happen.
    j_1(2) = 1.d0/denoms
    j_1(3) = 1.d0/denomT

    !j_2(1,1) = -1.d0*auxvar%ddden_dpp/denomp2
    !j_2(2,1) = -1.d0*phi*auxvar%dden_dp/denoms_sq
    !j_2(3,1) = -1.d0*phi*sat*auxvar%ddden_dT_dp/denomT_sq !could be cleaned up(cancellations)

    !j_2(1,2) = -1.d0/denomp2/auxvar%sat[oid]
    !j_2(2,2) = -1.d0*phi*auxvar%dden_dsat_oil/denoms_sq
    !j_2(3,2) = -1.d0*phi*auxvar%dden_dT/denomsT_sq

    !j_2(1,3) = -1.d0*auxvar%ddden/dpdT/denomp2
    !j_2(2,3) = -1.d0*phi*auxvar%dden_dT/denoms_sq
    !j_2(3,3) = -1.d0*phi*auxvar%ddden_dTT/denomsT_sq
  end if

end subroutine u1

#endif







end module TOilIms_derivs_module
