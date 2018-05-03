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

  public :: toil_accum_derivs_alyt, &
            TOilImsAverageDensity_derivs, &
            MoleFluxDerivs, &
            EnergyDrivenFluxDerivs, &
            DeltaPressureDerivs_up_and_down, &
            V_Darcy_Derivs, &
            InjectionEnergyPartDerivs, &
            Qsrc_mol_derivs

contains


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

function TOilImsAverageDensity_derivs(sat_up,sat_dn,density_up,density_dn, dden_up, dden_dn)
  ! 
  ! Modification of TOilImsAverageDensity which computes derivatives
  ! Daniel Stone, March-May 2018
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

  ! w.r.t. temperature:
  ddelta_pressure_dTup = dist_gravity * &
                         ddensity_ave_dden_up * &
                         dden_dt_up * fmw_use
  ddelta_pressure_dTdn = dist_gravity * &
                         ddensity_ave_dden_dn * &
                         dden_dt_dn * fmw_use


end subroutine DeltaPressureDerivs_up_and_down

! ************************************************************************** !

subroutine EnergyDrivenFluxDerivs(d_energy_flux, d_mole_flux, uH, updn_scale, mole_flux, &
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

end subroutine EnergyDrivenFluxDerivs

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

  !! first a sum over the two phases, with the term being
  !! 
  !!    (poro) (sat) (den) (U)

  !!  liquid phase
  !! 
  !! w.r.t pressure
  j(energy_id, 1) = j(energy_id, 1) + & 
                    porosity * &
                    toil_auxvar%sat(lid)* ( &
                    toil_auxvar%d%dden_dp(lid,1)*toil_auxvar%U(lid) + &
                    toil_auxvar%den(lid)*toil_auxvar%d%dU_dp(lid) )
  !! and density w.r.t. pressure:
  j(energy_id, 1) = j(energy_id, 1) + & 
                    toil_auxvar%d%dpor_dp*toil_auxvar%sat(lid)*toil_auxvar%den(lid)*toil_auxvar%u(lid)

  !! w.r.t oil sat:
  j(energy_id, 2) = j(energy_id, 2) - & !! note negative, next term is scaled by dsl/dso 
                    porosity*toil_auxvar%den(lid)*toil_auxvar%U(lid)

  !! w.r.t. temp
  j(energy_id,3) = j(energy_id,3) + &
                   porosity * &
                   toil_auxvar%sat(lid)* ( &
                   toil_auxvar%d%dden_dt(lid)*toil_auxvar%U(lid) + &
                   toil_auxvar%den(lid)*toil_auxvar%d%dU_dT(lid)  )

  !!  oil phase
  !!
  !! w.r.t pressure
  j(energy_id, 1) = j(energy_id, 1) + & 
                    porosity * &
                    toil_auxvar%sat(oid)* ( &
                    toil_auxvar%d%dden_dp(oid,1)*toil_auxvar%U(oid) + &
                    toil_auxvar%den(oid)*toil_auxvar%d%dU_dp(oid) )
  !! and density w.r.t. pressure:
  j(energy_id, 1) = j(energy_id, 1) + & 
                    toil_auxvar%d%dpor_dp*toil_auxvar%sat(oid)*toil_auxvar%den(oid)*toil_auxvar%u(oid)

  !! w.r.t oil sat:
  j(energy_id, 2) = j(energy_id, 2) + & 
                    porosity*toil_auxvar%den(oid)*toil_auxvar%U(oid)

  !! w.r.t. temp
  j(energy_id,3) = j(energy_id,3) + &
                   porosity * &
                   toil_auxvar%sat(oid)* ( &
                   toil_auxvar%d%dden_dt(oid)*toil_auxvar%U(oid) + &
                   toil_auxvar%den(oid)*toil_auxvar%d%dU_dT(oid)  )

  !! also the (1-por) ... term
  j(energy_id,1) = j(energy_id,1) - toil_auxvar%d%dpor_dp*material_auxvar%soil_particle_density*soil_heat_capacity*toil_auxvar%temp
  j(energy_id,3) = j(energy_id,3) + (1.d0 - porosity)*material_auxvar%soil_particle_density * &
                                                             soil_heat_capacity

  !! END ENERGY EQUATION


j = j*volume

end subroutine toil_accum_derivs_alyt

! ************************************************************************** !
end module TOilIms_derivs_module
