module Transport_NP_module
! This module handles specific calls from reactive_transport.f90 to
! take into account electromigration and specific
! diffusion terms (Nernst-Plank transport formulation)
! Author: Albert Nardi (Amphos21 - Barcelona Science)

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use Matrix_Block_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: TNPFlux, &
            TNPFluxBC, &
            TNPFluxDerivative, &
            TNPFluxDerivativeBC, &
            ComputeElectricPotentialTotalComponent

contains

! ************************************************************************** !
subroutine TNPFlux(reaction, &
                 rt_parameter, &
                 rt_auxvar_up, material_auxvar_up, global_auxvar_up, &
                 rt_auxvar_dn, material_auxvar_dn, global_auxvar_dn, &
                 dist, &
                 area, &
                 option, &
                 Res)
  !
  ! Computes flux term in residual function
  !
  ! Author: Albert Nardi (Amphos21 - Barcelona Science)
  ! Date: 29/06/20
  !
  use Option_module

  implicit none

  class(reaction_rt_type), pointer :: reaction
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_auxvar_up, rt_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(option_type) :: option

  PetscReal :: dist(-1:3)
  PetscReal :: area

  PetscReal :: Res(rt_parameter%ncomp)
  PetscReal :: res_cont(rt_parameter%ncomp) ! local contribution to residual
  PetscReal :: factor_up, factor_dn
  PetscReal :: factor_avg_em(rt_parameter%naqcomp)
  PetscReal :: tort_up, tort_dn
  PetscReal :: sat_up, sat_dn
  PetscReal :: harm_factor, cur_factor
  PetscReal :: cur_diff, curStoich
  PetscReal :: dist_up, dist_dn
  PetscReal :: potent_avg(rt_parameter%naqcomp)
  PetscReal :: sum_denom_avg, sum_transp
  PetscReal :: curcharge
  PetscReal :: curConc_up, curConc_dn

  PetscInt :: iphase, i
  PetscInt :: icplx, icomp, jcomp, lcomp, kcomp
  PetscInt :: ndof, ncomp

  iphase = option%liquid_phase
  ndof = rt_parameter%naqcomp

  dist_up = dist(0)*dist(-1)
  dist_dn = dist(0)-dist_up ! should avoid truncation error

  sat_up = global_auxvar_up%sat(iphase)
  sat_dn = global_auxvar_dn%sat(iphase)

  tort_up =  material_auxvar_up%tortuosity
  tort_dn =  material_auxvar_dn%tortuosity

  ! 1000 converts m^3 -> L
  factor_up = max(sat_up * material_auxvar_up%porosity * &
          tort_up * global_auxvar_up%den_kg(iphase)*1.d-3 * &
          area*1000.d0, 1d-40)
  factor_dn = max(sat_dn * material_auxvar_dn%porosity * &
          tort_dn * global_auxvar_dn%den_kg(iphase)*1.d-3* &
          area*1000.d0, 1d-40)

  harm_factor = factor_up*factor_dn/(factor_up*dist_dn + factor_dn*dist_up)

  res_cont = 0d0

  ! D_j \nabla C_j
  do icomp = 1, rt_parameter%naqcomp
    cur_factor = harm_factor*rt_parameter%pri_spec_diff_coef(icomp)
    res_cont(icomp) = cur_factor*(rt_auxvar_up%pri_molal(icomp) - rt_auxvar_dn%pri_molal(icomp))
  enddo

  ! summatory D_i \nabla C_i
  do icplx = 1, reaction%neqcplx
    cur_factor = harm_factor*rt_parameter%sec_spec_diff_coef(icplx)
    ncomp = reaction%eqcplxspecid(0,icplx)
    do i = 1, ncomp
      icomp = reaction%eqcplxspecid(i,icplx)
      res_cont(icomp) = res_cont(icomp) + &
                        cur_factor*reaction%eqcplxstoich(i,icplx)* &
                        (rt_auxvar_up%sec_molal(icplx) - rt_auxvar_dn%sec_molal(icplx))
    enddo
  enddo

  ! Electromigration term
  call ComputeElectricPotentialTotalComponent(reaction, &
                                            rt_parameter,  &
                                            rt_auxvar_up, &
                                            rt_auxvar_dn, &
                                            icomp, &
                                            potent_avg)
  sum_denom_avg = 1d-40
  do jcomp = 1, rt_parameter%naqcomp
    curCharge = reaction%primary_spec_Z(jcomp)
    sum_denom_avg = sum_denom_avg + potent_avg(jcomp)*curCharge
  enddo

  do jcomp = 1, rt_parameter%naqcomp
    factor_avg_em(jcomp) = potent_avg(jcomp) / sum_denom_avg
  enddo

  sum_transp = 0d0
  ! Primary l index
  do lcomp = 1, rt_parameter%naqcomp
    curCharge = reaction%primary_spec_Z(lcomp)
    cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)
    curConc_up = rt_auxvar_up%pri_molal(lcomp)
    curConc_dn = rt_auxvar_dn%pri_molal(lcomp)
    sum_transp = sum_transp + curCharge*cur_diff* &
                    (curConc_up - &
                     curConc_dn)
  enddo

  ! Secondary term k
  do icplx = 1, reaction%neqcplx
    ncomp = reaction%eqcplxspecid(0,icplx)
    do i = 1, ncomp
      kcomp = reaction%eqcplxspecid(i,icplx)
      curCharge = reaction%primary_spec_Z(kcomp)
      cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
      curStoich = reaction%eqcplxstoich(i,icplx)
      curConc_up = rt_auxvar_up%sec_molal(icplx)
      curConc_dn = rt_auxvar_dn%sec_molal(icplx)
      sum_transp = sum_transp + curCharge*curStoich*cur_diff* &
                    (curConc_up - &
                     curConc_dn)
    enddo
  enddo

  do icomp = 1, rt_parameter%naqcomp
    res_cont(icomp) = res_cont(icomp) - harm_factor*factor_avg_em(icomp)*sum_transp
  enddo

  ! Add to residual
  Res(1:ndof) = Res(1:ndof) + res_cont(1:ndof)

end subroutine TNPFlux

! ************************************************************************** !
subroutine TNPFluxBC( &
                 ibndtype, &
                 reaction, &
                 rt_parameter, &
                 rt_auxvar_up, &
                 global_auxvar_up, &
                 rt_auxvar_dn, &
                 material_auxvar_dn, &
                 global_auxvar_dn, &
                 dist_dn, &
                 area, &
                 option, &
                 Res)
  !
  ! Computes flux term in BC residual function
  !
  ! Author: Albert Nardi (Amphos21 - Barcelona Science)
  ! Date: 29/06/20
  !
  use Option_module

  implicit none

  PetscInt :: ibndtype
  class(reaction_rt_type), pointer :: reaction
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_auxvar_up, rt_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_dn, global_auxvar_up
  type(option_type) :: option


  PetscReal :: dist_dn(-1:3)
  PetscReal :: area

  PetscReal :: Res(rt_parameter%ncomp)
  PetscReal :: res_cont(rt_parameter%ncomp) ! local contribution to residual
  PetscReal :: factor_dn
  PetscReal :: factor_avg_em(rt_parameter%naqcomp)
  PetscReal :: tort_dn
  PetscReal :: sat_dn
  PetscReal :: harm_factor, cur_factor
  PetscReal :: cur_diff, curStoich
  PetscReal :: potent_avg(rt_parameter%naqcomp)
  PetscReal :: sum_denom_avg, sum_transp
  PetscReal :: curcharge
  PetscReal :: curConc_up, curConc_dn

  PetscInt :: iphase, i
  PetscInt :: icplx, icomp, jcomp, lcomp, kcomp
  PetscInt :: ndof, ncomp

  iphase = option%liquid_phase
  ndof = rt_parameter%naqcomp

  select case(ibndtype)
      case(DIRICHLET_BC, DIRICHLET_ZERO_GRADIENT_BC)
          ndof = rt_parameter%naqcomp
          sat_dn = global_auxvar_dn%sat(iphase)
          tort_dn =  material_auxvar_dn%tortuosity
          factor_dn = max(sat_dn * material_auxvar_dn%porosity * &
                  tort_dn * global_auxvar_dn%den_kg(iphase)*1.d-3* &
                  area*1000.d0, 1d-40) ! 1000 converts m^3 -> L
          harm_factor = factor_dn / dist_dn(0)
          res_cont = 0d0

          ! D_j \nabla C_j
          do icomp = 1, rt_parameter%naqcomp
            cur_factor = harm_factor*rt_parameter%pri_spec_diff_coef(icomp)
            res_cont(icomp) = cur_factor*(rt_auxvar_up%pri_molal(icomp) - rt_auxvar_dn%pri_molal(icomp))
          enddo

          ! sumatory D_i \nabla C_i
          do icplx = 1, reaction%neqcplx
            cur_factor = harm_factor*rt_parameter%sec_spec_diff_coef(icplx)
            ncomp = reaction%eqcplxspecid(0,icplx)
            do i = 1, ncomp
              icomp = reaction%eqcplxspecid(i,icplx)
              res_cont(icomp) = res_cont(icomp) + &
                                cur_factor*reaction%eqcplxstoich(i,icplx)* &
                                (rt_auxvar_up%sec_molal(icplx) - rt_auxvar_dn%sec_molal(icplx))
            enddo
          enddo

          ! Electromigration term
          call ComputeElectricPotentialTotalComponent(reaction, &
                                                    rt_parameter,  &
                                                    rt_auxvar_up, &
                                                    rt_auxvar_dn, &
                                                    icomp, &
                                                    potent_avg)

          sum_denom_avg = 1d-40
          do jcomp = 1, rt_parameter%naqcomp
            curCharge = reaction%primary_spec_Z(jcomp)
            sum_denom_avg = sum_denom_avg + potent_avg(jcomp)*curCharge
          enddo

          do jcomp = 1, rt_parameter%naqcomp
            factor_avg_em(jcomp) = potent_avg(jcomp) / sum_denom_avg
          enddo

          sum_transp = 0d0
          ! Primary l index
          do lcomp = 1, rt_parameter%naqcomp
            curCharge = reaction%primary_spec_Z(lcomp)
            cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)
            curConc_up = rt_auxvar_up%pri_molal(lcomp)
            curConc_dn = rt_auxvar_dn%pri_molal(lcomp)

            sum_transp = sum_transp + curCharge*cur_diff* &
                            (curConc_up - &
                             curConc_dn)
          enddo

          ! Secondary term k
          do icplx = 1, reaction%neqcplx
            ncomp = reaction%eqcplxspecid(0,icplx)
            do i = 1, ncomp
              kcomp = reaction%eqcplxspecid(i,icplx)
              curCharge = reaction%primary_spec_Z(kcomp)
              cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
              curStoich = reaction%eqcplxstoich(i,icplx)
              curConc_up = rt_auxvar_up%sec_molal(icplx)
              curConc_dn = rt_auxvar_dn%sec_molal(icplx)
              sum_transp = sum_transp + curCharge*curStoich*cur_diff* &
                            (curConc_up - &
                             curConc_dn)
            enddo
          enddo

          do icomp = 1, rt_parameter%naqcomp
            res_cont(icomp) = res_cont(icomp) - harm_factor*factor_avg_em(icomp)*sum_transp
          enddo

          ! Add to residual
          Res(1:ndof) = Res(1:ndof) + res_cont(1:ndof)

      case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)

    end select


end subroutine TNPFluxBC

! ************************************************************************** !
subroutine TNPFluxDerivative(reaction, &
                 rt_parameter, &
                 rt_auxvar_up, material_auxvar_up, global_auxvar_up, &
                 rt_auxvar_dn, material_auxvar_dn, global_auxvar_dn, &
                 dist, &
                 area, &
                 option, &
                 J_up, &
                 J_dn)
  !
  ! Computes function derivative term. The electromigration term
  !
  ! Author: Albert Nardi (Amphos21 - Barcelona Science)
  ! Date: 29/06/20
  !
  use Option_module

  implicit none

  class(reaction_rt_type), pointer :: reaction
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_auxvar_up, rt_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(option_type) :: option


  PetscReal :: dist(-1:3)
  PetscReal :: area

  PetscReal :: J_up(rt_parameter%ncomp,rt_parameter%ncomp), &
               J_dn(rt_parameter%ncomp,rt_parameter%ncomp)
  PetscReal :: factor_up, factor_dn
  PetscReal :: tort_up, tort_dn
  PetscReal :: sat_up, sat_dn
  PetscReal :: harm_factor, cur_factor
  PetscReal :: dist_up, dist_dn

  PetscInt :: iphase, i, j
  PetscInt :: icplx, icomp, jcomp
  PetscInt :: ndof, ncomp
  PetscInt :: istart
  PetscInt :: iend

  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: lnQK, tempreal_up, tempreal_dn
  PetscReal :: den_kg_per_L, xmass

  iphase = option%liquid_phase
  ndof = rt_parameter%naqcomp

  xmass = 1.d0
  if (associated(global_auxvar_up%xmass)) xmass = global_auxvar_up%xmass(iphase)
  den_kg_per_L = global_auxvar_up%den_kg(iphase)*xmass*1.d-3

  ln_conc = log(rt_auxvar_up%pri_molal)
  ln_act = ln_conc+log(rt_auxvar_up%pri_act_coef)
  lnQK = 0d0

  dist_up = dist(0)*dist(-1)
  dist_dn = dist(0)-dist_up ! should avoid truncation error

  sat_up = global_auxvar_up%sat(iphase)
  sat_dn = global_auxvar_dn%sat(iphase)

  tort_up =  material_auxvar_up%tortuosity
  tort_dn =  material_auxvar_dn%tortuosity

  factor_up = max(sat_up * material_auxvar_up%porosity * &
          tort_up * global_auxvar_up%den_kg(iphase)*1.d-3 * &
          area*1000.d0, 1d-40) ! 1000 converts m^3 -> L
  factor_dn = max(sat_dn * material_auxvar_dn%porosity * &
          tort_dn * global_auxvar_dn%den_kg(iphase)*1.d-3* &
          area*1000.d0, 1d-40) ! 1000 converts m^3 -> L

  harm_factor = factor_up*factor_dn/(factor_up*dist_dn + factor_dn*dist_up)

  istart = 1
  iend = rt_parameter%naqcomp

  ! d (D_j * C_j) / dC_j =  D_j
  do icomp = 1, rt_parameter%naqcomp
    cur_factor = harm_factor*rt_parameter%pri_spec_diff_coef(icomp)
    J_up(icomp,icomp) = &
      J_up(icomp,icomp) + &
      cur_factor
    J_dn(icomp,icomp) = &
      J_dn(icomp,icomp) - &
      cur_factor
  enddo

  ! d (sumatory ( frac * D_i * C_i )) / dC_j =
  ! sumatory ( frac * D_i * dC_i/dC_j )
  do icplx = 1, reaction%neqcplx

    lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

    ! activity of water
    if (reaction%eqcplxh2oid(icplx) > 0) then
      lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar_up%ln_act_h2o
    endif

    ncomp = reaction%eqcplxspecid(0,icplx)
    do i = 1, ncomp
      icomp = reaction%eqcplxspecid(i,icplx)
      lnQK = lnQK + reaction%eqcplxstoich(i,icplx)*ln_act(icomp)
    enddo

    cur_factor = harm_factor*rt_parameter%sec_spec_diff_coef(icplx)
    ncomp = reaction%eqcplxspecid(0,icplx)

    do j = 1, ncomp
      jcomp = reaction%eqcplxspecid(j,icplx)
      tempreal_up = reaction%eqcplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                 rt_auxvar_up%sec_act_coef(icplx)
      tempreal_dn = reaction%eqcplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                 rt_auxvar_dn%sec_act_coef(icplx)
      do i = 1, ncomp
        icomp = reaction%eqcplxspecid(i,icplx)
        J_up(icomp, jcomp) = J_up(icomp, jcomp) + &
                          cur_factor* &
                          reaction%eqcplxstoich(i,icplx)* &
                          tempreal_up
        J_dn(icomp, jcomp) = J_dn(icomp, jcomp) - &
                          cur_factor* &
                          reaction%eqcplxstoich(i,icplx)* &
                          tempreal_dn
      enddo
    enddo
  enddo


end subroutine TNPFluxDerivative

! ************************************************************************** !
subroutine TNPFluxDerivativeBC(&
                 ibndtype, &
                 reaction, &
                 rt_parameter, &
                 rt_auxvar_up, global_auxvar_up, &
                 rt_auxvar_dn, material_auxvar_dn, global_auxvar_dn, &
                 dist, &
                 area, &
                 option, &
                 J_up, &
                 J_dn)
  !
  ! Computes flux term in residual function
  !
  ! Author: Albert Nardi (Amphos21 - Barcelona Science)
  ! Date: 29/06/20
  !
  use Option_module

  implicit none

  PetscInt :: ibndtype
  class(reaction_rt_type), pointer :: reaction
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_auxvar_up, rt_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(option_type) :: option


  PetscReal :: dist(-1:3)
  PetscReal :: area

  PetscReal :: J_up(rt_parameter%ncomp,rt_parameter%ncomp), &
               J_dn(rt_parameter%ncomp,rt_parameter%ncomp)
  PetscReal :: factor_dn
  PetscReal :: tort_dn
  PetscReal :: sat_dn
  PetscReal :: harm_factor, cur_factor
  PetscReal :: dist_up, dist_dn

  PetscInt :: iphase, i, j
  PetscInt :: icplx, icomp, jcomp
  PetscInt :: ndof, ncomp
  PetscInt :: istart
  PetscInt :: iend

  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: lnQK, tempreal_up, tempreal_dn
  PetscReal :: den_kg_per_L, xmass

  ! print *, "TNPFluxDerivativeBC"

  iphase = option%liquid_phase
  ndof = rt_parameter%naqcomp

  xmass = 1.d0

  select case(ibndtype)
      case(DIRICHLET_BC, DIRICHLET_ZERO_GRADIENT_BC)
          if (associated(global_auxvar_up%xmass)) xmass = global_auxvar_up%xmass(iphase)
          den_kg_per_L = global_auxvar_up%den_kg(iphase)*xmass*1.d-3

          ln_conc = log(rt_auxvar_dn%pri_molal)
          ln_act = ln_conc+log(rt_auxvar_dn%pri_act_coef)

          do icplx = 1, reaction%neqcplx ! for each secondary species
            ! compute secondary species concentration
            lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

            ! activity of water
            if (reaction%eqcplxh2oid(icplx) > 0) then
              lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar_up%ln_act_h2o
            endif

            ncomp = reaction%eqcplxspecid(0,icplx)
            do i = 1, ncomp
              icomp = reaction%eqcplxspecid(i,icplx)
              lnQK = lnQK + reaction%eqcplxstoich(i,icplx)*ln_act(icomp)
            enddo
          enddo

          dist_up = dist(0)*dist(-1)
          dist_dn = dist(0)-dist_up ! should avoid truncation error

          sat_dn = global_auxvar_dn%sat(iphase)
          tort_dn =  material_auxvar_dn%tortuosity

          factor_dn = max(sat_dn * material_auxvar_dn%porosity * &
            tort_dn * global_auxvar_dn%den_kg(iphase)*1.d-3 * &
            area*1000.d0, 1d-40)

          ! 1000 converts m^3 -> L
          harm_factor = factor_dn/dist_dn

          istart = 1
          iend = rt_parameter%naqcomp

          ! d (D_j * C_j) / dC_j =  D_j
          do icomp = 1, rt_parameter%naqcomp
            cur_factor = harm_factor*rt_parameter%pri_spec_diff_coef(icomp)

            J_up(icomp,icomp) = &
              J_up(icomp,icomp) + &
              cur_factor
            J_dn(icomp,icomp) = &
              J_dn(icomp,icomp) - &
              cur_factor
          enddo

          ! d (sumatory ( frac * D_i * C_i )) / dC_j =
          ! sumatory ( frac * D_i * dC_i/dC_j )
          do icplx = 1, reaction%neqcplx

            cur_factor = harm_factor*rt_parameter%sec_spec_diff_coef(icplx)
            ncomp = reaction%eqcplxspecid(0,icplx)

            do j = 1, ncomp
              jcomp = reaction%eqcplxspecid(j,icplx)
              tempreal_up = reaction%eqcplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                         rt_auxvar_up%sec_act_coef(icplx)
              tempreal_dn = reaction%eqcplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                         rt_auxvar_dn%sec_act_coef(icplx)

              do i = 1, ncomp
                icomp = reaction%eqcplxspecid(i,icplx)
                J_up(icomp, jcomp) = J_up(icomp, jcomp) + &
                                  cur_factor* &
                                  reaction%eqcplxstoich(i,icplx)* &
                                  tempreal_up
                J_dn(icomp, jcomp) = J_dn(icomp, jcomp) - &
                                  cur_factor* &
                                  reaction%eqcplxstoich(i,icplx)* &
                                  tempreal_dn
              enddo
            enddo
          enddo


      case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
    end select

end subroutine TNPFluxDerivativeBC

! ************************************************************************** !
subroutine ComputeElectricPotentialTotalComponent(reaction, &
                                        rt_parameter, &
                                        rt_auxvar_up, &
                                        rt_auxvar_dn, &
                                        icomp, &
                                        potential_avg)
  !
  ! Computes the U electric potential vector term
  !
  ! Author: Albert Nardi (Amphos21 - Barcelona Science)
  ! Date: 29/06/20
  !
  class(reaction_rt_type), pointer :: reaction
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_auxvar_up, rt_auxvar_dn
  PetscInt :: icomp, jcomp
  PetscReal :: potential_avg(reaction%naqcomp)
  PetscReal :: curCoef, curCharge, curStoich
  PetscReal :: curConc_up, curConc_dn, curConc_avg
  PetscInt :: icplx, i, ncomp

  ! Primary species contribution
  do icomp = 1, reaction%naqcomp
    curCoef = rt_parameter%pri_spec_diff_coef(icomp)
    curCharge = reaction%primary_spec_Z(icomp)
    curConc_up = rt_auxvar_up%pri_molal(icomp)
    curConc_dn = rt_auxvar_dn%pri_molal(icomp)
    ! use of differential logarithmic avg
    if (log(curConc_up).eq.log(curConc_dn)) then
       curConc_avg = curConc_up
    else
       curConc_avg = (curConc_up-curConc_dn)/(log(curConc_up)-log(curConc_dn))
    endif
    potential_avg(icomp) = curCoef*curCharge*curConc_avg
  enddo

  ! Secondary species contribution
  do icplx = 1, reaction%neqcplx
    curCharge = reaction%eqcplx_Z(icplx)
    ncomp = reaction%eqcplxspecid(0,icplx)
    do i = 1, ncomp
      jcomp = reaction%eqcplxspecid(i,icplx)
      curCoef = rt_parameter%sec_spec_diff_coef(icplx)
      curStoich = reaction%eqcplxstoich(i,icplx)
      curConc_up = rt_auxvar_up%sec_molal(icplx)
      curConc_dn = rt_auxvar_dn%sec_molal(icplx)

      ! use of differential logarithmic avg
      if (log(curConc_up).eq.log(curConc_dn)) then
          curConc_avg = curConc_up
      else
        curConc_avg = (curConc_up-curConc_dn)/(log(curConc_up)-log(curConc_dn))
      endif
      potential_avg(jcomp) = potential_avg(jcomp) + curStoich * curCoef * curCharge * curConc_avg
    enddo
  enddo

end subroutine ComputeElectricPotentialTotalComponent

end module Transport_NP_module

