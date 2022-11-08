! Added by Amphos 21
! For technical details contact albert.nardi@amphos21.com
! ========================================================

module Secondary_Continuum_NP_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Secondary_Continuum_Aux_module

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal

  implicit none

  private

  PetscReal, parameter :: perturbation_tolerance = 1.d-5

  public :: SecondaryRTResJacMulti_NP, &
            SecondaryRTUpdateIterate_NP, &
            ComputeElectricPotentialTotalComponent_NP

contains

! ************************************************************************** !

subroutine SecondaryRTResJacMulti_NP(sec_transport_vars,auxvar, &
                                  global_auxvar,prim_vol, &
                                  reaction,rt_parameter,diffusion_coefficient, &
                                  porosity,tortuosity,option,res_transport)
  !
  ! RTSecondaryTransportMulti:  Calculates the source term contribution due to
  ! secondary continuum in the primary continuum residual for multicomponent
  ! system assuming only aqueous reaction
  !
  ! Author: Albert Nardi, Amphos 21
  ! Date: 8/31/2021
  !

  use Option_module
  use Global_Aux_module
  use Block_Solve_module
  use Block_Tridiag_module
  use Utility_module
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_module

  implicit none

  type(reactive_transport_param_type) :: rt_parameter
  type(sec_transport_type) :: sec_transport_vars
  type(reactive_transport_auxvar_type) :: auxvar
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type), pointer :: reaction
  type(option_type) :: option
  PetscReal :: coeff_left(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_diag(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right(reaction%naqcomp,reaction%naqcomp, &
                           sec_transport_vars%ncells)
  PetscReal :: res(sec_transport_vars%ncells*reaction%naqcomp)

  PetscReal :: curConc_up, curConc_dn
  PetscReal :: curCharge
  PetscReal :: cur_diff
  PetscReal :: curStoich
  PetscReal :: sum_transp_up, sum_transp_dn
  PetscReal :: sum_denom_up, sum_denom_dn
  PetscReal :: factor_up_em(reaction%naqcomp)
  PetscReal :: factor_dn_em(reaction%naqcomp)
  PetscReal :: potent_up(reaction%naqcomp)
  PetscReal :: potent_dn(reaction%naqcomp)

  PetscReal :: rhs(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: D_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: identity(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: b_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: sec_jac(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: inv_D_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: conc_upd(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: total_upd(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: total_prev(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: conc_current_M(reaction%naqcomp)
  PetscReal :: total_current_M(reaction%naqcomp)
  PetscReal :: res_transport(reaction%naqcomp)
  PetscReal :: total_primary_node(reaction%naqcomp)
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)
  PetscReal :: res_react(reaction%naqcomp)
  PetscReal :: jac_react(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: dtotal(reaction%naqcomp,reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: dtotal_prim(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: mc_pri_molal(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: mc_sec_molal(reaction%neqcplx,sec_transport_vars%ncells)
  PetscInt :: i, j, k, n, l
  PetscInt :: ngcells, ncomp
  PetscInt :: ni
  PetscInt :: nicomp, icomp
  PetscReal :: area_fm
  PetscReal :: diffusion_coefficient
  PetscReal :: porosity
  PetscReal :: tortuosity
  PetscReal :: dC_prim(reaction%naqcomp,reaction%naqcomp)

  PetscReal :: pordt, pordiff
  PetscReal :: pordiff_prim, pordiff_sec
  PetscReal :: portort
  PetscReal :: prim_vol ! volume of primary grid cell
  PetscReal :: dCsec_dCprim(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: dPsisec_dCprim(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: jcomp, lcomp, kcomp, icplx, ncompeq
  PetscReal :: sec_sec_molal_M(reaction%neqcplx)   ! secondary species molality of secondary continuum

  PetscInt :: pivot(reaction%naqcomp,sec_transport_vars%ncells)
  PetscInt :: indx(reaction%naqcomp)
  PetscInt :: d, ier
  PetscReal :: m

  ! Quantities for numerical jacobian
  PetscReal :: conc_prim(reaction%naqcomp)
  PetscReal :: conc_prim_pert(reaction%naqcomp)
  PetscReal :: sec_jac_num(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: conc_current_M_pert(reaction%naqcomp)
  PetscReal :: total_current_M_pert(reaction%naqcomp)
  PetscReal :: res_transport_pert(reaction%naqcomp)
  PetscReal :: total_primary_node_pert(reaction%naqcomp)
  PetscReal :: dtotal_prim_num(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: dPsisec_dCprim_num(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: pert
  PetscReal :: coeff_diag_dm(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_left_dm(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right_dm(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_left_pert(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_diag_pert(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right_pert(reaction%naqcomp,reaction%naqcomp, &
                           sec_transport_vars%ncells)
  PetscReal :: coeff_left_copy(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_diag_copy(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right_copy(reaction%naqcomp,reaction%naqcomp, &
                           sec_transport_vars%ncells)

  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: lnQK, tempreal
  PetscReal :: J_up, J_dn

  PetscReal :: total_sorb_upd(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: total_sorb_prev(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: dtotal_sorb_upd(reaction%naqcomp,reaction%naqcomp,sec_transport_vars%ncells)

  class(material_auxvar_type), allocatable :: material_auxvar

  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  vol = sec_transport_vars%vol
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus
  area_fm = sec_transport_vars%interfacial_area
  ncomp = reaction%naqcomp

  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  material_auxvar%porosity = porosity

  do j = 1, ncomp
    do i = 1, ngcells
      total_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total(j,1)
      if (reaction%neqsorb > 0) then
        total_sorb_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total_sorb_eq(j)
      endif
    enddo
  enddo
  conc_upd = sec_transport_vars%updated_conc


  ! Note that sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j) units are in mol/kg
  ! Need to convert to mol/L since the units of total. in the Thomas
  ! algorithm are in mol/L

  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  res = 0.d0
  res_transport = 0.d0
  res_transport_pert = 0.d0
  dC_prim = 0.0
  rhs = 0.d0
  D_M = 0.d0
  identity = 0.d0
  b_M = 0.d0
  inv_D_M = 0.d0
  total_current_M = 0.d0
  dPsisec_dCprim = 0.d0
  dCsec_dCprim = 0.d0

  total_primary_node = auxvar%total(:,1) ! in mol/L
  dtotal_prim = auxvar%aqueous%dtotal(:,:,1)

  ! Compute totals and its derivatives total_upd dtotal
  call RTAuxVarInit(rt_auxvar,reaction,option)
  do i = 1, ngcells
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i),option)
    rt_auxvar%pri_molal = conc_upd(:,i)
    mc_pri_molal(:, i) = conc_upd(:,i)

    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    if (reaction%neqsorb > 0) then
       call RTotalSorb(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                       reaction%isotherm%multicontinuum_isotherm_rxn,option)
    endif
    total_upd(:,i) = rt_auxvar%total(:,1)           ! phase 1 liquid
    dtotal(:,:,i) = rt_auxvar%aqueous%dtotal(:,:,1) ! phase 1 liquid
    mc_sec_molal(:, i) = rt_auxvar%sec_molal(:)
    if (reaction%neqsorb > 0) then
      total_sorb_upd(:,i) = rt_auxvar%total_sorb_eq(:)
      dtotal_sorb_upd(:,:,i) = rt_auxvar%dtotal_sorb_eq(:,:)
    endif
  enddo

!================ Calculate the secondary residual =============================

  pordt = porosity/option%tran_dt*1d3

  do i = 1, ngcells

    ! D_j \nabla C_j
    do icomp = 1, ncomp
      n = icomp + (i-1)*ncomp

      pordiff = porosity*diffusion_coefficient*tortuosity*global_auxvar%den_kg(1)
      pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(icomp)*&
              tortuosity*global_auxvar%den_kg(1)

      ! Accumulation
      res(n) = pordt*(total_upd(icomp,i) - total_prev(icomp,i))*vol(i)    ! in mol/L*m3/s
      if (reaction%neqsorb > 0) then
        res(n) = res(n) + vol(i)/option%tran_dt*(total_sorb_upd(icomp,i) - total_sorb_prev(icomp,i))
      endif

      ! Flux terms
      if (i.gt.1.and.i.lt.ngcells) then
        res(n) = res(n) - pordiff_prim*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                        (mc_pri_molal(icomp,i+1) - &
                        mc_pri_molal(icomp,i))
        res(n) = res(n) + pordiff_prim*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                        (mc_pri_molal(icomp,i) - &
                        mc_pri_molal(icomp,i-1))

      ! Apply boundary conditions
      ! Inner boundary
      else if (i.eq.1) then
        res(n) = res(n) - pordiff_prim*area(1)/(dm_minus(2) + dm_plus(1))* &
                        (mc_pri_molal(icomp,i+1) - &
                        mc_pri_molal(icomp,i))

      ! Outer boundary -- closest to primary node
      else !if (i.eq.ngcells) then
        res(n) = res(n) - &
                 pordiff_prim*area(ngcells)/dm_plus(ngcells)* &
                 (auxvar%pri_molal(icomp) - mc_pri_molal(icomp,i))
        res(n) = res(n) + &
                 pordiff_prim*area(ngcells-1)/(dm_minus(ngcells) &
                 + dm_plus(ngcells-1))* &
                 (mc_pri_molal(icomp,i) - mc_pri_molal(icomp,i-1))
      endif

    enddo


    ! summatory D_i \nabla C_i
    do icplx = 1, reaction%neqcplx

      pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)*&
                    tortuosity*global_auxvar%den_kg(1)
      nicomp = reaction%eqcplxspecid(0,icplx)

      do ni = 1, nicomp
        icomp = reaction%eqcplxspecid(ni,icplx)
        n = icomp + (i-1)*ncomp

        ! Flux terms
        if (i.gt.1.and.i.lt.ngcells) then
          res(n) = res(n) - pordiff_sec*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                   (mc_sec_molal(icplx,i+1) - &
                   mc_sec_molal(icplx,i))* &
                   reaction%eqcplxstoich(ni,icplx)
          res(n) = res(n) + pordiff_sec*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                   (mc_sec_molal(icplx,i) - &
                   mc_sec_molal(icplx,i-1))* &
                   reaction%eqcplxstoich(ni,icplx)

        ! Apply boundary conditions
        ! Inner boundary
        else if (i.eq.1) then
          res(n) = res(n) - pordiff_sec*area(1)/(dm_minus(2) + dm_plus(1))* &
                   (mc_sec_molal(icplx,i+1)- &
                   mc_sec_molal(icplx,i))* &
                   reaction%eqcplxstoich(ni,icplx)

        ! Outer boundary -- closest to primary node
        else !if (i.eq.ngcells) then
          res(n) = res(n) - &
                   pordiff_sec*area(ngcells)/dm_plus(ngcells)* &
                   (auxvar%sec_molal(icplx) - &
                   mc_sec_molal(icplx,i))* &
                   reaction%eqcplxstoich(ni,icplx)
          res(n) = res(n) + pordiff_sec* &
                   area(ngcells-1)/(dm_minus(ngcells) + dm_plus(ngcells-1))* &
                   (mc_sec_molal(icplx,i) - &
                   mc_sec_molal(icplx,i-1))* &
                   reaction%eqcplxstoich(ni,icplx)
        endif
      enddo
    enddo

    ! Electromigration term
    call ComputeElectricPotentialTotalComponent_NP(i,reaction, &
                                            rt_parameter,  &
                                            rt_auxvar, &
                                            sec_transport_vars, &
                                            icomp,auxvar, &
                                            potent_dn, &
                                            potent_up, &
                                            mc_pri_molal, &
                                            mc_sec_molal)

    sum_denom_up = 1d-40
    sum_denom_dn = 1d-40

    do jcomp = 1, rt_parameter%naqcomp
      curCharge = reaction%primary_spec_Z(jcomp)
      sum_denom_up = sum_denom_up + potent_up(jcomp)*curCharge
      sum_denom_dn = sum_denom_dn + potent_dn(jcomp)*curCharge
    enddo

    do jcomp = 1, rt_parameter%naqcomp
      factor_up_em(jcomp) = potent_up(jcomp) / sum_denom_up
      factor_dn_em(jcomp) = potent_dn(jcomp) / sum_denom_dn
    enddo

    portort = porosity*tortuosity*global_auxvar%den_kg(1)

    sum_transp_dn = 0d0
    sum_transp_up = 0d0

    ! Primary l index
    do lcomp = 1, rt_parameter%naqcomp
      curCharge = reaction%primary_spec_Z(lcomp)
      cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)

      ! Flux terms
      if (i.gt.1.and.i.lt.ngcells) then
        curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                     (mc_pri_molal(lcomp, i+1) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                     (mc_pri_molal(lcomp, i) - &
                     mc_pri_molal(lcomp, i-1))

      ! Apply boundary conditions
      ! Inner boundary
      else if (i.eq.1) then
        curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                     (mc_pri_molal(lcomp, i+1) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = 0.d0

      ! Outer boundary
      else !if (i.eq.ngcells) then
        curConc_up = area(ngcells)/dm_plus(ngcells)* &
                     (auxvar%pri_molal(lcomp) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                     + dm_plus(ngcells-1))* &
                     (mc_pri_molal(lcomp, i) - &
                     mc_pri_molal(lcomp, i-1))
      endif

      sum_transp_up = sum_transp_up - curCharge*cur_diff*(curConc_up)
      sum_transp_dn = sum_transp_dn + curCharge*cur_diff*(curConc_dn)

      ! Secondary term k
      do icplx = 1, reaction%neqcplx
        nicomp = reaction%eqcplxspecid(0,icplx)
        do ni = 1, nicomp
          kcomp = reaction%eqcplxspecid(ni,icplx)
          if (kcomp == lcomp) then
            curCharge = reaction%primary_spec_Z(kcomp)
            cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
            curStoich = reaction%eqcplxstoich(ni,icplx)

            ! Flux terms
            if (i.gt.1.and.i.lt.ngcells) then
              curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                          (mc_sec_molal(icplx, i+1) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                          (mc_sec_molal(icplx, i) - &
                          mc_sec_molal(icplx, i-1))

            ! Apply boundary conditions
            ! Inner boundary
            else if (i.eq.1) then
              curConc_up = area(1)/(dm_minus(2) + dm_plus(1))* &
                          (mc_sec_molal(icplx, i+1) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = 0.d0

            ! Outer boundary
            else !if (i.eq.ngcells) then
              curConc_up = area(ngcells)/dm_plus(ngcells)* &
                          (auxvar%sec_molal(icplx) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                          + dm_plus(ngcells-1))* &
                          (mc_sec_molal(icplx, i) - &
                          mc_sec_molal(icplx, i-1))
            endif

            sum_transp_up = sum_transp_up - curCharge*cur_diff*curStoich*(curConc_up)
            sum_transp_dn = sum_transp_dn + curCharge*cur_diff*curStoich*(curConc_dn)
          endif

        enddo
      enddo

    enddo

    do icomp = 1, rt_parameter%naqcomp
      n = icomp + (i-1)*ncomp
      res(n) = res(n) - factor_dn_em(icomp)*portort*sum_transp_dn
      res(n) = res(n) - factor_up_em(icomp)*portort*sum_transp_up
    enddo

  enddo

  !res = res*1.d3 ! Convert mol/L*m3/s to mol/s

!================ Calculate the secondary jacobian =============================

  do i = 1, ngcells

    ln_conc = log(rt_auxvar%pri_molal)
    ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
    lnQK = 0d0

    ! Accumulation
    do j = 1, ncomp
      do k = 1, ncomp
        coeff_diag(j,k,i) = coeff_diag(j,k,i) + pordt*vol(i)*dtotal(j,k,i)
        if (reaction%neqsorb > 0) then
          coeff_diag(j,k,i) = coeff_diag(j,k,i) + vol(i)/option%tran_dt*(dtotal_sorb_upd(j,k,i))
        endif
      enddo
    enddo


    ! d (D_j * C_j) / dC_j =  D_j
    !do j = 1, rt_parameter%naqcomp
    do j = 1, ncomp
      do k = 1, ncomp
        pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(j)*tortuosity &
                  *global_auxvar%den_kg(1)

        ! Flux terms
        if (i.gt.1.and.i.lt.ngcells) then
          coeff_diag(j,k,i) = coeff_diag(j,k,i) + &
                              pordiff_prim*area(i)/(dm_minus(i+1) + dm_plus(i)) + &
                              pordiff_prim*area(i-1)/(dm_minus(i) + dm_plus(i-1))
          coeff_left(j,k,i) = coeff_left(j,k,i) - &
                              pordiff_prim*area(i-1)/(dm_minus(i) + dm_plus(i-1))
          coeff_right(j,k,i) = coeff_right(j,k,i) - &
                              pordiff_prim*area(i)/(dm_minus(i+1) + dm_plus(i))

        ! Apply boundary conditions
        ! Inner boundary
        else if (i.eq.1) then
          coeff_diag(j,k,1) = coeff_diag(j,k,1) + &
                              pordiff_prim*area(1)/(dm_minus(2) + dm_plus(1))
          coeff_right(j,k,i) = coeff_right(j,k,1) - &
                              pordiff_prim*area(1)/(dm_minus(2) + dm_plus(1))

        ! Outer boundary -- closest to primary node
        else !if (i.eq.ngcells) then
          coeff_diag(j,k,ngcells) = coeff_diag(j,k,ngcells) + &
                                    pordiff_prim*area(ngcells-1)/(dm_minus(ngcells) &
                                    + dm_plus(ngcells-1)) &
                                    + pordiff_prim*area(ngcells)/dm_plus(ngcells)
          coeff_left(j,k,ngcells) = coeff_left(j,k,ngcells) - &
                                    pordiff_prim*area(ngcells-1)/(dm_minus(ngcells) &
                                    + dm_plus(ngcells-1))
        endif
      enddo
    enddo


    ! d (sumatory ( frac * D_i * C_i )) / dC_j =
    ! sumatory ( frac * D_i * dC_i/dC_j )
    do icplx = 1, reaction%neqcplx

      lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

      ! activity of water
      if (reaction%eqcplxh2oid(icplx) > 0) then
        lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
      endif

      nicomp = reaction%eqcplxspecid(0,icplx)
      do ni = 1, nicomp
        icomp = reaction%eqcplxspecid(ni,icplx)
        lnQK = lnQK + reaction%eqcplxstoich(ni,icplx)*ln_act(icomp)
      enddo

      pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)* &
                    tortuosity*global_auxvar%den_kg(1)

      nicomp = reaction%eqcplxspecid(0,icplx)
      do j = 1, nicomp
        jcomp = reaction%eqcplxspecid(j,icplx)
        tempreal = reaction%eqcplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                 rt_auxvar%sec_act_coef(icplx)

        do k = 1, nicomp
          icomp = reaction%eqcplxspecid(k,icplx)

          ! Flux terms
          if (i.gt.1.and.i.lt.ngcells) then

            J_up = pordiff_sec*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                   reaction%eqcplxstoich(k,icplx)*tempreal
            J_dn = pordiff_sec*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                   reaction%eqcplxstoich(k,icplx)*tempreal

            coeff_diag(icomp,jcomp,i) = coeff_diag(icomp,jcomp,i) + &
                                        J_up + J_dn
            coeff_left(icomp,jcomp,i) = coeff_left(icomp,jcomp,i) - J_dn
            coeff_right(icomp,jcomp,i) = coeff_right(icomp,jcomp,i) - J_up

          ! Apply boundary conditions
          ! Inner boundary
          else if (i.eq.1) then

            J_up = pordiff_sec*area(1)/(dm_minus(2) + dm_plus(1))* &
                   reaction%eqcplxstoich(k,icplx)*tempreal

            coeff_diag(j,k,1) = coeff_diag(j,k,1) + J_up
            coeff_right(j,k,1) = coeff_right(j,k,1) - J_up

          ! Outer boundary -- closest to primary node
          else !if (i.eq.ngcells) then

            J_up = pordiff_sec*area(ngcells)/dm_plus(ngcells)* &
                   reaction%eqcplxstoich(k,icplx)*tempreal
            J_dn = pordiff_sec*area(ngcells-1)/(dm_minus(ngcells) + dm_plus(ngcells-1))* &
                   reaction%eqcplxstoich(k,icplx)*tempreal

            coeff_diag(j,k,ngcells) = coeff_diag(j,k,ngcells) + &
                                      J_up + J_dn
            coeff_left(j,k,ngcells) = coeff_left(j,k,ngcells) - &
                                      J_dn
          endif

        enddo
      enddo
    enddo
  enddo


!====================== Add reaction contributions =============================

  ! Reaction
  do i = 1, ngcells
    res_react = 0.d0
    jac_react = 0.d0
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i), &
                      option)
    rt_auxvar%pri_molal = conc_upd(:,i) ! in mol/kg
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    material_auxvar%volume = vol(i)
    call RReaction(res_react,jac_react,PETSC_TRUE, &
                   rt_auxvar,global_auxvar,material_auxvar,reaction,option)
    do j = 1, ncomp
      res(j+(i-1)*ncomp) = res(j+(i-1)*ncomp) + res_react(j)
    enddo
    coeff_diag(:,:,i) = coeff_diag(:,:,i) + jac_react  ! in kg water/s
  enddo
  call MaterialAuxVarStrip(material_auxvar)
  deallocate(material_auxvar)

!============================== Forward solve ==================================

  rhs = -res

  if (reaction%use_log_formulation) then
  ! scale the jacobian by concentrations
    i = 1
    do k = 1, ncomp
      coeff_diag(:,k,i) = coeff_diag(:,k,i)*conc_upd(k,i) ! m3/s*kg/L
      coeff_right(:,k,i) = coeff_right(:,k,i)*conc_upd(k,i+1)
    enddo
    do i = 2, ngcells-1
      do k = 1, ncomp
        coeff_diag(:,k,i) = coeff_diag(:,k,i)*conc_upd(k,i) ! m3/s*kg/L
        coeff_left(:,k,i) = coeff_left(:,k,i)*conc_upd(k,i-1)
        coeff_right(:,k,i) = coeff_right(:,k,i)*conc_upd(k,i+1)
      enddo
    enddo
    i = ngcells
      do k = 1, ncomp
        coeff_diag(:,k,i) = coeff_diag(:,k,i)*conc_upd(k,i) ! m3/s*kg/L
        coeff_left(:,k,i) = coeff_left(:,k,i)*conc_upd(k,i-1)
      enddo
  endif

  ! First do an LU decomposition for calculating D_M matrix
  coeff_diag_dm = coeff_diag
  coeff_left_dm = coeff_left
  coeff_right_dm = coeff_right

  select case (option%secondary_continuum_solver)
    case(1)
      do i = 2, ngcells
        coeff_left_dm(:,:,i-1) = coeff_left_dm(:,:,i)
      enddo
      coeff_left_dm(:,:,ngcells) = 0.d0
      call bl3dfac(ngcells,ncomp,coeff_right_dm,coeff_diag_dm,coeff_left_dm,pivot)
    case(2)
      call decbt(ncomp,ngcells,ncomp,coeff_diag_dm,coeff_right_dm,coeff_left_dm,pivot,ier)
      if (ier /= 0) then
        print *,'error in matrix decbt: ier = ',ier
        stop
      endif
    case(3)
      ! Thomas algorithm for tridiagonal system
      ! Forward elimination
      if (ncomp /= 1) then
        option%io_buffer = 'THOMAS algorithm can be used only with single '// &
                           'component chemistry'
        call PrintErrMsg(option)
      endif
      do i = 2, ngcells
        m = coeff_left_dm(ncomp,ncomp,i)/coeff_diag_dm(ncomp,ncomp,i-1)
        coeff_diag_dm(ncomp,ncomp,i) = coeff_diag_dm(ncomp,ncomp,i) - &
                                    m*coeff_right_dm(ncomp,ncomp,i-1)
      enddo
    case default
      option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                         'HINDMARSH or KEARST. For single component'// &
                         'chemistry THOMAS can be used.'
      call PrintErrMsg(option)
  end select

  ! Set the values of D_M matrix and create identity matrix of size ncomp x ncomp
  do i = 1, ncomp
    do j = 1, ncomp
      D_M(i,j) = coeff_diag_dm(i,j,ngcells)
      if (j == i) then
        identity(i,j) = 1.d0
      else
        identity(i,j) = 0.d0
      endif
    enddo
  enddo

  ! Find the inverse of D_M
  call LUDecomposition(D_M,ncomp,indx,d)
  do j = 1, ncomp
    call LUBackSubstitution(D_M,ncomp,indx,identity(1,j))
  enddo
  inv_D_M = identity

  if (option%numerical_derivatives_multi_coupling) then
    ! Store the coeffs for numerical jacobian
    coeff_diag_copy = coeff_diag
    coeff_left_copy = coeff_left
    coeff_right_copy = coeff_right
  endif

  select case (option%secondary_continuum_solver)
    case(1)
      do i = 2, ngcells
        coeff_left(:,:,i-1) = coeff_left(:,:,i)
      enddo
      coeff_left(:,:,ngcells) = 0.d0
      call bl3dfac(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot)
      call bl3dsolf(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot, &
                    ONE_INTEGER,rhs)
    case(2)
      call decbt(ncomp,ngcells,ncomp,coeff_diag,coeff_right,coeff_left, &
                 pivot,ier)
      if (ier /= 0) then
        print *,'error in matrix decbt: ier = ',ier
        stop
      endif
      call solbtf(ncomp,ngcells,ncomp,coeff_diag,coeff_right,coeff_left, &
                  pivot,rhs)
    case(3)
      ! Thomas algorithm for tridiagonal system
      ! Forward elimination
      if (ncomp /= 1) then
        option%io_buffer = 'THOMAS algorithm can be used only with single '// &
                           'component chemistry'
        call PrintErrMsg(option)
      endif
      do i = 2, ngcells
        m = coeff_left(ncomp,ncomp,i)/coeff_diag(ncomp,ncomp,i-1)
        coeff_diag(ncomp,ncomp,i) = coeff_diag(ncomp,ncomp,i) - &
                                    m*coeff_right(ncomp,ncomp,i-1)
        rhs(i) = rhs(i) - m*rhs(i-1)
      enddo
      rhs(ngcells) = rhs(ngcells)/coeff_diag(ncomp,ncomp,ngcells)
    case default
      option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                         'HINDMARSH or KEARST. For single component'// &
                         'chemistry THOMAS can be used.'
      call PrintErrMsg(option)
  end select

  ! Update the secondary concentrations
  do i = 1, ncomp
    if (reaction%use_log_formulation) then
      ! convert log concentration to concentration
      rhs(i+(ngcells-1)*ncomp) = dsign(1.d0,rhs(i+(ngcells-1)*ncomp))* &
        min(dabs(rhs(i+(ngcells-1)*ncomp)),reaction%max_dlnC)
      conc_current_M(i) = conc_upd(i,ngcells)*exp(rhs(i+(ngcells-1)*ncomp))
    else
      conc_current_M(i) = conc_upd(i,ngcells) + rhs(i+(ngcells-1)*ncomp)
    endif
  enddo

  ! Update the secondary continuum totals at the outer matrix node
  call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(ngcells), &
                    option)
  rt_auxvar%pri_molal = conc_current_M ! in mol/kg
  call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
  total_current_M = rt_auxvar%total(:,1)
  if (reaction%neqcplx > 0) sec_sec_molal_M = rt_auxvar%sec_molal
  call RTAuxVarStrip(rt_auxvar)



  b_m = porosity*tortuosity/dm_plus(ngcells)*area(ngcells)*inv_D_M ! in m·s/mol

  dCsec_dCprim = b_m*dtotal_prim

  ! Calculate the dervative of outer matrix node total with respect to the
  ! primary node concentration
  if (reaction%use_log_formulation) then ! log formulation
    do j = 1, ncomp
      do l = 1, ncomp
        dPsisec_dCprim(j,l) = dCsec_dCprim(j,l)*conc_current_M(j)
      enddo
    enddo

    if (reaction%neqcplx > 0) then
      do icplx = 1, reaction%neqcplx
        ncompeq = reaction%eqcplxspecid(0,icplx)
        do j = 1, ncompeq
          jcomp = reaction%eqcplxspecid(j,icplx)
          do l = 1, ncompeq
            lcomp = reaction%eqcplxspecid(l,icplx)
            do k = 1, ncompeq
              kcomp = reaction%eqcplxspecid(k,icplx)
              dPsisec_dCprim(jcomp,lcomp) = dPsisec_dCprim(jcomp,lcomp) + &
                                            reaction%eqcplxstoich(j,icplx)* &
                                            reaction%eqcplxstoich(k,icplx)* &
                                            dCsec_dCprim(kcomp,lcomp)* &
                                            sec_sec_molal_M(icplx)
            enddo
          enddo
        enddo
      enddo
    endif

  else   ! linear case

    dPsisec_dCprim = dCsec_dCprim

    if (reaction%neqcplx > 0) then
      do icplx = 1, reaction%neqcplx
        ncompeq = reaction%eqcplxspecid(0,icplx)
        do j = 1, ncompeq
          jcomp = reaction%eqcplxspecid(j,icplx)
          do l = 1, ncompeq
            lcomp = reaction%eqcplxspecid(l,icplx)
            do k = 1, ncompeq
              kcomp = reaction%eqcplxspecid(k,icplx)

              dPsisec_dCprim(jcomp,lcomp) = dPsisec_dCprim(jcomp,lcomp) + &
                                            reaction%eqcplxstoich(j,icplx)* &
                                            reaction%eqcplxstoich(k,icplx)* &
                                            dCsec_dCprim(kcomp,lcomp)* &
                                            sec_sec_molal_M(icplx)/ &
                                            conc_current_M(kcomp)
            enddo
          enddo
        enddo
      enddo
    endif

  endif

  dPsisec_dCprim = dPsisec_dCprim*global_auxvar%den_kg(1)*1.d-3 ! in kg/L

  ! D_j \nabla C_j
  do j = 1, ncomp
    pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(j)* &
                   tortuosity*global_auxvar%den_kg(1)
    res_transport(j) = pordiff_prim*area_fm/dm_plus(ngcells)* &
                       (conc_current_M(j) - &
                       auxvar%pri_molal(j))*prim_vol ! in mol/s

  enddo

  ! summatory D_i \nabla C_i
  do icplx = 1, reaction%neqcplx
    pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)* &
                  tortuosity*global_auxvar%den_kg(1)
    nicomp = reaction%eqcplxspecid(0,icplx)
    do ni = 1, nicomp
      icomp = reaction%eqcplxspecid(ni,icplx)
      res_transport(icomp) = res_transport(icomp) + &
                       pordiff_sec*area_fm/dm_plus(ngcells)* &
                       (sec_sec_molal_M(icplx) - &
                       auxvar%sec_molal(icplx))*prim_vol* &
                       reaction%eqcplxstoich(ni,icplx) ! in mol/s
    enddo
  enddo

  ! Electromigration term
  call ComputeElectricPotentialTotalComponent_NP(ngcells,reaction, &
                                          rt_parameter,  &
                                          rt_auxvar, &
                                          sec_transport_vars, &
                                          icomp,auxvar, &
                                          potent_dn, &
                                          potent_up, &
                                          mc_pri_molal, &
                                          mc_sec_molal)

  sum_denom_up = 1d-40

  do jcomp = 1, rt_parameter%naqcomp
    curCharge = reaction%primary_spec_Z(jcomp)
    sum_denom_up = sum_denom_up + potent_up(jcomp)*curCharge
  enddo

  do jcomp = 1, rt_parameter%naqcomp
    factor_up_em(jcomp) = potent_up(jcomp) / sum_denom_up
  enddo

  portort = porosity*tortuosity*global_auxvar%den_kg(1)

  sum_transp_up = 0d0

  ! Primary l index
  do lcomp = 1, rt_parameter%naqcomp
    curCharge = reaction%primary_spec_Z(lcomp)
    cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)
    curConc_up = area_fm/dm_plus(ngcells)*prim_vol* &
              (auxvar%pri_molal(lcomp)-conc_current_M(lcomp))
    sum_transp_up = sum_transp_up - curCharge*cur_diff*curConc_up

    ! Secondary term k
    do icplx = 1, reaction%neqcplx
      nicomp = reaction%eqcplxspecid(0,icplx)
      do ni = 1, nicomp
        kcomp = reaction%eqcplxspecid(ni,icplx)
        if (kcomp == lcomp) then
          curCharge = reaction%primary_spec_Z(kcomp)
          cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
          curStoich = reaction%eqcplxstoich(ni,icplx)
          curConc_up = area_fm/dm_plus(ngcells)*prim_vol* &
                    (auxvar%sec_molal(icplx) - sec_sec_molal_M(icplx))
          sum_transp_up = sum_transp_up - curCharge*cur_diff*curStoich*curConc_up
        endif

      enddo
    enddo
  enddo

  do icomp = 1, rt_parameter%naqcomp
    res_transport(icomp) = res_transport(icomp) - portort*factor_up_em(icomp)*sum_transp_up ! in mol/s
  enddo

  ! Calculate the jacobian contribution due to coupling term
  sec_jac = area_fm*pordiff/dm_plus(ngcells)*(dPsisec_dCprim - dtotal_prim)* &
            prim_vol !*1.d3 ! in kg water/s

  ! Store the contribution to the primary jacobian term
  sec_transport_vars%sec_jac = sec_jac
  sec_transport_vars%sec_jac_update = PETSC_TRUE

  ! Store the coefficients from LU decomposition of the block tridiagonal
  ! sytem. These will be called later to perform backsolve to the get the
  ! updated secondary continuum concentrations at the end of the timestep
  sec_transport_vars%cxm = coeff_left
  sec_transport_vars%cxp = coeff_right
  sec_transport_vars%cdl = coeff_diag

  ! Store the solution of the forward solve
  sec_transport_vars%r = rhs


!============== Numerical jacobian for coupling term ===========================


  if (option%numerical_derivatives_multi_coupling) then

    call RTAuxVarInit(rt_auxvar,reaction,option)
    conc_prim = auxvar%pri_molal
    conc_prim_pert = conc_prim

    do l = 1, ncomp

      conc_prim_pert = conc_prim
      pert = conc_prim(l)*perturbation_tolerance
      conc_prim_pert(l) = conc_prim_pert(l) + pert

      res = 0.d0
      rhs = 0.d0

      coeff_diag_pert = coeff_diag_copy
      coeff_left_pert = coeff_left_copy
      coeff_right_pert = coeff_right_copy

      call RTAuxVarCopy(rt_auxvar,auxvar,option)
      rt_auxvar%pri_molal = conc_prim_pert ! in mol/kg
      call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
      total_primary_node_pert = rt_auxvar%total(:,1)

!================ Calculate the secondary residual =============================
      do i = 1, ngcells

        ! D_j \nabla C_j
        do icomp = 1, ncomp

          pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(icomp)*&
                         tortuosity*global_auxvar%den_kg(1)

          n = icomp + (i-1)*ncomp

          ! Accumulation
          res(n) = pordt*(total_upd(icomp,i) - total_prev(icomp,i))*vol(i)    ! in mol/L*m3/s
          if (reaction%neqsorb > 0) then
            res(n) = res(n) + vol(i)/option%tran_dt*(total_sorb_upd(icomp,i) - total_sorb_prev(icomp,i))
          endif

          ! Flux terms
          if (i.gt.1.and.i.lt.ngcells) then
            res(n) = res(n) - pordiff_prim*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                            (sec_transport_vars%sec_rt_auxvar(i+1)%pri_molal(icomp) - &
                            sec_transport_vars%sec_rt_auxvar(i)%pri_molal(icomp))
            res(n) = res(n) + pordiff_prim*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                            (sec_transport_vars%sec_rt_auxvar(i)%pri_molal(icomp) - &
                            sec_transport_vars%sec_rt_auxvar(i-1)%pri_molal(icomp))

          ! Apply boundary conditions
          ! Inner boundary
          else if (i.eq.1) then
            res(n) = res(n) - pordiff_prim*area(1)/(dm_minus(2) + dm_plus(1))* &
                            (sec_transport_vars%sec_rt_auxvar(2)%pri_molal(icomp) - &
                            sec_transport_vars%sec_rt_auxvar(1)%pri_molal(icomp))

          ! Outer boundary -- closest to primary node
          else !if (i.eq.ngcells) then
            res(n) = res(n) - &
                     pordiff_prim*area(ngcells)/dm_plus(ngcells)* &
                     (auxvar%pri_molal(icomp) - &
                     sec_transport_vars%sec_rt_auxvar(ngcells)%pri_molal(icomp))
            res(n) = res(n) + &
                     pordiff_prim*area(ngcells-1)/(dm_minus(ngcells) &
                     + dm_plus(ngcells-1))* &
                     (sec_transport_vars%sec_rt_auxvar(ngcells)%pri_molal(icomp) - &
                     sec_transport_vars%sec_rt_auxvar(ngcells-1)%pri_molal(icomp))
          endif

        enddo


        ! summatory D_i \nabla C_i
        do icplx = 1, reaction%neqcplx

          pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)*&
                        tortuosity*global_auxvar%den_kg(1)
          nicomp = reaction%eqcplxspecid(0,icplx)

          do ni = 1, nicomp
            icomp = reaction%eqcplxspecid(ni,icplx)

            n = icomp + (i-1)*ncomp

            ! Flux terms
            if (i.gt.1.and.i.lt.ngcells) then
              res(n) = res(n) - pordiff_sec*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                       (sec_transport_vars%sec_rt_auxvar(i+1)%sec_molal(icplx) - &
                       sec_transport_vars%sec_rt_auxvar(i)%sec_molal(icplx))* &
                       reaction%eqcplxstoich(ni,icplx)
              res(n) = res(n) + pordiff_sec*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                       (sec_transport_vars%sec_rt_auxvar(i)%sec_molal(icplx) - &
                       sec_transport_vars%sec_rt_auxvar(i-1)%sec_molal(icplx))* &
                       reaction%eqcplxstoich(ni,icplx)

            ! Apply boundary conditions
            ! Inner boundary
            else if (i.eq.1) then
              res(n) = res(n) - pordiff_sec*area(1)/(dm_minus(2) + dm_plus(1))* &
                       (sec_transport_vars%sec_rt_auxvar(2)%sec_molal(icplx) - &
                       sec_transport_vars%sec_rt_auxvar(1)%sec_molal(icplx))* &
                       reaction%eqcplxstoich(ni,icplx)

            ! Outer boundary -- closest to primary node
            else !if (i.eq.ngcells) then
              res(n) = res(n) - &
                       pordiff_sec*area(ngcells)/dm_plus(ngcells)* &
                       (auxvar%sec_molal(icplx) - &
                       sec_transport_vars%sec_rt_auxvar(ngcells)%sec_molal(icplx))* &
                       reaction%eqcplxstoich(ni,icplx)
              res(n) = res(n) + pordiff_sec* &
                       area(ngcells-1)/(dm_minus(ngcells) + dm_plus(ngcells-1))* &
                       (sec_transport_vars%sec_rt_auxvar(ngcells)%sec_molal(icplx) - &
                       sec_transport_vars%sec_rt_auxvar(ngcells-1)%sec_molal(icplx))* &
                       reaction%eqcplxstoich(ni,icplx)
            endif
          enddo
        enddo

        ! Electromigration term
        call ComputeElectricPotentialTotalComponent_NP(i,reaction, &
                                                rt_parameter,  &
                                                rt_auxvar, &
                                                sec_transport_vars, &
                                                icomp,auxvar, &
                                                potent_dn, &
                                                potent_up, &
                                                mc_pri_molal, &
                                                mc_sec_molal)

        sum_denom_up = 1d-40
        sum_denom_dn = 1d-40

        do jcomp = 1, rt_parameter%naqcomp
          curCharge = reaction%primary_spec_Z(jcomp)
          sum_denom_dn = sum_denom_dn + potent_dn(jcomp)*curCharge
          sum_denom_up = sum_denom_up + potent_up(jcomp)*curCharge
        enddo

        do jcomp = 1, rt_parameter%naqcomp
          factor_up_em(jcomp) = potent_up(jcomp) / sum_denom_up
          factor_dn_em(jcomp) = potent_dn(jcomp) / sum_denom_dn
        enddo

        portort = porosity*tortuosity* &
        global_auxvar%den_kg(1)

        ! Primary l index
        do lcomp = 1, rt_parameter%naqcomp
          curCharge = reaction%primary_spec_Z(lcomp)
          cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)

          ! Flux terms
          if (i.gt.1.and.i.lt.ngcells) then
            curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                         (sec_transport_vars%sec_rt_auxvar(i+1)%pri_molal(lcomp) - &
                         sec_transport_vars%sec_rt_auxvar(i)%pri_molal(lcomp))
            curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                         (sec_transport_vars%sec_rt_auxvar(i)%pri_molal(lcomp) - &
                         sec_transport_vars%sec_rt_auxvar(i-1)%pri_molal(lcomp))

          ! Apply boundary conditions
          ! Inner boundary
          else if (i.eq.1) then
            curConc_up = area(1)/(dm_minus(2) + dm_plus(1))* &
                         (sec_transport_vars%sec_rt_auxvar(2)%pri_molal(lcomp) - &
                         sec_transport_vars%sec_rt_auxvar(1)%pri_molal(lcomp))
            curConc_dn = 0.d0

          ! Outer boundary
          else !if (i.eq.ngcells) then
            curConc_up = area(ngcells)/dm_plus(ngcells)* &
                         (auxvar%pri_molal(lcomp) - &
                         sec_transport_vars%sec_rt_auxvar(ngcells)%pri_molal(lcomp))
            curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                         + dm_plus(ngcells-1))* &
                         (sec_transport_vars%sec_rt_auxvar(ngcells)%pri_molal(lcomp) - &
                         sec_transport_vars%sec_rt_auxvar(ngcells-1)%pri_molal(lcomp))
          endif

          sum_transp_up = sum_transp_up - curCharge*cur_diff*(curConc_up)
          sum_transp_dn = sum_transp_dn + curCharge*cur_diff*(curConc_dn)

          ! Secondary term k
          do icplx = 1, reaction%neqcplx
            nicomp = reaction%eqcplxspecid(0,icplx)
            do ni = 1, nicomp
              kcomp = reaction%eqcplxspecid(ni,icplx)
              if (kcomp == lcomp) then
                curCharge = reaction%primary_spec_Z(kcomp)
                cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
                curStoich = reaction%eqcplxstoich(ni,icplx)

                ! Flux terms
                if (i.gt.1.and.i.lt.ngcells) then
                  curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                              (sec_transport_vars%sec_rt_auxvar(i+1)%sec_molal(icplx) - &
                              sec_transport_vars%sec_rt_auxvar(i)%sec_molal(icplx))
                  curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                              (sec_transport_vars%sec_rt_auxvar(i)%sec_molal(icplx) - &
                              sec_transport_vars%sec_rt_auxvar(i-1)%sec_molal(icplx))

                ! Apply boundary conditions
                ! Inner boundary
                else if (i.eq.1) then
                  curConc_up = area(1)/(dm_minus(2) + dm_plus(1))* &
                              (sec_transport_vars%sec_rt_auxvar(2)%sec_molal(icplx) - &
                              sec_transport_vars%sec_rt_auxvar(1)%sec_molal(icplx))
                  curConc_dn = 0.d0

                ! Outer boundary
                else !if (i.eq.ngcells) then
                  curConc_up = area(ngcells)/dm_plus(ngcells)* &
                              (auxvar%sec_molal(icplx) - &
                              sec_transport_vars%sec_rt_auxvar(ngcells)%sec_molal(icplx))
                  curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                              + dm_plus(ngcells-1))* &
                              (sec_transport_vars%sec_rt_auxvar(ngcells)%sec_molal(icplx) - &
                              sec_transport_vars%sec_rt_auxvar(ngcells-1)%sec_molal(icplx))
                endif

                sum_transp_up = sum_transp_up - curCharge*cur_diff*curStoich*(curConc_up)
                sum_transp_dn = sum_transp_dn + curCharge*cur_diff*curStoich*(curConc_dn)
              endif
            enddo
          enddo
        enddo

        do icomp = 1, rt_parameter%naqcomp
          n = icomp + (i-1)*ncomp
          res(n) = res(n) - portort*factor_up_em(icomp)*sum_transp_up
          res(n) = res(n) - portort*factor_dn_em(icomp)*sum_transp_dn
        enddo

      enddo

!============================== Forward solve ==================================

    rhs = -res

    select case (option%secondary_continuum_solver)
      case(1)
        call bl3dfac(ngcells,ncomp,coeff_right_pert,coeff_diag_pert, &
                      coeff_left_pert,pivot)
        call bl3dsolf(ngcells,ncomp,coeff_right_pert,coeff_diag_pert, &
                       coeff_left_pert,pivot,ONE_INTEGER,rhs)
      case(2)
        call decbt(ncomp,ngcells,ncomp,coeff_diag_pert,coeff_right_pert, &
                    coeff_left_pert,pivot,ier)
        if (ier /= 0) then
          print *,'error in matrix decbt: ier = ',ier
          stop
        endif
        call solbtf(ncomp,ngcells,ncomp,coeff_diag_pert,coeff_right_pert, &
                     coeff_left_pert,pivot,rhs)
      case(3)
        ! Thomas algorithm for tridiagonal system
        ! Forward elimination
        if (ncomp /= 1) then
          option%io_buffer = 'THOMAS algorithm can be used only with '// &
                             'single component chemistry'
          call PrintErrMsg(option)
        endif
        do i = 2, ngcells
          m = coeff_left_pert(ncomp,ncomp,i)/coeff_diag_pert(ncomp,ncomp,i-1)
          coeff_diag_pert(ncomp,ncomp,i) = coeff_diag_pert(ncomp,ncomp,i) - &
                                      m*coeff_right_pert(ncomp,ncomp,i-1)
          rhs(i) = rhs(i) - m*rhs(i-1)
        enddo
        rhs(ngcells) = rhs(ngcells)/coeff_diag(ncomp,ncomp,ngcells)
      case default
        option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                           'HINDMARSH or KEARST. For single component'// &
                           'chemistry THOMAS can be used.'
        call PrintErrMsg(option)
      end select

      ! Update the secondary concentrations
      do i = 1, ncomp
        if (reaction%use_log_formulation) then
          ! convert log concentration to concentration
          rhs(i+(ngcells-1)*ncomp) = dsign(1.d0,rhs(i+(ngcells-1)*ncomp))* &
            min(dabs(rhs(i+(ngcells-1)*ncomp)),reaction%max_dlnC)
          conc_current_M_pert(i) = conc_upd(i,ngcells)* &
                                     exp(rhs(i+(ngcells-1)*ncomp))
        else
          conc_current_M_pert(i) = conc_upd(i,ngcells) + &
                                     rhs(i+(ngcells-1)*ncomp)
        endif
      enddo

      ! Update the secondary continuum totals at the outer matrix node
      call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(ngcells), &
                        option)
      rt_auxvar%pri_molal = conc_current_M_pert ! in mol/kg
      call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
      total_current_M_pert = rt_auxvar%total(:,1)

      ! ! Calculate the coupling term
      ! D_j \nabla C_j
      do j = 1, ncomp

        pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(j)* &
                   tortuosity*global_auxvar%den_kg(1)

        res_transport_pert(j) = res_transport_pert(j) + &
                       pordiff_prim*area_fm/dm_plus(ngcells)* &
                       (conc_current_M_pert(j) - &
                       auxvar%pri_molal(j))*prim_vol ! in mol/s

      enddo

      ! summatory D_i \nabla C_i
      do icplx = 1, reaction%neqcplx

        pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)* &
                  tortuosity*global_auxvar%den_kg(1)

        nicomp = reaction%eqcplxspecid(0,icplx)
        do ni = 1, nicomp
          icomp = reaction%eqcplxspecid(ni,icplx)

          res_transport_pert(icomp) = res_transport_pert(icomp) + &
                       pordiff_sec*area_fm/dm_plus(ngcells)* &
                       (rt_auxvar%sec_molal(icplx) - &
                       auxvar%sec_molal(icplx))*prim_vol* &
                       reaction%eqcplxstoich(ni,icplx) ! in mol/s
        enddo
      enddo

      ! Electromigration term
      call ComputeElectricPotentialTotalComponent_NP(ngcells,reaction, &
                                          rt_parameter,  &
                                          rt_auxvar, &
                                          sec_transport_vars, &
                                          icomp,auxvar, &
                                          potent_dn, &
                                          potent_up, &
                                          mc_pri_molal, &
                                          mc_sec_molal)

      sum_denom_up = 1d-40

      do jcomp = 1, rt_parameter%naqcomp
        curCharge = reaction%primary_spec_Z(jcomp)
        sum_denom_up = sum_denom_up + potent_up(jcomp)*curCharge
      enddo

      do jcomp = 1, rt_parameter%naqcomp
        factor_up_em(jcomp) = potent_up(jcomp) / sum_denom_up
      enddo

      portort = porosity*tortuosity*global_auxvar%den_kg(1)

      sum_transp_up = 0d0

      ! Primary l index
      do lcomp = 1, rt_parameter%naqcomp
        curCharge = reaction%primary_spec_Z(lcomp)
        cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)

        curConc_up = area_fm/dm_plus(ngcells)*prim_vol* &
              (auxvar%pri_molal(lcomp) - conc_current_M_pert(lcomp))
        sum_transp_up = sum_transp_up - curCharge*cur_diff*curConc_up

        ! Secondary term k
        do icplx = 1, reaction%neqcplx
          nicomp = reaction%eqcplxspecid(0,icplx)
          do ni = 1, nicomp
            kcomp = reaction%eqcplxspecid(ni,icplx)
            if (kcomp == lcomp) then
              curCharge = reaction%primary_spec_Z(kcomp)
              cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
              curStoich = reaction%eqcplxstoich(ni,icplx)

              curConc_up = area_fm/dm_plus(ngcells)*prim_vol* &
                  (auxvar%sec_molal(icplx) - rt_auxvar%sec_molal(icplx))
              sum_transp_up = sum_transp_up - curCharge*cur_diff*curStoich*curConc_up
            endif
          enddo
        enddo
      enddo

      do icomp = 1, rt_parameter%naqcomp
        res_transport_pert(icomp) = res_transport_pert(icomp) - &
                                  portort*factor_up_em(icomp)*sum_transp_up ! in mol/s
      enddo

      dtotal_prim_num(:,l) = (total_primary_node_pert(:) - &
                               total_primary_node(:))/pert

      dPsisec_dCprim_num(:,l) = (total_current_M_pert(:) - &
                                  total_current_M(:))/pert

      sec_jac_num(:,l) = (res_transport_pert(:) - res_transport(:))/pert

    enddo

    call RTAuxVarStrip(rt_auxvar)
    sec_transport_vars%sec_jac = sec_jac_num

  endif


end subroutine SecondaryRTResJacMulti_NP

! ************************************************************************** !

subroutine SecondaryRTUpdateIterate_NP(snes,P0,dP,P1,dX_changed, &
                                    X1_changed,realization,ierr)
  !
  ! Checks update after the update is done
  !
  ! Author: Satish Karra, LANL
  ! Date: 02/22/13
  !

  use Realization_Subsurface_class
  use Option_module
  use Grid_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module

  implicit none

  SNES :: snes
  Vec :: P0
  Vec :: dP
  Vec :: P1
  class(realization_subsurface_type) :: realization
  ! ignore changed flag for now.
  PetscBool :: dX_changed
  PetscBool :: X1_changed

  type(reactive_transport_param_type), pointer :: rt_parameter
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(reaction_rt_type), pointer :: reaction
  PetscInt :: local_id, ghosted_id
  PetscReal :: sec_diffusion_coefficient
  PetscReal :: sec_porosity
  PetscReal :: sec_tortuosity

  PetscErrorCode :: ierr
  PetscReal :: inf_norm_sec
  PetscReal :: max_inf_norm_sec



  option => realization%option
  grid => realization%patch%grid
  rt_auxvars => realization%patch%aux%RT%auxvars
  global_auxvars => realization%patch%aux%Global%auxvars
  reaction => realization%reaction
  rt_parameter => realization%patch%aux%RT%rt_parameter
  if (option%use_sc) then
    rt_sec_transport_vars => realization%patch%aux%SC_RT%sec_transport_vars
  endif

  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  max_inf_norm_sec = 0.d0

  if (option%use_sc) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (realization%patch%imat(ghosted_id) <= 0) cycle
      sec_diffusion_coefficient = realization%patch% &
                                  material_property_array(1)%ptr% &
                                  multicontinuum%diff_coeff(1)
      sec_porosity = realization%patch%material_property_array(1)%ptr% &
                    multicontinuum%porosity

      sec_tortuosity = realization%patch%material_property_array(1)%ptr% &
                    multicontinuum%tortuosity

      call SecondaryRTAuxVarComputeMulti(&
                                    rt_sec_transport_vars(ghosted_id), &
                                    reaction, &
                                    option)

      call SecondaryRTCheckResidual_np(rt_sec_transport_vars(ghosted_id), &
                                    rt_auxvars(ghosted_id), &
                                    global_auxvars(ghosted_id), &
                                    reaction,rt_parameter,sec_diffusion_coefficient, &
                                    sec_porosity, sec_tortuosity, &
                                    option,inf_norm_sec)

      max_inf_norm_sec = max(max_inf_norm_sec,inf_norm_sec)
    enddo
    call MPI_Allreduce(max_inf_norm_sec,option%infnorm_res_sec,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
  endif


end subroutine SecondaryRTUpdateIterate_NP

! ************************************************************************** !

subroutine SecondaryRTAuxVarComputeMulti(sec_transport_vars,reaction, &
                                         option)
  !
  ! Updates the secondary continuum
  ! concentrations at end of each time step for multicomponent system
  !
  ! Author: Satish Karra, LANL
  ! Date: 2/1/13
  !


  use Option_module
  use Reaction_Aux_module
  use Reaction_module
  use Reactive_Transport_Aux_module
  use Block_Solve_module
  use Block_Tridiag_module
  use Utility_module


  implicit none

  type(sec_transport_type) :: sec_transport_vars
  class(reaction_rt_type), pointer :: reaction
  type(option_type) :: option
  PetscReal :: coeff_left(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells)
  PetscReal :: coeff_diag(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells)
  PetscReal :: coeff_right(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells)
  PetscReal :: rhs(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: conc_upd(reaction%naqcomp,sec_transport_vars%ncells)
  PetscInt :: i, j, n
  PetscInt :: ngcells, ncomp
  PetscInt :: pivot(reaction%naqcomp,sec_transport_vars%ncells)

  ngcells = sec_transport_vars%ncells
  ncomp = reaction%naqcomp
  ! Note that sec_transport_vars%sec_conc units are in mol/kg
  ! Need to convert to mol/L since the units of conc. in the Thomas
  ! algorithm are in mol/L

  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0

  conc_upd = sec_transport_vars%updated_conc

  ! Use the stored coefficient matrices from LU decomposition of the
  ! block triagonal sytem
  coeff_left = sec_transport_vars%cxm
  coeff_right = sec_transport_vars%cxp
  coeff_diag = sec_transport_vars%cdl
  rhs = sec_transport_vars%r

  select case (option%secondary_continuum_solver)
    case(1)
      call bl3dsolb(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot, &
                    ONE_INTEGER,rhs)
    case(2)
      call solbtb(ncomp,ngcells,ncomp,coeff_diag,coeff_right,coeff_left, &
                  pivot,rhs)
    case(3)
      do i = ngcells-1, 1, -1
        rhs(i) = (rhs(i) - coeff_right(ncomp,ncomp,i)*rhs(i+1))/ &
                             coeff_diag(ncomp,ncomp,i)
      enddo
    case default
      option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                         'HINDMARSH or KEARST. For single component'// &
                         'chemistry THOMAS can be used.'
      call PrintErrMsg(option)
  end select

  do j = 1, ncomp
    do i = 1, ngcells
      n = j + (i - 1)*ncomp
      if (reaction%use_log_formulation) then
        ! convert log concentration to concentration
        rhs(n) = dsign(1.d0,rhs(n))*min(dabs(rhs(n)),reaction%max_dlnC)
        conc_upd(j,i) = exp(rhs(n))*conc_upd(j,i)
      else
        conc_upd(j,i) = rhs(n) + conc_upd(j,i)
      endif
      if (conc_upd(j,i) < 0.d0) conc_upd(j,i) = 1.d-8
    enddo
  enddo

  sec_transport_vars%updated_conc = conc_upd

end subroutine SecondaryRTAuxVarComputeMulti


! ************************************************************************** !


subroutine SecondaryRTCheckResidual_np(sec_transport_vars,auxvar, &
                                    global_auxvar, &
                                    reaction,rt_parameter, &
                                    diffusion_coefficient, &
                                    porosity, tortuosity, option, &
                                    inf_norm_sec)
  !
  ! The residual of the secondary domain are checked
  ! to ensure convergence
  !
  ! Author: Albert Nardi, Amphos 21
  ! Date: 8/31/2021
  !

  use Option_module
  use Global_Aux_module
  use Block_Solve_module
  use Block_Tridiag_module
  use Utility_module
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_module

  implicit none

  type(reactive_transport_param_type) :: rt_parameter
  type(sec_transport_type) :: sec_transport_vars
  type(reactive_transport_auxvar_type) :: auxvar
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type), pointer :: reaction
  type(option_type) :: option

  PetscReal :: res(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: curConc_up, curConc_dn
  PetscReal :: sum_transp_up, sum_transp_dn
  PetscReal :: curCharge
  PetscReal :: cur_diff
  PetscReal :: curStoich
  PetscReal :: sum_denom_up, sum_denom_dn
  PetscReal :: potent_up(reaction%naqcomp), potent_dn(reaction%naqcomp)
  PetscReal :: factor_up_em(reaction%naqcomp), factor_dn_em(reaction%naqcomp)
  PetscReal :: conc_upd(reaction%naqcomp, sec_transport_vars%ncells)
  PetscReal :: total_upd(reaction%naqcomp, sec_transport_vars%ncells)
  PetscReal :: total_prev(reaction%naqcomp, sec_transport_vars%ncells)
  PetscReal :: mc_pri_molal(reaction%naqcomp, sec_transport_vars%ncells)
  PetscReal :: mc_sec_molal(reaction%neqcplx, sec_transport_vars%ncells)
  PetscReal :: total_primary_node(reaction%naqcomp)
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)
  PetscReal :: res_react(reaction%naqcomp)
  PetscReal :: jac_react(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: i, j, n
  PetscInt :: ngcells, ncomp
  PetscInt :: ni
  PetscInt :: nicomp, icomp
  PetscReal :: area_fm
  PetscReal :: diffusion_coefficient
  PetscReal :: porosity
  PetscReal :: pordt
  PetscReal :: pordiff_prim, pordiff_sec
  PetscReal :: tortuosity
  PetscReal :: portort
  PetscInt :: jcomp, lcomp, kcomp, icplx
  PetscReal :: inf_norm_sec
  class(material_auxvar_type), allocatable :: material_auxvar

  PetscReal :: total_sorb_upd(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: total_sorb_prev(reaction%naqcomp,sec_transport_vars%ncells)

  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  vol = sec_transport_vars%vol
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus
  area_fm = sec_transport_vars%interfacial_area
  ncomp = reaction%naqcomp

  do j = 1, ncomp
    do i = 1, ngcells
      total_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total(j,1)
      if (reaction%neqsorb > 0) then
        total_sorb_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total_sorb_eq(j)
      endif
    enddo
  enddo
  conc_upd = sec_transport_vars%updated_conc

  ! Note that sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j) units are in mol/kg
  ! Need to convert to mol/L since the units of total. in the Thomas
  ! algorithm are in mol/L
  res = 0.d0

  total_primary_node = auxvar%total(:,1)                         ! in mol/L
  pordt = porosity/option%tran_dt
  !pordiff = porosity*diffusion_coefficient
  !pordiff = porosity*diffusion_coefficient*tortuosity

  call RTAuxVarInit(rt_auxvar,reaction,option)
  do i = 1, ngcells
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i),option)
    rt_auxvar%pri_molal = conc_upd(:,i)
    ! sec_transport_vars%sec_rt_auxvar(i)%pri_molal =  conc_upd(:,i)
    mc_pri_molal(:, i) =  conc_upd(:,i)
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    if (reaction%neqsorb > 0) then
      ! call SecondaryRTotalSorb(rt_auxvar,global_auxvar,material_auxvar,reaction,option)
    endif
    total_upd(:,i) = rt_auxvar%total(:,1)
    ! sec_transport_vars%sec_rt_auxvar(i)%sec_molal = rt_auxvar%sec_molal
    mc_sec_molal(:, i) =  rt_auxvar%sec_molal
    if (reaction%neqsorb > 0) then
      total_sorb_upd(:,i) = rt_auxvar%total_sorb_eq(:)
    endif
  enddo

!================ Calculate the secondary residual =============================

  do i = 1, ngcells

    ! D_j \nabla C_j
    do icomp = 1, ncomp

      pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(icomp)*&
                     tortuosity*global_auxvar%den_kg(1)/1000.d0
      n = icomp + (i-1)*ncomp

      ! Accumulation
      res(n) = pordt*(total_upd(icomp,i) - total_prev(icomp,i))*vol(i)    ! in mol/L*m3/s
      if (reaction%neqsorb > 0) then
        res(n) = res(n) + vol(i)/option%tran_dt*(total_sorb_upd(icomp,i) - total_sorb_prev(icomp,i))
      endif

      ! Flux terms
      if (i.gt.1.and.i.lt.ngcells) then
        res(n) = res(n) - pordiff_prim*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                        (mc_pri_molal(icomp, i+1) - &
                        mc_pri_molal(icomp, i))
        res(n) = res(n) + pordiff_prim*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                        (mc_pri_molal(icomp, i) - &
                        mc_pri_molal(icomp, i-1))

      ! Apply boundary conditions
      ! Inner boundary
      else if (i.eq.1) then
        res(n) = res(n) - pordiff_prim*area(1)/(dm_minus(2) + dm_plus(1))* &
                        (mc_pri_molal(icomp, i+1) - &
                        mc_pri_molal(icomp, i))

      ! Outer boundary -- closest to primary node
      else !if (i.eq.ngcells) then
        res(n) = res(n) - &
                 pordiff_prim*area(ngcells)/dm_plus(ngcells)* &
                 (auxvar%pri_molal(icomp) - &
                 mc_pri_molal(icomp, i))
        res(n) = res(n) + &
                 pordiff_prim*area(ngcells-1)/(dm_minus(ngcells) &
                 + dm_plus(ngcells-1))* &
                 (mc_pri_molal(icomp, i) - &
                 mc_pri_molal(icomp, i-1))
      endif

    enddo


    ! summatory D_i \nabla C_i
    do icplx = 1, reaction%neqcplx

      pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)*&
                    tortuosity*global_auxvar%den_kg(1)/1000
      nicomp = reaction%eqcplxspecid(0,icplx)

      do ni = 1, nicomp
        icomp = reaction%eqcplxspecid(ni,icplx)

        n = icomp + (i-1)*ncomp

        ! Flux terms
        if (i.gt.1.and.i.lt.ngcells) then
          res(n) = res(n) - pordiff_sec*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                   (mc_sec_molal(icplx, i+1) - &
                   mc_sec_molal(icplx, i))* &
                   reaction%eqcplxstoich(ni,icplx)
          res(n) = res(n) + pordiff_sec*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                   (mc_sec_molal(icplx, i) - &
                   mc_sec_molal(icplx, i-1))* &
                   reaction%eqcplxstoich(ni,icplx)

        ! Apply boundary conditions
        ! Inner boundary
        else if (i.eq.1) then
          res(n) = res(n) - pordiff_sec*area(1)/(dm_minus(2) + dm_plus(1))* &
                   (mc_sec_molal(icplx, i+1) - &
                   mc_sec_molal(icplx, i))* &
                   reaction%eqcplxstoich(ni,icplx)

        ! Outer boundary -- closest to primary node
        else !if (i.eq.ngcells) then
          res(n) = res(n) - &
                   pordiff_sec*area(ngcells)/dm_plus(ngcells)* &
                   (auxvar%sec_molal(icplx) - &
                   mc_sec_molal(icplx, i))* &
                   reaction%eqcplxstoich(ni,icplx)
          res(n) = res(n) + pordiff_sec* &
                   area(ngcells-1)/(dm_minus(ngcells) + dm_plus(ngcells-1))* &
                   (mc_sec_molal(icplx, i) - &
                   mc_sec_molal(icplx, i-1))* &
                   reaction%eqcplxstoich(ni,icplx)
        endif

      enddo
    enddo


    ! Electromigration term
    call ComputeElectricPotentialTotalComponent_NP(i,reaction, &
                                            rt_parameter,  &
                                            rt_auxvar, &
                                            sec_transport_vars, &
                                            icomp,auxvar, &
                                            potent_dn, &
                                            potent_up, &
                                            mc_pri_molal, &
                                            mc_sec_molal)

    sum_denom_dn = 1d-40
    sum_denom_up = 1d-40

    do jcomp = 1, rt_parameter%naqcomp
      curCharge = reaction%primary_spec_Z(jcomp)
      sum_denom_up = sum_denom_up + potent_up(jcomp)*curCharge
      sum_denom_dn = sum_denom_dn + potent_dn(jcomp)*curCharge
    enddo

    do jcomp = 1, rt_parameter%naqcomp
      factor_up_em(jcomp) = potent_up(jcomp) / sum_denom_up
      factor_dn_em(jcomp) = potent_dn(jcomp) / sum_denom_dn
    enddo

    portort = porosity*tortuosity*global_auxvar%den_kg(1)/1000

    sum_transp_dn = 0d0
    sum_transp_up = 0d0

    ! Primary l index
    do lcomp = 1, rt_parameter%naqcomp
      curCharge = reaction%primary_spec_Z(lcomp)
      cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)

      ! Flux terms
      if (i.gt.1.and.i.lt.ngcells) then
        curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                     (mc_pri_molal(lcomp, i+1) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                     (mc_pri_molal(lcomp, i) - &
                     mc_pri_molal(lcomp, i-1))

      ! Apply boundary conditions
      ! Inner boundary
      else if (i.eq.1) then
        curConc_up = area(1)/(dm_minus(2) + dm_plus(1))* &
                     (mc_pri_molal(lcomp, i+1) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = 0.d0

      ! Outer boundary
      else !if (i.eq.ngcells) then
        curConc_up = area(ngcells)/dm_plus(ngcells)* &
                     (auxvar%pri_molal(lcomp) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                     + dm_plus(ngcells-1))* &
                     (mc_pri_molal(lcomp, i) - &
                     mc_pri_molal(lcomp, i-1))
      endif

      sum_transp_up = sum_transp_up - curCharge*cur_diff*(curConc_up)
      sum_transp_dn = sum_transp_dn + curCharge*cur_diff*(curConc_dn)

      ! Secondary term k
      do icplx = 1, reaction%neqcplx
        nicomp = reaction%eqcplxspecid(0,icplx)
        do ni = 1, nicomp
          kcomp = reaction%eqcplxspecid(ni,icplx)
          if (kcomp == lcomp) then
            curCharge = reaction%primary_spec_Z(kcomp)
            cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
            curStoich = reaction%eqcplxstoich(ni,icplx)

            ! Flux terms
            if (i.gt.1.and.i.lt.ngcells) then
              curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                          (mc_sec_molal(icplx, i+1) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                          (mc_sec_molal(icplx, i) - &
                          mc_sec_molal(icplx, i-1))

            ! Apply boundary conditions
            ! Inner boundary
            else if (i.eq.1) then
              curConc_up = area(1)/(dm_minus(2) + dm_plus(1))* &
                          (mc_sec_molal(icplx, i+1) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = 0.d0

            ! Outer boundary
            else !if (i.eq.ngcells) then
              curConc_up = area(ngcells)/dm_plus(ngcells)* &
                          (auxvar%sec_molal(icplx) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                          + dm_plus(ngcells-1))* &
                          (mc_sec_molal(icplx, i) - &
                          mc_sec_molal(icplx, i-1))
            endif

            sum_transp_up = sum_transp_up - curCharge*cur_diff*curStoich*(curConc_up)
            sum_transp_dn = sum_transp_dn + curCharge*cur_diff*curStoich*(curConc_dn)
          endif
        enddo
      enddo

    enddo

    do icomp = 1, rt_parameter%naqcomp
      n = icomp + (i-1)*ncomp
      res(n) = res(n) - portort*factor_dn_em(icomp)*sum_transp_dn
      res(n) = res(n) - portort*factor_up_em(icomp)*sum_transp_up
    enddo

  enddo

  res = res*1.d3 ! Convert mol/L*m3/s to mol/s

!====================== Add reaction contributions =============================

  ! Reaction
  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  do i = 1, ngcells
    res_react = 0.d0
    jac_react = 0.d0
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i), &
                      option)
    rt_auxvar%pri_molal = conc_upd(:,i) ! in mol/kg
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    material_auxvar%porosity = porosity
    material_auxvar%volume = vol(i)
    call RReaction(res_react,jac_react,PETSC_FALSE, &
                   rt_auxvar,global_auxvar,material_auxvar,reaction,option)
    do j = 1, ncomp
      res(j+(i-1)*ncomp) = res(j+(i-1)*ncomp) + res_react(j)
    enddo
  enddo
  call MaterialAuxVarStrip(material_auxvar)
  deallocate(material_auxvar)

 ! Need to decide how to scale the residual with volumes
  do i = 1, ngcells
    do j = 1, ncomp
      res(j+(i-1)*ncomp) = res(j+(i-1)*ncomp)/vol(i)
    enddo
  enddo

  inf_norm_sec = maxval(abs(res))
  call RTAuxVarStrip(rt_auxvar)

end subroutine SecondaryRTCheckResidual_np

! ************************************************************************** !

subroutine ComputeElectricPotentialTotalComponent_NP(ic,reaction, &
                                                    rt_parameter, &
                                                    rt_auxvar, &
                                                    sec_transport_vars, &
                                                    icomp,auxvar, &
                                                    potent_dn, &
                                                    potent_up, &
                                                    mc_pri_molal, &
                                                    mc_sec_molal)
!
! Computes the U electric potential vector term
!
! Author: Amphos21 - Barcelona Science
! Date: 29/06/20
!

  use Reaction_Aux_module
  use Reactive_Transport_Aux_module

  class(reaction_rt_type), pointer :: reaction
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(reactive_transport_auxvar_type) :: auxvar
  type(sec_transport_type) :: sec_transport_vars
  PetscReal :: mc_pri_molal(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: mc_sec_molal(reaction%neqcplx,sec_transport_vars%ncells)

  PetscInt :: icomp, jcomp
  PetscReal :: potent_up(reaction%naqcomp), potent_dn(reaction%naqcomp)
  PetscReal :: curCoef, curCharge, curStoich
  PetscReal :: curConc_up, curConc_dn
  PetscInt :: icplx, i, ncomp
  PetscInt :: ic
  PetscReal :: a,b,c
  PetscInt :: ngcells
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)
  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus

  ! Primary species contribution
  do icomp = 1, reaction%naqcomp
    curCoef = rt_parameter%pri_spec_diff_coef(icomp)
    curCharge = reaction%primary_spec_Z(icomp)

    ! Flux terms
    if (ic.gt.1.and.ic.lt.ngcells) then
      a = mc_pri_molal(icomp, ic-1)
      b = mc_pri_molal(icomp, ic)
      c = mc_pri_molal(icomp, ic+1)

      ! use of differential logarithmic avg
      if (log(b).eq.log(c)) then
        curConc_up = b
      else
        curConc_up = (c-b)/ (log(c) - log(b))
      endif

      if (log(a).eq.log(b)) then
        curConc_dn = b
      else
        curConc_dn = (b-a)/ (log(b) - log(a))
      endif

    ! Apply boundary conditions
    ! Inner boundary
    else if (ic.eq.1) then
      b = mc_pri_molal(icomp, ic)
      c = mc_pri_molal(icomp, ic+1)

      ! use of differential logarithmic avg
      if (log(b).eq.log(c)) then
        curConc_up = b
      else
        curConc_up = (c-b)/ (log(c) - log(b))
      endif

      curConc_dn = b

    ! Outer boundary -- closest to primary node
    else !if (i.eq.ngcells) then
      a = mc_pri_molal(icomp, ic-1)
      b = mc_pri_molal(icomp, ic)
      c = auxvar%pri_molal(icomp)

      ! use of differential logarithmic avg
      if (log(a).eq.log(b)) then
        curConc_dn = b
      else
        curConc_dn = (b-a)/ (log(b) - log(a))
      endif

      if (log(b).eq.log(c)) then
        curConc_up = b
      else
        curConc_up = (c-b)/ (log(c) - log(b))
      endif

    endif

    potent_up(icomp) = curCoef*curCharge*curConc_up
    potent_dn(icomp) = curCoef*curCharge*curConc_dn
  enddo


  ! Secondary species contribution
  do icplx = 1, reaction%neqcplx
    curCharge = reaction%eqcplx_Z(icplx)
    ncomp = reaction%eqcplxspecid(0,icplx)
    do i = 1, ncomp
      jcomp = reaction%eqcplxspecid(i,icplx)
      curCoef = rt_parameter%sec_spec_diff_coef(icplx)
      curStoich = reaction%eqcplxstoich(i,icplx)

      ! Flux terms
      if (ic.gt.1.and.ic.lt.ngcells) then
        a = mc_sec_molal(icplx, ic-1)
        b = mc_sec_molal(icplx, ic)
        c = mc_sec_molal(icplx, ic+1)

        ! use of differential logarithmic avg
        if (log(b).eq.log(c)) then
          curConc_up = b
        else
          curConc_up = (c-b)/ (log(c) - log(b))
        endif

        if (log(a).eq.log(b)) then
          curConc_dn = b
        else
          curConc_dn = (b-a)/ (log(b) - log(a))
        endif

      ! Apply boundary conditions
      ! Inner boundary
      else if (ic.eq.1) then
        b = mc_sec_molal(icplx, ic)
        c = mc_sec_molal(icplx, ic+1)

        ! use of differential logarithmic avg
        if (log(b).eq.log(c)) then
          curConc_up = b
        else
          curConc_up = (c-b)/ (log(c) - log(b))
        endif

        curConc_dn = b

      ! Outer boundary -- closest to primary node
      else !if (i.eq.ngcells) then
        a = mc_sec_molal(icplx, ic-1)
        b = mc_sec_molal(icplx, ic)
        c = auxvar%sec_molal(icplx)

        ! use of differential logarithmic avg
        if (log(b).eq.log(c)) then
          curConc_up = b
        else
          curConc_up = (c-b)/ (log(c) - log(b))
        endif

        if (log(a).eq.log(b)) then
          curConc_dn = b
        else
          curConc_dn = (b-a)/ (log(b) - log(a))
        endif

      endif

      potent_up(jcomp) = potent_up(jcomp) + curStoich * curCoef * curCharge * curConc_up
      potent_dn(jcomp) = potent_dn(jcomp) + curStoich * curCoef * curCharge * curConc_dn
    enddo
  enddo

end subroutine ComputeElectricPotentialTotalComponent_NP

end module Secondary_Continuum_NP_module
