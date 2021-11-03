module NWT_Equilibrium_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use NW_Transport_Aux_module

  implicit none

  public :: NWTEquilibrateConstraint, &
            NWTEqDissPrecipSorb

contains

! ************************************************************************** !

subroutine NWTEquilibrateConstraint(reaction_nw,constraint,nwt_auxvar, &
                                    global_auxvar,material_auxvar, &
                                    option)
  ! 
  ! Calculates the transport constraints based on equilibrium conditions.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/30/2019
  ! 
  
  use Option_module
  use Global_Aux_module
  use Material_Aux_class
  use Transport_Constraint_NWT_module
  
  implicit none
  
  class(reaction_nw_type), pointer :: reaction_nw
  class(tran_constraint_nwt_type) :: constraint
  type(nw_transport_auxvar_type) :: nwt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  
  type(species_type), pointer :: cur_species
  type(nwt_species_constraint_type), pointer :: nwt_species
  PetscInt :: ispecies
  PetscInt :: c_type
  PetscBool :: dry_out
  PetscReal :: solubility(reaction_nw%params%nspecies)  ! [mol/m^3-liq]
  PetscReal :: mnrl_molar_density(reaction_nw%params%nspecies)  ! [mol/m^3-mnrl]
  PetscReal :: ele_kd(reaction_nw%params%nspecies)  ! [m^3-water/m^3-bulk]
  PetscReal :: aq_mass     ! [mol/m^3-liq]
  PetscReal :: ppt_mass    ! [mol/m^3-bulk]
  PetscReal :: sorb_mass   ! [mol/m^3-bulk]
  PetscReal :: sat, por
  PetscReal :: extra_mass  ! [mol/m^3-liq]

  nwt_species => constraint%nwt_species
  
  sat = max(MIN_LIQ_SAT,global_auxvar%sat(LIQUID_PHASE))
  por = material_auxvar%porosity

  cur_species => reaction_nw%species_list
  do 
    if (.not.associated(cur_species)) exit
    solubility(cur_species%id) = cur_species%solubility_limit
    mnrl_molar_density(cur_species%id) = cur_species%mnrl_molar_density
    ele_kd(cur_species%id) = cur_species%ele_kd
    cur_species => cur_species%next
  enddo
  
  if (sat > 0.d0) then
    dry_out = PETSC_FALSE
  else
    dry_out = PETSC_TRUE
  endif
  
  do ispecies = 1,reaction_nw%params%nspecies
  
    c_type = nwt_species%constraint_type(ispecies)
    select case(c_type)
    !---------------------------------------
      case(CONSTRAINT_T_EQUILIBRIUM)
        nwt_auxvar%total_bulk_conc(ispecies) = &
                              nwt_species%constraint_conc(ispecies)
        ! check aqueous concentration against solubility limit and update
        call NWTEqDissPrecipSorb(solubility(ispecies),material_auxvar, &
                                 global_auxvar,dry_out,ele_kd(ispecies), &
                                 nwt_auxvar%total_bulk_conc(ispecies), &
                                 aq_mass,ppt_mass,sorb_mass)   
        nwt_auxvar%aqueous_eq_conc(ispecies) = aq_mass      
        nwt_auxvar%sorb_eq_conc(ispecies) = sorb_mass
        nwt_auxvar%mnrl_eq_conc(ispecies) = ppt_mass
        nwt_auxvar%mnrl_vol_frac(:) = nwt_auxvar%mnrl_eq_conc(:)/ &
                                      (por*mnrl_molar_density(:))
    !---------------------------------------
      case(CONSTRAINT_AQ_EQUILIBRIUM)
        nwt_auxvar%aqueous_eq_conc(ispecies) = &
                                         nwt_species%constraint_conc(ispecies)
        !if (dry_out) then
        !  extra_mass = nwt_auxvar%aqueous_eq_conc(ispecies)
        !  nwt_auxvar%aqueous_eq_conc(ispecies) = 0.0d0
        !  nwt_auxvar%sorb_eq_conc(ispecies) = 0.d0
        !else
          if (nwt_auxvar%aqueous_eq_conc(ispecies) > solubility(ispecies)) then
            extra_mass = nwt_auxvar%aqueous_eq_conc(ispecies) &
                         - solubility(ispecies)
            nwt_auxvar%aqueous_eq_conc(ispecies) = solubility(ispecies)
          else
            extra_mass = 0.d0
          endif
          nwt_auxvar%sorb_eq_conc(ispecies) = &
                          nwt_auxvar%aqueous_eq_conc(ispecies)*ele_kd(ispecies)
        !endif
        nwt_auxvar%mnrl_eq_conc(ispecies) = extra_mass*(sat*por)
        nwt_auxvar%total_bulk_conc(ispecies) = &
                            (nwt_auxvar%aqueous_eq_conc(ispecies)*sat*por) + &
                            nwt_auxvar%mnrl_eq_conc(ispecies) + &
                            nwt_auxvar%sorb_eq_conc(ispecies)
    !---------------------------------------
      case(CONSTRAINT_PPT_EQUILIBRIUM)
        nwt_auxvar%mnrl_eq_conc(ispecies) = &
                              nwt_species%constraint_conc(ispecies)
        nwt_auxvar%mnrl_vol_frac(ispecies) = &
                                      nwt_auxvar%mnrl_eq_conc(ispecies)/ &
                                      (por*mnrl_molar_density(ispecies))
        !if (dry_out) then
        !  nwt_auxvar%aqueous_eq_conc(ispecies) = 0.0d0
        !  nwt_auxvar%sorb_eq_conc(ispecies) = 0.0d0
        !else
          nwt_auxvar%aqueous_eq_conc(ispecies) = solubility(ispecies)
          nwt_auxvar%sorb_eq_conc(ispecies) = solubility(ispecies)* &
                                              ele_kd(ispecies)
        !endif
        nwt_auxvar%total_bulk_conc(ispecies) = &
                            (nwt_auxvar%aqueous_eq_conc(ispecies)*sat*por) + &
                            nwt_auxvar%mnrl_eq_conc(ispecies) + &
                            nwt_auxvar%sorb_eq_conc(ispecies)
    !---------------------------------------
      case(CONSTRAINT_MNRL_VOL_FRAC_EQ)
        nwt_auxvar%mnrl_vol_frac(ispecies) = &
                              nwt_species%constraint_conc(ispecies)
        nwt_auxvar%mnrl_eq_conc(ispecies) = &
                          nwt_auxvar%mnrl_vol_frac(ispecies)* &
                          material_auxvar%porosity*mnrl_molar_density(ispecies)
        !if (dry_out) then
        !  nwt_auxvar%aqueous_eq_conc(ispecies) = 0.0d0
        !  nwt_auxvar%sorb_eq_conc(ispecies) = 0.0d0
        !else
          nwt_auxvar%aqueous_eq_conc(ispecies) = solubility(ispecies)  
          nwt_auxvar%sorb_eq_conc(ispecies) = solubility(ispecies)* &
                                              ele_kd(ispecies)
        !endif
        nwt_auxvar%total_bulk_conc(ispecies) = &
                            (nwt_auxvar%aqueous_eq_conc(ispecies)*sat*por) + &
                            nwt_auxvar%mnrl_eq_conc(ispecies) + &
                            nwt_auxvar%sorb_eq_conc(ispecies)
    !---------------------------------------
      case(CONSTRAINT_SB_EQUILIBRIUM)
        if (ele_kd(ispecies) == 0.d0) then
          option%io_buffer = 'No value given for elemental Kd, but the &
            &concentration is being constrained by SB.'
          call PrintErrMsg(option)
        endif
        ! use the sorbed constraint to get to total_bulk_conc
        nwt_auxvar%sorb_eq_conc(ispecies) = &
                              nwt_species%constraint_conc(ispecies)
        nwt_auxvar%aqueous_eq_conc(ispecies) = &
                            nwt_auxvar%sorb_eq_conc(ispecies)/ele_kd(ispecies)
        if (nwt_auxvar%aqueous_eq_conc(ispecies) > solubility(ispecies)) then
          extra_mass = nwt_auxvar%aqueous_eq_conc(ispecies) &
                       - solubility(ispecies)
          nwt_auxvar%aqueous_eq_conc(ispecies) = solubility(ispecies)
        else
          extra_mass = 0.d0
        endif
        nwt_auxvar%mnrl_eq_conc(ispecies) = extra_mass*(sat*por)
        nwt_auxvar%total_bulk_conc(ispecies) = &
                            (nwt_auxvar%aqueous_eq_conc(ispecies)*sat*por) + &
                            nwt_auxvar%mnrl_eq_conc(ispecies) + &
                            nwt_auxvar%sorb_eq_conc(ispecies)
        ! next, re-partition from the total_bulk_conc to correct the sorbed
        ! mass, in case the sorbed constraint caused super-saturated conditions,
        ! and/or to correct for dry-out 
        call NWTEqDissPrecipSorb(solubility(ispecies),material_auxvar, &
                                 global_auxvar,dry_out,ele_kd(ispecies), &
                                 nwt_auxvar%total_bulk_conc(ispecies), &
                                 aq_mass,ppt_mass,sorb_mass) 
        nwt_auxvar%aqueous_eq_conc(ispecies) = aq_mass
        nwt_auxvar%sorb_eq_conc(ispecies) = sorb_mass
        nwt_auxvar%mnrl_eq_conc(ispecies) = ppt_mass
        nwt_auxvar%mnrl_vol_frac(ispecies) = &
                                  nwt_auxvar%mnrl_eq_conc(ispecies)/ &
                                 (por*mnrl_molar_density(ispecies))
    !---------------------------------------
    end select
  
  enddo

end subroutine NWTEquilibrateConstraint

! ************************************************************************** !

subroutine NWTEqDissPrecipSorb(solubility,material_auxvar,global_auxvar, &
                               dry_out,ele_kd,total_bulk_conc,aqueous_eq_conc, &
                               ppt_mass_conc,sorb_mass_conc)
  ! 
  ! Computes the equilibrium state and partitions the total_bulk_conc.
  ! 
  ! Author: Jenn Frederick
  ! Date: 07/15/2019
  ! 
 
  use Material_Aux_class
  use Global_Aux_module
  
  implicit none

  PetscReal :: solubility       ! [mol/m^3-liq]
  class(material_auxvar_type) :: material_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscBool :: dry_out
  PetscReal :: ele_kd           ! [m^3-water/m^3-bulk]
  PetscReal :: total_bulk_conc  ! [mol/m^3-bulk]
  PetscReal :: aqueous_eq_conc  ! [mol/m^3-liq]
  PetscReal :: ppt_mass_conc    ! [mol/m^3-bulk]
  PetscReal :: sorb_mass_conc   ! [mol/m^3-bulk]
  
  PetscReal :: extra_mass_conc  ! [mol/m^3-liq]
  PetscReal :: por, sat
  PetscBool :: super_saturated

  por = material_auxvar%porosity
  sat = max(MIN_LIQ_SAT,global_auxvar%sat(LIQUID_PHASE))
  super_saturated = PETSC_FALSE

  !if (.not.dry_out) then
  !---- Cell is wet ----!
    aqueous_eq_conc = total_bulk_conc/(por*sat)             ! [mol/m^3-liq]
    if (aqueous_eq_conc > solubility) then
      extra_mass_conc = aqueous_eq_conc - solubility        ! [mol/m^3-liq]
      aqueous_eq_conc = solubility                          ! [mol/m^3-liq]
      ppt_mass_conc = extra_mass_conc*sat*por               ! [mol/m^3-bulk]
      super_saturated = PETSC_TRUE
    else
      ppt_mass_conc = 0.d0                                  ! [mol/m^3-bulk]
      super_saturated = PETSC_FALSE
    endif
    sorb_mass_conc = aqueous_eq_conc*ele_kd                 ! [mol/m^3-bulk]
    total_bulk_conc = (aqueous_eq_conc*sat*por) + &         ! [mol/m^3-bulk]
                      sorb_mass_conc + ppt_mass_conc        ! [mol/m^3-bulk]
  !else
  !!---- Cell is dry ---!
  !  ppt_mass_conc = total_bulk_conc  ! [mol/m^3-bulk]
  !  aqueous_eq_conc = 0.d0           ! [mol/m^3-liq]
  !  sorb_mass_conc = 0.d0            ! [mol/m^3-bulk]
  !  aqueous_eq_conc = 0.d0           ! [mol/m^3-liq]
  !endif

end subroutine NWTEqDissPrecipSorb

! ************************************************************************** !

end module NWT_Equilibrium_module
