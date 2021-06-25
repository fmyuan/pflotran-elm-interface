module Reaction_Redox_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Reaction_Gas_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: RRedoxCalcpH, &
            RRedoxCalcLogKEh, &
            RRedoxCalcEhpe, &
            RRedoxCalcLnFO2

contains

! ************************************************************************** !

subroutine RRedoxCalcpH(rt_auxvar,global_auxvar,reaction,ph,option)

  ! Calculates pH

  ! Author: Glenn Hammond
  !         Based on Peter Lichtner's previous implementation
  ! Date: 06/25/21

  use Option_module

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type) :: reaction
  PetscReal :: ph
  type(option_type) :: option

  ! pH
  if (reaction%species_idx%h_ion_id > 0) then
    ph = &
      -log10(rt_auxvar%pri_molal(reaction%species_idx%h_ion_id)* &
             rt_auxvar%pri_act_coef(reaction%species_idx%h_ion_id))
  else if (reaction%species_idx%h_ion_id < 0) then
    ph = &
      -log10(rt_auxvar%sec_molal(abs(reaction%species_idx%h_ion_id))* &
             rt_auxvar%sec_act_coef(abs(reaction%species_idx%h_ion_id)))
  endif

end subroutine RRedoxCalcpH

! ************************************************************************** !

subroutine RRedoxCalcEhpe(rt_auxvar,global_auxvar,reaction,eh,pe,option)

  ! Calculates pH, Eh and pe given an hydrogen ion and oxygen concentration

  ! Author: Glenn Hammond
  !         Based on Peter Lichtner's previous implementation
  ! Date: 06/25/21

  use Option_module

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type) :: reaction
  PetscReal :: eh
  PetscReal :: pe
  type(option_type) :: option

  PetscReal :: ph
  PetscReal :: lnQKo2
  PetscReal :: t_kelvin
  PetscReal :: ehfac

  ! pH
  call RRedoxCalcpH(rt_auxvar,global_auxvar,reaction,ph,option)
  call RRedoxCalcLnFO2(rt_auxvar,global_auxvar,reaction,lnQKo2,option)
  t_kelvin = global_auxvar%temp+273.15d0
  ehfac = IDEAL_GAS_CONSTANT*t_kelvin*LOG_TO_LN/FARADAY
  eh = ehfac*(-4.d0*ph+lnQKo2*LN_TO_LOG+ &
              RRedoxCalcLogKEh(t_kelvin))/4.d0
  pe = eh/ehfac

end subroutine RRedoxCalcEhpe

! ************************************************************************** !

subroutine RRedoxCalcLnFO2(rt_auxvar,global_auxvar,reaction,lnfo2,option)

  ! Calculates the natural log of O2 partial pressure

  ! Author: Glenn Hammond
  !         Based on Peter Lichtner's previous implementation
  ! Date: 06/25/21

  use Option_module

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type) :: reaction
  PetscReal :: lnfo2
  type(option_type) :: option

  PetscInt :: io2gas
  PetscReal :: lnQKo2
  PetscInt :: jcomp
  PetscInt :: comp_id

  io2gas = reaction%species_idx%o2_gas_id
    ! compute gas partial pressure
  lnQKo2 = -reaction%gas%paseqlogK(io2gas)*LOG_TO_LN
  ! activity of water
  if (reaction%gas%paseqh2oid(io2gas) > 0) then
    lnQKo2 = lnQKo2 + reaction%gas%paseqh2ostoich(io2gas) * &
                      rt_auxvar%ln_act_h2o
  endif
  do jcomp = 1, reaction%gas%paseqspecid(0,io2gas)
    comp_id = reaction%gas%paseqspecid(jcomp,io2gas)
    lnQKo2 = lnQKo2 + reaction%gas%paseqstoich(jcomp,io2gas)* &
                      log(rt_auxvar%pri_molal(comp_id)* &
                          rt_auxvar%pri_act_coef(comp_id))
  enddo
  lnfo2 = lnQKo2

end subroutine RRedoxCalcLnFO2

! ************************************************************************** !

function RRedoxCalcLogKEh(tk)
  !
  ! Function logKeh: Maier-Kelly fit to equilibrium constant half-cell
  ! reaction 2 H2O - 4 H+ - 4 e- = O2, to compute Eh and pe.
  !
  ! Author: Peter Lichtner
  ! Date: 04/27/13
  !
  PetscReal, intent(in) :: tk

  PetscReal :: RRedoxCalcLogKEh

  PetscReal, parameter :: cm1 = 6.745529048112373d0
  PetscReal, parameter :: c0 = -48.295936593543715d0
  PetscReal, parameter :: c1 = 0.0005578156078778505d0
  PetscReal, parameter :: c2 = 27780.749538022003d0
  PetscReal, parameter :: c3 = 4027.3376948579394d0

  RRedoxCalcLogKEh = cm1 * log(tk) + c0 + c1 * tk + c2 / tk + c3 / (tk * tk)

end function RRedoxCalcLogKEh

end module Reaction_Redox_module
