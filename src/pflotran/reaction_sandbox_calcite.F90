module Reaction_Sandbox_Calcite_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_calcite_type
    PetscInt :: auxiliary_offset
    PetscInt :: h_ion_id
    PetscInt :: calcium_id
    PetscInt :: bicarbonate_id
    PetscInt :: mineral_id
    PetscReal :: rate_constant1
    PetscReal :: rate_constant2
  contains
    procedure, public :: ReadInput => CalciteReadInput
    procedure, public :: Setup => CalciteSetup
    procedure, public :: AuxiliaryPlotVariables => CalciteAuxiliaryPlotVariables
    procedure, public :: Evaluate => CalciteEvaluate
    procedure, public :: UpdateKineticState => CalciteUpdateKineticState
  end type reaction_sandbox_calcite_type

  public :: CalciteCreate, &
            CalciteSetup

contains

! ************************************************************************** !

function CalciteCreate()
  !
  ! Allocates calcite reaction object.
  !
  implicit none

  class(reaction_sandbox_calcite_type), pointer :: CalciteCreate

  allocate(CalciteCreate)
  CalciteCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  CalciteCreate%h_ion_id = UNINITIALIZED_INTEGER
  CalciteCreate%calcium_id = UNINITIALIZED_INTEGER
  CalciteCreate%bicarbonate_id = UNINITIALIZED_INTEGER
  CalciteCreate%mineral_id = UNINITIALIZED_INTEGER
  CalciteCreate%rate_constant1 = UNINITIALIZED_DOUBLE
  CalciteCreate%rate_constant2 = UNINITIALIZED_DOUBLE
  nullify(CalciteCreate%next)

end function CalciteCreate

! ************************************************************************** !

subroutine CalciteReadInput(this,input,option)
  !
  ! Reads calcite reaction parameters
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(reaction_sandbox_calcite_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'CHEMISTRY,REACTION_SANDBOX,CALCITE'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    select case(word)
      case('RATE_CONSTANT1')
        call InputReadDouble(input,option,this%rate_constant1)
        call InputErrorMsg(input,option,word,error_string)
      case('RATE_CONSTANT2')
        call InputReadDouble(input,option,this%rate_constant2)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (Uninitialized(this%rate_constant1) .or. &
      Uninitialized(this%rate_constant2)) then
    option%io_buffer = 'RATE_CONSTANT1 and RATE_CONSTANT2 must be set for &
      REACTION_SANDBOX,CALCITE.'
    call PrintErrMsg(option)
  endif

end subroutine CalciteReadInput

! ************************************************************************** !

subroutine CalciteSetup(this,reaction,option)
  !
  ! Sets up the calcite reaction with hardwired parameters
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_calcite_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  ! rt_auxvar%auxiliary_data(:) is allocated to reaction%nauxiliary
  ! the offset points this sandbox to the correct entry for storing the rate
  this%auxiliary_offset = reaction%nauxiliary
  reaction%nauxiliary = reaction%nauxiliary + 1

  ! Aqueous species
  word = 'H+'
  this%h_ion_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Ca++'
  this%calcium_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'HCO3-'
  this%bicarbonate_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Calcite'
  this%mineral_id = &
    GetMineralIDFromName(word,reaction%mineral,option)

end subroutine CalciteSetup

! ************************************************************************** !

subroutine CalciteAuxiliaryPlotVariables(this,list,reaction,option)
  !
  ! Adds calcite auxiliary plot variables to output list
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY

  class(reaction_sandbox_calcite_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units

  word = 'Calcite Sandbox Rate'
  units = 'mol/m^3-sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)

end subroutine CalciteAuxiliaryPlotVariables

! ************************************************************************** !

subroutine CalciteEvaluate(this,Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !
  ! Evaluates calcite reaction storing residual but no Jacobian
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use Reaction_Mineral_Aux_module

  implicit none

  class(reaction_sandbox_calcite_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp) ! [mole / sec]
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  type(mineral_type), pointer :: mineral
  PetscReal :: volume               ! [m^3 bulk volume]
  PetscReal :: molality_to_molarity ! [kg water / L water]

  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: rate
  PetscReal :: drate_dQK
  PetscReal :: lnQK
  PetscReal :: QK
  PetscReal :: affinity_factor
  PetscReal :: dQK_dmj, dQK_dCj
  PetscReal :: sign_

  PetscInt :: ncomp
  PetscInt :: i, j
  PetscInt :: icomp, jcomp
  PetscInt :: imnrl
  PetscInt :: iauxiliary
  PetscBool :: calculate_rate

  mineral => reaction%mineral
  iauxiliary = this%auxiliary_offset + 1

  volume = material_auxvar%volume        ! den_kg [kg fluid / m^3 fluid]
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3  ! kg water/L water

  ! Reaction path #1
  ! the code in this block is very similar to the default mineral
  ! precipitation-dissolution capability in RKineticMineral(). this block
  ! illustrates how one can leverage the stoichiometries and logKs from the
  ! database and solely override the rate expression
  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  imnrl = this%mineral_id
  if (dabs(rt_auxvar%mnrl_rate(imnrl)) > 1.d-40) then
    option%io_buffer = 'For REACTION_SANDBOX CALCITE to function correctly, &
      &the RATE_CONSTANT in the default MINERAL_KINETICS block must be set &
      &to zero.'
    call PrintErrMsg(option)
  endif
  ! compute ion activity product
  lnQK = -mineral%kinmnrl_logK(imnrl)*LOG_TO_LN
  ! activity of water
  if (mineral%kinmnrlh2oid(imnrl) > 0) then
    lnQK = lnQK + mineral%kinmnrlh2ostoich(imnrl)* &
                  rt_auxvar%ln_act_h2o
  endif
  ! activity of other species
  ncomp = mineral%kinmnrlspecid(0,imnrl)
  do i = 1, ncomp
    icomp = mineral%kinmnrlspecid(i,imnrl)
    lnQK = lnQK + mineral%kinmnrlstoich(i,imnrl)*ln_act(icomp)
  enddo
  QK = exp(lnQK)
  affinity_factor = 1.d0-QK
  sign_ = sign(1.d0,affinity_factor)
  rate = 0.d0
  ! only calculate rate if mineral is present or precipitation (sign<0)
  calculate_rate = rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0
  if (calculate_rate) then
    ! mol/sec/m^3 bulk
    rate = -rt_auxvar%mnrl_area(imnrl) * &
           sign_ * abs(affinity_factor) * this%rate_constant1
  endif
  rt_auxvar%auxiliary_data(iauxiliary) = rate
  ! mol/sec
  rate = rate * material_auxvar%volume
  ncomp = mineral%kinmnrlspecid(0,imnrl)
  do i = 1, ncomp
    icomp = mineral%kinmnrlspecid(i,imnrl)
    Residual(icomp) = Residual(icomp) + mineral%kinmnrlstoich(i,imnrl)*rate
  enddo

  if (compute_derivative .and. calculate_rate) then
    ! mol/sec   ! m^2 mnrl/m^3 bulk          ! mol/m^2 mnrl/sec
    drate_dQK = rt_auxvar%mnrl_area(imnrl) * this%rate_constant1 * &
                ! m^3 bulk
                material_auxvar%volume
    do j = 1, ncomp
      jcomp = mineral%kinmnrlspecid(j,imnrl)
      ! units = kg water/mol
      dQK_dmj = mineral%kinmnrlstoich(j,imnrl)*QK*exp(-ln_conc(jcomp))* &
                molality_to_molarity
      do i = 1, ncomp
        icomp = mineral%kinmnrlspecid(i,imnrl)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jacobian(icomp,jcomp) = Jacobian(icomp,jcomp) + &
                                mineral%kinmnrlstoich(i,imnrl)* &
                                drate_dQK*dQK_dmj
      enddo
    enddo

  endif

  ! Reaction path #2
  ! this code block illustrate how to calculate the rate without the database

  ! QK = {Ca++}{HCO3-}/(Keq*{H+})
  ! 1.8487 is the pKeq from the database
  lnQK = -1.8487d0*LOG_TO_LN - ln_act(this%h_ion_id) + &
         ln_act(this%calcium_id) + ln_act(this%bicarbonate_id)
  affinity_factor = 1.d0-exp(lnQK)
  sign_ = sign(1.d0,affinity_factor)
  rate = 0.d0
  ! only calculate rate if mineral is present or precipitation (sign<0)
  calculate_rate = rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0
  if (calculate_rate) then
    ! mol/sec/m^3 bulk
    rate = -rt_auxvar%mnrl_area(imnrl) * &
           sign_ * abs(affinity_factor) * this%rate_constant2
  endif
  rt_auxvar%auxiliary_data(iauxiliary) = &
    rt_auxvar%auxiliary_data(iauxiliary) + rate
  ! mol/sec
  rate = rate * material_auxvar%volume
  Residual(this%h_ion_id) = Residual(this%h_ion_id) - rate
  Residual(this%calcium_id) = Residual(this%calcium_id) + rate
  Residual(this%bicarbonate_id) = Residual(this%bicarbonate_id) + rate

  if (compute_derivative .and. calculate_rate) then
    ! derivative of rate wrt affinity factor (1-QK)
    ! mol/sec   ! m^2 mnrl/m^3 bulk          ! mol/m^2 mnrl/sec
    drate_dQK = rt_auxvar%mnrl_area(imnrl) * this%rate_constant2 * &
                ! m^3 bulk
                material_auxvar%volume
    ! derivative wrt H+
    jcomp = this%h_ion_id
    ! stoich(H+)*QK/conc(H+)
    ! QK is a function of activities, dividing by conc leaves the
    ! activity coefficient in the derivative
              ! -1.d0 is the H+ stoichiometry in denominator
    dQK_dmj = -1.d0*exp(lnQK-ln_conc(jcomp))*molality_to_molarity
                                    ! subtract due to -1. H+ stoichiometry
    Jacobian(this%h_ion_id,jcomp) = &
      Jacobian(this%h_ion_id,jcomp) - drate_dQK * dQK_dmj
                                    ! add due to +1. Ca++ stoichiometry
    Jacobian(this%calcium_id,jcomp) = &
      Jacobian(this%calcium_id,jcomp) + drate_dQK * dQK_dmj
                                    ! add due to +1. HCO3- stoichiometry
    Jacobian(this%bicarbonate_id,jcomp) = &
      Jacobian(this%bicarbonate_id,jcomp) + drate_dQK * dQK_dmj
    ! derivative wrt Ca++
    jcomp = this%calcium_id
              ! 1.d0 is the Ca++ stoichiometry in denominator
    dQK_dmj = 1.d0*exp(lnQK-ln_conc(jcomp))*molality_to_molarity
    Jacobian(this%h_ion_id,jcomp) = &
      Jacobian(this%h_ion_id,jcomp) - drate_dQK * dQK_dmj
    Jacobian(this%calcium_id,jcomp) = &
      Jacobian(this%calcium_id,jcomp) + drate_dQK * dQK_dmj
    Jacobian(this%bicarbonate_id,jcomp) = &
      Jacobian(this%bicarbonate_id,jcomp) + drate_dQK * dQK_dmj
    ! derivative wrt HCO3-
    jcomp = this%bicarbonate_id
              ! 1.d0 is the HCO3- stoichiometry in denominator
    dQK_dmj = 1.d0*exp(lnQK-ln_conc(jcomp))*molality_to_molarity
    Jacobian(this%h_ion_id,jcomp) = &
      Jacobian(this%h_ion_id,jcomp) - drate_dQK * dQK_dmj
    Jacobian(this%calcium_id,jcomp) = &
      Jacobian(this%calcium_id,jcomp) + drate_dQK * dQK_dmj
    Jacobian(this%bicarbonate_id,jcomp) = &
      Jacobian(this%bicarbonate_id,jcomp) + drate_dQK * dQK_dmj
  endif

end subroutine CalciteEvaluate

! ************************************************************************** !

subroutine CalciteUpdateKineticState(this,rt_auxvar,global_auxvar, &
                                     material_auxvar,reaction,option)
  !
  ! Updates mineral volume fraction at end converged timestep based on latest
  ! rate
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_calcite_type) :: this
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscInt :: imnrl
  PetscReal :: delta_volfrac

  imnrl = this%mineral_id
  ! rate = mol/m^3/sec
  ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) *
  !                                mol_vol (m^3 mnrl/mol mnrl)
  delta_volfrac = rt_auxvar%auxiliary_data(this%auxiliary_offset+1)* &
                  reaction%mineral%kinmnrl_molar_vol(imnrl)* &
                  option%tran_dt
  ! m^3 mnrl/m^3 bulk
  rt_auxvar%mnrl_volfrac(imnrl) = rt_auxvar%mnrl_volfrac(imnrl) + &
                                  delta_volfrac
  ! zero to avoid negative volume fractions
  if (rt_auxvar%mnrl_volfrac(imnrl) < 0.d0) &
    rt_auxvar%mnrl_volfrac(imnrl) = 0.d0

end subroutine CalciteUpdateKineticState

end module Reaction_Sandbox_Calcite_class
