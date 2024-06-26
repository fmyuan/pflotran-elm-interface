module Reaction_Sandbox_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use Reaction_Sandbox_Calcite_class
  use Reaction_Sandbox_CLM_CN_class
  use Reaction_Sandbox_UFD_WP_class
  use Reaction_Sandbox_Example_class
  use Reaction_Sand_Equilibrate_class
  use Reaction_Sandbox_Simple_class
  use Reaction_Sandbox_Chromium_class
  use Reaction_Sandbox_Cyber_class
  use Reaction_Sandbox_Lambda_class
  use Reaction_Sandbox_Gas_class
  use Reaction_Sandbox_BioHill_class
  use Reaction_Sand_FlexBioHill_class
  use Reaction_Sandbox_BioTH_class
  use Reaction_Sandbox_Radon_class
  
  use Reaction_Sandbox_SomDec_class
  use Reaction_Sandbox_PlantN_class
  use Reaction_Sandbox_Langmuir_class
  use Reaction_Sandbox_Nitrif_class
  use Reaction_Sandbox_Denitr_class
  use Reaction_Sandbox_CNdegas_class

  ! Add new reacton sandbox classes here.

  use PFLOTRAN_Constants_module

  implicit none

  private

  class(reaction_sandbox_base_type), pointer, public :: rxn_sandbox_list

  interface RSandboxRead
    module procedure RSandboxRead1
    module procedure RSandboxRead2
  end interface

  interface RSandboxDestroy
    module procedure RSandboxDestroy1
    module procedure RSandboxDestroy2
  end interface

  public :: RSandboxInit, &
            RSandboxRead, &
            RSandboxSkipInput, &
            RSandboxSetup, &
            RSandboxEvaluate, &
            RSandboxUpdateKineticState, &
            RSandboxAuxiliaryPlotVariables, &
            RSandboxDestroy

contains

! ************************************************************************** !

subroutine RSandboxInit(option)
  !
  ! Initializes the sandbox list
  !
  ! Author: Glenn Hammond
  ! Date: 01/28/13
  !
  use Option_module
  implicit none
  type(option_type) :: option

  if (associated(rxn_sandbox_list)) then
    call RSandboxDestroy()
  endif
  nullify(rxn_sandbox_list)

end subroutine RSandboxInit

! ************************************************************************** !

subroutine RSandboxSetup(reaction,option)
  !
  ! Calls all the initialization routines for all reactions in
  ! the sandbox list
  !
  ! Author: Glenn Hammond
  ! Date: 01/28/13
  !

  use Option_module
  use Reaction_Aux_module, only : reaction_rt_type

  implicit none

  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  class(reaction_sandbox_base_type), pointer :: cur_sandbox

  ! sandbox reactions
  cur_sandbox => rxn_sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    call cur_sandbox%Setup(reaction,option)
    cur_sandbox => cur_sandbox%next
  enddo

end subroutine RSandboxSetup

! ************************************************************************** !

subroutine RSandboxRead1(input,option)
  !
  ! Reads input deck for reaction sandbox parameters
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/13
  !

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  call RSandboxRead(rxn_sandbox_list,input,option)

end subroutine RSandboxRead1

! ************************************************************************** !

subroutine RSandboxRead2(local_sandbox_list,input,option)
  !
  ! RSandboxRead: Reads input deck for reaction sandbox parameters
  !
  ! Author: Glenn Hammond
  ! Date: 11/08/12
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  class(reaction_sandbox_base_type), pointer :: local_sandbox_list
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  class(reaction_sandbox_base_type), pointer :: new_sandbox, cur_sandbox

  nullify(new_sandbox)
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,REACTION_SANDBOX')
    call StringToUpper(word)

    select case(trim(word))
      ! Add new cases statements for new reacton sandbox classes here.
      case('BIODEGRADATION_HILL')
        new_sandbox => BioHillCreate()
      case('BIOPARTICLE')
        new_sandbox => BioTH_Create()
      case('CALCITE')
        new_sandbox => CalciteCreate()
      case('CHROMIUM_REDUCTION')
        new_sandbox => ChromiumCreate()
      case('CLM-CN')
        new_sandbox => CLM_CN_Create()
      case('CYBERNETIC')
        new_sandbox => CyberCreate()
      case('EXAMPLE')
        new_sandbox => EXAMPLECreate()
      case('EQUILIBRATE')
        new_sandbox => EquilibrateCreate()
      case('FLEXIBLE_BIODEGRADATION_HILL')
        new_sandbox => FlexBioHillCreate()
      case('GAS')
        new_sandbox => GasCreate()
      case('LAMBDA')
        new_sandbox => LambdaCreate()
      case('RADON')
        new_sandbox => RadonCreate()
      case('SIMPLE')
        new_sandbox => SimpleCreate()
      case('UFD-WP')
        new_sandbox => WastePackageCreate()
      case('SOMDECOMP')
        new_sandbox => SomDecCreate()
      case('PLANTN')
        new_sandbox => PlantNCreate()
      case('NITRIFICATION')
        new_sandbox => NitrifCreate()
      case('DENITRIFICATION')
        new_sandbox => DenitrCreate()
      case('CNDEGAS')
        new_sandbox => CNdegasCreate()
      case('LANGMUIR')
        new_sandbox => LangmuirCreate()
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'CHEMISTRY,REACTION_SANDBOX',option)
    end select

    call new_sandbox%ReadInput(input,option)

    if (.not.associated(local_sandbox_list)) then
      local_sandbox_list => new_sandbox
    else
      cur_sandbox => local_sandbox_list
      do
        if (.not.associated(cur_sandbox%next)) exit
        cur_sandbox => cur_sandbox%next
      enddo
      cur_sandbox%next => new_sandbox
    endif
  enddo
  call InputPopBlock(input,option)

end subroutine RSandboxRead2

! ************************************************************************** !

subroutine RSandboxAuxiliaryPlotVariables(list,reaction,option)
  !
  ! Adds auxilairy plot variables to the list to be printed
  !
  ! Author: Glenn Hammond
  ! Date: 10/21/16
  !

  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module

  implicit none

  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction

  class(reaction_sandbox_base_type), pointer :: cur_reaction

  cur_reaction => rxn_sandbox_list
  do
    if (.not.associated(cur_reaction)) exit
    call cur_reaction%AuxiliaryPlotVariables(list,reaction,option)
    cur_reaction => cur_reaction%next
  enddo

end subroutine RSandboxAuxiliaryPlotVariables

! ************************************************************************** !

subroutine RSandboxSkipInput(input,option)
  !
  ! Intelligently skips over REACTION_SANDBOX block
  !
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  !

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  class(reaction_sandbox_base_type), pointer :: dummy_list

  nullify(dummy_list)
  call RSandboxRead(dummy_list,input,option)
  call RSandboxDestroy(dummy_list)

end subroutine RSandboxSkipInput

! ************************************************************************** !

subroutine RSandboxEvaluate(Residual,Jacobian,compute_derivative,rt_auxvar, &
                            global_auxvar,material_auxvar,reaction,option)
  !
  ! Evaluates reaction storing residual and/or Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 11/08/12
  !

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module, only: material_auxvar_type

  implicit none

  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  class(reaction_sandbox_base_type), pointer :: cur_reaction

  cur_reaction => rxn_sandbox_list
  do
    if (.not.associated(cur_reaction)) exit
    call cur_reaction%Evaluate(Residual,Jacobian,compute_derivative, &
                               rt_auxvar,global_auxvar,material_auxvar, &
                               reaction,option)
    cur_reaction => cur_reaction%next
  enddo

end subroutine RSandboxEvaluate

! ************************************************************************** !

subroutine RSandboxUpdateKineticState(rt_auxvar,global_auxvar, &
                                      material_auxvar,reaction,option)
  !
  ! Updates volume fractions, etc. at the end of a time step.
  !
  ! Author: Glenn Hammond
  ! Date: 09/06/16
  !

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module, only: material_auxvar_type

  implicit none

  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  class(reaction_sandbox_base_type), pointer :: cur_reaction

  cur_reaction => rxn_sandbox_list
  do
    if (.not.associated(cur_reaction)) exit
    call cur_reaction%UpdateKineticState(rt_auxvar,global_auxvar, &
                                         material_auxvar,reaction,option)
    cur_reaction => cur_reaction%next
  enddo

end subroutine RSandboxUpdateKineticState

! ************************************************************************** !

subroutine RSandboxDestroy1()
  !
  ! Destroys master sandbox list
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/13
  !

  implicit none

  call RSandboxDestroy(rxn_sandbox_list)

end subroutine RSandboxDestroy1

! ************************************************************************** !

subroutine RSandboxDestroy2(local_sandbox_list)
  !
  ! Destroys arbitrary sandbox list
  !
  ! Author: Glenn Hammond
  ! Date: 11/08/12
  !

  implicit none

  class(reaction_sandbox_base_type), pointer :: local_sandbox_list

  class(reaction_sandbox_base_type), pointer :: cur_sandbox, prev_sandbox

  ! sandbox reactions
  cur_sandbox => local_sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    prev_sandbox => cur_sandbox%next
    call cur_sandbox%Destroy()
    deallocate(cur_sandbox)
    cur_sandbox => prev_sandbox
  enddo
  nullify(local_sandbox_list)

end subroutine RSandboxDestroy2

end module Reaction_Sandbox_module
