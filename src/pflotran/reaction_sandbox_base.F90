module Reaction_Sandbox_Base_class
  
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
  type, abstract, public :: reaction_sandbox_base_type
    class(reaction_sandbox_base_type), pointer :: next
  contains
    procedure, public :: ReadInput => BaseRead
    procedure, public :: Setup => BaseSetup
    procedure, public :: Evaluate => BaseReact
    procedure, public :: UpdateKineticState => BaseUpdateKineticState
    procedure, public :: AuxiliaryPlotVariables => BaseAuxiliaryPlotVariables
    procedure, public :: Destroy => BaseDestroy    
  end type reaction_sandbox_base_type
  
contains

! ************************************************************************** !

  subroutine BaseSetup(this,reaction,option)
    
    use Option_module
    use Reaction_Aux_module
  
    implicit none
  
    class(reaction_sandbox_base_type) :: this
    class(reaction_rt_type) :: reaction
    type(option_type) :: option
  
  end subroutine BaseSetup 

! ************************************************************************** !

  subroutine BaseRead(this,input,option)
    
    use Option_module
    use Input_Aux_module
  
    implicit none
  
    class(reaction_sandbox_base_type) :: this
    type(input_type), pointer :: input
    type(option_type) :: option
  
  end subroutine BaseRead

! ************************************************************************** !

  subroutine BaseAuxiliaryPlotVariables(this,list,reaction,option)
    
    use Option_module
    use Reaction_Aux_module
    use Output_Aux_module
  
    implicit none
  
    class(reaction_sandbox_base_type) :: this
    type(output_variable_list_type), pointer :: list
    type(option_type) :: option
    class(reaction_rt_type) :: reaction
  
  end subroutine BaseAuxiliaryPlotVariables

! ************************************************************************** !

  subroutine BaseReact(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                        global_auxvar,material_auxvar,reaction,option)
    use Option_module
    use Reaction_Aux_module
    use Reactive_Transport_Aux_module
    use Global_Aux_module
    use Material_Aux_class
  
    implicit none
  
    class(reaction_sandbox_base_type) :: this
    type(option_type) :: option
    class(reaction_rt_type) :: reaction
    PetscBool :: compute_derivative
    ! the following arrays must be declared after reaction
    PetscReal :: Residual(reaction%ncomp)
    PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
    type(reactive_transport_auxvar_type) :: rt_auxvar
    type(global_auxvar_type) :: global_auxvar
    class(material_auxvar_type) :: material_auxvar
      
  end subroutine BaseReact

! ************************************************************************** !

  subroutine BaseUpdateKineticState(this,rt_auxvar,global_auxvar, &
                                     material_auxvar,reaction,option)
    use Option_module
    use Reaction_Aux_module
    use Reactive_Transport_Aux_module
    use Global_Aux_module
    use Material_Aux_class
  
    implicit none
  
    class(reaction_sandbox_base_type) :: this
    type(option_type) :: option
    class(reaction_rt_type) :: reaction
    type(reactive_transport_auxvar_type) :: rt_auxvar
    type(global_auxvar_type) :: global_auxvar
    class(material_auxvar_type) :: material_auxvar
      
  end subroutine BaseUpdateKineticState
  
! ************************************************************************** !

  subroutine BaseDestroy(this)

    implicit none
  
    class(reaction_sandbox_base_type) :: this

  end subroutine BaseDestroy  

end module Reaction_Sandbox_Base_class
