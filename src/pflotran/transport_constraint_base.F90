module Transport_Constraint_Base_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Global_Aux_module

  implicit none

  private
  
  type, public :: tran_constraint_base_type
    ! this class must be extended
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name         
    PetscBool :: equilibrate_at_each_cell
    class(tran_constraint_base_type), pointer :: next    
  contains
    procedure, public :: Strip => TranConstraintBaseStrip
  end type tran_constraint_base_type
  
  type, public :: tran_constraint_coupler_base_type
    ! this class must be extended
    character(len=MAXWORDLENGTH) :: constraint_name         
    PetscReal :: time
    PetscInt :: num_iterations
    PetscBool :: equilibrate_at_each_cell
    character(len=MAXWORDLENGTH) :: time_units
    class(tran_constraint_base_type), pointer :: constraint
    type(global_auxvar_type), pointer :: global_auxvar
    ! do not add PM auxvar here. they go in the daughter classes
    class(tran_constraint_coupler_base_type), pointer :: next   
  contains
    procedure, public :: Strip => TranConstraintCouplerBaseStrip
  end type tran_constraint_coupler_base_type
      
  public :: TranConstraintBaseInit, &
            TranConstraintCouplerBaseInit, &
            TranConstraintBaseRdSelectCase, &
            TranConstraintCouplerBaseStrip, &
            TranConstraintBaseStrip
    
contains

! ************************************************************************** !

subroutine TranConstraintBaseInit(this,option)
  ! 
  ! Initializes the base constraint class
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  ! 
  use Option_module
  
  implicit none
  
  class(tran_constraint_base_type) :: this
  type(option_type) :: option
  
  this%id = 0
  this%name = ''
  this%equilibrate_at_each_cell = PETSC_FALSE
  nullify(this%next)
  
end subroutine TranConstraintBaseInit

! ************************************************************************** !

subroutine TranConstraintCouplerBaseInit(this,option)
  ! 
  ! Initializes the base constraint coupler that ties a constraint to a
  ! transport condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  ! 

  use Option_module
  
  implicit none
  
  class(tran_constraint_coupler_base_type) :: this
  type(option_type) :: option
  
  this%constraint_name = ''
  this%time = 0.d0
  this%time_units = ''
  this%num_iterations = 0
  this%equilibrate_at_each_cell = PETSC_FALSE
  nullify(this%constraint)
  nullify(this%global_auxvar)
  nullify(this%next)
  
end subroutine TranConstraintCouplerBaseInit

! ************************************************************************** !

subroutine TranConstraintBaseRdSelectCase(this,input,keyword,found,option)
  ! 
  ! Reads a transport constraint from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  ! 
  use Option_module
  use Input_Aux_module
  use Units_module
  use String_module
  use Logging_module

  implicit none
  
  class(tran_constraint_base_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  type(option_type) :: option
  
  found = PETSC_TRUE

  select case(trim(keyword))
    case('EQUILIBRATE_AT_EACH_CELL')
      this%equilibrate_at_each_cell = PETSC_TRUE
    case('DO_NOT_EQUILIBRATE_AT_EACH_CELL')
      this%equilibrate_at_each_cell = PETSC_FALSE
    case default
      found = PETSC_FALSE
  end select 
  
end subroutine TranConstraintBaseRdSelectCase

! ************************************************************************** !

subroutine TranConstraintBaseStrip(this)
  ! 
  ! Strips any dynamically allocated members of base class
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  ! 

  implicit none
  
  class(tran_constraint_base_type) :: this
  
end subroutine TranConstraintBaseStrip

! ************************************************************************** !

subroutine TranConstraintCouplerBaseStrip(this)
  ! 
  ! Strips any dynamically allocated members of base class
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  ! 

  implicit none
  
  class(tran_constraint_coupler_base_type) :: this
  
  call GlobalAuxVarDestroy(this%global_auxvar)
  nullify(this%constraint)

end subroutine TranConstraintCouplerBaseStrip

end module Transport_Constraint_Base_module
