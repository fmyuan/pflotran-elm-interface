module Transport_Constraint_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  use Transport_Constraint_Base_module
  use Transport_Constraint_RT_module
  use Transport_Constraint_NWT_module

  implicit none

  private
  
  
  type, public :: tran_constraint_ptr_type
    type(tran_constraint_base_type), pointer :: ptr
  end type tran_constraint_ptr_type
  
  type, public :: tran_constraint_list_type
    PetscInt :: num_constraints
    class(tran_constraint_base_type), pointer :: first
    class(tran_constraint_base_type), pointer :: last
    class(tran_constraint_ptr_type), pointer :: array(:)    
  end type tran_constraint_list_type
  
  public :: TranConstraintAddToList, &
            TranConstraintInitList, &
            TranConstraintListDestroy, &
            TranConstraintGetPtrFromList, &
            TranConstraintDestroy, &
            TranConstraintCouplerDestroy
    
contains

! ************************************************************************** !

subroutine TranConstraintInitList(list)
  ! 
  ! Initializes a transport constraint list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  implicit none

  type(tran_constraint_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_constraints = 0

end subroutine TranConstraintInitList

! ************************************************************************** !

subroutine TranConstraintAddToList(new_constraint,list)
  ! 
  ! Adds a new constraint to a transport constraint
  ! list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  implicit none

  class(tran_constraint_base_type), pointer :: new_constraint
  type(tran_constraint_list_type) :: list

  list%num_constraints = list%num_constraints + 1
  new_constraint%id = list%num_constraints
  if (.not.associated(list%first)) list%first => new_constraint
  if (associated(list%last)) list%last%next => new_constraint
  list%last => new_constraint

end subroutine TranConstraintAddToList

! ************************************************************************** !

function TranConstraintGetPtrFromList(constraint_name,constraint_list)
  ! 
  ! Returns a pointer to the constraint matching
  ! constraint_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/08
  ! 

  use String_module

  implicit none

  class(tran_constraint_base_type), pointer :: TranConstraintGetPtrFromList
  character(len=MAXWORDLENGTH) :: constraint_name
  type(tran_constraint_list_type) :: constraint_list

  PetscInt :: length
  class(tran_constraint_base_type), pointer :: cur_constraint

  nullify(TranConstraintGetPtrFromList)
  cur_constraint => constraint_list%first

  do
    if (.not.associated(cur_constraint)) exit
    length = len_trim(constraint_name)
    if (length == len_trim(cur_constraint%name) .and. &
        StringCompare(cur_constraint%name,constraint_name, &
                        length)) then
      TranConstraintGetPtrFromList => cur_constraint
      return
    endif
    cur_constraint => cur_constraint%next
  enddo

end function TranConstraintGetPtrFromList

! ************************************************************************** !

subroutine TranConstraintListDestroy(constraint_list)
  ! 
  ! Deallocates a list of constraints
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  ! 

  implicit none

  type(tran_constraint_list_type), pointer :: constraint_list

  class(tran_constraint_base_type), pointer :: constraint, prev_constraint

  if (.not.associated(constraint_list)) return

  if (associated(constraint_list%first)) then
    call TranConstraintDestroy(constraint_list%first)
    nullify(constraint_list%first)
  endif

  constraint_list%num_constraints = 0
  nullify(constraint_list%last)

  deallocate(constraint_list)
  nullify(constraint_list)

end subroutine TranConstraintListDestroy

! ************************************************************************** !

recursive subroutine TranConstraintDestroy(constraint)
  ! 
  ! Strips any dynamically allocated members of base class
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  ! 
  implicit none

  class(tran_constraint_base_type), pointer :: constraint

  if (.not.associated(constraint)) return

  if (associated(constraint%next)) then
    call TranConstraintDestroy(constraint%next)
  endif

  call constraint%Strip()
  deallocate(constraint)
  nullify(constraint)

end subroutine TranConstraintDestroy

! ************************************************************************** !

recursive subroutine TranConstraintCouplerDestroy(coupler)
  ! 
  ! Destroys a constraint coupler linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/09
  ! 

  use Option_module

  implicit none

  class(tran_constraint_coupler_base_type), pointer :: coupler

  if (.not.associated(coupler)) return

  if (associated(coupler%next)) then
    call TranConstraintCouplerDestroy(coupler%next)
  endif

  call coupler%Strip()

  deallocate(coupler)
  nullify(coupler)

end subroutine TranConstraintCouplerDestroy

end module Transport_Constraint_module
