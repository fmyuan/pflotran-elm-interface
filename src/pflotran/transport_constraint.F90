module Transport_Constraint_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  use Transport_Constraint_Base_module
  use Transport_Aux_module
  use Global_Aux_module

  implicit none

  private
  
   ! concentration subcondition types
  PetscInt, parameter, public :: CONSTRAINT_NULL = 0
  PetscInt, parameter, public :: CONSTRAINT_AQUEOUS = 1
  PetscInt, parameter, public :: CONSTRAINT_GASEOUS = 2
  
  type, public :: tran_constraint_ptr_type
    class(tran_constraint_base_type), pointer :: ptr
  end type tran_constraint_ptr_type
  
  type, public :: tran_constraint_list_type
    PetscInt :: num_constraints
    class(tran_constraint_base_type), pointer :: first
    class(tran_constraint_base_type), pointer :: last
    class(tran_constraint_ptr_type), pointer :: array(:)    
  end type tran_constraint_list_type

  type, public, extends(tran_constraint_base_type) :: tran_constraint_type
    type(aq_species_constraint_type), pointer :: aqueous_species
    type(gas_constraint_type), pointer :: gaseous_species
  contains
    procedure, public :: Strip => TranConstraintStrip
  end type tran_constraint_type

  type, public, extends(tran_constraint_coupler_base_type) :: &
               tran_constraint_coupler_type
    type(transport_auxvar_type), pointer :: tran_auxvar
  end type tran_constraint_coupler_type

  
  public :: TranConstraintAddToList, &
            TranConstraintInitList, &
            TranConstraintListDestroy, &
            TranConstraintGetPtrFromList, &
            TranConstraintDestroy, &
            TranConstraintCouplerDestroy

  public :: TranConstraintCreate, &
            TranConstraintCast, &
            TranConstraintCouplerCast, &
            TranConstraintCouplerCreate, &
            TranConstraintGetAuxVar, &
            TranConstraintRead
    
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

! ************************************************************************** !

function TranConstraintCreate(option)
  !
  ! Creates a transport constraint (set of concentrations
  ! and constraints for setting boundary or initial
  ! condition).
  !
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  !
  use Option_module

  implicit none

  type(option_type) :: option
  class(tran_constraint_type), pointer :: TranConstraintCreate

  class(tran_constraint_type), pointer :: constraint

  allocate(constraint)
  call TranConstraintBaseInit(constraint,option)
  nullify(constraint%aqueous_species)
  nullify(constraint%gaseous_species)

  TranConstraintCreate => constraint

end function TranConstraintCreate

! ************************************************************************** !

function TranConstraintCouplerCreate(option)
  !
  ! Creates a coupler that ties a constraint to a
  ! transport condition
  !
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  !

  use Option_module

  implicit none

  type(option_type) :: option
  class(tran_constraint_coupler_type), pointer :: TranConstraintCouplerCreate

  class(tran_constraint_coupler_type), pointer :: coupler

  allocate(coupler)
  call TranConstraintCouplerBaseInit(coupler,option)
  nullify(coupler%rt_auxvar)

  TranConstraintCouplerCreate => coupler

end function TranConstraintCouplerCreate

! ************************************************************************** !

function TranConstraintCast(this)
  !
  ! Casts a tran_constraint_base_type to tran_constraint_type
  !
  ! Author: Glenn Hammond
  ! Date: 10/07/19
  !

  implicit none

  class(tran_constraint_base_type), pointer :: this

  class(tran_constraint_type), pointer :: TranConstraintCast

  nullify(TranConstraintCast)
  select type (this)
    class is (tran_constraint_type)
      TranConstraintCast => this
  end select

end function TranConstraintCast

! ************************************************************************** !

function TranConstraintCouplerCast(this)
  !
  ! Casts a tran_constraint_coupler_base_type to
  ! tran_constraint_coupler_type
  !
  ! Author: Glenn Hammond
  ! Date: 10/07/19
  !

  implicit none

  class(tran_constraint_coupler_base_type), pointer :: this

  class(tran_constraint_coupler_type), pointer :: TranConstraintCouplerCast

  nullify(TranConstraintCouplerCast)
  select type (this)
    class is (tran_constraint_coupler_type)
      TranConstraintCouplerCast => this
  end select

end function TranConstraintCouplerCast

! ************************************************************************** !

function TranConstraintGetAuxVar(this)
  !
  ! Returns the auxvar associated with tran_constraint_coupler_type
  !
  ! Author: Glenn Hammond
  ! Date: 10/07/19
  !

  implicit none

  class(tran_constraint_coupler_base_type), pointer :: this

  type(reactive_transport_auxvar_type), pointer :: TranConstraintGetAuxVar

  nullify(TranConstraintGetAuxVar)
  select type (coupler=>this)
    class is (tran_constraint_coupler_type)
      TranConstraintGetAuxVar => coupler%rt_auxvar
  end select

end function TranConstraintGetAuxVar

! ************************************************************************** !

subroutine TranConstraintRead(constraint,transport,input,option)
  !
  ! Reads a transport constraint from the input file
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  !
  use Option_module
  use Input_Aux_module
  use Units_module
  use String_module
  use Logging_module

  implicit none

  class(tran_constraint_type) :: constraint
  class(transport_type) :: transport
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: block_string
  PetscInt :: icomp, igas
  PetscInt :: length
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(gas_species_constraint_type), pointer :: gas_species_constraint
  PetscErrorCode :: ierr
  PetscReal :: tempreal
  PetscBool :: found

  call PetscLogEventBegin(logging%event_tran_constraint_read, &
                          ierr);CHKERRQ(ierr)

  ! read the constraint
  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONSTRAINT')

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CONSTRAINT')

    call TranConstraintBaseRdSelectCase(constraint,input,word,found,option)
    if (found) cycle

    select case(trim(word))

      case('CONC','CONCENTRATIONS','Aqu', 'Aqueous')

        aq_species_constraint => &
          AqueousSpeciesConstraintCreate(transport,option)

        block_string = 'AQUEOUS, CONSTRAINT, CONCENTRATIONS'
        icomp = 0
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)

          if (InputCheckExit(input,option)) exit

          icomp = icomp + 1

          if (icomp > transport%naqcomp) then
            option%io_buffer = 'Number of concentration constraints ' // &
                               'exceeds number of primary chemical ' // &
                               'components in constraint: ' // &
                                trim(constraint%name)
            call PrintErrMsg(option)
          endif

          call InputReadCard(input,option,aq_species_constraint%names(icomp))
          call InputErrorMsg(input,option,'aqueous species name',block_string)
          option%io_buffer = 'Constraint Species: ' // &
                             trim(aq_species_constraint%names(icomp))
          call PrintMsg(option)

          call InputReadDouble(input,option, &
                               aq_species_constraint%constraint_conc(icomp))
          call InputErrorMsg(input,option,'concentration',block_string)

          aq_species_constraint%constraint_type(icomp) = CONSTRAINT_AQUEOUS

          call InputReadCard(input,option,word,PETSC_FALSE)
          if (input%ierr == 0) then
            call StringToUpper(word)
            select case(word)
              case('DATASET')
                call InputReadWord(input,option,aq_species_constraint% &
                                 constraint_aux_string(icomp),PETSC_TRUE)
                call InputErrorMsg(input,option,'dataset name', &
                                       block_string)
                aq_species_constraint%external_dataset(icomp) = PETSC_TRUE
            end select
          endif

        enddo
        call InputPopBlock(input,option)

        if (icomp < reaction%naqcomp) then
          option%io_buffer = &
                   'Number of concentration constraints is less than ' // &
                   'number of primary species in aqueous constraint.'
          call PrintErrMsg(option)
        endif
        if (icomp > reaction%naqcomp) then
          option%io_buffer = &
                   'Number of concentration constraints is greater than ' // &
                   'number of primary species in aqueous constraint.'
          call PrintWrnMsg(option)
        endif

        if (associated(constraint%aqueous_species)) &
          call AqueousSpeciesConstraintDestroy(constraint%aqueous_species)
        constraint%aqueous_species => aq_species_constraint

      case('GAS','GASEOUS')

        gas_species_constraint => &
          GaseousSpeciesConstraintCreate(transport,option)

        block_string = 'GASEOUS, CONSTRAINT, CONCENTRATIONS'
        icomp = 0
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)

          if (InputCheckExit(input,option)) exit

          icomp = icomp + 1

          if (icomp > transport%ngascomp) then
            option%io_buffer = 'Number of concentration constraints ' // &
                               'exceeds number of primary chemical ' // &
                               'components in constraint: ' // &
                                trim(constraint%name)
            call PrintErrMsg(option)
          endif

          call InputReadCard(input,option,gas_species_constraint%names(icomp))
          call InputErrorMsg(input,option,'gaseous species name',block_string)
          option%io_buffer = 'Constraint Species: ' // &
                             trim(gas_species_constraint%names(icomp))
          call PrintMsg(option)

          call InputReadDouble(input,option, &
                               gas_species_constraint%constraint_conc(icomp))
          call InputErrorMsg(input,option,'concentration',block_string)

          gas_species_constraint%constraint_type(icomp) = CONSTRAINT_GASEOUS

          call InputReadCard(input,option,word,PETSC_FALSE)
          if (input%ierr == 0) then
            call StringToUpper(word)
            select case(word)
              case('DATASET')
                call InputReadWord(input,option,gas_species_constraint% &
                                 constraint_aux_string(icomp),PETSC_TRUE)
                call InputErrorMsg(input,option,'dataset name', &
                                       block_string)
                gas_species_constraint%external_dataset(icomp) = PETSC_TRUE
            end select
          endif

        enddo
        call InputPopBlock(input,option)

        if (icomp < reaction%ngascomp) then
          option%io_buffer = &
                   'Number of concentration constraints is less than ' // &
                   'number of primary species in aqueous constraint.'
          call PrintErrMsg(option)
        endif
        if (icomp > reaction%ngascomp) then
          option%io_buffer = &
                   'Number of concentration constraints is greater than ' // &
                   'number of gas species in gasueous constraint.'
          call PrintWrnMsg(option)
        endif

        if (associated(constraint%gaseous_species)) &
          call AqueousSpeciesConstraintDestroy(constraint%gaseous_species)
        constraint%gaseous_species => gas_species_constraint

      case default
        call InputKeywordUnrecognized(input,word,'CONSTRAINT',option)
    end select

  enddo
  call InputPopBlock(input,option)

  call PetscLogEventEnd(logging%event_tran_constraint_read,ierr);CHKERRQ(ierr)

end subroutine TranConstraintRead

! ************************************************************************** !

subroutine TranConstraintStrip(this)
  !
  ! Deallocates a constraint
  !
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  !

  implicit none

  class(tran_constraint_type) :: this

  call TranConstraintBaseStrip(this)

  if (associated(this%aqueous_species)) &
    call AqueousSpeciesConstraintDestroy(this%aqueous_species)
  nullify(this%aqueous_species)

  if (associated(this%gaseous_species)) &
    call GaseousSpeciesConstraintDestroy(this%gaseous_species)
  nullify(this%aqueous_species)

end subroutine TranConstraintStrip

! ************************************************************************** !

subroutine TranConstraintCouplerStrip(this)
  !
  ! Deallocate dynamic members of class
  !
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  !
  implicit none

  class(tran_constraint_coupler_type) :: this

  call TranConstraintCouplerBaseStrip(this)

  call AuxVarDestroy(this%rt_auxvar)

end subroutine TranConstraintCouplerStrip

end module Transport_Constraint_module


