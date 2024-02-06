module Reaction_Equation_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: reaction_equation_type
    PetscInt :: nspec
    character(len=MAXWORDLENGTH), pointer :: spec_name(:)
    PetscReal, pointer :: stoich(:)
    PetscInt, pointer :: specid(:)
  end type reaction_equation_type

  type, public :: reaction_equation_ptr_type
    type(reaction_equation_type), pointer :: reaction_equation
    type(reaction_equation_ptr_type), pointer :: next
  end type reaction_equation_ptr_type

  interface ReactionEquationCreate
    module procedure ReactionEquationCreate1
    module procedure ReactionEquationCreate2
  end interface ReactionEquationCreate

  public :: ReactionEquationCreate, &
            ReactionEquationCreateFromString, &
            ReactionEquationSubSpecInRxn, &
            ReactionEquationCreateRxnPtr, &
            ReactionEquationDestroyRxnPtr, &
            ReactionEquationDestroy

contains

! ************************************************************************** !

function ReactionEquationCreate1()
  !
  ! Allocate and initialize a reaction equation
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  type(reaction_equation_type), pointer :: ReactionEquationCreate1

  type(reaction_equation_type), pointer :: reaction_equation

  allocate(reaction_equation)
  reaction_equation%nspec = UNINITIALIZED_INTEGER
  nullify(reaction_equation%spec_name)
  nullify(reaction_equation%stoich)
  nullify(reaction_equation%specid)

  ReactionEquationCreate1 => reaction_equation

end function ReactionEquationCreate1

! ************************************************************************** !

function ReactionEquationCreate2(num_species)
  !
  ! Allocate and initialize a reaction equation
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  PetscInt :: num_species

  type(reaction_equation_type), pointer :: ReactionEquationCreate2

  type(reaction_equation_type), pointer :: reaction_equation

  reaction_equation => ReactionEquationCreate()
  reaction_equation%nspec = num_species
  allocate(reaction_equation%spec_name(num_species))
  reaction_equation%spec_name = ''
  allocate(reaction_equation%stoich(num_species))
  reaction_equation%stoich = UNINITIALIZED_DOUBLE
  allocate(reaction_equation%specid(num_species))
  reaction_equation%specid = UNINITIALIZED_INTEGER

  ReactionEquationCreate2 => reaction_equation

end function ReactionEquationCreate2

! ************************************************************************** !

function ReactionEquationCreateRxnPtr()
  !
  ! Allocate and initialize a pointer to the reaction equation
  !
  ! Author: Glenn Hammond
  ! Date: 09/01/08
  !

  implicit none

  type(reaction_equation_ptr_type), pointer :: ReactionEquationCreateRxnPtr

  type(reaction_equation_ptr_type), pointer :: reaction_equation_ptr

  allocate(reaction_equation_ptr)
  nullify(reaction_equation_ptr%reaction_equation)
  nullify(reaction_equation_ptr%next)

  ReactionEquationCreateRxnPtr => reaction_equation_ptr

end function ReactionEquationCreateRxnPtr

! ************************************************************************** !

function ReactionEquationCreateFromString(reaction_string, &
                                          naqcomp, aq_offset, &
                                          primary_aq_species_names, &
                                          nimcomp, im_offset, &
                                          primary_im_species_names, &
                                          consider_immobile_species,&
                                          option)
  !
  ! Creates a reaction eqaution given a reaction string
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Option_module
  use String_module
  use Input_Aux_module

  character(len=MAXSTRINGLENGTH) :: reaction_string
  PetscInt :: naqcomp ! mobile aqueoues species
  PetscInt :: aq_offset ! offset for aqueous species
  character(len=MAXWORDLENGTH) :: primary_aq_species_names(naqcomp)
  PetscInt :: nimcomp ! immobile primary speces (e.g. biomass)
  PetscInt :: im_offset ! offset for aqueous species
  character(len=MAXWORDLENGTH), pointer :: primary_im_species_names(:)
  PetscBool :: consider_immobile_species
  type(option_type) :: option

  type(reaction_equation_type), pointer :: ReactionEquationCreateFromString

  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word, word2
  PetscInt :: icount
  PetscInt :: midpoint
  PetscInt :: i, j, idum
  PetscReal :: value
  PetscReal :: tempreal
  PetscBool :: negative_flag
  PetscBool :: found
  PetscErrorCode :: ierr
  type(reaction_equation_type), pointer :: reaction_equation


  reaction_equation => ReactionEquationCreate()

  icount = 0
  ! Be sure to copy as words are removed when read.  Need to full string for
  ! later below
  string = reaction_string
  do
    ierr = 0
    call InputReadWord(string,word,PETSC_TRUE,ierr)
    if (InputError(ierr)) exit

    select case(word)
      case('+')
      case('-')
      case('=','<=>','<->','->','=>')
      case default
      ! try reading as double precision
      string2 = word
      if (.not.StringStartsWithAlpha(string2) .and. &
          StringIntegerDoubleOrWord(string2) /= STRING_IS_A_WORD) then
        ! the word is the stoichiometry value
      else
        ! check water
        word2 = 'H2O'
        if (.not.StringCompareIgnoreCase(word,word2)) then
          ! the word is the species name
          icount = icount + 1
        endif
      endif
    end select

  enddo

  ! load species into database format
  reaction_equation%nspec = icount
  allocate(reaction_equation%spec_name(icount))
  reaction_equation%spec_name = ''
  allocate(reaction_equation%stoich(icount))
  reaction_equation%stoich = UNINITIALIZED_DOUBLE
  allocate(reaction_equation%specid(icount))
  reaction_equation%specid = 0

  string = reaction_string
  icount = 1
  ! midpoint points to the first product species, as in
  ! reactant1 + reactant2 <-> product1 + product2
  midpoint = 0
  negative_flag = PETSC_FALSE
  do
    !geh: This conditional ensures that if water is at the end of
    !     the reaction expression, it is skipped.
    if (icount > reaction_equation%nspec) exit

    ierr = 0
    call InputReadWord(string,word,PETSC_TRUE,ierr)
    if (InputError(ierr)) exit

    select case(word)
      case('+')
      case('-')
        ! toggle negative flag
        if (negative_flag) then
          negative_flag = PETSC_FALSE
        else
          negative_flag = PETSC_TRUE
        endif
      case('=','<=>','<->','->','=>')
        midpoint = icount
      case default
        ! try reading as double precision
        string2 = word
        if (.not.StringStartsWithAlpha(string2) .and. &
            StringIntegerDoubleOrWord(string2) /= STRING_IS_A_WORD) then
          i = index(string2,'/')
          if (i > 0) then ! fraction exists
            string2 = word(:i-1)
            call InputReadDouble(string2,option,value,ierr)
            string2 = word(i+1:)
            call InputReadDouble(string2,option,tempreal,ierr)
            value = value / tempreal
          else
            call InputReadDouble(string2,option,value,ierr)
          endif
          if (ierr /= 0) then
            option%io_buffer = 'Keyword "' // trim(word) // &
               '" not recognized in reaction string "' // &
               trim(reaction_string) // '".'
            call PrintErrMsg(option)
          endif
          ! negate if negative stoichiometry
          if (negative_flag) value = -1.0*value
          reaction_equation%stoich(icount) = value
        else
          reaction_equation%spec_name(icount) = word
          if (negative_flag .and. &
              (reaction_equation%stoich(icount) + 999.d0) < 1.d-10) then
            reaction_equation%stoich(icount) = -1.d0
          endif

          ! set the primary aqueous species id
          found = PETSC_FALSE
          do i = 1, naqcomp
            if (StringCompare(word,primary_aq_species_names(i), &
                              MAXWORDLENGTH)) then
              reaction_equation%specid(icount) = i + aq_offset
              found = PETSC_TRUE
              exit
            endif
          enddo
          ! set the primary immobile species id
          if (.not.found .and. consider_immobile_species) then
            do i = 1, nimcomp
              if (StringCompare(word,primary_im_species_names(i), &
                                MAXWORDLENGTH)) then
                reaction_equation%specid(icount) = i + im_offset
                found = PETSC_TRUE
                exit
              endif
            enddo
          endif

          ! check water
          word2 = 'H2O'
          if (StringCompareIgnoreCase(word,word2)) then
            ! set stoichiometry back to uninitialized
            reaction_equation%stoich(icount) = UNINITIALIZED_DOUBLE
            ! don't increment icount
          else if (.not.found) then
            if (consider_immobile_species) then
              option%io_buffer = 'Species ' // trim(word) // &
                        ' in reaction not found among primary species list.'
            else
              option%io_buffer = 'Species ' // trim(word) // &
                ' in reaction not found among primary aqueous species list.'
            endif
            call PrintErrMsg(option)
          else
            icount = icount + 1
          endif
        endif
        negative_flag = PETSC_FALSE
    end select
  enddo

  ! if no stoichiometry specified, default = 1.
  do i = 1, reaction_equation%nspec
    if ((reaction_equation%stoich(i) + 999.d0) < 1.d-10) &
      reaction_equation%stoich(i) = 1.d0
  enddo
  if (midpoint > 0) then
    ! negate stoichiometries after midpoint
    do i = midpoint, reaction_equation%nspec
      reaction_equation%stoich(i) = -1.d0*reaction_equation%stoich(i)
    enddo
  endif
  ! now negate all stoichiometries to have - for reactants; + for products
  do i = 1, reaction_equation%nspec
    reaction_equation%stoich(i) = -1.d0*reaction_equation%stoich(i)
  enddo
  ! reorder species ids in ascending order
  do i = 1, reaction_equation%nspec
    do j = i+1, reaction_equation%nspec
      if (reaction_equation%specid(i) > reaction_equation%specid(j)) then
        ! swap ids
        idum = reaction_equation%specid(j)
        reaction_equation%specid(j) = reaction_equation%specid(i)
        reaction_equation%specid(i) = idum
        ! swap stoichiometry
        value = reaction_equation%stoich(j)
        reaction_equation%stoich(j) = reaction_equation%stoich(i)
        reaction_equation%stoich(i) = value
        ! swap names
        word = reaction_equation%spec_name(j)
        reaction_equation%spec_name(j) = reaction_equation%spec_name(i)
        reaction_equation%spec_name(i) = word
      endif
    enddo
  enddo

  ReactionEquationCreateFromString => reaction_equation

end function ReactionEquationCreateFromString

! ************************************************************************** !

subroutine ReactionEquationSubSpecInRxn(name1,reaction_equation1, &
                                        reaction_equation2,scale)
  !
  ! Swaps out a chemical species in a chemical reaction, replacing it with
  ! the species in a second reaction (swaps 1 into 2)
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use String_module

  character(len=MAXWORDLENGTH) :: name1
  type(reaction_equation_type) :: reaction_equation1
  type(reaction_equation_type) :: reaction_equation2
  PetscReal :: scale

  PetscInt :: i, j, tempcount, prevcount
  character(len=MAXWORDLENGTH) :: tempnames(20)
  PetscReal :: tempstoich(20)
  PetscBool :: found

  tempnames = ''
  tempstoich = 0.d0

  ! load species in reaction other than species 1 into new arrays
  scale = 1.d0
  tempcount = 0
  do i=1,reaction_equation2%nspec
    if (.not.StringCompare(name1, &
                           reaction_equation2%spec_name(i), &
                           MAXWORDLENGTH)) then
      tempcount = tempcount + 1
      tempnames(tempcount) = reaction_equation2%spec_name(i)
      tempstoich(tempcount) = reaction_equation2%stoich(i)
    else
      scale = reaction_equation2%stoich(i)
    endif
  enddo

  ! search for duplicate species and add stoichs or add new species
  ! if not duplicated
  do j=1,reaction_equation1%nspec
    found = PETSC_FALSE
    do i=1,tempcount
      if (StringCompare(tempnames(i), &
                        reaction_equation1%spec_name(j), &
                        MAXWORDLENGTH)) then
        tempstoich(i) = tempstoich(i) + scale*reaction_equation1%stoich(j)
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      tempcount = tempcount + 1
      tempnames(tempcount) = reaction_equation1%spec_name(j)
      tempstoich(tempcount) = scale*reaction_equation1%stoich(j)
    endif
  enddo

  ! deallocate arrays
  deallocate(reaction_equation2%spec_name)
  deallocate(reaction_equation2%stoich)

  ! check for zero stoichiometries due to cancelation
  prevcount = tempcount
  tempcount = 0
  do i=1,prevcount
    if (dabs(tempstoich(i)) > 1.d-10) then
      tempcount = tempcount + 1
      tempnames(tempcount) = tempnames(i)
      tempstoich(tempcount) = tempstoich(i)
    endif
  enddo

  tempnames(tempcount+1:) = ''
  tempstoich(tempcount+1:) = 0.d0

  ! reallocate
  allocate(reaction_equation2%spec_name(tempcount))
  allocate(reaction_equation2%stoich(tempcount))

  ! fill arrays in reaction_equation
  reaction_equation2%nspec = tempcount
  do i=1,tempcount
    reaction_equation2%spec_name(i) = tempnames(i)
    reaction_equation2%stoich(i) = tempstoich(i)
  enddo

end subroutine ReactionEquationSubSpecInRxn

! ************************************************************************** !

recursive subroutine ReactionEquationDestroyRxnPtr(reaction_equation_ptr)
  !
  ! Deallocates a reaction equation pointer
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  implicit none

  type(reaction_equation_ptr_type), pointer :: reaction_equation_ptr

  if (.not.associated(reaction_equation_ptr)) return

  call ReactionEquationDestroyRxnPtr(reaction_equation_ptr%next)
  call ReactionEquationDestroy(reaction_equation_ptr%reaction_equation)

  deallocate(reaction_equation_ptr)
  nullify(reaction_equation_ptr)

end subroutine ReactionEquationDestroyRxnPtr

! ************************************************************************** !

subroutine ReactionEquationDestroy(reaction_equation)
  !
  ! Deallocates a reaction equation
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Utility_module, only : DeallocateArray

  type(reaction_equation_type), pointer :: reaction_equation

  if (.not.associated(reaction_equation)) return

  if (associated(reaction_equation%spec_name)) &
    deallocate(reaction_equation%spec_name)
  nullify(reaction_equation%spec_name)
  call DeallocateArray(reaction_equation%specid)
  call DeallocateArray(reaction_equation%stoich)

  deallocate(reaction_equation)
  nullify(reaction_equation)

end subroutine ReactionEquationDestroy

end module Reaction_Equation_module
