module NWT_Constraint_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use NW_Transport_Aux_module
  use Global_Aux_module
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private
  
  ! concentration subcondition types
  ! constraint on aqueous concentration:
  PetscInt, parameter, public :: CONSTRAINT_AQ_EQUILIBRIUM = 1
  ! constraint on mineral (precipitated) concentration:
  PetscInt, parameter, public :: CONSTRAINT_PPT_EQUILIBRIUM = 2
  ! constraint on sorbed concentration:
  PetscInt, parameter, public :: CONSTRAINT_SB_EQUILIBRIUM = 3
  ! constraint on total bulk concentration:
  PetscInt, parameter, public :: CONSTRAINT_T_EQUILIBRIUM = 4
  ! constraint on mineral volume fraction:
  PetscInt, parameter, public :: CONSTRAINT_MNRL_VOL_FRAC_EQ = 5

  type, public :: nwt_constraint_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name         
    type(nwt_species_constraint_type), pointer :: nwt_species
    type(nwt_constraint_type), pointer :: next    
  end type nwt_constraint_type
  
  type, public :: nwt_constraint_ptr_type
    type(nwt_constraint_type), pointer :: ptr
  end type nwt_constraint_ptr_type
  
  type, public :: nwt_constraint_list_type
    PetscInt :: num_constraints
    type(nwt_constraint_type), pointer :: first
    type(nwt_constraint_type), pointer :: last
    type(nwt_constraint_ptr_type), pointer :: array(:)    
  end type nwt_constraint_list_type
  
  type, public :: nwt_constraint_coupler_type
    character(len=MAXWORDLENGTH) :: constraint_name  
    PetscReal :: time 
    character(len=MAXWORDLENGTH) :: time_units
    type(nwt_species_constraint_type), pointer :: nwt_species
    type(nw_transport_auxvar_type), pointer :: nwt_auxvar
    type(global_auxvar_type), pointer :: global_auxvar
    type(nwt_constraint_coupler_type), pointer :: next   
  end type nwt_constraint_coupler_type
  
  type, public :: nwt_species_constraint_type
    ! Any changes here must be incorporated within NWTProcessConstraint(),
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    PetscInt, pointer :: constraint_type(:)
  end type nwt_species_constraint_type
      
  public :: NWTConstraintAddToList, &
            NWTConstraintInitList, &
            NWTConstraintDestroyList, &
            NWTConstraintGetPtrFromList, &
            NWTConstraintCreate, &
            NWTConstraintDestroy, &
            NWTSpeciesConstraintCreate, &
            NWTSpeciesConstraintDestroy, &
            NWTConstraintCouplerCreate, &
            NWTConstraintCouplerDestroy, &
            NWTConstraintRead, &
            NWTConstraintProcess, &
            NWTConstraintMapToCoupler
            
    
contains

! ************************************************************************** !

function NWTConstraintCreate(option)
  ! 
  ! Creates an NW Transport constraint (set of concentrations
  ! and equilibrium constraints for setting boundary or initial
  ! condition).
  ! 
  ! Author: Jenn Frederick
  ! Date: 06/25/19
  ! 
  
  implicit none
  
  type(option_type) :: option
  type(nwt_constraint_type), pointer :: NWTConstraintCreate
  
  type(nwt_constraint_type), pointer :: constraint
  
  allocate(constraint)
  nullify(constraint%nwt_species)
  nullify(constraint%next)
  constraint%id = 0
  constraint%name = ''
  
  NWTConstraintCreate => constraint

end function NWTConstraintCreate

! ************************************************************************** !

function NWTConstraintCouplerCreate(option)
  ! 
  ! Creates a coupler that ties a constraint to a
  ! transport condition
  ! 
  ! Author: Jenn Frederick
  ! Date: 06/25/2019
  ! 
  
  implicit none
  
  type(option_type) :: option
  type(nwt_constraint_coupler_type), pointer :: NWTConstraintCouplerCreate
  
  type(nwt_constraint_coupler_type), pointer :: coupler
  
  allocate(coupler)
  nullify(coupler%nwt_species)
  nullify(coupler%nwt_auxvar)
  nullify(coupler%global_auxvar)
  
  coupler%time = 0.d0
  coupler%time_units = ''
  
  nullify(coupler%next)
  coupler%constraint_name = ''
  
  NWTConstraintCouplerCreate => coupler

end function NWTConstraintCouplerCreate

! ************************************************************************** !

function NWTSpeciesConstraintCreate(nw_trans,option)
  ! 
  ! Creates a nuclear waste transport species constraint object
  ! 
  ! Author: Jenn Frederick      
  ! Date: 03/21/2019
  ! 
  
  implicit none
  
  type(nw_trans_realization_type) :: nw_trans
  type(option_type) :: option
  type(nwt_species_constraint_type), pointer :: NWTSpeciesConstraintCreate

  type(nwt_species_constraint_type), pointer :: constraint
  
  allocate(constraint)
  allocate(constraint%names(nw_trans%params%nspecies))
  constraint%names = ''
  allocate(constraint%constraint_conc(nw_trans%params%nspecies))
  constraint%constraint_conc = 0.d0
  allocate(constraint%constraint_type(nw_trans%params%nspecies))
  constraint%constraint_type = 0

  NWTSpeciesConstraintCreate => constraint

end function NWTSpeciesConstraintCreate

! ************************************************************************** !

subroutine NWTConstraintRead(constraint,nw_trans,input,option)
  ! 
  ! Reads a transport constraint from the input file
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/29/2019
  ! 

  use Input_Aux_module
  use String_module
  use Logging_module

  implicit none
  
  type(nwt_constraint_type) :: constraint
  type(nw_trans_realization_type) :: nw_trans
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: block_string
  PetscInt :: icomp
  type(nwt_species_constraint_type), pointer :: nwt_species_constraint
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_tran_constraint_read, &
                          ierr);CHKERRQ(ierr)

  ! read the constraint
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONSTRAINT')
        
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONSTRAINT')   
      
    select case(trim(word))

      case('CONC','CONCENTRATIONS')

        nwt_species_constraint => &
          NWTSpeciesConstraintCreate(nw_trans,option)

        block_string = 'CONSTRAINT, CONCENTRATIONS'
        icomp = 0
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)
          
          if (InputCheckExit(input,option)) exit  
          
          icomp = icomp + 1        
          
          if (icomp > nw_trans%params%nspecies) then
            option%io_buffer = 'Number of concentration constraints exceeds &
                               &the number of species given in the &
                               &SUBSURFACE_NUCLEAR_WASTE_TRANSPORT block. &
                               &Error in constraint: ' // trim(constraint%name)
            call PrintErrMsg(option)
          endif
          
          call InputReadWord(input,option,nwt_species_constraint%names(icomp), &
                          PETSC_TRUE)
          call InputErrorMsg(input,option,'species name',block_string)
          option%io_buffer = 'Constraint Species: ' // &
                             trim(nwt_species_constraint%names(icomp))
          call PrintMsg(option)
          
          call InputReadDouble(input,option, &
                               nwt_species_constraint%constraint_conc(icomp))
          call InputErrorMsg(input,option,'concentration',block_string)
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputDefaultMsg(input,option,trim(block_string) // &
                               'constraint type')
          if (len_trim(word) > 0) then
            call StringToUpper(word)
            select case(word)
              case('AQ','AQUEOUS')
                nwt_species_constraint%constraint_type(icomp) = &
                                                    CONSTRAINT_AQ_EQUILIBRIUM
              case('PPT','PRECIPITATED')
                nwt_species_constraint%constraint_type(icomp) = &
                                                    CONSTRAINT_PPT_EQUILIBRIUM
              case('SB','SORBED')
                nwt_species_constraint%constraint_type(icomp) = &
                                                    CONSTRAINT_SB_EQUILIBRIUM
              case('T','TOTAL')
                nwt_species_constraint%constraint_type(icomp) = &
                                                    CONSTRAINT_T_EQUILIBRIUM
              case('VF','PRECIPITATED_VOLUME_FRACTION')
                nwt_species_constraint%constraint_type(icomp) = &
                                                    CONSTRAINT_MNRL_VOL_FRAC_EQ
              case default
                option%io_buffer = 'Error in constraint: ' // &
                  trim(constraint%name) // '. The constraint type given for &
                  &species ' // trim(nwt_species_constraint%names(icomp)) // &
                  ' is not recognized: ' // trim(word) // '. &
                  &Options include: VF, T, AQ, PPT, or SB only.'
                call PrintErrMsg(option)
            end select
          else
            option%io_buffer = 'Error in constraint: ' // &
              trim(constraint%name) // '. A constraint type was not specified &
              &for species ' // trim(nwt_species_constraint%names(icomp)) // '.'
            call PrintErrMsg(option)
          endif
        enddo  
        
        if (icomp < nw_trans%params%nspecies) then
          option%io_buffer = &
                   'Number of concentration constraints is less than ' // &
                   'number of species in species constraint.'
          call PrintErrMsg(option)
        endif
        if (icomp > nw_trans%params%nspecies) then
          option%io_buffer = &
                   'Number of concentration constraints is greater than ' // &
                   'number of species in species constraint.'
          call PrintWrnMsg(option)
        endif
        
        if (associated(constraint%nwt_species)) &
          call NWTSpeciesConstraintDestroy(constraint%nwt_species)
        constraint%nwt_species => nwt_species_constraint 
               
        
      case default
        call InputKeywordUnrecognized(word,'CONSTRAINT',option)
    end select 
  
  enddo  
  
  call PetscLogEventEnd(logging%event_tran_constraint_read,ierr);CHKERRQ(ierr)

end subroutine NWTConstraintRead

! ************************************************************************** ! 

subroutine NWTConstraintProcess(nw_trans,constraint_name, &
                                nwt_species_constraint,option)
  ! 
  ! Ensures ordering of species is consistant between the nw_trans object
  ! and the constraint object. 
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/22/2019
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use String_module
  use Utility_module
  
  implicit none
  
  type(nw_trans_realization_type), pointer :: nw_trans
  character(len=MAXWORDLENGTH) :: constraint_name
  type(nwt_species_constraint_type), pointer :: nwt_species_constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: ispecies, jspecies
  PetscReal :: constraint_conc(nw_trans%params%nspecies)
  PetscInt :: constraint_type(nw_trans%params%nspecies)
  character(len=MAXWORDLENGTH) :: constraint_species_names( &
                                                     nw_trans%params%nspecies)
  
  constraint_conc = 0.d0
  constraint_type = 0
  constraint_species_names = ''
  
  do ispecies = 1, nw_trans%params%nspecies
    found = PETSC_FALSE
    do jspecies = 1, nw_trans%params%nspecies
      if (StringCompare(nwt_species_constraint%names(ispecies), &
                        nw_trans%species_names(jspecies),MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
               'Species ' // trim(nwt_species_constraint%names(ispecies)) // &
               ' from CONSTRAINT ' // trim(constraint_name) // &
               ' not found among species.'
      call PrintErrMsg(option)
    else
      constraint_conc(jspecies) = &
                               nwt_species_constraint%constraint_conc(ispecies)
      constraint_type(jspecies) = &
                               nwt_species_constraint%constraint_type(ispecies)
      constraint_species_names(jspecies) = &
                                         nwt_species_constraint%names(ispecies)
    endif
  enddo
  
  ! place ordered constraint parameters back in original arrays
  nwt_species_constraint%constraint_conc = constraint_conc
  nwt_species_constraint%constraint_type = constraint_type
  nwt_species_constraint%names = constraint_species_names
    
end subroutine NWTConstraintProcess


! ************************************************************************** !

subroutine NWTConstraintInitList(list)
  ! 
  ! Initializes a transport constraint list
  ! 
  ! Author: Jenn Frederick
  ! Date: 06/25/2019
  ! 

  implicit none

  type(nwt_constraint_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_constraints = 0

end subroutine NWTConstraintInitList

! ************************************************************************** !

subroutine NWTConstraintAddToList(new_constraint,list)
  ! 
  ! Adds a new constraint to a transport constraint
  ! list
  ! 
  ! Author: Jenn Frederick
  ! Date: 06/25/2019
  ! 

  implicit none
  
  type(nwt_constraint_type), pointer :: new_constraint
  type(nwt_constraint_list_type) :: list
  
  list%num_constraints = list%num_constraints + 1
  new_constraint%id = list%num_constraints
  if (.not.associated(list%first)) list%first => new_constraint
  if (associated(list%last)) list%last%next => new_constraint
  list%last => new_constraint
  
end subroutine NWTConstraintAddToList

! ************************************************************************** !

function NWTConstraintGetPtrFromList(constraint_name,constraint_list)
  ! 
  ! Returns a pointer to the constraint matching
  ! constraint_name
  ! 
  ! Author: Jenn Frederick
  ! Date: 06/25/2019
  ! 

  use String_module

  implicit none
  
  type(nwt_constraint_type), pointer :: NWTConstraintGetPtrFromList
  character(len=MAXWORDLENGTH) :: constraint_name
  type(nwt_constraint_list_type) :: constraint_list
 
  PetscInt :: length
  type(nwt_constraint_type), pointer :: constraint
    
  nullify(NWTConstraintGetPtrFromList)
  constraint => constraint_list%first
  
  do 
    if (.not.associated(constraint)) exit
    length = len_trim(constraint_name)
    if (length == len_trim(constraint%name) .and. &
        StringCompare(constraint%name,constraint_name, &
                        length)) then
      NWTConstraintGetPtrFromList => constraint
      return
    endif
    constraint => constraint%next
  enddo
  
end function NWTConstraintGetPtrFromList

! ************************************************************************** !

subroutine NWTConstraintMapToCoupler(constraint_coupler,constraint)
  ! 
  ! Maps members of a constraint object to a constraint coupler.
  ! 
  ! Author: Jenn Frederick
  ! Date: 06/25/2019
  !  

  implicit none

  type(nwt_constraint_coupler_type) :: constraint_coupler
  type(nwt_constraint_type) :: constraint
  
  constraint_coupler%nwt_species => constraint%nwt_species

end subroutine NWTConstraintMapToCoupler

! ************************************************************************** !

subroutine NWTConstraintDestroy(constraint)
  ! 
  ! Deallocates a constraint
  ! 
  ! Author: Jenn Frederick
  ! Date: 06/25/2019
  ! 
  
  use Utility_module, only: DeallocateArray

  implicit none
  
  type(nwt_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return

  if (associated(constraint%nwt_species)) then
    call DeallocateArray(constraint%nwt_species%names)
    call DeallocateArray(constraint%nwt_species%constraint_conc)
    call DeallocateArray(constraint%nwt_species%constraint_type)
  endif
  nullify(constraint%nwt_species)
 
  deallocate(constraint)
  nullify(constraint)

end subroutine NWTConstraintDestroy

! ************************************************************************** !

subroutine NWTConstraintDestroyList(constraint_list)
  ! 
  ! Deallocates a list of constraints
  ! 
  ! Author: Jenn Frederick
  ! Date: 06/25/2019
  !  

  implicit none
  
  type(nwt_constraint_list_type), pointer :: constraint_list
  
  type(nwt_constraint_type), pointer :: constraint, prev_constraint
  
  if (.not.associated(constraint_list)) return
  
  constraint => constraint_list%first
  do 
    if (.not.associated(constraint)) exit
    prev_constraint => constraint
    constraint => constraint%next
    call NWTConstraintDestroy(prev_constraint)
  enddo
  
  constraint_list%num_constraints = 0
  nullify(constraint_list%first)
  nullify(constraint_list%last)
  if (associated(constraint_list%array)) deallocate(constraint_list%array)
  nullify(constraint_list%array)
  
  deallocate(constraint_list)
  nullify(constraint_list)

end subroutine NWTConstraintDestroyList

! ************************************************************************** !

subroutine NWTConstraintCouplerDestroy(coupler_list)
  ! 
  ! Destroys a constraint coupler linked list
  ! 
  ! Author: Jenn Frederick
  ! Date: 06/25/2019
  !  

  use Option_module
  
  implicit none
  
  type(nwt_constraint_coupler_type), pointer :: coupler_list
  
  type(nwt_constraint_coupler_type), pointer :: cur_coupler, prev_coupler
  
  cur_coupler => coupler_list
  
  do
    if (.not.associated(cur_coupler)) exit
    prev_coupler => cur_coupler
    cur_coupler => cur_coupler%next
    if (associated(prev_coupler%nwt_auxvar)) then
      call NWTAuxVarDestroy(prev_coupler%nwt_auxvar)
    endif
    nullify(prev_coupler%nwt_auxvar)
    if (associated(prev_coupler%global_auxvar)) then
      call GlobalAuxVarDestroy(prev_coupler%global_auxvar)
    endif
    nullify(prev_coupler%global_auxvar)
    nullify(prev_coupler%nwt_species)
    nullify(prev_coupler%next)
    deallocate(prev_coupler)
    nullify(prev_coupler)
  enddo
  
  nullify(coupler_list)
  
end subroutine NWTConstraintCouplerDestroy

! ************************************************************************** !

subroutine NWTSpeciesConstraintDestroy(constraint)
  ! 
  ! Deallocates a nuclear waste transport species constraint object
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/21/2019
  ! 

  use Utility_module, only: DeallocateArray
  
  implicit none
  
  type(nwt_species_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_conc)
  call DeallocateArray(constraint%constraint_type)

  deallocate(constraint)
  nullify(constraint)

end subroutine NWTSpeciesConstraintDestroy

! ************************************************************************** !

end module NWT_Constraint_module
