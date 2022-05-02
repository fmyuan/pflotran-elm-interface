module Transport_Constraint_NWT_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Transport_Constraint_Base_module
  use NW_Transport_Aux_module

  implicit none

  private

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

  type, public, extends(tran_constraint_base_type) :: tran_constraint_nwt_type
    type(nwt_species_constraint_type), pointer :: nwt_species
  contains
    procedure, public :: Strip => TranConstraintNWTStrip
  end type tran_constraint_nwt_type

  type, public, extends(tran_constraint_coupler_base_type) :: &
                                           tran_constraint_coupler_nwt_type
   type(nw_transport_auxvar_type), pointer :: nwt_auxvar
  end type tran_constraint_coupler_nwt_type

  public :: TranConstraintNWTCreate, &
            TranConstraintNWTCast, &
            TranConstraintCouplerNWTCast, &
            TranConstraintCouplerNWTCreate, &
            TranConstraintNWTGetAuxVar, &
            TranConstraintNWTRead, &
            NWTSpeciesConstraintCreate, &
            NWTConstraintProcess

contains

! ************************************************************************** !

function TranConstraintNWTCreate(option)
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
  class(tran_constraint_nwt_type), pointer :: TranConstraintNWTCreate

  class(tran_constraint_nwt_type), pointer :: constraint

  allocate(constraint)
  call TranConstraintBaseInit(constraint,option)
  nullify(constraint%nwt_species)

  TranConstraintNWTCreate => constraint

end function TranConstraintNWTCreate

! ************************************************************************** !

function TranConstraintCouplerNWTCreate(option)
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
  class(tran_constraint_coupler_nwt_type), pointer :: &
                                                TranConstraintCouplerNWTCreate

  class(tran_constraint_coupler_nwt_type), pointer :: coupler

  allocate(coupler)
  call TranConstraintCouplerBaseInit(coupler,option)
  nullify(coupler%nwt_auxvar)

  TranConstraintCouplerNWTCreate => coupler

end function TranConstraintCouplerNWTCreate


! ************************************************************************** !

function TranConstraintNWTCast(this)
  !
  ! Casts a tran_constraint_base_type to tran_constraint_nwt_type
  !
  ! Author: Glenn Hammond
  ! Date: 10/07/19
  !

  implicit none

  class(tran_constraint_base_type), pointer :: this

  class(tran_constraint_nwt_type), pointer :: TranConstraintNWTCast

  nullify(TranConstraintNWTCast)
  if (.not.associated(this)) return
  select type (this)
    class is (tran_constraint_nwt_type)
      TranConstraintNWTCast => this
  end select

end function TranConstraintNWTCast

! ************************************************************************** !

function TranConstraintCouplerNWTCast(this)
  !
  ! Casts a tran_constraint_coupler_base_type to
  ! tran_constraint_coupler_nwt_type
  !
  ! Author: Glenn Hammond
  ! Date: 10/07/19
  !

  implicit none

  class(tran_constraint_coupler_base_type), pointer :: this

  class(tran_constraint_coupler_nwt_type), pointer :: &
                     TranConstraintCouplerNWTCast

  nullify(TranConstraintCouplerNWTCast)
  if (.not.associated(this)) return
  select type (this)
    class is (tran_constraint_coupler_nwt_type)
      TranConstraintCouplerNWTCast => this
  end select

end function TranConstraintCouplerNWTCast

! ************************************************************************** !

function TranConstraintNWTGetAuxVar(this)
  !
  ! Returns the auxvar associated with tran_constraint_coupler_nwt_type
  !
  ! Author: Glenn Hammond
  ! Date: 10/07/19
  !

  implicit none

  class(tran_constraint_coupler_base_type), pointer :: this

  type(nw_transport_auxvar_type), pointer :: TranConstraintNWTGetAuxVar

  nullify(TranConstraintNWTGetAuxVar)
  if (.not.associated(this)) return
  select type (coupler=>this)
    class is (tran_constraint_coupler_nwt_type)
      TranConstraintNWTGetAuxVar => coupler%nwt_auxvar
  end select

end function TranConstraintNWTGetAuxVar

! ************************************************************************** !

subroutine TranConstraintNWTRead(constraint,reaction_nw,input,option)
  !
  ! Reads a transport constraint from the input file
  !
  ! Author: Jenn Frederick
  ! Date: 05/29/2019
  !
  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use NW_Transport_Aux_module

  implicit none

  class(tran_constraint_nwt_type) :: constraint
  class(reaction_nw_type) :: reaction_nw
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
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONSTRAINT')

    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONSTRAINT')

    select case(trim(word))

      case('CONC','CONCENTRATIONS')

        nwt_species_constraint => &
          NWTSpeciesConstraintCreate(reaction_nw,option)

        block_string = 'CONSTRAINT, CONCENTRATIONS'
        icomp = 0
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)

          if (InputCheckExit(input,option)) exit

          icomp = icomp + 1

          if (icomp > reaction_nw%params%nspecies) then
            option%io_buffer = 'Number of concentration constraints exceeds &
                               &the number of species given in the &
                               &NUCLEAR_WASTE_CHEMISTRY block. &
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

        if (icomp < reaction_nw%params%nspecies) then
          option%io_buffer = &
                   'Number of concentration constraints is less than ' // &
                   'number of species in species constraint.'
          call PrintErrMsg(option)
        endif
        if (icomp > reaction_nw%params%nspecies) then
          option%io_buffer = &
                   'Number of concentration constraints is greater than ' // &
                   'number of species in species constraint.'
          call PrintWrnMsg(option)
        endif

        if (associated(constraint%nwt_species)) &
          call NWTSpeciesConstraintDestroy(constraint%nwt_species)
        constraint%nwt_species => nwt_species_constraint


      case default
        call InputKeywordUnrecognized(input,word,'CONSTRAINT',option)
    end select

  enddo
  call InputPopBlock(input,option)

  call PetscLogEventEnd(logging%event_tran_constraint_read,ierr);CHKERRQ(ierr)

end subroutine TranConstraintNWTRead

! ************************************************************************** !

subroutine NWTConstraintProcess(reaction_nw,constraint,option)
  !
  ! Ensures ordering of species is consistant between the reaction_nw object
  ! and the constraint object.
  !
  ! Author: Jenn Frederick
  ! Date: 03/22/2019
  !
  use Option_module
  use String_module
  use Utility_module

  implicit none

  class(reaction_nw_type), pointer :: reaction_nw
  class(tran_constraint_nwt_type) :: constraint
  type(option_type) :: option

  PetscBool :: found
  PetscInt :: ispecies, jspecies
  PetscReal :: constraint_conc(reaction_nw%params%nspecies)
  PetscInt :: constraint_type(reaction_nw%params%nspecies)
  type(nwt_species_constraint_type), pointer :: nwt_species_constraint
  character(len=MAXWORDLENGTH) :: constraint_species_names( &
                                                     reaction_nw%params%nspecies)

  constraint_conc = 0.d0
  constraint_type = 0
  constraint_species_names = ''

  nwt_species_constraint => constraint%nwt_species

  do ispecies = 1, reaction_nw%params%nspecies
    found = PETSC_FALSE
    do jspecies = 1, reaction_nw%params%nspecies
      if (StringCompare(nwt_species_constraint%names(ispecies), &
                        reaction_nw%species_names(jspecies),MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
               'Species ' // trim(nwt_species_constraint%names(ispecies)) // &
               ' from CONSTRAINT ' // trim(constraint%name) // &
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

subroutine TranConstraintNWTStrip(this)
  !
  ! Deallocates a constraint
  !
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  !

  implicit none

  class(tran_constraint_nwt_type) :: this

  call TranConstraintBaseStrip(this)

  if (associated(this%nwt_species)) &
    call NWTSpeciesConstraintDestroy(this%nwt_species)

end subroutine TranConstraintNWTStrip

! ************************************************************************** !

subroutine TranConstraintCouplerNWTStrip(this)
  !
  ! Deallocate dynamic members of class
  !
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  !
  implicit none

  class(tran_constraint_coupler_nwt_type) :: this

  call TranConstraintCouplerBaseStrip(this)

  call NWTAuxVarDestroy(this%nwt_auxvar)

end subroutine TranConstraintCouplerNWTStrip

end module Transport_Constraint_NWT_module
