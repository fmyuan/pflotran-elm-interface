module Transport_Constraint_RT_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Transport_Constraint_Base_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module

  use Reaction_Surface_Complexation_Aux_module
  use Reaction_Mineral_Aux_module
  use Reaction_Immobile_Aux_module


  implicit none

  private

  ! concentration subcondition types
  PetscInt, parameter, public :: CONSTRAINT_NULL = 0
  PetscInt, parameter, public :: CONSTRAINT_FREE = 1
  PetscInt, parameter, public :: CONSTRAINT_TOTAL = 2
  PetscInt, parameter, public :: CONSTRAINT_LOG = 3
  PetscInt, parameter, public :: CONSTRAINT_PH = 4
  PetscInt, parameter, public :: CONSTRAINT_PE = 5
  PetscInt, parameter, public :: CONSTRAINT_EH = 6
  PetscInt, parameter, public :: CONSTRAINT_MINERAL = 7
  PetscInt, parameter, public :: CONSTRAINT_GAS = 8
  PetscInt, parameter, public :: CONSTRAINT_CHARGE_BAL = 9
  PetscInt, parameter, public :: CONSTRAINT_TOTAL_SORB = 10
  PetscInt, parameter, public :: CONSTRAINT_SUPERCRIT_CO2 = 11
  PetscInt, parameter, public :: CONSTRAINT_TOTAL_AQ_PLUS_SORB = 12

  type, public, extends(tran_constraint_base_type) :: tran_constraint_rt_type
    type(aq_species_constraint_type), pointer :: aqueous_species
    type(guess_constraint_type), pointer :: free_ion_guess
    type(mineral_constraint_type), pointer :: minerals
    type(srfcplx_constraint_type), pointer :: surface_complexes
    type(immobile_constraint_type), pointer :: immobile_species
  contains
    procedure, public :: Strip => TranConstraintRTStrip
  end type tran_constraint_rt_type

  type, public, extends(tran_constraint_coupler_base_type) :: &
                                           tran_constraint_coupler_rt_type
    type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  contains
    procedure, public :: Strip => TranConstraintCouplerRTStrip
  end type tran_constraint_coupler_rt_type

  public :: TranConstraintRTCreate, &
            TranConstraintRTCast, &
            TranConstraintCouplerRTCast, &
            TranConstraintCouplerRTCreate, &
            TranConstraintRTGetAuxVar, &
            TranConstraintRTRead

contains

! ************************************************************************** !

function TranConstraintRTCreate(option)
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
  class(tran_constraint_rt_type), pointer :: TranConstraintRTCreate

  class(tran_constraint_rt_type), pointer :: constraint

  allocate(constraint)
  call TranConstraintBaseInit(constraint,option)
  nullify(constraint%aqueous_species)
  nullify(constraint%free_ion_guess)
  nullify(constraint%minerals)
  nullify(constraint%surface_complexes)
  nullify(constraint%immobile_species)

  TranConstraintRTCreate => constraint

end function TranConstraintRTCreate

! ************************************************************************** !

function TranConstraintCouplerRTCreate(option)
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
  class(tran_constraint_coupler_rt_type), pointer :: &
                                                TranConstraintCouplerRTCreate

  class(tran_constraint_coupler_rt_type), pointer :: coupler

  allocate(coupler)
  call TranConstraintCouplerBaseInit(coupler,option)
  nullify(coupler%rt_auxvar)

  TranConstraintCouplerRTCreate => coupler

end function TranConstraintCouplerRTCreate

! ************************************************************************** !

function TranConstraintRTCast(this)
  !
  ! Casts a tran_constraint_base_type to tran_constraint_rt_type
  !
  ! Author: Glenn Hammond
  ! Date: 10/07/19
  !

  implicit none

  class(tran_constraint_base_type), pointer :: this

  class(tran_constraint_rt_type), pointer :: TranConstraintRTCast

  nullify(TranConstraintRTCast)
  if (.not.associated(this)) return
  select type (this)
    class is (tran_constraint_rt_type)
      TranConstraintRTCast => this
  end select

end function TranConstraintRTCast

! ************************************************************************** !

function TranConstraintCouplerRTCast(this)
  !
  ! Casts a tran_constraint_coupler_base_type to
  ! tran_constraint_coupler_rt_type
  !
  ! Author: Glenn Hammond
  ! Date: 10/07/19
  !

  implicit none

  class(tran_constraint_coupler_base_type), pointer :: this

  class(tran_constraint_coupler_rt_type), pointer :: &
                     TranConstraintCouplerRTCast

  nullify(TranConstraintCouplerRTCast)
  if (.not.associated(this)) return
  select type (this)
    class is (tran_constraint_coupler_rt_type)
      TranConstraintCouplerRTCast => this
  end select

end function TranConstraintCouplerRTCast

! ************************************************************************** !

function TranConstraintRTGetAuxVar(this)
  !
  ! Returns the auxvar associated with tran_constraint_coupler_rt_type
  !
  ! Author: Glenn Hammond
  ! Date: 10/07/19
  !

  implicit none

  class(tran_constraint_coupler_base_type), pointer :: this

  type(reactive_transport_auxvar_type), pointer :: TranConstraintRTGetAuxVar

  nullify(TranConstraintRTGetAuxVar)
  if (.not.associated(this)) return
  select type (coupler=>this)
    class is (tran_constraint_coupler_rt_type)
      TranConstraintRTGetAuxVar => coupler%rt_auxvar
  end select

end function TranConstraintRTGetAuxVar

! ************************************************************************** !

subroutine TranConstraintRTRead(constraint,reaction,input,option)
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

  class(tran_constraint_rt_type) :: constraint
  class(reaction_rt_type) :: reaction
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: block_string
  PetscInt :: icomp, imnrl, iimmobile
  PetscInt :: jcomp, jmnrl
  PetscInt :: tempint
  PetscInt :: isrfcplx
  PetscInt :: length
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(guess_constraint_type), pointer :: free_ion_guess_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint
  type(srfcplx_constraint_type), pointer :: srfcplx_constraint
  type(immobile_constraint_type), pointer :: immobile_constraint
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

      case('CONC','CONCENTRATIONS')

        aq_species_constraint => &
          ReactionAuxCreateSpecConstraint(reaction,option)

        block_string = 'CONSTRAINT, CONCENTRATIONS'
        icomp = 0
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)

          if (InputCheckExit(input,option)) exit

          icomp = icomp + 1

          if (icomp > reaction%naqcomp) then
            option%io_buffer = 'Number of concentration constraints exceeds &
              &number of primary chemical components in constraint: ' // &
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

          call InputReadCard(input,option,word)
          call InputDefaultMsg(input,option, &
                               trim(block_string) // ' constraint_type')
          length = len_trim(word)
          if (length > 0) then
            call StringToUpper(word)
            select case(word)
              case('F','FREE')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_FREE
              case('T','TOTAL')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_TOTAL
              case('TOTAL_SORB')
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_TOTAL_SORB
              case('TOTAL_AQ_PLUS_SORB')
                if (reaction%neqsorb == 0) then
                  option%io_buffer = 'TOTAL_AQ_PLUS_SORB constraint may not &
                    &be used unless equilibrium sorption is employed within &
                    &the CHEMISTRY block.'
                  call PrintErrMsg(option)
                endif
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_TOTAL_AQ_PLUS_SORB
              case('S')
                option%io_buffer = '"S" constraint type no longer &
                  &supported as of March 4, 2013.'
                call PrintErrMsg(option)
              case('P','PH')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_PH
              case('E','PE')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_PE
              case('L','LOG')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_LOG
              case('M','MINERAL','MNRL')
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_MINERAL
              case('G','GAS')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_GAS
              case('SC','CONSTRAINT_SUPERCRIT_CO2')
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_SUPERCRIT_CO2
              case('Z','CHARGE_BALANCE')
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_CHARGE_BAL
              case default
                call InputKeywordUnrecognized(input,word, &
                       'CONSTRAINT,CONCENTRATION,TYPE',option)
            end select

            if (aq_species_constraint%constraint_type(icomp) == &
                  CONSTRAINT_MINERAL .or. &
                aq_species_constraint%constraint_type(icomp) == &
                  CONSTRAINT_GAS .or.&
                aq_species_constraint%constraint_type(icomp) == &
                  CONSTRAINT_SUPERCRIT_CO2) then
              call InputReadCard(input,option,aq_species_constraint% &
                                 constraint_aux_string(icomp))
              call InputErrorMsg(input,option,'constraining species name', &
                                 block_string)
            else
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
            endif
          else
            option%io_buffer = 'A constraint type (e.g. T, F, P, etc.) is &
              &missing for primary species "' // &
              trim(aq_species_constraint%names(icomp)) // &
              '" in constraint "' // &
              trim(constraint%name) // &
              '".'
            call PrintErrMsg(option)
          endif

        enddo
        call InputPopBlock(input,option)

        if (icomp < reaction%naqcomp) then
          option%io_buffer = &
                   'Number of concentration constraints is less than &
                   &number of primary species in aqueous constraint.'
          call PrintErrMsg(option)
        endif
        if (icomp > reaction%naqcomp) then
          option%io_buffer = &
                   'Number of concentration constraints is greater than &
                   &number of primary species in aqueous constraint.'
          call PrintErrMsg(option)
        endif

        do icomp = 1, reaction%naqcomp
          tempint = 0
          do jcomp = 1, reaction%naqcomp
            if (StringCompare(aq_species_constraint%names(icomp), &
                              aq_species_constraint%names(jcomp), &
                              MAXWORDLENGTH)) tempint = tempint + 1
          enddo
          if (tempint > 1) then
            option%io_buffer = 'Duplicated primary species in concentration &
              &constraint: ' // trim(aq_species_constraint%names(icomp))
            call PrintErrMsg(option)
          endif
        enddo

        if (associated(constraint%aqueous_species)) &
          call ReactionAuxDestroySpecConstraint(constraint%aqueous_species)
        constraint%aqueous_species => aq_species_constraint

      case('FREE_ION_GUESS')

        free_ion_guess_constraint => &
          ReactionAuxCreateGuessConstraint(reaction,option)

        block_string = 'CONSTRAINT, FREE_ION_GUESS'
        icomp = 0
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)

          if (InputCheckExit(input,option)) exit

          icomp = icomp + 1

          if (icomp > reaction%naqcomp) then
            option%io_buffer = 'Number of free ion guess constraints &
                               &exceeds number of primary chemical &
                               &components in constraint: ' // &
                               trim(constraint%name)
            call PrintErrMsg(option)
          endif

          call InputReadCard(input,option, &
                             free_ion_guess_constraint%names(icomp))
          call InputErrorMsg(input,option,'free ion guess name',block_string)
          option%io_buffer = 'Constraint Species: ' // &
                             trim(free_ion_guess_constraint%names(icomp))
          call PrintMsg(option)

          call InputReadDouble(input,option,free_ion_guess_constraint%conc(icomp))
          call InputErrorMsg(input,option,'free ion guess',block_string)
        enddo
        call InputPopBlock(input,option)

        if (icomp < reaction%naqcomp) then
          option%io_buffer = &
                   'Number of free ion guess constraints is less than &
                   &number of primary species in aqueous constraint.'
          call PrintErrMsg(option)
        endif
        if (icomp > reaction%naqcomp) then
          option%io_buffer = &
                   'Number of free ion guess constraints is greater than &
                   &number of primary species in aqueous constraint.'
          call PrintErrMsg(option)
        endif

        do icomp = 1, reaction%naqcomp
          tempint = 0
          do jcomp = 1, reaction%naqcomp
            if (StringCompare(free_ion_guess_constraint%names(icomp), &
                              free_ion_guess_constraint%names(jcomp), &
                              MAXWORDLENGTH)) tempint = tempint + 1
          enddo
          if (tempint > 1) then
            option%io_buffer = 'Duplicated primary species in free ion &
              &guess constraint: ' // trim(aq_species_constraint%names(icomp))
            call PrintErrMsg(option)
          endif
        enddo

        if (associated(constraint%free_ion_guess)) &
          call ReactionAuxDestroyGuesConstraint(constraint%free_ion_guess)
        constraint%free_ion_guess => free_ion_guess_constraint
        nullify(free_ion_guess_constraint)

      case('MNRL','MINERALS')

        mineral_constraint => &
          ReactionMnrlCreateMnrlConstraint(reaction%mineral,option)

        block_string = 'CONSTRAINT, MINERALS'
        imnrl = 0
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)

          if (InputCheckExit(input,option)) exit

          imnrl = imnrl + 1

          if (imnrl > reaction%mineral%nkinmnrl) then
            option%io_buffer = &
                     'Number of mineral constraints exceeds number of &
                     &kinetic minerals in constraint: ' // &
                      trim(constraint%name)
            call PrintErrMsg(option)
          endif

          call InputReadCard(input,option,mineral_constraint%names(imnrl))
          call InputErrorMsg(input,option,'mineral name',block_string)
          option%io_buffer = 'Constraint Minerals: ' // &
                             trim(mineral_constraint%names(imnrl))
          call PrintMsg(option)

          ! volume fraction
          string = trim(input%buf)
          call InputReadWord(string,word,PETSC_TRUE,ierr)
          ! if a dataset
          if (StringCompareIgnoreCase(word,'DATASET')) then
            call InputPushCard(input,word,option)
            input%buf = trim(string)
            call InputReadWord(input,option,mineral_constraint% &
                                constraint_vol_frac_string(imnrl),PETSC_TRUE)
            call InputErrorMsg(input,option,'dataset name', &
                               trim(block_string) // ' VOL FRAC')
            mineral_constraint%external_voL_frac_dataset(imnrl) = PETSC_TRUE
          else
            call InputReadDouble(input,option, &
                                 mineral_constraint%constraint_vol_frac(imnrl))
            call InputErrorMsg(input,option,'volume fraction',block_string)
          endif

          string = trim(input%buf)
          call InputReadWord(string,word,PETSC_TRUE,ierr)
          ! if a dataset
          if (StringCompareIgnoreCase(word,'DATASET')) then
            call InputPushCard(input,word,option)
            input%buf = trim(string)
            call InputReadWord(input,option,mineral_constraint% &
                                constraint_area_string(imnrl),PETSC_TRUE)
            call InputErrorMsg(input,option,'dataset name', &
                               trim(block_string) // ' SURF AREA')
            mineral_constraint%external_area_dataset(imnrl) = PETSC_TRUE
          else
            ! specific surface area
            call InputReadDouble(input,option, &
                                 mineral_constraint%constraint_area(imnrl))
            call InputErrorMsg(input,option,'surface area',block_string)
          endif
          ! read units if they exist. conversion takes place later
          call InputReadWord(input,option,mineral_constraint% &
                             constraint_area_units(imnrl),PETSC_TRUE)
          if (InputError(input)) then
            mineral_constraint%constraint_area_units(imnrl) = 'm^2/m^3'
            input%err_buf = trim(mineral_constraint%names(imnrl)) // &
                             ' SPECIFIC SURFACE_AREA UNITS'
            call InputDefaultMsg(input,option)
          endif
        enddo
        call InputPopBlock(input,option)

        if (imnrl < reaction%mineral%nkinmnrl) then
          option%io_buffer = &
                   'Mineral lists in constraints must provide a volume &
                   &fraction and surface area for all kinetic minerals &
                   &(listed under MINERAL_KINETICS card in CHEMISTRY), &
                   &regardless of whether or not they are present (just &
                   &assign a zero volume fraction if not present).'
          call PrintErrMsg(option)
        endif

        do imnrl = 1, reaction%mineral%nkinmnrl
          tempint = 0
          do jmnrl = 1, reaction%mineral%nkinmnrl
            if (StringCompare(mineral_constraint%names(imnrl), &
                              mineral_constraint%names(jmnrl), &
                              MAXWORDLENGTH)) tempint = tempint + 1
          enddo
          if (tempint > 1) then
            option%io_buffer = 'Duplicated minerals in mineral &
              &constraint: ' // trim(mineral_constraint%names(imnrl))
            call PrintErrMsg(option)
          endif
        enddo

        if (associated(constraint%minerals)) then
          call ReactionMnrlDestMnrlConstraint(constraint%minerals)
        endif
        constraint%minerals => mineral_constraint

      case('SURFACE_COMPLEXES')

        srfcplx_constraint => &
          ReactionSrfCplxCreateConstraint(reaction%surface_complexation,option)

        block_string = 'CONSTRAINT, SURFACE_COMPLEXES'
        isrfcplx = 0
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)

          if (InputCheckExit(input,option)) exit

          isrfcplx = isrfcplx + 1

          if (isrfcplx > reaction%surface_complexation%nkinsrfcplx) then
            option%io_buffer = &
                     'Number of surface complex constraints exceeds &
                     &number of kinetic surface complexes in constraint: ' // &
                      trim(constraint%name)
            call PrintErrMsg(option)
          endif

          call InputReadCard(input,option,srfcplx_constraint%names(isrfcplx))
          call InputErrorMsg(input,option,'surface complex name',block_string)
          option%io_buffer = 'Constraint Surface Complex: ' // &
                             trim(srfcplx_constraint%names(isrfcplx))
          call PrintMsg(option)
          call InputReadDouble(input,option, &
                               srfcplx_constraint%constraint_conc(isrfcplx))
          call InputErrorMsg(input,option,'concentration',block_string)
        enddo
        call InputPopBlock(input,option)

        if (isrfcplx < reaction%surface_complexation%nkinsrfcplx) then
          option%io_buffer = &
                   'Number of surface complex constraints is less than &
                   &number of kinetic surface complexes in surface &
                   &complex constraint.'
          call PrintErrMsg(option)
        endif

        do icomp = 1, reaction%surface_complexation%nkinsrfcplx
          tempint = 0
          do jcomp = 1, reaction%surface_complexation%nkinsrfcplx
            if (StringCompare(srfcplx_constraint%names(icomp), &
                              srfcplx_constraint%names(jcomp), &
                              MAXWORDLENGTH)) tempint = tempint + 1
          enddo
          if (tempint > 1) then
            option%io_buffer = 'Duplicated surface complex in surface &
              &complex constraint: ' // &
              trim(srfcplx_constraint%names(icomp))
            call PrintErrMsg(option)
          endif
        enddo

        if (associated(constraint%surface_complexes)) then
          call ReactionSrfCplxDestroyConstraint(constraint%surface_complexes)
        endif
        constraint%surface_complexes => srfcplx_constraint

      case('IMMOBILE')

        immobile_constraint => &
          ReactionImConstraintCreate(reaction%immobile,option)

        block_string = 'CONSTRAINT, IMMOBILE'
        iimmobile = 0
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)

          if (InputCheckExit(input,option)) exit

          iimmobile = iimmobile + 1

          if (iimmobile > reaction%immobile%nimmobile) then
            option%io_buffer = &
                     'Number of immobile constraints exceeds number of ' // &
                     'immobile species in constraint: ' // &
                      trim(constraint%name)
            call PrintErrMsg(option)
          endif

          call InputReadCard(input,option, &
                             immobile_constraint%names(iimmobile))
          call InputErrorMsg(input,option,'immobile name',block_string)
          option%io_buffer = 'Constraint Immobile: ' // &
                             trim(immobile_constraint%names(iimmobile))
          call PrintMsg(option)

          ! concentration
          string = trim(input%buf)
          call InputReadWord(string,word,PETSC_TRUE,ierr)
          ! if a dataset
          if (StringCompareIgnoreCase(word,'DATASET')) then
            input%buf = trim(string)
            call InputReadWord(input,option,immobile_constraint% &
                                constraint_aux_string(iimmobile),PETSC_TRUE)
            call InputErrorMsg(input,option,'dataset name', &
                               trim(block_string) // ' concentration')
            immobile_constraint%external_dataset(iimmobile) = PETSC_TRUE
            ! set vol frac to NaN to catch bugs
            tempreal = -1.d0
            immobile_constraint%constraint_conc(iimmobile) = sqrt(tempreal)
          else
            call InputReadDouble(input,option, &
                                 immobile_constraint%constraint_conc(iimmobile))
            call InputErrorMsg(input,option,'concentration',block_string)
          endif

          ! read units if they exist
          internal_units = 'mol/m^3'
          call InputReadWord(input,option,word,PETSC_TRUE)
          string = trim(immobile_constraint%names(iimmobile)) // &
                             ' IMMOBILE CONCENTRATION UNITS'
          if (InputError(input)) then
            input%err_buf = string
            call InputDefaultMsg(input,option)
          else
            immobile_constraint%constraint_conc(iimmobile) = &
              immobile_constraint%constraint_conc(iimmobile) * &
              UnitsConvertToInternal(word,internal_units,string,option)
          endif
        enddo
        call InputPopBlock(input,option)

        if (iimmobile < reaction%immobile%nimmobile) then
          option%io_buffer = &
                   'Immobile lists in constraints must provide a &
                   &concentration for all immobile species &
                   &(listed under IMMOBILE card in CHEMISTRY), &
                   &regardless of whether or not they are present.'
          call PrintErrMsg(option)
        endif

        do icomp = 1, reaction%immobile%nimmobile
          tempint = 0
          do jcomp = 1, reaction%immobile%nimmobile
            if (StringCompare(immobile_constraint%names(icomp), &
                              immobile_constraint%names(jcomp), &
                              MAXWORDLENGTH)) tempint = tempint + 1
          enddo
          if (tempint > 1) then
            option%io_buffer = 'Duplicated immobile speies in immobile &
              &constraint: ' // &
              trim(immobile_constraint%names(icomp))
            call PrintErrMsg(option)
          endif
        enddo

        if (associated(constraint%immobile_species)) then
          call ReactionImConstraintDestroy(constraint%immobile_species)
        endif
        constraint%immobile_species => immobile_constraint

      case default
        call InputKeywordUnrecognized(input,word,'CONSTRAINT',option)
    end select

  enddo
  call InputPopBlock(input,option)

  if (.not.associated(constraint%aqueous_species)) then
    option%io_buffer = 'A CONCENTRATION block is missing in constraint "' // &
      trim(constraint%name) // '".'
    call PrintErrMsg(option)
  endif

  call PetscLogEventEnd(logging%event_tran_constraint_read, &
                        ierr);CHKERRQ(ierr)

end subroutine TranConstraintRTRead

! ************************************************************************** !

subroutine TranConstraintRTStrip(this)
  !
  ! Deallocates a constraint
  !
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  !

  implicit none

  class(tran_constraint_rt_type) :: this

  call TranConstraintBaseStrip(this)

  if (associated(this%aqueous_species)) &
    call ReactionAuxDestroySpecConstraint(this%aqueous_species)
  nullify(this%aqueous_species)
  if (associated(this%free_ion_guess)) &
    call ReactionAuxDestroyGuesConstraint(this%free_ion_guess)
  nullify(this%free_ion_guess)
  if (associated(this%minerals)) &
    call ReactionMnrlDestMnrlConstraint(this%minerals)
  nullify(this%minerals)
  if (associated(this%surface_complexes)) &
    call ReactionSrfCplxDestroyConstraint(this%surface_complexes)
  nullify(this%surface_complexes)
  if (associated(this%immobile_species)) &
    call ReactionImConstraintDestroy(this%immobile_species)
  nullify(this%immobile_species)

end subroutine TranConstraintRTStrip

! ************************************************************************** !

subroutine TranConstraintCouplerRTStrip(this)
  !
  ! Deallocate dynamic members of class
  !
  ! Author: Glenn Hammond
  ! Date: 10/04/19
  !
  implicit none

  class(tran_constraint_coupler_rt_type) :: this

  call TranConstraintCouplerBaseStrip(this)

  call RTAuxVarDestroy(this%rt_auxvar)

end subroutine TranConstraintCouplerRTStrip

end module Transport_Constraint_RT_module
