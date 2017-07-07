module PM_UFD_Biosphere_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PM_Base_class
  use Region_module
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: supported_rad_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: decay_rate                   ! [1/sec]
    PetscReal :: dcf                          ! [Sv/Bq]
    PetscReal :: kd                           ! (see note below on units)       
    PetscInt :: species_id
    PetscInt :: position_in_list
    type(supported_rad_type), pointer :: next
  end type supported_rad_type

  type, public :: unsupported_rad_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: supported_parent_name
    type(supported_rad_type), pointer :: supported_parent
    PetscReal :: decay_rate                   ! [1/sec]
    PetscReal :: dcf                          ! [Sv/Bq]
    PetscReal :: emanation_factor             ! [-]
    PetscReal :: kd                           ! (see note below on units) 
    PetscReal :: sorption_enhancement         ! [-]
    type(unsupported_rad_type), pointer :: next
  end type unsupported_rad_type
  
   ! NOTE: Kd units do not matter, but must be entered consistently in the 
   ! input file. They can be [L-water/kg-solid] or [kg-water/m3-bulk].
   ! The Kd refers to the material where the screen part of the well resides.

  type, public :: ERB_base_type
    class(ERB_base_type), pointer :: next
    character(len=MAXWORDLENGTH) :: name
    type(region_type), pointer :: region
    character(len=MAXWORDLENGTH) :: region_name
    PetscInt, pointer :: rank_list(:)
    PetscReal, pointer :: region_scaling_factor(:)
    PetscReal, pointer :: aqueous_conc_supp_rad_avg(:)          ! [mol/L]
    PetscReal, pointer :: aqueous_conc_supported_rad(:)         ! [Bq/L]
    PetscReal, pointer :: aqueous_conc_unsupported_rad(:)       ! [Bq/L]
    PetscReal, pointer :: annual_dose_supp_rad(:)               ! [Sv/yr]
    PetscReal, pointer :: annual_dose_unsupp_rad(:)             ! [Sv/yr]
    PetscReal, pointer :: annual_dose_supp_w_unsupp_rads(:)     ! [Sv/yr]
    character(len=MAXSTRINGLENGTH), pointer :: names_supp_w_unsupp_rads(:)
    PetscReal :: total_annual_dose                              ! [Sv/yr]
    PetscReal :: indv_consumption_rate                          ! [L/day]
    PetscBool :: incl_unsupported_rads
  contains
  end type ERB_base_type
  
  ! ERB_1A model assumes a pumping well (e.g., source/sink)
  type, public, extends(ERB_base_type) :: ERB_1A_type
  contains
  end type ERB_1A_type
  
  ! ERB_1B model assumes a hypothetical pumping well with a dilution factor
  type, public, extends(ERB_base_type) :: ERB_1B_type
    PetscReal :: dilution_factor
  contains
  end type ERB_1B_type

  type, public, extends(pm_base_type) :: pm_ufd_biosphere_type
    class(realization_subsurface_type), pointer :: realization
    class(ERB_base_type), pointer :: ERB_list
    type(supported_rad_type), pointer :: supported_rad_list
    type(unsupported_rad_type), pointer :: unsupported_rad_list
    PetscReal :: output_start_time
    PetscBool :: unsupp_rads_needed
  contains
    procedure, public :: PMUFDBSetRealization
    procedure, public :: Setup => PMUFDBSetup
    procedure, public :: Read => PMUFDBRead
    procedure, public :: InitializeRun => PMUFDBInitializeRun
    procedure, public :: InitializeTimestep => PMUFDBInitializeTimestep
    procedure, public :: FinalizeTimestep => PMUFDBFinalizeTimestep
    procedure, public :: Solve => PMUFDBSolve
    procedure, public :: Output => PMUFDBOutput
    procedure, public :: InputRecord => PMUFDBInputRecord
    procedure, public :: Destroy => PMUFDBDestroy
  end type pm_ufd_biosphere_type

  public :: PMUFDBCreate, &
            PMUFDB_ERB1ACreate, &
            PMUFDB_ERB1BCreate, &
            PMUFDBUnsuppRadCreate, &
            PMUFDBSupportedRadCreate
  
contains

! *************************************************************************** !

function PMUFDBCreate()
  !
  ! Creates and initializes the UFD Biosphere process model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none
  
  class(pm_ufd_biosphere_type), pointer :: PMUFDBCreate
  
  allocate(PMUFDBCreate)
  nullify(PMUFDBCreate%realization)
  nullify(PMUFDBCreate%ERB_list)
  nullify(PMUFDBCreate%unsupported_rad_list)
  nullify(PMUFDBCreate%supported_rad_list)
  PMUFDBCreate%output_start_time = 0.d0  ! [sec] default value
  PMUFDBCreate%unsupp_rads_needed = PETSC_FALSE
  PMUFDBCreate%name = 'ufd biosphere'

  call PMBaseInit(PMUFDBCreate)
  
end function PMUFDBCreate

! *************************************************************************** !

subroutine PMUFDB_ERBInit(ERB_model)
  !
  ! Initializes an ERB type model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !
  
  implicit none
  
  class(ERB_base_type) :: ERB_model
  
  ERB_model%name = ''
  ERB_model%region_name = ''
  ERB_model%indv_consumption_rate = UNINITIALIZED_DOUBLE
  ERB_model%total_annual_dose = UNINITIALIZED_DOUBLE
  ERB_model%incl_unsupported_rads = PETSC_FALSE
  nullify(ERB_model%region_scaling_factor)
  nullify(ERB_model%rank_list)
  nullify(ERB_model%aqueous_conc_supp_rad_avg)
  nullify(ERB_model%aqueous_conc_supported_rad)
  nullify(ERB_model%aqueous_conc_unsupported_rad)
  nullify(ERB_model%annual_dose_supp_rad)
  nullify(ERB_model%annual_dose_unsupp_rad)
  nullify(ERB_model%annual_dose_supp_w_unsupp_rads)
  nullify(ERB_model%names_supp_w_unsupp_rads)
  nullify(ERB_model%region)
  nullify(ERB_model%next)

end subroutine PMUFDB_ERBInit

! *************************************************************************** !

function PMUFDB_ERB1ACreate()
  !
  ! Creates and initializes an ERB_1A model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !
  
  implicit none
  
  type(ERB_1A_type), pointer :: PMUFDB_ERB1ACreate
  
  allocate(PMUFDB_ERB1ACreate)

  call PMUFDB_ERBInit(PMUFDB_ERB1ACreate)
  
end function PMUFDB_ERB1ACreate

! *************************************************************************** !

function PMUFDB_ERB1BCreate()
  !
  ! Creates and initializes an ERB_1B model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !
  
  implicit none
  
  type(ERB_1B_type), pointer :: PMUFDB_ERB1BCreate
  
  allocate(PMUFDB_ERB1BCreate)
  
  PMUFDB_ERB1BCreate%dilution_factor = UNINITIALIZED_DOUBLE
  call PMUFDB_ERBInit(PMUFDB_ERB1BCreate)
  
end function PMUFDB_ERB1BCreate

! *************************************************************************** !

function PMUFDBSupportedRadCreate()
  !
  ! Creates and initializes a supported radionuclide type.
  !
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !
  
  implicit none
  
  type(supported_rad_type), pointer :: PMUFDBSupportedRadCreate
  
  allocate(PMUFDBSupportedRadCreate)
  
  PMUFDBSupportedRadCreate%name = ''
  PMUFDBSupportedRadCreate%decay_rate = UNINITIALIZED_DOUBLE
  PMUFDBSupportedRadCreate%dcf = UNINITIALIZED_DOUBLE    
  PMUFDBSupportedRadCreate%kd = UNINITIALIZED_DOUBLE    
  PMUFDBSupportedRadCreate%species_id = 0
  PMUFDBSupportedRadCreate%position_in_list = 0
  nullify(PMUFDBSupportedRadCreate%next)
  
end function PMUFDBSupportedRadCreate

! *************************************************************************** !

function PMUFDBUnsuppRadCreate()
  !
  ! Creates and initializes an unsupported radionuclide type.
  !
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !
  
  implicit none
  
  type(unsupported_rad_type), pointer :: PMUFDBUnsuppRadCreate
  
  allocate(PMUFDBUnsuppRadCreate)
  
  PMUFDBUnsuppRadCreate%name = ''
  PMUFDBUnsuppRadCreate%supported_parent_name = ''
  PMUFDBUnsuppRadCreate%decay_rate = UNINITIALIZED_DOUBLE 
  PMUFDBUnsuppRadCreate%dcf = UNINITIALIZED_DOUBLE  
  PMUFDBUnsuppRadCreate%emanation_factor = 1.d0     
  PMUFDBUnsuppRadCreate%kd = UNINITIALIZED_DOUBLE  
  PMUFDBUnsuppRadCreate%sorption_enhancement = UNINITIALIZED_DOUBLE  
  nullify(PMUFDBUnsuppRadCreate%supported_parent)
  nullify(PMUFDBUnsuppRadCreate%next)
  
end function PMUFDBUnsuppRadCreate

! *************************************************************************** !

subroutine PMUFDBSetRealization(this,realization)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  use Realization_Subsurface_class

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMUFDBSetRealization

! *************************************************************************** !

subroutine PMUFDBRead(this,input)
  !
  ! Reads input file parameters for the UFD Biosphere process model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: double
  character(len=MAXSTRINGLENGTH) :: error_string
  type(ERB_1B_type), pointer :: new_ERB1B
  type(ERB_1A_type), pointer :: new_ERB1A
  class(ERB_base_type), pointer :: cur_ERB
  PetscBool :: added

  
  option => this%option
  input%ierr = 0
  option%io_buffer = 'pflotran card:: UFD_BIOSPHERE'
  call printMsg(option)
  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    error_string = 'UFD_BIOSPHERE'
    call StringToUpper(word)
    select case(trim(word))
    !-----------------------------------------
    !-----------------------------------------
      case('ERB_1A')
        error_string = trim(error_string) // ',ERB_1A'
        allocate(new_ERB1A)
        new_ERB1A => PMUFDB_ERB1ACreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        new_ERB1A%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_ERB1A%name)
        call PMUFDBReadERBmodel(this,input,option,new_ERB1A,error_string)
        ! add new ERB_1A model to ERB_list
        added = PETSC_FALSE
        if (.not.associated(this%ERB_list)) then
          this%ERB_list => new_ERB1A
        else
          cur_ERB => this%ERB_list
          do
            if (.not.associated(cur_ERB)) exit
            if (.not.associated(cur_ERB%next)) then
              cur_ERB%next => new_ERB1A
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_ERB => cur_ERB%next
          enddo
        endif
        nullify(new_ERB1A)
    !-----------------------------------------
    !-----------------------------------------
      case('ERB_1B')
        error_string = trim(error_string) // ',ERB_1B'
        allocate(new_ERB1B)
        new_ERB1B => PMUFDB_ERB1BCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        new_ERB1B%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_ERB1B%name)
        call PMUFDBReadERBmodel(this,input,option,new_ERB1B,error_string)
        ! add new ERB_1B model to ERB_list
        added = PETSC_FALSE
        if (.not.associated(this%ERB_list)) then
          this%ERB_list => new_ERB1B
        else
          cur_ERB => this%ERB_list
          do
            if (.not.associated(cur_ERB)) exit
            if (.not.associated(cur_ERB%next)) then
              cur_ERB%next => new_ERB1B
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_ERB => cur_ERB%next
          enddo
        endif
        nullify(new_ERB1B)
    !-----------------------------------------
    !-----------------------------------------
      case('SUPPORTED_RADIONUCLIDES')
        error_string = trim(error_string) // ',SUPPORTED_RADIONUCLIDES'
        call PMUFDBReadSupportedRad(this,input,option,error_string)  
    !-----------------------------------------
    !-----------------------------------------
      case('UNSUPPORTED_RADIONUCLIDES')
        error_string = trim(error_string) // ',UNSUPPORTED_RADIONUCLIDES'
        call PMUFDBReadUnsuppRad(this,input,option,error_string)      
    !-----------------------------------------
    !-----------------------------------------
      case('OUTPUT_START_TIME')
        call InputReadDouble(input,option,double)
        call InputErrorMsg(input,option,'OUTPUT_START_TIME',error_string)
        call InputReadAndConvertUnits(input,double,'sec',trim(error_string) &
                                      // ',OUTPUT_START_TIME units',option)
        this%output_start_time = double
    !-----------------------------------------
    !-----------------------------------------
      case default
        call InputKeywordUnrecognized(word,'UFD_BIOSPHERE',option)
    !-----------------------------------------
    end select  
  enddo
  
  ! error messages
  if (.not.associated(this%ERB_list)) then
    option%io_buffer = 'At least one ERB model must be specified &
                       &in the ' // trim(error_string) // ' block, with the &
                       &keyword ERB_1A or ERB_1B.'
    call printErrMsg(option)
  endif
  if (.not.associated(this%supported_rad_list)) then
    option%io_buffer = 'At least one supported radionuclide must be specified &
                       &in the ' // trim(error_string) // ' block, with the &
                       &keyword SUPPORTED_RADIONUCLIDES.'
    call printErrMsg(option)
  endif
  if (this%unsupp_rads_needed .and. &
     (.not.associated(this%unsupported_rad_list))) then
    option%io_buffer = 'At least one ERB model indicates that unsupported &
                       &radionuclides should be included, but no unsupported &
                       &radionuclides were specified. You must specify all &
                       &unsupported radionuclides using the keyword &
                       &UNSUPPORTED_RADIONUCLIDES block.'
    call printErrMsg(option)
  endif
  
  
end subroutine PMUFDBRead

! *************************************************************************** !

subroutine PMUFDBReadERBmodel(this,input,option,ERB_model,error_string)
  !
  ! Reads input file parameters for the UFD Biosphere process model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  use Input_Aux_module
  use Option_module
  use String_module

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(ERB_base_type) :: ERB_model
  character(len=MAXSTRINGLENGTH) :: error_string
  
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: double
  PetscInt :: num_errors
  
  num_errors = 0
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(trim(word))
    !-----------------------------------
      case('REGION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'region assignment',error_string)
        ERB_model%region_name = trim(word)
    !-----------------------------------
      case('DILUTION_FACTOR')
        select type(ERB_model)
          type is(ERB_1A_type)
            option%io_buffer = 'ERROR: DILUTION_FACTOR cannot be specified &
                               &for ' // trim(error_string)
            call printMsg(option)
            num_errors = num_errors + 1
          type is(ERB_1B_type)
            call InputReadDouble(input,option,ERB_model%dilution_factor)
            call InputErrorMsg(input,option,'DILUTION_FACTOR',error_string)
        end select
    !-----------------------------------
      case('INDIVIDUAL_CONSUMPTION_RATE')
          call InputReadDouble(input,option,double)
          call InputErrorMsg(input,option,'INDIVIDUAL_CONSUMPTION_RATE', &
                             error_string)
          call InputReadAndConvertUnits(input,double,'L/day', &
               trim(error_string) // ',INDIVIDUAL_CONSUMPTION_RATE units', &
               option)
          ERB_model%indv_consumption_rate = double
    !-----------------------------------
      case('INCLUDE_UNSUPPORTED_RADS')
        ERB_model%incl_unsupported_rads = PETSC_TRUE
        this%unsupp_rads_needed = PETSC_TRUE
    !-----------------------------------    
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    !----------------------------------- 
    end select
  enddo

  ! error messages
  if (ERB_model%region_name == '') then
    option%io_buffer = 'ERROR: REGION must be specified in the ' // &
                       trim(error_string) // ' block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(ERB_model%indv_consumption_rate)) then
    option%io_buffer = 'ERROR: INDIVIDUAL_CONSUMPTION_RATE must be specified &
                       &in the ' // trim(error_string) // ' block.'
    call printMsg(option)
    num_errors = num_errors + 1
  endif
  
  select type(ERB_model)
    type is(ERB_1B_type)
      if (Uninitialized(ERB_model%dilution_factor)) then
        option%io_buffer = 'ERROR: DILUTION_FACTOR must be specified in &
                           &the ' // trim(error_string) // ' block.'
        call printMsg(option)
        num_errors = num_errors + 1
      endif
  end select
  
  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the ' //trim(error_string) // ' block. See above.'
    call printErrMsg(option)
  endif
  
end subroutine PMUFDBReadERBmodel

! *************************************************************************** !

subroutine PMUFDBReadSupportedRad(this,input,option,error_string)
  !
  ! Reads input file parameters for the UFD Biosphere process model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  use Input_Aux_module
  use Option_module
  use String_module

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  type(supported_rad_type), pointer :: new_supp_rad
  type(supported_rad_type), pointer :: cur_supp_rad
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: double
  PetscInt :: num_errors
  PetscBool :: added
  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    num_errors = 0
    call StringToUpper(word)
    select case(trim(word))
    !-----------------------------------
      case('RADIONUCLIDE')
        error_string = trim(error_string) // ',RADIONUCLIDE'
        allocate(new_supp_rad)
        new_supp_rad => PMUFDBSupportedRadCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'radionuclide name',error_string)
        new_supp_rad%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_supp_rad%name)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
          !-----------------------------------
            case('ELEMENT_KD')
              call InputReadDouble(input,option,new_supp_rad%kd)
              call InputErrorMsg(input,option,'ELEMENT_KD',error_string)
          !-----------------------------------
            case('INGESTION_DOSE_COEF')
              call InputReadDouble(input,option,new_supp_rad%dcf)
              call InputErrorMsg(input,option,'INGESTION_DOSE_COEF', &
                                 error_string)
          !-----------------------------------
            case('DECAY_RATE')
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'DECAY_RATE',error_string)
              call InputReadAndConvertUnits(input,double,'1/sec', &
                     trim(error_string) // ',DECAY_RATE units',option)
              new_supp_rad%decay_rate = double
          !-----------------------------------  
            case default
              call InputKeywordUnrecognized(word,error_string,option)
          !----------------------------------- 
          end select
        enddo
        ! error messages
        if (Uninitialized(new_supp_rad%kd)) then
          option%io_buffer = 'ERROR: ELEMENT_KD must be specified &
                             &in the ' // trim(error_string) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_supp_rad%dcf)) then
          option%io_buffer = 'ERROR: INGESTION_DOSE_COEF must be specified &
                             &in the ' // trim(error_string) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_supp_rad%decay_rate)) then
          option%io_buffer = 'ERROR: DECAY_RATE must be specified &
                             &in the ' // trim(error_string) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (num_errors > 0) then
          write(option%io_buffer,*) num_errors
          option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the ' //trim(error_string) // ' block. See above.'
          call printErrMsg(option)
        endif
        ! add new supported radionuclide to list
        added = PETSC_FALSE
        if (.not.associated(this%supported_rad_list)) then
          this%supported_rad_list => new_supp_rad
        else
          cur_supp_rad => this%supported_rad_list
          do
            if (.not.associated(cur_supp_rad)) exit
            if (.not.associated(cur_supp_rad%next)) then
              cur_supp_rad%next => new_supp_rad
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_supp_rad => cur_supp_rad%next
          enddo
        endif
        nullify(new_supp_rad)
    !-----------------------------------
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    !-----------------------------------
    end select
  enddo
  
end subroutine PMUFDBReadSupportedRad


! *************************************************************************** !

subroutine PMUFDBReadUnsuppRad(this,input,option,error_string)
  !
  ! Reads input file parameters for the UFD Biosphere process model.
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  use Input_Aux_module
  use Option_module
  use String_module

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  type(unsupported_rad_type), pointer :: new_unsupp_rad
  type(unsupported_rad_type), pointer :: cur_unsupp_rad
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: double
  PetscInt :: num_errors
  PetscBool :: added
  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    num_errors = 0
    call StringToUpper(word)
    select case(trim(word))
    !-----------------------------------
      case('RADIONUCLIDE')
        error_string = trim(error_string) // ',RADIONUCLIDE'
        allocate(new_unsupp_rad)
        new_unsupp_rad => PMUFDBUnsuppRadCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'radionuclide name',error_string)
        new_unsupp_rad%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_unsupp_rad%name)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
          !-----------------------------------
            case('ELEMENT_KD')
              call InputReadDouble(input,option,new_unsupp_rad%kd)
              call InputErrorMsg(input,option,'ELEMENT_KD',error_string)
          !-----------------------------------
            case('SUPPORTED_PARENT')
              call InputReadWord(input,option, &
                              new_unsupp_rad%supported_parent_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'SUPPORTED_PARENT',error_string)
          !-----------------------------------
            case('INGESTION_DOSE_COEF')
              call InputReadDouble(input,option,new_unsupp_rad%dcf)
              call InputErrorMsg(input,option,'INGESTION_DOSE_COEF', &
                                 error_string)
          !-----------------------------------
            case('DECAY_RATE')
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'DECAY_RATE',error_string)
              call InputReadAndConvertUnits(input,double,'1/sec', &
                     trim(error_string) // ',DECAY_RATE units',option)
              new_unsupp_rad%decay_rate = double
          !-----------------------------------
            case('EMANATION_FACTOR')
              call InputReadDouble(input,option,new_unsupp_rad%emanation_factor)
              call InputErrorMsg(input,option,'EMANATION_FACTOR',error_string)
          !-----------------------------------  
            case default
              call InputKeywordUnrecognized(word,error_string,option)
          !----------------------------------- 
          end select
        enddo
        ! error messages
        if (Uninitialized(new_unsupp_rad%kd)) then
          option%io_buffer = 'ERROR: ELEMENT_KD must be specified in the ' // &
                             trim(error_string) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_unsupp_rad%dcf)) then
          option%io_buffer = 'ERROR: INGESTION_DOSE_COEF must be specified in &
                             &the ' // trim(error_string) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (Uninitialized(new_unsupp_rad%decay_rate)) then
          option%io_buffer = 'ERROR: DECAY_RATE must be specified &
                             &in the ' // trim(error_string) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (new_unsupp_rad%supported_parent_name == '') then
          option%io_buffer = 'ERROR: SUPPORTED_PARENT must be specified in &
                             &the ' // trim(error_string) // ' block.'
          call printMsg(option)
          num_errors = num_errors + 1
        endif
        if (num_errors > 0) then
          write(option%io_buffer,*) num_errors
          option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the ' //trim(error_string) // ' block. See above.'
          call printErrMsg(option)
        endif
        ! add new unsupported radionuclide to list
        added = PETSC_FALSE
        if (.not.associated(this%unsupported_rad_list)) then
          this%unsupported_rad_list => new_unsupp_rad
        else
          cur_unsupp_rad => this%unsupported_rad_list
          do
            if (.not.associated(cur_unsupp_rad)) exit
            if (.not.associated(cur_unsupp_rad%next)) then
              cur_unsupp_rad%next => new_unsupp_rad
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_unsupp_rad => cur_unsupp_rad%next
          enddo
        endif
        nullify(new_unsupp_rad)
    !-----------------------------------
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    !-----------------------------------
    end select
  enddo
  
end subroutine PMUFDBReadUnsuppRad

! *************************************************************************** !

subroutine PMUFDBAssociateRegion(this,region_list)
  ! 
  ! Associates the ERB model to its assigned region via the REGION keyword.
  ! And calculates the scaling factor by volume.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !

  use Region_module
  use Option_module
  use String_module
  use Material_Aux_class
  use Grid_module
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  type(region_list_type), pointer :: region_list
  
  type(region_type), pointer :: cur_region
  class(ERB_base_type), pointer :: cur_ERB
  type(option_type), pointer :: option
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type), pointer :: grid
  PetscInt :: ghosted_id
  PetscReal :: total_volume_local, total_volume_global
  PetscInt :: k 
  PetscErrorCode :: ierr
  
  option => this%option
  material_auxvars => this%realization%patch%aux%Material%auxvars
  grid => this%realization%patch%grid
  
  cur_ERB => this%ERB_list
  do
    if (.not.associated(cur_ERB)) exit
      cur_region => region_list%first     
      do
        if (.not.associated(cur_region)) exit
        if (StringCompare(cur_region%name, &
                          cur_ERB%region_name)) then
          cur_ERB%region => cur_region
          ! calculate scaling factor by cell volumes in region
          allocate(cur_ERB%region_scaling_factor(cur_ERB%region%num_cells))
          total_volume_global = 0.d0
          total_volume_local = 0.d0
          do k = 1,cur_ERB%region%num_cells
            ghosted_id = grid%nL2G(cur_ERB%region%cell_ids(k))
            cur_ERB%region_scaling_factor(k) = &
                                   material_auxvars(ghosted_id)%volume ! [m^3]
            total_volume_local = total_volume_local &
                                 + material_auxvars(ghosted_id)%volume ! [m^3]
          enddo
          call MPI_Allreduce(total_volume_local,total_volume_global, &
                             ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION,MPI_SUM, &
                             option%mycomm,ierr)
          cur_ERB%region_scaling_factor = cur_ERB%region_scaling_factor / &
                                          total_volume_global
          exit
        endif
        cur_region => cur_region%next
      enddo      
      if (.not.associated(cur_ERB%region)) then
        option%io_buffer = 'ERB model (' // trim(cur_ERB%name) // ') REGION ' &
                           // trim(cur_ERB%region_name) // ' not found among &
                           &defined regions.'
        call printErrMsg(option)
      endif
    cur_ERB => cur_ERB%next
  enddo
  
end subroutine PMUFDBAssociateRegion

! *************************************************************************** !

subroutine PMUFDBSetup(this)
  !
  ! Sets up the process model with external information.
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
  PetscMPIInt :: newcomm_size
  PetscInt, pointer :: ranks(:)
  PetscInt :: i, j
  PetscBool :: local
  PetscErrorCode :: ierr
  class(ERB_base_type), pointer :: cur_ERB
  class(ERB_base_type), pointer :: prev_ERB
  class(ERB_base_type), pointer :: next_ERB
  
  call PMUFDBAssociateRegion(this,this%realization%patch%region_list)
  call PMUFDBSupportedRadCheckRT(this)
  call PMUFDBAscUnsuppRadWithSuppRad(this)
  call PMUFDBCheckSrcSinkCouplers(this)
  call PMUFDBAllocateERBarrays(this)
  
  nullify(prev_ERB)
  allocate(ranks(this%option%mycommsize))
  cur_ERB => this%ERB_list
  do
    if (.not.associated(cur_ERB)) exit
    if (associated(cur_ERB%region)) then
      if (cur_ERB%region%num_cells > 0) then
          local = PETSC_TRUE
      endif
    endif
    ranks(:) = 0
    newcomm_size = 0
    if (local) ranks(this%option%myrank+1) = 1
    if (.not.local) ranks(this%option%myrank+1) = 0
    call MPI_Allreduce(MPI_IN_PLACE,ranks,this%option%mycommsize,MPI_INTEGER, &
                       MPI_SUM,this%option%mycomm,ierr)
    newcomm_size = sum(ranks)
    allocate(cur_ERB%rank_list(newcomm_size))
    j = 0
    do i = 1,this%option%mycommsize
      if (ranks(i) == 1) then
        j = j + 1
        cur_ERB%rank_list(j) = (i - 1)  ! (world ranks)
      endif
    enddo  
    if (local) then
      prev_ERB => cur_ERB
      cur_ERB => cur_ERB%next
    else 
      ! remove ERB model because it is not local
      next_ERB => cur_ERB%next
      if (associated(prev_ERB)) then
        prev_ERB%next => next_ERB
      else
        this%ERB_list => next_ERB
      endif
      call PMUFDBDestroyERB(cur_ERB)
      cur_ERB => next_ERB
    endif
    !cur_ERB => cur_ERB%next
  enddo
  deallocate(ranks)
  
end subroutine PMUFDBSetup

! *************************************************************************** !

subroutine PMUFDBSupportedRadCheckRT(this)
  ! 
  ! Checks whether the supported radionuclides are also defined as a
  ! primary or a secondary species in the CHEMISTRY block.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !

  use Option_module
  use String_module
  use Reaction_Aux_module
  
  implicit none
  
  type(pm_ufd_biosphere_type) :: this
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH), pointer :: pri_names(:)
  character(len=MAXWORDLENGTH), pointer :: sec_names(:)
  type(supported_rad_type), pointer :: cur_supp_rad
  character(len=MAXWORDLENGTH) :: rad_name
  PetscBool :: found
  PetscInt :: i
  
  option => this%option
  
  if (associated(this%realization%reaction)) then
    allocate(pri_names(GetPrimarySpeciesCount(this%realization%reaction)))
    pri_names => GetPrimarySpeciesNames(this%realization%reaction)
    allocate(sec_names(GetSecondarySpeciesCount(this%realization%reaction)))
    sec_names => GetSecondarySpeciesNames(this%realization%reaction)
  else
    option%io_buffer = 'The UFD_BIOSPHERE process model requires reactive &
                       &transport.'
    call printErrMsg(option)
  endif
  
  cur_supp_rad => this%supported_rad_list
  do
    if (.not.associated(cur_supp_rad)) exit
    rad_name = cur_supp_rad%name
    found = PETSC_FALSE
    do i = 1,len(pri_names)
      if (adjustl(trim(rad_name)) == adjustl(trim(pri_names(i)))) then
        found = PETSC_TRUE
        cur_supp_rad%species_id = GetPrimarySpeciesIDFromName(rad_name, &
                                              this%realization%reaction,option)
      endif
      if (found) exit
    enddo
    if (.not.found) then
      do i = 1,len(sec_names)
        if (adjustl(trim(rad_name)) == adjustl(trim(sec_names(i)))) then
          found = PETSC_TRUE
          cur_supp_rad%species_id = GetSecondarySpeciesIDFromName(rad_name, &
                                              this%realization%reaction,option)
        endif
        if (found) exit
      enddo
    endif
    if (.not.found) then
      option%io_buffer = 'SUPPORTED_RADIONUCLIDE ' // trim(rad_name) // &
        ' must be included as a primary or secondary species under the &
        &CHEMISTRY block.'
      call printErrMsg(option)
    endif
    
    cur_supp_rad => cur_supp_rad%next
  enddo
  
  deallocate(pri_names)
  deallocate(sec_names)
  
end subroutine PMUFDBSupportedRadCheckRT

! *************************************************************************** !

subroutine PMUFDBAscUnsuppRadWithSuppRad(this)
  ! 
  ! Associates the unsupported radionuclide with its support parent 
  ! radionuclide.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/22/2017
  !

  use Option_module
  use String_module
  use Material_Aux_class
  use Grid_module
  
  implicit none
  
  type(pm_ufd_biosphere_type) :: this
  
  type(supported_rad_type), pointer :: cur_supp_rad
  type(unsupported_rad_type), pointer :: cur_unsupp_rad
  class(ERB_base_type), pointer :: cur_ERB
  class(material_auxvar_type), pointer :: material_auxvar(:)
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  PetscInt :: ghosted_id
  PetscBool :: no_density
  
  option => this%option
  grid => this%realization%patch%grid
  material_auxvar => this%realization%patch%aux%Material%auxvars
  
  cur_unsupp_rad => this%unsupported_rad_list
  do
    if (.not.associated(cur_unsupp_rad)) exit
    cur_supp_rad => this%supported_rad_list     
    do
      if (.not.associated(cur_supp_rad)) exit
      if (StringCompare(cur_supp_rad%name, &
                        cur_unsupp_rad%supported_parent_name)) then
        cur_unsupp_rad%supported_parent => cur_supp_rad
        exit
      endif
      cur_supp_rad => cur_supp_rad%next
    enddo      
    if (.not.associated(cur_unsupp_rad%supported_parent)) then
      option%io_buffer = 'UNSUPPORTED_RADIONUCLIDE ' // &
        trim(cur_unsupp_rad%name) // "'s SUPPORTED_PARENT " &
        // trim(cur_unsupp_rad%supported_parent_name) // ' not found among &
        &defined supported radionuclides.'
      call printErrMsg(option)
    endif
    cur_unsupp_rad => cur_unsupp_rad%next
  enddo
  
  ! check that rock density has been specified because it is needed 
  ! for the Rfi/Rfu calculation later for all models that include
  ! unsupported radionuclides
  if (associated(this%unsupported_rad_list)) then
    cur_ERB => this%ERB_list
    do
      if (.not.associated(cur_ERB)) exit
      if (cur_ERB%region%num_cells > 0) then
        ghosted_id = grid%nL2G(cur_ERB%region%cell_ids(1))
        no_density = Uninitialized(material_auxvar(ghosted_id)% &
                                   soil_particle_density)
        if (cur_ERB%incl_unsupported_rads .and. no_density) then
          option%io_buffer = 'ROCK_DENSITY is not specified in the &
              &MATERIAL_PROPERTY that ERB Model ' // trim(cur_ERB%name) // &
              ' resides in.'
          call printErrMsgByRank(option)
        endif
      endif
    cur_ERB => cur_ERB%next
    enddo
  endif
  
end subroutine PMUFDBAscUnsuppRadWithSuppRad

! *************************************************************************** !

subroutine PMUFDBCheckSrcSinkCouplers(this)
  ! 
  ! Checks whether each ERB_1A model has a region that matches a defined
  ! src_sink coupler region, and that ERB_1B models do not have regions
  ! with src_sink couplers.
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/20/2017
  !
  
  use Coupler_module
  use Option_module
  use String_module
  
  implicit none
  
  type(pm_ufd_biosphere_type) :: this
  
  class(ERB_base_type), pointer :: cur_ERB
  type(coupler_type), pointer :: cur_ss
  type(option_type), pointer :: option
  PetscBool :: ERB1A_matched
  
  option => this%option
  
  cur_ERB => this%ERB_list
  do
    if (.not.associated(cur_ERB)) exit
    ERB1A_matched = PETSC_FALSE
    cur_ss => this%realization%patch%source_sink_list%first
    do
      if (.not.associated(cur_ss)) exit
      select type(cur_ERB)
      !------------------------------------
        type is(ERB_1A_type)
          if (StringCompare(cur_ss%region_name,cur_ERB%region_name)) then
            ERB1A_matched = PETSC_TRUE
            exit
          endif
      !------------------------------------
        type is(ERB_1B_type)
          if (StringCompare(cur_ss%region_name,cur_ERB%region_name)) then
            option%io_buffer = 'ERB_1B model "' // trim(cur_ERB%name) // &
              '" REGION "' // trim(cur_ERB%region_name) // '" cannot be &
              &associated with a SOURCE_SINK. Did you intend an ERB_1A model?' 
            call printErrMsg(option)
          endif
      !------------------------------------
      end select   
      cur_ss => cur_ss%next
    enddo
    select type(cur_ERB)
      type is(ERB_1A_type)
        if (.not.ERB1A_matched) then
          option%io_buffer = 'ERB_1A model "' // trim(cur_ERB%name) // &
                '" REGION "' // trim(cur_ERB%region_name) // '" is not &
                &associated with a SOURCE_SINK. ERB_1A model types must &
                &indicate a REGION that is associated with a SOURCE_SINK. &
                &Did you intend an ERB_1B model type?' 
          call printErrMsg(option)
        endif
    end select
    cur_ERB => cur_ERB%next
  enddo
  
end subroutine PMUFDBCheckSrcSinkCouplers

! *************************************************************************** !

subroutine PMUFDBAllocateERBarrays(this)
  ! 
  ! Allocates the concentraton and dose arrays in the ERB models. Also 
  ! allocates the radionuclide data arrays.
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/06/2017
  !

  use Option_module
  
  implicit none
  
  type(pm_ufd_biosphere_type) :: this
  
  class(ERB_base_type), pointer :: cur_ERB
  type(supported_rad_type), pointer :: cur_supp_rad
  type(unsupported_rad_type), pointer :: cur_unsupp_rad
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: temp_string
  PetscInt :: num_supp_rads
  PetscInt :: num_unsupp_rads
  PetscInt :: position
  PetscInt :: k
  
  option => this%option
  
  !----- Count the supported radionuclides --------------------------------
  num_supp_rads = 0
  cur_supp_rad => this%supported_rad_list
  do
    if (.not.associated(cur_supp_rad)) exit
    num_supp_rads = num_supp_rads + 1
    cur_supp_rad%position_in_list = num_supp_rads
    cur_supp_rad => cur_supp_rad%next
  enddo
  
  !----- Count the unsupported radionuclides ------------------------------
  num_unsupp_rads = 0
  cur_unsupp_rad => this%unsupported_rad_list
  do
    if (.not.associated(cur_unsupp_rad)) exit
    num_unsupp_rads = num_unsupp_rads + 1
    cur_unsupp_rad => cur_unsupp_rad%next
  enddo
  
  !----- Allocate the arrays ----------------------------------------------
  cur_ERB => this%ERB_list
  do
    if (.not.associated(cur_ERB)) exit
    allocate(cur_ERB%aqueous_conc_supp_rad_avg(num_supp_rads))
    allocate(cur_ERB%aqueous_conc_supported_rad(num_supp_rads))
    allocate(cur_ERB%aqueous_conc_unsupported_rad(num_unsupp_rads))
    allocate(cur_ERB%annual_dose_supp_rad(num_supp_rads))
    allocate(cur_ERB%annual_dose_unsupp_rad(num_unsupp_rads))
    allocate(cur_ERB%annual_dose_supp_w_unsupp_rads(num_supp_rads))
    allocate(cur_ERB%names_supp_w_unsupp_rads(num_supp_rads))
    cur_ERB%aqueous_conc_supp_rad_avg(:) = 0.d0
    cur_ERB%aqueous_conc_supported_rad(:) = 0.d0
    cur_ERB%aqueous_conc_unsupported_rad(:) = 0.d0
    cur_ERB%annual_dose_supp_rad(:) = 0.d0
    cur_ERB%annual_dose_unsupp_rad(:) = 0.d0
    cur_ERB%annual_dose_supp_w_unsupp_rads(:) = 0.d0
    cur_ERB%names_supp_w_unsupp_rads(:) = ''
    cur_ERB => cur_ERB%next
  enddo
  
  !-----Build-the-array-of-names-for-output--------------------------------
  cur_ERB => this%ERB_list
  do
    if (.not.associated(cur_ERB)) exit
    k = 0
    cur_supp_rad => this%supported_rad_list
    do
      if (.not.associated(cur_supp_rad)) exit
      k = k + 1
      cur_ERB%names_supp_w_unsupp_rads(k) = trim(cur_supp_rad%name)
      cur_supp_rad => cur_supp_rad%next
    enddo
    cur_unsupp_rad => this%unsupported_rad_list
    if (cur_ERB%incl_unsupported_rads) then
      do
        if (.not.associated(cur_unsupp_rad)) exit
        position = cur_unsupp_rad%supported_parent%position_in_list
        temp_string = trim(cur_ERB%names_supp_w_unsupp_rads(position)) // '+' &
                      // trim(cur_unsupp_rad%name) // '*'
        cur_ERB%names_supp_w_unsupp_rads(position) = temp_string
        cur_unsupp_rad => cur_unsupp_rad%next
      enddo
    endif
    cur_ERB => cur_ERB%next
  enddo
  
end subroutine PMUFDBAllocateERBarrays

! ************************************************************************** !

subroutine PMUFDBInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none

  class(pm_ufd_biosphere_type) :: this
  
  PetscErrorCode :: ierr
  
  ierr = 0
  
  call PMUFDBOutputHeader(this)
  call PMUFDBSolve(this,0.d0,ierr)
  
end subroutine PMUFDBInitializeRun

! *************************************************************************** !

subroutine PMUFDBInitializeTimestep(this)
  ! 
  ! Initializes the process model to take a time step in the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  implicit none

  class(pm_ufd_biosphere_type) :: this
  
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," UFD BIOSPHERE MODEL ",57("="))')
  endif
  
  if (this%option%time >= this%output_start_time) then
    call PMUFDBOutput(this)
  endif

  
end subroutine PMUFDBInitializeTimestep

! *************************************************************************** !

 subroutine PMUFDBSolve(this,time,ierr)
  ! 
  ! 
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  use Utility_module
  use Grid_module
  use Material_Aux_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  class(ERB_base_type), pointer :: cur_ERB
  type(supported_rad_type), pointer :: cur_supp_rad
  type(unsupported_rad_type), pointer :: cur_unsupp_rad
  PetscReal :: Rfi, Rfu
  PetscReal :: dry_bulk_density
  PetscReal :: water_content
  PetscReal :: den_sat_ratio
  PetscReal :: local_conc, global_conc
  PetscInt :: species_id
  PetscInt :: ghosted_id
  PetscInt :: position
  PetscInt :: i, k
  type(grid_type), pointer :: grid
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscReal, parameter :: avagadro = 6.0221409d23
  PetscInt, parameter :: iphase = 1  ! LIQUID_PHASE
  
  rt_auxvars => this%realization%patch%aux%RT%auxvars
  global_auxvars => this%realization%patch%aux%Global%auxvars
  material_auxvars => this%realization%patch%aux%Material%auxvars
  grid => this%realization%patch%grid
  
  ierr = 0
  
  cur_ERB => this%ERB_list
  do
    if (.not.associated(cur_ERB)) exit
    cur_supp_rad => this%supported_rad_list
    k = 0
    do
      if (.not.associated(cur_supp_rad)) exit
      k = k + 1
      species_id = cur_supp_rad%species_id
      !-----Calculate-average-molar-concentration-at-(potential)-well------   
      local_conc = 0.d0
      do i = 1,cur_ERB%region%num_cells
        ghosted_id = grid%nL2G(cur_ERB%region%cell_ids(i))
        local_conc = local_conc + &                              ! [mol/L]
            (rt_auxvars(ghosted_id)%total(species_id,iphase) * & ! [mol/L]
             cur_ERB%region_scaling_factor(i))                   ! [-]
      enddo
      call CalcParallelSum(this%option,cur_ERB%rank_list, &
                           local_conc,global_conc)
      cur_ERB%aqueous_conc_supp_rad_avg(k) = global_conc
      !-----Calculate-aqueous-concentration-(DF=1)-------------------------
      select type(cur_ERB)
        type is(ERB_1A_type)
          cur_ERB%aqueous_conc_supported_rad(k) = &           ! [Bq/L]
              cur_ERB%aqueous_conc_supp_rad_avg(k) * &        ! [mol/L]
              avagadro * &                                    ! [1/mol]
              cur_supp_rad%decay_rate                         ! [1/sec]
      !-----Calculate-aqueous-concentration-based-on-dilution-factor-------
        type is(ERB_1B_type)          
          cur_ERB%aqueous_conc_supported_rad(k) = &           ! [Bq/L]
              cur_ERB%aqueous_conc_supp_rad_avg(k) * &        ! [mol/L]
              avagadro * &                                    ! [1/mol]
              cur_supp_rad%decay_rate / &                     ! [1/sec] 
              cur_ERB%dilution_factor                         ! [-]       
      !--------------------------------------------------------------------
      end select
      !-----Calculate-dose:-supported-radionuclides------------------------
      cur_ERB%annual_dose_supp_rad(k) = &                     ! [Sv/yr]
          cur_ERB%aqueous_conc_supported_rad(k) * &           ! [Bq/L]
          cur_ERB%indv_consumption_rate * &                   ! [L/day]
          (365.d0/1.d0) * &                                   ! [day/L]
          cur_supp_rad%dcf                                    ! [Sv/Bq]
      !-----Initialize-dose-from-supp'd-+-unsupp'd-rads--------------------
      cur_ERB%annual_dose_supp_w_unsupp_rads(k) = &
          cur_ERB%annual_dose_supp_rad(k)  
      !--------------------------------------------------------------------
      cur_supp_rad => cur_supp_rad%next
    enddo
    
    !-----Unsupported-radionuclides----------------------------------------
    if (cur_ERB%incl_unsupported_rads) then
      dry_bulk_density = 0.d0
      do i = 1,cur_ERB%region%num_cells
        ghosted_id = grid%nL2G(cur_ERB%region%cell_ids(i))
        dry_bulk_density = dry_bulk_density + &   ! [kg/m3]
         ((material_auxvars(ghosted_id)%soil_particle_density * &
          (1.d0-material_auxvars(ghosted_id)%porosity)) * &
          cur_ERB%region_scaling_factor(i))
      enddo  ! get average dry_bulk_density over region:
      call CalcParallelSum(this%option,cur_ERB%rank_list, &
                           dry_bulk_density,dry_bulk_density)
      dry_bulk_density = dry_bulk_density/1.d3   ! [kg/L]
      water_content = 0.d0
      do i = 1,cur_ERB%region%num_cells
        ghosted_id = grid%nL2G(cur_ERB%region%cell_ids(i))
        water_content = water_content + &
         ((material_auxvars(ghosted_id)%porosity * &        
          global_auxvars(ghosted_id)%sat(LIQUID_PHASE)) * &
          cur_ERB%region_scaling_factor(i))
      enddo  ! get average water_content over region:
      call CalcParallelSum(this%option,cur_ERB%rank_list, &
                           water_content,water_content) 
      den_sat_ratio = (dry_bulk_density/water_content)
      cur_unsupp_rad => this%unsupported_rad_list
      k = 0
      do
        if (.not.associated(cur_unsupp_rad)) exit
        k = k + 1
    !-----Calculate-aqueous-concentration-of-unsupported-rads--------------
        Rfi = 1.d0 + cur_unsupp_rad%supported_parent%kd*(den_sat_ratio)
        Rfu = 1.d0 + cur_unsupp_rad%kd*(den_sat_ratio) 
        cur_unsupp_rad%sorption_enhancement = (Rfi/Rfu)
        position = cur_unsupp_rad%supported_parent%position_in_list
        cur_ERB%aqueous_conc_unsupported_rad(k) = &    
            cur_ERB%aqueous_conc_supported_rad(position) * & 
            (Rfi/Rfu) * cur_unsupp_rad%emanation_factor
    !-----Calculate-dose:-unsupported-radionuclides------------------------
        cur_ERB%annual_dose_unsupp_rad(k) = &                   ! [Sv/yr]
          cur_ERB%aqueous_conc_unsupported_rad(k) * &           ! [Bq/L]
          cur_ERB%indv_consumption_rate * &                     ! [L/day]
          (365.d0/1.d0) * &                                     ! [day/L]
          cur_unsupp_rad%dcf                                    ! [Sv/Bq]
    !-----Calculate-dose-from-supp'd-rads-with-their-unsupp'd-desc's-------
        cur_ERB%annual_dose_supp_w_unsupp_rads(position) = &    ! [Sv/yr]
          cur_ERB%annual_dose_supp_w_unsupp_rads(position) + &  ! [Sv/yr]
          cur_ERB%annual_dose_unsupp_rad(k)                     ! [Sv/yr]
        cur_unsupp_rad => cur_unsupp_rad%next
      enddo
    endif     
    !-----Calculate-total-annual-dose--------------------------------------
    cur_ERB%total_annual_dose = &                               ! [Sv/yr]
        sum(cur_ERB%annual_dose_supp_w_unsupp_rads)             ! [Sv/yr]
    !----------------------------------------------------------------------
    cur_ERB => cur_ERB%next
  enddo

end subroutine PMUFDBSolve

! ************************************************************************** !

subroutine PMUFDBFinalizeTimestep(this)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
end subroutine PMUFDBFinalizeTimestep

! *************************************************************************** !

subroutine PMUFDBOutput(this)
  ! 
  ! Sets up output for the process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  
  use Option_module
  use Output_Aux_module
  
  implicit none

  class(pm_ufd_biosphere_type) :: this
  
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  class(ERB_base_type), pointer :: cur_ERB
  type(supported_rad_type), pointer :: cur_supp_rad
  type(unsupported_rad_type), pointer :: cur_unsupp_rad
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: k
  PetscReal, parameter :: avagadro = 6.0221409d23
  
  if (.not.associated(this%ERB_list)) return
  
100 format(100es18.8)

  option => this%realization%option
  output_option => this%realization%output_option
  
  fid = 91
  filename = PMUFDBOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")
       
  ! this time is set at the end of the reactive transport step
  write(fid,100,advance="no") option%time / output_option%tconv
  
  cur_ERB => this%ERB_list
  do
    if (.not.associated(cur_ERB)) exit
    write(fid,'(a32)',advance="no") trim(cur_ERB%name)
    write(fid,100,advance="no") cur_ERB%total_annual_dose   ! [Sv/yr]
    
    do k = 1,size(cur_ERB%annual_dose_supp_w_unsupp_rads)
      write(fid,100,advance="no") cur_ERB%annual_dose_supp_w_unsupp_rads(k)
    enddo                                  ! [Sv/yr]
    
    k = 0
    cur_supp_rad => this%supported_rad_list
    do
      if (.not.associated(cur_supp_rad)) exit
      k = k + 1
      write(fid,100,advance="no") &
          cur_ERB%annual_dose_supp_rad(k)  ! [Sv/yr]
      cur_supp_rad => cur_supp_rad%next
    enddo
    
    if (cur_ERB%incl_unsupported_rads) then
      k = 0
      cur_unsupp_rad => this%unsupported_rad_list
      do
        if (.not.associated(cur_unsupp_rad)) exit
        k = k + 1
        write(fid,100,advance="no") cur_ERB%annual_dose_unsupp_rad(k)  ! [Sv/yr]
        cur_unsupp_rad => cur_unsupp_rad%next
      enddo
    endif
    
    k = 0
    cur_supp_rad => this%supported_rad_list
    do
      if (.not.associated(cur_supp_rad)) exit
      k = k + 1
      write(fid,100,advance="no") &
          (cur_ERB%aqueous_conc_supported_rad(k)/0.001d0), &     ! [Bq/m3]
          ((cur_ERB%aqueous_conc_supported_rad(k))/ &    ! [mol/L] 
           (avagadro*cur_supp_rad%decay_rate))
      cur_supp_rad => cur_supp_rad%next
    enddo
  
    if (cur_ERB%incl_unsupported_rads) then
      k = 0
      cur_unsupp_rad => this%unsupported_rad_list
      do
        if (.not.associated(cur_unsupp_rad)) exit
        k = k + 1
        write(fid,100,advance="no") &
             cur_unsupp_rad%sorption_enhancement, &                 ! [-]
            (cur_ERB%aqueous_conc_unsupported_rad(k)/0.001d0), &    ! [Bq/m3]
            ((cur_ERB%aqueous_conc_unsupported_rad(k))/ &           ! [mol/L] 
             (avagadro*cur_unsupp_rad%decay_rate))
        cur_unsupp_rad => cur_unsupp_rad%next
      enddo
    endif
    
    cur_ERB => cur_ERB%next
  enddo
  
  close(fid)
  
end subroutine PMUFDBOutput

! *************************************************************************** !

subroutine PMUFDBOutputHeader(this)
  !
  ! Opens the output file and writes the header line.
  !
  ! Author: Jenn Frederick
  ! Date: 04/06/2017
  !
  
  use Output_Aux_module
  use Utility_module
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
  type(output_option_type), pointer :: output_option
  character(len=MAXWORDLENGTH) :: units_string
  character(len=MAXWORDLENGTH) :: variable_string
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXSTRINGLENGTH) :: filename
  class(ERB_base_type), pointer :: cur_ERB
  type(supported_rad_type), pointer :: cur_supp_rad
  type(unsupported_rad_type), pointer :: cur_unsupp_rad
  PetscInt :: fid, i
  PetscInt :: icolumn
  PetscBool :: exist
  
  output_option => this%realization%output_option
  
  fid = 91
  filename = PMUFDBOutputFilename(this%option)
  exist = FileExists(trim(filename))
  if (this%option%restart_flag .and. exist) return
  open(unit=fid,file=filename,action="write",status="replace")  
  
  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif 
  
  write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'
  
  cur_ERB => this%ERB_list
  do
    if (.not.associated(cur_ERB)) exit
    variable_string = 'ERB Model'
    units_string = ''
    cell_string = '(' // trim(cur_ERB%region_name) // ')'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Total Dose'
    units_string = 'Sv/yr'
    cell_string = ''
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
                             
    do i = 1,size(cur_ERB%annual_dose_supp_w_unsupp_rads)
      variable_string = 'Dose'                   
      units_string = 'Sv/yr'
      cell_string = trim(cur_ERB%names_supp_w_unsupp_rads(i))
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
    enddo
    
    cur_supp_rad => this%supported_rad_list
    do
      if (.not.associated(cur_supp_rad)) exit
      variable_string = trim(cur_supp_rad%name) // ' Dose'                    
      units_string = 'Sv/yr'
      cell_string = '' 
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      cur_supp_rad =>cur_supp_rad%next
    enddo
    
    if (cur_ERB%incl_unsupported_rads) then
      cur_unsupp_rad => this%unsupported_rad_list
      do
        if (.not.associated(cur_unsupp_rad)) exit
        variable_string = trim(cur_unsupp_rad%name) // '* Dose'
        units_string = 'Sv/yr'
        cell_string = ''
        call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                                 icolumn)
        cur_unsupp_rad => cur_unsupp_rad%next
      enddo
    endif
    
    cur_supp_rad => this%supported_rad_list
    do
      if (.not.associated(cur_supp_rad)) exit
      variable_string = trim(cur_supp_rad%name) // ' Aq. Conc.'
      units_string = 'Bq/m3'
      cell_string = '' 
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      variable_string = trim(cur_supp_rad%name) // ' Aq. Conc.'
      units_string = 'mol/L'
      cell_string = '' 
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      cur_supp_rad =>cur_supp_rad%next
    enddo
    
    if (cur_ERB%incl_unsupported_rads) then
      cur_unsupp_rad => this%unsupported_rad_list
      do
        if (.not.associated(cur_unsupp_rad)) exit
        variable_string = trim(cur_unsupp_rad%name) // '* E_i'
        units_string = ''
        cell_string = ''
        call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                                 icolumn)
        variable_string = trim(cur_unsupp_rad%name) // '* Aq. Conc.'
        units_string = 'Bq/m3'
        cell_string = ''
        call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                                 icolumn)
        variable_string = trim(cur_unsupp_rad%name) // '* Aq. Conc.'
        units_string = 'mol/L'
        cell_string = ''
        call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                                 icolumn)
        cur_unsupp_rad => cur_unsupp_rad%next
      enddo
    endif
                                 
    cur_ERB => cur_ERB%next
  enddo
  
  close(fid)
  
end subroutine PMUFDBOutputHeader

! ************************************************************************** !

function PMUFDBOutputFilename(option)
  ! 
  ! Generates filename for biosphere output file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/06/2017
  !

  use Option_module

  implicit none
  
  type(option_type), pointer :: option
  
  character(len=MAXSTRINGLENGTH) :: PMUFDBOutputFilename
  character(len=MAXWORDLENGTH) :: word

  write(word,'(i6)') option%myrank
  PMUFDBOutputFilename = trim(option%global_prefix) // &
                         trim(option%group_prefix) // &
                         '-' // trim(adjustl(word)) // '.bio'
  
end function PMUFDBOutputFilename  

! *************************************************************************** !

subroutine PMUFDBInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  ! 
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT
  
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  
end subroutine PMUFDBInputRecord

! *************************************************************************** !

subroutine PMUFDBDestroyERB(ERB)
  ! 
  ! Strips and destroys a ERB model in the UFD Biosphere process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/24/2017
  !
  
  use Utility_module, only : DeallocateArray

  implicit none
  
  class(ERB_base_type), pointer :: ERB
  
  call DeallocateArray(ERB%rank_list)
  call DeallocateArray(ERB%region_scaling_factor)
  call DeallocateArray(ERB%aqueous_conc_supp_rad_avg)
  call DeallocateArray(ERB%aqueous_conc_supported_rad)
  call DeallocateArray(ERB%aqueous_conc_unsupported_rad)
  call DeallocateArray(ERB%annual_dose_supp_rad)
  call DeallocateArray(ERB%annual_dose_unsupp_rad)
  call DeallocateArray(ERB%annual_dose_supp_w_unsupp_rads)
  deallocate(ERB%names_supp_w_unsupp_rads)
  nullify(ERB%region)
  nullify(ERB%next)
  deallocate(ERB)
  nullify(ERB)
  
end subroutine PMUFDBDestroyERB

! *************************************************************************** !

subroutine PMUFDBStrip(this)
  ! 
  ! Strips the UFD Biosphere process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
  class(ERB_base_type), pointer :: cur_ERB 
  class(ERB_BASE_type), pointer :: prev_ERB
  type(supported_rad_type), pointer :: cur_supp_rad
  type(supported_rad_type), pointer :: prev_supp_rad
  type(unsupported_rad_type), pointer :: cur_unsupp_rad
  type(unsupported_rad_type), pointer :: prev_unsupp_rad
  
  nullify(this%realization)
  
  cur_ERB => this%ERB_list
  do
    if (.not.associated(cur_ERB)) exit
    prev_ERB => cur_ERB
    cur_ERB => cur_ERB%next
    call PMUFDBDestroyERB(prev_ERB)
  enddo
  nullify(this%ERB_list)

  cur_supp_rad => this%supported_rad_list
  do
    if (.not.associated(cur_supp_rad)) exit
    prev_supp_rad => cur_supp_rad
    cur_supp_rad => cur_supp_rad%next
    nullify(prev_supp_rad%next)
    deallocate(prev_supp_rad)
    nullify(prev_supp_rad)
  enddo
  nullify(this%supported_rad_list)
  
  cur_unsupp_rad => this%unsupported_rad_list
  do
    if (.not.associated(cur_unsupp_rad)) exit
    prev_unsupp_rad => cur_unsupp_rad
    cur_unsupp_rad => cur_unsupp_rad%next
    nullify(prev_unsupp_rad%supported_parent)
    nullify(prev_unsupp_rad%next)
    deallocate(prev_unsupp_rad)
    nullify(prev_unsupp_rad)
  enddo
  nullify(this%unsupported_rad_list)

end subroutine PMUFDBStrip

! ************************************************************************** !

subroutine PMUFDBDestroy(this)
  ! 
  ! Strips and destroys the UFD Biosphere process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !

  implicit none
  
  class(pm_ufd_biosphere_type) :: this
  
  call PMUFDBStrip(this)
  
end subroutine PMUFDBDestroy

! ************************************************************************** !

end module PM_UFD_Biosphere_class
