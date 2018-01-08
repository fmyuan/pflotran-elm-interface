module EOSData_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Lookup_Table_module

  implicit none

  private

  PetscInt, parameter, public :: EOS_DENSITY = 1
  PetscInt, parameter, public :: EOS_ENTHALPY = 2
  PetscInt, parameter, public :: EOS_VISCOSITY = 3
  PetscInt, parameter, public :: EOS_INTERNAL_ENERGY = 4
  PetscInt, parameter, public :: EOS_FVF = 5
  PetscInt, parameter, public :: EOS_RS = 6
  PetscInt, parameter, public :: EOS_COMPRESSIBILITY = 7
  PetscInt, parameter, public :: EOS_VISCOSIBILITY = 8
  !add here other properties and derivatives, then increase MAX_PROP_NUM
  !and add the default user and internal units in EOSDataBaseInit
  PetscInt, parameter, public :: MAX_PROP_NUM = 8

  !type, public :: eos_data_base_type
  type, abstract, public :: eos_data_base_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: id
    PetscInt :: num_p   ! number of pressure intervals
    PetscInt :: num_t   ! number of temperature points
    PetscInt :: num_prop ! number of properties in the database
    PetscInt :: data_to_prop_map(MAX_PROP_NUM) !map the data idx to the prop.
    PetscInt :: prop_to_data_map(MAX_PROP_NUM) !map the prop to the dat idx
    PetscReal :: press_unit_conv_factor
    character(len=MAXWORDLENGTH) :: press_internal_units
    character(len=MAXWORDLENGTH) :: press_user_units
    PetscReal :: temp_unit_conv_factor
    character(len=MAXWORDLENGTH) :: temp_internal_units
    character(len=MAXWORDLENGTH) :: temp_user_units
    PetscReal :: prop_unit_conv_factors(MAX_PROP_NUM)
    character(len=MAXWORDLENGTH) :: prop_internal_units(MAX_PROP_NUM)
    character(len=MAXWORDLENGTH) :: prop_user_units(MAX_PROP_NUM)
    PetscReal, pointer :: data(:,:)
  contains
    procedure :: EOSDataBaseInit
    procedure, public :: EOSPropPresent
    procedure :: ReadUserUnits
    procedure :: UnitConversionFactors
    procedure :: ConvertFVFtoMolarDensity
    procedure :: EOSDataBaseStrip
    procedure, public :: SetMetricUnits
  end type

  type, public, extends(eos_data_base_type) :: eos_database_type
    character(len=MAXWORDLENGTH) :: file_name
    PetscReal :: dp      ! uniform pressure interval
    PetscReal :: dt      ! uniform temperature interval
    class(lookup_table_uniform_type), pointer :: lookup_table_uni
  contains
    procedure, public :: Read => EOSDatabaseRead
    procedure, public :: EOSProp => EOSPropLinearInterp
  end type

  type, public, extends(eos_data_base_type) :: eos_table_type
    PetscReal :: temperature !temperature value for isothermal data
    PetscInt :: n_indices !number of indices to be saved for lookup
    PetscInt :: first_index !location of first index in auxvars
    class(lookup_table_general_type), pointer :: lookup_table_gen
    class(eos_table_type), pointer :: next => null()
  contains
    procedure, public :: Read => EOSTableRead
    procedure, public :: EOSProp => EOSPropTable
  end type

  type, public :: eos_table_list_type
    PetscInt :: num_eos_tables
    class(eos_table_type), pointer :: first
    class(eos_table_type), pointer :: last
    class(eos_table_type), pointer :: array(:)
  end type eos_table_list_type

  type(eos_table_list_type), pointer, public :: eos_table_list => null()
  !type(eos_table_list_type), publict :: eos_table_list

  public :: EOSDatabaseCreate, &
            EOSDatabaseDestroy, &
            EOSTableCreate, &
            EOSTableDestroy, &
            EOSTableInitList, &
            EOSTableAddToList, &
            EOSTableProcessList, &
            EOSTableDestroyList



contains

! ************************************************************************** !

subroutine EOSDataBaseInit(this)

  implicit none

  class(eos_data_base_type) :: this


  PetscInt :: i_prop
  this%name = ''
  this%id = UNINITIALIZED_INTEGER
  this%num_p = UNINITIALIZED_INTEGER
  this%num_t = UNINITIALIZED_INTEGER
  this%num_prop = UNINITIALIZED_INTEGER
  this%data_to_prop_map(1:MAX_PROP_NUM) = UNINITIALIZED_INTEGER
  this%prop_to_data_map(1:MAX_PROP_NUM) = UNINITIALIZED_INTEGER
  this%press_unit_conv_factor = 1.0d0
  this%temp_unit_conv_factor = 1.0d0
  this%prop_unit_conv_factors(1:MAX_PROP_NUM) = 1.0d0
  !set deafult internal units - they can be changed by internal tranformation
  this%press_internal_units = 'Pa'
  this%temp_internal_units = 'C'
  this%prop_internal_units(EOS_DENSITY) = 'kg/m^3'
  this%prop_internal_units(EOS_ENTHALPY) = 'J/kg'
  this%prop_internal_units(EOS_VISCOSITY) = 'Pa-s' !should replace with Pa.s
  this%prop_internal_units(EOS_INTERNAL_ENERGY) = 'J/kg'
  this%prop_internal_units(EOS_FVF) = 'm^3/m^3'
  this%prop_internal_units(EOS_RS) = 'm^3/m^3'
  this%prop_internal_units(EOS_COMPRESSIBILITY) = '1/Pa'
  this%prop_internal_units(EOS_VISCOSIBILITY) = '1/Pa'
  !set deafult user units - identical to internal units
  this%press_user_units = 'Pa'
  this%temp_user_units = 'C'
  this%prop_user_units(EOS_DENSITY) = 'kg/m^3'
  this%prop_user_units(EOS_ENTHALPY) = 'J/kg'
  this%prop_user_units(EOS_VISCOSITY) = 'Pa-s' !should replace with Pa.s
  this%prop_user_units(EOS_INTERNAL_ENERGY) = 'J/kg'
  this%prop_user_units(EOS_FVF) = 'm^3/m^3'
  this%prop_user_units(EOS_RS) = 'm^3/m^3'
  this%prop_user_units(EOS_COMPRESSIBILITY) = '1/Pa'
  this%prop_user_units(EOS_VISCOSIBILITY) = '1/Pa'

  do i_prop = 1,MAX_PROP_NUM
    this%prop_internal_units(i_prop) = trim(this%prop_internal_units(i_prop))
    this%prop_user_units(i_prop) = trim(this%prop_user_units(i_prop))
  end do

  nullify(this%data)

end subroutine EOSDataBaseInit

! ************************************************************************** !

function EOSPropPresent(this,prop_iname)
  !
  ! Author: Paolo Orsini
  ! Date: 12/18/15
  !
  ! Checks if a property is defined in the database

  implicit none

  class(eos_data_base_type) :: this
  PetscInt, intent(in) :: prop_iname
  PetscBool :: EOSPropPresent

  EOSPropPresent = Initialized(this%prop_to_data_map(prop_iname))

end function EOSPropPresent

! ************************************************************************** !

subroutine ReadUserUnits(this,input,option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/28/17
  !
  ! Read unit EOS data unit and compute convert them to internal

  use Input_Aux_module
  use Option_module

  implicit none

  class(eos_data_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word, internal_units, user_units
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: prop_idx
  prop_idx = 0

  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    select case(keyword)
      case('PRESSURE')
        !internal_units = 'Pa'
        !this%press_unit_conv_factor = UnitReadAndConversionFactor(input, &
        !                                      internal_units,keyword,option)
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        this%press_user_units = trim(user_units)
      case('TEMPERATURE')
        !only Celsius currently supported - currently unit must be defined
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          if (trim(word) /= "C") then
            option%io_buffer = "EOS data only Temperatures in Celcius " // &
                                "are supported"
            call printErrMsg(option)
          end if
        end if
        this%temp_unit_conv_factor = 1.0
        this%temp_user_units = 'C'
        !else !it not temperature unit defined
        !  string = "EOS data " // trim(keyword) // ' units'
        !  call InputDefaultMsg(input,option,string)
        !  this%temp_unit_conv_factor = 1.0
        !endif
      case('DENSITY')
        !prop_idx = this%prop_to_data_map(EOS_DENSITY)
        !internal_units = 'kg/m^3'
        !this%prop_unit_conv_factors(prop_idx) = &
        !       UnitReadAndConversionFactor(input,internal_units,keyword,option)
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        this%prop_user_units(EOS_DENSITY) = user_units
      case('ENTHALPY')
        !prop_idx = this%prop_to_data_map(EOS_ENTHALPY)
        !internal_units = 'J/kg'
        !this%prop_unit_conv_factors(prop_idx) = &
        !       UnitReadAndConversionFactor(input,internal_units,keyword,option)
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        this%prop_user_units(EOS_ENTHALPY) = trim(user_units)
      case('INTERNAL_ENERGY')
        !prop_idx = this%prop_to_data_map(EOS_INTERNAL_ENERGY)
        !internal_units = 'J/kg'
        !this%prop_unit_conv_factors(prop_idx) = &
        !       UnitReadAndConversionFactor(input,internal_units,keyword,option)
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        this%prop_user_units(EOS_INTERNAL_ENERGY) = trim(user_units)
      case('VISCOSITY')
        !prop_idx = this%prop_to_data_map(EOS_VISCOSITY)
        !internal_units = 'Pa-s'
        !this%prop_unit_conv_factors(prop_idx) = &
        !       UnitReadAndConversionFactor(input,internal_units,keyword,option)
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        this%prop_user_units(EOS_VISCOSITY) = trim(user_units)
      case('FVF')
        !prop_idx = this%prop_to_data_map(EOS_FVF)
        !internal_units = 'm^3/m^3'
        !this%prop_unit_conv_factors(prop_idx) = &
        !       UnitReadAndConversionFactor(input,internal_units,keyword,option)
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        this%prop_user_units(EOS_FVF) = trim(user_units)
      case('RS')
         !prop_idx = this%prop_to_data_map(EOS_RS)
         !internal_units = 'm^3/m^3'
         !this%prop_unit_conv_factors(prop_idx) = &
         !     UnitReadAndConversionFactor(input,internal_units,keyword,option)
         call InputReadWord(input,option,user_units,PETSC_TRUE)
         this%prop_user_units(EOS_RS) = trim(user_units)
      case('COMPRESSIBILITY')
         !prop_idx = this%prop_to_data_map(EOS_COMPRESSIBILITY)
         !internal_units = '1/Pa'
         !this%prop_unit_conv_factors(prop_idx) = &
        !      UnitReadAndConversionFactor(input,internal_units,keyword,option)
         call InputReadWord(input,option,user_units,PETSC_TRUE)
         this%prop_user_units(EOS_COMPRESSIBILITY) = trim(user_units)
      case('VISCOSIBILITY')
         !prop_idx = this%prop_to_data_map(EOS_VISCOSIBILITY)
         !internal_units = '1/Pa'
         !this%prop_unit_conv_factors(prop_idx) = &
         !     UnitReadAndConversionFactor(input,internal_units,keyword,option)
         call InputReadWord(input,option,user_units,PETSC_TRUE)
         this%prop_user_units(EOS_VISCOSIBILITY) = trim(user_units)
      case default
        error_string = trim(error_string) // ': ' // this%name // &
        ': EOS DATA units'
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  end do

  call this%UnitConversionFactors(option)

end subroutine ReadUserUnits

! ************************************************************************** !
subroutine UnitConversionFactors(this,option)
  !
  ! Author: Paolo Orsini
  ! Date: 11/08/17
  !
  ! Process user input units and compute properties conversion fator

  use Option_module
  use Units_module

  implicit none

  class(eos_data_base_type) :: this
  type(option_type) :: option

  PetscInt :: i_prop, data_idx

  !covert pressure
  this%press_unit_conv_factor = UnitsConvertToInternal(this%press_user_units, &
                                  this%press_internal_units, option)

  !no temperature conversion - currently allowed only C

  do i_prop = 1, size(this%prop_user_units(:))
    if ( Initialized( this%prop_to_data_map(i_prop) ) ) then
      data_idx = this%prop_to_data_map(i_prop)
      !prop_idx = this%prop_to_data_map(EOS_DENSITY)
      this%prop_unit_conv_factors(data_idx) = &
             UnitsConvertToInternal(this%prop_user_units(i_prop), &
                                    this%prop_internal_units(i_prop),option)
    end if
  end do

end subroutine UnitConversionFactors

! ************************************************************************** !
subroutine SetMetricUnits(this,option)
  !
  ! Author: Paolo Orsini
  ! Date: 11/08/17
  !
  ! Set METRIC units - to read pvt tables

  use Option_module

  implicit none

  class(eos_data_base_type) :: this
  type(option_type) :: option

  this%press_user_units = 'Bar'
  this%temp_user_units = 'C'
  this%prop_user_units(EOS_DENSITY) = 'kg/m^3'
  this%prop_user_units(EOS_ENTHALPY) = 'kJ/kg'
  this%prop_user_units(EOS_VISCOSITY) = 'cP' !should replace with Pa.s
  this%prop_user_units(EOS_INTERNAL_ENERGY) = 'kJ/kg'
  this%prop_user_units(EOS_FVF) = 'm^3/m^3'
  this%prop_user_units(EOS_RS) = 'm^3/m^3'
  this%prop_user_units(EOS_COMPRESSIBILITY) = '1/Bar'
  this%prop_user_units(EOS_VISCOSIBILITY) = '1/Bar'

  !set viscosity internal units to Pa.s - this is needed because not all
  !internal viscosity units are set to Pa.s
  this%prop_internal_units(EOS_VISCOSITY) = 'Pa.s'

  call this%UnitConversionFactors(option)

end subroutine SetMetricUnits

! ************************************************************************** !
subroutine ConvertFVFtoMolarDensity(this,FMW,reference_density_kg)
  !
  ! Author: Paolo Orsini
  ! Date: 10/28/17
  !
  ! Replaces FVF with molar density

  implicit none

  class(eos_data_base_type) :: this
  PetscReal, intent(in) :: FMW
  PetscReal, intent(in) :: reference_density_kg

  PetscInt :: i_data

  do i_data = 1,size(this%data,2)
    ! kg/sm^3 * kmol/kg * sm^3/rm^3 = kmol/rm^3
    this%data(this%prop_to_data_map(EOS_FVF),i_data) = &
       reference_density_kg / FMW / &
       this%data(this%prop_to_data_map(EOS_FVF),i_data)
  end do
  !change variable name and map
  this%data_to_prop_map(this%prop_to_data_map(EOS_FVF)) = EOS_DENSITY
  this%prop_to_data_map(EOS_DENSITY) = this%prop_to_data_map(EOS_FVF)
  ! change units after conversion
  this%prop_internal_units(EOS_DENSITY) = 'kmol/m^3'

end subroutine ConvertFVFtoMolarDensity

! ************************************************************************** !

subroutine EOSDataBaseStrip(this)

  use Utility_module

  implicit none

  class(eos_data_base_type) :: this

  call DeallocateArray(this%data)

end subroutine EOSDataBaseStrip

! ************************************************************************** !

function EOSDatabaseCreate(filename,dbase_name)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !

  implicit none

  class(eos_database_type), pointer :: EOSDatabaseCreate
  character(len=MAXWORDLENGTH) :: filename
  character(len=*) :: dbase_name

  allocate(EOSDatabaseCreate)
  call EOSDatabaseCreate%EOSDataBaseInit()
  EOSDatabaseCreate%name = trim(dbase_name)
  EOSDatabaseCreate%file_name = trim(filename)
  nullify(EOSDatabaseCreate%lookup_table_uni)

end function EOSDatabaseCreate

! ************************************************************************** !

subroutine EOSDatabaseRead(this,option)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  ! Reads the the an EOS database from a text file for one or more
  ! phase properties.
  !
  ! Database format (for look up table)
  ! Header: NUM_P, NUM_T and DATA_LIST_ORDER
  ! DATA
  ! column 1 = pressure
  ! column 2 = temperature
  ! column 2+1, column 2+n: properties in the order given in DATA_LIST_ORDER
  !
  ! The dataset must be followed by NUM_P * NUM_T lines, not interrupted
  ! by a commented line.
  ! Each line must include P,T, and all properties listed in DATA_LIST_ORDER
  !
  ! each property listed depends on P and T, prop(P,T).
  !
  ! temperature values must be equispaced (same dt in the entire table)
  ! pressure values must be equispaced (same dp in the entire table)
  ! The data must be ordered so for growing P and T, with T looping faster.
  !
  ! This format repeates unnecessary P and T values in favours o readibility
  ! and creation of small dataset by hand. Inefficient for for large datasets.
  !

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(eos_database_type) :: this
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word, internal_units
  character(len=MAXSTRINGLENGTH) :: error_string = 'EOS_DATABASE'
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: prop_idx, prop_count, i_idx, j_idx
  PetscInt :: data_size
  PetscReal :: tempreal
  !PetscReal :: press_conv_factor, temp_conv_factor
  !PetscReal, pointer :: prop_unit_conv_factor(:)
  PetscBool :: pres_present, temp_present

  type(input_type), pointer :: input_table


  if (len_trim(this%file_name) < 1) then
    option%io_buffer = 'FILENAME must be specified for EOS_DATABASE.'
    call printErrMsg(option)
  endif

  input_table => InputCreate(IUNIT_TEMP,this%file_name,option)
  input_table%ierr = 0
  input_table%force_units = PETSC_FALSE

  !if ( option%myrank == 0 ) then
    option%io_buffer = 'Reading database = ' // this%file_name
    call printMsg(option)
  !end if

  !allocate(prop_unit_conv_factor(MAX_PROP_NUM))
  !prop_unit_conv_factor = UNINITIALIZED_DOUBLE

  !reading the database file header
  do

    call InputReadPflotranString(input_table,option)
    !if (InputCheckExit(input_table,option)) exit
    if (InputError(input_table)) exit

    call InputReadWord(input_table,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input_table,option,'keyword',error_string)
    call StringToUpper(keyword)
    select case(keyword)
      case('NUM_P')
        call InputReadInt(input_table,option,this%num_p)
        call InputErrorMsg(input_table,option,'number of dp',error_string)
      case('NUM_T')
        call InputReadInt(input_table,option,this%num_t)
        call InputErrorMsg(input_table,option,'number of dt',error_string)
      case('DATA_LIST_ORDER')
        pres_present = PETSC_FALSE; temp_present = PETSC_FALSE;
        prop_idx = 0
        do
          call InputReadPflotranString(input_table,option)
          if (InputCheckExit(input_table,option)) exit
          call InputReadWord(input_table,option,word,PETSC_TRUE)
          select case(word)
            case('PRESSURE')
              pres_present = PETSC_TRUE
              !internal_units = 'Pa'
              !press_conv_factor = UnitReadAndConversionFactor(input_table, &
              !                                    internal_units,word,option)
            case('TEMPERATURE')
              temp_present = PETSC_TRUE
              !only C currently supported
              !internal_units = 'C'
              !temp_conv_factor = UnitReadAndConversionFactor(input_table, &
              !                                    internal_units,word,option)
            case('DENSITY')
              prop_idx = prop_idx + 1
              this%data_to_prop_map(prop_idx) = EOS_DENSITY
              this%prop_to_data_map(EOS_DENSITY) = prop_idx
              !internal_units = 'kg/m^3'
              !prop_unit_conv_factor(prop_idx) = &
              !                UnitReadAndConversionFactor(input_table, &
              !                               internal_units,word,option)
            case('ENTHALPY')
              prop_idx = prop_idx + 1
              this%data_to_prop_map(prop_idx) = EOS_ENTHALPY
              this%prop_to_data_map(EOS_ENTHALPY) = prop_idx
              !internal_units = 'J/kg'
              !prop_unit_conv_factor(prop_idx) = &
              !                UnitReadAndConversionFactor(input_table, &
              !                               internal_units,word,option)
            case('INTERNAL_ENERGY')
              prop_idx = prop_idx + 1
              this%data_to_prop_map(prop_idx) = EOS_INTERNAL_ENERGY
              this%prop_to_data_map(EOS_INTERNAL_ENERGY) = prop_idx
              !internal_units = 'J/kg'
              !prop_unit_conv_factor(prop_idx) = &
              !                UnitReadAndConversionFactor(input_table, &
              !                               internal_units,word,option)
            case('VISCOSITY')
              prop_idx = prop_idx + 1
              this%data_to_prop_map(prop_idx) = EOS_VISCOSITY
              this%prop_to_data_map(EOS_VISCOSITY) = prop_idx
              !internal_units = 'Pa-s'
              !prop_unit_conv_factor(prop_idx) = &
              !                UnitReadAndConversionFactor(input_table, &
              !                               internal_units,word,option)
            case default
              error_string = trim(error_string) // ': ' // this%file_name // &
              ': DATA_LIST_ORDER'
              call InputKeywordUnrecognized(keyword,error_string,option)
          end select
        end do
        this%num_prop = prop_idx
        ! go back to DATA_LIST_ORDER - read variable again to comput
        ! unit covnertion factors
        string = "DATA_LIST_ORDER"
        call InputFindStringInFile(input_table,option,string)
        call this%ReadUserUnits(input_table,option)

        if (.not.pres_present) then
          option%io_buffer = 'PRESSURE must be present in any EOS_DATABASE.'
          call printErrMsg(option)
        end if
        if (.not.temp_present) then
          option%io_buffer = 'TEMPERATURE must be present in any EOS_DATABASE.'
          call printErrMsg(option)
        end if
      case('DATA')
        exit
      case default
        error_string = trim(error_string) // ': ' // this%file_name
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select

  end do

  data_size = this%num_p * this%num_t

  this%num_prop = prop_idx
  allocate(this%data(prop_idx,data_size))

  !create lookup table
  this%lookup_table_uni => LookupTableCreateUniform(TWO_INTEGER)
  this%lookup_table_uni%dims(1) = this%num_t
  this%lookup_table_uni%dims(2) = this%num_p
  allocate(this%lookup_table_uni%axis1%values(this%num_t))
  this%lookup_table_uni%axis1%values(1:this%num_t) = UNINITIALIZED_DOUBLE
  allocate(this%lookup_table_uni%axis2%values(this%num_p))
  this%lookup_table_uni%axis2%values(1:this%num_p) = UNINITIALIZED_DOUBLE

  !TODO
  ! start loading data - at the moment using Input facility, however this file
  ! can be large. TODO. Implement more efficient solutions:
  ! - using read(,) without InputReadDouble filter
  ! - adding the option of reading a .h5 file where the database is defined

  ! go to data - first time to load axis1 and 2 values
  string = "DATA"
  call InputFindStringInFile(input_table,option,string)

  do j_idx = 1,this%num_p

    do i_idx = 1, this%num_t
      call InputReadPflotranString(input_table,option)

      call InputReadDouble(input_table,option,tempreal)
      call InputErrorMsg(input_table,option, &
                           'VALUE', 'EOS_DATABASE PRESS_VALUE')
      ! convert pressure to internal units - Pa
      this%lookup_table_uni%axis2%values(j_idx) = &
          tempreal * this%press_unit_conv_factor
      !reads temperature values - repeated this%num_p times - not efficient
      call InputReadDouble(input_table,option, &
                           this%lookup_table_uni%axis1%values(i_idx))
      call InputErrorMsg(input_table,option, &
                         'VALUE', 'EOS_DATABASE TEMP_VALUE')

      prop_count = i_idx + (j_idx-1) * this%num_t
      do prop_idx = 1,this%num_prop
        call InputReadDouble(input_table,option,this%data(prop_idx,prop_count))
        call InputErrorMsg(input_table,option,&
                           'VALUE','EOS_DATABASE PROP_VALUE')
        !convert to intenral units
        this%data(prop_idx,prop_count) = this%data(prop_idx,prop_count) * &
                                         this%prop_unit_conv_factors(prop_idx)
                                         !prop_unit_conv_factor(prop_idx)
      end do

    end do

  end do

  !deallocate(prop_unit_conv_factor)
  call InputDestroy(input_table)

end subroutine EOSDatabaseRead

! ************************************************************************** !

subroutine EOSPropLinearInterp(this,T,P,prop_iname,prop_value,ierr)
  !
  ! Author: Paolo Orsini
  ! Date: 12/12/15
  !
  ! interpolates a single EOS property from the EOS database
  ! Note: when more properties must be extracted from the same EOS database,
  !       i.e. properties listed for the same values and range of P and T,
  !       all propoertis should be extracted at the same time to perform
  !       only once the look up operations, which are:
  !       (a) LookupTableIndexUniform, i.e. axis look up
  !       (b) P,T location checks within a single 2D domain, (p,p+dp; t,t+dt)
  !           currently done within LookupTableInterpolate2DUniform
  !
  !       TODO: add a method lookup_table_uni to extract mulitdimensional data
  !             at the same time (data plus function)

  implicit none

  class(eos_database_type) :: this
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscInt, intent(in) :: prop_iname
  PetscReal, intent(out) :: prop_value ! database units (SI)
  PetscErrorCode, intent(out) :: ierr

  !PetscReal :: EOSEOSProp !property given in units specified in the database

  ierr = 0
  !insert check P and/or T out of range
  ! T smaller then minimum T in the database
  ! errors here a passed back to the EOSxxxPhase module
  ! and to ModeXXXComputeAuxVar where cells ids are available to check where
  ! P and T are out of range.
  ! TODO: (i) define error id identifiers for each EOSPhase, or for all EOSs.
  !       (ii) define a function that handles errors ids and print error message.
  !       (iii) can then remove the print statment below
  !       note: need to minimise the number if-checks for efficiency.
  if ( T < this%lookup_table_uni%axis1%values(1) ) then
    ierr = 101
    print*, "EOSEOSProp - T smaller than min val in EOSdatabase"
    print*, "Temp val [°C] = ", T
    stop
  end if
  if ( T > this%lookup_table_uni%axis1%values(this%num_t) ) then
    ierr = 102
    print*, "EOSEOSProp - T larger than max val in EOSdatabase"
    print*, "Temp val [°C] = ", T
    stop
  end if
  if ( P < this%lookup_table_uni%axis2%values(1) ) then
    ierr = 103
    print*, "EOSEOSProp - P smaller than min val in EOSdatabase"
    print*, "Press val [Mpa] = ", P*1.d-6
    stop
  end if
  if ( P > this%lookup_table_uni%axis2%values(this%num_p) ) then
    ierr = 104
    print*, "EOSEOSProp - P larger than max val in EOSdatabase"
    print*, "Press val [Mpa] = ", P*1.d-6
    stop
  end if

  this%lookup_table_uni%data => this%data(this%prop_to_data_map(prop_iname),:)

                             !T       P      optional
  !this%lookup_table_uni%Sample(lookup1,lookup2,lookup3)
  prop_value = this%lookup_table_uni%Sample(T,P)

  nullify(this%lookup_table_uni%data)

end subroutine EOSPropLinearInterp

! ************************************************************************** !

subroutine EOSDatabaseDestroy(eos_database)
  !
  ! Author: Paolo Orsini
  ! Date: 12/14/15
  !
  ! destroys EOS database

  use Utility_module

  implicit none

  class(eos_database_type), pointer :: eos_database

  if (.not.associated(eos_database)) return

  call eos_database%EOSDataBaseStrip()
  call LookupTableDestroy(eos_database%lookup_table_uni)

  deallocate(eos_database)
  nullify(eos_database)

end subroutine EOSDatabaseDestroy

! ************************************************************************** !

function EOSTableCreate(table_name,option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/18/17
  !
  ! Create EOS table

  use Option_module

  implicit none

  class(eos_table_type), pointer :: EOSTableCreate
  character(len=*) :: table_name
  type(option_type) :: option

  allocate(EOSTableCreate)
  call EOSTableCreate%EOSDataBaseInit()
  EOSTableCreate%name = trim(table_name)
  EOSTableCreate%n_indices = 0
  EOSTableCreate%first_index = 0
  nullify(EOSTableCreate%lookup_table_gen)

  !as default set up Metric units
  !call EOSTableCreate%SetMetricUnits(option)

end function EOSTableCreate

! ************************************************************************** !

subroutine EOSTableRead(this,input,option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/18/17
  !
  ! Reads the PVT table data

  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module

  implicit none

  class(eos_table_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  !character(len=MAXWORDLENGTH) :: keyword, word, input_unit, internal_units
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string = 'EOS_PVT_TABLE'
  !type(input_type), pointer :: input_table
  !type(input_type), pointer :: input_tmp
  PetscReal, pointer :: press_data_array(:,:)
  PetscInt :: n_temp_count
  PetscInt :: i_temp, i_press, i_prop
  PetscInt :: n_press_count, n_press_count_first, n_press_count_cur
  PetscInt :: temp_array_size, press_array_size
  PetscReal, pointer :: temp_array(:)

  n_press_count_first = 0
  n_temp_count = 0
  n_press_count = 0
  n_press_count_cur = 0
  temp_array_size = 100 !estimate for temperature points
  allocate(temp_array(temp_array_size))
  temp_array = UNINITIALIZED_DOUBLE
  press_array_size = 1000 !estimate of pressure points from all tables
  ! num_prop+1 = + 1 make space for the pressure
  allocate(press_data_array(this%num_prop+1,press_array_size))
  press_data_array = UNINITIALIZED_DOUBLE

  !reading pvt table
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    if (InputError(input)) then
      option%io_buffer = 'Error found as reading PVT table'
      call printErrMsg(option)
    end if
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    select case(keyword)
      case('DATA_UNITS')
        call this%ReadUserUnits(input,option)
      case('DATA')
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit
          if (InputError(input)) then
            option%io_buffer = 'Error found as reading PVT table - DATA block'
            call printErrMsg(option)
          end if
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'word',error_string)
          call StringToUpper(keyword)
          select case(word)
            case('TEMPERATURE')
              n_temp_count = n_temp_count + 1
              if (n_temp_count > temp_array_size) then
                call reallocateRealArray(temp_array,temp_array_size)
              end if
              call InputReadDouble(input,option,temp_array(n_temp_count))
              call InputErrorMsg(input,option,&
                                 'VALUE','EOS_PVT_TABLE - TEMPERATURE')
              !read pvt data for a given temperature
              !call ReadPressureTable(input,option,this%num_prop,&
              !                       n_press_count_cur,data_array)
              call ReadPressureTable(input,option,n_press_count, &
                                     n_press_count_cur,press_data_array)
              if (n_temp_count == 1) then
                n_press_count_first = n_press_count_cur
              else
                if ( n_press_count_cur /= n_press_count_first ) then
                  option%io_buffer =  'PVT Table = ' // 'this%name ' // &
                                      "PVT tables with pressure points that" // &
                                      "varies with T not currently supported"
                  call printErrMsg(option)
                end if
              end if
            case default
              option%io_buffer = 'PVT Table DATA block must contain ' // &
                                 'TEMPERATURE blocks'
              call printErrMsg(option)
          end select
        end do
      case default
        option%io_buffer = 'PVT Table = ' // 'this%name'
    end select
  end do

  this%num_t = n_temp_count
  this%num_p = n_press_count

  if (this%num_t ==1) then !isothermal press_data_array
    this%temperature = temp_array(1)
    ! create general 1D lookup table
    this%lookup_table_gen => LookupTableCreateGeneral(ONE_INTEGER)
    this%lookup_table_gen%dim = 1
    this%lookup_table_gen%dims(1) = this%num_p
    allocate(this%lookup_table_gen%axis1%values(this%num_p))
    do i_press = 1,this%num_p
      this%lookup_table_gen%axis1%values(i_press) = &
         press_data_array(1,i_press) * this%press_unit_conv_factor
    end do
    this%n_indices = 1
  else if (this%num_t > 1) then
    !creat general 2D lookip table and load the temperature values in axis1
    this%temperature = UNINITIALIZED_DOUBLE
    this%lookup_table_gen => LookupTableCreateGeneral(TWO_INTEGER)
    this%lookup_table_gen%dim = 2
    this%lookup_table_gen%dims(1) = this%num_t
    this%lookup_table_gen%dims(2) = n_press_count_first
    allocate(this%lookup_table_gen%axis1%values(this%num_t))
    do i_temp = 1,this%num_t
      this%lookup_table_gen%axis1%values(i_temp) = temp_array(i_temp) * &
                                                   this%temp_unit_conv_factor
    end do
    allocate(this%lookup_table_gen%axis2%values(this%num_p))
    do i_press = 1,this%num_p
      this%lookup_table_gen%axis2%values(i_press) = &
          press_data_array(1,i_press) * this%press_unit_conv_factor
    end do
    this%n_indices = 3
  end if
  !allocate general lookup Table
  !this%lookup_table_gen => LookupTableCreateGeneral(dim)
  allocate(this%data(this%num_prop,n_press_count))
  do i_press = 1,n_press_count
    do i_prop=1,this%num_prop
      this%data(i_prop,i_press) = press_data_array(i_prop+1,i_press) * &
                                  this%prop_unit_conv_factors(i_prop)
    end do
  end do

  call DeallocateArray(press_data_array)
  call DeallocateArray(temp_array)

  !if (len_trim(this%file_name) >= 1) then
  !  call InputDestroy(input_table)
  !end if
  !nullify(input_tmp)

end subroutine EOSTableRead

! ************************************************************************** !

subroutine ReadPressureTable(input,option,press_idx,n_press,press_data_array)
  !
  ! Author: Paolo Orsini
  ! Date: 10/18/17
  !
  ! Reads the PVT fields vs pressure
  !
  ! "*" entries not yet supported
  !PO TODO: when ecountering "*" as entry in the data field interpolates
  ! using previous and following data field.
  ! read the field data in each line using InputReadDouble,
  ! if getting an error, read as "word" and check if "*"
  ! compare each word vs "*" using function "StringCompare"
  ! Wehn finning "*", map location into an array to poste post-process later
  ! for the interpolation. "*" should not be entered in the first and last
  ! record
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module

  implicit none

  !class(eos_database_type) :: this
  !PetscInt, intent(in) :: num_fields
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscInt, intent(inout) :: press_idx
  PetscInt, intent(out) :: n_press
  PetscReal, pointer :: press_data_array(:,:)

  PetscInt :: i_data
  PetscInt :: size_rank2
  PetscInt :: num_fields !this is size_rank1 - include the pressure field

  num_fields = size(press_data_array,1)
  size_rank2 = size(press_data_array,2)

  n_press = 0
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    if (InputError(input)) then
      option%io_buffer = 'Error PVT table - reading Pressure Table'
      call printErrMsg(option)
    end if
    n_press = n_press + 1
    press_idx = press_idx + 1
    if ( press_idx > size_rank2 ) then
      !each time doubles the size of rank 2
      !tmp_array_size overwritten by new size
      call reallocateRealArray(press_data_array,size_rank2)
    end if
    do i_data = 1, num_fields
      call InputReadDouble(input,option,press_data_array(i_data,press_idx))
      call InputErrorMsg(input,option,&
                         'VALUE','EOS_PVT_TABLE - PRESSURE TABLE VALUE')
    end do
  end do

end subroutine ReadPressureTable

! ************************************************************************** !

!subroutine EOSPropTable(this,T,P,prop_iname,prop_value,ierr,iP1,iT,iP2)
subroutine EOSPropTable(this,T,P,prop_iname,prop_value,indices,ierr)
  !
  ! Author: Paolo Orsini
  ! Date: 10/20/17
  !
  ! interpolates a single EOS property from a PVT Table
  !
  ! if table%dim == 1, only pressure lookup
  !  indices(eos_table%first_index) = Pressure index
  ! else if table%dim == 2, one temperature lookup and two pressure lookups
  !  indices(eos_table%first_index+1) = iP1, i.e. Pressure_index_1
  !  indices(eos_table%first_index+2) = iP2, i.e. Pressure_index_2
  ! the case table%dim == 2, with one temperature lookup and
  !                       one pressure lookups not supported
  !
  ! Note: when more properties must be extracted from the same PVT table,
  !       i.e. properties listed for the same values and range of P and T,
  !       all propoertis should be extracted at the same time to perform
  !       only once the look up operations, which are:
  !       (a) LookupTableIndexGeneral, i.e. axis look up
  !       (b) P,T location checks within a single 2D domain, (p,p+dp; t,t+dt)
  !           currently done within LookupTableInterpolate2DGeneral
  !
  !    PO TODO: add a method lookup_table_gen to extract mulitdimensional data
  !             at the same time (data plus function)

  implicit none

  class(eos_table_type) :: this
  PetscReal, intent(in) :: T              ! temperature [C]
  PetscReal, intent(in) :: P              ! pressure [Pa]
  !PetscInt, intent(inout), optional :: iP1   ! pressure looup index1 - left T
  !PetscInt, intent(inout), optional :: iT    ! temperature lookup index
  !PetscInt, intent(inout), optional :: iP2   ! pressure looup index2 - right T
  PetscInt, intent(in) :: prop_iname
  PetscReal, intent(out) :: prop_value ! database units (SI)
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, intent(inout) :: indices(:)

  ierr = 0
  !no errors are assigned if P,T go out of range, since the values
  !are extrapolated using the edge values - however warnning are printed
  !on the screen to signal extrapolation

  this%lookup_table_gen%data => this%data(this%prop_to_data_map(prop_iname),:)

  if (this%lookup_table_gen%dim == ONE_INTEGER) then
    if ( P < this%lookup_table_gen%axis1%values(1) ) then
      print*, "EOSPropTable - P smaller than min val in table"
      print*, "table name = ", this%name
      print*, "Property extrapolated for Press val [Mpa] = ", P*1.d-6
    end if
    if ( P > this%lookup_table_gen%axis1%values(this%num_p) ) then
      print*, "EOSPropTable - P larger than max val in table"
      print*, "table name = ", this%name
      print*, "Property extrapolated for Press val [Mpa] = ", P*1.d-6
    end if
    this%lookup_table_gen%axis1%saved_index = indices(this%first_index)
    !T       P      optional
    !this%lookup_table_gen%Sample(lookup1,lookup2,lookup3)
    prop_value = this%lookup_table_gen%Sample(P)
    !must copy back iP1 because
    ! saved_index = indices(eos_table%first_index) is a copying operation not pointing
    ! saved_index should be a pointer to save an assigment operation
    indices(this%first_index) = this%lookup_table_gen%axis1%saved_index
  else if (this%lookup_table_gen%dim == TWO_INTEGER) then
    if ( T < this%lookup_table_gen%axis1%values(1) ) then
      print*, "EOSPropTable - T smaller than min val in table"
      print*, "table name = ", this%name
      print*, "Property extrapolated for Temp val [°C] = ", T
    end if
    if ( T > this%lookup_table_gen%axis1%values(this%num_t) ) then
      print*, "EOSPropTable - T larger than max val in table"
      print*, "table name = ", this%name
      print*, "Property extrapolated Temp val [°C] = ", T
    end if
    ! Due to the possible irregular shape of Pmax(T) and Pmin(T)
    ! a out of bound check requires the computation of Pmin & Pmax in the table
    ! the check below is valid if:
    ! Pmin = axis1%values(1) and Pmax=axis2%values(this%num_p)
    !if ( P < this%lookup_table_gen%axis1%values(1) ) then
    !  print*, "EOSEOSProp - P smaller than min val in EOSdatabase"
    !  print*, "Property extrapolated for Press val [Mpa] = ", P*1.d-6
    !end if
    !if ( P > this%lookup_table_gen%axis2%values(this%num_p) ) then
    !  print*, "EOSEOSProp - P larger than max val in EOSdatabase"
    !  print*, "Property extrapolated for Press val [Mpa] = ", P*1.d-6
    !end if
    !this%lookup_table_gen%axis1%saved_index = iT
    !this%lookup_table_gen%axis2%saved_index = iP1
    this%lookup_table_gen%axis1%saved_index = indices(this%first_index)
    this%lookup_table_gen%axis2%saved_index = indices(this%first_index+1)
    this%lookup_table_gen%axis2%saved_index2 = indices(this%first_index+2)
    !if (present(iP2)) then
    !  this%lookup_table_gen%axis2%saved_index2 = iP2
    !else
    !  this%lookup_table_gen%axis2%saved_index2 = iP1
    !end if
    !T       P      optional
    !this%lookup_table_gen%Sample(lookup1,lookup2,lookup3)
    prop_value = this%lookup_table_gen%Sample(T,P)
    !must copy back iT, iP1 and iP2 because
    ! saved_index = iT is a copying operation not a pointer
    ! saved_index should be a pointer to save an assigment operation
    !iT = this%lookup_table_gen%axis1%saved_index
    !iP1 = this%lookup_table_gen%axis2%saved_index
    !if (present(iP2)) then
    !  iP2 = this%lookup_table_gen%axis2%saved_index2
    !else
    !  iP2 = iP1
    !end if
    indices(this%first_index) = this%lookup_table_gen%axis1%saved_index
    indices(this%first_index+1)= this%lookup_table_gen%axis2%saved_index
    indices(this%first_index+2) = this%lookup_table_gen%axis2%saved_index2
  end if

  nullify(this%lookup_table_gen%data)

end subroutine EOSPropTable

! ************************************************************************** !

subroutine EOSTableDestroy(eos_table)
  !
  ! Author: Paolo Orsini
  ! Date: 10/21/17
  !
  ! destroys EOS Table

  use Utility_module

  implicit none

  class(eos_table_type), pointer :: eos_table

  if (.not.associated(eos_table)) return

  call eos_table%EOSDataBaseStrip()
  call LookupTableDestroy(eos_table%lookup_table_gen)

  deallocate(eos_table)
  nullify(eos_table)

end subroutine EOSTableDestroy

! ************************************************************************** !

subroutine EOSTableInitList()
  !
  ! Initializes eos_table_list
  !
  ! Author: Paolo Orsini
  ! Date: 10/24/17
  !

  implicit none

  !type(eos_table_list_type), pointer :: list

  allocate(eos_table_list)

  nullify(eos_table_list%first)
  nullify(eos_table_list%last)
  nullify(eos_table_list%array)
  eos_table_list%num_eos_tables = 0

end subroutine EOSTableInitList

! ************************************************************************** !

subroutine EOSTableAddToList(new_eos_table,list)
  !
  ! Adds a new eos_table to an eos table list
  !
  ! Author: Paolo Orsini
  ! Date: 10/24/17
  !
  implicit none

  class(eos_table_type), pointer :: new_eos_table
  type(eos_table_list_type) :: list

  list%num_eos_tables = list%num_eos_tables + 1
  new_eos_table%id = list%num_eos_tables
  if (.not.associated(list%first)) list%first => new_eos_table
  if (associated(list%last)) list%last%next => new_eos_table
  list%last => new_eos_table

end subroutine EOSTableAddToList

! ************************************************************************** !

subroutine EOSTableProcessList(option)
  !
  ! Loop through EOS table list and build indices
  !
  ! Author: Paolo Orsini
  ! Date: 10/24/17

  use Option_module

  implicit none

  type(option_type) :: option

  type(eos_table_list_type), pointer :: list
  class(eos_table_type), pointer :: eos_table

  list => eos_table_list

  if (.not.associated(list)) return

  option%neos_table_indices = 0
  !loop over  EOS tables and create a map to store save indices in auxvars
  eos_table => list%first
  do
    if (.not.(associated(eos_table))) exit
    eos_table%first_index = option%neos_table_indices + 1
    option%neos_table_indices = option%neos_table_indices + eos_table%n_indices
    ! add here other operation that must be performed on all eos tables
    ! move to next table if it is associated
    if( .not.(associated(eos_table%next))) exit
    eos_table => eos_table%next
  enddo

end subroutine EOSTableProcessList

! ************************************************************************** !

subroutine EOSTableDestroyList()
  !
  ! Destroy a list of EOS tables
  !
  ! Author: Paolo Orsini
  ! Date: 10/24/17
  !

  implicit none

  type(eos_table_list_type), pointer :: list

  class(eos_table_type), pointer :: eos_table, prev_eos_table

  if (.not.associated(list)) return

  list => eos_table_list

  eos_table => list%first
  do
    if (.not.associated(eos_table)) exit
    prev_eos_table => eos_table
    eos_table => eos_table%next
    call EOSTableDestroy(prev_eos_table)
    nullify(prev_eos_table)
  enddo

  list%num_eos_tables = 0
  nullify(list%first)
  nullify(list%last)
  if (associated(list%array)) deallocate(list%array)
  nullify(list%array)

  deallocate(list)
  nullify(list)

  nullify(eos_table_list)

end subroutine EOSTableDestroyList

! ************************************************************************** !


end module EOSData_module
