module Fracture_module
  
  use PFLOTRAN_Constants_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none
  
  private

  PetscInt, parameter, public :: frac_init_pres_index = 1
  PetscInt, parameter, public :: frac_alt_pres_index = 2
  PetscInt, parameter, public :: frac_max_poro_index = 3
  PetscInt, parameter, public :: frac_poro_exp_index = 4
  PetscInt, parameter, public :: frac_change_perm_x_index = 1
  PetscInt, parameter, public :: frac_change_perm_y_index = 2
  PetscInt, parameter, public :: frac_change_perm_z_index = 3
  
  type, public :: fracture_type
    PetscReal :: init_pressure
    PetscReal :: altered_pressure
    PetscReal :: maximum_porosity
    PetscReal :: porosity_exponent
    PetscReal :: change_perm_x
    PetscReal :: change_perm_y
    PetscReal :: change_perm_z
  contains
    procedure, public :: Read => FractureRead
  end type fracture_type
  
!  class(fracture_type), pointer, public :: fracture

  public :: FractureInit, &
            FractureCreate, &
            FractureSetInitialPressure, &
            FractureAuxVarInit, &
            FracturePropertytoAux, &
            FractureDestroy, &
            FracturePoroEvaluate, &
            FracturePermScale
  
  contains

! ************************************************************************** !

function FractureCreate()
  !
  ! Author: Heeho Park
  ! Date: 4/7/15
  !

  implicit none
  
  class(fracture_type), pointer :: FractureCreate
  class(fracture_type), pointer :: fracture
  
  allocate(fracture)
  call FractureInit(fracture)
  
  FractureCreate => fracture
  
end function FractureCreate

! ************************************************************************** !

subroutine FractureInit(this)
  !
  ! Author: Heeho Park
  ! Date: 4/7/2015
  !

  implicit none
  
  class(fracture_type), pointer :: this
  
  this%init_pressure = UNINITIALIZED_DOUBLE
  this%altered_pressure = UNINITIALIZED_DOUBLE
  this%maximum_porosity = UNINITIALIZED_DOUBLE
  this%porosity_exponent = UNINITIALIZED_DOUBLE
  this%change_perm_x = 0.d0
  this%change_perm_y = 0.d0
  this%change_perm_z = 0.d0

end subroutine FractureInit

! ************************************************************************** !

subroutine FractureAuxVarInit(auxvar)
  !
  ! Author: Heeho Park, Glenn Hammond
  ! Date: 7/8/2015, 6/15/17
  !

  use Material_Aux_class
  
  implicit none
  
  class(material_auxvar_type), intent(inout) :: auxvar

  call MaterialAuxVarFractureStrip(auxvar%fracture)
  allocate(auxvar%fracture)
  auxvar%fracture%initial_pressure = UNINITIALIZED_DOUBLE
  auxvar%fracture%properties = UNINITIALIZED_DOUBLE
  auxvar%fracture%vector = 0.d0

end subroutine FractureAuxVarInit

! ************************************************************************** !

subroutine FracturePropertytoAux(fracture_auxvar,fracture_property)
  !
  ! Author: Heeho Park
  ! Date: 7/8/2015
  !

  use Material_Aux_class
  
  implicit none

  type(fracture_auxvar_type), pointer :: fracture_auxvar
  class(fracture_type), pointer :: fracture_property

  if (associated(fracture_auxvar)) then
    if (associated(fracture_property)) then
      fracture_auxvar%fracture_is_on = PETSC_TRUE
      fracture_auxvar%properties(frac_init_pres_index) = &
        fracture_property%init_pressure
      fracture_auxvar%properties(frac_alt_pres_index) = &
        fracture_property%altered_pressure
      fracture_auxvar%properties(frac_max_poro_index) = &
        fracture_property%maximum_porosity
      fracture_auxvar%properties(frac_poro_exp_index) = &
        fracture_property%porosity_exponent
      fracture_auxvar%vector(frac_change_perm_x_index) = &
        fracture_property%change_perm_x
      fracture_auxvar%vector(frac_change_perm_y_index) = &
        fracture_property%change_perm_y
      fracture_auxvar%vector(frac_change_perm_z_index) = &
        fracture_property%change_perm_z
    else
      fracture_auxvar%fracture_is_on = PETSC_FALSE
      fracture_auxvar%properties = UNINITIALIZED_DOUBLE
      fracture_auxvar%vector = 0.d0
    endif
  endif

end subroutine FracturePropertytoAux

! ************************************************************************** !

subroutine FractureRead(this,input,option)
  ! 
  ! Author: Heeho Park
  ! Date: 4/7/15
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Input_Aux_module
  use String_module
  
  implicit none
  
  class(fracture_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word

  option%flow%fracture_on = PETSC_TRUE
  
  do
      call InputReadPflotranString(input,option)
      call InputReadStringErrorMsg(input,option, &
                                    'MATERIAL_PROPERTY,WIPP-FRACTURE')
          
      if (InputCheckExit(input,option)) exit
          
      if (InputError(input)) exit
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'keyword', &
                          'MATERIAL_PROPERTY,WIPP-FRACTURE')   
      select case(trim(word))
        case('INITIATING_PRESSURE')
          call InputReadDouble(input,option, &
                                this%init_pressure)
          call InputErrorMsg(input,option, &
                              'initiating pressure of fracturing', &
                              'MATERIAL_PROPERTY,WIPP-FRACTURE')
        case('ALTERED_PRESSURE')
          call InputReadDouble(input,option, &
                                this%altered_pressure)
          call InputErrorMsg(input,option, &
                              'altered pressure of fracturing', &
                              'MATERIAL_PROPERTY,WIPP-FRACTURE')
        case('MAXIMUM_FRACTURE_POROSITY')
          call InputReadDouble(input,option, &
                                this%maximum_porosity)
          call InputErrorMsg(input,option, &
                              'maximum fracture porosity', &
                              'MATERIAL_PROPERTY,WIPP-FRACTURE')
        case('FRACTURE_EXPONENT')
          call InputReadDouble(input,option, &
                              this%porosity_exponent)
          call InputErrorMsg(input,option, &
                          'dimensionless fracture exponent for porosity', &
                              'MATERIAL_PROPERTY,WIPP-FRACTURE')
        case('ALTER_PERM_X')
          this%change_perm_x = 1.d0
        case('ALTER_PERM_Y')
          this%change_perm_y = 1.d0
        case('ALTER_PERM_Z')
          this%change_perm_z = 1.d0
        case default
          call InputKeywordUnrecognized(word, &
                  'MATERIAL_PROPERTY,WIPP-FRACTURE',option)
      end select
    enddo

end subroutine FractureRead

! ************************************************************************** !

subroutine FractureSetInitialPressure(fracture,initial_cell_pressure)
  !
  ! Sets the pressure referenced in fracture
  !
  use Material_Aux_class

  implicit none
  
  type(fracture_auxvar_type) :: fracture
  PetscReal, intent(in) :: initial_cell_pressure
  
  fracture%initial_pressure = initial_cell_pressure
!  fracture%properties(frac_init_pres_index) = &
!    fracture%properties(frac_init_pres_index) + initial_cell_pressure
!  fracture%properties(frac_alt_pres_index) = &
!    fracture%properties(frac_alt_pres_index) + &
!    fracture%properties(frac_init_pres_index)

end subroutine FractureSetInitialPressure

! ************************************************************************** !

subroutine FracturePoroEvaluate(auxvar,pressure,compressed_porosity, &
                                dcompressed_porosity_dp)
  !
  ! Calculates porosity induced by fracture BRAGFLO_6.02_UM Eq. (136)
  ! 4.10 Pressure-Induced Fracture Treatment
  !
  ! Author: Heeho Park
  ! Date: 03/12/15
  !

  use Option_module
  use Material_Aux_class
  
  implicit none
  
  type(option_type) :: option
  
!  class(fracture_type) :: this
  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: Ci, Ca
  PetscReal :: P0, Pa, Pi
  PetscReal :: phia, phi0

  ! if fracture is off, still have to calculate soil compressiblity, if 
  ! soil_compressibility_index > 0
  if (.not.auxvar%fracture%fracture_is_on) then
    ! soil_compressibility_index is a file global in material_aux.F90
    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(auxvar,pressure, compressed_porosity, &
                                dcompressed_porosity_dp)
    endif
    return
  endif


  Ci = auxvar%soil_properties(soil_compressibility_index)
  if (associated(MaterialCompressSoilPtr,MaterialCompressSoilBRAGFLO)) then
    ! convert bulk compressibility to pore compressibility
    Ci = auxvar%soil_properties(soil_compressibility_index) / &
         auxvar%porosity_base
  endif
!  P0 = auxvar%soil_properties(soil_reference_pressure_index)
  P0 = auxvar%fracture%initial_pressure
  Pi = auxvar%fracture%properties(frac_init_pres_index) + P0
  Pa = auxvar%fracture%properties(frac_alt_pres_index) + Pi
  phia = auxvar%fracture%properties(frac_max_poro_index)
  phi0 = auxvar%porosity_base
  
  if (P0 < -998.d0) then ! not yet initialized
    compressed_porosity = phi0
    return
  endif
  
  if (pressure < Pi) then
!    call MaterialCompressSoil(auxvar,pressure, compressed_porosity, &
!                              dcompressed_porosity_dp)
    compressed_porosity = phi0 * exp(Ci*(pressure-P0))
  else if (pressure < Pa) then
    Ca = Ci*(1.d0 - 2.d0 * (Pa-P0)/(Pa-Pi)) + &
      2.d0/(Pa-Pi)*log(phia/phi0)
    compressed_porosity = phi0 * exp(Ci*(pressure-P0) + &
      ((Ca-Ci)*(pressure-Pi)**2.d0)/(2.d0*(Pa-Pi)))
    compressed_porosity=min(compressed_porosity, phia)
    !mathematica solution
    dcompressed_porosity_dp = exp(Ci*(pressure-P0) + &
      ((Ca-Ci)*(pressure-Pi)**2.d0) / (2.d0*(Pa-Pi))) * &
      phi0 * (Ci + ((Ca-Ci)*(pressure-Pi)) / (Pa-Pi))
  else if (pressure >= Pa) then
    compressed_porosity = phia
    dcompressed_porosity_dp = 0.d0
  endif

end subroutine FracturePoroEvaluate

! ************************************************************************** !
                                
subroutine FracturePermScale(auxvar,liquid_pressure,effective_porosity, &
                             scaling_factor)
  !
  ! Calculates permeability induced by fracture BRAGFLO_6.02_UM Eq. (136)
  ! 4.10 Pressure-Induced Fracture Treatment
  !
  ! Author: Glenn Hammond
  ! Date: 06/05/17
  !

  use Option_module
  use Material_Aux_class
  
  implicit none
  
  type(option_type) :: option
  
!  class(fracture_type) :: this
  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: liquid_pressure
  PetscReal, intent(in) :: effective_porosity
  PetscReal, intent(out) :: scaling_factor

  PetscReal :: n
  PetscReal :: phi
  PetscReal :: phi0
  PetscReal :: phii
  PetscReal :: phia
  PetscReal :: Ci ! pore compressibility
  PetscReal :: P0 ! initial pressure
  PetscReal :: Pi ! initiating pressure

  if (.not.auxvar%fracture%fracture_is_on) then
    scaling_factor = 1.d0
    return
  endif

  Ci = auxvar%soil_properties(soil_compressibility_index)
  if (associated(MaterialCompressSoilPtr,MaterialCompressSoilBRAGFLO)) then
    ! convert bulk compressibility to pore compressibility
    Ci = auxvar%soil_properties(soil_compressibility_index) / &
         auxvar%porosity_base
  endif
  P0 = auxvar%fracture%initial_pressure
  Pi = auxvar%fracture%properties(frac_init_pres_index) + P0

  phi = effective_porosity
  phi0 = auxvar%porosity_base

  n = auxvar%fracture%properties(frac_poro_exp_index)

  if (P0 < -998.d0) then ! not yet initialized
    scaling_factor = 1.d0
    return
  endif

  ! phi = altered porosity
  ! phii = porosity at initiating pressure
  ! phia = porosity at fully altered
  if (liquid_pressure < Pi) then
    scaling_factor = 1.d0
  else
    phii = phi0 * exp(Ci*(Pi-P0))
    scaling_factor = (phi/phii)**n
  endif

end subroutine FracturePermScale

! ************************************************************************** !

subroutine FractureDestroy(this)
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(fracture_type), pointer :: this
  
  if (.not.associated(this)) return

  deallocate(this)
  nullify(this)

end subroutine FractureDestroy

end module Fracture_module

! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
  
module Creep_Closure_module
  
  use PFLOTRAN_Constants_module
  use Lookup_Table_module

  implicit none
  
  private

#include "petsc/finclude/petscsys.h"

  type, public :: creep_closure_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: num_times
    PetscInt :: num_values_per_time
    PetscReal :: shutdown_pressure
    PetscReal :: time_closeoff
    PetscReal :: time_datamax
    PetscReal :: porosity_minimum
    class(lookup_table_general_type), pointer :: lookup_table
    
    class(creep_closure_type), pointer :: next
    
  contains
    procedure, public :: Read => CreepClosureRead
    procedure, public :: Evaluate => CreepClosureEvaluate
    procedure, public :: Test => CreepClosureTest
  end type creep_closure_type
  
  type, public :: creep_closure_ptr_type
    class(creep_closure_type), pointer :: ptr
  end type creep_closure_ptr_type

  public :: CreepClosureInit, &
            CreepClosureCreate, &
            CreepClosureDestroy, &
            CreepClosureArrayDestroy, &
            CreepClosureRead, &
            CreepClosureAddToList, &
            CreepClosureConvertListToArray, &
            CreepClosureGetID
  
  contains
  
! ************************************************************************** !

subroutine CreepClosureInit()
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
 
  
end subroutine CreepClosureInit

! ************************************************************************** !

function CreepClosureCreate()
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(creep_closure_type), pointer :: CreepClosureCreate
  
  allocate(CreepClosureCreate)
  CreepClosureCreate%name = ''
  CreepClosureCreate%num_times = UNINITIALIZED_INTEGER
  CreepClosureCreate%num_values_per_time = UNINITIALIZED_INTEGER
  CreepClosureCreate%shutdown_pressure = 5.d7 ! set to BRAGFLO default
  CreepClosureCreate%porosity_minimum = 1.d-2 
  CreepClosureCreate%time_closeoff = 1.d20 ! s
  CreepClosureCreate%time_datamax =  1.d20 ! s
  nullify(CreepClosureCreate%lookup_table)
  nullify(CreepClosureCreate%next)
  
end function CreepClosureCreate

! ************************************************************************** !

subroutine CreepClosureRead(this,input,option)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module
  use Units_module
  
  implicit none
  
  class(creep_closure_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: keyword, word, internal_units
  character(len=MAXSTRINGLENGTH) :: error_string = 'CREEP_CLOSURE'
  type(input_type), pointer :: input2
  PetscInt :: temp_int
  PetscReal :: time_units_conversion

  time_units_conversion = 1.d0
  filename = ''
  input%ierr = 0

  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('FILENAME') 
        call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'FILENAME',error_string)
      case('SHUTDOWN_PRESSURE')
        call InputReadDouble(input,option,this%shutdown_pressure)
        call InputErrorMsg(input,option,'shutdown pressure',error_string)
      case('TIME_CLOSEOFF')
        call InputReadDouble(input,option,this%time_closeoff)
        call InputErrorMsg(input,option,'time closeoff',error_string)
     case default
        call InputKeywordUnrecognized(keyword,'CREEP_CLOSURE',option)
    end select
  enddo
  
  if (len_trim(filename) < 1) then
    option%io_buffer = 'FILENAME must be specified for CREEP_CLOSURE.'
    call printErrMsg(option)
  endif
  
  this%lookup_table => LookupTableCreateGeneral(TWO_INTEGER)
  error_string = 'CREEP_CLOSURE file'
  input2 => InputCreate(IUNIT_TEMP,filename,option)
  input2%ierr = 0
  do
    call InputReadPflotranString(input2,option)
    if (InputError(input2)) exit

    call InputReadWord(input2,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input2,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('NUM_TIMES') 
        call InputReadInt(input2,option,this%num_times)
        call InputErrorMsg(input2,option,'number of times',error_string)
      case('NUM_VALUES_PER_TIME') 
        call InputReadInt(input2,option,this%num_values_per_time)
        call InputErrorMsg(input2,option,'number of pressure',error_string)
      case('TIME_UNITS') 
        internal_units = 'sec'
        call InputReadWord(input2,option,word,PETSC_TRUE) 
        call InputErrorMsg(input2,option,'UNITS','CONDITION')   
        call StringToLower(word)
        time_units_conversion = UnitsConvertToInternal(word, &
                                internal_units,option)
      case('TIME')
        if (Uninitialized(this%num_times) .or. &
            Uninitialized(this%num_values_per_time)) then
          option%io_buffer = 'NUM_TIMES and NUM_VALUES_PER_TIME must be ' // &
            'specified prior to reading the corresponding arrays.'
          call printErrMsg(option)
        endif
        this%lookup_table%dims(1) = this%num_times
        this%lookup_table%dims(2) = this%num_values_per_time
        temp_int = this%num_times*this%num_values_per_time
        allocate(this%lookup_table%axis1%values(this%num_times))
        allocate(this%lookup_table%axis2%values(temp_int))
        allocate(this%lookup_table%data(temp_int))
        string = 'TIME in CREEP_CLOSURE'
        call UtilityReadArray(this%lookup_table%axis1%values, &
                              NEG_ONE_INTEGER,string, &
                              input2,option)
        this%lookup_table%axis1%values = this%lookup_table%axis1%values * &
          time_units_conversion
      case('PRESSURE') 
        string = 'PRESSURE in CREEP_CLOSURE'
        call UtilityReadArray(this%lookup_table%axis2%values, &
                              NEG_ONE_INTEGER, &
                              string,input2,option)
      case('POROSITY') 
        string = 'POROSITY in CREEP_CLOSURE'
        call UtilityReadArray(this%lookup_table%data, &
                              NEG_ONE_INTEGER, &
                              string,input2,option)
     case default
        error_string = trim(error_string) // ': ' // filename
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  enddo
  call InputDestroy(input2)
  
  if (size(this%lookup_table%axis1%values) /= this%num_times) then
    option%io_buffer = 'Number of times does not match NUM_TIMES.'
    call printErrMsg(option)
  endif  
  if (size(this%lookup_table%axis2%values) /= &
    this%num_times*this%num_values_per_time) then
    option%io_buffer = 'Number of pressures does not match NUM_TIMES * ' // &
                       'NUM_VALUES_PER_TIME.'
    call printErrMsg(option)
  endif
  if (size(this%lookup_table%data) /= &
    this%num_times*this%num_values_per_time) then
    option%io_buffer = 'Number of porosities does not match NUM_TIMES * ' // &
                       'NUM_VALUES_PER_TIME.'
    call printErrMsg(option)
  endif

  ! set limits
  this%time_datamax = maxval(this%lookup_table%axis1%values)
  this%porosity_minimum = minval(this%lookup_table%data)

end subroutine CreepClosureRead


! ************************************************************************** !

subroutine CreepClosureAddToList(new_creep_closure,list)
  ! 
  ! Adds an object to linked list
  ! 

  implicit none
  
  class(creep_closure_type), pointer :: new_creep_closure
  class(creep_closure_type), pointer :: list

  class(creep_closure_type), pointer :: cur_creep_closure
  
  if (associated(list)) then
    cur_creep_closure => list
    ! loop to end of list
    do
      if (.not.associated(cur_creep_closure%next)) exit
      cur_creep_closure => cur_creep_closure%next
    enddo
    cur_creep_closure%next => new_creep_closure
  else
    list => new_creep_closure
  endif
  
end subroutine CreepClosureAddToList

! ************************************************************************** !

subroutine CreepClosureConvertListToArray(list,array,option)
  ! 
  ! Creates an array of pointers to the objects in the list
  ! 

  use String_module
  use Option_module
  
  implicit none
  
  class(creep_closure_type), pointer :: list
  class(creep_closure_ptr_type), pointer :: array(:)
  type(option_type) :: option
    
  class(creep_closure_type), pointer :: cur_creep_closure
  PetscInt :: count
  
  ! Start at 2
  ! The first element will be a null pointer (no creep closure)  
  count = 1
  cur_creep_closure => list
  do 
    if (.not.associated(cur_creep_closure)) exit
    count = count + 1
    cur_creep_closure => cur_creep_closure%next
  enddo
  
  if (associated(array)) deallocate(array)
  allocate(array(count))

  ! Start at 2
  ! The first element is a null pointer (no creep closure)  
  cur_creep_closure => list
  count = 1
  nullify(array(count)%ptr)
  do 
    if (.not.associated(cur_creep_closure)) exit
    count = count + 1
    array(count)%ptr => cur_creep_closure
    !if (cur_creep_closure%test .and. &
    !    option%myrank == option%io_rank) then
    !  call CreepClosureTest(cur_creep_closure,option)
    !endif
    cur_creep_closure => cur_creep_closure%next
  enddo

end subroutine CreepClosureConvertListToArray

! ************************************************************************** !

function CreepClosureGetID(creep_closure_array, &
                                   creep_closure_name, &
                                   material_property_name, option)
  ! 
  ! Returns the ID of the creep_closure object named
  ! "creep_closure_name"
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  use Option_module
  use String_module
  
  class(creep_closure_ptr_type), pointer :: &
    creep_closure_array(:)
  character(len=MAXWORDLENGTH) :: creep_closure_name
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: CreepClosureGetID

  ! CreepClosureGetID = 1 is a null pointer (no creep closure)
  CreepClosureGetID = 1

  if (len_trim(creep_closure_name) < 1) return
  
  do CreepClosureGetID = 2, size(creep_closure_array)
    if (StringCompare(creep_closure_name, &
                      creep_closure_array( &
                        CreepClosureGetID)%ptr%name)) then
      return
    endif
  enddo
  option%io_buffer = 'Creep closure "' // &
           trim(creep_closure_name) // &
           '" in material property "' // &
           trim(material_property_name) // &
           '" not found among available creep closure tables.'
  call printErrMsg(option)    

end function CreepClosureGetID

! ************************************************************************** !

function CreepClosureEvaluate(this,time,pressure)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(creep_closure_type) :: this
  PetscReal :: time
  PetscReal :: pressure
  
  PetscReal :: CreepClosureEvaluate
  
  CreepClosureEvaluate = this%lookup_table%Sample(time,pressure)
  
end function CreepClosureEvaluate


! ************************************************************************** !

subroutine CreepClosureTest(this,time,pressure)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(creep_closure_type) :: this
  PetscReal :: time
  PetscReal :: pressure
  
  print *, time, pressure, this%Evaluate(time,pressure)
  
end subroutine CreepClosureTest


! ************************************************************************** !

subroutine CreepClosureDestroy(creep_closure)
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(creep_closure_type), pointer :: creep_closure
  
  if (.not.associated(creep_closure)) return

  call LookupTableDestroy(creep_closure%lookup_table)
  deallocate(creep_closure)
  nullify(creep_closure)

end subroutine CreepClosureDestroy


! ************************************************************************** !

subroutine CreepClosureArrayDestroy(creep_closure_array)
  ! 
  ! Destroys an array of pointers
  ! 
  
  implicit none
  
  class(creep_closure_ptr_type), pointer :: &
    creep_closure_array(:)
  
  class(creep_closure_type), pointer :: cur_creep_closure
  PetscInt :: i
  
  if (.not. associated(creep_closure_array)) return

  ! The first element is a null pointer (no creep closure)
  do i=1,size(creep_closure_array)
    call CreepClosureDestroy(creep_closure_array(i)%ptr)
    nullify(creep_closure_array(i)%ptr)
  enddo
  
  deallocate(creep_closure_array)
  nullify(creep_closure_array)

end subroutine CreepClosureArrayDestroy

end module Creep_Closure_module

! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
  
module Klinkenberg_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private

#include "petsc/finclude/petscsys.h"

  type, public :: klinkenberg_type
    PetscReal :: a
    PetscReal :: b
  contains
    procedure, public :: Read => KlinkenbergRead
    procedure, public :: Evaluate => KlinkenbergEvaluate
    procedure, public :: Scale => KlinkenbergScale
    procedure, public :: Test => KlinkenbergTest
  end type klinkenberg_type
  
  class(klinkenberg_type), pointer, public :: klinkenberg
  
  interface KlinkenbergDestroy
    module procedure KlinkenbergDestroy1
    module procedure KlinkenbergDestroy2
  end interface

  public :: KlinkenbergInit, &
            KlinkenbergCreate, &
            KlinkenbergDestroy
  
  contains
  
! ************************************************************************** !

subroutine KlinkenbergInit()
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  if (associated(klinkenberg)) then
    call KlinkenbergDestroy(klinkenberg)
  endif
  nullify(klinkenberg)
  
end subroutine KlinkenbergInit

! ************************************************************************** !

function KlinkenbergCreate()
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(klinkenberg_type), pointer :: KlinkenbergCreate
  
  allocate(KlinkenbergCreate)
  KlinkenbergCreate%a = UNINITIALIZED_DOUBLE
  KlinkenbergCreate%b = UNINITIALIZED_DOUBLE
  
end function KlinkenbergCreate

! ************************************************************************** !

subroutine KlinkenbergRead(this,input,option)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module
  use Units_module
  
  implicit none
  
  class(klinkenberg_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string = 'KLINKENBERG_EFFECT'

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('A') 
        call InputReadDouble(input,option,this%a)
        call InputErrorMsg(input,option,'a',error_string)
      case('B') 
        call InputReadDouble(input,option,this%b)
        call InputErrorMsg(input,option,'b',error_string)
     case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  enddo
  
  if (Uninitialized(this%a)) then
    option%io_buffer = &
      'Parameter "a" must be included to model the Klinkenberg Effect.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%b)) then
    option%io_buffer = &
      'Parameter "b" must be included to model the Klinkenberg Effect.'
    call printErrMsg(option)
  endif
  
end subroutine KlinkenbergRead

! ************************************************************************** !

function KlinkenbergEvaluate(this,liquid_permeability,gas_pressure)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(klinkenberg_type) :: this
  PetscReal :: liquid_permeability(3)
  PetscReal :: gas_pressure
  
  PetscReal :: gas_permeability(3)
  PetscReal :: perm_scale(3)

  PetscReal :: KlinkenbergEvaluate(3)

  call this%Scale(liquid_permeability,gas_pressure,perm_scale)
  gas_permeability(:) = liquid_permeability(:) * perm_scale(:)
  KlinkenbergEvaluate(:) = gas_permeability(:)
  
end function KlinkenbergEvaluate

! ************************************************************************** !

subroutine KlinkenbergScale(this,liquid_permeability,gas_pressure, &
                          permeability_scale)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(klinkenberg_type) :: this
  PetscReal :: liquid_permeability(3) ! [m^2]
  PetscReal :: gas_pressure ! [Pa]
  PetscReal :: permeability_scale(3) ! [m^2]
  
  PetscInt :: i
  
  do i = 1, 3
    permeability_scale(i) = (1.d0 + &
                    this%b * (liquid_permeability(i)**this%a) / gas_pressure)
  enddo
  
end subroutine KlinkenbergScale

! ************************************************************************** !

subroutine KlinkenbergTest(this,liquid_permeability,gas_pressure)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(klinkenberg_type) :: this
  PetscReal :: liquid_permeability(3)
  PetscReal :: gas_pressure
  
  print *, liquid_permeability, gas_pressure, &
           this%Evaluate(liquid_permeability,gas_pressure)
  
end subroutine KlinkenbergTest

! ************************************************************************** !

subroutine KlinkenbergDestroy1()
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  call KlinkenbergDestroy(klinkenberg)

end subroutine KlinkenbergDestroy1

! ************************************************************************** !

subroutine KlinkenbergDestroy2(klinkenberg)
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(klinkenberg_type), pointer :: klinkenberg
  
  if (.not.associated(klinkenberg)) return

  deallocate(klinkenberg)
  nullify(klinkenberg)

end subroutine KlinkenbergDestroy2

end module Klinkenberg_module

! ************************************************************************** !

module WIPP_module
  
  use PFLOTRAN_Constants_module
  use Creep_Closure_module

  implicit none
  
  private

#include "petsc/finclude/petscsys.h"

  type :: wipp_type
    class(creep_closure_type), pointer :: creep_closure_tables
    class(creep_closure_ptr_type), pointer :: creep_closure_tables_array(:)
  end type wipp_type
  
  type(wipp_type), pointer, public :: wipp
  
  interface WIPPDestroy
    module procedure WIPPDestroy1
    module procedure WIPPDestroy2
  end interface
  
  public :: WIPPInit, &
            WIPPGetPtr, &
            WIPPRead, &
            WIPPDestroy

contains


! ************************************************************************** !

subroutine WIPPInit()
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/15
  !

  implicit none
  
  if (associated(wipp)) then
    call WIPPDestroy(wipp)
  endif
  nullify(wipp)  
  
end subroutine WIPPInit

! ************************************************************************** !

function WIPPGetPtr()
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/15
  !

  implicit none
  
  type(wipp_type), pointer :: WIPPGetPtr

  if (.not.associated(wipp)) then
    allocate(wipp)
    nullify(wipp%creep_closure_tables)
    nullify(wipp%creep_closure_tables_array)
  endif
  
  WIPPGetPtr => wipp
  
end function WIPPGetPtr

! ************************************************************************** !

subroutine WIPPRead(input,option)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Input_Aux_module
  use String_module
  use Creep_Closure_module
  
  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option
  
  type(wipp_type), pointer :: wipp
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string = 'WIPP'

  wipp => WIPPGetPtr()
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
!      case('CREEP_CLOSURE')
!        call CreepClosureInit()
!        creep_closure => CreepClosureCreate()
!        call creep_closure%Read(input,option)
!        option%flow%transient_porosity = PETSC_TRUE
!        wipp%creep_closure => creep_closure      
     case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  enddo
  
end subroutine WIPPRead

! ************************************************************************** !

subroutine WIPPDestroy1()
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/15
  !

  implicit none
  
  call WIPPDestroy(wipp)

end subroutine WIPPDestroy1

! ************************************************************************** !

subroutine WippDestroy2(wipp)
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/15
  !

  implicit none
  
  type(wipp_type), pointer :: wipp
  
  if (.not.associated(wipp)) return

  call CreepClosureArrayDestroy(wipp%creep_closure_tables_array)
  nullify(wipp%creep_closure_tables)
  nullify(wipp%creep_closure_tables_array)
  
  deallocate(wipp)
  nullify(wipp)

end subroutine WippDestroy2

end module WIPP_module
