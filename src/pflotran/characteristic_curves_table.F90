module Characteristic_Curves_Table_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Lookup_Table_module

  implicit none

  private
  !cc tables variable dictionary
  PetscInt, parameter, public :: CCT_SAT_WAT = 1
  PetscInt, parameter, public :: CCT_SAT_GAS = 2
  PetscInt, parameter, public :: CCT_SAT_OIL = 3
  PetscInt, parameter, public :: CCT_PCXW = 4
  PetscInt, parameter, public :: CCT_PCOG = 5
  PetscInt, parameter, public :: CCT_KRW = 6
  PetscInt, parameter, public :: CCT_KRG = 7
  PetscInt, parameter, public :: CCT_KROW = 8
  PetscInt, parameter, public :: CCT_KROG = 9
  !add here other characteristic curve table (CCT) variables, 
  !then increase CCT_MAX_VAR_NUM
  PetscInt, parameter, public :: CCT_MAX_VAR_NUM = 9
  
  !cc tables tepy
  PetscInt, parameter, public :: CCT_SWFN = 1
  PetscInt, parameter, public :: CCT_SGFN = 2
  PetscInt, parameter, public :: CCT_SOF2 = 3
  PetscInt, parameter, public :: CCT_SOF3 = 4
  
  !PO This can be defined also as extension of lookup_table_general_type
  ! could save some code, shorten the call signatures via lookup_table,
  ! and reduce the methods to access lookup_table memebers
  type, public :: char_curves_table_type
    character(len=MAXWORDLENGTH) :: name 
    PetscInt :: itype
    PetscInt :: num_fields
    PetscInt :: n_indices !number of indices to be saved for lookup
    PetscInt :: first_index !location of first index in auxvars
    PetscReal :: Swco,Swcr,Sgco,Sgcr,Sowcr,Sogcr !table end points
    PetscBool :: pc_inverse_available
    class(lookup_table_general_type), pointer :: lookup_table
    class(lookup_table_general_type), pointer :: pc_inv_lookup_table
    class(char_curves_table_type), pointer :: next
  contains
    procedure :: SetSWFNTable
    procedure :: SetSGFNTable
    procedure :: SetSOF3Table
    procedure :: SetSOF2Table
    procedure :: CCTVarIsMonotonic
    procedure, public :: CheckCCTVariableExists
    procedure, public :: CharCurveTableVarGrad
    procedure, public :: CharCurvePcInvTableVarGrad
  end type char_curves_table_type

  public :: CharCurvesTableCreate, &
            CharCurvesTableRead, &
            CharCurvesTableAddToList, &
            CharCurveTableGetPtrFromList, &
            SearchCCTVarInCCTableList, &
            CharCurvesTableDestroy
  
contains

! ************************************************************************** !

function CharCurvesTableCreate()
  ! 
  ! Creates a characteristic curve 
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/13/18
  ! 

  implicit none

  class(char_curves_table_type), pointer :: CharCurvesTableCreate
  
  class(char_curves_table_type), pointer :: CharCurvesTable
  
  allocate(CharCurvesTable)
  CharCurvesTable%name = ''
  CharCurvesTable%itype = UNINITIALIZED_INTEGER
  CharCurvesTable%num_fields = UNINITIALIZED_INTEGER
  CharCurvesTable%n_indices = UNINITIALIZED_INTEGER
  CharCurvesTable%first_index = UNINITIALIZED_INTEGER
  CharCurvesTable%Swco = UNINITIALIZED_DOUBLE
  CharCurvesTable%Swcr = UNINITIALIZED_DOUBLE
  CharCurvesTable%Sgco = UNINITIALIZED_DOUBLE
  CharCurvesTable%Sgcr = UNINITIALIZED_DOUBLE
  CharCurvesTable%Sowcr = UNINITIALIZED_DOUBLE
  CharCurvesTable%Sogcr = UNINITIALIZED_DOUBLE
  CharCurvesTable%pc_inverse_available = PETSC_FALSE
  nullify(CharCurvesTable%lookup_table)
  nullify(CharCurvesTable%pc_inv_lookup_table)
  
  !create gen one dim lookup tabel 
  CharCurvesTable%lookup_table => LookupTableCreateGeneral(ONE_INTEGER)
  call CharCurvesTable%lookup_table%LookupTableVarsInit(CCT_MAX_VAR_NUM)

  nullify(CharCurvesTable%next)

  CharCurvesTableCreate => CharCurvesTable

end function CharCurvesTableCreate


! ************************************************************************** !
  

subroutine CharCurvesTableRead(this,input,option)
  ! 
  ! Reads in contents of a characteristic curve table card
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/13/18
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(char_curves_table_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: press_unit
  character(len=MAXSTRINGLENGTH) :: error_string

  PetscReal :: krw,krg,krow,krog
  PetscInt :: s_idx
  PetscBool :: press_unit_found
  PetscBool :: table_found
  PetscReal, parameter :: kr_eps = 1.0d-20
  PetscBool, allocatable :: AxisIsSMInc(:)

  input%ierr = 0
  error_string = 'CHARACTERISTIC_CURVES_TABLE,' // trim(this%name)
  press_unit_found = PETSC_FALSE
  table_found = PETSC_FALSE
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    !-----------------------------------------------------------------------
      case('PRESSURE_UNITS')
        call InputReadWord(input,option,press_unit,PETSC_TRUE)
        call InputErrorMsg(input,option,'PRESSURE_UNITS',error_string)
        press_unit_found = PETSC_TRUE
      case('SWFN')
        this%itype = CCT_SWFN
        call this%SetSWFNTable(option)
        error_string = trim(error_string) // ', SWFN'
        call this%lookup_table%VarDataRead(input,this%num_fields,TWO_INTEGER, &
                                                        error_string,option)
        table_found = PETSC_TRUE
      case('SGFN')
        this%itype = CCT_SGFN
        call this%SetSGFNTable(option)
        error_string = trim(error_string) // ', SGFN'
        call this%lookup_table%VarDataRead(input,this%num_fields,TWO_INTEGER, &
                                                        error_string,option)
        table_found = PETSC_TRUE
      case('SOF2')
        this%itype = CCT_SOF2
        call this%SetSOF2Table(option)
        error_string = trim(error_string) // ', SOF2'
        call this%lookup_table%VarDataRead(input,this%num_fields,TWO_INTEGER, &
                                                        error_string,option)
        table_found = PETSC_TRUE
      case('SOF3')
        this%itype = CCT_SOF3
        call this%SetSOF3Table(option)
        error_string = trim(error_string) // ', SOF3'
        call this%lookup_table%VarDataRead(input,this%num_fields,TWO_INTEGER, &
                                                        error_string,option)
        table_found = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,error_string,option)  
   end select        
 end do
 call InputPopBlock(input,option)

 if ( .not. table_found ) then
   option%io_buffer = trim(error_string) // ', data block not found.'
   call PrintErrMsg(option)
 end if

 select case (this%itype)
   case(CCT_SWFN)
     if ( .not.press_unit_found ) then
       option%io_buffer = trim(error_string) //  &
                                ', SWFN table pressure unit not defined.'
       call PrintErrMsg(option)
     end if
     this%lookup_table%var_array(CCT_PCXW)%ptr%user_units = trim(press_unit)
   case(CCT_SGFN)
     if ( .not.press_unit_found ) then
       option%io_buffer = trim(error_string) //  &
                          ', SGFN table pressure unit not defined.'
       call PrintErrMsg(option)
     end if
     this%lookup_table%var_array(CCT_PCOG)%ptr%user_units = trim(press_unit)
 end select 

 call this%lookup_table%LookupTableVarConvFactors(option)
 call this%lookup_table%VarPointAndUnitConv(option)
 call this%lookup_table%SetupConstValExtrap(option)
 !allocate lookup var gradients 
 call this%lookup_table%LookupTableVarInitGradients(option)

 allocate(this%lookup_table%axis1)
 allocate(this%lookup_table%axis1%values(size( &
                                          this%lookup_table%var_data(1,:))) )
 this%lookup_table%axis1%values = this%lookup_table%var_data(1,:)
 this%lookup_table%dims(1) = size(this%lookup_table%axis1%values(:))
 
 !check saturation is monotonically growing
 allocate(AxisIsSMInc(this%lookup_table%dim))
 call this%lookup_table%AxesAreSMInc(AxisIsSMInc)

 if ( .not. AxisIsSMInc(ONE_INTEGER) ) then
   option%io_buffer = trim(error_string) //  &
                                   ', saturation not monotonically increasing.'
   call PrintErrMsg(option)
 end if

 deallocate(AxisIsSMInc)
 
 !load table end points
 select case (this%itype)
   case(CCT_SWFN)
     this%Swco = this%lookup_table%axis1%values(1)
     do s_idx =1, size(this%lookup_table%axis1%values(:))
       krw = this%lookup_table%var_array(CCT_KRW)%ptr%data(s_idx)
       if (krw < kr_eps ) then
         this%Swcr = this%lookup_table%axis1%values(s_idx)
       else
         exit
       end if   
     end do
     !create inverse table for capillary pressures Pcow
     if (this%CCTVarIsMonotonic(CCT_PCXW,option)) then
       this%pc_inverse_available = PETSC_TRUE
       this%pc_inv_lookup_table => &
           InverseLookupTableCreateGen(this%lookup_table,CCT_PCXW,option)
     end if     
   case(CCT_SGFN)
     this%Sgco = this%lookup_table%axis1%values(1)
     do s_idx =1, size(this%lookup_table%axis1%values(:))
       krg = this%lookup_table%var_array(CCT_KRG)%ptr%data(s_idx)
       if (krg < kr_eps ) then
         this%Sgcr = this%lookup_table%axis1%values(s_idx) 
       else
         exit
       end if       
     end do
     !create inverse table for capillary pressures Pcog
     if (this%CCTVarIsMonotonic(CCT_PCOG,option)) then
       this%pc_inverse_available = PETSC_TRUE
       this%pc_inv_lookup_table => &
          InverseLookupTableCreateGen(this%lookup_table,CCT_PCOG,option)
     end if    
   case(CCT_SOF2)
     do s_idx =1, size(this%lookup_table%axis1%values(:))
       krow = this%lookup_table%var_array(CCT_KROW)%ptr%data(s_idx)
       if (krow < kr_eps ) then
         this%Sowcr = this%lookup_table%axis1%values(s_idx)
       else
         exit
       end if
     end do      
   case(CCT_SOF3)
     do s_idx =1, size(this%lookup_table%axis1%values(:))
        krow = this%lookup_table%var_array(CCT_KROW)%ptr%data(s_idx)
        krog = this%lookup_table%var_array(CCT_KROG)%ptr%data(s_idx)
        if (krow < kr_eps ) then
          this%Sowcr = this%lookup_table%axis1%values(s_idx)
        end if
        if (krog < kr_eps ) then
          this%Sogcr = this%lookup_table%axis1%values(s_idx)
        end if        
      end do
 end select
 
 if ( associated(this%pc_inv_lookup_table) ) then
   this%n_indices = 2 
 else   
   this%n_indices = 1
 end if

end subroutine CharCurvesTableRead  

! ************************************************************************** !

subroutine SetSWFNTable(this,option)
  ! 
  ! Setup SWFN table
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/14/18
  ! 

  use Option_module
  use String_module

  implicit none

  class(char_curves_table_type) :: this
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: cct_var => null()
  character(len=MAXWORDLENGTH) :: internal_units, user_units
  PetscInt :: data_idx

  this%num_fields = 3

  !default units to none - only pressure for Pc is required
  internal_units = 'unitless' 
  user_units = 'unitless'
  
  data_idx = 1 !position in the table
  cct_var => CreateLookupTableVar(CCT_SAT_WAT,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var) 

  data_idx = 2 !position in the table
  cct_var => CreateLookupTableVar(CCT_KRW,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var)

  internal_units = 'Pa'
  user_units = 'unitless' !to be read from input
  data_idx = 3 !position in the table
  cct_var => CreateLookupTableVar(CCT_PCXW,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var)

end subroutine SetSWFNTable

! ************************************************************************** !

subroutine SetSGFNTable(this,option)
  ! 
  ! Setup SGFN table
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/17/18
  ! 

  use Option_module
  use String_module

  implicit none

  class(char_curves_table_type) :: this
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: cct_var => null()
  character(len=MAXWORDLENGTH) :: internal_units, user_units
  PetscInt :: data_idx

  this%num_fields = 3

  !default units to none - only pressure for Pc is required
  internal_units = 'unitless' 
  user_units = 'unitless'
  
  data_idx = 1 !position in the table
  cct_var => CreateLookupTableVar(CCT_SAT_GAS,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var) 

  data_idx = 2 !position in the table
  cct_var => CreateLookupTableVar(CCT_KRG,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var)

  internal_units = 'Pa'
  user_units = 'unitless' !to be read from input
  data_idx = 3 !position in the table
  cct_var => CreateLookupTableVar(CCT_PCOG,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var)

end subroutine SetSGFNTable

! ************************************************************************** !

subroutine SetSOF3Table(this,option)
  ! 
  ! Setup SOF3 table
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/17/18
  ! 

  use Option_module
  use String_module

  implicit none

  class(char_curves_table_type) :: this
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: cct_var => null()
  character(len=MAXWORDLENGTH) :: internal_units, user_units
  PetscInt :: data_idx

  this%num_fields = 3

  !default units to none - only pressure for Pc is required
  internal_units = 'unitless' 
  user_units = 'unitless'
  
  data_idx = 1 !position in the table
  cct_var => CreateLookupTableVar(CCT_SAT_OIL,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var) 

  data_idx = 2 !position in the table
  cct_var => CreateLookupTableVar(CCT_KROW,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var)

  data_idx = 3 !position in the table
  cct_var => CreateLookupTableVar(CCT_KROG,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var)

end subroutine SetSOF3Table

! ************************************************************************** !

subroutine SetSOF2Table(this,option)
  ! 
  ! Setup SOF2 table
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/17/18
  ! 

  use Option_module
  use String_module

  implicit none

  class(char_curves_table_type) :: this
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: cct_var => null()
  character(len=MAXWORDLENGTH) :: internal_units, user_units
  PetscInt :: data_idx

  this%num_fields = 2

  !default units to none - only pressure for Pc is required
  internal_units = 'unitless' 
  user_units = 'unitless'
  
  data_idx = 1 !position in the table
  cct_var => CreateLookupTableVar(CCT_SAT_OIL,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var) 

  data_idx = 2 !position in the table
  cct_var => CreateLookupTableVar(CCT_KROW,internal_units, &
                                                         user_units,data_idx)
  call LookupTableVarAddToList(cct_var,this%lookup_table%vars)
  this%lookup_table%var_array(cct_var%iname)%ptr => cct_var
  nullify(cct_var)

end subroutine SetSOF2Table

! ************************************************************************** !

function CCTVarIsMonotonic(this,var_iname,option)
  ! 
  ! Check if a CCT var (int name = var_iname) is monotonic or not
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/20/18
  ! 
  use Option_module

  implicit none

  PetscBool :: CCTVarIsMonotonic
  class(char_curves_table_type) :: this
  PetscInt, intent(in) :: var_iname
  type(option_type) :: option

  PetscReal, parameter :: diff_eps = 1.0d-6
  PetscInt :: data_idx, var_idx_max
  PetscInt :: i_tmp
  PetscReal :: data_tmp
  
  data_idx = this%lookup_table%var_array(var_iname)%ptr%data_idx
  var_idx_max = size(this%lookup_table%axis1%values(:))
  
  CCTVarIsMonotonic = PETSC_TRUE
  
  !check if a constant value has been given
  if ( dabs(this%lookup_table%var_data(data_idx,1) -  &
            this%lookup_table%var_data(data_idx,var_idx_max)  &
           )  < diff_eps &
     ) then
    CCTVarIsMonotonic = PETSC_FALSE
    return
  end if
  
  if (this%lookup_table%var_data(data_idx,1) > &
      this%lookup_table%var_data(data_idx,var_idx_max)) then
    do i_tmp = 1,var_idx_max
      data_tmp = this%lookup_table%var_data(data_idx,i_tmp)
      if ( data_tmp > this%lookup_table%var_data(data_idx,1) ) then
        CCTVarIsMonotonic = PETSC_FALSE
      end if
    end do  
  else if (this%lookup_table%var_data(data_idx,1) < &
        this%lookup_table%var_data(data_idx,var_idx_max)) then
    do i_tmp = 1,var_idx_max
      data_tmp = this%lookup_table%var_data(data_idx,i_tmp)
      if ( data_tmp < this%lookup_table%var_data(data_idx,1) ) then
        CCTVarIsMonotonic = PETSC_FALSE
      end if
    end do
  end if

end function CCTVarIsMonotonic

! ************************************************************************** !

subroutine CharCurvesTableAddToList(new_char_curves_table,list)
  !
  ! Adds a characteristic curves table object to the end of linked list
  !
  ! Author: Paolo Orsini
  ! Date: 08/13/18
  !

  implicit none

  class(char_curves_table_type), pointer :: new_char_curves_table
  class(char_curves_table_type), pointer :: list

  class(char_curves_table_type), pointer :: cur_char_curves_table

  if (associated(list)) then
    cur_char_curves_table => list
    ! loop to end of list
    do
      if (.not.associated(cur_char_curves_table%next)) exit
      cur_char_curves_table => cur_char_curves_table%next
    enddo
    cur_char_curves_table%next => new_char_curves_table
  else
    list => new_char_curves_table
  endif

end subroutine CharCurvesTableAddToList

! ************************************************************************** !

function CharCurveTableGetPtrFromList(cc_table_name,list,error_string,option)
  ! 
  ! Returns a pointer to the cc_table matching cc_table_name
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/16/18
  ! 

  use String_module
  use Option_module

  implicit none
  
  class(char_curves_table_type), pointer :: CharCurveTableGetPtrFromList
  character(len=MAXWORDLENGTH) :: cc_table_name
  class(char_curves_table_type), pointer :: list
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option  
  
  class(char_curves_table_type), pointer :: cur_char_curve_table
  PetscInt :: length
  
  nullify(CharCurveTableGetPtrFromList)
  !point the first object in the list 
  cur_char_curve_table => list
  do 
    if (.not.associated(cur_char_curve_table)) exit !retrun null pointer
    length = len_trim(cc_table_name)
    if (length == len_trim(cur_char_curve_table%name) .and. &
        StringCompare(cur_char_curve_table%name,cc_table_name,length)) then
      CharCurveTableGetPtrFromList => cur_char_curve_table
      return
    endif
    cur_char_curve_table => cur_char_curve_table%next
  enddo

  if (.not. associated(CharCurveTableGetPtrFromList) ) then
    error_string = trim(error_string) // "table name =" // &
                    trim(cc_table_name) // " not found"
    option%io_buffer = error_string
    call PrintErrMsg(option)
  end if
    
end function CharCurveTableGetPtrFromList

! ************************************************************************** !

subroutine SearchCCTVarInCCTableList(list,var_iname,cc_table_name, &
                                                        error_string,option)
  ! 
  ! Returns a pointer to the cc_table matching cc_table_name
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/16/18
  ! 

  use String_module
  use Option_module

  implicit none
  
  class(char_curves_table_type), pointer :: list
  PetscInt, intent(in) :: var_iname
  character(len=MAXWORDLENGTH), intent(out) :: cc_table_name
  character(len=MAXSTRINGLENGTH), intent(in) :: error_string
  type(option_type) :: option  
  
  class(char_curves_table_type), pointer :: cur_char_curves_table
  PetscInt :: num_occurrences
  
  num_occurrences = 0
  cc_table_name = ''
  
  !point the first object in the list 
  cur_char_curves_table => list
  do 
    if (.not.associated(cur_char_curves_table)) exit 
    !occurrence check
    if (cur_char_curves_table%lookup_table%LookupTableVarIsPresent( &
                                                          var_iname)) then
      num_occurrences = num_occurrences + 1
      cc_table_name = trim(cur_char_curves_table%name)
    end if
    cur_char_curves_table => cur_char_curves_table%next
  enddo

  if (num_occurrences == 1) then
    return
  else if (num_occurrences > 1) then
    option%io_buffer = trim(error_string) // &
            ' data found in multiple tables - please select one of the tables'
    call PrintErrMsg(option)
  else if ( num_occurrences < 1) then
    option%io_buffer = trim(error_string) // &
                               'data not found within the tables'
    call PrintErrMsg(option)
  end if
    
end subroutine SearchCCTVarInCCTableList

! ************************************************************************** !

subroutine CheckCCTVariableExists(this,var_iname,error_string,option)
  ! 
  ! Check if a CCT variable is present 
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/18/18
  ! 

  use Option_module
  use String_module
  
  implicit none
  
  class(char_curves_table_type) :: this
  PetscInt, intent(in) :: var_iname
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: error_string_lc

  if (.not.this%lookup_table%LookupTableVarIsPresent(var_iname)) then
    error_string_lc = 'data not found in table = ' // trim(this%name)
    option%io_buffer = error_string_lc
    call PrintErrMsg(option)
  end if

end subroutine CheckCCTVariableExists

! ************************************************************************** !

subroutine CharCurveTableVarGrad(this,sat,var_iname,sample,dSampledSat, &
                                 ierr,indices)
  !
  ! Author: Paolo Orsini
  ! Date: 08/17/18
  !
  ! interpolates a single var of cc table and computes its 
  ! gradient with respect to saturation
  !
  implicit none

  class(char_curves_table_type) :: this
  PetscReal, intent(in) :: sat              ! [-]
  PetscInt, intent(in) :: var_iname
  PetscReal, intent(out) :: sample ! PFLOTRAN intenral units (SI)
  PetscReal, intent(out) :: dSampledSat ! internal PFLOTRAN units ([Sample]/[-])
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer,optional, intent(inout) :: indices(:)

  ierr = 0

  if (present(indices) ) then
    this%lookup_table%axis1%saved_index = indices(this%first_index)
  else
    this%lookup_table%axis1%saved_index = 1
  endif

  call this%lookup_table%SampleAndGradient(var_iname,sat)
  sample = this%lookup_table%var_array(var_iname)%ptr%sample
  dSampledSat = this%lookup_table%var_array(var_iname)%ptr%sample_grad(1)
  !must copy back the index because
  ! saved_index = indices(eos_table%first_index) is a copying operation not pointing
  ! saved_index should be a pointer to save two assigment operations (in/out)
  if (present(indices) ) then
    indices(this%first_index) = this%lookup_table%axis1%saved_index
  endif

end subroutine CharCurveTableVarGrad

! ************************************************************************** !

subroutine CharCurvePcInvTableVarGrad(this,pc,var_iname,sample,dSampledPc, &
                                      ierr,indices)
  !
  ! Author: Paolo Orsini
  ! Date: 08/17/18
  !
  ! interpolates a single var of the Pc cc inverse table and computes its 
  ! gradient wrt Pc
  !
  implicit none

  class(char_curves_table_type) :: this
  PetscReal, intent(in) :: pc              ! [-]
  PetscInt, intent(in) :: var_iname
  PetscReal, intent(out) :: sample ! PFLOTRAN intenral units (SI)
  PetscReal, intent(out) :: dSampledPc ! internal PFLOTRAN units ([Sample]/[Pa])
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer,optional, intent(inout) :: indices(:)

  ierr = 0

  if (present(indices) ) then
    this%pc_inv_lookup_table%axis1%saved_index = indices(this%first_index+1)
  else
    this%pc_inv_lookup_table%axis1%saved_index = 1
  endif

  call this%pc_inv_lookup_table%SampleAndGradient(var_iname,pc)
  sample = this%pc_inv_lookup_table%var_array(var_iname)%ptr%sample
  dSampledPc = this%pc_inv_lookup_table%var_array(var_iname)%ptr%sample_grad(1)
  !must copy back the index because
  ! saved_index = indices(eos_table%first_index) is a copying operation not pointing
  ! saved_index should be a pointer to save two assigment operations (in/out)
  if (present(indices) ) then
    indices(this%first_index+1) = this%pc_inv_lookup_table%axis1%saved_index
  endif

end subroutine CharCurvePcInvTableVarGrad

! ************************************************************************** !


recursive subroutine CharCurvesTableDestroy(cc_table)
  ! 
  ! Destroys a characteristic curve
  ! 
  ! Author: Paolo Orsini
  ! Date: 08/14/18
  ! 

  implicit none
  
  class(char_curves_table_type), pointer :: cc_table
  
  if (.not.associated(cc_table)) return
  
  call CharCurvesTableDestroy(cc_table%next)
  
  call LookupTableDestroy(cc_table%pc_inv_lookup_table)
  call LookupTableDestroy(cc_table%lookup_table)

  nullify(cc_table%next)
  deallocate(cc_table)
  nullify(cc_table)
  
end subroutine CharCurvesTableDestroy

! ************************************************************************** !
  
end module Characteristic_Curves_Table_module
