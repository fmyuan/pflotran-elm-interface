module Characteristic_Curves_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_Common_module
  use Characteristic_Curves_OWG_module
  use Characteristic_Curves_WIPP_module

  implicit none

  private


  type, public :: characteristic_curves_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
    PetscBool :: test
    class(sat_func_base_type), pointer :: saturation_function
    class(sat_func_owg_base_type), pointer :: oil_wat_sat_func
    class(sat_func_owg_base_type), pointer :: oil_gas_sat_func
    class(rel_perm_func_base_type), pointer :: liq_rel_perm_function
    class(rel_perm_func_base_type), pointer :: gas_rel_perm_function
    class(rel_perm_func_base_type), pointer :: oil_rel_perm_function
    class(rel_perm_func_owg_base_type), pointer :: wat_rel_perm_func_owg
    class(rel_perm_func_owg_base_type), pointer :: oil_rel_perm_func_owg
    class(rel_perm_func_owg_base_type), pointer :: gas_rel_perm_func_owg
    class(characteristic_curves_type), pointer :: next
  end type characteristic_curves_type
  
  type, public :: characteristic_curves_ptr_type
    class(characteristic_curves_type), pointer :: ptr
  end type characteristic_curves_ptr_type 
  
  public :: CharacteristicCurvesCreate, &
            CharacteristicCurvesRead, &
            CharacteristicCurvesAddToList, &
            CharCurvesConvertListToArray, &
            CharacteristicCurvesGetID, &
            CharCurvesGetGetResidualSats, &
            CharacteristicCurvesDestroy, &
            CharCurvesInputRecord

contains

! ************************************************************************** !

function CharacteristicCurvesCreate()
  ! 
  ! Creates a characteristic curve object that holds parameters and pointers
  ! to functions for calculating saturation, capillary pressure, relative
  ! permeability, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/23/14
  ! 

  implicit none

  class(characteristic_curves_type), pointer :: CharacteristicCurvesCreate
  
  class(characteristic_curves_type), pointer :: characteristic_curves
  
  allocate(characteristic_curves)
  characteristic_curves%name = ''
  characteristic_curves%print_me = PETSC_FALSE
  characteristic_curves%test = PETSC_FALSE
  nullify(characteristic_curves%saturation_function)
  nullify(characteristic_curves%oil_wat_sat_func)
  nullify(characteristic_curves%oil_gas_sat_func)
  nullify(characteristic_curves%liq_rel_perm_function)
  nullify(characteristic_curves%gas_rel_perm_function)
  nullify(characteristic_curves%oil_rel_perm_function)
  nullify(characteristic_curves%wat_rel_perm_func_owg)
  nullify(characteristic_curves%oil_rel_perm_func_owg)
  nullify(characteristic_curves%gas_rel_perm_func_owg)
  nullify(characteristic_curves%next)

  CharacteristicCurvesCreate => characteristic_curves

end function CharacteristicCurvesCreate

! ************************************************************************** !

subroutine CharacteristicCurvesRead(this,input,option)
  ! 
  ! Reads in contents of a saturation_function card
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(characteristic_curves_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word, phase_keyword
  character(len=MAXWORDLENGTH) :: interface_keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  class(rel_perm_func_base_type), pointer :: rel_perm_function_ptr
  class(sat_func_owg_base_type), pointer :: sat_func_owg_ptr
  class(rel_perm_func_owg_base_type), pointer :: rel_perm_func_owg_ptr

  input%ierr = 0
  error_string = 'CHARACTERISTIC_CURVES'  
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    !-----------------------------------------------------------------------
      case('SATURATION_FUNCTION')
        call InputReadWordDbaseCompatible(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'SATURATION_FUNCTION',error_string)
        call StringToUpper(word)
        select case(word)
          case('CONSTANT')
            this%saturation_function => SF_Constant_Create()
          case('VAN_GENUCHTEN')
            this%saturation_function => SF_VG_Create()
          case('BROOKS_COREY')
            this%saturation_function => SF_BC_Create()
          case('LINEAR')
            this%saturation_function => SF_Linear_Create()
          case('MODIFIED_KOSUGI')
            this%saturation_function => SF_mK_Create()
          case('BRAGFLO_KRP1')
            this%saturation_function => SF_KRP1_Create()
          case('BRAGFLO_KRP2')
            this%saturation_function => SF_KRP2_Create()
          case('BRAGFLO_KRP3')
            this%saturation_function => SF_KRP3_Create()
          case('BRAGFLO_KRP4')
            this%saturation_function => SF_KRP4_Create()
          case('BRAGFLO_KRP5')
            this%saturation_function => SF_KRP5_Create()
          case('BRAGFLO_KRP8')
            this%saturation_function => SF_KRP8_Create()
          case('BRAGFLO_KRP9')
            this%saturation_function => SF_KRP9_Create()
          case('BRAGFLO_KRP11')
            this%saturation_function => SF_KRP11_Create()
          case('BRAGFLO_KRP12')
            this%saturation_function => SF_KRP12_Create()
          case default
            call InputKeywordUnrecognized(word,'SATURATION_FUNCTION',option)
        end select
        call SaturationFunctionRead(this%saturation_function,input,option)
    !-----------------------------------------------------------------------
      case('SATURATION_FUNCTION_OWG')
        call InputReadWordDbaseCompatible(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'SATURATION_FUNCTION_OWG',error_string)
        call StringToUpper(word)
        nullify(sat_func_owg_ptr)
        interface_keyword = 'NONE'
        select case(word)
          case('VAN_GENUCHTEN_OW')
            sat_func_owg_ptr => SF_OW_VG_Create()
            interface_keyword = 'OIL_WATER'
          case('BROOKS_COREY_OG')
            sat_func_owg_ptr => SF_OG_BC_Create()
            interface_keyword = 'OIL_GAS'
          case('VAN_GENUCHTEN_OG_SL')
            sat_func_owg_ptr => SF_OG_VG_SL_Create()
            interface_keyword = 'OIL_GAS'
          case default
            call InputKeywordUnrecognized(word,'SATURATION_FUNCTION_OWG',option)
        end select
        call SaturationFunctionOWGRead(sat_func_owg_ptr,input,option)
        select case(interface_keyword)
          case('OIL_WATER')
            this%oil_wat_sat_func => sat_func_owg_ptr
          case('OIL_GAS')
            this%oil_gas_sat_func => sat_func_owg_ptr
        end select
      case('PERMEABILITY_FUNCTION')
        nullify(rel_perm_function_ptr)
        phase_keyword = 'NONE'
        call InputReadWordDbaseCompatible(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'PERMEABILITY_FUNCTION',error_string)
        call StringToUpper(word)
        select case(word)
          case('MUALEM','MUALEM_VG_LIQ')
            rel_perm_function_ptr => RPF_Mualem_VG_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MUALEM_VG_GAS')
            rel_perm_function_ptr => RPF_Mualem_VG_Gas_Create()
            phase_keyword = 'GAS'
          case('BURDINE','BURDINE_BC_LIQ')
            rel_perm_function_ptr => RPF_Burdine_BC_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BURDINE_BC_GAS')
            rel_perm_function_ptr => RPF_Burdine_BC_Gas_Create()
            phase_keyword = 'GAS'
          case('TOUGH2_IRP7_LIQ')
            rel_perm_function_ptr => RPF_Mualem_VG_Liq_Create()
            phase_keyword = 'LIQUID'
          case('TOUGH2_IRP7_GAS')
            rel_perm_function_ptr => RPF_TOUGH2_IRP7_Gas_Create()
            phase_keyword = 'GAS'
          case('MUALEM_BC_LIQ')
            rel_perm_function_ptr => RPF_Mualem_BC_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MUALEM_BC_GAS')
            rel_perm_function_ptr => RPF_Mualem_BC_Gas_Create()
            phase_keyword = 'GAS'
          case('BURDINE_VG_LIQ')
            rel_perm_function_ptr => RPF_Burdine_VG_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BURDINE_VG_GAS')
            rel_perm_function_ptr => RPF_Burdine_VG_Gas_Create()
            phase_keyword = 'GAS'
          case('MUALEM_LINEAR_LIQ')
            rel_perm_function_ptr => RPF_Mualem_Linear_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MUALEM_LINEAR_GAS')
            rel_perm_function_ptr => RPF_Mualem_Linear_Gas_Create()
            phase_keyword = 'GAS'
          case('BURDINE_LINEAR_LIQ')
            rel_perm_function_ptr => RPF_Burdine_Linear_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BURDINE_LINEAR_GAS')
            rel_perm_function_ptr => RPF_Burdine_Linear_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP1_LIQ')
            rel_perm_function_ptr => RPF_KRP1_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP1_GAS')
            rel_perm_function_ptr => RPF_KRP1_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP2_LIQ')
            rel_perm_function_ptr => RPF_KRP2_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP2_GAS')
            rel_perm_function_ptr => RPF_KRP2_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP3_LIQ')
            rel_perm_function_ptr => RPF_KRP3_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP3_GAS')
            rel_perm_function_ptr => RPF_KRP3_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP4_LIQ')
            rel_perm_function_ptr => RPF_KRP4_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP4_GAS')
            rel_perm_function_ptr => RPF_KRP4_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP5_LIQ')
            rel_perm_function_ptr => RPF_KRP5_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP5_GAS')
            rel_perm_function_ptr => RPF_KRP5_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP8_LIQ')
            rel_perm_function_ptr => RPF_KRP8_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP8_GAS')
            rel_perm_function_ptr => RPF_KRP8_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP9_LIQ')
            rel_perm_function_ptr => RPF_KRP9_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP9_GAS')
            rel_perm_function_ptr => RPF_KRP9_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP11_LIQ')
            rel_perm_function_ptr => RPF_KRP11_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP11_GAS')
            rel_perm_function_ptr => RPF_KRP11_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP12_LIQ')
            rel_perm_function_ptr => RPF_KRP12_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP12_GAS')
            rel_perm_function_ptr => RPF_KRP12_Gas_Create()
            phase_keyword = 'GAS'
          case('MODIFIED_KOSUGI_LIQ')
            rel_perm_function_ptr => RPF_mK_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MODIFIED_KOSUGI_GAS')
            rel_perm_function_ptr => RPF_mK_Gas_Create()
            phase_keyword = 'GAS'
          case('TOUGH2_LINEAR_OIL')
            rel_perm_function_ptr => RPF_TOUGH2_Linear_Oil_Create()
            phase_keyword = 'OIL'
          case('MOD_BC_LIQ')
            rel_perm_function_ptr => RPF_Mod_BC_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MOD_BC_OIL')
            rel_perm_function_ptr => RPF_Mod_BC_Oil_Create()
            phase_keyword = 'OIL'
          case('CONSTANT')
            rel_perm_function_ptr => RPF_Constant_Create()
            ! phase_keyword = 'NONE'
          case default
            call InputKeywordUnrecognized(word,'PERMEABILITY_FUNCTION',option)
        end select
        call PermeabilityFunctionRead(rel_perm_function_ptr,phase_keyword, &
                                      input,option)
        ! if PHASE is specified, have to align correct pointer
        select case(phase_keyword)
          case('GAS')
            this%gas_rel_perm_function => rel_perm_function_ptr
          case('LIQUID')
            this%liq_rel_perm_function => rel_perm_function_ptr
          case('OIL')
            this%oil_rel_perm_function => rel_perm_function_ptr 
            ! PO: gas_rel_perm_fucntion initiated oil_rel_perm_function
            ! to pass the verification in CharacteristicCurvesVerify
            ! in case gas_rel_perm_function is not defined in the input
            ! We should change CharacteristicCurvesVerify instead
            ! this%gas_rel_perm_function => rel_perm_function_ptr
          case('NONE')
            option%io_buffer = 'PHASE has not been set for &
                               &CHARACTERISTIC_CURVES,PERMEABILITY_FUNCTION. &
                               &This is most likely a development issue, and &
                               &not an input deck mistake. Please e-mail &
                               &pflotran-dev@googlegroups.com.' 
            call printErrMsg(option)
          case default
            call InputKeywordUnrecognized(word, &
              'PERMEABILITY_FUNCTION,PHASE',option)
        end select
      case('PERMEABILITY_FUNCTION_OWG')
        nullify(rel_perm_func_owg_ptr)
        phase_keyword = 'NONE'
        call InputReadWordDbaseCompatible(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'PERMEABILITY_FUNCTION_OWG', &
                           error_string)
        call StringToUpper(word)
        select case(word)
        case('MOD_BROOKS_COREY_WATER')
            rel_perm_func_owg_ptr => RPF_wat_MBC_Create()
            phase_keyword = 'WATER'
          case('MOD_BROOKS_COREY_OIL')
            rel_perm_func_owg_ptr => RPF_oil_MBC_Create()
            phase_keyword = 'OIL'
          case('MOD_BROOKS_COREY_HYDROCARBON')
            rel_perm_func_owg_ptr => RPF_oil_MBC_Create()
            rel_perm_func_owg_ptr%So_is_Sh = PETSC_TRUE
            phase_keyword = 'OIL'
          case('MOD_BROOKS_COREY_GAS')
            rel_perm_func_owg_ptr => RPF_gas_MBC_Create()
            phase_keyword = 'GAS'
          case('ECLIPSE_OIL')
            rel_perm_func_owg_ptr => RPF_oil_ecl_Create()
            phase_keyword = 'OIL'
          case('MUALEM_VG_WAT')
            rel_perm_func_owg_ptr => RPF_OWG_Mualem_VG_wat_Create()
            phase_keyword = 'WATER'
          case('MUALEM_VG_GAS_SL')
            rel_perm_func_owg_ptr => RPF_OWG_Mualem_VG_gas_Create()
            rel_perm_func_owg_ptr%function_of_liquid_sat = PETSC_TRUE
            phase_keyword = 'GAS'
          case('TOUGH2_IRP7_GAS_SL')
            rel_perm_func_owg_ptr => RPF_OWG_TOUGH2_IRP7_gas_Create()
            rel_perm_func_owg_ptr%function_of_liquid_sat = PETSC_TRUE
            phase_keyword = 'GAS'
          case('BURDINE_VG_WAT')
            rel_perm_func_owg_ptr => RPF_OWG_Burdine_VG_wat_Create()
            phase_keyword = 'WATER'
          case('BURDINE_VG_GAS_SL')
            rel_perm_func_owg_ptr => RPF_OWG_Burdine_VG_gas_Create()
            rel_perm_func_owg_ptr%function_of_liquid_sat = PETSC_TRUE
            phase_keyword = 'GAS'
          case('BURDINE_BC_WAT')
            rel_perm_func_owg_ptr => RPF_OWG_Burdine_BC_wat_Create()
            phase_keyword = 'WATER'
          case('BURDINE_BC_GAS_SL')
            rel_perm_func_owg_ptr => RPF_OWG_Burdine_BC_gas_Create()
            rel_perm_func_owg_ptr%function_of_liquid_sat = PETSC_TRUE
            phase_keyword = 'GAS'
          case default
            call InputKeywordUnrecognized(word,'PERMEABILITY_FUNCTION_OWG', &
                                                option)
        end select
        call PermeabilityFunctionOWGRead(rel_perm_func_owg_ptr,phase_keyword, &
                                         input,option)
        ! align to correct pointer
        phase_keyword = trim(phase_keyword)
        select case(phase_keyword)
          case('WATER')
            this%wat_rel_perm_func_owg => rel_perm_func_owg_ptr
          case('OIL')
            this%oil_rel_perm_func_owg => rel_perm_func_owg_ptr
          case('GAS')
            this%gas_rel_perm_func_owg => rel_perm_func_owg_ptr
        end select
      case('TEST') 
        this%test = PETSC_TRUE
      case('DEFAULT')
        this%saturation_function => SF_Default_Create()
        this%liq_rel_perm_function => RPF_Default_Create()
        this%gas_rel_perm_function => this%liq_rel_perm_function
        ! PO TODO: adds default for OWG functions
      case default
        call InputKeywordUnrecognized(keyword,'CHARACTERISTIC_CURVES',option)
    end select 
  enddo
  
  select case(option%iflowmode)
    case(TOWG_MODE)
      call CharacteristicCurvesOWGVerify(this,option)
    case default
  call CharacteristicCurvesVerify(this,option)
   end select

end subroutine CharacteristicCurvesRead

! ************************************************************************** !

subroutine SaturationFunctionRead(saturation_function,input,option)
  !
  ! Reads in contents of a SATURATION_FUNCTION block
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(sat_func_base_type) :: saturation_function
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: smooth

  input%ierr = 0
  smooth = PETSC_FALSE
  error_string = 'CHARACTERISTIC_CURVES,SATURATION_FUNCTION,'
  select type(sf => saturation_function)
    class is(sat_func_constant_type)
      error_string = trim(error_string) // 'CONSTANT'
    class is(sat_func_VG_type)
      error_string = trim(error_string) // 'VAN_GENUCHTEN'
    class is(sat_func_BC_type)
      error_string = trim(error_string) // 'BROOKS_COREY'
    class is(sat_func_Linear_type)
      error_string = trim(error_string) // 'LINEAR'
    class is(sat_func_mK_type)
      error_string = trim(error_string) // 'MODIFIED_KOSUGI'
    class is(sat_func_KRP1_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP1'
    class is(sat_func_KRP2_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP2'
    class is(sat_func_KRP3_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP3'
    class is(sat_func_KRP4_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP4'
    class is(sat_func_KRP5_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP5'
    class is(sat_func_KRP8_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP8'
    class is(sat_func_KRP9_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP9'
    class is(sat_func_KRP11_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP11'
    class is(sat_func_KRP12_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP12'
  end select
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   

    ! base 
    found = PETSC_TRUE
    select case(keyword)
      case('LIQUID_RESIDUAL_SATURATION') 
        call InputReadDouble(input,option,saturation_function%Sr)
        call InputErrorMsg(input,option,'LIQUID_RESIDUAL_SATURATION', &
                           error_string)
      case('MAX_CAPILLARY_PRESSURE') 
        call InputReadDouble(input,option,saturation_function%pcmax)
        call InputErrorMsg(input,option,'MAX_CAPILLARY_PRESSURE', &
                            error_string)
      case('SMOOTH')
        smooth = PETSC_TRUE
      case default
        found = PETSC_FALSE
    end select
    
    if (found) cycle
    
    select type(sf => saturation_function)
    !------------------------------------------
      class is(sat_func_constant_type)
        select case(keyword)
          case('CONSTANT_CAPILLARY_PRESSURE') 
            call InputReadDouble(input,option,sf%constant_capillary_pressure)
            call InputErrorMsg(input,option,'constant capillary pressure', &
                               error_string)
          case('CONSTANT_SATURATION') 
            call InputReadDouble(input,option,sf%constant_saturation)
            call InputErrorMsg(input,option,'constant saturation', &
                                error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'constant saturation function',option)
        end select
    !------------------------------------------
      class is(sat_func_VG_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,sf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'van Genuchten saturation function',option)
        end select
    !------------------------------------------
      class is(sat_func_BC_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,sf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'Brooks-Corey saturation function',option)
        end select
    !------------------------------------------
      class is(sat_func_Linear_type)
        select case(keyword)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'Linear saturation function',option)
        end select
    !------------------------------------------
        class is(sat_func_mK_type)
          select case(keyword)
            case('SIGMAZ')
              call InputReadDouble(input,option,sf%sigmaz)
              call InputErrorMsg(input,option,'sigmaz',error_string)
            case('MUZ')
              call InputReadDouble(input,option,sf%muz)
              call InputErrorMsg(input,option,'muz',error_string)
            case('RMAX')
              call InputReadDouble(input,option,sf%rmax)
              call InputErrorMsg(input,option,'rmax',error_string)
            case('R0')
              call InputReadDouble(input,option,sf%r0)
              call InputErrorMsg(input,option,'r0',error_string)
            case('NPARAM')
              call InputReadInt(input,option,sf%nparam)
              call InputErrorMsg(input,option,'nparam',error_string)
            case default
              call InputKeywordUnrecognized(keyword, &
                   'MODIFIED_KOSUGI saturation function',option)
          end select
    !------------------------------------------
      class is(sat_func_KRP1_type)
        select case(keyword)
          case('KPC') 
            call InputReadInt(input,option,sf%kpc)
            call InputErrorMsg(input,option,'KPC',error_string)
          case('M') 
            call InputReadDouble(input,option,sf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case('PCT_A') 
            call InputReadDouble(input,option,sf%pct_a)
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP') 
            call InputReadDouble(input,option,sf%pct_exp)
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,sf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('IGNORE_PERMEABILITY') 
            sf%ignore_permeability = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP1',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP2_type)
        select case(keyword)
          case('KPC') 
            call InputReadInt(input,option,sf%kpc)
            call InputErrorMsg(input,option,'KPC',error_string)
          case('LAMBDA') 
            call InputReadDouble(input,option,sf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('PCT_A') 
            call InputReadDouble(input,option,sf%pct_a)
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP') 
            call InputReadDouble(input,option,sf%pct_exp)
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('IGNORE_PERMEABILITY') 
            sf%ignore_permeability = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP2',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP3_type)
        select case(keyword)
          case('KPC') 
            call InputReadInt(input,option,sf%kpc)
            call InputErrorMsg(input,option,'KPC',error_string)
          case('LAMBDA') 
            call InputReadDouble(input,option,sf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('PCT_A') 
            call InputReadDouble(input,option,sf%pct_a)
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP') 
            call InputReadDouble(input,option,sf%pct_exp)
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,sf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('IGNORE_PERMEABILITY') 
            sf%ignore_permeability = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP3',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP4_type)
        select case(keyword)
          case('KPC') 
            call InputReadInt(input,option,sf%kpc)
            call InputErrorMsg(input,option,'KPC',error_string)
          case('LAMBDA') 
            call InputReadDouble(input,option,sf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('PCT_A') 
            call InputReadDouble(input,option,sf%pct_a)
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP') 
            call InputReadDouble(input,option,sf%pct_exp)
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,sf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('IGNORE_PERMEABILITY') 
            sf%ignore_permeability = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP4',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP5_type)
        select case(keyword)
          case('KPC') 
            call InputReadInt(input,option,sf%kpc)
            call InputErrorMsg(input,option,'KPC',error_string)
          case('PCT_A') 
            call InputReadDouble(input,option,sf%pct_a)
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP') 
            call InputReadDouble(input,option,sf%pct_exp)
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,sf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('IGNORE_PERMEABILITY') 
            sf%ignore_permeability = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP5',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP8_type)
        select case(keyword)
          case('KPC') 
            call InputReadInt(input,option,sf%kpc)
            call InputErrorMsg(input,option,'KPC',error_string)
          case('M') 
            call InputReadDouble(input,option,sf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case('PCT_A') 
            call InputReadDouble(input,option,sf%pct_a)
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP') 
            call InputReadDouble(input,option,sf%pct_exp)
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,sf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('IGNORE_PERMEABILITY') 
            sf%ignore_permeability = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP8',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP9_type)
        select case(keyword)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP9',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP11_type)
        select case(keyword)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP11',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP12_type)
        select case(keyword)
          case('KPC') 
            call InputReadInt(input,option,sf%kpc)
            call InputErrorMsg(input,option,'KPC',error_string)
          case('PCT_A') 
            call InputReadDouble(input,option,sf%pct_a)
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP') 
            call InputReadDouble(input,option,sf%pct_exp)
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('LAMBDA') 
            call InputReadDouble(input,option,sf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('S_MIN') 
            call InputReadDouble(input,option,sf%s_min)
            call InputErrorMsg(input,option,'s_min',error_string)
          case('S_EFFMIN') 
            call InputReadDouble(input,option,sf%s_effmin)
            call InputErrorMsg(input,option,'s_effmin',error_string)
          case('IGNORE_PERMEABILITY') 
            sf%ignore_permeability = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP12',option)
        end select
    !------------------------------------------
      class default
        option%io_buffer = 'Read routine not implemented for ' &
                           // trim(error_string) // '.'
        call printErrMsg(option)
    !------------------------------------------
    end select
  enddo
  
  if (smooth) then
    call saturation_function%SetupPolynomials(option,error_string)
  endif

  select type(sf => saturation_function)
  !------------------------------------------
    class is(sat_func_constant_type)
      option%io_buffer = 'Constant saturation function is being used.'
      call printWrnMsg(option)
  !------------------------------------------
    class is(sat_func_VG_type)
  !------------------------------------------
    class is(sat_func_BC_type)
      if (.not.smooth) then
        option%io_buffer = 'Brooks-Corey saturation function is being used &
          &without SMOOTH option.'
        call printWrnMsg(option)
      endif
  !------------------------------------------
    class is(sat_func_Linear_type)
  !------------------------------------------
    class is(sat_func_WIPP_type)
      if (sf%ignore_permeability .and. Uninitialized(sf%alpha)) then
        option%io_buffer = 'If a WIPP capillary presure - saturation function &
          &is being used with the IGNORE_PERMEABILITY feature, you must &
          &specify ALPHA (inverse of the threshold capillary pressure). Do &
          &not specify PCT_A or PCT_EXP.'
        call printErrMsg(option)
      endif
  !------------------------------------------
  end select

end subroutine SaturationFunctionRead

! ************************************************************************** !

subroutine PermeabilityFunctionRead(permeability_function,phase_keyword, &
                                    input,option)
  !
  ! Reads in contents of a PERMEABILITY_FUNCTION block
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(rel_perm_func_base_type) :: permeability_function
  character(len=MAXWORDLENGTH) :: phase_keyword
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, new_phase_keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: smooth

  input%ierr = 0
  smooth = PETSC_FALSE
  new_phase_keyword = 'NONE'
  error_string = 'CHARACTERISTIC_CURVES,PERMEABILITY_FUNCTION,'
  select type(rpf => permeability_function)
    class is(rpf_Mualem_VG_liq_type)
      error_string = trim(error_string) // 'MUALEM_VG_LIQ'
    class is(rpf_Mualem_VG_gas_type)
      error_string = trim(error_string) // 'MUALEM_VG_GAS'
    class is(rpf_Burdine_BC_liq_type)
      error_string = trim(error_string) // 'BURDINE_BC_LIQ'
    class is(rpf_Burdine_BC_gas_type)
      error_string = trim(error_string) // 'BURDINE_BC_GAS'
    class is(rpf_TOUGH2_IRP7_gas_type)
      error_string = trim(error_string) // 'TOUGH2_IRP7_GAS'
    class is(rpf_Mualem_BC_liq_type)
      error_string = trim(error_string) // 'MUALEM_BC_LIQ'
    class is(rpf_Mualem_BC_gas_type)
      error_string = trim(error_string) // 'MUALEM_BC_GAS'
    class is(rpf_Burdine_VG_liq_type)
      error_string = trim(error_string) // 'BURDINE_VG_LIQ'
    class is(rpf_Burdine_VG_gas_type)
      error_string = trim(error_string) // 'BURDINE_VG_GAS'
    class is(rpf_Mualem_Linear_liq_type)
      error_string = trim(error_string) // 'MUALEM_Linear_LIQ'
    class is(rpf_Mualem_Linear_gas_type)
      error_string = trim(error_string) // 'MUALEM_Linear_GAS'
    class is(rpf_Burdine_Linear_liq_type)
      error_string = trim(error_string) // 'BURDINE_Linear_LIQ'
    class is(rpf_Burdine_Linear_gas_type)
      error_string = trim(error_string) // 'BURDINE_Linear_GAS'
    class is(rpf_KRP1_liq_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP1_LIQ'
    class is(rpf_KRP1_gas_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP1_GAS'
    class is(rpf_KRP2_liq_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP2_LIQ'
    class is(rpf_KRP2_gas_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP2_GAS'
    class is(rpf_KRP3_liq_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP3_LIQ'
    class is(rpf_KRP3_gas_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP3_GAS'
    class is(rpf_KRP4_liq_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP4_LIQ'
    class is(rpf_KRP4_gas_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP4_GAS'  
    class is(rpf_KRP5_liq_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP5_LIQ'
    class is(rpf_KRP5_gas_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP5_GAS'
    class is(rpf_KRP8_liq_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP8_LIQ'
    class is(rpf_KRP8_gas_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP8_GAS'
    class is(rpf_KRP9_liq_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP9_LIQ'
    class is(rpf_KRP9_gas_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP9_GAS'
    class is(rpf_KRP11_liq_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP11_LIQ'
    class is(rpf_KRP11_gas_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP11_GAS'
    class is(rpf_KRP12_liq_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP12_LIQ'
    class is(rpf_KRP12_gas_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP12_GAS'
    class is(rpf_mK_liq_type)
      error_string = trim(error_string) // 'MODIFIED_KOSUGI_LIQ'
    class is(rpf_mK_gas_type)
      error_string = trim(error_string) // 'MODIFIED_KOSUGI_GAS'
    class is(rpf_TOUGH2_Linear_oil_type)
      error_string = trim(error_string) // 'TOUGH2_Linear_OIL'
    class is(rpf_mod_BC_liq_type)
      error_string = trim(error_string) // 'Mod_BC_LIQ'
    class is(rpf_mod_BC_oil_type)
      error_string = trim(error_string) // 'Mod_BC_OIL'
    class is(rel_perm_func_constant_type)
      error_string = trim(error_string) // 'CONSTANT'
  end select

  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   

    ! base 
    found = PETSC_TRUE
    select case(keyword)
      case('LIQUID_RESIDUAL_SATURATION') 
        call InputReadDouble(input,option,permeability_function%Sr)
        call InputErrorMsg(input,option,'residual_saturation',error_string)
      case('PHASE')
        call InputReadWord(input,option,new_phase_keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'phase',error_string)
        call StringToUpper(phase_keyword) 
      case('SMOOTH')
        smooth = PETSC_TRUE
      case default
        found = PETSC_FALSE
    end select
    if (found) cycle

    select type(rpf => permeability_function)
    !------------------------------------------
      class is(rpf_Mualem_VG_liq_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem van Genuchten liquid relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Mualem_VG_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem van Genuchten gas relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Burdine_BC_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine Brooks-Corey liquid relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Burdine_BC_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine Brooks-Corey gas relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_TOUGH2_IRP7_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'TOUGH2 IRP7 gas relative permeability function',option)
        end select
    !------------------------------------------
      class is(rpf_Mualem_BC_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem Brooks-Corey liquid relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Mualem_BC_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem Brooks-Corey gas relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Burdine_VG_liq_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine van Genuchten liquid relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Burdine_VG_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine van Genuchten gas relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Mualem_Linear_liq_type)
        select case(keyword)
          case('MAX_CAPILLARY_PRESSURE') 
            call InputReadDouble(input,option,rpf%pcmax)
            call InputErrorMsg(input,option,'max_capillary_pressure',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,rpf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem Linear liquid relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Mualem_Linear_gas_type)
        select case(keyword)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case('MAX_CAPILLARY_PRESSURE') 
            call InputReadDouble(input,option,rpf%pcmax)
            call InputErrorMsg(input,option,'max_capillary_pressure',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,rpf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem Linear gas relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Burdine_Linear_liq_type)
        select case(keyword)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine Linear liquid relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Burdine_Linear_gas_type)
        select case(keyword)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine Linear gas relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP1_liq_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP1_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP1_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP1_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP2_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP2_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP2_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP2_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP3_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP3_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP3_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP3_GAS relative permeability function', &
              option)
        end select
        
    !------------------------------------------
      class is(rpf_KRP4_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP4_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP4_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP4_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP5_liq_type)
        select case(keyword)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP5_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP5_gas_type)
        select case(keyword)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP5_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP8_liq_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP8_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP8_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP8_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP9_liq_type)
        select case(keyword)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP9_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP9_gas_type)
        select case(keyword)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP9_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP11_liq_type)
        select case(keyword)
          case('TOLC') 
            call InputReadDouble(input,option,rpf%tolc)
            call InputErrorMsg(input,option,'TOLC',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                 error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP11_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP11_gas_type)
        select case(keyword)
          case('TOLC') 
            call InputReadDouble(input,option,rpf%tolc)
            call InputErrorMsg(input,option,'TOLC',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                 error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP11_GAS relative permeability function', &
              option)
        end select  
    !------------------------------------------
      class is(rpf_KRP12_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                 error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP12_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP12_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                 error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO_KRP12_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_mK_liq_type)
        select case(keyword)
          case('SIGMAZ')
            call InputReadDouble(input,option,rpf%sigmaz)
            call InputErrorMsg(input,option,'sigmaz',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                 'MODIFIED_KOSUGI liquid relative permeability '//&
                 &'function',option)
        end select
    !------------------------------------------
      class is(rpf_mK_gas_type)
        select case(keyword)
          case('SIGMAZ')
            call InputReadDouble(input,option,rpf%sigmaz)
            call InputErrorMsg(input,option,'sigmaz',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                 'MODIFIED_KOSUGI gas relative permeability '//&
                 &'function',option)
        end select
    !------------------------------------------
      class is(rpf_TOUGH2_Linear_oil_type)
        select case(keyword)
          case('OIL_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Sro)
            call InputErrorMsg(input,option,'Sro',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'TOUGH2 LINEAR oil relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_mod_BC_liq_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m - power',error_string)
          case('OIL_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Sro)
            call InputErrorMsg(input,option,'Sro',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case('LIQUID_MAX_REL_PERM') 
            call InputReadDouble(input,option,rpf%kr_max)
            call InputErrorMsg(input,option,'kr_max',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mod BC liq relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_mod_BC_oil_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m - power',error_string)
          case('OIL_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Sro)
            call InputErrorMsg(input,option,'Sro',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case('OIL_MAX_REL_PERM') 
            call InputReadDouble(input,option,rpf%kr_max)
            call InputErrorMsg(input,option,'kr_max',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mod BC oil relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rel_perm_func_constant_type)
        select case(keyword)
          case('RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Sr)
            call InputErrorMsg(input,option,'Sr',error_string)
          case('RELATIVE_PERMEABILITY') 
            call InputReadDouble(input,option,rpf%kr)
            call InputErrorMsg(input,option,'kr',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Constant relative permeability function', &
              option)
        end select
    !------------------------------------------
      class default
        option%io_buffer = 'Read routine not implemented for relative ' // &
                           'permeability function class.'
        call printErrMsg(option)
    !------------------------------------------
    end select
  enddo

  
  ! for functions that are not phase-specific, check if PHASE was given:
  if (StringCompare('NONE',phase_keyword)) then
    phase_keyword = new_phase_keyword
    ! a liq or gas phase should now be specified for the non-phase-specific
    ! functions, so check if it was:
    if (StringCompare('NONE',new_phase_keyword)) then
      ! entering means the new phase keyword was also NONE (the default), so
      ! throw an error and abort:
      option%io_buffer = 'PHASE is not specified for ' // trim(error_string) 
      call printErrMsg(option)
    endif
  endif
  
  ! liquid phase relative permeability function check:
  if (StringCompare('LIQUID',phase_keyword)) then
    if (StringCompare('GAS',new_phase_keyword)) then
      ! user is requesting a liquid relative perm func for a gas phase:
      option%io_buffer = 'A liquid-phase relative permeability function &
                         &is being requested for the gas phase under ' &
                         // trim(error_string) // '.'
      call printErrMsg(option)
    endif
  endif
  
  ! gas phase relative permeability function check:
  if (StringCompare('GAS',phase_keyword)) then
    if (StringCompare('LIQUID',new_phase_keyword)) then
      ! user is requesting a gas relative perm func for a liquid phase:
      option%io_buffer = 'A gas-phase relative permeability function &
                         &is being requested for the liquid phase under ' &
                         // trim(error_string) // '.'
      call printErrMsg(option)
    endif
  endif

  if (smooth) then
    call permeability_function%SetupPolynomials(option,error_string)
  endif
  
end subroutine PermeabilityFunctionRead

! ************************************************************************** !

subroutine CharacteristicCurvesAddToList(new_characteristic_curves,list)
  !
  ! Adds a characteristic curves object to linked list
  !
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  !

  implicit none

  class(characteristic_curves_type), pointer :: new_characteristic_curves
  class(characteristic_curves_type), pointer :: list

  class(characteristic_curves_type), pointer :: cur_characteristic_curves

  if (associated(list)) then
    cur_characteristic_curves => list
    ! loop to end of list
    do
      if (.not.associated(cur_characteristic_curves%next)) exit
      cur_characteristic_curves => cur_characteristic_curves%next
    enddo
    cur_characteristic_curves%next => new_characteristic_curves
  else
    list => new_characteristic_curves
  endif

end subroutine CharacteristicCurvesAddToList

! ************************************************************************** !

subroutine CharCurvesConvertListToArray(list,array,option)
  !
  ! Creates an array of pointers to the characteristic curves objects in the
  ! list
  !
  ! Author: Glenn Hammond
  ! Date: 12/11/07
  !

  use String_module
  use Option_module

  implicit none

  class(characteristic_curves_type), pointer :: list
  type(characteristic_curves_ptr_type), pointer :: array(:)
  type(option_type) :: option

  class(characteristic_curves_type), pointer :: cur_characteristic_curves
  PetscInt :: count

  count = 0
  cur_characteristic_curves => list
  do
    if (.not.associated(cur_characteristic_curves)) exit
    count = count + 1
    cur_characteristic_curves => cur_characteristic_curves%next
  enddo

  if (associated(array)) deallocate(array)
  allocate(array(count))

  count = 0
  cur_characteristic_curves => list
  do
    if (.not.associated(cur_characteristic_curves)) exit
    count = count + 1
    array(count)%ptr => cur_characteristic_curves
    if (cur_characteristic_curves%test .and. &
        option%myrank == option%io_rank) then
      call CharacteristicCurvesTest(cur_characteristic_curves,option)
    endif
    cur_characteristic_curves => cur_characteristic_curves%next
  enddo

end subroutine CharCurvesConvertListToArray

! ************************************************************************** !

function CharCurvesGetGetResidualSats(characteristic_curves,option)
  ! 
  ! Returns the residual saturations associated with a characteristic curves
  ! object
  ! 
  ! Author: Glenn Hammond, Paolo Orsini
  ! Date: 09/29/14, 03/27/17
  ! 

  use Option_module
  
  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option

  PetscReal :: CharCurvesGetGetResidualSats(option%nphase)

  select case(option%iflowmode)
    case(TOWG_MODE)
      CharCurvesGetGetResidualSats(option%liquid_phase) = &
        characteristic_curves%wat_rel_perm_func_owg%Swcr
      CharCurvesGetGetResidualSats(option%oil_phase) = &
        characteristic_curves%oil_rel_perm_func_owg%Socr
      if (option%iflow_sub_mode == TOWG_TODD_LONGSTAFF) then
        CharCurvesGetGetResidualSats(option%gas_phase) = 0.0d0
      else
        CharCurvesGetGetResidualSats(option%gas_phase) = &
          characteristic_curves%gas_rel_perm_func_owg%Sgcr
      end if
    case default
  CharCurvesGetGetResidualSats(1) = &
    characteristic_curves%liq_rel_perm_function%Sr
  !if (option%nphase > 1) then
  if ( (option%nphase > 1) .and. &
       associated(characteristic_curves%gas_rel_perm_function) ) then
    select type(rpf=>characteristic_curves%gas_rel_perm_function)
      class is(rpf_Mualem_VG_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_Mualem_VG_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_Burdine_BC_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_Burdine_BC_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_Mualem_BC_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_Mualem_BC_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_Burdine_VG_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_Burdine_VG_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_TOUGH2_IRP7_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_Mualem_Linear_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_Mualem_Linear_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_Burdine_Linear_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_Burdine_Linear_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_KRP1_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_KRP1_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_KRP2_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_KRP2_gas_type)
        ! KRP2 does not use a Srg, so return 0.d0
        CharCurvesGetGetResidualSats(option%gas_phase) = 0.d0
      class is(rpf_KRP3_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_KRP3_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_KRP4_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_KRP4_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_KRP5_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_KRP5_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_KRP8_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_KRP8_gas_type)
        ! KRP8 does not use a Srg, so return 0.d0
        CharCurvesGetGetResidualSats(option%gas_phase) = 0.d0
      class is(rpf_KRP9_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_KRP9_gas_type)
        ! KRP9 does not use a Srg, so return 0.d0
        CharCurvesGetGetResidualSats(option%gas_phase) = 0.d0
      class is(rpf_KRP11_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_KRP11_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_KRP12_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_KRP12_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_mK_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_mK_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_mod_BC_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rel_perm_func_constant_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rel_perm_func_default_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class default
        option%io_buffer = 'Relative permeability class not supported in &
              &CharCurvesGetGetResidualSats. &
              &Contact pflotran-dev@googlegroups.com'
        call printErrMsg(option)
    end select

  end if

  if ( (option%nphase > 1) .and. &
       associated(characteristic_curves%oil_rel_perm_function) ) then
    select type(rpf=>characteristic_curves%oil_rel_perm_function)
      class is(rpf_TOUGH2_Linear_oil_type)
        CharCurvesGetGetResidualSats(option%oil_phase) = rpf%Sro
      class is(rpf_mod_BC_oil_type)
        CharCurvesGetGetResidualSats(option%oil_phase) = rpf%Sro 
      class default
        option%io_buffer = 'Oil Relative permeability class ' // &
          'not supported in CharCurvesGetGetResidualSats.'
        call printErrMsg(option)
    end select
  endif

  end select ! end flow mode select

end function CharCurvesGetGetResidualSats

! ************************************************************************** !

function CharacteristicCurvesGetID(characteristic_curves_array, &
                                   characteristic_curves_name, &
                                   material_property_name, option)
  ! 
  ! Returns the ID of the characteristic curves object named
  ! "characteristic_curves_name"
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  use Option_module
  use String_module
  
  type(characteristic_curves_ptr_type), pointer :: &
    characteristic_curves_array(:)
  character(len=MAXWORDLENGTH) :: characteristic_curves_name
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: CharacteristicCurvesGetID

  CharacteristicCurvesGetID = 0
  do CharacteristicCurvesGetID = 1, size(characteristic_curves_array)
    if (StringCompare(characteristic_curves_name, &
                      characteristic_curves_array( &
                        CharacteristicCurvesGetID)%ptr%name)) then
      return
    endif
  enddo
  option%io_buffer = 'Characteristic curves "' // &
           trim(characteristic_curves_name) // &
           '" in material property "' // &
           trim(material_property_name) // &
           '" not found among available characteristic curves.'
  call printErrMsg(option)    

end function CharacteristicCurvesGetID

! ************************************************************************** !

subroutine CharacteristicCurvesTest(characteristic_curves,option)
  ! 
  ! Outputs values of characteristic curves over a range of values
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  use Option_module

  implicit none
  
  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: phase

  if (associated(characteristic_curves%saturation_function)) then
  call characteristic_curves%saturation_function%Test( &
                                                 characteristic_curves%name, &
                                                 option)
  end if

  if (associated(characteristic_curves%oil_wat_sat_func)) then
    call characteristic_curves%oil_wat_sat_func%Test( &
                                                 characteristic_curves%name, &
                                                 option)
  end if

  if (associated(characteristic_curves%oil_gas_sat_func)) then
    call characteristic_curves%oil_gas_sat_func%Test( &
                                                 characteristic_curves%name, &
                                                 option)
  end if

  if (associated(characteristic_curves%liq_rel_perm_function)) then
  phase = 'liquid'
  call characteristic_curves%liq_rel_perm_function%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
  end if
              
  if ( associated(characteristic_curves%gas_rel_perm_function) ) then
    phase = 'gas'
    call characteristic_curves%gas_rel_perm_function%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
  endif

  if ( associated(characteristic_curves%oil_rel_perm_function) ) then
    phase = 'oil'
    call characteristic_curves%oil_rel_perm_function%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
  end if
  
  if ( associated(characteristic_curves%wat_rel_perm_func_owg) ) then
    phase = 'water'
    call characteristic_curves%wat_rel_perm_func_owg%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
  end if

  if ( associated(characteristic_curves%oil_rel_perm_func_owg) ) then
    phase = 'oil'
    call characteristic_curves%oil_rel_perm_func_owg%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
  end if

  if ( associated(characteristic_curves%gas_rel_perm_func_owg) ) then
    phase = 'gas'
    call characteristic_curves%gas_rel_perm_func_owg%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
  end if

end subroutine CharacteristicCurvesTest

! ************************************************************************** !

subroutine CharacteristicCurvesVerify(characteristic_curves,option)
  ! 
  ! Checks if required parameters have been set for each curve type.
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  use Option_module

  implicit none
  
  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  string = 'CHARACTERISTIC_CURVES(' // trim(characteristic_curves%name) // &
           '),'

  if (associated(characteristic_curves%saturation_function)) then
    call characteristic_curves%saturation_function%Verify(string,option)
  else
    option%io_buffer = 'A saturation function has &
                       &not been set under CHARACTERISTIC_CURVES "' // &
                       trim(characteristic_curves%name) // '". A &
                       &PERMEABILITY_FUNCTION block must be specified &
                       &for the liquid phase.'
  endif
  
  if (associated(characteristic_curves%liq_rel_perm_function) ) then
    call characteristic_curves%liq_rel_perm_function%Verify(string,option)
  else
    option%io_buffer = 'A liquid phase relative permeability function has &
                       &not been set under CHARACTERISTIC_CURVES "' // &
                       trim(characteristic_curves%name) // '". A &
                       &PERMEABILITY_FUNCTION block must be specified &
                       &for the liquid phase.'
    call printErrMsg(option)
  end if

  if (associated(characteristic_curves%gas_rel_perm_function) ) then
    call characteristic_curves%gas_rel_perm_function%Verify(string,option)
  else
    if (option%iflowmode == G_MODE .or. option%iflowmode == TOWG_MODE .or. &
        option%iflowmode == WF_MODE) then
      option%io_buffer = 'A gas phase relative permeability function has &
                         &not been set under CHARACTERISTIC_CURVES "' // &
                         trim(characteristic_curves%name) // '". Another &
                         &PERMEABILITY_FUNCTION block must be specified &
                         &for the gas phase.'
      call printErrMsg(option)
    end if
  end if

  if ( associated(characteristic_curves%oil_rel_perm_function) ) then  
    call characteristic_curves%oil_rel_perm_function%Verify(string,option)
  else 
    if (option%iflowmode == TOIL_IMS_MODE .or. &
        option%iflowmode == TOWG_MODE  ) then
      option%io_buffer = 'An oil phase relative permeability function has &
                         &not been set under CHARACTERISTIC_CURVES "' // &
                         trim(characteristic_curves%name) // '". Another &
                         &PERMEABILITY_FUNCTION block must be specified &
                         &for the oil phase.'      
    end if
  end if
  
end subroutine CharacteristicCurvesVerify

! **************************************************************************** !

subroutine CharacteristicCurvesOWGVerify(characteristic_curves,option)
  !
  ! Checks if required parameters have been set for each curve type.
  !
  ! Author: Paolo Orsini
  ! Date: 11/16/17
  !
  use Option_module

  implicit none

  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  string = 'CHARACTERISTIC_CURVES(' // trim(characteristic_curves%name) // &
           '),'

  if (.not.(associated(characteristic_curves%oil_wat_sat_func))) then
    option%io_buffer = "A water/oil saturation function has  &
                        &not been set under CHARACTERISTIC_CURVES " // &
                        trim(characteristic_curves%name) // '". A &
                        &SATURATION_FUNCTION_OWG block must be &
                        &specified for the water/oil phase interface.'
    call printErrMsg(option)
  else
    call characteristic_curves%oil_wat_sat_func%verify(string,option)
  end if

  select case(option%iflow_sub_mode)
    case(TOWG_TODD_LONGSTAFF)
      !do nothing - for TL only one saturation function is required
    case default
      if (.not.(associated(characteristic_curves%oil_gas_sat_func))) then
        option%io_buffer = "An oil/gas saturation function has  &
                            &not been set under CHARACTERISTIC_CURVES " // &
                            trim(characteristic_curves%name) // '". A &
                            &SATURATION_FUNCTION_OWG block must be &
                            &specified for the oil/gas phase interface.'
        call printErrMsg(option)
      else
        call characteristic_curves%oil_gas_sat_func%verify(string,option)
      end if
  end select

  ! Verify relative permeabilities
  if (.not.(associated(characteristic_curves%wat_rel_perm_func_owg))) then
    option%io_buffer = "A water relative permeability function has &
                        &not been set under CHARACTERISTIC_CURVES " // &
                        trim(characteristic_curves%name) // '". A &
                        &PERMEABILITY_FUNCTION_OWG block must be &
                        &specified for the water phase.'
    call printErrMsg(option)
  else
    call characteristic_curves%wat_rel_perm_func_owg%verify(string,option)
  end if

  if (.not.(associated(characteristic_curves%oil_rel_perm_func_owg))) then
    option%io_buffer = "An oil relative permeability function has &
                        &not been set under CHARACTERISTIC_CURVES " // &
                        trim(characteristic_curves%name) // '". A &
                        &PERMEABILITY_FUNCTION_OWG block must be &
                        &specified for the oil phase.'
    call printErrMsg(option)
  else
    call characteristic_curves%oil_rel_perm_func_owg%verify(string,option)
  end if

  !check if end points in wat_rel_perm and oil_rel_perm are consistent
  !PO todo: add consistency checks for all functions in a seperate routine.
  !         Check Swcr vs Slcr: In some functions/tables only Slcr might
  !         be defined, thuse check Slcr == Socr + Swcr
  ! if (characteristic_curves%wat_rel_perm_func_owg%Swcr /= &
  !     characteristic_curves%oil_rel_perm_func_owg%Swcr
  !    ) then
  !    option%io_buffer = "Water critical saturation in the water and oil &
  !                       &relative permeability functions are different &
  !                        &for the CHARACTERISTIC_CURVES " // &
  !                        trim(characteristic_curves%name) // '". Please &
  !                        &ensure these to sse same values'
  ! end if

  select case(option%iflowmode)
    case(TOIL_IMS_MODE)
      ! do nothing - gas phase no defined in TOIL_IMS
    case(TOWG_MODE)
      if ( option%iflow_sub_mode /= TOWG_TODD_LONGSTAFF ) then
        if (.not.(associated(characteristic_curves%gas_rel_perm_func_owg)) &
           ) then
          option%io_buffer = "A gas relative permeability function has &
                              &not been set under CHARACTERISTIC_CURVES " // &
                              trim(characteristic_curves%name) // '". A &
                              &PERMEABILITY_FUNCTION_OWG block must be &
                              &specified for the gas phase.'
          call printErrMsg(option)
        else
          call characteristic_curves%gas_rel_perm_func_owg% &
                                                        verify(string,option)

        end if
      end if
  end select

end subroutine CharacteristicCurvesOWGVerify

! **************************************************************************** !

subroutine CharCurvesInputRecord(char_curve_list)
  ! 
  ! Prints ingested characteristic curves information to the input record file
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/11/2016
  ! 

  implicit none

  class(characteristic_curves_type), pointer :: char_curve_list
  
  class(characteristic_curves_type), pointer :: cur_ccurve
  character(len=MAXWORDLENGTH) :: word1
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'CHARACTERISTIC CURVES'
  
  cur_ccurve => char_curve_list
  do
    if (.not.associated(cur_ccurve)) exit
    
    write(id,'(a29)',advance='no') 'characteristic curve name: '
    write(id,'(a)') adjustl(trim(cur_ccurve%name))
    
    if (associated(cur_ccurve%saturation_function)) then
      write(id,'(a29)',advance='no') 'saturation function: '
      select type (sf => cur_ccurve%saturation_function)
      !---------------------------------
        class is (sat_func_VG_type)
          write(id,'(a)') 'van Genuchten'
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) sf%m
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'alpha: '
          write(word1,*) sf%alpha
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_BC_type)
          write(id,'(a)') 'Brooks Corey'
          write(id,'(a29)',advance='no') 'alpha: '
          write(word1,*) sf%alpha
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) sf%lambda
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_Linear_type)
          write(id,'(a)') 'linear'
          write(id,'(a29)',advance='no') 'alpha: '
          write(word1,*) sf%alpha
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_mK_type)
          write(id,'(a)') 'Modified Kosugi'
          write(id,'(a29)',advance='no') 'sigmaz: '
          write(word1,*) sf%sigmaz
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'muz: '
          write(word1,*) sf%muz
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'liquid residual sat.: '
          write(word1,*) sf%Sr
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'rmax: '
          write(word1,*) sf%rmax
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'r0: '
          write(word1,*) sf%r0
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_KRP1_type)
          write(id,'(a)') 'Bragflo KRP1 modified van Genuchten'
          write(id,'(a29)',advance='no') 'kpc: '
          write(word1,*) sf%kpc
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_a: '
          write(word1,*) sf%pct_a
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_exp: '
          write(word1,*) sf%pct_exp
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) sf%m
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) sf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_KRP2_type)
          write(id,'(a)') 'Bragflo KRP2 original Brooks Corey'
          write(id,'(a29)',advance='no') 'kpc: '
          write(word1,*) sf%kpc
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_a: '
          write(word1,*) sf%pct_a
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_exp: '
          write(word1,*) sf%pct_exp
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) sf%lambda
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_KRP3_type)
          write(id,'(a)') 'Bragflo KRP3 1st modified Brooks Corey'
          write(id,'(a29)',advance='no') 'kpc: '
          write(word1,*) sf%kpc
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_a: '
          write(word1,*) sf%pct_a
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_exp: '
          write(word1,*) sf%pct_exp
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) sf%lambda
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) sf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_KRP4_type)
          write(id,'(a)') 'Bragflo KRP4 2nd modified Brooks Corey'
          write(id,'(a29)',advance='no') 'kpc: '
          write(word1,*) sf%kpc
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_a: '
          write(word1,*) sf%pct_a
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_exp: '
          write(word1,*) sf%pct_exp
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) sf%lambda
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) sf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_KRP5_type)
          write(id,'(a)') 'Bragflo KRP5 modified Linear'
          write(id,'(a29)',advance='no') 'kpc: '
          write(word1,*) sf%kpc
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_a: '
          write(word1,*) sf%pct_a
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_exp: '
          write(word1,*) sf%pct_exp
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) sf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_KRP8_type)
          write(id,'(a)') 'Bragflo KRP8 original van Genuchten-Parker'
          write(id,'(a29)',advance='no') 'kpc: '
          write(word1,*) sf%kpc
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_a: '
          write(word1,*) sf%pct_a
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_exp: '
          write(word1,*) sf%pct_exp
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) sf%m
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) sf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_KRP9_type)
          write(id,'(a)') 'Bragflo KRP9 Vauchlin infiltration test'
      !---------------------------------
        class is (sat_func_KRP11_type)
          write(id,'(a)') 'Bragflo KRP11 open cavity modification'
      !---------------------------------
        class is (sat_func_KRP12_type)
          write(id,'(a)') 'Bragflo KRP12 waste area modification'
          write(id,'(a29)',advance='no') 'kpc: '
          write(word1,*) sf%kpc
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_a: '
          write(word1,*) sf%pct_a
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'pct_exp: '
          write(word1,*) sf%pct_exp
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) sf%lambda
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 's_min: '
          write(word1,*) sf%s_min
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 's_effmin: '
          write(word1,*) sf%s_effmin
          write(id,'(a)') adjustl(trim(word1))
      !---------------------------------
        class is (sat_func_default_type)
          write(id,'(a)') 'default'
      !---------------------------------
      end select
      write(id,'(a29)',advance='no') 'liquid residual sat.: '
      write(word1,*) cur_ccurve%saturation_function%Sr
      write(id,'(a)') adjustl(trim(word1))
      write(id,'(a29)',advance='no') 'max capillary pressure: '
      write(word1,*) cur_ccurve%saturation_function%pcmax
      write(id,'(a)') adjustl(trim(word1))
    endif
    
    if (associated(cur_ccurve%liq_rel_perm_function)) then
      write(id,'(a29)',advance='no') 'liq. relative perm. func.: '
      select type (rpf => cur_ccurve%liq_rel_perm_function)
      !------------------------------------
        class is (rel_perm_func_default_type)
          write(id,'(a)') 'default'
      !------------------------------------
        class is (rpf_Mualem_VG_liq_type)
          write(id,'(a)') 'mualem_vg_liq/tough2_irp7_liq'
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) rpf%m
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_Mualem_BC_liq_type)
          write(id,'(a)') 'mualem_bc_liq'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_Mualem_Linear_liq_type)
          write(id,'(a)') 'mualem_linear_liq'
          write(id,'(a29)',advance='no') 'alpha: '
          write(word1,*) rpf%alpha
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'max capillary pressure: '
          write(word1,*) rpf%pcmax
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_Burdine_VG_liq_type)
          write(id,'(a)') 'burdine_vg_liq'
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) rpf%m
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_Burdine_BC_liq_type)
          write(id,'(a)') 'burdine_bc_liq'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_Burdine_Linear_liq_type)
          write(id,'(a)') 'burdine_linear_liq'
      !------------------------------------
        class is (rpf_KRP1_liq_type)
          write(id,'(a)') 'Bragflo KRP1 liquid'
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) rpf%m
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP2_liq_type)
          write(id,'(a)') 'Bragflo KRP2 liquid'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP3_liq_type)
          write(id,'(a)') 'Bragflo KRP3 liquid'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP4_liq_type)
          write(id,'(a)') 'Bragflo KRP4 liquid'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP5_liq_type)
          write(id,'(a)') 'Bragflo KRP5 liquid'
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP8_liq_type)
          write(id,'(a)') 'Bragflo KRP1 liquid'
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) rpf%m
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP9_liq_type)
          write(id,'(a)') 'Bragflo KRP9 liquid'
      !------------------------------------
        class is (rpf_KRP11_liq_type)
          write(id,'(a)') 'Bragflo KRP11 liquid'
          write(id,'(a29)',advance='no') 'tolc: '
          write(word1,*) rpf%tolc
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP12_liq_type)
          write(id,'(a)') 'Bragflo KRP12 liquid'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_mK_liq_type)
          write(id,'(a)') 'modified_kosugi_liq'
          write(id,'(a29)',advance='no') 'sigmaz: '
          write(word1,*) rpf%sigmaz
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'liquid residual sat.: '
          write(word1,*) rpf%Sr
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class default
          write(id,'(a)') 'none'
      !------------------------------------
      end select
    endif
    
    if (associated(cur_ccurve%gas_rel_perm_function)) then
      write(id,'(a29)',advance='no') 'gas relative perm. func.: '
      select type (rpf => cur_ccurve%gas_rel_perm_function)
      !------------------------------------
        class is (rel_perm_func_default_type)
          write(id,'(a)') 'default'
      !------------------------------------
        class is (rpf_Mualem_VG_gas_type)
          write(id,'(a)') 'mualem_vg_gas'
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) rpf%m
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_Mualem_BC_gas_type)
          write(id,'(a)') 'mualem_bc_gas'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_Mualem_Linear_gas_type)
          write(id,'(a)') 'mualem_linear_gas'
          write(id,'(a29)',advance='no') 'alpha: '
          write(word1,*) rpf%alpha
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'max capillary pressure: '
          write(word1,*) rpf%pcmax
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_TOUGH2_IRP7_gas_type)
          write(id,'(a)') 'tough2_irp7_gas'
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) rpf%m
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_Burdine_VG_gas_type)
          write(id,'(a)') 'burdine_vg_gas'
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) rpf%m
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_Burdine_BC_gas_type)
          write(id,'(a)') 'burdine_bc_gas'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_Burdine_linear_gas_type)
          write(id,'(a)') 'burdine_linear_gas'
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP1_gas_type)
          write(id,'(a)') 'Bragflo KRP1 gas'
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) rpf%m
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP2_gas_type)
          write(id,'(a)') 'Bragflo KRP2 gas'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP3_gas_type)
          write(id,'(a)') 'Bragflo KRP3 gas'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP4_gas_type)
          write(id,'(a)') 'Bragflo KRP4 gas'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP5_gas_type)
          write(id,'(a)') 'Bragflo KRP5 gas'
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP8_gas_type)
          write(id,'(a)') 'Bragflo KRP1 gas'
          write(id,'(a29)',advance='no') 'm: '
          write(word1,*) rpf%m
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP9_gas_type)
          write(id,'(a)') 'Bragflo KRP9 gas'
      !------------------------------------
        class is (rpf_KRP11_gas_type)
          write(id,'(a)') 'Bragflo KRP11 gas'
          write(id,'(a29)',advance='no') 'tolc: '
          write(word1,*) rpf%tolc
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_KRP12_gas_type)
          write(id,'(a)') 'Bragflo KRP12 gas'
          write(id,'(a29)',advance='no') 'lambda: '
          write(word1,*) rpf%lambda
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class is (rpf_mK_gas_type)
          write(id,'(a)') 'modified_kosugi_gas'
          write(id,'(a29)',advance='no') 'sigmaz: '
          write(word1,*) rpf%sigmaz
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'gas residual sat.: '
          write(word1,*) rpf%Srg
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
        class default
          write(id,'(a)') 'none'
      !------------------------------------
      end select
    endif
    
    if (associated(cur_ccurve%oil_rel_perm_function)) then
      write(id,'(a29)',advance='no') 'oil relative perm. func.: '
      select type (rpf => cur_ccurve%oil_rel_perm_function)
      !------------------------------------
        class is (rel_perm_func_default_type)
          write(id,'(a)') 'default'
      !------------------------------------
        class is (rpf_TOUGH2_Linear_Oil_type)
          write(id,'(a)') 'tough2_linear_oil'
          write(id,'(a29)',advance='no') 'oil residual sat.: '
          write(word1,*) rpf%Sro
          write(id,'(a)') adjustl(trim(word1))
      !------------------------------------
      end select
    endif

    write(id,'(a29)') '---------------------------: '
    cur_ccurve => cur_ccurve%next
  enddo
  
end subroutine CharCurvesInputRecord


! ************************************************************************** !

recursive subroutine CharacteristicCurvesDestroy(cc)
  ! 
  ! Destroys a characteristic curve
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(characteristic_curves_type), pointer :: cc
  
  if (.not.associated(cc)) return
  
  call CharacteristicCurvesDestroy(cc%next)
  
  call SaturationFunctionDestroy(cc%saturation_function)

  call SaturationFunctionOWGDestroy(cc%oil_wat_sat_func)
  call SaturationFunctionOWGDestroy(cc%oil_gas_sat_func)

  ! the liquid and gas relative permeability pointers may pointer to the
  ! same address. if so, destroy one and nullify the other.
  if (associated(cc%liq_rel_perm_function,cc%gas_rel_perm_function)) then
    call PermeabilityFunctionDestroy(cc%liq_rel_perm_function)
    nullify(cc%gas_rel_perm_function)
  !PO how about avoiding xxx_rel_perm_function => aaa_rel_perm_function? 
  !   it should semplify code. It seems we do this only to pass verify 
  else if (associated(cc%oil_rel_perm_function,cc%gas_rel_perm_function)) then 
    call PermeabilityFunctionDestroy(cc%oil_rel_perm_function)
    nullify(cc%gas_rel_perm_function)
  else
    call PermeabilityFunctionDestroy(cc%liq_rel_perm_function)
    call PermeabilityFunctionDestroy(cc%gas_rel_perm_function)
    !call PermeabilityFunctionDestroy(cc%oil_rel_perm_function)
  endif

  call PermeabilityFunctionOWGDestroy(cc%wat_rel_perm_func_owg)
  call PermeabilityFunctionOWGDestroy(cc%oil_rel_perm_func_owg)
  call PermeabilityFunctionOWGDestroy(cc%gas_rel_perm_func_owg)

  deallocate(cc)
  nullify(cc)
  
end subroutine CharacteristicCurvesDestroy

! ************************************************************************** !

end module Characteristic_Curves_module
