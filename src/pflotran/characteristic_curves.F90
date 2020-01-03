module Characteristic_Curves_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_Common_module

  implicit none

  private

  type, public :: characteristic_curves_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
    PetscBool :: test
    class(sat_func_base_type), pointer :: saturation_function
    class(rel_perm_func_base_type), pointer :: liq_rel_perm_function
    class(rel_perm_func_base_type), pointer :: gas_rel_perm_function
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
  public :: CharacteristicCurvesCreateCopy
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
  nullify(characteristic_curves%liq_rel_perm_function)
  nullify(characteristic_curves%gas_rel_perm_function)
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
  character(len=MAXSTRINGLENGTH) :: error_string
  class(rel_perm_func_base_type), pointer :: rel_perm_function_ptr

  nullify(rel_perm_function_ptr)

  input%ierr = 0
  error_string = 'CHARACTERISTIC_CURVES' 
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    !-----------------------------------------------------------------------
      case('SATURATION_FUNCTION')
        call InputReadCardDbaseCompatible(input,option,word)
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
          case('IGHCC2_COMP')
            this%saturation_function => SF_IGHCC2_Comp_Create()
          case default
            call InputKeywordUnrecognized(input,word,'SATURATION_FUNCTION', &
                                          option)
        end select
        call SaturationFunctionRead(this%saturation_function,input,option)
    !-----------------------------------------------------------------------
      case('PERMEABILITY_FUNCTION')
        nullify(rel_perm_function_ptr)
        phase_keyword = 'NONE'
        call InputReadCardDbaseCompatible(input,option,word)
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
          case('MODIFIED_KOSUGI_LIQ')
            rel_perm_function_ptr => RPF_mK_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MODIFIED_KOSUGI_GAS')
            rel_perm_function_ptr => RPF_mK_Gas_Create()
            phase_keyword = 'GAS'
          case('IGHCC2_COMP_LIQ')
            rel_perm_function_ptr => RPF_IGHCC2_Comp_Liq_Create()
            phase_keyword = 'LIQUID'
          case('IGHCC2_COMP_GAS')
            rel_perm_function_ptr => RPF_IGHCC2_Comp_Gas_Create()
            phase_keyword = 'GAS'
          case('CONSTANT')
            rel_perm_function_ptr => RPF_Constant_Create()
            ! phase_keyword = 'NONE'
          case default
            call InputKeywordUnrecognized(input,word,'PERMEABILITY_FUNCTION', &
                                          option)
        end select
        call PermeabilityFunctionRead(rel_perm_function_ptr,phase_keyword, &
                                      input,option)
        ! if PHASE is specified, have to align correct pointer
        select case(phase_keyword)
          case('GAS')
            this%gas_rel_perm_function => rel_perm_function_ptr
          case('LIQUID')
            this%liq_rel_perm_function => rel_perm_function_ptr
          case('NONE')
            option%io_buffer = 'PHASE has not been set for &
                               &CHARACTERISTIC_CURVES,PERMEABILITY_FUNCTION. &
                               &This is most likely a development issue, and &
                               &not an input deck mistake. '
            call PrintErrMsgToDev(option,'')
          case default
            call InputKeywordUnrecognized(input,word, &
                                          'PERMEABILITY_FUNCTION,PHASE',option)
        end select
      case('PERMEABILITY_FUNCTION_WAT','KRW')
        call InputReadCardDbaseCompatible(input,option,word)
        call InputErrorMsg(input,option,'PERMEABILITY_FUNCTION_WAT - KRW ', &
                           error_string)
        call StringToUpper(word)
        select case(word)
          case('MOD_BROOKS_COREY')
            this%wat_rel_perm_func_owg => RPF_wat_owg_MBC_Create()
          case('MUALEM_VG')
            this%wat_rel_perm_func_owg => RPF_wat_owg_Mualem_VG_Create()
          case('BURDINE_VG')
            this%wat_rel_perm_func_owg => RPF_wat_owg_Burdine_VG_Create()
          case('BURDINE_BC')
            this%wat_rel_perm_func_owg => RPF_wat_owg_Burdine_BC_Create()
          ! case('TABLE')
          !   this%wat_rel_perm_func_owg => RPF_wat_owg_table_Create()
          !   call InputReadWord(input,option, &
          !                 this%wat_rel_perm_func_owg%table_name,PETSC_TRUE)
          !   call InputErrorMsg(input,option, &
          !                 'PERMEABILITY_FUNCTION_WAT/KRW,TABLE',error_string)
          case default
            call InputKeywordUnrecognized(input,word, &
                                      'PERMEABILITY_FUNCTION_WAT,KRW',option)
        end select
        call PermeabilityFunctionOWGRead(this%wat_rel_perm_func_owg, &
                                                   input,option)
      case('KRW_TABLE')
        this%wat_rel_perm_func_owg => RPF_wat_owg_table_Create()
        call InputReadWord(input,option, &
                      this%wat_rel_perm_func_owg%table_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'KRW_TABLE',error_string)        
      case('PERMEABILITY_FUNCTION_GAS','KRG')
        call InputReadCardDbaseCompatible(input,option,word)
        call InputErrorMsg(input,option,'PERMEABILITY_FUNCTION_GAS - KRG ', &
                           error_string)
        call StringToUpper(word)
        select case(word)
          case('MOD_BROOKS_COREY')
            this%gas_rel_perm_func_owg => RPF_gas_owg_MBC_Create()
          case('MUALEM_VG_SL')
            this%gas_rel_perm_func_owg => RPF_gas_owg_Mualem_VG_Create()
          case('BURDINE_VG_SL')
            this%gas_rel_perm_func_owg => RPF_gas_owg_Burdine_VG_Create()
          case('BURDINE_BC_SL')
            this%gas_rel_perm_func_owg => RPF_gas_owg_Burdine_BC_Create()
          case('TOUGH2_IRP7_SL')
            this%gas_rel_perm_func_owg => RPF_gas_owg_TOUGH2_IRP7_Create()
          ! case('TABLE')
          !   this%gas_rel_perm_func_owg => RPF_gas_owg_table_Create()
          !   call InputReadWord(input,option, &
          !                 this%gas_rel_perm_func_owg%table_name,PETSC_TRUE)
          !   call InputErrorMsg(input,option, &
          !                 'PERMEABILITY_FUNCTION_GAS/KRG,TABLE',error_string)
          case default
            call InputKeywordUnrecognized(input,word, &
                                  'PERMEABILITY_FUNCTION_GAS,KRG',option)
        end select
        call PermeabilityFunctionOWGRead(this%gas_rel_perm_func_owg, &
                                                   input,option)
      case('KRG_TABLE')
        this%gas_rel_perm_func_owg => RPF_gas_owg_table_Create()
        call InputReadWord(input,option, &
                      this%gas_rel_perm_func_owg%table_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'KRG_TABLE',error_string)
      case('PERMEABILITY_FUNCTION_OW','KROW', &
            'PERMEABILITY_FUNCTION_HC','KRH')
        !krh uses same class than krow    
        call InputReadCardDbaseCompatible(input,option,word)
        call InputErrorMsg(input,option,'PERMEABILITY_FUNCTION_OW/KROW/KRH ',&
                           error_string)
        call StringToUpper(word)
        select case(word)
          case('MOD_BROOKS_COREY')
            this%ow_rel_perm_func_owg => RPF_ow_owg_MBC_Create()
          case('TOUGH2_LINEAR')
            this%ow_rel_perm_func_owg => RPF_ow_owg_linear_Create()
          ! case('TABLE')
          !   this%ow_rel_perm_func_owg => RPF_ow_owg_table_Create()
          !   call InputReadWord(input,option, &
          !                 this%ow_rel_perm_func_owg%table_name,PETSC_TRUE)
          !   call InputErrorMsg(input,option, &
          !               'PERMEABILITY_FUNCTION_OW/KOW/KRH,TABLE',error_string)
          case default
            call InputKeywordUnrecognized(input,word, &
                              'PERMEABILITY_FUNCTION_OW/KROW/KRH',option)
        end select
        select case(word)
          case('PERMEABILITY_FUNCTION_HYDROCARBON','KRH')
            this%ow_rel_perm_func_owg%So_is_Sh = PETSC_TRUE
        end select
        call PermeabilityFunctionOWGRead(this%ow_rel_perm_func_owg, &
                                                   input,option)
      case('KROW_TABLE','KRH_TABLE')
        this%ow_rel_perm_func_owg => RPF_ow_owg_table_Create()
        call InputReadWord(input,option, &
                      this%ow_rel_perm_func_owg%table_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'KROW_TABLE/KRH_TABLE',error_string)        
      case('PERMEABILITY_FUNCTION_OG','KROG')
        call InputReadCardDbaseCompatible(input,option,word)
        call InputErrorMsg(input,option,'PERMEABILITY_FUNCTION_OG - KROG ', &
                           error_string)
        call StringToUpper(word)
        select case(word)
          case('MOD_BROOKS_COREY')
            this%og_rel_perm_func_owg => RPF_og_owg_MBC_Create()
          ! case('TABLE')
          !   this%og_rel_perm_func_owg => RPF_og_owg_table_Create()
          !   call InputReadWord(input,option, &
          !                 this%og_rel_perm_func_owg%table_name,PETSC_TRUE)
          !   call InputErrorMsg(input,option, &
          !                'PERMEABILITY_FUNCTION_OG/KROG,TABLE',error_string)
          case default
            call InputKeywordUnrecognized(input,word, &
                                  'PERMEABILITY_FUNCTION_OG,KROG',option)
        end select
        call PermeabilityFunctionOWGRead(this%og_rel_perm_func_owg, &
                                                   input,option)
      case('KROG_TABLE')
        this%og_rel_perm_func_owg => RPF_og_owg_table_Create()
        call InputReadWord(input,option, &
                      this%og_rel_perm_func_owg%table_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'KROG_TABLE',error_string)        
      case('PERMEABILITY_FUNCTION_OIL','KRO')
        call InputReadCardDbaseCompatible(input,option,word)
        call InputErrorMsg(input,option,'PERMEABILITY_FUNCTION_OIL - KRO ', &
                           error_string)
        call StringToUpper(word)
        select case(word)
        case('ECLIPSE')
            this%oil_rel_perm_func_owg => RPF_oil_ecl_Create()
          case default
            call InputKeywordUnrecognized(input,word, &
                                  'PERMEABILITY_FUNCTION_OIL,KRO',option)
        end select
      case('TEST')
        this%test = PETSC_TRUE
      case('DEFAULT')
        this%saturation_function => SF_Default_Create()
        this%liq_rel_perm_function => RPF_Default_Create()
        this%gas_rel_perm_function => this%liq_rel_perm_function
        ! PO TODO: adds default for OWG functions
      case default
        call InputKeywordUnrecognized(input,keyword,'CHARACTERISTIC_CURVES', &
                                      option)
    end select 
  enddo
  call InputPopBlock(input,option)
  
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
  end select
  
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadCard(input,option,keyword)
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
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputKeywordUnrecognized(input,keyword, &
                   'Brooks-Corey saturation function',option)
        end select
    !------------------------------------------
      class is(sat_func_Linear_type)
        select case(keyword)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
              call InputKeywordUnrecognized(input,keyword, &
                   'MODIFIED_KOSUGI saturation function',option)
          end select
    !------------------------------------------
      class default
        option%io_buffer = 'Read routine not implemented for ' &
                           // trim(error_string) // '.'
        call PrintErrMsg(option)
    !------------------------------------------
    end select
  enddo
  call InputPopBlock(input,option)
  
  if (smooth) then
    call saturation_function%SetupPolynomials(option,error_string)
  endif

  select type(sf => saturation_function)
  !------------------------------------------
    class is(sat_func_constant_type)
      option%io_buffer = 'Constant saturation function is being used.'
      call PrintWrnMsg(option)
  !------------------------------------------
    class is(sat_func_VG_type)
  !------------------------------------------
    class is(sat_func_BC_type)
      if (.not.smooth) then
        option%io_buffer = 'Brooks-Corey saturation function is being used &
          &without SMOOTH option.'
        call PrintWrnMsg(option)
      endif
  !------------------------------------------
    class is(sat_func_Linear_type)
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
    class is(rpf_mK_liq_type)
      error_string = trim(error_string) // 'MODIFIED_KOSUGI_LIQ'
    class is(rpf_mK_gas_type)
      error_string = trim(error_string) // 'MODIFIED_KOSUGI_GAS'
    class is(rel_perm_func_constant_type)
      error_string = trim(error_string) // 'CONSTANT'
  end select

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   

    ! base 
    found = PETSC_TRUE
    select case(keyword)
      case('LIQUID_RESIDUAL_SATURATION') 
        call InputReadDouble(input,option,permeability_function%Sr)
        call InputErrorMsg(input,option,'residual_saturation',error_string)
      case('PHASE')
        call InputReadCard(input,option,new_phase_keyword,PETSC_FALSE)
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
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputKeywordUnrecognized(input,keyword, &
              'Burdine Brooks-Corey gas relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Mualem_BC_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputKeywordUnrecognized(input,keyword, &
              'Burdine van Genuchten gas relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Mualem_Linear_liq_type)
        select case(keyword)
          case('MAX_CAPILLARY_PRESSURE') 
            call InputReadDouble(input,option,rpf%pcmax)
            call InputErrorMsg(input,option,'max_capillary_pressure', &
                               error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,rpf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputErrorMsg(input,option,'max_capillary_pressure', &
                               error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,rpf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'Mualem Linear gas relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Burdine_Linear_liq_type)
        select case(keyword)
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputKeywordUnrecognized(input,keyword, &
              'Burdine Linear gas relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_mK_liq_type)
        select case(keyword)
          case('SIGMAZ')
            call InputReadDouble(input,option,rpf%sigmaz)
            call InputErrorMsg(input,option,'sigmaz',error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
            call InputKeywordUnrecognized(input,keyword, &
                 'MODIFIED_KOSUGI gas relative permeability '//&
                 &'function',option)
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
            call InputKeywordUnrecognized(input,keyword, &
              'Constant relative permeability function', &
              option)
        end select
    !------------------------------------------
      class default
        option%io_buffer = 'Read routine not implemented for relative ' // &
                           'permeability function class.'
        call PrintErrMsg(option)
    !------------------------------------------
    end select
  enddo
  call InputPopBlock(input,option)
  
  ! for functions that are not phase-specific, check if PHASE was given:
  if (StringCompare('NONE',phase_keyword)) then
    phase_keyword = new_phase_keyword
    ! a liq or gas phase should now be specified for the non-phase-specific
    ! functions, so check if it was:
    if (StringCompare('NONE',new_phase_keyword)) then
      ! entering means the new phase keyword was also NONE (the default), so
      ! throw an error and abort:
      option%io_buffer = 'PHASE is not specified for ' // trim(error_string) 
      call PrintErrMsg(option)
    endif
  endif
  
  ! liquid phase relative permeability function check:
  if (StringCompare('LIQUID',phase_keyword)) then
    if (StringCompare('GAS',new_phase_keyword)) then
      ! user is requesting a liquid relative perm func for a gas phase:
      option%io_buffer = 'A liquid-phase relative permeability function &
                         &is being requested for the gas phase under ' &
                         // trim(error_string) // '.'
      call PrintErrMsg(option)
    endif
  endif
  
  ! gas phase relative permeability function check:
  if (StringCompare('GAS',phase_keyword)) then
    if (StringCompare('LIQUID',new_phase_keyword)) then
      ! user is requesting a gas relative perm func for a liquid phase:
      option%io_buffer = 'A gas-phase relative permeability function &
                         &is being requested for the liquid phase under ' &
                         // trim(error_string) // '.'
      call PrintErrMsg(option)
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
    if (cur_characteristic_curves%test) then
      call OptionSetBlocking(option,PETSC_FALSE)
      if (option%myrank == option%io_rank) then
        call CharacteristicCurvesTest(cur_characteristic_curves,option)
      endif
      call OptionSetBlocking(option,PETSC_TRUE)
      call OptionCheckNonBlockingError(option)
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
      class is(rpf_Mualem_Linear_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_Mualem_Linear_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_Burdine_Linear_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_Burdine_Linear_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rpf_mK_liq_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rpf_mK_gas_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Srg
      class is(rel_perm_func_constant_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class is(rel_perm_func_default_type)
        CharCurvesGetGetResidualSats(option%gas_phase) = rpf%Sr
      class default
        option%io_buffer = 'Relative permeability class not supported in &
              &CharCurvesGetGetResidualSats.'
        call PrintErrMsgToDev('',option)
    end select

  end if

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
  call PrintErrMsg(option)

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
    call PrintErrMsg(option)
  end if

  if (associated(characteristic_curves%gas_rel_perm_function) ) then
    call characteristic_curves%gas_rel_perm_function%Verify(string,option)
  else
    if (option%iflowmode == G_MODE .or. option%iflowmode == TOWG_MODE .or. &
        option%iflowmode == WF_MODE .or. option%iflowmode == H_MODE) then
      option%io_buffer = 'A gas phase relative permeability function has &
                         &not been set under CHARACTERISTIC_CURVES "' // &
                         trim(characteristic_curves%name) // '". Another &
                         &PERMEABILITY_FUNCTION block must be specified &
                         &for the gas phase.'
      call PrintErrMsg(option)
    end if
  end if

  
end subroutine CharacteristicCurvesVerify

! **************************************************************************** !

subroutine CharacteristicCurvesOWGVerify(characteristic_curves,option)
  !
  ! Checks if required parameters have been set for each curve type.
  ! Copy end points between classes where needed
  !
  ! Author: Paolo Orsini
  ! Date: 11/16/17 - 08/06/2018
  !
  use Option_module

  implicit none

  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: gas_present, oil_gas_interface_present 
  PetscBool :: wat_gas_interface_present
  PetscBool :: oil_perm_2ph_ow, oil_perm_3ph_owg
  PetscReal :: swcr, swco, sgco, sgcr, sowcr, sogcr
 
  swcr = 0.0
  swco = 0.0
  sgco = 0.0
  sgcr = 0.0
  sowcr = 0.0
  sogcr = 0.0
 
  call SetCCOWGPhaseFlags(option,oil_gas_interface_present, &
                          wat_gas_interface_present, gas_present, &
                          oil_perm_2ph_ow,oil_perm_3ph_owg)

  string = 'CHARACTERISTIC_CURVES(' // trim(characteristic_curves%name) // &
           '),'

  !check water rel perm - always
  if (.not.(associated(characteristic_curves%wat_rel_perm_func_owg))) then
    option%io_buffer = "A water relative permeability function has &
                        &not been set under CHARACTERISTIC_CURVES " // &
                        trim(characteristic_curves%name) // '". A &
                        &PERMEABILITY_FUNCTION_WAT/KRW block must be &
                        &specified for the water phase.'
    call PrintErrMsg(option)
  else
    !to avoid that verify for RPF_wat_MBC fails set sgcr and socr to zero
    !real values assigned later once RPF_gas and RPF_oil have been verirified
    select type(rpf => characteristic_curves%wat_rel_perm_func_owg)
      class is(RPF_wat_owg_MBC_type)
        call rpf%RPF_wat_owg_MBC_SetSowcr(0.d0,option)
    end select
    call characteristic_curves%wat_rel_perm_func_owg%verify(string,option)
    swco = characteristic_curves%wat_rel_perm_func_owg%GetConnateSaturation( &
                                                                      option)
    swcr = characteristic_curves%wat_rel_perm_func_owg%GetCriticalSaturation( &
                                                                       option)                                                                      
  end if  
  
  !check capillary pressures between oil and water - always
  if (.not.(associated(characteristic_curves%oil_wat_sat_func))) then
    option%io_buffer = "A water/oil saturation function has  &
                        &not been set under CHARACTERISTIC_CURVES " // &
                        trim(characteristic_curves%name) // '". A &
                        &CAP_PRESSURE_FUNCTION_OW/PC_OW block must be &
                        &specified for the water/oil phase interface.'
    call PrintErrMsg(option)
  else
    select type (sf => characteristic_curves%oil_wat_sat_func)
      class is(sat_func_xw_constant_type)
        call sf%SetConnateSaturation(swco,option)
        call sf%SetCriticalSaturation(swcr,option)
      class is(sat_func_xw_VG_type)
        call sf%SetConnateSaturation(swco,option)
      class is(sat_func_xw_table_type)
        !check if sat_func_of_pc_available
        if (sf%table%pc_inverse_available) then
          sf%sat_func_of_pc_available = PETSC_TRUE
        end if  
    end select
    call characteristic_curves%oil_wat_sat_func%verify(string,option)
    !Pm min is computed after function verification
    call characteristic_curves%oil_wat_sat_func%ComputePcMin(option)
    !check if swco and swcr defined in krw and pcw are the same
    if (characteristic_curves%oil_wat_sat_func%GetConnateSaturation(option) &
        /= swco ) then
      option%io_buffer = adjustl(trim(string)) // & 
                         'Swco defined in KRW and PC_XW differs- check input'
      call PrintErrMsg(option)
    end if
    if (characteristic_curves%oil_wat_sat_func%GetCriticalSaturation(option) &
        /= swcr ) then
      option%io_buffer = adjustl(trim(string)) // & 
                        'Swcr defined in KRW and PC_XW differs- check input'
      call PrintErrMsg(option)
    end if
  end if

  !check gas rel perm
  if (gas_present) then
    if (.not.(associated(characteristic_curves%gas_rel_perm_func_owg)) ) then
      option%io_buffer = "A gas relative permeability function has &
                          &not been set under CHARACTERISTIC_CURVES " // &
                          trim(characteristic_curves%name) // '". A &
                          &PERMEABILITY_FUNCTION_GAS/KRG block must be &
                          &specified for the gas phase.'
      call PrintErrMsg(option)
    else
      !to avoid that verify for RPF_gas_MBC fails set swcr and socr to zero
      !real values assigned later once RPF_wat and RPF_oil have been verirified      
      select type(rpf => characteristic_curves%gas_rel_perm_func_owg)
        class is(RPF_gas_owg_MBC_type)
          call rpf%RPF_gas_owg_MBC_SetSwcoSogcr(0.d0,0.d0,option)
      end select  
      call characteristic_curves%gas_rel_perm_func_owg%verify(string,option)
      sgco = &
           characteristic_curves%gas_rel_perm_func_owg%GetConnateSaturation( &
                                                                      option)
      sgcr = &
           characteristic_curves%gas_rel_perm_func_owg%GetCriticalSaturation( &
                                                                      option)
    end if
  end if ! end check gas rel perm

  !check oil/gas capillary pressure (Pcog)
  if (oil_gas_interface_present .and. gas_present) then
    if (.not.(associated(characteristic_curves%oil_gas_sat_func))) then
      option%io_buffer = "An oil/gas saturation function has  &
                          &not been set under CHARACTERISTIC_CURVES " // &
                          trim(characteristic_curves%name) // '". A &
                          &CAP_PRESSURE_FUNCTION_OG/PC_OG block must be &
                          &specified for the oil/gas phase interface.'
      call PrintErrMsg(option)
    else
      select type(sf => characteristic_curves%oil_gas_sat_func )
        class is(sat_func_og_constant_type)
          call sf%SetConnateSaturation(sgco,option)
          call sf%SetCriticalSaturation(sgcr,option)
        class is(sat_func_og_VG_SL_type)  
          call sf%SetConnateSaturation(sgco,option)
          call sf%SetCriticalSaturation(sgcr,option)
        class is(sat_func_og_table_type)
          !check if sat_func_of_pc_available
          if (sf%table%pc_inverse_available) then
            sf%sat_func_of_pc_available = PETSC_TRUE
          end if          
      end select
      call characteristic_curves%oil_gas_sat_func%verify(string,option)
      !Pcmin compute after fucntion vericifation
      call characteristic_curves%oil_gas_sat_func%ComputePcMin(option)
      !check if sgco and sgcr defined in KRG and PC_OG have the same values
      if (characteristic_curves%oil_gas_sat_func%GetConnateSaturation(option) &
          /= sgco ) then
        option%io_buffer = adjustl(trim(string)) // &
                          'Sgco in KRG and PC_OG differs - check input'
        call PrintErrMsg(option)
      end if
      if (characteristic_curves%oil_gas_sat_func%GetCriticalSaturation( &
          option) /= sgcr ) then
        option%io_buffer = adjustl(trim(string)) // &
                           'Sgcr in KRG and PC_OG differs - check input'
        call PrintErrMsg(option)
      end if    
    end if    
  end if  ! end check oil/gas capillary pressure (Pcog)

  if (oil_perm_2ph_ow) then
    if (.not.(associated(characteristic_curves%ow_rel_perm_func_owg)) ) then
      option%io_buffer = "A relative permeability function for oil in water &
                          &has not been set under CHARACTERISTIC_CURVES " // &
                          trim(characteristic_curves%name) // '". A &
                          &PERMEABILITY_FUNCTION_OW/KROW/KRH block must be &
                          &specified for the oil phase in water.'
      call PrintErrMsg(option)
    else
      select type(rpf => characteristic_curves%ow_rel_perm_func_owg)
        class is(rel_perm_ow_owg_MBC_type)
          call rpf%RPF_ow_owg_MBC_SetSwcr(swcr,option)
      end select
      call characteristic_curves%ow_rel_perm_func_owg%Verify(string,option)
      sowcr = characteristic_curves%ow_rel_perm_func_owg% &
                                              GetCriticalSaturation(option)
    end if
    if ( associated(characteristic_curves%oil_rel_perm_func_owg) ) then
      option%io_buffer = "Three phase (OWG) oil relative perability  &
                          &defined in CHARACTERISTIC_CURVES " // &
                          trim(characteristic_curves%name) // '". A &
                          &This is not supported in TOIL and TOWG_IMMISCIBLE'
      call PrintErrMsg(option)
    end if    
  end if !end oil phase check

  if (oil_perm_3ph_owg) then
    if (.not.(associated(characteristic_curves%oil_rel_perm_func_owg)) ) then
      option%io_buffer = "A oil relative permeability function has &
                          &not been set under CHARACTERISTIC_CURVES " // &
                          trim(characteristic_curves%name) // '". A &
                          &PERMEABILITY_FUNCTION_OIL/KRO block must be &
                          &specified for the oil phase.'
      call PrintErrMsg(option)
    else
      select type(rpf => characteristic_curves%oil_rel_perm_func_owg)
        class is(rel_perm_oil_owg_ecl_type)
          call rpf%RPF_oil_ecl_SetSwco(swco,option) 
          !check if KROW is defined and assign swcr if KROW_MBC
          if (.not.(associated(rpf%rel_perm_ow )) ) then
             option%io_buffer = "Within the oil Eclipse model krow has not &
                                 &been set under CHARACTERISTIC_CURVES " // &
                                 trim(characteristic_curves%name) // '". A &
                                 &PERMEABILITY_FUNCTION_OW/KROW block must be &
                                 &specified - KRO:ECLIPSE:KROW '         
          else
            select type (sub_rpf => rpf%rel_perm_ow)
              class is(rel_perm_ow_owg_MBC_type)
                call sub_rpf%RPF_ow_owg_MBC_SetSwcr(swcr,option)
            end select  
          end if
          !check if KROG is defined and assign swco and sgcr if KROG_MBC
          if (.not.(associated(rpf%rel_perm_og )) ) then
            option%io_buffer = "Within the oil Eclipse model krow has not &
                                &been set under CHARACTERISTIC_CURVES " // &
                                trim(characteristic_curves%name) // '". A &
                                &PERMEABILITY_FUNCTION_OG/KROG block must be &
                                &specified - KRO:ECLIPSE:KROG '         
          else
            select type (sub_rpf => rpf%rel_perm_og)
              class is(rel_perm_og_owg_MBC_type)
                call sub_rpf%RPF_og_owg_MBC_SetSwcoSgcr(swco,sgcr,option)
            end select  
          end if  
      end select
      call characteristic_curves%oil_rel_perm_func_owg%Verify(string,option)
      sowcr = characteristic_curves%oil_rel_perm_func_owg%GetSowcr(option)
      sogcr = characteristic_curves%oil_rel_perm_func_owg%GetSogcr(option)
    end if ! end if oil_rel_perm_func_owg
    if ( associated(characteristic_curves%ow_rel_perm_func_owg) ) then
      option%io_buffer = "KROW (Oil relative permeability in water)  &
                    &defined in CHARACTERISTIC_CURVES " // &
                    trim(characteristic_curves%name) // '". This is not &
                    &supported in TOWG:Black Oil,SOLVENT. KROW must be &
                    &defined within KRO'
      call PrintErrMsg(option)
    end if
    if ( associated(characteristic_curves%og_rel_perm_func_owg) ) then
      option%io_buffer = "KROG (Oil relative permeability in water) &
                    &defined in CHARACTERISTIC_CURVES " // &
                    trim(characteristic_curves%name) // '". This is not &
                    &supported in TOWG:Black Oil,SOLVENT. KROG must be &
                    &defined within KRO'
      call PrintErrMsg(option)
    end if        
  end if !end if oil_perm_3ph_owg

  !setup sgcr and socr for RPF_wat_MBC
  select type(rpf => characteristic_curves%wat_rel_perm_func_owg)
    class is(RPF_wat_owg_MBC_type)
      call rpf%RPF_wat_owg_MBC_SetSowcr(sowcr,option)
      call rpf%Verify(string,option)
  end select  

  !setup swcr and socr for RPF_gas_MBC
  if (gas_present) then
    select type(rpf => characteristic_curves%gas_rel_perm_func_owg)
      class is(RPF_gas_owg_MBC_type)
        call rpf%RPF_gas_owg_MBC_SetSwcoSogcr(swco,sogcr,option)
        call rpf%Verify(string,option)
    end select      
  end if


end subroutine CharacteristicCurvesOWGVerify

! **************************************************************************** !

subroutine CharCurvesOWGPostReadProcess(cc,option)
  !
  ! Process CharCurves:
  ! - associate cc curves with cc_tables
  ! - add here any other cc_curves post-read processing 
  !
  ! Author: Paolo Orsini
  ! Date: 08/06/18
  !
  use Option_module

  implicit none

  class(characteristic_curves_type) :: cc
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXSTRINGLENGTH) :: error_string_search
  character(len=MAXWORDLENGTH) :: table_name
  PetscBool :: gas_present, oil_gas_interface_present 
  PetscBool :: wat_gas_interface_present
  PetscBool :: oil_perm_2ph_ow, oil_perm_3ph_owg

  call SetCCOWGPhaseFlags(option,oil_gas_interface_present, &
                          wat_gas_interface_present, gas_present, &
                          oil_perm_2ph_ow,oil_perm_3ph_owg)

  error_string = 'CHARACTERISTIC CURVES,(' // trim(cc%name) // '),'
  error_string_search = ''
 
  if (associated(cc%oil_wat_sat_func)) then
    select type (sf => cc%oil_wat_sat_func)
      class is (sat_func_xw_table_type)
        call sf%ProcessTable(cc%char_curves_tables,error_string,option)
    end select
  else !attempt to create from list of cc_tables - always
    error_string_search = trim(error_string) // 'searching for PC_OW,'
    call SearchCCTVarInCCTableList(cc%char_curves_tables,CCT_PCXW, &
                                    table_name,error_string_search,option)
    cc%oil_wat_sat_func => SF_XW_table_Create()
    cc%oil_wat_sat_func%table_name = table_name
    call cc%oil_wat_sat_func%ProcessTable(cc%char_curves_tables, &
                                                error_string_search,option)
  end if
  
  if (associated(cc%gas_wat_sat_func)) then 
    select type (sf => cc%gas_wat_sat_func)
      class is (sat_func_xw_table_type)
        call sf%ProcessTable(cc%char_curves_tables,error_string,option)
    end select
  else if (wat_gas_interface_present) then !attempt to create from list of cc_tables
    error_string_search = trim(error_string) // 'searching for PC_GW,'
    call SearchCCTVarInCCTableList(cc%char_curves_tables,CCT_PCXW, &
                                    table_name,error_string_search,option)
    cc%gas_wat_sat_func => SF_XW_table_Create()
    cc%gas_wat_sat_func%table_name = table_name
    call cc%gas_wat_sat_func%ProcessTable(cc%char_curves_tables, &
                                                error_string_search,option)
  end if

  if (associated(cc%oil_gas_sat_func)) then 
    select type (sf => cc%oil_gas_sat_func)
      class is (sat_func_og_table_type)
        call sf%ProcessTable(cc%char_curves_tables,error_string,option)
    end select
  else if (oil_gas_interface_present) then !attempt to create from list of cc_tables
    error_string_search = trim(error_string) // 'searching for PC_OG,'
    call SearchCCTVarInCCTableList(cc%char_curves_tables,CCT_PCOG, &
                                     table_name,error_string_search,option)
    cc%oil_gas_sat_func => SF_OG_table_Create()
    cc%oil_gas_sat_func%table_name = table_name
    call cc%oil_gas_sat_func%ProcessTable(cc%char_curves_tables, &
                                                  error_string_search,option)
  end if

  if (associated(cc%wat_rel_perm_func_owg)) then
    select type (rpf => cc%wat_rel_perm_func_owg)
      class is (RPF_wat_owg_table_type)
        call rpf%ProcessTable(cc%char_curves_tables,error_string,option)
    end select
  else !attempt to create from list of cc_tables - always
    error_string_search = trim(error_string) // 'searching for KRW,'
    call SearchCCTVarInCCTableList(cc%char_curves_tables,CCT_KRW, &
                                    table_name,error_string_search,option)
    cc%wat_rel_perm_func_owg => RPF_wat_owg_table_Create()
    cc%wat_rel_perm_func_owg%table_name = table_name
    call cc%wat_rel_perm_func_owg%ProcessTable(cc%char_curves_tables, &
                                                  error_string_search,option)
  end if

  if (associated(cc%gas_rel_perm_func_owg)) then
    select type (rpf => cc%gas_rel_perm_func_owg)
     class is (RPF_gas_owg_table_type)
       call rpf%ProcessTable(cc%char_curves_tables,error_string,option)
    end select
  else if(gas_present) then !attempt to create from list of cc_tables
    error_string_search = trim(error_string) // 'searching for KRG,'
    call SearchCCTVarInCCTableList(cc%char_curves_tables,CCT_KRG, &
                                    table_name,error_string_search,option)
    cc%gas_rel_perm_func_owg => RPF_gas_owg_table_Create()
    cc%gas_rel_perm_func_owg%table_name = table_name
    call cc%gas_rel_perm_func_owg%ProcessTable(cc%char_curves_tables, &
                                                  error_string_search,option)
  end if

  if (associated(cc%ow_rel_perm_func_owg)) then
    select type (rpf => cc%ow_rel_perm_func_owg)
      class is (rel_perm_ow_owg_table_type)
        call rpf%ProcessTable(cc%char_curves_tables,error_string,option)
    end select
  else if (oil_perm_2ph_ow) then !attempt to create from list of cc_tables
    error_string_search = trim(error_string) // 'searching for KROW,'
    call SearchCCTVarInCCTableList(cc%char_curves_tables,CCT_KROW, &
                                    table_name,error_string_search,option)
    cc%ow_rel_perm_func_owg => RPF_ow_owg_table_Create()
    cc%ow_rel_perm_func_owg%table_name = table_name
    call cc%ow_rel_perm_func_owg%ProcessTable(cc%char_curves_tables, &
                                               error_string_search,option)
  end if

  if (associated(cc%oil_rel_perm_func_owg)) then
    if (associated(cc%oil_rel_perm_func_owg%rel_perm_ow)) then
      select type (rpf => cc%oil_rel_perm_func_owg%rel_perm_ow)
       class is(rel_perm_ow_owg_table_type)
         call rpf%ProcessTable(cc%char_curves_tables,error_string,option)
      end select
    end if
    if (associated(cc%oil_rel_perm_func_owg%rel_perm_og)) then
      select type (rpf => cc%oil_rel_perm_func_owg%rel_perm_og)
        class is(rel_perm_og_owg_table_type)
          call rpf%ProcessTable(cc%char_curves_tables,error_string,option)
       end select
    end if   
  else if(oil_perm_3ph_owg) then
    !default to eclipse - user must enter the KRO block to define different 
    !models when available
    cc%oil_rel_perm_func_owg => RPF_oil_ecl_Create()
    !attempt to create KROW from list of cc_tables
    error_string_search = trim(error_string) // 'searching for KROW,'
    call SearchCCTVarInCCTableList(cc%char_curves_tables,CCT_KROW, &
                                   table_name,error_string_search,option)
    cc%oil_rel_perm_func_owg%rel_perm_ow => RPF_ow_owg_table_Create()
    cc%oil_rel_perm_func_owg%rel_perm_ow%table_name = table_name
    call cc%oil_rel_perm_func_owg%rel_perm_ow%ProcessTable( &
                        cc%char_curves_tables,error_string_search,option)
    !attempt to create KROG from list of cc_tables
    error_string_search = trim(error_string) // 'searching for KROG,'
    call SearchCCTVarInCCTableList(cc%char_curves_tables,CCT_KROG, &
                                  table_name,error_string_search,option)
    cc%oil_rel_perm_func_owg%rel_perm_og => RPF_og_owg_table_Create()
    cc%oil_rel_perm_func_owg%rel_perm_og%table_name = table_name
    call cc%oil_rel_perm_func_owg%rel_perm_og%ProcessTable( &
                         cc%char_curves_tables,error_string_search,option)
  end if

end subroutine CharCurvesOWGPostReadProcess

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

  ! the liquid and gas relative permeability pointers may pointer to the
  ! same address. if so, destroy one and nullify the other.
  if (associated(cc%liq_rel_perm_function,cc%gas_rel_perm_function)) then
    call PermeabilityFunctionDestroy(cc%liq_rel_perm_function)
    nullify(cc%gas_rel_perm_function)
  else
    call PermeabilityFunctionDestroy(cc%liq_rel_perm_function)
    call PermeabilityFunctionDestroy(cc%gas_rel_perm_function)
  endif

  deallocate(cc)
  nullify(cc)
  
end subroutine CharacteristicCurvesDestroy

! ************************************************************************** !

function CharacteristicCurvesCreateCopy(cc, option)
  !
  ! Creates a SINGLE (i.e. NULL next) characteristic curve object that copy input 'cc'
  !
  ! Date: 10/18/2017
  ! Author: Fengming Yuan

  use Option_module
  use String_module

  implicit none

  type(option_type) :: option

  class(characteristic_curves_type), pointer :: CharacteristicCurvesCreateCopy

  class(characteristic_curves_type), pointer, intent(in) :: cc
  class(characteristic_curves_type), pointer :: cc2

  !
  cc2 => CharacteristicCurvesCreate()

  cc2%name = trim(adjustl(cc%name))//'Copy'
  !
  if (associated(cc%saturation_function)) then
    select type(sf => cc%saturation_function)
      class is(sat_func_VG_type)
        cc2%saturation_function => SF_VG_Create()
      select type(sf2=>cc2%saturation_function)
        class is(sat_func_VG_type)
          sf2%alpha  = sf%alpha
          sf2%m      = sf%m
      end select
      !
      class is(sat_func_BC_type)
        cc2%saturation_function => SF_BC_Create()
        select type(sf2=>cc2%saturation_function)
          class is(sat_func_BC_type)
            sf2%alpha  = sf%alpha
            sf2%lambda = sf%lambda
        end select
      !
      class default
        option%io_buffer = 'Currently NOT yet supported duplicating of saturation function type' // &
              ' when copying Characteristic_Curves'
        call printErrMsg(option)
    end select

    cc2%saturation_function%Sr     = &
      cc%saturation_function%Sr
    cc2%saturation_function%pcmax  = &
      cc%saturation_function%pcmax
    cc2%saturation_function%analytical_derivative_available = &
      cc%saturation_function%analytical_derivative_available
    !
    if(associated(cc%saturation_function%sat_poly)) then
        cc2%saturation_function%sat_poly => PolynomialCreate()
        cc2%saturation_function%sat_poly%low = &
          cc%saturation_function%sat_poly%low
        cc2%saturation_function%sat_poly%high = &
          cc%saturation_function%sat_poly%high
        cc2%saturation_function%sat_poly%coefficients(1:4) = &
          cc%saturation_function%sat_poly%coefficients(1:4)
    endif

    if(associated(cc%saturation_function%pres_poly)) then
        cc2%saturation_function%pres_poly => PolynomialCreate()
        cc2%saturation_function%pres_poly%low = &
          cc%saturation_function%pres_poly%low
        cc2%saturation_function%pres_poly%high = &
          cc%saturation_function%pres_poly%high
        cc2%saturation_function%pres_poly%coefficients(1:4) = &
          cc%saturation_function%pres_poly%coefficients(1:4)
    endif
    !
#ifdef SMOOTHING2
    if(associated(cc%saturation_function%sat_poly2)) then
      cc2%saturation_function%sat_poly2 => PolynomialCreate()

      cc2%saturation_function%sat_poly2%low = &
        cc%saturation_function%sat_poly2%low
      cc2%saturation_function%sat_poly2%high = &
        cc%saturation_function%sat_poly2%high
      cc2%saturation_function%sat_poly2%coefficients(1:4) = &
        cc%saturation_function%sat_poly2%coefficients(1:4)
    endif
    if(associated(cc%saturation_function%pres_poly2)) then
      cc2%saturation_function%pres_poly2 => PolynomialCreate()

      cc2%saturation_function%pres_poly2%low = &
        cc%saturation_function%pres_poly2%low
      cc2%saturation_function%pres_poly2%high = &
        cc%saturation_function%pres_poly2%high
      cc2%saturation_function%pres_poly2%coefficients(1:4) = &
        cc%saturation_function%pres_poly2%coefficients(1:4)
    endif
    !
#endif
  end if

  ! LIQ-RPF
  if(associated(cc%liq_rel_perm_function)) then
    !
    select type(rpf => cc%liq_rel_perm_function)
      class is(rpf_Mualem_VG_liq_type)
        cc2%liq_rel_perm_function => RPF_Mualem_VG_Liq_Create()
        select type(rpf2 => cc2%liq_rel_perm_function)
          class is(rpf_Mualem_VG_liq_type)
            rpf2%m     = rpf%m
        end select
      !
      class is(rpf_Burdine_BC_liq_type)
        cc2%liq_rel_perm_function => RPF_Burdine_BC_Liq_Create()
        select type(rpf2 => cc2%liq_rel_perm_function)
          class is(rpf_Burdine_BC_liq_type)
            rpf2%lambda = rpf%lambda
        end select
      !
      class is(rpf_Mualem_BC_liq_type)
        cc2%liq_rel_perm_function => RPF_Mualem_BC_Liq_Create()
        select type(rpf2 => cc2%liq_rel_perm_function)
          class is(rpf_Mualem_BC_liq_type)
            rpf2%lambda     = rpf%lambda
        end select
      !
      class default
        option%io_buffer = 'Currently NOT yet supported liq. ' // &
             ' permissivity function type when copying.'
        call printErrMsg(option)
    end select

    cc2%liq_rel_perm_function%Sr  = &
      cc%liq_rel_perm_function%Sr
    cc2%liq_rel_perm_function%analytical_derivative_available = &
      cc%liq_rel_perm_function%analytical_derivative_available
    !
    if(associated(cc%liq_rel_perm_function%poly)) then
      cc2%liq_rel_perm_function%poly => PolynomialCreate()

      cc2%liq_rel_perm_function%poly%low = &
        cc%liq_rel_perm_function%poly%low
      cc2%liq_rel_perm_function%poly%high = &
        cc%liq_rel_perm_function%poly%high
      cc2%liq_rel_perm_function%poly%coefficients(1:4) = &
        cc%liq_rel_perm_function%poly%coefficients(1:4)
    endif
#ifdef SMOOTHING2
    if(associated(cc%liq_rel_perm_function%poly2)) then
      cc2%liq_rel_perm_function%poly2 => PolynomialCreate()

      cc2%liq_rel_perm_function%poly2%low = &
        cc%liq_rel_perm_function%poly2%low
      cc2%liq_rel_perm_function%poly2%high = &
        cc%liq_rel_perm_function%poly2%high
      cc2%liq_rel_perm_function%poly2%coefficients(1:4) = &
        cc%liq_rel_perm_function%poly2%coefficients(1:4)
    endif
#endif
  endif


  ! GAS-RPF
  if(associated(cc%gas_rel_perm_function)) then
    !
    select type(rpf => cc%gas_rel_perm_function)
      class is(rpf_Mualem_VG_gas_type)
        cc2%gas_rel_perm_function => RPF_Mualem_VG_Gas_Create()
        select type(rpf2 => cc2%gas_rel_perm_function)
          class is(rpf_Mualem_VG_gas_type)
            rpf2%m       = rpf%m
            rpf2%Srg     = rpf%Srg
        end select
      !
      class is(rpf_Burdine_BC_gas_type)
        cc2%gas_rel_perm_function => RPF_Burdine_BC_Gas_Create()
        select type(rpf2 => cc2%gas_rel_perm_function)
          class is(rpf_Burdine_BC_gas_type)
            rpf2%lambda  = rpf%lambda
            rpf2%Srg     = rpf%Srg
        end select
      !
      class is(rpf_Mualem_BC_gas_type)
        cc2%gas_rel_perm_function => RPF_Mualem_BC_Gas_Create()
        select type(rpf2 => cc2%gas_rel_perm_function)
          class is(rpf_Mualem_BC_gas_type)
            rpf2%lambda  = rpf%lambda
            rpf2%Srg     = rpf%Srg
        end select
      !
      class default
        option%io_buffer = 'Currently NOT yet supported gas ' // &
             ' permissivity function type when copying.'
        call printErrMsg(option)
    end select

    cc2%gas_rel_perm_function%Sr  = &
      cc%gas_rel_perm_function%Sr
    cc2%gas_rel_perm_function%analytical_derivative_available = &
      cc%gas_rel_perm_function%analytical_derivative_available
    !
    if(associated(cc%gas_rel_perm_function%poly)) then
      cc2%gas_rel_perm_function%poly => PolynomialCreate()

      cc2%gas_rel_perm_function%poly%low = &
        cc%gas_rel_perm_function%poly%low
      cc2%gas_rel_perm_function%poly%high = &
        cc%gas_rel_perm_function%poly%high
      cc2%gas_rel_perm_function%poly%coefficients(1:4) = &
        cc%gas_rel_perm_function%poly%coefficients(1:4)
    endif
#ifdef SMOOTHING2
    if(associated(cc%gas_rel_perm_function%poly2)) then
      cc2%gas_rel_perm_function%poly2 => PolynomialCreate()

      cc2%gas_rel_perm_function%poly2%low = &
        cc%gas_rel_perm_function%poly2%low
      cc2%gas_rel_perm_function%poly2%high = &
        cc%gas_rel_perm_function%poly2%high
      cc2%gas_rel_perm_function%poly2%coefficients(1:4) = &
        cc%gas_rel_perm_function%poly2%coefficients(1:4)
    endif
#endif
  endif

  CharacteristicCurvesCreateCopy => cc2
end function CharacteristicCurvesCreateCopy
! ************************************************************************** !

end module Characteristic_Curves_module
