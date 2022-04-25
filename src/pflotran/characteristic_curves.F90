module Characteristic_Curves_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_Common_module
  use Characteristic_Curves_WIPP_module
  use Characteristic_Curves_loop_invariant_module
  use Characteristic_Curves_WIPP_Invariant_module
  use Characteristic_Curves_spline_module

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
  contains
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

  class(sat_func_base_type), pointer :: sf_swap
  class(rel_perm_func_base_type), pointer :: rpf_swap

  nullify(sf_swap)
  nullify(rpf_swap)
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
            this%saturation_function => SFConstantCreate()
          case('VAN_GENUCHTEN')
            this%saturation_function => SFVGCreate()
          case('BROOKS_COREY')
            this%saturation_function => SFBCCreate()
          case('LINEAR')
            this%saturation_function => SFLinearCreate()
          case('MODIFIED_KOSUGI')
            this%saturation_function => SFmKCreate()
          case('BRAGFLO_KRP1')
            this%saturation_function => SFKRP1Create()
          case('BRAGFLO_KRP2')
            this%saturation_function => SFKRP2Create()
          case('BRAGFLO_KRP3')
            this%saturation_function => SFKRP3Create()
          case('BRAGFLO_KRP4')
            this%saturation_function => SFKRP4Create()
          case('BRAGFLO_KRP5')
            this%saturation_function => SFKRP5Create()
          case('BRAGFLO_KRP8')
            this%saturation_function => SFKRP8Create()
          case('BRAGFLO_KRP9')
            this%saturation_function => SFKRP9Create()
          case('BRAGFLO_KRP11')
            this%saturation_function => SFKRP11Create()
          case('BRAGFLO_KRP12')
            this%saturation_function => SFKRP12Create()
          case('IGHCC2_COMP')
            this%saturation_function => SFIGHCC2CompCreate()
          case('LOOKUP_TABLE')
            this%saturation_function => SFTableCreate()
          case default
            call InputKeywordUnrecognized(input,word,'SATURATION_FUNCTION', &
                                          option)
        end select

        sf_swap => SaturationFunctionRead(this%saturation_function,input,option)
        ! If a constructor was used, reassign this%saturation_function to swap
        if (associated(sf_swap)) then
          deallocate(this%saturation_function)
          this%saturation_function => sf_swap
        end if
    !-----------------------------------------------------------------------
      case('SATURATION_FUNCTION_OWG')
        option%io_buffer = 'SATURATION_FUNCTION_OWG is not supported any more &
                           &in CHARACTERISTIC_CURVES. Please use either: &
                           &CAP_PRESSURE_FUNCTION_OW or PC_OW for Pcow; &
                           &CAP_PRESSURE_FUNCTION_WG or PC_WG for Pcwg or; &
                           &CAP_PRESSURE_FUNCTION_OG or PC_OG for Pcog'
        call PrintErrMsg(option)
      case('PERMEABILITY_FUNCTION')
        nullify(rel_perm_function_ptr)
        phase_keyword = 'NONE'
        call InputReadCardDbaseCompatible(input,option,word)
        call InputErrorMsg(input,option,'PERMEABILITY_FUNCTION',error_string)
        call StringToUpper(word)
        select case(word)
          case('MUALEM','MUALEM_VG_LIQ')
            rel_perm_function_ptr => RPFMualemVGLiqCreate()
            phase_keyword = 'LIQUID'
          case('MUALEM_VG_GAS')
            rel_perm_function_ptr => RPFMualemVGGasCreate()
            phase_keyword = 'GAS'
          case('BURDINE','BURDINE_BC_LIQ')
            rel_perm_function_ptr => RPFBurdineBCLiqCreate()
            phase_keyword = 'LIQUID'
          case('BURDINE_BC_GAS')
            rel_perm_function_ptr => RPFBurdineBCGasCreate()
            phase_keyword = 'GAS'
          case('TOUGH2_IRP7_LIQ')
            rel_perm_function_ptr => RPFMualemVGLiqCreate()
            phase_keyword = 'LIQUID'
          case('TOUGH2_IRP7_GAS')
            rel_perm_function_ptr => RPFTOUGH2IRP7GasCreate()
            phase_keyword = 'GAS'
          case('MUALEM_BC_LIQ')
            rel_perm_function_ptr => RPFMualemBCLiqCreate()
            phase_keyword = 'LIQUID'
          case('MUALEM_BC_GAS')
            rel_perm_function_ptr => RPFMualemBCGasCreate()
            phase_keyword = 'GAS'
          case('BURDINE_VG_LIQ')
            rel_perm_function_ptr => RPFBurdineVGLiqCreate()
            phase_keyword = 'LIQUID'
          case('BURDINE_VG_GAS')
            rel_perm_function_ptr => RPFBurdineVGGasCreate()
            phase_keyword = 'GAS'
          case('MUALEM_LINEAR_LIQ')
            rel_perm_function_ptr => RPFMualemLinearLiqCreate()
            phase_keyword = 'LIQUID'
          case('MUALEM_LINEAR_GAS')
            rel_perm_function_ptr => RPFMualemLinearGasCreate()
            phase_keyword = 'GAS'
          case('BURDINE_LINEAR_LIQ')
            rel_perm_function_ptr => RPFBurdineLinearLiqCreate()
            phase_keyword = 'LIQUID'
          case('BURDINE_LINEAR_GAS')
            rel_perm_function_ptr => RPFBurdineLinearGasCreate()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP1_LIQ')
            rel_perm_function_ptr => RPFKRP1LiqCreate()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP1_GAS')
            rel_perm_function_ptr => RPFKRP1GasCreate()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP2_LIQ')
            rel_perm_function_ptr => RPFKRP2LiqCreate()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP2_GAS')
            rel_perm_function_ptr => RPFKRP2GasCreate()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP3_LIQ')
            rel_perm_function_ptr => RPFKRP3LiqCreate()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP3_GAS')
            rel_perm_function_ptr => RPFKRP3GasCreate()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP4_LIQ')
            rel_perm_function_ptr => RPFKRP4LiqCreate()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP4_GAS')
            rel_perm_function_ptr => RPFKRP4GasCreate()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP5_LIQ')
            rel_perm_function_ptr => RPFKRP5LiqCreate()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP5_GAS')
            rel_perm_function_ptr => RPFKRP5GasCreate()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP8_LIQ')
            rel_perm_function_ptr => RPFKRP8LiqCreate()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP8_GAS')
            rel_perm_function_ptr => RPFKRP8GasCreate()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP9_LIQ')
            rel_perm_function_ptr => RPFKRP9LiqCreate()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP9_GAS')
            rel_perm_function_ptr => RPFKRP9GasCreate()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP11_LIQ')
            rel_perm_function_ptr => RPFKRP11LiqCreate()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP11_GAS')
            rel_perm_function_ptr => RPFKRP11GasCreate()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP12_LIQ')
            rel_perm_function_ptr => RPFKRP12LiqCreate()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP12_GAS')
            rel_perm_function_ptr => RPFKRP12GasCreate()
            phase_keyword = 'GAS'
          case('MODIFIED_KOSUGI_LIQ')
            rel_perm_function_ptr => RPFmKLiqCreate()
            phase_keyword = 'LIQUID'
          case('MODIFIED_KOSUGI_GAS')
            rel_perm_function_ptr => RPFmKGasCreate()
            phase_keyword = 'GAS'
          case('IGHCC2_COMP_LIQ')
            rel_perm_function_ptr => RPFIGHCC2CompLiqCreate()
            phase_keyword = 'LIQUID'
          case('IGHCC2_COMP_GAS')
            rel_perm_function_ptr => RPFIGHCC2CompGasCreate()
            phase_keyword = 'GAS'
          case('TABLE_LIQ')
            rel_perm_function_ptr => RPFTABLELiqCreate()
            phase_keyword = 'LIQUID'
          case('TABLE_GAS')
            rel_perm_function_ptr => RPFTABLEGasCreate()
            phase_keyword = 'GAS'
          case('CONSTANT')
            rel_perm_function_ptr => RPFConstantCreate()
            ! phase_keyword = 'NONE'
          case default
            call InputKeywordUnrecognized(input,word,'PERMEABILITY_FUNCTION', &
                                          option)
        end select

        rpf_swap => PermeabilityFunctionRead(rel_perm_function_ptr, &
                                             phase_keyword, input,option)
        ! If a constructor was used, redirect rel_perm_function_ptr to swap
        if (associated(rpf_swap)) then
          deallocate(rel_perm_function_ptr)
          rel_perm_function_ptr => rpf_swap
        end if

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
      case('TEST')
        this%test = PETSC_TRUE
      case('DEFAULT')
        this%saturation_function => SFDefaultCreate()
        this%liq_rel_perm_function => RPFDefaultCreate()
        this%gas_rel_perm_function => this%liq_rel_perm_function
        ! PO TODO: adds default for OWG functions
      case default
        call InputKeywordUnrecognized(input,keyword,'CHARACTERISTIC_CURVES', &
                                      option)
    end select
  enddo
  call InputPopBlock(input,option)

  call CharacteristicCurvesVerify(this,option)

end subroutine CharacteristicCurvesRead

! ************************************************************************** !

function SaturationFunctionRead(saturation_function,input,option) &
  result (sf_swap)
  !
  ! Reads in contents of a SATURATION_FUNCTION block
  !
  use Option_module
  use Input_Aux_module
  use String_module
  use Dataset_Ascii_class

  implicit none

  class(sat_func_base_type) :: saturation_function
  type(input_type), pointer :: input
  type(option_type) :: option
  class(sat_func_base_type), pointer :: sf_swap, sf_swap2

  character(len=MAXWORDLENGTH) :: keyword, internal_units
  character(len=MAXSTRINGLENGTH) :: error_string, table_name, temp_string
  PetscBool :: found
  PetscBool :: smooth
  PetscBool :: spline

  ! Lexicon of compiled parameters
  character(len=MAXWORDLENGTH) :: unsat_ext
  PetscBool :: loop_invariant, tension
  PetscInt :: vg_rpf_opt
  PetscReal :: alpha, m, Pcmax, Slj, Sr, Srg

  PetscInt :: wipp_krp, wipp_kpc
  PetscReal :: wipp_expon, wipp_pct_alpha, wipp_pct_expon
  PetscReal :: wipp_s_min, wipp_s_effmin
  PetscBool :: wipp_pct_ignore

  nullify(sf_swap)
  ! Default values for unspecified parameters
  loop_invariant = PETSC_FALSE
  tension = PETSC_FALSE
  unsat_ext = ''
  vg_rpf_opt = 1 ! Mualem. Burdine option in progress
  alpha = 0d0
  m = 0d0
  Pcmax = 1d9
  Slj = 0d0
  Sr = 0d0
  Srg = 0d0
  wipp_krp = 0
  wipp_kpc = 0
  wipp_expon = 0d0

  input%ierr = 0
  smooth = PETSC_FALSE
  spline = PETSC_FALSE
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
    class is(sat_func_IGHCC2_Comp_type)
      error_string = trim(error_string) // 'IGHCC2_COMP'
    class is (sat_func_Table_type)
      error_string = trim(error_string) // 'LOOKUP_TABLE'
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
        call InputReadDouble(input,option,Sr)
        saturation_function%Sr = Sr
        call InputErrorMsg(input,option,'LIQUID_RESIDUAL_SATURATION', &
                           error_string)
      case('MAX_CAPILLARY_PRESSURE')
        call InputReadDouble(input,option,Pcmax)
        saturation_function%Pcmax = Pcmax
        call InputErrorMsg(input,option,'MAX_CAPILLARY_PRESSURE', &
                            error_string)
      case('CALCULATE_INTERFACIAL_TENSION')
        tension = PETSC_TRUE
        saturation_function%calc_int_tension = PETSC_TRUE
      case('SMOOTH')
        smooth = PETSC_TRUE
      case('SPLINE')
        spline = PETSC_TRUE
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
            call InputReadDouble(input,option,m)
            sf%m = m
            call InputErrorMsg(input,option,'m',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,alpha)
            sf%alpha = alpha
            call InputErrorMsg(input,option,'alpha',error_string)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case('UNSATURATED_EXTENSION')
            call InputReadCard(input,option,unsat_ext)
            call InputErrorMsg(input,option,'unsaturated extension',error_string)
          case('LIQUID_JUNCTION_SATURATION')
            call InputReadDouble(input,option,Slj)
            call InputErrorMsg(input,option,'liquid junction saturation', &
                               error_string)
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
      class is(sat_func_KRP1_type)
        wipp_krp = 1
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case('KPC')
            call InputReadInt(input,option,sf%kpc)
            wipp_kpc = sf%kpc
            call InputErrorMsg(input,option,'KPC',error_string)
          case('M')
            call InputReadDouble(input,option,sf%m)
            wipp_expon = sf%m
            call InputErrorMsg(input,option,'M',error_string)
          case('PCT_A')
            call InputReadDouble(input,option,sf%pct_a)
            wipp_pct_alpha = sf%pct_a
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP')
            call InputReadDouble(input,option,sf%pct_exp)
            wipp_pct_expon = sf%pct_exp
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,sf%Srg)
            Srg = sf%Srg
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('IGNORE_PERMEABILITY')
            sf%ignore_permeability = PETSC_TRUE
            wipp_pct_ignore = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            alpha = sf%alpha
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case('LIQUID_JUNCTION_SATURATION')
            call InputReadDouble(input,option,Slj)
            call InputErrorMsg(input,option,'liquid junction saturation', &
                               error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP1',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP2_type)
        wipp_krp = 2
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case('KPC')
            call InputReadInt(input,option,sf%kpc)
            call InputErrorMsg(input,option,'KPC',error_string)
          case('LAMBDA')
            call InputReadDouble(input,option,sf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('PCT_A')
            call InputReadDouble(input,option,sf%pct_a)
            wipp_pct_alpha = sf%pct_a
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP')
            call InputReadDouble(input,option,sf%pct_exp)
            wipp_pct_expon = sf%pct_exp
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('IGNORE_PERMEABILITY')
            sf%ignore_permeability = PETSC_TRUE
            wipp_pct_ignore = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            alpha = sf%alpha
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case('LIQUID_JUNCTION_SATURATION')
            call InputReadDouble(input,option,Slj)
            call InputErrorMsg(input,option,'liquid junction saturation', &
                               error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP2',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP3_type)
        wipp_krp = 3
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case('KPC')
            call InputReadInt(input,option,sf%kpc)
            wipp_kpc = sf%kpc
            call InputErrorMsg(input,option,'KPC',error_string)
          case('LAMBDA')
            call InputReadDouble(input,option,sf%lambda)
            wipp_expon = sf%lambda
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('PCT_A')
            call InputReadDouble(input,option,sf%pct_a)
            wipp_pct_alpha = sf%pct_a
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP')
            call InputReadDouble(input,option,sf%pct_exp)
            wipp_pct_expon = sf%pct_exp
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,sf%Srg)
            Srg = sf%Srg
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('IGNORE_PERMEABILITY')
            sf%ignore_permeability = PETSC_TRUE
            wipp_pct_ignore = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            alpha = sf%alpha
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case('LIQUID_JUNCTION_SATURATION')
            call InputReadDouble(input,option,Slj)
            call InputErrorMsg(input,option,'liquid junction saturation', &
                               error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP3',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP4_type)
        wipp_krp = 4
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case('KPC')
            call InputReadInt(input,option,sf%kpc)
            wipp_kpc = sf%kpc
            call InputErrorMsg(input,option,'KPC',error_string)
          case('LAMBDA')
            call InputReadDouble(input,option,sf%lambda)
            wipp_expon = sf%lambda
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('PCT_A')
            call InputReadDouble(input,option,sf%pct_a)
            wipp_pct_alpha = sf%pct_a
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP')
            call InputReadDouble(input,option,sf%pct_exp)
            wipp_pct_expon = sf%pct_exp
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,sf%Srg)
            Srg = sf%Srg
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('IGNORE_PERMEABILITY')
            sf%ignore_permeability = PETSC_TRUE
            wipp_pct_ignore = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            alpha = sf%alpha
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case('LIQUID_JUNCTION_SATURATION')
            call InputReadDouble(input,option,Slj)
            call InputErrorMsg(input,option,'liquid junction saturation', &
                               error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP4',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP5_type)
        wipp_krp = 5
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case('KPC')
            call InputReadInt(input,option,sf%kpc)
            wipp_kpc = sf%kpc
            call InputErrorMsg(input,option,'KPC',error_string)
          case('PCT_A')
            call InputReadDouble(input,option,sf%pct_a)
            wipp_pct_alpha = sf%pct_a
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP')
            call InputReadDouble(input,option,sf%pct_exp)
            wipp_pct_expon = sf%pct_exp
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,sf%Srg)
            Srg = sf%Srg
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('IGNORE_PERMEABILITY')
            sf%ignore_permeability = PETSC_TRUE
            wipp_pct_ignore = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            alpha = sf%alpha
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case('LIQUID_JUNCTION_SATURATION')
            call InputReadDouble(input,option,Slj)
            call InputErrorMsg(input,option,'liquid junction saturation', &
                               error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP5',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP8_type)
        wipp_krp = 8
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case('KPC')
            call InputReadInt(input,option,sf%kpc)
            wipp_kpc = sf%kpc
            call InputErrorMsg(input,option,'KPC',error_string)
          case('M')
            call InputReadDouble(input,option,sf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case('PCT_A')
            call InputReadDouble(input,option,sf%pct_a)
            wipp_pct_alpha = sf%pct_a
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP')
            call InputReadDouble(input,option,sf%pct_exp)
            wipp_pct_expon = sf%pct_exp
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,sf%Srg)
            Srg = sf%Srg
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('IGNORE_PERMEABILITY')
            sf%ignore_permeability = PETSC_TRUE
            wipp_pct_ignore = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            alpha = sf%alpha
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case('LIQUID_JUNCTION_SATURATION')
            call InputReadDouble(input,option,Slj)
            call InputErrorMsg(input,option,'liquid junction saturation', &
                               error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP8',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP9_type)
        wipp_krp = 9
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP9',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP11_type)
        wipp_krp = 11
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP11',option)
        end select
    !------------------------------------------
      class is(sat_func_KRP12_type)
        wipp_krp = 12
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case('KPC')
            call InputReadInt(input,option,sf%kpc)
            wipp_kpc = sf%kpc
            call InputErrorMsg(input,option,'KPC',error_string)
          case('PCT_A')
            call InputReadDouble(input,option,sf%pct_a)
            wipp_pct_alpha = sf%pct_a
            call InputErrorMsg(input,option,'PCT_A',error_string)
          case('PCT_EXP')
            call InputReadDouble(input,option,sf%pct_exp)
            wipp_pct_expon = sf%pct_exp
            call InputErrorMsg(input,option,'PCT_EXP',error_string)
          case('LAMBDA')
            call InputReadDouble(input,option,sf%lambda)
            wipp_expon = sf%lambda
            call InputErrorMsg(input,option,'lambda',error_string)
          case('S_MIN')
            call InputReadDouble(input,option,sf%s_min)
            wipp_s_min = sf%s_min
            call InputErrorMsg(input,option,'s_min',error_string)
          case('S_EFFMIN')
            call InputReadDouble(input,option,sf%s_effmin)
            wipp_s_effmin = sf%s_effmin
            call InputErrorMsg(input,option,'s_effmin',error_string)
          case('IGNORE_PERMEABILITY')
            sf%ignore_permeability = PETSC_TRUE
            wipp_pct_ignore = PETSC_TRUE
            call InputErrorMsg(input,option,'IGNORE_PERMEABILITY',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            alpha = sf%alpha
            call InputErrorMsg(input,option,'ALPHA',error_string)
          case('LIQUID_JUNCTION_SATURATION')
            call InputReadDouble(input,option,Slj)
            call InputErrorMsg(input,option,'liquid junction saturation', &
                               error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'SATURATION_FUNCTION BRAGFLO_KRP12',option)
        end select
    !------------------------------------------
      class is(sat_func_IGHCC2_Comp_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,sf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'saturation function IGHCC2 Comparison',option)
        end select
      class is (sat_func_Table_type)
        select case(keyword)
          case('FILE')
            internal_units = 'unitless , Pa'
            call InputReadFilename(input,option,table_name)
            call DatasetAsciiReadFile(sf%pc_dataset,table_name, &
                                      temp_string, internal_units, &
                                      error_string,option)
        end select
      class default
        option%io_buffer = 'Read routine not implemented for ' &
                           // trim(error_string) // '.'
        call PrintErrMsg(option)
    !------------------------------------------
    end select
  enddo
  call InputPopBlock(input,option)

  ! At end of input block, call constructors if implemented
  ! Throw errors for invalid combinations of options or parametersa

  ! Error checking for wipp_pct_ignore option
  if (wipp_pct_ignore) then ! Check it is not overspecife
    if (alpha == 0d0) then
      option%io_buffer = 'Must specify ALPHA with IGNORE_PERMEABILITY option'
    else
      wipp_pct_alpha = alpha ! Copy to wipp_pct_alpha for constructor
    end if
  else
    if (alpha /= 0d0) then ! Error, pct_a must be specified
      option%io_buffer = 'CANNOT specify ALPHA without IGNORE_PERMEABILITY option'
    end if
  end if

  if (loop_invariant) then
    ! Use default junction saturation if not specified
    if (Slj == 0d0) Slj = Sr + 5d-2*(1d0-Srg-Sr)
    ! Call constructor
    if (wipp_krp /= 0) then ! WIPP invariants flagged by wipp_krp
      if (wipp_krp == 12) then ! wipp_s_min replaces Sr, wipp_s_effmin replaces Slj
        sf_swap => SFWIPPctor(wipp_krp, wipp_kpc, wipp_s_min, Srg, wipp_expon, &
                              wipp_pct_ignore, wipp_pct_alpha, wipp_pct_expon, &
                              Pcmax, wipp_s_effmin)
      else
        sf_swap => SFWIPPctor(wipp_krp, wipp_kpc, Sr, Srg, wipp_expon, &
                              wipp_pct_ignore, wipp_pct_alpha, wipp_pct_expon, &
                              Pcmax, Slj)
      end if
    else ! Old object type is used to identify common invariants
      select type (saturation_function)
      class is (sat_func_VG_type)
        call StringtoUpper(unsat_ext)
        sf_swap => SFVGctor(unsat_ext, alpha, m, Sr, vg_rpf_opt, Pcmax, Slj)
      class default
        option%io_buffer = 'Loop-invariant optimizations are not yet &
       & implemented for the designated saturation function type.'
        call PrintErrMsg(option)
      end select
    end if

    ! If successful, write tension option to the new object
    if (associated(sf_swap)) then
      sf_swap%calc_int_tension = tension
    else
    ! Throw an error the contructor failed. Most likley an invalid parameter
      option%io_buffer = 'Construction of the saturation function object &
      & failed.'
      call PrintErrMsg(option)
    end if
  else if (unsat_ext /= '') then
    ! Throw an error if unsaturated extensions are with loop_invariant
    option%io_buffer = 'Unsaturated extensions are unavailable without the &
    & loop-invariant optimization'
    call PrintErrMsg(option)
  end if


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
    class is(sat_func_WIPP_type)
      if (sf%ignore_permeability .and. Uninitialized(sf%alpha)) then
        option%io_buffer = 'If a WIPP capillary presure - saturation function &
          &is being used with the IGNORE_PERMEABILITY feature, you must &
          &specify ALPHA (inverse of the threshold capillary pressure). Do &
          &not specify PCT_A or PCT_EXP.'
        call PrintErrMsg(option)
      endif
  !------------------------------------------
  end select

  if (spline) then ! Create cubic approximation for any saturation function
    if (associated(sf_swap)) then ! Splining a loop-invariant replacement
      sf_swap2 => SFSplineCtor(sf_swap, 100)
      deallocate(sf_swap)
      sf_swap => sf_swap2
    else
      sf_swap => SFSplineCtor(saturation_function, 100)
    end if
    ! The calling CCRead will deallocated saturation_function
    sf_swap%calc_int_tension = tension
  end if

end function SaturationFunctionRead

! ************************************************************************** !

function PermeabilityFunctionRead(permeability_function,phase_keyword, &
                                    input,option) result (rpf_swap)
  !
  ! Reads in contents of a PERMEABILITY_FUNCTION block
  !
  use Option_module
  use Input_Aux_module
  use String_module
  use Dataset_Ascii_class

  implicit none

  class(rel_perm_func_base_type) :: permeability_function
  character(len=MAXWORDLENGTH) :: phase_keyword
  type(input_type), pointer :: input
  type(option_type) :: option
  class(rel_perm_func_base_type), pointer :: rpf_swap, rpf_swap2

  character(len=MAXWORDLENGTH) :: keyword, new_phase_keyword
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXSTRINGLENGTH) :: table_name, temp_string
  PetscBool :: found
  PetscBool :: smooth
  PetscBool :: spline

  ! Lexicon for compiled variables
  PetscBool :: loop_invariant
  PetscReal :: m, Srg, Sr

  nullify(rpf_swap)

  ! Default values for unspecified parameters
  spline = PETSC_FALSE
  loop_invariant = PETSC_FALSE
  m = 0d0
  Srg = 0d0
  Sr = 0d0

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
    class is(rpf_IGHCC2_Comp_liq_type)
      error_string = trim(error_string) // 'IGHCC2_COMP_LIQ'
    class is(rpf_IGHCC2_Comp_gas_type)
      error_string = trim(error_string) // 'IGHCC2_COMP_GAS'
    class is(rpf_Table_liq_type)
      error_string = trim(error_string) // 'LOOKUP_TABLE_LIQ'
    class is(rpf_Table_gas_type)
      error_string = trim(error_string) // 'LOOKUP_TABLE_GAS'
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
        call InputReadDouble(input,option,Sr)
        permeability_function%Sr = Sr
        call InputErrorMsg(input,option,'residual_saturation',error_string)
      case('PHASE')
        call InputReadCard(input,option,new_phase_keyword,PETSC_FALSE)
        call InputErrorMsg(input,option,'phase',error_string)
        call StringToUpper(phase_keyword)
      case('SMOOTH')
        smooth = PETSC_TRUE
      case('SPLINE')
        spline = PETSC_TRUE
      case default
        found = PETSC_FALSE
    end select
    if (found) cycle

    select type(rpf => permeability_function)
    !------------------------------------------
      class is(rpf_Mualem_VG_liq_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,m)
            rpf%m = m
            call InputErrorMsg(input,option,'m',error_string)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'Mualem van Genuchten liquid relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Mualem_VG_gas_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,m)
            rpf%m = m
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,Srg)
            rpf%Srg = Srg
            call InputErrorMsg(input,option,'Srg',error_string)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
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
      class is(rpf_TOUGH2_IRP7_gas_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,m)
            rpf%m = m
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
                   'TOUGH2 IRP7 gas relative permeability function',option)
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
            call InputReadDouble(input,option,m)
            rpf%m = m
            call InputErrorMsg(input,option,'m',error_string)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'Burdine van Genuchten liquid relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Burdine_VG_gas_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,m)
            rpf%m = m
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,Srg)
            rpf%Srg = Srg
            call InputErrorMsg(input,option,'Srg',error_string)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
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
      class is(rpf_KRP1_liq_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'BRAGFLO_KRP1_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP2_liq_type)
        select case(keyword)
          case('LAMBDA')
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'BRAGFLO_KRP2_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP2_gas_type)
        select case(keyword)
          case('LAMBDA')
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'BRAGFLO_KRP5_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP8_liq_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'BRAGFLO_KRP8_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP8_gas_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'BRAGFLO_KRP8_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP9_liq_type)
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'BRAGFLO_KRP9_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP9_gas_type)
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'BRAGFLO_KRP9_GAS relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP11_liq_type)
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case('TOLC')
            call InputReadDouble(input,option,rpf%tolc)
            call InputErrorMsg(input,option,'TOLC',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                 error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'BRAGFLO_KRP11_LIQ relative permeability function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_KRP11_gas_type)
        select case(keyword)
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case('TOLC')
            call InputReadDouble(input,option,rpf%tolc)
            call InputErrorMsg(input,option,'TOLC',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                 error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
          case('LOOP_INVARIANT')
            loop_invariant = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword, &
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
      class is(rpf_IGHCC2_Comp_liq_type)
        select case(keyword)
          case('LAMBDA')
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'IGHCC2 Comparison liq rel perm function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_IGHCC2_Comp_gas_type)
        select case(keyword)
          case('LAMBDA')
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'IGHCC2 Comparison gas rel perm function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Table_liq_type)
        select case(keyword)
          case('FILE')
            internal_units = 'unitless , unitless'
            call InputReadFilename(input,option,table_name)
            call DatasetAsciiReadFile(rpf%rpf_dataset,table_name, &
                                      temp_string, internal_units, &
                                      error_string,option)
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'Lookup Table liquid rel perm function', &
              option)
        end select
    !------------------------------------------
      class is(rpf_Table_gas_type)
        select case(keyword)
          case('FILE')
            internal_units = 'unitless , unitless'
            call InputReadFilename(input,option,table_name)
            call DatasetAsciiReadFile(rpf%rpf_dataset,table_name, &
                                      temp_string, internal_units, &
                                      error_string,option)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(input,keyword, &
              'Lookup Table gas rel perm function', &
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

  ! At the end of the input block, call constructors as applicable
  ! to replace with optimized relative permeability functions
  if (loop_invariant) then
    select type (rpf => permeability_function)
    class is (RPF_mualem_VG_liq_type)
      rpf_swap => RPFMVGliqCtor(m, Sr)
    class is (RPF_burdine_VG_liq_type)
      rpf_swap => RPFBVGliqCtor(m, Sr)
    class is (RPF_mualem_VG_gas_type)
      rpf_swap => RPFMVGgasCtor(m, Sr, Srg)
    class is (RPF_burdine_VG_gas_type)
      rpf_swap => RPFBVGgasCtor(m, Sr, Srg)
    class is (rpf_KRP1_liq_type)
      rpf_swap => RPFWIPPctor(PETSC_TRUE,1,rpf%Sr,rpf%Srg,rpf%m)
    class is (rpf_KRP1_gas_type)
      rpf_swap => RPFWIPPctor(PETSC_FALSE,1,rpf%Sr,rpf%Srg,rpf%m)
    class is (rpf_KRP2_liq_type)
      rpf_swap => RPFWIPPctor(PETSC_TRUE,2,rpf%Sr,rpf%Srg,rpf%lambda)
    class is (rpf_KRP2_gas_type)
      rpf_swap => RPFWIPPctor(PETSC_FALSE,2,rpf%Sr,rpf%Srg,rpf%lambda)
    class is (rpf_KRP3_liq_type)
      rpf_swap => RPFWIPPctor(PETSC_TRUE,3,rpf%Sr,rpf%Srg,rpf%lambda)
    class is (rpf_KRP3_gas_type)
      rpf_swap => RPFWIPPctor(PETSC_FALSE,3,rpf%Sr,rpf%Srg,rpf%lambda)
    class is (rpf_KRP4_liq_type)
      rpf_swap => RPFWIPPctor(PETSC_TRUE,4,rpf%Sr,rpf%Srg,rpf%lambda)
    class is (rpf_KRP4_gas_type)
      rpf_swap => RPFWIPPctor(PETSC_FALSE,4,rpf%Sr,rpf%Srg,rpf%lambda)
    class is (rpf_KRP5_liq_type)
      rpf_swap => RPFWIPPctor(PETSC_TRUE,5,rpf%Sr,rpf%Srg,0d0)
    class is (rpf_KRP5_gas_type)
      rpf_swap => RPFWIPPctor(PETSC_FALSE,5,rpf%Sr,rpf%Srg,0d0)
    class is (rpf_KRP8_liq_type)
      rpf_swap => RPFWIPPctor(PETSC_TRUE,8,rpf%Sr,rpf%Srg,rpf%m)
    class is (rpf_KRP8_gas_type)
      rpf_swap => RPFWIPPctor(PETSC_FALSE,8,rpf%Sr,rpf%Srg,rpf%m)
    class is (rpf_KRP9_liq_type)
      rpf_swap => RPFWIPPctor(PETSC_TRUE,9,rpf%Sr,rpf%Srg,0d0)
    class is (rpf_KRP9_gas_type)
      rpf_swap => RPFWIPPctor(PETSC_FALSE,9,rpf%Sr,rpf%Srg,0d0)
    class is (rpf_KRP11_liq_type)
      rpf_swap => RPFWIPPctor(PETSC_TRUE,11,rpf%Sr,rpf%Srg,rpf%tolc)
    class is (rpf_KRP11_gas_type)
      rpf_swap => RPFWIPPctor(PETSC_FALSE,11,rpf%Sr,rpf%Srg,rpf%tolc)
    class is (rpf_KRP12_liq_type)
      rpf_swap => RPFWIPPctor(PETSC_TRUE,12,rpf%Sr,rpf%Srg,rpf%lambda)
    class is (rpf_KRP12_gas_type)
      rpf_swap => RPFWIPPctor(PETSC_FALSE,12,rpf%Sr,rpf%Srg,rpf%lambda)
    end select
  end if

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

  if (spline) then ! Create cubic approximation for any saturation function
    if (associated(rpf_swap)) then ! Splining a loop-invariant replacement
      rpf_swap2 => RPFSplineCtor(permeability_function, 100)
      deallocate(rpf_swap)
      rpf_swap => rpf_swap2
    else
      rpf_swap => RPFSplineCtor(permeability_function, 100)
    end if
  end if

end function PermeabilityFunctionRead

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
  PetscInt :: i, ii

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
      if (OptionIsIORank(option)) then
        call CharacteristicCurvesTest(cur_characteristic_curves,option)
      endif
      call OptionSetBlocking(option,PETSC_TRUE)
      call OptionCheckNonBlockingError(option)
    endif
    cur_characteristic_curves => cur_characteristic_curves%next
  enddo

  do ii = 1, size(array)
    do i = ii+1, size(array)
      if (StringCompare(array(ii)%ptr%name,array(i)%ptr%name)) then
        option%io_buffer = 'Duplicate characteristic curves named "' // &
          trim(array(ii)%ptr%name) // '".'
        call PrintErrMsg(option)
      endif
    enddo
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
  PetscReal :: sgcr_dummy,sowcr_dummy, sogcr_dummy, swco_dummy
  PetscReal :: gas_res_sat

  select case(option%iflowmode)
    ! leave as case to accommodate OGS classes
    case default
      CharCurvesGetGetResidualSats(1) = &
        characteristic_curves%liq_rel_perm_function%Sr
      ! the Intel compiler on Windows complains about the use of a recursive
      ! subroutine when CharCurvesGetGetResidualSats(option%gas_phase) is
      ! set equal to the residual gas saturation below. Using gas_res_sat
      ! within the select type resolves the issue.
      if (option%nphase > 1 .and. &
          associated(characteristic_curves%gas_rel_perm_function) ) then
        select type(rpf=>characteristic_curves%gas_rel_perm_function)
          class is(rpf_Mualem_VG_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_Mualem_VG_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_Burdine_BC_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_Burdine_BC_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_Mualem_BC_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_Mualem_BC_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_Burdine_VG_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_Burdine_VG_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_TOUGH2_IRP7_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_Mualem_Linear_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_Mualem_Linear_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_Burdine_Linear_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_Burdine_Linear_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_KRP1_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_KRP1_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_KRP2_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_KRP2_gas_type)
            ! KRP2 does not use a Srg, so return 0.d0
            gas_res_sat = 0.d0
          class is(rpf_KRP3_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_KRP3_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_KRP4_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_KRP4_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_KRP5_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_KRP5_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_KRP8_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_KRP8_gas_type)
            ! KRP8 does not use a Srg, so return 0.d0
            gas_res_sat = 0.d0
          class is(rpf_KRP9_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_KRP9_gas_type)
            ! KRP9 does not use a Srg, so return 0.d0
            gas_res_sat = 0.d0
          class is(rpf_KRP11_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_KRP11_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_KRP12_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_KRP12_gas_type)
            gas_res_sat = rpf%Srg
          class is(rpf_mK_liq_type)
            gas_res_sat = rpf%Sr
          class is(rpf_mK_gas_type)
            gas_res_sat = rpf%Srg
          ! class is(rpf_mod_BC_liq_type)
          !   gas_res_sat = rpf%Sr
          class is(rel_perm_func_constant_type)
            gas_res_sat = rpf%Sr
          class is(rel_perm_func_default_type)
            gas_res_sat = rpf%Sr
          class is (rpf_IGHCC2_comp_liq_type)
            gas_res_sat = rpf%Sr
          class is (rpf_IGHCC2_comp_gas_type)
            gas_res_sat = rpf%Srg
          class default
            option%io_buffer = 'Relative permeability class not supported in &
                  &CharCurvesGetGetResidualSats.'
            call PrintErrMsgToDev(option,'')
        end select
        CharCurvesGetGetResidualSats(option%gas_phase) = gas_res_sat
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
    if (option%iflowmode == G_MODE .or. option%iflowmode == WF_MODE .or. &
        option%iflowmode == H_MODE) then
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

    !PO: todo - add cc_owg print out

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

end module Characteristic_Curves_module

