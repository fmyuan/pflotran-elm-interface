module Characteristic_Curves_OWG_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_Common_module
  use Characteristic_Curves_WIPP_module

  implicit none

  private
  
!-----------------------------------------------------------------------------
!-- OWG Saturation Functions -------------------------------------------------
!-----------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  type, public :: sat_func_owg_base_type
    !type(polynomial_type), pointer :: sat_poly
    !type(polynomial_type), pointer :: pres_poly
    PetscReal :: Swco !connate water saturation
    PetscReal :: Swcr !critical (residual) water saturation
    PetscReal :: Soco !connate oil saturation
    PetscReal :: Socr !critical (residual) oil saturation
    PetscReal :: Sgco !connate gas saturation
    PetscReal :: Sgcr !critical (residual) gas saturation
    PetscReal :: pcmax
    !lookup_table_general_type :: lookup_table
    class(sat_func_base_type), pointer :: sat_func_sl
    PetscBool :: analytical_derivative_available
  contains
    procedure, public :: Init => SFOWGBaseInit
    procedure, public :: Verify => SFOWGBaseVerify
    procedure, public :: Test => SFOWGBaseTest
    procedure, public :: SetupPolynomials => SFOWGBaseSetupPolynomials
    procedure, public :: CapillaryPressure => SFOWGBaseCapillaryPressure
    !procedure, public :: Saturation => SFOWGBaseSaturation
  end type sat_func_owg_base_type

   !--------------------------------------------------------------------------
  type, public, extends(sat_func_owg_base_type) :: sat_func_ow_VG_type
    PetscReal :: alpha
    PetscReal :: m
  contains
    procedure, public :: Init => SF_OW_VG_Init
    procedure, public :: Verify => SF_OW_VG_Verify
    procedure, public :: CapillaryPressure => SF_OW_VG_CapillaryPressure
    !procedure, public :: Saturation => SF_OW_VG_Saturation
  end type sat_func_ow_VG_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_owg_base_type) :: sat_func_og_BC_type
    PetscReal :: alpha
    PetscReal :: lambda
  contains
    procedure, public :: Init => SF_OG_BC_Init
    procedure, public :: Verify => SF_OG_BC_Verify
    procedure, public :: SetupPolynomials => SF_OG_BC_SetupPolynomials
    procedure, public :: CapillaryPressure => SF_OG_BC_CapillaryPressure
    !procedure, public :: Saturation => SF_OW_VG_Saturation
  end type sat_func_og_BC_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_owg_base_type) :: sat_func_og_VG_SL_type
    PetscReal :: alpha
    PetscReal :: m
    PetscReal :: Slcr
  contains
    procedure, public :: Init => SF_OG_VG_SL_Init
    procedure, public :: Verify => SF_OG_VG_SL_Verify
    procedure, public :: CapillaryPressure => SF_OG_VG_SL_CapillaryPressure
    !procedure, public :: Saturation => SF_OW_VG_Saturation
  end type sat_func_og_VG_SL_type
  
!-----------------------------------------------------------------------------
!-- Relative Permeability Functions ------------------------------------------
!-----------------------------------------------------------------------------  
  type, public, extends(rel_perm_func_base_type) :: RPF_Mod_BC_type
    PetscReal :: m   !exponential coeff. 
    PetscReal :: Srg 
    PetscReal :: Sro
    PetscReal :: kr_max
  contains
    procedure, public :: Init => RPF_Mod_BC_Init 
    procedure, public :: Verify => RPF_Mod_BC_Verify
    procedure, public :: SetupPolynomials => RPF_Mod_BC_SetupPolynomials
  end type RPF_Mod_BC_type
  !---------------------------------------------------------------------------
  type, public, extends(RPF_Mod_BC_type) :: RPF_Mod_BC_liq_type
  contains
    procedure, public :: RelativePermeability => RPF_Mod_BC_Liq_RelPerm
  end type RPF_Mod_BC_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(RPF_Mod_BC_type) :: RPF_Mod_BC_oil_type
  contains
    procedure, public :: RelativePermeability => RPF_Mod_BC_Oil_RelPerm
  end type RPF_Mod_BC_oil_type  
  
!-----------------------------------------------------------------------------
!-- OWG Relative Permeability Functions --------------------------------------
!-----------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_TOUGH2_Linear_oil_type
    PetscReal :: Sro !
  contains
    procedure, public :: Init => RPF_TOUGH2_Linear_Oil_Init 
    procedure, public :: Verify => RPF_TOUGH2_Linear_Oil_Verify
    procedure, public :: RelativePermeability => RPF_TOUGH2_Linear_Oil_RelPerm
  end type rpf_TOUGH2_Linear_Oil_type  
  !---------------------------------------------------------------------------
  type, public :: rel_perm_func_owg_base_type
    type(polynomial_type), pointer :: poly
    PetscReal :: Swco !connate water saturation
    PetscReal :: Swcr !critical (residual) water saturation
    PetscReal :: Soco !connate oil saturation
    PetscReal :: Socr !critical (residual) oil saturation
    PetscReal :: Sgco !connate gas saturation
    PetscReal :: Sgcr !critical (residual) gas saturation
    PetscReal :: Slcr !critical (residual) saturation of liqui (oil + water)
    PetscReal :: kr_max
    PetscBool :: function_of_liquid_sat
    PetscBool :: So_is_Sh
    !lookup_table_general_type :: lookup_table
    class(rel_perm_func_base_type), pointer :: rel_perm_func_sl
    PetscBool :: analytical_derivative_available
  contains
    procedure, public :: Init => RPFOWGBaseInit
    procedure, public :: Verify => RPFOWGBaseVerify
    procedure, public :: Test => RPFOWGBaseTest
    procedure, public :: SetupPolynomials => RPFOWGBaseSetupPolynomials
    procedure, public :: RelativePermeability => RPFOWGBaseRelPerm
  end type rel_perm_func_owg_base_type

  type, public, extends(rel_perm_func_owg_base_type) :: RPF_OWG_MBC_type
    PetscReal :: m   !exponential coeff.
    !PetscReal :: kr_max
  contains
    procedure, public :: Init => RPF_OWG_MBC_Init
    procedure, public :: Verify => RPF_OWG_MBC_Verify
    procedure, public :: SetupPolynomials => RPF_OWG_MBC_SetupPolynomials
    procedure  :: RPF_OWG_MBC_RelPerm_dkr_dSe
  end type RPF_OWG_MBC_type

  type, public, extends(RPF_OWG_MBC_type) :: RPF_OWG_MBC_wat_type
  contains
    procedure, public :: RelativePermeability => RPF_OWG_MBC_wat_RelPerm
  end type RPF_OWG_MBC_wat_type

  type, public, extends(RPF_OWG_MBC_type) :: RPF_OWG_MBC_oil_type
  contains
    procedure, public :: RelativePermeability => RPF_OWG_MBC_oil_RelPerm
  end type RPF_OWG_MBC_oil_type

  type, public, extends(RPF_OWG_MBC_type) :: RPF_OWG_MBC_gas_type
  contains
    procedure, public :: RelativePermeability => RPF_OWG_MBC_gas_RelPerm
  end type RPF_OWG_MBC_gas_type

  type, public, extends(rel_perm_func_owg_base_type) :: RPF_oil_ecl_type
    class(rel_perm_func_owg_base_type), pointer :: rel_perm_ow
    class(rel_perm_func_owg_base_type), pointer :: rel_perm_og
    !PetscReal :: kr_max
  contains
    procedure, public :: Init => RPF_oil_ecl_Init
    procedure, public :: Verify => RPF_oil_ecl_Verify
    procedure, public :: RelativePermeability => RPF_oil_ecl_RelPerm
  end type RPF_oil_ecl_type

  type, public, extends(rel_perm_func_owg_base_type) :: RPF_OWG_func_sl_type
    !PetscReal :: kr_max
  contains
    procedure, public :: Init => RPF_OWG_func_sl_Init
    procedure, public :: Verify => RPF_OWG_func_sl_Verify
    procedure, public :: RelativePermeability => RPF_OWG_func_sl_RelPerm
  end type RPF_OWG_func_sl_type

  type, public, extends(RPF_OWG_func_sl_type) :: RPF_OWG_func_sl_VG_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_OWG_func_sl_VG_Init
    procedure, public :: Verify => RPF_OWG_func_sl_VG_Verify
  end type RPF_OWG_func_sl_VG_type

  type, public, extends(RPF_OWG_func_sl_VG_type) :: RPF_OWG_Mualem_VG_wat_type
  contains
    procedure, public :: Verify => RPF_OWG_Mualem_VG_wat_Verify
  end type RPF_OWG_Mualem_VG_wat_type

  type, public, extends(RPF_OWG_func_sl_VG_type) :: RPF_OWG_Mualem_VG_gas_type
  contains
    procedure, public :: Verify => RPF_OWG_Mualem_VG_gas_Verify
  end type RPF_OWG_Mualem_VG_gas_type

  type, public, extends(RPF_OWG_func_sl_VG_type) :: &
                                                RPF_OWG_TOUGH2_IRP7_gas_type
  contains
    procedure, public :: Verify => RPF_OWG_TOUGH2_IRP7_gas_Verify
  end type RPF_OWG_TOUGH2_IRP7_gas_type

  type, public, extends(RPF_OWG_func_sl_VG_type) :: RPF_OWG_Burdine_VG_wat_type
  contains
    procedure, public :: Verify => RPF_OWG_Burdine_VG_wat_Verify
  end type RPF_OWG_Burdine_VG_wat_type

  type, public, extends(RPF_OWG_func_sl_VG_type) :: RPF_OWG_Burdine_VG_gas_type
  contains
    procedure, public :: Verify => RPF_OWG_Burdine_VG_gas_Verify
  end type RPF_OWG_Burdine_VG_gas_type

  type, public, extends(RPF_OWG_func_sl_type) :: RPF_OWG_func_sl_BC_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPF_OWG_func_sl_BC_Init
    procedure, public :: Verify => RPF_OWG_func_sl_BC_Verify
  end type RPF_OWG_func_sl_BC_type

  type, public, extends(RPF_OWG_func_sl_BC_type) :: RPF_OWG_Burdine_BC_wat_type
  contains
    procedure, public :: Verify => RPF_OWG_Burdine_BC_wat_Verify
  end type RPF_OWG_Burdine_BC_wat_type

  type, public, extends(RPF_OWG_func_sl_BC_type) :: RPF_OWG_Burdine_BC_gas_type
  contains
    procedure, public :: Verify => RPF_OWG_Burdine_BC_gas_Verify
  end type RPF_OWG_Burdine_BC_gas_type
  
  public :: SaturationFunctionOWGRead, &
            PermeabilityFunctionOWGRead, &
            ! OWG saturation functions
            SF_OW_VG_Create, &
            SF_OG_BC_Create, &
            SF_OG_VG_SL_Create, &
            RPF_TOUGH2_Linear_Oil_Create, &
            RPF_Mod_BC_Liq_Create, &
            RPF_Mod_BC_Oil_Create, &
            RPF_wat_MBC_Create, &
            RPF_oil_MBC_Create, &
            RPF_gas_MBC_Create, &
            RPF_oil_ecl_Create, &
            RPF_OWG_Mualem_VG_wat_Create, &
            RPF_OWG_Mualem_VG_gas_Create, &
            RPF_OWG_TOUGH2_IRP7_gas_Create, &
            RPF_OWG_Burdine_VG_wat_Create, &
            RPF_OWG_Burdine_VG_gas_Create, &
            RPF_OWG_Burdine_BC_wat_Create, &
            RPF_OWG_Burdine_BC_gas_Create, &
            SaturationFunctionOWGDestroy, &
            PermeabilityFunctionOWGDestroy
            
  
contains

! ************************************************************************** !

subroutine SaturationFunctionOWGRead(sat_func_owg,input,option)
  !
  ! Reads in contents of a SATURATION_FUNCTION_OWG block
  !
  ! Author: Paolo Orsini
  ! Date: 11/16/17
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(sat_func_owg_base_type) :: sat_func_owg
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: smooth

  PetscReal :: m_tmp

  input%ierr = 0
  smooth = PETSC_FALSE
  error_string = 'CHARACTERISTIC_CURVES,SATURATION_FUNCTION_OWG,'
  select type(sf_owg => sat_func_owg)
    class is(sat_func_ow_VG_type)
      error_string = trim(error_string) // 'VAN_GENUCHTEN_OW'
    class is(sat_func_og_BC_type)
      error_string = trim(error_string) // 'BROOKS_COREY_OG'
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
      case('MAX_CAPILLARY_PRESSURE')
        call InputReadDouble(input,option,sat_func_owg%pcmax)
        call InputErrorMsg(input,option,'MAX_CAPILLARY_PRESSURE', &
                            error_string)
        if ( associated(sat_func_owg%sat_func_sl) ) then
          sat_func_owg%sat_func_sl%pcmax = sat_func_owg%pcmax
        end if
      case('SMOOTH')
        smooth = PETSC_TRUE
      case default
        found = PETSC_FALSE
    end select

    if (found) cycle

    select type(sf_owg => sat_func_owg)
      class is(sat_func_ow_VG_type)
        select type(sf_sl => sf_owg%sat_func_sl)
          class is(sat_func_VG_type)
            select case(keyword)
              case('M')
                call InputReadDouble(input,option,sf_owg%m)
                call InputErrorMsg(input,option,'M',error_string)
                sf_sl%m = sf_owg%m
              case('ALPHA')
                call InputReadDouble(input,option,sf_owg%alpha)
                call InputErrorMsg(input,option,'ALPHA',error_string)
                sf_sl%alpha = sf_owg%alpha
              case('WATER_RESIDUAL_SATURATION')
                call InputReadDouble(input,option,sf_owg%Swcr)
                call InputErrorMsg(input,option,'WATER_RESIDUAL_SATURATION', &
                                   error_string)
                sf_sl%Sr = sf_owg%Swcr
              case default
                call InputKeywordUnrecognized(keyword, &
                     'Van Genuchten Oil-Water saturation function',option)
            end select
        end select
      class is(sat_func_og_BC_type)
        select type(sf_sl => sf_owg%sat_func_sl)
          class is(sat_func_BC_type)
            select case(keyword)
               case('LAMBDA')
                 call InputReadDouble(input,option,sf_owg%lambda)
                 call InputErrorMsg(input,option,'LAMBDA',error_string)
                 sf_sl%lambda = sf_owg%lambda
               case('ALPHA')
                 call InputReadDouble(input,option,sf_owg%alpha)
                 call InputErrorMsg(input,option,'ALPHA',error_string)
                 sf_sl%alpha = sf_owg%alpha
              case('OIL_RESIDUAL_SATURATION')
                call InputReadDouble(input,option,sf_owg%Socr)
                call InputErrorMsg(input,option,'OIL_RESIDUAL_SATURATION', &
                                   error_string)
                sf_sl%Sr = sf_owg%Socr
              case default
                call InputKeywordUnrecognized(keyword, &
                       'Brooks-Corey Oil-Gas saturation function',option)
            end select
          end select
      class is(sat_func_og_VG_SL_type)
        select type(sf_sl => sf_owg%sat_func_sl)
          class is(sat_func_VG_type)
            select case(keyword)
              case('M')
                call InputReadDouble(input,option,sf_owg%m)
                call InputErrorMsg(input,option,'M',error_string)
                sf_sl%m = sf_owg%m
              case('ALPHA')
                call InputReadDouble(input,option,sf_owg%alpha)
                call InputErrorMsg(input,option,'ALPHA',error_string)
                sf_sl%alpha = sf_owg%alpha
              case('LIQUID_RESIDUAL_SATURATION')
                call InputReadDouble(input,option,sf_owg%Slcr)
                call InputErrorMsg(input,option,'LIQUID_RESIDUAL_SATURATION', &
                                   error_string)
                sf_sl%Sr = sf_owg%Slcr
              case default
                call InputKeywordUnrecognized(keyword, &
                       'Van Genuchten Oil-Gas-SL saturation function',option)
            end select
        end select
      !add here table case
    end select
  !add reading instructions for other OWG saturation functions (tables etc)
  end do

  !pass on parameters to sl function if needed

  if ( smooth .and. associated(sat_func_owg%sat_func_sl) ) then
    call sat_func_owg%sat_func_sl%SetupPolynomials(option,error_string)
  end if

end subroutine SaturationFunctionOWGRead

! ************************************************************************** !

recursive subroutine PermeabilityFunctionOWGRead(permeability_function, &
                                                 phase_keyword,input,option)
  ! 
  ! Reads in contents of a PERMEABILITY_FUNCTION_OWG block
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  
  implicit none
  
  class(rel_perm_func_owg_base_type) :: permeability_function
  character(len=MAXWORDLENGTH) :: phase_keyword
  type(input_type), pointer :: input
  type(option_type) :: option
    
  class(rel_perm_func_owg_base_type), pointer :: perm_func_ptr => null()
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: perm_interface, perm_func_ch_type
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: smooth
  character(len=MAXWORDLENGTH) :: oilstring

  oilstring='OIL'
  
  input%ierr = 0
  smooth = PETSC_FALSE
  error_string = 'CHARACTERISTIC_CURVES_OWG,PERMEABILITY_FUNCTION,'
  select type(rpf => permeability_function)
    class is(RPF_OWG_MBC_wat_type)
      error_string = trim(error_string) // 'MOD_BROOKS_COREY_WAT'
    class is(RPF_OWG_MBC_oil_type)
      if (rpf%So_is_Sh) then
        error_string = trim(error_string) // 'MOD_BROOKS_COREY_HYDROCARBON'
      else
        error_string = trim(error_string) // 'MOD_BROOKS_COREY_OIL'
      end if
    class is(RPF_OWG_MBC_gas_type)
      error_string = trim(error_string) // 'MOD_BROOKS_COREY_GAS'
    class is(RPF_oil_ecl_type)
      error_string = trim(error_string) // 'ECLIPSE_OIL'
    class is(RPF_OWG_Mualem_VG_wat_type)
      error_string = trim(error_string) // 'MUALEM_VG_WAT'
    class is(RPF_OWG_Mualem_VG_gas_type)
      error_string = trim(error_string) // 'MUALEM_VG_GAS_SL'
    class is(RPF_OWG_TOUGH2_IRP7_gas_type)
        error_string = trim(error_string) // 'TOUGH2_IRP7_GAS_SL'
    class is(RPF_OWG_Burdine_VG_wat_type)
      error_string = trim(error_string) // 'BURDINE_VG_WAT'
    class is(RPF_OWG_Burdine_VG_gas_type)
        error_string = trim(error_string) // 'BURDINE_VG_GAS_SL'
    class is(RPF_OWG_Burdine_BC_wat_type)
      error_string = trim(error_string) // 'BURDINE_BC_WAT'
    class is(RPF_OWG_Burdine_BC_gas_type)
        error_string = trim(error_string) // 'BURDINE_BC_GAS_SL'
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
      case('WATER_CONNATE_SATURATION')
        call InputReadDouble(input,option,permeability_function%Swco)
        call InputErrorMsg(input,option,'WATER_CONNATE_SATURATION', &
                           error_string)
      case('WATER_RESIDUAL_SATURATION')
        call InputReadDouble(input,option,permeability_function%Swcr)
        call InputErrorMsg(input,option,'WATER_RESIDUAL_SATURATION', &
                           error_string)
        if(associated(permeability_function%rel_perm_func_sl) ) then
          permeability_function%rel_perm_func_sl%Sr = &
                                                    permeability_function%Swcr
        end if
      case('OIL_CONNATE_SATURATION')
        call InputReadDouble(input,option,permeability_function%Soco)
        call InputErrorMsg(input,option,'OIL_CONNATE_SATURATION', &
                          error_string)
      case('OIL_RESIDUAL_SATURATION')
        call InputReadDouble(input,option,permeability_function%Socr)
        call InputErrorMsg(input,option,'OIL_RESIDUAL_SATURATION', &
                           error_string)
      case('HYDROCARBON_RESIDUAL_SATURATION')
        call InputReadDouble(input,option,permeability_function%Socr)
        call InputErrorMsg(input,option,'HYDROCARBON_RESIDUAL_SATURATION', &
                          error_string)
        permeability_function%Sgcr = 0.0
      case('GAS_CONNATE_SATURATION')
        call InputReadDouble(input,option,permeability_function%Sgco)
        call InputErrorMsg(input,option,'GAS_CONNATE_SATURATION', &
                          error_string)
      case('GAS_RESIDUAL_SATURATION')
        call InputReadDouble(input,option,permeability_function%Sgcr)
        call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                           error_string)
        if(associated(permeability_function%rel_perm_func_sl) ) then
          select type(rpf_sl => permeability_function%rel_perm_func_sl)
            class is(RPF_Mualem_VG_gas_type)
              rpf_sl%Srg = permeability_function%Sgcr
            class is(RPF_TOUGH2_IRP7_gas_type)
              rpf_sl%Srg = permeability_function%Sgcr
            class is(RPF_Burdine_VG_gas_type)
              rpf_sl%Srg = permeability_function%Sgcr
            class is(RPF_Burdine_BC_gas_type)
              rpf_sl%Srg = permeability_function%Sgcr
          end select
        end if
      case('LIQUID_RESIDUAL_SATURATION')
        call InputReadDouble(input,option,permeability_function%Slcr)
        call InputErrorMsg(input,option,'LIQUID_RESIDUAL_SATURATION', &
                           error_string)
        if(associated(permeability_function%rel_perm_func_sl) ) then
          permeability_function%rel_perm_func_sl%Sr = &
                                                    permeability_function%Slcr
        end if
      case('MAX_RELATIVE_PERMEABILITY','MAX_REL_PERM')
        call InputReadDouble(input,option,permeability_function%kr_max)
        call InputErrorMsg(input,option,'MAX_RELATIVE_PERMEABILITY', &
                           error_string)
      case('SMOOTH')
        smooth = PETSC_TRUE
      case default
        found = PETSC_FALSE
    end select
    if (found) cycle

    select type(rpf => permeability_function)
    !------------------------------------------
      class is(RPF_OWG_MBC_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'M',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Modified Brooks-Corey relative permeability function', &
              option)
        end select
      class is(RPF_oil_ecl_type)
        perm_interface = 'NONE'
        select case(keyword)
          case('PERMEABILITY_FUNCTION_OW')
            perm_interface = 'OIL_WATER'
            call InputReadWord(input,option,perm_func_ch_type,PETSC_TRUE)
            call InputErrorMsg(input,option,'perm_func_ch_type',error_string)
          case('PERMEABILITY_FUNCTION_OG')
            perm_interface = 'OIL_GAS'
            call InputReadWord(input,option,perm_func_ch_type,PETSC_TRUE)
            call InputErrorMsg(input,option,'perm_func_ch_type',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                      'ECLIPSE_OIL relative permeability function',option)
          end select
          call StringToUpper(perm_func_ch_type)
          perm_func_ch_type = trim(perm_func_ch_type)
        select case(perm_func_ch_type)
          case('MOD_BROOKS_COREY_OIL')
            perm_func_ptr => RPF_oil_MBC_Create()
          case default
            call InputKeywordUnrecognized(perm_func_ch_type, &
                                          'PERMEABILITY_FUNCTION_OW/OG',option)
        end select
        call PermeabilityFunctionOWGRead(perm_func_ptr,oilstring,input,option)
        select case(perm_interface)
          case('OIL_WATER')
            rpf%rel_perm_ow => perm_func_ptr
          case('OIL_GAS')
            rpf%rel_perm_og => perm_func_ptr
        end select
        nullify(perm_func_ptr)
      class is(RPF_OWG_func_sl_VG_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'M', &
                               error_string)
            select type(rpf_sl => rpf%rel_perm_func_sl)
              class is(RPF_Mualem_VG_liq_type)
                rpf_sl%m = rpf%m
              class is(RPF_Mualem_VG_gas_type)
                rpf_sl%m = rpf%m
              class is(RPF_TOUGH2_IRP7_gas_type)
                rpf_sl%m = rpf%m
              class is(RPF_Burdine_VG_liq_type)
                rpf_sl%m = rpf%m
              class is(RPF_Burdine_VG_gas_type)
                rpf_sl%m = rpf%m
            end select
        end select
      class is(RPF_OWG_func_sl_BC_type)
        select case(keyword)
          case('LAMBDA')
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA', &
                               error_string)
             select type(rpf_sl => rpf%rel_perm_func_sl)
               class is(RPF_Burdine_BC_liq_type)
                 rpf_sl%lambda = rpf%lambda
               class is(RPF_Burdine_BC_gas_type)
                 rpf_sl%lambda = rpf%lambda
             end select
          end select
    end select
  end do

  !When a smoother for the sl function exists the OWG rel per smoother
  ! is not defined
  if ( associated(permeability_function%rel_perm_func_sl) ) then
    if (smooth) then
      call permeability_function%rel_perm_func_sl% &
                                 SetupPolynomials(option,error_string)
    end if
  else
    if (smooth) then
      call permeability_function%SetupPolynomials(option,error_string)
    endif
  end if

end subroutine PermeabilityFunctionOWGRead

! ************************************************************************** !
! ************************************************************************** !

function RPF_TOUGH2_Linear_Oil_Create()

  ! Creates the TOUGH2 Linear oil relative permeability function object
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/19/2015

  class(rpf_TOUGH2_Linear_oil_type), pointer :: RPF_TOUGH2_Linear_Oil_Create

  allocate(RPF_TOUGH2_Linear_Oil_Create)
  call RPF_TOUGH2_Linear_Oil_Create%Init()

end function RPF_TOUGH2_Linear_Oil_Create

! ************************************************************************** !

subroutine RPF_TOUGH2_Linear_Oil_Init(this)

  ! Initializes the TOUGH2 Linear Oil relative permeability function 
  ! object
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/19/2015

  implicit none
  
  class(rpf_TOUGH2_Linear_oil_type) :: this

  call RPFBaseInit(this)
  this%Sro = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_TOUGH2_Linear_Oil_Init

! ************************************************************************** !

subroutine RPF_TOUGH2_Linear_Oil_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_TOUGH2_Linear_oil_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,TOUGH2_LINEAR_OIL'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Sro)) then
    option%io_buffer = UninitializedMessage('OIL_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  
end subroutine RPF_TOUGH2_Linear_Oil_Verify

! ************************************************************************** !

subroutine RPF_TOUGH2_Linear_Oil_RelPerm(this,liquid_saturation, &
                                         relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/19/2015


  use Option_module
  use Utility_module, only : InitToNan
  
  implicit none

  class(rpf_TOUGH2_Linear_oil_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: So
  PetscReal :: Seo
  
  ! initialize to derivative to NaN so that not mistakenly used.
  dkr_sat = InitToNan()

  So = 1.d0 - liquid_saturation

  Seo = (So - this%Sro) / (1.d0 - this%Sro)

  !!! DS added
  !Seo = (So - this%Sro) / (1.d0 - this%Sro)
  !dSeo_so = 1.d0 / (1.d0 - this%Sro)
  !dkr_Seo = 1.d0
  !! return dkr_sat is derivative of relperm wrt liquid saturation to 
  !! be consistent 
  !dkr_sat = -1.d0*dSeo_so*dkr_Seo !! -1.d0 is dso/dsl
  !! just hard code it:
  dkr_sat = -1.d0 / (1.d0 - this%Sro)

  if (Seo >= 1.d0) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
    return
  else if (Seo <=  0.d0) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
    return
  endif

  relative_permeability = Seo

end subroutine RPF_TOUGH2_Linear_Oil_RelPerm

! ************************************************************************** !
! ************************************************************************** !

!Beginning RPF Modified Brooks-Corey for liq and oil phase (RPF_Mod_BC_Oil)

!  procedure, public :: Init => RPF_Mod_BC_Oil_Init 
!  procedure, public :: Verify => RPF_Mod_BC_Oil_Verify
!  procedure, public :: SetupPolynomials => RPF_Mod_BC_SetupPolynomials
!  procedure, public :: RelativePermeability => RPF_Mod_BC_Oil_RelPerm

function RPF_Mod_BC_Liq_Create()

  ! Creates the Modified BC Oil relative permeability function object
  ! Author: Paolo Orsini (OGS)
  ! Date: 02/20/2016

  class(rpf_mod_BC_liq_type), pointer :: RPF_Mod_BC_Liq_Create

  allocate(RPF_Mod_BC_Liq_Create)
  call RPF_Mod_BC_Liq_Create%Init()

end function RPF_Mod_BC_Liq_Create

! ************************************************************************** !

function RPF_Mod_BC_Oil_Create()

  ! Creates the Modified BC Oil relative permeability function object
  ! Author: Paolo Orsini (OGS)
  ! Date: 02/20/2016

  class(rpf_mod_BC_oil_type), pointer :: RPF_Mod_BC_Oil_Create

  allocate(RPF_Mod_BC_Oil_Create)
  call RPF_Mod_BC_Oil_Create%Init()

end function RPF_Mod_BC_Oil_Create

! ************************************************************************** !

!subroutine RPF_Mod_BC_Oil_Init(this)
subroutine RPF_Mod_BC_Init(this)

  ! Initializes the Modified BC Oil relative permeability function object 
  ! object
  ! Author: Paolo Orsini (OGS)
  ! Date: 02/20/2016

  implicit none
  
  !class(rpf_mod_BC_oil_type) :: this
  class(rpf_mod_BC_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  this%Sro = UNINITIALIZED_DOUBLE
  this%kr_max = 1.0d0
  
  this%analytical_derivative_available = PETSC_TRUE
   
!end subroutine RPF_Mod_BC_Oil_Init
end subroutine RPF_Mod_BC_Init

! ************************************************************************** !

!subroutine RPF_Mod_BC_Oil_Verify(this,name,option)
subroutine RPF_Mod_BC_Verify(this,name,option)

  use Option_module

  implicit none
  
  !class(rpf_mod_BC_oil_type) :: this
  class(rpf_mod_BC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    select type(rpf => this)
      class is(rpf_mod_BC_liq_type) 
        string = trim(name) // 'PERMEABILITY_FUNCTION,MOD_BC_LIQ'
      class is(rpf_mod_BC_oil_type)
        string = trim(name) // 'PERMEABILITY_FUNCTION,MOD_BC_OIL'
    end select
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Sro)) then
    option%io_buffer = UninitializedMessage('OIL_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  

  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('POWER EXPONENT',string)
    call printErrMsg(option)
  endif
  
!end subroutine RPF_Mod_BC_Oil_Verify
end subroutine RPF_Mod_BC_Verify

! ************************************************************************** !

subroutine RPF_Mod_BC_SetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing Modified BC permeability function

  use Option_module
  use Utility_module
  
  implicit none
  
  class(rpf_mod_BC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  PetscReal :: b(4)

  PetscReal :: Se_ph_low

  this%poly => PolynomialCreate()
  ! fill matix with values
  this%poly%low = 0.99d0  ! just below saturated
  !this%poly%low = 0.95d0  ! just below saturated 
  this%poly%high = 1.d0   ! saturated
  Se_ph_low = this%poly%low
  !select type(rpf => this)
  !  class is(rpf_mod_BC_liq_type) 
  !    Se_ph_low = ( this%poly%low - this%Sr ) / &
  !                (1.0 - this%Sro - this%Sr - this%Srg)
  !  class is(rpf_mod_BC_oil_type)
  !    Se_ph_low = ( this%poly%low - this%Sro ) / &
  !                (1.0 - this%Sro - this%Sr - this%Srg) 
  !end select 

  b(1) = this%kr_max
  b(2) = this%kr_max * (Se_ph_low ** this%m)
  b(3) = 0.d0
  b(4) = this%m * this%kr_max * Se_ph_low ** (this%m - 1.0 )
  
  call CubicPolynomialSetup(this%poly%high,this%poly%low,b)
  
  this%poly%coefficients(1:4) = b(1:4)
  
end subroutine RPF_Mod_BC_SetupPolynomials

! ************************************************************************** !

subroutine RPF_Mod_BC_Liq_RelPerm(this,liquid_saturation, &
                                  relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 02/21/2016

  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Mod_BC_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dSe_Sl, dkr_Se
  
  ! initialize to derivative to NaN so that not mistakenly used.
  dkr_sat = InitToNan()

  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sro - this%Sr - this%Srg )
  dSe_Sl = 1.d0 / (1.d0 - this%Sro - this%Sr - this%Srg )

  dkr_Se = this%m * this%kr_max * (Se ** (this%m - 1))
  dkr_sat = dSe_Sl * dkr_Se

  if (Se >= 1.d0) then
    relative_permeability = this%kr_max
    dkr_sat = 0.d0
    return
  else if (Se <=  0.d0) then
    dkr_sat = 0.d0
    relative_permeability = 0.d0
    return
  endif

  if (associated(this%poly)) then
    if (Se > this%poly%low) then
      call CubicPolynomialEvaluate(this%poly%coefficients, &
                                   Se,relative_permeability,dkr_Se)
      dkr_sat = dSe_Sl * dkr_Se
      return
    endif
  endif

  relative_permeability = this%kr_max * (Se ** this%m)

end subroutine RPF_Mod_BC_Liq_RelPerm

! ************************************************************************** !

subroutine RPF_Mod_BC_Oil_RelPerm(this,liquid_saturation, &
                                  relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 02/20/2016

  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Mod_BC_oil_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: So
  PetscReal :: Seo
  PetscReal :: dSe_So, dkr_Se
  
  ! initialize to derivative to NaN so that not mistakenly used.
  dkr_sat = InitToNan()

  So = 1.d0 - liquid_saturation

  Seo = (So - this%Sro) / (1.d0 - this%Sro - this%Sr - this%Srg ) 
  dSe_So = 1.d0 / (1.d0 - this%Sro - this%Sr - this%Srg )

  dkr_Se = this%m * this%kr_max * (Seo ** (this%m - 1))

  dkr_sat = -1.d0 * dSe_So * dkr_Se ! -1.d0 factor makes derivative w.r.t. Sl

  if (Seo >= 1.d0) then
    relative_permeability = this%kr_max
    dkr_sat = 0.d0
    return
  else if (Seo <=  0.d0) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
    return
  endif

  if (associated(this%poly)) then
    if (Seo > this%poly%low) then
      call CubicPolynomialEvaluate(this%poly%coefficients, &
                                   Seo,relative_permeability,dkr_Se)
      dkr_sat = -1.d0 * dSe_So * dkr_Se ! -1.d0 factor makes derivative w.r.t. Sl
      return
    endif
  endif

  relative_permeability = this%kr_max * (Seo ** this%m)

end subroutine RPF_Mod_BC_Oil_RelPerm

! ************************************************************************** !
! *********** OWG Saturaton functions  ************************************* !
! ************************************************************************** !

subroutine SFOWGBaseInit(this)

  implicit none
  
  class(sat_func_owg_base_type) :: this
  
  !nullify(this%sat_poly)
  !nullify(this%pres_poly)
  this%Swco = 0.0d0
  this%Swcr = UNINITIALIZED_DOUBLE
  this%Soco = 0.0
  this%Socr = UNINITIALIZED_DOUBLE
  this%Sgco = 0.0d0
  this%Sgcr = UNINITIALIZED_DOUBLE
  this%pcmax = DEFAULT_PCMAX
  this%analytical_derivative_available = PETSC_FALSE
  !nullify(lookup_table)
  nullify(this%sat_func_sl)

end subroutine SFOWGBaseInit

! ************************************************************************** !
subroutine SFOWGBaseVerify(this,name,option)

  use Option_module

  implicit none
  
  class(sat_func_owg_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  if ((.not.this%analytical_derivative_available) .and. &
      (.not.option%flow%numerical_derivatives)) then
    option%io_buffer = 'Analytical derivatives are not available for the &
      &capillary pressure - saturation function chosen: ' // &
      trim(name)
    call printErrMsg(option)
  endif

end subroutine SFOWGBaseVerify

! ************************************************************************** !

subroutine SFOWGBaseTest(this,cc_name,option)

  use Option_module
  use Material_Aux_class

  implicit none

  class(sat_func_owg_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string, sat_name
  PetscInt, parameter :: num_values = 101
  PetscReal :: capillary_pressure(num_values)
  PetscReal :: wat_saturation(num_values)
  PetscReal :: oil_saturation(num_values)
  PetscReal :: gas_saturation(num_values)
  PetscReal :: saturation(num_values)
  PetscReal :: dpc_dsato(num_values),dpc_dsatg(num_values)
  PetscInt :: i


 ! calculate capillary pressure as a function of saturation
  do i = 1, num_values
    wat_saturation(i) = dble(i-1)*0.01d0
    if (wat_saturation(i) < 1.d-7) then
      wat_saturation(i) = 1.d-7
    else if (wat_saturation(i) > (1.d0-1.d-7)) then
      wat_saturation(i) = 1.d0-1.d-7
    endif
  end do

  write(string,*) cc_name
  select type(sf => this)
    class is(sat_func_ow_VG_type)
      oil_saturation = 1.0 - wat_saturation
      gas_saturation = 0.0
      saturation = wat_saturation
      string = trim(cc_name) // '_pc_OW_wat_sat.dat'
      sat_name = 'wat_sat'
    class is(sat_func_og_BC_type)
      oil_saturation = wat_saturation
      gas_saturation = 0.0
      wat_saturation = 0.0
      saturation = oil_saturation
      string = trim(cc_name) // '_pc_OG_oil_sat.dat'
      sat_name = 'oil_sat'
    class is(sat_func_og_VG_SL_type)
      gas_saturation = 1.d0 - wat_saturation
      oil_saturation = 0.0
      saturation = wat_saturation
      string = trim(cc_name) // '_pc_OG_liq_sat.dat'
      sat_name = 'liq_sat'
  end select
  !sat_name = trim(sat_name)

  do i = 1, num_values
   call this%CapillaryPressure(oil_saturation(i), gas_saturation(i), &
                               capillary_pressure(i),dpc_dsato(i), &
                               dpc_dsatg(i),option)
    ! calculate numerical derivatives?
  enddo

  open(unit=86,file=string)
  write(86,*) '"',trim(sat_name), '"', ', "capillary pressure", "dpc/dsat0", &
              &dpc/dsatg"'
  do i = 1, num_values
    write(86,'(4es14.6)') saturation(i), capillary_pressure(i), &
                          dpc_dsato(i), dpc_dsatg(i)
  enddo
  close(86)

end subroutine SFOWGBaseTest

! ************************************************************************** !

subroutine SFOWGBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing saturation functions

  use Option_module

  implicit none

  class(sat_func_owg_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  option%io_buffer = 'SF OWG Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)

end subroutine SFOWGBaseSetupPolynomials

! ************************************************************************** !

subroutine SFOWGBaseCapillaryPressure(this,oil_saturation, gas_saturation, &
                                      capillary_pressure,dpc_dsato,dpc_dsatg, &
                                      option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_owg_base_type) :: this
  PetscReal, intent(in) :: oil_saturation
  PetscReal, intent(in) :: gas_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsato
  PetscReal, intent(out) :: dpc_dsatg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  option%io_buffer = 'SFOWGBaseCapillaryPressure must be extended.'
  call printErrMsg(option)

end subroutine SFOWGBaseCapillaryPressure

! ************************************************************************** !

function SF_OW_VG_Create()

  implicit none

  class(sat_func_ow_VG_type), pointer :: SF_OW_VG_Create

  allocate(SF_OW_VG_Create)

  call SFOWGBaseInit(SF_OW_VG_Create)

  SF_OW_VG_Create%sat_func_sl => SF_VG_Create()

  call SF_OW_VG_Create%Init()


end function SF_OW_VG_Create

! ************************************************************************** !

subroutine SF_OW_VG_Init(this)

  implicit none

  class(sat_func_ow_VG_type) :: this

  this%alpha = UNINITIALIZED_DOUBLE
  this%m = UNINITIALIZED_DOUBLE

  this%analytical_derivative_available = &
      this%sat_func_sl%analytical_derivative_available

end subroutine SF_OW_VG_Init

! ************************************************************************** !

subroutine SF_OW_VG_Verify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_ow_VG_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION_OWG,SF_OW_VG'
  endif

  call SFOWGBaseVerify(this,string,option)

  if (Uninitialized(this%Swcr)) then
    option%io_buffer = UninitializedMessage('WATER_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA = 1/Pcc',string)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call printErrMsg(option)
  endif

  if ( .not.associated(this%sat_func_sl) ) then
    option%io_buffer = 'Not analytical function associated with the capillary &
      & pressure between water and oil - saturation function chosen: ' // &
      trim(string)
    call printErrMsg(option)
  end if

  call this%sat_func_sl%verify(name,option)

end subroutine SF_OW_VG_Verify

! ************************************************************************** !

subroutine SF_OW_VG_CapillaryPressure(this,oil_saturation, gas_saturation, &
                                      capillary_pressure,dpc_dsato,dpc_dsatg, &
                                      option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_ow_VG_type) :: this
  PetscReal, intent(in) :: oil_saturation
  PetscReal, intent(in) :: gas_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsato
  PetscReal, intent(out) :: dpc_dsatg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: water_saturation, dpc_dsatw

  water_saturation = 1.0 - oil_saturation - gas_saturation

  call this%sat_func_sl%CapillaryPressure(water_saturation,capillary_pressure, &
                                          dpc_dsatw,option)

  dpc_dsato = - dpc_dsatw

  dpc_dsatg = - dpc_dsatw

end subroutine SF_OW_VG_CapillaryPressure

! ************************************************************************** !

function SF_OG_BC_Create()

  implicit none

  class(sat_func_og_BC_type), pointer :: SF_OG_BC_Create

  allocate(SF_OG_BC_Create)

  call SFOWGBaseInit(SF_OG_BC_Create)

  SF_OG_BC_Create%sat_func_sl => SF_BC_Create()

  call SF_OG_BC_Create%Init()

end function SF_OG_BC_Create

! ************************************************************************** !

subroutine SF_OG_BC_Init(this)

  implicit none

  class(sat_func_og_BC_type) :: this

  this%lambda =  UNINITIALIZED_DOUBLE
  this%alpha =  UNINITIALIZED_DOUBLE

  this%analytical_derivative_available = &
      this%sat_func_sl%analytical_derivative_available

end subroutine SF_OG_BC_Init

! ************************************************************************** !

subroutine SF_OG_BC_Verify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_og_BC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION_OWG,SF_OG_BC'
  endif

  if ( .not.associated(this%sat_func_sl) ) then
    option%io_buffer = 'Not analytical function associated with the capillary &
      & pressure between oil and gas - saturation function chosen: ' // &
      trim(string)
    call printErrMsg(option)
  end if

  call SFOWGBaseVerify(this,string,option)

  if (Uninitialized(this%Socr)) then
    option%io_buffer = UninitializedMessage('OIL_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('alpha = 1/Pcc',string)
    call printErrMsg(option)
  endif

  call this%sat_func_sl%verify(name,option)

end subroutine SF_OG_BC_Verify

! ************************************************************************** !

subroutine SF_OG_BC_SetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing Brooks-Corey saturation function

  use Option_module
  use Utility_module

  implicit none

  class(sat_func_og_BC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  call this%sat_func_sl%SetupPolynomials(option,error_string)

end subroutine SF_OG_BC_SetupPolynomials

! ************************************************************************** !

subroutine SF_OG_BC_CapillaryPressure(this,oil_saturation, gas_saturation, &
                                      capillary_pressure,dpc_dsato,dpc_dsatg, &
                                      option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_og_BC_type) :: this
  PetscReal, intent(in) :: oil_saturation
  PetscReal, intent(in) :: gas_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsato
  PetscReal, intent(out) :: dpc_dsatg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call this%sat_func_sl%CapillaryPressure(oil_saturation,capillary_pressure, &
                                          dpc_dsato,option)

  dpc_dsatg = 0.0d0

end subroutine SF_OG_BC_CapillaryPressure

! ************************************************************************** !

function SF_OG_VG_SL_Create()

  ! Creates the van Genutchten capillary pressure function object for use
  ! between the oil and gas phase treating OIl and Water joinly as liquid

  implicit none

  class(sat_func_og_VG_SL_type), pointer :: SF_OG_VG_SL_Create

  allocate(SF_OG_VG_SL_Create)

  call SFOWGBaseInit(SF_OG_VG_SL_Create)

  SF_OG_VG_SL_Create%sat_func_sl => SF_VG_Create()

  call SF_OG_VG_SL_Create%Init()


end function SF_OG_VG_SL_Create

! ************************************************************************** !

subroutine SF_OG_VG_SL_Init(this)

  implicit none

  class(sat_func_og_VG_SL_type) :: this

  this%alpha =  UNINITIALIZED_DOUBLE
  this%m =  UNINITIALIZED_DOUBLE
  this%Slcr = UNINITIALIZED_DOUBLE

  this%analytical_derivative_available = &
      this%sat_func_sl%analytical_derivative_available

end subroutine SF_OG_VG_SL_Init

! ************************************************************************** !

subroutine SF_OG_VG_SL_Verify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_og_VG_SL_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION_OWG,SF_OG_VG_SL'
  endif

  if ( .not.associated(this%sat_func_sl) ) then
    option%io_buffer = 'Not analytical function associated with the capillary &
      & pressure between oil and gas - saturation function chosen: ' // &
      trim(string)
    call printErrMsg(option)
  end if

  call SFOWGBaseVerify(this,string,option)

  if (Uninitialized(this%Slcr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('alpha = 1/Pcc',string)
    call printErrMsg(option)
  endif

  call this%sat_func_sl%verify(name,option)

end subroutine SF_OG_VG_SL_Verify

! ************************************************************************** !

subroutine SF_OG_VG_SL_CapillaryPressure(this,oil_saturation, gas_saturation, &
                                       capillary_pressure,dpc_dsato,dpc_dsatg, &
                                       option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_og_VG_SL_type) :: this
  PetscReal, intent(in) :: oil_saturation
  PetscReal, intent(in) :: gas_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsato
  PetscReal, intent(out) :: dpc_dsatg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: liq_saturation, dpc_dsatl

  liq_saturation = 1.0 - gas_saturation

  call this%sat_func_sl%CapillaryPressure(liq_saturation,capillary_pressure, &
                                          dpc_dsatl,option)

  dpc_dsato = dpc_dsatl
  dpc_dsatg = - dpc_dsatl


end subroutine SF_OG_VG_SL_CapillaryPressure

! ************************************************************************** !

! ************************************************************************** !
! *********** END OWG Saturation functions    ****************************** !
! ************************************************************************** !

! ************************************************************************** !
! *********** OWG Relative Permeability functions  ************************* !
! ************************************************************************** !

! ************************************************************************** !

subroutine RPFOWGBaseInit(this)

  implicit none

  class(rel_perm_func_owg_base_type) :: this

  nullify(this%poly)

  this%Swco = 0.0d0
  this%Swcr = UNINITIALIZED_DOUBLE
  this%Soco = 0.0d0
  this%Socr = UNINITIALIZED_DOUBLE
  this%Sgco = 0.0d0
  this%Sgcr = UNINITIALIZED_DOUBLE
  this%Slcr = UNINITIALIZED_DOUBLE
  this%kr_max = 1.0d0

  this%analytical_derivative_available = PETSC_FALSE
  this%function_of_liquid_sat = PETSC_FALSE
  this%So_is_Sh = PETSC_FALSE

  nullify(this%rel_perm_func_sl)

end subroutine RPFOWGBaseInit

! ************************************************************************** !

subroutine RPFOWGBaseVerify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_func_owg_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  ! by default kr_max = 1.0 - if entered the value must be bound beween 0 and 1
  if ( this%kr_max < 0.0 .or. this%kr_max > 1.0 ) then
    option%io_buffer = adjustl(trim(name)) // ' MAX_REL_PERM entered &
                                   &not valid, must be 0 <= Kr_max <= 1'
    call printErrMsg(option)
  end if

  if ((.not.this%analytical_derivative_available) .and. &
      (.not.option%flow%numerical_derivatives)) then
    option%io_buffer = 'Analytical derivatives are not available for the &
      &relative permeability function chosen: ' // trim(name)
    call printErrMsg(option)
  endif

end subroutine RPFOWGBaseVerify

! ************************************************************************** !

subroutine RPFOWGBaseTest(this,cc_name,phase,option)

  use Option_module

  implicit none

  class(rel_perm_func_owg_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  character(len=MAXWORDLENGTH) :: phase
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, parameter :: num_values = 101
  PetscReal :: saturation(num_values)
  PetscReal :: wat_saturation(num_values)
  PetscReal :: oil_saturation(num_values)
  PetscReal :: gas_saturation(num_values)
  PetscReal :: kr(num_values), dkr_sato(num_values), dkr_satg(num_values)
  PetscReal :: krow(num_values), dkrow_sato(num_values),dkrow_satg(num_values)
  PetscReal :: krog(num_values), dkrog_sato(num_values),dkrog_satg(num_values)

  saturation = 0.0
  wat_saturation = 0.0
  gas_saturation = 0.0
  oil_saturation = 0.0
  kr = 0.0
  dkr_sato = 0.0
  dkr_satg = 0.0
  krow = 0.0
  dkrow_sato = 0.0
  dkrow_satg = 0.0
  krog = 0.0
  dkrog_sato = 0.0
  dkrog_satg = 0.0

  do i = 1, num_values
    saturation(i) = dble(i-1)*0.01d0
  enddo

  select type(this)
    class is(RPF_OWG_MBC_wat_type)
      oil_saturation = 1.0 - saturation
      gas_saturation = 0.0
    class is(RPF_OWG_MBC_oil_type)
      oil_saturation = saturation
      gas_saturation = 0.0
    class is(RPF_oil_ecl_type)
      oil_saturation = saturation
      gas_saturation = 0.0
    class is(RPF_OWG_MBC_gas_type)
      gas_saturation = saturation
      oil_saturation = 0.0
    class is(RPF_OWG_func_sl_type)
      if (phase == 'water') then
        gas_saturation = 1.0 - saturation
        oil_saturation = 0.0
      else if (phase == 'gas') then
        gas_saturation = saturation
        oil_saturation = 0.0
      end if
  end select

  write(string,*) cc_name
  select type(this)
    class is(RPF_oil_ecl_type)
      do i = 1, num_values
        call this%rel_perm_ow%RelativePermeability(oil_saturation(i), &
                                   gas_saturation(i),krow(i),dkrow_sato(i), &
                                   dkrow_satg(i),option)
        call this%rel_perm_og%RelativePermeability(oil_saturation(i), &
                                   gas_saturation(i),krog(i),dkrog_sato(i), &
                                   dkrog_satg(i),option)
      end do
      !print ow rel perm
      string = trim(cc_name) // '_' // trim(phase) // '_ow_rel_perm.dat'
      open(unit=86,file=string)
      write(86,*) '"oil_saturation", "' // trim(phase) // '_ow_rel_perm", &
                  &"dkrow/dsato", " dkrow/dsatg"'
      do i = 1, size(saturation)
        write(86,'(4es14.6)') oil_saturation(i), krow(i), dkrow_sato(i), &
                              dkrow_satg(i)
      enddo
      close(86)
      !print og rel perm
      string = trim(cc_name) // '_' // trim(phase) // '_og_rel_perm.dat'
      open(unit=86,file=string)
      write(86,*) '"oil_saturation", "' // trim(phase) // '_og_rel_perm", &
                  &"dkrog/dsato", " dkrog/dsatg"'
      do i = 1, size(saturation)
        write(86,'(4es14.6)') oil_saturation(i), krog(i), dkrog_sato(i), &
                              dkrog_satg(i)
      enddo
      close(86)
    class default
      do i = 1, num_values
        call this%RelativePermeability(oil_saturation(i),gas_saturation(i), &
                                       kr(i),dkr_sato(i),dkr_satg(i),option)
      enddo
      !print any other owg rel perms
      string = trim(cc_name) // '_' // trim(phase) // '_rel_perm.dat'
      open(unit=86,file=string)
      write(86,*) '"saturation", "' // trim(phase) // '_rel_perm", &
                  &"dkr/dsato", " dkr/dsatg"'
      do i = 1, size(saturation)
        write(86,'(4es14.6)') saturation(i), kr(i), dkr_sato(i), &
                              dkr_satg(i)
      enddo
      close(86)
  end select

end subroutine RPFOWGBaseTest

! ************************************************************************** !

! ************************************************************************** !

subroutine RPFOWGBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing OWG relative permeability functions

  use Option_module

  implicit none

  class(rel_perm_func_owg_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  option%io_buffer = 'RPF Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)

end subroutine RPFOWGBaseSetupPolynomials

! ************************************************************************** !

subroutine RPFOWGBaseRelPerm(this,oil_sat,gas_sat,rel_perm, &
                             dkr_sato,dkr_satg,option,table_idxs)
  use Option_module

  implicit none

  class(rel_perm_func_owg_base_type) :: this
  PetscReal, intent(in) :: oil_sat
  PetscReal, intent(in) :: gas_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkr_sato
  PetscReal, intent(out) :: dkr_satg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  option%io_buffer = 'RPFOWGBaseRelPerm must be extended.'
  call printErrMsg(option)

end subroutine RPFOWGBaseRelPerm

! ************************************************************************** !

subroutine RPF_OWG_MBC_Init(this)

  implicit none

  class(RPF_OWG_MBC_type) :: this

  call RPFOWGBaseInit(this)

  this%analytical_derivative_available = PETSC_TRUE

  this%m = UNINITIALIZED_DOUBLE
  !this%kr_max = 1.0d0

end subroutine RPF_OWG_MBC_Init

! ************************************************************************** !

subroutine RPF_OWG_MBC_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_MBC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OWG,RPF_OWG_MBC'
  endif

  call RPFOWGBaseVerify(this,string,option)

  if (Uninitialized(this%Swcr)) then
    option%io_buffer = UninitializedMessage('WATER_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif

  if (this%So_is_Sh) then
    if (Uninitialized(this%Socr)) then
      option%io_buffer = &
            UninitializedMessage('HYDROCARBON_RESIDUAL_SATURATION',string)
      call printErrMsg(option)
    endif
  else
    if (Uninitialized(this%Socr)) then
      option%io_buffer = UninitializedMessage('OIL_RESIDUAL_SATURATION',string)
      call printErrMsg(option)
    endif
    if (Uninitialized(this%Sgcr)) then
      option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
      call printErrMsg(option)
    endif
  end if

  ! by default kr_max = 1.0 - if entered the value must be bound beween 0 and 1
  !if ( this%kr_max < 0.0 .or. this%kr_max > 1.0 ) then
  !  option%io_buffer = string &
  !    // ' MAX_REL_PERM entered not valid, must be 0 <= Kr_max <= 1'
  !  call printErrMsg(option)
  !end if
  !if (Uninitialized(this%kr_max)) then
  !  option%io_buffer = UninitializedMessage('M',string)
  !  call printErrMsg(option)
  !endif

end subroutine RPF_OWG_MBC_Verify

! ************************************************************************** !

subroutine RPF_OWG_MBC_SetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing Modified BC permeability function

  use Option_module
  use Utility_module

  implicit none

  class(RPF_OWG_MBC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  PetscReal :: b(4)

  PetscReal :: Se_ph_low

  this%poly => PolynomialCreate()
  ! fill matix with values
  this%poly%low = 0.99d0  ! just below saturated
  !this%poly%low = 0.95d0  ! just below saturated
  this%poly%high = 1.d0   ! saturated
  Se_ph_low = this%poly%low

  b(1) = this%kr_max
  b(2) = this%kr_max * (Se_ph_low ** this%m)
  b(3) = 0.d0
  b(4) = this%m * this%kr_max * Se_ph_low ** (this%m - 1.0 )

  call CubicPolynomialSetup(this%poly%high,this%poly%low,b)

  this%poly%coefficients(1:4) = b(1:4)

end subroutine RPF_OWG_MBC_SetupPolynomials

! ************************************************************************** !

subroutine RPF_OWG_MBC_RelPerm_dkr_dSe(this,effective_sat,rel_perm,&
                                       dkr_Se,option)
 !
 ! Computes the relative permeability and its derivative WRT Se
 ! for the MBC RPF given the Se value
 !
 ! Author: Paolo Orsini (OGS)
 ! Date: 11/18/2017

  use Option_module
  use Utility_module

  implicit none

  class(RPF_OWG_MBC_type) :: this
  PetscReal, intent(in) :: effective_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkr_Se
  type(option_type), intent(inout) :: option

  PetscReal :: Se

  Se = effective_sat

  dkr_Se = 0.0d0

  if (Se >= 1.d0) then
    rel_perm = this%kr_max
    return
  else if (Se <=  0.d0) then
    rel_perm = 0.d0
    return
  endif

  if (associated(this%poly)) then
    if (Se > this%poly%low) then
      call CubicPolynomialEvaluate(this%poly%coefficients, &
                                   Se,rel_perm,dkr_Se)
      return
    endif
  endif

  rel_perm = this%kr_max * (Se ** this%m)

  dkr_Se = this%kr_max * this%m * (Se ** (this%m-1.0))


end subroutine RPF_OWG_MBC_RelPerm_dkr_dSe

! ************************************************************************** !

function RPF_wat_MBC_Create()

  implicit none

  class(RPF_OWG_MBC_wat_type), pointer :: RPF_wat_MBC_Create

  allocate(RPF_wat_MBC_Create)

  call RPF_wat_MBC_Create%Init()

end function RPF_wat_MBC_Create

! ************************************************************************** !

subroutine RPF_OWG_MBC_wat_RelPerm(this,oil_sat,gas_sat,rel_perm, &
                             dkr_sato,dkr_satg,option,table_idxs)
 !
 ! Computes the relative permeability (and associated derivatives) as a
 ! function of water saturation
 !
 ! Author: Paolo Orsini (OGS)
 ! Date: 11/18/2017

  use Option_module
  use Utility_module

  implicit none

  class(RPF_OWG_MBC_wat_type) :: this
  PetscReal, intent(in) :: oil_sat
  PetscReal, intent(in) :: gas_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkr_sato
  PetscReal, intent(out) :: dkr_satg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: wat_sat, Se, dSe_dSw, dkr_dSe


  wat_sat = 1.0 - oil_sat - gas_sat

  Se = (wat_sat - this%Swcr) / (1.d0 - this%Socr - this%Swcr - this%Sgcr )

  dSe_dSw = 1.0 / (1.d0 - this%Socr - this%Swcr - this%Sgcr )

  ! scaling WRT Kr_max occurs within RPF_OWG_MBC_RelPerm_dkr_dSe
  call this%RPF_OWG_MBC_RelPerm_dkr_dSe(Se,rel_perm,dkr_dSe,option)

  ! dkr_sato = dkr_dSe * dSe_dSw * dSw_dSo , with dSw_dSo = -1
  dkr_sato = - dkr_dSe * dSe_dSw

  ! dkr_satg = dkr_dSe * dSe_dSw * dSw_dSg , with dSw_dSg = -1
  dkr_satg = - dkr_dSe * dSe_dSw

end subroutine RPF_OWG_MBC_wat_RelPerm

! ************************************************************************** !

function RPF_oil_MBC_Create()

  implicit none

  class(RPF_OWG_MBC_oil_type), pointer :: RPF_oil_MBC_Create

  allocate(RPF_oil_MBC_Create)

  call RPF_oil_MBC_Create%Init()

end function RPF_oil_MBC_Create

! ************************************************************************** !

subroutine RPF_OWG_MBC_oil_RelPerm(this,oil_sat,gas_sat,rel_perm, &
                             dkr_sato,dkr_satg,option,table_idxs)
 !
 ! Computes the relative permeability (and associated derivatives) as a
 ! function of water saturation
 !
 ! Author: Paolo Orsini (OGS)
 ! Date: 11/18/2017

  use Option_module
  use Utility_module

  implicit none

  class(RPF_OWG_MBC_oil_type) :: this
  PetscReal, intent(in) :: oil_sat
  PetscReal, intent(in) :: gas_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkr_sato
  PetscReal, intent(out) :: dkr_satg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: Se, dkr_dSe, dSe_dSo

  Se = (oil_sat - this%Socr) / (1.d0 - this%Socr - this%Swcr - this%Sgcr )

  dSe_dSo = 1.0 / (1.d0 - this%Socr - this%Swcr - this%Sgcr )

  ! scaling WRT Kr_max occurs within RPF_OWG_MBC_RelPerm_dkr_dSe
  call this%RPF_OWG_MBC_RelPerm_dkr_dSe(Se,rel_perm,dkr_dSe,option)

  ! dkr_sato = dkr_dSe * dSe_dSo
  dkr_sato = dkr_dSe * dSe_dSo

  ! dkr_satg = dkr_dSe * dSe_dSo * dSodSg, with dSodSg = 1
  dkr_satg = - dkr_dSe * dSe_dSo

end subroutine RPF_OWG_MBC_oil_RelPerm

! ************************************************************************** !

function RPF_gas_MBC_Create()

  implicit none

  class(RPF_OWG_MBC_gas_type), pointer :: RPF_gas_MBC_Create

  allocate(RPF_gas_MBC_Create)

  call RPF_gas_MBC_Create%Init()

end function RPF_gas_MBC_Create

! ************************************************************************** !

subroutine RPF_OWG_MBC_gas_RelPerm(this,oil_sat,gas_sat,rel_perm, &
                               dkr_sato,dkr_satg,option,table_idxs)
 !
 ! Computes the relative permeability with a modified Brooks and Cory law
 ! (and associated derivatives) as function of water saturation
 !
 ! Author: Paolo Orsini (OGS)
 ! Date: 11/18/2017

  use Option_module
  use Utility_module

  implicit none

  class(RPF_OWG_MBC_gas_type) :: this
  PetscReal, intent(in) :: oil_sat
  PetscReal, intent(in) :: gas_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkr_sato
  PetscReal, intent(out) :: dkr_satg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: Se, dkr_dSe, dSe_dSg

  Se = (gas_sat - this%Sgcr) / (1.d0 - this%Socr - this%Swcr - this%Sgcr )

  dSe_dSg = 1.0 / (1.d0 - this%Socr - this%Swcr - this%Sgcr )

  ! scaling WRT Kr_max occurs within RPF_OWG_MBC_RelPerm_dkr_dSe
  call this%RPF_OWG_MBC_RelPerm_dkr_dSe(Se,rel_perm,dkr_dSe,option)

  ! dkr_sato = dkr_dSe * dSe_dSg * dSgdSo, with dSgdSo = -1
  dkr_sato = - dkr_dSe * dSe_dSg

  ! dkr_satg = dkr_dSe * dSe_dSg
  dkr_satg = dkr_dSe * dSe_dSg

end subroutine RPF_OWG_MBC_gas_RelPerm

! ************************************************************************** !

subroutine RPF_OWG_func_sl_Init(this)

  implicit none

  class(RPF_OWG_func_sl_type) :: this

  call RPFOWGBaseInit(this)

  this%analytical_derivative_available = PETSC_TRUE

  !this%kr_max = 1.0d0

end subroutine RPF_OWG_func_sl_Init

! ************************************************************************** !

subroutine RPF_OWG_func_sl_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_func_sl_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  call RPFOWGBaseVerify(this,string,option)

  !if ( this%kr_max < 0.0 .or. this%kr_max > 1.0 ) then
  !  option%io_buffer = name &
  !    // ' MAX_REL_PERM entered not valid, must be 0 <= Kr_max <= 1'
  !  call printErrMsg(option)
  !end if

  if (.not.associated(this%rel_perm_func_sl) ) then
    option%io_buffer = name &
      // ' Sub SL analytical model not defined'
    call printErrMsg(option)
  else
    call this%rel_perm_func_sl%verify(name,option)
  end if

end subroutine RPF_OWG_func_sl_Verify

! ************************************************************************** !

subroutine RPF_OWG_func_sl_RelPerm(this,oil_sat,gas_sat,rel_perm, &
                                         dkr_sato,dkr_satg,option,table_idxs)
!
! Computes the relative permeability (and associated derivatives) using
! analytical models implemented as function of sl
!
! Author: Paolo Orsini (OGS)
! Date: 11/20/2017

 use Option_module
 use Utility_module

 implicit none

 class(RPF_OWG_func_sl_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(in) :: gas_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkr_sato
 PetscReal, intent(out) :: dkr_satg
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 PetscReal :: func_sat, dkr_func_sat

 if ( this%function_of_liquid_sat ) then
   func_sat = 1.0 - gas_sat
 else
   func_sat = 1.0 - oil_sat - gas_sat
 end if

 call this%rel_perm_func_sl%RelativePermeability(func_sat,rel_perm, &
                                                 dkr_func_sat,option)

 rel_perm = this%kr_max * rel_perm

 dkr_satg = - this%kr_max *dkr_func_sat

 if ( this%function_of_liquid_sat ) then
   dkr_sato = this%kr_max * dkr_func_sat
 else
   dkr_sato = - this%kr_max * dkr_func_sat
 end if

end subroutine RPF_OWG_func_sl_relPerm

! ************************************************************************** !

subroutine RPF_OWG_func_sl_VG_Init(this)

  implicit none

  class(RPF_OWG_func_sl_VG_type) :: this

  call RPF_OWG_func_sl_Init(this)

  this%m = UNINITIALIZED_DOUBLE

end subroutine RPF_OWG_func_sl_VG_Init

! ************************************************************************** !

subroutine RPF_OWG_func_sl_VG_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_func_sl_VG_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  call RPF_OWG_func_sl_Verify(this,name,option)

  if (this%function_of_liquid_sat) then
    if (Uninitialized(this%Slcr)) then
      option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                               name)
      call printErrMsg(option)
    endif
  else
    if (Uninitialized(this%Swcr)) then
      option%io_buffer = UninitializedMessage('WATER_RESIDUAL_SATURATION',name)
      call printErrMsg(option)
    endif
  end if

  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',name)
    call printErrMsg(option)
  endif

end subroutine RPF_OWG_func_sl_VG_Verify

! ************************************************************************** !

function RPF_OWG_Mualem_VG_wat_Create()

  implicit none

  class(RPF_OWG_Mualem_VG_wat_type), pointer :: RPF_OWG_Mualem_VG_wat_Create

  allocate(RPF_OWG_Mualem_VG_wat_Create)

  call RPF_OWG_Mualem_VG_wat_Create%Init()

  RPF_OWG_Mualem_VG_wat_Create%rel_perm_func_sl => RPF_Mualem_VG_liq_create()

end function RPF_OWG_Mualem_VG_wat_Create

! ************************************************************************** !

subroutine RPF_OWG_Mualem_VG_wat_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_Mualem_VG_wat_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OWG,RPF_OWG_Mualem_VG_wat'
  endif

  call RPF_OWG_func_sl_VG_Verify(this,string,option)

end subroutine RPF_OWG_Mualem_VG_wat_Verify

! ************************************************************************** !

function RPF_OWG_Mualem_VG_gas_Create()

  implicit none

  class(RPF_OWG_Mualem_VG_gas_type), pointer :: RPF_OWG_Mualem_VG_gas_Create

  allocate(RPF_OWG_Mualem_VG_gas_Create)

  call RPF_OWG_Mualem_VG_gas_Create%Init()

  RPF_OWG_Mualem_VG_gas_Create%rel_perm_func_sl => RPF_Mualem_VG_gas_create()

end function RPF_OWG_Mualem_VG_gas_Create

! ************************************************************************** !

subroutine RPF_OWG_Mualem_VG_gas_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_Mualem_VG_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OWG,RPF_OWG_Mualem_VG_gas'
  endif

  if ( .not. this%function_of_liquid_sat) then
    option%io_buffer = string // ' only supported as function of Liquid Sat'
  end if

  call RPF_OWG_func_sl_VG_Verify(this,string,option)

  if (Uninitialized(this%Sgcr)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif

end subroutine RPF_OWG_Mualem_VG_gas_Verify

! ************************************************************************** !

function RPF_OWG_TOUGH2_IRP7_gas_Create()

  implicit none

  class(RPF_OWG_TOUGH2_IRP7_gas_type), pointer :: &
                                            RPF_OWG_TOUGH2_IRP7_gas_Create

  allocate(RPF_OWG_TOUGH2_IRP7_gas_Create)

  call RPF_OWG_TOUGH2_IRP7_gas_Create%Init()

  RPF_OWG_TOUGH2_IRP7_gas_Create%rel_perm_func_sl => &
                                      RPF_TOUGH2_IRP7_gas_create()

end function RPF_OWG_TOUGH2_IRP7_gas_Create

! ************************************************************************** !

subroutine RPF_OWG_TOUGH2_IRP7_gas_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_TOUGH2_IRP7_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OWG,RPF_OWG_TOUGH2_IRP7_gas'
  endif

  if ( .not. this%function_of_liquid_sat) then
    option%io_buffer = string // ' only supported as function of Liquid Sat'
  end if

  call RPF_OWG_func_sl_VG_Verify(this,string,option)

  if (Uninitialized(this%Sgcr)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif

end subroutine RPF_OWG_TOUGH2_IRP7_gas_Verify

! ************************************************************************** !

function RPF_OWG_Burdine_VG_wat_Create()

  implicit none

  class(RPF_OWG_Burdine_VG_wat_type), pointer :: RPF_OWG_Burdine_VG_wat_Create

  allocate(RPF_OWG_Burdine_VG_wat_Create)

  call RPF_OWG_Burdine_VG_wat_Create%Init()

  RPF_OWG_Burdine_VG_wat_Create%rel_perm_func_sl => RPF_Burdine_VG_liq_create()

end function RPF_OWG_Burdine_VG_wat_Create

! ************************************************************************** !

subroutine RPF_OWG_Burdine_VG_wat_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_Burdine_VG_wat_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OWG,RPF_OWG_Burdine_VG_wat'
  endif

  call RPF_OWG_func_sl_VG_Verify(this,string,option)

end subroutine RPF_OWG_Burdine_VG_wat_Verify

! ************************************************************************** !

function RPF_OWG_Burdine_VG_gas_Create()

  implicit none

  class(RPF_OWG_Burdine_VG_gas_type), pointer :: RPF_OWG_Burdine_VG_gas_Create

  allocate(RPF_OWG_Burdine_VG_gas_Create)

  call RPF_OWG_Burdine_VG_gas_Create%Init()

  RPF_OWG_Burdine_VG_gas_Create%rel_perm_func_sl => RPF_Burdine_VG_gas_create()

end function RPF_OWG_Burdine_VG_gas_Create

! ************************************************************************** !

subroutine RPF_OWG_Burdine_VG_gas_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_Burdine_VG_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OWG,RPF_OWG_Burdine_VG_gas'
  endif

  if ( .not. this%function_of_liquid_sat) then
    option%io_buffer = string // ' only supported as function of Liquid Sat'
  end if

  call RPF_OWG_func_sl_VG_Verify(this,string,option)

  if (Uninitialized(this%Sgcr)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif

end subroutine RPF_OWG_Burdine_VG_gas_Verify

! ************************************************************************** !

subroutine RPF_OWG_func_sl_BC_Init(this)

  implicit none

  class(RPF_OWG_func_sl_BC_type) :: this

  call RPF_OWG_func_sl_Init(this)

  this%lambda = UNINITIALIZED_DOUBLE

end subroutine RPF_OWG_func_sl_BC_Init

! ************************************************************************** !

subroutine RPF_OWG_func_sl_BC_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_func_sl_BC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  call RPF_OWG_func_sl_Verify(this,name,option)

  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',name)
    call printErrMsg(option)
  endif

  if (this%function_of_liquid_sat) then
    if (Uninitialized(this%Slcr)) then
      option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                               name)
      call printErrMsg(option)
    endif
  else
    if (Uninitialized(this%Swcr)) then
      option%io_buffer = UninitializedMessage('WATER_RESIDUAL_SATURATION',name)
      call printErrMsg(option)
    endif
  end if

end subroutine RPF_OWG_func_sl_BC_Verify

! ************************************************************************** !

function RPF_OWG_Burdine_BC_wat_Create()

  implicit none

  class(RPF_OWG_Burdine_BC_wat_type), pointer :: RPF_OWG_Burdine_BC_wat_Create

  allocate(RPF_OWG_Burdine_BC_wat_Create)

  call RPF_OWG_Burdine_BC_wat_Create%Init()

  RPF_OWG_Burdine_BC_wat_Create%rel_perm_func_sl => RPF_Burdine_BC_liq_create()

end function RPF_OWG_Burdine_BC_wat_Create

! ************************************************************************** !

subroutine RPF_OWG_Burdine_BC_wat_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_Burdine_BC_wat_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // &
         'PERMEABILITY_FUNCTION_OWG,RPF_OWG_Burdine_BC_wat_Verify'
  endif

  call RPF_OWG_func_sl_BC_Verify(this,string,option)

end subroutine RPF_OWG_Burdine_BC_wat_Verify

! ************************************************************************** !

function RPF_OWG_Burdine_BC_gas_Create()

  implicit none

  class(RPF_OWG_Burdine_BC_gas_type), pointer :: RPF_OWG_Burdine_BC_gas_Create

  allocate(RPF_OWG_Burdine_BC_gas_Create)

  call RPF_OWG_Burdine_BC_gas_Create%Init()

  RPF_OWG_Burdine_BC_gas_Create%rel_perm_func_sl => RPF_Burdine_BC_gas_create()

end function RPF_OWG_Burdine_BC_gas_Create

! ************************************************************************** !

subroutine RPF_OWG_Burdine_BC_gas_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_OWG_Burdine_BC_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // &
         'PERMEABILITY_FUNCTION_OWG,RPF_OWG_Burdine_BC_gas_Verify'
  endif

  if ( .not. this%function_of_liquid_sat) then
    option%io_buffer = string // ' only supported as function of Liquid Sat'
  end if

  call RPF_OWG_func_sl_BC_Verify(this,string,option)

  if (Uninitialized(this%Sgcr)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif

end subroutine RPF_OWG_Burdine_BC_gas_Verify

! ************************************************************************** !

function RPF_oil_ecl_Create()

  implicit none

  class(RPF_oil_ecl_type), pointer :: RPF_oil_ecl_Create

  allocate(RPF_oil_ecl_Create)

  call RPF_oil_ecl_Create%Init()

end function RPF_oil_ecl_Create

! ************************************************************************** !

subroutine RPF_oil_ecl_Init(this)

  implicit none

  class(RPF_oil_ecl_type) :: this

  call RPFOWGBaseInit(this)

  nullify(this%rel_perm_ow)
  nullify(this%rel_perm_og)
  !this%kr_max = 1.0

end subroutine RPF_oil_ecl_Init

! ************************************************************************** !

subroutine RPF_oil_ecl_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_oil_ecl_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string, string_ow, string_og

  if (index(name,'PERMEABILITY_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OWG,RPF_OIL_ECL'
  endif

  call RPFOWGBaseVerify(this,string,option)

  if ( .not.associated(this%rel_perm_ow) ) then
    string_ow = string // ' ,RPF_OIL_WATER not defined'
    option%io_buffer = string_ow
    call printErrMsg(option)
  end if

  if ( .not.associated(this%rel_perm_og) ) then
    string_og = string // ' ,RPF_OIL_GAS not defined'
    option%io_buffer = string_og
    call printErrMsg(option)
  end if

  ! PO: TODO if ow and/or og function crtitical saturation not defined
  ! take the values defined in the ECLIPSE oil rpf
  if (Uninitialized(this%Swcr)) then
    this%Swcr=0.0
    !option%io_buffer = UninitializedMessage('WATER_RESIDUAL_SATURATION',string)
    !call printErrMsg(option)
  else
    if (Uninitialized(this%rel_perm_ow%Swcr)) then
       this%rel_perm_ow%Swcr = this%Swcr
    end if
    if (Uninitialized(this%rel_perm_og%Swcr)) then
      this%rel_perm_og%Swcr = this%Swcr
    end if
  end if

  if (Uninitialized(this%Socr)) then
    this%Socr = 0.0
    !option%io_buffer = UninitializedMessage('OIL_RESIDUAL_SATURATION',string)
    !call printErrMsg(option)
  else
    if (Uninitialized(this%rel_perm_ow%Socr)) then
       this%rel_perm_ow%Socr = this%Socr
    end if
    if (Uninitialized(this%rel_perm_og%Socr)) then
      this%rel_perm_og%Socr = this%Socr
    end if
  end if

  if (Uninitialized(this%Sgcr)) then
    this%Sgcr = 0.0
    !option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    !call printErrMsg(option)
  else
    if (Uninitialized(this%rel_perm_ow%Sgcr)) then
       this%rel_perm_ow%Sgcr = this%Sgcr
    end if
    if (Uninitialized(this%rel_perm_og%Sgcr)) then
      this%rel_perm_og%Sgcr = this%Sgcr
    end if
  end if

  !string_ow = string //'RPF_OIL_WATER, choose a kr function that supports krmax'
  !string_og = string //'RPF_OIL_GAS, choose a kr function that supports krmax'

  !if ( this%kr_max < 0.0 .or. this%kr_max > 1.0 ) then
  !  option%io_buffer = string &
  !    // ' MAX_REL_PERM entered not valid, must be 0 <= Kr_max <= 1'
  !  call printErrMsg(option)
  !else
  !  this%rel_perm_ow%kr_max = this%kr_max
  !  this%rel_perm_og%kr_max = this%kr_max
    ! select type(rel_perm_ow => this%rel_perm_ow)
    !   class is(RPF_OWG_MBC_oil_type)
    !     rel_perm_ow%kr_max = this%kr_max
    !   class default
    !     option%io_buffer = string_ow
    !     call printErrMsg(option)
    ! end select
    ! select type(rel_perm_og => this%rel_perm_og)
    !   class is(RPF_OWG_MBC_oil_type)
    !     rel_perm_og%kr_max = this%kr_max
    !   class default
    !     option%io_buffer = string_og
    !     call printErrMsg(option)
    ! end select
  !end if

  ! pass connate water, gas and oil to ow and og functions
  this%rel_perm_ow%Swco = this%Swco
  this%rel_perm_ow%Soco = this%Soco
  this%rel_perm_ow%Sgco = this%Sgco
  this%rel_perm_og%Swco = this%Swco
  this%rel_perm_og%Soco = this%Soco
  this%rel_perm_og%Sgco = this%Sgco

  string_ow = string // ',RPF_OIL_WATER'

  call this%rel_perm_ow%Verify(string_ow,option)

  string_og = string // ',RPF_OIL_GAS'

  call this%rel_perm_og%Verify(string_og,option)

  if (this%rel_perm_ow%kr_max /= this%rel_perm_og%kr_max) then
    option%io_buffer = trim(string) &
      // ' MAX_REL_PERM different in Krow and krog'
    call printErrMsg(option)
  end if

end subroutine RPF_oil_ecl_Verify

! ************************************************************************** !

subroutine RPF_oil_ecl_RelPerm(this,oil_sat,gas_sat,rel_perm, &
                               dkr_sato,dkr_satg,option,table_idxs)
 !
 ! Computes the oil relative permeability and associated derivatives
 ! using the eclipse default model
 !
 ! Author: Paolo Orsini (OGS)
 ! Date: 11/18/2017

  use Option_module
  use Utility_module

  implicit none

  class(RPF_oil_ecl_type) :: this
  PetscReal, intent(in) :: oil_sat
  PetscReal, intent(in) :: gas_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkr_sato
  PetscReal, intent(out) :: dkr_satg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: Krow, Krog
  PetscReal :: dKrow_sato, dKrow_satg
  PetscReal :: dKrog_sato, dKrog_satg
  PetscReal :: wat_sat
  PetscReal :: sg_plus_sw_minus_swco
  !PetscReal :: den_eps = 1.0d-8

  wat_sat = 1.0 - oil_sat - gas_sat

  call this%rel_perm_ow%RelativePermeability(oil_sat,gas_sat,Krow,dKrow_sato, &
                                             dKrow_satg,option,table_idxs)

  call this%rel_perm_og%RelativePermeability(oil_sat,gas_sat,Krog,dKrog_sato, &
                                             dKrog_satg,option,table_idxs)

  sg_plus_sw_minus_swco = gas_sat + wat_sat - this%Swco

  if ( sg_plus_sw_minus_swco > 0.0 ) then
    rel_perm = (gas_sat * Krog +  (wat_sat - this%Swco) * Krow ) / &
               sg_plus_sw_minus_swco
  else
    rel_perm = this%kr_max
  end if

end subroutine RPF_oil_ecl_RelPerm

! ************************************************************************** !
! *********** END OWG Relative Permeability functions  ********************* !
! ************************************************************************** !

! ************************************************************************** !

subroutine SaturationFunctionOWGDestroy(sf_owg)
  !
  ! Destroys an OWG saturuation function
  !
  ! Author: Paolo Orsini
  ! Date: 11/16/17
  !

  implicit none

  class(sat_func_owg_base_type), pointer :: sf_owg

  if (.not.associated(sf_owg)) return

  call SaturationFunctionDestroy(sf_owg%sat_func_sl)

  deallocate(sf_owg)
  nullify(sf_owg)

end subroutine SaturationFunctionOWGDestroy

! ************************************************************************** !

recursive subroutine PermeabilityFunctionOWGDestroy(rpf)
  !
  ! Destroys an OWG permeability function
  !
  ! Author: Paolo Orsini
  ! Date: 11/18/17
  !

  implicit none

  class(rel_perm_func_owg_base_type), pointer :: rpf

  if (.not.associated(rpf)) return

  call PolynomialDestroy(rpf%poly)

  call PermeabilityFunctionDestroy(rpf%rel_perm_func_sl)

  select type(rpf)
    class is(RPF_oil_ecl_type)
      call PermeabilityFunctionOWGDestroy(rpf%rel_perm_ow)
      call PermeabilityFunctionOWGDestroy(rpf%rel_perm_og)
  end select

  deallocate(rpf)
  nullify(rpf)

end subroutine PermeabilityFunctionOWGDestroy

end module Characteristic_Curves_OWG_module
