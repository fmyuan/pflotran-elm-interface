module Characteristic_Curves_OWG_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_Common_module
  use Characteristic_Curves_WIPP_module
  use Characteristic_Curves_Table_module

  implicit none

  private
  
!-----------------------------------------------------------------------------
!-- OWG Saturation Functions Base class  -------------------------------------
!-----------------------------------------------------------------------------

  type, abstract, public :: sat_func_owg_base_type
    PetscReal :: pcmax
    PetscReal :: pcmin
    class(sat_func_base_type), pointer :: sat_func_sl
    PetscBool :: analytical_derivative_available
    PetscBool :: sat_func_of_pc_available
    character(len=MAXWORDLENGTH) :: table_name
    class(char_curves_table_type), pointer :: table
  contains
    procedure, public :: Init => SFOWGBaseInit !only for pcmax, nullify pointer (table and sl), init table name
    procedure, public :: Verify => SFOWGBaseVerify ! to be extended
    procedure, public :: Test => SFOWGBaseTest     ! to be extended
    procedure, public :: SetupPolynomials => SFOWGBaseSetupPolynomials
    procedure, public :: ProcessTable => SFOWGBaseProcessTable
    procedure, public :: GetPcMax
    procedure, public :: GetPcMin
  end type sat_func_owg_base_type

  !-----------------------------------------------------------------------------
  !-- XW Saturation Functions  -------------------------------------------------
  !-- XW = x/water - any two-phase interfaces where water is the wetting phase -
  !-----------------------------------------------------------------------------

  !type, public :: sat_func_xw_base_type
  type,abstract,public,extends(sat_func_owg_base_type) :: sat_func_xw_base_type
    ! xw =  x/water - any two-phase interfaces where water is wetting
    PetscReal :: Swco !connate water saturation
    PetscReal :: Swcr !critical (residual) water saturation
  contains
    procedure, public :: Init => SFXWBaseInit
    procedure, public :: Verify => SFXWBaseVerify
    procedure, public :: Test => SFXWBaseTest
    !procedure, public :: SetupPolynomials => SFXWBaseSetupPolynomials
    procedure, public :: CapillaryPressure => SFXWBaseCapillaryPressure
    procedure, public :: Saturation => SFXWBaseSaturation
    procedure, public :: SetConnateSaturation => SFXWSetConnateSatBase
    procedure, public :: GetConnateSaturation => SFXWGetConnateSatBase
    procedure, public :: SetCriticalSaturation => SFXWSetCriticalSatBase
    procedure, public :: GetCriticalSaturation => SFXWGetCriticalSatBase
    procedure, public :: ComputePcMin => SFXWComputePcMin
  end type sat_func_xw_base_type  
  
  type, public, extends(sat_func_xw_base_type) :: sat_func_xw_constant_type
    PetscReal :: constant_capillary_pressure
  contains
    procedure, public :: Init => SF_XW_constant_Init
    procedure, public :: Verify => SF_XW_constant_Verify
    procedure, public :: CapillaryPressure => SF_XW_const_CapillaryPressure
    procedure, public :: Saturation => SF_XW_const_Saturation
  end type sat_func_xw_constant_type

  ! add other analytical model extensions  
  type, public, extends(sat_func_xw_base_type) :: sat_func_xw_VG_type
    PetscReal :: alpha
    PetscReal :: m
  contains
    procedure, public :: Init => SF_XW_VG_Init
    procedure, public :: Verify => SF_XW_VG_Verify
    procedure, public :: CapillaryPressure => SF_XW_VG_CapillaryPressure
    procedure, public :: Saturation => SF_XW_VG_Saturation
  end type sat_func_xw_VG_type

  type, public, extends(sat_func_xw_base_type) :: sat_func_xw_BC_type
    PetscReal :: alpha
    PetscReal :: lambda
  contains
    procedure, public :: Init => SF_XW_BC_Init
    procedure, public :: Verify => SF_XW_BC_Verify
    procedure, public :: CapillaryPressure => SF_XW_BC_CapillaryPressure
    procedure, public :: Saturation => SF_XW_BC_Saturation
  end type sat_func_xw_BC_type

   type, public, extends(sat_func_xw_base_type) :: sat_func_xw_table_type
   contains
     procedure, public :: Init => SF_XW_table_Init
     !procedure, public :: Verify => SF_XW_table_Verify
     procedure, public :: CapillaryPressure => SF_XW_table_CapillaryPressure
     procedure, public :: Saturation => SF_XW_table_Saturation
   end type sat_func_xw_table_type

!-----------------------------------------------------------------------------
!-- OG Saturation Functions --------------------------------------------------
!-----------------------------------------------------------------------------

  type,abstract,public,extends(sat_func_owg_base_type) :: sat_func_og_base_type
    PetscReal :: Sgco !connate oil saturation
    PetscReal :: Sgcr !critical (residual) gas saturation
  contains
    procedure, public :: Init => SFOGBaseInit
    procedure, public :: Verify => SFOGBaseVerify
    procedure, public :: Test => SFOGBaseTest
    procedure, public :: SetupPolynomials => SFOGBaseSetupPolynomials
    procedure, public :: CapillaryPressure => SFOGBaseCapillaryPressure
    procedure, public :: Saturation => SFOGBaseSaturation
    procedure, public :: SetConnateSaturation => SFOGSetConnateSatBase
    procedure, public :: GetConnateSaturation => SFOGGetConnateSatBase
    procedure, public :: SetCriticalSaturation => SFOGSetCriticalSatBase
    procedure, public :: GetCriticalSaturation => SFOGGetCriticalSatBase
    procedure, public :: ComputePcMin => SFOGComputePcMin
  end type sat_func_og_base_type

  type,public,extends(sat_func_og_base_type) :: sat_func_og_constant_type
    PetscReal :: constant_capillary_pressure
  contains
    procedure, public :: Init => SF_OG_constant_Init
    procedure, public :: Verify => SF_OG_constant_Verify
    procedure, public :: CapillaryPressure => SF_OG_Const_CapillaryPressure
  end type sat_func_og_constant_type

  type, public, extends(sat_func_og_base_type) :: sat_func_og_VG_SL_type
    PetscReal :: alpha
    PetscReal :: m
    PetscReal :: Slcr
  contains
    procedure, public :: Init => SF_OG_VG_SL_Init
    procedure, public :: Verify => SF_OG_VG_SL_Verify
    procedure, public :: CapillaryPressure => SF_OG_VG_SL_CapillaryPressure
    !procedure, public :: Saturation => SF_OW_VG_Saturation
  end type sat_func_og_VG_SL_type

  type, public, extends(sat_func_og_base_type) :: sat_func_og_table_type
  contains
    procedure, public :: Init => SF_OG_table_Init
    !procedure, public :: Verify => SF_OG_table_Verify
    procedure, public :: CapillaryPressure => SF_OG_table_CapillaryPressure
    procedure, public :: Saturation => SF_OG_table_Saturation
  end type sat_func_og_table_type

  !-----------------------------------------------------------------------------
  !-- OWG Relative Permeability Functions --------------------------------------
  !-----------------------------------------------------------------------------  
  
  type,abstract,public :: rel_perm_owg_base_type
    PetscReal :: m
    PetscReal :: kr_max
    class(rel_perm_func_base_type), pointer :: rel_perm_func_sl
    type(polynomial_type), pointer :: poly
    PetscBool :: analytical_derivative_available
    character(len=MAXWORDLENGTH) :: table_name
    class(char_curves_table_type), pointer :: table
  contains
    procedure, public :: Init => RPFOWGBaseInit
    procedure, public :: Verify => RPFOWGBaseVerify
    procedure, public :: SetupPolynomials => RPFOWGBaseSetupPolynomials !to be extended
    procedure, public :: Strip => PermFunctionOWGBaseStrip
    procedure, public :: ProcessTable => RPFOWGBaseProcessTable
  end type rel_perm_owg_base_type

  !-----------------------------------------------------------------------------
  !-- Water OWG Relative Permeability Functions --------------------------------
  !-----------------------------------------------------------------------------

  type,abstract,public,extends(rel_perm_owg_base_type) :: &
                                                  rel_perm_wat_owg_base_type
    PetscReal :: Swco !connate water saturation
    PetscReal :: Swcr !critical water saturation
  contains
    procedure, public :: Init => RPFWatOWGBaseInit
    procedure, public :: Verify => RPFWatOWGBaseVerify
    procedure, public :: Test => RPFWatOWGBaseTest
    procedure, public :: SetupPolynomials => RPFWatOWGBaseSetupPolynomials
    procedure, public :: RelativePermeability => RPFWatOWGBaseRelPerm !defines argument template
    procedure, public :: GetConnateSaturation => RPFWatGetConnateSatOWGBase
    procedure, public :: GetCriticalSaturation => RPFWatGetCriticalSatOWGBase
  end type rel_perm_wat_owg_base_type

  type,public,extends(rel_perm_wat_owg_base_type) :: RPF_wat_owg_MBC_type
    PetscReal :: Sowcr !critical oil saturation
  contains
    procedure, public :: Init => RPF_wat_owg_MBC_Init
    procedure, public :: Verify => RPF_wat_owg_MBC_Verify
    !procedure, public :: Test => RPFOWGWatBaseTest
    procedure, public :: SetupPolynomials => RPF_wat_owg_MBC_SetupPoly
    procedure, public :: RPF_wat_owg_MBC_SetSowcr
    procedure, public :: RelativePermeability => RPF_wat_owg_MBC_RelPerm !define argument template
  end type RPF_wat_owg_MBC_type

  type,abstract,public,extends(rel_perm_wat_owg_base_type) :: &
                                                      RPF_wat_owg_func_sl_type
  contains  
    procedure, public :: Init => RPF_wat_owg_func_sl_Init
    procedure, public :: Verify => RPF_wat_owg_func_sl_Verify
    procedure, public :: RelativePermeability => RPF_wat_owg_func_sl_RelPerm
  end type RPF_wat_owg_func_sl_type

  type, public, extends(RPF_wat_owg_func_sl_type) :: &
                                                  RPF_wat_owg_Mualem_VG_type
  contains                                                
  end type RPF_wat_owg_Mualem_VG_type

  type, public, extends(RPF_wat_owg_func_sl_type) :: &
                                                  RPF_wat_owg_Burdine_VG_type                
  contains
  end type RPF_wat_owg_Burdine_VG_type

  type, public, extends(RPF_wat_owg_func_sl_type) :: &
                                                RPF_wat_owg_Burdine_BC_type
    PetscReal :: lambda                                            
  contains
  end type RPF_wat_owg_Burdine_BC_type

  type,public,extends(rel_perm_wat_owg_base_type) :: RPF_wat_owg_table_type
  contains
    procedure, public :: Init => RPF_wat_owg_table_Init
    !procedure, public :: Verify => RPF_wat_owg_table_Verify
    procedure, public :: RelativePermeability => RPF_wat_owg_table_RelPerm
  end type RPF_wat_owg_table_type

  !-----------------------------------------------------------------------------
  !-- Gas OWG Relative Permeability Functions ----------------------------------
  !-----------------------------------------------------------------------------
  type,abstract,public,extends(rel_perm_owg_base_type) :: &
                                                  rel_perm_gas_owg_base_type
    PetscReal :: Sgco !Gas connate saturation
    PetscReal :: Sgcr !critical gas saturation
  contains
    procedure, public :: Init => RPFGasOWGBaseInit
    procedure, public :: Verify => RPFGasOWGBaseVerify
    procedure, public :: Test => RPFGasOWGBaseTest
    procedure, public :: SetupPolynomials => RPFGasOWGBaseSetupPolynomials
    procedure, public :: RelativePermeability => RPFGasOWGBaseRelPerm !defines argument template
    procedure, public :: GetConnateSaturation => RPFGasGetConnateSatOWGBase
    procedure, public :: GetCriticalSaturation => RPFGasGetCriticalSatOWGBase
  end type rel_perm_gas_owg_base_type

  type,public,extends(rel_perm_gas_owg_base_type) :: RPF_gas_owg_MBC_type
    PetscReal :: Swco !connate water saturation
    PetscReal :: Sogcr !critical saturation of oil in gas
  contains
    procedure, public :: Init => RPF_gas_owg_MBC_Init
    procedure, public :: Verify => RPF_gas_owg_MBC_Verify
    procedure, public :: SetupPolynomials => RPF_gas_owg_MBC_SetupPoly 
    procedure, public :: RPF_gas_owg_MBC_SetSwcoSogcr
    procedure, public :: RelativePermeability => RPF_gas_owg_MBC_RelPerm !define argument template
  end type RPF_gas_owg_MBC_type

  type,abstract,public,extends(rel_perm_gas_owg_base_type) :: &
                                                      RPF_gas_owg_func_sl_type
    PetscReal :: Slcr !critical saturation of liquid
  contains  
    procedure, public :: Init => RPF_gas_owg_func_sl_Init
    procedure, public :: Verify => RPF_gas_owg_func_sl_Verify
    procedure, public :: RelativePermeability => RPF_gas_owg_func_sl_RelPerm
  end type RPF_gas_owg_func_sl_type
 
  type, public, extends(RPF_gas_owg_func_sl_type) :: RPF_gas_owg_Mualem_VG_type
  contains
  end type RPF_gas_owg_Mualem_VG_type

  type, public, extends(RPF_gas_owg_func_sl_type) :: &
                                                 RPF_gas_owg_TOUGH2_IRP7_type
  contains
  end type RPF_gas_owg_TOUGH2_IRP7_type

  type, public, extends(RPF_gas_owg_func_sl_type) :: &
                                                 RPF_gas_owg_Burdine_VG_type
  contains
  end type RPF_gas_owg_Burdine_VG_type

  type, public, extends(RPF_gas_owg_func_sl_type) :: &
                                                 RPF_gas_owg_Burdine_BC_type
    PetscReal :: lambda                                             
  contains
  end type RPF_gas_owg_Burdine_BC_type

  type,public,extends(rel_perm_gas_owg_base_type) :: &
                                                      RPF_gas_owg_table_type
  contains  
    procedure, public :: Init => RPF_gas_owg_table_Init
    !procedure, public :: Verify => RPF_gas_owg_table_Verify
    procedure, public :: RelativePermeability => RPF_gas_owg_table_RelPerm
  end type RPF_gas_owg_table_type

  !-----------------------------------------------------------------------------
  !-- Oil-Water OWG Relative Permeability Functions ----------------------------
  !-----------------------------------------------------------------------------
  ! two options: MBC and table
  
  type,abstract,public,extends(rel_perm_owg_base_type) :: &
                                                  rel_perm_ow_owg_base_type 
    PetscReal :: Soco !oil connate saturation
    PetscReal :: Sowcr !critical oil-water saturation
    PetscBool :: So_is_Sh
  contains
    procedure, public :: Init => RPFOWOWGBaseInit
    procedure, public :: Verify => RPFOWOWGBaseVerify
    procedure, public :: Test => RPFOWOWGBaseTest
    procedure, public :: SetupPolynomials => RPFOWOWGBaseSetupPoly
    procedure, public :: RelativePermeability => RPFOWOWGBaseRelPerm !defines argument template
    procedure, public :: GetCriticalSaturation => RPFOWGetCriticalSatOWGBase
    procedure, public :: GetConnateSaturation => RPFOWGetConnateSaturation
  end type rel_perm_ow_owg_base_type  

  type,public,extends(rel_perm_ow_owg_base_type) :: rel_perm_ow_owg_linear_type
  contains
    procedure, public :: Init => RPF_ow_owg_linear_Init
    procedure, public :: Verify => RPF_ow_owg_linear_Verify
    procedure, public :: RelativePermeability => RPF_ow_owg_linear_RelPerm
  end type rel_perm_ow_owg_linear_type
    
  type,public,extends(rel_perm_ow_owg_base_type) :: rel_perm_ow_owg_MBC_type
    PetscReal :: Swcr
  contains
    procedure, public :: Init => RPF_ow_owg_MBC_Init
    procedure, public :: Verify => RPF_ow_owg_MBC_Verify
    procedure, public :: SetupPolynomials => RPF_ow_owg_MBC_SetupPoly
    procedure, public :: RelativePermeability => RPF_ow_owg_MBC_RelPerm
    procedure, public :: RPF_ow_owg_MBC_SetSwcr
  end type rel_perm_ow_owg_MBC_type

  type,public,extends(rel_perm_ow_owg_base_type) :: rel_perm_ow_owg_table_type
  contains
    procedure, public :: Init => RPF_ow_owg_table_Init
    !procedure, public :: Verify => RPF_ow_owg_table_Verify
    procedure, public :: RelativePermeability => RPF_ow_owg_table_RelPerm
  end type rel_perm_ow_owg_table_type
 
  !-----------------------------------------------------------------------------
  !-- Oil-Gas OWG Relative Permeability Functions ------------------------------
  !-----------------------------------------------------------------------------
  ! two options: MBC and table
  type,abstract,public,extends(rel_perm_owg_base_type) :: &
                                                  rel_perm_og_owg_base_type
    PetscReal :: Soco !oil connate saturation
    PetscReal :: Sogcr !critical oil-water saturation    
  contains
    procedure, public :: Init => RPFOGOWGBaseInit
    procedure, public :: Verify => RPFOGOWGBaseVerify
    procedure, public :: Test => RPFOGOWGBaseTest
    procedure, public :: SetupPolynomials => RPFOGOWGBaseSetupPoly
    procedure, public :: RelativePermeability => RPFOGOWGBaseRelPerm !defines argument template
    procedure, public :: GetCriticalSaturation => GetCriticalSatOGOWGBase
  end type rel_perm_og_owg_base_type

  type,public,extends(rel_perm_og_owg_base_type) :: rel_perm_og_owg_MBC_type
    PetscReal :: Sgcr !gas residual saturation
    PetscReal :: Swco !water connate saturation
  contains
    procedure, public :: Init => RPF_OG_OWG_MBC_Init
    procedure, public :: Verify => RPF_OG_OWG_MBC_Verify
    procedure, public :: SetupPolynomials => RPF_OG_OWG_MBC_SetupPoly
    procedure, public :: RelativePermeability => RPF_OG_OWG_MBC_RelPerm
    procedure, public :: RPF_og_owg_MBC_SetSwcoSgcr
  end type rel_perm_og_owg_MBC_type

  type,public,extends(rel_perm_og_owg_base_type) :: rel_perm_og_owg_table_type
  contains
    procedure, public :: Init => RPF_OG_OWG_table_Init
  !   procedure, public :: Verify => RPF_OG_OWG_table_Verify
     procedure, public :: RelativePermeability => RPF_OG_OWG_table_RelPerm
  end type rel_perm_og_owg_table_type

  !-----------------------------------------------------------------------------
  !-- Oil OWG Relative Permeability Functions ------------------------------
  !-----------------------------------------------------------------------------
  ! currently only one option: ECLIPSE model (shoud add here Stone, etc)
  type,abstract,public,extends(rel_perm_owg_base_type) :: &
                                                  rel_perm_oil_owg_base_type
    PetscReal :: Soco !oil connate saturation
    PetscReal :: Socr !oil residual saturation
    class(rel_perm_ow_owg_base_type), pointer :: rel_perm_ow
    class(rel_perm_og_owg_base_type), pointer :: rel_perm_og    
  contains
    procedure, public :: Init => RPFOilOWGBaseInit
    procedure, public :: Verify => RPFOilOWGBaseVerify
    procedure, public :: Test => RPFOilOWGBaseTest
    !procedure, public :: SetupPolynomials => RPFOilOWGBaseSetupPoly
    procedure, public :: RelPermOW
    procedure, public :: RelPermOG
    procedure, public :: GetSowcr => RPFOilOWGBaseGetSowcr
    procedure, public :: GetSogcr => RPFOilOWGBaseGetSogcr
    procedure, public :: RelativePermeability => RPFOilOWGBaseRelPerm !defines argument template
    procedure, public :: GetCriticalSaturation => GetCriticalSatOilOWGBase
    procedure, public :: GetConnateSaturation => GetConnateSatOilOWGBase
    procedure, public :: Strip => RPFOilOWGBaseStrip
  end type rel_perm_oil_owg_base_type

  type,public,extends(rel_perm_oil_owg_base_type) :: rel_perm_oil_owg_ecl_type
    PetscReal :: Swco !oil residual saturation      
    !class(rel_perm_ow_owg_base_type), public, pointer :: rel_perm_ow
    !class(rel_perm_og_owg_base_type), public, pointer :: rel_perm_og    
  contains
    procedure, public :: Init => RPF_oil_ecl_Init
    procedure, public :: Verify => RPF_oil_ecl_Verify
    !procedure, public :: Test => RPF_oil_ecl_Test
    procedure, public :: RelativePermeability => RPF_oil_ecl_RelPerm !defines argument template
    procedure, public :: RPF_oil_ecl_SetSwco
  end type rel_perm_oil_owg_ecl_type

  
  public :: SaturationFunctionOWGRead, &
            PermeabilityFunctionOWGRead, &
            SF_XW_VG_Create, &
            SF_XW_BC_Create, &
            SF_XW_constant_Create, &
            SF_XW_table_Create, &
            SF_OG_VG_SL_Create, &
            SF_OG_constant_Create, &
            SF_OG_table_Create, &
            RPF_ow_owg_linear_Create, &
            RPF_wat_owg_MBC_Create, &
            RPF_gas_owg_MBC_Create, &
            RPF_og_owg_MBC_Create, &
            RPF_ow_owg_MBC_Create, &
            RPF_oil_ecl_Create, &
            RPF_wat_owg_Mualem_VG_Create, &
            RPF_gas_owg_Mualem_VG_Create, &
            RPF_gas_owg_TOUGH2_IRP7_Create, &
            RPF_wat_owg_Burdine_VG_Create, &
            RPF_gas_owg_Burdine_VG_Create, &
            RPF_wat_owg_Burdine_BC_Create, &
            RPF_gas_owg_Burdine_BC_Create, &
            RPF_wat_owg_table_Create, &
            RPF_gas_owg_table_Create, &
            RPF_ow_owg_table_Create, &
            RPF_og_owg_table_Create, &
            SaturationFunctionXWDestroy, &
            SaturationFunctionOGDestroy, &
            WatPermFunctionOWGDestroy, &
            GasPermFunctionOWGDestroy, &
            OGPermFunctionOWGDestroy, &
            OWPermFunctionOWGDestroy, &
            OilPermFunctionOWGDestroy, &
            SetCCOWGPhaseFlags
  
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


  input%ierr = 0
  smooth = PETSC_FALSE
  error_string = 'CHARACTERISTIC_CURVES,SATURATION_FUNCTION_OWG,'
  select type(sf_owg => sat_func_owg)
    class is(sat_func_xw_VG_type)
      error_string = trim(error_string) // 'VAN_GENUCHTEN_XW'
    class is(sat_func_xw_constant_type)
      error_string = trim(error_string) // 'CONSTANT_PRESSURE_XW'
    class is(sat_func_og_constant_type)
      error_string = trim(error_string) // 'CONSTANT_PRESSURE_OG'
  end select

  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    ! read sat_func_owg_base
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
      case('TABLE_NAME')
        !read table name
      case default
        found = PETSC_FALSE
    end select

    if (found) cycle

    !read sat_func_xw_base/sat_func_og_base
    found = PETSC_TRUE
    select type(sf_owg => sat_func_owg)
      class is(sat_func_xw_base_type)
        select case(keyword)
          case('WATER_CONNATE_SATURATION')
            call InputReadDouble(input,option,sf_owg%Swco)
            call InputErrorMsg(input,option,'CONNATE_WATER_SATURATION', &
                               error_string)
          case default
            found = PETSC_FALSE
        end select  
      class is(sat_func_og_base_type)  
        select case(keyword)
          case('GAS_CONNATE_SATURATION')
            call InputReadDouble(input,option,sf_owg%Sgco)
            call InputErrorMsg(input,option,'CONNATE_WATER_SATURATION', &
                               error_string)
          case default
            found = PETSC_FALSE
        end select
    end select    

    if (found) cycle

    select type(sf_owg => sat_func_owg)
      class is(sat_func_xw_VG_type)
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
      class is(sat_func_xw_BC_type)
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
      class is(sat_func_xw_constant_type)
        select case(keyword)                
          case('CONSTANT_PRESSURE')
            call InputReadDouble(input,option, &
                                  sf_owg%constant_capillary_pressure)
            call InputErrorMsg(input,option,'CONSTANT_CAPILLARY_PRESSURE', &
                                error_string)
            sf_owg%pcmax = sf_owg%constant_capillary_pressure
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
      class is(sat_func_og_constant_type)
        select case(keyword)
          case('CONSTANT_PRESSURE')
            call InputReadDouble(input,option, &
                                  sf_owg%constant_capillary_pressure)
            call InputErrorMsg(input,option,'CONSTANT_CAPILLARY_PRESSURE', &
                                error_string)
            sf_owg%pcmax = sf_owg%constant_capillary_pressure
        end select  
    end select
  !add reading instructions for other OWG saturation functions (tables etc)
  end do

  if ( smooth .and. associated(sat_func_owg%sat_func_sl) ) then
    call sat_func_owg%sat_func_sl%SetupPolynomials(option,error_string)
  end if

end subroutine SaturationFunctionOWGRead

! ************************************************************************** !

recursive subroutine PermeabilityFunctionOWGRead(permeability_function, &
                                                 input,option)
  ! 
  ! Reads in contents of a PERMEABILITY_FUNCTION_OWG block
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  
  implicit none
  
  class(rel_perm_owg_base_type) :: permeability_function
  type(input_type), pointer :: input
  type(option_type) :: option
    
  character(len=MAXWORDLENGTH) :: keyword
  
  character(len=MAXWORDLENGTH) :: perm_func_ch_type
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: smooth
  
  
  input%ierr = 0
  smooth = PETSC_FALSE
  error_string = 'CHARACTERISTIC_CURVES,'
  
  select type(rpf => permeability_function)
    class is(rel_perm_wat_owg_base_type)
      error_string = trim(error_string) // 'PERMEABILITY_FUNCTION_WAT,KRW,'
    class is(rel_perm_gas_owg_base_type)
      error_string = trim(error_string) // 'PERMEABILITY_FUNCTION_GAS,KRG,'
    class is(rel_perm_ow_owg_base_type)
      if (rpf%So_is_Sh) then
        error_string = trim(error_string) // &
                                     'PERMEABILITY_FUNCTION_HYDROCARBON,KRH,'
      else  
        error_string = trim(error_string) // 'PERMEABILITY_FUNCTION_OW,KROW,'
      end if  
    class is(rel_perm_og_owg_base_type)
      error_string = trim(error_string) // 'PERMEABILITY_FUNCTION_OG,KROG,'
    class is(rel_perm_oil_owg_base_type)
      error_string = trim(error_string) // 'PERMEABILITY_FUNCTION_OIL,KRO,'
  end select    
      
  select type(rpf => permeability_function)
    class is(RPF_wat_owg_MBC_type)
      error_string = trim(error_string) // 'MOD_BROOKS_COREY'
    class is(RPF_gas_owg_MBC_type)
      error_string = trim(error_string) // 'MOD_BROOKS_COREY'
    class is(rel_perm_ow_owg_MBC_type)
      error_string = trim(error_string) // 'MOD_BROOKS_COREY'
    class is(rel_perm_ow_owg_linear_type)
      error_string = trim(error_string) // 'TOUGH2_LINEAR'
    class is(rel_perm_og_owg_MBC_type)
      error_string = trim(error_string) // 'MOD_BROOKS_COREY'            
    class is(RPF_wat_owg_Mualem_VG_type)
      error_string = trim(error_string) // 'MUALEM_VG'
    class is(RPF_gas_owg_Mualem_VG_type)
      error_string = trim(error_string) // 'MUALEM_VG_SL'
    class is(RPF_wat_owg_Burdine_VG_type)
      error_string = trim(error_string) // 'BURDINE_VG'
    class is(RPF_gas_owg_Burdine_VG_type)
      error_string = trim(error_string) // 'BURDINE_VG_SL'
    class is(RPF_wat_owg_Burdine_BC_type)
      error_string = trim(error_string) // 'BURDINE_BC'
    class is(RPF_gas_owg_Burdine_BC_type)
      error_string = trim(error_string) // 'BURDINE_BC_SL'
    class is(RPF_gas_owg_TOUGH2_IRP7_type)
      error_string = trim(error_string) // 'TOUGH2_IRP7'
    class is(rel_perm_oil_owg_ecl_type)
      error_string = trim(error_string) // 'ECLIPSE'  
  end select
  
  do 
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    !rel perm owg base
    found = PETSC_TRUE
    select case(keyword)
      case('M')
        call InputReadDouble(input,option,permeability_function%m)
        call InputErrorMsg(input,option,'M',error_string)    
      case('MAX_RELATIVE_PERMEABILITY','MAX_REL_PERM')
        call InputReadDouble(input,option,permeability_function%kr_max)
        call InputErrorMsg(input,option,'MAX_RELATIVE_PERMEABILITY', &
                           error_string)
      case('TABLE_NAME')
        !read table name
      case('SMOOTH')
        smooth = PETSC_TRUE
      case default
        found = PETSC_FALSE
    end select
    if (found) cycle

    ! wat base
    found = PETSC_FALSE
    select type(rpf => permeability_function)
      class is(rel_perm_wat_owg_base_type)
        select case(keyword)
        case('WATER_CONNATE_SATURATION')
          call InputReadDouble(input,option,rpf%Swco)
          call InputErrorMsg(input,option,'WATER_CONNATE_SATURATION', &
                             error_string)
          found = PETSC_TRUE          
          case('WATER_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Swcr)
            call InputErrorMsg(input,option,'WATER_RESIDUAL_SATURATION', &
                               error_string)
            found = PETSC_TRUE
          !case default
          !  found = PETSC_FALSE
        end select                       
    end select
    if (found) cycle


    ! wat Burding_BC function
    found = PETSC_FALSE
    select type(rpf =>permeability_function)
      class is(RPF_wat_owg_Burdine_BC_type)
        select case(keyword)
          case('LAMBDA')
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
            found = PETSC_TRUE
        end select    
    end select
    if (found) cycle
    
    ! gas base
    found = PETSC_FALSE
    select type(rpf => permeability_function)
      class is(rel_perm_gas_owg_base_type)
        select case(keyword)
        case('GAS_CONNATE_SATURATION')
            call InputReadDouble(input,option,rpf%Sgco)
            call InputErrorMsg(input,option,'GAS_CONNATE_SATURATION', &
                             error_string)
            found = PETSC_TRUE          
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Sgcr)
            call InputErrorMsg(input,option,'GAS_RESIDUAL_SATURATION', &
                               error_string)
            found = PETSC_TRUE
        end select  
    end select
    if (found) cycle


    ! gas sl function
    found = PETSC_FALSE
    select type(rpf => permeability_function)
      class is(RPF_gas_owg_func_sl_type)
        select case(keyword)
          case('LIQUID_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Slcr)
            call InputErrorMsg(input,option,'LIQUID_RESIDUAL_SATURATION', &
                               error_string)
            found = PETSC_TRUE
        end select    
    end select
    if (found) cycle    

    ! gas Burding_BC function
    found = PETSC_FALSE
    select type(rpf =>permeability_function)
      class is(RPF_gas_owg_Burdine_BC_type)
        select case(keyword)
          case('LAMBDA')
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'LAMBDA',error_string)
            found = PETSC_TRUE
        end select    
    end select
    if (found) cycle
    
    ! ow base
    found = PETSC_FALSE
    select type(rpf => permeability_function)
      class is(rel_perm_ow_owg_base_type)
        select case(keyword)
          case('OIL_RESIDUAL_SATURATION')
            if (rpf%So_is_Sh) then
              option%io_buffer = 'RELATIVE_PERMEABILITY_HYDROCARVON &
                             &requires an hydrocarbon residual sarturation, &
                             &instead of an oil residual saturation: &
                             & Shcr = Socr + Sgcr + Sscr must be entered' 
              call printErrMsg(option)              
            end if  
            call InputReadDouble(input,option,rpf%Sowcr)
            call InputErrorMsg(input,option,'OIL_RESIDUAL_SATURATION', &
                               error_string)
            found = PETSC_TRUE
          case('HYDROCARBON_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Sowcr)
            call InputErrorMsg(input,option,'HYDROCARBON_RESIDUAL_SATURATION',&
                               error_string)
            found = PETSC_TRUE            
        end select    
    end select
    if (found) cycle

    ! og base
    found = PETSC_FALSE
    select type(rpf => permeability_function)
      class is(rel_perm_og_owg_base_type)
        select case(keyword)
          case('OIL_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Sogcr)
            call InputErrorMsg(input,option,'OIL_RESIDUAL_SATURATION', &
                               error_string)
            found = PETSC_TRUE
          !case default
          !  found = PETSC_FALSE
        end select    
    end select
    if (found) cycle

    select type(rpf => permeability_function)    
      class is(rel_perm_oil_owg_ecl_type)
        select case(keyword)
          case('PERMEABILITY_FUNCTION_OW','KROW')
            call InputReadWord(input,option,perm_func_ch_type,PETSC_TRUE)
            call InputErrorMsg(input,option,'perm_func_ch_type',error_string)
            perm_func_ch_type = trim(perm_func_ch_type)
            call StringToUpper(perm_func_ch_type)
            select case(perm_func_ch_type)
              case('MOD_BROOKS_COREY')
                rpf%rel_perm_ow => RPF_ow_owg_MBC_Create()
              case default
                call InputKeywordUnrecognized(perm_func_ch_type, &
                                            'PERMEABILITY_FUNCTION_OW',option)
            end select
            call PermeabilityFunctionOWGRead(rpf%rel_perm_ow,input,option)
          case('KROW_TABLE')
            rpf%rel_perm_ow => RPF_ow_owg_table_Create()
            call InputReadWord(input,option, &
                               rpf%rel_perm_ow%table_name,PETSC_TRUE)
            call InputErrorMsg(input,option,'KROW_TABLE',error_string)
          case('PERMEABILITY_FUNCTION_OG','KROG')
            call InputReadWord(input,option,perm_func_ch_type,PETSC_TRUE)
            call InputErrorMsg(input,option,'perm_func_ch_type',error_string)
            perm_func_ch_type = trim(perm_func_ch_type)
            call StringToUpper(perm_func_ch_type)
            select case(perm_func_ch_type)
              case('MOD_BROOKS_COREY')
                rpf%rel_perm_og => RPF_og_owg_MBC_Create()
              case default
                call InputKeywordUnrecognized(perm_func_ch_type, &
                                            'PERMEABILITY_FUNCTION_OG',option)
            end select
            call PermeabilityFunctionOWGRead(rpf%rel_perm_og,input,option)
          case('KROG_TABLE')
            rpf%rel_perm_og => RPF_og_owg_table_Create()
            call InputReadWord(input,option, &
                               rpf%rel_perm_og%table_name,PETSC_TRUE)
            call InputErrorMsg(input,option,'KROG_TABLE',error_string)
          case('KRO_TABLE')
            rpf%rel_perm_ow => RPF_ow_owg_table_Create()
            rpf%rel_perm_og => RPF_og_owg_table_Create()
            call InputReadWord(input,option, &
                               rpf%rel_perm_ow%table_name,PETSC_TRUE)
            call InputErrorMsg(input,option,'KRO_TABLE',error_string)
            rpf%rel_perm_og%table_name = rpf%rel_perm_ow%table_name
          case default
            call InputKeywordUnrecognized(keyword, &
                      'ECLIPSE relative permeability function',option)
          end select
    end select
  end do !end loop line within rel perm funct def

  !pass parameter to sl_functions
  select type(rpf => permeability_function)
    class is(RPF_wat_owg_func_sl_type)
      select type(rpf_sl => rpf%rel_perm_func_sl)
        class is(RPF_Mualem_VG_liq_type)
          rpf_sl%Sr = rpf%Swcr
          rpf_sl%m = rpf%m
        class is(RPF_Burdine_VG_liq_type)
          rpf_sl%Sr = rpf%Swcr
          rpf_sl%m = rpf%m
      end select  
    class is(RPF_gas_owg_func_sl_type)
      select type(rpf_sl => rpf%rel_perm_func_sl)
        class is(RPF_Mualem_VG_gas_type)
          rpf_sl%Sr = rpf%Slcr
          rpf_sl%Srg = rpf%Sgcr
          rpf_sl%m = rpf%m
        class is(RPF_Burdine_VG_gas_type)
          rpf_sl%Sr = rpf%Slcr
          rpf_sl%Srg = rpf%Sgcr
          rpf_sl%m = rpf%m
        class is(rpf_TOUGH2_IRP7_gas_type)
          rpf_sl%Sr = rpf%Slcr
          rpf_sl%Srg = rpf%Sgcr
          rpf_sl%m = rpf%m          
      end select
    class is(RPF_wat_owg_Burdine_BC_type)
      select type(rpf_sl => rpf%rel_perm_func_sl)
        class is(RPF_Burdine_BC_liq_type)
          rpf_sl%Sr = rpf%Swcr
          rpf_sl%lambda = rpf%lambda
      end select    
    class is(RPF_gas_owg_Burdine_BC_type)
      select type(rpf_sl => rpf%rel_perm_func_sl)
        class is(RPF_Burdine_BC_gas_type)
          rpf_sl%Sr = rpf%Slcr
          rpf_sl%Srg = rpf%Sgcr
          rpf_sl%lambda = rpf%lambda
      end select
  end select  

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
! *********** OWG Saturaton functions base class memebers  ***************** !
! ************************************************************************** !

subroutine SFOWGBaseInit(this)

  implicit none
  
  class(sat_func_owg_base_type) :: this

  this%pcmax = DEFAULT_PCMAX
  this%pcmin = UNINITIALIZED_DOUBLE
  this%analytical_derivative_available = PETSC_FALSE
  this%sat_func_of_pc_available = PETSC_FALSE
  
  this%table_name =''
  nullify(this%table)
  nullify(this%sat_func_sl)

end subroutine SFOWGBaseInit


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


subroutine SFOWGBaseTest(this,cc_name,option)

  use Option_module
  use Material_Aux_class

  implicit none

  class(sat_func_owg_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option

  option%io_buffer = 'SFOWGBaseTest must be extended.'
  call printErrMsg(option)

end subroutine SFOWGBaseTest


subroutine SFOWGBaseSetupPolynomials(this,option,error_string)

  use Option_module

  implicit none

  class(sat_func_owg_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  option%io_buffer = 'SF OWG Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)

end subroutine SFOWGBaseSetupPolynomials


subroutine SFOWGBaseProcessTable(this,char_curves_tables,error_string,option)

  use Option_module

  implicit none

  class(sat_func_owg_base_type) :: this
  class(char_curves_table_type), pointer :: char_curves_tables
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: error_string_lc
  PetscInt :: sat_i_max

  select type(sf => this)
    class is(sat_func_xw_base_type)
      error_string_lc = trim(error_string) // 'PCOW/PCWG,'
    class is(sat_func_og_base_type)
      error_string_lc = trim(error_string) // 'PCOG,'
  end select

  this%table =>  CharCurveTableGetPtrFromList(this%table_name, &
                           char_curves_tables,error_string_lc,option)

  !load cc curve end points from table
 select type(sf => this)
   class is(sat_func_xw_base_type)
     call this%table%CheckCCTVariableExists(CCT_PCXW,error_string_lc,option)
     sf%Swco = sf%table%Swco
     sf%Swcr = sf%table%Swcr
     sf%pcmax = this%table%lookup_table%var_array(CCT_PCXW)%ptr%data(1)
     sat_i_max = this%table%lookup_table%dims(1)
     sf%pcmin = this%table%lookup_table%var_array(CCT_PCXW)%ptr%data(sat_i_max)
   class is(sat_func_og_base_type)
     call this%table%CheckCCTVariableExists(CCT_PCOG,error_string_lc,option)
     sf%Sgco = sf%table%Sgco
     sf%Sgcr = sf%table%Sgcr
     sat_i_max = this%table%lookup_table%dims(1)
     sf%pcmax = this%table%lookup_table%var_array(CCT_PCOG)%ptr%data(sat_i_max)
     sf%pcmin = this%table%lookup_table%var_array(CCT_PCOG)%ptr%data(1)
end select

end subroutine SFOWGBaseProcessTable


function GetPcMax(this)

  implicit none

  PetscReal :: GetPcMax
  class(sat_func_owg_base_type) :: this

  GetPcMax = this%pcmax

end function GetPcMax


function GetPcMin(this)

  implicit none

  PetscReal :: GetPcMin
  class(sat_func_owg_base_type) :: this

  GetPcMin = this%pcmin

end function GetPcMin

! ************************************************************************** !
! *********** XW Saturaton functions  ************************************** !
! ************************************************************************** !

subroutine SFXWBaseInit(this)

  implicit none
  

  class(sat_func_xw_base_type) :: this
  
  call SFOWGBaseInit(this)
  
  this%Swco = UNINITIALIZED_DOUBLE
  this%Swcr = UNINITIALIZED_DOUBLE

end subroutine SFXWBaseInit

! ************************************************************************** !
subroutine SFXWBaseVerify(this,name,option)

  use Option_module

  implicit none
  
  class(sat_func_xw_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  call SFOWGBaseVerify(this,name,option)

  if (Uninitialized(this%Swco)) then
    option%io_buffer = UninitializedMessage('WATER_CONNATE_SATURATION',name)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%Swcr)) then
    option%io_buffer = UninitializedMessage('WATER_RESIDUAL_SATURATION',name)
    call printErrMsg(option)
  endif

end subroutine SFXWBaseVerify

! ************************************************************************** !

subroutine SFXWBaseTest(this,cc_name,option)

  use Option_module
  use Material_Aux_class

  implicit none

  !class(sat_func_owg_base_type) :: this
  class(sat_func_xw_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: num_values = 101
  PetscReal :: capillary_pressure(num_values)
  PetscReal :: wat_saturation(num_values)
  PetscReal :: dpc_dsatw(num_values)
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

  do i = 1, num_values
   call this%CapillaryPressure(wat_saturation(i), capillary_pressure(i), &
                               dpc_dsatw(i),option)
    ! calculate numerical derivatives?
  enddo

  string = trim(cc_name) // '_pcxw_sw.dat'

  open(unit=86,file=string)
  write(86,*) ' "wat_saturation", "capillary pressure", "dpc/dsatw" '
  do i = 1, num_values
    write(86,'(3es14.6)') wat_saturation(i),capillary_pressure(i),dpc_dsatw(i)
  enddo
  close(86)

end subroutine SFXWBaseTest


! ************************************************************************** !

subroutine SFXWBaseCapillaryPressure(this,wat_saturation,capillary_pressure,&
                                       dpc_dsatw,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_xw_base_type) :: this
  PetscReal, intent(in) :: wat_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatw
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  option%io_buffer = 'SFXWBaseCapillaryPressure must be extended.'
  call printErrMsg(option)

end subroutine SFXWBaseCapillaryPressure


subroutine SFXWBaseSaturation(this,capillary_pressure,wat_saturation, &
                                      dsatw_dpc,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_xw_base_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: wat_saturation
  PetscReal, intent(out) :: dsatw_dpc
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  option%io_buffer = 'SFXWBaseSaturation must be extended.'
  call printErrMsg(option)

end subroutine SFXWBaseSaturation


subroutine SFXWSetConnateSatBase(this,swco,option)

use Option_module

implicit none

class(sat_func_xw_base_type) :: this
PetscReal, intent(in) :: swco
type(option_type), intent(inout) :: option

this%Swco = swco

end subroutine SFXWSetConnateSatBase


function SFXWGetConnateSatBase(this,option)

use Option_module

implicit none

PetscReal :: SFXWGetConnateSatBase
class(sat_func_xw_base_type) :: this
type(option_type), intent(inout) :: option

SFXWGetConnateSatBase = this%Swco

end function SFXWGetConnateSatBase


subroutine SFXWSetCriticalSatBase(this,swcr,option)

use Option_module

implicit none

class(sat_func_xw_base_type) :: this
PetscReal, intent(in) :: swcr
type(option_type), intent(inout) :: option

this%Swcr = swcr

end subroutine SFXWSetCriticalSatBase


function SFXWGetCriticalSatBase(this,option)

use Option_module

implicit none

PetscReal :: SFXWGetCriticalSatBase
class(sat_func_xw_base_type) :: this
type(option_type), intent(inout) :: option

SFXWGetCriticalSatBase = this%Swcr

end function SFXWGetCriticalSatBase


subroutine SFXWComputePcMin(this,option)
  
  use Option_module

  implicit none

  class(sat_func_xw_base_type) :: this
  type(option_type), intent(inout) :: option
  
  PetscReal :: sw_max, dpc_dsatw

  sw_max = 1.0d0
  
  select type(sf => this)
    class is(sat_func_xw_table_type)
      !do not do anything - value loaded from table
    class default
       call sf%CapillaryPressure(sw_max,this%pcmin,dpc_dsatw,option)
  end select

end subroutine SFXWComputePcMin

! ************************************************************************** !

function SF_XW_VG_Create()

  implicit none

  class(sat_func_xw_VG_type), pointer :: SF_XW_VG_Create

  allocate(SF_XW_VG_Create)

  call SFXWBaseInit(SF_XW_VG_Create)

  SF_XW_VG_Create%sat_func_sl => SF_VG_Create()

  call SF_XW_VG_Create%Init()

end function SF_XW_VG_Create

! ************************************************************************** !

subroutine SF_XW_VG_Init(this)

  implicit none

  class(sat_func_xw_VG_type) :: this

  this%alpha = UNINITIALIZED_DOUBLE
  this%m = UNINITIALIZED_DOUBLE

  this%analytical_derivative_available = &
      this%sat_func_sl%analytical_derivative_available
  this%sat_func_of_pc_available = PETSC_TRUE

end subroutine SF_XW_VG_Init

! ************************************************************************** !

subroutine SF_XW_VG_Verify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_xw_VG_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION_OWG,SF_OW_VG'
  endif

  call SFXWBaseVerify(this,string,option)

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
      & pressure between X and water - saturation function chosen: ' // &
      trim(string)
    call printErrMsg(option)
  end if

  call this%sat_func_sl%verify(name,option)

end subroutine SF_XW_VG_Verify

! ************************************************************************** !

subroutine SF_XW_VG_CapillaryPressure(this,wat_saturation,capillary_pressure, &
                                      dpc_dsatw,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_xw_VG_type) :: this
  PetscReal, intent(in) :: wat_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatw
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal, parameter :: eps_wat=1.0d-5
  PetscReal :: swa
  PetscReal :: dum

  if (wat_saturation < 0.0) then
    swa = eps_wat
  else   
    swa = wat_saturation
  end if

    if (wat_saturation > 1.D0-1.D-8) then
    call this%sat_func_sl%CapillaryPressure(swa,capillary_pressure, &
                                            dum,option)
    call this%sat_func_sl%CapillaryPressure(1.D0-1.D-8,dum, &
                                            dpc_dsatw,option)
  else
    call this%sat_func_sl%CapillaryPressure(swa,capillary_pressure, &
                                            dpc_dsatw,option)
  endif

end subroutine SF_XW_VG_CapillaryPressure

! ************************************************************************** !

subroutine SF_XW_VG_Saturation(this,capillary_pressure,wat_saturation, &
                                      dsatw_dpc,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_xw_VG_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: wat_saturation
  PetscReal, intent(out) :: dsatw_dpc
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call this%sat_func_sl%Saturation(capillary_pressure, &
                                     wat_saturation,dsatw_dpc,option)

end subroutine SF_XW_VG_Saturation

! ************************************************************************** !
! *****************SF x-Water BC function ********************************** !

function SF_XW_BC_Create()

  implicit none

  class(sat_func_xw_BC_type), pointer :: SF_XW_BC_Create

  allocate(SF_XW_BC_Create)

  call SFXWBaseInit(SF_XW_BC_Create)

  SF_XW_BC_Create%sat_func_sl => SF_BC_Create()

  call SF_XW_BC_Create%Init()

end function SF_XW_BC_Create


subroutine SF_XW_BC_Init(this)

  implicit none

  class(sat_func_xw_BC_type) :: this

  this%alpha = UNINITIALIZED_DOUBLE
  this%lambda = UNINITIALIZED_DOUBLE

  this%analytical_derivative_available = &
      this%sat_func_sl%analytical_derivative_available
  this%sat_func_of_pc_available = PETSC_TRUE

end subroutine SF_XW_BC_Init


subroutine SF_XW_BC_Verify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_xw_BC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION_OWG,SF_OW_BC'
  endif

  call SFXWBaseVerify(this,string,option)

  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA = 1/Pcc',string)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call printErrMsg(option)
  endif

  if ( .not.associated(this%sat_func_sl) ) then
    option%io_buffer = 'Not analytical function associated with the capillary &
      & pressure between X and water - saturation function chosen: ' // &
      trim(string)
    call printErrMsg(option)
  end if

  call this%sat_func_sl%verify(name,option)

end subroutine SF_XW_BC_Verify


subroutine SF_XW_BC_CapillaryPressure(this,wat_saturation,capillary_pressure, &
                                      dpc_dsatw,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_xw_BC_type) :: this
  PetscReal, intent(in) :: wat_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatw
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal, parameter :: eps_wat=1.0d-5
  PetscReal :: swa
  PetscReal :: dum

  if (wat_saturation < 0.0) then
    swa = eps_wat
  else   
    swa = wat_saturation
  end if

  call this%sat_func_sl%CapillaryPressure(swa,capillary_pressure, &
                                          dpc_dsatw,option)
  if (wat_saturation > 1.D0-1.D-8) then
    call this%sat_func_sl%CapillaryPressure(swa,capillary_pressure, &
                                            dum,option)
    call this%sat_func_sl%CapillaryPressure(1.D0-1.D-8,dum, &
                                            dpc_dsatw,option)
  else
    call this%sat_func_sl%CapillaryPressure(swa,capillary_pressure, &
                                            dpc_dsatw,option)
  endif

end subroutine SF_XW_BC_CapillaryPressure


subroutine SF_XW_BC_Saturation(this,capillary_pressure,wat_saturation, &
                                      dsatw_dpc,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_xw_BC_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: wat_saturation
  PetscReal, intent(out) :: dsatw_dpc
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call this%sat_func_sl%Saturation(capillary_pressure, &
                                     wat_saturation,dsatw_dpc,option)

end subroutine SF_XW_BC_Saturation

! ************************************************************************** !

function SF_XW_table_Create()

  implicit none

  class(sat_func_xw_table_type), pointer :: SF_XW_table_Create

  allocate(SF_XW_table_Create)

  call SFXWBaseInit(SF_XW_table_Create)

  call SF_XW_table_Create%Init()

end function SF_XW_table_Create


subroutine SF_XW_table_Init(this)

  implicit none

  class(sat_func_xw_table_type) :: this

  this%analytical_derivative_available = PETSC_TRUE

end subroutine SF_XW_table_Init


subroutine SF_XW_table_CapillaryPressure(this,wat_saturation, &
                              capillary_pressure,dpc_dsatw,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_xw_table_type) :: this
  PetscReal, intent(in) :: wat_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatw
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscErrorCode :: ierr

  call this%table%CharCurveTableVarGrad(wat_saturation,CCT_PCXW, &
                              capillary_pressure,dpc_dsatw,ierr,table_idxs)

end subroutine SF_XW_table_CapillaryPressure


subroutine SF_XW_table_Saturation(this,capillary_pressure,wat_saturation, &
                                       dsatw_dpc,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_xw_table_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: wat_saturation
  PetscReal, intent(out) :: dsatw_dpc
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscErrorCode :: ierr

  call this%table%CharCurvePcInvTableVarGrad(capillary_pressure,CCT_SAT_WAT, &
                                     wat_saturation,dsatw_dpc,ierr,table_idxs)

end subroutine SF_XW_table_Saturation

!************************************************************************** !

function SF_XW_constant_Create()

  implicit none

  class(sat_func_xw_constant_type), pointer :: SF_XW_constant_Create

  allocate(SF_XW_constant_Create)

  call SFXWBaseInit(SF_XW_constant_Create)

  call SF_XW_constant_Create%Init()

end function SF_XW_constant_Create

! ************************************************************************** !

subroutine SF_XW_constant_Init(this)

  implicit none

  class(sat_func_xw_constant_type) :: this

  this%constant_capillary_pressure = UNINITIALIZED_DOUBLE

  this%analytical_derivative_available = PETSC_TRUE

end subroutine SF_XW_constant_Init

! ************************************************************************** !

subroutine SF_XW_constant_Verify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_xw_constant_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  !PO check if this make sense
  if (index(name,'SATURATION_FUNCTION_XW') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION_XW,SF_XW_Constant'
  endif

  call SFXWBaseVerify(this,string,option)

  if (Uninitialized(this%constant_capillary_pressure)) then
    option%io_buffer = 'CONSTANT_CAPILLARY_PRESSURE must be specified &
      &for ' // trim(string) // '.'
    call printErrMsg(option)
  endif

end subroutine SF_XW_constant_Verify

! ************************************************************************** !

subroutine SF_XW_const_CapillaryPressure(this,wat_saturation, &
                               capillary_pressure,dpc_dsatw,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_xw_constant_type) :: this
  PetscReal, intent(in) :: wat_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatw
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  capillary_pressure = this%constant_capillary_pressure
  dpc_dsatw = 0.0d0

end subroutine SF_XW_const_CapillaryPressure


subroutine SF_XW_const_Saturation(this,capillary_pressure,wat_saturation, &
                                      dsatw_dpc,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_xw_constant_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: wat_saturation
  PetscReal, intent(out) :: dsatw_dpc
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  option%io_buffer = 'SF_XW_const_Saturation not supported.'
  call printErrMsg(option)

end subroutine SF_XW_const_Saturation

! ************************************************************************** !
! *********** OG Saturaton functions  ************************************** !
! ************************************************************************** !

subroutine SFOGBaseInit(this)

  implicit none
  
  class(sat_func_og_base_type) :: this
  
  call SFOWGBaseInit(this)
    
  this%Sgco = 0.0d0
  this%Sgcr = UNINITIALIZED_DOUBLE

end subroutine SFOGBaseInit


subroutine SFOGBaseVerify(this,name,option)

  use Option_module

  implicit none
  
  class(sat_func_og_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string

  call SFOWGBaseVerify(this,name,option)

  if (index(name,'SATURATION_FUNCTION_OWG') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION_OWG,SF_OG'
  endif  
  
  if (Uninitialized(this%Sgcr)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%Sgco)) then
    option%io_buffer = UninitializedMessage('GAS_CONNATE_SATURATION',string)
    call printErrMsg(option)
  endif

end subroutine SFOGBaseVerify


subroutine SFOGBaseTest(this,cc_name,option)

  use Option_module
  use Material_Aux_class

  implicit none

  class(sat_func_og_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: num_values = 101
  PetscReal :: capillary_pressure(num_values)
  PetscReal :: gas_saturation(num_values)
  PetscReal :: dpc_dsatg(num_values)
  PetscInt :: i

 ! calculate capillary pressure as a function of saturation
  do i = 1, num_values
    gas_saturation(i) = dble(i-1)*0.01d0
    if (gas_saturation(i) < 1.d-7) then
      gas_saturation(i) = 1.d-7
    else if (gas_saturation(i) > (1.d0-1.d-7)) then
      gas_saturation(i) = 1.d0-1.d-7
    endif
  end do

  write(string,*) cc_name

  do i = 1, num_values
   call this%CapillaryPressure(gas_saturation(i),capillary_pressure(i), &
                               dpc_dsatg(i),option)
  enddo

  string = trim(cc_name) // '_pcog_sg.dat'

  open(unit=86,file=string)
  write(86,*) ' "gas_saturation", "capillary pressure", "dpc/dsatg" '
  do i = 1, num_values
    write(86,'(3es14.6)') gas_saturation(i),capillary_pressure(i),dpc_dsatg(i)
  enddo
  close(86)

end subroutine SFOGBaseTest


subroutine SFOGBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing saturation functions

  use Option_module

  implicit none

  class(sat_func_og_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  option%io_buffer = 'SF OG Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)

end subroutine SFOGBaseSetupPolynomials


subroutine SFOGBaseCapillaryPressure(this,gas_saturation,capillary_pressure, &
                                     dpc_dsatg,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_og_base_type) :: this
  PetscReal, intent(in) :: gas_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  option%io_buffer = 'SFOGBaseCapillaryPressure must be extended.'
  call printErrMsg(option)

end subroutine SFOGBaseCapillaryPressure


subroutine SFOGBaseSaturation(this,capillary_pressure,gas_saturation, &
                                       dsatg_dpc,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_og_base_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: gas_saturation
  PetscReal, intent(out) :: dsatg_dpc
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscErrorCode :: ierr

  option%io_buffer = 'SFOGBaseSaturation must be extended.'
  call printErrMsg(option)

end subroutine SFOGBaseSaturation


subroutine SFOGSetConnateSatBase(this,sgco,option)

use Option_module

implicit none

class(sat_func_og_base_type) :: this
PetscReal :: sgco
type(option_type), intent(inout) :: option

this%Sgco = sgco

end subroutine SFOGSetConnateSatBase


function SFOGGetConnateSatBase(this,option)

use Option_module

implicit none

PetscReal :: SFOGGetConnateSatBase
class(sat_func_og_base_type) :: this
type(option_type), intent(inout) :: option

SFOGGetConnateSatBase = this%Sgco

end function SFOGGetConnateSatBase


subroutine SFOGSetCriticalSatBase(this,sgcr,option)

use Option_module

implicit none

class(sat_func_og_base_type) :: this
PetscReal :: sgcr
type(option_type), intent(inout) :: option

this%Sgcr = sgcr

end subroutine SFOGSetCriticalSatBase


function SFOGGetCriticalSatBase(this,option)

use Option_module

implicit none

PetscReal :: SFOGGetCriticalSatBase
class(sat_func_og_base_type) :: this
type(option_type), intent(inout) :: option

SFOGGetCriticalSatBase = this%Sgcr

end function SFOGGetCriticalSatBase


subroutine SFOGComputePcMin(this,option)

  use Option_module

  implicit none
  
  class(sat_func_og_base_type) :: this
  type(option_type), intent(inout) :: option
  
  PetscReal :: sg_min, dpc_dsatg
  sg_min = 0.0d0
  
  select type(sf => this)
    class is(sat_func_og_table_type)
      !do not do anything - loaded in table post processing
    class default
      call this%CapillaryPressure(sg_min,sf%pcmin,dpc_dsatg,option)
  end select  
  
end subroutine SFOGComputePcMin

! ************************************************************************** !

function SF_OG_constant_Create()

  implicit none

  class(sat_func_og_constant_type), pointer :: SF_OG_constant_Create

  allocate(SF_OG_constant_Create)

  call SFOGBaseInit(SF_OG_constant_Create)

  call SF_OG_constant_Create%Init()

end function SF_OG_constant_Create

! ************************************************************************** !

subroutine SF_OG_constant_Init(this)

  implicit none

  class(sat_func_og_constant_type) :: this

  this%constant_capillary_pressure =  UNINITIALIZED_DOUBLE
  this%analytical_derivative_available = PETSC_TRUE

end subroutine SF_OG_constant_Init

! ************************************************************************** !

subroutine SF_OG_constant_Verify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_og_constant_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION_OG') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION_OG,SF_OG_constant'
  endif

  call SFOGBaseVerify(this,string,option)

  if (Uninitialized(this%constant_capillary_pressure)) then
    option%io_buffer = UninitializedMessage('Pcog_const',string)
    call printErrMsg(option)
  endif

end subroutine SF_OG_constant_Verify

! ************************************************************************** !

subroutine SF_OG_Const_CapillaryPressure(this, gas_saturation, &
                                         capillary_pressure,dpc_dsatg, &
                                         option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_og_constant_type) :: this
  PetscReal, intent(in) :: gas_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  capillary_pressure = this%constant_capillary_pressure
  dpc_dsatg = 0.0d0

end subroutine SF_OG_Const_CapillaryPressure

! ************************************************************************** !


function SF_OG_VG_SL_Create()

  ! Creates the van Genutchten capillary pressure function object for use
  ! between the oil and gas phase treating OIl and Water joinly as liquid

  implicit none

  class(sat_func_og_VG_SL_type), pointer :: SF_OG_VG_SL_Create

  allocate(SF_OG_VG_SL_Create)

  call SFOGBaseInit(SF_OG_VG_SL_Create)

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

  !call SFOWGBaseVerify(this,string,option)
  call SFOGBaseVerify(this,string,option)

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

subroutine SF_OG_VG_SL_CapillaryPressure(this, gas_saturation, &
                                         capillary_pressure,dpc_dsatg, &
                                         option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_og_VG_SL_type) :: this
  PetscReal, intent(in) :: gas_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal, parameter :: eps_liq = 1.0d-5
  PetscReal :: liq_saturation, dpc_dsatl

  liq_saturation = 1.0 - gas_saturation

  if ( liq_saturation < 0.0 ) then
    liq_saturation = eps_liq
  end if

  call this%sat_func_sl%CapillaryPressure(liq_saturation,capillary_pressure, &
                                          dpc_dsatl,option)

  dpc_dsatg = - dpc_dsatl

end subroutine SF_OG_VG_SL_CapillaryPressure

! ************************************************************************** !

function SF_OG_table_Create()

  implicit none

  class(sat_func_og_table_type), pointer :: SF_OG_table_Create

  allocate(SF_OG_table_Create)

  call SFOGBaseInit(SF_OG_table_Create)

  call SF_OG_table_Create%Init()

end function SF_OG_table_Create


subroutine SF_OG_table_Init(this)

  implicit none

  class(sat_func_og_table_type) :: this

  this%analytical_derivative_available = PETSC_TRUE

end subroutine SF_OG_table_Init


subroutine SF_OG_table_CapillaryPressure(this,gas_saturation, &
                                         capillary_pressure,dpc_dsatg, &
                                         option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_og_table_type) :: this
  PetscReal, intent(in) :: gas_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscErrorCode :: ierr

  call this%table%CharCurveTableVarGrad(gas_saturation,CCT_PCOG, &
                                capillary_pressure,dpc_dsatg,ierr,table_idxs)

end subroutine SF_OG_table_CapillaryPressure


subroutine SF_OG_table_Saturation(this,capillary_pressure,gas_saturation, &
                                       dsatg_dpc,option,table_idxs)
  use Option_module

  implicit none

  class(sat_func_og_table_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: gas_saturation
  PetscReal, intent(out) :: dsatg_dpc
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscErrorCode :: ierr

  call this%table%CharCurvePcInvTableVarGrad(capillary_pressure,CCT_SAT_GAS, &
                                     gas_saturation,dsatg_dpc,ierr,table_idxs)

end subroutine SF_OG_table_Saturation


! ************************************************************************** !
! *********** END OWG Saturation functions    ****************************** !
! ************************************************************************** !

! ************************************************************************** !
! *********** OWG Relative Permeability functions  ************************* !
! ************************************************************************** !

subroutine RPFOWGBaseInit(this)

  implicit none

  class(rel_perm_owg_base_type) :: this

  this%kr_max = 1.0d0
  this%m = UNINITIALIZED_DOUBLE

  this%analytical_derivative_available = PETSC_FALSE

  this%table_name =''
  nullify(this%table)
  nullify(this%rel_perm_func_sl)
  nullify(this%poly)

end subroutine RPFOWGBaseInit

! ************************************************************************** !

subroutine RPFOWGBaseVerify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_owg_base_type) :: this
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

subroutine RPFOWGBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing OWG relative permeability functions

  use Option_module

  implicit none

  class(rel_perm_owg_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  option%io_buffer = 'RPF Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)

end subroutine RPFOWGBaseSetupPolynomials


subroutine PermFunctionOWGBaseStrip(this)
  !
  ! Destroys an OWG permeability function
  !
  ! Author: Paolo Orsini
  ! Date: 11/18/17
  !

  implicit none
        
  class(rel_perm_owg_base_type) :: this

  call PolynomialDestroy(this%poly)

  call PermeabilityFunctionDestroy(this%rel_perm_func_sl)

  nullify(this%table)


end subroutine PermFunctionOWGBaseStrip


subroutine RPFOWGBaseProcessTable(this,char_curves_tables,error_string,option)

  use Option_module

  implicit none

  class(rel_perm_owg_base_type) :: this
  class(char_curves_table_type), pointer :: char_curves_tables
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: error_string_lc

  select type(rpf => this)
    class is(rel_perm_wat_owg_base_type)
      error_string_lc = trim(error_string) // 'KRW,'
    class is(rel_perm_gas_owg_base_type)
      error_string_lc = trim(error_string) // 'KRG,'
    class is(rel_perm_ow_owg_base_type)
      error_string_lc = trim(error_string) // 'KROW,'
    class is(rel_perm_og_owg_base_type)
      error_string_lc = trim(error_string) // 'KROG,'
  end select

  this%table =>  CharCurveTableGetPtrFromList(this%table_name, &
                           char_curves_tables,error_string_lc,option)
 !load cc curve end points from table
 select type(rpf => this)
   class is(rel_perm_wat_owg_base_type)
     call this%table%CheckCCTVariableExists(CCT_KRW,error_string_lc,option)
     rpf%Swco = rpf%table%Swco
     rpf%Swcr = rpf%table%Swcr
   class is(rel_perm_gas_owg_base_type)
     call this%table%CheckCCTVariableExists(CCT_KRG,error_string_lc,option)
     rpf%Sgco = rpf%table%Sgco
     rpf%Sgcr = rpf%table%Sgcr
   class is(rel_perm_ow_owg_base_type)
     call this%table%CheckCCTVariableExists(CCT_KROW,error_string_lc,option)
     rpf%Sowcr = rpf%table%Sowcr
   class is(rel_perm_og_owg_base_type)
     call this%table%CheckCCTVariableExists(CCT_KROG,error_string_lc,option)
     rpf%Sogcr = rpf%table%Sogcr
end select

end subroutine RPFOWGBaseProcessTable

! ************************************************************************** !
! **************RPF MBC common functions ************************************ !
! ************************************************************************** !
subroutine MBC_SetupPolynomials(poly,m,kr_max,option,error_string)

  ! Sets up polynomials for smoothing Modified BC permeability function

  use Option_module
  use Utility_module

  implicit none

  type(polynomial_type), pointer :: poly
  PetscReal, intent(in) :: m
  PetscReal, intent(in) :: kr_max
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  PetscReal :: b(4)

  PetscReal :: Se_ph_low

  poly => PolynomialCreate()
  ! fill matix with values
  poly%low = 0.99d0  ! just below saturated
  !poly%low = 0.95d0  ! just below saturated
  poly%high = 1.d0   ! saturated
  Se_ph_low = poly%low

  b(1) = kr_max
  b(2) = kr_max * (Se_ph_low ** m)
  b(3) = 0.d0
  b(4) = m * kr_max * Se_ph_low ** (m - 1.0 )

  call CubicPolynomialSetup(poly%high,poly%low,b)

  poly%coefficients(1:4) = b(1:4)

end subroutine MBC_SetupPolynomials


subroutine MBC_RelPerm_dkr_dSe(poly,m,kr_max,effective_sat,rel_perm,&
                                       dkr_Se,option)
 !
 ! Computes the relative permeability and its derivative WRT Se
 ! for the MBC RPF given the Se value
 !
 ! Author: Paolo Orsini (OGS)
 ! Date: 11/18/2017 - 07/24/2018

  use Option_module
  use Utility_module

  implicit none

  type(polynomial_type), pointer :: poly
  PetscReal, intent(in) :: m
  PetscReal, intent(in) :: kr_max
  PetscReal, intent(in) :: effective_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkr_Se
  type(option_type), intent(inout) :: option

  PetscReal :: Se

  Se = effective_sat

  dkr_Se = 0.0d0

  if (Se >= 1.d0) then
    rel_perm = kr_max
    !!!! EXPERIMENTAL - DS
    ! return slope as though at se = 1 for derivative
    if (Se == 1.d0) then
      dkr_Se = kr_max*m
    endif
    !!!! 
    return
  else if (Se <=  0.d0) then
    rel_perm = 0.d0
    return
  endif

  if (associated(poly)) then
    if (Se > poly%low) then
      call CubicPolynomialEvaluate(poly%coefficients,Se,rel_perm,dkr_Se)
      return
    endif
  endif

  rel_perm = kr_max * (Se ** m)

  dkr_Se = kr_max * m * (Se ** (m-1.0))

end subroutine MBC_RelPerm_dkr_dSe

! ************************************************************************** !

subroutine RPFWatOWGBaseInit(this)

  implicit none

  class(rel_perm_wat_owg_base_type) :: this

  call RPFOWGBaseInit(this)

  this%Swco = UNINITIALIZED_DOUBLE
  this%Swcr = UNINITIALIZED_DOUBLE

end subroutine RPFWatOWGBaseInit


subroutine RPFWatOWGBaseVerify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_wat_owg_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option


  call RPFOWGBaseVerify(this,name,option)

  if (Uninitialized(this%Swco)) then
    option%io_buffer = UninitializedMessage('WATER_CONNATE_SATURATION',name)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%Swcr)) then
    option%io_buffer = UninitializedMessage('WATER_RESIDUAL_SATURATION',name)
    call printErrMsg(option)
  endif

  if (this%Swco > this%Swcr ) then
    option%io_buffer = adjustl(trim(name)) // &
                               ',water connate sat > water critical sat'
    call printErrMsg(option)
  end if  

end subroutine RPFWatOWGBaseVerify


subroutine RPFWatOWGBaseTest(this,cc_name,option)

  use Option_module

  implicit none

  class(rel_perm_wat_owg_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, parameter :: num_values = 101
  PetscReal :: wat_saturation(num_values)
  PetscReal :: krw(num_values), dkrow_satw(num_values)


  wat_saturation = 0.0
  krw = 0.0
  dkrow_satw = 0.0

  do i = 1, num_values
    wat_saturation(i) = dble(i-1)*0.01d0
  enddo

  do i = 1, num_values
    call this%RelativePermeability(wat_saturation(i),krw(i), &
                                   dkrow_satw(i),option)
  end do                                   
 
  write(string,*) cc_name
  
  string = trim(cc_name) // '_wat_rel_perm.dat'
  open(unit=86,file=string)
  write(86,*) '"wat_saturation", "wat_rel_perm", "dkrw/dsatw"'  
  do i = 1, size(wat_saturation)
    write(86,'(4es14.6)') wat_saturation(i), krw(i), dkrow_satw(i)
  enddo
  close(86)

end subroutine RPFWatOWGBaseTest


subroutine RPFWatOWGBaseSetupPolynomials(this,option,error_string)

  use Option_module

  implicit none

  class(rel_perm_wat_owg_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  option%io_buffer = 'RPF Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)

end subroutine RPFWatOWGBaseSetupPolynomials


subroutine RPFWatOWGBaseRelPerm(this,wat_sat,rel_perm, &
                                dkrw_satw,option,table_idxs)
  use Option_module

  implicit none

  class(rel_perm_wat_owg_base_type) :: this
  PetscReal, intent(in) :: wat_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkrw_satw
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  option%io_buffer = 'RPFWatOWGBaseRelPerm must be extended.'
  call printErrMsg(option)

end subroutine RPFWatOWGBaseRelPerm


function RPFWatGetConnateSatOWGBase(this,option)

  use Option_module

  implicit none

  PetscReal :: RPFWatGetConnateSatOWGBase
  class(rel_perm_wat_owg_base_type) :: this
  
  type(option_type), intent(inout) :: option

  RPFWatGetConnateSatOWGBase = this%Swco
  
end function RPFWatGetConnateSatOWGBase


function RPFWatGetCriticalSatOWGBase(this,option)

  use Option_module

  implicit none

  PetscReal :: RPFWatGetCriticalSatOWGBase
  class(rel_perm_wat_owg_base_type) :: this
  
  type(option_type), intent(inout) :: option

  RPFWatGetCriticalSatOWGBase = this%Swcr
  
end function RPFWatGetCriticalSatOWGBase

! ************************************************************************** !
! **************RPF wat MBC function *************************************** !
! ************************************************************************** !
function RPF_wat_owg_MBC_Create()

  implicit none

  class(RPF_wat_owg_MBC_type), pointer :: RPF_wat_owg_MBC_Create

  allocate(RPF_wat_owg_MBC_Create)

  call RPF_wat_owg_MBC_Create%Init()

end function RPF_wat_owg_MBC_Create


subroutine RPF_wat_owg_MBC_Init(this)
  
  implicit none

  class(RPF_wat_owg_MBC_type) :: this

  call RPFWatOWGBaseInit(this)

  this%Sowcr = UNINITIALIZED_DOUBLE

  this%analytical_derivative_available = PETSC_TRUE

end subroutine RPF_wat_owg_MBC_Init


subroutine RPF_wat_owg_MBC_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_wat_owg_MBC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_WAT') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_WAT,RPF_WAT_OWG_MBC'
  endif

  call RPFWatOWGBaseVerify(this,string,option)

  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%Sowcr)) then
    option%io_buffer = UninitializedMessage('OIL_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif

end subroutine RPF_wat_owg_MBC_Verify


subroutine RPF_wat_owg_MBC_SetupPoly(this,option,error_string)

  ! Sets up polynomials for smoothing Modified BC permeability function

  use Option_module
  use Utility_module

  implicit none

  class(RPF_wat_owg_MBC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  !call RPF_OWG_MBC_SetupPolynomials(this,option,error_string)
  call MBC_SetupPolynomials(this%poly,this%m,this%kr_max,option,error_string)

end subroutine RPF_wat_owg_MBC_SetupPoly


subroutine RPF_wat_owg_MBC_SetSowcr(this,sowcr,option)

  ! Sets up Sgcr and Socr needed for Wat MBC

  use Option_module

  implicit none

  class(RPF_wat_owg_MBC_type) :: this
  PetscReal, intent(in) :: sowcr  
  type(option_type) :: option

  this%Sowcr = sowcr

end subroutine RPF_wat_owg_MBC_SetSowcr


subroutine RPF_wat_owg_MBC_RelPerm(this,wat_sat,rel_perm, &
                                       dkrw_satw,option,table_idxs)
  use Option_module

  implicit none

  class(RPF_wat_owg_MBC_type) :: this
  PetscReal, intent(in) :: wat_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkrw_satw
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: Se, dSe_dSw, dkr_dSe


  Se = (wat_sat - this%Swcr) / (1.d0 - this%Swcr - this%Sowcr )

  dSe_dSw = 1.0 / (1.d0 - this%Swcr - this%Sowcr )

  ! scaling WRT Kr_max occurs within RPF_OWG_MBC_RelPerm_dkr_dSe
  call MBC_RelPerm_dkr_dSe(this%poly,this%m,this%kr_max,Se,rel_perm, &
                                         dkr_dSe,option)
                                        
  dkrw_satw = dkr_dSe * dSe_dSw

end subroutine RPF_wat_owg_MBC_RelPerm

! ************************************************************************** !
! **************RPF wat func_sl commons function *************************** !
! ************************************************************************** !
subroutine RPF_wat_owg_func_sl_Init(this)

  implicit none

  class(RPF_wat_owg_func_sl_type) :: this
  
  call RPFWatOWGBaseInit(this)

  this%analytical_derivative_available = PETSC_TRUE
  this%m = UNINITIALIZED_DOUBLE

  select type(rpf => this)
    class is(RPF_wat_owg_Burdine_BC_type)
       rpf%lambda = UNINITIALIZED_DOUBLE   
  end select  
     
end subroutine RPF_wat_owg_func_sl_Init


subroutine RPF_wat_owg_func_sl_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_wat_owg_func_sl_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_WAT') > 0) then
    string = name
  else
    !add select case for messages
    string = trim(name) // 'PERMEABILITY_FUNCTION_WAT,RPF_WAT_OWG_func_sl'
  endif

  call RPFWatOWGBaseVerify(this,string,option)

  if (.not.associated(this%rel_perm_func_sl) ) then
    option%io_buffer = trim(string) // ' Sub SL analytical model not defined'
    call printErrMsg(option)
  else
    call this%rel_perm_func_sl%verify(string,option)
  end if

  select type(rpf => this)
    class is(RPF_wat_owg_Burdine_BC_type)
      if (Uninitialized(rpf%lambda)) then
        option%io_buffer = UninitializedMessage('LAMBDA',string)
        call printErrMsg(option)      
      end if
    class default
      if (Uninitialized(this%m)) then
        option%io_buffer = UninitializedMessage('M',string)
        call printErrMsg(option)
      endif
  end select  

end subroutine RPF_wat_owg_func_sl_Verify


subroutine RPF_wat_owg_func_sl_RelPerm(this,wat_sat,rel_perm, &
                                          dkrw_satw,option,table_idxs)
!
! Author: Paolo Orsini (OGS)
! Date: 07/24/2018

 use Option_module
 use Utility_module

 implicit none

 class(RPF_wat_owg_func_sl_type) :: this
 PetscReal, intent(in) :: wat_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkrw_satw
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 call this%rel_perm_func_sl%RelativePermeability(wat_sat,rel_perm, &
                                                        dkrw_satw,option)

 rel_perm = this%kr_max * rel_perm
 dkrw_satw = this%kr_max * dkrw_satw

end subroutine RPF_wat_owg_func_sl_RelPerm

! ************************************************************************** !
! **************RPF wat Mualem VG function ********************************* !
! ************************************************************************** !
function RPF_wat_owg_Mualem_VG_Create()

  implicit none

  class(RPF_wat_owg_Mualem_VG_type), pointer :: RPF_wat_owg_Mualem_VG_Create

  allocate(RPF_wat_owg_Mualem_VG_Create)

  call RPF_wat_owg_Mualem_VG_Create%Init()
  
  RPF_wat_owg_Mualem_VG_Create%rel_perm_func_sl => RPF_Mualem_VG_liq_create()

end function RPF_wat_owg_Mualem_VG_Create

! ************************************************************************** !
! **************RPF wat Burdine VG function ******************************** !
! ************************************************************************** !
function RPF_wat_owg_Burdine_VG_Create()

  implicit none

  class(RPF_wat_owg_Burdine_VG_type), pointer :: RPF_wat_owg_Burdine_VG_Create

  allocate(RPF_wat_owg_Burdine_VG_Create)

  call RPF_wat_owg_Burdine_VG_Create%Init()

  RPF_wat_owg_Burdine_VG_Create%rel_perm_func_sl => RPF_Burdine_VG_liq_create()

end function RPF_wat_owg_Burdine_VG_Create

! ************************************************************************** !
! **************RPF wat Burdine BC function ******************************** !
! ************************************************************************** !
function RPF_wat_owg_Burdine_BC_Create()

  implicit none

  class(RPF_wat_owg_Burdine_BC_type), pointer :: RPF_wat_owg_Burdine_BC_Create

  allocate(RPF_wat_owg_Burdine_BC_Create)

  call RPF_wat_owg_Burdine_BC_Create%Init()

  RPF_wat_owg_Burdine_BC_Create%rel_perm_func_sl => RPF_Burdine_BC_liq_create()

end function RPF_wat_owg_Burdine_BC_Create

! ************************************************************************** !
! **************RPF wat table function ******************************** !
! ************************************************************************** !
function RPF_wat_owg_table_Create()

  implicit none

  class(RPF_wat_owg_table_type), pointer :: RPF_wat_owg_table_Create

  allocate(RPF_wat_owg_table_Create)

  call RPF_wat_owg_table_Create%Init()

end function RPF_wat_owg_table_Create


subroutine RPF_wat_owg_table_Init(this)

  implicit none

  class(RPF_wat_owg_table_type) :: this
  
  call RPFWatOWGBaseInit(this)

  this%analytical_derivative_available = PETSC_TRUE
     
end subroutine RPF_wat_owg_table_Init


subroutine RPF_wat_owg_table_RelPerm(this,wat_sat,rel_perm, &
                                          dkrw_satw,option,table_idxs)
!
! Author: Paolo Orsini (OGS)
! Date: 08/17/2018

 use Option_module
 use Utility_module

 implicit none

 class(RPF_wat_owg_table_type) :: this
 PetscReal, intent(in) :: wat_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkrw_satw
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 PetscErrorCode:: ierr

 call this%table%CharCurveTableVarGrad(wat_sat,CCT_KRW, &
                                      rel_perm,dkrw_satw,ierr,table_idxs)

end subroutine RPF_wat_owg_table_RelPerm

! ************************************************************************** !
! **************RPF gas owg base functions ********************************* !
! ************************************************************************** !

subroutine RPFGasOWGBaseInit(this)

  implicit none

  class(rel_perm_gas_owg_base_type) :: this

  call RPFOWGBaseInit(this)

  this%Sgco = UNINITIALIZED_DOUBLE
  this%Sgcr = UNINITIALIZED_DOUBLE

end subroutine RPFGasOWGBaseInit


subroutine RPFGasOWGBaseVerify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_gas_owg_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option


  call RPFOWGBaseVerify(this,name,option)

  if (Uninitialized(this%Sgco)) then
    option%io_buffer = UninitializedMessage('GAS_CONNATE_SATURATION',name)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%Sgcr)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',name)
    call printErrMsg(option)
  endif

  if (this%Sgco > 0.0 ) then
    option%io_buffer = trim(name) // ', WARNING Sgco > 0 has been set'
    call printMsg(option)
end if

end subroutine RPFGasOWGBaseVerify


subroutine RPFGasOWGBaseTest(this,cc_name,option)

  use Option_module

  implicit none

  class(rel_perm_gas_owg_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, parameter :: num_values = 101
  PetscReal :: gas_saturation(num_values)
  PetscReal :: krg(num_values), dkrg_satg(num_values)

  gas_saturation = 0.0
  krg = 0.0
  dkrg_satg = 0.0

  do i = 1, num_values
    gas_saturation(i) = dble(i-1)*0.01d0
  enddo

  write(string,*) cc_name

  do i = 1, num_values
    call this%RelativePermeability(gas_saturation(i), &
                                   krg(i),dkrg_satg(i),option)
  enddo
  
  !print any other owg rel perms
  string = trim(cc_name) // '_gas_rel_perm.dat'
  open(unit=86,file=string)
  write(86,*) '"gas_saturation", "gas_rel_perm", "dkrg/dsatg" '
  do i = 1, size(gas_saturation)
    write(86,'(4es14.6)') gas_saturation(i), krg(i), dkrg_satg(i)
  enddo
  close(86)  

end subroutine RPFGasOWGBaseTest


subroutine RPFGasOWGBaseSetupPolynomials(this,option,error_string)

  use Option_module

  implicit none

  class(rel_perm_gas_owg_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  option%io_buffer = 'RPF Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)

end subroutine RPFGasOWGBaseSetupPolynomials


subroutine RPFGasOWGBaseRelPerm(this,gas_sat,rel_perm, &
                                dkrg_satg,option,table_idxs)
  use Option_module

  implicit none

  class(rel_perm_gas_owg_base_type) :: this
  PetscReal, intent(in) :: gas_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkrg_satg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  option%io_buffer = 'RPFGasOWGBaseRelPerm must be extended.'
  call printErrMsg(option)

end subroutine RPFGasOWGBaseRelPerm


function RPFGasGetConnateSatOWGBase(this,option)

  use Option_module

  implicit none

  PetscReal :: RPFGasGetConnateSatOWGBase
  class(rel_perm_gas_owg_base_type) :: this
  type(option_type), intent(inout) :: option

  RPFGasGetConnateSatOWGBase = this%Sgco
  
end function RPFGasGetConnateSatOWGBase


function RPFGasGetCriticalSatOWGBase(this,option)

  use Option_module

  implicit none

  PetscReal :: RPFGasGetCriticalSatOWGBase
  class(rel_perm_gas_owg_base_type) :: this
  type(option_type), intent(inout) :: option

  RPFGasGetCriticalSatOWGBase = this%Sgcr
  
end function RPFGasGetCriticalSatOWGBase


! ************************************************************************** !
! *********** RPF gas MBC ************************************************** !
! ************************************************************************** !

function RPF_gas_owg_MBC_Create()

  implicit none

  class(RPF_gas_owg_MBC_type), pointer :: RPF_gas_owg_MBC_Create

  allocate(RPF_gas_owg_MBC_Create)

  call RPF_gas_owg_MBC_Create%Init()

end function RPF_gas_owg_MBC_Create


subroutine RPF_gas_owg_MBC_Init(this)
  
  implicit none

  class(RPF_gas_owg_MBC_type) :: this

  call RPFGasOWGBaseInit(this)

  this%Swco = UNINITIALIZED_DOUBLE
  this%Sogcr = UNINITIALIZED_DOUBLE
 
  this%analytical_derivative_available = PETSC_TRUE

end subroutine RPF_gas_owg_MBC_Init


subroutine RPF_gas_owg_MBC_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_gas_owg_MBC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_GAS') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_GAS,RPF_GAS_OWG_MBC'
  endif

  call RPFGasOWGBaseVerify(this,string,option)

  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call printErrMsg(option)
  endif

  if (Uninitialized(this%Sogcr)) then
    option%io_buffer = UninitializedMessage('OIL_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif
  if (Uninitialized(this%Swco)) then
    option%io_buffer = UninitializedMessage('WAT_CONNATE_SATURATION',string)
    call printErrMsg(option)
  endif

end subroutine RPF_gas_owg_MBC_Verify


subroutine RPF_gas_owg_MBC_SetupPoly(this,option,error_string)

  ! Sets up polynomials for smoothing Modified BC permeability function

  use Option_module
  use Utility_module

  implicit none

  class(RPF_gas_owg_MBC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  call MBC_SetupPolynomials(this%poly,this%m,this%kr_max,option,error_string)

end subroutine RPF_gas_owg_MBC_SetupPoly


subroutine RPF_gas_owg_MBC_SetSwcoSogcr(this,swco,sogcr,option)

  ! Sets up Swcr and Socr needed for Gas MBC

  use Option_module

  implicit none

  class(RPF_gas_owg_MBC_type) :: this
  PetscReal, intent(in) :: swco,sogcr  
  type(option_type) :: option
  
  this%Swco = swco
  this%Sogcr = sogcr
  
end subroutine RPF_gas_owg_MBC_SetSwcoSogcr


subroutine RPF_gas_owg_MBC_RelPerm(this,gas_sat,rel_perm, &
                                       dkrg_satg,option,table_idxs)
  use Option_module

  implicit none

  class(RPF_gas_owg_MBC_type) :: this
  PetscReal, intent(in) :: gas_sat
  PetscReal, intent(out) :: rel_perm
  PetscReal, intent(out) :: dkrg_satg
  type(option_type), intent(inout) :: option
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: Se, dkr_dSe, dSe_dSg

  Se = (gas_sat - this%Sgcr) / (1.d0 - this%Sgcr - this%Sogcr - this%Swco )

  dSe_dSg = 1.0 / (1.d0 - this%Sgcr - this%Sogcr - this%Swco )

  ! scaling WRT Kr_max occurs within _MBC_RelPerm_dkr_dSe
  call MBC_RelPerm_dkr_dSe(this%poly,this%m,this%kr_max,Se,rel_perm, &
                                         dkr_dSe,option)

  dkrg_satg = dkr_dSe * dSe_dSg

end subroutine RPF_gas_owg_MBC_RelPerm

! ************************************************************************** !
! **************RPF gas func_sl commons function *************************** !
! ************************************************************************** !
subroutine RPF_gas_owg_func_sl_Init(this)

  implicit none

  class(RPF_gas_owg_func_sl_type) :: this
  
  call RPFGasOWGBaseInit(this)

  ! base overwrite
  this%analytical_derivative_available = PETSC_TRUE
  this%m = UNINITIALIZED_DOUBLE
  
  this%Slcr= UNINITIALIZED_DOUBLE
  
  select type(rpf => this)
    class is(RPF_gas_owg_Burdine_BC_type)
       rpf%lambda = UNINITIALIZED_DOUBLE   
  end select  
     
end subroutine RPF_gas_owg_func_sl_Init


subroutine RPF_gas_owg_func_sl_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_gas_owg_func_sl_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_GAS') > 0) then
    string = name
  else
    !add select case for messages
    string = trim(name) // 'PERMEABILITY_FUNCTION_GAS,RPF_GAS_OWG_func_sl:'
  endif

  call RPFGasOWGBaseVerify(this,string,option)

  if (.not.associated(this%rel_perm_func_sl) ) then
    option%io_buffer = trim(string) // ' Sub SL analytical model not defined'
    call printErrMsg(option)
  else
    call this%rel_perm_func_sl%verify(string,option)
  end if

  if (Uninitialized(this%Slcr)) then
    option%io_buffer = UninitializedMessage('LIQUID_CRITICAL_SATURATION', &
                                                                       string)
    call printErrMsg(option)          
  end if

  select type(rpf => this)
    class is(RPF_gas_owg_Burdine_BC_type)
      if (Uninitialized(rpf%lambda)) then
        option%io_buffer = UninitializedMessage('LAMBDA',string)
        call printErrMsg(option)      
      end if  
    class default
      if (Uninitialized(this%m)) then
        option%io_buffer = UninitializedMessage('M',string)
        call printErrMsg(option)
      end if
  end select  

end subroutine RPF_gas_owg_func_sl_Verify


subroutine RPF_gas_owg_func_sl_RelPerm(this,gas_sat,rel_perm, &
                                          dkrg_satg,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(RPF_gas_owg_func_sl_type) :: this
 PetscReal, intent(in) :: gas_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkrg_satg
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 PetscReal :: liq_sat, dkrg_satl

 liq_sat = 1.0 - gas_sat
 dkrg_satl = 0.0

 if (gas_sat <= this%Sgcr ) then
   rel_perm = 0.0
   dkrg_satl = 0.0
   return
 end if

 call this%rel_perm_func_sl%RelativePermeability(liq_sat,rel_perm, &
                                                        dkrg_satl,option)

 rel_perm = this%kr_max * rel_perm
 
 !Sg + Sl = 1  => dSl/dSg -1; dkrg/dSg = dkrg/dSl * dSl/dSg = - dkrg/dSl
 !and scaling to the maximum value
 dkrg_satg = - this%kr_max * dkrg_satl

end subroutine RPF_gas_owg_func_sl_RelPerm

! ************************************************************************** !
! *********** RPF gas Mualem VG ******************************************** !
! ************************************************************************** !

function RPF_gas_owg_Mualem_VG_Create()

  implicit none

  class(RPF_gas_owg_Mualem_VG_type), pointer :: RPF_gas_owg_Mualem_VG_Create

  allocate(RPF_gas_owg_Mualem_VG_Create)

  call RPF_gas_owg_Mualem_VG_Create%Init()

  RPF_gas_owg_Mualem_VG_Create%rel_perm_func_sl => RPF_Mualem_VG_gas_create()

end function RPF_gas_owg_Mualem_VG_Create

! ************************************************************************** !
! *********** RPF gas TOUGH2_IRP7  ***************************************** !
! ************************************************************************** !

function RPF_gas_owg_TOUGH2_IRP7_Create()

  implicit none

  class(RPF_gas_owg_TOUGH2_IRP7_type), pointer :: &
                                            RPF_gas_owg_TOUGH2_IRP7_Create

  allocate(RPF_gas_owg_TOUGH2_IRP7_Create)

  call RPF_gas_owg_TOUGH2_IRP7_Create%Init()

  RPF_gas_owg_TOUGH2_IRP7_Create%rel_perm_func_sl => &
                                      RPF_TOUGH2_IRP7_gas_create()

end function RPF_gas_owg_TOUGH2_IRP7_Create

! ************************************************************************** !
! *********** RPF gas Burding VG ******************************************* !
! ************************************************************************** !

function RPF_gas_owg_Burdine_VG_Create()

  implicit none

  class(RPF_gas_owg_Burdine_VG_type), pointer :: RPF_gas_owg_Burdine_VG_Create

  allocate(RPF_gas_owg_Burdine_VG_Create)

  call RPF_gas_owg_Burdine_VG_Create%Init()

  RPF_gas_owg_Burdine_VG_Create%rel_perm_func_sl => RPF_Burdine_VG_gas_create()

end function RPF_gas_owg_Burdine_VG_Create

! ************************************************************************** !
! *********** RPF gas Burding BC ******************************************* !
! ************************************************************************** !

function RPF_gas_owg_Burdine_BC_Create()

  implicit none

  class(RPF_gas_owg_Burdine_BC_type), pointer :: RPF_gas_owg_Burdine_BC_Create

  allocate(RPF_gas_owg_Burdine_BC_Create)

  call RPF_gas_owg_Burdine_BC_Create%Init()

  RPF_gas_owg_Burdine_BC_Create%rel_perm_func_sl => RPF_Burdine_BC_gas_create()

end function RPF_gas_owg_Burdine_BC_Create

! ************************************************************************** !
! *********** RPF gas table ******************************************* !
! ************************************************************************** !

function RPF_gas_owg_table_Create()

  implicit none

  class(RPF_gas_owg_table_type), pointer :: RPF_gas_owg_table_Create

  allocate(RPF_gas_owg_table_Create)

  call RPF_gas_owg_table_Create%Init()

end function RPF_gas_owg_table_Create


subroutine RPF_gas_owg_table_Init(this)

  implicit none

  class(RPF_gas_owg_table_type) :: this
  
  call RPFGasOWGBaseInit(this)

  this%analytical_derivative_available = PETSC_TRUE
     
end subroutine RPF_gas_owg_table_Init


subroutine RPF_gas_owg_table_RelPerm(this,gas_sat,rel_perm, &
                                          dkrg_satg,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(RPF_gas_owg_table_type) :: this
 PetscReal, intent(in) :: gas_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkrg_satg
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 PetscErrorCode:: ierr

 call this%table%CharCurveTableVarGrad(gas_sat,CCT_KRG, &
                                      rel_perm,dkrg_satg,ierr,table_idxs)

end subroutine RPF_gas_owg_table_RelPerm

!-----------------------------------------------------------------------------
!-- RPF Oil-Water base ------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine RPFOWOWGBaseInit(this)

  implicit none

  class(rel_perm_ow_owg_base_type) :: this
  
  call RPFOWGBaseInit(this)

  this%Soco = 0.0 !different than zero only for oil-wetting problems
  this%Sowcr = UNINITIALIZED_DOUBLE
  this%So_is_Sh = PETSC_FALSE
     
end subroutine RPFOWOWGBaseInit


subroutine RPFOWOWGBaseVerify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_ow_owg_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  call RPFOWGBaseVerify(this,name,option)


   if (this%Soco > 0.0) then
     option%io_buffer = 'Oil Connate Sat /=0; Oil-wet system not supported'
     call printErrMsg(option)
   end if

   if (this%So_is_Sh) then
     if (Uninitialized(this%Sowcr)) then
       option%io_buffer = &
             UninitializedMessage('HYDROCARBON_RESIDUAL_SATURATION',name)
       call printErrMsg(option)
     endif
   else
     if (Uninitialized(this%Sowcr)) then
       option%io_buffer = UninitializedMessage('OW_RESIDUAL_SATURATION',name)
       call printErrMsg(option)
     endif
   end if

end subroutine RPFOWOWGBaseVerify


subroutine RPFOWOWGBaseTest(this,cc_name,option)

  use Option_module

  implicit none

  class(rel_perm_ow_owg_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, parameter :: num_values = 101
  PetscReal :: oil_saturation(num_values)
  PetscReal :: krow(num_values), dkrow_sato(num_values)

  oil_saturation = 0.0
  krow = 0.0
  dkrow_sato = 0.0

  do i = 1, num_values
    oil_saturation(i) = dble(i-1)*0.01d0
  enddo

  write(string,*) cc_name

  do i = 1, num_values
    call this%RelativePermeability(oil_saturation(i), &
                                   krow(i),dkrow_sato(i),option)
  enddo

  if (this%So_is_Sh) then
    string = trim(cc_name) // '_hw_rel_perm.dat'
  else  
    string = trim(cc_name) // '_ow_rel_perm.dat'
  end if  
  
  open(unit=86,file=string)
  
  if (this%So_is_Sh) then
    write(86,*) '"Sh", "hw_rel_perm, "dkrhw/dsath" '
  else  
    write(86,*) '"oil_saturation", "ow_rel_perm, "dkrow/dsato" '
  end if  
  
  do i = 1, size(oil_saturation)
    write(86,'(3es14.6)') oil_saturation(i), krow(i), dkrow_sato(i)
  enddo
  close(86)

end subroutine RPFOWOWGBaseTest


subroutine RPFOWOWGBaseSetupPoly(this,option,error_string)

  use Option_module

  implicit none

  class(rel_perm_ow_owg_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  option%io_buffer ='RPF OW Smoothing not supported for' // trim(error_string)
  call printErrMsg(option)

end subroutine RPFOWOWGBaseSetupPoly


subroutine RPFOWOWGBaseRelPerm(this,oil_sat,rel_perm, &
                                          dkr_sato,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_ow_owg_base_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkr_sato
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 option%io_buffer ='RPFOWOWGBaseRelPerm must be extended.'
 call printErrMsg(option)

end subroutine RPFOWOWGBaseRelPerm


function RPFOWGetCriticalSatOWGBase(this,option)

  use Option_module

  implicit none

  PetscReal :: RPFOWGetCriticalSatOWGBase
  class(rel_perm_ow_owg_base_type) :: this
  type(option_type), intent(inout) :: option

  RPFOWGetCriticalSatOWGBase = this%Sowcr
  
end function RPFOWGetCriticalSatOWGBase


function RPFOWGetConnateSaturation(this)
  
  implicit none
  
  PetscReal :: RPFOWGetConnateSaturation
  class(rel_perm_ow_owg_base_type) :: this
  
  RPFOWGetConnateSaturation = this%Soco
  
end function  

!-----------------------------------------------------------------------------
!-- RPF Oil-Water linear -----------------------------------------------------
!-----------------------------------------------------------------------------
function RPF_ow_owg_linear_Create()

  implicit none

  class(rel_perm_ow_owg_linear_type), pointer :: RPF_ow_owg_linear_Create

  allocate(RPF_ow_owg_linear_Create)

  call RPF_ow_owg_linear_Create%Init()

end function RPF_ow_owg_linear_Create


subroutine RPF_ow_owg_linear_Init(this)
 
  implicit none
   
  class(rel_perm_ow_owg_linear_type) :: this
   
  call RPFOWOWGBaseInit(this)
   
  this%analytical_derivative_available = PETSC_TRUE
 
end subroutine RPF_ow_owg_linear_Init


subroutine RPF_ow_owg_linear_Verify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_ow_owg_linear_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OW') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OW,TOUGH2_LINEAR.'
  endif

  call RPFOWOWGBaseVerify(this,string,option)

end subroutine RPF_ow_owg_linear_Verify


subroutine RPF_ow_owg_linear_RelPerm(this,oil_sat,rel_perm, &
                                          dkr_sato,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_ow_owg_linear_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkr_sato
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 PetscReal :: Seo, swa, soa

 !PO done to get same truncation error than previous Kro implementation
 !   in order to pass regression tests
 swa = 1.0d0 - oil_sat
 soa = 1.0d0 - swa
 
 !Seo = (oil_sat - this%Sowcr) / (1.d0 - this%Sowcr)
 Seo = (soa - this%Sowcr) / (1.d0 - this%Sowcr)

 dkr_sato = 1.d0 / (1.d0 - this%Sowcr)
 
 if (Seo >= 1.d0) then
   rel_perm = 1.d0
   dkr_sato = 0.d0
   return
 else if (Seo <=  0.d0) then
   rel_perm = 0.d0
   dkr_sato = 0.d0
   return
 endif

 rel_perm = Seo


end subroutine RPF_ow_owg_linear_RelPerm

!-----------------------------------------------------------------------------
!-- RPF Oil-Water MBC ------------------------------------------------------
!-----------------------------------------------------------------------------
function RPF_ow_owg_MBC_Create()

  implicit none

  class(rel_perm_ow_owg_MBC_type), pointer :: RPF_ow_owg_MBC_Create

  allocate(RPF_ow_owg_MBC_Create)

  call RPF_ow_owg_MBC_Create%Init()

end function RPF_ow_owg_MBC_Create


subroutine RPF_ow_owg_MBC_Init(this)
 
  implicit none
   
  class(rel_perm_ow_owg_MBC_type) :: this
   
  call RPFOWOWGBaseInit(this)
   
  this%Swcr = UNINITIALIZED_DOUBLE
   
  this%analytical_derivative_available = PETSC_TRUE
 
end subroutine RPF_ow_owg_MBC_Init


subroutine RPF_ow_owg_MBC_Verify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_ow_owg_MBC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OW') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OW, MBC'
  endif

  call RPFOWOWGBaseVerify(this,string,option)

  if (Uninitialized(this%Swcr)) then
    option%io_buffer = UninitializedMessage('WATER_CRITICAL_SATURATION',string)
    call printErrMsg(option)          
  end if

  if (Uninitialized(this%m)) then
   option%io_buffer = UninitializedMessage('M',string)
   call printErrMsg(option)          
  end if

end subroutine RPF_ow_owg_MBC_Verify


subroutine RPF_ow_owg_MBC_SetupPoly(this,option,error_string)

  use Option_module

  implicit none

  class(rel_perm_ow_owg_MBC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  call MBC_SetupPolynomials(this%poly,this%m,this%kr_max,option,error_string)

end subroutine RPF_ow_owg_MBC_SetupPoly


subroutine RPF_ow_owg_MBC_RelPerm(this,oil_sat,rel_perm, &
                                          dkr_sato,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_ow_owg_MBC_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkr_sato
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 PetscReal :: Se, dSe_dSo, dkr_dSe

 Se = (oil_sat - this%Sowcr) / (1.d0 - this%Sowcr - this%Swcr )

 dSe_dSo = 1.0 / (1.d0 - this%Sowcr - this%Swcr)
 
 ! scaling WRT Kr_max occurs within MBC_RelPerm_dkr_dSe
 call MBC_RelPerm_dkr_dSe(this%poly,this%m,this%kr_max,Se,rel_perm, &
                                        dkr_dSe,option)
                                       
 dkr_sato = dkr_dSe * dSe_dSo

end subroutine RPF_ow_owg_MBC_RelPerm


subroutine RPF_ow_owg_MBC_SetSwcr(this,swcr,option)

  use Option_module

  implicit none

  class(rel_perm_ow_owg_MBC_type) :: this
  PetscReal :: swcr
  type(option_type), intent(inout) :: option

  this%Swcr = swcr
  
end subroutine RPF_ow_owg_MBC_SetSwcr

! ************************************************************************** !
! *********** RPF oil-water table ****************************************** !
! ************************************************************************** !

function RPF_ow_owg_table_Create()

  implicit none

  class(rel_perm_ow_owg_table_type), pointer :: RPF_ow_owg_table_Create

  allocate(RPF_ow_owg_table_Create)

  call RPF_ow_owg_table_Create%Init()

end function RPF_ow_owg_table_Create


subroutine RPF_ow_owg_table_Init(this)

  implicit none

  class(rel_perm_ow_owg_table_type) :: this
  
  call RPFOWOWGBaseInit(this)

  this%analytical_derivative_available = PETSC_TRUE
     
end subroutine RPF_ow_owg_table_Init


subroutine RPF_ow_owg_table_RelPerm(this,oil_sat,rel_perm, &
                                          dkr_sato,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_ow_owg_table_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkr_sato
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 PetscErrorCode:: ierr

 call this%table%CharCurveTableVarGrad(oil_sat,CCT_KROW, &
                                      rel_perm,dkr_sato,ierr,table_idxs)

end subroutine RPF_ow_owg_table_RelPerm


!-----------------------------------------------------------------------------
!-- RPF Oil-Gas base ---------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine RPFOGOWGBaseInit(this)

  implicit none

  class(rel_perm_og_owg_base_type) :: this
  
  call RPFOWGBaseInit(this)

  this%Soco = 0.0d0 !different than zero only for non-wetting system
  this%Sogcr = UNINITIALIZED_DOUBLE
     
end subroutine RPFOGOWGBaseInit


subroutine RPFOGOWGBaseVerify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_og_owg_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  call RPFOWGBaseVerify(this,name,option)

   if ( this%Soco > 0.0d0 ) then
     option%io_buffer = 'Oil Connate Sat /=0; Oil-wet system not supported'
     call printErrMsg(option)    
   end if

   if (Uninitialized(this%Sogcr)) then
     option%io_buffer = UninitializedMessage('OG_CRITICAL_SATURATION',name)
     call printErrMsg(option)          
   end if

end subroutine RPFOGOWGBaseVerify


subroutine RPFOGOWGBaseTest(this,cc_name,option)

  use Option_module

  implicit none

  class(rel_perm_og_owg_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, parameter :: num_values = 101
  PetscReal :: oil_saturation(num_values)
  PetscReal :: krog(num_values), dkrog_sato(num_values)

  oil_saturation = 0.0
  krog = 0.0
  dkrog_sato = 0.0

  do i = 1, num_values
    oil_saturation(i) = dble(i-1)*0.01d0
  enddo

  write(string,*) cc_name

  do i = 1, num_values
    call this%RelativePermeability(oil_saturation(i), &
                                   krog(i),dkrog_sato(i),option)
  enddo

  string = trim(cc_name) // '_og_rel_perm.dat'
  open(unit=86,file=string)
  write(86,*) '"oil_saturation", "og_rel_perm, "dkrog/dsato" '
  do i = 1, size(oil_saturation)
    write(86,'(3es14.6)') oil_saturation(i), krog(i), dkrog_sato(i)
  enddo
  close(86)

end subroutine RPFOGOWGBaseTest


subroutine RPFOGOWGBaseSetupPoly(this,option,error_string)

  use Option_module

  implicit none

  class(rel_perm_og_owg_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  option%io_buffer ='RPF OG Smoothing not supported for' // trim(error_string)
  call printErrMsg(option)

end subroutine RPFOGOWGBaseSetupPoly


subroutine RPFOGOWGBaseRelPerm(this,oil_sat,rel_perm, &
                                          dkr_sato,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_og_owg_base_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkr_sato
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 option%io_buffer ='RPFOWOWGBaseRelPerm must be extended.'
 call printErrMsg(option)

end subroutine RPFOGOWGBaseRelPerm


function GetCriticalSatOGOWGBase(this,option)

  use Option_module

  implicit none

  PetscReal :: GetCriticalSatOGOWGBase
  class(rel_perm_og_owg_base_type) :: this
  type(option_type), intent(inout) :: option

  GetCriticalSatOGOWGBase = this%Sogcr
  
end function GetCriticalSatOGOWGBase

!-----------------------------------------------------------------------------
!-- RPF Oil-Gas MBC ------------------------------------------------------
!-----------------------------------------------------------------------------
function RPF_og_owg_MBC_Create()

  implicit none

  class(rel_perm_og_owg_MBC_type), pointer :: RPF_og_owg_MBC_Create

  allocate(RPF_og_owg_MBC_Create)

  call RPF_og_owg_MBC_Create%Init()

end function RPF_og_owg_MBC_Create


subroutine RPF_og_owg_MBC_Init(this)
 
  implicit none
   
  class(rel_perm_og_owg_MBC_type) :: this
   
  call RPFOGOWGBaseInit(this)
  
  this%Sgcr = UNINITIALIZED_DOUBLE
  this%Swco = UNINITIALIZED_DOUBLE
   
  this%analytical_derivative_available = PETSC_TRUE
 
end subroutine RPF_og_owg_MBC_Init


subroutine RPF_og_owg_MBC_Verify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_og_owg_MBC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION_OG') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OG, MBC:'
  endif

  call RPFOGOWGBaseVerify(this,string,option)

  if (Uninitialized(this%Sgcr)) then
    option%io_buffer = UninitializedMessage('GAS_CRITICAL_SATURATION',string)
    call printErrMsg(option)          
  end if

  if (Uninitialized(this%Swco)) then
    option%io_buffer = UninitializedMessage('WAT_CONNATE_SATURATION',string)
    call printErrMsg(option)          
  end if

  if (Uninitialized(this%m)) then
   option%io_buffer = UninitializedMessage('M',string)
   call printErrMsg(option)          
  end if

end subroutine RPF_og_owg_MBC_Verify


subroutine RPF_og_owg_MBC_SetupPoly(this,option,error_string)

  use Option_module

  implicit none

  class(rel_perm_og_owg_MBC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  call MBC_SetupPolynomials(this%poly,this%m,this%kr_max,option,error_string)

end subroutine RPF_og_owg_MBC_SetupPoly


subroutine RPF_og_owg_MBC_RelPerm(this,oil_sat,rel_perm, &
                                          dkr_sato,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_og_owg_MBC_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkr_sato
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 PetscReal :: Se, dSe_dSo, dkr_dSe


 Se = (oil_sat - this%Sogcr) / (1.d0 - this%Sogcr - this%Sgcr - this%Swco )

 dSe_dSo = 1.0 / (1.d0 - this%Sogcr - this%Sgcr - this%Swco )
 
 ! scaling WRT Kr_max occurs within MBC_RelPerm_dkr_dSe
 call MBC_RelPerm_dkr_dSe(this%poly,this%m,this%kr_max,Se,rel_perm, &
                                        dkr_dSe,option)
                                       
 dkr_sato = dkr_dSe * dSe_dSo

end subroutine RPF_og_owg_MBC_RelPerm


subroutine RPF_og_owg_MBC_SetSwcoSgcr(this,swco,sgcr,option)

  use Option_module

  implicit none

  class(rel_perm_og_owg_MBC_type) :: this
  PetscReal :: swco,sgcr
  type(option_type), intent(inout) :: option

  this%Swco = swco
  this%Sgcr = sgcr
  
end subroutine RPF_og_owg_MBC_SetSwcoSgcr

! ************************************************************************** !
! *********** RPF oil-gas table ****************************************** !
! ************************************************************************** !

function RPF_og_owg_table_Create()

  implicit none

  class(rel_perm_og_owg_table_type), pointer :: RPF_og_owg_table_Create

  allocate(RPF_og_owg_table_Create)

  call RPF_og_owg_table_Create%Init()

end function RPF_og_owg_table_Create


subroutine RPF_og_owg_table_Init(this)

  implicit none

  class(rel_perm_og_owg_table_type) :: this
  
  call RPFOGOWGBaseInit(this)

  this%analytical_derivative_available = PETSC_TRUE
     
end subroutine RPF_og_owg_table_Init


subroutine RPF_og_owg_table_RelPerm(this,oil_sat,rel_perm, &
                                          dkr_sato,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_og_owg_table_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkr_sato
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 PetscErrorCode:: ierr

 call this%table%CharCurveTableVarGrad(oil_sat,CCT_KROG, &
                                      rel_perm,dkr_sato,ierr,table_idxs)

end subroutine RPF_og_owg_table_RelPerm

!-----------------------------------------------------------------------------
!-- RPF Oil base -------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine RPFOilOWGBaseInit(this)

  implicit none

  class(rel_perm_oil_owg_base_type) :: this
  
  call RPFOWGBaseInit(this)

  this%Soco = 0.0d0
  this%Socr = UNINITIALIZED_DOUBLE
  
  nullify(this%rel_perm_ow)
  nullify(this%rel_perm_og)  
     
end subroutine RPFOilOWGBaseInit


subroutine RPFOilOWGBaseVerify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_oil_owg_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  call RPFOWGBaseVerify(this,name,option)

   if (Uninitialized(this%Socr)) then
     option%io_buffer = UninitializedMessage('OIL_CRITICAL_SATURATION',name)
     call printErrMsg(option)          
   end if

end subroutine RPFOilOWGBaseVerify


subroutine RPFOilOWGBaseTest(this,cc_name,option)

  use Option_module

  implicit none

  class(rel_perm_oil_owg_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i, i_oil, i_gas
  PetscInt, parameter :: num_values = 11
  PetscReal :: oil_saturation(num_values)
  PetscReal :: gas_saturation(num_values)
  PetscReal :: kr(num_values,num_values)
  PetscReal :: dkr_sato(num_values,num_values)
  PetscReal :: dkr_satg(num_values,num_values)  

  gas_saturation = 0.0
  oil_saturation = 0.0
  kr = 0.0
  dkr_sato = 0.0
  dkr_satg = 0.0

  do i = 1, num_values
    oil_saturation(i) = dble(i-1)*0.1d0
  enddo
  gas_saturation = oil_saturation

  !two-dimensional function
  do i_oil = 1, num_values
    do i_gas = 1, num_values
      call this%RelativePermeability(oil_saturation(i_oil), & 
                                     gas_saturation(i_gas), &
                                     kr(i_oil,i_gas),dkr_sato(i_oil,i_gas), &
                                     dkr_satg(i_oil,i_gas),option)
    end do 
  end do

  write(string,*) trim(cc_name) // '_oil_rel_perm.dat'  

  open(unit=86,file=string)
  write(86,*) '"oil_saturation", "gas_saturation", "oil_rel_perm", &
              &"dkro/dsato", " dkro/dsatg"'
              
  do i_oil = 1, size(oil_saturation)
    do i_gas = 1, size(gas_saturation)
        write(86,'(5es14.6)') oil_saturation(i_oil), gas_saturation(i_gas), &
                              kr(i_oil,i_gas),dkr_sato(i_oil,i_gas), &
                              dkr_satg(i_oil,i_gas)
    end do                          
  enddo
  close(86)

  !print the krow and kwog associated with kro_eclipse 
  select type (this)
    class is(rel_perm_oil_owg_ecl_type)
      call this%rel_perm_ow%Test(cc_name,option)
      call this%rel_perm_og%Test(cc_name,option)
  end select

end subroutine RPFOilOWGBaseTest


subroutine RPFOilOWGBaseRelPerm(this,oil_sat,gas_sat,rel_perm, &
                                dkro_sato,dkro_satg,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_oil_owg_base_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(in) :: gas_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkro_sato
 PetscReal, intent(out) :: dkro_satg
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 option%io_buffer ='RPFOilOWGBaseRelPerm must be extended.'
 call printErrMsg(option)

end subroutine RPFOilOWGBaseRelPerm


function GetCriticalSatOilOWGBase(this,option)

  use Option_module

  implicit none

  PetscReal :: GetCriticalSatOilOWGBase
  class(rel_perm_oil_owg_base_type) :: this
  type(option_type), intent(inout) :: option

  GetCriticalSatOilOWGBase = this%Socr
  
end function GetCriticalSatOilOWGBase


function GetConnateSatOilOWGBase(this)

  implicit none

  PetscReal :: GetConnateSatOilOWGBase
  class(rel_perm_oil_owg_base_type) :: this
  
  GetConnateSatOilOWGBase = this%Soco

end function


subroutine RPFOilOWGBaseStrip(this)
  !
  ! Author: Paolo Orsini
  ! Date: 11/18/17
  !

  implicit none
        
  class(rel_perm_oil_owg_base_type) :: this

  call PermFunctionOWGBaseStrip(this)
 
  call OWPermFunctionOWGDestroy(this%rel_perm_ow)
  call OGPermFunctionOWGDestroy(this%rel_perm_og)
  
end subroutine RPFOilOWGBaseStrip

!-----------------------------------------------------------------------------
!-- RPF Oil Eclipse model ----------------------------------------------------
!-----------------------------------------------------------------------------
function RPF_oil_ecl_Create()

  implicit none

  class(rel_perm_oil_owg_ecl_type), pointer :: RPF_oil_ecl_Create

  allocate(RPF_oil_ecl_Create)

  call RPF_oil_ecl_Create%Init()

end function RPF_oil_ecl_Create


subroutine RPF_oil_ecl_Init(this)

  implicit none

  class(rel_perm_oil_owg_ecl_type) :: this
  
  call RPFOilOWGBaseInit(this)
  
  this%Swco = 0.0
     
end subroutine RPF_oil_ecl_Init


subroutine RPF_oil_ecl_Verify(this,name,option)

  use Option_module

  implicit none

  class(rel_perm_oil_owg_ecl_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: Soil_max, krow_max,krog_max, dkr_dummy

  if (index(name,'PERMEABILITY_FUNCTION_OIL') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION_OIL,Eclipse.'
  endif

   if (Uninitialized(this%Swco)) then
     option%io_buffer = UninitializedMessage('WAT_CONNATE_SATURATION',string)
     call printErrMsg(option)          
   end if

   if ( .not.associated(this%rel_perm_ow) ) then
     option%io_buffer = trim(string) // ' ,RPF_OIL_WATER not defined'
     call printErrMsg(option)
   else
     call this%rel_perm_ow%Verify(string,option)
   end if

   if ( .not.associated(this%rel_perm_og) ) then
     option%io_buffer = trim(string) // ' ,RPF_OIL_GAS not defined'
     call printErrMsg(option)
   else
     call this%rel_perm_og%Verify(string,option)
   end if   

   Soil_max = 1.0 - this%Swco   

   call this%rel_perm_ow%RelativePermeability(Soil_max,krow_max, &
                                              dkr_dummy,option)
                                              
   call this%rel_perm_og%RelativePermeability(Soil_max,krog_max, &
                                               dkr_dummy,option)
                                               
   if ( krow_max /= krog_max ) then
     option%io_buffer = trim(string) // ' krow_max /= krog_max '
     call printErrMsg(option)
   end if 

   if ( this%rel_perm_ow%analytical_derivative_available .and. &
        this%rel_perm_og%analytical_derivative_available ) then
      this%analytical_derivative_available = PETSC_TRUE  
   end if

   !check what's the smallest between Sowcr and Sogcr and assign to Socr
   if ( this%rel_perm_ow%GetCriticalSaturation(option) <= &
        this%rel_perm_og%GetCriticalSaturation(option) ) then
     this%Socr = this%rel_perm_ow%GetCriticalSaturation(option)
   else 
     this%Socr = this%rel_perm_og%GetCriticalSaturation(option)
   end if      

   ! this is called at the end after Socr is intialised
   call RPFOilOWGBaseVerify(this,string,option)
   
end subroutine RPF_oil_ecl_Verify


subroutine RPF_oil_ecl_RelPerm(this,oil_sat,gas_sat,rel_perm, &
                                dkro_sato,dkro_satg,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_oil_owg_ecl_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(in) :: gas_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkro_sato
 PetscReal, intent(out) :: dkro_satg
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 PetscReal :: wat_sat, krow, dkrow_sato, krog, dkrog_sato
 PetscReal :: den, den2

 ! Form the Eclipse three-phase Kro expression

 call this%rel_perm_ow%RelativePermeability(oil_sat,krow,dkrow_sato,option, &
                                                                   table_idxs)
                                            
 call this%rel_perm_og%RelativePermeability(oil_sat,krog,dkrog_sato,option, &
                                                                    table_idxs)

 wat_sat =  1.0 - oil_sat - gas_sat
 
 den=gas_sat+ wat_sat - this%Swco
 if( den>0.0 ) then
   rel_perm=(gas_sat*krog+(wat_sat-this%Swco)*krow)/den
   den2 = den * den
   dkro_sato = gas_sat/den * dkrog_sato + &
               gas_sat * krog / den2 + & 
               (wat_sat - this%Swco) / den * dkrow_sato - &
               krow / den + &
               (wat_sat - this%Swco) * krow / den2
   dkro_satg = (krog - 1.0*krow ) / den 
 else
   rel_perm=0.5*(krog+krow)
   dkro_sato = 0.5*(dkrog_sato + dkrow_sato)
   dkro_satg = 0.0d0   
 endif

end subroutine RPF_oil_ecl_RelPerm


subroutine RelPermOW(this,oil_sat,rel_perm,dkr_sato,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_oil_owg_base_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkr_sato
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 call this%rel_perm_ow%RelativePermeability(oil_sat,rel_perm,dkr_sato,option)

end subroutine RelPermOW


subroutine RelPermOG(this,oil_sat,rel_perm,dkr_sato,option,table_idxs)

 use Option_module
 use Utility_module

 implicit none

 class(rel_perm_oil_owg_base_type) :: this
 PetscReal, intent(in) :: oil_sat
 PetscReal, intent(out) :: rel_perm
 PetscReal, intent(out) :: dkr_sato
 type(option_type), intent(inout) :: option
 PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

 call this%rel_perm_og%RelativePermeability(oil_sat,rel_perm,dkr_sato,option)

end subroutine RelPermOG

function RPFOilOWGBaseGetSowcr(this,option)

 use Option_module
 use Utility_module

 implicit none

 PetscReal :: RPFOilOWGBaseGetSowcr
 class(rel_perm_oil_owg_base_type) :: this
 type(option_type), intent(inout) :: option

 RPFOilOWGBaseGetSowcr = this%rel_perm_ow%GetCriticalSaturation(option)

end function RPFOilOWGBaseGetSowcr


function RPFOilOWGBaseGetSogcr(this,option)

 use Option_module
 use Utility_module

 implicit none

 PetscReal :: RPFOilOWGBaseGetSogcr
 class(rel_perm_oil_owg_base_type) :: this
 type(option_type), intent(inout) :: option

 RPFOilOWGBaseGetSogcr = this%rel_perm_og%GetCriticalSaturation(option)

end function RPFOilOWGBaseGetSogcr


subroutine RPF_oil_ecl_SetSwco(this,swco,option)

  use Option_module

  implicit none

  class(rel_perm_oil_owg_ecl_type) :: this
  PetscReal :: swco
  type(option_type), intent(inout) :: option

  this%Swco = swco
  
end subroutine RPF_oil_ecl_SetSwco

! ************************************************************************** !
! *********** END OWG Relative Permeability functions  ********************* !
! ************************************************************************** !

! ************************************************************************** !
! *********** Destroy Saturation and Relative Permeability functions  ****** !
! ************************************************************************** !

subroutine SaturationFunctionXWDestroy(sf_xw)
  !
  ! Destroys an XW saturuation function
  !
  ! Author: Paolo Orsini
  ! Date: 11/16/17 - 07/20/18
  !

  implicit none

  class(sat_func_xw_base_type), pointer :: sf_xw

  if (.not.associated(sf_xw)) return

  call SaturationFunctionDestroy(sf_xw%sat_func_sl)
  nullify(sf_xw%table)

  deallocate(sf_xw)
  nullify(sf_xw)

end subroutine SaturationFunctionXWDestroy

! ************************************************************************** !

subroutine SaturationFunctionOGDestroy(sf_og)
  !
  ! Destroys an OG saturuation function
  !
  ! Author: Paolo Orsini
  ! Date: 07/20/18
  !

  implicit none

  class(sat_func_og_base_type), pointer :: sf_og

  if (.not.associated(sf_og)) return

  call SaturationFunctionDestroy(sf_og%sat_func_sl)
  nullify(sf_og%table)

  deallocate(sf_og)
  nullify(sf_og)

end subroutine SaturationFunctionOGDestroy

! ************************************************************************** !


subroutine WatPermFunctionOWGDestroy(rpf)
  !
  ! Author: Paolo Orsini
  ! Date: 11/18/17
  !
  implicit none

  class(rel_perm_wat_owg_base_type), pointer :: rpf

  if (.not.associated(rpf)) return

  call rpf%Strip()

  deallocate(rpf)
  nullify(rpf)

end subroutine WatPermFunctionOWGDestroy


subroutine GasPermFunctionOWGDestroy(rpf)
  !
  ! Author: Paolo Orsini
  ! Date: 11/18/17
  !
  implicit none

  class(rel_perm_gas_owg_base_type), pointer :: rpf

  if (.not.associated(rpf)) return

  call rpf%Strip()

  deallocate(rpf)
  nullify(rpf)

end subroutine GasPermFunctionOWGDestroy


subroutine OWPermFunctionOWGDestroy(rpf)
  !
  ! Author: Paolo Orsini
  ! Date: 11/18/17
  !
  implicit none

  class(rel_perm_ow_owg_base_type), pointer :: rpf

  if (.not.associated(rpf)) return

  call rpf%Strip()

  deallocate(rpf)
  nullify(rpf)

end subroutine OWPermFunctionOWGDestroy


subroutine OGPermFunctionOWGDestroy(rpf)
  !
  ! Author: Paolo Orsini
  ! Date: 11/18/17
  !
  implicit none

  class(rel_perm_og_owg_base_type), pointer :: rpf

  if (.not.associated(rpf)) return

  call rpf%Strip()

  deallocate(rpf)
  nullify(rpf)

end subroutine OGPermFunctionOWGDestroy


subroutine OilPermFunctionOWGDestroy(rpf)
  !
  ! Author: Paolo Orsini
  ! Date: 11/18/17
  !
  implicit none

  class(rel_perm_oil_owg_base_type), pointer :: rpf

  if (.not.associated(rpf)) return

  call rpf%Strip()

  deallocate(rpf)
  nullify(rpf)

end subroutine OilPermFunctionOWGDestroy

! ************************************************************************** !

subroutine SetCCOWGPhaseFlags(option,oil_gas_interface_present, &
                              wat_gas_interface_present, gas_present, &
                              oil_perm_2ph_ow,oil_perm_3ph_owg)
  !
  ! Set phase flags for error checks based on PM moeds
  !
  ! Author: Paolo Orsini
  ! Date: 08/18/18
  !
  use Option_module

  implicit none

  type(option_type) :: option
  PetscBool, intent(out) :: oil_gas_interface_present
  PetscBool, intent(out) :: wat_gas_interface_present
  PetscBool, intent(out) :: gas_present
  PetscBool, intent(out) :: oil_perm_2ph_ow
  PetscBool, intent(out) :: oil_perm_3ph_owg

  oil_perm_2ph_ow = PETSC_FALSE
  oil_perm_3ph_owg = PETSC_FALSE
  gas_present = PETSC_FALSE
  oil_gas_interface_present = PETSC_FALSE
  wat_gas_interface_present = PETSC_FALSE

  select case(option%iflowmode)
    case(TOIL_IMS_MODE)
      oil_perm_2ph_ow = PETSC_TRUE
    case(TOWG_MODE)
      select case(option%iflow_sub_mode)
        case(TOWG_TODD_LONGSTAFF)
          oil_perm_2ph_ow = PETSC_TRUE
        case(TOWG_IMMISCIBLE)
          oil_perm_2ph_ow = PETSC_TRUE !oil_perm_3ph_owg no yet supported
          gas_present = PETSC_TRUE
          oil_gas_interface_present = PETSC_TRUE
        case default ! Black Oil and Solvent models
          oil_perm_3ph_owg = PETSC_TRUE
          gas_present = PETSC_TRUE
          oil_gas_interface_present = PETSC_TRUE
    end select
  end select

end subroutine SetCCOWGPhaseFlags

! ************************************************************************** !

end module Characteristic_Curves_OWG_module
