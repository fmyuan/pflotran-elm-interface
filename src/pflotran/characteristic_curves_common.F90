module Characteristic_Curves_Common_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Characteristic_Curves_Base_module
  use Dataset_Ascii_class

  implicit none

  private
  
!-----------------------------------------------------------------------------
!-- Saturation Functions -----------------------------------------------------
!-----------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_default_type
  contains
    procedure, public :: Verify => SFDefaultVerify
    procedure, public :: CapillaryPressure => SFDefaultCapillaryPressure
    procedure, public :: Saturation => SFDefaultSaturation
    procedure, public :: D2SatDP2 => SFDefaultD2SatDP2
  end type sat_func_default_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_constant_type
    PetscReal :: constant_capillary_pressure
    PetscReal :: constant_saturation
  contains
    procedure, public :: Verify => SFConstantVerify
    procedure, public :: CapillaryPressure => SFConstantCapillaryPressure
    procedure, public :: Saturation => SFConstantSaturation
    procedure, public :: D2SatDP2 => SFConstantD2SatDP2
  end type sat_func_constant_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_BC_type
    PetscReal :: alpha
    PetscReal :: lambda
  contains
    procedure, public :: Init => SFBCInit
    procedure, public :: Verify => SFBCVerify
    procedure, public :: SetupPolynomials => SFBCSetupPolynomials
    procedure, public :: CapillaryPressure => SFBCCapillaryPressure
    procedure, public :: Saturation => SFBCSaturation
    procedure, public :: D2SatDP2 => SFBCD2SatDP2
  end type sat_func_BC_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_Linear_type
    PetscReal :: alpha
  contains
    procedure, public :: Init => SFLinearInit
    procedure, public :: Verify => SFLinearVerify
    procedure, public :: CapillaryPressure => SFLinearCapillaryPressure
    procedure, public :: Saturation => SFLinearSaturation
    procedure, public :: D2SatDP2 => SFLinearD2SatDP2
  end type sat_func_Linear_type
  type, public, extends(sat_func_base_type) :: sat_func_mK_type
    PetscReal :: sigmaz, muz
    PetscReal :: rmax, r0
    PetscInt :: nparam
  contains
    procedure, public :: Init => SFmKInit
    procedure, public :: Verify => SFmKVerify
    procedure, public :: CapillaryPressure => SFmKCapillaryPressure
    procedure, public :: Saturation => SFmKSaturation
  end type sat_func_mK_type
!---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_IGHCC2_Comp_type
    PetscReal :: alpha
    PetscReal :: m
  contains
    procedure, public :: Init => SFIGHCC2CompInit
    procedure, public :: Verify => SFIGHCC2CompVerify
    procedure, public :: CapillaryPressure => SFIGHCC2CompCapillaryPressure
  end type sat_func_IGHCC2_Comp_type
!---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_Table_type
    class(dataset_ascii_type), pointer :: pc_dataset
  contains
    procedure, public :: Init => SFTableInit
    procedure, public :: Verify => SFTableVerify
    procedure, public :: CapillaryPressure => SFTableCapillaryPressure
    procedure, public :: Saturation => SFTableSaturation
  end type sat_func_Table_type
  
!-----------------------------------------------------------------------------
!-- Relative Permeability Functions ------------------------------------------
!-----------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rel_perm_func_default_type
  contains
    procedure, public :: Verify => RPFDefaultVerify
    procedure, public :: RelativePermeability => RPFDefaultRelPerm
  end type rel_perm_func_default_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_BC_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPFBurdineBCLiqInit
    procedure, public :: Verify => RPFBurdineBCLiqVerify
    procedure, public :: RelativePermeability => RPFBurdineBCLiqRelPerm
  end type rpf_Burdine_BC_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_BC_gas_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPFBurdineBCGasInit
    procedure, public :: Verify => RPFBurdineBCGasVerify
    procedure, public :: RelativePermeability => RPFBurdineBCGasRelPerm
  end type rpf_Burdine_BC_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_BC_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPFMualemBCLiqInit
    procedure, public :: Verify => RPFMualemBCLiqVerify
    procedure, public :: RelativePermeability => RPFMualemBCLiqRelPerm
  end type rpf_MUALEM_BC_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_BC_gas_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPFMualemBCGasInit
    procedure, public :: Verify => RPFMualemBCGasVerify
    procedure, public :: RelativePermeability => RPFMualemBCGasRelPerm
  end type rpf_Mualem_BC_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_Linear_liq_type
    PetscReal :: pcmax
    PetscReal :: alpha
  contains
    procedure, public :: Init => RPFMualemLinearLiqInit
    procedure, public :: Verify => RPFMualemLinearLiqVerify
    procedure, public :: RelativePermeability => RPFMualemLinearLiqRelPerm
  end type rpf_Mualem_Linear_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_Mualem_Linear_liq_type) :: & 
                        rpf_Mualem_Linear_gas_type
  contains
    procedure, public :: Init => RPFMualemLinearGasInit
    procedure, public :: Verify => RPFMualemLinearGasVerify
    procedure, public :: RelativePermeability => RPFMualemLinearGasRelPerm
  end type rpf_Mualem_Linear_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_Linear_liq_type
  contains
    procedure, public :: Init => RPFBurdineLinearLiqInit
    procedure, public :: Verify => RPFBurdineLinearLiqVerify
    procedure, public :: RelativePermeability => RPFBurdineLinearLiqRelPerm
  end type rpf_Burdine_Linear_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: & 
                        rpf_Burdine_Linear_gas_type
  contains
    procedure, public :: Init => RPFBurdineLinearGasInit
    procedure, public :: Verify => RPFBurdineLinearGasVerify
    procedure, public :: RelativePermeability => RPFBurdineLinearGasRelPerm
  end type rpf_Burdine_Linear_gas_type  
  !---------------------------------------------------------------------------
  ! Constant: for running tests with a fixed relative permeability
  type, public, extends(rel_perm_func_base_type) :: rel_perm_func_constant_type
    PetscReal :: kr
  contains
    procedure, public :: Verify => RPFConstantVerify
    procedure, public :: RelativePermeability => RPFConstantRelPerm
  end type rel_perm_func_constant_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_mK_liq_type
    PetscReal :: sigmaz
  contains
    procedure, public :: Init => RPFmKLiqInit
    procedure, public :: Verify => RPFmKLiqVerify
    procedure, public :: RelativePermeability => RPFmKLiqRelPerm
  end type rpf_mK_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_mK_gas_type
    PetscReal :: sigmaz
  contains
    procedure, public :: Verify => RPFmKGasVerify
    procedure, public :: RelativePermeability => RPFmKGasRelPerm
  end type rpf_mK_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: &
                                     rpf_IGHCC2_Comp_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPFIGHCC2CompLiqInit
    procedure, public :: Verify => RPFIGHCC2CompLiqVerify
    procedure, public :: RelativePermeability => &
                                  RPFIGHCC2CompLiqRelPerm
  end type rpf_IGHCC2_Comp_liq_type  
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: &
                                       rpf_IGHCC2_Comp_gas_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPFIGHCC2CompGasInit
    procedure, public :: Verify => RPFIGHCC2CompGasVerify
    procedure, public :: RelativePermeability => &
                                  RPFIGHCC2CompGasRelPerm
  end type rpf_IGHCC2_Comp_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: &
                                     rpf_Table_liq_type
    class(dataset_ascii_type), pointer :: rpf_dataset
  contains
    procedure, public :: Init => RPFTableLiqInit
    procedure, public :: Verify => RPFTableLiqVerify
    procedure, public :: RelativePermeability => &
                                  RPFTableLiqRelPerm
  end type rpf_Table_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: &
                                       rpf_Table_gas_type
    class(dataset_ascii_type), pointer :: rpf_dataset
  contains
    procedure, public :: Init => RPFTableGasInit
    procedure, public :: Verify => RPFTableGasVerify
    procedure, public :: RelativePermeability => &
                                  RPFTableGasRelPerm
  end type rpf_Table_gas_type

 
  public :: &! standard char. curves:
            SFDefaultCreate, &
            SFConstantCreate, &
            SFBCCreate, &
            SFLinearCreate, &
            SFmKCreate, &
            SFIGHCC2CompCreate, &
            SFTableCreate, &
            ! standard rel. perm. curves:
            RPFDefaultCreate, &
            RPFConstantCreate, &  
            RPFBurdineBCLiqCreate, &
            RPFBurdineBCGasCreate, &
            RPFMualemBCLiqCreate, &
            RPFMualemBCGasCreate, &
            RPFMualemLinearLiqCreate, &
            RPFMualemLinearGasCreate, &
            RPFBurdineLinearLiqCreate, &
            RPFBurdineLinearGasCreate, &
            RPFmKLiqCreate, &
            RPFmKGasCreate, &
            RPFIGHCC2CompLiqCreate, &
            RPFIGHCC2CompGasCreate, &
            RPFTableLiqCreate, &
            RPFTableGasCreate
  
contains

! ************************************************************************** !
! ************************************************************************** !

function SFDefaultCreate()

  ! Creates the default saturation function object

  implicit none
  
  class(sat_func_default_type), pointer :: SFDefaultCreate
  
  allocate(SFDefaultCreate)
  call SFBaseInit(SFDefaultCreate)
  SFDefaultCreate%Sr = 0.d0
  
  SFDefaultCreate%analytical_derivative_available = PETSC_TRUE
  
end function SFDefaultCreate

! ************************************************************************** !

subroutine SFDefaultVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_default_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  option%io_buffer = 'A default Saturation Function has been chosen in ' // &
    trim(name) // '.'
  call PrintWrnMsg(option)
  
end subroutine SFDefaultVerify

! ************************************************************************** !

subroutine SFDefaultCapillaryPressure(this,liquid_saturation, &
                                      capillary_pressure,dpc_dsatl,option)
  use Option_module
  
  implicit none
  
  class(sat_func_default_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation < 1.d0) then
    option%io_buffer = 'SFDefaultCapillaryPressure is a dummy routine used &
      &for saturated flow only.  The user must specify a valid &
      &SATURATION_FUNCTION.'
    call PrintErrMsgByRank(option)
  endif

end subroutine SFDefaultCapillaryPressure

! ************************************************************************** !

subroutine SFDefaultSaturation(this,capillary_pressure, &
                               liquid_saturation,dsat_dpres,option)
  use Option_module

  implicit none
  
  class(sat_func_default_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFDefaultSaturation is a dummy routine used &
    &for saturated flow only.  The user must specify a valid &
    &SATURATION_FUNCTION.'
  call PrintErrMsgByRank(option)

end subroutine SFDefaultSaturation

! ************************************************************************** !

subroutine SFDefaultD2SatDP2(this,pc, &
                               d2s_dp2,option)
  use Option_module

  implicit none
  
  class(sat_func_default_type) :: this
  PetscReal, intent(in) :: pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFDefaultD2SatDP2 is a dummy routine used &
    &for saturated flow only.  The user must specify a valid &
    &SATURATION_FUNCTION.'
  call PrintErrMsgByRank(option)

end subroutine SFDefaultD2SatDP2

! ************************************************************************** !

function RPFDefaultCreate()

  ! Creates the default relative permeability function object

  implicit none
  
  class(rel_perm_func_default_type), pointer :: RPFDefaultCreate
  
  allocate(RPFDefaultCreate)
  call RPFBaseInit(RPFDefaultCreate)
  RPFDefaultCreate%Sr = 0.d0
  
end function RPFDefaultCreate

! ************************************************************************** !

subroutine RPFDefaultVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rel_perm_func_default_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  option%io_buffer = 'A default Relative Permeability Function has been ' // &
    'chosen in ' // trim(name) // '.'
  call PrintWrnMsg(option)

end subroutine RPFDefaultVerify

! ************************************************************************** !

subroutine RPFDefaultRelPerm(this,liquid_saturation,relative_permeability, &
                            dkr_sat,option)
  use Option_module

  implicit none
  
  class(rel_perm_func_default_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation < 1.d0) then
    option%io_buffer = 'RPFDefaultRelPerm is a dummy routine used &
      &for saturated flow only.  The user must specify a valid &
      &PERMEABILITY_FUNCTION.'
    call PrintErrMsgByRank(option)
  endif
  relative_permeability = 1.d0
  
end subroutine RPFDefaultRelPerm

! ************************************************************************** !
! ************************************************************************** !

function SFConstantCreate()

  ! Creates the default saturation function object

  implicit none
  
  class(sat_func_constant_type), pointer :: SFConstantCreate
  
  allocate(SFConstantCreate)
  call SFBaseInit(SFConstantCreate)
  ! set Sr to zero as it doesn't matter, but must be initialized
  SFConstantCreate%Sr = 0.d0 
  SFConstantCreate%constant_capillary_pressure = UNINITIALIZED_DOUBLE
  SFConstantCreate%constant_saturation = UNINITIALIZED_DOUBLE
  
  SFConstantCreate%analytical_derivative_available = PETSC_TRUE
  
end function SFConstantCreate

! ************************************************************************** !

subroutine SFConstantVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_constant_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  character(len=MAXSTRINGLENGTH) :: string  

  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,CONSTANT'
  endif
  call SFBaseVerify(this,string,option)
  select case(option%iflowmode)
    case(RICHARDS_MODE,RICHARDS_TS_MODE,TH_MODE,TH_TS_MODE)
      if (Initialized(this%constant_capillary_pressure)) then
        option%io_buffer = 'CONSTANT_CAPILLARY_PRESSURE is not supported for &
          &Richards or TH flow modes as CONSTANT_SATURATION must be applied. &
          &See ' // trim(string) // '.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%constant_saturation)) then
        option%io_buffer = 'CONSTANT_SATURATION must be specified for ' // &
          trim(string) // '.'
        call PrintErrMsg(option)
      endif
    case(WF_MODE,G_MODE,MPH_MODE,H_MODE)
      if (Initialized(this%constant_saturation)) then
        option%io_buffer = 'CONSTANT_SATURATION is not supported for &
          &multiphase flow modes as CONSTANT_CAPILLARY_PRESSURE must be &
          &applied. Saturation is a primary dependent variables. &
          &See ' // trim(string) // '.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%constant_capillary_pressure)) then
        option%io_buffer = 'CONSTANT_CAPILLARY_PRESSURE must be specified &
          &for ' // trim(string) // '.'
        call PrintErrMsg(option)
      endif
    case default
  end select

end subroutine SFConstantVerify

! ************************************************************************** !

subroutine SFConstantCapillaryPressure(this,liquid_saturation, &
                                       capillary_pressure,dpc_dsatl,option)
  use Option_module
  
  implicit none
  
  class(sat_func_constant_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  dpc_dsatl = 0.d0
  capillary_pressure = this%constant_capillary_pressure

end subroutine SFConstantCapillaryPressure

! ************************************************************************** !

subroutine SFConstantSaturation(this,capillary_pressure, &
                                liquid_saturation,dsat_dpres,option)
  use Option_module

  implicit none
  
  class(sat_func_constant_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option

  liquid_saturation = this%constant_saturation
  dsat_dpres = 0.d0

end subroutine SFConstantSaturation

! ************************************************************************** !

subroutine SFConstantD2SatDP2(this,pc,d2s_dp2,option)
  use Option_module

  implicit none
  
  class(sat_func_constant_type) :: this
  PetscReal, intent(in) :: pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option
  
  d2s_dp2 = 0.d0

end subroutine SFConstantD2SatDP2

! ************************************************************************** !

function RPFConstantCreate()

  ! Creates the constant relative permeability function object

  implicit none
  
  class(rel_perm_func_constant_type), pointer :: RPFConstantCreate
  
  allocate(RPFConstantCreate)
  call RPFBaseInit(RPFConstantCreate)
  ! set Sr = 0. to avoid uninitialized failure
  RPFConstantCreate%Sr = 0.d0
  RPFConstantCreate%kr = 0.d0
  
end function RPFConstantCreate

! ************************************************************************** !

subroutine RPFConstantVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rel_perm_func_constant_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,CONSTANT'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%kr)) then
    option%io_buffer = UninitializedMessage('RELATIVE_PERMEABILITY',string)
    call PrintErrMsg(option)
  endif   

end subroutine RPFConstantVerify

! ************************************************************************** !

subroutine RPFConstantRelPerm(this,liquid_saturation,relative_permeability, &
                            dkr_sat,option)
  use Option_module

  implicit none
  
  class(rel_perm_func_constant_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  relative_permeability = this%kr
  dkr_sat = 0.d0
  
end subroutine RPFConstantRelPerm

! ************************************************************************** !

function SFIGHCC2CompCreate()

  ! Creates the IGHCC2 Comparison capillary pressure function object

  implicit none

  class(sat_func_IGHCC2_Comp_type), pointer :: &
                              SFIGHCC2CompCreate

  allocate(SFIGHCC2CompCreate)
  call SFIGHCC2CompCreate%Init()

end function SFIGHCC2CompCreate

! ************************************************************************** !

subroutine SFIGHCC2CompInit(this)

  ! Creates the IGHCC2 Comparison capillary pressure function object

  implicit none

  class(sat_func_IGHCC2_Comp_type) :: this

  call SFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%m = UNINITIALIZED_DOUBLE

  this%analytical_derivative_available = PETSC_TRUE

end subroutine SFIGHCC2CompInit


! ************************************************************************** !

subroutine SFIGHCC2CompVerify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_IGHCC2_Comp_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,IGHCC2 Comp'
  endif
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call PrintErrMsg(option)
  endif

end subroutine SFIGHCC2CompVerify

! ************************************************************************** !

subroutine SFIGHCC2CompCapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,dpc_dsatl,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation, adapted to
  ! benchmark against the IGHCC2 study.
  ! 
  ! Author: Michael Nole
  ! Date: 05/16/19
  !
  use Option_module

  implicit none

  class(sat_func_IGHCC2_Comp_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option

  PetscReal :: n
  PetscReal :: Se

  PetscReal :: neg_one_over_m
  PetscReal :: one_over_n
  PetscReal :: dSe_dsatl
  PetscReal :: Se_sup_neg_one_over_m
  PetscReal :: Se_sup_neg_one_over_m_minus_one

  dpc_dsatl = 0.d0

  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif

  n = 1.d0/(1.d0-this%m)
  neg_one_over_m = -1.d0/this%m
  one_over_n = 1.d0/n
  dSe_dsatl = 1.d0 / (1.d0-this%Sr)
  Se = (liquid_saturation-this%Sr)*dSe_dsatl
  Se_sup_neg_one_over_m = Se**neg_one_over_m
  Se_sup_neg_one_over_m_minus_one = Se_sup_neg_one_over_m - 1.d0
  
  capillary_pressure = (Se_sup_neg_one_over_m_minus_one**this%m)/this%alpha
  
  dpc_dsatl = capillary_pressure/Se_sup_neg_one_over_m_minus_one * &
              one_over_n * neg_one_over_m * Se_sup_neg_one_over_m / Se * &
              dSe_dsatl

  if (capillary_pressure > this%pcmax) then
    capillary_pressure = this%pcmax
    dpc_dsatl = 0.d0
  endif

end subroutine SFIGHCC2CompCapillaryPressure

! ************************************************************************** !
! ************************************************************************** !

function SFTableCreate()

  ! Creates the Lookup Table capillary pressure function object

  implicit none

  class(sat_func_Table_type), pointer :: &
                              SFTableCreate

  allocate(SFTableCreate)
  call SFTableCreate%Init()

end function SFTableCreate

! ************************************************************************** !

subroutine SFTableInit(this)

  ! Creates the Lookup Table capillary pressure function object

  implicit none

  class(sat_func_Table_type) :: this

  call SFBaseInit(this)
  this%pc_dataset => DatasetAsciiCreate()

  this%analytical_derivative_available = PETSC_TRUE

end subroutine SFTableInit


! ************************************************************************** !

subroutine SFTableVerify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_Table_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION, Lookup Table'
  endif
  call SFBaseVerify(this,string,option)
  if (.not.associated(this%pc_dataset)) then
    option%io_buffer = UninitializedMessage('TABLE',string)
    call PrintErrMsg(option)
  endif

end subroutine SFTableVerify

! ************************************************************************** !

subroutine SFTableCapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,dpc_dsatl,option)
  ! 
  ! Computes the capillary pressure as a function of saturation.
  ! 
  ! Author: Michael Nole
  ! Date: 01/15/20
  !
  use Option_module

  implicit none

  class(sat_func_Table_type) :: this
  PetscReal, intent(in) :: liquid_saturation

  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option

  class(dataset_ascii_type), pointer :: dataset
  PetscReal, pointer :: times(:)
  PetscInt :: i, j, num_entries

  dataset => this%pc_dataset
  times => dataset%time_storage%times
  num_entries = 0
  j = 0
  do i = 1,size(times)
    if (times(i) <= dataset%time_storage%cur_time) then
      if (i > 1 .and. times(i) > times(i-1)) then
        j = 0
        num_entries = 0
      endif
      if (j==0) j = i
      num_entries = num_entries + 1
    endif
  enddo


  if (liquid_saturation < dataset%rbuffer(2*j-1)) then
    capillary_pressure = dataset%rbuffer(2*j)
    dpc_dsatl = 0.d0
  elseif (liquid_saturation > dataset%rbuffer(2*(j-1+num_entries)-1)) then
    dpc_dsatl = (capillary_pressure - dataset%rbuffer(2*(j-1+num_entries))) / &
               (liquid_saturation - dataset%rbuffer(2*(j-1+num_entries)-1))
    capillary_pressure = (liquid_saturation - dataset% &
                         rbuffer(2*(j-1+num_entries)-1)) * dpc_dsatl + &
                         dataset%rbuffer(2*(j-1+num_entries))
  else
    do i = j+1, j+num_entries-1
      if (liquid_saturation <= dataset%rbuffer(2*i-1)) then
        dpc_dsatl = (dataset%rbuffer(2*i) - dataset%rbuffer(2*i-2)) / &
                    (dataset%rbuffer(2*i-1) - dataset%rbuffer(2*i-3))
        capillary_pressure = (liquid_saturation - dataset%rbuffer(2*i-3)) * &
                         dpc_dsatl + dataset%rbuffer(2*i-2)
        exit
      endif
    enddo
  endif

  if (capillary_pressure > this%pcmax) then
    capillary_pressure = this%pcmax
    dpc_dsatl = 0.d0
  endif

end subroutine SFTableCapillaryPressure

! ************************************************************************** !

subroutine SFTableSaturation(this,capillary_pressure, &
                            liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes saturation as a function of capillary pressure.
  ! 
  ! Author: Michael Nole
  ! Date: 02/04/20
  !
  use Option_module

  implicit none

  class(sat_func_Table_type) :: this
  PetscReal, intent(in) :: capillary_pressure

  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option

  class(dataset_ascii_type), pointer :: dataset
  PetscReal, pointer :: times(:)
  PetscInt :: i, j, num_entries

  dataset => this%pc_dataset
  times => dataset%time_storage%times
  num_entries = 0
  j = 0
  do i = 1,size(times)
    if (times(i) <= dataset%time_storage%cur_time) then
      if (i > 1 .and. times(i) > times(i-1)) then
        j = 0
        num_entries = 0
      endif
      if (j==0) j = i
      num_entries = num_entries + 1
    endif
  enddo

  if (capillary_pressure >  dataset%rbuffer(2*j)) then
    dsat_dpres = (dataset%rbuffer(2*j+1) - dataset%rbuffer(2*j-1)) / &
                 (dataset%rbuffer(2*j+2) - dataset%rbuffer(2*j))
    liquid_saturation = dataset%rbuffer(2*j-1) - dsat_dpres * &
                        (capillary_pressure - dataset%rbuffer(2*j))
  elseif (capillary_pressure < dataset%rbuffer(2*(j-1+num_entries))) then
    dsat_dpres = (liquid_saturation - dataset%rbuffer(2*(j-1+num_entries)-1))/ &
                 (capillary_pressure - dataset%rbuffer(2*(j-1+num_entries)))
    liquid_saturation = dsat_dpres * (0.d0 - &
                        dataset%rbuffer(2*(j-1+num_entries))) + &
                        dataset%rbuffer(2*(j-1+num_entries))
  else
    do i = j+1, j+num_entries-1
      if (capillary_pressure >= dataset%rbuffer(2*i)) then
        dsat_dpres = (dataset%rbuffer(2*i-1) - dataset%rbuffer(2*i-3)) / &
                    (dataset%rbuffer(2*i) - dataset%rbuffer(2*i-2))
        liquid_saturation = (dataset%rbuffer(2*i) - capillary_pressure) * &
                            dsat_dpres  + dataset%rbuffer(2*i-1)
        exit
      endif
    enddo
  endif

  liquid_saturation = maxval([0.d0,liquid_saturation])
  liquid_saturation = minval([1.d0,liquid_saturation])

end subroutine SFTableSaturation

! ************************************************************************** !


function SFBCCreate()

  ! Creates the Brooks Corey capillary pressure function object

  implicit none
  
  class(sat_func_BC_type), pointer :: SFBCCreate
  
  allocate(SFBCCreate)
  call SFBCCreate%Init()
  
end function SFBCCreate

! ************************************************************************** !

subroutine SFBCInit(this)

  use Option_module

  implicit none
  
  class(sat_func_BC_type) :: this

  call SFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine SFBCInit

! ************************************************************************** !

subroutine SFBCVerify(this,name,option)

  use Option_module

  implicit none
  
  class(sat_func_BC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BROOKS_COREY'
  endif  
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call PrintErrMsg(option)
  endif 
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif 
  
end subroutine SFBCVerify

! ************************************************************************** !

subroutine SFBCSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing Brooks-Corey saturation function

  use Option_module
  use Utility_module
  
  implicit none
  
  class(sat_func_BC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  PetscReal :: b(4)

  ! polynomial fitting pc as a function of saturation
  ! 1.05 is essentially pc*alpha (i.e. pc = 1.05/alpha)
  this%sat_poly => PolynomialCreate()
  this%sat_poly%low = 1.05d0**(-this%lambda)
  this%sat_poly%high = 1.d0
  
  b = 0.d0
  ! fill right hand side
  ! capillary pressure at 1
  b(1) = 1.05d0/this%alpha 
  ! capillary pressure at 2
  b(2) = 0.d0
  ! derivative of pressure at saturation_1
  ! pc = Se**(-1/lambda)/alpha
  ! dpc_dSe = -1/lambda*Se**(-1/lambda-1)/alpha
  b(3) = -1.d0/this%lambda* &
          this%sat_poly%low**(-1.d0/this%lambda-1.d0)/ &
          this%alpha

  call QuadraticPolynomialSetup(this%sat_poly%low,this%sat_poly%high,b(1:3), &
                                ! indicates derivative given at 1
                                PETSC_TRUE) 
      
  this%sat_poly%coefficients(1:3) = b(1:3)

  ! polynomial fitting saturation as a function of pc
  !geh: cannot invert the pressure/saturation relationship above
  !     since it can result in saturations > 1 with both
  !     quadratic and cubic polynomials
  ! fill matix with values
  this%pres_poly => PolynomialCreate()
  this%pres_poly%low = 0.95/this%alpha
  this%pres_poly%high = 1.05/this%alpha
  
  b = 0.d0
  ! Se at 1
  b(1) = 1.d0
  ! Se at 2
  b(2) = (this%pres_poly%high*this%alpha)** &
          (-this%lambda)
  ! derivative of Se at 1
  b(3) = 0.d0 
  ! derivative of Se at 2
  b(4) = -this%lambda/this%pres_poly%high* &
            (this%pres_poly%high*this%alpha)** &
              (-this%lambda)

  call CubicPolynomialSetup(this%pres_poly%low,this%pres_poly%high,b)

  this%pres_poly%coefficients(1:4) = b(1:4)
  
  
end subroutine SFBCSetupPolynomials

! ************************************************************************** !

subroutine SFBCCapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,dpc_dsatl,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation using the
  ! Brooks-Corey formulation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  !
  use Option_module
  use Utility_module
  
  implicit none
  
  class(sat_func_BC_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dSe_dsatl
  PetscReal :: dpc_dSe
  PetscReal :: neg_one_over_lambda
  PetscReal :: Pcmax_copy

  dpc_dsatl = 0.d0
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  dSe_dsatl = 1.d0 / (1.d0-this%Sr)
  Se = (liquid_saturation-this%Sr)*dSe_dsatl
  if (associated(this%sat_poly)) then
    if (Se > this%sat_poly%low) then
      call QuadraticPolynomialEvaluate(this%sat_poly%coefficients(1:3), &
                                       Se,capillary_pressure,dpc_dSe)
      dpc_dsatl = dpc_dSe*dSe_dsatl
      return
    endif
  endif
  neg_one_over_lambda = -1.d0/this%lambda
  capillary_pressure = (Se**neg_one_over_lambda)/this%alpha
  dpc_dsatl = capillary_pressure/Se*neg_one_over_lambda*dSe_dsatl

#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
    dpc_dsatl = dpc_satl*(1.d0-liquid_saturation)/0.001d0 + &
                capillary_pressure*(-1.d0/0.001d0)
  endif
#endif  

  if (capillary_pressure > this%pcmax) then
    capillary_pressure = this%pcmax
    dpc_dsatl = 0.d0
  endif
  
end subroutine SFBCCapillaryPressure

! ************************************************************************** !

subroutine SFBCSaturation(this,capillary_pressure, &
                            liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
    
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_BC_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: pc_alpha_neg_lambda
  PetscReal :: Se
  PetscReal :: dSe_dpc
  PetscReal, parameter :: dpc_dpres = -1.d0
  
  dsat_dpres = 0.d0
  
  ! reference #1
  if (associated(this%pres_poly)) then
    if (capillary_pressure < this%pres_poly%low) then
      liquid_saturation = 1.d0
      return
    else if (capillary_pressure < this%pres_poly%high) then
      call CubicPolynomialEvaluate(this%pres_poly%coefficients, &
                                   capillary_pressure,Se,dSe_dpc)
      liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
      dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
      return
    endif
  else
    if (capillary_pressure < 1.d0/this%alpha) then
      liquid_saturation = 1.d0
      dsat_dpres = 0.d0
      return
    endif
  endif

  pc_alpha_neg_lambda = (capillary_pressure*this%alpha)**(-this%lambda)
  Se = pc_alpha_neg_lambda
  dSe_dpc = -this%lambda/capillary_pressure*pc_alpha_neg_lambda
  liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
  dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
  
end subroutine SFBCSaturation

! ************************************************************************** !

subroutine SFBCD2SatDP2(this,pc,d2s_dp2,option)

  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_BC_type) :: this
  PetscReal, intent(in) :: pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option
  
  PetscReal :: pc_alpha_neg_lambda
  PetscReal :: d2Se_dpc2
  PetscReal, parameter :: dpc_dpres = -1.d0
  
  ! reference #1
  if (associated(this%pres_poly)) then
    if (pc < this%pres_poly%low) then
      d2s_dp2 = 0.d0
      return
    else if (pc < this%pres_poly%high) then
      d2Se_dpc2 = this%pres_poly%coefficients(3)*2.d0 + &
                  this%pres_poly%coefficients(4)*6.d0*pc
      d2s_dp2 = (1.d0-this%Sr)*d2Se_dpc2*dpc_dpres*dpc_dpres
      return
    endif
  else
    if (pc < 1.d0/this%alpha) then
      d2s_dp2 = 0.d0
      return
    endif
  endif

  pc_alpha_neg_lambda = (pc*this%alpha)**(-this%lambda)
  d2Se_dpc2 = (this%lambda*this%lambda + this%lambda)/(pc*2.d0)*pc_alpha_neg_lambda
  d2s_dp2 = (1.d0-this%Sr)*d2Se_dpc2*dpc_dpres*dpc_dpres
  
end subroutine SFBCD2SatDP2

! ************************************************************************** !

function SFLinearCreate()

  ! Creates the Linear capillary pressure function object

  implicit none
  
  class(sat_func_Linear_type), pointer :: SFLinearCreate
  
  allocate(SFLinearCreate)
  call SFLinearCreate%Init()
  
end function SFLinearCreate

! ************************************************************************** !

subroutine SFLinearInit(this)

  ! Creates the Linear capillary pressure function object

  implicit none
  
  class(sat_func_Linear_type) :: this

  call SFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine SFLinearInit

! ************************************************************************** !

subroutine SFLinearVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_Linear_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,LINEAR'
  endif
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call PrintErrMsg(option)
  endif   

end subroutine SFLinearVerify

! ************************************************************************** !

subroutine SFLinearCapillaryPressure(this,liquid_saturation, &
                                       capillary_pressure,dpc_dsatl,option)
  ! 
  ! Computes the capillary pressure as a function of saturation.
  !
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  !
  use Option_module
  use Material_Aux_class
  
  implicit none
  
  class(sat_func_Linear_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dSe_dsatl
  PetscReal :: one_over_alpha_minus_pcmax

  dpc_dsatl = 0.d0
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  dSe_dsatl = 1.d0/(1.d0-this%Sr)
  Se = (liquid_saturation-this%Sr)*dSe_dsatl
  one_over_alpha_minus_pcmax = 1.d0/this%alpha-this%pcmax
  capillary_pressure = one_over_alpha_minus_pcmax*Se + this%pcmax
  dpc_dsatl = one_over_alpha_minus_pcmax*dSe_dsatl

#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
    dpc_dsatl = dpc_satl*(1.d0-liquid_saturation)/0.001d0 + &
                capillary_pressure*(-1.d0/0.001d0)
  endif
#endif  

  if (capillary_pressure > this%pcmax) then
    capillary_pressure = this%pcmax
    dpc_dsatl = 0.d0
  endif
  
end subroutine SFLinearCapillaryPressure

! ************************************************************************** !

subroutine SFLinearSaturation(this,capillary_pressure, &
                                liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
  !   
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_Linear_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dSe_dpc
  PetscReal, parameter :: dpc_dpres = -1.d0
  
  dsat_dpres = 0.d0

  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  else
    Se = (this%pcmax-capillary_pressure) / (this%pcmax-1.d0/this%alpha)
    dSe_dpc = -1.d0/(this%pcmax-1.d0/this%alpha)
    liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
    dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
  endif 

end subroutine SFLinearSaturation

! ************************************************************************** !

subroutine SFLinearD2SatDP2(this,pc,d2s_dp2,option)

  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_Linear_type) :: this
  PetscReal, intent(in) :: pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dSe_dpc
  PetscReal, parameter :: dpc_dpres = -1.d0
  
  d2s_dp2 = 0.d0

end subroutine SFLinearD2SatDP2

! ************************************************************************** !

function SFmKCreate()

  ! Creates the modified Kosugi saturation function object

  implicit none

  class(sat_func_mK_type), pointer :: SFmKCreate

  allocate(SFmKCreate)
  call SFmKCreate%Init()

end function SFmKCreate

! ************************************************************************** !

subroutine SFmKInit(this)

  ! Initializes modified Kosugi saturation function object

  implicit none
  
  class(sat_func_mK_type) :: this

  call SFBaseInit(this)
  this%sigmaz = UNINITIALIZED_DOUBLE
  this%muz = UNINITIALIZED_DOUBLE
  this%rmax = UNINITIALIZED_DOUBLE
  this%r0 = UNINITIALIZED_DOUBLE
  this%nparam = UNINITIALIZED_INTEGER
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine SFmKInit

! ************************************************************************** !

subroutine SFmKVerify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_mK_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,MODIFIED_KOSUGI'
  endif
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%sigmaz)) then
    option%io_buffer = UninitializedMessage('SIGMAZ',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%muz)) then
    option%io_buffer = UninitializedMessage('MUZ',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%nparam)) then
    option%io_buffer = UninitializedMessage('NPARAM',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%rmax)) then
    ! rmax is used for both nparam 3 and 4
    option%io_buffer = UninitializedMessage('RMAX',string)
    call PrintErrMsg(option)
  endif
  select case(this%nparam)
    case(4)
      ! r0 is only used for nparam 4
      if (Uninitialized(this%r0)) then
        option%io_buffer = UninitializedMessage('R0',string)
        call PrintErrMsg(option)
      endif
      if (this%r0 >= this%rmax) then
        option%io_buffer = trim(string) // ' requires RMAX > R0'
        call PrintErrMsg(option)
      end if
    case(3)
      continue ! rmax handled above
    case default
      option%io_buffer = 'invalid NPARAM value in' // &
        trim(string) // '. Only NPARAM=(3,4) supported.'
      call PrintErrMsg(option)
  end select

end subroutine SFmKVerify

! ************************************************************************** !

subroutine SFmKCapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,dpc_dsatl,option)
  !
  ! Computes the capillary_pressure as a function of saturation
  ! for modified Kosugi model.
  !
  ! Malama, B. & K.L. Kuhlman, 2015. Unsaturated Hydraulic Conductivity
  ! Models Based on Truncated Lognormal Pore-size Distributions, Groundwater,
  ! 53(3):498-502. http://dx.doi.org/10.1111/gwat.12220
  !
  ! Author: Kris Kuhlman
  ! Date: 2017
  !
  use Option_module
  use Utility_module, only : InverseNorm

  implicit none

  PetscReal, parameter :: KAPPA = 1.49D-1 !  water in glass tube
  PetscReal, parameter :: LNKAP = log(KAPPA)
  PetscReal, parameter :: UNIT_CONVERSION = 9.982D+2*9.81d0/1.0D+2
  
  class(sat_func_mK_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option

  PetscReal :: Se
  PetscReal :: inverse, exparg
  PetscReal :: hc, hmaxinv
  PetscReal :: dinverse_dSe
  PetscReal :: dSe_dsatl, dexparg_dinverse, dpc_dexparg
  PetscReal :: one_over_pc
  PetscReal :: tempreal

  dpc_dsatl = 0.d0

  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif

  dSe_dsatl = 1.d0/(1.d0 - this%Sr)
  Se = (liquid_saturation - this%Sr)*dSe_dsatl
!  inverse = -InverseNorm(Se)
  call InverseNorm(Se,inverse,PETSC_TRUE,dinverse_dSe)
  inverse = -1.d0*inverse
  dinverse_dSe = -1.d0*dinverse_dSe
  exparg = this%sigmaz*inverse + LNKAP - this%muz
  dexparg_dinverse = this%sigmaz

  hc = KAPPA/this%rmax
  dpc_dexparg = exp(exparg)
  capillary_pressure = dpc_dexparg + hc
  dpc_dsatl = dpc_dexparg*dexparg_dinverse*dinverse_dSe*dSe_dsatl
  if (this%nparam == 4) then
    hmaxinv = this%r0/KAPPA
    one_over_pc = 1.d0/capillary_pressure
    tempreal = 1.d0/(one_over_pc + hmaxinv)
    capillary_pressure = tempreal
    dpc_dsatl = capillary_pressure*tempreal*one_over_pc*one_over_pc*dpc_dsatl
  end if

  capillary_pressure = capillary_pressure*UNIT_CONVERSION
  dpc_dsatl = dpc_dsatl*UNIT_CONVERSION
  if (capillary_pressure > this%pcmax) then
    capillary_pressure = this%pcmax
    dpc_dsatl = 0.d0
  endif

end subroutine SFmKCapillaryPressure

! ************************************************************************** !

subroutine SFmKSaturation(this,capillary_pressure, &
                            liquid_saturation,dsat_dpres,option)
  !
  ! Computes the saturation (and associated derivatives) as a function of
  ! capillary pressure for modified Kosugi model
  !
  ! Malama, B. & K.L. Kuhlman, 2015. Unsaturated Hydraulic Conductivity
  ! Models Based on Truncated Lognormal Pore-size Distributions, Groundwater,
  ! 53(3):498-502. http://dx.doi.org/10.1111/gwat.12220
  !
  ! Author: Kris Kuhlman
  ! Date: 2017
  !
  use Option_module
  use Utility_module, only : InverseNorm

  implicit none

  ! gnu & intel extension and required in f2008
  intrinsic :: erfc

  PetscReal, parameter :: KAPPA = 1.49D-1 ! water in glass tube
  PetscReal, parameter :: LNKAP = log(KAPPA)
  PetscReal, parameter :: SQRT2 = sqrt(2.0d0)
  PetscReal, parameter :: SQRTPI = sqrt(4.0d0*atan(1.0d0))
  PetscReal, parameter :: UNIT_CONVERSION = 9.982D+2*9.81d0/1.0D+2
  
  class(sat_func_mK_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option

  PetscReal :: hc, hmax, cap_press_scaled
  PetscReal :: rt2sz
  PetscReal :: lnArg, erfcArg

  dsat_dpres = 0.0d0
  cap_press_scaled = capillary_pressure/UNIT_CONVERSION
  
  hc = KAPPA/this%rmax
  if (cap_press_scaled <= hc) then
    liquid_saturation = 1.d0
    return
  end if

  if (this%nparam == 3) then
    lnArg = cap_press_scaled - hc
  else ! nparam == 4 
    hmax = KAPPA/this%r0
    if (cap_press_scaled >= hmax) then
      liquid_saturation = this%Sr
      return
    end if
    lnArg = 1.d0/(1.d0/cap_press_scaled - 1.d0/hmax) - hc
  end if

  rt2sz = SQRT2*this%sigmaz
  erfcArg = (log(lnArg) - LNKAP + this%muz)/rt2sz
  liquid_saturation = this%Sr + (1.0d0-this%Sr)*5.0D-1*erfc(erfcArg)
  dsat_dpres = exp(-erfcArg**2)/(SQRTPI*rt2sz*lnArg)/UNIT_CONVERSION

end subroutine SFmKSaturation

! ************************************************************************** !

function RPFBurdineBCLiqCreate()

  ! Creates the Brooks-Corey Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_BC_liq_type), pointer :: RPFBurdineBCLiqCreate
  
  allocate(RPFBurdineBCLiqCreate)
  call RPFBurdineBCLiqCreate%Init()
  
end function RPFBurdineBCLiqCreate

! ************************************************************************** !

subroutine RPFBurdineBCLiqInit(this)

  ! Initializes the Brooks-Corey Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_BC_liq_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFBurdineBCLiqInit

! ************************************************************************** !

subroutine RPFBurdineBCLiqVerify(this,name,option)

  ! Initializes the Brooks-Corey Burdine relative permeability function object

  use Option_module
  
  implicit none
  
  class(rpf_Burdine_BC_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE'
  endif    
  call RPFBaseVerify(this,name,option)
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif
  
end subroutine RPFBurdineBCLiqVerify

! ************************************************************************** !

subroutine RPFBurdineBCLiqRelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Burdine_BC_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: power
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  ! reference #1
  power = 3.d0+2.d0/this%lambda
  relative_permeability = Se**power
  dkr_Se = power*relative_permeability/Se          
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPFBurdineBCLiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFIGHCC2CompLiqCreate()

  ! Creates the IGHCC2 Comparison relative permeability function object

  implicit none

  class(rpf_IGHCC2_Comp_liq_type), pointer :: &
                        RPFIGHCC2CompLiqCreate

  allocate(RPFIGHCC2CompLiqCreate)
  call RPFIGHCC2CompLiqCreate%Init()

end function RPFIGHCC2CompLiqCreate

! ************************************************************************** !

subroutine RPFIGHCC2CompLiqInit(this)

  ! Initializes the IGHCC2 Comparison relative permeability function object

  implicit none

  class(rpf_IGHCC2_Comp_liq_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE

  this%analytical_derivative_available = PETSC_TRUE

end subroutine RPFIGHCC2CompLiqInit

! ************************************************************************** !

subroutine RPFIGHCC2CompLiqVerify(this,name,option)

  ! Initializes the IGHCC2 Comparison relative permeability function object

  use Option_module

  implicit none

  class(rpf_IGHCC2_Comp_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,IGHCC2_COMP'
  endif
  call RPFBaseVerify(this,name,option)
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif

end subroutine RPFIGHCC2CompLiqVerify

! ************************************************************************** !

subroutine RPFIGHCC2CompLiqRelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation, to benchmark against IGHCC2 study.
  ! 
  ! Author: Michael Nole
  ! Date: 05/16/19
  ! 
  use Option_module

  implicit none

  class(rpf_IGHCC2_Comp_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option

  PetscReal :: Se
  PetscReal :: power
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0

  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif

  power = this%lambda
  relative_permeability = Se**power
  dkr_Se = power*relative_permeability/Se
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat

end subroutine RPFIGHCC2CompLiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFTableLiqCreate()

  ! Creates the Lookup Table relative permeability function object

  implicit none

  class(rpf_Table_liq_type), pointer :: &
                        RPFTableLiqCreate

  allocate(RPFTableLiqCreate)
  call RPFTableLiqCreate%Init()

end function RPFTableLiqCreate

! ************************************************************************** !

subroutine RPFTableLiqInit(this)

  implicit none

  class(rpf_Table_liq_type) :: this

  call RPFBaseInit(this)
  this%rpf_dataset => DatasetAsciiCreate()

  this%analytical_derivative_available = PETSC_TRUE

end subroutine RPFTableLiqInit

! ************************************************************************** !

subroutine RPFTableLiqVerify(this,name,option)

  use Option_module

  implicit none

  class(rpf_Table_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,LOOKUP_TABLE'
  endif
  call RPFBaseVerify(this,name,option)
  if (.not.associated(this%rpf_dataset)) then
    option%io_buffer = UninitializedMessage('TABLE',string)
    call PrintErrMsg(option)
  endif

end subroutine RPFTableLiqVerify

! ************************************************************************** !

subroutine RPFTableLiqRelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)

  use Option_module
  
  implicit none

  class(rpf_Table_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option

  class(dataset_ascii_type), pointer :: dataset
  PetscReal, pointer :: times(:)
  PetscInt :: i, j, num_entries

  dataset => this%rpf_dataset
  times => dataset%time_storage%times
  num_entries = 0
  j = 0
  do i = 1,size(times)
    if (times(i) <= dataset%time_storage%cur_time) then
      if (i > 1 .and. times(i) > times(i-1)) then
        j = 0
        num_entries = 0
      endif
      if (j==0) j = i
      num_entries = num_entries + 1
    endif
  enddo

  if (liquid_saturation <= this%sr) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0 !Not exactly true
    return
  endif

  if (liquid_saturation < dataset%rbuffer(2*j-1)) then
    relative_permeability = dataset%rbuffer(2*j)
    dkr_sat = 0.d0
  elseif (liquid_saturation > dataset%rbuffer(2*(j-1+num_entries)-1)) then
    dkr_sat = (relative_permeability - dataset%rbuffer(2*(j-1+num_entries))) / &
              (liquid_saturation - dataset%rbuffer(2*(j-1+num_entries)-1))
    relative_permeability = (liquid_saturation - dataset% &
                         rbuffer(2*(j-1+num_entries)-1)) * dkr_sat + & 
                         dataset%rbuffer(2*(j-1+num_entries))
  else
    do i = j+1, j+num_entries-1
      if (liquid_saturation <= dataset%rbuffer(2*i-1)) then
        dkr_sat = (dataset%rbuffer(2*i) - dataset%rbuffer(2*i-2)) / &
                  (dataset%rbuffer(2*i-1) - dataset%rbuffer(2*i-3))
        relative_permeability = (liquid_saturation - dataset%rbuffer(2*i-3)) * &
                         dkr_sat + dataset%rbuffer(2*i-2)
        exit
      endif
    enddo
  endif
  
  relative_permeability = maxval([relative_permeability, 0.d0])
  relative_permeability = minval([relative_permeability, 1.d0])

end subroutine RPFTableLiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFBurdineBCGasCreate()

  ! Creates the Brooks-Corey Burdine gas relative permeability function
  ! object

  implicit none
  
  class(rpf_Burdine_BC_gas_type), pointer :: RPFBurdineBCGasCreate
  
  allocate(RPFBurdineBCGasCreate)
  call RPFBurdineBCGasCreate%Init()
  
end function RPFBurdineBCGasCreate

! ************************************************************************** !

subroutine RPFBurdineBCGasInit(this)

  ! Initializes the Brooks-Corey Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Burdine_BC_gas_type) :: this

  call RPFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFBurdineBCGasInit

! ************************************************************************** !

subroutine RPFBurdineBCGasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Burdine_BC_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE_BC_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif  
  
end subroutine RPFBurdineBCGasVerify

! ************************************************************************** !

subroutine RPFBurdineBCGasRelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Burdine_BC_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  ! reference #1
  relative_permeability = Seg*Seg*(1.d0-Se**(1.d0+2.d0/this%lambda))
  ! Mathematica analytical solution (Heeho Park)
  dkr_Se = -(1.d0+2.d0/this%lambda)*Seg**2.d0*Se**(2.d0/this%lambda) &
           - 2.d0*Seg*(1.d0-Se**(1.d0+2.d0/this%lambda))
  dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPFBurdineBCGasRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFIGHCC2CompGasCreate()

  ! Creates the IGHCC2 Comparison gas relative permeability function
  ! object

  implicit none

  class(rpf_IGHCC2_Comp_gas_type), pointer :: &
                        RPFIGHCC2CompGasCreate

  allocate(RPFIGHCC2CompGasCreate)
  call RPFIGHCC2CompGasCreate%Init()

end function RPFIGHCC2CompGasCreate

! ************************************************************************** !

subroutine RPFIGHCC2CompGasInit(this)

  ! Initializes the IGHCC2 Comparison gas relative permeability function 
  ! object

  implicit none

  class(rpf_IGHCC2_Comp_gas_type) :: this

  call RPFBaseInit(this)

  this%analytical_derivative_available = PETSC_TRUE

end subroutine RPFIGHCC2CompGasInit

! ************************************************************************** !

subroutine RPFIGHCC2CompGasVerify(this,name,option)

  use Option_module

  implicit none

  class(rpf_IGHCC2_Comp_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,IGHCC2_Comp'
  endif
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif

end subroutine RPFIGHCC2CompGasVerify

! ************************************************************************** !
subroutine RPFIGHCC2CompGasRelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation, to benchmark against IGHCC2 study.
  ! 
  ! Author: Michael Nole
  ! Date: 05/16/19
  ! 
  use Option_module

  implicit none

  class(rpf_IGHCC2_Comp_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option

  PetscReal :: Se
  PetscReal :: power
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0

  Se = (1.d0 - liquid_saturation - this%Srg) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif

  power = this%lambda
  relative_permeability = Se**power
  dkr_Se = power*relative_permeability/Se
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat

end subroutine RPFIGHCC2CompGasRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFTableGasCreate()

  ! Creates the Lookup Table relative permeability function object

  implicit none

  class(rpf_Table_gas_type), pointer :: &
                        RPFTableGasCreate

  allocate(RPFTableGasCreate)
  call RPFTableGasCreate%Init()

end function RPFTableGasCreate

! ************************************************************************** !

subroutine RPFTableGasInit(this)

  implicit none

  class(rpf_Table_gas_type) :: this

  call RPFBaseInit(this)
  this%rpf_dataset => DatasetAsciiCreate()

  this%analytical_derivative_available = PETSC_TRUE

end subroutine RPFTableGasInit

! ************************************************************************** !

subroutine RPFTableGasVerify(this,name,option)

  use Option_module

  implicit none

  class(rpf_Table_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,LOOKUP_TABLE'
  endif
  call RPFBaseVerify(this,name,option)
  if (.not.associated(this%rpf_dataset)) then
    option%io_buffer = UninitializedMessage('TABLE',string)
    call PrintErrMsg(option)
  endif

end subroutine RPFTableGasVerify

! ************************************************************************** !

subroutine RPFTableGasRelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  use Option_module

  implicit none

  class(rpf_Table_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option

  class(dataset_ascii_type), pointer :: dataset
  PetscReal, pointer :: times(:)
  PetscInt :: i, j, num_entries

  dataset => this%rpf_dataset
  times => dataset%time_storage%times
  num_entries = 0
  j = 0
  do i = 1,size(times)
    if (times(i) <= dataset%time_storage%cur_time) then
      if (i > 1 .and. times(i) > times(i-1)) then
        j = 0
        num_entries = 0
      endif
      if (j==0) j = i
      num_entries = num_entries + 1
    endif
  enddo

  if (1.d0 - liquid_saturation <= this%srg) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0 !Not exactly true
    return
  endif

  if (liquid_saturation < dataset%rbuffer(2*j-1)) then
    relative_permeability = dataset%rbuffer(2*j)
    dkr_sat = 0.d0
  elseif (liquid_saturation > dataset%rbuffer(2*(j-1+num_entries)-1)) then
    dkr_sat = (relative_permeability - dataset%rbuffer(2*(j-1+num_entries))) / &
              (liquid_saturation - dataset%rbuffer(2*(j-1+num_entries)-1))
    relative_permeability = (liquid_saturation - dataset% &
                         rbuffer(2*(j-1+num_entries)-1)) * dkr_sat + &
                         dataset%rbuffer(2*(j-1+num_entries))
  else
    do i = j+1, j+num_entries-1
      if (liquid_saturation <= dataset%rbuffer(2*i-1)) then
        dkr_sat = (dataset%rbuffer(2*i) - dataset%rbuffer(2*i-2)) / &
                  (dataset%rbuffer(2*i-1) - dataset%rbuffer(2*i-3))
        relative_permeability = (liquid_saturation - dataset%rbuffer(2*i-3)) * &
                         dkr_sat + dataset%rbuffer(2*i-2)
        exit
      endif
    enddo
  endif

  relative_permeability = maxval([relative_permeability, 0.d0])
  relative_permeability = minval([relative_permeability, 1.d0])

end subroutine RPFTableGasRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFMualemBCLiqCreate()

  ! Creates the Brooks-Corey Mualem liquid relative permeability function object

  implicit none
  
  class(rpf_Mualem_BC_liq_type), pointer :: RPFMualemBCLiqCreate
  
  allocate(RPFMualemBCLiqCreate)
  call RPFMualemBCLiqCreate%Init()
  
end function RPFMualemBCLiqCreate

! ************************************************************************** !

subroutine RPFMualemBCLiqInit(this)

  ! Initializes the Brooks-Corey Mualem liquid relative permeability function 
  ! object

  implicit none
  
  class(rpf_Mualem_BC_liq_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE

end subroutine RPFMualemBCLiqInit

! ************************************************************************** !

subroutine RPFMualemBCLiqVerify(this,name,option)

  ! Initializes the Brooks-Corey Mualem liquid relative permeability function object

  use Option_module
  
  implicit none
  
  class(rpf_Mualem_BC_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM'
  endif    
  call RPFBaseVerify(this,name,option)
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif
  
end subroutine RPFMualemBCLiqVerify

! ************************************************************************** !

subroutine RPFMualemBCLiqRelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Mualem_BC_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: power
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  ! reference #1
  power = 2.5d0+2.d0/this%lambda
  relative_permeability = Se**power
  dkr_Se = power*relative_permeability/Se          
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat 

end subroutine RPFMualemBCLiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFMualemBCGasCreate()

  ! Creates the Brooks-Corey Mualem gas relative permeability function object

  implicit none
  
  class(rpf_Mualem_BC_gas_type), pointer :: RPFMualemBCGasCreate
  
  allocate(RPFMualemBCGasCreate)
  call RPFMualemBCGasCreate%Init()
  
end function RPFMualemBCGasCreate

! ************************************************************************** !

subroutine RPFMualemBCGasInit(this)

  ! Initializes the Brooks-Corey Mualem gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Mualem_BC_gas_type) :: this

  call RPFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFMualemBCGasInit

! ************************************************************************** !

subroutine RPFMualemBCGasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Mualem_BC_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM_BC_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif  
  
end subroutine RPFMualemBCGasVerify

! ************************************************************************** !

subroutine RPFMualemBCGasRelPerm(this,liquid_saturation, &
                                       relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Mualem_BC_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  ! reference Table 2
  relative_permeability = sqrt(Seg)* &
                             (1.d0-Se**(1.d0+1.d0/this%lambda))**2.d0
  ! Mathematica analytical solution (Heeho Park)
  dkr_Se = -2.d0*(1.d0+1.d0/this%lambda)*sqrt(Seg)*Se**(1.d0/this%lambda) &
          * (1.d0-Se**(1.d0+1.d0/this%lambda)) &
          - (1.d0-Se**(1.d0+1.d0/this%lambda))**2.d0/(2.d0*sqrt(Seg))
  dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPFMualemBCGasRelPerm

! ************************************************************************** !

function RPFMualemLinearLiqCreate()

  ! Creates the Linear Mualem relative permeability function object

  implicit none
  
  class(rpf_Mualem_linear_liq_type), pointer :: RPFMualemLinearLiqCreate
  
  allocate(RPFMualemLinearLiqCreate)
  call RPFMualemLinearLiqCreate%Init()
  
end function RPFMualemLinearLiqCreate

! ************************************************************************** !

subroutine RPFMualemLinearLiqInit(this)

  ! Initializes the Linear Mualem relative permeability function object

  implicit none
  
  class(rpf_Mualem_Linear_liq_type) :: this

  call RPFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%pcmax = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFMualemLinearLiqInit

! ************************************************************************** !

subroutine RPFMualemLinearLiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Mualem_Linear_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call PrintErrMsg(option)
  endif  
  if (Uninitialized(this%pcmax)) then
    option%io_buffer = UninitializedMessage('MAX_CAPILLARY_PRESSURE',string)
    call PrintErrMsg(option)
  endif
  
end subroutine RPFMualemLinearLiqVerify

! ************************************************************************** !

subroutine RPFMualemLinearLiqRelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  !   
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Mualem_Linear_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: one_over_alpha
  PetscReal :: pct_over_pcmax
  PetscReal :: pc_over_pcmax
  PetscReal :: pc_log_ratio
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  one_over_alpha = 1.d0/this%alpha
  pct_over_pcmax = one_over_alpha/this%pcmax
  pc_over_pcmax = 1.d0-(1.d0-pct_over_pcmax)*Se
  pc_log_ratio = log(pc_over_pcmax) / log(pct_over_pcmax)
  relative_permeability = (Se**0.5d0)*(pc_log_ratio**2.d0)
  ! ***used Mathematica to verify***
  ! In[3]:
  ! D[Se^(1/2)*(Log[1 - (1 - pctoverpcmax)*Se]/Log[pctoverpcmax])^2, Se]
  ! Out[3]:
  ! (2 (-1 + pctoverpcmax) Sqrt[Se]
  !  Log[1 - (1 - pctoverpcmax) Se])/((1 - (1 - pctoverpcmax) Se) Log[
  !  pctoverpcmax]^2) + Log[1 - (1 - pctoverpcmax) Se]^2/(
  ! 2 Sqrt[Se] Log[pctoverpcmax]^2)
  dkr_Se = 2.d0*(-1.d0+pct_over_pcmax)*sqrt(Se)* log(pc_over_pcmax) / &
    (pc_over_pcmax*log(pct_over_pcmax)**2.d0) + &
    log(pc_over_pcmax)**2.d0 / (2.d0*sqrt(Se)*log(pct_over_pcmax)**2.d0)
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat
         
end subroutine RPFMualemLinearLiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFMualemLinearGasCreate()

  ! Creates the Linear Mualem gas relative permeability function object

  implicit none
  
  class(rpf_Mualem_Linear_gas_type), pointer :: RPFMualemLinearGasCreate
  
  allocate(RPFMualemLinearGasCreate)
  call RPFMualemLinearGasCreate%Init()
  
end function RPFMualemLinearGasCreate
 
! ************************************************************************** !

subroutine RPFMualemLinearGasInit(this)

  ! Initializes the Linear Mualem gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Mualem_Linear_gas_type) :: this

  call RPFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%pcmax = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFMualemLinearGasInit

! ************************************************************************** !

subroutine RPFMualemLinearGasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Mualem_Linear_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM_LINEAR_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif  
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call PrintErrMsg(option)
  endif  
  if (Uninitialized(this%pcmax)) then
    option%io_buffer = UninitializedMessage('MAX_CAPILLARY_PRESSURE',string)
    call PrintErrMsg(option)
  endif
  
end subroutine RPFMualemLinearGasVerify

! ************************************************************************** !

subroutine RPFMualemLinearGasRelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14

  use Option_module
  
  implicit none

  class(rpf_Mualem_Linear_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: liquid_relative_permeability
  PetscReal :: liquid_dkr_sat  
  PetscReal :: dkr_dSe
  PetscReal :: dSe_dsat
  
  call RPFMualemLinearLiqRelPerm(this,liquid_saturation, &
                                     liquid_relative_permeability, &
                                     liquid_dkr_sat,option)
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0

  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  ! reference Table 2
  relative_permeability = Seg**0.5d0 * &
                 (1.d0-sqrt(liquid_relative_permeability*Se**(-0.5d0)))**2.d0
  ! Python analytical derivative (Jenn Frederick)
  dkr_dSe = 0.5d0*1.d0/Se*sqrt(Se**(-0.5d0)*liquid_relative_permeability)* &
    sqrt(1.d0-Se)*(1.d0-sqrt(Se**(-0.5d0)*liquid_relative_permeability))**1.0 &
    - (1.d0-sqrt(Se**(-0.5d0)*liquid_relative_permeability))**2.d0 &
    /(2.d0*sqrt(1.d0-Se))
  !one_over_apcm = 1.d0/(1.d-7)/(1.d9)
  !dkr_dSe = -2.0*Se**0.5*sqrt(Se**(-0.5)*log(-Se*(-one_over_apcm + 1.0) + 1.0)/log(one_over_apcm))*sqrt(-Se + 1.0)*(-0.25*Se**(-1.5)*log(-Se*(-one_over_apcm + 1.0) + 1.0)/log(one_over_apcm) + Se**(-0.5)*(one_over_apcm - 1.0)/(2*(-Se*(-one_over_apcm + 1.0) + 1.0)*log(one_over_apcm)))*(-sqrt(Se**(-0.5)*log(-Se*(-one_over_apcm + 1.0) + 1.0)/log(one_over_apcm)) + 1.0)**1.0*log(one_over_apcm)/log(-Se*(-one_over_apcm + 1.0) + 1.0) - (-sqrt(Se**(-0.5)*log(-Se*(-one_over_apcm + 1.0) + 1.0)/log(one_over_apcm)) + 1.0)**2.0/(2*sqrt(-Se + 1.0))
  dSe_dsat = 1.d0/(1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_dSe*dSe_dsat

end subroutine RPFMualemLinearGasRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFBurdineLinearLiqCreate()

  ! Creates the Linear Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_Linear_liq_type), pointer :: RPFBurdineLinearLiqCreate
  
  allocate(RPFBurdineLinearLiqCreate)
  call RPFBurdineLinearLiqCreate%Init()
  
end function RPFBurdineLinearLiqCreate

! ************************************************************************** !

subroutine RPFBurdineLinearLiqInit(this)

  ! Initializes the Linear Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_Linear_liq_type) :: this

  call RPFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFBurdineLinearLiqInit

! ************************************************************************** !

subroutine RPFBurdineLinearLiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Burdine_Linear_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE'
  endif  
  call RPFBaseVerify(this,string,option)
  
end subroutine RPFBurdineLinearLiqVerify

! ************************************************************************** !

subroutine RPFBurdineLinearLiqRelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !  
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  ! 
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Burdine_Linear_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  relative_permeability = Se
  dkr_sat = 1.d0 / (1.d0 - this%Sr)
  
end subroutine RPFBurdineLinearLiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFBurdineLinearGasCreate()

  ! Creates the Linear Burdine gas relative permeability function object

  implicit none
  
  class(rpf_Burdine_Linear_gas_type), pointer :: RPFBurdineLinearGasCreate
  
  allocate(RPFBurdineLinearGasCreate)
  call RPFBurdineLinearGasCreate%Init()
  
end function RPFBurdineLinearGasCreate

! ************************************************************************** !

subroutine RPFBurdineLinearGasInit(this)

  ! Initializes the Linear Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Burdine_Linear_gas_type) :: this

  call RPFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFBurdineLinearGasInit

! ************************************************************************** !

subroutine RPFBurdineLinearGasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Burdine_Linear_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE_LINEAR_GAS&
             &/BRAGFLO_ KRP5'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif  
  
end subroutine RPFBurdineLinearGasVerify

! ************************************************************************** !

subroutine RPFBurdineLinearGasRelPerm(this,liquid_saturation, &
                                          relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  !

  use Option_module
  
  implicit none

  class(rpf_Burdine_Linear_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  relative_permeability = Seg
  dkr_Se = -1.d0
  dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPFBurdineLinearGasRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFmKLiqCreate()

  ! Creates the modified Kosugi liq relative permeability function object

  implicit none

  class(rpf_mK_liq_type), pointer :: RPFmKLiqCreate

  allocate(RPFmKLiqCreate)
  call RPFmKLiqCreate%Init()

end function RPFmKLiqCreate

! ************************************************************************** !

subroutine RPFmKLiqInit(this)

  ! Initializes modified Kosugi saturation function object

  implicit none
  
  class(rpf_mK_liq_type) :: this

  call RPFBaseInit(this)
  this%sigmaz = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFmKLiqInit

! ************************************************************************** !

subroutine RPFmKLiqVerify(this,name,option)

  use Option_module

  implicit none

  class(rpf_mK_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'LIQUID_RELATIVE_PERM') > 0) then
    string = name
  else
    string = trim(name) // 'LIQUID_RELATIVE_PERM,MODIFIED_KOSUGI'
  endif
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%sigmaz)) then
    option%io_buffer = UninitializedMessage('SIGMAZ',string)
    call PrintErrMsg(option)
  endif

end subroutine RPFmKLiqVerify

! ************************************************************************** !

subroutine RPFmKLiqRelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  !
  ! Computes the relative permeability (and associated derivatives) as a
  ! function of saturation for modified Kosugi model
  !
  ! Malama, B. & K.L. Kuhlman, 2015. Unsaturated Hydraulic Conductivity
  ! Models Based on Truncated Lognormal Pore-size Distributions, Groundwater,
  ! 53(3):498�502. http://dx.doi.org/10.1111/gwat.12220
  !
  ! Author: Kris Kuhlman
  ! Date: 2017
  !
  use Option_module
  use Utility_module

  implicit none

  ! gnu & intel extension and required in f2008
  intrinsic :: erfc

  PetscReal, parameter :: SQRT2 = sqrt(2.0d0)

  class(rpf_mK_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option

  PetscReal :: Se, dkr_Se
  PetscReal :: InvSatRange
  PetscReal :: erfcArg, erfcRes
  PetscReal :: invErfcRes
  PetscReal :: sqrtSe, expArg
  PetscReal :: dinvErfcRes_dSe

  relative_permeability = 0.d0
  dkr_sat = 0.d0

  InvSatRange = 1.0d0/(1.0d0 - this%Sr)
  Se = (liquid_saturation - this%Sr)*InvSatRange
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif

!  invErfcRes = InverseNorm(Se)
  call InverseNorm(Se,invErfcRes,PETSC_TRUE,dinvErfcRes_dSe)
  erfcArg = (this%sigmaz - invErfcRes)/SQRT2
  erfcRes = erfc(erfcArg)
  sqrtSe = sqrt(Se)
  relative_permeability = sqrtSe*erfcRes*5.0D-1

  ! from Wolfram Alpha (x -> Se)
  ! (InverseErfc[x] -> -1/Sqrt[x] InverseNorm[x/2])
  !
  ! D[(Sqrt[x] Erfc[sigmaz/Sqrt[2] + InverseErfc[2 x]])/2, x] =
  ! E^(InverseErfc[2 x]^2 - (simgaz/Sqrt[2] + InverseErfc[2 x])^2) * ...
  ! Sqrt[x] + Erfc[sigmaz/Sqrt[2] + InverseErfc[2 x]]/(4 Sqrt[x])
  expArg = 5.0D-1*invErfcRes**2 - erfcArg**2
  dkr_Se = erfcres/(4.0D0*sqrtSe) + sqrtSe*exp(expArg)

  ! InvSatRange = dSe/dsat
  dkr_sat = dkr_Se * InvSatRange 

end subroutine RPFmKLiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFmKGasCreate()

  ! Creates the modified Kosugi gas relative permeability function object

  implicit none

  class(rpf_mK_gas_type), pointer :: RPFmKGasCreate

  allocate(RPFmKGasCreate)
  call RPFmKGasCreate%Init()

end function RPFmKGasCreate

! ************************************************************************** !

subroutine RPFmKGasInit(this)

  ! Initializes modified Kosugi saturation function object

  implicit none
  
  class(rpf_mK_gas_type) :: this

  call RPFBaseInit(this)
  this%sigmaz = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFmKGasInit

! ************************************************************************** !

subroutine RPFmKGasVerify(this,name,option)

  use Option_module

  implicit none

  class(rpf_mK_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'GAS_RELATIVE_PERM') > 0) then
    string = name
  else
    string = trim(name) // 'GAS_RELATIVE_PERM,MODIFIED_KOSUGI'
  endif
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%sigmaz)) then
    option%io_buffer = UninitializedMessage('SIGMAZ',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%srg)) then
    option%io_buffer = UninitializedMessage('SRG',string)
    call PrintErrMsg(option)
  endif

end subroutine RPFmKGasVerify

! ************************************************************************** !

subroutine RPFmKGasRelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  !
  ! Computes the relative permeability (and associated derivatives) as a
  ! function of saturation for modified Kosugi model
  !
  ! Malama, B. & K.L. Kuhlman, 2015. Unsaturated Hydraulic Conductivity
  ! Models Based on Truncated Lognormal Pore-size Distributions, Groundwater,
  ! 53(3):498�502. http://dx.doi.org/10.1111/gwat.12220
  !
  ! Author: Kris Kuhlman
  ! Date: 2017
  !
  use Option_module
  use Utility_module, only : InverseNorm

  implicit none

  ! gnu & intel extension and required in f2008
  intrinsic :: erfc

  PetscReal, parameter :: SQRT2 = sqrt(2.0d0)

  class(rpf_mK_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option

  PetscReal :: Se, Seg, InvSatRange
  PetscReal :: dkr_Se
  PetscReal :: erfcArg, erfcRes
  PetscReal :: invErfcRes
  PetscReal :: sqrtSe, expArg
  PetscReal :: dinvErfcRes_dSeg

  InvSatRange = 1.d0/(1.d0 - this%Sr - this%Srg)
  Se = (liquid_saturation - this%Sr)*InvSatRange

  relative_permeability = 0.d0
  dkr_sat = 0.d0
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif

  Seg = 1.d0 - Se

!  invErfcRes = InverseNorm(Seg)
  call InverseNorm(Seg,invErfcRes,PETSC_TRUE,dinvErfcRes_dSeg)
  erfcArg = (this%sigmaz - invErfcRes)/SQRT2
  erfcRes = erfc(erfcArg)
  sqrtSe = sqrt(Seg)
  relative_permeability = sqrtSe*erfcRes*5.0D-1

  ! from Wolfram Alpha (x -> Seg)
  ! (InverseErfc[x] -> -1/Sqrt[x] InverseNorm[x/2])
  !
  ! D[(Sqrt[x] Erfc[sigmaz/Sqrt[2] + InverseErfc[2 x]])/2, x] =
  ! E^(InverseErfc[2 x]^2 - (simgaz/Sqrt[2] + InverseErfc[2 x])^2) * ...
  ! Sqrt[x] + Erfc[sigmaz/Sqrt[2] + InverseErfc[2 x]]/(4 Sqrt[x])
  expArg = 5.0D-1*invErfcRes**2 - erfcArg**2
  dkr_Se = erfcres/(4.0D0*sqrtSe) + sqrtSe*exp(expArg)

  ! -1 = dSeg/dSe
  ! InvSatRange = dSe/dsat
  dkr_sat = -1.d0 * dkr_Se * InvSatRange 

end subroutine RPFmKGasRelPerm

 
end module Characteristic_Curves_Common_module
