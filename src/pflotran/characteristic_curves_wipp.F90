module Characteristic_Curves_WIPP_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_Common_module

  implicit none

  private
  
!-----------------------------------------------------------------------------
!-- Saturation Functions -----------------------------------------------------
!-----------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_WIPP_type
    PetscInt :: kpc
    PetscReal :: pct_a
    PetscReal :: pct_exp
    PetscReal :: pct
    PetscReal :: alpha
    PetscBool :: ignore_permeability
  contains
    procedure, public :: Init => SFWIPPInit
    procedure, public :: Verify => SFWIPPVerify
    procedure, public :: CapillaryPressure => SFWIPPCapillaryPressure
    procedure, public :: Saturation => SFWIPPSaturation
  end type sat_func_WIPP_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP1_type
    PetscReal :: Srg
    PetscReal :: m
  contains
    procedure, public :: Init => SFKRP1Init
    procedure, public :: Verify => SFKRP1Verify
    procedure, public :: CapillaryPressure => SFKRP1CapillaryPressure
    procedure, public :: Saturation => SFKRP1Saturation
  end type sat_func_KRP1_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP2_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => SFKRP2Init
    procedure, public :: Verify => SFKRP2Verify
    procedure, public :: CapillaryPressure => SFKRP2CapillaryPressure
    procedure, public :: Saturation => SFKRP2Saturation
  end type sat_func_KRP2_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP3_type
    PetscReal :: Srg
    PetscReal :: lambda
  contains
    procedure, public :: Init => SFKRP3Init
    procedure, public :: Verify => SFKRP3Verify
    procedure, public :: CapillaryPressure => SFKRP3CapillaryPressure
    procedure, public :: Saturation => SFKRP3Saturation
  end type sat_func_KRP3_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP4_type
    PetscReal :: Srg
    PetscReal :: lambda
  contains
    procedure, public :: Init => SFKRP4Init
    procedure, public :: Verify => SFKRP4Verify
    procedure, public :: CapillaryPressure => SFKRP4CapillaryPressure
    procedure, public :: Saturation => SFKRP4Saturation
  end type sat_func_KRP4_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP5_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => SFKRP5Init
    procedure, public :: Verify => SFKRP5Verify
    procedure, public :: CapillaryPressure => SFKRP5CapillaryPressure
    procedure, public :: Saturation => SFKRP5Saturation
  end type sat_func_KRP5_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP8_type
    PetscReal :: Srg
    PetscReal :: m
  contains
    procedure, public :: Init => SFKRP8Init
    procedure, public :: Verify => SFKRP8Verify
    procedure, public :: CapillaryPressure => SFKRP8CapillaryPressure
    procedure, public :: Saturation => SFKRP8Saturation
  end type sat_func_KRP8_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_KRP9_type
  contains
    procedure, public :: Init => SFKRP9Init
    procedure, public :: Verify => SFKRP9Verify
    procedure, public :: CapillaryPressure => SFKRP9CapillaryPressure
    procedure, public :: Saturation => SFKRP9Saturation
  end type sat_func_KRP9_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_KRP11_type
  contains
    procedure, public :: Init => SFKRP11Init
    procedure, public :: Verify => SFKRP11Verify
    procedure, public :: CapillaryPressure => SFKRP11CapillaryPressure
    procedure, public :: Saturation => SFKRP11Saturation
  end type sat_func_KRP11_type 
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP12_type
    PetscReal :: lambda
    PetscReal :: s_min
    PetscReal :: s_effmin
  contains
    procedure, public :: Init => SFKRP12Init
    procedure, public :: Verify => SFKRP12Verify
    procedure, public :: CapillaryPressure => SFKRP12CapillaryPressure
    procedure, public :: Saturation => SFKRP12Saturation
  end type sat_func_KRP12_type 
  
!-----------------------------------------------------------------------------
!-- Relative Permeability Functions ------------------------------------------
!-----------------------------------------------------------------------------  
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP1_liq_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPFKRP1LiqInit
    procedure, public :: Verify => RPFKRP1LiqVerify
    procedure, public :: RelativePermeability => RPFKRP1LiqRelPerm
  end type rpf_KRP1_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP1_gas_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPFKRP1GasInit
    procedure, public :: Verify => RPFKRP1GasVerify
    procedure, public :: RelativePermeability => RPFKRP1GasRelPerm
  end type rpf_KRP1_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP2_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPFKRP2LiqInit
    procedure, public :: Verify => RPFKRP2LiqVerify
    procedure, public :: RelativePermeability => RPFKRP2LiqRelPerm
  end type rpf_KRP2_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP2_gas_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPFKRP2GasInit
    procedure, public :: Verify => RPFKRP2GasVerify
    procedure, public :: RelativePermeability => RPFKRP2GasRelPerm
  end type rpf_KRP2_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP3_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPFKRP3LiqInit
    procedure, public :: Verify => RPFKRP3LiqVerify
    procedure, public :: RelativePermeability => RPFKRP3LiqRelPerm
  end type rpf_KRP3_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP3_gas_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPFKRP3GasInit
    procedure, public :: Verify => RPFKRP3GasVerify
    procedure, public :: RelativePermeability => RPFKRP3GasRelPerm
  end type rpf_KRP3_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP3_liq_type) :: rpf_KRP4_liq_type
  contains
    procedure, public :: Verify => RPFKRP4LiqVerify
    procedure, public :: RelativePermeability => RPFKRP4LiqRelPerm
  end type rpf_KRP4_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP3_gas_type) :: rpf_KRP4_gas_type
  contains
    procedure, public :: Verify => RPFKRP4GasVerify
  end type rpf_KRP4_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP5_liq_type
  contains
    procedure, public :: Init => RPFKRP5LiqInit
    procedure, public :: Verify => RPFKRP5LiqVerify
    procedure, public :: RelativePermeability => RPFKRP5LiqRelPerm
  end type rpf_KRP5_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP5_gas_type
  contains
    procedure, public :: Init => RPFKRP5GasInit
    procedure, public :: Verify => RPFKRP5GasVerify
    procedure, public :: RelativePermeability => RPFKRP5GasRelPerm
  end type rpf_KRP5_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP8_liq_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPFKRP8LiqInit
    procedure, public :: Verify => RPFKRP8LiqVerify
    procedure, public :: RelativePermeability => RPFKRP8LiqRelPerm
  end type rpf_KRP8_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP8_gas_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPFKRP8GasInit
    procedure, public :: Verify => RPFKRP8GasVerify
    procedure, public :: RelativePermeability => RPFKRP8GasRelPerm
  end type rpf_KRP8_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP9_liq_type
  contains
    procedure, public :: Init => RPFKRP9LiqInit
    procedure, public :: Verify => RPFKRP9LiqVerify
    procedure, public :: RelativePermeability => RPFKRP9LiqRelPerm
  end type rpf_KRP9_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP9_liq_type) :: rpf_KRP9_gas_type
  contains
    procedure, public :: Init => RPFKRP9GasInit
    procedure, public :: Verify => RPFKRP9GasVerify
    procedure, public :: RelativePermeability => RPFKRP9GasRelPerm
  end type rpf_KRP9_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP11_liq_type
    PetscReal :: tolc
  contains
    procedure, public :: Init => RPFKRP11LiqInit
    procedure, public :: Verify => RPFKRP11LiqVerify
    procedure, public :: RelativePermeability => RPFKRP11LiqRelPerm
  end type rpf_KRP11_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP11_liq_type) :: rpf_KRP11_gas_type
  contains
    procedure, public :: Verify => RPFKRP11GasVerify
    procedure, public :: RelativePermeability => RPFKRP11GasRelPerm
  end type rpf_KRP11_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP4_liq_type) :: rpf_KRP12_liq_type
  contains
    procedure, public :: Verify => RPFKRP12LiqVerify
    procedure, public :: RelativePermeability => RPFKRP12LiqRelPerm
  end type rpf_KRP12_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP3_gas_type) :: rpf_KRP12_gas_type
  contains
    procedure, public :: Verify => RPFKRP12GasVerify
    procedure, public :: RelativePermeability => RPFKRP12GasRelPerm
  end type rpf_KRP12_gas_type
  !---------------------------------------------------------------------------
  ! since the TOUGH2_Corey relative permeability function (IRP=7 in 
  ! TOUGH2 manual) calculates relative perm as a function of the 
  ! Mualem-based  liquid relative permeability when Srg = 0., we extend 
  ! the rpf_Mualem_type to save code
  type, public, extends(rpf_Mualem_VG_liq_type) :: rpf_TOUGH2_IRP7_gas_type
  contains
    procedure, public :: Init => RPFTOUGH2IRP7GasInit
    procedure, public :: Verify => RPFTOUGH2IRP7GasVerify
    procedure, public :: RelativePermeability => RPFTOUGH2IRP7GasRelPerm
  end type rpf_TOUGH2_IRP7_gas_type
  
  public :: &! WIPP saturation functions:
            SFKRP1Create, &
            SFKRP2Create, &
            SFKRP3Create, &
            SFKRP4Create, &
            SFKRP5Create, &
            SFKRP8Create, &
            SFKRP9Create, &
            SFKRP11Create, &
            SFKRP12Create, &
            ! WIPP rel. perm. curves:
            RPFKRP1LiqCreate, &
            RPFKRP1GasCreate, &
            RPFKRP2LiqCreate, &
            RPFKRP2GasCreate, &
            RPFKRP3LiqCreate, &
            RPFKRP3GasCreate, &
            RPFKRP4LiqCreate, &
            RPFKRP4GasCreate, &
            RPFKRP5LiqCreate, &
            RPFKRP5GasCreate, &
            RPFKRP8LiqCreate, &
            RPFKRP8GasCreate, &
            RPFKRP9LiqCreate, &
            RPFKRP9GasCreate, &
            RPFKRP11LiqCreate, &
            RPFKRP11GasCreate, &
            RPFKRP12LiqCreate, &
            RPFKRP12GasCreate, &
            RPFTOUGH2IRP7GasCreate
  
contains

! ************************************************************************** !

subroutine SFWIPPInit(this)

  ! Initializes a sat_func_WIPP_type object.

  implicit none
  
  class(sat_func_WIPP_type) :: this

  call SFBaseInit(this)
  this%kpc = UNINITIALIZED_INTEGER
  this%pct_a = UNINITIALIZED_DOUBLE
  this%pct_exp = UNINITIALIZED_DOUBLE
  this%ignore_permeability = PETSC_FALSE
  this%alpha = UNINITIALIZED_DOUBLE
  
end subroutine SFWIPPInit

! ************************************************************************** !

subroutine SFWIPPVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_WIPP_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: num_errors
  
  num_errors = 0
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION'
  endif
  call SFBaseVerify(this,string,option)
  
  if (.not.this%ignore_permeability) then
    if (Uninitialized(this%pct_a)) then
      option%io_buffer = UninitializedMessage('PCT_A',string)
      call PrintMsg(option)
      num_errors = num_errors + 1
    endif   
    if (Uninitialized(this%pct_exp)) then
      option%io_buffer = UninitializedMessage('PCT_EXP',string)
      call PrintMsg(option)
      num_errors = num_errors + 1
    endif 
  else
    if (Uninitialized(this%alpha)) then
      option%io_buffer = UninitializedMessage('ALPHA',string)
      call PrintMsg(option)
      num_errors = num_errors + 1
    endif   
  endif
  if (Uninitialized(this%kpc)) then
    option%io_buffer = UninitializedMessage('KPC',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif 
   
  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the ' // trim(string) // ' block. See above.'
    call PrintErrMsg(option)
  endif

end subroutine SFWIPPVerify

! ************************************************************************** !

subroutine SFWIPPCapillaryPressure(this,liquid_saturation, & 
                                     capillary_pressure,dpc_dsatl,option)
  use Option_module
  
  implicit none
  
  class(sat_func_WIPP_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFWIPPCapillaryPressure must be extended.'
  call PrintErrMsg(option)
  
end subroutine SFWIPPCapillaryPressure

! ************************************************************************** !

subroutine SFWIPPSaturation(this,capillary_pressure, &
                              liquid_saturation,dsat_dpres,option)
  use Option_module

  implicit none
  
  class(sat_func_WIPP_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFWIPPSaturation must be extended.'
  call PrintErrMsg(option)
  
end subroutine SFWIPPSaturation

! ************************************************************************** !

subroutine SFWIPPKPC(this,lambda,PT,Se,capillary_pressure)
  !
  ! Calculates the SEMIN value that is used for truncation of capillary 
  ! pressure in the CapillaryPressure functions if KPC=2.
  !
  ! Author: Jennifer Frederick
  ! Date: 06/12/2017
  !
  
  implicit none
  
  class(sat_func_WIPP_type) :: this
  PetscReal :: lambda
  PetscReal :: PT
  PetscReal :: Se
  PetscReal :: capillary_pressure
  
  PetscReal :: SEMIN
  PetscBool :: BC     ! Brooks-Corey type
  PetscBool :: VG     ! van Genuchten type
  
  BC = PETSC_FALSE
  VG = PETSC_FALSE
  
  if (this%kpc /= 2) return
  
  ! calculate SEMIN:
  if (PT < 0.d0) then
    SEMIN = -1.d0
  else
    select type(this)
      type is(sat_func_KRP1_type)
        VG = PETSC_TRUE
      type is(sat_func_KRP2_type)
        BC = PETSC_TRUE
      type is(sat_func_KRP3_type)
        BC = PETSC_TRUE
      type is(sat_func_KRP4_type)
        BC = PETSC_TRUE
      type is(sat_func_KRP8_type)
        VG = PETSC_TRUE
      type is(sat_func_KRP12_type)
        SEMIN = 0.d0
        ! BC = PETSC_TRUE
      class default
        return
    end select
    if (BC) then
      SEMIN = (PT/this%pcmax)**lambda
    else if (VG) then
      SEMIN = (1.d0+((this%pcmax/PT)**(lambda+1.d0))) &
               **(-1.d0*lambda/(1.d0+lambda))
    endif
  endif
  
  if (Se <= SEMIN) then
    capillary_pressure = this%pcmax
  endif
  
end subroutine SFWIPPKPC

! ************************************************************************** !
! ************************************************************************** !

function SFKRP1Create()

  ! Creates the BRAGFLO KRP1 capillary pressure function object

  implicit none
  
  class(sat_func_KRP1_type), pointer :: SFKRP1Create
  
  allocate(SFKRP1Create)
  call SFKRP1Create%Init()
  
end function SFKRP1Create

! ************************************************************************** !

subroutine SFKRP1Init(this)

  ! Creates the BRAGFLO KRP1 capillary pressure function object

  implicit none
  
  class(sat_func_KRP1_type) :: this

  call SFBaseInit(this)
  call SFWIPPInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%m = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SFKRP1Init

! ************************************************************************** !

subroutine SFKRP1Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_KRP1_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: num_errors
  
  num_errors = 0
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP1'
  endif
  call SFBaseVerify(this,string,option)
  call SFWIPPVerify(this,string,option)
  
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif   
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif   
  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the ' // trim(string) // ' block. See above.'
    call PrintErrMsg(option)
  endif

end subroutine SFKRP1Verify

! ************************************************************************** !

subroutine SFKRP1CapillaryPressure(this,liquid_saturation, &
                                     capillary_pressure,dpc_dsatl,option)     
  ! 
  ! Computes the capillary pressure as a function of saturation using the
  ! modified van Genuchten-Parker formulation.
  ! BRAGFLO UM 6.02 pg 41, 42; eq.120
  ! Modified according to KRP=1 option of BRAGFLO:
  ! Effective saturation includes residual gas saturation.
  ! Threshold pressure is a function of permeability.
  ! P0 parameter.
  !
  ! Author: Heeho Park; Modified by Jennifer Frederick
  ! Date: 11/17/16; Modified 04/26/2017
  !
  use Option_module
  
  implicit none
  
  class(sat_func_KRP1_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: lambda
  PetscReal :: Se2
  PetscReal :: P0
  
  dpc_dsatl = 0.d0
  dpc_dsatl = 1.d0/dpc_dsatl
  dpc_dsatl = 0.d0*dpc_dsatl
   
  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP1_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  lambda = this%m/(1.d0-this%m)
  ! Derivation for P0 is in Appendix PA:
  ! It is derived by setting Se2 in KRP4 and KRP1 to the value 0.5, and then 
  ! equating Pc(KRP1,Se2=0.5) = Pc(KRP4,Se2=0.5) and solving for P0.
  P0 = this%pct * (2.d0**(1.d0/lambda)) * &
       (((0.5d0**(-1.d0/this%m))-1.d0)**(this%m-1.d0))
  Se2 = (liquid_saturation-this%Sr)/(1.d0-this%Sr-this%Srg)
  Se2 = min(Se2,1.d0)
  
  if ((liquid_saturation > this%Sr) .or. &
      ((1.d0 - liquid_saturation) <= this%Srg)) then
    capillary_pressure = P0*(Se2**(-1.d0/this%m)-1.d0)**(1.d0-this%m)
  else
    capillary_pressure = 0.d0
  endif

#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
  endif
#endif

  call SFWIPPKPC(this,lambda,P0,Se2,capillary_pressure)

end subroutine SFKRP1CapillaryPressure

! ************************************************************************** !

subroutine SFKRP1Saturation(this,capillary_pressure, &
                              liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the liquid saturation as a function of capillary pressure using
  ! the modified van Genuchten-Parker formulation.
  ! BRAGFLO UM 6.02 pg 41, 42; eq.120
  ! Modified according to KRP=1 option of BRAGFLO:
  ! Effective saturation includes residual gas saturation.
  ! Threshold pressure is a function of permeability.
  ! P0 parameter.
  !
  ! Author: Heeho Park; Modified by Jenn Frederick
  ! Date: 11/17/16; Modified on 05/03/2017
  !
  use Option_module
  
  implicit none

  class(sat_func_KRP1_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: lambda
  PetscReal :: Se2
  PetscReal :: P0
  PetscReal :: n
  
  dsat_dpres = 0.d0
  dsat_dpres = 1.d0/dsat_dpres
  dsat_dpres = 0.d0*dsat_dpres
  
  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP1_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  lambda = this%m/(1.d0-this%m)
  ! Derivation for P0 is in Appendix PA:
  ! It is derived by setting Se2 in KRP4 and KRP1 to the value 0.5, and then 
  ! equating Pc(KRP1,Se2=0.5) = Pc(KRP4,Se2=0.5) and solving for P0.
  P0 = this%pct * (2.d0**(1.d0/lambda)) * &
       (((0.5d0**(-1.d0/this%m))-1.d0)**(this%m-1.d0))
    
  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  else
    n = 1.d0/(1.d0-this%m)
    Se2 = (((capillary_pressure**n)/P0) + 1.d0)**(-1.d0*this%m)
    liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se2
  endif
  
end subroutine SFKRP1Saturation

! ************************************************************************** !
! ************************************************************************** !

function SFKRP2Create()

  ! Creates the BRAGFLO KRP2 capillary pressure function object

  implicit none
  
  class(sat_func_KRP2_type), pointer :: SFKRP2Create
  
  allocate(SFKRP2Create)
  call SFKRP2Create%Init()
  
end function SFKRP2Create

! ************************************************************************** !

subroutine SFKRP2Init(this)

  ! Creates the BRAGFLO KRP2 capillary pressure function object

  implicit none
  
  class(sat_func_KRP2_type) :: this

  call SFBaseInit(this)
  call SFWIPPInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SFKRP2Init

! ************************************************************************** !

subroutine SFKRP2Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_KRP2_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP2'
  endif
  call SFBaseVerify(this,string,option)
  call SFWIPPVerify(this,string,option)
  
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif    

end subroutine SFKRP2Verify

! ************************************************************************** !

subroutine SFKRP2CapillaryPressure(this,liquid_saturation, &
                                     capillary_pressure,dpc_dsatl,option)     
  ! 
  ! Computes the capillary pressure as a function of saturation using the
  ! original Brooks Corey formulation.
  ! BRAGFLO UM 6.02 pg 41, 42; eq.126
  ! Modified according to KRP=2 option of BRAGFLO:
  ! Threshold pressure is a function of permeability.
  !
  ! Author: Jennifer Frederick
  ! Date: 05/02/2017
  !
  use Option_module
  
  implicit none
  
  class(sat_func_KRP2_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1

  dpc_dsatl = 0.d0
  dpc_dsatl = 1.d0/dpc_dsatl
  dpc_dsatl = 0.d0*dpc_dsatl
  
  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP2_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  Se1 = (liquid_saturation-this%Sr)/(1.d0-this%Sr)
      
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = 0.d0
  else 
    capillary_pressure = this%pct/(Se1**(1.d0/this%lambda))
  endif
  
  call SFWIPPKPC(this,this%lambda,this%pct,Se1,capillary_pressure)
  
end subroutine SFKRP2CapillaryPressure

! ************************************************************************** !

subroutine SFKRP2Saturation(this,capillary_pressure, &
                              liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the liquid saturation as a function of capillary pressure using 
  ! the original Brooks Corey formulation.
  ! BRAGFLO UM 6.02 pg 41, 42; eq.126
  ! Modified according to KRP=2 option of BRAGFLO:
  ! Threshold pressure is a function of permeability.
  !
  ! Author: Jennifer Frederick
  ! Date: 05/02/2017
  !
  use Option_module
  
  implicit none

  class(sat_func_KRP2_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1     
  
  dsat_dpres = 0.d0
  dsat_dpres = 1.d0/dsat_dpres
  dsat_dpres = 0.d0*dsat_dpres

  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP2_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  if (capillary_pressure < this%pct) then
    liquid_saturation = 1.d0
    dsat_dpres = 0.d0
    return
  endif
  
  Se1 = (capillary_pressure/this%pct)**(-1.d0*this%lambda)
  liquid_saturation = this%Sr + (1.d0-this%Sr)*Se1
  
end subroutine SFKRP2Saturation

! ************************************************************************** !
! ************************************************************************** !

function SFKRP3Create()

  ! Creates the BRAGFLO KRP3 capillary pressure function object

  implicit none
  
  class(sat_func_KRP3_type), pointer :: SFKRP3Create
  
  allocate(SFKRP3Create)
  call SFKRP3Create%Init()
  
end function SFKRP3Create

! ************************************************************************** !

subroutine SFKRP3Init(this)

  ! Creates the BRAGFLO KRP3 capillary pressure function object

  implicit none
  
  class(sat_func_KRP3_type) :: this

  call SFBaseInit(this)
  call SFWIPPInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SFKRP3Init

! ************************************************************************** !

subroutine SFKRP3Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_KRP3_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: num_errors

  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP3'
  endif
  call SFBaseVerify(this,string,option)
  call SFWIPPVerify(this,string,option)
  
  num_errors = 0
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif   
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif   
  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the ' // trim(string) // ' block. See above.'
    call PrintErrMsg(option)
  endif   

end subroutine SFKRP3Verify

! ************************************************************************** !

subroutine SFKRP3CapillaryPressure(this,liquid_saturation, &
                                     capillary_pressure,dpc_dsatl,option)     
  ! 
  ! Computes the capillary pressure as a function of saturation using the
  ! modified Brooks Corey formulation.
  ! BRAGFLO UM 6.02 pg 41, 42; eq.126, eq. 129
  ! Modified according to KRP=3 option of BRAGFLO:
  ! Threshold pressure is a function of permeability.
  ! Effective saturation is a function of residual gas saturation.
  !
  ! Author: Jennifer Frederick
  ! Date: 05/04/2017
  !
  use Option_module
  
  implicit none
  
  class(sat_func_KRP3_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  
  dpc_dsatl = 0.d0
  dpc_dsatl = 1.d0/dpc_dsatl
  dpc_dsatl = 0.d0*dpc_dsatl
  
  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP3_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  Se2 = (liquid_saturation-this%Sr)/(1.d0-this%Sr-this%Srg)
  
  if ((1.d0-liquid_saturation) <= this%Srg) then
    capillary_pressure = this%pct
  elseif (liquid_saturation > this%Sr) then
    capillary_pressure = this%pct/(Se2**(1.d0/this%lambda))
  else
    capillary_pressure = 0.d0
  endif
  
  call SFWIPPKPC(this,this%lambda,this%pct,Se2,capillary_pressure)
  
end subroutine SFKRP3CapillaryPressure

! ************************************************************************** !

subroutine SFKRP3Saturation(this,capillary_pressure, &
                              liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the liquid saturation as a function of capillary pressure using
  ! the modified Brooks Corey formulation.
  ! BRAGFLO UM 6.02 pg 41, 42; eq.126, eq. 129
  ! Modified according to KRP=3 option of BRAGFLO:
  ! Threshold pressure is a function of permeability.
  ! Effective saturation is a function of residual gas saturation.
  !
  ! Author: Jennifer Frederick
  ! Date: 05/04/2017
  !
  use Option_module
  
  implicit none

  class(sat_func_KRP3_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2    
  PetscReal :: term
  
  dsat_dpres = 0.d0
  dsat_dpres = 1.d0/dsat_dpres
  dsat_dpres = 0.d0*dsat_dpres

  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP3_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  term = this%pct*((1.d0-this%Sr)/(1.d0-this%Sr-this%Srg))**(-1.d0/this%lambda)

  if ((capillary_pressure) < term) then
    liquid_saturation = 1.d0
    dsat_dpres = 0.d0
    return
  endif
  
  Se2 = (this%pct/capillary_pressure)**this%lambda
  liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se2
  
end subroutine SFKRP3Saturation

! ************************************************************************** !
! ************************************************************************** !

function SFKRP4Create()

  ! Creates the BRAGFLO KRP4 capillary pressure function object

  implicit none
  
  class(sat_func_KRP4_type), pointer :: SFKRP4Create
  
  allocate(SFKRP4Create)
  call SFKRP4Create%Init()
  
end function SFKRP4Create

! ************************************************************************** !

subroutine SFKRP4Init(this)

  ! Creates the BRAGFLO KRP4 capillary pressure function object

  implicit none
  
  class(sat_func_KRP4_type) :: this

  call SFBaseInit(this)
  call SFWIPPInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SFKRP4Init

! ************************************************************************** !

subroutine SFKRP4Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_KRP4_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: num_errors

  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP4'
  endif
  call SFBaseVerify(this,string,option)
  call SFWIPPVerify(this,string,option)
  
  num_errors = 0
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif   
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif   
  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the ' // trim(string) // ' block. See above.'
    call PrintErrMsg(option)
  endif   

end subroutine SFKRP4Verify

! ************************************************************************** !

subroutine SFKRP4CapillaryPressure(this,liquid_saturation, &
                                     capillary_pressure,dpc_dsatl,option)     
  ! 
  ! Computes the capillary pressure as a function of saturation using the
  ! modified Brooks Corey formulation.
  ! BRAGFLO UM 6.02 pg 41, 42; eq.126, eq. 129
  ! Modified according to KRP=4 option of BRAGFLO:
  ! Threshold pressure is a function of permeability.
  ! Effective saturation is a function of residual gas saturation.
  !
  ! Author: Jennifer Frederick
  ! Date: 06/01/2017
  !
  use Option_module
  
  implicit none
  
  class(sat_func_KRP4_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  
  dpc_dsatl = 0.d0
  dpc_dsatl = 1.d0/dpc_dsatl
  dpc_dsatl = 0.d0*dpc_dsatl
  
  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP4_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  Se2 = (liquid_saturation-this%Sr)/(1.d0-this%Sr-this%Srg)
  
  capillary_pressure = 0.d0
  if ((1.d0-liquid_saturation) <= this%Srg) then
    capillary_pressure = this%pct/(Se2**(1.d0/this%lambda))
  endif
  if (liquid_saturation > this%Sr) then
    capillary_pressure = this%pct/(Se2**(1.d0/this%lambda))
  endif
  
  call SFWIPPKPC(this,this%lambda,this%pct,Se2,capillary_pressure)
  
end subroutine SFKRP4CapillaryPressure

! ************************************************************************** !

subroutine SFKRP4Saturation(this,capillary_pressure, &
                              liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the liquid saturation as a function of capillary pressure using 
  ! the modified Brooks Corey formulation.
  ! BRAGFLO UM 6.02 pg 41, 42; eq.126, eq. 129
  ! Modified according to KRP=4 option of BRAGFLO:
  ! Threshold pressure is a function of permeability.
  ! Effective saturation is a function of residual gas saturation.
  !
  ! Author: Jennifer Frederick
  ! Date: 06/01/2017
  !
  use Option_module
  
  implicit none

  class(sat_func_KRP4_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2  
  PetscReal :: term
  
  dsat_dpres = 0.d0
  dsat_dpres = 1.d0/dsat_dpres
  dsat_dpres = 0.d0*dsat_dpres

  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP1_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  term = this%pct*((1.d0-this%Sr)/(1.d0-this%Sr-this%Srg))**(-1.d0/this%lambda)

  if ((capillary_pressure) < term) then
    liquid_saturation = 1.d0
    dsat_dpres = 0.d0
    return
  endif
  
  Se2 = (this%pct/capillary_pressure)**this%lambda
  liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se2
  
end subroutine SFKRP4Saturation

! ************************************************************************** !
! ************************************************************************** !

function SFKRP5Create()

  ! Creates the BRAGFLO KRP5 capillary pressure function object

  implicit none
  
  class(sat_func_KRP5_type), pointer :: SFKRP5Create
  
  allocate(SFKRP5Create)
  call SFKRP5Create%Init()
  
end function SFKRP5Create

! ************************************************************************** !

subroutine SFKRP5Init(this)

  ! Creates the BRAGFLO KRP5 capillary pressure function object

  implicit none
  
  class(sat_func_KRP5_type) :: this

  call SFBaseInit(this)
  call SFWIPPInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%pcmax = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SFKRP5Init

! ************************************************************************** !

subroutine SFKRP5Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_KRP5_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: num_errors

  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP5'
  endif
  call SFBaseVerify(this,string,option)
  call SFWIPPVerify(this,string,option)
  
  num_errors = 0
  if (Uninitialized(this%pcmax)) then
    option%io_buffer = UninitializedMessage('MAX_CAPILLARY_PRESSURE',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif   
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif   
  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the ' // trim(string) // ' block. See above.'
    call PrintErrMsg(option)
  endif   

end subroutine SFKRP5Verify

! ************************************************************************** !

subroutine SFKRP5CapillaryPressure(this,liquid_saturation, &
                                     capillary_pressure,dpc_dsatl,option)  
  ! 
  ! Computes the capillary pressure as a function of saturation linearly.
  ! BRAGFLO UM 6.02 pg 45; Fig. 21
  ! Modified according to KRP=5 option of BRAGFLO.
  ! Threshold pressure is a function of permeability.
  ! Effective saturation is a function of residual gas pressure.
  !
  ! Author: Jennifer Frederick
  ! Date: 06/01/2017
  !
  use Option_module
  
  implicit none
  
  class(sat_func_KRP5_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  PetscReal :: dummy_lambda
  
  dpc_dsatl = 0.d0
  dpc_dsatl = 1.d0/dpc_dsatl
  dpc_dsatl = 0.d0*dpc_dsatl
  
  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP5_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  Se2 = (liquid_saturation-this%Sr)/(1.d0-this%Sr-this%Srg)
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
  else if ((1.d0 - liquid_saturation) <= this%Srg) then
    capillary_pressure = this%pct
  else 
    capillary_pressure = (this%pct-this%pcmax)*Se2 + this%pcmax
  endif
  
end subroutine SFKRP5CapillaryPressure

! ************************************************************************** !

subroutine SFKRP5Saturation(this,capillary_pressure, &
                              liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the liquid saturation as a function of capillary pressure linearly.
  ! BRAGFLO UM 6.02 pg 45; Fig. 21
  ! Modified according to KRP=5 option of BRAGFLO.
  ! Threshold pressure is a function of permeability.
  ! Effective saturation is a function of residual gas pressure.
  !
  ! Author: Jennifer Frederick
  ! Date: 06/01/2017
  !
  use Option_module
  
  implicit none

  class(sat_func_KRP5_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  
  dsat_dpres = 0.d0
  dsat_dpres = 1.d0/dsat_dpres
  dsat_dpres = 0.d0*dsat_dpres

  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP5_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif

  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  endif
  
  Se2 = (capillary_pressure-this%pcmax)/(this%pct-this%pcmax)
  liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se2

end subroutine SFKRP5Saturation
                            
! ************************************************************************** !
! ************************************************************************** !

function SFKRP8Create()

  ! Creates the BRAGFLO KRP8 capillary pressure function object

  implicit none
  
  class(sat_func_KRP8_type), pointer :: SFKRP8Create
  
  allocate(SFKRP8Create)
  call SFKRP8Create%Init()
  
end function SFKRP8Create

! ************************************************************************** !

subroutine SFKRP8Init(this)

  ! Creates the BRAGFLO KRP8 capillary pressure function object

  implicit none
  
  class(sat_func_KRP8_type) :: this

  call SFBaseInit(this)
  call SFWIPPInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%m = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SFKRP8Init

! ************************************************************************** !

subroutine SFKRP8Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_KRP8_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: num_errors
  
  num_errors = 0
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP8'
  endif
  call SFBaseVerify(this,string,option)
  call SFWIPPVerify(this,string,option)
  
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif    
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif 
  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the ' // trim(string) // ' block. See above.'
    call PrintErrMsg(option)
  endif

end subroutine SFKRP8Verify

! ************************************************************************** !

subroutine SFKRP8CapillaryPressure(this,liquid_saturation, &
                                     capillary_pressure,dpc_dsatl,option)     
  ! 
  ! Computes the capillary pressure as a function of saturation using the
  ! original van Genuchten-Parker formulation.
  ! BRAGFLO UM 6.02 pg 41, 42; eq.120
  ! Modified according to KRP=8 option of BRAGFLO:
  ! Threshold pressure is a function of permeability.
  ! P0 parameter.
  !
  ! Author: Jennifer Frederick
  ! Date: 06/01/2017
  !
  use Option_module
  
  implicit none
  
  class(sat_func_KRP8_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: lambda
  PetscReal :: Se1
  PetscReal :: P0
  
  dpc_dsatl = 0.d0
  dpc_dsatl = 1.d0/dpc_dsatl
  dpc_dsatl = 0.d0*dpc_dsatl
   
  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP8_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  lambda = this%m/(1.d0-this%m)
  ! Derivation for P0 is in Appendix PA:
  ! It is derived by setting Se2 in KRP4 and Se1 KRP8 to the value 0.5, and then 
  ! equating Pc(KRP4,Se2=0.5) = Pc(KRP8,Se1=0.5) and solving for P0.
  P0 = this%pct * (2.d0**(1.d0/lambda)) * &
       (((((0.5d0*(1.d0-this%Srg-this%Sr))/(1.d0-this%Sr)) &
       **(-1.d0/this%m))-1.d0)**(this%m-1.d0))
  Se1 = (liquid_saturation-this%Sr)/(1.d0-this%Sr)
  
  if ((Se1 < 1.d0) .and. (liquid_saturation > this%Sr)) then
    capillary_pressure = P0*(Se1**(-1.d0/this%m)-1.d0)**(1.d0-this%m)
  else 
    capillary_pressure = 0.d0
  endif

  call SFWIPPKPC(this,lambda,P0,Se1,capillary_pressure)
  
end subroutine SFKRP8CapillaryPressure

! ************************************************************************** !

subroutine SFKRP8Saturation(this,capillary_pressure, &
                              liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the liquid saturation as a function of capillary pressure using
  ! the original van Genuchten-Parker formulation.
  ! BRAGFLO UM 6.02 pg 41, 42; eq.120
  ! Modified according to KRP=8 option of BRAGFLO:
  ! Threshold pressure is a function of permeability.
  ! P0 parameter.
  !
  ! Author: Jennifer Frederick
  ! Date: 06/01/2017
  !
  use Option_module
  
  implicit none

  class(sat_func_KRP8_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: lambda
  PetscReal :: Se1
  PetscReal :: P0
  PetscReal :: n
  
  dsat_dpres = 0.d0
  dsat_dpres = 1.d0/dsat_dpres
  dsat_dpres = 0.d0*dsat_dpres
  
  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP8_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  lambda = this%m/(1.d0-this%m)
  ! Derivation for P0 is in Appendix PA:
  ! It is derived by setting Se2 in KRP4 and Se1 KRP8 to the value 0.5, and then 
  ! equating Pc(KRP4,Se2=0.5) = Pc(KRP8,Se1=0.5) and solving for P0.
  P0 = this%pct * (2.d0**(1.d0/lambda)) * &
       (((((0.5d0*(1.d0-this%Srg-this%Sr))/(1.d0-this%Sr)) &
       **(-1.d0/this%m))-1.d0)**(this%m-1.d0))
  Se1 = (liquid_saturation-this%Sr)/(1.d0-this%Sr)
    
  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  else
    n = 1.d0/(1.d0-this%m)
    Se1 = (((capillary_pressure**n)/P0) + 1.d0)**(-1.d0*this%m)
    liquid_saturation = this%Sr + (1.d0-this%Sr)*Se1
  endif
  
end subroutine SFKRP8Saturation

! ************************************************************************** !
! ************************************************************************** !

function SFKRP9Create()

  ! Creates the BRAGFLO KRP9 capillary pressure function object

  implicit none
  
  class(sat_func_KRP9_type), pointer :: SFKRP9Create
  
  allocate(SFKRP9Create)
  call SFKRP9Create%Init()
  
end function SFKRP9Create

! ************************************************************************** !

subroutine SFKRP9Init(this)

  ! Creates the BRAGFLO KRP9 capillary pressure function object

  implicit none
  
  class(sat_func_KRP9_type) :: this

  call SFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine SFKRP9Init

! ************************************************************************** !

subroutine SFKRP9Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_KRP9_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP9'
  endif
  call SFBaseVerify(this,string,option)

end subroutine SFKRP9Verify

! ************************************************************************** !

subroutine SFKRP9CapillaryPressure(this,liquid_saturation, &
                                     capillary_pressure,dpc_dsatl,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation
  ! based on experimental measurements and analyses done by Vauclin et al.
  ! as discussed by Moridis and Pruess, and the BRAGFLO V6.02 Requirements
  ! Document and Verification and Validation Plan, Sandia National Laboratories,
  ! Carlsbad, NM. ERMS #558659.  
  ! Moridis, G. J., and K. Pruess.  1992.  TOUGH Simulations of 
  ! Updegraff's Set of Fluid and Heat Flow Problems.  LBL-32611, ERMS# 138458. 
  ! Berkeley, CA:  Lawrence Berkeley Laboratory.
  !
  ! Author: Heeho Park
  ! Date: 03/26/15
  !
  use Option_module
  
  implicit none
  
  class(sat_func_KRP9_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1
  PetscReal, parameter :: a = 3783.0145d0
  PetscReal, parameter :: b = 2.9d0
  
  dpc_dsatl = capillary_pressure / &
              (liquid_saturation*b*(liquid_saturation - 1.d0))
              
  Se1 = (1.d0-liquid_saturation)/(liquid_saturation)
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = 0.d0
  else 
    capillary_pressure = a*Se1**(1.d0/b)
  endif
  
#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
  endif
#endif
  
end subroutine SFKRP9CapillaryPressure

! ************************************************************************** !

subroutine SFKRP9Saturation(this,capillary_pressure, &
                              liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the liquid saturation as a function of capillary pressure
  ! based on experimental measurements and analyses done by Vauclin et al.
  ! as discussed by Moridis and Pruess, and the BRAGFLO V6.02 Requirements
  ! Document and Verification and Validation Plan, Sandia National Laboratories,
  ! Carlsbad, NM. ERMS #558659.  
  ! Moridis, G. J., and K. Pruess.  1992.  TOUGH Simulations of 
  ! Updegraff's Set of Fluid and Heat Flow Problems.  LBL-32611, ERMS# 138458. 
  ! Berkeley, CA:  Lawrence Berkeley Laboratory.
  !
  ! Author: Heeho Park
  ! Date: 03/26/15
  ! 
  use Option_module
  
  implicit none

  class(sat_func_KRP9_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1
  PetscReal :: dS_dSe
  PetscReal :: dSe_dpc
  PetscReal, parameter :: dpc_dpres = -1.d0
  PetscReal, parameter :: a = 3783.0145d0
  PetscReal, parameter :: b = 2.9d0
  
  dsat_dpres = 0.d0
  
  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  else
    Se1 = (capillary_pressure/a)**(b)
    liquid_saturation = 1.d0 / (Se1+1.d0)
    ! Python analytical derivative (Jenn Frederick)
    dS_dSe = -1.d0/(Se1 + 1.d0)**2
    dSe_dpc = b*(capillary_pressure/a)**b/capillary_pressure
    dsat_dpres = dS_dSe*dSe_dpc*dpc_dpres
  endif 

end subroutine SFKRP9Saturation

! ************************************************************************** !
! ************************************************************************** !

function SFKRP11Create()

  ! Creates the BRAGFLO KRP11 capillary pressure function object

  implicit none
  
  class(sat_func_KRP11_type), pointer :: SFKRP11Create
  
  allocate(SFKRP11Create)
  call SFKRP11Create%Init()
  
end function SFKRP11Create

! ************************************************************************** !

subroutine SFKRP11Init(this)

  ! Creates the BRAGFLO KRP11 capillary pressure function object

  implicit none
  
  class(sat_func_KRP11_type) :: this

  call SFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine SFKRP11Init

! ************************************************************************** !

subroutine SFKRP11Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_KRP11_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: tempreal
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP11'
  endif
  ! SFBaseVerify checks whether residual saturation is initialized, and it
  ! will not be with KRP11, therefore, set to dummy value and back to 
  ! uninitialized.
  tempreal = this%Sr
  if (Uninitialized(this%Sr)) this%Sr = 0.d0
  call SFBaseVerify(this,string,option)
  this%Sr = tempreal

end subroutine SFKRP11Verify

! ************************************************************************** !

subroutine SFKRP11CapillaryPressure(this,liquid_saturation, &
                                      capillary_pressure,dpc_dsatl,option)
  ! 
  ! Computes the capillary pressure as a function of saturation using the
  ! open cavity modification logic.
  ! BRAGFLO UM 6.02 pg 48
  ! Modified according to KRP=11 option of BRAGFLO:
  ! Capillary pressure is zero at all saturations.
  !
  ! Author: Heeho Park
  ! Date: 03/26/15
  !
  use Option_module
  
  implicit none
  
  class(sat_func_KRP11_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option

  dpc_dsatl = 0.d0
  capillary_pressure = 0.0d0
  
end subroutine SFKRP11CapillaryPressure

! ************************************************************************** !

subroutine SFKRP11Saturation(this,capillary_pressure, &
                               liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the liquid saturation as a function of capillary pressure using
  ! the open cavity modification logic.
  ! BRAGFLO UM 6.02 pg 48
  ! Modified according to KRP=11 option of BRAGFLO:
  ! Saturation is 1.0 for any capillary pressure.
  !   
  ! Author: Heeho Park
  ! Date: 03/26/15
  !
  use Option_module
  
  implicit none

  class(sat_func_KRP11_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  dsat_dpres = 0.d0
  liquid_saturation = 1.d0

end subroutine SFKRP11Saturation

! ************************************************************************** !
! ************************************************************************** !

function SFKRP12Create()

  ! Creates the BRAGFLO KRP12 capillary pressure function object

  implicit none
  
  class(sat_func_KRP12_type), pointer :: SFKRP12Create
  
  allocate(SFKRP12Create)
  call SFKRP12Create%Init()
  
end function SFKRP12Create

! ************************************************************************** !

subroutine SFKRP12Init(this)

  ! Creates the BRAGFLO KRP12 capillary pressure function object

  implicit none
  
  class(sat_func_KRP12_type) :: this

  call SFBaseInit(this)
  call SFWIPPInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  this%s_min = UNINITIALIZED_DOUBLE
  this%s_effmin = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SFKRP12Init

! ************************************************************************** !

subroutine SFKRP12Verify(this,name,option)

  use Option_module

  implicit none
  
  class(sat_func_KRP12_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 
  PetscInt :: num_errors
  
  num_errors = 0
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP12'
  endif  
  call SFBaseVerify(this,string,option)
  call SFWIPPVerify(this,string,option)
 
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%s_min)) then
    option%io_buffer = UninitializedMessage('S_MIN',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif
  if (Uninitialized(this%s_effmin)) then
    option%io_buffer = UninitializedMessage('S_EFFMIN',string)
    call PrintMsg(option)
    num_errors = num_errors + 1
  endif
  if (num_errors > 0) then
    write(option%io_buffer,*) num_errors
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                       &the ' // trim(string) // ' block. See above.'
    call PrintErrMsg(option)
  endif
  
end subroutine SFKRP12Verify

! ************************************************************************** !

subroutine SFKRP12CapillaryPressure(this,liquid_saturation, &
                                      capillary_pressure,dpc_dsatl,option)
  ! 
  ! Computes the capillary pressure as a function of saturation using the
  ! modified Brooks Corey formulation for a waste area.
  ! BRAGFLO UM 6.02 pg 47; eq. 130
  ! Modified according to KRP=12 option of BRAGFLO:
  ! Threshold pressure is a function of permeability.
  ! Effective saturation is a function of the parameters s_min and s_effmin.
  !
  ! Author: Jennifer Frederick
  ! Date: 06/02/2017
  !
  use Option_module
  
  implicit none
  
  class(sat_func_KRP12_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se21
  PetscReal :: Se1
  PetscReal :: Se
  
  dpc_dsatl = 0.d0
  dpc_dsatl = 1.d0/dpc_dsatl
  dpc_dsatl = 0.d0*dpc_dsatl

  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP12_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif

  Se = (liquid_saturation - (this%s_min - this%s_effmin)) / &
       (1.d0 - (this%s_min - this%s_effmin))
  Se21 = max(min(Se,1.d0),this%s_effmin)
  
  capillary_pressure = this%pct/(Se21**(1.d0/this%lambda))
  
  ! do not pass in Se21 into the following function, it needs Se1:
  Se1 = (liquid_saturation - this%Sr)/(1.d0 - this%Sr)
  call SFWIPPKPC(this,this%lambda,this%pct,Se1,capillary_pressure)

end subroutine SFKRP12CapillaryPressure

! ************************************************************************** !

subroutine SFKRP12Saturation(this,capillary_pressure, &
                               liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the liquid saturation as a function of capillary pressure using 
  ! the modified Brooks Corey formulation for a waste area.
  ! BRAGFLO UM 6.02 pg 47; eq. 130
  ! Modified according to KRP=12 option of BRAGFLO:
  ! Threshold pressure is a function of permeability.
  ! Effective saturation is a function of the parameters s_min and s_effmin.
  !
  ! Author: Jennifer Frederick
  ! Date: 06/02/2017
  !
  use Option_module
  
  implicit none      
  
  class(sat_func_KRP12_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se21
  
  dsat_dpres = 0.d0
  dsat_dpres = 1.d0/dsat_dpres
  dsat_dpres = 0.d0*dsat_dpres

  if (this%ignore_permeability) then
    this%pct = 1.d0/this%alpha
  else
    ! check if pct has been updated before using
    if (.not.option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP12_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%flow%pct_updated = PETSC_FALSE
  endif
  
  if (capillary_pressure < this%pct) then
    liquid_saturation = 1.d0
    dsat_dpres = 0.d0
    return
  endif
  
  Se21 = (this%pct/capillary_pressure)**this%lambda
  liquid_saturation = Se21*(1.d0 - this%s_min - this%s_effmin) + &
                      this%s_min + this%s_effmin
                               
end subroutine SFKRP12Saturation

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP1LiqCreate()

  ! Creates the BRAGFLO_KRP1_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP1_liq_type), pointer :: RPFKRP1LiqCreate
  
  allocate(RPFKRP1LiqCreate)
  call RPFKRP1LiqCreate%Init()
  
end function RPFKRP1LiqCreate

! ************************************************************************** !

subroutine RPFKRP1LiqInit(this)

  ! Initializes the BRAGFLO_KRP1_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP1_liq_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP1LiqInit

! ************************************************************************** !

subroutine RPFKRP1LiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP1_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP1_LIQ'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call PrintErrMsg(option)
  endif  
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif
  
end subroutine RPFKRP1LiqVerify

! ************************************************************************** !

subroutine RPFKRP1LiqRelPerm(this,liquid_saturation, &
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
  ! Author: Glenn Hammond, Modified by Jennifer Frederick
  ! Date: 12/11/07, 09/23/14, 06/12/2017
  ! 
  use Option_module
  
  implicit none

  class(rpf_KRP1_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1
  PetscReal :: one_over_m
  PetscReal :: Se1_one_over_m
  PetscReal :: dkr_Se1
  PetscReal :: dSe1_sat
  
  Se1 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  Se1 = min(Se1,1.d0)
  
  if (((1.d0-liquid_saturation) <= this%Srg) .or. &
      (liquid_saturation > this%Sr)) then
    one_over_m = 1.d0/this%m
    Se1_one_over_m = Se1**one_over_m
    relative_permeability = sqrt(Se1)*(1.d0-(1.d0-Se1_one_over_m)**this%m)**2.d0
    dkr_Se1 = 0.5d0*relative_permeability/Se1+ &
              2.d0*Se1**(one_over_m-0.5d0)* &
              (1.d0-Se1_one_over_m)**(this%m-1.d0)* &
              (1.d0-(1.d0-Se1_one_over_m)**this%m)
    dSe1_sat = 1.d0 / (1.d0 - this%Sr)
    dkr_sat = dkr_Se1 * dSe1_sat
  else
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  endif
  
end subroutine RPFKRP1LiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP1GasCreate()

  ! Creates the BRAGFLO_KRP1_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP1_gas_type), pointer :: RPFKRP1GasCreate
  
  allocate(RPFKRP1GasCreate)
  call RPFKRP1GasCreate%Init() 
  
end function RPFKRP1GasCreate

! ************************************************************************** !

subroutine RPFKRP1GasInit(this)

  ! Initializes the BRAGFLO_KRP1_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP1_gas_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP1GasInit

! ************************************************************************** !

subroutine RPFKRP1GasVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rpf_KRP1_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP1_GAS'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif 
  
end subroutine RPFKRP1GasVerify

! ************************************************************************** !

subroutine RPFKRP1GasRelPerm(this,liquid_saturation, &
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
  ! Author: Glenn Hammond, Modified by Jennifer Frederick
  ! Date: 12/11/07, 09/23/14, 06/12/2017
  ! 
  use Option_module
  
  implicit none

  class(rpf_KRP1_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  PetscReal :: Seg
  PetscReal :: dkr_Se2
  PetscReal :: dSe2_sat
  
  Se2 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  Se2 = min(Se2,1.d0)
  
  if ((1.d0-liquid_saturation) <= this%Srg) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  else if (liquid_saturation > this%Sr) then
    Seg = 1.d0 - Se2
    relative_permeability = sqrt(Seg)*(1.d0-Se2**(1.d0/this%m))**(2.d0*this%m)
    ! Mathematica analytical solution (Heeho Park)
    dkr_Se2 = -(1.d0-Se2**(1.d0/this%m))**(2.d0*this%m)/(2.d0*sqrt(Seg)) &
              - 2.d0*sqrt(Seg)*Se2**(1.d0/this%m-1.d0) &
              * (1.d0-Se2**(1.d0/this%m))**(2.d0*this%m-1.d0)
    dSe2_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
    dkr_sat = dkr_Se2 * dSe2_sat
  else
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  endif
  
end subroutine RPFKRP1GasRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP2LiqCreate()

  ! Creates the BRAGFLO_KRP2_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP2_liq_type), pointer :: RPFKRP2LiqCreate
  
  allocate(RPFKRP2LiqCreate)
  call RPFKRP2LiqCreate%Init()
  
end function RPFKRP2LiqCreate

! ************************************************************************** !

subroutine RPFKRP2LiqInit(this)

  ! Initializes the BRAGFLO_KRP2_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP2_liq_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP2LiqInit

! ************************************************************************** !

subroutine RPFKRP2LiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP2_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP2_LIQ'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif   
  
end subroutine RPFKRP2LiqVerify

! ************************************************************************** !

subroutine RPFKRP2LiqRelPerm(this,liquid_saturation, &
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
  ! Author: Glenn Hammond, Modified by Jennifer Frederick
  ! Date: 12/11/07, 09/23/14, 06/12/2017
  ! 
  use Option_module
  
  implicit none

  class(rpf_KRP2_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1
  PetscReal :: power
  PetscReal :: dkr_Se1
  PetscReal :: dSe1_sat
  
  Se1 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  else 
    power = 3.d0+2.d0/this%lambda
    relative_permeability = Se1**power
    dkr_Se1 = power*relative_permeability/Se1          
    dSe1_sat = 1.d0 / (1.d0 - this%Sr)
    dkr_sat = dkr_Se1 * dSe1_sat
  endif
  
end subroutine RPFKRP2LiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP2GasCreate()

  ! Creates the BRAGFLO_KRP2_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP2_gas_type), pointer :: RPFKRP2GasCreate
  
  allocate(RPFKRP2GasCreate)
  call RPFKRP2GasCreate%Init()
  
end function RPFKRP2GasCreate

! ************************************************************************** !

subroutine RPFKRP2GasInit(this)

  ! Initializes the BRAGFLO_KRP2_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP2_gas_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP2GasInit

! ************************************************************************** !

subroutine RPFKRP2GasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP2_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP2_GAS'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif   
  
end subroutine RPFKRP2GasVerify

! ************************************************************************** !

subroutine RPFKRP2GasRelPerm(this,liquid_saturation, &
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
  ! Author: Glenn Hammond, Modified by Jennifer Frederick
  ! Date: 12/11/07, 09/23/14, 06/12/2017
  ! 
  use Option_module
  
  implicit none

  class(rpf_KRP2_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1
  PetscReal :: Seg
  PetscReal :: dkr_Se1
  PetscReal :: dSe1_sat
  
  Se1 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else 
    Seg = 1.d0 - Se1
    relative_permeability = Seg*Seg*(1.d0-Se1**(1.d0+2.d0/this%lambda))
    ! Mathematica analytical solution (Heeho Park)
    dkr_Se1 = -(1.d0+2.d0/this%lambda)*Seg**2.d0*Se1**(2.d0/this%lambda) &
              - 2.d0*Seg*(1.d0-Se1**(1.d0+2.d0/this%lambda))
    dSe1_sat = 1.d0 / (1.d0 - this%Sr)
    dkr_sat = dkr_Se1 * dSe1_sat
  endif
    
end subroutine RPFKRP2GasRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP3LiqCreate()

  ! Creates the BRAGFLO_KRP3_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP3_liq_type), pointer :: RPFKRP3LiqCreate
  
  allocate(RPFKRP3LiqCreate)
  call RPFKRP3LiqCreate%Init() 
  
end function RPFKRP3LiqCreate

! ************************************************************************** !

subroutine RPFKRP3LiqInit(this)

  ! Initializes the BRAGFLO_KRP3_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP3_liq_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP3LiqInit

! ************************************************************************** !

subroutine RPFKRP3LiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP3_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP3_LIQ'
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
  
end subroutine RPFKRP3LiqVerify

! ************************************************************************** !

subroutine RPFKRP3LiqRelPerm(this,liquid_saturation, &
                                relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Jennifer Frederick
  ! Date: 06/06/2017
  ! 
  use Option_module
  
  implicit none
  
  class(rpf_KRP3_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  PetscReal :: lambda_exp
  PetscReal :: dkr_dSe2
  PetscReal :: dSe2_dsat
  
  if ((1.d0-liquid_saturation) <= this%Srg) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else if (liquid_saturation > this%Sr) then
    Se2 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
    lambda_exp = (2.d0+(3.d0*this%lambda))/this%lambda
    relative_permeability = Se2**(lambda_exp)
    dSe2_dsat = 1.d0 / (1.d0 - this%Sr - this%Srg)
    dkr_dSe2 = lambda_exp*relative_permeability/Se2 
    dkr_sat = dkr_dSe2 * dSe2_dsat
  else
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  endif
  
end subroutine RPFKRP3LiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP3GasCreate()

  ! Creates the BRAGFLO_KRP3_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP3_gas_type), pointer :: RPFKRP3GasCreate
  
  allocate(RPFKRP3GasCreate)
  call RPFKRP3GasCreate%Init() 
  
end function RPFKRP3GasCreate

! ************************************************************************** !

subroutine RPFKRP3GasInit(this)

  ! Initializes the BRAGFLO_KRP3_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP3_gas_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP3GasInit

! ************************************************************************** !

subroutine RPFKRP3GasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP3_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP3_GAS'
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
  
end subroutine RPFKRP3GasVerify

! ************************************************************************** !

subroutine RPFKRP3GasRelPerm(this,liquid_saturation, &
                                relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Jennifer Frederick
  ! Date: 06/06/2017
  ! 
  use Option_module
  
  implicit none
  
  class(rpf_KRP3_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  PetscReal :: lambda_exp
  PetscReal :: dkr_dSe2
  PetscReal :: dSe2_dsat
  
  if ((1.d0-liquid_saturation) <= this%Srg) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  else if (liquid_saturation > this%Sr) then
    Se2 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
    lambda_exp = (2.d0+this%lambda)/this%lambda
    relative_permeability = ((1.d0-Se2)**2.d0) * (1.d0-(Se2**lambda_exp))
    dSe2_dsat = 1.d0 / (1.d0 - this%Sr - this%Srg)
    ! Python analytical derivative (Jenn Frederick)
    dkr_dSe2 = -1.d0*(Se2-1.d0)*(lambda_exp*(Se2**lambda_exp)*(Se2-1.d0) + &
               2.d0*Se2*(Se2**lambda_exp-1.d0))/Se2
    dkr_sat = dkr_dSe2 * dSe2_dsat
  else
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  endif
   
end subroutine RPFKRP3GasRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP4LiqCreate()

  ! Creates the BRAGFLO_KRP4_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP4_liq_type), pointer :: RPFKRP4LiqCreate
  
  allocate(RPFKRP4LiqCreate)
  call RPFKRP4LiqCreate%Init() ! Calls KRP3_Liq's Init()
  
end function RPFKRP4LiqCreate

! ************************************************************************** !

subroutine RPFKRP4LiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP4_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP4_LIQ'
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
  
end subroutine RPFKRP4LiqVerify

! ************************************************************************** !

subroutine RPFKRP4LiqRelPerm(this,liquid_saturation, &
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
  ! Author: Jennifer Frederick
  ! Date: 06/06/2017
  ! 
  use Option_module
  
  implicit none

  class(rpf_KRP4_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1
  PetscReal :: power
  PetscReal :: dkr_dSe1
  PetscReal :: dSe1_dsat
  
  if (((1.d0-liquid_saturation) <= this%Srg) .or. &
      (liquid_saturation > this%Sr)) then
    Se1 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
    power = (2.d0+(3.d0*this%lambda))/this%lambda
    relative_permeability = Se1**power
    dkr_dSe1 = power*relative_permeability/Se1          
    dSe1_dsat = 1.d0 / (1.d0 - this%Sr)
    dkr_sat = dkr_dSe1 * dSe1_dsat
  else
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  endif
  
end subroutine RPFKRP4LiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP4GasCreate()

  ! Creates the BRAGFLO_KRP4_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP4_gas_type), pointer :: RPFKRP4GasCreate
  
  allocate(RPFKRP4GasCreate)
  call RPFKRP4GasCreate%Init() ! calls KRP3_Gas's Init()
  
end function RPFKRP4GasCreate

! ************************************************************************** !

subroutine RPFKRP4GasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP4_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP4_GAS'
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
  
end subroutine RPFKRP4GasVerify

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP5LiqCreate()

  ! Creates the BRAGFLO_KRP5_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP5_liq_type), pointer :: RPFKRP5LiqCreate
  
  allocate(RPFKRP5LiqCreate)
  call RPFKRP5LiqCreate%Init()
  
end function RPFKRP5LiqCreate

! ************************************************************************** !

subroutine RPFKRP5LiqInit(this)

  ! Initializes the BRAGFLO_KRP5_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP5_liq_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP5LiqInit

! ************************************************************************** !

subroutine RPFKRP5LiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP5_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP5_LIQ'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif  
    
end subroutine RPFKRP5LiqVerify

! ************************************************************************** !

subroutine RPFKRP5LiqRelPerm(this,liquid_saturation, &
                                relative_permeability,dkr_sat,option)
  !
  ! Author: Heeho Park; Modified by Jennifer Frederick
  ! Date: 11/18/16; 06/06/2017
  !
  use Option_module
  
  implicit none

  class(rpf_KRP5_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  else if ((1.d0-liquid_saturation) <= this%Srg) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else
    Se2 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
    relative_permeability = Se2
    dkr_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  endif
   
end subroutine RPFKRP5LiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP5GasCreate()

  ! Creates the BRAGFLO_KRP5_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP5_gas_type), pointer :: RPFKRP5GasCreate
  
  allocate(RPFKRP5GasCreate)
  call RPFKRP5GasCreate%Init()
  
end function RPFKRP5GasCreate

! ************************************************************************** !

subroutine RPFKRP5GasInit(this)

  ! Initializes the BRAGFLO_KRP5_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP5_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP5GasInit

! ************************************************************************** !

subroutine RPFKRP5GasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP5_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP5_GAS'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif  
    
end subroutine RPFKRP5GasVerify

! ************************************************************************** !

subroutine RPFKRP5GasRelPerm(this,liquid_saturation, &
                                relative_permeability,dkr_sat,option)
  !
  ! Author: Jennifer Frederick
  ! Date: 06/12/2017
  !
  use Option_module
  
  implicit none

  class(rpf_KRP5_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  PetscReal :: dkr_Se2
  PetscReal :: dSe2_sat
  
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else if ((1.d0-liquid_saturation) <= this%Srg) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  else
    Se2 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
    relative_permeability = 1.d0 - Se2
    dkr_Se2 = -1.d0
    dSe2_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
    dkr_sat = dkr_Se2 * dSe2_sat
  endif
   
end subroutine RPFKRP5GasRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP8LiqCreate()

  ! Creates the BRAGFLO_KRP8_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP8_liq_type), pointer :: RPFKRP8LiqCreate
  
  allocate(RPFKRP8LiqCreate)
  call RPFKRP8LiqCreate%Init()
  
end function RPFKRP8LiqCreate

! ************************************************************************** !

subroutine RPFKRP8LiqInit(this)

  ! Initializes the BRAGFLO_KRP8_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP8_liq_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP8LiqInit

! ************************************************************************** !

subroutine RPFKRP8LiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP8_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP8_LIQ'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call PrintErrMsg(option)
  endif   
  
end subroutine RPFKRP8LiqVerify

! ************************************************************************** !

subroutine RPFKRP8LiqRelPerm(this,liquid_saturation, &
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
  ! Author: Jennifer Frederick
  ! Date: 06/07/2017
  ! 
  use Option_module
  
  implicit none

  class(rpf_KRP8_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1
  PetscReal :: dkr_dSe1
  PetscReal :: dSe1_dsat
  
  Se1 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  
  if (liquid_saturation > this%Sr) then
    if (Se1 < 1.d0) then
      relative_permeability = sqrt(Se1)* &
                              (1.d0-((1.d0-(Se1**(1.d0/this%m)))**this%m))**2.d0
      dkr_dSe1 = 0.5d0*relative_permeability/Se1+ &
                 2.d0*Se1**((1.d0/this%m)-0.5d0)* &
                 (1.d0-(Se1**(1.d0/this%m)))**(this%m-1.d0)* &
                 (1.d0-(1.d0-(Se1**(1.d0/this%m)))**this%m)
      dSe1_dsat = 1.d0 / (1.d0 - this%Sr)
      dkr_sat = dkr_dSe1 * dSe1_dsat
    else
      relative_permeability = 1.d0
      dkr_sat = 0.d0
    endif
  else
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  endif
  
end subroutine RPFKRP8LiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP8GasCreate()

  ! Creates the BRAGFLO_KRP8_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP8_gas_type), pointer :: RPFKRP8GasCreate
  
  allocate(RPFKRP8GasCreate)
  call RPFKRP8GasCreate%Init()
  
end function RPFKRP8GasCreate

! ************************************************************************** !

subroutine RPFKRP8GasInit(this)

  ! Initializes the BRAGFLO_KRP8_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP8_gas_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP8GasInit

! ************************************************************************** !

subroutine RPFKRP8GasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP8_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP8_GAS'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call PrintErrMsg(option)
  endif   
  
end subroutine RPFKRP8GasVerify

! ************************************************************************** !

subroutine RPFKRP8GasRelPerm(this,liquid_saturation, &
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
  ! Author: Jennifer Frederick
  ! Date: 06/07/2017
  ! 
  use Option_module
  
  implicit none

  class(rpf_KRP8_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1
  PetscReal :: dkr_dSe1
  PetscReal :: dSe1_dsat
  
  Se1 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  
  if (liquid_saturation > this%Sr) then
    if (Se1 < 1.d0) then
      relative_permeability = sqrt(1.d0-Se1)* &
                              (1.d0-Se1**(1.d0/this%m))**(2.d0*this%m)
      ! Mathematica analytical derivative (Heeho Park)
      dkr_dSe1 = -(1.d0-Se1**(1.d0/this%m))**(2.d0*this%m)/ &
          (2.d0*sqrt(1.d0-Se1)) - 2.d0*sqrt(1.d0-Se1)*Se1**(1.d0/this%m-1.d0) &
          * (1.d0-Se1**(1.d0/this%m))**(2.d0*this%m-1.d0)
      dSe1_dsat = 1.d0 / (1.d0 - this%Sr)
      dkr_sat = dkr_dSe1 * dSe1_dsat
    else
      relative_permeability = 0.d0
      dkr_sat = 0.d0
    endif
  else
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  endif
  
end subroutine RPFKRP8GasRelPerm

! ************************************************************************** !
! ************************************************************************** !
                                     
function RPFKRP9LiqCreate()

  ! Creates the BRAGFLO_KRP9_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP9_liq_type), pointer :: RPFKRP9LiqCreate
  
  allocate(RPFKRP9LiqCreate)
  call RPFKRP9LiqCreate%Init()
  
end function RPFKRP9LiqCreate

! ************************************************************************** !

subroutine RPFKRP9LiqInit(this)

  ! Initializes the BRAGFLO_KRP9_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP9_liq_type) :: this

  call RPFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP9LiqInit

! ************************************************************************** !

subroutine RPFKRP9LiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP9_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP9_LIQ'
  endif  
  call RPFBaseVerify(this,string,option)
  
end subroutine RPFKRP9LiqVerify

! ************************************************************************** !

subroutine RPFKRP9LiqRelPerm(this,liquid_saturation, &
                                relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation based on experimental measurements and analyses 
  ! done by Vauclin et al. as discussed by Moridis and Pruess. 
  ! Moridis, G. J., and K. Pruess.  1992.  TOUGH Simulations of 
  ! Updegraff's Set of Fluid and Heat Flow Problems.  LBL-32611, ERMS# 138458. 
  ! Berkeley, CA:  Lawrence Berkeley Laboratory.
  ! 
  ! Author: Heeho Park
  ! Date: 03/26/15
  ! 
  use Option_module
  
  implicit none

  class(rpf_KRP9_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal, parameter :: a = 28.768353d0
  PetscReal, parameter :: b = 1.7241379d0
  PetscReal :: Se
  PetscReal :: dkr_dSe
  PetscReal :: dSe_dsat
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (1.d0-liquid_saturation)/(liquid_saturation)
  
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 0.d0
    return
  endif
  
  relative_permeability = 1.d0/(1.d0+a*Se**b)
  ! Python analytical derivative (Jenn Frederick)
  dkr_dSe = -1.d0*Se**(b-1.d0)*a*b/(Se**b*a + 1.d0)**2.d0
  dSe_dsat = -1.d0/(liquid_saturation**2.d0)
  dkr_sat = dkr_dSe * dSe_dsat
  
end subroutine RPFKRP9LiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP9GasCreate()

  ! Creates the BRAGFLO_KRP9_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP9_gas_type), pointer :: RPFKRP9GasCreate
  
  allocate(RPFKRP9GasCreate)
  call RPFKRP9GasCreate%Init()
  
end function RPFKRP9GasCreate

! ************************************************************************** !

subroutine RPFKRP9GasInit(this)

  ! Initializes the BRAGFLO_KRP9_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP9_gas_type) :: this

  call RPFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP9GasInit

! ************************************************************************** !

subroutine RPFKRP9GasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP9_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP9_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  
end subroutine RPFKRP9GasVerify

! ************************************************************************** !

subroutine RPFKRP9GasRelPerm(this,liquid_saturation, &
                                relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation based on experimental measurements and analyses 
  ! done by Vauclin et al. as discussed by Moridis and Pruess. 
  ! Moridis, G. J., and K. Pruess.  1992.  TOUGH Simulations of 
  ! Updegraff's Set of Fluid and Heat Flow Problems.  LBL-32611, ERMS# 138458. 
  ! Berkeley, CA:  Lawrence Berkeley Laboratory.
  !
  ! Author: Heeho Park
  ! Date: 03/26/15
  ! 

  use Option_module
  
  implicit none

  class(rpf_KRP9_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: liquid_relative_permeability
  PetscReal :: liquid_dkr_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0

  Se = (1.d0-liquid_saturation)/(liquid_saturation)
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 1.d0
    return
  endif
  
  call RPFKRP9LiqRelPerm(this,liquid_saturation, &
                            liquid_relative_permeability, &
                            liquid_dkr_sat,option)
  
  relative_permeability = 1.d0 - liquid_relative_permeability
  dkr_sat = -1.d0 * liquid_dkr_sat
  
end subroutine RPFKRP9GasRelPerm

! ************************************************************************** !
! ************************************************************************** !
 
function RPFKRP11LiqCreate()

  ! Creates the BRAGFLO_KRP11_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP11_liq_type), pointer :: RPFKRP11LiqCreate
  
  allocate(RPFKRP11LiqCreate)
  call RPFKRP11LiqCreate%Init()
  
end function RPFKRP11LiqCreate

! ************************************************************************** !

subroutine RPFKRP11LiqInit(this)

  ! Initializes the BRAGFLO_KRP11_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP11_liq_type) :: this

  call RPFBaseInit(this)
  this%tolc = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFKRP11LiqInit

! ************************************************************************** !

subroutine RPFKRP11LiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP11_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP11_LIQ'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%tolc)) then
    option%io_buffer = UninitializedMessage('TOLC',string)
    call PrintErrMsg(option)
  endif

end subroutine RPFKRP11LiqVerify

! ************************************************************************** !

subroutine RPFKRP11LiqRelPerm(this,liquid_saturation, &
                                 relative_permeability,dkr_sat,option)
  !
  ! Author: Heeho Park
  ! Date: 03/26/15
  ! 
  use Option_module
  
  implicit none

  class(rpf_KRP11_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: gas_saturation
  PetscReal :: tol
  
  gas_saturation = 1.d0 - liquid_saturation
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  tol = this%tolc * (1 - this%Sr - this%Srg)
  
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  else if (gas_saturation <= this%Srg) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else if (liquid_saturation <= (this%Sr+tol)) then
    relative_permeability = (liquid_saturation - this%Sr)/tol
    dkr_sat = 1.d0/tol
  else if (gas_saturation <= (this%Srg+tol)) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  endif
    
end subroutine RPFKRP11LiqRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP11GasCreate()

  ! Creates the BRAGFLO_KRP11_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP11_gas_type), pointer :: RPFKRP11GasCreate
  
  allocate(RPFKRP11GasCreate)
  call RPFKRP11GasCreate%Init() ! calls KRP11_LIQ's Init()
  
end function RPFKRP11GasCreate

! ************************************************************************** !

subroutine RPFKRP11GasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP11_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP11_GAS'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%tolc)) then
    option%io_buffer = UninitializedMessage('TOLC',string)
    call PrintErrMsg(option)
  endif

end subroutine RPFKRP11GasVerify

! ************************************************************************** !

subroutine RPFKRP11GasRelPerm(this,liquid_saturation, &
                                 relative_permeability,dkr_sat,option)
  !
  ! Author: Heeho Park
  ! Date: 03/26/15
  ! 

  use Option_module
  
  implicit none

  class(rpf_KRP11_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: gas_saturation
  PetscReal :: tol
  
  gas_saturation = 1.d0 - liquid_saturation
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  tol = this%tolc * (1 - this%Sr - this%Srg)
  
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else if (gas_saturation <= this%Srg) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  else if (liquid_saturation <= (this%Sr+tol)) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else if (gas_saturation <= (this%Srg+tol)) then
    relative_permeability = (gas_saturation - this%Srg)/tol
    dkr_sat = -1.d0/tol
  else
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  endif
  
  end subroutine RPFKRP11GasRelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPFKRP12LiqCreate()

  ! Creates the BRAGFLO_KRP12_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP12_liq_type), pointer :: RPFKRP12LiqCreate
  
  allocate(RPFKRP12LiqCreate)
  call RPFKRP12LiqCreate%Init()
  
end function RPFKRP12LiqCreate

! ************************************************************************** !

subroutine RPFKRP12LiqVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP12_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP12_LIQ'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif

end subroutine RPFKRP12LiqVerify

! ************************************************************************** !

subroutine RPFKRP12LiqRelPerm(this,liquid_saturation, &
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
  ! Author: Jennifer Frederick
  ! Date: 06/07/2017
  ! 
  use Option_module
  
  implicit none

  class(rpf_KRP12_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se1
  PetscReal :: power
  PetscReal :: dkr_dSe1
  PetscReal :: dSe1_dsat
  
  if (((1.d0-liquid_saturation) <= this%Srg) .or. &
      (liquid_saturation > this%Sr)) then
    Se1 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
    Se1 = max(min(Se1,1.d0),0.d0)
    power = (2.d0+(3.d0*this%lambda))/this%lambda
    relative_permeability = Se1**power
    dkr_dSe1 = power*relative_permeability/Se1          
    dSe1_dsat = 1.d0 / (1.d0 - this%Sr)
    dkr_sat = dkr_dSe1 * dSe1_dsat
  else
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  endif
  
end subroutine RPFKRP12LiqRelPerm
  
! ************************************************************************** !
! ************************************************************************** !

function RPFKRP12GasCreate()

  ! Creates the BRAGFLO_KRP12_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP12_gas_type), pointer :: RPFKRP12GasCreate
  
  allocate(RPFKRP12GasCreate)
  call RPFKRP12GasCreate%Init()
  
end function RPFKRP12GasCreate

! ************************************************************************** !

subroutine RPFKRP12GasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_KRP12_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP12_GAS'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif

end subroutine RPFKRP12GasVerify

! ************************************************************************** !

subroutine RPFKRP12GasRelPerm(this,liquid_saturation, &
                                 relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Jennifer Frederick
  ! Date: 06/07/2017
  ! 
  use Option_module
  
  implicit none
  
  class(rpf_KRP12_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  PetscReal :: lambda_exp
  PetscReal :: dkr_dSe2
  PetscReal :: dSe2_dsat
  
  if ((1.d0-liquid_saturation) <= this%Srg) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  else if (liquid_saturation > this%Sr) then
    Se2 = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
    Se2 = max(min(Se2,1.d0),0.d0)
    lambda_exp = (2.d0+this%lambda)/this%lambda
    relative_permeability = ((1.d0-Se2)**2.d0) * (1.d0-(Se2**lambda_exp))
    dSe2_dsat = 1.d0 / (1.d0 - this%Sr - this%Srg)
    ! Python analytical derivative (Jenn Frederick)
    dkr_dSe2 = -1.d0*(Se2-1.d0)*(lambda_exp*(Se2**lambda_exp)*(Se2-1.d0) + &
               2.d0*Se2*(Se2**lambda_exp-1.d0))/Se2
    dkr_sat = dkr_dSe2 * dSe2_dsat
  else
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  endif
  
end subroutine RPFKRP12GasRelPerm
  
! ************************************************************************** !
! ************************************************************************** !

function RPFTOUGH2IRP7GasCreate()

  ! Creates the Brooks-Corey Burdine gas relative permeability function object

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type), pointer :: RPFTOUGH2IRP7GasCreate
  
  allocate(RPFTOUGH2IRP7GasCreate)
  call RPFTOUGH2IRP7GasCreate%Init()
  
end function RPFTOUGH2IRP7GasCreate

! ************************************************************************** !

subroutine RPFTOUGH2IRP7GasInit(this)

  ! Initializes the Brooks-Corey Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPFTOUGH2IRP7GasInit

! ************************************************************************** !

subroutine RPFTOUGH2IRP7GasVerify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,TOUGH2_IRP7_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif  
  
end subroutine RPFTOUGH2IRP7GasVerify

! ************************************************************************** !

subroutine RPFTOUGH2IRP7GasRelPerm(this,liquid_saturation, &
                                       relative_permeability,dkr_sat,option)
  ! 
  ! TOUGH2 IRP(7) equations from Appendix G of TOUGH2 user manual
  !
  use Option_module
  
  implicit none

  class(rpf_TOUGH2_IRP7_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: liquid_relative_permeability
  PetscReal :: liquid_dkr_sat
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
                 ! essentially zero
  if (this%Srg <= 0.d0) then
    call RPFMualemVGLiqRelPerm(this,liquid_saturation, &
                               liquid_relative_permeability, &
                               liquid_dkr_sat,option)
    relative_permeability = 1.d0 - liquid_relative_permeability 
  else if ((1.d0 - liquid_saturation) <= this%Srg) then
    relative_permeability = 0.d0
  else if (this%Srg > 0.d0) then
    if (liquid_saturation < this%Sr) then
      relative_permeability = 1.d0
    else
      Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
      Seg = 1.d0 - Se
      relative_permeability = Seg**2*(1.d0-Se*Se)
      ! Mathematica Analytical solution (Heeho Park)
      dkr_Se = -2.d0*Seg**2.d0*Se - 2.d0*Seg*(1.d0-Se**2.d0)
      dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
      dkr_sat = dkr_Se * dSe_sat
    end if
  end if
    
end subroutine RPFTOUGH2IRP7GasRelPerm

end module Characteristic_Curves_WIPP_module

