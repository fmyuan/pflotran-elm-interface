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
    procedure, public :: Init => SF_WIPP_Init
    procedure, public :: Verify => SF_WIPP_Verify
    procedure, public :: CapillaryPressure => SF_WIPP_CapillaryPressure
    procedure, public :: Saturation => SF_WIPP_Saturation
  end type sat_func_WIPP_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP1_type
    PetscReal :: Srg
    PetscReal :: m
  contains
    procedure, public :: Init => SF_KRP1_Init
    procedure, public :: Verify => SF_KRP1_Verify
    procedure, public :: CapillaryPressure => SF_KRP1_CapillaryPressure
    procedure, public :: Saturation => SF_KRP1_Saturation
  end type sat_func_KRP1_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP2_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => SF_KRP2_Init
    procedure, public :: Verify => SF_KRP2_Verify
    procedure, public :: CapillaryPressure => SF_KRP2_CapillaryPressure
    procedure, public :: Saturation => SF_KRP2_Saturation
  end type sat_func_KRP2_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP3_type
    PetscReal :: Srg
    PetscReal :: lambda
  contains
    procedure, public :: Init => SF_KRP3_Init
    procedure, public :: Verify => SF_KRP3_Verify
    procedure, public :: CapillaryPressure => SF_KRP3_CapillaryPressure
    procedure, public :: Saturation => SF_KRP3_Saturation
  end type sat_func_KRP3_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP4_type
    PetscReal :: Srg
    PetscReal :: lambda
  contains
    procedure, public :: Verify => SF_KRP4_Verify
    procedure, public :: CapillaryPressure => SF_KRP4_CapillaryPressure
    procedure, public :: Saturation => SF_KRP4_Saturation
  end type sat_func_KRP4_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP5_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => SF_KRP5_Init
    procedure, public :: Verify => SF_KRP5_Verify
    procedure, public :: CapillaryPressure => SF_KRP5_CapillaryPressure
    procedure, public :: Saturation => SF_KRP5_Saturation
  end type sat_func_KRP5_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP8_type
    PetscReal :: Srg
    PetscReal :: m
  contains
    procedure, public :: Init => SF_KRP8_Init
    procedure, public :: Verify => SF_KRP8_Verify
    procedure, public :: CapillaryPressure => SF_KRP8_CapillaryPressure
    procedure, public :: Saturation => SF_KRP8_Saturation
  end type sat_func_KRP8_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_KRP9_type
  contains
    procedure, public :: Init => SF_KRP9_Init
    procedure, public :: Verify => SF_KRP9_Verify
    procedure, public :: CapillaryPressure => SF_KRP9_CapillaryPressure
    procedure, public :: Saturation => SF_KRP9_Saturation
  end type sat_func_KRP9_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_KRP11_type
  contains
    procedure, public :: Init => SF_KRP11_Init
    procedure, public :: Verify => SF_KRP11_Verify
    procedure, public :: CapillaryPressure => SF_KRP11_CapillaryPressure
    procedure, public :: Saturation => SF_KRP11_Saturation
  end type sat_func_KRP11_type 
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_WIPP_type) :: sat_func_KRP12_type
    PetscReal :: lambda
    PetscReal :: s_min
    PetscReal :: s_effmin
  contains
    procedure, public :: Init => SF_KRP12_Init
    procedure, public :: Verify => SF_KRP12_Verify
    procedure, public :: CapillaryPressure => SF_KRP12_CapillaryPressure
    procedure, public :: Saturation => SF_KRP12_Saturation
  end type sat_func_KRP12_type 
  
!-----------------------------------------------------------------------------
!-- Relative Permeability Functions ------------------------------------------
!-----------------------------------------------------------------------------  
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP1_liq_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_KRP1_Liq_Init
    procedure, public :: Verify => RPF_KRP1_Liq_Verify
    procedure, public :: RelativePermeability => RPF_KRP1_Liq_RelPerm
  end type rpf_KRP1_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP1_gas_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_KRP1_Gas_Init
    procedure, public :: Verify => RPF_KRP1_Gas_Verify
    procedure, public :: RelativePermeability => RPF_KRP1_Gas_RelPerm
  end type rpf_KRP1_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP2_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPF_KRP2_Liq_Init
    procedure, public :: Verify => RPF_KRP2_Liq_Verify
    procedure, public :: RelativePermeability => RPF_KRP2_Liq_RelPerm
  end type rpf_KRP2_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP2_gas_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPF_KRP2_Gas_Init
    procedure, public :: Verify => RPF_KRP2_Gas_Verify
    procedure, public :: RelativePermeability => RPF_KRP2_Gas_RelPerm
  end type rpf_KRP2_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP3_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPF_KRP3_Liq_Init
    procedure, public :: Verify => RPF_KRP3_Liq_Verify
    procedure, public :: RelativePermeability => RPF_KRP3_Liq_RelPerm
  end type rpf_KRP3_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP3_gas_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPF_KRP3_Gas_Init
    procedure, public :: Verify => RPF_KRP3_Gas_Verify
    procedure, public :: RelativePermeability => RPF_KRP3_Gas_RelPerm
  end type rpf_KRP3_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP3_liq_type) :: rpf_KRP4_liq_type
  contains
    procedure, public :: Verify => RPF_KRP4_Liq_Verify
    procedure, public :: RelativePermeability => RPF_KRP4_Liq_RelPerm
  end type rpf_KRP4_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP3_gas_type) :: rpf_KRP4_gas_type
  contains
    procedure, public :: Verify => RPF_KRP4_Gas_Verify
  end type rpf_KRP4_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP5_liq_type
  contains
    procedure, public :: Init => RPF_KRP5_Liq_Init
    procedure, public :: Verify => RPF_KRP5_Liq_Verify
    procedure, public :: RelativePermeability => RPF_KRP5_Liq_RelPerm
  end type rpf_KRP5_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP5_gas_type
  contains
    procedure, public :: Init => RPF_KRP5_Gas_Init
    procedure, public :: Verify => RPF_KRP5_Gas_Verify
    procedure, public :: RelativePermeability => RPF_KRP5_Gas_RelPerm
  end type rpf_KRP5_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP8_liq_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_KRP8_Liq_Init
    procedure, public :: Verify => RPF_KRP8_Liq_Verify
    procedure, public :: RelativePermeability => RPF_KRP8_Liq_RelPerm
  end type rpf_KRP8_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP8_gas_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_KRP8_Gas_Init
    procedure, public :: Verify => RPF_KRP8_Gas_Verify
    procedure, public :: RelativePermeability => RPF_KRP8_Gas_RelPerm
  end type rpf_KRP8_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP9_liq_type
  contains
    procedure, public :: Init => RPF_KRP9_Liq_Init
    procedure, public :: Verify => RPF_KRP9_Liq_Verify
    procedure, public :: RelativePermeability => RPF_KRP9_Liq_RelPerm
  end type rpf_KRP9_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP9_liq_type) :: rpf_KRP9_gas_type
  contains
    procedure, public :: Init => RPF_KRP9_Gas_Init
    procedure, public :: Verify => RPF_KRP9_Gas_Verify
    procedure, public :: RelativePermeability => RPF_KRP9_Gas_RelPerm
  end type rpf_KRP9_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_KRP11_liq_type
    PetscReal :: tolc
  contains
    procedure, public :: Init => RPF_KRP11_Liq_Init
    procedure, public :: Verify => RPF_KRP11_Liq_Verify
    procedure, public :: RelativePermeability => RPF_KRP11_Liq_RelPerm
  end type rpf_KRP11_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP11_liq_type) :: rpf_KRP11_gas_type
  contains
    procedure, public :: Verify => RPF_KRP11_Gas_Verify
    procedure, public :: RelativePermeability => RPF_KRP11_Gas_RelPerm
  end type rpf_KRP11_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP4_liq_type) :: rpf_KRP12_liq_type
  contains
    procedure, public :: Verify => RPF_KRP12_Liq_Verify
    procedure, public :: RelativePermeability => RPF_KRP12_Liq_RelPerm
  end type rpf_KRP12_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_KRP3_gas_type) :: rpf_KRP12_gas_type
  contains
    procedure, public :: Verify => RPF_KRP12_Gas_Verify
    procedure, public :: RelativePermeability => RPF_KRP12_Gas_RelPerm
  end type rpf_KRP12_gas_type
  !---------------------------------------------------------------------------
  ! since the TOUGH2_Corey relative permeability function (IRP=7 in 
  ! TOUGH2 manual) calculates relative perm as a function of the 
  ! Mualem-based  liquid relative permeability when Srg = 0., we extend 
  ! the rpf_Mualem_type to save code
  type, public, extends(rpf_Mualem_VG_liq_type) :: rpf_TOUGH2_IRP7_gas_type
  contains
    procedure, public :: Init => RPF_TOUGH2_IRP7_Gas_Init
    procedure, public :: Verify => RPF_TOUGH2_IRP7_Gas_Verify
    procedure, public :: RelativePermeability => RPF_TOUGH2_IRP7_Gas_RelPerm
  end type rpf_TOUGH2_IRP7_gas_type
  
  public :: &! WIPP saturation functions:
            SF_KRP1_Create, &
            SF_KRP2_Create, &
            SF_KRP3_Create, &
            SF_KRP4_Create, &
            SF_KRP5_Create, &
            SF_KRP8_Create, &
            SF_KRP9_Create, &
            SF_KRP11_Create, &
            SF_KRP12_Create, &
            ! WIPP rel. perm. curves:
            RPF_KRP1_Liq_Create, &
            RPF_KRP1_Gas_Create, &
            RPF_KRP2_Liq_Create, &
            RPF_KRP2_Gas_Create, &
            RPF_KRP3_Liq_Create, &
            RPF_KRP3_Gas_Create, &
            RPF_KRP4_Liq_Create, &
            RPF_KRP4_Gas_Create, &
            RPF_KRP5_Liq_Create, &
            RPF_KRP5_Gas_Create, &
            RPF_KRP8_Liq_Create, &
            RPF_KRP8_Gas_Create, &
            RPF_KRP9_Liq_Create, &
            RPF_KRP9_Gas_Create, &
            RPF_KRP11_Liq_Create, &
            RPF_KRP11_Gas_Create, &
            RPF_KRP12_Liq_Create, &
            RPF_KRP12_Gas_Create, &
            RPF_TOUGH2_IRP7_Gas_Create
  
contains

! ************************************************************************** !

subroutine SF_WIPP_Init(this)

  ! Initializes a sat_func_WIPP_type object.

  implicit none
  
  class(sat_func_WIPP_type) :: this

  call SFBaseInit(this)
  this%kpc = UNINITIALIZED_INTEGER
  this%pct_a = UNINITIALIZED_DOUBLE
  this%pct_exp = UNINITIALIZED_DOUBLE
  this%ignore_permeability = PETSC_FALSE
  this%alpha = UNINITIALIZED_DOUBLE
  
end subroutine SF_WIPP_Init

! ************************************************************************** !

subroutine SF_WIPP_Verify(this,name,option)

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

end subroutine SF_WIPP_Verify

! ************************************************************************** !

subroutine SF_WIPP_CapillaryPressure(this,liquid_saturation, & 
                                     capillary_pressure,dpc_dsatl,option)
  use Option_module
  
  implicit none
  
  class(sat_func_WIPP_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SF_WIPP_CapillaryPressure must be extended.'
  call PrintErrMsg(option)
  
end subroutine SF_WIPP_CapillaryPressure

! ************************************************************************** !

subroutine SF_WIPP_Saturation(this,capillary_pressure, &
                              liquid_saturation,dsat_dpres,option)
  use Option_module

  implicit none
  
  class(sat_func_WIPP_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SF_WIPP_Saturation must be extended.'
  call PrintErrMsg(option)
  
end subroutine SF_WIPP_Saturation

! ************************************************************************** !

subroutine SF_WIPP_KPC(this,lambda,PT,Se,capillary_pressure)
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
  
end subroutine SF_WIPP_KPC

! ************************************************************************** !
! ************************************************************************** !

function SF_KRP1_Create()

  ! Creates the BRAGFLO KRP1 capillary pressure function object

  implicit none
  
  class(sat_func_KRP1_type), pointer :: SF_KRP1_Create
  
  allocate(SF_KRP1_Create)
  call SF_KRP1_Create%Init()
  
end function SF_KRP1_Create

! ************************************************************************** !

subroutine SF_KRP1_Init(this)

  ! Creates the BRAGFLO KRP1 capillary pressure function object

  implicit none
  
  class(sat_func_KRP1_type) :: this

  call SFBaseInit(this)
  call SF_WIPP_Init(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%m = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SF_KRP1_Init

! ************************************************************************** !

subroutine SF_KRP1_Verify(this,name,option)

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
  call SF_WIPP_Verify(this,string,option)
  
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

end subroutine SF_KRP1_Verify

! ************************************************************************** !

subroutine SF_KRP1_CapillaryPressure(this,liquid_saturation, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP1_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
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

  call SF_WIPP_KPC(this,lambda,P0,Se2,capillary_pressure)

end subroutine SF_KRP1_CapillaryPressure

! ************************************************************************** !

subroutine SF_KRP1_Saturation(this,capillary_pressure, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP1_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
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
  
end subroutine SF_KRP1_Saturation

! ************************************************************************** !
! ************************************************************************** !

function SF_KRP2_Create()

  ! Creates the BRAGFLO KRP2 capillary pressure function object

  implicit none
  
  class(sat_func_KRP2_type), pointer :: SF_KRP2_Create
  
  allocate(SF_KRP2_Create)
  call SF_KRP2_Create%Init()
  
end function SF_KRP2_Create

! ************************************************************************** !

subroutine SF_KRP2_Init(this)

  ! Creates the BRAGFLO KRP2 capillary pressure function object

  implicit none
  
  class(sat_func_KRP2_type) :: this

  call SFBaseInit(this)
  call SF_WIPP_Init(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SF_KRP2_Init

! ************************************************************************** !

subroutine SF_KRP2_Verify(this,name,option)

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
  call SF_WIPP_Verify(this,string,option)
  
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call PrintErrMsg(option)
  endif    

end subroutine SF_KRP2_Verify

! ************************************************************************** !

subroutine SF_KRP2_CapillaryPressure(this,liquid_saturation, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP2_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
  endif
  
  Se1 = (liquid_saturation-this%Sr)/(1.d0-this%Sr)
      
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = 0.d0
  else 
    capillary_pressure = this%pct/(Se1**(1.d0/this%lambda))
  endif
  
  call SF_WIPP_KPC(this,this%lambda,this%pct,Se1,capillary_pressure)
  
end subroutine SF_KRP2_CapillaryPressure

! ************************************************************************** !

subroutine SF_KRP2_Saturation(this,capillary_pressure, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP2_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
  endif
  
  if (capillary_pressure < this%pct) then
    liquid_saturation = 1.d0
    dsat_dpres = 0.d0
    return
  endif
  
  Se1 = (capillary_pressure/this%pct)**(-1.d0*this%lambda)
  liquid_saturation = this%Sr + (1.d0-this%Sr)*Se1
  
end subroutine SF_KRP2_Saturation

! ************************************************************************** !
! ************************************************************************** !

function SF_KRP3_Create()

  ! Creates the BRAGFLO KRP3 capillary pressure function object

  implicit none
  
  class(sat_func_KRP3_type), pointer :: SF_KRP3_Create
  
  allocate(SF_KRP3_Create)
  call SF_KRP3_Create%Init()
  
end function SF_KRP3_Create

! ************************************************************************** !

subroutine SF_KRP3_Init(this)

  ! Creates the BRAGFLO KRP3 capillary pressure function object

  implicit none
  
  class(sat_func_KRP3_type) :: this

  call SFBaseInit(this)
  call SF_WIPP_Init(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SF_KRP3_Init

! ************************************************************************** !

subroutine SF_KRP3_Verify(this,name,option)

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
  call SF_WIPP_Verify(this,string,option)
  
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

end subroutine SF_KRP3_Verify

! ************************************************************************** !

subroutine SF_KRP3_CapillaryPressure(this,liquid_saturation, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP3_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
  endif
  
  Se2 = (liquid_saturation-this%Sr)/(1.d0-this%Sr-this%Srg)
  
  if ((1.d0-liquid_saturation) <= this%Srg) then
    capillary_pressure = this%pct
  elseif (liquid_saturation > this%Sr) then
    capillary_pressure = this%pct/(Se2**(1.d0/this%lambda))
  else
    capillary_pressure = 0.d0
  endif
  
  call SF_WIPP_KPC(this,this%lambda,this%pct,Se2,capillary_pressure)
  
end subroutine SF_KRP3_CapillaryPressure

! ************************************************************************** !

subroutine SF_KRP3_Saturation(this,capillary_pressure, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP3_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
  endif
  
  term = this%pct*((1.d0-this%Sr)/(1.d0-this%Sr-this%Srg))**(-1.d0/this%lambda)

  if ((capillary_pressure) < term) then
    liquid_saturation = 1.d0
    dsat_dpres = 0.d0
    return
  endif
  
  Se2 = (this%pct/capillary_pressure)**this%lambda
  liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se2
  
end subroutine SF_KRP3_Saturation

! ************************************************************************** !
! ************************************************************************** !

function SF_KRP4_Create()

  ! Creates the BRAGFLO KRP4 capillary pressure function object

  implicit none
  
  class(sat_func_KRP4_type), pointer :: SF_KRP4_Create
  
  allocate(SF_KRP4_Create)
  call SF_KRP4_Create%Init()
  
end function SF_KRP4_Create

! ************************************************************************** !

subroutine SF_KRP4_Init(this)

  ! Creates the BRAGFLO KRP4 capillary pressure function object

  implicit none
  
  class(sat_func_KRP4_type) :: this

  call SFBaseInit(this)
  call SF_WIPP_Init(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SF_KRP4_Init

! ************************************************************************** !

subroutine SF_KRP4_Verify(this,name,option)

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
  call SF_WIPP_Verify(this,string,option)
  
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

end subroutine SF_KRP4_Verify

! ************************************************************************** !

subroutine SF_KRP4_CapillaryPressure(this,liquid_saturation, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP4_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
  endif
  
  Se2 = (liquid_saturation-this%Sr)/(1.d0-this%Sr-this%Srg)
  
  capillary_pressure = 0.d0
  if ((1.d0-liquid_saturation) <= this%Srg) then
    capillary_pressure = this%pct/(Se2**(1.d0/this%lambda))
  endif
  if (liquid_saturation > this%Sr) then
    capillary_pressure = this%pct/(Se2**(1.d0/this%lambda))
  endif
  
  call SF_WIPP_KPC(this,this%lambda,this%pct,Se2,capillary_pressure)
  
end subroutine SF_KRP4_CapillaryPressure

! ************************************************************************** !

subroutine SF_KRP4_Saturation(this,capillary_pressure, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP1_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
  endif
  
  term = this%pct*((1.d0-this%Sr)/(1.d0-this%Sr-this%Srg))**(-1.d0/this%lambda)

  if ((capillary_pressure) < term) then
    liquid_saturation = 1.d0
    dsat_dpres = 0.d0
    return
  endif
  
  Se2 = (this%pct/capillary_pressure)**this%lambda
  liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se2
  
end subroutine SF_KRP4_Saturation

! ************************************************************************** !
! ************************************************************************** !

function SF_KRP5_Create()

  ! Creates the BRAGFLO KRP5 capillary pressure function object

  implicit none
  
  class(sat_func_KRP5_type), pointer :: SF_KRP5_Create
  
  allocate(SF_KRP5_Create)
  call SF_KRP5_Create%Init()
  
end function SF_KRP5_Create

! ************************************************************************** !

subroutine SF_KRP5_Init(this)

  ! Creates the BRAGFLO KRP5 capillary pressure function object

  implicit none
  
  class(sat_func_KRP5_type) :: this

  call SFBaseInit(this)
  call SF_WIPP_Init(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%pcmax = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SF_KRP5_Init

! ************************************************************************** !

subroutine SF_KRP5_Verify(this,name,option)

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
  call SF_WIPP_Verify(this,string,option)
  
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

end subroutine SF_KRP5_Verify

! ************************************************************************** !

subroutine SF_KRP5_CapillaryPressure(this,liquid_saturation, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP5_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
  endif
  
  Se2 = (liquid_saturation-this%Sr)/(1.d0-this%Sr-this%Srg)
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
  else if ((1.d0 - liquid_saturation) <= this%Srg) then
    capillary_pressure = this%pct
  else 
    capillary_pressure = (this%pct-this%pcmax)*Se2 + this%pcmax
  endif
  
end subroutine SF_KRP5_CapillaryPressure

! ************************************************************************** !

subroutine SF_KRP5_Saturation(this,capillary_pressure, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP5_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
  endif

  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  endif
  
  Se2 = (capillary_pressure-this%pcmax)/(this%pct-this%pcmax)
  liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se2

end subroutine SF_KRP5_Saturation
                            
! ************************************************************************** !
! ************************************************************************** !

function SF_KRP8_Create()

  ! Creates the BRAGFLO KRP8 capillary pressure function object

  implicit none
  
  class(sat_func_KRP8_type), pointer :: SF_KRP8_Create
  
  allocate(SF_KRP8_Create)
  call SF_KRP8_Create%Init()
  
end function SF_KRP8_Create

! ************************************************************************** !

subroutine SF_KRP8_Init(this)

  ! Creates the BRAGFLO KRP8 capillary pressure function object

  implicit none
  
  class(sat_func_KRP8_type) :: this

  call SFBaseInit(this)
  call SF_WIPP_Init(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%m = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SF_KRP8_Init

! ************************************************************************** !

subroutine SF_KRP8_Verify(this,name,option)

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
  call SF_WIPP_Verify(this,string,option)
  
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

end subroutine SF_KRP8_Verify

! ************************************************************************** !

subroutine SF_KRP8_CapillaryPressure(this,liquid_saturation, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP8_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
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

  call SF_WIPP_KPC(this,lambda,P0,Se1,capillary_pressure)
  
end subroutine SF_KRP8_CapillaryPressure

! ************************************************************************** !

subroutine SF_KRP8_Saturation(this,capillary_pressure, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP8_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
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
  
end subroutine SF_KRP8_Saturation

! ************************************************************************** !
! ************************************************************************** !

function SF_KRP9_Create()

  ! Creates the BRAGFLO KRP9 capillary pressure function object

  implicit none
  
  class(sat_func_KRP9_type), pointer :: SF_KRP9_Create
  
  allocate(SF_KRP9_Create)
  call SF_KRP9_Create%Init()
  
end function SF_KRP9_Create

! ************************************************************************** !

subroutine SF_KRP9_Init(this)

  ! Creates the BRAGFLO KRP9 capillary pressure function object

  implicit none
  
  class(sat_func_KRP9_type) :: this

  call SFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine SF_KRP9_Init

! ************************************************************************** !

subroutine SF_KRP9_Verify(this,name,option)

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

end subroutine SF_KRP9_Verify

! ************************************************************************** !

subroutine SF_KRP9_CapillaryPressure(this,liquid_saturation, &
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
  
end subroutine SF_KRP9_CapillaryPressure

! ************************************************************************** !

subroutine SF_KRP9_Saturation(this,capillary_pressure, &
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

end subroutine SF_KRP9_Saturation

! ************************************************************************** !
! ************************************************************************** !

function SF_KRP11_Create()

  ! Creates the BRAGFLO KRP11 capillary pressure function object

  implicit none
  
  class(sat_func_KRP11_type), pointer :: SF_KRP11_Create
  
  allocate(SF_KRP11_Create)
  call SF_KRP11_Create%Init()
  
end function SF_KRP11_Create

! ************************************************************************** !

subroutine SF_KRP11_Init(this)

  ! Creates the BRAGFLO KRP11 capillary pressure function object

  implicit none
  
  class(sat_func_KRP11_type) :: this

  call SFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine SF_KRP11_Init

! ************************************************************************** !

subroutine SF_KRP11_Verify(this,name,option)

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

end subroutine SF_KRP11_Verify

! ************************************************************************** !

subroutine SF_KRP11_CapillaryPressure(this,liquid_saturation, &
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
  
end subroutine SF_KRP11_CapillaryPressure

! ************************************************************************** !

subroutine SF_KRP11_Saturation(this,capillary_pressure, &
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

end subroutine SF_KRP11_Saturation

! ************************************************************************** !
! ************************************************************************** !

function SF_KRP12_Create()

  ! Creates the BRAGFLO KRP12 capillary pressure function object

  implicit none
  
  class(sat_func_KRP12_type), pointer :: SF_KRP12_Create
  
  allocate(SF_KRP12_Create)
  call SF_KRP12_Create%Init()
  
end function SF_KRP12_Create

! ************************************************************************** !

subroutine SF_KRP12_Init(this)

  ! Creates the BRAGFLO KRP12 capillary pressure function object

  implicit none
  
  class(sat_func_KRP12_type) :: this

  call SFBaseInit(this)
  call SF_WIPP_Init(this)
  this%lambda = UNINITIALIZED_DOUBLE
  this%s_min = UNINITIALIZED_DOUBLE
  this%s_effmin = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SF_KRP12_Init

! ************************************************************************** !

subroutine SF_KRP12_Verify(this,name,option)

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
  call SF_WIPP_Verify(this,string,option)
 
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
  
end subroutine SF_KRP12_Verify

! ************************************************************************** !

subroutine SF_KRP12_CapillaryPressure(this,liquid_saturation, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP12_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
  endif

  Se = (liquid_saturation - (this%s_min - this%s_effmin)) / &
       (1.d0 - (this%s_min - this%s_effmin))
  Se21 = max(min(Se,1.d0),this%s_effmin)
  
  capillary_pressure = this%pct/(Se21**(1.d0/this%lambda))
  
  ! do not pass in Se21 into the following function, it needs Se1:
  Se1 = (liquid_saturation - this%Sr)/(1.d0 - this%Sr)
  call SF_WIPP_KPC(this,this%lambda,this%pct,Se1,capillary_pressure)

end subroutine SF_KRP12_CapillaryPressure

! ************************************************************************** !

subroutine SF_KRP12_Saturation(this,capillary_pressure, &
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
    if (.not.option%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sat_func_KRP12_type. STOPPING.'
      call PrintErrMsg(option)
    endif
    option%pct_updated = PETSC_FALSE
  endif
  
  if (capillary_pressure < this%pct) then
    liquid_saturation = 1.d0
    dsat_dpres = 0.d0
    return
  endif
  
  Se21 = (this%pct/capillary_pressure)**this%lambda
  liquid_saturation = Se21*(1.d0 - this%s_min - this%s_effmin) + &
                      this%s_min + this%s_effmin
                               
end subroutine SF_KRP12_Saturation

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP1_Liq_Create()

  ! Creates the BRAGFLO_KRP1_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP1_liq_type), pointer :: RPF_KRP1_Liq_Create
  
  allocate(RPF_KRP1_Liq_Create)
  call RPF_KRP1_Liq_Create%Init()
  
end function RPF_KRP1_Liq_Create

! ************************************************************************** !

subroutine RPF_KRP1_Liq_Init(this)

  ! Initializes the BRAGFLO_KRP1_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP1_liq_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP1_Liq_Init

! ************************************************************************** !

subroutine RPF_KRP1_Liq_Verify(this,name,option)

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
  
end subroutine RPF_KRP1_Liq_Verify

! ************************************************************************** !

subroutine RPF_KRP1_Liq_RelPerm(this,liquid_saturation, &
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
  
end subroutine RPF_KRP1_Liq_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP1_Gas_Create()

  ! Creates the BRAGFLO_KRP1_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP1_gas_type), pointer :: RPF_KRP1_Gas_Create
  
  allocate(RPF_KRP1_Gas_Create)
  call RPF_KRP1_Gas_Create%Init() 
  
end function RPF_KRP1_Gas_Create

! ************************************************************************** !

subroutine RPF_KRP1_Gas_Init(this)

  ! Initializes the BRAGFLO_KRP1_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP1_gas_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP1_Gas_Init

! ************************************************************************** !

subroutine RPF_KRP1_Gas_Verify(this,name,option)

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
  
end subroutine RPF_KRP1_Gas_Verify

! ************************************************************************** !

subroutine RPF_KRP1_Gas_RelPerm(this,liquid_saturation, &
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
  
end subroutine RPF_KRP1_Gas_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP2_Liq_Create()

  ! Creates the BRAGFLO_KRP2_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP2_liq_type), pointer :: RPF_KRP2_Liq_Create
  
  allocate(RPF_KRP2_Liq_Create)
  call RPF_KRP2_Liq_Create%Init()
  
end function RPF_KRP2_Liq_Create

! ************************************************************************** !

subroutine RPF_KRP2_Liq_Init(this)

  ! Initializes the BRAGFLO_KRP2_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP2_liq_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP2_Liq_Init

! ************************************************************************** !

subroutine RPF_KRP2_Liq_Verify(this,name,option)

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
  
end subroutine RPF_KRP2_Liq_Verify

! ************************************************************************** !

subroutine RPF_KRP2_Liq_RelPerm(this,liquid_saturation, &
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
  
end subroutine RPF_KRP2_Liq_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP2_Gas_Create()

  ! Creates the BRAGFLO_KRP2_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP2_gas_type), pointer :: RPF_KRP2_Gas_Create
  
  allocate(RPF_KRP2_Gas_Create)
  call RPF_KRP2_Gas_Create%Init()
  
end function RPF_KRP2_Gas_Create

! ************************************************************************** !

subroutine RPF_KRP2_Gas_Init(this)

  ! Initializes the BRAGFLO_KRP2_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP2_gas_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP2_Gas_Init

! ************************************************************************** !

subroutine RPF_KRP2_Gas_Verify(this,name,option)

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
  
end subroutine RPF_KRP2_Gas_Verify

! ************************************************************************** !

subroutine RPF_KRP2_Gas_RelPerm(this,liquid_saturation, &
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
    
end subroutine RPF_KRP2_Gas_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP3_Liq_Create()

  ! Creates the BRAGFLO_KRP3_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP3_liq_type), pointer :: RPF_KRP3_Liq_Create
  
  allocate(RPF_KRP3_Liq_Create)
  call RPF_KRP3_Liq_Create%Init() 
  
end function RPF_KRP3_Liq_Create

! ************************************************************************** !

subroutine RPF_KRP3_Liq_Init(this)

  ! Initializes the BRAGFLO_KRP3_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP3_liq_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP3_Liq_Init

! ************************************************************************** !

subroutine RPF_KRP3_Liq_Verify(this,name,option)

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
  
end subroutine RPF_KRP3_Liq_Verify

! ************************************************************************** !

subroutine RPF_KRP3_Liq_RelPerm(this,liquid_saturation, &
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
  
end subroutine RPF_KRP3_Liq_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP3_Gas_Create()

  ! Creates the BRAGFLO_KRP3_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP3_gas_type), pointer :: RPF_KRP3_Gas_Create
  
  allocate(RPF_KRP3_Gas_Create)
  call RPF_KRP3_Gas_Create%Init() 
  
end function RPF_KRP3_Gas_Create

! ************************************************************************** !

subroutine RPF_KRP3_Gas_Init(this)

  ! Initializes the BRAGFLO_KRP3_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP3_gas_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP3_Gas_Init

! ************************************************************************** !

subroutine RPF_KRP3_Gas_Verify(this,name,option)

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
  
end subroutine RPF_KRP3_Gas_Verify

! ************************************************************************** !

subroutine RPF_KRP3_Gas_RelPerm(this,liquid_saturation, &
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
   
end subroutine RPF_KRP3_Gas_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP4_Liq_Create()

  ! Creates the BRAGFLO_KRP4_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP4_liq_type), pointer :: RPF_KRP4_Liq_Create
  
  allocate(RPF_KRP4_Liq_Create)
  call RPF_KRP4_Liq_Create%Init() ! Calls KRP3_Liq's Init()
  
end function RPF_KRP4_Liq_Create

! ************************************************************************** !

subroutine RPF_KRP4_Liq_Verify(this,name,option)

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
  
end subroutine RPF_KRP4_Liq_Verify

! ************************************************************************** !

subroutine RPF_KRP4_Liq_RelPerm(this,liquid_saturation, &
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
  
end subroutine RPF_KRP4_Liq_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP4_Gas_Create()

  ! Creates the BRAGFLO_KRP4_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP4_gas_type), pointer :: RPF_KRP4_Gas_Create
  
  allocate(RPF_KRP4_Gas_Create)
  call RPF_KRP4_Gas_Create%Init() ! calls KRP3_Gas's Init()
  
end function RPF_KRP4_Gas_Create

! ************************************************************************** !

subroutine RPF_KRP4_Gas_Verify(this,name,option)

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
  
end subroutine RPF_KRP4_Gas_Verify

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP5_Liq_Create()

  ! Creates the BRAGFLO_KRP5_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP5_liq_type), pointer :: RPF_KRP5_Liq_Create
  
  allocate(RPF_KRP5_Liq_Create)
  call RPF_KRP5_Liq_Create%Init()
  
end function RPF_KRP5_Liq_Create

! ************************************************************************** !

subroutine RPF_KRP5_Liq_Init(this)

  ! Initializes the BRAGFLO_KRP5_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP5_liq_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP5_Liq_Init

! ************************************************************************** !

subroutine RPF_KRP5_Liq_Verify(this,name,option)

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
    
end subroutine RPF_KRP5_Liq_Verify

! ************************************************************************** !

subroutine RPF_KRP5_Liq_RelPerm(this,liquid_saturation, &
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
   
end subroutine RPF_KRP5_Liq_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP5_Gas_Create()

  ! Creates the BRAGFLO_KRP5_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP5_gas_type), pointer :: RPF_KRP5_Gas_Create
  
  allocate(RPF_KRP5_Gas_Create)
  call RPF_KRP5_Gas_Create%Init()
  
end function RPF_KRP5_Gas_Create

! ************************************************************************** !

subroutine RPF_KRP5_Gas_Init(this)

  ! Initializes the BRAGFLO_KRP5_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP5_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP5_Gas_Init

! ************************************************************************** !

subroutine RPF_KRP5_Gas_Verify(this,name,option)

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
    
end subroutine RPF_KRP5_Gas_Verify

! ************************************************************************** !

subroutine RPF_KRP5_Gas_RelPerm(this,liquid_saturation, &
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
   
end subroutine RPF_KRP5_Gas_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP8_Liq_Create()

  ! Creates the BRAGFLO_KRP8_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP8_liq_type), pointer :: RPF_KRP8_Liq_Create
  
  allocate(RPF_KRP8_Liq_Create)
  call RPF_KRP8_Liq_Create%Init()
  
end function RPF_KRP8_Liq_Create

! ************************************************************************** !

subroutine RPF_KRP8_Liq_Init(this)

  ! Initializes the BRAGFLO_KRP8_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP8_liq_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP8_Liq_Init

! ************************************************************************** !

subroutine RPF_KRP8_Liq_Verify(this,name,option)

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
  
end subroutine RPF_KRP8_Liq_Verify

! ************************************************************************** !

subroutine RPF_KRP8_Liq_RelPerm(this,liquid_saturation, &
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
  
end subroutine RPF_KRP8_Liq_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP8_Gas_Create()

  ! Creates the BRAGFLO_KRP8_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP8_gas_type), pointer :: RPF_KRP8_Gas_Create
  
  allocate(RPF_KRP8_Gas_Create)
  call RPF_KRP8_Gas_Create%Init()
  
end function RPF_KRP8_Gas_Create

! ************************************************************************** !

subroutine RPF_KRP8_Gas_Init(this)

  ! Initializes the BRAGFLO_KRP8_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP8_gas_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP8_Gas_Init

! ************************************************************************** !

subroutine RPF_KRP8_Gas_Verify(this,name,option)

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
  
end subroutine RPF_KRP8_Gas_Verify

! ************************************************************************** !

subroutine RPF_KRP8_Gas_RelPerm(this,liquid_saturation, &
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
  
end subroutine RPF_KRP8_Gas_RelPerm

! ************************************************************************** !
! ************************************************************************** !
                                     
function RPF_KRP9_Liq_Create()

  ! Creates the BRAGFLO_KRP9_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP9_liq_type), pointer :: RPF_KRP9_Liq_Create
  
  allocate(RPF_KRP9_Liq_Create)
  call RPF_KRP9_Liq_Create%Init()
  
end function RPF_KRP9_Liq_Create

! ************************************************************************** !

subroutine RPF_KRP9_Liq_Init(this)

  ! Initializes the BRAGFLO_KRP9_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP9_liq_type) :: this

  call RPFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP9_Liq_Init

! ************************************************************************** !

subroutine RPF_KRP9_Liq_Verify(this,name,option)

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
  
end subroutine RPF_KRP9_Liq_Verify

! ************************************************************************** !

subroutine RPF_KRP9_Liq_RelPerm(this,liquid_saturation, &
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
  
end subroutine RPF_KRP9_Liq_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP9_Gas_Create()

  ! Creates the BRAGFLO_KRP9_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP9_gas_type), pointer :: RPF_KRP9_Gas_Create
  
  allocate(RPF_KRP9_Gas_Create)
  call RPF_KRP9_Gas_Create%Init()
  
end function RPF_KRP9_Gas_Create

! ************************************************************************** !

subroutine RPF_KRP9_Gas_Init(this)

  ! Initializes the BRAGFLO_KRP9_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP9_gas_type) :: this

  call RPFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP9_Gas_Init

! ************************************************************************** !

subroutine RPF_KRP9_Gas_Verify(this,name,option)

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
  
end subroutine RPF_KRP9_Gas_Verify

! ************************************************************************** !

subroutine RPF_KRP9_Gas_RelPerm(this,liquid_saturation, &
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
  
  call RPF_KRP9_Liq_RelPerm(this,liquid_saturation, &
                            liquid_relative_permeability, &
                            liquid_dkr_sat,option)
  
  relative_permeability = 1.d0 - liquid_relative_permeability
  dkr_sat = -1.d0 * liquid_dkr_sat
  
end subroutine RPF_KRP9_Gas_RelPerm

! ************************************************************************** !
! ************************************************************************** !
 
function RPF_KRP11_Liq_Create()

  ! Creates the BRAGFLO_KRP11_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP11_liq_type), pointer :: RPF_KRP11_Liq_Create
  
  allocate(RPF_KRP11_Liq_Create)
  call RPF_KRP11_Liq_Create%Init()
  
end function RPF_KRP11_Liq_Create

! ************************************************************************** !

subroutine RPF_KRP11_Liq_Init(this)

  ! Initializes the BRAGFLO_KRP11_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP11_liq_type) :: this

  call RPFBaseInit(this)
  this%tolc = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_KRP11_Liq_Init

! ************************************************************************** !

subroutine RPF_KRP11_Liq_Verify(this,name,option)

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

end subroutine RPF_KRP11_Liq_Verify

! ************************************************************************** !

subroutine RPF_KRP11_Liq_RelPerm(this,liquid_saturation, &
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
    
end subroutine RPF_KRP11_Liq_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP11_Gas_Create()

  ! Creates the BRAGFLO_KRP11_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP11_gas_type), pointer :: RPF_KRP11_Gas_Create
  
  allocate(RPF_KRP11_Gas_Create)
  call RPF_KRP11_Gas_Create%Init() ! calls KRP11_LIQ's Init()
  
end function RPF_KRP11_Gas_Create

! ************************************************************************** !

subroutine RPF_KRP11_Gas_Verify(this,name,option)

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

end subroutine RPF_KRP11_Gas_Verify

! ************************************************************************** !

subroutine RPF_KRP11_Gas_RelPerm(this,liquid_saturation, &
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
  
  end subroutine RPF_KRP11_Gas_RelPerm

! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP12_Liq_Create()

  ! Creates the BRAGFLO_KRP12_LIQ relative permeability function object

  implicit none
  
  class(rpf_KRP12_liq_type), pointer :: RPF_KRP12_Liq_Create
  
  allocate(RPF_KRP12_Liq_Create)
  call RPF_KRP12_Liq_Create%Init()
  
end function RPF_KRP12_Liq_Create

! ************************************************************************** !

subroutine RPF_KRP12_Liq_Verify(this,name,option)

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

end subroutine RPF_KRP12_Liq_Verify

! ************************************************************************** !

subroutine RPF_KRP12_Liq_RelPerm(this,liquid_saturation, &
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
  
end subroutine RPF_KRP12_Liq_RelPerm
  
! ************************************************************************** !
! ************************************************************************** !

function RPF_KRP12_Gas_Create()

  ! Creates the BRAGFLO_KRP12_GAS relative permeability function object

  implicit none
  
  class(rpf_KRP12_gas_type), pointer :: RPF_KRP12_Gas_Create
  
  allocate(RPF_KRP12_Gas_Create)
  call RPF_KRP12_Gas_Create%Init()
  
end function RPF_KRP12_Gas_Create

! ************************************************************************** !

subroutine RPF_KRP12_Gas_Verify(this,name,option)

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

end subroutine RPF_KRP12_Gas_Verify

! ************************************************************************** !

subroutine RPF_KRP12_Gas_RelPerm(this,liquid_saturation, &
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
  
end subroutine RPF_KRP12_Gas_RelPerm
  
! ************************************************************************** !
! ************************************************************************** !

function RPF_TOUGH2_IRP7_Gas_Create()

  ! Creates the Brooks-Corey Burdine gas relative permeability function object

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type), pointer :: RPF_TOUGH2_IRP7_Gas_Create
  
  allocate(RPF_TOUGH2_IRP7_Gas_Create)
  call RPF_TOUGH2_IRP7_Gas_Create%Init()
  
end function RPF_TOUGH2_IRP7_Gas_Create

! ************************************************************************** !

subroutine RPF_TOUGH2_IRP7_Gas_Init(this)

  ! Initializes the Brooks-Corey Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_TOUGH2_IRP7_Gas_Init

! ************************************************************************** !

subroutine RPF_TOUGH2_IRP7_Gas_Verify(this,name,option)

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
  
end subroutine RPF_TOUGH2_IRP7_Gas_Verify

! ************************************************************************** !

subroutine RPF_TOUGH2_IRP7_Gas_RelPerm(this,liquid_saturation, &
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
  dkr_sat = dkr_sat / 0.d0
  dkr_sat = dkr_sat * 0.d0
  
                 ! essentially zero
  if (this%Srg <= 0.d0) then
    call RPF_Mualem_VG_Liq_RelPerm(this,liquid_saturation, &
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
    
end subroutine RPF_TOUGH2_IRP7_Gas_RelPerm

end module Characteristic_Curves_WIPP_module
