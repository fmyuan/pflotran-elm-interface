module characteristic_curves_WIPP_invariant_module
#include "petsc/finclude/petscsys.h"

use petscsys ! Necessary for PETSC_TRUE / PETSC_FALSE
use option_module ! Necessary for characterisic curve base class arguments
use characteristic_curves_base_module ! Needed to define base type
implicit none

! **************************************************************************** !
! van Genuchten vs Brooks-Corey enumerations
! **************************************************************************** !
! Not applicable to KRP 6,7,9,10,11
#define KRP_VG_range      1,        8
#define KRP_BC_range        2,3,4,       12
#define KRP_lin_range             5,   11
! KRP 4 is peculiar in that effective saturation for Pc may exceed unity
#define KRP_BC1_range       2,3,         12
#define KRP_BC2_range           4

! **************************************************************************** !
! Effective saturation enumerations
! **************************************************************************** !
! Se1 has only a liquid residual (Sgr = 0)
#define KRP_Se1_Pc_range    2,      8
#define KRP_Se1_Krw_range 1,2,  4,  8,   12
#define KRP_Se1_Krg_range   2,      8

! Se2 has both liquid and gas residuals
#define KRP_Se2_Pc_range  1,  3,4,5
#define KRP_Se2_Krw_range     3,  5,  11
#define KRP_Se2_Krg_range 1,  3,4,5,  11,12

! **************************************************************************** !
! Common WIPP Characteristic Curve Types
! **************************************************************************** !

type, public, extends(sat_func_base_type) :: sf_WIPP_type
  private
    procedure(set_k_type), public, pointer :: setK
    procedure(set_Swj_type), pointer :: setSwj

    procedure(calc_Pc_type), pointer :: KPCPc
    procedure(calc_Sw_type), pointer :: KPCSw

    procedure(calc_Pc_type), pointer :: KRPPc
    procedure(calc_Sw_type), pointer :: KRPSw

!   PFLOTRAN object parameter             BRAGFLO equivalent

!   Effective Saturation Parameters
    PetscReal :: Swr                      ! SWR or SOCZRO
    PetscReal :: Sgr_comp                 ! 1 - SGR
    PetscReal :: Sw_span                  ! 1/SETERM or 1/SETERM2
    PetscReal :: dSe_dSw                  ! SETERM or SETERM2
    PetscReal :: Semin                    ! SOCEFFMIN for KRP12

!   Brooks-Corey Parameters
    PetscReal :: lambda                   ! 1/XLAM1
    PetscReal :: lambda_nrec              ! -XLAM1

!   Van Genuchten Parameters
    PetscReal :: m                        ! XLAM4
    PetscReal :: m_nrec                   ! XLAM6
    PetscReal :: m_comp                   ! XLAM7
    PetscReal :: n                        ! 1/XLAM7

!   Coefficient for analytical derivatives
    PetscReal :: k_dSe_dSw                ! -m*n*dSe_dSw or -dSe_dSw/lambda

!   Unsaturated Extension Parameters
!   PetscReal :: Pcmax                    ! PCFIX          Defined in base
    PetscReal :: Swj, Pcj, dPcj_dSwj

!   Threshold Pressure Parameters
    PetscReal :: Pct                      ! Calculated Pct
    PetscReal :: Pct_a
    PetscReal :: Pct_exp
    PetscReal :: permeability
    PetscReal :: Pcm_Pct                  ! Cached PCM:PCT ratio
  contains
! Overridden function pointers from the PFLOTRAN base class
    procedure, public :: CapillaryPressure => SFWIPPCapillaryPressure
    procedure, public :: Saturation        => SFWIPPSaturation
end type

! **************************************************************************** !

type, public, extends(rel_perm_func_base_type) :: rpf_WIPP_type
  private
    procedure(calc_Kr_type), pointer :: KRPKr

!   Effective Saturation Parameters
    PetscReal :: Swr                      ! SWR
    PetscReal :: Sgr_comp                 ! 1 - SGR
    PetscReal :: Sw_span                  ! TOLC*(1-SWR-SGR)
    PetscReal :: dSe_dSw                  ! SETERM or SETERM2

!   Brooks-Corey / van Genuchten Parameters
!   PetscReal :: lambda                   ! 1/XLAM1
    PetscReal :: m                        ! XLAM4
    PetscReal :: m_rec                    ! -XLAM6
    PetscReal :: expon                    ! 3+2/lambda OR 1+2/lambda OR 2*m

!   Cached Saturated Limits
    PetscReal :: Kr_Swr, dKr_dSwr         ! Residual water limit values
    PetscReal :: Kr_Sgr, dKr_dSgr         ! Residual   gas limit values
  contains
! Overridden function pointers from the PFLOTRAN base class
    procedure, public :: RelativePermeability => RPFWIPPRelativePermeability
end type

! **************************************************************************** !
! Function prototypes for function pointers
! **************************************************************************** !

abstract interface
  function set_Swj_type(this, Swj) result (error)
    import sf_WIPP_type
    class(sf_WIPP_type), intent(inout) :: this
    PetscReal, intent(in)  :: Swj
    PetscInt :: error
  end function
  pure subroutine calc_Pc_type(this, Sw, Pc, dPc_dSw)
    import sf_WIPP_type
    class(sf_WIPP_type), intent(in) :: this
    PetscReal, intent(in)  :: Sw
    PetscReal, intent(out) :: Pc, dPc_dSw
  end subroutine
  pure subroutine calc_Kr_type(this, Sw, Kr, dKr_dSw)
    import rpf_WIPP_type
    class(rpf_WIPP_type), intent(in) :: this
    PetscReal, intent(in)  :: Sw
    PetscReal, intent(out) :: Kr, dKr_dSw
  end subroutine
  pure subroutine calc_Sw_type(this, Pc, Sw)
    import sf_WIPP_type
    class(sf_WIPP_type), intent(in) :: this
    PetscReal, intent(in)  :: Pc
    PetscReal, intent(out) :: Sw
  end subroutine
  subroutine set_k_type(this, k)
    import sf_WIPP_type
    class(sf_WIPP_type), intent(inout) :: this
    PetscReal, intent(in)  :: k
  end subroutine
end interface

! **************************************************************************** !
! Public and Private Procedure Declarations
! **************************************************************************** !

! WIPP constructors
public  :: SFWIPPctor, RPFWIPPctor

! Implemented WIPP KPC Pc procedures
private :: SFWIPPKPC1Pc , SFWIPPKPC1Sw , SFWIPPKPC1Swj, &
           SFWIPPKPC2Pc , SFWIPPKPC2Sw , SFWIPPKPC2Swj, &
           SFWIPPKPC6Pc , SFWIPPKPC6Sw , SFWIPPKPC6Swj, &
                                         SFWIPPKRP12Swj

! Implemented WIPP KRP Pc procedures
private :: SFWIPPVGPc    , SFWIPPVGSw , &    ! van Genuchten
           SFWIPPBCPc    , SFWIPPBCSw , &    ! Brooks-Corey
           SFWIPPKRP4Pc  , SFWIPPKRP4Sw , &  ! Extended Brooks-Corey
           SFWIPPKRP5Pc  , SFWIPPKRP5Sw , &
           SFWIPPKRP9Pc  , SFWIPPKRP9Sw , &
           SFWIPPKRP11Pc , SFWIPPKRP11Sw

! Implemented WIPP Permeability Pct procedures
private :: SFWIPPSetK, SFWIPPIgnoreK

! Implemented WIPP KRP Kr procedures
private :: RPFWIPPMVGKrw , RPFWIPPMVGKrg,  & ! Mualem - van Genuchten
           RPFWIPPBBCKrw , RPFWIPPBBCKrg,  & ! Burdine - Brooks-Corey
           RPFWIPPlinKrw , RPFWIPPlinKrg,  & ! Linear
           RPFWIPPKRP9Krw  , RPFWIPPKRP9Krg

contains 

! **************************************************************************** !
! WIPP Constructors
! **************************************************************************** !

function SFWIPPctor(KRP, KPC, Swr, Sgr, expon, Pct_ignore, Pct_alpha, &
                    Pct_expon, Pcmax, Swj) result (new)
  class(sf_WIPP_type), pointer :: new
  PetscInt, intent(in)  :: KRP, KPC
  PetscReal, intent(in) :: Swr, Sgr, expon, Pct_alpha, Pct_expon, Pcmax
  PetscReal, intent(inout) ::  Swj ! KRP12 sends Semin and returns Swj
  PetscBool, intent(in) :: Pct_ignore
  PetscInt :: error

  ! Memory allocation
  allocate(new)
  if (.not. associated(new)) return ! Memory allocation failed, abort

  ! Derivatives have been defined to enable smooth unsaturated extensions
  new%analytical_derivative_available = PETSC_TRUE

  ! Data validation
                                        error = 0
  !  For KRP 11 - the cavity model - Pc is always zero.
  !  Except for KRP 9, if PCT_A is zero, Pc is also always zero.
  if (KRP == 11 .OR. (KRP /= 9 .AND. Pct_alpha == 0d0)) then
    new%setK   => SFWIPPIgnoreK
    new%setSwj => SFWIPPKPC1Swj

    new%KRPPc  => SFWIPPKRP11Pc
    new%KRPSw  => SFWIPPKRP11Sw

    new%KPCPc  => new%KRPPc
    new%KPCSw  => new%KRPSW
    return ! No need to check anything else, the Pc/Sw are are fully defined
  end if

  ! If KRP is in branch table, set function pointers, else flag error
  select case(KRP)
  case (KRP_VG_range)
    new%KRPPc => SFWIPPVGPc
    new%KRPSw => SFWIPPVGSw
  case (KRP_BC1_range)
    new%KRPPc => SFWIPPBCPc
    new%KRPSw => SFWIPPBCSw
  case (KRP_BC2_range)
    new%KRPPc => SFWIPPKRP4Pc
    new%KRPSw => SFWIPPKRP4Sw
  case (5)
    new%KRPPc => SFWIPPKRP5Pc
    new%KRPSw => SFWIPPKRP5Sw
  case (9)
    new%KRPPc => SFWIPPKRP9Pc
    new%KRPSw => SFWIPPKRP9Sw
  case default
                                        error = error + 1
  end select

  ! If KPC is in branch table, set function pointers, else flag error
  select case(KPC)
  case (1) ! 0 at or below residual
    new%KPCPc  => SFWIPPKPC1Pc 
    new%KPCSw  => SFWIPPKPC1Sw
    new%setSwj => SFWIPPKPC1Swj
  case (2) ! Pcmax at or below residual
    new%KPCPc  => SFWIPPKPC2Pc
    new%KPCSw  => SFWIPPKPC2Sw
    new%setSwj => SFWIPPKPC2Swj
  case (6) ! Linear at or below junction
    new%KPCPc  => SFWIPPKPC6Pc
    new%KPCSw  => SFWIPPKPC6Sw
    new%setSwj => SFWIPPKPC6Swj
  case default
                                        error = error + 2
  end select

  ! Check the residual saturations, and their sum, are bound between 0 and 1
  if (Swr < 0d0 .or. Swr >= 1d0)        error = error + 4
  if (Sgr < 0d0 .or. Sgr >= 1d0)        error = error + 8
  if (Swr + Sgr > 1d0)                  error = error + 12

  ! Check the exponent parameter is valid for the KRP option chosen
  select case(KRP)
  case (KRP_BC_range) ! Brooks-Corey types
    if (expon <= 0d0)                   error = error + 16
  case (KRP_VG_range) ! Van Genuchten types
    if (expon <= 0d0 .or. expon >= 1d0) error = error + 16
  case default ! Other KRP functions do not use lambda or m
  end select

  ! Abort if errors caught in data validation
  if (error /= 0) then
    ! TODO: Print error code(s) to error stream
    deallocate(new)
    nullify(new)
    return
  end if

  new%Pcmax = Pcmax

  ! Assign residual saturation parameters
  new%Swr = Swr
  select case(KRP)
  case (KRP_Se1_Pc_range)
    new%Sgr_comp = 1d0
  case (KRP_Se2_Pc_range)
    new%Sgr_comp = 1d0 - Sgr
  case (12) ! KRP 12 is peculiar
    new%Sgr_comp = 1d0
    new%Swr = Swr - Swj ! SOCZRO = SOCMIN - SOCEFFMIN
  ! Warning: BRAGFLO UM 6.02 indicates +, but BRAGFLO code is -
  ! I.e. 1 - (Smin - Seffmin) /= 1 - Smin - Seffmin
  ! Matching the behaivor of the code, not the documentation
  end select
  new%Sw_span = new%Sgr_comp - new%Swr 
  new%dSe_dSw = 1d0 / new%Sw_span

  ! Assign exponential BC and VG parameters
  select case(KRP)
  case (KRP_BC_range) ! Brooks-Corey types
    new%lambda = expon
    new%lambda_nrec = -1d0/new%lambda
    new%k_dSe_dSw = -new%dSe_dSw / new%lambda
    new%pcm_pct = 0d0
  case (KRP_VG_range) ! Van Genuchten types
  ! In BRAGFLO, VG functions use the closed form Mualem condition n = 1/(1-m)
    new%m = expon
    new%m_nrec = -1d0 / new%m
    new%m_comp =  1d0 - new%m
    new%n      =  1d0 / new%m_comp
    new%k_dSe_dSw = new%m_comp / (new%m * new%Sw_span)
    new%pcm_pct = 2d0**(1d0/new%m-1d0)
  end select

  ! Assign Pct/Pcm model parameters
  new%permeability = -1d0   ! Initial value ensures Pct calculation on first pass
  if (pct_ignore) then
    new%setK    => SFWIPPIgnoreK
    new%pct_a   = 0d0       ! Not used
    new%pct_exp = 0d0       ! Not used
    new%pct     = pct_alpha ! Pct permanently set to alpha
  else
    new%setK    => SFWIPPSetK
    new%pct_a   = pct_alpha
    new%pct_exp = pct_expon
    new%pct     = 1E6       ! Arbitrary 1 MPa to enable test output
  end if

  if (KRP == 12) then ! KRP 12 is peculiar, superceding KPC options
    Swj = Swj * new%Sw_span + new%Swr
    new%Swj = Swj
    new%KPCPc  => SFWIPPKPC2Pc
    new%KPCSw  => SFWIPPKPC2Sw
    new%setSwj => SFWIPPKRP12Swj
    ! TODO print warning or error if KRP == 12 and KPC /= 1
  end if

  ! Initialize unsaturated extensions 
  error = new%setSwj(Swj)
  ! Abort if error in unsaturated extension (e.g. Swj < Swr)
  if (error /= 0) then
    ! TODO: Print error code(s) to error stream
    deallocate(new)
    nullify(new)
    return
  end if
end function

! **************************************************************************** !
! WIPP Characteristic Curve Roots
! **************************************************************************** !

subroutine SFWIPPCapillaryPressure(this, liquid_saturation, capillary_pressure,&
                                         dpc_dsatl, option)
  class(sf_WIPP_type)              :: this
  PetscReal, intent(in)            :: liquid_saturation
  PetscReal, intent(out)           :: capillary_pressure, dpc_dsatl
  type(option_type), intent(inout) :: option

  if (liquid_saturation <= this%Swj) then
    call this%KPCPc(liquid_saturation, capillary_pressure, dpc_dsatl)
  else
    call this%KRPPc(liquid_saturation, capillary_pressure, dpc_dsatl)
  end if
end subroutine

! **************************************************************************** !

subroutine SFWIPPSaturation(this, capillary_pressure, liquid_saturation, &
                                  dsat_dpres, option)
  class(sf_WIPP_type)              :: this
  PetscReal, intent(in)            :: capillary_pressure
  PetscReal, intent(out)           :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  if (capillary_pressure >= this%Pcj) then
    call this%KPCSw(capillary_pressure, liquid_saturation)
  else
    call this%KRPSw(capillary_pressure, liquid_saturation)
  end if
  dsat_dpres = 0d0 ! TODO analytic derivatives for Richard's mode
end subroutine

! **************************************************************************** !
! WIPP PCT Subroutines
! **************************************************************************** !

subroutine SFWIPPSetK(this, permeability)
  class(sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: permeability
  PetscInt :: error
 
  ! Permeability changes with material region, but because this occurs in 
  ! blocks, frequently permeability has not changed from call to call.
  ! Because exponentiaton is expensive, avoid this if possible.
  if (permeability /= this%permeability) then
    this%permeability = permeability
    this%pct = this%pct_a * permeability ** this%pct_exp
    ! Update unsaturated extension to reflect new Pct
    error = this%setSwj(this%Swj)
  end if
end subroutine

! **************************************************************************** !

subroutine SFWIPPIgnoreK(this, permeability)
  class(sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: permeability
end subroutine

! **************************************************************************** !
! WIPP KPC Subroutines
! **************************************************************************** !

function SFWIPPKPC1Swj(this,Swj) result (error)
  class (sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Swj
  PetscInt :: error

  ! Ignore input, set Swj to Swr
  this%Swj = this%Swr
  this%Pcj = huge(this%Pcj)
  error = 0
end function

! **************************************************************************** !

pure subroutine SFWIPPKPC1Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  Pc = 0d0
  dPc_dSw = 0d0
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKPC1Sw(this, Pc, Sw)
 class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw

  Sw = this%Swr
end subroutine

! **************************************************************************** !

function SFWIPPKPC2Swj(this,Swj) result (error)
  class (sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Swj
  PetscInt :: error

  ! Ignore input, calculate Swj based on Pcmax
  call this%KRPSw(this%Pcmax, this%Swj)
  this%Pcj = this%Pcmax
  this%dPcj_dSwj = 0d0 ! LHS derivative TODO consider RHS derivative

  error = 0
end function

! **************************************************************************** !

pure subroutine SFWIPPKPC2Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  Pc = this%Pcmax
  dPc_dSw = this%dPcj_dSwj
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKPC2Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw

  Sw = this%Swj
end subroutine

! **************************************************************************** !

function SFWIPPKPC6Swj(this,Swj) result (error)
  class(sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Swj
  PetscInt :: error

  if (Swj > this%Swr) then ! Linearly extrapolate from the valid Swj
    this%Swj = Swj
    call this%KRPPc(Swj, this%Pcj, this%dPcj_dSwj)
    this%Pcmax = this%Pcj - this%dPcj_dSwj * Swj
    error = 0
  else ! Invalid Swj
    error = 1
  end if
end function

! **************************************************************************** !

pure subroutine SFWIPPKPC6Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  if (Sw > 0d0) then                                  ! Linear interpolation
    Pc = this%Pcmax + this%dPcj_dSwj*Sw
  else                                                ! y-intercept
    Pc = this%Pcmax
  end if
  dPc_dSw = this%dPcj_dSwj
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKPC6Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw

  if (Pc < this%Pcmax) then                           ! Linear interpolation
    Sw = (Pc - this%Pcmax) / this%dPcj_dSwj
  else                                                ! y-intercept
    Sw = 0d0
  end if
end subroutine

! **************************************************************************** !

function SFWIPPKRP12Swj(this, Swj) result (error)
  class(sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Swj
  PetscInt :: error
  ! Update Pcmax and dPcj_dSwj for a fixed Swj due to changes in Pct
  call this%KRPPc(this%Swj, this%Pcmax, this%dPcj_dSwj)
  error = 0
end function

! **************************************************************************** !
! WIPP KRP Subroutines
! **************************************************************************** !

pure subroutine SFWIPPVGPc(this, Sw, Pc, dPc_dSw)
! Author: Heeho Park; Modified by Jennifer Frederick; Refactored matpaul
! Date: 11/17/16; Modified 04/26/2017 ; Refactored 12/1/2021
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se, Se_mrtrec, aPc_n

  if (Sw >= this%Sgr_comp) then ! Saturated limit
    Pc = 0d0
    dPc_dSw = -huge(dPc_dSw) ! TODO calculate finite difference limit
  else
    Se = (Sw-this%Swr) * this%dSe_dSw
    Se_mrtrec = Se**this%m_nrec
    aPc_n = Se_mrtrec - 1d0
    Pc = this%Pcm_Pct * this%Pct * aPc_n**this%m_comp
    dPc_dSw = this%k_dSe_dSw * Pc * Se_mrtrec/(aPc_n*Se)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPVGSw(this, Pc, Sw)
! Author: Heeho Park; Modified by Jennifer Frederick; Refactored matpaul
! Date: 11/17/16; Modified 04/26/2017; Refactored 12/1/2021
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: aPc_n, Se_mrtrec, Se

  if (Pc <= 0d0) then
    Sw = 1d0
  else
    aPc_n = (Pc / (this%Pcm_Pct * this%Pct))**this%n
    Se_mrtrec = aPc_n + 1d0
    Se = Se_mrtrec**(-this%m)
    Sw = this%Swr + this%Sw_span * Se
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPBCPc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se

  if (Sw >= this%Sgr_comp) then
    Pc = this%Pct
    dPc_dSw = this%k_dSe_dSw * Pc
  else
    Se      = (Sw - this%Swr) * this%dSe_dSw
    Pc      = this%Pct * Se**this%lambda_nrec
    dPc_dSw = this%k_dSe_dSw * Pc / Se
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPBCSw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se

  if (Pc <= this%Pct) then
    Sw = this%Sgr_comp
  else
    Se = (this%Pct/Pc)**this%lambda
    Sw = this%Swr + this%Sw_span*Se
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP4Pc(this, Sw, Pc, dPc_dSw)
! Brooks-Corey extending past Sgr
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se

  Se = (Sw - this%Swr) * this%dSe_dSw
  Pc = this%Pct * Se**this%lambda_nrec
  dPc_dSw = this%k_dSe_dSw * Pc / Se
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP4Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se

  Se = (this%Pct/Pc)**this%lambda
  Sw = this%Swr + this%Sw_span*Se
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP5Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
!------- Linear model (A)
  if (Sw <= this%Swr) then
    Pc = this%Pcmax
    dPc_dSw = 0d0
  else if (Sw <= this%Sgr_comp) then
    Pc = this%Pct
    dPc_dSw = 0d0
  else
    dPc_dSw = (this%Pct-this%Pcmax)*this%dSe_dSw
    Pc = dPc_dSw*(Sw-this%Swr) + this%Pcmax
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP5Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se

  if (Pc <= 0d0) then
   Sw = 1d0
  else
    Se = (Pc-this%Pcmax)/(this%Pct-this%Pcmax)
    Sw = this%Swr + this%Sw_span*Se
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP9Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
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
  PetscReal :: Sg, Se9

  PetscReal, parameter :: a = 3783.0145d0
  PetscReal, parameter :: b = 2.9d0
  PetscReal, parameter :: b_rec = 1d0/b

  if (Sw <= this%Swr) then
    Pc = 0d0
    dPc_dSw = 0d0
  else
    Sg = 1d0 - Sw
    Se9 = 1d0/Sw - 1d0
    Pc = a*Se9**b_rec
    dPc_dSw = -Pc/(b*Sw*Sg)
  endif
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP9Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se9

  PetscReal, parameter :: a = 3783.0145d0
  PetscReal, parameter :: b = 2.9d0
  PetscReal, parameter :: b_rec = 1d0/b

  if (Pc <= 0d0) then
    Sw = 1d0
  else
    Se9 = (Pc/a)**b
    Sw = 1d0 / (Se9+1d0)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP11Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  Pc = 0d0
  dPc_dSw = 0d0
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP11Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  Sw = 1d0
end subroutine

! **************************************************************************** !
! Relative Permeability Functions
! **************************************************************************** !

function RPFWIPPctor(liquid, KRP, Swr, Sgr, expon) result (new)
  class(rpf_WIPP_type), pointer :: new
  PetscBool, intent(in) :: liquid
  PetscInt , intent(in) :: KRP
  PetscReal, intent(in) :: Swr, Sgr, expon
  PetscInt :: error

  ! Memory allocation
  allocate(new)
  if (.not. associated(new)) return ! Memory allocation failed, abort

  ! Derivatives have been defined to enable smooth unsaturated extensions
  new%analytical_derivative_available = .TRUE.

  ! Data validation
                                        error = 0
  ! If KRP is in branch table, set function pointers, else flag error
  if (liquid) then
    select case(KRP)
    case (KRP_VG_range)
      new%KRPKr => RPFWIPPMVGKrw
    case (KRP_BC_range)
      new%KRPKr => RPFWIPPBBCKrw
    case (KRP_lin_range)
      new%KRPKr => RPFWIPPlinKrw
    case (9)
      new%KRPKr => RPFWIPPKRP9Krw
    case default
                                        error = error + 1
    end select
  else
    select case(KRP)
    case (KRP_VG_range)
      new%KRPKr => RPFWIPPMVGKrg
    case (KRP_BC_range)
      new%KRPKr => RPFWIPPBBCKrg
    case (KRP_lin_range)
      new%KRPKr => RPFWIPPlinKrg
    case (9)
      new%KRPKr => RPFWIPPKRP9Krg
    case default
                                        error = error + 2
    end select
  end if

  ! Check the residual saturations, and their sum, are bound between 0 and 1
  if (Swr < 0d0 .or. Swr >= 1d0)        error = error + 4
  if (Sgr < 0d0 .or. Sgr >= 1d0)        error = error + 8
  if (Swr + Sgr > 1d0)                  error = error + 12

  ! Check the exponent parameter is valid for the KRP option chosen
  select case(KRP)
  case (KRP_BC_range) ! Brooks-Corey types
    if (expon <= 0d0)                   error = error + 16
  case (KRP_VG_range) ! Van Genuchten types
    if (expon <= 0d0 .or. expon >= 1d0) error = error + 16
  case default ! Other KRP functions do not use lambda or m
  end select

  ! Abort if errors caught in data validation
  if (error /= 0) then
    ! Potentally write which parameters were invalid to error stream
    deallocate(new)
    nullify(new)
    return
  end if

  ! Assign residual saturation parameters
  new%Swr = Swr
  if (liquid) then
    select case(KRP)
    case (KRP_Se1_Krw_range)
      new%Sgr_comp = 1d0
    case (KRP_Se2_Krw_range)
      new%Sgr_comp = 1d0 - Sgr
    end select
  else
    select case(KRP)
    case (KRP_Se1_Krg_range)
      new%Sgr_comp = 1d0
    case (KRP_Se2_Krg_range)
      new%Sgr_comp = 1d0 - Sgr
    end select
  end if

  ! Calculate span and derivative
  new%Sw_span = new%Sgr_comp - new%Swr
  if (KRP == 11) then ! Reduce linear span for KRP 11
    new%Sw_span = expon * new%Sw_span
    if (liquid) then
      new%Sgr_comp = new%Swr + new%Sw_span
    else
      new%Swr = new%Sgr_comp - new%Sw_span
    end if
  end if
  new%dSe_dSw = 1d0 / new%Sw_span

  ! Assign exponential BC and VG parameters 
  select case(KRP)
  case (KRP_BC_range) ! Brooks-Corey types
    new%m     = 0d0
    new%m_rec = 0d0
    if (liquid) then 
      new%expon = 3d0 + 2d0/expon
    else
      new%expon = 1d0 + 2d0/expon
    end if
  case (KRP_VG_range) ! Van Genuchten types
  ! VG functions embed the closed form Mualem condition n = 1/(1-m)
    new%m     = expon
    new%m_rec = 1d0 / new%m
    new%expon = 2d0 * new%m
  case default        ! All others
    new%m     = 0d0
    new%m_rec = 0d0
    new%expon = 0d0
  end select

! Cache values at saturated limits
! TODO check this behaves as expected at limits
  call new%KRPKr(new%Swr, new%Kr_Swr, new%dKr_dSwr)
  call new%KRPKr(new%Sgr_comp, new%Kr_Sgr, new%dKr_dSgr)

end function

! **************************************************************************** !

subroutine RPFWIPPRelativePermeability(this, liquid_saturation, &
                  relative_permeability, dkr_sat, option)
  class(rpf_WIPP_type)             :: this
  PetscReal, intent(in)            :: liquid_saturation
  PetscReal, intent(out)           :: relative_permeability, dkr_sat
  type(option_type), intent(inout) :: option

  ! Use cached values at saturated limits
  if (liquid_saturation <= this%Swr) then
    relative_permeability = this%Kr_Swr
    dkr_sat = this%dKr_dSwr
  else if (liquid_saturation >= this%Sgr_comp) then
    relative_permeability = this%Kr_Sgr
    dkr_sat = this%dKr_dSgr
  else
    call this%KRPKr(liquid_saturation, relative_permeability, dkr_sat)
  end if
end subroutine

! **************************************************************************** !

pure subroutine RPFWIPPMVGKrw(this, Sw, Kr, dKr_dSw)
  class(rpf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Kr, dKr_dSw
  PetscReal :: Se, Se_mrt, Se_mrt_comp, Se_mrt_comp_m, f

  Se = (Sw - this%Swr) * this%dSe_dSw
  Se_mrt = Se**this%m_rec
  Se_mrt_comp = 1d0 - Se_mrt
  Se_mrt_comp_m = Se_mrt_comp**this%m
  f = 1d0 - Se_mrt_comp_m

  Kr = sqrt(Se)*f*f
  dKr_dSw = this%dSe_dSw * Kr / Se * &
          (0.5d0 + 2d0*Se_mrt*Se_mrt_comp_m/(f*Se_mrt_comp))
end subroutine

! **************************************************************************** !

pure subroutine RPFWIPPMVGKrg(this, Sw, Kr, dKr_dSw)
  class(rpf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Kr, dKr_dSw
  PetscReal :: Se, Se_comp, Se_mrt, Se_mrt_comp

  Se = (Sw - this%Swr) * this%dSe_dSw
  Se_comp = 1.d0 - Se
  Se_mrt = Se**this%m_rec
  Se_mrt_comp = 1d0 - Se_mrt

  Kr = sqrt(Se_comp)*Se_mrt_comp**this%expon
  dKr_dSw = -this%dSe_dSw * Kr * (0.5d0/Se_comp + 2d0*Se_mrt/(Se*Se_mrt_comp))
end subroutine

! **************************************************************************** !

pure subroutine RPFWIPPBBCKrw(this, Sw, Kr, dKr_dSw)
  class(rpf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Kr, dKr_dSw
  PetscReal :: Se

  Se = (Sw - this%Swr) * this%dSe_dSw

  Kr = Se**this%expon
  dKr_dSw = this%dSe_dSw*this%expon*Kr/Se
end subroutine

! **************************************************************************** !

pure subroutine RPFWIPPBBCKrg(this, Sw, Kr, dKr_dSw)
  class(rpf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Kr, dKr_dSw
  PetscReal :: Se, Se_comp, Se_expon, Se_expon_comp

  Se = (Sw - this%Swr) * this%dSe_dSw
  Se_comp = 1d0 - Se
  Se_expon_comp = 1d0 - Se**this%expon

  Kr = Se_comp*Se_comp*Se_expon_comp
  dKr_dSw = -this%dSe_dSw * Kr * (2d0/Se_comp + this%expon/Se)
end subroutine

! **************************************************************************** !

pure subroutine RPFWIPPlinKrw(this, Sw, Kr, dKr_dSw)
  class(rpf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Kr, dKr_dSw

  Kr = (Sw - this%Swr) * this%dSe_dSw
  dKr_dSw = this%dSe_dSw
end subroutine

! **************************************************************************** !

pure subroutine RPFWIPPlinKrg(this, Sw, Kr, dKr_dSw)
  class(rpf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Kr, dKr_dSw

  Kr = 1d0 - (Sw - this%Swr) * this%dSe_dSw
  dKr_dSw = -this%dSe_dSw
end subroutine

! **************************************************************************** !

pure subroutine RPFWIPPKRP9Krw(this, Sw, Kr, dKr_dSw)
  class(rpf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Kr, dKr_dSw

  PetscReal, parameter :: a = 28.768353d0
  PetscReal, parameter :: b = 1.7241379d0

  PetscReal :: Se9, dKr_dSe, dSe_dsat

  Se9 = 1d0/Sw - 1d0
  
  Kr = 1d0/(1d0 + a*Se9**b)

  dKr_dSe = -Se9**(b-1d0)*a*b/(Se9**b*a + 1d0)**2
  dSe_dsat = 1d0/(Se9*Se9)
  dKr_dSw = dKr_dSe * dSe_dsat
end subroutine

! **************************************************************************** !

pure subroutine RPFWIPPKRP9Krg(this, Sw, Kr, dKr_dSw)
  class(rpf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Kr, dKr_dSw
  PetscReal, parameter :: a = 28.768353d0
  PetscReal, parameter :: b = 1.7241379d0
  PetscReal :: Se9, dKr_dSe, dSe_dsat

  Se9 = 1d0/Sw - 1d0

  Kr = 1d0 -1d0/(1d0 + a*Se9**b)

  dKr_dSe = -Se9**(b-1d0)*a*b/(Se9**b*a + 1d0)**2
  dSe_dsat = 1d0/(Se9*Se9)
  dKr_dSw = -dKr_dSe * dSe_dsat
end subroutine

! **************************************************************************** !

end module
