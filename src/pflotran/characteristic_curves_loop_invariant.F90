module Characteristic_Curves_loop_invariant_module
#include "petsc/finclude/petscsys.h"

use petscsys ! Necessary for PETSC_TRUE / PETSC_FALSE
use Characteristic_Curves_Base_module   ! Needed to define base type
use Characteristic_Curves_Common_module ! Needed to inherit VG type
use Option_module ! Needed for Verify any unused arguments in Pc and Sl

implicit none

private

! **************************************************************************** !
!
! References:
!
! Chen, J, JW Hopmans, and ME Grismer (1999) "Parameter estimation of two-fluid
! capillary pressure-saturation and pemeability functions", Adv Water Resour
! 22(5):479-493. doi: 10.1016/s0309-1708(98)00025-6
!
! Sun, Y, et al. (2010) "Modeling Thermal-Hydrologic Processes for a Heated
! Fractured Rock System: Impact of a Capillary-Pressure Maximum",
! Transp Porous Med 83:501-523. doi: 10.1007/s11242-009-9459-1
!
! Van Genuchten, MT (1980) "A Closed-form Equation for Predicting the
! Hydraulic ! Conductivity of Unsaturated Soils", Soil Sci Soc Am J
! 44(5):892-898. doi: 10.2136/sssaj1980.03615995004400050002x
!
! Note: As the VG capillary pressure function has an infinite derivative
! approaching saturation, derivatives above Sl_max are approximated using a
! backwards finite difference method using machine epsilon.
!
! **************************************************************************** !

PetscReal, private, parameter :: Sl_max = 1d0 - epsilon(Sl_max) ! Saturated limita

! **************************************************************************** !
!
! VG Saturation Function Type Declarations
!
! sat_func_base_type        External base type, quasi-abstract
! |
! |-->sat_func_VG_type      External common type, VG constant extension
!     |
!     |-->sf_VG_type        VG type with loop-invariant parameters, no extension
!         |
!         |-->sf_VG_extn_type*     Abstract VG type for extensions
!             |
!             |-->sf_VG_cons_type  Constant    extension type
!             |-->sf_VG_expn_type  Exponential extension type
!             |-->sf_VG_line_type  Linear      extension type
!             |-->sf_VG_quad_type  Quadratic   extension type
!
! Objects of sf_VG_type are created and initialized with the constructor:
! sf_VG_type       => SF_VG_ctor(unsat_ext,alpha,m,Sr,rpf,Pcmax,Sj)
!
! Warning, Pcmax and Sr are unprotected as Fortran will not extend private
! scope to child classed in separate modules.
!
! **************************************************************************** !
! VG Saturation Function constructors
! **************************************************************************** !

  public  :: SFVGCtor         ! Generic constructor
  private :: SFVGNEVGCtor, & ! No capped capillary pressure
             SFVGFCPCCtor, & ! Flat, specified maximum
             SFVGFNOCCtor, & ! Flat, specified junction
             SFVGECPCCtor, & ! Exponential, specified maximum
             SFVGENOCCtor, & ! Exponential, specified junction
             SFVGLCPCCtor, & ! Linear, specified maximum
             SFVGLNOCCtor, & ! Linear, specified junction
             SFVGQUADCtor    ! Quadratic, specified maximum and junction

! **************************************************************************** !
! Van Genuchten Saturation Function type defintions
! **************************************************************************** !
  type, public, extends(sat_func_VG_type) :: sf_VG_type
!   public                           ! Unprotected parameters from parent types
!     PetscReal :: Sr                ! Base   - Residual saturation
!     PetscReal :: Pcmax             ! Base   - Maximum capillary pressure
!     PetscReal :: alpha             ! Common - van Genuchten coefficient * Pa
!     PetscReal :: m                 ! Common - van Genuchten exponent
    private                          ! Loop-invariant parameters
      PetscReal :: a_rec             ! Alpha reciprocal
      PetscReal :: m_nrec, m_a2      ! Negative reciprocal M, M add 2
      PetscReal :: n, n_rec, n_m1    ! N, reciprocal N, and N minus 1
      PetscReal :: mn_a1             ! MN product add 1
      PetscReal :: Sl_span, dSe_dSl  ! Effective saturation span and reciprocal
      PetscReal :: dSe_mndSl         ! Coefficent in VG derivative
      PetscReal :: dPc_dSl_max       ! Finite difference derivative at saturation
      PetscReal :: dSl_dPcmin        ! Finite difference derivative at saturation
      PetscReal :: d2Sl_dPc2min      ! Finite difference derivative at saturation
      PetscInt  :: rpf               ! Mualem/Burdine model flag TODO add custom
  contains
! Definition of base type methods
    procedure, public  :: Init                  => SFVGInit
    procedure, public  :: CapillaryPressure     => SFVGCapillaryPressure
    procedure, public  :: Saturation            => SFVGSaturation
    procedure, public  :: D2SatDP2              => SFVGD2SatDP2
!   procedure, public  :: Verify
! Common VG methods
    procedure, private :: Configure             => SFVGConfigure
    procedure, private :: PcInline              => SFVGPcInline
    procedure, private :: SlInline              => SFVGSlInline
    procedure, private :: D2SlDPc2Inline        => SFVGD2SlDPc2Inline
    procedure, private :: SlInflection          => SFVGSlInflection
! No-extension mutator methods
    procedure, private :: Set_alpha             => SFVGSetAlpha
    procedure, private :: Set_m                 => SFVGSetM
! Internal no-extension pure methods
    procedure, private :: Pc                    => SFVGPc
    procedure, private :: Sl                    => SFVGSl
    procedure, private :: D2SlDPc2              => SFVGD2SlDPc2
  end type

! VG Unsaturated Extension Abstract type
    type, abstract, extends(sf_VG_type) :: sf_VG_extn_type
      private
        PetscReal :: Pj, Sj           ! Piecewise junction point
        PetscBool :: Pcmax_designated ! Flag for designated Pcmax
    contains
      ! Common mutator methods
      procedure, public :: SetAlpha             => SFVGextnSetAlpha
      procedure, public :: SetM                 => SFVGextnSetM
      ! Mutator prototype
      procedure(set_Pcmax_type), deferred, public    :: SetPcmax
      procedure(set_Sj_type   ), deferred, public    :: SetSj
    end type

! VG Constant Unsaturated Extension
      type, public, extends(sf_VG_extn_type) :: sf_VG_cons_type
        private
          PetscReal :: dSj_dPj   ! Non-zero derivative at cusp
          PetscReal :: d2Sj_dPj2 ! Non-zero 2nd derivative at cusp
      contains
        ! Public mutator methods
        procedure, public  :: SetPcmax          => SFVGconsSetPcmax
        procedure, public  :: SetSj             => SFVGconsSetSj
        ! Internal methods for constant extensions
        procedure, private :: Pc                => SFVGconsPc
        procedure, private :: Sl                => SFVGconsSl
        procedure, private :: D2SlDPc2          => SFVGconsD2SlDPc2
      end type

! VG Exponential Unsaturated Extension
      type, public, extends(sf_VG_extn_type) :: sf_VG_expn_type
        private
          PetscReal :: beta, beta_rec   ! Exponential coefficient and reciprocal
          PetscReal :: dPcmax_dSl       ! Unsaturated limits
          PetscReal :: dSl_dPcmax       ! Derivative at unsaturated limit
          PetscReal :: d2Sl_dPc2max     ! 2nd derivative at unsaturated limit
      contains
        ! Public mutator methods
        procedure, public  :: SetPcmax          => SFVGexpnSetPcmax
        procedure, public  :: SetSj             => SFVGexpnSetSj
        ! Internal methods for exponential extensions
        procedure, private :: Pc                => SFVGexpnPc
        procedure, private :: Sl                => SFVGexpnSl
        procedure, private :: D2SlDPc2          => SFVGexpnD2SlDPc2
      end type

! VG Linear Unsaturated Extension
      type, public, extends(sf_VG_extn_type) :: sf_VG_line_type
        private
          PetscReal :: dPj_dSj, dSj_dPj ! Linear coefficients
      contains
        ! Public mutator methods
        procedure, public  :: SetPcmax          => SFVGlineSetPcmax
        procedure, public  :: SetSj             => SFVGlineSetSj
        ! Internal methods for linear extensions
        procedure, private :: Pc                => SFVGlinePc
        procedure, private :: Sl                => SFVGlineSl
        procedure, private :: D2SlDPc2          => SFVGlineD2SlDPc2
      end type

! VG Quadratic Unsaturated Extension
      type, public, extends(sf_VG_extn_type) :: sf_VG_quad_type
        private
          PetscReal :: A, B            ! Quadratic coefficients
          PetscReal :: dSl_dPcmax
          PetscReal :: d2Sl_dPc2max
          PetscBool :: qa_branch       ! Select inverse solution
      contains
        ! Public mutator methods
        procedure, public  :: SetPcmax          => SFVGquadSetPcmax
        procedure, public  :: SetSj             => SFVGquadSetSj
        procedure, public  :: SetQuad           => SFVGquadSetQuad
        ! Internal methods for quadratic extensions
        procedure, private :: Pc                => SFVGquadPc
        procedure, private :: Sl                => SFVGquadSl
        procedure, private :: D2SlDPc2          => SFVGquadD2SlDPc2
        procedure, private :: SlRoot            => SFVGquadSlRoot
      end type

! **************************************************************************** !
! Saturation Function Prototype Block
! **************************************************************************** !

abstract interface
  function set_Pcmax_type(this, Pcmax) result (error)
    import :: sf_VG_extn_type
    class(sf_VG_extn_type), intent(inout) :: this
    PetscReal, intent(in) :: Pcmax
    PetscInt :: error
  end function

! **************************************************************************** !

  function set_Sj_type(this, Sj) result (error)
    import :: sf_VG_extn_type
    class(sf_VG_extn_type), intent(inout) :: this
    PetscReal, intent(in) :: Sj
    PetscInt :: error
  end function
end interface

! **************************************************************************** !
! VG Relative Permeability Function Type Declarations
! **************************************************************************** !
!
! rel_perm_func_base_type   External definition in characteristic_curves_base
! |
! |-->rpf_VG_type*                   Abstract base type
!    |
!    |-->rpf_VG_liq_type*            Abstract liquid base type
!    |   |
!    |   |-->rpf_MVG_liq_type        Mualem  - van Genuchten - Liquid
!    |   |-->rpf_BVG_liq_type        Burdine - van Genuchten - Liquid
!    |
!    |-->rpf_VG_gas_type*            Abstract gas base type
!        |
!        |-->rpf_MVG_gas_type        Mualem  - van Genuchten - Gas
!        |-->rpf_BVG_gas_type        Burdine - van Genuchten - Gas
!
! Caution, Sr is unprotected as Fortran will not extend private
! scope to child classed in separate modules. The rel_perm_base_type would
! need to be modified.
!
! Due to Fortran scoping rules, the loop-invariant types do not extend the
! loop-variant types. Should loop-invariant replace the original types, the
! classis() statements throughout the code must be updated.
!
! **************************************************************************** !
! RPF VG Contructors
! **************************************************************************** !

  public :: RPFMVGliqCtor, &
            RPFMVGgasCtor, &
            RPFBVGliqCtor, &
            RPFBVGgasCtor

! **************************************************************************** !
! RPF type definitions
! **************************************************************************** !

type, abstract, extends(rel_perm_func_base_type) :: rpf_VG_type
! public
!   PetscReal :: Sr         ! Base   - Residual saturation
  private
    PetscReal :: m, m_rec   ! Van Genuchten m exponent and reciprocal
    PetscReal :: dSe_dSl    ! Effective saturation span reciprocal
contains
  ! Overloaded base methods
  procedure, public :: Init                 => RPFVGInit
  procedure, public :: RelativePermeability => RPFVGRelativePermeability
  ! Public accessor/mutator methods
  procedure, private                         :: GetM  => RPFVGGetM
  procedure(set_m_type)  , deferred, public  :: SetM
  ! Private deferred methods called by RelativePermeability
  procedure(calc_Kr_type), deferred, private :: Kr
  procedure(calc_Kr_type), deferred, private :: KrInline
end type

! **************************************************************************** !

  type, abstract, extends(rpf_VG_type) :: rpf_VG_liq_type
    private
      PetscReal :: dKr_dSl_max ! Finite difference derivative at saturation
      PetscReal :: Sl_min      ! Unsaturated limit to avoid NaN
  contains
    procedure, public  :: SetM                  => RPFVGliqSetM
    procedure, private :: Configure             => RPFVGliqConfigure
    procedure, private :: Kr                    => RPFVGliqKr
  end type

! **************************************************************************** !

    type, public, extends(rpf_VG_liq_type) :: rpf_MVG_liq_type
    contains
      procedure, private :: KrInline            => RPFMVGliqKrInline
    end type

! **************************************************************************** !

    type, public, extends(rpf_VG_liq_type) :: rpf_BVG_liq_type
    contains
      procedure, private :: KrInline            => RPFBVGliqKrInline
    end type

! **************************************************************************** !

  type, abstract, extends(rpf_VG_type) :: rpf_VG_gas_type
  private
    PetscReal :: mx2        ! m times 2
    PetscReal :: Sgr        ! Gas residual saturation
    PetscReal :: Sl_max     ! Saturated limit reduced by Sgr
  contains
    procedure, public  :: SetM                   => RPFVGgasSetM
    procedure, private :: Configure              => RPFVGgasConfigure
    procedure, private :: Kr                     => RPFVGgasKr
  end type

! **************************************************************************** !

    type, public, extends(rpf_VG_gas_type) :: rpf_MVG_gas_type
    contains
      procedure, private :: KrInline            => RPFMVGgasKrinline
    end type

! **************************************************************************** !

    type, public, extends(rpf_VG_gas_type) :: rpf_BVG_gas_type
    contains
      procedure, private :: KrInline            => RPFBVGgasKrinline
    end type

! **************************************************************************** !
! Relative Permeability Function Prototype Block
! **************************************************************************** !

abstract interface
  function set_m_type(this, m) result (error)
    import :: rpf_VG_type
    class(rpf_VG_type), intent(inout) :: this
    PetscReal, intent(in) :: m
    PetscInt :: error
  end function

! **************************************************************************** !

  pure subroutine calc_Kr_type(this, Sl, Kr, dKr_dSl)
    import :: rpf_VG_type
    class(rpf_VG_type), intent(in) :: this
    PetscReal, intent(in)  :: Sl
    PetscReal, intent(out) :: Kr, dKr_dSl
  end subroutine
end interface

contains

! **************************************************************************** !
! VG Saturation Function Methods
! **************************************************************************** !

function SFVGCtor(unsat_ext, alpha, m, Sr, Sgt, vg_rpf_opt, Pcmax, Sj) &
         result (new)

  implicit none

  class(sf_VG_type), pointer :: new
  character(*), intent(in) :: unsat_ext
  PetscReal, intent(in) :: alpha, m, Sr, Sgt, Pcmax, Sj
  PetscInt, intent(in) :: vg_rpf_opt

! This function returns a the van Genuchten saturation function object using
! the correct extension constructor method

  select case (unsat_ext)
  case ('NONE') ! No extension
    new => SFVGNEVGCtor(alpha,m,Sr,Sgt,vg_rpf_opt)
  case ('FCPC') ! Flat specified cap
    new => SFVGFCPCCtor(alpha,m,Sr,Sgt,vg_rpf_opt,Pcmax)
  case ('FNOC') ! Flat specificed junction
    new => SFVGFNOCCtor(alpha,m,Sr,Sgt,vg_rpf_opt,Sj)
  case ('ECPC') ! Exponential specified cap
    new => SFVGECPCCtor(alpha,m,Sr,Sgt,vg_rpf_opt,Pcmax)
  case ('ENOC') ! Exponential specified junction
    new => SFVGENOCCtor(alpha,m,Sr,Sgt,vg_rpf_opt,Sj)
  case ('LCPC') ! Linear specified cap
    new => SFVGLCPCCtor(alpha,m,Sr,Sgt,vg_rpf_opt,Pcmax)
  case ('LNOC') ! Linear specified junction
    new => SFVGLNOCCtor(alpha,m,Sr,Sgt,vg_rpf_opt,Sj)
  case ('QUAD') ! Quadratic specified cap and junction
    new => SFVGquadCtor(alpha,m,Sr,Sgt,vg_rpf_opt,Pcmax,Sj)
  case default
    nullify(new)
  end select
end function SFVGCtor

! **************************************************************************** !

subroutine SFVGInit(this)

  implicit none

  class(sf_VG_type) :: this
  ! This method is intentionally left blank.
end subroutine SFVGInit

! **************************************************************************** !

subroutine SFVGCapillaryPressure(this, liquid_saturation, &
                                   capillary_pressure, dPc_dSatl, option)

  implicit none

  class(sf_VG_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option

  call this%Pc(liquid_saturation, capillary_pressure, dPc_dSatl)
end subroutine SFVGCapillaryPressure

! **************************************************************************** !

subroutine SFVGSaturation(this, capillary_pressure, &
                            liquid_saturation, dsat_dpres, option)

  implicit none

  class(sf_VG_type) :: this
  PetscReal, intent(in)  :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  call this%Sl(capillary_pressure, liquid_saturation, dsat_dpres)
  dsat_dpres = -dsat_dpres ! Replicating existing reversed signed behaivor
end subroutine SFVGSaturation

! **************************************************************************** !

subroutine SFVGD2SatDP2(this,Pc, d2s_dp2, option)

  implicit none

  class(sf_VG_type) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option

  call this%D2SlDPc2(Pc, d2s_dp2)
end subroutine SFVGD2SatDP2

! **************************************************************************** !

function SFVGConfigure(this, alpha, m, Sr, Sgt, rpf) result (error)

  implicit none

  ! Configure loop-invariant parameters common to all loop-invariant types
  class(sf_VG_type) :: this
  PetscReal, intent(in) :: alpha,m,Sr,Sgt
  PetscInt, intent(in) :: rpf
  PetscInt :: error
  PetscReal :: Pc1, Pc2, dPc_dSl

  ! Can eliminate the pointers if legacy smoothing implementation is removed
  nullify(this%sat_poly)
  nullify(this%pres_poly)
  this%analytical_derivative_available = PETSC_TRUE

  error = 0                  ! Using bitwise errorcodes
  if (alpha <= 0d0)            error = error + 1
  if (m <= 0d0 .OR. m >= 1d0)  error = error + 2
  if (Sr < 0d0 .OR. Sr >= 1d0) error = error + 4
  if (rpf < 1 .OR. rpf > 2)    error = error + 8

  if (error == 0) then
    ! Calculate loop-invariant parameters if data validation is successful
    this%alpha = alpha
    this%a_rec = 1d0 / alpha
    this%m = m
    this%m_a2 = m + 2d0
    this%m_nrec = -1d0 / m
    this%rpf = rpf
    this%Sgt_max = Sgt

    ! While the Mualem is more common, the Burdine assumption is also supported
    select case (rpf)
    case (1) ! Mualem
      this%n_rec = 1d0 - m
      this%n = 1d0 / this%n_rec
    case (2) ! Burdine
      this%n_rec = (1d0 - m) / 2d0
      this%n = 2d0 / (1d0 - m)
    end select
    this%n_m1 = this%n - 1d0
    this%mn_a1 = m*this%n + 1d0
    this%Sr = Sr
    this%Sl_span = 1d0 - Sr
    this%dSe_dSl = 1d0 / this%Sl_span
    this%dSe_mndSl = this%n_rec / (this%m*this%Sl_span)

    ! Estimate saturated limits using backwards finite difference on Sl
    ! Pc(1.0) = 0; Pc1 = Pc(1.0-eps); Pc2 = Pc(1.0-2.0*eps);
    call this%PcInline(Sl_max, Pc1, dPc_dSl)
    call this%PcInline(Sl_max-epsilon(Sl_max), Pc2, dPc_dSl)
    this%dPc_dSl_max = -Pc1 / epsilon(Sl_max)
    this%dSl_dPcmin = -epsilon(Sl_max) / Pc1
    this%d2Sl_dPc2min = epsilon(Sl_max)*(2d0-Pc2/Pc1)/Pc1**2
  end if
end function SFVGConfigure

! **************************************************************************** !

pure subroutine SFVGPcInline(this, Sl, Pc, dPc_dSl)

  implicit none

  class(sf_VG_type), intent(in) :: this
  PetscReal, intent(in) :: Sl
  PetscReal, intent(out) :: Pc, dPc_dSl
  PetscReal :: Se, Se_mrtrec, aPc_n

  Se = (Sl - this%Sr) * this%dSe_dSl
  Se_mrtrec = Se**this%m_nrec
  aPc_n = Se_mrtrec - 1d0

  Pc = this%a_rec * aPc_n**this%n_rec
  dPc_dSl = (this%dSe_mndSl*Pc) / (Se/Se_mrtrec-Se)
end subroutine SFVGPcInline

! **************************************************************************** !

pure subroutine SFVGSlInline(this, Pc, Sl, dSl_dPc)

  implicit none

  class(sf_VG_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: Sl, dSl_dPc
  PetscReal :: aPc_n, Se_mrtrec, Se

  aPc_n = (this%alpha*Pc)**this%n
  Se_mrtrec = aPc_n + 1d0
  Se = Se_mrtrec**(-this%m)

  Sl = this%Sr + Se * this%Sl_span
  dSl_dPc = (Se/Se_mrtrec-Se) / (this%dSe_mndSl*Pc)
end subroutine SFVGSlInline

! **************************************************************************** !

pure subroutine SFVGD2SlDPc2Inline(this, Pc, d2Sl_dPc2)

  implicit none

  class(sf_VG_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2
  PetscReal :: aPc_n, Se_mrtrec, d2Se_mndPc2

  aPc_n = (this%alpha*Pc)**this%n
  Se_mrtrec = aPc_n + 1d0
  d2Se_mndPc2 =  aPc_n*(this%mn_a1*aPc_n-this%n_m1)/(Se_mrtrec**this%m_a2*Pc**2)
  ! Note, the exponentiation of Se_mrtrec could be reduced if Se is provived
  ! Se_mrtrec**(m+2) == Se_mrtrec**2/Se
  ! As it is not, the terms are combined

  d2Sl_dPc2 = d2Se_mndPc2 / this%dSe_mndSl
end subroutine SFVGD2SlDPc2Inline

! **************************************************************************** !

pure function SFVGSlInflection(this) result (Sl)

  implicit none

  class(sf_VG_type), intent(in) :: this
  PetscReal :: Se_mrtrec, Se, Sl

  Se_mrtrec = (1d0+this%m)/(this%n_rec+this%m)
  Se = Se_mrtrec**(-this%m)
  Sl = this%Sr + Se*this%Sl_span
end function SFVGSlInflection

! **************************************************************************** !
! Van Genuchten
! **************************************************************************** !

function SFVGNEVGCtor(alpha,m,Sr,Sgt,rpf) result (new)

  implicit none

  class(sf_VG_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sgt
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,Sgt,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  end if
  new%Pcmax = huge(new%Pcmax)
end function SFVGNEVGCtor

! **************************************************************************** !

function SFVGSetAlpha(this,alpha) result (error)

  implicit none

  class(sf_VG_type), intent(inout) :: this
  PetscReal, intent(in) :: alpha
  PetscInt :: error
  error = this%Configure(alpha,this%m,this%Sr,this%Sgt_max,this%rpf)
end function SFVGSetAlpha

! **************************************************************************** !

function SFVGSetM(this,m) result (error)

  implicit none

  class(sf_VG_type), intent(inout) :: this
  PetscReal, intent(in) :: m
  PetscInt :: error
  error = this%Configure(this%alpha,m,this%Sr,this%Sgt_max,this%rpf)
end function SFVGSetM

! **************************************************************************** !

pure subroutine SFVGPc(this, Sl, Pc, dPc_dSl)

  implicit none

  class(sf_VG_type), intent(in) :: this
  PetscReal, intent(in)   :: Sl
  PetscReal, intent(out)  :: Pc, dPc_dSl

  if (Sl <= this%Sr) then               ! Unsaturated limit
    Pc = huge(Pc)
    dPc_dSl = -huge(Pc)
  else if (Sl > Sl_max) then            ! Saturated limit
    Pc= 0d0
    dPc_dSl = this%dPc_dSl_max
  else                                  ! Ordinary VG domain
    call this%PcInline(Sl,Pc,dPc_dSl)
  end if
end subroutine SFVGPc

! **************************************************************************** !

pure subroutine SFVGSl(this, Pc, Sl, dSl_dPc)

  implicit none

  class(sf_VG_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sl, dSl_dPc

  if (Pc <= 0d0) then                   ! Saturated limit
    Sl = 1d0
    dSl_dPc = this%dSl_dPcmin
  else                                  ! Ordinary VG domain
    call this%SlInline(Pc,Sl,dSl_dPc)
  end if
end subroutine SFVGSl

! **************************************************************************** !

pure subroutine SFVGD2SlDPc2(this, Pc, d2Sl_dPc2)

  implicit none

  class(sf_VG_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2

  if (Pc <= 0d0) then                   ! Saturated limit
    d2Sl_dPc2 = this%d2Sl_dPc2min
  else                                  ! Ordinary VG domain
    call this%D2SlDPc2Inline(Pc,d2Sl_dPc2)
  end if
end subroutine SFVGD2SlDPc2

! **************************************************************************** !
! Abstract Van Genuchten Extension
! **************************************************************************** !

function SFVGextnSetAlpha(this,alpha) result (error)

  implicit none

  class(sf_VG_extn_type), intent(inout) :: this
  PetscReal, intent(in) :: alpha
  PetscInt :: error
  PetscReal :: alpha_old

  alpha_old = this%alpha
  error = this%Configure(alpha,this%m,this%Sr,this%Sgt_max,this%rpf)
  if (error == 0) then                  ! Update unsaturated extension
    if (this%Pcmax_designated) then
      error = this%SetPcmax(this%Pcmax)
    else
      error = this%SetSj(this%Sj)
    end if
    if (error /= 0) then                ! Restore previous state upon error
      error = error + this%Configure(alpha_old,this%m,this%Sr,this%Sgt_max, &
                                     this%rpf)
    end if
  end if
end function SFVGextnSetAlpha

! **************************************************************************** !

function SFVGextnSetM(this,m) result (error)

  implicit none

  class(sf_VG_extn_type), intent(inout) :: this
  PetscReal, intent(in) :: m
  PetscInt :: error
  PetscReal :: m_old

  m_old = this%m
  error = this%Configure(this%alpha,m,this%Sr,this%Sgt_max,this%rpf)
  if (error == 0) then                   ! Update unsaturated extension
    if (this%Pcmax_designated) then
      error = this%SetPcmax(this%Pcmax)
    else
      error = this%SetSj(this%Sj)
    end if
    if (error /= 0) then                ! Restore previous state upon error
      error = error + this%Configure(this%alpha,m_old,this%Sr,this%Sgt_max, &
                                     this%rpf)
    end if
  end if
end function SFVGextnSetM

! **************************************************************************** !
! Constant Van Genuchten Extension
! **************************************************************************** !

function SFVGFCPCCtor(alpha,m,Sr,Sgt,rpf,Pcmax) result (new)

  implicit none

  class(sf_VG_cons_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sgt, Pcmax
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,Sgt,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%SetPcmax(Pcmax) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function SFVGFCPCCtor

! **************************************************************************** !

function SFVGconsSetPcmax(this,Pcmax) result (error)

  implicit none

  class(sf_VG_cons_type), intent(inout) :: this
  PetscReal, intent(in) :: Pcmax
  PetscInt :: error

  if (Pcmax > 0d0) then
    error = 0
    this%Pcmax_designated = PETSC_TRUE
    this%Pcmax = Pcmax
    call this%SlInline(Pcmax,this%Sj,this%dSj_dPj)
    call this%D2SlDPc2Inline(Pcmax,this%d2Sj_dPj2)
  else
    error = 1
  end if
end function SFVGconsSetPcmax

! **************************************************************************** !

function SFVGFNOCCtor(alpha,m,Sr,Sgt,rpf,Sj) result (new)

  implicit none

  class(sf_VG_cons_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sgt, Sj
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,Sgt,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%SetSj(Sj) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function SFVGFNOCCtor

! **************************************************************************** !

function SFVGconsSetSj(this,Sj) result (error)

  implicit none

  class(sf_VG_cons_type), intent(inout) :: this
  PetscReal, intent(in) :: Sj
  PetscInt :: error
  PetscReal :: dPj_dSj

  if (Sj > this%Sr) then
    error = 0
    this%Pcmax_designated = PETSC_FALSE
    this%Sj = Sj
    call this%PcInline(Sj,this%Pcmax,dPj_dSj)
    this%dSj_dPj = 1d0 / dPj_dSj
    call this%D2SlDPc2Inline(this%Pcmax,this%d2Sj_dPj2)
  else
    error = 1
  end if
end function SFVGconsSetSj

! **************************************************************************** !

pure subroutine SFVGconsPc(this, Sl, Pc, dPc_dSl)

  implicit none

  class(sf_VG_cons_type), intent(in) :: this
  PetscReal, intent(in)  :: Sl
  PetscReal, intent(out) :: Pc, dPc_dSl

  if (Sl < this%Sj) then                ! Unsaturated limit
    Pc = this%Pcmax
    dPc_dSl = 0d0
  else if (Sl > Sl_max) then            ! Saturated limit
    Pc = 0d0
    dPc_dSl = this%dPc_dSl_max
  else                                  ! Ordinary VG domain
    call this%PcInline(Sl,Pc,dPc_dSl)
  end if
end subroutine SFVGconsPc

! **************************************************************************** !

pure subroutine SFVGconsSl(this, Pc, Sl, dSl_dPc)

  implicit none

  class(sf_VG_cons_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sl, dSl_dPc

  if (Pc <= 0d0) then                   ! Saturated limit
    Sl = 1d0
    dSl_dPc = this%dSl_dPcmin
  else if (Pc >= this%Pcmax) then       ! Unsaturated limit
    Sl = this%Sj
    dSl_dPc = this%dSj_dPj
  else                                  ! Ordinary VG domain
    call this%SlInline(Pc, Sl, dSl_dPc)
  end if
end subroutine SFVGconsSl

! **************************************************************************** !

pure subroutine SFVGconsD2SlDPc2(this,Pc,d2Sl_dPc2)

  implicit none

  class(sf_VG_cons_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2

  if (Pc <= 0d0) then                   ! Saturated limit
    d2Sl_dPc2 = this%d2Sl_dPc2min
  else if (Pc >= this%Pcmax) then       ! Unsaturated limit
    d2Sl_dPc2 = this%d2Sj_dPj2
  else                                  ! Ordinary VG domain
    call this%D2SlDPc2Inline(Pc, d2Sl_dPc2)
  end if
end subroutine SFVGconsD2SlDPc2

! **************************************************************************** !
! Exponential Van Genuchten Extension
! **************************************************************************** !

function SFVGECPCctor(alpha,m,Sr,Sgt,rpf,Pcmax) result (new)

  implicit none

  class(sf_VG_expn_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sgt, Pcmax
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,Sgt,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%SetPcmax(Pcmax) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function SFVGECPCctor

! **************************************************************************** !

function SFVGexpnSetPcmax(this,Pcmax) result (error)

  implicit none

  class(sf_VG_expn_type), intent(inout) :: this
  PetscReal, intent(in) :: Pcmax
  PetscInt :: error
  PetscReal :: Sa, Sb, Pe, dPj_dSj

  ! Solve the nonlinear system to satisfy C1 continuity using bisection method

  ! Set lower saturation limit (Sa) to be intercept (P(Sa) = Pcmax)
  call this%SlInline(Pcmax,Sa,dPj_dSj) ! dPj_dSj is inverted, but not used

  ! Set upper saturation limit (Sb) to be the the inflection point in log space
  Sb = (1d0+1d0/this%m)**(-this%m)
  Sb = this%Sr + Sb*this%Sl_span

  ! Confirm Pcmax is above minimum extrapolating from inflection point
  call this%PcInline(Sb,Pe,dPj_dSj)
  if (Pcmax >= Pe*exp(-dPj_dSj*Sb/Pe)) then
    ! Set Pcmax and begin iteration loop
    error = 0
    this%Pcmax = Pcmax
    this%Pcmax_designated = PETSC_TRUE

    do while (Sb/Sa > 1d0 + epsilon(this%Sj)) ! Tolerance interval epsilon
      this%Sj = (Sa+Sb)/2d0                ! Bisect bracket
      call this%PcInline(this%Sj,this%Pj,dPj_dSj)
      this%beta = dPj_dSj / this%Pj

      ! Residual error = Pcmax*exp(beta*Sj) - Pf(Sj)
      Pe = Pcmax*exp(this%beta*this%Sj) - this%Pj
      if (Pe < 0d0) then ! Error on side a is always negative on LHS for VG
        Sa = this%Sj
      else
        Sb = this%Sj
      end if
    end do
    this%dPcmax_dSl = this%beta*Pcmax
    this%beta_rec = 1d0 / this%beta
    this%dSl_dPcmax = this%beta_rec/Pcmax
    this%d2Sl_dPc2max = -this%beta_rec/Pcmax**2
  else
    error = 1
  end if
end function SFVGexpnSetPcmax

! **************************************************************************** !

function SFVGENOCCtor(alpha,m,Sr,Sgt,rpf,Sj) result (new)

  implicit none

  class(sf_VG_expn_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sgt, Sj
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,Sgt,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%SetSj(Sj) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function SFVGENOCCtor

! **************************************************************************** !

function SFVGexpnSetSj(this,Sj) result (error)

  implicit none

  class(sf_VG_expn_type), intent(inout) :: this
  PetscReal, intent(in) :: Sj
  PetscInt :: error
  PetscReal :: dPj_dSj

  if (Sj > this%Sr) then
    error = 0
    this%Sj = Sj
    this%Pcmax_designated = PETSC_FALSE

    ! Find the slope at the piecewise junction point and extrapolate
    call this%PcInline(Sj,this%Pj,dPj_dSj)
    this%beta = dPj_dSj / this%Pj
    this%beta_rec = this%Pj / dPj_dSj
    this%Pcmax = this%Pj*exp(-this%beta*Sj)
    this%dPcmax_dSl = this%beta*this%Pcmax
    this%dSl_dPcmax = this%beta_rec/this%Pcmax
    this%d2Sl_dPc2max = -this%beta_rec/this%Pcmax**2
  else
    error = 1
  end if
end function SFVGexpnSetSj

! **************************************************************************** !

pure subroutine SFVGexpnPc(this, Sl, Pc, dPc_dSl)

  implicit none

  class(sf_VG_expn_type), intent(in) :: this
  PetscReal, intent(in)   :: Sl
  PetscReal, intent(out)  :: Pc, dPc_dSl

  if (Sl < this%Sj) then                ! Unsaturated domain
    if (Sl <= 0d0) then                 ! Unsaturated limit
      Pc = this%Pcmax
      dPc_dSl = this%dPcmax_dSl
    else                                ! Exponential extension
      Pc = this%Pcmax*exp(this%beta*Sl)
      dPc_dSl = this%beta*Pc
    end if
  else if (Sl >= 1d0) then              ! Saturated limit
    Pc = 0d0
    dPc_dSl = this%dPc_dSl_max
  else                                  ! Ordinary VG domain
    call this%PcInline(Sl, Pc, dPc_dSl)
  end if
end subroutine SFVGexpnPc

! **************************************************************************** !

pure subroutine SFVGexpnSl(this, Pc, Sl, dSl_dPc)

  implicit none

  class(sf_VG_expn_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sl, dSl_dPc

  if (Pc <= 0d0) then                   ! Saturated limit
    Sl = 1d0
    dSl_dPc = this%dSl_dPcmin
  else if (Pc >= this%Pj) then          ! Unsaturated domain
    if (Pc >= this%Pcmax) then          ! Unsaturated limit
      Sl = 0d0
      dSl_dPc = this%dSl_dPcmax
    else                                ! Exponential extension
      Sl = this%beta_rec*log(Pc/this%Pcmax)
      dSl_dPc = this%beta_rec/Pc
    end if
  else                                  ! Ordinary VG domain
    call this%SlInline(Pc, Sl, dSl_dPc)
  end if
end subroutine SFVGexpnSl

! **************************************************************************** !

pure subroutine SFVGexpnD2SlDPc2(this, Pc, d2Sl_dPc2)

  implicit none

  class(sf_VG_expn_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2

  if (pc <= 0d0) then                   ! Saturated limit
    d2Sl_dPc2 = this%d2Sl_dPc2min
  else if (Pc >= this%Pj) then          ! Unsaturated domain
    if (Pc >= this%Pcmax) then          ! Unsaturated limit
      d2Sl_dPc2 = this%d2Sl_dPc2max
    else                                ! Exponential extension
      d2Sl_dPc2 = -this%beta_rec/Pc**2
    end if
  else                                  ! Ordinary VG domain
    call this%D2SlDPc2Inline(Pc, d2Sl_dPc2)
  end if
end subroutine SFVGexpnD2SlDPc2

! **************************************************************************** !
! Linear Van Genuchten Extension
! **************************************************************************** !

function SFVGLCPCCtor(alpha,m,Sr,Sgt,rpf,Pcmax) result (new)

  implicit none

  class(sf_VG_line_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sgt, Pcmax
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,Sgt,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%SetPcmax(Pcmax) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function SFVGLCPCCtor

! **************************************************************************** !

function SFVGlineSetPcmax(this,Pcmax) result (error)

  implicit none

  class(sf_VG_line_type), intent(inout) :: this
  PetscReal, intent(in) :: Pcmax
  PetscInt :: error
  PetscReal :: Sa, Sb, Pe, dPj_dSj

  ! Solve the nonlinear system to satisfy C1 continuity using bisection method

  ! Set lower saturation bracket Set where P(Sa) = Pcmax
  call this%SlInline(Pcmax,Sa,dPj_dSj)

  ! Set upper saturation limit (Sb) to be the the inflection point
  Sb = this%SlInflection()

  ! Confirm Pcmax is above minimum extrapolating from inflection point
  call this%PcInline(Sb,Pe,dPj_dSj)
  if (Pcmax > Pe - dPj_dSj*Sb) then
    error = 0
    this%Pcmax_designated = PETSC_TRUE
    this%Pcmax = Pcmax

    do while (Sb-Sa > epsilon(this%Sj)) ! Tolerance interval epsilon
      this%Sj = (Sa+Sb)/2d0             ! Bisect bracket
      call this%PcInline(this%Sj,this%Pj,this%dPj_dSj)

      ! Residual error = Pcmax + dPj_dSj * Sj - Pf(Sj)
      Pe = Pcmax + this%dPj_dSj*this%Sj - this%Pj
      if (Pe < 0d0) then ! Error on side a is always negative on LHS of VG
        Sa = this%Sj
      else
        Sb = this%Sj
      end if
    end do
    this%dSj_dPj = 1d0 / this%dPj_dSj
  else
    error = 1
  end if
end function SFVGlineSetPcmax

! **************************************************************************** !

function SFVGLNOCCtor(alpha,m,Sr,Sgt,rpf,Sj) result (new)

  implicit none

  class(sf_VG_line_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sgt, Sj
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,Sgt,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%SetSj(Sj) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function SFVGLNOCCtor

! **************************************************************************** !

function SFVGlineSetSj(this,Sj) result (error)
  class(sf_VG_line_type), intent(inout) :: this
  PetscReal, intent(in) :: Sj
  PetscInt :: error

  if (Sj > this%Sr .AND. Sj < 1d0) then
    error = 0
    this%Pcmax_designated = PETSC_FALSE
    this%Sj = Sj

    ! Find the value and slope at the piecewise junction point and extrapolate
    call this%PcInline(Sj,this%Pj,this%dPj_dSj)
    this%Pcmax = this%Pj - this%dPj_dSj * Sj
    this%dSj_dPj = 1d0 / this%dPj_dSj
  else
    error = 1
  end if
end function SFVGlineSetSj

! **************************************************************************** !

pure subroutine SFVGlinePc(this, Sl, Pc, dPc_dSl)
  class(sf_VG_line_type), intent(in) :: this
  PetscReal, intent(in)   :: Sl
  PetscReal, intent(out)  :: Pc, dPc_dSl

  if (Sl < this%Sj) then                ! Unsaturated domain
    dPc_dSl = this%dPj_dSj
    if (Sl <= 0d0) then                 ! Unsaturated limit
      Pc = this%Pcmax
    else                                ! Linear extension
      Pc = this%Pcmax + this%dPj_dSj*Sl
    end if
  else if (Sl > Sl_max) then            ! Saturated limit
    Pc = 0d0
    dPc_dSl = this%dPc_dSl_max
  else                                  ! Ordinary VG domain
    call this%PcInline(Sl, Pc, dPc_dSl)
  end if
end subroutine SFVGlinePc

! **************************************************************************** !

pure subroutine SFVGlineSl(this, Pc, Sl, dSl_dPc)

  implicit none

  class(sf_VG_line_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sl, dSl_dPc

  if (Pc <= 0d0) then                   ! Saturated limit
    Sl = 1d0
    dSl_dPc = this%dSl_dPcmin
  else if (Pc >= this%Pj) then          ! Unsaturated domain
    dSl_dPc = this%dSj_dPj
    if (Pc >= this%Pcmax) then          ! Unsaturated limit
      Sl = 0d0
    else                                ! Linear extensions
      Sl = (Pc - this%Pcmax) * this%dSj_dPj
    end if
  else                                  ! Ordinary VG domain
    call this%SlInline(Pc, Sl, dSl_dPc)
  end if
end subroutine SFVGlineSl

! **************************************************************************** !

pure subroutine SFVGlineD2SlDPc2(this, Pc, d2Sl_dPc2)
  class(sf_VG_line_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2

  if (pc <= 0d0) then                   ! Saturated limit
    d2Sl_dPc2 = this%d2Sl_dPc2min
  else if (Pc >= this%Pj) then          ! Unsaturated domain
    d2Sl_dPc2 = 0d0
  else                                  ! Ordinary VG domain
    call this%D2SlDPc2Inline(Pc, d2Sl_dPc2)
  end if
end subroutine SFVGlineD2SlDPc2

! **************************************************************************** !
! Quadratic Van Genuchten Extension
! **************************************************************************** !

function SFVGquadCtor(alpha, m, Sr, Sgt, rpf, Pcmax, Sj) result (new)
  class(sf_VG_quad_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Pcmax, Sgt, Sj
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,Sgt,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%SetQuad(Pcmax, Sj) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function SFVGquadCtor

! **************************************************************************** !

function SFVGquadSetPcmax(this, Pcmax) result (error)

  implicit none

  class(sf_VG_quad_type), intent(inout) :: this
  PetscReal, intent(in) :: Pcmax
  PetscInt :: error

  error = this%SetQuad(Pcmax, this%Sj)
end function SFVGquadSetPcmax

! **************************************************************************** !

function SFVGquadSetSj(this, Sj) result (error)

  implicit none

  class(sf_VG_quad_type), intent(inout) :: this
  PetscReal, intent(in) :: Sj
  PetscInt :: error

  error = this%SetQuad(this%Pcmax, this%Sj)
end function SFVGquadSetSj

! **************************************************************************** !

function SFVGquadSetQuad(this, Pcmax, Sj) result (error)

  implicit none

  class(sf_VG_quad_type), intent(inout) :: this
  PetscReal, intent(in) :: Pcmax, Sj
  PetscInt :: error
  PetscReal :: Pj, dPj_dSj, dPc_dSj, A, B, Sl_extremum, Sl_qa, Sl_cq

  ! Must satisfy Pc(0) = Pcmax and C1 continuity at Sj

  !                  C = Pc(0) = Pcmax
  ! A*Sj**2 + B*Sj + C = Pc(Sj)
  ! A*Sj*2  + B        = dPc_dSl(Sj)

  ! Let dPc_dSj = (Pj - Pcmax) / Sj

  ! A = (dPj_dSj - dPc_dSj) / Sj
  ! B = 2*dPc_dSj - dPj_dSj

  if (Sj <= this%Sr) then ! Sj is below residual saturation
    error = 2
  else
    call this%PcInline(Sj, Pj, dPj_dSj)
    dPc_dSj = (Pj - Pcmax)/Sj
    A = (dPj_dSj - dPc_dSj) / Sj
    B = 2d0*dPc_dSj - dPj_dSj
    Sl_extremum = -B/(2d0*A)

    if (Sl_extremum > 0d0 .and. Sl_extremum < Sj) then
      error = 1 ! Parabola is not monotonic over extension. Return error
    else
      error = 0
      this%Pcmax_designated = PETSC_TRUE
      this%Pcmax = Pcmax ! Pcmax = C
      this%Pj = Pj
      this%Sj = Sj

      this%A = A
      this%B = B
      this%dSl_dPcmax = 1d0 / B

      ! Select the correct branch for inverse solutions
      this%qa_branch = PETSC_TRUE
      Sl_qa = this%SlRoot(Pj)
      this%qa_branch = PETSC_FALSE
      Sl_cq = this%SlRoot(Pj)

      this%qa_branch = (abs(Sl_qa - Sj) < abs(Sl_cq - Sj))
    end if
  end if
end function SFVGquadSetQuad

! **************************************************************************** !

pure subroutine SFVGquadPc(this, Sl, Pc, dPc_dSl)

  implicit none

  class(sf_VG_quad_type), intent(in) :: this
  PetscReal, intent(in)   :: Sl
  PetscReal, intent(out)  :: Pc, dPc_dSl

  if (Sl < this%Sj) then                ! Unsaturated domain
    if (Sl <= 0d0) then                 ! Unsaturated limit
      Pc = this%Pcmax
      dPc_dSl = this%B
    else                                ! Quadratic extension
      Pc = (this%A*Sl + this%B)*Sl + this%Pcmax
      dPc_dSl = 2d0*this%A*Sl + this%B
    end if
  else if (Sl > Sl_max) then            ! Saturated limit
    Pc = 0d0
    dPc_dSl = this%dPc_dSl_max
  else                                  ! Ordinary VG domain
    call this%PcInline(Sl, Pc, dPc_dSl)
  end if
end subroutine SFVGquadPc

! **************************************************************************** !

pure subroutine SFVGquadSl(this, Pc, Sl, dSl_dPc)

  implicit none

  class(sf_VG_quad_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sl, dSl_dPc

  if (Pc <= 0d0) then                   ! Saturated limit
    Sl = 1d0
    dSl_dPc = this%dSl_dPcmin
  else if (Pc >= this%Pj) then          ! Unsaturated domain
    if (Pc >= this%Pcmax) then          ! Unsaturated limit
      Sl = 0d0
      dSl_dPc = this%dSl_dPcmax
    else                                ! Quadratic extension
      Sl = this%SlRoot(Pc)
      dSl_dPc = 1d0 / (2d0*this%A*Sl + this%B)
    end if
  else                                  ! Ordinary VG domain
    call this%SlInline(Pc, Sl, dSl_dPc)
  end if
end subroutine SFVGquadSl

! **************************************************************************** !

pure subroutine SFVGquadD2SlDPc2(this, Pc, d2Sl_dPc2)

  implicit none

  class(sf_VG_quad_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2
  PetscReal :: Sl

  if (pc <= 0d0) then                   ! Saturated limit
    d2Sl_dPc2 = this%d2Sl_dPc2min
  else if (Pc >= this%Pj) then          ! Unsaturated domain
    if (Pc >= this%Pcmax) then          ! Unsaturated limit
      d2Sl_dPc2 = -1d0 / this%B**3
    else                                ! Quadratic extension
      Sl = this%SlRoot(Pc)
      d2Sl_dPc2 = -2d0*this%A / (2d0*this%A*Sl + this%B)**3
    end if
  else                                  ! Ordinary VG domain
    call this%D2SlDPc2Inline(Pc, d2Sl_dPc2)
  end if
end subroutine SFVGquadD2SlDPc2

! **************************************************************************** !

pure function SFVGquadSlRoot(this,Pc) result (Sl)

  implicit none

  class(sf_VG_quad_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal :: Sl
  PetscReal :: C, Q

  C = this%Pcmax - Pc
  Q = -0.5d0*(this%B+sign(sqrt(this%B**2-4d0*this%A*C), this%B))
  if (this%qa_branch) then
    Sl = Q/this%A
  else
    Sl = C/Q
  end if
end function SFVGquadSlRoot

! **************************************************************************** !
! Abstract van Genuchten Relative Permeability Methods
! **************************************************************************** !

pure function RPFVGGetM(this) result (m)

  implicit none

  class(rpf_VG_type), intent(in) :: this
  PetscReal :: m
  m = this%m
end function RPFVGGetM

! **************************************************************************** !

subroutine RPFVGInit(this)

  implicit none

  class(rpf_VG_type) :: this
  ! This method is intentionally left blank.
end subroutine RPFVGInit

! **************************************************************************** !

subroutine RPFVGRelativePermeability(this, liquid_saturation, &
                                       relative_permeability, dkr_sat, option)

  implicit none

  class(rpf_VG_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability, dkr_sat
  type(option_type), intent(inout) :: option

  call this%Kr(liquid_saturation, relative_permeability, dkr_sat)
end subroutine RPFVGRelativePermeability

! **************************************************************************** !
! Abstract Liquid van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPFVGliqConfigure(this, m, Sr) result (error)

  implicit none

  class(rpf_VG_liq_type), intent(inout) :: this
  PetscReal, intent(in) :: m, Sr
  PetscInt :: error
  PetscReal :: Kr

  error = 0
  if (m <= 0d0 .OR. m >= 1d0) error = error + 1
  if (Sr < 0d0 .OR. Sr >= 1d0) error = error + 2

  if (error == 0) then
    this%m = m
    this%m_rec = 1d0 / m

    this%Sr = Sr
    this%dSe_dSl = 1d0 / (1d0 - Sr)

    call this%KrInline(Sl_max,Kr,this%dKr_dSl_max)
    this%dKr_dSl_max = (1d0 - Kr) / epsilon(Sl_max)

    this%Sl_min = Sr + epsilon(Sr)

    this%analytical_derivative_available = PETSC_TRUE
  end if
end function RPFVGliqConfigure

! **************************************************************************** !

function RPFVGliqSetM(this, m) result (error)

  implicit none

  class(rpf_VG_liq_type), intent(inout) :: this
  PetscReal, intent(in) :: m
  PetscInt :: error
  error = this%Configure(m, this%Sr)
end function RPFVGliqSetM

! **************************************************************************** !

pure subroutine RPFVGliqKr(this, Sl, Kr, dKr_dSl)

  implicit none

  class(rpf_VG_liq_type), intent(in) :: this
  PetscReal, intent(in) :: Sl
  PetscReal, intent(out) :: Kr, dKr_dSl

  if (Sl > Sl_max) then                 ! Saturated limit
    Kr = 1d0
    dKr_dSl = this%dKr_dSl_max
  else if (Sl < this%Sl_min) then       ! Unsaturated limit
    Kr = 0d0
    dKr_dSl = 0d0
  else                                  ! Ordinary domain
    call this%KrInline(Sl, Kr, dKr_dSl)
  end if
end subroutine RPFVGliqKr

! **************************************************************************** !
! Liquid Mualem - van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPFMVGliqCtor(m, Sr) result (new)

  implicit none

  class(rpf_MVG_liq_type), pointer :: new
  PetscReal, intent(in) :: m, Sr
  PetscInt :: error

  allocate(new)
  error = new%Configure(m,Sr)
  if (error == 0) then
    nullify(new%poly)
  else
    deallocate(new)
    nullify(new)
  end if
end function RPFMVGliqCtor

! **************************************************************************** !

pure subroutine RPFMVGliqKrInline(this, Sl, Kr, dKr_dSl)

  implicit none

  class(rpf_MVG_liq_type), intent(in) :: this
  PetscReal, intent(in) :: Sl
  PetscReal, intent(out) :: Kr, dKr_dSl
  PetscReal :: Se, Se_mrt, Se_mrt_comp, Se_mrt_comp_m, f

  Se = (Sl- this%Sr) * this%dSe_dSl
  Se_mrt  = Se**this%m_rec
  Se_mrt_comp = 1d0 - Se_mrt
  Se_mrt_comp_m = Se_mrt_comp**this%m
  f = 1d0 - Se_mrt_comp_m

  Kr = sqrt(Se)*f*f
  if (f > 0.d0) then
    dKr_dSl = this%dSe_dSl * Kr / Se * &
             (0.5d0 + 2d0*Se_mrt*Se_mrt_comp_m/(f*Se_mrt_comp))
  else
    dKr_dSl = 0.d0
  endif
end subroutine RPFMVGliqKrInline

! **************************************************************************** !
! Liquid Burdine - van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPFBVGliqCtor(m, Sr) result (new)

  implicit none

  class(rpf_BVG_liq_type), pointer :: new
  PetscReal, intent(in) :: m, Sr
  PetscInt :: error

  allocate(new)
  error = new%Configure(m, Sr)

  if (error == 0) then
    nullify(new%poly)
  else
    deallocate(new)
    nullify(new)
  end if
end function RPFBVGliqCtor

! **************************************************************************** !

pure subroutine RPFBVGliqKrInline(this, Sl, Kr, dKr_dSl)

  implicit none

  class(rpf_BVG_liq_type), intent(in) :: this
  PetscReal, intent(in) :: Sl
  PetscReal, intent(out) :: Kr, dKr_dSl
  PetscReal :: Se, Se_mrt, Se_mrt_comp, Se_mrt_comp_m, Kr_Se

  Se = (Sl - this%Sr) * this%dSe_dSl
  Se_mrt = Se**this%m_rec
  Se_mrt_comp = 1d0 - Se_mrt
  Se_mrt_comp_m = Se_mrt_comp**this%m
  Kr_Se = Se*(1d0 - Se_mrt_comp_m)

  Kr = Kr_Se*Se
  dKr_dSl = this%dSe_dSl * (2d0*Kr_Se + Se*Se_mrt*Se_mrt_comp_m/Se_mrt_comp)
end subroutine RPFBVGliqKrInline

! **************************************************************************** !
! Abstract Gas van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPFVGgasConfigure(this, m, Sr, Sgr) result (error)

  implicit none

  class(rpf_VG_gas_type), intent(inout) :: this
  PetscReal, intent(in) :: m, Sr, Sgr
  PetscInt :: error

  error = 0
  if (m <= 0d0 .OR. m >= 1d0) error = error + 1
  if (Sr < 0d0 .OR. Sgr < 0d0 .OR. Sr + Sgr >= 1d0) error = error + 2

  if (error == 0) then
    this%m = m
    this%m_rec = 1d0 / m
    this%mx2 = 2d0 * m

    this%Sr = Sr
    this%Sgr = Sgr
    this%dSe_dSl = 1d0 / (1d0 - Sr - Sgr)

    this%Sl_max = 1d0 - Sgr - epsilon(Sgr)

    this%analytical_derivative_available = PETSC_TRUE
  end if
end function RPFVGgasConfigure

! **************************************************************************** !

function RPFVGgasSetM(this, m) result (error)

  implicit none

  class(rpf_VG_gas_type), intent(inout) :: this
  PetscReal, intent(in) :: m
  PetscInt :: error
  error = this%Configure(m, this%Sr, this%Sgr)
end function RPFVGgasSetM

! **************************************************************************** !

pure subroutine RPFVGgasKr(this, Sl, Kr, dKr_dSl)

  implicit none

  class(rpf_VG_gas_type), intent(in) :: this
  PetscReal, intent(in) :: Sl
  PetscReal, intent(out) :: Kr, dKr_dSl

  if (Sl > this%Sl_max) then            ! Saturated limit
    Kr = 0d0
    dKr_dSl = 0d0
  else if (Sl < this%Sr) then           ! Unsaturated limit
    Kr = 1d0
    dKr_dSl = 0d0                       ! Cusp at limit
  else                                  ! Ordinary domain
    call this%KrInline(Sl, Kr, dKr_dSl)
  end if
end subroutine RPFVGgasKr

! **************************************************************************** !
! Gas Mualem - van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPFMVGgasCtor(m, Sr, Sgr) result (new)

  implicit none

  class(rpf_MVG_gas_type), pointer :: new
  PetscReal, intent(in) :: m, Sr, Sgr
  PetscInt :: error

  allocate(new)
  error = new%Configure(m, Sr, Sgr)

  if (error == 0) then
    nullify(new%poly)
  else
    deallocate(new)
    nullify(new)
  end if
end function RPFMVGgasCtor

! **************************************************************************** !

pure subroutine RPFMVGgasKrInline(this, Sl, Kr, dKr_dSl)

  implicit none

  class(rpf_MVG_gas_type), intent(in) :: this
  PetscReal, intent(in) :: Sl
  PetscReal, intent(out) :: Kr, dKr_dSl
  PetscReal :: Se, Se_comp, Se_mrt, Se_mrt_comp

  Se = (Sl - this%Sr) * this%dSe_dSl
  Se_comp = 1.d0 - Se
  Se_mrt = Se**this%m_rec
  Se_mrt_comp = 1d0 - Se_mrt

  Kr = sqrt(Se_comp)*Se_mrt_comp**this%mx2
  dKr_dSl = - this%dSe_dSl * Kr * (0.5d0/Se_comp + 2d0*Se_mrt/(Se*Se_mrt_comp))
end subroutine RPFMVGgasKrInline

! **************************************************************************** !
! Gas Burdine - van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPFBVGgasCtor(m, Sr, Sgr) result (new)

  implicit none

  class(rpf_BVG_gas_type), pointer :: new
  PetscReal, intent(in) :: m, Sr, Sgr
  PetscInt :: error

  allocate(new)
  error = new%Configure(m, Sr, Sgr)

  if (error == 0) then
    nullify(new%poly)
  else
    deallocate(new)
    nullify(new)
  end if
end function RPFBVGgasCtor

! **************************************************************************** !

pure subroutine RPFBVGgasKrInline(this, Sl, Kr, dKr_dSl)

  implicit none

  class(rpf_BVG_gas_type), intent(in) :: this
  PetscReal, intent(in) :: Sl
  PetscReal, intent(out) :: Kr, dKr_dSl
  PetscReal :: Se, Se_comp, Se_mrt, Se_mrt_comp, Kr_Se_comp

  Se = (Sl - this%Sr) * this%dSe_dSl
  Se_comp = 1d0 - Se
  Se_mrt = Se**this%m_rec
  Se_mrt_comp = 1d0 - Se_mrt
  Kr_Se_comp = Se_comp*Se_mrt_comp**this%m

  Kr = Kr_Se_comp*Se_comp
  dKr_dSl = -this%dSe_dSl * (2d0*Kr_Se_comp + Kr*Se_mrt/(Se*Se_mrt_comp))
end subroutine RPFBVGgasKrInline

! **************************************************************************** !

end module

