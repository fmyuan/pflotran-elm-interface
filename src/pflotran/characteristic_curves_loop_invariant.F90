module Characteristic_Curves_loop_invariant_module
#include "petsc/finclude/petscsys.h"

use petscsys ! Necessary for PETSC_TRUE / PETSC_FALSE
use Characteristic_Curves_Base_module   ! Needed to define base type
use Characteristic_Curves_Common_module ! Needed to inherent VG type
use Option_module ! Needed for Verify and unused arguments in Pc and Sl
use PFLOTRAN_Constants_module ! Needed for overridden keyword "name" in verify
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
!     |-->SF_VG_type        VG type with loop-invariant parameters, no extension
!         |
!         |-->SF_VG_extn_type*     Abstract VG type for extensions
!             |
!             |-->SF_VG_cons_type  Constant    extension type
!             |-->SF_VG_expn_type  Exponential extension type
!             |-->SF_VG_line_type  Linear      extension type
!             |-->SF_VG_quad_type  Quadratic   extension type
!
! Objects of SF_VG_type are created and initialized with the constructor:
! SF_VG_type       => SF_VG_ctor(unsat_ext,alpha,m,Sr,rpf,Pcmax,Sj)
!
! Warning, Pcmax and Sr are unprotected as Fortran will not extend private
! scope to child classed in separate modules. 
!
! **************************************************************************** !
! VG Saturation Function constructors
! **************************************************************************** !

  public  :: SF_VG_ctor         ! Generic constructor
  private :: SF_VG_NEVG_ctor, & ! No capped capillary pressure
             SF_VG_FCPC_ctor, & ! Flat, specified maximum
             SF_VG_FNOC_ctor, & ! Flat, specified junction
             SF_VG_ECPC_ctor, & ! Exponential, specified maximum
             SF_VG_ENOC_ctor, & ! Exponential, specified junction
             SF_VG_LCPC_ctor, & ! Linear, specified maximum
             SF_VG_LNOC_ctor, & ! Linear, specified junction
             SF_VG_QUAD_ctor    ! Quadratic, specified maximum and junction

! **************************************************************************** !
! Van Genuchten Saturation Function type defintions
! **************************************************************************** !
  type, public, extends(sat_func_VG_type) :: SF_VG_type
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
    procedure, public  :: Init                  => SF_VG_Init
    procedure, public  :: CapillaryPressure     => SF_VG_CapillaryPressure
    procedure, public  :: Saturation            => SF_VG_Saturation
    procedure, public  :: D2SatDP2              => SF_VG_D2SatDP2
!   procedure, public  :: Verify
! Common VG methods
    procedure, private :: Configure             => SF_VG_Configure
    procedure, private :: Pc_inline             => SF_VG_Pc_inline
    procedure, private :: Sl_inline             => SF_VG_Sl_inline
    procedure, private :: d2Sl_dPc2_inline      => SF_VG_d2Sl_dPc2_inline
    procedure, private :: Sl_inflection         => SF_VG_Sl_inflection
! No-extension mutator methods
    procedure, private :: Set_alpha             => SF_VG_Set_alpha
    procedure, private :: Set_m                 => SF_VG_Set_m
! Internal no-extension pure methods
    procedure, private :: Pc                    => SF_VG_Pc
    procedure, private :: Sl                    => SF_VG_Sl
    procedure, private :: d2Sl_dPc2             => SF_VG_d2Sl_dPc2
  end type

! VG Unsaturated Extension Abstract type
    type, abstract, extends(SF_VG_type) :: SF_VG_extn_type
      private
        PetscReal :: Pj, Sj           ! Piecewise junction point
        PetscBool :: Pcmax_designated ! Flag for designated Pcmax
    contains
      ! Common mutator methods
      procedure, public :: Set_alpha            => SF_VG_extn_Set_alpha
      procedure, public :: Set_m                => SF_VG_extn_Set_m
      ! Mutator prototype
      procedure(Set_Pcmax), deferred, public    :: Set_Pcmax
      procedure(Set_Sj   ), deferred, public    :: Set_Sj
    end type

! VG Constant Unsaturated Extension
      type, public, extends(SF_VG_extn_type) :: SF_VG_cons_type
        private
          PetscReal :: dSj_dPj   ! Non-zero derivative at cusp
          PetscReal :: d2Sj_dPj2 ! Non-zero 2nd derivative at cusp
      contains
        ! Public mutator methods
        procedure, public  :: Set_Pcmax         => SF_VG_cons_Set_Pcmax
        procedure, public  :: Set_Sj            => SF_VG_cons_Set_Sj
        ! Internal methods for constant extensions
        procedure, private :: Pc                => SF_VG_cons_Pc
        procedure, private :: Sl                => SF_VG_cons_Sl
        procedure, private :: d2Sl_dPc2         => SF_VG_cons_d2Sl_dPc2
      end type

! VG Exponential Unsaturated Extension
      type, public, extends(SF_VG_extn_type) :: SF_VG_expn_type
        private
          PetscReal :: beta, beta_rec   ! Exponential coefficient and reciprocal
          PetscReal :: dPcmax_dSl       ! Unsaturated limits
          PetscReal :: dSl_dPcmax       ! Derivative at unsaturated limit
          PetscReal :: d2Sl_dPc2max     ! 2nd derivative at unsaturated limit
      contains
        ! Public mutator methods
        procedure, public  :: Set_Pcmax         => SF_VG_expn_Set_Pcmax
        procedure, public  :: Set_Sj            => SF_VG_expn_Set_Sj
        ! Internal methods for exponential extensions
        procedure, private :: Pc                => SF_VG_expn_Pc
        procedure, private :: Sl                => SF_VG_expn_Sl
        procedure, private :: d2Sl_dPc2         => SF_VG_expn_d2Sl_dPc2
      end type

! VG Linear Unsaturated Extension
      type, public, extends(SF_VG_extn_type) :: SF_VG_line_type
        private
          PetscReal :: dPj_dSj, dSj_dPj ! Linear coefficients
      contains
        ! Public mutator methods
        procedure, public  :: Set_Pcmax         => SF_VG_line_Set_Pcmax
        procedure, public  :: Set_Sj            => SF_VG_line_Set_Sj
        ! Internal methods for linear extensions
        procedure, private :: Pc                => SF_VG_line_Pc
        procedure, private :: Sl                => SF_VG_line_Sl
        procedure, private :: d2Sl_dPc2         => SF_VG_line_d2Sl_dPc2
      end type

! VG Quadratic Unsaturated Extension
      type, public, extends(SF_VG_extn_type) :: SF_VG_quad_type
        private
          PetscReal :: A, B            ! Quadratic coefficients
          PetscReal :: dSl_dPcmax
          PetscReal :: d2Sl_dPc2max
          PetscBool :: qa_branch       ! Select inverse solution
      contains
        ! Public mutator methods
        procedure, public  :: Set_Pcmax         => SF_VG_quad_Set_Pcmax
        procedure, public  :: Set_Sj            => SF_VG_quad_Set_Sj
        procedure, public  :: Set_quad          => SF_VG_quad_Set_quad
        ! Internal methods for quadratic extensions
        procedure, private :: Pc                => SF_VG_quad_Pc
        procedure, private :: Sl                => SF_VG_quad_Sl
        procedure, private :: d2Sl_dPc2         => SF_VG_quad_d2Sl_dPc2
        procedure, private :: Sl_root           => SF_VG_quad_Sl_root
      end type

! **************************************************************************** !
! Saturation Function Prototype Block
! **************************************************************************** !

interface
  function Set_Pcmax(this, Pcmax) result (error)
    import :: SF_VG_extn_type
    class(SF_VG_extn_type), intent(inout) :: this
    PetscReal, intent(in) :: Pcmax
    PetscInt :: error
  end function

! **************************************************************************** !

  function Set_Sj(this, Sj) result (error)
    import :: SF_VG_extn_type
    class(SF_VG_extn_type), intent(inout) :: this
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
! |-->RPF_VG_type*                   Abstract base type
!    |
!    |-->RPF_VG_liq_type*            Abstract liquid base type
!    |   |
!    |   |-->RPF_MVG_liq_type        Mualem  - van Genuchten - Liquid
!    |   |-->RPF_BVG_liq_type        Burdine - van Genuchten - Liquid
!    |
!    |-->RPF_VG_gas_type*            Abstract gas base type
!        |
!        |-->RPF_MVG_gas_type        Mualem  - van Genuchten - Gas
!        |-->RPF_BVG_gas_type        Burdine - van Genuchten - Gas
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
! VG RPF Contructors
! **************************************************************************** !

! Relative permeability function constructors
  public :: RPF_MVG_Liq_ctor, &
            RPF_MVG_Gas_ctor, &
            RPF_BVG_Liq_ctor, &
            RPF_BVG_Gas_ctor

type, abstract, extends(rel_perm_func_base_type) :: RPF_VG_type
! public
!   PetscReal :: Sr         ! Base   - Residual saturation
  private
    PetscReal :: m, m_rec   ! Van Genuchten m exponent and reciprocal
    PetscReal :: dSe_dSl    ! Effective saturation span reciprocal
contains
  ! Overloaded base methods
  procedure, public :: Init                 => RPF_VG_Init
  procedure, public :: RelativePermeability => RPF_VG_RelativePermeability
  ! Public accessor/mutator methods
  procedure, public                     :: Get_m                => RPF_VG_Get_m
  procedure(Set_m)  , deferred, public  :: Set_m
  ! Private deferred methods called by RelativePermeability
  procedure(Calc_Kr), deferred, private :: Kr
  procedure(Calc_Kr), deferred, private :: Kr_inline
end type

! **************************************************************************** !

  type, abstract, extends(RPF_VG_type) :: RPF_VG_liq_type
    private
      PetscReal :: dKr_dSl_max ! Finite difference derivative at saturation
      PetscReal :: Sl_min      ! Unsaturated limit to avoid NaN
  contains
    procedure, public  :: Set_m                 => RPF_VG_liq_Set_m
    procedure, private :: Configure             => RPF_VG_liq_Configure
    procedure, private :: Kr                    => RPF_VG_liq_Kr
  end type

! **************************************************************************** !

    type, public, extends(RPF_VG_liq_type) :: RPF_MVG_liq_type
    contains
      procedure, private :: Kr_inline           => RPF_MVG_liq_Kr_inline
    end type

! **************************************************************************** !

    type, public, extends(RPF_VG_liq_type) :: RPF_BVG_liq_type
    contains
      procedure, private :: Kr_inline           => RPF_BVG_liq_Kr_inline
    end type

! **************************************************************************** !

  type, abstract, extends(RPF_VG_type) :: RPF_VG_gas_type
  private
    PetscReal :: mx2        ! m times 2
    PetscReal :: Sgr        ! Gas residual saturation
    PetscReal :: Sl_max     ! Saturated limit reduced by Sgr
  contains
    procedure, public  :: Set_m                  => RPF_VG_gas_Set_m
    procedure, private :: Configure              => RPF_VG_gas_Configure
    procedure, private :: Kr                     => RPF_VG_gas_Kr
  end type

! **************************************************************************** !

    type, public, extends(RPF_VG_gas_type) :: RPF_MVG_gas_type
    contains
      procedure, private :: Kr_inline           => RPF_MVG_gas_Kr_inline
    end type

! **************************************************************************** !

    type, public, extends(RPF_VG_gas_type) :: RPF_BVG_gas_type
    contains
      procedure, private :: Kr_inline           => RPF_BVG_gas_Kr_inline
    end type

! **************************************************************************** !
! Relative Permeability Function Prototype Block
! **************************************************************************** !

interface
  function Set_m(this, m) result (error)
    import :: RPF_VG_type
    class(RPF_VG_type), intent(inout) :: this
    PetscReal, intent(in) :: m
    PetscInt :: error
  end function

! **************************************************************************** !

  pure subroutine Calc_Kr(this, Sl, Kr, dKr_dSl)
    import :: RPF_VG_type
    class(RPF_VG_type), intent(in) :: this
    PetscReal, intent(in)  :: Sl
    PetscReal, intent(out) :: Kr, dKr_dSl
  end subroutine
end interface

contains

! **************************************************************************** !
! VG Saturation Function Methods
! **************************************************************************** !

function SF_VG_ctor(unsat_ext, alpha, m, Sr, vg_rpf_opt, Pcmax, Sj) result (new)
 class(SF_VG_type), pointer :: new
 character(*), intent(in) :: unsat_ext
 PetscReal, intent(in) :: alpha, m, Sr, Pcmax, Sj
 PetscInt, intent(in) :: vg_rpf_opt

! This function returns a the van Genuchten saturation function object using
! the correct constructor method

  select case (unsat_ext)
  case ('NONE') ! No extension
    new => SF_VG_NEVG_ctor(alpha,m,Sr,vg_rpf_opt)
  case ('FCPC') ! Flat specified cap
    new => SF_VG_FCPC_ctor(alpha,m,Sr,vg_rpf_opt,Pcmax)
  case ('FNOC') ! Flat specificed junction
    new => SF_VG_FNOC_ctor(alpha,m,Sr,vg_rpf_opt,Sj)
  case ('ECPC') ! Exponential specified cap
    new => SF_VG_ECPC_ctor(alpha,m,Sr,vg_rpf_opt,Pcmax)
  case ('ENOC') ! Exponential specified junction
    new => SF_VG_ENOC_ctor(alpha,m,Sr,vg_rpf_opt,Sj)
  case ('LCPC') ! Linear specified cap
    new => SF_VG_LCPC_ctor(alpha,m,Sr,vg_rpf_opt,Pcmax)
  case ('LNOC') ! Linear specified junction
    new => SF_VG_LNOC_ctor(alpha,m,Sr,vg_rpf_opt,Sj)
  case ('QUAD') ! Quadratic specified cap and junction
    new => SF_VG_quad_ctor(alpha,m,Sr,vg_rpf_opt,Pcmax,Sj)
  case default
    nullify(new)
  end select
end function

! **************************************************************************** !

subroutine SF_VG_init(this)
  class(SF_VG_type) :: this
  ! This method is intentionally left blank.
end subroutine

! **************************************************************************** !

subroutine SF_VG_CapillaryPressure(this, liquid_saturation, &
                                   capillary_pressure, dPc_dSatl, option)
  class(SF_VG_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option
  
  call this%Pc(liquid_saturation, capillary_pressure, dPc_dSatl)
end subroutine

! **************************************************************************** !

subroutine SF_VG_Saturation(this, capillary_pressure, &
                            liquid_saturation, dsat_dpres, option)
  class(SF_VG_type) :: this
  PetscReal, intent(in)  :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  call this%Sl(capillary_pressure, liquid_saturation, dsat_dpres)
end subroutine

! **************************************************************************** !

subroutine SF_VG_D2SatDP2(this,Pc, d2s_dp2, option)
  class(SF_VG_type) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option

  call this%d2Sl_dPc2(Pc, d2s_dp2)
end subroutine

! **************************************************************************** !

function SF_VG_Configure(this, alpha, m, Sr, rpf) result (error)
  ! Configure loop-invariant parameters common to all loop-invariant types
  class(SF_VG_type) :: this
  PetscReal, intent(in) :: alpha,m,Sr
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
    call this%Pc_inline(Sl_max, Pc1, dPc_dSl)
    call this%Pc_inline(Sl_max-epsilon(Sl_max), Pc2, dPc_dSl)
    this%dPc_dSl_max = -Pc1 / epsilon(Sl_max)
    this%dSl_dPcmin = -epsilon(Sl_max) / Pc1
    this%d2Sl_dPc2min = epsilon(Sl_max)*(2d0-Pc2/Pc1)/Pc1**2
  end if
end function

! **************************************************************************** !

pure subroutine SF_VG_Pc_inline(this, Sl, Pc, dPc_dSl)
  class(SF_VG_type), intent(in) :: this
  PetscReal, intent(in) :: Sl
  PetscReal, intent(out) :: Pc, dPc_dSl
  PetscReal :: Se, Se_mrtrec, aPc_n

  Se = (Sl - this%Sr) * this%dSe_dSl
  Se_mrtrec = Se**this%m_nrec
  aPc_n = Se_mrtrec - 1d0

  Pc = this%a_rec * aPc_n**this%n_rec
  dPc_dSl = (this%dSe_mndSl*Pc) / (Se/Se_mrtrec-Se)
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_Sl_inline(this, Pc, Sl, dSl_dPc)
  class(SF_VG_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: Sl, dSl_dPc
  PetscReal :: aPc_n, Se_mrtrec, Se

  aPc_n = (this%alpha*Pc)**this%n
  Se_mrtrec = aPc_n + 1d0
  Se = Se_mrtrec**(-this%m)

  Sl = this%Sr + Se * this%Sl_span
  dSl_dPc = (Se/Se_mrtrec-Se) / (this%dSe_mndSl*Pc)
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_d2Sl_dPc2_inline(this, Pc, d2Sl_dPc2)
  class(SF_VG_type), intent(in) :: this
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
end subroutine

! **************************************************************************** !

pure function SF_VG_Sl_inflection(this) result (Sl)
  class(SF_VG_type), intent(in) :: this
  PetscReal :: Se_mrtrec, Se, Sl

  Se_mrtrec = (1d0+this%m)/(this%n_rec+this%m)
  Se = Se_mrtrec**(-this%m)
  Sl = this%Sr + Se*this%Sl_span
end function

! **************************************************************************** !
! No-extension VG Saturation Function Methods
! **************************************************************************** !

function SF_VG_NEVG_ctor(alpha,m,Sr,rpf) result (new)
  class(SF_VG_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  end if
  new%Pcmax = huge(new%Pcmax)
end function

! **************************************************************************** !

function SF_VG_Set_alpha(this,alpha) result (error)
  class(SF_VG_type), intent(inout) :: this
  PetscReal, intent(in) :: alpha
  PetscInt :: error
  error = this%Configure(alpha,this%m,this%Sr,this%rpf)
end function

! **************************************************************************** !

function SF_VG_Set_m(this,m) result (error)
  class(SF_VG_type), intent(inout) :: this
  PetscReal, intent(in) :: m
  PetscInt :: error
  error = this%Configure(this%alpha,m,this%Sr,this%rpf)
end function

! **************************************************************************** !

pure subroutine SF_VG_Pc(this, Sl, Pc, dPc_dSl)
  class(SF_VG_type), intent(in) :: this
  PetscReal, intent(in)   :: Sl
  PetscReal, intent(out)  :: Pc, dPc_dSl

  if (Sl <= this%Sr) then               ! Unsaturated limit
    Pc = huge(Pc)
    dPc_dSl = -huge(Pc)
  else if (Sl > Sl_max) then            ! Saturated limit
    Pc= 0d0
    dPc_dSl = this%dPc_dSl_max
  else                                  ! Ordinary VG domain
    call this%Pc_inline(Sl,Pc,dPc_dSl)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_Sl(this, Pc, Sl, dSl_dPc)
  class(SF_VG_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sl, dSl_dPc

  if (Pc <= 0d0) then                   ! Saturated limit
    Sl = 1d0
    dSl_dPc = this%dSl_dPcmin
  else                                  ! Ordinary VG domain
    call this%Sl_inline(Pc,Sl,dSl_dPc)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_d2Sl_dPc2(this, Pc, d2Sl_dPc2)
  class(SF_VG_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2

  if (Pc <= 0d0) then                   ! Saturated limit
    d2Sl_dPc2 = this%d2Sl_dPc2min
  else                                  ! Ordinary VG domain
    call this%d2Sl_dPc2_inline(Pc,d2Sl_dPc2)
  end if
end subroutine

! **************************************************************************** !
! SF VG Extension Methods
! **************************************************************************** !

function SF_VG_extn_Set_alpha(this,alpha) result (error)
  class(SF_VG_extn_type), intent(inout) :: this
  PetscReal, intent(in) :: alpha
  PetscInt :: error
  PetscReal :: alpha_old

  alpha_old = this%alpha
  error = this%Configure(alpha,this%m,this%Sr,this%rpf)
  if (error == 0) then                  ! Update unsaturated extension
    if (this%Pcmax_designated) then
      error = this%Set_Pcmax(this%Pcmax)
    else
      error = this%Set_Sj(this%Sj)
    end if
    if (error /= 0) then                ! Restore previous state upon error
      error = error + this%Configure(alpha_old,this%m,this%Sr,this%rpf)
    end if
  end if
end function

! **************************************************************************** !

function SF_VG_extn_Set_m(this,m) result (error)
  class(SF_VG_extn_type), intent(inout) :: this
  PetscReal, intent(in) :: m
  PetscInt :: error
  PetscReal :: m_old

  m_old = this%m
  error = this%Configure(this%alpha,m,this%Sr,this%rpf)
  if (error == 0) then                   ! Update unsaturated extension
    if (this%Pcmax_designated) then
      error = this%Set_Pcmax(this%Pcmax)
    else
      error = this%Set_Sj(this%Sj)
    end if
    if (error /= 0) then                ! Restore previous state upon error
      error = error + this%Configure(this%alpha,m_old,this%Sr,this%rpf)
    end if
  end if
end function

! **************************************************************************** !
! VG Constant Extension Type
! **************************************************************************** !

function SF_VG_FCPC_ctor(alpha,m,Sr,rpf,Pcmax) result (new)
  class(SF_VG_cons_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Pcmax
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%Set_Pcmax(Pcmax) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function

! **************************************************************************** !

function SF_VG_cons_Set_Pcmax(this,Pcmax) result (error)
  class(SF_VG_cons_type), intent(inout) :: this
  PetscReal, intent(in) :: Pcmax
  PetscInt :: error

  if (Pcmax > 0d0) then
    error = 0
    this%Pcmax_designated = PETSC_TRUE
    this%Pcmax = Pcmax
    call this%Sl_inline(Pcmax,this%Sj,this%dSj_dPj)
    call this%d2Sl_dPc2_inline(Pcmax,this%d2Sj_dPj2)
  else
    error = 1
  end if
end function

! **************************************************************************** !

function SF_VG_FNOC_ctor(alpha,m,Sr,rpf,Sj) result (new)
  class(SF_VG_cons_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sj
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%Set_Sj(Sj) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function

! **************************************************************************** !

function SF_VG_cons_Set_Sj(this,Sj) result (error)
  class(SF_VG_cons_type), intent(inout) :: this
  PetscReal, intent(in) :: Sj
  PetscInt :: error
  PetscReal :: dPj_dSj

  if (Sj > this%Sr) then
    error = 0
    this%Pcmax_designated = PETSC_FALSE
    this%Sj = Sj
    call this%Pc_inline(Sj,this%Pcmax,dPj_dSj)
    this%dSj_dPj = 1d0 / dPj_dSj
    call this%d2Sl_dPc2_inline(this%Pcmax,this%d2Sj_dPj2)
  else
    error = 1
  end if
end function

! **************************************************************************** !

pure subroutine SF_VG_cons_Pc(this, Sl, Pc, dPc_dSl)
  class(SF_VG_cons_type), intent(in) :: this
  PetscReal, intent(in)  :: Sl
  PetscReal, intent(out) :: Pc, dPc_dSl

  if (Sl < this%Sj) then                ! Unsaturated limit
    Pc = this%Pcmax
    dPc_dSl = 0d0
  else if (Sl > Sl_max) then            ! Saturated limit
    Pc = 0d0
    dPc_dSl = this%dPc_dSl_max
  else                                  ! Ordinary VG domain
    call this%Pc_inline(Sl,Pc,dPc_dSl)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_cons_Sl(this, Pc, Sl, dSl_dPc)
  class(SF_VG_cons_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sl, dSl_dPc

  if (Pc <= 0d0) then                   ! Saturated limit
    Sl = 1d0
    dSl_dPc = this%dSl_dPcmin
  else if (Pc >= this%Pcmax) then       ! Unsaturated limit
    Sl = this%Sj
    dSl_dPc = this%dSj_dPj
  else                                  ! Ordinary VG domain
    call this%Sl_inline(Pc, Sl, dSl_dPc)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_cons_d2Sl_dPc2(this,Pc,d2Sl_dPc2)
  class(SF_VG_cons_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2

  if (Pc <= 0d0) then                   ! Saturated limit
    d2Sl_dPc2 = this%d2Sl_dPc2min
  else if (Pc >= this%Pcmax) then       ! Unsaturated limit
    d2Sl_dPc2 = this%d2Sj_dPj2
  else                                  ! Ordinary VG domain
    call this%d2Sl_dPc2_inline(Pc, d2Sl_dPc2)
  end if
end subroutine

! **************************************************************************** !
! Exponential Capillary Pressure Extension
! **************************************************************************** !

function SF_VG_ECPC_ctor(alpha,m,Sr,rpf,Pcmax) result (new)
  class(SF_VG_expn_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Pcmax
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%Set_Pcmax(Pcmax) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function

! **************************************************************************** !

function SF_VG_expn_Set_Pcmax(this,Pcmax) result (error)
  class(SF_VG_expn_type), intent(inout) :: this
  PetscReal, intent(in) :: Pcmax
  PetscInt :: error
  PetscReal :: Sa, Sb, Pe, dPj_dSj

  ! Solve the nonlinear system to satisfy C1 continuity using bisection method

  ! Set lower saturation limit (Sa) to be intercept (P(Sa) = Pcmax)
  call this%Sl_inline(Pcmax,Sa,dPj_dSj) ! dPj_dSj is inverted, but not used

  ! Set upper saturation limit (Sb) to be the the inflection point
  Sb = this%Sl_inflection()

  ! Confirm Pcmax is above minimum extrapolating from inflection point
  call this%Pc_inline(Sb,Pe,dPj_dSj)
  if (Pcmax >= Pe*exp(-dPj_dSj/Pe)) then
    ! Set Pcmax and begin iteration loop
    error = 0
    this%Pcmax = Pcmax
    this%Pcmax_designated = PETSC_TRUE

    do while (Sb-Sa > epsilon(this%Sj)) ! Tolerance interval epsilon
      this%Sj = (Sa+Sb)/2d0             ! Bisect bracket
      call this%Pc_inline(this%Sj,this%Pj,dPj_dSj)
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
end function

! **************************************************************************** !

function SF_VG_ENOC_ctor(alpha,m,Sr,rpf,Sj) result (new)
  class(SF_VG_expn_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sj
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%Set_Sj(Sj) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function

! **************************************************************************** !

function SF_VG_expn_Set_Sj(this,Sj) result (error)
  class(SF_VG_expn_type), intent(inout) :: this
  PetscReal, intent(in) :: Sj
  PetscInt :: error
  PetscReal :: dPj_dSj

  if (Sj > this%Sr) then
    error = 0
    this%Sj = Sj
    this%Pcmax_designated = PETSC_FALSE

    ! Find the slope at the piecewise junction point and extrapolate
    call this%Pc_inline(Sj,this%Pj,dPj_dSj)
    this%beta = dPj_dSj / this%Pj
    this%beta_rec = this%Pj / dPj_dSj
    this%Pcmax = this%Pj*exp(-this%beta*Sj)
    this%dPcmax_dSl = this%beta*this%Pcmax
    this%dSl_dPcmax = this%beta_rec/this%Pcmax
    this%d2Sl_dPc2max = -this%beta_rec/this%Pcmax**2
  else
    error = 1
  end if
end function

! **************************************************************************** !

pure subroutine SF_VG_expn_Pc(this, Sl, Pc, dPc_dSl)
  class(SF_VG_expn_type), intent(in) :: this
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
    call this%Pc_inline(Sl, Pc, dPc_dSl)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_expn_Sl(this, Pc, Sl, dSl_dPc)
  class(SF_VG_expn_type), intent(in) :: this
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
    call this%Sl_inline(Pc, Sl, dSl_dPc)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_expn_d2Sl_dPc2(this, Pc, d2Sl_dPc2)
  class(SF_VG_expn_type), intent(in) :: this
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
    call this%d2Sl_dPc2_inline(Pc, d2Sl_dPc2)
  end if
end subroutine

! **************************************************************************** !
! Linear Capillary Pressure Extension
! **************************************************************************** !

function SF_VG_LCPC_ctor(alpha,m,Sr,rpf,Pcmax) result (new)
  class(SF_VG_line_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Pcmax
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%Set_Pcmax(Pcmax) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function

! **************************************************************************** !

function SF_VG_line_Set_Pcmax(this,Pcmax) result (error)
  class(SF_VG_line_type), intent(inout) :: this
  PetscReal, intent(in) :: Pcmax
  PetscInt :: error
  PetscReal :: Sa, Sb, Pe, dPj_dSj

  ! Solve the nonlinear system to satisfy C1 continuity using bisection method

  ! Set lower saturation bracket Set where P(Sa) = Pcmax
  call this%Sl_inline(Pcmax,Sa,dPj_dSj)

  ! Set upper saturation limit (Sb) to be the the inflection point
  Sb = this%Sl_inflection()

  ! Confirm Pcmax is above minimum extrapolating from inflection point
  call this%Pc_inline(Sb,Pe,dPj_dSj)
  if (Pcmax > Pe - dPj_dSj*Sb) then
    error = 0
    this%Pcmax_designated = PETSC_TRUE
    this%Pcmax = Pcmax

    do while (Sb-Sa > epsilon(this%Sj)) ! Tolerance interval epsilon
      this%Sj = (Sa+Sb)/2d0             ! Bisect bracket
      call this%Pc_inline(this%Sj,this%Pj,this%dPj_dSj)

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
end function

! **************************************************************************** !

function SF_VG_LNOC_ctor(alpha,m,Sr,rpf,Sj) result (new)
  class(SF_VG_line_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sj
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%Set_Sj(Sj) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function

! **************************************************************************** !

function SF_VG_line_Set_Sj(this,Sj) result (error)
  class(SF_VG_line_type), intent(inout) :: this
  PetscReal, intent(in) :: Sj
  PetscInt :: error

  if (Sj > this%Sr .AND. Sj < 1d0) then
    error = 0
    this%Pcmax_designated = PETSC_FALSE
    this%Sj = Sj

    ! Find the value and slope at the piecewise junction point and extrapolate
    call this%Pc_inline(Sj,this%Pj,this%dPj_dSj)
    this%Pcmax = this%Pj - this%dPj_dSj * Sj
    this%dSj_dPj = 1d0 / this%dPj_dSj
  else
    error = 1
  end if
end function


! **************************************************************************** !

pure subroutine SF_VG_line_Pc(this, Sl, Pc, dPc_dSl)
  class(SF_VG_line_type), intent(in) :: this
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
    call this%Pc_inline(Sl, Pc, dPc_dSl)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_line_Sl(this, Pc, Sl, dSl_dPc)
  class(SF_VG_line_type), intent(in) :: this
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
    call this%Sl_inline(Pc, Sl, dSl_dPc)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_line_d2Sl_dPc2(this, Pc, d2Sl_dPc2)
  class(SF_VG_line_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2

  if (pc <= 0d0) then                   ! Saturated limit
    d2Sl_dPc2 = this%d2Sl_dPc2min
  else if (Pc >= this%Pj) then          ! Unsaturated domain
    d2Sl_dPc2 = 0d0
  else                                  ! Ordinary VG domain
    call this%d2Sl_dPc2_inline(Pc, d2Sl_dPc2)
  end if
end subroutine

! **************************************************************************** !
! Quadratic Capillary Pressure Extension
! **************************************************************************** !

function SF_VG_quad_ctor(alpha, m, Sr, rpf, Pcmax, Sj) result (new)
  class(SF_VG_quad_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Pcmax, Sj
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%Configure(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else if (new%Set_quad(Pcmax, Sj) /= 0) then
    deallocate(new)
    nullify(new)
  end if
end function

! **************************************************************************** !

function SF_VG_quad_Set_Pcmax(this, Pcmax) result (error)
  class(SF_VG_quad_type), intent(inout) :: this
  PetscReal, intent(in) :: Pcmax
  PetscInt :: error

  error = this%Set_quad(Pcmax, this%Sj)
end function

! **************************************************************************** !

function SF_VG_quad_Set_Sj(this, Sj) result (error)
  class(SF_VG_quad_type), intent(inout) :: this
  PetscReal, intent(in) :: Sj
  PetscInt :: error

  error = this%Set_quad(this%Pcmax, this%Sj)
end function
  
! **************************************************************************** !

function SF_VG_quad_Set_quad(this, Pcmax, Sj) result (error)
  class(SF_VG_quad_type), intent(inout) :: this
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
    call this%Pc_inline(Sj, Pj, dPj_dSj)
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
      Sl_qa = this%Sl_root(Pj)
      this%qa_branch = PETSC_FALSE
      Sl_cq = this%Sl_root(Pj)

      this%qa_branch = (abs(Sl_qa - Sj) < abs(Sl_cq - Sj))
    end if
  end if
end function

! **************************************************************************** !

pure subroutine SF_VG_quad_Pc(this, Sl, Pc, dPc_dSl)
  class(SF_VG_quad_type), intent(in) :: this
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
    call this%Pc_inline(Sl, Pc, dPc_dSl)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_quad_Sl(this, Pc, Sl, dSl_dPc)
  class(SF_VG_quad_type), intent(in) :: this
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
      Sl = this%Sl_root(Pc)
      dSl_dPc = 1d0 / (2d0*this%A*Sl + this%B)
    end if
  else                                  ! Ordinary VG domain
    call this%Sl_inline(Pc, Sl, dSl_dPc)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SF_VG_quad_d2Sl_dPc2(this, Pc, d2Sl_dPc2)
  class(SF_VG_quad_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2
  PetscReal :: Sl

  if (pc <= 0d0) then                   ! Saturated limit
    d2Sl_dPc2 = this%d2Sl_dPc2min
  else if (Pc >= this%Pj) then          ! Unsaturated domain
    if (Pc >= this%Pcmax) then          ! Unsaturated limit
      d2Sl_dPc2 = -1d0 / this%B**3
    else                                ! Quadratic extension
      Sl = this%Sl_root(Pc)
      d2Sl_dPc2 = -2d0*this%A / (2d0*this%A*Sl + this%B)**3
    end if
  else                                  ! Ordinary VG domain
    call this%d2Sl_dPc2_inline(Pc, d2Sl_dPc2)
  end if
end subroutine

! **************************************************************************** !

pure function SF_VG_quad_Sl_root(this,Pc) result (Sl)
  class(SF_VG_quad_type), intent(in) :: this
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
end function

! **************************************************************************** !
! van Genuchten Relative Permeability Methods
! **************************************************************************** !

pure function RPF_VG_Get_m(this) result (m)
  class(RPF_VG_type), intent(in) :: this
  PetscReal :: m
  m = this%m
end function

! **************************************************************************** !

subroutine RPF_VG_Init(this)
  class(RPF_VG_type) :: this
  ! This method is intentionally left blank.
end subroutine

! **************************************************************************** !

subroutine RPF_VG_RelativePermeability(this, liquid_saturation, &
                                       relative_permeability, dkr_sat, option)
  class(RPF_VG_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability, dkr_sat
  type(option_type), intent(inout) :: option

  call this%Kr(liquid_saturation, relative_permeability, dkr_sat)
end subroutine

! **************************************************************************** !
! Abstract Liquid van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPF_VG_liq_Configure(this, m, Sr) result (error)
  class(RPF_VG_liq_type), intent(inout) :: this
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

    call this%Kr_inline(Sl_max,Kr,this%dKr_dSl_max)
    this%dKr_dSl_max = (1d0 - Kr) / epsilon(Sl_max)

    this%Sl_min = Sr + epsilon(Sr)

    this%analytical_derivative_available = PETSC_TRUE
  end if
end function

! **************************************************************************** !

function RPF_VG_liq_Set_m(this, m) result (error)
  class(RPF_VG_liq_type), intent(inout) :: this
  PetscReal, intent(in) :: m
  PetscInt :: error
  error = this%Configure(m, this%Sr)
end function

! **************************************************************************** !

pure subroutine RPF_VG_liq_Kr(this, Sl, Kr, dKr_dSl)
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
    call this%Kr_inline(Sl, Kr, dKr_dSl)
  end if
end subroutine

! **************************************************************************** !
! Loop-invariant Liquid Mualem - van Genuchten Relative Permeability Methods
! **************************************************************************** !

function  RPF_MVG_liq_ctor(m, Sr) result (new)
  class(RPF_MVG_liq_type), pointer :: new
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
end function

! **************************************************************************** !

pure subroutine RPF_MVG_liq_Kr_inline(this, Sl, Kr, dKr_dSl)
  class(RPF_MVG_liq_type), intent(in) :: this
  PetscReal, intent(in) :: Sl
  PetscReal, intent(out) :: Kr, dKr_dSl
  PetscReal :: Se, Se_mrt, Se_mrt_comp, Se_mrt_comp_m, f

  Se = (Sl- this%Sr) * this%dSe_dSl
  Se_mrt  = Se**this%m_rec
  Se_mrt_comp = 1d0 - Se_mrt
  Se_mrt_comp_m = Se_mrt_comp**this%m
  f = 1d0 - Se_mrt_comp_m

  Kr = sqrt(Se)*f*f
  dKr_dSl = this%dSe_dSl * Kr / Se * &
           (0.5d0 + 2d0*Se_mrt*Se_mrt_comp_m/(f*Se_mrt_comp))
end subroutine

! **************************************************************************** !
! Loop-invariant Liquid Burdine- van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPF_BVG_liq_ctor(m, Sr) result (new)
  class(RPF_BVG_liq_type), pointer :: new
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
end function

! **************************************************************************** !

pure subroutine RPF_BVG_liq_Kr_inline(this, Sl, Kr, dKr_dSl)
  class(RPF_BVG_liq_type), intent(in) :: this
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
end subroutine

! **************************************************************************** !
! Abstract Gas van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPF_VG_gas_Configure(this, m, Sr, Sgr) result (error)
  class(RPF_VG_gas_type), intent(inout) :: this
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
end function

! **************************************************************************** !

function RPF_VG_gas_Set_m(this, m) result (error)
  class(RPF_VG_gas_type), intent(inout) :: this
  PetscReal, intent(in) :: m
  PetscInt :: error
  error = this%Configure(m, this%Sr, this%Sgr)
end function

! **************************************************************************** !

pure subroutine RPF_VG_gas_Kr(this, Sl, Kr, dKr_dSl)
  class(RPF_VG_gas_type), intent(in) :: this
  PetscReal, intent(in) :: Sl
  PetscReal, intent(out) :: Kr, dKr_dSl

  if (Sl > this%Sl_max) then            ! Saturated limit
    Kr = 0d0
    dKr_dSl = 0d0
  else if (Sl < this%Sr) then           ! Unsaturated limit
    Kr = 1d0
    dKr_dSl = 0d0                       ! Cusp at limit
  else                                  ! Ordinary domain
    call this%Kr_inline(Sl, Kr, dKr_dSl)
  end if
end subroutine

! **************************************************************************** !
! Loop-invariant Gas Mualem - van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPF_MVG_gas_ctor(m, Sr, Sgr) result (new)
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
end function

! **************************************************************************** !

pure subroutine RPF_MVG_gas_Kr_inline(this, Sl, Kr, dKr_dSl)
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
end subroutine

! **************************************************************************** !
! Gas Burdine - van Genuchten Relative Permeability Methods
! **************************************************************************** !

function RPF_BVG_gas_ctor(m, Sr, Sgr) result (new)
  class(RPF_BVG_gas_type), pointer :: new
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
end function

! **************************************************************************** !

pure subroutine RPF_BVG_gas_Kr_inline(this, Sl, Kr, dKr_dSl)
  class(RPF_BVG_gas_type), intent(in) :: this
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
end subroutine

! **************************************************************************** !

end module

