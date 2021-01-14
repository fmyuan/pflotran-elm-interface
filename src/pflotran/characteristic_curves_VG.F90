module Characteristic_Curves_VG_module
#include "petsc/finclude/petscsys.h"

use Characteristic_Curves_Base_module
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
! Child classes are used to override parent functions to avoid additional
! virtual function calls or branching statements in the inner loop.
! Loop invariant parameters are calculated upon construction for memoization.
!
! As the VG capillary pressure function has an infinite derivative approaching
! saturation, derivatives above Slmax are approximated using a backwards finite
! difference method using machine epsilon
!
! **************************************************************************** !

PetscReal, private, parameter :: Slmax = 1d0 - epsilon(Slmax) ! Saturated limit

! **************************************************************************** !
!
! VG Saturation Function Type Declarations
!
! sat_func_base_type        External definition in characteristic_curves_base
! |
! |-->sat_func_VG_type      Common parent, no extension
!     |                     Capillary pressure extensions below residual
!     |-->SF_VG_cons_type   Constant    extension
!     |-->SF_VG_expn_type   Exponential extension
!     |-->SF_VG_line_type   Linear      extension
!
! The objects are intended to be immutable and constructed using one of the
! provided constructors:
! 
! sat_func_VG_type => SF_VG_NEVG_ctor(alpha,m,Sr,rpf)
! SF_VG_cons_type => SF_VG_FCPC_ctor(alpha,m,Sr,rpf,Pcmax)
! SF_VG_cons_type => SF_VG_FNOC_ctor(alpha,m,Sr,rpf,Sj)
! SF_VG_expn_type => SF_VG_ECPC_ctop(alpha,m,Sr,rpf,Pcmax)
! SF_VG_expn_type => SF_VG_ENOC_ctop(alpha,m,Sr,rpf,Sj)
! SF_VG_line_type => SF_VG_LCPC_ctor(alpha,m,Sr,rpf,Pcmax)
! SF_VG_line_type => SF_VG_LNOC_ctor(alpha,m,Sr,rpf,Sj)
!
! **************************************************************************** !

! Parent Saturation Function VG Type
type, public, extends(sat_func_base_type) :: sat_func_VG_type
  private                          ! Immutable loop-invariant parameters
    PetscReal :: alpha, a_rec      ! Alpha and alpha reciprocal
    PetscReal :: m, m_nrec, m_a2   ! M, negative reciprocal M, M add 2
    PetscReal :: n, n_rec, n_m1    ! N, reciprocal N, and N minus 1
    PetscReal :: mn_a1             ! MN product add 1
    PetscReal :: Sl_span, dSe_dSl  ! Effective saturation span and reciprocal
    PetscReal :: dSe_mndSl         ! Coefficent of VG derivative
    PetscReal :: dPc_dSlmax        ! Finite difference derivative at saturation
    PetscReal :: dSl_dPcmin
    PetscReal :: d2Sl_dPc2min
contains
  ! Public compatibility methods
  procedure, public :: Verify            => SF_VG_Verify
  ! Common parameter data validation and initialization 
  procedure, private :: SF_VG_init
  ! Common member methods for ordinary VG domain
  procedure, private :: SF_VG_Pc_inline
  procedure, private :: SF_VG_Sl_inline
  procedure, private :: SF_VG_d2Sl_dPc2_inline
  ! Public accessor methods
  procedure, public :: get_alpha         => SF_VG_get_alpha
  procedure, public :: get_m             => SF_VG_get_m
  ! Van Genuchten functions without unsaturated extension
  procedure, public :: CapillaryPressure => SF_VG_Pc
  procedure, public :: Saturation        => SF_VG_Sl
  procedure, public :: D2SatDP2          => SF_VG_d2Sl_dPc2
end type

! Child VG Constant Capillary Pressure Type
type, public, extends(sat_func_VG_type) :: SF_VG_cons_type
  private
    PetscReal :: Sj        ! Piecewise junction point
    PetscReal :: dSj_dPj   ! Non-zero derivative at cusp
    PetscReal :: d2Sj_dPj2 ! Non-zero 2nd derivative at cusp
contains
  procedure, public :: CapillaryPressure => SF_VG_cons_Pc
  procedure, public :: Saturation        => SF_VG_cons_Sl
  procedure, public :: D2SatDP2          => SF_VG_cons_d2Sl_dPc2
end type

! Child VG Exponential Capillary Pressure Type
type, public, extends(sat_func_VG_type) :: SF_VG_expn_type
  private
    PetscReal :: Sj, Pj           ! Piecewise junction point
    PetscReal :: beta, beta_rec   ! Exponential coefficient and reciprocal
    PetscReal :: dPcmax_dSl       ! Unsaturated limits
    PetscReal :: dSl_dPcmax
    PetscReal :: d2Sl_dPc2max
contains
  procedure, public :: CapillaryPressure => SF_VG_expn_Pc
  procedure, public :: Saturation        => SF_VG_expn_Sl
  procedure, public :: D2SatDP2          => SF_VG_expn_d2Sl_dPc2
end type

! Child VG Linear Capillary Pressure Type
type, public, extends(sat_func_VG_type) :: SF_VG_line_type
  private
    PetscReal :: Sj, Pj           ! Piecewise junction point
    PetscReal :: dPj_dSj, dSj_dPj ! Linear coefficients
contains
  procedure, public :: CapillaryPressure => SF_VG_line_Pc
  procedure, public :: Saturation        => SF_VG_line_Sl
  procedure, public :: D2SatDP2          => SF_VG_line_d2Sl_dPc2
end type

! **************************************************************************** !
!
! VG Relative Permeability Function Type Declarations
!
! rel_perm_func_base_type   External definition in characteristic_curves_base
! |
! |-->RPF_VG_type           Common van Genuchten parent
!     |
!     |-->RPF_MVG_liq_type  Mualem  - van Genuchten liquid 
!     |-->RPF_MVG_gas_type  Mualem  - van Genuchten gas
!     |-->RPF_BVG_liq_type  Burdine - van Genuchten liquid 
!     |-->RPF_BVG_gas_type  Burdine - van Genuchten gas
!
! **************************************************************************** !

! Parent VG Relative Permeability Function Type
type, public, extends(rel_perm_func_base_type) :: RPF_VG_type
 ! private TODO make private when data encapsulation is complete
    PetscReal :: m, m_rec   ! Van Genuchten m exponent and reciprocal
    PetscReal :: dSe_dSl    ! Effective saturation span reciprocal
    PetscReal :: dKr_dSlmax ! Finite difference derivative at saturation
  contains
    procedure, public :: get_m             => RPF_VG_get_m
end type

! Child Mualem VG Relative Permeability Function Types
type, public, extends(RPF_VG_type) :: RPF_MVG_liq_type
contains
  procedure, private :: RPF_MVG_liq_inline
  procedure, public :: Verify               => RPF_MVG_liq_verify
  procedure, public :: RelativePermeability => RPF_MVG_liq_relperm
end type

type, public, extends(RPF_VG_type) :: RPF_MVG_gas_type
  private
    PetscReal :: mx2        ! m exponent times 2
    PetscReal :: Sgr        ! Gas residual saturation
    PetscReal :: Slmax      ! Saturated limit
    PetscReal :: Slmin      ! Unsaturated limit
    PetscReal :: dKr_dSlmin ! Finite difference derivative at residual
contains
  procedure, private :: RPF_MVG_gas_inline
  procedure, public :: Verify               => RPF_MVG_gas_verify
  procedure, public :: RelativePermeability => RPF_MVG_gas_relperm
end type

! Child Burdine VG Relative Permeability Function Types
type, public, extends(RPF_VG_type) :: RPF_BVG_liq_type
contains
  procedure, private :: RPF_BVG_liq_inline
  procedure, public :: Verify               => RPF_BVG_liq_verify
  procedure, public :: RelativePermeability => RPF_BVG_liq_relperm
end type

type, public, extends(RPF_VG_type) :: RPF_BVG_gas_type
  private 
    PetscReal :: Sgr        ! Gas residual saturation
    PetscReal :: Slmax      ! Saturated limit
    PetscReal :: Slmin      ! Unsaturated limit
    PetscReal :: dKr_dSlmin ! Finite difference derivative at residual
contains
  procedure, private :: RPF_BVG_gas_inline
  procedure, public :: Verify               => RPF_BVG_gas_verify
  procedure, public :: RelativePermeability => RPF_BVG_gas_relperm
end type

! **************************************************************************** !
! Public constuctors and procedures
! **************************************************************************** !

! Public VG Saturation Function constructors
  public :: SF_VG_Create, &    ! For backwards compatibility
            SF_VG_NEVG_ctor, & ! No capped capillary pressure
            SF_VG_FCPC_ctor, & ! Flat, specified maximum
            SF_VG_FNOC_ctor, & ! Flat, specified junction
            SF_VG_ECPC_ctor, & ! Exponential, specified maximum
            SF_VG_ENOC_ctor, & ! Exponential, specified junction
            SF_VG_LCPC_ctor, & ! Linear, specified maximum
            SF_VG_LNOC_ctor    ! Linear, specified junction

! Public VG Relative Permeability constructors
  public :: RPFMualemVGLiqCreate, &
            RPF_MVG_Liq_ctor, &
            RPFMualemVGGasCreate, &
            RPF_MVG_Gas_ctor, &
            RPFBurdineVGLiqCreate, &
            RPF_BVG_Liq_ctor, &
            RPFBurdineVGGasCreate, &
            RPF_BVG_Gas_ctor

! Public VG Relative Permeability method accessed by child class
  public :: RPF_MVG_liq_relperm

contains
! **************************************************************************** !
! Common VG Saturation Function Methods
! **************************************************************************** !

function SF_VG_Create() result (new)
  ! Creates a dummy van Genutchten saturation function object
  ! The CharacteristicCurvesRead function is passed an object with a generic type
  ! Once the block is parsed, a new object is created that replaces this object.
  ! If object creation by the parser is delayed, this method will not be needed
  implicit none
  class(sat_func_VG_type), pointer :: new
  allocate(new)
end function

! **************************************************************************** !

subroutine SF_VG_Verify(this,name,option)
  implicit none
  class(sat_func_VG_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,VAN_GENUCHTEN'
  endif
end subroutine

! **************************************************************************** !

function SF_VG_init(this,alpha,m,Sr,rpf) result (error)
  ! Initalize loop-invariant parameters common to all extensions
  class(sat_func_VG_type) :: this
  PetscReal, intent(in) :: alpha,m,Sr
  PetscInt, intent(in) :: rpf
  PetscInt :: error
  PetscReal :: Pc1, Pc2, dPc_dSl

  ! Can eliminate the pointers if smoothing implementation is replaced
  nullify(this%sat_poly)
  nullify(this%pres_poly) 
  this%analytical_derivative_available = .TRUE.
  this%calc_int_tension = .FALSE.
 
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
    ! Pc0 = 0 = Pc(1); Pc1 = Pc(1-eps); Pc2 = Pc(1-2*eps); 
    call this%SF_VG_Pc_inline(Slmax,Pc1,dPc_dSl)
    call this%SF_VG_Pc_inline(Slmax-epsilon(Slmax),Pc2,dPc_dSl)
    this%dPc_dSlmax = -Pc1 / epsilon(Slmax)
    this%dSl_dPcmin = -epsilon(Slmax) / Pc1
    this%d2Sl_dPc2min = epsilon(Slmax)*(2d0-Pc2/Pc1)/Pc1**2
  end if
end function

! **************************************************************************** !

pure subroutine SF_VG_Pc_inline(this,Sl,Pc,dPc_dSl)
  class(sat_func_VG_type), intent(in) :: this
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

pure subroutine SF_VG_Sl_inline(this,Pc,Sl,dSl_dPc)
  class(sat_func_VG_type), intent(in) :: this
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

pure subroutine SF_VG_d2Sl_dPc2_inline(this,Pc,d2Sl_dPc2)
  class(sat_func_VG_type), intent(in) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2Sl_dPc2
  PetscReal :: aPc_n, Se_mrtrec, d2Se_mndPc2

  aPc_n = (this%alpha*Pc)**this%n
  Se_mrtrec = aPc_n + 1d0

  ! Note, the exponentiation of Se_mrtrec could be reduced if Se is provived
  ! Se_mrtrec**(m+2) == Se_mrtrec**2/Se
  ! As it is not, the terms are combined

  d2Se_mndPc2 =  aPc_n*(this%mn_a1*aPc_n-this%n_m1)/(Se_mrtrec**this%m_a2*Pc**2)
  d2Sl_dPc2 = d2Se_mndPc2 / this%dSe_mndSl
end subroutine

! **************************************************************************** !

pure function SF_VG_get_m(this) result (m)
  implicit none
  class(sat_func_VG_type), intent(in) :: this
  PetscReal :: m
  m = this%m
end function

! **************************************************************************** !

pure function SF_VG_get_alpha(this) result (alpha)
  class(sat_func_VG_type), intent(in) :: this
  PetscReal :: alpha
  alpha = this%alpha
end function

! **************************************************************************** !
! No-extension VG Saturation Function Methods
! **************************************************************************** !
function SF_VG_NEVG_ctor(alpha,m,Sr,rpf) result (new)
  class(sat_func_VG_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr
  PetscInt,  intent(in) :: rpf

  allocate(new)

  if (new%SF_VG_init(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  end if
  new%Pcmax = huge(new%Pcmax)
end function
! **************************************************************************** !

subroutine SF_VG_Pc(this, liquid_saturation, &
                         capillary_pressure, dPc_dSatl, option)
  implicit none
  class(sat_func_VG_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation <= this%Sr) then
    ! Unsaturated limit
    capillary_pressure = huge(capillary_pressure)
    dPc_dSatl = -huge(capillary_pressure)
  else if (liquid_saturation > Slmax) then
    ! Saturated limit
    capillary_pressure = 0d0
    dPc_dSatl = this%dPc_dSlmax
  else
    ! Ordinary VG domain
    call this%SF_VG_Pc_inline(liquid_saturation,capillary_pressure,dPc_dSatl)
  end if
end subroutine

! **************************************************************************** !

subroutine SF_VG_Sl(this, capillary_pressure, &
                         liquid_saturation, dsat_dpres, option)
  implicit none
  class(sat_func_VG_type) :: this
  PetscReal, intent(in)  :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  if (capillary_pressure <= 0d0) then
    ! Saturated limit
    liquid_saturation = 1d0
    dsat_dpres = this%dSl_dPcmin
  else
    ! Ordinary VG domain
    call this%SF_VG_Sl_inline(capillary_pressure,liquid_saturation,dsat_dpres)
  end if
end subroutine


! **************************************************************************** !

subroutine SF_VG_d2Sl_dPc2(this,Pc,d2s_dp2,option)
  implicit none
  class(sat_func_VG_type) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option

  if (pc <= 0d0) then
    ! Saturated limit
    d2s_dp2 = this%d2Sl_dPc2min
  else
    ! Ordinary VG domain
    call this%SF_VG_d2Sl_dPc2_inline(Pc,d2s_dp2)
  end if
end subroutine

! **************************************************************************** !
! Child Constant Extension Type
! **************************************************************************** !
function SF_VG_FCPC_ctor(alpha,m,Sr,rpf,Pcmax) result (new)
  class(SF_VG_cons_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Pcmax
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%SF_VG_init(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else
    if (Pcmax <= 0d0) then
      deallocate(new)
      nullify(new)
    else
      ! Find the piecewise junction point saturation (Pc(Sj) = Pcmax)
      new%Pcmax = Pcmax
      call new%SF_VG_Sl_inline(Pcmax,new%Sj,new%dSj_dPj)
      call new%SF_VG_d2Sl_dPc2_inline(Pcmax,new%d2Sj_dPj2)
    end if
  end if
end function

! **************************************************************************** !

function SF_VG_FNOC_ctor(alpha,m,Sr,rpf,Sj) result (new)
  class(SF_VG_cons_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sj
  PetscInt,  intent(in) :: rpf
  PetscReal :: dPj_dSj

  allocate(new)
  if (new%SF_VG_init(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else
    if (Sj <= Sr) then
      deallocate(new)
      nullify(new)
    else
      ! Find the capillary pressure at piecewise junction point Sj
      new%Sj = Sj
      call new%SF_VG_Pc_inline(Sj,new%Pcmax,dPj_dSj)
      new%dSj_dPj = 1d0 / dPj_dSj
      call new%SF_VG_d2Sl_dPc2_inline(new%Pcmax,new%d2Sj_dPj2)
    end if
  end if
end function

! **************************************************************************** !

subroutine SF_VG_cons_Pc(this, liquid_saturation, &
                         capillary_pressure, dPc_dSatl, option)
  implicit none
  class(SF_VG_cons_type) :: this
  PetscReal, intent(in)  :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation < this%Sj) then
    ! Unsaturated limit
    capillary_pressure = this%Pcmax
    dPc_dSatl = 0d0
  else if (liquid_saturation > Slmax) then
    ! Saturated limit
    capillary_pressure = 0d0
    dPc_dSatl = this%dPc_dSlmax
  else
    ! Ordinary VG domain
    call this%SF_VG_Pc_inline(liquid_saturation,capillary_pressure,dPc_dSatl)
  end if

end subroutine

! **************************************************************************** !

subroutine SF_VG_cons_Sl(this, capillary_pressure, &
                         liquid_saturation, dsat_dpres, option)
  implicit none
  class(SF_VG_cons_type) :: this
  PetscReal, intent(in)  :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  if (capillary_pressure <= 0d0) then
    ! Saturated limit
    liquid_saturation = 1d0
    dsat_dpres = this%dSl_dPcmin
  else if (capillary_pressure >= this%Pcmax) then
    ! Unsaturated limit
    liquid_saturation = this%Sj
    dsat_dpres = this%dSj_dPj
  else
    ! Ordinary VG domain
    call this%SF_VG_Sl_inline(capillary_pressure,liquid_saturation,dsat_dpres)
  end if

end subroutine

! **************************************************************************** !

subroutine SF_VG_cons_d2Sl_dPc2(this,Pc,d2s_dp2,option)
  implicit none
  class(SF_VG_cons_type) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option

  if (pc <= 0d0) then
    ! Saturated limit
    d2s_dp2 = this%d2Sl_dPc2min
  else if (Pc >= this%Pcmax) then
    ! Unsaturated limit
    d2s_dp2 = this%d2Sj_dPj2
  else
    ! Ordinary VG domain
    call this%SF_VG_d2Sl_dPc2_inline(Pc,d2s_dp2)
  end if
end subroutine

! **************************************************************************** !
! Child Exponential Capillary Pressure Extension
! **************************************************************************** !
function SF_VG_ECPC_ctor(alpha,m,Sr,rpf,Pcmax) result (new)
  class(SF_VG_expn_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Pcmax
  PetscInt,  intent(in) :: rpf
  PetscReal :: Se, Se_mrtrec, Sa, Sb, Pe, dPj_dSj

  allocate(new)
  if (new%SF_VG_init(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else
    ! Solve the nonlinear system to ensure function is C1 continuous
    ! Set upper limit of the bracket to be the the inflection point
    Se_mrtrec = (1d0+new%m)/(new%n_rec+new%m)
    Se = Se_mrtrec**(-new%m)
    Pe = new%a_rec * (Se_mrtrec-1d0)**new%n_rec
    Sb = new%Sr + Se*new%Sl_span
    ! Confirm Pcmax is greater than P(Sb)
    if (Pcmax <= Pe) then
      deallocate(new)
      nullify(new)
    else
      new%Pcmax = Pcmax
      ! Set lower saturation bracket set where P(Sa) = Pcmax
      Se = (1d0 + (new%alpha*new%Pcmax)**new%n)**(-new%m)
      Sa = new%Sr + Se*new%Sl_span
      do while (Sb-Sa > epsilon(Sa)) ! Use intrinsic epsilon() as Sl is O(1)
        ! Bisect bracket bound between (Sa,Sb):
        new%Sj = (Sa+Sb)/2d0
        call new%SF_VG_Pc_inline(new%Sj,new%Pj,dPj_dSj)
        new%beta = dPj_dSj / new%Pj
        ! Residual error = Pcmax*exp(beta*Sj) - Pf(Sj)
        Pe = new%Pcmax*exp(new%beta*new%Sj) - new%Pj
        if (Pe < 0d0) then ! Error on side a is always negative for VG
          Sa = new%Sj
        else
          Sb = new%Sj
        end if
      end do
      new%dPcmax_dSl = new%beta*new%Pcmax
      new%beta_rec = 1d0 / new%beta
      new%dSl_dPcmax = new%beta_rec/new%Pcmax
      new%d2Sl_dPc2max = -new%beta_rec/new%Pcmax**2
    end if
  end if
end function

! **************************************************************************** !

function SF_VG_ENOC_ctor(alpha,m,Sr,rpf,Sj) result (new)
  class(SF_VG_expn_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sj
  PetscInt,  intent(in) :: rpf
  PetscReal :: dPj_dSj

  allocate(new)
  if (new%SF_VG_init(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else
    if (Sj <= Sr) then
      deallocate(new)
      nullify(new)
    else
      ! Find the slope at the piecewise junction point and extrapolate
      new%Sj = Sj
      call new%SF_VG_Pc_inline(Sj,new%Pj,dPj_dSj)
      new%beta = dPj_dSj / new%Pj
      new%beta_rec = new%Pj / dPj_dSj
      new%Pcmax = new%Pj*exp(-new%beta*Sj)
      new%dPcmax_dSl = new%beta*new%Pcmax
      new%dSl_dPcmax = new%beta_rec/new%Pcmax
      new%d2Sl_dPc2max = -new%beta_rec/new%Pcmax**2
    end if
  end if
end function

! **************************************************************************** !

subroutine SF_VG_expn_Pc(this, liquid_saturation, &
                         capillary_pressure, dPc_dSatl, option)
  implicit none
  class(SF_VG_expn_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation < this%Sj) then
    ! Exponential domain
    if (liquid_saturation < 0d0) then
      ! Unsaturated limit
      capillary_pressure = this%Pcmax
      dPc_dSatl = this%dPcmax_dSl
    else
      ! Exponential extension
      capillary_pressure = this%Pcmax*exp(this%beta*liquid_saturation)
      dPc_dSatl = this%beta*capillary_pressure
    end if
  else if (liquid_saturation >= 1d0) then
    ! Saturated limit
    capillary_pressure = 0d0
    dPc_dSatl = this%dPc_dSlmax
  else
    ! Ordinary VG domain
    call this%SF_VG_Pc_inline(liquid_saturation,capillary_pressure,dPc_dSatl)
  end if

end subroutine

! **************************************************************************** !

subroutine SF_VG_expn_Sl(this, capillary_pressure, &
                         liquid_saturation, dsat_dpres, option)
  implicit none
  class(SF_VG_expn_type) :: this
  PetscReal, intent(in)  :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  if (capillary_pressure <= 0d0) then
    ! Saturated limit
    liquid_saturation = 1d0
    dsat_dpres = this%dSl_dPcmin
  else if (capillary_pressure >= this%Pj) then
    ! Exponential domain
    if (capillary_pressure >= this%Pcmax) then
      ! Unsaturated limit
      liquid_saturation = 0d0
      dsat_dpres = this%dSl_dPcmax
    else
      ! Exponential extension
      liquid_saturation = this%beta_rec*log(capillary_pressure/this%Pcmax)
      dsat_dpres = this%beta_rec/capillary_pressure
    end if
  else
    ! Ordinary VG domain
    call this%SF_VG_Sl_inline(capillary_pressure,liquid_saturation,dsat_dpres)
  end if

end subroutine

! **************************************************************************** !

subroutine SF_VG_expn_d2Sl_dPc2(this,Pc,d2s_dp2,option)
  implicit none
  class(SF_VG_expn_type) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option

  if (pc <= 0d0) then
    ! Saturated limit
    d2s_dp2 = this%d2Sl_dPc2min
  else if (Pc >= this%Pj) then
    ! Exponential domain
    if (Pc >= this%Pcmax) then
      ! Unsaturated limit
      d2s_dp2 = this%d2Sl_dPc2max
    else
      ! Exponential extension
      d2s_dp2 = -this%beta_rec/Pc**2
    end if
  else
    ! Ordinary VG domain
    call this%SF_VG_d2Sl_dPc2_inline(Pc,d2s_dp2)
  end if
end subroutine

! **************************************************************************** !
! Child Linear Capillary Pressure Extension
! **************************************************************************** !
function SF_VG_LNOC_ctor(alpha,m,Sr,rpf,Sj) result (new)
  class(SF_VG_line_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Sj
  PetscInt,  intent(in) :: rpf

  allocate(new)
  if (new%SF_VG_init(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else
    if (Sj <= Sr .OR. Sj >= 1d0) then
      deallocate(new)
      nullify(new)
    else
      ! Find the slope at the piecewise junction point and extrapolate
      new%Sj = Sj
      call new%SF_VG_Pc_inline(Sj,new%Pj,new%dPj_dSj)
      new%Pcmax = new%Pj - new%dPj_dSj * new%Sj
      new%dSj_dPj = 1d0 / new%dPj_dSj
    end if
  end if
end function

! **************************************************************************** !

function SF_VG_LCPC_ctor(alpha,m,Sr,rpf,Pcmax) result (new)
  class(SF_VG_line_type), pointer :: new
  PetscReal, intent(in) :: alpha, m, Sr, Pcmax
  PetscInt,  intent(in) :: rpf
  PetscReal :: Se, Se_mrtrec, Sa, Sb, Pe

  allocate(new)
  if (new%SF_VG_init(alpha,m,Sr,rpf) /= 0) then
    deallocate(new)
    nullify(new)
  else
    ! Solve the nonlinear system to ensure function is C1 continuous
    ! Set upper limit of the bracket to be the the inflection point
    Se_mrtrec = (1d0+new%m)/(new%n_rec+new%m)
    Se = Se_mrtrec**(-new%m)
    Pe = new%a_rec * (Se_mrtrec-1d0)**new%n_rec
    Sb = new%Sr + Se*new%Sl_span
    ! Confirm Pcmax is greater than P(Sb)
    if (Pcmax <= Pe) then
      deallocate(new)
      nullify(new)
    else
      new%Pcmax = Pcmax
      ! Set lower saturation bracket set where P(Sa) = Pcmax
      Se = (1d0 + (new%alpha*new%Pcmax)**new%n)**(-new%m)
      Sa = new%Sr + Se*new%Sl_span
      do while (Sb-Sa > epsilon(Sa)) ! Use intrinsic epsilon() as Sl is O(1)
        ! Bisect bracket bound between (Sa,Sb):
        new%Sj = (Sa+Sb)/2d0
        call new%SF_VG_Pc_inline(new%Sj,new%Pj,new%dPj_dSj)
        ! Residual error = Pcmax + beta * Sj - Pf(Sj)
        Pe = new%Pcmax + new%dPj_dSj*new%Sj - new%Pj
        if (Pe < 0d0) then ! Error on side a is always negative for VG
          Sa = new%Sj
        else
          Sb = new%Sj
        end if
      end do
      new%dSj_dPj = 1d0 / new%dPj_dSj
    end if
  end if
end function

! **************************************************************************** !

subroutine SF_VG_line_Pc(this, liquid_saturation, &
                         capillary_pressure, dPc_dSatl, option)
  implicit none
  class(SF_VG_line_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure, dPc_dSatl
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation < this%Sj) then
    ! Linear domain
    dPc_dSatl = this%dPj_dSj
    if (liquid_saturation < 0d0) then
      ! Unsaturated limit
      capillary_pressure = this%Pcmax
    else
      ! Linear extension
      capillary_pressure = this%Pcmax + this%dPj_dSj*liquid_saturation
    end if
  else if (liquid_saturation > Slmax) then
    ! Saturated limit
    capillary_pressure = 0d0
    dPc_dSatl = this%dPc_dSlmax
  else
    ! Ordinary VG domain
    call this%SF_VG_Pc_inline(liquid_saturation,capillary_pressure,dPc_dSatl)
  end if

end subroutine

! **************************************************************************** !

subroutine SF_VG_line_Sl(this, capillary_pressure, &
                         liquid_saturation, dsat_dpres, option)
  implicit none
  class(SF_VG_line_type) :: this
  PetscReal, intent(in)  :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  if (capillary_pressure <= 0d0) then
    ! Saturated limit
    liquid_saturation = 1d0
    dsat_dpres = this%dSl_dPcmin
  else if (capillary_pressure >= this%Pj) then
    ! Linear domain
    dsat_dpres = this%dSj_dPj
    if (capillary_pressure >= this%Pcmax) then
      ! Unsaturated limit
      liquid_saturation = 0d0
    else
      ! Linear extensions
      liquid_saturation = (capillary_pressure - this%Pcmax) * this%dSj_dPj
    end if
  else
    ! Ordinary VG domain
    call this%SF_VG_Sl_inline(capillary_pressure,liquid_saturation,dsat_dpres)
  end if

end subroutine

! **************************************************************************** !

subroutine SF_VG_line_d2Sl_dPc2(this,Pc,d2s_dp2,option)
  implicit none
  class(SF_VG_line_type) :: this
  PetscReal, intent(in) :: Pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option

  if (pc <= 0d0) then
    ! Saturated limit
    d2s_dp2 = this%d2Sl_dPc2min
  else if (Pc >= this%Pj) then
    ! Linear domain
    d2s_dp2 = 0d0
  else
    ! Ordinary VG domain
    call this%SF_VG_d2Sl_dPc2_inline(Pc,d2s_dp2)
  end if
end subroutine

! **************************************************************************** !
! Common van Genuchten Relative Permeability Function Methods
! **************************************************************************** !

pure function RPF_VG_get_m(this) result (m)
  implicit none
  class(rpf_vg_type), intent(in) :: this
  PetscReal :: m
  m = this%m
end function

! **************************************************************************** !
! Mualem - van Genuchten Liquid Relative Permeability Function Methods
! **************************************************************************** !

function RPFMualemVGLiqCreate() result (new)
  ! Creates a dummy Mualem - van Genutchten relative permeability object
  ! The CharacteristicCurvesRead function is passed an object with a generic type
  ! Once the block is parsed, a new object is created that replaces this object.
  ! If object creation by the parser is delayed, this method will not be needed
  implicit none
  class(rpf_MVG_liq_type), pointer :: new
  allocate(new)
end function

! **************************************************************************** !

subroutine RPF_MVG_liq_verify(this,name,option)
  implicit none
  class(rpf_MVG_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MVG_LIQ'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call PrintErrMsg(option)
  end if
end subroutine

! **************************************************************************** !
  
function  RPF_MVG_liq_ctor(m, Sr) result (new)
  class(rpf_MVG_liq_type), pointer :: new
  PetscReal, intent(in) :: m, Sr
  PetscInt :: error
  PetscReal :: Kr

  error = 0
  if (m <= 0d0 .OR. m >= 1d0) error = error + 1
  if (Sr < 0d0 .OR. Sr >= 1d0) error = error + 2

  if (error == 0) then
    allocate(new)
    new%analytical_derivative_available = .TRUE.
    nullify(new%poly)

    new%m = m
    new%m_rec = 1d0 / m
    new%Sr = Sr
    new%dSe_dSl = 1d0 / (1d0 - Sr)

    call new%RPF_MVG_liq_inline(Slmax, Kr, new%dKr_dSlmax)
    new%dKr_dSlmax = (1d0 - Kr) / epsilon(Slmax)
  else
    nullify(new)
  end if
end function

! **************************************************************************** !

pure subroutine RPF_MVG_liq_inline(this, Sl, Kr, dKr_dSl)
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
  dKr_dSl = this%dSe_dSl * Kr / Se * &
           (0.5d0 + 2d0*Se_mrt*Se_mrt_comp_m/(f*Se_mrt_comp))
end subroutine

! **************************************************************************** !

subroutine RPF_MVG_liq_relperm(this, liquid_saturation, relative_permeability, &
                               dkr_sat, option)
  implicit none
  class(rpf_MVG_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability, dkr_sat
  type(option_type), intent(inout) :: option

  if (liquid_saturation > Slmax) then           ! Saturated limit
    relative_permeability = 1d0
    dkr_sat = this%dKr_dSlmax
  else if (liquid_saturation <= this%Sr) then ! Unsaturated limit
    relative_permeability = 0d0
    dkr_sat = 0d0
  else                                           ! Ordinary MVG
    call this%RPF_MVG_liq_inline(liquid_saturation, relative_permeability, &
                                 dkr_sat)
  end if
end subroutine

! **************************************************************************** !
! Mualem - van Genuchten Gas Relative Permeability Function Methods
! **************************************************************************** !

function RPFMualemVGGasCreate() result (new)
  ! The CharacteristicCurvesRead function is passed an object with a generic type
  ! Once the block is parsed, a new object is created that replaces this object.
  ! If object creation by the parser is delayed, this method will not be needed
  implicit none
  class(rpf_MVG_gas_type), pointer :: new
  allocate(new)
end function

! **************************************************************************** !

subroutine RPF_MVG_gas_verify(this,name,option)
  implicit none
  class(rpf_MVG_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MVG_GAS'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif 
end subroutine

! **************************************************************************** !

function  RPF_MVG_gas_ctor(m, Slr, Sgr) result (new)
  class(rpf_MVG_gas_type), pointer :: new
  PetscReal, intent(in) :: m, Slr, Sgr
  PetscInt :: error
  PetscReal :: Kr

  error = 0
  if (m <= 0d0 .OR. m >= 1d0) error = error + 1
  if (Slr < 0d0 .OR. Sgr < 0d0 .OR. Slr + Sgr >= 1d0) error = error + 2

  if (error == 0) then
    allocate(new)
    new%analytical_derivative_available = .TRUE.
    nullify(new%poly)

    new%m = m
    new%m_rec = 1d0 / m
    new%mx2 = 2d0 * m

    new%Sr = Slr
    new%Sgr = Sgr
    new%dSe_dSl = 1d0 / (1d0 - Slr - Sgr)

    new%Slmax = 1d0 - Sgr - epsilon(Sgr)
    call new%RPF_MVG_gas_inline(new%Slmax, Kr, new%dKr_dSlmax)
    new%dKr_dSlmax = - Kr / epsilon(Sgr)

    new%Slmin = new%Sr + epsilon(Slr)
    new%dKr_dSlmin = -0.5d0 * new%dSe_dSl
  else
    nullify(new)
  end if
end function

! **************************************************************************** !

pure subroutine RPF_MVG_gas_inline(this, Sl, Kr, dKr_dSl)
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
end subroutine

! **************************************************************************** !

subroutine RPF_MVG_gas_relperm(this,liquid_saturation, relative_permeability, &
                               dkr_sat, option)
  implicit none
  class(rpf_MVG_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability, dkr_sat
  type(option_type), intent(inout) :: option

  if (liquid_saturation > this%Slmax) then
    relative_permeability = 0d0
    dkr_sat = this%dKr_dSlmax
  else if (liquid_saturation < this%Slmin) then
    relative_permeability = 1d0
    dkr_sat = this%dKr_dSlmin
  else
    ! Ordinary MVG
    call this%RPF_MVG_gas_inline(liquid_saturation, relative_permeability, &
                                 dkr_sat)
  end if
end subroutine

! **************************************************************************** !
! Burdine - van Genuchten Liquid Relative Permeability Function Methods
! **************************************************************************** !

function RPFBurdineVGLiqCreate() result (new)
  ! Creates the van Genutchten Mualem relative permeability function object
  implicit none
  class(RPF_BVG_liq_type), pointer :: new
  allocate(new)
end function

! **************************************************************************** !

function RPF_BVG_liq_ctor(m, Sr) result (new)
  class(RPF_BVG_liq_type), pointer :: new
  PetscReal, intent(in) :: m, Sr
  PetscInt :: error
  PetscReal :: Kr

  error = 0
  if (m <= 0d0 .OR. m >= 1d0) error = error + 1
  if (Sr < 0d0 .OR. Sr >= 1d0) error = error + 2

  if (error == 0) then
    allocate(new)
    new%analytical_derivative_available = .TRUE.
    nullify(new%poly)

    new%m = m
    new%m_rec = 1d0 / m
    new%Sr = Sr
    new%dSe_dSl = 1d0 / (1d0 - Sr)

    call new%RPF_BVG_liq_inline(Slmax,Kr,new%dKr_dSlmax)
    new%dKr_dSlmax = (1d0 - Kr) / epsilon(Slmax)
  else
    nullify(new)
  end if
end function

! **************************************************************************** !

subroutine RPF_BVG_liq_verify(this,name,option)
  implicit none
  class(RPF_BVG_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call PrintErrMsg(option)
  endif   
end subroutine

! **************************************************************************** !

pure subroutine RPF_BVG_liq_inline(this, Sl, Kr, dKr_dSl)
  implicit none
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

subroutine RPF_BVG_liq_relperm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  implicit none
  class(RPF_BVG_liq_type) :: this
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: relative_permeability, dkr_sat
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation > Slmax) then
    relative_permeability = 1d0
    dkr_sat = this%dKr_dSlmax
  else if (liquid_saturation <= this%Sr) then
    relative_permeability = 0d0
    dkr_sat = 0d0
  else
    call this%RPF_BVG_liq_inline(liquid_saturation, relative_permeability, &
                                 dkr_sat)
  end if
end subroutine

! **************************************************************************** !
! Burdine - van Genuchten Gas Relative Permeability Function Methods
! **************************************************************************** !

function RPFBurdineVGGasCreate() result (new)
  implicit none
  class(RPF_BVG_gas_type), pointer :: new
  allocate(new)
end function

! **************************************************************************** !

function RPF_BVG_gas_ctor(m, Slr, Sgr) result (new)
  class(RPF_BVG_gas_type), pointer :: new
  PetscReal, intent(in) :: m, Slr, Sgr
  PetscInt :: error
  PetscReal :: Kr

  error = 0
  if (m <= 0d0 .OR. m >= 1d0) error = error + 1
  if (Slr < 0d0 .OR. Sgr < 0d0 .OR. Slr + Sgr >= 1d0) error = error + 2

  if (error == 0) then
    allocate(new)
    new%analytical_derivative_available = .TRUE.
    nullify(new%poly)

    new%m = m
    new%m_rec = 1d0 / m

    new%Sr = Slr
    new%Sgr = Sgr
    new%dSe_dSl = 1d0 / (1d0 - Slr - Sgr)

    new%Slmax = 1d0 - Sgr - epsilon(Sgr)
    call new%RPF_BVG_gas_inline(new%Slmax, Kr, new%dKr_dSlmax)
    new%dKr_dSlmax = - Kr / epsilon(Sgr)

    new%Slmin = new%Sr + epsilon(Slr)
    new%dKr_dSlmin = -2d0 * new%dSe_dSl
  else
    nullify(new)
  end if
end function

! **************************************************************************** !

subroutine RPF_BVG_gas_verify(this,name,option)
  implicit none
  class(RPF_BVG_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE_VG_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif  
end subroutine

! **************************************************************************** !

pure subroutine RPF_BVG_gas_inline(this, Sl, Kr, dKr_dSl)
  implicit none
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

subroutine RPF_BVG_gas_relperm(this,liquid_saturation, relative_permeability, &
                               dkr_sat,option)
  implicit none
  class(RPF_BVG_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability, dkr_sat
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation > this%Slmax) then
    relative_permeability = 0d0
    dkr_sat = this%dKr_dSlmax
  else if (liquid_saturation < this%Slmin) then
    relative_permeability = 1d0
    dKr_sat = this%dKr_dSlmin
  else
    call this%RPF_BVG_gas_inline(liquid_saturation, relative_permeability, &
                               dkr_sat)
  end if
end subroutine

end module
