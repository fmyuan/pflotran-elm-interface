module Characteristic_Curves_VG_module
#include "petsc/finclude/petscsys.h"

use Characteristic_Curves_Base_module
use Option_module ! Needed for Verify and unused arguments in Pc and Sl
use PFLOTRAN_Constants_module ! Needed for overridden keyword "name" in verify
implicit none
private

! **************************************************************************** !
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
! Child classes are used to override parent functions to avoid additional
! virtual function calls or branching statements in the inner loop.
! Loop invariant parameters are calculated upon construction for memoization.
!
! References:
! Sun, Y., et al. (2010) "Modeling Thermal-Hydrologic Processes for a Heated
! Fractured Rock System: Impact of a Capillary-Pressure Maximum",
! Transp Porous Med 83:501-523. doi: 10.1007/s11242-009-9459-1
!
! Van Genuchten (1980) "A Closed-form Equation for Predicting the Hydraulic
! Conductivity of Unsaturated Soils", Soil Sci Soc Am J 44(5):892-898.
! doi: 10.2136/sssaj1980.03615995004400050002x 
!
! **************************************************************************** !

! As the VG capillary pressure function has an infinite derivative approaching
! saturation, derivatives above Sl_max are approximated using a backwards finite
! difference method using machine epsilon
  PetscReal, private, parameter :: Sl_max = 1d0 - epsilon(Sl_max)

! Parent Saturation Function VG Type
  type, public, extends(sat_func_base_type) :: sat_func_VG_type
    private                          ! Immutable loop-invariant parameters
      PetscReal :: alpha, a_rec      ! Alpha and alpha reciprocal
      PetscReal :: m, m_nrec, m_a2   ! M, negative reciprocal M, M add 2
      PetscReal :: n, n_rec, n_m1    ! N, reciprocal N, and N minus 1
      PetscReal :: mn_a1             ! MN product add 1
      PetscReal :: Sl_span, dSe_dSl  ! Effective saturation span and reciprocal
      PetscReal :: dSe_mndSl         ! Coefficent of VG derivative
      PetscReal :: dPc_dSlmax        ! Finite difference limits at saturation
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

! VG Relative Permeability Function Types
! TODO presently, does not account for unsaturated extensions nor accurately
! reflects the Burdine case. I.e. the saturation functions assume m = 1 - 1/n
! In contrast, the Burdine closed-form solution requires m = 1 - 2/n

! Omitting liquid saturation below the residual, there is still *some* deviation
! in the area under the saturation curve from Sr to Sj. The total area under
! the curve should be found at construction. Integrate 1/h(Sl) from Sr to Sj and
! 1/h(Sl from Sj to 1. This gives the correct denominator.

! Save the value of Int(1/h(Sl)) from Sr to Sj. There is *probably* an easy
! integral for the unsaturated extension. Then also, find the value for the VG integral from Sr to Sj, so that it can be subtracted from the normal MVG relation

! if (Sl < Sr) then
!    Kr = 0d0
! else if (Sl < Sj) then
!    f(Sl) = integrate(1/extension,Sr,Sl)
!    Kr = sqrt(Sl)*(f(Sl) / f(1d0))**2
! else if (Sl > 1d0) then
!    Kr = 1d0
! else
!    f(Sl)= integrate(1/vangenuchten,Sr,Sl) - integrate(1/vangenuchten,Sr,Sj) + integrate(1/extension,Sr,Sj)
!    Kr = sqrt(Sl)*(f(Sl) / f(1d0))**2
! end if
   
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_VG_liq_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPFMualemVGLiqInit
    procedure, public :: Verify => RPFMualemVGLiqVerify
    procedure, public :: SetupPolynomials => RPFMualemVGSetupPolynomials
    procedure, public :: RelativePermeability => RPFMualemVGLiqRelPerm
  end type

  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_VG_gas_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPFMualemVGGasInit
    procedure, public :: Verify => RPFMualemVGGasVerify
    procedure, public :: RelativePermeability => RPFMualemVGGasRelPerm
  end type

  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_VG_liq_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPFBurdineVGLiqInit
    procedure, public :: Verify => RPFBurdineVGLiqVerify
    procedure, public :: SetupPolynomials => RPFBurdineVGSetupPolynomials
    procedure, public :: RelativePermeability => RPFBurdineVGLiqRelPerm
  end type

  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_VG_gas_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPFBurdineVGGasInit
    procedure, public :: Verify => RPFBurdineVGGasVerify
    procedure, public :: RelativePermeability => RPFBurdineVGGasRelPerm
  end type

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
            RPFMualemVGGasCreate, &
            RPFBurdineVGLiqCreate, &
            RPFBurdineVGGasCreate

! Public VG Relative Permeability method accessed by child class
  public :: RPFMualemVGLiqRelPerm

contains
! **************************************************************************** !
! Common VG Saturation Function Methods
! **************************************************************************** !
function SF_VG_Create() result (this)
  ! Creates a dummy van Genutchten saturation function object
  ! The CharacteristicCurvesRead function is passed an object with a generic type
  ! Once the block is parsed, a new object is created that replaces this object.
  ! If object creation by the parser is delayed, this method will not be needed
  implicit none
  class(sat_func_VG_type), pointer :: this
  allocate(this)
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
  PetscReal :: Pc1, Pc2
  PetscReal :: dPc_dSl

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
    call this%SF_VG_Pc_inline(Sl_max,Pc1,dPc_dSl)
    call this%SF_VG_Pc_inline(Sl_max-epsilon(Sl_max),Pc2,dPc_dSl)
    this%dPc_dSlmax = -Pc1 / epsilon(Sl_max)
    this%dSl_dPcmin = -epsilon(Sl_max) / Pc1
    this%d2Sl_dPc2min = epsilon(Sl_max)*(2d0-Pc2/Pc1)/Pc1**2
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
  PetscReal :: aPc_n
  PetscReal :: Se_mrtrec
  PetscReal :: d2Se_mndPc2

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
  PetscReal, intent(out)  :: capillary_pressure
  PetscReal, intent(out)  :: dPc_dSatl
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation <= this%Sr) then
    ! Unsaturated limit
    capillary_pressure = huge(capillary_pressure)
    dPc_dSatl = -huge(capillary_pressure)
  else if (liquid_saturation > Sl_max) then
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
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
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
  PetscReal, intent(in)   :: liquid_saturation
  PetscReal, intent(out)  :: capillary_pressure
  PetscReal, intent(out)  :: dPc_dSatl
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation < this%Sj) then
    ! Unsaturated limit
    capillary_pressure = this%Pcmax
    dPc_dSatl = 0d0
  else if (liquid_saturation > Sl_max) then
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
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
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
  PetscReal :: Se, Se_mrtrec
  PetscReal :: Sa, Sb, Pe, dPj_dSj

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
  PetscReal, intent(out)  :: capillary_pressure
  PetscReal, intent(out)  :: dPc_dSatl
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
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
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
  PetscReal :: Se, Se_mrtrec
  PetscReal :: Sa, Sb, Pe

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
  PetscReal, intent(out)  :: capillary_pressure
  PetscReal, intent(out)  :: dPc_dSatl
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
  else if (liquid_saturation > Sl_max) then
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
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
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
! VG Mualem Liquid Relative Permeability Function Methods
! **************************************************************************** !

function RPFMualemVGLiqCreate()
  ! Creates the van Genutchten Mualem relative permeability function object
  implicit none
  class(rpf_Mualem_vg_liq_type), pointer :: RPFMualemVGLiqCreate
  allocate(RPFMualemVGLiqCreate)
  call RPFMualemVGLiqCreate%Init()
end function

! **************************************************************************** !

subroutine RPFMualemVGLiqInit(this)
  ! Initializes the van Genutchten Mualem relative permeability function object
  implicit none
  class(rpf_Mualem_VG_liq_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  this%analytical_derivative_available = .TRUE.
end subroutine

! **************************************************************************** !
subroutine RPFMualemVGLiqVerify(this,name,option)
  implicit none
  class(rpf_Mualem_VG_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM_VG_LIQ'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call PrintErrMsg(option)
  end if
end subroutine

! **************************************************************************** !
  
subroutine RPFMualemVGSetupPolynomials(this,option,error_string)
  ! Sets up polynomials for smoothing Mualem - van Genuchten relative 
  ! permeability function
  use Utility_module
  implicit none
  class(rpf_Mualem_VG_liq_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  PetscReal :: b(4)
  PetscReal :: one_over_m, Se_one_over_m, m
  
  this%poly => PolynomialCreate()
  ! fill matix with values
  this%poly%low = 0.99d0  ! just below saturated
  this%poly%high = 1.d0   ! saturated
  
  m = this%m
  one_over_m = 1.d0/m
  Se_one_over_m = this%poly%low**one_over_m
  b(1) = 1.d0
  b(2) = sqrt(this%poly%low)*(1.d0-(1.d0-Se_one_over_m)**m)**2.d0
  b(3) = 0.d0
  b(4) = 0.5d0*b(2)/this%poly%low+ &
          2.d0*this%poly%low**(one_over_m-0.5d0)* &
          (1.d0-Se_one_over_m)**(m-1.d0)* &
          (1.d0-(1.d0-Se_one_over_m)**m)
  
  call CubicPolynomialSetup(this%poly%high,this%poly%low,b)
  this%poly%coefficients(1:4) = b(1:4)
end subroutine

! **************************************************************************** !

subroutine RPFMualemVGLiqRelPerm(this,liquid_saturation, &
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
  use Utility_module

  implicit none

  class(rpf_Mualem_VG_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option

  PetscReal :: Se
  PetscReal :: one_over_m
  PetscReal :: Se_one_over_m
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
  
  if (associated(this%poly)) then
    if (Se > this%poly%low) then
      call CubicPolynomialEvaluate(this%poly%coefficients, &
                                   Se,relative_permeability,dkr_Se)
      return 
    endif
  endif
  
  one_over_m = 1.d0/this%m
  Se_one_over_m = Se**one_over_m
  relative_permeability = sqrt(Se)*(1.d0-(1.d0-Se_one_over_m)**this%m)**2.d0
  dkr_Se = 0.5d0*relative_permeability/Se+ &
            2.d0*Se**(one_over_m-0.5d0)* &
                (1.d0-Se_one_over_m)**(this%m-1.d0)* &
                (1.d0-(1.d0-Se_one_over_m)**this%m)

  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat
end subroutine

! **************************************************************************** !
! VG Mualem Gas Relative Permeability Function Methods
! **************************************************************************** !

function RPFMualemVGGasCreate()
  ! Creates the van Genutchten Mualem gas relative permeability function object
  implicit none
  class(rpf_Mualem_VG_gas_type), pointer :: RPFMualemVGGasCreate
  
  allocate(RPFMualemVGGasCreate)
  call RPFMualemVGGasCreate%Init()
end function

! **************************************************************************** !

subroutine RPFMualemVGGasInit(this)
  ! Initializes the van Genutchten Mualem gas relative permeability function 
  ! object
  implicit none
  class(rpf_Mualem_VG_gas_type) :: this

  call RPFBaseInit(this)
  this%analytical_derivative_available = .TRUE.
end subroutine

! **************************************************************************** !

subroutine RPFMualemVGGasVerify(this,name,option)
  implicit none
  class(rpf_Mualem_VG_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM_VG_GAS'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call PrintErrMsg(option)
  endif 
end subroutine
                          
! **************************************************************************** !

subroutine RPFMualemVGGasRelPerm(this,liquid_saturation, &
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
  implicit none
  class(rpf_Mualem_VG_gas_type) :: this
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
  relative_permeability = sqrt(Seg)*(1.d0-Se**(1.d0/this%m))**(2.d0*this%m)
  ! Mathematica analytical solution (Heeho Park)
  dkr_Se = -(1.d0-Se**(1.d0/this%m))**(2.d0*this%m)/(2.d0*sqrt(Seg)) &
          - 2.d0*sqrt(Seg)*Se**(1.d0/this%m-1.d0) &
          * (1.d0-Se**(1.d0/this%m))**(2.d0*this%m-1.d0)
  dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_Se * dSe_sat
end subroutine

! **************************************************************************** !
! VG Burdine Liquid Relative Permeability Function Methods
! **************************************************************************** !
function RPFBurdineVGLiqCreate()
  ! Creates the van Genutchten Mualem relative permeability function object
  implicit none
  class(rpf_burdine_vg_liq_type), pointer :: RPFBurdineVGLiqCreate
  
  allocate(RPFBurdineVGLiqCreate)
  call RPFBurdineVGLiqCreate%Init()
end function

! **************************************************************************** !

subroutine RPFBurdineVGLiqInit(this)
  ! Initializes the van Genutchten Mualem relative permeability function object
  implicit none
  class(rpf_Burdine_VG_liq_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  this%analytical_derivative_available = .TRUE.
end subroutine

! **************************************************************************** !

subroutine RPFBurdineVGLiqVerify(this,name,option)
  implicit none
  class(rpf_Burdine_VG_liq_type) :: this
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

subroutine RPFBurdineVGSetupPolynomials(this,option,error_string)
  ! Sets up polynomials for smoothing Burdine - van Genuchten relative 
  ! permeability function
  use Utility_module
  implicit none
  class(rpf_Burdine_VG_liq_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  PetscReal :: b(4)
  PetscReal :: one_over_m, Se_one_over_m, m

  this%poly => PolynomialCreate()
  ! fill matix with values
  this%poly%low = 0.99d0  ! just below saturated
  this%poly%high = 1.d0   ! saturated
  
  m = this%m
  one_over_m = 1.d0/m
  Se_one_over_m = this%poly%low**one_over_m
  b(1) = 1.d0
  b(2) = this%poly%low*this%poly%low*(1.d0-(1.d0-Se_one_over_m)**this%m)
  b(3) = 0.d0
  b(4) = 2.d0*b(2)/this%poly%low + &
         this%poly%low*Se_one_over_m*(1.d0-Se_one_over_m)**(this%m-1.d0)

  call CubicPolynomialSetup(this%poly%high,this%poly%low,b)
  this%poly%coefficients(1:4) = b(1:4)
end subroutine

! **************************************************************************** !

subroutine RPFBurdineVGLiqRelPerm(this,liquid_saturation, &
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
  use Utility_module
  implicit none
  class(rpf_Burdine_VG_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: one_over_m
  PetscReal :: Se_one_over_m
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
  
  if (associated(this%poly)) then
    if (Se > this%poly%low) then
      call CubicPolynomialEvaluate(this%poly%coefficients, &
                                   Se,relative_permeability,dkr_Se)
      return
    endif
  endif

  one_over_m = 1.d0/this%m
  Se_one_over_m = Se**one_over_m
  relative_permeability = Se*Se*(1.d0-(1.d0-Se_one_over_m)**this%m)
  dkr_Se = 2.d0*relative_permeability/Se + &
                 Se*Se_one_over_m*(1.d0-Se_one_over_m)**(this%m-1.d0)
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat
end subroutine

! **************************************************************************** !
! VG Burdine Liquid Relative Permeability Function Methods
! **************************************************************************** !

function RPFBurdineVGGasCreate()
  ! Creates the Brooks-Corey Burdine gas relative permeability function object
  implicit none
  class(rpf_Burdine_VG_gas_type), pointer :: RPFBurdineVGGasCreate
  
  allocate(RPFBurdineVGGasCreate)
  call RPFBurdineVGGasCreate%Init()
end function

! **************************************************************************** !

subroutine RPFBurdineVGGasInit(this)
  ! Initializes the Brooks-Corey Burdine gas relative permeability function 
  ! object
  implicit none
  class(rpf_Burdine_VG_gas_type) :: this

  call RPFBaseInit(this)
  this%analytical_derivative_available = .TRUE.
end subroutine

! **************************************************************************** !

subroutine RPFBurdineVGGasVerify(this,name,option)
  implicit none
  class(rpf_Burdine_VG_gas_type) :: this
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

subroutine RPFBurdineVGGasRelPerm(this,liquid_saturation, &
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
  implicit none
  class(rpf_Burdine_VG_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
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
  relative_permeability = Seg*Seg*(1.d0-Se**(1.d0/this%m))**this%m
  dkr_Se = -Seg**2.d0*Se**(1.d0/this%m-1.d0) &
          *(1.d0-Se**(1.d0/this%m))**(this%m-1.d0) &
          - 2.d0*Seg*(1.d0-Se**(1.d0/this%m))**this%m
  dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_Se * dSe_sat
end subroutine

end module
