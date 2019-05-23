module Characteristic_Curves_Base_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private
  
  PetscReal, parameter, public :: DEFAULT_PCMAX = 1.d9
  
  type, public :: polynomial_type
    PetscReal :: low
    PetscReal :: high
    PetscReal :: coefficients(4)
  end type polynomial_type  
  
!-----------------------------------------------------------------------------
!-- Saturation Functions -----------------------------------------------------
!-----------------------------------------------------------------------------
  type, public :: sat_func_base_type
    type(polynomial_type), pointer :: sat_poly
    type(polynomial_type), pointer :: pres_poly
    PetscReal :: Sr
    PetscReal :: pcmax
    PetscBool :: analytical_derivative_available
#ifdef SMOOTHING2
    type(polynomial_type), pointer :: sat_poly2      ! dry-end of the curve  ! added by F.-M. Yuan (2017-03-10)
    type(polynomial_type), pointer :: pres_poly2     ! dry-end of the curve  ! added by F.-M. Yuan (2017-03-10)
#endif
  contains
    procedure, public :: Init => SFBaseInit
    procedure, public :: Verify => SFBaseVerify
    procedure, public :: Test => SFBaseTest
    procedure, public :: SetupPolynomials => SFBaseSetupPolynomials
    procedure, public :: CapillaryPressure => SFBaseCapillaryPressure
    procedure, public :: Saturation => SFBaseSaturation
    procedure, public :: D2SatDP2 => SFBaseD2SatDP2
    !added pc function if ice exists
    procedure, public :: IceCapillaryPressure => SFBaseIceCapillaryPressure
  end type sat_func_base_type



!-----------------------------------------------------------------------------
!-- Relative Permeability Functions ------------------------------------------
!-----------------------------------------------------------------------------  
  type, public :: rel_perm_func_base_type
    type(polynomial_type), pointer :: poly
#ifdef SMOOTHING2
    type(polynomial_type), pointer :: poly2     ! dry-end of the curve  ! added by F.-M. Yuan (2017-03-10)
#endif
    PetscReal :: Sr
    PetscBool :: analytical_derivative_available
  contains
    procedure, public :: Init => RPFBaseInit
    procedure, public :: Verify => RPFBaseVerify
    procedure, public :: Test => RPF_Base_Test
    procedure, public :: SetupPolynomials => RPFBaseSetupPolynomials
    procedure, public :: RelativePermeability => RPF_Base_RelPerm
  end type rel_perm_func_base_type
  
  public :: PolynomialCreate, &
            PolynomialDestroy, &
            SFBaseInit, &
            SFBaseVerify, &
            SFBaseTest, &
            SFBaseCapillaryPressure, &
            SFBaseSaturation, &
            RPFBaseInit, &
            RPFBaseVerify, &
            RPF_Base_Test, &
            RPF_Base_RelPerm, &
            SaturationFunctionDestroy, &
            PermeabilityFunctionDestroy

contains

! ************************************************************************** !

function PolynomialCreate()

  implicit none
  
  type(polynomial_type), pointer :: PolynomialCreate  

  allocate(PolynomialCreate)
  PolynomialCreate%low = 0.d0
  PolynomialCreate%high = 0.d0
  PolynomialCreate%coefficients(:) = 0.d0
  
end function PolynomialCreate

! ************************************************************************** !
! ************************************************************************** !

subroutine SFBaseInit(this)

  implicit none
  
  class(sat_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%sat_poly)
  nullify(this%pres_poly)
#ifdef SMOOTHING2
  nullify(this%sat_poly2)
  nullify(this%pres_poly2)
#endif
  this%Sr = UNINITIALIZED_DOUBLE
  this%pcmax = DEFAULT_PCMAX
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine SFBaseInit

! ************************************************************************** !

subroutine SFBaseVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  
  
  if (Uninitialized(this%Sr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                            name)
    call printErrMsg(option)
  endif
  
  if ((.not.this%analytical_derivative_available) .and. &
      (.not.option%flow%numerical_derivatives)) then
    option%io_buffer = 'Analytical derivatives are not available for the &
      &capillary pressure - saturation function chosen: ' // &
      trim(name)
    call printErrMsg(option)
  endif
  
end subroutine SFBaseVerify

! ************************************************************************** !

subroutine RPFBaseInit(this)

  implicit none
  
  class(rel_perm_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%poly)
#ifdef SMOOTHING2
  nullify(this%poly2)
#endif
  this%Sr = UNINITIALIZED_DOUBLE
  this%analytical_derivative_available = PETSC_FALSE
  
end subroutine RPFBaseInit

! ************************************************************************** !

subroutine RPFBaseVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rel_perm_func_base_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  if (Uninitialized(this%Sr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                            name)
    call printErrMsg(option)
  endif
  
  if ((.not.this%analytical_derivative_available) .and. &
      (.not.option%flow%numerical_derivatives)) then
    option%io_buffer = 'Analytical derivatives are not available for the &
      &relative permeability function chosen: ' // trim(name)
    call printErrMsg(option)
  endif
  
end subroutine RPFBaseVerify

! ************************************************************************** !

subroutine SFBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing saturation functions

  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  option%io_buffer = 'SF Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)
  
end subroutine SFBaseSetupPolynomials

! ************************************************************************** !

subroutine RPFBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing relative permeability functions

  use Option_module
  
  implicit none
  
  class(rel_perm_func_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  option%io_buffer = 'RPF Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)
  
end subroutine RPFBaseSetupPolynomials

! ************************************************************************** !

subroutine SFBaseCapillaryPressure(this,liquid_saturation, & 
                                   capillary_pressure,dpc_dsatl,option)
  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseCapillaryPressure must be extended.'
  call printErrMsg(option)
  
end subroutine SFBaseCapillaryPressure

! ************************************************************************** !

subroutine SFBaseIceCapillaryPressure(this, pres_l, tc, &
                           ice_pc, dice_pc_dpres, dice_pc_dt, option)

  ! Ice capillary_pressure as a function of Pres_l, Tc

  use Option_module

  implicit none

  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: pres_l      ! liquid water pressure head (atm. P adjusted), in Pa
  PetscReal, intent(in) :: tc          ! in oC
  PetscReal, intent(out) :: ice_pc, dice_pc_dt, dice_pc_dpres     ! in Pa
  type(option_type), intent(inout) :: option

  call SF_Ice_CapillaryPressure(this, pres_l, tc, &
                           ice_pc, dice_pc_dpres, dice_pc_dt, option)

end subroutine SFBaseIceCapillaryPressure

! ************************************************************************** !

subroutine SFBaseSaturation(this,capillary_pressure, &
                            liquid_saturation,dsat_dpres,option)
  use Option_module

  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseSaturation must be extended.'
  call printErrMsg(option)
  
end subroutine SFBaseSaturation

! ************************************************************************** !

subroutine SFBaseD2SatDP2(this,pc,d2s_dp2,option)

  use Option_module

  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseD2SatDP2 must be extended.'
  call printErrMsg(option)
  
end subroutine SFBaseD2SatDP2

! ************************************************************************** !

subroutine SFBaseTest(this,cc_name,option)

  use Option_module
  use Material_Aux_class

  implicit none
  
  class(sat_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: num_values = 101
  PetscReal :: pc, pc_increment
  PetscReal :: capillary_pressure(num_values)
  PetscReal :: liquid_saturation(num_values)
  PetscReal :: dpc_dsatl(num_values)
  PetscReal :: dpc_dsatl_numerical(num_values)
  PetscReal :: dsat_dpres(num_values)
  PetscReal :: dsat_dpres_numerical(num_values)
  PetscReal :: capillary_pressure_pert
  PetscReal :: liquid_saturation_pert
  PetscReal :: perturbation
  PetscReal :: pert
  PetscReal :: dummy_real
  PetscInt :: count, i
 
  ! calculate saturation as a function of capillary pressure
  ! start at 1 Pa up to maximum capillary pressure
  pc = 1.d0
  pc_increment = 1.d0
  perturbation = 1.d-6
  count = 0
  do
    if (pc > this%pcmax) exit
    count = count + 1
    call this%Saturation(pc,liquid_saturation(count),dsat_dpres(count),option)
    capillary_pressure(count) = pc
    ! calculate numerical derivative dsat_dpres_numerical
    capillary_pressure_pert = pc + pc*perturbation
    call this%Saturation(capillary_pressure_pert,liquid_saturation_pert, &
                         dummy_real,option)
    dsat_dpres_numerical(count) = (liquid_saturation_pert - &
         & liquid_saturation(count))/(pc*perturbation)*(-1.d0) ! dPc/dPres
    ! get next value for pc
    if (pc > 0.99d0*pc_increment*10.d0) pc_increment = pc_increment*10.d0
    pc = pc + pc_increment
  enddo

  write(string,*) cc_name
  string = trim(cc_name) // '_sat_from_pc.dat'
  open(unit=86,file=string)
  write(86,*) '"capillary pressure", "saturation", "dsat/dpres", &
              &"dsat/dpres_numerical"'
  do i = 1, count
    write(86,'(4es14.6)') capillary_pressure(i), liquid_saturation(i), &
                          dsat_dpres(i), dsat_dpres_numerical(i)
  enddo
  close(86)

 ! calculate capillary pressure as a function of saturation
  do i = 1, num_values
    liquid_saturation(i) = dble(i-1)*0.01d0
    if (liquid_saturation(i) < 1.d-7) then
      liquid_saturation(i) = 1.d-7
    else if (liquid_saturation(i) > (1.d0-1.d-7)) then
      liquid_saturation(i) = 1.d0-1.d-7
    endif
    call this%CapillaryPressure(liquid_saturation(i), &
                                capillary_pressure(i),dpc_dsatl(i),option)
    ! calculate numerical derivative dpc_dsatl_numerical
    pert = liquid_saturation(i) * perturbation
    if (liquid_saturation(i) > 0.5d0) then
      pert = -1.d0 * pert
    endif
    liquid_saturation_pert = liquid_saturation(i) + pert
    call this%CapillaryPressure(liquid_saturation_pert, &
                                capillary_pressure_pert,dummy_real,option)
    dpc_dsatl_numerical(i) = (capillary_pressure_pert - &
         & capillary_pressure(i))/pert 
  enddo
  count = num_values

  write(string,*) cc_name
  string = trim(cc_name) // '_pc_from_sat.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "capillary pressure", "dpc/dsat", &
              &dpc_dsat_numerical"'
  do i = 1, count
    write(86,'(4es14.6)') liquid_saturation(i), capillary_pressure(i), &
                          dpc_dsatl(i), dpc_dsatl_numerical(i)
  enddo
  close(86)

end subroutine SFBaseTest

! ************************************************************************** !
! ************************************************************************** !

subroutine RPF_Base_RelPerm(this,liquid_saturation,relative_permeability, &
                            dkr_sat,option)
  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'RPF_Base_RelPerm must be extended.'
  call printErrMsg(option)
  
end subroutine RPF_Base_RelPerm

! ************************************************************************** !

subroutine RPF_Base_Test(this,cc_name,phase,option)

  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  character(len=MAXWORDLENGTH) :: phase
  type(option_type), intent(inout) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, parameter :: num_values = 101
  PetscReal :: perturbation
  PetscReal :: liquid_saturation(num_values)
  PetscReal :: liquid_saturation_pert(num_values)
  PetscReal :: kr(num_values)
  PetscReal :: kr_pert(num_values)
  PetscReal :: dkr_dsat(num_values)
  PetscReal :: dkr_dsat_numerical(num_values)
  PetscReal :: dummy_real(num_values)
  
  perturbation = 1.d-6

  do i = 1, num_values
    liquid_saturation(i) = dble(i-1)*0.01d0
    call this%RelativePermeability(liquid_saturation(i),kr(i),dkr_dsat(i), &
                                   option)
    ! calculate numerical derivative dkr_dsat_numerical
    liquid_saturation_pert(i) = liquid_saturation(i) &
                                + liquid_saturation(i)*perturbation
    call this%RelativePermeability(liquid_saturation_pert(i),kr_pert(i), &
                                   dummy_real(i),option)
    if( i>1 ) then
      dkr_dsat_numerical(i) = (kr_pert(i) - kr(i))/ &
                              (liquid_saturation(i)*perturbation)
    else
! Trap case of i=0 as liquid_saturation is 0 and will otherwise divide by zero
      dkr_dsat_numerical(i) = 0.0
    endif
  enddo

  write(string,*) cc_name
  string = trim(cc_name) // '_' // trim(phase) // '_rel_perm.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "' // trim(phase) // ' relative permeability", "' &
              // trim(phase) // ' dkr/dsat", "' // trim(phase) // &
              ' dkr/dsat_numerical"'
  do i = 1, size(liquid_saturation)
    write(86,'(4es14.6)') liquid_saturation(i), kr(i), dkr_dsat(i), &
                          dkr_dsat_numerical(i)
  enddo
  close(86)

end subroutine RPF_Base_Test

! ************************************************************************** !

subroutine SF_Ice_CapillaryPressure(this, pres_l, tc, &
                                ice_pc, dice_pc_dpres, dice_pc_dt, option)
!
! Computes the ice capillary_pressure as a function of Pres_l, Tc
! Mainly from Painter et al. (2011), Painter and Karra (2014)
!
! Author: Fengming Yuan
!         Based on relevant saturation_functions by Satish K.
!         Revised with freezing-thawing zone smoothing following ATS algorithm
! Date: 06/01/2016
!
! Two Cautions here: (1) when w.r.t. '_dp', means to 'pc'; while w.r.t. '_dpl' or '_dpres', it's to 'presure head'
!                 (2) 'p' or 'pc' or 'cp' is capillary potential and always positive; while 'pres' or 'presl' or 'pw' NOT.


  use Option_module
  use EOS_Water_module
  use Utility_module

  implicit none

  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: pres_l      ! liquid water pressure head (refrence atm. P added), in Pa
  PetscReal, intent(in) :: tc          ! in oC
  PetscReal, intent(out) :: ice_pc     ! in Pa
  PetscReal, intent(out) :: dice_pc_dt, dice_pc_dpres     ! w.r.t. 'oC' or 'pressure head'
  type(option_type) :: option

  PetscReal :: pcgl, pw, tk

  PetscReal, parameter :: beta = 2.33d0           ! dimensionless -- ratio of soil ice surf. tension
  PetscReal, parameter :: T0   = TC2TK            ! freezing-point at standard pressure: in K
  PetscReal, parameter :: Lf   = HEAT_OF_FUSION   ! fusion heat (in J/kg)
  PetscReal :: gammar, alpha, dalpha_drhol
  PetscReal :: rhol, rhol_mol, drhol_dp, drhol_dt
  PetscReal :: rhoi, rhoi_mol, drhoi_dp, drhoi_dt

  PetscReal :: Tf, dTf_dt, dTf_dp
  PetscReal :: tftheta, dtftheta_dt, dtftheta_dp

  PetscReal :: Hfunc, dHfunc, tempreal
  PetscReal :: deltaTf, xTf, a, b, c
  PetscReal :: beta2, T_star, T_star_th, T_star_min, T_star_max
  PetscReal :: dice_pc_dp
  PetscReal, parameter :: dpc_dpres = -1.d0
  PetscReal :: PCGL_MAX_FRZ = 1.d8   ! max. soil matrical potential (Pa) for liq. water exists under super-cooling

  PetscErrorCode :: ierr

  !---------------------
  !
  pcgl = max(0.d0, option%reference_pressure - pres_l)       ! always non-negative (0 = saturated)
  if (pcgl > abs(this%pcmax)) pcgl = this%pcmax

  Tk = tc + T0

  ! if ice module turns on, 2-phase saturation recalculated (liq. and ice) under common 'pw' and 'tc'
  pw = max(option%reference_pressure, pres_l)

  ! --------------------
  ! if not constant rhol
  call EOSWaterDensity(tc, pw, &
                       rhol, rhol_mol, drhol_dp, drhol_dt,ierr)
  rhol = rhol_mol * FMWH2O   ! in eos_water.F90, the conversion from mol->kg may not be consistent for all methods.
  drhol_dp = drhol_dp * FMWH2O
  drhol_dt = drhol_dt * FMWH2O

  call EOSWaterDensityIce(tc, pw, &
                          rhoi_mol, drhoi_dt, drhoi_dp, ierr)
  rhoi = rhoi_mol * FMWH2O   ! in eos_water.F90, the conversion from mol->kg may not be consistent for all methods.
  drhoi_dp = drhoi_dp * FMWH2O
  drhoi_dt = drhoi_dt * FMWH2O

  ! --------------------

    select case (option%ice_model)

      case (PAINTER_EXPLICIT)

        ! explicit model from Painter (Comp. Geosci, 2011)

        ! @T0+deltaTF, ice_pc
        ice_pc = pcgl
        dice_pc_dt = 0.d0
        dice_pc_dp = 1.d0

        if(Tk<=T0) then

          ! @deltaTF, ice_pc = 10*pcgl (reverse function of 'ice_pc' below)
          ! So that, ice_pc will not less than pcgl, i.e. freezing will not allow more ice than total water available
          deltaTf = 10.d0*pcgl*T0/(rhoi*beta*Lf)
          if(Tk>deltaTF+T0) Tk=deltaTF+T0

          ice_pc = - rhoi*beta*Lf*(Tk-T0)/T0               ! in positive value
          dice_pc_dt = - beta*Lf/T0*(rhoi+(Tk-T0)*drhoi_dt)
          dice_pc_dp = - beta*Lf/T0*((Tk-T0)*drhoi_dp)

          ! smoothing from T0 ~ deltaTf
          call HFunctionSmooth(tc+T0, T0, deltaTf+T0, Hfunc, dHfunc)
          dice_pc_dt = dice_pc_dt * Hfunc + (ice_pc - pcgl) * dHfunc
          dice_pc_dp = (dice_pc_dp - 1.0d0) * Hfunc + 1.0d0          ! dHfunc_dp = 0, dpcgl_dp = 1
          ice_pc = (ice_pc - pcgl) * Hfunc + pcgl                    ! do this after derivatives

        endif

      case (PAINTER_KARRA_EXPLICIT, PAINTER_KARRA_EXPLICIT_SMOOTH)

        gammar      = beta*Lf
        alpha       = gammar/T0*rhol
        dalpha_drhol= gammar/T0

        ! freezing-point depression due to capilliary pressure
        PCGL_MAX_FRZ = min(10.d0*alpha, this%pcmax)                 ! 'Tf' of max. 10oK below T0, otherwise, it could be crazy.
        Tf     = T0 - 1.d0/alpha*pcgl !min(PCGL_MAX_FRZ,pcgl)             ! P.-K. Eq.(10), omiga=1/beta. Unit: Kelvin
        dTf_dt = pcgl/alpha/alpha*(dalpha_drhol*drhol_dt)
        dTf_dp = (pcgl*dalpha_drhol*drhol_dp - alpha)/alpha/alpha   ! dpcgl_dp = 1.0
        !if(pcgl>PCGL_MAX_FRZ) then
        !  dTf_dt = 0.d0
        !  dTf_dp = 0.d0
        !endif

        !--------------------------------------------
        deltaTf = 1.0d-50              ! half-width of smoothing zone (by default, nearly NO smoothing)
        if(option%frzthw_halfwidth /= UNINITIALIZED_DOUBLE) deltaTf = option%frzthw_halfwidth
        if (option%ice_model==PAINTER_KARRA_EXPLICIT_SMOOTH .and. deltaTf>1.0d-50) then
          ! if smoothing, only need to calculate the two ends
          if (Tk<T0+deltaTf .and. Tk>=Tf-deltaTf) then
            Tk=Tf-deltaTf
          endif
        endif
        !--------------------------------------------


         ! The following is a slightly-modified version from PKE in saturation_function module
         ! without smoothing of freezing-thawing zone

         ! explicit model from Painter & Karra, VJZ (2014)
         tftheta = -(Tk-T0)/T0                              ! P.-K. Eq.(18): theta: (Tk-T0)/T0, assuming Tf~T0 (ignored FP depression) in Eq. (12)
         dtftheta_dt = -1.0d0/T0
         dtftheta_dp = 0.d0

         ice_pc     = gammar * rhol*tftheta                  ! P.-K. Eq.(18), first term (i.e. ice only), but in positive value as all in CC equations
         tempreal   = rhol*dtftheta_dt+tftheta*drhol_dt
         dice_pc_dt = gammar * tempreal
         tempreal   = rhol*dtftheta_dp+tftheta*drhol_dp
         dice_pc_dp = gammar * tempreal

         ! Heaviside function to truncate Eq. (18) at 'Tf', pcgl@>Tf and ice_pc@<Tf
         ! Tf = T0-1/alpha*pcgl, i.e. -alpha*(Tf-T0)=pcgl, the left side IS exactly ice_pc @Tf
         Hfunc = sign(0.5d0, -(Tk-Tf))+0.5d0
         dice_pc_dt =  dice_pc_dt*Hfunc                          ! no need to adjust dice_pc_dt due to dpcgl_dt = 0
         dice_pc_dp =  dice_pc_dp*Hfunc + (1.d0-Hfunc)           ! dpcgl_dp = 1.0
         ice_pc     =  ice_pc*Hfunc + pcgl*(1.d0-Hfunc)

        !--------------------------------------------
         if (option%ice_model==PAINTER_KARRA_EXPLICIT_SMOOTH .and. deltaTf>1.0d-50) then
            ! using the new Heaveside Smoothing function:
            ! ice_pc@Tf-deltaTf --> ice_pc = pcgl @T0+deltaTf
            !  note: @Tf, ice_pc = pcgl by PKE, so this smoothing actually span over a larger range around Tf~T0)
            call HFunctionSmooth(tc+T0, Tf-deltaTf, Tf+deltaTf, Hfunc, dHfunc)
            dice_pc_dt = dice_pc_dt * Hfunc + (ice_pc - pcgl) * dHfunc
            dice_pc_dp = (dice_pc_dp - 1.0d0) * Hfunc + 1.0d0          ! dHfunc_dp = 0, dpcgl_dp = 1
            ice_pc = (ice_pc - pcgl) * Hfunc + pcgl                    ! do this after derivatives
          endif !(option%ice_model==PAINTER_KARRA_EXPLICIT_SMOOTH)
        !--------------------------------------------

       case (DALL_AMICO)

         ! Model from Dall'Amico (2010) and Dall' Amico et al. (2011)
         ! rewritten following 'saturation_function.F90:SatFuncComputeIceDallAmico()'
         ! NOTE: here only calculate Pc1 and its derivatives
         ice_pc = pcgl
         dice_pc_dt = 0.d0
         dice_pc_dp = 1.d0

         !
         T_star_th  = 5.d-1                       ! unit: Kevin
         beta2      = beta !1.d0                  ! NOT sure why changed to 1. in saturation_function.F90.

         T_star = T0-1.d0/beta2/Lf/rhol*pcgl      ! 'T_star' shall be Tf same as other ice_model, when using same CONSTANTS
         T_star_min = T_star - T_star_th
         T_star_max = T_star

         ! H function
         if (Tk<T_star_min) then
           Hfunc = 1.0d0
           dHfunc= 0.d0
         else if (Tk>T_star_max) then
           Hfunc = 0.0d0
           dHfunc= 0.d0
         else
           tempreal = (Tk-T_star_min)/(T_star_max-T_star_min)
           tempreal = 1.d0 - tempreal*tempreal
           Hfunc = tempreal * tempreal
           dHfunc= -4.0d0*tempreal &
                   *(Tk-T_star_min)/(T_star_max-T_star_min)/(T_star_max-T_star_min)
         endif

         !
         ice_pc     = pcgl - beta2 *  (Tk-T_star)/T_star *Lf*rhol * Hfunc  ! L2030 in saturation_function.F90

         ! saturation_function.F90: L2053 -
         ! dsl_dT = dS1 * (-beta*L_f*rhol_l*H/T_star-beta*theta*L_f*rhol_l*dH_dT)
         ! i.e., dS1* [(-beta*L_f*rhol_l)*(H/T_star+theta*dH_dT)], with theta=(T-T_star)/T_star
         !       for second term:
         !           -beta/T_star*L_f*rhol_l*(H+(T-T_star)*dH_dT)
         !               (because dT_star/dt=0)
         dice_pc_dt = -beta2/T_star*Lf*rhol*        &
                      (Hfunc + (Tk-T_star)*dHfunc)

         ! saturation_function.F90: L2050 -
         !  dsl_dpl = -dS1*(1-T*T_0/T_star/T_star*H), in which
         !   (1) '-dS1' is w.r.t. '-Pc1' actually, because in L2045 it was negatived once
         !   (2) (1-T*T_0/T_star/T_star*H) is 'dp_fh2o_dP' in L.2034.
         !dice_pc_dp = -(1.0d0 - Tk*T0/T_star/T_star*Hfunc)   ! '-1' is for reverse '_dP' to '_dpc'

         ! HOWever, if we do derivative here for 'ice_pc' above, 'ice_pc' shall be
         !   ice_pc = pcgl - beta2*Lf*rhol*Hfunc* (Tk/T_star-1)
         !  w.r.t. 'pc', the second term can be simplified as:
         !              -beta2*Lf*rhol*Hfunc*Tk * (1/T_star)
         !              in which, because: dT_star_dp = -1.d0/beta2/Lf/rhol,
         !                        then:  d(1/T_star)_dp = dT_star_dp/T_star/T_star
         ! i.e.,
         dice_pc_dp = 1.0d0 + Tk/T_star/T_star*Hfunc

       case default

         ! no ice-model defined
         ice_pc = pcgl
         dice_pc_dt = 0.d0
         dice_pc_dp = 1.d0

       end select ! select case (option%ice_model)
     !

  ! --------------------
  dice_pc_dpres = dice_pc_dp*dpc_dpres  ! convert 'pc_dp' to 'pc_dpresl'
  !(NOTE: this is very confusing, so if something wrong better checking this one in all CC equations)

end subroutine SF_Ice_CapillaryPressure


! ************************************************************************** !

subroutine PolynomialDestroy(poly)
  !
  ! Destroys a polynomial smoother
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  !

  implicit none

  type(polynomial_type), pointer :: poly

  if (.not.associated(poly)) return

  deallocate(poly)
  nullify(poly)

end subroutine PolynomialDestroy
  
! ************************************************************************** !

subroutine SaturationFunctionDestroy(sf)
  !
  ! Destroys a saturuation function
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  !

  implicit none

  class(sat_func_base_type), pointer :: sf

  if (.not.associated(sf)) return
  
  call PolynomialDestroy(sf%sat_poly)
  call PolynomialDestroy(sf%pres_poly)
#ifdef SMOOTHING2
  call PolynomialDestroy(sf%sat_poly2)
  call PolynomialDestroy(sf%pres_poly2)
#endif
  deallocate(sf)
  nullify(sf)

end subroutine SaturationFunctionDestroy

! ************************************************************************** !

subroutine PermeabilityFunctionDestroy(rpf)
  ! 
  ! Destroys a saturuation function
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(rel_perm_func_base_type), pointer :: rpf
  
  if (.not.associated(rpf)) return
  
  call PolynomialDestroy(rpf%poly)
#ifdef SMOOTHING2
  call PolynomialDestroy(rpf%poly2)
#endif
  deallocate(rpf)
  nullify(rpf)

end subroutine PermeabilityFunctionDestroy

end module Characteristic_Curves_Base_module
