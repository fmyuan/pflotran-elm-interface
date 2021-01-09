module AuxVars_TOWG_module

  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module
  use AuxVars_FlowEnergy_module
  use petscsys

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscBool, public :: towg_analytical_derivatives = PETSC_FALSE
  PetscBool, public :: towg_analytical_derivatives_compare = PETSC_FALSE

  PetscReal, public :: towg_dcomp_tol = 1.d-1
  PetscReal, public :: towg_dcomp_reltol = 1.d-1


  type, public, extends(auxvar_flow_energy_type) :: auxvar_towg_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    PetscBool :: hastl3p_test_object
    type(tl_auxvar_type), pointer :: tl=>null()
    type(bo_auxvar_type), pointer :: bo=>null()
    type(tl3_test_type), pointer  :: tl3TEST=>null()

    PetscBool :: has_TL_test_object
    type(tl_auxvar_testing_type), pointer :: tlT=>null()
  contains
    procedure, public :: Init => AuxVarTOWGInit
    procedure, public :: Strip => AuxVarTOWGStrip
    procedure, public :: InitTL
    procedure, public :: StripTL
    procedure, public :: InitBO
    procedure, public :: StripBO
  end type auxvar_towg_type

  type, public ::  tl_auxvar_type
    PetscReal :: den_oil_eff_kg
    PetscReal :: den_gas_eff_kg
    PetscBool :: has_derivs
    PetscReal, pointer :: D_den_oil_eff_kg(:)   ! (idof)
    PetscReal, pointer :: D_den_gas_eff_kg(:)   ! (idof)
  end type tl_auxvar_type

  type, public ::  bo_auxvar_type
    PetscReal :: bubble_point
    PetscReal :: xo
    PetscReal :: xg

    PetscBool :: has_derivs
    ! derivatives:
    PetscReal, pointer :: D_xo(:)   ! (idof)
    PetscReal, pointer :: D_xg(:)   ! (idof)
  end type bo_auxvar_type

  type, public :: tl3_test_type
    PetscReal :: denotl,dengtl,viscotl,viscgtl,krotl,krgtl,krh
    PetscReal,pointer :: D_denotl(:),D_dengtl(:),D_viscotl(:),D_viscgtl(:),D_krotl(:),D_krgtl(:) 
    PetscReal,pointer :: D_krh(:)
  end type tl3_test_type

    type, public ::  tl_auxvar_testing_type
    !!! hack for testing: store all the intermediate tl variables
    !!! so we can test analytical against numerical derivs
    PetscReal :: krotl,krgtl,viscotl,viscgtl,denotl,dengtl
    PetscReal :: krstl,viscstl,denstl
    PetscReal :: krom,krgm,krsm,krvm,kroi,krog,krow
    PetscReal :: fm,viso,viss,visg
    PetscReal :: uoil,uvap
    PetscReal :: cellpres
    PetscReal  :: krsi,krgi
    PetscReal  :: denos,dengs,denogs,denos_pre,denog

    PetscReal, pointer :: D_krotl(:),D_krgtl(:),D_viscotl(:),D_viscgtl(:)
    PetscReal, pointer :: D_denotl(:),D_dengtl(:)
    PetscReal, pointer :: D_krstl(:),D_viscstl(:),D_denstl(:)
    PetscReal, pointer :: D_krom(:),D_krgm(:),D_krsm(:),D_krvm(:),D_kroi(:)
    PetscReal, pointer :: D_krog(:),D_krow(:)
    PetscReal, pointer :: D_fm(:),D_viso(:),D_visg(:),D_viss(:)
    PetscReal, pointer :: D_uoil(:),D_uvap(:)
    PetscReal, pointer :: D_cellpres(:)
    PetscReal, pointer :: D_krsi(:),D_krgi(:)
    PetscReal, pointer :: D_denos(:),D_dengs(:),D_denogs(:),D_denos_pre(:),D_denog(:)

  end type tl_auxvar_testing_type

  public :: AuxVarTOWGStrip

contains

! ************************************************************************** !

subroutine AuxVarTOWGInit(this,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: PAolo Orsini
  ! Date: 11/07/16
  ! 

  use Option_module

  implicit none
  
  class(auxvar_towg_type) :: this
  type(option_type) :: option

  this%effective_porosity = 0.d0
  this%pert = 0.d0

  !this%has_TL_test_object = option%flow%numerical_derivatives_compare_TL_intermediates
  this%has_TL_test_object = PETSC_TRUE

  if (this%has_TL_test_object) then
    allocate(this%tlT)

    this%tlT%krotl = 0.d0
    this%tlT%krgtl = 0.d0
    this%tlT%viscotl = 0.d0
    this%tlT%viscgtl = 0.d0
    this%tlT%denotl = 0.d0
    this%tlT%dengtl = 0.d0

    this%tlT%krstl = 0.d0
    this%tlT%viscstl = 0.d0
    this%tlT%denstl = 0.d0

    this%tlT%krom = 0.d0
    this%tlT%krgm = 0.d0
    this%tlT%krsm = 0.d0
    this%tlT%krvm = 0.d0

    this%tlT%fm = 0.d0
    this%tlT%viso = 0.d0

    allocate(this%tlT%D_krotl(option%nflowdof))
    this%tlT%D_krotl = 0.d0
    allocate(this%tlT%D_krgtl(option%nflowdof))
    this%tlT%D_krgtl = 0.d0
    allocate(this%tlT%D_viscotl(option%nflowdof))
    this%tlT%D_viscotl = 0.d0
    allocate(this%tlT%D_viscgtl(option%nflowdof))
    this%tlT%D_viscgtl = 0.d0
    allocate(this%tlT%D_denotl(option%nflowdof))
    this%tlT%D_denotl = 0.d0
    allocate(this%tlT%D_dengtl(option%nflowdof))
    this%tlT%D_dengtl = 0.d0

    allocate(this%tlT%D_krstl(option%nflowdof))
    this%tlT%D_krstl = 0.d0
    allocate(this%tlT%D_viscstl(option%nflowdof))
    this%tlT%D_viscstl = 0.d0
    allocate(this%tlT%D_denstl(option%nflowdof))
    this%tlT%D_denstl = 0.d0

    allocate(this%tlT%D_krom(option%nflowdof))
    this%tlT%D_krom = 0.d0
    allocate(this%tlT%D_krgm(option%nflowdof))
    this%tlT%D_krgm = 0.d0
    allocate(this%tlT%D_krsm(option%nflowdof))
    this%tlT%D_krsm = 0.d0
    allocate(this%tlT%D_krvm(option%nflowdof))
    this%tlT%D_krvm = 0.d0
    allocate(this%tlT%D_kroi(option%nflowdof))
    this%tlT%D_kroi = 0.d0
    allocate(this%tlT%D_krog(option%nflowdof))
    this%tlT%D_krog = 0.d0
    allocate(this%tlT%D_krow(option%nflowdof))
    this%tlT%D_krow = 0.d0

    allocate(this%tlT%D_fm(option%nflowdof))
    this%tlT%D_fm = 0.d0
    allocate(this%tlT%D_viso(option%nflowdof))
    this%tlT%D_viso = 0.d0
    allocate(this%tlT%D_visg(option%nflowdof))
    this%tlT%D_visg = 0.d0
    allocate(this%tlT%D_viss(option%nflowdof))
    this%tlT%D_viss = 0.d0

    allocate(this%tlT%D_uoil(option%nflowdof))
    this%tlT%D_uoil = 0.d0
    allocate(this%tlT%D_uvap(option%nflowdof))
    this%tlT%D_uvap = 0.d0

    allocate(this%tlT%D_cellpres(option%nflowdof))
    this%tlT%D_cellpres= 0.d0

    allocate(this%tlT%D_krsi(option%nflowdof))
    this%tlT%D_krsi= 0.d0

    allocate(this%tlT%D_krgi(option%nflowdof))
    this%tlT%D_krgi= 0.d0

    allocate(this%tlT%D_denos(option%nflowdof))
    this%tlT%D_denos= 0.d0
    allocate(this%tlT%D_dengs(option%nflowdof))
    this%tlT%D_dengs= 0.d0
    allocate(this%tlT%D_denogs(option%nflowdof))
    this%tlT%D_denogs= 0.d0

    allocate(this%tlT%D_denos_pre(option%nflowdof))
    this%tlT%D_denos_pre= 0.d0

    allocate(this%tlT%D_denog(option%nflowdof))
    this%tlT%D_denog= 0.d0
  endif

  call AuxVarFlowInit(this,option)

  call AuxVarFlowEnergyInit(this,option)

  this%istate_store = 0

end subroutine AuxVarTOWGInit

! ************************************************************************** !

subroutine InitTL(this,option)
  ! 
  ! Initialize auxiliary object or Todd Longstaff model
  ! 
  ! Author: Paolo Orsini
  ! Date: 07/06/17
  ! 

  use Option_module

  implicit none
  
  class(auxvar_towg_type) :: this
  type(option_type) :: option

  allocate(this%tl)

  nullify(this%tl%D_den_oil_eff_kg)
  nullify(this%tl%D_den_gas_eff_kg)


  this%tl%den_oil_eff_kg = 0.0
  this%tl%den_gas_eff_kg = 0.0

  if (.NOT. option%flow%numerical_derivatives) then
    this%tl%has_derivs = PETSC_TRUE

    allocate(this%tl%D_den_oil_eff_kg(option%nflowdof))
    this%tl%D_den_oil_eff_kg= 0.d0
    allocate(this%tl%D_den_gas_eff_kg(option%nflowdof))
    this%tl%D_den_gas_eff_kg= 0.d0
  endif

    this%hastl3p_test_object = PETSC_FALSE
  if (option%flow%numerical_derivatives_compare) then
    this%hastl3p_test_object = PETSC_TRUE

    allocate(this%tl3TEST)

    nullify(this%tl3TEST%D_denotl)
    nullify(this%tl3TEST%D_dengtl)
    nullify(this%tl3TEST%D_viscotl)
    nullify(this%tl3TEST%D_viscgtl)
    nullify(this%tl3TEST%D_krotl)
    nullify(this%tl3TEST%D_krgtl)
    nullify(this%tl3TEST%D_krh)

    allocate(this%tl3TEST%D_denotl(option%nflowdof))
    this%tl3TEST%D_denotl= 0.d0
    allocate(this%tl3TEST%D_dengtl(option%nflowdof))
    this%tl3TEST%D_dengtl= 0.d0
    allocate(this%tl3TEST%D_viscotl(option%nflowdof))
    this%tl3TEST%D_viscotl= 0.d0
    allocate(this%tl3TEST%D_viscgtl(option%nflowdof))
    this%tl3TEST%D_viscgtl= 0.d0
    allocate(this%tl3TEST%D_krotl(option%nflowdof))
    this%tl3TEST%D_krotl= 0.d0
    allocate(this%tl3TEST%D_krgtl(option%nflowdof))
    this%tl3TEST%D_krgtl= 0.d0
    allocate(this%tl3TEST%D_krh(option%nflowdof))
    this%tl3TEST%D_krh= 0.d0
  endif


end subroutine InitTL

!--Routine to initialise the black oil substructure----------------------------

subroutine InitBO(this,option)

!------------------------------------------------------------------------------
! Used in TOWG_BLACK_OIL and TOWG_SOLVENT_TL, initialises the bo sub-structure
! of auxvars. Contains the bubble point and the oil mole fractions.
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Sep 2017
!------------------------------------------------------------------------------

  use Option_module

  implicit none

  class(auxvar_towg_type) :: this
  type(option_type) :: option

  allocate(this%bo)

  this%bo%bubble_point=0.0
  this%bo%xg          =0.0
  this%bo%xo          =0.0

  nullify(this%bo%D_xo)
  nullify(this%bo%D_xg)

  !if (towg_analytical_derivatives) then
  if (.NOT. option%flow%numerical_derivatives) then
    this%bo%has_derivs= PETSC_TRUE

    allocate(this%bo%D_xo(option%nflowdof))
    this%bo%D_xo = 0.d0
    allocate(this%bo%D_xg(option%nflowdof))
    this%bo%D_xg = 0.d0

  else
    this%bo%has_derivs = PETSC_FALSE
  endif

end subroutine InitBO

!--Routine to initialise the TOWG substructure-------------------------------

subroutine AuxVarTOWGStrip(this)
  ! 
  ! AuxVarTOWGStrip: Deallocates a towg auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/30/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_towg_type) :: this

  call AuxVarFlowStrip(this)

  call AuxVarFlowEnergyStrip(this)

  if (associated(this%tl)) call this%StripTL()
  if (associated(this%bo)) call this%StripBO()

  if (this%has_TL_test_object) then
    call DeallocateArray(this%tlT%D_krotl)
    call DeallocateArray(this%tlT%D_krgtl)
    call DeallocateArray(this%tlT%D_viscotl)
    call DeallocateArray(this%tlT%D_viscgtl)
    call DeallocateArray(this%tlT%D_denotl)
    call DeallocateArray(this%tlT%D_dengtl)

    call DeallocateArray(this%tlT%D_krstl)
    call DeallocateArray(this%tlT%D_viscstl)
    call DeallocateArray(this%tlT%D_denstl)

    call DeallocateArray(this%tlT%D_krom)
    call DeallocateArray(this%tlT%D_krgm)
    call DeallocateArray(this%tlT%D_krsm)
    call DeallocateArray(this%tlT%D_krvm)
    call DeallocateArray(this%tlT%D_kroi)
    call DeallocateArray(this%tlT%D_krog)
    call DeallocateArray(this%tlT%D_krow)

    call DeallocateArray(this%tlT%D_fm)
    call DeallocateArray(this%tlT%D_viso)
    call DeallocateArray(this%tlT%D_visg)
    call DeallocateArray(this%tlT%D_viss)

    call DeallocateArray(this%tlT%D_uoil)
    call DeallocateArray(this%tlT%D_uvap)

    call DeallocateArray(this%tlT%D_cellpres)

    call DeallocateArray(this%tlT%D_krsi)
    call DeallocateArray(this%tlT%D_krgi)

    call DeallocateArray(this%tlT%D_denos)
    call DeallocateArray(this%tlT%D_dengs)
    call DeallocateArray(this%tlT%D_denogs)
    call DeallocateArray(this%tlT%D_denos_pre)
    call DeallocateArray(this%tlT%D_denog)
    deallocate(this%tlT)
  endif

end subroutine AuxVarTOWGStrip


! ************************************************************************** !

subroutine StripTL(this)
  ! 
  ! StripTL: Deallocates a the Todd Longstaff component of 
  !          the towg auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 07/06/17
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_towg_type) :: this

  if (this%has_derivs) then
    call DeallocateArray(this%tl%D_den_oil_eff_kg)
    call DeallocateArray(this%tl%D_den_gas_eff_kg)
  endif

  if (this%hastl3p_test_object) then
    call DeallocateArray(this%tl3TEST%D_denotl)
    call DeallocateArray(this%tl3TEST%D_dengtl)
    call DeallocateArray(this%tl3TEST%D_viscotl)
    call DeallocateArray(this%tl3TEST%D_viscotl)
    call DeallocateArray(this%tl3TEST%D_krotl)
    call DeallocateArray(this%tl3TEST%D_krgtl)
    call DeallocateArray(this%tl3TEST%D_krh)
  endif


  deallocate(this%tl)

end subroutine StripTL

!--Routine to strip the black oil substructure---------------------------------

subroutine StripBO(this)

!------------------------------------------------------------------------------
! Used in TOWG_BLACK_OIL and TOWG_SOLVENT_TL.
! Deallocate the bo sub-structure of auxvars
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Sep 2017
!------------------------------------------------------------------------------

  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_towg_type) :: this

  if (this%has_derivs) then
    call DeallocateArray(this%bo%D_xo)
    call DeallocateArray(this%bo%D_xg)
  endif

  deallocate(this%bo)

end subroutine StripBO

! ************************************************************************** !

end module AuxVars_TOWG_module

