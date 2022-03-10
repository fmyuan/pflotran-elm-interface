module Material_Transform_module
  !
  ! Models to transform material properties
  !
  ! Author: Alex Salazar III
  ! Date: 02/25/21
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  !---------------------------------------------------------------------------
  type, public :: illitization_base_type
    PetscReal :: ilt_threshold ! temperature threshold to begin illitization
    PetscReal :: ilt_fs0 ! initial fraction of smectite in material
    class(ilt_perm_effects_type), pointer :: ilt_shift_perm ! shift permeability
    class(ilt_kd_effects_type), pointer :: ilt_shift_kd_list ! shift sorption
  contains
    procedure, public :: Verify => ILTBaseVerify
    procedure, public :: Test => ILTBaseTest
    procedure, public :: CalculateILT => ILTBaseIllitization
    procedure, public :: ShiftKd => ILTBaseShiftSorption
    procedure, public :: ShiftPerm => ILTBaseShiftPerm
    procedure, public :: CheckElements => ILTBaseCheckElements
  end type illitization_base_type
  !---------------------------------------------------------------------------
  type, public, extends(illitization_base_type) :: ILT_default_type
    ! Model by W-L Huang, J.M. Longo, & D.R. Pevear, 1993
    PetscReal :: ilt_ea   ! activation energy in J/mol
    PetscReal :: ilt_freq ! frequency term (in L/mol-sec for default)
    PetscReal :: ilt_K_conc ! molar concentration of potassium
  contains
    procedure, public :: Verify => ILTDefaultVerify
    procedure, public :: CalculateILT => ILTDefaultIllitization
    procedure, public :: ShiftKd => ILTShiftSorption
    procedure, public :: ShiftPerm => ILTShiftPerm
    procedure, public :: CheckElements => ILTCheckElements
  end type ILT_default_type
  !---------------------------------------------------------------------------
  type, public, extends(ILT_default_type) :: ILT_general_type
    ! Generalized model by J. Cuadros and J. Linares, 1996
    PetscReal :: ilt_K_exp ! exponent of potassium concentration
    PetscReal :: ilt_exp   ! exponent of smectite fraction
  contains
    procedure, public :: Verify => ILTGeneralVerify
    procedure, public :: CalculateILT => ILTGeneralIllitization
  end type ILT_general_type
  !---------------------------------------------------------------------------
  type :: ilt_perm_effects_type
    PetscReal, pointer :: f_perm(:) ! factors for modifying the permeability
    character(len=MAXWORDLENGTH), pointer :: f_perm_mode ! function type
  end type
  !---------------------------------------------------------------------------
  type :: ilt_kd_effects_type
    PetscInt :: num_elements
    PetscReal, pointer :: f_kd(:,:) ! factors for modifying the kd value
    character(len=MAXWORDLENGTH), pointer :: f_kd_mode(:) ! function type
    character(len=MAXWORDLENGTH), pointer :: f_kd_element(:) ! element affected
  end type
  !---------------------------------------------------------------------------
  type, public :: illitization_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
    PetscBool :: test
    class(illitization_base_type), pointer :: illitization_function
  end type illitization_type
  !---------------------------------------------------------------------------
  type, public :: illitization_auxvar_type
    PetscReal :: fs0    ! initial fraction of smectite in material
    PetscReal :: fs     ! fraction of smectite in material
    PetscReal :: fi     ! fraction of illite in material
    PetscReal :: ts     ! track time of last change in smectite
    PetscBool :: qperm0 ! save initial permeability
    PetscReal :: scale  ! scale factor
    PetscReal, allocatable :: perm0(:) ! intiial permeability
  end type illitization_auxvar_type
  !---------------------------------------------------------------------------
  type, public :: buffer_erosion_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
    PetscBool :: test
    ! class(buffer_erosion_base_type), pointer :: buffer_erosion_model
  end type buffer_erosion_type
  !---------------------------------------------------------------------------
  type, public :: buffer_erosion_auxvar_type
    ! Placeholder for buffer erosion model auxvars
  end type buffer_erosion_auxvar_type
  !---------------------------------------------------------------------------
  type, public :: material_transform_auxvar_type
    class(illitization_auxvar_type), pointer :: il_aux ! auxvars for illitization class
    class(buffer_erosion_auxvar_type), pointer :: be_aux ! auxvars for buffer erosion class
  end type material_transform_auxvar_type
  !---------------------------------------------------------------------------
  type, public :: material_transform_type
    character(len=MAXWORDLENGTH) :: name ! name of material transform
    PetscBool :: test ! not yet implemented
    
    ! Auxiliary variables
    PetscBool :: auxvars_up_to_date
    PetscInt :: num_aux
    type(material_transform_auxvar_type), pointer :: auxvars(:)
    
    ! Classes for material transformations
    class(illitization_type), pointer :: illitization
    class(buffer_erosion_type), pointer :: buffer_erosion
    
    ! Linked list
    type(material_transform_type), pointer :: next
  end type material_transform_type
  !---------------------------------------------------------------------------
  type, public :: material_transform_ptr_type
    class(material_transform_type), pointer :: ptr
  end type material_transform_ptr_type
  !---------------------------------------------------------------------------

  public :: MaterialTransformCreate, &
            MaterialTransformGetID, &
            MaterialTransformCheckILT, &
            MaterialTransformCheckBE, &
            MaterialTransformAddToList, &
            MaterialTransformConvertListToArray, &
            MaterialTransformDestroy, &
            MaterialTransformInputRecord, &
            MaterialTransformRead, &
            MaterialTransformAuxVarInit, &
            IllitizationCreate, &
            IllitizationAuxVarInit, &
            BufferErosionCreate, &
            BufferErosionAuxVarInit

contains

! ************************************************************************** !

subroutine ILTBaseVerify(this,name,option)
  ! 
  ! Checks parameters in the illitization_base_type class
  ! 
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  use Option_module

  implicit none

  class(illitization_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  if (Uninitialized(this%ilt_threshold)) then
    option%io_buffer = 'Illitization temperature threshold must be specified ' &
                     //'for function "'//trim(name)//'".'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%ilt_fs0)) then
    option%io_buffer = 'Initial smectite fraction must be specified ' &
                     //'for function "'//trim(name)//'".'
    call PrintErrMsg(option)
  else
    if (this%ilt_fs0 <= 0.0d0 .or. this%ilt_fs0 > 1.0d0) then
      option%io_buffer = 'Initial smectite fraction for function "' &
        //trim(name)//'" must be nonzero positive number up to 1.'
      call PrintErrMsg(option)
    endif
  endif

end subroutine ILTBaseVerify

! ************************************************************************** !

subroutine ILTBaseIllitization(this,fs,temperature,dt,fi,scale,option)

  use Option_module

  implicit none

  class(illitization_base_type) :: this
  PetscReal, intent(inout) :: fs
  PetscReal, intent(in) :: temperature
  PetscReal, intent(in) :: dt
  PetscReal, intent(out) :: fi
  PetscReal, intent(out) :: scale
  type(option_type), intent(inout) :: option

  fi = 0.0d+0
  scale = 0.0d+0

end subroutine ILTBaseIllitization

! ************************************************************************** !

subroutine ILTBaseTest(this,name,option)
  !
  ! Tests illitization functions using a range of initial smectite contents and
  !   temperatures over geological time
  !
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  use Option_module

  implicit none

  class(illitization_base_type) :: this
  character(len=MAXWORDLENGTH) :: name
  type(option_type), intent(inout) :: option

  ! Test with pertubrations to initial smectite and temperature over time
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: ns = 10
  PetscInt, parameter :: nt = 61
  PetscInt, parameter :: np = 27
  PetscReal, parameter :: perturbation = 1.0d-6
  PetscReal :: deltaSmec
  PetscReal :: deltaTemp
  PetscReal :: smec_vec(ns)
  PetscReal :: temp_vec(nt)
  PetscReal :: time_vec(np)
  PetscReal :: fi(ns,nt,np)
  PetscReal :: dfi_dtemp_numerical(ns,nt,np)
  PetscReal :: perturbed_temp
  PetscReal :: fi_temp_pert
  PetscReal :: smec_min, smec_max
  PetscReal :: temp_min, temp_max
  PetscReal :: sc(ns,nt,np), sc_temp_pert
  PetscReal :: dt,fs0_original
  PetscReal :: fs, fsp
  PetscInt :: i, j, k

  ! thermal conductivity as a function of temp. and liq. sat.
  smec_min = 1.0d-1 ! Minimum fraction smectite
  smec_max = 1.0d+0 ! Maximum fraction smectite
  temp_min = 2.0d+1 ! Minimum temperature in Celcius
  temp_max = 2.6d+2 ! Maximum temperature in Celcius

  deltaSmec = (smec_max - smec_min)/(ns - 1)
  deltaTemp = (temp_max - temp_min)/(nt - 1)

  smec_vec = [(smec_min + i*deltaSmec, i=0,ns-1)]
  temp_vec = [(temp_min + i*deltaTemp, i=0,nt-1)]
  time_vec = (/0.,1.,2.5,5.,7.5,10.,25.,50.,75.,100.,250.,500.,750.,1000., &
               2500.,5000.,7500.,10000.,20000.,30000.,40000.,50000.,60000.,&
               70000.,80000.,90000.,100000./)

  fs0_original = this%ilt_fs0

  do i = 1,ns
    do j = 1,nt
      ! reset base variables to initial
      this%ilt_fs0 = smec_vec(i)
      fs  = smec_vec(i)
      do k = 2,np

        ! get change in time
        dt = time_vec(k) - time_vec(k-1) ! years
        dt = dt*(365.25*24*60*60)        ! convert to seconds

        ! base case with analytical derivatives
        fsp = fs
        call this%CalculateILT(fs,temp_vec(j),dt,fi(i,j,k),sc(i,j,k),option)

        ! calculate numerical derivatives via finite differences
        perturbed_temp = temp_vec(j) * (1.d0 + perturbation)
        call this%CalculateILT(fsp,perturbed_temp,dt,fi_temp_pert, &
                               sc_temp_pert,option)

        dfi_dtemp_numerical(i,j,k) = (fi_temp_pert - fi(i,j,k))/ &
                                      (temp_vec(j)*perturbation)
      enddo
    enddo
  enddo

  write(string,*) name
  string = trim(name) // '_ilt_vs_time_and_temp.dat'
  open(unit=86,file=string)
  write(86,*) '"initial smectite [-]", "temperature [C]", &
               "time [yr]", "illite [-]", "dillite/dT [1/yr]", "scale [-]"'
  do i = 1,ns
    do j = 1,nt
      do k = 2,np
        write(86,'(6(ES14.6))') smec_vec(i), temp_vec(j), time_vec(k), &
             fi(i,j,k),dfi_dtemp_numerical(i,j,k),sc(i,j,k)
      enddo
    enddo
  enddo
  close(86)

  ! reset to original values
  this%ilt_fs0 = fs0_original

end subroutine ILTBaseTest

! ************************************************************************** !

subroutine ILTDestroy(ilf)

  implicit none

  class(illitization_base_type), pointer :: ilf

  if (.not.associated(ilf)) return
  deallocate(ilf)
  nullify(ilf)

end subroutine ILTDestroy

! ************************************************************************** !

function ILTBaseCreate()

  implicit none

  class(illitization_base_type), pointer :: ILTBaseCreate

  allocate(ILTBaseCreate)

  ILTBaseCreate%ilt_threshold  = 0.0d0
  ILTBaseCreate%ilt_fs0        = 0.0d0
  nullify(ILTBaseCreate%ilt_shift_perm)
  nullify(ILTBaseCreate%ilt_shift_kd_list)

end function ILTBaseCreate

! ************************************************************************** !

function ILTDefaultCreate()

  implicit none

  class(ILT_default_type), pointer :: ILTDefaultCreate

  allocate(ILTDefaultCreate)

  ILTDefaultCreate%ilt_threshold  = 0.0d0
  ILTDefaultCreate%ilt_fs0        = 1.0d0
  ILTDefaultCreate%ilt_ea     = UNINITIALIZED_DOUBLE
  ILTDefaultCreate%ilt_freq   = UNINITIALIZED_DOUBLE
  ILTDefaultCreate%ilt_K_conc = UNINITIALIZED_DOUBLE
  nullify(ILTDefaultCreate%ilt_shift_perm)
  nullify(ILTDefaultCreate%ilt_shift_kd_list)

end function ILTDefaultCreate

! ************************************************************************** !

function ILTGeneralCreate()

  implicit none

  class(ILT_general_type), pointer :: ILTGeneralCreate

  allocate(ILTGeneralCreate)

  ILTGeneralCreate%ilt_threshold  = 0.0d0
  ILTGeneralCreate%ilt_fs0        = 1.0d0
  ILTGeneralCreate%ilt_freq       = 1.0d0 ! Default of 1.0 in general model
  ILTGeneralCreate%ilt_ea     = UNINITIALIZED_DOUBLE
  ILTGeneralCreate%ilt_K_conc = UNINITIALIZED_DOUBLE
  ILTGeneralCreate%ilt_K_exp  = UNINITIALIZED_DOUBLE
  ILTGeneralCreate%ilt_exp    = UNINITIALIZED_DOUBLE
  nullify(ILTGeneralCreate%ilt_shift_perm)
  nullify(ILTGeneralCreate%ilt_shift_kd_list)

end function ILTGeneralCreate

! ************************************************************************** !

function ILTPermEffectsCreate()
  ! 
  ! Creates object for modifying permeability from smectite/illite transition
  ! 
  ! Author: Alex Salazar III
  ! Date: 11/11/2021

  implicit none

  class(ilt_perm_effects_type), pointer :: ILTPermEffectsCreate

  allocate(ILTPermEffectsCreate)
  
  nullify(ILTPermEffectsCreate%f_perm)
  nullify(ILTPermEffectsCreate%f_perm_mode)

end function ILTPermEffectsCreate

! ************************************************************************** !

function ILTKdEffectsCreate()
  ! 
  ! Creates object for modifying Kd values from smectite/illite transition
  ! 
  ! Author: Alex Salazar III
  ! Date: 10/06/2021

  implicit none

  class(ilt_kd_effects_type), pointer :: ILTKdEffectsCreate

  allocate(ILTKdEffectsCreate)
  
  ILTKdEffectsCreate%num_elements = UNINITIALIZED_INTEGER
  nullify(ILTKdEffectsCreate%f_kd)
  nullify(ILTKdEffectsCreate%f_kd_mode)
  nullify(ILTKdEffectsCreate%f_kd_element)

end function ILTKdEffectsCreate

! ************************************************************************** !

subroutine ILTDefaultVerify(this,name,option)
  ! 
  ! Checks parameters in the ILT_default_type class
  ! 
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  use Option_module

  implicit none

  class(ILT_default_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'ILLITIZATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'ILLITIZATION_FUNCTION, DEFAULT'
  endif
  call ILTBaseVerify(this,string,option)
  if (Uninitialized(this%ilt_ea)) then
    option%io_buffer = 'Illitization activation energy must be specified in' &
                     //' function "'//trim(string)//'".'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%ilt_freq)) then
    option%io_buffer = 'Illitization frequency term must be specified in' &
                     //' function "'//trim(string)//'".'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%ilt_K_conc)) then
    option%io_buffer = 'Illitization potassium concentration must be ' &
                     //'specified in function "'//trim(string)//'".'
    call PrintErrMsg(option)
  endif

end subroutine ILTDefaultVerify

! ************************************************************************** !

subroutine ILTGeneralVerify(this,name,option)
  ! 
  ! Checks parameters in the ILT_general_type class
  ! 
  ! Author: Alex Salazar III
  ! Date: 06/16/2021
  !
  use Option_module

  implicit none

  class(ILT_general_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'ILLITIZATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'ILLITIZATION_FUNCTION, GENERAL'
  endif
  call ILTBaseVerify(this,string,option)
  call ILTDefaultVerify(this,string,option)
  if (Uninitialized(this%ilt_exp)) then
    option%io_buffer = 'Illitization smectite exponent must be specified in' &
                     //' function "'//trim(string)//'".'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%ilt_K_exp)) then
    option%io_buffer = 'Illitization postassium exponent must be specified in' &
                     //' function "'//trim(string)//'".'
    call PrintErrMsg(option)
  endif

end subroutine ILTGeneralVerify

! ************************************************************************** !

subroutine ILTDefaultIllitization(this,fs,temperature,dt,fi,scale,option)

  use Option_module

  implicit none

  class(ILT_default_type) :: this
  PetscReal, intent(inout) :: fs
  PetscReal, intent(in) :: temperature
  PetscReal, intent(in) :: dt
  PetscReal, intent(out) :: fi
  PetscReal, intent(out) :: scale
  type(option_type), intent(inout) :: option

  PetscReal :: ds   ! change in smectite
  PetscReal :: T    ! temperature in Kelvin
  PetscReal :: rate ! temperature-dependent illitization rate in sec^-1

  ! Model based on W-L Huang, J.M. Longo, & D.R. Pevear,
  !  "An Experimentally Derived Kinetic Model for Smectite-to-Illite Conversion
  !  and its Use as a Geothermometer," Clay and Clay Minerals, vol. 41, no 2., 
  !  pp. 162-177, 1993

  ! Use Kelvin to calculate rate
  T = temperature + 273.15d0

  ! Check if temperature is above threshold for illitization
  if (temperature >= this%ilt_threshold) then
    ! Negative of illitization rate [L/mol-s]
    rate = this%ilt_K_conc * this%ilt_freq * &
      exp(-1.0d0 * this%ilt_ea / (IDEAL_GAS_CONSTANT * T))
  else
    rate = 0.0d0
  endif

  ! Log change in smectite as time proceeds
  ds = rate * dt

  ! Fraction smectite
  fs = fs / (1.0d0 + (fs * ds))

  if (fs > 1.0d0) then
    fs = 1.0d0
  elseif (fs < 0.0d0) then
    fs = 0.0d0
  endif

  ! Fraction illite
  fi = 1.0d0 - fs
  
  ! Calculate scale factor
  scale = ((fi - (1.0d+0 - this%ilt_fs0)) / this%ilt_fs0)

end subroutine ILTDefaultIllitization

! ************************************************************************** !

subroutine ILTGeneralIllitization(this,fs,temperature,dt,fi,scale,option)

  use Option_module

  implicit none

  class(ILT_general_type) :: this
  PetscReal, intent(inout) :: fs
  PetscReal, intent(in) :: temperature
  PetscReal, intent(in) :: dt
  PetscReal, intent(out) :: fi
  PetscReal, intent(out) :: scale
  type(option_type), intent(inout) :: option

  PetscReal :: T    ! temperature in Kelvin
  PetscReal :: rate ! temperature-dependent illitization rate in sec^-1

  ! Model based on J. Cuadros & J. Linares, "Experimental Kinetic Study of the
  !   Smectite-to-Illite Transformation," Geochimica et Cosmochimica Acta, 
  !   vol. 60, no. 3, pp. 439-453, 1996

  ! Use Kelvin to calculate rate
  T = temperature + 273.15d0

  ! Check if temperature is above threshold for illitization
  if (temperature >= this%ilt_threshold) then
    ! Negative of illitization rate [L/mol-s]
    rate = this%ilt_freq * &
      exp(-1.0d0 * this%ilt_ea / (IDEAL_GAS_CONSTANT * T))
  else
    rate = 0.0d0
  endif

  ! Fraction smectite - pivot solution based on choice of exponent
  if (this%ilt_exp == 1.0d0) then
    ! n = 1
    fs = fs * exp(-1.0d0 * rate * (this%ilt_K_conc**this%ilt_K_exp) * dt)
  else
    ! n != 1
    fs = (rate * (this%ilt_K_conc**this%ilt_K_exp) * &
         (this%ilt_exp - 1.0d0) * dt + &
         fs**(1.0d0 - this%ilt_exp))**(1.0d0/(1.0d0 - this%ilt_exp))
  endif

  if (fs > 1.0d0) then
    fs = 1.0d0
  elseif (fs < 0.0d0) then
    fs = 0.0d0
  endif

  ! Fraction illite
  fi = 1.0d0 - fs
  
  ! Calculate scale factor
  scale = ((fi - (1.0d+0 - this%ilt_fs0)) / this%ilt_fs0)

end subroutine ILTGeneralIllitization

! ************************************************************************** !

subroutine ILTBaseShiftSorption(this,kd0,ele,auxvar,option)

  use Option_module
  use Material_Aux_module

  implicit none

  class(illitization_base_type) :: this
  PetscReal, intent(inout) :: kd0
  character(len=MAXWORDLENGTH), intent(in) :: ele
  class(illitization_auxvar_type), intent(in) :: auxvar
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'Illitization function must be extended to modify ' &
                   //'the kd values of elements in UFD Decay.'
  call PrintErrMsgByRank(option)
  
end subroutine ILTBaseShiftSorption

! ************************************************************************** !

subroutine ILTShiftSorption(this,kd0,ele,auxvar,option)
  !
  ! Modifies the kd of selected elements using results from the
  !   illitization model and a user-specified functional form.
  !
  ! Author: Alex Salazar III
  ! Date: 10/21/2021
  !
  use Option_module
  use Material_Aux_module

  implicit none

  class(ILT_default_type) :: this
  PetscReal, intent(inout) :: kd0
  character(len=MAXWORDLENGTH), intent(in) :: ele
  class(illitization_auxvar_type), intent(in) :: auxvar
  type(option_type), intent(inout) :: option
  
  class(ilt_kd_effects_type), pointer :: kdl
  character(len=MAXWORDLENGTH) :: fkdele
  character(len=MAXWORDLENGTH) :: fkdmode
  PetscReal, allocatable :: fkd(:)
  PetscInt :: i, j, k
  PetscReal :: scale, factor
  
  if (.not. associated(this%ilt_shift_kd_list)) return

  ! Find element and functional properties
  j = 0
  kdl => this%ilt_shift_kd_list
  do i = 1, kdl%num_elements
    fkdele = kdl%f_kd_element(i)
    ! If elements match, proceed
    if (trim(fkdele) == trim(ele)) then
      ! Identify function
      fkdmode = kdl%f_kd_mode(i)
      ! Allocate vector of function values
      select case(fkdmode)
        case ('DEFAULT','LINEAR')
          j = 1
        case ('QUADRATIC')
          j = 2
        case ('POWER')
          j = 2
        case ('EXPONENTIAL')
          j = 1
        case default
          option%io_buffer = 'Sorption modification function "' &
                           // trim(fkdmode) &
                           //'" was not found among the available options.'
          call PrintErrMsgByRank(option)
      end select
      allocate(fkd(j))
      ! Populate local vector of function values
      do k = 1, j
        fkd(k) = kdl%f_kd(i,k)
      enddo
      ! Done
      exit
    endif
  enddo
  
  ! Apply function to modify kd
  scale = auxvar%scale
  factor = 1.0d0
  select case(fkdmode)
    case ('DEFAULT','LINEAR')
      factor = 1.0d0 + fkd(1)*scale
    case ('QUADRATIC')
      factor = 1.0d0 + fkd(1)*scale + fkd(2)*(scale**2)
    case ('POWER')
      factor = 1.0d0 + fkd(1)*(scale**fkd(2))
    case ('EXPONENTIAL')
      factor = exp(fkd(1) * scale)
    case default
      option%io_buffer = 'No analytical expression available for sorption ' &
                       //'modification function "'// trim(fkdmode)//'".'
      call PrintErrMsgByRank(option)
  end select
  
  kd0 = kd0 * factor
  
  if (allocated(fkd)) deallocate(fkd)

end subroutine ILTShiftSorption

! ************************************************************************** !

subroutine ILTBaseCheckElements(this,pm_ufd_elements,num,option)

  use Option_module

  implicit none

  class(illitization_base_type) :: this
  PetscInt, intent(in) :: num
  character(len=MAXWORDLENGTH), intent(in) :: pm_ufd_elements(num)
  type(option_type), intent(inout) :: option
  
  return
  
end subroutine ILTBaseCheckElements

! ************************************************************************** !

subroutine ILTCheckElements(this,pm_ufd_elements,num,option)
  !
  ! Ensures that the elements specified for kd modification in ILLITIZATION are
  !   present in the reference list (UFD Decay)
  !
  ! Author: Alex Salazar III
  ! Date: 11/01/2021
  !
  use Option_module

  implicit none

  class(ILT_default_type) :: this
  PetscInt, intent(in) :: num
  character(len=MAXWORDLENGTH), intent(in) :: pm_ufd_elements(num)
  type(option_type), intent(inout) :: option
  
  class(ilt_kd_effects_type), pointer :: kdl
  character(len=MAXWORDLENGTH) :: fkdele1, fkdele2
  PetscInt :: i, j
  PetscBool :: found
  
  if (.not. associated(this%ilt_shift_kd_list)) return
  
  kdl => this%ilt_shift_kd_list
  ! Check for duplicates in list
  do i = 1, kdl%num_elements
    fkdele1 = kdl%f_kd_element(i)
    do j = 1, kdl%num_elements
      if (i == j) cycle
      fkdele2 = kdl%f_kd_element(j)
      if (trim(fkdele1) == trim(fkdele2)) then
        option%io_buffer = 'Duplicate element "'// trim(fkdele1) &
                         //'" has been detected in SHIFT_KD.'
        call PrintErrMsgByRank(option)
      endif
    enddo
  enddo
  ! Check if present in UFD Decay
  do i = 1, kdl%num_elements
    ! Element specified in illitization function
    fkdele1 = kdl%f_kd_element(i)
    found = PETSC_FALSE
    do j = 1, num
      ! Element specified in UFD Decay
      fkdele2 = pm_ufd_elements(j)
      ! If elements match, proceed to next in list
      if (trim(fkdele1) == trim(fkdele2)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not. found) then
      option%io_buffer = 'Element "'// trim(fkdele1) &
                       //'" listed for kd modification was not found among ' &
                       //'the elements in UFD Decay.'
      call PrintErrMsgByRank(option)
    endif
  enddo
  
end subroutine ILTCheckElements

! ************************************************************************** !

subroutine ILTBaseShiftPerm(this,material_auxvar,auxvar,option)

  use Option_module
  use Material_Aux_module

  implicit none

  class(illitization_base_type), intent(inout) :: this
  class(material_auxvar_type), intent(inout) :: material_auxvar
  class(illitization_auxvar_type), intent(inout) :: auxvar
  class(option_type), intent(inout) :: option

  option%io_buffer = 'Illitization function must be extended to modify ' &
                   //'the permeability.'
  call PrintErrMsgByRank(option)

end subroutine ILTBaseShiftPerm

! ************************************************************************** !

subroutine ILTShiftPerm(this,material_auxvar,auxvar,option)
  !
  ! Modifies the permeability tensor using results from the
  !   illitization model and a user-specified functional form.
  !
  ! Author: Alex Salazar III
  ! Date: 11/11/2021
  !

  use Option_module
  use Material_Aux_module

  implicit none

  class(ILT_default_type), intent(inout) :: this
  class(material_auxvar_type), intent(inout) :: material_auxvar
  class(illitization_auxvar_type), intent(inout) :: auxvar
  class(option_type), intent(inout) :: option

  PetscInt  :: ps, i, j, k
  PetscReal :: scale, factor
  PetscReal, allocatable :: fperm(:)
  character(len=MAXWORDLENGTH) :: fpermmode
  class(ilt_perm_effects_type), pointer :: perm

  ! Check whether illitization and permeability modification are active
  if (.not. associated(this%ilt_shift_perm)) return

  ! Assess whether original permeability was saved in the auxvar
  ps = size(material_auxvar%permeability)
  if (.not. auxvar%qperm0) then
    ! allocate(auxvar%perm0(ps))
    auxvar%perm0 = UNINITIALIZED_DOUBLE
    do i = 1, ps
      auxvar%perm0(i) = material_auxvar%permeability(i)
    enddo
    auxvar%qperm0 = PETSC_TRUE
  endif

  ! Find functional properties
  j = 0
  perm => this%ilt_shift_perm
  fpermmode = perm%f_perm_mode
  select case(fpermmode)
    case ('DEFAULT','LINEAR')
      j = 1
    case ('QUADRATIC')
      j = 2
    case ('POWER')
      j = 2
    case ('EXPONENTIAL')
      j = 1
    case default
      option%io_buffer = 'Permeability modification function "' &
                       // trim(fpermmode) &
                       //'" was not found among the available options.'
      call PrintErrMsgByRank(option)
  end select
  allocate(fperm(j))
  do k = 1, j
    fperm(k) = perm%f_perm(k)
  enddo

  ! Apply function to modify permeability tensor
  scale = auxvar%scale
  factor = 1.0d0
  select case(fpermmode)
    case ('DEFAULT','LINEAR')
      factor = 1.0d0 + fperm(1)*scale
    case ('QUADRATIC')
      factor = 1.0d0 + fperm(1)*scale + fperm(2)*(scale**2)
    case ('POWER')
      factor = 1.0d0 + fperm(1)*(scale**fperm(2))
    case ('EXPONENTIAL')
      factor = exp(fperm(1) * scale)
    case default
      option%io_buffer = 'No analytical expression available for permeability '&
                       //'modification function "'// trim(fpermmode)//'".'
      call PrintErrMsgByRank(option)
  end select
  
  do i = 1, ps
    material_auxvar%permeability(i) = auxvar%perm0(i) * factor
  enddo
  
  if (allocated(fperm)) deallocate(fperm)

  ! Save time of permeability modification
  auxvar%ts = option%time

end subroutine ILTShiftPerm

! ************************************************************************** !

function MaterialTransformCreate()
  !
  ! Creates a material transform object
  !
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  implicit none

  class(material_transform_type), pointer :: MaterialTransformCreate
  class(material_transform_type), pointer :: material_transform

  allocate(material_transform)
  material_transform%name = ''
  material_transform%test = PETSC_FALSE
  material_transform%num_aux = 0
  nullify(material_transform%auxvars)
  nullify(material_transform%illitization)
  nullify(material_transform%buffer_erosion)
  nullify(material_transform%next)

  MaterialTransformCreate => material_transform

end function MaterialTransformCreate

! ************************************************************************** !

function IllitizationCreate()
  !
  ! Creates an illitization object
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022
  !
  implicit none

  class(illitization_type), pointer :: IllitizationCreate
  class(illitization_type), pointer :: Illitization

  allocate(Illitization)
  Illitization%name = ''
  Illitization%print_me = PETSC_FALSE
  Illitization%test = PETSC_FALSE
  nullify(Illitization%illitization_function)

  IllitizationCreate => Illitization

end function IllitizationCreate

! ************************************************************************** !

function BufferErosionCreate()
  !
  ! Creates an object for a buffer erosion model
  !
  ! Author: Alex Salazar III
  ! Date: 01/20/2022
  !
  implicit none

  class(buffer_erosion_type), pointer :: BufferErosionCreate
  class(buffer_erosion_type), pointer :: BufferErosion

  allocate(BufferErosion)
  BufferErosion%name = ''
  BufferErosion%print_me = PETSC_FALSE
  BufferErosion%test = PETSC_FALSE
  ! nullify(BufferErosion%buffer_erosion_model)

  BufferErosionCreate => BufferErosion

end function BufferErosionCreate

! ************************************************************************** !

subroutine MaterialTransformRead(this,input,option)
  !
  ! Reads in components of a MATERIAL_TRANSFORM block
  !
  ! Author: Alex Salazar III
  ! Date: 11/09/2021
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(material_transform_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string

  input%ierr = 0
  error_string = 'MATERIAL_TRANSFORM "'//trim(this%name)//'"'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      !------------------------------------------
      case('ILLITIZATION')
        this%illitization => IllitizationCreate()
        call IllitizationRead(this%illitization,input,option)
      !------------------------------------------
      case('BUFFER_EROSION')
        this%buffer_erosion=> BufferErosionCreate()
        call BufferErosionRead(this%buffer_erosion,input,option)
      !------------------------------------------
      case default
        call InputKeywordUnrecognized(input,keyword, &
               'MATERIAL_TRANSFORM "'//trim(this%name)//'"',option)
    end select
    
  enddo
  
  call InputPopBlock(input,option)

end subroutine MaterialTransformRead

! ************************************************************************** !

subroutine IllitizationRead(this,input,option)
  !
  ! Reads in contents of an ILLITIZATION block from MATERIAL_TRANSFORM
  !
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(illitization_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string, verify_string
  class(illitization_base_type), pointer :: illitization_function_ptr

  nullify(illitization_function_ptr)

  input%ierr = 0
  error_string = 'ILLITIZATION'
  
  if (associated(this%illitization_function)) then
    option%io_buffer = 'There may only be one instance of '// &
                       'ILLITIZATION_FUNCTION in ILLITIZATION "'// &
                       trim(this%name)//'".'
    call PrintErrMsg(option)
  endif
  
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      !------------------------------------------
      case('ILLITIZATION_FUNCTION')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option, &
             'ILLITIZATION_FUNCTION',error_string)
        call StringToUpper(word)
        select case(word)
          !-------------------------------------
          case('DEFAULT','HUANG')
            this%illitization_function => ILTDefaultCreate()
          !-------------------------------------
          case('GENERAL','CUADROS_AND_LINARES')
            this%illitization_function => ILTGeneralCreate()
          !-------------------------------------
          case default
            call InputKeywordUnrecognized(input,word, &
                 'ILLITIZATION_FUNCTION',option)
        end select
        call ILTRead(this%illitization_function,input,option)
      !------------------------------------------
      case('TEST')
        this%test = PETSC_TRUE
      !------------------------------------------
      case default
        call InputKeywordUnrecognized(input,keyword,'ILLITIZATION',option)
    end select
  enddo
  call InputPopBlock(input,option)

  verify_string = 'ILLITIZATION(' // trim(this%name) // '),'

  if (associated(this%illitization_function)) then
    call this%illitization_function%Verify(verify_string,option)
  else
    option%io_buffer = 'A illitization function has &
         &not been set under ILLITIZATION "' // &
         trim(this%name) // '". An ILLITIZATION_FUNCTION &
         &block must be specified.'
  endif

end subroutine IllitizationRead

! ************************************************************************** !

subroutine BufferErosionRead(this,input,option)
  !
  ! Reads in contents of a BUFFER_EROSION block from MATERIAL_TRANSFORM
  !
  ! Author: Alex Salazar III
  ! Date: 01/20/2022
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(buffer_erosion_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string, verify_string
  ! class(buffer_erosion_base_type), pointer :: buffer_erosion_model_ptr

  ! nullify(buffer_erosion_model_ptr)

  input%ierr = 0
  error_string = 'BUFFER_EROSION'
  
  ! if (associated(this%buffer_erosion_model)) then
  !   option%io_buffer = 'There may only be one instance of '// &
  !                      'BUFFER_EROSION_MODEL in BUFFER_EROSION "'// &
  !                      trim(this%name)//'".'
  !   call PrintErrMsg(option)
  ! endif
  
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      !------------------------------------------
      case('BUFFER_EROSION_MODEL')
        ! Placeholder for erosion models
      !------------------------------------------
      case('TEST')
        this%test = PETSC_TRUE
      !------------------------------------------
      case default
        call InputKeywordUnrecognized(input,keyword,'BUFFER_EROSION',option)
    end select
  enddo
  call InputPopBlock(input,option)

  verify_string = 'BUFFER_EROSION(' // trim(this%name) // '),'

  ! if (associated(this%buffer_erosion_model)) then
  !   call this%buffer_erosion_model%Verify(verify_string,option)
  ! else
  !   option%io_buffer = 'A buffer erosion model has &
  !        &not been set under BUFFER_EROSION "' // &
  !        trim(this%name) // '". A BUFFER_EROSION_MODEL &
  !        &block must be specified.'
  ! endif

end subroutine BufferErosionRead

! ************************************************************************** !

subroutine ILTRead(illitization_function,input,option)
  !
  ! Reads in contents of a ILLITIZATION_FUNCTION block
  !
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(illitization_base_type) :: illitization_function
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string

  input%ierr = 0
  error_string = 'ILLITIZATION_FUNCTION,'
  select type(ilf => illitization_function)
    class is(ILT_default_type)
      error_string = trim(error_string) // 'DEFAULT'
    class is(ILT_general_type)
      error_string = trim(error_string) // 'GENERAL'
  end select
  
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select type(ilf => illitization_function)
      !------------------------------------------
      class is(ILT_default_type)
        select case(trim(keyword))
          case default
            call ILTDefaultRead(ilf,input,keyword,error_string,'DEFAULT',option)
        end select
      !------------------------------------------
      class is(ILT_general_type)
        select case(trim(keyword))
          case('K_EXP')
            ! Exponent of potassium cation concentration
            call InputReadDouble(input,option,ilf%ilt_K_exp)
            call InputErrorMsg(input,option,'potassium concentration exponent',&
                               'ILLITIZATION, GENERAL')
          case('SMECTITE_EXP')
            ! Exponent of smectite fraction
            call InputReadDouble(input,option,ilf%ilt_exp)
            call InputErrorMsg(input,option,'smectite exponent', &
                               'ILLITIZATION, GENERAL')
          case default
            call ILTDefaultRead(ilf,input,keyword,error_string,'GENERAL',option)
        end select
      !------------------------------------------
      class default
        option%io_buffer = 'Read routine not implemented for ' &
             // trim(error_string) // '.'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine ILTRead

! ************************************************************************** !

subroutine ILTBaseRead(ilf,input,keyword,error_string,kind,option)
  !
  ! Reads in contents of ILLITIZATION_FUNCTION block for the illitization 
  !   base class
  !
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  use Option_module
  use Input_Aux_module
  use String_module

  class(illitization_base_type) :: ilf
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH)   :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=*)  :: kind
  type(option_type) :: option
  
  PetscInt :: i, j
  PetscReal :: a, b, v1, v2, v3, r1, r2
  PetscInt, parameter :: MAX_KD_SIZE = 100
  character(len=MAXWORDLENGTH) :: word
  class(ilt_kd_effects_type), pointer :: shift_kd_list
  class(ilt_perm_effects_type), pointer :: shift_perm
  PetscReal :: f_kd(MAX_KD_SIZE,10), f_perm(10)
  PetscInt :: f_kd_mode_size(MAX_KD_SIZE), f_perm_mode_size
  character(len=MAXWORDLENGTH) :: f_kd_element(MAX_KD_SIZE)
  character(len=MAXWORDLENGTH) :: f_kd_mode(MAX_KD_SIZE), f_perm_mode

  select case(keyword)
    case('SMECTITE_INITIAL')
      ! Initial fraction of smectite in the smectite/illite mixture
      call InputReadDouble(input,option,ilf%ilt_fs0)
      call InputErrorMsg(input,option,'initial smectite fraction', &
                         'ILLITIZATION, '//trim(kind)//'')
    case('THRESHOLD_TEMPERATURE')
      ! Specifies the temperature threshold for activating illitization
      call InputReadDouble(input,option, &
                           ilf%ilt_threshold)
      call InputErrorMsg(input,option,'temperature threshold', &
                         'ILLITIZATION, '//trim(kind)//'')
      call InputReadAndConvertUnits(input,ilf%ilt_threshold,'C', &
                                    'ILLITIZATION, '//trim(kind)// &
                                    ', temperature threshold',option)
    case('SHIFT_PERM')
      ! Functions and parameters modifying the permeability using elements
      !   from the illitization model
      shift_perm => ILTPermEffectsCreate()
      f_perm_mode = ''
      f_perm_mode_size = 0
      f_perm(:) = UNINITIALIZED_DOUBLE
      
      ! Function type
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'f_perm function type', &
                         error_string)
      f_perm_mode = word
      
      ! Function parameters
      select case(f_perm_mode)
        case ('DEFAULT','LINEAR')
          f_perm_mode_size = 1
          call InputReadDouble(input,option,f_perm(1))
          call InputErrorMsg(input,option,'f_perm(1), DEFAULT/LINEAR', &
                             error_string)
          ! Check user values
          if (f_perm(1) < -1.0d+0) then
            option%io_buffer = 'Function parameter #1 in "' &
                             // trim(f_perm_mode) &
                             //'" must not be less than -1 for ' &
                             //'SHIFT_PERM in ILLITIZATION, '//trim(kind)//'.'
            call PrintErrMsg(option)
          endif
        case ('QUADRATIC')
          f_perm_mode_size = 2
          call InputReadDouble(input,option,f_perm(1))
          call InputErrorMsg(input,option,'f_perm(1), QUADRATIC',error_string)
          call InputReadDouble(input,option,f_perm(2))
          call InputErrorMsg(input,option,'f_perm(2), QUADRATIC',error_string)
          ! Check user values
          a = f_perm(1)
          b = f_perm(2)
          v1 = 1 + a + b    ! value at x = 1
          v2 = a + 2*b      ! slope at x = 1
          v3 = a**2 - 4*a*b ! check for real roots
          if (v1 < 0.d0) then ! negative values
            option%io_buffer = 'Function parameters in "' &
                             // trim(f_perm_mode) //'" cannot result ' &
                             //'in a negative value at 100% illite for ' &
                             //'SHIFT_PERM in ILLITIZATION, '//trim(kind)//'.'
            call PrintErrMsg(option)
          endif
          if (a > 0.d0 .and. v2 < 0.d0) then ! positive monotonic
            option%io_buffer = 'Function parameters in "'//trim(f_perm_mode) &
                             //'" must provide monotonic results for ' &
                             //'SHIFT_PERM in ILLITIZATION, '//trim(kind)//'.'
            call PrintErrMsg(option)
          endif
          if (a < 0.d0 .and. v2 > 0.d0) then ! negative monotonic
            option%io_buffer = 'Function parameters in "'//trim(f_perm_mode) &
                             //'" must provide monotonic results for ' &
                             //'SHIFT_PERM in ILLITIZATION, '//trim(kind)//'.'
            call PrintErrMsg(option)
          endif
          if (v3 > 0) then ! roots between 0 and 1
            r1 = (-1*a - sqrt(v3))/(2*b)
            r2 = (-1*a + sqrt(v3))/(2*b)
            if ((r1 >= 0.d0 .and. r1 < 1.d0) .or. & 
                (r2 >= 0.d0 .and. r2 < 1.d0)) then
              option%io_buffer = 'Function parameters in "'//trim(f_perm_mode) &
                               //'" must not have roots between 0 and 1 for ' &
                               //'SHIFT_PERM in ILLITIZATION, '//trim(kind)//'.'
              call PrintErrMsg(option)
            endif
          endif
        case ('POWER')
          f_perm_mode_size = 2
          call InputReadDouble(input,option,f_perm(1))
          call InputErrorMsg(input,option,'f_perm(1), POWER',error_string)
          call InputReadDouble(input,option,f_perm(2))
          call InputErrorMsg(input,option,'f_perm(2), POWER',error_string)
          ! Check user values
          a = f_perm(1)
          b = f_perm(2)
          v1 = 1 + a ! value at x = 1
          if (v1 < 0.d0) then
            option%io_buffer = 'Function parameters in "' &
                             // trim(f_perm_mode) //'" cannot result ' &
                             //'in a negative value at 100% illite for ' &
                             //'SHIFT_PERM in ILLITIZATION, '//trim(kind)//'.'
            call PrintErrMsg(option)
          endif
          if (b <= 0.d0) then
            option%io_buffer = 'Function parameter #2 in "' &
                             // trim(f_perm_mode) //'" must be greater than ' &
                             //'zero for SHIFT_PERM in ILLITIZATION, ' &
                             //trim(kind)//'.'
            call PrintErrMsg(option)
          endif
        case ('EXPONENTIAL')
          f_perm_mode_size = 1
          call InputReadDouble(input,option,f_perm(1))
          call InputErrorMsg(input,option,'f_perm(1), EXPONENTIAL',error_string)
        case default
          option%io_buffer = 'Permeability modification function "' &
                           // trim(f_perm_mode) &
                           //'" was not found among the available options ' &
                           //'for SHIFT_PERM in ILLITIZATION, '//trim(kind)//'.'
          call PrintErrMsg(option)
      end select
      
      if (f_perm_mode_size == 0) then
        option%io_buffer = 'No function parameters were specified &
          &for SHIFT_PERM in ' // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif
      
      allocate(shift_perm%f_perm(f_perm_mode_size))
      shift_perm%f_perm = f_perm(1:f_perm_mode_size)
      allocate(shift_perm%f_perm_mode)
      shift_perm%f_perm_mode = f_perm_mode
      
      ilf%ilt_shift_perm => shift_perm
      
      nullify(shift_perm)
      
    case('SHIFT_KD')
      ! Functions and parameters modifying selected kd values using elements
      !   from the illitization model
      shift_kd_list => ILTKdEffectsCreate()
      i = 0
      f_kd_mode_size(:) = 0
      f_kd(:,:) = UNINITIALIZED_DOUBLE
      f_kd_mode(:) = ''
      f_kd_element(:) = ''
      
      call InputPushBlock(input,option)
      
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        i = i + 1
        if (i > MAX_KD_SIZE) then
          write(word,*) i-1
          option%io_buffer = 'f_kd array in ILLITIZATION must be' &
            //'allocated larger than ' // trim(adjustl(word)) &
            //' under SHIFT_KD in' // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        
        ! Element
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'f_kd element symbol', &
                           error_string)
        f_kd_element(i) = word
        
        ! Function type
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'f_kd function type', &
                           error_string)
        f_kd_mode(i) = word
        
        ! Function parameters
        select case(f_kd_mode(i))
          case ('DEFAULT','LINEAR')
            f_kd_mode_size(i) = 1
            call InputReadDouble(input,option,f_kd(i,1))
            call InputErrorMsg(input,option,'f_kd(*,1), DEFAULT/LINEAR', &
                               error_string)
            
            ! Check user values
            if (f_kd(i,1) < -1.0d+0) then
              option%io_buffer = 'Function parameter #1 in "' &
                               // trim(f_kd_mode(i)) // '" for element "'&
                               // trim(f_kd_element(i)) &
                               //'" must not be less than -1 for ' &
                               //'SHIFT_KD in ILLITIZATION, '//trim(kind)//'.'
              call PrintErrMsg(option)
              
            endif
          case ('QUADRATIC')
            f_kd_mode_size(i) = 2
            call InputReadDouble(input,option,f_kd(i,1))
            call InputErrorMsg(input,option,'f_kd(*,1), QUADRATIC',error_string)
            call InputReadDouble(input,option,f_kd(i,2))
            call InputErrorMsg(input,option,'f_kd(*,2), QUADRATIC',error_string)
            ! Check user values
            a = f_kd(i,1)
            b = f_kd(i,2)
            v1 = 1 + a + b    ! value at x = 1
            v2 = a + 2*b      ! slope at x = 1
            v3 = a**2 - 4*a*b ! check for real roots
            if (v1 < 0.d0) then ! negative values
              option%io_buffer = 'Function parameters in "' &
                               // trim(f_kd_mode(i)) //'" cannot result ' &
                               //'in a negative value at 100% illite for ' &
                               //'SHIFT_KD in ILLITIZATION, '//trim(kind)//'.'
              call PrintErrMsg(option)
            endif
            if (a > 0.d0 .and. v2 < 0.d0) then ! positive monotonic
              option%io_buffer = 'Function parameters in "' &
                               // trim(f_kd_mode(i)) &
                               //'" must provide monotonic results for ' &
                               //'SHIFT_KD in ILLITIZATION, '//trim(kind)//'.'
              call PrintErrMsg(option)
            endif
            if (a < 0.d0 .and. v2 > 0.d0) then ! negative monotonic
              option%io_buffer = 'Function parameters in "' &
                               // trim(f_kd_mode(i)) &
                               //'" must provide monotonic results for ' &
                               //'SHIFT_KD in ILLITIZATION, '//trim(kind)//'.'
              call PrintErrMsg(option)
            endif
            if (v3 > 0) then ! roots between 0 and 1
              r1 = (-1*a - sqrt(v3))/(2*b)
              r2 = (-1*a + sqrt(v3))/(2*b)
              if ((r1 >= 0.d0 .and. r1 < 1.d0) .or. & 
                  (r2 >= 0.d0 .and. r2 < 1.d0)) then
                option%io_buffer = 'Function parameters in "' &
                                 // trim(f_kd_mode(i)) &
                                 //'" must not have roots between 0 and 1 for '&
                                 //'SHIFT_KD in ILLITIZATION, '//trim(kind)//'.'
                call PrintErrMsg(option)
              endif
            endif
          case ('POWER')
            f_kd_mode_size(i) = 2
            call InputReadDouble(input,option,f_kd(i,1))
            call InputErrorMsg(input,option,'f_kd(*,1), POWER',error_string)
            call InputReadDouble(input,option,f_kd(i,2))
            call InputErrorMsg(input,option,'f_kd(*,2), POWER',error_string)
            ! Check user values
            a = f_kd(i,1)
            b = f_kd(i,2)
            v1 = 1 + a ! value at x = 1
            if (v1 < 0.d0) then
              option%io_buffer = 'Function parameters in "' &
                               // trim(f_kd_mode(i)) //'" cannot result ' &
                               //'in a negative value at 100% illite for ' &
                               //'SHIFT_KD in ILLITIZATION, '//trim(kind)//'.'
              call PrintErrMsg(option)
            endif
            if (b <= 0.d0) then
              option%io_buffer = 'Function parameter #2 in "' &
                               // trim(f_kd_mode(i)) //'" must be greater ' &
                               //'than zero for SHIFT_KD in ILLITIZATION, ' &
                               //trim(kind)//'.'
              call PrintErrMsg(option)
            endif
          case ('EXPONENTIAL')
            f_kd_mode_size(i) = 1
            call InputReadDouble(input,option,f_kd(i,1))
            call InputErrorMsg(input,option,'f_kd(*,1), EXPONENTIAL', &
                               error_string)
          case default
            option%io_buffer = 'Sorption modification function "' &
                             // trim(f_kd_mode(i)) // '" for element "'&
                             // trim(f_kd_element(i)) &
                             //'" was not found among the available options ' &
                             //'for SHIFT_KD in ILLITIZATION, '//trim(kind)//'.'
            call PrintErrMsg(option)
        end select
      enddo
      
      call InputPopBlock(input,option)
      
      if (i == 0) then
        option%io_buffer = 'No element/function parameter combinations &
          &specified under SHIFT_KD in ' // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif
      
      j = maxval(f_kd_mode_size)
      
      if (j == 0) then
        option%io_buffer = 'No function parameters were &
          &specified under SHIFT_KD in ' // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif
      
      allocate(shift_kd_list%f_kd(i,j))
      shift_kd_list%f_kd = f_kd(1:i,1:j)
      allocate(shift_kd_list%f_kd_element(i))
      shift_kd_list%f_kd_element = f_kd_element(1:i)
      allocate(shift_kd_list%f_kd_mode(i))
      shift_kd_list%f_kd_mode = f_kd_mode(1:i)
      shift_kd_list%num_elements = i
      
      ilf%ilt_shift_kd_list => shift_kd_list
      
      nullify(shift_kd_list)
      
    case default
      call InputKeywordUnrecognized(input,keyword, &
           'illitization function ('//trim(kind)//')',option)
  end select

end subroutine ILTBaseRead

! ************************************************************************** !

subroutine ILTDefaultRead(ilf,input,keyword,error_string,kind,option)
  !
  ! Reads in contents of ILLITIZATION_FUNCTION block for illitization
  !   default class
  !
  ! Author: Alex Salazar III
  ! Date: 10/12/2021
  !
  use Option_module
  use Input_Aux_module
  use String_module

  class(ILT_default_type) :: ilf
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH)   :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=*)  :: kind
  type(option_type) :: option
  
  select case(keyword)
    case('EA')
      ! Activation energy in Arrhenius term
      call InputReadDouble(input,option,ilf%ilt_ea)
      call InputErrorMsg(input,option,'activation energy', &
                         'ILLITIZATION, '//trim(kind)//'')
      call InputReadAndConvertUnits(input,ilf%ilt_ea, &
                                    'J/mol','ILLITIZATION, '//trim(kind)// &
                                    ', activation energy',option)
    case('FREQ')
      ! Frequency factor (scaling constant of Arrhenius term)
      call InputReadDouble(input,option,ilf%ilt_freq)
      call InputErrorMsg(input,option,'frequency term', &
                         'ILLITIZATION, '//trim(kind)//'')
      call InputReadAndConvertUnits(input,ilf%ilt_freq, &
                                    'L/s-mol','ILLITIZATION, '//trim(kind)// &
                                    ', frequency term',option)
    case('K_CONC')
      ! Concentration of potassium cation
      call InputReadDouble(input,option,ilf%ilt_K_conc)
      call InputErrorMsg(input,option,'potassium concentration', &
                         'ILLITIZATION, '//trim(kind)//'')
      call InputReadAndConvertUnits(input,ilf%ilt_K_conc,'M',&
                                    'ILLITIZATION, ' //trim(kind)// &
                                    ', potassium concentration',option)
    case default
      call ILTBaseRead(ilf,input,keyword,error_string,kind,option)
  end select

end subroutine ILTDefaultRead

! ************************************************************************** !

subroutine MaterialTransformAddToList(new_mtf,list)
  !
  ! Populates the next pointer with a new material transform
  !
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  implicit none

  type(material_transform_type), pointer :: new_mtf
  type(material_transform_type), pointer :: list

  class(material_transform_type), pointer :: cur_mtf

  if (associated(list)) then
    cur_mtf => list
    ! loop to end of list
    do
      if (.not.associated(cur_mtf%next)) exit
      cur_mtf => cur_mtf%next
    enddo
    cur_mtf%next => new_mtf
  else
    list => new_mtf
  endif

end subroutine MaterialTransformAddToList

! ************************************************************************** !

subroutine MaterialTransformConvertListToArray(list,array,option)
  !
  ! Populates the material transform pointer type
  !
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  use String_module
  use Option_module

  implicit none

  class(material_transform_type), pointer :: list
  type(material_transform_ptr_type), pointer :: array(:)
  type(option_type) :: option

  class(material_transform_type), pointer :: cur_mtf
  PetscInt :: count

  count = 0
  cur_mtf => list
  do
    if (.not.associated(cur_mtf)) exit
    count = count + 1
    cur_mtf => cur_mtf%next
  enddo

  if (associated(array)) deallocate(array)
  allocate(array(count))

  count = 0
  cur_mtf => list
  do
    if (.not.associated(cur_mtf)) exit
    count = count + 1
    array(count)%ptr => cur_mtf
    call OptionSetBlocking(option,PETSC_FALSE)
    if (OptionIsIORank(option)) then
      if (associated(cur_mtf%illitization%illitization_function)) then
        if (cur_mtf%illitization%test) then
          call cur_mtf%illitization%illitization_function%Test( &
            cur_mtf%illitization%name,option)
        endif
      endif
      ! if (associated(cur_mtf%buffer_erosion%buffer_erosion_model)) then
      !   if (cur_mtf%buffer_erosion%test) then
      !     call cur_mtf%buffer_erosion%buffer_erosion_model%Test( &
      !       cur_mtf%buffer_erosion_model%name,option)
      !   endif
      ! endif
    endif
    call OptionSetBlocking(option,PETSC_TRUE)
    call OptionCheckNonBlockingError(option)
    cur_mtf => cur_mtf%next
  enddo

end subroutine MaterialTransformConvertListToArray

! ************************************************************************** !

function MaterialTransformGetID(material_transform_array, &
           material_transform_name, material_property_name, option)
  !
  ! Obtains the id number of the material transform from the list
  !
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  use Option_module
  use String_module

  type(material_transform_ptr_type), pointer :: material_transform_array(:)
  character(len=MAXWORDLENGTH) :: material_transform_name
  character(len=MAXWORDLENGTH) :: test1, test2
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: iid, MaterialTransformGetID
  PetscInt :: i, j

  do i = 1, size(material_transform_array)
      test1 = material_transform_array(i)%ptr%name
      do j = 1, size(material_transform_array)
        if (i == j) cycle
        test2 = material_transform_array(j)%ptr%name
        if (test1 == test2) then
          option%io_buffer = 'Duplicate material transform function '//&
                             trim(test2)//&
                             ' has been detected.'
          call PrintErrMsg(option)
        endif
      enddo
  enddo

  MaterialTransformGetID = 0
  do iid = 1, size(material_transform_array)
    if (StringCompare(material_transform_name, &
                      material_transform_array(iid)%ptr%name)) then
      MaterialTransformGetID = iid
      return
    endif
  enddo

  ! MaterialTransformGetID = UNINITIALIZED_INTEGER
  option%io_buffer = 'Material transform function "' // &
                     trim(material_transform_name) // &
                     '" specified in material property "' // &
                     trim(material_property_name) // &
                     '" not found among available functions.'
  call PrintErrMsg(option)

end function MaterialTransformGetID

! ************************************************************************** !

function MaterialTransformCheckILT(material_transform_array,id)
  !
  ! A logical check to determine whether an illitization function is extended
  !
  ! Author: Alex Salazar III
  ! Date: 11/10/2021
  !
  type(material_transform_ptr_type), pointer :: material_transform_array(:)
  PetscInt, intent(in) :: id

  PetscBool :: MaterialTransformCheckILT
  type(illitization_base_type), pointer :: ilt
  
  MaterialTransformCheckILT = PETSC_FALSE
  
  if (associated(material_transform_array(id)%ptr%illitization &
        %illitization_function)) then
    select type(ilt => material_transform_array(id)%ptr%illitization &
                  %illitization_function)
      ! Type must be extended
      class is(ILT_default_type)
        MaterialTransformCheckILT = PETSC_TRUE
    end select
  endif

end function MaterialTransformCheckILT

! ************************************************************************** !

function MaterialTransformCheckBE(material_transform_array,id)
  !
  ! A logical check to determine whether a buffer erosion model is extended
  !
  ! Author: Alex Salazar III
  ! Date: 03/03/2022
  !
  type(material_transform_ptr_type), pointer :: material_transform_array(:)
  PetscInt, intent(in) :: id

  PetscBool :: MaterialTransformCheckBE
  ! type(buffer_erosion_base_type), pointer :: bem
  
  MaterialTransformCheckBE = PETSC_FALSE
  
  ! if (associated(material_transform_array(id)%ptr%buffer_erosion &
  !       %buffer_erosion_model)) then
  !   select type(bem => material_transform_array(id)%ptr%buffer_erosion &
  !                 %buffer_erosion_model)
  !     ! Type must be extended
  !     class is(BE_default_type)
  !       MaterialTransformCheckBE = PETSC_TRUE
  !   end select
  ! endif

end function MaterialTransformCheckBE

! ************************************************************************** !

subroutine MaterialTransformInputRecord(material_transform_list)
  !
  ! Adds details on material transform functions to the input record file
  !
  ! Author: Alex Salazar III
  ! Date: 02/26/2021
  !
  implicit none

  class(material_transform_type), pointer :: material_transform_list

  class(material_transform_type), pointer :: cur_mtf
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
  
  class(illitization_base_type), pointer :: ilf
  class(ilt_kd_effects_type), pointer :: kdl
  class(ilt_perm_effects_type), pointer :: perm
  PetscInt :: i, j, k
  PetscBool :: inactive
  
  inactive = PETSC_TRUE
  
  cur_mtf => material_transform_list
  do
    if (.not.associated(cur_mtf)) exit
    if (associated(cur_mtf%illitization%illitization_function)) then
      select type (ilf => cur_mtf%illitization%illitization_function)
        class is (ILT_default_type)
          inactive = PETSC_FALSE
          exit
        end select
    endif
    cur_mtf => cur_mtf%next
  enddo

  if (inactive) return

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
       &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'MATERIAL TRANSFORM FUNCTIONS'

  cur_mtf => material_transform_list
  do
    if (.not.associated(cur_mtf)) exit
    
    if (associated(cur_mtf%illitization%illitization_function)) then
      select type (ilf => cur_mtf%illitization%illitization_function)
        type is (illitization_base_type)
          exit
        end select
    endif

    write(id,'(a29)',advance='no') 'material transform name: '
    write(id,'(a)') adjustl(trim(cur_mtf%name))
    
    if (associated(cur_mtf%illitization%illitization_function)) then
      write(id,'(a29)') '--------------: '
      write(id,'(a29)',advance='no') 'illitization model: '
      select type (ilf => cur_mtf%illitization%illitization_function)
        !---------------------------------
        class is (ILT_default_type)
          write(id,'(a)') 'Huang et al., 1993'
          write(id,'(a29)',advance='no') 'initial smectite: '
          write(word,'(es12.5)') ilf%ilt_fs0
          write(id,'(a)') adjustl(trim(word))
          write(id,'(a29)',advance='no') 'frequency: '
          write(word,'(es12.5)') ilf%ilt_freq
          write(id,'(a)') adjustl(trim(word))//' L/mol-s'
          write(id,'(a29)',advance='no') 'activation energy: '
          write(word,'(es12.5)') ilf%ilt_ea
          write(id,'(a)') adjustl(trim(word))//' J/mol'
          write(id,'(a29)',advance='no') 'K+ concentration: '
          write(word,'(es12.5)') ilf%ilt_K_conc
          write(id,'(a)') adjustl(trim(word))//' M'
          write(id,'(a29)',advance='no') 'temperature threshold: '
          write(word,'(es12.5)') ilf%ilt_threshold
          write(id,'(a)') adjustl(trim(word))//' C'
          call MaterialTransformPrintPerm(ilf)
          call MaterialTransformPrintKd(ilf)
        !---------------------------------
        class is (ILT_general_type)
          write(id,'(a)') 'Cuadros and Linares, 1996'
          write(id,'(a29)',advance='no') 'initial smectite: '
          write(word,'(es12.5)') ilf%ilt_fs0
          write(id,'(a)') adjustl(trim(word))
          write(id,'(a29)',advance='no') 'smectite exponent: '
          write(word,'(es12.5)') ilf%ilt_exp
          write(id,'(a)') adjustl(trim(word))
          write(id,'(a29)',advance='no') 'frequency: '
          write(word,'(es12.5)') ilf%ilt_freq
          write(id,'(a)') adjustl(trim(word))//' '
          write(id,'(a29)',advance='no') 'activation energy: '
          write(word,'(es12.5)') ilf%ilt_ea
          write(id,'(a)') adjustl(trim(word))//' J/mol'
          write(id,'(a29)',advance='no') 'K+ concentration: '
          write(word,'(es12.5)') ilf%ilt_K_conc
          write(id,'(a)') adjustl(trim(word))//' M'
          write(id,'(a29)',advance='no') 'K+ conc. exponent: '
          write(word,'(es12.5)') ilf%ilt_K_exp
          write(id,'(a)') adjustl(trim(word))
          write(id,'(a29)',advance='no') 'temperature threshold: '
          write(word,'(es12.5)') ilf%ilt_threshold
          write(id,'(a)') adjustl(trim(word))//' C'
          call MaterialTransformPrintPerm(ilf)
          call MaterialTransformPrintKd(ilf)
      end select
    endif

    write(id,'(a29)') '---------------------------: '
    cur_mtf => cur_mtf%next
  enddo

end subroutine MaterialTransformInputRecord

! ************************************************************************** !

subroutine MaterialTransformPrintKd(ilf)
  !
  ! Adds details on kd parameters to the input record file
  !
  ! Author: Alex Salazar III
  ! Date: 11/15/2021
  !
  implicit none

  class(illitization_base_type) :: ilf
  
  class(ilt_kd_effects_type), pointer :: kdl
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
  PetscInt :: i, j, k
  
  if (.not. associated(ilf%ilt_shift_kd_list)) return
  
  j = 0
  kdl => ilf%ilt_shift_kd_list
  
  write(id,'(a29)',advance='no') 'shift (kd): '
  do i = 1, kdl%num_elements
    if (.not. i == 1) then
      write(id,'(a29)',advance='no') ""
    endif
    write(word,'(a)') kdl%f_kd_element(i)
    write(id,'(a)',advance='no') adjustl(trim(word))//" "
    write(word,'(a)') kdl%f_kd_mode(i)
    write(id,'(a)',advance='no') adjustl(trim(word))
    select case(kdl%f_kd_mode(i))
      case ('DEFAULT','LINEAR')
        j = 1
      case ('QUADRATIC')
        j = 2
      case ('POWER')
        j = 2
      case ('EXPONENTIAL')
        j = 1
    end select
    do k = 1, j
      write(word,'(es12.5)') kdl%f_kd(i,k)
      write(id,'(a)',advance='no') " "//adjustl(trim(word))
    enddo
    write(id,'(a)')
  enddo

end subroutine MaterialTransformPrintKd

! ************************************************************************** !

subroutine MaterialTransformPrintPerm(ilf)
  !
  ! Adds details on permeability parameters to the input record file
  !
  ! Author: Alex Salazar III
  ! Date: 11/15/2021
  !
  implicit none

  class(illitization_base_type) :: ilf
  
  class(ilt_perm_effects_type), pointer :: perm
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
  PetscInt :: j, k
  
  if (.not. associated(ilf%ilt_shift_perm)) return
  
  j = 0
  perm => ilf%ilt_shift_perm
  
  write(id,'(a29)',advance='no') 'shift (permeability): '
  perm => ilf%ilt_shift_perm
  write(word,'(a)') perm%f_perm_mode
  write(id,'(a)',advance='no') adjustl(trim(word))
  select case(perm%f_perm_mode)
    case ('DEFAULT','LINEAR')
      j = 1
    case ('QUADRATIC')
      j = 2
    case ('POWER')
      j = 2
    case ('EXPONENTIAL')
      j = 1
  end select
  do k = 1, j
    write(word,'(es12.5)') perm%f_perm(k)
    write(id,'(a)',advance='no') " "//adjustl(trim(word))
  enddo
  write(id,'(a)')

end subroutine MaterialTransformPrintPerm

! ************************************************************************** !

subroutine MaterialTransformAuxVarInit(auxvar)
  !
  ! Initializes a material transform auxiliary object
  !
  ! Author: Alex Salazar III
  ! Date: 02/10/2022
  !

  implicit none

  type(material_transform_auxvar_type) :: auxvar

  nullify(auxvar%il_aux)
  nullify(auxvar%be_aux)


end subroutine MaterialTransformAuxVarInit

! ************************************************************************** !

function IllitizationAuxVarInit(option)
  !
  ! Initializes an illitization auxiliary object
  !
  ! Author: Alex Salazar III
  ! Date: 02/10/2022
  !

  use Option_module

  implicit none

  class(illitization_auxvar_type), pointer :: IllitizationAuxVarInit
  class(illitization_auxvar_type), pointer :: auxvar
  type(option_type) :: option

  allocate(auxvar)
  ! auxvar%il_aux%initial_pressure = UNINITIALIZED_DOUBLE
  auxvar%fs0    = 1.0d+0               ! initial fraction of smectite in material
  auxvar%fs     = UNINITIALIZED_DOUBLE ! fraction of smectite in material
  auxvar%fi     = UNINITIALIZED_DOUBLE ! fraction of illite in material
  auxvar%ts     = UNINITIALIZED_DOUBLE ! track time of last change in smectite
  auxvar%scale  = UNINITIALIZED_DOUBLE ! scale factor
  auxvar%qperm0 = PETSC_FALSE          ! save initial permeability

  if (option%iflowmode /= NULL_MODE) then
    if (option%flow%full_perm_tensor) then
      allocate(auxvar%perm0(6))
    else
      allocate(auxvar%perm0(3))
    endif
    auxvar%perm0 = UNINITIALIZED_DOUBLE
  else
    ! nullify(auxvar%perm0)
  endif

  IllitizationAuxVarInit => auxvar

end function IllitizationAuxVarInit

! ************************************************************************** !

subroutine IllitizationAuxVarStrip(auxvar)
  !
  ! Deallocates an illitization auxiliary object
  !
  ! Author: Alex Salazar III
  ! Date: 02/10/2022
  !

  use Utility_module, only : DeallocateArray

  implicit none

  class(illitization_auxvar_type), pointer :: auxvar

  if (.not.associated(auxvar)) return
  
  ! call DeallocateArray(auxvar%perm0)

  deallocate(auxvar)
  nullify(auxvar)

end subroutine IllitizationAuxVarStrip

! ************************************************************************** !

function BufferErosionAuxVarInit()
  !
  ! Initializes a buffer erosion auxiliary object
  !
  ! Author: Alex Salazar III
  ! Date: 02/10/2022
  !

  implicit none

  class(buffer_erosion_auxvar_type), pointer :: BufferErosionAuxVarInit
  class(buffer_erosion_auxvar_type), pointer :: auxvar

  allocate(auxvar)

  BufferErosionAuxVarInit => auxvar

end function BufferErosionAuxVarInit

! ************************************************************************** !

subroutine BufferErosionAuxVarStrip(auxvar)
  !
  ! Deallocates a buffer erosion auxiliary object
  !
  ! Author: Alex Salazar III
  ! Date: 02/10/2022
  !

  implicit none

  class(buffer_erosion_auxvar_type), pointer :: auxvar

  if (.not.associated(auxvar)) return

  deallocate(auxvar)
  nullify(auxvar)

end subroutine BufferErosionAuxVarStrip

! ************************************************************************** !

subroutine MaterialTransformAuxVarStrip(auxvar)
  ! 
  ! Deallocates a material transform auxiliary object
  ! 
  ! Author: Alex Salazar
  ! Date: 01/20/2022
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  type(material_transform_auxvar_type) :: auxvar
  
  if (associated(auxvar%il_aux)) then
    call IllitizationAuxVarStrip(auxvar%il_aux)
  endif
  if (associated(auxvar%be_aux)) then
    call BufferErosionAuxVarStrip(auxvar%be_aux)
  endif
  
  ! call DeallocateArray(auxvar%variable)  
  
end subroutine MaterialTransformAuxVarStrip

! ************************************************************************** !

recursive subroutine MaterialTransformDestroy(mtf)

  implicit none

  type(material_transform_type), pointer :: mtf
  
  PetscInt :: i

  if (.not. associated(mtf)) return

  call MaterialTransformDestroy(mtf%next)

  if (associated(mtf%auxvars)) then
    do i = 1, size(mtf%auxvars)
      call MaterialTransformAuxVarStrip(mtf%auxvars(i))
    enddo
    deallocate(mtf%auxvars)
  endif
  nullify(mtf%auxvars)

  if (associated(mtf%illitization)) then
    call IllitizationDestroy(mtf%illitization)
  endif
  
  if (associated(mtf%buffer_erosion)) then
    call BufferErosionDestroy(mtf%buffer_erosion)
  endif

  nullify(mtf)

end subroutine MaterialTransformDestroy

! ************************************************************************** !

recursive subroutine IllitizationDestroy(illitization)

  implicit none

  class(illitization_type), pointer :: illitization

  if (.not. associated(illitization)) return

  if (associated(illitization%illitization_function)) then
    call ILTDestroy(illitization%illitization_function)
  endif

  deallocate(illitization)
  nullify(illitization)

end subroutine IllitizationDestroy

! ************************************************************************** !

recursive subroutine BufferErosionDestroy(buffer_erosion)

  implicit none

  class(buffer_erosion_type), pointer :: buffer_erosion

  if (.not. associated(buffer_erosion)) return

  deallocate(buffer_erosion)
  nullify(buffer_erosion)

end subroutine BufferErosionDestroy

end module Material_Transform_module
