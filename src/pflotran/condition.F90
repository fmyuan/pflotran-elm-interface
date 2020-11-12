module Condition_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Global_Aux_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Time_Storage_module
  use Lookup_Table_module
  use Transport_Constraint_Base_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: flow_condition_type
    PetscInt :: id                          ! id from which condition can be referenced
    PetscBool :: is_transient
    PetscBool :: sync_time_with_update
    character(len=MAXWORDLENGTH) :: name    ! name of condition (e.g. initial, recharge)
    PetscInt :: num_sub_conditions
    PetscInt :: iphase
    PetscInt, pointer :: itype(:)
    character(len=MAXWORDLENGTH) :: time_units
    character(len=MAXWORDLENGTH) :: length_units
    type(time_storage_type), pointer :: default_time_storage
    class(dataset_base_type), pointer :: datum
    type(flow_sub_condition_type), pointer :: datum_z
    type(flow_sub_condition_type), pointer :: pressure
    type(flow_sub_condition_type), pointer :: saturation
    type(flow_sub_condition_type), pointer :: rate
    type(flow_sub_condition_type), pointer :: well
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: concentration
    type(flow_sub_condition_type), pointer :: enthalpy
    type(flow_sub_condition_type), pointer :: energy_rate
    type(flow_sub_condition_type), pointer :: energy_flux
    type(flow_general_condition_type), pointer :: general
    type(flow_hydrate_condition_type), pointer :: hydrate
    type(flow_toil_ims_condition_type), pointer :: toil_ims
    type(flow_towg_condition_type), pointer :: towg
    class(lookup_table_general_type), pointer :: rtempvz_table  !temperature variation over z
    ! any new sub conditions must be added to FlowConditionIsTransient
    type(sub_condition_ptr_type), pointer :: sub_condition_ptr(:)
    type(flow_condition_type), pointer :: next ! pointer to next condition_type for linked-lists
  end type flow_condition_type

  ! data structure for general phase
  type, public :: flow_general_condition_type
    type(flow_sub_condition_type), pointer :: liquid_pressure
    type(flow_sub_condition_type), pointer :: gas_pressure
    type(flow_sub_condition_type), pointer :: gas_saturation
    type(flow_sub_condition_type), pointer :: mole_fraction
    type(flow_sub_condition_type), pointer :: relative_humidity
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: rate
    type(flow_sub_condition_type), pointer :: liquid_flux
    type(flow_sub_condition_type), pointer :: gas_flux
    type(flow_sub_condition_type), pointer :: energy_flux
    ! any new sub conditions must be added to FlowConditionIsTransient
  end type flow_general_condition_type

  ! data structure for hydrate
  type, public :: flow_hydrate_condition_type
    type(flow_sub_condition_type), pointer :: liquid_pressure
    type(flow_sub_condition_type), pointer :: gas_pressure
    type(flow_sub_condition_type), pointer :: liquid_saturation
    type(flow_sub_condition_type), pointer :: gas_saturation
    type(flow_sub_condition_type), pointer :: hydrate_saturation
    type(flow_sub_condition_type), pointer :: ice_saturation
    type(flow_sub_condition_type), pointer :: mole_fraction
    type(flow_sub_condition_type), pointer :: relative_humidity
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: rate
    type(flow_sub_condition_type), pointer :: liquid_flux
    type(flow_sub_condition_type), pointer :: gas_flux
    type(flow_sub_condition_type), pointer :: energy_flux
    ! any new sub conditions must be added to FlowConditionIsTransient
  end type flow_hydrate_condition_type

  ! data structure for toil_ims
  type, public :: flow_toil_ims_condition_type
    type(flow_sub_condition_type), pointer :: pressure
    type(flow_sub_condition_type), pointer :: saturation
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: enthalpy
    type(flow_sub_condition_type), pointer :: rate
    type(flow_sub_condition_type), pointer :: liquid_flux
    type(flow_sub_condition_type), pointer :: oil_flux
    type(flow_sub_condition_type), pointer :: energy_flux
    type(flow_sub_condition_type), pointer :: owc   ! oil water contact
    type(flow_sub_condition_type), pointer :: owc_z
    type(flow_sub_condition_type), pointer :: pcow_owc
    type(flow_sub_condition_type), pointer :: liq_press_grad ! water piezometric head gradient
    ! any new sub conditions must be added to FlowConditionIsTransient
  end type flow_toil_ims_condition_type

  ! data structure for towg
  ! some of these variables depend on primary variable choice
  type, public :: flow_towg_condition_type
    PetscBool :: is_wg_equilibration
    type(flow_sub_condition_type), pointer :: oil_pressure
    type(flow_sub_condition_type), pointer :: gas_pressure
    type(flow_sub_condition_type), pointer :: oil_saturation
    type(flow_sub_condition_type), pointer :: gas_saturation
    type(flow_sub_condition_type), pointer :: solvent_saturation
    type(flow_sub_condition_type), pointer :: bubble_point
    type(flow_sub_condition_type), pointer :: liquid_flux
    type(flow_sub_condition_type), pointer :: oil_flux
    type(flow_sub_condition_type), pointer :: gas_flux
    type(flow_sub_condition_type), pointer :: solvent_flux
    type(flow_sub_condition_type), pointer :: energy_flux
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: enthalpy
    type(flow_sub_condition_type), pointer :: rate
    type(flow_sub_condition_type), pointer :: bhp_pressure
    type(flow_sub_condition_type), pointer :: owc_z
    type(flow_sub_condition_type), pointer :: pcow_owc
    type(flow_sub_condition_type), pointer :: ogc_z
    type(flow_sub_condition_type), pointer :: pcog_ogc
    class(lookup_table_general_type), pointer :: pbvz_table
    ! any new sub conditions must be added to FlowConditionIsTransient
  end type flow_towg_condition_type

  type, public :: flow_sub_condition_type
    PetscInt :: itype                  ! integer describing type of condition
    PetscInt :: isubtype
    character(len=MAXWORDLENGTH) :: ctype ! character string describing type of condition
    character(len=MAXWORDLENGTH) :: units      ! units
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: aux_real(2)
    class(dataset_base_type), pointer :: gradient
    class(dataset_base_type), pointer :: dataset
  end type flow_sub_condition_type

  type, public :: sub_condition_ptr_type
    type(flow_sub_condition_type), pointer :: ptr
  end type sub_condition_ptr_type

  type, public :: condition_ptr_type
    type(flow_condition_type), pointer :: ptr
  end type condition_ptr_type

  type, public :: condition_list_type
    PetscInt :: num_conditions
    type(flow_condition_type), pointer :: first
    type(flow_condition_type), pointer :: last
    type(flow_condition_type), pointer :: array(:)
  end type condition_list_type

  type, public :: tran_condition_type
    PetscInt :: id                     ! id from which condition can be referenced
    PetscInt :: itype                  ! integer describing type of condition
    PetscBool :: is_transient
    character(len=MAXWORDLENGTH) :: name  ! name of condition (e.g. initial, recharge)
    class(tran_constraint_coupler_base_type), pointer :: constraint_coupler_list
    class(tran_constraint_coupler_base_type), pointer :: cur_constraint_coupler
    type(tran_condition_type), pointer :: next
  end type tran_condition_type

  type, public :: tran_condition_ptr_type
    type(tran_condition_type), pointer :: ptr
  end type tran_condition_ptr_type

  type, public :: tran_condition_list_type
    PetscInt :: num_conditions
    type(tran_condition_type), pointer :: first
    type(tran_condition_type), pointer :: last
    type(tran_condition_ptr_type), pointer :: array(:)
  end type tran_condition_list_type

  public :: FlowConditionCreate, FlowConditionDestroy, FlowConditionRead, &
            FlowConditionGeneralRead, FlowConditionTOilImsRead, &
            FlowConditionHydrateRead, FlowConditionTOWGRead, &
            FlowConditionAddToList, FlowConditionInitList, &
            FlowConditionDestroyList, &
            FlowConditionGetPtrFromList, FlowConditionUpdate, &
            FlowConditionPrint, &
            TranConditionCreate, &
            TranConditionAddToList, TranConditionInitList, &
            TranConditionDestroyList, TranConditionGetPtrFromList, &
            TranConditionRead, &
            TranConditionUpdate, &
            FlowConditionIsTransient, &
            ConditionReadValues, &
            GetSubConditionName, &
            FlowConditionUnknownItype, &
            FlowCondInputRecord, &
            TranCondInputRecord

contains

! ************************************************************************** !

function FlowConditionCreate(option)
  !
  ! Creates a condition
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(flow_condition_type), pointer :: FlowConditionCreate

  type(flow_condition_type), pointer :: condition

  allocate(condition)
  nullify(condition%pressure)
  nullify(condition%saturation)
  nullify(condition%rate)
  nullify(condition%energy_rate)
  nullify(condition%energy_flux)
  nullify(condition%well)
  nullify(condition%temperature)
  nullify(condition%concentration)
  nullify(condition%enthalpy)
  nullify(condition%sub_condition_ptr)
  nullify(condition%general)
  nullify(condition%hydrate)
  nullify(condition%toil_ims)
  nullify(condition%towg)
  nullify(condition%rtempvz_table)
  nullify(condition%itype)
  nullify(condition%next)
  nullify(condition%datum)
  nullify(condition%datum_z)
  nullify(condition%default_time_storage)
  condition%is_transient = PETSC_FALSE
  condition%sync_time_with_update = PETSC_FALSE
  condition%time_units = ''
  condition%length_units = ''
  condition%id = 0
  condition%iphase = 0
  condition%num_sub_conditions = 0
  condition%name = ''

  FlowConditionCreate => condition

end function FlowConditionCreate

! ************************************************************************** !

function TranConditionCreate(option)
  !
  ! Creates a transport condition
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(tran_condition_type), pointer :: TranConditionCreate

  type(tran_condition_type), pointer :: condition

  allocate(condition)
  nullify(condition%constraint_coupler_list)
  nullify(condition%cur_constraint_coupler)
  nullify(condition%next)
  condition%id = 0
  condition%itype = 0
  condition%name = ''

  TranConditionCreate => condition

end function TranConditionCreate

! ************************************************************************** !

function FlowGeneralConditionCreate(option)
  !
  ! Creates a condition for general mode
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/11
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(flow_general_condition_type), pointer :: FlowGeneralConditionCreate

  type(flow_general_condition_type), pointer :: general_condition

  allocate(general_condition)
  nullify(general_condition%liquid_pressure)
  nullify(general_condition%gas_pressure)
  nullify(general_condition%gas_saturation)
  nullify(general_condition%relative_humidity)
  nullify(general_condition%mole_fraction)
  nullify(general_condition%temperature)
  nullify(general_condition%liquid_flux)
  nullify(general_condition%gas_flux)
  nullify(general_condition%energy_flux)
  nullify(general_condition%rate)

  FlowGeneralConditionCreate => general_condition

end function FlowGeneralConditionCreate

! ************************************************************************** !
function FlowHydrateConditionCreate(option)
  !
  ! Creates a condition for hydrate mode
  !
  ! Author: Michael Nole
  ! Date: 07/22/19
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(flow_hydrate_condition_type), pointer :: FlowHydrateConditionCreate

  type(flow_hydrate_condition_type), pointer :: hydrate_condition

  allocate(hydrate_condition)
  nullify(hydrate_condition%liquid_pressure)
  nullify(hydrate_condition%gas_pressure)
  nullify(hydrate_condition%liquid_saturation)
  nullify(hydrate_condition%gas_saturation)
  nullify(hydrate_condition%hydrate_saturation)
  nullify(hydrate_condition%ice_saturation)
  nullify(hydrate_condition%relative_humidity)
  nullify(hydrate_condition%mole_fraction)
  nullify(hydrate_condition%temperature)
  nullify(hydrate_condition%liquid_flux)
  nullify(hydrate_condition%gas_flux)
  nullify(hydrate_condition%energy_flux)
  nullify(hydrate_condition%rate)

  FlowHydrateConditionCreate => hydrate_condition

end function FlowHydrateConditionCreate

! ************************************************************************** !

function FlowTOilImsConditionCreate(option)
  !
  ! Creates a condition for toil_ims mode
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 09/9/2015
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(flow_toil_ims_condition_type), pointer :: FlowTOilImsConditionCreate

  type(flow_toil_ims_condition_type), pointer :: toil_ims_condition

  allocate(toil_ims_condition)
  nullify(toil_ims_condition%pressure)
  nullify(toil_ims_condition%saturation)
  nullify(toil_ims_condition%temperature)
  nullify(toil_ims_condition%enthalpy)
  nullify(toil_ims_condition%rate)
  nullify(toil_ims_condition%liquid_flux)
  nullify(toil_ims_condition%oil_flux)
  nullify(toil_ims_condition%energy_flux)
  nullify(toil_ims_condition%owc)
  nullify(toil_ims_condition%owc_z)
  nullify(toil_ims_condition%pcow_owc)
  nullify(toil_ims_condition%liq_press_grad)

  FlowTOilImsConditionCreate => toil_ims_condition

end function FlowTOilImsConditionCreate

! ************************************************************************** !

function FlowTOWGConditionCreate(option)
  !
  ! Creates a condition for TOWG mode
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/16/2016
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(flow_towg_condition_type), pointer :: FlowTOWGConditionCreate

  type(flow_towg_condition_type), pointer :: towg_condition

  allocate(towg_condition)

  towg_condition%is_wg_equilibration = PETSC_FALSE
  nullify(towg_condition%oil_pressure)
  nullify(towg_condition%gas_pressure)
  nullify(towg_condition%oil_saturation)
  nullify(towg_condition%gas_saturation)
  nullify(towg_condition%solvent_saturation)
  nullify(towg_condition%bubble_point)
  nullify(towg_condition%liquid_flux)
  nullify(towg_condition%oil_flux)
  nullify(towg_condition%gas_flux)
  nullify(towg_condition%solvent_flux)
  nullify(towg_condition%energy_flux)
  nullify(towg_condition%temperature)
  nullify(towg_condition%enthalpy)
  nullify(towg_condition%rate)
  nullify(towg_condition%bhp_pressure)
  nullify(towg_condition%owc_z)
  nullify(towg_condition%pcow_owc)
  nullify(towg_condition%ogc_z)
  nullify(towg_condition%pcog_ogc)
  nullify(towg_condition%pbvz_table)

  FlowTOWGConditionCreate => towg_condition

end function FlowTOWGConditionCreate

! ************************************************************************** !

function FlowGeneralSubConditionPtr(input,sub_condition_name,general, &
                                    option)
  !
  ! Returns a pointer to a subcondition, creating
  ! them if necessary
  !
  ! Author: Glenn Hammond
  ! Date: 06/09/11
  !

  use Option_module
  use Input_Aux_module

  implicit none

  type(input_type) :: input
  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_general_condition_type) :: general
  type(option_type) :: option

  type(flow_sub_condition_type), pointer :: FlowGeneralSubConditionPtr
  type(flow_sub_condition_type), pointer :: sub_condition_ptr

  select case(sub_condition_name)
    case('LIQUID_PRESSURE')
      if (associated(general%liquid_pressure)) then
        sub_condition_ptr => general%liquid_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%liquid_pressure => sub_condition_ptr
      endif
    case('GAS_PRESSURE')
      if (associated(general%gas_pressure)) then
        sub_condition_ptr => general%gas_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%gas_pressure => sub_condition_ptr
      endif
    case('LIQUID_SATURATION','GAS_SATURATION')
      if (associated(general%gas_saturation)) then
        sub_condition_ptr => general%gas_saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%gas_saturation => sub_condition_ptr
      endif
    case('TEMPERATURE')
      if (associated(general%temperature)) then
        sub_condition_ptr => general%temperature
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%temperature => sub_condition_ptr
      endif
    case('RELATIVE_HUMIDITY')
      if (associated(general%relative_humidity)) then
        sub_condition_ptr => general%relative_humidity
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%relative_humidity => sub_condition_ptr
      endif
    case('MOLE_FRACTION')
      if (associated(general%mole_fraction)) then
        sub_condition_ptr => general%mole_fraction
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%mole_fraction => sub_condition_ptr
      endif
    case('LIQUID_FLUX')
      if (associated(general%liquid_flux)) then
        sub_condition_ptr => general%liquid_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%liquid_flux => sub_condition_ptr
      endif
    case('GAS_FLUX')
      if (associated(general%gas_flux)) then
        sub_condition_ptr => general%gas_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%gas_flux => sub_condition_ptr
      endif
    case('ENERGY_FLUX')
      if (associated(general%energy_flux)) then
        sub_condition_ptr => general%energy_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%energy_flux => sub_condition_ptr
      endif
    case('RATE')
      if (associated(general%rate)) then
        sub_condition_ptr => general%rate
      else
        sub_condition_ptr => FlowSubConditionCreate(option%nflowdof)
        general%rate => sub_condition_ptr
      endif
    case default
      call InputKeywordUnrecognized(input,sub_condition_name, &
                                    'general condition,type',option)
  end select

  FlowGeneralSubConditionPtr => sub_condition_ptr

end function FlowGeneralSubConditionPtr

! ************************************************************************** !

function FlowHydrateSubConditionPtr(input,sub_condition_name,hydrate, &
                                    option)
  !
  ! Returns a pointer to a subcondition, creating
  ! them if necessary
  !
  ! Author: Michael Nole
  ! Date: 07/22/19
  !

  use Option_module
  use Input_Aux_module

  implicit none

  type(input_type) :: input
  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_hydrate_condition_type) :: hydrate
  type(option_type) :: option

  type(flow_sub_condition_type), pointer :: FlowHydrateSubConditionPtr
  type(flow_sub_condition_type), pointer :: sub_condition_ptr

  select case(sub_condition_name)
    case('LIQUID_PRESSURE')
      if (associated(hydrate%liquid_pressure)) then
        sub_condition_ptr => hydrate%liquid_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%liquid_pressure => sub_condition_ptr
      endif
    case('GAS_PRESSURE')
      if (associated(hydrate%gas_pressure)) then
        sub_condition_ptr => hydrate%gas_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%gas_pressure => sub_condition_ptr
      endif
    case('LIQUID_SATURATION','GAS_SATURATION') !MAN should split this
      if (associated(hydrate%gas_saturation)) then
        sub_condition_ptr => hydrate%gas_saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%gas_saturation => sub_condition_ptr
      endif
!    case('LIQUID_SATURATION')
!      if (associated(hydrate%gas_saturation)) then
!         sub_condition_ptr => hydrate%liquid_saturation
!      else
!         sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
!         hydrate%liquid_saturation => sub_condition_ptr
!      endif
    case('HYDRATE_SATURATION')
      if (associated(hydrate%hydrate_saturation)) then
        sub_condition_ptr => hydrate%hydrate_saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%hydrate_saturation => sub_condition_ptr
      endif
    case('ICE_SATURATION')
      if (associated(hydrate%ice_saturation)) then
        sub_condition_ptr => hydrate%ice_saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%ice_saturation => sub_condition_ptr
      endif
    case('TEMPERATURE')
      if (associated(hydrate%temperature)) then
        sub_condition_ptr => hydrate%temperature
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%temperature => sub_condition_ptr
      endif
    case('RELATIVE_HUMIDITY')
      if (associated(hydrate%relative_humidity)) then
        sub_condition_ptr => hydrate%relative_humidity
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%relative_humidity => sub_condition_ptr
      endif
    case('MOLE_FRACTION')
      if (associated(hydrate%mole_fraction)) then
        sub_condition_ptr => hydrate%mole_fraction
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%mole_fraction => sub_condition_ptr
      endif
    case('LIQUID_FLUX')
      if (associated(hydrate%liquid_flux)) then
        sub_condition_ptr => hydrate%liquid_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%liquid_flux => sub_condition_ptr
      endif
    case('GAS_FLUX')
      if (associated(hydrate%gas_flux)) then
        sub_condition_ptr => hydrate%gas_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%gas_flux => sub_condition_ptr
      endif
    case('ENERGY_FLUX')
      if (associated(hydrate%energy_flux)) then
        sub_condition_ptr => hydrate%energy_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        hydrate%energy_flux => sub_condition_ptr
      endif
    case('RATE')
      if (associated(hydrate%rate)) then
        sub_condition_ptr => hydrate%rate
      else
        sub_condition_ptr => FlowSubConditionCreate(option%nflowdof)
        hydrate%rate => sub_condition_ptr
      endif
    case default
      call InputKeywordUnrecognized(input,sub_condition_name, &
                                    'hydrate condition,type',option)
  end select

  FlowHydrateSubConditionPtr => sub_condition_ptr

end function FlowHydrateSubConditionPtr

! ************************************************************************** !

function FlowTOilImsSubConditionPtr(input,sub_condition_name,toil_ims, &
                                    option)
  !
  ! Returns a pointer to a subcondition, creating
  ! them if necessary for toil_ims flow mode
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 09/09/2015
  !

  use Option_module
  use Input_Aux_module

  implicit none

  type(input_type) :: input
  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_toil_ims_condition_type) :: toil_ims
  type(option_type) :: option

  type(flow_sub_condition_type), pointer :: FlowTOilImsSubConditionPtr
  type(flow_sub_condition_type), pointer :: sub_condition_ptr

  select case(sub_condition_name)
    case('PRESSURE','OIL_PRESSURE','WATER_PRESSURE')
      if (associated(toil_ims%pressure)) then
        sub_condition_ptr => toil_ims%pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%pressure => sub_condition_ptr
      endif
    case('LIQUID_SATURATION','OIL_SATURATION')
      if (associated(toil_ims%saturation)) then
        sub_condition_ptr => toil_ims%saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%saturation => sub_condition_ptr
      endif
    case('TEMPERATURE','RTEMP','TEMPERATURE_AT_DATUM')
      if (associated(toil_ims%temperature)) then
        sub_condition_ptr => toil_ims%temperature
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%temperature => sub_condition_ptr
      endif
    case('ENTHALPY')
      if (associated(toil_ims%enthalpy)) then
        sub_condition_ptr => toil_ims%enthalpy
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%enthalpy => sub_condition_ptr
      endif
    case('LIQUID_FLUX')
      if (associated(toil_ims%liquid_flux)) then
        sub_condition_ptr => toil_ims%liquid_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%liquid_flux => sub_condition_ptr
      endif
    case('OIL_FLUX')
      if (associated(toil_ims%oil_flux)) then
        sub_condition_ptr => toil_ims%oil_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%oil_flux => sub_condition_ptr
      endif
    case('ENERGY_FLUX')
      if (associated(toil_ims%energy_flux)) then
        sub_condition_ptr => toil_ims%energy_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%energy_flux => sub_condition_ptr
      endif
    case('OWC')
      if (associated(toil_ims%owc)) then
        sub_condition_ptr => toil_ims%owc
      else
        sub_condition_ptr => FlowSubConditionCreate(THREE_INTEGER)
        toil_ims%owc => sub_condition_ptr
      endif
    case('WATER_PRESSURE_GRAD')
      if (associated(toil_ims%liq_press_grad)) then
        sub_condition_ptr => toil_ims%liq_press_grad
      else
        sub_condition_ptr => FlowSubConditionCreate(THREE_INTEGER)
        toil_ims%liq_press_grad => sub_condition_ptr
      endif
    !a condition can have either RATE or WELL_RATE (target/limits)
    case('RATE')
      if (associated(toil_ims%rate)) then
        sub_condition_ptr => toil_ims%rate
      else
        ! energy rate is loaded in the third record
        sub_condition_ptr => FlowSubConditionCreate(THREE_INTEGER)
        toil_ims%rate => sub_condition_ptr
      endif
    case('OWC_Z','OWC_D')
      if (associated(toil_ims%owc_z)) then
        sub_condition_ptr => toil_ims%owc_z
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%owc_z => sub_condition_ptr
      endif
    case('PCOW_OWC')
      if (associated(toil_ims%pcow_owc)) then
        sub_condition_ptr => toil_ims%pcow_owc
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%pcow_owc => sub_condition_ptr
      endif
    case default
      call InputKeywordUnrecognized(input,sub_condition_name, &
                                    'toil_ims condition,type',option)
  end select

  FlowTOilImsSubConditionPtr => sub_condition_ptr

end function FlowTOilImsSubConditionPtr

! ************************************************************************** !

function FlowTOWGSubConditionPtr(input,sub_condition_name,towg,option)
  !
  ! Returns a pointer to a subcondition, creating
  ! them if necessary for TOWG mode
  !
  ! Author: Paolo Orsini
  ! Date: 10/17/16
  !

  use Option_module
  use Input_Aux_module

  implicit none

  type(input_type) :: input
  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_towg_condition_type) :: towg
  type(option_type) :: option

  type(flow_sub_condition_type), pointer :: FlowTOWGSubConditionPtr
  type(flow_sub_condition_type), pointer :: sub_condition_ptr

  select case(sub_condition_name)
  case('OIL_PRESSURE','PRESSURE')
      if (associated(towg%oil_pressure)) then
        sub_condition_ptr => towg%oil_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%oil_pressure => sub_condition_ptr
      endif
    case('GAS_PRESSURE')
      if (associated(towg%gas_pressure)) then
        sub_condition_ptr => towg%gas_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%gas_pressure => sub_condition_ptr
      endif
    case('OIL_SATURATION')
      if (associated(towg%oil_saturation)) then
        sub_condition_ptr => towg%oil_saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%oil_saturation => sub_condition_ptr
      endif
    case('GAS_SATURATION')
      if (associated(towg%gas_saturation)) then
        sub_condition_ptr => towg%gas_saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%gas_saturation => sub_condition_ptr
      endif
    case('SOLVENT_SATURATION')
      if (associated(towg%solvent_saturation)) then
        sub_condition_ptr => towg%solvent_saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%solvent_saturation => sub_condition_ptr
      endif
    case('BUBBLE_POINT')
      if (associated(towg%bubble_point)) then
        sub_condition_ptr => towg%bubble_point
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%bubble_point => sub_condition_ptr
      endif
    case('TEMPERATURE','RTEMP','TEMPERATURE_AT_DATUM')
      if (associated(towg%temperature)) then
        sub_condition_ptr => towg%temperature
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%temperature => sub_condition_ptr
      endif
    case('ENTHALPY')
      if (associated(towg%enthalpy)) then
        sub_condition_ptr => towg%enthalpy
      else
        sub_condition_ptr => FlowSubConditionCreate(option%nphase)
        towg%enthalpy => sub_condition_ptr
      endif
    case('LIQUID_FLUX')
      if (associated(towg%liquid_flux)) then
        sub_condition_ptr => towg%liquid_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%liquid_flux => sub_condition_ptr
      endif
    case('OIL_FLUX')
      if (associated(towg%oil_flux)) then
        sub_condition_ptr => towg%oil_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%oil_flux => sub_condition_ptr
      endif
    case('GAS_FLUX')
      if (associated(towg%gas_flux)) then
        sub_condition_ptr => towg%gas_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%gas_flux => sub_condition_ptr
      endif
    case('SOLVENT_FLUX')
      if (associated(towg%solvent_flux)) then
        sub_condition_ptr => towg%solvent_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%solvent_flux => sub_condition_ptr
      endif
    case('ENERGY_FLUX')
      if (associated(towg%energy_flux)) then
        sub_condition_ptr => towg%energy_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%energy_flux => sub_condition_ptr
      endif
    case('RATE')
      if (associated(towg%rate)) then
        sub_condition_ptr => towg%rate
      else
        sub_condition_ptr => FlowSubConditionCreate(option%nflowdof)
        towg%rate => sub_condition_ptr
      endif
    case('BHP_PRESSURE')
      if (associated(towg%bhp_pressure)) then
        sub_condition_ptr => towg%bhp_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%bhp_pressure => sub_condition_ptr
      endif
    case('OWC_Z','OWC_D','WGC_Z','WGC_D')
      if (associated(towg%owc_z)) then
        sub_condition_ptr => towg%owc_z
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%owc_z => sub_condition_ptr
      endif
    case('PCOW_OWC','PCWG_WGC')
      if (associated(towg%pcow_owc)) then
        sub_condition_ptr => towg%pcow_owc
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%pcow_owc => sub_condition_ptr
      endif      
    case('OGC_Z','OGC_D')
      if (associated(towg%ogc_z)) then
        sub_condition_ptr => towg%ogc_z
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%ogc_z => sub_condition_ptr
      endif
    case('PCOG_OGC')
      if (associated(towg%pcog_ogc)) then
        sub_condition_ptr => towg%pcog_ogc
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        towg%pcog_ogc => sub_condition_ptr
      endif
    case default
      call InputKeywordUnrecognized(input,sub_condition_name, &
                                    'towg condition,type',option)
  end select

  FlowTOWGSubConditionPtr => sub_condition_ptr

end function FlowTOWGSubConditionPtr


! ************************************************************************** !

function FlowSubConditionCreate(ndof)
  !
  ! Creates a sub_condition
  !
  ! Author: Glenn Hammond
  ! Date: 02/04/08
  !

  use Dataset_Ascii_class
  use Option_module

  implicit none

  type(flow_sub_condition_type), pointer :: FlowSubConditionCreate

  PetscInt :: ndof

  type(flow_sub_condition_type), pointer :: sub_condition
  class(dataset_ascii_type), pointer :: dataset_ascii

  allocate(sub_condition)
  sub_condition%units = ''
  sub_condition%itype = NULL_CONDITION
  sub_condition%isubtype = 0
  sub_condition%ctype = ''
  sub_condition%name = ''
  sub_condition%aux_real = UNINITIALIZED_DOUBLE
  nullify(sub_condition%gradient)
  nullify(sub_condition%dataset)

  ! by default, all dataset are of type dataset_ascii_type, unless overwritten
  dataset_ascii => DatasetAsciiCreate()
  call DatasetAsciiInit(dataset_ascii)
  dataset_ascii%array_width = ndof
  dataset_ascii%data_type = DATASET_REAL
  sub_condition%dataset => dataset_ascii
  nullify(dataset_ascii)

  FlowSubConditionCreate => sub_condition

end function FlowSubConditionCreate

! ************************************************************************** !

function GetFlowSubCondFromArrayByName(sub_condition_ptr_list,name)
  !
  ! returns a pointer to a subcondition with
  ! matching name
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/08
  !

  use Input_Aux_module
  use String_module

  implicit none

  type(flow_sub_condition_type), pointer :: GetFlowSubCondFromArrayByName
  type(sub_condition_ptr_type), pointer :: sub_condition_ptr_list(:)
  character(len=MAXWORDLENGTH) :: name

  PetscInt :: idof
  PetscInt :: length

  nullify(GetFlowSubCondFromArrayByName)
  length = len_trim(name)
  do idof = 1, size(sub_condition_ptr_list)
    if (length == len_trim(sub_condition_ptr_list(idof)%ptr%name) .and. &
        StringCompare(name,sub_condition_ptr_list(idof)%ptr%name,length)) then
      GetFlowSubCondFromArrayByName => sub_condition_ptr_list(idof)%ptr
      return
    endif
  enddo

  print *, 'GetFlowSubCondFromArrayByName() needs to be updated to include &
           &the general_condition_type.'
  stop

end function GetFlowSubCondFromArrayByName

! ************************************************************************** !

subroutine FlowSubConditionVerify(option, condition, sub_condition_name, &
                                  sub_condition, default_time_storage, &
                                  destroy_if_null)
  !
  ! Verifies the data in a subcondition
  !
  ! Author: Glenn Hammond
  ! Date: 02/04/08
  !
  use Time_Storage_module
  use Option_module
  use Dataset_module

  implicit none

  type(option_type) :: option
  type(flow_condition_type) :: condition
  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_sub_condition_type), pointer :: sub_condition
  type(time_storage_type), pointer :: default_time_storage
  PetscBool :: destroy_if_null

  character(len=MAXSTRINGLENGTH) :: header

  if (.not.associated(sub_condition)) return

  ! dataset is not optional
  if (.not.(associated(sub_condition%dataset%rarray) .or. &
            associated(sub_condition%dataset%rbuffer) .or. &
            ! if a dataset name is read, instead of data at this point
            len_trim(sub_condition%dataset%name) > 0 .or. &
            sub_condition%itype /= NULL_CONDITION)) then
    if (destroy_if_null) call FlowSubConditionDestroy(sub_condition)
    return
  endif

!  if (len_trim(sub_condition%ctype) == NULL_CONDITION) then
  if (sub_condition%itype == NULL_CONDITION) then
    option%io_buffer = 'TYPE of condition ' // trim(condition%name) // &
      ' ' // trim(sub_condition_name) // ' dataset not defined.'
    call PrintErrMsg(option)
  endif

  header = 'SUBSURFACE/FLOW_CONDITION/' // &
           trim(condition%name) // '/' // &
           trim(sub_condition_name) // '/Value(s)'
  call DatasetVerify(sub_condition%dataset,default_time_storage, &
                     header,option)
  header = 'SUBSURFACE/FLOW_CONDITION/' // &
           trim(condition%name) // '/' // &
           trim(sub_condition_name) // '/Gradient'
  call DatasetVerify(sub_condition%gradient,default_time_storage, &
                     header,option)

end subroutine FlowSubConditionVerify

! ************************************************************************** !

subroutine FlowConditionRead(condition,input,option)
  !
  ! Reads a condition from the input file
  !
  ! Author: Glenn Hammond
  ! Date: 10/31/07
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Time_Storage_module
  use Dataset_module

  implicit none

  type(flow_condition_type) :: condition
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: rate_unit_string
  character(len=MAXWORDLENGTH) :: energy_rate_unit_string
  character(len=MAXWORDLENGTH) :: internal_units
  type(flow_sub_condition_type), pointer :: pressure, flux, temperature, &
                                       concentration, enthalpy, rate, well,&
                                       sub_condition_ptr, saturation, &
                                       energy_rate, energy_flux
  PetscInt :: default_iphase
  PetscInt :: idof
  type(time_storage_type), pointer :: default_time_storage
  class(dataset_ascii_type), pointer :: dataset_ascii

  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read, &
                          ierr);CHKERRQ(ierr)

  default_iphase = 0

  default_time_storage => TimeStorageCreate()
  default_time_storage%is_cyclic = PETSC_FALSE
  default_time_storage%time_interpolation_method = INTERPOLATION_STEP

  rate_unit_string = 'not_assigned'
  energy_rate_unit_string = 'not_assigned'
  internal_units = 'not_assigned'

  pressure => FlowSubConditionCreate(option%nphase)
  pressure%name = 'pressure'
  flux => pressure
  rate => FlowSubConditionCreate(option%nflowspec)
  rate%name = 'rate'
  energy_rate => FlowSubConditionCreate(ONE_INTEGER)
  energy_rate%name = 'energy_rate'
  energy_flux => FlowSubConditionCreate(ONE_INTEGER)
  energy_flux%name = 'energy_flux'
  well => FlowSubConditionCreate(7 + option%nflowspec)
  well%name = 'well'
  saturation => FlowSubConditionCreate(option%nphase)
  saturation%name = 'saturation'
  temperature => FlowSubConditionCreate(ONE_INTEGER)
  temperature%name = 'temperature'
  concentration => FlowSubConditionCreate(ONE_INTEGER)
  concentration%name = 'concentration'
  enthalpy => FlowSubConditionCreate(option%nphase)
  enthalpy%name = 'enthalpy'

  condition%time_units = 'yr'
  condition%length_units = 'm'
  pressure%units = 'Pa'
  rate%units = 'kg/s'
  energy_rate%units = 'W'
  energy_flux%units = 'W/m^2'
  well%units = 'Pa'
  saturation%units = ' '
  temperature%units = 'C'
  concentration%units = 'M'
  enthalpy%units = 'kJ/mol'

  ! read the condition
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    internal_units = 'not_assigned'

    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CONDITION')

    select case(trim(word))

      case('UNITS') ! read default units for condition arguments
        do
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (InputError(input)) exit
          select case(trim(word))
            case('s','sec','min','hr','d','day','y','yr')
              condition%time_units = trim(word)
            case('mm','cm','m','met','meter','dm','km')
              condition%length_units = trim(word)
            case('Pa','KPa')
              pressure%units = trim(word)
            case('kg/s','kg/yr')
              rate%units = trim(word)
            case('W','J/yr')
              energy_rate%units = trim(word)
            case('W/m^2','J/m^2/yr')
              energy_flux%units = trim(word)
            case('m/s','m/yr')
              flux%units = trim(word)
            case('C','K')
              temperature%units = trim(word)
            case('M','mol/L')
              concentration%units = trim(word)
            case('kJ/mol')
              enthalpy%units = trim(word)
            case default
              call InputKeywordUnrecognized(input,word,'condition,units', &
                                            option)
          end select
        enddo
      case('CYCLIC')
        ! by default, is_cyclic is set to PETSC_FALSE
        default_time_storage%is_cyclic = PETSC_TRUE
      case('SYNC_TIMESTEP_WITH_UPDATE')
        condition%sync_time_with_update = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')
        call StringToUpper(word)
        select case(word)
          case('STEP')
            default_time_storage%time_interpolation_method = &
              INTERPOLATION_STEP
          case('LINEAR')
            default_time_storage%time_interpolation_method = &
              INTERPOLATION_LINEAR
          case default
            call InputKeywordUnrecognized(input,word, &
                                          'condition,interpolation', &
                                          option)
        end select
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          call StringToUpper(word)
          select case(trim(word))
            case('PRESSURE')
              sub_condition_ptr => pressure
            case('RATE')
              sub_condition_ptr => rate
            case('ENERGY_RATE')
              sub_condition_ptr => energy_rate
            case('WELL')
              sub_condition_ptr => well
            case('FLUX')
              sub_condition_ptr => flux
            case('ENERGY_FLUX')
              sub_condition_ptr => energy_flux
            case('SATURATION')
              sub_condition_ptr => saturation
            case('TEMPERATURE')
              sub_condition_ptr => temperature
            case('CONCENTRATION')
              sub_condition_ptr => concentration
            case('ENTHALPY')
              sub_condition_ptr => enthalpy
            case default
              call InputKeywordUnrecognized(input,word,'condition,type',option)
          end select
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'TYPE','CONDITION')
          call StringToUpper(word)
          sub_condition_ptr%ctype = word
          select case(word)
            case('DIRICHLET')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('NEUMANN')
              sub_condition_ptr%itype = NEUMANN_BC
            case('MASS_RATE')
              sub_condition_ptr%itype = MASS_RATE_SS
              rate_unit_string = 'kg/sec'
            case('TOTAL_MASS_RATE')
              sub_condition_ptr%itype = TOTAL_MASS_RATE_SS
              rate_unit_string = 'kg/sec'
            case('ENERGY_RATE')
              sub_condition_ptr%itype = ENERGY_RATE_SS
              energy_rate_unit_string = 'MJ/sec|MW'
            case('HETEROGENEOUS_ENERGY_RATE')
              sub_condition_ptr%itype = HET_ENERGY_RATE_SS
              energy_rate_unit_string = 'MJ/sec|MW'
            case('SCALED_MASS_RATE','SCALED_VOLUMETRIC_RATE', &
                 'SCALED_ENERGY_RATE')
              select case(word)
                case('SCALED_MASS_RATE')
                  sub_condition_ptr%itype = SCALED_MASS_RATE_SS
                  rate_unit_string = 'kg/sec'
                case('SCALED_VOLUMETRIC_RATE')
                  sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
                  rate_unit_string = 'm^3/sec'
                case('SCALED_ENERGY_RATE')
                  sub_condition_ptr%itype = SCALED_ENERGY_RATE_SS
                  energy_rate_unit_string = 'MW|MJ/sec'
              end select
              ! store name of type for error messaging below.
              string = word
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call InputPushCard(input,word,option)
                call StringToUpper(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('NEIGHBOR_PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('VOLUME')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" ' // trim(string)
                    call InputKeywordUnrecognized(input,word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in flow condition "' // &
                  trim(condition%name) // '" ' // trim(string)
                call PrintErrMsg(option)
                endif
            case('HYDROSTATIC')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('CONDUCTANCE')
              sub_condition_ptr%itype = HYDROSTATIC_CONDUCTANCE_BC
            case('ZERO_GRADIENT')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case('WELL','PRODUCTION_WELL', 'INJECTION_WELL')
              sub_condition_ptr%itype = WELL_SS
            case('SEEPAGE')
              sub_condition_ptr%itype = HYDROSTATIC_SEEPAGE_BC
            case('DIRICHLET_SEEPAGE')
              sub_condition_ptr%itype = DIRICHLET_SEEPAGE_BC
            case('DIRICHLET_CONDUCTANCE')
              sub_condition_ptr%itype = DIRICHLET_CONDUCTANCE_BC
            case('VOLUMETRIC_RATE')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
              rate_unit_string = 'm^3/sec'
            case('EQUILIBRIUM')
              sub_condition_ptr%itype = EQUILIBRIUM_SS
            case('UNIT_GRADIENT')
              if (.not.associated(sub_condition_ptr,pressure)) then
                option%io_buffer = 'unit_gradient flow condition type may &
                  &only be associated with a PRESSURE flow condition.'
                call PrintErrMsg(option)
              endif
              sub_condition_ptr%itype = UNIT_GRADIENT_BC
            case('HETEROGENEOUS_VOLUMETRIC_RATE')
              sub_condition_ptr%itype = HET_VOL_RATE_SS
              rate_unit_string = 'm^3/sec'
            case('HETEROGENEOUS_MASS_RATE')
              sub_condition_ptr%itype = HET_MASS_RATE_SS
              rate_unit_string = 'kg/sec'
            case('HETEROGENEOUS_DIRICHLET')
              sub_condition_ptr%itype = HET_DIRICHLET_BC
            case('HETEROGENEOUS_SEEPAGE')
              sub_condition_ptr%itype = HET_HYDROSTATIC_SEEPAGE_BC
            case('HETEROGENEOUS_CONDUCTANCE')
              sub_condition_ptr%itype = HET_HYDROSTATIC_CONDUCTANCE_BC
            case('HETEROGENEOUS_SURFACE_SEEPAGE')
              sub_condition_ptr%itype = HET_SURF_HYDROSTATIC_SEEPAGE_BC
            case('SPILLOVER')
              sub_condition_ptr%itype = SPILLOVER_BC
            case('SURFACE_DIRICHLET')
              sub_condition_ptr%itype = SURFACE_DIRICHLET
            case('SURFACE_ZERO_GRADHEIGHT')
              sub_condition_ptr%itype = SURFACE_ZERO_GRADHEIGHT
            case('SURFACE_SPILLOVER')
              sub_condition_ptr%itype = SURFACE_SPILLOVER
            case default
              call InputKeywordUnrecognized(input,word,'condition bc type',option)
          end select
        enddo
        call InputPopBlock(input,option)
      case('IPHASE')
        call InputReadInt(input,option,default_iphase)
        call InputErrorMsg(input,option,'IPHASE','CONDITION')
      case('DATUM')
        dataset_ascii => DatasetAsciiCreate()
        call DatasetAsciiInit(dataset_ascii)
        dataset_ascii%array_width = 3
        dataset_ascii%data_type = DATASET_REAL
        condition%datum => dataset_ascii
        nullify(dataset_ascii)
        internal_units = 'meter'
        call ConditionReadValues(input,option,word, &
                                 condition%datum,word,internal_units)
      case('GRADIENT','GRAD')
        call InputPushBlock(input,option)
        do
          internal_units = 'not_assigned'
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          select case(trim(word))
            case('PRES','PRESS','PRESSURE')
              sub_condition_ptr => pressure
              internal_units = 'Pa/meter'
            case('RATE')
              sub_condition_ptr => rate
              internal_units = 'kg/sec-meter'
            case('ENERGY_RATE')
              sub_condition_ptr => energy_rate
              internal_units = 'MW/meter|MJ/sec-meter'
            case('WELL')
              sub_condition_ptr => well
              internal_units = 'Pa/meter'
            case('FLUX')
              sub_condition_ptr => flux
              internal_units = 'm/sec-m|unitless/sec'
            case('SATURATION')
              sub_condition_ptr => saturation
              internal_units = 'unitless/meter'
            case('TEMP','TEMPERATURE')
              sub_condition_ptr => temperature
              internal_units = 'temperature/m'
            case('CONC','CONCENTRATION')
              sub_condition_ptr => concentration
              internal_units = 'unitless'
            case('H','ENTHALPY')
              sub_condition_ptr => enthalpy
              internal_units = 'kJ/mol-meter'
            case default
              call InputKeywordUnrecognized(input,word, &
                     'FLOW CONDITION,GRADIENT,TYPE',option)
          end select
          dataset_ascii => DatasetAsciiCreate()
          call DatasetAsciiInit(dataset_ascii)
          dataset_ascii%array_width = 3
          dataset_ascii%data_type = DATASET_REAL
          sub_condition_ptr%gradient => dataset_ascii
          nullify(dataset_ascii)
          call ConditionReadValues(input,option,word, &
                                   sub_condition_ptr%gradient, &
                                   word,internal_units)
          nullify(sub_condition_ptr)
        enddo
        call InputPopBlock(input,option)
      case('TEMPERATURE','TEMP')
        internal_units = 'C'
        call ConditionReadValues(input,option,word, &
                                 temperature%dataset, &
                                 temperature%units,internal_units)
      case('ENTHALPY','H')
        internal_units = 'kJ/mol'
        call ConditionReadValues(input,option,word, &
                                 enthalpy%dataset, &
                                 enthalpy%units,internal_units)
      case('PRESSURE','PRES','PRESS')
        internal_units = 'Pa'
        call ConditionReadValues(input,option,word, &
                                 pressure%dataset, &
                                 pressure%units,internal_units)
      case('RATE')
        internal_units = rate_unit_string
        call ConditionReadValues(input,option,word, &
                                 rate%dataset, &
                                 rate%units,internal_units)
      case('ENERGY_FLUX')
        input%force_units = PETSC_TRUE
        internal_units = 'MW/m^2|MJ/m^2-sec'
        call ConditionReadValues(input,option,word, &
                                 energy_flux%dataset, &
                                 energy_flux%units,internal_units)
        input%force_units = PETSC_FALSE
      case('ENERGY_RATE')
        input%force_units = PETSC_TRUE
        internal_units = energy_rate_unit_string
        input%err_buf = word
        call ConditionReadValues(input,option,word, &
                                 energy_rate%dataset, &
                                 energy_rate%units,internal_units)
        input%force_units = PETSC_FALSE
      case('WELL')
        internal_units = 'Pa'
        call ConditionReadValues(input,option,word, &
                                 well%dataset, &
                                 well%units,internal_units)
      case('FLUX','VELOCITY','VEL')
        internal_units = 'meter/sec'
        call ConditionReadValues(input,option,word, &
                                 pressure%dataset, &
                                 pressure%units,internal_units)
      case('CONC','CONCENTRATION')
        internal_units = 'unitless'
        call ConditionReadValues(input,option,word, &
                                 concentration%dataset, &
                                 concentration%units,internal_units)
      case('SAT','SATURATION')
        internal_units = 'unitless'
        call ConditionReadValues(input,option,word, &
                                 saturation%dataset, &
                                 saturation%units,internal_units)
      case('CONDUCTANCE')
        call InputReadDouble(input,option,pressure%aux_real(1))
        call InputErrorMsg(input,option,'CONDUCTANCE','CONDITION')
      case default
        call InputKeywordUnrecognized(input,word,'flow condition',option)
    end select

  enddo
  call InputPopBlock(input,option)

  ! check whether
  if (default_iphase == 0) then
    condition%iphase = 1
  else
    condition%iphase = default_iphase
  endif

  !geh: simple check to ensure that DIRICHLET_SEEPAGE and 
  !     DIRICHLET_CONDUCTANCE_BC are only used in TH and RICHARDS
  select case(option%iflowmode)
    case(RICHARDS_MODE,TH_MODE)
    case default
      if (pressure%itype == DIRICHLET_SEEPAGE_BC .or. &
          pressure%itype == DIRICHLET_CONDUCTANCE_BC) then
        option%io_buffer = 'DIRICHLET_SEEPAGE_BC and DIRICHLET_CONDUCTANCE_BC &
          &only supported for RICHARDS and TH.'
        call PrintErrMsg(option)
      endif
  end select

  ! datum is not required
  string = trim(condition%name) // '/' // 'Datum'
  call DatasetVerify(condition%datum,default_time_storage,string,option)

  ! check to ensure that a rate condition is not of type pressure
  if (associated(rate)) then
    select case(rate%itype)
      case(DIRICHLET_BC,NEUMANN_BC,HYDROSTATIC_BC,UNIT_GRADIENT_BC, &
           HYDROSTATIC_CONDUCTANCE_BC,ZERO_GRADIENT_BC,HYDROSTATIC_SEEPAGE_BC, &
           DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC,SURFACE_DIRICHLET, &
           SURFACE_SPILLOVER,HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC,&
           HET_HYDROSTATIC_CONDUCTANCE_BC)
        option%io_buffer = 'RATE condition must not be of type: dirichlet, &
          &neumann, zero_gradient, dirichlet_zero_gradient, hydrostatic, &
          &seepage, or conductance".'
        call PrintErrMsg(option)
    end select
  endif
  ! check to ensure that a pressure condition is not of type rate
  if (associated(pressure)) then
    select case(pressure%itype)
      case(MASS_RATE_SS,SCALED_MASS_RATE_SS,VOLUMETRIC_RATE_SS, &
           SCALED_VOLUMETRIC_RATE_SS,EQUILIBRIUM_SS)
        option%io_buffer = 'PRESSURE or FLUX condition must not be of type: &
          &mass_rate, scaled_mass_rate, volumetric_rate, &
          &scaled_volumetric_rate, equilibrium, or production_well.'
        call PrintErrMsg(option)
    end select
  endif

  ! verify the datasets
  word = 'pressure/flux'
  call FlowSubConditionVerify(option,condition,word,pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,rate, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy_flux'
  call FlowSubConditionVerify(option,condition,word,energy_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy_rate'
  call FlowSubConditionVerify(option,condition,word,energy_rate, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'well'
  call FlowSubConditionVerify(option,condition,word,well, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,temperature, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'saturation'
  call FlowSubConditionVerify(option,condition,word,saturation, &
                              default_time_storage, &
                              PETSC_TRUE)

  word = 'concentration'
  call FlowSubConditionVerify(option,condition,word,concentration, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'enthalpy'
  call FlowSubConditionVerify(option,condition,word,enthalpy, &
                              default_time_storage, &
                              PETSC_TRUE)

  select case(option%iflowmode)
    case default
      option%io_buffer = 'The flow mode not supported in original &
        &FlowConditionRead.'
      call PrintMsg(option)
    case(G_MODE)
      option%io_buffer = 'General mode not supported in original &
        &FlowConditionRead.'
      call PrintMsg(option)
    case(H_MODE)
      option%io_buffer = 'Hydrate mode not supported in original &
        &FlowConditionRead.'
      call PrintMsg(option)
    case(WF_MODE)
      option%io_buffer = 'WIPP Flow mode not supported in original &
        &FlowConditionRead.'
      call PrintMsg(option)
    case(TOIL_IMS_MODE)
      option%io_buffer = 'TOilIms mode not supported in original &
        &FlowConditionRead.'
      call PrintMsg(option)
    case(MPH_MODE,IMS_MODE,FLASH2_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)&
           .and. .not.associated(well) .and. .not.associated(saturation)) then
        option%io_buffer = 'pressure, rate and saturation condition null in &
                           &condition: ' // trim(condition%name)
        call PrintErrMsg(option)
      endif

      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif
      if (associated(saturation)) then
        condition%saturation => saturation
      endif


      if (.not.associated(temperature) .and. .not.associated(energy_rate)) then
        option%io_buffer = 'temperature and energy rate condition null &
          &in condition: ' // trim(condition%name)
        call PrintErrMsg(option)
      endif
      if (associated(temperature)) then
        condition%temperature => temperature
      endif
      if (associated(energy_flux)) then
        condition%energy_flux => energy_flux
      endif
      if (associated(energy_rate)) then
        condition%energy_rate => energy_rate
      endif

      if (.not.associated(concentration)) then
        option%io_buffer = 'concentration condition null in condition: ' // &
                            trim(condition%name)
        call PrintErrMsg(option)
      endif
      condition%concentration => concentration

      if (.not.associated(enthalpy)) then
        option%io_buffer = 'enthalpy condition null in condition: ' // &
                            trim(condition%name)
        call PrintErrMsg(option)
      endif
      condition%enthalpy => enthalpy

      condition%num_sub_conditions = 4
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, 4
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo

      ! must be in this order, which matches the dofs i problem
      if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      if (associated(well)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
      if (associated(saturation)) condition%sub_condition_ptr(ONE_INTEGER)%ptr &
                                  => saturation
      condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      condition%sub_condition_ptr(THREE_INTEGER)%ptr => concentration
      if (associated(enthalpy)) condition%sub_condition_ptr(FOUR_INTEGER)%ptr => enthalpy
      if (associated(energy_rate)) &
        condition%sub_condition_ptr(FOUR_INTEGER)%ptr => energy_rate

      allocate(condition%itype(FIVE_INTEGER))
      condition%itype = 0
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
      if (associated(well)) condition%itype(ONE_INTEGER) = well%itype
      if (associated(saturation)) condition%itype(ONE_INTEGER) = &
                                    saturation%itype
      condition%itype(TWO_INTEGER) = temperature%itype
      condition%itype(THREE_INTEGER) = concentration%itype
      if (associated(enthalpy)) condition%itype(FOUR_INTEGER) = concentration%itype
      if (associated(energy_rate)) condition%itype(FOUR_INTEGER) = energy_rate%itype

    case(TH_MODE,TH_TS_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)&
           .and. .not.associated(well) .and. .not.associated(saturation)) then
        option%io_buffer = 'pressure, rate and saturation condition null in &
          &condition: ' // trim(condition%name)
        call PrintErrMsg(option)
      endif

      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif
      if (associated(saturation)) then
        condition%saturation => saturation
      endif

      if (.not.associated(temperature) .and. .not.associated(energy_rate) &
          .and. .not.associated(energy_flux)) then
        option%io_buffer = 'temperature, energy_flux, and energy_rate &
          &condition null in condition: ' // trim(condition%name)
        call PrintErrMsg(option)
      endif
      if (associated(temperature) .and. associated(energy_rate) ) then
        option%io_buffer = 'Both, temperature and energy_rate cannot be &
                           &specified in condition: ' // trim(condition%name)
        call PrintErrMsg(option)
      endif
      if (associated(temperature)) condition%temperature => temperature
      if (associated(energy_flux)) condition%energy_flux => energy_flux
      if (associated(energy_rate)) condition%energy_rate => energy_rate

      if (associated(enthalpy)) then
        option%io_buffer = 'enthalpy condition not supported in TH mode: ' // &
                            trim(condition%name)
        call PrintErrMsg(option)
      endif
      if (associated(enthalpy)) condition%enthalpy => enthalpy

      condition%num_sub_conditions = TWO_INTEGER
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, 2
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo

      ! must be in this order, which matches the dofs i problem
      if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      if (associated(well)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
      if (associated(saturation)) condition%sub_condition_ptr(ONE_INTEGER)%ptr &
                                  => saturation
      if ( associated(temperature)) &
        condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      if (associated(energy_flux)) condition%sub_condition_ptr(TWO_INTEGER)%ptr => energy_flux
      if (associated(energy_rate)) condition%sub_condition_ptr(TWO_INTEGER)%ptr => energy_rate

      allocate(condition%itype(TWO_INTEGER))
      condition%itype = 0
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
      if (associated(well)) condition%itype(ONE_INTEGER) = well%itype
      if (associated(saturation)) condition%itype(ONE_INTEGER) = &
                                    saturation%itype
      if (associated(temperature)) condition%itype(TWO_INTEGER) = temperature%itype
      if (associated(energy_flux)) condition%itype(TWO_INTEGER) = energy_flux%itype
      if (associated(energy_rate)) condition%itype(TWO_INTEGER) = energy_rate%itype

!#if 0
    case(MIS_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)&
           .and. .not.associated(well)) then
        option%io_buffer = 'pressure and rate condition null in &
                           &condition: ' // trim(condition%name)
        call PrintErrMsg(option)
      endif

      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif

      if (.not.associated(concentration)) then
        option%io_buffer = 'concentration condition null in condition: ' // &
                            trim(condition%name)
        call PrintErrMsg(option)
      endif
      condition%concentration => concentration

#if 0
      if (.not.associated(temperature)) then
        option%io_buffer = 'temperature condition null in condition: ' // &
                            trim(condition%name)
        call PrintErrMsg(option)
      endif
      condition%temperature => temperature

      if (.not.associated(enthalpy)) then
        option%io_buffer = 'enthalpy condition null in condition: ' // &
                            trim(condition%name)
        call PrintErrMsg(option)
      endif
      condition%enthalpy => enthalpy
#endif

      condition%num_sub_conditions = 2
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, 2
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo

      ! must be in this order, which matches the dofs in problem
      if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      if (associated(well)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
!     condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      condition%sub_condition_ptr(TWO_INTEGER)%ptr => concentration
!     if (associated(enthalpy)) condition%sub_condition_ptr(FOUR_INTEGER)%ptr => enthalpy

      allocate(condition%itype(TWO_INTEGER))
      condition%itype = 0
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
      if (associated(well)) condition%itype(ONE_INTEGER) = well%itype
      condition%itype(TWO_INTEGER) = concentration%itype
!#endif

    case(RICHARDS_MODE,RICHARDS_TS_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate) .and. &
          .not.associated(saturation) .and. .not.associated(well)) then
        option%io_buffer = 'pressure, rate and saturation condition null in &
                           &condition: ' // trim(condition%name)
        call PrintErrMsg(option)
      endif

      if (associated(saturation)) then
        condition%saturation => saturation
      endif
      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif

      condition%num_sub_conditions = 1
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      if (associated(pressure)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      elseif (associated(saturation)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => saturation
      elseif (associated(rate)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      elseif (associated(well)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
      endif

      allocate(condition%itype(ONE_INTEGER))
      if (associated(pressure)) then
        condition%itype(ONE_INTEGER) = pressure%itype
      else if (associated(saturation)) then
        condition%itype(ONE_INTEGER) = saturation%itype
      else if (associated(rate)) then
        condition%itype(ONE_INTEGER) = rate%itype
      else if (associated(well)) then
        condition%itype(ONE_INTEGER) = well%itype
      endif

      ! these are not used with richards
      if (associated(temperature)) call FlowSubConditionDestroy(temperature)
      if (associated(enthalpy)) call FlowSubConditionDestroy(enthalpy)

  end select

  condition%default_time_storage => default_time_storage

  call PetscLogEventEnd(logging%event_flow_condition_read,ierr);CHKERRQ(ierr)

end subroutine FlowConditionRead

! ************************************************************************** !

subroutine FlowConditionGeneralRead(condition,input,option)
  !
  ! Reads a condition from the input file for
  ! general mode
  !
  ! Author: Glenn Hammond
  ! Date: 09/14/11
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Time_Storage_module
  use Dataset_module

  ! needed for STATES
  use General_Aux_module

  implicit none

  type(flow_condition_type) :: condition
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: rate_string, internal_units
  character(len=MAXWORDLENGTH) :: word
  type(flow_general_condition_type), pointer :: general
  type(flow_sub_condition_type), pointer :: sub_condition_ptr
  PetscReal :: default_time
  PetscInt :: default_iphase
  class(dataset_base_type), pointer :: default_flow_dataset
  class(dataset_base_type), pointer :: default_gradient
  PetscInt :: idof, i
  PetscBool :: default_is_cyclic
  type(time_storage_type), pointer :: default_time_storage
  class(dataset_ascii_type), pointer :: dataset_ascii
  character(len=MAXWORDLENGTH) :: flow_mode_chars
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read, &
                          ierr);CHKERRQ(ierr)

  select case(option%iflowmode)
    case(G_MODE)
      flow_mode_chars = 'General Mode'
    case(WF_MODE)
      flow_mode_chars = 'WIPP Flow Mode'
  end select

  rate_string = 'not_assigned'
  internal_units = 'not_assigned'

  default_time = 0.d0
  default_iphase = 0

  default_time_storage => TimeStorageCreate()
  default_time_storage%is_cyclic = PETSC_FALSE
  default_time_storage%time_interpolation_method = INTERPOLATION_STEP

  select case(option%iflowmode)
    case(G_MODE,WF_MODE)
      general => FlowGeneralConditionCreate(option)
      condition%general => general
  end select

  ! read the condition
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    internal_units = 'not_assigned'

    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CONDITION')

    select case(trim(word))

      case('CYCLIC')
        ! by default, is_cyclic is set to PETSC_FALSE
        default_time_storage%is_cyclic = PETSC_TRUE
      case('SYNC_TIMESTEP_WITH_UPDATE')
        condition%sync_time_with_update = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')
        call StringToUpper(word)
        select case(word)
          case('STEP')
            default_time_storage%time_interpolation_method = &
              INTERPOLATION_STEP
          case('LINEAR')
            default_time_storage%time_interpolation_method = &
              INTERPOLATION_LINEAR
        end select
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          call StringToUpper(word)
          select case(option%iflowmode)
            case(G_MODE,WF_MODE)
              sub_condition_ptr => &
                FlowGeneralSubConditionPtr(input,word,general,option)
          end select
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'TYPE','CONDITION')
          call StringToUpper(word)
          sub_condition_ptr%ctype = word
          select case(word)
            case('DIRICHLET')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('NEUMANN')
              sub_condition_ptr%itype = NEUMANN_BC
            case('HYDROSTATIC')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('CONDUCTANCE')
              sub_condition_ptr%itype = HYDROSTATIC_CONDUCTANCE_BC
            case('SEEPAGE')
              sub_condition_ptr%itype = HYDROSTATIC_SEEPAGE_BC
            case('MASS_RATE')
              sub_condition_ptr%itype = MASS_RATE_SS
              rate_string = 'kg/sec'
            case('TOTAL_MASS_RATE')
              sub_condition_ptr%itype = TOTAL_MASS_RATE_SS
              rate_string = 'kg/sec'
            case('SCALED_MASS_RATE')
              sub_condition_ptr%itype = SCALED_MASS_RATE_SS
              rate_string = 'kg/sec'
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call InputPushCard(input,word,option)
                call StringToUpper(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('NEIGHBOR_PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('VOLUME')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" scaled_mass_rate type'
                    call InputKeywordUnrecognized(input,word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in &
                  &flow condition "' // trim(condition%name) // &
                  '" scaled_mass_rate type'
                call PrintErrMsg(option)
              endif
            case('VOLUMETRIC_RATE')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'
            case('SCALED_VOLUMETRIC_RATE')
              sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call InputPushCard(input,word,option)
                call StringToUpper(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('NEIGHBOR_PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('VOLUME')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" scaled_volumetric_rate type'
                    call InputKeywordUnrecognized(input,word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in &
                  &flow condition "' // trim(condition%name) // &
                  '" scaled_volumetric_rate type'
                call PrintErrMsg(option)
              endif
            case('HETEROGENEOUS_VOLUMETRIC_RATE')
              sub_condition_ptr%itype = HET_VOL_RATE_SS
              rate_string = 'm^3/sec'
            case('HETEROGENEOUS_MASS_RATE')
              sub_condition_ptr%itype = HET_MASS_RATE_SS
              rate_string = 'kg/sec'
            case('HETEROGENEOUS_DIRICHLET')
              sub_condition_ptr%itype = HET_DIRICHLET_BC
            case('HETEROGENEOUS_SURFACE_SEEPAGE')
              sub_condition_ptr%itype = HET_SURF_HYDROSTATIC_SEEPAGE_BC
            case default
              call InputKeywordUnrecognized(input,word, &
                                            'flow condition,type',option)
          end select
        enddo
        call InputPopBlock(input,option)
      case('DATUM')
        dataset_ascii => DatasetAsciiCreate()
        call DatasetAsciiInit(dataset_ascii)
        dataset_ascii%array_width = 3
        dataset_ascii%data_type = DATASET_REAL
        condition%datum => dataset_ascii
        nullify(dataset_ascii)
        internal_units = 'meter'
        call ConditionReadValues(input,option,word,condition%datum, &
                                 word,internal_units)
      case('GRADIENT')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          call StringToUpper(word)
          select case(option%iflowmode)
            case(G_MODE,WF_MODE)
              sub_condition_ptr => &
                FlowGeneralSubConditionPtr(input,word,general,option)
          end select
          dataset_ascii => DatasetAsciiCreate()
          call DatasetAsciiInit(dataset_ascii)
          dataset_ascii%array_width = 3
          dataset_ascii%data_type = DATASET_REAL
          sub_condition_ptr%gradient => dataset_ascii
          nullify(dataset_ascii)
          internal_units = 'unitless/meter'
          call ConditionReadValues(input,option,word, &
                                   sub_condition_ptr%gradient, &
                                   word,internal_units)
          nullify(sub_condition_ptr)
        enddo
        call InputPopBlock(input,option)
      case('CONDUCTANCE')
        word = 'LIQUID_PRESSURE'
        select case(option%iflowmode)
          case(G_MODE,WF_MODE)
            sub_condition_ptr => &
                FlowGeneralSubConditionPtr(input,word,general,option)
        end select
        call InputReadDouble(input,option,sub_condition_ptr%aux_real(1))
        call InputErrorMsg(input,option,'LIQUID_CONDUCTANCE','CONDITION')
      case('LIQUID_PRESSURE','GAS_PRESSURE','LIQUID_SATURATION', &
           'GAS_SATURATION', 'TEMPERATURE','MOLE_FRACTION','RATE', &
           'LIQUID_FLUX','GAS_FLUX','ENERGY_FLUX','RELATIVE_HUMIDITY')
        select case(option%iflowmode)
          case(G_MODE,WF_MODE)
            sub_condition_ptr => &
                FlowGeneralSubConditionPtr(input,word,general,option)
        end select
        internal_units = 'not_assigned'
        select case(trim(word))
          case('LIQUID_PRESSURE','GAS_PRESSURE')
            internal_units = 'Pa'
          case('LIQUID_SATURATION','GAS_SATURATION','MOLE_FRACTION', &
                'RELATIVE_HUMIDITY')
            internal_units = 'unitless'
          case('TEMPERATURE')
            internal_units = 'C'
          case('RATE')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            select case(option%iflowmode)
              case(WF_MODE)
                internal_units = trim(rate_string) // ',' // trim(rate_string)
              case(G_MODE)
                internal_units = trim(rate_string) // ',' // &
                  trim(rate_string) // ',MJ/sec|MW'
            end select
          case('LIQUID_FLUX','GAS_FLUX')
            internal_units = 'meter/sec'
          case('ENERGY_FLUX')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            internal_units = 'MW/m^2|MJ/m^2-sec'
        end select
        call ConditionReadValues(input,option,word, &
                                 sub_condition_ptr%dataset, &
                                 sub_condition_ptr%units,internal_units)
        input%force_units = PETSC_FALSE
        select case(word)
          case('LIQUID_SATURATION') ! convert to gas saturation
            if (associated(sub_condition_ptr%dataset%rbuffer)) then
              sub_condition_ptr%dataset%rbuffer(:) = 1.d0 - &
                sub_condition_ptr%dataset%rbuffer(:)
            endif
            sub_condition_ptr%dataset%rarray(:) = 1.d0 - &
              sub_condition_ptr%dataset%rarray(:)
        end select
      case default
        call InputKeywordUnrecognized(input,word,'flow condition',option)
    end select

  enddo
  call InputPopBlock(input,option)

  ! datum is not required
  string = 'SUBSURFACE/FLOW_CONDITION' // trim(condition%name) // '/Datum'
  call DatasetVerify(condition%datum,default_time_storage,string,option)

  ! need mole fraction and some sort of saturation
  if (.false.) then
    ! neumann or mass/volumetric flux
    ! need temperature
    if (.not.associated(general%mole_fraction) .and. &
        .not.associated(general%gas_saturation)) then
      option%io_buffer = trim(flow_mode_chars) // ' flux condition must &
        &include a MOLE_FRACTION or GAS/LIQUID_SATURATION.'
      call PrintErrMsg(option)
    endif
    if (associated(general%mole_fraction) .and. &
        associated(general%gas_saturation)) then
      option%io_buffer = trim(flow_mode_chars) // ' flux condition must &
        &include only a MOLE_FRACTION or GAS/LIQUID_SATURATION, not both.'
      call PrintErrMsg(option)
    endif
    if (.not.associated(general%temperature)) then
      option%io_buffer = trim(flow_mode_chars) // ' flux condition must &
        &include a temperature.'
      call PrintErrMsg(option)
    endif
  else
    if (associated(general%rate)) then
      condition%iphase = ANY_STATE
    elseif (associated(general%liquid_flux) .and. &
            associated(general%gas_flux) .and. &
            (option%iflowmode == WF_MODE .or. &
             associated(general%energy_flux) .or. &
             associated(general%temperature))) then
      condition%iphase = ANY_STATE
    else
      ! some sort of dirichlet-based pressure, temperature, etc.
      if (option%iflowmode == G_MODE) then
        if (.not.associated(general%liquid_pressure) .and. &
            .not.associated(general%gas_pressure)) then
          option%io_buffer = 'General Mode non-rate condition must include &
            &a liquid or gas pressure'
          call PrintErrMsg(option)
        endif
        if (.not.associated(general%mole_fraction) .and. &
            .not.associated(general%relative_humidity) .and. &
            .not.associated(general%gas_saturation)) then
          option%io_buffer = 'General Mode non-rate condition must include &
            &a mole fraction, relative humidity, or gas/liquid saturation'
          call PrintErrMsg(option)
        endif
        if (.not.associated(general%temperature)) then
          option%io_buffer = 'General Mode non-rate condition must include &
            &a temperature'
          call PrintErrMsg(option)
        endif
        if ( associated(general%gas_pressure) .and. &
             associated(general%gas_saturation) .and. &
             associated(general%liquid_pressure) .and. &
             (associated(general%mole_fraction) .or. &
              associated(general%relative_humidity)) ) then
          ! multiphase condition
          condition%iphase = MULTI_STATE
        else if (associated(general%gas_pressure) .and. &
                associated(general%gas_saturation)) then
          ! two phase condition
          condition%iphase = TWO_PHASE_STATE
        else if (associated(general%liquid_pressure) .and. &
                 associated(general%mole_fraction)) then
          ! liquid phase condition
          condition%iphase = LIQUID_STATE
        else if (associated(general%gas_pressure) .and. &
                 (associated(general%mole_fraction) .or. &
                  associated(general%relative_humidity))) then
          ! gas phase condition
          condition%iphase = GAS_STATE
        endif
      else
        if (.not.associated(general%liquid_pressure)) then
          option%io_buffer = 'WIPP Flow Mode non-rate condition must include &
            &a liquid pressure'
          call PrintErrMsg(option)
        endif
        if (.not.associated(general%gas_saturation)) then
          option%io_buffer = 'WIPP Flow Mode non-rate condition must include &
            &a gas saturation'
          call PrintErrMsg(option)
        endif
        condition%iphase = TWO_PHASE_STATE
      endif
    endif
    if (condition%iphase == NULL_STATE) then
      option%io_buffer = 'General Phase non-rate/flux condition contains &
        &an unsupported combination of primary dependent variables.'
      call PrintErrMsg(option)
    endif
  endif

  ! verify the datasets
  word = 'liquid pressure'
  call FlowSubConditionVerify(option,condition,word,general%liquid_pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas pressure'
  call FlowSubConditionVerify(option,condition,word,general%gas_pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas saturation'
  call FlowSubConditionVerify(option,condition,word,general%gas_saturation, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'relative humidity'
  call FlowSubConditionVerify(option,condition,word,general%relative_humidity, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'mole fraction'
  call FlowSubConditionVerify(option,condition,word,general%mole_fraction, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,general%temperature, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'liquid flux'
  call FlowSubConditionVerify(option,condition,word,general%liquid_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas flux'
  call FlowSubConditionVerify(option,condition,word,general%gas_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy flux'
  call FlowSubConditionVerify(option,condition,word,general%energy_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,general%rate, &
                              default_time_storage, &
                              PETSC_TRUE)

  condition%num_sub_conditions = 0
  i = 0
  if (associated(general%liquid_pressure)) &
    i = i + 1
  if (associated(general%gas_pressure)) &
    i = i + 1
  if (associated(general%gas_saturation)) &
    i = i + 1
  if (associated(general%relative_humidity)) &
    i = i + 1
  if (associated(general%mole_fraction)) &
    i = i + 1
  if (associated(general%temperature)) &
    i = i + 1
  if (associated(general%liquid_flux)) &
    i = i + 1
  if (associated(general%gas_flux)) &
    i = i + 1
  if (associated(general%energy_flux)) &
    i = i + 1
  if (associated(general%rate)) &
    i = i + 1
  condition%num_sub_conditions = i
  allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    nullify(condition%sub_condition_ptr(idof)%ptr)
  enddo
  i = 0
  if (associated(general%liquid_pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%liquid_pressure
  endif
  if (associated(general%gas_pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%gas_pressure
  endif
  if (associated(general%gas_saturation)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%gas_saturation
  endif
  if (associated(general%relative_humidity)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%relative_humidity
  endif
  if (associated(general%mole_fraction)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%mole_fraction
  endif
  if (associated(general%temperature)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%temperature
  endif
  if (associated(general%liquid_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%liquid_flux
  endif
  if (associated(general%gas_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%gas_flux
  endif
  if (associated(general%energy_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%energy_flux
  endif
  if (associated(general%rate)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%rate
  endif

  ! set condition types
  allocate(condition%itype(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    condition%itype(idof) = condition%sub_condition_ptr(idof)%ptr%itype
  enddo

  condition%default_time_storage => default_time_storage

  call PetscLogEventEnd(logging%event_flow_condition_read,ierr);CHKERRQ(ierr)

end subroutine FlowConditionGeneralRead

! ************************************************************************** !

subroutine FlowConditionHydrateRead(condition,input,option)

  !
  ! Reads a condition from the input file for
  ! hydrate mode
  !
  ! Author: Michael Nole
  ! Date: 07/22/19
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Time_Storage_module
  use Dataset_module

  ! needed for STATES
  use Hydrate_Aux_module

  implicit none

  type(flow_condition_type) :: condition
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: rate_string, internal_units
  character(len=MAXWORDLENGTH) :: word
  type(flow_hydrate_condition_type), pointer :: hydrate
  type(flow_sub_condition_type), pointer :: sub_condition_ptr
  PetscReal :: default_time
  PetscInt :: default_iphase
  class(dataset_base_type), pointer :: default_flow_dataset
  class(dataset_base_type), pointer :: default_gradient
  PetscInt :: idof, i
  PetscBool :: default_is_cyclic
  type(time_storage_type), pointer :: default_time_storage
  class(dataset_ascii_type), pointer :: dataset_ascii
  character(len=MAXWORDLENGTH) :: flow_mode_chars
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read, &
                          ierr);CHKERRQ(ierr)

  select case(option%iflowmode)
    case(H_MODE)
      flow_mode_chars = 'Hydrate Mode'
  end select

  rate_string = 'not_assigned'
  internal_units = 'not_assigned'

  default_time = 0.d0
  default_iphase = 0

  default_time_storage => TimeStorageCreate()
  default_time_storage%is_cyclic = PETSC_FALSE
  default_time_storage%time_interpolation_method = INTERPOLATION_STEP

  select case(option%iflowmode)
    case(H_MODE)
      hydrate => FlowHydrateConditionCreate(option)
      condition%hydrate => hydrate
  end select

  ! read the condition
  input%ierr = 0
  call InputPushBlock(input,option)
  do

    internal_units = 'not_assigned'

    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CONDITION')

    select case(trim(word))

      case('CYCLIC')
        ! by default, is_cyclic is set to PETSC_FALSE
        default_time_storage%is_cyclic = PETSC_TRUE
      case('SYNC_TIMESTEP_WITH_UPDATE')
        condition%sync_time_with_update = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')
        call StringToUpper(word)
        select case(word)
          case('STEP')
            default_time_storage%time_interpolation_method = &
              INTERPOLATION_STEP
          case('LINEAR')
            default_time_storage%time_interpolation_method = &
              INTERPOLATION_LINEAR
        end select
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          call StringToUpper(word)
          select case(option%iflowmode)
            case(H_MODE)
              sub_condition_ptr => &
                FlowHydrateSubConditionPtr(input,word,hydrate,option)
          end select
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'TYPE','CONDITION')
          call StringToUpper(word)
          sub_condition_ptr%ctype = word
          select case(word)
            case('DIRICHLET')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('NEUMANN')
              sub_condition_ptr%itype = NEUMANN_BC
            case('HYDROSTATIC')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('CONDUCTANCE')
              sub_condition_ptr%itype = HYDROSTATIC_CONDUCTANCE_BC
            case('SEEPAGE')
              sub_condition_ptr%itype = HYDROSTATIC_SEEPAGE_BC
            case('MASS_RATE')
              sub_condition_ptr%itype = MASS_RATE_SS
              rate_string = 'kg/sec'
            case('TOTAL_MASS_RATE')
              sub_condition_ptr%itype = TOTAL_MASS_RATE_SS
              rate_string = 'kg/sec'
            case('SCALED_MASS_RATE')
              sub_condition_ptr%itype = SCALED_MASS_RATE_SS
              rate_string = 'kg/sec'
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call InputPushCard(input,word,option)
                call StringToUpper(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('NEIGHBOR_PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('VOLUME')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" scaled_mass_rate type'
                    call InputKeywordUnrecognized(input,word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in &
                  &flow condition "' // trim(condition%name) // &
                  '" scaled_mass_rate type'
                call printErrMsg(option)
              endif
            case('VOLUMETRIC_RATE')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'
            case('SCALED_VOLUMETRIC_RATE')
              sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call InputPushCard(input,word,option)
                call StringToUpper(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('NEIGHBOR_PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('VOLUME')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" scaled_volumetric_rate type'
                    call InputKeywordUnrecognized(input,word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in &
                  &flow condition "' // trim(condition%name) // &
                  '" scaled_volumetric_rate type'
                call printErrMsg(option)
              endif
            case('HETEROGENEOUS_VOLUMETRIC_RATE')
              sub_condition_ptr%itype = HET_VOL_RATE_SS
              rate_string = 'm^3/sec'
            case('HETEROGENEOUS_MASS_RATE')
              sub_condition_ptr%itype = HET_MASS_RATE_SS
              rate_string = 'kg/sec'
            case('HETEROGENEOUS_DIRICHLET')
              sub_condition_ptr%itype = HET_DIRICHLET_BC
            case('HETEROGENEOUS_SURFACE_SEEPAGE')
              sub_condition_ptr%itype = HET_SURF_HYDROSTATIC_SEEPAGE_BC
            case default
              call InputKeywordUnrecognized(input,word, &
                                            'flow condition,type',option)
          end select
        enddo
        call InputPopBlock(input,option)
      case('DATUM')
        dataset_ascii => DatasetAsciiCreate()
        call DatasetAsciiInit(dataset_ascii)
        dataset_ascii%array_width = 3
        dataset_ascii%data_type = DATASET_REAL
        condition%datum => dataset_ascii
        nullify(dataset_ascii)
        internal_units = 'meter'
        call ConditionReadValues(input,option,word,condition%datum, &
                                 word,internal_units)
      case('GRADIENT')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          call StringToUpper(word)
          select case(option%iflowmode)
            case(H_MODE)
              sub_condition_ptr => &
                FlowHydrateSubConditionPtr(input,word,hydrate,option)
          end select
          dataset_ascii => DatasetAsciiCreate()
          call DatasetAsciiInit(dataset_ascii)
          dataset_ascii%array_width = 3
          dataset_ascii%data_type = DATASET_REAL
          sub_condition_ptr%gradient => dataset_ascii
          nullify(dataset_ascii)
          internal_units = 'unitless/meter'
          call ConditionReadValues(input,option,word, &
                                   sub_condition_ptr%gradient, &
                                   word,internal_units)
          nullify(sub_condition_ptr)
        enddo
        call InputPopBlock(input,option)
      case('CONDUCTANCE')
        word = 'LIQUID_PRESSURE'
        select case(option%iflowmode)
          case(H_MODE)
            sub_condition_ptr => &
              FlowHydrateSubConditionPtr(input,word,hydrate,option)
        end select
        call InputReadDouble(input,option,sub_condition_ptr%aux_real(1))
        call InputErrorMsg(input,option,'LIQUID_CONDUCTANCE','CONDITION')
      case('LIQUID_PRESSURE','GAS_PRESSURE','LIQUID_SATURATION', &
           'ICE_SATURATION','GAS_SATURATION','HYDRATE_SATURATION', &
           'TEMPERATURE','MOLE_FRACTION','RATE','LIQUID_FLUX','GAS_FLUX', &
           'ENERGY_FLUX','RELATIVE_HUMIDITY')
        select case(option%iflowmode)
          case(H_MODE)
            sub_condition_ptr => &
              FlowHydrateSubConditionPtr(input,word,hydrate,option)
        end select
        internal_units = 'not_assigned'
        select case(trim(word))
          case('LIQUID_PRESSURE','GAS_PRESSURE')
            internal_units = 'Pa'
          case('LIQUID_SATURATION','GAS_SATURATION','HYDRATE_SATURATION', &
               'ICE_SATURATION','MOLE_FRACTION','RELATIVE_HUMIDITY')
            internal_units = 'unitless'
          case('TEMPERATURE')
            internal_units = 'C'
          case('RATE')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            select case(option%iflowmode)
              case(H_MODE)
                internal_units = trim(rate_string) // ',' // &
                  trim(rate_string) // ',MJ/sec|MW'
            end select
          case('LIQUID_FLUX','GAS_FLUX')
            internal_units = 'meter/sec'
          case('ENERGY_FLUX')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            internal_units = 'MW/m^2|MJ/m^2-sec'
        end select
        call ConditionReadValues(input,option,word, &
                                 sub_condition_ptr%dataset, &
                                 sub_condition_ptr%units,internal_units)
        input%force_units = PETSC_FALSE
        select case(word)
          case('LIQUID_SATURATION') ! convert to gas saturation
            if (associated(sub_condition_ptr%dataset%rbuffer)) then
              sub_condition_ptr%dataset%rbuffer(:) = 1.d0 - &
                sub_condition_ptr%dataset%rbuffer(:)
            endif
            sub_condition_ptr%dataset%rarray(:) = 1.d0 - &
              sub_condition_ptr%dataset%rarray(:)
        end select
      case default
        call InputKeywordUnrecognized(input,word,'flow condition',option)
    end select

  enddo
  call InputPopBlock(input,option)

  ! datum is not required
  string = 'SUBSURFACE/FLOW_CONDITION' // trim(condition%name) // '/Datum'
  call DatasetVerify(condition%datum,default_time_storage,string,option)

  ! need mole fraction and some sort of saturation
  if (.false.) then
    ! neumann or mass/volumetric flux
    ! need temperature
    if (.not.associated(hydrate%mole_fraction) .and. &
        .not.associated(hydrate%gas_saturation)) then
      option%io_buffer = trim(flow_mode_chars) // ' flux condition must &
        &include a MOLE_FRACTION or GAS/LIQUID_SATURATION.'
      call printErrMsg(option)
    endif
    if (associated(hydrate%mole_fraction) .and. &
        associated(hydrate%gas_saturation)) then
      option%io_buffer = trim(flow_mode_chars) // ' flux condition must &
        &include only a MOLE_FRACTION or GAS/LIQUID_SATURATION, not both.'
      call printErrMsg(option)
    endif
    if (.not.associated(hydrate%temperature)) then
      option%io_buffer = trim(flow_mode_chars) // ' flux condition must &
        &include a temperature.'
      call printErrMsg(option)
    endif
  else
    if (associated(hydrate%rate)) then
      condition%iphase = HYD_ANY_STATE
    elseif (associated(hydrate%liquid_flux) .and. &
            associated(hydrate%gas_flux) .and. &
            (associated(hydrate%energy_flux) .or. &
             associated(hydrate%temperature))) then
      condition%iphase = HYD_ANY_STATE
    else
      ! some sort of dirichlet-based pressure, temperature, etc.
      if (.not.associated(hydrate%liquid_pressure) .and. &
          .not.associated(hydrate%gas_pressure)) then
        option%io_buffer = 'Hydrate Mode non-rate condition must include &
          &a liquid or gas pressure'
        call printErrMsg(option)
      endif
      if (.not.associated(hydrate%mole_fraction) .and. &
          .not.associated(hydrate%relative_humidity) .and. &
          .not.associated(hydrate%gas_saturation)) then
        if (.not.associated(hydrate%hydrate_saturation) .and. &
                .not.associated(hydrate%ice_saturation)) then
          option%io_buffer = 'Hydrate Mode non-rate condition must &
                  &include a mole fraction, relative humidity, or &
                  &gas/liquid/hydrate/ice saturation'
          call printErrMsg(option)
        endif
      endif
      if (.not.associated(hydrate%temperature)) then
        option%io_buffer = 'Hydrate Mode non-rate condition must include &
          &a temperature, for now...'
        call printErrMsg(option)
      endif
      if ( associated(hydrate%gas_pressure) .and. &
           associated(hydrate%gas_saturation) .and. &
           associated(hydrate%liquid_pressure) .and. &
           (associated(hydrate%mole_fraction) .or. &
            associated(hydrate%relative_humidity)) ) then
        ! multiphase condition
        condition%iphase = HYD_MULTI_STATE
      elseif (associated(hydrate%liquid_pressure) .and. &
              associated(hydrate%mole_fraction) .and. &
              associated(hydrate%temperature)) then
        ! liquid phase condition
        condition%iphase = L_STATE
      elseif (associated(hydrate%gas_pressure) .and. &
               (associated(hydrate%mole_fraction) .or. &
                associated(hydrate%relative_humidity)) .and. &
                associated(hydrate%temperature)) then
        ! gas phase condition
        condition%iphase = G_STATE
      elseif (associated(hydrate%gas_pressure) .and. &
              associated(hydrate%gas_saturation) .and. &
              associated(hydrate%temperature)) then
        ! two phase condition
        condition%iphase = GA_STATE
      elseif ((associated(hydrate%gas_pressure) .or. &
          associated(hydrate%liquid_pressure)) .and. &
          associated(hydrate%hydrate_saturation)) then
        condition%iphase = HA_STATE
      elseif (associated(hydrate%gas_pressure) .and. &
              associated(hydrate%ice_saturation) .and. &
              associated(hydrate%temperature)) then
        condition%iphase = GI_STATE
      elseif (associated(hydrate%gas_pressure) .and. &
              associated(hydrate%mole_fraction) .and. &
              associated(hydrate%liquid_saturation)) then
        condition%iphase = AI_STATE
      elseif (associated(hydrate%gas_saturation) .and. &
              associated(hydrate%hydrate_saturation) .and. &
              associated(hydrate%temperature)) then
        condition%iphase = HGA_STATE
      elseif (associated(hydrate%gas_pressure) .and. &
              associated(hydrate%liquid_saturation) .and. &
              associated(hydrate%ice_saturation)) then
        condition%iphase = HAI_STATE
      elseif (associated(hydrate%ice_saturation) .and. &
              associated(hydrate%hydrate_saturation) .and. &
              associated(hydrate%temperature)) then
        condition%iphase = HGI_STATE
      elseif (associated(hydrate%liquid_saturation) .and. &
              associated(hydrate%gas_saturation) .and. &
              associated(hydrate%ice_saturation)) then
        condition%iphase = QUAD_STATE
      endif
      
    endif
    if (condition%iphase == NULL_STATE) then
      option%io_buffer = 'General Phase non-rate/flux condition contains &
        &an unsupported combination of primary dependent variables.'
      call printErrMsg(option)
    endif
  endif
  
   ! verify the datasets
  word = 'liquid pressure'
  call FlowSubConditionVerify(option,condition,word,hydrate%liquid_pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas pressure'
  call FlowSubConditionVerify(option,condition,word,hydrate%gas_pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas saturation'
  call FlowSubConditionVerify(option,condition,word,hydrate%gas_saturation, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'hydrate saturation'
  call FlowSubConditionVerify(option,condition,word,hydrate%hydrate_saturation,&
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'ice saturation'
  call FlowSubConditionVerify(option,condition,word,hydrate%ice_saturation, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'relative humidity'
  call FlowSubConditionVerify(option,condition,word,hydrate%relative_humidity,&
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'mole fraction'
  call FlowSubConditionVerify(option,condition,word,hydrate%mole_fraction, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,hydrate%temperature, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'liquid flux'
  call FlowSubConditionVerify(option,condition,word,hydrate%liquid_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas flux'
  call FlowSubConditionVerify(option,condition,word,hydrate%gas_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy flux'
  call FlowSubConditionVerify(option,condition,word,hydrate%energy_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,hydrate%rate, &
                              default_time_storage, &
                              PETSC_TRUE)

  condition%num_sub_conditions = 0
  i = 0
  if (associated(hydrate%liquid_pressure)) &
    i = i + 1
  if (associated(hydrate%gas_pressure)) &
    i = i + 1
  if (associated(hydrate%gas_saturation)) &
    i = i + 1
  if (associated(hydrate%hydrate_saturation)) &
    i = i + 1
  if (associated(hydrate%ice_saturation)) &
    i = i + 1
  if (associated(hydrate%relative_humidity)) &
    i = i + 1
  if (associated(hydrate%mole_fraction)) &
    i = i + 1
  if (associated(hydrate%temperature)) &
    i = i + 1
  if (associated(hydrate%liquid_flux)) &
    i = i + 1
  if (associated(hydrate%gas_flux)) &
    i = i + 1
  if (associated(hydrate%energy_flux)) &
    i = i + 1
  if (associated(hydrate%rate)) &
    i = i + 1
  condition%num_sub_conditions = i
  allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    nullify(condition%sub_condition_ptr(idof)%ptr)
  enddo
  i = 0
  if (associated(hydrate%liquid_pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%liquid_pressure
  endif
  if (associated(hydrate%gas_pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%gas_pressure
  endif
  if (associated(hydrate%gas_saturation)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%gas_saturation
  endif
  if (associated(hydrate%hydrate_saturation)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%hydrate_saturation
  endif
  if (associated(hydrate%ice_saturation)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%ice_saturation
  endif
  if (associated(hydrate%relative_humidity)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%relative_humidity
  endif
  if (associated(hydrate%mole_fraction)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%mole_fraction
  endif
  if (associated(hydrate%temperature)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%temperature
  endif
  if (associated(hydrate%liquid_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%liquid_flux
  endif
  if (associated(hydrate%gas_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%gas_flux
  endif
  if (associated(hydrate%energy_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%energy_flux
  endif
  if (associated(hydrate%rate)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => hydrate%rate
  endif

  ! set condition types
  allocate(condition%itype(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    condition%itype(idof) = condition%sub_condition_ptr(idof)%ptr%itype
  enddo

  condition%default_time_storage => default_time_storage

  call PetscLogEventEnd(logging%event_flow_condition_read,ierr);CHKERRQ(ierr)

end subroutine FlowConditionHydrateRead

! ************************************************************************** !

subroutine FlowConditionTOilImsRead(condition,input,option)
  !
  ! Reads a condition from the input file for
  ! toil_ims mode
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 9/9/2015
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Time_Storage_module
  use Dataset_module

  !use TOilIms_Aux_module
  use PM_TOilIms_Aux_module

  implicit none

  type(flow_condition_type) :: condition
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: rate_string, internal_units
  character(len=MAXWORDLENGTH) :: word, sub_word
  type(flow_toil_ims_condition_type), pointer :: toil_ims
  !type(flow_well_condition_type), pointer :: flow_well
  type(flow_sub_condition_type), pointer :: sub_condition_ptr
  PetscReal :: default_time
  PetscInt :: default_iphase
  PetscBool :: comm_card_found
  class(dataset_base_type), pointer :: default_flow_dataset
  class(dataset_base_type), pointer :: default_gradient
  PetscInt :: idof, i
  PetscBool :: default_is_cyclic
  type(time_storage_type), pointer :: default_time_storage
  class(dataset_ascii_type), pointer :: dataset_ascii
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read, &
                          ierr);CHKERRQ(ierr)

  rate_string = 'not_assigned'
  internal_units = 'not_assigned'

  default_time = 0.d0
  default_iphase = 0

  default_time_storage => TimeStorageCreate()
  default_time_storage%is_cyclic = PETSC_FALSE
  default_time_storage%time_interpolation_method = INTERPOLATION_STEP

  toil_ims => FlowTOilImsConditionCreate(option)
  condition%toil_ims => toil_ims

  ! read the condition
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    internal_units = 'not_assigned'

    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CONDITION')

    !reads cards common to all modes
    call FlowConditionCommonRead(condition,input,word,default_time_storage, &
                                 comm_card_found,option)
    if (comm_card_found) cycle

    select case(trim(word))

      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          call StringToUpper(word)

          select case(word)
            case('PRESSURE','OIL_PRESSURE','WATER_PRESSURE', &
                 'LIQUID_SATURATION', 'OIL_SATURATION','TEMPERATURE','RATE', &
                 'LIQUID_FLUX','OIL_FLUX', 'ENERGY_FLUX','ENTHALPY','OWC', &
                 'WATER_PRESSURE_GRAD')

              sub_condition_ptr => &
                FlowTOilImsSubConditionPtr(input,word,toil_ims,option)
            case default
              call InputKeywordUnrecognized(input,word,'flow condition',option)
          end select

          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'TYPE','CONDITION')
          call StringToUpper(word)

          sub_condition_ptr%ctype = word
          select case(word)
            case('DIRICHLET')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('NEUMANN')
              sub_condition_ptr%itype = NEUMANN_BC
            case('HYDROSTATIC')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('CONDUCTANCE')
              sub_condition_ptr%itype = HYDROSTATIC_CONDUCTANCE_BC
            case('SEEPAGE')
              sub_condition_ptr%itype = HYDROSTATIC_SEEPAGE_BC
            case('ZERO_GRADIENT')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case('MASS_RATE')
              sub_condition_ptr%itype = MASS_RATE_SS
              rate_string = 'kg/sec'
            case('SCALED_MASS_RATE')
              sub_condition_ptr%itype = SCALED_MASS_RATE_SS
              rate_string = 'kg/sec'
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call InputPushCard(input,word,option)
                call StringToUpper(word)
                sub_condition_ptr%ctype = &
                      trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('NEIGHBOR_PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('VOLUME')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                case('PERM')
                  sub_condition_ptr%isubtype = SCALE_BY_PERM
                case default
                  string = 'flow condition "' // trim(condition%name) // &
                    '" scaled_mass_rate type'
                  call InputKeywordUnrecognized(input,word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in flow condition "' // &
                  trim(condition%name) // &
                  '" scaled_mass_rate type'
                call PrintErrMsg(option)
              endif
            case('VOLUMETRIC_RATE')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'
            case('SCALED_VOLUMETRIC_RATE')
              sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call InputPushCard(input,word,option)
                call StringToUpper(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('NEIGHBOR_PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('VOLUME')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" scaled_volumetric_rate type'
                    call InputKeywordUnrecognized(input,word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in flow condition "' // &
                  trim(condition%name) // '" scaled_volumetric_rate type'
                call PrintErrMsg(option)
              endif
            case('HETEROGENEOUS_VOLUMETRIC_RATE')
              sub_condition_ptr%itype = HET_VOL_RATE_SS
              rate_string = 'm^3/sec'
            case('HETEROGENEOUS_MASS_RATE')
              sub_condition_ptr%itype = HET_MASS_RATE_SS
              rate_string = 'kg/sec'
            case('HETEROGENEOUS_DIRICHLET')
              sub_condition_ptr%itype = HET_DIRICHLET_BC
            case('HETEROGENEOUS_SURFACE_SEEPAGE')
              sub_condition_ptr%itype = HET_SURF_HYDROSTATIC_SEEPAGE_BC
            case('BHP')
              sub_condition_ptr%itype = WELL_BHP
            case('BHP_MIN')
              sub_condition_ptr%itype = WELL_BHP_MIN
            case('BHP_MAX')
              sub_condition_ptr%itype = WELL_BHP_MAX
            case default
              call InputKeywordUnrecognized(input,word, &
                                            'flow condition,type',option)
          end select
        enddo
        call InputPopBlock(input,option)
      
      case('GRADIENT','GRADIENT_D')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,sub_word)
          call InputErrorMsg(input,option,'keyword','GRADIENT,TYPE')
          call StringToUpper(sub_word)
          sub_condition_ptr => &
            FlowTOilImsSubConditionPtr(input,sub_word,toil_ims,option)
          dataset_ascii => DatasetAsciiCreate()
          call DatasetAsciiInit(dataset_ascii)
          dataset_ascii%array_width = 3
          dataset_ascii%data_type = DATASET_REAL
          sub_condition_ptr%gradient => dataset_ascii
          nullify(dataset_ascii)
          internal_units = 'unitless/meter'
          call ConditionReadValues(input,option,sub_word, &
                                   sub_condition_ptr%gradient, &
                                   sub_word, internal_units)
          select case(trim(word))
            case('GRADIENT_D')
              sub_condition_ptr%gradient%rarray(1:3) =  - &
                                       sub_condition_ptr%gradient%rarray(1:3)
          end select
          nullify(sub_condition_ptr)
        enddo
        call InputPopBlock(input,option)
      case('CONDUCTANCE')
        word = 'PRESSURE'
        select case(option%iflowmode)
          case(TOIL_IMS_MODE)
            sub_condition_ptr => &
              FlowTOilImsSubConditionPtr(input,word,toil_ims,option)
        end select
        call InputReadDouble(input,option,sub_condition_ptr%aux_real(1))
        call InputErrorMsg(input,option,'LIQUID_CONDUCTANCE','CONDITION')
  
      case('PRESSURE','OIL_PRESSURE','WATER_PRESSURE','LIQUID_SATURATION', &
           'OIL_SATURATION','TEMPERATURE','RATE', 'LIQUID_FLUX','OIL_FLUX', &
           'ENERGY_FLUX','ENTHALPY','WATER_PRESSURE_GRAD','OWC','OWC_Z', &
            'OWC_D','RTEMP','TEMPERATURE_AT_DATUM','PCOW_OWC')
        sub_condition_ptr => &
          FlowTOilImsSubConditionPtr(input,word,toil_ims,option)

        select case(trim(word))
        !give a type to pass FlowSubConditionVerify.
          case('OWC','WATER_PRESSURE_GRAD','OWC_Z','OWC_D','PCOW_OWC', &
               'RTEMP','TEMPERATURE_AT_DATUM')
            sub_condition_ptr%itype = DIRICHLET_BC
            sub_condition_ptr%ctype = 'dirichlet'
        end select

        internal_units = 'not_assigned'
        select case(trim(word))
          case('PRESSURE','OIL_PRESSURE','WATER_PRESSURE','PCOW_OWC')
            input%force_units = PETSC_TRUE
            internal_units = 'Pa'
          case('LIQUID_SATURATION','OIL_SATURATION')
            internal_units = 'unitless'
          case('TEMPERATURE','RTEMP','TEMPERATURE_AT_DATUM')
            input%force_units = PETSC_TRUE
            internal_units = 'C'
          case('OWC','OWC_Z','OWC_D')
            input%force_units = PETSC_TRUE
            internal_units = 'meter'
          case('WATER_PRESSURE_GRAD')
            internal_units = 'Pa/meter'
          case('RATE')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            internal_units = trim(rate_string) // ',' // trim(rate_string) //&
                             ',MJ/sec|MW'
          case('LIQUID_FLUX','OIL_FLUX')
            internal_units = 'meter/sec'
          case('ENERGY_FLUX')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            internal_units = 'MW/m^2|MJ/sec-m^2'
          case('ENTHALPY')
            input%force_units = PETSC_TRUE
            internal_units = 'J/kg'
            !internal_units = 'MJ/mol'
        end select
        call ConditionReadValues(input,option,word, &
                                 sub_condition_ptr%dataset, &
                                 sub_condition_ptr%units,internal_units)
        input%force_units = PETSC_FALSE
        select case(word)
          case('LIQUID_SATURATION') ! convert to oil saturation
            if (associated(sub_condition_ptr%dataset%rbuffer)) then
              sub_condition_ptr%dataset%rbuffer(:) = 1.d0 - &
                sub_condition_ptr%dataset%rbuffer(:)
            endif
            sub_condition_ptr%dataset%rarray(:) = 1.d0 - &
              sub_condition_ptr%dataset%rarray(:)
          case('OWC_D')
            sub_condition_ptr%dataset%rarray(:) = - &
                                     sub_condition_ptr%dataset%rarray(:)
        end select
      case default
        call InputKeywordUnrecognized(input,word,'flow condition',option)
    end select

  enddo
  call InputPopBlock(input,option)

  ! phase condition should never be used in TOilIms
  condition%iphase = ZERO_INTEGER

  ! unless the coondtion is a rate or pressure bhp (i.e. a bhp controlled well)
  ! - pressure is required
  ! - unless hydrostatic condition (sat and temp computed in hydrostatic equil):
  !    - water or oil saturation is required
  !    - temperature required (temp input checked in hydrostatic equil)
  if (.not.associated(toil_ims%rate)) then
    ! this branch is executed for sub_conditions that are not a rate or well
    ! some sort of dirichlet-based pressure, temperature, etc.
    if (.not.associated(toil_ims%pressure)) then
      option%io_buffer = 'TOilIms Phase non-rate condition must &
        &include a pressure'
      call PrintErrMsg(option)
    endif
    if ( toil_ims%pressure%itype /= HYDROSTATIC_BC ) then
      if (.not.associated(toil_ims%saturation) ) then
        option%io_buffer = 'TOilIms Phase non-rate condition must &
          &include liquid or oil saturation'
        call PrintErrMsg(option)
      endif
      if (.not.associated(toil_ims%temperature)) then
        option%io_buffer = 'TOilIms Phase non-rate condition must &
          &include temperature'
        call PrintErrMsg(option)
      endif
    end if 
  endif

  ! control that enthalpy is used for src/sink only
  if ( (.not.associated(toil_ims%rate)) .and. &
        associated(toil_ims%enthalpy)  ) then
      option%io_buffer = 'TOilIms Enthlapy condition is not &
       &currently supported for boundary & initial conditions'
      call PrintErrMsg(option)
  end if
  ! within a src/sink either temp or enthalpy can be defined
  if (associated(toil_ims%rate)) then
    if ( associated(toil_ims%temperature).and. &
        associated(toil_ims%enthalpy) &
       ) then
      option%io_buffer = 'TOilIms Rate condition can &
       &have either temp or enthalpy'
      call PrintErrMsg(option)
    end if
    ! only dirich condition supported for src/sink temp or enthalpy
    if ( associated(toil_ims%temperature) ) then
      if (toil_ims%temperature%itype /= DIRICHLET_BC) then
        option%io_buffer = 'TOilIms Src/Sink; only dirichlet type &
         &is supported for temperature conditions'
        call PrintErrMsg(option)
      end if
    end if

    if ( associated(toil_ims%enthalpy) ) then
      if (toil_ims%enthalpy%itype /= DIRICHLET_BC) then
        option%io_buffer = 'TOilIms Src/Sink; only dirichlet type &
         &is supported for enthalpy conditions'
        call PrintErrMsg(option)
      end if
    end if

  end if ! end if rate


  ! verify the datasets
  word = 'pressure'
  call FlowSubConditionVerify(option,condition,word,toil_ims%pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'oil saturation'
  call FlowSubConditionVerify(option,condition,word,toil_ims%saturation, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,toil_ims%temperature, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'enthalpy'
  call FlowSubConditionVerify(option,condition,word,toil_ims%enthalpy, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'liquid flux'
  call FlowSubConditionVerify(option,condition,word,toil_ims%liquid_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'oil flux'
  call FlowSubConditionVerify(option,condition,word,toil_ims%oil_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy flux'
  call FlowSubConditionVerify(option,condition,word,toil_ims%energy_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'oil water contact'
  call FlowSubConditionVerify(option,condition,word,toil_ims%owc, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'water pressure gradient'
  call FlowSubConditionVerify(option,condition,word,toil_ims%liq_press_grad, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,toil_ims%rate, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'datum_z'
  call FlowSubConditionVerify(option,condition,word,condition%datum_z, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'owc_z'
  call FlowSubConditionVerify(option,condition,word,toil_ims%owc_z, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'pcow_owc'
  call FlowSubConditionVerify(option,condition,word,toil_ims%pcow_owc, &
                              default_time_storage, &
                              PETSC_TRUE)

  condition%num_sub_conditions = 0
  i = 0
  if (associated(toil_ims%pressure)) &
    i = i + 1
  if (associated(toil_ims%saturation)) &
    i = i + 1
  if (associated(toil_ims%temperature)) &
    i = i + 1
  if (associated(toil_ims%enthalpy)) &
    i = i + 1
  if (associated(toil_ims%liquid_flux)) &
    i = i + 1
  if (associated(toil_ims%oil_flux)) &
    i = i + 1
  if (associated(toil_ims%energy_flux)) &
    i = i + 1
  if (associated(toil_ims%owc)) &
    i = i + 1
  if (associated(toil_ims%liq_press_grad)) &
    i = i + 1
  if (associated(toil_ims%rate)) &
    i = i + 1
  if (associated(toil_ims%owc_z)) &
    i = i + 1
  if (associated(toil_ims%pcow_owc)) &
    i = i + 1        

  ! assing number of sub_condition
  condition%num_sub_conditions = i
  allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    nullify(condition%sub_condition_ptr(idof)%ptr)
  enddo
  i = 0
  if (associated(toil_ims%pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%pressure
  endif
  if (associated(toil_ims%saturation)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%saturation
  endif
  if (associated(toil_ims%temperature)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%temperature
  endif
  if (associated(toil_ims%enthalpy)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%enthalpy
  endif
  if (associated(toil_ims%liquid_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%liquid_flux
  endif
  if (associated(toil_ims%oil_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%oil_flux
  endif
  if (associated(toil_ims%energy_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%energy_flux
  endif
  if (associated(toil_ims%owc)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%owc
  endif
  if (associated(toil_ims%liq_press_grad)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%liq_press_grad
  endif
  if (associated(toil_ims%rate)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%rate
  endif
  if (associated(toil_ims%owc_z)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%owc_z
  end if  
  if (associated(toil_ims%pcow_owc)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%pcow_owc
  end if

  ! set condition types
  allocate(condition%itype(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    condition%itype(idof) = condition%sub_condition_ptr(idof)%ptr%itype
  enddo

  condition%default_time_storage => default_time_storage

  call PetscLogEventEnd(logging%event_flow_condition_read,ierr);CHKERRQ(ierr)

end subroutine FlowConditionTOilImsRead

! ************************************************************************** !

subroutine FlowConditionTOWGRead(condition,input,option)
  !
  ! Reads a condition from the input file for TOWG mode
  !
  ! Note: consider refactoring to avoid repeating code. Modularise common
  !       contents of: FlowConditionGeneralRead, FlowConditionTOilImsRead
  !                    and FlowConditionTOWGRead
  !       And refactor to use these modules within FlowConditionPMRead
  !
  ! Author: Paolo Orsini
  ! Date: 10/16/16
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Time_Storage_module
  use Dataset_module

  use PM_TOWG_Aux_module

  implicit none

  type(flow_condition_type) :: condition
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: rate_string
  character(len=MAXSTRINGLENGTH) :: internal_units_string
  character(len=MAXWORDLENGTH) :: word, sub_word
  character(len=MAXWORDLENGTH) :: internal_units_word
  character(len=MAXWORDLENGTH) :: usr_tbl_len_units, usr_tbl_press_units
  type(flow_towg_condition_type), pointer :: towg
  type(flow_sub_condition_type), pointer :: sub_condition_ptr
  PetscReal :: default_time
  PetscInt :: default_iphase
  class(dataset_base_type), pointer :: default_flow_dataset
  class(dataset_base_type), pointer :: default_gradient
  PetscInt :: idof, i
  PetscInt :: data_idx
  PetscBool :: default_is_cyclic
  PetscBool :: phase_state_found
  PetscBool :: comm_card_found
  PetscBool :: usr_tbl_press_units_found
  PetscBool :: usr_tbl_z_units_found
  PetscBool :: pbvz_found
  type(time_storage_type), pointer :: default_time_storage
  class(dataset_ascii_type), pointer :: dataset_ascii
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read, &
                          ierr);CHKERRQ(ierr)

  default_time = 0.d0
  default_iphase = 0

  default_time_storage => TimeStorageCreate()
  default_time_storage%is_cyclic = PETSC_FALSE
  default_time_storage%time_interpolation_method = INTERPOLATION_STEP

  towg => FlowTOWGConditionCreate(option)
  condition%towg => towg

  rate_string = 'not_assigned'

  ! read the condition
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    internal_units_string = 'not_assigned'

    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CONDITION')

    !reads cards common to all modes
    call FlowConditionCommonRead(condition,input,word,default_time_storage, &
                                 comm_card_found,option)
    if (comm_card_found) cycle
    
    select case(trim(word))

      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          call StringToUpper(word)
          sub_condition_ptr => &
            FlowTOWGSubConditionPtr(input,word,towg,option)
          !when refactoring
          !sub_condition_ptr => FlowPMSubConditionPtr(word,condition,option)
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'TYPE','CONDITION')
          call StringToUpper(word)
          sub_condition_ptr%ctype = word
          select case(word)
            case('DIRICHLET')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('NEUMANN')
              sub_condition_ptr%itype = NEUMANN_BC
            case('HYDROSTATIC')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('CONDUCTANCE')
              sub_condition_ptr%itype = HYDROSTATIC_CONDUCTANCE_BC
            case('SEEPAGE')
              sub_condition_ptr%itype = HYDROSTATIC_SEEPAGE_BC
            case('MASS_RATE')
              sub_condition_ptr%itype = MASS_RATE_SS
              rate_string = 'kg/sec'
            case('SCALED_MASS_RATE')
              sub_condition_ptr%itype = SCALED_MASS_RATE_SS
              rate_string = 'kg/sec'
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call InputPushCard(input,word,option)
                call StringToUpper(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('NEIGHBOR_PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('VOLUME')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" scaled_mass_rate type'
                    call InputKeywordUnrecognized(input,word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in flow condition "' // &
                  trim(condition%name) // &
                  '" scaled_mass_rate type'
                call PrintErrMsg(option)
              endif
            case('VOLUMETRIC_RATE')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'
            case('SCALED_VOLUMETRIC_RATE')
              sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call InputPushCard(input,word,option)
                call StringToUpper(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('NEIGHBOR_PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('VOLUME')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('PERM')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" scaled_volumetric_rate type'
                    call InputKeywordUnrecognized(input,word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in flow condition "' // &
                  trim(condition%name) // &
                  '" scaled_volumetric_rate type'
                call PrintErrMsg(option)
              endif
            case('HETEROGENEOUS_VOLUMETRIC_RATE')
              sub_condition_ptr%itype = HET_VOL_RATE_SS
              rate_string = 'm^3/sec'
            case('HETEROGENEOUS_MASS_RATE')
              sub_condition_ptr%itype = HET_MASS_RATE_SS
              rate_string = 'kg/sec'
            case('HETEROGENEOUS_DIRICHLET')
              sub_condition_ptr%itype = HET_DIRICHLET_BC
            case default
              call InputKeywordUnrecognized(input,word, &
                                            'flow condition,type',option)
          end select
        enddo
        call InputPopBlock(input,option)

      case('GRADIENT','GRADIENT_D')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,sub_word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          call StringToUpper(sub_word)
          sub_condition_ptr => &
            FlowTOWGSubConditionPtr(input,sub_word,towg,option)
          dataset_ascii => DatasetAsciiCreate()
          call DatasetAsciiInit(dataset_ascii)
          dataset_ascii%array_width = 3
          dataset_ascii%data_type = DATASET_REAL
          sub_condition_ptr%gradient => dataset_ascii
          nullify(dataset_ascii)
          internal_units_string = 'unitless/meter'
          call ConditionReadValues(input,option,sub_word, &
                                   sub_condition_ptr%gradient, &
                                   sub_word,internal_units_string)
          select case(trim(word))
            case ('GRADIENT_D')
              sub_condition_ptr%gradient%rarray(1:3) = - &
                                 sub_condition_ptr%gradient%rarray(1:3)
          end select  
          nullify(sub_condition_ptr)
        enddo
        call InputPopBlock(input,option)
      case('CONDUCTANCE')
        word = 'LIQUID_PRESSURE'
        sub_condition_ptr => &
          FlowTOWGSubConditionPtr(input,word,towg,option)
        call InputReadDouble(input,option,sub_condition_ptr%aux_real(1))
        call InputErrorMsg(input,option,'LIQUID_CONDUCTANCE','CONDITION')
      case('WATER_GAS_EQUILIBRATION')
        towg%is_wg_equilibration = PETSC_TRUE
      case('PRESSURE','OIL_PRESSURE','GAS_PRESSURE','OIL_SATURATION', &
          'GAS_SATURATION','SOLVENT_SATURATION','BUBBLE_POINT','TEMPERATURE', &
          'RATE','BHP_PRESSURE', 'LIQUID_FLUX','OIL_FLUX','GAS_FLUX', &
          'SOLVENT_FLUX','ENERGY_FLUX','ENTHALPY', &
          'OWC_Z','OWC_D','PCOW_OWC', 'OGC_Z','OGC_D', 'PCOG_OGC', &
          'WGC_Z','WGC_D','PCWG_WGC','RTEMP','TEMPERATURE_AT_DATUM')
        sub_condition_ptr => &
          FlowTOWGSubConditionPtr(input,word,towg,option)
        select case(trim(word))
          case('PRESSURE','OIL_PRESSURE','GAS_PRESSURE','BHP_PRESSURE', &
               'BUBBLE_POINT','PCOW_OWC','PCOG_OGC','PCWG_WGC')
            internal_units_string = 'Pa'
            input%force_units = PETSC_TRUE
          case('OWC_Z','OWC_D','OGC_Z','OGC_D','WGC_Z','WGC_D', &
               'DATUM_Z','DATUM_D')
            internal_units_string = 'meter'
            input%force_units = PETSC_TRUE 
          case('OIL_SATURATION','GAS_SATURATION','SOLVENT_SATURATION')
            internal_units_string = 'unitless'
          case('TEMPERATURE','RTEMP','TEMPERATURE_AT_DATUM')
            internal_units_string = 'C'
            input%force_units = PETSC_TRUE
          case('RATE')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            if (towg_miscibility_model == TOWG_SOLVENT_TL) then
              internal_units_string= trim(rate_string) // ',' // &
                                     trim(rate_string) // ',' // &
                                     trim(rate_string) // ',' // &
                                     trim(rate_string) // ',MJ/sec|MW'
            else
              internal_units_string = trim(rate_string) // ',' // &
                                      trim(rate_string) // ',' // &
                                      trim(rate_string) // ',MJ/sec|MW'
            end if
          case('ENTHALPY')
            input%force_units = PETSC_TRUE
            if (towg_miscibility_model == TOWG_SOLVENT_TL) then
              internal_units_string = 'J/kg' // ',' // 'J/kg' // ',' // &
                                      'J/kg' // ',' // 'J/kg'
            else
              internal_units_string = 'J/kg' // ',' // 'J/kg' // ',' // 'J/kg'
            end if   
          case('LIQUID_FLUX','OIL_FLUX','GAS_FLUX')
            internal_units_string = 'meter/sec'
            input%force_units = PETSC_TRUE
          case('ENERGY_FLUX')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            internal_units_string = 'MW/m^2|MJ/m^2-sec'
        end select
        select case(trim(word))
        !give a type to pass FlowSubConditionVerify.
          case('OWC_Z','OWC_D','OGC_Z','OGC_D','PCOW_OWC','PCOG_OGC', &
                'PCWG_WGC','RTEMP','TEMPERATURE_AT_DATUM')
            sub_condition_ptr%itype = DIRICHLET_BC
            sub_condition_ptr%ctype = 'dirichlet'
        end select
        call ConditionReadValues(input,option,word, &
                            sub_condition_ptr%dataset, &
                            sub_condition_ptr%units,internal_units_string)
        input%force_units = PETSC_FALSE
        select case(trim(word))
          case('OWC_D','OGC_D','WGC_D')
            sub_condition_ptr%dataset%rarray(:) =  - &
                                    sub_condition_ptr%dataset%rarray(:)
        end select
      !PO move BUBBLE_POINT_TABLE to new routine to improve 
      !  FlowConditionTOWGRead readibility
      case('BUBBLE_POINT_TABLE')
        towg%pbvz_table => LookupTableCreateGeneral(ONE_INTEGER)
        call towg%pbvz_table%LookupTableVarsInit(TWO_INTEGER)
        !set up default units
        internal_units_word = 'not_assigned'
        usr_tbl_len_units = 'm'
        usr_tbl_press_units = 'Pa'
        usr_tbl_press_units_found = PETSC_FALSE
        usr_tbl_z_units_found = PETSC_FALSE
        pbvz_found = PETSC_FALSE
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option, &
                                                'CONDITION,BUBBLE_POINT_TABLE')
          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword', &
                                               'CONDITION,BUBBLE_POINT_TABLE')
          call StringToUpper(word)
          select case (trim(word))
            case('Z_UNITS','D_UNITS')
              call InputReadWord(input,option,usr_tbl_len_units,PETSC_TRUE)
              call InputErrorMsg(input,option,'BUBBLE_POINT_TABLE','Z/D_UNITS')
              usr_tbl_z_units_found = PETSC_TRUE
            case('PRESSURE_UNITS')
              call InputReadWord(input,option,usr_tbl_press_units,PETSC_TRUE)
              call InputErrorMsg(input,option,'BUBBLE_POINT_TABLE', &
                                                            'PRESSURE_UNITS')
              usr_tbl_press_units_found = PETSC_TRUE
            case('PBVZ','PBVD')
              data_idx = 1 !elevation/depth in the first column
              internal_units_word = 'm'
              call towg%pbvz_table%CreateAddLookupTableVar(ONE_INTEGER, &
                        internal_units_word,usr_tbl_len_units,data_idx,option)
              data_idx = 2 !Bubble point pressure in the second column
              internal_units_word = 'Pa'
              call towg%pbvz_table%CreateAddLookupTableVar(TWO_INTEGER, &
                      internal_units_word,usr_tbl_press_units,data_idx,option)
              string = 'SUBSURFACE/FLOW_CONDITION' // trim(condition%name)  &
                                         // '/reading PBVZ or PBVD table'
              call towg%pbvz_table%VarDataRead(input,TWO_INTEGER,TWO_INTEGER, &
                                               string,option)
              !table input for decreasing z (or decreasing depth) - needs reversing
              call towg%pbvz_table%VarDataReverse(option)
              ! PO todo: move this error to table general reads
              if ( size(towg%pbvz_table%var_data(1,:)) < 2 ) then
                option%io_buffer = 'PBVZ, PBVD tables require at least two &
                                    &entries'
                call PrintErrMsg(option)
              end if

              select case(trim(word))
                case('PBVD')
                  !convert depth in elevation z for internal PFLOTRAN use
                  towg%pbvz_table%var_data(ONE_INTEGER,:) = - &
                                       towg%pbvz_table%var_data(ONE_INTEGER,:)
              end select

              pbvz_found = PETSC_TRUE
            case default
              call InputKeywordUnrecognized(input,word, &
                                  'flow condition,BUBBLE_POINT_TABLE',option)
          end select
        end do
        call InputPopBlock(input,option)
        if ( .not. usr_tbl_z_units_found ) then
          option%io_buffer = 'TOWG condition - BUBBLE_POINT_TABLE: &
            &lenght units must be entered for the z/depths'
          call PrintErrMsg(option)
        else if ( usr_tbl_z_units_found .and. pbvz_found ) then
          call towg%pbvz_table%SetupVarUserUnits(ONE_INTEGER, &
                                               usr_tbl_len_units,option)
        end if
        if ( .not. usr_tbl_press_units_found) then
          option%io_buffer = 'TOWG condition - BUBBLE_POINT_TABLE: &
            &pressure units must be entered for the bubble point'
          call PrintErrMsg(option)
        else if ( usr_tbl_press_units_found .and. pbvz_found ) then
          call towg%pbvz_table%SetupVarUserUnits(TWO_INTEGER, &
                                              usr_tbl_press_units,option)
        end if
        ! LookupTable unit conversion after reading units to make the table input
        ! independent from the order the table instructions are given
        if (pbvz_found) then
          call towg%pbvz_table%LookupTableVarConvFactors(option)
          call towg%pbvz_table%VarPointAndUnitConv(option)
          call towg%pbvz_table%SetupConstValExtrap(option)
          call towg%pbvz_table%LookupTableVarInitGradients(option)
          ! define table axis
          !PO: to move into table%SetUpIndependentVars(var1,var2), var 2 is optional
          allocate(towg%pbvz_table%axis1)
          allocate(towg%pbvz_table%axis1%values( &
                                    size(towg%pbvz_table%var_data(1,:))))
          towg%pbvz_table%axis1%values = towg%pbvz_table%var_data(1,:)
          towg%pbvz_table%dims(1) = size(towg%pbvz_table%axis1%values(:))
        else
          option%io_buffer = 'TOWG condition - BUBBLE_POINT_TABLE: &
                              &defined but pressure table (PBVZ or PBVD) &
                              &not found.If you are trying to run without &
                              &PBVD table, please remove the entire &
                              BUBBLE_POINT_TABLE card'
          call PrintErrMsg(option)
        end if
      case default
        call InputKeywordUnrecognized(input,word,'flow condition',option)
    end select

  enddo
  call InputPopBlock(input,option)

  !initialise phase_state to null
  condition%iphase = TOWG_NULL_STATE
  !assign phase state
  if (associated(towg%rate)) then
    condition%iphase = TOWG_ANY_STATE
  !Neumann condition for all phases - IMMISCIBLE, LONGOSTAFF, BLACK_OIL
  else if (associated(towg%liquid_flux) .and. &
           associated(towg%oil_flux) .and. &
           associated(towg%gas_flux) .and. &
           (associated(towg%energy_flux) .or. &
            associated(towg%temperature)) &
          ) then
    condition%iphase = TOWG_ANY_STATE
  !Neumann condition for all phases - SOLVENT_TL
  else if (associated(towg%liquid_flux) .and. &
           associated(towg%oil_flux) .and. &
           associated(towg%gas_flux) .and. &
           associated(towg%solvent_flux) .and. &
           (associated(towg%energy_flux) .or. &
            associated(towg%temperature)) &
          ) then
    condition%iphase = TOWG_ANY_STATE
  else
    !check oil/gas phase state
    phase_state_found=PETSC_FALSE
    if ( associated(towg%oil_pressure  ) .and. &
         associated(towg%oil_saturation) ) then
      if (     (towg_miscibility_model == TOWG_BLACK_OIL ) &
           .or.(towg_miscibility_model == TOWG_SOLVENT_TL) ) then
! Setup for black oil case - needs gas saturation and/or bubble point
        if(      associated(towg%gas_saturation) &
           .and. associated(towg%bubble_point  ) ) then
! Both found - can be saturated or undersaturated - i.e. any state
          condition%iphase  = TOWG_ANY_STATE
          phase_state_found = PETSC_TRUE
        else if( associated(towg%gas_saturation) ) then
! Only gas saturation only - is saturated three-phase
          condition%iphase  = TOWG_THREE_PHASE_STATE
          phase_state_found = PETSC_TRUE
        else if( associated(towg%bubble_point  ) ) then
! Only bubble point found - is undersaturated water-oil
          condition%iphase  = TOWG_LIQ_OIL_STATE
          phase_state_found = PETSC_TRUE
        end if
      else
! Setup for all but black oil case - just needs gas saturation
        if( associated(towg%gas_saturation) ) then
          condition%iphase  = TOWG_THREE_PHASE_STATE
          phase_state_found = PETSC_TRUE
        end if
      end if
    !only oil pressuse associated: this is a hydrostatic equilibration
    !initialise to TOWG_ANY_STATE. Phase state computed in HydrostaticMPUpdateCoupler
    else if ( associated(towg%oil_pressure) .and. &
              towg%oil_pressure%itype == HYDROSTATIC_BC ) then
       condition%iphase  = TOWG_ANY_STATE
       phase_state_found = PETSC_TRUE
    end if

    !check that a valid phase state has been found
    if( .not.phase_state_found ) then
      option%io_buffer = 'TOWG condition - phase state  &
        &Currently only THREE_PHASE_STATE, LIQ_OIL_STATE and ANY_STATE implemented'
      call PrintErrMsg(option)
    endif

    !check if conditions are compatible with miscibility model
    if ( (towg_miscibility_model == TOWG_IMMISCIBLE .or. &
          towg_miscibility_model == TOWG_TODD_LONGSTAFF) .and. &
          condition%iphase /= TOWG_THREE_PHASE_STATE .and. &
          .not.( (associated(towg%oil_pressure) .and. &
                 towg%oil_pressure%itype == HYDROSTATIC_BC ) )  &
       ) then
      option%io_buffer = 'FlowConditionTOWGRead: For TOWG_IMMISCIBLE and &
         &TOWG_TODD_LONGSTAFF only three phase state conditions &
         &are supported other than rate and flux conditions '
      call PrintErrMsg(option)
    end if

  end if

  ! control that enthalpy is used for src/sink only
  if ( (.not.associated(towg%rate)) .and. &
        associated(towg%enthalpy)  ) then
      option%io_buffer = 'TOWG Enthlapy condition is not &
       &currently supported for boundary & initial conditions'
      call PrintErrMsg(option)
  end if
  ! within a src/sink either temp or enthalpy can be defined
  if (associated(towg%rate)) then
    if ( associated(towg%temperature).and. &
        associated(towg%enthalpy) &
       ) then
      option%io_buffer = 'TOWG Rate condition can &
       &have either temp or enthalpy'
      call PrintErrMsg(option)
    end if
    ! only dirich condition supported for src/sink temp or enthalpy
    if (associated(towg%temperature)) then
      if (towg%temperature%itype /= DIRICHLET_BC) then
         option%io_buffer = 'TOWG Src/Sink; only dirichlet type &
                             &is supported for temperature conditions'
         call PrintErrMsg(option)
      end if
    end if
    if (associated(towg%enthalpy)) then
      if (towg%enthalpy%itype /= DIRICHLET_BC ) then
        option%io_buffer = 'TOWG Src/Sink; only dirichlet type &
                            &is supported for enthalpy conditions'
        call PrintErrMsg(option)
      end if
    end if
  end if ! end if rate

  if (condition%iphase == TOWG_NULL_STATE) then
    option%io_buffer = 'TOWG Phase non-rate/flux condition contains &
      &an unsupported combination of primary dependent variables.'
    call PrintErrMsg(option)
  endif

  ! verify the datasets
  word = 'oil pressure'
  call FlowSubConditionVerify(option,condition,word,towg%oil_pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas pressure'
  call FlowSubConditionVerify(option,condition,word,towg%gas_pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'oil saturation'
  call FlowSubConditionVerify(option,condition,word,towg%oil_saturation, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas saturation'
  call FlowSubConditionVerify(option,condition,word,towg%gas_saturation, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'solvent saturation'
  call FlowSubConditionVerify(option,condition,word,towg%solvent_saturation, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'bubble point'
  call FlowSubConditionVerify(option,condition,word,towg%bubble_point, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,towg%temperature, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'enthalpy'
  call FlowSubConditionVerify(option,condition,word,towg%enthalpy, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'liquid flux'
  call FlowSubConditionVerify(option,condition,word,towg%liquid_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'oil flux'
  call FlowSubConditionVerify(option,condition,word,towg%oil_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas flux'
  call FlowSubConditionVerify(option,condition,word,towg%gas_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'solvent flux'
  call FlowSubConditionVerify(option,condition,word,towg%solvent_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy flux'
  call FlowSubConditionVerify(option,condition,word,towg%energy_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,towg%rate, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'bhp pressure'
  call FlowSubConditionVerify(option,condition,word,towg%bhp_pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'owc_z'
  call FlowSubConditionVerify(option,condition,word,towg%owc_z, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'pcow_owc'
  call FlowSubConditionVerify(option,condition,word,towg%pcow_owc, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'ogc_z'
  call FlowSubConditionVerify(option,condition,word,towg%ogc_z, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'pcog_ogc'
  call FlowSubConditionVerify(option,condition,word,towg%pcog_ogc, &
                              default_time_storage, &
                              PETSC_TRUE)

  condition%num_sub_conditions = 0
  i = 0
  if (associated(towg%oil_pressure)) &
    i = i + 1
  if (associated(towg%gas_pressure)) &
    i = i + 1
  if (associated(towg%oil_saturation)) &
    i = i + 1
  if (associated(towg%gas_saturation)) &
    i = i + 1
  if (associated(towg%solvent_saturation)) &
    i = i + 1
  if (associated(towg%bubble_point)) &
    i = i + 1
  if (associated(towg%temperature)) &
    i = i + 1
  if (associated(towg%enthalpy)) &
    i = i + 1
  if (associated(towg%liquid_flux)) &
    i = i + 1
  if (associated(towg%oil_flux)) &
    i = i + 1
  if (associated(towg%gas_flux)) &
    i = i + 1
  if (associated(towg%solvent_flux)) &
    i = i + 1
  if (associated(towg%energy_flux)) &
    i = i + 1
  if (associated(towg%rate)) &
    i = i + 1
  if (associated(towg%bhp_pressure)) &
    i = i + 1
  if (associated(towg%owc_z)) &
    i = i + 1
  if (associated(towg%pcow_owc)) &
    i = i + 1
  if (associated(towg%ogc_z)) &
    i = i + 1
  if (associated(towg%pcog_ogc)) &
    i = i + 1

  condition%num_sub_conditions = i
  allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    nullify(condition%sub_condition_ptr(idof)%ptr)
  enddo
  i = 0
  if (associated(towg%oil_pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%oil_pressure
  endif
  if (associated(towg%gas_pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%gas_pressure
  endif
  if (associated(towg%oil_saturation)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%oil_saturation
  endif
  if (associated(towg%gas_saturation)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%gas_saturation
  endif
  if (associated(towg%solvent_saturation)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%solvent_saturation
  endif
  if (associated(towg%bubble_point)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%bubble_point
  endif
  if (associated(towg%temperature)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%temperature
  endif
  if (associated(towg%enthalpy)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%enthalpy
  endif
  if (associated(towg%liquid_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%liquid_flux
  endif
  if (associated(towg%oil_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%oil_flux
  endif
  if (associated(towg%gas_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%gas_flux
  endif
  if (associated(towg%solvent_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%solvent_flux
  endif
  if (associated(towg%energy_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%energy_flux
  endif
  if (associated(towg%rate)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%rate
  endif
  if (associated(towg%bhp_pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%bhp_pressure
  endif
  if (associated(towg%owc_z)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%owc_z
  endif
  if (associated(towg%pcow_owc)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%pcow_owc
  endif
  if (associated(towg%ogc_z)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%ogc_z
  endif
  if (associated(towg%pcog_ogc)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => towg%pcog_ogc
  endif

  ! set condition types
  allocate(condition%itype(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    condition%itype(idof) = condition%sub_condition_ptr(idof)%ptr%itype
  enddo

  condition%default_time_storage => default_time_storage

  call PetscLogEventEnd(logging%event_flow_condition_read,ierr);CHKERRQ(ierr)

end subroutine FlowConditionTOWGRead

! ************************************************************************** !

subroutine FlowConditionCommonRead(condition,input,word,default_time_storage, &
                                   card_found,option)
  !
  ! Reads flow conditions block common to all modes:
  !  these are sflow_sub_conditions defined in flow_condition_type
  !  for each commom card must add card_found = PETSC_TRUE
  !
  ! Author: Paolo Orsini
  ! Date: 01/17/19
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Time_Storage_module
  use Dataset_module
  use Lookup_Table_module

  implicit none


  type(flow_condition_type) :: condition
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: word
  type(time_storage_type), pointer :: default_time_storage
  PetscBool, intent(out) :: card_found
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: internal_units_string
  character(len=MAXWORDLENGTH) :: internal_units_word
  character(len=MAXWORDLENGTH) :: usr_lenght_units
  character(len=MAXWORDLENGTH) :: usr_temp_units
  PetscInt :: data_idx
  PetscBool :: rtempvz_found
  PetscBool :: rtempvz_z_units_found
  PetscBool :: rtempvz_temp_units_found

  class(dataset_ascii_type), pointer :: dataset_ascii
!  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: sub_word
  class(lookup_table_general_type), pointer :: lkp_table => null()
  !type(lookup_table_var_type), pointer :: lkp_var => null()


  card_found = PETSC_FALSE
  select case (trim(word))

    case('CYCLIC')
      card_found = PETSC_TRUE
      ! by default, is_cyclic is set to PETSC_FALSE
      default_time_storage%is_cyclic = PETSC_TRUE
      
    case('SYNC_TIMESTEP_WITH_UPDATE')
      card_found = PETSC_TRUE
      condition%sync_time_with_update = PETSC_TRUE
      
    case('INTERPOLATION')
      card_found = PETSC_TRUE
      call InputReadCard(input,option,word)
      call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')
      call StringToUpper(word)
      select case(word)
        case('STEP')
          default_time_storage%time_interpolation_method = &
            INTERPOLATION_STEP
        case('LINEAR')
          default_time_storage%time_interpolation_method = &
            INTERPOLATION_LINEAR
      end select
      
    case('DATUM')
      card_found = PETSC_TRUE
      dataset_ascii => DatasetAsciiCreate()
      call DatasetAsciiInit(dataset_ascii)
      dataset_ascii%array_width = 3
      dataset_ascii%data_type = DATASET_REAL
      condition%datum => dataset_ascii
      nullify(dataset_ascii)
      internal_units_string = 'meter'
      call ConditionReadValues(input,option,word,condition%datum, &
                               word,internal_units_string)
      !
      string = 'SUBSURFACE/FLOW_CONDITION' // trim(condition%name) // '/Datum'
      call DatasetVerify(condition%datum,default_time_storage,string,option)      
    
    case('DATUM_Z','DATUM_D')
      card_found = PETSC_TRUE
      condition%datum_z => FlowSubConditionCreate(ONE_INTEGER)
      internal_units_string = 'meter'
      input%force_units = PETSC_TRUE
      call ConditionReadValues(input,option,word, &
                               condition%datum_z%dataset, &
                               condition%datum_z%units,internal_units_string)
      input%force_units = PETSC_FALSE
      !give a condition type to pass verify
      condition%datum_z%itype = DIRICHLET_BC
      condition%datum_z%ctype = 'dirichlet'
      select case(word)
        case('DATUM_D')
          condition%datum_z%dataset%rarray(:) = - &
                                        condition%datum_z%dataset%rarray(:)
      end select
      call FlowSubConditionVerify(option,condition,word,condition%datum_z, &
                                  default_time_storage, &
                                  PETSC_TRUE)
    case ('TEMPERATURE_TABLE')
      card_found = PETSC_TRUE
      rtempvz_found = PETSC_FALSE
      rtempvz_z_units_found = PETSC_FALSE
      rtempvz_temp_units_found = PETSC_FALSE
      condition%rtempvz_table => LookupTableCreateGeneral(ONE_INTEGER)
      call condition%rtempvz_table%LookupTableVarsInit(TWO_INTEGER)
      lkp_table => condition%rtempvz_table
      !set up default units
      usr_lenght_units = 'm'
      usr_temp_units = 'C'
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option, &
                                               'CONDITION,TEMPERATURE_TABLE')
        if (InputCheckExit(input,option)) exit

        if (InputError(input)) exit
        call InputReadCard(input,option,sub_word,PETSC_FALSE)
        call InputErrorMsg(input,option,'keyword', &
                                                'CONDITION,TEMPERATURE_TABLE')
        call StringToUpper(sub_word)
        select case (trim(sub_word))
          case('Z_UNITS','D_UNITS')
            call InputReadWord(input,option,usr_lenght_units,PETSC_TRUE)
            call InputErrorMsg(input,option,'TEMPERATURE_TABLE','Z/D_UNITS')
            rtempvz_z_units_found = PETSC_TRUE
          case('TEMPERATURE_UNITS')
            call InputReadWord(input,option,usr_temp_units,PETSC_TRUE)
            call InputErrorMsg(input,option,'TEMPERATURE_TABLE', &
                                                          'TEMPERATURE_UNITS')
            if ( trim(usr_temp_units) /= 'C' ) then
              option%io_buffer = 'TEMPERATURE_TABLE supports only &
                                  &degree celcius, C' 
               call PrintErrMsg(option)
            end if
            rtempvz_temp_units_found = PETSC_TRUE
          case('RTEMPVZ','RTEMPVD')
            data_idx = 1 !elevation/depth in the first column
            internal_units_word = 'm'
            call lkp_table%CreateAddLookupTableVar(ONE_INTEGER, &
                         internal_units_word,usr_lenght_units,data_idx,option)
            data_idx = 2 !temperature in the second column
            internal_units_word = 'C'
            call lkp_table%CreateAddLookupTableVar(TWO_INTEGER, &
                            internal_units_word,usr_temp_units,data_idx,option)
            string = 'SUBSURFACE/FLOW_CONDITION' // trim(condition%name)  &
                      // '/reading RTEMPVZ or RTEMPVD table'
            call lkp_table%VarDataRead(input,TWO_INTEGER,TWO_INTEGER, &
                                       string,option)
            !table input for decreasing z (or decreasing depth) - needs reversing
            call lkp_table%VarDataReverse(option)
            ! PO todo: move this error to table general reads
            if ( size(lkp_table%var_data(1,:)) < 2 ) then
              option%io_buffer = 'RTEMPVZ, RTEMPVD tables require at least two &
                                  &entries'
              call PrintErrMsg(option)
            end if

            select case(trim(sub_word))
              case('RTEMPVD')
                !convert depth in elevation z for internal PFLOTRAN use
                lkp_table%var_data(ONE_INTEGER,:) = - &
                                         lkp_table%var_data(ONE_INTEGER,:)
            end select

            rtempvz_found = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,sub_word, &
                                    'flow condition,TEMPERATURE_TABLE',option)
        end select
      end do
      call InputPopBlock(input,option)
      !PO: consider to include a force unit check in lookup table
      if ( .not. rtempvz_z_units_found ) then
        option%io_buffer = 'TOWG condition - RTEMPVZ/RTEMPVD: &
          &z/depth units must be entered for temperature table'
        call PrintErrMsg(option)
      else if ( rtempvz_z_units_found .and. rtempvz_found ) then
        call lkp_table%SetupVarUserUnits(ONE_INTEGER,usr_lenght_units,option)
      end if
      if ( .not. rtempvz_temp_units_found ) then
        option%io_buffer = 'TOWG condition - RTEMPVZ/RTEMPVD: &
          &temperature units must be entered for temperature table'
        call PrintErrMsg(option)
      else if ( rtempvz_temp_units_found .and. rtempvz_found ) then
        call lkp_table%SetupVarUserUnits(TWO_INTEGER,usr_temp_units,option)
      end if
      ! LookupTable unit conversion after reading units to make the table input
      ! independent from order the table instructions are given
      if ( rtempvz_found ) then
        call lkp_table%LookupTableVarConvFactors(option)
        call lkp_table%VarPointAndUnitConv(option)
        call lkp_table%SetupConstValExtrap(option)
        call lkp_table%LookupTableVarInitGradients(option)
        ! define table axis
        !PO: to move into table%SetUpIndependentVars(var1,var2), var 2 is optional
        allocate(lkp_table%axis1)
        allocate(lkp_table%axis1%values(size(lkp_table%var_data(1,:))))
        lkp_table%axis1%values = lkp_table%var_data(1,:)
        lkp_table%dims(1) = size(lkp_table%axis1%values(:))
        if ( .not. lkp_table%LookupTableVarIsSMInc(ONE_INTEGER) ) then
          option%io_buffer = 'TEMPERATURE_TABLE: temperature values must &
                              &be entered for strictly growing depths '
          call PrintErrMsg(option)
        end if
        nullify(lkp_table)
      else
        option%io_buffer = 'Flow condition - TEMPERATURE_TABLE: &
                            &RTEMPVZ or RTEMPVD tables not found'
        call PrintErrMsg(option)
      end if
   case default
     ! do nothing - do not throw error as other cards might be found
     ! in the mode-specific reading routine
 end select

end subroutine FlowConditionCommonRead

! ************************************************************************** !

subroutine TranConditionRead(condition,constraint_list, &
                             reaction_base,input,option)
  !
  ! Reads a transport condition from the input file
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Units_module

  use Transport_Constraint_Base_module
  use Transport_Constraint_RT_module
  use Transport_Constraint_NWT_module
  use Transport_Constraint_module
  use Reaction_Base_module
  use Reaction_Aux_module
  use NW_Transport_Aux_module

  implicit none

  type(tran_condition_type) :: condition
  type(tran_constraint_list_type) :: constraint_list
  class(reaction_base_type), pointer :: reaction_base
  type(input_type), pointer :: input
  type(option_type) :: option

  class(tran_constraint_base_type), pointer :: constraint
  class(tran_constraint_coupler_base_type), pointer :: constraint_coupler 
  class(tran_constraint_coupler_base_type), pointer :: cur_constraint_coupler
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word, internal_units
  character(len=MAXWORDLENGTH) :: default_time_units
  class(reaction_rt_type), pointer :: reaction
  class(reaction_nw_type), pointer :: reaction_nw
  PetscInt :: default_itype
  PetscBool :: found
  PetscInt :: icomp
  PetscBool :: minerals_exist
  PetscErrorCode :: ierr
  PetscReal :: conversion

  call PetscLogEventBegin(logging%event_tran_condition_read, &
                          ierr);CHKERRQ(ierr)

  select type(r=>reaction_base)
    class is(reaction_rt_type)
      reaction => r
    class is(reaction_nw_type)
      reaction_nw => r
  end select

  default_time_units = ''

  ! read the condition
  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CONDITION')

    select case(trim(word))

      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'TYPE','CONDITION')
        call StringToUpper(word)
        select case(word)
            case('DIRICHLET')
              condition%itype = DIRICHLET_BC
            case('DIRICHLET_ZERO_GRADIENT')
              condition%itype = DIRICHLET_ZERO_GRADIENT_BC
            case('EQUILIBRIUM')
              condition%itype = EQUILIBRIUM_SS
            case('NEUMANN')
              condition%itype = NEUMANN_BC
            case('MOLE','MOLE_RATE')
              condition%itype = MASS_RATE_SS
            case('ZERO_GRADIENT')
              condition%itype = ZERO_GRADIENT_BC
            case default
              call InputKeywordUnrecognized(input,word, &
                                            'transport condition type', &
                                            option)
        end select
      case('TIME_UNITS')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'UNITS','CONDITION')
        select case(trim(word))
          case('s','sec','min','m','hr','h','d','day','y','yr')
            default_time_units = trim(word)
          case default
            option%io_buffer = 'Units "' // trim(word) // '" not recognized.'
            call PrintErrMsg(option)
        end select
      case('CONSTRAINT_LIST')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONSTRAINT')

          if (InputCheckExit(input,option)) exit

          select case(option%itranmode)
            case(RT_MODE)
              constraint_coupler => TranConstraintCouplerRTCreate(option)
            case(NWT_MODE)
              constraint_coupler => TranConstraintCouplerNWTCreate(option)
          end select
          call InputReadDouble(input,option,constraint_coupler%time)
          call InputErrorMsg(input,option,'time','CONSTRAINT_LIST')
          ! time units are optional
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'constraint name','CONSTRAINT_LIST')
          ! read constraint name
          call InputReadWord(input,option, &
                             constraint_coupler%constraint_name, &
                             PETSC_TRUE)
          if (InputError(input)) then
            constraint_coupler%time_units = default_time_units
            constraint_coupler%constraint_name = trim(word)
          else
            constraint_coupler%time_units = word
          endif
          ! convert time units
          if (len_trim(constraint_coupler%time_units) > 0) then
            internal_units = 'sec'
            constraint_coupler%time = constraint_coupler%time* &
              UnitsConvertToInternal(constraint_coupler%time_units, &
                                     internal_units,option)
          endif
          ! add to end of list
          if (.not.associated(condition%constraint_coupler_list)) then
            condition%constraint_coupler_list => constraint_coupler
          else
            cur_constraint_coupler => condition%constraint_coupler_list
            do
              if (.not.associated(cur_constraint_coupler%next)) exit
              cur_constraint_coupler => cur_constraint_coupler%next
            enddo
            cur_constraint_coupler%next => constraint_coupler
          endif
        enddo
        call InputPopBlock(input,option)
        
      case('CONSTRAINT')
        select case(option%itranmode)
          case(RT_MODE)
            constraint_coupler => TranConstraintCouplerRTCreate(option)
            constraint => TranConstraintRTCreate(option)
          case(NWT_MODE)
            constraint_coupler => TranConstraintCouplerNWTCreate(option)
            constraint => TranConstraintNWTCreate(option)
        end select
        call InputReadWord(input,option,constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'constraint','name')
        option%io_buffer = 'Constraint: ' // trim(constraint%name)
        call PrintMsg(option)
        ! have to leave this select type as the argument lists differ
        select type(c=>constraint)
          class is(tran_constraint_rt_type)
            call TranConstraintRTRead(c,reaction,input,option)
          class is(tran_constraint_nwt_type)
            call TranConstraintNWTRead(c,reaction_nw,input,option)
        end select
        call TranConstraintAddToList(constraint,constraint_list)
        constraint_coupler%constraint => constraint
        constraint_coupler%time = 0.d0
        ! add to end of coupler list
        if (.not.associated(condition%constraint_coupler_list)) then
          condition%constraint_coupler_list => constraint_coupler
        else
          cur_constraint_coupler => condition%constraint_coupler_list
          do
            if (.not.associated(cur_constraint_coupler%next)) exit
            cur_constraint_coupler => cur_constraint_coupler%next
          enddo
          cur_constraint_coupler%next => constraint_coupler
        endif
      case default
        call InputKeywordUnrecognized(input,word,'transport condition',option)
    end select

  enddo
  call InputPopBlock(input,option)

  if (.not.associated(condition%constraint_coupler_list)) then
    option%io_buffer = 'No CONSTRAINT or CONSTRAINT_LIST defined in &
                       &TRANSPORT_CONDITION "' // trim(condition%name) // '".'
    call PrintErrMsg(option)
  endif
  
  if (len_trim(default_time_units) > 0) then
    internal_units = 'sec'
    conversion = UnitsConvertToInternal(default_time_units,internal_units, &
                                        option)
    cur_constraint_coupler => condition%constraint_coupler_list
    do
      if (.not.associated(cur_constraint_coupler)) exit
      if (len_trim(cur_constraint_coupler%time_units) == 0) then
        cur_constraint_coupler%time = cur_constraint_coupler%time*conversion
      endif
      cur_constraint_coupler => cur_constraint_coupler%next
    enddo
  endif

  call PetscLogEventEnd(logging%event_tran_condition_read,ierr);CHKERRQ(ierr)

end subroutine TranConditionRead

! ************************************************************************** !

subroutine ConditionReadValues(input,option,keyword,dataset_base, &
                               data_external_units,data_internal_units)
  ! 
  ! Read the value(s) of a condition variable
  !
  ! Author: Glenn Hammond
  ! Date: 10/31/07
  !

  use Input_Aux_module
  use String_module
  use Option_module
  use Logging_module
  use HDF5_Aux_module
  use Units_module
  use Dataset_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use hdf5

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  class(dataset_base_type), pointer :: dataset_base
  character(len=*) :: data_external_units
  character(len=*) :: data_internal_units

  character(len=MAXSTRINGLENGTH), pointer :: internal_unit_strings(:)
  class(dataset_ascii_type), pointer :: dataset_ascii
  character(len=MAXSTRINGLENGTH) :: string2, filename, hdf5_path
  character(len=MAXWORDLENGTH) :: word, realization_word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: length, i, icount
  PetscInt :: icol
  PetscInt :: ndims
  PetscInt, pointer :: dims(:)
  PetscReal, pointer :: real_buffer(:)
  PetscErrorCode :: ierr

  integer(HID_T) :: file_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: hdf5_err

  call PetscLogEventBegin(logging%event_flow_condition_read_values, &
                          ierr);CHKERRQ(ierr)

  ! dataset_base, though of type dataset_base_type, should always be created
  ! as dataset_ascii_type.
  dataset_ascii => DatasetAsciiCast(dataset_base)
  if (.not.associated(dataset_ascii)) then
    ! The dataset was not of type dataset_asci and was likely set to a different
    ! type.  There is a bug in the input file.
    option%io_buffer = 'Dataset associated with ' // trim(keyword) // &
      ' in the input file is already associated with a different dataset &
      &type.  Check for duplicate definitions of ' // trim(keyword) // '.'
    call PrintErrMsg(option)
  endif

  filename = ''
  realization_word = ''
  hdf5_path = ''

  internal_unit_strings => StringSplit(data_internal_units,',')

  input%ierr = 0
  string2 = trim(input%buf)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'file or value','CONDITION')
  call StringToUpper(word)
  length = len_trim(word)
  if (StringStartsWithAlpha(word)) then
    call InputPushCard(input,word,option)
    if (length == FOUR_INTEGER .and. &
        StringCompare(word,'FILE',FOUR_INTEGER)) then 
      input%err_buf2 = trim(keyword) // ', FILE'
      input%err_buf = 'keyword'
      call InputReadFilename(input,option,string2)
      if (input%ierr == 0) then
        filename = string2
      else
        option%io_buffer = 'The ability to read realization dependent &
          &datasets outside the DATASET block is no longer supported'
        call PrintErrMsg(option)
      endif

      if (len_trim(filename) < 2) then
        option%io_buffer = 'No filename listed under Flow_Condition: ' // &
                           trim(keyword)
        call PrintErrMsg(option)
      endif

      if (index(filename,'.h5') > 0) then
        write(option%io_buffer,'("Reading of HDF5 datasets for flow ", &
                                 &"conditions not currently supported.")')
        call PrintErrMsg(option)
#if 0
        if (len_trim(hdf5_path) < 1) then
          option%io_buffer = 'No hdf5 path listed under Flow_Condition: ' // &
                             trim(keyword)
          call PrintErrMsg(option)
        endif

        option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
        call PrintMsg(option)
        call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
        call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
        call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
        call h5pclose_f(prop_id,hdf5_err)

        hdf5_path = trim(hdf5_path) // trim(realization_word)
        call HDF5ReadNDimRealArray(option,file_id,hdf5_path,ndims,dims, &
                                   real_buffer)
        option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
        call PrintMsg(option)
        call h5fclose_f(file_id,hdf5_err)

        ! dims(1) = size of array
        ! dims(2) = number of data point in time
        if (dims(1)-1 == flow_dataset%time_series%rank) then
          ! alright, the 2d data is layed out in C-style.  now place it in
          ! the appropriate arrays
          allocate(flow_dataset%time_series%times(dims(2)))
          flow_dataset%time_series%times = UNINITIALIZED_DOUBLE
          allocate(flow_dataset%time_series%values(flow_dataset% &
                                                     time_series%rank,dims(2))) 
          flow_dataset%time_series%values = UNINITIALIZED_DOUBLE
          icount = 1
          do i = 1, dims(2)
            flow_dataset%time_series%times(i) = real_buffer(icount)
            icount = icount + 1
            do icol = 1, flow_dataset%time_series%rank
              flow_dataset%time_series%values(icol,i) = real_buffer(icount)
              icount = icount + 1
            enddo
          enddo
        else
          option%io_buffer = 'HDF condition data set rank does not match &
            &rank of internal data set.  Email Glenn for additions'
          call PrintErrMsg(option)
        endif
        if (associated(dims)) deallocate(dims)
        nullify(dims)
        if (associated(real_buffer)) deallocate(real_buffer)
        nullify(real_buffer)
#endif
! if 0
      else
        i = index(filename,'.',PETSC_TRUE)
        if (i > 2) then
          filename = filename(1:i-1) // trim(realization_word) // filename(i:)
        else
          filename = trim(filename) // trim(realization_word)
        endif
        error_string = 'CONDITION,' // trim(keyword) // ',FILE'
        call DatasetAsciiReadFile(dataset_ascii,filename,data_external_units, &
                                  data_internal_units,error_string,option)
        dataset_ascii%filename = filename
      endif
    else if (StringCompare(word,'DATASET')) then
      call InputReadWord(input,option,word,PETSC_TRUE)
      input%err_buf2 = trim(keyword) // ', DATASET'
      input%err_buf = 'dataset name'
      call InputErrorMsg(input,option)
      call DatasetDestroy(dataset_base)
      dataset_base => DatasetBaseCreate()
      dataset_base%name = word
    else if (length==FOUR_INTEGER .and. StringCompare(word,'LIST',length)) then 
      error_string = 'CONDITION,' // trim(keyword) // ',LIST'
      call DatasetAsciiReadList(dataset_ascii,input,data_external_units, &
                                data_internal_units,error_string,option)
    else if (StringCompare(word,'DBASE_VALUE')) then
      input%buf = trim(string2)
      error_string = 'CONDITION,' // trim(keyword) // ',SINGLE'
      call DatasetAsciiReadSingle(dataset_ascii,input,data_external_units, &
                                  data_internal_units,error_string,option)
    else
      option%io_buffer = 'Keyword "' // trim(word) // &
        '" not recognized in when reading condition values for "' // &
        trim(keyword) // '".'
      call PrintErrMsg(option)
    endif
  else
    input%buf = trim(string2)
    error_string = 'CONDITION,' // trim(keyword) // ',SINGLE'
    call DatasetAsciiReadSingle(dataset_ascii,input,data_external_units, &
                                data_internal_units,error_string,option)
  endif

  deallocate(internal_unit_strings)
  nullify(internal_unit_strings)

  call PetscLogEventEnd(logging%event_flow_condition_read_values, &
                        ierr);CHKERRQ(ierr)

end subroutine ConditionReadValues

! ************************************************************************** !

subroutine FlowConditionPrint(condition,option)
  !
  ! Prints flow condition info
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  !

  use Option_module
  use Dataset_module

  implicit none

  type(flow_condition_type) :: condition
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i

99 format(/,80('-'))

  write(option%fid_out,'(/,2x,''Flow Condition: '',a)') trim(condition%name)

  if (condition%sync_time_with_update) then
    string = 'yes'
  else
    string = 'no'
  endif
  write(option%fid_out,'(4x,''Synchronize time with update: '', a)') trim(string)
  write(option%fid_out,'(4x,''Time units: '', a)') trim(condition%time_units)
  write(option%fid_out,'(4x,''Length units: '', a)') trim(condition%length_units)

100 format(6x,a)
  write(option%fid_out,100) 'Datum:'
  if (associated(condition%datum)) then
    call DatasetPrint(condition%datum,option)
  endif

  do i=1, condition%num_sub_conditions
    call FlowConditionPrintSubCondition(condition%sub_condition_ptr(i)%ptr, &
                                        option)
  enddo
  write(option%fid_out,99)

end subroutine FlowConditionPrint

! ************************************************************************** !

subroutine FlowConditionPrintSubCondition(subcondition,option)
  !
  ! Prints flow subcondition info
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  !

  use Option_module
  use Dataset_module

  implicit none

  type(flow_sub_condition_type) :: subcondition
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  write(option%fid_out,'(/,4x,''Sub Condition: '',a)') trim(subcondition%name)
  string = GetSubConditionName(subcondition%itype)

  105 format(6x,'Type: ',a)
  write(option%fid_out,105) trim(string)

  110 format(6x,a)

  write(option%fid_out,110) 'Gradient:'
  if (associated(subcondition%gradient)) then
    call DatasetPrint(subcondition%gradient,option)
  endif

  write(option%fid_out,110) 'Data:'
  if (associated(subcondition%dataset)) then
    call DatasetPrint(subcondition%dataset,option)
  endif

end subroutine FlowConditionPrintSubCondition

! ************************************************************************** !

function GetSubConditionName(subcon_itype)
  !
  ! SubConditionName: Return name of subcondition
  !
  ! Author: Gautam Bisht
  ! Date: 10/16/13
  !

  implicit none

  PetscInt :: subcon_itype

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: GetSubConditionName

  select case(subcon_itype)
    case(DIRICHLET_BC)
      string = 'dirichlet'
    case(NEUMANN_BC)
      string = 'neumann'
    case(DIRICHLET_ZERO_GRADIENT_BC)
      string = 'dirichlet-zero gradient'
    case(MASS_RATE_SS)
      string = 'mass_rate'
    case(WELL_SS)
      string = 'well'
    case(HYDROSTATIC_BC)
      string = 'hydrostatic'
    case(HYDROSTATIC_CONDUCTANCE_BC)
      string = 'conductance'
    case(ZERO_GRADIENT_BC)
      string = 'zero gradient'
    case(HYDROSTATIC_SEEPAGE_BC)
      string = 'seepage'
    case(DIRICHLET_SEEPAGE_BC)
      string = 'dirichlet seepage'
    case(DIRICHLET_CONDUCTANCE_BC)
      string = 'dirichlet conductance'
    case(VOLUMETRIC_RATE_SS)
      string = 'volumetric rate'
    case(EQUILIBRIUM_SS)
      string = 'equilibrium'
    case(UNIT_GRADIENT_BC)
      string = 'unit gradient'
    case(SCALED_MASS_RATE_SS)
      string = 'scaled mass rate'
    case(SCALED_VOLUMETRIC_RATE_SS)
      string = 'scaled volumetric rate'
    case(HET_VOL_RATE_SS)
      string = 'heterogeneous volumetric rate'
    case(HET_MASS_RATE_SS)
      string = 'heterogeneous mass rate'
    case(HET_DIRICHLET_BC)
      string = 'heterogeneous dirichlet'
    case(HET_HYDROSTATIC_SEEPAGE_BC)
      string = 'heterogeneous seepage'
    case(HET_HYDROSTATIC_CONDUCTANCE_BC)
      string = 'heterogeneous conductance'
    case(ENERGY_RATE_SS)
      string = 'energy rate'
    case(SCALED_ENERGY_RATE_SS)
      string = 'scaled energy rate'
    case(HET_ENERGY_RATE_SS)
      string = 'heterogeneous energy rate'
    case(HET_SURF_HYDROSTATIC_SEEPAGE_BC)
      string = 'heterogeneous surface seepage'
    case(SPILLOVER_BC)
      string = 'spillover'
    case(WELL_MASS_RATE_TARGET)
      string = 'well target mass rate'
    case(WELL_MASS_RATE_MAX)
      string = 'well maximum mass rate'
    case(WELL_MASS_RATE_MIN)
      string = 'well minimum mass rate'
    case(WELL_VOL_RATE_TARGET)
      string = 'well target volumetric rate'
    case(WELL_VOL_RATE_MAX)
      string = 'well maximum volumetric rate'
    case(WELL_VOL_RATE_MIN)
      string = 'well minimum volumetric rate'
    case(WELL_BHP)
      string = 'well bottom hole pressure'
    case(WELL_BHP_MIN)
      string = 'well minimum bottom hole pressure'
    case(WELL_BHP_MAX)
      string = 'well maximum bottom hole pressure'
    case(SURFACE_DIRICHLET)
      string = 'surface_dirichlet'
    case(SURFACE_ZERO_GRADHEIGHT)
      string = 'surface_zero_gradheight'
    case(SURFACE_SPILLOVER)
      string = 'surface_spillover'
    end select

  GetSubConditionName = trim(string)

end function GetSubConditionName

! ************************************************************************** !

subroutine FlowConditionUpdate(condition_list,option)
  ! 
  ! Updates a transient condition
  !
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  !

  use Option_module
  use Dataset_module

  implicit none

  type(condition_list_type) :: condition_list
  type(option_type) :: option

  type(flow_condition_type), pointer :: condition
  type(flow_sub_condition_type), pointer :: sub_condition
  PetscInt :: isub_condition

  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    call DatasetUpdate(condition%datum,option)
    do isub_condition = 1, condition%num_sub_conditions

      sub_condition => condition%sub_condition_ptr(isub_condition)%ptr

      if (associated(sub_condition)) then
        call DatasetUpdate(sub_condition%dataset,option)
        call DatasetUpdate(sub_condition%gradient,option)
      endif

    enddo

    condition => condition%next

  enddo

end subroutine FlowConditionUpdate

! ************************************************************************** !

subroutine TranConditionUpdate(condition_list,option)
  ! 
  ! Updates a transient transport condition
  !
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  !
  use Option_module

  implicit none

  type(tran_condition_list_type) :: condition_list
  type(option_type) :: option

  type(tran_condition_type), pointer :: condition

  condition => condition_list%first
  do
    if (.not.associated(condition)) exit

    do
      if (associated(condition%cur_constraint_coupler%next)) then
        if (option%time >= condition%cur_constraint_coupler%next%time) then
          condition%cur_constraint_coupler => &
            condition%cur_constraint_coupler%next
        else
          exit
        endif
      else
        exit
      endif
    enddo
    condition => condition%next

  enddo

end subroutine TranConditionUpdate

! ************************************************************************** !

subroutine FlowConditionInitList(list)
  !
  ! Initializes a condition list
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  implicit none

  type(condition_list_type) :: list

  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine FlowConditionInitList

! ************************************************************************** !

subroutine FlowConditionAddToList(new_condition,list)
  !
  ! Adds a new condition to a condition list
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  implicit none

  type(flow_condition_type), pointer :: new_condition
  type(condition_list_type) :: list

  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition

end subroutine FlowConditionAddToList

! ************************************************************************** !

function FlowConditionGetPtrFromList(condition_name,condition_list)
  !
  ! Returns a pointer to the condition matching &
  ! condition_name
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  use String_module

  implicit none

  type(flow_condition_type), pointer :: FlowConditionGetPtrFromList
  character(len=MAXWORDLENGTH) :: condition_name
  type(condition_list_type) :: condition_list

  PetscInt :: length
  type(flow_condition_type), pointer :: condition

  nullify(FlowConditionGetPtrFromList)
  condition => condition_list%first

  do
    if (.not.associated(condition)) exit
    length = len_trim(condition_name)
    if (length == len_trim(condition%name) .and. &
        StringCompare(condition%name,condition_name, &
                      length)) then
      FlowConditionGetPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo

end function FlowConditionGetPtrFromList

! ************************************************************************** !

subroutine TranConditionInitList(list)
  !
  ! Initializes a transport condition list
  !
  ! Author: Glenn Hammond
  ! Date: 10/13/08
  !

  implicit none

  type(tran_condition_list_type) :: list

  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine TranConditionInitList

! ************************************************************************** !

subroutine TranConditionAddToList(new_condition,list)
  !
  ! Adds a new condition to a transport condition list
  !
  ! Author: Glenn Hammond
  ! Date: 10/13/08
  !

  implicit none

  type(tran_condition_type), pointer :: new_condition
  type(tran_condition_list_type) :: list

  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition

end subroutine TranConditionAddToList

! ************************************************************************** !

function TranConditionGetPtrFromList(condition_name,condition_list)
  !
  ! Returns a pointer to the condition matching
  ! condition_name
  !
  ! Author: Glenn Hammond
  ! Date: 10/13/08
  !

  use String_module

  implicit none

  type(tran_condition_type), pointer :: TranConditionGetPtrFromList
  character(len=MAXWORDLENGTH) :: condition_name
  type(tran_condition_list_type) :: condition_list

  PetscInt :: length
  type(tran_condition_type), pointer :: condition

  nullify(TranConditionGetPtrFromList)
  condition => condition_list%first

  do
    if (.not.associated(condition)) exit
    length = len_trim(condition_name)
    if (length == len_trim(condition%name) .and. &
        StringCompare(condition%name,condition_name, &
                        length)) then
      TranConditionGetPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo

end function TranConditionGetPtrFromList

! ************************************************************************** !

function FlowConditionIsTransient(condition)
  !
  ! Returns PETSC_TRUE
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  !

  use Dataset_module

  implicit none

  type(flow_condition_type) :: condition

  PetscBool :: FlowConditionIsTransient

  FlowConditionIsTransient = PETSC_FALSE

  if (DatasetIsTransient(condition%datum) .or. &
      FlowSubConditionIsTransient(condition%pressure) .or. &
      FlowSubConditionIsTransient(condition%temperature) .or. &
      FlowSubConditionIsTransient(condition%concentration) .or. &
      FlowSubConditionIsTransient(condition%saturation) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%well) .or. &
      FlowSubConditionIsTransient(condition%enthalpy) .or. &
      FlowSubConditionIsTransient(condition%energy_rate) .or. &
      FlowSubConditionIsTransient(condition%energy_flux) .or. &
      FlowConditionTOilImsIsTransient(condition%toil_ims) .or. &
      FlowConditionTOWGIsTransient(condition%towg) .or. &
      FlowConditionGeneralIsTransient(condition%general) .or. &
      FlowConditionHydrateIsTransient(condition%hydrate)) then
    FlowConditionIsTransient = PETSC_TRUE
  endif

end function FlowConditionIsTransient

! ************************************************************************** !

function FlowConditionGeneralIsTransient(condition)
  !
  ! Returns PETSC_TRUE
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  !

  use Dataset_module

  implicit none

  type(flow_general_condition_type), pointer :: condition

  PetscBool :: FlowConditionGeneralIsTransient

  FlowConditionGeneralIsTransient = PETSC_FALSE

  if (.not.associated(condition)) return

  if (FlowSubConditionIsTransient(condition%liquid_pressure) .or. &
      FlowSubConditionIsTransient(condition%gas_pressure) .or. &
      FlowSubConditionIsTransient(condition%gas_saturation) .or. &
      FlowSubConditionIsTransient(condition%relative_humidity) .or. &
      FlowSubConditionIsTransient(condition%mole_fraction) .or. &
      FlowSubConditionIsTransient(condition%temperature) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%liquid_flux) .or. &
      FlowSubConditionIsTransient(condition%gas_flux) .or. &
      FlowSubConditionIsTransient(condition%energy_flux)) then
    FlowConditionGeneralIsTransient = PETSC_TRUE
  endif

end function FlowConditionGeneralIsTransient

! ************************************************************************** !

function FlowConditionHydrateIsTransient(condition)
  !
  ! Returns PETSC_TRUE
  !
  ! Author: Michael Nole
  ! Date: 11/04/20
  !

  use Dataset_module

  implicit none

  type(flow_hydrate_condition_type), pointer :: condition

  PetscBool :: FlowConditionHydrateIsTransient

  FlowConditionHydrateIsTransient = PETSC_FALSE

  if (.not.associated(condition)) return

  if (FlowSubConditionIsTransient(condition%liquid_pressure) .or. &
      FlowSubConditionIsTransient(condition%gas_pressure) .or. &
      FlowSubConditionIsTransient(condition%gas_saturation) .or. &
      FlowSubConditionIsTransient(condition%hydrate_saturation) .or. &
      FlowSubConditionIsTransient(condition%liquid_saturation) .or. &
      FlowSubConditionIsTransient(condition%ice_saturation) .or. &
      FlowSubConditionIsTransient(condition%relative_humidity) .or. &
      FlowSubConditionIsTransient(condition%mole_fraction) .or. &
      FlowSubConditionIsTransient(condition%temperature) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%liquid_flux) .or. &
      FlowSubConditionIsTransient(condition%gas_flux) .or. &
      FlowSubConditionIsTransient(condition%energy_flux)) then
    FlowConditionHydrateIsTransient = PETSC_TRUE
  endif

end function FlowConditionHydrateIsTransient

! ************************************************************************** !

function FlowConditionTOilImsIsTransient(condition)
  !
  ! Returns PETSC_TRUE
  !
  ! Author: Paolo Orsini
  ! Date: 10/20/15
  !

  use Dataset_module

  implicit none

  type(flow_toil_ims_condition_type), pointer :: condition

  PetscBool :: FlowConditionTOilImsIsTransient

  FlowConditionTOilImsIsTransient = PETSC_FALSE

  if (.not.associated(condition)) return

  if (FlowSubConditionIsTransient(condition%pressure) .or. &
      FlowSubConditionIsTransient(condition%saturation) .or. &
      FlowSubConditionIsTransient(condition%temperature) .or. &
      FlowSubConditionIsTransient(condition%enthalpy) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%liquid_flux) .or. &
      FlowSubConditionIsTransient(condition%oil_flux) .or. &
      FlowSubConditionIsTransient(condition%energy_flux) .or. &
      FlowSubConditionIsTransient(condition%owc) .or. &
      FlowSubConditionIsTransient(condition%liq_press_grad)) then
    FlowConditionTOilImsIsTransient = PETSC_TRUE
  endif

end function FlowConditionTOilImsIsTransient

! ************************************************************************** !

function FlowConditionTOWGIsTransient(condition)
  !
  ! Returns PETSC_TRUE - if it founds a trannsient sub condition
  !
  ! Author: Paolo Orsini
  ! Date: 11/09/17
  !

  use Dataset_module

  implicit none

  type(flow_towg_condition_type), pointer :: condition

  PetscBool :: FlowConditionTOWGIsTransient

  FlowConditionTOWGIsTransient = PETSC_FALSE

  if (.not.associated(condition)) return

  if (FlowSubConditionIsTransient(condition%oil_pressure) .or. &
      FlowSubConditionIsTransient(condition%gas_pressure) .or. &
      FlowSubConditionIsTransient(condition%oil_saturation) .or. &
      FlowSubConditionIsTransient(condition%gas_saturation) .or. &
      FlowSubConditionIsTransient(condition%solvent_saturation) .or. &
      FlowSubConditionIsTransient(condition%bubble_point) .or. &
      FlowSubConditionIsTransient(condition%liquid_flux) .or. &
      FlowSubConditionIsTransient(condition%oil_flux) .or. &
      FlowSubConditionIsTransient(condition%gas_flux) .or. &
      FlowSubConditionIsTransient(condition%solvent_flux) .or. &
      FlowSubConditionIsTransient(condition%energy_flux) .or. &
      FlowSubConditionIsTransient(condition%temperature) .or. &
      FlowSubConditionIsTransient(condition%enthalpy) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%bhp_pressure) ) then
    FlowConditionTOWGIsTransient = PETSC_TRUE
  endif

end function FlowConditionTOWGIsTransient

! ************************************************************************** !

function FlowSubConditionIsTransient(sub_condition)
  !
  ! Returns PETSC_TRUE
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  !

  use Dataset_module

  implicit none

  type(flow_sub_condition_type), pointer :: sub_condition

  PetscBool :: FlowSubConditionIsTransient

  FlowSubConditionIsTransient = PETSC_FALSE

  if (associated(sub_condition)) then
    if (DatasetIsTransient(sub_condition%dataset) .or. &
        DatasetIsTransient(sub_condition%gradient)) then
      FlowSubConditionIsTransient = PETSC_TRUE
    endif
  endif

end function FlowSubConditionIsTransient

! ************************************************************************** !

function FlowConditionUnknownItype(condition,message,type_name)
  !
  ! Returns a string indicating which flow condition has a wrong type.
  !
  ! Author: Glenn Hammond
  ! Date: 03/21/16
  !
  implicit none

  type(flow_condition_type) :: condition
  character(len=*) :: message
  character(len=*) :: type_name

  character(len=MAXSTRINGLENGTH) :: FlowConditionUnknownItype

  FlowConditionUnknownItype = 'Unknown TYPE (' // trim(type_name) // &
    ') for ' // trim(message) // ' within FLOW_CONDITION "' // &
    trim(condition%name) // '".'

end function FlowConditionUnknownItype

! **************************************************************************** !

subroutine FlowCondInputRecord(flow_condition_list,option)
  !
  ! Prints ingested flow condition information to
  ! the input record file.
  !
  ! Author: Jenn Frederick
  ! Date: 04/19/2016
  !
  use Option_module
  use Dataset_Base_class

  implicit none

  type(condition_list_type), pointer :: flow_condition_list
  type(option_type), pointer :: option

  type(flow_condition_type), pointer :: cur_fc
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'FLOW CONDITIONS'

  cur_fc => flow_condition_list%first
  do
    if (.not.associated(cur_fc)) exit
    write(id,'(a29)',advance='no') 'flow condition name: '
    write(id,'(a)') adjustl(trim(cur_fc%name))
    if (cur_fc%num_sub_conditions > 0) then
      do k = 1,cur_fc%num_sub_conditions
        write(id,'(a29)',advance='no') 'sub condition name: '
        write(id,'(a)') adjustl(trim(cur_fc%sub_condition_ptr(k)%ptr%name))
        write(id,'(a29)',advance='no') 'sub condition type: '
        write(id,'(a)') adjustl(trim(cur_fc%sub_condition_ptr(k)%ptr%ctype))
        if (associated(cur_fc%sub_condition_ptr(k)%ptr%dataset)) then
          call DatasetBasePrint(cur_fc%sub_condition_ptr(k)%ptr%dataset,option)
          ! DatasetBasePrint doesn't seem to do anything?
        endif
        if (associated(cur_fc%sub_condition_ptr(k)%ptr%gradient)) then
          call DatasetBasePrint(cur_fc%sub_condition_ptr(k)%ptr%gradient,option)
        endif
      enddo
    endif

    write(id,'(a29)') '---------------------------: '
    cur_fc => cur_fc%next
  enddo

end subroutine FlowCondInputRecord

! **************************************************************************** !

subroutine TranCondInputRecord(tran_condition_list,option)
  !
  ! Prints ingested transport condition information to
  ! the input record file.
  !
  ! Author: Jenn Frederick
  ! Date: 04/19/2016
  !
  use Option_module
  use Dataset_Base_class
  Use Transport_Constraint_module
  Use Transport_Constraint_RT_module

  implicit none

  type(tran_condition_list_type), pointer :: tran_condition_list
  type(option_type), pointer :: option

  class(tran_constraint_coupler_base_type), pointer :: cur_constraint_coupler
  class(tran_constraint_base_type), pointer :: constraint
  type(tran_condition_type), pointer :: cur_condition
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'TRANSPORT CONDITIONS'

  cur_condition => tran_condition_list%first
  do
    if (.not.associated(cur_condition)) exit
    write(id,'(a29)',advance='no') 'transport condition name: '
    write(id,'(a)') adjustl(trim(cur_condition%name))
    write(id,'(a29)',advance='no') 'transport condition type: '
    select case (cur_condition%itype)
      case(DIRICHLET_BC)
        write(id,'(a)') 'dirichlet'
      case(DIRICHLET_ZERO_GRADIENT_BC)
        write(id,'(a)') 'dirichlet_zero_gradient'
      case(EQUILIBRIUM_SS)
        write(id,'(a)') 'equilibrium'
      case(NEUMANN_BC)
        write(id,'(a)') 'neumann'
      case(MASS_RATE_SS)
        write(id,'(a)') 'mole_rate'
      case(ZERO_GRADIENT_BC)
        write(id,'(a)') 'zero_gradient'
    end select
    write(id,'(a29)',advance='no') 'is transient?: '
    if (cur_condition%is_transient) then
      write(id,'(a)') 'YES'
    else
      write(id,'(a)') 'NO'
    endif
    cur_constraint_coupler => cur_condition%constraint_coupler_list
    do
      if (.not.associated(cur_constraint_coupler)) exit
      write(id,'(a29)',advance='no') 'transport constraint name: '
      write(id,'(a)') adjustl(trim(cur_constraint_coupler%constraint_name))

      constraint => cur_constraint_coupler%constraint

      select type(c=>constraint)
        class is(tran_constraint_rt_type)
          ! aqueous species concentraion constraint
          if (associated(c%aqueous_species)) then
            do k = 1,size(c%aqueous_species%names)
              write(id,'(a29)',advance='no') 'aqueous species constraint: '
              write(string,*) trim(c%aqueous_species%names(k))
              select case (c%aqueous_species%constraint_type(k))
                case(CONSTRAINT_FREE)
                  string = trim(string) // ', free'
                case(CONSTRAINT_TOTAL)
                  string = trim(string) // ', total'
                case(CONSTRAINT_TOTAL_SORB)
                  string = trim(string) // ', total_sorb'
                case(CONSTRAINT_PH)
                  string = trim(string) // ', ph'
                case(CONSTRAINT_LOG)
                  string = trim(string) // ', log'
                case(CONSTRAINT_MINERAL)
                  string = trim(string) // ', mineral'
                case(CONSTRAINT_GAS)
                  string = trim(string) // ', gas'
                case(CONSTRAINT_SUPERCRIT_CO2)
                  string = trim(string) // ', super critical CO2'
                case(CONSTRAINT_CHARGE_BAL)
                  string = trim(string) // ', charge balance'
              end select
              write(word1,*) c%aqueous_species%constraint_conc(k)
              write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                              // ' mol'
            enddo
          endif

          ! free-ion guess constraint
          if (associated(c%free_ion_guess)) then
            do k = 1,size(c%free_ion_guess%names)
              write(id,'(a29)',advance='no') 'free ion guess constraint: '
              write(string,*) trim(c%free_ion_guess%names(k))
              write(word1,*) c%free_ion_guess%conc(k)
              write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                              // ' mol'
            enddo
          endif

          ! mineral constraint
          if (associated(c%minerals)) then
            do k = 1,size(c%minerals%names)
              write(id,'(a29)',advance='no') 'mineral vol. frac. constraint: '
              write(string,*) trim(c%minerals%names(k))
              write(word1,*) c%minerals%constraint_vol_frac(k)
              write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                              // ' m^3/m^3'
              write(id,'(a29)',advance='no') 'mineral area constraint: '
              write(string,*) trim(c%minerals%names(k))
              write(word1,*) c%minerals%constraint_area(k)
              write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                              // ' m^2/m^3'
            enddo
          endif

          ! surface complexes constraint
          if (associated(c%surface_complexes)) then
            do k = 1,size(c%surface_complexes%names)
              write(id,'(a29)',advance='no') 'surface complex constraint: '
              write(string,*) trim(c%surface_complexes%names(k))
              write(word1,*) c%surface_complexes% &
                               constraint_conc(k)
              write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                              // ' mol/m^3'
            enddo
          endif

          ! colloids constraint
          if (associated(c%colloids)) then
            do k = 1,size(c%colloids%names)
              write(id,'(a29)',advance='no') 'colloid mobile constraint: '
              write(string,*) trim(c%colloids%names(k))
              write(word1,*) c%colloids%constraint_conc_mob(k)
              write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                              // ' mol'
              write(id,'(a29)',advance='no') 'colloid immobile constraint: '
              write(string,*) trim(c%colloids%names(k))
              write(word1,*) c%colloids%constraint_conc_imb(k)
              write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                              // ' mol'
            enddo
          endif

          ! immobile species constraint
          if (associated(c%immobile_species)) then
            do k = 1,size(c%immobile_species%names)
              write(id,'(a29)',advance='no') 'immobile species constraint: '
              write(string,*) trim(c%immobile_species%names(k))
              write(word1,*) c%immobile_species% &
                               constraint_conc(k)
              write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                              // ' mol/m^3'
            enddo
          endif
      end select

      cur_constraint_coupler => cur_constraint_coupler%next
    enddo

    write(id,'(a29)') '---------------------------: '
    cur_condition => cur_condition%next
  enddo

end subroutine TranCondInputRecord

! ************************************************************************** !

subroutine FlowConditionDestroyList(condition_list)
  !
  ! Deallocates a list of conditions
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  implicit none

  type(condition_list_type), pointer :: condition_list

  type(flow_condition_type), pointer :: condition, prev_condition

  if (.not.associated(condition_list)) return

  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call FlowConditionDestroy(prev_condition)
  enddo

  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  if (associated(condition_list%array)) deallocate(condition_list%array)
  nullify(condition_list%array)

  deallocate(condition_list)
  nullify(condition_list)

end subroutine FlowConditionDestroyList

! ************************************************************************** !

subroutine FlowConditionDestroy(condition)
  !
  ! Deallocates a condition
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !

  use Dataset_module
  use Dataset_Ascii_class
  use Utility_module

  implicit none

  type(flow_condition_type), pointer :: condition

  class(dataset_ascii_type), pointer :: dataset_ascii
  PetscInt :: i

  if (.not.associated(condition)) return

  ! if dataset_ascii_type, destroy.  Otherwise, they are in another list
  dataset_ascii => DatasetAsciiCast(condition%datum)
  ! dataset_ascii will be NULL if not dataset_ascii_type
  call DatasetAsciiDestroy(dataset_ascii)

  if (associated(condition%sub_condition_ptr)) then
    ! only nullify the pointers; don't destroy as they are destroyed below
    do i=1,condition%num_sub_conditions
      nullify(condition%sub_condition_ptr(i)%ptr)
    enddo
    deallocate(condition%sub_condition_ptr)
    nullify(condition%sub_condition_ptr)
  endif

  call DeallocateArray(condition%itype)

  call FlowSubConditionDestroy(condition%datum_z)
  call FlowSubConditionDestroy(condition%pressure)
  call FlowSubConditionDestroy(condition%saturation)
  call FlowSubConditionDestroy(condition%rate)
  call FlowSubConditionDestroy(condition%well)
  call FlowSubConditionDestroy(condition%temperature)
  call FlowSubConditionDestroy(condition%concentration)
  call FlowSubConditionDestroy(condition%enthalpy)
  call FlowSubConditionDestroy(condition%energy_rate)

  call TimeStorageDestroy(condition%default_time_storage)
  call FlowGeneralConditionDestroy(condition%general)
  call FlowTOilConditionDestroy(condition%toil_ims)
  call FlowTOWGConditionDestroy(condition%towg)
  call LookupTableDestroy(condition%rtempvz_table)

  nullify(condition%next)

  deallocate(condition)
  nullify(condition)

end subroutine FlowConditionDestroy

! ************************************************************************** !

subroutine FlowGeneralConditionDestroy(general_condition)
  !
  ! Destroys a general mode condition
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/11
  !

  use Option_module

  implicit none

  type(flow_general_condition_type), pointer :: general_condition

  if (.not.associated(general_condition)) return

  call FlowSubConditionDestroy(general_condition%liquid_pressure)
  call FlowSubConditionDestroy(general_condition%gas_pressure)
  call FlowSubConditionDestroy(general_condition%gas_saturation)
  call FlowSubConditionDestroy(general_condition%relative_humidity)
  call FlowSubConditionDestroy(general_condition%mole_fraction)
  call FlowSubConditionDestroy(general_condition%temperature)
  call FlowSubConditionDestroy(general_condition%liquid_flux)
  call FlowSubConditionDestroy(general_condition%gas_flux)
  call FlowSubConditionDestroy(general_condition%energy_flux)
  call FlowSubConditionDestroy(general_condition%rate)

  deallocate(general_condition)
  nullify(general_condition)

end subroutine FlowGeneralConditionDestroy

! ************************************************************************** !

subroutine FlowTOilConditionDestroy(toil_ims_condition)
  !
  ! Destroys a toil_ims mode condition
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/06/15
  !

  use Option_module

  implicit none

  type(flow_toil_ims_condition_type), pointer :: toil_ims_condition

  if (.not.associated(toil_ims_condition)) return

  call FlowSubConditionDestroy(toil_ims_condition%pressure)
  call FlowSubConditionDestroy(toil_ims_condition%saturation)
  call FlowSubConditionDestroy(toil_ims_condition%temperature)
  call FlowSubConditionDestroy(toil_ims_condition%enthalpy)
  call FlowSubConditionDestroy(toil_ims_condition%liquid_flux)
  call FlowSubConditionDestroy(toil_ims_condition%oil_flux)
  call FlowSubConditionDestroy(toil_ims_condition%energy_flux)
  call FlowSubConditionDestroy(toil_ims_condition%rate)
  call FlowSubConditionDestroy(toil_ims_condition%owc)
  call FlowSubConditionDestroy(toil_ims_condition%owc_z)
  call FlowSubConditionDestroy(toil_ims_condition%pcow_owc)
  call FlowSubConditionDestroy(toil_ims_condition%liq_press_grad)

  deallocate(toil_ims_condition)
  nullify(toil_ims_condition)

end subroutine FlowTOilConditionDestroy

! ************************************************************************** !

subroutine FlowTOWGConditionDestroy(towg_condition)
  !
  ! Destroys a towg mode condition
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/06/15
  !

  use Option_module

  implicit none

  type(flow_towg_condition_type), pointer :: towg_condition

  if (.not.associated(towg_condition)) return

  call FlowSubConditionDestroy(towg_condition%oil_pressure)
  call FlowSubConditionDestroy(towg_condition%gas_pressure)
  call FlowSubConditionDestroy(towg_condition%oil_saturation)
  call FlowSubConditionDestroy(towg_condition%gas_saturation)
  call FlowSubConditionDestroy(towg_condition%bubble_point)
  call FlowSubConditionDestroy(towg_condition%solvent_saturation)
  call FlowSubConditionDestroy(towg_condition%liquid_flux)
  call FlowSubConditionDestroy(towg_condition%oil_flux)
  call FlowSubConditionDestroy(towg_condition%gas_flux)
  call FlowSubConditionDestroy(towg_condition%energy_flux)
  call FlowSubConditionDestroy(towg_condition%temperature)
  call FlowSubConditionDestroy(towg_condition%enthalpy)
  call FlowSubConditionDestroy(towg_condition%rate)
  call FlowSubConditionDestroy(towg_condition%owc_z)
  call FlowSubConditionDestroy(towg_condition%pcow_owc)
  call FlowSubConditionDestroy(towg_condition%ogc_z)
  call FlowSubConditionDestroy(towg_condition%pcog_ogc)
  call LookupTableDestroy(towg_condition%pbvz_table)

  deallocate(towg_condition)
  nullify(towg_condition)

end subroutine FlowTOWGConditionDestroy

! ************************************************************************** !

subroutine FlowSubConditionDestroy(sub_condition)
  !
  ! Destroys a sub_condition
  !
  ! Author: Glenn Hammond
  ! Date: 02/04/08
  !

  use Dataset_module
  use Dataset_Ascii_class

  implicit none

  type(flow_sub_condition_type), pointer :: sub_condition

  class(dataset_ascii_type), pointer :: dataset_ascii

  if (.not.associated(sub_condition)) return

  ! if dataset_ascii_type, destroy.  Otherwise, they are in another list
  dataset_ascii => DatasetAsciiCast(sub_condition%dataset)
  ! dataset_ascii will be NULL if not dataset_ascii_type
  call DatasetAsciiDestroy(dataset_ascii)
  dataset_ascii => DatasetAsciiCast(sub_condition%gradient)
  call DatasetAsciiDestroy(dataset_ascii)

  deallocate(sub_condition)
  nullify(sub_condition)

end subroutine FlowSubConditionDestroy

! ************************************************************************** !

subroutine TranConditionDestroyList(condition_list)
  !
  ! Deallocates a list of conditions
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  implicit none

  type(tran_condition_list_type), pointer :: condition_list

  type(tran_condition_type), pointer :: condition, prev_condition

  if (.not.associated(condition_list)) return

  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call TranConditionDestroy(prev_condition)
  enddo

  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  if (associated(condition_list%array)) deallocate(condition_list%array)
  nullify(condition_list%array)

  deallocate(condition_list)
  nullify(condition_list)

end subroutine TranConditionDestroyList

! ************************************************************************** !

subroutine TranConditionDestroy(condition)
  !
  ! Deallocates a condition
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !
  use Transport_Constraint_module

  implicit none

  type(tran_condition_type), pointer :: condition

  if (.not.associated(condition)) return

  if (associated(condition%constraint_coupler_list)) &
    call TranConstraintCouplerDestroy(condition%constraint_coupler_list)

  deallocate(condition)
  nullify(condition)

end subroutine TranConditionDestroy

end module Condition_module
