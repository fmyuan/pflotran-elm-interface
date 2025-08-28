module PM_Ponded_Water_class

#include "petsc/finclude/petscvec.h"
  use petscsys
  use petscvec
  use PM_Base_class
  use Realization_Subsurface_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter :: METRIC_ELEVATION = 1
  PetscInt, parameter :: METRIC_AREA = 2
  PetscInt, parameter :: METRIC_UNDERLYING_VOLUME = 3

  type, public, extends(pm_base_type) :: pm_ponded_water_type
    class(realization_subsurface_type), pointer :: realization
    type(ponded_water_aux_type), pointer :: aux
    procedure(PMPondedWaterUpdateDummy), pointer :: Update => null()
  contains
    procedure, public :: Setup => PMPondedWaterSetup
    procedure, public :: SetRealization => PMPondedWaterSetRealization
    procedure, public :: InitializeRun => PMPondedWaterInitializeRun
    procedure, public :: FinalizeRun => PMPondedWaterFinalizeRun
    procedure, public :: Destroy => PMPondedWaterDestroy
  end type pm_ponded_water_type

  type :: ponded_water_aux_type
    character(len=MAXWORDLENGTH) :: region_name
    PetscReal :: pond_elevation
    PetscReal, pointer :: metric_by_elevation(:,:)
    type(ponded_water_aux_type), pointer :: next
  end type ponded_water_aux_type

  ! interface blocks
  interface
    subroutine PMPondedWaterUpdateDummy(this,time,ierr)
      use petscsys
      import :: pm_ponded_water_type
      implicit none
      class(pm_ponded_water_type) :: this
      PetscReal :: time
      PetscErrorCode :: ierr
    end subroutine PMPondedWaterUpdateDummy
  end interface

  public :: PMPondedWaterCreate, &
            PMPondedWaterInit, &
            PMPondedWaterCast, &
            PMPondedWaterRead

contains

! ************************************************************************** !

function PMPondedWaterCreate()
  !
  ! Creates a ponded water process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/25
  !

  class(pm_ponded_water_type), pointer :: PMPondedWaterCreate

  class(pm_ponded_water_type), pointer :: pm

  allocate(pm)
  call PMPondedWaterInit(pm)

  PMPondedWaterCreate => pm

end function PMPondedWaterCreate

! ************************************************************************** !

subroutine PMPondedWaterInit(this)
  !
  ! Initializes the ponded water process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/25

  class(pm_ponded_water_type) :: this

  call PMBaseInit(this)
  this%header = 'Ponded Water'
  nullify(this%aux)
  nullify(this%realization)
  this%Update => null()

end subroutine PMPondedWaterInit

! ************************************************************************** !

function PondedWaterAuxCreate()
  !
  ! Creates a ponded water auxiliary data object
  !
  ! Author: Glenn Hammond
  ! Date: 02/19/25
  !

  class(ponded_water_aux_type), pointer :: PondedWaterAuxCreate

  class(ponded_water_aux_type), pointer :: aux

  allocate(aux)
  aux%region_name = ''
  aux%pond_elevation = UNINITIALIZED_DOUBLE
  nullify(aux%metric_by_elevation)
  nullify(aux%next)

  PondedWaterAuxCreate => aux

end function PondedWaterAuxCreate

! ************************************************************************** !

subroutine PMPondedWaterSetup(this)
  !
  ! Sets up ponded water process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/25

  class(pm_ponded_water_type) :: this

end subroutine PMPondedWaterSetup

! ************************************************************************** !

function PMPondedWaterCast(this)
  !
  ! Casts a base process model to ponded water
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/25

  use Option_module

  class(pm_base_type), pointer :: this

  class(pm_ponded_water_type), pointer :: PMPondedWaterCast

  nullify(PMPondedWaterCast)
  if (.not.associated(this)) return
  select type (this)
    class is (pm_ponded_water_type)
      PMPondedWaterCast => this
    class default
      this%option%io_buffer = 'Cannot cast pm_base_type to &
        &pm_ponded_water_type in PMPondedWaterCast.'
      call PrintErrMsg(this%option)
  end select

end function PMPondedWaterCast

! ************************************************************************** !

subroutine PMPondedWaterRead(input,option,this)
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/25
  !
  use Input_Aux_module
  use Option_module
  use String_module

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_ponded_water_type), pointer :: this

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_str

  error_str = 'SIMULATION,PROCESS_MODELS,PONDED_WATER'

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_str)
    call StringToUpper(keyword)

    select case(trim(keyword))
      case default
        call InputKeywordUnrecognized(input,keyword,error_str,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine PMPondedWaterRead

! ************************************************************************** !

subroutine PMPondedWaterSetRealization(this,realization)
  !
  ! Sets the realization pointer
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/25

  use Realization_Subsurface_class

  class(pm_ponded_water_type) :: this
  class(realization_subsurface_type), pointer :: realization

  this%realization => realization

end subroutine PMPondedWaterSetRealization

! ************************************************************************** !

recursive subroutine PMPondedWaterInitializeRun(this)
  !
  ! Initializes the process model for the simulation
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/25

  use Option_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Dataset_Ascii_class
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_module
  use Grid_module
  use ZFlow_Aux_module, only : ZFLOW_COND_WATER_INDEX
  use Utility_module, only : DotProduct

  class(pm_ponded_water_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(ponded_water_aux_type), pointer :: cur_ponded_water_aux
  type(ponded_water_aux_type), pointer :: prev_ponded_water_aux
  PetscInt :: water_index
  PetscInt :: num_connections
  PetscInt :: num_ponded_water_bc
  PetscInt :: i, iface, iconn, iconn_prev_elev
  PetscReal :: sum_area
  PetscInt, allocatable :: cell_id(:)
  PetscInt, allocatable :: face_id(:)
  PetscInt, allocatable :: permutation(:)
  PetscReal, allocatable :: face_projected_area(:,:)
  PetscReal, allocatable :: face_elevation(:)
  PetscReal, pointer :: metric_by_elevation(:,:)
  PetscReal :: gravity_unit_vector(3)
  PetscReal :: cell_z, face_z, delta_z
  PetscErrorCode :: ierr

  option => this%realization%option
  patch => this%realization%patch
  grid => patch%grid

  if (this%option%iflowmode /= ZFLOW_MODE) then
  endif

  this%Update => PMPondedWaterUpdate1

  gravity_unit_vector = dabs(option%gravity/ &
                             sqrt(DotProduct(option%gravity,option%gravity)))

  ! count number of ponded water boundary conditions
  num_ponded_water_bc = 0
  boundary_condition => patch%boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    if (boundary_condition%flow_bc_type(1) == PONDED_WATER_BC) then
      num_ponded_water_bc = num_ponded_water_bc + 1
    endif
    boundary_condition => boundary_condition%next
  enddo

  ! count number of ponded water boundary conditions
  boundary_condition => patch%boundary_condition_list%first
  nullify(prev_ponded_water_aux)
  do
    if (.not.associated(boundary_condition)) exit
    if (boundary_condition%flow_bc_type(1) == PONDED_WATER_BC) then
      cur_ponded_water_aux => PondedWaterAuxCreate()
      cur_ponded_water_aux%region_name = boundary_condition%region_name
      cur_connection_set => boundary_condition%connection_set
      num_connections = cur_connection_set%num_connections
      allocate(permutation(num_connections))
      allocate(face_elevation(num_connections))
      allocate(face_projected_area(2,num_connections))
      allocate(cell_id(num_connections))
      allocate(face_id(num_connections))
      cell_id = UNINITIALIZED_INTEGER
      face_id = UNINITIALIZED_INTEGER
      face_elevation = UNINITIALIZED_DOUBLE
      face_projected_area = UNINITIALIZED_DOUBLE
      do iconn = 1, num_connections
        permutation(iconn) = iconn
        cell_id(iconn) = grid%nG2A(grid%nL2G(cur_connection_set%id_dn(iconn)))
        if (dabs(cur_connection_set%dist(1,iconn)) > 0.1) then
          if (cur_connection_set%dist(1,iconn) > 0.1) then
            face_id(iconn) = 1
          elseif (cur_connection_set%dist(1,iconn) < 0.1) then
            face_id(iconn) = 2
          endif
        elseif (dabs(cur_connection_set%dist(2,iconn)) > 0.1) then
          if (cur_connection_set%dist(2,iconn) > 0.1) then
            face_id(iconn) = 3
          elseif (cur_connection_set%dist(2,iconn) < 0.1) then
            face_id(iconn) = 4
          endif
        elseif (dabs(cur_connection_set%dist(3,iconn)) > 0.1) then
          if (cur_connection_set%dist(3,iconn) > 0.1) then
            face_id(iconn) = 5
          else
            face_id(iconn) = 6
          endif
        endif
        cell_z = grid%z(grid%nL2G(cur_connection_set%id_dn(iconn)))
        delta_z = cur_connection_set%dist(0,iconn)* &
                  cur_connection_set%dist(3,iconn)
        face_z = cell_z-delta_z
        face_elevation(iconn) = face_z
        face_projected_area(1,iconn) = cur_connection_set%area(iconn)* &
          dabs(DotProduct(cur_connection_set%dist(1:3,iconn), &
                          gravity_unit_vector))
      enddo
      permutation = permutation - 1
      call PetscSortRealWithPermutation(num_connections,face_elevation, &
                                        permutation,ierr);CHKERRQ(ierr)
      permutation = permutation + 1
      iconn = 1
      iface = permutation(iconn)
      iconn_prev_elev = iconn
      do
        if (iconn > num_connections) exit
        iface = permutation(iconn)
        if (face_elevation(iface) > &
            face_elevation(permutation(iconn_prev_elev))) then
          do i = iconn_prev_elev, iconn-1
            face_projected_area(2,permutation(i)) = sum_area
          enddo
          iconn_prev_elev = iconn
        endif
        sum_area = sum_area + face_projected_area(1,iface)
        iconn = iconn + 1
      enddo
      do i = iconn_prev_elev, num_connections
        face_projected_area(2,permutation(i)) = sum_area
      enddo
      i = 1
      do iconn = 2, num_connections
        if (face_projected_area(2,permutation(iconn)) > &
            face_projected_area(2,permutation(iconn-1))) then
          i = i + 1
        endif
      enddo
      allocate(metric_by_elevation(METRIC_UNDERLYING_VOLUME,i+1))
      metric_by_elevation = UNINITIALIZED_DOUBLE
      i = 1
      metric_by_elevation(METRIC_ELEVATION,i) = grid%z_min_global
      metric_by_elevation(METRIC_AREA,i) = 0.d0
      metric_by_elevation(METRIC_UNDERLYING_VOLUME,i) = 0.d0
      i = i + 1
      metric_by_elevation(METRIC_ELEVATION,i) = &
        face_elevation(permutation(1))
      metric_by_elevation(METRIC_AREA,i) = &
        face_projected_area(2,permutation(1))
      metric_by_elevation(METRIC_UNDERLYING_VOLUME,i) = &
        metric_by_elevation(METRIC_UNDERLYING_VOLUME,i-1) + &
        (metric_by_elevation(METRIC_ELEVATION,i) - &
         metric_by_elevation(METRIC_ELEVATION,i-1)) * &
        metric_by_elevation(METRIC_AREA,i-1)
      do iconn = 2, num_connections
        if (face_projected_area(2,permutation(iconn)) > &
            face_projected_area(2,permutation(iconn-1))) then
          i = i + 1
          metric_by_elevation(METRIC_ELEVATION,i) = &
            face_elevation(permutation(iconn))
          metric_by_elevation(METRIC_AREA,i) = &
            face_projected_area(2,permutation(iconn))
          metric_by_elevation(METRIC_UNDERLYING_VOLUME,i) = &
            metric_by_elevation(METRIC_UNDERLYING_VOLUME,i-1) + &
            (metric_by_elevation(METRIC_ELEVATION,i) - &
            metric_by_elevation(METRIC_ELEVATION,i-1)) * &
            metric_by_elevation(METRIC_AREA,i-1)
        endif
      enddo
      do i = 1, size(metric_by_elevation,2)
        iface = permutation(iconn)
        write(*,'(3f9.5)') metric_by_elevation(:,i)
      enddo
      do iconn = 1, num_connections
        iface = permutation(iconn)
        write(*,'(i4,i2,3f9.5)') cell_id(iface), face_id(iface), &
                face_elevation(iface), face_projected_area(:,iface)
      enddo
      deallocate(cell_id)
      deallocate(face_id)
      deallocate(face_projected_area)
      deallocate(face_elevation)
      deallocate(permutation)
      cur_ponded_water_aux%metric_by_elevation => metric_by_elevation
      nullify(metric_by_elevation)
      if (associated(prev_ponded_water_aux)) then
        prev_ponded_water_aux%next => cur_ponded_water_aux
      else
        this%aux => cur_ponded_water_aux
      endif
      prev_ponded_water_aux => cur_ponded_water_aux
      nullify(cur_ponded_water_aux)
    endif
    boundary_condition => boundary_condition%next
  enddo

  ! update ponded water depth
  boundary_condition => patch%boundary_condition_list%first
  cur_ponded_water_aux => this%aux
  do
    if (.not.associated(boundary_condition)) exit
    if (boundary_condition%flow_bc_type(1) == PONDED_WATER_BC) then
      select case(option%iflowmode)
        case(ZFLOW_MODE)
          water_index = &
            boundary_condition%flow_aux_mapping(ZFLOW_COND_WATER_INDEX)
        case default
          this%option%io_buffer = 'The ponded water process model is &
            &not currently supported by ' // trim(option%flowmode) // '.'
          call PrintErrMsg(this%option)
      end select
      cur_connection_set => boundary_condition%connection_set
      num_connections = cur_connection_set%num_connections
      select type(dataset => &
                    boundary_condition%flow_condition%pressure%dataset)
        class is(dataset_ascii_type)
          cur_ponded_water_aux%pond_elevation = dataset%rarray(1)
          do iconn = 1, num_connections
            cell_z = grid%z(grid%nL2G(cur_connection_set%id_dn(iconn)))
            delta_z = cur_connection_set%dist(0,iconn)* &
                      cur_connection_set%dist(3,iconn)
            face_z = cell_z-delta_z
            boundary_condition%flow_aux_real_var(water_index,iconn) = &
              max(cur_ponded_water_aux%pond_elevation-face_z,0.d0)
          enddo
#if 0
        class is(dataset_gridded_hdf5_type)
          call PatchUpdateCouplerGridDataset(boundary_condition, &
                                             option,patch%grid, &
                                             dataset,water_index)
        class is(dataset_common_hdf5_type)
          ! skip cell indexed datasets used in initial conditions
#endif
        class default
          call PrintMsg(option,'pressure%itype,PONDED_WATER')
          call DatasetUnknownClass(dataset,option, &
                                   'PMPondedWaterInitializeRun')
      end select
      cur_ponded_water_aux => cur_ponded_water_aux%next
    endif
    boundary_condition => boundary_condition%next
  enddo

end subroutine PMPondedWaterInitializeRun

! ************************************************************************** !

subroutine PMPondedWaterUpdate1(this,time,ierr)
  !
  ! Updates the depth of each ponded water boundary condition by column
  !
  ! Author: Glenn Hammond
  ! Date: 02/07/25

  use Patch_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use ZFlow_Aux_module

  class(pm_ponded_water_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(ponded_water_aux_type), pointer :: cur_ponded_water_aux
  PetscInt :: iconn
  PetscInt :: water_index
  PetscInt :: sum_connection

  ierr = 0

  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option

  ! update ponded water depth
  boundary_condition => patch%boundary_condition_list%first
  cur_ponded_water_aux => this%aux
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    if (boundary_condition%flow_bc_type(1) /= PONDED_WATER_BC) then
      sum_connection = sum_connection + cur_connection_set%num_connections
      boundary_condition => boundary_condition%next
      cycle
    endif

    water_index = boundary_condition%flow_aux_mapping(ZFLOW_COND_WATER_INDEX)
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      boundary_condition%flow_aux_real_var(water_index,iconn) = &
        boundary_condition%flow_aux_real_var(water_index,iconn) - &
        patch%boundary_velocities(1,sum_connection)*option%flow_dt
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine PMPondedWaterUpdate1

! ************************************************************************** !

subroutine PMPondedWaterUpdate2(this,time,ierr)
  !
  ! Updates the depth of each ponded water boundary condition, redistributing
  ! mass among columns
  !
  ! Author: Glenn Hammond
  ! Date: 02/07/25

  use Patch_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use ZFlow_Aux_module

  class(pm_ponded_water_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(ponded_water_aux_type), pointer :: cur_ponded_water_aux
  PetscInt :: iconn
  PetscInt :: water_index
  PetscInt :: sum_connection
  PetscReal :: cell_z
  PetscReal :: delta_z
  PetscReal :: face_z
  PetscReal :: sum_flux
  PetscReal :: water_volume
  PetscReal :: pond_elevation

  ierr = 0

  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option

  ! update ponded water depth
  boundary_condition => patch%boundary_condition_list%first
  cur_ponded_water_aux => this%aux
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    if (boundary_condition%flow_bc_type(1) /= PONDED_WATER_BC) then
      sum_connection = sum_connection + cur_connection_set%num_connections
      boundary_condition => boundary_condition%next
      cycle
    endif
    water_volume = PMPondedWaterCalcTotalVolume(cur_ponded_water_aux, &
                                                cur_ponded_water_aux%pond_elevation)
    water_index = boundary_condition%flow_aux_mapping(ZFLOW_COND_WATER_INDEX)
    sum_flux = 0.d0
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      sum_flux = sum_flux + &
                 patch%boundary_velocities(1,sum_connection)*option%flow_dt
    enddo
    water_volume = water_volume + sum_flux
    pond_elevation = &
      PMPondedWaterElevationFromVolume(cur_ponded_water_aux,water_volume)
    do iconn = 1, cur_connection_set%num_connections
      cell_z = grid%z(grid%nL2G(cur_connection_set%id_dn(iconn)))
      delta_z = cur_connection_set%dist(0,iconn)* &
                cur_connection_set%dist(3,iconn)
      face_z = cell_z-delta_z
      boundary_condition%flow_aux_real_var(water_index,iconn) = &
        pond_elevation - face_z
    enddo
    cur_ponded_water_aux%pond_elevation = pond_elevation
    cur_ponded_water_aux => cur_ponded_water_aux%next
    boundary_condition => boundary_condition%next
  enddo

end subroutine PMPondedWaterUpdate2

! ************************************************************************** !

function PMPondedWaterCalcTotalVolume(aux,pond_elevation)
  !
  ! Calculates the current volume of water in the pond
  !
  ! Author: Glenn Hammond
  ! Date: 02/20/25

  type(ponded_water_aux_type) :: aux
  PetscReal :: pond_elevation

  PetscReal :: PMPondedWaterCalcTotalVolume

  PetscInt :: num_elevation
  PetscInt :: ielevation
  PetscReal :: volume

  num_elevation = size(aux%metric_by_elevation,2)

  volume = 0.d0
  do ielevation = 2, num_elevation
    volume = volume + &
      PMPondedWaterCalcIncrVolume(aux%metric_by_elevation,ielevation,pond_elevation)
  enddo

#if 0
  ielevation = aux%ielevation
  if (pond_elevation < &
      aux%metric_by_elevation(METRIC_ELEVATION,ielevation)) then
    do
      ielevation = ielevation - 1
      if (ielevation < 1) then
        ielevation = 1
        exit
      else if (pond_elevation > &
               aux%metric_by_elevation(METRIC_ELEVATION,ielevation)) then
        exit
      endif
    enddo
  else if (pond_elevation > &
           aux%metric_by_elevation(METRIC_ELEVATION,ielevation+1)) then
    do
      ielevation = ielevation + 1
      if (ielevation > num_elevation) then
        ielevation = num_elevation
        exit
      else if (pond_elevation < &
               aux%metric_by_elevation(METRIC_ELEVATION,ielevation)) then
        ielevation = ielevation - 1
        exit
      endif
    enddo
  endif

  PMPondedWaterCalcTotalVolume = &
    aux%metric_by_elevation(METRIC_UNDERLYING_VOLUME,ielevation) + &
    PMPondedWaterCalcIncrVolume(aux%metric_by_elevation,ielevation,elevation)




    else if (pond_elevation > &
             aux%metric_by_elevation(METRIC_ELEVATION, &
                                     min(ielevation+1,num_elevation)) then
      ielevation = min(ielevation+1,num_elevation)
    endif
  enddo

  do
    if (pond_elevation < &
        aux%metric_by_elevation(METRIC_ELEVATION,ielevation)) then
      ielevation = ielevation - 1
    endif
  enddo
  volume = 0.d0
  do i = 2, size(aux%metric_by_elevation,2)
    volume = volume + aux%metric_by_elevation(METRIC_AREA)
  enddo
#endif

  PMPondedWaterCalcTotalVolume = volume

end function PMPondedWaterCalcTotalVolume

! ************************************************************************** !

function PMPondedWaterElevationFromVolume(aux,pond_volume)
  !
  ! Calculates the current volume of water in the pond
  !
  ! Author: Glenn Hammond
  ! Date: 02/20/25

  type(ponded_water_aux_type) :: aux
  PetscReal :: pond_volume

  PetscReal :: PMPondedWaterElevationFromVolume

  PetscInt :: num_elevation
  PetscInt :: ielevation
  PetscReal :: elevation
  PetscReal :: delta_elevation
  PetscReal :: volume

  num_elevation = size(aux%metric_by_elevation,2)

  ielevation = 2
  do
    if (ielevation > num_elevation) exit
    if (pond_volume > &
        aux%metric_by_elevation(METRIC_UNDERLYING_VOLUME,ielevation-1) .and. &
        pond_volume <= &
        aux%metric_by_elevation(METRIC_UNDERLYING_VOLUME,ielevation)) then
      ielevation = ielevation - 1
      exit
    endif
    ielevation = ielevation + 1
  enddo
  ielevation = min(ielevation,num_elevation)
  volume = aux%metric_by_elevation(METRIC_UNDERLYING_VOLUME,ielevation)
  volume = pond_volume - volume
  delta_elevation = volume / aux%metric_by_elevation(METRIC_AREA,ielevation)
  elevation = aux%metric_by_elevation(METRIC_ELEVATION,ielevation) + delta_elevation

  PMPondedWaterElevationFromVolume = elevation

end function PMPondedWaterElevationFromVolume

! ************************************************************************** !

function PMPondedWaterCalcIncrVolume(metric_by_elevation,i,elevation)
  !
  ! Calculates the current volume of water in the pond
  !
  ! Author: Glenn Hammond
  ! Date: 02/20/25

  PetscReal :: metric_by_elevation(:,:)
  PetscInt :: i
  PetscReal :: elevation

  PetscReal :: PMPondedWaterCalcIncrVolume

  PetscReal :: delta_z

  delta_z = max(min(metric_by_elevation(METRIC_ELEVATION,i),elevation) - &
                metric_by_elevation(METRIC_ELEVATION,i-1),0.d0)
  PMPondedWaterCalcIncrVolume = delta_z * metric_by_elevation(METRIC_AREA,i)

end function PMPondedWaterCalcIncrVolume

! ************************************************************************** !

recursive subroutine PMPondedWaterFinalizeRun(this)
  !
  ! Finalizes the run
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/25
  !

  class(pm_ponded_water_type) :: this

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMPondedWaterFinalizeRun

! ************************************************************************** !

recursive subroutine PMPondedWaterAuxDestroy(aux)
  !
  ! Destroys ponded water process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/25

  use Utility_module

  type(ponded_water_aux_type), pointer :: aux

  if (.not.associated(aux)) return
  call PMPondedWaterAuxDestroy(aux%next)
  call DeallocateArray(aux%metric_by_elevation)

end subroutine PMPondedWaterAuxDestroy

! ************************************************************************** !

subroutine PMPondedWaterDestroy(this)
  !
  ! Destroys ponded water process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/25

  class(pm_ponded_water_type) :: this

  call PMBaseDestroy(this)
  call PMPondedWaterAuxDestroy(this%aux)
  ! destroyed in realization
  nullify(this%realization)
  nullify(this%option)

end subroutine PMPondedWaterDestroy

end module PM_Ponded_Water_class
