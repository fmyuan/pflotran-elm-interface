module SrcSink_Sandbox_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Geometry_module

  implicit none

  private

  type, abstract, public :: srcsink_sandbox_base_type
    PetscInt, pointer :: local_cell_ids(:)
    PetscInt, pointer :: natural_cell_ids(:)
    type(point3d_type) :: coordinate
    PetscReal, pointer :: instantaneous_mass_rate(:)
    PetscReal, pointer :: cumulative_mass(:)
    class(srcsink_sandbox_base_type), pointer :: next
  contains
    procedure, public :: ReadInput => SSSandboxBaseRead
    procedure, public :: Setup => SSSandboxBaseSetup
    procedure, public :: Update => SSSandboxBaseUpdate
    procedure, public :: Evaluate => SSSandboxBaseEvaluate
    procedure, public :: Destroy => SSSandboxBaseDestroy
  end type srcsink_sandbox_base_type

  public :: SSSandboxBaseInit, &
            SSSandboxBaseSetup, &
            SSSandboxBaseRead, &
            SSSandboxBaseSelectCase, &
            SSSandboxBaseDestroy

contains

! ************************************************************************** !

subroutine SSSandboxBaseInit(this)

  implicit none

  class(srcsink_sandbox_base_type) :: this

  this%coordinate%x = UNINITIALIZED_DOUBLE
  this%coordinate%y = UNINITIALIZED_DOUBLE
  this%coordinate%z = UNINITIALIZED_DOUBLE
  nullify(this%local_cell_ids)
  nullify(this%natural_cell_ids)
  nullify(this%instantaneous_mass_rate)
  nullify(this%cumulative_mass)
  nullify(this%next)

end subroutine SSSandboxBaseInit

! ************************************************************************** !

subroutine SSSandboxBaseSetup(this,grid,material_auxvars,option)

  use Option_module
  use Grid_module
  use Material_Aux_module, only: material_auxvar_type

  implicit none

  class(srcsink_sandbox_base_type) :: this
  type(grid_type) :: grid
  type(material_auxvar_type) :: material_auxvars(:)
  type(option_type) :: option

  PetscInt :: local_id, natural_id
  PetscInt, allocatable :: local_cell_ids(:)
  PetscInt :: icell, num_local, num_global
  PetscInt :: max_natural, min_natural
  PetscErrorCode :: ierr

  num_local = 0
  num_global = 0
  if (Initialized(this%coordinate%x)) then
    call GridGetLocalIDFromCoordinate(grid,this%coordinate,option,local_id)
    num_local = 1
    allocate(local_cell_ids(num_local))
    local_cell_ids(num_local) = local_id
    num_global = 1
  else if (associated(this%natural_cell_ids)) then
    allocate(local_cell_ids(size(this%natural_cell_ids)))
    local_cell_ids = UNINITIALIZED_INTEGER
    num_global = size(this%natural_cell_ids)
    max_natural = maxval(this%natural_cell_ids)
    min_natural = minval(this%natural_cell_ids)
    do local_id = 1, grid%nlmax
      natural_id = grid%nG2A(grid%nL2G(local_id))
      ! speed up search
      if (natural_id >= min_natural .and. natural_id <= max_natural) then
        do icell = 1, num_global
          if (natural_id == this%natural_cell_ids(icell)) then
            num_local = num_local + 1
            local_cell_ids(num_local) = local_id
            exit
          endif
        enddo
      endif
    enddo
  else
    option%io_buffer = 'Source/sink in SSSandbox not associate with the &
      &domain through either a CELL_ID or COORDINATE.'
    call PrintErrMsg(option)
  endif
  if (num_local > 0) then
    allocate(this%local_cell_ids(num_local))
    this%local_cell_ids = local_cell_ids
  endif
  deallocate(local_cell_ids)

  ! check to ensure that each cell is mapped once
  call MPI_Allreduce(num_local,icell,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)

  if (icell == 0) then
    option%io_buffer = 'No grid cells mapped in SSSandboxBaseSetup.'
    call PrintErrMsg(option)
  else if (icell > num_global) then
    option%io_buffer = 'More grid cells than those listed were mapped &
                       &in SSSandboxBaseSetup.'
    call PrintErrMsg(option)
  else if (icell < num_global) then
    option%io_buffer = 'Less grid cells than those listed were mapped &
                       &in SSSandboxBaseSetup.'
    call PrintErrMsg(option)
  endif

  if (associated(this%natural_cell_ids)) deallocate(this%natural_cell_ids)
  if (num_local > 0) then
    allocate(this%natural_cell_ids(num_local))
    do icell = 1, size(this%local_cell_ids)
      this%natural_cell_ids(icell) = &
        grid%nG2A(grid%nL2G(this%local_cell_ids(icell)))
    enddo
  endif

end subroutine SSSandboxBaseSetup

! ************************************************************************** !

subroutine SSSandboxBaseRead(this,input,option)

  use Option_module
  use Input_Aux_module

  implicit none

  class(srcsink_sandbox_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

end subroutine SSSandboxBaseRead

! ************************************************************************** !

subroutine SSSandboxBaseSelectCase(this,input,option,keyword,found)

  use Option_module
  use Input_Aux_module
  use Geometry_module
  use Utility_module

  implicit none

  class(srcsink_sandbox_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found

  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SOURCE_SINK_SANDBOX'

  found = PETSC_TRUE
  select case(trim(keyword))
    case('COORDINATE')
      call GeometryReadCoordinate(input,option,this%coordinate,error_string)
    case('CELL_ID','CELL_IDS')
      call UtilityReadArray(this%natural_cell_ids,UNINITIALIZED_INTEGER, &
                            trim(error_string)//','//trim(keyword), &
                            input,option)
    case default
      found = PETSC_FALSE
  end select

end subroutine SSSandboxBaseSelectCase

! ************************************************************************** !

subroutine SSSandboxBaseUpdate(this,option)

  use Option_module

  implicit none

  class(srcsink_sandbox_base_type) :: this
  type(option_type) :: option

  if (associated(this%cumulative_mass)) then
    this%cumulative_mass(:) = this%cumulative_mass(:) + &
      option%flow_dt*this%instantaneous_mass_rate(:)
  endif

end subroutine SSSandboxBaseUpdate

! ************************************************************************** !

subroutine SSSandboxBaseEvaluate(this,Residual,Jacobian,compute_derivative, &
                                 material_auxvar,aux_real,option)

  use Option_module
  use Material_Aux_module

  implicit none

  class(srcsink_sandbox_base_type) :: this
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: aux_real(:)

end subroutine SSSandboxBaseEvaluate

! ************************************************************************** !

subroutine SSSandboxBaseDestroy(this)

  use Utility_module

  implicit none

  class(srcsink_sandbox_base_type) :: this

  call DeallocateArray(this%natural_cell_ids)
  call DeallocateArray(this%local_cell_ids)
  call DeallocateArray(this%instantaneous_mass_rate)
  call DeallocateArray(this%cumulative_mass)

end subroutine SSSandboxBaseDestroy

end module SrcSink_Sandbox_Base_class
