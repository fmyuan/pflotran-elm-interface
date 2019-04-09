module Criticality_module
! 
! Author: Michael Nole
! Date: 11/01/18

! MODULE DESCRIPTION:
! ===========================================================================
! This module calculates the radionuclide and energy source terms due 
! to waste form criticality.
! ===========================================================================
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Data_Mediator_Vec_class
  use Region_module

  implicit none
  
  private
  
  ! Stores variables relevant to criticality calculations
  type, public :: criticality_mechanism_type
    character(len=MAXWORDLENGTH) :: mech_name
    PetscReal :: heat_released
    PetscReal :: sw
    PetscReal :: rho_w
    PetscReal :: temperature
    PetscReal :: k_effective
    type(criticality_mechanism_type), pointer :: next
  end type criticality_mechanism_type
  
  ! Stores information regarding the criticality event
  type, public :: criticality_event_type
    character(len=MAXWORDLENGTH) :: mech_name
    PetscBool :: steady_state
    PetscReal :: crit_time
    PetscBool :: crit_flag
  end type criticality_event_type

  ! Criticality process model object. Extends the pm_waste_form type to
  ! include relevant variables for criticality consequence calculations.
  type, public :: criticality_type
    type(criticality_event_type), pointer :: crit_event
    type(criticality_mechanism_type), pointer :: crit_mech
    type(region_type), pointer :: region
    type(criticality_type), pointer :: next
  end type criticality_type
  
  type, public :: criticality_mediator_type
    type(criticality_type), pointer :: criticality_list
    type(criticality_mechanism_type), pointer :: crit_mech_list
    class(data_mediator_vec_type), pointer :: data_mediator
    PetscInt :: total_num_cells
  end type criticality_mediator_type
  
! -------------------------------------------------------------------
  public :: CriticalityMediatorCreate, &
            CriticalityMechCreate, &
            CriticalityCreate, &
            ReadCriticalityMech, &
            CriticalityCalc, &
            CriticalityInitializeRun, &
            AssignCritMech, &
            CriticalitySolve

! -------------------------------------------------------------------

contains

subroutine CriticalityInit(this)
  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none
  
  type(criticality_type), pointer :: this
  
  allocate(this)
  allocate(this%crit_event)
  allocate(this%crit_mech)
  nullify(this%region)
  nullify(this%next)
  
  this%crit_event%steady_state = PETSC_FALSE
  this%crit_event%crit_flag = PETSC_FALSE
  this%crit_event%crit_time = UNINITIALIZED_DOUBLE
  
  
  
end subroutine CriticalityInit

! ************************************************************************** !
subroutine CriticalityMechInit(this)
  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none
  
  type(criticality_mechanism_type), pointer :: this
  
  allocate(this)
  nullify(this%next)
  
  this%heat_released = UNINITIALIZED_DOUBLE
  this%sw = UNINITIALIZED_DOUBLE
  this%rho_w = UNINITIALIZED_DOUBLE
  this%temperature = UNINITIALIZED_DOUBLE
  this%k_effective = UNINITIALIZED_DOUBLE
  
  
  
end subroutine CriticalityMechInit
! ************************************************************************** !

subroutine CriticalityMediatorInit(this)
  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18

  implicit none
  
  type(criticality_mediator_type), pointer :: this
  
  nullify(this%data_mediator)
  nullify(this%criticality_list)
  nullify(this%crit_mech_list)
  
  
end subroutine CriticalityMediatorInit

! ************************************************************************** !

function CriticalityMediatorCreate()

  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18
  
  implicit none
  
  type(criticality_mediator_type), pointer :: CriticalityMediatorCreate
  type(criticality_mediator_type), pointer :: crit
  
  allocate(crit)
  call CriticalityMediatorInit(crit)
  
  CriticalityMediatorCreate => crit
  
end function CriticalityMediatorCreate


! ************************************************************************** !

function CriticalityCreate()

  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18
  
  implicit none
  
  type(criticality_type), pointer :: CriticalityCreate
  type(criticality_type), pointer :: crit
  
  allocate(crit)
  call CriticalityInit(crit)
  
  CriticalityCreate => crit
  
end function CriticalityCreate

! ************************************************************************** !

function CriticalityMechCreate()

  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18
  
  implicit none
  
  type(criticality_mechanism_type), pointer :: CriticalityMechCreate
  type(criticality_mechanism_type), pointer :: crit
  
  allocate(crit)
  call CriticalityMechInit(crit)
  
  CriticalityMechCreate => crit
  
end function CriticalityMechCreate

! ************************************************************************** !

subroutine ReadCriticalityMech(this,input,option,keyword,error_string,found)

  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18
  
  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  type(criticality_mediator_type), pointer :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  
  character(len=MAXWORDLENGTH) :: word
  type(criticality_mechanism_type), pointer :: new_crit_mech, cur_crit_mech
  
  error_string = trim(error_string) // ',CRITICALITY'
  found = PETSC_TRUE
  
  select case(trim(keyword))
    case('CRITICALITY_MECH')
      allocate(new_crit_mech)
      new_crit_mech => CriticalityMechCreate()
      do
        call InputReadPflotranString(input, option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case (trim(word))
          case('NAME')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option, &
                  'criticality mechanism assignment',error_string)
            call StringToUpper(word)
            new_crit_mech%mech_name = trim(word)
          case('HEAT_RELEASED')
            call InputReadDouble(input,option,new_crit_mech%heat_released)
            call InputErrorMsg(input,option,'HEAT_RELEASED',error_string)
        end select
      enddo      
      if (.not. associated(this%crit_mech_list)) then
        this%crit_mech_list => new_crit_mech
      else
        cur_crit_mech => this%crit_mech_list
        do
          if (.not. associated(cur_crit_mech)) exit
          if (.not. associated(cur_crit_mech%next)) then
            cur_crit_mech => new_crit_mech
          endif
          cur_crit_mech => cur_crit_mech%next
        enddo
      endif
      nullify(cur_crit_mech)
    case default
      found = PETSC_FALSE
  end select   

end subroutine ReadCriticalityMech

! ************************************************************************** !

subroutine CriticalityCalc(this,time,ierr)

  ! Calculate mass and heat source terms as a function of time.
  ! Author: Michael Nole
  ! Date: 11/01/18
  
  implicit none
  
  type(criticality_mechanism_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  this%heat_released = 1.e-4
  
end subroutine CriticalityCalc

! ************************************************************************** !

subroutine CriticalityInitializeRun(this, realization, option)

  ! Author: Michael Nole
  ! Date: 11/01/18
  
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Subsurface_class
  use Realization_Base_class
  use Option_module
  
  implicit none
  
  type(criticality_mediator_type), pointer :: this
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  
  type(criticality_type), pointer :: cur_criticality
  PetscInt :: vec_size, i, j
  PetscInt, allocatable :: energy_indices_in_residual(:)
  PetscErrorCode :: ierr
  IS :: is

  call RealizCreateFlowMassTransferVec(realization)
  this%data_mediator => DataMediatorVecCreate()
  call this%data_mediator%AddToList(realization%flow_data_mediator_list)
  
  cur_criticality => this%criticality_list
  vec_size = 0
  
  do
    if (.not. associated(cur_criticality)) exit
    vec_size = vec_size + cur_criticality%region%num_cells
    cur_criticality => cur_criticality%next
  enddo
  
  call VecCreateSeq(PETSC_COMM_SELF, vec_size,this%data_mediator%vec, &
                    ierr);CHKERRQ(ierr)
  call VecSetFromOptions(this%data_mediator%vec,ierr); CHKERRQ(ierr)
  
  cur_criticality => this%criticality_list
  allocate(energy_indices_in_residual(vec_size))
  j = 0
  do
    if (.not. associated(cur_criticality)) exit
      do i = 1, cur_criticality%region%num_cells
        j = j + 1
        energy_indices_in_residual(j) = (cur_criticality%region% &
                                        cell_ids(i) - 1) * option%nflowdof + 3
      enddo
    cur_criticality => cur_criticality%next
  enddo
  
  this%total_num_cells = j
  
  call ISCreateGeneral(option%mycomm,vec_size, &
                       energy_indices_in_residual, &
                       PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)
  call VecScatterCreate(this%data_mediator%vec,PETSC_NULL_IS, &
                        realization%field%flow_r, is, &
                        this%data_mediator%scatter_ctx, ierr); CHKERRQ(ierr)
  if (allocated(energy_indices_in_residual)) then 
      deallocate(energy_indices_in_residual)
  endif
  
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  
end subroutine CriticalityInitializeRun

! ************************************************************************** !
subroutine AssignCritMech(this)
  
  use String_module
  
  implicit none
  
  type(criticality_mediator_type), pointer :: this
  
  type(criticality_mechanism_type), pointer :: cur_mechanism
  type(criticality_type), pointer :: cur_criticality
  
  cur_criticality => this%criticality_list
  do
    if(.not. associated(cur_criticality)) exit
    cur_mechanism => this%crit_mech_list
    do
      if (.not. associated(cur_mechanism)) exit
      if (StringCompare(cur_criticality%crit_event%mech_name, &
          cur_mechanism%mech_name)) then
        cur_criticality%crit_mech => cur_mechanism
        exit
      endif
    enddo
    cur_criticality => cur_criticality%next
  enddo
  
end subroutine AssignCritMech

! ************************************************************************** !


subroutine CriticalityStrip(this)

  ! 
  ! Author: Michael Nole
  ! Date: 11/01/18
  implicit none
  
  type(criticality_mediator_type), pointer :: this
  nullify(this)
  

end subroutine CriticalityStrip

! ************************************************************************** !

subroutine CriticalitySolve(this,realization,time,ierr)
  !
  !Author: Michael Nole
  !Date: 11/05/18
  !
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Subsurface_class
  
  implicit none
  
  type(criticality_mediator_type), pointer :: this
  class(realization_subsurface_type), pointer :: realization
  PetscReal :: time
  PetscErrorCode :: ierr
  
  PetscInt :: i,j
  type(criticality_type), pointer :: cur_criticality
  PetscReal, pointer :: heat_source(:)
  
  call VecGetArrayF90(this%data_mediator%vec,heat_source, &
                      ierr);CHKERRQ(ierr)
  
  cur_criticality => this%criticality_list
  j = 0
  do
    if (.not. associated(cur_criticality)) exit
    do i = 1, cur_criticality%region%num_cells
      j = j + 1
      heat_source(j) = -cur_criticality%crit_mech%heat_released
    enddo
    cur_criticality => cur_criticality%next
  enddo
  
  call VecRestoreArrayF90(this%data_mediator%vec,heat_source, &
                          ierr);CHKERRQ(ierr)
  
end subroutine CriticalitySolve

! ************************************************************************** !

end module Criticality_module
