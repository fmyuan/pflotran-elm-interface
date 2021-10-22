module Inversion_Perturbation_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Inversion_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_subsurface_type) :: &
                                            inversion_perturbation_type
    Mat :: Jsensitivity_pert
    Vec :: quantity_of_interest_base
    PetscInt :: ndof
    PetscInt :: idof_pert
    PetscReal :: pert
    PetscReal :: perturbation_tolerance
  contains
    procedure, public :: Init => InversionPerturbationInit
    procedure, public :: ReadBlock => InversionPerturbationReadBlock
    procedure, public :: Initialize => InversionPerturbationInitialize
    procedure, public :: Step => InversionPerturbationStep
    procedure, public :: ConnectToForwardRun => InvPerturbationConnectForwardRun
    procedure, public :: Finalize => InversionPerturbationFinalize
    procedure, public :: Strip => InversionPerturbationStrip
  end type inversion_perturbation_type

  public :: InversionPerturbationCreate, &
            InversionPerturbationFinalize, &
            InversionPerturbationStrip, &
            InversionPerturbationDestroy

contains

! ************************************************************************** !

function InversionPerturbationCreate(driver)
  !
  ! Allocates and initializes a new perturbation inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Driver_module

  class(driver_type), pointer :: driver

  class(inversion_perturbation_type), pointer :: InversionPerturbationCreate

  allocate(InversionPerturbationCreate)
  call InversionPerturbationCreate%Init(driver)

end function InversionPerturbationCreate

! ************************************************************************** !

subroutine InversionPerturbationInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Driver_module
  use ZFlow_Aux_module, only : zflow_calc_adjoint

  class(inversion_perturbation_type) :: this
  class(driver_type), pointer :: driver

  call InversionSubsurfaceInit(this,driver)

  this%quantity_of_interest_base = PETSC_NULL_VEC
  this%Jsensitivity_pert = PETSC_NULL_MAT

  this%ndof = 0
  this%idof_pert = 0
  this%pert = 0.d0
  this%perturbation_tolerance = 1.d-6

  nullify(this%measurement)
  nullify(this%imeasurement)

  zflow_calc_adjoint = PETSC_FALSE

end subroutine InversionPerturbationInit

! ************************************************************************** !

subroutine InversionPerturbationReadBlock(this,input,option)

  use Input_Aux_module
  use Option_module
  use String_module

  class(inversion_perturbation_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'Perturbation Inversion'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_TRUE
    call InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('PERTURBATION_TOLERANCE')
        call InputReadDouble(input,option,this%perturbation_tolerance)
        call InputErrorMsg(input,option,keyword,error_string)
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine InversionPerturbationReadBlock

! ************************************************************************** !

subroutine InversionPerturbationInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  class(inversion_perturbation_type) :: this

  PetscInt :: num_measurement

  call InversionBaseInitialize(this)

  if (Uninitialized(this%iqoi)) then
    call this%driver%PrintErrMsg('Quantity of interest not specified in &
      InversionPerturbationInitialize.')
  endif

end subroutine InversionPerturbationInitialize

! ************************************************************************** !

subroutine InversionPerturbationStep(this)
  !
  ! Execute a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21

  use Option_module
  use Factory_Forward_module
  use String_module

  class(inversion_perturbation_type) :: this

  type(option_type), pointer :: option
  PetscInt :: iteration

  iteration = 0
  do
    this%idof_pert = iteration
    option => OptionCreate()
    option%group_prefix = 'Run' // trim(StringWrite(this%iteration)) // &
                          '_' // StringWrite(iteration)
    call OptionSetDriver(option,this%driver)
    call OptionSetInversionOption(option,this%inversion_option)
    call FactoryForwardInitialize(this%forward_simulation, &
                                  this%forward_simulation_filename,option)
    this%realization => this%forward_simulation%realization
    call this%Initialize()
    call this%forward_simulation%InitializeRun()
    call this%ConnectToForwardRun()
    if (option%status == PROCEED) then
      call this%forward_simulation%ExecuteRun()
    endif
    call InversionPerturbationFillColumn(this,iteration)
    if (iteration == this%ndof) then
      call this%OutputSensitivity('pert')
    endif
    call this%forward_simulation%FinalizeRun()
    call this%forward_simulation%Strip()
    deallocate(this%forward_simulation)
    nullify(this%forward_simulation)
    iteration = iteration + 1
    if (iteration > this%ndof) exit
  enddo

  this%converg_flag = PETSC_FALSE
  if (this%iteration > this%maximum_iteration) this%converg_flag = PETSC_TRUE

end subroutine InversionPerturbationStep

! ************************************************************************** !

subroutine InvPerturbationConnectForwardRun(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Discretization_module
  use Material_module

  class(inversion_perturbation_type) :: this

  PetscReal, pointer :: vec_ptr(:)
  PetscInt, parameter :: isubdof = ZERO_INTEGER
  PetscErrorCode :: ierr

  call InvSubsurfConnectToForwardRun(this)

  ! on first pass, store and set thereafter
  if (this%idof_pert == 0) then
    call MatDuplicate(this%inversion_aux%Jsensitivity, &
                      MAT_SHARE_NONZERO_PATTERN, &
                      this%Jsensitivity_pert,ierr);CHKERRQ(ierr)
    call VecGetSize(this%realization%field%work,this%ndof,ierr);CHKERRQ(ierr)
    call VecDuplicate(this%quantity_of_interest, &
                      this%quantity_of_interest_base,ierr);CHKERRQ(ierr)
    call VecCopy(this%quantity_of_interest,this%quantity_of_interest_base, &
                                     ierr);CHKERRQ(ierr)
  else
    call VecCopy(this%quantity_of_interest_base,this%quantity_of_interest, &
                 ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%quantity_of_interest,vec_ptr, &
                        ierr);CHKERRQ(ierr)
    this%pert = this%perturbation_tolerance*vec_ptr(this%idof_pert)
    vec_ptr(this%idof_pert) = vec_ptr(this%idof_pert) + this%pert
    call VecRestoreArrayF90(this%quantity_of_interest,vec_ptr, &
                            ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(this%realization%discretization, &
                                     this%quantity_of_interest, &
                                     this%realization%field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi,isubdof)
  endif

end subroutine InvPerturbationConnectForwardRun

! ************************************************************************** !

subroutine InversionPerturbationFillColumn(this,iteration)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Debug_module
  use Variables_module, only : LIQUID_PRESSURE
  use Realization_Base_class
  use String_module

  class(inversion_perturbation_type) :: this
  PetscInt :: iteration

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: soln_ptr(:)
  PetscReal, allocatable :: temp_array(:)
  PetscInt, allocatable :: rows(:)
  PetscInt :: cols(1)
  PetscInt :: num_measurement
  PetscInt :: i
  Mat :: M 
  PetscErrorCode :: ierr

  num_measurement = size(this%measurement)
  call RealizationGetVariable(this%realization, &
                              this%realization%field%work, &
                              LIQUID_PRESSURE,ZERO_INTEGER)

  call VecGetArrayReadF90(this%realization%field%work,soln_ptr,ierr);CHKERRQ(ierr)
  if (iteration == 0) then
    do i = 1, num_measurement
      this%measurement(i) = soln_ptr(this%imeasurement(i))
    enddo
    call MatZeroEntries(this%Jsensitivity_pert,ierr);CHKERRQ(ierr)
  else
    allocate(temp_array(num_measurement))
    print *, 'pert: ', this%pert
    print *, 'meas: ', this%measurement(:)
    print *, 'soln: ', soln_ptr(:)
    temp_array = 0.d0
    do i = 1, num_measurement
      temp_array(i) = soln_ptr(this%imeasurement(i)) - this%measurement(i)
    enddo
    temp_array = temp_array / this%pert
  endif
  call VecRestoreArrayReadF90(this%realization%field%work,soln_ptr,ierr);CHKERRQ(ierr)

  if (iteration == 0) return

  allocate(rows(num_measurement))
  do i = 1, num_measurement
    rows(i) = i-1
  enddo
  cols = iteration-1
  call MatSetValues(this%Jsensitivity_pert, &
                    num_measurement,rows,1,cols, &
                    temp_array,INSERT_VALUES,ierr);CHKERRQ(ierr)
  deallocate(rows)
  deallocate(temp_array)

  if (iteration == this%ndof) then
    call MatAssemblyBegin(this%Jsensitivity_pert,MAT_FINAL_ASSEMBLY, &
                          ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(this%Jsensitivity_pert,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
    ! on last iteration, have to point sensitivity matrix to this one in
    ! order to print the matrix at the end; avoids having  so swap
    M = this%Jsensitivity_pert
    this%Jsensitivity_pert = this%inversion_aux%Jsensitivity
    this%inversion_aux%Jsensitivity = M
  endif

end subroutine InversionPerturbationFillColumn

! ************************************************************************** !

subroutine InversionPerturbationFinalize(this)
  !
  ! Finalizes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  class(inversion_perturbation_type) :: this

  call InversionBaseStrip(this)

end subroutine InversionPerturbationFinalize

! ************************************************************************** !

subroutine InversionPerturbationStrip(this)
  !
  ! Deallocates members of inversion perturbation
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  class(inversion_perturbation_type) :: this

  PetscErrorCode :: ierr

   call InversionSubsurfaceStrip(this)

   call MatDestroy(this%Jsensitivity_pert,ierr);CHKERRQ(ierr)
   call VecDestroy(this%quantity_of_interest,ierr);CHKERRQ(ierr)
   call VecDestroy(this%quantity_of_interest_base,ierr);CHKERRQ(ierr)

end subroutine InversionPerturbationStrip

! ************************************************************************** !

subroutine InversionPerturbationDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  class(inversion_perturbation_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionPerturbationDestroy

end module Inversion_Perturbation_class
