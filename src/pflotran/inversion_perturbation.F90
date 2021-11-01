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
    Vec :: quantity_of_interest_base
    Vec :: natural_vec
    PetscInt :: ndof
    PetscInt :: idof_pert
    PetscReal :: pert
    PetscReal :: perturbation_tolerance
    PetscReal, pointer :: base_solution_measurement(:)
  contains
    procedure, public :: Init => InversionPerturbationInit
    procedure, public :: ReadBlock => InversionPerturbationReadBlock
    procedure, public :: Initialize => InversionPerturbationInitialize
    procedure, public :: Step => InversionPerturbationStep
    procedure, public :: ConnectToForwardRun => &
                           InvPerturbationConnectForwardRun
    procedure, public :: Strip => InversionPerturbationStrip
  end type inversion_perturbation_type

  public :: InversionPerturbationCreate, &
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
  this%natural_vec = PETSC_NULL_VEC

  this%ndof = 0
  this%idof_pert = 0
  this%pert = 0.d0
  this%perturbation_tolerance = 1.d-6

  nullify(this%base_solution_measurement)

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

  call InversionSubsurfInitialize(this)

  if (this%idof_pert == 0) then
    allocate(this%base_solution_measurement(size(this%imeasurement)))
    this%base_solution_measurement = UNINITIALIZED_DOUBLE
    this%ndof = this%realization%patch%grid%nmax
  endif

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
      call this%OutputSensitivity('')
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

  type(discretization_type), pointer :: discretization
  Vec :: work
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, parameter :: isubdof = ZERO_INTEGER
  PetscReal :: rmin, rmax
  PetscErrorCode :: ierr

  call InvSubsurfConnectToForwardRun(this)

  ! on first pass, store and set thereafter
  if (this%quantity_of_interest_base == PETSC_NULL_VEC) then
    call VecDuplicate(this%quantity_of_interest, &
                      this%quantity_of_interest_base,ierr);CHKERRQ(ierr)
  endif
  if (this%idof_pert == 0) then
    call VecCopy(this%quantity_of_interest,this%quantity_of_interest_base, &
                                     ierr);CHKERRQ(ierr)
  else
    discretization => this%realization%discretization
    work = this%realization%field%work
    call DiscretizationCreateVector(discretization, &
                                    ONEDOF,this%natural_vec, &
                                    NATURAL,this%realization%option)
    call VecCopy(this%quantity_of_interest_base,this%quantity_of_interest, &
                 ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%quantity_of_interest,vec_ptr, &
                        ierr);CHKERRQ(ierr)
    call VecZeroEntries(this%natural_vec,ierr);CHKERRQ(ierr)
    if (this%driver%comm%myrank == 0) then
      call VecSetValue(this%natural_vec,this%idof_pert-1, &
                       this%perturbation_tolerance,INSERT_VALUES, &
                       ierr);CHKERRQ(ierr)
    endif
    call VecAssemblyBegin(this%natural_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(this%natural_vec,ierr);CHKERRQ(ierr)
    call DiscretizationNaturalToGlobal(discretization,this%natural_vec, &
                                       work,ONEDOF)
    call VecPointwiseMult(work,work,this%quantity_of_interest,ierr);CHKERRQ(ierr)
    call VecMax(work,PETSC_NULL_INTEGER,rmax,ierr);CHKERRQ(ierr)
    call VecMin(work,PETSC_NULL_INTEGER,rmin,ierr);CHKERRQ(ierr)
    if (rmax > 0.d0) then
      this%pert = rmax
    else
      this%pert = rmin
    endif
    call VecAXPY(this%quantity_of_interest,1.d0,work,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(this%realization%discretization, &
                                     this%quantity_of_interest, &
                                     this%realization%field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi,isubdof)
    call VecDestroy(this%natural_vec,ierr);CHKERRQ(ierr)
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
  PetscInt, allocatable :: cols(:)
  PetscInt :: rows(1)
  PetscInt :: num_measurement
  PetscInt :: i
  Mat :: M
  PetscErrorCode :: ierr

  num_measurement = size(this%imeasurement)
  call RealizationGetVariable(this%realization, &
                              this%realization%field%work, &
                              LIQUID_PRESSURE,ZERO_INTEGER)

  call VecGetArrayReadF90(this%realization%field%work,soln_ptr, &
                          ierr);CHKERRQ(ierr)
  if (iteration == 0) then
    do i = 1, num_measurement
      this%base_solution_measurement(i) = soln_ptr(this%imeasurement(i))
    enddo
    call MatZeroEntries(this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
  else
    allocate(temp_array(num_measurement))
    print *, 'pert: ', this%pert
    print *, 'meas: ', this%measurement(:)
    print *, 'soln: ', soln_ptr(:)
    temp_array = 0.d0
    do i = 1, num_measurement
      temp_array(i) = soln_ptr(this%imeasurement(i)) - &
                      this%base_solution_measurement(i)
    enddo
    temp_array = temp_array / this%pert
  endif
  call VecRestoreArrayReadF90(this%realization%field%work,soln_ptr, &
                              ierr);CHKERRQ(ierr)

  if (iteration == 0) return

  allocate(cols(num_measurement))
  do i = 1, num_measurement
    cols(i) = i-1
  enddo
  rows = iteration-1
  ! remember the Jacobian is a transpose
  ! icol = measurement
  ! irow = parameter
  call MatSetValues(this%inversion_aux%JsensitivityT, &
                    1,rows,num_measurement,cols, &
                    temp_array,INSERT_VALUES,ierr);CHKERRQ(ierr)
  deallocate(cols)
  deallocate(temp_array)

  if (iteration == this%ndof) then
    call MatAssemblyBegin(this%inversion_aux%JsensitivityT,MAT_FINAL_ASSEMBLY, &
                          ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(this%inversion_aux%JsensitivityT,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  endif

end subroutine InversionPerturbationFillColumn

! ************************************************************************** !

subroutine InversionPerturbationStrip(this)
  !
  ! Deallocates members of inversion perturbation
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Utility_module

  class(inversion_perturbation_type) :: this

  PetscErrorCode :: ierr

   call InversionSubsurfaceStrip(this)

   call DeallocateArray(this%base_solution_measurement)
   if (this%quantity_of_interest_base /= PETSC_NULL_VEC) then
    call VecDestroy(this%quantity_of_interest_base,ierr);CHKERRQ(ierr)
   endif
   if (this%natural_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%natural_vec,ierr);CHKERRQ(ierr)
   endif

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
