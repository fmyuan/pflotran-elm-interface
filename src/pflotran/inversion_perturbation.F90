module Inversion_Perturbation_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Inversion_Subsurface_class

  implicit none

  private

  type, public :: inversion_abc_type
    Vec :: quantity_of_interest_base
    Vec :: base_measurement_vec
    PetscInt :: ndof
    PetscInt :: idof_pert
    PetscReal :: pert
    PetscReal :: base_value
    PetscReal :: tolerance
    PetscInt, pointer :: select_cells(:)
  end type inversion_abc_type

  type, public, extends(inversion_subsurface_type) :: &
                                            inversion_perturbation_type
    type(inversion_abc_type), pointer :: perturbation
  contains
    procedure, public :: Init => InversionPerturbationInit
    procedure, public :: ReadBlock => InversionPerturbationReadBlock
    procedure, public :: Initialize => InversionPerturbationInitialize
    procedure, public :: Step => InversionPerturbationStep
    procedure, public :: ConnectToForwardRun => &
                           InvPerturbationConnectForwardRun
    procedure, public :: CalculateSensitivity => &
                           InvPerturbationCalculateSensitivity
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

  allocate(this%perturbation)
  this%perturbation%quantity_of_interest_base = PETSC_NULL_VEC
  this%perturbation%base_measurement_vec = PETSC_NULL_VEC

  this%perturbation%ndof = 0
  this%perturbation%idof_pert = 0
  this%perturbation%pert = 0.d0
  this%perturbation%base_value = 0.d0
  this%perturbation%tolerance = 1.d-6
  nullify(this%perturbation%select_cells)

  zflow_calc_adjoint = PETSC_FALSE

end subroutine InversionPerturbationInit

! ************************************************************************** !

subroutine InversionPerturbationReadBlock(this,input,option)

  use Input_Aux_module
  use Option_module
  use String_module
  use Utility_module

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
        call InputReadDouble(input,option,this%perturbation%tolerance)
        call InputErrorMsg(input,option,keyword,error_string)
      case('SELECT_CELLS')
        call UtilityReadArray(this%perturbation%select_cells,ZERO_INTEGER, &
                              error_string,input,option)
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
  use Inversion_TS_Aux_module
  use String_module

  class(inversion_perturbation_type) :: this

  PetscErrorCode :: ierr

  call InversionSubsurfInitialize(this)
  call InvForwardAuxDestroyList(this%inversion_aux%inversion_forward_aux, &
                                PETSC_FALSE)
  this%inversion_aux%inversion_forward_aux%store_adjoint = PETSC_FALSE

  if (this%perturbation%idof_pert == 0) then
    if (this%qoi_is_full_vector) then
      if (associated(this%perturbation%select_cells)) then
        this%perturbation%ndof = size(this%perturbation%select_cells)
        if (this%perturbation%ndof > this%realization%patch%grid%nmax) then
          call this%driver%PrintErrMsg('Number of SELECT_CELLS is larger than &
                                      &the problem size: '// &
                    trim(StringWrite(this%perturbation%ndof))//' '// &
                    trim(StringWrite(this%realization%patch%grid%nmax)))
        endif
      else
        this%perturbation%ndof = this%realization%patch%grid%nmax
      endif
    else
      this%perturbation%ndof = size(this%parameters)
    endif
    call VecDuplicate(this%measurement_vec, &
                      this%perturbation%base_measurement_vec, &
                      ierr);CHKERRQ(ierr)
  else
    if (this%qoi_is_full_vector) then
      if (this%perturbation%idof_pert > this%realization%patch%grid%nmax) then
        call this%driver%PrintErrMsg('SELECT_CELLS ID is larger than &
                                    &the problem size: '// &
                      trim(StringWrite(this%perturbation%idof_pert))//' '// &
                      trim(StringWrite(this%realization%patch%grid%nmax)))
      endif
    endif
  endif

  if (Uninitialized(this%iqoi(1)) .and. this%qoi_is_full_vector) then
    call this%driver%PrintErrMsg('Quantity of interest not specified in &
      &InversionPerturbationInitialize.')
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

    this%perturbation%idof_pert = 0
    call this%InitializeForwardRun(option)
    call this%Initialize()
    call this%ConnectToForwardRun()
    call this%ExecuteForwardRun()
    call this%CalculateSensitivity()
    call this%OutputSensitivity('')
    call this%DestroyForwardRun()

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
  use Init_Subsurface_module
  use Material_module
  use Region_module
  use Strata_module
  use String_module

  class(inversion_perturbation_type) :: this

  type(discretization_type), pointer :: discretization
  Vec :: work
  Vec :: natural_vec
  PetscReal, pointer :: vec_ptr(:)
  PetscReal :: rmin, rmax
  PetscInt :: i
  type(material_property_type), pointer :: material_property
  PetscErrorCode :: ierr

  call InvSubsurfConnectToForwardRun(this)

  ! on first pass, store and set thereafter
  if (this%perturbation%quantity_of_interest_base == PETSC_NULL_VEC) then
    call VecDuplicate(this%quantity_of_interest, &
                      this%perturbation%quantity_of_interest_base, &
                      ierr);CHKERRQ(ierr)
  endif
  if (this%perturbation%idof_pert == 0) then
    call VecCopy(this%quantity_of_interest, &
                 this%perturbation%quantity_of_interest_base, &
                 ierr);CHKERRQ(ierr)
  else
    if (this%qoi_is_full_vector) then
      discretization => this%realization%discretization
      work = this%realization%field%work
      call DiscretizationCreateVector(discretization, &
                                      ONEDOF,natural_vec, &
                                      NATURAL,this%realization%option)
      call VecCopy(this%perturbation%quantity_of_interest_base, &
                   this%quantity_of_interest, &
                   ierr);CHKERRQ(ierr)
      call VecGetArrayF90(this%quantity_of_interest,vec_ptr, &
                          ierr);CHKERRQ(ierr)
      call VecZeroEntries(natural_vec,ierr);CHKERRQ(ierr)
      if (this%driver%comm%myrank == 0) then
        call VecSetValue(natural_vec,this%perturbation%idof_pert-1, &
                        this%perturbation%tolerance,INSERT_VALUES, &
                        ierr);CHKERRQ(ierr)
      endif
      call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
      call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
      call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                        work,ONEDOF)
      call VecPointwiseMult(work,work,this%quantity_of_interest, &
                            ierr);CHKERRQ(ierr)
      call VecMax(work,PETSC_NULL_INTEGER,rmax,ierr);CHKERRQ(ierr)
      call VecMin(work,PETSC_NULL_INTEGER,rmin,ierr);CHKERRQ(ierr)
      if (rmax > 0.d0) then
        this%perturbation%pert = rmax
      else
        this%perturbation%pert = rmin
      endif
      call VecAXPY(this%quantity_of_interest,1.d0,work,ierr);CHKERRQ(ierr)
      call DiscretizationGlobalToLocal(this%realization%discretization, &
                                      this%quantity_of_interest, &
                                      this%realization%field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                  this%realization%field%work_loc, &
                                  this%iqoi(1),this%iqoi(2))
      call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
    else
      this%perturbation%base_value = &
        this%parameters(this%perturbation%idof_pert)%value
      this%perturbation%pert = this%perturbation%base_value* &
                               this%perturbation%tolerance
      this%parameters(this%perturbation%idof_pert)%value = &
        this%perturbation%base_value + this%perturbation%pert
      do i = 1, size(this%parameters)
        material_property => &
          MaterialPropGetPtrFromArray(this%parameters(i)%material_name, &
                            this%realization%patch%material_property_array)
        if (this%first_iteration) then
          this%parameters(i)%imat = abs(material_property%internal_id)
          this%parameters(i)%value = material_property%permeability(1,1)
        else
          material_property%permeability(1,1) = this%parameters(i)%value
          material_property%permeability(2,2) = this%parameters(i)%value
          material_property%permeability(3,3) = this%parameters(i)%value
        endif
      enddo
      call InitSubsurfAssignMatIDsToRegns(this%realization)
      call InitSubsurfAssignMatProperties(this%realization)
    endif
  endif

end subroutine InvPerturbationConnectForwardRun

! ************************************************************************** !

subroutine InvPerturbationCalculateSensitivity(this)
  !
  ! Calculates sensitivity matrix Jsensitivity
  !
  ! Author: Glenn Hammond
  ! Date: 03/21/22
  !
  use Option_module

  class(inversion_perturbation_type) :: this

  type(option_type), pointer :: option
  PetscInt :: iteration

  ! destroy non-perturbed forward run
  iteration = 0
  ! InversionPerturbationFillRow performs setup on iteration 0
  call InversionPerturbationFillRow(this,iteration)
  call this%DestroyForwardRun()
  iteration = 1
  do
    if (associated(this%perturbation%select_cells)) then
      this%perturbation%idof_pert = this%perturbation%select_cells(iteration)
    else
      this%perturbation%idof_pert = iteration
    endif
    call this%InitializeForwardRun(option)
    call this%Initialize()
    call this%ConnectToForwardRun()
    call this%ExecuteForwardRun()
    call InversionPerturbationFillRow(this,iteration)
    iteration = iteration + 1
    if (iteration > this%perturbation%ndof) exit
    ! the last forward run will be destroyed after any output of
    ! sensitivity matrices
    call this%DestroyForwardRun()
  enddo

end subroutine InvPerturbationCalculateSensitivity

! ************************************************************************** !

subroutine InversionPerturbationFillRow(this,iteration)
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
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscErrorCode :: ierr


  call VecGetArrayF90(this%measurement_vec,vec_ptr,ierr)
  do i = 1, size(this%measurements)
    vec_ptr(i) = this%measurements(i)%simulated_value
  enddo
  call VecRestoreArrayF90(this%measurement_vec,vec_ptr,ierr)

  if (this%perturbation%idof_pert == 0) then
    call VecCopy(this%measurement_vec,this%perturbation%base_measurement_vec, &
                 ierr);CHKERRQ(ierr)
    call MatZeroEntries(this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
  else
    call VecAXPY(this%measurement_vec,-1.d0, &
                 this%perturbation%base_measurement_vec, &
                 ierr);CHKERRQ(ierr)
    call VecScale(this%measurement_vec,1.d0/this%perturbation%pert, &
                  ierr);CHKERRQ(ierr)
  endif

  if (this%perturbation%idof_pert == 0) return

  ! don't need to use the distributed vec, but why not
  call VecScatterBegin(this%scatter_measure_to_dist_measure, &
                       this%measurement_vec,this%dist_measurement_vec, &
                       INSERT_VALUES,SCATTER_FORWARD_LOCAL, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(this%scatter_measure_to_dist_measure, &
                     this%measurement_vec,this%dist_measurement_vec, &
                     INSERT_VALUES,SCATTER_FORWARD_LOCAL, &
                     ierr);CHKERRQ(ierr)
  call VecGetArrayF90(this%dist_measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)
  do i = 1, size(vec_ptr)
    call MatSetValue(this%inversion_aux%JsensitivityT,this%perturbation%idof_pert-1, &
                     this%dist_measurement_offset+i-1,vec_ptr(i), &
                     INSERT_VALUES,ierr);CHKERRQ(ierr)
  enddo
  call VecRestoreArrayF90(this%dist_measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)

  if (iteration == this%perturbation%ndof) then
    call MatAssemblyBegin(this%inversion_aux%JsensitivityT,MAT_FINAL_ASSEMBLY, &
                          ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(this%inversion_aux%JsensitivityT,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  endif

  if (.not.this%qoi_is_full_vector) then
    ! revert back to base value
    this%parameters(this%perturbation%idof_pert)%value = this%perturbation%base_value
  endif

end subroutine InversionPerturbationFillRow

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

  call DeallocateArray(this%perturbation%select_cells)
  if (this%perturbation%quantity_of_interest_base /= PETSC_NULL_VEC) then
    call VecDestroy(this%perturbation%quantity_of_interest_base, &
                    ierr);CHKERRQ(ierr)
  endif
  if (this%perturbation%base_measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%perturbation%base_measurement_vec,ierr);CHKERRQ(ierr)
  endif
  deallocate(this%perturbation)
  nullify(this%perturbation)

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
