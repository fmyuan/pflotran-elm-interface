module Inversion_Test_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Inversion_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_subsurface_type) :: inversion_test_type
    PetscBool :: store_adjoint_matrices
    Mat :: Jsensitivity_test
  contains
  procedure, public :: Init => InversionTestInit
  procedure, public :: ReadBlock => InversionTestReadBlock
    procedure, public :: CalculateSensitivity => InvTestCalculateSensitivity
!    procedure, public :: Invert => InversionTestInvert
  end type inversion_test_type

  public :: InversionTestCreate, &
            InversionTestDestroy

contains

! ************************************************************************** !

function InversionTestCreate(driver)
  !
  ! Allocates and initializes a new test inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Driver_module

  class(driver_type), pointer :: driver

  class(inversion_test_type), pointer :: InversionTestCreate

  allocate(InversionTestCreate)
  call InversionTestCreate%Init(driver)

end function InversionTestCreate

! ************************************************************************** !

subroutine InversionTestInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  use Driver_module

  class(inversion_test_type) :: this
  class(driver_type), pointer :: driver

  call InversionSubsurfaceInit(this,driver)

  this%store_adjoint_matrices = PETSC_FALSE

end subroutine InversionTestInit

! ************************************************************************** !

subroutine InversionTestReadBlock(this,input,option)

  use Input_Aux_module
  use Option_module
  use String_module

  class(inversion_test_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'Test Inversion'

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
      case('STORE_ADJOINT_MATRICES')
        this%store_adjoint_matrices = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine InversionTestReadBlock

! ************************************************************************** !

subroutine InversionTestInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Inversion_Aux_module

  class(inversion_test_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr

  if (this%store_adjoint_matrices) then
    allocate(this%inversion_aux%dMdK(this%realization%patch%grid%nmax))
    this%inversion_aux%dMdK = PETSC_NULL_MAT
    allocate(this%inversion_aux%dbdK(this%realization%patch%grid%nmax))
    this%inversion_aux%dbdK = PETSC_NULL_VEC
    do i = 1, size(this%inversion_aux%dMdK)
      call MatDuplicate(this%forward_simulation%flow_process_model_coupler% &
                          timestepper%solver%M, &
                        MAT_SHARE_NONZERO_PATTERN, &
                        this%inversion_aux%dMdK(i),ierr);CHKERRQ(ierr)
    enddo
    do i = 1, size(this%inversion_aux%dbdK)
      call VecDuplicate(this%realization%field%work, &
                        this%inversion_aux%dbdK(i),ierr);CHKERRQ(ierr)
    enddo
  endif

end subroutine InversionTestInitialize

! ************************************************************************** !

subroutine InvTestCalculateSensitivity(this)
  !
  ! Calculates sensitivity matrix Jsensitivity
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Connection_module
  use Debug_module
  use Grid_module
  use Inversion_Aux_module
  use Option_module
  use Patch_module
  use Realization_Base_class
  use Solver_module
  use String_module

  use PM_Base_class
  use PM_ZFlow_class

  class(inversion_test_type) :: this

  type(connection_set_type), pointer :: connection_set
  type(grid_type), pointer :: grid
  type(inversion_aux_type), pointer :: inversion_aux
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(solver_type), pointer :: solver
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, allocatable :: flux_coef(:)
  PetscReal :: tempreal
  PetscInt :: i, j, k
  PetscInt :: icell_measurement
  PetscInt :: local_id, local_id_up, local_id_dn
  PetscInt :: iconn
  PetscReal :: hTdMdKTlambda, dbdKTlambda
  Vec :: work
  Vec :: solution
  Vec :: p, lambda
  Mat :: M, Pmat
  Mat :: dMdK
  Vec :: dbdK
  PetscErrorCode :: ierr

  solver => this%forward_simulation%flow_process_model_coupler% &
              timestepper%solver
  option => this%realization%option
  patch => this%realization%patch
  grid => patch%grid
  connection_set => grid%internal_connection_set_list%first
  inversion_aux => this%inversion_aux

  allocate(flux_coef(size(inversion_aux%dFluxdIntConn,1)))
  work = this%realization%field%work ! DO NOT DESTROY!
  solution = this%realization%field%flow_xx ! DO NOT DESTROY!
  call VecDuplicate(work,p,ierr);CHKERRQ(ierr)
  call VecDuplicate(work,lambda,ierr);CHKERRQ(ierr)

  call MatDuplicate(solver%M,MAT_SHARE_NONZERO_PATTERN,dMdK, &
                    ierr);CHKERRQ(ierr)
  call VecDuplicate(work,dbdK,ierr);CHKERRQ(ierr)

  call MatZeroEntries(inversion_aux%Jsensitivity,ierr);CHKERRQ(ierr)
  if (this%debug_adjoint) then
    print *, 'solution'
    call VecView(solution,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(solution,vec_ptr,ierr);CHKERRQ(ierr)
    print *, vec_ptr(:)
    call VecRestoreArrayF90(solution,vec_ptr,ierr);CHKERRQ(ierr)
    print *, 'residual'
    call VecView(this%realization%field%flow_r,PETSC_VIEWER_STDOUT_WORLD, &
                  ierr);CHKERRQ(ierr)
    print *, 'J'
    call MatView(solver%M,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
    call KSPGetOperators(solver%ksp,M,Pmat,ierr);CHKERRQ(ierr)
    print *, 'M'
    call MatView(M,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
  endif
  do j = 1, size(this%imeasurement)
    call VecZeroEntries(p,ierr);CHKERRQ(ierr)
    tempreal = 1.d0
    icell_measurement = this%imeasurement(j)
    call VecSetValue(p,icell_measurement-1,tempreal, &
                      INSERT_VALUES,ierr);CHKERRQ(ierr)
    if (this%debug_adjoint) then
      print *, 'p'
      call VecView(p,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
    endif
    call KSPSolveTranspose(solver%ksp,p,lambda,ierr);CHKERRQ(ierr)
    if (this%debug_adjoint) then
      print *, 'lambda'
      call VecView(lambda,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(lambda,vec_ptr,ierr);CHKERRQ(ierr)
      print *, vec_ptr(:)
      call VecRestoreArrayF90(lambda,vec_ptr,ierr);CHKERRQ(ierr)
    endif

    do i = 1, grid%nmax

      call MatZeroEntries(dMdK,ierr);CHKERRQ(ierr)
      call VecZeroEntries(dbdK,ierr);CHKERRQ(ierr)

      local_id = i
      do k = 1, inversion_aux%cell_to_internal_connection(0,local_id)
        iconn = inversion_aux%cell_to_internal_connection(k,local_id)
        local_id_up = grid%nG2L(connection_set%id_up(iconn))
        local_id_dn = grid%nG2L(connection_set%id_dn(iconn))
        flux_coef = inversion_aux%dFluxdIntConn(:,iconn)
        if (local_id_up == local_id) then
          call MatSetValue(dMdK,local_id_up-1,local_id_up-1, &
                            flux_coef(1), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK,local_id_up-1,local_id_dn-1, &
                            flux_coef(3), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK,local_id_dn-1,local_id_dn-1, &
                            -1.d0*flux_coef(3), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK,local_id_dn-1,local_id_up-1, &
                            -1.d0*flux_coef(1), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call VecSetValue(dbdK,local_id_up-1, &
                            -1.d0*flux_coef(5), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call VecSetValue(dbdK,local_id_dn-1, &
                            flux_coef(5), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
        elseif (local_id_dn == local_id) then
          call MatSetValue(dMdK,local_id_up-1,local_id_up-1, &
                            flux_coef(2), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK,local_id_up-1,local_id_dn-1, &
                            flux_coef(4), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK,local_id_dn-1,local_id_dn-1, &
                            -1.d0*flux_coef(4), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK,local_id_dn-1,local_id_up-1, &
                            -1.d0*flux_coef(2), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call VecSetValue(dbdK,local_id_up-1, &
                            -1.d0*flux_coef(6), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call VecSetValue(dbdK,local_id_dn-1, &
                            flux_coef(6), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
        else
          option%io_buffer = 'Incorrect mapping of connection'
          call PrintErrMsg(option)
        endif
      enddo
      flux_coef = 0.d0
      do k = 1, inversion_aux%cell_to_bc_connection(0,local_id)
        iconn = inversion_aux%cell_to_bc_connection(k,local_id)
        flux_coef(1:2) = inversion_aux%dFluxdBCConn(:,iconn)
        call MatSetValue(dMdK,local_id-1,local_id-1, &
                          flux_coef(1),ADD_VALUES,ierr);CHKERRQ(ierr)
        call VecSetValue(dbdK,local_id-1,flux_coef(2),ADD_VALUES, &
                          ierr);CHKERRQ(ierr)
      enddo
      call MatAssemblyBegin(dMdK,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      call MatAssemblyEnd(dMdK,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      call VecAssemblyBegin(dbdK,ierr);CHKERRQ(ierr)
      call VecAssemblyEnd(dbdK,ierr);CHKERRQ(ierr)

      if (this%debug_adjoint) then
        if (associated(inversion_aux%dMdK)) then
          print *, 'dMdK_petsc ', i
          call MatView(inversion_aux%dMdK(i),PETSC_VIEWER_STDOUT_WORLD, &
                       ierr);CHKERRQ(ierr)
        endif
        print *, 'dMdK ', i
        call MatView(dMdK,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
        if (associated(inversion_aux%dbdK)) then
          print *, 'dbdK_petsc ', i
          call VecView(inversion_aux%dbdK(i),PETSC_VIEWER_STDOUT_WORLD, &
                       ierr);CHKERRQ(ierr)
        endif
        print *, 'dbdK ', i
        call VecView(dbdK,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      endif
      if (associated(inversion_aux%dMdK)) then
        call MatAXPY(dMdK,-1.d0,inversion_aux%dMdK(i),SAME_NONZERO_PATTERN, &
                     ierr);CHKERRQ(ierr)
        call MatNorm(dMdK,NORM_FROBENIUS,tempreal,ierr);CHKERRQ(ierr)
        print *, 'dMdK diff norm: ', tempreal
      endif
      call MatMultTranspose(dMdK,lambda,work,ierr);CHKERRQ(ierr)
      if (this%debug_adjoint) then
        print *, 'dMdK^T * lambda'
        call VecView(work,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      endif
      call VecDot(solution,work,hTdMdKTlambda,ierr);CHKERRQ(ierr)
      if (this%debug_adjoint) then
        print *, '-h^T * dMdK^T * lambda'
        print *, i, ' : ', -hTdMdKTlambda
      endif
      call VecDot(dbdK,lambda,dbdKTlambda,ierr);CHKERRQ(ierr)
      if (this%debug_adjoint) then
        print *, 'dbdK^T * lambda'
        print *, i, ' : ', dbdKTlambda
      endif
      tempreal = dbdKTlambda-hTdMdKTlambda
      if (this%debug_adjoint) then
        print *, '(dbdK^T - h^T * dMdK^T) * lambda'
        print *, i, ' : ', tempreal
      endif
      call MatSetValue(inversion_aux%Jsensitivity,j-1,i-1,-tempreal, &
                        INSERT_VALUES,ierr);CHKERRQ(ierr)
    enddo
  enddo
  call MatDestroy(dMdK,ierr);CHKERRQ(ierr)
  call VecDestroy(dbdK,ierr);CHKERRQ(ierr)
  call VecDestroy(p,ierr);CHKERRQ(ierr)
  call VecDestroy(lambda,ierr);CHKERRQ(ierr)
  call MatAssemblyBegin(inversion_aux%Jsensitivity, &
                        MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(inversion_aux%Jsensitivity, &
                      MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  deallocate(flux_coef)

end subroutine InvTestCalculateSensitivity

! ************************************************************************** !

subroutine InversionTestDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  class(inversion_test_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionTestDestroy

end module Inversion_Test_class
