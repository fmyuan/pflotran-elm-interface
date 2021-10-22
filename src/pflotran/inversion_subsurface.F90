module Inversion_Subsurface_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Aux_module
  use Inversion_Base_class
  use Realization_Subsurface_class
  use Simulation_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_base_type) :: inversion_subsurface_type
    character(len=MAXSTRINGLENGTH) :: forward_simulation_filename
    class(simulation_subsurface_type), pointer :: forward_simulation
    class(realization_subsurface_type), pointer :: realization
    type(inversion_aux_type), pointer :: inversion_aux
    PetscReal, pointer :: measurement(:)
    PetscInt, pointer :: imeasurement(:)
    Vec :: quantity_of_interest
    PetscInt :: iqoi
    Vec :: ref_quantity_of_interest
    character(len=MAXWORDLENGTH) :: ref_qoi_dataset_name
    PetscBool :: print_sensitivity_jacobian
    PetscBool :: debug_adjoint
  contains
    procedure, public :: Init => InversionSubsurfaceInit
    procedure, public :: Initialize => InversionSubsurfInitialize
    procedure, public :: Step => InversionSubsurfaceStep
    procedure, public :: ConnectToForwardRun => InvSubsurfConnectToForwardRun
    procedure, public :: CalculateSensitivity => InvSubsurfCalculateSensitivity
    procedure, public :: OutputSensitivity => InvSubsurfOutputSensitivity
    procedure, public :: Invert => InversionSubsurfaceInvert
  end type inversion_subsurface_type

  public :: InversionSubsurfaceInit, &
            InversionSubsurfReadSelectCase, &
            InversionSubsurfInitialize, &
            InvSubsurfConnectToForwardRun, &
            InvSubsurfOutputSensitivity, &
            InversionSubsurfaceStrip

contains

! ************************************************************************** !

subroutine InversionSubsurfaceInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  use Driver_module
  use ZFlow_Aux_module, only : zflow_calc_adjoint

  class(inversion_subsurface_type) :: this
  class(driver_type), pointer :: driver

  call InversionBaseInit(this,driver)

  this%quantity_of_interest = PETSC_NULL_VEC
  this%iqoi = UNINITIALIZED_INTEGER
  this%ref_quantity_of_interest = PETSC_NULL_VEC
  this%ref_qoi_dataset_name = ''
  this%forward_simulation_filename = ''
  this%print_sensitivity_jacobian = PETSC_FALSE
  this%debug_adjoint = PETSC_FALSE

  nullify(this%measurement)
  nullify(this%imeasurement)

  nullify(this%forward_simulation)
  nullify(this%realization)
  nullify(this%inversion_aux)

  zflow_calc_adjoint = PETSC_TRUE

end subroutine InversionSubsurfaceInit

! ************************************************************************** !

subroutine InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                          error_string,option)

  use Input_Aux_module
  use Option_module
  use String_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, POROSITY
  use Utility_module

  class(inversion_subsurface_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, pointer :: tempint(:)
  PetscReal, pointer :: tempreal(:)

  found = PETSC_TRUE
  call InversionBaseReadSelectCase(this,input,keyword,found, &
                                   error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('FORWARD_SIMULATION_FILENAME')
      call InputReadFilename(input,option,this%forward_simulation_filename)
      call InputErrorMsg(input,option,keyword,error_string)
    case('QUANTITY_OF_INTEREST')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,keyword,error_string)
      call StringToUpper(word)
      select case(word)
        case('ELECTRICAL_CONDUCTIVITY')
          this%iqoi = ELECTRICAL_CONDUCTIVITY
        case('PERMEABILITY')
          this%iqoi = PERMEABILITY
        case('POROSITY')
          this%iqoi = POROSITY
        case default
          call InputKeywordUnrecognized(input,word,trim(error_string)// &
                                        & ',QUANTITY_OF_INTEREST',option)
      end select
    case('REFERENCE_QUANTITY_OF_INTEREST')
      call InputReadNChars(input,option,this%ref_qoi_dataset_name, &
                           MAXWORDLENGTH,PETSC_TRUE)
      call InputErrorMsg(input,option,'DATASET NAME', &
                         keyword)
    case('MEASUREMENTS')
      string = trim(error_string)//keyword
      i = 10
      allocate(tempint(i))
      tempint = UNINITIALIZED_INTEGER
      allocate(tempreal(i))
      tempreal = UNINITIALIZED_DOUBLE
      i = 0
      do
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,error_string)
        if (InputCheckExit(input,option)) exit
        i = i + 1
        if (i > size(tempint)) then
          call ReallocateArray(tempint)
          call ReallocateArray(tempreal)
        endif
        call InputReadInt(input,option,tempint(i))
        call InputErrorMsg(input,option,'cell id',string)
        call InputReadDouble(input,option,tempreal(i))
        call InputErrorMsg(input,option,'measurement',string)
      enddo
      allocate(this%imeasurement(i))
      this%imeasurement(:) = tempint(1:i)
      allocate(this%measurement(i))
      this%measurement(:) = tempreal(1:i)
      call DeallocateArray(tempint)
      call DeallocateArray(tempreal)
    case('PRINT_SENSITIVITY_JACOBIAN')
      this%print_sensitivity_jacobian = PETSC_TRUE
    case('DEBUG_ADJOINT')
      this%debug_adjoint = PETSC_TRUE
    case default
      found = PETSC_FALSE
  end select

end subroutine InversionSubsurfReadSelectCase

! ************************************************************************** !

subroutine InversionSubsurfInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Connection_module
  use Coupler_module
  use Grid_module
  use Patch_module

  class(inversion_subsurface_type) :: this

  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(inversion_aux_type), pointer :: inversion_aux
  type(patch_type), pointer :: patch
  PetscInt :: local_id
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscInt, allocatable :: int_array(:)

  !TODO(geh): compress the mappings

  if (.not.associated(this%inversion_aux)) then
    patch => this%realization%patch
    inversion_aux => InversionAuxCreate()
    call GridMapCellsToConnections(patch%grid, &
                                   inversion_aux%cell_to_internal_connection)
    sum_connection = &
      ConnectionGetNumberInList(patch%grid%internal_connection_set_list)
    allocate(inversion_aux%dFluxdIntConn(6,sum_connection))
    inversion_aux%dFluxdIntConn = 0.d0
    sum_connection = &
      CouplerGetNumConnectionsInList(patch%boundary_condition_list)
    allocate(inversion_aux%dFluxdBCConn(2,sum_connection))
    inversion_aux%dFluxdBCConn = 0.d0

    allocate(int_array(patch%grid%nlmax))
    int_array = 0
    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0
    do
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        int_array(local_id) = int_array(local_id) + 1
      enddo
      boundary_condition => boundary_condition%next
    enddo
    iconn = maxval(int_array)
    deallocate(int_array)
    allocate(inversion_aux% &
               cell_to_bc_connection(0:iconn,patch%grid%nlmax))
    inversion_aux%cell_to_bc_connection = 0
    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0
    do
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        inversion_aux%cell_to_bc_connection(0,local_id) = &
          inversion_aux%cell_to_bc_connection(0,local_id) + 1
          inversion_aux%cell_to_bc_connection( &
            inversion_aux%cell_to_bc_connection(0,local_id), &
            local_id) = sum_connection
      enddo
      boundary_condition => boundary_condition%next
    enddo
    this%inversion_aux => inversion_aux
  endif

end subroutine InversionSubsurfInitialize

! ************************************************************************** !

subroutine InversionSubsurfaceStep(this)
  !
  ! Execute a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Option_module
  use Factory_Forward_module

  class(inversion_subsurface_type) :: this

  type(option_type), pointer :: option

  option => OptionCreate()
  write(option%group_prefix,'(i6)') this%iteration
  option%group_prefix = 'Run' // trim(adjustl(option%group_prefix))
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
  call this%CalculateSensitivity()
  call this%OutputSensitivity('adjoint')
  nullify(this%realization)
  call this%forward_simulation%FinalizeRun()
  call this%forward_simulation%Strip()
  deallocate(this%forward_simulation)
  nullify(this%forward_simulation)
  call this%Invert()

  this%converg_flag = PETSC_FALSE
  if (this%iteration > this%maximum_iteration) this%converg_flag = PETSC_TRUE

end subroutine InversionSubsurfaceStep

! ************************************************************************** !

subroutine InvSubsurfConnectToForwardRun(this)
  !
  ! Sets up the interface between a single forward run and the outer
  ! inversion wrapper
  !
  ! Author: Glenn Hammond
  ! Date: 10/18/21
  !
  use Discretization_module
  use Material_module

  class(inversion_subsurface_type) :: this

  PetscInt :: num_measurements
  PetscErrorCode :: ierr

  this%realization%patch%aux%inversion_aux => this%inversion_aux

  ! on first pass, store and set thereafter
  if (this%quantity_of_interest == PETSC_NULL_VEC) then
    num_measurements = size(this%imeasurement)
    call VecDuplicate(this%realization%field%work, &
                      this%quantity_of_interest,ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(this%realization%discretization, &
                                     this%realization%field%work_loc, &
                                     this%quantity_of_interest,ONEDOF)
    call MatCreateDense(this%driver%comm%mycomm, &
                        num_measurements, &
                        this%realization%patch%grid%nlmax, &
                        num_measurements, &
                        this%realization%patch%grid%nmax, &
                        PETSC_NULL_SCALAR, &
                        this%inversion_aux%Jsensitivity,ierr);CHKERRQ(ierr)
  endif

  call DiscretizationGlobalToLocal(this%realization%discretization, &
                                   this%quantity_of_interest, &
                                   this%realization%field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                               this%realization%field%work_loc, &
                               this%iqoi,ZERO_INTEGER)

end subroutine InvSubsurfConnectToForwardRun

! ************************************************************************** !

subroutine InvSubsurfCalculateSensitivity(this)
  !
  ! Calculates sensitivity matrix Jsensitivity
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Connection_module
  use Debug_module
  use Grid_module
  use Option_module
  use Patch_module
  use Solver_module
  use String_module

  use PM_Base_class
  use PM_ZFlow_class

  class(inversion_subsurface_type) :: this

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
  Mat :: dMdK_
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

  call MatDuplicate(solver%M,MAT_SHARE_NONZERO_PATTERN,dMdK_, &
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

      call MatZeroEntries(dMdK_,ierr);CHKERRQ(ierr)
      call VecZeroEntries(dbdK,ierr);CHKERRQ(ierr)

      local_id = i
      do k = 1, inversion_aux%cell_to_internal_connection(0,local_id)
        iconn = inversion_aux%cell_to_internal_connection(k,local_id)
        local_id_up = grid%nG2L(connection_set%id_up(iconn))
        local_id_dn = grid%nG2L(connection_set%id_dn(iconn))
        flux_coef = inversion_aux%dFluxdIntConn(:,iconn)
        if (local_id_up == local_id) then
          call MatSetValue(dMdK_,local_id_up-1,local_id_up-1, &
                            flux_coef(1), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK_,local_id_up-1,local_id_dn-1, &
                            flux_coef(3), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK_,local_id_dn-1,local_id_dn-1, &
                            -1.d0*flux_coef(3), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK_,local_id_dn-1,local_id_up-1, &
                            -1.d0*flux_coef(1), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call VecSetValue(dbdK,local_id_up-1, &
                            -1.d0*flux_coef(5), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call VecSetValue(dbdK,local_id_dn-1, &
                            flux_coef(5), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
        elseif (local_id_dn == local_id) then
          call MatSetValue(dMdK_,local_id_up-1,local_id_up-1, &
                            flux_coef(2), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK_,local_id_up-1,local_id_dn-1, &
                            flux_coef(4), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK_,local_id_dn-1,local_id_dn-1, &
                            -1.d0*flux_coef(4), &
                            ADD_VALUES,ierr);CHKERRQ(ierr)
          call MatSetValue(dMdK_,local_id_dn-1,local_id_up-1, &
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
        call MatSetValue(dMdK_,local_id-1,local_id-1, &
                          flux_coef(1),ADD_VALUES,ierr);CHKERRQ(ierr)
        call VecSetValue(dbdK,local_id-1,flux_coef(2),ADD_VALUES, &
                          ierr);CHKERRQ(ierr)
      enddo
      call MatAssemblyBegin(dMdK_,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      call MatAssemblyEnd(dMdK_,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      call VecAssemblyBegin(dbdK,ierr);CHKERRQ(ierr)
      call VecAssemblyEnd(dbdK,ierr);CHKERRQ(ierr)

      if (this%debug_adjoint) then
        print *, 'dMdK ', i
        call MatView(dMdK_,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
        print *, 'dbdK ', i
        call VecView(dbdK,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      endif
      call MatMultTranspose(dMdK_,lambda,work,ierr);CHKERRQ(ierr)
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
  call MatDestroy(dMdK_,ierr);CHKERRQ(ierr)
  call VecDestroy(dbdK,ierr);CHKERRQ(ierr)
  call VecDestroy(p,ierr);CHKERRQ(ierr)
  call VecDestroy(lambda,ierr);CHKERRQ(ierr)
  call MatAssemblyBegin(inversion_aux%Jsensitivity, &
                        MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(inversion_aux%Jsensitivity, &
                      MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  deallocate(flux_coef)

end subroutine InvSubsurfCalculateSensitivity

! ************************************************************************** !

subroutine InvSubsurfScaleSensitivity(this,Jsensitivity)
  !
  ! Writes sensitivity Jacobian to an ASCII output file
  !
  ! Author: Glenn hammond
  ! Date: 10/11/21
  !
  use Realization_Base_class
  use Variables_module, only : PERMEABILITY

  class(inversion_subsurface_type) :: this
  Mat :: Jsensitivity

  PetscErrorCode :: ierr

  call RealizationGetVariable(this%realization, &
                              this%realization%field%work, &
                              PERMEABILITY,ZERO_INTEGER)
  call MatDiagonalScale(Jsensitivity,PETSC_NULL_VEC, &
                        this%realization%field%work, &
                        ierr);CHKERRQ(ierr)

end subroutine InvSubsurfScaleSensitivity

! ************************************************************************** !

subroutine InvSubsurfOutputSensitivity(this,suffix)
  !
  ! Writes sensitivity Jacobian
  !
  ! Author: Glenn hammond
  ! Date: 10/11/21
  !

  class(inversion_subsurface_type) :: this
  character(len=*) :: suffix

  character(len=MAXSTRINGLENGTH) :: filename_prefix

  filename_prefix = 'Jsensitivity'
  if (len_trim(suffix) > 0) filename_prefix = trim(filename_prefix) // '_' // &
                            suffix
  call InvSubsurfOutputSensitivityASCII(this,this%inversion_aux%Jsensitivity, &
                                        filename_prefix)
  call InvSubsurfOutputSensitivityHDF5(this,this%inversion_aux%Jsensitivity, &
                                       filename_prefix)

end subroutine InvSubsurfOutputSensitivity

! ************************************************************************** !

subroutine InvSubsurfOutputSensitivityASCII(this,Jsensitivity,filename_prefix)
  !
  ! Writes sensitivity Jacobian to an ASCII output file
  !
  ! Author: Glenn hammond
  ! Date: 10/11/21
  !
  use Realization_Subsurface_class
  use Variables_module, only : PERMEABILITY

  class(inversion_subsurface_type) :: this
  Mat :: Jsensitivity
  character(len=*) :: filename_prefix

  character(len=MAXSTRINGLENGTH) :: string
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  if (.not.associated(this%realization)) then
    call this%driver%PrintErrMsg('InvSubsurfOutputSensitivityASCII must be &
           &called before the forward simulation is destroyed.')
  endif

  string = trim(filename_prefix) // '.txt'
  call PetscViewerASCIIOpen(this%realization%option%mycomm,string,viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(Jsensitivity,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

end subroutine InvSubsurfOutputSensitivityASCII

! ************************************************************************** !

subroutine InvSubsurfOutputSensitivityHDF5(this,Jsensitivity,filename_prefix)
  !
  ! Writes sensitivity Jacobian to an HDF5 output file
  !
  ! Author: Glenn hammond
  ! Date: 10/18/21
  !
  use hdf5
  use HDF5_module
  use Output_HDF5_module
  use String_module

  class(inversion_subsurface_type) :: this
  Mat :: Jsensitivity
  character(len=*) :: filename_prefix

  Vec :: row_vec
  PetscReal, pointer :: row_ptr(:)
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscInt :: imeasurement
  PetscInt :: num_measurement
  PetscErrorCode :: ierr

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
  character(len=MAXSTRINGLENGTH) :: string
  integer :: hdf5_err

  if (.not.associated(this%realization)) then
    call this%driver%PrintErrMsg('InvSubsurfOutputSensitivityHDF5 must be &
           &called before the forward simulation is destroyed.')
  endif

  !HDF5 formatted output
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,this%realization%option%mycomm, &
                          MPI_INFO_NULL,hdf5_err)
#endif
  string = trim(filename_prefix) // '.h5'
  call h5fcreate_f(string,H5F_ACC_TRUNC_F,file_id,hdf5_err, &
                   H5P_DEFAULT_F,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  call OutputHDF5WriteStructCoordGroup(file_id, &
                                       this%realization%discretization, &
                                       this%realization%patch%grid, &
                                       this%realization%option)
  ! create a group for the data set
  this%realization%option%time = 0.d0
  write(string,'(''Time:'',es13.5,x,a1)') &
        this%realization%option%time/this%realization%output_option%tconv, &
        this%realization%output_option%tunit
  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)

  num_measurement = size(this%measurement)
  call VecCreateMPI(this%realization%option%mycomm,num_measurement, &
                    num_measurement, &
                    row_vec,ierr);CHKERRQ(ierr)
  do imeasurement = 1, num_measurement
    call VecZeroEntries(row_vec,ierr);CHKERRQ(ierr)
    call VecSetValue(row_vec,imeasurement-1,1.d0,INSERT_VALUES, &
                     ierr);CHKERRQ(ierr)
    call VecAssemblyBegin(row_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(row_vec,ierr);CHKERRQ(ierr)
    call MatMultTranspose(Jsensitivity,row_vec, &
                          this%realization%field%work, &
                          ierr);CHKERRQ(ierr)
    string = 'Measurement ' // StringWrite(this%imeasurement(imeasurement))
    call HDF5WriteStructDataSetFromVec(string,this%realization, &
                                       this%realization%field%work,grp_id, &
                                       H5T_NATIVE_DOUBLE)
  enddo
  call VecDestroy(row_vec,ierr);CHKERRQ(ierr)
  call h5gclose_f(grp_id,hdf5_err)
  call OutputHDF5CloseFile(this%realization%option,file_id)

end subroutine InvSubsurfOutputSensitivityHDF5

! ************************************************************************** !

subroutine InversionSubsurfaceInvert(this)
  !
  ! Inverts for a parameter update
  !
  ! Author: Glenn hammond
  ! Date: 09/20/21
  !

  class(inversion_subsurface_type) :: this

end subroutine InversionSubsurfaceInvert

! ************************************************************************** !

subroutine InversionSubsurfaceStrip(this)
  !
  ! Deallocates members of inversion Subsurface
  !
  ! Author: Glenn hammond
  ! Date: 09/20/21
  !
  use Utility_module

  class(inversion_subsurface_type) :: this

  PetscErrorCode :: ierr

  call InversionBaseStrip(this)

  nullify(this%realization)
  call InversionAuxDestroy(this%inversion_aux)
  if (associated(this%forward_simulation)) then
    print *, 'Why is forward simulation still associated in &
             &InversionSubSurfStrip?'
    stop
  endif
  nullify(this%forward_simulation)
  if (this%quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%quantity_of_interest,ierr);CHKERRQ(ierr)
  endif
  if (this%ref_quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%ref_quantity_of_interest,ierr);CHKERRQ(ierr)
  endif

  call DeallocateArray(this%imeasurement)
  call DeallocateArray(this%measurement)

end subroutine InversionSubsurfaceStrip

end module Inversion_Subsurface_class
