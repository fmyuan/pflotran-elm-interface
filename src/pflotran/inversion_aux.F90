module Inversion_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module
  use Characteristic_Curves_module
  use Communicator_Aux_module
  use Driver_class
  use Inversion_Coupled_Aux_module
  use Inversion_Measurement_Aux_module
  use Inversion_Parameter_module
  use Inversion_TS_Aux_module
  use Material_module
  use Option_Inversion_module

  implicit none

  private

  PetscInt, parameter, public :: INVAUX_GET_MATERIAL_VALUE = 0
  PetscInt, parameter, public :: INVAUX_OVERWRITE_MATERIAL_VALUE = 1
  PetscInt, parameter, public :: INVAUX_COPY_TO_VEC = 2
  PetscInt, parameter, public :: INVAUX_COPY_FROM_VEC = 3
  PetscInt, parameter, public :: INVAUX_PARAMETER_VALUE = 4
  PetscInt, parameter, public :: INVAUX_PARAMETER_UPDATE = 5

  PetscInt, parameter, public :: INVAUX_SCATFORWARD = 0
  PetscInt, parameter, public :: INVAUX_SCATREVERSE = 1

  type, public :: inversion_aux_type
    class(driver_type), pointer :: driver
    type(material_property_ptr_type), pointer :: material_property_array(:)
    type(characteristic_curves_ptr_type), pointer :: cc_array(:)
    PetscBool :: qoi_is_full_vector
    PetscBool :: startup_phase
    Vec :: solution ! solely a pointer
    PetscInt :: isync_time              ! current index of sync_times
    PetscReal, pointer :: sync_times(:) ! an array with all measurement times
    type(inversion_coupled_aux_type), pointer :: coupled_aux
    type(inversion_measurement_aux_type), pointer :: measurements(:)
    type(inversion_parameter_type), pointer :: parameters(:)
    Vec :: measurement_vec
    Vec :: dist_measurement_vec
    Vec :: parameter_vec
    Vec :: dist_parameter_vec
    Vec :: del_parameter_vec
    VecScatter :: scatter_measure_to_dist_measure
    VecScatter :: scatter_param_to_dist_param
    VecScatter :: scatter_global_to_dist_param
    Mat :: JsensitivityT
    ! adjoint data structures
    PetscBool :: store_adjoint
    Mat :: M_ptr
    type(inversion_forward_ts_aux_type), pointer :: first_forward_ts_aux
    type(inversion_forward_ts_aux_type), pointer :: last_forward_ts_aux
    PetscReal, pointer :: local_measurement_values_ptr(:)
    PetscReal, pointer :: local_dobs_dunknown_values_ptr(:)
    PetscReal, pointer :: local_dobs_dparam_values_ptr(:)
    type(inversion_perturbation_type), pointer :: perturbation
  end type inversion_aux_type

  type, public :: inversion_perturbation_type
    Vec :: base_parameter_vec
    Vec :: base_measurement_vec
    PetscInt :: ndof
    PetscInt :: idof_pert
    PetscReal :: pert
    PetscReal :: base_value
    PetscReal :: tolerance
    PetscInt, pointer :: select_cells(:)
  end type inversion_perturbation_type

  public :: InversionAuxCreate, &
            InversionAuxPerturbationCreate, &
            InversionAuxResetMeasurements, &
            InversionAuxAdjointRecordTS, &
            InvAuxAdjCleanupAfterForwardRun, &
            InversionAuxDestroy

  public :: InvAuxCopyParameterValue, &
            InvAuxCopyParamToFromParamVec, &
            InvAuxGetParamValueByCell, &
            InvAuxGetSetParamValueByMat, &
            InvAuxScatMeasToDistMeas, &
            InvAuxCopyMeasToFromMeasVec, &
            InvAuxScatParamToDistParam, &
            InvAuxScatGlobalToDistParam, &
            InvAuxBCastVecForCommI, &
            InvAuxParamVecToMaterial, &
            InvAuxMaterialToParamVec

contains

! ************************************************************************** !

function InversionAuxCreate(driver)
  !
  ! Allocate and initialize auxiliary inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  class(driver_type), pointer :: driver

  type(inversion_aux_type), pointer :: InversionAuxCreate

  type(inversion_aux_type), pointer :: aux

  allocate(aux)

  aux%solution = PETSC_NULL_VEC

  aux%driver => driver
  nullify(aux%material_property_array)
  nullify(aux%cc_array)
  aux%qoi_is_full_vector = PETSC_FALSE
  aux%startup_phase = PETSC_TRUE
  aux%isync_time = 1
  nullify(aux%sync_times)
  nullify(aux%coupled_aux)
  nullify(aux%measurements)
  nullify(aux%parameters)
  aux%measurement_vec = PETSC_NULL_VEC
  aux%dist_measurement_vec = PETSC_NULL_VEC
  aux%parameter_vec = PETSC_NULL_VEC
  aux%dist_parameter_vec = PETSC_NULL_VEC
  aux%del_parameter_vec = PETSC_NULL_VEC
  aux%scatter_measure_to_dist_measure = PETSC_NULL_VECSCATTER
  aux%scatter_param_to_dist_param = PETSC_NULL_VECSCATTER
  aux%scatter_global_to_dist_param = PETSC_NULL_VECSCATTER
  aux%JsensitivityT = PETSC_NULL_MAT
  nullify(aux%local_measurement_values_ptr)
  nullify(aux%local_dobs_dunknown_values_ptr)
  nullify(aux%local_dobs_dparam_values_ptr)
  ! adjoint
  call InversionAuxInitAdjoint(aux)
  ! perturbation
  nullify(aux%perturbation)

  InversionAuxCreate => aux

end function InversionAuxCreate

! ************************************************************************** !

subroutine InversionAuxInitAdjoint(aux)
  !
  ! Initializes adjoint portion of object
  !
  ! Author: Glenn Hammond
  ! Date: 11/28/22

  type(inversion_aux_type) :: aux

  aux%store_adjoint = PETSC_TRUE
  aux%M_ptr = PETSC_NULL_MAT
  nullify(aux%first_forward_ts_aux)
  nullify(aux%last_forward_ts_aux)

end subroutine InversionAuxInitAdjoint

! ************************************************************************** !

function InversionAuxPerturbationCreate()
  !
  ! Allocates and initializes a new perturbation object
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  type(inversion_perturbation_type), pointer :: InversionAuxPerturbationCreate

  allocate(InversionAuxPerturbationCreate)
  InversionAuxPerturbationCreate%base_parameter_vec = PETSC_NULL_VEC
  InversionAuxPerturbationCreate%base_measurement_vec = PETSC_NULL_VEC

  InversionAuxPerturbationCreate%ndof = 0
  InversionAuxPerturbationCreate%idof_pert = 0
  InversionAuxPerturbationCreate%pert = 0.d0
  InversionAuxPerturbationCreate%base_value = 0.d0
  InversionAuxPerturbationCreate%tolerance = 1.d-6
  nullify(InversionAuxPerturbationCreate%select_cells)

end function InversionAuxPerturbationCreate

! ************************************************************************** !

subroutine InversionAuxResetMeasurements(aux)
  !
  ! Resets flags for forward run back to original settings.
  !
  ! Author: Glenn Hammond
  ! Date: 02/21/22

  type(inversion_aux_type), pointer :: aux

  PetscInt :: imeasurement

  aux%isync_time = 1
  do imeasurement = 1, size(aux%measurements)
    call InversionMeasurementAuxReset(aux%measurements(imeasurement))
  enddo

end subroutine InversionAuxResetMeasurements

! ************************************************************************** !

subroutine InversionAuxAdjointRecordTS(aux,time)
  !
  ! Appends a time step to the linked list
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/22, 11/28/22

  use Utility_module

  implicit none

  type(inversion_aux_type), pointer :: aux
  PetscReal :: time

  if (associated(aux%last_forward_ts_aux)) then
    aux%last_forward_ts_aux%time = time
    ! store the solution
    call InvForTSAuxDupForwardJacobian(aux%M_ptr,aux%last_forward_ts_aux)
    ! append next time step
    aux%last_forward_ts_aux => &
      InversionForwardTSAuxCreate(aux%last_forward_ts_aux)
  endif

end subroutine InversionAuxAdjointRecordTS

! ************************************************************************** !

subroutine InvAuxAdjCleanupAfterForwardRun(aux)
  !
  ! Destroys the linked list of adjoint objects
  !
  ! Author: Glenn Hammond
  ! Date: 11/28/22

  implicit none

  type(inversion_aux_type) :: aux

  ! destroy the lists
  call InvForwardTSAuxDestroyList(aux%first_forward_ts_aux,PETSC_FALSE)
  ! initialize everything else
  call InversionAuxInitAdjoint(aux)

end subroutine InvAuxAdjCleanupAfterForwardRun

! ************************************************************************** !

subroutine InvAuxCopyParameterValue(aux,iparam,iflag)
  !
  ! Copies parameter values back and forth
  !
  ! Author: Glenn Hammond
  ! Date: 03/30/22

  type(inversion_aux_type) :: aux
  PetscInt :: iparam
  PetscInt :: iflag

  PetscReal :: tempreal

  if (iflag /= INVAUX_GET_MATERIAL_VALUE) then
    ! everything else is implicit OVERWRITE_MATERIAL_VALUE
    tempreal = aux%parameters(iparam)%value
  endif

  call InvAuxGetSetParamValueByMat(aux,tempreal, &
                                   aux%parameters(iparam)%itype, &
                                   aux%parameters(iparam)%imat,iflag)

  if (iflag == INVAUX_GET_MATERIAL_VALUE) then
    aux%parameters(iparam)%value = tempreal
  endif

end subroutine InvAuxCopyParameterValue

! ************************************************************************** !

subroutine InvAuxGetSetParamValueByMat(aux,value,iparameter_type,imat,iflag)
  !
  ! Copies parameter values back and forth
  !
  ! Author: Glenn Hammond
  ! Date: 03/30/22

  use String_module
  use Utility_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, PERMEABILITY, &
                               POROSITY, VG_ALPHA, VG_SR, VG_M, &
                               ARCHIE_CEMENTATION_EXPONENT, &
                               ARCHIE_SATURATION_EXPONENT, &
                               ARCHIE_TORTUOSITY_CONSTANT

  type(inversion_aux_type) :: aux
  PetscReal :: value
  PetscInt :: iparameter_type
  PetscInt :: imat
  PetscInt :: iflag

  type(material_property_type), pointer :: material_property
  class(characteristic_curves_type), pointer :: cc
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: tempreal

  material_property => aux%material_property_array(imat)%ptr
  select case(iparameter_type)
    case(ELECTRICAL_CONDUCTIVITY)
      if (iflag == INVAUX_GET_MATERIAL_VALUE) then
        value = material_property%electrical_conductivity
      else
        material_property%electrical_conductivity = value
      endif
    case(PERMEABILITY)
      if (iflag == INVAUX_GET_MATERIAL_VALUE) then
        value = material_property%permeability(1,1)
      else
        material_property%permeability(1,1) = value
        material_property%permeability(2,2) = value
        if (Initialized(material_property%vertical_anisotropy_ratio)) then
          value = value * material_property%vertical_anisotropy_ratio
        endif
        material_property%permeability(3,3) = value
      endif
    case(POROSITY)
      if (iflag == INVAUX_GET_MATERIAL_VALUE) then
        value = material_property%porosity
      else
        material_property%porosity = value
      endif
    case(VG_ALPHA,VG_SR,VG_M)
      cc => aux%cc_array(material_property%saturation_function_id)%ptr
      select case(iparameter_type)
        case(VG_ALPHA)
          if (iflag == INVAUX_GET_MATERIAL_VALUE) then
            value = cc%saturation_function%GetAlpha_()
          else
            call cc%saturation_function%SetAlpha_(value)
          endif
        case(VG_M)
          if (iflag == INVAUX_GET_MATERIAL_VALUE) then
            value = cc%saturation_function%GetM_()
            tempreal = cc%liq_rel_perm_function%GetM_()
            if (.not.Equal(value,tempreal)) then
              string = 'For inversion, saturation and relative permeability &
                &function van Genuchten "m" values must match in &
                &characteristic curve "' // trim(cc%name)
              call aux%driver%PrintErrMsg(string)
            endif
          else
            call cc%saturation_function%SetM_(value)
            call cc%liq_rel_perm_function%SetM_(value)
          endif
        case(VG_SR)
          if (iflag == INVAUX_GET_MATERIAL_VALUE) then
            value = cc%saturation_function%GetResidualSaturation()
            tempreal = cc%liq_rel_perm_function%GetResidualSaturation()
            if (.not.Equal(value,tempreal)) then
              string = 'For inversion, saturation and relative permeability &
                &function  saturations must match in characteristic &
                &curve "' // trim(cc%name)
              call aux%driver%PrintErrMsg(string)
            endif
          else
            call cc%saturation_function%SetResidualSaturation(value)
            call cc%liq_rel_perm_function%SetResidualSaturation(value)
          endif
      end select
    case(ARCHIE_CEMENTATION_EXPONENT)
      if (iflag == INVAUX_GET_MATERIAL_VALUE) then
        value = material_property%archie_cementation_exponent
      else
        material_property%archie_cementation_exponent = value
      endif
    case(ARCHIE_SATURATION_EXPONENT)
      if (iflag == INVAUX_GET_MATERIAL_VALUE) then
        value = material_property%archie_saturation_exponent
      else
        material_property%archie_saturation_exponent = value
      endif
    case(ARCHIE_TORTUOSITY_CONSTANT)
      if (iflag == INVAUX_GET_MATERIAL_VALUE) then
        value = material_property%archie_tortuosity_constant
      else
        material_property%archie_tortuosity_constant = value
      endif
    case default
      string = 'Unrecognized variable in &
        &InvAuxGetSetParamValueByMat: ' // &
        trim(StringWrite(iparameter_type))
      call aux%driver%PrintErrMsg(string)
  end select

end subroutine InvAuxGetSetParamValueByMat

! ************************************************************************** !

subroutine InvAuxCopyParamToFromParamVec(aux,itype,idirection)
  !
  ! Copies parameter values back and forth
  !
  ! Author: Glenn Hammond
  ! Date: 03/30/22

  class(inversion_aux_type) :: aux
  PetscInt :: itype
  PetscInt :: idirection

  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscErrorCode :: ierr

  select case(itype)
    case(INVAUX_PARAMETER_VALUE)
      call VecGetArrayF90(aux%parameter_vec,vec_ptr,ierr);CHKERRQ(ierr)
      select case(idirection)
        case(INVAUX_COPY_FROM_VEC)
          do i = 1, size(aux%parameters)
            aux%parameters(i)%value = vec_ptr(i)
          enddo
        case(INVAUX_COPY_TO_VEC)
          do i = 1, size(aux%parameters)
            vec_ptr(i) = aux%parameters(i)%value
          enddo
        case default
          stop 'Error: idirection in InvAuxCopyParamToFromParamVec,Value'
      end select
      call VecRestoreArrayF90(aux%parameter_vec,vec_ptr,ierr);CHKERRQ(ierr)
    case(INVAUX_PARAMETER_UPDATE)
      call VecGetArrayF90(aux%del_parameter_vec,vec_ptr,ierr);CHKERRQ(ierr)
      select case(idirection)
        case(INVAUX_COPY_FROM_VEC)
          do i = 1, size(aux%parameters)
            aux%parameters(i)%update = vec_ptr(i)
          enddo
        case(INVAUX_COPY_TO_VEC)
          do i = 1, size(aux%parameters)
            vec_ptr(i) = aux%parameters(i)%update
          enddo
        case default
          stop 'Error: idirection in InvAuxCopyParamToFromParamVec,Update'
      end select
      call VecRestoreArrayF90(aux%del_parameter_vec,vec_ptr,ierr);CHKERRQ(ierr)
    case default
      stop 'Error: itype in InvAuxCopyParamToFromParamVec'
  end select

end subroutine InvAuxCopyParamToFromParamVec

! ************************************************************************** !

subroutine InvAuxParamVecToMaterial(aux)
  !
  ! Copies parameter values from parameter vec to material properties
  !
  ! Author: Glenn Hammond
  ! Date: 03/10/23

  class(inversion_aux_type) :: aux

  PetscInt :: i

  call InvAuxCopyParamToFromParamVec(aux,INVAUX_PARAMETER_VALUE, &
                                     INVAUX_COPY_FROM_VEC)
  do i = 1, size(aux%parameters)
    call InvAuxCopyParameterValue(aux,i,INVAUX_OVERWRITE_MATERIAL_VALUE)
  enddo

end subroutine InvAuxParamVecToMaterial

! ************************************************************************** !

subroutine InvAuxMaterialToParamVec(aux)
  !
  ! Copies parameter values from material properties to parameter vec
  !
  ! Author: Glenn Hammond
  ! Date: 03/10/23

  class(inversion_aux_type) :: aux

  PetscInt :: i

  do i = 1, size(aux%parameters)
    call InvAuxCopyParameterValue(aux,i,INVAUX_GET_MATERIAL_VALUE)
  enddo
  call InvAuxCopyParamToFromParamVec(aux,INVAUX_PARAMETER_VALUE, &
                                     INVAUX_COPY_TO_VEC)

end subroutine InvAuxMaterialToParamVec

! ************************************************************************** !

subroutine InvAuxGetParamValueByCell(aux,value,iparameter_type,imat, &
                                     material_auxvar)
  !
  ! Returns the parameter value at the cell
  !
  ! Author: Glenn Hammond
  ! Date: 11/11/22

  use Material_Aux_module, only : material_auxvar_type, &
                                  MaterialAuxVarGetValue
  use String_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, PERMEABILITY_X, &
                               POROSITY, BASE_POROSITY, &
                               VG_ALPHA, VG_SR, VG_M, &
                               ARCHIE_CEMENTATION_EXPONENT, &
                               ARCHIE_SATURATION_EXPONENT, &
                               ARCHIE_TORTUOSITY_CONSTANT

  class(inversion_aux_type) :: aux
  PetscReal :: value
  PetscInt :: iparameter_type
  PetscInt :: imat
  type(material_auxvar_type) :: material_auxvar

  type(material_property_type), pointer :: material_property
  type(characteristic_curves_type), pointer :: cc

  select case(iparameter_type)
    case(ELECTRICAL_CONDUCTIVITY,ARCHIE_CEMENTATION_EXPONENT, &
         ARCHIE_SATURATION_EXPONENT,ARCHIE_TORTUOSITY_CONSTANT)
      value = MaterialAuxVarGetValue(material_auxvar,iparameter)
    case(ELECTRICAL_CONDUCTIVITY)
      value = MaterialAuxVarGetValue(material_auxvar,ELECTRICAL_CONDUCTIVITY)
    case(PERMEABILITY)
      value = MaterialAuxVarGetValue(material_auxvar,PERMEABILITY_X)
    case(POROSITY)
      value = MaterialAuxVarGetValue(material_auxvar,BASE_POROSITY)
    case(VG_ALPHA,VG_SR,VG_M)
      material_property => aux%material_property_array(imat)%ptr
      cc => aux%cc_array(material_property%saturation_function_id)%ptr
      select case(iparameter_type)
        case(VG_ALPHA)
          value = cc%saturation_function%GetAlpha_()
        case(VG_M)
          value = cc%saturation_function%GetM_()
        case(VG_SR)
          value = cc%saturation_function%GetResidualSaturation()
      end select
    case default
      call aux%driver%PrintErrMsg('Unrecognized variable in &
                                   &InvAuxGetParamValueByCell: ' // &
                                   trim(StringWrite(iparameter_type)))
  end select

end subroutine InvAuxGetParamValueByCell

! ************************************************************************** !

subroutine InvAuxScatGlobalToDistParam(inversion_aux,global_, &
                                       dist_parameter_vec,direction)
  !
  ! Scatters from work to dist_parameter_vec
  !
  ! Author: Glenn Hammond
  ! Date: 04/01/22
  !
  type(inversion_aux_type) :: inversion_aux
  Vec :: global_
  Vec :: dist_parameter_vec
  PetscInt :: direction

  PetscErrorCode :: ierr

  if (direction == INVAUX_SCATFORWARD) then
    call VecScatterBegin(inversion_aux%scatter_global_to_dist_param, &
                         global_,dist_parameter_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(inversion_aux%scatter_global_to_dist_param, &
                       global_,dist_parameter_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  else ! INVAUX_SCATREVERSE
    call VecScatterBegin(inversion_aux%scatter_global_to_dist_param, &
                         dist_parameter_vec,global_, &
                         INSERT_VALUES,SCATTER_REVERSE,ierr);CHKERRQ(ierr)
    call VecScatterEnd(inversion_aux%scatter_global_to_dist_param, &
                       dist_parameter_vec,global_, &
                       INSERT_VALUES,SCATTER_REVERSE,ierr);CHKERRQ(ierr)
  endif

end subroutine InvAuxScatGlobalToDistParam

! ************************************************************************** !

subroutine InvAuxScatParamToDistParam(inversion_aux,parameter_vec, &
                                      dist_parameter_vec,direction)
  !
  ! Scatters from parameter_vec to dist_parameter_vec
  !
  ! Author: Glenn Hammond
  ! Date: 04/11/22
  !
  type(inversion_aux_type) :: inversion_aux
  Vec :: parameter_vec
  Vec :: dist_parameter_vec
  PetscInt :: direction

  PetscErrorCode :: ierr

  if (direction == INVAUX_SCATFORWARD) then
    ! the parameter_vec is full on each process
    call VecScatterBegin(inversion_aux%scatter_param_to_dist_param, &
                         parameter_vec, &
                         dist_parameter_vec,INSERT_VALUES, &
                         SCATTER_FORWARD_LOCAL,ierr);CHKERRQ(ierr)
    call VecScatterEnd(inversion_aux%scatter_param_to_dist_param, &
                       parameter_vec, &
                       dist_parameter_vec,INSERT_VALUES, &
                       SCATTER_FORWARD_LOCAL,ierr);CHKERRQ(ierr)
  else ! INVAUX_SCATREVERSE
    call VecScatterBegin(inversion_aux%scatter_param_to_dist_param, &
                         dist_parameter_vec, &
                         parameter_vec,INSERT_VALUES, &
                         SCATTER_REVERSE,ierr);CHKERRQ(ierr)
    call VecScatterEnd(inversion_aux%scatter_param_to_dist_param, &
                       dist_parameter_vec, &
                       parameter_vec,INSERT_VALUES, &
                       SCATTER_REVERSE,ierr);CHKERRQ(ierr)
  endif

end subroutine InvAuxScatParamToDistParam

! ************************************************************************** !

subroutine InvAuxScatMeasToDistMeas(inversion_aux,measurement_vec, &
                                    dist_measurement_vec,direction)
  !
  ! Scatters from measurement_vec to dist_measurement_vec
  !
  ! Author: Glenn Hammond
  ! Date: 04/01/22
  !
  type(inversion_aux_type) :: inversion_aux
  Vec :: measurement_vec
  Vec :: dist_measurement_vec
  PetscInt :: direction

  PetscErrorCode :: ierr

  if (direction == INVAUX_SCATFORWARD) then
    ! the measurement_vec is full on each process
    call VecScatterBegin(inversion_aux%scatter_measure_to_dist_measure, &
                         measurement_vec, &
                         dist_measurement_vec,INSERT_VALUES, &
                         SCATTER_FORWARD_LOCAL,ierr);CHKERRQ(ierr)
    call VecScatterEnd(inversion_aux%scatter_measure_to_dist_measure, &
                       measurement_vec, &
                       dist_measurement_vec,INSERT_VALUES, &
                       SCATTER_FORWARD_LOCAL,ierr);CHKERRQ(ierr)
  else ! INVAUX_SCATREVERSE
    call VecScatterBegin(inversion_aux%scatter_measure_to_dist_measure, &
                         dist_measurement_vec, &
                         measurement_vec,INSERT_VALUES, &
                         SCATTER_REVERSE,ierr);CHKERRQ(ierr)
    call VecScatterEnd(inversion_aux%scatter_measure_to_dist_measure, &
                       dist_measurement_vec, &
                       measurement_vec,INSERT_VALUES, &
                       SCATTER_REVERSE,ierr);CHKERRQ(ierr)
  endif

end subroutine InvAuxScatMeasToDistMeas

! ************************************************************************** !

subroutine InvAuxCopyMeasToFromMeasVec(aux,idirection)
  !
  ! Copies parameter values back and forth
  !
  ! Author: Glenn Hammond
  ! Date: 03/30/22
  use Inversion_Measurement_Aux_module

  class(inversion_aux_type) :: aux
  PetscInt :: idirection

  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscErrorCode :: ierr

  call VecGetArrayF90(aux%measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)
  select case(idirection)
    case(INVAUX_COPY_FROM_VEC)
      do i = 1, size(aux%measurements)
        aux%measurements(i)%simulated_value = vec_ptr(i)
      enddo
    case(INVAUX_COPY_TO_VEC)
      do i = 1, size(aux%measurements)
        vec_ptr(i) = aux%measurements(i)%simulated_value
      enddo
    case default
      stop 'Error: idirection in InvAuxCopyMeasToFromMeasVec,Update'
  end select
  call VecRestoreArrayF90(aux%measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine InvAuxCopyMeasToFromMeasVec

! ************************************************************************** !

subroutine InvAuxBCastVecForCommI(comm,vec,driver)
  !
  ! Broadcasts the contents of a Vec segment to the perturbation ranks
  !
  ! Author: Glenn Hammond
  ! Date: 03/06/23
  !
  type(comm_type) :: comm
  Vec :: vec
  type(driver_type) :: driver

  PetscInt :: vec_size
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  call VecGetLocalSize(vec,vec_size,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  call MPI_BCast(vec_ptr,vec_size,MPI_DOUBLE_PRECISION,ZERO_INTEGER_MPI, &
                 comm%communicator,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine InvAuxBCastVecForCommI

! ************************************************************************** !

subroutine InversionAuxPerturbationStrip(perturbation)
  !
  ! Deallocates members of inversion perturbation
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Utility_module

  type(inversion_perturbation_type), pointer :: perturbation

  PetscErrorCode :: ierr

  if (.not.associated(perturbation)) return

  call DeallocateArray(perturbation%select_cells)
  if (perturbation%base_parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(perturbation%base_parameter_vec,ierr);CHKERRQ(ierr)
  endif
  if (perturbation%base_measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(perturbation%base_measurement_vec,ierr);CHKERRQ(ierr)
  endif
  deallocate(perturbation)
  nullify(perturbation)

end subroutine InversionAuxPerturbationStrip

! ************************************************************************** !

subroutine InversionAuxDestroy(aux)
  !
  ! Deallocates a inversion auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  use Utility_module, only : DeallocateArray

  type(inversion_aux_type), pointer :: aux

  PetscInt :: i
  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

  call DeallocateArray(aux%sync_times)

  ! these are owned and must be destroyed
  if (associated(aux%coupled_aux)) then
    call InversionCoupledAuxDestroy(aux%coupled_aux)
  endif
  if (associated(aux%measurements)) then
    do i = 1, size(aux%measurements)
      call InversionMeasurementAuxStrip(aux%measurements(i))
    enddo
    deallocate(aux%measurements)
  endif
  nullify(aux%measurements)
  if (associated(aux%parameters)) then
    do i = 1, size(aux%parameters)
      call InversionParameterStrip(aux%parameters(i))
    enddo
    deallocate(aux%parameters)
  endif
  nullify(aux%parameters)
  if (aux%measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(aux%measurement_vec,ierr);CHKERRQ(ierr)
  endif
  if (aux%dist_measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(aux%dist_measurement_vec,ierr);CHKERRQ(ierr)
  endif
  if (aux%parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(aux%parameter_vec,ierr);CHKERRQ(ierr)
  endif
  if (aux%dist_parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(aux%dist_parameter_vec,ierr);CHKERRQ(ierr)
  endif
  if (aux%del_parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(aux%del_parameter_vec,ierr);CHKERRQ(ierr)
  endif
  if (aux%scatter_measure_to_dist_measure /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(aux%scatter_measure_to_dist_measure, &
                           ierr);CHKERRQ(ierr)
  endif
  if (aux%scatter_param_to_dist_param /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(aux%scatter_measure_to_dist_measure, &
                           ierr);CHKERRQ(ierr)
  endif
  if (aux%scatter_global_to_dist_param /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(aux%scatter_measure_to_dist_measure, &
                           ierr);CHKERRQ(ierr)
  endif
  if (aux%JsensitivityT /= PETSC_NULL_MAT) then
    call MatDestroy(aux%JsensitivityT,ierr);CHKERRQ(ierr)
  endif

  ! nullify objects owned by other objects
  nullify(aux%driver)
  nullify(aux%material_property_array)
  nullify(aux%cc_array)
  aux%solution = PETSC_NULL_VEC
  nullify(aux%local_measurement_values_ptr)
  nullify(aux%local_dobs_dunknown_values_ptr)
  nullify(aux%local_dobs_dparam_values_ptr)
  ! adjoints
  call InvAuxAdjCleanupAfterForwardRun(aux)
  ! perturbation
  call InversionAuxPerturbationStrip(aux%perturbation)

  deallocate(aux)
  nullify(aux)

end subroutine InversionAuxDestroy

end module Inversion_Aux_module
