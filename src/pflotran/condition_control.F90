module Condition_Control_module

  ! This module store routines that operate on conditions from a level above
  ! that of the realization_module.  This is necessary to access capability
  ! such as HDF5 which is unavailable from within the realization object
  ! and below.  Routines in this module will loop over realization, levels,
  ! and patches without calling underlying level/patch versions of the
  ! subroutines, which is common in realization.F90 - GEH
#include "petsc/finclude/petscvec.h"
  use petscvec

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: CondControlAssignFlowInitCond, &
            CondControlAssignRTTranInitCond, &
            CondControlAssignNWTranInitCond, &
            CondControlScaleSourceSink

contains

! ************************************************************************** !

subroutine CondControlAssignFlowInitCond(realization)
  !
  ! Assigns flow initial conditions to model
  !
  ! Author: Glenn Hammond
  ! Date: 11/02/07, 10/18/11
  !
  use Realization_Subsurface_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Dataset_Base_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Common_HDF5_class
  use Dataset_module
  use Grid_module
  use Patch_module
  use EOS_Water_module

  use Global_module
  use Variables_module, only : STATE
  use Global_Aux_module
  use General_Aux_module, gen_dof_to_primary_variable => dof_to_primary_variable
  use WIPP_Flow_Aux_module, wf_dof_to_primary_variable => dof_to_primary_variable
  use Hydrate_Aux_module, hyd_dof_to_primary_variable => dof_to_primary_variable

  implicit none

  class(realization_subsurface_type) :: realization

  PetscInt :: icell, iconn, idof
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscReal, pointer :: xx_p(:)
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(flow_general_condition_type), pointer :: general
  type(flow_hydrate_condition_type), pointer :: hydrate
  class(dataset_base_type), pointer :: dataset
  PetscBool :: dataset_flag(realization%option%nflowdof)
  PetscInt :: num_connections
  PetscInt, pointer :: conn_id_ptr(:)
  PetscInt :: offset, istate
  PetscReal :: temperature, p_sat
  PetscReal :: tempreal

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  ! to catch uninitialized grid cells.  see VecMin check at bottom.
  call VecSet(field%work_loc,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
  call GlobalSetAuxVarVecLoc(realization,field%work_loc,STATE,ZERO_INTEGER)

  ! TODO(patch_list): fix indentation

    select case(option%iflowmode)

      case(WF_MODE)

        call VecGetArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

        xx_p = UNINITIALIZED_DOUBLE

        initial_condition => patch%initial_condition_list%first
        do

          if (.not.associated(initial_condition)) exit

          dataset_flag = PETSC_FALSE
          do idof = 1, option%nflowdof
            dataset =>  initial_condition%flow_condition% &
                              sub_condition_ptr(idof)%ptr%dataset
            select type(dataset_ptr => dataset)
              class is(dataset_gridded_hdf5_type)
                ! already mapped to flow_aux_real_var
                if (.not.associated(initial_condition%flow_aux_real_var)) then
                  option%io_buffer = 'A gridded dataset is being &
                    &used with WIPP_FLOW, yet flow_aux_real_var is not &
                    &allocated.'
                  call PrintErrMsgToDev(option,'')
                endif
              class is(dataset_common_hdf5_type)
                dataset_flag(idof) = PETSC_TRUE
                call VecRestoreArrayF90(field%flow_xx,xx_p, &
                                        ierr);CHKERRQ(ierr)
                call ConditionControlMapDatasetToVec(realization, &
                        initial_condition%flow_condition% &
                          sub_condition_ptr(idof)%ptr%dataset,idof, &
                        field%flow_xx,GLOBAL)
                call VecGetArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
            end select
          enddo

          if (.not.associated(initial_condition%flow_aux_real_var)) then
            if (.not.associated(initial_condition%flow_condition)) then
              option%io_buffer = 'Flow condition is NULL in initial condition'
              call PrintErrMsg(option)
            endif

            general => initial_condition%flow_condition%general

            string = 'in flow condition "' // &
              trim(initial_condition%flow_condition%name) // &
              '" within initial condition "' // &
              trim(initial_condition%flow_condition%name) // &
              '" must be of type Dirichlet or Hydrostatic'
            ! error checking.  the data must match the state
            if (.not. &
                (general%liquid_pressure%itype == DIRICHLET_BC .or. &
                  general%liquid_pressure%itype == HYDROSTATIC_BC)) then
              option%io_buffer = 'Liquid pressure ' // trim(string)
              call PrintErrMsg(option)
            endif
            if (.not. &
                (general%gas_saturation%itype == DIRICHLET_BC .or. &
                  general%gas_saturation%itype == HYDROSTATIC_BC)) then
              option%io_buffer = 'Gas saturation ' // trim(string)
              call PrintErrMsg(option)
            endif

            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = grid%nL2G(local_id)
              iend = local_id*option%nflowdof
              ibegin = iend-option%nflowdof+1
              if (patch%imat(ghosted_id) <= 0) then
                xx_p(ibegin:iend) = 0.d0
                patch%aux%Global%auxvars(ghosted_id)%istate = 0
                cycle
              endif
              ! decrement ibegin to give a local offset of 0
              ibegin = ibegin - 1
              if (.not.dataset_flag(WIPPFLO_LIQUID_PRESSURE_DOF)) then
                xx_p(ibegin+WIPPFLO_LIQUID_PRESSURE_DOF) = &
                  general%liquid_pressure%dataset%rarray(1)
              endif
              if (.not.dataset_flag(WIPPFLO_GAS_SATURATION_DOF)) then
                xx_p(ibegin+WIPPFLO_GAS_SATURATION_DOF) = &
                  general%gas_saturation%dataset%rarray(1)
              endif
              patch%aux%Global%auxvars(ghosted_id)%istate = &
                initial_condition%flow_condition%iphase
            enddo
          else
            do iconn=1,initial_condition%connection_set%num_connections
              local_id = initial_condition%connection_set%id_dn(iconn)
              ghosted_id = grid%nL2G(local_id)
              if (patch%imat(ghosted_id) <= 0) then
                iend = local_id*option%nflowdof
                ibegin = iend-option%nflowdof+1
                xx_p(ibegin:iend) = 0.d0
                patch%aux%Global%auxvars(ghosted_id)%istate = 0
                cycle
              endif
              offset = (local_id-1)*option%nflowdof
              istate = initial_condition%flow_aux_int_var(1,iconn)
              do idof = 1, option%nflowdof
                if (dataset_flag(idof)) cycle
                xx_p(offset+idof) = &
                  initial_condition%flow_aux_real_var( &
                    initial_condition%flow_aux_mapping( &
                      wf_dof_to_primary_variable(idof)),iconn)
              enddo
              patch%aux%Global%auxvars(ghosted_id)%istate = istate
            enddo
          endif
          initial_condition => initial_condition%next
        enddo

        call VecRestoreArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

      case(G_MODE) ! general phase mode

        call VecGetArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

        xx_p = UNINITIALIZED_DOUBLE

        initial_condition => patch%initial_condition_list%first
        do

          if (.not.associated(initial_condition)) exit

          if (.not.associated(initial_condition%flow_aux_real_var)) then
            if (.not.associated(initial_condition%flow_condition)) then
              option%io_buffer = 'Flow condition is NULL in initial condition'
              call PrintErrMsg(option)
            endif

            general => initial_condition%flow_condition%general

            string = 'in flow condition "' // &
              trim(initial_condition%flow_condition%name) // &
              '" within initial condition "' // &
              trim(initial_condition%flow_condition%name) // &
              '" must be of type Dirichlet or Hydrostatic'
            ! error checking.  the data must match the state
            select case(initial_condition%flow_condition%iphase)
              case(TWO_PHASE_STATE)
                if (.not. &
                    (general%gas_pressure%itype == DIRICHLET_BC .or. &
                      general%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (general%gas_saturation%itype == DIRICHLET_BC .or. &
                      general%gas_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
              case(LGP_STATE)
                if (.not. &
                     (general%gas_pressure%itype == DIRICHLET_BC .or. &
                     general%gas_pressure%itype == HYDROSTATIC_BC)) then
                   option%io_buffer = 'Gas pressure ' // trim(string)
                   call PrintErrMsg(option)
                endif
                if (.not. &
                     (general%gas_saturation%itype == DIRICHLET_BC .or. &
                     general%gas_saturation%itype == HYDROSTATIC_BC)) then
                   option%io_buffer = 'Gas saturation ' // trim(string)
                   call PrintErrMsg(option)
                endif
                if (general_salt .and. general_soluble_matrix) then
                   if (.not. &
                        (general%porosity%itype == DIRICHLET_BC .or. &
                        general%porosity%itype == HYDROSTATIC_BC)) then
                      option%io_buffer = 'Porosity ' // trim(string)
                      call PrintErrMsg(option)
                   endif
                endif
                if (general_salt) then
                   if (.not. &
                        (general%precipitate_saturation%itype == DIRICHLET_BC .or. &
                        general%precipitate_saturation%itype == HYDROSTATIC_BC)) then
                      option%io_buffer = 'Precipitate saturation ' // trim(string)
                      call PrintErrMsg(option)
                   endif
                endif
              case(LIQUID_STATE)
                if (.not. &
                    (general%liquid_pressure%itype == DIRICHLET_BC .or. &
                      general%liquid_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Liquid pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (general%mole_fraction%itype == DIRICHLET_BC .or. &
                      general%mole_fraction%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Mole fraction ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (general_salt .and. general_soluble_matrix) then
                  if (.not. &
                       (general%porosity%itype == DIRICHLET_BC .or. &
                         general%porosity%itype == HYDROSTATIC_BC)) then
                    option%io_buffer = 'Porosity ' // trim(string)
                    call PrintErrMsg(option)
                  endif
                endif
                if (general_salt) then
                  if (.not. &
                       (general%salt_mole_fraction%itype == DIRICHLET_BC .or. &
                         general%salt_mole_fraction%itype == HYDROSTATIC_BC)) then
                    option%io_buffer = 'Salt mole fraction ' // trim(string)
                    call PrintErrMsg(option)
                  endif
                endif
              case(LP_STATE)
                if (.not. &
                    (general%liquid_pressure%itype == DIRICHLET_BC .or. &
                      general%liquid_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Liquid pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (general%mole_fraction%itype == DIRICHLET_BC .or. &
                      general%mole_fraction%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Mole fraction ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (general_salt .and. general_soluble_matrix) then
                  if (.not. &
                       (general%porosity%itype == DIRICHLET_BC .or. &
                         general%porosity%itype == HYDROSTATIC_BC)) then
                    option%io_buffer = 'Porosity ' // trim(string)
                    call PrintErrMsg(option)
                  endif
                endif
                if (general_salt .and. .not. general_soluble_matrix) then
                  if (.not. &
                       (general%precipitate_saturation%itype == DIRICHLET_BC .or. &
                         general%precipitate_saturation%itype == HYDROSTATIC_BC)) then
                    option%io_buffer = 'Precipitate saturation ' // trim(string)
                    call PrintErrMsg(option)
                  endif
                endif
              case(GAS_STATE)
                if (.not. &
                    (general%gas_pressure%itype == DIRICHLET_BC .or. &
                      general%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (general%mole_fraction%itype == DIRICHLET_BC .or. &
                      general%mole_fraction%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
             case(GP_STATE)
                if (.not. &
                     (general%gas_pressure%itype == DIRICHLET_BC .or. &
                     general%gas_pressure%itype == HYDROSTATIC_BC)) then
                   option%io_buffer = 'Gas pressure ' // trim(string)
                   call PrintErrMsg(option)
                endif
                if (general_salt .and. .not. general_soluble_matrix) then
                  if (.not. &
                       (general%mole_fraction%itype == DIRICHLET_BC .or. &
                        general%mole_fraction%itype == HYDROSTATIC_BC)) then
                       option%io_buffer = 'Gas saturation ' // trim(string)
                     call PrintErrMsg(option)
                  endif
                elseif (general_salt .and. general_soluble_matrix) then
                  if (.not. &
                       (general%porosity%itype == DIRICHLET_BC .or. &
                        general%porosity%itype == HYDROSTATIC_BC)) then
                        option%io_buffer = 'Porosity ' // trim(string)
                    call PrintErrMsg(option)
                  endif
                endif
            end select
            if (.not. &
                (general%temperature%itype == DIRICHLET_BC .or. &
                  general%temperature%itype == HYDROSTATIC_BC)) then
              option%io_buffer = 'Temperature ' // trim(string)
              call PrintErrMsg(option)
            endif


            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = grid%nL2G(local_id)
              iend = local_id*option%nflowdof
              ibegin = iend-option%nflowdof+1
              if (patch%imat(ghosted_id) <= 0) then
                xx_p(ibegin:iend) = 0.d0
                patch%aux%Global%auxvars(ghosted_id)%istate = 0
                cycle
              endif
              ! decrement ibegin to give a local offset of 0
              ibegin = ibegin - 1
              select case(initial_condition%flow_condition%iphase)
                case(TWO_PHASE_STATE)
                  xx_p(ibegin+GENERAL_GAS_PRESSURE_DOF) = &
                    general%gas_pressure%dataset%rarray(1)
                  xx_p(ibegin+GENERAL_GAS_SATURATION_DOF) = &
                    general%gas_saturation%dataset%rarray(1)
                  temperature = general%temperature%dataset%rarray(1)
                  if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
                    xx_p(ibegin+GENERAL_ENERGY_DOF) = temperature
                  else
                    call EOSWaterSaturationPressure(temperature,p_sat,ierr)
                    ! p_a = p_g - p_s(T)
                    xx_p(ibegin+GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
                      general%gas_pressure%dataset%rarray(1) - &
                      p_sat
                  endif
                case(LGP_STATE)
                  xx_p(ibegin+GENERAL_GAS_PRESSURE_DOF) = &
                       general%gas_pressure%dataset%rarray(1)
                  xx_p(ibegin+GENERAL_GAS_SATURATION_DOF) = &
                       general%gas_saturation%dataset%rarray(1)
                  temperature = general%temperature%dataset%rarray(1)
                  if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
                     xx_p(ibegin+GENERAL_ENERGY_DOF) = temperature
                  else
                     call EOSWaterSaturationPressure(temperature,p_sat,ierr)
                     ! p_a = p_g - p_s(T)
                     xx_p(ibegin+GENERAL_2PH_STATE_AIR_PRESSURE_DOF) = &
                          general%gas_pressure%dataset%rarray(1) - &
                          p_sat
                  endif
                  if (general_salt .and. .not. general_soluble_matrix) then
                    xx_p(ibegin+GENERAL_PRECIPITATE_SAT_DOF) = &
                         general%precipitate_saturation%dataset%rarray(1)
                  elseif (general_salt .and. general_soluble_matrix) then
                    xx_p(ibegin+GENERAL_POROSITY_DOF) = &
                         general%porosity%dataset%rarray(1)
                  endif
                case(LIQUID_STATE)
                  xx_p(ibegin+GENERAL_LIQUID_PRESSURE_DOF) = &
                    general%liquid_pressure%dataset%rarray(1)
                  xx_p(ibegin+GENERAL_LIQUID_STATE_X_MOLE_DOF) = &
                    general%mole_fraction%dataset%rarray(1)
                  xx_p(ibegin+GENERAL_ENERGY_DOF) = &
                    general%temperature%dataset%rarray(1)
                  if (general_salt) then
                    xx_p(ibegin+GENERAL_LIQUID_STATE_S_MOLE_DOF) = &
                      general%salt_mole_fraction%dataset%rarray(1)
                  endif
                case(GAS_STATE)
                  xx_p(ibegin+GENERAL_GAS_PRESSURE_DOF) = &
                    general%gas_pressure%dataset%rarray(1)
                  if (general_gas_air_mass_dof == &
                      GENERAL_AIR_PRESSURE_INDEX) then
                    xx_p(ibegin+GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
                      general%gas_pressure%dataset%rarray(1) * &
                      general%mole_fraction%dataset%rarray(1)
                  else
                    xx_p(ibegin+GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
                      general%mole_fraction%dataset%rarray(1)
                  endif
                  xx_p(ibegin+GENERAL_ENERGY_DOF) = &
                    general%temperature%dataset%rarray(1)
                case(LP_STATE)
                  xx_p(ibegin+GENERAL_LIQUID_PRESSURE_DOF) = &
                       general%liquid_pressure%dataset%rarray(1)
                  if (general_salt .and. .not. general_soluble_matrix) then
                    xx_p(ibegin+GENERAL_PRECIPITATE_SAT_DOF) = &
                         general%precipitate_saturation%dataset%rarray(1)
                  elseif (general_salt .and. general_soluble_matrix) then
                     xx_p(ibegin+GENERAL_POROSITY_DOF) = &
                          general%porosity%dataset%rarray(1)
                  endif
                  xx_p(ibegin+GENERAL_ENERGY_DOF) = &
                       general%temperature%dataset%rarray(1)
                  xx_p(ibegin+GENERAL_LIQUID_STATE_X_MOLE_DOF) = &
                     general%mole_fraction%dataset%rarray(1)
                case(GP_STATE)
                  xx_p(ibegin+GENERAL_GAS_PRESSURE_DOF) = &
                       general%gas_pressure%dataset%rarray(1)
                  if (general_gas_air_mass_dof == &
                       GENERAL_AIR_PRESSURE_INDEX) then
                     xx_p(ibegin+GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
                          general%gas_pressure%dataset%rarray(1) * &
                          general%mole_fraction%dataset%rarray(1)
                  else
                     xx_p(ibegin+GENERAL_GAS_STATE_AIR_PRESSURE_DOF) = &
                          general%mole_fraction%dataset%rarray(1)
                  endif
                  xx_p(ibegin+GENERAL_ENERGY_DOF) = &
                       general%temperature%dataset%rarray(1)
                  if (general_salt .and. .not. general_soluble_matrix) then
                    xx_p(ibegin+GENERAL_PRECIPITATE_SAT_DOF) = &
                         general%precipitate_saturation%dataset%rarray(1)
                  elseif (general_salt .and. general_soluble_matrix) then
                    xx_p(ibegin+GENERAL_POROSITY_DOF) = &
                         general%porosity%dataset%rarray(1)
                  endif
              end select
              patch%aux%Global%auxvars(ghosted_id)%istate = &
                initial_condition%flow_condition%iphase
            enddo
          else
            do iconn=1,initial_condition%connection_set%num_connections
              local_id = initial_condition%connection_set%id_dn(iconn)
              ghosted_id = grid%nL2G(local_id)
              if (patch%imat(ghosted_id) <= 0) then
                iend = local_id*option%nflowdof
                ibegin = iend-option%nflowdof+1
                xx_p(ibegin:iend) = 0.d0
                patch%aux%Global%auxvars(ghosted_id)%istate = 0
                cycle
              endif
              offset = (local_id-1)*option%nflowdof
              istate = initial_condition%flow_aux_int_var(1,iconn)
              do idof = 1, option%nflowdof
                xx_p(offset+idof) = &
                  initial_condition%flow_aux_real_var( &
                    initial_condition%flow_aux_mapping( &
                      gen_dof_to_primary_variable(idof,istate)),iconn)
              enddo
              patch%aux%Global%auxvars(ghosted_id)%istate = istate
            enddo
          endif
          initial_condition => initial_condition%next
        enddo

        call VecRestoreArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

      case(H_MODE)

        call VecGetArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

        xx_p = UNINITIALIZED_DOUBLE

        initial_condition => patch%initial_condition_list%first
        do

          if (.not.associated(initial_condition)) exit

          if (.not.associated(initial_condition%flow_aux_real_var)) then
            if (.not.associated(initial_condition%flow_condition)) then
              option%io_buffer = 'Flow condition is NULL in initial condition'
              call PrintErrMsg(option)
            endif

            hydrate => initial_condition%flow_condition%hydrate

            string = 'in flow condition "' // &
              trim(initial_condition%flow_condition%name) // &
              '" within initial condition "' // &
              trim(initial_condition%flow_condition%name) // &
              '" must be of type Dirichlet or Hydrostatic'
            ! error checking.  the data must match the state
            select case(initial_condition%flow_condition%iphase)
              case(L_STATE)
                if (.not. &
                    (hydrate%liquid_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%liquid_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Liquid pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%mole_fraction%itype == DIRICHLET_BC .or. &
                      hydrate%mole_fraction%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Mole fraction ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%temperature%itype == DIRICHLET_BC .or. &
                  hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif

              case(G_STATE)
                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%mole_fraction%itype == DIRICHLET_BC .or. &
                      hydrate%mole_fraction%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%temperature%itype == DIRICHLET_BC .or. &
                  hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif
              case(H_STATE)
                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%temperature%itype == DIRICHLET_BC .or. &
                  hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif
              case(I_STATE)
                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%temperature%itype == DIRICHLET_BC .or. &
                  hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif
              case(GA_STATE)
                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%gas_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%gas_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%temperature%itype == DIRICHLET_BC .or. &
                  hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif

              case(HG_STATE)
                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%gas_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%gas_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%temperature%itype == DIRICHLET_BC .or. &
                  hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif

              case(HA_STATE)
                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%hydrate_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%hydrate_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Hydrate saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%temperature%itype == DIRICHLET_BC .or. &
                  hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif

              case(HI_STATE)
                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%hydrate_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%hydrate_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Hydrate saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%temperature%itype == DIRICHLET_BC .or. &
                  hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif

              case(GI_STATE)
                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%ice_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%ice_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Ice saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%temperature%itype == DIRICHLET_BC .or. &
                  hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif

              case(AI_STATE)
                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%mole_fraction%itype == DIRICHLET_BC .or. &
                      hydrate%mole_fraction%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Liquid Mole Fraction ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%liquid_saturation%itype == DIRICHLET_BC .or. &
                  hydrate%liquid_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Liquid saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif

              case(HGA_STATE)
                if (.not. &
                    (hydrate%liquid_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%liquid_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Liquid saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%hydrate_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%hydrate_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Hydrate saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                   (hydrate%temperature%itype == DIRICHLET_BC .or. &
                  hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif

              case(HAI_STATE)

                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%liquid_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%liquid_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Liquid saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%ice_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%ice_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Ice saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif

              case(HGI_STATE)

                if (.not. &
                    (hydrate%ice_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%ice_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Ice saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%hydrate_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%hydrate_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Hydrate saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%temperature%itype == DIRICHLET_BC .or. &
                      hydrate%temperature%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Temperature ' // trim(string)
                  call PrintErrMsg(option)
                endif


              case(GAI_STATE)

                if (.not. &
                    (hydrate%gas_pressure%itype == DIRICHLET_BC .or. &
                      hydrate%gas_pressure%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas pressure ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%liquid_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%liquid_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Liquid saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%ice_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%ice_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Ice saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif


              case(HGAI_STATE)

                if (.not. &
                    (hydrate%gas_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%gas_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Gas saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%liquid_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%liquid_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Liquid saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif
                if (.not. &
                    (hydrate%ice_saturation%itype == DIRICHLET_BC .or. &
                      hydrate%ice_saturation%itype == HYDROSTATIC_BC)) then
                  option%io_buffer = 'Ice saturation ' // trim(string)
                  call PrintErrMsg(option)
                endif


            end select


            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = grid%nL2G(local_id)
              iend = local_id*option%nflowdof
              ibegin = iend-option%nflowdof+1
              if (patch%imat(ghosted_id) <= 0) then
                xx_p(ibegin:iend) = 0.d0
                patch%aux%Global%auxvars(ghosted_id)%istate = 0
                cycle
              endif
              ! decrement ibegin to give a local offset of 0
              ibegin = ibegin - 1
              select case(initial_condition%flow_condition%iphase)
                case(L_STATE)
                  xx_p(ibegin+HYDRATE_LIQUID_PRESSURE_DOF) = &
                    hydrate%liquid_pressure%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_L_STATE_X_MOLE_DOF) = &
                    hydrate%mole_fraction%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%temperature%dataset%rarray(1)
                case(G_STATE)
                  xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                    hydrate%gas_pressure%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_G_STATE_AIR_PRESSURE_DOF) = &
                    hydrate%gas_pressure%dataset%rarray(1) * &
                    hydrate%mole_fraction%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%temperature%dataset%rarray(1)
                case(H_STATE)
                  xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                    hydrate%gas_pressure%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    MOL_RATIO_METH
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%temperature%dataset%rarray(1)
                case(I_STATE)
                  xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                    hydrate%gas_pressure%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    0.d0
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%temperature%dataset%rarray(1)
                case(GA_STATE)
                  if (associated(hydrate%gas_pressure)) then
                    xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%gas_pressure%dataset%rarray(1)
                  else
                    xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%liquid_pressure%dataset%rarray(1)
                  endif
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%gas_saturation%dataset%rarray(1)
                  temperature = hydrate%temperature%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = temperature
                case(HG_STATE)
                  xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                    hydrate%gas_pressure%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%gas_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%temperature%dataset%rarray(1)
                case(HA_STATE)
                  if (associated(hydrate%gas_pressure)) then
                    xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%gas_pressure%dataset%rarray(1)
                  else
                    xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%liquid_pressure%dataset%rarray(1)
                  endif
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%hydrate_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%temperature%dataset%rarray(1)
                case(HI_STATE)
                  xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                    hydrate%gas_pressure%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%hydrate_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%temperature%dataset%rarray(1)
                case(GI_STATE)
                  xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                    hydrate%gas_pressure%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%ice_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%temperature%dataset%rarray(1)
                case(AI_STATE)
                  if (associated(hydrate%gas_pressure)) then
                    xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%gas_pressure%dataset%rarray(1)
                  else
                    xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%liquid_pressure%dataset%rarray(1)
                  endif
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%mole_fraction%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%liquid_saturation%dataset%rarray(1)
                case(HGA_STATE)
                  xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%liquid_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%hydrate_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%temperature%dataset%rarray(1)
                case(HAI_STATE)
                  if (associated(hydrate%gas_pressure)) then
                    xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%gas_pressure%dataset%rarray(1)
                  else
                    xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%liquid_pressure%dataset%rarray(1)
                  endif
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%liquid_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%ice_saturation%dataset%rarray(1)
                case(HGI_STATE)
                  xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%ice_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%hydrate_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%temperature%dataset%rarray(1)
                case(GAI_STATE)
                  if (associated(hydrate%gas_pressure)) then
                    xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%gas_pressure%dataset%rarray(1)
                  else
                    xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%liquid_pressure%dataset%rarray(1)
                  endif
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%liquid_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%ice_saturation%dataset%rarray(1)
                case(HGAI_STATE)
                  xx_p(ibegin+HYDRATE_GAS_PRESSURE_DOF) = &
                      hydrate%liquid_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_GAS_SATURATION_DOF) = &
                    hydrate%gas_saturation%dataset%rarray(1)
                  xx_p(ibegin+HYDRATE_ENERGY_DOF) = &
                    hydrate%ice_saturation%dataset%rarray(1)
              end select
              patch%aux%Global%auxvars(ghosted_id)%istate = &
                initial_condition%flow_condition%iphase
            enddo
          else
            do iconn=1,initial_condition%connection_set%num_connections
              local_id = initial_condition%connection_set%id_dn(iconn)
              ghosted_id = grid%nL2G(local_id)
              if (patch%imat(ghosted_id) <= 0) then
                iend = local_id*option%nflowdof
                ibegin = iend-option%nflowdof+1
                xx_p(ibegin:iend) = 0.d0
                patch%aux%Global%auxvars(ghosted_id)%istate = 0
                cycle
              endif
              offset = (local_id-1)*option%nflowdof
              istate = initial_condition%flow_aux_int_var(1,iconn)
              do idof = 1, option%nflowdof
                xx_p(offset+idof) = &
                  initial_condition%flow_aux_real_var( &
                    initial_condition%flow_aux_mapping( &
                      hyd_dof_to_primary_variable(idof,istate)),iconn)
              enddo
              patch%aux%Global%auxvars(ghosted_id)%istate = istate
            enddo
          endif
          initial_condition => initial_condition%next
        enddo

        call VecRestoreArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

      case default
        ! assign initial conditions values to domain
        call VecGetArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

        xx_p = UNINITIALIZED_DOUBLE

        initial_condition => patch%initial_condition_list%first
        do

          if (.not.associated(initial_condition)) exit

          dataset_flag = PETSC_FALSE
          do idof = 1, option%nflowdof
            dataset =>  initial_condition%flow_condition% &
                              sub_condition_ptr(idof)%ptr%dataset
            select type(dataset_ptr => dataset)
              class is(dataset_gridded_hdf5_type)
                ! already mapped to flow_aux_real_var
                if (.not.associated(initial_condition%flow_aux_real_var)) then
                  option%io_buffer = 'A gridded dataset is being &
                    &used with ' // trim(option%flowmode) // &
                    ', yet flow_aux_real_var is not allocated.'
                  call PrintErrMsgToDev(option,'')
                endif
              class is(dataset_common_hdf5_type)
                dataset_flag(idof) = PETSC_TRUE
                call VecRestoreArrayF90(field%flow_xx,xx_p, &
                                        ierr);CHKERRQ(ierr)
                call ConditionControlMapDatasetToVec(realization, &
                        initial_condition%flow_condition% &
                          sub_condition_ptr(idof)%ptr%dataset,idof, &
                        field%flow_xx,GLOBAL)
                call VecGetArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
            end select
          enddo
          if (.not.associated(initial_condition%flow_aux_real_var) .and. &
              .not.associated(initial_condition%flow_condition)) then
            option%io_buffer = 'Flow condition is NULL in initial condition'
            call PrintErrMsg(option)
          endif
          if (associated(initial_condition%flow_aux_real_var)) then
            num_connections = &
              initial_condition%connection_set%num_connections
            conn_id_ptr => initial_condition%connection_set%id_dn
          else
            num_connections = initial_condition%region%num_cells
            conn_id_ptr => initial_condition%region%cell_ids
          endif
          do iconn=1, num_connections
            local_id = conn_id_ptr(iconn)
            ghosted_id = grid%nL2G(local_id)
            iend = local_id*option%nflowdof
            ibegin = iend-option%nflowdof+1
            if (patch%imat(ghosted_id) <= 0) then
              xx_p(ibegin:iend) = 0.d0
              patch%aux%Global%auxvars(ghosted_id)%istate = 0
              cycle
            endif
            if (associated(initial_condition%flow_aux_real_var)) then
              do idof = 1, option%nflowdof
                if (.not.dataset_flag(idof)) then
                  xx_p(ibegin+idof-1) =  &
                    initial_condition%flow_aux_real_var(idof,iconn)
                endif
              enddo
            else
              do idof = 1, option%nflowdof
                if (.not.dataset_flag(idof)) then
                  xx_p(ibegin+idof-1) = &
                    initial_condition%flow_condition% &
                      sub_condition_ptr(idof)%ptr%dataset%rarray(1)
                endif
              enddo
            endif
            patch%aux%Global%auxvars(ghosted_id)%istate = &
              initial_condition%flow_condition%iphase
          enddo
          initial_condition => initial_condition%next
        enddo
        call VecRestoreArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

    end select

  select case(option%iflowmode)
    case(RICHARDS_MODE,RICHARDS_TS_MODE,ZFLOW_MODE,PNF_MODE)
    case default
      call GlobalUpdateState(realization)
  end select

  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                   field%flow_xx_loc,NFLOWDOF)

  call VecCopy(field%flow_xx,field%flow_yy,ierr);CHKERRQ(ierr)

  ! cannot perform VecMin on local vector as the ghosted corner values are not
  ! updated during the local to local update.
  call GlobalGetAuxVarVecLoc(realization,field%work_loc,STATE)
  call DiscretizationLocalToGlobal(discretization,field%work_loc,field%work, &
                                   ONEDOF)
  call VecMin(field%work,PETSC_NULL_INTEGER,tempreal,ierr);CHKERRQ(ierr)
  if (tempreal < 0.d0) then
!    print *, tempreal
    option%io_buffer = 'Uninitialized cells in domain.'
    call PrintErrMsg(option)
  endif

end subroutine CondControlAssignFlowInitCond

! ************************************************************************** !

subroutine CondControlAssignRTTranInitCond(realization)
  !
  ! Assigns transport initial conditions to model
  !
  ! Author: Glenn Hammond
  ! Date: 11/02/07, 10/18/11
  !

  use Realization_Subsurface_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Transport_Constraint_Base_module
  use Transport_Constraint_RT_module
  use Grid_module
  use Dataset_Base_class
  use Patch_module
  use Reactive_Transport_module, only : RTUpdateAuxVars, &
                                        RTUpdateActivityCoefficients
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use Reaction_module
  use HDF5_module
  use Secondary_Continuum_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization

  PetscInt :: icell, idof, temp_int, iimmobile, cell
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscInt :: irxn, isite, imnrl, ikinrxn
  PetscReal, pointer :: xx_p(:), xx_loc_p(:), vec_p(:)
  Vec :: vec1_loc
  Vec :: vec2_loc
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  class(reaction_rt_type), pointer :: reaction
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(tran_constraint_coupler_rt_type), pointer :: constraint_coupler
  class(tran_constraint_coupler_rt_type), pointer :: sec_constraint_coupler
  class(tran_constraint_rt_type), pointer :: constraint
  type(material_auxvar_type), pointer :: material_auxvars(:)
  class(tran_constraint_rt_type), pointer :: sec_tran_constraint
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)

  PetscInt :: iphase
  PetscInt :: offset
  PetscBool :: equilibrate_at_each_cell
  character(len=MAXSTRINGLENGTH) :: string, string2
  class(dataset_base_type), pointer :: dataset
  PetscInt :: aq_dataset_to_idof(realization%reaction%naqcomp)
  PetscInt :: iaqdataset, num_aq_datasets
  PetscBool :: use_aq_dataset
  PetscReal :: ave_num_iterations
  PetscReal :: tempreal
  PetscReal, parameter :: epsilon = 1.d-16
  PetscInt :: prev_equilibrated_ghosted_id
  PetscReal, pointer :: flow_xx_p(:)
  PetscLogDouble :: tstart, tend

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  reaction => realization%reaction

  iphase = 1
  vec1_loc = PETSC_NULL_VEC
  vec2_loc = PETSC_NULL_VEC

  ! TODO(patch_list) fix indentation

    rt_auxvars => patch%aux%RT%auxvars
    global_auxvars => patch%aux%Global%auxvars
    material_auxvars => patch%aux%Material%auxvars

    ! assign initial conditions values to domain
    call VecGetArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
    select case(option%iflowmode)
      case(MPH_MODE)
        call VecGetArrayF90(field%flow_xx,flow_xx_p,ierr);CHKERRQ(ierr)
    end select

    xx_p = UNINITIALIZED_DOUBLE

    initial_condition => patch%initial_condition_list%first
    do

      if (.not.associated(initial_condition)) exit

      constraint_coupler => &
        TranConstraintCouplerRTCast(initial_condition%tran_condition% &
                                      cur_constraint_coupler)
      constraint => TranConstraintRTCast(constraint_coupler%constraint)
      if (option%use_sc) then
        rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars
        if (associated(initial_condition%tran_condition%sec_constraint_coupler)) then
          sec_constraint_coupler => TranConstraintCouplerRTCast(initial_condition%tran_condition% &
                                                                  sec_constraint_coupler)
          sec_tran_constraint => TranConstraintRTCast(sec_constraint_coupler%constraint)
        else
          sec_constraint_coupler => constraint_coupler
          sec_tran_constraint => constraint
        endif
      endif

      equilibrate_at_each_cell = constraint_coupler%equilibrate_at_each_cell
      use_aq_dataset = PETSC_FALSE
      num_aq_datasets = 0
      aq_dataset_to_idof = 0
      do idof = 1, reaction%naqcomp ! primary aqueous concentrations
        if (constraint%aqueous_species%external_dataset(idof)) then
          num_aq_datasets = num_aq_datasets + 1
          aq_dataset_to_idof(num_aq_datasets) = idof
          equilibrate_at_each_cell = PETSC_TRUE
          use_aq_dataset = PETSC_TRUE
          string = 'constraint ' // trim(constraint%name)
          dataset => DatasetBaseGetPointer(realization%datasets, &
                        constraint%aqueous_species%constraint_aux_string(idof), &
                        string,option)
          call ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                                field%tran_xx_loc,LOCAL)
        endif
      enddo

      ! read in heterogeneous mineral volume fractions
      if (associated(constraint%minerals)) then
        do imnrl = 1, reaction%mineral%nkinmnrl
          if (constraint%minerals%external_vol_frac_dataset(imnrl)) then
            equilibrate_at_each_cell = PETSC_TRUE
            string = 'constraint ' // trim(constraint%name)
            dataset => DatasetBaseGetPointer(realization%datasets, &
                          constraint%minerals% &
                            constraint_vol_frac_string(imnrl), &
                          string,option)
            if (vec1_loc == PETSC_NULL_VEC) then
              ! cannot use field%work_loc as it is used within ConditionCo...
              call VecDuplicate(field%work_loc,vec1_loc,ierr);CHKERRQ(ierr)
            endif
            idof = ONE_INTEGER
            call ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                                 vec1_loc,LOCAL)
            call VecGetArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = grid%nL2G(local_id)
              rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl) = vec_p(ghosted_id)
              rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl) = vec_p(ghosted_id)
            enddo
            call VecRestoreArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
          endif
        enddo
      endif

      ! read in heterogeneous mineral surface area
      if (associated(constraint%minerals)) then
        do imnrl = 1, reaction%mineral%nkinmnrl
          if (constraint%minerals%external_area_dataset(imnrl)) then
            equilibrate_at_each_cell = PETSC_TRUE
            string = 'constraint ' // trim(constraint%name)
            dataset => DatasetBaseGetPointer(realization%datasets, &
                          constraint%minerals% &
                          constraint_area_string(imnrl), &
                          string,option)
            if (vec1_loc == PETSC_NULL_VEC) then
              ! cannot use field%work_loc as it is used within ConditionCo...
              call VecDuplicate(field%work_loc,vec1_loc,ierr);CHKERRQ(ierr)
            endif
            idof = ONE_INTEGER
            call ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                                 vec1_loc,LOCAL)
            call VecScale(vec1_loc, &
                          constraint%minerals% &
                            constraint_area_conv_factor(imnrl), &
                          ierr);CHKERRQ(ierr)
            if (constraint%minerals%area_per_unit_mass(imnrl)) then
              if (constraint%minerals% &
                    external_vol_frac_dataset(imnrl)) then
                dataset => DatasetBaseGetPointer(realization%datasets, &
                              constraint%minerals% &
                                constraint_vol_frac_string(imnrl), &
                              string,option)
                if (vec2_loc == PETSC_NULL_VEC) then
                  call VecDuplicate(vec1_loc,vec2_loc,ierr);CHKERRQ(ierr)
                endif
                idof = ONE_INTEGER
                call ConditionControlMapDatasetToVec(realization,dataset, &
                                                     idof,vec2_loc,LOCAL)
                call VecPointwiseMult(vec1_loc,vec1_loc,vec2_loc, &
                                      ierr);CHKERRQ(ierr)
              else
                call VecScale(vec1_loc, &
                              constraint%minerals%constraint_vol_frac(imnrl), &
                              ierr);CHKERRQ(ierr)
              endif
              call DiscretizationLocalToGlobal(discretization,vec1_loc, &
                                               field%work,ONEDOF)
              call VecMin(field%work,PETSC_NULL_INTEGER,tempreal, &
                          ierr);CHKERRQ(ierr)
              if (tempreal < epsilon) then
                option%io_buffer = 'A zero volume fraction assigned to &
                  &mineral "' // &
                  trim(reaction%mineral%kinmnrl_names(imnrl)) // &
                  '" in constraint "' // trim(constraint%name) // &
                  '" prevents the use of a mass-based surface area in the &
                  &constraint.'
                call PrintErrMsg(option)
              endif
            endif
            call VecGetArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = grid%nL2G(local_id)
              rt_auxvars(ghosted_id)%mnrl_area0(imnrl) = vec_p(ghosted_id)
              rt_auxvars(ghosted_id)%mnrl_area(imnrl) = vec_p(ghosted_id)
            enddo
            call VecRestoreArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
          endif
        enddo
      endif

      ! read in heterogeneous immobile
      if (associated(constraint%immobile_species)) then
        do iimmobile = 1, reaction%immobile%nimmobile
          if (constraint%immobile_species%external_dataset(iimmobile)) then
            ! no need to requilibrate at each cell
            string = 'constraint ' // trim(constraint%name)
            dataset => DatasetBaseGetPointer(realization%datasets, &
                constraint%immobile_species%constraint_aux_string(iimmobile), &
                string,option)
            if (vec1_loc == PETSC_NULL_VEC) then
              ! cannot use field%work_loc as it is used within ConditionCo...
              call VecDuplicate(field%work_loc,vec1_loc,ierr);CHKERRQ(ierr)
            endif
            idof = ONE_INTEGER
            call ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                                 vec1_loc,LOCAL)
            call VecGetArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = grid%nL2G(local_id)
              rt_auxvars(ghosted_id)%immobile(iimmobile) = vec_p(ghosted_id)
            enddo
            call VecRestoreArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
          endif
        enddo
      endif

      if (.not.option%use_isothermal) then
        equilibrate_at_each_cell = PETSC_TRUE
      endif

      if (use_aq_dataset) then
        call VecGetArrayF90(field%tran_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
        call PetscTime(tstart,ierr);CHKERRQ(ierr)
      endif

      ave_num_iterations = 0.d0
      prev_equilibrated_ghosted_id = 0
      do icell=1,initial_condition%region%num_cells
        local_id = initial_condition%region%cell_ids(icell)
        ghosted_id = grid%nL2G(local_id)
        iend = local_id*option%ntrandof
        ibegin = iend-option%ntrandof+1
        if (patch%imat(ghosted_id) <= 0) then
          xx_p(ibegin:iend) = 1.d-200
          cycle
        endif
        if (equilibrate_at_each_cell) then
          if (use_aq_dataset) then
            offset = (ghosted_id-1)*option%ntrandof
            do iaqdataset = 1, num_aq_datasets
              ! remember that xx_loc_p holds the data set values that
              ! were read in
              temp_int = aq_dataset_to_idof(iaqdataset)
              constraint%aqueous_species%constraint_conc(temp_int) = &
                xx_loc_p(offset+temp_int)
            enddo
          endif
          option%iflag = grid%nG2A(grid%nL2G(local_id))
          if (prev_equilibrated_ghosted_id > 0) then
            ! copy molalities from previous equilibrated auxvar as initial guess
            call RTAuxVarCopyInitialGuess( &
                         rt_auxvars(prev_equilibrated_ghosted_id), &
                         rt_auxvars(ghosted_id),option)
          endif
          call ReactionEquilibrateConstraint(rt_auxvars(ghosted_id), &
            global_auxvars(ghosted_id),material_auxvars(ghosted_id), &
            reaction,constraint, &
            constraint_coupler%num_iterations, &
            (prev_equilibrated_ghosted_id > 0),option)
          option%iflag = 0
          ave_num_iterations = ave_num_iterations + &
            constraint_coupler%num_iterations
          ! update CO2 mole fraction for CO2 modes
#if 0
          ! TODO(geh): ideally, the intermingling of the flow process model
          ! with transport is not ideal.  Peter should be looking into whether
          ! we can remove this code in favor of a slighly less accurate
          ! solution.
          select case(option%iflowmode)
            case(MPH_MODE)
              if (global_auxvars(ghosted_id)%istate == 1) then
                tempreal = &
                  RCO2MoleFraction(rt_auxvars(ghosted_id), &
                                   global_auxvars(ghosted_id),reaction,option)
                ! concentration dof in flow solution vector
                flow_xx_p(local_id*option%nflowdof) = tempreal
              endif
          end select
#endif
          ! prev_eq_ghosted_id is only updated to the prev active cell
          prev_equilibrated_ghosted_id = ghosted_id
        endif
        ! ibegin is the local non-ghosted offset: (local_id-1)*option%ntrandof+1
        offset = ibegin + reaction%offset_aqueous - 1
        ! primary aqueous concentrations
        do idof = 1, reaction%naqcomp
          xx_p(offset+idof) = &
            constraint%aqueous_species%basis_molarity(idof) / &
            global_auxvars(ghosted_id)%den_kg(iphase)*1000.d0 ! convert molarity -> molality
        enddo
        ! mineral volume fractions
        if (associated(constraint%minerals)) then
          do imnrl = 1, reaction%mineral%nkinmnrl
            ! if read from a dataset, the vol frac was set above.  Don't want to
            ! overwrite
            if (.not.constraint%minerals% &
                  external_vol_frac_dataset(imnrl)) then
              rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl) = &
                constraint%minerals%constraint_vol_frac(imnrl)
              rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl) = &
                constraint%minerals%constraint_vol_frac(imnrl)
            endif
            if (.not.constraint%minerals% &
                  external_area_dataset(imnrl)) then
              rt_auxvars(ghosted_id)%mnrl_area0(imnrl) = &
                constraint%minerals%constraint_area(imnrl)
              rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
                constraint%minerals%constraint_area(imnrl)
            endif
          enddo
        endif
        ! kinetic surface complexes
        if (associated(constraint%surface_complexes)) then
          do idof = 1, reaction%surface_complexation%nkinsrfcplx
            rt_auxvars(ghosted_id)%kinsrfcplx_conc(idof,-1) = & !geh: to catch bug
              constraint%surface_complexes%constraint_conc(idof)
          enddo
          do ikinrxn = 1, reaction%surface_complexation%nkinsrfcplxrxn
            irxn = reaction%surface_complexation%kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
            isite = reaction%surface_complexation%srfcplxrxn_to_surf(irxn)
            rt_auxvars(ghosted_id)%kinsrfcplx_free_site_conc(isite) = &
              constraint%surface_complexes%basis_free_site_conc(isite)
          enddo
        endif
        ! this is for the multi-rate surface complexation model
        if (reaction%surface_complexation%nkinmrsrfcplxrxn > 0 .and. &
          ! geh: if we re-equilibrate at each grid cell, we do not want to
          ! overwrite the reequilibrated values with those from the constraint
            .not. equilibrate_at_each_cell) then
          ! copy over total sorbed concentration
          rt_auxvars(ghosted_id)%kinmr_total_sorb = &
            constraint_coupler%rt_auxvar%kinmr_total_sorb
          ! copy over free site concentration
          rt_auxvars(ghosted_id)%srfcplxrxn_free_site_conc = &
            constraint_coupler%rt_auxvar%srfcplxrxn_free_site_conc
        endif
        ! immobile
        if (associated(constraint%immobile_species)) then
          offset = ibegin + reaction%offset_immobile - 1
          do iimmobile = 1, reaction%immobile%nimmobile
            if (constraint%immobile_species%external_dataset(iimmobile)) then
              ! already read into rt_auxvars above.
              xx_p(offset+iimmobile) = &
                rt_auxvars(ghosted_id)%immobile(iimmobile)
            else
              xx_p(offset+iimmobile) = &
                constraint%immobile_species%constraint_conc(iimmobile)
              rt_auxvars(ghosted_id)%immobile(iimmobile) = &
                constraint%immobile_species%constraint_conc(iimmobile)
            endif
          enddo
        endif
        if (option%use_sc) then
          reaction%mc_flag = 1
          do cell = 1, rt_sec_transport_vars(ghosted_id)%ncells
            call ReactionEquilibrateConstraint( &
                                     rt_sec_transport_vars(ghosted_id)% &
                                       sec_rt_auxvar(cell), &
                                     global_auxvars(ghosted_id), &
                                     material_auxvars(ghosted_id),reaction, &
                                     sec_tran_constraint, &
                                     sec_constraint_coupler%num_iterations, &
                                     PETSC_FALSE,option)
            if (option%transport%sc_fixed_water_density) then
              ! convert molality to molarity so that secondary continuum is
              ! solely molarity (independent of fluid density)
              rt_sec_transport_vars(ghosted_id)% &
                      sec_rt_auxvar(cell)%pri_molal = &
                global_auxvars(ghosted_id)%den_kg(1) * 1.d-3 * &
                rt_sec_transport_vars(ghosted_id)%sec_rt_auxvar(cell)%pri_molal
            endif
            rt_sec_transport_vars(ghosted_id)%updated_conc(:,cell) =  &
              rt_sec_transport_vars(ghosted_id)%sec_rt_auxvar(cell)%pri_molal
          enddo
          reaction%mc_flag = 0
        endif
      enddo ! icell=1,initial_condition%region%num_cells
      if (use_aq_dataset) then
        call PetscTime(tend,ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(field%tran_xx_loc,xx_loc_p, &
                                ierr);CHKERRQ(ierr)
        ave_num_iterations = ave_num_iterations / &
          initial_condition%region%num_cells
        write(option%io_buffer,&
              '("Average number of iterations in ReactionEquilibrateConstraint():", &
              & f5.1)') ave_num_iterations
        call PrintMsg(option)
        write(option%io_buffer,'(f10.2," Seconds to equilibrate constraints")') &
          tend-tstart
        call PrintMsg(option)
      endif
      initial_condition => initial_condition%next
    enddo

    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
    select case(option%iflowmode)
      case(MPH_MODE)
        call VecRestoreArrayF90(field%flow_xx,flow_xx_p,ierr);CHKERRQ(ierr)
    end select

  ! check to ensure that minimum concentration is not less than or equal
  ! to zero
  call VecMin(field%tran_xx,PETSC_NULL_INTEGER,tempreal,ierr);CHKERRQ(ierr)
  if (tempreal <= 0.d0) then
    option%io_buffer = 'ERROR: Zero concentrations found in initial ' // &
      'transport solution.'
    call PrintMsg(option)
    ! now figure out which species have zero concentrations
    do idof = 1, option%ntrandof
      call VecStrideMin(field%tran_xx,idof-1,offset,tempreal, &
                        ierr);CHKERRQ(ierr)
      if (tempreal <= 0.d0) then
        write(string,*) tempreal
        if (idof <= reaction%naqcomp) then
          string2 = '  Aqueous species "' // &
            trim(reaction%primary_species_names(idof))
        else
          string2 = '  Immobile species "' // &
            trim(reaction%immobile%names(idof-reaction%offset_immobile))
        endif
          string2 = trim(string2) // &
            '" has zero concentration (' // &
            trim(adjustl(string)) // ').'
        call PrintMsg(option,string2)
      endif
    enddo
    option%io_buffer = ''
    call PrintMsg(option)
    option%io_buffer = '*** Begin Note'
    call PrintMsg(option)
    option%io_buffer = 'If concentrations = -999., they have not ' // &
              'been initialized properly.'
    call PrintMsg(option)
    option%io_buffer = '*** End Note'
    call PrintMsg(option)
    option%io_buffer = 'Free ion concentrations must be positive.  Try ' // &
      'using a small value such as 1.e-20 or 1.e-40 instead of zero.'
    call PrintErrMsg(option)
  endif

  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)
  call VecCopy(field%tran_xx,field%tran_yy,ierr);CHKERRQ(ierr)

  ! override initial conditions if they are to be read from a file
  if (len_trim(option%initialize_transport_filename) > 1) then
    call CondControlReadTransportIC(realization, &
                                    option%initialize_transport_filename)
  endif
  call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)
  ! at this point the auxvars have been computed with activity coef = 1.d0
  ! to use intitial condition with activity coefs /= 1.d0, must update
  ! activity coefs and recompute auxvars
  if (realization%reaction%act_coef_update_frequency /= &
      ACT_COEF_FREQUENCY_OFF) then
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_TRUE)
    !geh: you may ask, why call this twice....  We need to iterate at least
    !     once to ensure that the activity coefficients are more accurate.
    !     Otherwise, the total component concentrations can be quite
    !     different from what is defined in the input file.
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_TRUE)
  endif

  if (vec1_loc /= PETSC_NULL_VEC) then
    call VecDestroy(vec1_loc,ierr);CHKERRQ(ierr)
  endif
  if (vec2_loc /= PETSC_NULL_VEC) then
    call VecDestroy(vec2_loc,ierr);CHKERRQ(ierr)
  endif

end subroutine CondControlAssignRTTranInitCond

! ************************************************************************** !

subroutine CondControlAssignNWTranInitCond(realization)
  !
  ! Assigns transport initial conditions to model, and equilibrates the
  ! initial conditions according to the constraint types.
  !
  ! Author: Jenn Frederick
  ! Date: 04/02/2019
  !

  use Realization_Subsurface_class
  use Discretization_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Transport_Constraint_NWT_module
  use Grid_module
  use Global_Aux_module
  use Patch_module
  use NW_Transport_module
  use NW_Transport_Aux_module
  use NWT_Equilibrium_module
  use Material_Aux_module
  use HDF5_module

  implicit none

  class(realization_subsurface_type) :: realization

  PetscInt :: icell, idof
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscReal, pointer :: xx_p(:)
  Vec :: vec1_loc
  Vec :: vec2_loc
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(patch_type), pointer :: patch
  class(reaction_nw_type), pointer :: reaction_nw
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  class(tran_constraint_coupler_nwt_type), pointer :: constraint_coupler
  class(tran_constraint_nwt_type), pointer :: constraint

  PetscInt :: iphase
  PetscInt :: offset
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscReal :: tempreal

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  reaction_nw => realization%reaction_nw
  patch => realization%patch
  grid => patch%grid

  iphase = 1
  vec1_loc = PETSC_NULL_VEC
  vec2_loc = PETSC_NULL_VEC

  !TODO(jenn) Do not allow MPH_MODE with NW Transport.

  ! TODO(patch_list) fix indentation

    material_auxvars => patch%aux%Material%auxvars
    global_auxvars => patch%aux%Global%auxvars
    nwt_auxvars => patch%aux%NWT%auxvars

    call VecGetArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
    xx_p = UNINITIALIZED_DOUBLE

    initial_condition => patch%initial_condition_list%first
    do
      if (.not.associated(initial_condition)) exit

      constraint_coupler => &
        TranConstraintCouplerNWTCast(initial_condition%tran_condition% &
                                       cur_constraint_coupler)
      constraint => TranConstraintNWTCast(constraint_coupler%constraint)

      do icell=1,initial_condition%region%num_cells

        local_id = initial_condition%region%cell_ids(icell)
        ghosted_id = grid%nL2G(local_id)

        iend = local_id*option%ntrandof
        ibegin = iend-option%ntrandof+1

        if (patch%imat(ghosted_id) <= 0) then
          xx_p(ibegin:iend) = 1.d-200
          cycle
        endif

        call NWTEquilibrateConstraint(reaction_nw,constraint, &
                                      nwt_auxvars(ghosted_id), &
                                      global_auxvars(ghosted_id), &
                                      material_auxvars(ghosted_id), &
                                      option)


        ! ibegin is the local non-ghosted offset: (local_id-1)*option%ntrandof+1
        offset = ibegin - 1

        ! species concentrations
        do idof = 1, reaction_nw%params%nspecies
          xx_p(offset+idof) = nwt_auxvars(ghosted_id)%total_bulk_conc(idof)
        enddo

      enddo ! icell=1,initial_condition%region%num_cells
      initial_condition => initial_condition%next
    enddo

    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)

  ! check to ensure that minimum concentration is not less than or equal
  ! to zero
  call VecMin(field%tran_xx,PETSC_NULL_INTEGER,tempreal,ierr);CHKERRQ(ierr)
  if (tempreal <= 0.d0) then
    option%io_buffer = 'ERROR: Zero concentrations found in initial &
                       &transport solution.'
    call PrintMsg(option)
    ! now figure out which species have zero concentrations
    do idof = 1, option%ntrandof
      call VecStrideMin(field%tran_xx,idof-1,offset,tempreal, &
                        ierr);CHKERRQ(ierr)
      if (tempreal <= 0.d0) then
        write(string,*) tempreal
        string2 = '  Species "' // trim(reaction_nw%species_names(idof))
        string2 = trim(string2) // '" has zero concentration (' // &
                  trim(adjustl(string)) // ').'
        call PrintMsg(option,string2)
      endif
    enddo
    option%io_buffer = ''
    call PrintMsg(option)
    option%io_buffer = '*** Begin Note'
    call PrintMsg(option)
    option%io_buffer = 'If concentrations = -999., they have not ' // &
              'been initialized properly.'
    call PrintMsg(option)
    option%io_buffer = '*** End Note'
    call PrintMsg(option)
    option%io_buffer = 'Species concentrations must be positive.  Try ' // &
      'using a small value such as 1.e-20 or 1.e-40 instead of zero.'
    call PrintErrMsg(option)
  endif

  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)
  call VecCopy(field%tran_xx,field%tran_yy,ierr);CHKERRQ(ierr)

  ! override initial conditions if they are to be read from a file
  if (len_trim(option%initialize_transport_filename) > 1) then
    call CondControlReadTransportIC(realization, &
                                    option%initialize_transport_filename)
  endif

  call NWTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE)

  if (vec1_loc /= PETSC_NULL_VEC) then
    call VecDestroy(vec1_loc,ierr);CHKERRQ(ierr)
  endif
  if (vec2_loc /= PETSC_NULL_VEC) then
    call VecDestroy(vec2_loc,ierr);CHKERRQ(ierr)
  endif

end subroutine CondControlAssignNWTranInitCond

! ************************************************************************** !

subroutine ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                           mdof_vec,vec_type)
  !
  ! maps an external dataset to a PETSc vec
  ! representing values at each grid cell
  !
  ! Author: Glenn Hammond
  ! Date: 03/23/12
  !
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Dataset_Common_HDF5_class
  use Dataset_Base_class
  use HDF5_module
  use Discretization_module

  implicit none


  class(realization_subsurface_type) :: realization
  class(dataset_base_type), pointer :: dataset
  PetscInt :: idof
  Vec :: mdof_vec
  PetscInt :: vec_type

  type(field_type), pointer :: field
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  field => realization%field
  option => realization%option

  call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
  if (associated(dataset)) then
    select type(dataset)
      class is (dataset_common_hdf5_type)
        string = '' ! group name
        ! have to copy to string2 due to mismatch in string size
        string2 = dataset%hdf5_dataset_name
        call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                          dataset%filename, &
                                          string,string2, &
                                          dataset%realization_dependent)
        if (vec_type == GLOBAL) then
          call VecStrideScatter(field%work,idof-1,mdof_vec,INSERT_VALUES, &
                                ierr);CHKERRQ(ierr)
        else
          call DiscretizationGlobalToLocal(realization%discretization, &
                                           field%work, &
                                           field%work_loc,ONEDOF)
          call VecStrideScatter(field%work_loc,idof-1,mdof_vec,INSERT_VALUES, &
                                ierr);CHKERRQ(ierr)
        endif
      class default
        option%io_buffer = 'Dataset "' // trim(dataset%name) // &
          '" not supported in ConditionControlMapDatasetToVec.'
        call PrintErrMsg(option)
    end select
  endif

end subroutine ConditionControlMapDatasetToVec

! ************************************************************************** !

subroutine CondControlScaleSourceSink(realization)
  !
  ! Scales select source/sinks based on perms
  !
  ! Author: Glenn Hammond
  ! Date: 09/03/08, 10/18/11
  !
#include "petsc/finclude/petscdmda.h"
  use petscdmda

  use Realization_Subsurface_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  use Patch_module
  use Material_Aux_module
  use Variables_module, only : PERMEABILITY_X

  implicit none

  class(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: cur_source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id, neighbor_ghosted_id
  PetscInt :: iconn
  PetscReal :: scale, sum
  PetscInt :: icount
  PetscInt :: x_count, y_count, z_count
  PetscInt, parameter :: x_width = 1, y_width = 1, z_width = 0

  PetscInt :: ghosted_neighbors(0:27)

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  material_auxvars => realization%patch%aux%Material%auxvars

  select case(option%iflowmode)
    case(TH_MODE,TH_TS_MODE,MPH_MODE,PNF_MODE)
      option%io_buffer = 'Flow mode ' // trim(option%flowmode) // ' not &
        &supported in CondControlScaleSourceSink().'
      call PrintErrMsg(option)
  end select

  ! TODO(patch_list) fix indentation

    cur_source_sink => patch%source_sink_list%first
    do
      if (.not.associated(cur_source_sink)) exit

      call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

      cur_connection_set => cur_source_sink%connection_set

      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)

        select case(option%iflowmode)
          case(RICHARDS_MODE,RICHARDS_TS_MODE,G_MODE,WF_MODE,H_MODE,&
               ZFLOW_MODE)
              call GridGetGhostedNeighbors(grid,ghosted_id,DMDA_STENCIL_STAR, &
                                          x_width,y_width,z_width, &
                                          x_count,y_count,z_count, &
                                          ghosted_neighbors,option)
              ! ghosted neighbors is ordered first in x, then, y, then z
              icount = 0
              sum = 0.d0
              ! x-direction
              do while (icount < x_count)
                icount = icount + 1
                neighbor_ghosted_id = ghosted_neighbors(icount)
                sum = sum + MaterialAuxVarGetValue(material_auxvars( &
                              neighbor_ghosted_id),PERMEABILITY_X) * &
                            grid%structured_grid%dy(neighbor_ghosted_id)* &
                            grid%structured_grid%dz(neighbor_ghosted_id)

              enddo
              ! y-direction
              do while (icount < x_count + y_count)
                icount = icount + 1
                neighbor_ghosted_id = ghosted_neighbors(icount)
                sum = sum + MaterialAuxVarGetValue(material_auxvars( &
                              neighbor_ghosted_id),PERMEABILITY_X) * &
                            grid%structured_grid%dx(neighbor_ghosted_id)* &
                            grid%structured_grid%dz(neighbor_ghosted_id)

              enddo
              ! z-direction
              do while (icount < x_count + y_count + z_count)
                icount = icount + 1
                neighbor_ghosted_id = ghosted_neighbors(icount)
                sum = sum + MaterialAuxVarGetValue(material_auxvars( &
                              neighbor_ghosted_id),PERMEABILITY_X) * &
                            grid%structured_grid%dx(neighbor_ghosted_id)* &
                            grid%structured_grid%dy(neighbor_ghosted_id)
              enddo
              vec_ptr(local_id) = vec_ptr(local_id) + sum
          case(TH_MODE,TH_TS_MODE)
          case(MPH_MODE)
        end select

      enddo

      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      call VecNorm(field%work,NORM_1,scale,ierr);CHKERRQ(ierr)
      scale = 1.d0/scale
      call VecScale(field%work,scale,ierr);CHKERRQ(ierr)

      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        select case(option%iflowmode)
          case(RICHARDS_MODE,RICHARDS_TS_MODE,G_MODE,WF_MODE,H_MODE, &
               ZFLOW_MODE)
            cur_source_sink%flow_aux_real_var(ONE_INTEGER,iconn) = &
              vec_ptr(local_id)
        end select

      enddo
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

      cur_source_sink => cur_source_sink%next
    enddo

end subroutine CondControlScaleSourceSink

! ************************************************************************** !

subroutine CondControlReadTransportIC(realization,filename)
  !
  ! Assigns transport initial condition from
  ! HDF5 file
  !
  ! Author: Glenn Hammond
  ! Date: 03/05/10
  !

  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
  use Reactive_Transport_module
  use Reaction_Aux_module
  use Discretization_module
  use HDF5_module

  implicit none

  class(realization_subsurface_type) :: realization
  character(len=MAXSTRINGLENGTH) :: filename

  PetscInt :: local_id, idx, offset, idof
  PetscReal, pointer :: xx_p(:)
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscReal, pointer :: vec_p(:)
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  class(reaction_rt_type), pointer :: reaction

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  reaction => realization%reaction

  ! TODO(patch_list) fix indentation

      ! assign initial conditions values to domain
    call VecGetArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)

    ! Primary species concentrations for all modes
    do idof = 1, option%ntrandof ! primary aqueous concentrations
      offset = idof
      group_name = ''
      if (associated(reaction)) &
        dataset_name = reaction%primary_species_names(idof)
      if (associated(realization%reaction_nw)) &
        dataset_name = realization%reaction_nw%species_names(idof)
      call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                        filename,group_name, &
                                        dataset_name,option%id>0)
      call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
      do local_id=1, grid%nlmax
        if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
        if (vec_p(local_id) < 1.d-40) then
          print *,  option%myrank, grid%nG2A(grid%nL2G(local_id)), &
            ': Zero free-ion concentration in Initial Condition read from file.'
        endif
        idx = (local_id-1)*option%ntrandof + offset
        xx_p(idx) = vec_p(local_id)
      enddo
      call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)

    enddo

    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)

  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)
  call VecCopy(field%tran_xx,field%tran_yy,ierr);CHKERRQ(ierr)

end subroutine CondControlReadTransportIC

end module Condition_Control_module
