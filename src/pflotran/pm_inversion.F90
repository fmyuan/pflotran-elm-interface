module PM_Inversion_class

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_class
  use Option_Inversion_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_inversion_type
    class(realization_subsurface_type), pointer :: realization
    character(len=MAXWORDLENGTH) :: ctype
    procedure(PMInversionEvaluate), pointer :: Evaluate => null()
  contains
    procedure, public :: Setup => PMInversionSetup
    procedure, public :: InitializeRun => PMInversionInitializeRun
    procedure, public :: FinalizeRun => PMInversionFinalizeRun
    procedure, public :: Destroy => PMInversionDestroy
  end type pm_inversion_type

  ! interface blocks
  interface
    subroutine PMInversionEvaluate(this,time,ierr)
      import :: pm_inversion_type
      implicit none
      class(pm_inversion_type) :: this
      PetscReal :: time
      PetscErrorCode :: ierr
    end subroutine PMInversionEvaluate
  end interface

  public :: PMInversionCreate, &
            PMInversionInit, &
            PMInversionSetFunctionPointer

contains

! ************************************************************************** !

function PMInversionCreate()
  !
  ! Creates an inversion process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22
  !

  implicit none

  class(pm_inversion_type), pointer :: PMInversionCreate

  class(pm_inversion_type), pointer :: pm

  allocate(pm)
  call PMInversionInit(pm)

  PMInversionCreate => pm

end function PMInversionCreate

! ************************************************************************** !

subroutine PMInversionInit(this)
  !
  ! Initializes inversion process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  implicit none

  class(pm_inversion_type) :: this

  nullify(this%realization)
  this%ctype = ''
  this%name = ''

  call PMBaseInit(this)
  ! restart not currently supported for inversion pm's, and not needed.
  this%skip_restart = PETSC_TRUE

end subroutine PMInversionInit

! ************************************************************************** !

subroutine PMInversionSetup(this)
  !
  ! Dummy setup routine
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  implicit none

  class(pm_inversion_type) :: this

end subroutine PMInversionSetup

! ************************************************************************** !

subroutine PMInversionSetFunctionPointer(this,string)
  !
  ! Initializes inversion process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  use Option_module

  implicit none

  class(pm_inversion_type) :: this
  character(len=*) :: string

  this%ctype = trim(string)
  select case(string)
    case('INVERSION_MEASUREMENT')
      this%Evaluate => PMInversionInversionMeasurement
      this%header = 'INVERSION MEASUREMENT'
      this%name = 'Inversion Measurement'
    case('INVERSION_ADJOINT')
      this%Evaluate => PMInversionInversionAdjoint
      this%header = 'INVERSION ADJOINT'
      this%name = 'Inversion Adjoint'
    case default
      this%option%io_buffer = 'Function pointer "' // trim(string) // '" not &
        &found among available functions in PMInversionSetFunctionPointer.'
      call PrintErrMsg(this%option)
  end select

end subroutine PMInversionSetFunctionPointer

! ************************************************************************** !

recursive subroutine PMInversionInitializeRun(this)
  !
  ! Initializes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  use Condition_module
  use Coupler_module
  use Option_module

  implicit none

  class(pm_inversion_type) :: this

end subroutine PMInversionInitializeRun

! ************************************************************************** !

subroutine PMInversionInversionMeasurement(this,time,ierr)
  !
  ! Takes measurements for inversion calculation
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  use Inversion_Aux_module
  use Inversion_Coupled_Aux_module
  use Inversion_Measurement_Aux_module
  use Option_module
  use Realization_Base_class
  use Realization_Subsurface_class
  use String_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_SATURATION, &
                               SOLUTE_CONCENTRATION, DERIVATIVE, POROSITY
  use ZFlow_Aux_module
  use Utility_module

  implicit none

  class(pm_inversion_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(inversion_aux_type), pointer :: inversion_aux
  type(inversion_measurement_aux_type), pointer :: measurements(:)
  type(inversion_coupled_soln_type), pointer :: solutions(:)
  PetscInt :: imeasurement
  PetscInt :: iert_measurement
  PetscInt :: icount
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscReal :: tempreal
  PetscBool :: iflag
  PetscBool :: no_measurement_flag

  option => this%option
  inversion_aux => this%realization%patch%aux%inversion_aux
  nullify(measurements)

  ierr = 0
  no_measurement_flag = PETSC_TRUE
  if (associated(inversion_aux)) then
    if (option%inversion%record_measurements .and. &
        inversion_aux%isync_time <= &
        size(inversion_aux%sync_times)) then
      measurements => inversion_aux%measurements
      if (Equal(inversion_aux%sync_times( &
                  inversion_aux%isync_time),time)) then
        inversion_aux%isync_time = inversion_aux%isync_time + 1
        icount = 0
        do imeasurement = 1, size(measurements)
          local_id = measurements(imeasurement)%local_id
          if (Initialized(local_id)) then
            icount = icount + 1
            if (.not.measurements(imeasurement)%measured .and. &
                Equal(measurements(imeasurement)%time,time)) then
              tempreal = UNINITIALIZED_DOUBLE
              ghosted_id = this%realization%patch%grid%nL2G(local_id)
              select case(measurements(imeasurement)%iobs_var)
                case(OBS_ERT_MEASUREMENT)
                  iert_measurement = measurements(imeasurement)%cell_id
                  tempreal = this%realization%survey%dsim(iert_measurement)
                case default
                  select case(measurements(imeasurement)%iobs_var)
                    case(OBS_LIQUID_PRESSURE)
                      ivar = LIQUID_PRESSURE
                    case(OBS_LIQUID_SATURATION)
                      ivar = LIQUID_SATURATION
                    case(OBS_SOLUTE_CONCENTRATION)
                      ivar = SOLUTE_CONCENTRATION
                    case default
                      option%io_buffer = 'Unrecognized observed variable in &
                        &PMInversionInversionMeasurement: ' // &
                        trim(StringWrite(measurements(imeasurement)%iobs_var))
                      call PrintErrMsgByRank(option)
                  end select
                  if (Initialized(inversion_aux% &
                                    local_measurement_values_ptr(icount))) then
                    option%io_buffer = 'Duplicate measurement in &
                      &PMInversionInversionMeasurement for measurement: ' // &
                      trim(StringWrite(imeasurement))
                    call PrintErrMsgByRank(option)
                  endif
                  tempreal = &
                    RealizGetVariableValueAtCell(this%realization,ghosted_id, &
                                                ivar,ZERO_INTEGER)
              end select
              no_measurement_flag = PETSC_FALSE
              inversion_aux%local_measurement_values_ptr(icount) = &
                tempreal
              option%io_buffer = &
                InvMeasAnnounceToString(measurements(imeasurement), &
                                        tempreal,option)
              call PrintMsgByRank(option)
              if (inversion_aux%store_adjoint) then
                ! store derivative of observation wrt unknown
                if (associated(inversion_aux% &
                                local_dobs_dunknown_values_ptr)) then
                  isubvar = UNINITIALIZED_INTEGER
                  select case(measurements(imeasurement)%iobs_var)
                    case(OBS_LIQUID_SATURATION)
                      isubvar = ZFLOW_LIQ_SAT_WRT_LIQ_PRES
                  end select
                  if (Initialized(isubvar)) then
                    tempreal = &
                      RealizGetVariableValueAtCell(this%realization, &
                            ghosted_id,DERIVATIVE,isubvar)
                    inversion_aux% &
                      local_dobs_dunknown_values_ptr(icount) = tempreal
                  endif
                endif
                ! store derivative of observation wrt parameter
                if (associated(inversion_aux% &
                                local_dobs_dparam_values_ptr)) then
                  isubvar = UNINITIALIZED_INTEGER
                  select case(measurements(imeasurement)%iobs_var)
                    case(OBS_LIQUID_PRESSURE)
                      select case(inversion_aux%parameters(1)%itype)
                        case(POROSITY)
!                          isubvar = ZFLOW_LIQ_PRES_WRT_POROS
                      end select
                  end select
                  if (Initialized(isubvar)) then
                    tempreal = &
                      RealizGetVariableValueAtCell(this%realization, &
                            ghosted_id,DERIVATIVE,isubvar)
                    inversion_aux% &
                      local_dobs_dparam_values_ptr(icount) = tempreal
                  endif
                endif
              endif ! store adjoint
            endif ! not measured
          endif ! initialized
        enddo
      endif
    endif
    if (option%inversion%coupled_flow_ert) then
      solutions => inversion_aux%coupled_aux%solutions
      do icount = 1, size(solutions)
        if (.not.solutions(icount)%measured .and. &
            Equal(solutions(icount)%time,time)) then
          call RealizationGetVariable(this%realization, &
                                      solutions(icount)% &
                                        perturbed_saturation_solution, &
                                      LIQUID_SATURATION,ZERO_INTEGER)
          iflag = PETSC_FALSE
          no_measurement_flag = PETSC_FALSE
          if (solutions(icount)%perturbed_solute_solution /= &
              PETSC_NULL_VEC) then
            iflag = PETSC_TRUE
            call RealizationGetVariable(this%realization, &
                                        solutions(icount)% &
                                          perturbed_solute_solution, &
                                        SOLUTE_CONCENTRATION,ZERO_INTEGER)
          endif
          option%io_buffer = 'Full saturation'
          if (iflag) then
            option%io_buffer = trim(option%io_buffer) // &
              ' and solute concentration'
          endif
          option%io_buffer = trim(option%io_buffer) // &
              ' solutions measured.'
          call PrintMsg(option)
        endif
      enddo
    endif
    if (no_measurement_flag) then
      call PrintMsg(option,'  No measurement requested at this time.')
    endif
  endif

end subroutine PMInversionInversionMeasurement

! ************************************************************************** !

subroutine PMInversionInversionAdjoint(this,time,ierr)
  !
  ! Initializes inversion process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  use Inversion_Aux_module
  use Realization_Subsurface_class

  implicit none

  class(pm_inversion_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  ierr = 0
  if (associated(this%realization%patch%aux%inversion_aux)) then
    call InversionAuxAdjointRecordTS(this%realization%patch%aux% &
                                       inversion_aux,time)
  endif

end subroutine PMInversionInversionAdjoint

! ************************************************************************** !

recursive subroutine PMInversionFinalizeRun(this)
  !
  ! Finalizes the run
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22
  !

  implicit none

  class(pm_inversion_type) :: this

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMInversionFinalizeRun

! ************************************************************************** !

subroutine PMInversionDestroy(this)
  !
  ! Destroys inversion process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  implicit none

  class(pm_inversion_type) :: this

  call PMBaseDestroy(this)

  ! destroyed in realization
  nullify(this%realization)
  nullify(this%option)

end subroutine PMInversionDestroy

end module PM_Inversion_class
