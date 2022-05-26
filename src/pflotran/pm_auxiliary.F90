module PM_Auxiliary_class

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_auxiliary_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    character(len=MAXWORDLENGTH) :: ctype
    type(pm_auxiliary_salinity_type), pointer :: salinity
    PetscBool :: evaluate_at_end_of_simulation
    procedure(PMAuxliaryEvaluate), pointer :: Evaluate => null()
  contains
    procedure, public :: Setup => PMAuxiliarySetup
    procedure, public :: InitializeRun => PMAuxiliaryInitializeRun
    procedure, public :: FinalizeRun => PMAuxiliaryFinalizeRun
    procedure, public :: InputRecord => PMAuxiliaryInputRecord
    procedure, public :: Destroy => PMAuxiliaryDestroy
  end type pm_auxiliary_type

  type :: pm_auxiliary_salinity_type
    PetscInt :: nspecies
    character(len=MAXWORDLENGTH) :: species_names(6)
    PetscInt :: ispecies(6)
    PetscReal :: molecular_weights(6)
  end type pm_auxiliary_salinity_type

  ! interface blocks
  interface
    subroutine PMAuxliaryEvaluate(this,time,ierr)
      import :: pm_auxiliary_type
      implicit none
      class(pm_auxiliary_type) :: this
      PetscReal :: time
      PetscErrorCode :: ierr
    end subroutine PMAuxliaryEvaluate
  end interface

  public :: PMAuxiliaryCreate, &
            PMAuxiliaryInit, &
            PMAuxiliaryCast, &
            PMAuxiliaryRead, &
            PMAuxiliaryInputRecord, &
            PMAuxiliarySetFunctionPointer

contains

! ************************************************************************** !

function PMAuxiliaryCreate()
  !
  ! Creates reactive transport process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  implicit none

  class(pm_auxiliary_type), pointer :: PMAuxiliaryCreate

  class(pm_auxiliary_type), pointer :: pm

  allocate(pm)
  call PMAuxiliaryInit(pm)

  PMAuxiliaryCreate => pm

end function PMAuxiliaryCreate

! ************************************************************************** !

subroutine PMAuxiliaryInit(this)
  !
  ! Initializes auxiliary process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none

  class(pm_auxiliary_type) :: this

  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%salinity)
  this%ctype = ''
  this%name = ''
  this%evaluate_at_end_of_simulation = PETSC_TRUE

  call PMBaseInit(this)
  ! restart not currently supported for auxiliary pm's, and not needed.
  this%skip_restart = PETSC_TRUE

end subroutine PMAuxiliaryInit

! ************************************************************************** !

subroutine PMAuxiliarySetup(this)
  !
  ! Sets up auxiliary process model
  !
  ! Author: Glenn Hammond
  ! Date: 03/02/16

  implicit none

  class(pm_auxiliary_type) :: this

end subroutine PMAuxiliarySetup

! ************************************************************************** !

function PMAuxiliaryCast(this)
  !
  ! Initializes auxiliary process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none

  class(pm_base_type), pointer :: this

  class(pm_auxiliary_type), pointer :: PMAuxiliaryCast

  nullify(PMAuxiliaryCast)
  if (.not.associated(this)) return
  select type (this)
    class is (pm_auxiliary_type)
      PMAuxiliaryCast => this
    class default
      !geh: have default here to pass a null pointer if not of type ascii
  end select

end function PMAuxiliaryCast

! ************************************************************************** !

subroutine PMAuxiliaryRead(input, option, this)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_auxiliary_type), pointer :: this

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: i
  PetscReal :: tempreal

  error_string = 'SIMULATION,PROCESS_MODELS,AUXILIARY'
  call InputReadCard(input,option,word,PETSC_FALSE)
  call InputErrorMsg(input,option,'type',error_string)
  call StringToUpper(word)
  error_string = trim(error_string) // ',' // trim(word)

  this%ctype = word
  select case(word)
    case('SALINITY')
      option%flow%density_depends_on_salinity = PETSC_TRUE
      allocate(this%salinity)
      this%salinity%nspecies = 0
      this%salinity%species_names = ''
      this%salinity%ispecies = UNINITIALIZED_INTEGER
      this%salinity%molecular_weights = UNINITIALIZED_DOUBLE
      i = 0
      word = ''
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(word)
          case('SPECIES')
            i = i + 1
            if (i > 6) then
              option%io_buffer = 'Number of salinity species exceeded.'
              call PrintErrMsgToDev(option, &
                                    'ask for the maximum number of salinity &
                                    &species to be increased')
            endif
            call InputReadWord(input,option,this%salinity% &
                                 species_names(i),PETSC_TRUE)
            call InputErrorMsg(input,option,'species_name',error_string)
            call InputReadDouble(input,option,tempreal)
            if (input%ierr == 0) then
              this%salinity%molecular_weights(i) = tempreal
            else
              ! for now let's print an error message.  Decide on whether to
              ! read from database later.
              call InputErrorMsg(input,option,'molecular weight',error_string)
            endif
          case default
            error_string = trim(error_string) // 'SALINITY'
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
      enddo
      call InputPopBlock(input,option)
      this%salinity%nspecies = i
    case default
      call InputKeywordUnrecognized(input,word,error_string,option)
  end select

  call PMAuxiliarySetFunctionPointer(this,this%ctype)

end subroutine PMAuxiliaryRead

! ************************************************************************** !

subroutine PMAuxiliarySetFunctionPointer(this,string)
  !
  ! Initializes auxiliary process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  use Option_module

  implicit none

  class(pm_auxiliary_type) :: this
  character(len=*) :: string

  this%ctype = trim(string)
  select case(string)
    case('EVOLVING_STRATA')
      this%Evaluate => PMAuxiliaryEvolvingStrata
      this%name = 'auxiliary evolving strata'
      this%header = 'AUXILIARY EVOLVING STRATA'
      this%evaluate_at_end_of_simulation = PETSC_FALSE
    case('INVERSION_MEASUREMENT')
      this%Evaluate => PMAuxiliaryInversionMeasurement
      this%header = 'AUXILIARY INVERSION MEASUREMENT'
    case('INVERSION_ADJOINT')
      this%Evaluate => PMAuxiliaryInversionAdjoint
      this%header = 'AUXILIARY INVERSION ADJOINT'
    case('SALINITY')
      this%Evaluate => PMAuxiliarySalinity
      this%header = 'AUXILIARY SALINITY'
    case default
      this%option%io_buffer = 'Function pointer "' // trim(string) // '" not &
        &found among available functions in PMAuxiliarySetFunctionPointer.'
      call PrintErrMsg(this%option)
  end select

end subroutine PMAuxiliarySetFunctionPointer

! ************************************************************************** !

recursive subroutine PMAuxiliaryInitializeRun(this)
  !
  ! Initializes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  use Condition_module
  use Coupler_module
  use Option_module
  use Reaction_Aux_module

  implicit none

  class(pm_auxiliary_type) :: this

  type(coupler_type), pointer :: boundary_condition
  PetscReal :: time
  PetscInt :: i
  PetscErrorCode :: ierr

  ierr = 0
  time = 0.d0
  select case(this%ctype)
    case('EVOLVING_STRATA')
!      call MatSetOption(Jacobian,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE, &
!                        ierr);CHKERRQ(ierr)
!      call MatSetOption(Jacobian,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE, &
!                        ierr);CHKERRQ(ierr)
    case('INVERSION')
    case('SALINITY')
      ! check if any boundary conditions are hydrostatic, as hydrostatic are
      ! currently not supported
      boundary_condition => &
        this%realization%patch%boundary_condition_list%first
      do
        if (.not.associated(boundary_condition)) exit
        if (FlowConditionIsHydrostatic(boundary_condition%flow_condition)) then
          this%option%io_buffer = 'Hydrostatic flow conditions are currently &
            &not supported by the SALINITY process model.'
          call PrintErrMsg(this%option)
        endif
        boundary_condition => boundary_condition%next
      enddo

      ! set up species names
      do i =1, this%salinity%nspecies
        this%salinity%ispecies(i) = &
          GetPrimarySpeciesIDFromName(this%salinity%species_names(i), &
                                      this%realization%patch%reaction, &
                                      this%option)
        if (Uninitialized(this%salinity%molecular_weights(i))) then
          this%salinity%molecular_weights(i) = this%realization%patch% &
            reaction%primary_spec_molar_wt(this%salinity%ispecies(i))
        endif
      enddo
      ! execute twice to initiaize both m_nacl entries
      call this%Evaluate(time,ierr)
      call this%Evaluate(time,ierr)
  end select

end subroutine PMAuxiliaryInitializeRun

! ************************************************************************** !

recursive subroutine PMAuxiliaryFinalizeRun(this)
  !
  ! Finalizes the run
  !
  ! Author: Glenn Hammond
  ! Date: 10/22/18
  !

  implicit none

  class(pm_auxiliary_type) :: this

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMAuxiliaryFinalizeRun

! ************************************************************************** !

subroutine PMAuxiliaryEvolvingStrata(this,time,ierr)
  !
  ! Initializes auxiliary process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  use Init_Subsurface_module

  implicit none

  class(pm_auxiliary_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  ierr = 0
  call InitSubsurfAssignMatIDsToRegns(this%realization)
  call InitSubsurfAssignMatProperties(this%realization)
  call InitSubsurfaceSetupZeroArrays(this%realization)

end subroutine PMAuxiliaryEvolvingStrata

! ************************************************************************** !

subroutine PMAuxiliaryInversionMeasurement(this,time,ierr)
  !
  ! Takes measurements for inversion calculation
  !
  ! Author: Glenn Hammond
  ! Date: 05/25/22

  use Inversion_TS_Aux_module
  use Option_module
  use Realization_Subsurface_class
  use Utility_module

  implicit none

  class(pm_auxiliary_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(inversion_forward_aux_type), pointer :: inversion_forward_aux

  ierr = 0
  inversion_forward_aux => this%realization%patch%aux%inversion_forward_aux
  if (associated(inversion_forward_aux)) then
    if (Equal(inversion_forward_aux%sync_times( &
                inversion_forward_aux%isync_time),time)) then
      call RealizationGetObservedVariables(this%realization)
      call InversionForwardAuxMeasure(inversion_forward_aux,time,this%option)
    else
      call PrintMsg(this%option,'  No measurement requested at this time.')
    endif
  endif

end subroutine PMAuxiliaryInversionMeasurement

! ************************************************************************** !

subroutine PMAuxiliaryInversionAdjoint(this,time,ierr)
  !
  ! Initializes auxiliary process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  use Inversion_TS_Aux_module
  use Realization_Subsurface_class

  implicit none

  class(pm_auxiliary_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  ierr = 0
  if (associated(this%realization%patch%aux%inversion_forward_aux)) then
    call InversionForwardAuxStep(this%realization%patch%aux% &
                                   inversion_forward_aux,time)
  endif

end subroutine PMAuxiliaryInversionAdjoint

! ************************************************************************** !

subroutine PMAuxiliarySalinity(this,time,ierr)
  !
  ! Initializes auxiliary process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/16
  !
  use Reactive_Transport_Aux_module
  use EOS_Water_module
  use Global_Aux_module

  implicit none

  class(pm_auxiliary_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  PetscInt :: ghosted_id, i, j, ispecies, num_auxvars
  PetscReal :: sum_mass_species, xnacl, mass_h2o
  PetscReal :: dw_mol, dw_dp, dw_dt
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscInt, parameter :: iphase = 1

  ierr = 0
  do j = 1, 2
    if (j == 1) then
      rt_auxvars => this%realization%patch%aux%RT%auxvars
      global_auxvars => this%realization%patch%aux%Global%auxvars
      num_auxvars = this%realization%patch%aux%RT%num_aux
    else
      rt_auxvars => this%realization%patch%aux%RT%auxvars_bc
      global_auxvars => this%realization%patch%aux%Global%auxvars_bc
      num_auxvars = this%realization%patch%aux%RT%num_aux_bc
    endif
    do ghosted_id = 1, num_auxvars
      sum_mass_species = 0.d0
      do i = 1, this%salinity%nspecies
        ispecies = this%salinity%ispecies(i)
        sum_mass_species = sum_mass_species + &
          rt_auxvars(ghosted_id)%total(ispecies,iphase)* &
          this%salinity%molecular_weights(i) ! mol/L * g/mol = g/L and
                                             !   g/L => kg/m^3
      enddo

      ! Save NaCl from pervious timestep
      global_auxvars(ghosted_id)%m_nacl(TWO_INTEGER) = &
        global_auxvars(ghosted_id)%m_nacl(ONE_INTEGER)

      ! Compute NaCl for new timestep
      call EOSWaterDensity(global_auxvars(ghosted_id)%temp, &
                           global_auxvars(ghosted_id)%pres(1), &
                           mass_h2o,dw_mol,dw_dp,dw_dt,ierr)
      xnacl = sum_mass_species / (sum_mass_species + mass_h2o)
      global_auxvars(ghosted_id)%m_nacl(iphase) = xnacl
    enddo
  enddo

end subroutine PMAuxiliarySalinity

! ************************************************************************** !

subroutine PMAuxiliaryInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 04/21/2016
  !

  implicit none

  class(pm_auxiliary_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMAuxiliaryInputRecord

! ************************************************************************** !

subroutine PMAuxiliaryDestroy(this)
  !
  ! Destroys auxiliary process model
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none

  class(pm_auxiliary_type) :: this

  call PMBaseDestroy(this)

  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%option)
  nullify(this%output_option)

  if (associated(this%salinity)) then
    deallocate(this%salinity)
    nullify(this%salinity)
  endif

end subroutine PMAuxiliaryDestroy

end module PM_Auxiliary_class
