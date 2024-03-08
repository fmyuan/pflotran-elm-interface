module Option_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Communicator_Aux_module
  use Driver_class
  use Option_Checkpoint_module
  use Option_Flow_module
  use Option_Transport_module
  use Option_Geophysics_module
  use Option_Inversion_module
  use Option_Parameter_module
  use Print_module

  implicit none

  private

  type, public :: option_type

    type(print_flags_type), pointer :: print_flags
    type(flow_option_type), pointer :: flow
    type(transport_option_type), pointer :: transport
    type(geophysics_option_type), pointer :: geophysics
    type(checkpoint_option_type), pointer :: checkpoint
    type(inversion_option_type), pointer :: inversion
    type(parameter_option_type), pointer :: parameter
    type(comm_type), pointer :: comm
    class(driver_type), pointer :: driver

    PetscInt :: id                         ! id of realization

    PetscMPIInt :: mycomm                  ! PETSC_COMM_WORLD
    PetscMPIInt :: myrank                  ! rank in PETSC_COMM_WORLD

    PetscMPIInt :: hdf5_read_group_size, hdf5_write_group_size
    PetscBool :: broadcast_read
    PetscBool :: blocking
    PetscBool :: error_while_nonblocking

    character(len=MAXSTRINGLENGTH) :: io_buffer

    PetscInt :: fid_out
    PetscInt :: fid_inputrecord

    ! defines the mode (e.g. mph, richards, vadose, etc.
    character(len=MAXWORDLENGTH) :: flowmode
    PetscInt :: iflowmode
    PetscInt :: iflow_sub_mode
    character(len=MAXWORDLENGTH) :: tranmode
    PetscInt :: itranmode
    character(len=MAXWORDLENGTH) :: geopmode
    PetscInt :: igeopmode

    PetscInt :: nphase
    PetscInt :: liquid_phase
    PetscInt :: gas_phase
    PetscInt :: hydrate_phase
    PetscInt :: ice_phase
    PetscInt :: precipitate_phase
    PetscInt :: pure_water_phase ! for storing pure water properties
    PetscInt :: pure_brine_phase ! for storing pure brine properties
    PetscInt :: trapped_gas_phase
    PetscInt :: nflowdof
    PetscInt :: nflowspec
    PetscInt :: nmechdof
    PetscInt :: nsec_cells
    PetscInt :: num_table_indices

    PetscBool :: geomech_on
    PetscBool :: geomech_initial
    PetscInt :: ngeomechdof
    PetscInt :: n_stress_strain_dof
    PetscReal :: geomech_time
    PetscInt :: geomech_subsurf_coupling
    PetscReal :: geomech_gravity(3)
    PetscBool :: sec_vars_update

    PetscInt :: air_pressure_id
    PetscInt :: co2_pressure_id
    PetscInt :: capillary_pressure_id
    PetscInt :: vapor_pressure_id
    PetscInt :: reduced_vapor_pressure_id
    PetscInt :: saturation_pressure_id
    PetscInt :: water_id  ! index of water component dof
    PetscInt :: air_id  ! index of air component dof
    PetscInt :: co2_id ! index of co2 component dof
    PetscInt :: energy_id  ! index of energy dof
    PetscInt :: salt_id ! index of salt dof

    PetscInt :: ntrandof

    PetscInt :: ngeopdof ! geophysics # of dof

    PetscBool :: coupled_well

    PetscInt :: iflag
    PetscInt :: ierror
    PetscInt :: status
    PetscBool :: input_record
    PetscInt :: verbosity  ! Values >0 indicate additional console output.
    PetscBool :: keyword_logging
    PetscBool :: keyword_logging_screen_output
    character(len=MAXSTRINGLENGTH) :: keyword_log
    character(len=MAXSTRINGLENGTH) :: keyword_buf
    PetscInt :: keyword_block_map(20)
    PetscInt :: keyword_block_count

    ! Program options
    PetscBool :: use_matrix_free  ! If true, do not form the Jacobian.

    PetscBool :: use_isothermal
    PetscBool :: use_sc           ! If true, multiple continuum formulation is used.
    PetscReal :: flow_time, tran_time, time  ! The time elapsed in the simulation.
    PetscReal :: flow_dt ! The size of the time step.
    PetscReal :: tran_dt
    PetscReal :: dt
    PetscBool :: match_waypoint

    PetscReal :: gravity(3)

    PetscReal :: scale

    PetscReal :: m_nacl

    PetscInt :: ideriv
    PetscInt :: idt_switch

    PetscBool :: converged
    PetscInt :: convergence

    PetscReal :: infnorm_res_sec  ! inf. norm of secondary continuum rt residual

!   table lookup
    PetscInt :: itable
    PetscInt :: co2eos
    character(len=MAXSTRINGLENGTH) :: co2_database_filename

    PetscBool :: restart_flag
    PetscReal :: restart_time
    character(len=MAXSTRINGLENGTH) :: restart_filename
    character(len=MAXSTRINGLENGTH) :: input_filename

    PetscBool :: wallclock_stop_flag
    PetscLogDouble :: wallclock_stop_time

    PetscInt :: log_stage(10)

    PetscBool :: numerical_derivatives_multi_coupling
    PetscBool :: compute_statistics
    PetscBool :: compute_mass_balance_new
    PetscBool :: mass_bal_detailed
    PetscBool :: use_touch_options
    PetscInt :: io_handshake_buffer_size

    character(len=MAXSTRINGLENGTH) :: initialize_flow_filename
    character(len=MAXSTRINGLENGTH) :: initialize_transport_filename

    character(len=MAXSTRINGLENGTH) :: global_prefix
    character(len=MAXWORDLENGTH) :: group_prefix

    PetscBool :: use_matrix_buffer
    PetscBool :: force_newton_iteration
    PetscBool :: out_of_table

    ! Specify secondary continuum solver
    PetscInt :: secondary_continuum_solver     ! Specify secondary continuum solver

    ! when the scaling factor is too small, stop in reactive transport
    PetscReal :: min_allowable_scale

    PetscBool :: print_ekg

  end type option_type

  interface OptionCreate
    module procedure OptionCreate1
    module procedure OptionCreate2
  end interface

  interface PrintMsg
    module procedure PrintMsg1
    module procedure PrintMsg2
    module procedure PrintMsg3
  end interface

  interface PrintMsgNoAdvance
    module procedure PrintMsgNoAdvance1
    module procedure PrintMsgNoAdvance2
  end interface

  interface PrintMsgAnyRank
    module procedure PrintMsgAnyRank1
    module procedure PrintMsgAnyRank2
  end interface

  interface PrintMsgByRank
    module procedure PrintMsgByRank1
    module procedure PrintMsgByRank2
  end interface

  interface PrintErrMsgByRank
    module procedure PrintErrMsgByRank1
    module procedure PrintErrMsgByRank2
  end interface

  interface PrintErrMsgNoStopByRank
    module procedure PrintErrMsgNoStopByRank1
    module procedure PrintErrMsgNoStopByRank2
  end interface

  interface PrintErrMsg
    module procedure PrintErrMsg1
    module procedure PrintErrMsg2
  end interface

  interface PrintWrnMsg
    module procedure PrintWrnMsg1
    module procedure PrintWrnMsg2
  end interface

  interface OptionIsIORank
    module procedure OptionIsIORank1
    module procedure OptionIsIORank2
  end interface

  public :: OptionCreate, &
            OptionSetDriver, &
            OptionSetComm, &
            OptionSetInversionOption, &
            OptionCheckCommandLine, &
            PrintErrMsg, &
            PrintErrMsgToDev, &
            PrintErrMsgByRank, &
            PrintErrMsgByRankToDev, &
            PrintWrnMsg, &
            PrintMsg, &
            PrintMsgNoAdvance, &
            PrintMsgAnyRank, &
            PrintMsgByRank, &
            PrintMsgByCell, &
            PrintErrMsgNoStopByRank, &
            PrintVerboseMsg, &
            OptionCheckTouch, &
            OptionPrintToScreen, &
            OptionPrintToFile, &
            OptionInitRealization, &
            OptionMeanVariance, &
            OptionMaxMinMeanVariance, &
            OptionPrintPFLOTRANHeader, &
            OptionSetBlocking, &
            OptionCheckNonBlockingError, &
            OptionIsIORank, &
            OptionCreatePrintHandler, &
            OptionDestroy

contains

! ************************************************************************** !

function OptionCreate1()
  !
  ! Allocates and initializes a new Option object
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !
  implicit none

  type(option_type), pointer :: OptionCreate1

  type(option_type), pointer :: option

  allocate(option)
  option%print_flags => PrintCreateFlags()
  option%flow => OptionFlowCreate()
  option%transport => OptionTransportCreate()
  option%geophysics => OptionGeophysicsCreate()
  option%parameter => OptionParameterCreate()
  nullify(option%checkpoint)
  nullify(option%inversion)
  nullify(option%driver)
  nullify(option%comm)

  ! DO NOT initialize members of the option type here.  One must decide
  ! whether the member needs initialization once for all stochastic
  ! simulations or initialization for every realization (e.g. within multiple
  ! stochastic simulations).  This is done in OptionInitAll() and
  ! OptionInitRealization()
  call OptionInitAll(option)
  OptionCreate1 => option

end function OptionCreate1

! ************************************************************************** !

function OptionCreate2(outer_option)
  !
  ! Same as OptionCreate() but increments file pointers
  !
  ! Author: Glenn Hammond
  ! Date: 11/21/22
  !
  use String_module

  implicit none

  type(option_type), pointer :: outer_option

  type(option_type), pointer :: OptionCreate2

  type(option_type), pointer :: option
  PetscInt :: fid_out
  PetscErrorCode :: ierr

  option => OptionCreate()
  if (associated(outer_option)) then
    call PrintInitFlags(option%print_flags,outer_option%print_flags)
    call MPI_Allreduce(outer_option%fid_out,fid_out,ONE_INTEGER_MPI, &
                       MPI_INTEGER,MPI_MAX,outer_option%mycomm, &
                       ierr);CHKERRQ(ierr)
    if (fid_out <= 0) then
      outer_option%io_buffer = 'outer_option%fid_out not set properly in &
        &OptionCreate2: (' // StringWrite(outer_option%fid_out) // ').'
      call PrintErrMsg(outer_option)
    endif
    if (fid_out + 1 > MAX_OUT_UNIT) then
      option%io_buffer = 'The maximum output file id (MAX_OUT_UNIT) has been &
        &exceeded. Please increase MAX_OUT_UNIT in pflotran_constants.F90.'
      call PrintErrMsg(option)
    endif
    if (outer_option%fid_out > 0) then
      option%fid_out = outer_option%fid_out + 1
    endif
    option%group_prefix = outer_option%group_prefix
  else
    option%io_buffer = 'outer_option not associated in OptionCreate2.'
    call PrintErrMsg(option)
  endif
  OptionCreate2 => option

end function OptionCreate2

! ************************************************************************** !

subroutine OptionSetDriver(option,driver)

  implicit none

  type(option_type) :: option
  class(driver_type), pointer :: driver

  option%driver => driver
  call PrintInitFlags(option%print_flags,driver%print_flags)

end subroutine OptionSetDriver

! ************************************************************************** !

subroutine OptionSetComm(option,comm)

  ! If the MPI communicator is split, we need to update the values local
  ! values in option

  use Communicator_Aux_module

  implicit none

  type(option_type) :: option
  type(comm_type), pointer :: comm

  call CommResetStartTime(comm)
  option%comm => comm
  option%mycomm = comm%communicator
  option%myrank = comm%rank

end subroutine OptionSetComm

! ************************************************************************** !

subroutine OptionSetInversionOption(option,inversion_option)

  implicit none

  type(option_type) :: option
  type(inversion_option_type), pointer :: inversion_option

  option%inversion => inversion_option

end subroutine OptionSetInversionOption

! ************************************************************************** !

subroutine OptionInitAll(option)
  !
  ! Initializes all option variables
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

  implicit none

  type(option_type) :: option

  ! These variables should only be initialized once at the beginning of a
  ! PFLOTRAN run (regardless of whether stochastic)

  call PrintInitFlags(option%print_flags)
  call OptionFlowInitAll(option%flow)
  call OptionTransportInitAll(option%transport)

  option%id = 0

  option%mycomm = 0
  option%myrank = 0

  option%group_prefix = ''
  option%global_prefix = ''

  option%broadcast_read = PETSC_FALSE
  option%hdf5_read_group_size = 0
  option%hdf5_write_group_size = 0
  option%blocking = PETSC_TRUE
  option%error_while_nonblocking = PETSC_FALSE

  option%input_record = PETSC_FALSE
  option%verbosity = 0
  option%keyword_logging = PETSC_TRUE
  option%keyword_logging_screen_output = PETSC_FALSE
  option%keyword_log = ''
  option%keyword_buf = ''
  option%keyword_block_map(:) = 0
  option%keyword_block_count = 0

  option%input_filename = ''

  option%out_of_table = PETSC_FALSE

  call OptionInitRealization(option)

end subroutine OptionInitAll

! ************************************************************************** !

subroutine OptionInitRealization(option)
  !
  ! Initializes option variables specific to a single
  ! realization
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

  implicit none

  type(option_type) :: option

  ! These variables should be initialized once at the beginning of every
  ! PFLOTRAN realization or simulation of a single realization
  call OptionFlowInitRealization(option%flow)
  call OptionTransportInitRealization(option%transport)

  option%fid_out = 0
  option%fid_inputrecord = INPUT_RECORD_UNIT

  option%iflag = 0
  option%ierror = 0
  option%io_buffer = ''

  option%use_isothermal = PETSC_FALSE
  option%use_matrix_free = PETSC_FALSE
  option%use_sc = PETSC_FALSE

  option%flowmode = ""
  option%iflowmode = NULL_MODE
  option%iflow_sub_mode = NULL_MODE
  option%nflowdof = 0
  option%nflowspec = 0
  option%nmechdof = 0
  option%nsec_cells = 0
  option%num_table_indices = 0

  option%geomech_on = PETSC_FALSE
  option%geomech_initial = PETSC_FALSE
  option%ngeomechdof = 0
  option%n_stress_strain_dof = 0
  option%geomech_time = 0.d0
  option%geomech_subsurf_coupling = 0
  option%geomech_gravity(:) = 0.d0
  option%geomech_gravity(3) = -1.d0*EARTH_GRAVITY    ! m/s^2

  option%tranmode = ""
  option%itranmode = NULL_MODE
  option%ntrandof = 0

  option%geopmode = ""
  option%igeopmode = NULL_MODE
  option%ngeopdof = 0

  option%nphase = 0

  option%liquid_phase  = UNINITIALIZED_INTEGER
  option%gas_phase     = UNINITIALIZED_INTEGER
  option%hydrate_phase = UNINITIALIZED_INTEGER
  option%ice_phase = UNINITIALIZED_INTEGER
  option%precipitate_phase = UNINITIALIZED_INTEGER
  option%pure_water_phase = UNINITIALIZED_INTEGER
  option%pure_brine_phase = UNINITIALIZED_INTEGER

  option%air_pressure_id = 0
  option%co2_pressure_id = 0
  option%capillary_pressure_id = 0
  option%vapor_pressure_id = 0
  option%reduced_vapor_pressure_id = 0
  option%saturation_pressure_id = 0

  option%water_id = 0
  option%air_id = 0
  option%co2_id = 0
  option%energy_id = 0
  option%salt_id = 0

!-----------------------------------------------------------------------
      ! Initialize some parameters to sensible values.  These are parameters
      ! which should be set via the command line or the input file, but it
      ! seems good practice to set them to sensible values when a pflowGrid
      ! is created.
!-----------------------------------------------------------------------
  option%converged = PETSC_FALSE
  option%convergence = CONVERGENCE_OFF

  option%infnorm_res_sec = 0.d0

  !set scale factor for heat equation, i.e. use units of MJ for energy
  option%scale = 1.d-6

  option%ideriv = 1

  option%gravity(:) = 0.d0
  option%gravity(3) = -1.d0*EARTH_GRAVITY ! m/s^2

  !physical constants and defult variables
!  option%difaq = 1.d-9 ! m^2/s read from input file
!  option%difaq = 0.d0
!  option%delhaq = 12.6d0 ! kJ/mol read from input file
!  option%eqkair = 1.d10 ! Henry's constant for air: Xl = eqkair * pa

  ! default brine concentrations
  option%m_nacl = 0.d0

!  option%disp = 0.d0

  option%restart_flag = PETSC_FALSE
  option%restart_filename = ""
  option%restart_time = UNINITIALIZED_DOUBLE

  option%wallclock_stop_flag = PETSC_FALSE
  option%wallclock_stop_time = 0.d0

  option%log_stage = 0

  option%numerical_derivatives_multi_coupling = PETSC_FALSE
  option%compute_statistics = PETSC_FALSE
  option%compute_mass_balance_new = PETSC_FALSE
  option%mass_bal_detailed = PETSC_FALSE

  option%use_touch_options = PETSC_FALSE

  option%time = 0.d0
  option%flow_dt = 0.d0
  option%tran_dt = 0.d0
  option%dt = 0.d0
  option%match_waypoint = PETSC_FALSE

  option%io_handshake_buffer_size = 0

  option%initialize_flow_filename = ''
  option%initialize_transport_filename = ''

  option%itable = 0
  option%co2eos = EOS_SPAN_WAGNER
  option%co2_database_filename = ''

! option%idt_switch = 1
  option%idt_switch = -1

  option%use_matrix_buffer = PETSC_FALSE
  option%status = PROCEED
  option%force_newton_iteration = PETSC_FALSE
  option%secondary_continuum_solver = 1

  ! when the scaling factor is too small, stop in reactive transport
  option%min_allowable_scale = 1.0d-10

  option%print_ekg = PETSC_FALSE

  option%coupled_well = PETSC_FALSE

end subroutine OptionInitRealization

! ************************************************************************** !

subroutine OptionCheckCommandLine(option)
  !
  ! Checks all PETSc options on input
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  implicit none

  type(option_type) :: option

  PetscBool :: option_found
  PetscErrorCode :: ierr

  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           "-buffer_matrix",option%use_matrix_buffer, &
                           ierr);CHKERRQ(ierr)
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-snes_mf", &
                           option%use_matrix_free,ierr);CHKERRQ(ierr)
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           "-use_isothermal",option%use_isothermal, &
                           ierr);CHKERRQ(ierr)
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-use_sc", &
                           option%use_sc,ierr);CHKERRQ(ierr)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                             '-restart',option%restart_filename, &
                             option%restart_flag,ierr);CHKERRQ(ierr)
  ! check on possible modes
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           "-use_richards",option_found,ierr);CHKERRQ(ierr)
  if (option_found) option%flowmode = "richards"
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-use_thc", &
                           option_found,ierr);CHKERRQ(ierr)
  if (option_found) option%flowmode = "thc"
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-use_mph", &
                           option_found,ierr);CHKERRQ(ierr)
  if (option_found) option%flowmode = "mph"
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           "-use_flash2",option_found,ierr);CHKERRQ(ierr)
  if (option_found) option%flowmode = "flash2"

end subroutine OptionCheckCommandLine

! ************************************************************************** !

subroutine PrintErrMsg1(option)
  !
  ! Prints the error message from p0 and stops
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  implicit none

  type(option_type) :: option

  call PrintErrMsg2(option,option%io_buffer)

end subroutine PrintErrMsg1

! ************************************************************************** !

subroutine PrintErrMsg2(option,string)
  !
  ! Prints the error message from p0 and stops
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !
  use Print_module

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  print_handler%byrank = PETSC_FALSE
  call PrintErrorMessage(print_handler,string)
  call PrintDestroyHandler(print_handler)
  if (.not.option%blocking) then
    option%error_while_nonblocking = PETSC_TRUE
  endif

end subroutine PrintErrMsg2

! ************************************************************************** !

subroutine OptionCheckNonBlockingError(option)
  !
  ! Checks whether error_while_nonblocking was set and stops if TRUE
  !
  ! Author: Glenn Hammond
  ! Date: 06/2/19
  !

  implicit none

  type(option_type) :: option

  PetscBool :: petsc_initialized
  PetscMPIInt :: mpi_int
  PetscErrorCode :: ierr

  mpi_int = 1
  call MPI_Allreduce(MPI_IN_PLACE,option%error_while_nonblocking,mpi_int, &
                     MPI_LOGICAL,MPI_LOR,option%mycomm,ierr);CHKERRQ(ierr)
  if (option%error_while_nonblocking) then
    call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
    call PetscInitialized(petsc_initialized,ierr);CHKERRQ(ierr)
    if (petsc_initialized) then
      call PetscFinalize(ierr);CHKERRQ(ierr)
    endif
    call exit(EXIT_USER_ERROR)
  endif

end subroutine OptionCheckNonBlockingError

! ************************************************************************** !

subroutine PrintErrMsgByRank1(option)
  !
  ! Prints the error message from processor with error along
  ! with rank
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/11
  !

  implicit none

  type(option_type) :: option

  call PrintErrMsgByRank2(option,option%io_buffer)

end subroutine PrintErrMsgByRank1

! ************************************************************************** !

subroutine PrintErrMsgByRank2(option,string)
  !
  ! Prints the error message from processor with error along
  ! with rank
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/11
  !
  use Print_module

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  print_handler%byrank = PETSC_TRUE
  call PrintErrorMessage(print_handler,string)
  call PrintDestroyHandler(print_handler)

end subroutine PrintErrMsgByRank2

! ************************************************************************** !

subroutine PrintErrMsgNoStopByRank1(option)
  !
  ! Prints the error message from processor with error along
  ! with rank
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/11
  !

  implicit none

  type(option_type) :: option

  call PrintErrMsgNoStopByRank2(option,option%io_buffer)

end subroutine PrintErrMsgNoStopByRank1

! ************************************************************************** !

subroutine PrintErrMsgNoStopByRank2(option,string)
  !
  ! Prints the error message from processor with error along
  ! with rank
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/11
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  type(print_handler_type), pointer :: print_handler

  print_handler => OptionCreatePrintHandler(option)
  print_handler%blocking = PETSC_FALSE ! keeps it from stopping
  print_handler%byrank = PETSC_TRUE
  call PrintErrorMessage(print_handler,string)
  call PrintDestroyHandler(print_handler)

end subroutine PrintErrMsgNoStopByRank2

! ************************************************************************** !

subroutine PrintErrMsgToDev(option,string)
  !
  ! Prints the error message from p0, appends a request to submit input
  ! deck to pflotran-dev, and stops.  The reverse order of arguments is
  ! to avoid conflict with variants of PrintErrMsg()
  !
  ! Author: Glenn Hammond
  ! Date: 07/26/18
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  if (len_trim(string) > 0) then
    option%io_buffer = trim(option%io_buffer) // &
      ' Please email pflotran-dev@googlegroups.com and ' // &
      trim(adjustl(string)) // '.'
  else
    option%io_buffer = trim(option%io_buffer) // &
      ' Please email pflotran-dev@googlegroups.com.'
  endif
  call PrintErrMsg(option)

end subroutine PrintErrMsgToDev

! ************************************************************************** !

subroutine PrintErrMsgByRankToDev(option,string)
  !
  ! Prints the error message from processor with error along
  ! with rank. The reverse order of arguments is to avoid conflict with
  ! variants of PrintErrMsg()
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/11
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  if (len_trim(string) > 0) then
    option%io_buffer = trim(option%io_buffer) // &
      ' Please email pflotran-dev@googlegroups.com and ' // &
      trim(adjustl(string)) // '.'
  else
    option%io_buffer = trim(option%io_buffer) // &
      ' Please email pflotran-dev@googlegroups.com.'
  endif
  call PrintErrMsgByRank(option)

end subroutine PrintErrMsgByRankToDev

! ************************************************************************** !

subroutine PrintWrnMsg1(option)
  !
  ! Prints the warning message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  implicit none

  type(option_type) :: option

  call PrintWrnMsg2(option,option%io_buffer)

end subroutine PrintWrnMsg1

! ************************************************************************** !

subroutine PrintWrnMsg2(option,string)
  !
  ! Prints the warning message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  character(len=:), allocatable :: local_string
  PetscBool, parameter :: advance_ = PETSC_TRUE
  PetscBool, parameter :: byrank = PETSC_FALSE

  allocate(local_string,source = ' WARNING: ' // trim(string))
  call PrintMessage(option%print_flags,option%comm,option%fid_out, &
                    local_string,advance_,byrank)

end subroutine PrintWrnMsg2

! ************************************************************************** !

subroutine PrintMsgNoAdvance1(option)
  !
  ! Prints output to the screen and/or output file
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22
  !
  implicit none

  type(option_type) :: option


  call PrintMsgNoAdvance2(option,option%io_buffer)

end subroutine PrintMsgNoAdvance1

! ************************************************************************** !

subroutine PrintMsgNoAdvance2(option,string)
  !
  ! Prints output to the screen and/or output file
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22
  !
  implicit none

  type(option_type) :: option
  character(len=*) :: string

  PetscBool, parameter :: advance_ = PETSC_FALSE
  PetscBool, parameter :: byrank = PETSC_FALSE

  call PrintMessage(option%print_flags,option%comm,option%fid_out, &
                    string,advance_,byrank)

end subroutine PrintMsgNoAdvance2

! ************************************************************************** !

subroutine PrintMsg1(option)
  !
  ! Prints output to the screen and/or output file
  !
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  !
  implicit none

  type(option_type) :: option

  call PrintMsg2(option,option%io_buffer)

end subroutine PrintMsg1

! ************************************************************************** !

subroutine PrintMsg2(option,string)
  !
  ! Prints output to the screen and/or output file
  !
  ! Author: Glenn Hammond
  ! Date: 11/14/07, 06/03/22
  !
  implicit none

  type(option_type) :: option
  character(len=*) :: string

  PetscBool, parameter :: advance_ = PETSC_TRUE

  call PrintMsg3(option,string,advance_)

end subroutine PrintMsg2

! ************************************************************************** !

subroutine PrintMsg3(option,string,advance_)
  !
  ! Prints output to the screen and/or output file
  !
  ! Author: Glenn Hammond
  ! Date: 11/14/07, 06/03/22
  !
  implicit none

  type(option_type) :: option
  character(len=*) :: string
  PetscBool :: advance_

  PetscBool, parameter :: byrank = PETSC_FALSE

  call PrintMessage(option%print_flags,option%comm,option%fid_out, &
                    string,advance_,byrank)

end subroutine PrintMsg3

! ************************************************************************** !

subroutine PrintMsgAnyRank1(option)
  !
  ! Prints the message from any processor core
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  !
  implicit none

  type(option_type) :: option

  if (option%print_flags%print_to_screen) then
    call PrintMsgAnyRank2(option%io_buffer)
  endif

end subroutine PrintMsgAnyRank1

! ************************************************************************** !

subroutine PrintMsgAnyRank2(string)
  !
  ! Prints the message from any processor core
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  !
  implicit none

  character(len=*) :: string

  print *, trim(string)

end subroutine PrintMsgAnyRank2

! ************************************************************************** !

subroutine PrintMsgByRank1(option)
  !
  ! Prints a message from processor along with rank
  !
  ! Author: Glenn Hammond
  ! Date: 03/27/12
  !

  implicit none

  type(option_type) :: option

  call PrintMsgByRank2(option,option%io_buffer)

end subroutine PrintMsgByRank1

! ************************************************************************** !

subroutine PrintMsgByRank2(option,string)
  !
  ! Prints a message from processor along with rank
  !
  ! Author: Glenn Hammond
  ! Date: 03/27/12
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  PetscBool, parameter :: advance_ = PETSC_TRUE
  PetscBool, parameter :: byrank = PETSC_TRUE

  call PrintMessage(option%print_flags,option%comm,option%fid_out, &
                    string,advance_,byrank)

end subroutine PrintMsgByRank2

! ************************************************************************** !

subroutine PrintMsgByCell(option,cell_id,string)
  !
  ! Prints the message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  !
  use String_module

  implicit none

  type(option_type) :: option
  PetscInt :: cell_id
  character(len=*) :: string

  option%io_buffer = trim(string) // ' for cell ' // &
                     StringWrite(cell_id) // '.'
  call PrintMsgByRank(option)

end subroutine PrintMsgByCell

! ************************************************************************** !

subroutine PrintVerboseMsg(option)
  !
  ! Prints the message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  !

  implicit none

  type(option_type) :: option

  if (option%verbosity > 0) then
    call PrintMsg(option,option%io_buffer)
  endif

end subroutine PrintVerboseMsg

! ************************************************************************** !

function OptionCheckTouch(option,filename)
  !
  ! Users can steer the code by touching files.
  !
  ! Author: Glenn Hammond
  ! Date: 03/04/08
  !

  implicit none

  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: filename

  PetscBool :: OptionCheckTouch

  PetscInt :: ios
  PetscInt :: fid = 86
  PetscBool :: is_io_rank
  PetscErrorCode :: ierr

  OptionCheckTouch = PETSC_FALSE

  is_io_rank = CommIsIORank(option%comm)

  if (is_io_rank) open(unit=fid,file=trim(filename),status='old',iostat=ios)
  call MPI_Bcast(ios,ONE_INTEGER_MPI,MPIU_INTEGER,option%comm%io_rank, &
                 option%mycomm,ierr);CHKERRQ(ierr)

  if (ios == 0) then
    if (is_io_rank) close(fid,status='delete')
    OptionCheckTouch = PETSC_TRUE
  endif

end function OptionCheckTouch

! ************************************************************************** !

function OptionPrintToScreen(option)
  !
  ! Determines whether printing should occur
  !
  ! Author: Glenn Hammond
  ! Date: 12/09/08
  !

  implicit none

  type(option_type) :: option

  PetscBool :: OptionPrintToScreen

  OptionPrintToScreen = PrintToScreen(option%print_flags,option%comm)

end function OptionPrintToScreen

! ************************************************************************** !

function OptionPrintToFile(option)
  !
  ! Determines whether printing to file should occur
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/09
  !

  implicit none

  type(option_type) :: option

  PetscBool :: OptionPrintToFile

  OptionPrintToFile = PrintToFile(option%print_flags,option%comm)

end function OptionPrintToFile

! ************************************************************************** !

subroutine OptionMaxMinMeanVariance(value,max,min,mean,variance, &
                                    calculate_variance,option)
  !
  ! Calculates the maximum, minumum, mean and
  ! optionally variance of a number across processor
  ! cores
  !
  ! Author: Glenn Hammond
  ! Date: 06/01/10
  !

  implicit none

  type(option_type) :: option
  PetscReal :: value
  PetscReal :: max
  PetscReal :: min
  PetscReal :: mean
  PetscReal :: variance
  PetscBool :: calculate_variance

  PetscReal :: temp_real_in(2), temp_real_out(2)
  PetscErrorCode :: ierr

  temp_real_in(1) = value
  temp_real_in(2) = -1.d0*value
  call MPI_Allreduce(temp_real_in,temp_real_out,TWO_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                     ierr);CHKERRQ(ierr)
  max = temp_real_out(1)
  min = -1.d0*temp_real_out(2)

  call OptionMeanVariance(value,mean,variance,calculate_variance,option)

end subroutine OptionMaxMinMeanVariance

! ************************************************************************** !

subroutine OptionMeanVariance(value,mean,variance,calculate_variance,option)
  !
  ! Calculates the mean and optionally variance of a number
  ! across processor cores
  !
  ! Author: Glenn Hammond
  ! Date: 05/29/10
  !

  implicit none

  type(option_type) :: option
  PetscReal :: value
  PetscReal :: mean
  PetscReal :: variance
  PetscBool :: calculate_variance

  PetscReal :: temp_real
  PetscErrorCode :: ierr

  call MPI_Allreduce(value,temp_real,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  mean = temp_real / dble(option%comm%size)

  if (calculate_variance) then
    temp_real = value-mean
    temp_real = temp_real*temp_real
    call MPI_Allreduce(temp_real,variance,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                       ierr);CHKERRQ(ierr)
    variance = variance / dble(option%comm%size)
  endif

end subroutine OptionMeanVariance

! ************************************************************************** !

function OptionIsIORank1(option)
  !
  ! Returns PETSC_TRUE if I/O rank
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/21
  !
  implicit none

  type(option_type) :: option

  PetscBool :: OptionIsIORank1

  OptionIsIORank1 = CommIsIORank(option%comm)

end function OptionIsIORank1

! ************************************************************************** !

function OptionIsIORank2(option,irank)
  !
  ! Returns PETSC_TRUE if I/O rank
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/21
  !

  implicit none

  type(option_type) :: option
  PetscInt :: irank

  PetscBool :: OptionIsIORank2

  OptionIsIORank2 = CommIsIORank(option%comm,irank)

end function OptionIsIORank2

! ************************************************************************** !

subroutine OptionPrintPFLOTRANHeader(option)
  !
  ! Start outer timing.
  !
  ! Author: Glenn Hammond
  ! Date: 04/20/20
  !

  implicit none

  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: version
  character(len=MAXSTRINGLENGTH) :: string

  version = GetVersion()
  write(string,*) len_trim(version)+4
  string = trim(adjustl(string)) // '("=")'
  string = '(/,' // trim(string) // ',/,"  '// &
           trim(version) // &
           '",/,' // trim(string) // ',/)'
  if (OptionPrintToScreen(option)) then
    write(*,string)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,string)
  endif

end subroutine OptionPrintPFLOTRANHeader

! ************************************************************************** !
subroutine OptionSetBlocking(option,flag)
  !
  ! Sets blocking flag
  !
  ! Author: Glenn Hammond
  ! Date: 06/24/19
  !

  implicit none

  type(option_type) :: option
  PetscBool :: flag

  option%blocking = flag

end subroutine OptionSetBlocking

! ************************************************************************** !

function OptionCreatePrintHandler(option)
  !
  ! Creates a print handler with flags for file and screen io
  !
  ! Author: Glenn Hammond
  ! Date: 04/13/23
  !
  use Print_module

  implicit none

  type(option_type) :: option

  type(print_handler_type), pointer :: OptionCreatePrintHandler

  OptionCreatePrintHandler => &
    PrintCreateHandler(option%print_flags, &
                       option%comm, &
                       option%fid_out, &
                       option%driver%exit_code, &
                       PETSC_TRUE, & ! advance
                       option%blocking, &
                       PETSC_FALSE)  ! byrank

end function OptionCreatePrintHandler

! ************************************************************************** !

subroutine OptionDestroy(option)
  !
  ! Deallocates an option
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !
  use Communicator_Aux_module
  use Driver_class

  implicit none

  type(option_type), pointer :: option

  call PrintDestroyFlags(option%print_flags)
  call OptionFlowDestroy(option%flow)
  call OptionTransportDestroy(option%transport)
  call OptionGeophysicsDestroy(option%geophysics)
  call OptionCheckpointDestroy(option%checkpoint)
  call OptionParameterDestroy(option%parameter)
  nullify(option%comm)
  nullify(option%inversion)
  ! never destroy the driver as it was created elsewhere
  nullify(option%driver)

  ! all the below should be placed somewhere other than option.F90

  deallocate(option)
  nullify(option)

end subroutine OptionDestroy

end module Option_module
