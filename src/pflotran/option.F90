module Option_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Communicator_Aux_module
  use Driver_module
  use Option_Flow_module
  use Option_Transport_module
  use Option_Geophysics_module
  use Option_Inversion_module

  implicit none

  private

  type, public :: option_type

    type(flow_option_type), pointer :: flow
    type(transport_option_type), pointer :: transport
    type(geophysics_option_type), pointer :: geophysics

    type(inversion_option_type), pointer :: inversion
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
    PetscInt :: phase_map(MAX_PHASE)
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
    PetscInt :: capillary_pressure_id
    PetscInt :: vapor_pressure_id
    PetscInt :: saturation_pressure_id
    PetscInt :: water_id  ! index of water component dof
    PetscInt :: air_id  ! index of air component dof
    PetscInt :: energy_id  ! index of energy dof

    PetscInt :: ntrandof

    PetscInt :: ngeopdof ! geophysics # of dof

    PetscInt :: iflag
    PetscInt :: ierror
    PetscInt :: status
    PetscBool :: input_record
    ! these flags are for printing outside of time step loop
    ! these flags are for printing within time step loop where printing may
    ! need to be temporarily turned off to accommodate periodic screen outout.
    PetscBool :: print_screen_flag
    PetscBool :: print_file_flag
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
    PetscBool :: use_mc           ! If true, multiple continuum formulation is used.
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
    PetscBool :: use_upwinding
    PetscBool :: out_of_table

    ! Specify secondary continuum solver
    PetscInt :: secondary_continuum_solver     ! Specify secondary continuum solver

    ! when the scaling factor is too small, stop in reactive transport
    PetscReal :: min_allowable_scale

    PetscBool :: print_ekg

  end type option_type

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
            OptionSetInversionOption, &
            OptionUpdateComm, &
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
            OptionDestroy

contains

! ************************************************************************** !

function OptionCreate()
  !
  ! Allocates and initializes a new Option object
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !
  implicit none

  type(option_type), pointer :: OptionCreate

  type(option_type), pointer :: option

  allocate(option)
  option%flow => OptionFlowCreate()
  option%transport => OptionTransportCreate()
  option%geophysics => OptionGeophysicsCreate()
  nullify(option%inversion)
  nullify(option%driver)
  nullify(option%comm)

  ! DO NOT initialize members of the option type here.  One must decide
  ! whether the member needs initialization once for all stochastic
  ! simulations or initialization for every realization (e.g. within multiple
  ! stochastic simulations).  This is done in OptionInitAll() and
  ! OptionInitRealization()
  call OptionInitAll(option)
  OptionCreate => option

end function OptionCreate

! ************************************************************************** !

subroutine OptionSetDriver(option,driver)

  implicit none

  type(option_type) :: option
  class(driver_type), pointer :: driver

  option%driver => driver
  option%comm => driver%comm
  call OptionUpdateComm(option)
  if (option%comm%start_time < 1.d-40) then
    option%io_buffer = 'option%comm%start_time not set. WALLCLOCK_STOP &
      &will not function properly.'
    call PrintErrMsg(option)
  endif
  option%print_screen_flag = driver%PrintToScreen()
  option%print_file_flag = PETSC_FALSE

end subroutine OptionSetDriver

! ************************************************************************** !

subroutine OptionUpdateComm(option)

  ! If the MPI communicator is split, we need to update the values local
  ! values in option

  use Communicator_Aux_module

  implicit none

  type(option_type) :: option

  option%mycomm          = option%comm%mycomm
  option%myrank          = option%comm%myrank

end subroutine OptionUpdateComm

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
  option%print_screen_flag = PETSC_FALSE
  option%print_file_flag = PETSC_FALSE
  option%verbosity = 0
  option%keyword_logging = PETSC_TRUE
  option%keyword_logging_screen_output = PETSC_FALSE
  option%keyword_log = ''
  option%keyword_buf = ''
  option%keyword_block_map(:) = 0
  option%keyword_block_count = 0

  option%input_filename = ''

  option%use_upwinding = PETSC_TRUE

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


  option%fid_out = FORWARD_OUT_UNIT
  option%fid_inputrecord = INPUT_RECORD_UNIT

  option%iflag = 0
  option%ierror = 0
  option%io_buffer = ''

  option%use_isothermal = PETSC_FALSE
  option%use_matrix_free = PETSC_FALSE
  option%use_mc = PETSC_FALSE

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

  option%phase_map = UNINITIALIZED_INTEGER

  option%nphase = 0

  option%liquid_phase  = UNINITIALIZED_INTEGER
  option%gas_phase     = UNINITIALIZED_INTEGER
  option%hydrate_phase = UNINITIALIZED_INTEGER
  option%ice_phase = UNINITIALIZED_INTEGER

  option%air_pressure_id = 0
  option%capillary_pressure_id = 0
  option%vapor_pressure_id = 0
  option%saturation_pressure_id = 0

  option%water_id = 0
  option%air_id = 0
  option%energy_id = 0

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
  PetscInt :: temp_int
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string

  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           "-buffer_matrix",option%use_matrix_buffer, &
                           ierr);CHKERRQ(ierr)
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-snes_mf", &
                           option%use_matrix_free,ierr);CHKERRQ(ierr)
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           "-use_isothermal",option%use_isothermal, &
                           ierr);CHKERRQ(ierr)
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-use_mc", &
                           option%use_mc,ierr);CHKERRQ(ierr)

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

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  PetscBool :: petsc_initialized
  PetscErrorCode :: ierr

  if (OptionPrintToScreen(option)) then
    print *
    print *, 'ERROR: ' // trim(string)
    print *
    print *, 'Stopping!'
  endif
  if (option%blocking) then
    call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
    call PetscInitialized(petsc_initialized,ierr);CHKERRQ(ierr)
    if (petsc_initialized) then
      call PetscFinalize(ierr);CHKERRQ(ierr)
    endif
    select case(option%driver%exit_code)
      case(EXIT_FAILURE)
        call exit(option%driver%exit_code)
      case default
        call exit(EXIT_USER_ERROR)
    end select
  else
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

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  character(len=MAXWORDLENGTH) :: word

  if (option%driver%PrintToScreen()) then
    write(word,*) option%myrank
    print *
    print *, 'ERROR(' // trim(adjustl(word)) // '): ' // trim(string)
    print *
    print *, 'Stopping!'
  endif
  select case(option%driver%exit_code)
    case(EXIT_FAILURE)
      call exit(option%driver%exit_code)
    case default
      call exit(EXIT_USER_ERROR)
  end select

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

  character(len=MAXWORDLENGTH) :: word

  if (option%driver%PrintToScreen()) then
    write(word,*) option%myrank
    print *
    print *, 'ERROR(' // trim(adjustl(word)) // '): ' // trim(string)
    print *
  endif

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

  if (OptionPrintToScreen(option)) print *, 'WARNING: ' // trim(string)

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

  call PrintMsg3(option,string,PETSC_FALSE)

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

  call PrintMsg3(option,option%io_buffer,PETSC_TRUE)

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

  call PrintMsg3(option,string,PETSC_TRUE)

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

  ! note that these flags can be toggled off specific time steps
  if (option%print_screen_flag) then
    if (.not.advance_) then
      write(STDOUT_UNIT,'(a)',advance='no') trim(string)
    else
      write(STDOUT_UNIT,'(a)') trim(string)
    endif
  endif
  if (option%print_file_flag) then
    if (.not.advance_) then
      write(option%fid_out,'(a)',advance='no') trim(string)
    else
      write(option%fid_out,'(a)') trim(string)
    endif
  endif

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

  if (option%driver%PrintToScreen()) call PrintMsgAnyRank2(option%io_buffer)

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

  character(len=MAXWORDLENGTH) :: word

  if (option%driver%PrintToScreen()) then
    write(word,*) option%myrank
    print *, '(' // trim(adjustl(word)) // '): ' // trim(string)
  endif

end subroutine PrintMsgByRank2

! ************************************************************************** !

subroutine PrintMsgByCell(option,cell_id,string)
  !
  ! Prints the message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  !

  implicit none

  type(option_type) :: option
  PetscInt :: cell_id
  character(len=*) :: string

  character(len=MAXWORDLENGTH) :: word

  write(word,*) cell_id
  word = adjustl(word)
  option%io_buffer = trim(string) // ' for cell ' // trim(word) // '.'
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

  is_io_rank = option%driver%IsIORank()

  if (is_io_rank) open(unit=fid,file=trim(filename),status='old',iostat=ios)
  call MPI_Bcast(ios,ONE_INTEGER_MPI,MPIU_INTEGER,option%driver%io_rank, &
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

  OptionPrintToScreen = option%driver%PrintToScreen()

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

  OptionPrintToFile = option%driver%PrintToFile()

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
  mean = temp_real / dble(option%comm%mycommsize)

  if (calculate_variance) then
    temp_real = value-mean
    temp_real = temp_real*temp_real
    call MPI_Allreduce(temp_real,variance,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                       ierr);CHKERRQ(ierr)
    variance = variance / dble(option%comm%mycommsize)
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

  OptionIsIORank1 = option%driver%IsIORank()

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

  OptionIsIORank2 = (irank == option%driver%io_rank)

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
  if (option%driver%IsIORank()) then
    write(string,*) len_trim(version)+4
    string = trim(adjustl(string)) // '("=")'
    string = '(/,' // trim(string) // ',/,"  '// &
             trim(version) // &
             '",/,' // trim(string) // ',/)'
    if (option%driver%PrintToScreen()) then
      write(*,string)
    endif
    if (option%driver%PrintToFile()) then
      write(option%fid_out,string)
    endif
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

subroutine OptionDestroy(option)
  !
  ! Deallocates an option
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !
  use Communicator_Aux_module
  use Driver_module

  implicit none

  type(option_type), pointer :: option

  call OptionFlowDestroy(option%flow)
  call OptionTransportDestroy(option%transport)
  call OptionGeophysicsDestroy(option%geophysics)
  ! never destroy the driver as it was created elsewhere
  nullify(option%inversion)
  nullify(option%driver)
  nullify(option%comm)

  ! all the below should be placed somewhere other than option.F90

  deallocate(option)
  nullify(option)

end subroutine OptionDestroy

end module Option_module
