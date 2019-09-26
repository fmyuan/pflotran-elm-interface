module Well_Data_class

! This is container for information associated with the WELL_DATA type wells

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Well_Type_class

  implicit none

  private

  ! Event types

  PetscInt, parameter, public :: EVENT_OPEN = 1
  PetscInt, parameter, public :: EVENT_SHUT = 2
  PetscInt, parameter, public :: EVENT_STOP = 3
  PetscInt, parameter, public :: EVENT_PINJ = 4
  PetscInt, parameter, public :: EVENT_TINJ = 5
  PetscInt, parameter, public :: EVENT_TARG = 6
  PetscInt, parameter, public :: EVENT_TYPE = 7

  PetscReal, parameter :: w_event_time_eps = 1.0 ! One second

  ! Module variable (number of phases)

  PetscInt :: n_phase = 0

  ! Well target rate array types: target, actual and total flags

  PetscInt, parameter, public :: VALTYPE_TARGET    = 1
  PetscInt, parameter, public :: VALTYPE_ACTUAL    = 2
  PetscInt, parameter, public :: VALTYPE_ACTUALP   = 3
  PetscInt, parameter, public :: VALTYPE_ACTUALI   = 4
  PetscInt, parameter, public :: VALTYPE_TOTALP    = 5
  PetscInt, parameter, public :: VALTYPE_TOTALI    = 6

  ! Global field pressure

  PetscReal :: f_hpav = 0.0

  ! Default wellbore radius

  PetscReal, parameter :: w_vrdef  =  0.25

  ! Flag on use of well_data wells

  PetscBool :: use_wells = PETSC_FALSE

  ! Flag indicating global information set

  PetscBool :: global_info_set = PETSC_FALSE

  ! Warning counts

  PetscInt :: w_unconv_t = 0
  PetscInt :: w_unconv_w = 0
  PetscInt :: w_unconv_b = 0
  PetscInt :: w_mpierr   = 0

  ! This module contains a list of instances of the well_data_type type
  ! The definition of well_data_type follows

  type, public :: well_data_type

      private

      character(len = MAXSTRINGLENGTH), public :: w_name ! well name
      PetscInt  :: w_index                            ! well index
      PetscInt  :: w_itype                            ! well integer type

      PetscBool :: w_radius_set       ! well radius      set
      PetscBool :: w_skin_factor_set  ! well skin factor set
      PetscBool :: w_theta_frac_set   ! well mult factor set
      ! true if the well reference elevation has been read from input
      PetscBool :: w_z_ref_set
      ! true if well drilling direction has been read from input
      PetscBool :: w_const_drill_dir_set

      PetscReal :: w_radius            ! [m] well radius
      PetscReal :: w_skin_factor       ! [-] well skin factor
      PetscReal :: w_theta_frac        ! [-] well mult factor
      PetscReal :: w_z_ref             ! [m] well reference elevation
      PetscInt  :: w_const_drill_dir   ! cmpl drilling direction

      PetscReal :: w_pw                ! Well (bhp) pressure
      PetscReal :: w_pb                ! Well bubble point pressure
      PetscReal , pointer :: w_sp(:)   ! Well saturations
      PetscBool :: w_issat             ! Well state (sat or un-sat)
      PetscReal :: w_trel              ! Well temperature
      PetscBool :: well_solution_set   ! Indicates well solution set

      PetscReal :: w_injection_t       ! Injection temperature
      PetscReal :: w_injection_p       ! Injection pressure

      PetscInt , pointer :: c_ci(:) => null() ! cmpl i-location
      PetscInt , pointer :: c_cj(:) => null() ! cmpl j-location
      PetscInt , pointer :: c_ck(:) => null() ! cmpl k-location

 ! local addresses   for local completions
      PetscInt , pointer :: c_local_id  (:) => null()
 ! ghosted addresses for local completions
      PetscInt , pointer :: c_ghosted_id(:) => null()
 ! global cmpl for each local completion
      PetscInt , pointer :: c_to_cg     (:) => null()
 ! ghosted addresses for all completions
      PetscInt , pointer :: cg_ghosted_id(:) => null()

      PetscBool, pointer :: c_onproc(:) => null() ! on proc flags

      PetscReal, pointer :: c_dx(:) => null()  ! [m] cmpl dx value
      PetscReal, pointer :: c_dy(:) => null()  ! [m] cmpl dy value
      PetscReal, pointer :: c_dz(:) => null()  ! [m] cmpl dz value

      PetscReal, pointer :: c_radius(:) => null()      ! [m] cmpl radius
      PetscReal, pointer :: c_skin_factor(:) => null() ! [-] cmpl skin factor
      PetscReal, pointer :: c_theta_frac(:)  => null() ! [-] cmpl mult factor
      PetscInt , pointer :: c_drill_dir(:)   => null() ! cmpl drilling dirn
      PetscReal, pointer :: c_z(:)           => null() ! cmpl elevation

      PetscReal, pointer  :: c_ccf(:) => null() ! CCF

      PetscInt :: w_status   ! well status (can be open or closed)
      PetscInt :: w_TT       ! well target type

      PetscReal, pointer :: w_targets(:) => null() ! well targets
 ! well actual rates, signed
      PetscReal, pointer :: w_actuals(:) => null()
 ! global-sum well act rates, signed
      PetscReal, pointer :: w_actualsG(:)=> null()
 ! global-sum well prd tots, unsigned
      PetscReal, pointer :: w_totalsPG(:)=> null()
 ! global-sum well inj tots, unsigned
      PetscReal, pointer :: w_totalsIG(:)=> null()

      PetscInt  :: w_nrankw  ! number of ranks on which this well appears
      PetscInt  :: w_ncmpl   ! number of completions (this proc)
      PetscInt  :: w_mcmpl   ! Size of completion list (at least w_ncmpl)
      PetscInt  :: w_ncmplg  ! number of completions (all procs)
      PetscInt  :: w_mcmplg  ! number of all completion list

      PetscBool :: w_ismp   ! Multi-processor well flag

      MPI_Group :: w_group  ! MPI group for this well
      MPI_Comm  :: w_comm   ! MPI communicator for this well

 ! Completions marked for deletion
      PetscBool , pointer :: c_mfd(:) => null()

      PetscInt             :: w_mevent        ! Size of event list
      PetscInt             :: w_nevent        ! Number of events
      PetscInt  , pointer  :: w_event_code(:) =>null()! Event code
      PetscReal , pointer  :: w_event_time(:) =>null()! Event time
      PetscInt  , pointer  :: w_event_ival(:) =>null()! Event integer value
      PetscReal , pointer  :: w_event_rval(:) =>null()! EVent real value
      PetscBool , pointer  :: w_event_used(:) =>null()! Event code used flag

      PetscReal            :: w_readtime

      PetscInt :: w_ncompe ! number of components+energy in flows

      ! cmpl flows
      PetscReal, pointer :: w_cmplflows(:,:) => null()
      ! flag indicating cmpl flows allocated
      PetscBool          :: w_cmplflows_Allocated

      PetscReal, pointer :: w_cmplflowsX(:,:,:,:) => null()
      PetscBool          :: w_cmplflowsX_Allocated

      ! points to next link in list
      class(well_data_type), pointer, public :: next

  contains

  ! Type-bound procedures held by well_data_type

    procedure, public :: Init => WellDataInit
    procedure, public :: Read => WellDataRead
    procedure, public :: Clear => WellDataClear
    procedure, public :: GetCmplLocation => GetCmplLocationInList
    procedure, public :: GetCmplLocationG => GetCmplLocationGInList
    procedure, public :: GetCmplGlobalLoc => GetCmplGlobalLocInList
    procedure, public :: FillCmplData => FillCmplDataInList
    procedure, public :: GetType => GetTypeInList
    procedure, public :: GetTargets => GetTargetsInList
    procedure, public :: SetActuals => SetActualsInList
    procedure, public :: ZeroActuals => ZeroActualsInList
    procedure, public :: DoUpdate => DoUpdateInList
    procedure, public :: DoIncrJac => DoIncrJacInList
    procedure, public :: GetNCmpl => GetNCmplInList
    procedure, public :: GetNCmplG => GetNCmplGInList
    procedure, public :: GetWellComm => GetCommInList
    procedure, public :: GetCmplDrillingDirection => &
                         GetCmplDrillingDirectionInList
    procedure, public :: GetCmplRadius => GetCmplRadiusInList
    procedure, public :: GetCmplThetaFactor => GetCmplThetaFactorInList
    procedure, public :: GetCmplSkinFactor => GetCmplSkinFactorInList
    procedure, public :: GetCmplDx => GetCmplDxInList
    procedure, public :: GetCmplDy => GetCmplDyInList
    procedure, public :: GetCmplDz => GetCmplDzInList
    procedure, public :: GetWellName => GetNameInList
    procedure, public :: GetWellTTVal => GetWellTTValInList
    procedure, public :: SetCCF => SetCCFInList
    procedure, public :: SetCmplZ => SetZInList
    procedure, public :: GetCCF => GetCCFInList
    procedure, public :: GetCmplZ => GetZInList
    procedure, public :: GetTT => GetTTInList
    procedure, public :: SetTT => SetTTInList
    procedure, public :: SetNcomp => SetNCompInList
    procedure, public :: SetCmplFlows => SetCmplFlowsInList
    procedure, public :: ZeroCmplFlows => ZeroCmplFlowsInList
    procedure, public :: MarkCmplForDeletion => MarkCmplForDeletionInList
    procedure, public :: DeleteMarkedCompletions => &
                         DeleteMarkedCompletionsInList
    procedure, public :: SetZRef => SetZRefInList
    procedure, public :: GetZRef => GetZRefInList
    procedure, public :: SetWellSolution => SetWellSolutionInList
    procedure, public :: GetWellSolution => GetWellSolutionInList
    procedure, public :: SetWellSolutionSet => SetWellSolutionSetInList
    procedure, public :: GetWellSolutionSet => GetWellSolutionSetInList
    procedure, public :: GetWellInjectionPAndT => GetWellInjectionPAndTInList
    procedure, public :: GetWellStatus => GetWellStatusInList

  end type well_data_type

  ! And now the actual list type

  type, public :: well_data_list_type
    PetscInt  :: num_well
    PetscReal :: f_actualP(N_WELL_TT) ! group prod rates
    PetscReal :: f_actualI(N_WELL_TT) ! group inj  totals
    PetscReal :: f_totalP (N_WELL_TT) ! group prod rates
    PetscReal :: f_totalI (N_WELL_TT) ! group inj  totals
    class(well_data_type), pointer :: first
    class(well_data_type), pointer :: last
    class(well_data_type), pointer :: array(:)
  end type well_data_list_type

  ! The public interface to this module

  public ::  WellDataSetFlag, WellDataGetFlag, &
             WellDataCreate, WellDataInitList, WellDataDestroyList, &
             WellDataAddToList, GetWellNCmpl, &
             WellSetGlobalInfo, WellSetGlobalInfoSet, &
             GetTargetUnitType, getnwell, getWellNameI, &
             GetWellTTValI, GetWellTypeI, &
             GetFieldData, SetFieldData, GetWellNCmplGI, GetCmplGlobalLocI, &
             FindGroupRates, FindGroupTotals, GetFieldTTVal, &
             IncrementWellWarningCount

  ! Private routines within this module

  private :: GetNCmplInList, &
             GetTargetsInList, SetActualsInList, GetTypeInList, &
             GetTTInList, SetTTInList, &
             DoUpdateInList, &
             GetCmplLocationInlist, GetCmplLocationGInlist, &
             FillCmplDataInList, GetCmplRadiusInList, &
             GetCmplThetaFactorInList, GetCmplSkinFactorInList,  &
             GetCmplDxInList, GetCmplDYInList, GetCmplDZInList, &
             SetCCFInList, GetCCFInList, &
             SetZInList, GetZInList, &
             SetWellSolutionInList, GetWellSolutionInList   , &
             SetWellSolutionSetInList, GetWellSolutionSetInList, &
             FindWellInList, GetWellInjectionPAndTInList

contains

! *************************************************************************** !

subroutine WellDataSetFlag()
  !
  ! Set flag indicating that well_data wells are used in this run
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  use_wells = PETSC_TRUE

end subroutine WellDataSetFlag

!*****************************************************************************!

function WellDataGetFlag()
  !
  ! Return status of flag indicating that well_data wells are used in this run
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscBool :: WellDataGetFlag

  WellDataGetFlag = use_wells

end function WellDataGetFlag

! *************************************************************************** !

function WellDataCreate()
  !
  ! Create a well_data instance and return a pointer to it
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type), pointer :: WellDataCreate
  class(well_data_type), pointer :: well_data

  allocate(well_data)
  call well_data%Init()

  WellDataCreate => well_data

end function WellDataCreate

! *************************************************************************** !

subroutine WellDataInit(this)
  !
  ! Initialise the data in a well_data_type instance
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this

  this%w_name  = ""
  this%w_index = 0
  this%w_itype = 0

  this%w_radius_set      = PETSC_FALSE
  this%w_skin_factor_set = PETSC_FALSE
  this%w_theta_frac_set  = PETSC_FALSE

  this%w_radius          = w_vrdef
  this%w_skin_factor     = 0.0
  this%w_theta_frac      = 1.0

  this%w_injection_t     = 15.0
  this%w_injection_p     = 1.01325d5

  this%w_pw = 1.01325d5
  this%w_pb = this%w_pw/2.0
  allocate(this%w_sp(n_phase));this%w_sp = 0.0
  this%w_issat  = PETSC_FALSE
  this%w_trel   = 0.0
  this%well_solution_set = PETSC_FALSE

  this%w_TT        = W_BHP_LIMIT
  this%w_status    = W_STATUS_OPEN

  this%w_injection_p = 1.01325d5
  this%w_injection_t = 15.0

  this%c_ci         => null()
  this%c_cj         => null()
  this%c_ck         => null()

  this%c_local_id    => null()
  this%c_ghosted_id  => null()
  this%c_to_cg       => null()

  this%cg_ghosted_id => null()

  this%c_onproc      => null()

  this%c_dx          => null()
  this%c_dy          => null()
  this%c_dz          => null()

  this%c_radius      => null()
  this%c_skin_factor => null()
  this%c_theta_frac  => null()
  this%c_drill_dir   => null()
  this%c_z           => null()

  this%c_ccf         => null()

  allocate(this%w_targets (N_WELL_TT));this%w_targets  = -1.0
  allocate(this%w_actuals (N_WELL_TT));this%w_actuals  =  0.0
  allocate(this%w_actualsG(N_WELL_TT));this%w_actualsG =  0.0
  allocate(this%w_totalsPG(N_WELL_TT));this%w_totalsPG =  0.0
  allocate(this%w_totalsIG(N_WELL_TT));this%w_totalsIG =  0.0

  this%w_z_ref_set = PETSC_FALSE
  this%w_z_ref     = 0.0

  this%w_const_drill_dir_set = PETSC_FALSE
  this%w_const_drill_dir     = Z_DIRECTION

  this%w_nrankw = 0

  this%w_ncmpl  = 0
  this%w_mcmpl  = 0
  this%w_ncmplg = 0
  this%w_mcmplg = 0

  this%w_ismp  = PETSC_FALSE
  this%w_group = 0
  this%w_comm  = 0

  this%w_nevent   = 0
  this%w_mevent   = 0
  this%w_readtime = 0.0

  this%w_ncompe = 0

  this%w_cmplflows => null()
  this%w_cmplflows_Allocated = PETSC_FALSE

  this%w_cmplflowsX => null()
  this%w_cmplflowsX_Allocated = PETSC_FALSE

  nullify(this%next)

end subroutine WellDataInit

! ************************************************************************** !

subroutine WellDataInitList(list, nphase)
  !
  ! Set up the initial state of a newly created well data list
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_list_type) :: list
  PetscInt, intent(in) :: nphase

  n_phase = nphase

  nullify(list%first)
  nullify(list%last)
  nullify(list%array)

  list%num_well = 0

  list%f_actualP = 0.0
  list%f_actualI = 0.0
  list%f_totalP  = 0.0
  list%f_totalI  = 0.0

end subroutine WellDataInitList

! ************************************************************************** !

subroutine WellDataRead(this, input, option, waytime, nwaytime, mwaytime)
  !
  ! Reads the data for a well from a WELL_DATA block in the input file
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  use Input_Aux_module
  use String_module
  use Option_module
  use Units_module
  use Utility_module, only : ReallocateArray

  implicit none

  type(option_type) :: option
  class(well_data_type) :: this
  type(input_type), pointer :: input
  PetscReal, pointer :: waytime(:)
  PetscInt, intent(inout) :: nwaytime, mwaytime

  character(len = MAXWORDLENGTH) :: keyword, word, units
  character(len = MAXWORDLENGTH) :: internal_units
  PetscInt :: ci, cj, ckl, cku, ival
  PetscReal :: v = 0.0

  ! Initialise

  internal_units = 'not_assigned'
  ci  = 1
  cj  = 1
  ckl = 1
  cku = 1

  input%ierr = 0

  ! Start reading data items
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input, option)
    if (InputError(input)) exit
    if (InputCheckExit(input, option)) exit

    call InputReadCard(input, option, keyword)
    call InputErrorMsg(input, option, 'keyword', 'WELL_DATA')
    call StringToUpper(keyword)

    select case(trim(keyword))
      case('RADIUS')  ! this has a conversion factor to account for
        ! Read value
        v = 0.0
        call InputReadDouble(input, option, v)
        call InputErrorMsg(input, option, 'well-radius' , 'WELL_DATA')
        call InputReadWord(input, option, word, PETSC_TRUE)
        ! Sort out units and convert if required
        internal_units = 'meter'
        if (InputError(input)) then
          word = trim(keyword) // ' UNITS'
          call InputDefaultMsg(input, option, word)
        else
          units = trim(word)
          v = v*UnitsConvertToInternal(units, internal_units, option)
        endif
        ! Store value
        this%w_radius     = v
        this%w_radius_set = PETSC_TRUE
      case('SKIN_FACTOR')
        call InputReadDouble(input, option, this%w_skin_factor)
        this%w_skin_factor_set = PETSC_TRUE
        call InputErrorMsg(input, option, 'skin factor', 'WELL_DATA')
        call InputReadWord(input, option, word, PETSC_TRUE)
      case('THETA_FRACTION')
        call InputReadDouble(input, option, this%w_theta_frac)
        this%w_theta_frac_set = PETSC_TRUE
        call InputErrorMsg(input, option, 'theta angle fraction', 'WELL_DATA')
        call InputReadWord(input, option, word, PETSC_TRUE)
      case('INJECTION_ENTHALPY_P')
        !  Read value
        v = 101325.0
        call InputReadDouble(input, option, v)
        call InputErrorMsg(input, option, &
                           'Injection enthalpy pressure', 'WELL_DATA')
        call InputReadWord(input, option, word, PETSC_TRUE)
        ! Sort out units and convert if required
        internal_units = 'Pa'
        if (InputError(input)) then
          word = trim(keyword) // ' UNITS'
          call InputDefaultMsg(input, option, word)
        else
          units = trim(word)
          v = v*UnitsConvertToInternal(units, internal_units, option)
        endif
        ! Store value
        if (this%w_readtime > 0.0) then
          call StoreEvent(this, EVENT_PINJ, ival, v)
        else
          this%w_injection_p = v
        endif
      case('INJECTION_ENTHALPY_T')
        v = 15.0
        call InputReadDouble(input, option, v)
        call InputErrorMsg(input, option, 'injection enthalpy temp', &
                           'WELL_DATA')
        call InputReadWord(input, option, word, PETSC_TRUE)
        if (this%w_readtime>0.0) then
          call StoreEvent(this, EVENT_TINJ, ival, v)
        else
          this%w_injection_t = v
        endif
      case('WELL_TYPE')
        call InputReadCard(input, option,  keyword , PETSC_FALSE)
        call InputErrorMsg(input, option, 'keyword', 'WELL_TYPE')
        call StringToUpper(keyword)

        select case(trim(keyword))
          case('PRODUCER', 'PROD')
            call SetWellType(this, PROD_WELL_TYPE)
          case('OIL_INJECTOR', 'OIL_INJ')
            call SetWellType(this, OIL_INJ_WELL_TYPE)
          case('GAS_INJECTOR', 'GAS_INJ')
            call SetWellType(this, GAS_INJ_WELL_TYPE)
          case('WATER_INJECTOR', 'WAT_INJ')
            call SetWellType(this, WAT_INJ_WELL_TYPE)
          case('SOLVENT_INJECTOR', 'SLV_INJ')
            call SetWellType(this, SLV_INJ_WELL_TYPE)
          case default
            option%io_buffer = 'WELL_DATA keyword: ' &
                               // trim(keyword) // ' not recognized'
            call PrintErrMsg(option)
        end select
      case('CONST_DRILL_DIR')
        this%w_const_drill_dir_set = PETSC_TRUE
        call InputReadCard(input, option, keyword, PETSC_FALSE)
        call InputErrorMsg(input, option, 'keyword', 'CONST_DRILL_DIR')
        call StringToUpper(keyword)

        select case(trim(keyword))
          case('DIR_X')
            this%w_const_drill_dir = X_DIRECTION
          case('DIR_Y')
            this%w_const_drill_dir = Y_DIRECTION
          case('DIR_Z')
            this%w_const_drill_dir = Z_DIRECTION
          case default
            option%io_buffer = &
            'Drilling direction ' // trim(keyword) // ' not recognized'
            call PrintErrMsg(option)
        end select
      case('Z_REF','D_REF')
        ! Read value
        v = 0.0
        call InputReadDouble(input, option, v)
        call InputErrorMsg(input, option, 'z_ref', 'WELL_DATA')
        call InputReadWord(input, option, word, PETSC_TRUE)
        ! Sort out units and convert if required
        internal_units = 'meter'
        if (InputError(input)) then
          word = trim(keyword) // ' UNITS'
          call InputDefaultMsg(input, option, word)
        else
          units = trim(word)
          v = v*UnitsConvertToInternal(units, internal_units, option)
        endif
        ! Case of DREF - convert depth to elevation
        if( trim(keyword) == 'D_REF' ) v = -v
        ! Store value
        this%w_z_ref     = v
        this%w_z_ref_set = PETSC_TRUE
      case('CIJK', 'CIJK_Z', 'CIJK_D')
      case('BHPL')
  ! Read a well bhp limit (will be max for injector, min for producer)
        v = 101325.0
        call InputReadDouble(input, option, v)
        call InputErrorMsg(input, option, 'BHPL', 'WELL_DATA')
        call InputReadWord(input, option, word, PETSC_TRUE)
        internal_units = 'Pa'
        if (InputError(input)) then
          word = trim(keyword) // ' UNITS'
          call InputDefaultMsg(input, option, word)
        else
          units = trim(word)
          v = v*UnitsConvertToInternal(units, internal_units, option)
        endif
        if (this%w_readtime > 0.0) then
          call StoreEvent(this, EVENT_TARG, W_BHP_LIMIT, v)
        else
          this%w_targets(W_BHP_LIMIT) = v
        endif
      case('OPEN')
          if (this%w_readtime>0.0) then
            call StoreEvent(this, EVENT_OPEN, W_BHP_LIMIT, v)
          else
            this%w_status    = W_STATUS_OPEN
          endif
      case('SHUT')
          if (this%w_readtime>0.0) then
            call StoreEvent(this, EVENT_SHUT, W_BHP_LIMIT, v)
          else
            this%w_status    = W_STATUS_SHUT
          endif
      case('TIME')
        ! Read the time associated with susbsequent instructions
        v = 0.0
        call InputReadDouble(input, option, v)
        call InputErrorMsg(input, option, 'TIME', 'WELL_DATA')
        call InputReadWord(input, option, word, PETSC_TRUE)
        ! Sort out units and convert if required
        internal_units = 'sec'
        if (InputError(input)) then
          word = trim(keyword) // ' UNITS'
          call InputDefaultMsg(input, option, word)
        else
          units = trim(word)
          v = v*UnitsConvertToInternal(units, internal_units, option)
        endif
        ! Store value
        this%w_readtime = v
        if (nwaytime+1 > mwaytime) call ReallocateArray(waytime, mwaytime)
        waytime(nwaytime+1) = v
        nwaytime = nwaytime + 1
  ! Read a surface volume rate (several options)
      case('TARG_OSV')
        call readWellTarget(this, input, option, 'TARG_OSV', word, W_TARG_OSV)
      case('TARG_GSV')
        call readWellTarget(this, input, option, 'TARG_GSV', word, W_TARG_GSV)
      case('TARG_WSV')
        call readWellTarget(this, input, option, 'TARG_WSV', word, W_TARG_WSV)
      case('TARG_SSV')
        call readWellTarget(this, input, option, 'TARG_SSV', word, W_TARG_SSV)
      case('TARG_LSV')
        call readWellTarget(this, input, option, 'TARG_LSV', word, W_TARG_LSV)
  ! Read a mass rate (several options)
      case('TARG_OM' )
        call readWellTarget(this, input, option, 'TARG_OM' , word, W_TARG_OM)
      case('TARG_GM' )
        call readWellTarget(this, input, option, 'TARG_GM' , word, W_TARG_GM)
      case('TARG_WM' )
        call readWellTarget(this, input, option, 'TARG_WM' , word, W_TARG_WM)
      case('TARG_SM' )
        call readWellTarget(this, input, option, 'TARG_SM' , word, W_TARG_SM)
      case('TARG_RV')
        call readWellTarget(this, input, option, 'TARG_RV' , word, W_TARG_RV)
    case default
        option%io_buffer = 'WELL_DATA keyword: ' &
                           // trim(keyword) // ' not recognized'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine WellDataRead

! ************************************************************************** !

subroutine WellDataAddToList(new_well_data, list)
  !
  ! Add a new well_data_type item to the list (ie data for a new well)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type), pointer :: new_well_data
  type(well_data_list_type), pointer :: list

  list%num_well   = list%num_well + 1
  new_well_data%w_index = list%num_well
  if (.not.associated(list%first)) list%first => new_well_data
  if (associated(list%last)) list%last%next => new_well_data
  list%last => new_well_data

end subroutine WellDataAddToList

! ************************************************************************** !

subroutine WellDataDestroyList(well_data_list, option)
  !
  ! De-allocates a whole list of well_data_type items
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  use Option_module

  implicit none

  type(well_data_list_type), pointer :: well_data_list
  type(option_type) :: option

  class(well_data_type), pointer :: well_data, prev_well_data

  PetscInt :: nwarn, unconv_tg, unconv_wg, unconv_bg, mpierr_g
  PetscMPIInt :: ibufl(4), ibufg(4), ierr

  ! If wells used, check for warnings

  if (use_wells) then

    ! Globalise the warnings

    ibufl(1) = w_unconv_t
    ibufl(2) = w_unconv_w
    ibufl(3) = w_unconv_b
    ibufl(4) = w_mpierr
    ibufg    = 0
    ierr     = 0

    call MPI_Reduce(ibufl, ibufg, FOUR_INTEGER_MPI, &
                    MPI_INTEGER, MPI_SUM, &
                    option%io_rank, option%mycomm, ierr)

    unconv_tg = ibufg(1)
    unconv_wg = ibufg(2)
    unconv_bg = ibufg(3)
    mpierr_g  = ibufg(4)

    ! Report warnings

    if (option%io_rank == option%myrank) then
      nwarn = unconv_tg+unconv_wg+unconv_bg+mpierr_g
      if (nwarn>0) then
        print *, 'Well model convergence failure counts:'
        if (unconv_tg > 0) print *, unconv_tg, ' mode selection'
        if (unconv_wg > 0) print *, unconv_wg, ' bhp solution'
        if (unconv_bg > 0) print *, unconv_bg, ' wellbore composition'
        if (mpierr_g  > 0) print *, mpierr_g , ' well solver MPI errors'
      endif
    endif

  endif

  ! Skip if list not allocated

  if (.not.associated(well_data_list)) return

  ! Start at the first item and loop through list destroying items

  well_data => well_data_list%first
  do
    if (.not.associated(well_data)) exit
    prev_well_data => well_data
    well_data => well_data%next
    call WellDataDestroy(prev_well_data)
  enddo

  ! Nullify and deallocate the list itself

  well_data_list%num_well = 0
  nullify(well_data_list%first)
  nullify(well_data_list%last)
  if (associated(well_data_list%array)) deallocate(well_data_list%array)
  nullify(well_data_list%array)

  deallocate(well_data_list)
  nullify(well_data_list)

end subroutine WellDataDestroyList

! ************************************************************************** !

subroutine WellDataDestroy(well_data)
  !
  ! De-allocates a single well_data_type item
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type), pointer :: well_data

  ! Skip if item not allocated

  if (.not.associated(well_data)) return

  ! Clear out and deallocate

  call well_data%Clear()

  deallocate(well_data)
  nullify(well_data)

end subroutine WellDataDestroy

! ************************************************************************** !

subroutine WellDataClear(this)
  !
  ! Clear out data contained in a single well_data_type item
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  use Utility_module, only : DeallocateArray

  implicit none

  class(well_data_type) :: this

  call DeallocateArray(this%w_targets )
  call DeallocateArray(this%w_actuals )
  call DeallocateArray(this%w_actualsG)
  call DeallocateArray(this%w_totalsPG)
  call DeallocateArray(this%w_totalsIG)

  call DeallocateArray(this%c_ci)
  call DeallocateArray(this%c_cj)
  call DeallocateArray(this%c_ck)

  call DeallocateArray(this%c_local_id  )
  call DeallocateArray(this%c_ghosted_id)
  call DeallocateArray(this%c_to_cg     )

  call DeallocateArray(this%cg_ghosted_id)

  call DeallocateArray(this%c_onproc)

  call DeallocateArray(this%c_dx)
  call DeallocateArray(this%c_dy)
  call DeallocateArray(this%c_dz)

  call DeallocateArray(this%c_radius     )
  call DeallocateArray(this%c_skin_factor)
  call DeallocateArray(this%c_theta_frac )
  call DeallocateArray(this%c_drill_dir  )
  call DeallocateArray(this%c_z          )

  call DeallocateArray(this%c_ccf        )

  call DeallocateArray(this%c_mfd        )

  call DeallocateArray(this%w_event_code)
  call DeallocateArray(this%w_event_time)
  call DeallocateArray(this%w_event_ival)
  call DeallocateArray(this%w_event_rval)
  call DeallocateArray(this%w_event_used)

  if (this%w_cmplFlows_allocated) then
    call DeallocateArray(this%w_cmplFlows)
    this%w_cmplFlows_allocated = PETSC_FALSE
  endif

  if (this%w_cmplFlowsX_allocated) then
    call DeallocateArray(this%w_cmplFlowsX)
    this%w_cmplFlowsX_allocated = PETSC_FALSE
  endif

  nullify(this%next)

end subroutine WellDataClear

! ************************************************************************** !

function GetWellNCmpl(iwell, list)
  !
  ! Get the number of completions (on this proc) on a given well
  ! iwell - in - The index of the well (1..m_numwell)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt :: GetWellNCmpl
  PetscInt, intent(in) :: iwell
  type(well_data_list_type), pointer :: list
  class(well_data_type), pointer :: well_data

  PetscBool :: found

  GetWellNCmpl = 0

  ! Find well with index iwell
  found = FindWellInList(iwell, well_data, list)

  ! If found, return number of completions
  if (found) then
    GetWellNCmpl = well_data%w_ncmpl
  endif

end function GetWellNCmpl

! ************************************************************************** !

function GetWellNCmplGI(iwell, list)
  !
  ! Get the number of completions (on this proc) on a given well
  ! iwell - in - The index of the well (1..m_numwell)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt :: GetWellNCmplGI
  PetscInt, intent(in) :: iwell
  type(well_data_list_type), pointer :: list
  class(well_data_type), pointer :: well_data

  PetscBool :: found

  GetWellNCmplGI = 0

  ! Find well with index iwell
  found = FindWellInList(iwell, well_data, list)

  ! If found, return number of completions
  if (found) then
    GetWellNCmplGI = well_data%w_ncmplg
  endif

end function GetWellNCmplGI

! ************************************************************************** !

subroutine getWellNameI(iwell, list, name)
  !
  ! Get the name of the (iwell)th well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt, intent(in) :: iwell
  type(well_data_list_type), pointer :: list
  character(len = MAXSTRINGLENGTH), intent(out) :: name

  class(well_data_type), pointer :: well_data
  PetscBool :: found

  name = ''

  ! Find well with index iwell
  found = FindWellInList(iwell, well_data, list)

  ! If found, return name (otherwise leave at default)
  if (found) then
    call well_data%GetWellName(name)
  endif

end subroutine getWellNameI

! ************************************************************************** !

function getWellTypeI(iwell, list)
  !
  ! Get the type of the (iwell)th well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt, intent(in) :: iwell
  type(well_data_list_type), pointer :: list
  PetscInt :: getWellTypeI

  class(well_data_type), pointer :: well_data
  PetscBool :: found

  ! Default to producer
  getWellTypeI = PROD_WELL_TYPE

  ! Find well with index iwell
  found = FindWellInList(iwell, well_data, list)

  ! If found, return name (otherwise leave at default)
  if (found) then
    getWellTypeI = well_data%GetType()
  endif

end function getWellTypeI

! ************************************************************************** !

subroutine GetCmplGlobalLocI(iwell, icmplG, ci, cj, ck, cdd, list)
  !
  ! Get the location for the icmplG'th completion of the iw'th well
  ! Note this is the global completion count, not that on the proc
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt, intent(in) :: iwell
  PetscInt, intent(in) :: icmplG
  PetscInt, intent(out) :: ci, cj, ck, cdd
  type(well_data_list_type), pointer :: list
  class(well_data_type), pointer :: well_data
  PetscBool :: found

  ! Find well with index iwell
  found = FindWellInList(iwell, well_data, list)

  ! If found, return completion location
  if (found) then
    call well_data%GetCmplGlobalLoc(icmplG, ci, cj, ck, cdd)
  else
    ci  = 1
    cj  = 1
    ck  = icmplG
    cdd = 3
  endif

end subroutine GetCmplGlobalLocI

! ************************************************************************** !

subroutine GetCmplGlobalLocInList(this, icmplG, ci, cj, ck, cdd)
  !
  ! Given the list elemetn for a well, get the completion details
  ! Note this is the global completion index, not that on the proc
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this
  PetscInt, intent(in)  :: icmplG
  PetscInt, intent(out) :: ci, cj, ck, cdd

  ci  = 1
  cj  = 1
  ck  = 1
  cdd = 3

  if (icmplG .le. this%w_ncmplg) then
    ci  = this%c_ci(icmplG)
    cj  = this%c_cj(icmplG)
    ck  = this%c_ck(icmplG)
    cdd = this%c_drill_dir(icmplG)
  endif

end subroutine GetCmplGlobalLocInList

! ************************************************************************** !

function getFieldTTVal(itt, valType, list)
  !
  ! Get value of the (itt)th target type of the (0)th group
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal :: getFieldTTVal
  PetscInt, intent(in) :: itt
  PetscInt, intent(in) :: valType
  type(well_data_list_type), pointer :: list

  getFieldTTVal = 0.0

  if (valType == VALTYPE_ACTUALP) getFieldTTVal = list%f_actualP(itt)
  if (valType == VALTYPE_ACTUALI) getFieldTTVal = list%f_actualI(itt)
  if (valType == VALTYPE_TOTALP ) getFieldTTVal = list%f_totalP (itt)
  if (valType == VALTYPE_TOTALI ) getFieldTTVal = list%f_totalI (itt)

end function getFieldTTVal

! ************************************************************************** !

function GetWellTTValI(iwell, itt, valType, list)
  !
  ! Get value of the (itt)th target type of the (iwell)th well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal :: getWellTTValI
  PetscInt, intent(in) :: iwell
  PetscInt, intent(in) :: itt
  PetscInt, intent(in) :: valType
  type(well_data_list_type), pointer :: list

  class(well_data_type), pointer :: well_data

  PetscBool :: found

  GetWellTTValI = 0.0

  ! Find well with index iwell
  found = FindWellInList(iwell, well_data, list)

  ! Find value of target type itt
  if (found) then
    GetWellTTValI = well_data%GetWellTTVal(itt, valType)
  endif

end function GetWellTTValI

! ************************************************************************** !

function GetWellTTValInList(this, itt, valType)
  !
  ! Get value of the (itt)th target type for a given well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal :: GetWellTTValInList
  class(well_data_type) :: this
  PetscInt, intent(in) :: itt
  PetscInt, intent(in) :: valType

  ! May be the required value, the actual value or the time-total value
  GetWellTTValInList = 0.0
  if (valType == VALTYPE_TARGET) GetWellTTValInList = this%w_targets (itt)
  if (valType == VALTYPE_ACTUAL) GetWellTTValInList = this%w_actualsG(itt)
  if (valType == VALTYPE_TOTALP) GetWellTTValInList = this%w_totalsPG(itt)
  if (valType == VALTYPE_TOTALI) GetWellTTValInList = this%w_totalsIG(itt)

end function GetWellTTValInList

! ************************************************************************** !

subroutine WellSetGlobalInfo(iwell, nrankw, ncmplg, ismp, group, comm, list)
  !
  ! Set a package of useful global information for the (iwell) the well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !
  ! nrankw - number of ranks on which this well appears
  ! ncmplg - number of global completions (ie over all procs)
  ! ismp   - bool indicating this is a cross-processor well
  ! group  - MPI group index
  ! comm   - MPI communicator index

  implicit none

  PetscInt   , intent(in) :: iwell, ncmplg
  PetscBool  , intent(in) :: ismp
  MPI_Group, intent(in) :: group
  MPI_Comm , intent(in) :: comm
  type(well_data_list_type), pointer :: list

  class(well_data_type), pointer :: well_data
  PetscBool :: found
  PetscInt :: nrankw

  ! Find well in list
  found = FindWellInList(iwell, well_data, list)

  ! Set values if found
  if (found) then
    well_data%w_nrankw = nrankw
    well_data%w_ncmplg = ncmplg
    well_data%w_ismp   = ismp
    well_data%w_group  = group
    well_data%w_comm   = comm
  else
    call throwWellDataException('Unable to find well')
  endif

end subroutine WellSetGlobalInfo

! ************************************************************************** !

subroutine WellSetGlobalInfoSet
  !
  ! Set flag indicating that global information has been set
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  global_info_set = PETSC_TRUE

end subroutine WellSetGlobalInfoSet

! ************************************************************************** !

subroutine GetCmplLocationInList(this, icmpl, local_id, &
                                 ghosted_id, onproc, icmplg)
  !
  ! Get the well location (structured grid only)
  ! Error will occur if unstructured grid used with this code
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this

  PetscInt , intent(in ) :: icmpl
  PetscInt , intent(out) :: local_id
  PetscInt , intent(out) :: ghosted_id
  PetscBool, intent(out) :: onproc
  PetscInt , intent(out) :: icmplg

  PetscInt :: ncmpl

  ncmpl = this%GetNCmpl()

  if (icmpl > 0 .and. icmpl <= ncmpl) then
    local_id   = this%c_local_id  (icmpl)
    ghosted_id = this%c_ghosted_id(icmpl)
    onproc     = this%c_onproc    (icmpl)
    icmplg     = this%c_to_cg     (icmpl)
  endif

end subroutine GetCmplLocationInList

! ************************************************************************** !

subroutine GetCmplLocationGInList(this, icmplg, ghosted_id)
  !
  ! Get the well location wrt global completion index
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this

  PetscInt, intent(in ) :: icmplg
  PetscInt, intent(out) :: ghosted_id

  PetscInt :: ncmplg

  ncmplg = this%GetNCmplG()

  if (icmplg > 0 .and. icmplg <= ncmplg) then
    ghosted_id = this%cg_ghosted_id(icmplg)
  endif

end subroutine GetCmplLocationGInList

! ************************************************************************** !

subroutine FillCmplDataInList(this, iw, grid)
  !
  ! Fill in any user-supplied completion data
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  use Grid_Grdecl_module, only : GetGrdNCmpl, GetCmplData
  use Grid_module

  implicit none

  class(well_data_type) :: this
  PetscInt, intent(in) :: iw
  type(grid_type), pointer :: grid

  PetscInt  :: icmpl, ci, cj, ck, ia, local_id, ghosted_id
  PetscReal :: dx, dy, dz, z
  PetscBool :: onproc

  local_id   = 1
  ghosted_id = 1

  this%w_ncmpl = GetGrdNCmpl(iw)

  call CheckCompletionCount(this)

  do icmpl = 1, this%w_ncmpl

    call GetCmplData(iw, icmpl, ci, cj, ck, ia, dx, dy, dz, z)

    call getLocalAndGhostedIDs(ia, local_id, ghosted_id, onproc, grid)

    this%c_ci        (icmpl) = ci
    this%c_cj        (icmpl) = cj
    this%c_ck        (icmpl) = ck
    this%c_local_id  (icmpl) = local_id
    this%c_ghosted_id(icmpl) = ghosted_id
    this%c_onproc    (icmpl) = onproc

    if (this%w_radius_set         ) &
      this%c_radius     (icmpl) = this%w_radius

    if (this%w_skin_factor_set    ) &
      this%c_skin_factor(icmpl) = this%w_skin_factor

    if (this%w_theta_frac_set     ) &
      this%c_theta_frac (icmpl) = this%w_theta_frac

    if (this%w_const_drill_dir_set) &
      this%c_drill_dir  (icmpl) = this%w_const_drill_dir

    this%c_dx        (icmpl)  = dx
    this%c_dy        (icmpl)  = dy
    this%c_dz        (icmpl)  = dz

    this%c_z         (icmpl)  = z

  enddo

end subroutine FillCmplDataInList

! ************************************************************************** !

subroutine SetZRefInList(this, option)
  !
  ! Set the reference elevation for this well (if the user has not set a value)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  use Option_module

  implicit none

  class(well_data_type) :: this
  type(option_type) :: option

  PetscInt :: icmpl, ierr
  PetscReal :: ztop, zcmpl
  PetscReal :: zl(1)
  PetscReal :: zg(1)

  zl = 0.0
  zg = 0.0

  ! Check that not user-set

  if (.not.(this%w_z_ref_set)) then

  ! Set to very low elevation and reset to max elevation
  ! of completion on this proc (may be none)

    ztop = -1.0D6
    do icmpl = 1, this%w_ncmpl
      zcmpl = this%c_Z(icmpl)
      if (zcmpl > ztop) ztop = zcmpl
    enddo

  ! Take global max to find final value

    zl(1) = ztop
    call MPI_Allreduce(zl, zg, ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, MPI_MAX, option%mycomm, ierr)
    ztop = zg(1)

  ! If still original value (no completions anywhere) reset to zero

    if (ztop<-0.5D5) ztop = 0.0

  ! Store calculated value

    this%w_z_ref = ztop
  endif

end subroutine SetZRefInList

! ************************************************************************** !

function GetZRefInList(this)
  !
  ! Get the reference elevation for this well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal :: GetZRefInList
  class(well_data_type) :: this

  GetZRefInList = this%w_z_ref

end function GetZRefInList

! ************************************************************************** !

subroutine SetWellSolutionInList(this, pw, pb, sp, issat, trel)
  !
  ! Set the wellbore solution for this well (as obtained by well_solver)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this
  PetscReal, intent(in) :: pw
  PetscReal, intent(in) :: pb
  PetscReal, intent(in) :: sp(:)
  PetscBool, intent(in) :: issat
  PetscReal, intent(in) :: trel

  this%w_pw     = pw
  this%w_pb     = pb
  this%w_sp     = sp
  this%w_issat  = issat
  this%w_trel   = trel

end subroutine SetWellSolutionInList

! ************************************************************************** !

subroutine GetWellSolutionInList(this, pw, pb, sp, issat, trel)
  !
  ! Set the wellbore solution for this well (as predictor for well_solver)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this
  PetscReal, intent(out) :: pw
  PetscReal, intent(out) :: pb
  PetscReal, intent(out) :: sp(:)
  PetscBool, intent(out) :: issat
  PetscReal, intent(out) :: trel

  pw     = this%w_pw
  pb     = this%w_pb
  sp     = this%w_sp
  issat  = this%w_issat
  trel   = this%w_trel

end subroutine GetWellSolutionInList

! ************************************************************************** !

subroutine SetWellSolutionSetInList(this)
  !
  ! Set the flag indicating that a wellbore solution stored for this well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this

  this%well_solution_set = PETSC_TRUE

end subroutine SetWellSolutionSetInList

! ************************************************************************** !

function GetWellSolutionSetInList(this)
  !
  ! Get the flag indicating that a wellbore solution stored for this well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscBool :: GetWellSolutionSetInList
  class(well_data_type) :: this

  GetWellSolutionSetInList = this%well_solution_set

end function GetWellSolutionSetInList

! ************************************************************************** !

function GetTargetsInList(this, targets)
  !
  ! Get the current well targets
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscBool :: GetTargetsInList
  PetscInt  :: iTarget
  class(well_data_type) :: this
  PetscReal, intent(out) :: targets(N_WELL_TT)

  do iTarget = 1, N_WELL_TT
    targets(iTarget) = this%w_targets(iTarget)
  enddo
  GetTargetsInList = PETSC_TRUE

end function GetTargetsInList

! ************************************************************************** !

subroutine SetActualsInList(this, actuals)
  !
  ! Get the actual values for each well target
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt  :: iTarget
  class(well_data_type) :: this
  PetscReal, intent(in) :: actuals(N_WELL_TT)

  do iTarget = 1, N_WELL_TT
    this%w_actuals(iTarget) = actuals(iTarget)
  enddo

end subroutine SetActualsInList

! ************************************************************************** !

subroutine ZeroActualsInList(this)
  !
  ! Zero the actual values for each well target
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this

  this%w_actuals = 0.0

end subroutine ZeroActualsInList

! ************************************************************************** !

function GetTypeInList(this)
  !
  ! Get the current well type
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt :: GetTypeInList
  class(well_data_type) :: this

  GetTypeInList = this%w_itype

end function GetTypeInList

! ************************************************************************** !

subroutine SetTTInList(this, w_TT)
  !
  ! Set the well target type
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt :: w_TT
  class(well_data_type) :: this

  this%w_TT = w_TT

end subroutine SetTTInList

! ************************************************************************** !

subroutine SetNCompInList(this, ncomp)
  !
  ! Set the number of components in the problem
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this
  PetscInt, intent(in)   :: nComp

  this%w_ncompe = ncomp+1

end subroutine SetNCompInList

! ************************************************************************** !

subroutine SetCmplFlowsInList(this, cmplflows, cmplflowsX, &
                              ncmpl, ncompe, ncmplg, ndof, isothermal)
  !
  ! Set the completion flows
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this
  PetscReal             :: cmplflows(:,:)
  PetscReal             :: cmplflowsX(:,:,:,:)
  PetscInt, intent(in)  :: ncmpl
  PetscInt, intent(in)  :: ncompe
  PetscInt, intent(in)  :: ncmplg
  PetscInt, intent(in)  :: ndof
  PetscBool, intent(in) :: isothermal
  PetscInt              :: nicmplArg, ncompeArg, njcmplArg, ndofArg, &
                           icmpl, icompe, jcmpl, jdof, &
                           nicmpl, nicompe, njcmpl, njdof

  wd_isothermal = isothermal

  ! If the completion flows have not been allocated, allocate them now

  if (.not.(this%w_cmplflows_Allocated)) then
    allocate(this%w_cmplflows(ncmpl, ncompe))
    this%w_cmplflows_allocated = PETSC_TRUE
    this%w_cmplflows = 0.0
  endif

  ! If the completion flow derivatives have not been allocated, allocate now

  if (.not.(this%w_cmplflowsX_Allocated)) then
    allocate(this%w_cmplflowsX(ncmpl, ncompe, ncmplg, ndof))
    this%w_cmplflowsX_allocated = PETSC_TRUE
    this%w_cmplflowsX = 0.0
  endif

  ! Copy over the completion flows (use min. dim. of stored and arg arrays)

  nicmplArg = size(cmplflows, 1)
  ncompeArg = size(cmplflows, 2)

  if (this%w_cmplflows_allocated) then
    do icmpl = 1, min(nicmplArg, this%w_ncmpl)
      do icompe = 1, min(ncompeArg, this%w_ncompe)
        this%w_cmplflows(icmpl, icompe) = cmplflows(icmpl, icompe)
      enddo
    enddo
  endif

  ! Copy over the completion flow derivatives
  ! Use minimum dimension of stored and argument arrays

  nicmplArg = size(cmplflowsX, 1)
  ncompeArg = size(cmplflowsX, 2)
  njcmplArg = size(cmplflowsX, 3)
  ndofArg   = size(cmplflowsX, 4)

  if (this%w_cmplflowsX_allocated) then
    nicmpl  = min(nicmplArg, this%w_ncmpl )
    nicompe = min(ncompeArg, this%w_ncompe)
    njcmpl  = min(njcmplArg, this%w_ncmplg)
    njdof   = min(ndofArg, ndof)
    do icmpl = 1, nicmpl
      do icompe = 1, nicompe
        do jcmpl = 1, njcmpl
          do jdof = 1, njdof
             this%w_cmplflowsX(icmpl, icompe, jcmpl, jdof) = &
                    cmplflowsX(icmpl, icompe, jcmpl, jdof)
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine SetCmplFlowsInList

! ************************************************************************** !

subroutine ZeroCmplFlowsInList(this, ncmpl, ncompe, ncmplg, ndof)
  !
  ! Set the completion flows
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this
  PetscInt, intent(in)   :: ncmpl
  PetscInt, intent(in)   :: ncompe
  PetscInt, intent(in)   :: ncmplg
  PetscInt, intent(in)   :: ndof

  ! If the completion flows have not been allocated, allocate them now

  if (.not.(this%w_cmplflows_Allocated)) then
    allocate(this%w_cmplflows(ncmpl, ncompe))
    this%w_cmplflows_allocated = PETSC_TRUE
  endif

  ! If the completion flow derivatives have not been allocated, allocate now

  if (.not.(this%w_cmplflowsX_Allocated)) then
    allocate(this%w_cmplflowsX(ncmpl, ncompe, ncmplg, ndof))
    this%w_cmplflowsX_allocated = PETSC_TRUE
  endif

  ! Zero the completion flows

  if (this%w_cmplflows_allocated) then
   this%w_cmplflows = 0.0
  endif

  ! Zero the completion flow derivatives

  if (this%w_cmplflowsX_allocated) then
    this%w_cmplflowsX = 0.0
  endif

end subroutine ZeroCmplFlowsInList

! ************************************************************************** !

subroutine DoUpdateInList(this, dt, option)
  !
  ! Update this well at the end of a timestep
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  use Option_module

  implicit none

  class(well_data_type) :: this
  PetscReal, intent(in) :: dt
  type(option_type) :: option

  PetscInt :: itt, vmtype, ierr, ievent, event_code, itarg, itype
  PetscReal :: den, deni, time, event_time, vtarg

  ! Check that simulation has started

  if (global_info_set) then

  ! Globalise well flows over all procs-even if well has no cmpls on this proc

  ! Do collective sum back into actuals 
  ! Note we do all procs to sync all the wells on all the procs

    call MPI_Allreduce( this%w_actuals, this%w_actualsG, N_WELL_TT &
                       , MPI_DOUBLE_PRECISION, MPI_SUM, option%mycomm, ierr)

    den = this%w_nrankw
    deni = 0.0
    if (den>0.0) deni = 1.0/den

  ! Convert sum to average for pressures

    do itt = 1, N_WELL_TT
      vmtype = GetTargetUnitType(itt)
      if (vmtype == TT_P) then
        this%w_actualsG(itt) = this%w_actualsG(itt)*deni
      endif
    enddo

  ! Increment the well totals and zero actuals for next step

    do itt = 1, N_WELL_TT
      vmtype = GetTargetUnitType(itt)
      if (vmtype /= TT_P) then
        if (this%w_itype == PROD_WELL_TYPE) then
          this%w_totalsPG(itt) = this%w_totalsPG(itt)+this%w_actualsG(itt)*dt
        else
  ! Note conversion to unsigned
          this%w_totalsIG(itt) = this%w_totalsIG(itt)-this%w_actualsG(itt)*dt
        endif
      endif
    enddo

  endif

  ! Go through event list and extract any instructions

  time = option%time
  do ievent = 1, this%w_nevent
    event_time = this%w_event_time(ievent)
    if (.not.this%w_event_used(ievent)) then
      if (event_time<(time+w_event_time_eps)) then
        event_code = this%w_event_code(ievent)
        if (event_code == EVENT_TARG) then
          itarg = this%w_event_ival(ievent)
          vtarg = this%w_event_rval(ievent)
          this%w_targets(itarg) = vtarg
        endif
        if (event_code == EVENT_TYPE) then
          itype = this%w_event_ival(ievent)
          call UpdateWellType(this, itype)
        endif
        if (event_code == EVENT_OPEN) then
          this%w_status = W_STATUS_OPEN
        endif
        if (event_code == EVENT_SHUT) then
          this%w_status = W_STATUS_SHUT
        endif
        this%w_event_used(ievent) = PETSC_TRUE
      endif
    endif
  enddo

end subroutine DoUpdateInList

! ************************************************************************** !

subroutine DoIncrJacInList(this, option, nflowdof, Jup, A)
  !
  ! Increment Jacobian
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module

  implicit none

  class(well_data_type) :: this
  type(option_type) :: option
  PetscInt , intent(in ) :: nflowdof
  PetscReal, intent(out) :: Jup(nflowdof, nflowdof)
  Mat :: A

  PetscInt :: ghosted_id, ghosted_jd, ierr, icmpl, jcmplg, &
              ncmpl, ncmplg, i, j, status

  ncmpl  = this%GetNCmpl()
  ncmplg = this%GetNCmplG()

  status = this%GetWellStatus()

  if (ncmpl > 0 .and. status == W_STATUS_OPEN) then

    ierr = 0

  ! Loop over all connection terms d(flow(icmpl)/d(Xc(jcmpl))

    do icmpl = 1, ncmpl

      ghosted_id = this%c_ghosted_id(icmpl)

      do jcmplg = 1, ncmplg

        ghosted_jd = this%cg_ghosted_id(jcmplg)

       Jup = 0.0

        do i = 1, nflowdof
          do j = 1, nflowdof
            Jup(i, j) = this%w_cmplflowsX(icmpl, i, jcmplg, j)
          enddo
        enddo

        if (wd_isothermal) then
          Jup(option%energy_id, :) = 0.d0
          Jup(:, option%energy_id) = 0.d0
        endif

        call MatSetValuesBlockedLocal(A, 1, ghosted_id-1, &
                                         1, ghosted_jd-1, Jup, &
                                         ADD_VALUES, ierr);CHKERRQ(ierr)

      enddo
    enddo
  endif ! Well open and completed on this rank

end subroutine DoIncrJacInList

! ************************************************************************** !

function GetNCmplInList(this)
  !
  ! Returns number of completions for this well (on this proc)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt :: GetNCmplInList

  class(well_data_type) :: this

  GetNCmplInList = this%w_ncmpl

end function GetNCmplInList

! *************************************************************************** !

function GetNCmplGInList(this)
  !
  ! Returns number of completions for this well (on all procs)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt :: GetNCmplGInList

  class(well_data_type) :: this

  GetNCmplGInList = this%w_ncmplG

end function GetNCmplGInList

! *************************************************************************** !

function GetCommInList(this,ismp)
  !
  ! Returns MPI communicator index for this well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  MPI_Comm :: GetCommInList
  PetscBool,intent(out) :: ismp

  class(well_data_type) :: this

  GetCommInList = this%w_comm
  ismp          = this%w_ismp

end function GetCommInList

! *************************************************************************** !

function GetCmplDrillingDirectionInList(this, icmpl)
  !
  ! Get drilling direction for (icmpl) th completion on this well on this proc
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt :: GetCmplDrillingDirectionInList
  PetscInt, intent(in) :: icmpl

  class(well_data_type) :: this

  GetCmplDrillingDirectionInList = this%c_drill_dir(icmpl)

end function GetCmplDrillingDirectionInList

! *************************************************************************** !

function GetCmplRadiusInList(this, icmpl)
  !
  ! Get wellbore radius for (icmpl) th completion on this well on this proc
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal :: GetCmplRadiusInList
  PetscInt, intent(in) :: icmpl
  class(well_data_type) :: this

  GetCmplRadiusInList = this%c_radius(icmpl)

end function GetCmplRadiusInList

! *************************************************************************** !

function GetCmplSkinFactorInList(this, icmpl)
  !
  ! Get skin factor for (icmpl) th completion on this well on this proc
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal :: GetCmplSkinFactorInList
  PetscInt, intent(in) :: icmpl
  class(well_data_type) :: this

  GetCmplSkinFactorInList = this%c_skin_factor(icmpl)

end function GetCmplSkinFactorInList

! *************************************************************************** !

function GetCmplDxInList(this, icmpl)
  !
  ! Get cell dx size for (icmpl) th completion on this well on this proc
  ! This is usually -1, indicating that grid value should be used
  ! However, the value can be set by the user
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal             :: GetCmplDxInList
  class(well_data_type) :: this
  PetscInt, intent(in)  :: icmpl

  GetCmplDxInList = this%c_dx(icmpl)

end function GetCmplDxInList

! *************************************************************************** !

function GetCmplDyInList(this, icmpl)
  !
  ! Get cell dy size for (icmpl) th completion on this well on this proc
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal             :: GetCmplDyInList
  class(well_data_type) :: this
  PetscInt, intent(in)   :: icmpl

  GetCmplDyInList = this%c_dy(icmpl)

end function GetCmplDyInList

! *************************************************************************** !

function GetCmplDzInList(this, icmpl)
  !
  ! Get cell dz size for (icmpl) th completion on this well on this proc
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal             :: GetCmplDzInList
  class(well_data_type) :: this
  PetscInt, intent(in)  :: icmpl

  GetCmplDzInList = this%c_dz(icmpl)

end function GetCmplDzInList

! *************************************************************************** !

function GetCmplThetaFactorInList(this, icmpl)
  !
  ! Get theta factor for (icmpl) th completion on this well on this proc
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal             :: GetCmplThetaFactorInList
  class(well_data_type) :: this
  PetscInt, intent(in)  :: icmpl

  GetCmplThetaFactorInList = this%c_theta_frac(icmpl)

end function GetCmplThetaFactorInList

! *************************************************************************** !

subroutine GetNameInList(this, name)
  !
  ! Get name of this well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  character(len=MAXSTRINGLENGTH), intent(out) :: name
  class(well_data_type) :: this

  name = this%w_name

end subroutine GetNameInList

! *************************************************************************** !

subroutine SetCCFInList(this, icmpl, ccf)
  !
  ! Set completion connection factor for the (icmpl)th completion on this well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal :: ccf
  PetscInt, intent(in) :: icmpl
  class(well_data_type) :: this

  this%c_ccf(icmpl) = ccf

end subroutine SetCCFInList

! ************************************************************************** !

function GetCCFInList(this, icmpl)
  !
  ! Get completion connection factor for the (icmpl)th completion on this well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal :: GetCCFInList
  PetscInt, intent(in) :: icmpl
  class(well_data_type) :: this

  GetCCFInList = this%c_ccf(icmpl)

end function GetCCFInList

! *************************************************************************** !

subroutine SetZInList(this, icmpl, z)
  !
  ! Set elevation for the (icmpl)th completion on this well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal :: z
  PetscInt, intent(in) :: icmpl
  class(well_data_type) :: this

  this%c_z(icmpl) = z

end subroutine SetZInList

! *************************************************************************** !

function GetZInList(this, icmpl)
  !
  ! Get elevation for the (icmpl)th completion on this well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscReal :: GetZInList
  PetscInt, intent(in) :: icmpl
  class(well_data_type) :: this

  GetZInList = this%c_z(icmpl)

end function GetZInList

! *************************************************************************** !

function GetTTInList(this)
  !
  ! Get current target type this well
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt :: GetTTInList
  class(well_data_type) :: this

  GetTTInList = this%w_TT

end function GetTTInList

! *************************************************************************** !

subroutine GetWellInjectionPAndTInList(this, injection_p, injection_t)
  !
  ! Get current target injection pressure and temperature
  ! (for specific enthalpy calculation)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none
  class(well_data_type) :: this

  PetscReal, intent(out) :: injection_p
  PetscReal, intent(out) :: injection_t

  injection_p = this%w_injection_p
  injection_t = this%w_injection_t

end subroutine GetWellInjectionPAndTInList

! *************************************************************************** !

subroutine readWellTarget(this, input, option, keyword, word, target_type)
  !
  ! Read a general well target from the input file
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  use Input_Aux_module
  use String_module
  use Option_module
  use Units_module

  implicit none

  class(well_data_type) :: this
  type(option_type)            :: option
  type(input_type), pointer    :: input
  character(len=*), intent(in) :: keyword
  PetscInt, intent(in)         :: target_type
  PetscInt                     :: vmtype

  character(len=MAXWORDLENGTH) :: internal_units, word, units
  PetscReal                    :: v

  ! Read a well surface target rate int v
  v = 0.0
  call InputReadDouble(input, option, v)
  call InputErrorMsg(input, option, keyword, 'WELL_DATA')
  ! Read the units into word
  call InputReadWord(input, option, word, PETSC_TRUE)
  vmtype=GetTargetUnitType(target_type)

  internal_units = 'unitless'
  if (vmtype == TT_P) internal_units = 'Pa'
  if (vmtype == TT_V) internal_units = 'm^3/sec'
  if (vmtype == TT_M) internal_units = 'Kg/sec'

  if (InputError(input)) then
    option%io_buffer = 'Keyword ' // trim(keyword) // ' units not found'
    call PrintErrMsg(option)
  else
    ! All OK, convert units and store
    units = trim(word)
    v = v*UnitsConvertToInternal(units, internal_units, option)
    if (this%w_readtime>0.0) then
      call StoreEvent(this, EVENT_TARG, target_type, v)
    else
      this%w_targets(target_type) = v
    endif
  endif

end subroutine readWellTarget

! *************************************************************************** !

subroutine SetWellType(this, well_type)
  !
  ! Store a well type
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this
  PetscInt, intent(in)  :: well_type

  PetscReal :: rval

  if (this%w_readtime>0.0) then
   rval = 0.0
   call StoreEvent(this, EVENT_TYPE, well_type, rval)
  else
    call UpdateWellType(this, well_type)
  endif

end subroutine SetWellType

! *************************************************************************** !

subroutine MarkCmplForDeletionInList(this, icmpl)
  !
  ! Mark completion icmpl for deletion
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this
  PetscInt, intent(in) :: icmpl

  this%c_mfd(icmpl) = PETSC_TRUE

end subroutine MarkCmplForDeletionInList

! *************************************************************************** !

subroutine DeleteMarkedCompletionsInList(this)
  !
  ! Delete marked completions (they are on other procs)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  class(well_data_type) :: this
  PetscInt :: icmpl, jcmpl

  this%w_ncmplg = this%w_ncmpl
  this%w_mcmplg = max(1, this%w_ncmplg)

  allocate(this%cg_ghosted_id(this%w_mcmplg))

  jcmpl = 0
  do icmpl = 1, this%w_ncmpl

    this%cg_ghosted_id(icmpl) = this%c_ghosted_id(icmpl)

    if (.not.this%c_mfd(icmpl)) then
      jcmpl = jcmpl+1

      this%c_to_cg      (jcmpl) = icmpl

      this%c_ci         (jcmpl) = this%c_ci         (icmpl)
      this%c_cj         (jcmpl) = this%c_cj         (icmpl)
      this%c_ck         (jcmpl) = this%c_ck         (icmpl)

      this%c_local_id   (jcmpl) = this%c_local_id   (icmpl)
      this%c_ghosted_id (jcmpl) = this%c_ghosted_id (icmpl)
      this%c_onproc     (jcmpl) = this%c_onproc     (icmpl)

      this%c_dx         (jcmpl) = this%c_dx         (icmpl)
      this%c_dy         (jcmpl) = this%c_dy         (icmpl)
      this%c_dz         (jcmpl) = this%c_dz         (icmpl)

      this%c_radius     (jcmpl) = this%c_radius     (icmpl)
      this%c_skin_factor(jcmpl) = this%c_skin_factor(icmpl)
      this%c_theta_frac (jcmpl) = this%c_theta_frac (icmpl)
      this%c_drill_dir  (jcmpl) = this%c_drill_dir  (icmpl)
      this%c_z          (jcmpl) = this%c_z          (icmpl)

      this%c_ccf        (jcmpl) = this%c_ccf        (icmpl)

    endif
  enddo
  this%w_ncmpl = jcmpl
  this%c_mfd   = PETSC_FALSE

end subroutine DeleteMarkedCompletionsInList

! *************************************************************************** !

function FindWellInList(iwell, well_data, list)
  !
  ! Find the well_data item for a given well index
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscBool :: FindWellInList
  PetscInt, intent(in) :: iwell
  class(well_data_type), pointer :: well_data
  type (well_data_list_type), pointer :: list

  FindWellInList = PETSC_FALSE
  well_data => list%first
  do
    if (.not.associated(well_data)) exit
    if (iwell == well_data%w_index) then
      FindWellInList = PETSC_TRUE
      exit
    endif
    well_data => well_data%next
  enddo

end function FindWellInList

! *************************************************************************** !

function FindWellInListZ(wname, well_data, list)
  !
  ! Find the well_data item for a given well index
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscBool :: FindWellInListZ
  character(len=*), intent(in) :: wname
  class(well_data_type), pointer :: well_data
  type (well_data_list_type), pointer :: list

  FindWellInListZ = PETSC_FALSE
  well_data => list%first
  do
    if (.not.associated(well_data)) exit
    if (wname == well_data%w_name) then
      FindWellInListZ = PETSC_TRUE
      exit
    endif
    well_data => well_data%next
  enddo

end function FindWellInListZ

! *************************************************************************** !

function GetWellStatusInList(this)
  !
  ! Get well status (open or shut)
  !
  ! Author: Dave Ponting
  ! Date: 10/25/18

  implicit none

  PetscInt :: GetWellStatusInList
  class(well_data_type) :: this

  GetWellStatusInList = this%w_status

end function GetWellStatusInList

! *************************************************************************** !

function getnwell(list)
  !
  ! Return the number of wells in the list of wells
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt :: getnwell
  type (well_data_list_type), pointer :: list

  getnwell = list%num_well

end function getnwell

! *************************************************************************** !

function GetTargetUnitType(itt)
  !
  ! Return the type (pressure, surface volume, mass) for a given target type
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  PetscInt :: GetTargetUnitType
  PetscInt, intent(in) :: itt

  GetTargetUnitType = -1
  if (itt == W_BHP_LIMIT) GetTargetUnitType = TT_P
  if (itt == W_TARG_OSV ) GetTargetUnitType = TT_V
  if (itt == W_TARG_GSV ) GetTargetUnitType = TT_V
  if (itt == W_TARG_WSV ) GetTargetUnitType = TT_V
  if (itt == W_TARG_SSV ) GetTargetUnitType = TT_V
  if (itt == W_TARG_LSV ) GetTargetUnitType = TT_V
  if (itt == W_TARG_OM  ) GetTargetUnitType = TT_M
  if (itt == W_TARG_GM  ) GetTargetUnitType = TT_M
  if (itt == W_TARG_WM  ) GetTargetUnitType = TT_M
  if (itt == W_TARG_SM  ) GetTargetUnitType = TT_M
  if (itt == W_TARG_RV  ) GetTargetUnitType = TT_M

  if (itt == -1) then
    call  throwWellDataException('Unable to identify target unit type')
  endif

end function getTargetUnitType

! *************************************************************************** !

subroutine StoreEvent(this, code, ival, rval)
  !
  ! Store an event
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  implicit none

  class(well_data_type)  :: this
  PetscInt,  intent(in)  :: code
  PetscInt,  intent(in)  :: ival
  PetscReal, intent(in)  :: rval

  call checkEventCount(this)

  this%w_nevent = this%w_nevent+1
  this%w_event_time(this%w_nevent) = this%w_readtime
  this%w_event_code(this%w_nevent) = code
  this%w_event_ival(this%w_nevent) = ival
  this%w_event_rval(this%w_nevent) = rval

end subroutine StoreEvent

! *************************************************************************** !

subroutine CheckEventCount(this)
  !
  ! Check number of elements in the event arrays
  ! Allocate if none, and extend if full
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  implicit none

  class(well_data_type) :: this

  PetscInt :: mold, mnew, nreq
  PetscReal :: vzero

  vzero = 0.0

  if (this%w_mevent == 0 .or. (this%w_nevent >= (this%w_mevent-1) ) ) then

    nreq = this%w_nevent+1
    mold = this%w_mevent
    mnew = mold

    call AllocOrReallocI(this%w_event_code, -1         , mold, mnew, nreq)
    call AllocOrReallocR(this%w_event_time, vzero      , mold, mnew, nreq)
    call AllocOrReallocI(this%w_event_ival, 0          , mold, mnew, nreq)
    call AllocOrReallocR(this%w_event_rval, vzero      , mold, mnew, nreq)
    call AllocOrReallocB(this%w_event_used, PETSC_FALSE, mold, mnew, nreq)

    this%w_mevent = mnew
  endif

end subroutine CheckEventCount

! *************************************************************************** !

subroutine CheckCompletionCount(this)
  !
  ! Check number of elements in the completion arrays
  ! Allocate if none, and extend if full
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  implicit none

  class(well_data_type) :: this

  PetscInt :: mold, mnew, nreq
  PetscReal :: vneg, vzero, vunity

  vneg   = -1.0
  vzero  =  0.0
  vunity =  1.0

  if (this%w_mCmpl == 0 .or. (this%w_ncmpl >= this%w_mcmpl)) then

    mold = this%w_mCmpl
    mnew = mold
    nreq = this%w_ncmpl

    call AllocOrReallocI(this%c_ci         , 0          , mold, mnew, nreq)
    call AllocOrReallocI(this%c_cj         , 0          , mold, mnew, nreq)
    call AllocOrReallocI(this%c_ck         , 0          , mold, mnew, nreq)

    call AllocOrReallocI(this%c_local_id   , 0          , mold, mnew, nreq)
    call AllocOrReallocI(this%c_ghosted_id , 0          , mold, mnew, nreq)
    call AllocOrReallocB(this%c_onproc     , PETSC_FALSE, mold, mnew, nreq)
    call AllocOrReallocI(this%c_to_cg      , 0          , mold, mnew, nreq)

    call AllocOrReallocR(this%c_dx         , vneg       , mold, mnew, nreq)
    call AllocOrReallocR(this%c_dy         , vneg       , mold, mnew, nreq)
    call AllocOrReallocR(this%c_dz         , vneg       , mold, mnew, nreq)

    call AllocOrReallocR(this%c_radius     , w_vrdef    , mold, mnew, nreq)
    call AllocOrReallocR(this%c_skin_factor, vzero      , mold, mnew, nreq)
    call AllocOrReallocR(this%c_theta_frac , vunity     , mold, mnew, nreq)
    call AllocOrReallocI(this%c_drill_dir  , Z_DIRECTION, mold, mnew, nreq)
    call AllocOrReallocR(this%c_z          , vzero      , mold, mnew, nreq)
    call AllocOrReallocR(this%c_ccf        , vzero      , mold, mnew, nreq)

    call AllocOrReallocB(this%c_mfd        , PETSC_FALSE, mold, mnew, nreq)

    this%w_mcmpl = mnew

  endif

end subroutine CheckCompletionCount

! ************************************************************************** !

subroutine AllocOrReallocI(pi, idef, mold, mnew, nreq)
  !
  ! Allocate or extend a integer array
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  use Utility_module, only : ReallocateArray

  Implicit none

  PetscInt, pointer, intent(inout) :: pi(:)
  PetscInt,          intent(in   ) :: idef
  PetscInt,          intent(in   ) :: mold
  PetscInt,          intent(out  ) :: mnew
  PetscInt,          intent(in   ) :: nreq

  mnew = mold

  if (mold == 0) then
    mnew = nreq
    allocate(pi(mnew))
    pi = idef
  else
    call ReallocateArray(pi, mnew)
    pi(mold+1:mnew) = idef
  endif

end subroutine AllocOrReallocI

! ************************************************************************** !

subroutine AllocOrReallocR(pr, rdef, mold, mnew, nreq)
  !
  ! Allocate or extend a real array
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  use Utility_module, only : ReallocateArray

  Implicit none

  PetscReal, pointer, intent(inout) :: pr(:)
  PetscReal,          intent(in   ) :: rdef
  PetscInt ,          intent(in   ) :: mold
  PetscInt ,          intent(out  ) :: mnew
  PetscInt ,          intent(in   ) :: nreq

  mnew = mold

  if (mold == 0) then
    mnew = nreq
    allocate(pr(mnew))
    pr = rdef
  else
    call ReallocateArray(pr, mnew)
    pr(mold+1:mnew) = rdef
  endif

end subroutine AllocOrReallocR

! ************************************************************************** !

subroutine AllocOrReallocB(pb, bdef, mold, mnew, nreq)
  !
  ! Allocate or extend a bool array
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  use Utility_module, only : ReallocateArray

  Implicit none

  PetscBool, pointer, intent(inout) :: pb(:)
  PetscBool,          intent(in   ) :: bdef
  PetscInt ,          intent(in   ) :: mold
  PetscInt ,          intent(out  ) :: mnew
  PetscInt ,          intent(in   ) :: nreq

  mnew = mold

  if (mold == 0) then
    mnew = nreq
    allocate(pb(mnew))
    pb = bdef
  else
    call ReallocateArray(pb, mnew)
    pb(mold+1:mnew) = bdef
  endif

end subroutine AllocOrReallocB

! *************************************************************************** !

subroutine UpdateWellType(this, itype)
  !
  ! Update the type of a well
  ! This needs care, as the well will generally need a fresh solve.
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  implicit none

  class(well_data_type) :: this
  PetscInt, intent(in)   :: itype
  PetscBool :: isbhp, wasbhp
  PetscReal :: old_bhp_limit

  ! Set up old and new target type values

  wasbhp = PETSC_FALSE
  isbhp  = PETSC_FALSE

  if (this%w_TT == W_BHP_LIMIT) wasbhp = PETSC_TRUE
  if (this%w_TT == W_BHP_LIMIT) isbhp  = PETSC_TRUE

  old_bhp_limit = this%w_targets(W_BHP_LIMIT)

  ! Make the change

  this%well_solution_set = PETSC_FALSE ! Unset solution to get clean predictor

  this%w_targets         = -1.0        ! Clear old targets

  ! if prod/inj status has not changed, leave the bhp limit alone

  if ((wasbhp .and. isbhp) .or. ((.not.wasbhp) .and. (.not.isbhp))) then
    this%w_targets(W_BHP_LIMIT) = old_bhp_limit
  endif

  ! Set to bhp control for first attempt on new type

  this%w_TT    = W_BHP_LIMIT

  ! Store new type

  this%w_itype = itype

end subroutine UpdateWellType

! *************************************************************************** !

subroutine getLocalAndGhostedIDs(natural_id, local_id, &
                                 ghosted_id, onproc, grid)
  !
  ! Find a completion from its natural id
  ! This needs a search
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

use Grid_module

  implicit none

  PetscInt,  intent(in)  :: natural_id
  PetscInt,  intent(out) :: local_id
  PetscInt,  intent(out) :: ghosted_id
  PetscBool, intent(out) :: onproc
  type(grid_type), pointer :: grid

  PetscInt :: il, ig, ia, ighost

  ! Initial values

  local_id   = 1
  ghosted_id = 1
  onproc     = PETSC_FALSE

  ! Search local cells

  do il = 1, grid%nlmax
    ig = grid%nL2G(il)
    ia = grid%nG2A(ig)
    if (ia == natural_id) then
      local_id   = il
      ghosted_id = ig
      onproc = PETSC_TRUE
      exit
    endif
  enddo

  ! Search other cells

  if (.not.onproc) then
    do ighost = 1, grid%ngmax
      ia = grid%nG2A(ighost)
      if (ia == natural_id) then
        local_id   = -1
        ghosted_id = ighost
        exit
       endif
    enddo
  endif

end subroutine getLocalAndGhostedIDs

! ************************************************************************** !

subroutine GetFieldData(fhpav)
  !
  ! Get field data not formed from well totals
  ! Currently just one item, the average hydrocarbon volume pressure 
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  implicit none

  PetscReal, intent(out) :: fhpav

  fhpav = f_hpav

end subroutine GetFieldData

! ************************************************************************** !

subroutine SetFieldData(fhpav)
  !
  ! Set field data not formed from well totals
  ! Currently just one item, the average hydrocarbon volume pressure 
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  implicit none

  PetscReal, intent(in) :: fhpav

  f_hpav = fhpav

end subroutine SetFieldData

! *************************************************************************** !

subroutine FindGroupRates(list)
  !
  ! Find well group rates
  ! Currently just one group, the whole field
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  implicit none

   type(well_data_list_type), pointer :: list

  PetscInt :: iw, nw, itt, itype
  PetscReal :: gactp, gacti, wact

  list%f_actualP = 0.0
  list%f_actualI = 0.0

  nw = list%num_well

  do itt = 1, N_WELL_TT

    gactp = 0.0
    gacti = 0.0

    do iw = 1, nw

      itype = getWellTypeI(iw, list)

      wact = GetWellTTValI(iw, itt, VALTYPE_ACTUAL, list)

      if (itype == PROD_WELL_TYPE) then
        gactp = gactp+wact
      else
        ! Note conversion to unsigned
        gacti = gacti-wact
      endif

    enddo

    if (itt /= W_BHP_LIMIT) then
      list%f_actualP(itt) = list%f_actualP(itt)+gactp
      list%f_actualI(itt) = list%f_actualI(itt)+gacti
    endif

  enddo

end subroutine FindGroupRates

! *************************************************************************** !

subroutine FindGroupTotals(list)
  !
  ! Find well group totals (cumulatives over time)
  ! Currently just one group, the whole field
  !
  ! Author: Dave Ponting
  ! Date: 10/23/18

  implicit none

  type(well_data_list_type), pointer :: list

  PetscInt :: iw, nw, itt
  PetscReal :: gtotp, gtoti, wtotp, wtoti

  list%f_totalP = 0.0
  list%f_totalI = 0.0

  nw = list%num_well

  do itt = 1, N_WELL_TT

    gtotp = 0.0
    gtoti = 0.0

    do iw = 1, nw

      wtotp = GetWellTTValI(iw, itt, VALTYPE_TOTALP, list)
      wtoti = GetWellTTValI(iw, itt, VALTYPE_TOTALI, list)

      gtotp = gtotp + wtotp
      gtoti = gtoti + wtoti

    enddo

    if (itt /= W_BHP_LIMIT) then
      list%f_totalP(itt) = list%f_totalP (itt)+gtotp
      list%f_totalI(itt) = list%f_totalI (itt)+gtoti
    endif

  enddo

end subroutine FindGroupTotals

! *************************************************************************** !

subroutine IncrementWellWarningCount( unconv_t, &
                                      unconv_w, &
                                      unconv_b, &
                                      mpierr )
  !
  ! Increment the well iteratin convergence failure counts
  !
  ! Author: Dave Ponting
  ! Date: 01/24/19

  implicit none

  PetscInt, intent(in) :: unconv_t, unconv_w, unconv_b,mpierr

  w_unconv_t = w_unconv_t + unconv_t
  w_unconv_w = w_unconv_w + unconv_w
  w_unconv_b = w_unconv_b + unconv_b
  w_mpierr   = w_mpierr   + mpierr

end subroutine IncrementWellWarningCount

! *************************************************************************** !

subroutine throwWellDataException(message)
  !
  ! Throw a general well_data error
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  character(len = *) :: message
  print *, message
  stop
end subroutine throwWellDataException

end module Well_Data_class
