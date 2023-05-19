
module PM_NWT_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_class
  use Option_module
  use PFLOTRAN_Constants_module
  use NW_Transport_module
  use NW_Transport_Aux_module

  implicit none

  private

  type, public :: pm_nwt_controls_type
    PetscReal, pointer :: itol_rel_update(:)
    PetscReal, pointer :: itol_abs_update(:)
    PetscReal, pointer :: itol_scaled_res(:)
    PetscReal, pointer :: itol_abs_res(:)
    PetscReal, pointer :: cnvg_criteria_value(:)
    PetscInt, pointer :: i_mapping(:)
    character(len=MAXWORDLENGTH), pointer :: names_itol_rel_update(:)
    character(len=MAXWORDLENGTH), pointer :: names_itol_abs_update(:)
    character(len=MAXWORDLENGTH), pointer :: names_itol_scaled_res(:)
    character(len=MAXWORDLENGTH), pointer :: names_itol_abs_res(:)
    PetscInt :: max_newton_iterations
    PetscReal :: dt_cut
    PetscBool :: check_post_converged
    PetscBool :: check_post_convergence
    PetscBool :: check_update
    PetscBool :: verbose_newton
    PetscBool :: scaling_cut_dt
#ifdef OS_STATISTICS
! use PetscReal for large counts
    PetscInt :: newton_call_count
    PetscReal :: sum_newton_call_count
    PetscInt :: newton_iterations
    PetscReal :: sum_newton_iterations
    PetscInt :: overall_max_newton_iterations
#endif
  end type pm_nwt_controls_type

  type, public :: pm_nwt_params_type
    PetscInt :: nphase
    PetscInt :: ncomp
    PetscInt :: nspecies
    PetscInt :: nauxiliary
    PetscReal :: tran_weight_t0
    PetscReal :: tran_weight_t1
    PetscBool :: calculate_transverse_dispersion
    PetscBool :: temperature_dependent_diffusion
    PetscBool :: transient_porosity
    PetscBool :: steady_flow
    PetscBool :: zero_out_borehole
    PetscReal :: bh_zero_value
    PetscReal :: init_total_mass_conc
    PetscReal :: wm_start_time
    PetscReal :: wm_end_time
    PetscReal :: wm_value
    character(len=MAXWORDLENGTH), pointer :: bh_material_names(:)
    character(len=MAXWORDLENGTH), pointer :: dirichlet_material_names(:)
    PetscInt, pointer :: bh_material_ids(:)
    PetscInt, pointer :: dirichlet_material_ids(:)

  end type pm_nwt_params_type

  type, public, extends(pm_base_type) :: pm_nwt_type
  ! realization_base_type has the nwt object (equivalent to reaction)
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    type(pm_nwt_controls_type), pointer :: controls
    type(pm_nwt_params_type), pointer :: params
  contains
    procedure, public :: Setup => PMNWTSetup
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMNWTReadSimOptionsBlock
    procedure, public :: ReadTSBlock => PMNWTReadTSSelectCase
    procedure, public :: ReadNewtonBlock => PMNWTReadNewtonSelectCase
    procedure, public :: SetRealization => PMNWTSetRealization
    procedure, public :: InitializeRun => PMNWTInitializeRun
    procedure, public :: FinalizeRun => PMNWTFinalizeRun
    procedure, public :: InitializeTimestep => PMNWTInitializeTimestep
    procedure, public :: FinalizeTimestep => PMNWTFinalizeTimestep
    procedure, public :: UpdateTimestep => PMNWTUpdateTimestep
    procedure, public :: Residual => PMNWTResidual
    procedure, public :: Jacobian => PMNWTJacobian
    procedure, public :: PreSolve => PMNWTPreSolve
    procedure, public :: PostSolve => PMNWTPostSolve
    procedure, public :: UpdateSolution => PMNWTUpdateSolution
    procedure, public :: CheckConvergence => PMNWTCheckConvergence
    procedure, public :: CheckUpdatePre => PMNWTCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMNWTCheckUpdatePost
    procedure, public :: AcceptSolution => PMNWTAcceptSolution
    procedure, public :: TimeCut => PMNWTTimeCut
    procedure, public :: SetTranWeights => PMNWTSetTranWeights
    procedure, public :: UpdateAuxVars => PMNWTUpdateAuxVars
    procedure, public :: MaxChange => PMNWTMaxChange
    procedure, public :: CheckpointBinary => PMNWTCheckpointBinary
    procedure, public :: CheckpointHDF5 => PMNWTCheckpointHDF5
    procedure, public :: RestartBinary => PMNWTRestartBinary
    procedure, public :: RestartHDF5 => PMNWTRestartHDF5
    procedure, public :: InputRecord => PMNWTInputRecord
    procedure, public :: Destroy => PMNWTDestroy
  end type pm_nwt_type

  public :: PMNWTCreate


contains

! ************************************************************************** !

function PMNWTCreate()
  !
  ! Creates the nuclear waste transport process model shell.
  !
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  !

  implicit none

  class(pm_nwt_type), pointer :: PMNWTCreate

  class(pm_nwt_type), pointer :: nwt_pm

  allocate(nwt_pm)
  nullify(nwt_pm%option)
  nullify(nwt_pm%output_option)
  nullify(nwt_pm%realization)
  nullify(nwt_pm%comm1)

  allocate(nwt_pm%controls)
  nullify(nwt_pm%controls%itol_rel_update)
  nullify(nwt_pm%controls%itol_abs_update)
  nullify(nwt_pm%controls%itol_scaled_res)
  nullify(nwt_pm%controls%itol_abs_res)
  nullify(nwt_pm%controls%cnvg_criteria_value)
  nullify(nwt_pm%controls%i_mapping)
  nullify(nwt_pm%controls%names_itol_rel_update)
  nullify(nwt_pm%controls%names_itol_abs_update)
  nullify(nwt_pm%controls%names_itol_scaled_res)
  nullify(nwt_pm%controls%names_itol_abs_res)
  nwt_pm%controls%max_newton_iterations = 10
  nwt_pm%controls%dt_cut = 0.5d0
  nwt_pm%controls%scaling_cut_dt = PETSC_FALSE
  nwt_pm%controls%check_post_converged = PETSC_FALSE
  nwt_pm%controls%check_post_convergence = PETSC_FALSE
  nwt_pm%controls%check_update = PETSC_TRUE
  nwt_pm%controls%verbose_newton = PETSC_FALSE
#ifdef OS_STATISTICS
  nwt_pm%controls%newton_call_count = 0
  nwt_pm%controls%sum_newton_call_count = 0.d0
  nwt_pm%controls%newton_iterations = 0
  nwt_pm%controls%sum_newton_iterations = 0.d0
  nwt_pm%controls%overall_max_newton_iterations = 0
#endif

  allocate(nwt_pm%params)
  nwt_pm%params%ncomp = 0
  nwt_pm%params%nphase = 0
  nwt_pm%params%nspecies = 0
  nwt_pm%params%nauxiliary = 0
  nwt_pm%params%tran_weight_t0 = 0.d0
  nwt_pm%params%tran_weight_t1 = 0.d0
  nwt_pm%params%transient_porosity = PETSC_FALSE
  nwt_pm%params%calculate_transverse_dispersion = PETSC_FALSE
  nwt_pm%params%temperature_dependent_diffusion = PETSC_FALSE
  nwt_pm%params%steady_flow = PETSC_FALSE
  nwt_pm%params%zero_out_borehole = PETSC_FALSE
  nwt_pm%params%bh_zero_value = 1.0d-40  ! [mol/m3-bulk]
  nwt_pm%params%wm_start_time = UNINITIALIZED_DOUBLE  ! [sec]
  nwt_pm%params%wm_end_time = UNINITIALIZED_DOUBLE  ! [sec]
  nwt_pm%params%wm_value = 1.0d-40  ! [mol/m3-bulk]
  nwt_pm%params%init_total_mass_conc = UNINITIALIZED_DOUBLE
  nullify(nwt_pm%params%bh_material_names)
  nullify(nwt_pm%params%dirichlet_material_names)
  nullify(nwt_pm%params%bh_material_ids)
  nullify(nwt_pm%params%dirichlet_material_ids)


  call PMBaseInit(nwt_pm)
  nwt_pm%name = 'Nuclear Waste Transport'
  nwt_pm%header = 'NUCLEAR WASTE TRANSPORT'

  PMNWTCreate => nwt_pm

end function PMNWTCreate

! ************************************************************************** !

subroutine PMNWTReadSimOptionsBlock(this,input)
  !
  ! Reads input file parameters associated with the nuclear waste transport
  ! process model in the SIMULATION block.
  !
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  !
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_nwt_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option

  error_string = 'NUCLEAR_WASTE_TRANSPORT OPTIONS'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMBaseReadSimOptionsSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))
!geh: yet to be implemented
!      case('TEMPERATURE_DEPENDENT_DIFFUSION')
!        this%temperature_dependent_diffusion = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine PMNWTReadSimOptionsBlock

! ************************************************************************** !

subroutine PMNWTReadTSSelectCase(this,input,keyword,found, &
                                 error_string,option)
  !
  ! Read timestepper settings specific to this process model
  !
  ! Author: Glenn Hammond
  ! Date: 03/23/20

  use Input_Aux_module
  use Option_module

  implicit none

  class(pm_nwt_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

!  found = PETSC_TRUE
!  call PMBaseReadSelectCase(this,input,keyword,found,error_string,option)
!  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case default
      found = PETSC_FALSE
  end select

end subroutine PMNWTReadTSSelectCase

! ************************************************************************** !

subroutine PMNWTReadNewtonSelectCase(this,input,keyword,found, &
                                     error_string,option)
  !
  ! Reads input file parameters associated with the NWT process model
  ! Newton solver convergence
  !
  ! Author: Glenn Hammond
  ! Date: 03/25/20

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_nwt_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  character(len=MAXWORDLENGTH), pointer :: temp_species_names(:)
  PetscReal, pointer :: temp_itol_abs_update(:)
  PetscReal, pointer :: temp_itol_rel_update(:)
  PetscReal, pointer :: temp_itol_scaled_res(:)
  PetscReal, pointer :: temp_itol_abs_res(:)
  PetscInt :: k
  character(len=MAXSTRINGLENGTH) :: error_string_ex
  character(len=MAXWORDLENGTH) :: word

  option => this%option

  allocate(temp_species_names(50))

  found = PETSC_TRUE
  select case(trim(keyword))
    !------------------------------------------------------------------------
    case('MAXIMUM_NUMBER_OF_ITERATIONS')
      error_string_ex = trim(error_string) // ',MAXIMUM_NUMBER_OF_ITERATIONS'
      call InputReadInt(input,option,this%controls%max_newton_iterations)
      call InputErrorMsg(input,option,'VALUE',error_string)
      this%solver%newton_max_iterations = this%controls%max_newton_iterations
    !------------------------------------------------------------------------
    case('VERBOSE_NEWTON')
      error_string_ex = trim(error_string) // ',VERBOSE_NEWTON'
      this%controls%verbose_newton = PETSC_TRUE
    !------------------------------------------------------------------------
    case('NWT_ITOL_RELATIVE_UPDATE')
      k = 0
      temp_species_names = ''
      error_string_ex = trim(error_string) // ',NWT_ITOL_RELATIVE_UPDATE'
      allocate(temp_itol_rel_update(50))
      temp_itol_rel_update = -999.99999
      do
        k = k + 1
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'SPECIES NAME',error_string_ex)
        call StringToUpper(word)
        temp_species_names(k) = trim(word)
        call InputReadDouble(input,option,temp_itol_rel_update(k))
        call InputErrorMsg(input,option,'SPECIES TOLERANCE VALUE', &
                           error_string_ex)
      enddo
      allocate(this%controls%itol_rel_update(k-1))
      this%controls%itol_rel_update(:) = temp_itol_rel_update(1:k-1)
      allocate(this%controls%names_itol_rel_update(k-1))
      this%controls%names_itol_rel_update(:) = temp_species_names(1:k-1)
      deallocate(temp_itol_rel_update)
    !------------------------------------------------------------------------
    case('NWT_ITOL_ABSOLUTE_UPDATE')
      k = 0
      temp_species_names = ''
      error_string_ex = trim(error_string) // ',NWT_ITOL_ABSOLUTE_UPDATE'
      allocate(temp_itol_abs_update(50))
      temp_itol_abs_update = -999.99999
      do
        k = k + 1
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'SPECIES NAME',error_string_ex)
        call StringToUpper(word)
        temp_species_names(k) = trim(word)
        call InputReadDouble(input,option,temp_itol_abs_update(k))
        call InputErrorMsg(input,option,'SPECIES TOLERANCE VALUE', &
                           error_string_ex)
      enddo
      allocate(this%controls%itol_abs_update(k-1))
      this%controls%itol_abs_update(:) = temp_itol_abs_update(1:k-1)
      allocate(this%controls%names_itol_abs_update(k-1))
      this%controls%names_itol_abs_update(:) = temp_species_names(1:k-1)
      deallocate(temp_itol_abs_update)
    !------------------------------------------------------------------------
    case('NWT_ITOL_SCALED_RESIDUAL')
      k = 0
      temp_species_names = ''
      error_string_ex = trim(error_string) // ',NWT_ITOL_SCALED_RESIDUAL'
      allocate(temp_itol_scaled_res(50))
      temp_itol_scaled_res = -999.99999
      do
        k = k + 1
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'SPECIES NAME',error_string_ex)
        call StringToUpper(word)
        temp_species_names(k) = trim(word)
        call InputReadDouble(input,option,temp_itol_scaled_res(k))
        call InputErrorMsg(input,option,'SPECIES TOLERANCE VALUE', &
                           error_string_ex)
      enddo
      allocate(this%controls%itol_scaled_res(k-1))
      this%controls%itol_scaled_res(:) = temp_itol_scaled_res(1:k-1)
      allocate(this%controls%names_itol_scaled_res(k-1))
      this%controls%names_itol_scaled_res(:) = temp_species_names(1:k-1)
      deallocate(temp_itol_scaled_res)
    !------------------------------------------------------------------------
    case('NWT_ITOL_ABSOLUTE_RESIDUAL')
      k = 0
      temp_species_names = ''
      error_string_ex = trim(error_string) // ',ITOL_ABSOLUTE_RESIDUAL'
      allocate(temp_itol_abs_res(50))
      temp_itol_abs_res = -999.99999
      do
        k = k + 1
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'SPECIES NAME',error_string_ex)
        call StringToUpper(word)
        temp_species_names(k) = trim(word)
        call InputReadDouble(input,option,temp_itol_abs_res(k))
        call InputErrorMsg(input,option,'SPECIES TOLERANCE VALUE', &
                           error_string_ex)
      enddo
      allocate(this%controls%itol_abs_res(k-1))
      this%controls%itol_abs_res(:) = temp_itol_abs_res(1:k-1)
      allocate(this%controls%names_itol_abs_res(k-1))
      this%controls%names_itol_abs_res(:) = temp_species_names(1:k-1)
      deallocate(temp_itol_abs_res)
    !------------------------------------------------------------------------
  case default
      found = PETSC_FALSE
  end select

  deallocate(temp_species_names)

end subroutine PMNWTReadNewtonSelectCase

! ************************************************************************** !

subroutine PMNWTSetup(this)
  !
  ! Initializes variables associated with nuclear waste transport.
  !
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  !

  use String_module

  implicit none

  class(pm_nwt_type) :: this

  class(reaction_nw_type), pointer :: reaction_nw
  PetscInt :: k, j
  character(len=MAXWORDLENGTH) :: name_species
  character(len=MAXWORDLENGTH) :: name_cnvgcrit

  reaction_nw => this%realization%reaction_nw

  if (this%option%nflowdof > 0) then
    this%option%flow%store_state_variables_in_global = PETSC_TRUE
  endif
  this%params%nphase = reaction_nw%params%nphase
  this%params%ncomp = reaction_nw%params%ncomp
  this%params%nspecies = reaction_nw%params%nspecies
  this%params%nauxiliary = reaction_nw%params%nauxiliary
  this%params%calculate_transverse_dispersion = &
                           reaction_nw%params%calculate_transverse_dispersion
  this%params%temperature_dependent_diffusion = &
                           reaction_nw%params%temperature_dependent_diffusion
  if(associated(reaction_nw%params%bh_material_names)) then
    this%params%bh_zero_value = reaction_nw%params%bh_zero_value
    allocate(this%params%bh_material_names(&
                       size(reaction_nw%params%bh_material_names)))
    this%params%bh_material_names = reaction_nw%params%bh_material_names
    allocate(this%params%bh_material_ids(&
                       size(reaction_nw%params%bh_material_names)))
    this%params%bh_material_ids = UNINITIALIZED_INTEGER
  endif
  if(reaction_nw%screening_run) then
    allocate(this%params%dirichlet_material_names(&
                       size(reaction_nw%params%dirichlet_material_names)))
    this%params%dirichlet_material_names = &
                               reaction_nw%params%dirichlet_material_names
    allocate(this%params%dirichlet_material_ids(&
                       size(reaction_nw%params%dirichlet_material_names)))
    this%params%dirichlet_material_ids = UNINITIALIZED_INTEGER
    this%params%init_total_mass_conc = reaction_nw%params%init_total_mass_conc
  endif
  this%params%wm_start_time = reaction_nw%params%wm_start_time
  this%params%wm_end_time = reaction_nw%params%wm_end_time
  this%params%wm_value = reaction_nw%params%wm_value

  ! set the communicator
  this%comm1 => this%realization%comm1

  allocate(this%controls%cnvg_criteria_value(this%params%nspecies))
  allocate(this%controls%i_mapping(this%params%nspecies))

  if (.not.associated(this%controls%itol_abs_res)) then
    this%option%io_buffer = 'A NWT_ITOL_ABSOLUTE_RESIDUAL block is required &
      &in the NUMERICAL_METHODS,NEWTON_SOLVER block for SUBSURFACE_TRANSPORT &
      &MODE NWT.'
    call PrintErrMsg(this%option)
  endif
  if (.not.associated(this%controls%itol_scaled_res)) then
    this%option%io_buffer = 'A NWT_ITOL_SCALED_RESIDUAL block is required &
      &in the NUMERICAL_METHODS,NEWTON_SOLVER block for SUBSURFACE_TRANSPORT &
      &MODE NWT.'
    call PrintErrMsg(this%option)
  endif
  if (.not.associated(this%controls%itol_rel_update)) then
    this%option%io_buffer = 'A NWT_ITOL_RELATIVE_UPDATE block is required &
      &in the NUMERICAL_METHODS,NEWTON_SOLVER block for SUBSURFACE_TRANSPORT &
      &MODE NWT.'
    call PrintErrMsg(this%option)
  endif
  if (.not.associated(this%controls%itol_abs_update)) then
    this%option%io_buffer = 'A NWT_ITOL_ABSOLUTE_UPDATE block is required &
      &in the NUMERICAL_METHODS,NEWTON_SOLVER block for SUBSURFACE_TRANSPORT &
      &MODE NWT.'
    call PrintErrMsg(this%option)
  endif

  if (associated(this%controls%itol_abs_res)) then
    if (size(this%controls%itol_abs_res) /= this%params%nspecies) then
      this%option%io_buffer = 'The number of species-dependent values of &
        &ITOL_ABSOLUTE_RESIDUAL provided in the NUMERICAL_METHODS,NEWTON_SOLVER &
        &block for SUBSURFACE_TRANSPORT MODE NWT does not match the number of &
        &species provided in the NUCLEAR_WASTE_CHEMISTRY block.'
      call PrintErrMsg(this%option)
    endif
  endif
  if (associated(this%controls%itol_scaled_res)) then
    if (size(this%controls%itol_scaled_res) /= this%params%nspecies) then
      this%option%io_buffer = 'The number of species-dependent values of &
        &ITOL_SCALED_RESIDUAL provided in the NUMERICAL_METHODS,NEWTON_SOLVER &
        &block for SUBSURFACE_TRANSPORT MODE NWT does not match the number of &
        &species provided in the NUCLEAR_WASTE_CHEMISTRY block.'
      call PrintErrMsg(this%option)
    endif
  endif
  if (associated(this%controls%itol_rel_update)) then
    if (size(this%controls%itol_rel_update) /= this%params%nspecies) then
      this%option%io_buffer = 'The number of species-dependent values of &
        &ITOL_RELATIVE_UPDATE provided in the NUMERICAL_METHODS,NEWTON_SOLVER &
        &block for SUBSURFACE_TRANSPORT MODE NWT does not match the number of &
        &species provided in the NUCLEAR_WASTE_CHEMISTRY block.'
      call PrintErrMsg(this%option)
    endif
  endif
  if (associated(this%controls%itol_abs_update)) then
    if (size(this%controls%itol_abs_update) /= this%params%nspecies) then
      this%option%io_buffer = 'The number of species-dependent values of &
        &ITOL_ABSOLUTE_UPDATE provided in the NUMERICAL_METHODS,NEWTON_SOLVER &
        &block for SUBSURFACE_TRANSPORT MODE NWT does not match the number of &
        &species provided in the NUCLEAR_WASTE_CHEMISTRY block.'
      call PrintErrMsg(this%option)
    endif
  endif

  ! reorder the convergence based on the reaction_nw%species_names(:) order
  !---------------------------------------------------------------------------
  do k = 1,this%params%nspecies
    this%controls%i_mapping(k) = 0
    name_species = trim(reaction_nw%species_names(k))
    call StringToUpper(name_species)
    do j = 1,this%params%nspecies
      name_cnvgcrit = trim(this%controls%names_itol_abs_res(j))
      call StringToUpper(name_cnvgcrit)
      if (name_cnvgcrit == name_species) then
        this%controls%i_mapping(k) = j ! gives the index in controls order
      endif
    enddo
    if (this%controls%i_mapping(k) == 0) then
      this%option%io_buffer = 'Species ' // trim(name_species) // ' not &
        &found among species included in the NUMERICAL_METHODS,NEWTON_SOLVER,&
        &NWT_ITOL_ABSOLUTE_RESIDUAL block.'
      call PrintErrMsg(this%option)
    endif
    this%controls%cnvg_criteria_value(k) = &
                    this%controls%itol_abs_res(this%controls%i_mapping(k))
  enddo
  this%controls%itol_abs_res = this%controls%cnvg_criteria_value
  !---------------------------------------------------------------------------
  do k = 1,this%params%nspecies
    this%controls%i_mapping(k) = 0
    name_species = trim(reaction_nw%species_names(k))
    call StringToUpper(name_species)
    do j = 1,this%params%nspecies
      name_cnvgcrit = trim(this%controls%names_itol_scaled_res(j))
      call StringToUpper(name_cnvgcrit)
      if (name_cnvgcrit == name_species) then
        this%controls%i_mapping(k) = j ! gives the index in controls order
      endif
    enddo
    if (this%controls%i_mapping(k) == 0) then
      this%option%io_buffer = 'Species ' // trim(name_species) // ' not &
        &found among species included in the NUMERICAL_METHODS,NEWTON_SOLVER,&
        &NWT_ITOL_SCALED_RESIDUAL block.'
      call PrintErrMsg(this%option)
    endif
    this%controls%cnvg_criteria_value(k) = &
                    this%controls%itol_scaled_res(this%controls%i_mapping(k))
  enddo
  this%controls%itol_scaled_res = this%controls%cnvg_criteria_value
  !---------------------------------------------------------------------------
  do k = 1,this%params%nspecies
    this%controls%i_mapping(k) = 0
    name_species = trim(reaction_nw%species_names(k))
    call StringToUpper(name_species)
    do j = 1,this%params%nspecies
      name_cnvgcrit = trim(this%controls%names_itol_rel_update(j))
      call StringToUpper(name_cnvgcrit)
      if (name_cnvgcrit == name_species) then
        this%controls%i_mapping(k) = j ! gives the index in controls order
      endif
    enddo
    if (this%controls%i_mapping(k) == 0) then
      this%option%io_buffer = 'Species ' // trim(name_species) // ' not &
        &found among species included in the NUMERICAL_METHODS,NEWTON_SOLVER,&
        &NWT_ITOL_RELATIVE_UPDATE block.'
      call PrintErrMsg(this%option)
    endif
    this%controls%cnvg_criteria_value(k) = &
                    this%controls%itol_rel_update(this%controls%i_mapping(k))
  enddo
  this%controls%itol_rel_update = this%controls%cnvg_criteria_value
  !---------------------------------------------------------------------------
  do k = 1,this%params%nspecies
    this%controls%i_mapping(k) = 0
    name_species = trim(reaction_nw%species_names(k))
    call StringToUpper(name_species)
    do j = 1,this%params%nspecies
      name_cnvgcrit = trim(this%controls%names_itol_abs_update(j))
      call StringToUpper(name_cnvgcrit)
      if (name_cnvgcrit == name_species) then
        this%controls%i_mapping(k) = j ! gives the index in controls order
      endif
    enddo
    if (this%controls%i_mapping(k) == 0) then
      this%option%io_buffer = 'Species ' // trim(name_species) // ' not &
        &found among species included in the NUMERICAL_METHODS,NEWTON_SOLVER,&
        &NWT_ITOL_ABSOLUTE_UPDATE block.'
      call PrintErrMsg(this%option)
    endif
    this%controls%cnvg_criteria_value(k) = &
                    this%controls%itol_abs_update(this%controls%i_mapping(k))
  enddo
  this%controls%itol_abs_update = this%controls%cnvg_criteria_value

end subroutine PMNWTSetup

! ************************************************************************** !

subroutine PMNWTSetRealization(this,realization)
  !
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  !

  implicit none

  class(pm_nwt_type) :: this
  class(realization_subsurface_type), pointer :: realization

  this%realization => realization
  this%realization_base => realization

  if (this%realization%reaction_nw%use_log_formulation) then
    this%solution_vec = realization%field%tran_log_xx
  else
    this%solution_vec = realization%field%tran_xx
  endif
  this%residual_vec = realization%field%tran_r

end subroutine PMNWTSetRealization

! ************************************************************************** !

subroutine PMNWTInitializeRun(this)
  !
  ! Initializes process model time stepping.
  !
  ! Author: Jenn Frederick
  ! Date: 04/02/2019
  !

  use Material_module
  use Coupler_module
  use Region_module

  implicit none

  class(pm_nwt_type) :: this

  PetscInt :: p
  type(material_property_type), pointer :: material_property
  type(region_type), pointer :: region
  character(len=MAXSTRINGLENGTH) :: hack_region_name

  ! check for uninitialized flow variables
  call RealizUnInitializedVarsTran(this%realization)

  ! update the boundary conditions
  !geh: need to update cells also, as the flow solution may have changed
  !     during restart and transport may have been skipped
  call NWTUpdateAuxVars(this%realization,PETSC_TRUE,PETSC_TRUE)

  this%realization%patch%aux%NWT%truncate_output = &
    this%realization%reaction_nw%truncate_output

  ! set borehole material ids
  if (associated(this%params%bh_material_names)) then
    this%params%zero_out_borehole = PETSC_TRUE
    do p = 1, size(this%params%bh_material_names)
      material_property => &
        MaterialPropGetPtrFromList(this%params%bh_material_names(p), &
                                   this%realization%patch%material_properties)
      if (.not.associated(material_property)) then
        this%option%io_buffer = 'Borehole material "' // &
          trim(this%params%bh_material_names(p)) // &
          '" not found among material properties.'
        call PrintErrMsg(this%option)
      endif
      this%params%bh_material_ids(p) = material_property%internal_id
    enddo
  endif

  ! set screening run dirichlet material ids
  if (this%realization%reaction_nw%screening_run) then
    do p = 1, size(this%params%dirichlet_material_names)
      material_property => &
        MaterialPropGetPtrFromList(this%params%dirichlet_material_names(p), &
                                   this%realization%patch%material_properties)
      if (.not.associated(material_property)) then
        this%option%io_buffer = 'Dirichlet material "' // &
          trim(this%params%dirichlet_material_names(p)) // &
          '" not found among material properties.'
        call PrintErrMsg(this%option)
      endif
      this%params%dirichlet_material_ids(p) = material_property%internal_id
    enddo
  endif

  !===========================================================================!
  ! Jenn's hack for upper borehole
  if (associated(this%params%bh_material_names)) then
    hack_region_name = 'rBH_OPEN_UPPER'
    region => RegionGetPtrFromList(hack_region_name, &
                                   this%realization%patch%region_list)
    !if (.not.associated(region)) then
    !  this%option%io_buffer = 'Borehole region "rBH_OPEN_UPPER" not found &
    !                          &among regions. This is a hack!'
    !  call PrintErrMsg(this%option)
    !endif
  endif
  !===========================================================================!

  call PMNWTUpdateSolution(this)

end subroutine PMNWTInitializeRun

! ************************************************************************** !

subroutine PMNWTInitializeTimestep(this)
  !
  ! Author: Jenn Frederick
  ! Date: 04/03/2019
  !

  use Global_module
  use Material_module
  use Patch_module
  use Field_module

  implicit none

  class(pm_nwt_type) :: this

  PetscErrorCode :: ierr

  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  type(material_property_type), pointer :: material_property
  PetscInt :: local_id, ghosted_id
  PetscInt :: i, idof
  PetscInt :: offset, index
  PetscInt :: istart, iend
  PetscReal, pointer :: xx_p(:)
  character(len=MAXSTRINGLENGTH) :: hack_material_name

  this%option%tran_dt = this%option%dt
  field => this%realization%field
  patch => this%realization%patch
  nwt_auxvars => this%realization%patch%aux%NWT%auxvars

  !call PMBasePrintHeader(this)

  ! If a material change to borehole materials has happened, remove all
  ! NWT species mass from the borehole material region.
  if (this%params%zero_out_borehole) then
    do local_id = 1, patch%grid%nlmax
      ghosted_id = patch%grid%nL2G(local_id)
      ! loop through the borehole material IDs:
      do i = 1, size(this%params%bh_material_ids)
        if (patch%imat(ghosted_id) == this%params%bh_material_ids(i)) then
          ! Means the current grid cell is one of the borehole materials
          ! We now want to zero out all the mass from NWT process model. But,
          ! it should only happen at the very first timestep that we find
          ! this borehole material. So after we get here, we should set the
          ! flag zero_out_borehole to FALSE.
          ! NOTE: Will this work if there are two intrusions? No!

          call VecGetArrayReadF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
          ! compute offset in solution vector for first dof in grid cell
          offset = (local_id-1)*this%params%nspecies
          do idof = 1, this%params%nspecies
            index = idof + offset
            xx_p(index) = this%params%bh_zero_value
          enddo
          ! calculate range of species
          istart = offset + 1
          iend = offset + this%params%nspecies
          nwt_auxvars(ghosted_id)%total_bulk_conc = xx_p(istart:iend)
          call VecRestoreArrayReadF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)

          write(this%option%io_buffer,*) local_id
          this%option%io_buffer = "All NWT mass zero'd at cell " // &
            trim(adjustl(this%option%io_buffer)) // &
            ' due to borehole material.'
          call PrintMsg(this%option)
          write(this%option%io_buffer,*) this%params%bh_zero_value
          this%option%io_buffer = 'The "zero" value is set to ' // &
            trim(adjustl(this%option%io_buffer)) // ' mol/m3.'
          call PrintMsg(this%option)

          this%params%zero_out_borehole = PETSC_FALSE
        endif
      enddo
    enddo
  endif

  ! If this is a screening run, then hold the mass in the list of Dirichlet
  ! materials constant (at initial value).
  if (this%realization%reaction_nw%screening_run) then
    do local_id = 1, patch%grid%nlmax
      ghosted_id = patch%grid%nL2G(local_id)
      ! loop through the dirichlet material IDs:
      do i = 1, size(this%params%dirichlet_material_ids)
        if (patch%imat(ghosted_id) == this%params%dirichlet_material_ids(i)) then
          ! Means the current grid cell is one of the dirichlet materials
          call VecGetArrayReadF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
          ! compute offset in solution vector for first dof in grid cell
          offset = (local_id-1)*this%params%nspecies
          do idof = 1, this%params%nspecies
            index = idof + offset
            xx_p(index) = this%params%init_total_mass_conc
          enddo
          ! calculate range of species
          istart = offset + 1
          iend = offset + this%params%nspecies
          nwt_auxvars(ghosted_id)%total_bulk_conc = xx_p(istart:iend)
          call VecRestoreArrayReadF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
        endif
      enddo
    enddo
  endif

  !===========================================================================!
  ! Jenn's hack for upper borehole:
  ! You must include a material property called BH_OPEN_UPPER, and if you
  ! do, then the concentration there is kept at background levels.
  hack_material_name = 'BH_OPEN_UPPER'
  material_property => &
       MaterialPropGetPtrFromList(hack_material_name, &
                                  this%realization%patch%material_properties)
  if (associated(material_property) .and. &
      Uninitialized(this%params%wm_start_time)) then
    this%option%io_buffer = 'START_TIME was not provided in the &
      &WASHING_MACHINE block for SUBSURFACE_TRANSPORT MODE NWT.'
    call PrintErrMsg(this%option)
  endif
  if (associated(material_property) .and. &
      Uninitialized(this%params%wm_end_time)) then
    this%option%io_buffer = 'END_TIME was not provided in the &
      &WASHING_MACHINE block for SUBSURFACE_TRANSPORT MODE NWT.'
    call PrintErrMsg(this%option)
  endif
  if (.not.associated(material_property) .and. &
      ( Initialized(this%params%wm_end_time) .or. &
        Initialized(this%params%wm_start_time) ) ) then
    this%option%io_buffer = 'WASHING_MACHINE block for SUBSURFACE_TRANSPORT &
      &MODE NWT was provided, but region and material BH_OPEN_UPPER not found.'
    call PrintErrMsg(this%option)
  endif
  if (associated(material_property)) then
    do local_id = 1, patch%grid%nlmax
      ghosted_id = patch%grid%nL2G(local_id)
      ! if within BH_OPEN time interval:
      if (this%option%time >= this%params%wm_start_time .and. &
          this%option%time < this%params%wm_end_time) then
        if (patch%imat(ghosted_id) == material_property%internal_id) then
          ! Means the current grid cell is in BH_OPEN_UPPER
          call VecGetArrayReadF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
          ! compute offset in solution vector for first dof in grid cell
          offset = (local_id-1)*this%params%nspecies
          do idof = 1, this%params%nspecies
            index = idof + offset
            xx_p(index) = this%params%wm_value  ! mol/m3-bulk
          enddo
          ! calculate range of species
          istart = offset + 1
          iend = offset + this%params%nspecies
          nwt_auxvars(ghosted_id)%total_bulk_conc = xx_p(istart:iend)
          call VecRestoreArrayReadF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
        endif
      endif
    enddo
  endif
  !===========================================================================!

  ! interpolate flow parameters/data
  ! this must remain here as these weighted values are used by both
  ! NWTInitializeTimestep and NWTTimeCut (which calls NWTInitializeTimestep)
  if ((this%option%nflowdof > 0) .and. (.not.this%params%steady_flow)) then
    call this%SetTranWeights()
    if (this%option%flow%transient_porosity) then
      ! weight material properties (e.g. porosity)
      call MaterialWeightAuxVars(this%realization%patch%aux%Material, &
                                 this%params%tran_weight_t0, &
                                 this%realization%field,this%comm1)
    endif
    ! set densities and saturations to t
    call GlobalWeightAuxVars(this%realization,this%params%tran_weight_t0)
  else if (this%params%transient_porosity) then
    this%params%tran_weight_t0 = 0.d0
    call MaterialWeightAuxVars(this%realization%patch%aux%Material, &
                               this%params%tran_weight_t0, &
                               this%realization%field,this%comm1)
  endif

  call NWTInitializeTimestep(this%realization)
  ! This will eventually call NWTEqDissPrecipSorb() so that the aqueous,
  ! precipitated, and sorbed amounts of the total bulk mass get adjusted.

end subroutine PMNWTInitializeTimestep

! ************************************************************************** !

subroutine PMNWTFinalizeTimestep(this)
  !
  ! Author: Jenn Frederick
  ! Date: 04/18/2019
  !

  use Variables_module, only : POROSITY
  use Material_module, only : MaterialGetAuxVarVecLoc
  use Material_Aux_module, only : POROSITY_BASE
  use Global_module

  implicit none

  class(pm_nwt_type) :: this
  PetscErrorCode :: ierr

  if (this%params%transient_porosity) then
    call VecCopy(this%realization%field%porosity_tpdt, &
                 this%realization%field%porosity_t,ierr);CHKERRQ(ierr)
    call RealizationUpdatePropertiesTS(this%realization)
    call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 POROSITY,POROSITY_BASE)
    call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                  this%realization%field%porosity_tpdt)
  endif

end subroutine PMNWTFinalizeTimestep

! ************************************************************************** !

subroutine PMNWTUpdateTimestep(this,update_dt, &
                               dt,dt_min,dt_max,iacceleration, &
                               num_newton_iterations,tfac, &
                               time_step_max_growth_factor)
  !
  ! Author: Jenn Frederick, Glenn Hammond
  ! Date: 05/27/2019
  !

  implicit none

  class(pm_nwt_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  PetscReal :: dtt
  PetscReal, parameter :: pert = 1.d-20

  if (update_dt .and. iacceleration /= 0) then
  
    ! this overrides the default setting of iacceleration = 5
    ! by default size(tfac) is 13, so this is safe
    dtt = dt
    iacceleration = size(tfac)
    if (num_newton_iterations <= iacceleration) then
      dtt = tfac(num_newton_iterations) * dt
    else
      dtt = this%controls%dt_cut * dt
    endif

    dtt = min(time_step_max_growth_factor*dt,dtt)
    if (dtt > dt_max) dtt = dt_max
    ! geh: see comment above under flow stepper
    dtt = max(dtt,dt_min)
    dt = dtt
  endif

end subroutine PMNWTUpdateTimestep

! ************************************************************************** !

subroutine PMNWTPreSolve(this)
  !
  ! Author: Jenn Frederick
  ! Date: 04/18/2019
  !

  use Global_module
  use Material_module
  use Data_Mediator_module

  implicit none

  class(pm_nwt_type) :: this

  PetscErrorCode :: ierr

  ! set densities and saturations to t+dt
  if (this%option%nflowdof > 0 .and. .not. this%params%steady_flow) then
    if (this%option%flow%transient_porosity) then
      ! weight material properties (e.g. porosity)
      call MaterialWeightAuxVars(this%realization%patch%aux%Material, &
                                 this%params%tran_weight_t1, &
                                 this%realization%field,this%comm1)
    endif
    call GlobalWeightAuxVars(this%realization,this%params%tran_weight_t1)
  else if (this%params%transient_porosity) then
    this%params%tran_weight_t1 = 1.d0
    call MaterialWeightAuxVars(this%realization%patch%aux%Material, &
                               this%params%tran_weight_t1, &
                               this%realization%field,this%comm1)
  endif

  if (this%realization%reaction_nw%use_log_formulation) then
    call VecCopy(this%realization%field%tran_xx, &
                 this%realization%field%tran_log_xx,ierr);CHKERRQ(ierr)
    call VecLog(this%realization%field%tran_log_xx,ierr);CHKERRQ(ierr)
  endif

  call DataMediatorUpdate(this%realization%tran_data_mediator_list, &
                          this%realization%field%tran_mass_transfer, &
                          this%realization%option)

end subroutine PMNWTPreSolve

! ************************************************************************** !

subroutine PMNWTPostSolve(this)
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  implicit none

  class(pm_nwt_type) :: this

  ! placeholder for now, nothing to do at the moment

end subroutine PMNWTPostSolve

! ************************************************************************** !

subroutine PMNWTUpdateSolution(this)
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  use Condition_module
  use Integral_Flux_module

  implicit none

  class(pm_nwt_type) :: this

  call TranConditionUpdate(this%realization%transport_conditions, &
                           this%realization%option)

  if (associated(this%realization%uniform_velocity_dataset)) then
    call RealizUpdateUniformVelocity(this%realization)
  endif

  if (this%realization%option%compute_mass_balance_new) then
    call NWTUpdateMassBalance(this%realization)
  endif

  if (this%option%transport%store_fluxes) then
    call IntegralFluxUpdate(this%realization%patch%integral_flux_list, &
                            this%realization%patch%internal_tran_fluxes, &
                            this%realization%patch%boundary_tran_fluxes, &
                            INTEGRATE_TRANSPORT,this%option)
  endif

end subroutine PMNWTUpdateSolution

! ************************************************************************** !

subroutine PMNWTResidual(this,snes,xx,r,ierr)
  !
  ! Author: Jenn Frederick
  ! Date: 03/14/13
  !

  implicit none

  class(pm_nwt_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  call NWTResidual(snes,xx,r,this%realization,ierr)

end subroutine PMNWTResidual

! ************************************************************************** !

subroutine PMNWTJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Jenn Frederick
  ! Date: 05/14/2019
  !

  implicit none

  class(pm_nwt_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  call NWTJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMNWTJacobian

! ************************************************************************** !

subroutine PMNWTCheckConvergence(this,snes,it,xnorm,unorm,fnorm,reason,ierr)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 05/14/2019, refactored 01/10/2023
  !
  use Convergence_module
  use Field_module
  use Option_module
  use Grid_module

  implicit none

  class(pm_nwt_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: out_string, string, rsn_string
  character(len=6) :: aR_str, sR_str, rUP_str, aUP_str
  character(len=4) :: rsn
  character(len=2) :: sign_str
  Vec :: dX
  PetscReal, pointer :: dX_p(:)    ! solution update
  PetscReal, pointer :: X0_p(:)    ! previous solution
  PetscReal, pointer :: X1_p(:)    ! current solution
  PetscReal, pointer :: r_p(:)     ! current residual R
  PetscReal, pointer :: accum_p(:) ! current accumulation term in R
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  PetscReal :: tempreal
  PetscReal :: max_relative_change
  PetscReal :: max_absolute_change
  PetscReal :: max_scaled_residual
  PetscReal :: max_absolute_residual
  PetscReal :: species_max_relative_update(this%option%ntrandof)
  PetscReal :: species_max_absolute_update(this%option%ntrandof)
  PetscReal :: species_max_scaled_residual(this%option%ntrandof)
  PetscReal :: species_max_absolute_residual(this%option%ntrandof)
  PetscReal :: species_max(this%option%ntrandof,4)
  PetscInt :: loc_max_scaled_residual
  PetscInt :: loc_max_abs_residual
  PetscInt :: loc_max_abs_update
  PetscInt :: loc_max_rel_update
  PetscInt :: i 
  PetscInt :: local_id, offset, idof, index
  PetscBool :: converged_flag
  
  PetscBool :: idof_is_positive(this%option%ntrandof, &
                                         this%realization%patch%grid%nlmax)
  PetscBool :: idof_cnvgd_due_to_rel_update(this%option%ntrandof, &
                                            this%realization%patch%grid%nlmax)
  PetscBool :: idof_cnvgd_due_to_abs_update(this%option%ntrandof, &
                                            this%realization%patch%grid%nlmax)
  PetscBool :: idof_cnvgd_due_to_update(this%option%ntrandof)
  PetscBool :: idof_cnvgd_due_to_scaled_res(this%option%ntrandof, &
                                            this%realization%patch%grid%nlmax)
  PetscBool :: idof_cnvgd_due_to_abs_res(this%option%ntrandof, &
                                         this%realization%patch%grid%nlmax)  
  PetscBool :: idof_cnvgd_due_to_residual(this%option%ntrandof)
  PetscBool :: rsn_rUP(this%option%ntrandof)
  PetscBool :: rsn_aUP(this%option%ntrandof)
  PetscBool :: rsn_aR(this%option%ntrandof)
  PetscBool :: rsn_sR(this%option%ntrandof)
  PetscBool :: rsn_pos(this%option%ntrandof)
  PetscBool :: idof_cnvgd(this%option%ntrandof,5)

  PetscMPIInt :: mpi_int1, mpi_int2
  PetscInt, parameter :: MAX_INDEX = 5

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field

  ! check if scaling from PMNWTCheckUpdatePre() requires ts cut
  if (this%controls%scaling_cut_dt) then
    option%convergence = CONVERGENCE_CUT_TIMESTEP
    reason = -88
    this%controls%scaling_cut_dt = PETSC_FALSE
    return
  endif

  if (it == 0) then
    option%converged = PETSC_FALSE
    option%convergence = CONVERGENCE_KEEP_ITERATING
    reason = 0
    return
  endif

  converged_flag = PETSC_FALSE
  idof_is_positive = PETSC_FALSE
  idof_cnvgd_due_to_update = PETSC_FALSE
  idof_cnvgd_due_to_rel_update = PETSC_FALSE
  idof_cnvgd_due_to_abs_update = PETSC_FALSE
  idof_cnvgd_due_to_residual = PETSC_FALSE
  idof_cnvgd_due_to_scaled_res = PETSC_FALSE
  idof_cnvgd_due_to_abs_res = PETSC_FALSE
  idof_cnvgd = PETSC_FALSE
  rsn_aR = PETSC_FALSE
  rsn_sR = PETSC_FALSE
  rsn_rUP = PETSC_FALSE
  rsn_aUP = PETSC_FALSE
  rsn_pos = PETSC_FALSE
  species_max_relative_update = 0.d0
  species_max_absolute_update = 0.d0
  species_max_absolute_residual = 0.d0
  species_max_scaled_residual = 0.d0
  species_max = 0.d0

  !call ConvergenceTest(snes,it,xnorm,unorm,fnorm,reason, &
  !                     this%realization%patch%grid, &
  !                     this%option,this%solver,ierr)
  ierr = 0
  call SNESGetSolutionUpdate(this%solver%snes,dX,ierr);CHKERRQ(ierr)

  call VecGetArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%tran_yy,X0_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%tran_xx,X1_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%tran_r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%tran_accum,accum_p,ierr);CHKERRQ(ierr)

  max_relative_change = maxval(dabs(dX_p(:)/X0_p(:)))
  max_absolute_change = maxval(dabs(dX_p(:)))
  max_absolute_residual = maxval(dabs(r_p(:)))
  max_scaled_residual = maxval(dabs(r_p(:)/accum_p(:)))

  loc_max_scaled_residual = maxloc(dabs(r_p(:)/accum_p(:)),1)
  loc_max_abs_residual = maxloc(dabs(r_p(:)),1)
  loc_max_rel_update = maxloc(dabs(dX_p(:)/X0_p(:)),1)
  loc_max_abs_update = maxloc(dabs(dX_p(:)),1)

  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%ntrandof
    do idof = 1, option%ntrandof
      index = idof + offset
    !-----------------------------------------------------------------
    !---solution-is-positive------------------------------------------
      tempreal = X1_p(index)
      if (tempreal > 0.d0) then
        idof_is_positive(idof,local_id) = PETSC_TRUE
      endif
    !-----------------------------------------------------------------
    !---relative-update-----------------------------------------------
      tempreal = dabs((dX_p(index))/X0_p(index))
      if (tempreal < this%controls%itol_rel_update(idof)) then
        idof_cnvgd_due_to_rel_update(idof,local_id) = PETSC_TRUE
      endif 
      species_max_relative_update(idof) = &
          max(species_max_relative_update(idof),tempreal)
    !-----------------------------------------------------------------
    !---absolute-update-----------------------------------------------
      tempreal = dabs(dX_p(index))
      if (tempreal < this%controls%itol_abs_update(idof)) then
        idof_cnvgd_due_to_abs_update(idof,local_id) = PETSC_TRUE
      endif 
      species_max_absolute_update(idof) = &
          max(species_max_absolute_update(idof),tempreal)
    !-----------------------------------------------------------------
    !---absolute-residual---------------------------------------------
      tempreal = dabs(r_p(index))
      if (tempreal < this%controls%itol_abs_res(idof)) then
        idof_cnvgd_due_to_abs_res(idof,local_id) = PETSC_TRUE
      endif
      species_max_absolute_residual(idof) = &
          max(species_max_absolute_residual(idof),tempreal)
    !-----------------------------------------------------------------
    !---scaled-residual-----------------------------------------------
      tempreal = dabs(r_p(index)/accum_p(index))
      if (tempreal < this%controls%itol_scaled_res(idof)) then
        idof_cnvgd_due_to_scaled_res(idof,local_id) = PETSC_TRUE
      endif
      species_max_scaled_residual(idof) = &
          max(species_max_scaled_residual(idof),tempreal)
    !-----------------------------------------------------------------
    enddo 
  enddo
  
  call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%tran_yy,X0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%tran_xx,X1_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%tran_r,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%tran_accum,accum_p,ierr);CHKERRQ(ierr)


  do idof = 1, option%ntrandof
    species_max(idof,1) = species_max_absolute_update(idof)
    if (all(idof_cnvgd_due_to_abs_update(idof,:))) then
      rsn_aUP(idof) = PETSC_TRUE
    endif
    species_max(idof,2) = species_max_relative_update(idof)
    if (all(idof_cnvgd_due_to_rel_update(idof,:))) then
      rsn_rUP(idof) = PETSC_TRUE
      rUP_str = ' *rUP:'
    endif
    species_max(idof,3) = species_max_absolute_residual(idof)
    if (all(idof_cnvgd_due_to_abs_res(idof,:))) then
      rsn_aR(idof) = PETSC_TRUE
      aR_str  = '  *aR:'
    endif
    species_max(idof,4) = species_max_scaled_residual(idof)
    if (all(idof_cnvgd_due_to_scaled_res(idof,:))) then
      rsn_sR(idof) = PETSC_TRUE
      sR_str  = '  *sR:'
    endif
    if (all(idof_is_positive(idof,:))) then
      rsn_pos(idof) = PETSC_TRUE
      sign_str  = ' +'
    endif

    if (rsn_aUP(idof) .or. rsn_rUP(idof)) then
      idof_cnvgd_due_to_update(idof) = PETSC_TRUE
    endif
    if (rsn_aR(idof) .or. rsn_sR(idof)) then
      idof_cnvgd_due_to_residual(idof) = PETSC_TRUE
    endif
    idof_cnvgd(idof,1)=rsn_aUP(idof)
    idof_cnvgd(idof,2)=rsn_rUP(idof)
    idof_cnvgd(idof,3)=rsn_aR(idof)
    idof_cnvgd(idof,4)=rsn_sR(idof)
    idof_cnvgd(idof,5)=rsn_pos(idof)
  enddo

  if (all(idof_cnvgd_due_to_update) .and. all(idof_cnvgd_due_to_residual) &
      .and. all(idof_is_positive)) then
    converged_flag = PETSC_TRUE
  else
    converged_flag = PETSC_FALSE
  endif

  mpi_int1 = option%ntrandof*5
  mpi_int2 = option%ntrandof*4
  call MPI_Allreduce(MPI_IN_PLACE,converged_flag,ONE_INTEGER, &
                     MPI_LOGICAL,MPI_LAND,option%mycomm,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,idof_cnvgd,mpi_int1, &
                     MPI_LOGICAL,MPI_LAND,option%mycomm,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,species_max,mpi_int2, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)

  rUP_str = '  rUP:'
  aUP_str = '  aUP:'
  aR_str  = '   aR:'
  sR_str  = '   sR:'
  sign_str = ' -'
  rsn_string = ''
  
  if (all(idof_cnvgd(:,1))) then
    aUP_str = ' *aUP:'
    rsn_string = trim(rsn_string) // ' aUP,'
  endif
  if (all(idof_cnvgd(:,2))) then
    rUP_str = ' *rUP:'
    rsn_string = trim(rsn_string) // ' rUP,' 
  endif
  if (all(idof_cnvgd(:,3))) then
    aR_str = '  *aR:'
    rsn_string = trim(rsn_string) // ' aR,'
  endif
  if (all(idof_cnvgd(:,4))) then
    sR_str = '  *sR:'
    rsn_string = trim(rsn_string) // ' sR,'
  endif  
  if (all(idof_cnvgd(:,5))) then
    sign_str = ' +'
    rsn_string = trim(rsn_string) // ' +'
  endif
  
  if (converged_flag .eqv. PETSC_TRUE) then
    option%converged = PETSC_TRUE
    option%convergence = CONVERGENCE_CONVERGED
    reason = 999
  else  ! means ITOL_* tolerances were satisfied
    ! means ITOL_* tolerances were not satisfied:
    option%converged = PETSC_FALSE
    option%convergence = CONVERGENCE_KEEP_ITERATING
    reason = 0
    if (it >= this%controls%max_newton_iterations) then
      option%convergence = CONVERGENCE_CUT_TIMESTEP
      reason = -88
    endif
  endif

  rsn = 'rsn:'
  write(out_string,'(i4,a6,es10.3,a6,es10.3,a6,es10.3,a6,es10.3,a5,i4,a2)') &
                    it, aR_str, maxval(species_max(:,3)), &
                    sR_str, maxval(species_max(:,4)), &
                    rUP_str, maxval(species_max(:,2)), &
                    aUP_str, maxval(species_max(:,1)), rsn, reason, &
                    sign_str
  call PrintMsg(option,out_string)

  if (this%controls%verbose_newton) then
    ! only works for single core
    write(out_string,'(a4,a6,i10,a6,i10,a6,i10,a6,i10)') &
        'cell', aR_str, loc_max_abs_residual, sR_str, &
        loc_max_scaled_residual, rUP_str, loc_max_rel_update, aUP_str, &
        loc_max_abs_update
    call PrintMsg(option,out_string)

    out_string = "     update converged (T/F) = "
    do i = 1,option%ntrandof
      write(string,'(L2)') idof_cnvgd_due_to_update(i)
      out_string = trim(out_string) // trim(string)
    enddo
    call PrintMsg(option,out_string)

    out_string = "     residual converged (T/F) = "
    do i = 1,option%ntrandof
      write(string,'(L2)') idof_cnvgd_due_to_residual(i)
      out_string = trim(out_string) // trim(string)
    enddo
    call PrintMsg(option,out_string)
  endif
  
  ! finalize convergence and reset option flags:
  if (option%convergence == CONVERGENCE_CONVERGED) then
    out_string = ' NW TRANS converged!   Rsn ->' // trim(rsn_string)
    call PrintMsg(option,out_string)
    option%convergence = CONVERGENCE_OFF
    option%converged = PETSC_FALSE
  endif 

end subroutine PMNWTCheckConvergence

! ************************************************************************** !

subroutine PMNWTCheckUpdatePre(this,snes,X,dX,changed,ierr)
  !
  ! Checks the solution update vector prior to updating the solution.
  !
  ! Author: Glenn Hammond, Jenn Frederick
  ! Date: 05/28/2019
  !

  use Realization_Subsurface_class
  use Grid_module
  use Option_module

  implicit none

  class(pm_nwt_type) :: this
  SNES :: snes
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr

  PetscReal, pointer :: X_p(:)  ! CURRENT SOLUTION
  PetscReal, pointer :: dX_p(:) ! SOLUTION UPDATE STEP
  type(grid_type), pointer :: grid
  class(reaction_nw_type), pointer :: reaction_nw
  type(option_type), pointer :: option
  PetscReal :: ratio, min_ratio, a, b
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i, k
  PetscReal,parameter :: TOL=1.0d-20

  grid => this%realization%patch%grid
  reaction_nw => this%realization%reaction_nw
  option => this%realization%option

  changed = PETSC_FALSE

  call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)

  if (Initialized(reaction_nw%params%truncated_concentration)) then
    do k = 1,size(X_p)
      if (X_p(k) < reaction_nw%params%truncated_concentration) then
        X_p(k) = reaction_nw%params%truncated_concentration
        dX_p(k) = 0.0d0
      else
        dX_p(k) = min(dX_p(k),X_p(k) - &
                      reaction_nw%params%truncated_concentration)
      endif
    enddo
  endif

  call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)

end subroutine PMNWTCheckUpdatePre

! ************************************************************************** !

subroutine PMNWTCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                X1_changed,ierr)
  !
  ! Checks convergence after the solution update
  !
  ! Author: Jenn Frederick, Glenn Hammond
  ! Date: 05/27/2019
  !

  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Output_EKG_module

  implicit none

  class(pm_nwt_type) :: this
  SNES :: snes
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr


end subroutine PMNWTCheckUpdatePost

! ************************************************************************** !

function PMNWTAcceptSolution(this)
  !
  ! Right now this does nothing, but if it did spit out PETSC_FALSE, it would
  ! prompt a time step cut and start the solution process over again.
  !
  ! Author: Jenn Frederick
  ! Date: 03/14/13
  !

  implicit none

  class(pm_nwt_type) :: this

  PetscBool :: PMNWTAcceptSolution

  ! no nothing
  PMNWTAcceptSolution = PETSC_TRUE

end function PMNWTAcceptSolution

! ************************************************************************** !

recursive subroutine PMNWTFinalizeRun(this)
  !
  ! Finalizes the time stepping
  !
  ! Author: Jenn Frederick
  ! Date: 05/28/2019
  !

  implicit none

  class(pm_nwt_type) :: this

  ! placeholder for now, doesn't do anything at the moment.

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMNWTFinalizeRun

! ************************************************************************** !

subroutine PMNWTTimeCut(this)
  !
  ! Resets arrays for a time step cut.
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  use Option_module
  use Field_module
  use Global_module

  implicit none

  class(pm_nwt_type) :: this

  class(realization_subsurface_type), pointer :: realization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  PetscErrorCode :: ierr

  realization => this%realization
  field => realization%field
  option => realization%option

  this%option%tran_dt = this%option%dt
  if (this%option%nflowdof > 0 .and. .not. this%params%steady_flow) then
    call this%SetTranWeights()
  endif

  ! copy previous solution back to current solution
  call VecCopy(field%tran_yy,field%tran_xx,ierr);CHKERRQ(ierr)

  ! set densities and saturations to t+dt
  if (realization%option%nflowdof > 0) then
    call GlobalWeightAuxVars(realization, &
                             realization%option%transport%tran_weight_t1)
  endif

end subroutine PMNWTTimeCut

! ************************************************************************** !

subroutine PMNWTSetTranWeights(this)
  !
  ! Sets the weights at t0 or t1 for transport
  !
  ! Author: Glenn Hammond
  ! Date: 01/17/11; 04/03/13
  !

  implicit none

  class(pm_nwt_type) :: this

  PetscReal :: flow_dt
  PetscReal :: flow_t0
  PetscReal :: flow_t1

  ! option%tran_time is the time at beginning of transport step
  flow_t0 = this%realization%patch%aux%Global%time_t
  flow_t1 = this%realization%patch%aux%Global%time_tpdt
  flow_dt = flow_t1-flow_t0
  this%params%tran_weight_t0 = max(0.d0,(this%option%time-flow_t0)/flow_dt)
  this%params%tran_weight_t1 = min(1.d0, &
                            (this%option%time+this%option%tran_dt-flow_t0)/ &
                            flow_dt)

end subroutine PMNWTSetTranWeights

! ************************************************************************** !

subroutine PMNWTComputeMassBalance(this,mass_balance_array)
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  implicit none

  class(pm_nwt_type) :: this
  PetscReal :: mass_balance_array(:)

#ifndef SIMPLIFY
  ! passing in dummy -999 for max_size doesn't do anything, so why even call it?
  !call NWTComputeMassBalance(this%realization,-999,mass_balance_array)
#endif

end subroutine PMNWTComputeMassBalance

! ************************************************************************** !

subroutine PMNWTUpdateAuxVars(this)
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  implicit none

  class(pm_nwt_type) :: this
                                      !  cells      bcs
  call NWTUpdateAuxVars(this%realization,PETSC_TRUE,PETSC_FALSE)

end subroutine PMNWTUpdateAuxVars

! ************************************************************************** !

subroutine PMNWTMaxChange(this)
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  implicit none

  class(pm_nwt_type) :: this

  !TODO(jenn)
  print *, 'PMNWTMaxChange not implemented.'
  stop

end subroutine PMNWTMaxChange

! ************************************************************************** !

subroutine PMNWTCheckpointBinary(this,viewer)
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  use Option_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Field_module
  use Discretization_module
  use Variables_module, only : NWT_AUXILIARY

  implicit none

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: pm_base_header_type
      implicit none
      PetscBag :: bag
      class(pm_base_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData

  class(pm_nwt_type) :: this
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  Vec :: global_vec
  PetscInt :: i

  class(pm_base_header_type), pointer :: header
  type(pm_base_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize

  realization => this%realization
  option => realization%option
  field => realization%field
  discretization => realization%discretization

  global_vec = PETSC_NULL_VEC

  bagsize = size(transfer(dummy_header,dummy_char))

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%ndof,0,"ndof","",ierr);CHKERRQ(ierr)
  header%ndof = option%ntrandof
  call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

  if (option%ntrandof > 0) then
    call VecView(field%tran_xx,viewer,ierr);CHKERRQ(ierr)

    if (global_vec == PETSC_NULL_VEC) then
      call DiscretizationCreateVector(discretization,ONEDOF, &
                                      global_vec,GLOBAL,option)
    endif
    ! auxiliary data for reactions (e.g. cumulative mass)
    if (realization%reaction_nw%params%nauxiliary > 0) then
      do i = 1, realization%reaction_nw%params%nauxiliary
        call RealizationGetVariable(realization,global_vec, &
                                    NWT_AUXILIARY,i)
        call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
      enddo
    endif
  endif

  if (global_vec /= PETSC_NULL_VEC) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif

end subroutine PMNWTCheckpointBinary

! ************************************************************************** !

subroutine PMNWTCheckpointHDF5(this,pm_grp_id)
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  use hdf5
  use Option_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Field_module
  use Discretization_module
  use Variables_module, only : NWT_AUXILIARY
  use Checkpoint_module, only: CheckPointWriteIntDatasetHDF5
  use HDF5_module, only : HDF5WriteDataSetFromVec

  implicit none

  class(pm_nwt_type) :: this
  integer(HID_T) :: pm_grp_id

  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)

  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: i
  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  discretization => realization%discretization

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  dataset_name = "NDOF" // CHAR(0)
  int_array(1) = option%ntrandof
  call CheckPointWriteIntDatasetHDF5(pm_grp_id, dataset_name, dataset_rank, &
                                     dims, start, length, stride, &
                                     int_array, option)

  if (option%ntrandof > 0) then
    call DiscretizationCreateVector(discretization, NTRANDOF, &
                                    natural_vec, NATURAL, option)
    call DiscretizationGlobalToNatural(discretization, &
                                       field%tran_xx, &
                                       natural_vec, NTRANDOF)
    dataset_name = "Primary_Variable" // CHAR(0)
    call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
                                 pm_grp_id, H5T_NATIVE_DOUBLE)
    call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)

    ! auxiliary data for reactions (e.g. cumulative mass)
    if (realization%reaction_nw%params%nauxiliary> 0) then
      call DiscretizationCreateVector(discretization, ONEDOF, &
                                      global_vec, GLOBAL, option)
      call DiscretizationCreateVector(discretization, ONEDOF, &
                                      natural_vec, NATURAL, option)
      do i = 1, realization%reaction_nw%params%nauxiliary
        call RealizationGetVariable(realization,global_vec, &
                                    NWT_AUXILIARY,i)
        call DiscretizationGlobalToNatural(discretization, &
                                           global_vec, natural_vec, ONEDOF)
        write(dataset_name,*) i
        dataset_name = 'NWT_auxiliary_' // trim(adjustl(dataset_name))
        call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
                                     pm_grp_id, H5T_NATIVE_DOUBLE)
      enddo
      call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
    endif
  endif

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)

end subroutine PMNWTCheckpointHDF5

! ************************************************************************** !

subroutine PMNWTRestartBinary(this,viewer)
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  use Option_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Patch_module
  use Variables_module, only : NWT_AUXILIARY

  implicit none

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: pm_base_header_type
      implicit none
      PetscBag :: bag
      class(pm_base_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData

  class(pm_nwt_type) :: this
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  Vec :: global_vec
  PetscInt :: i

  class(pm_base_header_type), pointer :: header
  type(pm_base_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize

  realization => this%realization
  option => realization%option
  field => realization%field
  discretization => realization%discretization

  global_vec = PETSC_NULL_VEC

  bagsize = size(transfer(dummy_header,dummy_char))

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%ndof,0,"ndof","",ierr);CHKERRQ(ierr)
  call PetscBagLoad(viewer,bag,ierr);CHKERRQ(ierr)
  option%ntrandof = header%ndof

  call VecLoad(field%tran_xx,viewer,ierr);CHKERRQ(ierr)
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)
  call VecCopy(field%tran_xx,field%tran_yy,ierr);CHKERRQ(ierr)

  if (global_vec == PETSC_NULL_VEC) then
    call DiscretizationCreateVector(discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
  endif

  ! auxiliary data for reactions (e.g. cumulative mass)
  if (realization%reaction_nw%params%nauxiliary> 0) then
    do i = 1, realization%reaction_nw%params%nauxiliary
      call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
      call RealizationSetVariable(realization,global_vec, &
                                    GLOBAL,NWT_AUXILIARY,i)
    enddo
  endif

  if (global_vec /= PETSC_NULL_VEC) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif

  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)
  call NWTUpdateAuxVars(realization,PETSC_FALSE,PETSC_TRUE)
  call PMNWTUpdateSolution(this)

end subroutine PMNWTRestartBinary

! ************************************************************************** !

subroutine PMNWTRestartHDF5(this,pm_grp_id)
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  use hdf5
  use Option_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Field_module
  use Discretization_module
  use Variables_module, only : NWT_AUXILIARY
  use Checkpoint_module, only: CheckPointReadIntDatasetHDF5
  use HDF5_module, only : HDF5ReadDataSetInVec

  implicit none

  class(pm_nwt_type) :: this
  integer(HID_T) :: pm_grp_id

  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)

  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  Vec :: global_vec
  Vec :: natural_vec
  Vec :: local_vec
  PetscInt :: i
  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  discretization => realization%discretization

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  dataset_name = "NDOF" // CHAR(0)
  int_array(1) = option%ntrandof
  call CheckPointReadIntDatasetHDF5(pm_grp_id, dataset_name, dataset_rank, &
                                    dims, start, length, stride, &
                                    int_array, option)
  option%ntrandof = int_array(1)

  if (option%ntrandof > 0) then
    call DiscretizationCreateVector(discretization, NTRANDOF, &
                                    natural_vec, NATURAL, option)
    dataset_name = "Primary_Variable" // CHAR(0)
    call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
                              pm_grp_id, H5T_NATIVE_DOUBLE)
    call DiscretizationNaturalToGlobal(discretization, &
                                       natural_vec, &
                                       field%tran_xx,NTRANDOF)
    call DiscretizationGlobalToNatural(discretization, &
                                       field%tran_xx, &
                                       field%tran_xx_loc, NTRANDOF)
    call VecCopy(field%tran_xx,field%tran_yy,ierr);CHKERRQ(ierr)
    call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)

    call DiscretizationCreateVector(discretization, ONEDOF, &
                                    global_vec, GLOBAL, option)
    call DiscretizationCreateVector(discretization, ONEDOF, &
                                    natural_vec, NATURAL, option)
    call DiscretizationCreateVector(discretization, ONEDOF, &
                                    local_vec, LOCAL, option)
    ! auxiliary data for reactions (e.g. cumulative mass)
    if (realization%reaction_nw%params%nauxiliary> 0) then
       do i = 1, realization%reaction_nw%params%nauxiliary
         write(dataset_name,*) i
         dataset_name = 'NWT_auxiliary_' // trim(adjustl(dataset_name))
         call HDF5ReadDataSetInVec(dataset_name,option,natural_vec, &
                                   pm_grp_id,H5T_NATIVE_DOUBLE)
         call DiscretizationNaturaltoGlobal(discretization,natural_vec, &
                                            global_vec,ONEDOF)
         call DiscretizationGlobaltoLocal(discretization,global_vec, &
                                          local_vec, ONEDOF)
         call RealizationSetVariable(realization,local_vec, &
                                     LOCAL,NWT_AUXILIARY,i)
      enddo
    endif
  endif

  call NWTUpdateAuxVars(realization,PETSC_FALSE,PETSC_TRUE)
  call PMNWTUpdateSolution(this)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)

end subroutine PMNWTRestartHDF5

! ************************************************************************** !

subroutine PMNWTInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  implicit none

  class(pm_nwt_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMNWTInputRecord

! ************************************************************************** !

subroutine PMNWTDestroy(this)
  !
  ! Destroys objects in the NW Transport process model.
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !

  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_nwt_type) :: this

  call DeallocateArray(this%controls%itol_rel_update)
  call DeallocateArray(this%controls%itol_scaled_res)
  call DeallocateArray(this%controls%itol_abs_res)
  call DeallocateArray(this%controls%cnvg_criteria_value)
  call DeallocateArray(this%controls%i_mapping)
  call DeallocateArray(this%controls%names_itol_rel_update)
  call DeallocateArray(this%controls%names_itol_scaled_res)
  call DeallocateArray(this%controls%names_itol_abs_res)

  if (associated(this%params%bh_material_ids)) then
    call DeallocateArray(this%params%bh_material_ids)
    call DeallocateArray(this%params%bh_material_names)
  endif
  if (associated(this%params%dirichlet_material_ids)) then
    call DeallocateArray(this%params%dirichlet_material_ids)
    call DeallocateArray(this%params%dirichlet_material_names)
  endif

  call PMBaseDestroy(this)
  call NWTDestroy(this%realization)

  nullify(this%comm1) ! already destroyed in realization

end subroutine PMNWTDestroy

! ************************************************************************** !

end module PM_NWT_class
