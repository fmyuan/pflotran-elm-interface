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
    PetscReal, pointer :: max_concentration_change(:)
    PetscReal, pointer :: max_volfrac_change(:)
    PetscReal :: volfrac_change_governor
    PetscReal :: cfl_governor
    PetscReal, pointer :: itol_rel_update(:)
    PetscReal, pointer :: itol_scaled_res(:)
    PetscReal, pointer :: itol_abs_res(:)
    PetscReal, pointer :: cnvg_criteria_value(:) 
    PetscInt, pointer :: i_mapping(:) 
    character(len=MAXWORDLENGTH), pointer :: names_itol_rel_update(:)
    character(len=MAXWORDLENGTH), pointer :: names_itol_scaled_res(:)
    character(len=MAXWORDLENGTH), pointer :: names_itol_abs_res(:)
    PetscInt :: max_newton_iterations
    PetscReal :: max_dlnC
    PetscReal :: dt_cut
    PetscBool :: check_post_converged
    PetscBool :: check_post_convergence
    PetscBool :: check_update
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
  nullify(nwt_pm%controls%max_concentration_change)
  nullify(nwt_pm%controls%max_volfrac_change)
  nwt_pm%controls%volfrac_change_governor = 1.d0
  nwt_pm%controls%cfl_governor = UNINITIALIZED_DOUBLE
  nullify(nwt_pm%controls%itol_rel_update)
  nullify(nwt_pm%controls%itol_scaled_res)
  nullify(nwt_pm%controls%itol_abs_res)
  nullify(nwt_pm%controls%cnvg_criteria_value)
  nullify(nwt_pm%controls%i_mapping)
  nullify(nwt_pm%controls%names_itol_rel_update)
  nullify(nwt_pm%controls%names_itol_scaled_res)
  nullify(nwt_pm%controls%names_itol_abs_res)
  nwt_pm%controls%max_newton_iterations = 10
  nwt_pm%controls%dt_cut = 0.5d0
  nwt_pm%controls%max_dlnC = 5.0d0
  nwt_pm%controls%check_post_converged = PETSC_FALSE
  ! set to true so that itol_relative_update always gets checked:
  nwt_pm%controls%check_post_convergence = PETSC_TRUE
  nwt_pm%controls%check_update = PETSC_TRUE
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
    case('CFL_GOVERNOR')
      call InputReadDouble(input,option,this%controls%cfl_governor)
      call InputErrorMsg(input,option,keyword,error_string)
    case('VOLUME_FRACTION_CHANGE_GOVERNOR')
      call InputReadDouble(input,option,this%controls%volfrac_change_governor)
      call InputErrorMsg(input,option,keyword,error_string)
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
  PetscReal, pointer :: temp_itol_rel_update(:)
  PetscReal, pointer :: temp_itol_scaled_res(:)
  PetscReal, pointer :: temp_itol_abs_res(:)
  PetscInt :: k, j
  character(len=MAXSTRINGLENGTH) :: error_string_ex
  character(len=MAXWORDLENGTH) :: word

  option => this%option

  allocate(temp_species_names(50))
  
  !found = PETSC_TRUE
  !call PMBaseReadSelectCase(this,input,keyword,found,error_string,option)
  !if (found) return
    
  found = PETSC_TRUE
  select case(trim(keyword))
!geh: these have not been implemented
!    case('NUMERICAL_JACOBIAN')
!      option%transport%numerical_derivatives = PETSC_TRUE
    !------------------------------------------------------------------------
    case('MAXIMUM_NUMBER_OF_ITERATIONS')
      error_string_ex = trim(error_string) // ',MAXIMUM_NUMBER_OF_ITERATIONS'
      call InputReadInt(input,option,this%controls%max_newton_iterations)
      call InputErrorMsg(input,option,'VALUE',error_string)
      this%solver%newton_max_iterations = this%controls%max_newton_iterations
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
  
  allocate(this%controls%max_concentration_change(this%params%nspecies))
  allocate(this%controls%max_volfrac_change(this%params%nspecies))
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
  type(coupler_type), pointer :: init_condition
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
  use Material_Aux_class, only : POROSITY_BASE 
  use Global_module

  implicit none
  
  class(pm_nwt_type) :: this
  PetscReal :: time  
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
  
  call NWTMaxChange(this%realization,this%controls%max_concentration_change)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dcmx= ",1pe12.4,"  dc/dt= ",1pe12.4, &
            &" [mol/s]")') &
      maxval(this%controls%max_concentration_change), &
      maxval(this%controls%max_concentration_change)/this%option%tran_dt
  endif
  if (this%option%print_file_flag) then  
    write(this%option%fid_out,&
            '("  --> max chng: dcmx= ",1pe12.4,"  dc/dt= ",1pe12.4, &
            &" [mol/s]")') &
      maxval(this%controls%max_concentration_change), &
      maxval(this%controls%max_concentration_change)/this%option%tran_dt
  endif
  
end subroutine PMNWTFinalizeTimestep

! ************************************************************************** !

subroutine PMNWTUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                               num_newton_iterations,tfac, &
                               time_step_max_growth_factor)
  ! 
  ! Author: Jenn Frederick, Glenn Hammond
  ! Date: 05/27/2019
  ! 

  implicit none
  
  class(pm_nwt_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  
  PetscReal :: dtt, uvf, dt_vf, dt_tfac, fac
  PetscInt :: ifac
  PetscReal, parameter :: pert = 1.d-20
    
  if (this%controls%volfrac_change_governor < 1.d0) then
    ! with volume fraction potentially scaling the time step
    if (iacceleration > 0) then
      fac = 0.5d0
      if (num_newton_iterations >= iacceleration) then
        fac = 0.33d0
        uvf = 0.d0
      else
        uvf = this%controls%volfrac_change_governor/ &
              (maxval(this%controls%max_volfrac_change)+pert)
      endif
      dtt = fac * dt * (1.d0 + uvf)
    else
      ifac = max(min(num_newton_iterations,size(tfac)),1)
      dt_tfac = tfac(ifac) * dt

      fac = 0.5d0
      uvf= this%controls%volfrac_change_governor/ &
           (maxval(this%controls%max_volfrac_change)+pert)
      dt_vf = fac * dt * (1.d0 + uvf)

      dtt = min(dt_tfac,dt_vf)
    endif
  else
    ! original implementation
    ! this overrides the default setting of iacceleration = 5
    ! by default size(tfac) is 13, so this is safe
    dtt = dt
    iacceleration = size(tfac)
    if (num_newton_iterations <= iacceleration) then
      dtt = tfac(num_newton_iterations) * dt
    else
      dtt = this%controls%dt_cut * dt
    endif
  endif

  dtt = min(time_step_max_growth_factor*dt,dtt)
  if (dtt > dt_max) dtt = dt_max
  ! geh: see comment above under flow stepper
  dtt = max(dtt,dt_min)
  dt = dtt

  call RealizationLimitDTByCFL(this%realization,this%controls%cfl_governor, &
                               dt,dt_max)
  
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
  ! Author: Jenn Frederick
  ! Date: 05/14/2019
  ! 
  use Convergence_module

  implicit none

  class(pm_nwt_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  PetscReal, pointer :: dC_p(:) ! SOLUTION UPDATE STEP
  PetscReal, pointer :: C_p(:) ! CURRENT SOLUTION 
  Vec :: update_vec, curr_solution_vec
  PetscReal :: max_relative_change
  PetscReal :: max_update
  PetscBool :: converged_due_to_rel_update
  PetscInt :: converged_flag, temp_int

  call ConvergenceTest(snes,it,xnorm,unorm,fnorm,reason, &
                       this%realization%patch%grid, &
                       this%option,this%solver,ierr)

end subroutine PMNWTCheckConvergence

! ************************************************************************** !

subroutine PMNWTCheckUpdatePre(this,snes,X,dX,changed,ierr)
  ! 
  ! In the case of the log formulation, ensures that the update
  ! vector does not exceed a prescribed tolerance
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
  
  PetscReal, pointer :: C_p(:)  ! CURRENT SOLUTION
  PetscReal, pointer :: dC_p(:) ! SOLUTION UPDATE STEP
  type(grid_type), pointer :: grid
  class(reaction_nw_type), pointer :: reaction_nw
  PetscReal :: ratio, min_ratio
  PetscReal, parameter :: min_allowable_scale = 1.d-10
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i, n, k 
  
  grid => this%realization%patch%grid
  reaction_nw => this%realization%reaction_nw
  
  call VecGetArrayF90(X,C_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(dX,dC_p,ierr);CHKERRQ(ierr)
  
  if (reaction_nw%use_log_formulation) then
    ! C and dC are actually lnC and dlnC
    dC_p = dsign(1.d0,dC_p)*min(dabs(dC_p),this%controls%max_dlnC)
    ! at this point, it does not matter whether "changed" is set to true, 
    ! since it is not checkied in PETSc.  Thus, I don't want to spend 
    ! time checking for changes and performing an allreduce for log 
    ! formulation.
    if (Initialized(reaction_nw%params%truncated_concentration)) then
      dC_p = min(C_p-log(reaction_nw%params%truncated_concentration),dC_p)
    endif
  else
        
    if (Initialized(reaction_nw%params%truncated_concentration)) then
      do k = 1,size(C_p)
        if (C_p(k) < reaction_nw%params%truncated_concentration) then
          C_p(k) = reaction_nw%params%truncated_concentration
          dC_p(k) = 0.0d0 
        else
          dC_p(k) = min(dC_p(k),C_p(k) - &
                        reaction_nw%params%truncated_concentration)
        endif 
      enddo 
    else
      ! C^p+1 = C^p - dC^p
      ! if dC is positive and abs(dC) larger than C
      ! we need to scale the update
      
      ! compute smallest ratio of C to dC
#if 0
      min_ratio = 1.d0/maxval(dC_p/C_p)
#else
      min_ratio = MAX_DOUBLE ! large number
      k = 0
      do i = 1, n
        if (C_p(i) <= dC_p(i)) then
          WRITE(*,*)  ' i =', i, '  C_p(i) =', C_p(i), '  dC_p(i) =', dC_p(i)
          ratio = abs(C_p(i)/dC_p(i))
          if (ratio < min_ratio) then
            min_ratio = ratio
            k = i
          endif
        endif
      enddo
#endif
      ratio = min_ratio
      !WRITE(*,*)  ' location of min_ratio =', k, '   min_ratio =', min_ratio
    
      ! get global minimum
      call MPI_Allreduce(ratio,min_ratio,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                         MPI_MIN,this%realization%option%mycomm,ierr)
                       
      ! scale if necessary
      if (min_ratio < 1.d0) then
        if (min_ratio < this%realization%option%min_allowable_scale) then
          write(string,'(es10.3)') min_ratio
          string = 'The update of primary species concentration is being ' // &
            'scaled by a very small value (i.e. ' // &
            trim(adjustl(string)) // &
            ') to prevent negative concentrations.  This value is too ' // &
            'small and will likely cause the solver to mistakenly ' // &
            'converge based on the infinity norm of the update vector. ' // &
            'In this case, it is recommended that you use the ' // &
            'LOG_FORMULATION keyword or TRUNCATE_CONCENTRATION <float> in &
            &the NUCLEAR_WASTE_TRANSPORT block).'
          this%realization%option%io_buffer = string
          call PrintErrMsgToDev(this%realization%option, 'send your input &
                                &deck if that does not work')
        endif
        ! scale by 0.99 to make the update slightly smaller than the min_ratio
        dC_p = dC_p*min_ratio*0.99d0
        changed = PETSC_TRUE
      endif
    endif
  endif

  call VecRestoreArrayF90(X,C_p,ierr);CHKERRQ(ierr) 
  call VecRestoreArrayF90(dX,dC_p,ierr);CHKERRQ(ierr)

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
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  PetscReal, pointer :: C0_p(:)
  PetscReal, pointer :: dC_p(:)
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:)
  PetscBool :: idof_cnvgd_due_to_update(this%option%ntrandof)
  PetscBool :: idof_cnvgd_due_to_residual(this%option%ntrandof)
  PetscBool :: idof_cnvgd_due_to_rel_update(this%option%ntrandof,this%realization%patch%grid%nlmax)
  PetscBool :: idof_cnvgd_due_to_scaled_res(this%option%ntrandof,this%realization%patch%grid%nlmax)
  PetscBool :: idof_cnvgd_due_to_abs_res(this%option%ntrandof,this%realization%patch%grid%nlmax)
  PetscReal :: max_relative_change
  PetscReal :: max_scaled_residual
  PetscReal :: max_absolute_change
  PetscReal :: max_absolute_residual
  PetscReal :: min_C0, min_C_prev
  PetscInt :: loc_max_scaled_residual
  PetscInt :: loc_max_abs_residual
  PetscInt :: loc_max_rel_update
  PetscReal :: residual_at_max
  PetscReal :: accum_at_max
  PetscReal :: update_at_max
  PetscReal :: soln_at_max
  PetscReal :: species_max_relative_change(this%option%ntrandof)
  PetscReal :: species_max_scaled_residual(this%option%ntrandof)
  PetscReal :: species_max_absolute_change(this%option%ntrandof)
  PetscReal :: species_max_absolute_residual(this%option%ntrandof)
  PetscInt :: converged_flag
  PetscInt :: temp_int, i
  PetscInt :: newton_iter_number
  PetscReal :: max_relative_change_by_dof(this%option%ntrandof)
  PetscReal :: global_max_rel_change_by_dof(this%option%ntrandof)
  PetscMPIInt :: mpi_int
  PetscInt :: local_id, offset, idof, index
  PetscReal :: tempreal
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  call SNESGetIterationNumber(this%solver%snes,newton_iter_number,ierr)
  
  converged_flag = 0
  if (this%controls%check_post_convergence) then
    idof_cnvgd_due_to_rel_update = PETSC_FALSE
    idof_cnvgd_due_to_scaled_res = PETSC_FALSE
    idof_cnvgd_due_to_abs_res = PETSC_FALSE
    idof_cnvgd_due_to_update = PETSC_FALSE
    idof_cnvgd_due_to_residual = PETSC_FALSE

    call VecGetArrayReadF90(dX,dC_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(X0,C0_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%tran_r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%tran_accum,accum_p,ierr);CHKERRQ(ierr)

    min_C0 = minval(dabs(C0_p(:)))

    max_relative_change = maxval(dabs(dC_p(:)/C0_p(:)))  
    max_scaled_residual = maxval(dabs(r_p(:)/accum_p(:)))
    max_absolute_change = maxval(dabs(dC_p(:)))
    max_absolute_residual = maxval(dabs(r_p(:)))

    loc_max_scaled_residual = maxloc(dabs(r_p(:)/accum_p(:)),1)
    loc_max_abs_residual = maxloc(dabs(r_p(:)),1)
    loc_max_rel_update = maxloc(dabs(dC_p(:)/C0_p(:)),1)

    residual_at_max = dabs(r_p(loc_max_scaled_residual))
    accum_at_max = dabs(accum_p(loc_max_scaled_residual))
    update_at_max = dabs(dC_p(loc_max_rel_update))
    soln_at_max = dabs(C0_p(loc_max_rel_update))

    do local_id = 1, grid%nlmax
      offset = (local_id-1)*option%ntrandof
      do idof = 1, option%ntrandof
        index = idof + offset
      !-----------------------------------------------------------------
        idof_cnvgd_due_to_rel_update(idof,local_id) = PETSC_FALSE
        tempreal = dabs((dC_p(index))/C0_p(index))
        if (tempreal < this%controls%itol_rel_update(idof)) then
          idof_cnvgd_due_to_rel_update(idof,local_id) = PETSC_TRUE
        else
          if (dC_p(index) < 1.d-18) then
            idof_cnvgd_due_to_rel_update(idof,local_id) = PETSC_TRUE
          endif
        endif
        species_max_relative_change(idof) = &
          max(species_max_relative_change(idof),tempreal)
      !-----------------------------------------------------------------
        idof_cnvgd_due_to_scaled_res(idof,local_id) = PETSC_FALSE
        tempreal = dabs(r_p(index)/accum_p(index))
        if (tempreal < this%controls%itol_scaled_res(idof)) then
          idof_cnvgd_due_to_scaled_res(idof,local_id) = PETSC_TRUE
        else
          if (r_p(index) < this%controls%itol_abs_res(idof)) then
            idof_cnvgd_due_to_scaled_res(idof,local_id) = PETSC_TRUE
          endif
        endif
        species_max_scaled_residual(idof) = &
          max(species_max_scaled_residual(idof),tempreal)
      !-----------------------------------------------------------------
        tempreal = dabs(dC_p(index))
        species_max_absolute_change(idof) = &
          max(species_max_absolute_change(idof),tempreal)
      !-----------------------------------------------------------------
        idof_cnvgd_due_to_abs_res(idof,local_id) = PETSC_FALSE
        tempreal = dabs(r_p(index))
        if (tempreal < this%controls%itol_abs_res(idof)) then
          idof_cnvgd_due_to_abs_res(idof,local_id) = PETSC_TRUE
        endif
        species_max_absolute_residual(idof) = &
          max(species_max_absolute_residual(idof),tempreal)
      !-----------------------------------------------------------------
      enddo
    enddo
    call VecRestoreArrayReadF90(dX,dC_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(X0,C0_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%tran_r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%tran_accum,accum_p,ierr);CHKERRQ(ierr)

    do idof = 1, option%ntrandof
      if (all(idof_cnvgd_due_to_rel_update(idof,:))) then
        idof_cnvgd_due_to_update(idof) = PETSC_TRUE
      endif
      if (all(idof_cnvgd_due_to_scaled_res(idof,:)) .or. &
          all(idof_cnvgd_due_to_abs_res(idof,:))) then
        idof_cnvgd_due_to_residual(idof) = PETSC_TRUE
      endif
    enddo

    if (all(idof_cnvgd_due_to_residual) .and. all(idof_cnvgd_due_to_update)) then
      converged_flag = 1
    else
      converged_flag = 0
    endif

  endif

  WRITE(*,*)  ' --------------------------------------------------------------'
  WRITE(*,*)  '          max scaled residual = ', max_scaled_residual
  WRITE(*,*)  '                     location = ', loc_max_scaled_residual
  WRITE(*,*)  '               residual @ max = ', residual_at_max
  WRITE(*,*)  '                  accum @ max = ', accum_at_max
  WRITE(*,*)  ' --------------------------------------------------------------'
  WRITE(*,*)  '     max absolute residual = ', max_absolute_residual
  WRITE(*,*)  '                  location = ', loc_max_abs_residual
  WRITE(*,*)  ' --------------------------------------------------------------'
  WRITE(*,*)  ' idof_cnvgd_due_to_residual = ', idof_cnvgd_due_to_residual
  WRITE(*,*)  ' --------------------------------------------------------------'
  WRITE(*,*)  '      max relative update = ', max_relative_change
  WRITE(*,*)  '                 location = ', loc_max_rel_update
  WRITE(*,*)  '             update @ max = ', update_at_max
  WRITE(*,*)  '               soln @ max = ', soln_at_max
  WRITE(*,*)  '                   min C0 = ', min_C0
  WRITE(*,*)  ' --------------------------------------------------------------'
  WRITE(*,*)  ' idof_cnvgd_due_to_update = ', idof_cnvgd_due_to_update
  WRITE(*,*)  ' --------------------------------------------------------------'
  WRITE(*,*)  ' --------------------------------------------------------------'
  WRITE(*,*)  ' ITOL converged_flag = ', converged_flag
  WRITE(*,*)  ' --------------------------------------------------------------'

  
  ! get global minimum
  call MPI_Allreduce(converged_flag,temp_int,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_MIN,this%realization%option%mycomm,ierr)
  
  ! this will override all previous convergence criteria to keep iterating
    if (temp_int /= 1) then  ! means ITOL_* tolerances were not satisfied:
      this%realization%option%converged = PETSC_FALSE
      !!this%realization%option%convergence = CONVERGENCE_CUT_TIMESTEP
      this%realization%option%convergence = CONVERGENCE_KEEP_ITERATING
      if (newton_iter_number >= this%controls%max_newton_iterations) then
        this%realization%option%convergence = CONVERGENCE_CUT_TIMESTEP
      endif
    else  ! means ITOL_* tolerances were satisfied, but the previous
          ! criteria were maybe not met
      ! do nothing - let the instruction proceed based on previous criteria

      ! test:
      this%realization%option%converged = PETSC_TRUE
      this%realization%option%convergence = CONVERGENCE_CONVERGED
    endif

  if (this%print_ekg) then
    call VecGetArrayReadF90(dX,dC_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(X0,C0_p,ierr);CHKERRQ(ierr)
    max_relative_change_by_dof = -MAX_DOUBLE
    do local_id = 1, grid%nlmax
      offset = (local_id-1)*option%ntrandof
      do idof = 1, option%ntrandof
        index = idof + offset
        tempreal = dabs(dC_p(index)/C0_p(index))
        max_relative_change_by_dof(idof) = &
          max(max_relative_change_by_dof(idof),tempreal)
      enddo
    enddo
    call VecRestoreArrayReadF90(dX,dC_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(X0,C0_p,ierr);CHKERRQ(ierr)
    mpi_int = option%ntrandof
    call MPI_Allreduce(MPI_IN_PLACE,max_relative_change_by_dof,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,this%option%mycomm,ierr)
    if (OptionPrintToFile(option)) then
100 format("NUCLEAR WASTE TRANSPORT  NEWTON_ITERATION ",30es16.8)
      write(IUNIT_EKG,100) max_relative_change_by_dof(:)
    endif    
  endif

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

  implicit none
  
  class(pm_nwt_type) :: this
  PetscViewer :: viewer

  !TODO(jenn)
  print *, 'PMNWTCheckpointBinary not yet implemented.'
  stop

end subroutine PMNWTCheckpointBinary

! ************************************************************************** !

subroutine PMNWTCheckpointHDF5(this,pm_grp_id)
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  ! 

  use hdf5
  
  implicit none
  
  class(pm_nwt_type) :: this
  integer(HID_T) :: pm_grp_id

  !TODO(jenn)
  print *, 'PMNWTCheckpointHDF5 not yet implemented.'
  stop

end subroutine PMNWTCheckpointHDF5

! ************************************************************************** !

subroutine PMNWTRestartBinary(this,viewer)
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  ! 

  implicit none
  
  class(pm_nwt_type) :: this
  PetscViewer :: viewer

  !TODO(jenn)
  print *, 'PMNWTRestartBinary not yet implemented.'
  stop

end subroutine PMNWTRestartBinary

! ************************************************************************** !

subroutine PMNWTRestartHDF5(this,pm_grp_id)
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  ! 
  
  use hdf5

  implicit none
  
  class(pm_nwt_type) :: this
  integer(HID_T) :: pm_grp_id

  !TODO(jenn)
  print *, 'PMNWTRestartHDF5 not yet implemented.'
  stop

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

  character(len=MAXWORDLENGTH) :: word
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

  call DeallocateArray(this%controls%max_concentration_change)
  call DeallocateArray(this%controls%max_volfrac_change)
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
