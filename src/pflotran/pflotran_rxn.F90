!=======================================================================
! Please see LICENSE and COPYRIGHT files at top of repository.
!=======================================================================
module BatchChem

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: BatchChemInitializeReactions, &
            BatchChemProcessConstraints

contains

! ************************************************************************** !

subroutine BatchChemInitializeReactions(option, input, reaction)

  use Reaction_module
  use Reaction_Aux_module
  use Reaction_Database_module
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  class(reaction_rt_type), pointer :: reaction
  character(len=MAXSTRINGLENGTH) :: string

  ! check for a chemistry block in the  input file
  string = "CHEMISTRY"
  call InputFindStringInFile(input, option, string)
  if (.not.InputError(input)) then
    ! found a chemistry block, initialize the chemistry.

    ! NOTE(bja): ReactionInit() only does a first pass through the
    ! input file to check for a select items
    call ReactionInit(reaction, input, option)
    ! rewind the input file to prepare for the second pass
    call InputFindStringInFile(input, option, string)
    ! the second pass through the input file to read the remaining blocks
    call ReactionReadPass2(reaction, input, option)
  else
     ! TODO(bja): no chemistry block --> fatal error
  endif

  if (associated(reaction)) then
    if (reaction%use_full_geochemistry) then
       call DatabaseRead(reaction, option)
       call BasisInit(reaction, option)
    else
      ! NOTE(bja): do we need this for the batch chemistry driver?

      ! turn off activity coefficients since the database has not been read
      reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
      allocate(reaction%primary_species_print(option%ntrandof))
      reaction%primary_species_print = PETSC_TRUE
    endif
  endif

end subroutine BatchChemInitializeReactions

! ************************************************************************** !

subroutine BatchChemProcessConstraints(option, input, reaction, &
     global_auxvars, rt_auxvars, material_auxvars, transport_constraints, &
     constraint_coupler)

  use Reaction_module
  use Reaction_Aux_module
  use Reaction_Database_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use Transport_Constraint_module
  use Transport_Constraint_RT_module
  use Transport_Constraint_Base_module
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none


  type(option_type), pointer :: option
  type(input_type), pointer :: input
  class(reaction_rt_type), pointer :: reaction
  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars
  type(material_auxvar_type), pointer :: material_auxvars
  class(tran_constraint_coupler_rt_type), pointer :: constraint_coupler
  type(tran_constraint_list_type), pointer :: transport_constraints

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  class(tran_constraint_base_type), pointer :: base_constraint
  class(tran_constraint_rt_type), pointer :: tran_constraint
  class(tran_constraint_rt_type), pointer :: constraint
  PetscBool :: use_prev_soln_as_guess
  PetscInt :: num_iterations


  !
  ! read the constraints...
  !

  ! look through the input file
  call InputRewind(input)
  do
    call InputReadPflotranString(input, option)
    if (InputError(input)) exit

    call InputReadWord(input, option, word, PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    option%io_buffer = 'pflotran card:: ' // trim(card)
    call PrintMsg(option)

    select case(trim(card))
      case('CONSTRAINT')
        if (.not.associated(reaction)) then
          option%io_buffer = 'CONSTRAINTs not supported without CHEMISTRY.'
          call PrintErrMsg(option)
        endif
        tran_constraint => TranConstraintRTCreate(option)
        base_constraint => tran_constraint
        call InputReadWord(input, option, tran_constraint%name, PETSC_TRUE)
        call InputErrorMsg(input, option, 'constraint', 'name')
        call PrintMsg(option, tran_constraint%name)
        call TranConstraintRTRead(tran_constraint, reaction, input, option)
        call TranConstraintAddToList(base_constraint, transport_constraints)
        nullify(tran_constraint)

      case default
         ! do nothing
    end select
  enddo

  !
  ! process constraints
  !
  num_iterations = 0
  use_prev_soln_as_guess = PETSC_FALSE
  base_constraint => transport_constraints%first
  ! NOTE(bja): we only created one set of global and rt auxvars, so if
  ! there is more than one constratint in the input file, they will be
  ! over written.
  do
     if (.not. associated(base_constraint)) exit
     tran_constraint => TranConstraintRTCast(base_constraint)
     ! initialize constraints
     option%io_buffer = "initializing constraint : " // tran_constraint%name
     call PrintMsg(option)
     call ReactionProcessConstraint(reaction, &
                                    tran_constraint, &
                                    option)

     ! link the constraint to the constraint coupler
     constraint => TranConstraintRTCast(constraint_coupler%constraint)
     constraint_coupler%constraint_name = tran_constraint%name
     constraint%aqueous_species => tran_constraint%aqueous_species
     constraint%free_ion_guess => tran_constraint%free_ion_guess
     constraint%minerals => tran_constraint%minerals
     constraint%surface_complexes => tran_constraint%surface_complexes
     constraint%colloids => tran_constraint%colloids
     constraint_coupler%global_auxvar => global_auxvars
     constraint_coupler%rt_auxvar => rt_auxvars

     ! equilibrate
     option%io_buffer = "equilibrate constraint : " // tran_constraint%name
     call PrintMsg(option)
     call ReactionEquilibrateConstraint(rt_auxvars, global_auxvars, &
                                        material_auxvars, reaction, &
                                        tran_constraint, &
                                        num_iterations, &
                                        use_prev_soln_as_guess, &
                                        option)
     call ReactionPrintConstraint(constraint_coupler, reaction, option)
     base_constraint => base_constraint%next
  enddo

end subroutine BatchChemProcessConstraints


end module BatchChem


! ************************************************************************** !
program pflotran_rxn

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use Reaction_Database_module
  use Communicator_Aux_module
  use Option_module
  use Driver_module
  use Input_Aux_module
  use String_module

  use Transport_Constraint_module
  use Transport_Constraint_RT_module
  use Transport_Constraint_Base_module
  use PFLOTRAN_Constants_module

  use BatchChem

  implicit none


  PetscErrorCode :: ierr
  PetscBool :: option_found
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: filename_out
  class(reaction_rt_type), pointer :: reaction
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  class(driver_type), pointer :: driver

  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars
  type(material_auxvar_type), pointer :: material_auxvars

  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  type(tran_constraint_list_type), pointer :: transport_constraints
  class(tran_constraint_coupler_base_type), pointer :: constraint_coupler

  driver => DriverCreate()
  call MPI_Init(ierr);CHKERRQ(ierr)
  call CommInitPetsc(driver%comm,MPI_COMM_WORLD)

  option => OptionCreate()
  option%driver => driver
  option%fid_out = FORWARD_OUT_UNIT

  ! check for non-default input filename
  option%input_filename = "pflotran.in"
  string = '-pflotranin'
  call InputGetCommandLineString(string, option%input_filename, option_found, option)

  string = '-output_prefix'
  call InputGetCommandLineString(string, option%global_prefix, option_found, option)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)

  input => InputCreate(IN_UNIT, option%input_filename, option)

  filename_out = trim(option%global_prefix) // trim(option%group_prefix) // &
                 '.out'

  if (driver%IsIORank() .and. driver%PrintToFile()) then
    open(option%fid_out, file=filename_out, action="write", status="unknown")
  endif

  !
  ! manual initialization...
  !
  option%nphase = 1
  option%liquid_phase = 1
  option%flow%reference_density(option%liquid_phase) = 998.2

  call BatchChemInitializeReactions(option, input, reaction)

  !
  ! create the storage containers
  !
  ! NOTE(bja) : batch chem --> one cell

  ! global_auxvars --> cell by cell temperature, pressure, saturation, density
  allocate(global_auxvars)
  call GlobalAuxVarInit(global_auxvars, option)

  ! rt_auxvars --> cell by cell chemistry data
  allocate(rt_auxvars)
  call RTAuxVarInit(rt_auxvars, reaction, option)

  ! material_auxvars --> cell by cell material property data
  allocate(material_auxvars)
  call MaterialAuxVarInit(material_auxvars, option)
  material_auxvars%porosity = option%flow%reference_porosity

  ! assign default state values
  global_auxvars%pres = option%flow%reference_pressure
  global_auxvars%temp = option%flow%reference_temperature
  ! global_auxvars%den_kg = option%flow%reference_water_density
  ! NOTE(bja): option%ref_density = 0.0, so we set it manually. This is a Bad Thing(TM)
  global_auxvars%den_kg = 998.2
  global_auxvars%sat = option%flow%reference_saturation

  ! create the constraint list
  allocate(transport_constraints)
  call TranConstraintInitList(transport_constraints)
  allocate(constraint_coupler)
  constraint_coupler => TranConstraintCouplerRTCreate(option)

  call BatchChemProcessConstraints(option, input, reaction, &
                                   global_auxvars, rt_auxvars, &
                                   material_auxvars, transport_constraints, &
                               TranConstraintCouplerRTCast(constraint_coupler))

  ! cleanup
  call TranConstraintCouplerDestroy(constraint_coupler)
  call TranConstraintListDestroy(transport_constraints)
  ! FIXME(bja) : causes error freeing memory.
  !call RTAuxVarDestroy(rt_auxvars)
  !call GlobalAuxVarDestroy(global_auxvars)
  call ReactionDestroy(reaction,option)
  call GlobalAuxVarStrip(global_auxvars)
  deallocate(global_auxvars)
  nullify(global_auxvars)
  call RTAuxVarStrip(rt_auxvars)
  deallocate(rt_auxvars)
  nullify(rt_auxvars)
  call MaterialAuxVarStrip(material_auxvars)
  deallocate(material_auxvars)
  nullify(material_auxvars)
  call InputDestroy(input)
  call OptionDestroy(option)
  call DriverDestroy(driver)
  call PetscFinalize(ierr);CHKERRQ(ierr)
  call MPI_Finalize(ierr);CHKERRQ(ierr)

end program pflotran_rxn

