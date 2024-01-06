module Factory_Subsurface_Read_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Simulation_Subsurface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: &
            FactorySubsurfReadFlowPM, &
            FactorySubsurfReadTransportPM, &
            FactorySubsurfReadWasteFormPM, &
            FactorySubsurfReadUFDDecayPM, &
            FactorySubsurfReadUFDBiospherePM, &
            FactorySubsurfReadWellPM, &
            FactorySubsurfReadMTPM, &
            FactorySubsurfReadGeophysicsPM, &
            FactorySubsurfReadFracturePM, &
            FactorySubsurfReadRequiredCards, &
            FactorySubsurfReadInput

contains

! ************************************************************************** !

subroutine FactorySubsurfReadFlowPM(input,option,pm)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_General_class
  use PM_Hydrate_class
  use PM_WIPP_Flow_class
  use PM_Mphase_class
  use PM_Richards_class
  use PM_TH_class
  use PM_Richards_TS_class
  use PM_TH_TS_class
  use PM_ZFlow_class
  use PM_PNF_class
  use Init_Common_module
  use General_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,SUBSURFACE_FLOW'

  nullify(pm)
  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('MODE')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'mode',error_string)
        call StringToUpper(word)
        select case(word)
          case('GENERAL','HYDRATE','WIPP_FLOW')
          ! In OptionFlowInitRealization(), numerical_derivatives is set to
          ! PETSC_FALSE, but the default for GENERAL needs to be PETSC_TRUE.
          ! This is will eventually affect all flow modes with numerical
          ! derivatives as default if analytical derivatives are available
          ! and we are keying off this flag.
          option%flow%numerical_derivatives = PETSC_TRUE
        end select
        select case(word)
          case('GENERAL')
            pm => PMGeneralCreate()
          case('HYDRATE')
            pm => PMHydrateCreate()
          case('WIPP_FLOW')
            pm => PMWIPPFloCreate()
          case('BRAGFLO')
            option%io_buffer = 'BRAGFLO mode has been merged with WIPP_FLOW. &
              &Please use WIPP_FLOW instead.'
            call PrintErrMsg(option)
          case('MPHASE')
            pm => PMMphaseCreate()
          case('RICHARDS')
            pm => PMRichardsCreate()
          case('TH')
            pm => PMTHCreate()
          case ('RICHARDS_TS')
            pm => PMRichardsTSCreate()
          case ('TH_TS')
            pm => PMTHTSCreate()
          case ('ZFLOW')
            pm => PMZFlowCreate()
          case ('PORE_FLOW')
            pm => PMPNFCreate()
          case default
            error_string = trim(error_string) // ',MODE'
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
        pm%option => option
      case('OPTIONS')
        if (.not.associated(pm)) then
          option%io_buffer = 'MODE keyword must be read first under ' // &
                             trim(error_string)
          call PrintErrMsg(option)
        endif
        call pm%ReadSimulationOptionsBlock(input)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (.not.associated(pm)) then
    option%io_buffer = 'A flow MODE (card) must be included in the &
      &SUBSURFACE_FLOW block in ' // trim(error_string) // '.'
    call PrintErrMsg(option)
  endif

end subroutine FactorySubsurfReadFlowPM

! ************************************************************************** !

subroutine FactorySubsurfReadTransportPM(input,option,pm)
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/19
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_NWT_class
  use PM_OSRT_class
  use PM_RT_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: print_refactor_msg

  error_string = 'SIMULATION,PROCESS_MODELS,SUBSURFACE_TRANSPORT'

  print_refactor_msg = PETSC_FALSE

  nullify(pm)
  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('MODE')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'mode',error_string)
        call StringToUpper(word)
        select case(word)
          case('GIRT')
            pm => PMRTCreate()
            option%itranmode = RT_MODE
          case('OSRT')
            pm => PMOSRTCreate()
            option%itranmode = RT_MODE
            option%transport%reactive_transport_coupling = OPERATOR_SPLIT
          case('NWT')
            pm => PMNWTCreate()
            option%itranmode = NWT_MODE
          case default
            error_string = trim(error_string) // ',MODE'
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
        pm%option => option
      case('OPTIONS')
        if (.not.associated(pm)) then
          option%io_buffer = 'MODE keyword must be read first under ' // &
                             trim(error_string)
          call PrintErrMsg(option)
        endif
        call pm%ReadSimulationOptionsBlock(input)
      case('GLOBAL_IMPLICIT')
        print_refactor_msg = PETSC_TRUE
        exit
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (print_refactor_msg .or. .not.associated(pm)) then
    if (OptionPrintToScreen(option)) then
      print *
      print *, 'SIMULATION'
      print *, '  SIMULATION_TYPE SUBSURFACE'
      print *, '  PROCESS_MODELS'
      print *, '    SUBSURFACE_TRANSPORT'
      print *, '      MODE GIRT'
      print *, '      OPTIONS'
      print *, '      /'
      print *, '    /'
      print *, '  /'
      print *, 'END'
      print *
    endif
    option%io_buffer = "PFLOTRAN's SUBSURFACE_TRANSPORT &
      &process model has been refactored to use the &
      &combination of the SUBSURFACE_TRANSPORT and 'MODE &
      &GIRT' keywords and an (optional) OPTIONS block. &
      &Please use the keywords above in reformatting the &
      &SIMULATION block."
    call PrintErrMsg(option)
  endif

  if (.not.associated(pm)) then
    option%io_buffer = 'A transport MODE (card) must be included in the &
      &SUBSURFACE_TRANSPORT block in ' // trim(error_string) // '.'
    call PrintErrMsg(option)
  endif

end subroutine FactorySubsurfReadTransportPM

! ************************************************************************** !

subroutine FactorySubsurfReadWasteFormPM(input,option,pm)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_Waste_Form_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,WASTE_FORM'

  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)

    select case(word)
      case('TYPE')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'mode',error_string)
        call StringToUpper(word)
        select case(word)
          case('GENERAL')
            pm => PMWFCreate()
          case default
            option%io_buffer = 'WASTE FORM type ' // trim(word) // &
              ' not recognized. Only TYPE GENERAL currently supported. &
              & TYPE GLASS or TYPE FMDM no longer supported.'
            call PrintErrMsg(option)
        end select
        pm%option => option
      case('SKIP_RESTART')
        if (.not.associated(pm)) then
          option%io_buffer = 'TYPE keyword must be read first under ' // &
                             trim(error_string)
          call PrintErrMsg(option)
        endif
        pm%skip_restart = PETSC_TRUE
      case('STEADY_STATE')
        if (.not.associated(pm)) then
          option%io_buffer = 'TYPE keyword must be read first under ' // &
                             trim(error_string)
          call PrintErrMsg(option)
        endif
        pm%steady_state = PETSC_TRUE
      case('OPTIONS')
        if (.not.associated(pm)) then
          option%io_buffer = 'TYPE keyword must be read first under ' // &
                             trim(error_string)
          call PrintErrMsg(option)
        endif
        call pm%ReadSimulationOptionsBlock(input)
      case default
        option%io_buffer = 'Keyword ' // trim(word) // &
              ' not recognized for the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (.not.associated(pm)) then
    option%io_buffer = 'TYPE card missing in ' // trim(error_string)
    call PrintErrMsg(option)
  endif

  pm%option => option

end subroutine FactorySubsurfReadWasteFormPM

! ************************************************************************** !

subroutine FactorySubsurfReadUFDDecayPM(input,option,pm)
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_UFD_Decay_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'SIMULATION,PROCESS_MODELS,UFD_DECAY'

  pm => PMUFDDecayCreate()
  pm%option => option

  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)

    found = PETSC_FALSE
    call PMBaseReadSimOptionsSelectCase(pm,input,word,found, &
                                        error_string,option)
    if (found) cycle

    select case(word)
      case default
        option%io_buffer = 'Keyword ' // trim(word) // &
              ' not recognized for the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine FactorySubsurfReadUFDDecayPM

! ************************************************************************** !

subroutine FactorySubsurfReadUFDBiospherePM(input,option,pm)
  !
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_UFD_Biosphere_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,UFD_BIOSPHERE'

  pm => PMUFDBCreate()
  pm%option => option

  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case default
        option%io_buffer = 'Keyword ' // trim(word) // &
              ' not recognized for the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine FactorySubsurfReadUFDBiospherePM

! ************************************************************************** !

subroutine FactorySubsurfReadGeophysicsPM(input,option,pm)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/25/21
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_ERT_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,SUBSURFACE_GEOPHYSICS'

  nullify(pm)
  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('MODE')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'mode',error_string)
        call StringToUpper(word)
        select case(word)
          case('ERT')
            pm => PMERTCreate()
            option%igeopmode = ERT_MODE
          case default
            option%io_buffer = 'MODE ' // trim(word) // &
              ' not recognized. Only MODE ERT currently supported for &
              & SUBSURFACE_GEOPHYSICS process models.'
            call PrintErrMsg(option)
        end select
        pm%option => option
      case('OPTIONS')
        if (.not.associated(pm)) then
          option%io_buffer = 'MODE keyword must be read first under ' // &
                             trim(error_string)
          call PrintErrMsg(option)
        endif
        call pm%ReadSimulationOptionsBlock(input)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (.not.associated(pm)) then
    option%io_buffer = 'A geophysics MODE (card) must be included in the &
      &SUBSURFACE_GEOPHYSICS block in ' // trim(error_string) // '.'
    call PrintErrMsg(option)
  endif

end subroutine FactorySubsurfReadGeophysicsPM

! ************************************************************************** !

subroutine FactorySubsurfReadWellPM(input,option,pm)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_Well_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,WELL_MODEL'

  pm => PMWellCreate()
  pm%option => option

  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case default
        option%io_buffer = 'Keyword ' // trim(word) // &
              ' not recognized for the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine FactorySubsurfReadWellPM

! ************************************************************************** !

subroutine FactorySubsurfReadFracturePM(input,option,pm)
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_Fracture_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,FRACTURE_MODEL'

  pm => PMFracCreate()
  pm%option => option

  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case default
        option%io_buffer = 'Keyword ' // trim(word) // &
              ' not recognized for the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine FactorySubsurfReadFracturePM

! ************************************************************************** !

subroutine FactorySubsurfReadMTPM(input, option, pm)
  !
  ! Author: Alex Salazar III
  ! Date: 01/19/2022
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_Material_Transform_class

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'SIMULATION,PROCESS_MODELS,MATERIAL_TRANSFORM'

  pm => PMMaterialTransformCreate()
  pm%option => option

  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)

    found = PETSC_FALSE
    call PMBaseReadSimOptionsSelectCase(pm,input,word,found, &
                                        error_string,option)
    if (found) cycle

    select case(word)
      case default
        option%io_buffer = 'Keyword ' // trim(word) // &
              ' not recognized for the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine FactorySubsurfReadMTPM

! ************************************************************************** !

subroutine FactorySubsurfReadRequiredCards(simulation,input)
  !
  ! Reads required cards from input file
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07, refactored 08/20/14, refactored 12/10/14
  !

  use Option_module
  use Discretization_module
  use Grid_module
  use Input_Aux_module
  use String_module
  use Patch_module
  use Realization_Subsurface_class
  use HDF5_Aux_module

  use Simulation_Subsurface_class
  use General_module
  use Reaction_module
  use Reaction_Aux_module
  use NW_Transport_Aux_module
  use Init_Common_module

  implicit none

  class(simulation_subsurface_type) :: simulation

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  class(realization_subsurface_type), pointer :: realization
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  PetscBool :: found
  PetscBool :: qerr

  character(len = MAXSTRINGLENGTH) :: wname

  realization => simulation%realization
  patch => realization%patch
  option => realization%option
  discretization => realization%discretization

  qerr  = PETSC_FALSE
  wname = '<missing>'
  found = PETSC_FALSE

  call InputPushBlock(input,'SUBSURFACE',option)

  ! GRID information - GRID is a required card for every simulation
  string = "GRID"
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call InputPushBlock(input,'GRID',option)
  call DiscretizationReadRequiredCards(discretization,input,option)

  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      patch => PatchCreate()
      patch%grid => discretization%grid
      realization%patch => patch
  end select
  call InputPopBlock(input,option)

  ! optional required cards - yes, an oxymoron, but we need to know if
  ! these exist before we can go any further.
  call InputRewind(input)
  call InputPushBlock(input,'REQUIRED_CARDS',option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    ! do not use InputReadCard here as this is a search operation
    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    select case(trim(card))

!....................
      case('DBASE_FILENAME')
        call InputPushCard(input,card,option)
        call InputReadFilename(input,option,string)
        call InputErrorMsg(input,option,'filename','DBASE_FILENAME')
        if (index(string,'.h5') > 0) then
          call HDF5ReadDbase(string,option)
        else
          call InputReadASCIIDbase(string,option)
        endif

!....................
      case('HDF5_WRITE_GROUP_SIZE')
        call InputPushCard(input,card,option)
        call InputReadInt(input,option,option%hdf5_write_group_size)
        call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')
        call InputSkipToEnd(input,option,'HDF5_WRITE_GROUP_SIZE')

      case('HDF5_READ_GROUP_SIZE')
        call InputPushCard(input,card,option)
        call InputReadInt(input,option,option%hdf5_read_group_size)
        call InputErrorMsg(input,option,'HDF5_READ_GROUP_SIZE','Group size')

!....................
      case('PROC')
        call InputPushCard(input,card,option)
        ! processor decomposition
        if (realization%discretization%itype == STRUCTURED_GRID) then
          grid => realization%patch%grid
          ! strip card from front of string
          call InputReadInt(input,option,grid%structured_grid%npx)
          call InputDefaultMsg(input,option,'npx')
          call InputReadInt(input,option,grid%structured_grid%npy)
          call InputDefaultMsg(input,option,'npy')
          call InputReadInt(input,option,grid%structured_grid%npz)
          call InputDefaultMsg(input,option,'npz')

          if (OptionIsIORank(option) .and. OptionPrintToScreen(option)) then
            option%io_buffer = ' Processor Decomposition:'
            call PrintMsg(option)
            write(option%io_buffer,'("  npx   = ",3x,i4)') &
              grid%structured_grid%npx
            call PrintMsg(option)
            write(option%io_buffer,'("  npy   = ",3x,i4)') &
              grid%structured_grid%npy
            call PrintMsg(option)
            write(option%io_buffer,'("  npz   = ",3x,i4)') &
              grid%structured_grid%npz
            call PrintMsg(option)
          endif

          if (option%comm%size /= grid%structured_grid%npx * &
                                 grid%structured_grid%npy * &
                                 grid%structured_grid%npz) then
            write(option%io_buffer,*) 'Incorrect number of processors &
              &specified: ',grid%structured_grid%npx*grid%structured_grid%npy* &
              grid%structured_grid%npz,' commsize = ',option%comm%size
            call PrintErrMsg(option)
          endif
        endif

!....................
      case('CHEMISTRY')
        call InputPushCard(input,card,option)
        if (option%itranmode /= RT_MODE) then
          option%io_buffer = 'CHEMISTRY card is included, but &
            &SUBSURFACE_TRANSPORT MODE GIRT/OSRT was not specified in the &
            &SIMULATION block.'
          call PrintErrMsg(option)
        endif
        !geh: for some reason, we need this with CHEMISTRY read for
        !     multicontinuum
 !       option%use_sc = PETSC_TRUE
        call ReactionInit(realization%reaction,input,option)
        realization%reaction_base => realization%reaction

!....................
      case('NUCLEAR_WASTE_CHEMISTRY')
        call InputPushCard(input,card,option)
        if (option%itranmode /= NWT_MODE) then
          option%io_buffer = 'NUCLEAR_WASTE_CHEMISTRY card is included, but &
            &SUBSURFACE_TRANSPORT MODE NWT was not specified in the &
            &SIMULATION block.'
          call PrintErrMsg(option)
        endif
        realization%reaction_nw => NWTReactionCreate()
        realization%reaction_base => realization%reaction_nw
        call NWTRead(realization%reaction_nw,input,option)

    end select
  enddo
  call InputPopBlock(input,option) ! REQUIRED_CARDS
  call InputPopBlock(input,option) ! SUBSURFACE

end subroutine FactorySubsurfReadRequiredCards

! ************************************************************************** !

subroutine FactorySubsurfReadInput(simulation,input)
  !
  ! Reads pflow input file
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !

  use Option_module
  use Field_module
  use Grid_module
  use Grid_Unstructured_Aux_module
  use Grid_Structured_module
  use Solver_module
  use Material_module
  use Saturation_Function_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Creep_Closure_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Dataset_module
  use Dataset_Common_HDF5_class
  use Fluid_module
  use Realization_Common_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Region_module
  use Condition_module
  use Transport_Constraint_Base_module
  use Transport_Constraint_RT_module
  use Transport_Constraint_NWT_module
  use Transport_Constraint_module
  use Coupler_module
  use Strata_module
  use Observation_module
  use Integral_Flux_module
  use Waypoint_module
  use Debug_module
  use Patch_module
  use Reaction_module
  use Reaction_Aux_module
  use NW_Transport_module
  use NW_Transport_Aux_module
  use Discretization_module
  use Input_Aux_module
  use String_module
  use Units_module
  use Reaction_Mineral_module
  use Regression_module
  use Output_Aux_module
  use Output_module
  use Output_Tecplot_module
  use Data_Mediator_Dataset_class
  use EOS_module
  use EOS_Water_module
  use SrcSink_Sandbox_module
  use Klinkenberg_module
  use WIPP_module
  use Utility_module
  use Checkpoint_module
  use Simulation_Subsurface_class
  use PMC_Base_class
  use PMC_Subsurface_class
  use PMC_Geophysics_class
  use PMC_Subsurface_OSRT_class
  use PM_Base_class
  use PM_RT_class
  use PM_NWT_class
  use PM_Well_class
  use PM_Hydrate_class
  use PM_Fracture_class
  use PM_Base_class
  use Print_module
  use Timestepper_Base_class
  use Timestepper_KSP_class
  use Timestepper_SNES_class
  use Timestepper_Steady_class
  use Timestepper_TS_class
  use Time_Storage_module
  use TH_Aux_module
  use Survey_module

#ifdef SOLID_SOLUTION
  use Reaction_Solid_Solution_module, only : SolidSolutionReadFromInputFile
#endif

  implicit none

  class(simulation_subsurface_type) :: simulation

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXSTRINGLENGTH) :: string, temp_string
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: error_string

  character(len=1) :: backslash
  PetscReal :: temp_real, temp_real2
  PetscReal, pointer :: temp_real_array(:)
  PetscInt :: temp_int

  PetscBool :: vel_cent
  PetscBool :: vel_face
  PetscBool :: fluxes
  PetscBool :: mass_flowrate
  PetscBool :: energy_flowrate
  PetscBool :: aveg_mass_flowrate
  PetscBool :: aveg_energy_flowrate

  PetscInt :: flag1

  type(region_type), pointer :: region
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(geop_condition_type), pointer :: geop_condition
  class(tran_constraint_base_type), pointer :: tran_constraint
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  type(observation_type), pointer :: observation
  type(integral_flux_type), pointer :: integral_flux
  class(pmc_base_type), pointer :: master_pmc

  type(waypoint_type), pointer :: waypoint

  type(material_property_type), pointer :: material_property
  type(fluid_property_type), pointer :: fluid_property
  type(saturation_function_type), pointer :: saturation_function
  class(characteristic_curves_type), pointer :: characteristic_curves
  class(cc_thermal_type), pointer :: characteristic_curves_thermal
  class(creep_closure_type), pointer :: creep_closure

  class(realization_subsurface_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  class(reaction_rt_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  class(dataset_base_type), pointer :: dataset
  class(dataset_ascii_type), pointer :: dataset_ascii
  type(time_storage_type), pointer :: default_time_storage
  class(data_mediator_dataset_type), pointer :: flow_data_mediator
  class(data_mediator_dataset_type), pointer :: rt_data_mediator
  type(waypoint_list_type), pointer :: waypoint_list
  type(waypoint_list_type), pointer :: waypoint_list_time_card
  type(input_type), pointer :: input
  type(survey_type), pointer :: survey

  PetscReal :: dt_init
  PetscReal :: dt_min
  PetscReal :: units_conversion

  class(timestepper_base_type), pointer :: temp_timestepper

  PetscReal :: msfsalt, msfwatr, mlfsalt, mlfwatr

  class(pm_base_type), pointer :: pm_flow

  internal_units = 'not_assigned'
  nullify(default_time_storage)
  nullify(waypoint_list_time_card)

  realization => simulation%realization
  output_option => simulation%output_option
  waypoint_list => simulation%waypoint_list_subsurface
  patch => realization%patch

  if (associated(patch)) grid => patch%grid

  option => realization%option
  field => realization%field
  reaction => realization%reaction

  master_pmc => simulation%flow_process_model_coupler
  if (associated(simulation%tran_process_model_coupler)) then
    if (.not.associated(master_pmc)) then
      master_pmc => simulation%tran_process_model_coupler
    endif
  endif

  if (associated(simulation%geop_process_model_coupler)) then
    if (.not.associated(master_pmc)) then
      master_pmc => simulation%geop_process_model_coupler
    endif
  endif

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++

  call InputRewind(input)
  string = 'SUBSURFACE'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call InputPushBlock(input,'SUBSURFACE',option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    option%io_buffer = 'pflotran card:: ' // trim(card)
    call PrintMsg(option)

    select case(trim(card))

!....................
      case ('GRID')
        call DiscretizationRead(realization%discretization,input,option)

!....................
      case ('CHEMISTRY')
        call ReactionReadPass2(reaction,input,option)

!....................
      case('NUCLEAR_WASTE_CHEMISTRY')
        call NWTReadPass2(realization%reaction_nw,input,option)

!....................
      case ('SPECIFIED_VELOCITY')
        if (option%nflowdof > 0) then
          option%io_buffer = 'SPECIFIED_VELOCITY fields may not be used &
            &with a SUBSURFACE_FLOW mode.'
          call PrintErrMsg(option)
        endif
        internal_units = 'm/sec'
        flag1 = UNINITIALIZED_INTEGER ! uniform?
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','SPECIFIED_VELOCITY')
          call StringToUpper(word)
          select case(trim(word))
            case('UNIFORM?')
              flag1 = StringYesNoOther(input%buf)
            case('DATASET')
              if (flag1 == STRING_OTHER) then
                option%io_buffer = 'SPECIFIED_VELOCITY card "UNIFORM?" &
                  &must be answered with "YES"/"NO" before velocity data &
                  &can can be read.'
                call PrintErrMsg(option)
              endif
              if (flag1 == STRING_YES) then
                error_string = 'SPECIFIED_VELOCITY,UNIFORM,DATASET'
                dataset_ascii => DatasetAsciiCreate()
                dataset_ascii%data_type = DATASET_REAL
                dataset_ascii%array_width = 3 * &
                  max(option%nphase,option%transport%nphase)
                realization%uniform_velocity_dataset => dataset_ascii

                string = input%buf
                call InputReadDouble(input,option,temp_real)
                if (.not.InputError(input)) then
                  input%buf = string
                  call DatasetAsciiReadSingle(dataset_ascii,input, &
                                              temp_string,internal_units, &
                                              error_string,option)
                else
                  input%buf = string
                  input%ierr = 0
                  call InputReadCard(input,option,word)
                  call InputErrorMsg(input,option,'keyword',error_string)
                  call StringToUpper(word)
                  select case(word)
                    case('FILE')
                      error_string = trim(error_string) // ',FILE'
                      call InputReadNChars(input,option,string, &
                                           MAXSTRINGLENGTH,PETSC_TRUE)
                      call InputErrorMsg(input,option,'filename',error_string)
                      call DatasetAsciiReadFile(dataset_ascii,string, &
                                                temp_string,internal_units, &
                                                error_string,option)
                    case('LIST')
                      error_string = trim(error_string) // ',LIST'
                      call DatasetAsciiReadList(dataset_ascii,input, &
                                                temp_string,internal_units, &
                                                error_string,option)
                    case default
                      call InputKeywordUnrecognized(input,word, &
                                                    error_string,option)
                  end select
                  if (dataset_ascii%time_storage%time_interpolation_method == &
                      INTERPOLATION_NULL) then
                    option%io_buffer = 'An INTERPOLATION method (LINEAR or &
                      &STEP) must be specified for: ' // trim(error_string)
                    call PrintErrMsg(option)
                  endif
                endif
                string = 'SPECIFIED_VELOCITY,UNIFORM,DATASET'
                ! have to pass in dataset_base_type
                dataset => dataset_ascii
                call DatasetVerify(dataset,default_time_storage,string,option)
              else
! Add interface for non-uniform dataset
                call InputReadFilename(input,option, &
                                       realization%nonuniform_velocity_filename)
                call InputErrorMsg(input,option,'filename', &
                                   'SPECIFIED_VELOCITY,NONUNIFORM,DATASET')
              endif
          end select
        enddo
        call InputPopBlock(input,option)

!....................
      case ('DEBUG')
        call DebugRead(realization%debug,input,option)

!....................
      case ('PROC')

!....................
      case ('REGION')
        region => RegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','REGION')
        call PrintMsg(option,region%name)
        call RegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call RegionAddToList(region,realization%region_list)
        nullify(region)

!....................
      case ('FLOW_CONDITION')
        flow_condition => FlowConditionCreate(option)
        call InputReadWord(input,option,flow_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'FLOW_CONDITION','name')
        call PrintMsg(option,flow_condition%name)
        select case(option%iflowmode)
          case(G_MODE,WF_MODE)
            call FlowConditionGeneralRead(flow_condition,input,option)
          case(H_MODE)
            call FlowConditionHydrateRead(flow_condition,input,option)
          case default
            call FlowConditionRead(flow_condition,input,option)
        end select
        call FlowConditionAddToList(flow_condition,realization%flow_conditions)
        nullify(flow_condition)

!....................
      case ('TRANSPORT_CONDITION')
        if (option%itranmode == NULL_MODE) then
          option%io_buffer = 'TRANSPORT_CONDITIONs are not supported without &
                             &a SUBSURFACE_TRANSPORT PROCESS_MODEL.'
          call PrintErrMsg(option)
        endif
        tran_condition => TranConditionCreate(option)
        call InputReadWord(input,option,tran_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'TRANSPORT_CONDITION','name')
        call PrintMsg(option,tran_condition%name)
        call TranConditionRead(tran_condition, &
                               realization%transport_constraints, &
                               realization%reaction_base,input,option)
        call TranConditionAddToList(tran_condition, &
                                    realization%transport_conditions)
        nullify(tran_condition)

!....................
      case ('GEOPHYSICS_CONDITION')
        if (option%igeopmode == NULL_MODE) then
          option%io_buffer = 'GEOPHYSICS_CONDITIONs are not supported without &
                             &a SUBSURFACE_GEOPHYSICS PROCESS_MODEL.'
          call PrintErrMsg(option)
        endif
        geop_condition => GeopConditionCreate(option)
        call InputReadWord(input,option,geop_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'GEOPHYSICS_CONDITION','name')
        call PrintMsg(option,geop_condition%name)
        call GeopConditionRead(geop_condition,input,option)
        call GeopConditionAddToList(geop_condition, &
                                    realization%geophysics_conditions)
        nullify(geop_condition)

!....................
      case('CONSTRAINT')
        select case(option%itranmode)
          case(RT_MODE)
            tran_constraint => TranConstraintRTCreate(option)
          case(NWT_MODE)
            tran_constraint => TranConstraintNWTCreate(option)
          case default
            option%io_buffer = 'CONSTRAINTs are not supported without &
                               &a SUBSURFACE_TRANSPORT PROCESS_MODEL.'
          call PrintErrMsg(option)
        end select
        call InputReadWord(input,option,tran_constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'constraint','name')
        call PrintMsg(option,tran_constraint%name)
        select type(tc=>tran_constraint)
          class is(tran_constraint_rt_type)
            call TranConstraintRTRead(tc,reaction,input,option)
          class is(tran_constraint_nwt_type)
            call TranConstraintNWTRead(tc,realization%reaction_nw,input,option)
        end select
        call TranConstraintAddToList(tran_constraint, &
                                     realization%transport_constraints)
        nullify(tran_constraint)

!....................
      case ('BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Boundary Condition name')
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization%patch,coupler)
        nullify(coupler)

!....................
      case ('INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Initial Condition name')
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization%patch,coupler)
        nullify(coupler)

!....................
      case ('SOURCE_SINK')
        coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name')
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization%patch,coupler)
        nullify(coupler)

!....................
      case ('SOURCE_SINK_SANDBOX')
        call SSSandboxRead(input,option)

!....................
      case ('FLOW_MASS_TRANSFER')
        flow_data_mediator => DataMediatorDatasetCreate()
        call InputReadWord(input,option,flow_data_mediator%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Flow Mass Transfer name')
        call DataMediatorDatasetRead(flow_data_mediator,input,option)
        call flow_data_mediator%AddToList(realization%flow_data_mediator_list)
        nullify(flow_data_mediator)

!....................
      case ('RT_MASS_TRANSFER')
        rt_data_mediator => DataMediatorDatasetCreate()
        call InputReadWord(input,option,rt_data_mediator%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'RT Mass Transfer name')
        call DataMediatorDatasetRead(rt_data_mediator,input,option)
        call rt_data_mediator%AddToList(realization%tran_data_mediator_list)
        nullify(rt_data_mediator)

!....................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,input,option)
        call RealizationAddStrata(realization%patch,strata)
        nullify(strata)

!.....................
      case ('DATASET')
        nullify(dataset)
        call DatasetRead(input,dataset,option)
        call DatasetBaseAddToList(dataset,realization%datasets)
        nullify(dataset)

!....................

      case('REFERENCE_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%flow%reference_pressure)
        call InputErrorMsg(input,option,'Reference Pressure','value')
        call InputReadAndConvertUnits(input,option%flow%reference_pressure, &
                                      'Pa','Reference Pressure',option)
!....................

      case('REFERENCE_LIQUID_DENSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option, &
                             option%flow%reference_density(option%liquid_phase))
        call InputErrorMsg(input,option,'Reference Liquid Density','value')
        call InputReadAndConvertUnits(input, &
                           option%flow%reference_density(option%liquid_phase), &
                              'kg/m^3','Reference Density',option)
!....................

      case('REFERENCE_GAS_DENSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option, &
                             option%flow%reference_density(option%gas_phase))
        call InputErrorMsg(input,option,'Reference Gas Density','value')
        call InputReadAndConvertUnits(input, &
                              option%flow%reference_density(option%gas_phase), &
                              'kg/m^3','Reference Density',option)
!....................

      case('MINIMUM_HYDROSTATIC_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option, &
                             option%flow%minimum_hydrostatic_pressure)
        call InputErrorMsg(input,option,'Minimum Hydrostatic Pressure','value')
        call InputReadAndConvertUnits(input, &
                                    option%flow%minimum_hydrostatic_pressure, &
                                    'Pa','Minimum Hydrostatic Pressure',option)
!......................

      case('REFERENCE_TEMPERATURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%flow%reference_temperature)
        call InputErrorMsg(input,option,'Reference Temperature','value')

!......................

      case('REFERENCE_POROSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%flow%reference_porosity)
        call InputErrorMsg(input,option,'Reference Porosity','value')

!......................

      case('REFERENCE_SATURATION')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%flow%reference_saturation)
        call InputErrorMsg(input,option,'Reference Saturation','value')

!......................

      case('NONISOTHERMAL')
        option%use_isothermal = PETSC_FALSE

!......................

      case('ISOTHERMAL')
        option%use_isothermal = PETSC_TRUE

!......................

      case('UPDATE_FLOW_PERMEABILITY')
        option%flow%update_flow_perm = PETSC_TRUE

!......................

      case('DFN')
        grid%unstructured_grid%grid_type = TWO_DIM_GRID

!......................

      case("MULTIPLE_CONTINUUM")
        option%io_buffer = 'MULTIPLE_CONTINUUM must be entered under the &
          &SUBSURFACE_TRANSPORT block within the SIMULATION block.'
        call PrintErrMsg(option)

!......................

      case('SECONDARY_CONTINUUM_SOLVER')
        if (.not.option%use_sc) then
          option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can only be used &
                             &with MULTIPLE_CONTINUUM keyword.'
          call PrintErrMsg(option)
        endif
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('KEARST')
            option%secondary_continuum_solver = 1
          case('HINDMARSH')
            option%secondary_continuum_solver = 2
          case('THOMAS')
            option%secondary_continuum_solver = 3
          case default
            option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only &
                               &HINDMARSH or KEARST. For single component &
                               &chemistry THOMAS can be used.'
          call PrintErrMsg(option)
        end select
!....................

      case('SECONDARY_CONSTRAINT')
        option%io_buffer = 'SECONDARY_CONSTRAINT has been moved to the &
                            &TRANSPORT_CONDITION block. If a SECONDARY_CONSTRAINT &
                            &is not specified the secondary initial conditions &
                            &are copied from the primary.'
        call PrintErrMsg(option)


!......................

      case('BRIN','BRINE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%m_nacl)
        call InputDefaultMsg(input,option,'NaCl Concentration')

        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word(1:len_trim(word)))
          case('MOLAL')
          case('MASS')
            msfsalt = option%m_nacl
            msfwatr = 1.0 - msfsalt
            ! Convert mass salt (kg)/mass water(kg) to Kg-moles = kmol
            ! by dividing mass of salt by molecular weight salt
            ! Then factor of 1000 to convert Kg-mol to mol
            option%m_nacl = (1000.0*msfsalt/FMWNACL)/msfwatr
          case('MOLE')
            ! Convert Kg-mole salt (kg)/Kg mole water(kg) to Kg-mol/Kg water
            ! by multiplying mass of water by molecular weight water
            ! Then factor of 1000 to convert Kg-mol to mol
            mlfsalt = option%m_nacl
            mlfwatr = 1.0 - mlfsalt
            option%m_nacl = 1000.0*mlfsalt/(mlfwatr*FMWH2O)
          case default
            print *, 'Wrong unit: ', word(1:len_trim(word))
            stop
         end select

         ! Report to screen with units
         if (OptionPrintToScreen(option)) then
           print *,'Molality is ',option%m_nacl,' mol/Kg'
         endif

         ! Saturated molality is ~6.16, so above 10 is suspicious
         ! May be units error, so warn user
         if (option%m_nacl > 10.0) then
           option%io_buffer = &
           'More that 10 mols/Kg ~ 584 gms/Kg '// &
           'is an unusually high brine molality'
           call PrintWrnMsg(option)
         endif

!......................

      case ('TIMESTEPPER','NEWTON_SOLVER','LINEAR_SOLVER')
        option%io_buffer = 'TIMESTEPPER, NEWTON_SOLVER and LINEAR_SOLVER &
          &have been moved inside a NUMERICAL_METHODS block in the input &
          &file. Please see "Numerical Methods Refactor" under &
          &Announcements in the online documentation.'
        call PrintErrMsg(option)

!......................

      case ('NUMERICAL_JACOBIAN_MULTI_COUPLE')
        option%numerical_derivatives_multi_coupling = PETSC_TRUE

!......................

      case ('COMPUTE_STATISTICS')
        option%compute_statistics = PETSC_TRUE

!....................

      case ('CO2_DATABASE')
        call InputReadFilename(input,option,option%co2_database_filename)
        call InputErrorMsg(input,option,'CO2_DATABASE','filename')

!....................

      case ('NUMERICAL_METHODS')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            if (associated(simulation%flow_process_model_coupler)) then
              call simulation%flow_process_model_coupler% &
                     ReadNumericalMethods(input,word)
            else
              option%io_buffer = 'A SUBSURFACE_FLOW process model must &
                &be defined to read NUMERICAL_METHODS for FLOW.'
              call PrintErrMsg(option)
            endif
          case('TRAN','TRANSPORT')
            if (associated(simulation%tran_process_model_coupler)) then
              call simulation%tran_process_model_coupler% &
                     ReadNumericalMethods(input,word)
            else
              option%io_buffer = 'A SUBSURFACE_TRANSPORT process model must &
                &be defined to read NUMERICAL_METHODS for TRANSPORT.'
              call PrintErrMsg(option)
            endif
          case('GEOPHYSICS','GEOP')
            if (associated(simulation%geop_process_model_coupler)) then
              call simulation%geop_process_model_coupler% &
                     ReadNumericalMethods(input,word)
            else
              option%io_buffer = 'A SUBSURFACE_GEOPHYSICS process model must &
                &be defined to read NUMERICAL_METHODS for GEOPHYSICS.'
              call PrintErrMsg(option)
            endif
          case default
            option%io_buffer = 'NUMERICAL_METHODS must specify FLOW or &
                               &TRANSPORT.'
            call PrintErrMsg(option)
        end select

!....................

      case ('FLUID_PROPERTY')

        fluid_property => FluidPropertyCreate()
        call FluidPropertyRead(fluid_property,input,option)
        call FluidPropertyAddToList(fluid_property,realization%fluid_properties)
        nullify(fluid_property)

!....................

      case ('EOS')
        call EOSRead(input,option)

!....................

      case ('SATURATION_FUNCTION')
        flag1 = 0
        select case(option%iflowmode)
          case(TH_MODE)
            if (.not.option%flow%th_freezing) flag1 = 1
          case(MPH_MODE)
          case default
            flag1 = 1
        end select
        if (flag1 == 1) then
          option%io_buffer = &
            'Must compile with legacy_saturation_function=1 to use the &
            &SATURATION_FUNCTION keyword.  Otherwise, use &
            &CHARACTERISTIC_CURVES.'
          call PrintErrMsg(option)
        endif
        saturation_function => SaturationFunctionCreate(option)
        call InputReadWord(input,option,saturation_function%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','SATURATION_FUNCTION')
        call SaturationFunctionRead(saturation_function,input,option)
        call SatFunctionComputePolynomial(option,saturation_function)
        call PermFunctionComputePolynomial(option,saturation_function)
        call SaturationFunctionAddToList(saturation_function, &
                                         realization%saturation_functions)
        nullify(saturation_function)

!....................

      case ('CHARACTERISTIC_CURVES')
        select case(option%iflowmode)
          case(TH_MODE)
            if (option%flow%th_freezing) then
              option%io_buffer = 'CHARACTERISTIC_CURVES not supported in &
                    &flow mode TH with freezing. Use SATURATION_FUNCTION.'
              call PrintErrMsg(option)
            endif
          case(MPH_MODE)
              option%io_buffer = 'CHARACTERISTIC_CURVES not supported in &
                    &flow mode MPH. Use SATURATION_FUNCTION.'
              call PrintErrMsg(option)
          case(PNF_MODE)
            option%io_buffer = 'Variably-saturated flow not supported &
              &in PORE mode.'
            call PrintErrMsg(option)
        end select
        characteristic_curves => CharacteristicCurvesCreate()
        call InputReadWord(input,option,characteristic_curves%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','CHARACTERISTIC_CURVES')
        option%io_buffer = '  Name :: ' // &
          trim(characteristic_curves%name)
        call PrintMsg(option)
        call CharacteristicCurvesRead(characteristic_curves,input,option)
!        call SatFunctionComputePolynomial(option,saturation_function)
!        call PermFunctionComputePolynomial(option,saturation_function)
        call CharacteristicCurvesAddToList(characteristic_curves, &
                                          realization%characteristic_curves)
        nullify(characteristic_curves)

!....................

      case ('THERMAL_CHARACTERISTIC_CURVES')
        characteristic_curves_thermal => CharCurvesThermalCreate()
        call InputReadWord(input,option, &
             characteristic_curves_thermal%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','THERMAL_CHARACTERISTIC_CURVES')
        option%io_buffer = '  Name :: ' // &
             trim(characteristic_curves_thermal%name)
        call PrintMsg(option)
        call CharCurvesThermalRead( &
             characteristic_curves_thermal,input,option)
        call CharCurvesThermalAddToList( &
             characteristic_curves_thermal, &
             realization%characteristic_curves_thermal)
        nullify(characteristic_curves_thermal)

!....................

      case('CREEP_CLOSURE_TABLE')
        wipp => WIPPGetPtr()
        option%flow%transient_porosity = PETSC_TRUE
        option%flow%creep_closure_on = PETSC_TRUE
        creep_closure => CreepClosureCreate()
        call InputReadWord(input,option,creep_closure%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','CREEP_CLOSURE_TABLE')
        option%io_buffer = '  Name :: ' // trim(creep_closure%name)
        call PrintMsg(option)
        call creep_closure%Read(input,option)
        call CreepClosureAddToList(creep_closure, &
             wipp%creep_closure_tables)
        nullify(creep_closure)

!....................

      case ('MATERIAL_PROPERTY')

        material_property => MaterialPropertyCreate(option)
        call InputReadWord(input,option,material_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')
        option%io_buffer = '  Name :: ' // trim(material_property%name)
        call PrintMsg(option)
        call MaterialPropertyRead(material_property,input,option)
        call MaterialPropertyAddToList(material_property, &
             realization%material_properties)
        nullify(material_property)

!....................

      case ('USE_TOUCH_OPTIONS')
        option%use_touch_options = PETSC_TRUE

      case ('MPI_IO')
!        call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
!                                       '-viewer_binary_mpiio')

      case ('HANDSHAKE_IO')
        call InputReadInt(input,option,option%io_handshake_buffer_size)
        call InputErrorMsg(input,option,'io_handshake_buffer_size', &
                           'HANDSHAKE_IO')

      case ('OVERWRITE_RESTART_TRANSPORT')
        option%io_buffer = 'OVERWRITE_RESTART_TRANSPORT no longer &
          &supported. Please use SKIP_RESTART in the SUBSURFACE_TRANSPORT &
          &process model options block.'
        call PrintErrMsg(option)

      case ('OVERWRITE_RESTART_FLOW_PARAMS')
        option%io_buffer = 'OVERWRITE_RESTART_FLOW_PARAMS no longer &
          &supported. Please use REVERT_PARAMETERS_ON_RESTART in &
          &the SUBSURFACE_FLOW process model options block.'
        call PrintErrMsg(option)

      case ('INITIALIZE_FLOW_FROM_FILE')
        call InputReadFilename(input,option,option%initialize_flow_filename)
        call InputErrorMsg(input,option,'filename','INITIALIZE_FLOW_FROM_FILE')

      case ('INITIALIZE_TRANSPORT_FROM_FILE')
        call InputReadFilename(input,option, &
                               option%initialize_transport_filename)
        call InputErrorMsg(input,option,'filename', &
                           'INITIALIZE_TRANSPORT_FROM_FILE')

!....................
      case ('OBSERVATION')
        observation => ObservationCreate()
        call ObservationRead(observation,input,option)
        call ObservationAddToList(observation, &
                                  realization%patch%observation_list)
        nullify(observation)

!....................
      case ('INTEGRAL_FLUX')
        integral_flux => IntegralFluxCreate()
        call InputReadWord(input,option,integral_flux%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Integral Flux name')
        call IntegralFluxRead(integral_flux,input,option)
        call IntegralFluxAddToList(integral_flux, &
                                   realization%patch%integral_flux_list)
        nullify(integral_flux)

!.....................
      case ('WALLCLOCK_STOP')
        option%wallclock_stop_flag = PETSC_TRUE
        call InputReadDouble(input,option,option%wallclock_stop_time)
        call InputErrorMsg(input,option,'stop time',card)
        internal_units = 'sec'
        call InputReadAndConvertUnits(input,option%wallclock_stop_time, &
                                      internal_units,card,option)
        ! append to start time
        option%wallclock_stop_time = option%comm%start_time + &
                                     option%wallclock_stop_time

!....................
      case ('OUTPUT')
        if (output_option%output_read) then
          option%io_buffer = 'Only one OUTPUT block may be included in a &
            &PFLOTRAN input deck.'
          call PrintErrMsg(option)
        endif
        output_option%output_read = PETSC_TRUE
        vel_cent = PETSC_FALSE
        vel_face = PETSC_FALSE
        fluxes = PETSC_FALSE
        mass_flowrate = PETSC_FALSE
        energy_flowrate = PETSC_FALSE
        aveg_mass_flowrate = PETSC_FALSE
        aveg_energy_flowrate = PETSC_FALSE
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','OUTPUT')
          call StringToUpper(word)
        !----------------------------------------------------------------------
        !----- NEW INPUT FORMAT: ----------------------------------------------
        !----------------------------------------------------------------------
          select case(trim(word))
            case('OBSERVATION_FILE')
              call OutputFileRead(input,realization,output_option, &
                                  waypoint_list,trim(word))
            case('SNAPSHOT_FILE')
              call OutputFileRead(input,realization,output_option, &
                                  waypoint_list,trim(word))
            case('MASS_BALANCE_FILE')
              call OutputFileRead(input,realization,output_option, &
                                  waypoint_list,trim(word))
            case('TIME_UNITS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Output Time Units','OUTPUT')
              output_option%tunit = trim(word)
              internal_units = 'sec'
              output_option%tconv = &
                UnitsConvertToInternal(word,internal_units, &
                                       'OUTPUT,TIME_UNITS',option)
            case('VARIABLES')
              call OutputVariableRead(input,option, &
                                      output_option%output_variable_list)
            case('AVERAGE_VARIABLES')
              call OutputVariableRead(input,option, &
                                      output_option%aveg_output_variable_list)
            case('UNFILTER_NON_STATE_VARIABLES')
              output_option%filter_non_state_variables = PETSC_FALSE
            case('NO_SYNCHRONIZED_OUTPUT')
              output_option%force_synchronized_output = PETSC_FALSE


        !----------------------------------------------------------------------
        !----- SUPPORT FOR OLD INPUT FORMAT: ----------------------------------
        !----------------------------------------------------------------------
            case('NO_FINAL','NO_PRINT_FINAL')
              output_option%print_final_obs = PETSC_FALSE
              output_option%print_final_snap = PETSC_FALSE
              output_option%print_final_massbal = PETSC_FALSE
            case('NO_INITIAL','NO_PRINT_INITIAL')
              output_option%print_initial_obs = PETSC_FALSE
              output_option%print_initial_snap = PETSC_FALSE
              output_option%print_initial_massbal = PETSC_FALSE
            case('PRINT_INITIAL')
              output_option%print_final_obs = PETSC_TRUE
              output_option%print_final_snap = PETSC_TRUE
              output_option%print_final_massbal = PETSC_TRUE
            case('MASS_BALANCE')
              option%compute_mass_balance_new = PETSC_TRUE
              output_option%periodic_msbl_output_ts_imod = 1
              call InputReadCard(input,option,word)
              call InputDefaultMsg(input,option, &
                                   'OUTPUT,MASS_BALANCE,DETAILED')
              if (len_trim(word) > 0) then
                call StringToUpper(word)
                select case(trim(word))
                  case('DETAILED')
                    option%mass_bal_detailed = PETSC_TRUE
                  case default
                    call InputKeywordUnrecognized(input,word, &
                           'OUTPUT,MASS_BALANCE',option)
                end select
              endif
            case('PRINT_COLUMN_IDS')
              output_option%print_column_ids = PETSC_TRUE

            case ('PRINT_PRIMAL_GRID')
              output_option%print_explicit_primal_grid = PETSC_TRUE

            !out_mesh_type defaults for primal_explicit grid is vetex_centered
            case ('EXPLICIT_GRID_PRIMAL_GRID_TYPE')
              if (associated(grid%unstructured_grid)) then
                if (associated(grid%unstructured_grid%explicit_grid)) then
                  call InputReadCard(input,option,word)
                  call InputErrorMsg(input,option,word, &
                        'EXPLICIT_GRID_PRIMAL_GRID_TYPE')
                  call PrintMsg(option,word)
                  call StringToUpper(word)
                    select case (trim(word))
                      case ('VERTEX_CENTERED')
                        grid%unstructured_grid%explicit_grid% &
                           output_mesh_type = VERTEX_CENTERED_OUTPUT_MESH
                      case ('CELL_CENTERED')
                        grid%unstructured_grid%explicit_grid% &
                           output_mesh_type = CELL_CENTERED_OUTPUT_MESH
                        call OptionSetBlocking(option,PETSC_FALSE)
                        if (OptionIsIORank(option)) then
                          if (grid%unstructured_grid% &
                                explicit_grid%num_elems /= &
                              grid%unstructured_grid% &
                                explicit_grid%num_cells_global) then
                            option%io_buffer = &
                              'EXPLICIT_GRID_PRIMAL_GRID_TYPE &
                              &if CELL_CENTERED option, the number of cells &
                              &of the grid to print and those &
                              &of the computational grid must be equal.'
                            call PrintErrMsg(option)
                          end if
                        end if
                        call OptionSetBlocking(option,PETSC_TRUE)
                        call OptionCheckNonBlockingError(option)
                      case default
                        option%io_buffer ='EXPLICIT_GRID_PRIMAL_GRID_TYPE &
                                   &only VERTEX_CENTERED and CELL_CENTERED &
                                   &are supported.'
                        call PrintErrMsg(option)
                    end select
                endif
              endif

            case ('PRINT_DUAL_GRID')
              output_option%print_explicit_dual_grid = PETSC_TRUE

            case('TIMES')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'units','OUTPUT,TIMES')
              internal_units = 'sec'
              string = 'OUTPUT,TIMES'
              units_conversion = &
                UnitsConvertToInternal(word,internal_units,string,option)
              nullify(temp_real_array)
              call UtilityReadArray(temp_real_array,NEG_ONE_INTEGER, &
                                    string,input,option)
              do temp_int = 1, size(temp_real_array)
                waypoint => WaypointCreate()
                waypoint%time = temp_real_array(temp_int)*units_conversion
                waypoint%print_snap_output = PETSC_TRUE
                call WaypointInsertInList(waypoint,waypoint_list)
              enddo
              call DeallocateArray(temp_real_array)
            case('OUTPUT_FILE')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'time increment', &
                                 'OUTPUT,OUTPUT_FILE')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  call PrintSetPrintToFileFlag(option%print_flags, &
                                               PETSC_FALSE)
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%output_file_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,OUTPUT_FILE,PERIODIC')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'OUTPUT,OUTPUT_FILE',option)
              end select
            case('SCREEN')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'time increment','OUTPUT,SCREEN')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  call PrintSetPrintToScreenFlag(option%print_flags, &
                                                 PETSC_FALSE)
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%screen_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,SCREEN')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'OUTPUT,SCREEN',option)
              end select
            case('PERIODIC')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'time increment', &
                                 'OUTPUT,PERIODIC')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  string = 'OUTPUT,PERIODIC,TIME'
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment',string)
                  internal_units = 'sec'
                  call InputReadAndConvertUnits(input,temp_real, &
                                                internal_units,string,option)
                  output_option%periodic_snap_output_time_incr = temp_real
                  call InputReadCard(input,option,word)
                  if (input%ierr == 0) then
                    if (StringCompareIgnoreCase(word,'between')) then
                      call InputReadDouble(input,option,temp_real)
                      call InputErrorMsg(input,option,'start time',string)
                      internal_units = 'sec'
                      call InputReadAndConvertUnits(input,temp_real, &
                                                internal_units, &
                                                trim(string)//',START TIME', &
                                                option)
                      call InputReadCard(input,option,word)
                      if (.not.StringCompareIgnoreCase(word,'and')) then
                        input%ierr = 1
                      endif
                      call InputErrorMsg(input,option,'AND',string)
                      call InputReadDouble(input,option,temp_real2)
                      call InputErrorMsg(input,option,'end time',string)
                      internal_units = 'sec'
                      call InputReadAndConvertUnits(input,temp_real2, &
                                                internal_units, &
                                                trim(string)//',END TIME', &
                                                option)
                      do
                        waypoint => WaypointCreate()
                        waypoint%time = temp_real
                        waypoint%print_snap_output = PETSC_TRUE
                        call WaypointInsertInList(waypoint,waypoint_list)
                        temp_real = temp_real + &
                          output_option%periodic_snap_output_time_incr
                        if (temp_real > temp_real2) exit
                      enddo
                      output_option%periodic_snap_output_time_incr = 0.d0
                    else
                      input%ierr = 1
                      call InputErrorMsg(input,option,'BETWEEN', &
                                          'OUTPUT,PERIODIC,TIME')
                    endif
                  endif
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_snap_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,TIMESTEP')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'OUTPUT,PERIODIC',option)
              end select
            case('OBSERVATION_TIMES')
              output_option%print_observation = PETSC_TRUE
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time units', &
                   'OUTPUT,OBSERVATION_TIMES')
              internal_units = 'sec'
              units_conversion = &
                UnitsConvertToInternal(word,internal_units, &
                                       'OUTPUT,OBSERVATION_TIMES,TIME_UNITS', &
                                       option)
              string = 'OUTPUT,OBSERVATION_TIMES,TIMES'
              nullify(temp_real_array)
              call UtilityReadArray(temp_real_array,NEG_ONE_INTEGER, &
                                    string,input,option)
              do temp_int = 1, size(temp_real_array)
                waypoint => WaypointCreate()
                waypoint%time = temp_real_array(temp_int)*units_conversion
                waypoint%print_obs_output = PETSC_TRUE
                call WaypointInsertInList(waypoint,waypoint_list)
                waypoint => WaypointCreate()
                waypoint%time = temp_real_array(temp_int)*units_conversion
                waypoint%print_msbl_output = PETSC_TRUE
                call WaypointInsertInList(waypoint,waypoint_list)
              enddo
              call DeallocateArray(temp_real_array)
            case('PERIODIC_OBSERVATION')
              output_option%print_observation = PETSC_TRUE
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'periodic increment type', &
                'OUTPUT,PERIODIC_OBSERVATION')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  string = 'OUTPUT,PERIODIC_OBSERVATION,TIME'
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment',string)
                  internal_units = 'sec'
                  call InputReadAndConvertUnits(input,temp_real, &
                                                internal_units,string,option)
                  output_option%periodic_obs_output_time_incr = temp_real
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_obs_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIMESTEP')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'OUTPUT,PERIODIC_OBSERVATION',option)
              end select
            case('FORMAT')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'keyword','OUTPUT,FORMAT')
              call StringToUpper(word)
              select case(trim(word))
                case ('HDF5')
                  output_option%print_hdf5 = PETSC_TRUE
                  call InputReadCard(input,option,word)
                  call InputDefaultMsg(input,option, &
                                       'OUTPUT,FORMAT,HDF5,# FILES')
                  if (len_trim(word) > 0) then
                    call StringToUpper(word)
                    select case(trim(word))
                      case('SINGLE_FILE')
                        output_option%print_single_h5_file = PETSC_TRUE
                      case('MULTIPLE_FILES')
                        output_option%print_single_h5_file = PETSC_FALSE
                        output_option%times_per_h5_file = 1
                        call InputReadCard(input,option,word)
                        if (len_trim(word)>0) then
                          select case(trim(word))
                            case('TIMES_PER_FILE')
                              call InputReadInt(input,option, &
                                              output_option%times_per_h5_file)
                              call InputErrorMsg(input,option, &
                                'timestep increment', &
                                'OUTPUT,FORMAT,HDF5,MULTIPLE_FILES,&
                                &TIMES_PER_FILE')
                            case default
                              call InputKeywordUnrecognized(input,word, &
                                    'OUTPUT,FORMAT,HDF5,MULTIPLE_FILES',option)
                          end select
                        endif
                      case default
                        call InputKeywordUnrecognized(input,word, &
                               'OUTPUT,FORMAT,HDF5',option)
                    end select
                  endif
                case ('TECPLOT')
                  output_option%print_tecplot = PETSC_TRUE
                  call InputReadCard(input,option,word)
                  call InputErrorMsg(input,option,'TECPLOT','OUTPUT,FORMAT')
                  call StringToUpper(word)
                  select case(trim(word))
                    case('POINT')
                      output_option%tecplot_format = TECPLOT_POINT_FORMAT
                    case('BLOCK')
                      output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
                    case('FEBRICK')
                      output_option%tecplot_format = TECPLOT_FEBRICK_FORMAT
                    case default
                      call InputKeywordUnrecognized(input,word, &
                               'OUTPUT,FORMAT,TECPLOT',option)
                  end select
                  if (output_option%tecplot_format == TECPLOT_POINT_FORMAT &
                      .and. option%comm%size > 1) then
                    output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
                  endif
                  if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
                    output_option%tecplot_format = TECPLOT_FEBRICK_FORMAT
                  endif
                case ('VTK')
                  output_option%print_vtk = PETSC_TRUE
                case default
                  call InputKeywordUnrecognized(input,word, &
                                                'OUTPUT,FORMAT',option)
              end select
            case('VELOCITY_AT_CENTER')
              vel_cent = PETSC_TRUE
            case('VELOCITY_AT_FACE')
              vel_face = PETSC_TRUE
            case('FLUXES')
              fluxes = PETSC_TRUE
            case('FLOWRATES','FLOWRATE')
              mass_flowrate = PETSC_TRUE
              energy_flowrate = PETSC_TRUE
            case('MASS_FLOWRATE')
              mass_flowrate = PETSC_TRUE
            case('ENERGY_FLOWRATE')
              energy_flowrate = PETSC_TRUE
            case('AVERAGE_FLOWRATES','AVERAGE_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
              aveg_energy_flowrate = PETSC_TRUE
            case('AVERAGE_MASS_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
            case('AVERAGE_ENERGY_FLOWRATE')
              aveg_energy_flowrate = PETSC_TRUE
            case ('HDF5_WRITE_GROUP_SIZE')
              call InputReadInt(input,option,option%hdf5_write_group_size)
              call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE', &
                                 'Group size')
            case('EXTEND_HDF5_TIME_FORMAT')
              output_option%extend_hdf5_time_format = PETSC_TRUE
            case ('ACKNOWLEDGE_VTK_FLAW')
              output_option%vtk_acknowledgment = PETSC_TRUE
            case default
              call InputKeywordUnrecognized(input,word,'OUTPUT',option)
          end select

        enddo
        call InputPopBlock(input,option)

  ! If VARIABLES were not specified within the *_FILE blocks, point their
  ! variable lists to the master variable list, which can be specified within
  ! the OUTPUT block. If no VARIABLES are specified for the master list, the
  ! defaults will be populated.
          if (.not.associated(output_option%output_snap_variable_list%first) &
              .and.(output_option%output_snap_variable_list%flow_vars .and. &
                    output_option%output_snap_variable_list%energy_vars)) then
            call OutputVariableListDestroy( &
                 output_option%output_snap_variable_list)
            output_option%output_snap_variable_list => &
                 output_option%output_variable_list
          endif
          if (.not.associated(output_option%output_obs_variable_list%first) &
              .and.(output_option%output_obs_variable_list%flow_vars .and. &
                    output_option%output_obs_variable_list%energy_vars)) then
            call OutputVariableListDestroy( &
                 output_option%output_obs_variable_list)
            output_option%output_obs_variable_list => &
                output_option%output_variable_list
          endif

        if (vel_cent) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_vel_cent = PETSC_TRUE
          if (output_option%print_hdf5) &
            output_option%print_hdf5_vel_cent = PETSC_TRUE
          if (output_option%print_vtk) &
            output_option%print_vtk_vel_cent = PETSC_TRUE
        endif
        if (vel_face) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_vel_face = PETSC_TRUE
          if (output_option%print_hdf5) &
           output_option%print_hdf5_vel_face = PETSC_TRUE
        endif
        if (fluxes) then
          output_option%print_fluxes = PETSC_TRUE
        endif
        if (output_option%aveg_output_variable_list%nvars>0) then
          if (Equal(output_option%periodic_snap_output_time_incr,0.d0)) then
            option%io_buffer = 'Keyword: AVERAGE_VARIABLES defined without &
                               &PERIODIC TIME being set.'
            call PrintErrMsg(option)
          endif
          if (.not.output_option%print_hdf5) then
            option%io_buffer = 'Keyword: AVERAGE_VARIABLES only defined for FORMAT HDF5'
            call PrintErrMsg(option)
          endif
        endif
        if (mass_flowrate .or. energy_flowrate .or. aveg_mass_flowrate .or. &
            aveg_energy_flowrate .or. fluxes) then
          if (output_option%print_hdf5) then
            output_option%print_hdf5_mass_flowrate = mass_flowrate
            output_option%print_hdf5_energy_flowrate = energy_flowrate
            output_option%print_hdf5_aveg_mass_flowrate = aveg_mass_flowrate
            output_option%print_hdf5_aveg_energy_flowrate = aveg_energy_flowrate
            if (aveg_mass_flowrate.or.aveg_energy_flowrate) then
              if (Equal(output_option%periodic_snap_output_time_incr,0.d0)) then
                option%io_buffer = 'Keyword: AVEGRAGE_FLOWRATES/ &
                  &AVEGRAGE_MASS_FLOWRATE/ENERGY_FLOWRATE defined without &
                  &PERIODIC TIME being set.'
                call PrintErrMsg(option)
              endif
            endif
          endif
          option%flow%store_fluxes = PETSC_TRUE
          if (associated(grid%unstructured_grid)) then
            if (associated(grid%unstructured_grid%explicit_grid)) then
              option%flow%store_fluxes = PETSC_TRUE
              output_option%print_explicit_flowrate = mass_flowrate
            endif
          endif
        endif
        if (associated(grid%unstructured_grid)) then
          if (associated(grid%unstructured_grid%explicit_grid)) then
            if (.not.output_option%print_hdf5.and.  &
                (grid%unstructured_grid%explicit_grid%output_mesh_type == &
                 CELL_CENTERED_OUTPUT_MESH)) then
                option%io_buffer = 'unstructured explicit grid &
                  &output_mesh_type = CELL_CENTERED supported for hdf5 only'
                call PrintErrMsg(option)
            end if
          end if
        end if

!.....................
      case ('REGRESSION')
        call RegressionRead(simulation%regression,input,option)

!.....................
      case ('TIME')
        dt_init = UNINITIALIZED_DOUBLE
        dt_min = UNINITIALIZED_DOUBLE
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'word','TIME')
          select case(trim(word))
            case('SCREEN_UNITS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Screen Units','TIME')
              internal_units = 'sec'
              temp_real2 = UnitsConvertToInternal(word,internal_units, &
                                                  'TIME,SCREEN_UNITS',option)
              output_option%tunit = trim(word)
              output_option%tconv = temp_real2
            case('STEADY_STATE')
              option%io_buffer = 'STEADY_STATE no longer supported under &
                &TIME card. Please enter under process model OPTIONS.'
              call PrintErrMsg(option)
            case('FINAL_TIME')
              ! cannot use InputReadAndConvertUnits here because we need to
              ! store the units if output_option%tunit is not set
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Final Time',card)
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Final Time Units',card)
              internal_units = 'sec'
              temp_real2 = UnitsConvertToInternal(word,internal_units, &
                                                  'TIME,FINAL_TIME',option)
              if (len_trim(output_option%tunit) == 0) then
                output_option%tunit = trim(word)
                output_option%tconv = temp_real2
              endif
              waypoint => WaypointCreate()
              waypoint%final = PETSC_TRUE
              waypoint%time = temp_real*temp_real2
              waypoint%print_snap_output = PETSC_TRUE
              ! do not place final time in waypoint_list_time_card
              call WaypointInsertInList(waypoint,waypoint_list)
            case('INITIAL_TIMESTEP_SIZE')
              call InputReadDouble(input,option,dt_init)
              call InputErrorMsg(input,option,'INITIAL_TIMESTEP_SIZE',card)
              internal_units = 'sec'
              call InputReadAndConvertUnits(input,dt_init,internal_units, &
                                            'TIME,INITIAL_TIMESTEP_SIZE', &
                                            option)
            case('MINIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,dt_min)
              call InputErrorMsg(input,option,'MINIMUM_TIMESTEP_SIZE',card)
              internal_units = 'sec'
              call InputReadAndConvertUnits(input,dt_min,internal_units, &
                                            'TIME,MINIMUM_TIMESTEP_SIZE', &
                                            option)
            case('MAXIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Minimum Timestep Size',card)
              internal_units = 'sec'
              call InputReadAndConvertUnits(input,temp_real,internal_units, &
                                            'TIME,MAXIMUM_TIMESTEP_SIZE', &
                                            option)
              waypoint => WaypointCreate()
              waypoint%dt_max = temp_real
              call InputReadCard(input,option,word)
              if (input%ierr == 0) then
                call StringToUpper(word)
                if (StringCompare(word,'AT',TWO_INTEGER)) then
                  call InputReadDouble(input,option,waypoint%time)
                  call InputErrorMsg(input,option,'MAXIMUM_TIMESTEP_SIZE &
                                                  &Update Time',card)
                  internal_units = 'sec'
                  call InputReadAndConvertUnits(input,waypoint%time, &
                                                internal_units, &
                                                'TIME,MAXIMUM_TIMESTEP_SIZE,&
                                                &Update Time',option)
                else
                  option%io_buffer = 'Keyword under "MAXIMUM_TIMESTEP_SIZE" &
                                     &after maximum timestep size should &
                                     &be "AT".'
                  call PrintErrMsg(option)
                endif
              else
                waypoint%time = 0.d0
              endif
              if (.not.associated(waypoint_list_time_card)) then
                waypoint_list_time_card => WaypointListCreate()
              endif
              call WaypointInsertInList(waypoint, &
                                        waypoint_list_time_card)
            case default
              call InputKeywordUnrecognized(input,word,'TIME',option)
          end select
        enddo
        call InputPopBlock(input,option)

        ! we store dt_init and dt_min in local variables so that they
        ! cannot overwrite what has previously been set in the respective
        ! timestepper object member variable
        if (Initialized(dt_init)) then
          if (Initialized(master_pmc%timestepper%dt_init)) then
            option%io_buffer = 'INITIAL_TIMESTEP_SIZE may be included &
              &under either the TIME or TIMESTEPPER ' // &
              trim(master_pmc%timestepper%name) // ' card, but not both.'
            call PrintErrMsg(option)
          endif
          if (associated(simulation%flow_process_model_coupler)) then
            temp_timestepper => &
              simulation%flow_process_model_coupler%timestepper
            if (associated(temp_timestepper)) then
              if (Uninitialized(temp_timestepper%dt_init)) then
                temp_timestepper%dt_init = dt_init
              endif
            endif
          endif
          if (associated(simulation%tran_process_model_coupler)) then
            temp_timestepper => &
              simulation%tran_process_model_coupler%timestepper
            if (associated(temp_timestepper)) then
              if (Uninitialized(temp_timestepper%dt_init)) then
                temp_timestepper%dt_init = dt_init
              endif
            endif
          endif
        endif
        if (Initialized(dt_min)) then
          if (Initialized(master_pmc%timestepper%dt_min)) then
            option%io_buffer = 'MINIMUM_TIMESTEP_SIZE may be included &
              &under either the TIME or TIMESTEPPER ' // &
              trim(master_pmc%timestepper%name) // ' card, but not both.'
            call PrintErrMsg(option)
          endif
          if (associated(simulation%flow_process_model_coupler)) then
            temp_timestepper => &
              simulation%flow_process_model_coupler%timestepper
            if (associated(temp_timestepper)) then
              if (Uninitialized(temp_timestepper%dt_min)) then
                temp_timestepper%dt_min = dt_min
              endif
            endif
          endif
          if (associated(simulation%tran_process_model_coupler)) then
            temp_timestepper => &
              simulation%tran_process_model_coupler%timestepper
            if (associated(temp_timestepper)) then
              if (Uninitialized(temp_timestepper%dt_min)) then
                temp_timestepper%dt_min = dt_min
              endif
            endif
          endif
        endif

!......................
      case ('HDF5_READ_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_read_group_size)
        call InputErrorMsg(input,option,'HDF5_READ_GROUP_SIZE','Group size')

!......................
      case ('HDF5_WRITE_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_write_group_size)
        call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')

!....................
      case('WIPP')
        wipp => WIPPGetPtr()
        call WIPPRead(input,option)

!....................
      case('KLINKENBERG_EFFECT')
        wipp => WIPPGetPtr()
        call KlinkenbergInit()
        klinkenberg => KlinkenbergCreate()
        call Klinkenberg%Read(input,option)

!....................
      case ('ONLY_VERTICAL_FLOW')
        option%flow%only_vertical_flow = PETSC_TRUE
        if (option%iflowmode /= TH_MODE .and. &
            option%iflowmode /= TH_TS_MODE .and. &
            option%iflowmode /= RICHARDS_MODE .and. &
            option%iflowmode /= RICHARDS_TS_MODE) then
          option%io_buffer = 'ONLY_VERTICAL_FLOW implemented in RICHARDS, &
                              &RICHARDS_TS and TH, TH_TS modes.'
          call PrintErrMsg(option)
        endif

!....................
      case ('QUASI_3D')
        option%flow%quasi_3d = PETSC_TRUE
        option%flow%only_vertical_flow = PETSC_TRUE
        if (option%iflowmode /= RICHARDS_MODE .and. &
            option%iflowmode /= RICHARDS_TS_MODE) then
          option%io_buffer = 'QUASI_3D implemented in RICHARDS and &
                              &RICHARDS_TS modes.'
          call PrintErrMsg(option)
        endif

!....................
      case ('ONLY_ENERGY_EQ')
        option%flow%only_energy_eq = PETSC_TRUE
        if (option%iflowmode /= TH_MODE .and. &
            option%iflowmode /= TH_TS_MODE) then
          option%io_buffer = 'ONLY_ENERGY_EQ applicable only in TH and &
                              &TH_TS modes.'
          call PrintErrMsg(option)
        endif

!....................
      case ('RELATIVE_PERMEABILITY_AVERAGE')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (trim(word))
          case ('UPWIND')
            option%flow%rel_perm_aveg = UPWIND
          case ('HARMONIC')
            option%flow%rel_perm_aveg = HARMONIC
          case ('DYNAMIC_HARMONIC')
            option%flow%rel_perm_aveg = DYNAMIC_HARMONIC
          case default
            option%io_buffer = 'Cannot identify the specificed &
              &RELATIVE_PERMEABILITY_AVERAGE.'
            call PrintErrMsg(option)
          end select

!....................
      case ('DBASE_FILENAME')

!....................
      case ('MIN_ALLOWABLE_SCALE')
        call InputReadDouble(input,option,option%min_allowable_scale)
        call InputErrorMsg(input,option,'minimium allowable scaling factor', &
                           'InitSubsurface')

!....................
      case ('HYDRATE')
        pm_flow => simulation%flow_process_model_coupler%pm_list
        select type (pm_flow)
          class is (pm_hydrate_type)
            call PMHydrateReadParameters(input,pm_flow,option)
          class default
            option%io_buffer = 'Keyword HYDRATE not recognized for the ' // &
                               trim(option%flowmode) // ' flow process model.'
            call PrintErrMsg(option)
        end select

!....................
      case ('SURVEY')
        survey => SurveyCreate()
        call SurveyRead(survey,input,option)
        realization%survey => survey
        nullify(survey)

!....................
      case ('WELLBORE_MODEL')
        call PMWellReadPass2(input,option)

!....................
      case ('GEOTHERMAL_FRACTURE_MODEL')
        call PMFracReadPass2(input,option)

!....................
      case ('END_SUBSURFACE')
        exit

      case default
        call InputKeywordUnrecognized(input,word, &
                                      'SubsurfaceReadInput()',option)
    end select

  enddo
  call InputPopBlock(input,option) ! SUBSURFACE

  call PrintInitFlags(option%print_flags,option%driver%print_flags)

  if (associated(simulation%flow_process_model_coupler)) then
    select case(option%iflowmode)
      case(MPH_MODE,G_MODE,TH_MODE,WF_MODE,RICHARDS_TS_MODE,TH_TS_MODE, &
           H_MODE)
        if (option%flow%steady_state) then
          option%io_buffer = 'Steady state solution is not supported with &
            &the current flow mode.'
          call PrintErrMsg(option)
        endif
    end select
  endif
  if (associated(simulation%tran_process_model_coupler)) then
    select case(option%itranmode)
      case(NULL_MODE)
        if (option%transport%steady_state) then
          option%io_buffer = 'Steady state solution is not supported with &
                             &the current transport mode.'
          call PrintErrMsg(option)
        endif
    end select
  endif

  ! must come after setup of timestepper steady above. otherwise, the
  ! destruction of the waypoint lists will fail with to pointer to the
  ! same list.
  if (associated(master_pmc%timestepper%local_waypoint_list) .and. &
      associated(waypoint_list_time_card)) then
    option%io_buffer = 'MAXIMUM_TIMESTEP_SIZE may be included under either &
      &the TIME or TIMESTEPPER ' // trim(master_pmc%timestepper%name) // &
      ' card, but not both.'
    call PrintErrMsg(option)
  endif
  if (associated(waypoint_list_time_card)) then
    call WaypointListMerge(simulation%waypoint_list_subsurface, &
                           waypoint_list_time_card,option)
    ! DO NOT destroy as both pointer point to the same list
    nullify(waypoint_list_time_card)
  else
    call WaypointListMerge(simulation%waypoint_list_subsurface, &
                           master_pmc%timestepper%local_waypoint_list,option)
  endif

end subroutine FactorySubsurfReadInput

end module Factory_Subsurface_Read_module
