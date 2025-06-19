module Reaction_Mineral_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Mineral_Aux_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscReal, parameter :: perturbation_tolerance = 1.d-5

  public :: ReactionMnrlRead, &
            ReactionMnrlReadKinetics, &
            ReactionMnrlReadNucleation, &
            ReactionMnrlReadFromDatabase, &
            ReactionMnrlReadMassActOverride, &
            ReactionMnrlProcessConstraint, &
            ReactionMnrlSetup, &
            ReactionMnrlKinetics, &
            ReactionMnrlSaturationIndex, &
            ReactionMnrlUpdateTempDepCoefs, &
            ReactionMnrlUpdateSpecSurfArea, &
            ReactionMnrlUpdateKineticState, &
            ReactionMnrlReportZeroSurfArea

contains

! ************************************************************************** !

subroutine ReactionMnrlRead(mineral,input,option)
  !
  ! Reads chemical species
  !
  ! Author: Glenn Hammond
  ! Date: 08/16/12
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  type(mineral_type) :: mineral
  type(input_type), pointer :: input
  type(option_type) :: option

  type(mineral_rxn_type), pointer :: cur_mineral, prev_mineral

  nullify(prev_mineral)
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    mineral%nmnrl = mineral%nmnrl + 1

    cur_mineral => ReactionMnrlCreateMineralRxn()
    call InputReadCard(input,option,cur_mineral%name)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,MINERALS')
    if (.not.associated(mineral%mineral_list)) then
      mineral%mineral_list => cur_mineral
      cur_mineral%id = 1
    endif
    if (associated(prev_mineral)) then
      prev_mineral%next => cur_mineral
      cur_mineral%id = prev_mineral%id + 1
    endif
    prev_mineral => cur_mineral
    nullify(cur_mineral)
  enddo
  call InputPopBlock(input,option)

end subroutine ReactionMnrlRead

! ************************************************************************** !

subroutine ReactionMnrlReadKinetics(mineral,input,option)
  !
  ! Reads mineral kinetics
  !
  ! Author: Glenn Hammond
  ! Date: 10/16/08
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use Units_module
  use Reaction_Database_Aux_module
  use Utility_module

  implicit none

  type(mineral_type) :: mineral
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name

  type(mineral_rxn_type), pointer :: cur_mineral
  type(transition_state_rxn_type), pointer :: tstrxn, cur_tstrxn
  type(transition_state_prefactor_type), pointer :: prefactor, &
                                                    cur_prefactor
  type(ts_prefactor_species_type), pointer :: prefactor_species, &
                                              cur_prefactor_species
  PetscBool :: found
  PetscInt :: imnrl,icount
  PetscReal :: rate_constant

  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = -1*abs(cur_mineral%id)
    cur_mineral => cur_mineral%next
  enddo

  input%ierr = INPUT_ERROR_NONE
  icount = 0
  call InputPushBlock(input,'MINERAL_KINETICS',option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,name,PETSC_TRUE)
    call InputErrorMsg(input,option,'mineral name', &
                       'CHEMISTRY,MINERAL_KINETICS')

    cur_mineral => mineral%mineral_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_mineral)) exit
      if (StringCompare(cur_mineral%name,name,MAXWORDLENGTH)) then
        found = PETSC_TRUE
        cur_mineral%itype = MINERAL_KINETIC
        tstrxn => ReactionMnrlCreateTSTRxn()
        ! initialize to UNINITIALIZED_INTEGER to ensure that it is set
        tstrxn%precipitation_rate_constant = UNINITIALIZED_DOUBLE
        call InputPushBlock(input,name,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,error_string)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,keyword)
          error_string = 'CHEMISTRY,MINERAL_KINETICS,'//name
          call InputErrorMsg(input,option,'keyword',error_string)

          select case(trim(keyword))
            case('RATE_CONSTANT','PRECIPITATION_RATE_CONSTANT', &
                 'DISSOLUTION_RATE_CONSTANT')
              rate_constant = UNINITIALIZED_DOUBLE
              call ReactionMnrlReadRateConstant(input,cur_mineral%name, &
                                                rate_constant,error_string, &
                                                '',option)
              select case(keyword)
                case('RATE_CONSTANT')
                  tstrxn%precipitation_rate_constant = rate_constant
                  tstrxn%dissolution_rate_constant = rate_constant
                case('PRECIPITATION_RATE_CONSTANT')
                  tstrxn%precipitation_rate_constant = rate_constant
                case('DISSOLUTION_RATE_CONSTANT')
                  tstrxn%dissolution_rate_constant = rate_constant
              end select
            case('ACTIVATION_ENERGY')
!             read activation energy for Arrhenius law
              call InputReadDouble(input,option,tstrxn%activation_energy)
              call InputErrorMsg(input,option,keyword,error_string)
              call InputReadAndConvertUnits(input,tstrxn%activation_energy, &
                                            'J/mol', &
                              trim(error_string)//',ACTIVATION_ENERGY',option)
            case('AFFINITY_THRESHOLD')
!             read affinity threshold for precipitation
              call InputReadDouble(input,option,tstrxn%affinity_threshold)
              call InputErrorMsg(input,option,keyword,error_string)
            case('AFFINITY_POWER')
!             reads exponent on affinity term
              call InputReadDouble(input,option,tstrxn%affinity_factor_beta)
              call InputErrorMsg(input,option,keyword,error_string)
            case('MINERAL_SCALE_FACTOR')
!             read mineral scale factor term
              call InputReadDouble(input,option,tstrxn%mnrl_scale_factor)
              call InputErrorMsg(input,option,keyword,error_string)
            case('TEMKIN_CONSTANT')
!             reads exponent on affinity term
              call InputReadDouble(input,option,tstrxn%affinity_factor_sigma)
              call InputErrorMsg(input,option,keyword,error_string)
            case('SURFACE_AREA_POROSITY_POWER')
              call InputReadDouble(input,option,tstrxn%surf_area_porosity_pwr)
              call InputErrorMsg(input,option,keyword,error_string)
            case('SURFACE_AREA_VOL_FRAC_POWER')
              call InputReadDouble(input,option,tstrxn%surf_area_vol_frac_pwr)
              call InputErrorMsg(input,option,keyword,error_string)
            case('RATE_LIMITER')
!             read rate limiter for precipitation
              call InputReadDouble(input,option,tstrxn%rate_limiter)
              call InputErrorMsg(input,option,keyword,error_string)
            case('ARMOR_MINERAL')
                    ! read mineral name
              call InputReadWord(input,option,tstrxn%armor_min_name,PETSC_TRUE)
              call InputErrorMsg(input,option,keyword,error_string)
            case('ARMOR_PWR')
                    ! read power law exponent
              call InputReadDouble(input,option,tstrxn%armor_pwr)
              call InputErrorMsg(input,option,keyword,error_string)
            case('ARMOR_CRIT_VOL_FRAC')
              call InputReadDouble(input,option,tstrxn%armor_crit_vol_frac)
              call InputErrorMsg(input,option,keyword,error_string)
            case('SPECIFIC_SURFACE_AREA_EPSILON')
              call InputReadDouble(input,option,tstrxn%surf_area_epsilon)
              call InputErrorMsg(input,option,keyword,error_string)
            case('VOLUME_FRACTION_EPSILON')
              call InputReadDouble(input,option,tstrxn%vol_frac_epsilon)
              call InputErrorMsg(input,option,keyword,error_string)
            case('PREFACTOR')
              prefactor => ReactionMnrlCreateTSTPrefactor()
              ! Initialize to UNINITIALIZED_DOUBLE to check later whether
              ! they were set
              prefactor%activation_energy = UNINITIALIZED_DOUBLE
              call InputPushBlock(input,'PREFACTOR',option)
              do
                error_string = 'CHEMISTRY,MINERAL_KINETICS,'// &
                               trim(name)//',PREFACTOR'
                call InputReadPflotranString(input,option)
                call InputReadStringErrorMsg(input,option,error_string)
                if (InputCheckExit(input,option)) exit
                call InputReadCard(input,option,keyword)
                call InputErrorMsg(input,option,'keyword',error_string)
                select case(trim(keyword))
                  case('RATE_CONSTANT','PRECIPITATION_RATE_CONSTANT', &
                       'DISSOLUTION_RATE_CONSTANT')
                    rate_constant = UNINITIALIZED_DOUBLE
                    call ReactionMnrlReadRateConstant(input, &
                                                      cur_mineral%name, &
                                                      rate_constant, &
                                                      error_string, &
                                                      'PREFACTOR',option)
                    select case(keyword)
                      case('RATE_CONSTANT')
                        prefactor%precipitation_rate_constant = rate_constant
                        prefactor%dissolution_rate_constant = rate_constant
                      case('PRECIPITATION_RATE_CONSTANT')
                        prefactor%precipitation_rate_constant = rate_constant
                      case('DISSOLUTION_RATE_CONSTANT')
                        prefactor%dissolution_rate_constant = rate_constant
                    end select
                  case('ACTIVATION_ENERGY')
                    ! read activation energy for Arrhenius law
                    call InputReadDouble(input,option, &
                                         prefactor%activation_energy)
                    call InputErrorMsg(input,option,keyword,error_string)
                    call InputReadAndConvertUnits(input, &
                            prefactor%activation_energy,'J/mol', &
                            trim(error_string)//',ACTIVATION_ENERGY',option)
                  case('PREFACTOR_SPECIES')
                    prefactor_species => ReactionMnrlCreateTSTPrefSpec()
                    call InputReadCard(input,option,prefactor_species%name, &
                                       PETSC_TRUE)
                    call InputErrorMsg(input,option,keyword,error_string)
                    call InputPushBlock(input,'PREFACTOR_SPECIES',option)
                    do
                      error_string = 'CHEMISTRY,MINERAL_KINETICS,PREFACTOR,&
                                    &SPECIES'
                      call InputReadPflotranString(input,option)
                      call InputReadStringErrorMsg(input,option,error_string)
                      if (InputCheckExit(input,option)) exit
                      call InputReadCard(input,option,keyword)
                      call InputErrorMsg(input,option,'keyword',error_string)
                      select case(trim(keyword))
                        case('ALPHA')
                          call InputReadDouble(input,option, &
                                               prefactor_species%alpha)
                          call InputErrorMsg(input,option,keyword,error_string)
                        case('BETA')
                          call InputReadDouble(input,option, &
                                               prefactor_species%beta)
                          call InputErrorMsg(input,option,keyword,error_string)
                        case('ATTENUATION_COEF')
                          call InputReadDouble(input,option, &
                                            prefactor_species%attenuation_coef)
                          call InputErrorMsg(input,option,keyword,error_string)
                        case default
                          call InputKeywordUnrecognized(input,keyword,word, &
                                                        error_string,option)
                      end select
                    enddo
                    call InputPopBlock(input,option)
                    ! add prefactor species
                    if (.not.associated(prefactor%species)) then
                      prefactor%species => prefactor_species
                    else ! append to end of list
                      cur_prefactor_species => prefactor%species
                      do
                        if (.not.associated(cur_prefactor_species%next)) then
                          cur_prefactor_species%next => prefactor_species
                          exit
                        else
                          cur_prefactor_species => cur_prefactor_species%next
                        endif
                      enddo
                    endif
                  case default
                    call InputKeywordUnrecognized(input,keyword,error_string, &
                                                  option)
                end select
              enddo
              call InputPopBlock(input,option)
              ! add prefactor
              if (.not.associated(tstrxn%prefactor)) then
                tstrxn%prefactor => prefactor
              else ! append to end of list
                cur_prefactor => tstrxn%prefactor
                do
                  if (.not.associated(cur_prefactor%next)) then
                    cur_prefactor%next => prefactor
                    exit
                  else
                    cur_prefactor => cur_prefactor%next
                  endif
                enddo
              endif
            case('SURFACE_AREA_FUNCTION')
              mineral%update_surface_area = PETSC_TRUE
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,keyword,error_string)
              call StringToUpper(word)
              select case(word)
                case('CONSTANT')
                  tstrxn%surf_area_function = MINERAL_SURF_AREA_F_NULL
                case('POROSITY_RATIO')
                  tstrxn%surf_area_function = MINERAL_SURF_AREA_F_POR_RATIO
                case('VOLUME_FRACTION_RATIO')
                  tstrxn%surf_area_function = MINERAL_SURF_AREA_F_VF_RATIO
                case('POROSITY_VOLUME_FRACTION_RATIO')
                  tstrxn%surf_area_function = MINERAL_SURF_AREA_F_POR_VF_RATIO
                case('MINERAL_MASS')
                  tstrxn%surf_area_function = MINERAL_SURF_AREA_F_MNRL_MASS
                case default
                  error_string = error_string // ',' // &
                                 trim(keyword) // ',' // word
                  call InputKeywordUnrecognized(input,word, &
                                                error_string,option)
              end select
            case('SPECIFIC_SURFACE_AREA')
              call InputReadDouble(input,option,tstrxn%spec_surf_area)
              call InputErrorMsg(input,option,keyword,error_string)
              call InputReadAndConvertUnits(input,tstrxn%spec_surf_area, &
                        'm^2/kg',trim(error_string)//','//trim(keyword),option)
            case('NUCLEATION_KINETICS')
              tstrxn%nucleation => ReactionMnrlCreateNucleation()
              call InputReadWord(input,option,tstrxn%nucleation%name,PETSC_TRUE)
              call InputErrorMsg(input,option,keyword,error_string)
            case default
              call InputKeywordUnrecognized(input,keyword,error_string,option)
          end select
        enddo
        call InputPopBlock(input,option)
        ! Ensure that both forward and reverse rate constants have been set
        if ((Uninitialized(tstrxn%precipitation_rate_constant) .and. &
             Initialized(tstrxn%dissolution_rate_constant)) .or. &
            (Initialized(tstrxn%precipitation_rate_constant) .and. &
             Uninitialized(tstrxn%dissolution_rate_constant))) then
          option%io_buffer = 'Both forward and reverse rate constants must &
            &be specified if either is specified for kinetic mineral "' // &
            trim(cur_mineral%name) // '".'
          call PrintErrMsg(option)
        endif
        ! Ensure that correct surface area parameters are set
        if (tstrxn%surf_area_porosity_pwr > 0.d0) then
          select case(tstrxn%surf_area_function)
            case(MINERAL_SURF_AREA_F_POR_RATIO)
            case(MINERAL_SURF_AREA_F_POR_VF_RATIO)
            case(MINERAL_SURF_AREA_F_MNRL_MASS)
              option%io_buffer = 'SURFACE_AREA_FUNCTION MINERAL_MASS is &
                &not compatible with SURFACE_AREA_POROSITY_POWER for &
                &mineral "' // trim(cur_mineral%name) // '".'
              call PrintErrMsg(option)
            case default
              option%io_buffer = 'A SURFACE_AREA_FUNCTION must be set for &
                &mineral "' // trim(cur_mineral%name) // '" since &
                &SURFACE_AREA_POROSITY_POWER is specified.'
              call PrintErrMsg(option)
          end select
        endif
        if (tstrxn%surf_area_vol_frac_pwr > 0.d0) then
          select case(tstrxn%surf_area_function)
            case(MINERAL_SURF_AREA_F_VF_RATIO)
            case(MINERAL_SURF_AREA_F_POR_VF_RATIO)
            case(MINERAL_SURF_AREA_F_MNRL_MASS)
              option%io_buffer = 'SURFACE_AREA_FUNCTION MINERAL_MASS is &
                &not compatible with SURFACE_AREA_VOL_FRAC_POWER for &
                &mineral "' // trim(cur_mineral%name) // '".'
              call PrintErrMsg(option)
            case default
              option%io_buffer = 'A SURFACE_AREA_FUNCTION must be set for &
                &mineral "' // trim(cur_mineral%name) // '" since &
                &SURFACE_AREA_VOL_FRAC_POWER is specified.'
              call PrintErrMsg(option)
          end select
        endif
        select case(tstrxn%surf_area_function)
          case(MINERAL_SURF_AREA_F_POR_RATIO)
            if (Uninitialized(tstrxn%surf_area_porosity_pwr)) then
              option%io_buffer = 'A SURFACE_AREA_POROSITY_POWER must be &
                &specified for mineral "' // trim(cur_mineral%name) // '".'
              call PrintErrMsg(option)
            endif
          case(MINERAL_SURF_AREA_F_VF_RATIO)
            if (Uninitialized(tstrxn%surf_area_vol_frac_pwr)) then
              option%io_buffer = 'A SURFACE_AREA_VOL_FRAC_POWER must be &
                &specified for mineral "' // trim(cur_mineral%name) // '".'
              call PrintErrMsg(option)
            endif
          case(MINERAL_SURF_AREA_F_POR_VF_RATIO)
            if (Uninitialized(tstrxn%surf_area_porosity_pwr)) then
              option%io_buffer = 'A SURFACE_AREA_POROSITY_POWER must be &
                &specified for mineral "' // trim(cur_mineral%name) // '".'
              call PrintErrMsg(option)
            endif
            if (Uninitialized(tstrxn%surf_area_vol_frac_pwr)) then
              option%io_buffer = 'A SURFACE_AREA_VOL_FRAC_POWER must be &
                &specified for mineral "' // trim(cur_mineral%name) // '".'
              call PrintErrMsg(option)
            endif
          case(MINERAL_SURF_AREA_F_MNRL_MASS)
            if (Uninitialized(tstrxn%spec_surf_area)) then
              option%io_buffer = 'A SPECIFIC_SURFACE_AREA must be &
                &specified for mineral "' // trim(cur_mineral%name) // '".'
              call PrintErrMsg(option)
            endif
        end select
        ! Loop over prefactors and set kinetic rates and activation energies
        ! equal to the "outer" values if zero.
        cur_prefactor => tstrxn%prefactor
        do
          if (.not.associated(cur_prefactor)) exit
          ! if not initialized
          if (Uninitialized(cur_prefactor%precipitation_rate_constant)) then
            cur_prefactor%precipitation_rate_constant = &
              tstrxn%precipitation_rate_constant
            if (Uninitialized(cur_prefactor%precipitation_rate_constant)) then
              option%io_buffer = 'Both outer and inner prefactor forward &
                &rate constants uninitialized for kinetic mineral ' // &
                trim(cur_mineral%name) // '.'
              call PrintErrMsg(option)
            endif
          endif
          if (Uninitialized(cur_prefactor%dissolution_rate_constant)) then
            cur_prefactor%dissolution_rate_constant = &
              tstrxn%dissolution_rate_constant
            if (Uninitialized(cur_prefactor%dissolution_rate_constant)) then
              option%io_buffer = 'Both outer and inner prefactor reverse &
                &rate constants uninitialized for kinetic mineral ' // &
                trim(cur_mineral%name) // '.'
              call PrintErrMsg(option)
            endif
          endif
          if (Uninitialized(cur_prefactor%activation_energy)) then
            cur_prefactor%activation_energy = tstrxn%activation_energy
          endif
          cur_prefactor => cur_prefactor%next
        enddo
        if (.not.associated(tstrxn%prefactor) .and. &
            Uninitialized(tstrxn%precipitation_rate_constant)) then
          option%io_buffer = 'A rate constant must be defined for kinetic &
            &mineral "' // trim(cur_mineral%name) // '".'
          call PrintErrMsg(option)
        endif
        ! add tst rxn
        if (.not.associated(cur_mineral%tstrxn)) then
          cur_mineral%tstrxn => tstrxn
        else ! append to end of list
          cur_tstrxn => cur_mineral%tstrxn
          do
            if (.not.associated(cur_tstrxn%next)) then
              cur_tstrxn%next => tstrxn
              exit
            else
              cur_tstrxn => cur_tstrxn%next
            endif
          enddo
        endif
        cur_mineral%id = abs(cur_mineral%id)
        exit
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (.not.found) then
      option%io_buffer = 'Mineral "' // trim(name) // '" specified under &
        &CHEMISTRY,MINERAL_KINETICS not found in list of available minerals.'
      call PrintErrMsg(option)
    endif
  enddo
  call InputPopBlock(input,option)

  cur_mineral => mineral%mineral_list
  imnrl = 0
  do
    if (.not.associated(cur_mineral)) exit
    if (cur_mineral%id < 0 .and. &
        cur_mineral%itype == MINERAL_KINETIC) then
      option%io_buffer = 'No rate provided in input file for mineral: ' // &
               trim(cur_mineral%name) // '.'
      call PrintErrMsg(option)
    endif
    if (associated(cur_mineral%tstrxn)) then
      imnrl = imnrl + 1
!geh  reaction%kinmnrl_names(imnrl) = cur_mineral%name
    endif
    cur_mineral => cur_mineral%next
  enddo

  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = abs(cur_mineral%id)
    cur_mineral => cur_mineral%next
  enddo

end subroutine ReactionMnrlReadKinetics

! ************************************************************************** !

subroutine ReactionMnrlReadRateConstant(input,mineral_name, &
                                        rate_constant,error_string, &
                                        error_string2,option)
  !
  ! Reads mineral kinetics
  !
  ! Author: Glenn Hammond
  ! Date: 10/16/08
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use Units_module
  use Reaction_Database_Aux_module
  use Utility_module

  implicit none

  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: mineral_name
  PetscReal :: rate_constant
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=*) :: error_string2
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: internal_units

  call InputReadDouble(input,option,rate_constant)
  if (rate_constant < 0.d0) then
    rate_constant = 10.d0**rate_constant
  endif
  call InputErrorMsg(input,option,'RATE_CONSTANT',error_string)
  ! read units if they exist
  internal_units = 'mol/m^2-sec'
  call InputReadWord(input,option,word,PETSC_TRUE)
  if (InputError(input)) then
    input%err_buf = trim(mineral_name) // ' ' // trim(error_string2) // &
      ' RATE_CONSTANT UNITS'
    call InputDefaultMsg(input,option)
  else
    rate_constant = rate_constant * &
      UnitsConvertToInternal(word,internal_units, &
                             trim(error_string)//','//trim(error_string2)// &
                             ',RATE_CONSTANT',option)
  endif

end subroutine ReactionMnrlReadRateConstant

! ************************************************************************** !

subroutine ReactionMnrlReadMassActOverride(mineral,input,option)
  !
  ! Reads parameters for overriding mineral mass action and log Ks
  !
  ! Author: Glenn Hammond
  ! Date: 11/18/24
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use Reaction_Database_Aux_module
  use Utility_module

  implicit none

  type(mineral_type) :: mineral
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: name

  type(mineral_rxn_type), pointer :: cur_mineral
  type(mass_action_override_type), pointer :: mass_action_override
  PetscBool :: found
  PetscInt :: icount

  input%ierr = INPUT_ERROR_NONE
  icount = 0
  call InputPushBlock(input,'OVERRIDE_MINERAL_MASS_ACTION',option)
  do

    error_string = 'CHEMISTRY,OVERRIDE_MINERAL_MASS_ACTION'
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,name,PETSC_TRUE)
    call InputErrorMsg(input,option,'mineral name',error_string)

    cur_mineral => mineral%mineral_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_mineral)) exit
      if (StringCompare(cur_mineral%name,name,MAXWORDLENGTH)) then
        found = PETSC_TRUE
        mass_action_override => ReactionDBCreateMassActOverride()
        cur_mineral%mass_action_override => mass_action_override
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,error_string)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,keyword)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(keyword)

          select case(keyword)
            case('REACTION')
              mass_action_override%reaction_string = trim(input%buf)
            case('LOGK')
              call DeallocateArray(mass_action_override%logK)
              allocate(mass_action_override%logK(100))
              mass_action_override%logK = UNINITIALIZED_DOUBLE
              icount = 0
              do
                icount = icount + 1
                call InputReadDouble(input,option, &
                                      mass_action_override%logK(icount))
                if (icount == 1) then ! have to read at least 1
                  call InputErrorMsg(input,option,keyword,error_string)
                endif
                if (InputError(input)) then
                  icount = icount-1 ! decrement
                  exit
                endif
              enddo
              mass_action_override%logK => &
                TruncateArray(mass_action_override%logK,icount)
            case default
              call InputKeywordUnrecognized(input,keyword,error_string, &
                                            option)
          end select
        enddo
        call InputPopBlock(input,option)
        exit
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (.not.found) then
      option%io_buffer = 'Mineral "' // trim(name) // '" specified under ' // &
        trim(error_string) // ' not found in list of available minerals.'
      call PrintErrMsg(option)
    endif
  enddo
  call InputPopBlock(input,option)

end subroutine ReactionMnrlReadMassActOverride

! ************************************************************************** !

subroutine ReactionMnrlReadNucleation(mineral,input,option)
  !
  ! Reads parameters for overriding mineral mass action and log Ks
  !
  ! Author: Glenn Hammond
  ! Date: 11/18/24
  !
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(mineral_type) :: mineral
  type(input_type), pointer :: input
  type(option_type) :: option

  type(nucleation_type), pointer :: nucleation
  type(nucleation_type), pointer :: last_nucleation
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: keyword

  ! find last in list, if the list exists
  last_nucleation => mineral%nucleation_list
  if (associated(last_nucleation)) then
    do
      if (associated(last_nucleation%next)) then
        last_nucleation => last_nucleation%next
      else
        exit
      endif
    enddo
  endif

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,'MINERAL_NUCLEATION_KINETICS',option)
  do

    error_string = 'CHEMISTRY,MINERAL_NUCLEATION_KINETICS'
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    nucleation => ReactionMnrlCreateNucleation()
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,keyword,error_string)
    call StringToUpper(keyword)
    select case(keyword)
      case('CLASSICAL')
        nucleation%itype = MINERAL_NUCLEATION_CLASSICAL
      case('SIMPLIFIED')
        nucleation%itype = MINERAL_NUCLEATION_SIMPLIFIED
      case default
        error_string = trim(error_string)//',TYPE'
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
    call InputReadWord(input,option,nucleation%name,PETSC_TRUE)
    call InputErrorMsg(input,option,'nucleation name',error_string)

    call InputPushBlock(input,nucleation%name,option)
    do
      error_string = trim(error_string) // ',' // nucleation%name
      call InputReadPflotranString(input,option)
      call InputReadStringErrorMsg(input,option,error_string)
      if (InputCheckExit(input,option)) exit
      call InputReadCard(input,option,keyword)
      call InputErrorMsg(input,option,'keyword',error_string)
      call StringToUpper(keyword)

      select case(nucleation%itype)
        case(MINERAL_NUCLEATION_CLASSICAL)
          select case(keyword)
            case('RATE_CONSTANT')
              call ReactionMnrlReadRateConstant(input,nucleation%name, &
                                                nucleation%rate_constant, &
                                                error_string, &
                                                'NUCLEATION',option)
            case('GEOMETRIC_SHAPE_FACTOR')
              call InputReadDouble(input,option, &
                                  nucleation%geometric_shape_factor)
              call InputErrorMsg(input,option,keyword,error_string)
            case('HETEROGENEOUS_CORRECTION_FACTOR')
              call InputReadDouble(input,option, &
                                  nucleation%heterogenous_correction_factor)
              call InputErrorMsg(input,option,keyword,error_string)
            case('SURFACE_TENSION')
              call InputReadDouble(input,option,nucleation%surface_tension)
              call InputErrorMsg(input,option,keyword,error_string)
            case default
              call InputKeywordUnrecognized(input,keyword,error_string,option)
          end select
        case(MINERAL_NUCLEATION_SIMPLIFIED)
          select case(keyword)
            case('RATE_CONSTANT')
              call ReactionMnrlReadRateConstant(input,nucleation%name, &
                                                nucleation%rate_constant, &
                                                error_string, &
                                                'NUCLEATION',option)
            case('GAMMA')
              call InputReadDouble(input,option,nucleation%gamma)
              call InputErrorMsg(input,option,keyword,error_string)
            case default
              call InputKeywordUnrecognized(input,keyword,error_string,option)
          end select
      end select
    enddo
    call InputPopBlock(input,option)

    ! error checking
    select case(nucleation%itype)
      case(MINERAL_NUCLEATION_CLASSICAL)
        if (len_trim(nucleation%name) == 0 .or. &
            Uninitialized(nucleation%rate_constant) .or. &
            Uninitialized(nucleation%geometric_shape_factor) .or. &
            Uninitialized(nucleation%heterogenous_correction_factor) .or. &
            Uninitialized(nucleation%surface_tension)) then
          option%io_buffer = 'Uninitialized values in classical &
            &nucleation reaction "' // trim(nucleation%name) // '".'
          call PrintErrMsg(option)
        endif
      case(MINERAL_NUCLEATION_SIMPLIFIED)
        if (len_trim(nucleation%name) == 0 .or. &
            Uninitialized(nucleation%rate_constant) .or. &
            Uninitialized(nucleation%gamma)) then
          option%io_buffer = 'Uninitialized values in simplified &
            &nucleation reaction "' // trim(nucleation%name) // '".'
          call PrintErrMsg(option)
        endif
    end select

    if (associated(last_nucleation)) then
      last_nucleation%next => nucleation
    else
      mineral%nucleation_list => nucleation
    endif
    last_nucleation => nucleation
    nullify(nucleation)

  enddo
  call InputPopBlock(input,option)

end subroutine ReactionMnrlReadNucleation

! ************************************************************************** !

subroutine ReactionMnrlReadFromDatabase(mineral,num_dbase_temperatures, &
                                        input,option)
  !
  ! Reads mineral from database
  !
  ! Author: Glenn Hammond
  ! Date: 10/16/08
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use Reaction_Database_Aux_module

  implicit none

  type(mineral_rxn_type) :: mineral
  PetscInt :: num_dbase_temperatures
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: ispec
  PetscInt :: itemp
  PetscInt :: num_species_in_rxn

  ! read the molar volume
  call InputReadDouble(input,option,mineral%molar_volume)
  call InputErrorMsg(input,option,'MINERAL molar volume','DATABASE')
  ! convert from cm^3/mol to m^3/mol
  mineral%molar_volume = mineral%molar_volume*1.d-6
  ! create mineral reaction
  if (.not.associated(mineral%tstrxn)) then
    mineral%tstrxn => ReactionMnrlCreateTSTRxn()
  endif
  ! read the number of aqueous species in mineral rxn
  call InputReadInt(input,option,num_species_in_rxn)
  call InputErrorMsg(input,option,'Number of species in mineral reaction', &
                  'DATABASE')
  mineral%dbaserxn => &
    ReactionDBCreateRxn(num_species_in_rxn,num_dbase_temperatures)
  ! read in species and stoichiometries
  do ispec = 1, mineral%dbaserxn%reaction_equation%nspec
    call InputReadDouble(input,option,mineral%dbaserxn% &
                           reaction_equation%stoich(ispec))
    call InputErrorMsg(input,option,'MINERAL species stoichiometry','DATABASE')
    call InputReadQuotedWord(input,option,mineral%dbaserxn% &
                               reaction_equation%spec_name(ispec),PETSC_TRUE)
    call InputErrorMsg(input,option,'MINERAL species name','DATABASE')
  enddo
  !note: logKs read are pK so that K is in the denominator (i.e. Q/K)
  do itemp = 1, num_dbase_temperatures
    call InputReadDouble(input,option,mineral%dbaserxn%logK(itemp))
    call InputErrorMsg(input,option,'MINERAL logKs','DATABASE')
  enddo
  ! read the molar weight
  call InputReadDouble(input,option,mineral%molar_weight)
  call InputErrorMsg(input,option,'MINERAL molar weight','DATABASE')
  ! convert from g/mol to kg/mol
  mineral%molar_weight = mineral%molar_weight*1.d-3

end subroutine ReactionMnrlReadFromDatabase

! ************************************************************************** !

subroutine ReactionMnrlProcessConstraint(mineral,constraint_name, &
                                         constraint,option)
  !
  ! Initializes constraints based on mineral
  ! species in system
  !
  ! Author: Glenn Hammond
  ! Date: 01/07/13
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module
  use Units_module

  implicit none

  type(mineral_type), pointer :: mineral
  character(len=MAXWORDLENGTH) :: constraint_name
  type(mineral_constraint_type), pointer :: constraint
  type(option_type) :: option

  PetscBool :: found
  PetscInt :: imnrl, jmnrl
  PetscReal, parameter :: epsilon = 1.d-16

  character(len=MAXWORDLENGTH) :: mineral_names(mineral%nkinmnrl)
  character(len=MAXWORDLENGTH) :: constraint_vol_frac_string(mineral%nkinmnrl)
  character(len=MAXWORDLENGTH) :: constraint_area_string(mineral%nkinmnrl)
  character(len=MAXWORDLENGTH) :: constraint_area_units(mineral%nkinmnrl)
  PetscReal :: constraint_vol_frac(mineral%nkinmnrl)
  PetscReal :: constraint_area(mineral%nkinmnrl)
  PetscBool :: external_vol_frac_dataset(mineral%nkinmnrl)
  PetscBool :: external_area_dataset(mineral%nkinmnrl)
  character(len=MAXWORDLENGTH) :: units
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXWORDLENGTH) :: mineral_name
  PetscReal :: tempreal
  PetscReal :: molar_volume
  PetscReal :: molar_weight
  PetscReal :: specific_surface_area
  PetscReal :: constraint_surface_area

  if (.not.associated(constraint)) return

  mineral_names = ''
  constraint_vol_frac_string = ''
  constraint_area_string = ''
  external_vol_frac_dataset = PETSC_FALSE
  external_area_dataset = PETSC_FALSE
  do imnrl = 1, mineral%nkinmnrl
    found = PETSC_FALSE
    do jmnrl = 1, mineral%nkinmnrl
      if (StringCompare(constraint%names(imnrl), &
                        mineral%kinmnrl_names(jmnrl), &
                        MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
                'Mineral ' // trim(constraint%names(imnrl)) // &
                'from CONSTRAINT ' // trim(constraint_name) // &
                ' not found among kinetic minerals.'
      call PrintErrMsg(option)
    else
      constraint_vol_frac(jmnrl) = &
        constraint%constraint_vol_frac(imnrl)
      constraint_area(jmnrl) = &
        constraint%constraint_area(imnrl)
      mineral_names(jmnrl) = constraint%names(imnrl)
      constraint_vol_frac_string(jmnrl) = &
        constraint%constraint_vol_frac_string(imnrl)
      constraint_area_string(jmnrl) = &
        constraint%constraint_area_string(imnrl)
      constraint_area_units(jmnrl) = &
        constraint%constraint_area_units(imnrl)
      external_vol_frac_dataset(jmnrl) = &
        constraint%external_vol_frac_dataset(imnrl)
      external_area_dataset(jmnrl) = &
        constraint%external_area_dataset(imnrl)
    endif
  enddo
  constraint%names = mineral_names
  constraint%constraint_vol_frac = constraint_vol_frac
  constraint%constraint_area = constraint_area
  constraint%constraint_vol_frac_string = constraint_vol_frac_string
  constraint%constraint_area_string = constraint_area_string
  constraint%constraint_area_units = constraint_area_units
  constraint%external_vol_frac_dataset = external_vol_frac_dataset
  constraint%external_area_dataset = external_area_dataset

  ! set up constraint specific surface area conversion factor
  do imnrl = 1, mineral%nkinmnrl
    mineral_name = mineral%kinmnrl_names(imnrl)
    molar_volume = mineral%kinmnrl_molar_vol(imnrl)
    molar_weight = mineral%kinmnrl_molar_wt(imnrl)
    specific_surface_area = mineral%kinmnrl_spec_surf_area(imnrl)

    units = constraint%constraint_area_units(imnrl)
    call StringToLower(units)
    if (StringEndsWith(units,'m^3_mnrl')) then
      units = units(1:index(units,'m^3_mnrl')+2)
      constraint%area_units_type(imnrl) = MINERAL_SURF_AREA_PER_MNRL_VOL
      internal_units = 'm^2/m^3' ! m^2 mnrl/m^3 mnrl
    elseif (StringEndsWith(units,'g')) then
      constraint%area_units_type(imnrl) = MINERAL_SURF_AREA_PER_MNRL_MASS
      internal_units = 'm^2/kg' ! m^2 mnrl/kg mnrl
    elseif (StringEndsWith(units,'m^3')) then
      constraint%area_units_type(imnrl) = MINERAL_SURF_AREA_PER_BULK_VOL
      internal_units = 'm^2/m^3' ! m^2 mnrl/m^3 bulk
    else
      option%io_buffer = 'Unrecognized mineral specific surface area units &
        &for mineral ' // constraint%names(imnrl) // '.'
      call PrintErrMsg(option)
    endif
    ! convert to m^2/m^3 (mnrl or bulk) or m^2/kg mnrl
    !   tempreal is solely a conversion factor at this point
    !   DO NOT scale the constraint surface area
    tempreal = UnitsConvertToInternal(units,internal_units, &
                         trim(constraint%names(imnrl))// &
                         ',specific surface area', &
                         option)
    if (constraint%area_units_type(imnrl) == &
        MINERAL_SURF_AREA_PER_MNRL_MASS) then
      ! m^2 mnrl/kg mnrl -> m^2 mnrl/m^3 mnrl
      ! if a different specific surface area is specified in constraint, scale
      ! to be wrt the specific surface are specified under mineral kinetics
      if (Initialized(specific_surface_area)) then
        ! solely use in scaling below
        constraint_surface_area = tempreal * constraint%constraint_area(imnrl)
        if (.not.Equal(specific_surface_area,constraint_surface_area)) then
          option%io_buffer = 'Use of differing specific surface areas &
            &(mineral kinetics vs. constraint) is currently not &
            &supported: ' // StringWrite(specific_surface_area) // ' vs ' // &
            &StringWrite(constraint_surface_area)
          call PrintErrMsg(option)
        endif
        ! rescale if they differ
        tempreal = tempreal * constraint_surface_area / specific_surface_area
      endif

      if (molar_weight < epsilon .or. Equal(molar_weight,500.d0)) then
        option%io_buffer = 'Zero or undefined molar weight for mineral "' // &
          trim(mineral_name) // '" prevents specifying mineral specific &
          &surface area per mass mineral in constraint "' // &
          trim(constraint_name) // '".'
        call PrintErrMsg(option)
      endif
      if (molar_volume < epsilon .or. Equal(molar_volume,500.d0)) then
        option%io_buffer = 'Zero or undefined molar volume for mineral "' // &
          trim(mineral_name) // '" prevents specifying mineral specific &
          &surface area per mass mineral in constraint "' // &
          trim(constraint_name) // '".'
        call PrintErrMsg(option)
      endif
      ! m^2 mnrl/m^3 mnrl = m^2 mnrl/kg mnrl * kg mnrl/mol mnrl /
      !                     m^3 mnrl/mol mnrl
      tempreal = tempreal * molar_weight / molar_volume
    endif
    constraint%constraint_area_conv_factor(imnrl) = tempreal
    constraint%constraint_area_units(imnrl) = internal_units
    if (Initialized(constraint%constraint_vol_frac(imnrl))) then
      if (constraint%area_units_type(imnrl) == &
          MINERAL_SURF_AREA_PER_MNRL_MASS) then
        ! m^2 mnrl/m^3 mnrl -> m^2 mnrl/m^3 bulk
        tempreal = tempreal * constraint%constraint_vol_frac(imnrl)
      endif
    endif
    ! this is where we scale the constraint surface area (if initialized)
    if (Initialized(constraint%constraint_area(imnrl))) then
      ! tempreal converts input value to m^2 mnrl/m^3 bulk
      constraint%constraint_area(imnrl) = tempreal * &
        constraint%constraint_area(imnrl)

      ! check for zero initial surface areas
      call ReactionMnrlReportZeroSurfArea(imnrl, &
                                          constraint%constraint_area(imnrl), &
                                          mineral,constraint_name, &
                                          constraint,option)
    endif

  enddo

end subroutine ReactionMnrlProcessConstraint

! ************************************************************************** !

subroutine ReactionMnrlSetup(reaction,option)
  !
  ! Ensure that the mineral object is properly initialized.
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/25
  !
  use Option_module

  implicit none

  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  type(mineral_type), pointer :: mineral
  PetscBool :: error_found
  PetscInt :: imnrl

  mineral => reaction%mineral

  error_found = PETSC_FALSE
  do imnrl = 1, mineral%nkinmnrl
    select case(mineral%kinmnrl_surf_area_function(imnrl))
      case(MINERAL_SURF_AREA_F_MNRL_MASS)
        if (Uninitialized(mineral%kinmnrl_spec_surf_area(imnrl))) then
          error_found = PETSC_TRUE
          option%io_buffer = 'A specific surface area must be defined in &
            &in the MINERAL_KINETICS block for mineral "' // &
            trim(mineral%mineral_names(imnrl)) // '".'
          call PrintMsg(option)
        endif
    end select
  enddo

  if (error_found) then
    call PrintErrMsg(option,'Mineral kinetics missing parameters.')
  endif

end subroutine ReactionMnrlSetup


! ************************************************************************** !

subroutine ReactionMnrlKinetics(Res,Jac, &
                                compute_analytical_derivative,store_rate, &
                                rt_auxvar,global_auxvar,material_auxvar, &
                                reaction,option)
  !
  ! Wrapper routine for mineral kinetic reactions
  !
  ! Author: Glenn Hammond
  ! Date: 01/24/25

  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_analytical_derivative
  PetscBool :: store_rate
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: ln_sec(reaction%neqcplx)
  PetscReal :: ln_sec_act(reaction%neqcplx)
  PetscInt :: imnrl

  ! zero rates
  rt_auxvar%mnrl_rate(:) = 0.d0

  ! log concentration
  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  if (reaction%neqcplx > 0) then
    ln_sec = log(rt_auxvar%sec_molal)
    ln_sec_act = ln_sec+log(rt_auxvar%sec_act_coef)
  endif

  do imnrl = 1, reaction%mineral%nkinmnrl

    select case(reaction%mineral%kinmnrl_tst_itype(imnrl))
      case(MINERAL_KINETICS_TST_SIMPLE)
        call ReactionMnrlKineticRateTSTSimple(Res,Jac, &
                                compute_analytical_derivative,store_rate, &
                                imnrl,ln_conc,ln_act, &
                                rt_auxvar,global_auxvar,material_auxvar, &
                                reaction,option)
      case(MINERAL_KINETICS_TST_COMPLEX)
        call ReactionMnrlKineticRateTST(Res,Jac, &
                                compute_analytical_derivative,store_rate, &
                                imnrl,ln_conc,ln_act,ln_sec,ln_sec_act, &
                                rt_auxvar,global_auxvar,material_auxvar, &
                                reaction,option)
    end select
  enddo

  ! nucleation reactions
  if (associated(reaction%mineral%nucleation_array)) then
    call ReactionMnrlNucleationKinetics(Res,Jac, &
                                compute_analytical_derivative,store_rate, &
                                rt_auxvar,global_auxvar,material_auxvar, &
                                reaction,option)
  endif

end subroutine ReactionMnrlKinetics

! ************************************************************************** !

subroutine ReactionMnrlKineticRateTSTSimple(Res,Jac, &
                                      compute_derivative,store_rate, &
                                      imnrl,ln_conc,ln_act, &
                                      rt_auxvar,global_auxvar, &
                                      material_auxvar,reaction,option)
  !
  ! Computes the kinetic mineral precipitation/dissolution
  ! rates
  !
  ! Author: Glenn Hammond
  ! Date: 09/04/08
  !

  use Option_module
  use Material_Aux_module
  use Utility_module, only : Arrhenius

  implicit none

  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscBool :: store_rate
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscInt :: imnrl
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt :: i, j, icomp, jcomp
  PetscReal :: affinity_factor, sign_
  PetscReal :: Im, Im_const, dIm_dQK
  PetscReal :: QK, lnQK, dQK_dmj
  PetscBool :: precipitation
  PetscReal :: rate_constant
  type(mineral_type), pointer :: mineral

  mineral => reaction%mineral

  ! compute ion activity product
  lnQK = -mineral%kinmnrl_logK(imnrl)*LOG_TO_LN

  ! activity of water
  if (mineral%kinmnrlh2oid(imnrl) > 0) then
    lnQK = lnQK + mineral%kinmnrlh2ostoich(imnrl)* &
                  rt_auxvar%ln_act_h2o
  endif

  do i = 1, mineral%kinmnrlspecid(0,imnrl)
    icomp = mineral%kinmnrlspecid(i,imnrl)
    lnQK = lnQK + mineral%kinmnrlstoich(i,imnrl)*ln_act(icomp)
  enddo

  QK = exp(lnQK)

  affinity_factor = 1.d0-QK

  sign_ = sign(1.d0,affinity_factor) ! sign_ > 0 = dissolution

  if (rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0) then

    precipitation = (sign_ < 0.d0)

    if (precipitation) then
      rate_constant = mineral%kinmnrl_precip_rate_constant(imnrl)
    else
      rate_constant = mineral%kinmnrl_dissol_rate_constant(imnrl)
    endif
    if (.not.(rate_constant > 0.d0)) return

    ! compute rate
    ! rate: mol/m^2 mnrl/sec
    ! area: m^2 mnrl/m^3 bulk
    ! volume: m^3 bulk
    Im_const = -rt_auxvar%mnrl_area(imnrl)

    ! units: mol/sec/m^3 bulk
    Im = Im_const*sign_*dabs(affinity_factor)*rate_constant
    ! store volumetric rate to be used in updating mineral volume fractions
    ! at end of time step
    if (store_rate) then
      ! mol/sec/m^3 bulk
      rt_auxvar%mnrl_rate(imnrl) = rt_auxvar%mnrl_rate(imnrl) + Im
    endif
  else ! rate is already zero by default; move on to next mineral
    return
  endif

  ! scale Im_const by volume for calculating derivatives below
  ! units: m^2 mnrl
  Im_const = Im_const*material_auxvar%volume

  ! convert rate from volumetric (mol/sec/m^3 bulk) to mol/sec
  ! units: mol/sec
  Im = Im*material_auxvar%volume

  do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
    icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
    Res(icomp) = Res(icomp) + mineral%kinmnrlstoich_in_residual(i,imnrl)*Im
  enddo

  if (.not. compute_derivative) return

  ! calculate derivatives of rate with respect to free
  ! units = mol/sec
  dIm_dQK = -Im_const*rate_constant

  ! derivatives with respect to primary species in reaction quotient
  do j = 1, mineral%kinmnrlspecid(0,imnrl)
    jcomp = mineral%kinmnrlspecid(j,imnrl)
    ! unit = kg water/mol
    dQK_dmj = mineral%kinmnrlstoich(j,imnrl)*QK*exp(-ln_conc(jcomp))
    do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
      icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
      ! units = (mol/sec)*(kg water/mol) = kg water/sec
      Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                          mineral%kinmnrlstoich_in_residual(i,imnrl)* &
                          dIm_dQK*dQK_dmj
    enddo
  enddo

end subroutine ReactionMnrlKineticRateTSTSimple

! ************************************************************************** !

subroutine ReactionMnrlKineticRateTST(Res,Jac, &
                                      compute_derivative,store_rate, &
                                      imnrl,ln_conc,ln_act,ln_sec,ln_sec_act, &
                                      rt_auxvar,global_auxvar, &
                                      material_auxvar,reaction,option)
  !
  ! Computes the kinetic mineral precipitation/dissolution
  ! rates
  !
  ! Author: Glenn Hammond
  ! Date: 09/04/08
  !

  use Option_module
  use Material_Aux_module
  use Utility_module, only : Arrhenius

  implicit none

  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscBool :: store_rate
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  PetscInt :: imnrl
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: ln_sec(reaction%neqcplx)
  PetscReal :: ln_sec_act(reaction%neqcplx)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt :: i, j, icomp, jcomp
  PetscInt :: ipref, ipref_species
  ! I am assuming a maximum of 10 prefactors and 5 species per prefactor
  PetscReal :: tempreal
  PetscReal :: dspec_dprimary
  PetscReal :: affinity_factor, sign_
  PetscReal :: Im, Im_const, dIm_dQK
  PetscReal :: QK, lnQK, lnQK2, dQK_dmj
  PetscReal :: den
  PetscReal :: ln_spec_act, spec_act_coef
  PetscReal :: ln_prefactor, ln_numerator, ln_denominator
  PetscReal :: prefactor(10), ln_prefactor_spec(5,10)
  PetscReal :: sum_prefactor_rate
  PetscReal :: dIm_dsum_prefactor_rate, dIm_dspec
  PetscReal :: dprefactor_dprefactor_spec, dprefactor_spec_dspec
  PetscReal :: dprefactor_spec_dspec_numerator
  PetscReal :: dprefactor_spec_dspec_denominator
  PetscReal :: denominator
  PetscInt ::  icplx
  PetscReal :: ln_gam_m_beta
  PetscReal, parameter :: TREF = 25.d0
  PetscBool :: precipitation
  PetscReal :: rate_constant
  PetscReal :: arrhenius_factor
  type(mineral_type), pointer :: mineral

  mineral => reaction%mineral

  ! compute ion activity product
  lnQK = -mineral%kinmnrl_logK(imnrl)*LOG_TO_LN

  ! activity of water
  if (mineral%kinmnrlh2oid(imnrl) > 0) then
    lnQK = lnQK + mineral%kinmnrlh2ostoich(imnrl)* &
                  rt_auxvar%ln_act_h2o
  endif

  do i = 1, mineral%kinmnrlspecid(0,imnrl)
    icomp = mineral%kinmnrlspecid(i,imnrl)
    lnQK = lnQK + mineral%kinmnrlstoich(i,imnrl)*ln_act(icomp)
  enddo

  QK = exp(lnQK)

  if (associated(mineral%kinmnrl_Temkin_const)) then
    if (associated(mineral%kinmnrl_mnrl_scale_factor)) then
      affinity_factor = 1.d0-QK**(1.d0/ &
        (mineral%kinmnrl_mnrl_scale_factor(imnrl)* &
          mineral%kinmnrl_Temkin_const(imnrl)))
    else
      affinity_factor = 1.d0-QK**(1.d0/ &
                                mineral%kinmnrl_Temkin_const(imnrl))
    endif
  else if (associated(mineral%kinmnrl_mnrl_scale_factor)) then
      affinity_factor = 1.d0-QK**(1.d0/ &
        mineral%kinmnrl_mnrl_scale_factor(imnrl))
  else
    affinity_factor = 1.d0-QK
  endif

  sign_ = sign(1.d0,affinity_factor) ! sign_ > 0 = dissolution

  if (rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0) then

    precipitation = (sign_ < 0.d0)

!     check for supersaturation threshold for precipitation
!     if (associated(mineral%kinmnrl_affinity_threshold)) then
    if (mineral%kinmnrl_affinity_threshold(imnrl) > 0.d0) then
      if (precipitation .and. &
          QK < mineral%kinmnrl_affinity_threshold(imnrl)) return
    endif

!     check for rate limiter for precipitation
    if (mineral%kinmnrl_rate_limiter(imnrl) > 0.d0) then
      affinity_factor = affinity_factor/(1.d0+(1.d0-affinity_factor) &
        /mineral%kinmnrl_rate_limiter(imnrl))
    endif

    ! compute prefactor
    if (mineral%kinmnrl_num_prefactors(imnrl) > 0) then
      sum_prefactor_rate = 0.d0
      prefactor = 0.d0
      ln_prefactor_spec = 0.d0
      ! sum over parallel prefactors
      do ipref = 1, mineral%kinmnrl_num_prefactors(imnrl)
        ln_prefactor = 0.d0
        ! product of "monod" equations
        do ipref_species = 1, mineral%kinmnrl_prefactor_id(0,ipref,imnrl)
          icomp = mineral%kinmnrl_prefactor_id(ipref_species,ipref,imnrl)
          if (icomp > 0) then ! primary species
            ln_spec_act = ln_act(icomp)
          else ! secondary species (given a negative id to differentiate)
            ln_spec_act = ln_sec_act(-icomp)
          endif
          ln_numerator = &
            mineral%kinmnrl_pref_alpha(ipref_species,ipref,imnrl)* &
            ln_spec_act
          ln_denominator = log(1.d0 + &
            exp(log(mineral%kinmnrl_pref_atten_coef(ipref_species, &
                                                    ipref,imnrl)) + &
                mineral%kinmnrl_pref_beta(ipref_species,ipref,imnrl)* &
                ln_spec_act))
          ln_prefactor = ln_prefactor + ln_numerator
          ln_prefactor = ln_prefactor - ln_denominator
          ln_prefactor_spec(ipref_species,ipref) = &
            ln_numerator - ln_denominator
        enddo
        prefactor(ipref) = exp(ln_prefactor)
        arrhenius_factor = 1.d0
        if (mineral%kinmnrl_pref_activation_energy(ipref,imnrl) > 0.d0) then
          arrhenius_factor = &
            Arrhenius(mineral%kinmnrl_pref_activation_energy(ipref,imnrl),&
                      global_auxvar%temp,TREF)
        endif
        if (precipitation) then
          rate_constant = mineral%kinmnrl_pref_precip_rate_const(ipref,imnrl)
        else
          rate_constant = mineral%kinmnrl_pref_dissol_rate_const(ipref,imnrl)
        endif
        sum_prefactor_rate = sum_prefactor_rate + &
          prefactor(ipref)*rate_constant*arrhenius_factor
      enddo
    else
      arrhenius_factor = 1.d0
      if (mineral%kinmnrl_activation_energy(imnrl) > 0.d0) then
        arrhenius_factor = &
          Arrhenius(mineral%kinmnrl_activation_energy(imnrl),&
                    global_auxvar%temp,TREF)
      endif
      if (precipitation) then
        rate_constant = mineral%kinmnrl_precip_rate_constant(imnrl)
      else
        rate_constant = mineral%kinmnrl_dissol_rate_constant(imnrl)
      endif
      sum_prefactor_rate = rate_constant*arrhenius_factor
    endif
    if (.not.(sum_prefactor_rate > 0.d0)) return

    ! compute rate
    ! rate: mol/m^2 mnrl/sec
    ! area: m^2 mnrl/m^3 bulk
    ! volume: m^3 bulk
    Im_const = -rt_auxvar%mnrl_area(imnrl)
    if (associated(mineral%kinmnrl_mnrl_scale_factor)) then
      Im_const = Im_const/mineral%kinmnrl_mnrl_scale_factor(imnrl)
    endif

    ! units: mol/sec/m^3 bulk
    if (associated(mineral%kinmnrl_affinity_power)) then
      ! Im_const: m^2 mnrl/m^3 bulk
      ! sum_prefactor_rate: mol/m^2 mnrl/sec
      Im = Im_const*sign_* &
            dabs(affinity_factor)**mineral%kinmnrl_affinity_power(imnrl)* &
            sum_prefactor_rate
    else
      Im = Im_const*sign_*dabs(affinity_factor)*sum_prefactor_rate
    endif
    ! store volumetric rate to be used in updating mineral volume fractions
    ! at end of time step
    if (store_rate) then
      ! mol/sec/m^3 bulk
      rt_auxvar%mnrl_rate(imnrl) = rt_auxvar%mnrl_rate(imnrl) + Im
    endif
  else ! rate is already zero by default; move on to next mineral
    return
  endif

  ! scale Im_const by volume for calculating derivatives below
  ! units: m^2 mnrl
  Im_const = Im_const*material_auxvar%volume

  ! convert rate from volumetric (mol/sec/m^3 bulk) to mol/sec
  ! units: mol/sec
  Im = Im*material_auxvar%volume

  do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
    icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
    Res(icomp) = Res(icomp) + mineral%kinmnrlstoich_in_residual(i,imnrl)*Im
  enddo

  if (.not. compute_derivative) return

  ! calculate derivatives of rate with respect to free
  ! units = mol/sec
  if (associated(mineral%kinmnrl_affinity_power)) then
    tempreal = mineral%kinmnrl_affinity_power(imnrl)
    dIm_dQK = Im_const*tempreal* &
              dabs(affinity_factor)**(tempreal-1.d0)* &
              sum_prefactor_rate
    ! The separation of sign(affinity_factor) and dabs(affinity_factor)
    ! results in erroneous derivatives when the sign is negative. The
    ! call to sign(dIm_dQK,-Im_const*sum_prefactor_rate)
    ! corrects this issue.
    dIm_dQK = sign(dIm_dQK,-Im_const*sum_prefactor_rate)
  else
    dIm_dQK = -Im_const*sum_prefactor_rate
  endif

  if (associated(mineral%kinmnrl_Temkin_const)) then
    if (associated(mineral%kinmnrl_mnrl_scale_factor)) then
      dIm_dQK = dIm_dQK*(1.d0/(mineral%kinmnrl_mnrl_scale_factor(imnrl)* &
                          mineral%kinmnrl_Temkin_const(imnrl))) / &
                QK*(1.d0-affinity_factor)
    else
      dIm_dQK = dIm_dQK*(1.d0/mineral%kinmnrl_Temkin_const(imnrl))/QK* &
                (1.d0-affinity_factor)
    endif
  else if (associated(mineral%kinmnrl_mnrl_scale_factor)) then
    dIm_dQK = dIm_dQK*(1.d0/mineral%kinmnrl_mnrl_scale_factor(imnrl))/QK* &
              (1.d0-affinity_factor)
  endif

  ! derivatives with respect to primary species in reaction quotient
  if (mineral%kinmnrl_rate_limiter(imnrl) <= 0.d0) then
    do j = 1, mineral%kinmnrlspecid(0,imnrl)
      jcomp = mineral%kinmnrlspecid(j,imnrl)
      ! unit = kg water/mol
      dQK_dmj = mineral%kinmnrlstoich(j,imnrl)*QK*exp(-ln_conc(jcomp))
      do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
        icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                            mineral%kinmnrlstoich_in_residual(i,imnrl)* &
                            dIm_dQK*dQK_dmj
      enddo
    enddo

  else

    den = 1.d0+(1.d0-affinity_factor)/mineral%kinmnrl_rate_limiter(imnrl)
    do j = 1, mineral%kinmnrlspecid(0,imnrl)
      jcomp = mineral%kinmnrlspecid(j,imnrl)
      ! unit = kg water/mol
      dQK_dmj = mineral%kinmnrlstoich(j,imnrl)*QK*exp(-ln_conc(jcomp))
      do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
        icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
          mineral%kinmnrlstoich_in_residual(i,imnrl)*dIm_dQK*  &
          (1.d0 + QK/mineral%kinmnrl_rate_limiter(imnrl)/den)*dQK_dmj/den
      enddo
    enddo
  endif

  if (mineral%kinmnrl_num_prefactors(imnrl) > 0) then
    ! add contribution of derivative in prefactor - messy
    dIm_dsum_prefactor_rate = Im/sum_prefactor_rate
    ! summation over parallel reactions (prefactors)
    do ipref = 1, mineral%kinmnrl_num_prefactors(imnrl)
      arrhenius_factor = 1.d0
      if (mineral%kinmnrl_pref_activation_energy(ipref,imnrl) > 0.d0) then
        arrhenius_factor = &
          Arrhenius(mineral%kinmnrl_pref_activation_energy(ipref,imnrl),&
                    global_auxvar%temp,TREF)
      endif
      ! prefactor() saved in residual calc above
      ln_prefactor = log(prefactor(ipref))
      ! product of "monod" equations
      do ipref_species = 1, mineral%kinmnrl_prefactor_id(0,ipref,imnrl)
        ! derivative of 54 with respect to a single "monod" equation
        ! ln_prefactor_spec(,) saved in residual calc above
        dprefactor_dprefactor_spec = &
          exp(ln_prefactor-ln_prefactor_spec(ipref_species,ipref))
        icomp = mineral%kinmnrl_prefactor_id(ipref_species,ipref,imnrl)
        if (icomp > 0) then ! primary species
          ln_spec_act = ln_act(icomp)
          spec_act_coef = rt_auxvar%pri_act_coef(icomp)
        else ! secondary species
          ln_spec_act = ln_sec_act(-icomp)
          spec_act_coef = rt_auxvar%sec_act_coef(-icomp)
        endif
        ! derivative of numerator in eq. 54 with respect to species activity
        dprefactor_spec_dspec_numerator = &
          mineral%kinmnrl_pref_alpha(ipref_species,ipref,imnrl) * &
          exp(ln_prefactor_spec(ipref_species,ipref) - ln_spec_act)
        ln_gam_m_beta = mineral%kinmnrl_pref_beta(ipref_species, &
                                                  ipref,imnrl) * &
                        ln_spec_act
        ! denominator
        denominator = 1.d0 + &
            exp(log(mineral%kinmnrl_pref_atten_coef(ipref_species, &
                                                    ipref,imnrl)) + &
                ln_gam_m_beta)
        ! derivative of denominator in eq. 54 with respect to species
        ! activity
        dprefactor_spec_dspec_denominator = -1.d0 * &
          exp(ln_prefactor_spec(ipref_species,ipref)) / denominator * &
          mineral%kinmnrl_pref_atten_coef(ipref_species,ipref,imnrl) * &
          mineral%kinmnrl_pref_beta(ipref_species,ipref,imnrl) * &
          exp(ln_gam_m_beta - ln_spec_act)

        ! chain rule for derivative of "monod" equation
        dprefactor_spec_dspec = dprefactor_spec_dspec_numerator + &
          dprefactor_spec_dspec_denominator

        ! thus far the derivative is with respect to the activity, convert
        ! to with respect to molality
        dprefactor_spec_dspec = dprefactor_spec_dspec * spec_act_coef

        if (precipitation) then
          rate_constant = mineral%kinmnrl_pref_precip_rate_const(ipref,imnrl)
        else
          rate_constant = mineral%kinmnrl_pref_dissol_rate_const(ipref,imnrl)
        endif
        dIm_dspec = dIm_dsum_prefactor_rate * dprefactor_dprefactor_spec * &
                    dprefactor_spec_dspec * rate_constant * arrhenius_factor

        if (icomp > 0) then
          ! add derivative for primary species
          do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
            jcomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
            ! units = (mol/sec)*(kg water/mol) = kg water/sec
            Jac(jcomp,icomp) = Jac(jcomp,icomp) + &
                                mineral%kinmnrlstoich_in_residual(i,imnrl)* &
                                dIm_dspec
          enddo
        else ! secondary species -- have to calculate the derivative
          ! have to recalculate the reaction quotient (QK) for secondary
          ! species
          icplx = -icomp

          ! compute secondary species concentration
          lnQK2 = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

          ! activity of water
          if (reaction%eqcplxh2oid(icplx) > 0) then
            lnQK2 = lnQK2 + reaction%eqcplxh2ostoich(icplx) * &
                            rt_auxvar%ln_act_h2o
          endif

          do i = 1, reaction%eqcplxspecid(0,icplx)
            icomp = reaction%eqcplxspecid(i,icplx)
            lnQK2 = lnQK2 + reaction%eqcplxstoich(i,icplx)*ln_act(icomp)
          enddo
          ! add contribution to derivatives secondary prefactor with
          ! respect to free
          do j = 1, reaction%eqcplxspecid(0,icplx)
            jcomp = reaction%eqcplxspecid(j,icplx)
            dspec_dprimary = reaction%eqcplxstoich(j,icplx) * &
                              exp(lnQK2-ln_conc(jcomp)) / &
                              rt_auxvar%sec_act_coef(icplx)
            do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
              icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
              Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                mineral%kinmnrlstoich_in_residual(i,imnrl)* &
                dIm_dspec*dspec_dprimary
            enddo
          enddo
        endif
      enddo
    enddo  ! loop over prefactors
  endif

end subroutine ReactionMnrlKineticRateTST

! ************************************************************************** !

subroutine ReactionMnrlNucleationKinetics(Res,Jac, &
                                          compute_derivative,store_rate, &
                                          rt_auxvar,global_auxvar, &
                                          material_auxvar,reaction,option)
  !
  ! Calculates mineral precipitation nucleation rates.
  !
  ! Author: Glenn Hammond
  ! Date: 01/23/25
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscBool :: store_rate
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  type(mineral_type), pointer :: mineral
  PetscInt :: nid
  PetscInt :: i, j, imnrl, icomp, jcomp
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: QK, lnQK
  PetscReal :: dQK_dmj, dlnQK_dQK
  PetscReal :: nucleation_rate, dnucleation_rate_dlnQK

  mineral => reaction%mineral

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  do imnrl = 1, mineral%nkinmnrl ! for each mineral

    nid = mineral%kinmnrl_nucleation_id(imnrl)
    if (nid == 0) cycle

    ! compute ion activity product
    lnQK = -mineral%kinmnrl_logK(imnrl)*LOG_TO_LN

    ! activity of water
    if (mineral%kinmnrlh2oid(imnrl) > 0) then
      lnQK = lnQK + mineral%kinmnrlh2ostoich(imnrl)* &
                    rt_auxvar%ln_act_h2o
    endif

    do i = 1, mineral%kinmnrlspecid(0,imnrl)
      icomp = mineral%kinmnrlspecid(i,imnrl)
      lnQK = lnQK + mineral%kinmnrlstoich(i,imnrl)*ln_act(icomp)
    enddo

    QK = exp(lnQK)

    call ReactionMnrlNucleation(mineral%nucleation_array(nid), &
                                mineral%kinmnrl_molar_vol(imnrl), &
                                global_auxvar%temp,lnQK, &
                                nucleation_rate,dnucleation_rate_dlnQK, &
                                option)
!    print *, 'nucleation rate: ', nucleation_rate
#if 0
print *, 'analytical derivative: ', dnucleation_rate_dlnQK
    call ReactionMnrlNucleation(mineral%nucleation_array(nid), &
                                mineral%kinmnrl_molar_vol(imnrl), &
                                global_auxvar%temp,lnQK+1.d-6*lnQK, &
                                dlnQK_dQK,dQK_dmj, &
                                option)
print *, 'numerical derivative: ',  (dlnQK_dQK-nucleation_rate)/(1.d-6*lnQK)
!stop
#endif
    if (store_rate) then
      rt_auxvar%mnrl_rate(imnrl) = rt_auxvar%mnrl_rate(imnrl) + &
                                   nucleation_rate
    endif
    do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
      icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
      Res(icomp) = Res(icomp) + &
                   mineral%kinmnrlstoich_in_residual(i,imnrl) * &
                   nucleation_rate
    enddo

    if (.not. compute_derivative) cycle

    dlnQK_dQK = 1.d0/QK
    do j = 1, mineral%kinmnrlspecid(0,imnrl)
      jcomp = mineral%kinmnrlspecid(j,imnrl)
      ! unit = kg water/mol
      dQK_dmj = mineral%kinmnrlstoich(j,imnrl)*QK*exp(-ln_conc(jcomp))
      do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
        icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
          mineral%kinmnrlstoich_in_residual(i,imnrl) * &
          dnucleation_rate_dlnQK * dlnQK_dQK * dQK_dmj
      enddo
    enddo

  enddo

end subroutine ReactionMnrlNucleationKinetics

! ************************************************************************** !

subroutine ReactionMnrlNucleation(nucleation,molar_volume, &
                                  temperature_C,lnQK, &
                                  rate,drate_dlnQK,option)
  !
  ! Calculates individual mineral nucleation rates
  !
  ! Author: Glenn Hammond
  ! Date: Date: 01/23/25
  !
  use Option_module

  implicit none

  type(nucleation_type) :: nucleation
  PetscReal :: molar_volume
  PetscReal :: temperature_C
  PetscReal :: lnQK
  PetscReal :: rate
  PetscReal :: drate_dlnQK
  type(option_type) :: option

  PetscReal, parameter :: Avogadros_number = 6.02d23
  PetscReal :: T_K
  PetscReal :: exp_term, exp_term1, exp_term2

  if (lnQK <= 0.d0) then
    rate = 0.d0
    drate_dlnQK = 0.d0
    return
  endif

  T_K = temperature_C + T273K
  select case(nucleation%itype)
    case(MINERAL_NUCLEATION_CLASSICAL)
      exp_term1 = -1.d0 * nucleation%geometric_shape_factor * &
                  Avogadros_number * &
                  nucleation%heterogenous_correction_factor * &
                  molar_volume**2 * &
                  nucleation%surface_tension**3 / &
                  ((IDEAL_GAS_CONSTANT * T_K)**3)
      exp_term2 = 1.d0/lnQK
      exp_term = exp_term1 * exp_term2**2
      rate = nucleation%rate_constant * exp(exp_term)
      drate_dlnQK = rate*exp_term*(-2.d0*exp_term2)
    case(MINERAL_NUCLEATION_SIMPLIFIED)
      exp_term1 = -1.d0 * nucleation%gamma / (T_K**3)
      exp_term2 = 1.d0/lnQK
      exp_term = exp_term1 * exp_term2**2
      rate = nucleation%rate_constant * exp(exp_term)
      drate_dlnQK = rate*exp_term*(-2.d0*exp_term2)
  end select
#if 0
  print *, rate, lnQK
  print *, 'exp: ', exp_term, exp_term1, exp_term2
  print *, 'rate_constant: ', nucleation%rate_constant
  print *, 'lnQK: ', lnQK
  print *, 'geom: ', nucleation%geometric_shape_factor
  print *, 'hetero: ', nucleation%heterogenous_correction_factor
  print *, 'molvol: ', molar_volume
  print *, 'surf_tens: ', nucleation%surface_tension
  print *, 'TC: ', temperature_C
#endif

end subroutine ReactionMnrlNucleation

! ************************************************************************** !

function ReactionMnrlSaturationIndex(imnrl,rt_auxvar,global_auxvar, &
                                     reaction,option)
  !
  ! Calculates the mineral saturation index
  !
  ! Author: Glenn Hammond
  ! Date: 08/29/11
  !
  use Option_module

  implicit none

  type(option_type) :: option
  PetscInt :: imnrl
  class(reaction_rt_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  PetscReal :: ReactionMnrlSaturationIndex
  PetscInt :: i, icomp
  PetscReal :: lnQK
  PetscInt, parameter :: iphase = 1
  type(mineral_type), pointer :: mineral

  mineral => reaction%mineral

  if (.not.option%transport%isothermal_reaction) then
    call ReactionMnrlUpdateTempDepCoefs(global_auxvar%temp, &
                                   global_auxvar%pres(iphase), &
                                   reaction%mineral, &
                                   reaction%use_geothermal_hpt, &
                                   PETSC_TRUE,option)
  endif

  ! compute saturation
  lnQK = -mineral%mnrl_logK(imnrl)*LOG_TO_LN
  if (mineral%mnrlh2oid(imnrl) > 0) then
    lnQK = lnQK + mineral%mnrlh2ostoich(imnrl)*rt_auxvar%ln_act_h2o
  endif
  do i = 1, mineral%mnrlspecid(0,imnrl)
    icomp = mineral%mnrlspecid(i,imnrl)
    lnQK = lnQK + mineral%mnrlstoich(i,imnrl)* &
           log(rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp))
  enddo
  ReactionMnrlSaturationIndex = exp(lnQK)

end function ReactionMnrlSaturationIndex

! ************************************************************************** !

subroutine ReactionMnrlUpdateTempDepCoefs(temp,pres,mineral, &
                                          use_geothermal_hpt, &
                                          update_mnrl,option)
  !
  ! Updates temperature dependent coefficients for
  ! anisothermal simulations
  !
  ! Author: Glenn Hammond
  ! Date: 01/25/13
  !

  use Option_module

  implicit none

  PetscReal :: temp
  PetscReal :: pres
  type(mineral_type) :: mineral
  PetscBool :: use_geothermal_hpt
  PetscBool :: update_mnrl
  type(option_type) :: option

  if (.not.use_geothermal_hpt) then
    if (associated(mineral%kinmnrl_logKcoef)) then
      call ReactionAuxInterpolateLogK(mineral%kinmnrl_logKcoef, &
                                      mineral%kinmnrl_logK, &
                                      temp, &
                                      mineral%nkinmnrl)
    endif
    if (update_mnrl .and. associated(mineral%mnrl_logKcoef)) then
      call ReactionAuxInterpolateLogK(mineral%mnrl_logKcoef, &
                                      mineral%mnrl_logK, &
                                      temp, &
                                      mineral%nmnrl)
    endif
  else
    if (associated(mineral%kinmnrl_logKcoef)) then
      call ReactionAuxInterpolateLogK_hpt(mineral%kinmnrl_logKcoef, &
                                          mineral%kinmnrl_logK, &
                                          temp,pres,mineral%nkinmnrl)
    endif
    if (update_mnrl .and. associated(mineral%mnrl_logKcoef)) then
      call ReactionAuxInterpolateLogK_hpt(mineral%mnrl_logKcoef, &
                                          mineral%mnrl_logK, &
                                          temp,pres,mineral%nmnrl)
    endif
  endif

end subroutine ReactionMnrlUpdateTempDepCoefs

! ************************************************************************** !

subroutine ReactionMnrlUpdateSpecSurfArea(reaction,rt_auxvar, &
                                          material_auxvar, &
                                          porosity0,option)
  !
  ! Updates mineral specific surface area
  !
  ! Author: Glenn Hammond
  ! Date: 03/04/21
  !
  use Material_Aux_module
  use Option_module

  implicit none

  class(reaction_rt_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: porosity0
  type(option_type) :: option

  type(mineral_type), pointer :: mineral
  PetscInt :: imnrl
  PetscInt :: imnrl1
  PetscInt :: imnrl_armor
  PetscReal :: porosity_scale
  PetscReal :: volfrac_scale
  PetscReal :: mnrl_volfrac0
  PetscReal :: mnrl_volfrac
  PetscBool :: calc_porosity
  PetscBool :: calc_volfrac

  mineral => reaction%mineral

  do imnrl = 1, mineral%nkinmnrl

    select case(mineral%kinmnrl_surf_area_function(imnrl))
      case(MINERAL_SURF_AREA_F_POR_VF_RATIO, &
           MINERAL_SURF_AREA_F_POR_RATIO, &
           MINERAL_SURF_AREA_F_VF_RATIO)
        select case(mineral%kinmnrl_surf_area_function(imnrl))
          case(MINERAL_SURF_AREA_F_POR_RATIO)
            calc_porosity = PETSC_TRUE
            calc_volfrac = PETSC_FALSE
            volfrac_scale = 1.d0
          case(MINERAL_SURF_AREA_F_VF_RATIO)
            calc_porosity = PETSC_FALSE
            calc_volfrac = PETSC_TRUE
            porosity_scale = 1.d0
          case(MINERAL_SURF_AREA_F_POR_VF_RATIO)
            calc_porosity = PETSC_TRUE
            calc_volfrac = PETSC_TRUE
        end select
        if (calc_porosity) then
          porosity_scale = &
              ((1.d0-material_auxvar%porosity_base) / &
              (1.d0-porosity0))** &
              mineral%kinmnrl_surf_area_porosity_pwr(imnrl)
        endif
        if (calc_volfrac) then
          mnrl_volfrac0 = max(rt_auxvar%mnrl_volfrac0(imnrl), &
                              mineral%kinmnrl_vol_frac_epsilon(imnrl))
          mnrl_volfrac = max(rt_auxvar%mnrl_volfrac(imnrl), &
                            mineral%kinmnrl_vol_frac_epsilon(imnrl))
          if (mnrl_volfrac0 > 0.d0) then
            volfrac_scale = (mnrl_volfrac/mnrl_volfrac0)** &
                            mineral%kinmnrl_surf_area_vol_frac_pwr(imnrl)
          endif
        endif

        rt_auxvar%mnrl_area(imnrl) = &
            max(rt_auxvar%mnrl_area0(imnrl), &
                mineral%kinmnrl_surf_area_epsilon(imnrl)) * &
            porosity_scale*volfrac_scale

        if (reaction%update_armor_mineral_surface .and. &
            mineral%kinmnrl_armor_crit_vol_frac(imnrl) > 0.d0) then
          imnrl_armor = imnrl
          do imnrl1 = 1, mineral%nkinmnrl
            if (mineral%kinmnrl_armor_min_names(imnrl) == &
                mineral%kinmnrl_names(imnrl1)) then
              imnrl_armor = imnrl1
              exit
            endif
          enddo

          ! check for negative surface area armoring correction
          if (mineral%kinmnrl_armor_crit_vol_frac(imnrl) > &
              rt_auxvar%mnrl_volfrac(imnrl_armor)) then

            if (reaction%update_armor_mineral_surface_flag == 0) then
              ! surface unarmored
              rt_auxvar%mnrl_area(imnrl) = &
                rt_auxvar%mnrl_area(imnrl) * &
                ((mineral%kinmnrl_armor_crit_vol_frac(imnrl) &
                - rt_auxvar%mnrl_volfrac(imnrl_armor))/ &
                mineral%kinmnrl_armor_crit_vol_frac(imnrl))** &
                mineral%kinmnrl_surf_area_vol_frac_pwr(imnrl)
            else
              rt_auxvar%mnrl_area(imnrl) = &
                rt_auxvar%mnrl_area0(imnrl)
              reaction%update_armor_mineral_surface_flag = 0
            endif
          else
            rt_auxvar%mnrl_area(imnrl) = 0.d0
            reaction%update_armor_mineral_surface_flag = 1 ! surface armored
          endif
        endif
      case(MINERAL_SURF_AREA_F_MNRL_MASS)
        ! m^2 mnrl/m^3 bulk = m^2 mnrl/kg mnrl *    [specific surface area]
        !                     kg mnrl/mol mnrl /    [formula weight]
        !                     m^3 mnrl/mol mnrl *   [molar volume]
        !                     m^3 mnrl/m^3 bulk     [volume fraction]
        rt_auxvar%mnrl_area(imnrl) = &
          mineral%kinmnrl_spec_surf_area(imnrl) * &
          mineral%kinmnrl_molar_wt(imnrl) / &
          mineral%kinmnrl_molar_vol(imnrl) * &
          rt_auxvar%mnrl_volfrac(imnrl)
      case default
    end select
  enddo

end subroutine ReactionMnrlUpdateSpecSurfArea

! ************************************************************************** !

subroutine ReactionMnrlUpdateKineticState(rt_auxvar,global_auxvar, &
                                          material_auxvar,reaction, &
                                          kinetic_state_updated,option)
  !
  ! Update the mineral volume fraction and mass exchange due to mineral
  ! precipitation-dissolution
  !
  ! Author: Glenn Hammond
  ! Date: 08/18/23
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(reaction_rt_type) :: reaction
  PetscBool :: kinetic_state_updated
  type(option_type) :: option

  PetscInt, save :: icount = 0
  PetscInt :: imnrl, iaqspec, icomp
  PetscReal :: delta_volfrac
  PetscReal :: res(reaction%ncomp) ! has to be sized accurately
  PetscReal :: jac(1,1) ! strictly a dummy array
  PetscBool, parameter :: store_rate = PETSC_TRUE
  PetscBool, parameter :: compute_analytical_derivative = PETSC_FALSE

  if (reaction%mineral%nkinmnrl == 0) return

  kinetic_state_updated = PETSC_TRUE

  ! Updates the mineral rates, res is not needed
  call ReactionMnrlKinetics(res,jac, &
                            compute_analytical_derivative,store_rate, &
                            rt_auxvar,global_auxvar, &
                            material_auxvar,reaction,option)

  do imnrl = 1, reaction%mineral%nkinmnrl
    ! rate = mol/m^3/sec
    ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) *
    !                                mol_vol (m^3 mnrl/mol mnrl)
    delta_volfrac = rt_auxvar%mnrl_rate(imnrl)* &
                    reaction%mineral%kinmnrl_molar_vol(imnrl)* &
                    option%tran_dt
    rt_auxvar%mnrl_volfrac(imnrl) = rt_auxvar%mnrl_volfrac(imnrl) + &
                                    delta_volfrac
    if (rt_auxvar%mnrl_volfrac(imnrl) < 0.d0) &
      rt_auxvar%mnrl_volfrac(imnrl) = 0.d0

    ! CO2-specific
    if (option%transport%couple_co2) then
      do iaqspec = 1, reaction%mineral%kinmnrlspecid_in_residual(0,imnrl)
        icomp = reaction%mineral%kinmnrlspecid_in_residual(iaqspec,imnrl)
        !geh: this is problematic. once back to the flow side, the co2 src/sink
        !     rate depends on the previous flow time step size
        if (icomp == reaction%species_idx%pri_co2_id) then
          if (icount < 20) then
            icount = icount + 1
            call PrintMsg(option,'problematic mineral rate for co2 coupling')
          endif
          global_auxvar%reaction_rate(2) = &
            global_auxvar%reaction_rate(2) + &
            rt_auxvar%mnrl_rate(imnrl)*option%tran_dt * &
            reaction%mineral%kinmnrlstoich_in_residual(iaqspec,imnrl)
          cycle
        endif
      enddo
      ! water rate
      if (reaction%mineral%kinmnrlh2ostoich_in_residual(imnrl) /= 0) then
        global_auxvar%reaction_rate(1) = &
          global_auxvar%reaction_rate(1) + &
          rt_auxvar%mnrl_rate(imnrl)*option%tran_dt * &
          reaction%mineral%kinmnrlh2ostoich_in_residual(imnrl)
      endif
      if (option%iflowmode /= SCO2_MODE) global_auxvar%reaction_rate = &
                                   global_auxvar%reaction_rate / option%flow_dt
    endif
  enddo

end subroutine ReactionMnrlUpdateKineticState

! ************************************************************************** !

subroutine ReactionMnrlReportZeroSurfArea(imnrl,surface_area,mineral, &
                                          constraint_name, &
                                          mineral_constraint,option)
  !
  ! Reports an error is the initial mineral reactive surface area is zero
  ! and settings are not configure to overcome that
  !
  ! Author: Glenn Hammond
  ! Date: 06/19/25
  !

  use Option_module

  implicit none

  PetscInt :: imnrl
  PetscReal :: surface_area
  type(mineral_type), pointer :: mineral
  character(len=MAXWORDLENGTH) :: constraint_name
  type(mineral_constraint_type), pointer :: mineral_constraint
  type(option_type) :: option

  PetscReal, parameter :: epsilon = 1.d-40

  if (mineral_constraint%acknowledge_zero_surface_area) return

  ! to not throw this error, the initial surface area must be > 0 or
  ! nucleation must be specified for the mineral.

  if (surface_area < epsilon .and. &
      mineral%kinmnrl_nucleation_id(imnrl) == 0) then
    option%io_buffer = 'The volume fraction and/or specific &
      &surface area assigned to mineral "' // &
      trim(mineral%kinmnrl_names(imnrl)) // &
      '" in constraint "' // trim(constraint_name) // &
      '" results in a zero reactive mineral surface area. Please add &
      &ACKNOWLEDGE_ZERO_SURFACE_AREA to the CONSTRAINT block or &
      &assign a nucleation reaction to the mineral &
      &to avoid this error message.'
    call PrintErrMsg(option)
  endif

end subroutine ReactionMnrlReportZeroSurfArea

end module Reaction_Mineral_module
