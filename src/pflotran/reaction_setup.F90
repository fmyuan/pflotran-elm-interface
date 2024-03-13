module Reaction_Setup_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Reaction_module
  use Reaction_Aux_module
  use Reaction_Database_Aux_module

  implicit none

  private

  public :: ReactionSetupKinetics, &
            ReactionSetupSpecificSpecies, &
            ReactionSetupSpeciesSummary, &
            ReactionSetupPrimaryPrint

contains

! ************************************************************************** !

subroutine ReactionSetupKinetics(reaction,option)
  !
  ! Initializes non database related kinetic reactions
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/24
  !
  use CLM_Rxn_module
  use Option_module
  use Reaction_Equation_module
  use Reaction_Immobile_Aux_module
  use Reaction_Inhibition_Aux_module
  use Reaction_Isotherm_Aux_module
  use Reaction_Microbial_Aux_module
  use Reaction_Mineral_Aux_module
  use Reaction_Sandbox_module
  use String_module
  use Utility_module

  implicit none

  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  type(microbial_type), pointer :: microbial
  type(immobile_type), pointer :: immobile

  type(radioactive_decay_rxn_type), pointer :: cur_radiodecay_rxn
  type(microbial_rxn_type), pointer :: cur_microbial_rxn
  type(immobile_decay_rxn_type), pointer :: cur_immobile_decay_rxn
  type(dynamic_kd_rxn_type), pointer :: cur_dynamic_kd_rxn
  type(isotherm_link_type), pointer :: cur_isotherm_rxn, &
                                       sec_cont_cur_isotherm_rxn
  type(immobile_species_type), pointer :: cur_immobile_spec
  type(general_rxn_type), pointer :: cur_general_rxn
  type(monod_type), pointer :: cur_monod
  type(inhibition_type), pointer :: cur_inhibition
  type(reaction_equation_type), pointer :: reaction_equation

  PetscInt :: i, irxn, icomp
  PetscInt :: icplx, ipri_spec
  PetscInt :: temp_int
  PetscInt :: max_species_count
  PetscInt :: monod_count, inhibition_count, activation_energy_count
  PetscInt :: max_monod_count, max_inhibition_count
  PetscInt :: forward_count, max_forward_count
  PetscInt :: backward_count, max_backward_count
  PetscInt :: species_count
  PetscBool :: found

  microbial => reaction%microbial
  immobile => reaction%immobile

  ! immobile species
  immobile%nimmobile = ReactionImGetCount(immobile)
  if (immobile%nimmobile > 0) then
    allocate(immobile%names(immobile%nimmobile))
    immobile%names = ''
    allocate(immobile%print_me(immobile%nimmobile))
    immobile%print_me = PETSC_FALSE

    cur_immobile_spec => immobile%list
    temp_int = 0
    do
      if (.not.associated(cur_immobile_spec)) exit
      temp_int = temp_int + 1
      immobile%names(temp_int) = cur_immobile_spec%name
      immobile%print_me(temp_int) = cur_immobile_spec%print_me .or. &
                                   immobile%print_all
      cur_immobile_spec => cur_immobile_spec%next
    enddo
  endif
  reaction%offset_immobile = reaction%offset_aqueous + reaction%naqcomp

  ! radioactive decay reaction

  if (reaction%nradiodecay_rxn > 0) then

    ! process reaction equation into the database format
    cur_radiodecay_rxn => reaction%radioactive_decay_rxn_list
    do
      if (.not.associated(cur_radiodecay_rxn)) exit
      cur_radiodecay_rxn%reaction_equation => &
        ReactionEquationCreateFromString(cur_radiodecay_rxn%reaction, &
                                         reaction%naqcomp, &
                                         reaction%offset_aqueous, &
                                         reaction%primary_species_names, &
                                         reaction%nimcomp, &
                                         reaction%offset_immobile, &
                                         reaction%immobile%names, &
                                         PETSC_FALSE,option)
      cur_radiodecay_rxn => cur_radiodecay_rxn%next
    enddo
    nullify(cur_radiodecay_rxn)

    ! determine max # species for a given radiodecay rxn
    max_species_count = 0
    cur_radiodecay_rxn => reaction%radioactive_decay_rxn_list
    do
      if (.not.associated(cur_radiodecay_rxn)) exit

      ! zero count
      forward_count = 0

      ! max species in reaction
      species_count = cur_radiodecay_rxn%reaction_equation%nspec

      ! sum forward and reverse species
      reaction_equation => cur_radiodecay_rxn%reaction_equation
      do i = 1, reaction_equation%nspec
        if (reaction_equation%stoich(i) < 0.d0) then
          forward_count = forward_count + 1
        endif
      enddo

      if (forward_count > 1) then ! currently cannot have more than one species
        option%io_buffer = 'Cannot have more than one reactant in &
                           &radioactive decay reaction: (' // &
                           trim(cur_radiodecay_rxn%reaction) // ').'
        call PrintErrMsg(option)
      endif

      ! calculate maximum
      if (species_count > max_species_count) max_species_count = species_count

      cur_radiodecay_rxn => cur_radiodecay_rxn%next

    enddo
    nullify(cur_radiodecay_rxn)

    allocate(reaction%radiodecayspecid(0:max_species_count, &
                                       reaction%nradiodecay_rxn))
    reaction%radiodecayspecid = 0
    allocate(reaction%radiodecaystoich(max_species_count, &
                                       reaction%nradiodecay_rxn))
    reaction%radiodecaystoich = 0.d0
    allocate(reaction%radiodecayforwardspecid(reaction%nradiodecay_rxn))
    reaction%radiodecayforwardspecid = 0
    allocate(reaction%radiodecay_kf(reaction%nradiodecay_rxn))
    reaction%radiodecay_kf = 0.d0

    ! load the data into the compressed arrays
    irxn = 0
    cur_radiodecay_rxn => reaction%radioactive_decay_rxn_list
    do
      if (.not.associated(cur_radiodecay_rxn)) exit

      reaction_equation => cur_radiodecay_rxn%reaction_equation

      irxn = irxn + 1

      forward_count = 0
      backward_count = 0
      do i = 1, reaction_equation%nspec
        reaction%radiodecayspecid(i,irxn) = reaction_equation%specid(i)
        reaction%radiodecaystoich(i,irxn) = reaction_equation%stoich(i)
        if (reaction_equation%stoich(i) < 0.d0) then
          reaction%radiodecayforwardspecid(irxn) = reaction_equation%specid(i)
        endif
      enddo
      reaction%radiodecayspecid(0,irxn) = reaction_equation%nspec
      reaction%radiodecay_kf(irxn) = cur_radiodecay_rxn%rate_constant

      cur_radiodecay_rxn => cur_radiodecay_rxn%next

    enddo

  endif

  ! general reaction

  if (reaction%ngeneral_rxn > 0) then

    ! process reaction equation into the database format
    cur_general_rxn => reaction%general_rxn_list
    do
      if (.not.associated(cur_general_rxn)) exit
      cur_general_rxn%reaction_equation => &
        ReactionEquationCreateFromString(cur_general_rxn%reaction, &
                                         reaction%naqcomp, &
                                         reaction%offset_aqueous, &
                                         reaction%primary_species_names, &
                                         reaction%nimcomp, &
                                         reaction%offset_immobile, &
                                         reaction%immobile%names, &
                                         PETSC_FALSE,option)
      cur_general_rxn => cur_general_rxn%next
    enddo
    nullify(cur_general_rxn)

    ! determine max # species, forward species and backward species
    !  for a given general rxn
    max_species_count = 0
    max_forward_count = 0
    max_backward_count = 0
    cur_general_rxn => reaction%general_rxn_list
    do
      if (.not.associated(cur_general_rxn)) exit

      ! zero count
      forward_count = 0
      backward_count = 0

      ! max species in reaction
      species_count = cur_general_rxn%reaction_equation%nspec

      ! sum forward and reverse species
      reaction_equation => cur_general_rxn%reaction_equation
      do i = 1, reaction_equation%nspec
        if (reaction_equation%stoich(i) < 0.d0) then
          forward_count = forward_count + 1
        else if (reaction_equation%stoich(i) > 0.d0) then
          backward_count = backward_count + 1
        endif
      enddo

      ! calculate maximum
      if (forward_count > max_forward_count) max_forward_count = forward_count
      if (backward_count > max_backward_count) &
        max_backward_count = backward_count
      if (species_count > max_species_count) max_species_count = species_count

      cur_general_rxn => cur_general_rxn%next

    enddo
    nullify(cur_general_rxn)

    allocate(reaction%generalspecid(0:max_species_count,reaction%ngeneral_rxn))
    reaction%generalspecid = 0
    allocate(reaction%generalstoich(max_species_count,reaction%ngeneral_rxn))
    reaction%generalstoich = 0.d0
    allocate(reaction%generalforwardspecid(0:max_forward_count, &
                                           reaction%ngeneral_rxn))
    reaction%generalforwardspecid = 0
    allocate(reaction%generalforwardstoich(max_forward_count, &
                                           reaction%ngeneral_rxn))
    reaction%generalforwardstoich = 0.d0
    allocate(reaction%generalbackwardspecid(0:max_backward_count, &
                                            reaction%ngeneral_rxn))
    reaction%generalbackwardspecid = 0
    allocate(reaction%generalbackwardstoich(max_backward_count, &
                                            reaction%ngeneral_rxn))
    reaction%generalbackwardstoich = 0.d0
    allocate(reaction%generalh2oid(reaction%ngeneral_rxn))
    reaction%generalh2oid = 0
    allocate(reaction%generalh2ostoich(reaction%ngeneral_rxn))
    reaction%generalh2ostoich = 0.d0
    allocate(reaction%general_kf(reaction%ngeneral_rxn))
    reaction%general_kf = 0.d0
    allocate(reaction%general_kr(reaction%ngeneral_rxn))
    reaction%general_kr = 0.d0

    ! load the data into the compressed arrays
    irxn = 0
    cur_general_rxn => reaction%general_rxn_list
    do
      if (.not.associated(cur_general_rxn)) exit

      reaction_equation => cur_general_rxn%reaction_equation

      irxn = irxn + 1

      forward_count = 0
      backward_count = 0
      do i = 1, reaction_equation%nspec
        reaction%generalspecid(i,irxn) = reaction_equation%specid(i)
        reaction%generalstoich(i,irxn) = reaction_equation%stoich(i)
        if (reaction_equation%stoich(i) < 0.d0) then
          forward_count = forward_count + 1
          reaction%generalforwardspecid(forward_count,irxn) = &
            reaction_equation%specid(i)
          ! ensure that forward stoich is positive for rate expression
          reaction%generalforwardstoich(forward_count,irxn) = &
            dabs(reaction_equation%stoich(i))
        else if (reaction_equation%stoich(i) > 0.d0) then
          backward_count = backward_count + 1
          reaction%generalbackwardspecid(backward_count,irxn) = &
            reaction_equation%specid(i)
          reaction%generalbackwardstoich(backward_count,irxn) = &
            reaction_equation%stoich(i)
        endif
      enddo
      reaction%generalspecid(0,irxn) = reaction_equation%nspec
      reaction%generalforwardspecid(0,irxn) = forward_count
      reaction%generalbackwardspecid(0,irxn) = backward_count

      reaction%general_kf(irxn) = cur_general_rxn%forward_rate
      reaction%general_kr(irxn) = cur_general_rxn%backward_rate

      cur_general_rxn => cur_general_rxn%next

    enddo

  endif

  ! microbial reaction
  if (microbial%nrxn > 0) then

    ! process reaction equation into the database format
    max_species_count = 0
    max_monod_count = 0
    max_inhibition_count = 0
    monod_count = 0
    inhibition_count = 0
    activation_energy_count = 0
    cur_microbial_rxn => microbial%microbial_rxn_list
    do
      if (.not.associated(cur_microbial_rxn)) exit
      cur_microbial_rxn%reaction_equation => &
        ReactionEquationCreateFromString(cur_microbial_rxn%reaction, &
                                         reaction%naqcomp, &
                                         reaction%offset_aqueous, &
                                         reaction%primary_species_names, &
                                         reaction%nimcomp, &
                                         reaction%offset_immobile, &
                                         reaction%immobile%names, &
                                         PETSC_TRUE,option)
      if (cur_microbial_rxn%activation_energy > 0.d0) then
        activation_energy_count = activation_energy_count + 1
      endif
      temp_int = cur_microbial_rxn%reaction_equation%nspec
      if (temp_int > max_species_count) max_species_count = temp_int
      temp_int = ReactionMicrobGetMonodCount(cur_microbial_rxn)
      monod_count = monod_count + temp_int
      if (temp_int > max_monod_count) max_monod_count = temp_int
      temp_int = ReactionMicrobGetInhibtionCount(cur_microbial_rxn)
      inhibition_count = inhibition_count + temp_int
      if (temp_int > max_inhibition_count) max_inhibition_count = temp_int
      cur_microbial_rxn => cur_microbial_rxn%next
    enddo
    nullify(cur_microbial_rxn)
    if (max_inhibition_count > MAX_NUM_INHIBITION_TERMS) then
      option%io_buffer = 'The number of microbial inhibition terms (' // &
        StringWrite(max_inhibition_count) // ') exceeds &
        &MAX_NUM_INHIBITION_TERMS defined in reaction_microbial_aux.F90. &
        &Please increase the values.'
      call PrintErrMsg(option)
    endif

    ! rate constant
    allocate(microbial%rate_constant(microbial%nrxn))
    microbial%rate_constant = 0.d0

    ! activation_energy
    if (activation_energy_count > 0) then
      allocate(microbial%activation_energy(microbial%nrxn))
      microbial%activation_energy = 0.d0
    endif

    ! species ids and stoichiometry
    allocate(microbial%specid(0:max_species_count,microbial%nrxn))
    microbial%specid = 0
    allocate(microbial%stoich(max_species_count,microbial%nrxn))
    microbial%stoich = 0.d0

    ! biomass id and yield
    allocate(microbial%biomassid(microbial%nrxn))
    microbial%biomassid = 0
    allocate(microbial%biomass_yield(microbial%nrxn))
    microbial%biomass_yield = 0.d0

    ! linkage to monod and inhibition terms
    allocate(microbial%monodid(0:max_monod_count,microbial%nrxn))
    microbial%monodid = 0
    allocate(microbial%inhibitionid(0:max_inhibition_count, &
                                    microbial%nrxn))
    microbial%inhibitionid = 0

    ! monod
    allocate(microbial%monod_specid(monod_count))
    microbial%monod_specid = 0
    allocate(microbial%monod_K(monod_count))
    microbial%monod_K = 0.d0
    allocate(microbial%monod_Cth(monod_count))
    microbial%monod_Cth = 0.d0

    ! inhibition
    allocate(microbial%inhibition_specid(inhibition_count))
    microbial%inhibition_specid = 0
    allocate(microbial%inhibition_type(inhibition_count))
    microbial%inhibition_type = 0
    allocate(microbial%inhibition_C(inhibition_count))
    microbial%inhibition_C = 0.d0
    allocate(microbial%inhibition_C2(inhibition_count))
    microbial%inhibition_C2 = 0.d0

    ! load the data into the compressed arrays
    irxn = 0
    monod_count = 0
    inhibition_count = 0
    cur_microbial_rxn => microbial%microbial_rxn_list
    do
      if (.not.associated(cur_microbial_rxn)) exit

      reaction_equation => cur_microbial_rxn%reaction_equation

      irxn = irxn + 1

      microbial%rate_constant(irxn) = cur_microbial_rxn%rate_constant
      if (associated(microbial%activation_energy)) then
        microbial%activation_energy(irxn) = cur_microbial_rxn%activation_energy
      endif

      microbial%specid(0,irxn) = reaction_equation%nspec
      do i = 1, reaction_equation%nspec
        microbial%specid(i,irxn) = reaction_equation%specid(i)
        microbial%stoich(i,irxn) = reaction_equation%stoich(i)
      enddo

      if (associated(cur_microbial_rxn%biomass)) then
        ! try aqueous
        temp_int = &
          ReactionAuxGetPriSpecIDFromName(cur_microbial_rxn% &
                                            biomass%species_name, &
                                          reaction,PETSC_FALSE,option)
        ! temp_int will be UNINITIALIZED_INTEGER if not found
        if (Uninitialized(temp_int)) then
          ! check for biomass species in global immobile list
          temp_int = &
            StringFindEntryInList(cur_microbial_rxn%biomass%species_name, &
                                  immobile%names)
          ! temp_int will be zero if not found
          temp_int = -temp_int ! toggle to negative index for immobile
        endif
        if (temp_int == 0) then
          option%io_buffer = 'Biomass species "' // &
            trim(cur_microbial_rxn%biomass%species_name) // &
            ' not found among the primary aqueous or immobile species.'
          call PrintErrMsg(option)
        endif
        microbial%biomassid(irxn) = temp_int
        microbial%biomass_yield(irxn) = &
          cur_microbial_rxn%biomass%yield
        ! check for biomass species in microbial reaction
        temp_int = &
          StringFindEntryInList(cur_microbial_rxn%biomass%species_name, &
                                reaction_equation%spec_name)
        if (temp_int /= 0) then
          option%io_buffer = 'Biomass species "' // &
            trim(cur_microbial_rxn%biomass%species_name) // &
            ' should not be included in microbial reaction mass action &
            &expression.'
          if (microbial%biomassid(irxn) > 0) then
            option%io_buffer = trim(option%io_buffer) // ' Mobile biomass &
              &growth and decay is specified through a BIOMASS &
              &YIELD and a first-order aqueous GENERAL_REACTION, respectively.'
          else
            option%io_buffer = trim(option%io_buffer) // ' Immobile biomass &
              &growth and decay is specified through a BIOMASS &
              &YIELD and an IMMOBLE_DECAY_REACTION, respectively.'
          endif
          call PrintErrMsg(option)
        endif
      endif

      cur_monod => cur_microbial_rxn%monod
      do
        if (.not.associated(cur_monod)) exit
        monod_count = monod_count + 1

        ! increment # of monod reactions in microbial reaction
        microbial%monodid(0,irxn) = microbial%monodid(0,irxn) + 1
        ! set global id of this monod reaction
        microbial%monodid(microbial%monodid(0,irxn),irxn) = monod_count

        ! ensure that monod species exists in reaction expression
        temp_int = StringFindEntryInList(cur_monod%species_name, &
                                         reaction_equation%spec_name)
        if (temp_int == 0) then
          option%io_buffer = 'Monod species "' // &
            trim(cur_monod%species_name) // ' not found in microbial reaction.'
          call PrintErrMsg(option)
        endif
        ! if species stoichiometry is > 0., it is a product and cannot be
        ! used in a monod expression.
        if (reaction_equation%stoich(temp_int) > 0.d0) then
          option%io_buffer = 'Monod species "' // &
            trim(cur_monod%species_name) // ' must be a reactant and not ' // &
            'a product in microbial reaction.'
          call PrintErrMsg(option)
        endif

        microbial%monod_specid(monod_count) = &
          ReactionAuxGetPriSpecIDFromName(cur_monod%species_name, &
                                          reaction,option)
        microbial%monod_K(monod_count) = cur_monod%half_saturation_constant
        microbial%monod_Cth(monod_count) = cur_monod%threshold_concentration
        cur_monod => cur_monod%next
      enddo

      cur_inhibition => cur_microbial_rxn%inhibition
      do
        if (.not.associated(cur_inhibition)) exit
        inhibition_count = inhibition_count + 1

        ! increment # of inhibition reactions in microbial reaction
        microbial%inhibitionid(0,irxn) = microbial%inhibitionid(0,irxn) + 1
        ! set global id of this inhibition reaction
        microbial%inhibitionid(microbial%inhibitionid(0,irxn),irxn) = &
          inhibition_count

        ! Check whether inhibition species exists in reaction expression
        ! If no, print warning.
        temp_int = StringFindEntryInList(cur_inhibition%species_name, &
                                         reaction_equation%spec_name)
        if (temp_int == 0) then
          option%io_buffer = 'Inhibition species "' // &
            trim(cur_inhibition%species_name) // &
            ' not found in microbial reaction.'
          call PrintWrnMsg(option)
        endif

        microbial%inhibition_specid(inhibition_count) = &
          ReactionAuxGetPriSpecIDFromName(cur_inhibition%species_name, &
                                          reaction,option)
        microbial%inhibition_type(inhibition_count) = &
          cur_inhibition%itype
        microbial%inhibition_C(inhibition_count) = &
          cur_inhibition%inhibition_constant
        microbial%inhibition_C2(inhibition_count) = &
          cur_inhibition%inhibition_constant2
        cur_inhibition => cur_inhibition%next
      enddo

      cur_microbial_rxn => cur_microbial_rxn%next

    enddo

  endif

  ! immobile decay reaction

  if (reaction%immobile%ndecay_rxn > 0) then

    allocate(reaction%immobile%decayspecid(reaction%immobile%ndecay_rxn))
    allocate(reaction%immobile%decay_rate_constant(reaction%immobile%ndecay_rxn))

    cur_immobile_decay_rxn => reaction%immobile%decay_rxn_list
    irxn = 0
    do
      if (.not.associated(cur_immobile_decay_rxn)) exit

      irxn = irxn + 1

      found = PETSC_FALSE
      do i = 1, reaction%immobile%nimmobile
        if (StringCompare(cur_immobile_decay_rxn%species_name, &
                          reaction%immobile%names(i), &
                          MAXWORDLENGTH)) then
          reaction%immobile%decayspecid(irxn) = i
          found = PETSC_TRUE
          exit
        endif
      enddo
      if (.not.found) then
        option%io_buffer = 'Species "' // &
        trim(cur_immobile_decay_rxn%species_name) // &
        '" in immobile decay reaction not found among immobile species.'
        call PrintErrMsg(option)
      endif
      reaction%immobile%decay_rate_constant(irxn) = &
        cur_immobile_decay_rxn%rate_constant
      cur_immobile_decay_rxn => cur_immobile_decay_rxn%next
    enddo
    nullify(cur_immobile_decay_rxn)

  endif

  ! Smart Kd reactions

  if (reaction%neqdynamickdrxn > 0) then

    ! allocate arrays
    allocate(reaction%eqdynamickdspecid(reaction%neqdynamickdrxn))
    allocate(reaction%eqdynamickdrefspecid(reaction%neqdynamickdrxn))
    allocate(reaction%eqdynamickdrefspechigh(reaction%neqdynamickdrxn))
    allocate(reaction%eqdynamickdlow(reaction%neqdynamickdrxn))
    allocate(reaction%eqdynamickdhigh(reaction%neqdynamickdrxn))
    allocate(reaction%eqdynamickdpower(reaction%neqdynamickdrxn))
    reaction%eqdynamickdspecid = 0
    reaction%eqdynamickdrefspecid = 0
    reaction%eqdynamickdrefspechigh = 0.d0
    reaction%eqdynamickdlow = 0.d0
    reaction%eqdynamickdhigh = 0.d0
    reaction%eqdynamickdpower = 0.d0

    cur_dynamic_kd_rxn => reaction%dynamic_kd_rxn_list
    irxn = 0
    do
      if (.not.associated(cur_dynamic_kd_rxn)) exit

      irxn = irxn + 1

      found = PETSC_FALSE
      do i = 1, reaction%naqcomp
        if (StringCompare(cur_dynamic_kd_rxn%kd_species_name, &
                          reaction%primary_species_names(i), &
                          MAXWORDLENGTH)) then
          reaction%eqdynamickdspecid(irxn) = i
          found = PETSC_TRUE
          exit
        endif
      enddo
      if (.not.found) then
        option%io_buffer = 'KD species ' // &
             trim(cur_dynamic_kd_rxn%kd_species_name) // &
             ' in dynamic kd reaction not found among primary species list.'
        call PrintErrMsg(option)
      endif
      found = PETSC_FALSE
      do i = 1, reaction%naqcomp
        if (StringCompare(cur_dynamic_kd_rxn%ref_species_name, &
                          reaction%primary_species_names(i), &
                          MAXWORDLENGTH)) then
          reaction%eqdynamickdrefspecid(irxn) = i
          found = PETSC_TRUE
          exit
        endif
      enddo
      if (.not.found) then
        option%io_buffer = 'Reference species ' // &
           trim(cur_dynamic_kd_rxn%ref_species_name) // &
           ' in dynamic kd reaction not found among primary species list.'
        call PrintErrMsg(option)
      endif
      reaction%eqdynamickdrefspechigh(irxn) = &
         cur_dynamic_kd_rxn%ref_species_high
      reaction%eqdynamickdlow(irxn) = cur_dynamic_kd_rxn%KD_low
      reaction%eqdynamickdhigh(irxn) = cur_dynamic_kd_rxn%KD_high
      reaction%eqdynamickdpower(irxn) = cur_dynamic_kd_rxn%KD_power
      cur_dynamic_kd_rxn => cur_dynamic_kd_rxn%next
    enddo
  endif

  ! Kd reactions

  if (reaction%isotherm%neqkdrxn > 0) then

    call ReactionIsothermCreateRxn(reaction%isotherm%isotherm_rxn, &
      reaction%isotherm)
    ! allocate arrays
    allocate(reaction%isotherm%eqkdspecid(reaction%isotherm%neqkdrxn))
    reaction%isotherm%eqkdspecid = 0
    allocate(reaction%isotherm%eqisothermtype(reaction%isotherm%neqkdrxn))
    reaction%isotherm%eqisothermtype = 0
    allocate(reaction%isotherm%eqkdmineral(reaction%isotherm%neqkdrxn))
    reaction%isotherm%eqkdmineral = 0

    cur_isotherm_rxn => reaction%isotherm%isotherm_list

    if (option%use_sc) then
      call ReactionIsothermCreateRxn(reaction%isotherm%multicontinuum_isotherm_rxn, &
                             reaction%isotherm)
      sec_cont_cur_isotherm_rxn => &
        reaction%isotherm%multicontinuum_isotherm_list
    endif

    irxn = 0
    do
      if (.not.associated(cur_isotherm_rxn)) exit

      irxn = irxn + 1

      found = PETSC_FALSE
      do i = 1, reaction%naqcomp
        if (StringCompare(cur_isotherm_rxn%species_name, &
                          reaction%primary_species_names(i), &
                          MAXWORDLENGTH)) then
          reaction%isotherm%eqkdspecid(irxn) = i
          found = PETSC_TRUE
          exit
        endif
      enddo
      if (.not.found) then
        option%io_buffer = 'Species ' // trim(cur_isotherm_rxn%species_name) // &
                 ' in kd reaction &
                 & not found among primary species list.'
        call PrintErrMsg(option)
      endif
      reaction%isotherm%eqisothermtype(irxn) = cur_isotherm_rxn%itype
      ! associate mineral id
      if (len_trim(cur_isotherm_rxn%kd_mineral_name) > 1) then
        reaction%isotherm%eqkdmineral(irxn) = &
          ReactionMnrlGetKinMnrlIDFromName(cur_isotherm_rxn%kd_mineral_name, &
                                           reaction%mineral,option)
        if (reaction%isotherm%eqkdmineral(irxn) < 0) then
          option%io_buffer = 'Mineral ' // &
                             trim(cur_isotherm_rxn%kd_mineral_name) // &
                             ' listed in kd (linear sorption) &
                             &reaction not found in mineral list'
          call PrintErrMsg(option)
        endif
      endif
      reaction%isotherm%isotherm_rxn%eqisothermcoeff(irxn) = cur_isotherm_rxn%Kd
      reaction%isotherm%isotherm_rxn%eqisothermlangmuirb(irxn) = cur_isotherm_rxn%Langmuir_b
      reaction%isotherm%isotherm_rxn%eqisothermfreundlichn(irxn) = cur_isotherm_rxn%Freundlich_n

      cur_isotherm_rxn => cur_isotherm_rxn%next

      if (option%use_sc) then
        reaction%isotherm%multicontinuum_isotherm_rxn%eqisothermcoeff(irxn) = &
          sec_cont_cur_isotherm_rxn%Kd
        reaction%isotherm%multicontinuum_isotherm_rxn%eqisothermlangmuirb(irxn) = &
          sec_cont_cur_isotherm_rxn%Langmuir_b
        reaction%isotherm%multicontinuum_isotherm_rxn%eqisothermfreundlichn(irxn) = &
          sec_cont_cur_isotherm_rxn%Freundlich_n
        sec_cont_cur_isotherm_rxn => sec_cont_cur_isotherm_rxn%next
      endif

    enddo

    ! check for isotherm reaction using species with complexes
    found = PETSC_FALSE
    do icplx = 1, reaction%neqcplx
      do icomp = 1, reaction%eqcplxspecid(0,icplx)
        ipri_spec = reaction%eqcplxspecid(icomp,icplx)
        do irxn = 1, reaction%isotherm%neqkdrxn
          if (reaction%isotherm%eqkdspecid(irxn) == ipri_spec) then
            found = PETSC_TRUE
            option%io_buffer = 'Primary aqueous species "' // &
              trim(reaction%primary_species_names(ipri_spec)) // &
              '" is referenced in a sorption isotherm reaction and &
             &is associated with secondary aqueous complex "' // &
              trim(reaction%secondary_species_names(icplx)) // '".'
            call PrintMsg(option)
          endif
        enddo
      enddo
    enddo
    if (found) then
      option%io_buffer = 'Isotherm reactions can only be simulated for &
        &species without secondary aqueous complexes.  See comments above.'
      call PrintErrMsg(option)
    endif
  endif

  ! sandbox reactions
  call RSandboxSetup(reaction,option)
  call ReactionCLMRxnSetup(reaction,option)

end subroutine ReactionSetupKinetics

! ************************************************************************** !

subroutine ReactionSetupSpecificSpecies(reaction,option)
!
! Sets up integer pointers to specific species
!
! Author: Glenn Hammond
! Date: 01/29/24
!
  use Option_module
  use String_module

  implicit none

  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: ispec

  ! locate specific species
  reaction%species_idx => ReactionAuxCreateAqSpeciesIndex()
  do ispec = 1, reaction%naqcomp
    if (reaction%species_idx%h_ion_id == 0) then
      word = 'H+'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%h_ion_id = ispec
      endif
    endif
    if (reaction%species_idx%na_ion_id == 0) then
      word = 'Na+'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%na_ion_id = ispec
      endif
    endif
    if (reaction%species_idx%cl_ion_id == 0) then
      word = 'Cl-'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%cl_ion_id = ispec
      endif
    endif
    if (reaction%species_idx%co2_aq_id == 0) then
      word = 'CO2(aq)'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%co2_aq_id = ispec
      endif
    endif
    if (reaction%species_idx%tracer_aq_id == 0) then
      word = 'Tracer'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%tracer_aq_id = ispec
      endif
    endif
    if (reaction%species_idx%h2o_aq_id == 0) then
      word = 'H2O'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%h2o_aq_id = ispec
      endif
    endif
    if (reaction%species_idx%tracer_age_id == 0) then
      word = 'Tracer_Age'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%tracer_age_id = ispec
        reaction%calculate_tracer_age = PETSC_TRUE
      endif
    endif
    if (reaction%species_idx%water_age_id == 0) then
      word = 'Water_Age'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%water_age_id = ispec
        reaction%calculate_water_age = PETSC_TRUE
      endif
    endif
  enddo

  do ispec = 1, reaction%neqcplx
    if (reaction%species_idx%h_ion_id == 0) then
      word = 'H+'
      if (StringCompareIgnoreCase(reaction%secondary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%h_ion_id = -ispec
      endif
    endif
    if (reaction%species_idx%na_ion_id == 0) then
      word = 'Na+'
      if (StringCompareIgnoreCase(reaction%secondary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%na_ion_id = -ispec
      endif
    endif
    if (reaction%species_idx%cl_ion_id == 0) then
      word = 'Cl-'
      if (StringCompareIgnoreCase(reaction%secondary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%cl_ion_id = -ispec
      endif
    endif
    if (reaction%species_idx%co2_aq_id == 0) then
      word = 'CO2(aq)'
      if (StringCompareIgnoreCase(reaction%secondary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%co2_aq_id = -ispec
      endif
    endif
  enddo

  ! these are passive gases used in CONSTRAINTS
  do ispec = 1, reaction%gas%npassive_gas
    if (reaction%species_idx%o2_gas_id == 0) then
      word = 'O2(g)'
      if (StringCompareIgnoreCase(reaction%gas%passive_names(ispec), &
                                  word)) then
        reaction%species_idx%o2_gas_id = ispec
      endif
    endif
    if (reaction%species_idx%co2_gas_id == 0) then
      word = 'CO2(g)'
      if (StringCompareIgnoreCase(reaction%gas%passive_names(ispec), &
                                  word)) then
        reaction%species_idx%co2_gas_id = ispec
      endif
      word = 'CO2(g)*'
      if (StringCompareIgnoreCase(reaction%gas%passive_names(ispec), &
                                  word)) then
        reaction%species_idx%co2_gas_id = ispec
      endif

    endif

  enddo

end subroutine ReactionSetupSpecificSpecies

! ************************************************************************** !

subroutine ReactionSetupSpeciesSummary(reaction,option)
!
! Prints a summary of reaction species to .out file
!
! Author: Glenn Hammond
! Date: 01/29/24
!
  use Option_module
  use Reaction_Mineral_Aux_module
  use Reaction_Surface_Complexation_Aux_module

  implicit none

  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscInt :: i
  type(mineral_type), pointer :: mineral
  type(surface_complexation_type), pointer :: surface_complexation

  mineral => reaction%mineral
  surface_complexation => reaction%surface_complexation

90 format(80('-'))
100 format(/,2x,i4,2x,a)
110 format(100(/,14x,3(a20,2x)))

  if (OptionPrintToFile(option)) then
    write(option%fid_out,90)
    write(option%fid_out,100) reaction%naqcomp, 'Primary Species'
    write(option%fid_out,110) &
      (reaction%primary_species_names(i),i=1,reaction%naqcomp)

    write(option%fid_out,100) reaction%neqcplx, 'Secondary Complex Species'
    write(option%fid_out,110) &
      (reaction%secondary_species_names(i),i=1,reaction%neqcplx)

    write(option%fid_out,100) reaction%gas%nactive_gas, 'Active Gas Species'
    write(option%fid_out,110) (reaction%gas%active_names(i),i=1, &
                               reaction%gas%nactive_gas)

    write(option%fid_out,100) reaction%gas%npassive_gas, 'Passive Gas Species'
    write(option%fid_out,110) (reaction%gas%passive_names(i),i=1, &
                               reaction%gas%npassive_gas)

    write(option%fid_out,100) mineral%nmnrl, 'Reference Minerals'
    write(option%fid_out,110) (mineral%mineral_names(i),i=1,mineral%nmnrl)

    write(option%fid_out,100) mineral%nkinmnrl, 'Kinetic Mineral Reactions'
    write(option%fid_out,110) (mineral%kinmnrl_names(i),i=1,mineral%nkinmnrl)

    if (surface_complexation%nsrfcplxrxn > 0) then
      write(option%fid_out,100) surface_complexation%nsrfcplxrxn, &
        'Surface Complexation Reaction Sites'
      write(option%fid_out,110) &
        (surface_complexation%srfcplxrxn_site_names(i), &
         i=1,surface_complexation%nsrfcplxrxn)
      write(option%fid_out,100) surface_complexation%nsrfcplx, &
        'Surface Complexes'
      write(option%fid_out,110) (surface_complexation%srfcplx_names(i), &
        i=1,surface_complexation%nsrfcplx)
    endif

    write(option%fid_out,100) reaction%neqionxrxn, 'Ion Exchange Reactions'
    write(option%fid_out,100) reaction%neqionxcation, 'Ion Exchange Cations'
    write(option%fid_out,90)
  endif

#ifdef AMANZI_BGD
  ! output reaction in amanzi "bgd" formatted file
  if (OptionPrintToFile(option)) then
    string = trim(option%global_prefix) // '.bgd'
    open(unit=86,file=trim(string))

    write(86,'("# pflotran database preprocessing :")')
    call date_and_time(date=word,time=word2)
    write(86,'("#        date : ",a,"   ",a)') trim(word), trim(word2)
    write(86,'("#       input : ",a)') trim(option%input_filename)

    write(86,'(/,"<Primary Aqueous Species")')
    do icomp = 1, reaction%naqcomp
      write(86,'(a,x,3(" ; ",f6.2))') &
        trim(reaction%primary_species_names(icomp)), &
        reaction%primary_spec_a0(icomp), &
        reaction%primary_spec_Z(icomp), &
        reaction%primary_spec_molar_wt(icomp)
    enddo

    write(86,'(/,"<Aqueous Equilibrium Complexes")')
    do icplx = 1, reaction%neqcplx
      write(86,'(a," = ")',advance='no') &
        trim(reaction%secondary_species_names(icplx))
      if (reaction%eqcplxh2oid(icplx) > 0) then
        write(86,'(f6.2," H2O ")',advance='no') reaction%eqcplxh2ostoich(icplx)
      endif

      do i = 1,reaction%eqcplxspecid(0,icplx)
        temp_tin = reaction%eqcplxspecid(i,icplx)
        write(86,'(f6.2,x,a,x)',advance='no') reaction%eqcplxstoich(i,icplx), &
                                 trim(reaction%primary_species_names(temp_int))
      enddo
      write(86,'(4(" ; ",f10.5))') reaction%eqcplx_logK(icplx), &
                                   reaction%eqcplx_a0(icplx), &
                                   reaction%eqcplx_Z(icplx), &
                                   reaction%eqcplx_molar_wt(icplx)
    enddo

    write(86,'(/,"<General Kinetics")')
    do irxn = 1, reaction%ngeneral_rxn
      do i = 1, reaction%generalforwardspecid(0,irxn)
        temp_int = reaction%generalforwardspecid(i,irxn)
        write(86,'(f6.2,x,a)',advance='no') &
          reaction%generalforwardstoich(i,irxn), &
          trim(reaction%primary_species_names(temp_int))
        if (i /= reaction%generalforwardspecid(0,irxn)) then
          write(86,'(" + ")',advance='no')
        endif
      enddo
      write(86,'(" <-> ")',advance='no')
      do i = 1, reaction%generalbackwardspecid(0,irxn)
        temp_int = reaction%generalbackwardspecid(i,irxn)
        write(86,'(f6.2,x,a)',advance='no') &
          reaction%generalbackwardstoich(i,irxn), &
          trim(reaction%primary_species_names(temp_int))
        if (i /= reaction%generalbackwardspecid(0,irxn)) then
          write(86,'(" + ")',advance='no')
        endif
      enddo
      write(86,'(" ; ")',advance='no')
      do i = 1, reaction%generalforwardspecid(0,irxn)
        temp_int = reaction%generalforwardspecid(i,irxn)
        write(86,'(f6.2,x,a)',advance='no') &
          reaction%generalforwardstoich(i,irxn), &
          trim(reaction%primary_species_names(temp_int))
      enddo
      write(86,'(" ; ")',advance='no')
      write(86,'(1es13.5)',advance='no') reaction%general_kf(irxn)
      write(86,'(" ; ")',advance='no')
      do i = 1, reaction%generalbackwardspecid(0,irxn)
        temp_int = reaction%generalbackwardspecid(i,irxn)
        write(86,'(f6.2,x,a)',advance='no') &
          reaction%generalbackwardstoich(i,irxn), &
          trim(reaction%primary_species_names(temp_int))
      enddo
      write(86,'(" ; ")',advance='no')
      write(86,'(1es13.5)') reaction%general_kr(irxn)
      !write(86,'(" ; ")',advance='no')
      !write(86,'(f6.2)',advance='no') reaction%generalh2ostoich(irxn)
    enddo

    write(86,'(/,"<Minerals")')

    do imnrl = 1, mineral%nkinmnrl
      write(86,'(a," = ")',advance='no') trim(mineral%kinmnrl_names(imnrl))
      if (mineral%kinmnrlh2oid(imnrl) > 0) then
        write(86,'(f6.2," H2O ")',advance='no') mineral%kinmnrlh2ostoich(imnrl)
      endif
      do i = 1, mineral%kinmnrlspecid(0,imnrl)
        temp_tin = mineral%kinmnrlspecid(i,imnrl)
        write(86,'(f6.2,x,a,x)',advance='no') mineral%kinmnrlstoich(i,imnrl), &
                                 trim(reaction%primary_species_names(temp_int))
      enddo
      !molar volume has been converted to m^3/mol!
      write(86,'(4(" ; ",1es13.5))') mineral%kinmnrl_logK(imnrl), &
                                     mineral%kinmnrl_molar_wt(imnrl), &
                                     mineral%kinmnrl_molar_vol(imnrl)*1.d6, 1.0
    enddo

    write(86,'(/,"<Mineral Kinetics")')
    do imnrl = 1, mineral%nkinmnrl
      write(86,'(a," ; TST ; log10_rate_constant ")',advance='no') &
        trim(mineral%kinmnrl_names(imnrl))
      write(86,'(1es13.5," moles/cm^2/sec ")',advance='no') &
        log10(mineral%kinmnrl_rate_constant(imnrl))
      if (mineral%kinmnrl_num_prefactors(imnrl) /= 0) then
        write(86,'(" ; ")',advance='no')
        do i = 1, mineral%kinmnrl_num_prefactors(imnrl)
          ! number of prefactor species stored in
          ! kinmnrl_prefactor_id(0,i,imnrl)
          do j = 1, mineral%kinmnrl_prefactor_id(0,i,imnrl)
            temp_int = mineral%kinmnrl_prefactor_id(j,i,imnrl)
            if (temp_int > 0) then
              write(86,'(a)',advance='no') &
                trim(reaction%primary_species_names(temp_int))
            else
              write(86,'(a)',advance='no') &
                trim(reaction%secondary_species_names(-temp_int))
            endif
            write(86,'(x,1es13.5,x)',advance='no') &
              mineral%kinmnrl_pref_alpha(j,i,imnrl)
          enddo
        enddo
      endif
      write(86,*)
    enddo

    write(86,'(/,"<Ion Exchange Sites")')
    do irxn = 1, reaction%neqionxrxn
      write(86,'("X- ; -1.0 ; ",a)') &
        trim(reaction%ion_exchange_rxn_list%mineral_name)
    enddo

    write(86,'(/,"<Ion Exchange Complexes")')
    do irxn = 1, reaction%neqionxrxn
      do i = 1, reaction%neqionxcation
        temp_int = reaction%eqionx_rxn_cationid(i,irxn)
        write(86,'(a,"X = 1.0 ",a)',advance='no') &
          trim(reaction%primary_species_names(temp_int)), &
          trim(reaction%primary_species_names(temp_int))
        write(86,'(f6.2," X- ")',advance='no') reaction%primary_spec_Z(temp_int)
        write(86,'(" ; ",1es13.5)') reaction%eqionx_rxn_k(i,irxn)
      enddo
    enddo

    write(86,'(/,"<Surface Complex Sites")')
    do ieqrxn = 1, surface_complexation%neqsrfcplxrxn
      irxn = surface_complexation%eqsrfcplxrxn_to_srfcplxrxn(ieqrxn)
      write(86,'(a, " ; ")',advance='no') &
        trim(surface_complexation%srfcplxrxn_site_names(irxn))
      write(86,'(1es13.5)') surface_complexation%srfcplxrxn_site_density(irxn)
    enddo

    write(86,'(/,"<Surface Complexes")')
    do ieqrxn = 1, surface_complexation%neqsrfcplxrxn
      irxn = surface_complexation%eqsrfcplxrxn_to_srfcplxrxn(ieqrxn)
      do i = 1, surface_complexation%srfcplxrxn_to_complex(0,irxn)
        icplx = surface_complexation%srfcplxrxn_to_complex(i,irxn)
        write(86,'(a, " = ")',advance='no') &
          trim(surface_complexation%srfcplx_names(icplx))
        write(86,'(f6.2,x,a)',advance='no') &
          surface_complexation%srfcplx_free_site_stoich(icplx), &
          trim(surface_complexation%srfcplxrxn_site_names(irxn))

        if (surface_complexation%srfcplxh2oid(icplx) > 0) then
          write(86,'(f6.2," H2O ")',advance='no') &
            surface_complexation%srfcplxh2ostoich(icplx)
        endif
        do j = 1, surface_complexation%srfcplxspecid(0,icplx)
          temp_int = surface_complexation%srfcplxspecid(j,icplx)
          write(86,'(f6.2,x,a)',advance='no') &
            surface_complexation%srfcplxstoich(j,icplx), &
            trim(reaction%primary_species_names(temp_int))
        enddo
        write(86,'(" ; ",1es13.5," ; ",f6.2)') &
          surface_complexation%srfcplx_logK(icplx), &
          surface_complexation%srfcplx_Z(icplx)

      enddo
    enddo

    write(86,'(/,"<Isotherms")')
    do irxn = 1, reaction%isotherm%neqkdrxn
       write(86,'(a," ; ")',advance='no') &
         trim(reaction%primary_species_names(reaction%isotherm%eqkdspecid(irxn)))
      select case (reaction%isotherm%eqisothermtype(irxn))
        case(SORPTION_LINEAR)
           write(86,'("linear ; ",es13.5)',advance='no') &
             reaction%isotherm%isotherm_rxn%eqisothermcoeff(irxn)
           write(86,'()')
        case(SORPTION_LANGMUIR)
           write(86,'("langmuir ; ",es13.5)',advance='no') &
             reaction%isotherm%isotherm_rxn%eqisothermcoeff(irxn)
           write(86,'(es13.5)') reaction%isotherm%isotherm_rxn%eqisothermlangmuirb(irxn)
        case(SORPTION_FREUNDLICH)
           write(86,'("freundlich ; ",es13.5)',advance='no') &
             reaction%isotherm%isotherm_rxn%eqisothermcoeff(irxn)
           write(86,'(es13.5)') reaction%isotherm%isotherm_rxn%eqisothermfreundlichn(irxn)
      end select
    enddo

    close(86)
  endif
#endif
! AMANZI_BGD

#if 0
  ! output for ASCEM reactions
  if (OptionPrintToFile(option)) then
    open(unit=86,file='reaction.dat')
    write(86,'(10i4)') reaction%naqcomp, reaction%neqcplx, &
                       reaction%ngeneral_rxn, &
                       reaction%neqsrfcplxrxn, mineral%nkinmnrl
    do icomp = 1, reaction%naqcomp
      write(86,'(a12,f6.2,f6.2)') reaction%primary_species_names(icomp), &
                                  reaction%primary_spec_Z(icomp), &
                                  reaction%primary_spec_a0(icomp)
    enddo
    do icplx = 1, reaction%neqcplx
      write(86,'(a32,f6.2,f6.2)') reaction%secondary_species_names(icplx), &
                                  reaction%eqcplx_Z(icplx), &
                                  reaction%eqcplx_a0(icplx)
      write(86,'(40i4)') reaction%eqcplxspecid(:,icplx)
      write(86,'(40f6.2)') reaction%eqcplxstoich(:,icplx)
      write(86,'(i4)') reaction%eqcplxh2oid(icplx)
      write(86,'(f6.2)') reaction%eqcplxh2ostoich(icplx)
      write(86,'(1es13.5)') reaction%eqcplx_logK(icplx)
    enddo
    do irxn = 1, reaction%ngeneral_rxn
      write(86,'(40i4)') reaction%generalspecid(:,irxn)
      write(86,'(40f6.2)') reaction%generalstoich(:,irxn)
      write(86,'(40i4)') reaction%generalforwardspecid(:,irxn)
      write(86,'(40f6.2)') reaction%generalforwardstoich(:,irxn)
      write(86,'(40i4)') reaction%generalbackwardspecid(:,irxn)
      write(86,'(40f6.2)') reaction%generalbackwardstoich(:,irxn)
      write(86,'(f6.2)') reaction%generalh2ostoich(irxn)
      write(86,'(1es13.5)') reaction%general_kf(irxn)
      write(86,'(1es13.5)') reaction%general_kr(irxn)
    enddo
    do irxn = 1, reaction%neqsrfcplxrxn
      write(86,'(a32)')reaction%eqsrfcplx_site_names(irxn)
      write(86,'(1es13.5)') reaction%eqsrfcplx_rxn_site_density(irxn)
      write(86,'(i4)') reaction%srfcplxrxn_to_complex(0,irxn) ! # complexes
      do i = 1, reaction%srfcplxrxn_to_complex(0,irxn)
        icplx = reaction%srfcplxrxn_to_complex(i,irxn)
        write(86,'(a32,f6.2)') reaction%eqsrfcplx_names(icplx), &
                               reaction%eqsrfcplx_Z(icplx)
        write(86,'(40i4)') reaction%srfcplxspecid(:,icplx)
        write(86,'(40f6.2)') reaction%eqsrfcplxstoich(:,icplx)
        write(86,'(i4)') reaction%eqsrfcplxh2oid(icplx)
        write(86,'(f6.2)') reaction%eqsrfcplxh2ostoich(icplx)
        write(86,'(i4)') reaction%eqsrfcplx_free_site_id(icplx)
        write(86,'(f6.2)') reaction%eqsrfcplx_free_site_stoich(icplx)
        write(86,'(1es13.5)') reaction%eqsrfcplx_logK(icplx)

      enddo
    enddo
    do imnrl = 1, mineral%nkinmnrl
      write(86,'(a32)') mineral%kinmnrl_names(imnrl)
      write(86,'(40i4)') mineral%kinmnrlspecid(:,imnrl)
      write(86,'(40f6.2)') mineral%kinmnrlstoich(:,imnrl)
      write(86,'(i4)') mineral%kinmnrlh2oid(imnrl)
      write(86,'(f6.2)') mineral%kinmnrlh2ostoich(imnrl)
      write(86,'(1es13.5)') mineral%kinmnrl_logK(imnrl)
      write(86,'(1es13.5)') mineral%kinmnrl_molar_vol(imnrl)
      write(86,'(1es13.5)') mineral%kinmnrl_molar_wt(imnrl)
      write(86,'(1es13.5)') mineral%kinmnrl_rate_constant(1,imnrl)
      write(86,'(1es13.5)') 1.d0 ! specific surface area 1 cm^2 / cm^3
    enddo
        close(86)
  endif
#endif

end subroutine ReactionSetupSpeciesSummary

! ************************************************************************** !

subroutine ReactionSetupPrimaryPrint(reaction,option)
!
! Prints a summary of reaction species to .out file
!
! Author: Glenn Hammond
! Date: 01/29/24
!
  use Option_module
  use Reaction_Mineral_Aux_module
  use Reaction_Surface_Complexation_Aux_module

  implicit none

  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  type(aq_species_type), pointer :: cur_pri_aq_spec
  PetscInt :: ispec

  ! primary species
  allocate(reaction%primary_species_print(reaction%naqcomp))
  reaction%primary_species_print = PETSC_FALSE

  allocate(reaction%kd_print(reaction%naqcomp))
  reaction%kd_print = PETSC_FALSE
  if (reaction%nsorb > 0) then
    allocate(reaction%total_sorb_print(reaction%naqcomp))
    reaction%total_sorb_print = PETSC_FALSE
  endif

  cur_pri_aq_spec => reaction%primary_species_list
  ispec = 1
  do
    if (.not.associated(cur_pri_aq_spec)) exit
    reaction%primary_species_print(ispec) = cur_pri_aq_spec%print_me .or. &
                                            reaction%print_all_primary_species
    reaction%kd_print(ispec) = (cur_pri_aq_spec%print_me .or. &
                                reaction%print_all_primary_species) .and. &
                                reaction%print_kd
    if (reaction%nsorb > 0) then
      reaction%total_sorb_print(ispec) = (cur_pri_aq_spec%print_me .or. &
                                  reaction%print_all_primary_species) .and. &
                                  reaction%print_total_sorb
    endif
    ispec = ispec + 1
    cur_pri_aq_spec => cur_pri_aq_spec%next
  enddo
  nullify(cur_pri_aq_spec)

end subroutine ReactionSetupPrimaryPrint

end module Reaction_Setup_module
