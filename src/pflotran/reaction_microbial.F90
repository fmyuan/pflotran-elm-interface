module Reaction_Microbial_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Microbial_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  public :: MicrobialRead, &
            RMicrobial

contains

! ************************************************************************** !

subroutine MicrobialRead(microbial,input,option)
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
  
  type(microbial_type) :: microbial
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  type(microbial_rxn_type), pointer :: microbial_rxn, cur_microbial_rxn
  type(monod_type), pointer :: monod, prev_monod
  type(microbial_biomass_type), pointer :: microbial_biomass
  type(inhibition_type), pointer :: inhibition, prev_inhibition
  PetscInt :: temp_int
  
  microbial%nrxn = microbial%nrxn + 1
        
  microbial_rxn => MicrobialRxnCreate()
  nullify(prev_monod)
  nullify(prev_inhibition)
  nullify(microbial_biomass)
  temp_int = UNINITIALIZED_INTEGER
  call InputPushBlock(input,option)
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,MICROBIAL_REACTION')
    call StringToUpper(word)   

    select case(trim(word))
      case('REACTION')
        ! remainder of string should be the reaction equation
        microbial_rxn%reaction = trim(adjustl(input%buf))
        ! set flag for error message
        if (len_trim(microbial_rxn%reaction) < 2) input%ierr = 1
        call InputErrorMsg(input,option,'reaction string', &
                            'CHEMISTRY,MICROBIAL_REACTION,REACTION')     
      case('CONCENTRATION_UNITS')
        call InputReadWord(input,option,word,PETSC_TRUE)  
        call InputErrorMsg(input,option,'concentration units', &
                           'CHEMISTRY,MICROBIAL_REACTION') 
        call StringToUpper(word)   
        select case(word)
          case('MOLALITY')
            temp_int = MICROBIAL_MOLALITY
          case('ACTIVITY')
            temp_int = MICROBIAL_ACTIVITY
          case('MOLARITY')
            temp_int = MICROBIAL_MOLARITY
          case default
            option%io_buffer = 'Unrecognized concentration units in a &
              &microbial reaction.'
            call PrintErrMsg(option)
        end select
        if (Initialized(microbial%concentration_units)) then
          if (temp_int /= microbial%concentration_units) then
            option%io_buffer = 'Concentration units must be consistent for &
              &all microbial reactions. Default is MOLARITY if not specified.'
            call PrintErrMsg(option)
          endif
        else
          microbial%concentration_units = temp_int
        endif
      case('RATE_CONSTANT')
        call InputReadDouble(input,option,microbial_rxn%rate_constant)  
        call InputErrorMsg(input,option,'rate constant', &
                           'CHEMISTRY,MICROBIAL_REACTION') 
      case('ACTIVATION_ENERGY')
        call InputReadDouble(input,option,microbial_rxn%activation_energy)  
        call InputErrorMsg(input,option,'activation energy', &
                           'CHEMISTRY,MICROBIAL_REACTION') 
        call InputReadAndConvertUnits(input,microbial_rxn%activation_energy, &
                     'J/mol', &
                     'CHEMISTRY,MICROBIAL_REACTION,ACTIVATION_ENERGY',option)
      case('MONOD')
        monod => MicrobialMonodCreate()
        call InputPushBlock(input,option)
        do 
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,MICROBIAL_REACTION,MONOD')
          call StringToUpper(word)   
          select case(trim(word))
            case('SPECIES_NAME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'species name', &
                           'CHEMISTRY,MICROBIAL_REACTION,MONOD')
              monod%species_name = word
            case('HALF_SATURATION_CONSTANT')
              call InputReadDouble(input,option,monod%half_saturation_constant)
              call InputErrorMsg(input,option,'half saturation constant', &
                                 'CHEMISTRY,MICROBIAL_REACTION,MONOD')
            case('THRESHOLD_CONCENTRATION')
              call InputReadDouble(input,option,monod%threshold_concentration)  
              call InputErrorMsg(input,option,'threshold concdntration', &
                                 'CHEMISTRY,MICROBIAL_REACTION,MONOD')
            case default
              call InputKeywordUnrecognized(input,word, &
                                     'CHEMISTRY,MICROBIAL_REACTION,MONOD', &
                                     option)
          end select
        enddo
        call InputPopBlock(input,option)
        ! append to list
        if (.not.associated(microbial_rxn%monod)) then
          microbial_rxn%monod => monod
        else
          prev_monod%next => monod
        endif
        prev_monod => monod
        nullify(monod)
      case('INHIBITION')
        inhibition => MicrobialInhibitionCreate()
        call InputPushBlock(input,option)
        do 
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,MICROBIAL_REACTION,INHIBITION')
          call StringToUpper(word)   
          select case(trim(word))
            case('SPECIES_NAME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'species name', &
                                 'CHEMISTRY,MICROBIAL_REACTION,INHIBITION')
              inhibition%species_name = word
            case('TYPE')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'inhibition type', &
                                 'CHEMISTRY,MICROBIAL_REACTION,INHIBITION')
              call StringToUpper(word)   
              select case(word)
                case('MONOD')
                  inhibition%itype = INHIBITION_MONOD
                case('INVERSE_MONOD')
                  inhibition%itype = INHIBITION_INVERSE_MONOD
                case('THRESHOLD')
                  inhibition%itype = INHIBITION_THRESHOLD
                  call InputReadDouble(input,option, &
                                       inhibition%inhibition_constant2)  
                  call InputErrorMsg(input,option,'scaling factor', &
                                     'CHEMISTRY,MICROBIAL_REACTION,&
                                     &INHIBITION,THRESHOLD_INHIBITION')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'CHEMISTRY,MICROBIAL_REACTION,INHIBITION,TYPE',option)
              end select
            case('INHIBITION_CONSTANT')
              call InputReadDouble(input,option,inhibition%inhibition_constant)  
              call InputErrorMsg(input,option,'inhibition constant', &
                                 'CHEMISTRY,MICROBIAL_REACTION,INHIBITION')
            case default
              call InputKeywordUnrecognized(input,word, &
                      'CHEMISTRY,MICROBIAL_REACTION,INHIBITION',option)
          end select        
        enddo
        call InputPopBlock(input,option)
        if (len_trim(inhibition%species_name) < 2 .or. &
            inhibition%itype == 0 .or. &
            Uninitialized(inhibition%inhibition_constant)) then
          option%io_buffer = 'A SPECIES_NAME, TYPE, and INHIBITION_CON' // &
            'STANT must be defined for INHIBITION in MICROBIAL_REACTION ' // &
            'with REACTION "' // trim(microbial_rxn%reaction) // '".'
          call PrintErrMsg(option)
        endif
        ! append to list
        if (.not.associated(microbial_rxn%inhibition)) then
          microbial_rxn%inhibition => inhibition
        else
          prev_inhibition%next => inhibition
        endif
        prev_inhibition => inhibition
        nullify(inhibition)
      case('BIOMASS')
        if (associated(microbial_biomass)) then
          call MicrobialBiomassDestroy(microbial_biomass)
        endif
        microbial_biomass => MicrobialBiomassCreate()
        call InputPushBlock(input,option)
        do 
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,MICROBIAL_REACTION,BIOMASS')
          call StringToUpper(word)   
          select case(trim(word))
            case('SPECIES_NAME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'species name', &
                                 'CHEMISTRY,MICROBIAL_REACTION,BIOMASS')
              microbial_biomass%species_name = word
            case('YIELD')
              call InputReadDouble(input,option,microbial_biomass%yield)
              call InputErrorMsg(input,option,'yield', &
                                 'CHEMISTRY,MICROBIAL_REACTION,BIOMASS')
            case default
              call InputKeywordUnrecognized(input,word, &
                                      'CHEMISTRY,MICROBIAL_REACTION,BIOMASS', &
                                            option)
          end select
        enddo
        call InputPopBlock(input,option)
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'CHEMISTRY,MICROBIAL_REACTION', &
                                      option)
    end select
  enddo
  call InputPopBlock(input,option)
  
  ! add linkage to biomass if exists
  microbial_rxn%biomass => microbial_biomass
  
  ! add to list
  if (.not.associated(microbial%microbial_rxn_list)) then
    microbial%microbial_rxn_list => microbial_rxn
    microbial_rxn%id = 1
  else
    cur_microbial_rxn => microbial%microbial_rxn_list
    do
      if (.not.associated(cur_microbial_rxn%next)) then
        cur_microbial_rxn%next => microbial_rxn
        microbial_rxn%id = cur_microbial_rxn%id + 1
        exit
      endif
      cur_microbial_rxn => cur_microbial_rxn%next
    enddo
  endif
  nullify(microbial_rxn)

  if (Uninitialized(microbial%concentration_units)) then
    microbial%concentration_units = MICROBIAL_MOLARITY
  endif

end subroutine MicrobialRead

! ************************************************************************** !

subroutine RMicrobial(Res,Jac,compute_derivative,rt_auxvar, &
                      global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Computes the microbial reaction
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/31/12
  ! 

  use Option_module, only : option_type
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type
  use Global_Aux_module, only : global_auxvar_type
  use Material_Aux_class, only : material_auxvar_type
  use Reaction_Aux_module, only : reaction_rt_type
  use Reaction_Immobile_Aux_module, only : immobile_type
  
  implicit none
  
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  
  PetscInt, parameter :: iphase = 1
  PetscInt :: irxn, i, ii, icomp, jcomp, ncomp
  PetscInt :: imonod, iinhibition, ibiomass
  PetscReal :: Im
  PetscReal :: rate_constant
  PetscReal :: concentration(reaction%naqcomp)
  PetscReal :: dconcentration_dmolal(reaction%naqcomp)
  PetscReal :: conc
  PetscReal :: dconc_dmolal
  PetscReal :: monod(10)
  PetscReal :: inhibition(10)
  PetscReal :: biomass_conc, yield
  PetscInt :: immobile_id
  PetscReal :: denominator, dR_dX, dX_dc, dR_dc, dR_dbiomass
  PetscReal :: tempreal
  PetscReal :: L_water
  type(microbial_type), pointer :: microbial
  type(immobile_type), pointer :: immobile
  
  microbial => reaction%microbial
  immobile => reaction%immobile

  select case(microbial%concentration_units)
    case(MICROBIAL_MOLALITY)
      dconcentration_dmolal = 1.d0
    case(MICROBIAL_ACTIVITY)
      dconcentration_dmolal = rt_auxvar%pri_act_coef
    case(MICROBIAL_MOLARITY)
      dconcentration_dmolal = global_auxvar%den_kg(iphase)*1.d-3
  end select
  concentration = rt_auxvar%pri_molal*dconcentration_dmolal

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  ! units:
  ! Residual: mol/sec
  ! Jacobian: (mol/sec)*(kg water/mol) = kg water/sec

  do irxn = 1, microbial%nrxn
  
    ! units:
    !   without biomass: mol/L-sec
    !   with biomass: mol/L-sec * (m^3 bulk / mol biomass)

    ncomp = microbial%specid(0,irxn)
    rate_constant = microbial%rate_constant(irxn)
    Im = rate_constant
    if (associated(microbial%activation_energy)) then
      ! ideal gas constant units: J/mol-K
      Im = Im * exp(microbial%activation_energy(irxn)/IDEAL_GAS_CONSTANT* &
                    (1.d0/298.15d0-1.d0/(global_auxvar%temp+273.15d0)))
    endif
    yield = 0.d0
    biomass_conc = 0.d0

    ! monod expressions
    do ii = 1, microbial%monodid(0,irxn)
      imonod = microbial%monodid(ii,irxn)
      icomp = microbial%monod_specid(imonod)
      conc = concentration(icomp)
      monod(ii) = (conc - microbial%monod_Cth(imonod)) / &
                  (microbial%monod_K(imonod) + conc - &
                   microbial%monod_Cth(imonod))
      Im = Im*monod(ii)
    enddo

    ! inhibition expressions
    do ii = 1, microbial%inhibitionid(0,irxn)
      iinhibition = microbial%inhibitionid(ii,irxn)
      icomp = microbial%inhibition_specid(iinhibition)
      conc = concentration(icomp)
      select case(microbial%inhibition_type(iinhibition))
        case(INHIBITION_MONOD)
          inhibition(ii) = microbial%inhibition_C(iinhibition) / &
                          (microbial%inhibition_C(iinhibition) + conc)
        case(INHIBITION_INVERSE_MONOD)
          inhibition(ii) = conc / &
                          (microbial%inhibition_C(iinhibition) + conc)
        case(INHIBITION_THRESHOLD)
          ! if microbial%inhibition_C2 is negative, inhibition kicks in 
          ! above the concentration
          inhibition(ii) = 0.5d0 + &
                           atan((conc - &
                                 microbial%inhibition_C(iinhibition)) * &
                                microbial%inhibition_C2(iinhibition)) / PI
      end select        
      Im = Im*inhibition(ii)
    enddo
    
    ! biomass term
    ibiomass = microbial%biomassid(irxn)
    if (ibiomass > 0) then
      immobile_id = reaction%offset_immobile + ibiomass
      biomass_conc = rt_auxvar%immobile(ibiomass)
      yield = microbial%biomass_yield(irxn)
      ! Im units (before): mol/mol biomass-sec
      Im = Im*biomass_conc*material_auxvar%volume
    else
      ! Im units (before): mol/L-sec
      Im = Im * L_water
    endif
    ! Im units (after): mol/sec
    
    do i = 1, ncomp
      icomp = microbial%specid(i,irxn)
      Res(icomp) = Res(icomp) - microbial%stoich(i,irxn)*Im
    enddo

    if (ibiomass > 0) then
      Res(immobile_id) = Res(immobile_id) - yield*Im
    endif
    
    if (.not. compute_derivative) cycle
    
    ! monod expressions
    do ii = 1, microbial%monodid(0,irxn)
      imonod = microbial%monodid(ii,irxn)
      jcomp = microbial%monod_specid(imonod)
      conc = concentration(jcomp)
      dconc_dmolal = dconcentration_dmolal(jcomp)
        
      dR_dX = Im / monod(ii)
        
      denominator = microbial%monod_K(imonod) + conc - &
                    microbial%monod_Cth(imonod)
        
      dX_dc = dconc_dmolal / denominator - &
              dconc_dmolal * (conc - microbial%monod_Cth(imonod)) / &
              (denominator*denominator)
        
      dR_dc = -1.d0*dR_dX*dX_dc
      do i = 1, ncomp
        icomp = microbial%specid(i,irxn)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                            microbial%stoich(i,irxn)*dR_dc
      enddo
      if (ibiomass > 0) then
        Jac(immobile_id,jcomp) = Jac(immobile_id,jcomp) + yield*dR_dc      
      endif      
    enddo

    ! inhibition expressions
    do ii = 1, microbial%inhibitionid(0,irxn)
      iinhibition = microbial%inhibitionid(ii,irxn)
      jcomp = microbial%inhibition_specid(iinhibition)
      conc = concentration(jcomp)
      dconc_dmolal = dconcentration_dmolal(jcomp)

      dR_dX = Im / inhibition(ii)
        
      select case(microbial%inhibition_type(iinhibition))
        case(INHIBITION_MONOD)
          denominator = microbial%inhibition_C(iinhibition) + conc
          dX_dc = -1.d0 * dconc_dmolal *microbial%inhibition_C(iinhibition) / &
                  (denominator*denominator)
        case(INHIBITION_INVERSE_MONOD)
          denominator = microbial%inhibition_C(iinhibition) + conc
          dX_dc = dconc_dmolal / denominator - &
                  dconc_dmolal * conc / &
                  (denominator*denominator)
        case(INHIBITION_THRESHOLD)
          ! derivative of atan(X) = 1 / (1 + X^2) dX
          tempreal = (conc - microbial%inhibition_C(iinhibition)) * &
                     microbial%inhibition_C2(iinhibition)
          dX_dc = (microbial%inhibition_C2(iinhibition) * dconc_dmolal / &
                   (1.d0 + tempreal*tempreal)) / PI
      end select        
      
      dR_dc = -1.d0*dR_dX*dX_dc
      do i = 1, ncomp
        icomp = microbial%specid(i,irxn)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                            microbial%stoich(i,irxn)*dR_dc
      enddo
      if (ibiomass > 0) then
        Jac(immobile_id,jcomp) = Jac(immobile_id,jcomp) + yield*dR_dc        
      endif
    enddo

    ! biomass expression
    if (ibiomass > 0) then
      dR_dbiomass = -1.d0*Im / biomass_conc
      do i = 1, ncomp
        icomp = microbial%specid(i,irxn)
        ! units = (mol/sec)*(m^3/mol) = m^3/sec
        Jac(icomp,immobile_id) = Jac(icomp,immobile_id) + &
                            microbial%stoich(i,irxn)*dR_dbiomass
      enddo
      Jac(immobile_id,immobile_id) = Jac(immobile_id,immobile_id) + &
        yield*dR_dbiomass
    endif
    
  enddo
    
end subroutine RMicrobial

end module Reaction_Microbial_module
