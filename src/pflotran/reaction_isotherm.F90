module Reaction_Isotherm_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  
  use PFLOTRAN_Constants_module
  use Reaction_Isotherm_Aux_module
  use Reactive_Transport_Aux_module

  implicit none

  public :: IsothermRead, &
            RTotalSorbKD

contains

! ************************************************************************** !

subroutine IsothermRead(isotherm,input,option)
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

  type(isotherm_type) :: isotherm
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word, word2
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXWORDLENGTH) :: kd_units
  character(len=MAXWORDLENGTH) :: multi_kd_units
  PetscInt :: ikd_units, imulti_kd_units
  
  type(isotherm_link_type), pointer :: isotherm_rxn, prev_isotherm_rxn
  type(isotherm_link_type), pointer :: sec_cont_isotherm_rxn, &
                                       sec_cont_prev_isotherm_rxn

  PetscErrorCode :: ierr

  if (associated(isotherm%isotherm_list)) then
    option%io_buffer = 'SORPTION ISOTHERM_REACTIONS must be entered in a &
      &single block.'
    call PrintErrMsg(option)
  endif

  nullify(prev_isotherm_rxn)
  nullify(sec_cont_isotherm_rxn)
  nullify(sec_cont_prev_isotherm_rxn)
 
  ! set default units
  isotherm%ikd_units = UNINITIALIZED_INTEGER

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    ! first string is species name
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'species name', &
                       'CHEMISTRY,ISOTHERM_REACTIONS')
    kd_units = ''
    multi_kd_units = ''
    isotherm%neqkdrxn = isotherm%neqkdrxn + 1
    isotherm_rxn => IsothermLinkCreate()
    isotherm_rxn%species_name = trim(word)
    if (option%use_mc) then
      sec_cont_isotherm_rxn => IsothermLinkCreate()
      sec_cont_isotherm_rxn%species_name = isotherm_rxn%species_name
      sec_cont_isotherm_rxn%Kd = UNINITIALIZED_DOUBLE
    endif

    call InputPushBlock(input,option)
    do 
      call InputReadPflotranString(input,option)
      if (InputError(input)) exit
      if (InputCheckExit(input,option)) exit

      call InputReadCard(input,option,word)
      call InputErrorMsg(input,option,'keyword', &
                         'CHEMISTRY,ISOTHERM_REACTIONS')
      call StringToUpper(word)
                  
      ! default type is linear
      isotherm_rxn%itype = SORPTION_LINEAR
      select case(trim(word))
        case('TYPE')
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'type', &
                             'CHEMISTRY,ISOTHERM_REACTIONS')
          select case(word)
            case('LINEAR')
              isotherm_rxn%itype = SORPTION_LINEAR
            case('LANGMUIR')
              isotherm_rxn%itype = SORPTION_LANGMUIR
            case('FREUNDLICH')
              isotherm_rxn%itype = SORPTION_FREUNDLICH
            case default
              call InputKeywordUnrecognized(input,word, &
                    'CHEMISTRY,SORPTION,ISOTHERM_REACTIONS,TYPE', &
                    option)
          end select
          if (option%use_mc) then
            sec_cont_isotherm_rxn%itype = isotherm_rxn%itype
          endif
        case('KD')
          call InputKeywordDeprecated('KD', &
                               'DISTRIBUTION_COEFFICIENT',option)
        case('DISTRIBUTION_COEFFICIENT')
          call InputReadDouble(input,option,isotherm_rxn%Kd)
          call InputErrorMsg(input,option, &
                             'DISTRIBUTION_COEFFICIENT', &
                             'CHEMISTRY,ISOTHERM_REACTIONS')
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (input%ierr == 0) kd_units = trim(word)
        ! S.Karra, 02/20/2014
        case('SEC_CONT_DISTRIBUTION_COEFFICIENT', &
              'SEC_CONT_KD')
          if (.not.associated(sec_cont_isotherm_rxn)) then
            option%io_buffer = 'Make sure MULTIPLE_CONTINUUM &
              &keyword is set, SECONDARY_CONTINUUM_KD.'
            call PrintErrMsg(option)
          endif
          call InputReadDouble(input,option,sec_cont_isotherm_rxn%Kd)
          call InputErrorMsg(input,option, &
            'SECONDARY_CONTINUUM_DISTRIBUTION_COEFFICIENT', &
            'CHEMISTRY,ISOTHERM_REACTIONS')
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (input%ierr == 0) multi_kd_units = trim(word)
        case('LANGMUIR_B')
          call InputReadDouble(input,option,isotherm_rxn%Langmuir_B)
          call InputErrorMsg(input,option,'Langmuir_B', &
                                         'CHEMISTRY,ISOTHERM_REACTIONS')
          isotherm_rxn%itype = SORPTION_LANGMUIR
        case('FREUNDLICH_N')
          call InputReadDouble(input,option,isotherm_rxn%Freundlich_N)
          call InputErrorMsg(input,option,'Freundlich_N', &
                             'CHEMISTRY,ISOTHERM_REACTIONS')
          isotherm_rxn%itype = SORPTION_FREUNDLICH
        case('KD_MINERAL_NAME')
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'KD_MINERAL_NAME', &
                             'ISOTHERM_REACTIONS,KD_MINERAL_NAME')
          isotherm_rxn%kd_mineral_name = word                      
        case default
          call InputKeywordUnrecognized(input,word, &
                 'CHEMISTRY,SORPTION,ISOTHERM_REACTIONS',option)
      end select
    enddo
    call InputPopBlock(input,option)

    if (associated(sec_cont_isotherm_rxn)) then
      if (UnInitialized(sec_cont_isotherm_rxn%Kd)) then
        option%io_buffer = 'SEC_CONT_DISTRIBUTION_COEFFICIENT must be &
          &specified for all ISOTHERM_REACTIONs when using &
          &MULTIPLE_CONTINUUM.'
        call PrintErrMsg(option)
      endif
    endif

    call IsothermConvertKDUnits(isotherm_rxn%Kd,kd_units,ikd_units,option)
    if (associated(sec_cont_isotherm_rxn)) then
      call IsothermConvertKDUnits(sec_cont_isotherm_rxn%Kd,multi_kd_units, &
                                  imulti_kd_units,option)
      if (ikd_units /= imulti_kd_units) then
        if (len_trim(kd_units) == 0) kd_units = 'default'
        if (len_trim(multi_kd_units) == 0) multi_kd_units = 'default'
        option%io_buffer = 'kd unit must match between primary &
                            &and secondary continuum: ' // &
                            trim(kd_units) // ' vs ' // trim(multi_kd_units)
        call PrintErrMsg(option)
      endif
    endif
    if (Initialized(isotherm%ikd_units) .and. &
        isotherm%ikd_units /= ikd_units) then
      if (ikd_units == KD_UNIT_KG_M3_BULK) then
        word = 'kg/m^3'
      else
        word = 'L/kg'
      endif
      if (isotherm%ikd_units == KD_UNIT_KG_M3_BULK) then
        word2 = 'kg/m^3'
      else
        word2 = 'L/kg'
      endif
      option%io_buffer = 'Isotherm KD units must be consistent throughout &
        &input file: ' // trim(word) // ' vs ' // trim(word2)
    else
      isotherm%ikd_units = ikd_units
    endif
      
    ! add to list
    if (.not.associated(isotherm%isotherm_list)) then
      isotherm%isotherm_list => isotherm_rxn
      isotherm_rxn%id = 1
    endif
    if (associated(prev_isotherm_rxn)) then
      prev_isotherm_rxn%next => isotherm_rxn
      isotherm_rxn%id = prev_isotherm_rxn%id + 1
    endif
    prev_isotherm_rxn => isotherm_rxn
    nullify(isotherm_rxn)               

    if (associated(sec_cont_isotherm_rxn)) then
      ! add to list
      if (.not.associated(isotherm%multicontinuum_isotherm_list)) then
        isotherm%multicontinuum_isotherm_list => sec_cont_isotherm_rxn
        sec_cont_isotherm_rxn%id = 1
      endif
      if (associated(sec_cont_prev_isotherm_rxn)) then
        sec_cont_prev_isotherm_rxn%next => sec_cont_isotherm_rxn
        sec_cont_isotherm_rxn%id = sec_cont_prev_isotherm_rxn%id + 1
      endif
      sec_cont_prev_isotherm_rxn => sec_cont_isotherm_rxn
      nullify(sec_cont_isotherm_rxn)
    endif
  enddo
  call InputPopBlock(input,option)

end subroutine IsothermRead

! ************************************************************************** !
subroutine IsothermConvertKDUnits(kd,kd_units,ikd_units,option)

  ! Converts units of isotherm reaction

  ! Author: Glenn Hammond
  ! Date: 03/16/21

  use Option_module
  use Units_module
   
  implicit none

  PetscReal :: kd
  character(MAXWORDLENGTH) :: kd_units
  PetscInt, intent(out) :: ikd_units
  type(option_type) :: option

  character(MAXWORDLENGTH) :: internal_units
  PetscErrorCode :: ierr

  ikd_units = UNINITIALIZED_INTEGER
  if (len_trim(kd_units) > 0) then
    ierr = 1
    internal_units = 'kg/m^3'
    Kd = Kd * UnitsConvertToInternal(kd_units,internal_units,option,ierr)
    if (ierr < 0) then
      ierr = 1
      internal_units = 'L/kg'
      Kd = Kd * UnitsConvertToInternal(kd_units,internal_units,option,ierr)
      if (ierr < 0) then
        option%io_buffer = 'Unrecognized kd_units: ' // trim(kd_units)
        call PrintErrMsg(option)
      else
        ikd_units = KD_UNIT_MLW_GSOIL
      endif
    else
      ikd_units = KD_UNIT_KG_M3_BULK
    endif
  endif

end subroutine IsothermConvertKDUnits

! ************************************************************************** !
subroutine RTotalSorbKD(rt_auxvar,global_auxvar,material_auxvar,isotherm, &
                        isotherm_rxn,option)
  ! 
  ! Computes the total sorbed component concentrations and
  ! derivative with respect to free-ion for the linear
  ! K_D model
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/30/2010
  ! 

  use Option_module
   
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(isotherm_type) :: isotherm
  type(isotherm_rxn_type) :: isotherm_rxn
  type(option_type) :: option
  
  PetscInt :: irxn
  PetscInt :: icomp
  PetscReal :: res
  PetscReal :: dres_dc
  PetscReal :: molality
  PetscReal :: tempreal
  PetscReal :: one_over_n
  PetscReal :: molality_one_over_n
  PetscReal :: kd_kgw_m3b  

  PetscInt, parameter :: iphase = 1

  do irxn = 1, isotherm%neqkdrxn
    icomp = isotherm%eqkdspecid(irxn)
    molality = rt_auxvar%pri_molal(icomp)
    if (isotherm%ikd_units == KD_UNIT_MLW_GSOIL) then
                   !KD units [mL water/g soil]
      kd_kgw_m3b = isotherm_rxn%eqisothermcoeff(irxn) * & 
                   global_auxvar%den_kg(iphase) * &
                   (1.d0-material_auxvar%porosity) * &
                   material_auxvar%soil_particle_density * &
                   1.d-3 ! convert mL water/g soil to m^3 water/kg soil

    else
      ! kd_unit = KD_UNIT_KG_M3_BULK
      kd_kgw_m3b = isotherm_rxn%eqisothermcoeff(irxn)              
    endif
    if (isotherm%eqkdmineral(irxn) > 0) then
      ! NOTE: mineral volume fraction here is solely a scaling factor.  It has 
      ! nothing to do with the soil volume; that is calculated through as a 
      ! function of porosity.
      kd_kgw_m3b = kd_kgw_m3b * &
                   (rt_auxvar%mnrl_volfrac(isotherm%eqkdmineral(irxn)))
    endif
    select case(isotherm%eqisothermtype(irxn))
      case(SORPTION_LINEAR)
        ! Csorb = Kd*Caq
        res = kd_kgw_m3b*molality
        dres_dc = kd_kgw_m3b
      case(SORPTION_LANGMUIR)
        ! Csorb = K*Caq*b/(1+K*Caq)
        tempreal = kd_kgw_m3b*molality
        res = tempreal*isotherm_rxn%eqisothermlangmuirb(irxn) / &
              (1.d0 + tempreal)
        dres_dc = res/molality - &
                  res / (1.d0 + tempreal) * tempreal / molality
      case(SORPTION_FREUNDLICH)
        ! Csorb = Kd*Caq**(1/n)
        one_over_n = 1.d0/isotherm_rxn%eqisothermfreundlichn(irxn)
        molality_one_over_n = molality**one_over_n
        res = kd_kgw_m3b*molality**one_over_n
        dres_dc = res/molality*one_over_n
      case default
        res = 0.d0
        dres_dc = 0.d0
    end select
    rt_auxvar%total_sorb_eq(icomp) = rt_auxvar%total_sorb_eq(icomp) + res
    rt_auxvar%dtotal_sorb_eq(icomp,icomp) = &
      rt_auxvar%dtotal_sorb_eq(icomp,icomp) + dres_dc 
  enddo

end subroutine RTotalSorbKD

end module Reaction_Isotherm_module
