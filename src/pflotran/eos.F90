module EOS_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use EOS_Water_module
  use EOS_Gas_module
  use EOSData_module

  implicit none

  private


  public :: EOSInit, &
            EOSRead, &
            EOSReferenceDensity, &
            EOSProcess, &
            EOSInputRecord, &
            AllEOSDBaseDestroy

contains

! ************************************************************************** !

subroutine EOSInit()

  implicit none

  call EOSWaterInit()
  call EOSGasInit()

  call EOSTableInitList()

end subroutine EOSInit

! ************************************************************************** !

subroutine EOSRead(input,option)

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word, subkeyword
  character(len=MAXWORDLENGTH) :: test_filename
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: tempreal
  PetscReal :: rks_tc = UNINITIALIZED_DOUBLE
  PetscReal :: rks_pc = UNINITIALIZED_DOUBLE
  PetscReal :: rks_acen = UNINITIALIZED_DOUBLE
  PetscReal :: rks_omegaa = UNINITIALIZED_DOUBLE
  PetscReal :: rks_omegab = UNINITIALIZED_DOUBLE
  PetscBool :: rks_hydrogen = PETSC_TRUE
  PetscBool :: rks_use_effective_properties = PETSC_TRUE
  PetscBool :: rks_use_cubic_root_solution = PETSC_FALSE
  PetscReal :: temparray(10)
  PetscReal :: test_t_high, test_t_low, test_p_high, test_p_low
  PetscInt :: test_n_temp, test_n_pres
  PetscBool :: test_uniform_temp, test_uniform_pres
  PetscErrorCode :: ierr

  input%ierr = 0

  call InputReadWord(input,option,keyword,PETSC_TRUE)
  call InputErrorMsg(input,option,'keyword','EOS')
  call StringToUpper(keyword)

  select case(trim(keyword))
    case('WATER')
      do
        temparray = 0.d0
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','EOS,WATER')
        call StringToUpper(keyword)
        select case(trim(keyword))
          case('DENSITY')
            temparray = 0.d0
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'DENSITY','EOS,WATER')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,DENSITY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                               'kg/m^3','EOS,WATER,DENSITY,CONSTANT',option)
              case('EXPONENTIAL','BRAGFLO')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'REFERENCE_DENSITY', &
                                   'EOS,WATER,DENSITY,EXPONENTIAL')
                call InputReadDouble(input,option,temparray(2))
                call InputErrorMsg(input,option,'REFERENCE_PRESSURE', &
                                   'EOS,WATER,DENSITY,EXPONENTIAL')
                call InputReadDouble(input,option,temparray(3))
                call InputErrorMsg(input,option,'WATER_COMPRESSIBILITY', &
                                   'EOS,WATER,DENSITY,EXPONENTIAL')
              case('LINEAR')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'REFERENCE_DENSITY', &
                                   'EOS,WATER,DENSITY,LINEAR')
                call InputReadDouble(input,option,temparray(2))
                call InputErrorMsg(input,option,'REFERENCE_PRESSURE', &
                                   'EOS,WATER,DENSITY,LINEAR')
                call InputReadDouble(input,option,temparray(3))
                call InputErrorMsg(input,option,'WATER_COMPRESSIBILITY', &
                                   'EOS,WATER,DENSITY,LINEAR')
              case('QUADRATIC')
                do
                  call InputReadPflotranString(input,option)
                  call InputReadStringErrorMsg(input,option, &
                                               'EOS,WATER,DENSITY,QUADRATIC')
                  if (InputCheckExit(input,option)) exit
                  if (InputError(input)) exit
                  call InputReadWord(input,option,subkeyword,PETSC_TRUE)
                  call InputErrorMsg(input,option,'subkeyword', &
                                       'EOS,WATER,DENSITY,QUADRATIC')
                  select case(trim(subkeyword))
                    case('REFERENCE_DENSITY')
                      call InputReadDouble(input,option,temparray(1))
                      call InputErrorMsg(input,option,'REFERENCE_DENSITY', &
                                         'EOS,WATER,DENSITY,QUADRATIC')
                    case('REFERENCE_PRESSURE')
                      call InputReadDouble(input,option,temparray(2))
                      call InputErrorMsg(input,option,'REFERENCE_PRESSURE', &
                                         'EOS,WATER,DENSITY,QUADRATIC')
                    case('WATER_COMPRESSIBILITY')
                      call InputReadDouble(input,option,temparray(3))
                      call InputErrorMsg(input,option,'WATER_COMPRESSIBILITY', &
                                         'EOS,WATER,DENSITY,QUADRATIC')
                    case default
                      call InputKeywordUnrecognized(subkeyword, &
                                'EOS,WATER,DENSITY,QUADRATIC',option)
                  end select
                enddo
              case('IFC67','DEFAULT','BATZLE_AND_WANG','TGDPB01','PLANAR', &
                   'TRANGENSTEIN','PAINTER')
              case default
                call InputKeywordUnrecognized(word,'EOS,WATER,DENSITY',option)
            end select
            call EOSWaterSetDensity(word,temparray)
          case('ENTHALPY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'ENTHALPY','EOS,WATER')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,ENTHALPY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                               'J/kmol','EOS,WATER,ENTHALPY,CONSTANT',option)
               case('IFC67','PAINTER','DEFAULT','PLANAR')
              case default
                call InputKeywordUnrecognized(word,'EOS,WATER,ENTHALPY',option)
            end select
            call EOSWaterSetEnthalpy(word,temparray)
          case('VISCOSITY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'VISCOSITY','EOS,WATER')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,VISCOSITY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                              'Pa-s','EOS,WATER,VISCOSITY,CONSTANT',option)
              case('DEFAULT','BATZLE_AND_WANG','GRABOWSKI')
              case default
                call InputKeywordUnrecognized(word,'EOS,WATER,VISCOSITY', &
                                              option)
            end select
            call EOSWaterSetViscosity(word,temparray)
          case('STEAM_DENSITY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'STEAM_DENSITY','EOS,WATER')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,STEAM_DENSITY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                           'kg/m^3','EOS,WATER,STEAM_DENSITY,CONSTANT',option)
              case('IFC67','DEFAULT','PLANAR','IF97')
              case default
                call InputKeywordUnrecognized(word,'EOS,WATER,STEAM_DENSITY', &
                                              option)
            end select
            call EOSWaterSetSteamDensity(word,temparray)
          case('STEAM_ENTHALPY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'STEAM_ENTHALPY','EOS,WATER')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,STEAM_ENTHALPY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                        'J/kmol','EOS,WATER,STEAM_ENTHALPY,CONSTANT',option)
              case('IFC67','DEFAULT','PLANAR','IF97')
              case default
                call InputKeywordUnrecognized(word, &
                       'EOS,WATER,STEAM_ENTHALPY',option)
            end select
            call EOSWaterSetSteamEnthalpy(word,temparray)
          case('TEST')
            if (option%global_rank == 0) then
              call InputReadDouble(input,option,test_t_low)
              call InputErrorMsg(input,option,'T_low', &
                                 'EOS,WATER,TEST')
              call InputReadDouble(input,option,test_t_high)
              call InputErrorMsg(input,option,'T_high', &
                                 'EOS,WATER,TEST')
              call InputReadDouble(input,option,test_p_low)
              call InputErrorMsg(input,option,'P_low', &
                                 'EOS,WATER,TEST')
              call InputReadDouble(input,option,test_p_high)
              call InputErrorMsg(input,option,'P_high', &
                                 'EOS,WATER,TEST')
              call InputReadInt(input,option,test_n_temp)
              call InputErrorMsg(input,option,'num_temperatures', &
                                 'EOS,WATER,TEST')
              call InputReadInt(input,option,test_n_pres)
              call InputErrorMsg(input,option,'num_pressures', &
                                 'EOS,WATER,TEST')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'temperature distribution type', &
                                 'EOS,WATER,TEST')
              if (StringCompareIgnoreCase(word,'uniform')) then
                test_uniform_temp = PETSC_TRUE
              else if (StringCompareIgnoreCase(word,'log')) then
                test_uniform_temp = PETSC_FALSE
              else
                option%io_buffer = 'Temperature distribution type "' // &
                  trim(word) // '" for EOS Water not recognized.'
                call printErrMsg(option)
              endif
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'pressure distribution type', &
                                 'EOS,WATER,TEST,')
              if (StringCompareIgnoreCase(word,'uniform')) then
                test_uniform_pres = PETSC_TRUE
              else if (StringCompareIgnoreCase(word,'log')) then
                test_uniform_pres = PETSC_FALSE
              else
                option%io_buffer = 'Pressure distribution type "' // &
                  trim(word) // '" for EOS Water not recognized.'
                call printErrMsg(option)
              endif
              call InputReadWord(input,option,word,PETSC_TRUE)
              test_filename = ''
              if (input%ierr == 0) then
                test_filename = word
              endif
              call EOSWaterTest(test_t_low,test_t_high,test_p_low,test_p_high, &
                                test_n_temp, test_n_pres, &
                                test_uniform_temp, test_uniform_pres, &
                                test_filename)
            endif
          case default
            call InputKeywordUnrecognized(keyword,'EOS,WATER',option)
        end select
      enddo
      string = ''
      call EOSWaterVerify(ierr,string)
      if (ierr /= 0) then
        option%io_buffer = 'Error in Water EOS'
        if (len_trim(string) > 1) then
          option%io_buffer = trim(option%io_buffer) // ': ' // trim(string)
        endif
        call printErrMsg(option)
      endif
    !
    case('ICE')
      do
        temparray = 0.d0
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','EOS,ICE')
        call StringToUpper(keyword)
        select case(trim(keyword))
          case('DENSITY')
            temparray = 0.d0
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'DENSITY','EOS,ICE')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,ICE,DENSITY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                               'kg/m^3','EOS,ICE,DENSITY,CONSTANT',option)
              case('DEFAULT','PAINTER')
                ! nothing
              case default
                call InputKeywordUnrecognized(word, &
                       'EOS,ICE,DENSITY',option)

            end select
            call EOSWaterSetDensityIce(word,temparray)
         end select
      end do
    !
    case('GAS')
      do
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','EOS,GAS')
        call StringToUpper(keyword)
        select case(trim(keyword))
          case('SURFACE_DENSITY','STANDARD_DENSITY')
            call InputReadDouble(input,option,tempreal)
            call InputErrorMsg(input,option,'VALUE', &
                               'EOS,GAS,REFERENCE_DENSITY')
            call InputReadAndConvertUnits(input,tempreal, &
                           'kg/m^3','EOS,GAS,REFERENCE_DENSITY',option)
            call EOSGasSetSurfaceDensity(tempreal)
          case('DENSITY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'DENSITY','EOS,GAS')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,GAS,DENSITY,CONSTANT')
                call InputReadAndConvertUnits(input,tempreal, &
                            'kmol/m^3','EOS,GAS,DENSITY,CONSTANT',option)
                call EOSGasSetDensityConstant(tempreal)
              case('RKS')
                ! if nothing is entered, it will calculate as hydrogen gas
                  do
                    call InputReadPflotranString(input,option)
                    call InputReadStringErrorMsg(input,option, &
                                                 'EOS GAS,RKS')
                    if (InputCheckExit(input,option)) exit
                    if (InputError(input)) exit
                    call InputReadWord(input,option,word,PETSC_TRUE)
                    call InputErrorMsg(input,option,'keyword', &
                                       'EOS GAS, RKS')
                    select case(trim(word))
                      case('HYDROGEN')
                        rks_hydrogen = PETSC_TRUE
                      case('NON-HYDROGEN')
                        rks_hydrogen = PETSC_FALSE
                      case('USE_EFFECTIVE_PROPERTIES')
                        rks_use_effective_properties = PETSC_TRUE
                      case('DONT_USE_EFFECTIVE_PROPERTIES')
                        rks_use_effective_properties = PETSC_FALSE
                      case('USE_CUBIC_ROOT_SOLUTION')
                        rks_use_cubic_root_solution = PETSC_TRUE
                      case('DONT_USE_CUBIC_ROOT_SOLUTION')
                        rks_use_cubic_root_solution = PETSC_FALSE
                      case('CRITICAL_TEMPERATURE','TC')
                        call InputReadDouble(input,option,rks_tc)
                        call InputErrorMsg(input,option, &
                                            'critical temperature for RKS', &
                                            'EOS GAS,RKS')
                      case('CRITICAL_PRESSURE','PC')
                        call InputReadDouble(input,option,rks_pc)
                        call InputErrorMsg(input,option, &
                                            'critical pressure for RKS', &
                                            'EOS GAS,RKS')
                      case('ACENTRIC,ACENTRIC_FACTOR','ACEN','AC')
                        ! acentric factor is only used for non-hydrogen gas
                        call InputReadDouble(input,option,rks_acen)
                        call InputErrorMsg(input,option, &
                                            'accentric factor for RKS', &
                                            'EOS GAS,RKS')
                      case('OMEGAA','A')
                        call InputReadDouble(input,option,rks_omegaa)
                        call InputErrorMsg(input,option, &
                                        'omega_a factor for RKS', &
                                            'EOS GAS,RKS')
                      case('OMEGAB','B')
                        call InputReadDouble(input,option,rks_omegab)
                        call InputErrorMsg(input,option, &
                                        'omega_b factor for RKS', &
                                            'EOS GAS,RKS')
                      case default
                        call InputKeywordUnrecognized(word, &
                                'EOS GAS,RKS',option)
                    end select
                enddo
                call EOSGasSetDensityRKS(rks_hydrogen, &
                                         rks_use_effective_properties, &
                                         rks_use_cubic_root_solution, &
                                         rks_tc,rks_pc,rks_acen, &
                                         rks_omegaa,rks_omegab)
              case('PR_METHANE')
                call EOSGasSetDensityPRMethane()
              case('IDEAL','DEFAULT')
                call EOSGasSetDensityIdeal()
              case default
                call InputKeywordUnrecognized(word,'EOS,GAS,DENSITY',option)
            end select
          case('ENTHALPY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'ENTHALPY','EOS,GAS')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,GAS,ENTHALPY,CONSTANT')
                call InputReadAndConvertUnits(input,tempreal, &
                                 'J/kmol','EOS,GAS,ENTHALPY,CONSTANT',option)
                call EOSGasSetEnergyConstant(tempreal)
              case('IDEAL_METHANE')
                call EOSGasSetEnergyIdealMethane()
              case('IDEAL','DEFAULT')
                call EOSGasSetEnergyIdeal()
              case default
                call InputKeywordUnrecognized(word,'EOS,GAS,ENTHALPY',option)
            end select
          case('VISCOSITY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'VISCOSITY','EOS,GAS')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,GAS,VISCOSITY,CONSTANT')
                call InputReadAndConvertUnits(input,tempreal, &
                                 'Pa-s','EOS,GAS,VISCOSITY,CONSTANT',option)
                call EOSGasSetViscosityConstant(tempreal)
              case('DEFAULT')
              case default
                call InputKeywordUnrecognized(word,'EOS,GAS,VISCOSITY',option)
            end select
          case('HENRYS_CONSTANT')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'HENRYS_CONSTANT','EOS,GAS')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,GAS,HENRYS_CONSTANT,CONSTANT')
                call EOSGasSetHenryConstant(tempreal)
              case('DEFAULT')
                call EOSGasSetHenry()
              case default
                call InputKeywordUnrecognized(word,'EOS,GAS,HENRYS_CONSTANT', &
                                              option)
            end select
          case('TEST')
            if (option%global_rank == 0) then
              call InputReadDouble(input,option,test_t_low)
              call InputErrorMsg(input,option,'T_low', &
                                 'EOS,GAS,TEST')
              call InputReadDouble(input,option,test_t_high)
              call InputErrorMsg(input,option,'T_high', &
                                 'EOS,GAS,TEST')
              call InputReadDouble(input,option,test_p_low)
              call InputErrorMsg(input,option,'P_low', &
                                 'EOS,GAS,TEST')
              call InputReadDouble(input,option,test_p_high)
              call InputErrorMsg(input,option,'P_high', &
                                 'EOS,GAS,TEST')
              call InputReadInt(input,option,test_n_temp)
              call InputErrorMsg(input,option,'num_temperatures', &
                                 'EOS,GAS,TEST')
              call InputReadInt(input,option,test_n_pres)
              call InputErrorMsg(input,option,'num_pressures', &
                                 'EOS,GAS,TEST')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'temperature distribution type', &
                                 'EOS,GAS,TEST')
              if (StringCompareIgnoreCase(word,'uniform')) then
                test_uniform_temp = PETSC_TRUE
              else if (StringCompareIgnoreCase(word,'log')) then
                test_uniform_temp = PETSC_FALSE
              else
                option%io_buffer = 'Temperature distribution type "' // &
                  trim(word) // '" for EOS Gas not recognized.'
                call printErrMsg(option)
              endif
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'pressure distribution type', &
                                 'EOS,GAS,TEST,')
              if (StringCompareIgnoreCase(word,'uniform')) then
                test_uniform_pres = PETSC_TRUE
              else if (StringCompareIgnoreCase(word,'log')) then
                test_uniform_pres = PETSC_FALSE
              else
                option%io_buffer = 'Pressure distribution type "' // &
                  trim(word) // '" for EOS Gas not recognized.'
                call printErrMsg(option)
              endif
              call InputReadWord(input,option,word,PETSC_TRUE)
              test_filename = ''
              if (input%ierr == 0) then
                test_filename = word
              endif
              call EOSGasTest(test_t_low,test_t_high,test_p_low,test_p_high, &
                              test_n_temp, test_n_pres, &
                              test_uniform_temp, test_uniform_pres, &
                              test_filename)
            endif
          case('FORMULA_WEIGHT')
            call InputReadDouble(input,option,tempreal)
            call InputErrorMsg(input,option,'VALUE','EOS,GAS,FORMULA_WEIGHT')
            call InputReadAndConvertUnits(input,tempreal, &
                             'g/mol','EOS,GAS,FORMULA_WEIGHT',option)
            call EOSGasSetFMWConstant(tempreal)

          case('DATABASE')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'EOS,GAS','DATABASE filename')
            call EOSGasSetEOSDBase(word,option)
          case('PVDG')
            call EOSGasSetPVDG(input,option)
          case default
            call InputKeywordUnrecognized(keyword,'EOS,GAS',option)
        end select
      enddo
      string = ''
      call EOSGasVerify(ierr,string)
      if (ierr == 5) then
        option%io_buffer = 'set to default value for RKS hydrogen'
        if (len_trim(string) > 1) then
          option%io_buffer =  trim(string) // ': ' // trim(option%io_buffer)
        endif
        call printMsg(option)
      !else if (ierr == 6) then
      !  option%io_buffer = 'set as default value for gas fmw'
      !  if (len_trim(string) > 1) then
      !    option%io_buffer =  trim(string) // ': ' // trim(option%io_buffer)
      !  endif
      !  call printMsg(option)
      else if (ierr /= 0) then
        option%io_buffer = 'Error in Gas EOS'
        if (len_trim(string) > 1) then
          option%io_buffer = trim(option%io_buffer) // ': ' // trim(string)
        endif
        call printErrMsg(option)
      endif

    case default
      call InputKeywordUnrecognized(keyword,'EOS',option)
  end select

end subroutine EOSRead

! **************************************************************************** !

subroutine EOSReferenceDensity(option)
  ! 
  ! Calculates reference densities for phases if not specified in input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/18
  ! 
  use Option_module
  
  implicit none
  type(option_type) :: option

  PetscReal :: vapor_saturation_pressure
  PetscReal :: vapor_density_kg
  PetscReal :: air_density_kmol
  PetscReal :: dum1, dum2
  PetscErrorCode :: ierr

  if (Initialized(option%liquid_phase)) then
    if (option%reference_density(option%liquid_phase) < 1.d-40) then
      call EOSWaterDensity(option%reference_temperature, &
                           option%reference_pressure, &
                           option%reference_density(option%liquid_phase), &
                           dum1,ierr)
    endif
  endif
  if (Initialized(option%gas_phase)) then
    if (option%reference_density(option%gas_phase) < 1.d-40) then
      ! assume saturated vapor pressure
      call EOSWaterSaturationPressure(option%reference_temperature, &
                                      vapor_saturation_pressure,dum1,ierr)
      call EOSWaterSteamDensityEnthalpy(option%reference_temperature, &
                                        vapor_saturation_pressure, &
                                        vapor_density_kg,dum1,dum2,ierr)
      ! call no-derivative version of EOSGasDensity
      call EOSGasDensity(option%reference_temperature, &
                         option%reference_pressure-vapor_saturation_pressure, &
                         air_density_kmol,ierr)
      option%reference_density(option%gas_phase) = &
        vapor_density_kg + air_density_kmol*FMWAIR
    endif
  endif

end subroutine EOSReferenceDensity

! **************************************************************************** !

subroutine EOSProcess(option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/28/17
  !
  ! Process EOS data after all input data have been read

  use Option_module

  implicit none

  type(option_type) :: option

  call EOSGasTableProcess(option)

  call EOSTableProcessList(option)

end subroutine EOSProcess

! **************************************************************************** !

subroutine EOSInputRecord()
  !
  ! Prints ingested equation of state information to the input record file.
  !
  ! Author: Jenn Frederick
  ! Date: 05/04/2016
  !
  use EOS_Water_module
  use EOS_Gas_module

  implicit none

  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'EQUATIONS OF STATE (EOS)'

  call EOSWaterInputRecord()
  call EOSGasInputRecord()

end subroutine EOSInputRecord

! ************************************************************************** !

subroutine AllEOSDBaseDestroy()
  !
  ! Author: Paolo Orsini
  ! Date: 12/14/15
  !

  implicit none

  call EOSTableDestroyList()

  call EOSGasDBaseDestroy()

end subroutine AllEOSDBaseDestroy

! ************************************************************************** !

end module EOS_module
