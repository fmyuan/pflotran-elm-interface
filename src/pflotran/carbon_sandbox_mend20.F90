module Carbon_Sandbox_MEND20_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Carbon_Sandbox_Base_class
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, extends(carbon_sandbox_base_type), public :: carbon_sandbox_mend20_type
    character(len=MAXWORDLENGTH) :: PCO_species_name
    character(len=MAXWORDLENGTH) :: PNO_species_name
    character(len=MAXWORDLENGTH) :: PCH_species_name
    character(len=MAXWORDLENGTH) :: PNH_species_name
    character(len=MAXWORDLENGTH) :: MC_species_name
    character(len=MAXWORDLENGTH) :: MN_species_name
    character(len=MAXWORDLENGTH) :: DC_species_name
    character(len=MAXWORDLENGTH) :: DN_species_name
    character(len=MAXWORDLENGTH) :: QC_species_name
    character(len=MAXWORDLENGTH) :: QN_species_name
    character(len=MAXWORDLENGTH) :: BAC_species_name
    character(len=MAXWORDLENGTH) :: BAN_species_name
    character(len=MAXWORDLENGTH) :: BDC_species_name
    character(len=MAXWORDLENGTH) :: BDN_species_name
    character(len=MAXWORDLENGTH) :: EPCO_species_name
    character(len=MAXWORDLENGTH) :: EPNO_species_name !
    character(len=MAXWORDLENGTH) :: EPCH_species_name
    character(len=MAXWORDLENGTH) :: EPNH_species_name
    character(len=MAXWORDLENGTH) :: EMC_species_name
    character(len=MAXWORDLENGTH) :: EMN_species_name
    character(len=MAXWORDLENGTH) :: NH4_species_name
    character(len=MAXWORDLENGTH) :: NO3_species_name
    ! C = carbon pool, N = nitrogen pool
    ! P = particulate, M = mineral, D = dissolved, Q = sorbing
    ! BA = active biomass, BD = dormant biomass
    ! O = oxidative, H = hydrolytic
    ! E = enzyme
    PetscInt :: PCO_species_index  ! oxidative particulate organic matter (POM)
    PetscInt :: PNO_species_index
    PetscInt :: PCH_species_index  ! hydrolytic POM
    PetscInt :: PNH_species_index
    PetscInt :: MC_species_index   ! mineral associated organic matter (MOM)
    PetscInt :: MN_species_index
    PetscInt :: DC_species_index   ! dissolved organic matter (DOC)
    PetscInt :: DN_species_index
    PetscInt :: QC_species_index   ! active MOM interacting with DOC (sorbing)
    PetscInt :: QN_species_index
    PetscInt :: BAC_species_index  ! active biomass
    PetscInt :: BAN_species_index
    PetscInt :: BDC_species_index  ! dormant biomass
    PetscInt :: BDN_species_index
    PetscInt :: EPCO_species_index ! oxidative enzymes decomposing POM
    PetscInt :: EPNO_species_index
    PetscInt :: EPCH_species_index ! hydroylitic enzymes decomposing POM
    PetscInt :: EPNH_species_index
    PetscInt :: EMC_species_index  ! enzymes decomposing MOM
    PetscInt :: EMN_species_index
    PetscInt :: NH4_species_index  ! ammonium
    PetscInt :: NO3_species_index  ! nitrate
    PetscReal :: LF_0
    PetscReal :: r_0
    PetscReal :: fINP
    PetscReal :: Vd_PO
    PetscReal :: Vd_PH
    PetscReal :: Vd_M
    PetscReal :: K_PO
    PetscReal :: K_PH
    PetscReal :: K_M
    PetscReal :: QC_max
    PetscReal :: K_ba
    PetscReal :: k_des
    PetscReal :: r_E
    PetscReal :: p_EP
    PetscReal :: fp_EM
    PetscReal :: f_D
    PetscReal :: g_D
    PetscReal :: g_PO
    PetscReal :: V_g
    PetscReal :: alpha
    PetscReal :: K_D
    PetscReal :: Y_g_T_ref
    PetscReal :: k_Yg
    PetscReal :: Q10
    PetscReal :: gamma
    PetscReal :: beta
    PetscReal :: psi_A2D
    PetscReal :: tau
    PetscReal :: omega
    PetscReal :: VN_imNH4
    PetscReal :: VN_imNO3
    PetscReal :: KS_NH4
    PetscReal :: KS_NO3
    PetscReal :: VN_nit
    PetscReal :: VN_denit
    PetscReal :: reference_temperature
  contains
    procedure, public :: ReadInput => CarbonMEND20ReadInput
    procedure, public :: Setup => CarbonMEND20Setup
    procedure, public :: Evaluate => CarbonMEND20Evaluate
    procedure, public :: Strip => CarbonMEND20Strip
  end type carbon_sandbox_mend20_type

  type, extends(carbon_sandbox_rxn_base_type), public :: &
                                                 carbon_sandbox_rxn_mend_type
  contains
    procedure, public :: ReadInput => CarbonRxnMEND20ReadInput
    procedure, public :: Setup => CarbonRxnMEND20Setup
    procedure, public :: Evaluate => CarbonRxnMEND20Evaluate
    procedure, public :: Strip => CarbonRxnMEND20Strip
  end type carbon_sandbox_rxn_mend_type

  public :: CarbonMEND20Create, &
            CarbonRxnMEND20Create

contains

! ************************************************************************** !

function CarbonMEND20Create()
  !
  ! Allocates and initializes the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_mend20_type), pointer :: CarbonMEND20Create

  class(carbon_sandbox_mend20_type), pointer :: this

  allocate(this)
  call CarbonBaseInit(this)
  this%PCO_species_name = ''
  this%PNO_species_name = ''
  this%PCH_species_name = ''
  this%PNH_species_name = ''
  this%MC_species_name = ''
  this%MN_species_name = ''
  this%DC_species_name = ''
  this%DN_species_name = ''
  this%QC_species_name = ''
  this%QN_species_name = ''
  this%BAC_species_name = ''
  this%BAN_species_name = ''
  this%BDC_species_name = ''
  this%BDN_species_name = ''
  this%EPCO_species_name = ''
  this%EPNO_species_name = ''
  this%EPCH_species_name = ''
  this%EPNH_species_name = ''
  this%EMC_species_name = ''
  this%EMN_species_name = ''
  this%NH4_species_name = ''
  this%NO3_species_name = ''

  this%PCO_species_index = UNINITIALIZED_INTEGER
  this%PNO_species_index = UNINITIALIZED_INTEGER
  this%PCH_species_index = UNINITIALIZED_INTEGER
  this%PNH_species_index = UNINITIALIZED_INTEGER
  this%MC_species_index = UNINITIALIZED_INTEGER
  this%MN_species_index = UNINITIALIZED_INTEGER
  this%DC_species_index = UNINITIALIZED_INTEGER
  this%DN_species_index = UNINITIALIZED_INTEGER
  this%QC_species_index = UNINITIALIZED_INTEGER
  this%QN_species_index = UNINITIALIZED_INTEGER
  this%BAC_species_index = UNINITIALIZED_INTEGER
  this%BAN_species_index = UNINITIALIZED_INTEGER
  this%BDC_species_index = UNINITIALIZED_INTEGER
  this%BDN_species_index = UNINITIALIZED_INTEGER
  this%EPCO_species_index = UNINITIALIZED_INTEGER
  this%EPNO_species_index = UNINITIALIZED_INTEGER
  this%EPCH_species_index = UNINITIALIZED_INTEGER
  this%EPNH_species_index = UNINITIALIZED_INTEGER
  this%EMC_species_index = UNINITIALIZED_INTEGER
  this%EMN_species_index = UNINITIALIZED_INTEGER
  this%NH4_species_index = UNINITIALIZED_INTEGER
  this%NO3_species_index = UNINITIALIZED_INTEGER

  this%LF_0 = UNINITIALIZED_DOUBLE
  this%r_0 = UNINITIALIZED_DOUBLE
  this%fINP = UNINITIALIZED_DOUBLE
  this%Vd_PO = UNINITIALIZED_DOUBLE
  this%Vd_PH = UNINITIALIZED_DOUBLE
  this%Vd_M = UNINITIALIZED_DOUBLE
  this%K_PO = UNINITIALIZED_DOUBLE
  this%K_PH = UNINITIALIZED_DOUBLE
  this%K_M = UNINITIALIZED_DOUBLE
  this%QC_max = UNINITIALIZED_DOUBLE
  this%K_ba = UNINITIALIZED_DOUBLE
  this%k_des = UNINITIALIZED_DOUBLE
  this%r_E = UNINITIALIZED_DOUBLE
  this%p_EP = UNINITIALIZED_DOUBLE
  this%fp_EM = UNINITIALIZED_DOUBLE
  this%f_D = UNINITIALIZED_DOUBLE
  this%g_D = UNINITIALIZED_DOUBLE
  this%g_PO = UNINITIALIZED_DOUBLE
  this%V_g = UNINITIALIZED_DOUBLE
  this%alpha = UNINITIALIZED_DOUBLE
  this%K_D = UNINITIALIZED_DOUBLE
  this%Y_g_T_ref = UNINITIALIZED_DOUBLE
  this%k_Yg = UNINITIALIZED_DOUBLE
  this%Q10 = UNINITIALIZED_DOUBLE
  this%gamma = UNINITIALIZED_DOUBLE
  this%beta = UNINITIALIZED_DOUBLE
  this%psi_A2D = UNINITIALIZED_DOUBLE
  this%tau = UNINITIALIZED_DOUBLE
  this%omega = UNINITIALIZED_DOUBLE
  this%VN_imNH4 = UNINITIALIZED_DOUBLE
  this%VN_imNO3 = UNINITIALIZED_DOUBLE
  this%KS_NH4 = UNINITIALIZED_DOUBLE
  this%KS_NO3 = UNINITIALIZED_DOUBLE
  this%VN_nit = UNINITIALIZED_DOUBLE
  this%VN_denit = UNINITIALIZED_DOUBLE
  this%reference_temperature = UNINITIALIZED_DOUBLE

  CarbonMEND20Create => this

end function CarbonMEND20Create

! ************************************************************************** !

subroutine CarbonMEND20ReadInput(this,input,option)
  !
  ! Reads parameters from input deck
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_mend20_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: err_string
  PetscReal :: tempreal
  PetscBool :: found

  err_string = 'CHEMISTRY,CARBON_SANDBOX,MEND20'
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',err_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call CarbonBaseReadSelectCase(this,input,keyword,found,err_string,option)
    if (found) cycle

    select case(keyword)
      case('POMC_O_C_SPECIES_NAME')
        call InputReadWord(input,option,this%PCO_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('POMC_O_N_SPECIES_NAME')
        call InputReadWord(input,option,this%PNO_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('POMC_H_C_SPECIES_NAME')
        call InputReadWord(input,option,this%PCH_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('POMC_H_N_SPECIES_NAME')
        call InputReadWord(input,option,this%PNH_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('MOMC_C_SPECIES_NAME')
        call InputReadWord(input,option,this%MC_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('MOMC_N_SPECIES_NAME')
        call InputReadWord(input,option,this%MN_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('DOMC_C_SPECIES_NAME')
        call InputReadWord(input,option,this%DC_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('DOMC_N_SPECIES_NAME')
        call InputReadWord(input,option,this%DN_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('QOMC_C_SPECIES_NAME')
        call InputReadWord(input,option,this%QC_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('QOMC_N_SPECIES_NAME')
        call InputReadWord(input,option,this%QN_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('MBAC_C_SPECIES_NAME')
        call InputReadWord(input,option,this%BAC_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('MBAC_N_SPECIES_NAME')
        call InputReadWord(input,option,this%BAN_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('MBDC_C_SPECIES_NAME')
        call InputReadWord(input,option,this%BDC_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('MBDC_N_SPECIES_NAME')
        call InputReadWord(input,option,this%BDN_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('EPCO_C_SPECIES_NAME')
        call InputReadWord(input,option,this%EPCO_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('EPCO_N_SPECIES_NAME')
        call InputReadWord(input,option,this%EPNO_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('EPCH_C_SPECIES_NAME')
        call InputReadWord(input,option,this%EPCH_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('EPCH_N_SPECIES_NAME')
        call InputReadWord(input,option,this%EPNH_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('EMC_C_SPECIES_NAME')
        call InputReadWord(input,option,this%EMC_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('EMC_N_SPECIES_NAME')
        call InputReadWord(input,option,this%EMN_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('NH4_SPECIES_NAME')
        call InputReadWord(input,option,this%NH4_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('NO3_SPECIES_NAME')
        call InputReadWord(input,option,this%NO3_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('K_YG')
        call InputReadDouble(input,option,this%k_Yg)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = '1/C' ! 1/C
        call InputReadAndConvertUnits(input,this%k_Yg,internal_units, &
                                      keyword,option)
      case('PSI_A2D')
        call InputReadDouble(input,option,this%psi_A2D)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'Pa' ! MPa
        call InputReadAndConvertUnits(input,this%psi_A2D,internal_units, &
                                      keyword,option)
      case('K_DES')
        call InputReadDouble(input,option,this%k_des)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'kg/m^3-s' ! mg C / cm^3 soil / h
        call InputReadAndConvertUnits(input,this%k_des,internal_units, &
                                      keyword,option)
      case('VN_NIT','VN_DENIT')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = '1/s' ! 1 / h
        call InputReadAndConvertUnits(input,tempreal,internal_units, &
                                      keyword,option)
        select case(keyword)
          case('VN_NIT')
            this%VN_nit = tempreal
          case('VN_DENIT')
            this%VN_denit = tempreal
        end select
      case('VDC_PO','VDC_PH','VDC_M','R_E','V_G','VN_IMNH4','VN_IMNO3')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'kg/kg-s' ! mg C / mg C / h
        call InputReadAndConvertUnits(input,tempreal,internal_units, &
                                      keyword,option)
        select case(keyword)
          case('VDC_PO')
            this%Vd_PO = tempreal
          case('VDC_PH')
            this%Vd_PH = tempreal
          case('VDC_M')
            this%Vd_M = tempreal
          case('R_E')
            this%r_E = tempreal
          case('V_G')
            this%V_g = tempreal
          case('VN_IMNH4')
            this%VN_imNH4 = tempreal
          case('VN_IMNO3')
            this%VN_imNO3 = tempreal
        end select
      case('K_PO','K_PH','K_M','K_D','KS_NH4','KS_NO3','QC_MAX')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'kg/m^3' ! mg C/cm^3 soil
        call InputReadAndConvertUnits(input,tempreal,internal_units, &
                                      keyword,option)
        select case(keyword)
          case('K_PO')
            this%K_PO = tempreal
          case('K_PH')
            this%K_PH = tempreal
          case('K_M')
            this%K_M = tempreal
          case('K_D')
            this%K_D = tempreal
          case('KS_NH4')
            this%KS_NH4 = tempreal
          case('KS_NO3')
            this%KS_NO3 = tempreal
          case('QC_MAX')
            this%QC_max = tempreal
        end select
      case('K_BA')
        call InputReadDouble(input,option,this%K_ba)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'm^3/kg' ! cm^3 soil / mg C
        call InputReadAndConvertUnits(input,this%K_ba,internal_units, &
                                      keyword,option)
      case('LF_0','R_0','FINP','P_EP','FP_EM','F_D','G_D','G_PO','Q10', &
           'ALPHA','Y_G_T_REF','GAMMA','BETA','TAU','OMEGA') ! unitless
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,err_string)
        select case(keyword)
          case('LF_0')
            this%LF_0 = tempreal
          case('R_0')
            this%r_0 = tempreal
          case('FINP')
            this%fINP = tempreal
          case('P_EP')
            this%p_EP = tempreal
          case('FP_EM')
            this%fp_EM = tempreal
          case('F_D')
            this%f_D = tempreal
          case('G_D')
            this%g_D = tempreal
          case('G_PO')
            this%g_PO = tempreal
          case('Q10')
            this%g_PO = tempreal
          case('ALPHA')
            this%alpha = tempreal
          case('Y_G_T_REF')
            this%Y_g_T_ref = tempreal
          case('GAMMA')
            this%gamma = tempreal
          case('BETA')
            this%beta = tempreal
          case('TAU')
            this%tau = tempreal
          case('OMEGA')
            this%omega = tempreal
        end select
      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,this%reference_temperature)
        call InputErrorMsg(input,option,keyword,err_string)
      case default
        call InputKeywordUnrecognized(input,keyword,err_string,option)
    end select
  enddo

end subroutine CarbonMEND20ReadInput

! ************************************************************************** !

subroutine CarbonMEND20Setup(this,reaction,option)
  !
  ! Configures the reaction and associated data structures
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module

  class(carbon_sandbox_mend20_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  call CarbonBaseSetup(this,reaction,option)

  call SetupMobile(this%PCO_species_index,this%PCO_species_name)
  call SetupMobile(this%PNO_species_index,this%PNO_species_name)
  call SetupMobile(this%PCH_species_index,this%PCH_species_name)
  call SetupMobile(this%PNH_species_index,this%PNH_species_name)
  call SetupMobile(this%DC_species_index,this%DC_species_name)
  call SetupMobile(this%DN_species_index,this%DN_species_name)
  call SetupMobile(this%BAC_species_index,this%BAC_species_name)
  call SetupMobile(this%BAN_species_index,this%BAN_species_name)
  call SetupMobile(this%BDC_species_index,this%BDC_species_name)
  call SetupMobile(this%BDN_species_index,this%BDN_species_name)
  call SetupMobile(this%EPNO_species_index,this%EPNO_species_name)
  call SetupMobile(this%EPCH_species_index,this%EPCH_species_name)
  call SetupMobile(this%EPNH_species_index,this%EPNH_species_name)
  call SetupMobile(this%EMC_species_index,this%EMC_species_name)
  call SetupMobile(this%EMN_species_index,this%EMN_species_name)
  call SetupMobile(this%NH4_species_index,this%NH4_species_name)
  call SetupMobile(this%NO3_species_index,this%NO3_species_name)

  call SetupImmobile(this%MC_species_index,this%MC_species_name)
  call SetupImmobile(this%MN_species_index,this%MN_species_name)
  call SetupImmobile(this%QC_species_index,this%QC_species_name)
  call SetupImmobile(this%QN_species_index,this%QN_species_name)

#if 0
  if (this%MC_species_index <= 0) then
    option%io_buffer = 'Species "' // trim(this%MC_species_name) // &
      '" not found among the available immobile species. Mineral associated &
      &biomass must be an immobile species.'
    call PrintErrMsg(option)
  endif
  this%P_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%P_species_name,reaction,option)
  this%QC_species_index = &
    ReactionImGetSpeciesIDFromName(this%QC_species_name,reaction%immobile, &
                                   PETSC_FALSE,option) + &
    reaction%offset_immobile
  if (this%MC_species_index <= 0) then
    option%io_buffer = 'Species "' // trim(this%MC_species_name) // &
      '" not found among the available immobile species. Adsorbed biomass &
      &must be an immobile species.'
    call PrintErrMsg(option)
  endif
#endif

  call CheckUninitialized(this%LF_0,'LF_0')
  call CheckUninitialized(this%r_0,'R_0')
  call CheckUninitialized(this%fINP,'FINP')
  call CheckUninitialized(this%Vd_PO,'VDC_PO')
  call CheckUninitialized(this%Vd_PH,'VDC_PH')
  call CheckUninitialized(this%Vd_M,'VDC_M')
  call CheckUninitialized(this%K_PO,'K_PO')
  call CheckUninitialized(this%K_PH,'K_PH')
  call CheckUninitialized(this%K_M,'K_M')
  call CheckUninitialized(this%QC_max,'QC_MAX')
  call CheckUninitialized(this%K_ba,'K_BA')
  call CheckUninitialized(this%k_des,'K_DES')
  call CheckUninitialized(this%r_E,'R_E')
  call CheckUninitialized(this%k_des,'K_des')
  call CheckUninitialized(this%p_EP,'P_EP')
  call CheckUninitialized(this%fp_EM,'FP_EM')
  call CheckUninitialized(this%f_D,'F_D')
  call CheckUninitialized(this%g_D,'G_D')
  call CheckUninitialized(this%g_PO,'G_PO')
  call CheckUninitialized(this%V_g,'V_G')
  call CheckUninitialized(this%alpha,'ALPHA')
  call CheckUninitialized(this%K_D,'K_D')
  call CheckUninitialized(this%Y_g_T_ref,'Y_G_T_REF')
  call CheckUninitialized(this%k_Yg,'K_YG')
  call CheckUninitialized(this%Q10,'Q10')
  call CheckUninitialized(this%gamma,'GAMMA')
  call CheckUninitialized(this%beta,'BETA')
  call CheckUninitialized(this%psi_A2D,'PSI_A2D')
  call CheckUninitialized(this%tau,'TAU')
  call CheckUninitialized(this%omega,'OMEGA')
  call CheckUninitialized(this%VN_imNH4,'VN_IMNH4')
  call CheckUninitialized(this%VN_imNO3,'VN_IMNO3')
  call CheckUninitialized(this%KS_NH4,'KS_NH4')
  call CheckUninitialized(this%KS_NO3,'KS_NO3')
  call CheckUninitialized(this%VN_nit,'VN_NIT')
  call CheckUninitialized(this%VN_denit,'VN_DENIT')

  if (Uninitialized(this%reference_temperature)) then
    this%reference_temperature = option%flow%reference_temperature
  endif

  call this%EnforceConcentrationUnits(CARBON_UNITS_MOLE_PER_KG_SOIL,option)

contains

subroutine CheckUninitialized(var,varname)

  PetscReal :: var
  character(len=*) :: varname

  if (Uninitialized(var)) then
    option%io_buffer = trim(varname) // 'not defined in MEND20 carbon sandbox.'
    call PrintErrMsg(option)
  endif

end subroutine CheckUninitialized

subroutine SetupMobile(ivar,varname)

  PetscInt :: ivar
  character(len=*) :: varname

  ivar = ReactionAuxGetPriSpecIDFromName(varname,reaction,option)
  if (Uninitialized(ivar)) then
    option%io_buffer = trim(varname) // 'not defined in MEND20 carbon sandbox.'
    call PrintErrMsg(option)
  endif

end subroutine SetupMobile

subroutine SetupImmobile(ivar,varname)

  PetscInt :: ivar
  character(len=*) :: varname

  ivar = ReactionImGetSpeciesIDFromName(varname,reaction%immobile, &
                                        PETSC_FALSE,option) + &
    reaction%offset_immobile

  if (ivar <= 0) then
    option%io_buffer = 'Species "' // trim(varname) // '" in MEND20 carbon &
      &sandbox not found among the available immobile species. It must be &
      &an immobile species.'
    call PrintErrMsg(option)
  endif

end subroutine SetupImmobile

end subroutine CarbonMEND20Setup

! ************************************************************************** !

subroutine CarbonMEND20Evaluate(this,Residual,Jacobian,compute_derivative, &
                              rt_auxvar,global_auxvar,material_auxvar, &
                              reaction,option)
  !
  ! Evaluates the rate expression
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Global_Aux_module
  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Inhibition_Aux_module
  use Reactive_Transport_Aux_module
  use Utility_module

  class(carbon_sandbox_mend20_type) :: this
  class(reaction_rt_type) :: reaction
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(:)
  PetscReal :: Jacobian(:,:)
  PetscBool :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option

  PetscReal :: PCO, PNO, PCH, PNH, MC, MN, DC, DN, QC, QN, BAC, BAN, BDC, BDN
  PetscReal :: EPCO, EPNO, EPCH, EPNH, EMC, EMN, NH4, NO3
  PetscReal :: F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14
  PetscReal :: F15, F15_EPCO, F15_EPCH, F15_EMC
  PetscReal :: F16, F16_EPCO, F16_EPCH, F16_EMC
  PetscReal :: FN_im_NH4toBA, FN_im_NO3toBA
!  PetscReal :: FN_mn_BA
  PetscReal :: dPCO, dPCH, dMC, dQC, dDC, dBAC, dBDC, dEPCO, dEPCH, dEMC, dCO2
  PetscReal :: YN_g, one_over_Y_g
  PetscReal :: four_thirds
  PetscReal :: FN_nit
  PetscReal :: FN_nit_denit
  PetscReal :: f_O2
  PetscReal :: WFP
  PetscReal :: FN_denit
  PetscReal :: V_m, k_ads
  PetscReal :: conc_units_conversion
  PetscReal :: t
  PetscReal :: CN_BA, CN_BAmin, CN_BAmax
  PetscReal :: DC_monod

  PetscReal :: I_PO, I_PH, I_D

  I_PO = 0.d0
  I_PO = (1.d0/I_PO)*0.d0
  I_PH = I_PO
  I_D = I_PO

  CN_BA = 0.5d0
  CN_BAmin = 0.d0
  CN_BAmax = 1.d0

  call this%MapStateVariables(rt_auxvar,global_auxvar,material_auxvar, &
                              reaction,option)
  t = this%aux%temperature
  WFP = 1.d0

  k_ads = this%k_des * this%K_BA
  V_m = this%alpha * this%V_g / (1.d0 - this%alpha)

  ! mgC/gsoil = molC/kgsoil*gC/molC*mgC/gC*kgsoil/gsoil
  conc_units_conversion = 30.d0*1000.d0*1.d-3
  PCO = this%aux%conc(this%PCO_species_index)*conc_units_conversion
  PNO = this%aux%conc(this%PNO_species_index)*conc_units_conversion
  PCH = this%aux%conc(this%PCH_species_index)*conc_units_conversion
  PNH = this%aux%conc(this%PNH_species_index)*conc_units_conversion
  MC = this%aux%conc(this%MC_species_index)*conc_units_conversion
  MN = this%aux%conc(this%MN_species_index)*conc_units_conversion
  DC = this%aux%conc(this%DC_species_index)*conc_units_conversion
  DN = this%aux%conc(this%DN_species_index)*conc_units_conversion
  QC = this%aux%conc(this%QC_species_index)*conc_units_conversion
  QN = this%aux%conc(this%QN_species_index)*conc_units_conversion
  BAC = this%aux%conc(this%BAC_species_index)*conc_units_conversion
  BAN = this%aux%conc(this%BAN_species_index)*conc_units_conversion
  BDC = this%aux%conc(this%BDC_species_index)*conc_units_conversion
  BDN = this%aux%conc(this%BDN_species_index)*conc_units_conversion
  EPCO = this%aux%conc(this%EPCO_species_index)*conc_units_conversion
  EPNO = this%aux%conc(this%EPNO_species_index)*conc_units_conversion
  EPCH = this%aux%conc(this%EPCH_species_index)*conc_units_conversion
  EPNH = this%aux%conc(this%EPNH_species_index)*conc_units_conversion
  EMC = this%aux%conc(this%EMC_species_index)*conc_units_conversion
  EMN = this%aux%conc(this%EMN_species_index)*conc_units_conversion
  NH4 = this%aux%conc(this%NH4_species_index)*conc_units_conversion
  NO3 = this%aux%conc(this%NO3_species_index)*conc_units_conversion

  one_over_Y_g = 1.d0 / (this%Y_g_T_ref - &
                         this%k_Yg * (t - this%reference_temperature))
  YN_g = ((CN_BA-CN_BAmin)/(CN_BAmax-CN_BAmin))**this%omega
  four_thirds = 4.d0 / 3.d0
  f_O2 = (((1.d0 - WFP)**four_thirds) / &
          ((0.5**four_thirds)+((1.d0-WFP)**four_thirds)))

  DC_monod = DC / (this%K_D+DC)

  F1 = this%Vd_PO * EPCO * PCO / (this%K_PO+PCO)
  F2 = this%Vd_PH * EPCH * PCH / (this%K_PH+PCH)
  F3 = this%Vd_M * EMC * MC / (this%K_M+MC)
  F4 = k_ads * (1.d0 - QC/this%QC_max) * DC
  F5 = this%k_des * (QC/this%QC_max)
  F6 = one_over_Y_g * (this%V_g + V_m) * BAC * DC_monod
  F7 = (1.d0 - DC_monod) * V_m * BAC
  F8 = DC_monod * V_m * BDC
  F9 = (one_over_Y_g - 1.d0) * this%V_g * BAC * DC_monod
  F10 = (one_over_Y_g - 1.d0) * V_m * BAC * DC_monod
  F11 = max(0.d0,BAC-BAN*CN_BAmax)
  F12 = this%beta * V_m * BDC
  F13 = max(0.d0,BDC-BDN*CN_BAmax)
  F14 = this%gamma * V_m * BDC
  F15_EPCO = PCO / (PCO+PCH) * this%p_EP * V_m * BAC
  F15_EPCH = PCH / (PCO+PCH) * this%p_EP * V_m * BAC
  F15_EMC = this%fp_EM * this%p_EP * V_m * BAC
  F15 = F15_EPCO + F15_EPCH + F15_EMC
  F16_EPCO = this%r_E * EPCO
  F16_EPCH = this%r_E * EPCH
  F16_EMC = this%r_E * EMC
  F16 = F16_EPCO + F16_EPCH + F16_EMC
  FN_im_NH4toBA = this%VN_imNH4 * YN_g * BAC * NH4 / &
                  (this%KS_NH4 * (1.d0 + NH4 / this%KS_NH4 + &
                                  NO3 / this%KS_NO3 + BAC / this%KS_NH4))
  FN_im_NO3toBA = this%VN_imNO3 * YN_g * BAC * NO3 / &
                  (this%KS_NO3 * (1.d0 + NH4 / this%KS_NH4 + &
                                  NO3 / this%KS_NO3 + BAC / this%KS_NO3))
!  FN_mn_BA = (1.d0 - YN_g) * FN6
  FN_nit = this%VN_nit * NH4
  FN_nit_denit = FN_nit * (1.d0 - f_O2)
  FN_denit = this%VN_denit * NO3

  dPCO = I_PO + (1.d0-this%g_D)*this%g_PO*F12 - F1
  dPCH = I_PH + (1.d0-this%g_D)*(1.d0-this%g_PO)*F12 - F2
  dMC = (1.d0-this%f_D)*(F1+F2) - F3
  dQC = F4 - F5
  dDC = I_D + this%f_D*(F1+F2) + F3 + this%g_D*F12 + F14 - F6 - (F4-F5)
  dBAC = F6 - (F7-F8) - (F9+F10+F11) - F14 - F15
  dBDC = (F7-F8) - (F12+F13)
  dEPCO = F15_EPCO - F16_EPCO
  dEPCH = F15_EPCH - F16_EPCH
  dEMC = F15_EMC - F16_EMC
  dCO2 = (F9 + F10 + F11) + (F12 + F13)

!  dPCH = I_PH + (1.d0-this%g_D)*this%g_PO*F12 - F1


end subroutine CarbonMEND20Evaluate

! ************************************************************************** !

subroutine CarbonMEND20Strip(this)
  !
  ! Destroys members of the carbon sandbox MEND20 object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Utility_module

  class(carbon_sandbox_mend20_type) :: this

  call CarbonBaseStrip(this)

end subroutine CarbonMEND20Strip

! ************************************************************************** !

function CarbonRxnMEND20Create()
  !
  ! Allocates and initializes the carbon sandbox rxn object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_rxn_mend_type), pointer :: CarbonRxnMEND20Create

  class(carbon_sandbox_rxn_mend_type), pointer :: this

  allocate(this)
  this%rate_constant = UNINITIALIZED_DOUBLE
  this%reaction_string = ''
  nullify(this%aux)
  nullify(this%reaction_equation)
  nullify(this%next)

  CarbonRxnMEND20Create => this

end function CarbonRxnMEND20Create

! ************************************************************************** !

subroutine CarbonRxnMEND20ReadInput(this,input,option)
  !
  ! Reads reaction parameters from input deck
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_rxn_mend_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: err_string

  err_string = 'CHEMISTRY,CARBON_SANDBOX,REACTION_NETWORK'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',err_string)

    call StringToUpper(keyword)
    select case(trim(keyword))
      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,keyword,err_string)
      case('REACTION_EQUATION')
        this%reaction_string = trim(input%buf)
      case default
        call InputKeywordUnrecognized(input,keyword,err_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT not defined in ' // trim(err_string)
    call PrintErrMsg(option)
  endif
  if (len_trim(this%reaction_string) <= 1) then
    option%io_buffer = 'A REACTION_EQUATION is not defined in ' // &
      trim(err_string)
    call PrintErrMsg(option)
  endif

end subroutine CarbonRxnMEND20ReadInput

! ************************************************************************** !

subroutine CarbonRxnMEND20Setup(this,aux,reaction,option)
  !
  ! Configures the reaction and associated data structures
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Reaction_Aux_module

  class(carbon_sandbox_rxn_mend_type) :: this
  type(carbon_sandbox_aux_type), pointer :: aux
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  call CarbonRxnBaseSetup(this,aux,reaction,option)

end subroutine CarbonRxnMEND20Setup

! ************************************************************************** !

subroutine CarbonRxnMEND20Evaluate(this,Residual,Jacobian,option)
  !
  ! Evaluates the rate expression
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Option_module
  use Reaction_Inhibition_Aux_module

  class(carbon_sandbox_rxn_mend_type) :: this
  PetscReal :: Residual(:)
  PetscReal :: Jacobian(:,:)
  type(option_type) :: option

  PetscInt, pointer :: specid(:)
  PetscReal, pointer :: stoich(:)
  PetscReal :: effective_rate
  PetscReal :: conc,inhibition,inhibition_factor
  PetscReal :: dummy
  PetscInt :: i, icomp, ncomp

  ncomp = this%reaction_equation%nspec
  specid => this%reaction_equation%specid
  stoich => this%reaction_equation%stoich

  effective_rate = log(this%rate_constant)
  inhibition = 0.d0
  do i = 1, ncomp
    if (stoich(i) > 0.d0) cycle
    icomp = specid(i)
    ! subtract due to negative stoichiometry
    effective_rate = effective_rate - stoich(i) * this%aux%ln_conc(icomp)
    conc = this%aux%conc(icomp) / this%aux%liter_water
    call ReactionInhibitionSmoothStep(conc,1.d-20,inhibition_factor,dummy)
    inhibition = inhibition + log(inhibition_factor)
  enddo
  effective_rate = exp(effective_rate+inhibition)
  do i = 1, ncomp
    icomp = specid(i)
    Residual(icomp) = Residual(icomp) - stoich(i)*effective_rate
  enddo

end subroutine CarbonRxnMEND20Evaluate

! ************************************************************************** !

recursive subroutine CarbonRxnMEND20Strip(this)
  !
  ! Destroys the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_rxn_mend_type) :: this

  call CarbonRxnBaseStrip(this)

end subroutine CarbonRxnMEND20Strip

end module Carbon_Sandbox_MEND20_class
