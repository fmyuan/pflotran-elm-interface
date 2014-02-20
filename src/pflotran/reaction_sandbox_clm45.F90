module Reaction_Sandbox_CLM45_class

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

                          ! 14.00674d0 / 12.011d0
  PetscReal, parameter :: CN_ratio_mass_to_mol = 1.16616d0 
  PetscInt, parameter :: CARBON_INDEX = 1
  PetscInt, parameter :: NITROGEN_INDEX = 2
  PetscInt, parameter :: SOM_INDEX = 1
  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_CLM4 = 1 
  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_Q10 = 2 

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_clm45_type

    PetscInt :: temperature_response_function
    PetscReal :: Q10
    PetscReal :: half_saturation_o2
    PetscReal :: half_saturation_nh3
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh3_no3
    PetscReal :: k_nitr_max    ! nitrification rate

    PetscInt :: nrxn
    PetscInt :: npool
    PetscReal, pointer :: CN_ratio(:)
    PetscReal, pointer :: rate_constant(:)
    PetscReal, pointer :: respiration_fraction(:)
    PetscReal, pointer :: inhibition_constant(:)
    PetscInt,  pointer :: upstream_pool_id(:)
    PetscInt,  pointer :: downstream_pool_id(:)
    PetscInt,  pointer :: pool_id_to_species_id(:,:)

    PetscInt :: species_id_co2
    PetscInt :: species_id_nh3
    PetscInt :: species_id_no3
    PetscInt :: species_id_n2o
    PetscInt :: species_id_o2
    PetscInt :: species_id_n2
    PetscInt :: species_id_nplant

    type(pool_type), pointer :: pools
    type(clm45_reaction_type), pointer :: reactions
  contains
    procedure, public :: ReadInput => CLM45_Read
    procedure, public :: Setup => CLM45_Setup
    procedure, public :: Evaluate => CLM45_React
    procedure, public :: Destroy => CLM45_Destroy
  end type reaction_sandbox_clm45_type
  
  type :: pool_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: CN_ratio
    type(pool_type), pointer :: next
  end type pool_type
  
  type :: clm45_reaction_type
    character(len=MAXWORDLENGTH) :: upstream_pool_name
    character(len=MAXWORDLENGTH) :: downstream_pool_name
    PetscReal :: rate_constant
    PetscReal :: respiration_fraction
    PetscReal :: inhibition_constant
    type(clm45_reaction_type), pointer :: next
  end type clm45_reaction_type
  
  public :: CLM45_Create

contains

! ************************************************************************** !

function CLM45_Create()
  ! 
  ! Allocates CLM45 reaction object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  implicit none
  
  type(reaction_sandbox_clm45_type), pointer :: CLM45_Create
  
  allocate(CLM45_Create)

  CLM45_Create%temperature_response_function=TEMPERATURE_RESPONSE_FUNCTION_CLM4
  CLM45_Create%Q10 = 1.5d0
  CLM45_Create%half_saturation_o2 = 1.0d-20
  CLM45_Create%half_saturation_nh3 = 1.0d-20
  CLM45_Create%half_saturation_no3 = 1.0d-20
  CLM45_Create%inhibition_nh3_no3 = 1.0d-10 
  CLM45_Create%k_nitr_max = 0.1d0/86400.0d0   ! nitrification rate

  CLM45_Create%nrxn = 0
  CLM45_Create%npool = 0
  nullify(CLM45_Create%CN_ratio)
  nullify(CLM45_Create%rate_constant)
  nullify(CLM45_Create%respiration_fraction)
  nullify(CLM45_Create%inhibition_constant)
  nullify(CLM45_Create%pool_id_to_species_id)
  nullify(CLM45_Create%upstream_pool_id)
  nullify(CLM45_Create%downstream_pool_id)
  CLM45_Create%species_id_co2 = 0
  CLM45_Create%species_id_nh3 = 0
  CLM45_Create%species_id_no3 = 0
  CLM45_Create%species_id_n2o = 0
  CLM45_Create%species_id_o2 = 0
  CLM45_Create%species_id_n2 = 0
  CLM45_Create%species_id_nplant = 0
  nullify(CLM45_Create%next)
  nullify(CLM45_Create%pools)
  nullify(CLM45_Create%reactions)

end function CLM45_Create

! ************************************************************************** !

subroutine CLM45_Read(this,input,option)
  ! 
  ! Reads input deck for reaction sandbox parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_clm45_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  
  type(pool_type), pointer :: new_pool, prev_pool
  type(clm45_reaction_type), pointer :: new_reaction, prev_reaction
  
  PetscReal :: rate_constant, turnover_time
  
  nullify(new_pool)
  nullify(prev_pool)
  
  nullify(new_reaction)
  nullify(prev_reaction)
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CLM45')
    call StringToUpper(word)   

    select case(trim(word))
      case('TEMPERATURE_RESPONSE_FUNCTION')
        do
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit

         call InputReadWord(input,option,word,PETSC_TRUE)
         call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CLM45P,TEMPERATURE RESPONSE FUNCTION')
         call StringToUpper(word)   

            select case(trim(word))
              case('CLM4')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4    
              case('Q10')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_Q10    
                  call InputReadDouble(input,option,this%Q10)  
                  call InputErrorMsg(input,option,'Q10', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM45P,TEMPERATURE RESPONSE FUNCTION')
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM45P,TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
         enddo 

     case('AMMONIA_INHIBITION_COEF')
         call InputReadDouble(input,option,this%inhibition_nh3_no3)
         call InputErrorMsg(input,option,'ammonia inhibition coefficient', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM45,REACTION')

     case('NITRIFICATION_RATE_COEF')
         call InputReadDouble(input,option,this%k_nitr_max)
         call InputErrorMsg(input,option,'nitrification rate coefficient', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM45,REACTION')
     case('POOLS')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit   

          allocate(new_pool)
          new_pool%name = ''
          new_pool%CN_ratio = 0.d0
          nullify(new_pool%next)

          call InputReadWord(input,option,new_pool%name,PETSC_TRUE)
          call InputErrorMsg(input,option,'pool name', &
            'CHEMISTRY,REACTION_SANDBOX,CLM45,POOLS')
          call InputReadDouble(input,option,new_pool%CN_ratio)
          if (InputError(input)) then
            new_pool%CN_ratio = -999.d0
          else
            ! convert CN ratio from mass C/mass N to mol C/mol N
            new_pool%CN_ratio = new_pool%CN_ratio * CN_ratio_mass_to_mol
          endif
          if (associated(this%pools)) then
            prev_pool%next => new_pool
          else
            this%pools => new_pool
          endif
          prev_pool => new_pool
          nullify(new_pool)
        enddo
      case('REACTION')
      
        allocate(new_reaction)
        new_reaction%upstream_pool_name = ''
        new_reaction%downstream_pool_name = ''
        new_reaction%rate_constant = -999.d0
        new_reaction%respiration_fraction = -999.d0
        new_reaction%inhibition_constant = 0.d0
        nullify(new_reaction%next)
        
        ! need to set these temporarily in order to check that they
        ! are not both set.
        turnover_time = 0.d0
        rate_constant = 0.d0
        
        do 
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,REACTION_SANDBOX,CLM45,REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('UPSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%upstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'upstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM45,REACTION')
            case('DOWNSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%downstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'downstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM45,REACTION')
            case('RATE_CONSTANT')
              call InputReadDouble(input,option,rate_constant)
              call InputErrorMsg(input,option,'rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM45,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLM45 RATE CONSTANT UNITS'
                call InputDefaultMsg(input,option)
              else              
                rate_constant = rate_constant * &
                  UnitsConvertToInternal(word,option)
              endif
            case('TURNOVER_TIME')
              call InputReadDouble(input,option,turnover_time)
              call InputErrorMsg(input,option,'turnover time', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM45,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLM45 TURNOVER TIME UNITS'
                call InputDefaultMsg(input,option)
              else              
                turnover_time = turnover_time * &
                  UnitsConvertToInternal(word,option)
              endif
            case('RESPIRATION_FRACTION')
              call InputReadDouble(input,option,new_reaction%respiration_fraction)
              call InputErrorMsg(input,option,'respiration fraction', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM45,REACTION')
            case('N_INHIBITION')
              call InputReadDouble(input,option,new_reaction%inhibition_constant)
              call InputErrorMsg(input,option,'inhibition constant', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM45,REACTION')
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM45,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
              call printErrMsg(option)
          end select
        enddo
        
        ! check to ensure that one of turnover time or rate constant is set.
        if (turnover_time > 0.d0 .and. rate_constant > 0.d0) then
          option%io_buffer = 'Only TURNOVER_TIME or RATE_CONSTANT may ' // &
            'be included in a CLM45 reaction definition, but not both. ' // &
            'See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        else if (turnover_time > 0.d0) then
          new_reaction%rate_constant = 1.d0 / turnover_time
        else
          new_reaction%rate_constant = rate_constant
        endif
        ! ensure that respiration fraction is 0-1.
        if (new_reaction%respiration_fraction < 0.d0 .or. &
                  new_reaction%respiration_fraction > 1.d0) then
          option%io_buffer = 'Respiratory fraction (rf) must be between ' // &
            'zero and one (i.e. 0. <= rf <= 1.) in a CLM45 reaction ' // &
            'definition. See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        endif
        ! If no downstream pool exists, ensure that respiration fraction = 1
        if (len_trim(new_reaction%downstream_pool_name) < 1 .and. &
            (1.d0 - new_reaction%respiration_fraction) > 1.d-40) then
          option%io_buffer = 'Respiratory fraction (rf) must be set to ' // &
            '1.0 if no downstream pool is specified in a CLM45 reaction ' // &
            'definition. See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        endif
        if (associated(this%reactions)) then
          prev_reaction%next => new_reaction
        else
          this%reactions => new_reaction
        endif
        prev_reaction => new_reaction
        nullify(new_reaction)        
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM45 keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine CLM45_Read

! ************************************************************************** !

subroutine CLM45_Setup(this,reaction,option)
  ! 
  ! Sets up CLM45 reaction after it has been read from input
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  use Reaction_Aux_module
  use Option_module
  use String_module
  use Immobile_Aux_module
  
  implicit none

  class(reaction_sandbox_clm45_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  
  character(len=MAXWORDLENGTH), allocatable :: pool_names(:)
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt :: icount

  type(pool_type), pointer :: cur_pool
  type(clm45_reaction_type), pointer :: cur_rxn
  
  ! count # pools
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    cur_pool => cur_pool%next
  enddo
  this%npool = icount
  
  ! count # reactions
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1
    cur_rxn => cur_rxn%next
  enddo
  this%nrxn = icount
  
  ! allocate and initialize arrays
  allocate(this%CN_ratio(this%npool))
  allocate(this%pool_id_to_species_id(0:2,this%npool))
  allocate(this%upstream_pool_id(this%nrxn))
  allocate(this%downstream_pool_id(this%nrxn))
  allocate(this%rate_constant(this%nrxn))
  allocate(this%respiration_fraction(this%nrxn))
  allocate(this%inhibition_constant(this%nrxn))
  this%CN_ratio = 0.d0
  this%pool_id_to_species_id = 0
  this%upstream_pool_id = 0
  this%downstream_pool_id = 0
  this%rate_constant = 0.d0
  this%respiration_fraction = 0.d0
  this%inhibition_constant = 0.d0
  
  ! temporary array for mapping pools in reactions
  allocate(pool_names(this%npool))
  
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    this%CN_ratio(icount) = cur_pool%CN_ratio
    pool_names(icount) = cur_pool%name
    if (cur_pool%CN_ratio < 0.d0) then
      ! Since no CN ratio provided, must provide two species with the
      ! same name as the pool with C or N appended.
      word = trim(cur_pool%name) // 'C'
      this%pool_id_to_species_id(CARBON_INDEX,icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      word = trim(cur_pool%name) // 'N'
      this%pool_id_to_species_id(NITROGEN_INDEX,icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      this%pool_id_to_species_id(0,icount) = 2
      if (minval(this%pool_id_to_species_id(:,icount)) <= 0) then
        option%io_buffer = 'For CLM45 pools with no CN ratio defined, ' // &
          'the user must define two immobile species with the same root ' // &
          'name as the pool with "C" or "N" appended, respectively.'
        call printErrMsg(option)
      endif
    else ! only one species (e.g. SOMX)
      this%pool_id_to_species_id(SOM_INDEX,icount) = &
        GetImmobileSpeciesIDFromName(cur_pool%name,reaction%immobile, &
                                     PETSC_TRUE,option)
      this%pool_id_to_species_id(0,icount) = 1
    endif
    cur_pool => cur_pool%next
  enddo
  
  ! map C and N species (solid phase for now)
  word = 'CO2(aq)'
  this%species_id_co2 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)

  word = 'Ammonia'
  this%species_id_nh3 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)

  word = 'Nitrate'
  this%species_id_no3 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)
  
  word = 'Nitrousoxide'
  this%species_id_n2o = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)
  
  word = 'O2(aq)'
  this%species_id_o2 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)

  word = 'N2(aq)'
  this%species_id_n2 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)

  word = 'PlantN'
  this%species_id_nplant = GetImmobileSpeciesIDFromName(word, &
                               reaction%immobile, PETSC_FALSE,option)

  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1
    this%upstream_pool_id(icount) = &
      StringFindEntryInList(cur_rxn%upstream_pool_name,pool_names)
    if (this%upstream_pool_id(icount) == 0) then
      option%io_buffer = 'Upstream pool "' // &
        trim(cur_rxn%upstream_pool_name) // &
        '" in reaction with downstream pool "' // &
        trim(cur_rxn%downstream_pool_name) // '" not found in list of pools.'
      call printErrMsg(option)
    endif
    if (len_trim(cur_rxn%downstream_pool_name) > 0) then
      this%downstream_pool_id(icount) = &
        StringFindEntryInList(cur_rxn%downstream_pool_name,pool_names)
      if (this%downstream_pool_id(icount) == 0) then
        option%io_buffer = 'Downstream pool "' // &
          trim(cur_rxn%downstream_pool_name) // &
          '" in reaction with upstream pool "' // &
          trim(cur_rxn%upstream_pool_name) // '" not found in list of pools.'
        call printErrMsg(option)
      endif
      if (this%CN_ratio(this%downstream_pool_id(icount)) < 0.d0) then
        option%io_buffer = 'For CLM45 reactions, downstream pools ' // &
          'must have a constant C:N ratio (i.e. C and N are not tracked ' // &
          ' individually.  Therefore, pool "' // &
          trim(cur_rxn%downstream_pool_name) // &
          '" may not be used as a downstream pool.'
        call printErrMsg(option)
      endif
    endif
    this%rate_constant(icount) = cur_rxn%rate_constant
    this%respiration_fraction(icount) = cur_rxn%respiration_fraction
    this%inhibition_constant(icount) = cur_rxn%inhibition_constant
    cur_rxn => cur_rxn%next
  enddo 
  
  deallocate(pool_names)
  
end subroutine CLM45_Setup

! ************************************************************************** !

subroutine CLM45_React(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                        global_auxvar,porosity,volume,reaction,option, local_id)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  use Option_module
  use Reaction_Aux_module

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif

  implicit none

#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif

  
  class(reaction_sandbox_clm45_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  
! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  PetscReal :: saturation, tc, theta

  PetscInt, parameter :: iphase = 1
  PetscInt :: ipool_up, ipool_down
  PetscInt :: ispec_d
  PetscInt :: ispec_uc, ispec_un
  PetscInt :: ires_d, ires_co2, ires_nh3, ires_o2
  PetscInt :: ires_uc, ires_un, ires_no3, ires_nplant
  PetscReal :: scaled_rate_const, rate

  PetscReal :: drate_uc   ! d Rate / d upstream c
  PetscReal :: drate_nh3  ! d Rate / d nh3 ammonia limitation
  PetscReal :: drate_o2   ! d Rate / d o2  oxygen limitation
  PetscReal :: d_nh3  ! d inhibition / d nh3 ammonia limitation
  PetscReal :: d_o2   ! d inhibition / d o2  oxygen limitation

! for N immobilization reactions with NO3 as N source 
  PetscReal :: rate_no3 
  PetscReal :: d_no3         ! d inhibition / d no3 = half_saturation/(c + half_saturation)^2 
  PetscReal :: drate_no3     ! d Rate_no3 / d no3 
  PetscReal :: drate_o2_no3  ! d Rate_no3 / d o2 
  PetscReal :: drate_uc_no3  ! d Rate_no3 / d uc 

  PetscInt :: i, j
  PetscInt :: irxn
  PetscInt :: local_id
  PetscErrorCode :: ierr

  PetscReal :: f_t    ! temperature response function
  PetscReal :: f_w    ! moisture response function
  PetscReal :: f_o2   ! o2 limitation
  PetscReal :: f_nh3  ! ammonia limitation
  PetscReal :: f_no3  ! ammonia limitation
  PetscReal :: f_nh3_inhibit  ! when ammonia available, immobilization pefers ammonia to nitrate

  PetscReal :: CN_ratio_up, CN_ratio_down
  PetscBool :: constant_CN_ratio_up
  PetscReal :: resp_frac
  PetscReal :: stoich_mineraln
  PetscReal :: stoich_co2
  PetscReal :: stoich_dc
  PetscReal :: stoich_uc, stoich_un
  
! /i PetscReal :: N_inhibition, d_N_inhibition
!  PetscReal :: drate_dN_inhibition
!  PetscBool :: use_N_inhibition
  PetscReal :: temp_real
  
  PetscReal :: c_nh3, c_no3, c_o2, c_uc, c_un
!  PetscReal :: o_inhibition, d_o_inhibition

! nitrification, denitrification 
  PetscReal :: respiration_rate ! CO2 production / O2 uptake rate from decomposition 
  PetscReal :: d_respiration_rate_dno3
  PetscReal :: f_ph, rate_nitri, drate_nitri, anaerobic_frac

  PetscReal :: diffus
  PetscReal :: soil_bulkdensity  ! clm4.5.35 denitrification, dry bulk density + water
  PetscReal :: bulkdensity_dry_cell

  PetscScalar, pointer :: bulkdensity_dry(:)
  PetscReal :: ratio_k1, ratio_no3_co2
  PetscReal :: M_to_ug_per_gsoil, kg_water_m3
  PetscReal :: fr_WFPS, wfps_vr
  PetscReal :: n2_n2o_ratio, stoich_n2o, stoich_n2
  PetscReal :: d_stoich_n2o_d_no3, d_stoich_n2_d_no3
  PetscReal :: rate_deni_resp, rate_deni_no3, rate_deni
  PetscReal :: drate_deni_resp, drate_deni
  PetscReal, parameter :: N_molecular_weight = 14.0067d0
  PetscReal, parameter :: nitrif_n2o_loss_frac = 6.d-4 !fraction of N lost as N2O in nitrification (Li et al., 2000)
!  !  real(r8), parameter :: nitrif_n2o_loss_frac = 0.02_r8   !fraction of N lost as N2O in nitrification (Parton et al., 2001)
!  real(r8), parameter :: nitrif_n2o_loss_frac = 6.e-4_r8   !fraction of N lost as N2O in nitrification (Li et al., 2000)
  PetscInt :: ires_n2, ires_n2o

  PetscReal :: rate_nplant, drate_nplant
  PetscReal :: rate_nplant_no3, drate_nplant_no3

  PetscScalar, pointer :: rate_plantnuptake_pf_loc(:)   !

  c_nh3     = rt_auxvar%total(this%species_id_nh3, iphase)
  temp_real = c_nh3 + this%half_saturation_nh3
  f_nh3     = c_nh3 / temp_real 
  d_nh3     = this%half_saturation_nh3 / temp_real / temp_real

  if (this%species_id_o2 > 0) then
      c_o2 = rt_auxvar%total(this%species_id_o2, iphase)
      temp_real = c_o2 + this%half_saturation_o2
      f_o2 = c_o2 / temp_real 
      d_o2 = this%half_saturation_o2 / temp_real / temp_real 
  endif 

  if (this%species_id_no3 > 0) then
      c_no3 = rt_auxvar%total(this%species_id_no3, iphase)
      temp_real = c_no3 + this%half_saturation_no3
      f_no3 = c_no3 / temp_real 
      d_no3 = this%half_saturation_no3 / temp_real / temp_real 

      temp_real = this%inhibition_nh3_no3 + c_nh3 
      f_nh3_inhibit = this%inhibition_nh3_no3/temp_real
  endif 

  ! indices for C and N species
  ires_co2 = this%species_id_co2
  ires_nh3 = this%species_id_nh3
  ires_o2  = this%species_id_o2
  ires_no3 = this%species_id_no3
  ires_n2  = this%species_id_n2
  ires_n2o = this%species_id_n2o
  ires_nplant = this%species_id_nplant + reaction%offset_immobile 

! plant n uptake
#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%rate_plantnuptake_pf, rate_plantnuptake_pf_loc, ierr)

  rate_nplant = rate_plantnuptake_pf_loc(local_id) * volume ! mol/m3/s * m3

  call VecRestoreArrayReadF90(clm_pf_idata%rate_plantnuptake_pf, rate_plantnuptake_pf_loc, ierr)

  rate_nplant_no3 = rate_nplant * f_nh3_inhibit
  drate_nplant = rate_nplant * d_nh3
  rate_nplant = rate_nplant * f_nh3 

  Residual(ires_nh3) = Residual(ires_nh3) + rate_nplant
  Residual(ires_nplant) = Residual(ires_nplant) - rate_nplant

  if (compute_derivative) then 
     Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) + drate_nplant
     Jacobian(ires_nplant,ires_nh3) = Jacobian(ires_nplant,ires_nh3) - drate_nplant
  endif

  if(this%species_id_no3 > 0) then
    drate_nplant_no3 = rate_nplant_no3 *d_no3
    rate_nplant_no3  = rate_nplant_no3 *f_no3

    Residual(ires_no3) = Residual(ires_no3) + rate_nplant_no3
    Residual(ires_nplant) = Residual(ires_nplant) - rate_nplant_no3

    if (compute_derivative) then 
     Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) + drate_nplant_no3
     Jacobian(ires_nplant,ires_no3) = Jacobian(ires_nplant,ires_no3) - drate_nplant_no3
    endif
  endif
#endif

 ! temperature response function 
  tc = global_auxvar%temp(1) 
  f_t = GetTemperatureResponse(tc, this%temperature_response_function, this%Q10) 

  saturation = global_auxvar%sat(1)
  theta = saturation * porosity 

  ! moisture response function 
  f_w = GetMoistureResponse(theta, local_id)

  if(f_t < 1.0d-20 .or. f_w < 1.0d-20) then
     return
  endif

! decomposition reactions
  respiration_rate = 0.d0
  d_respiration_rate_dno3 = 0.d0

  do irxn = 1, this%nrxn
  
    ! scaled_rate_const units: (m^3 bulk / s) = (1/s) * (m^3 bulk)
    scaled_rate_const = this%rate_constant(irxn)*volume*f_t*f_w

    resp_frac = this%respiration_fraction(irxn)

    ! upstream pool    
    ipool_up = this%upstream_pool_id(irxn)
    constant_CN_ratio_up = (this%pool_id_to_species_id(0,ipool_up) == 1)

    stoich_uc = 1.d0
    if (.not. constant_CN_ratio_up) then
      ! upstream pool is Litter pool with two species (C,N)
      !
      !  a LitC_i + b LitN_i -> c SOM_i + d C + e N
      !
      ispec_uc = this%pool_id_to_species_id(CARBON_INDEX,ipool_up)
      ispec_un = this%pool_id_to_species_id(NITROGEN_INDEX,ipool_up)
      c_uc = rt_auxvar%immobile(ispec_uc)
      c_un = rt_auxvar%immobile(ispec_un)
      stoich_un = c_un / c_uc
    else
      ! upstream pool is an SOM pool with one species
      !
      !  a SOM_i -> c SOM_(i+1) + d C + e N
      !
      ispec_uc    = this%pool_id_to_species_id(SOM_INDEX,ipool_up)
      c_uc        = rt_auxvar%immobile(ispec_uc)
      CN_ratio_up = this%CN_ratio(ipool_up)
      stoich_un   = stoich_uc / CN_ratio_up
    endif

    ! downstream pool
    ipool_down = this%downstream_pool_id(irxn)
    ! a downstream pool need not exist if last in succession and
    ! respiration fraction = 1.
    if (ipool_down > 0) then
      ! downstream pool is always SOM (i.e. not split between C and N
      ! as a litter would be).
      ispec_d = this%pool_id_to_species_id(SOM_INDEX,ipool_down)
      CN_ratio_down = this%CN_ratio(ipool_down)
      ! c = (1-resp_frac) * a
      stoich_dc = (1.d0-resp_frac) * stoich_uc
    else    
      ispec_d = 0
      stoich_dc = 0.d0
      CN_ratio_down = 1.d0 ! to prevent divide by zero below.
    endif
      
    ! d = resp_frac * a
    stoich_co2 = resp_frac * stoich_uc

    ! e = b - c / CN_ratio_dn
    stoich_mineraln = stoich_un - stoich_dc / CN_ratio_down
 
    ! residual units: (mol/sec) = (m^3 bulk/s) * (mol/m^3 bulk)
    rate     = c_uc * scaled_rate_const
    drate_uc = scaled_rate_const 

!   NH3 limiting
    if(stoich_mineraln < 0.0d0) then
       drate_nh3 = rate * d_nh3 
       rate      = rate * f_nh3
       drate_uc  = drate_uc * f_nh3
    endif 

!   O2 limiting
    if (this%species_id_o2 > 0) then
       drate_o2 = rate * d_o2
       rate     = rate * f_o2
       drate_nh3 = drate_nh3 * f_o2 
       drate_uc  = drate_uc * f_o2
    endif 
    
    ! calculation of residual
    ! carbon
    Residual(ires_co2) = Residual(ires_co2) - stoich_co2 * rate
    
    respiration_rate = respiration_rate + stoich_co2 * rate
 
    ! nitrogen
    Residual(ires_nh3) = Residual(ires_nh3) - stoich_mineraln * rate

    ! C species in upstream pool (litter or SOM)
    ires_uc = reaction%offset_immobile + ispec_uc
    ! scaled by negative one since it is a reactant 
    Residual(ires_uc) = Residual(ires_uc) - (-1.d0) * stoich_uc * rate
    
    ! N species in upstream pool (litter only)
    if (.not.constant_CN_ratio_up) then
      ! N species in upstream pool
      ires_un = reaction%offset_immobile + ispec_un
      ! scaled by negative one since it is a reactant 
      Residual(ires_un) = Residual(ires_un) - (-1.d0) * stoich_un * rate
    endif
    
    if (ispec_d > 0) then
      ! downstream pool
      ires_d = reaction%offset_immobile + ispec_d
      Residual(ires_d) = Residual(ires_d) - stoich_dc * rate
    endif
    
    if (this%species_id_o2 > 0) then
      Residual(ires_o2) = Residual(ires_o2) - (-1.d0) * stoich_co2 * rate
    endif

!   start residual calculation for N immobilization reaction with NO3 uptake
!   if nitrate is available, N immobilization decomposition reactions occurs
!   with rate depending on NH3, with reduced rate if NH3 is abundent  
    if(this%species_id_no3 > 0 .and. stoich_mineraln < 0.d0) then
       ires_no3 = this%species_id_no3
       rate_no3     = c_uc * scaled_rate_const * f_nh3_inhibit
       drate_uc_no3 = scaled_rate_const * f_nh3_inhibit

!      NO3 limiting
       if(stoich_mineraln < 0.0d0) then   ! same stoichiometric coefficient as nh3
          drate_no3     = rate_no3 * d_no3 
          rate_no3      = rate_no3 * f_no3
          drate_uc_no3  = drate_uc_no3 * f_no3
       endif 

!      O2 limiting
       if (this%species_id_o2 > 0) then
          drate_o2_no3 = rate_no3 * d_o2
          rate_no3     = rate_no3 * f_o2
          drate_no3    = drate_no3 * f_o2 
          drate_uc_no3 = drate_uc_no3 * f_o2
       endif 
    
    ! residuals
    ! carbon
       Residual(ires_co2) = Residual(ires_co2) - stoich_co2 * rate_no3
    
       respiration_rate = respiration_rate + stoich_co2 * rate_no3

       d_respiration_rate_dno3 = d_respiration_rate_dno3 + stoich_co2 * drate_no3

    ! nitrogen
       Residual(ires_no3) = Residual(ires_no3) - stoich_mineraln * rate_no3

    ! C species in upstream pool (litter or SOM)
       Residual(ires_uc) = Residual(ires_uc) - (-1.d0) * stoich_uc * rate_no3
    
    ! N species in upstream pool (litter only)
       if (.not.constant_CN_ratio_up) then
         Residual(ires_un) = Residual(ires_un) - (-1.d0) * stoich_un * rate_no3
       endif

    ! downstream pool 
       if (ispec_d > 0) then
         Residual(ires_d) = Residual(ires_d) - stoich_dc * rate_no3
       endif
    
       if (this%species_id_o2 > 0) then
         Residual(ires_o2) = Residual(ires_o2) - (-1.d0) * stoich_co2 * rate_no3
       endif
    endif 
!   end residual calculation for N immobilization reaction with NO3 uptake

    if (compute_derivative) then
! with respect to upstream c    
      ! carbon
      Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
        stoich_co2 * drate_uc

      ! nitrogen
         ! dR_N/dC_u = d nR/dC_u = dn/dC_u R + n dR/dC_u
         ! first, n dR/dC_u
      Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - &
        stoich_mineraln * drate_uc

      ! upstream N pool
      if (.not.constant_CN_ratio_up) then
         ! second, dn/dC_u R, n = u - (1 - f)d, dn/dC_u = du/dC_u = -N_u/C_u^2 = -u/C_u
         ! dn/dC_u R = -u R/C_u = -u rate_uc 
        Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - &
          (-1.d0) *stoich_un * drate_uc
      endif

      ! upstream C pool
      Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) - &
        (-1.d0) * stoich_uc * drate_uc
        ! scaled by negative one since it is a reactant 

      ! upstream N pool
      if (.not.constant_CN_ratio_up) then
        ! R_Nu = Nu/Cu * R_Cu, dR_Nu/dCu = Nu/Cu dR_Cu/dCu - Nu/Cu^2 R_Cu
        ! R_Cu = k Cu,         dR_Cu/dCu = k
        ! dR_Nu/dCu = Nu/Cu k - Nu/Cu kCu/Cu = 0
!        Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) - &
!          (-1.d0) * stoich_un * drate_uc 
      endif

      ! downstream pool
      if (ispec_d > 0) then
        Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - &
          stoich_dc * drate_uc
      endif

      if (this%species_id_o2 > 0) then
        Jacobian(ires_o2,ires_uc) = Jacobian(ires_o2,ires_uc) -  &
          (-1.d0) * stoich_co2 * drate_uc
      endif

!with respect to upstream n
      if (.not.constant_CN_ratio_up) then
        ! upstream n, dR_Nu/dNu = d uR/dNu = R/Cu = dR/dCu
         Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) -(-1.d0)*drate_uc

        ! mineral N, dR_N/dNu = d nR/du = Rd (u - (1-f)d)/dNu = dR/dCu 
         Jacobian(ires_nh3,ires_un) = Jacobian(ires_nh3,ires_un) - drate_uc
      endif

!with respect to O2
      if (this%species_id_o2 > 0) then
        ! carbon
        Jacobian(ires_co2,ires_o2) = Jacobian(ires_co2,ires_o2) - &
          stoich_co2 * drate_o2 * &
          rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_o2,iphase)
  
        ! nitrogen
        Jacobian(ires_nh3,ires_o2) = Jacobian(ires_nh3,ires_o2) - &
          stoich_mineraln * drate_o2 * &
          rt_auxvar%aqueous%dtotal(this%species_id_nh3,this%species_id_o2,iphase)
  
        ! upstream C pool
        Jacobian(ires_uc,ires_o2) = Jacobian(ires_uc,ires_o2) - &
          (-1.d0) * stoich_uc * drate_o2
  
        ! upstream N pool
        if (.not.constant_CN_ratio_up) then
          Jacobian(ires_un,ires_o2) = Jacobian(ires_un,ires_o2) - &
            (-1.d0) * stoich_un * drate_o2 
        endif
  
        ! downstream pool
        if (ispec_d > 0) then
          Jacobian(ires_d,ires_o2) = Jacobian(ires_d,ires_o2) - &
            stoich_dc * drate_o2
        endif

        Jacobian(ires_o2,ires_o2) = Jacobian(ires_o2,ires_o2) - & 
          (-1.0) * stoich_co2 * drate_o2 * &
          rt_auxvar%aqueous%dtotal(this%species_id_o2,this%species_id_o2,iphase)
      endif

!with respect to nh3
      if(stoich_mineraln < 0.0d0) then
        ! carbon
        Jacobian(ires_co2,ires_nh3) = Jacobian(ires_co2,ires_nh3) - &
          stoich_co2 * drate_nh3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_nh3,iphase)
  
        ! nitrogen
        Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) - &
          stoich_mineraln * drate_nh3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_nh3,this%species_id_nh3,iphase)
  
        ! upstream C pool
        Jacobian(ires_uc,ires_nh3) = Jacobian(ires_uc,ires_nh3) - &
          (-1.d0) * stoich_uc * drate_nh3
  
        ! upstream N pool
        if (.not.constant_CN_ratio_up) then
          Jacobian(ires_un,ires_nh3) = Jacobian(ires_un,ires_nh3) - &
            (-1.d0) * stoich_un * drate_nh3 
        endif
  
        ! downstream pool
        if (ispec_d > 0) then
          Jacobian(ires_d,ires_nh3) = Jacobian(ires_d,ires_nh3) - &
            stoich_dc * drate_nh3
        endif

        Jacobian(ires_o2,ires_nh3) = Jacobian(ires_o2,ires_nh3) - &
          (-1.d0) * stoich_co2 * drate_nh3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_o2,this%species_id_nh3,iphase)
      endif

!   start jacobian calculation for N immobilization reaction with NO3 uptake
      if(this%species_id_no3 > 0 .and. stoich_mineraln < 0.d0) then
! with respect to upstream c    
      ! carbon
        Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
          stoich_co2 * drate_uc_no3

      ! nitrogen
        Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - &
          stoich_mineraln * drate_uc_no3

        if (.not.constant_CN_ratio_up) then
          Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - &
            (-1.d0) *stoich_un * drate_uc_no3
        endif

      ! upstream C pool
        Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) - &
          (-1.d0) * stoich_uc * drate_uc_no3

      ! upstream N pool
        if (.not.constant_CN_ratio_up) then
         ! 0 for first order rate
        endif

      ! downstream pool
        if (ispec_d > 0) then
          Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - &
            stoich_dc * drate_uc_no3
        endif

        if (this%species_id_o2 > 0) then
          Jacobian(ires_o2,ires_uc) = Jacobian(ires_o2,ires_uc) -  &
            (-1.d0) * stoich_co2 * drate_uc_no3
        endif

!with respect to upstream n
        if (.not.constant_CN_ratio_up) then
        ! upstream n, dR_Nu/dNu = d uR/dNu = R/Cu = dR/dCu
           Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) -(-1.d0)*drate_uc_no3

        ! mineral N, dR_N/dNu = d nR/du = Rd (u - (1-f)d)/dNu = dR/dCu 
           Jacobian(ires_no3,ires_un) = Jacobian(ires_no3,ires_un) - drate_uc_no3
        endif

!with respect to O2
        if (this%species_id_o2 > 0) then
        ! carbon
          Jacobian(ires_co2,ires_o2) = Jacobian(ires_co2,ires_o2) - &
            stoich_co2 * drate_o2_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_o2,iphase)
  
        ! nitrogen
          Jacobian(ires_no3,ires_o2) = Jacobian(ires_no3,ires_o2) - &
            stoich_mineraln * drate_o2_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_o2,iphase)
  
        ! upstream C pool
          Jacobian(ires_uc,ires_o2) = Jacobian(ires_uc,ires_o2) - &
            (-1.d0) * stoich_uc * drate_o2_no3
  
        ! upstream N pool
          if (.not.constant_CN_ratio_up) then
            Jacobian(ires_un,ires_o2) = Jacobian(ires_un,ires_o2) - &
              (-1.d0) * stoich_un * drate_o2_no3 
          endif
  
        ! downstream pool
          if (ispec_d > 0) then
            Jacobian(ires_d,ires_o2) = Jacobian(ires_d,ires_o2) - &
              stoich_dc * drate_o2_no3
          endif

          Jacobian(ires_o2,ires_o2) = Jacobian(ires_o2,ires_o2) - & 
            (-1.0) * stoich_co2 * drate_o2_no3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_o2,this%species_id_o2,iphase)
        endif

!with respect to no3
        if(stoich_mineraln < 0.0d0) then
        ! carbon
          Jacobian(ires_co2,ires_no3) = Jacobian(ires_co2,ires_no3) - &
            stoich_co2 * drate_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_no3,iphase)
  
        ! nitrogen
          Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) - &
            stoich_mineraln * drate_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_no3,iphase)
  
        ! upstream C pool
          Jacobian(ires_uc,ires_no3) = Jacobian(ires_uc,ires_no3) - &
            (-1.d0) * stoich_uc * drate_no3
  
        ! upstream N pool
          if (.not.constant_CN_ratio_up) then
            Jacobian(ires_un,ires_no3) = Jacobian(ires_un,ires_no3) - &
              (-1.d0) * stoich_un * drate_no3 
          endif
  
        ! downstream pool
          if (ispec_d > 0) then
            Jacobian(ires_d,ires_no3) = Jacobian(ires_d,ires_no3) - &
              stoich_dc * drate_no3
          endif

          Jacobian(ires_o2,ires_no3) = Jacobian(ires_o2,ires_no3) - &
            (-1.d0) * stoich_co2 * drate_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_o2,this%species_id_no3,iphase)
        endif
      endif
!     end jacobian calculation for N immobilization reaction with NO3 uptake
    endif
!   end of jacobian calculation
  enddo
   
!  do i = 1, 5
!     write(*, *) (Jacobian(i, j), j = 1, 5)
!  enddo
  if(this%species_id_no3 < 0) then
    return ! no no3, then no nitrification or denitrification
  endif
! nitrification   
! moisture function-- assume the same moisture function as limits heterotrophic respiration
! Parton et al. base their nitrification- soil moisture rate constants based on heterotrophic rates-- can we do the same?
! k_nitr_h2o_vr(c,j) = w_scalar(c,j)

  ! Set maximum nitrification rate constant 
!  k_nitr_max =  0.1 / 86400.0   ! [1/sec] 10%/day  Parton et al., 2001 
  ! Todo:  SPM - the explicit divide gives different results than when that
  ! value is placed in the parameters netcdf file.  To get bfb, keep the 
  ! divide in source.
  !k_nitr_max = CNNitrifDenitrifParamsInst%k_nitr_max

!  k_nitr_max =  0.1 / 86400.0   ! [1/sec] 10%/day  Parton et al., 2001 
! nitrification constant is a set scalar * temp, moisture, and ph scalars
!         k_nitr_vr(c,j) = k_nitr_max * k_nitr_t_vr(c,j) * k_nitr_h2o_vr(c,j) * k_nitr_ph_vr(c,j)

  call GetAnaerobicFraction(saturation, porosity, tc, f_t, f_w, local_id, &
         c_o2*theta*1.d3, respiration_rate/volume, anaerobic_frac, diffus)

  f_ph = GetpHResponse(6.5d0)

  rate_nitri = this%k_nitr_max * f_t * f_ph * f_w * (1.0 - anaerobic_frac) * c_nh3 ! do we need to * water_volume?  

!  if(this%half_saturation_nh3 > 1.0d-20) then
!    f = k N N/(s + N)
!    df/dN = k (2N(s + N) - N^2)/(s + N)^2 = k (2Ns + N^2)/(s + N)^2
!                                          = k N (2s + N)/(S + N)^2
  temp_real   = this%half_saturation_nh3 + c_nh3
  drate_nitri = rate_nitri  * (2.0 * this%half_saturation_nh3 + c_nh3) / temp_real / temp_real
  rate_nitri  = rate_nitri * c_nh3 / temp_real

!         ! limit to non-frozen soil layers
!         if ( t_soisno(c,j) <= SHR_CONST_TKFRZ .and. no_frozen_nitrif_denitrif) then
!            pot_f_nit_vr(c,j) = 0._r8
!         endif

  Residual(ires_nh3) = Residual(ires_nh3) - (-1.0) * rate_nitri
 
  if(this%species_id_n2o > 0) then
     Residual(ires_no3) = Residual(ires_no3) - rate_nitri * (1.0d0 - nitrif_n2o_loss_frac)
     Residual(ires_n2o) = Residual(ires_n2o) - rate_nitri * nitrif_n2o_loss_frac
  else
     Residual(ires_no3) = Residual(ires_no3) - rate_nitri
  endif
  
  if (compute_derivative) then
     Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) + drate_nitri * &
        rt_auxvar%aqueous%dtotal(this%species_id_nh3,this%species_id_nh3,iphase)

     if(this%species_id_n2o > 0) then
       Jacobian(ires_no3,ires_nh3) = Jacobian(ires_no3,ires_nh3) - drate_nitri * &
        (1.0d0 - nitrif_n2o_loss_frac) * &   
        rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_nh3,iphase)
     else
       Jacobian(ires_no3,ires_nh3) = Jacobian(ires_no3,ires_nh3) - drate_nitri * &
        nitrif_n2o_loss_frac * &
        rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_nh3,iphase)
     endif
  endif

  if(this%species_id_n2o < 0 .or. this%species_id_n2 < 0) then
    return ! no n2o or n2, then no denitrification
  endif
!---------------- denitrification
! follows CENTURY denitrification scheme (Parton et al., (2001, 1996))
#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%bulkdensity_dry_pf, bulkdensity_dry, ierr)
  bulkdensity_dry_cell = bulkdensity_dry(local_id) ! kg/m3
  call VecRestoreArrayReadF90(clm_pf_idata%bulkdensity_dry_pf, bulkdensity_dry, ierr)
#else
  bulkdensity_dry_cell = 1.0d0
#endif

  kg_water_m3 = porosity*global_auxvar%sat(iphase)*global_auxvar%den_kg(iphase)
  soil_bulkdensity = bulkdensity_dry_cell + kg_water_m3

  M_to_ug_per_gsoil = 1.d6 *theta / soil_bulkdensity * N_molecular_weight

!! maximum potential denitrification rates based on heterotrophic respiration rates or nitrate concentrations, 
!! from (del Grosso et al., 2000)
  rate_deni_no3 = 1.15d0 * (c_no3 * M_to_ug_per_gsoil)**0.57d0  ! ug/(g soil day)
  rate_deni_no3 = rate_deni_no3/M_to_ug_per_gsoil/86400.0d0       ! mole / L s
  rate_deni_no3 = rate_deni_no3 * theta * volume * 1.d3  ! mole/s
! rate_deni = 1.15 * M_to_ug_per_gsoil**(-0.43)/86400.0 * conc_no3**1.57/(k + c)
! d(x^n/(k+x))/dx = [nx^(n-1)(k+x) - x^n)]/(k + x)^2

  rate_deni_no3 = rate_deni_no3 * f_no3 * anaerobic_frac

  rate_deni_resp=0.1d0*(respiration_rate*M_to_ug_per_gsoil*86400.0d0)**1.3d0/M_to_ug_per_gsoil/86400.0d0

  rate_deni_resp = rate_deni_resp * anaerobic_frac

  if(rate_deni_no3 < rate_deni_resp) then  ! smaller rate
    rate_deni = rate_deni_no3
    temp_real = c_no3 + this%half_saturation_no3
    drate_deni = 1.15d0 * M_to_ug_per_gsoil**(-0.43d0)/86400.0d0 * &
              (1.57d0 * c_no3**0.57 * temp_real - c_no3**1.57d0)/temp_real/temp_real * &
              theta * volume * 1.d3
  else
    rate_deni = rate_deni_resp
! R = 0.1 (M_to_ug_per_gsoil*86400)^0.3 Rr^1.3 
! dR = 0.1 (M_to_ug_per_gsoil*86400)^0.3 1.3 * Rr^0.3 * dRr
    drate_deni = 0.1d0 * (M_to_ug_per_gsoil*86400.0d0)**0.3d0 * 1.3d0 * respiration_rate**0.3d0 * d_respiration_rate_dno3
  endif

! now calculate the ratio of N2O to N2 from denitrifictaion, following Del Grosso et al., 2000
! diffusivity constant (figure 6b)
  ratio_k1 = max(1.7d0, 38.4d0 - 350.0d0 * diffus)

! ratio function (figure 7c)
  if ( respiration_rate > 0.0 ) then
       ratio_no3_co2 = c_no3 / respiration_rate / 86400.0    !  not unitless?
  else
! fucntion saturates at large no3/co2 ratios, so set as some nominally large number
       ratio_no3_co2 = 100.0
  endif

! total water limitation function (Del Grosso et al., 2000, figure 7a)
!  wfps_vr = max(min(h2osoi_vol/porosity, 1.0), 0.0) * 100.0
  wfps_vr = max(min(saturation, 1.0d0), 0.0d0) * 100.0d0
  fr_WFPS = max(0.1d0, 0.015d0 * wfps_vr - 0.32d0)
!         if (use_lch4) then
!            if (anoxia_wtsat) then
!               fr_WFPS(c,j) = fr_WFPS(c,j)*(1._r8 - finundated(c)) + finundated(c)*1.18_r8
!            end if
!         end if

         ! final ratio expression 
  n2_n2o_ratio = ratio_k1 * max(0.16d0, exp(-0.8d0 * ratio_no3_co2)) * fr_WFPS
  stoich_n2 = n2_n2o_ratio/(1.0 + n2_n2o_ratio)
  stoich_n2o = 1.0 - stoich_n2

  if(respiration_rate > 0.0) then
     d_stoich_n2_d_no3 = 1.0 / (1.0 + n2_n2o_ratio) / (1.0 + n2_n2o_ratio) * &
            ratio_k1 * fr_WFPS * (-0.8 / respiration_rate) * exp (-0.8*ratio_no3_co2)
     d_stoich_n2o_d_no3 = 1.0 - d_stoich_n2_d_no3
  endif

! residual
  Residual(ires_no3) = Residual(ires_no3) - (-1.0) * rate_deni
  Residual(ires_n2o) = Residual(ires_n2o) - rate_deni * stoich_n2o
  Residual(ires_n2) = Residual(ires_n2) - rate_deni * stoich_n2
 
  if (compute_derivative) then
     Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) + drate_deni * &
        rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_no3,iphase)
      
     Jacobian(ires_n2o,ires_no3) = Jacobian(ires_n2o,ires_no3) + drate_deni * &
        stoich_n2o * &      
        rt_auxvar%aqueous%dtotal(this%species_id_n2o,this%species_id_no3,iphase)
      
     Jacobian(ires_n2,ires_no3) = Jacobian(ires_n2,ires_no3) + drate_deni * &
        stoich_n2 * &      
        rt_auxvar%aqueous%dtotal(this%species_id_n2,this%species_id_no3,iphase)
  endif

end subroutine CLM45_React

! ************************************************************************** !

subroutine CLM45_Destroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_clm45_type) :: this
  
  type(pool_type), pointer :: cur_pool, prev_pool
  type(clm45_reaction_type), pointer :: cur_reaction, prev_reaction
  
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    prev_pool => cur_pool
    cur_pool => cur_pool%next
    deallocate(prev_pool)
    nullify(prev_pool)
  enddo
  
  cur_reaction => this%reactions
  do
    if (.not.associated(cur_reaction)) exit
    prev_reaction => cur_reaction
    cur_reaction => cur_reaction%next
    deallocate(prev_reaction)
    nullify(prev_reaction)
  enddo
  
  call DeallocateArray(this%CN_ratio)
  call DeallocateArray(this%rate_constant)
  call DeallocateArray(this%respiration_fraction)
  call DeallocateArray(this%inhibition_constant)
  call DeallocateArray(this%upstream_pool_id)
  call DeallocateArray(this%downstream_pool_id)
  call DeallocateArray(this%pool_id_to_species_id)
  
end subroutine CLM45_Destroy

! ************************************************************************** !
! temperature response function 

function GetTemperatureResponse(tc, itype, Q10)

  implicit none
  
  PetscInt  :: itype
  PetscReal :: Ft, tc, Q10, tk

  PetscReal, parameter :: one_over_71_02 = 1.408054069d-2
  PetscReal :: GetTemperatureResponse

  select case(itype)
!     CLM4.5 temperature response function
      case(TEMPERATURE_RESPONSE_FUNCTION_Q10)
          Ft = Q10 ** ((tc - 25.0d0) / 10.0d0)

!     CLM-CN
!     Equation: F_t = exp(308.56*(1/17.02 - 1/(T - 227.13)))
      case(TEMPERATURE_RESPONSE_FUNCTION_CLM4) 
  
          tk = tc + 273.15d0

          if(tk > 227.15d0) then
              Ft = exp(308.56d0*(one_over_71_02 - 1.d0/(tk - 227.13d0)))
          else
              Ft = 0.d0
          endif
  end select
  GetTemperatureResponse = Ft 

end function GetTemperatureResponse
  
! ************************************************************************** !

Function GetMoistureResponse(theta, local_id)

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif
  
  implicit none
  
#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif

  PetscReal :: F_theta, theta
  PetscReal :: GetMoistureResponse
  PetscReal, parameter :: theta_min = 0.01d0     ! 1/nat log(0.01d0)
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1
  PetscReal, parameter :: twelve_over_14 = 0.857142857143d0

  PetscInt :: local_id
  PetscReal :: maxpsi, psi, tc
  PetscReal, parameter :: minpsi = -10.0d0  ! MPa

#ifdef CLM_PFLOTRAN
  PetscScalar, pointer :: sucsat_pf_loc(:)   !
  PetscScalar, pointer :: soilpsi_pf_loc(:)   !
#endif

  PetscErrorCode :: ierr


#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)
  call VecGetArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi_pf_loc, ierr)

  maxpsi = sucsat_pf_loc(local_id) * (-9.8d-6)
  psi = min(soilpsi_pf_loc(local_id), maxpsi)

  if(psi > minpsi) then
     F_theta = log(minpsi/psi)/log(minpsi/maxpsi)
  else
     F_theta = 0.0d0
     call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi_pf_loc, ierr)
  endif

  call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi_pf_loc, ierr)
#else

  ! inhibition due to moisture content
  ! Equation: F_theta = log(theta_min/theta) / log(theta_min/theta_max)
  ! Assumptions: theta is saturation
  !              theta_min = 0.01, theta_max = 1.
  F_theta = 1.0d0 !log(theta_min/max(theta_min,theta)) * one_over_log_theta_min 
  
#endif
  GetMoistureResponse = F_theta

end function GetMoistureResponse

! ************************************************************************** !
subroutine GetAnaerobicFraction(saturation, porosity, tc, f_t, f_w, local_id, conc_o2, rate_o2, anaerobic_frac, diffus)

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif
  
  implicit none

#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif
  
  PetscReal :: volume, porosity, saturation
  PetscInt :: local_id, ires_nh3, ires_no3
  PetscErrorCode     :: ierr
  PetscReal :: theta, tc, tk, conc_o2, rate_o2

  ! inhibition variables
  PetscReal :: f_t
  PetscReal :: f_w
  PetscReal :: f_ph
  PetscReal :: tmp_real 
  PetscReal :: anaerobic_frac, diffus
!  PetscReal :: GetAnaerobicFraction
  PetscReal, parameter :: WT_Saturation = 0.95 !volumetric soil water defining top of water table clm4.5

#ifdef CLM_PFLOTRAN
  PetscInt, parameter :: iphase = 1
  PetscReal, parameter :: spval = 1.0d36
 
  PetscReal :: surface_tension_water
  PetscReal :: rij_kro_a                  !  Arah and Vinten 1995
  PetscReal :: rij_kro_alpha              !  Arah and Vinten 1995
  PetscReal :: rij_kro_beta               !  Arah and Vinten 1995
  PetscReal :: rij_kro_gamma              !  Arah and Vinten 1995
  PetscReal :: rij_kro_delta              !  Arah and Vinten 1995
  PetscReal :: organic_max                ! organic matter content (kg/m3) where soil is assumed to act like peat
  PetscReal :: pH                         ! placeholder
  PetscReal :: co2diff_con1, co2diff_con2 ! diffusion constants for CO2
 
  PetscReal :: d_con_w_1_o2, d_con_w_2_o2, d_con_w_3_o2           ! water diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-4)
  PetscReal :: d_con_g_1_o2, d_con_g_2_o2           ! gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)

  PetscReal :: f_a   !
  PetscReal :: e_a   ! air filled fraction of total soil volume
  PetscReal :: om_frac
  PetscReal :: r_min, r_max, r_psi
  PetscReal :: rho_w, grav
  PetscReal :: ratio_diffusivity_water_gas
!  PetscReal :: h2osoi_vol

  PetscScalar, pointer :: sucsat(:)
  PetscScalar, pointer :: watfc(:)
  PetscScalar, pointer :: cellorg(:)
  PetscScalar, pointer :: bsw(:)
  PetscScalar, pointer :: soilpsi(:)

  tk = tc + 273.15d0
! anaerobic fraction
!  saturation = global_auxvar%sat(iphase) 
!  h2osoi_vol = porosity * saturation  !? 
! what is the difference among h2osoi_vol (for dust model?), watsat, and volumetric water content?
! here I assume that they are identical     
  theta = porosity * saturation

  surface_tension_water = 73.d-3   ! (J/m^2), Arah and Vinten 1995

  ! Set parameters from simple-structure model to calculate anoxic fratction (Arah and Vinten 1995)
  rij_kro_a = 1.5d-10              !  Arah and Vinten 1995
  rij_kro_alpha = 1.26             !  Arah and Vinten 1995
  rij_kro_beta = 0.6               !  Arah and Vinten 1995
  rij_kro_gamma = 0.6              !  Arah and Vinten 1995
  rij_kro_delta = 0.85             !  Arah and Vinten 1995

  organic_max  = 130.0             ! organic matter content (kg/m3) where soil is assumed to act like peat
  ! for diffusion. Very large values will lead to all soil being treated as mineral. Negative values will lead   ! to all soil being treated as peat.


!  pH = 6.5  !!! set all soils with the same pH as placeholder here
  co2diff_con1 =   0.1325
  co2diff_con2 =   0.0009

!  data (d_con_w(1,i),i=1,3) /0.9798_r8, 0.02986_r8, 0.0004381_r8/ ! CH4
!  data (d_con_w(2,i),i=1,3) /1.172_r8, 0.03443_r8, 0.0005048_r8/ ! O2
!  data (d_con_w(3,i),i=1,3) /0.939_r8, 0.02671_r8, 0.0004095_r8/ ! CO2

  
  d_con_g_1_o2 = 0.1759
  d_con_g_2_o2 = 0.00117

  d_con_w_1_o2 = 1.172
  d_con_w_2_o2 = 0.03443
  d_con_w_3_o2 = 0.0005048

  rho_w  = 1.d3                   ! (kg/m3)
  grav   = 9.80616                ! acceleration of gravity ~ m/s^2
 

!  call VecGetArrayReadF90(clm_pf_idata%watsat_pf, watsat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%sucsat_pf, sucsat, ierr)
  call VecGetArrayReadF90(clm_pf_idata%watfc_pf, watfc, ierr)
  call VecGetArrayReadF90(clm_pf_idata%cellorg_pf, cellorg, ierr)
  call VecGetArrayReadF90(clm_pf_idata%bsw_pf, bsw, ierr)
  call VecGetArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi, ierr)
!  call VecGetArrayReadF90(clm_pf_idata%o2_decomp_depth_unsat_pf, o2_decomp_depth_unsat, ierr)
!  call VecGetArrayReadF90(clm_pf_idata%o2_decomp_depth_sat_pf, o2_decomp_depth_sat, ierr)
!  call VecGetArrayReadF90(clm_pf_idata%conc_o2_unsat_pf, conc_o2_unsat, ierr)
!  call VecGetArrayReadF90(clm_pf_idata%conc_o2_sat_pf, conc_o2_sat, ierr)

!  f_a = 1.0 - watfc(local_id) / !watsat(local_id)
!  e_a = watsat(local_id) - watfc(local_id)
  f_a = 1.0 - watfc(local_id) / porosity
  e_a = porosity - watfc(local_id)

!  if (clm_pf_idata%use_lch4) then
  if (organic_max > 0.0) then
     om_frac = min(cellorg(local_id)/organic_max, 1.0)
  else
     om_frac = 1.0
  end if

  diffus = (d_con_g_1_o2 + d_con_g_2_o2*tk) * 1.d-4 * &
        (om_frac * f_a**(10.0/3.0) / porosity**2 + &
        (1.0 - om_frac) * e_a**2 * f_a**(3.0 / bsw(local_id)))
!        (om_frac * f_a**(10.0/3.0) / watsat(local_id)**2 + &

  ! calculate anoxic fraction of soils
  ! use rijtema and kroess model after Riley et al., 2000
  ! caclulated r_psi as a function of psi

  if(saturation < WT_saturation) then
     tmp_real = soilpsi(local_id)
     r_min = 2.0 * surface_tension_water / (rho_w * grav * abs(soilpsi(local_id)))
     r_max = 2.0 * surface_tension_water / (rho_w * grav * 0.1)
     r_psi = sqrt(r_min * r_max)

     ratio_diffusivity_water_gas = (d_con_g_1_o2 + d_con_g_2_o2*tk) * 1.d-4 / &
          ((d_con_w_1_o2 + d_con_w_2_o2*tk + d_con_w_3_o2*tk**2) * 1.d-9)

!     if (o2_decomp_depth_unsat(local_id) .ne. spval .and. &
!        conc_o2_unsat(local_id) .ne. spval .and. &
!        o2_decomp_depth_unsat(local_id) > 0.0) then
     if(rate_o2 > 0.0) then
        anaerobic_frac = exp(-rij_kro_a * r_psi**(-rij_kro_alpha) * &
            rate_o2**(-rij_kro_beta) * &
            conc_o2**rij_kro_gamma * (theta + &
            ratio_diffusivity_water_gas * &
            porosity)**rij_kro_delta)
!            o2_decomp_depth_unsat(local_id)**(-rij_kro_beta) * &
!            watsat(local_id))**rij_kro_delta)
     else
        anaerobic_frac = 0.0
     endif
  else
     anaerobic_frac = 1.0d0   ! saturated zone, all anaerobic
! anoxia_wtsat = .false by default, NaN in o2_decomp_depth_sat
!     if (anoxia_wtsat) then ! Average saturated fraction values into anaerobic_frac(c,j).
!         r_min = 2.0 * surface_tension_water / (rho_w * grav * abs(grav * 1.e-6 * sucsat(local_id)))
!         r_max = 2.0 * surface_tension_water / (rho_w * grav * 0.1)
!         r_psi = sqrt(r_min * r_max)
!         ratio_diffusivity_water_gas = (d_con_g_1_o2 + d_con_g_2_o2*tk) * 1.d-4 / &
!             ((d_con_w_1_o2 + d_con_w_2_o2*tk + d_con_w_3_o2*tk**2) * 1.d-9)

 !        if (o2_decomp_depth_sat(local_id) .ne. spval .and. &
 !            conc_o2_sat(local_id) .ne. spval .and. &
 !            o2_decomp_depth_sat(local_id) > 0.0) then
 !            anaerobic_frac = exp(-rij_kro_a * r_psi**(-rij_kro_alpha) * &
 !                      o2_decomp_depth_sat(local_id)**(-rij_kro_beta) * &
 !                      conc_o2_sat(local_id)**rij_kro_gamma * (watsat(local_id) +  &
 !                      ratio_diffusivity_water_gas * watsat(local_id))**rij_kro_delta)
!             anaerobic_frac = 0.0
!         else
!             anaerobic_frac = 0.0
!         endif
!               anaerobic_frac(c,j) = (1._r8 - finundated(c))*anaerobic_frac(c,j) + finundated(c)*anaerobic_frac_sat
!     end if
  endif

!  call VecRestoreArrayReadF90(clm_pf_idata%watsat_pf, watsat, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%watfc_pf, watfc, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%cellorg_pf, cellorg, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%bsw_pf, bsw, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi, ierr)
!  call VecRestoreArrayReadF90(clm_pf_idata%o2_decomp_depth_unsat_pf, o2_decomp_depth_unsat, ierr)
!  call VecRestoreArrayReadF90(clm_pf_idata%o2_decomp_depth_sat_pf, o2_decomp_depth_sat, ierr)
!  call VecRestoreArrayReadF90(clm_pf_idata%conc_o2_unsat_pf, conc_o2_unsat, ierr)
!  call VecRestoreArrayReadF90(clm_pf_idata%conc_o2_sat_pf, conc_o2_sat, ierr)
#else
  if(saturation < WT_saturation) then
     anaerobic_frac = 0.0d0
  else
     anaerobic_frac = 1.0d0
  endif
  diffus = 1.0d-6  ! not sure about this 
#endif

!  GetAnaerobicFraction = anaerobic_frac

end subroutine GetAnaerobicFraction
 
Function GetpHResponse(pH)

  PetscReal, parameter :: rpi = 3.14159265358979323846
  PetscReal :: f_ph, pH, GetpHResponse

! ph function from Parton et al., (2001, 1996)
!  k_nitr_ph_vr(c,j) = 0.56 + atan(rpi * 0.45 * (-5.+ pH(c)))/rpi
#ifdef CLM_PFLOTRAN
  f_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + pH))/rpi
#else
  f_ph = 1.0
#endif

  GetpHResponse = f_ph

end function GetpHResponse

end module Reaction_Sandbox_CLM45_class
