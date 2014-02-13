module Reaction_Sandbox_CLMCND_class

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module

! -----------------------------------------------------------------------------
! description
! Extend reaction_sandbox_clm_cn.F90 to include nitrification & denitrification,
! and plant uptake and immobilization of ammonia and nitrate
! CLM-CN follows Thornton and Rosenbloom 2005
! de/nitrification follows Parton et al. 1996 and Dickinson et al. 2002 
! competition is implemented by regulating the rate by N/(S + N)
! preference of ammonia over nitrate by plants and microbes is implemented by
! regulating the rate by I/(I + ammonia) for uptake of nitrate 
! Suggested by Fengming Yuan, immplemented by Guoping Tang 02/04/2014
! -----------------------------------------------------------------------------

  implicit none
  
  private
  
#include "finclude/petscsys.h"

                          ! 14.00674d0 / 12.011d0
  PetscReal, parameter :: CN_ratio_mass_to_mol = 1.16616d0 
  PetscInt, parameter :: CARBON_INDEX = 1
  PetscInt, parameter :: NITROGEN_INDEX = 2
  PetscInt, parameter :: SOM_INDEX = 1
  PetscInt, parameter :: iphase = 1
  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_CLM4 = 1 
  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_Q10 = 2 
  PetscReal, parameter :: rpi = 3.14159265358979323846
  PetscReal, parameter :: N_molecular_weight = 14.0067d0

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_clmcnd_type

    PetscInt :: temperature_response_function
    PetscReal :: Q10
    PetscReal :: half_saturation_nh3
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh3_no3
    PetscReal :: n2o_frac_mineralization ! fraction of n2o from net N mineralization
    PetscReal :: k_nitr_max    ! nitrification rate
    PetscReal :: k_nitr_n2o    ! N2O production rate with nitrification
    PetscReal :: k_deni_max    ! denitrification rate
    PetscReal :: x0eps

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
    PetscInt :: species_id_n2
    PetscInt :: species_id_nplant

    type(pool_type), pointer :: pools
    type(clmcnd_reaction_type), pointer :: reactions
  contains
    procedure, public :: ReadInput => CLMCND_Read
    procedure, public :: Setup => CLMCND_Setup
    procedure, public :: Evaluate => CLMCND_React
    procedure, public :: Destroy => CLMCND_Destroy
  end type reaction_sandbox_clmcnd_type
  
  type :: pool_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: CN_ratio
    type(pool_type), pointer :: next
  end type pool_type
  
  type :: clmcnd_reaction_type
    character(len=MAXWORDLENGTH) :: upstream_pool_name
    character(len=MAXWORDLENGTH) :: downstream_pool_name
    PetscReal :: rate_constant
    PetscReal :: respiration_fraction
    PetscReal :: inhibition_constant
    type(clmcnd_reaction_type), pointer :: next
  end type clmcnd_reaction_type
  
  public :: CLMCND_Create

contains

! ************************************************************************** !

function CLMCND_Create()
  ! 
  ! Allocates CLMCND reaction object.
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
  ! 

  implicit none
  
  type(reaction_sandbox_clmcnd_type), pointer :: CLMCND_Create
  
  allocate(CLMCND_Create)

  CLMCND_Create%temperature_response_function=TEMPERATURE_RESPONSE_FUNCTION_CLM4
  CLMCND_Create%Q10 = 1.5d0
  CLMCND_Create%half_saturation_nh3 = 1.0d-15
  CLMCND_Create%half_saturation_no3 = 1.0d-15
  CLMCND_Create%inhibition_nh3_no3 = 1.0d-15 
  CLMCND_Create%k_nitr_max = 1.0d-6   ! nitrification rate
  CLMCND_Create%k_nitr_n2o = 3.5d-8   
  CLMCND_Create%k_deni_max = 2.5d-5   ! denitrification rate
  CLMCND_Create%n2o_frac_mineralization = 0.02d0  ! Parton et al. 2001
  CLMCND_Create%x0eps = 1.0d-20

  CLMCND_Create%nrxn = 0
  CLMCND_Create%npool = 0
  nullify(CLMCND_Create%CN_ratio)
  nullify(CLMCND_Create%rate_constant)
  nullify(CLMCND_Create%respiration_fraction)
  nullify(CLMCND_Create%inhibition_constant)
  nullify(CLMCND_Create%pool_id_to_species_id)
  nullify(CLMCND_Create%upstream_pool_id)
  nullify(CLMCND_Create%downstream_pool_id)
  CLMCND_Create%species_id_co2 = 0
  CLMCND_Create%species_id_nh3 = 0
  CLMCND_Create%species_id_no3 = 0
  CLMCND_Create%species_id_n2o = 0
  CLMCND_Create%species_id_n2 = 0
  CLMCND_Create%species_id_nplant = 0
  nullify(CLMCND_Create%next)
  nullify(CLMCND_Create%pools)
  nullify(CLMCND_Create%reactions)

end function CLMCND_Create

! ************************************************************************** !

subroutine CLMCND_Read(this,input,option)
  ! 
  ! Reads input deck for reaction sandbox parameters
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_clmcnd_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  
  type(pool_type), pointer :: new_pool, prev_pool
  type(clmcnd_reaction_type), pointer :: new_reaction, prev_reaction
  
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
                       'CHEMISTRY,REACTION_SANDBOX,CLMCND')
    call StringToUpper(word)   

    select case(trim(word))
      case('TEMPERATURE_RESPONSE_FUNCTION')
        do
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit

         call InputReadWord(input,option,word,PETSC_TRUE)
         call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,REACTION_SANDBOX,CLMCNDP,TEMPERATURE RESPONSE FUNCTION')
         call StringToUpper(word)   

            select case(trim(word))
              case('CLM4')
                  this%temperature_response_function = &
                       TEMPERATURE_RESPONSE_FUNCTION_CLM4    
              case('Q10') 
                  this%temperature_response_function = &
                       TEMPERATURE_RESPONSE_FUNCTION_Q10    
                  call InputReadDouble(input,option,this%Q10)  
                  call InputErrorMsg(input,option,'Q10', 'CHEMISTRY,' // &
                       'REACTION_SANDBOX_CLMCND,TEMPERATURE RESPONSE FUNCTION')
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMCND,' // &
                                'TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
         enddo 

     case('X0EPS')
         call InputReadDouble(input,option,this%x0eps)
         call InputErrorMsg(input,option,'x0eps', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')

     case('AMMONIA_HALF_SATURATION')
         call InputReadDouble(input,option,this%half_saturation_nh3)
         call InputErrorMsg(input,option,'ammonia half saturation', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')

     case('NITRATE_HALF_SATURATION')
         call InputReadDouble(input,option,this%half_saturation_no3)
         call InputErrorMsg(input,option,'nitrate half saturation', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')


     case('AMMONIA_INHIBITION_COEF')
         call InputReadDouble(input,option,this%inhibition_nh3_no3)
         call InputErrorMsg(input,option,'ammonia inhibition coefficient', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')

     case('NITRIFICATION_RATE_COEF')
         call InputReadDouble(input,option,this%k_nitr_max)
         call InputErrorMsg(input,option,'nitrification rate coefficient', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')

     case('DENITRIFICATION_RATE_COEF')
         call InputReadDouble(input,option,this%k_deni_max)
         call InputErrorMsg(input,option,'denitrification rate coefficient', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')

     case('N2O_FRAC_MINERALIZATION')
         call InputReadDouble(input,option,this%n2o_frac_mineralization)
         call InputErrorMsg(input,option,'n2o fraction from mineralization', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')

     case('N2O_RATE_COEF_NITRIFICATION')
         call InputReadDouble(input,option,this%k_nitr_n2o)
         call InputErrorMsg(input,option,'N2O rate coefficient from nirification', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')

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
            'CHEMISTRY,REACTION_SANDBOX,CLMCND,POOLS')
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
                             'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('UPSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%upstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'upstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')
            case('DOWNSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%downstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'downstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')
            case('RATE_CONSTANT')
              call InputReadDouble(input,option,rate_constant)
              call InputErrorMsg(input,option,'rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLMCND RATE CONSTANT UNITS'
                call InputDefaultMsg(input,option)
              else              
                rate_constant = rate_constant * &
                  UnitsConvertToInternal(word,option)
              endif
            case('TURNOVER_TIME')
              call InputReadDouble(input,option,turnover_time)
              call InputErrorMsg(input,option,'turnover time', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLMCND TURNOVER TIME UNITS'
                call InputDefaultMsg(input,option)
              else              
                turnover_time = turnover_time * &
                  UnitsConvertToInternal(word,option)
              endif
            case('RESPIRATION_FRACTION')
              call InputReadDouble(input,option,new_reaction%respiration_fraction)
              call InputErrorMsg(input,option,'respiration fraction', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')
            case('N_INHIBITION')
              call InputReadDouble(input,option,new_reaction%inhibition_constant)
              call InputErrorMsg(input,option,'inhibition constant', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMCND,REACTION')
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMCND,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
              call printErrMsg(option)
          end select
        enddo
        
        ! check to ensure that one of turnover time or rate constant is set.
        if (turnover_time > 0.d0 .and. rate_constant > 0.d0) then
          option%io_buffer = 'Only TURNOVER_TIME or RATE_CONSTANT may ' // &
            'be included in a CLMCND reaction definition, but not both. ' // &
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
            'zero and one (i.e. 0. <= rf <= 1.) in a CLMCND reaction ' // &
            'definition. See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        endif
        ! If no downstream pool exists, ensure that respiration fraction = 1
        if (len_trim(new_reaction%downstream_pool_name) < 1 .and. &
            (1.d0 - new_reaction%respiration_fraction) > 1.d-40) then
          option%io_buffer = 'Respiratory fraction (rf) must be set to ' // &
            '1.0 if no downstream pool is specified in a CLMCND reaction ' // &
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
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMCND keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine CLMCND_Read

! ************************************************************************** !

subroutine CLMCND_Setup(this,reaction,option)
  ! 
  ! Sets up CLMCND reaction after it has been read from input
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
  ! 

  use Reaction_Aux_module
  use Option_module
  use String_module
  use Immobile_Aux_module
  
  implicit none

  class(reaction_sandbox_clmcnd_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  
  character(len=MAXWORDLENGTH), allocatable :: pool_names(:)
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt :: icount

  type(pool_type), pointer :: cur_pool
  type(clmcnd_reaction_type), pointer :: cur_rxn
  
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
        option%io_buffer = 'For CLMCND pools with no CN ratio defined, ' // &
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
  
  word = 'CO2(aq)'
  this%species_id_co2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'Ammonia'
  this%species_id_nh3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'Nitrate'
  this%species_id_no3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  
  word = 'N2O(aq)'
  this%species_id_n2o = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  
  word = 'N2(aq)'
  this%species_id_n2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

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
        option%io_buffer = 'For CLMCND reactions, downstream pools ' // &
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
  
end subroutine CLMCND_Setup

! ************************************************************************** !

subroutine CLMCND_React(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                       global_auxvar,porosity,volume,reaction,option, local_id)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
  ! 

  use Option_module
  use Reaction_Aux_module

  implicit none

  class(reaction_sandbox_clmcnd_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscInt :: local_id

  call CLMCND_PlantNUptake(this,Residual,Jacobian,compute_derivative,&
           rt_auxvar,global_auxvar,porosity,volume,reaction,option, local_id)

  call CLMCND_Decomposition(this,Residual,Jacobian,compute_derivative,&
           rt_auxvar,global_auxvar,porosity,volume,reaction,option, local_id)

  call CLMCND_Nitrification(this,Residual,Jacobian,compute_derivative,&
           rt_auxvar,global_auxvar,porosity,volume,reaction,option, local_id)

  call CLMCND_Denitrification(this,Residual,Jacobian,compute_derivative,&
           rt_auxvar,global_auxvar,porosity,volume,reaction,option, local_id)

end subroutine CLMCND_React

! ************************************************************************** !

subroutine CLMCND_PlantNUptake(this,Residual,Jacobian,compute_derivative,&
           rt_auxvar,global_auxvar,porosity,volume,reaction,option, local_id)

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

  class(reaction_sandbox_clmcnd_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscInt :: local_id
  PetscErrorCode :: ierr

  PetscReal :: c_nh3         ! concentration (mole/L)
  PetscReal :: f_nh3         ! nh3 / (half_saturation + nh3)
  PetscReal :: d_nh3         ! half_saturation / (half_saturation + nh3)^2
  PetscReal :: f_nh3_inhibit ! inhibition_coef/(inhibition_coef + nh3)
  PetscReal :: c_no3         ! concentration (mole/L)
  PetscReal :: f_no3         ! no3 / (half_saturation + no3)
  PetscReal :: d_no3         ! half_saturation/(no3 + half_saturation)^2 
  PetscReal :: temp_real

  PetscScalar, pointer :: rate_plantnuptake_pf_loc(:)   !
  PetscReal :: rate_nplant
  PetscReal :: rate_nplant_no3
  PetscReal :: drate_nplant       !drate/dnh4+ 
  PetscReal :: drate_nplant_no3   !drate/dno3-

  PetscInt :: ires_nh3, ires_no3, ires_nplant

  c_nh3     = rt_auxvar%total(this%species_id_nh3, iphase)
  temp_real = c_nh3 - this%x0eps + this%half_saturation_nh3
  f_nh3     = (c_nh3 - this%x0eps) / temp_real 
  d_nh3     = this%half_saturation_nh3 / temp_real / temp_real

  if (this%species_id_no3 > 0) then
      c_no3 = rt_auxvar%total(this%species_id_no3, iphase)
      temp_real = c_no3 -this%x0eps + this%half_saturation_no3
      f_no3 = (c_no3 - this%x0eps) / temp_real 
      d_no3 = this%half_saturation_no3 / temp_real / temp_real 

      temp_real = this%inhibition_nh3_no3 + c_nh3 
      f_nh3_inhibit = this%inhibition_nh3_no3/temp_real
  endif 

  ! indices for C and N species
  ires_nh3 = this%species_id_nh3
  ires_no3 = this%species_id_no3
  ires_nplant = this%species_id_nplant + reaction%offset_immobile 

#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%rate_plantnuptake_pf, &
       rate_plantnuptake_pf_loc, ierr)

  rate_nplant = rate_plantnuptake_pf_loc(local_id) * volume ! mol/m3/s * m3

  call VecRestoreArrayReadF90(clm_pf_idata%rate_plantnuptake_pf, &
       rate_plantnuptake_pf_loc, ierr)
#else
  rate_nplant = 0.0d0
#endif

  rate_nplant_no3 = rate_nplant * f_nh3_inhibit   ! NO3 uptake rate

  drate_nplant = rate_nplant * d_nh3
  rate_nplant = rate_nplant * f_nh3 

  Residual(ires_nh3) = Residual(ires_nh3) + rate_nplant
  Residual(ires_nplant) = Residual(ires_nplant) - rate_nplant

  if (compute_derivative) then 
     Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3)+drate_nplant * &
       rt_auxvar%aqueous%dtotal(this%species_id_nh3,this%species_id_nh3,iphase)

     Jacobian(ires_nplant,ires_nh3)=Jacobian(ires_nplant,ires_nh3)-drate_nplant
  endif

  if(this%species_id_no3 > 0) then
    drate_nplant_no3 = rate_nplant_no3 * d_no3
    rate_nplant_no3  = rate_nplant_no3 * f_no3

    Residual(ires_no3) = Residual(ires_no3) + rate_nplant_no3
    Residual(ires_nplant) = Residual(ires_nplant) - rate_nplant_no3

    if (compute_derivative) then 
     Jacobian(ires_no3,ires_no3)=Jacobian(ires_no3,ires_no3)+drate_nplant_no3* &
       rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_no3,iphase)

     Jacobian(ires_nplant,ires_no3)=Jacobian(ires_nplant,ires_no3)-drate_nplant_no3
    endif
  endif

end subroutine CLMCND_PlantNUptake

! ************************************************************************** !

subroutine CLMCND_Decomposition(this,Residual,Jacobian,compute_derivative,&
           rt_auxvar,global_auxvar,porosity,volume,reaction,option, local_id)

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

  class(reaction_sandbox_clmcnd_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscReal :: saturation
  PetscReal :: theta
  PetscInt :: local_id
  PetscErrorCode :: ierr

  PetscReal :: c_nh3         ! concentration (mole/L)
  PetscReal :: f_nh3         ! nh3 / (half_saturation + nh3)
  PetscReal :: d_nh3         ! half_saturation / (half_saturation + nh3)^2
  PetscReal :: f_nh3_inhibit ! inhibition_coef/(inhibition_coef + nh3)
  PetscReal :: d_nh3_inhibit ! -inhibition_coef/(inhibition_coef + nh3)^2
  PetscReal :: c_no3         ! concentration (mole/L)
  PetscReal :: f_no3         ! no3 / (half_saturation + no3)
  PetscReal :: d_no3         ! half_saturation/(no3 + half_saturation)^2 
  PetscReal :: temp_real

  PetscInt :: irxn
  PetscInt :: ipool_up, ipool_down
  PetscReal :: CN_ratio_up, CN_ratio_down
  PetscBool :: constant_CN_ratio_up
  PetscReal :: resp_frac

  PetscInt :: ispec_uc, ispec_un, ispec_d   ! species id for upstream C, N, and downstream
  PetscInt :: ires_uc, ires_un, ires_d      ! id used for residual and Jacobian
  PetscInt :: ires_co2, ires_nh3, ires_n2o, ires_no3
  PetscReal :: stoich_uc, stoich_un, stoich_dc, stoich_mineraln, stoich_co2

  PetscReal :: scaled_rate_const
  PetscReal :: rate       ! mole/s 
  PetscReal :: drate_uc   ! d Rate / d upstream c
  PetscReal :: drate_nh3  ! d Rate / d nh3 ammonia limitation

! for N immobilization reactions with NO3 as N source 
  PetscReal :: rate_no3      ! mole/s
  PetscReal :: drate_no3     ! d Rate_no3 / d no3 
  PetscReal :: drate_uc_no3  ! d Rate_no3 / d uc 
  PetscReal :: drate_nh3_no3 ! d Rate_no3 / d nh3 

  PetscInt :: i, j
  PetscReal :: tc     ! temperature in C
  PetscReal :: f_t    ! temperature response function
  PetscReal :: f_w    ! moisture response function

! save mineral N fraction and decomposition rate for net N mineralization and N2O calculation 
  PetscReal :: net_n_mineralization_rate
  PetscReal :: ph, f_ph
  PetscReal :: rate_n2o, drate_n2o

  PetscReal, parameter :: rpi = 3.14159265358979323846

  PetscReal :: c_uc, c_un

  c_nh3     = rt_auxvar%total(this%species_id_nh3, iphase)
  temp_real = c_nh3 - this%x0eps + this%half_saturation_nh3
  f_nh3     = (c_nh3 - this%x0eps) / temp_real 
  d_nh3     = this%half_saturation_nh3 / temp_real / temp_real

  if (this%species_id_no3 > 0) then
      c_no3 = rt_auxvar%total(this%species_id_no3, iphase)
      temp_real = c_no3 -this%x0eps + this%half_saturation_no3
      f_no3 = (c_no3 - this%x0eps)/ temp_real 
      d_no3 = this%half_saturation_no3 / temp_real / temp_real 

      temp_real = this%inhibition_nh3_no3 + c_nh3 
      f_nh3_inhibit = this%inhibition_nh3_no3/temp_real
      d_nh3_inhibit = -1.0d0*this%inhibition_nh3_no3/temp_real/temp_real
  endif 

  ires_co2 = this%species_id_co2
  ires_nh3 = this%species_id_nh3
  ires_no3 = this%species_id_no3
  ires_n2o = this%species_id_n2o

! temperature response function 
  tc = global_auxvar%temp(1) 

#ifdef CLM_PFLOTRAN 
  f_t = GetTemperatureResponse(tc, this%temperature_response_function, this%Q10) 
#else
  f_t = 1.0d0
#endif

  saturation = global_auxvar%sat(1)
  theta = saturation * porosity 

  ! moisture response function 
#ifdef CLM_PFLOTRAN
  f_w = GetMoistureResponse(theta, local_id)
#else
  f_w = 1.0d0
#endif

  if(f_t < 1.0d-20 .or. f_w < 1.0d-20) then
     return
  endif

  net_n_mineralization_rate = 0.0d0

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

!    if(c_uc < 1.0d-20) cycle

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
    rate     = scaled_rate_const * (c_uc - this%x0eps)
    drate_uc = scaled_rate_const 

!   NH3 limiting
    if(stoich_mineraln < 0.0d0) then
!       if(c_nh3 < 1.0d-20) then
!          rate = 0.0d0
!          drate_uc = 0.0d0
!       endif

       drate_nh3 = rate * d_nh3 
       rate      = rate * f_nh3
       drate_uc  = drate_uc * f_nh3
    endif 

!    write(*,*)(Residual(i), i = 1, reaction%ncomp)
!    write(*, *) rate, drate_nh3, c_nh3, stoich_mineraln
    ! calculation of residual
    ! carbon
    Residual(ires_co2) = Residual(ires_co2) - stoich_co2 * rate
    
    ! nitrogen
    Residual(ires_nh3) = Residual(ires_nh3) - stoich_mineraln * rate

    net_n_mineralization_rate=net_n_mineralization_rate+stoich_mineraln*rate

    ! C species in upstream pool (litter or SOM)
    ires_uc = reaction%offset_immobile + ispec_uc
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
    
!    write(*,*)(Residual(i), i = 1, reaction%ncomp)
!   start residual calculation for N immobilization reaction with NO3 uptake
!   if nitrate is available, N immobilization decomposition reactions occurs
!   with rate depending on NH3, with reduced rate if NH3 is abundent  
    if(this%species_id_no3 > 0 .and. stoich_mineraln < 0.d0) then
!       rate_no3     = c_uc * scaled_rate_const * f_nh3_inhibit
!       drate_uc_no3 = scaled_rate_const * f_nh3_inhibit

!      NO3 limiting
!       drate_no3     = rate_no3 * d_no3 
!       rate_no3      = rate_no3 * f_no3
!       drate_uc_no3  = drate_uc_no3 * f_no3
!       drate_nh3_no3 = rate_no3 * 

       rate_no3     = scaled_rate_const * (c_uc - this%x0eps) * f_no3 * f_nh3_inhibit
       drate_uc_no3 = scaled_rate_const * 1.0d0 * f_no3 * f_nh3_inhibit
       drate_no3    = scaled_rate_const * c_uc  * d_no3 * f_nh3_inhibit
       drate_nh3_no3= scaled_rate_const * c_uc  * f_no3 * d_nh3_inhibit

!       write(*, *) rate_no3, drate_no3, stoich_mineraln, stoich_un, c_nh3, c_no3, f_nh3_inhibit
    ! residuals
    ! carbon
       Residual(ires_co2) = Residual(ires_co2) - stoich_co2 * rate_no3
    
    ! nitrogen
       Residual(ires_no3) = Residual(ires_no3) - stoich_mineraln * rate_no3
       net_n_mineralization_rate=net_n_mineralization_rate+stoich_mineraln*rate_no3

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
    
!  write(*, *)'----------'
!    write(*,*)(Residual(i), i = 1, reaction%ncomp)
    endif 
!   end residual calculation for N immobilization reaction with NO3 uptake

!    return

    if (compute_derivative) then
!  do i = 1, reaction%ncomp
!     write(*,*)(Jacobian(i, j), j = 1, reaction%ncomp)
!  enddo
 
! with respect to upstream c    
      ! carbon
      Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
        stoich_co2 * drate_uc

      ! nitrogen
         ! dR_N/dC_u = d nR/dC_u = dn/dC_u R + n dR/dC_u
         ! first, n dR/dC_u
      Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - &
        stoich_mineraln * drate_uc

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

!with respect to upstream n
      if (.not.constant_CN_ratio_up) then
        ! upstream n, dR_Nu/dNu = d uR/dNu = R/Cu = dR/dCu
         Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) -(-1.d0)*drate_uc

        ! mineral N, dR_N/dNu = d nR/du = Rd (u - (1-f)d)/dNu = dR/dCu 
         Jacobian(ires_nh3,ires_un) = Jacobian(ires_nh3,ires_un) - drate_uc
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

      endif

!  do i = 1, reaction%ncomp
!     write(*,*)(Jacobian(i, j), j = 1, reaction%ncomp)
!  enddo
!      return 
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

!with respect to upstream n
        if (.not.constant_CN_ratio_up) then
        ! upstream n, dR_Nu/dNu = d uR/dNu = R/Cu = dR/dCu
           Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) -(-1.d0)*drate_uc_no3

        ! mineral N, dR_N/dNu = d nR/du = Rd (u - (1-f)d)/dNu = dR/dCu 
           Jacobian(ires_no3,ires_un) = Jacobian(ires_no3,ires_un) - drate_uc_no3
        endif

!with respect to no3
!        if(stoich_mineraln < 0.0d0) then
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

! with respect to nh3
        ! carbon
          Jacobian(ires_co2,ires_nh3) = Jacobian(ires_co2,ires_nh3) - &
            stoich_co2 * drate_nh3_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_nh3,iphase)
  
        ! nitrogen
          Jacobian(ires_no3,ires_nh3) = Jacobian(ires_no3,ires_nh3) - &
            stoich_mineraln * drate_nh3_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_nh3,iphase)
  
        ! upstream C pool
          Jacobian(ires_uc,ires_nh3) = Jacobian(ires_uc,ires_nh3) - &
            (-1.d0) * stoich_uc * drate_nh3_no3
  
        ! upstream N pool
          if (.not.constant_CN_ratio_up) then
            Jacobian(ires_un,ires_nh3) = Jacobian(ires_un,ires_nh3) - &
              (-1.d0) * stoich_un * drate_nh3_no3 
          endif
  
        ! downstream pool
          if (ispec_d > 0) then
            Jacobian(ires_d,ires_nh3) = Jacobian(ires_d,ires_nh3) - &
              stoich_dc * drate_nh3_no3
          endif

!  write(*, *)'----------'
!  do i = 1, reaction%ncomp
!     write(*,*)(Jacobian(i, j), j = 1, reaction%ncomp)
!  enddo
!        endif
      endif  ! jacobian calculation for N immobilization reaction with NO3 uptake
    endif    ! jacobian calculation
 
! save net N mineralization and N2O calculation 

  enddo      !reactions

!  return

  if(this%species_id_n2o > 0 .and. net_n_mineralization_rate > 0.0d0) then
  ! temperature response function (Parton et al. 1996)
    f_t = -0.06d0 + 0.13d0 * exp( 0.07d0 * tc )
    f_w = ((1.27d0 - saturation)/0.67d0)**(3.1777d0) * &
        ((saturation - 0.0012d0)/0.5988d0)**2.84d0  
    ph = 6.5d0
    f_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + ph))/rpi

    if(f_t > 0.0d0 .and. f_w > 0.0d0 .and. f_ph > 0.0d0) then 
      temp_real = f_t * f_w * f_ph

      if(temp_real > 1.0d0) then
         temp_real = 1.0d0
      endif

      temp_real = temp_real * this%n2o_frac_mineralization 
    
    ! residuals 
      rate_n2o = temp_real * net_n_mineralization_rate * f_nh3
 
      Residual(ires_nh3) = Residual(ires_nh3) + rate_n2o 

      Residual(ires_n2o) = Residual(ires_n2o) - 0.5d0 * rate_n2o

      if (compute_derivative) then
        drate_n2o = temp_real * net_n_mineralization_rate * d_nh3
        Jacobian(ires_nh3,ires_n2o) = Jacobian(ires_nh3,ires_n2o) + drate_n2o
        Jacobian(ires_n2o,ires_n2o) = Jacobian(ires_n2o,ires_n2o) - &
                                      0.5d0 * drate_n2o
      endif
    endif
  endif

end subroutine CLMCND_Decomposition

! ************************************************************************** !

subroutine CLMCND_Nitrification(this,Residual,Jacobian,compute_derivative,&
           rt_auxvar,global_auxvar,porosity,volume,reaction,option, local_id)

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

  class(reaction_sandbox_clmcnd_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscInt :: local_id
  PetscErrorCode :: ierr

  PetscReal :: temp_real

  PetscInt :: ires_nh3, ires_no3, ires_n2o

  PetscScalar, pointer :: bulkdensity(:)
  PetscReal :: rho_b
  PetscReal :: theta
  PetscReal :: c_nh3      ! mole/L
  PetscReal :: c_nh3_ugg  ! ug ammonia N / g soil
  PetscReal :: ph
  PetscReal :: rate_n2o, drate_n2o
  PetscReal :: rate_nitri, drate_nitri
  PetscReal :: f_t, f_w, f_ph
  PetscReal :: dfw_dnh3
  PetscReal :: saturation
  PetscReal :: tc

  ! indices for C and N species
  ires_nh3 = this%species_id_nh3
  ires_no3 = this%species_id_no3
  ires_n2o = this%species_id_n2o 

  saturation = global_auxvar%sat(1)
  theta = saturation * porosity

  tc = global_auxvar%temp(1) 

! nitrification (Dickinson et al. 2002)
  if(this%species_id_no3 > 0) then
    f_t = exp(0.08d0 * (tc - 298.0d0 + 273.15d0))
    f_w = saturation * (1.0d0 - saturation)

    rate_nitri = f_t * f_w * this%k_nitr_max * c_nh3
    Residual(ires_nh3) = Residual(ires_nh3) + rate_nitri
    Residual(ires_no3) = Residual(ires_no3) - rate_nitri
   
    if (compute_derivative) then
     drate_nitri = f_t*f_w*this%k_nitr_max*(0.25d0*c_nh3*c_nh3+2.0d0*c_nh3) &
                 / (0.25d0*c_nh3+1.0d0) / (0.25d0 * c_nh3 + 1.0d0) 

     Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) + drate_nitri * &
        rt_auxvar%aqueous%dtotal(this%species_id_nh3,this%species_id_nh3,iphase)

     Jacobian(ires_no3,ires_nh3) = Jacobian(ires_no3,ires_nh3) - drate_nitri * &
        rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_nh3,iphase)
    endif
  endif

! N2O production from nitrification (Parton et al. 1996)
#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%bulkdensity_dry_pf, bulkdensity, ierr)
  rho_b = bulkdensity(local_id) ! kg/m3
  call VecRestoreArrayReadF90(clm_pf_idata%bulkdensity_dry_pf, bulkdensity, ierr)
#else
  rho_b = 1.0d0
#endif
!             mole/L * 1000 L/m3 * g/mol / kg/m3 = g/kg = mg/g = 1000 ug/g  
!  c_nh3_ugg = c_nh3 / theta *1000.0d0 * N_molecular_weight / rho_b * 1000.0d0
  temp_real  = 1.0d0 / theta *1000.0d0 * N_molecular_weight / rho_b * 1000.0d0
  c_nh3     = rt_auxvar%total(this%species_id_nh3, iphase)
  c_nh3_ugg = c_nh3 * temp_real

  if(this%species_id_n2o > 0.0d0 .and. c_nh3_ugg > 3.0d0 ) then 
  ! temperature response function (Parton et al. 1996)
    f_t = -0.06d0 + 0.13d0 * exp( 0.07d0 * tc )

    f_w = ((1.27d0 - saturation)/0.67d0)**(3.1777d0) * &
        ((saturation - 0.0012d0)/0.5988d0)**2.84d0  

    ph = 6.5d0
    f_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + ph))/rpi

    if(f_t > 0.0d0 .and. f_w > 0.0d0 .and. f_ph > 0.0d0) then 
       if(f_w > 1.0d0) then
          f_w = 1.0d0
       endif

       if(f_ph > 1.0d0) then
          f_ph = 1.0d0
       endif

       rate_n2o = 1.0 - exp(-0.0105d0 * c_nh3_ugg)  ! need to change units 
       rate_n2o = rate_n2o * f_t * f_w * f_ph * this%k_nitr_n2o

       Residual(ires_nh3) = Residual(ires_nh3) + rate_n2o
       Residual(ires_n2o) = Residual(ires_n2o) - 0.5d0 * rate_n2o

       if (compute_derivative) then
           drate_n2o = 0.0105d0 * exp(-0.0105d0 * c_nh3_ugg) * temp_real
           drate_n2o = drate_n2o * f_t * f_w * f_ph * this%k_nitr_n2o

           Jacobian(ires_nh3,ires_nh3)=Jacobian(ires_nh3,ires_nh3)+drate_n2o * &
           rt_auxvar%aqueous%dtotal(this%species_id_nh3,this%species_id_nh3,iphase)

           Jacobian(ires_n2o,ires_nh3)=Jacobian(ires_n2o,ires_nh3)- &
           0.5d0 * drate_n2o * &
           rt_auxvar%aqueous%dtotal(this%species_id_n2o,this%species_id_nh3,iphase)
       endif
     endif
  endif

end subroutine CLMCND_Nitrification


! ************************************************************************** !

subroutine CLMCND_Denitrification(this,Residual,Jacobian,compute_derivative,&
           rt_auxvar,global_auxvar,porosity,volume,reaction,option, local_id)

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

  class(reaction_sandbox_clmcnd_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscInt :: local_id
  PetscErrorCode :: ierr

  PetscReal :: temp_real

  PetscInt :: ires_no3, ires_n2o, ires_n2

  PetscScalar, pointer :: bsw(:)
  PetscScalar, pointer :: bulkdensity(:)

  PetscReal :: s_min 
  PetscReal :: tc
  PetscReal :: f_t, f_w

  PetscReal :: c_no3      ! mole/L
  PetscReal :: rate_deni, drate_deni
  PetscReal :: saturation

  ! indices for C and N species
  ires_no3 = this%species_id_no3
  ires_n2o = this%species_id_n2o 
  ires_n2 = this%species_id_n2

! denitrification (Dickinson et al. 2002)
  if(this%species_id_n2 < 0) return 

#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%bsw_pf, bsw, ierr)
  temp_real = bsw(local_id)
  call VecRestoreArrayReadF90(clm_pf_idata%bsw_pf, bsw, ierr)
#else
  temp_real = 1.0d0
#endif

  tc = global_auxvar%temp(1) 
  f_t = exp(0.08d0 * (tc - 298.0d0 + 273.15d0))

  saturation = global_auxvar%sat(1)
  s_min = 0.6d0
  if(saturation > s_min) then
     f_w = (saturation - s_min)/(1.0d0 - s_min)
     f_w = f_w ** temp_real
  else
     return
  endif

  c_no3     = rt_auxvar%total(this%species_id_no3, iphase)

  if(f_w > 0.0d0) then
     rate_deni = this%k_deni_max * c_no3 * f_t * f_w 

     Residual(ires_no3) = Residual(ires_no3) + rate_deni
     Residual(ires_n2) = Residual(ires_n2) - 0.5d0*rate_deni
     
    if (compute_derivative) then
     
       drate_deni = f_t*f_w*this%k_deni_max

       Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) + drate_deni * &
        rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_no3,iphase)

       Jacobian(ires_n2,ires_no3)=Jacobian(ires_n2,ires_no3)-0.5d0*drate_deni * &
        rt_auxvar%aqueous%dtotal(this%species_id_n2,this%species_id_no3,iphase)

    endif
  endif

end subroutine CLMCND_Denitrification

! ************************************************************************** !

subroutine CLMCND_Destroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_clmcnd_type) :: this
  
  type(pool_type), pointer :: cur_pool, prev_pool
  type(clmcnd_reaction_type), pointer :: cur_reaction, prev_reaction
  
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
  
end subroutine CLMCND_Destroy

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

end module Reaction_Sandbox_CLMCND_class
