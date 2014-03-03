module Reaction_Sandbox_CLM_Decomp_class

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  use CLM_BGC_module

! -----------------------------------------------------------------------------
! description
! -----------------------------------------------------------------------------

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  PetscInt, parameter :: LITTER_DECOMP_CLMCN = 1 
  PetscInt, parameter :: LITTER_DECOMP_CLMMICROBE = 2 
  PetscReal, parameter :: CN_ratio_microbe = 9.32928d0   ! 8.0d0 
  PetscReal, parameter :: CN_ratio_bacteria = 5.8038d0   ! 5.0d0 
  PetscReal, parameter :: CN_ratio_fungi = 17.4924d0     !15.0d0 ! or 10.0 
  PetscReal, parameter :: CUE_max = 0.6d0
  PetscReal, parameter :: fraction_bacteria = 0.340927d0 ! 5.0**0.6/(5.0**0.6+15.0**0.6) 


  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_clm_decomp_type

    PetscInt :: temperature_response_function
    PetscReal :: Q10
    PetscInt :: litter_decomp_type          ! CLM-CN or CLM-Microbe
    PetscReal :: half_saturation_nh3
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh3_no3
    PetscReal :: n2o_frac_mineralization    ! fraction of n2o from net N mineralization
    PetscReal :: k_nitr_max                 ! nitrification rate
    PetscReal :: k_nitr_n2o                 ! N2O production rate with nitrification
    PetscReal :: k_deni_max                 ! denitrification rate
    PetscReal :: x0eps

    PetscInt :: npool                       ! litter or variable CN ration pools
    PetscReal, pointer :: pool_nc_ratio(:)         ! NC ratio in mole  

    PetscInt :: nrxn
    PetscReal, pointer :: rate_constant(:)         !nrxn

    PetscBool, pointer :: is_litter_decomp(:)
    PetscInt,  pointer :: upstream_c_id(:)         !nrxn
    PetscInt,  pointer :: upstream_n_id(:)         !nrxn
    PetscReal, pointer :: upstream_nc(:)           !nrxn
    PetscBool, pointer :: upstream_is_aqueous(:)   !nrxn

    PetscInt,  pointer :: n_downstream_pools(:)   !nrxn by maximum # of downstream pools
    PetscInt,  pointer :: downstream_id(:,:) !nrxn by maximum # of downstream pools
    PetscBool, pointer :: downstream_is_aqueous(:,:) !nrxn by maximum # of downstream pools
    PetscReal, pointer :: downstream_stoich(:,:)  !nrxn by maximum # of downstream pools
    PetscReal, pointer :: downstream_nc(:,:)  !nrxn by maximum # of downstream pools
    PetscReal, pointer :: mineral_c_stoich(:)     !nrxn by maximum # of downstream pools
    PetscReal, pointer :: mineral_n_stoich(:)     !nrxn by maximum # of downstream pools

    PetscInt :: species_id_co2
    PetscInt :: species_id_nh3
    PetscInt :: species_id_no3
    PetscInt :: species_id_n2o
    PetscInt :: species_id_dom
    PetscInt :: species_id_bacteria
    PetscInt :: species_id_fungi

    type(pool_type), pointer :: pools
    type(clm_decomp_reaction_type), pointer :: reactions
  contains
    procedure, public :: ReadInput => CLM_Decomp_Read
    procedure, public :: Setup => CLM_Decomp_Setup
    procedure, public :: Evaluate => CLM_Decomp_React
    procedure, public :: Destroy => CLM_Decomp_Destroy
  end type reaction_sandbox_clm_decomp_type
  
  type :: pool_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: stoich
    PetscReal :: nc_ratio
    type(pool_type), pointer :: next
  end type pool_type
  
  type :: clm_decomp_reaction_type
    character(len=MAXWORDLENGTH) :: upstream_pool_name
    type(pool_type), pointer :: downstream_pools
    PetscReal :: rate_constant
    type(clm_decomp_reaction_type), pointer :: next
  end type clm_decomp_reaction_type
  
  public :: CLM_Decomp_Create

contains

! ************************************************************************** !

function CLM_Decomp_Create()
  ! 
  ! Allocates CLM_Decomp reaction object.
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
  ! 

  implicit none
  
  type(reaction_sandbox_clm_decomp_type), pointer :: CLM_Decomp_Create
  
  allocate(CLM_Decomp_Create)

#ifdef CLM_PFLOTRAN
  CLM_Decomp_Create%temperature_response_function=TEMPERATURE_RESPONSE_FUNCTION_CLM4
#endif

  CLM_Decomp_Create%Q10 = 1.5d0
  CLM_Decomp_Create%litter_decomp_type=LITTER_DECOMP_CLMCN
  CLM_Decomp_Create%half_saturation_nh3 = 1.0d-15
  CLM_Decomp_Create%half_saturation_no3 = 1.0d-15
  CLM_Decomp_Create%inhibition_nh3_no3 = 1.0d-15 
  CLM_Decomp_Create%k_nitr_max = 1.0d-6   ! nitrification rate
  CLM_Decomp_Create%k_nitr_n2o = 3.5d-8   
  CLM_Decomp_Create%k_deni_max = 2.5d-5   ! denitrification rate
  CLM_Decomp_Create%n2o_frac_mineralization = 0.02d0  ! Parton et al. 2001
  CLM_Decomp_Create%x0eps = 1.0d-20

  CLM_Decomp_Create%npool = 0
  nullify(CLM_Decomp_Create%pool_nc_ratio)

  CLM_Decomp_Create%nrxn = 0
  nullify(CLM_Decomp_Create%rate_constant)
  nullify(CLM_Decomp_Create%is_litter_decomp)
  nullify(CLM_Decomp_Create%upstream_c_id)
  nullify(CLM_Decomp_Create%upstream_n_id)
  nullify(CLM_Decomp_Create%upstream_nc)
  nullify(CLM_Decomp_Create%upstream_is_aqueous)
  
  nullify(CLM_Decomp_Create%n_downstream_pools)
  nullify(CLM_Decomp_Create%downstream_id)
  nullify(CLM_Decomp_Create%downstream_is_aqueous)
  nullify(CLM_Decomp_Create%downstream_stoich)
  nullify(CLM_Decomp_Create%mineral_c_stoich)
  nullify(CLM_Decomp_Create%mineral_n_stoich)

  CLM_Decomp_Create%species_id_co2 = 0
  CLM_Decomp_Create%species_id_nh3 = 0
  CLM_Decomp_Create%species_id_no3 = 0
  CLM_Decomp_Create%species_id_n2o = 0
!  CLM_Decomp_Create%species_id_n2 = 0
  CLM_Decomp_Create%species_id_dom = 0
  CLM_Decomp_Create%species_id_bacteria = 0
  CLM_Decomp_Create%species_id_fungi = 0

  nullify(CLM_Decomp_Create%next)

  nullify(CLM_Decomp_Create%pools)
  nullify(CLM_Decomp_Create%reactions)

end function CLM_Decomp_Create

! ************************************************************************** !

subroutine CLM_Decomp_Read(this,input,option)
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
!  use CLM_BGC_module
 
  implicit none
  
  class(reaction_sandbox_clm_decomp_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  
  type(pool_type), pointer :: new_pool, prev_pool
  type(pool_type), pointer :: new_pool_rxn, prev_pool_rxn
  type(clm_decomp_reaction_type), pointer :: new_reaction, prev_reaction
  
  PetscReal :: rate_constant, turnover_time
  PetscReal :: temp_real
  
  nullify(new_pool)
  nullify(prev_pool)

  nullify(new_pool_rxn)
  nullify(prev_pool_rxn)

  nullify(new_reaction)
  nullify(prev_reaction)
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp')
    call StringToUpper(word)   
    select case(trim(word))
#ifdef CLM_PFLOTRAN
      case('TEMPERATURE_RESPONSE_FUNCTION')
        do
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit

         call InputReadWord(input,option,word,PETSC_TRUE)
         call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,TEMPERATURE RESPONSE FUNCTION')
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
                       'REACTION_SANDBOX_CLM_Decomp,TEMPERATURE RESPONSE FUNCTION')
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,' // &
                                'TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
         enddo 
#endif

     case('CLM-MICROBE-LITTER-DECOMPOSITION')
          this%litter_decomp_type = LITTER_DECOMP_CLMMICROBE    
!     case('LITTER_DECOMPOSITION')
!        do
!         call InputReadPflotranString(input,option)
!         if (InputError(input)) exit
!         if (InputCheckExit(input,option)) exit
!
!         call InputReadWord(input,option,word,PETSC_TRUE)
!         call InputErrorMsg(input,option,'keyword', &
!            'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,LITTER DECOMPOSITION TYPE')
!         call StringToUpper(word)   
!
!            select case(trim(word))
!              case('CLM-CN')
!                  this%litter_decomp_type = LITTER_DECOMP_CLMCN    
!1              case('CLM-MICROBE') 
!                  this%litter_decomp_type = LITTER_DECOMP_CLMMICROBE    
!              case default
!                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,' // &
!                                'LITTER DECOMPOSITION TYPE keyword: ' // &
!                                     trim(word) // ' not recognized.'
!                  call printErrMsg(option)
!            end select
!         enddo 

     case('X0EPS')
         call InputReadDouble(input,option,this%x0eps)
         call InputErrorMsg(input,option,'x0eps', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')

     case('AMMONIA_HALF_SATURATION')
         call InputReadDouble(input,option,this%half_saturation_nh3)
         call InputErrorMsg(input,option,'ammonia half saturation', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')

     case('NITRATE_HALF_SATURATION')
         call InputReadDouble(input,option,this%half_saturation_no3)
         call InputErrorMsg(input,option,'nitrate half saturation', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')


     case('AMMONIA_INHIBITION_COEF')
         call InputReadDouble(input,option,this%inhibition_nh3_no3)
         call InputErrorMsg(input,option,'ammonia inhibition coefficient', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')

     case('NITRIFICATION_RATE_COEF')
         call InputReadDouble(input,option,this%k_nitr_max)
         call InputErrorMsg(input,option,'nitrification rate coefficient', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')

     case('N2O_FRAC_MINERALIZATION')
         call InputReadDouble(input,option,this%n2o_frac_mineralization)
         call InputErrorMsg(input,option,'n2o fraction from mineralization', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')

     case('POOLS')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit   

          allocate(new_pool)
          new_pool%name = ''
          new_pool%nc_ratio = -999.d0
          nullify(new_pool%next)

          call InputReadWord(input,option,new_pool%name,PETSC_TRUE)
          call InputErrorMsg(input,option,'pool name', &
            'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,POOLS')
          call InputReadDouble(input,option,temp_real)
          if (InputError(input)) then
            new_pool%nc_ratio = -999.d0
          else
            ! convert CN ratio from mass C/mass N to mol C/mol N
             if(temp_real > 0.0d0 ) then
                new_pool%nc_ratio = 1.0d0/temp_real/CN_ratio_mass_to_mol
             endif
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
        new_reaction%rate_constant = -999.d0
        nullify(new_reaction%downstream_pools)
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
                             'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('UPSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%upstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'upstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')
            case('DOWNSTREAM_POOL')
              allocate(new_pool_rxn)
              new_pool_rxn%name = ''
              new_pool_rxn%stoich = 0.d0
              nullify(new_pool_rxn%next)

              call InputReadWord(input,option, &
                                 new_pool_rxn%name,PETSC_TRUE)
              call InputErrorMsg(input,option,'downstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')
              call InputReadDouble(input,option,new_pool_rxn%stoich)
              call InputErrorMsg(input,option,'Downstream pool stoich', 'CHEMISTRY,' // &
                  'REACTION_SANDBOX_CLM_Decomp,TEMPERATURE RESPONSE FUNCTION')

              if (associated(new_reaction%downstream_pools)) then
                  prev_pool_rxn%next => new_pool_rxn
              else
                  new_reaction%downstream_pools => new_pool_rxn
              endif
              prev_pool_rxn => new_pool_rxn
              nullify(new_pool_rxn)

            case('RATE_CONSTANT')
              call InputReadDouble(input,option,rate_constant)
              call InputErrorMsg(input,option,'rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLM_Decomp RATE CONSTANT UNITS'
                call InputDefaultMsg(input,option)
              else              
                rate_constant = rate_constant * &
                  UnitsConvertToInternal(word,option)
              endif
            case('TURNOVER_TIME')
              call InputReadDouble(input,option,turnover_time)
              call InputErrorMsg(input,option,'turnover time', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLM_Decomp TURNOVER TIME UNITS'
                call InputDefaultMsg(input,option)
              else              
                turnover_time = turnover_time * &
                  UnitsConvertToInternal(word,option)
              endif
!            case('RESPIRATION_FRACTION')
!              call InputReadDouble(input,option,new_reaction%resp_frac)
!              call InputErrorMsg(input,option,'respiration fraction', &
!                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')
!            case('N_INHIBITION')
!              call InputReadDouble(input,option,new_reaction%inhibition_constant)
!              call InputErrorMsg(input,option,'inhibition constant', &
!                     'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,REACTION')
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
              call printErrMsg(option)
          end select
        enddo
        
        ! check to ensure that one of turnover time or rate constant is set.
        if (turnover_time > 0.d0 .and. rate_constant > 0.d0) then
          option%io_buffer = 'Only TURNOVER_TIME or RATE_CONSTANT may ' // &
            'be included in a CLM_Decomp reaction definition, but not both. ' // &
            'See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        else if (turnover_time > 0.d0) then
          new_reaction%rate_constant = 1.d0 / turnover_time
        else
          new_reaction%rate_constant = rate_constant
        endif
        ! ensure that respiration fraction is 0-1.
 !       if (new_reaction%resp_frac < 0.d0 .or. &
 !                 new_reaction%resp_frac > 1.d0) then
 !         option%io_buffer = 'Respiratory fraction (rf) must be between ' // &
 !           'zero and one (i.e. 0. <= rf <= 1.) in a CLM_Decomp reaction ' // &
 !           'definition. See reaction with upstream pool "' // &
 !           trim(new_reaction%upstream_pool_name) // '".'
 !         call printErrMsg(option)
 !       endif
        ! If no downstream pool exists, ensure that respiration fraction = 1
!        if (len_trim(new_reaction%downstream_pool_name) < 1 .and. &
!            (1.d0 - new_reaction%resp_frac) > 1.d-40) then
!          option%io_buffer = 'Respiratory fraction (rf) must be set to ' // &
!            '1.0 if no downstream pool is specified in a CLM_Decomp reaction ' // &
!            'definition. See reaction with upstream pool "' // &
!            trim(new_reaction%upstream_pool_name) // '".'
!          call printErrMsg(option)
!        endif
        if (associated(this%reactions)) then
          prev_reaction%next => new_reaction
        else
          this%reactions => new_reaction
        endif
        prev_reaction => new_reaction
        nullify(new_reaction)        
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_Decomp keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine CLM_Decomp_Read

! ************************************************************************** !

subroutine CLM_Decomp_Setup(this,reaction,option)
  ! 
  ! Sets up CLM_Decomp reaction after it has been read from input
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
  ! 

  use Reaction_Aux_module
  use Option_module
  use String_module
  use Immobile_Aux_module
  use Utility_module, only : DeallocateArray
  
  implicit none

  class(reaction_sandbox_clm_decomp_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  
  character(len=MAXWORDLENGTH), allocatable :: pool_names(:)
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt, pointer :: species_id_pool_c(:)
  PetscInt, pointer :: species_id_pool_n(:)
  PetscBool, pointer :: pool_is_aqueous(:)

  PetscInt :: icount, jcount, max_downstream_pools, ipool
  PetscReal :: stoich_c, stoich_n

  type(pool_type), pointer :: cur_pool
  type(clm_decomp_reaction_type), pointer :: cur_rxn
  
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
 
  allocate(this%n_downstream_pools(this%nrxn))
 
  ! count # downstream pools in each reaction
  max_downstream_pools = -1
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1

    jcount = 0
    cur_pool => cur_rxn%downstream_pools

    do
      if (.not.associated(cur_pool)) exit
      jcount = jcount + 1
      cur_pool => cur_pool%next
    enddo

    this%n_downstream_pools(icount) = jcount

    if(max_downstream_pools < jcount) then
      max_downstream_pools = jcount
    endif 

    cur_rxn => cur_rxn%next
  enddo

  ! allocate and initialize arrays
  allocate(this%pool_nc_ratio(this%npool))

  allocate(this%rate_constant(this%nrxn))
  allocate(this%is_litter_decomp(this%nrxn))
  allocate(this%upstream_c_id(this%nrxn))
  allocate(this%upstream_n_id(this%nrxn))
  allocate(this%upstream_nc(this%nrxn))
  allocate(this%upstream_is_aqueous(this%nrxn))
  
  allocate(this%downstream_id(this%nrxn,max_downstream_pools))
  allocate(this%downstream_stoich(this%nrxn,max_downstream_pools))
  allocate(this%downstream_nc(this%nrxn,max_downstream_pools))
  allocate(this%downstream_is_aqueous(this%nrxn,max_downstream_pools))
  allocate(this%mineral_c_stoich(this%nrxn))
  allocate(this%mineral_n_stoich(this%nrxn))

  this%pool_nc_ratio = 0.d0
  this%rate_constant = 0.d0
  this%is_litter_decomp = PETSC_FALSE
  this%upstream_c_id = 0
  this%upstream_n_id = 0
  this%upstream_nc = -999.9
  this%upstream_is_aqueous = PETSC_FALSE

  this%downstream_id = 0
  this%downstream_is_aqueous = PETSC_FALSE
  this%downstream_stoich = 0.d0
  this%mineral_c_stoich = 0.d0
  this%mineral_n_stoich = 0.d0
  
! temporary array for mapping pools in reactions
  allocate(pool_names(this%npool))
  allocate(pool_is_aqueous(this%npool))
  allocate(species_id_pool_c(this%npool))
  allocate(species_id_pool_n(this%npool))

  pool_names = ''
  pool_is_aqueous = PETSC_FALSE
  species_id_pool_c = -999 
  species_id_pool_n = -999 

! pools
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    this%pool_nc_ratio(icount) = cur_pool%nc_ratio
    pool_names(icount) = cur_pool%name

    if (cur_pool%nc_ratio < 0.d0) then
      ! Since no CN ratio provided, must provide two species with the
      ! same name as the pool with C or N appended.
      word = trim(cur_pool%name) // 'C'
      species_id_pool_c(icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      word = trim(cur_pool%name) // 'N'
      species_id_pool_n(icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      if (species_id_pool_c(icount)<=0 .or. species_id_pool_n(icount)<=0) then
        option%io_buffer = 'For CLM_Decomp pools with no CN ratio defined, ' // &
          'the user must define two immobile species with the same root ' // &
          'name as the pool with "C" or "N" appended, respectively.'
        call printErrMsg(option)
      endif
    else ! only one species (e.g. SOMX)
      species_id_pool_c(icount) = &
        GetImmobileSpeciesIDFromName(cur_pool%name,reaction%immobile, &
                                     PETSC_FALSE,option)
      if(species_id_pool_c(icount) <= 0) then
          species_id_pool_c(icount) = GetPrimarySpeciesIDFromName(cur_pool%name, &
                reaction, PETSC_FALSE,option)
          if(species_id_pool_c(icount) <= 0) then
             option%io_buffer = 'CLM_Decomp pool: ' // cur_pool%name // &
             'is not specified either in the IMMOBILE_SPECIES or PRIMARY_SPECIES!'
             call printErrMsg(option)
          else
             pool_is_aqueous(icount) = PETSC_TRUE
          endif
      endif
    endif
    cur_pool => cur_pool%next
  enddo
 
 
! reactions
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
!   upstream pools
    icount = icount + 1
    ipool = StringFindEntryInList(cur_rxn%upstream_pool_name,pool_names)
    if (ipool == 0) then
      option%io_buffer = 'Upstream pool ' // &
        trim(cur_rxn%upstream_pool_name) // &
        'in reaction not found in list of pools.'
      call printErrMsg(option)
    else
      this%upstream_c_id(icount) = species_id_pool_c(ipool)
      this%upstream_n_id(icount) = species_id_pool_n(ipool)
      this%upstream_nc(icount) = this%pool_nc_ratio(ipool) 
      this%upstream_is_aqueous(icount) = pool_is_aqueous(ipool) 
      if(this%upstream_n_id(icount) > 0) then
         this%is_litter_decomp(icount) = PETSC_TRUE
      else
         if(this%upstream_nc(icount) < 0.0d0) then
            option%io_buffer = 'SOM decomposition reaction with upstream pool ' // &
              trim(cur_rxn%upstream_pool_name) // &
              'has negative C:N ratio in upstream pool.'
            call printErrMsg(option)
         endif
      endif
    endif

!   downstream pools
    jcount = 0
    cur_pool => cur_rxn%downstream_pools

    do
      if (.not.associated(cur_pool)) exit
      jcount = jcount + 1

      if (len_trim(cur_pool%name) > 0) then
        ipool = StringFindEntryInList(cur_pool%name,pool_names)
        if (ipool == 0) then
          option%io_buffer = 'Downstream pool "' // &
          trim(cur_pool%name) // &
          '" in reaction with upstream pool "' // &
          trim(cur_rxn%upstream_pool_name) // '" not found in list of pools.'
          call printErrMsg(option)
        else
          this%downstream_id(icount, jcount) = species_id_pool_c(ipool)
          this%downstream_stoich(icount, jcount) = cur_pool%stoich 
          this%downstream_nc(icount, jcount) = this%pool_nc_ratio(ipool) 
          this%downstream_is_aqueous(icount, jcount) = pool_is_aqueous(ipool) 

          if (this%downstream_nc(icount,jcount) < 0.d0) then
            option%io_buffer = 'For CLM_Decomp reactions, downstream pools ' // &
             'must have a constant C:N ratio (i.e. C and N are not tracked ' // &
             ' individually.  Therefore, pool "' // &
             trim(cur_pool%name) // &
            '" may not be used as a downstream pool.'
            call printErrMsg(option)
          endif
        endif
      endif

      cur_pool => cur_pool%next

    enddo

    this%rate_constant(icount) = cur_rxn%rate_constant
    cur_rxn => cur_rxn%next
  enddo 
  
  deallocate(pool_names)
  call DeallocateArray(pool_is_aqueous)
  call DeallocateArray(species_id_pool_c)
  call DeallocateArray(species_id_pool_n)

! set stoichiometric coefficients for som decomposition reactions  
  do icount = 1, this%nrxn
     if(this%is_litter_decomp(icount)) then
        cycle
     else
!    calculate respiration factor
       stoich_c = 1.0d0
       stoich_n = this%upstream_nc(icount)

       do jcount = 1, this%n_downstream_pools(icount)
          stoich_c = stoich_c - this%downstream_stoich(icount, jcount)
          stoich_n = stoich_n - this%downstream_stoich(icount, jcount) * &
                              this%downstream_nc(icount, jcount)
       enddo

       if(stoich_c < -1.0d-10) then
         option%io_buffer = 'CLM_Decomp SOM decomposition reaction has negative respiration fraction!'
         call printErrMsg(option)
       endif

       this%mineral_c_stoich(icount) = stoich_c
       this%mineral_n_stoich(icount) = stoich_n

     endif
  enddo

  word = 'CO2(aq)*'
  this%species_id_co2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'NH4+'
  this%species_id_nh3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if(this%species_id_nh3 <= 0) then
    option%io_buffer = 'NH4+ is not defined in the database for CLM_Decomp'
    call printErrMsg(option)
  endif

  word = 'NO3-'
  this%species_id_no3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  
  word = 'N2O(aq)'
  this%species_id_n2o = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
 
  word = 'Bacteria'
  this%species_id_bacteria = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)

  word = 'Fungi'
  this%species_id_fungi = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
end subroutine CLM_Decomp_Setup

! ************************************************************************** !

subroutine CLM_Decomp_React(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                       global_auxvar,porosity,volume,reaction,option, local_id)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
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

  class(reaction_sandbox_clm_decomp_type) :: this
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
  PetscInt, parameter :: iphase = 1

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
!  PetscReal :: stoich_un, stoich_dc, stoich_mineraln, stoich_co2
  PetscReal :: stoich_c, stoich_n

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

!  PetscReal, parameter :: rpi = 3.14159265358979323846

  PetscReal :: c_uc, c_un
  PetscReal :: nc_bacteria, nc_fungi

  c_nh3     = rt_auxvar%total(this%species_id_nh3, iphase)
  temp_real = c_nh3 - this%x0eps + this%half_saturation_nh3
  f_nh3     = (c_nh3 - this%x0eps) / temp_real 
  d_nh3     = this%half_saturation_nh3 / temp_real / temp_real

  if (this%species_id_no3 > 0) then
      c_no3 = rt_auxvar%total(this%species_id_no3, iphase)
      temp_real = c_no3 - this%x0eps + this%half_saturation_no3
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

    ispec_uc = this%upstream_c_id(irxn)
    if(this%upstream_is_aqueous(irxn)) then
      c_uc = rt_auxvar%total(ispec_uc, iphase)
      c_uc = theta * 1000.0d0 * c_uc    ! from mol/L -> mol/m3
      ires_uc = ispec_uc
    else
      c_uc = rt_auxvar%immobile(ispec_uc)
      ires_uc = reaction%offset_immobile + ispec_uc
    endif

    if(c_uc < 1.0d-20) cycle

!   for litter decomposition reactions, stoich needs to be calculated on the fly
    if(this%is_litter_decomp(irxn)) then

      ispec_un = this%upstream_n_id(irxn)
      ires_un = ispec_un + reaction%offset_immobile

      if(ispec_un > 0) then
         c_un = rt_auxvar%immobile(ispec_un)
         this%upstream_nc(irxn) = c_un / c_uc
      endif

      if(this%litter_decomp_type == LITTER_DECOMP_CLMCN) then

!    calculate respiration factor (CO2 stoichiometry)
        stoich_c = 1.0d0

        do j = 1, this%n_downstream_pools(irxn)
           stoich_c = stoich_c - this%downstream_stoich(irxn, j)
        enddo

        if(stoich_c < 0.0d0) then
           option%io_buffer = 'CLM_Decomp litter decomposition reaction has negative respiration fraction!'
           call printErrMsg(option)
        endif

        this%mineral_c_stoich(irxn) = stoich_c

!     calculate N stoichiometry 
        stoich_n = this%upstream_nc(irxn)

        do j = 1, this%n_downstream_pools(irxn)
           stoich_n = stoich_n - this%downstream_stoich(irxn, j) * &
                              this%downstream_nc(irxn, j)
        enddo

        this%mineral_n_stoich(irxn) = stoich_n

      elseif(this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then

!    Sinsabaugh et al. 2013 Ecology Letters, 16, 930-939
        resp_frac = CN_ratio_microbe * this%upstream_nc(irxn) !c_un/c_uc

        if(resp_frac > CUE_max) then
          resp_frac = CUE_max
        endif

!    c pools     
        this%mineral_c_stoich(irxn) = resp_frac
       
        if(this%n_downstream_pools(irxn) .ne. 2) then
           option%io_buffer = 'CLM_Microbe litter decomposition reaction more than 2 (bacteria and fungi pools)!'
           call printErrMsg(option)
        endif

        do i = 1, this%n_downstream_pools(irxn)
          if(this%downstream_id(irxn, i) == this%species_id_bacteria) then
            this%downstream_stoich(irxn, i) = fraction_bacteria * (1.0d0  - resp_frac)
          else
            this%downstream_stoich(irxn, i) = (1.0d0 - fraction_bacteria) * (1.0d0  - resp_frac)
          endif
        enddo

        stoich_n = this%upstream_nc(irxn)

        do j = 1, this%n_downstream_pools(irxn)
           stoich_n = stoich_n - this%downstream_stoich(irxn, j) * &
                              this%downstream_nc(irxn, j)
        enddo

        this%mineral_n_stoich(irxn) = stoich_n

      endif
    endif

    ! residual units: (mol/sec) = (m^3 bulk/s) * (mol/m^3 bulk)
    rate     = scaled_rate_const * (c_uc - this%x0eps)
    drate_uc = scaled_rate_const 

!   NH3 limiting
    if(this%mineral_n_stoich(irxn) < 0.0d0) then
       drate_nh3 = rate * d_nh3 
       rate      = rate * f_nh3
       drate_uc  = drate_uc * f_nh3
    endif 

    ! calculation of residual
    ! carbon
    Residual(ires_co2) = Residual(ires_co2) - this%mineral_c_stoich(irxn) * rate
    
    ! nitrogen
    Residual(ires_nh3) = Residual(ires_nh3) - this%mineral_n_stoich(irxn) * rate

    ! upstream c
    Residual(ires_uc) = Residual(ires_uc) - (-1.d0) * rate

    ! upstream n
    if(this%is_litter_decomp(irxn)) then
      Residual(ires_un) = Residual(ires_un) - (-1.d0) * this%upstream_nc(irxn) * rate
    endif
    
!   downstream pools
    do j = 1, this%n_downstream_pools(irxn)
       ispec_d = this%downstream_id(irxn, j)
       if(this%downstream_is_aqueous(irxn, j)) then
         ires_d = ispec_d
       else
         ires_d = reaction%offset_immobile + ispec_d
       endif
       if(ispec_d > 0) then
          Residual(ires_d) = Residual(ires_d) - this%downstream_stoich(irxn, j) * rate
       endif
    enddo

    net_n_mineralization_rate=net_n_mineralization_rate+this%mineral_n_stoich(irxn)*rate

!   start residual calculation for N immobilization reaction with NO3 uptake
!   if nitrate is available, N immobilization decomposition reactions occurs
!   with rate depending on NH3, with reduced rate if NH3 is abundent  
    if(this%species_id_no3 > 0 .and. this%mineral_n_stoich(irxn) < 0.d0) then

       rate_no3     = scaled_rate_const * (c_uc - this%x0eps) * f_no3 * f_nh3_inhibit
       drate_uc_no3 = scaled_rate_const * 1.0d0 * f_no3 * f_nh3_inhibit
       drate_no3    = scaled_rate_const * c_uc  * d_no3 * f_nh3_inhibit
       drate_nh3_no3= scaled_rate_const * c_uc  * f_no3 * d_nh3_inhibit

    ! residuals
    ! carbon
       Residual(ires_co2) = Residual(ires_co2) - this%mineral_c_stoich(irxn) * rate_no3
    
    ! nitrogen
       Residual(ires_no3) = Residual(ires_no3) - this%mineral_n_stoich(irxn) * rate_no3

    ! upstream c 
       Residual(ires_uc) = Residual(ires_uc) - (-1.d0) * rate_no3
    
    ! upstream n
       if(this%is_litter_decomp(irxn)) then
         Residual(ires_un) = Residual(ires_un)+this%upstream_nc(irxn)*rate_no3
       endif
    
!   downstream pools
       do j = 1, this%n_downstream_pools(irxn)
         ispec_d = this%downstream_id(irxn, j)
         if(this%downstream_is_aqueous(irxn, j)) then
           ires_d = ispec_d
         else
           ires_d = reaction%offset_immobile + ispec_d
         endif
         if(ispec_d > 0) then
            Residual(ires_d) = Residual(ires_d) - this%downstream_stoich(irxn, j) * rate_no3
         endif
       enddo

       net_n_mineralization_rate=net_n_mineralization_rate+this%mineral_n_stoich(irxn)*rate_no3
    endif 
!   end residual calculation for N immobilization reaction with NO3 uptake

    if (compute_derivative) then
 
! with respect to upstream c    
      ! carbon
      if(this%upstream_is_aqueous(irxn)) then
        Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
          this%mineral_c_stoich(irxn) * drate_uc * &
          rt_auxvar%aqueous%dtotal(this%species_id_co2,ispec_uc,iphase)

      else
        Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
          this%mineral_c_stoich(irxn) * drate_uc

        if(this%is_litter_decomp(irxn) .and. &
           this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE .and. &
           resp_frac < CUE_max) then
! dRco2/dLit1C = dcR/dLit1C = cdR/dLit1C + R dc/dLit1C 
! c = min(CUEmax, Lit1N/Lit1C*CN_ratio_microbe)
! dc/dLit1C =  - Lit1N/Lit1C^2 CN_ratio_microbe  
! R dc/dLit1C =  - R Lit1N/Lit1C^2 CN_ratio_microbe  = -dR/dLit1C u CN_ratio_microbe
           Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) + &
            drate_uc * this%upstream_nc(irxn) * CN_ratio_microbe 
        endif
      endif

      ! nitrogen
         ! dR_N/dC_u = d nR/dC_u = dn/dC_u R + n dR/dC_u
         ! first, n dR/dC_u
      if(this%upstream_is_aqueous(irxn)) then
        Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - &
          this%mineral_n_stoich(irxn) * drate_uc * &
          rt_auxvar%aqueous%dtotal(this%species_id_nh3,ispec_uc,iphase)
      else
        Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - &
          this%mineral_n_stoich(irxn) * drate_uc
      endif

      if(this%is_litter_decomp(irxn)) then
         ! litter pool is immobile
         ! second, dn/dC_u R, n = u - (1 - f)d, dn/dC_u = du/dC_u = -N_u/C_u^2 = -u/C_u
         ! dn/dC_u R = -u R/C_u = -u drate_uc 
         if(this%litter_decomp_type == LITTER_DECOMP_CLMCN) then
           Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_uc
         else  ! litter decomposition CLM-MICROBE
           if(resp_frac < CUE_max) then
! dRN/dLit1C = dnR/dLit1C = ndR/dLit1C + R dn/dLit1C 
! n = u - b nb - f nf = u - g(1-c)n_b - (1 - g)(1 - c) n_f
! dn/dLit1C =  - Lit1N/Lit1C^2 + (g n_b + (1-g)n_f)dc/dLit1C
! R dn/dLit1C =  - dR/dLit1C u - (g n_b + (1 - g) n_f) dR/dLit1C u CN_ratio_microbe
! R dn/dLit1C =  - dR/dLit1C u (1 + (g n_b + (1 - g) n_f) CN_ratio_microbe)
             do i = 1, this%n_downstream_pools(irxn)
                if(this%downstream_id(irxn, i) == this%species_id_bacteria) then
                  nc_bacteria = this%downstream_nc(irxn, i)
                else
                  nc_fungi = this%downstream_nc(irxn, i)
                endif
             enddo

             Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) + &
               this%upstream_nc(irxn) * drate_uc * (1.0d0 + &
               (fraction_bacteria * nc_bacteria + &
               (1.0d0 - fraction_bacteria) * nc_fungi) * CN_ratio_microbe)

           else
             Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - &
              (-1.d0) * this%upstream_nc(irxn) * drate_uc
           endif
         endif
      endif

      ! upstream C pool
      if(this%upstream_is_aqueous(irxn)) then
        Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) - (-1.d0) * drate_uc * &
          rt_auxvar%aqueous%dtotal(ispec_uc,ispec_uc,iphase)
      else
        Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) - (-1.d0) * drate_uc
      endif

      ! upstream N pool
      if(this%is_litter_decomp(irxn)) then
        ! litter pools are immobile
        ! R_Nu = Nu/Cu * R_Cu, dR_Nu/dCu = Nu/Cu dR_Cu/dCu - Nu/Cu^2 R_Cu
        ! R_Cu = k Cu,         dR_Cu/dCu = k
        ! dR_Nu/dCu = Nu/Cu k - Nu/Cu kCu/Cu = 0
         if(this%litter_decomp_type == LITTER_DECOMP_CLMCN) then
!        Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) - &
!          (-1.d0) * stoich_un * drate_uc 
         else   ! litter decomposition CLM-MICROBE

         endif 
      endif

!   downstream pools
      do j = 1, this%n_downstream_pools(irxn)
         ispec_d = this%downstream_id(irxn, j)
         if(ispec_d < 0) then
           cycle
         endif

         if(this%downstream_is_aqueous(irxn, j)) then
           ires_d = ispec_d
         else
           ires_d = reaction%offset_immobile + ispec_d
         endif
         
         if(this%upstream_is_aqueous(irxn) .and. &
            this%downstream_is_aqueous(irxn, j)) then
               Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - &
                 this%downstream_stoich(irxn, j) * drate_uc * &
                 rt_auxvar%aqueous%dtotal(ispec_d,ispec_uc,iphase)
         else
            Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - &
               this%downstream_stoich(irxn, j) * drate_uc

            if(this%is_litter_decomp(irxn) .and. &
               this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE .and. &
               resp_frac < CUE_max) then
               if(ispec_d == this%species_id_bacteria) then
                 Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc)-drate_uc* &
                 fraction_bacteria * this%upstream_nc(irxn) * CN_ratio_microbe
               elseif(ispec_d == this%species_id_fungi) then
                 Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc)-drate_uc* &
               (1.0d0-fraction_bacteria)*this%upstream_nc(irxn)*CN_ratio_microbe
               else
               endif
            endif
         endif
      enddo

!with respect to upstream n (due to variable CN ratio)
      if(this%is_litter_decomp(irxn)) then
         if(this%litter_decomp_type == LITTER_DECOMP_CLMCN) then
        ! upstream n, dR_Nu/dNu = d uR/dNu = R/Cu = dR/dCu
           Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) -(-1.d0)*drate_uc
        ! mineral N, dR_N/dNu = d nR/du = Rd (u - (1-f)d)/dNu = dR/dCu 
           Jacobian(ires_nh3,ires_un) = Jacobian(ires_nh3,ires_un) - drate_uc
         else   ! litter decomposition CLM-MICROBE

         endif 
      endif

!with respect to nh3
      if(this%mineral_n_stoich(irxn) < 0.0d0) then
        ! carbon
        Jacobian(ires_co2,ires_nh3) = Jacobian(ires_co2,ires_nh3) - &
          this%mineral_c_stoich(irxn) * drate_nh3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_nh3,iphase)
  
        ! nitrogen
        Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) - &
          this%mineral_n_stoich(irxn) * drate_nh3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_nh3,this%species_id_nh3,iphase)
  
        ! upstream C pool
        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_uc,ires_nh3) = Jacobian(ires_uc,ires_nh3) - &
            (-1.d0) * drate_nh3 * &
          rt_auxvar%aqueous%dtotal(ispec_uc,this%species_id_nh3,iphase)
        else
          Jacobian(ires_uc,ires_nh3) = Jacobian(ires_uc,ires_nh3) - &
            (-1.d0) * drate_nh3
        endif
 
        ! upstream N pool
        if(this%is_litter_decomp(irxn)) then
          Jacobian(ires_un,ires_nh3) = Jacobian(ires_un,ires_nh3) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_nh3 
        endif

!   downstream pools
        do j = 1, this%n_downstream_pools(irxn)
           ispec_d = this%downstream_id(irxn, j)
           if(ispec_d < 0) then
             cycle
           endif
           if(this%downstream_is_aqueous(irxn, j)) then
              ires_d = ispec_d
              Jacobian(ires_d,ires_nh3) = Jacobian(ires_d,ires_nh3) - &
                  this%downstream_stoich(irxn, j) * drate_nh3 * &
              rt_auxvar%aqueous%dtotal(ispec_d,this%species_id_nh3,iphase)
           else
              ires_d = reaction%offset_immobile + ispec_d
              Jacobian(ires_d,ires_nh3) = Jacobian(ires_d,ires_nh3) - &
                  this%downstream_stoich(irxn, j) * drate_nh3
           endif
        enddo
      endif

!   start jacobian calculation for N immobilization reaction with NO3 uptake
      if(this%species_id_no3 > 0 .and. this%mineral_n_stoich(irxn) < 0.d0) then
! with respect to upstream c    
      ! carbon
        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
            this%mineral_c_stoich(irxn) * drate_uc_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_co2,ispec_uc,iphase)
        else
          Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
            this%mineral_c_stoich(irxn) * drate_uc_no3
        endif

      ! nitrogen
        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - &
            this%mineral_n_stoich(irxn) * drate_uc_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_no3,ispec_uc,iphase)
        else
          Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - &
            this%mineral_n_stoich(irxn) * drate_uc_no3
        endif

        if(this%is_litter_decomp(irxn)) then
         ! litter pool is immobile
         ! second, dn/dC_u R, n = u - (1 - f)d, dn/dC_u = du/dC_u = -N_u/C_u^2 = -u/C_u
         ! dn/dC_u R = -u R/C_u = -u drate_uc 
          if(this%litter_decomp_type == LITTER_DECOMP_CLMCN) then
            Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - &
             (-1.d0) * this%upstream_nc(irxn) * drate_uc_no3
          else  ! litter decomposition CLM-MICROBE

          endif
        endif

      ! upstream C pool
        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) + drate_uc_no3 * &
            rt_auxvar%aqueous%dtotal(ispec_uc,ispec_uc,iphase)
        else
          Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) + drate_uc_no3
        endif

      ! upstream N pool
        if(this%is_litter_decomp(irxn)) then
        ! litter pools are immobile
        ! R_Nu = Nu/Cu * R_Cu, dR_Nu/dCu = Nu/Cu dR_Cu/dCu - Nu/Cu^2 R_Cu
        ! R_Cu = k Cu,         dR_Cu/dCu = k
        ! dR_Nu/dCu = Nu/Cu k - Nu/Cu kCu/Cu = 0
          if(this%litter_decomp_type == LITTER_DECOMP_CLMCN) then
!        Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) - &
!          (-1.d0) * stoich_un * drate_uc 
          else   ! litter decomposition CLM-MICROBE

          endif 
        endif

!   downstream pools
        do j = 1, this%n_downstream_pools(irxn)
           ispec_d = this%downstream_id(irxn, j)
           if(ispec_d < 0) then
             cycle
           endif

           if(this%downstream_is_aqueous(irxn, j)) then
             ires_d = ispec_d
           else
             ires_d = reaction%offset_immobile + ispec_d
           endif
         
           if(this%upstream_is_aqueous(irxn) .and. &
              this%downstream_is_aqueous(irxn, j)) then
               Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - &
                 this%downstream_stoich(irxn, j) * drate_uc_no3 * &
                 rt_auxvar%aqueous%dtotal(ispec_d,ispec_uc,iphase)
           else
              Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - &
                 this%downstream_stoich(irxn, j) * drate_uc_no3
           endif
        enddo

!with respect to upstream n (due to variable CN ratio)
        if(this%is_litter_decomp(irxn)) then
          if(this%litter_decomp_type == LITTER_DECOMP_CLMCN) then
        ! upstream n, dR_Nu/dNu = d uR/dNu = R/Cu = dR/dCu
            Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) -(-1.d0)*drate_uc_no3
        ! mineral N, dR_N/dNu = d nR/du = Rd (u - (1-f)d)/dNu = dR/dCu 
            Jacobian(ires_no3,ires_un) = Jacobian(ires_no3,ires_un) - drate_uc_no3
          else   ! litter decomposition CLM-MICROBE

          endif 
        endif

!with respect to no3
!        if(stoich_mineraln < 0.0d0) then
        ! carbon
        Jacobian(ires_co2,ires_no3) = Jacobian(ires_co2,ires_no3) - &
            this%mineral_c_stoich(irxn) * drate_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_no3,iphase)
  
        ! nitrogen
        Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) - &
            this%mineral_n_stoich(irxn) * drate_no3 * &
            rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_no3,iphase)
  
        ! upstream C pool
        if(this%upstream_is_aqueous(irxn)) then
           Jacobian(ires_uc,ires_no3) = Jacobian(ires_uc,ires_no3) - &
            (-1.d0) * drate_no3 * &
          rt_auxvar%aqueous%dtotal(ispec_uc,this%species_id_no3,iphase)
        else
           Jacobian(ires_uc,ires_no3) = Jacobian(ires_uc,ires_no3) - &
            (-1.d0) * drate_no3
        endif
 
        ! upstream N pool
        if(this%is_litter_decomp(irxn)) then
           Jacobian(ires_un,ires_no3) = Jacobian(ires_un,ires_no3) - &
             (-1.d0) * this%upstream_nc(irxn) * drate_no3
        endif
  
!   downstream pools
        do j = 1, this%n_downstream_pools(irxn)

           ispec_d = this%downstream_id(irxn, j)

           if(ispec_d < 0) then
             cycle
           endif

           if(this%downstream_is_aqueous(irxn, j)) then
              ires_d = ispec_d
              Jacobian(ires_d,ires_no3) = Jacobian(ires_d,ires_no3) - &
                  this%downstream_stoich(irxn, j) * drate_no3 * &
              rt_auxvar%aqueous%dtotal(ispec_d,this%species_id_no3,iphase)
           else
              ires_d = reaction%offset_immobile + ispec_d
              Jacobian(ires_d,ires_no3) = Jacobian(ires_d,ires_no3) - &
                  this%downstream_stoich(irxn, j) * drate_no3
           endif
        enddo

! with respect to nh3
        ! carbon
        Jacobian(ires_co2,ires_nh3) = Jacobian(ires_co2,ires_nh3) - &
          this%mineral_c_stoich(irxn) * drate_nh3_no3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_nh3,iphase)
  
        ! nitrogen
        Jacobian(ires_no3,ires_nh3) = Jacobian(ires_no3,ires_nh3) - &
          this%mineral_n_stoich(irxn) * drate_nh3_no3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_nh3,iphase)
  
        ! upstream C pool
        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_uc,ires_nh3) = Jacobian(ires_uc,ires_nh3) - &
          (-1.d0) * drate_nh3_no3 * &
          rt_auxvar%aqueous%dtotal(ispec_uc,this%species_id_nh3,iphase)
        else
          Jacobian(ires_uc,ires_nh3) = Jacobian(ires_uc,ires_nh3) - &
          (-1.d0) * drate_nh3_no3
        endif
  
        ! upstream N pool
        if(this%is_litter_decomp(irxn)) then
           Jacobian(ires_un,ires_nh3) = Jacobian(ires_un,ires_nh3) - &
             (-1.d0) * this%upstream_nc(irxn) * drate_nh3_no3 
        endif
  
!   downstream pools
        do j = 1, this%n_downstream_pools(irxn)

           ispec_d = this%downstream_id(irxn, j)

           if(ispec_d < 0) then
             cycle
           endif

           if(this%downstream_is_aqueous(irxn, j)) then
              ires_d = ispec_d
              Jacobian(ires_d,ires_nh3) = Jacobian(ires_d,ires_nh3) - &
                  this%downstream_stoich(irxn, j) * drate_nh3_no3 * &
              rt_auxvar%aqueous%dtotal(ispec_d,this%species_id_nh3,iphase)
           else
              ires_d = reaction%offset_immobile + ispec_d
              Jacobian(ires_d,ires_nh3) = Jacobian(ires_d,ires_nh3) - &
                  this%downstream_stoich(irxn, j) * drate_nh3_no3
           endif
        enddo
      endif  ! jacobian calculation for N immobilization reaction with NO3 uptake

!  write(*,*)(Residual(j), j = 1, reaction%ncomp)
!  do i = 1, reaction%ncomp
!     write(*,*)(Jacobian(i, j), j = 1, reaction%ncomp)
!  enddo

    endif    ! jacobian calculation
 
  enddo      !reactions

  !if(this%no_n2o_decomp) then
!     return
 ! endif

  if(this%species_id_n2o > 0 .and. net_n_mineralization_rate > 0.0d0) then
  ! temperature response function (Parton et al. 1996)
#ifdef CLM_PFLOTRAN
    f_t = -0.06d0 + 0.13d0 * exp( 0.07d0 * tc )
    f_w = ((1.27d0 - saturation)/0.67d0)**(3.1777d0) * &
        ((saturation - 0.0012d0)/0.5988d0)**2.84d0  
    ph = 6.5d0
    f_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + ph))/rpi
#else
    f_t = 1.0d0
    f_w = 1.0d0
    f_ph = 1.0d0
#endif

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
        Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) + drate_n2o
        Jacobian(ires_n2o,ires_nh3) = Jacobian(ires_n2o,ires_nh3) - &
                                      0.5d0 * drate_n2o
      endif
    endif
  endif

end subroutine CLM_Decomp_React

! ************************************************************************** !
!
! CLM_DecompDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Guoping Tang
! ************************************************************************** !
subroutine CLM_Decomp_Destroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_clm_decomp_type) :: this
  
  type(pool_type), pointer :: cur_pool, prev_pool
  type(clm_decomp_reaction_type), pointer :: cur_reaction, prev_reaction
  
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

    cur_pool => cur_reaction%downstream_pools  
    do
      if (.not.associated(cur_pool)) exit
      prev_pool => cur_pool
      cur_pool => cur_pool%next
      deallocate(prev_pool)
      nullify(prev_pool)
    enddo

    prev_reaction => cur_reaction
    cur_reaction => cur_reaction%next

    deallocate(prev_reaction)
    nullify(prev_reaction)
  enddo
  
  call DeallocateArray(this%pool_nc_ratio)
  call DeallocateArray(this%rate_constant)
  call DeallocateArray(this%is_litter_decomp)
  call DeallocateArray(this%upstream_c_id)
  call DeallocateArray(this%upstream_n_id)
  call DeallocateArray(this%upstream_nc)
  call DeallocateArray(this%upstream_is_aqueous)
  call DeallocateArray(this%downstream_id)
  call DeallocateArray(this%downstream_stoich)
  call DeallocateArray(this%downstream_is_aqueous)
  call DeallocateArray(this%mineral_c_stoich) 
  call DeallocateArray(this%mineral_n_stoich) 
 
end subroutine CLM_Decomp_Destroy

end module Reaction_Sandbox_CLM_Decomp_class
