module Reaction_Sandbox_SomDec_class

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  use CLM_RspFuncs_module

! -----------------------------------------------------------------------------
! description
! -----------------------------------------------------------------------------

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_somdec_type

    PetscInt :: temperature_response_function
    PetscInt :: moisture_response_function
    PetscReal :: Q10
    PetscReal :: decomp_depth_efolding
    PetscReal :: half_saturation_nh4
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh4_no3
    PetscReal :: n2o_frac_mineralization    ! fraction of n2o from net N mineralization
    PetscReal :: x0eps

    PetscInt :: npool                       ! litter or variable CN ration pools
    PetscReal, pointer :: pool_nc_ratio(:)         ! NC ratio in mole  

    PetscInt :: nrxn
    PetscReal, pointer :: rate_constant(:)         !nrxn
    PetscReal, pointer :: rate_decomposition(:)    !nrxn
    PetscReal, pointer :: rate_ad_factor(:)        !nrxn

    PetscInt,  pointer :: upstream_c_id(:)         !nrxn
    PetscInt,  pointer :: upstream_n_id(:)         !nrxn
    PetscReal, pointer :: upstream_nc(:)           !nrxn
    PetscBool, pointer :: upstream_is_aqueous(:)   !nrxn
    PetscBool, pointer :: upstream_is_varycn(:)

    PetscInt,  pointer :: n_downstream_pools(:)   !nrxn by maximum # of downstream pools
    PetscInt,  pointer :: downstream_id(:,:)      !nrxn by maximum # of downstream pools
    PetscBool, pointer :: downstream_is_aqueous(:,:) !nrxn by maximum # of downstream pools
    PetscReal, pointer :: downstream_stoich(:,:)  !nrxn by maximum # of downstream pools
    PetscReal, pointer :: downstream_nc(:,:)      !nrxn by maximum # of downstream pools
    PetscBool, pointer :: downstream_is_varycn(:,:)
    PetscReal, pointer :: mineral_c_stoich(:)     !nrxn by maximum # of downstream pools
    PetscReal, pointer :: mineral_n_stoich(:)     !nrxn by maximum # of downstream pools

    PetscInt :: species_id_co2
    PetscInt :: species_id_nh4
    PetscInt :: species_id_no3
    PetscInt :: species_id_n2o

    PetscInt :: species_id_hrimm
    PetscInt :: species_id_nmin
    PetscInt :: species_id_nimm
    PetscInt :: species_id_ngasmin
    PetscInt :: species_id_proton

    type(pool_type), pointer :: pools
    type(somdec_reaction_type), pointer :: reactions
  contains
    procedure, public :: ReadInput => SomDecRead
    procedure, public :: Setup => SomDecSetup
    procedure, public :: Evaluate => SomDecReact
    procedure, public :: Destroy => SomDecDestroy
  end type reaction_sandbox_somdec_type
  
  type :: pool_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: stoich
    PetscReal :: nc_ratio
    type(pool_type), pointer :: next
  end type pool_type
  
  type :: somdec_reaction_type
    character(len=MAXWORDLENGTH) :: upstream_pool_name
    type(pool_type), pointer :: downstream_pools
    PetscReal :: rate_constant
    PetscReal :: rate_decomposition
    PetscReal :: rate_ad_factor
    type(somdec_reaction_type), pointer :: next
  end type somdec_reaction_type
  
  public :: SomDecCreate

contains

! ************************************************************************** !

function SomDecCreate()
  ! 
  ! Allocates SomDec reaction object.
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
  ! 

  implicit none
  
  type(reaction_sandbox_somdec_type), pointer :: SomDecCreate
  
  allocate(SomDecCreate)

#ifdef CLM_PFLOTRAN
  SomDecCreate%temperature_response_function=TEMPERATURE_RESPONSE_FUNCTION_CLM4
  SomDecCreate%moisture_response_function = MOISTURE_RESPONSE_FUNCTION_CLM4
#endif

  SomDecCreate%Q10 = 1.5d0
  SomDecCreate%decomp_depth_efolding = 0.d0      ! non-positive value will turn off this option
  SomDecCreate%half_saturation_nh4 = 1.0d-15
  SomDecCreate%half_saturation_no3 = 1.0d-15
  SomDecCreate%inhibition_nh4_no3 = 1.0d0
  SomDecCreate%n2o_frac_mineralization = 0.02d0  ! Parton et al. 2001
  SomDecCreate%x0eps = 1.0d-20

  SomDecCreate%npool = 0
  nullify(SomDecCreate%pool_nc_ratio)

  SomDecCreate%nrxn = 0
  nullify(SomDecCreate%rate_constant)
  nullify(SomDecCreate%rate_decomposition)
  nullify(SomDecCreate%rate_ad_factor)
  nullify(SomDecCreate%upstream_is_varycn)
  nullify(SomDecCreate%upstream_c_id)
  nullify(SomDecCreate%upstream_n_id)
  nullify(SomDecCreate%upstream_nc)
  nullify(SomDecCreate%upstream_is_aqueous)
  
  nullify(SomDecCreate%n_downstream_pools)
  nullify(SomDecCreate%downstream_id)
  nullify(SomDecCreate%downstream_is_aqueous)
  nullify(SomDecCreate%downstream_stoich)
  nullify(SomDecCreate%downstream_is_varycn)
  nullify(SomDecCreate%mineral_c_stoich)
  nullify(SomDecCreate%mineral_n_stoich)

  SomDecCreate%species_id_co2 = 0
  SomDecCreate%species_id_nh4 = 0
  SomDecCreate%species_id_no3 = 0
  SomDecCreate%species_id_n2o = 0
  SomDecCreate%species_id_hrimm = 0
  SomDecCreate%species_id_nmin = 0
  SomDecCreate%species_id_nimm = 0
  SomDecCreate%species_id_ngasmin = 0

  nullify(SomDecCreate%next)

  nullify(SomDecCreate%pools)
  nullify(SomDecCreate%reactions)

end function SomDecCreate

! ************************************************************************** !

subroutine SomDecRead(this,input,option)
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
  
  class(reaction_sandbox_somdec_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  
  type(pool_type), pointer :: new_pool, prev_pool
  type(pool_type), pointer :: new_pool_rxn, prev_pool_rxn
  type(somdec_reaction_type), pointer :: new_reaction, prev_reaction
  
  PetscReal :: rate_constant, turnover_time
  PetscReal :: rate_decomposition, rate_ad_factor
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
                       'CHEMISTRY,REACTION_SANDBOX,SomDecomp')
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
            'CHEMISTRY,REACTION_SANDBOX,SomDecomp,TEMPERATURE RESPONSE FUNCTION')
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
                       'REACTION_SANDBOX_SomDecomp,TEMPERATURE RESPONSE FUNCTION')
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,SomDec,' // &
                                'TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
         enddo 
      case('MOISTURE_RESPONSE_FUNCTION')
        do
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit

         call InputReadWord(input,option,word,PETSC_TRUE)
         call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,SomDecomp,MOISTURE RESPONSE FUNCTION')
         call StringToUpper(word)   

            select case(trim(word))
              case('CLM4')
                  this%moisture_response_function = MOISTURE_RESPONSE_FUNCTION_CLM4    
              case('DLEM')
                  this%moisture_response_function = MOISTURE_RESPONSE_FUNCTION_DLEM    
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,SomDecomp,TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
         enddo 

#endif

     case('X0EPS')
         call InputReadDouble(input,option,this%x0eps)
         call InputErrorMsg(input,option,'x0eps', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDecomp,REACTION')

     case('AMMONIUM_HALF_SATURATION')
         call InputReadDouble(input,option,this%half_saturation_nh4)
         call InputErrorMsg(input,option,'ammonium half saturation', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDecomp,REACTION')

     case('NITRATE_HALF_SATURATION')
         call InputReadDouble(input,option,this%half_saturation_no3)
         call InputErrorMsg(input,option,'nitrate half saturation', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDecomp,REACTION')

     case('AMMONIUM_INHIBITION_NITRATE')
         call InputReadDouble(input,option,this%inhibition_nh4_no3)
         call InputErrorMsg(input,option,'ammonium inhibition on nitrate immobilization', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDecomp,REACTION')
         if (this%inhibition_nh4_no3<0.d-20) then
           option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,SomDecomp,' // &
             'AMMONIUM_INHIBITION_NITRATE cannot be too small to close to 0'
           call printErrMsg(option)
         endif

     case('N2O_FRAC_MINERALIZATION')
         call InputReadDouble(input,option,this%n2o_frac_mineralization)
         call InputErrorMsg(input,option,'n2o fraction from mineralization', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDecomp,REACTION')

     case('DECOMP_DEPTH_EFOLDING')
         call InputReadDouble(input,option,this%decomp_depth_efolding)
         call InputErrorMsg(input,option,'decomp_depth_efolding', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDecomp,REACTION')

     case('POOLS')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit   

          allocate(new_pool)
          new_pool%name = ''
          new_pool%nc_ratio = UNINITIALIZED_DOUBLE
          nullify(new_pool%next)

          call InputReadWord(input,option,new_pool%name,PETSC_TRUE)
          call InputErrorMsg(input,option,'pool name', &
            'CHEMISTRY,REACTION_SANDBOX,SomDecomp,POOLS')
          call InputReadDouble(input,option,temp_real)
          if (InputError(input)) then
            new_pool%nc_ratio = UNINITIALIZED_DOUBLE
          else
            ! convert CN ratio from mass C/mass N to mol N/mol C
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
        new_reaction%rate_constant = UNINITIALIZED_DOUBLE
        new_reaction%rate_decomposition = UNINITIALIZED_DOUBLE
        new_reaction%rate_ad_factor = 1.0d0
        nullify(new_reaction%downstream_pools)
        nullify(new_reaction%next)
        
        ! need to set these temporarily in order to check that they
        ! are not both set.
        turnover_time = 0.d0
        rate_constant = -1.d0
        rate_decomposition = -1.d0
        rate_ad_factor = 1.d0
        
        do 
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,REACTION_SANDBOX,SomDec,REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('UPSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%upstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'upstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDec,REACTION')
            case('DOWNSTREAM_POOL')
              allocate(new_pool_rxn)
              new_pool_rxn%name = ''
              new_pool_rxn%stoich = 0.d0
              nullify(new_pool_rxn%next)

              call InputReadWord(input,option, &
                                 new_pool_rxn%name,PETSC_TRUE)
              call InputErrorMsg(input,option,'downstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDec,REACTION')
              call InputReadDouble(input,option,new_pool_rxn%stoich)
              call InputErrorMsg(input,option,'Downstream pool stoich', 'CHEMISTRY,' // &
                  'REACTION_SANDBOX_SomDec,REACTION')

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
                     'CHEMISTRY,REACTION_SANDBOX,SomDec,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'SomDec RATE CONSTANT UNITS'
                call InputDefaultMsg(input,option)
              else              
                rate_constant = rate_constant * &
                  UnitsConvertToInternal(word,option)
              endif
            case('TURNOVER_TIME')
              call InputReadDouble(input,option,turnover_time)
              call InputErrorMsg(input,option,'turnover time', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDec,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'SomDec TURNOVER TIME UNITS'
                call InputDefaultMsg(input,option)
              else              
                turnover_time = turnover_time * &
                  UnitsConvertToInternal(word,option)
              endif
            case('RATE_DECOMPOSITION')
              call InputReadDouble(input,option,rate_decomposition)
              call InputErrorMsg(input,option,'rate for decomposition', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDec,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'SomDec RATE DECOMPOSITION UNITS'
                call InputDefaultMsg(input,option)
              else
                rate_decomposition = rate_decomposition * &
                  UnitsConvertToInternal(word,option)
              endif
            case('RATE_AD_FACTOR')
              call InputReadDouble(input,option,rate_ad_factor)
              call InputErrorMsg(input,option,'Accelerated rate factor', &
                     'CHEMISTRY,REACTION_SANDBOX,SomDec,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,SomDec,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
              call printErrMsg(option)
          end select
        enddo
        
        ! check to ensure that one of turnover time or rate constant is set.
        if ( (turnover_time > 0.d0 .and. rate_constant >= 0.d0) .or. &
             (turnover_time > 0.d0 .and. rate_decomposition >= 0.d0) .or. &
             (rate_decomposition >= 0.d0 .and. rate_constant >= 0.d0) ) then
          option%io_buffer = 'Only TURNOVER_TIME or RATE_CONSTANT or RATE_DECOMPOSITION ' // &
            'may be included in a SomDec reaction definition, but not any two. ' // &
            'See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        else if (turnover_time > 0.d0) then
          new_reaction%rate_constant = 1.d0 / turnover_time
          new_reaction%rate_decomposition = -1.d0
        else if (rate_constant >= 0.d0) then
          new_reaction%rate_constant = rate_constant
          new_reaction%rate_decomposition = -1.d0
        else if (rate_decomposition >= 0.d0) then
          new_reaction%rate_constant = -1.d0
          new_reaction%rate_decomposition = rate_decomposition
        else
          option%io_buffer = 'ONE of TURNOVER_TIME or RATE_CONSTANT or RATE_DECOMPOSITION ' // &
            'must be correctly set in a SomDec reaction definition ' // &
            'See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        endif
        if(rate_ad_factor .ne. 1.d0) then
          new_reaction%rate_ad_factor = rate_ad_factor
        endif
        if (associated(this%reactions)) then
          prev_reaction%next => new_reaction
        else
          this%reactions => new_reaction
        endif
        prev_reaction => new_reaction
        nullify(new_reaction)        
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,SomDec keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine SomDecRead

! ************************************************************************** !

subroutine SomDecSetup(this,reaction,option)
  ! 
  ! Sets up SomDec reaction after it has been read from input
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
  ! 

  use Reaction_Aux_module
  use Option_module
  use String_module
  use Reaction_Immobile_Aux_module
  use Utility_module, only : DeallocateArray
  
  implicit none

  class(reaction_sandbox_somdec_type) :: this
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
  type(somdec_reaction_type), pointer :: cur_rxn
  
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
  allocate(this%rate_decomposition(this%nrxn))
  allocate(this%rate_ad_factor(this%nrxn))
  allocate(this%upstream_is_varycn(this%nrxn))
  allocate(this%upstream_c_id(this%nrxn))
  allocate(this%upstream_n_id(this%nrxn))
  allocate(this%upstream_nc(this%nrxn))
  allocate(this%upstream_is_aqueous(this%nrxn))
  
  allocate(this%downstream_id(this%nrxn,max_downstream_pools))
  allocate(this%downstream_stoich(this%nrxn,max_downstream_pools))
  allocate(this%downstream_nc(this%nrxn,max_downstream_pools))
  allocate(this%downstream_is_varycn(this%nrxn,max_downstream_pools))
  allocate(this%downstream_is_aqueous(this%nrxn,max_downstream_pools))
  allocate(this%mineral_c_stoich(this%nrxn))
  allocate(this%mineral_n_stoich(this%nrxn))

  this%pool_nc_ratio = 0.d0
  this%rate_constant = 0.d0
  this%rate_decomposition = 0.d0
  this%rate_ad_factor = 1.d0
  this%upstream_is_varycn = PETSC_FALSE
  this%upstream_c_id = 0
  this%upstream_n_id = 0
  this%upstream_nc = UNINITIALIZED_DOUBLE
  this%upstream_is_aqueous = PETSC_FALSE

  this%downstream_id = 0
  this%downstream_is_aqueous = PETSC_FALSE
  this%downstream_stoich = 0.d0
  this%downstream_is_varycn = PETSC_FALSE
  this%mineral_c_stoich = 0.d0
  this%mineral_n_stoich = 0.d0
  
! temporary array for mapping pools in reactions
  allocate(pool_names(this%npool))
  allocate(pool_is_aqueous(this%npool))
  allocate(species_id_pool_c(this%npool))
  allocate(species_id_pool_n(this%npool))

  pool_names = ''
  pool_is_aqueous = PETSC_FALSE
  species_id_pool_c = UNINITIALIZED_INTEGER
  species_id_pool_n = UNINITIALIZED_INTEGER

! pools
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    this%pool_nc_ratio(icount) = cur_pool%nc_ratio
    pool_names(icount) = cur_pool%name

    if (cur_pool%nc_ratio <= 0.d0) then
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
        option%io_buffer = 'For SomDec pools with no CN ratio defined, ' // &
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
             option%io_buffer = 'SomDec pool: ' // cur_pool%name // &
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
         this%upstream_is_varycn(icount) = PETSC_TRUE
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
            this%downstream_is_varycn = PETSC_TRUE
            ! (TODO) may need modify this to have more general applications
            option%io_buffer = 'For SomDec reactions, downstream pool ' // &
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
    this%rate_decomposition(icount) = cur_rxn%rate_decomposition
    this%rate_ad_factor(icount) = cur_rxn%rate_ad_factor

    cur_rxn => cur_rxn%next
  enddo 
  
  deallocate(pool_names)
  call DeallocateArray(pool_is_aqueous)
  call DeallocateArray(species_id_pool_c)
  call DeallocateArray(species_id_pool_n)

! set stoichiometric coefficients for fixed-CNratio SOM decomposition reactions
  do icount = 1, this%nrxn
     if(this%upstream_is_varycn(icount)) then
       cycle
     else
       ! calculate respiration C fraction from upstream C decomposition
       stoich_c = 1.0d0
       stoich_n = this%upstream_nc(icount)

       do jcount = 1, this%n_downstream_pools(icount)
          stoich_c = stoich_c - this%downstream_stoich(icount, jcount)
          stoich_n = stoich_n - this%downstream_stoich(icount, jcount) * &
                              this%downstream_nc(icount, jcount)
       enddo

       if(stoich_c < 0.0d0) then
         option%io_buffer = 'SomDec fixedCN C pool decomposition '// &
          'has a negative respiration fraction!' // &
            'Please check reactions of upstream pool: ' // &
             trim(cur_pool%name)
         call printErrMsg(option)
       endif

       if(stoich_n < 0.0d0) then
         option%io_buffer = 'SomDec fixedCN C pool N mineralization is negative,' // &
           'which (i.e. N immobilization) is not supported currently.' // &
           ' VariableCN ratio suggested.'
         call printErrMsg(option)
       endif

       this%mineral_c_stoich(icount) = stoich_c  ! this will be updated, if variable C/N Organic C pool(s) exists.
       this%mineral_n_stoich(icount) = stoich_n  ! this will be updated, if variable C/N Organic C pool(s) exists.

     endif
  enddo

  word = 'CO2(aq)'
  this%species_id_co2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  if(this%species_id_co2 <= 0) then
    option%io_buffer = 'CO2(aq) is not specified in the input file!'
    call printErrMsg(option)
  endif

  word = 'NH4+'
  this%species_id_nh4 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  if(this%species_id_nh4 <= 0) then
    option%io_buffer = 'NH4+ is not specified in the input file!'
    call printErrMsg(option)
  endif

  word = 'NO3-'
  this%species_id_no3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  
  word = 'N2O(aq)'
  this%species_id_n2o = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'H+'
  this%species_id_proton = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'HRimm'
  this%species_id_hrimm = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
#ifdef CLM_PFLOTRAN
  if(this%species_id_hrimm < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMdec: ' // &
       ' HRimm is not specified as immobile species in the input file. ' // &
       ' It is required when coupled with CLM.'
     call printErrMsg(option)
  endif
#endif

  word = 'Nmin'
  this%species_id_nmin = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
#ifdef CLM_PFLOTRAN
  if(this%species_id_nmin < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMdec: ' // &
       ' Nmin is not specified as immobile species in the input file. ' // &
       ' It is required when coupled with CLM.'
     call printErrMsg(option)
  endif
#endif
 
  word = 'Nimm'
  this%species_id_nimm = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
#ifdef CLM_PFLOTRAN
  if(this%species_id_nimm < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,SomDec: ' // &
       ' Nimm is not specified as immobile species in the input file. ' // &
       ' It is required when coupled with CLM.'
     call printErrMsg(option)
  endif
#endif
 
  word = 'NGASmin'
  this%species_id_ngasmin = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
#ifdef CLM_PFLOTRAN
  if(this%species_id_ngasmin < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,SomDec: ' // &
       ' NGASmin is not specified as immobile species in the input file. ' // &
       ' It is required when coupled with CLM.'
     call printErrMsg(option)
  endif
#endif

end subroutine SomDecSetup

! ************************************************************************** !
subroutine SomDecReact(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                       global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
  !
  ! Rewritten by Fengming Yuan @Aug-14-2014. The orginal was totally messed up,
  ! which caused a lot of issues.
  ! 
! ----------------------------------------------------------------------------!
  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only : material_auxvar_type

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif

  implicit none

#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif

  class(reaction_sandbox_somdec_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscReal :: saturation
  PetscReal :: theta
  PetscReal :: psi
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr
  PetscInt, parameter :: iphase = 1
  PetscScalar, pointer :: zsoi_pf_loc(:)    !

  PetscReal :: temp_real
  PetscReal :: c_uc, c_un    ! concentration (mole/m3 or mole/L, upon species type) => mole/m3bulk if needed
  PetscReal :: c_nh4         ! concentration (mole/L): substrate OR product,  => mole/m3bulk
  PetscReal :: c_no3         ! concentration (mole/L): substrate only,  => mole/m3bulk

  PetscInt :: irxn
  PetscInt :: ipool_up, ipool_down

  PetscInt :: ispec_uc, ispec_un, ispec_d   ! species id for upstream C, N, and downstream
  PetscInt :: ires_uc, ires_un, ires_d      ! id used for residual and Jacobian
  PetscInt :: ires_co2, ires_nh4, ires_n2o, ires_no3
  PetscInt :: ires_hrimm, ires_nmin, ires_nimm, ires_ngasmin
  PetscReal :: stoich_c, stoich_n

  PetscReal :: scaled_crate_const
  PetscReal :: tc     ! temperature in degC
  PetscReal :: f_t    ! temperature response function
  PetscReal :: f_w    ! moisture response function
  PetscReal :: f_depth     ! reduction factor due to deep soil

  PetscReal :: crate       ! moleC/s: overall = crate_uc, or,
                           !                  = crate_uc*crate_nh4, or,
                           !                  = crate_uc*(crate_nh4*(1-fnh4_inhibit_no3)+crate_no3*fnh4_inhibit_no3)
  PetscReal :: dcrate_dx   ! d(crate)/d(x), x = uc, un, dc, dn, nh4, or no3, respectively

  PetscReal :: crate_uc        ! crate function of upstream C ('uc')
  PetscReal :: dcrate_uc_duc   ! d crate_uc / d uc

  PetscReal :: crate_nh4        ! crate function of c_nh4 (for nh4 immobilization)
  PetscReal :: dcrate_nh4_dnh4  ! d(crate_nh4)/d(nh4) (for nh4 immobilization)
  PetscReal :: fnh4             ! = nh4/(half_saturation + nh4)  ( N 'resource' limitation on immobilization)
  PetscReal :: dfnh4_dnh4       ! d(fnh4)/d(nh4)

  PetscReal :: crate_no3        ! crate function of c_no3 (for no3 immobilization)
  PetscReal :: dcrate_no3_dno3  ! d(crate_no3)/d(no3)
  PetscReal :: fno3             ! = no3/(half_saturation + no3)  ( N 'resource' limitation on immobilization)
  PetscReal :: dfno3_dno3       ! d(fnh4)/d(nh4)

  PetscReal :: nratecap, dtmin  ! max. n immobilization rate within allowable min. timestep
  PetscReal :: fnratecap        ! max. nratecap as function of c_nh4/c_no3 consumption rate vs. potential immobilization
  PetscReal :: dfnratecap_dnh4, dfnratecap_dno3

  !nh4 inhibition on no3 immobilization, or microbial N immobilization preference btw nh4 and no3
  ! crate_nh4 = fnh4*fnh4_inhibit_no3, while crate_no3 = 1.-fnh4*fnh4_inhibition_no3
  ! by the following eq., if 'inhibition=1', uptake will be equal for both (if same conc.);
  ! higher coef, strong NH4 inhibition on NO3 (i.e., more NH4 immobilization over NO3)
  PetscReal :: fnh4_inhibit_no3 ! inhibition_coef/(inhibition_coef + no3/nh4):
  PetscReal :: dfnh4_inhibit_no3_dnh4 ! d(fnh4_inhibit_no3)/dnh4
  PetscReal :: dfnh4_inhibit_no3_dno3 ! d(fnh4_inhibit_no3)/dno3

  ! save mineral N fraction and decomposition rate for net N mineralization and associated N2O calculation
  PetscReal :: net_nmin_rate
  PetscReal :: dnet_nmin_rate_dx
  PetscReal :: ph, f_ph
  PetscReal :: rate_n2o, drate_n2o_dx

! other local variables
  ! UC + u*UN --> (1-di)*DCi +(u-ni-n2)*DNi+ d*CO2 + n*N[H3] [+ n2*NO3], and
  ! the reaction (decomposition) rate = crate, and independent unkowns: d, n, and possibly u, n2. ('di','ni' imply multi-downpools)
  ! n1 may be negative, i.e., n1*nh4 should be in the left side of the Eq. (N immobilization).
  ! n2 must be negative, i.e., n2*NO3 should be in the left side of the Eq. (N immobilization).
  ! u = (UC/UN), uc/un for upstream C/N (UC/UN here), n for nh4, c for CO2
  ! So, total 7 variables, and 6 possible independent variables.
  ! but, (1) 'u-ni-n2' is implicitly completely dependent on others currently, so ignored in the code. So 6 variables;
  !      (2) any rates NOW are not functions of down-stream C-N pools, So 4 possible independent variables.
  PetscReal :: dc_duc, duc_duc, ddc_duc, dn_duc, dno3_duc, du_duc        ! 'uc' as independent variable
  PetscReal :: dc_dun, duc_dun, ddc_dun, dn_dun, dno3_dun, du_dun        ! 'un' as independent variable (NOT YET either - TODO)

  PetscReal :: dc_dn, duc_dn, ddc_dn, dn_dn, dno3_dn, du_dn              ! 'nh4' as independent variable
  PetscReal :: dc_dno3, duc_dno3, ddc_dno3, dn_dno3, dno3_dno3, du_dno3  ! 'no3' as independent variable

! misc. local variables
  PetscInt :: i, j, ires
  PetscReal:: feps0, dfeps0_dx
  PetscReal:: k_decomp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !------------------------------------------------------------------------------------------------
  porosity = material_auxvar%porosity
  volume = material_auxvar%volume
  saturation = global_auxvar%sat(1)

  theta = saturation * porosity
  psi = min(global_auxvar%pres(1) - option%reference_pressure, -1.d-20)     ! if positive, saturated soil's psi is nearly zero

  ! moisture response function
#ifdef CLM_PFLOTRAN
  ghosted_id = option%iflag
  if(this%moisture_response_function == MOISTURE_RESPONSE_FUNCTION_CLM4) then
     f_w = GetMoistureResponse(theta, ghosted_id, this%moisture_response_function)
  elseif(this%moisture_response_function == MOISTURE_RESPONSE_FUNCTION_DLEM) then
     f_w = GetMoistureResponse(theta, ghosted_id, this%moisture_response_function)
  endif
#else
  f_w = 1.0d0
#endif

  !----------------------------------------------------------------------------------------
  ! temperature response function
  tc = global_auxvar%temp

#ifdef CLM_PFLOTRAN
  f_t = GetTemperatureResponse(tc,this%temperature_response_function, this%Q10)
#else
  f_t = 1.0d0
#endif

  if(f_t < 1.0d-20 .or. f_w < 1.0d-20) then
     return
  endif

#ifdef CLM_PFLOTRAN
  if (this%decomp_depth_efolding > 0.d0) then
    call VecGetArrayReadF90(clm_pf_idata%zsoi_pfs, zsoi_pf_loc, ierr)
    CHKERRQ(ierr)
    f_depth = exp(-zsoi_pf_loc(ghosted_id)/this%decomp_depth_efolding)
    f_depth = min(1.0d0, max(1.d-20, f_depth))
    call VecRestoreArrayReadF90(clm_pf_idata%zsoi_pfs, zsoi_pf_loc, ierr)
    CHKERRQ(ierr)
  else
    f_depth = 1.0d0
  endif
#else
  f_depth = 1.0d0
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! -----------------------------------------------------------------------------------------------
  ! for 'nh4' dependent rate function (i.e. immobilization, will be dynamically changing for 'irxn')
  ires_nh4  = this%species_id_nh4 + reaction%offset_aqueous       ! as aq. species
  c_nh4     = rt_auxvar%total(ires_nh4, iphase)*theta*1000.0d0    ! from mol/L -> mol/m3 bulk

  fnh4_inhibit_no3 = 1.d0
  dfnh4_inhibit_no3_dnh4 = 0.d0
  dfnh4_inhibit_no3_dno3 = 0.d0
  if (this%species_id_no3 > 0) then
    ! for 'no3' dependent rate function (i.e. immobilization), if needed
    ires_no3  = this%species_id_no3 + reaction%offset_aqueous        ! as aq. species
    c_no3 = rt_auxvar%total(this%species_id_no3, iphase)*theta*1000.0d0    ! from mol/L -> mol/m3 bulk

    ! nh4 inhibition on no3 immobilization, if any ('this%inhibition_nh4_no3')
    ! this is for quantifying microbial N immobilization preference btw NH4 and NO3.
    ! (DON'T change the 'rate' and 'derivatives' after this)
    if(c_nh4>this%x0eps .and. c_no3>this%x0eps .and. this%inhibition_nh4_no3>0.d0) then
      temp_real = this%inhibition_nh4_no3 + c_no3/c_nh4
      fnh4_inhibit_no3 = this%inhibition_nh4_no3/temp_real
      dfnh4_inhibit_no3_dnh4 = this%inhibition_nh4_no3*(c_no3/c_nh4/c_nh4)  &
                               /temp_real/temp_real     ! over 'd_nh4'
      dfnh4_inhibit_no3_dno3 = -this%inhibition_nh4_no3/c_nh4  &
                               /temp_real/temp_real     ! over 'd_no3'
    else
      if (c_nh4>this%x0eps .and. c_no3<=this%x0eps) then
        fnh4_inhibit_no3 = 0.999d0
      elseif (c_nh4<=this%x0eps .and. c_no3>this%x0eps) then
        fnh4_inhibit_no3 = 0.001d0
      else
        fnh4_inhibit_no3 = 0.50d0
      endif
      dfnh4_inhibit_no3_dnh4 = 0.d0
      dfnh4_inhibit_no3_dno3 = 0.d0
    endif

  endif ! if (this%species_id_no3>0)

  !----------------------------------------------------------------------------------------------
  ires_co2 = this%species_id_co2 + reaction%offset_aqueous       ! as aq. species
  if(this%species_id_hrimm > 0) then
     ires_hrimm = this%species_id_hrimm + reaction%offset_immobile
  endif

  if(this%species_id_nmin > 0) then
     ires_nmin = this%species_id_nmin + reaction%offset_immobile
  endif

  if(this%species_id_nimm > 0) then
     ires_nimm = this%species_id_nimm + reaction%offset_immobile
  endif

  net_nmin_rate     = 0.0d0
  dnet_nmin_rate_dx = 0.0d0

  do irxn = 1, this%nrxn
  
    !-----------------------------------------------------------------------------------------------------
    ! (i) calculate C rate ( and, N rate is derived by Crate*N/C)

    ! rate in unit of 1/second
    if (this%rate_constant(irxn) >= 0.d0) then
      k_decomp = this%rate_constant(irxn)
    elseif (this%rate_decomposition(irxn) >= 0.d0) then
      k_decomp = 1.0d0-exp(-this%rate_decomposition(irxn)*option%tran_dt)
      k_decomp = k_decomp/option%tran_dt
    endif
    k_decomp = this%rate_ad_factor(irxn)*k_decomp
    k_decomp = min(k_decomp, 1.0d0/volume/option%tran_dt)        ! make sure of NO over-decomposition rate (just in case, maybe not needed)

    ! scaled_rate_const units: (m^3 bulk / s) = (1/s) * (m^3 bulk)
    scaled_crate_const = k_decomp*volume*f_t*f_w*f_depth

    ! (ii) substrates
    ispec_uc = this%upstream_c_id(irxn)
    if(this%upstream_is_aqueous(irxn)) then
      c_uc = rt_auxvar%total(ispec_uc, iphase)
      c_uc = theta * 1000.0d0 * c_uc    ! from mol/L -> mol/m3
      ires_uc = reaction%offset_aqueous + ispec_uc
    else
      c_uc = rt_auxvar%immobile(ispec_uc)
      ires_uc = reaction%offset_immobile + ispec_uc
    endif

    if(this%x0eps>0.d0) then
      ! GP's cut-off approach (sort of Heaviside function)
      feps0     = funcTrailersmooth(c_uc, this%x0eps*10.d0, this%x0eps, PETSC_FALSE)
      dfeps0_dx = funcTrailersmooth(c_uc, this%x0eps*10.d0, this%x0eps, PETSC_TRUE)

    else
      feps0 = 1.0d0
      dfeps0_dx = 0.d0
      if(c_uc <= this%x0eps) cycle     ! this may bring in 'oscillation' around 'this%x0eps'
    endif

    ! (iii) C substrate only dependent rate/derivative  (DON'T change after this)
    crate_uc  = scaled_crate_const * c_uc * feps0    ! moles/s: (m3 bulk/s)* (moles/m3 bulk) * (-)
    dcrate_uc_duc = scaled_crate_const * (feps0 + c_uc*dfeps0_dx)

    !-----------------------------------------------------------------------------------------------------

    ! for varying N/C ratio C pool decomposition reactions, N/C stoichiometry needs to be calculated on the fly
    if(this%upstream_is_varycn(irxn)) then

      ispec_un = this%upstream_n_id(irxn)
      if (ispec_un>0) then  ! just in case
        if(this%upstream_is_aqueous(irxn)) then
          c_un = rt_auxvar%total(ispec_un, iphase)
          c_un = theta * 1000.0d0 * c_uc    ! from mol/L -> mol/m3
          ires_un = ispec_un + reaction%offset_aqueous
        else
          c_un = rt_auxvar%immobile(ispec_un)
          ires_un = ispec_un + reaction%offset_immobile
        endif
        this%upstream_nc(irxn) = c_un / c_uc

        ! calculate respiration factor (CO2 stoichiometry)
        stoich_c = 1.0d0
        do j = 1, this%n_downstream_pools(irxn)
          stoich_c = stoich_c - this%downstream_stoich(irxn, j)
        enddo

        if(stoich_c < 0.0d0) then
          option%io_buffer = 'SomDec variableCN ratio C pool decomposition has' // &
                              'a negative respiration fraction!'
          call printErrMsg(option)
        endif
        this%mineral_c_stoich(irxn) = stoich_c

        ! calculate N (nh4) stoichiometry
        stoich_n = this%upstream_nc(irxn)
        do j = 1, this%n_downstream_pools(irxn)
          stoich_n = stoich_n - this%downstream_stoich(irxn, j) * &
                              this%downstream_nc(irxn, j)
        enddo
        this%mineral_n_stoich(irxn) = stoich_n

      endif

    endif ! end of 'upstream_is_varycn' (NOTE: 'mineral_n_stoich' for non-litter species IS constant, which are calc. in 'setup')

    !-----------------------------------------------------------------------------------------------------
    ! calculation of residuals
    ! residual units: (mol/sec) = (m^3 bulk/s) * (mol/m^3 bulk)

    ! overall C rates
    crate     = crate_uc

    ! NH4/NO3 limiting on N-immobolization involved C rates, if any
    ! NOTE: don't change the 'fnh4' ('fno3') and 'dfnh4_dnh4'('dfno3_dno3') after this block of codes
    ! it will be used for each possible 'irxn' species (and for this reason, initializing it here)
    fnh4 = 1.d0
    dfnh4_dnh4 = 0.d0
    fno3 = 1.d0
    dfno3_dno3 = 0.d0
    crate_nh4  = 1.d0
    crate_no3  = 1.d0

    if(this%mineral_n_stoich(irxn) < 0.0d0) then

      ! half-saturation regulated rate (by input, it may be representing competetiveness)
      fnh4       = funcMonod(c_nh4, this%half_saturation_nh4, PETSC_FALSE)
      dfnh4_dnh4 = funcMonod(c_nh4, this%half_saturation_nh4, PETSC_TRUE)
      if(this%species_id_no3 > 0) then
        fno3       = funcMonod(c_no3, this%half_saturation_no3, PETSC_FALSE)
        dfno3_dno3 = funcMonod(c_no3, this%half_saturation_no3, PETSC_TRUE)
      endif

      ! using the following for trailer smoothing
      ! note: 'x0eps' is different from 'half_saturation_nh4' above.
      !       'x0eps' is for mathematic reason in the code;
      !       'half_saturation_nh4' is for using 'monod-type' function to quantify microbial
      !    NH4 immobilization dependence on resources (NH4). So, physiologically it may be
      !    used as a method to quantify competition over resources.
      if(this%x0eps>0.d0) then
        ! GP's cut-off approach (sort of Heaviside function)
        feps0     = funcTrailersmooth(c_nh4, this%x0eps*10.d0, this%x0eps, PETSC_FALSE)
        dfeps0_dx = funcTrailersmooth(c_nh4, this%x0eps*10.d0, this%x0eps, PETSC_TRUE)
      else
        feps0 = 1.0d0
        dfeps0_dx = 0.d0
      endif
      dfnh4_dnh4 = dfnh4_dnh4*feps0 + fnh4 * dfeps0_dx
      fnh4 = fnh4 * feps0
      !
      if(this%species_id_no3 > 0) then
        if(this%x0eps>0.d0) then
          ! GP's cut-off approach
          feps0     = funcTrailersmooth(c_no3, this%x0eps*10.d0, this%x0eps, PETSC_FALSE)
          dfeps0_dx = funcTrailersmooth(c_no3, this%x0eps*10.d0, this%x0eps, PETSC_TRUE)
        else
          feps0 = 1.0d0
          dfeps0_dx = 0.d0
        endif
        dfno3_dno3 = dfno3_dno3*feps0 + fno3 * dfeps0_dx
        fno3 = fno3 * feps0
      endif

      ! constraining 'N immobilization rate' locally if too high compared to available within the allowable min. time-step
      ! It can be achieved by cutting time-step, but it may be taking a very small timestep finally
      ! - implying tiny timestep in model, which potentially crashes model

      ! in the following, '2.1' multiplier is chosen because that is slightly larger(5% to avoid numerical issue) than '2.0',
      ! which should be the previous-timestep before reaching the 'option%dt_min'
      ! (the time-cut used in PF is like dt=0.5*dt, when cutting)
      dtmin = 2.1d0*option%dt_min
      !dtmin = max(option%tran_dt, 2.1d0*option%dt_min)   ! this 'dtmin' may be accelerating the timing, but may not be appropriate to mulitple consummers

      if (crate_uc*this%mineral_n_stoich(irxn) < 0.d0) then
        !
        nratecap = -crate_uc*this%mineral_n_stoich(irxn)*dtmin
        if (nratecap*fnh4_inhibit_no3 > c_nh4) then
          ! assuming that monod type reduction used and 'nratecap' must be reduced to 'c_nh4', i.e.
          ! c_nh4/dt = nratecap*fnh4_inhibit_no3/dt * (c_nh4/(c_nh4+c0))
          ! then, c0 = nratecap*fnh4_inhibit_no3 - c_nh4 (this is the 'half-saturation' term used above)
          fnratecap = funcMonod(c_nh4, nratecap*fnh4_inhibit_no3-c_nh4, PETSC_FALSE)
          dfnratecap_dnh4 = funcMonod(c_nh4, nratecap*fnh4_inhibit_no3-c_nh4, PETSC_TRUE)
        else
          fnratecap       = 1.d0
          dfnratecap_dnh4 = 0.d0
        endif
        dfnh4_dnh4 = dfnh4_dnh4*fnratecap + fnh4 * dfnratecap_dnh4
        fnh4 = fnh4 * fnratecap
        !
        if (this%species_id_no3 > 0) then
          if (nratecap*(1.d0-fnh4_inhibit_no3) > c_no3) then
            fnratecap = funcMonod(c_no3, nratecap*(1.d0-fnh4_inhibit_no3)-c_no3, PETSC_FALSE)
            dfnratecap_dno3 = funcMonod(c_no3, nratecap*(1.d0-fnh4_inhibit_no3)-c_no3, PETSC_TRUE)
          else
            fnratecap       = 1.d0
            dfnratecap_dno3 = 0.d0
          endif
          dfno3_dno3 = dfno3_dno3*fnratecap + fno3 * dfnratecap_dno3
          fno3 = fno3 * fnratecap
        endif
        !

      endif

      !----------- overall rate ------
      crate_nh4 = fnh4
      ! overall C rates
      crate = crate_uc * crate_nh4

      if (this%species_id_no3 > 0) then
        ! 'f_nh4_inhibit_no3' is somehow a spliting fraction for NH4:NO3 immoblization,
        crate_nh4 = fnh4 * fnh4_inhibit_no3

        ! f_no3 should be adjusted by R reduced by 'crate_nh4' so that it is inhibited by nh4
        crate_no3 = fno3 * (1.0d0-fnh4_inhibit_no3)

        ! overall C rates
        crate = crate_uc * (crate_nh4 + crate_no3)
      endif

!#ifdef DEBUG
      if(option%tran_dt<=option%dt_min .and. option%print_file_flag) then
        if (-crate_uc*crate_nh4*this%mineral_n_stoich(irxn)*option%dt_min>=c_nh4-this%x0eps .or.  &
            -crate_uc*crate_no3*this%mineral_n_stoich(irxn)*option%dt_min>=c_no3-this%x0eps) then
          write(option%fid_out, *) '----------------------------------------------'
          write(option%fid_out, *) 'Reaction Sandbox: SOMDEC'
          write(option%fid_out, *) 'dt=',option%tran_dt, ' dt_min=',option%dt_min
          write(option%fid_out, *) 'ghosted_id=',ghosted_id, ' irxn=', irxn, &
            ' c_nh4=',c_nh4, ' c_no3=',c_no3, &
            ' fnh4=',fnh4,'crate_nh4=',crate_uc*crate_nh4*this%mineral_n_stoich(irxn)*option%dt, &
            'crate0_nh4=',crate_uc*this%mineral_n_stoich(irxn)*fnh4_inhibit_no3*option%dt, &
            ' fno3=',fno3,'crate_no3=',crate_uc*crate_no3*this%mineral_n_stoich(irxn)*option%dt, &
             'crate0_n03=',crate_uc*this%mineral_n_stoich(irxn)*(1.d0-fnh4_inhibit_no3)*option%dt
        endif
      endif
!#endif

    endif  !if(this%mineral_n_stoich(irxn) < 0.0d0) (i.e. immobilization)

    !------------------------------------------------------------------------------------------
    ! -- Residuals for all C-N pools
    ! CO2 [1]
    Residual(ires_co2) = Residual(ires_co2) - this%mineral_c_stoich(irxn) * crate
    if (this%species_id_hrimm > 0) then    ! for tracking
       Residual(ires_hrimm) = Residual(ires_hrimm) - this%mineral_c_stoich(irxn) * crate
    endif
    
    ! upstream c [2]
    Residual(ires_uc) = Residual(ires_uc) + crate

    ! downstream non-mineral pools (fixed C/N ratio only now - TODO: will modify) [3]
    do j = 1, this%n_downstream_pools(irxn)
       ispec_d = this%downstream_id(irxn, j)
       if(this%downstream_is_aqueous(irxn, j)) then
         ires_d = ispec_d
       else
         ires_d = reaction%offset_immobile + ispec_d
       endif
       if(ispec_d > 0) then
          Residual(ires_d) = Residual(ires_d) - this%downstream_stoich(irxn, j) * crate
       endif
    enddo

    ! inorg. nitrogen (NH4+/nh4(aq)) [4], [and, no3 [5] ]
    if(this%mineral_n_stoich(irxn) >= 0.0d0) then        ! mineralization
      Residual(ires_nh4) = Residual(ires_nh4) - this%mineral_n_stoich(irxn) * crate

      if(this%species_id_nmin > 0) then    ! for tracking
        Residual(ires_nmin) = Residual(ires_nmin) - this%mineral_n_stoich(irxn) * crate
      endif

    else                                                 ! immobilization

      Residual(ires_nh4) = Residual(ires_nh4) - &
        this%mineral_n_stoich(irxn) * crate_uc * crate_nh4
      if (this%species_id_no3 > 0) then
        Residual(ires_no3) = Residual(ires_no3) - &
          this%mineral_n_stoich(irxn) * crate_uc * crate_no3
      endif

      if(this%species_id_nimm > 0) then    ! for tracking
         Residual(ires_nimm) = Residual(ires_nimm) + this%mineral_n_stoich(irxn) * crate
      endif

    endif

    net_nmin_rate = net_nmin_rate + &
        this%mineral_n_stoich(irxn) * crate

    ! upstream n [6]
    if(this%upstream_is_varycn(irxn)) then
      Residual(ires_un) = Residual(ires_un) + this%upstream_nc(irxn) * crate
    endif
    

    !-----------------------------------------------------------------------------------------------------
    ! calculate jacobians
    if (compute_derivative) then

      ! ----- SOM decomposition network --------------------------------------------------------

! -- derivatives with fixed C/N decomposing pools AND no reduction of nh4 (positive 'n')

        ! SOMc + u SOMn -> (1-di) SOMci + (u-n) SOMni + d C[O2] + n N[H3]
        ! [2]    [6]       [3]            [ignored]     [1]       [4]
        ! d - minearl_c_stoich, n - mineral_n_stoich, u - upstream_nc, 1-di - (multi-)downstream_stoich
        ! NOTE: (1) '(u-n)' is implicitly calc'ed (totally ignored in the codes, so NO 'ddn_d??' in the codes)
        !       (2) 'SOMci' will be caclucated below, because multiple-downstream pools('ddc_d??') are allowed
        !       (3) 'u' here is implicitly calc'ed, but not for variale-CN pools
        !       (4) currently, rates NOT depends upon 'un' or 'u', so,
        dc_dun  = 0.d0     ! [6-1] 'C rates'
        duc_dun = 0.d0     ! [6-2]
        ddc_dun = 0.d0     ! [6-3 - will be done in Jacobians below)

        dn_dun  = 0.d0     ! [6-4] 'N rates'
        dno3_dun= 0.d0     ! [6-5 - not associated with SOM N mineralization]
        du_dun  = 0.d0     ! [6-6]

        ! so, 'uc' is the only independent variable, by this point
        dcrate_dx = dcrate_uc_duc
        dc_duc  = dcrate_dx * this%mineral_c_stoich(irxn)              ! [6-1]
        duc_duc = -1.0d0*dcrate_dx                                     ! [6-2]
        ddc_duc = 0.d0                                                 ! [6-3 - will be done in Jacobians below)

        dn_duc  = dc_duc * this%mineral_n_stoich(irxn)                 ! [6-4]
        dno3_duc= 0.d0                                                 ! [6-5 - not associated with SOM decomposition]
        du_duc  = this%upstream_nc(irxn) * duc_duc                     ! [6-6]

! -- derivatives with variable C/N decompsing pools
        ! LitC + u LitN -> (1-di) SOMci + (u-n-n2) SOMni + d C[O2] + n N[H3] [+ n2 NO3]
        ! [2]    [6]       [3]            [ignored]          [1]       [4]        [5]
        ! 'n+n2' IS mineral_n_stoich
        ! Note: (1) only 'u-n-n2' is implicitly calc'ed (i.e., ignored in the codes)

        if(this%upstream_is_varycn(irxn)) then

! -- derivatives with variable C/N decompsing pools and N[H3] as reactants
           ! 'un' is limited, so 'n' is negative (nh4 as reactant, i.e. N immobilization)
           if (this%mineral_n_stoich(irxn) < 0.d0) then

             ! -- depending on 'uc'
             ! crate = crate_uc * crate_nh4 = crate_uc*fnh4,
             ! So, d(crate)/d(uc) = fnh4*dcrate_uc_duc + crate_uc * dfnh4_duc
             !     in which 'dfnh4_duc' = 0, because crate_nh4 = fnh4, only dependent on 'nh4'
             dcrate_dx = dcrate_uc_duc*fnh4
             dc_duc  =  dcrate_dx * this%mineral_c_stoich(irxn)                         ! [6-1]
             duc_duc = -1.0d0 * dcrate_dx                                               ! [6-2]

             dn_duc  = dc_duc * this%mineral_n_stoich(irxn)                             ! [6-4]
             du_duc  = this%upstream_nc(irxn) * duc_duc                                 ! [6-6]

             ! -- depending on 'nh4'
             ! d(crate)/d(n) = d(crate_nh4)*crate_uc = dfnh4_dnh4*crate_uc
             dcrate_dx = dfnh4_dnh4*crate_uc
             dc_dn  = dcrate_dx * this%mineral_c_stoich(irxn)                           ! [6-1]
             duc_dn = -1.0d0*dcrate_dx                                                  ! [6-2]

             dn_dn  = dc_dn * this%mineral_n_stoich(irxn)                               ! [6-4]
             du_dn  = this%upstream_nc(irxn) * duc_dn                                   ! [6-6]

! -- derivatives with variable C/N decompsing pools and N[H3]+NO3 as reactants
             if (this%species_id_no3 > 0) then
               ! -- depending on 'uc'
               ! crate = crate_uc * (crate_nh4 + crate_no3)
               !       = crate_uc * (fnh4*fnh4_inhibit_no3+fno3*(1-fnh4_inhibit_no3)),
               ! So, d(crate)/d(uc) = (fnh4*fnh4_inhibit_no3+fno3-fno3*fnh4_inhibit_no3)*dcrate_uc_duc
               dcrate_dx  = dcrate_uc_duc * &
                         (fnh4*fnh4_inhibit_no3+fno3-fno3*fnh4_inhibit_no3)

               dc_duc  = dcrate_dx * this%mineral_c_stoich(irxn)                         ! [6-1]

               duc_duc = -1.0d0 * dcrate_dx                                              ! [6-2]

               ! nrate = crate_uc*mineral_n_stoich*fnh4*fnh4_inhibit_no3: 'mineral_n_stoich = u-sum(ni)'
               dn_duc  = dcrate_uc_duc* this%mineral_n_stoich(irxn) * &
                         fnh4*fnh4_inhibit_no3                                           ! [6-4]

               ! no3rate = crate_uc*mineral_n_stoich*fno3*(1-fnh4_inhibit_no3): 'mineral_n_stoich = u-sum(ni)'
               dno3_duc= dcrate_uc_duc* this%mineral_n_stoich(irxn) * &
                         fno3*(1.0d0-fnh4_inhibit_no3)                                   ! [6-5]

               du_duc  = this%upstream_nc(irxn) * duc_duc                                ! [6-6]

               ! -- depending on 'nh4'
               ! crate = crate_uc * (crate_nh4 + crate_no3)
               ! d(crate)/d(n) = d(crate_nh4+crate_no3)*crate_uc/dnh4
               !               = d(fnh4*fnh4_inhibit_no3+fno3*(1.-fnh4_inhibit_no3))/dnh4*crate_uc
               !               = d(fnh4_inhibit_no3*(fnh4-fno3) + fno3)/dnh4*crate_uc
               !               = d(fnh4_inhibit_no3*(fnh4-fno3))/dnh4*crate_uc                     ! dfno3_dnh4 = 0
               !               = (fnh4-fno3)*dfnh4_inhibit_no3_dnh4
               !                  +dfnh4_dnh4*fnh4_inhibit_no3)
               !                 *crate_uc
               dcrate_dx  = ( dfnh4_dnh4*fnh4_inhibit_no3  &
                             +(fnh4-fno3)*dfnh4_inhibit_no3_dnh4)
               dcrate_dx  = dcrate_dx*crate_uc

               dc_dn  = dcrate_dx * this%mineral_c_stoich(irxn)                          ! [6-1]

               duc_dn = -1.0d0 * dcrate_dx                                               ! [6-2]

               ! nrate = crate_uc*mineral_n_stoich*fnh4*fnh4_inhibit_no3: 'mineral_n_stoich = u-sum(ni)'
               ! d(nrate)/d(n) = ( fnh4*dfnh4_inhibit_no3_dnh4
               !                  +dfnh4_dnh4*fnh4_inhibit_no3
               !                 ) * crate_uc * mineral_n_stoich
               dn_dn  = fnh4*dfnh4_inhibit_no3_dnh4 &
                       +dfnh4_dnh4*fnh4_inhibit_no3
               dn_dn  = dn_dn * crate_uc*this%mineral_n_stoich(irxn)                     ! [6-4]

               ! no3rate = crate_uc*mineral_n_stoich*fno3*(1-fnh4_inhibit_no3): 'mineral_n_stoich = u-sum(ni)'
               ! d(no3rate)/d(n) = -1.0*(fno3*dfnh4_inhibit_no3_dnh4
               !                    +dfno3_dnh4*fnh4_inhibit_no3)         ! 'dfno3_dnh4' = 0
               !                  * crate_uc * mineral_n_stoich
               dno3_dn= -1.d0*fno3*dfnh4_inhibit_no3_dnh4
               dno3_dn= dno3_dn * crate_uc*this%mineral_n_stoich(irxn)                   ! [6-5]

               du_dn  = this%upstream_nc(irxn) * duc_dn                                  ! [6-6]

               ! -- depending on 'no3'
               ! d(crate)/d(no3) = d(crate_nh4+crate_no3)*crate_uc/dno3
               !               = d(fnh4*fnh4_inhibit_no3+fno3*(1-fnh4_inhibit_no3))/dno3*crate_uc
               !               = d((fnh4-fno3)*fnh4_inhibit_no3+fno3)/dno3*crate_uc
               !               = (-dfno3_dno3*fnh4_inhibit_no3           ! 'dfnh4_dno3' = 0
               !                  +(fnh4-fno3)*dfnh4_inhibit_no3_dno3
               !                  +dfno3_dno3
               !                 ) * crate_uc
               dcrate_dx = (fnh4-fno3)*dfnh4_inhibit_no3_dno3    &
                           + dfno3_dno3*(1.d0-fnh4*fnh4_inhibit_no3)
               dcrate_dx = dcrate_dx * crate_uc

               dc_dno3  = dcrate_dx * this%mineral_c_stoich(irxn)                        ! [6-1]

               duc_dno3 = -1.0d0*dcrate_dx                                               ! [6-2]

               ! nrate = crate_uc*mineral_n_stoich*fnh4*fnh4_inhibit_no3: 'mineral_n_stoich = u-sum(ni)'
               ! d(nrate)/d(no3) = ( fnh4*dfnh4_inhibit_no3_dno3
               !                    +dfnh4_dno3*fnh4_inhibit_no3               ! 'dfnh4_dno3' = 0
               !                   ) * crate_uc * mineral_n_stoich
               dn_dno3  = fnh4*dfnh4_inhibit_no3_dno3  &
                         * crate_uc * this%mineral_n_stoich(irxn)                        ! [6-4]

               ! no3rate = crate_uc*mineral_n_stoich*fno3*(1-fnh4_inhibit_no3): 'mineral_n_stoich = u-sum(ni)'
               ! d(no3rate)/d(no3) = ( fno3*(-dfnh4_inhibit_no3_dno3)     ! 'dfnh4_dno3' = 0
               !                    +dfno3_dno3*(1-fnh4_inhibit_no3)
               !                 ) * crate_uc * mineral_n_stoich
               dno3_dno3= -1.0d0*fno3*dfnh4_inhibit_no3_dno3 &
                          + dfno3_dno3 * (1.d0-fnh4_inhibit_no3)
               dno3_dno3= dno3_dno3 * crate_uc * this%mineral_n_stoich(irxn)             ! [6-5]

               du_dno3  = this%upstream_nc(irxn) * duc_dno3                              ! [6-6]

             endif  ! this%species_id_no3 > 0

           endif !this%mineral_n_stoich(irxn) < 0.d0

        endif  !this%upstream_is_varycn(irxn)

!--- mixed derivatives for 'net_nmin_rate' (TODO - need more thinking here?)
        dnet_nmin_rate_dx = dnet_nmin_rate_dx + dn_duc + dn_dun
        if(this%mineral_n_stoich(irxn) < 0.d0) then
          dnet_nmin_rate_dx = dnet_nmin_rate_dx + dn_dn
          if(this%species_id_no3>0) then
            dnet_nmin_rate_dx = dnet_nmin_rate_dx + dn_dno3
          endif
        endif

! -- Jacobians with respect to upstream C ('uc')
      ! CO2 [6-1]
      if(this%upstream_is_aqueous(irxn)) then   ! aqueous c pools must have fixed CN ratio
        Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - dc_duc*  &
          rt_auxvar%aqueous%dtotal(this%species_id_co2,ispec_uc,iphase)
      else
        Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - dc_duc
      endif
      !if(this%species_id_hrimm > 0) then    ! for tracking
      !  Jacobian(ires_hrimm,ires_uc) = Jacobian(ires_hrimm,ires_uc) - dc_duc
      !endif

      ! upstream C pool [6-2]
      if(this%upstream_is_aqueous(irxn)) then
        Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) - duc_duc * &
          rt_auxvar%aqueous%dtotal(ispec_uc,ispec_uc,iphase)
      else
        Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) - duc_duc
      endif

      ! downstream C pools [6-3]
      do j = 1, this%n_downstream_pools(irxn)
         ispec_d = this%downstream_id(irxn, j)
         if(ispec_d < 0) then
           option%io_buffer = 'Downstream pool species not specified!'
           call printErrMsg(option)
         endif

         if(this%downstream_is_aqueous(irxn, j)) then
           ires_d = ispec_d
         else
           ires_d = reaction%offset_immobile + ispec_d
         endif

         ! if given that: 'mineral_c_stoich + sum(di) = 1',
         !                 [6-1] dc_dx=dcrate_dx*mineral_c_stoich,
         !             and,[6-2] duc_dx=-dcrate_dx
         ! then, d(dci)/dx = -duc_dx * di
         ddc_duc = this%downstream_stoich(irxn, j) * (-1.d0*duc_duc)
         
         if(this%upstream_is_aqueous(irxn) .and. &
           this%downstream_is_aqueous(irxn, j)) then
             Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - ddc_duc * &
                 rt_auxvar%aqueous%dtotal(ispec_d,ispec_uc,iphase)
         else
             Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - ddc_duc
         endif
      enddo

      ! N[H3] [6-4]
      if(this%upstream_is_aqueous(irxn)) then
        Jacobian(ires_nh4,ires_uc) = Jacobian(ires_nh4,ires_uc) - dn_duc*  &
          rt_auxvar%aqueous%dtotal(this%species_id_nh4,ispec_uc,iphase)
      else
        Jacobian(ires_nh4,ires_uc) = Jacobian(ires_nh4,ires_uc) - dn_duc
      endif
      !!for tracking
      !if(this%species_id_nmin > 0 .and. this%mineral_n_stoich(irxn) > 0.0d0) then
      !  Jacobian(ires_nmin,ires_uc) = Jacobian(ires_nmin,ires_uc) - dn_duc
      !endif
      !if(this%species_id_nimm > 0 .and. this%mineral_n_stoich(irxn) < 0.0d0) then
      !  Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) + dn_duc
      !endif

      ! NO3 [6-5], if any
      if (this%species_id_no3>0 .and. &
          this%upstream_is_varycn(irxn) .and. &
          this%mineral_n_stoich(irxn) < 0.0d0) then

        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - &
            dno3_duc*  &
            rt_auxvar%aqueous%dtotal(this%species_id_no3,ispec_uc,iphase)
        else
          Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - &
            dno3_duc
        endif
        !!for tracking
        !if(this%species_id_nimm >0) then
        !  Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) + dn_duc
        !endif
      endif

      ! upstream N pool [6-6], if any
      if (this%upstream_is_varycn(irxn)) then
        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) - du_duc * &
            rt_auxvar%aqueous%dtotal(ispec_un,ispec_uc,iphase)
        else
          Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) - du_duc
        endif
      endif

! -- Jacobians with respect to upstream n ('un' or'u', due to variable CN ratio reactant)
      ! (TODO) currently not-yet-supported function.
      !

! -- Jacobians with respect to nh4, if any (nh4 as a reactant for N immoblization)
      if(this%upstream_is_varycn(irxn) .and. &
         this%mineral_n_stoich(irxn) < 0.0d0) then

        ! CO2 [6-1]
        Jacobian(ires_co2,ires_nh4) = Jacobian(ires_co2,ires_nh4) - dc_dn * &
          rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_nh4,iphase)
        !if(this%species_id_hrimm > 0) then   ! for tracking
        !  Jacobian(ires_hrimm,ires_nh4) = Jacobian(ires_hrimm,ires_nh4) - dc_dn
        !endif

        ! upstream C [6-2]
        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_uc,ires_nh4) = Jacobian(ires_uc,ires_nh4) - duc_dn* &
            rt_auxvar%aqueous%dtotal(ispec_uc,this%species_id_nh4,iphase)
        else
          Jacobian(ires_uc,ires_nh4) = Jacobian(ires_uc,ires_nh4) - duc_dn
        endif
 
        ! downstream C pools [6-3]
        do j = 1, this%n_downstream_pools(irxn)
           ispec_d = this%downstream_id(irxn, j)
           if(ispec_d < 0) then
             option%io_buffer = 'Downstream pool species not specified!'
             call printErrMsg(option)
           endif

           ! if given that: 'mineral_c_stoich + sum(di) = 1',
           !                 [6-1] dc_dx=dcrate_dx*mineral_c_stoich,
           !             and,[6-2] duc_dx=-dcrate_dx
           ! then, d(dci)/dx = -duc_dx * di
           ddc_dn = this%downstream_stoich(irxn, j) * (-1.d0*duc_dn)

           if(this%downstream_is_aqueous(irxn, j)) then
              ires_d = ispec_d
              Jacobian(ires_d,ires_nh4) = Jacobian(ires_d,ires_nh4) - ddc_dn* &
              rt_auxvar%aqueous%dtotal(ispec_d,this%species_id_nh4,iphase)
           else
              ires_d = reaction%offset_immobile + ispec_d
              Jacobian(ires_d,ires_nh4) = Jacobian(ires_d,ires_nh4) - ddc_dn
           endif
        enddo

        ! N[H3] [6-4]
        Jacobian(ires_nh4,ires_nh4) = Jacobian(ires_nh4,ires_nh4) - dn_dn * &
          rt_auxvar%aqueous%dtotal(this%species_id_nh4,this%species_id_nh4,iphase)
        !if(this%species_id_nimm > 0) then  ! for tracking
        !  Jacobian(ires_nimm,ires_nh4) = Jacobian(ires_nimm,ires_nh4) + dn_dn
        !endif

        ! NO3 [6-5], if any
        if (this%species_id_no3>0) then
          Jacobian(ires_no3,ires_nh4) = Jacobian(ires_no3,ires_nh4) - dno3_dn*  &
            rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_nh4,iphase)
          !!for tracking
          !if(this%species_id_nimm >0) then
          !  Jacobian(ires_nimm,ires_nh4) = Jacobian(ires_nimm,ires_nh4) + dno3_dn
          !endif
        endif

        ! upstream N pool [6-6], if any
        if (this%upstream_is_varycn(irxn)) then
          if(this%upstream_is_aqueous(irxn)) then
            Jacobian(ires_un,ires_nh4) = Jacobian(ires_un,ires_nh4) - du_dn * &
              rt_auxvar%aqueous%dtotal(ispec_un,this%species_id_nh4,iphase)
          else
            Jacobian(ires_un,ires_nh4) = Jacobian(ires_un,ires_nh4) - du_dn
          endif
        endif

      endif !if(this%upstream_is_varycn(irxn) .and. &
            !   this%mineral_n_stoich(irxn) < 0.0d0)

! -- Jacobians with respect to no3 (NO3 as a reactant for N immobilization)
      if(this%upstream_is_varycn(irxn) .and. &
         this%mineral_n_stoich(irxn) < 0.0d0 .and. &
         this%species_id_no3>0) then

        ! CO2 [6-1]
        Jacobian(ires_co2,ires_no3) = Jacobian(ires_co2,ires_no3) - dc_dno3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_no3,iphase)
        !if(this%species_id_hrimm > 0) then  ! for tracking
        !  Jacobian(ires_hrimm,ires_no3) = Jacobian(ires_hrimm,ires_no3) - dc_dno3
        !endif

        ! upstream C pool [6-2]
        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_uc,ires_no3) = Jacobian(ires_uc,ires_no3) - duc_dno3 * &
            rt_auxvar%aqueous%dtotal(ispec_uc,this%species_id_no3,iphase)
        else
          Jacobian(ires_uc,ires_no3) = Jacobian(ires_uc,ires_no3) - duc_dno3
        endif
 
        ! downstream C pools [6-3]
        do j = 1, this%n_downstream_pools(irxn)

           ispec_d = this%downstream_id(irxn, j)
           if(ispec_d < 0) then
             option%io_buffer = 'Downstream pool species not specified!'
             call printErrMsg(option)
           endif

           ! if given that: 'mineral_c_stoich + sum(di) = 1',
           !                 [6-1] dc_dx=dcrate_dx*mineral_c_stoich,
           !             and,[6-2] duc_dx=-dcrate_dx
           ! then, d(dci)/dx = -duc_dx * di
           ddc_dno3 = this%downstream_stoich(irxn, j) * (-1.d0*duc_dno3)

           if(this%downstream_is_aqueous(irxn, j)) then
              ires_d = ispec_d
              Jacobian(ires_d,ires_no3) = Jacobian(ires_d,ires_no3)-ddc_dno3 * &
                rt_auxvar%aqueous%dtotal(ispec_d,this%species_id_no3,iphase)
           else
              ires_d = reaction%offset_immobile + ispec_d
              Jacobian(ires_d,ires_no3) = Jacobian(ires_d,ires_no3)-ddc_dno3
           endif
        enddo

        ! N[H3] [6-4]
        Jacobian(ires_nh4,ires_no3) = Jacobian(ires_nh4,ires_no3) - dn_dno3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_nh4,this%species_id_no3,iphase)
        !if(this%species_id_nimm > 0) then  ! for tracking
        !  Jacobian(ires_nimm,ires_no3) = Jacobian(ires_nimm,ires_no3) + dn_dno3
        !endif

        ! NO3 [6-5]
        Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) - dno3_dno3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_no3,iphase)
        !if(this%species_id_nimm > 0) then  ! for tracking
        !  Jacobian(ires_nimm,ires_no3) = Jacobian(ires_nimm,ires_no3) + dno3_dno3
        !endif

        ! upstream N pool [6-6]
        Jacobian(ires_un,ires_no3) = Jacobian(ires_un,ires_no3) - du_dno3

      endif       !if(this%upstream_is_varycn(irxn) .and. &
                  !   this%mineral_n_stoich(irxn) < 0.0d0 .and. &
                  !   this%species_id_no3>0) then

    endif ! if(compute_derivative) then -- end of jacobian calculations

  enddo ! reactions loop

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(this%species_id_n2o>0) then
    ires_n2o = this%species_id_n2o + reaction%offset_aqueous       ! as aq. species
  endif
  if(this%species_id_ngasmin > 0) then
     ires_ngasmin = this%species_id_ngasmin + reaction%offset_immobile
  endif

  if(this%species_id_n2o>0 .and. net_nmin_rate>this%x0eps) then

    ! temperature/moisture/pH response functions (Parton et al. 1996)
    f_t = -0.06d0 + 0.13d0 * exp( 0.07d0 * tc )
    f_w = ((1.27d0 - saturation)/0.67d0)**(3.1777d0) * &
        ((saturation - 0.0012d0)/0.5988d0)**2.84d0  

    ph = 6.5d0       ! default
    if (this%species_id_proton > 0) then
      if (reaction%species_idx%h_ion_id > 0) then
        ph = &
          -log10(rt_auxvar%pri_molal(reaction%species_idx%h_ion_id)* &
                 rt_auxvar%pri_act_coef(reaction%species_idx%h_ion_id))
      else if (reaction%species_idx%h_ion_id < 0) then
        ph = &
          -log10(rt_auxvar%sec_molal(abs(reaction%species_idx%h_ion_id))* &
                 rt_auxvar%sec_act_coef(abs(reaction%species_idx%h_ion_id)))
      endif
    endif
    f_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + ph))/rpi

    if(f_t > this%x0eps .and. f_w > this%x0eps .and. f_ph > this%x0eps) then
      f_t = min(f_t, 1.0d0)
      f_w = min(f_w, 1.0d0)
      f_ph= min(f_ph, 1.0d0)
      temp_real = f_t * f_w * f_ph

      if(this%x0eps>0.d0) then
        ! GP's cut-off approach (sort of Heaviside function)
        feps0     = funcTrailersmooth(c_nh4, this%x0eps*10.d0, this%x0eps, PETSC_FALSE)
        dfeps0_dx = funcTrailersmooth(c_nh4, this%x0eps*10.d0, this%x0eps, PETSC_TRUE)
      else
        feps0 = 1.0d0
        dfeps0_dx = 0.d0
      endif
      ! constraining 'N potential rate' if too high compared to available within the allowable min. time-step
      ! It can be achieved by cutting time-step, but it may be taking a very small timestep finally
      ! - implying tiny timestep in model, which potentially crashes model

      ! in the following, '2.1' multiplier is chosen because that is slightly larger(5% to avoid numerical issue) than '2.0',
      ! which should be the previous-timestep before reaching the 'option%dt_min'
      ! (the time-cut used in PF is like dt=0.5*dt, when cutting)
      dtmin = 2.1d0*option%dt_min
      !dtmin = max(option%tran_dt, 2.1d0*option%dt_min)   ! this 'dt_min' may be accelerating the timing, but may not be good to mulitple consummers

      nratecap = temp_real*this%n2o_frac_mineralization*dtmin*net_nmin_rate
      if (nratecap > c_nh4) then
         ! assuming that monod type reduction used and 'nratecap' must be reduced to 'c_nh4', i.e.
         ! c_nh4/dt = nratecap*fnh4_inhibit_no3/dt * (c_nh4/(c_nh4+c0))
         ! then, c0 = nratecap*fnh4_inhibit_no3 - c_nh4 (this is the 'half-saturation' term used above)
         fnratecap = funcMonod(c_nh4, nratecap*fnh4_inhibit_no3-c_nh4, PETSC_FALSE)
         dfnratecap_dnh4 = funcMonod(c_nh4, nratecap*fnh4_inhibit_no3-c_nh4, PETSC_TRUE)
      else
         fnratecap       = 1.d0
         dfnratecap_dnh4 = 0.d0
      endif
      ! modifying the 'feps0' and 'dfeps0_dx' so that NO need to modify the major codes below
      dfeps0_dx = dfeps0_dx*fnratecap + feps0 * dfnratecap_dnh4   ! do the derivative first
      feps0 = feps0 * fnratecap

      !--------------------------------------------------------------
      ! residuals
      rate_n2o = temp_real*this%n2o_frac_mineralization*net_nmin_rate * feps0
 
      Residual(ires_nh4) = Residual(ires_nh4) + rate_n2o

      Residual(ires_n2o) = Residual(ires_n2o) - 0.5d0 * rate_n2o

      if(this%species_id_ngasmin > 0) then
         Residual(ires_ngasmin) = Residual(ires_ngasmin) - rate_n2o
      endif

     !Jacobians
      if (compute_derivative) then
        drate_n2o_dx = temp_real*this%n2o_frac_mineralization * &       ! 'constant' portion of 'rate_n2o'
                  (dnet_nmin_rate_dx * feps0                    &       ! d(net_nmin_rate*feps0)/dnh4
                  + net_nmin_rate * dfeps0_dx)

        Jacobian(ires_nh4,ires_nh4) = Jacobian(ires_nh4,ires_nh4) + drate_n2o_dx* &
           rt_auxvar%aqueous%dtotal(this%species_id_nh4,this%species_id_nh4,iphase)
        Jacobian(ires_n2o,ires_nh4) = Jacobian(ires_n2o,ires_nh4) - 0.5d0*drate_n2o_dx* &
           rt_auxvar%aqueous%dtotal(this%species_id_n2o,this%species_id_nh4,iphase)

! The following IS not needed and causes issue - breaking  whole reacton network if going wrong!
        !if(this%species_id_ngasmin > 0) then
        !   Jacobian(ires_ngasmin,ires_nh4) = Jacobian(ires_ngasmin,ires_nh4) - &
        !                              drate_n2o_dx
        !endif
      endif

    endif  ! end of 'f_t/f_w/f_ph > 1.0d-20'

  endif ! end of 'species_id_n2o > 0'

#ifdef DEBUG
  do ires=1, reaction%ncomp
    temp_real = Residual(ires)

    if (abs(temp_real) > huge(temp_real)) then
      write(option%fid_out, *) 'infinity of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: DECOMPOSITION'
      option%io_buffer = ' checking infinity of Residuals matrix  @ SomDecReact '
      call printErrMsg(option)
    endif
  enddo
#endif

end subroutine SomDecReact

! ************************************************************************** !
!
! SomDecDestroy: Destroys allocatable or pointer objects created in this module
! author: Guoping Tang
! ************************************************************************** !
subroutine SomDecDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_somdec_type) :: this
  
  type(pool_type), pointer :: cur_pool, prev_pool
  type(somdec_reaction_type), pointer :: cur_reaction, prev_reaction
  
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
  call DeallocateArray(this%rate_decomposition)
  call DeallocateArray(this%rate_ad_factor)
  call DeallocateArray(this%upstream_is_varycn)
  call DeallocateArray(this%upstream_c_id)
  call DeallocateArray(this%upstream_n_id)
  call DeallocateArray(this%upstream_nc)
  call DeallocateArray(this%upstream_is_aqueous)
  call DeallocateArray(this%downstream_id)
  call DeallocateArray(this%downstream_stoich)
  call DeallocateArray(this%downstream_is_aqueous)
  call DeallocateArray(this%downstream_is_varycn)
  call DeallocateArray(this%mineral_c_stoich) 
  call DeallocateArray(this%mineral_n_stoich) 
 
end subroutine SomDecDestroy

end module Reaction_Sandbox_SomDec_class
