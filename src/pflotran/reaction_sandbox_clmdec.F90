module Reaction_Sandbox_CLMDec_class

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

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_clm_decomp_type

    PetscInt :: temperature_response_function
    PetscInt :: moisture_response_function
    PetscReal :: Q10
    PetscReal :: half_saturation_nh3
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh3_no3
    PetscReal :: n2o_frac_mineralization    ! fraction of n2o from net N mineralization
    PetscReal :: x0eps

    PetscInt :: npool                       ! litter or variable CN ration pools
    PetscReal, pointer :: pool_nc_ratio(:)         ! NC ratio in mole  

    PetscInt :: nrxn
    PetscReal, pointer :: rate_constant(:)         !nrxn

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
    PetscInt :: species_id_nh3
    PetscInt :: species_id_no3
    PetscInt :: species_id_n2o
    PetscInt :: species_id_dom

    PetscInt :: species_id_hrimm
    PetscInt :: species_id_nmin
    PetscInt :: species_id_nimm
    PetscInt :: species_id_ngasmin
    PetscInt :: species_id_proton

    type(pool_type), pointer :: pools
    type(clm_decomp_reaction_type), pointer :: reactions
  contains
    procedure, public :: ReadInput => CLMDec_Read
    procedure, public :: Setup => CLMDec_Setup
    procedure, public :: Evaluate => CLMDec_React
    procedure, public :: Destroy => CLMDec_Destroy
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
  
  public :: CLMDec_Create

contains

! ************************************************************************** !

function CLMDec_Create()
  ! 
  ! Allocates CLMDec reaction object.
  ! 
  ! Author: Guoping Tang
  ! Date: 02/04/14
  ! 

  implicit none
  
  type(reaction_sandbox_clm_decomp_type), pointer :: CLMDec_Create
  
  allocate(CLMDec_Create)

#ifdef CLM_PFLOTRAN
  CLMDec_Create%temperature_response_function=TEMPERATURE_RESPONSE_FUNCTION_CLM4
  CLMDec_Create%moisture_response_function = MOISTURE_RESPONSE_FUNCTION_CLM4
#endif

  CLMDec_Create%Q10 = 1.5d0
  CLMDec_Create%half_saturation_nh3 = 1.0d-15
  CLMDec_Create%half_saturation_no3 = 1.0d-15
  CLMDec_Create%inhibition_nh3_no3 = 1.0d-15
  CLMDec_Create%n2o_frac_mineralization = 0.02d0  ! Parton et al. 2001
  CLMDec_Create%x0eps = 1.0d-20

  CLMDec_Create%npool = 0
  nullify(CLMDec_Create%pool_nc_ratio)

  CLMDec_Create%nrxn = 0
  nullify(CLMDec_Create%rate_constant)
  nullify(CLMDec_Create%upstream_is_varycn)
  nullify(CLMDec_Create%upstream_c_id)
  nullify(CLMDec_Create%upstream_n_id)
  nullify(CLMDec_Create%upstream_nc)
  nullify(CLMDec_Create%upstream_is_aqueous)
  
  nullify(CLMDec_Create%n_downstream_pools)
  nullify(CLMDec_Create%downstream_id)
  nullify(CLMDec_Create%downstream_is_aqueous)
  nullify(CLMDec_Create%downstream_stoich)
  nullify(CLMDec_Create%downstream_is_varycn)
  nullify(CLMDec_Create%mineral_c_stoich)
  nullify(CLMDec_Create%mineral_n_stoich)

  CLMDec_Create%species_id_co2 = 0
  CLMDec_Create%species_id_nh3 = 0
  CLMDec_Create%species_id_no3 = 0
  CLMDec_Create%species_id_n2o = 0
  CLMDec_Create%species_id_dom = 0
  CLMDec_Create%species_id_hrimm = 0
  CLMDec_Create%species_id_nmin = 0
  CLMDec_Create%species_id_nimm = 0
  CLMDec_Create%species_id_ngasmin = 0

  nullify(CLMDec_Create%next)

  nullify(CLMDec_Create%pools)
  nullify(CLMDec_Create%reactions)

end function CLMDec_Create

! ************************************************************************** !

subroutine CLMDec_Read(this,input,option)
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
                       'CHEMISTRY,REACTION_SANDBOX,CLMDec')
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
            'CHEMISTRY,REACTION_SANDBOX,CLMDec,TEMPERATURE RESPONSE FUNCTION')
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
                       'REACTION_SANDBOX_CLMDec,TEMPERATURE RESPONSE FUNCTION')
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMDec,' // &
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
                       'CHEMISTRY,REACTION_SANDBOX,CLMDec,MOISTURE RESPONSE FUNCTION')
         call StringToUpper(word)   

            select case(trim(word))
              case('CLM4')
                  this%moisture_response_function = MOISTURE_RESPONSE_FUNCTION_CLM4    
              case('DLEM')
                  this%moisture_response_function = MOISTURE_RESPONSE_FUNCTION_DLEM    
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMDec,TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
         enddo 

#endif

     case('X0EPS')
         call InputReadDouble(input,option,this%x0eps)
         call InputErrorMsg(input,option,'x0eps', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

     case('AMMONIA_HALF_SATURATION')
         call InputReadDouble(input,option,this%half_saturation_nh3)
         call InputErrorMsg(input,option,'ammonia half saturation', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

     case('NITRATE_HALF_SATURATION')
         call InputReadDouble(input,option,this%half_saturation_no3)
         call InputErrorMsg(input,option,'nitrate half saturation', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

     case('AMMONIA_INHIBITION_NITRATE')
         call InputReadDouble(input,option,this%inhibition_nh3_no3)
         call InputErrorMsg(input,option,'ammonia inhibition on nitrate immobilization', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

     case('N2O_FRAC_MINERALIZATION')
         call InputReadDouble(input,option,this%n2o_frac_mineralization)
         call InputErrorMsg(input,option,'n2o fraction from mineralization', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')

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
            'CHEMISTRY,REACTION_SANDBOX,CLMDec,POOLS')
          call InputReadDouble(input,option,temp_real)
          if (InputError(input)) then
            new_pool%nc_ratio = -999.d0
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
                             'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('UPSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%upstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'upstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
            case('DOWNSTREAM_POOL')
              allocate(new_pool_rxn)
              new_pool_rxn%name = ''
              new_pool_rxn%stoich = 0.d0
              nullify(new_pool_rxn%next)

              call InputReadWord(input,option, &
                                 new_pool_rxn%name,PETSC_TRUE)
              call InputErrorMsg(input,option,'downstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
              call InputReadDouble(input,option,new_pool_rxn%stoich)
              call InputErrorMsg(input,option,'Downstream pool stoich', 'CHEMISTRY,' // &
                  'REACTION_SANDBOX_CLMDec,REACTION')

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
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLMDec RATE CONSTANT UNITS'
                call InputDefaultMsg(input,option)
              else              
                rate_constant = rate_constant * &
                  UnitsConvertToInternal(word,option)
              endif
            case('TURNOVER_TIME')
              call InputReadDouble(input,option,turnover_time)
              call InputErrorMsg(input,option,'turnover time', &
                     'CHEMISTRY,REACTION_SANDBOX,CLMDec,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLMDec TURNOVER TIME UNITS'
                call InputDefaultMsg(input,option)
              else              
                turnover_time = turnover_time * &
                  UnitsConvertToInternal(word,option)
              endif
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMDec,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
              call printErrMsg(option)
          end select
        enddo
        
        ! check to ensure that one of turnover time or rate constant is set.
        if (turnover_time > 0.d0 .and. rate_constant > 0.d0) then
          option%io_buffer = 'Only TURNOVER_TIME or RATE_CONSTANT may ' // &
            'be included in a CLMDec reaction definition, but not both. ' // &
            'See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        else if (turnover_time > 0.d0) then
          new_reaction%rate_constant = 1.d0 / turnover_time
        else
          new_reaction%rate_constant = rate_constant
        endif
        if (associated(this%reactions)) then
          prev_reaction%next => new_reaction
        else
          this%reactions => new_reaction
        endif
        prev_reaction => new_reaction
        nullify(new_reaction)        
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLMDec keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine CLMDec_Read

! ************************************************************************** !

subroutine CLMDec_Setup(this,reaction,option)
  ! 
  ! Sets up CLMDec reaction after it has been read from input
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
  this%upstream_is_varycn = PETSC_FALSE
  this%upstream_c_id = 0
  this%upstream_n_id = 0
  this%upstream_nc = -999.9
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
        option%io_buffer = 'For CLMDec pools with no CN ratio defined, ' // &
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
             option%io_buffer = 'CLMDec pool: ' // cur_pool%name // &
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
            option%io_buffer = 'For CLMDec reactions, downstream pool ' // &
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
         option%io_buffer = 'CLMDec fixed-CN C pool decomposition '// &
          'has a negative respiration fraction!' // &
            'Please check reactions of upstream pool: ' // &
             trim(cur_pool%name)
         call printErrMsg(option)
       endif

       if(stoich_n < 0.0d0) then
         option%io_buffer = 'CLMDec fixed-CN SOM N mineralization is negative,' // &
           'i.e., fixed-CN SOM-N immobilization is not supposed currently,' // &
           'and, variable-CN ratio suggested for SOM pool: ' // &
             trim(cur_pool%name)
         call printErrMsg(option)
       endif

       this%mineral_c_stoich(icount) = stoich_c  ! this will be updated, if variable C/N Organic C pool(s) exists.
       this%mineral_n_stoich(icount) = stoich_n  ! this will be updated, if variable C/N Organic C pool(s) exists.

     endif
  enddo

  word = 'HCO3-'
  this%species_id_co2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if(this%species_id_co2 < 0) then
     word = 'CO2(aq)'
     this%species_id_co2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  endif

  word = 'NH4+'
  this%species_id_nh3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if(this%species_id_nh3 < 0) then
    word = 'NH3(aq)'
    this%species_id_nh3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  endif

  if(this%species_id_nh3 <= 0) then
    option%io_buffer = 'NH4+ or NH3(aq) is not specified in the input file!'
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

  word = 'Nmin'
  this%species_id_nmin = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  word = 'Nimm'
  this%species_id_nimm = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  word = 'NGASmin'
  this%species_id_ngasmin = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)

end subroutine CLMDec_Setup

! ************************************************************************** !
subroutine CLMDec_React(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
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

  class(reaction_sandbox_clm_decomp_type) :: this
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

  PetscReal :: temp_real
  PetscReal :: c_uc, c_un    ! concentration (mole/m3 or mole/L, upon species type)
  PetscReal :: c_nh3         ! concentration (mole/L): substrate OR product
  PetscReal :: c_no3         ! concentration (mole/L): substrate only

  PetscInt :: irxn
  PetscInt :: ipool_up, ipool_down

  PetscInt :: ispec_uc, ispec_un, ispec_d   ! species id for upstream C, N, and downstream
  PetscInt :: ires_uc, ires_un, ires_d      ! id used for residual and Jacobian
  PetscInt :: ires_co2, ires_nh3, ires_n2o, ires_no3
  PetscInt :: ires_hrimm, ires_nmin, ires_nimm, ires_ngasmin
  PetscReal :: stoich_c, stoich_n

  PetscReal :: scaled_crate_const
  PetscReal :: tc     ! temperature in degC
  PetscReal :: f_t    ! temperature response function
  PetscReal :: f_w    ! moisture response function

  PetscReal :: crate       ! moleC/s: overall = crate_uc, or,
                           !                  = crate_uc*crate_nh3, or,
                           !                  = crate_uc*(crate_nh3*(1-fnh3_inhibit_no3)+crate_no3*fnh3_inhibit_no3)
  PetscReal :: dcrate_dx   ! d(crate)/d(x), x = uc, un, dc, dn, nh3, or no3, respectively

  PetscReal :: crate_uc        ! crate function of upstream C ('uc')
  PetscReal :: dcrate_uc_duc   ! d crate_uc / d uc

  PetscReal :: crate_nh3        ! crate function of c_nh3 (for nh3 immobilization)
  PetscReal :: dcrate_nh3_dnh3  ! d(crate_nh3)/d(nh3) (for nh3 immobilization)
  PetscReal :: fnh3             ! = nh3/(half_saturation + nh3)  ( N 'resource' limitation on immobilization)
  PetscReal :: dfnh3_dnh3       ! d(fnh3)/d(nh3)

  PetscReal :: crate_no3        ! crate function of c_no3 (for no3 immobilization)
  PetscReal :: dcrate_no3_dno3  ! d(crate_no3)/d(no3)
  PetscReal :: fno3             ! = no3/(half_saturation + no3)  ( N 'resource' limitation on immobilization)
  PetscReal :: dfno3_dno3       ! d(fnh3)/d(nh3)

  !nh3 inhibition on no3 immobilization, or microbial N immobilization preference btw nh3 and no3
  PetscReal :: fnh3_inhibit_no3 ! inhibition_coef/(inhibition_coef + nh3):
  PetscReal :: dfnh3_inhibit_no3_dnh3 ! d(fnh3_inhibit_no3)/dnh3

  ! save mineral N fraction and decomposition rate for net N mineralization and associated N2O calculation
  PetscReal :: net_nmin_rate
  PetscReal :: dnet_nmin_rate_dx
  PetscReal :: ph, f_ph
  PetscReal :: rate_n2o, drate_n2o_dx

! other local variables
  ! UC + u*UN --> (1-di)*DCi +(u-ni-n2)*DNi+ d*CO2 + n*N[H3] [+ n2*NO3], and
  ! the reaction (decomposition) rate = crate, and independent unkowns: d, n, and possibly u, n2. ('di','ni' imply multi-downpools)
  ! n1 may be negative, i.e., n1*NH3 should be in the left side of the Eq. (N immobilization).
  ! n2 must be negative, i.e., n2*NO3 should be in the left side of the Eq. (N immobilization).
  ! u = (UC/UN), uc/un for upstream C/N (UC/UN here), n for NH3, c for CO2
  ! So, total 7 variables, and 6 possible independent variables.
  ! but, (1) 'u-ni-n2' is implicitly completely dependent on others currently, so ignored in the code. So 6 variables;
  !      (2) any rates NOW are not functions of down-stream C-N pools, So 4 possible independent variables.
  PetscReal :: dc_duc, duc_duc, ddc_duc, dn_duc, dno3_duc, du_duc        ! 'uc' as independent variable
  PetscReal :: dc_dun, duc_dun, ddc_dun, dn_dun, dno3_dun, du_dun        ! 'un' as independent variable (NOT YET either - TODO)

  PetscReal :: dc_dn, duc_dn, ddc_dn, dn_dn, dno3_dn, du_dn              ! 'nh3' as independent variable
  PetscReal :: dc_dno3, duc_dno3, ddc_dno3, dn_dno3, dno3_dno3, du_dno3  ! 'no3' as independent variable

! misc. local variables
  PetscInt :: i, j, ires
  PetscReal:: feps0, dfeps0_dx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !------------------------------------------------------------------------------------------------
  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  ! -----------------------------------------------------------------------------------------------
  ! 'nh3' dependent rate function (DON'T change the 'rate' and 'derivatives' after this)
  c_nh3     = rt_auxvar%total(this%species_id_nh3, iphase)
  temp_real = c_nh3 + this%half_saturation_nh3
  fnh3      = c_nh3 / temp_real
  dfnh3_dnh3 = this%half_saturation_nh3 / temp_real / temp_real

  if (this%species_id_no3 > 0) then
    ! 'no3' dependent rate function (DON'T change the 'rate' and 'derivatives' after this)
    c_no3 = rt_auxvar%total(this%species_id_no3, iphase)
    temp_real = c_no3 + this%half_saturation_no3
    fno3      = c_no3 / temp_real
    dfno3_dno3= this%half_saturation_no3 / temp_real / temp_real

    ! nh3 inhibition on no3 immobilization, if any ('this%inhibition_nh3_no3')
    ! (DON'T change the 'rate' and 'derivatives' after this)
    temp_real = this%inhibition_nh3_no3 + c_nh3
    fnh3_inhibit_no3 = this%inhibition_nh3_no3/temp_real
    dfnh3_inhibit_no3_dnh3 = -this%inhibition_nh3_no3/temp_real/temp_real     ! over 'd_nh3'

  endif ! if (this%species_id_no3>0)

  !----------------------------------------------------------------------------------------------
  ires_co2 = this%species_id_co2
  ires_nh3 = this%species_id_nh3
  ires_no3 = this%species_id_no3
  ires_n2o = this%species_id_n2o

  if(this%species_id_hrimm > 0) then
     ires_hrimm = this%species_id_hrimm + reaction%offset_immobile
  endif

  if(this%species_id_nmin > 0) then
     ires_nmin = this%species_id_nmin + reaction%offset_immobile
  endif

  if(this%species_id_nimm > 0) then
     ires_nimm = this%species_id_nimm + reaction%offset_immobile
  endif

  if(this%species_id_ngasmin > 0) then
     ires_ngasmin = this%species_id_ngasmin + reaction%offset_immobile
  endif

  !----------------------------------------------------------------------------------------
  ! temperature response function
  tc = global_auxvar%temp

#ifdef CLM_PFLOTRAN 
  f_t = GetTemperatureResponse(tc,this%temperature_response_function, this%Q10) 
#else
  f_t = 1.0d0
#endif

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

  if(f_t < 1.0d-20 .or. f_w < 1.0d-20) then
     return
  endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  net_nmin_rate     = 0.0d0
  dnet_nmin_rate_dx = 0.0d0

  do irxn = 1, this%nrxn
  
    !-----------------------------------------------------------------------------------------------------
    ! calculate C rate (N rate is derived by Crate*N/C)

    ! scaled_rate_const units: (m^3 bulk / s) = (1/s) * (m^3 bulk)
    scaled_crate_const = this%rate_constant(irxn)*volume*f_t*f_w

    ! substrates
    ispec_uc = this%upstream_c_id(irxn)
    if(this%upstream_is_aqueous(irxn)) then
      c_uc = rt_auxvar%total(ispec_uc, iphase)
      c_uc = theta * 1000.0d0 * c_uc    ! from mol/L -> mol/m3
      ires_uc = ispec_uc
    else
      c_uc = rt_auxvar%immobile(ispec_uc)
      ires_uc = reaction%offset_immobile + ispec_uc
    endif

    !if(c_uc <= this%x0eps) cycle     ! this may bring in 'oscillation' around 'this%x0eps'
    if(this%x0eps>0.d0) then
      feps0 = c_uc/(c_uc+this%x0eps)    ! using these two for trailer smoothing, alternatively
      dfeps0_dx = this%x0eps/(c_uc+this%x0eps)/(c_uc+this%x0eps)
    else
      feps0 = 1.0d0
      dfeps0_dx = 0.d0
    endif

    ! C substrate only dependent rate/derivative  (DON'T change after this)
    crate_uc  = scaled_crate_const * c_uc * feps0
    dcrate_uc_duc = scaled_crate_const * (feps0 + c_uc*dfeps0_dx)

    !-----------------------------------------------------------------------------------------------------

    ! for litter decomposition reactions, N/C stoichiometry needs to be calculated on the fly
    if(this%upstream_is_varycn(irxn)) then

      ispec_un = this%upstream_n_id(irxn)
      ires_un = ispec_un + reaction%offset_immobile

      if(ispec_un > 0) then
         c_un = rt_auxvar%immobile(ispec_un)
         this%upstream_nc(irxn) = c_un / c_uc
      endif

      ! calculate respiration factor (CO2 stoichiometry)
      stoich_c = 1.0d0
      do j = 1, this%n_downstream_pools(irxn)
        stoich_c = stoich_c - this%downstream_stoich(irxn, j)
      enddo

      if(stoich_c < 0.0d0) then
        option%io_buffer = 'CLMDec variable-CNratio C pool decomposition has' // &
                              'a negative respiration fraction!'
        call printErrMsg(option)
      endif
      this%mineral_c_stoich(irxn) = stoich_c

      ! calculate N (NH3) stoichiometry
      stoich_n = this%upstream_nc(irxn)
      do j = 1, this%n_downstream_pools(irxn)
        stoich_n = stoich_n - this%downstream_stoich(irxn, j) * &
                              this%downstream_nc(irxn, j)
      enddo
      this%mineral_n_stoich(irxn) = stoich_n

    endif ! end of 'upstream_is_varycn' (NOTE: 'mineral_n_stoich' for non-litter species IS constant, which are calc. in 'setup')

    !-----------------------------------------------------------------------------------------------------
    ! calculation of residuals
    ! residual units: (mol/sec) = (m^3 bulk/s) * (mol/m^3 bulk)

    ! overall C rates
    crate     = crate_uc

    ! NH3/NO3 limiting on N-immobolization involved C rates, if any
    if(this%mineral_n_stoich(irxn) < 0.0d0) then
      crate_nh3 = fnh3
      ! overall C rates
      crate = crate_uc * crate_nh3

      if (this%species_id_no3 > 0) then
        ! since 'f_nh3_inhibit_no3' is somehow a spliting fraction for NO3 immoblization,
        ! '1.-f_nh3_inhibit_no3' should be those for NH4 immobilization
        crate_nh3 = fnh3*(1.d0-fnh3_inhibit_no3)

        ! f_no3 should be adjusted by f_nh3_inhibit_no3 so that it is inhibited by nh3
        crate_no3 = fno3*fnh3_inhibit_no3

        ! overall C rates
        crate = crate_uc * (crate_nh3 + crate_no3)
      endif
    endif

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

    ! inorg. nitrogen (NH4+/NH3(aq)) [4], [and, no3 [5] ]
    if(this%mineral_n_stoich(irxn) >= 0.0d0) then        ! mineralization
      Residual(ires_nh3) = Residual(ires_nh3) - this%mineral_n_stoich(irxn) * crate

      if(this%species_id_nmin > 0) then    ! for tracking
        Residual(ires_nmin) = Residual(ires_nmin) - this%mineral_n_stoich(irxn) * crate
      endif

    else                                                 ! immobilization

      Residual(ires_nh3) = Residual(ires_nh3) - &
        this%mineral_n_stoich(irxn) * crate_uc * crate_nh3
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

! -- derivatives with fixed C/N decomposing pools AND no reduction of NH3 (positive 'n')

        ! SOMc + u SOMn -> (1-di) SOMci + (u-n) SOMni + d C[O2] + n N[H3]
        ! [2]    [6]       [3]            [ignored]     [1]       [4]
        ! d - minearl_c_stoich, n - mineral_n_stoich, u - upstream_nc, 1-di - (multi-)downstream_stoich
        ! NOTE: (1) '(u-n)' is implicitly calc'ed (totally ignored in the codes, so NO 'ddn_d??' in the codes)
        !       (2) 'SOMci' will be caclucated below, because multiple-downstream pools('ddc_d??') are allowed
        !       (3) 'u' here is implicitly calc'ed, but not for variale-CN pools
        !       (4) currently, rates NOT depends upon 'un' or 'u', so,
        dc_dun  = 0.d0     ! [6-1] 'C rates'
        duc_dun = 0.d0     ! [6-2]
        ddc_dun = 0.d0     ! [6-3 - will be done in Jacobians)

        dn_dun  = 0.d0     ! [6-4] 'N rates'
        dno3_dun= 0.d0     ! [6-5]
        du_dun  = 0.d0     ! [6-6]

        ! so, 'uc' is the only independent variable, by this point
        dcrate_dx = dcrate_uc_duc
        dc_duc  = dcrate_dx * this%mineral_c_stoich(irxn)              ! [6-1]
        duc_duc = -1.0d0*dcrate_dx                                     ! [6-2]

        dn_duc  = dc_duc * this%mineral_n_stoich(irxn)                 ! [6-4]

! -- derivatives with variable C/N decompsing pools
        ! LitC + u LitN -> (1-di) SOMci + (u-n-n2) SOMni + d C[O2] + n N[H3] [+ n2 NO3]
        ! [2]    [6]       [3]            [ignored]          [1]       [4]        [5]
        ! 'n+n2' IS mineral_n_stoich
        ! Note: (1) only 'u-n-n2' is implicitly calc'ed (i.e., ignored in the codes)

        if(this%upstream_is_varycn(irxn)) then

           ! first, considering variable 'u', but not limited for 'crate_uc' (i.e., 'n' is still non-negative)
           du_duc = this%upstream_nc(irxn) * duc_duc                                    ! [6-6]

! -- derivatives with variable C/N decompsing pools and N[H3] as reactants
           ! 'un' is limited, so 'n' is negative (NH3 as reactant, i.e. N immobilization)
           if (this%mineral_n_stoich(irxn) < 0.d0) then

             ! -- depending on 'uc'
             ! crate = crate_uc * crate_nh3 = crate_uc*fnh3,
             ! So, d(crate)/d(uc) = fnh3*dcrate_uc_duc + crate_uc * dfnh3_duc
             !     in which 'dfnh3_duc' = 0, because crate_nh3 = fnh3, only dependent on 'nh3'
             dcrate_dx = dcrate_uc_duc*fnh3
             dc_duc  =  dcrate_dx * this%mineral_c_stoich(irxn)                         ! [6-1]
             duc_duc = -1.0d0 * dcrate_dx                                               ! [6-2]

             dn_duc  = dc_duc * this%mineral_n_stoich(irxn)                             ! [6-4]
             du_duc  = this%upstream_nc(irxn) * duc_duc                                 ! [6-6]

             ! -- depending on 'nh3'
             ! d(crate)/d(n) = d(crate_nh3)*crate_uc = dfnh3_dnh3*crate_uc
             dcrate_dx = dfnh3_dnh3*crate_uc
             dc_dn  = dcrate_dx * this%mineral_c_stoich(irxn)                           ! [6-1]
             duc_dn = -1.0d0*dcrate_dx                                                  ! [6-2]

             dn_dn  = dc_dn * this%mineral_n_stoich(irxn)                               ! [6-4]
             du_dn  = this%upstream_nc(irxn) * duc_dn                                   ! [6-6]

! -- derivatives with variable C/N decompsing pools and N[H3]+NO3 as reactants
             if (this%species_id_no3 > 0) then
               ! -- depending on 'uc'
               ! crate = crate_uc * (crate_nh3 + crate_no3)
               !       = crate_uc * (fnh3*(1-fnh3_inhibit_no3)+fno3*fnh3_inhibit_no3),
               ! So, d(crate)/d(uc) = (fnh3*(1-fnh3_inhibit_no3)+fno3*fnh3_inhibit_no3)*dcrate_uc_duc
               dcrate_dx  = dcrate_uc_duc * &
                         (fnh3*(1.d0-fnh3_inhibit_no3)+fno3*fnh3_inhibit_no3)

               dc_duc  = dcrate_dx * this%mineral_c_stoich(irxn)                         ! [6-1]

               duc_duc = -1.0d0 * dcrate_dx                                              ! [6-2]

               ! nrate = crate_uc*fnh3*(1-fnh3_inhibit_no3)*mineral_n_stoich: 'mineral_n_stoich = u-sum(ni)'
               dn_duc  = dcrate_uc_duc*fnh3*(1.d0-fnh3_inhibit_no3) &
                           * this%mineral_n_stoich(irxn)                                 ! [6-4]

               ! no3rate = crate_uc*fno3*fnh3_inhibit_no3*mineral_n_stoich: 'mineral_n_stoich = u-sum(ni)'
               dno3_duc= dcrate_uc_duc*fno3*fnh3_inhibit_no3 &
                           * this%mineral_n_stoich(irxn)                                 ! [6-5]

               du_duc  = this%upstream_nc(irxn) * duc_duc                                ! [6-6]

               ! -- depending on 'nh3'
               ! d(crate)/d(n) = d(crate_nh3+dcrate_no3)*crate_uc/dnh3
               !               = d(fnh3*(1-fnh3_inhibit_no3)+fno3*fnh3_inhibit_no3)/dnh3*crate_uc
               !               = ( dfnh3_dnh3*(1-fnh3_inhibit_no3)
               !                  +fnh3*(-dfnh3_inhibit_no3_dnh3)
               !                  +fno3*dfnh3_inhibit_no3_dnh3
               !                 ) * crate_uc
               dcrate_dx  =  dfnh3_dnh3*(1.d0-fnh3_inhibit_no3)  &
                           + fnh3*(-dfnh3_inhibit_no3_dnh3)      &
                           + fno3*dfnh3_inhibit_no3_dnh3
               dcrate_dx  = dcrate_dx*crate_uc

               dc_dn  = dcrate_dx * this%mineral_c_stoich(irxn)                          ! [6-1]

               duc_dn = -1.0d0 * dcrate_dx                                               ! [6-2]

               ! nrate = crate_uc*fnh3*(1-fnh3_inhibit_no3)*mineral_n_stoich: 'mineral_n_stoich = u-sum(ni)'
               ! d(nrate)/d(n) = ( fnh3*(-dfnh3_inhibit_no3_dnh3)
               !                  +dfnh3_dnh3*(1-fnh3_inhibit_no3)
               !                 ) * crate_uc * mineral_n_stoich
               dn_dn  = fnh3*(-dfnh3_inhibit_no3_dnh3) &
                       +dfnh3_dnh3*(1.d0-fnh3_inhibit_no3)
               dn_dn  = dn_dn * crate_uc*this%mineral_n_stoich(irxn)                     ! [6-4]

               ! no3rate = crate_uc*fno3*fnh3_inhibit_no3*mineral_n_stoich: 'mineral_n_stoich = u-sum(ni)'
               ! d(no3rate)/d(n) = ( fno3*dfnh3_inhibit_no3_dnh3)
               !                    +dfno3_dnh3*fnh3_inhibit_no3         ! 'dfno3_dnh3' = 0
               !                 ) * crate_uc * mineral_n_stoich
               dno3_dn= fno3*dfnh3_inhibit_no3_dnh3  &
                         * crate_uc*this%mineral_n_stoich(irxn)                          ! [6-5]

               du_dn  = this%upstream_nc(irxn) * duc_dn                                  ! [6-6]

               ! -- depending on 'no3'
               ! d(crate)/d(no3) = d(crate_nh3+dcrate_no3)*crate_uc/dno3
               !               = d(fnh3*(1-fnh3_inhibit_no3)+fno3*fnh3_inhibit_no3)/dno3*crate_uc
               !               = ( dfnh3_dno3*(1-fnh3_inhibit_no3)       ! 'dfnh3_dno3' = 0
               !                  +fnh3*(-dfnh3_inhibit_no3_dno3)        ! 'dfnh3_inhibit_no3_dno3' = 0
               !                  +fno3*dfnh3_inhibit_no3_dno3           ! 'dfnh3_inhibit_no3_dno3' = 0
               !                  +dfno3_dno3*fnh3_inhibit_no3
               !                 ) * crate_uc
               dcrate_dx = dfno3_dno3*fnh3_inhibit_no3 * crate_uc

               dc_dno3  = dcrate_dx * this%mineral_c_stoich(irxn)                        ! [6-1]

               duc_dno3 = -1.0d0*dcrate_dx                                               ! [6-2]

               ! nrate = crate_uc*fnh3*(1-fnh3_inhibit_no3)*mineral_n_stoich: 'mineral_n_stoich = u-sum(ni)'
               ! d(nrate)/d(no3) = ( fnh3*(-dfnh3_inhibit_no3_dno3)      ! 'dfnh3_inhibit_no3_dno3' = 0
               !                    +dfnh3_dno3*(1-fnh3_inhibit_no3)     ! 'dfnh3_dno3' = 0
               !                   ) * crate_uc * mineral_n_stoich
               dn_dno3  = 0.d0                                                           ! [6-4]

               ! no3rate = crate_uc*fno3*fnh3_inhibit_no3*mineral_n_stoich: 'mineral_n_stoich = u-sum(ni)'
               ! d(no3rate)/d(no3) = ( fno3*dfnh3_inhibit_no3_dno3)      ! 'dfnh3_inhibit_no3_dno3' = 0
               !                    +dfno3_dno3*fnh3_inhibit_no3
               !                 ) * crate_uc * mineral_n_stoich
               dno3_dno3= dfno3_dno3 * fnh3_inhibit_no3 &
                         * crate_uc * this%mineral_n_stoich(irxn)                        ! [6-5]

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
      if(this%species_id_hrimm > 0) then    ! for tracking
        Jacobian(ires_hrimm,ires_uc) = Jacobian(ires_hrimm,ires_uc) - dc_duc
      endif

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
        Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - dn_duc*  &
          rt_auxvar%aqueous%dtotal(this%species_id_nh3,ispec_uc,iphase)
      else
        Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - dn_duc
      endif
      !for tracking
      if(this%species_id_nmin > 0 .and. this%mineral_n_stoich(irxn) > 0.0d0) then
        Jacobian(ires_nmin,ires_uc) = Jacobian(ires_nmin,ires_uc) - dn_duc
      endif
      if(this%species_id_nimm > 0 .and. this%mineral_n_stoich(irxn) < 0.0d0) then
        Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) + dn_duc
      endif

      ! NO3 [6-5], if any
      if (this%species_id_no3>0 .and. &
          this%upstream_is_varycn(irxn) .and. &
          this%mineral_n_stoich(irxn) < 0.0d0) then

        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - &
            dno3_duc*  &
            rt_auxvar%aqueous%dtotal(this%species_id_no3,ispec_uc,iphase)
        else
          Jacobian(ires_nh3,ires_uc) = Jacobian(ires_nh3,ires_uc) - &
            dno3_duc
        endif
        !for tracking
        if(this%species_id_nimm >0) then
          Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) + dn_duc
        endif
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

! -- Jacobians with respect to nh3, if any (nh3 as a reactant for N immoblization)
      if(this%upstream_is_varycn(irxn) .and. &
         this%mineral_n_stoich(irxn) < 0.0d0) then

        ! CO2 [6-1]
        Jacobian(ires_co2,ires_nh3) = Jacobian(ires_co2,ires_nh3) - dc_dn * &
          rt_auxvar%aqueous%dtotal(this%species_id_co2,this%species_id_nh3,iphase)
        if(this%species_id_hrimm > 0) then   ! for tracking
          Jacobian(ires_hrimm,ires_nh3) = Jacobian(ires_hrimm,ires_nh3) - dc_dn
        endif

        ! upstream C [6-2]
        if(this%upstream_is_aqueous(irxn)) then
          Jacobian(ires_uc,ires_nh3) = Jacobian(ires_uc,ires_nh3) - duc_dn* &
            rt_auxvar%aqueous%dtotal(ispec_uc,this%species_id_nh3,iphase)
        else
          Jacobian(ires_uc,ires_nh3) = Jacobian(ires_uc,ires_nh3) - duc_dn
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
              Jacobian(ires_d,ires_nh3) = Jacobian(ires_d,ires_nh3) - ddc_dn* &
              rt_auxvar%aqueous%dtotal(ispec_d,this%species_id_nh3,iphase)
           else
              ires_d = reaction%offset_immobile + ispec_d
              Jacobian(ires_d,ires_nh3) = Jacobian(ires_d,ires_nh3) - ddc_dn
           endif
        enddo

        ! N[H3] [6-4]
        Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) - dn_dn * &
          rt_auxvar%aqueous%dtotal(this%species_id_nh3,this%species_id_nh3,iphase)
        if(this%species_id_nimm > 0) then  ! for tracking
          Jacobian(ires_nimm,ires_nh3) = Jacobian(ires_nimm,ires_nh3) + dn_dn
        endif

        ! NO3 [6-5], if any
        if (this%species_id_no3>0) then
          Jacobian(ires_no3,ires_nh3) = Jacobian(ires_no3,ires_nh3) - dno3_dn*  &
            rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_nh3,iphase)
          !for tracking
          if(this%species_id_nimm >0) then
            Jacobian(ires_nimm,ires_nh3) = Jacobian(ires_nimm,ires_nh3) + dno3_dn
          endif
        endif

        ! upstream N pool [6-6], if any
        if (this%upstream_is_varycn(irxn)) then
          if(this%upstream_is_aqueous(irxn)) then
            Jacobian(ires_un,ires_nh3) = Jacobian(ires_un,ires_nh3) - du_dn * &
              rt_auxvar%aqueous%dtotal(ispec_un,this%species_id_nh3,iphase)
          else
            Jacobian(ires_un,ires_nh3) = Jacobian(ires_un,ires_nh3) - du_dn
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
        if(this%species_id_hrimm > 0) then  ! for tracking
          Jacobian(ires_hrimm,ires_no3) = Jacobian(ires_hrimm,ires_no3) - dc_dno3
        endif

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
        Jacobian(ires_nh3,ires_no3) = Jacobian(ires_nh3,ires_no3) - dn_dno3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_nh3,this%species_id_no3,iphase)
        if(this%species_id_nimm > 0) then  ! for tracking
          Jacobian(ires_nimm,ires_no3) = Jacobian(ires_nimm,ires_no3) + dn_dno3
        endif

        ! NO3 [6-5]
        Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) - dno3_dno3 * &
          rt_auxvar%aqueous%dtotal(this%species_id_no3,this%species_id_no3,iphase)
        if(this%species_id_nimm > 0) then  ! for tracking
          Jacobian(ires_nimm,ires_no3) = Jacobian(ires_nimm,ires_no3) + dno3_dno3
        endif

        ! upstream N pool [6-6]
        Jacobian(ires_un,ires_no3) = Jacobian(ires_un,ires_no3) - du_dno3

      endif       !if(this%upstream_is_varycn(irxn) .and. &
                  !   this%mineral_n_stoich(irxn) < 0.0d0 .and. &
                  !   this%species_id_no3>0) then

    endif ! if(compute_derivative) then -- end of jacobian calculations

  enddo ! reactions loop

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  if(this%species_id_n2o>0 .and. net_nmin_rate>this%x0eps) then

#ifdef CLM_PFLOTRAN
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
#else
    f_t = 1.0d0
    f_w = 1.0d0
    f_ph = 1.0d0
#endif

    if(f_t > this%x0eps .and. f_w > this%x0eps .and. f_ph > this%x0eps) then
      f_t = min(f_t, 1.0d0)
      f_w = min(f_w, 1.0d0)
      f_ph= min(f_ph, 1.0d0)

      temp_real = f_t * f_w * f_ph
      temp_real = temp_real * this%n2o_frac_mineralization 
    
      ! residuals
      rate_n2o = temp_real * net_nmin_rate * fnh3
 
      Residual(ires_nh3) = Residual(ires_nh3) + rate_n2o 

      Residual(ires_n2o) = Residual(ires_n2o) - 0.5d0 * rate_n2o

      if(this%species_id_ngasmin > 0) then
         Residual(ires_ngasmin) = Residual(ires_ngasmin) - 0.5d0 * rate_n2o
      endif

     !Jacobians
      if (compute_derivative) then
        drate_n2o_dx = temp_real * &
                  (dnet_nmin_rate_dx * fnh3 &
                  + net_nmin_rate * dfnh3_dnh3)

        Jacobian(ires_nh3,ires_nh3) = Jacobian(ires_nh3,ires_nh3) + drate_n2o_dx* &
           rt_auxvar%aqueous%dtotal(this%species_id_nh3,this%species_id_nh3,iphase)
        Jacobian(ires_n2o,ires_nh3) = Jacobian(ires_n2o,ires_nh3) - 0.5d0*drate_n2o_dx* &
           rt_auxvar%aqueous%dtotal(this%species_id_n2o,this%species_id_nh3,iphase)

        if(this%species_id_ngasmin > 0) then
           Jacobian(ires_ngasmin,ires_nh3) = Jacobian(ires_ngasmin,ires_nh3) - &
                                      0.5d0 * drate_n2o_dx
        endif
      endif

    endif  ! end of 'f_t/f_w/f_ph > 1.0d-20'

  endif ! end of 'species_id_n2o > 0'

#ifdef DEBUG
  do ires=1, reaction%ncomp
    temp_real = Residual(ires)

    if (abs(temp_real) > huge(temp_real)) then
      write(option%fid_out, *) 'infinity of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: DECOMPOSITION'
      option%io_buffer = ' checking infinity of Residuals matrix  @ CLMDec_React '
      call printErrMsg(option)
    endif
  enddo
#endif

end subroutine CLMDec_React

! ************************************************************************** !
!
! CLMDecDestroy: Destroys allocatable or pointer objects created in this
!                  module
! author: Guoping Tang
! ************************************************************************** !
subroutine CLMDec_Destroy(this)
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
 
end subroutine CLMDec_Destroy

end module Reaction_Sandbox_CLMDec_class
