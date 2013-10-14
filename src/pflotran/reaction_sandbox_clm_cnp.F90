module Reaction_Sandbox_CLM_CNP_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
!  use Microbial_Aux_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
!#include "definitions.h"
#include "finclude/petscsys.h"

                          ! 14.00674d0 / 12.011d0
  PetscReal, parameter :: CN_ratio_mass_to_mol = 1.16616d0 
                          ! 30.973762 /12.011
  PetscReal, parameter :: CP_ratio_mass_to_mol = 2.578783d0 
  
  type, public :: pool_type
    character(len=MAXWORDLENGTH) :: name_c
    character(len=MAXWORDLENGTH) :: name_n
    character(len=MAXWORDLENGTH) :: name_p
    PetscReal                    :: stoich
    PetscReal                    :: ratio_cn
    PetscReal                    :: ratio_cp
    type(pool_type), pointer     :: next
  end type pool_type

  type, public :: rate_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal                    :: value
    type(rate_type), pointer     :: next
  end type rate_type

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_CLM_CNP_type

! pools

    type(pool_type), pointer  :: Upstream

    PetscInt :: upstream_ispec_c
    PetscInt :: upstream_ispec_n
    PetscInt :: upstream_ispec_p
   
    PetscReal :: upstream_stoich_c
    PetscReal :: upstream_stoich_n
    PetscReal :: upstream_stoich_p

    PetscBool :: bfixed_upstream_cn
    PetscBool :: bfixed_upstream_cp

    PetscInt :: nDownstream
    type(pool_type), pointer  :: Downstream

    PetscInt, pointer :: downstream_ispec_c(:)
    PetscInt, pointer :: downstream_ispec_n(:)
    PetscInt, pointer :: downstream_ispec_p(:)

    PetscReal, pointer :: downstream_stoich_c(:)
    PetscReal, pointer :: downstream_stoich_n(:)
    PetscReal, pointer :: downstream_stoich_p(:)

    PetscReal, pointer :: downstream_cn(:)
    PetscReal, pointer :: downstream_cp(:)

    PetscBool, pointer :: bfixed_downstream_cn(:)
    PetscBool, pointer :: bfixed_downstream_cp(:)

    PetscBool :: bNEnabled
    PetscBool :: bPEnabled

    PetscReal :: NInhibitionCoef
    PetscReal :: PInhibitionCoef

    PetscInt  :: mineral_c_ispec
    PetscInt  :: mineral_n_ispec
    PetscInt  :: mineral_p_ispec

    PetscReal :: mineral_c_stoich  !respiration fraction
    PetscReal :: mineral_n_stoich
    PetscReal :: mineral_p_stoich

! rate terms

    PetscReal :: rate_constant

    PetscInt :: nFirstOrder
    PetscInt :: nMonod
    PetscInt :: nInhibition

    type(rate_type), pointer :: FirstOrder
    type(rate_type), pointer :: Monod
    type(rate_type), pointer :: Inhibition

    PetscInt, pointer :: ispec_1st(:)
    PetscInt, pointer :: ispec_mnd(:)
    PetscInt, pointer :: ispec_inh(:)

    PetscReal, pointer :: half_saturation(:)
    PetscReal, pointer :: inhibition_coef(:)

  contains
    procedure, public :: ReadInput => CLM_CNPRead
    procedure, public :: Setup => CLM_CNPSetup
    procedure, public :: Evaluate => CLM_CNPReact
    procedure, public :: Destroy => CLM_CNPDestroy
  end type reaction_sandbox_CLM_CNP_type

  public :: CLM_CNPCreate

contains

! ************************************************************************** !
!
! CLM_CNPCreate: Allocates CLM_CNP reaction object.
! author: Guoping Tang
! date: 80/15/13
!
! ************************************************************************** !
function CLM_CNPCreate()

  implicit none
  
  class(reaction_sandbox_CLM_CNP_type), pointer :: CLM_CNPCreate

  allocate(CLM_CNPCreate)

  nullify(CLM_CNPCreate%Upstream)

  CLM_CNPCreate%upstream_ispec_c = -1
  CLM_CNPCreate%upstream_ispec_n = -1
  CLM_CNPCreate%upstream_ispec_p = -1
  
  CLM_CNPCreate%upstream_stoich_c = -1.d0
  CLM_CNPCreate%upstream_stoich_n = -1.d0
  CLM_CNPCreate%upstream_stoich_p = -1.d0

  CLM_CNPCreate%bfixed_upstream_cn = PETSC_FALSE
  CLM_CNPCreate%bfixed_upstream_cp = PETSC_FALSE
   
  CLM_CNPCreate%nDownstream = 0
  
nullify(CLM_CNPCreate%Downstream)

  nullify(CLM_CNPCreate%downstream_ispec_c)
  nullify(CLM_CNPCreate%downstream_ispec_n)
  nullify(CLM_CNPCreate%downstream_ispec_p)
  nullify(CLM_CNPCreate%downstream_stoich_c)
  nullify(CLM_CNPCreate%downstream_stoich_n)
  nullify(CLM_CNPCreate%downstream_stoich_p)
  nullify(CLM_CNPCreate%downstream_cn)
  nullify(CLM_CNPCreate%downstream_cp)
  nullify(CLM_CNPCreate%bfixed_downstream_cn)
  nullify(CLM_CNPCreate%bfixed_downstream_cp)

  CLM_CNPCreate%bNEnabled = PETSC_FALSE
  CLM_CNPCreate%bPEnabled = PETSC_FALSE
  CLM_CNPCreate%NInhibitionCoef = -1.0d0
  CLM_CNPCreate%PInhibitionCoef = -1.0d0
  
  CLM_CNPCreate%mineral_c_ispec = -1
  CLM_CNPCreate%mineral_n_ispec = -1
  CLM_CNPCreate%mineral_p_ispec = -1

  CLM_CNPCreate%mineral_c_stoich = -1.0d0
  CLM_CNPCreate%mineral_n_stoich = -1.0d0
  CLM_CNPCreate%mineral_p_stoich = -1.0d0

  CLM_CNPCreate%rate_constant = 0.d0

  CLM_CNPCreate%nFirstOrder = 0
  CLM_CNPCreate%nMonod = 0
  CLM_CNPCreate%nInhibition = 0

  nullify(CLM_CNPCreate%FirstOrder)
  nullify(CLM_CNPCreate%Monod)
  nullify(CLM_CNPCreate%Inhibition)
 
  nullify(CLM_CNPCreate%ispec_1st)
  nullify(CLM_CNPCreate%ispec_mnd)
  nullify(CLM_CNPCreate%ispec_inh)

  nullify(CLM_CNPCreate%half_saturation)
  nullify(CLM_CNPCreate%inhibition_coef)

  nullify(CLM_CNPCreate%next)  

end function CLM_CNPCreate

function PoolCreate()
  implicit none
  type(pool_type), pointer :: PoolCreate
  allocate(PoolCreate)  
  PoolCreate%name_c = ''
  PoolCreate%name_n = ''
  PoolCreate%name_p = ''
  PoolCreate%stoich = -1.d0
  PoolCreate%ratio_cn = -1.d0
  PoolCreate%ratio_cp = -1.d0
  nullify(PoolCreate%next)
end function PoolCreate

recursive subroutine PoolDestroy(pool)
  implicit none
  type(pool_type), pointer :: pool
  if (.not.associated(pool)) return
  call PoolDestroy(pool%next)
  deallocate(pool)
  nullify(pool)
end subroutine PoolDestroy

function RateCreate()
  implicit none
  type(rate_type), pointer :: RateCreate
  allocate(RateCreate)  
  RateCreate%name = ''
  RateCreate%value = -1.0d0
  nullify(RateCreate%next)
end function RateCreate

recursive subroutine RateDestroy(rate)
  implicit none
  type(rate_type), pointer :: rate
  if (.not.associated(rate)) return
  call RateDestroy(rate%next)
  deallocate(rate)
  nullify(rate)
end subroutine RateDestroy

! ************************************************************************** !
!
! CLM_CNPRead:
! author: Guoping Tang
! date: 08/15/13
!
! ************************************************************************** !
subroutine CLM_CNPRead(this,input,option)

  use Option_module
  use String_module
  use Input_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_CLM_CNP_type) :: this
  type(input_type)                     :: input
  type(option_type)                    :: option

  type(pool_type), pointer :: pool
  type(pool_type), pointer :: downstream_pool_prev

  type(rate_type), pointer :: firstorder,firstorder_prev
  type(rate_type), pointer :: monod, monod_prev
  type(rate_type), pointer :: inhibition, inhibition_prev

  PetscReal :: tmp_real, rate_constant, turnover_time

  character(len=MAXWORDLENGTH) :: word

  nullify(downstream_pool_prev)
  nullify(firstorder_prev)
  nullify(monod_prev)
  nullify(inhibition_prev)

  turnover_time = 0.d0
  rate_constant = 0.d0

  do 
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CLM_CNP')
    call StringToUpper(word)   

    select case(trim(word))
      case('UPSTREAM')
        if (.not.associated(this%Upstream)) then
            pool => PoolCreate()
            this%Upstream => pool
            nullify(pool)
        endif

        do
         call InputReadFlotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit

         call InputReadWord(input,option,word,PETSC_TRUE)
         call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CLM_CNP,UPSTREAM')
         call StringToUpper(word)   

            select case(trim(word))
              case('CPOOL')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'name', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,UPSTREAM,CPOOL')
                  this%Upstream%name_c = word
                  this%Upstream%stoich = 1.0
              case('NPOOL')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'name', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,UPSTREAM,NPOOL')
                  this%Upstream%name_n = word
              case('PPOOL')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'name', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,UPSTREAM,PPOOL')
                 this%Upstream%name_p = word
              case('CNRATIO')
                  call InputReadDouble(input,option,tmp_real)  
                  call InputErrorMsg(input,option,'CN ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,UPSTREAMCPOOL')
                  this%Upstream%ratio_cn = tmp_real*CN_ratio_mass_to_mol
              case('CPRATIO')
                  call InputReadDouble(input,option,tmp_real)  
                  call InputErrorMsg(input,option,'CP ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,UPSTREAMCPOOL')
                  this%Upstream%ratio_cp = tmp_real*CP_ratio_mass_to_mol
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CNP,UPSTREAM keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
        enddo

      case('DOWNSTREAM')
        pool => PoolCreate()

        do
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CLM_CNP,DOWNSTREAM')
          call StringToUpper(word)   

          select case(trim(word))
            case('CPOOL')
               call InputReadWord(input,option,word,PETSC_TRUE)
               call InputErrorMsg(input,option,'name', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,DOWNSTREAM,CPOOL')
               call InputReadDouble(input,option,tmp_real)  
               call InputErrorMsg(input,option,'stoichoimetric coefficient', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,DOWNSTREAMCPOOL')
               pool%name_c = word
               pool%stoich = tmp_real
            case('NPOOL')
               call InputReadWord(input,option,word,PETSC_TRUE)
               call InputErrorMsg(input,option,'name', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,DOWNSTREAM,NPOOL')
               pool%name_n = word
            case('PPOOL')
               call InputReadWord(input,option,word,PETSC_TRUE)
               call InputErrorMsg(input,option,'name', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,DOWNSTREAM,PPOOL')
               pool%name_p = word
            case('CNRATIO')
               call InputReadDouble(input,option,tmp_real)  
               call InputErrorMsg(input,option,'CN ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,DOWNSTREAMCPOOL')
               pool%ratio_cn = tmp_real*CN_ratio_mass_to_mol
            case('CPRATIO')
               call InputReadDouble(input,option,tmp_real)  
               call InputErrorMsg(input,option,'CP ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,DOWNSTREAMCPOOL')
               pool%ratio_cp = tmp_real*CP_ratio_mass_to_mol
            case default
               option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CNP,DOWNSTREAM keyword: ' // &
                                  trim(word) // ' not recognized.'
               call printErrMsg(option)
          end select
        enddo

        if (.not.associated(this%Downstream)) then
           this%Downstream => pool
        else
           downstream_pool_prev%next => pool
        endif

        downstream_pool_prev => pool
        nullify(pool)
        this%nDownstream = this%nDownstream + 1

      case('RATE_CONSTANT')
        call InputReadDouble(input,option,rate_constant)
        call InputErrorMsg(input,option,'rate constant', &
               'CHEMISTRY,REACTION_SANDBOX,CLM-CNP,RATE_CONSTANT')
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) then
            input%err_buf = 'CLM-CNP RATE CONSTANT UNITS'
            call InputDefaultMsg(input,option)
        else              
            this%rate_constant = rate_constant * &
            UnitsConvertToInternal(word,option)
        endif
      case('TURNOVER_TIME')
        call InputReadDouble(input,option,turnover_time)
        call InputErrorMsg(input,option,'turnover time', &
               'CHEMISTRY,REACTION_SANDBOX,CLM-CNP,REACTION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) then
            input%err_buf = 'CLM-CN TURNOVER TIME UNITS'
            call InputDefaultMsg(input,option)
        else              
            turnover_time = turnover_time * &
                  UnitsConvertToInternal(word,option)
            this%rate_constant = 1.d0 / turnover_time 
        endif
      case('FIRSTORDER')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'firstorder rate name', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,FIRST ORDER')
        firstorder => RateCreate()
        firstorder%name = word
        if (.not.associated(this%FirstOrder)) then
          this%FirstOrder => firstorder
        else
          firstorder_prev%next => firstorder
        endif
        firstorder_prev => firstorder
        nullify(firstorder)
        this%nFirstOrder = this%nFirstOrder + 1

      case('MONOD')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'species name', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,MONOD')
        call InputReadDouble(input,option,tmp_real)  
        call InputErrorMsg(input,option,'half saturation constant', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,MONOD')
        monod => RateCreate()
        monod%name = word
        monod%value = tmp_real
        if (.not.associated(this%Monod)) then
          this%Monod => monod
        else
          monod_prev%next => monod
        endif
        monod_prev => monod
        nullify(monod)
        this%nMonod = this%nMonod + 1

      case('INHIBITION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'species name', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,INHIBITION')
        call InputReadDouble(input,option,tmp_real)  
        call InputErrorMsg(input,option,'inhibition constant', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,INHIBITION')
        inhibition => RateCreate()
        inhibition%name = word
        inhibition%value = tmp_real

        if (.not.associated(this%inhibition)) then
          this%inhibition => inhibition
        else
          inhibition_prev%next => inhibition
        endif
        inhibition_prev => inhibition
        nullify(inhibition)
        this%nInhibition = this%nInhibition + 1

      case('NINHIBITION')
        call InputReadDouble(input,option,tmp_real)  
        call InputErrorMsg(input,option,'N inhibition constant', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,NINHIBITION')
        this%NInhibitionCoef = tmp_real

      case('PINHIBITION')
        call InputReadDouble(input,option,tmp_real)  
        call InputErrorMsg(input,option,'P inhibition constant', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CNP,NINHIBITION')
        this%PInhibitionCoef = tmp_real

      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CNP keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo

  ! check to ensure that one of turnover time or rate constant is set.
  if (turnover_time > 0.d0 .and. rate_constant > 0.d0) then
       option%io_buffer = 'Only TURNOVER_TIME or RATE_CONSTANT may ' // &
        'be included in a CLM-CNP reaction definition, but not both. ' // &
        'See reaction with upstream pool "' // &
            trim(this%Upstream%name_c) // '".'
          call printErrMsg(option)
  endif  
end subroutine CLM_CNPRead

! ************************************************************************** !
!
! CLM_CNPSetup: 
! author: Guoping Tang
! date: 08/15/13
!
! ************************************************************************** !
subroutine CLM_CNPSetup(this,reaction,option)

  use Option_module
  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Immobile_Aux_module

  implicit none

  class(reaction_sandbox_CLM_CNP_type) :: this
  type(reaction_type)                  :: reaction
  type(option_type)                    :: option

  type(pool_type), pointer             :: cur_pool
  type(rate_type), pointer             :: cur_rate

  PetscInt :: i, icount
  PetscReal :: sum_stoich_c_prod

  character(len=MAXWORDLENGTH) :: word

  word = 'C'
  this%mineral_c_ispec = GetImmobileSpeciesIDFromName( &
         word,reaction%immobile,PETSC_FALSE,option)  

! upstream
  if(.not.associated(this%Upstream)) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CNP check: ' // &
       ' upstream pool is not specified.'
     call printErrMsg(option)
  else
     this%upstream_ispec_c = GetImmobileSpeciesIDFromName( &
         this%Upstream%name_c,reaction%immobile,PETSC_FALSE,option)  
     if(trim(this%Upstream%name_n) .ne. '') then
       this%upstream_ispec_n = GetImmobileSpeciesIDFromName( &
         this%Upstream%name_n,reaction%immobile,PETSC_FALSE,option)  
     else
         if(this%Upstream%ratio_cn .GT. 1.0d-10) then
            this%bfixed_upstream_cn = PETSC_TRUE
!           option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CNP check:' // &
!            ' upstream CN ratio is zero.'
!          call printErrMsg(option)
!         else
            this%upstream_stoich_n = 1.0d0/this%Upstream%ratio_cn
         endif 
     endif
     if(trim(this%Upstream%name_p) .ne. '') then
       this%upstream_ispec_p = GetImmobileSpeciesIDFromName( &
         this%Upstream%name_p,reaction%immobile,PETSC_FALSE,option)  
     else
         if(this%Upstream%ratio_cp .GT. 1.0d-10) then
            this%bfixed_upstream_cp = PETSC_TRUE
            this%upstream_stoich_p = 1.0d0/this%Upstream%ratio_cp
         endif
     endif
  endif   

! downstream
  sum_stoich_c_prod = 0.0d0
  if(this%nDownstream .GE. 1) then
     allocate(this%downstream_ispec_c(this%nDownstream))
     allocate(this%downstream_ispec_n(this%nDownstream))
     allocate(this%downstream_ispec_p(this%nDownstream))

     allocate(this%downstream_stoich_c(this%nDownstream))
     allocate(this%downstream_stoich_n(this%nDownstream))
     allocate(this%downstream_stoich_p(this%nDownstream))

     allocate(this%downstream_cn(this%nDownstream))
     allocate(this%downstream_cp(this%nDownstream))

     allocate(this%bfixed_downstream_cn(this%nDownstream))
     allocate(this%bfixed_downstream_cp(this%nDownstream))

     cur_pool => this%Downstream
     icount = 1
     do
        if(.not.associated(cur_pool)) exit
        this%downstream_ispec_c(icount) = GetImmobileSpeciesIDFromName( &
                    cur_pool%name_c,reaction%immobile,PETSC_FALSE,option)
        this%downstream_stoich_c(icount) = cur_pool%stoich
        sum_stoich_c_prod = sum_stoich_c_prod + cur_pool%stoich
        if(trim(cur_pool%name_n) .ne. '') then
           this%downstream_ispec_n(icount) = GetImmobileSpeciesIDFromName( &
                    cur_pool%name_n,reaction%immobile,PETSC_FALSE,option)  
        else
           this%downstream_cn(icount) = cur_pool%ratio_cn
           if(this%downstream_cn(icount) .GT. 1.0d-10) then
!             option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CNP check:' // &
!               ' downstream CN ratio is zero.'
!             call printErrMsg(option)
!           else
              this%downstream_stoich_n(icount) = 1.0d0/this%downstream_cn(icount)
              this%bfixed_downstream_cn(icount) = PETSC_TRUE
           endif 
        endif
        if(trim(cur_pool%name_p) .ne. '') then
           this%downstream_ispec_p(icount) = GetImmobileSpeciesIDFromName( &
                    cur_pool%name_p,reaction%immobile,PETSC_FALSE,option)  
        else
           this%downstream_cp(icount) = cur_pool%ratio_cp
           if(this%downstream_cp(icount) .GT. 1.0d-10) then
              this%bfixed_downstream_cp(icount) = PETSC_TRUE
              this%downstream_stoich_p(icount) = 1.0d0/this%downstream_cp(icount)
           endif
        endif

        cur_pool => cur_pool%next
        icount = icount + 1
     enddo 
     if(sum_stoich_c_prod .GT. 1.0d0) then
         option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CNP check:' // &
            ' downstream C pools stoichoimetric coefficient sum > 1.'
         call printErrMsg(option)
     else
         this%mineral_c_stoich = 1.0d0 - sum_stoich_c_prod
     endif
  else
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CNP check:' // &
            ' downstream pool not specified.'
     call printErrMsg(option)
  endif   

! check if N is included
  if((this%Upstream%ratio_cn .LT. 0) .and. (this%upstream_ispec_n .LT. 0)) then
      write(option%fid_out,*) 'Neither CN ratio nor valid upstream N pool is specified. N is not considered.'
  else
     do i = 1, this%nDownstream
        if((this%downstream_cn(i) .LT. 0.0d0) .and. (this%downstream_ispec_n(i) .LT. 0)) then
          write(option%fid_out,*) 'Neither CN ratio nor valid N pool is specified for downstream. N is not considered.'
          exit
        endif
     enddo
     this%bNEnabled = PETSC_TRUE  
     word = 'N'
     this%mineral_n_ispec = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)  
  endif

! check if P is included
  if((this%Upstream%ratio_cp .LT. 0) .and. (this%upstream_ispec_p .LT. 0)) then
      write(option%fid_out,*) 'Neither CP ratio nor valid upstream P pool is specified. P is not considered.'
  else
     do i = 1, this%nDownstream
        if((this%downstream_cp(i) .LT. 0.0d0) .and. (this%downstream_ispec_p(i) .LT. 0)) then
          write(option%fid_out,*) 'Neither CP ratio nor valid P pool is specified for downstream. P is not considered.'
          exit
        endif
     enddo
     this%bPEnabled = PETSC_TRUE  
     word = 'P'
     this%mineral_p_ispec = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)  
  endif

! first order rate terms
  if(this%nFirstOrder .GE. 1) then
     allocate(this%ispec_1st(this%nFirstOrder))

     cur_rate => this%FirstOrder
     icount = 1
     do
        if(.not.associated(cur_rate)) exit
        this%ispec_1st(icount) = GetImmobileSpeciesIDFromName( &
                    cur_rate%name,reaction%immobile,PETSC_FALSE,option)
        cur_rate => cur_rate%next
        icount = icount + 1
     enddo 
  endif   

! monod rate terms
  if(this%nMonod .GE. 1) then
     allocate(this%ispec_mnd(this%nMonod))
     allocate(this%half_saturation(this%nMonod))

     cur_rate => this%Monod
     icount = 1
     do
        if(.not.associated(cur_rate)) exit
        this%ispec_mnd(icount) = GetImmobileSpeciesIDFromName( &
                    cur_rate%name,reaction%immobile,PETSC_FALSE,option)
        this%half_saturation(icount) = cur_rate%value
        cur_rate => cur_rate%next
        icount = icount + 1
     enddo 
  endif   

! inhibition rate terms
  if(this%nInhibition .GE. 1) then
     allocate(this%ispec_inh(this%nInhibition))
     allocate(this%inhibition_coef(this%nInhibition))

     cur_rate => this%Inhibition
     icount = 1
     do
        if(.not.associated(cur_rate)) exit
        this%ispec_inh(icount) = GetImmobileSpeciesIDFromName( &
                    cur_rate%name,reaction%immobile,PETSC_FALSE,option)
        this%inhibition_coef(icount) = cur_rate%value
        cur_rate => cur_rate%next
        icount = icount + 1
     enddo 
  endif   

end subroutine CLM_CNPSetup

! ************************************************************************** !
!
! CLM_CNPReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 08/16/13
!
! ************************************************************************** !
subroutine CLM_CNPReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,porosity,volume,reaction, &
                         option,local_id)

  use Option_module
  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  
#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data !, only : rate_plantnuptake_pf 
#endif
  
  implicit none
  
#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif
  
  class(reaction_sandbox_CLM_CNP_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  PetscInt, parameter :: iphase = 1
  PetscInt :: offset 
  PetscInt :: i, j, ires, ires_j, ires_c, ires_n, ires_p
  PetscInt :: ires_uc, ires_un, ires_up, ires_dc, ires_dn, ires_dp
  PetscReal :: conc, c_up, c_N
  PetscReal :: L_water
  PetscReal :: rate, drate, drate_uc 
  PetscReal :: tmp_real
!  PetscReal :: stoich_mineral_n_sign 

  PetscReal :: F_t
  PetscReal :: F_theta
  PetscReal :: temp_K, tc
  PetscReal, parameter :: eps = 1.0d-150
  PetscReal, parameter :: one_over_71_02 = 1.408054069d-2
  PetscReal, parameter :: theta_min = 0.01d0     ! 1/nat log(0.01d0)
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1
  PetscReal, parameter :: twelve_over_14 = 0.857142857143d0
  PetscInt :: local_id
  PetscReal :: maxpsi, psi
  PetscReal, parameter :: minpsi = -10.0d0  ! MPa
  PetscScalar, pointer :: sucsat_pf_loc(:)   !
  PetscScalar, pointer :: soilpsi_pf_loc(:)   !
  PetscErrorCode :: ierr

! CLM4.5 temperature response function
  tc = global_auxvar%temp(1)
  F_t = 1.5d0 ** ((tc - 25.0d0) / 10.0d0)
 
! CLM-CN temperature response function
  ! inhibition due to temperature
  ! Equation: F_t = exp(308.56*(1/17.02 - 1/(T - 227.13)))

!  temp_K = global_auxvar%temp(1) + 273.15d0

!  if(temp_K .GT. 227.15d0) then
!     F_t = exp(308.56d0*(one_over_71_02 - 1.d0/(temp_K - 227.13d0)))
!  else
!     F_t = 0.0d0
!  endif
  ! inhibition due to moisture content
  ! Equation: F_theta = log(theta_min/theta) / log(theta_min/theta_max)
  ! Assumptions: theta is saturation
  !              theta_min = 0.01, theta_max = 1.

#ifdef CLM_PFLOTRAN
  call VecGetArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)
  call VecGetArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi_pf_loc, ierr)

  maxpsi = sucsat_pf_loc(local_id) * (-9.8d-6)
  psi = min(soilpsi_pf_loc(local_id), maxpsi)
  F_theta = log(minpsi/psi)/log(minpsi/maxpsi)

  call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi_pf_loc, ierr)
#else
  F_theta = 1.0d0 
#endif

!  F_theta = log(theta_min/max(theta_min,global_auxvar%sat(1))) &  ! theta could be zero, which will cause math issue
!          * one_over_log_theta_min

  offset = reaction%offset_immobile

  rate = this%rate_constant * volume * F_t * F_theta

  drate_uc = rate

  c_up = rt_auxvar%immobile(this%upstream_ispec_c) 

!  if(c_up .LE. 0.0d0) then
!     rate = 0.0d0
!     return 
!  endif
       
! first order term
  do i = 1, this%nFirstOrder
     conc = rt_auxvar%immobile(this%ispec_1st(i))  
     rate = rate * conc 
     if(this%ispec_1st(i) .ne. this%upstream_ispec_c) then
       drate_uc = drate_uc * conc 
     endif
  enddo 

! monod term
  do i = 1, this%nMonod
     conc = rt_auxvar%immobile(this%ispec_mnd(i))  
     rate = rate * conc/(conc + this%half_saturation(i))  
  enddo 

! inhibition term
  do i = 1, this%nInhibition
     conc = rt_auxvar%immobile(this%ispec_inh(i))  
     rate = rate * this%inhibition_coef(i)/(conc + this%inhibition_coef(i))  
  enddo 

! N limiting
  if((this%bNEnabled)) then
! upstream CN
    if(.not.this%bfixed_upstream_cn) then
       if(rt_auxvar%immobile(this%upstream_ispec_c) .LT. eps) then
          write(option%fid_out,*) 'Upstream C concentration is 0 in CN ratio calculation.'
       endif 
       this%upstream_stoich_n = rt_auxvar%immobile(this%upstream_ispec_n) & 
                  /rt_auxvar%immobile(this%upstream_ispec_c)
    endif

! downstream CN
    this%mineral_n_stoich = this%upstream_stoich_n   ! start from upstream N, substract downstream N
    do i = 1, this%nDownstream
       if(.not.this%bfixed_downstream_cn(i)) then
          if(rt_auxvar%immobile(this%downstream_ispec_c(i)) .LT. eps) then
             write(option%fid_out,*) 'Downtream C concentration 0 in CN ratio calculation.'
          endif 
          this%downstream_stoich_n(i)=rt_auxvar%immobile(this%downstream_ispec_n(i)) & 
                 / rt_auxvar%immobile(this%downstream_ispec_c(i))
       endif
       this%mineral_n_stoich = this%mineral_n_stoich - this%downstream_stoich_n(i) &
                 * this%downstream_stoich_c(i)
    enddo         

    if(this%mineral_n_stoich .LT. 0.0d0 .and. this%NInhibitionCoef .GT. 0.0d0) then
       conc = rt_auxvar%immobile(this%mineral_n_ispec)  
       rate = rate * conc/(conc + this%NInhibitionCoef)  ! N limiting  
    endif
  endif

! P limiting
  if((this%bPEnabled)) then
! upstream CP
    if(.not.this%bfixed_upstream_cp) then
       if(rt_auxvar%immobile(this%upstream_ispec_c) .LT. eps) then
          write(option%fid_out,*) 'Upstream C concentration is 0 in CP ratio calculation.'
       endif 
       this%upstream_stoich_p = rt_auxvar%immobile(this%upstream_ispec_p) & 
                  / rt_auxvar%immobile(this%upstream_ispec_c)
    endif

! downstream CP
    this%mineral_p_stoich = this%upstream_stoich_p   ! start from upstream P, substract downstream P
    do i = 1, this%nDownstream
       if(.not.this%bfixed_downstream_cp(i)) then
          if(rt_auxvar%immobile(this%downstream_ispec_c(i)) .LT. eps) then
             write(option%fid_out,*) 'Downtream C concentration 0 in CP ratio calculation.'
          endif 
          this%downstream_stoich_p(i)=rt_auxvar%immobile(this%downstream_ispec_p(i)) & 
                 / rt_auxvar%immobile(this%downstream_ispec_c(i))
       endif
       this%mineral_p_stoich = this%mineral_p_stoich - this%downstream_stoich_p(i) &
                 * this%downstream_stoich_c(i)
    enddo         

    if(this%mineral_p_stoich .LT. 0.0d0 .and. this%PInhibitionCoef .GT. 0.0d0) then
       conc = rt_auxvar%immobile(this%mineral_p_ispec)  
       rate = rate * conc/(conc + this%PInhibitionCoef)  ! P limiting  
    endif
  endif

! residuals
! upstream C
  ires = reaction%offset_immobile + this%upstream_ispec_c      
  Residual(ires) = Residual(ires) + rate       

! mineral C
  ires_c = reaction%offset_immobile + this%mineral_c_ispec
  Residual(ires_c) = Residual(ires_c) - rate * this%mineral_c_stoich

! downstream C
  do i = 1, this%nDownstream
     ires = reaction%offset_immobile + this%downstream_ispec_c(i)
     Residual(ires) = Residual(ires) - rate * this%downstream_stoich_c(i)
  enddo 

  if(this%bNEnabled) then
! upstream N
     if(.not.this%bfixed_upstream_cn) then
        ires = reaction%offset_immobile + this%upstream_ispec_n
        Residual(ires) = Residual(ires) + rate * this%upstream_stoich_n      
     endif
! mineral N
     ires_n = reaction%offset_immobile + this%mineral_n_ispec      
     Residual(ires_n) = Residual(ires_n) - rate * this%mineral_n_stoich      
! downstream N
     do i = 1, this%nDownstream
       if(.not.this%bfixed_downstream_cn(i)) then
          ires = reaction%offset_immobile + this%downstream_ispec_n(i)      
          Residual(ires) = Residual(ires) - rate * this%downstream_stoich_n(i)
       endif
     enddo 
  endif

  if(this%bPEnabled) then
! upstream P
     if(.not.this%bfixed_upstream_cp) then
        ires = reaction%offset_immobile + this%upstream_ispec_p      
        Residual(ires) = Residual(ires) + rate * this%upstream_stoich_p      
     endif
! mineral P
     ires_p = reaction%offset_immobile + this%mineral_p_ispec      
     Residual(ires_p) = Residual(ires_p) - rate * this%mineral_p_stoich      
! downstream P
     do i = 1, this%nDownstream
       if(.not.this%bfixed_downstream_cp(i)) then
          ires = reaction%offset_immobile + this%downstream_ispec_p(i)      
          Residual(ires) = Residual(ires) - rate * this%downstream_stoich_p(i)      
       endif
     enddo 
  endif

! Jacobian
  if (.not.compute_derivative) return 

  if(this%bNEnabled) then
     if(.not.this%bfixed_upstream_cn) then
        ires_uc = reaction%offset_immobile + this%upstream_ispec_c
        ires_un = reaction%offset_immobile + this%upstream_ispec_n
        tmp_real = rt_auxvar%immobile(this%upstream_ispec_n) &
                 / (rt_auxvar%immobile(this%upstream_ispec_c)) !&
!                 / (rt_auxvar%immobile(this%upstream_ispec_c))
        Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) - drate_uc * tmp_real !rate * tmp_real
        Jacobian(ires_n, ires_uc) = Jacobian(ires_n, ires_uc) + drate_uc * tmp_real !rate * tmp_real
        Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) + drate_uc !rate / & 
                                    !rt_auxvar%immobile(this%upstream_ispec_c)
        Jacobian(ires_n,ires_un) = Jacobian(ires_n,ires_un) - drate_uc   !rate / & 
                                    !rt_auxvar%immobile(this%upstream_ispec_c)
     endif

     do i = 1, this%nDownstream
        if(.not.this%bfixed_downstream_cn(i)) then
           ires_dc = reaction%offset_immobile + this%downstream_ispec_c(i)
           ires_dn = reaction%offset_immobile + this%downstream_ispec_n(i)
           tmp_real = rt_auxvar%immobile(this%downstream_ispec_n(i)) &
                    / rt_auxvar%immobile(this%downstream_ispec_c(i)) &
                    / rt_auxvar%immobile(this%downstream_ispec_c(i))
           Jacobian(ires_dn,ires_dc) = Jacobian(ires_dn,ires_dc) + rate * tmp_real
           Jacobian(ires_n, ires_dc) = Jacobian(ires_n, ires_dc) - rate * tmp_real
           Jacobian(ires_dn,ires_dn) = Jacobian(ires_dn,ires_dn) - rate / & 
                                    rt_auxvar%immobile(this%downstream_ispec_c(i))
           Jacobian(ires_n, ires_dn) = Jacobian(ires_n, ires_dn) + rate / &
                                    rt_auxvar%immobile(this%downstream_ispec_c(i))
        endif
     enddo 
  endif

  if(this%bPEnabled) then
     if(.not.this%bfixed_upstream_cp) then
        ires_uc = reaction%offset_immobile + this%upstream_ispec_c
        ires_up = reaction%offset_immobile + this%upstream_ispec_p
        tmp_real = rt_auxvar%immobile(this%upstream_ispec_p) &
                 / rt_auxvar%immobile(this%upstream_ispec_c) !&
!                 / rt_auxvar%immobile(this%upstream_ispec_c)
        Jacobian(ires_up,ires_uc) = Jacobian(ires_up,ires_uc) - drate_uc * tmp_real !rate * tmp_real
        Jacobian(ires_p, ires_uc) = Jacobian(ires_p, ires_uc) + drate_uc * tmp_real !rate * tmp_real
        Jacobian(ires_up,ires_up) = Jacobian(ires_up,ires_up) + drate_uc !rate / & 
!                                    rt_auxvar%immobile(this%upstream_ispec_c)
        Jacobian(ires_p, ires_up) = Jacobian(ires_p, ires_up) - drate_uc !rate / &
!                                    rt_auxvar%immobile(this%upstream_ispec_c)
     endif

     do i = 1, this%nDownstream
        if(.not.this%bfixed_downstream_cp(i)) then
           ires_dc = reaction%offset_immobile + this%downstream_ispec_c(i)
           ires_dp = reaction%offset_immobile + this%downstream_ispec_p(i)
           tmp_real = rt_auxvar%immobile(this%downstream_ispec_p(i)) &
                    / rt_auxvar%immobile(this%downstream_ispec_c(i)) &
                    / rt_auxvar%immobile(this%downstream_ispec_c(i))
           Jacobian(ires_dp,ires_dc) = Jacobian(ires_dp,ires_dc) + rate * tmp_real
           Jacobian(ires_p, ires_dc) = Jacobian(ires_p, ires_dc) - rate * tmp_real
           Jacobian(ires_dp,ires_dp) = Jacobian(ires_dp,ires_dp) - rate / & 
                                    rt_auxvar%immobile(this%downstream_ispec_c(i))
           Jacobian(ires_p, ires_dp) = Jacobian(ires_p, ires_dp) + rate / &
                                    rt_auxvar%immobile(this%downstream_ispec_c(i))
        endif
     enddo 
  endif

!  first order terms
  do j = 1, this%nFirstOrder
     conc = rt_auxvar%immobile(this%ispec_1st(j))  
     if(this%ispec_1st(i) .eq. this%upstream_ispec_c) then
       drate = drate_uc 
     else
       drate = rate / conc 
     endif

     ires_j = reaction%offset_immobile + this%ispec_1st(j) 

! upstream C
     ires = reaction%offset_immobile + this%upstream_ispec_c
     Jacobian(ires,ires_j) = Jacobian(ires,ires_j) + drate
       
! mineral C
     Jacobian(ires_c, ires_j) = Jacobian(ires_c, ires_j) - drate * this%mineral_c_stoich

! downstream C
     do i = 1, this%nDownstream
        ires = reaction%offset_immobile + this%downstream_ispec_c(i)
        Jacobian(ires,ires_j) = Jacobian(ires,ires_j) - drate * this%downstream_stoich_c(i)
     enddo 

     if(this%bNEnabled) then
! upstream N
        if(.not.this%bfixed_upstream_cn) then
           ires = reaction%offset_immobile + this%upstream_ispec_n      
           Jacobian(ires,ires_j) = Jacobian(ires,ires_j) + drate * this%upstream_stoich_n      
        endif 
! mineral N
        Jacobian(ires_n,ires_j) = Jacobian(ires_n,ires_j) - drate * this%mineral_n_stoich
! downstream N      
        do i = 1, this%nDownstream
           if(.not.this%bfixed_downstream_cn(i)) then
             ires = reaction%offset_immobile + this%downstream_ispec_n(i)      
             Jacobian(ires,ires_j) = Jacobian(ires,ires_j) - drate * this%downstream_stoich_n(i)
           endif
        enddo 
     endif

     if(this%bPEnabled) then
! upstream P
        if(.not.this%bfixed_upstream_cp) then
           ires = reaction%offset_immobile + this%upstream_ispec_p      
           Jacobian(ires,ires_j) = Jacobian(ires,ires_j) + drate * this%upstream_stoich_p      
        endif 
! mineral P
        Jacobian(ires_p,ires_j) = Jacobian(ires_p,ires_j) - drate * this%mineral_p_stoich
! downstream P      
        do i = 1, this%nDownstream
           if(.not.this%bfixed_downstream_cp(i)) then
             ires = reaction%offset_immobile + this%downstream_ispec_p(i)      
             Jacobian(ires,ires_j) = Jacobian(ires,ires_j) - drate * this%downstream_stoich_p(i)
           endif
        enddo 
     endif
  enddo 

!  monod terms
  do j = 1, this%nMonod
! r     = r0 * x / (k + x)     
! dr/dx = r0 * k / (k + x) / (k + x) = r * k / x / (k + x)
     conc = rt_auxvar%immobile(this%ispec_mnd(j))  
     drate = rate * this%half_saturation(j) / conc / (this%half_saturation(j) + conc) 
     ires_j = reaction%offset_immobile + this%ispec_mnd(j) 

! upstream C
     ires = reaction%offset_immobile + this%upstream_ispec_c      
     Jacobian(ires, ires_j) = Jacobian(ires, ires_j) + drate       

! mineral C
     Jacobian(ires_c, ires_j) = Jacobian(ires, ires_j) - drate * this%mineral_c_stoich

! downstream C
     do i = 1, this%nDownstream
        ires = reaction%offset_immobile + this%downstream_ispec_c(i)      
        Jacobian(ires,ires_j) = Jacobian(ires,ires_j) - drate * this%downstream_stoich_c(i)
     enddo 

     if(this%bNEnabled) then
! upstream N
       if(.not.this%bfixed_upstream_cn) then
          ires = reaction%offset_immobile + this%upstream_ispec_n      
          Jacobian(ires,ires_j) = Jacobian(ires,ires_j) + drate * this%upstream_stoich_n
       endif
 ! mineral N
       Jacobian(ires_n,ires_j) = Jacobian(ires_n,ires_j) - drate * this%mineral_n_stoich      

! downstream N       
       do i = 1, this%nDownstream
        if(.not.this%bfixed_downstream_cn(i)) then
          ires = reaction%offset_immobile + this%downstream_ispec_n(i)      
          Jacobian(ires,ires_j) = Jacobian(ires,ires_j) - drate * this%downstream_stoich_n(i)      
        endif
       enddo 
     endif

     if(this%bPEnabled) then
! upstream P
       if(.not.this%bfixed_upstream_cp) then
          ires = reaction%offset_immobile + this%upstream_ispec_p      
          Jacobian(ires,ires_j) = Jacobian(ires,ires_j) + drate * this%upstream_stoich_p
       endif
 ! mineral P
       Jacobian(ires_p,ires_j) = Jacobian(ires_p,ires_j) - drate * this%mineral_p_stoich      

! downstream P       
       do i = 1, this%nDownstream
        if(.not.this%bfixed_downstream_cp(i)) then
          ires = reaction%offset_immobile + this%downstream_ispec_p(i)      
          Jacobian(ires,ires_j) = Jacobian(ires,ires_j) - drate * this%downstream_stoich_p(i)      
        endif
       enddo 
     endif
   enddo 

!  inhibition terms
   do j = 1, this%nInhibition
! r     = r0 * k / (k + x)     
! dr/dx = - r0 * k / (k + x) / (k + x) = - r / (k + x)
     conc = rt_auxvar%immobile(this%ispec_inh(j))  
     drate = -rate / (this%inhibition_coef(j) + conc) 
     ires_j = reaction%offset_immobile + this%ispec_inh(j) 

! upstream C
     ires = reaction%offset_immobile + this%upstream_ispec_c      
     Jacobian(ires, ires_j) = Jacobian(ires, ires_j) + drate       

! mineral C
     ires_c = reaction%offset_immobile + this%mineral_c_ispec      
     Jacobian(ires_c, ires_j) = Jacobian(ires_c, ires_j) - drate * this%mineral_c_stoich      

! downstream C
     do i = 1, this%nDownstream
        ires = reaction%offset_immobile + this%downstream_ispec_c(i)      
        Jacobian(ires,ires_j) = Jacobian(ires,ires_j) - drate * this%downstream_stoich_c(i)
     enddo 

     if(this%bNEnabled) then
! upstream N
       if(.not.this%bfixed_upstream_cn) then
          ires = reaction%offset_immobile + this%upstream_ispec_n      
          Jacobian(ires,ires_j) = Jacobian(ires,ires_j) + drate * this%upstream_stoich_n      
       endif

! mineral N 
       ires_n = reaction%offset_immobile + this%mineral_n_ispec      
       Jacobian(ires_n,ires_j) = Jacobian(ires_n,ires_j) - drate * this%mineral_n_stoich      

! downstream N       
       do i = 1, this%nDownstream
        if(.not.this%bfixed_downstream_cn(i)) then
          ires = reaction%offset_immobile + this%downstream_ispec_n(i)      
          Jacobian(ires,ires_j) = Jacobian(ires,ires_j) - drate * this%downstream_stoich_n(i)      
        endif
       enddo 
     endif

     if(this%bPEnabled) then
! upstream P
       if(.not.this%bfixed_upstream_cp) then
          ires = reaction%offset_immobile + this%upstream_ispec_p      
          Jacobian(ires,ires_j) = Jacobian(ires,ires_j) + drate * this%upstream_stoich_p      
       endif

! mineral P 
       ires_p = reaction%offset_immobile + this%mineral_p_ispec      
       Jacobian(ires_p,ires_j) = Jacobian(ires_p,ires_j) - drate * this%mineral_p_stoich      

! downstream P       
       do i = 1, this%nDownstream
        if(.not.this%bfixed_downstream_cp(i)) then
          ires = reaction%offset_immobile + this%downstream_ispec_p(i)      
          Jacobian(ires,ires_j) = Jacobian(ires,ires_j) - drate * this%downstream_stoich_p(i)      
        endif
       enddo 
     endif
   enddo 

! N limiting
   if(this%bNEnabled .and. this%mineral_n_stoich .LT. 0.0d0 .and. this%NInhibitionCoef .GT. 0.0d0) then
! r     = r0 * x / (k + x)     
! dr/dx = r0 * k / (k + x) / (k + x) = r * k / x / (k + x)
     conc = rt_auxvar%immobile(this%mineral_n_ispec)  
     drate = rate * this%NInhibitionCoef / conc / (this%NInhibitionCoef + conc)
     ires_n = reaction%offset_immobile + this%mineral_n_ispec      

! upstream C
     ires = reaction%offset_immobile + this%upstream_ispec_c
     Jacobian(ires, ires_n) = Jacobian(ires, ires_n) + drate

! mineral C
     ires_c = reaction%offset_immobile + this%mineral_c_ispec      
     Jacobian(ires_c, ires_n) = Jacobian(ires_c, ires_n) - drate * this%mineral_c_stoich

! downstream C
     do i = 1, this%nDownstream
        ires = reaction%offset_immobile + this%downstream_ispec_c(i)      
        Jacobian(ires,ires_n) = Jacobian(ires,ires_n) - drate * this%downstream_stoich_c(i)
     enddo 

     if(.not.this%bfixed_upstream_cn) then
! upstream N
       ires = reaction%offset_immobile + this%upstream_ispec_n      
       Jacobian(ires,ires_n) = Jacobian(ires,ires_n) + drate * this%upstream_stoich_n      

! mineral N       
       Jacobian(ires_n,ires_n) = Jacobian(ires_n,ires_n) - drate * this%mineral_n_stoich      
     endif 

! downstream N
     do i = 1, this%nDownstream
        if(.not.this%bfixed_downstream_cn(i)) then
          ires = reaction%offset_immobile + this%downstream_ispec_n(i)      
          Jacobian(ires,ires_n) = Jacobian(ires,ires_n) - drate * this%downstream_stoich_n(i)
        endif
     enddo 
     
     if(this%bPEnabled) then
! upstream P
        if(.not.this%bfixed_upstream_cp) then
          ires = reaction%offset_immobile + this%upstream_ispec_p      
          Jacobian(ires,ires_n) = Jacobian(ires,ires_n) + drate * this%upstream_stoich_p      

! mineral P       
          Jacobian(ires_p,ires_n) = Jacobian(ires_p,ires_n) - drate * this%mineral_p_stoich
        endif 

! downstream P
        do i = 1, this%nDownstream
           if(.not.this%bfixed_downstream_cp(i)) then
               ires = reaction%offset_immobile + this%downstream_ispec_p(i)
               Jacobian(ires,ires_n) = Jacobian(ires,ires_n) - drate * this%downstream_stoich_p(i)
           endif
        enddo 
     endif
   endif  ! if N limiting

! P limiting
   if(this%bPEnabled .and. this%mineral_p_stoich .LT. 0.0d0 .and. this%PInhibitionCoef .GT. 0.0d0) then
! r     = r0 * x / (k + x)     
! dr/dx = r0 * k / (k + x) / (k + x) = r * k / x / (k + x)
     conc = rt_auxvar%immobile(this%mineral_p_ispec)  
     drate = rate * this%PInhibitionCoef / conc / (this%PInhibitionCoef + conc)
     ires_p = reaction%offset_immobile + this%mineral_p_ispec      

! upstream C
     ires = reaction%offset_immobile + this%upstream_ispec_c
     Jacobian(ires, ires_p) = Jacobian(ires, ires_p) + drate

! mineral C
     ires_c = reaction%offset_immobile + this%mineral_c_ispec      
     Jacobian(ires_c, ires_p) = Jacobian(ires_c, ires_p) - drate * this%mineral_c_stoich

! downstream C
     do i = 1, this%nDownstream
        ires = reaction%offset_immobile + this%downstream_ispec_c(i)      
        Jacobian(ires,ires_p) = Jacobian(ires,ires_p) - drate * this%downstream_stoich_c(i)
     enddo 

     if(.not.this%bfixed_upstream_cp) then
! upstream P
       ires = reaction%offset_immobile + this%upstream_ispec_p      
       Jacobian(ires,ires_p) = Jacobian(ires,ires_p) + drate * this%upstream_stoich_p      

! mineral P       
       Jacobian(ires_p,ires_p) = Jacobian(ires_p,ires_p) - drate * this%mineral_p_stoich      
     endif 

! downstream P
     do i = 1, this%nDownstream
        if(.not.this%bfixed_downstream_cp(i)) then
          ires = reaction%offset_immobile + this%downstream_ispec_p(i)      
          Jacobian(ires,ires_p) = Jacobian(ires,ires_p) - drate * this%downstream_stoich_p(i)
        endif
     enddo 

     if(this%bNEnabled) then
       if(.not.this%bfixed_upstream_cn) then
! upstream N
          ires = reaction%offset_immobile + this%upstream_ispec_n      
          Jacobian(ires,ires_p) = Jacobian(ires,ires_p) + drate * this%upstream_stoich_n      

! mineral N       
          Jacobian(ires_n,ires_p) = Jacobian(ires_n,ires_p) - drate * this%mineral_n_stoich      
        endif 

! downstream N
        do i = 1, this%nDownstream
           if(.not.this%bfixed_downstream_cn(i)) then
             ires = reaction%offset_immobile + this%downstream_ispec_n(i)      
             Jacobian(ires,ires_p) = Jacobian(ires,ires_p) - drate * this%downstream_stoich_n(i)
           endif
        enddo 
     endif
   endif  ! if P limiting
end subroutine CLM_CNPReact

! ************************************************************************** !
!
! CLM_CNPDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Guoping Tang
! date: 00/00/00
!
! ************************************************************************** !
subroutine CLM_CNPDestroy(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_CLM_CNP_type) :: this  

  type(pool_type), pointer :: cur_pool, prev_pool
  type(rate_type), pointer :: cur_rate, prev_rate

  cur_pool => this%Upstream
  deallocate(cur_pool)

  cur_pool => this%Downstream
  do 
    if(.not.associated(cur_pool)) exit
       prev_pool => cur_pool
       cur_pool => cur_pool%next
       deallocate(prev_pool)
       nullify(prev_pool)
  enddo

  cur_rate => this%FirstOrder
  do 
    if(.not.associated(cur_rate)) exit
       prev_rate => cur_rate
       cur_rate => cur_rate%next
       deallocate(prev_rate)
       nullify(prev_rate)
  enddo

  cur_rate => this%Monod
  do 
    if(.not.associated(cur_rate)) exit
       prev_rate => cur_rate
       cur_rate => cur_rate%next
       deallocate(prev_rate)
       nullify(prev_rate)
  enddo

  cur_rate => this%Inhibition
  do 
    if(.not.associated(cur_rate)) exit
       prev_rate => cur_rate
       cur_rate => cur_rate%next
       deallocate(prev_rate)
       nullify(prev_rate)
  enddo

  call DeallocateArray(this%downstream_ispec_c) 
  call DeallocateArray(this%downstream_ispec_n) 
  call DeallocateArray(this%downstream_ispec_p) 

  call DeallocateArray(this%downstream_stoich_c) 
  call DeallocateArray(this%downstream_stoich_n) 
  call DeallocateArray(this%downstream_stoich_p) 

  call DeallocateArray(this%downstream_cn) 
  call DeallocateArray(this%downstream_cp)

  call DeallocateArray(this%bfixed_downstream_cn) 
  call DeallocateArray(this%bfixed_downstream_cp)

  call DeallocateArray(this%ispec_1st) 
  call DeallocateArray(this%ispec_mnd) 
  call DeallocateArray(this%ispec_inh) 

  call DeallocateArray(this%half_saturation) 
  call DeallocateArray(this%inhibition_coef) 
 
end subroutine CLM_CNPDestroy

end module Reaction_Sandbox_CLM_CNP_class
