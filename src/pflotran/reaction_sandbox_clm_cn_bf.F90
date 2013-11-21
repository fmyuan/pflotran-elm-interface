module Reaction_Sandbox_CLM_CN_BF_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
!  use Microbial_Aux_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

                          ! 14.00674d0 / 12.011d0
  PetscReal, parameter :: CN_ratio_mass_to_mol = 1.16616d0 
                          ! 30.973762 /12.011
  PetscReal, parameter :: CP_ratio_mass_to_mol = 2.578783d0 

! added for CLM-CN-MICROBE  11/20/2013 
!-------------------------------------------------------------------------------
  PetscBool, parameter :: CLM_CN_MICROBE = PETSC_TRUE 
  PetscReal, parameter :: CN_ratio_microbe = 9.32928d0   ! 8.0d0 
  PetscReal, parameter :: CN_ratio_bacteria = 5.8038d0   ! 5.0d0 
  PetscReal, parameter :: CN_ratio_fungi = 17.4924d0     !15.0d0 ! or 10.0 
  PetscReal, parameter :: CUE_max = 0.6d0 
  PetscReal, parameter :: fraction_bacteria = 0.340927d0 ! 5.0**0.6/(5.0**0.6+15.0**0.6) 
!-------------------------------------------------------------------------------

  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_CLM4 = 1 
  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_Q10 = 2 

  type, public :: pool_type
    character(len=MAXWORDLENGTH) :: name_c
    character(len=MAXWORDLENGTH) :: name_n
    character(len=MAXWORDLENGTH) :: name_p
    PetscReal                    :: stoich
    PetscReal                    :: ratio_cn
    PetscReal                    :: ratio_cp
    PetscReal                    :: ratio_nc
    PetscReal                    :: ratio_pc
    PetscBool                    :: bratio_cn_follow_upstream
    PetscBool                    :: bratio_cp_follow_upstream
    type(pool_type), pointer     :: next
  end type pool_type

  type, public :: rate_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal                    :: value
    type(rate_type), pointer     :: next
  end type rate_type

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_CLM_CN_BF_type

    PetscInt :: temperature_response_function
    PetscReal :: Q10

! pools

    type(pool_type), pointer  :: Upstream

    PetscInt :: upstream_ispec_c
    PetscInt :: upstream_ispec_n
    PetscInt :: upstream_ispec_p
   
    PetscInt :: upstream_ires_c
    PetscInt :: upstream_ires_n
    PetscInt :: upstream_ires_p
   
    PetscBool :: upstream_isaqueous_c
    PetscBool :: upstream_isaqueous_n
    PetscBool :: upstream_isaqueous_p

    character(len=MAXWORDLENGTH) :: name_mc
    character(len=MAXWORDLENGTH) :: name_mn
    character(len=MAXWORDLENGTH) :: name_mp

    PetscBool :: isaqueous_mc
    PetscBool :: isaqueous_mn
    PetscBool :: isaqueous_mp
   
    PetscReal :: upstream_stoich_c    !cu
    PetscReal :: upstream_stoich_n    !nu
    PetscReal :: upstream_stoich_p    !pu

    PetscBool :: bfixed_upstream_cn
    PetscBool :: bfixed_upstream_cp

! added for CLM-CN-MICROBE  11/20/2013 
!-------------------------------------------------------------------------------
    PetscBool :: is_upstream_litter
!-------------------------------------------------------------------------------

    PetscInt :: nDownstream
    type(pool_type), pointer  :: Downstream

    PetscInt, pointer :: downstream_ispec_c(:)
    PetscInt, pointer :: downstream_ispec_n(:)
    PetscInt, pointer :: downstream_ispec_p(:)

    PetscInt, pointer :: downstream_ires_c(:)
    PetscInt, pointer :: downstream_ires_n(:)
    PetscInt, pointer :: downstream_ires_p(:)

    PetscBool, pointer :: downstream_isaqueous_c(:)
    PetscBool, pointer :: downstream_isaqueous_n(:)
    PetscBool, pointer :: downstream_isaqueous_p(:)

    PetscReal, pointer :: downstream_stoich_c(:) !ci
    PetscReal, pointer :: downstream_stoich_n(:) !ni
    PetscReal, pointer :: downstream_stoich_p(:) !pi

!    PetscReal, pointer :: downstream_cn(:)
!    PetscReal, pointer :: downstream_cp(:)

    PetscBool, pointer :: b_downstream_cn_fixed(:)
    PetscBool, pointer :: b_downstream_cp_fixed(:)

    ! downstream CN or CP ratio follow upstream?
    PetscBool, pointer :: b_downstream_cn_follow_upstream(:)
    PetscBool, pointer :: b_downstream_cp_follow_upstream(:)

    PetscBool :: bNEnabled
    PetscBool :: bPEnabled

    PetscReal :: NInhibitionCoef
    PetscReal :: PInhibitionCoef

    PetscInt  :: ispec_mc
    PetscInt  :: ispec_mn
    PetscInt  :: ispec_mp

    PetscInt  :: ires_mc
    PetscInt  :: ires_mn
    PetscInt  :: ires_mp

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

    PetscBool, pointer :: isaqueous_1st(:)
    PetscBool, pointer :: isaqueous_mnd(:)
    PetscBool, pointer :: isaqueous_inh(:)

    PetscReal, pointer :: half_saturation(:)
    PetscReal, pointer :: inhibition_coef(:)

  contains
    procedure, public :: ReadInput => CLM_CN_BFRead
    procedure, public :: Setup => CLM_CN_BFSetup
    procedure, public :: Evaluate => CLM_CN_BFReact
    procedure, public :: Destroy => CLM_CN_BFDestroy
  end type reaction_sandbox_CLM_CN_BF_type

  public :: CLM_CN_BFCreate

contains

! ************************************************************************** !
!
! CLM_CN_BFCreate: Allocates CLM_CN_BF reaction object.
! author: Guoping Tang
! date: 80/15/13
!
! ************************************************************************** !
function CLM_CN_BFCreate()

  implicit none
  
  class(reaction_sandbox_CLM_CN_BF_type), pointer :: CLM_CN_BFCreate

  allocate(CLM_CN_BFCreate)

  nullify(CLM_CN_BFCreate%Upstream)

  CLM_CN_BFCreate%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_Q10
  CLM_CN_BFCreate%Q10 = 1.5d0

  CLM_CN_BFCreate%upstream_ispec_c = -1
  CLM_CN_BFCreate%upstream_ispec_n = -1
  CLM_CN_BFCreate%upstream_ispec_p = -1
  
  CLM_CN_BFCreate%upstream_ires_c = -1
  CLM_CN_BFCreate%upstream_ires_n = -1
  CLM_CN_BFCreate%upstream_ires_p = -1
  
  CLM_CN_BFCreate%upstream_isaqueous_c = PETSC_FALSE
  CLM_CN_BFCreate%upstream_isaqueous_n = PETSC_FALSE
  CLM_CN_BFCreate%upstream_isaqueous_p = PETSC_FALSE

  CLM_CN_BFCreate%isaqueous_mc = PETSC_FALSE
  CLM_CN_BFCreate%isaqueous_mn = PETSC_FALSE
  CLM_CN_BFCreate%isaqueous_mp = PETSC_FALSE
  
  CLM_CN_BFCreate%upstream_stoich_c = -1.d0
  CLM_CN_BFCreate%upstream_stoich_n = -1.d0
  CLM_CN_BFCreate%upstream_stoich_p = -1.d0

  CLM_CN_BFCreate%bfixed_upstream_cn = PETSC_FALSE
  CLM_CN_BFCreate%bfixed_upstream_cp = PETSC_FALSE
   
! added for CLM-CN-MICROBE  11/20/2013 
!-------------------------------------------------------------------------------
  CLM_CN_BFCreate%is_upstream_litter = PETSC_FALSE
!-------------------------------------------------------------------------------

  CLM_CN_BFCreate%nDownstream = 0
  
  nullify(CLM_CN_BFCreate%Downstream)

  nullify(CLM_CN_BFCreate%downstream_ispec_c)
  nullify(CLM_CN_BFCreate%downstream_ispec_n)
  nullify(CLM_CN_BFCreate%downstream_ispec_p)

  nullify(CLM_CN_BFCreate%downstream_ires_c)
  nullify(CLM_CN_BFCreate%downstream_ires_n)
  nullify(CLM_CN_BFCreate%downstream_ires_p)

  nullify(CLM_CN_BFCreate%downstream_isaqueous_c)
  nullify(CLM_CN_BFCreate%downstream_isaqueous_n)
  nullify(CLM_CN_BFCreate%downstream_isaqueous_p)

  nullify(CLM_CN_BFCreate%downstream_stoich_c)
  nullify(CLM_CN_BFCreate%downstream_stoich_n)
  nullify(CLM_CN_BFCreate%downstream_stoich_p)
!  nullify(CLM_CN_BFCreate%downstream_cn)
!  nullify(CLM_CN_BFCreate%downstream_cp)
  nullify(CLM_CN_BFCreate%b_downstream_cn_fixed)
  nullify(CLM_CN_BFCreate%b_downstream_cp_fixed)
  nullify(CLM_CN_BFCreate%b_downstream_cn_follow_upstream)
  nullify(CLM_CN_BFCreate%b_downstream_cp_follow_upstream)

  CLM_CN_BFCreate%bNEnabled = PETSC_FALSE
  CLM_CN_BFCreate%bPEnabled = PETSC_FALSE
  CLM_CN_BFCreate%NInhibitionCoef = -1.0d0
  CLM_CN_BFCreate%PInhibitionCoef = -1.0d0
  
  CLM_CN_BFCreate%name_mc = 'C'
  CLM_CN_BFCreate%name_mn = 'N'
  CLM_CN_BFCreate%name_mp = 'P'

  CLM_CN_BFCreate%ispec_mc = -1
  CLM_CN_BFCreate%ispec_mn = -1
  CLM_CN_BFCreate%ispec_mp = -1

  CLM_CN_BFCreate%ires_mc = -1
  CLM_CN_BFCreate%ires_mn = -1
  CLM_CN_BFCreate%ires_mp = -1

  CLM_CN_BFCreate%mineral_c_stoich = -1.0d0
  CLM_CN_BFCreate%mineral_n_stoich = -1.0d0
  CLM_CN_BFCreate%mineral_p_stoich = -1.0d0

  CLM_CN_BFCreate%rate_constant = 0.d0

  CLM_CN_BFCreate%nFirstOrder = 0
  CLM_CN_BFCreate%nMonod = 0
  CLM_CN_BFCreate%nInhibition = 0

  nullify(CLM_CN_BFCreate%FirstOrder)
  nullify(CLM_CN_BFCreate%Monod)
  nullify(CLM_CN_BFCreate%Inhibition)
 
  nullify(CLM_CN_BFCreate%ispec_1st)
  nullify(CLM_CN_BFCreate%ispec_mnd)
  nullify(CLM_CN_BFCreate%ispec_inh)

  nullify(CLM_CN_BFCreate%isaqueous_1st)
  nullify(CLM_CN_BFCreate%isaqueous_mnd)
  nullify(CLM_CN_BFCreate%isaqueous_inh)

  nullify(CLM_CN_BFCreate%half_saturation)
  nullify(CLM_CN_BFCreate%inhibition_coef)

  nullify(CLM_CN_BFCreate%next)  

end function CLM_CN_BFCreate

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
  PoolCreate%ratio_nc = -1.d0
  PoolCreate%ratio_pc = -1.d0
  PoolCreate%bratio_cn_follow_upstream = PETSC_FALSE
  PoolCreate%bratio_cp_follow_upstream = PETSC_FALSE
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
! CLM_CN_BFRead:
! author: Guoping Tang
! date: 08/15/13
!
! ************************************************************************** !
subroutine CLM_CN_BFRead(this,input,option)

  use Option_module
  use String_module
  use Input_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_CLM_CN_BF_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF')
    call StringToUpper(word)   

    select case(trim(word))
      case('TEMPERATURE_RESPONSE_FUNCTION')
        do
         call InputReadFlotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit

         call InputReadWord(input,option,word,PETSC_TRUE)
         call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF,TEMPERATURE RESPONSE FUNCTION')
         call StringToUpper(word)   

            select case(trim(word))
              case('CLM4')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4    
              case('Q10')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_Q10    
                  call InputReadDouble(input,option,tmp_real)  
                  call InputErrorMsg(input,option,'Q10', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,TEMPERATURE RESPONSE FUNCTION')
                  this%Q10 = tmp_real
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF,TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
         enddo 

      case('MINERALC')
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'mineral c', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,UPSTREAM,MINERALC')
          this%name_mc = word

      case('MINERALN')
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'mineral n', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,UPSTREAM,MINERALN')
          this%name_mn = word

      case('MINERALP')
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'mineral p', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,UPSTREAM,MINERALP')
          this%name_mp = word

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
                       'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF,UPSTREAM')
         call StringToUpper(word)   

            select case(trim(word))
              case('CPOOL')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'name', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,UPSTREAM,CPOOL')
                  this%Upstream%name_c = word
                  this%Upstream%stoich = 1.0
              case('NPOOL')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'name', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,UPSTREAM,NPOOL')
                  this%Upstream%name_n = word
              case('PPOOL')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'name', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,UPSTREAM,PPOOL')
                 this%Upstream%name_p = word
              case('CNRATIO')
                  call InputReadDouble(input,option,tmp_real)  
                  call InputErrorMsg(input,option,'CN ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,UPSTREAMCPOOL')
                  this%Upstream%ratio_cn = tmp_real*CN_ratio_mass_to_mol
              case('CPRATIO')
                  call InputReadDouble(input,option,tmp_real)  
                  call InputErrorMsg(input,option,'CP ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,UPSTREAMCPOOL')
                  this%Upstream%ratio_cp = tmp_real*CP_ratio_mass_to_mol
              case('NCRATIO')
                  call InputReadDouble(input,option,tmp_real)  
                  call InputErrorMsg(input,option,'NC ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,UPSTREAMCPOOL')
                  this%Upstream%ratio_nc = tmp_real/CN_ratio_mass_to_mol
              case('PCRATIO')
                  call InputReadDouble(input,option,tmp_real)  
                  call InputErrorMsg(input,option,'CP ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,UPSTREAMCPOOL')
                  this%Upstream%ratio_pc = tmp_real/CP_ratio_mass_to_mol
              case('LITTER')
! added for CLM-CN-MICROBE  11/20/2013 
!-------------------------------------------------------------------------------
                   this%is_upstream_litter = PETSC_TRUE
!-------------------------------------------------------------------------------
! added for CLM-CN-MICROBE  11/20/2013 
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF,UPSTREAM keyword: ' // &
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
                       'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF,DOWNSTREAM')
          call StringToUpper(word)   

          select case(trim(word))
            case('CPOOL')
               call InputReadWord(input,option,word,PETSC_TRUE)
               call InputErrorMsg(input,option,'name', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,DOWNSTREAM,CPOOL')
               call InputReadDouble(input,option,tmp_real)  
               call InputErrorMsg(input,option,'stoichoimetric coefficient', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,DOWNSTREAMCPOOL')
               pool%name_c = word
               pool%stoich = tmp_real
            case('NPOOL')
               call InputReadWord(input,option,word,PETSC_TRUE)
               call InputErrorMsg(input,option,'name', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,DOWNSTREAM,NPOOL')
               pool%name_n = word
            case('PPOOL')
               call InputReadWord(input,option,word,PETSC_TRUE)
               call InputErrorMsg(input,option,'name', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,DOWNSTREAM,PPOOL')
               pool%name_p = word
            case('CNRATIO')
               call InputReadDouble(input,option,tmp_real)  
               call InputErrorMsg(input,option,'CN ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,DOWNSTREAMCPOOL')
               pool%ratio_cn = tmp_real*CN_ratio_mass_to_mol
            case('CPRATIO')
               call InputReadDouble(input,option,tmp_real)  
               call InputErrorMsg(input,option,'CP ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,DOWNSTREAMCPOOL')
               pool%ratio_cp = tmp_real*CP_ratio_mass_to_mol
            case('NCRATIO')
                call InputReadDouble(input,option,tmp_real)  
                call InputErrorMsg(input,option,'NC ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,DOWNSTREAMCPOOL')
                pool%ratio_nc = tmp_real/CN_ratio_mass_to_mol
            case('PCRATIO')
                call InputReadDouble(input,option,tmp_real)  
                call InputErrorMsg(input,option,'CP ratio', &
                        'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,DOWNSTREAMCPOOL')
                pool%ratio_pc = tmp_real/CP_ratio_mass_to_mol
            case('CNRATIO_FOLLOW_UPSTREAM')
               pool%bratio_cn_follow_upstream = PETSC_TRUE
            case('CPRATIO_FOLLOW_UPSTREAM')
               pool%bratio_cp_follow_upstream = PETSC_TRUE
            case default
               option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF,DOWNSTREAM keyword: ' // &
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
               'CHEMISTRY,REACTION_SANDBOX,CLM-CN_BF,RATE_CONSTANT')
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) then
            input%err_buf = 'CLM-CN_BF RATE CONSTANT UNITS'
            call InputDefaultMsg(input,option)
        else              
            this%rate_constant = rate_constant * &
            UnitsConvertToInternal(word,option)
        endif
      case('TURNOVER_TIME')
        call InputReadDouble(input,option,turnover_time)
        call InputErrorMsg(input,option,'turnover time', &
               'CHEMISTRY,REACTION_SANDBOX,CLM-CN_BF,REACTION')
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
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,FIRST ORDER')
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
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,MONOD')
        call InputReadDouble(input,option,tmp_real)  
        call InputErrorMsg(input,option,'half saturation constant', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,MONOD')
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
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,INHIBITION')
        call InputReadDouble(input,option,tmp_real)  
        call InputErrorMsg(input,option,'inhibition constant', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,INHIBITION')
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
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,NINHIBITION')
        this%NInhibitionCoef = tmp_real

      case('PINHIBITION')
        call InputReadDouble(input,option,tmp_real)  
        call InputErrorMsg(input,option,'P inhibition constant', &
                           'CHEMISTRY,REACTION_SANDBOX_CLM_CN_BF,NINHIBITION')
        this%PInhibitionCoef = tmp_real

      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo

  ! check to ensure that one of turnover time or rate constant is set.
  if (turnover_time > 0.d0 .and. rate_constant > 0.d0) then
       option%io_buffer = 'Only TURNOVER_TIME or RATE_CONSTANT may ' // &
        'be included in a CLM-CN_BF reaction definition, but not both. ' // &
        'See reaction with upstream pool "' // &
            trim(this%Upstream%name_c) // '".'
          call printErrMsg(option)
  endif  
end subroutine CLM_CN_BFRead

! ************************************************************************** !
!
! CLM_CN_BFSetup: 
! author: Guoping Tang
! date: 08/15/13
!
! ************************************************************************** !
subroutine CLM_CN_BFSetup(this,reaction,option)

  use Option_module
  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Immobile_Aux_module

  implicit none

  class(reaction_sandbox_CLM_CN_BF_type) :: this
  type(reaction_type)                  :: reaction
  type(option_type)                    :: option

  type(pool_type), pointer             :: cur_pool
  type(rate_type), pointer             :: cur_rate

  PetscInt :: i, icount, ispec
  PetscReal :: sum_stoich_c_prod
  PetscReal, parameter :: eps = 1.0d-150

  character(len=MAXWORDLENGTH) :: word

! mineral C  
  ispec = GetImmobileSpeciesIDFromName( &
          this%name_mc,reaction%immobile,PETSC_FALSE,option)  

  if(ispec > 0) then
       this%ispec_mc = ispec
       this%ires_mc = ispec + reaction%offset_immobile
       this%isaqueous_mc = PETSC_FALSE
  else
       ispec = GetPrimarySpeciesIDFromName(this%name_mc,reaction,PETSC_FALSE,option)
       if(ispec > 0) then
          this%ispec_mc = ispec
          this%ires_mc = ispec
          this%isaqueous_mc = PETSC_TRUE
       else
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
           'mineral c species ' // trim(this%name_mc)// &
           ' is neither an immobile species nor an aqueous species.'
          call printErrMsg(option)
       endif
  endif 

! upstream
  if(.not.associated(this%Upstream)) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
       ' upstream pool is not specified.'
     call printErrMsg(option)
  else
!     this%upstream_ispec_c = GetImmobileSpeciesIDFromName( &
!         this%Upstream%name_c,reaction%immobile,PETSC_FALSE,option)  
     ispec = GetImmobileSpeciesIDFromName( &
         this%Upstream%name_c,reaction%immobile,PETSC_FALSE,option)
     
     if(ispec > 0) then
          this%upstream_ispec_c = ispec
          this%upstream_ires_c = ispec + reaction%offset_immobile
          this%upstream_isaqueous_c = PETSC_FALSE
     else
          ispec = GetPrimarySpeciesIDFromName(this%Upstream%name_c,reaction,PETSC_FALSE,option)
          if(ispec > 0) then
             this%upstream_ispec_c = ispec
             this%upstream_ires_c = ispec
             this%upstream_isaqueous_c = PETSC_TRUE
          else
             option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
              'upstream c pool ' // trim(this%Upstream%name_c)// &
              ' is neither an immobile species nor an aqueous species.'
             call printErrMsg(option)
          endif
     endif

     if(trim(this%Upstream%name_n) /= '') then
!       this%upstream_ispec_n = GetImmobileSpeciesIDFromName( &
!         this%Upstream%name_n,reaction%immobile,PETSC_FALSE,option)  
       ispec = GetImmobileSpeciesIDFromName( &
         this%Upstream%name_n,reaction%immobile,PETSC_FALSE,option)  
       if(ispec > 0) then
          this%upstream_ispec_n = ispec
          this%upstream_ires_n = ispec + reaction%offset_immobile
          this%upstream_isaqueous_n = PETSC_FALSE
       else
          ispec = GetPrimarySpeciesIDFromName(this%Upstream%name_n,reaction,PETSC_FALSE,option)
          if(ispec > 0) then
             this%upstream_ispec_n = ispec
             this%upstream_ires_n = ispec
             this%upstream_isaqueous_n = PETSC_TRUE
          else
             option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
              'upstream N pool ' // trim(this%Upstream%name_n)// &
              ' is neither an immobile species nor an aqueous species.'
             call printErrMsg(option)
          endif
       endif
     else
         this%bfixed_upstream_cn = PETSC_TRUE
         if(this%Upstream%ratio_cn > eps) then
            this%upstream_stoich_n = 1.0d0/this%Upstream%ratio_cn
         elseif (this%Upstream%ratio_nc > eps) then
            this%upstream_stoich_n = this%Upstream%ratio_nc
         else
            this%upstream_stoich_n = 0.0d0 
         endif 
     endif
     if(trim(this%Upstream%name_p) /= '') then
!       this%upstream_ispec_p = GetImmobileSpeciesIDFromName( &
!         this%Upstream%name_p,reaction%immobile,PETSC_FALSE,option)  
       ispec = GetImmobileSpeciesIDFromName( &
         this%Upstream%name_p,reaction%immobile,PETSC_FALSE,option)  
       if(ispec > 0) then
          this%upstream_ispec_p = ispec
          this%upstream_ires_p = ispec + reaction%offset_immobile
          this%upstream_isaqueous_p = PETSC_FALSE
       else
          ispec = GetPrimarySpeciesIDFromName(this%Upstream%name_p,reaction,PETSC_FALSE,option)
          if(ispec > 0) then
             this%upstream_ispec_p = ispec
             this%upstream_ires_p = ispec
             this%upstream_isaqueous_p = PETSC_TRUE
          else
             option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
              'upstream P pool ' // trim(this%Upstream%name_p)// &
              ' is neither an immobile species nor an aqueous species.'
             call printErrMsg(option)
          endif
       endif
     else
         this%bfixed_upstream_cp = PETSC_TRUE
         if(this%Upstream%ratio_cp > eps) then
            this%upstream_stoich_p = 1.0d0/this%Upstream%ratio_cp
         elseif(this%Upstream%ratio_pc > eps) then
            this%upstream_stoich_p = this%Upstream%ratio_pc
         else
            this%upstream_stoich_p = 0.0d0 
         endif
     endif
  endif   

! downstream
  sum_stoich_c_prod = 0.0d0
  if(this%nDownstream >= 1) then
     allocate(this%downstream_ispec_c(this%nDownstream))
     allocate(this%downstream_ispec_n(this%nDownstream))
     allocate(this%downstream_ispec_p(this%nDownstream))

     allocate(this%downstream_ires_c(this%nDownstream))
     allocate(this%downstream_ires_n(this%nDownstream))
     allocate(this%downstream_ires_p(this%nDownstream))

     allocate(this%downstream_isaqueous_c(this%nDownstream))
     allocate(this%downstream_isaqueous_n(this%nDownstream))
     allocate(this%downstream_isaqueous_p(this%nDownstream))

     allocate(this%downstream_stoich_c(this%nDownstream))
     allocate(this%downstream_stoich_n(this%nDownstream))
     allocate(this%downstream_stoich_p(this%nDownstream))

!     allocate(this%downstream_cn(this%nDownstream))
!     allocate(this%downstream_cp(this%nDownstream))

     allocate(this%b_downstream_cn_fixed(this%nDownstream))
     allocate(this%b_downstream_cp_fixed(this%nDownstream))

     allocate(this%b_downstream_cn_follow_upstream(this%nDownstream))
     allocate(this%b_downstream_cp_follow_upstream(this%nDownstream))

     cur_pool => this%Downstream
     icount = 1
     do
       if(.not.associated(cur_pool)) exit
!        this%downstream_ispec_c(icount) = GetImmobileSpeciesIDFromName( &
!                    cur_pool%name_c,reaction%immobile,PETSC_FALSE,option)
       ispec = GetImmobileSpeciesIDFromName( &
          cur_pool%name_c,reaction%immobile,PETSC_FALSE,option) 
       if(ispec > 0) then
          this%downstream_ispec_c(icount) = ispec
          this%downstream_ires_c(icount) = ispec + reaction%offset_immobile
          this%downstream_isaqueous_c(icount) = PETSC_FALSE
       else
          ispec = GetPrimarySpeciesIDFromName(cur_pool%name_c,reaction,PETSC_FALSE,option)
          if(ispec > 0) then
             this%downstream_ispec_c(icount) = ispec
             this%downstream_ires_c(icount) = ispec
             this%downstream_isaqueous_c(icount) = PETSC_TRUE
          else
             option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
              'downstream C pool ' // trim(cur_pool%name_c)// &
              ' is neither an immobile species nor an aqueous species.'
             call printErrMsg(option)
          endif
       endif
       this%downstream_stoich_c(icount) = cur_pool%stoich
       sum_stoich_c_prod = sum_stoich_c_prod + cur_pool%stoich

       if(cur_pool%bratio_cn_follow_upstream) then
           this%b_downstream_cn_follow_upstream(icount) = PETSC_TRUE
           if(this%bfixed_upstream_cn) then
              this%b_downstream_cn_fixed(icount) = PETSC_TRUE
!              this%downstream_cn(icount) = this%Upstream%ratio_cn 
!              this%downstream_stoich_n(icount) = 1.0d0/this%downstream_cn(icount)
              this%downstream_stoich_n(icount) = this%upstream_stoich_n  !1.0d0/this%Upstream%ratio_cn
           endif
        else
           this%b_downstream_cn_follow_upstream(icount) = PETSC_FALSE
        endif

        if(trim(cur_pool%name_n) /= '') then
!           this%downstream_ispec_n(icount) = GetImmobileSpeciesIDFromName( &
!                    cur_pool%name_n,reaction%immobile,PETSC_FALSE,option)  
          ispec = GetImmobileSpeciesIDFromName( &
          cur_pool%name_n,reaction%immobile,PETSC_FALSE,option) 
          if(ispec > 0) then
             this%downstream_ispec_n(icount) = ispec
             this%downstream_ires_n(icount) = ispec + reaction%offset_immobile
             this%downstream_isaqueous_n(icount) = PETSC_FALSE
          else
             ispec = GetPrimarySpeciesIDFromName(cur_pool%name_n,reaction,PETSC_FALSE,option)
             if(ispec > 0) then
                this%downstream_ispec_n(icount) = ispec
                this%downstream_ires_n(icount) = ispec
                this%downstream_isaqueous_n(icount) = PETSC_TRUE
             else
                option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
                 'downstream N pool ' // trim(cur_pool%name_n)// &
                 ' is neither an immobile species nor an aqueous species.'
                call printErrMsg(option)
             endif
          endif
        else
           if(.not. (this%b_downstream_cn_follow_upstream(icount))) then
!             this%downstream_cn(icount) = cur_pool%ratio_cn
!             if(this%downstream_cn(icount) .GT. 1.0d-10) then
             this%b_downstream_cn_fixed(icount) = PETSC_TRUE
             if(cur_pool%ratio_cn > eps) then
              this%downstream_stoich_n(icount) = 1.0d0/cur_pool%ratio_cn
             elseif(cur_pool%ratio_nc > eps) then
              this%downstream_stoich_n(icount) = cur_pool%ratio_nc 
             else
              this%downstream_stoich_n(icount) = 0.0d0
             endif 
           endif
        endif

        if(cur_pool%bratio_cp_follow_upstream) then
           this%b_downstream_cp_follow_upstream(icount) = PETSC_TRUE
           if(this%bfixed_upstream_cp) then
              this%b_downstream_cp_fixed(icount) = PETSC_TRUE
!              this%downstream_cp(icount) = this%Upstream%ratio_cp 
!              this%downstream_stoich_p(icount) = 1.0d0/this%downstream_cp(icount)
              this%downstream_stoich_p(icount) = this%upstream_stoich_p  !1.0d0/this%Upstream%ratio_cp
           endif
        else
           this%b_downstream_cp_follow_upstream(icount) = PETSC_FALSE
        endif

        if(trim(cur_pool%name_p) /= '') then
!           this%downstream_ispec_p(icount) = GetImmobileSpeciesIDFromName( &
!                    cur_pool%name_p,reaction%immobile,PETSC_FALSE,option)  
          ispec = GetImmobileSpeciesIDFromName( &
          cur_pool%name_p,reaction%immobile,PETSC_FALSE,option) 
          if(ispec > 0) then
             this%downstream_ispec_p(icount) = ispec
             this%downstream_ires_p(icount) = ispec + reaction%offset_immobile
             this%downstream_isaqueous_p(icount) = PETSC_FALSE
          else
             ispec = GetPrimarySpeciesIDFromName(cur_pool%name_p,reaction,PETSC_FALSE,option)
             if(ispec > 0) then
                this%downstream_ispec_p(icount) = ispec
                this%downstream_ires_p(icount) = ispec
                this%downstream_isaqueous_p(icount) = PETSC_TRUE
             else
                option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
                 'downstream P pool ' // trim(cur_pool%name_p)// &
                 ' is neither an immobile species nor an aqueous species.'
                call printErrMsg(option)
             endif
          endif
        else
           if(.not. (this%b_downstream_cp_follow_upstream(icount))) then
!             this%downstream_cp(icount) = cur_pool%ratio_cp
!             if(this%downstream_cp(icount) .GT. 1.0d-10) then
             this%b_downstream_cp_fixed(icount) = PETSC_TRUE
             if(cur_pool%ratio_cp > eps) then
              this%b_downstream_cp_fixed(icount) = PETSC_TRUE
!              this%downstream_stoich_p(icount) = 1.0d0/this%downstream_cp(icount)
              this%downstream_stoich_p(icount) = 1.0d0/cur_pool%ratio_cp
             elseif(cur_pool%ratio_pc > eps) then
              this%downstream_stoich_p(icount) = cur_pool%ratio_pc
             else
              this%downstream_stoich_p(icount) = 0.0d0
             endif
           endif
        endif

        cur_pool => cur_pool%next
        icount = icount + 1
     enddo 
     if(sum_stoich_c_prod > 1.0d0) then
         option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check:' // &
            ' downstream C pools stoichoimetric coefficient sum > 1.'
         call printErrMsg(option)
     else
         this%mineral_c_stoich = 1.0d0 - sum_stoich_c_prod
     endif
  else
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check:' // &
            ' downstream pool not specified.'
     call printErrMsg(option)
  endif   

! check if N is included
  if((this%Upstream%ratio_cn < 0) .AND. (this%upstream_ispec_n < 0)) then
      write(option%fid_out,*) 'Neither CN ratio nor valid upstream N pool is specified. N is not considered.'
  else
     do i = 1, this%nDownstream
!        if((this%downstream_cn(i) .LT. 0.0d0) .and. (this%downstream_ispec_n(i) .LT. 0) &
        if((this%downstream_stoich_n(i) < 0.0d0) .AND. (this%downstream_ispec_n(i) < 0) &
          .AND. (.NOT.this%b_downstream_cn_follow_upstream(i))) then
          write(option%fid_out,*) 'Neither CN ratio nor valid N pool is specified for downstream. N is not considered.'
          exit
        endif
     enddo
     this%bNEnabled = PETSC_TRUE  
     
!     word = 'N'
!     this%ispec_mn = GetImmobileSpeciesIDFromName( &
!            word,reaction%immobile,PETSC_FALSE,option)  
!     word = 'NH3(aq)'
!     this%ispec_mn = GetPrimarySpeciesIDFromName(word,reaction,option)
! mineral N  
     ispec = GetImmobileSpeciesIDFromName( &
             this%name_mn,reaction%immobile,PETSC_FALSE,option)  

     if(ispec > 0) then
       this%ispec_mn = ispec
       this%ires_mn = ispec + reaction%offset_immobile
       this%isaqueous_mn = PETSC_FALSE
     else
       ispec = GetPrimarySpeciesIDFromName(this%name_mn,reaction,PETSC_FALSE,option)
       if(ispec > 0) then
          this%ispec_mn = ispec
          this%ires_mn = ispec
          this%isaqueous_mn = PETSC_TRUE
       else
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
           'mineral n species ' // trim(this%name_mc)// &
          ' is neither an immobile species nor an aqueous species.'
          call printErrMsg(option)
       endif
     endif 
  endif

! check if P is included
  if((this%Upstream%ratio_cp < 0) .and. (this%upstream_ispec_p < 0)) then
      write(option%fid_out,*) 'Neither CP ratio nor valid upstream P pool is specified. P is not considered.'
  else
     do i = 1, this%nDownstream
!        if((this%downstream_cp(i) .LT. 0.0d0) .and. (this%downstream_ispec_p(i) .LT. 0) &
        if((this%downstream_stoich_p(i) < 0.0d0) .and. (this%downstream_ispec_p(i) < 0) &
          .and. (.not.this%b_downstream_cn_follow_upstream(i))) then
          write(option%fid_out,*) 'Neither CP ratio nor valid P pool is specified for downstream. P is not considered.'
          exit
        endif
     enddo
     this%bPEnabled = PETSC_TRUE  
!     word = 'P'
!     this%ispec_mp = GetImmobileSpeciesIDFromName( &
!            word,reaction%immobile,PETSC_FALSE,option)  
     ispec = GetImmobileSpeciesIDFromName( &
             this%name_mp,reaction%immobile,PETSC_FALSE,option)  

     if(ispec > 0) then
       this%ispec_mp = ispec
       this%ires_mp = ispec + reaction%offset_immobile
       this%isaqueous_mp = PETSC_FALSE
     else
       ispec = GetPrimarySpeciesIDFromName(this%name_mp,reaction,PETSC_FALSE,option)
       if(ispec > 0) then
          this%ispec_mp = ispec
          this%ires_mp = ispec
          this%isaqueous_mp = PETSC_TRUE
       else
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
           'mineral p species ' // trim(this%name_mc)// &
           ' is neither an immobile species nor an aqueous species.'
          call printErrMsg(option)
       endif
     endif 
  endif

! first order rate terms
  if(this%nFirstOrder >= 1) then
     allocate(this%ispec_1st(this%nFirstOrder))
     allocate(this%isaqueous_1st(this%nFirstOrder))

     cur_rate => this%FirstOrder
     icount = 1
     do
        if(.not.associated(cur_rate)) exit
!        this%ispec_1st(icount) = GetImmobileSpeciesIDFromName( &
!                    cur_rate%name,reaction%immobile,PETSC_FALSE,option)
        ispec = GetImmobileSpeciesIDFromName( &
            cur_rate%name,reaction%immobile,PETSC_FALSE,option) 
        if(ispec > 0) then
             this%ispec_1st(icount) = ispec
             this%isaqueous_1st(icount) = PETSC_FALSE
        else
             ispec = GetPrimarySpeciesIDFromName(cur_rate%name,reaction,PETSC_FALSE,option)
             if(ispec > 0) then
                this%ispec_1st(icount) = ispec
                this%isaqueous_1st(icount) = PETSC_TRUE
             else
                option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
                 'First order term C pool ' // trim(cur_rate%name)// &
                 ' is neither an immobile species nor an aqueous species.'
                call printErrMsg(option)
             endif
        endif
        cur_rate => cur_rate%next
        icount = icount + 1
     enddo 
  endif   

! monod rate terms
  if(this%nMonod >= 1) then
     allocate(this%ispec_mnd(this%nMonod))
     allocate(this%half_saturation(this%nMonod))
     allocate(this%isaqueous_mnd(this%nMonod))

     cur_rate => this%Monod
     icount = 1
     do
        if(.not.associated(cur_rate)) exit
!        this%ispec_mnd(icount) = GetImmobileSpeciesIDFromName( &
!                    cur_rate%name,reaction%immobile,PETSC_FALSE,option)
        ispec = GetImmobileSpeciesIDFromName( &
                cur_rate%name,reaction%immobile,PETSC_FALSE,option) 
        if(ispec > 0) then
             this%ispec_mnd(icount) = ispec
             this%isaqueous_mnd(icount) = PETSC_FALSE
        else
             ispec = GetPrimarySpeciesIDFromName(cur_rate%name,reaction,PETSC_FALSE,option)
             if(ispec > 0) then
                this%ispec_mnd(icount) = ispec
                this%isaqueous_mnd(icount) = PETSC_TRUE
             else
                option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
                 'Monod term C pool ' // trim(cur_rate%name)// &
                 ' is neither an immobile species nor an aqueous species.'
                call printErrMsg(option)
             endif
        endif
        this%half_saturation(icount) = cur_rate%value
        cur_rate => cur_rate%next
        icount = icount + 1
     enddo 
  endif   

! inhibition rate terms
  if(this%nInhibition >= 1) then
     allocate(this%ispec_inh(this%nInhibition))
     allocate(this%inhibition_coef(this%nInhibition))
     allocate(this%isaqueous_inh(this%nInhibition))

     cur_rate => this%Inhibition
     icount = 1
     do
        if(.not.associated(cur_rate)) exit
!        this%ispec_inh(icount) = GetImmobileSpeciesIDFromName( &
!                    cur_rate%name,reaction%immobile,PETSC_FALSE,option)
        ispec = GetImmobileSpeciesIDFromName( &
                cur_rate%name,reaction%immobile,PETSC_FALSE,option) 
        if(ispec > 0) then
             this%ispec_inh(icount) = ispec
             this%isaqueous_inh(icount) = PETSC_FALSE
        else
             ispec = GetPrimarySpeciesIDFromName(cur_rate%name,reaction,PETSC_FALSE,option)
             if(ispec > 0) then
                this%ispec_inh(icount) = ispec
                this%isaqueous_inh(icount) = PETSC_TRUE
             else
                option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
                 'Inhibition term C pool ' // trim(cur_rate%name)// &
                 ' is neither an immobile species nor an aqueous species.'
                call printErrMsg(option)
             endif
        endif
        this%inhibition_coef(icount) = cur_rate%value
        cur_rate => cur_rate%next
        icount = icount + 1
     enddo 
  endif   

end subroutine CLM_CN_BFSetup

! ************************************************************************** !
!
! CLM_CN_BFReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 08/16/13
!
! ************************************************************************** !
subroutine CLM_CN_BFReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,porosity,volume,reaction, &
                         option,local_id)

  use Option_module
  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  
#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data !, only : rate_plantnuptake_pf 
#endif
  
  implicit none
  
#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif
  
  class(reaction_sandbox_CLM_CN_BF_type) :: this  
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
  PetscInt :: i, j, ires, ires_j, ires_mc, ires_mn, ires_mp
  PetscInt :: ires_uc, ires_un, ires_up, ires_dc, ires_dn, ires_dp
  PetscInt :: ispec_bacteria, ispec_fungi
  PetscReal :: conc, c_uc, c_un, c_up, c_mc, c_mn, c_mp, c_dc, c_dn, c_dp
  PetscReal :: L_water
  PetscReal :: rate, drate, drate_uc, drate_n, drate_p 
  PetscReal :: tmp_real, respiration_fraction

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
  PetscBool :: is_pool_in_rxn
  PetscBool :: isaqueous
  PetscErrorCode :: ierr
  PetscReal :: Lwater_m3
  character(len=MAXWORDLENGTH) :: word

  ! temperature response function 
  select case(this%temperature_response_function)
      case(TEMPERATURE_RESPONSE_FUNCTION_Q10)
! CLM4.5 temperature response function
          tc = global_auxvar%temp(1)
          F_t = this%Q10 ** ((tc - 25.0d0) / 10.0d0)
      case(TEMPERATURE_RESPONSE_FUNCTION_CLM4) 
! CLM-CN temperature response function
  ! inhibition due to temperature
  ! Equation: F_t = exp(308.56*(1/17.02 - 1/(T - 227.13)))

          temp_K = global_auxvar%temp(1) + 273.15d0

          if(temp_K .GT. 227.15d0) then
            F_t = exp(308.56d0*(one_over_71_02 - 1.d0/(temp_K - 227.13d0)))
          else
            F_t = 0.0d0
            return
          endif
  end select     
  ! inhibition due to moisture content
  ! Equation: F_theta = log(theta_min/theta) / log(theta_min/theta_max)
  ! Assumptions: theta is saturation
  !              theta_min = 0.01, theta_max = 1.

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
     return
  endif

  call VecRestoreArrayReadF90(clm_pf_idata%sucsat_pf, sucsat_pf_loc, ierr)
  call VecRestoreArrayReadF90(clm_pf_idata%soilpsi_pf, soilpsi_pf_loc, ierr)
#else
  F_theta = 1.0d0 
#endif

!  F_theta = log(theta_min/max(theta_min,global_auxvar%sat(1))) &  ! theta could be zero, which will cause math issue
!          * one_over_log_theta_min

  offset = reaction%offset_immobile

  Lwater_m3 = porosity * global_auxvar%sat(iphase)*1.d3

  rate = this%rate_constant * volume * F_t * F_theta

  drate_uc = rate    ! to avoid rate/upstream c when upstream c goes to 0

  if(this%upstream_isaqueous_c) then
     c_uc = rt_auxvar%total(this%upstream_ispec_c,iphase)
     c_uc = c_uc * Lwater_m3 
  else 
     c_uc = rt_auxvar%immobile(this%upstream_ispec_c) 
  endif
       
! first order term
  do i = 1, this%nFirstOrder
     if(this%isaqueous_1st(i)) then
        conc = rt_auxvar%total(this%ispec_1st(i), iphase) 
        conc = conc * Lwater_m3 
     else
        conc = rt_auxvar%immobile(this%ispec_1st(i))  
     endif
     rate = rate * conc 
     if(this%ispec_1st(i) /= this%upstream_ispec_c) then
       drate_uc = drate_uc * conc 
     endif
  enddo 

! monod term
  do i = 1, this%nMonod
     if(this%isaqueous_mnd(i)) then
        conc = rt_auxvar%total(this%ispec_mnd(i), iphase) 
        conc = conc * Lwater_m3 
     else
        conc = rt_auxvar%immobile(this%ispec_mnd(i))  
     endif
     if(this%ispec_mnd(i) /= this%upstream_ispec_c) then
       drate_uc = drate_uc * conc/(conc + this%half_saturation(i))  
     endif
     rate = rate * conc/(conc + this%half_saturation(i))  
  enddo 

! inhibition term
  do i = 1, this%nInhibition
     if(this%isaqueous_inh(i)) then
        conc = rt_auxvar%total(this%ispec_inh(i),iphase) 
        conc = conc * Lwater_m3 
     else
        conc = rt_auxvar%immobile(this%ispec_inh(i))
     endif
     rate = rate * this%inhibition_coef(i)/(conc + this%inhibition_coef(i))  
  enddo 

! N limiting
  if((this%bNEnabled)) then
! upstream CN
    if(.not.this%bfixed_upstream_cn) then
!       if(rt_auxvar%immobile(this%upstream_ispec_c) < eps) then
       if(this%upstream_isaqueous_n) then
          c_un = rt_auxvar%total(this%upstream_ispec_n, iphase)
          c_un = c_un * Lwater_m3 
       else 
          c_un = rt_auxvar%immobile(this%upstream_ispec_n) 
       endif
       if(c_uc < eps) then
          write(option%fid_out,*) 'Upstream C concentration is 0 in CN ratio calculation.'
       endif 
!       this%upstream_stoich_n = rt_auxvar%immobile(this%upstream_ispec_n) & 
!                  /rt_auxvar%immobile(this%upstream_ispec_c)
       this%upstream_stoich_n = c_un/c_uc 
    endif

! mineral N
    if(this%isaqueous_mn) then
       c_mn = rt_auxvar%total(this%ispec_mn, iphase)
       c_mn = c_mn * Lwater_m3 
    else
       c_mn = rt_auxvar%immobile(this%ispec_mn) 
    endif

! added for CLM-CN-MICROBE  11/20/2013 
!-------------------------------------------------------------------------------
!  if (CLM_CN_MICROBE) then
  if (this%is_upstream_litter) then

     if( this%nDownstream /= 2 ) then
             option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
              'There should be 2 downstream pools (bacteria and fungi) for CLM-CN-Microbe.'
             call printErrMsg(option)
     endif

     word = 'Bacteria'
     ispec_bacteria = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option) 
     if(ispec_bacteria <= 0) then
          ispec_bacteria = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)
          if(ispec_bacteria <= 0) then
             option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
              'Bacteria pool ' // trim(word)// &
              ' not specified for CLM-CN-Microbe.'
             call printErrMsg(option)
          endif
     endif   

     word = 'Fungi'
     ispec_fungi = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option) 
     if(ispec_fungi <= 0) then
          ispec_fungi = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)
          if(ispec_fungi <= 0) then
             option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,CLM_CN_BF check: ' // &
              'Fungi pool ' // trim(word)// &
              ' not specified for CLM-CN-Microbe.'
             call printErrMsg(option)
          endif
     endif   

!     call CLM_CN_BF_Update_Stoich(this, ispec_bacteria, ispec_fungi) 
! upstream pool is litter
     if(this%is_upstream_litter) then
!    Sinsabaugh et al. 2013 Ecology Letters, 16, 930-939
       respiration_fraction = CN_ratio_microbe * this%upstream_stoich_n

       if(respiration_fraction > CUE_max) then
          respiration_fraction = CUE_max
       endif

!    c pools     
       this%mineral_c_stoich = respiration_fraction

       do i = 1, this%nDownstream
        if(this%downstream_ispec_c(i) == ispec_bacteria) then
           this%downstream_stoich_c(i) = fraction_bacteria * (1.0d0  - respiration_fraction)
        else
           this%downstream_stoich_c(i) = (1.0d0 - fraction_bacteria) * (1.0d0  - respiration_fraction) 
        endif
       enddo
     endif
  endif
!-------------------------------------------------------------------------------

! downstream CN
    this%mineral_n_stoich = this%upstream_stoich_n   ! start from upstream N, substract downstream N
    do i = 1, this%nDownstream
       if(.not.this%b_downstream_cn_fixed(i)) then
         if(this%b_downstream_cn_follow_upstream(i)) then
          this%downstream_stoich_n(i) = this%upstream_stoich_n  
         else
          if(this%downstream_isaqueous_c(i)) then
             c_dc = rt_auxvar%total(this%downstream_ispec_c(i), iphase)
             c_dc = c_dc * Lwater_m3 
          else 
             c_dc = rt_auxvar%immobile(this%downstream_ispec_c(i)) 
          endif
          if(this%downstream_isaqueous_n(i)) then
             c_dn = rt_auxvar%total(this%downstream_ispec_n(i), iphase)
             c_dn = c_dn * Lwater_m3 
          else 
             c_dn = rt_auxvar%immobile(this%downstream_ispec_n(i)) 
          endif

!          if(rt_auxvar%immobile(this%downstream_ispec_c(i)) < eps) then
          if(c_dc < eps) then
             write(option%fid_out,*) 'Downtream C concentration 0 in CN ratio calculation.'
          endif 
!          this%downstream_stoich_n(i)=rt_auxvar%immobile(this%downstream_ispec_n(i)) & 
!                 / rt_auxvar%immobile(this%downstream_ispec_c(i))
          this%downstream_stoich_n(i) = c_dn / c_dc
         endif
       endif
       this%mineral_n_stoich = this%mineral_n_stoich - this%downstream_stoich_n(i) &
                 * this%downstream_stoich_c(i)
    enddo         

    if(this%mineral_n_stoich < 0.0d0 .and. this%NInhibitionCoef > 0.0d0) then
!       conc = rt_auxvar%immobile(this%ispec_mn)  
!       drate_n = rate / (conc + this%NInhibitionCoef) 
!       drate_uc = drate_uc * conc/(conc + this%NInhibitionCoef)  
!       rate = rate * conc/(conc + this%NInhibitionCoef)  ! N limiting  
       drate_n = rate / (c_mn + this%NInhibitionCoef) 
       drate_uc = drate_uc * c_mn/(c_mn + this%NInhibitionCoef)  
       rate = rate * c_mn/(c_mn + this%NInhibitionCoef)  ! N limiting  
    endif
  endif

! P limiting
  if((this%bPEnabled)) then
! upstream CP
    if(.not.this%bfixed_upstream_cp) then
!       if(rt_auxvar%immobile(this%upstream_ispec_c) < eps) then
       if(this%upstream_isaqueous_p) then
          c_up = rt_auxvar%total(this%upstream_ispec_p, iphase)
          c_up = c_up * Lwater_m3 
       else 
          c_up = rt_auxvar%immobile(this%upstream_ispec_p) 
       endif
       if(c_uc < eps) then
          write(option%fid_out,*) 'Upstream C concentration is 0 in CP ratio calculation.'
       endif 
!       this%upstream_stoich_p = rt_auxvar%immobile(this%upstream_ispec_p) & 
!                  / rt_auxvar%immobile(this%upstream_ispec_c)
       this%upstream_stoich_p = c_up / c_uc
    endif

! mineral P
    if(this%isaqueous_mp) then
       c_mp = rt_auxvar%total(this%ispec_mp, iphase)
       c_mp = c_mp * Lwater_m3 
    else
       c_mp = rt_auxvar%immobile(this%ispec_mp) 
    endif

! downstream CP
    this%mineral_p_stoich = this%upstream_stoich_p   ! start from upstream P, substract downstream P
    do i = 1, this%nDownstream
       if(.not.this%b_downstream_cp_fixed(i)) then
         if(this%b_downstream_cp_follow_upstream(i)) then
          this%downstream_stoich_p(i) = this%upstream_stoich_p  
         else
          if(this%downstream_isaqueous_c(i)) then
             c_dc = rt_auxvar%total(this%downstream_ispec_c(i), iphase)
             c_dc = c_dc * Lwater_m3 
          else 
             c_dc = rt_auxvar%immobile(this%downstream_ispec_c(i)) 
          endif
          if(this%downstream_isaqueous_p(i)) then
             c_dp = rt_auxvar%total(this%downstream_ispec_p(i), iphase)
             c_dp = c_dp * Lwater_m3 
          else 
             c_dp = rt_auxvar%immobile(this%downstream_ispec_p(i)) 
          endif

!          if(rt_auxvar%immobile(this%downstream_ispec_c(i)) < eps) then
          if(c_dc < eps) then
             write(option%fid_out,*) 'Downtream C concentration 0 in CP ratio calculation.'
          endif 
!          this%downstream_stoich_p(i)=rt_auxvar%immobile(this%downstream_ispec_p(i)) & 
!                 / rt_auxvar%immobile(this%downstream_ispec_c(i))
          this%downstream_stoich_p(i) = c_dp / c_dc
         endif
       endif
       this%mineral_p_stoich = this%mineral_p_stoich - this%downstream_stoich_p(i) &
                 * this%downstream_stoich_c(i)
    enddo         

    if(this%mineral_p_stoich < 0.0d0 .and. this%PInhibitionCoef > 0.0d0) then
!       conc = rt_auxvar%immobile(this%ispec_mp)  
!       drate_p = rate / (conc + this%PInhibitionCoef) 
!       rate = rate * conc/(conc + this%PInhibitionCoef)  ! P limiting  
!       drate_n = drate_n * conc / (conc + this%PInhibitionCoef) 
!       drate_uc = drate_uc * conc / (conc + this%PInhibitionCoef)  
       drate_p = rate / (c_mp + this%PInhibitionCoef) 
       rate = rate * c_mp/(c_mp + this%PInhibitionCoef)  ! P limiting  
       drate_n = drate_n * c_mp / (c_mp + this%PInhibitionCoef) 
       drate_uc = drate_uc * c_mp / (c_mp + this%PInhibitionCoef)  
    endif
  endif

! residuals

  ires_uc = this%upstream_ires_c
  Residual(ires_uc) = Residual(ires_uc) + rate

! mineral C
  ires_mc = this%ires_mc 
  Residual(ires_mc) = Residual(ires_mc) - rate * this%mineral_c_stoich

! downstream C
  do i = 1, this%nDownstream
     ires_dc = this%downstream_ires_c(i)
     Residual(ires_dc) = Residual(ires_dc) - rate * this%downstream_stoich_c(i)
  enddo 

  if(this%bNEnabled) then
! upstream N
     if(.not.this%bfixed_upstream_cn) then
        ires_un = this%upstream_ires_n
        Residual(ires_un) = Residual(ires_un) + rate * this%upstream_stoich_n
     endif

! mineral N
     ires_mn = this%ires_mn

     Residual(ires_mn) = Residual(ires_mn) - rate * this%mineral_n_stoich

! downstream N
     do i = 1, this%nDownstream
       if(.not.this%b_downstream_cn_fixed(i)) then
          ires_dn = this%downstream_ires_n(i)
          Residual(ires_dn) = Residual(ires_dn) - rate * &
                               this%downstream_stoich_n(i) * &
                               this%downstream_stoich_c(i)
       endif
     enddo 
  endif

  if(this%bPEnabled) then

! upstream P
     if(.not.this%bfixed_upstream_cp) then
        ires_up = this%upstream_ires_p
        Residual(ires_up) = Residual(ires_up) + rate * this%upstream_stoich_p
     endif

     ires_mp = this%ires_mp
      
     Residual(ires_mp) = Residual(ires_mp) - rate * this%mineral_p_stoich

! downstream P
     do i = 1, this%nDownstream
       if(.not.this%b_downstream_cp_fixed(i)) then
          ires_dp = this%downstream_ires_p(i)
          Residual(ires_dp) = Residual(ires_dp) - rate * &
                               this%downstream_stoich_p(i) * & 
                               this%downstream_stoich_c(i)      
       endif
     enddo 
  endif

! Jacobian
  if (.not.compute_derivative) return 

  if(this%bNEnabled) then
     if(.not.this%bfixed_upstream_cn) then
        tmp_real = c_un / c_uc
        if(this%upstream_isaqueous_c) then
           Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) - drate_uc * tmp_real * & 
                       rt_auxvar%aqueous%dtotal(ires_un,ires_uc,iphase)
           Jacobian(ires_mn,ires_uc) = Jacobian(ires_mn,ires_uc) + drate_uc * tmp_real * &
                       rt_auxvar%aqueous%dtotal(ires_mn,ires_uc,iphase)
        else
           Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) - drate_uc * tmp_real !rate * tmp_real  dnu/dunRuc
           Jacobian(ires_mn,ires_uc) = Jacobian(ires_mn, ires_uc) + drate_uc * tmp_real !rate * tmp_real  dnu/dmnRuc 
        endif

        if(this%upstream_isaqueous_n) then
           Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) + drate_uc * & 
                       rt_auxvar%aqueous%dtotal(ires_un,ires_un,iphase)
           Jacobian(ires_mn,ires_un) = Jacobian(ires_mn,ires_un) - drate_uc  * &
                       rt_auxvar%aqueous%dtotal(ires_mn,ires_un,iphase)
        else
           Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) + drate_uc !rate / &                    dnu/dunRun   
                                    !rt_auxvar%immobile(this%upstream_ispec_c)                        
           Jacobian(ires_mn,ires_un) = Jacobian(ires_mn,ires_un) - drate_uc   !rate / &                    dnu/dmnRun  
                                    !rt_auxvar%immobile(this%upstream_ispec_c)
        endif

        do i = 1, this%nDownstream
           if(this%b_downstream_cn_follow_upstream(i)) then
              ires_dn = this%downstream_ires_n(i)
            if(this%upstream_isaqueous_c) then
              Jacobian(ires_dn,ires_uc) = Jacobian(ires_dn,ires_uc) + &
                       drate_uc * tmp_real * this%downstream_stoich_c(i) * & 
                       rt_auxvar%aqueous%dtotal(ires_dn,ires_uc,iphase)
              Jacobian(ires_mn, ires_uc) = Jacobian(ires_mn,  ires_uc) - &
                       drate_uc * tmp_real * this%downstream_stoich_c(i) * &
                       rt_auxvar%aqueous%dtotal(ires_mn,ires_uc,iphase)
            else
              Jacobian(ires_dn,ires_uc) = Jacobian(ires_dn,ires_uc) + & 
                       drate_uc * tmp_real * this%downstream_stoich_c(i) 
              Jacobian(ires_mn,ires_uc) = Jacobian(ires_mn,ires_uc) - &
                       drate_uc * tmp_real * this%downstream_stoich_c(i)
            endif

            if(this%upstream_isaqueous_n) then
              Jacobian(ires_dn,ires_un) = Jacobian(ires_dn,ires_un) - &
                       drate_uc * this%downstream_stoich_c(i) * & 
                       rt_auxvar%aqueous%dtotal(ires_dn,ires_un,iphase)
              Jacobian(ires_mn,ires_un) = Jacobian(ires_mn,ires_un) + &
                       drate_uc * this%downstream_stoich_c(i) * &
                       rt_auxvar%aqueous%dtotal(ires_mn,ires_un,iphase)
            else
              Jacobian(ires_dn,ires_un) = Jacobian(ires_dn,ires_un) - &
                       drate_uc * this%downstream_stoich_c(i) 
              Jacobian(ires_mn,ires_un) = Jacobian(ires_mn,ires_un) + &
                       drate_uc * this%downstream_stoich_c(i)  
            endif
           endif
        enddo
     endif

     do i = 1, this%nDownstream
        if((.not.this%b_downstream_cn_fixed(i)).and. &
           (.not.this%b_downstream_cn_follow_upstream(i))) then

          if(this%downstream_isaqueous_c(i)) then
             ires_dc = this%downstream_ispec_c(i)
             c_dc = rt_auxvar%total(this%downstream_ispec_c(i), iphase)
             c_dc = c_dc * Lwater_m3 
          else 
             ires_dc = this%downstream_ires_c(i) !reaction%offset_immobile + this%downstream_ispec_c(i)
             c_dc = rt_auxvar%immobile(this%downstream_ispec_c(i)) 
          endif

          if(this%downstream_isaqueous_n(i)) then
             ires_dn = this%downstream_ispec_n(i)
             c_dn = rt_auxvar%total(this%downstream_ispec_n(i), iphase)
             c_dn = c_dn * Lwater_m3 
          else 
             ires_dn = this%downstream_ires_n(i)  !reaction%offset_immobile + this%downstream_ispec_n(i)
             c_dn = rt_auxvar%immobile(this%downstream_ispec_n(i)) 
          endif

           tmp_real = c_dn / c_dc / c_dc
           
          if(this%downstream_isaqueous_c(i)) then
           Jacobian(ires_dn,ires_dc) = Jacobian(ires_dn,ires_dc) + &
                    rate * tmp_real * this%downstream_stoich_c(i) * &
                    rt_auxvar%aqueous%dtotal(ires_dn,ires_dc,iphase)
           Jacobian(ires_mn,ires_dc) = Jacobian(ires_mn, ires_dc) - &
                    rate * tmp_real * this%downstream_stoich_c(i) * &
                    rt_auxvar%aqueous%dtotal(ires_mn,ires_dc,iphase)
          else
           Jacobian(ires_dn,ires_dc) = Jacobian(ires_dn,ires_dc) + &
                    rate * tmp_real * this%downstream_stoich_c(i)
           Jacobian(ires_mn, ires_dc) = Jacobian(ires_mn, ires_dc) - &
                    rate * tmp_real * this%downstream_stoich_c(i)
          endif

          if(this%downstream_isaqueous_n(i)) then
           Jacobian(ires_dn,ires_dn) = Jacobian(ires_dn,ires_dn) - &
                    rate / c_dc * this%downstream_stoich_c(i) * &
                    rt_auxvar%aqueous%dtotal(ires_dn,ires_dn,iphase)
           Jacobian(ires_mn,ires_dn) = Jacobian(ires_mn,ires_dn) + &
                    rate / c_dc * this%downstream_stoich_c(i) * &
                    rt_auxvar%aqueous%dtotal(ires_mn,ires_dn,iphase)
          else
           Jacobian(ires_dn,ires_dn) = Jacobian(ires_dn,ires_dn) - &
                    rate / c_dc * this%downstream_stoich_c(i)
           Jacobian(ires_mn,ires_dn) = Jacobian(ires_mn,ires_dn) + &
                    rate / c_dc * this%downstream_stoich_c(i)
          endif
        endif
     enddo 
  endif

  if(this%bPEnabled) then
     if(.not.this%bfixed_upstream_cp) then
        tmp_real = c_up /c_uc
       if(this%upstream_isaqueous_c) then
        Jacobian(ires_up,ires_uc) = Jacobian(ires_up,ires_uc) - &
                 drate_uc * tmp_real * rt_auxvar%aqueous%dtotal(ires_up,ires_uc,iphase)
        Jacobian(ires_mp, ires_uc) = Jacobian(ires_mp, ires_uc) + &
                 drate_uc * tmp_real * rt_auxvar%aqueous%dtotal(ires_mp,ires_uc,iphase) 
       else
        Jacobian(ires_up,ires_uc) = Jacobian(ires_up,ires_uc) - drate_uc * tmp_real !rate * tmp_real
        Jacobian(ires_mp, ires_uc) = Jacobian(ires_mp, ires_uc) + drate_uc * tmp_real !rate * tmp_real
       endif

       if(this%upstream_isaqueous_p) then
        Jacobian(ires_up,ires_up) = Jacobian(ires_up,ires_up) + &
                 drate_uc * rt_auxvar%aqueous%dtotal(ires_up,ires_up,iphase) 
        Jacobian(ires_mp,ires_up) = Jacobian(ires_mp,ires_up) - &
                 drate_uc * rt_auxvar%aqueous%dtotal(ires_mp,ires_up,iphase)
       else
        Jacobian(ires_up,ires_up) = Jacobian(ires_up,ires_up) + drate_uc !rate / & 
!                                    rt_auxvar%immobile(this%upstream_ispec_c)
        Jacobian(ires_mp,ires_up) = Jacobian(ires_mp,ires_up) - drate_uc !rate / &
!                                    rt_auxvar%immobile(this%upstream_ispec_c)
       endif
        do i = 1, this%nDownstream
           if(this%b_downstream_cp_follow_upstream(i)) then
               ires_dp = this%downstream_ires_p(i)

             if(this%upstream_isaqueous_c) then
              Jacobian(ires_dp,ires_uc) = Jacobian(ires_dp,ires_uc) + &
                       drate_uc * tmp_real * this%downstream_stoich_c(i) * &
                       rt_auxvar%aqueous%dtotal(ires_dp, ires_uc, iphase) 
              Jacobian(ires_mp,ires_uc) = Jacobian(ires_mp,ires_uc) - &
                       drate_uc * tmp_real * this%downstream_stoich_c(i) * &
                       rt_auxvar%aqueous%dtotal(ires_mp, ires_uc, iphase) 
             else
              Jacobian(ires_dp,ires_uc) = Jacobian(ires_dp,ires_uc) + &
                       drate_uc * tmp_real * this%downstream_stoich_c(i) 
              Jacobian(ires_mp,ires_uc) = Jacobian(ires_mp,ires_uc) - &
                       drate_uc * tmp_real * this%downstream_stoich_c(i)
             endif

             if(this%upstream_isaqueous_p) then
              Jacobian(ires_dp,ires_up) = Jacobian(ires_dp,ires_up) - &
                       drate_uc * this%downstream_stoich_c(i) * &
                       rt_auxvar%aqueous%dtotal(ires_dp, ires_up, iphase) 
              Jacobian(ires_mp,ires_up) = Jacobian(ires_mp,ires_up) + &
                       drate_uc * this%downstream_stoich_c(i) * & 
                       rt_auxvar%aqueous%dtotal(ires_mp, ires_up, iphase) 
             else
              Jacobian(ires_dp,ires_up) = Jacobian(ires_dp,ires_up) - &
                       drate_uc * this%downstream_stoich_c(i) 
              Jacobian(ires_mp,ires_up) = Jacobian(ires_mp,ires_up) + &
                       drate_uc * this%downstream_stoich_c(i)  
             endif
           endif
        enddo
     endif

     do i = 1, this%nDownstream
        if((.not.this%b_downstream_cp_fixed(i)) .and. &
           (.not.this%b_downstream_cn_follow_upstream(i))) then
          if(this%downstream_isaqueous_c(i)) then
             ires_dc = this%downstream_ispec_c(i)
             c_dc = rt_auxvar%total(this%downstream_ispec_c(i), iphase)
             c_dc = c_dc * Lwater_m3 
          else 
             ires_dc = this%downstream_ires_c(i) !reaction%offset_immobile + this%downstream_ispec_c(i)
             c_dc = rt_auxvar%immobile(this%downstream_ispec_c(i)) 
          endif

          if(this%downstream_isaqueous_p(i)) then
             ires_dp = this%downstream_ispec_p(i)
             c_dp = rt_auxvar%total(this%downstream_ispec_p(i), iphase)
             c_dp = c_dp * Lwater_m3 
          else 
             ires_dp = this%downstream_ires_p(i) !reaction%offset_immobile + this%downstream_ispec_p(i)
             c_dp = rt_auxvar%immobile(this%downstream_ispec_p(i)) 
          endif
          tmp_real = c_dp / c_dc / c_dc

          if(this%downstream_isaqueous_c(i)) then
           Jacobian(ires_dp,ires_dc) = Jacobian(ires_dp,ires_dc) + &
                    rate * tmp_real * this%downstream_stoich_c(i) * &
                    rt_auxvar%aqueous%dtotal(ires_dp, ires_dc, iphase) 
           Jacobian(ires_mp,ires_dc) = Jacobian(ires_mp,ires_dc) - &
                    rate * tmp_real * this%downstream_stoich_c(i) * &
                    rt_auxvar%aqueous%dtotal(ires_mp, ires_dc, iphase) 
          else
           Jacobian(ires_dp,ires_dc) = Jacobian(ires_dp,ires_dc) + &
                    rate * tmp_real * this%downstream_stoich_c(i)
           Jacobian(ires_mp,ires_dc) = Jacobian(ires_mp,ires_dc) - &
                    rate * tmp_real * this%downstream_stoich_c(i)
          endif

          if(this%downstream_isaqueous_p(i)) then
           Jacobian(ires_dp,ires_dp) = Jacobian(ires_dp,ires_dp) - &
                    rate / c_dc * this%downstream_stoich_c(i) * &
                    rt_auxvar%aqueous%dtotal(ires_dp, ires_dp, iphase) 
           Jacobian(ires_mp, ires_dp) = Jacobian(ires_mp, ires_dp) + &
                    rate / c_dc * this%downstream_stoich_c(i) * &
                    rt_auxvar%aqueous%dtotal(ires_mp, ires_dp, iphase)
          else
           Jacobian(ires_dp,ires_dp) = Jacobian(ires_dp,ires_dp) - &
                    rate / c_dc * this%downstream_stoich_c(i)
           Jacobian(ires_mp, ires_dp) = Jacobian(ires_mp, ires_dp) + &
                    rate / c_dc * this%downstream_stoich_c(i)
          endif 
        endif
     enddo 
  endif

!  first order terms
  do j = 1, this%nFirstOrder
     is_pool_in_rxn = PETSC_FALSE
     isaqueous = PETSC_FALSE
     if(this%ispec_1st(j) == this%upstream_ispec_c) then
       drate = drate_uc
       ires_j = ires_uc
       is_pool_in_rxn = PETSC_TRUE
       isaqueous = this%upstream_isaqueous_c
     else
       do i = 1, this%nDownstream
            if(this%ispec_1st(j) == this%downstream_ispec_c(i)) then
!               isaqueous = this%downstream_isaqueous_c(this%ispec_1st(j))
               isaqueous = this%downstream_isaqueous_c(i)
               if(isaqueous) then
                  ires_j = this%ispec_1st(j)
                  conc = rt_auxvar%total(this%ispec_1st(j), iphase)
                  conc = conc * Lwater_m3  
               else
                  ires_j = this%ispec_1st(j) + reaction%offset_immobile
                  conc = rt_auxvar%immobile(this%ispec_1st(j))  
               endif
               drate = rate / conc
               is_pool_in_rxn = PETSC_TRUE
               exit
            endif 
       enddo
     endif

     if(.not.is_pool_in_rxn) then
       cycle
     endif

! upstream C
     if(isaqueous) then
        Jacobian(ires_uc,ires_j) = Jacobian(ires_uc,ires_j) + drate * &
                 rt_auxvar%aqueous%dtotal(ires_uc,ires_j, iphase)
     else
        Jacobian(ires_uc,ires_j) = Jacobian(ires_uc,ires_j) + drate
     endif
       
! mineral C
     if(isaqueous) then
        Jacobian(ires_mc, ires_j) = Jacobian(ires_mc, ires_j) - &
                 drate * this%mineral_c_stoich * &
                 rt_auxvar%aqueous%dtotal(ires_mc,ires_j, iphase)
     else
        Jacobian(ires_mc, ires_j) = Jacobian(ires_mc, ires_j) - &
                 drate * this%mineral_c_stoich
     endif

! downstream C
     do i = 1, this%nDownstream
        ires_dc = this%downstream_ires_c(i)
       if(isaqueous) then
        Jacobian(ires_dc,ires_j) = Jacobian(ires_dc,ires_j) - &
                 drate * this%downstream_stoich_c(i) * &
                 rt_auxvar%aqueous%dtotal(ires_dc,ires_j,iphase)
       else
        Jacobian(ires_dc,ires_j) = Jacobian(ires_dc,ires_j) - &
                 drate * this%downstream_stoich_c(i)
       endif
     enddo 

     if(this%bNEnabled) then
! upstream N
        if(.not.this%bfixed_upstream_cn) then
          if(isaqueous) then
           Jacobian(ires_un,ires_j) = Jacobian(ires_un,ires_j) + &
                    drate * this%upstream_stoich_n * &
                    rt_auxvar%aqueous%dtotal(ires_un,ires_j,iphase)
          else
           Jacobian(ires_un,ires_j) = Jacobian(ires_un,ires_j) + &
                    drate * this%upstream_stoich_n      
          endif
        endif 
! mineral N
        if(isaqueous) then
          Jacobian(ires_mn,ires_j) = Jacobian(ires_mn,ires_j) - &
                   drate * this%mineral_n_stoich * &
                   rt_auxvar%aqueous%dtotal(ires_mn,ires_j,iphase)
        else
          Jacobian(ires_mn,ires_j) = Jacobian(ires_mn,ires_j) - &
                   drate * this%mineral_n_stoich
        endif
! downstream N      
        do i = 1, this%nDownstream
           if(.not.this%b_downstream_cn_fixed(i)) then
             ires_dn = this%downstream_ires_n(i)      
             if(isaqueous) then
                Jacobian(ires_dn,ires_j) = Jacobian(ires_dn,ires_j) - &
                      drate * this%downstream_stoich_n(i) * &
                      this%downstream_stoich_c(i) * &
                      rt_auxvar%aqueous%dtotal(ires_dn,ires_j,iphase)
             else
                Jacobian(ires_dn,ires_j) = Jacobian(ires_dn,ires_j) - &
                      drate * this%downstream_stoich_n(i) * &
                      this%downstream_stoich_c(i)
             endif
           endif
        enddo 
     endif

     if(this%bPEnabled) then
! upstream P
        if(.not.this%bfixed_upstream_cp) then
          if(isaqueous) then
             Jacobian(ires_up,ires_j) = Jacobian(ires_up,ires_j) + &
                      drate * this%upstream_stoich_p * &
                      rt_auxvar%aqueous%dtotal(ires_up,ires_j,iphase)
          else
             Jacobian(ires_up,ires_j) = Jacobian(ires_up,ires_j) + &
                      drate * this%upstream_stoich_p      
          endif
        endif 
! mineral P
        if(isaqueous) then
           Jacobian(ires_mp,ires_j) = Jacobian(ires_mp,ires_j) - &
                    drate * this%mineral_p_stoich * &
                    rt_auxvar%aqueous%dtotal(ires_mp,ires_j,iphase)
        else
           Jacobian(ires_mp,ires_j) = Jacobian(ires_mp,ires_j) - &
                    drate * this%mineral_p_stoich
        endif
! downstream P      
        do i = 1, this%nDownstream
           if(.not.this%b_downstream_cp_fixed(i)) then
             ires_dp = this%downstream_ires_p(i)
             if(isaqueous) then
                Jacobian(ires_dp,ires_j) = Jacobian(ires_dp,ires_j) - &
                         drate * this%downstream_stoich_p(i) * &
                         this%downstream_stoich_c(i) * &
                         rt_auxvar%aqueous%dtotal(ires_dp,ires_j,iphase)
             else
                Jacobian(ires_dp,ires_j) = Jacobian(ires_dp,ires_j) - &
                         drate * this%downstream_stoich_p(i) * &
                         this%downstream_stoich_c(i)
             endif
           endif
        enddo 
     endif
  enddo 

!  monod terms
  do j = 1, this%nMonod
     is_pool_in_rxn = PETSC_FALSE
     isaqueous = PETSC_FALSE
     if(this%ispec_mnd(j) == this%upstream_ispec_c) then
       is_pool_in_rxn = PETSC_TRUE
       conc = c_uc
       ires_j = ires_uc
       isaqueous = this%upstream_isaqueous_c
     else
       do i = 1, this%nDownstream
            if(this%ispec_mnd(j) == this%downstream_ispec_c(i)) then
               isaqueous = this%downstream_isaqueous_c(i)
               if(isaqueous) then
                   conc = rt_auxvar%total(this%ispec_mnd(j),iphase)
                   conc = conc * Lwater_m3  
               else
                   conc = rt_auxvar%immobile(this%ispec_mnd(j))  
               endif    
               is_pool_in_rxn = PETSC_TRUE
               ires_j = this%downstream_ires_c(j)
               exit
            endif 
       enddo
     endif

     if(.not.is_pool_in_rxn) then
        cycle
     endif

! r     = r0 * x / (k + x)     
! dr/dx = r0 * k / (k + x) / (k + x) = r * k / x / (k + x)
     drate = rate * this%half_saturation(j) / conc / (this%half_saturation(j) + conc) 

! upstream C
     if(isaqueous) then 
        Jacobian(ires_uc, ires_j) = Jacobian(ires_uc, ires_j) + drate * &
                 rt_auxvar%aqueous%dtotal(ires_uc,ires_j,iphase)
     else
        Jacobian(ires_uc, ires_j) = Jacobian(ires_uc, ires_j) + drate       
     endif

! mineral C
     if(isaqueous) then 
        Jacobian(ires_mc, ires_j) = Jacobian(ires_mc, ires_j) - &
                 drate * this%mineral_c_stoich * &
                 rt_auxvar%aqueous%dtotal(ires_mc, ires_j,iphase)
     else
        Jacobian(ires_mc, ires_j) = Jacobian(ires_mc, ires_j) - &
                 drate * this%mineral_c_stoich
     endif

! downstream C
     do i = 1, this%nDownstream
        ires_dc = this%downstream_ires_c(i)
       if(isaqueous) then
        Jacobian(ires_dc,ires_j) = Jacobian(ires_dc,ires_j) - &
                 drate * this%downstream_stoich_c(i) * &
                 rt_auxvar%aqueous%dtotal(ires_dc, ires_j,iphase)
       else
        Jacobian(ires_dc,ires_j) = Jacobian(ires_dc,ires_j) - &
                 drate * this%downstream_stoich_c(i)
       endif
     enddo 

     if(this%bNEnabled) then
! upstream N
       if(.not.this%bfixed_upstream_cn) then
         if(isaqueous) then
          Jacobian(ires_un,ires_j) = Jacobian(ires_un,ires_j) + &
                   drate * this%upstream_stoich_n * &
                   rt_auxvar%aqueous%dtotal(ires_un,ires_j,iphase)
         else
          Jacobian(ires_un,ires_j) = Jacobian(ires_un,ires_j) + &
                   drate * this%upstream_stoich_n
         endif
       endif

 ! mineral N
       if(isaqueous) then
          Jacobian(ires_mn,ires_j) = Jacobian(ires_mn,ires_j) - &
                   drate * this%mineral_n_stoich * &     
                   rt_auxvar%aqueous%dtotal(ires_mn,ires_j,iphase)
       else
          Jacobian(ires_mn,ires_j) = Jacobian(ires_mn,ires_j) - &
                   drate * this%mineral_n_stoich
       endif

! downstream N       
       do i = 1, this%nDownstream
        if(.not.this%b_downstream_cn_fixed(i)) then
          ires_dn = this%downstream_ires_n(i)
          if(isaqueous) then
             Jacobian(ires_dn,ires_j) = Jacobian(ires_dn,ires_j) - &
                      drate * this%downstream_stoich_n(i) * &
                      this%downstream_stoich_c(i) * & 
                      rt_auxvar%aqueous%dtotal(ires_dn,ires_j,iphase)
          else
             Jacobian(ires_dn,ires_j) = Jacobian(ires_dn,ires_j) - &
                      drate * this%downstream_stoich_n(i) * &
                      this%downstream_stoich_c(i)      
          endif
        endif
       enddo 
     endif

     if(this%bPEnabled) then
! upstream P
       if(.not.this%bfixed_upstream_cp) then
          if(isaqueous) then
             Jacobian(ires_up,ires_j) = Jacobian(ires_up,ires_j) + &
                      drate * this%upstream_stoich_p * &
                      rt_auxvar%aqueous%dtotal(ires_up,ires_j,iphase)
          else
             Jacobian(ires_up,ires_j) = Jacobian(ires_up,ires_j) + &
                      drate * this%upstream_stoich_p
          endif
       endif

 ! mineral P
       if(isaqueous) then
          Jacobian(ires_mp,ires_j) = Jacobian(ires_mp,ires_j) - &
                   drate * this%mineral_p_stoich * &
                   rt_auxvar%aqueous%dtotal(ires_mp,ires_j,iphase)
       else
          Jacobian(ires_mp,ires_j) = Jacobian(ires_mp,ires_j) - &
                   drate * this%mineral_p_stoich      
       endif

! downstream P       
       do i = 1, this%nDownstream
        if(.not.this%b_downstream_cp_fixed(i)) then
          ires_dp = this%downstream_ires_p(i)
          if(isaqueous) then
             Jacobian(ires_dp,ires_j) = Jacobian(ires_dp,ires_j) - &
                   drate * this%downstream_stoich_p(i) * &     
                   this%downstream_stoich_c(i) * &
                   rt_auxvar%aqueous%dtotal(ires_dp,ires_j,iphase)
          else
             Jacobian(ires_dp,ires_j) = Jacobian(ires_dp,ires_j) - &
                   drate * this%downstream_stoich_p(i) * &     
                   this%downstream_stoich_c(i)      
          endif
        endif
       enddo 
     endif
   enddo 

!  inhibition terms
   do j = 1, this%nInhibition
     is_pool_in_rxn = PETSC_FALSE
     isaqueous = PETSC_FALSE
     if(this%ispec_inh(j) == this%upstream_ispec_c) then
       conc = c_uc
       ires_j = ires_uc
       is_pool_in_rxn = PETSC_TRUE
       isaqueous = this%upstream_isaqueous_c
     else
       do i = 1, this%nDownstream
            if(this%ispec_inh(j) == this%downstream_ispec_c(i)) then
               isaqueous = this%downstream_isaqueous_c(i)
               if(isaqueous) then
                   conc = rt_auxvar%total(this%ispec_inh(j),iphase)
                   conc = conc * Lwater_m3  
               else
                   conc = rt_auxvar%immobile(this%ispec_inh(j))  
               endif    
               ires_j = this%downstream_ires_c(i)
               is_pool_in_rxn = PETSC_TRUE
               exit
            endif 
       enddo
     endif

     if(.not.is_pool_in_rxn) then
       cycle !continue
     endif

! r     = r0 * k / (k + x)     
! dr/dx = - r0 * k / (k + x) / (k + x) = - r / (k + x)
     drate = -rate / (this%inhibition_coef(j) + conc) 

! upstream C
     if(isaqueous) then
        Jacobian(ires_uc, ires_j) = Jacobian(ires_uc, ires_j) + drate * &
                 rt_auxvar%aqueous%dtotal(ires_uc,ires_j,iphase)
     else
        Jacobian(ires_uc, ires_j) = Jacobian(ires_uc, ires_j) + drate       
     endif

! mineral C
     if(isaqueous) then
        Jacobian(ires_mc, ires_j) = Jacobian(ires_mc, ires_j) - &
                 drate * this%mineral_c_stoich * &     
                 rt_auxvar%aqueous%dtotal(ires_mc,ires_j,iphase)
     else
        Jacobian(ires_mc, ires_j) = Jacobian(ires_mc, ires_j) - &
                 drate * this%mineral_c_stoich      
     endif

! downstream C
     do i = 1, this%nDownstream
        ires_dc = this%downstream_ires_c(i)
        if(isaqueous) then
           Jacobian(ires_dc,ires_j) = Jacobian(ires_dc,ires_j) - &
                    drate * this%downstream_stoich_c(i) * &
                    rt_auxvar%aqueous%dtotal(ires_dc,ires_j,iphase)
        else
           Jacobian(ires_dc,ires_j) = Jacobian(ires_dc,ires_j) - &
                    drate * this%downstream_stoich_c(i)
        endif
     enddo 

     if(this%bNEnabled) then
! upstream N
       if(.not.this%bfixed_upstream_cn) then
         if(isaqueous) then
          Jacobian(ires_un,ires_j) = Jacobian(ires_un,ires_j) + &
                   drate * this%upstream_stoich_n * &
                   rt_auxvar%aqueous%dtotal(ires_un,ires_j,iphase)
         else
          Jacobian(ires_un,ires_j) = Jacobian(ires_un,ires_j) + &
                   drate * this%upstream_stoich_n      
         endif
       endif

! mineral N 
       if(isaqueous) then
          Jacobian(ires_mn,ires_j) = Jacobian(ires_mn,ires_j) - &
                   drate * this%mineral_n_stoich * &     
                   rt_auxvar%aqueous%dtotal(ires_mn,ires_j,iphase)
       else
          Jacobian(ires_mn,ires_j) = Jacobian(ires_mn,ires_j) - &
                   drate * this%mineral_n_stoich      
       endif

! downstream N       
       do i = 1, this%nDownstream
        if(.not.this%b_downstream_cn_fixed(i)) then
          ires_dn = this%downstream_ires_n(i)
          if(isaqueous) then
             Jacobian(ires_dn,ires_j) = Jacobian(ires_dn,ires_j) - &
                      drate * this%downstream_stoich_n(i) * &
                      this%downstream_stoich_c(i) * &      
                      rt_auxvar%aqueous%dtotal(ires_dn,ires_j,iphase)
          else
             Jacobian(ires_dn,ires_j) = Jacobian(ires_dn,ires_j) - &
                      drate * this%downstream_stoich_n(i) * &
                      this%downstream_stoich_c(i)      
          endif
        endif
       enddo 
     endif

     if(this%bPEnabled) then
! upstream P
       if(.not.this%bfixed_upstream_cp) then
         if(isaqueous) then
          Jacobian(ires_up,ires_j) = Jacobian(ires_up,ires_j) + &
                   drate * this%upstream_stoich_p * &     
                   rt_auxvar%aqueous%dtotal(ires_up,ires_j,iphase)
         else
          Jacobian(ires_up,ires_j) = Jacobian(ires_up,ires_j) + &
                   drate * this%upstream_stoich_p      
         endif
       endif

! mineral P 
       if(isaqueous) then
         Jacobian(ires_mp,ires_j) = Jacobian(ires_mp,ires_j) - &
                  drate * this%mineral_p_stoich * &     
                  rt_auxvar%aqueous%dtotal(ires_mp,ires_j,iphase)
       else
         Jacobian(ires_mp,ires_j) = Jacobian(ires_mp,ires_j) - &
                drate * this%mineral_p_stoich      
       endif

! downstream P       
       do i = 1, this%nDownstream
        if(.not.this%b_downstream_cp_fixed(i)) then
          ires_dp = this%downstream_ires_p(i)
          Jacobian(ires_dp,ires_j) = Jacobian(ires_dp,ires_j) - drate * this%downstream_stoich_p(i) * &      
                                  this%downstream_stoich_c(i)      
        endif
       enddo 
     endif
   enddo 

! N limiting
   if(this%bNEnabled .and. this%mineral_n_stoich < 0.0d0 .and. this%NInhibitionCoef > 0.0d0) then
! r     = r0 * x / (k + x)     
! dr/dx = r0 * k / (k + x) / (k + x) = r * k / x / (k + x)
     conc = c_mn
     drate = drate_n * this%NInhibitionCoef / (this%NInhibitionCoef + conc)

     if(this%isaqueous_mn) then
! upstream C
       Jacobian(ires_uc, ires_mn) = Jacobian(ires_uc, ires_mn) + drate * &
                rt_auxvar%aqueous%dtotal(ires_uc,ires_mn,iphase)

! mineral C
       Jacobian(ires_mc, ires_mn) = Jacobian(ires_mc, ires_mn) - &
                drate * this%mineral_c_stoich * &
                rt_auxvar%aqueous%dtotal(ires_mc,ires_mn,iphase)

! downstream C
       do i = 1, this%nDownstream
        ires_dc = this%downstream_ires_c(i)
        Jacobian(ires_dc,ires_mn) = Jacobian(ires_dc,ires_mn) - &
                 drate * this%downstream_stoich_c(i) * &
                 rt_auxvar%aqueous%dtotal(ires_dc,ires_mn,iphase)
       enddo 

       if(.not.this%bfixed_upstream_cn) then
! upstream N
         Jacobian(ires_un,ires_mn) = Jacobian(ires_un,ires_mn) + &
                  drate * this%upstream_stoich_n * &     
                  rt_auxvar%aqueous%dtotal(ires_un,ires_mn,iphase)

! mineral N       
         Jacobian(ires_mn,ires_mn) = Jacobian(ires_mn,ires_mn) - &
                  drate * this%mineral_n_stoich * &
                  rt_auxvar%aqueous%dtotal(ires_mn,ires_mn,iphase)
       endif 

! downstream N
       do i = 1, this%nDownstream
        if(.not.this%b_downstream_cn_fixed(i)) then
          ires_dn = this%downstream_ires_n(i) !reaction%offset_immobile + this%downstream_ispec_n(i)      
          Jacobian(ires_dn,ires_mn) = Jacobian(ires_dn,ires_mn) - &
                   drate * this%downstream_stoich_n(i) * &
                   this%downstream_stoich_c(i) * &
                   rt_auxvar%aqueous%dtotal(ires_dn,ires_mn,iphase)
        endif
       enddo 
     
       if(this%bPEnabled) then
! upstream P
        if(.not.this%bfixed_upstream_cp) then
          Jacobian(ires_up,ires_mn) = Jacobian(ires_up,ires_mn) + &
                   drate * this%upstream_stoich_p * &     
                   rt_auxvar%aqueous%dtotal(ires_up,ires_mn,iphase)

! mineral P       
          Jacobian(ires_mp,ires_mn) = Jacobian(ires_mp,ires_mn) - &
                   drate * this%mineral_p_stoich * &
                   rt_auxvar%aqueous%dtotal(ires_mp,ires_mn,iphase)
        endif 

! downstream P
        do i = 1, this%nDownstream
           if(.not.this%b_downstream_cp_fixed(i)) then
!               ires = reaction%offset_immobile + this%downstream_ispec_p(i)
               ires_dp = this%downstream_ires_p(i)
               Jacobian(ires_dp,ires_mn) = Jacobian(ires_dp,ires_mn) - &
                        drate * this%downstream_stoich_p(i) * &
                        this%downstream_stoich_c(i) * &
                        rt_auxvar%aqueous%dtotal(ires_dp,ires_mn,iphase)
           endif
        enddo 
       endif
     else  ! aqueous
! upstream C
       Jacobian(ires_uc, ires_mn) = Jacobian(ires_uc, ires_mn) + drate

! mineral C
       Jacobian(ires_mc, ires_mn) = Jacobian(ires_mc, ires_mn) - &
                drate * this%mineral_c_stoich

! downstream C
       do i = 1, this%nDownstream
        ires_dc = this%downstream_ires_c(i)
        Jacobian(ires_dc,ires_mn) = Jacobian(ires_dc,ires_mn) - &
                 drate * this%downstream_stoich_c(i)
       enddo 

       if(.not.this%bfixed_upstream_cn) then
! upstream N
         Jacobian(ires_un,ires_mn) = Jacobian(ires_un,ires_mn) + &
                  drate * this%upstream_stoich_n

! mineral N       
         Jacobian(ires_mn,ires_mn) = Jacobian(ires_mn,ires_mn) - &
                  drate * this%mineral_n_stoich
       endif 

! downstream N
       do i = 1, this%nDownstream
        if(.not.this%b_downstream_cn_fixed(i)) then
          ires_dn = this%downstream_ires_n(i) !reaction%offset_immobile + this%downstream_ispec_n(i)      
          Jacobian(ires_dn,ires_mn) = Jacobian(ires_dn,ires_mn) - &
                   drate * this%downstream_stoich_n(i) * &
                   this%downstream_stoich_c(i)
        endif
       enddo 
     
       if(this%bPEnabled) then
! upstream P
        if(.not.this%bfixed_upstream_cp) then
          Jacobian(ires_up,ires_mn) = Jacobian(ires_up,ires_mn) + &
                   drate * this%upstream_stoich_p

! mineral P       
          Jacobian(ires_mp,ires_mn) = Jacobian(ires_mp,ires_mn) - &
                   drate * this%mineral_p_stoich
        endif 

! downstream P
        do i = 1, this%nDownstream
           if(.not.this%b_downstream_cp_fixed(i)) then
               ires_dp = this%downstream_ires_p(i)
               Jacobian(ires_dp,ires_mn) = Jacobian(ires_dp,ires_mn) - &
                        drate * this%downstream_stoich_p(i) * &
                        this%downstream_stoich_c(i)
           endif
        enddo 
       endif
     endif  ! if aqueous
   endif  ! if N limiting

! P limiting
   if(this%bPEnabled .and. this%mineral_p_stoich < 0.0d0 .and. this%PInhibitionCoef > 0.0d0) then
! r     = r0 * x / (k + x)     
! dr/dx = r0 * k / (k + x) / (k + x) = r * k / x / (k + x)
     conc = c_mp
     drate = drate_p * this%PInhibitionCoef / (this%PInhibitionCoef + conc)

     if(this%isaqueous_mp) then
! upstream C
        Jacobian(ires_uc, ires_mp) = Jacobian(ires_uc, ires_mp) + drate * &
                rt_auxvar%aqueous%dtotal(ires_uc,ires_mp,iphase)

! mineral C
!     ires_c = reaction%offset_immobile + this%ispec_mc      
        Jacobian(ires_mc, ires_mp) = Jacobian(ires_mc, ires_mp) - &
                 drate * this%mineral_c_stoich * &
                 rt_auxvar%aqueous%dtotal(ires_mc,ires_mp,iphase)

! downstream C
        do i = 1, this%nDownstream
!        ires = reaction%offset_immobile + this%downstream_ispec_c(i)      
           ires_dc = this%downstream_ires_c(i)
           Jacobian(ires_dc,ires_mp) = Jacobian(ires_dc,ires_mp) - &
                    drate * this%downstream_stoich_c(i) * &
                    rt_auxvar%aqueous%dtotal(ires_dc,ires_mp,iphase)
        enddo 

        if(.not.this%bfixed_upstream_cp) then
! upstream P
           Jacobian(ires_up,ires_mp) = Jacobian(ires_up,ires_mp) + &
                    drate * this%upstream_stoich_p * &
                    rt_auxvar%aqueous%dtotal(ires_up,ires_mp,iphase)

! mineral P       
           Jacobian(ires_mp,ires_mp) = Jacobian(ires_mp,ires_mp) - &
                    drate * this%mineral_p_stoich * &
                    rt_auxvar%aqueous%dtotal(ires_mp,ires_mp,iphase)
        endif 

! downstream P
        do i = 1, this%nDownstream
           if(.not.this%b_downstream_cp_fixed(i)) then
             ires_dp = this%downstream_ires_p(i)
             Jacobian(ires_dp,ires_mp) = Jacobian(ires_dp,ires_mp) - &
                      drate * this%downstream_stoich_p(i) * &
                      this%downstream_stoich_c(i) * &
                      rt_auxvar%aqueous%dtotal(ires_dp,ires_mp,iphase)
           endif
        enddo 

        if(this%bNEnabled) then
          if(.not.this%bfixed_upstream_cn) then
! upstream N
             Jacobian(ires_un,ires_mp) = Jacobian(ires_un,ires_mp) + &
                      drate * this%upstream_stoich_n * &     
                      rt_auxvar%aqueous%dtotal(ires_un,ires_mp,iphase)

! mineral N       
             Jacobian(ires_mn,ires_mp) = Jacobian(ires_mn,ires_mp) - &
                      drate * this%mineral_n_stoich * &
                      rt_auxvar%aqueous%dtotal(ires_mn,ires_mp,iphase)
          endif 

! downstream N
           do i = 1, this%nDownstream
              if(.not.this%b_downstream_cn_fixed(i)) then
                ires_dn = this%downstream_ires_n(i)
!             ires = reaction%offset_immobile + this%downstream_ispec_n(i)      
                Jacobian(ires_dn,ires_mp) = Jacobian(ires_dn,ires_mp) - &
                         drate * this%downstream_stoich_n(i) * &
                         this%downstream_stoich_c(i) * &
                         rt_auxvar%aqueous%dtotal(ires_dn,ires_mp,iphase)
              endif
           enddo 
        endif
      else  !if aqueous
        Jacobian(ires_uc, ires_mp) = Jacobian(ires_uc, ires_mp) + drate

! mineral C
        Jacobian(ires_mc, ires_mp) = Jacobian(ires_mc, ires_mp) - &
                 drate * this%mineral_c_stoich

! downstream C
        do i = 1, this%nDownstream
           ires_dc = this%downstream_ires_c(i)
           Jacobian(ires_dc,ires_mp) = Jacobian(ires_dc,ires_mp) - &
                    drate * this%downstream_stoich_c(i)
        enddo 

        if(.not.this%bfixed_upstream_cp) then
! upstream P
           Jacobian(ires_up,ires_mp) = Jacobian(ires_up,ires_mp) + &
                    drate * this%upstream_stoich_p      

! mineral P       
           Jacobian(ires_mp,ires_mp) = Jacobian(ires_mp,ires_mp) - &
                    drate * this%mineral_p_stoich      
        endif 

! downstream P
        do i = 1, this%nDownstream
           if(.not.this%b_downstream_cp_fixed(i)) then
             ires_dp = this%downstream_ires_p(i)
             Jacobian(ires_dp,ires_mp) = Jacobian(ires_dp,ires_mp) - &
                      drate * this%downstream_stoich_p(i) * &
                      this%downstream_stoich_c(i)
           endif
        enddo 

        if(this%bNEnabled) then
          if(.not.this%bfixed_upstream_cn) then
! upstream N
             Jacobian(ires_un,ires_mp) = Jacobian(ires_un,ires_mp) + &
                      drate * this%upstream_stoich_n      

! mineral N       
             Jacobian(ires_mn,ires_mp) = Jacobian(ires_mn,ires_mp) - &
                      drate * this%mineral_n_stoich      
          endif 

! downstream N
           do i = 1, this%nDownstream
              if(.not.this%b_downstream_cn_fixed(i)) then
                ires_dn = this%downstream_ires_n(i)
!             ires = reaction%offset_immobile + this%downstream_ispec_n(i)      
                Jacobian(ires_dn,ires_mp) = Jacobian(ires_dn,ires_mp) - &
                         drate * this%downstream_stoich_n(i) * &
                         this%downstream_stoich_c(i)
              endif
           enddo 
        endif
      endif !if aqueous
   endif  ! if P limiting
end subroutine CLM_CN_BFReact

! ************************************************************************** !
!
! CLM_CN_BF_Update_Stoich: Update stoichiometric coefficients for CLM_CN_Mirobe  
! author: Guoping Tang
! date: 00/00/00
! added for CLM-CN-MICROBE  11/20/2013 
!
! ************************************************************************** !
subroutine CLM_CN_BF_Update_Stoich(this,ispec_bacteria,ispec_fungi) 

  implicit none
  
  class(reaction_sandbox_CLM_CN_BF_type) :: this  

  PetscReal :: respiration_fraction, c_uc, c_un
  PetscInt :: i, ispec_bacteria, ispec_fungi
  character(len=MAXWORDLENGTH) :: word

! upstream pool is litter
  if(this%is_upstream_litter) then
!    Sinsabaugh et al. 2013 Ecology Letters, 16, 930-939
     respiration_fraction = CN_ratio_microbe * this%upstream_stoich_n !c_un/c_uc

     if(respiration_fraction > CUE_max) then
        respiration_fraction = CUE_max
     endif

!    c pools     
     this%mineral_c_stoich = respiration_fraction

     do i = 1, this%nDownstream
        if(this%downstream_ispec_c(i) == ispec_bacteria) then
           this%downstream_stoich_c(i) = fraction_bacteria * (1.0d0  - respiration_fraction)
        else
           this%downstream_stoich_c(i) = (1.0d0 - fraction_bacteria) * (1.0d0  - respiration_fraction) 
        endif
     enddo
  endif

end subroutine CLM_CN_BF_Update_Stoich 
!-------------------------------------------------------------------------------

! ************************************************************************** !
!
! CLM_CN_BFDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Guoping Tang
! date: 00/00/00
!
! ************************************************************************** !
subroutine CLM_CN_BFDestroy(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_CLM_CN_BF_type) :: this  

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

  call DeallocateArray(this%downstream_ires_c) 
  call DeallocateArray(this%downstream_ires_n) 
  call DeallocateArray(this%downstream_ires_p) 

  call DeallocateArray(this%downstream_isaqueous_c) 
  call DeallocateArray(this%downstream_isaqueous_n) 
  call DeallocateArray(this%downstream_isaqueous_p) 

  call DeallocateArray(this%downstream_stoich_c) 
  call DeallocateArray(this%downstream_stoich_n) 
  call DeallocateArray(this%downstream_stoich_p) 

!  call DeallocateArray(this%downstream_cn) 
!  call DeallocateArray(this%downstream_cp)

  call DeallocateArray(this%b_downstream_cn_fixed) 
  call DeallocateArray(this%b_downstream_cp_fixed)

  call DeallocateArray(this%b_downstream_cn_follow_upstream) 
  call DeallocateArray(this%b_downstream_cp_follow_upstream)

  call DeallocateArray(this%ispec_1st) 
  call DeallocateArray(this%ispec_mnd) 
  call DeallocateArray(this%ispec_inh) 

  call DeallocateArray(this%isaqueous_1st) 
  call DeallocateArray(this%isaqueous_mnd) 
  call DeallocateArray(this%isaqueous_inh) 

  call DeallocateArray(this%half_saturation) 
  call DeallocateArray(this%inhibition_coef) 
 
end subroutine CLM_CN_BFDestroy

end module Reaction_Sandbox_CLM_CN_BF_class
