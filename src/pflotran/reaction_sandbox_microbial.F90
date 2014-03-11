module Reaction_Sandbox_Microbial_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  use Database_Aux_module

#ifdef CLM_PFLOTRAN
  use CLM_BGC_module
#endif

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public :: rate_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal                    :: value
    type(rate_type), pointer     :: next
  end type rate_type

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_Microbial_type

    PetscInt :: temperature_response_function
    PetscReal :: Q10
    PetscInt :: moisture_response_function
    PetscInt :: ph_response_function
    PetscReal :: fixed_ph

    character(len=MAXSTRINGLENGTH) :: str_reaction
    type(database_rxn_type), pointer ::dbase_rxn

! rate terms

    PetscReal :: rate_constant
    PetscReal :: x0eps

    character(len=MAXWORDLENGTH) :: electron_donor_name
    PetscReal :: electron_donor_half_saturation
    PetscInt :: electron_donor_ispec

    character(len=MAXWORDLENGTH) :: electron_acceptor_name
    PetscReal :: electron_acceptor_half_saturation
    PetscInt :: electron_acceptor_ispec

    character(len=MAXWORDLENGTH) :: microbial_mass_name
    PetscInt :: microbial_mass_ispec

    PetscInt :: nInhibition

    type(rate_type), pointer :: Inhibition

    PetscInt, pointer :: ispec_inh(:)
    PetscReal, pointer :: inhibition_coef(:)

  contains
    procedure, public :: ReadInput => MicrobialRead
    procedure, public :: Setup => MicrobialSetup
    procedure, public :: Evaluate => MicrobialReact
    procedure, public :: Destroy => MicrobialDestroy
  end type reaction_sandbox_Microbial_type

  public :: MicrobialCreate

contains

! ************************************************************************** !
!
! MicrobialCreate: Allocates Microbial reaction object.
! author: Guoping Tang
! date: 11/25/13
!
! ************************************************************************** !
function MicrobialCreate()

  implicit none
  
  class(reaction_sandbox_Microbial_type), pointer :: MicrobialCreate

  allocate(MicrobialCreate)

#ifdef CLM_PFLOTRAN
  MicrobialCreate%temperature_response_function = -1
  MicrobialCreate%moisture_response_function = -1
  MicrobialCreate%ph_response_function = -1
#endif

  MicrobialCreate%Q10 = 1.5d0

  MicrobialCreate%str_reaction = ''

  nullify(MicrobialCreate%dbase_rxn)

  MicrobialCreate%rate_constant = 0.d0
  MicrobialCreate%x0eps = 1.d-20
  MicrobialCreate%fixed_ph = -1.0d0

  MicrobialCreate%electron_donor_name = ''
  MicrobialCreate%electron_donor_half_saturation = 1.0d-6
  MicrobialCreate%electron_donor_ispec = -1

  MicrobialCreate%electron_acceptor_name = ''
  MicrobialCreate%electron_acceptor_half_saturation = 1.0d-6
  MicrobialCreate%electron_acceptor_ispec = -1

  MicrobialCreate%microbial_mass_name = ''
  MicrobialCreate%microbial_mass_ispec = -1

  MicrobialCreate%nInhibition = 0

  nullify(MicrobialCreate%Inhibition)
  nullify(MicrobialCreate%ispec_inh)
  nullify(MicrobialCreate%inhibition_coef)
  nullify(MicrobialCreate%next)  

end function MicrobialCreate

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
! MicrobialRead:
! author: Guoping Tang
! date: 11/25/13
!
! ************************************************************************** !
subroutine MicrobialRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_Microbial_type) :: this
  type(input_type)                     :: input
  type(option_type)                    :: option

  type(rate_type), pointer :: firstorder,firstorder_prev
  type(rate_type), pointer :: monod, monod_prev
  type(rate_type), pointer :: inhibition, inhibition_prev

  PetscReal :: tmp_real, rate_constant, turnover_time

  character(len=MAXWORDLENGTH) :: word

  nullify(firstorder_prev)
  nullify(monod_prev)
  nullify(inhibition_prev)

  turnover_time = 0.d0
  rate_constant = 0.d0

  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,Microbial')
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
                       'CHEMISTRY,REACTION_SANDBOX,Microbial,TEMPERATURE RESPONSE FUNCTION')
         call StringToUpper(word)   

            select case(trim(word))
              case('CLM4')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4    
              case('Q10')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_Q10    
                  call InputReadDouble(input,option,tmp_real)  
                  call InputErrorMsg(input,option,'Q10', &
                        'CHEMISTRY,REACTION_SANDBOX_Microbial,TEMPERATURE RESPONSE FUNCTION')
                  this%Q10 = tmp_real
              case('DLEM')
                  this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_DLEM    
                  call InputReadDouble(input,option,tmp_real)  
                  call InputErrorMsg(input,option,'Q10', &
                        'CHEMISTRY,REACTION_SANDBOX_Microbial,TEMPERATURE RESPONSE FUNCTION')
                  this%Q10 = tmp_real
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,Microbial,TEMPERATURE RESPONSE FUNCTION keyword: ' // &
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
                       'CHEMISTRY,REACTION_SANDBOX,Microbial,MOISTURE RESPONSE FUNCTION')
         call StringToUpper(word)   

            select case(trim(word))
              case('CLM4')
                  this%moisture_response_function = MOISTURE_RESPONSE_FUNCTION_CLM4    
              case('DLEM')
                  this%moisture_response_function = MOISTURE_RESPONSE_FUNCTION_DLEM    
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,Microbial,TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
         enddo 

      case('PH_RESPONSE_FUNCTION')
        do
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit

         call InputReadWord(input,option,word,PETSC_TRUE)
         call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,Microbial,PH RESPONSE FUNCTION')
         call StringToUpper(word)   

            select case(trim(word))
              case('CENTURY')
                  this%ph_response_function = PH_RESPONSE_FUNCTION_CENTURY    
              case('DLEM')
                  this%ph_response_function = PH_RESPONSE_FUNCTION_DLEM    
              case default
                  option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,Microbial,TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                                     trim(word) // ' not recognized.'
                  call printErrMsg(option)
            end select
         enddo 

      case('FIXED_PH')
        call InputReadDouble(input,option,this%fixed_ph)
        call InputErrorMsg(input,option,'fixed ph', &
               'CHEMISTRY,REACTION_SANDBOX,Microbial,RATE_CONSTANT')
#endif

      case('REACTION')
        ! remainder of string should be the reaction equation
        this%str_reaction = trim(adjustl(input%buf))
        ! set flag for error message
        if (len_trim(this%str_reaction) < 2) input%ierr = 1
        call InputErrorMsg(input,option,'reaction string', &
                            'CHEMISTRY,REACTION_SANDBOX_MICROBIAL_REACTION,REACTION')     

      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'rate constant', &
               'CHEMISTRY,REACTION_SANDBOX,Microbial,RATE_CONSTANT')
!        call InputReadWord(input,option,word,PETSC_TRUE)
!        if (InputError(input)) then
!            input%err_buf = 'Microbial RATE CONSTANT UNITS'
!            call InputDefaultMsg(input,option)
!        else              
!            this%rate_constant = rate_constant * &
!            UnitsConvertToInternal(word,option)
!            write(*, *) this%rate_constant, rate_constant
!        endif

      case('X0EPS')
        call InputReadDouble(input,option,this%x0eps)
        call InputErrorMsg(input,option,'x0eps', &
               'CHEMISTRY,REACTION_SANDBOX,Microbial,X0EPS')
      case('ELECTRON_DONOR')
        call InputReadWord(input,option,this%electron_donor_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'Electron donor species name', &
                           'CHEMISTRY,REACTION_SANDBOX_Microbial,ELECTRON_DONOR')
        call InputReadDouble(input,option,this%electron_donor_half_saturation)  
        call InputErrorMsg(input,option,'Electron donor half saturation constant', &
                           'CHEMISTRY,REACTION_SANDBOX_Microbial,ELECTRON_DONOR')

      case('ELECTRON_ACCEPTOR')
        call InputReadWord(input,option,this%electron_acceptor_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'Electron acceptor species name', &
                           'CHEMISTRY,REACTION_SANDBOX_Microbial,ELECTRON_ACCEPTOR')
        call InputReadDouble(input,option,this%electron_acceptor_half_saturation)  
        call InputErrorMsg(input,option,'Electron acceptor half saturation constant', &
                           'CHEMISTRY,REACTION_SANDBOX_Microbial,ELECTRON_ACCEPTOR')
      
      case('MICROBIAL_MASS')
        call InputReadWord(input,option,this%microbial_mass_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'Microbial mass species name', &
                           'CHEMISTRY,REACTION_SANDBOX_Microbial,MICROBIAL_MASS')

      case('INHIBITION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'species name', &
                           'CHEMISTRY,REACTION_SANDBOX_Microbial,INHIBITION')
        call InputReadDouble(input,option,tmp_real)  
        call InputErrorMsg(input,option,'inhibition constant', &
                           'CHEMISTRY,REACTION_SANDBOX_Microbial,INHIBITION')
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

      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,Microbial keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo

end subroutine MicrobialRead

! ************************************************************************** !
!
! MicrobialSetup: 
! author: Guoping Tang
! date: 11/25/13
!
! ************************************************************************** !
subroutine MicrobialSetup(this,reaction,option)

  use Option_module
  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Immobile_Aux_module
  use Database_Aux_module

  implicit none

  class(reaction_sandbox_Microbial_type) :: this
  type(reaction_type)                  :: reaction
  type(option_type)                    :: option

  type(rate_type), pointer             :: cur_rate

  PetscInt :: i, icount, ispec
  PetscReal :: sum_stoich_c_prod

  character(len=MAXWORDLENGTH) :: word

! parse reaction
  this%dbase_rxn => DatabaseRxnCreateFromRxnString(this%str_reaction, &
                                     reaction%naqcomp, &
                                     reaction%offset_aqueous, &
                                     reaction%primary_species_names, &
                                     reaction%nimcomp, &
                                     reaction%offset_immobile, &
                                     reaction%immobile%names, &
                                     option)

!  write(*, *) this%str_reaction
!  write(*, *) this%dbase_rxn%nspec
!  do i = 1, this%dbase_rxn%nspec
!     write(*, *) this%dbase_rxn%spec_name(i), this%dbase_rxn%spec_ids(i), &
!                 this%dbase_rxn%stoich(i) 
!  enddo

! electron donor
  if(trim(this%electron_donor_name) /= '')then
    this%electron_donor_ispec = GetPrimarySpeciesIDFromName(  &
         this%electron_donor_name,reaction,PETSC_FALSE,option)
!     write(*,*) 'electron donor', this%electron_donor_name, this%electron_donor_ispec
  endif

  if(trim(this%electron_acceptor_name) /= '')then
    this%electron_acceptor_ispec = GetPrimarySpeciesIDFromName(  &
         this%electron_acceptor_name,reaction,PETSC_FALSE,option)
!     write(*,*) 'electron acceptor', this%electron_acceptor_name, this%electron_acceptor_ispec
  endif

  if(trim(this%microbial_mass_name) /= '')then
!       this%microbial_mass_ispec = GetPrimarySpeciesIDFromName(  &
!               this%microbial_mass_name,reaction,PETSC_FALSE,option)
    this%microbial_mass_ispec = GetImmobileSpeciesIDFromName( &
         this%microbial_mass_name,reaction%immobile,PETSC_FALSE,option) 
  endif

! inhibition rate terms
  if(this%nInhibition >= 1) then
     allocate(this%ispec_inh(this%nInhibition))
     allocate(this%inhibition_coef(this%nInhibition))
 
     cur_rate => this%Inhibition
     icount = 1
     do
        if(.not.associated(cur_rate)) exit
        ispec = GetPrimarySpeciesIDFromName(cur_rate%name,reaction,PETSC_FALSE,option)
        if(ispec > 0) then
           this%ispec_inh(icount) = ispec
        else
           option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,Microbial check: ' // &
            'Inhibition term ' // trim(cur_rate%name)// &
           ' is not an aqueous species.'
           call printErrMsg(option)
        endif
        this%inhibition_coef(icount) = cur_rate%value
        cur_rate => cur_rate%next
        icount = icount + 1
     enddo 
  endif   

end subroutine MicrobialSetup

! ************************************************************************** !
!
! MicrobialReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 11/25/13
!
! ************************************************************************** !
subroutine MicrobialReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,porosity,volume,reaction, &
                         option,local_id)

  use Option_module
  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  
#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif
  
  implicit none
  
#ifdef CLM_PFLOTRAN
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#endif
  
  class(reaction_sandbox_Microbial_type) :: this  
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
  PetscReal :: Lwater
  PetscReal :: rate, drate, drate_donor, drate_acceptor, drate_biomass
  PetscReal :: tmp_real, conc, activity 

  PetscReal :: f_t
  PetscReal :: f_w
  PetscReal :: f_ph, ph
#ifdef CLM_PFLOTRAN
  PetscReal :: tc, theta
#endif

  PetscInt :: local_id, i, icomp, ires
  PetscErrorCode :: ierr

! temperature response function 
#ifdef CLM_PFLOTRAN 
  tc = global_auxvar%temp(1) 
  f_t = GetTemperatureResponse(tc, this%temperature_response_function, this%Q10) 
#else
  f_t = 1.0d0
#endif

  ! moisture response function 
#ifdef CLM_PFLOTRAN
  theta = global_auxvar%sat(1) * porosity 
  f_w = GetMoistureResponse(theta, local_id, this%moisture_response_function)
#else
  f_w = 1.0d0
#endif

#ifdef CLM_PFLOTRAN
  if(this%ph_response_function > 0) then
    if(this%fixed_ph > 0.0d0) then
       ph = this%fixed_ph
    else
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
      f_ph = GetpHResponse(ph, this%ph_response_function)  
  endif 
#else
  f_ph = 1.0d0
#endif 
  
  if(f_t < 1.0d-20 .or. f_w < 1.0d-20 .or. f_ph < 1.0d-20) then
     return
  endif

! residual
  Lwater = porosity * global_auxvar%sat(iphase)*volume*1.d3

  rate = this%rate_constant * Lwater * f_t * f_w * f_ph

!  write(*, *) this%rate_constant

  drate_donor    = rate
  drate_acceptor = rate
  drate_biomass  = rate

! biomass
  icomp = this%microbial_mass_ispec
  if(icomp > 0) then
     conc = rt_auxvar%immobile(icomp)
     rate = rate * conc 
     if (compute_derivative) then 
!        drate_biomass = drate_biomass * 1.0
         drate_donor = drate_donor * conc
         drate_acceptor = drate_acceptor * conc
     endif
  endif

! electron donor
  icomp = this%electron_donor_ispec
  if(icomp > 0) then
     activity = (rt_auxvar%pri_molal(icomp) - this%x0eps)* rt_auxvar%pri_act_coef(icomp)
     tmp_real = activity /(activity + this%electron_donor_half_saturation) 
     rate = rate * tmp_real
     if (compute_derivative) then 
        drate_biomass = drate_biomass * tmp_real
        drate_donor = drate_donor * this%electron_donor_half_saturation / &
                   (activity + this%electron_donor_half_saturation) / &
                   (activity + this%electron_donor_half_saturation) * &
                   rt_auxvar%pri_act_coef(icomp)
        drate_acceptor = drate_acceptor * tmp_real
     endif
  endif

! electron acceptor
  icomp = this%electron_acceptor_ispec
  if(icomp > 0) then
     activity = (rt_auxvar%pri_molal(icomp) - this%x0eps) * rt_auxvar%pri_act_coef(icomp)
     tmp_real = activity /(activity + this%electron_acceptor_half_saturation) 
     rate = rate * tmp_real
     if (compute_derivative) then 
        drate_biomass = drate_biomass * tmp_real
        drate_donor = drate_donor * tmp_real
        drate_acceptor = drate_acceptor * this%electron_acceptor_half_saturation / &
                   (activity + this%electron_acceptor_half_saturation) / &
                   (activity + this%electron_acceptor_half_saturation) * &
                   rt_auxvar%pri_act_coef(icomp)
     endif 
  endif

! inhibition term
  do i = 1, this%nInhibition
     icomp = this%ispec_inh(i)
     activity = rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp) 
     tmp_real = this%inhibition_coef(i)/(activity + this%inhibition_coef(i))  
     rate = rate * tmp_real  
     drate_biomass = drate_biomass * tmp_real
  enddo
  
  do i = 1, this%dbase_rxn%nspec
     if(this%dbase_rxn%spec_ids(i) == this%microbial_mass_ispec) then
       ires = this%dbase_rxn%spec_ids(i) + reaction%offset_immobile
     else
       ires = this%dbase_rxn%spec_ids(i)
     endif
!     Residual(this%dbase_rxn%spec_ids(i)) = Residual(this%dbase_rxn%spec_ids(i)) &
     Residual(this%dbase_rxn%spec_ids(i)) = Residual(ires) &
                                     - this%dbase_rxn%stoich(i) * rate 
  enddo

! Jacobian
  if (.not.compute_derivative) return 

! biomass
  icomp = this%microbial_mass_ispec
  if(icomp > 0) then
     do i = 1, this%dbase_rxn%nspec
        if(this%dbase_rxn%spec_ids(i) == this%microbial_mass_ispec) then
          ires = this%dbase_rxn%spec_ids(i) + reaction%offset_immobile
        else
          ires = this%dbase_rxn%spec_ids(i)
        endif
       Jacobian(ires, icomp + reaction%offset_immobile) = &
           Jacobian(ires, icomp + reaction%offset_immobile) &
           - this%dbase_rxn%stoich(i) * drate_biomass 
     enddo
  endif

! electron donor     
  icomp = this%electron_donor_ispec
  if(icomp > 0) then
     do i = 1, this%dbase_rxn%nspec
        if(this%dbase_rxn%spec_ids(i) == this%microbial_mass_ispec) then
          ires = this%dbase_rxn%spec_ids(i) + reaction%offset_immobile
        else
          ires = this%dbase_rxn%spec_ids(i)
        endif
        Jacobian(ires, icomp) = &
                Jacobian(ires,icomp) &
              - this%dbase_rxn%stoich(i) * drate_donor 
     enddo
  endif
! electron acceptor
  icomp = this%electron_acceptor_ispec
  if(icomp > 0) then
     do i = 1, this%dbase_rxn%nspec
        if(this%dbase_rxn%spec_ids(i) == this%microbial_mass_ispec) then
          ires = this%dbase_rxn%spec_ids(i) + reaction%offset_immobile
        else
          ires = this%dbase_rxn%spec_ids(i)
        endif
        Jacobian(ires, icomp) = Jacobian(ires,icomp) &
              - this%dbase_rxn%stoich(i) * drate_acceptor 
     enddo
  endif

!  write(*,*)(Residual(i), i = 1, reaction%ncomp)
!  do i = 1, reaction%ncomp
!     write(*,*)(Jacobian(icomp, i), icomp = 1, reaction%ncomp)
!  enddo

end subroutine MicrobialReact

! ************************************************************************** !
!
! MicrobialDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Guoping Tang
! date: 11/25/13
!
! ************************************************************************** !
subroutine MicrobialDestroy(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_Microbial_type) :: this  

  type(rate_type), pointer :: cur_rate, prev_rate

  call DatabaseRxnDestroy(this%dbase_rxn)    

  cur_rate => this%Inhibition
  do 
    if(.not.associated(cur_rate)) exit
       prev_rate => cur_rate
       cur_rate => cur_rate%next
       deallocate(prev_rate)
       nullify(prev_rate)
  enddo

  call DeallocateArray(this%ispec_inh) 

  call DeallocateArray(this%inhibition_coef) 
 
end subroutine MicrobialDestroy

end module Reaction_Sandbox_Microbial_class
