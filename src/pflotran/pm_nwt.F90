module PM_NWT_class
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class 
  use Realization_Subsurface_class
  use Communicator_Base_module  
  use Option_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
  type, public, extends(pm_base_type) :: pm_nwt_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    ! local variables
    PetscBool :: check_post_convergence
    ! governs the size of subsequent time steps
    PetscReal, pointer :: max_concentration_change(:)
    PetscReal, pointer :: max_volfrac_change(:)
    PetscReal :: volfrac_change_governor
    PetscReal :: cfl_governor
    PetscBool :: temperature_dependent_diffusion
  contains
    procedure, public :: Setup => PMNWTSetup 
    procedure, public :: Read => PMNWTRead
    procedure, public :: SetRealization => PMNWTSetRealization   
  end type pm_nwt_type
  
  public :: PMNWTCreate
  
  
contains

! ************************************************************************** !

function PMNWTCreate()
  ! 
  ! Creates the nuclear waste transport process model shell.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  ! 

  implicit none
  
  class(pm_nwt_type), pointer :: PMNWTCreate

  class(pm_nwt_type), pointer :: nwt_pm
    
  allocate(nwt_pm)
  nullify(nwt_pm%option)
  nullify(nwt_pm%output_option)
  nullify(nwt_pm%realization)
  nullify(nwt_pm%comm1)
  
  ! local variables
  nwt_pm%check_post_convergence = PETSC_FALSE
  nullify(nwt_pm%max_concentration_change)
  nullify(nwt_pm%max_volfrac_change)
  nwt_pm%volfrac_change_governor = 1.d0
  nwt_pm%cfl_governor = UNINITIALIZED_DOUBLE
  nwt_pm%temperature_dependent_diffusion = PETSC_FALSE

  call PMBaseInit(nwt_pm)
  nwt_pm%name = 'Nuclear Waste Transport'
  nwt_pm%header = 'NUCLEAR WASTE TRANSPORT'
  
  PMNWTCreate => nwt_pm
  
end function PMNWTCreate

! ************************************************************************** !
  
subroutine PMNWTSetup(this)
  ! 
  ! Initializes variables associated with nuclear waste transport.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  ! 
  
  use NW_Transport_Aux_module

  implicit none
  
  class(pm_nwt_type) :: this
  
  type(nw_transport_param_type), pointer :: nwt_parameter
  
  ! jenn:todo get NWT in patch%aux pointed to correct place
  nwt_parameter => this%realization%patch%aux%NWT%nwt_parameter
  
  ! pass down flags from PMNWT class
  nwt_parameter%temperature_dependent_diffusion = &
    this%temperature_dependent_diffusion
  
  ! set the communicator
  this%comm1 => this%realization%comm1
  
  allocate(this%max_concentration_change(nwt_parameter%ncomp))
  allocate(this%max_volfrac_change(nwt_parameter%ncomp))

end subroutine PMNWTSetup

! ************************************************************************** !

subroutine PMNWTRead(this,input)
  ! 
  ! Reads input file parameters associated with the nuclear waste transport 
  ! process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use NW_Transport_Aux_module
 
  implicit none
  
  class(pm_nwt_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option
  
  error_string = 'NUCLEAR_WASTE_TRANSPORT OPTIONS'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    found = PETSC_FALSE
    call PMBaseReadSelectCase(this,input,keyword,found,error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('GLOBAL_IMPLICIT','OPERATOR_SPLIT','OPERATOR_SPLITTING')
      case('MAX_VOLUME_FRACTION_CHANGE')
        call InputReadDouble(input,option,this%volfrac_change_governor)
        call InputErrorMsg(input,option,'MAX_VOLUME_FRACTION_CHANGE', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
      case('ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,nwt_itol_rel_update)
        call InputErrorMsg(input,option,'ITOL_RELATIVE_UPDATE', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
        this%check_post_convergence = PETSC_TRUE
      case('MINIMUM_SATURATION')
        call InputReadDouble(input,option,nwt_min_saturation)
        call InputErrorMsg(input,option,'MINIMUM_SATURATION', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
      case('NUMERICAL_JACOBIAN')
        option%transport%numerical_derivatives = PETSC_TRUE
      case('TEMPERATURE_DEPENDENT_DIFFUSION')
        this%temperature_dependent_diffusion = PETSC_TRUE
      case('MAX_CFL')
        call InputReadDouble(input,option,this%cfl_governor)
        call InputErrorMsg(input,option,'MAX_CFL', &
                           'NUCLEAR_WASTE_TRANSPORT OPTIONS')
      case('MULTIPLE_CONTINUUM')
        option%use_mc = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  enddo
  
end subroutine PMNWTRead

! ************************************************************************** !

subroutine PMNWTSetRealization(this,realization)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/08/2019
  ! 

  use Realization_Subsurface_class  

  implicit none
  
  class(pm_nwt_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization
  
  if (realization%reaction%use_log_formulation) then
    this%solution_vec = realization%field%tran_log_xx
  else
    this%solution_vec = realization%field%tran_xx
  endif
  this%residual_vec = realization%field%tran_r
  
end subroutine PMNWTSetRealization

! ************************************************************************** !
  

end module PM_NWT_class