module PM_Well_class

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscsnes.h"
  use petscsys
  use petscsnes
  use PM_Base_class
  use Option_module
  use Realization_Subsurface_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type :: well_grid_type
    ! number of well segments
    PetscInt :: nsegments      
    ! delta h discretization of each segment center       
    PetscReal, pointer :: dh(:)      
    ! h coordinate of each segment center 
    PetscReal, pointer :: h(:)        
    ! gravity vector magnitude
    PetscReal :: g
  end type well_grid_type

  type :: well_type
    ! cross-sectional area of each well segment
    PetscReal, pointer :: area(:)         
    ! diameter of each well segment      
    PetscReal, pointer :: diameter(:) 
    ! volume of each well segment
    PetscReal, pointer :: volume(:) 
    ! friction ceofficient of each well segment        
    PetscReal, pointer :: f(:)      
    ! well index of each well segment [0,1]  0 = cased; 1 = open     
    PetscReal, pointer :: WI(:)    
  end type well_type

  type, public, extends(pm_base_type) :: pm_well_type
    class(realization_subsurface_type), pointer :: realization
    type(well_grid_type), pointer :: grid
    type(well_type), pointer :: well
  contains
    procedure, public :: Setup => PMWellSetup
    procedure, public :: ReadPMBlock => PMWellReadPMBlock
    procedure, public :: SetRealization => PMWellSetRealization
    procedure, public :: InitializeRun => PMWellInitializeRun
    procedure, public :: FinalizeRun => PMWellFinalizeRun
    procedure, public :: InitializeTimestep => PMWellInitializeTimestep
    procedure, public :: UpdateTimestep => PMWellUpdateTimestep
    procedure, public :: FinalizeTimestep => PMWellFinalizeTimestep
    procedure, public :: Residual => PMWellResidual
    procedure, public :: Jacobian => PMWellJacobian
    procedure, public :: PreSolve => PMWellPreSolve
    procedure, public :: PostSolve => PMWellPostSolve
    procedure, public :: Destroy => PMWellDestroy
  end type pm_well_type

  public :: PMWellCreate

  contains

! ************************************************************************** !

function PMWellCreate()
  ! 
  ! Creates the well process model.
  ! 
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type), pointer :: PMWellCreate
  
  allocate(PMWellCreate)
  call PMBaseInit(PMWellCreate)

  nullify(PMWellCreate%realization)

  ! create the grid object:
  nullify(PMWellCreate%grid)
  PMWellCreate%grid%nsegments = UNINITIALIZED_INTEGER
  nullify(PMWellCreate%grid%dh)
  nullify(PMWellCreate%grid%h)
  PMWellCreate%grid%g = 9.81

  ! create the well object:
  nullify(PMWellCreate%well)
  nullify(PMWellCreate%well%area)
  nullify(PMWellCreate%well%diameter)
  nullify(PMWellCreate%well%volume)
  nullify(PMWellCreate%well%f)
  nullify(PMWellCreate%well%WI)


end function PMWellCreate

! ************************************************************************** !
  
subroutine PMWellSetup(this)
  ! 
  ! Initializes variables associated with the well process model.
  ! 
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/04/2021
  ! 

  implicit none
  
  class(pm_well_type) :: this
  
  ! placeholder

end subroutine PMWellSetup

! ************************************************************************** !

subroutine PMWellReadPMBlock(this,input)
  ! 
  ! Reads input file parameters associated with the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  use Input_Aux_module
  use String_module
  use Utility_module
  
  implicit none

  class(pm_well_type) :: this
  type(input_type), pointer :: input

  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  option => this%option
  
  option%io_buffer = 'pflotran card:: WELL_MODEL'
  call PrintMsg(option)

end subroutine PMWellReadPMBlock

! ************************************************************************** !

subroutine PMWellSetRealization(this,realization)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  use Realization_Subsurface_class

  implicit none
  
  class(pm_well_type) :: this
  class(realization_subsurface_type), pointer :: realization

  ! There will probably need to be a realization_well created because we need
  ! our own grid and solution vectors and all that. But we will still need to
  ! be connected to the subsurface realization too.
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMWellSetRealization

! ************************************************************************** !

recursive subroutine PMWellInitializeRun(this)
  ! 
  ! Initializes the simulation run for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  
  implicit none

  class(pm_well_type) :: this
  
  ! placeholder
  
end subroutine PMWellInitializeRun

! ************************************************************************** !

recursive subroutine PMWellFinalizeRun(this)
  ! 
  ! Finalizes the simulation run for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  
  ! placeholder 
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMWellFinalizeRun

! ************************************************************************** !

subroutine PMWellInitializeTimestep(this)
  !
  ! Initializes the time step for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  
  implicit none
  
  class(pm_well_type) :: this

  ! placeholder

end subroutine PMWellInitializeTimestep

! ************************************************************************** !

subroutine PMWellUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                               num_newton_iterations,tfac, &
                               time_step_max_growth_factor)
  !
  ! Updates the time step for the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  
  implicit none
  
  class(pm_well_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  ! placeholder

end subroutine PMWellUpdateTimestep

! ************************************************************************** !

subroutine PMWellFinalizeTimestep(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this

  ! placeholder
  
end subroutine PMWellFinalizeTimestep

! ************************************************************************** !

subroutine PMWellResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
    
  ! placeholder

end subroutine PMWellResidual

! ************************************************************************** !

subroutine PMWellJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  ! placeholder

end subroutine PMWellJacobian

! ************************************************************************** !

subroutine PMWellPreSolve(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  
  ! placeholder
  
end subroutine PMWellPreSolve

! ************************************************************************** !

subroutine PMWellPostSolve(this)
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021

  implicit none
  
  class(pm_well_type) :: this
  
  ! placeholder
  
end subroutine PMWellPostSolve

! ************************************************************************** !

subroutine PMWellDestroy(this)
  ! 
  ! Destroys the well process model.
  ! 
  ! Author: Jennifer M. Frederick
  ! Date: 08/04/2021
  !
  use Utility_module, only : DeallocateArray
  use Option_module

  implicit none
  
  class(pm_well_type) :: this
    
  call PMBaseDestroy(this)

  call DeallocateArray(this%grid%h)
  
end subroutine PMWellDestroy

! ************************************************************************** !

end module PM_Well_class