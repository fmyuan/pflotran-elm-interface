module PM_Well_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PM_Base_class
  use Option_module
  use Realization_Subsurface_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type :: well_grid_type
    PetscInt :: nsegments       ! number of well segments
    PetscReal :: dh             ! delta h discretization of each segment center
    PetscReal, pointer :: h(:)  ! h coordinate of each segment center
  end type well_grid_type

  type, public, extends(pm_base_type) :: pm_well_type
    class(realization_subsurface_type), pointer :: realization
    type(well_grid_type), pointer :: grid
  contains
    procedure, public :: Setup => PMWellSetup
    procedure, public :: ReadPMBlock => PMWellReadPMBlock
    procedure, public :: SetRealization => PMWellSetRealization
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

  nullify(PMWellCreate%grid)
  PMWellCreate%grid%nsegments = UNINITIALIZED_INTEGER
  PMWellCreate%grid%dh = UNINITIALIZED_DOUBLE
  nullify(PMWellCreate%grid%h)

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
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMWellSetRealization

! ************************************************************************** !

end module PM_Well_class