module PM_Fracture_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module
  use PM_Base_class
  use Option_module
  use String_module
  use Input_Aux_module
  use Geometry_module

  implicit none

  private

  type :: fracture_type 
    ! fracture ID number
    PetscInt :: id 
    ! cell permeability, initial [m]
    PetscReal, pointer :: kx0(:)
    ! cell permeability, initial [m]
    PetscReal, pointer :: ky0(:)
    ! cell permeability, initial [m]
    PetscReal, pointer :: kz0(:)
    ! fracture hydraulic aperture, initial [m]
    PetscReal :: hap0
    ! fracture hydraulic aperture, dynamic [m]
    PetscReal, pointer :: hap(:)
    ! fracture radius - x [m]
    PetscReal :: radx
    ! fracture radius - y [m]
    PetscReal :: rady
    ! fracture center coordinates [m]
    type(point3d_type) :: center
    ! fracture normal vector [m]
    type(point3d_type) :: normal
    ! fracture parallel vector [m]
    type(point3d_type) :: parallel
    ! rotation matrix
    PetscReal :: RM(3,3)
    ! list of global cell ids the fracture occupies
    PetscInt, pointer :: cell_ids(:)
    ! list of previous temperatures in the fracture cells
    PetscReal, pointer :: prev_temperature(:)
    ! list of change in temperatures in the fracture cells
    PetscReal, pointer :: dT(:)
    ! list of change in length in the fracture cells
    PetscReal, pointer :: dL(:)
    ! total number of cells this fracture occupies
    PetscInt :: ncells
    ! linked list next object
    class(fracture_type), pointer :: next 
  end type fracture_type

  type, public, extends(pm_base_type) :: pm_fracture_type
    class(realization_subsurface_type), pointer :: realization
    class(fracture_type), pointer :: fracture_list
    ! maximum distance grid cell center can be from the fracture plane in 
    ! order to mark the cell as being fractured # [m]
    PetscReal :: max_distance
    ! thermal expansion coefficient # [1/C]
    PetscReal :: t_coeff
    ! number of fractures in the domain
    PetscInt :: nfrac 
  contains
    procedure, public :: Setup => PMFracSetup
    procedure, public :: ReadPMBlock => PMFracReadPMBlock
    procedure, public :: SetRealization => PMFracSetRealization
    procedure, public :: InitializeRun => PMFracInitializeRun
    procedure, public :: FinalizeRun => PMFracFinalizeRun
    procedure, public :: InitializeTimestep => PMFracInitializeTimestep
    !procedure, public :: UpdateTimestep => PMFracUpdateTimestep
    procedure, public :: FinalizeTimestep => PMFracFinalizeTimestep
    !procedure, public :: PreSolve => PMFracPreSolve
    procedure, public :: Solve => PMFracSolve
    !procedure, public :: PostSolve => PMFracPostSolve
    procedure, public :: Destroy => PMFracDestroy
  end type pm_fracture_type

  public :: PMFracCreate, &
            PMFracReadPMBlock, &
            PMFracReadPass2

  contains

! ************************************************************************** !

function PMFracCreate()
  !
  ! Creates the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023

  implicit none

  class(pm_fracture_type), pointer :: PMFracCreate
  class(pm_fracture_type), pointer :: this

  allocate(this)
  call PMBaseInit(this)

  this%header = 'GEOTHERMAL FRACTURE MODEL'

  nullify(this%realization)
  nullify(this%fracture_list)
  this%nfrac = 0  
  this%max_distance = 5.0  
  this%t_coeff = 40.d-6

  PMFracCreate => this

end function PMFracCreate

! ************************************************************************** !

function PMFractureCreate()
  !
  ! Creates a fracture object.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023

  implicit none

  class(fracture_type), pointer :: PMFractureCreate

  allocate(PMFractureCreate)
  call PMFractureInit(PMFractureCreate)

end function PMFractureCreate

! ************************************************************************** !

subroutine PMFractureInit(this)
  !
  ! Initializes variables associated with a fracture object.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  implicit none

  class(fracture_type) :: this

  nullify(this%next)
  nullify(this%cell_ids)
  nullify(this%prev_temperature)
  nullify(this%dT)
  nullify(this%dL)
  nullify(this%hap)
  nullify(this%kx0); nullify(this%ky0); nullify(this%kz0)
  this%ncells = UNINITIALIZED_INTEGER
  this%id = UNINITIALIZED_INTEGER
  this%hap0 = UNINITIALIZED_DOUBLE
  this%radx = UNINITIALIZED_DOUBLE
  this%rady = UNINITIALIZED_DOUBLE
  this%center%x = UNINITIALIZED_DOUBLE
  this%center%y = UNINITIALIZED_DOUBLE
  this%center%z = UNINITIALIZED_DOUBLE
  this%normal%x = UNINITIALIZED_DOUBLE
  this%normal%y = UNINITIALIZED_DOUBLE
  this%normal%z = UNINITIALIZED_DOUBLE
  this%parallel%x = UNINITIALIZED_DOUBLE
  this%parallel%y = UNINITIALIZED_DOUBLE
  this%parallel%z = UNINITIALIZED_DOUBLE

end subroutine PMFractureInit

! ************************************************************************** !

subroutine PMFracSetup(this)
  !
  ! Initializes variables associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  use Grid_module

  implicit none

  class(pm_fracture_type) :: this

  type(option_type), pointer :: option
  type(fracture_type), pointer :: cur_fracture
  type(point3d_type) :: coordinate
  type(grid_type), pointer :: res_grid
  character(len=MAXWORDLENGTH) :: word
  PetscInt, pointer :: temp_cell_ids(:)
  PetscInt :: k,j,local_id,local_id_center
  PetscInt :: nf,read_max
  PetscReal :: D,distance 
  PetscReal :: a1,a2,a3,b1,b2,b3,c1,c2,c3,vmag
  PetscReal :: ra, sinra, cosra

  option => this%option
  res_grid => this%realization%patch%grid
  nf = 0
  read_max = 5000 ! Set to a large but realistic number

  allocate(temp_cell_ids(read_max))

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit
    nf = nf + 1
    cur_fracture => cur_fracture%next
  enddo
  write(word,'(i4)') nf

  option%io_buffer = ' '
  call PrintMsg(option)
  option%io_buffer = 'GEOTHERMAL_FRACTURE_MODEL: Mapping [' // trim(word) // &
    '] fractures onto equivalent porous media domain:'
  call PrintMsg(option)

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    write(word,'(i4)') cur_fracture%id 
    option%io_buffer = 'GEOTHERMAL_FRACTURE_MODEL: Mapping fracture ID# [' &
                       // trim(word) // '].'
    call PrintMsg(option)
    call GridGetLocalIDFromCoordinate(res_grid,cur_fracture%center,option, &
    	                              local_id_center)
    if (Uninitialized(local_id_center)) then
      option%io_buffer = 'Fracture ID# [' // trim(word) // '] CENTER &
        &coordinate is not within the reservoir domain.'
      call PrintErrMsg(option)
    endif
    ! Equation of a plane normal to vector (A,B,C) is
    ! Ax + By + Cz + D = 0
    ! Get D by knowing the plane must contain the center point (x,y,z)
    D = -1.d0*(cur_fracture%normal%x*cur_fracture%center%x + &
    	       cur_fracture%normal%y*cur_fracture%center%y + &
    	       cur_fracture%normal%z*cur_fracture%center%z)
    j = 0
    do k = 1,res_grid%nmax
      distance = cur_fracture%normal%x*res_grid%x(k) + &
                 cur_fracture%normal%y*res_grid%y(k) + &
                 cur_fracture%normal%z*res_grid%z(k) + D
      if (abs(distance) < this%max_distance) then
        j = j + 1
        temp_cell_ids(j) = k
      endif
    enddo
    cur_fracture%ncells = j
    allocate(cur_fracture%cell_ids(cur_fracture%ncells))
    allocate(cur_fracture%prev_temperature(cur_fracture%ncells))
    allocate(cur_fracture%dT(cur_fracture%ncells))
    allocate(cur_fracture%dL(cur_fracture%ncells))
    allocate(cur_fracture%hap(cur_fracture%ncells))
    allocate(cur_fracture%kx0(cur_fracture%ncells))
    allocate(cur_fracture%ky0(cur_fracture%ncells))
    allocate(cur_fracture%kz0(cur_fracture%ncells))

    cur_fracture%hap(:) = cur_fracture%hap0
    cur_fracture%cell_ids(1:cur_fracture%ncells) = &
                                          temp_cell_ids(1:cur_fracture%ncells)

    ! get the unit vector parallel to the fracture plane
    ! a X b = c, where c is perpendicular to a and b (a and b can't be parallel)
    ! take a as the unit normal, and b as slightly rotated unit normal
    ! then c will be perpendicular to a, and parallel to fracture plane
    a1 = cur_fracture%normal%x; b1 = cur_fracture%normal%x
    a2 = cur_fracture%normal%y; b2 = cur_fracture%normal%y 
    a3 = cur_fracture%normal%z; b3 = cur_fracture%normal%z + 1
    c1 = a2*b3 - a3*b2
    c2 = -1.d0*(a1*b3 - a3*b1)
    c3 = a1*b2 - a2*b1
    cur_fracture%parallel%x = c1
    cur_fracture%parallel%y = c2
    cur_fracture%parallel%z = c3
    vmag = sqrt((cur_fracture%parallel%x)**2.d0 + &
        	    (cur_fracture%parallel%y)**2.d0 + &
        	    (cur_fracture%parallel%z)**2.d0)
    cur_fracture%parallel%x = cur_fracture%parallel%x/vmag
    cur_fracture%parallel%y = cur_fracture%parallel%y/vmag
    cur_fracture%parallel%z = cur_fracture%parallel%z/vmag
    c1 = cur_fracture%parallel%x
    c2 = cur_fracture%parallel%y
    c3 = cur_fracture%parallel%z
    
    ! calculate the rotation angle (ra) between domain/fracture orientations
    ra = acos(a3)

    ! calculate the rotation matrix between domain/fracture orientations
    ! this is done knowing the rotation angle, and fracture parallel direction
    sinra = sin(ra); cosra = cos(ra)
    cur_fracture%RM(1,1) = cosra + c1*c1*(1.d0-cosra)
    cur_fracture%RM(1,2) = c1*c2*(1.d0-cosra) - c3*sinra
    cur_fracture%RM(1,3) = c1*c3*(1.d0-cosra) + c2*sinra
    cur_fracture%RM(2,1) = c2*c1*(1.d0-cosra) + c3*sinra
    cur_fracture%RM(2,2) = cosra + c2*c2*(1.d0-cosra)
    cur_fracture%RM(2,3) = c2*c3*(1.d0-cosra) -c1*sinra
    cur_fracture%RM(3,1) = c3*c1*(1.d0-cosra) -c2*sinra
    cur_fracture%RM(3,2) = c3*c2*(1.d0-cosra) + c1*sinra
    cur_fracture%RM(3,3) = cosra + c3*c3*(1.d0-cosra)

    cur_fracture => cur_fracture%next
  enddo

  deallocate(temp_cell_ids)

end subroutine PMFracSetup

! ************************************************************************** !

subroutine PMFracSetRealization(this,realization)
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023

  implicit none

  class(pm_fracture_type) :: this
  class(realization_subsurface_type), pointer :: realization

  this%realization => realization
  this%realization_base => realization

end subroutine PMFracSetRealization

! ************************************************************************** !

subroutine PMFracInitializeRun(this)
  !
  ! Initializes the run associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  use Material_Aux_module
  use Global_Aux_module

  implicit none

  class(pm_fracture_type) :: this

  type(material_auxvar_type), pointer :: material_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(option_type), pointer :: option
  type(fracture_type), pointer :: cur_fracture
  PetscInt :: icell,k
  PetscReal :: kx,ky,kz
  PetscReal :: L 
  PetscErrorCode :: ierr

  option => this%option

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    do k = 1,cur_fracture%ncells
      icell = cur_fracture%cell_ids(k)
      material_auxvar => &
        this%realization%patch%aux%material%auxvars(icell) ! the input here
                                                           ! should be ghosted_id
      ! if cells aren't too pancaked, L should be a good estimate of length
      L = (material_auxvar%volume)**(1.d0/3.d0) 
      cur_fracture%kx0(k) = material_auxvar%permeability(1)
      cur_fracture%ky0(k) = material_auxvar%permeability(2)
      cur_fracture%kz0(k) = material_auxvar%permeability(3)

      call PMFracCalcK(L,cur_fracture%hap0,cur_fracture%RM,kx,ky,kz)

      call PMFracUpdateK(material_auxvar,cur_fracture%kx0(k), &
      	                 cur_fracture%ky0(k),cur_fracture%kz0(k), &
      	                 kx,ky,kz)

      global_auxvar => this%realization%patch%aux%Global%auxvars(icell)! the input here
                                                                       ! should be ghosted_id
      cur_fracture%prev_temperature(k) = global_auxvar%temp
    enddo

    cur_fracture => cur_fracture%next
  enddo

end subroutine PMFracInitializeRun

! ************************************************************************** !

subroutine PMFracFinalizeRun(this)
  !
  ! Finalizes the run associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  implicit none

  class(pm_fracture_type) :: this

end subroutine PMFracFinalizeRun

! ************************************************************************** !

subroutine PMFracInitializeTimestep(this)
  !
  ! Initializes the timestep associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  implicit none

  class(pm_fracture_type) :: this


end subroutine PMFracInitializeTimestep

! ************************************************************************** !

subroutine PMFracFinalizeTimestep(this)
  !
  ! Finalizes the timestep associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  implicit none

  class(pm_fracture_type) :: this

end subroutine PMFracFinalizeTimestep

! ************************************************************************** !

subroutine PMFracReadPMBlock(this,input)
  !
  ! Reads input file parameters associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023

  implicit none

  class(pm_fracture_type) :: this
  type(input_type), pointer :: input

  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option
  input%ierr = 0
  error_string = 'GEOTHERMAL_FRACTURE_MODEL'

  option%io_buffer = 'pflotran card:: GEOTHERMAL_FRACTURE_MODEL'
  call PrintMsg(option)

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    found = PETSC_FALSE

    ! Read keywords within GEOTHERMAL_FRACTURE_MODEL block:
    select case(trim(word))
      case('MAX_DISTANCE')
        call InputReadDouble(input,option,this%max_distance)
        call InputErrorMsg(input,option,'MAX_DISTANCE',error_string)
        cycle
    !-------------------------------------
      case('THERMAL_EXPANSION_COEFFICIENT')
        call InputReadDouble(input,option,this%t_coeff)
        call InputErrorMsg(input,option,'THERMAL_EXPANSION_COEFFICIENT', &
        	               error_string)
        cycle
    !-------------------------------------
    end select

    ! Read sub-blocks within GEOTHERMAL_FRACTURE_MODEL block:
    error_string = 'GEOTHERMAL_FRACTURE_MODEL'
    call PMFracReadFracture(this,input,option,word,error_string,found)
    if (found) cycle

    if (.not. found) then
      option%io_buffer = 'Keyword "' // trim(word) // &
                         '" does not exist for GEOTHERMAL_FRACTURE_MODEL.'
      call PrintErrMsg(option)
    endif

  enddo
  call InputPopBlock(input,option)

end subroutine PMFracReadPMBlock

! ************************************************************************** !

subroutine PMFracReadFracture(pm_fracture,input,option,keyword,error_string, &
	                          found)
  !
  ! Reads input file parameters associated with the fracture model grid.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  implicit none

  class(pm_fracture_type) :: pm_fracture
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  character(len=MAXWORDLENGTH) :: word
  class(fracture_type), pointer :: new_fracture,cur_fracture
  PetscReal :: vmag
  PetscInt :: num_errors
  PetscBool :: added

  error_string = trim(error_string) // ',FRACTURE'
  found = PETSC_TRUE
  added = PETSC_FALSE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('FRACTURE')
      allocate(new_fracture)
      new_fracture => PMFractureCreate()
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !-----------------------------
          case('ID')
            call InputReadInt(input,option,new_fracture%id) 
            call InputErrorMsg(input,option,'ID',error_string)
        !-----------------------------
          case('HYDRAULIC_APERTURE')
            call InputReadDouble(input,option,new_fracture%hap0)
            call InputErrorMsg(input,option,'HYDRAULIC_APERTURE',error_string)
        !-----------------------------
          case('RADIUS_X')
            call InputReadDouble(input,option,new_fracture%radx)
            call InputErrorMsg(input,option,'RADIUS_X',error_string)
        !-----------------------------
          case('RADIUS_Y')
            call InputReadDouble(input,option,new_fracture%rady)
            call InputErrorMsg(input,option,'RADIUS_Y',error_string)
        !-----------------------------
          case('CENTER')
            call InputReadDouble(input,option,new_fracture%center%x)
            call InputErrorMsg(input,option,'CENTER, X COORDINATE', &
            	               error_string)
            call InputReadDouble(input,option,new_fracture%center%y)
            call InputErrorMsg(input,option,'CENTER, Y COORDINATE', &
            	               error_string)
            call InputReadDouble(input,option,new_fracture%center%z)
            call InputErrorMsg(input,option,'CENTER, Z COORDINATE', &
            	               error_string)
        !-----------------------------
          case('NORMAL_VECTOR')
            call InputReadDouble(input,option,new_fracture%normal%x)
            call InputErrorMsg(input,option,'NORMAL_VECTOR, X COORDINATE', &
            	               error_string)
            call InputReadDouble(input,option,new_fracture%normal%y)
            call InputErrorMsg(input,option,'NORMAL_VECTOR, Y COORDINATE', &
            	               error_string)
            call InputReadDouble(input,option,new_fracture%normal%z)
            call InputErrorMsg(input,option,'NORMAL_VECTOR, Z COORDINATE', &
            	               error_string)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)
      !------ error messaging ----------------------------------------
      if (Uninitialized(new_fracture%id)) then
        option%io_buffer = 'ERROR: ID must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%hap0)) then
        option%io_buffer = 'ERROR: HYDRAULIC_APERTURE must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%radx)) then
        option%io_buffer = 'ERROR: RADIUS_X must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%rady)) then
        option%io_buffer = 'ERROR: RADIUS_Y must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%center%x) .or. &
      	  Uninitialized(new_fracture%center%y) .or. &
      	  Uninitialized(new_fracture%center%z)) then
        option%io_buffer = 'ERROR: CENTER coordinate must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracture%normal%x) .or. &
      	  Uninitialized(new_fracture%normal%y) .or. &
      	  Uninitialized(new_fracture%normal%z)) then
        option%io_buffer = 'ERROR: NORMAL_VECTOR coordinate must be specified &
                            &in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      else
        vmag = sqrt((new_fracture%normal%x)**2.d0 + &
        	        (new_fracture%normal%y)**2.d0 + &
        	        (new_fracture%normal%z)**2.d0)
        new_fracture%normal%x = new_fracture%normal%x/vmag
        new_fracture%normal%y = new_fracture%normal%y/vmag
        new_fracture%normal%z = new_fracture%normal%z/vmag
      endif
      if (num_errors > 0) then
        write(option%io_buffer,*) num_errors
        option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                           &the FRACTURE block. See above to fix.'
        call PrintErrMsg(option)
      endif
      !------ add fracture to list -----------------------------------
      if (.not.associated(pm_fracture%fracture_list)) then
        pm_fracture%fracture_list => new_fracture
      else
        cur_fracture => pm_fracture%fracture_list
        do
          if (.not.associated(cur_fracture)) exit
          if (.not.associated(cur_fracture%next)) then
            cur_fracture%next => new_fracture
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_fracture => cur_fracture%next
        enddo
      endif
      nullify(new_fracture)
  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

end subroutine PMFracReadFracture

! ************************************************************************** !

subroutine PMFracReadPass2(input,option)
  !
  ! Reads input file parameters associated with the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: card

  error_string = 'SUBSURFACE,GEOTHERMAL_FRACTURE_MODEL'

  input%ierr = 0
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,card)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      !--------------------
      case('FRACTURE')
        call InputSkipToEND(input,option,card)
      !--------------------
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine PMFracReadPass2

! ************************************************************************** !

subroutine PMFracCalcK(L,hap,RM,kx,ky,kz)
  !
  ! Calculates cell anisotropic permeability from fracture permeability.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/15/2023
  !

  implicit none

  PetscReal :: L 
  PetscReal :: hap
  PetscReal :: RM(3,3)
  PetscReal :: kx,ky,kz

  PetscReal :: K_frac(3,3)
  PetscReal :: K_frac_rot(3,3)
  PetscReal :: K_domain(3,3)

  K_frac(:,:) = 0.d0
  K_frac(1,1) = (1.d0/12.d0)*(hap**3.0d0)/L ! kxx
  K_frac(2,2) = (1.d0/12.d0)*(hap**3.0d0)/L ! kyy
  ! note: zero in kzz

  K_frac_rot = MATMUL(RM,K_frac)
  K_domain = MATMUL(K_frac_rot,TRANSPOSE(RM))

  kx = K_domain(1,1)
  ky = K_domain(2,2)
  kz = K_domain(3,3)

end subroutine PMFracCalcK

! ************************************************************************** !

subroutine PMFracUpdateK(material_auxvar,kx0,ky0,kz0,kx,ky,kz)
  !
  ! Calculates cell permeability from fracture permeability.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/15/2023
  !
  use Material_Aux_module

  implicit none

  type(material_auxvar_type), pointer :: material_auxvar
  PetscReal :: kx0,ky0,kz0
  PetscReal :: kx,ky,kz

  material_auxvar%permeability(1) = kx0+kx ! to fix: this will keep adding background over and over again
  material_auxvar%permeability(2) = ky0+ky ! but maybe it doesn't matter b/c background is really tiny
  material_auxvar%permeability(3) = kz0+kz

end subroutine PMFracUpdateK

! ************************************************************************** !

subroutine PMFracSolve(this,time,ierr)
  !
  ! Main solve step for the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  use Material_Aux_module
  use Global_Aux_module

  implicit none

  class(pm_fracture_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(material_auxvar_type), pointer :: material_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(option_type), pointer :: option
  type(fracture_type), pointer :: cur_fracture
  PetscReal :: cur_temperature,L
  PetscReal :: kx,ky,kz
  PetscInt :: icell,k

  option => this%option

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    do k = 1,cur_fracture%ncells
      icell = cur_fracture%cell_ids(k)
      material_auxvar => &
        this%realization%patch%aux%material%auxvars(icell) ! the input here
                                                           ! should be ghosted_id
      global_auxvar => this%realization%patch%aux%Global%auxvars(icell)! the input here
                                                                       ! should be ghosted_id
      cur_temperature = global_auxvar%temp
      cur_fracture%dT(k) = cur_temperature - cur_fracture%prev_temperature(k)

      ! if cells aren't too pancaked, L should be a good estimate of length
      L = (material_auxvar%volume)**(1.d0/3.d0)
      cur_fracture%dL(k) = L*cur_fracture%dT(k)*this%t_coeff

      ! update the aperture with dL and calculate new perm
      cur_fracture%hap(k) = cur_fracture%hap(k) - cur_fracture%dL(k)

      ! get anisotropic domain permeability from rotation transformation
      call PMFracCalcK(L,cur_fracture%hap(k),cur_fracture%RM,kx,ky,kz)

      ! update domain material permeability of cells with fractures
      call PMFracUpdateK(material_auxvar,cur_fracture%kx0(k), &
      	                 cur_fracture%ky0(k),cur_fracture%kz0(k), &
      	                 kx,ky,kz)

      ! reset previous temperature for next time step
      cur_fracture%prev_temperature(k) = cur_temperature
    enddo

    cur_fracture => cur_fracture%next
  enddo

  ierr = 0 ! If this is not set to zero, TS_STOP_FAILURE occurs

end subroutine PMFracSolve

! ************************************************************************** !

subroutine PMFracDestroy(this)
  !
  ! Destroys objects in the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !

  implicit none

  class(pm_fracture_type) :: this

end subroutine PMFracDestroy

! ************************************************************************** !

end module PM_Fracture_class