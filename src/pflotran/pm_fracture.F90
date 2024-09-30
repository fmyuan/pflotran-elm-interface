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
    ! fracture radius - z [m]
    PetscReal :: radz
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
    ! maximum distance grid cell center can be from the fracture plane in
    ! order to mark the cell as being fractured # [m]
    PetscReal :: max_distance
    ! linked list next object
    type(fracture_type), pointer :: next
  end type fracture_type

  type :: fracfam_type
    ! fracture family ID number
    PetscInt :: id
    ! number of fractures in the family
    PetscInt :: nfrac_in_fam
    ! fracture intensity of fracture family [num fractures]
    PetscInt :: intensity_fracfam
    ! total surface area of fracture family [m^2] (includes both sides of fracture)
    PetscReal :: surface_area_fracfam
    ! fracture family hydraulic aperture, initial [m]
    PetscReal :: hap0
    ! fracture family hydraulic aperture, initial maximum [m]
    PetscReal :: hap0_max
    ! fracture family hydraulic standard deviation, initial [m]
    PetscReal :: hap0_stdev
    ! fracture family hydraulic aperture seed [-]
    PetscInt :: hap0_seed
    ! fracture family center coordinates [m]
    type(point3d_type) :: center
    ! fracture family center coordinates standard deviation [m]
    type(point3d_type) :: center_stdev
    ! fracture family center seed [-]
    PetscInt :: center_seed
    ! fracture family normal vector [m]
    type(point3d_type) :: normal
    ! fracture family normal vector coordinates standard deviation [m]
    type(point3d_type) :: normal_stdev
    ! fracture family normal vector seed [-]
    PetscInt :: normal_seed
    ! fracture family radius lengths [m]
    type(point3d_type) :: radius
    ! fracture family radius lengths standard deviation [m]
    type(point3d_type) :: radius_stdev
    ! fracture family radius seed [-]
    PetscInt :: radius_seed
    ! maximum distance grid cell center can be from the fracture plane in
    ! order to mark the cell as being fractured # [m]
    PetscReal :: max_distance
    ! list of fractures in this family
    type(fracture_type), pointer :: fracture_list
    ! linked list next object
    type(fracfam_type), pointer :: next
  end type fracfam_type

  type, public, extends(pm_base_type) :: pm_fracture_type
    class(realization_subsurface_type), pointer :: realization
    type(fracture_type), pointer :: fracture_list
    type(fracfam_type), pointer :: fracfam_list
    ! list of all global cell ids that contain fractures (some repeated)
    PetscInt, pointer :: allfrac_cell_ids(:)
    ! list of all global cell ids that contain fracture intersections
    PetscInt, pointer :: frac_common_cell_ids(:)
    ! thermal expansion coefficient # [1/C]
    PetscReal :: t_coeff
    ! number of fractures in the domain
    PetscInt :: nfrac
    ! maximum number of fractures in the domain
    PetscInt :: max_frac
    ! flag to update material_auxvar (permeability) object
    PetscBool :: update_material_auxvar_perm
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

  public :: PMFractureCreate, &
            FractureCreate, &
            FractureFamCreate, &
            PMFracReadPMBlock, &
            PMFracReadPass2

  contains

! ************************************************************************** !

function PMFractureCreate()
  !
  ! Creates the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023

  implicit none

  class(pm_fracture_type), pointer :: PMFractureCreate
  class(pm_fracture_type), pointer :: this

  allocate(this)
  call PMBaseInit(this)

  this%header = 'GEOTHERMAL FRACTURE MODEL'

  nullify(this%realization)
  nullify(this%fracture_list)
  nullify(this%fracfam_list)
  nullify(this%allfrac_cell_ids)
  nullify(this%frac_common_cell_ids)
  this%nfrac = 0
  this%max_frac = UNINITIALIZED_INTEGER
  this%t_coeff = UNINITIALIZED_DOUBLE ! 40.d-6 is value for granite
  this%update_material_auxvar_perm = PETSC_FALSE

  PMFractureCreate => this

end function PMFractureCreate

! ************************************************************************** !

function FractureCreate()
  !
  ! Creates a fracture object.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023

  implicit none

  type(fracture_type), pointer :: FractureCreate

  allocate(FractureCreate)
  call FractureInit(FractureCreate)

end function FractureCreate

! ************************************************************************** !

function FractureFamCreate()
  !
  ! Creates a fracture family object.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/19/2024

  implicit none

  type(fracfam_type), pointer :: FractureFamCreate

  allocate(FractureFamCreate)
  call FractureFamInit(FractureFamCreate)

end function FractureFamCreate

! ************************************************************************** !

subroutine FractureInit(this)
  !
  ! Initializes variables associated with a fracture object.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  implicit none

  type(fracture_type) :: this

  nullify(this%next)
  nullify(this%cell_ids)
  nullify(this%prev_temperature)
  nullify(this%dT)
  nullify(this%dL)
  nullify(this%hap)
  nullify(this%kx0); nullify(this%ky0); nullify(this%kz0)
  this%max_distance = 5.0
  this%ncells = UNINITIALIZED_INTEGER
  this%id = UNINITIALIZED_INTEGER
  this%hap0 = UNINITIALIZED_DOUBLE
  this%radx = UNINITIALIZED_DOUBLE
  this%rady = UNINITIALIZED_DOUBLE
  this%radz = UNINITIALIZED_DOUBLE
  this%center%x = UNINITIALIZED_DOUBLE
  this%center%y = UNINITIALIZED_DOUBLE
  this%center%z = UNINITIALIZED_DOUBLE
  this%normal%x = UNINITIALIZED_DOUBLE
  this%normal%y = UNINITIALIZED_DOUBLE
  this%normal%z = UNINITIALIZED_DOUBLE
  this%parallel%x = UNINITIALIZED_DOUBLE
  this%parallel%y = UNINITIALIZED_DOUBLE
  this%parallel%z = UNINITIALIZED_DOUBLE

end subroutine FractureInit

! ************************************************************************** !

subroutine FractureFamInit(this)
  !
  ! Initializes variables associated with a fracture family object.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/19/2024
  !
  implicit none

  type(fracfam_type) :: this

  nullify(this%next)
  nullify(this%fracture_list)
  this%id = UNINITIALIZED_INTEGER
  this%nfrac_in_fam = 0
  this%max_distance = 5.0
  this%intensity_fracfam= UNINITIALIZED_DOUBLE
  this%surface_area_fracfam = 0.d0
  this%hap0 = UNINITIALIZED_DOUBLE
  this%hap0_stdev = UNINITIALIZED_DOUBLE
  this%hap0_seed = 5
  this%hap0_max = UNINITIALIZED_DOUBLE
  this%center%x = UNINITIALIZED_DOUBLE
  this%center%y = UNINITIALIZED_DOUBLE
  this%center%z = UNINITIALIZED_DOUBLE
  this%center_stdev%x = UNINITIALIZED_DOUBLE
  this%center_stdev%y = UNINITIALIZED_DOUBLE
  this%center_stdev%z = UNINITIALIZED_DOUBLE
  this%center_seed = 10
  this%normal%x = UNINITIALIZED_DOUBLE
  this%normal%y = UNINITIALIZED_DOUBLE
  this%normal%z = UNINITIALIZED_DOUBLE
  this%normal_stdev%x = UNINITIALIZED_DOUBLE
  this%normal_stdev%y = UNINITIALIZED_DOUBLE
  this%normal_stdev%z = UNINITIALIZED_DOUBLE
  this%normal_seed = 20
  this%radius%x = UNINITIALIZED_DOUBLE
  this%radius%y = UNINITIALIZED_DOUBLE
  this%radius%z = UNINITIALIZED_DOUBLE
  this%radius_stdev%x = UNINITIALIZED_DOUBLE
  this%radius_stdev%y = UNINITIALIZED_DOUBLE
  this%radius_stdev%z = UNINITIALIZED_DOUBLE
  this%radius_seed = 30

end subroutine FractureFamInit

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
  type(fracture_type), pointer :: cur_fracture,new_fracture_list
  type(fracfam_type), pointer :: cur_fracfam
  type(grid_type), pointer :: res_grid
  character(len=MAXWORDLENGTH) :: word
  PetscInt, pointer :: temp_cell_ids(:)
  PetscInt, pointer :: temp_allfrac_cell_ids(:)
  PetscInt, pointer :: temp_frac_common_cell_ids(:)
  PetscInt :: k,j
  PetscInt :: nf,tfc,read_max
  PetscReal :: D,distance
  PetscReal :: min_x,max_x,min_y,max_y,min_z,max_z
  PetscReal :: a1,a2,a3,b1,b2,b3,c1,c2,c3,vmag
  PetscReal :: ra, sinra, cosra
  PetscBool :: within_x,within_y,within_z
  PetscBool :: added

  call this%SetRealization()

  option => this%option
  res_grid => this%realization%patch%grid
  nf = 0
  read_max = 5000 ! Set to a large but realistic number

  allocate(temp_cell_ids(read_max))
  allocate(temp_allfrac_cell_ids(read_max*10))
  allocate(temp_frac_common_cell_ids(read_max))

  cur_fracfam => this%fracfam_list
  do
    if (.not.associated(cur_fracfam)) exit

    ! generate fractures in this fracture family
    call PMFracGenerateFracFam(cur_fracfam)

    cur_fracfam => cur_fracfam%next
  enddo

  added = PETSC_FALSE
  cur_fracfam => this%fracfam_list
  do
    if (.not.associated(cur_fracfam)) exit
    new_fracture_list => cur_fracfam%fracture_list

    ! add fracture family fractures to the individual fracture list
    if (.not.associated(this%fracture_list)) then
      this%fracture_list => new_fracture_list
    else
      cur_fracture => this%fracture_list
      do
        if (.not.associated(cur_fracture)) exit
        if (.not.associated(cur_fracture%next)) then
          cur_fracture%next => new_fracture_list
          added = PETSC_TRUE
        endif
        if (added) exit
        cur_fracture => cur_fracture%next
      enddo
    endif

    cur_fracfam => cur_fracfam%next
    added = PETSC_FALSE
  enddo

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

  tfc = 0
  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    write(word,'(i4)') cur_fracture%id
    option%io_buffer = 'GEOTHERMAL_FRACTURE_MODEL: Mapping fracture ID# [' &
                       // trim(word) // '].'
    call PrintMsg(option)

    ! Equation of a plane normal to vector (A,B,C) is
    ! Ax + By + Cz + D = 0
    ! Get D by knowing the plane must contain the center point (x,y,z)
    D = -1.d0*(cur_fracture%normal%x*cur_fracture%center%x + &
             cur_fracture%normal%y*cur_fracture%center%y + &
             cur_fracture%normal%z*cur_fracture%center%z)
    j = 0
    do k = 1,res_grid%nlmax
      distance = cur_fracture%normal%x*res_grid%x(res_grid%nL2G(k)) + &
                 cur_fracture%normal%y*res_grid%y(res_grid%nL2G(k)) + &
                 cur_fracture%normal%z*res_grid%z(res_grid%nL2G(k)) + D
      if (abs(distance) < cur_fracture%max_distance) then
        ! entered here if the center of the grid cell is within a tolerance
        ! of max_distance of the plane.
        ! next check if the (x,y) position of the grid cell lies within the
        ! requested rad(x,y) of the plane from the center point.
        min_x = cur_fracture%center%x - cur_fracture%radx
        max_x = cur_fracture%center%x + cur_fracture%radx
        min_y = cur_fracture%center%y - cur_fracture%rady
        max_y = cur_fracture%center%y + cur_fracture%rady
        min_z = cur_fracture%center%z - cur_fracture%radz
        max_z = cur_fracture%center%z + cur_fracture%radz
        within_x = PETSC_FALSE; within_y = PETSC_FALSE
        within_z = PETSC_FALSE
        if ((res_grid%x(res_grid%nL2G(k)) <= max_x) .and. &
            (res_grid%x(res_grid%nL2G(k)) >= min_x)) then
          within_x = PETSC_TRUE
        endif
        if ((res_grid%y(res_grid%nL2G(k)) <= max_y) .and. &
            (res_grid%y(res_grid%nL2G(k)) >= min_y)) then
          within_y = PETSC_TRUE
        endif
        if ((res_grid%z(res_grid%nL2G(k)) <= max_z) .and. &
            (res_grid%z(res_grid%nL2G(k)) >= min_z)) then
          within_z = PETSC_TRUE
        endif
        ! if you want the fracture oriented exactly along an axis, then allow
        ! that axis to pass the following logic statement by setting the
        ! within_* booleans to TRUE
        if (abs(max_x-min_x) < 1.d-5) within_x = PETSC_TRUE
        if (abs(max_y-min_y) < 1.d-5) within_y = PETSC_TRUE
        if (abs(max_z-min_z) < 1.d-5) within_z = PETSC_TRUE
        if (within_x .and. within_y .and. within_z) then
          j = j + 1
          temp_cell_ids(j) = k
        endif
      endif
    enddo
    cur_fracture%ncells = j
    if (j > read_max) then
      option%io_buffer = 'The number of grid cells that FRACTURE with &
        &ID# [' // word // '] occupies exceeds the set limit of 5000 cells. If &
        &you are certain the fracture should be this large, increase the size &
        &of read_max in PMFracSetup() or email the PFLOTRAN developers group.'
      call PrintErrMsg(option)
    endif
    if (cur_fracture%ncells > 0) then
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
      temp_allfrac_cell_ids(tfc+1:tfc+cur_fracture%ncells) = &
                                  cur_fracture%cell_ids(1:cur_fracture%ncells)
    endif

    tfc = tfc + cur_fracture%ncells
    if (tfc > (read_max*10)) then
      option%io_buffer = 'tfc exceed maximum. Email the PFLOTRAN developers &
        &group, or increase the size of read_max*10 in PMFracSetup().'
      call PrintErrMsg(option)
    endif

    ! get the unit vector parallel to the fracture plane
    ! a X b = c, where c is perpendicular to a and b (a and b can't be parallel)
    ! take a as the unit normal, and b as slightly rotated unit normal
    ! then c will be perpendicular to a, and parallel to fracture plane
    vmag = sqrt((cur_fracture%normal%x)**2.d0 + &
                (cur_fracture%normal%y)**2.d0 + &
                (cur_fracture%normal%z)**2.d0)
    cur_fracture%normal%x = cur_fracture%normal%x/vmag
    cur_fracture%normal%y = cur_fracture%normal%y/vmag
    cur_fracture%normal%z = cur_fracture%normal%z/vmag
    a1 = cur_fracture%normal%x; b1 = cur_fracture%normal%x + 0.25
    a2 = cur_fracture%normal%y; b2 = cur_fracture%normal%y + 0.5
    a3 = cur_fracture%normal%z; b3 = cur_fracture%normal%z + 1.0
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

  if (tfc > 0) then
    allocate(this%allfrac_cell_ids(tfc))
    this%allfrac_cell_ids(1:tfc) = temp_allfrac_cell_ids(1:tfc)
    deallocate(temp_allfrac_cell_ids)
    allocate(temp_allfrac_cell_ids(tfc))
    temp_allfrac_cell_ids(:) = UNINITIALIZED_INTEGER

    deallocate(temp_cell_ids)
    allocate(temp_cell_ids(tfc))
    temp_cell_ids(:) = UNINITIALIZED_INTEGER
    j = 1
    do k = 1,tfc
      if (any(temp_allfrac_cell_ids == this%allfrac_cell_ids(k))) then
        temp_cell_ids(j) = this%allfrac_cell_ids(k)
        j = j + 1
      endif
      temp_allfrac_cell_ids(k) = this%allfrac_cell_ids(k)
    enddo
    allocate(this%frac_common_cell_ids(j))
    this%frac_common_cell_ids(1:j) = temp_cell_ids(1:j)
  deallocate(temp_frac_common_cell_ids)
  endif

  if (.not.associated(this%realization%reaction)) then
    this%update_material_auxvar_perm = PETSC_TRUE
  else
    if (.not.this%realization%reaction%update_permeability) then
      this%update_material_auxvar_perm = PETSC_TRUE
    endif
  endif

  deallocate(temp_cell_ids)
  deallocate(temp_allfrac_cell_ids)

end subroutine PMFracSetup

! ************************************************************************** !

subroutine PMFracSetRealization(this)
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  use Realization_Subsurface_class

  implicit none

  class(pm_fracture_type) :: this

  this%realization => RealizationCast(this%realization_base)

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
  use General_Aux_module
  use Field_module
  use Grid_module

  implicit none

  class(pm_fracture_type) :: this

  type(material_auxvar_type), pointer :: material_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(general_auxvar_type), pointer :: general_auxvar
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(fracture_type), pointer :: cur_fracture
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscInt :: icell,k
  PetscReal :: kx,ky,kz,L
  PetscErrorCode :: ierr

  option => this%option
  field => this%realization%field
  grid => this%realization%patch%grid

  call VecGetArrayReadF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    do k = 1,cur_fracture%ncells
      icell = cur_fracture%cell_ids(k)
      material_auxvar => &
        this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
      ! if cells aren't too pancaked, L should be a good estimate of length
      L = (material_auxvar%volume)**(1.d0/3.d0)

      cur_fracture%kx0(k) = perm0_xx_p(icell)
      cur_fracture%ky0(k) = perm0_yy_p(icell)
      cur_fracture%kz0(k) = perm0_zz_p(icell)

      call PMFracCalcK(L,cur_fracture%hap0,cur_fracture%RM,kx,ky,kz)

      material_auxvar%permeability(1) = material_auxvar%permeability(1) + kx
      material_auxvar%permeability(2) = material_auxvar%permeability(2) + ky
      material_auxvar%permeability(3) = material_auxvar%permeability(3) + kz

      material_auxvar%porosity_0 = material_auxvar%porosity_0 + &
                                   (cur_fracture%hap0/L)
      material_auxvar%porosity_base = material_auxvar%porosity_base + &
                                      (cur_fracture%hap0/L)
      material_auxvar%porosity = material_auxvar%porosity + &
                                 (cur_fracture%hap0/L)

      if (option%iflowmode == G_MODE) then
        general_auxvar => &
          this%realization%patch%aux%General%auxvars(1,grid%nL2G(icell))
        cur_fracture%prev_temperature(k) = general_auxvar%temp
      elseif (option%iflowmode == TH_MODE) then
        global_auxvar => &
          this%realization%patch%aux%Global%auxvars(grid%nL2G(icell))
        cur_fracture%prev_temperature(k) = global_auxvar%temp
      endif
    enddo

    cur_fracture => cur_fracture%next
  enddo

  if (associated(this%allfrac_cell_ids)) then
    do k = 1,size(this%allfrac_cell_ids)
      icell = this%allfrac_cell_ids(k)
      material_auxvar => &
          this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
      perm0_xx_p(icell) = material_auxvar%permeability(1)
      perm0_yy_p(icell) = material_auxvar%permeability(2)
      perm0_zz_p(icell) = material_auxvar%permeability(3)
    enddo
  endif

  call VecRestoreArrayReadF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)

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
  use Material_Aux_module
  use Field_module
  use Grid_module

  implicit none

  class(pm_fracture_type) :: this

  type(material_auxvar_type), pointer :: material_auxvar
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(fracture_type), pointer :: cur_fracture
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscInt :: icell,k
  PetscErrorCode :: ierr

  field => this%realization%field
  grid => this%realization%patch%grid

  call VecGetArrayReadF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    do k = 1,cur_fracture%ncells
      icell = cur_fracture%cell_ids(k)
      material_auxvar => &
        this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
      ! reset permeability back to initial permeability
      perm0_xx_p(icell) = cur_fracture%kx0(k)
      perm0_yy_p(icell) = cur_fracture%ky0(k)
      perm0_zz_p(icell) = cur_fracture%kz0(k)
    enddo

    cur_fracture => cur_fracture%next
  enddo

  call VecRestoreArrayReadF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)

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
  PetscBool :: found,fracfam_given,frac_given

  option => this%option
  fracfam_given = PETSC_FALSE
  frac_given = PETSC_FALSE
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
    !-------------------------------------
      case('THERMAL_EXPANSION_COEFFICIENT')
        call InputReadDouble(input,option,this%t_coeff)
        call InputErrorMsg(input,option,'THERMAL_EXPANSION_COEFFICIENT', &
                         error_string)
        cycle
    !-------------------------------------
      case('MAXIMUM_NUMBER_OF_FRACTURES')
        call InputReadInt(input,option,this%max_frac)
        call InputErrorMsg(input,option,'MAXIMUM_NUMBER_OF_FRACTURES', &
                         error_string)
        cycle
    !-------------------------------------

    end select

    ! Read sub-blocks within GEOTHERMAL_FRACTURE_MODEL block:
    error_string = 'GEOTHERMAL_FRACTURE_MODEL'
    call PMFracReadFracture(this,input,option,word,error_string,found)
    if (found) cycle

    error_string = 'GEOTHERMAL_FRACTURE_MODEL'
    call PMFracReadFracFam(this,input,option,word,error_string,found)
    if (found) cycle

    if (.not. found) then
      option%io_buffer = 'Keyword "' // trim(word) // &
                         '" does not exist for GEOTHERMAL_FRACTURE_MODEL.'
      call PrintErrMsg(option)
    endif

  enddo
  call InputPopBlock(input,option)

if (Uninitialized(this%t_coeff)) then
  option%io_buffer = 'THERMAL_EXPANSION_COEFFICIENT not provided in &
    &the GEOTHERMAL_FRACTURE_MODEL block.'
  call PrintErrMsg(option)
endif

if (Uninitialized(this%max_frac)) then
  this%max_frac = 100000000  ! set to a very high number so the max
                             ! is never reached
endif

if (associated(this%fracfam_list)) fracfam_given = PETSC_TRUE
if (associated(this%fracture_list)) frac_given = PETSC_TRUE

if (.not.fracfam_given .and. .not.frac_given) then
  option%io_buffer = 'At least one FRACTURE or FRACTURE_FAMILY block must &
    &be given in the GEOTHERMAL_FRACTURE_MODEL block.'
  call PrintErrMsg(option)
endif

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
  type(fracture_type), pointer :: new_fracture,cur_fracture
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
      new_fracture => FractureCreate()
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
          case('MAX_DISTANCE')
            call InputReadDouble(input,option,new_fracture%max_distance)
            call InputErrorMsg(input,option,'MAX_DISTANCE',error_string)
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
          case('RADIUS_Z')
            call InputReadDouble(input,option,new_fracture%radz)
            call InputErrorMsg(input,option,'RADIUS_Z',error_string)
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
      if (Uninitialized(new_fracture%radz)) then
        option%io_buffer = 'ERROR: RADIUS_Z must be specified in ' // &
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

subroutine PMFracReadFracFam(pm_fracture,input,option,keyword,error_string, &
                            found)
  !
  ! Reads input file parameters associated with the fracture model grid.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/19/2024
  !

  implicit none

  class(pm_fracture_type) :: pm_fracture
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  type(fracfam_type), pointer :: new_fracfam, cur_fracfam
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: added
  PetscInt :: num_errors

  error_string = trim(error_string) // ',FRACTURE_FAMILY'
  found = PETSC_TRUE
  added = PETSC_FALSE
  num_errors = 0

  select case(trim(keyword))
  !-------------------------------------
    case('FRACTURE_FAMILY')
      allocate(new_fracfam)
      new_fracfam => FractureFamCreate()
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
            call InputReadInt(input,option,new_fracfam%id)
            call InputErrorMsg(input,option,'ID',error_string)
        !-----------------------------
          case('MAX_DISTANCE')
            call InputReadDouble(input,option,new_fracfam%max_distance)
            call InputErrorMsg(input,option,'MAX_DISTANCE',error_string)
        !-----------------------------
          case('NUMBER_OF_FRACTURES')
            call InputReadInt(input,option,new_fracfam%intensity_fracfam)
            call InputErrorMsg(input,option,'NUMBER_OF_FRACTURES',error_string)
        !-----------------------------
          case('HYDRAULIC_APERTURE')
            call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'keyword',error_string)
                call StringToUpper(word)
                select case(trim(word))
                !----------------------------------
                  case('HYDRAULIC_APERTURE_VALUE')
                    call InputReadDouble(input,option,new_fracfam%hap0)
                    call InputErrorMsg(input,option, &
                                       'HYDRAULIC_APERTURE_VALUE',error_string)
                !----------------------------------
                  case('HYDRAULIC_APERTURE_MAX')
                    call InputReadDouble(input,option,new_fracfam%hap0_max)
                    call InputErrorMsg(input,option, &
                                       'HYDRAULIC_APERTURE_MAX',error_string)
                !----------------------------------
                  case('HYDRAULIC_APERTURE_STDEV')
                    call InputReadDouble(input,option,new_fracfam%hap0_stdev)
                    call InputErrorMsg(input,option, &
                                       'HYDRAULIC_APERTURE_STDEV',error_string)
                !-----------------------------
                  case('HYDRAULIC_APERTURE_SEED')
                    call InputReadInt(input,option,new_fracfam%hap0_seed)
                    call InputErrorMsg(input,option,'HYDRAULIC_APERTURE_SEED', &
                                       error_string)
                !-----------------------------------
                  case default
                    call InputKeywordUnrecognized(input,word, &
                                                  error_string,option)
                !-----------------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case('CENTER')
            call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'keyword',error_string)
                call StringToUpper(word)
                select case(trim(word))
                !----------------------------------
                  case('COORDINATE')
                    call InputReadDouble(input,option,new_fracfam%center%x)
                    call InputErrorMsg(input,option,'CENTER,X COORDINATE', &
                                       error_string)
                    call InputReadDouble(input,option,new_fracfam%center%y)
                    call InputErrorMsg(input,option,'CENTER,Y COORDINATE', &
                                       error_string)
                    call InputReadDouble(input,option,new_fracfam%center%z)
                    call InputErrorMsg(input,option,'CENTER,Z COORDINATE', &
                                       error_string)
                !----------------------------------
                  case('XCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%center_stdev%x)
                    call InputErrorMsg(input,option,'CENTER,XCOORD_STDEV', &
                                       error_string)
                !----------------------------------
                  case('YCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%center_stdev%y)
                    call InputErrorMsg(input,option,'CENTER,YCOORD_STDEV', &
                                       error_string)
                !----------------------------------
                  case('ZCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%center_stdev%z)
                    call InputErrorMsg(input,option,'CENTER,ZCOORD_STDEV', &
                                       error_string)
                !-----------------------------
                  case('CENTER_SEED')
                    call InputReadInt(input,option,new_fracfam%center_seed)
                    call InputErrorMsg(input,option,'CENTER_SEED',error_string)
                !-----------------------------------
                  case default
                    call InputKeywordUnrecognized(input,word, &
                                                  error_string,option)
                !-----------------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case('NORMAL_VECTOR')
            call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'keyword',error_string)
                call StringToUpper(word)
                select case(trim(word))
                !----------------------------------
                  case('VECTOR_COORDINATES')
                    call InputReadDouble(input,option,new_fracfam%normal%x)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,X COORDINATE',error_string)
                    call InputReadDouble(input,option,new_fracfam%normal%y)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,Y COORDINATE',error_string)
                    call InputReadDouble(input,option,new_fracfam%normal%z)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,Z COORDINATE',error_string)
                !----------------------------------
                  case('XCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%normal_stdev%x)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,XCOORD_STDEV',error_string)
                !----------------------------------
                  case('YCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%normal_stdev%y)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,YCOORD_STDEV',error_string)
                !----------------------------------
                  case('ZCOORD_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%normal_stdev%z)
                    call InputErrorMsg(input,option, &
                      'NORMAL_VECTOR,ZCOORD_STDEV',error_string)
                !-----------------------------
                  case('NORMAL_SEED')
                    call InputReadInt(input,option,new_fracfam%normal_seed)
                    call InputErrorMsg(input,option,'NORMAL_SEED',error_string)
                !-----------------------------------
                  case default
                    call InputKeywordUnrecognized(input,word, &
                                                  error_string,option)
                !-----------------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case('RADIUS')
            call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'keyword',error_string)
                call StringToUpper(word)
                select case(trim(word))
                !----------------------------------
                  case('RADIUS_XYZ')
                    call InputReadDouble(input,option,new_fracfam%radius%x)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,X',error_string)
                    call InputReadDouble(input,option,new_fracfam%radius%y)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,Y',error_string)
                    call InputReadDouble(input,option,new_fracfam%radius%z)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,Z',error_string)
                !----------------------------------
                  case('RAD_X_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%radius_stdev%x)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,RAD_X_STDEV', &
                                       error_string)
                !----------------------------------
                  case('RAD_Y_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%radius_stdev%y)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,RAD_Y_STDEV', &
                                       error_string)
                !----------------------------------
                  case('RAD_Z_STDEV')
                    call InputReadDouble(input,option, &
                                         new_fracfam%radius_stdev%z)
                    call InputErrorMsg(input,option,'RADIUS_XYZ,RAD_Z_STDEV', &
                                       error_string)
                !-----------------------------
                  case('RADIUS_SEED')
                    call InputReadInt(input,option,new_fracfam%radius_seed)
                    call InputErrorMsg(input,option,'RADIUS_SEED',error_string)
                !-----------------------------------
                  case default
                    call InputKeywordUnrecognized(input,word, &
                                                  error_string,option)
                !-----------------------------------
                end select
              enddo
            call InputPopBlock(input,option)
        !-----------------------------
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        !-----------------------------
        end select
      enddo
      call InputPopBlock(input,option)

      !------ error messaging ----------------------------------------
      if (Uninitialized(new_fracfam%id)) then
        option%io_buffer = 'ERROR: ID must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%intensity_fracfam)) then
        option%io_buffer = 'ERROR: FRACTURE_INTENSITY must be specified &
          &in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%hap0)) then
        option%io_buffer = 'ERROR: HYDRAULIC_APERTURE,HYDRAULIC_APERTURE_VALUE &
        &must be specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%hap0_max)) then
        option%io_buffer = 'ERROR: HYDRAULIC_APERTURE,HYDRAULIC_APERTURE_MAX &
        &must be specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%hap0_stdev)) then
        option%io_buffer = 'ERROR: HYDRAULIC_APERTURE,HYDRAULIC_APERTURE_STDEV &
        &must be specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%center%x) .or. &
          Uninitialized(new_fracfam%center%y) .or. &
          Uninitialized(new_fracfam%center%z)) then
        option%io_buffer = 'ERROR: CENTER coordinate must be specified in ' // &
                           trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%center_stdev%x) .or. &
          Uninitialized(new_fracfam%center_stdev%y) .or. &
          Uninitialized(new_fracfam%center_stdev%z)) then
        option%io_buffer = 'ERROR: CENTER,(X/Y/Z)COORD_STDEV must be &
          &specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%normal%x) .or. &
          Uninitialized(new_fracfam%normal%y) .or. &
          Uninitialized(new_fracfam%normal%z)) then
        option%io_buffer = 'ERROR: NORMAL_VECTOR coordinate must be &
                           &specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%normal_stdev%x) .or. &
          Uninitialized(new_fracfam%normal_stdev%y) .or. &
          Uninitialized(new_fracfam%normal_stdev%z)) then
        option%io_buffer = 'ERROR: NORMAL,(X/Y/Z)COORD_STDEV must be &
          &specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%radius%x) .or. &
          Uninitialized(new_fracfam%radius%y) .or. &
          Uninitialized(new_fracfam%radius%z)) then
        option%io_buffer = 'ERROR: RADIUS,RADIUS_XYZ must be &
                           &specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif
      if (Uninitialized(new_fracfam%radius_stdev%x) .or. &
          Uninitialized(new_fracfam%radius_stdev%y) .or. &
          Uninitialized(new_fracfam%radius_stdev%z)) then
        option%io_buffer = 'ERROR: RADIUS,RAD_(X/Y/Z)_STDEV must be &
          &specified in ' // trim(error_string) // ' block.'
        call PrintMsg(option); num_errors = num_errors + 1
      endif

      if (num_errors > 0) then
        write(option%io_buffer,*) num_errors
        option%io_buffer = trim(adjustl(option%io_buffer)) // ' errors in &
                           &the FRACTURE_FAMILY block. See above to fix.'
        call PrintErrMsg(option)
      endif
      !------ add fracture family to list -----------------------------------
      if (.not.associated(pm_fracture%fracfam_list)) then
        pm_fracture%fracfam_list => new_fracfam
      else
        cur_fracfam => pm_fracture%fracfam_list
        do
          if (.not.associated(cur_fracfam)) exit
          if (.not.associated(cur_fracfam%next)) then
            cur_fracfam%next => new_fracfam
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_fracfam => cur_fracfam%next
        enddo
      endif
      nullify(new_fracfam)
  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

end subroutine PMFracReadFracFam

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
      case('FRACTURE_FAMILY')
        call InputSkipToEND(input,option,card)
        call InputSkipToEND(input,option,card)
        call InputSkipToEND(input,option,card)
        call InputSkipToEND(input,option,card)
        call InputSkipToEND(input,option,card)
      !--------------------
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine PMFracReadPass2

! ************************************************************************** !

subroutine PMFracGenerateFracFam(fracfam)
  !
  ! Generates fractures based on the information about a fracture family.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 08/20/2024
  !
  use Utility_module

  implicit none

  type(fracfam_type) :: fracfam

  type(fracture_type), pointer :: new_fracture,cur_fracture
  PetscBool :: added

  do while (fracfam%nfrac_in_fam < fracfam%intensity_fracfam)
    added = PETSC_FALSE
    allocate(new_fracture)
    fracfam%nfrac_in_fam = fracfam%nfrac_in_fam + 1
    new_fracture => FractureCreate()
    new_fracture%id = fracfam%id*100 + (fracfam%nfrac_in_fam)
    new_fracture%max_distance = fracfam%max_distance

    call GetRndNumFromNormalDist(fracfam%center%x,fracfam%center_stdev%x, &
      new_fracture%center%x,(fracfam%center_seed + fracfam%nfrac_in_fam))
    call GetRndNumFromNormalDist(fracfam%center%y,fracfam%center_stdev%y, &
      new_fracture%center%y,(fracfam%center_seed + fracfam%nfrac_in_fam))
    call GetRndNumFromNormalDist(fracfam%center%z,fracfam%center_stdev%z, &
      new_fracture%center%z,(fracfam%center_seed + fracfam%nfrac_in_fam))

    call GetRndNumFromNormalDist(fracfam%normal%x,fracfam%normal_stdev%x, &
      new_fracture%normal%x,(fracfam%normal_seed + fracfam%nfrac_in_fam))
    call GetRndNumFromNormalDist(fracfam%normal%y,fracfam%normal_stdev%y, &
      new_fracture%normal%y,(fracfam%normal_seed + fracfam%nfrac_in_fam))
    call GetRndNumFromNormalDist(fracfam%normal%z,fracfam%normal_stdev%z, &
      new_fracture%normal%z,(fracfam%normal_seed + fracfam%nfrac_in_fam))

    call GetRndNumFromNormalDist(fracfam%radius%x,fracfam%radius_stdev%x, &
      new_fracture%radx,(fracfam%radius_seed + fracfam%nfrac_in_fam))
    call GetRndNumFromNormalDist(fracfam%radius%y,fracfam%radius_stdev%y, &
      new_fracture%rady,(fracfam%normal_seed + fracfam%nfrac_in_fam))
    call GetRndNumFromNormalDist(fracfam%radius%z,fracfam%radius_stdev%z, &
      new_fracture%radz,(fracfam%normal_seed + fracfam%nfrac_in_fam))
    new_fracture%radx = abs(new_fracture%radx)  ! ensures radius is positive
    new_fracture%rady = abs(new_fracture%rady)  ! ensures radius is positive
    new_fracture%radz = abs(new_fracture%radz)  ! ensures radius is positive

    call GetRndNumFromNormalDist(fracfam%hap0,fracfam%hap0_stdev, &
      new_fracture%hap0,(fracfam%hap0_seed + fracfam%nfrac_in_fam))
    new_fracture%hap0 = abs(new_fracture%hap0)  ! ensures hap0 is positive
    new_fracture%hap0 = min(new_fracture%hap0,fracfam%hap0_max)

    !------ add fracture to list -------------------------
    if (.not.associated(fracfam%fracture_list)) then
      fracfam%fracture_list => new_fracture
    else
      cur_fracture => fracfam%fracture_list
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
  enddo

end subroutine PMFracGenerateFracFam

! ************************************************************************** !

subroutine PMFracCalc_dK(L,d_hap,RM,kx,ky,kz)
  !
  ! Calculates change in cell anisotropic permeability from fracture
  ! permeability.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/18/2023
  !

  implicit none

  PetscReal :: L
  PetscReal :: d_hap ! change in hap
  PetscReal :: RM(3,3)
  PetscReal :: kx,ky,kz

  PetscReal :: K_frac(3,3)
  PetscReal :: K_frac_rot(3,3)
  PetscReal :: K_domain(3,3)

  ! the math:
  ! K = (1.d0/12.d0)*(hap**3.0d0)/L
  ! d_K/d_hap = (3.d0/12.d0)*(hap**2.0d0)/L

  K_frac(:,:) = 0.d0
  K_frac(1,1) = (3.d0/12.d0)*(d_hap**2.0d0)/L ! kxx
  K_frac(2,2) = (3.d0/12.d0)*(d_hap**2.0d0)/L ! kyy
  ! note: zero in kzz

  K_frac_rot = MATMUL(RM,K_frac)
  K_domain = MATMUL(K_frac_rot,TRANSPOSE(RM))

  kx = K_domain(1,1)
  ky = K_domain(2,2)
  kz = K_domain(3,3)

end subroutine PMFracCalc_dK

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

subroutine PMFracSolve(this,time,ierr)
  !
  ! Main solve step for the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/13/2023
  !
  use Material_Aux_module
  use Global_Aux_module
  use General_Aux_module
  use Field_module
  use Grid_module

  implicit none

  class(pm_fracture_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(material_auxvar_type), pointer :: material_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(general_auxvar_type), pointer :: general_auxvar
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(fracture_type), pointer :: cur_fracture
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscReal :: cur_temperature,L
  PetscReal :: kx,ky,kz
  PetscReal :: dkx,dky,dkz
  PetscInt :: icell,k

  option => this%option
  field => this%realization%field
  grid => this%realization%patch%grid

  call VecGetArrayReadF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    do k = 1,cur_fracture%ncells
      icell = cur_fracture%cell_ids(k)
      material_auxvar => &
        this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
      if (option%iflowmode == G_MODE) then
        general_auxvar => &
          this%realization%patch%aux%General%auxvars(1,grid%nL2G(icell))
        cur_temperature = general_auxvar%temp
      elseif (option%iflowmode == TH_MODE) then
        global_auxvar => &
          this%realization%patch%aux%Global%auxvars(grid%nL2G(icell))
        cur_temperature = global_auxvar%temp
      endif

      if (cur_fracture%prev_temperature(k) == 0.d0) then
        ! for whatever reason, in general_mode the auxvar temperature is not
        ! initialized before the fracture pm calls initialize run.
        cur_fracture%dT(k) = 0.d0
      else
        cur_fracture%dT(k) = cur_temperature - cur_fracture%prev_temperature(k)
      endif

      ! if cells aren't too pancaked, L should be a good estimate of length
      L = (material_auxvar%volume)**(1.d0/3.d0)
      cur_fracture%dL(k) = L*cur_fracture%dT(k)*this%t_coeff

      ! update the aperture with dL and calculate new perm
      cur_fracture%hap(k) = cur_fracture%hap(k) - cur_fracture%dL(k)
      ! make sure aperature doesn't go negative
      cur_fracture%hap(k) = max(0.d0,cur_fracture%hap(k))

      ! get change in anisotropic domain permeability from rotation transform
      call PMFracCalc_dK(L,-1.d0*cur_fracture%dL(k),cur_fracture%RM, &
                         dkx,dky,dkz)

      ! get anisotropic domain permeability from rotation transformation
      call PMFracCalcK(L,cur_fracture%hap(k),cur_fracture%RM,kx,ky,kz)

      ! update domain material permeability of cells with fractures
      perm0_xx_p(icell) = perm0_xx_p(icell) + kx
      perm0_yy_p(icell) = perm0_yy_p(icell) + ky
      perm0_zz_p(icell) = perm0_zz_p(icell) + kz

      ! reset previous temperature for next time step
      cur_fracture%prev_temperature(k) = cur_temperature
    enddo

    cur_fracture => cur_fracture%next
  enddo

  if (this%update_material_auxvar_perm) then
    if (associated(this%allfrac_cell_ids)) then
      do k = 1,size(this%allfrac_cell_ids)
        icell = this%allfrac_cell_ids(k)
        material_auxvar => &
          this%realization%patch%aux%material%auxvars(grid%nL2G(icell))
        material_auxvar%permeability(1) = perm0_xx_p(icell)
        material_auxvar%permeability(2) = perm0_yy_p(icell)
        material_auxvar%permeability(3) = perm0_zz_p(icell)
      enddo
    endif
  endif

  call VecRestoreArrayReadF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)

  ierr = 0 ! If this is not set to zero, TS_STOP_FAILURE occurs

end subroutine PMFracSolve

! ************************************************************************** !

subroutine PMFracDestroy(this)
  !
  ! Destroys objects in the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/18/2023
  !
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_fracture_type) :: this

  type(fracture_type), pointer :: cur_fracture

  if (associated(this%allfrac_cell_ids)) then
    call DeallocateArray(this%allfrac_cell_ids)
  endif
  if (associated(this%frac_common_cell_ids)) then
    call DeallocateArray(this%frac_common_cell_ids)
  endif

  cur_fracture => this%fracture_list
  do
    if (.not.associated(cur_fracture)) exit

    if (cur_fracture%ncells > 0) call PMFracStripFrac(cur_fracture)

    cur_fracture => cur_fracture%next
  enddo

end subroutine PMFracDestroy

! ************************************************************************** !

subroutine PMFracStripFrac(frac)
  !
  ! Destroys objects in the fracture process model.
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 12/18/2023
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(fracture_type), pointer :: frac

  call DeallocateArray(frac%kx0)
  call DeallocateArray(frac%ky0)
  call DeallocateArray(frac%kx0)
  call DeallocateArray(frac%hap)
  call DeallocateArray(frac%cell_ids)
  call DeallocateArray(frac%prev_temperature)
  call DeallocateArray(frac%dT)
  call DeallocateArray(frac%dL)

end subroutine PMFracStripFrac

! ************************************************************************** !

end module PM_Fracture_class
