module Geomechanics_Subsurface_Properties_module
  
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  private


  PetscInt, parameter, public :: Bandis_A_index = 1
  PetscInt, parameter, public :: Bandis_B_index = 2
  PetscInt, parameter, public :: maximum_aperture_index = 3
  PetscInt, parameter, public :: model_index = 4
  PetscInt, parameter, public :: normal_vector_x_index = 5
  PetscInt, parameter, public :: normal_vector_y_index = 6
  PetscInt, parameter, public :: normal_vector_z_index = 7
  
  PetscInt, parameter, public :: BANDIS_MODEL = 8
  PetscInt, parameter, public :: TURNER_MODEL = 9
  PetscInt, parameter, public :: LINEAR_MODEL = 10

  type, public :: geomechanics_subsurface_properties_type
    character(len=MAXWORDLENGTH) :: geomechanical_compressibility_function
    PetscReal :: Bandis_A
    PetscReal :: Bandis_B
    PetscReal :: maximum_aperture
    PetscReal :: normal_vector_x
    PetscReal :: normal_vector_y
    PetscReal :: normal_vector_z
  contains
    procedure, public :: Read => GeomechanicsSubsurfacePropsRead
  end type geomechanics_subsurface_properties_type
  
  public :: GeomechanicsSubsurfacePropsInit, &
            GeomechanicsSubsurfacePropsCreate, &
            GeomechanicsSubsurfacePropsAuxvarInit, &
            GeomechanicsSubsurfacePropsPropertytoAux, &
            GeomechanicsSubsurfacePropsDestroy, &
            GeomechanicsSubsurfacePropsPoroEvaluate, &
            GeomechanicsSubsurfacePropsPermEvaluate
  
  contains

! ************************************************************************** !

function GeomechanicsSubsurfacePropsCreate()
  !
  ! Author: Satish Karra
  ! Date: 07/29/16
  !

  implicit none
  
  class(geomechanics_subsurface_properties_type), pointer :: &  
    GeomechanicsSubsurfacePropsCreate
  class(geomechanics_subsurface_properties_type), pointer :: &
    this 
  
  allocate(this)
  call GeomechanicsSubsurfacePropsInit(this)
  
  GeomechanicsSubsurfacePropsCreate => this 
  
end function GeomechanicsSubsurfacePropsCreate

! ************************************************************************** !

subroutine GeomechanicsSubsurfacePropsInit(this)
  !
  ! Author: Satish Karra
  ! Date: 07/29/16
  !

  implicit none
  
  class(geomechanics_subsurface_properties_type), pointer :: this 

  this%geomechanical_compressibility_function = ''  
  this%Bandis_A = UNINITIALIZED_DOUBLE
  this%Bandis_B = UNINITIALIZED_DOUBLE
  this%maximum_aperture = UNINITIALIZED_DOUBLE
  this%normal_vector_x = UNINITIALIZED_DOUBLE
  this%normal_vector_y = UNINITIALIZED_DOUBLE
  this%normal_vector_z = UNINITIALIZED_DOUBLE

end subroutine GeomechanicsSubsurfacePropsInit

! ************************************************************************** !

subroutine GeomechanicsSubsurfacePropsAuxvarInit( &
                                      this,auxvar)
  !
  ! Author: Satish Karra
  ! Date: 07/29/16
  !

  use Material_Aux_class
  
  implicit none
  
  class(geomechanics_subsurface_properties_type), pointer :: &
    this 
  class(material_auxvar_type), intent(inout) :: auxvar    
   
  allocate(auxvar%geomechanics_subsurface_prop(7))
  auxvar%geomechanics_subsurface_prop = 0.d0

end subroutine GeomechanicsSubsurfacePropsAuxvarInit

! ************************************************************************** !

subroutine GeomechanicsSubsurfacePropsPropertytoAux(auxvar,this)
  !
  ! Author: Satish Karra
  ! Date: 07/29/16
  !

  use Material_Aux_class
  use String_module
  use Option_module
  
  implicit none

  class(material_auxvar_type), intent(inout) :: auxvar
  class(geomechanics_subsurface_properties_type), pointer :: &
    this 
  type(option_type) :: option
 
  auxvar%geomechanics_subsurface_prop(Bandis_A_index) = &
    this%Bandis_A  
  auxvar%geomechanics_subsurface_prop(Bandis_B_index) = &
    this%Bandis_B
  auxvar%geomechanics_subsurface_prop(maximum_aperture_index) = &
    this%maximum_aperture

  ! Normal vector to the fracture/fault plane
  auxvar%geomechanics_subsurface_prop(normal_vector_x_index) = &
    this%normal_vector_x
  auxvar%geomechanics_subsurface_prop(normal_vector_y_index) = &
    this%normal_vector_y
  auxvar%geomechanics_subsurface_prop(normal_vector_z_index) = &
    this%normal_vector_z
  
  ! set the model index
  call StringToUpper(this%geomechanical_compressibility_function)
  select case(this%geomechanical_compressibility_function)
    case ('BANDIS')
      auxvar%geomechanics_subsurface_prop(model_index) = &
        BANDIS_MODEL
    case ('LINEAR')
      auxvar%geomechanics_subsurface_prop(model_index) = &
        LINEAR_MODEL
    case ('TURNER')
      auxvar%geomechanics_subsurface_prop(model_index) = &
        TURNER_MODEL
    case default
      auxvar%geomechanics_subsurface_prop(model_index) = &
        LINEAR_MODEL
  end select
 
end subroutine GeomechanicsSubsurfacePropsPropertytoAux

! ************************************************************************** !

subroutine GeomechanicsSubsurfacePropsRead(this,input,option)
  ! 
  ! Author: Satish Karra
  ! Date: 07/29/16
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  
  implicit none
  
  class(geomechanics_subsurface_properties_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  
  do
      call InputReadPflotranString(input,option)
      call InputReadStringErrorMsg(input,option, &
                        'MATERIAL_PROPERTY,GEOMECHANICS_SUBSURFACE_PROPS')
          
      if (InputCheckExit(input,option)) exit
          
      if (InputError(input)) exit
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'keyword', &
                          'MATERIAL_PROPERTY,GEOMECHANICS_SUBSURFACE_PROPS')   
      select case(trim(word))
        case('COMPRESSIBILITY_FUNCTION')
          call InputReadWord(input,option, &
                             this%geomechanical_compressibility_function, &
                             PETSC_TRUE)
          call InputErrorMsg(input,option, &
                             'geomechanical compressibility function', &
                             'GEOMECHANICS_SUBSURFACE_PROPS')
        case('BANDIS_A') 
          call InputReadDouble(input,option, &
                               this%Bandis_A)
          call InputErrorMsg(input,option,'Bandis A parameter', &
                             'GEOMECHANICS_SUBSURFACE_PROPS')
        case('BANDIS_B') 
          call InputReadDouble(input,option, &
                               this%Bandis_B)
          call InputErrorMsg(input,option,'Bandis B parameter', &
                             'GEOMECHANICS_SUBSURFACE_PROPS')
        case('MAXIMUM_APERTURE') 
          call InputReadDouble(input,option, &
                               this%maximum_aperture)
          call InputErrorMsg(input,option,'max aperture for Bandis Model', &
                             'GEOMECHANICS_SUBSURFACE_PROPS')
        case('NORMAL_VECTOR')
          call InputReadDouble(input,option,this%normal_vector_x)
          call InputErrorMsg(input,option,'x-direction','NORMAL_VECTOR')
          call InputReadDouble(input,option,this%normal_vector_y)
          call InputErrorMsg(input,option,'y-direction','NORMAL_VECTOR')
          call InputReadDouble(input,option,this%normal_vector_z)
          call InputErrorMsg(input,option,'z-direction','NORMAL_VECTOR')
        case default
          call InputKeywordUnrecognized(word, &
                  'MATERIAL_PROPERTY,GEOMECHANICS_SUBSURFACE_PROPS', &
                  option)
      end select
    enddo
    
end subroutine GeomechanicsSubsurfacePropsRead

! ************************************************************************** !

subroutine GeomechanicsSubsurfacePropsPoroEvaluate(grid, &
                                                   auxvar,porosity_before, &
                                                   local_stress, &
                                                   local_strain, &
                                                   local_pressure, &
                                                   porosity_after)
  !
  ! Calculates the change in porosity due to geomechanical strains
  !
  ! Author: Satish Karra
  ! Date: 07/29/16
  !

  use Option_module
  use Material_Aux_class
  use Grid_module

  
  implicit none
  
  type(option_type) :: option
  
  class(material_auxvar_type), intent(inout) :: auxvar
  type(grid_type), pointer, intent(inout) :: grid
  PetscReal, intent(in) :: porosity_before
  PetscReal, intent(in) :: local_stress(6), local_strain(6), local_pressure
  PetscReal, intent(out) :: porosity_after
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscReal :: Bandis_A, Bandis_B, maximum_aperture
  PetscReal :: normal_vector_x, normal_vector_y, normal_vector_z
  PetscInt :: model_id

  model_id = int(auxvar%geomechanics_subsurface_prop(model_index))

  if (model_id == 0) then
    model_id = LINEAR_MODEL ! set linear model by default if nothing is specified in the input file
  endif
        
  select case(model_id)
    case(BANDIS_MODEL)
      Bandis_A = auxvar%geomechanics_subsurface_prop(Bandis_A_index)
      Bandis_B = auxvar%geomechanics_subsurface_prop(Bandis_B_index)
      maximum_aperture = auxvar%geomechanics_subsurface_prop(maximum_aperture_index)
      normal_vector_x = auxvar%geomechanics_subsurface_prop(normal_vector_x_index)
      normal_vector_y = auxvar%geomechanics_subsurface_prop(normal_vector_y_index)
      normal_vector_z = auxvar%geomechanics_subsurface_prop(normal_vector_z_index)
      call GeomechanicsSubsurfaceBandisPoroEvaluate(grid,porosity_before, &
        local_stress,local_strain,local_pressure, &
        Bandis_A,Bandis_B,maximum_aperture,normal_vector_x,normal_vector_y, &
        normal_vector_z,porosity_after) 
    case(LINEAR_MODEL)
      call GeomechanicsSubsurfaceLinearPoroEvaluate(porosity_before, &
        local_stress,local_strain,local_pressure,porosity_after)    
    case(TURNER_MODEL)
      call GeomechanicsSubsurfaceTurnerPoroEvaluate(porosity_before, &
        local_stress,local_strain,local_pressure,porosity_after)
    case default
      write(string,*) model_id
      option%io_buffer = 'geomechanical compressibility model "' // &
        trim(string) // '" not recognized.'
      call printErrMsg(option)
    end select
          
end subroutine GeomechanicsSubsurfacePropsPoroEvaluate

! ************************************************************************** !

subroutine GeomechanicsSubsurfaceBandisPoroEvaluate(grid,porosity_before, &
                                  local_stress,local_strain,local_pressure, &
                                  Bandis_A,Bandis_B, &
                                  maximum_aperture,normal_vector_x, &
                                  normal_vector_y, &
                                  normal_vector_z,porosity_after) 
  ! 
  ! Calculates soil matrix compression for based on Bandis model (1983) 
  ! Citation: Bandis, S.C., Lumsden, A.C. and Barton, N.R., Fundamentals
  ! of Rock Joint Deformation, Int. J. Rock. Mech. Min. Sci. & Geomech. Abstr.
  ! Vol. 20, No. 6, pp. 249--268, 1983.
  !
  ! Author: Satish Karra
  ! Date: 07/29/16
  !

  use Grid_module

  implicit none
  
  type(grid_type), pointer, intent(inout) :: grid
  PetscReal, intent(in) :: porosity_before
  PetscReal, intent(in) :: local_stress(6), local_strain(6), local_pressure !DANNY-local stress is effective stress
  PetscReal, intent(in) :: Bandis_A, Bandis_B, maximum_aperture
  PetscReal, intent(in) :: normal_vector_x, normal_vector_y, normal_vector_z
  PetscReal, intent(out) :: porosity_after
  PetscReal :: effective_stress, b_p, b
  PetscInt :: icount

  effective_stress = local_stress(3)
  b_p = 1.d0  ! need to extract delta_x, delta_y, delta_z to get b_p
  b = maximum_aperture + Bandis_A*effective_stress/ &
            (1.d0 - Bandis_B*effective_stress)
  porosity_after = b/b_p

end subroutine GeomechanicsSubsurfaceBandisPoroEvaluate

! ************************************************************************** !

subroutine GeomechanicsSubsurfaceLinearPoroEvaluate(porosity_before, &
                                            local_stress,local_strain, &
                                            local_pressure,porosity_after)
  !
  ! Calculates the change in porosity due to geomechanical strains
  ! following linear model (needs citation)
  !
  ! Author: Satish Karra
  ! Date: 07/29/16
  !
  
  implicit none
  
  PetscReal, intent(in) :: porosity_before
  PetscReal, intent(in) :: local_stress(6), local_strain(6), local_pressure
  PetscReal, intent(out) :: porosity_after
  PetscReal :: volumetric_strain
  
  volumetric_strain = local_strain(1) + local_strain(2) + local_strain(3)
  porosity_after = porosity_before + volumetric_strain
  
end subroutine GeomechanicsSubsurfaceLinearPoroEvaluate

! ************************************************************************** !

subroutine GeomechanicsSubsurfaceTurnerPoroEvaluate(porosity_before, &
                                            local_stress,local_strain, &
                                            local_pressure,porosity_after)
  !
  ! Calculates the change in porosity due to geomechanical strains
  ! following Turner model (needs citation)
  !
  ! Author: Satish Karra
  ! Date: 07/29/16
  !
  
  implicit none
  
  PetscReal, intent(in) :: porosity_before
  PetscReal, intent(in) :: local_stress(6), local_strain(6), local_pressure
  PetscReal, intent(out) :: porosity_after
  PetscReal :: volumetric_strain
  
  volumetric_strain = local_strain(1) + local_strain(2) + local_strain(3)

  porosity_after = porosity_before/ &
      (1.d0 + (1.d0 - porosity_before)*volumetric_strain)

  
end subroutine GeomechanicsSubsurfaceTurnerPoroEvaluate

! ************************************************************************** !

subroutine GeomechanicsSubsurfacePropsPermEvaluate(grid, &
                                                   auxvar,permeability_before, &
                                                   local_stress, &
                                                   local_strain, &
                                                   local_pressure, &
                                                   permeability_after)
  !
  ! Calculates the change in permeability due to geomechanical strains
  ! according to Bandis calculation of b
  ! 
  ! Author: Daniel Birdsell, Satish Karra
  ! Date: 10/4/16
  !

  use Option_module
  use Material_Aux_class
  use Grid_module

  
  implicit none
  
  type(option_type) :: option
  
  class(material_auxvar_type), intent(inout) :: auxvar
  type(grid_type), pointer, intent(inout) :: grid
  PetscReal, intent(in) :: permeability_before
  PetscReal, intent(in) :: local_stress(6), local_strain(6), local_pressure
  PetscReal, intent(out) :: permeability_after
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscReal :: Bandis_A, Bandis_B, maximum_aperture
  PetscReal :: normal_vector_x, normal_vector_y, normal_vector_z
  PetscInt :: model_id
  
  model_id = int(auxvar%geomechanics_subsurface_prop(model_index))

  if (model_id == 0) then
    model_id = LINEAR_MODEL ! set linear model by default if nothing is specified in the input file
  endif
    
  select case(model_id)
    case(BANDIS_MODEL)
      Bandis_A = auxvar%geomechanics_subsurface_prop(Bandis_A_index)
      Bandis_B = auxvar%geomechanics_subsurface_prop(Bandis_B_index)
      maximum_aperture = auxvar%geomechanics_subsurface_prop(maximum_aperture_index)
      normal_vector_x = auxvar%geomechanics_subsurface_prop(normal_vector_x_index)
      normal_vector_y = auxvar%geomechanics_subsurface_prop(normal_vector_y_index)
      normal_vector_z = auxvar%geomechanics_subsurface_prop(normal_vector_z_index)
      call GeomechanicsSubsurfaceBandisPermEvaluate(grid,permeability_before, &
        local_stress,local_strain,local_pressure, &
        Bandis_A,Bandis_B,maximum_aperture,normal_vector_x,normal_vector_y, &
        normal_vector_z,permeability_after) 
    case(LINEAR_MODEL)
      call GeomechanicsSubsurfaceLinearPermEvaluate(permeability_before, &
        local_stress,local_strain,local_pressure,permeability_after)    
    case default
      write(string,*) model_id
      option%io_buffer = 'geomechanical perm model "' // &
      trim(string) // '" not recognized.'
      call printErrMsg(option)
  end select

end subroutine GeomechanicsSubsurfacePropsPermEvaluate

! ************************************************************************** !

subroutine GeomechanicsSubsurfaceBandisPermEvaluate(grid,permeability_before, &
                                  local_stress,local_strain,local_pressure, &
                                  Bandis_A,Bandis_B, &
                                  maximum_aperture,normal_vector_x, &
                                  normal_vector_y, &
                                  normal_vector_z,permeability_after) 
  ! 
  ! Calculates soil matrix compression for based on Bandis model (1983) 
  ! Citation: Bandis, S.C., Lumsden, A.C. and Barton, N.R., Fundamentals
  ! of Rock Joint Deformation, Int. J. Rock. Mech. Min. Sci. & Geomech. Abstr.
  ! Vol. 20, No. 6, pp. 249--268, 1983.
  !
  ! Author: Satish Karra
  ! Date: 07/29/16
  !

  use Grid_module

  implicit none
  
  type(grid_type), pointer, intent(inout) :: grid
  PetscReal, intent(in) :: permeability_before
  PetscReal, intent(in) :: local_stress(6), local_strain(6), local_pressure
  PetscReal, intent(in) :: Bandis_A, Bandis_B, maximum_aperture
  PetscReal, intent(in) :: normal_vector_x, normal_vector_y, normal_vector_z
  PetscReal, intent(out) :: permeability_after
  PetscReal :: effective_stress, b_p, b
  PetscInt :: icount

  effective_stress = local_stress(3)
  b_p = 1.d0  ! need to extract delta_x, delta_y, delta_x to get b_p
  b = maximum_aperture + Bandis_A*effective_stress/ &
            (1.d0 - Bandis_B*effective_stress)
  permeability_after = b**3/b_p/12.

end subroutine GeomechanicsSubsurfaceBandisPermEvaluate

! ************************************************************************** !

subroutine GeomechanicsSubsurfaceLinearPermEvaluate(permeability_before, &
                                            local_stress,local_strain, &
                                            local_pressure,permeability_after)
  !
  ! Calculates the change in porosity due to geomechanical strains
  ! following linear model (needs citation)
  !
  ! Author: Satish Karra
  ! Date: 07/29/16
  !
  
  implicit none
  PetscReal, intent(in) :: permeability_before
  PetscReal, intent(in) :: local_stress(6), local_strain(6), local_pressure
  PetscReal, intent(out) :: permeability_after
  PetscReal :: volumetric_strain
  
  volumetric_strain = local_strain(1) + local_strain(2) + local_strain(3)

  ! do nothing 
  permeability_after = permeability_before
  

end subroutine GeomechanicsSubsurfaceLinearPermEvaluate

! ************************************************************************** !

subroutine GeomechanicsSubsurfacePropsDestroy(this)
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Satish Karra
  ! Date: 07/29/16
  !

  implicit none
  
  class(geomechanics_subsurface_properties_type), pointer :: this
  
  if (.not.associated(this)) return

  deallocate(this)
  nullify(this)

end subroutine GeomechanicsSubsurfacePropsDestroy

end module Geomechanics_Subsurface_Properties_module
