module ERT_module

#include "petsc/finclude/petscksp.h"
  use petscksp

  use ERT_Aux_module
  use Material_Aux_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: ERTSetup, &
            ERTCalculateMatrix

contains

! ************************************************************************** !

subroutine ERTSetup(realization)
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/29/21
  ! 

  use Realization_Subsurface_class
  use Patch_module

  type(realization_subsurface_type) :: realization
  
  call ERTSetupPatch(realization)
 
  ! Setup other ERT requirements e.g. plot, output, etc.

end subroutine ERTSetup

! ************************************************************************** !

subroutine ERTSetupPatch(realization)
  ! 
  ! Creates arrays for ERT auxiliary variables
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/29/21
  ! 

  use Realization_Subsurface_class
  use Option_module 
  use Patch_module
  use Grid_module
  use Survey_module   

  implicit none

  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(survey_type), pointer :: survey

  class(material_auxvar_type), pointer :: material_auxvars(:) 
  type(ert_auxvar_type), pointer :: ert_auxvars(:) 

  PetscInt :: flag(2)
  PetscInt :: local_id, ghosted_id
  PetscBool :: error_found
  PetscReal :: tempreal
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%ERT => ERTAuxCreate()

  ! ensure that material properties specific to this module are properly
  ! initialized i.e. electrical_conductivity is initialized
  material_auxvars => patch%aux%Material%auxvars
  error_found = PETSC_FALSE
  flag = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle       
    tempreal = minval(material_auxvars(ghosted_id)%electrical_conductivity)   
    if (Uninitialized(tempreal) .and. flag(2) == 0) then
      option%io_buffer = 'ERROR: Non-initialized electrical conductivity.'
      call PrintMsgByRank(option)
      flag(2) = 1
    endif      
  enddo  

  error_found = error_found .or. (maxval(flag) > 0)
  call MPI_Allreduce(MPI_IN_PLACE,error_found,ONE_INTEGER_MPI,MPI_LOGICAL, &
                     MPI_LOR,option%mycomm,ierr)
  if (error_found) then
    option%io_buffer = 'Material property errors found in ERTSetup.'
    call PrintErrMsg(option)
  endif

  survey => realization%survey

 ! allocate auxvars data structures for all grid cells  
  allocate(ert_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call ERTAuxVarInit(ert_auxvars(ghosted_id),survey,option)   
  enddo
  patch%aux%ERT%auxvars => ert_auxvars
  patch%aux%ERT%num_aux = grid%ngmax

end subroutine ERTSetupPatch

! ************************************************************************** !

subroutine ERTCalculateMatrix(realization,M)
  ! 
  ! Calculate System matrix for ERT
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/26/21
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Grid_module
  use Patch_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Debug_module

  implicit none

  type(realization_subsurface_type) :: realization
  Mat :: M

  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid 
  type(field_type), pointer :: field 
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, ghosted_id_up
  PetscInt :: local_id_dn, ghosted_id_dn

  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  ! Matrix coefficients
  PetscReal :: coef_up, coef_dn
  PetscReal :: area
  ! Electrical conductivity
  PetscReal :: cond_up, cond_dn, cond_avg
  PetscReal :: dist_up, dist_dn, dist_0
  PetscReal :: up_frac 
  PetscViewer :: viewer
  PetscErrorCode :: ierr 
  character(len=MAXSTRINGLENGTH) :: string

  option => realization%option
  field => realization%field
  patch => realization%patch
  material_auxvars => patch%aux%Material%auxvars
  grid => patch%grid

  ! Pre-set Matrix to zeros
  call MatZeroEntries(M,ierr); CHKERRQ(ierr)

  ! Setting matrix enteries for Internal Flux terms/connections
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do
    if (.not. associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ! get ghosted ids of up and down
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
       
      ! Ghosted to local id mapping. Local id is zero/-1 for
      ! ghosted cells
      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      
      ! cycle if material is negative for any cell -> inactive cell
      if (patch%imat(ghosted_id_up) <= 0 .or.   & 
          patch%imat(ghosted_id_dn) <=0 ) cycle

      cond_up = material_auxvars(ghosted_id_up)%electrical_conductivity(1)
      cond_dn = material_auxvars(ghosted_id_dn)%electrical_conductivity(1)
      
      !dist(-1) -> scalar fractional distance up
      !dist(-1) = d_up/d_0
      up_frac = cur_connection_set%dist(-1,iconn)
      dist_0  = cur_connection_set%dist( 0,iconn)
      dist_up = up_frac * dist_0
      dist_dn = dist_0 - dist_up
      
      ! get harmonic averaged conductivity at the face
      ! NB: cond_avg is actually cond_avg/dist_0 
      cond_avg = (cond_up * cond_dn) / (dist_up*cond_dn + dist_dn*cond_up)

      area = cur_connection_set%area(iconn)

      ! get matrix coefficients for up cell
      coef_up = - cond_avg * area
      coef_dn =   cond_avg * area

      if (local_id_up > 0) then
        ! set matrix coefficients for the upwind cell
        call MatSetValuesLocal(M,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                               coef_up,ADD_VALUES,ierr); CHKERRQ(ierr)
        call MatSetValuesLocal(M,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                               coef_dn,ADD_VALUES,ierr); CHKERRQ(ierr)                       
      endif

      if (local_id_dn > 0) then
        ! set matrix coefficients for the downwind cell
        coef_up = - coef_up
        coef_dn = - coef_dn
        call MatSetValuesLocal(M,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                               coef_dn,ADD_VALUES,ierr); CHKERRQ(ierr)
        call MatSetValuesLocal(M,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                               coef_up,ADD_VALUES,ierr);CHKERRQ(ierr)        
      endif

    enddo
    cur_connection_set => cur_connection_set%next      
  enddo
  
  ! Add Dirichley Boundary condition -> potential at boundaries = 0
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit

    cur_connection_set => boundary_condition%connection_set

    do iconn = 1, cur_connection_set%num_connections

      sum_connection = sum_connection + 1
  
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      cond_dn = material_auxvars(ghosted_id)%electrical_conductivity(1)

      dist_0  = cur_connection_set%dist( 0,iconn)

      ! get harmonic averaged conductivity at the face
      ! NB: cond_avg is actually cond_avg/dist_0 
      ! Use just the same value of down cell
      ! also note dist(0) = distance from center to boundary/face
      cond_avg = cond_dn / dist_0

      area = cur_connection_set%area(iconn)

      ! get matrix coefficients for up cell -> NO up cell since it's the boundary
      ! down cell for it is the interior cell
      coef_dn =   cond_avg * area

      ! We need matrix coeff only for down cell so
      coef_dn = - coef_dn

      call MatSetValuesLocal(M,1,ghosted_id-1,1,ghosted_id-1,coef_dn, &
                             ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
  call MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)

  if (realization%debug%matview_Matrix) then
    string = 'Mmatrix'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(M,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

end subroutine ERTCalculateMatrix

end module ERT_module