module ERT_module

#include "petsc/finclude/petscksp.h"
  use petscksp

  use ERT_Aux_module
  use Material_Aux_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: ERTSetup, &
            ERTCalculateMatrix, &
            ERTCalculateAnalyticPotential, &
            ERTCalculateAverageConductivity

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

  PetscInt :: flag(1)
  PetscInt :: local_id, ghosted_id
  PetscBool :: error_found
  PetscReal :: tempreal
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%ERT => ERTAuxCreate()

  ! ensure mapping of local cell ids to neighboring ghosted ids exits
  call GridSetupCellNeighbors(grid,option)

  ! ensure that material properties specific to this module are properly
  ! initialized i.e. electrical_conductivity is initialized
  material_auxvars => patch%aux%Material%auxvars
  error_found = PETSC_FALSE
  flag = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    tempreal = minval(material_auxvars(ghosted_id)%electrical_conductivity)
    if (Uninitialized(tempreal) .and. flag(1) == 0) then
      option%io_buffer = 'ERROR: Non-initialized electrical conductivity.'
      call PrintMsgByRank(option)
      flag(1) = 1
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

subroutine ERTCalculateMatrix(realization,M,compute_delM)
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
  PetscBool :: compute_delM

  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(ert_auxvar_type), pointer :: ert_auxvars(:)
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
  PetscInt :: num_neighbors_up, num_neighbors_dn
  PetscInt :: ineighbor
  ! Matrix coefficients
  PetscReal :: coef_up, coef_dn
  PetscReal :: dcoef_up, dcoef_dn
  PetscReal :: area
  PetscReal :: factor,factor2
  ! Electrical conductivity
  PetscReal :: cond_up, cond_dn, cond_avg
  PetscReal :: dcond_avg_up, dcond_avg_dn
  PetscReal :: dist_up, dist_dn, dist_0
  PetscReal :: up_frac
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string

  option => realization%option
  field => realization%field
  patch => realization%patch
  material_auxvars => patch%aux%Material%auxvars
  ert_auxvars => patch%aux%ERT%auxvars
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
      factor = dist_up*cond_dn + dist_dn*cond_up
      cond_avg = (cond_up * cond_dn) / factor

      if (compute_delM) then
        ! For dM/dcond matrix
        ! dcond_avg_* = dcond_avg/dcond_*
        ! NB: dcond_avg_* is acutally dcond_avg_*/dist_0 
        factor2 = factor * factor
        dcond_avg_up = (dist_up*cond_dn*cond_dn) / factor2
        dcond_avg_dn = (dist_dn*cond_up*cond_up) / factor2
      endif

      area = cur_connection_set%area(iconn)

      ! get matrix coefficients for up cell
      coef_up = - cond_avg * area
      coef_dn = - coef_up                   ! cond_avg * area

      if (compute_delM) then
        ! For dM/dcond matrix wrt cond_up
        dcoef_up = - dcond_avg_up * area
        dcoef_dn = - dcoef_up                 ! dcond_avg_up * area
      endif

      if (local_id_up > 0) then
        ! set matrix coefficients for the upwind cell
        call MatSetValuesLocal(M,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                               coef_up,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesLocal(M,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                               coef_dn,ADD_VALUES,ierr);CHKERRQ(ierr)

        if (compute_delM) then
          ! For dM/dcond_up matrix
          num_neighbors_up = grid%cell_neighbors_local_ghosted(0,local_id_up)
          ineighbor = FindLocNeighbor(grid%cell_neighbors_local_ghosted      &
                                      (1:num_neighbors_up,local_id_up),      &
                                      num_neighbors_up,ghosted_id_dn)
        
          if (.not.associated(ert_auxvars(ghosted_id_up)%delM)) then
            allocate(ert_auxvars(ghosted_id_up)%delM(num_neighbors_up + 1))
            ert_auxvars(ghosted_id_up)%delM = 0.d0
          endif
        
          ! Fill values to dM/dcond_up matrix for up cell
          call FillValuesToDelM(dcoef_up,dcoef_dn,num_neighbors_up,          &
                                ineighbor,ert_auxvars(ghosted_id_up)%delM)
        endif

      endif

      if (local_id_dn > 0) then
        ! set matrix coefficients for the downwind cell
        coef_up = - coef_up
        coef_dn = - coef_dn

        ! For dM/dcond matrix
        dcoef_up =   dcond_avg_dn * area
        dcoef_dn = - dcoef_up               ! - dcond_avg_dn * area

        call MatSetValuesLocal(M,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                               coef_dn,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesLocal(M,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                               coef_up,ADD_VALUES,ierr);CHKERRQ(ierr)

        if (compute_delM) then                      
          ! For dM/dcond_dn matrix
          num_neighbors_dn = grid%cell_neighbors_local_ghosted(0,local_id_dn)
          ineighbor = FindLocNeighbor(grid%cell_neighbors_local_ghosted      &
                                      (1:num_neighbors_dn,local_id_dn),      &
                                      num_neighbors_dn,ghosted_id_up)
        
          if (.not.associated(ert_auxvars(ghosted_id_dn)%delM)) then
            allocate(ert_auxvars(ghosted_id_dn)%delM(num_neighbors_dn + 1))
            ert_auxvars(ghosted_id_dn)%delM = 0.d0
          endif

          ! Fill values to dM/dcond_dn matrix for up cell
          call FillValuesToDelM(dcoef_dn,dcoef_up,num_neighbors_dn,          &
                                ineighbor,ert_auxvars(ghosted_id_dn)%delM)
        endif

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

      ! get matrix coefficients for up cell -> NO up cell since it's
      ! the boundary
      ! down cell for it is the interior cell
      coef_dn =   cond_avg * area

      ! We need matrix coeff only for down cell so
      coef_dn = - coef_dn

      call MatSetValuesLocal(M,1,ghosted_id-1,1,ghosted_id-1,coef_dn, &
                             ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  if (realization%debug%matview_Matrix) then
    string = 'Mmatrix'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(M,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

contains
  subroutine FillValuesToDelM(dcoef_self,dcoef_neighbor,num_neighbors, &
                              ineighbor,delM)
    ! 
    ! Fills out upper traingle part of the dM/dcond matrix for each cell
    !   Storing only first rows as other rows can easily be retrieved
    !   from off-diagonal elements of the first row.
    !
    ! Author: Piyoosh Jaysaval
    ! Date: 03/08/21

    implicit none

    PetscReal :: dcoef_self, dcoef_neighbor
    PetscInt :: num_neighbors
    PetscInt :: ineighbor
    PetscReal :: delM(num_neighbors + 1)

    ! Add values for self
    delM(1) = delM(1) + dcoef_self
    ! insert values for neighbor
    delM(1+ineighbor) = dcoef_neighbor

    ! Not storing following as these can be retrieved from ineighbor's value
    !delM(1 + num_neighbor + 2*ineighbor -1) = dcoef_neighbor
    !delM(1 + num_neighbor + 2*ineighbor) = - dcoef_neighbor

  end subroutine FillValuesToDelM

  function FindLocNeighbor(neighbors,num_neighbors,id_neighbor) &
                                                            result(ineighbor)
    implicit none

    PetscInt :: num_neighbors
    PetscInt :: neighbors(num_neighbors)
    PetscInt :: id_neighbor
    PetscInt :: ineighbor

    do ineighbor = 1,num_neighbors
      if(abs(neighbors(ineighbor)) == id_neighbor) exit
    enddo

    if (ineighbor > num_neighbors) then
      option%io_buffer = 'ERTCalculateMatrix: There is something wrong ' // &
        'in finding neighbor location. '
      call PrintErrMsg(option)
    endif    

  end function FindLocNeighbor
  
end subroutine ERTCalculateMatrix

! ************************************************************************** !

subroutine ERTConductivityFromEmpiricalEqs(global_auxvar, material_auxvar)
  !
  ! Calculates conductivity using petrophysical empirical relations
  ! using Archie's law or Waxman-Smits equation
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 03/18/21
  !

  use Global_Aux_module

  implicit none

  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  ! Archie's law parameters
  ! TODO: Should read from input file
  PetscReal, parameter :: a = 1.d0        ! Tortuosity factor constant
  PetscReal, parameter :: m = 1.9d0       ! Cementation exponent
  PetscReal, parameter :: n = 2.d0        ! Saturation exponent
  PetscReal, parameter :: cond_w = 0.01d0 ! Water conductivity
  ! Waxman-Smits additional paramters
  PetscReal, parameter :: cond_c = 0.03d0 ! Clay conductivity
  PetscReal, parameter :: Vc = 0.2d0      ! Clay/Shale volume 
  ! calculated total resistivity
  PetscReal :: cond
  PetscReal :: porosity, saturation

  porosity = material_auxvar%porosity
  saturation = global_auxvar%sat(1)

  ! Archie's law
  cond = cond_w * (porosity**m) * (saturation**n) / a

  ! Waxmax-Smits equations/Dual-Water model 
  !cond = cond + cond_c * Vc * (1-porosity) * saturation**(n-1)

  material_auxvar%electrical_conductivity(1) = cond

end subroutine ERTConductivityFromEmpiricalEqs

! ************************************************************************** !

subroutine ERTCalculateAnalyticPotential(realization,ielec,average_conductivity)
  !
  ! Calculates Analytic potential for all electrodes
  ! for a given apparent/average conductivity model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/05/21
  !
  use Realization_Subsurface_class
  use Option_module
  use Grid_module
  use Patch_module
  use Survey_module

  implicit none

  type(realization_subsurface_type) :: realization
  PetscInt :: ielec
  PetscReal, optional :: average_conductivity

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(survey_type), pointer :: survey
  type(ert_auxvar_type), pointer :: ert_auxvars(:)

  PetscReal :: cond
  PetscReal :: r,epos(3),cell_center(3)

  PetscInt :: local_id
  PetscInt :: ghosted_id

  survey => realization%survey
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  ert_auxvars => patch%aux%ERT%auxvars

  if (present(average_conductivity)) then
    cond = average_conductivity
  else
    if (Initialized(survey%apparent_conductivity)) then
      cond = survey%apparent_conductivity
    else
      option%io_buffer = "ERT potential can't be computed analytically &
        &without given average conductivity or survey's &
        &apparent conductivity."
      call PrintErrMsg(option)
    endif
  endif

  epos = survey%pos_electrode(:,ielec)

  ! get & store potentials for each electrode
  do local_id=1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    cell_center(1) = grid%x(ghosted_id)
    cell_center(2) = grid%y(ghosted_id)
    cell_center(3) = grid%z(ghosted_id)

    r = norm2(epos - cell_center)
    ! Add small value to avoid overshooting at electrode position
    r = r + 1.0d-15
    ert_auxvars(ghosted_id)%potential(ielec) = 1 / (2*pi*r*cond)
  enddo

end subroutine ERTCalculateAnalyticPotential

! ************************************************************************** !

subroutine ERTCalculateAverageConductivity(realization)
  !
  ! Calculates Average conductivity of a given conductivity
  ! model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/05/21
  !
  use Realization_Subsurface_class
  use Option_module
  use Grid_module
  use Patch_module
  use Survey_module

  implicit none

  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(survey_type), pointer :: survey
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: local_average_cond
  PetscErrorCode :: ierr

  survey => realization%survey
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars

  local_average_cond = 0.d0
  ! Get part of average conductivity locally
  do local_id=1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    local_average_cond = local_average_cond + &
                         material_auxvars(ghosted_id)%electrical_conductivity(1)
  enddo
  local_average_cond = local_average_cond / grid%nmax

  ! get the average conductivity
  call MPI_Allreduce(MPI_IN_PLACE,local_average_cond,ONE_INTEGER_MPI,   &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

  survey%average_conductivity = local_average_cond

end subroutine ERTCalculateAverageConductivity

end module ERT_module