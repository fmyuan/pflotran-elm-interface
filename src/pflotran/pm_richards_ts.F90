module PM_Richards_TS_class

#include "petsc/finclude/petscts.h"
#include "petsc/finclude/petscvec.h"
  use petscts
  use petscvec
  use Richards_module
  use Richards_Aux_module
  use Richards_Common_module
  use Global_Aux_module
  use Material_Aux_class
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_Richards_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_richards_type) :: pm_richards_ts_type
  contains
    procedure, public :: UpdateAuxVars => PMRichardsTSUpdateAuxVars
    procedure, public :: IFunction => PMRichardsTSIFunction
    procedure, public :: IJacobian => PMRichardsTSIJacobian
    procedure, public :: InitializeTimestep => PMRichardsTSInitializeTimestep
    procedure, public :: Destroy => PMRichardsTSDestroy
    procedure, public :: CheckConvergence => PMRichardsTSCheckConvergence
  end type pm_richards_ts_type

  public :: PMRichardsTSCreate, &
            PMRichardsTSUpdateAuxVarsPatch

contains

! ************************************************************************** !

function PMRichardsTSCreate()
  ! 
  ! Creates Richards TS process models shell
  ! 
  ! Author: Gautam Bisht
  ! Date: 06/07/18
  ! 

  implicit none
  
  class(pm_richards_ts_type), pointer :: PMRichardsTSCreate

  class(pm_richards_ts_type), pointer :: richards_ts_pm
  
  allocate(richards_ts_pm)
  call PMRichardsInit(richards_ts_pm)
  richards_ts_pm%name = 'Richards TS Flow'
  richards_ts_pm%header = 'RICHARDS TS FLOW'

  PMRichardsTSCreate => richards_ts_pm
  
end function PMRichardsTSCreate

! ************************************************************************** !

subroutine PMRichardsTSUpdateAuxVars(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 06/07/18

  implicit none
  
  class(pm_richards_ts_type) :: this

  call PMRichardsTSUpdateAuxVarsPatch(this%realization)

end subroutine PMRichardsTSUpdateAuxVars

! ************************************************************************** !

subroutine PMRichardsTSUpdateAuxVarsPatch(realization)
  ! 
  ! Author: Gautam Bisht
  ! Date: 06/07/18

  use Realization_Subsurface_class
  use Field_module
  use Grid_module
  use Patch_module
  use Option_module
  use Richards_Aux_module, only : RichardsAuxVarCompute2ndOrderDeriv

  implicit none
  
  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(richards_auxvar_type), pointer :: rich_auxvars(:) 
  PetscInt :: ghosted_id
  PetscReal, pointer :: xx_loc_p(:),xxdot_loc_p(:)
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  rich_auxvars => patch%aux%Richards%auxvars

  ! 1. Update auxvars based on new values of pressure
  call RichardsUpdateAuxVars(realization)

  ! 2. Update auxvars based on new value of dpressure_dtime, mass, and 
  !    dmass_dtime
  call VecGetArrayReadF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_xxdot_loc,xxdot_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle 
    rich_auxvars(ghosted_id)%dpres_dtime = xxdot_loc_p(ghosted_id)
    call RichardsAuxVarCompute2ndOrderDeriv(rich_auxvars(ghosted_id), &
                                       patch%characteristic_curves_array( &
                                         patch%sat_func_id(ghosted_id))%ptr, &
                                       option)
  enddo

  call VecRestoreArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_xxdot_loc,xxdot_loc_p, &
                              ierr);CHKERRQ(ierr)

end subroutine PMRichardsTSUpdateAuxVarsPatch

! ************************************************************************** !

subroutine PMRichardsTSIFunction(this,ts,time,U,Udot,F,ierr)
  !
  !
  !
  ! Author: Gautam Bisht
  ! Date: 06/07/18
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module

  implicit none

  class(pm_richards_ts_type) :: this

  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  Vec :: F
  PetscErrorCode :: ierr

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  class(realization_subsurface_type), pointer :: realization
  PetscInt :: skip_conn_type

  realization => this%realization
  field => realization%field
  discretization => realization%discretization

  call VecZeroEntries(F, ierr); CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(discretization,U,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationGlobalToLocal(discretization,Udot,field%flow_xxdot_loc, &
                                   NFLOWDOF)

  call PMRichardsTSUpdateAuxVarsPatch(realization)

  skip_conn_type = 1

  call RichardsResidualInternalConn(F,realization,skip_conn_type,ierr)
  call RichardsResidualBoundaryConn(F,realization,ierr)
  call RichardsResidualSourceSink(F,realization,ierr)
  call IFunctionAccumulation(F,realization,ierr)

end subroutine PMRichardsTSIFunction
! ************************************************************************** !
subroutine IFunctionAccumulation(F,realization,ierr)
  !
  !
  !
  ! Author: Gautam Bisht
  ! Date: 06/07/18
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Patch_module
  use Option_module
  use Material_Aux_class

  implicit none

  class(realization_subsurface_type), pointer :: realization
  Vec :: F
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: ibeg
  PetscReal, pointer :: f_p(:)
  PetscReal :: compressed_porosity, dcompressed_porosity_dp
  PetscReal :: por, dpor_dP, dmass_dP

  option => realization%option
  grid => realization%patch%grid
  patch => realization%patch
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(F, f_p, ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)

    if (patch%imat(ghosted_id) <= 0) cycle

    ibeg = (ghosted_id-1)*option%nflowdof + 1

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvars(ghosted_id), &
             global_auxvars(ghosted_id)%pres(1), &
             compressed_porosity, dcompressed_porosity_dp)
      por = compressed_porosity
      dpor_dP = dcompressed_porosity_dp
    else
      por = material_auxvars(ghosted_id)%porosity
      dpor_dP = 0.d0
    endif

    ! F = d(rho*phi*s)/dP * dP_dtime * Vol
    dmass_dP = global_auxvars(ghosted_id)%den(1)*dpor_dP* &
                             global_auxvars(ghosted_id)%sat(1) + &
               rich_auxvars(ghosted_id)%dden_dp *por    * &
                             global_auxvars(ghosted_id)%sat(1) + &
               global_auxvars(ghosted_id)%den(1)*por    * &
                             rich_auxvars(ghosted_id)%dsat_dp

    f_p(ibeg) = f_p(ibeg) + &
                  dmass_dP*rich_auxvars(ghosted_id)%dpres_dtime * &
                  material_auxvars(ghosted_id)%volume

  enddo

  call VecRestoreArrayF90(F, f_p, ierr);CHKERRQ(ierr)

end subroutine IFunctionAccumulation

! ************************************************************************** !

subroutine PMRichardsTSIJacobian(this,ts,time,U,Udot,shift,A,B,ierr)
  !
  !
  !
  ! Author: Gautam Bisht
  ! Date: 06/07/18
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module

  implicit none

  class(pm_richards_ts_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  PetscReal :: shift
  Mat :: A, B
  PetscErrorCode :: ierr

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  class(realization_subsurface_type), pointer :: realization
  Mat :: J

  realization => this%realization
  field => realization%field
  discretization => realization%discretization

  J = B

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  call RichardsJacobianInternalConn(J,realization,ierr)
  call RichardsJacobianBoundaryConn(J,realization,ierr)
  call RichardsJacobianSourceSink(J,realization,ierr)
  call IJacobianAccumulation(J,shift,realization,ierr)

  if (A /= B) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr);
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr);
  endif

end subroutine PMRichardsTSIJacobian

! ************************************************************************** !
subroutine IJacobianAccumulation(J,shift,realization,ierr)
  !
  !
  !
  ! Author: Gautam Bisht
  ! Date: 06/07/18
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Patch_module
  use Option_module
  use Material_Aux_class

  implicit none

  class(realization_subsurface_type), pointer :: realization
  PetscReal :: shift
  Mat :: J
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: ibeg
  PetscInt :: row, col
  PetscReal:: val
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof)
  PetscReal :: compressed_porosity, dcompressed_porosity_dp
  PetscReal :: por, dpor_dP, d2por_dP2
  PetscReal :: sat, dsat_dP, d2sat_dP2
  PetscReal :: den, dden_dP, d2den_dP2
  PetscReal :: dmass_dP, d2mass_dP2

  option => realization%option
  grid => realization%patch%grid
  patch => realization%patch
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)

    if (patch%imat(ghosted_id) <= 0) cycle

    ibeg = (ghosted_id-1)*option%nflowdof + 1

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvars(ghosted_id), &
             global_auxvars(ghosted_id)%pres(1), &
             compressed_porosity, dcompressed_porosity_dp)
      por = compressed_porosity
      dpor_dP = dcompressed_porosity_dp
    else
      por = material_auxvars(ghosted_id)%porosity
    endif

    ! F = d(rho*phi*s)/dP * dP_dtime * Vol

    ! J(1,1) = shift*d(F)/d(Pdot) + d(F)/d(P)
    !        = shift*d(rho*phi*s)/dP        * Vol +
    !          d2(rho*phi*s)/dP2 * dP_dtime * Vol

    den = global_auxvars(ghosted_id)%den(1)
    sat = global_auxvars(ghosted_id)%sat(1)
    dden_dP = rich_auxvars(ghosted_id)%dden_dp
    dsat_dP = rich_auxvars(ghosted_id)%dsat_dp
    d2den_dP2 = rich_auxvars(ghosted_id)%d2den_dp2
    d2sat_dP2 = rich_auxvars(ghosted_id)%d2sat_dp2
    d2por_dP2 = 0.d0

    row = ibeg - 1; col = ibeg - 1
    dmass_dP = ( &
      sat     * dden_dP * por     + &
      dsat_dP * den     * por     + &
      sat     * den     * dpor_dP &
      )

    d2mass_dP2 = ( &
      dsat_dP   * dden_dP   * por       + &
      sat       * d2den_dP2 * por       + &
      sat       * dden_dP   * dpor_dP   + &
      d2sat_dP2 * den       * por       + &
      dsat_dP   * dden_dP   * por       + &
      dsat_dP   * den       * dpor_dP   + &
      dsat_dP   * den       * dpor_dP   + &
      sat       * dden_dP   * dpor_dP   + &
      sat       * den       * d2por_dP2 &
       )

    val = (shift* dmass_dP + d2mass_dP2 * &
            rich_auxvars(ghosted_id)%dpres_dtime)* &
            material_auxvars(ghosted_id)%volume

    call MatSetValuesLocal(J,1,row,1,col,val,ADD_VALUES,ierr);CHKERRQ(ierr)

  enddo

  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

end subroutine IJacobianAccumulation

! ************************************************************************** !

subroutine PMRichardsTSInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  use Richards_module, only : RichardsInitializeTimestep
  
  implicit none
  
  class(pm_richards_ts_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," RICHARDS TS FLOW ",60("="))')
  endif

  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMRichardsTSInitializeTimestep

! ************************************************************************** !

subroutine PMRichardsTSCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)
  ! 
  ! Adds a convergence check for the nonlinear problem
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  ! 

#include "petsc/finclude/petscsnes.h"
  use petscsnes

  implicit none
  
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm ! 2-norm of updated solution
  PetscReal :: unorm ! 2-norm of update. PETSc refers to this as snorm
  PetscReal :: fnorm ! 2-norm of updated residual
  SNESConvergedReason :: reason
  class(pm_richards_ts_type) :: this
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string

  call SNESConvergedDefault(snes,it,xnorm,unorm,fnorm,reason, &
                            0,ierr);CHKERRQ(ierr)

  if (this%option%print_screen_flag) then
    select case(int(reason))
      case(2)
        string = 'atol'
      case(3)
        string = 'rtol'
      case(4)
        string = 'stol'
      case default
        write(string,'(i3)') reason
    end select

    write(*,'(i3," 2r:",es9.2, &
            & " 2x:",es9.2, &
            & " 2u:",es9.2, &
            & " rsn: ",a)') &
            it, fnorm, xnorm, unorm, trim(string)
  endif

end subroutine PMRichardsTSCheckConvergence

! ************************************************************************** !

subroutine PMRichardsTSDestroy(this)
  ! 
  ! Destroys Richards process model
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  ! 
  use Richards_module, only : RichardsDestroy

  implicit none
  
  class(pm_richards_ts_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call PMRichardsDestroy(this)
  
end subroutine PMRichardsTSDestroy

end module PM_Richards_TS_class
