module PM_Richards_2DOFs_class

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
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_subsurface_flow_type) :: pm_richards_2DOFs_type
  contains
    procedure, public :: PreSolve => PMRichards2DOFsPreSolve
    procedure, public :: UpdateAuxVars => PMRichards2DOFsUpdateAuxVars
    procedure, public :: UpdateSolution => PMRichards2DOFsUpdateSolution
    procedure, public :: IFunction => PMRichards2DOFsIFunction
    procedure, public :: IJacobian => PMRichards2DOFsIJacobian
    procedure, public :: PostSolve => PMRichards2DOFsPostSolve
    procedure, public :: UpdateTimestep => PMRichards2DOFsUpdateTimestep
    procedure, public :: MaxChange => PMRichards2DOFsMaxChange
    procedure, public :: InitializeTimestep => PMRichards2DOFsInitializeTimestep
    procedure, public :: Destroy => PMRichards2DOFsDestroy
    procedure, public :: CheckConvergence => PMRichards2DOFsCheckConvergence
  end type pm_richards_2DOFs_type

  public :: PMRichards2DOFsCreate, &
            PMRichards2DOFsUpdateAuxVarsPatch

contains

! ************************************************************************** !

function PMRichards2DOFsCreate()
  ! 
  ! Creates Richards 2DOFs process models shell
  ! 
  ! Author: Gautam Bisht
  ! Date: 06/07/18
  ! 

  implicit none
  
  class(pm_richards_2DOFs_type), pointer :: PMRichards2DOFsCreate

  class(pm_richards_2DOFs_type), pointer :: richards_2dofs_pm
  
  allocate(richards_2dofs_pm)

  call PMSubsurfaceFlowCreate(richards_2dofs_pm)

  richards_2dofs_pm%name = 'Richards 2DOFs Flow'

  PMRichards2DOFsCreate => richards_2dofs_pm
  
end function PMRichards2DOFsCreate

! ************************************************************************** !

subroutine PMRichards2DOFsPreSolve(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 06/07/18

  implicit none
  
  class(pm_richards_2DOFs_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMRichards2DOFsPreSolve

! ************************************************************************** !

subroutine PMRichards2DOFsUpdateAuxVars(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 06/07/18

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_module
  use Discretization_module

  implicit none
  
  class(pm_richards_2DOFs_type) :: this

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: rich_auxvars(:) 
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt  :: ghosted_id, iend
  PetscReal :: compressed_porosity, dcompressed_porosity_dp
  PetscReal :: por
  PetscReal, pointer :: xx_loc_p(:)
  PetscErrorCode :: ierr

  field => this%realization%field
  option => this%realization%option
  patch => this%realization%patch
  grid => patch%grid
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  ! Given P, update auxvars
  call PMRichards2DOFsUpdateAuxVarsPatch(this%realization)

  ! Given P, compute mass = rho(P)*sat(P)*por(P)
  call VecGetArrayReadF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    if (patch%imat(ghosted_id) <= 0) cycle

    iend = ghosted_id*option%nflowdof

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvars(ghosted_id), &
             global_auxvars(ghosted_id)%pres(1), &
             compressed_porosity, dcompressed_porosity_dp)
      por = compressed_porosity
    else
      por = material_auxvars(ghosted_id)%porosity
    endif

    xx_loc_p(iend) = global_auxvars(ghosted_id)%sat(1) * &
                     global_auxvars(ghosted_id)%den(1) * por
  enddo

  call VecRestoreArrayReadF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)

  call DiscretizationLocalToGlobal(this%realization%discretization, &
         field%flow_xx_loc, field%flow_xx, NFLOWDOF)

end subroutine PMRichards2DOFsUpdateAuxVars

! ************************************************************************** !

subroutine PMRichards2DOFsUpdateAuxVarsPatch(realization)
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

  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: rich_auxvars(:) 
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  PetscInt :: ghosted_id
  PetscInt :: ibeg, iend
  PetscReal, pointer :: xx_loc_p(:),xxdot_loc_p(:)
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  grid => realization%patch%grid
  patch => realization%patch
  rich_auxvars => patch%aux%Richards%auxvars

  ! 1. Update auxvars based on new values of pressure
  call RichardsUpdateAuxVars(realization)

  ! 2. Update auxvars based on new value of dpressure_dtime, mass, and dmass_dtime
  call VecGetArrayReadF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_xxdot_loc,xxdot_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax

    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    if (patch%imat(ghosted_id) <= 0) cycle !Ignore inactive cells with inactive materials

     iend = ghosted_id*option%nflowdof
     ibeg = iend - option%nflowdof + 1

     rich_auxvars(ghosted_id)%mass = xx_loc_p(iend)
     rich_auxvars(ghosted_id)%dpres_dtime = xxdot_loc_p(ibeg)
     rich_auxvars(ghosted_id)%dmass_dtime = xxdot_loc_p(iend)

  enddo

  call VecRestoreArrayReadF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_xxdot_loc,xxdot_loc_p, ierr);CHKERRQ(ierr)

end subroutine PMRichards2DOFsUpdateAuxVarsPatch

! ************************************************************************** !

subroutine PMRichards2DOFsUpdateSolution(this)
  !
  !
  !
  ! Author: Gautam Bisht
  ! Date: 06/07/18
  !
  implicit none
  
  class(pm_richards_2DOFs_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)

  call RichardsUpdateSolution(this%realization)

end subroutine PMRichards2DOFsUpdateSolution

! ************************************************************************** !
subroutine PMRichards2DOFsIFunction(this,ts,time,U,Udot,F,ierr)
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

  class(pm_richards_2DOFs_type) :: this

  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  Vec :: F
  PetscErrorCode :: ierr

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(realization_subsurface_type), pointer :: realization
  PetscInt :: skip_conn_type

  realization => this%realization
  field => realization%field
  discretization => realization%discretization

  call VecZeroEntries(F, ierr); CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(discretization,U,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationGlobalToLocal(discretization,Udot,field%flow_xxdot_loc,NFLOWDOF)

  call PMRichards2DOFsUpdateAuxVarsPatch(realization)

  skip_conn_type = 1

  call RichardsResidualInternalConn(F,realization,skip_conn_type,ierr)
  call RichardsResidualBoundaryConn(F,realization,ierr)
  call IFunctionAccumulation(F,realization,ierr)

end subroutine PMRichards2DOFsIFunction

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

  type(realization_subsurface_type), pointer :: realization
  Vec :: F
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: ibeg, iend
  PetscReal, pointer :: f_p(:)
  PetscReal :: compressed_porosity, dcompressed_porosity_dp
  PetscReal :: por

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

    iend = ghosted_id*option%nflowdof
    ibeg = iend - option%nflowdof + 1

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvars(ghosted_id), &
             global_auxvars(ghosted_id)%pres(1), &
             compressed_porosity, dcompressed_porosity_dp)
      por = compressed_porosity
    else
      por = material_auxvars(ghosted_id)%porosity
    endif

    ! dmass_dtime * Vol
    f_p(ibeg) = f_p(ibeg) + &
                  rich_auxvars(ghosted_id)%dmass_dtime * &
                  material_auxvars(ghosted_id)%volume

    ! mass - sat(P)*rho(P)*por(P)
    f_p(iend) = f_p(iend) + &
                  rich_auxvars(ghosted_id)%mass - &
                  global_auxvars(ghosted_id)%sat(1) * &
                  global_auxvars(ghosted_id)%den(1) * por
  enddo

  call VecRestoreArrayF90(F, f_p, ierr);CHKERRQ(ierr)

end subroutine IFunctionAccumulation

! ************************************************************************** !

subroutine PMRichards2DOFsIJacobian(this,ts,time,U,Udot,shift,A,B,ierr)
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

  class(pm_richards_2DOFs_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  PetscReal :: shift
  Mat :: A, B
  PetscErrorCode :: ierr

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(realization_subsurface_type), pointer :: realization
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

end subroutine PMRichards2DOFsIJacobian

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

  type(realization_subsurface_type), pointer :: realization
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
  PetscInt :: ibeg, iend
  PetscInt :: row, col
  PetscReal:: val
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof)
  PetscReal :: compressed_porosity, dcompressed_porosity_dp
  PetscReal :: por, dpor_dP

  option => realization%option
  grid => realization%patch%grid
  patch => realization%patch
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)

    if (patch%imat(ghosted_id) <= 0) cycle

    iend = ghosted_id*option%nflowdof
    ibeg = iend - option%nflowdof + 1

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvars(ghosted_id), &
             global_auxvars(ghosted_id)%pres(1), &
             compressed_porosity, dcompressed_porosity_dp)
      por = compressed_porosity
      dpor_dP = dcompressed_porosity_dp
    else
      por = material_auxvars(ghosted_id)%porosity
    endif

    ! F1 = dmass_dtime * Vol
    ! F2 = mass - sat(P)*rho(P)*por(P)

    ! J(1,2) = shift*d(F1)/d(mdot) + d(F1)/d(m)
    row = ibeg - 1; col = iend - 1
    val = shift*material_auxvars(ghosted_id)%volume
    call MatSetValuesLocal(J,1,row,1,col,val,ADD_VALUES,ierr);CHKERRQ(ierr)

    ! J(2,1) = shift*d(F2)/d(Pdot) + d(F2)/d(P)
    row = iend - 1; col = ibeg - 1
    val = -(global_auxvars(ghosted_id)%sat(1) * rich_auxvars(ghosted_id)%dden_dp  * por     + &
            rich_auxvars(ghosted_id)%dsat_dp  * global_auxvars(ghosted_id)%den(1) * por     + &
            global_auxvars(ghosted_id)%sat(1) * global_auxvars(ghosted_id)%den(1) * dpor_dP)
    call MatSetValuesLocal(J,1,row,1,col,val,ADD_VALUES,ierr);CHKERRQ(ierr)

    ! J(2,2) = shift*d(F2)/d(mdot) + d(F2)/d(m)
    row = iend - 1; col = iend - 1
    val = 1.d0
    call MatSetValuesLocal(J,1,row,1,col,val,ADD_VALUES,ierr);CHKERRQ(ierr)

  enddo

  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

end subroutine IJacobianAccumulation

! ************************************************************************** !

subroutine PMRichards2DOFsPostSolve(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 06/07/18

  implicit none
  
  class(pm_richards_2DOFs_type) :: this

end subroutine PMRichards2DOFsPostSolve

! ************************************************************************** !
subroutine PMRichards2DOFsUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac, &
                                    time_step_max_growth_factor)
  ! 
  ! Author: Gautam Bisht
  ! Date: 06/07/18
  !
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL

  implicit none
  
  class(pm_richards_2DOFs_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: dtt
  PetscReal :: dt_p
  PetscReal :: dt_tfac
  PetscInt :: ifac

  if (iacceleration > 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = this%pressure_change_governor/(this%max_pressure_change+0.1)
      ut = up
    endif
    dtt = fac * dt * (1.d0 + ut)
  else
    ifac = max(min(num_newton_iterations,size(tfac)),1)
    dt_tfac = tfac(ifac) * dt

    fac = 0.5d0
    up = this%pressure_change_governor/(this%max_pressure_change+0.1)
    dt_p = fac * dt * (1.d0 + up)

    dtt = min(dt_tfac,dt_p)
  endif

  dtt = min(time_step_max_growth_factor*dt,dtt)
  if (dtt > dt_max) dtt = dt_max
  ! geh: There used to be code here that cut the time step if it is too
  !      large relative to the simulation time.  This has been removed.
  dtt = max(dtt,dt_min)
  dt = dtt

  call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)

end subroutine PMRichards2DOFsUpdateTimestep

! ************************************************************************** !

subroutine PMRichards2DOFsMaxChange(this)
  ! 
  ! Not needed given RichardsMaxChange is called in PostSolve
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  use Richards_module, only : RichardsMaxChange

  implicit none
  
  class(pm_richards_2DOFs_type) :: this

  call RichardsMaxChange(this%realization,this%max_pressure_change)

  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4)') this%max_pressure_change
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4)') &
      this%max_pressure_change
  endif    

end subroutine PMRichards2DOFsMaxChange

! ************************************************************************** !

subroutine PMRichards2DOFsInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  use Richards_module, only : RichardsInitializeTimestep
  
  implicit none
  
  class(pm_richards_2DOFs_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," RICHARDS FLOW ",63("="))')
  endif

  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMRichards2DOFsInitializeTimestep

! ************************************************************************** !

subroutine PMRichards2DOFsCheckConvergence(this,snes,it,xnorm,unorm,fnorm,reason,ierr)
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
  class(pm_richards_2DOFs_type) :: this
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

end subroutine PMRichards2DOFsCheckConvergence

! ************************************************************************** !

subroutine PMRichards2DOFsDestroy(this)
  ! 
  ! Destroys Richards process model
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  ! 

  use Richards_module, only : RichardsDestroy

  implicit none
  
  class(pm_richards_2DOFs_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call RichardsDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMRichards2DOFsDestroy

end module PM_Richards_2DOFs_class
