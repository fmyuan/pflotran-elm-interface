module PM_TH_TS_class

#include "petsc/finclude/petscts.h"
#include "petsc/finclude/petscvec.h"
  use petscts
  use petscvec
  use TH_module
  use TH_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_TH_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_th_type) :: pm_th_ts_type
  contains
    procedure, public :: UpdateAuxVars => PMTHTSUpdateAuxVars
    procedure, public :: IFunction => PMTHTSIFunction
    procedure, public :: IJacobian => PMTHTSIJacobian
    procedure, public :: InitializeTimestep => PMTHTSInitializeTimestep
    procedure, public :: Destroy => PMTHTSDestroy
    procedure, public :: CheckConvergence => PMTHTSCheckConvergence
  end type pm_th_ts_type

  public :: PMTHTSCreate, &
            PMTHTSUpdateAuxVarsPatch

  
contains

! ************************************************************************** !

function PMTHTSCreate()
  ! 
  ! Creates TH TS process models shell
  ! 
  ! Author: Satish Karra
  ! Date: 05/08/19
  ! 

  implicit none
  
  class(pm_th_ts_type), pointer :: PMTHTSCreate

  class(pm_th_ts_type), pointer :: this

  PetscReal, parameter :: pres_abs_inf_tol = 1.d0
  PetscReal, parameter :: temp_abs_inf_tol = 1.d-5
  PetscReal, parameter :: abs_update_inf_tol(2) = &
                            [pres_abs_inf_tol,temp_abs_inf_tol]
  PetscReal, parameter :: pres_rel_inf_tol = 1.d-5
  PetscReal, parameter :: temp_rel_inf_tol = 1.d-5
  PetscReal, parameter :: rel_update_inf_tol(2) = &
                            [pres_rel_inf_tol,temp_rel_inf_tol]
  PetscReal, parameter :: residual_abs_inf_tol(2) = 1.d-5
  PetscReal, parameter :: residual_scaled_inf_tol(2) = 1.d-5
  
#ifdef PM_TH_TS_DEBUG
  print *, 'PMTHTSCreate()'
#endif  

  allocate(this)

  nullify(this%commN)

  call PMSubsurfaceFlowCreate(this)
  this%name = 'TH_TS Flow'
  this%header = 'TH_TS FLOW'

  this%residual_abs_inf_tol = residual_abs_inf_tol
  this%residual_scaled_inf_tol = residual_scaled_inf_tol
  this%abs_update_inf_tol = abs_update_inf_tol
  this%rel_update_inf_tol = rel_update_inf_tol

  PMTHTSCreate => this
  
end function PMTHTSCreate

! ************************************************************************** !

subroutine PMTHTSUpdateAuxVars(this)
  ! 
  ! Author: Satish Karra
  ! Date: 05/08/19

  implicit none
  
  class(pm_th_ts_type) :: this

  call PMTHTSUpdateAuxVarsPatch(this%realization)

end subroutine PMTHTSUpdateAuxVars


! ************************************************************************** !

subroutine PMTHTSUpdateAuxVarsPatch(realization)
  ! 
  ! Author: Satish Karra
  ! Date: 05/08/19


  use Realization_Subsurface_class
  use Field_module
  use Grid_module
  use Patch_module
  use Option_module
  
  implicit none

  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(TH_auxvar_type), pointer :: TH_auxvars(:) 
  type(TH_parameter_type), pointer :: TH_parameter
  PetscInt :: ghosted_id,local_id,istart,iend
  PetscReal, pointer :: xx_loc_p(:),xxdot_loc_p(:)
  PetscReal, pointer :: icap_loc_p(:),ithrm_loc_p(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: ithrm,icap
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  TH_auxvars => patch%aux%TH%auxvars
  TH_parameter => patch%aux%TH%TH_parameter
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  ! 1. Update auxvars based on new values of pressure, temperature
  call THUpdateAuxVars(realization)

  ! 2. Update auxvars based on new value of dpressure_dtime, mass, and 
  !    dmass_dtime
  call VecGetArrayReadF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_xxdot_loc,xxdot_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  
  
  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle 
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    icap = int(icap_loc_p(ghosted_id))
    ithrm = int(ithrm_loc_p(ghosted_id))

    th_auxvars(ghosted_id)%dpres_dtime = & 
      xxdot_loc_p((ghosted_id-1)*option%nflowdof+1)
    th_auxvars(ghosted_id)%dtemp_dtime = &
      xxdot_loc_p((ghosted_id-1)*option%nflowdof+2)
    call THAuxVarCompute2ndOrderDeriv(TH_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id), &
                                  material_auxvars(ghosted_id), &
                                  TH_parameter,ithrm, &
                                  patch%characteristic_curves_array(icap)%ptr, &
                                  option)
  enddo

  call VecRestoreArrayReadF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_xxdot_loc,xxdot_loc_p, &
                              ierr);CHKERRQ(ierr)
  

end subroutine PMTHTSUpdateAuxVarsPatch

! ************************************************************************** !

subroutine PMTHTSIFunction(this,ts,time,U,Udot,F,ierr)
  !
  !
  ! Author: Satish Karra
  ! Date: 05/08/19
  
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module

  implicit none
  
  class(pm_th_ts_type) :: this

  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  Vec :: F
  PetscErrorCode :: ierr

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  class(realization_subsurface_type), pointer :: realization

  realization => this%realization
  field => realization%field
  discretization => realization%discretization

  call VecZeroEntries(F, ierr); CHKERRQ(ierr)

  call THUpdateLocalVecs(U,realization,ierr)

  call DiscretizationGlobalToLocal(discretization,Udot,field%flow_xxdot_loc, &
                                   NFLOWDOF)

  call PMTHTSUpdateAuxVarsPatch(realization)


  call THResidualInternalConn(F,realization,ierr)
  call THResidualBoundaryConn(F,realization,ierr)
  call THResidualSourceSink(F,realization,ierr)
  call IFunctionAccumulation(F,realization,ierr)


end subroutine PMTHTSIFunction


! ************************************************************************** !
subroutine IFunctionAccumulation(F,realization,ierr)
  !
  !
  !
  ! Author: Satish Karra
  ! Date: 06/07/19
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Patch_module
  use Option_module
  use Material_Aux_class
  use Field_module

  implicit none
  
  class(realization_subsurface_type), pointer :: realization
  Vec :: F
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(TH_auxvar_type), pointer :: TH_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(TH_parameter_type), pointer :: TH_parameter
  PetscInt :: local_id, ghosted_id

  PetscInt :: istart, iend
  PetscReal, pointer :: f_p(:)
  PetscReal :: compressed_porosity, dcompressed_porosity_dp
  PetscReal :: por, den, sat, u, temp
  PetscReal :: dpor_dP, dpor_dt
  PetscReal :: dmass_dP, dmass_dt
  PetscReal :: denergy_dP, denergy_dt
  PetscReal :: dden_dP, dden_dt
  PetscReal :: dsat_dP, dsat_dt
  PetscReal  :: du_dP, du_dt
  PetscReal :: rock_dencpr
  PetscReal, pointer :: ithrm_loc_p(:)
  
  option => realization%option
  grid => realization%patch%grid
  patch => realization%patch
  field => realization%field
  TH_auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  TH_parameter => patch%aux%TH%TH_parameter

  call VecGetArrayF90(F, f_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1

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
    
    den = global_auxvars(ghosted_id)%den(1)
    sat = global_auxvars(ghosted_id)%sat(1)
    temp = global_auxvars(ghosted_id)%temp
    u = TH_auxvars(ghosted_id)%u
    dden_dP = TH_auxvars(ghosted_id)%dden_dp
    dsat_dP = TH_auxvars(ghosted_id)%dsat_dp
    dpor_dt = 0.d0
    dden_dt = TH_auxvars(ghosted_id)%dden_dt
    dsat_dt = TH_auxvars(ghosted_id)%dsat_dt
    du_dP = TH_auxvars(ghosted_id)%du_dp
    du_dt = TH_auxvars(ghosted_id)%du_dt
    rock_dencpr = TH_parameter%dencpr(int(ithrm_loc_p(ghosted_id)))

    ! A_M = d(rho*phi*s)/dP * dP_dtime *Vol + d(rho*phi*s)/dT * dT_dtime *Vol
    dmass_dP = den     * dpor_dP * sat     + &
               dden_dP * por     * sat     + &
               den     * por     * dsat_dP

    dmass_dt = den     * dpor_dt * sat + &
               dden_dt * por     * sat + &
               den     * por     * dsat_dt

    f_p(istart) = f_p(istart) + &
                    dmass_dP*TH_auxvars(ghosted_id)%dpres_dtime * &
                    material_auxvars(ghosted_id)%volume + &
                    dmass_dt*TH_auxvars(ghosted_id)%dtemp_dtime * &
                    material_auxvars(ghosted_id)%volume
    
    
    ! A_E = [d(rho*phi*s*U)/dP + d(rho*(1-phi)*T)/dP] * dP_dtime *Vol + 
    !       [d(rho*phi*s*U)/dT + d(rho*(1-phi)*T)/dT] * dT_dtime *Vol
    denergy_dP = dden_dP     * por        * sat     * u     + &
                 den         * dpor_dP    * sat     * u     + &
                 den         * por        * dsat_dP * u     + &
                 den         * por        * sat     * du_dP + &
                 rock_dencpr * (-dpor_dP) * temp

    denergy_dt = dden_dt     * por        * sat     * u     + &
                 den         * dpor_dt    * sat     * u     + &
                 den         * por        * dsat_dt * u     + &
                 den         * por        * sat     * du_dt + &
                 rock_dencpr * (-dpor_dt) * temp            + &
                 rock_dencpr * (1-por)
                 
                                  
    f_p(iend) = f_p(iend) + &
                  denergy_dP*TH_auxvars(ghosted_id)%dpres_dtime * &
                  material_auxvars(ghosted_id)%volume + &
                  denergy_dt*TH_auxvars(ghosted_id)%dtemp_dtime * &
                  material_auxvars(ghosted_id)%volume
    
    
  enddo
    
  call VecRestoreArrayF90(F, f_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
    
    
end subroutine IFunctionAccumulation

! ************************************************************************** !

subroutine PMTHTSIJacobian(this,ts,time,U,Udot,shift,A,B,ierr)
  !
  !
  !
  ! Author: Satish Karra
  ! Date: 06/07/19
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Debug_module
  use Option_module
  
  implicit none
  
  class(pm_th_ts_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  PetscReal :: shift
  Mat :: A, B
  PetscErrorCode :: ierr  
  PetscViewer :: viewer

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  class(realization_subsurface_type), pointer :: realization
  Mat :: J

  realization => this%realization
  field => realization%field
  discretization => realization%discretization
  option => realization%option

  J = B

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  call THJacobianInternalConn(J,realization,ierr)
  call THJacobianBoundaryConn(J,realization,ierr)
  call THJacobianSourceSink(J,realization,ierr)
  call IJacobianAccumulation(J,shift,realization,ierr)

  if (A /= B) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr);
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr);
  endif  
  
  
  call PetscViewerASCIIOpen(option%mycomm,'THTSjacobian.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(J,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

 
  
end subroutine PMTHTSIJacobian


! ************************************************************************** !
subroutine IJacobianAccumulation(J,shift,realization,ierr)
  !
  !
  !
  ! Author: Satish Karra
  ! Date: 06/07/19
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Patch_module
  use Option_module
  use Material_Aux_class
  use Field_module

  implicit none
  
  class(realization_subsurface_type), pointer :: realization
  PetscReal :: shift
  Mat :: J
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(TH_auxvar_type), pointer :: TH_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(TH_parameter_type), pointer :: TH_parameter

  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend
  PetscReal :: compressed_porosity, dcompressed_porosity_dp
  PetscReal :: por, dpor_dP, d2por_dP2, dpor_dt, d2por_dt2
  PetscReal :: d2por_dtdP, d2por_dPdt
  PetscReal :: sat, dsat_dP, d2sat_dP2, dsat_dt, d2sat_dt2
  PetscReal :: d2sat_dtdP, d2sat_dPdt
  PetscReal :: den, dden_dP, d2den_dP2, dden_dt, d2den_dt2
  PetscReal :: d2den_dtdP, d2den_dPdt
  PetscReal :: u, du_dP, du_dt, d2u_dP2, d2u_dt2, d2u_dtdP, d2u_dPdt
  PetscReal :: dmass_dP, d2mass_dP2, d2mass_dPdt
  PetscReal :: dmass_dt, d2mass_dt2, d2mass_dtdP
  PetscReal :: denergy_dP, d2energy_dP2, d2energy_dPdt
  PetscReal :: denergy_dt, d2energy_dt2, d2energy_dtdP
  PetscReal :: Jlocal(2,2)
  PetscReal :: temp
  PetscReal :: rock_dencpr
  PetscReal, pointer :: ithrm_loc_p(:)

  option => realization%option
  grid => realization%patch%grid
  patch => realization%patch
  field => realization%field
  TH_auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars  
  TH_parameter => patch%aux%TH%TH_parameter


  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1

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
    
    den = global_auxvars(ghosted_id)%den(1)
    sat = global_auxvars(ghosted_id)%sat(1)
    temp = global_auxvars(ghosted_id)%temp
    u = TH_auxvars(ghosted_id)%u
    rock_dencpr = TH_parameter%dencpr(int(ithrm_loc_p(ghosted_id)))
    dden_dP = TH_auxvars(ghosted_id)%dden_dp
    dsat_dP = TH_auxvars(ghosted_id)%dsat_dp
    dpor_dt = 0.d0
    dden_dt = TH_auxvars(ghosted_id)%dden_dt
    dsat_dt = TH_auxvars(ghosted_id)%dsat_dt
    du_dP = TH_auxvars(ghosted_id)%du_dp
    du_dt = TH_auxvars(ghosted_id)%du_dt
    d2den_dP2 = TH_auxvars(ghosted_id)%d2den_dp2
    d2den_dtdP = TH_auxvars(ghosted_id)%d2den_dtdp
    d2den_dt2 = TH_auxvars(ghosted_id)%d2den_dt2
    d2sat_dP2 = TH_auxvars(ghosted_id)%d2sat_dp2
    d2sat_dtdP = TH_auxvars(ghosted_id)%d2sat_dtdp
    d2sat_dt2 = TH_auxvars(ghosted_id)%d2sat_dt2
    d2u_dP2 = TH_auxvars(ghosted_id)%d2u_dp2
    d2u_dtdP = TH_auxvars(ghosted_id)%d2u_dtdp
    d2u_dt2 = TH_auxvars(ghosted_id)%d2u_dt2
    d2por_dP2 = 0.d0
    d2por_dtdP = 0.d0
    d2por_dt2 = 0.d0
    d2den_dPdt = d2den_dtdP
    d2sat_dPdt = d2sat_dtdP
    d2por_dPdt = d2por_dtdP
    d2u_dPdt = d2u_dtdP
    
    ! A_M = d(rho*phi*s)/dP * dP_dtime * Vol + d(rho*phi*s)/dT * dT_dtime * Vol
    
    ! Jlocal(1,1) = shift*d(A_M)/d(Pdot) + d(A_M)/d(P)
    !         		= shift*d(rho*phi*s)/dP*Vol + d2(rho*phi*s)/dP2*dP_dtime*Vol +
    !               d2(rho*phi*s)/dTdP*dT_dtime*Vol
    
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
    
    d2mass_dPdt = ( &
      dsat_dt    * dden_dP    * por       + &
      sat        * d2den_dtdP * por       + &
      sat        * dden_dP    * dpor_dt   + &
      d2sat_dtdP * den        * por       + &
      dsat_dP    * dden_dt    * por       + &
      dsat_dP    * den        * dpor_dt   + &
      dsat_dt    * den        * dpor_dP   + &
      sat        * dden_dt    * dpor_dP   + &
      sat        * den        * d2por_dtdP &
      )
      
    d2mass_dtdP = d2mass_dPdt
     
    Jlocal(1,1) = (shift*dmass_dP + &
                   d2mass_dP2*TH_auxvars(ghosted_id)%dpres_dtime + &
                   d2mass_dtdP*TH_auxvars(ghosted_id)%dtemp_dtime)* &
                   material_auxvars(ghosted_id)%volume

    
    ! Jlocal(1,2) = shift*d(A_M)/d(Tdot) + d(A_M)/d(T)
    !             = shift*d(rho*phi*s)/dT*Vol + d2(rho*phi*s)/dT2*dT_dtime*Vol +
    !               d2(rho*phi*s)/dPdT*dP_dtime*Vol  
    
    dmass_dt = ( &
      sat     * dden_dt * por     + &
      dsat_dt * den     * por     + &
      sat     * den     * dpor_dt &
      )   
    
    d2mass_dt2 = ( &
      dsat_dt   * dden_dt   * por       + &
      sat       * d2den_dt2 * por       + &
      sat       * dden_dt   * dpor_dt   + &
      d2sat_dt2 * den       * por       + &
      dsat_dt   * dden_dt   * por       + &
      dsat_dt   * den       * dpor_dt   + &
      dsat_dt   * den       * dpor_dt   + &
      sat       * dden_dt   * dpor_dt   + &
      sat       * den       * d2por_dt2 &
      )    
      
    Jlocal(1,2) = (shift*dmass_dt + &
                   d2mass_dt2*TH_auxvars(ghosted_id)%dtemp_dtime + &
                   d2mass_dPdt*TH_auxvars(ghosted_id)%dpres_dtime)* &
                   material_auxvars(ghosted_id)%volume
    
    
    ! A_E = [d(rho*phi*s*U)/dP + d(rock_dencpr*(1-phi)*T)/dP] * dP_dtime *Vol + 
    !       [d(rho*phi*s*U)/dT + d(rock_dencpr*(1-phi)*T)/dT] * dT_dtime *Vol
    
    
    ! Jlocal(2,1) = shift*d(A_E)/d(Pdot) + d(A_E)/d(P)
    !             = shift*[d(rho*phi*s*U)/dP + d(rock_dencpr*(1-phi)*T)/dP]*Vol + 
    !               [d2(rho*phi*s*U)/dP2 + d2(rock_dencpr*(1-phi)*T)/dP2]*dP_dtime*Vol +
    !               [d2(rho*phi*s*U)/dTdP + d2(rock_dencpr*(1-phi)*T)/dTdP]*dT_dtime*Vol
    
    denergy_dP = dden_dP     * por        * sat     * u     + &
                 den         * dpor_dP    * sat     * u     + &
                 den         * por        * dsat_dP * u     + &
                 den         * por        * sat     * du_dp + &
                 rock_dencpr * (-dpor_dP) * temp 
    
    d2energy_dP2 = d2den_dP2   * por          * sat       * u       + &
                   dden_dP     * dpor_dP      * sat       * u       + &
                   dden_dP     * por          * dsat_dP   * u       + &
                   dden_dP     * por          * sat       * du_dP   + &
                   dden_dP     * dpor_dP      * sat       * u       + &
                   den         * d2por_dP2    * sat       * u       + &
                   den         * dpor_dP      * dsat_dP   * u       + &
                   den         * dpor_dP      * sat       * du_dP   + &
                   dden_dP     * por          * dsat_dP   * u       + &
                   den         * dpor_dP      * dsat_dP   * u       + &
                   den         * por          * d2sat_dP2 * u       + &
                   den         * por          * dsat_dP   * du_dP   + &
                   dden_dP     * por          * sat       * du_dP   + &
                   den         * dpor_dP      * sat       * du_dP   + &
                   den         * por          * dsat_dP   * du_dP   + &
                   den         * por          * sat       * d2u_dP2 + &
                   rock_dencpr * (-d2por_dP2) * temp 

    d2energy_dPdt = d2den_dPdt  * por           * sat        * u        + &
                    dden_dP     * dpor_dt       * sat        * u        + &
                    dden_dP     * por           * dsat_dt    * u        + &
                    dden_dP     * por           * sat        * du_dt    + &
                    dden_dt     * dpor_dP       * sat        * u        + &
                    den         * d2por_dPdt    * sat        * u        + &
                    den         * dpor_dP       * dsat_dt    * u        + &
                    den         * dpor_dP       * sat        * du_dt    + &
                    dden_dt     * por           * dsat_dP    * u        + &
                    den         * dpor_dt       * dsat_dP    * u        + &
                    den         * por           * d2sat_dPdt * u        + &
                    den         * por           * dsat_dP    * du_dt    + &
                    dden_dt     * por           * sat        * du_dP    + &
                    den         * dpor_dt       * sat        * du_dP    + &
                    den         * por           * dsat_dt    * du_dP    + &
                    den         * por           * sat        * d2u_dPdt + &
                    rock_dencpr * (-d2por_dPdt) * temp                 + &
                    rock_dencpr * (-dpor_dP)  

    d2energy_dtdP = d2energy_dPdt


    Jlocal(2,1) = (shift*denergy_dP + &
                   d2energy_dP2*TH_auxvars(ghosted_id)%dpres_dtime + &
                   d2energy_dtdP*TH_auxvars(ghosted_id)%dtemp_dtime)* &
                   material_auxvars(ghosted_id)%volume


    ! Jlocal(2,2) = shift*d(A_E)/d(Tdot) + d(A_E)/d(T) 
    !             = shift*[d(rho*phi*s*U)/dT + d(rock_dencpr*(1-phi)*T)/dT]*Vol + 
    !               [d2(rho*phi*s*U)/dPdt + d2(rock_dencpr*(1-phi)*T)/dPdt]*dP_dtime*Vol +
    !               [d2(rho*phi*s*U)/dT2 + d2(rock_dencpr*(1-phi)*T)/dT2]*dT_dtime*Vol

    denergy_dt = dden_dt     * por        * sat     * u     + &
                 den         * dpor_dt    * sat     * u     + &
                 den         * por        * dsat_dt * u     + &
                 den         * por        * sat     * du_dt + &                          
                 rock_dencpr * (-dpor_dt) * temp            + &
                 rock_dencpr * (1-por)     
                 
    d2energy_dt2 = d2den_dt2   * por          * sat       * u       + &
                   dden_dt     * dpor_dt      * sat       * u       + &
                   dden_dt     * por          * dsat_dt   * u       + &
                   dden_dt     * por          * sat       * du_dt   + &
                   dden_dt     * dpor_dt      * sat       * u       + &
                   den         * d2por_dt2    * sat       * u       + &
                   den         * dpor_dt      * dsat_dt   * u       + &
                   den         * dpor_dt      * sat       * du_dt   + &
                   dden_dt     * por          * dsat_dt   * u       + &
                   den         * dpor_dt      * dsat_dt   * u       + &
                   den         * por          * d2sat_dt2 * u       + &
                   den         * por          * dsat_dt   * du_dt   + &
                   dden_dt     * por          * sat       * du_dt   + &
                   den         * dpor_dt      * sat       * du_dt   + &
                   den         * por          * dsat_dt   * du_dt   + &
                   den         * por          * sat       * d2u_dt2 + &
                   rock_dencpr * (-d2por_dt2) * temp              + &
                   rock_dencpr * (-dpor_dt)                       + &
                   rock_dencpr * (-dpor_dt)


    Jlocal(2,2) = (shift*denergy_dt + &
                   d2energy_dt2*TH_auxvars(ghosted_id)%dtemp_dtime + &
                   d2energy_dPdt*TH_auxvars(ghosted_id)%dpres_dtime)* &
                   material_auxvars(ghosted_id)%volume

    call MatSetValuesBlockedLocal(J,1,ghosted_id-1,1,ghosted_id-1,Jlocal, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)    
    
    
  enddo
    
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
  
  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

      
end subroutine IJacobianAccumulation

! ************************************************************************** !

subroutine PMTHTSInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/08/19
  !

  use TH_module, only : THInitializeTimestep
  
  implicit none
  
  class(pm_th_ts_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMTHTSInitializeTimestep

! ************************************************************************** !

subroutine PMTHTSCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)
  ! 
  ! Adds a convergence check for the nonlinear problem
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/08/19
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
  class(pm_th_ts_type) :: this
  PetscErrorCode :: ierr

  call PMTHCheckConvergence(this,snes,it,xnorm,unorm,fnorm,reason,ierr)

end subroutine PMTHTSCheckConvergence

! ************************************************************************** !

subroutine PMTHTSDestroy(this)
  ! 
  ! Destroys TH process model
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/08/19
  ! 
  use TH_module, only : THDestroy

  implicit none
  
  class(pm_th_ts_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call PMTHDestroy(this)
  
end subroutine PMTHTSDestroy

end module PM_TH_TS_class

