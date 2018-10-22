module Well_Solver_module

!  This well solver finds the current state of a well, automatically switching
!  to the controlling target mode.
!  The calculated well flows are place in the Pflotran residual
!  Well data is obtained from  the Well_Data class
!  Solution results are stored the Well_Data class

#include "petsc/finclude/petscsys.h"

  use petscsys
  use PFLOTRAN_Constants_module
  use Well_Data_class
  use PM_TOWG_Aux_module
  use Option_module

  implicit none

  private

    PetscInt,parameter :: max_itertt = 20 ! Max iterations over target types
    PetscInt,parameter :: max_iterws = 20 ! Max iterations for well solution

    PetscReal,parameter :: ws_eps_p = 1.0E-3  ! Pressure small shift value

    PetscReal :: ws_gravity                  ! Vertical gravity constant

    PetscBool :: ws_oil ! Bools indicating phases are present
    PetscBool :: ws_gas
    PetscBool :: ws_wat
    PetscBool :: ws_slv

    PetscInt :: ws_nphase                    ! Number of phases and components
    PetscInt :: ws_ncomp
    PetscInt :: ws_ncompe                    ! Number of components plus energy

    PetscInt :: w_type                       ! Well type (producer, gas inj...)

    PetscBool :: ws_isproducer = PETSC_TRUE ! Flags indicating producer or inj.
    PetscBool :: ws_isinjector = PETSC_FALSE

    PetscInt :: w_ncmpl =  0            ! Number of completions (this proc)
    PetscInt :: w_ncmplg = 0            ! Number of completions (all procs)

    PetscBool :: w_crossproc = PETSC_FALSE  ! Flag indicating cross-proc well

    PetscReal :: w_sign   ! Sign flag for well (+ve prod,-ve inje)

    PetscReal :: w_pw                           ! Well pressure at reference
    PetscReal :: w_pb                           ! Well bubble point pressure
    PetscReal,dimension(:),allocatable :: w_sp  ! Well saturations by phase
    PetscReal,dimension(:),allocatable :: w_bp  ! Well molar densities by phase
    PetscReal,dimension(:),allocatable :: w_gdp ! Well densities by phase

    PetscReal :: ws_pwbf ! Well pressure such that all completions backflow

    PetscReal :: w_zref  ! Reference elevation for well
    PetscReal :: w_gd    ! Gravity density (G*wellbore fluid density)

    PetscReal :: ws_targets(N_WELL_TT)         ! Well targets by target type
    PetscReal :: ws_actuals(N_WELL_TT)         ! Well actuals by target type

    PetscReal,dimension(:),allocatable :: ws_svpm ! Well surf. vol./mole by comp
    PetscReal,dimension(:),allocatable :: ws_mspm ! Well mass/mole for each comp

    PetscInt ,dimension(:),allocatable :: c_ghosted_id ! Cell ghosted id by compl
    PetscInt ,dimension(:),allocatable :: c_local_id   ! Cell local   id by compl
    PetscReal,dimension(:),allocatable :: c_ccf  ! Compl. conn. factor by compl.
    PetscReal,dimension(:),allocatable :: c_z    ! Cell elevation by completion

! Molar comp. and enegy flows by compl (this proc)
    PetscReal,dimension(:,:),allocatable :: c_flows

! For well,this proc:
    PetscReal,dimension(  :),allocatable :: w_flows    ! Molar component flows
    PetscReal,dimension(  :),allocatable :: w_flows_ms ! Mass  component flows
    PetscReal,dimension(  :),allocatable :: w_flows_sv ! Surf. vol. comp. flows

! For well, all procs
    PetscReal,dimension(  :),allocatable :: w_flowsG    ! Molar component flows
    PetscReal,dimension(  :),allocatable :: w_flowsG_ms ! Mass component flows
    PetscReal,dimension(  :),allocatable :: w_flowsG_sv ! Surf.vol. comp. flows

! Cell values (this proc)
    PetscReal,dimension(  :),allocatable :: c_p    ! Pressures
    PetscReal,dimension(  :),allocatable :: c_pb   ! Bubble points
    PetscReal,dimension(:,:),allocatable :: c_mdp  ! Molar densities
    PetscReal,dimension(:,:),allocatable :: c_sp   ! Saturations
    PetscReal,dimension(:,:),allocatable :: c_mob  ! Mobilities K/vsc
    PetscReal,dimension(  :),allocatable :: c_xo   ! Oil mole fraction)
    PetscReal,dimension(  :),allocatable :: c_xg   ! Gas mole fraction
    PetscReal,dimension(:,:),allocatable :: c_kgdp ! Mass density
    PetscReal,dimension(:,:),allocatable :: c_hp   ! Enthalpy density

    character(len=MAXSTRINGLENGTH) :: ws_name    ! Well name

    PetscInt :: w_comm  ! Well MPI communicator

  public ::  SolveWell,InitialiseWell,doWellMPISetup

contains

! ************************************************************************** !

subroutine InitialiseWell(well_data,grid,material_auxvars,option)
  !
  ! Do well initialisation, including completion connection factors
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  use Well_Data_class
  use Grid_module
  use Material_Aux_class

  implicit none

  class(well_data_type) , pointer :: well_data
  type (grid_type)      , pointer :: grid
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option

  PetscInt  :: icmpl,ncmpl,ixuser,iyuser,izuser,ix,iy,iz
  PetscInt  :: drilling_direction
  PetscReal :: dx,dy,dz,z,ccf,dx1,dx2,dh,k1,k2,radius,skinfactor,thetafactor,r0

  PetscInt :: lx,ly,lz,nlx,nly,nlz,ixl,ixu,iyl,iyu,izl,izu
  PetscBool :: onproc
  PetscInt :: local_id,ghosted_id

! Initialise

  ixuser = 0
  iyuser = 0
  izuser = 0

  nCmpl = well_data%GetNCmpl()

!  Store gravity constant

   ws_gravity = abs(option%gravity(3))

!  Do completion connection factor calculations

!  Get grid details

  lx = grid%structured_grid%lxs
  ly = grid%structured_grid%lys
  lz = grid%structured_grid%lzs

  nlx = grid%structured_grid%nlx
  nly = grid%structured_grid%nly
  nlz = grid%structured_grid%nlz

  ixl = lx+1
  iyl = ly+1
  izl = lz+1

  ixu = ixl+nlx-1
  iyu = iyl+nly-1
  izu = izl+nlz-1

!  Fill in completion data (from values set for wells etc.)

  call well_data%FillCmplData()

!  Loop over completions

  do icmpl = 1,ncmpl

!  Get the completion location

    call well_data%GetCmplLocation(icmpl,ixuser,iyuser,izuser, &
                                   local_id,ghosted_id)

!  Check if on-proc

    onproc = PETSC_FALSE
    if (      (ixuser >= ixl) .and. (ixuser <= ixu) &
        .and. (iyuser >= iyl) .and. (iyuser <= iyu) &
        .and. (izuser >= izl) .and. (izuser <= izu)  ) then
      onproc=PETSC_TRUE
      ix = ixuser-ixl+1
      iy = iyuser-iyl+1
      iz = izuser-izl+1
      local_id   = ix+(iy-1)*nlx+(iz-1)*nlx*nly
      ghosted_id = grid%nL2G(local_id)
      call well_data%SetCmplIndices(icmpl,local_id,ghosted_id)
    endif

!  if not on-proc then mark the completion for deletion

    if (.not.onproc) then
      call well_data%MarkCmplForDeletion(icmpl)
    endif
 enddo

!  Now delete the marked completions

 call well_data%DeleteMarkedCompletions()

!  Get new number of completions and generate ccfs

  ncmpl = well_data%GetNCmpl()
  do icmpl = 1,ncmpl

! Get completion location

    call well_data%GetCmplLocation(icmpl,ixuser,iyuser,izuser, &
                                   local_id,ghosted_id)

!  Get cell location

    ix = ixuser-ixl+1
    iy = iyuser-iyl+1
    iz = izuser-izl+1
    local_id   = ix+(iy-1)*nlx+(iz-1)*nlx*nly
    ghosted_id = grid%nL2G(local_id)

! Get cell dimensions

    dx = well_data%GetCmplDx(icmpl)
    dy = well_data%GetCmplDy(icmpl)
    dz = well_data%GetCmplDz(icmpl)

    if (dx<0.0) dx = grid%structured_grid%dx(ghosted_id)
    if (dy<0.0) dy = grid%structured_grid%dy(ghosted_id)
    if (dz<0.0) dz = grid%structured_grid%dz(ghosted_id)

    z = grid%z(ghosted_id)

!  Get cell directional properties

    drilling_direction = well_data%GetCmplDrillingDirection(icmpl)
    radius             = well_data%GetCmplRadius           (icmpl)
    skinfactor         = well_data%GetCmplSkinFactor       (icmpl)
    thetafactor        = well_data%GetCmplThetaFactor      (icmpl)

    select case(drilling_direction)
      case(X_DIRECTION)
        dx1 = dy
        dx2 = dz
        dh = dx
        k1 = material_auxvars(ghosted_id)%permeability(perm_yy_index)
        k2 = material_auxvars(ghosted_id)%permeability(perm_zz_index)
      case(Y_DIRECTION)
        dx1 = dx
        dx2 = dz
        dh = dy
        k1 = material_auxvars(ghosted_id)%permeability(perm_xx_index)
        k2 = material_auxvars(ghosted_id)%permeability(perm_zz_index)
      case(Z_DIRECTION)
        dx1 = dx
        dx2 = dy
        dh = dz
        k1 = material_auxvars(ghosted_id)%permeability(perm_xx_index)
        k2 = material_auxvars(ghosted_id)%permeability(perm_yy_index)
     end select

      r0 = (dx1**2.d0 * (k2/k1)**0.5d0 + dx2**2.d0 * (k1/k2)**0.5d0)**0.5d0 * &
           0.28d0 / ((k2/k1)**0.25d0 + (k1/k2)**0.25d0)

      ccf = 2.0d0 * PI * dh * dsqrt(k1*k2) * &
                      thetaFactor / &
                      ( dlog(r0/radius) + skinfactor )

      call well_data%SetCCF  (icmpl,ccf)
      call well_data%SetCmplZ(icmpl,z  )

  enddo

!  Fill in well reference z value

  call well_data%SetZRef(option)

end subroutine InitialiseWell

! *************************************************************************** !

subroutine doWellMPISetup(option,num_well,well_data_list)
  !
  ! Setup the communicators for the wells
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  type(option_type), pointer :: option

  PetscInt,intent(in) :: num_well
  type(well_data_list_type),pointer :: well_data_list

  PetscInt,allocatable :: ncpr   (:)
  PetscInt,allocatable :: ncprall(:)
  PetscInt,allocatable :: rfortw (:)

  PetscInt :: comm_size,iwell,ncmpl,nipr,irankp
  PetscMPIInt :: myrank,ierr,nrankw,nCmplG
  PetscBool :: ismp
  PetscMPIInt :: comm,group

  comm_size = option%mycommsize
  myrank    = option%myrank

!  Allocate the completion count buffers

  allocate(ncpr   (comm_size))
  allocate(ncprall(comm_size))
  allocate(rfortw (comm_size))

! Items per rank

  nipr = 1

!  Loop over wells

  do iwell = 1,num_well

!  Find which procs this well exists on

    comm  =0
    group =0

!  First set all rank completion counts to zero

    ncpr    = 0
    ncprall = 0
    rfortw  = 0

! Now increment the completion count for this rank

    ncmpl=GetWellNCmpl(iwell,well_data_list)
    ncpr(myrank+1) = ncmpl

! Do a gather to find out where the other rank completions are

    ierr = 0
    call MPI_ALLREDUCE(ncpr, ncprall,comm_size, &
                       MPI_INTEGER,MPI_SUM,option%mycomm,ierr)

!  Is this a multi=proc well?

    ncmplg = 0
    nrankw = 0

    ismp = PETSC_FALSE
    do irankp = 1,comm_size
      if ( ncprall(irankp)>0 ) then
        ncmplg = ncmplg+ncprall(irankp)
        nrankw = nrankw+1
        rfortw(nrankw) = irankp-1
      endif
    enddo

! If more than two ranks see the well, is multi-proc

    if ( nrankw >= 2 ) ismp = PETSC_TRUE

! If multi-proc, create the group and communicator

    if ( ismp ) then

      call MPI_GROUP_INCL (option%mygroup,nrankw,rfortw,group,ierr)
      call MPI_COMM_CREATE(option%mycomm,group,comm,ierr)

    endif

! Set global number of completions, multi-proc flag, group and communicator

    call WellSetGlobalInfo(iwell,nrankw,ncmplg,ismp,group,comm,well_data_list)

  enddo

  call WellSetGlobalInfoSet()

!  Free completion count buffers

  deallocate(ncpr   )
  deallocate(ncprall)
  deallocate(rfortw )

end subroutine doWellMPISetup

! *************************************************************************** !

subroutine SolveWell(aux,option,well_data,r_p)
  !
  ! Supervise well solution to obtain residual and Jacobian contributions
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Well_Data_class
  use EOS_Oil_module
  use EOS_Gas_module
  use EOS_Slv_module
  use Auxiliary_module

  implicit none

  type(auxiliary_type) :: aux
  type(option_type), pointer :: option
  class(well_data_type), pointer :: well_data
  PetscReal, pointer :: r_p(:)

  PetscInt :: ci = 1
  PetscInt :: cj = 1
  PetscInt :: ck = 1
  PetscInt :: local_id = 1
  PetscInt :: ghosted_id = 1
  PetscInt :: icmpl
  PetscBool :: finished,well_is_dead

  PetscBool :: welllocationfound
  PetscBool :: welltargetsfound
  PetscReal :: mw,sd
  PetscInt :: itertt,itt,jtt,rank,iphase
  PetscBool :: possible,well_solution_set

!  Basic setup

  welllocationfound = PETSC_FALSE
  welltargetsfound  = PETSC_FALSE
  possible          = PETSC_FALSE

  ws_nphase = option%nphase
  ws_ncomp  = option%nphase
  ws_ncompe = option%nphase+1

  ws_oil = PETSC_FALSE
  ws_gas = PETSC_FALSE
  ws_wat = PETSC_FALSE
  ws_slv = PETSC_FALSE

  do iphase = 1,ws_nphase
    if ( iphase .eq. option%oil_phase     ) ws_oil = PETSC_TRUE
    if ( iphase .eq. option%gas_phase     ) ws_gas = PETSC_TRUE
    if ( iphase .eq. option%liquid_phase  ) ws_wat = PETSC_TRUE
    if ( iphase .eq. option%solvent_phase ) ws_slv = PETSC_TRUE
  enddo

  w_type = well_data%GetType()
  if ( w_type .eq. PROD_WELL_TYPE ) then
    ws_isproducer = PETSC_TRUE
    ws_isinjector = PETSC_FALSE
  else
    ws_isproducer = PETSC_FALSE
    ws_isinjector = PETSC_TRUE
  endif

  if ( w_type .eq. PROD_WELL_TYPE ) then
    w_sign = 1.0
  else
    w_sign =-1.0
  endif

  w_zref = well_data%GetZRef()

!  Solve well if required

  w_ncmpl  = well_data%GetNCmpl ()
  w_ncmplg = well_data%GetNCmplG()
  w_comm   = well_data%GetWellComm()

  if ( (w_ncmpl > 0) .and. (w_ncmplg > w_ncmpl) ) then
    w_crossproc = PETSC_TRUE
  else
    w_crossproc = PETSC_FALSE
  endif

  if ( w_ncmpl > 0 ) then

!  Debug

    call well_data%GetWellName(ws_name)
    rank = option%myrank

!  Allocate the work arrays

    call allocateWorkArrays()

!  Extract completion data

    do icmpl =1,w_ncmpl
      call well_data%GetCmplLocation(icmpl,ci,cj,ck,local_id,ghosted_id)
      c_local_id  (icmpl) = local_id
      c_ghosted_id(icmpl) = ghosted_id
      c_ccf       (icmpl) = well_data%GetCCF  (icmpl)
      c_z         (icmpl) = well_data%GetCmplZ(icmpl)
    enddo

!  Load the cell properties into arrays indexed by completion

    if ( option%iflowmode == TOWG_MODE       ) then
      call wellSolverLoaderTOWG(aux,option)
    elseif ( option%iflowmode == TOIL_IMS_MODE ) then
      call wellSolverLoaderTOIL(aux,option)
    else
      call throwWellSolverException('This mode not yet supported by SolveWell')
    endif

!  Setup well solution

    well_solution_set = well_data%GetWellSolutionSet()
    if ( well_solution_set ) then
      call well_data%GetWellSolution(w_pw,w_pb,w_sp)
    else
!  Take initial solution from first completion
      w_pw = c_p (1)
      w_pb = c_pb(1)
      do iphase = 1,ws_nphase
        w_sp(iphase) = c_sp(1,iphase)
      enddo
    endif

!  Set up surface volume/mole values

    do iphase = 1,ws_nphase

      mw = 1.0
      sd = 1.0

      if ( iphase .eq. option%oil_phase     ) mw = EOSOilGetFMW()
      if ( iphase .eq. option%gas_phase     ) mw = EOSGasGetFMW()
      if ( iphase .eq. option%liquid_phase  ) mw = FMWH2O
      if ( iphase .eq. option%solvent_phase ) mw = EOSSlvGetFMW()

      if ( iphase .eq. option%oil_phase     ) sd = EOSOilGetReferenceDensity()
      if ( iphase .eq. option%gas_phase     ) sd = EOSGasGetReferenceDensity()
      if ( iphase .eq. option%liquid_phase  ) then
        sd = option%reference_density(option%liquid_phase)
      endif
      if ( iphase.eq.option%solvent_phase ) sd = EOSSlvGetReferenceDensity()

      ws_svpm(iphase) = mw/sd
      ws_mspm(iphase) = mw
    enddo

!  Extract well status data

    welltargetsfound = well_data%GetTargets(ws_targets)
    itt=well_data%GetTT()

!  Find well gravity density based on cell properties

    call findWellboreGravityDensity()

!  Find backflow well pressure

   call findBackflowWellPressure()
   well_is_dead = PETSC_FALSE
   if( ws_isProducer ) then
    if( ws_targets(W_BHP_LIMIT) > ws_pwbf ) well_is_dead = PETSC_TRUE
   endif
   if( ws_isInjector ) then
    if( ws_targets(W_BHP_LIMIT) < ws_pwbf ) well_is_dead = PETSC_TRUE
   endif

!  Either solve well or zero rates

   if( .not.well_is_dead ) then

!  Loop over attempts to find the well target type

      do itertt = 1,max_itertt
        finished = solveForWellTarget(option,itertt,itt)
        if ( finished ) exit
      enddo
      if ( .not.finished ) then
        call throwWellSolverException('Unable to converge target iteration')
      endif

!  Update residual

      call updateMainResidual(option,r_p)

!  Get actual rates for each target for this proc (will be globalised in well_data)

      do jtt = 1,N_WELL_TT
        ws_actuals(jtt) = getActualFlowForTargetType(option,jtt,possible, &
                                                     w_flows_sv,w_flows_ms)
      enddo

!  Update stored well status

      call well_data%SetTT(itt)
      call well_data%SetNComp(ws_ncomp)
      call well_data%SetWellFlows(w_flows        ,ws_ncomp)
      call well_data%SetCmplFlows(c_flows,w_ncmpl,ws_ncomp)
      call well_data%SetActuals(ws_actuals)

    else

      itt = W_BHP_LIMIT
      call well_data%SetTT(itt)
      call well_data%SetNComp(ws_ncomp)
      call well_data%ZeroWellFlows(ws_ncomp)
      call well_data%ZeroCmplFlows(w_ncmpl,ws_ncomp)
      call well_data%ZeroActuals()

    endif

!  Update stored well solution

    call well_data%SetWellSolution(w_pw,w_pb,w_sp)
    call well_data%SetWellSolutionSet()

!  Free the work arrays

    call freeWorkArrays()

  endif

end subroutine SolveWell

! *************************************************************************** !

function solveForWellTarget(option,itertt,itt)
  !
  ! Solve for a given well target
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  PetscBool :: solveForWellTarget
  type(option_type), pointer :: option
  PetscInt,intent(in   ) :: itertt
  PetscInt,intent(inout) :: itt

  PetscInt :: iterws,jtt,nconv
  PetscBool :: finished,possible
! fpsc: flow with positive sign convention
  PetscReal::conv_crit,rw,jww,jwwi,fpsc,ftarg,pwtarg,bhp_target

  conv_crit = 1.0E-6
  rw        = 0.0
  jww       = 1.0

!  Do the inner iteration for the well solution when on this target

  finished = PETSC_FALSE
  nconv = 0
  do iterws = 1,max_iterws

!  Check that solution is within range

    if(itt .ne. W_BHP_LIMIT) then

      bhp_target = ws_targets(W_BHP_LIMIT)
      if (ws_isproducer) then
        if (w_pw < bhp_target) w_pw = bhp_target
! Leave space for differencing
        if (w_pw > ws_pwbf) w_pw = ws_pwbf-2.0*ws_eps_p
      endif
      if (ws_isinjector) then
        if (w_pw > bhp_target) w_pw = bhp_target
        if (w_pw < ws_pwbf) w_pw = ws_pwbf
      endif
    endif

!  Get the residual and Jacobian, exit if converged

    call getWellResidualAndJacobian(option,iterws,itt,rw,jww)
    if (abs(rw) < conv_crit) nconv = nconv+1
    if (nconv > 3) then
      finished = PETSC_TRUE
      exit
    endif

!  Find inverse Jacobian

    if (abs(jww) > 0.0) then
      jwwi = 1.0/jww
    else
      jwwi = 0.0
    endif

!  Update solution

    w_pw = w_pw - rw*jwwi

  enddo

  if (.not.finished) then
    call throwWellSolverException('Unable to converge well iteration')
  endif

!  Check if all other targets satisfied

  pwtarg = ws_targets(W_BHP_LIMIT)

  solveForWellTarget = PETSC_TRUE
  do jtt = 1,N_WELL_TT
! No need to check current type
    if (jtt .eq. itt) cycle
! Case of bhp control
    if (jtt .eq. W_BHP_LIMIT) then
      if (     ( ws_isproducer .and. ( w_pw<pwtarg )) &
          .or. ( ws_isinjector .and. ( w_pw>pwtarg ))) then
!  Switch to bhp control
        itt = jtt
        solveForWellTarget = PETSC_FALSE
      endif
    else
! Case of rate control - get positive convention flow
      fpsc = w_sign*getActualFlowForTargetType(option,jtt,possible, &
                                             w_flowsG_sv,w_flowsG_ms)
      if (possible) then
        ftarg = ws_targets(jtt)
        if (ftarg > -0.5) then
          if (fpsc > ftarg) then
! Non-defaulted flow target exceeded - switch to this target
            itt = jtt
            solveForWellTarget = PETSC_FALSE
          endif
        endif
      endif
    endif
  enddo

end function solveForWellTarget

! *************************************************************************** !

subroutine getWellResidualAndJacobian(option,iterws,itt,rw,jww)
  !
  ! Do one iteration of well solution
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  type(option_type), pointer :: option
  PetscInt,intent(in)::iterws
  PetscInt,intent(in)::itt
  PetscReal,intent(out)::rw,jww
  PetscReal::rwe

!  Build up the residual for this target and solution

  rwe = 0.0
  rw  = 0.0

  w_pw = w_pw + ws_eps_p
  call findWellFlows(option)
  rwe = extractResidual(option,itt)
  w_pw = w_pw - ws_eps_p

  call findWellFlows(option)
  rw = extractResidual(option,itt)

  jww = (rwe-rw) / ws_eps_p

end subroutine getWellResidualAndJacobian

!***************************************************************************************!

function extractResidual(option,itt)
  !
  ! Get residual for well solver
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  PetscReal :: extractResidual
  type(option_type), pointer :: option
  PetscInt,intent(in) :: itt
  PetscBool :: possible

  PetscReal :: tact,treq,rw

  rw = 0.0
  possible = PETSC_FALSE

  if ( itt .eq. W_BHP_LIMIT ) then
    treq = ws_targets(itt)
    rw = w_pw- treq
  else
    tact= getActualFlowForTargetType(option,itt,possible, &
                                     w_flowsG_sv,w_flowsG_ms)
    treq = ws_targets(itt)
    rw = w_sign*tact - treq
    extractResidual = rw
  endif

  extractResidual = rw

end function extractResidual

! *************************************************************************** !

function getActualFlowForTargetType(option,itt,possible,flows_sv,flows_ms)
  !
  ! Get actual flow for a given target type
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  PetscReal :: getActualFlowForTargetType
  type(option_type), pointer :: option
  PetscInt,intent(in) :: itt
  PetscBool,intent(out) :: possible
  PetscReal,dimension(:) :: flows_sv
  PetscReal,dimension(:) :: flows_ms

  PetscReal :: f

  f = 0.0
  possible = PETSC_FALSE

  if( itt .eq. W_BHP_LIMIT ) then
    f = w_pw
  endif

  if ((itt .eq. W_TARG_OSV) .and. ws_oil) then
    f = flows_sv(option%oil_phase    )
    possible = PETSC_TRUE
  endif

  if ((itt.eq. W_TARG_GSV) .and. ws_gas) then
    f = flows_sv(option%gas_phase    )
    possible = PETSC_TRUE
  endif

  if ((itt .eq. W_TARG_WSV).and. ws_wat) then
    f = flows_sv(option%liquid_phase )
    possible = PETSC_TRUE
  endif

  if ((itt .eq. W_TARG_SSV).and. ws_slv) then
    f = flows_sv(option%solvent_phase)
    possible = PETSC_TRUE
  endif

  if ((itt .eq. W_TARG_LSV).and. ws_oil .and. ws_wat) then
    f =  flows_sv(option%oil_phase    ) &
        +flows_sv(option%liquid_phase )
    possible = PETSC_TRUE
  endif

  if ((itt .eq. W_TARG_OM).and. ws_oil) then
    f = flows_ms(option%oil_phase    )
    possible = PETSC_TRUE
  endif

  if ((itt .eq. W_TARG_GM).and. ws_gas) then
    f = flows_ms(option%gas_phase    )
    possible = PETSC_TRUE
  endif

  if ((itt .eq. W_TARG_WM).and. ws_wat) then
    f = flows_ms(option%liquid_phase )
    possible = PETSC_TRUE
  endif

  if ((itt .eq. W_TARG_SM).and. ws_slv) then
    f = flows_ms(option%solvent_phase)
    possible = PETSC_TRUE
  endif

  getActualFlowForTargetType = f

end function getActualFlowForTargetType

! ************************************************************************** !

subroutine updateMainResidual(option,r_p)
  !
  ! Insert the well flow contributiosn into the mail Pflotran residual
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  type(option_type), pointer :: option
  PetscReal, pointer :: r_p(:)

  PetscReal :: res(option%nflowdof)

  PetscInt :: icmpl,icompe,local_id,istart,iend

!  Loop over completions

  do icmpl = 1,w_ncmpl

    res = 0.0

    local_id =c_local_id(icmpl)
    do icompe = 1,ws_ncompe
      res(icompe) = c_flows(icmpl,icompe)
    enddo

    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1

    r_p(istart:iend) = r_p(istart:iend) + res(1:option%nflowdof)
  enddo

end subroutine updateMainResidual

! *************************************************************************** !

subroutine allocateWorkArrays()
  !
  ! Allocate well solve work arrays
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  allocate(c_local_id  (w_ncmpl          ))
  allocate(c_ghosted_id(w_ncmpl          ))
  allocate(c_ccf       (w_ncmpl          ))
  allocate(c_z         (w_ncmpl          ))

  allocate(ws_svpm     (        ws_ncomp ))
  allocate(ws_mspm     (        ws_ncomp ))

  allocate(c_flows     (w_ncmpl,ws_ncompe)) ! Note extra energy location
  allocate(w_flows     (        ws_ncomp ))
  allocate(w_flows_ms  (        ws_ncomp ))
  allocate(w_flows_sv  (        ws_ncomp ))
  allocate(w_flowsG    (        ws_ncomp ))
  allocate(w_flowsG_ms (        ws_ncomp ))
  allocate(w_flowsG_sv (        ws_ncomp ))

  allocate(c_p         (w_ncmpl          ))
  allocate(c_pb        (w_ncmpl          ))
  allocate(c_mdp       (w_ncmpl,ws_nphase))
  allocate(c_sp        (w_ncmpl,ws_nphase))
  allocate(c_mob       (w_ncmpl,ws_nphase))
  allocate(c_xo        (w_ncmpl          ))
  allocate(c_xg        (w_ncmpl          ))
  allocate(c_kgdp      (w_ncmpl,ws_nphase))
  allocate(c_hp        (w_ncmpl,ws_nphase))

  allocate(w_sp        (        ws_nphase))
  allocate(w_bp        (        ws_nphase))
  allocate(w_gdp       (        ws_nphase))

end subroutine allocateWorkArrays

! *************************************************************************** !

subroutine freeWorkArrays()
  !
  ! De-allocate well solve work arrays
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  deallocate(c_local_id  )
  deallocate(c_ghosted_id)
  deallocate(c_ccf       )
  deallocate(c_z         )

  deallocate(ws_svpm     )
  deallocate(ws_mspm     )

  deallocate(c_flows     )
  deallocate(w_flows     )
  deallocate(w_flows_ms  )
  deallocate(w_flows_sv  )
  deallocate(w_flowsG    )
  deallocate(w_flowsG_ms )
  deallocate(w_flowsG_sv )

  deallocate(c_P         )
  deallocate(c_pb        )
  deallocate(c_mdp       )
  deallocate(c_sp        )
  deallocate(c_mob       )
  deallocate(c_xo        )
  deallocate(c_xg        )
  deallocate(c_kgdp      )
  deallocate(c_hp        )

  deallocate(w_sp        )
  deallocate(w_bp        )
  deallocate(w_gdp       )

end subroutine freeWorkArrays

! *************************************************************************** !

subroutine findWellFlows(option)
  !
  ! Find well flows
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  type(option_type), pointer :: option
  PetscInt :: ghosted_id,icomp,icmpl,ierr

!  Initialise the output arrays

  c_flows     = 0.0
  w_flows     = 0.0
  w_flows_ms  = 0.0
  w_flows_sv  = 0.0
  w_flowsG    = 0.0
  w_flowsG_ms = 0.0
  w_flowsG_sv = 0.0

!  Get completion flows

  do icmpl = 1,w_ncmpl
   ghosted_id = c_ghosted_id(icmpl)
   call findCompletionFlows(option,icmpl)
  enddo

!  Build well flows

  do icomp = 1,ws_ncomp
    do icmpl = 1,w_ncmpl
      w_flows   (icomp) = w_flows   (icomp)+               c_flows(icmpl,icomp)
      w_flows_sv(icomp) = w_flows_sv(icomp)+ws_svpm(icomp)*c_flows(icmpl,icomp)
      w_flows_ms(icomp) = w_flows_ms(icomp)+ws_mspm(icomp)*c_flows(icmpl,icomp)
    enddo
  enddo

!  Build well flows across procs if required

  if (w_crossproc) then
    call MPI_Allreduce(w_flows   ,w_flowsG   ,ws_ncomp, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,w_comm,ierr)
    call MPI_Allreduce(w_flows_sv,w_flowsG_sv,ws_ncomp, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,w_comm,ierr)
    call MPI_Allreduce(w_flows_ms,w_flowsG_ms,ws_ncomp, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,w_comm,ierr)
  else
    w_flowsG    = w_flows
    w_flowsG_sv = w_flows_sv
    w_flowsG_ms = w_flows_ms
  endif

end subroutine findWellFlows

! *************************************************************************** !

subroutine findCompletionFlows(option,icmpl)
  !
  ! Find completion flows
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  type(option_type), pointer :: option
  PetscInt ,intent(in) :: icmpl

  PetscBool :: is_black_oil,componentInPhase,is_oil_in_oil,is_gas_in_oil
  PetscInt :: iphase,icomp,eid
  PetscReal :: xmf,ccf,z,dhh,pdd,mob,mobt,mdp,flow

!  Initialise

  eid = option%energy_id

!  Extract completion connection factor and drawdown

  ccf = c_ccf(icmpl)
  z   = c_z  (icmpl)
  dhh = z-w_zref
  pdd = c_P  (icmpl)-w_pw+w_gd*dhh

!  Check for black oil cases (with dissolved gas)

  is_black_oil = PETSC_FALSE
  if (     ( towg_miscibility_model == TOWG_BLACK_OIL) &
      .or. ( towg_miscibility_model == TOWG_SOLVENT_TL)) is_black_oil = PETSC_TRUE

  mobt = 0.0
  do iphase = 1, ws_nphase
    mobt = mobt + c_mob(icmpl,iphase)
  enddo

!  Loop over phases and components building up component flows

  if (is_black_oil) then
    do iphase = 1, ws_nphase
      do icomp = 1, ws_ncomp

        componentInPhase = PETSC_FALSE
        is_oil_in_oil    = PETSC_FALSE
        is_gas_in_oil    = PETSC_FALSE

! disgas is special case

        if (iphase == option%oil_phase) then
          if (icomp == option%oil_phase) is_oil_in_oil = PETSC_TRUE
          if (icomp == option%gas_phase) is_gas_in_oil = PETSC_TRUE
        endif

! OK if phase and component match or dissolved gas

        if (iphase == icomp .or. (ws_isproducer .and. is_gas_in_oil)) then
          componentInPhase = PETSC_TRUE
        endif

        if (componentInPhase) then
          xmf = 1.0d0
          if (ws_isproducer) then
            if (is_oil_in_oil) xmf = c_xo(icmpl)
            if (is_gas_in_oil) xmf = c_xg(icmpl)
          endif

!  Initialise mobility and get molar density of phase

          mob =0.0
          mdp =c_mdp(icmpl,iphase)

! Check for backflow (leave mob at zero otherwise)

          if (     ( ws_isproducer .and. (Pdd > 0.0)) &
              .or. ( ws_isinjector .and. (Pdd < 0.0)) ) then

!  Set up mobility estimate

            if (ws_isproducer) mob = c_mob(icmpl,iphase)
            if (      (w_type .eq. OIL_INJ_WELL_TYPE) &
                .and. (icomp.eq.option%oil_phase    )) mob = mobt
            if (      (w_type .eq. GAS_INJ_WELL_TYPE) &
                .and. (icomp.eq.option%gas_phase    )) mob = mobt
            if (      (w_type .eq. WAT_INJ_WELL_TYPE) &
                .and. (icomp.eq.option%liquid_phase )) mob = mobt
            if (      (w_type .eq. SLV_INJ_WELL_TYPE) &
                .and. (icomp.eq.option%solvent_phase)) mob = mobt
          endif

!  Find flow

          flow = ccf*xmf*mob*mdp*pdd
          c_flows(icmpl,icomp) = c_flows(icmpl,icomp)+flow

!  Increment energy flow by component

          if(iphase == icomp) then
            c_flows(icmpl,eid)=c_flows(icmpl,eid) + flow*c_hp(icmpl,iphase)
          endif
        endif

      enddo
    enddo

  else

!  Non-black-oil case

    do iphase = 1, ws_nphase

      icomp = iphase ! Component and phase indices are the same in non-black-oil case

      if (ws_isproducer) mob = c_mob(icmpl,iphase)
      if (      w_type .eq. OIL_INJ_WELL_TYPE &
          .and. iphase.eq.option%oil_phase    ) mob = mobt
      if (      w_type .eq. GAS_INJ_WELL_TYPE &
          .and. iphase.eq.option%gas_phase    ) mob = mobt
      if (      w_type .eq. WAT_INJ_WELL_TYPE &
          .and. iphase.eq.option%liquid_phase ) mob = mobt
      if (      w_type .eq. SLV_INJ_WELL_TYPE &
          .and. iphase.eq.option%solvent_phase) mob = mobt

      mdp = c_mdp(icmpl,iphase)

      flow = ccf*mob*mdp*pdd
      c_flows(icmpl,icomp) = c_flows(icmpl,icomp)+flow

    enddo

  endif

end subroutine findCompletionFlows

! *************************************************************************** !

subroutine findWellboreGravityDensity()
  !
  ! Find gravity*(wellbore mass density)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  PetscReal :: wd,sumgdm,sumgds,summ,sums,sp,gdp,mp
  PetscInt :: icmpl,iphase,ierr
  PetscReal,dimension(FOUR_INTEGER) :: sl
  PetscReal,dimension(FOUR_INTEGER) :: rg

  wd = 0.0

!  Form saturation and mobility weighted sums

  sumgdm = 0.0
  sumgds = 0.0

  summ = 0.0
  sums = 0.0

  do icmpl=1,w_ncmpl
    do iphase=1,ws_nphase
      sp  = c_sp  (icmpl,iphase)
      mp  = c_mob (icmpl,iphase)
      gdp = c_kgdp(icmpl,iphase)
      sumgdm = sumgdm+mp*gdp
      sumgds = sumgds+sp*gdp
      summ   = summ  +mp
      sums   = sums  +sp
    enddo
  enddo

!  In case of cross-proc well continue summation over processors

  if (w_crossproc) then
    ierr= 0
    rg = 0.0
    sl(1) = sumgdm
    sl(2) = sumgds
    sl(3) = summ
    sl(4) = sums
    call MPI_ALLREDUCE(sl,rg,FOUR_INTEGER, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,w_comm,ierr)
    sumgdm = rg(1)
    sumgds = rg(2)
    summ   = rg(3)
    sums   = rg(4)
  endif

!  For density in Kg/rm3

  if (summ > 0.0) then
    wd = sumgdm/summ
  else
    wd = sumgds/sums
  endif

!--Store result with gravity constant included

  w_gd = wd*ws_gravity

end subroutine findWellboreGravityDensity

! *************************************************************************** !

subroutine findBackflowWellPressure()
  !
  ! Find the pressure at which the well will attempt to backflow
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  implicit none

  PetscInt :: icmpl,ierr
  PetscReal :: pwbf,z,dhh,pwbfc
  PetscReal,dimension(1) :: sl
  PetscReal,dimension(1) :: rg

  ierr = 0

!  For producer, find maximum of pw values at which injection starts

  if (ws_isproducer) then
    pwbf = 0.0
    do icmpl = 1,w_ncmpl
      z   = c_z  (icmpl)
      dhh = z-w_zref
      pwbfc =c_P(icmpl)+w_gd*dhh
      if( pwbfc > pwbf ) pwbf = pwbfc
    enddo

!  Take cross-proc value if required

    if (w_crossproc) then
      rg =0.0
      sl(1) = pwbf
      call MPI_ALLREDUCE(sl,rg,ONE_INTEGER, &
                         MPI_DOUBLE_PRECISION,MPI_MAX,w_comm,ierr)
      pwbf = rg(1)
    endif

  endif

!  For injector, find minimum of pw values at which production start

  if (ws_isinjector) then
    pwbf = 1.0D10
    do icmpl = 1,w_ncmpl
      z   = c_z  (icmpl)
      dhh = z-w_zref
      pwbfc = c_P(icmpl)+w_gd*dhh
      if( pwbfc < pwbf ) pwbf = pwbfc
    enddo

!  Take cross-proc value if required

    if (w_crossproc) then
      rg = 0.0
      sl(1) = pwbf
      call MPI_ALLREDUCE(sl,rg,ONE_INTEGER, &
                         MPI_DOUBLE_PRECISION,MPI_MIN,w_comm,ierr)
      pwbf = rg(1)
    endif

  endif

  ws_pwbf = pwbf

end subroutine findBackflowWellPressure

! *************************************************************************** !

subroutine wellSolverLoaderTOWG(aux,option)
  !
  ! Load completion solution arrays for TOWG mode
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Auxiliary_module
  use PM_TOWG_Aux_module

  implicit none

  type(auxiliary_type) :: aux
  type(option_type), pointer :: option

  class(pm_towg_aux_type), pointer :: TOWG

  PetscInt :: icmpl,ghosted_id

  towg=>aux%TOWG

  do icmpl = 1,w_ncmpl
   ghosted_id = c_ghosted_id(icmpl)
   call loadCellDataTOWG(towg%auxvars(ZERO_INTEGER,ghosted_id),option,icmpl)
  enddo

end subroutine wellSolverLoaderTOWG

!*****************************************************************************!

subroutine loadCellDataTOWG(auxvar,option,icmpl)
  !
  ! Load completion cell arrays for TOWG mode (for one cell)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

  use AuxVars_TOWG_module

  implicit none

  class(auxvar_towg_type) :: auxvar
  type(option_type), pointer :: option
  PetscInt ,intent(in) :: icmpl
  PetscInt :: iphase

  c_p (icmpl) = auxvar%pres(option%oil_phase)
  c_pb(icmpl) = auxvar%bo%bubble_point

  do iphase = 1, ws_nphase
    c_mdp (icmpl,iphase) = auxvar%den     (iphase)
    c_sp  (icmpl,iphase) = auxvar%sat     (iphase)
    c_mob (icmpl,iphase) = auxvar%mobility(iphase)
    c_kgdp(icmpl,iphase) = auxvar%den_kg  (iphase)
    c_hp  (icmpl,iphase) = auxvar%H       (iphase)
  enddo

  c_xo(icmpl) = auxvar%bo%xo
  c_xg(icmpl) = auxvar%bo%xg

end subroutine loadCellDataTOWG

! **************************************************************************** !

subroutine wellSolverLoaderTOIL(aux,option)
  !
  ! Well solver loader for TOIl mode (not yet supported)
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Auxiliary_module
  use PM_TOilIms_Aux_module

  implicit none

  type(auxiliary_type) :: aux
  type(option_type), pointer :: option

  class(pm_toil_ims_aux_type), pointer :: TOil_ims

  TOil_ims => aux%TOil_ims

  call throwWellSolverException('TOIL mode not yet supported by SolveWell')

end subroutine wellSolverLoaderTOIL

! *************************************************************************** !

subroutine throwWellSolverException(message)
  !
  ! Issue a general well solver error
  !
  ! Author: Dave Ponting
  ! Date: 08/15/18

  implicit none

  character(len=*) :: message
  print *,message

end subroutine throwWellSolverException

end module Well_Solver_module

